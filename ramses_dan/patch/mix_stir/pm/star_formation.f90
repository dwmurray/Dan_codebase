!################################################################
!################################################################
!################################################################
!################################################################
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, twopi
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled. 
  ! It modifies hydrodynamic variables according to mass conservation 
  ! and assumes an isothermal transformation... 
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,e0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer ::ntot,ntot_all,info,nstar_corrected,ncell
!#ifdef SOLVERhydro
!  integer ::imetal=6
!#endif
!#ifdef SOLVERmhd
!  integer ::imetal=9
!#endif
  logical ::ok_free,ok_all
  real(dp)::d,x,y,z,u,v,w,e,zg,ppf,pgz,vdisp,dgas,feg,oxg
  real(dp)::mstar,dstar,tstar,nISM,nCOM
  real(dp)::tot_nlw_photons, tot_nHgas, tot_nHstar_this_step ! RS - h2_frac decay
  real(dp)::tot_nHgas_all, tot_nHstar_this_step_all ! MPI tots
  real(dp)::velc,uc,vc,wc,mass_load,ul,vl,wl,ur,vr,wr,divv,curlv,alpha
  real(dp)::vxgauss,vygauss,vzgauss,birth_epoch,factG
  real(kind=8)::mlost,mtot,mlost_all,mtot_all
  real(kind=8)::RandNum,GaussNum,PoissMean   
  real(dp)::vsn,costheta,sintheta,phi,cosphi,sinphi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,T2,nH,T_poly
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  integer ,dimension(1:nvector,0:twondim)      :: ind_nbor
!  real(dp),dimension(1:nchem)::chem1=0.0
  
  if(numbtot(1,ilevel)==0) return ! Number of grids in level (I beleive, RS)
  if(.not. hydro)return
  if(ndim.ne.3)return


  !if(mstar_tot.gt.0) return
  if(verbose)write(*,*)' Entering star_formation'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only
  t0=t_star*(1d9*365.*24.*3600.)/scale_t

  ! ISM density threshold from H/cc to code units
  ! R. Sarmento
  ! Here's where we look at n_star and del_star, star
  ! formation threshold is the larger of the two.
  nISM = n_star
  nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
  ! So nCOM is the threshold overdensity at time "aexp" in terms of H/cc
  ! Note the use of omega_b vs omega_m in the units.f90 package
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH

  ! Initial star particle mass
  mstar=n_star/(scale_nH*aexp**3)*vol_min
  dstar=mstar/vol_loc
  factG = 1d0
  if (cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  !mstar=MAX(del_star*omega_b*rhoc*XH/mH,n_star)/(scale_nH*aexp**3)*vol_min
  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

#if NDIM==3
  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)/d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e-0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e
        end do
        ! Rick Sarmento
        ! Converting our variables -- dividing out the density
        ! Looks like this is done so we can use the actual
        ! values of the uold scalar array after 'this' pt in the code... 
        ! nvar is the total number of variables - so we're
        ! starting with the metals, and continuing up with
        ! the new scalars I've created
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

  ! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
!  ndebris_tot=0
  tot_nHgas=0.0D0  ! R. Sarmento -- Total # gas particles
  tot_nHstar_this_step = 0.0D0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Density criterion
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           if(d<=d0)ok(i)=.false. 
        end do
        ! Temperature criterion
        do i=1,ngrid
           T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0)
           nH=uold(ind_cell(i),1)*scale_nH
           T_poly=T2_star*(nH/nISM)**(g_star-1.0)
           T2=T2-T_poly
           if(T2>1.1*T_poly)ok(i)=.false. 
        end do
        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           ! Compute total gas mass... not just for leaves. Moved outside
           ! of ok check
           d=uold(ind_cell(i),1)
           mcell=d*vol_loc 
           tot_nHgas = tot_nHgas + (mcell * scale_nH * scale_l**3)! RS - tot gas mass as # particles
           if(ok(i))then
              ! Compute mean number of events
!              print *,"Cell mass ", mcell
!              print *,"Cell mass [g]", mcell * scale_d * scale_l**3
!              print *,"Cell mass [#H]", mcell * scale_nH * scale_l**3
              tstar  = t0*sqrt(d0/d)
              PoissMean = dtnew(ilevel)/tstar*mcell/mstar
              ! If catastrophic star formation (massive star cluster) wants to occur, we need to limit the 
              ! maximal mass of the star particle we want to create in a cell for numerical reasons (gravity solver).
              PoissMean = min(PoissMean,10.0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas=nstar(i)*mstar
              ! Security to prevent more than 90% of gas depletion
              if (mgas > 0.9*mcell) then
                 nstar_corrected=int(0.9*mcell/mstar)
                 mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                 nstar(i)=nstar_corrected
                 mgas = 0.9 * mcell ! RS - correct mgas needed for gas loss
              endif
              tot_nHgas = tot_nHgas - (mgas * scale_nH * scale_l**3)! Account for star loss
              ! Compute new stars local statistics
              tot_nHstar_this_step = tot_nHstar_this_step + (nstar(i)*mstar * scale_nH * scale_l**3) ! Just stars created 'now'
              mstar_tot=mstar_tot+nstar(i)*mstar
              ! Print info for this level - RS
!!$              print *,"ilevel ", ilevel
!!$              print *,"aexp ",aexp
!!$              print *,"dx ",dx
!!$              print *,"dx (pc)",dx * scale_l / 3.085D18
!!$              print *,"dx loc ",dx_loc
              if(nstar(i)>0)then
                 ntot=ntot+1
              endif
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ! R. Sarmento - need to get other processors totals
  call MPI_ALLREDUCE(tot_nHgas,tot_nHgas_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tot_nHstar_this_step,tot_nHstar_this_step_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF4.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/mtot_all*100.
        ! R. Sarmento - Update the H2 fraction based on the # of star
        ! particles created in this time-step
        if (h2_frac > 0.0d0 .and. tot_nHgas_all > 0.0d0 .and. tot_nHstar_this_step_all > 0.0d0)then
           ! Assume every photon dissociates H2 -- May want 
           ! to change this to 0.15 since only 15% of the time 
           ! does absorption result in dissociation
           ! Compute: 
           !      h2_frac = (h2_frac * tot_nHgas - tot_nlw_photons)/tot_nHgas
           tot_nlw_photons = tot_nHstar_this_step_all * lw_photon
           print *,"H2 frac BEFORE ()", h2_frac
           h2_frac = h2_frac - (tot_nlw_photons/tot_nHgas_all)
           print *,"Total stellar numB ()", tot_nHstar_this_step_all
           ! Convert to # of baryons... assume we've taken care of He fraction
           print *,"Total stellar mass (Msol)", tot_nHstar_this_step_all / scale_nH * scale_d / (1.988d33)
           print *,"Total gas numB ()", tot_nHgas_all
           print *,"Total gas mass (Msol)", tot_nHgas_all / scale_nH * scale_d / (1.988d33)
           print *,"Total lw photons ()", tot_nlw_photons
           if (h2_frac < 0.0d0)h2_frac = 0.0d0
           print *,"H2 frac AFTER ()", h2_frac
        end if ! End h2_frac        
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0 ! If flag2 > 0 we have sp's in grid i
        end do

        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i)) ! n is number of sp's in the grid
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale ! xg is grid positions
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
!!$           print *,"star form: position - ", x,y,z
           if(metal)then
              zg=uold(ind_cell_new(i),imetal)
!              do ichem=1,nchem
!                 chem1(ichem)=uold(ind_cell_new(i),imetal+ichem)
!              enddo
              ! Rick Sarmento - Jul 2014
              ! Also capture the pristine gas fraction for this cell
              ! if we are tracking pristine fraction. 
              if(prist_gas_fraction) then
                 ppf=uold(ind_cell_new(i),iprist)
                 pgz=uold(ind_cell_new(i),iprimordz)
              endif
           endif

           ! Set star particle variables
           tp(ind_part(i))=birth_epoch  ! Birth epoch
           mp(ind_part(i))=n*mstar      ! Mass
           levelp(ind_part(i))=ilevel   ! Level
           idp(ind_part(i))=-index_star  ! Star identity
!!$           print *, "star form ind_part, i", ind_part(i), i
           xp(ind_part(i),1)=x
           xp(ind_part(i),2)=y
           xp(ind_part(i),3)=z
           vp(ind_part(i),1)=u
           vp(ind_part(i),2)=v
           vp(ind_part(i),3)=w
           if(metal)then
              zp(ind_part(i))=zg  ! Initial star metallicity
!              do ichem=1,nchem
!                 chp(ind_part(i),ichem)=chem1(ichem)  ! Initial chemical abudance
!              enddo
              ! Track the fraction of stars that are pristine - Rick Sarmento
              if(prist_gas_fraction)then
                 ! Since we just lost gas from the cell, we need to adjust the
                 ! pristine gas density ... it reflects the pristine fraction at
                 ! the density before the star was formed.
                 pfp(ind_part(i))=ppf ! RS - The star gets the pristine fraction of the gas
                 z3p(ind_part(i))=pgz ! RS - The star gets the gas primord Z
!!$                 print *,"star_formation: Setting particle Z = ",zg
!!$                 print *,"star_formation: Setting particle ppf = ",ppf
!!$                 print *,"star_formation: Setting particle z3p = ",pgz
              endif
           endif
        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           if(flag2(ind_cell(i))>0)then ! flag2 > 0 indicates SNe in the cell
              n=flag2(ind_cell(i)) ! n is the number of SN in the cell...
              d=uold(ind_cell(i),1)
              uold(ind_cell(i),1)=d-n*dstar
           endif
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids
  
  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=d*e
        end do
        ! Rick Sarmento - this is scaling the uold by density...
        ! See comment above - undoing what was done above
        ! This will adjust the scalars (Z, iprist, iturbvel, iprimordz)
        ! by the new cell density.
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

#endif

end subroutine star_formation 
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the input cell.
  !-----------------------------------------------------------------
  integer::nxny,i,idim,j,iok,ind
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:nvector,1:3),save::ix
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok

  if(ilevel==1)then
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do

  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do

  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do

  ! Get neighboring father grids  
  call getnborgrids(ind_grid_father,igridn,ncell)

  ! Loop over position
  do ind=1,twotondim

     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do

     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
              else
                 ind_father(i,j)=ind_cell(i)
              end if
           end if
        end do
     end do

  end do

end subroutine getnbor
