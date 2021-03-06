!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  ! Sim scales not used in this routine... 
  !  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return ! Check number of grids at ilevel
  if(verbose)write(*,111)ilevel

  ! Gather star particles only.
#if NDIM==3
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0.and.tp(ipart).ne.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(idp(ipart).lt.0.and.tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

#endif

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons ! This is where h2_frac is defined
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine thermal_feedback. Each debris
  ! particle dumps mass, momentum and energy in the nearest grid cell
  ! using array uold.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,mstar,dx_min,vol_min, ctv
  real(dp)::xxx,mmm,t0,ESN,mejecta,zloss,z3loss
  real(dp)::ERAD,RAD_BOOST,tauIR,eta_sig,current_time
  real(dp)::sigma_d,delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,vol_loc,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,mz3loss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
!  integer ::idelay
!#ifdef SOLVERhydro
!  integer ::imetal=6
!#endif
!#ifdef SOLVERmhd
!  integer ::imetal=9
!#endif

! RS - handled in read_hyro_params.f90
!  if(metal)then
!     idelay=imetal+1
!  else
!     idelay=6
!  endif
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

  ! Massive star lifetime from Myr to code units
  t0=10.*1d6*(365.*24.*3600.)/scale_t
  current_time=t

  ! Minimum star particle mass
  mstar=n_star/(scale_nH*aexp**3)*vol_min

  ! Compute stochastic boost to account for target GMC mass
  SN_BOOST=MAX(mass_gmc*2d33/(scale_d*scale_l**3)/mstar,1d0)

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10.*2d33)/scale_v**2 ! RS - We'll assume E/mass is the same for regular and PISN 

  ! Life time radiation specific energy from cgs to code units
  ERAD=1d53/(10.*2d33)/scale_v**2

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sn2'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     if(ok(j))then
        if(levelp(ind_part(j))>=ilevel)then
           dteff(j)=dtnew(levelp(ind_part(j)))
        else
           dteff(j)=dtold(levelp(ind_part(j)))
        end if
     endif
  end do

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     if(ok(j))then
        mloss(j)=0d0
        mzloss(j)=0d0
        if (prist_gas_fraction) mz3loss(j)=0d0 ! RS - needed for iprimordz
        ethermal(j)=0d0
     endif
  end do

  ! Compute stellar mass loss and thermal feedback due to supernovae
  if(f_w==0)then ! Mean's we're not carrying any of the ISM gas along with the SNe - RS
     do j=1,np
        if(ok(j))then
           birth_time=tp(ind_part(j))
           ! Make sure that we don't count feedback twice
           if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j)))then
              ! Stellar mass loss
!!$              mejecta=eta_sn*mp(ind_part(j))
              ! Correct mass loss for two types of SN - PopIII and regular... 
              mejecta=eta_sn*mp(ind_part(j)) * (1.0-pfp(ind_part(j))) + eta_sn3 * mp(ind_part(j)) * pfp(ind_part(j))
              mloss(j)=mloss(j)+mejecta/vol_loc ! RS- Summing up mloss (in units of density)
              ! Thermal energy
              ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc
              ! Metallicity
              if(metal)then
                 ! yield - fraction of star particle that is transformed into Z
                 ! (1d0-yield)*zp(ind_part(j)) - part of the star particle that is "already" Z
                 zloss=yield+(1d0-yield)*zp(ind_part(j)) ! zloss is the total Z of ejecta
                 mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc ! RS - Sum up mzloss (in density units)
                 if(prist_gas_fraction) then
                    ! Rick Sarmento
                    ! We're tracking iprist and iprimordz: pristine fraction & primordial Z
                    ! For iprist, we'll let the increase in the cell's density be the factor
                    ! by which the pristine fraction is reduced. iprist tracks the cell's
                    ! pristine density, so when we update the cell's density with the
                    ! mass of the ejecta, and later use that density to convert iprist
                    ! to a mass-fraction, it will be scaled by the new density: rho_old/rho_new
                    ! where rho_new = rho_old + mzloss
!!$                    print *,"feedbk -- NO GAS CARRIED. SN Mass only"

                    ! For iprimordz, only count the fraction of metals created by primordial stars
                    ! So, there is only a yield of Z due to the SN... no old Z mass to carry along
                    z3loss = yield 
                    ! Sum up the Z mass due to POP III stars
                    mz3loss(j)=mz3loss(j)+mejecta*z3loss*pfp(ind_part(j))/vol_loc 
                    ! For testing - let the mz3loss exactly track mzloss ... ensure
                    ! that it is tracking exactly... 
!!$                    mz3loss(j)=mz3loss(j)+mejecta*zloss/vol_loc ! FOR TEST ONLY - should track mzloss
                 endif
              endif
              ! Reduce star particle mass
              mp(ind_part(j))=mp(ind_part(j))-mejecta
              idp(ind_part(j))=-idp(ind_part(j)) ! Reverse index to tag exploded stars
              ! Boost SNII energy and depopulate accordingly
              if(SN_BOOST>1d0)then
                 call ranf(localseed,RandNum)
                 if(RandNum<1d0/SN_BOOST)then
                    mloss(j)=SN_BOOST*mloss(j)
                    mzloss(j)=SN_BOOST*mzloss(j)
                    if(prist_gas_fraction) then
                       mz3loss(j)=SN_BOOST*mz3loss(j) ! Do I want to do this? RS
                       print *,"feedback: BOOSTING mz3loss. SN_BOOST",SN_BOOST
                    endif
                    ethermal(j)=SN_BOOST*ethermal(j)
                 else
                    mloss(j)=0d0
                    mzloss(j)=0d0
                    if(prist_gas_fraction) then
                       mz3loss(j)=0d0
                       print *,"feedbk: - reseting mz3loss = 0 because of RandNum!"
                    endif
                    ethermal(j)=0d0
                 endif
              endif
           endif
        end if
     end do ! end loop over all the star particles
  endif ! end if f_w = 0

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use the cell resolution to estimate the column density of gas
  delta_x=200*3d18
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN
  do j=1,np
     if(ok(j))then

        ! Infrared photon trapping boost
        if(metal)then
           tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
        else
           tauIR=tau_factor*max(uold(indp(j),1),smallr)
        endif
        RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
        
        ! Specific kinetic energy of the star
        ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
             &          +vp(ind_part(j),2)**2 &
             &          +vp(ind_part(j),3)**2)
        ! Update hydro variable in NGP cell
        ! Rick Sarmento - Note that we're updating the cell's
        ! density and other parms BEFORE making the pristine
        ! fraction computation (called in amr_step).
        ! So when iprist (which was scaled with the pre-SNe cell
        ! density) is updated, it'll be reduced because of the addition
        ! of the new mass to the cell. This is an approximation for the
        ! change in pristine fraction.

        ! Ensure the v_t isn't artificially descreased with the increase in rho
        if (turbulent_velocity) then
!!$           print *, "feedback v_t before rho update ", uold(indp(j),iturbvel)/uold(indp(j),1) * scale_v/1.0d5
           ctv = uold(indp(j),iturbvel)/uold(indp(j),1)
        end if
        uold(indp(j),1)=uold(indp(j),1)+mloss(j)
        if (turbulent_velocity) then
           uold(indp(j),iturbvel) = ctv * uold(indp(j),1) ! Using the updated rho
!!$           print *, "feedback v_t after rho update ", uold(indp(j),iturbvel)/uold(indp(j),1) * scale_v/1.0d5
        end if
        uold(indp(j),2)=uold(indp(j),2)+mloss(j)*vp(ind_part(j),1)
        uold(indp(j),3)=uold(indp(j),3)+mloss(j)*vp(ind_part(j),2)
        uold(indp(j),4)=uold(indp(j),4)+mloss(j)*vp(ind_part(j),3)
        uold(indp(j),5)=uold(indp(j),5)+mloss(j)*ekinetic(j)+ &
             & ethermal(j)*(1d0+RAD_BOOST)
     endif
  end do

  ! Add metals
  ! RS - this is only for the f_w=0 case -- no ISM gas carried along.
  ! I think this should be in the if block for f_w == 0.0 since
  ! that's the only way mz(3)loss can be updated! - ask Yohan
  ! Although, the mloss array is initialized to 0.0 so this won't
  ! hurt if it is executed
  if(metal)then
     do j=1,np
        if(ok(j))then
           uold(indp(j),imetal)=uold(indp(j),imetal)+mzloss(j) ! RS - Update the Z 'density' for the cell
           if(prist_gas_fraction) then
              ! Nothing to be done for iprist... Since iprist is based on
              ! the original cell's density, and that density has increase due
              ! to SN ejecta, the pristine fraction will naturally decrease when
              ! scaled by the cell's new rho. (rho updated above)
!!$              print *,"feedbk: f_w=0, imetal adjusted by mzloss:",mzloss(j)
!!$              print *,"feedbk: f_w=0, PGF before:",uold(indp(j),iprist)/(uold(indp(j),1)-mloss(j)) 
!!$              print *,"feedbk: f_w=0, PGF after :",uold(indp(j),iprist)/(uold(indp(j),1)) 
              ! iprimordz - Z from primordial (POPIII) stars
              ! mz3loss is density of 'Z' generated by POPIII (primordial) stars
              uold(indp(j),iprimordz) = uold(indp(j),iprimordz) + mz3loss(j) 
!!$              if (uold(indp(j),imetal) > 1d-10) print *,"feedbk: Z/(1-P) ",uold(indp(j),imetal)/(uold(indp(j),1)-uold(indp(j),iprist))
           endif
        endif
     end do
  endif

  ! Add delayed cooling switch variable
  if(delayed_cooling)then
     do j=1,np
        if(ok(j))then
           uold(indp(j),idelay)=uold(indp(j),idelay)+mloss(j)
        endif
     end do
  endif

#endif
  
end subroutine feedbk
!################################################################
!################################################################
! Called before thermal_feedback (above) in amr_step...
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all,PSN_all ! RS - arrays for particle pristine frac
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine checks SN events in cells where a
  ! star particle has been spawned.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_tot,info,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  logical ::ok_free
!#ifdef SOLVERhydro
!  integer ::imetal=6
!#endif
!#ifdef SOLVERmhd
!  integer ::imetal=9
!#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0,mass_load
  integer ,dimension(:),allocatable::indSN,iSN_myid
  ! Rick Sarmento - Adding PSN array tracking the part pristine fraction and Z3loadSN for gas xfer
  real(dp),dimension(:),allocatable::mSN,ZSN,PSN,vol_gas,ekBlast,mloadSN,ZloadSN, Z3loadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,vloadSN
  integer ,dimension(:),allocatable::indSN_tot,itemp
  real(dp),dimension(:),allocatable::mSN_tot,ZSN_tot,PSN_tot ! RS - arrays for particle prist frac 
  real(dp),dimension(:,:),allocatable::xSN_tot,vSN_tot
  integer::isort

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering kinetic_feedback'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Time delay for SN explosion from Myr to code units
  t0=t_delay*(1d6*365.*24.*3600.)/scale_t ! RS - t_delay is 10 Myr in amr_parameters.f90

  !------------------------------------------------------
  ! Gather star particles eligible for a SN event
  !------------------------------------------------------
  nSN_tot=0
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0        
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then ! RS - sim times counts up toward 0 from neg. Ensure we're t_delay past
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        nSN_tot=nSN_tot+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  enddo

  nSN_icpu=0
  nSN_icpu(myid)=nSN_tot  
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif
  nSN_tot=sum(nSN_icpu(1:ncpu))

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  if (nSN_tot .eq. 0) return
  
  ! Allocate arrays for the position and the mass of the SN
  ! Rick Sarmento - Adding PSN_tot to track the pristine fraction of the particle
  allocate(xSN_tot(1:nSN_tot,1:3),vSN_tot(1:nSN_tot,1:3),mSN_tot(1:nSN_tot),ZSN_tot(1:nSN_tot),PSN_tot(1:nSN_tot),itemp(1:nSN_tot))

  xSN_tot=0.;vSN_tot=0.;mSN_tot=0.;ZSN_tot=0.;PSN_tot=0

  !------------------------------------------------------
  ! Give position and mass of the star to the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid) ! RS - pointer to head of particle list
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart) ! Pointer to particle right after ipart
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then
                 iSN=iSN+1
                 xSN_tot(iSN,1)=xp(ipart,1)
                 xSN_tot(iSN,2)=xp(ipart,2)
                 xSN_tot(iSN,3)=xp(ipart,3)
                 vSN_tot(iSN,1)=vp(ipart,1)
                 vSN_tot(iSN,2)=vp(ipart,2)
                 vSN_tot(iSN,3)=vp(ipart,3)

                 if (metal) ZSN_tot(iSN)=zp(ipart)
                 ! Rick Sarmento -
                 ! We need to keep track of the particles' pristine fraction
                 ! When we do mass loading we'll have to adjust the primordial Z of the
                 ! star particle.
                 if (prist_gas_fraction) PSN_tot(iSN)=pfp(ipart) ! pfp: sp's pristing frac

                 ! **********************************************************
                 ! RS - 20 July 2016
                 ! Change mass yield based on pristine fraction of the sp.
                 ! For Pop III use eta_sn3 as the mass fraction returned to 
                 ! the ism.
                 ! **********************************************************
                 mSN_tot(iSN) = eta_sn*mp(ipart)*(1.0-PSN_tot(iSN)) + eta_sn3*mp(ipart)*PSN_tot(iSN)
!!$                 mSN_tot(iSN)  =eta_sn*mp(ipart)       ! Orig code... 
                 ! Remove the mass ejected by the SN
                 mp(ipart) =mp(ipart)-mSN_tot(iSN)
!!$                 mp(ipart) =mp(ipart)-eta_sn*mp(ipart) ! Orig code... 
                 ! Also compute energy, its dependent on SN type: 
                 ! NO... assume energy per unit mass in IMF is the same... (per discussion with Evan)
!!$                 ESN_tot(iSN) = ((1d51/(10d0*2d33))/scale_v**2)*(1.0-PSN_tot(iSN)) + & ! Regular SN
!!$                              &  10d0*(1d51/(10d0*2d33))/scale_v**2 * PSN_tot(iSN)   ! PISN SN, 10x Reg SN Energy
!!$                 print *,"old, new sn ejecta E",((1d51/(10d0*2d33))/scale_v**2),ESN_tot(iSN)
                 idp(ipart) =-idp(ipart) ! Negate the id, this particle has gone SN
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
     
        igrid=next(igrid)   ! Go to next grid
     end do
  end do 
  ! End loop over levels

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),ZSN_all(1:nSN_tot),PSN_all(1:nSN_tot)) 
  call MPI_ALLREDUCE(xSN_tot,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN_tot,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN_tot,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN_tot,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(PSN_tot,PSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info) ! Particle pristine fraction
!!$  call MPI_ALLREDUCE(ESN_tot,ESN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info) ! Particle SN energy
  ! RS - July 20 2016: added ESN for Sedov_Blast so it can use SN-type specific energy

  ! Copy back sums into local 'tot'
  xSN_tot=xSN_all
  vSN_tot=vSN_all
  mSN_tot=mSN_all
  ZSN_tot=ZSN_all
  PSN_tot=PSN_all
!!$  ESN_tot=ESN_all 
  deallocate(xSN_all,vSN_all,mSN_all,ZSN_all,PSN_all)
#endif

  call getSNonmyid(itemp,nSN,xSN_tot,nSN_tot)

  ! Allocate the arrays for the position and the mass of the SN
  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),PSN(1:nSN),iSN_myid(1:nSN)) ! RS - Adding PSN
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; PSN=0d0; iSN_myid=0

  do iSN=1,nSN
     isort=itemp(iSN)
     iSN_myid(iSN)=isort
     xSN(iSN,1)=xSN_tot(isort,1)
     xSN(iSN,2)=xSN_tot(isort,2)
     xSN(iSN,3)=xSN_tot(isort,3)
     vSN(iSN,1)=vSN_tot(isort,1)
     vSN(iSN,2)=vSN_tot(isort,2)
     vSN(iSN,3)=vSN_tot(isort,3)
     mSN(iSN)  =mSN_tot(isort)
     ZSN(iSN)  =ZSN_tot(isort) ! ZSN tracks the Z of the SNe at index iSN - RS
     PSN(iSN)  =PSN_tot(isort) ! Rick Sarmento - Adding an array with the particle's pristine fraction 
  enddo
  deallocate(xSN_tot,vSN_tot,mSN_tot,ZSN_tot,PSN_tot,itemp)

  allocate(vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN),indSN(1:nSN))
  allocate(mloadSN(1:nSN),ZloadSN(1:nSN),Z3loadSN(1:nSN),vloadSN(1:nSN,1:3))

  ! Compute the grid discretization effects
  ! Rick Sarmento - PSN is an array of the star particles' pristine fraction. 
  ! Z3loadSN generated by average_SN ... needed in Sedov_blast
  call average_SN(xSN,vSN,vol_gas,dq,ekBlast,indSN,nSN,nSN_tot,iSN_myid,mSN,mloadSN,ZSN,PSN,ZloadSN,Z3loadSN,vloadSN) 
  ! The above populates ZloadSN, need to add Z3loadSN
  
  deallocate(PSN)
  ! Modify hydro quantities to account for a Sedov blast wave
  ! Z3loadSN needed when blast mass xfer
  call Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,Z3loadSN,vloadSN) ! RS - Added Z3loadSN

  deallocate(xSN,vSN,mSN,ZSN,iSN_myid)
  deallocate(indSN,vol_gas,dq,ekBlast)
  deallocate(mloadSN,ZloadSN,Z3loadSN,vloadSN) ! Release Z3loadSN

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vSN,vol_gas,dq,ekBlast,ind_blast,nSN,nSN_tot,iSN_myid, &
     & mSN,mloadSN,ZSN,PSN,ZloadSN,Z3loadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  ! and do the mass loading process
  ! Called by kinetic_feedback
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,nSN_tot,j,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
!#ifdef SOLVERhydro
!  integer ::imetal=6
!#endif
!#ifdef SOLVERmhd
!  integer ::imetal=9
!#endif
  real(dp)::x,y,z,dr_SN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::eint,ekk,ekk1,ekk2,heat,mload,Zload,cpgf,ctv,Z3load 
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  ! RS - Adding PSN: the particles' pristine fraction, Z3loadSN for primordially generated Z gas xfer
  real(dp),dimension(1:nSN)::mSN,mloadSN,ZSN, PSN, ZloadSN, Z3loadSN,vol_gas,ekBlast 
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq,u2Blast,vloadSN
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN_tot)::vol_gas_mpi
  real(dp),dimension(1:nSN_tot)::vol_gas_all
  real(dp),dimension(1:nSN_tot)::mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all,Z3loadSN_mpi,Z3loadSN_all ! RS - adding Z3load...
  real(dp),dimension(1:nSN_tot,1:3)::dq_mpi,u2Blast_mpi,vloadSN_mpi
  real(dp),dimension(1:nSN_tot,1:3)::dq_all,u2Blast_all,vloadSN_all
#endif
  logical ,dimension(1:nvector),save::ok
  integer ,dimension(1:nSN)::iSN_myid
  integer::ind_SN

  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1;mloadSN=0.0;ZloadSN=0.0;Z3loadSN=0.0;vloadSN=0.0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel     ! ilevel 1->0.5, 2->0.25 ... dx is fractional box len
     dx_loc=dx*scale      ! Fraction of boxlen (which starts at 1 for ilevel 0)
     vol_loc=dx_loc**ndim ! So this is a volume as a fraction of total boxlen units
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 ! xc is cell center, xg grids position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Distance from cell center to SN: 
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2 ! Distance to SN squared - RS
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz)) ! Largest delta distance dim - RS
                    if(dr_SN.lt.rmax2)then ! Compare distance^2 to bubble radius^2... - RS
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc ! Add up gas for all cells in bubble
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax ! Distance as fraction of max radius of SN - RS
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then ! dx_loc is size of cell
                       ! We get in here if the SN went off in THIS cell
                       ind_blast(iSN)=ind_cell(i) ! This is the cell with the blast!
                       ekBlast  (iSN)=vol_loc     ! Not sure what's going on here??
                       d=uold(ind_blast(iSN),1)
                       if (d < smallr) print *,"ave_SN d < smallr but using d:", d, smallr
                       u=uold(ind_blast(iSN),2)/d
                       v=uold(ind_blast(iSN),3)/d
                       w=uold(ind_blast(iSN),4)/d
                       ekk=0.5d0*d*(u*u+v*v+w*w)       ! 0.5 rho v^2 - KE of gas in cell
                       eint=uold(ind_blast(iSN),5)-ekk ! Internal gas energy... 
                       ! Mass loading factor of the Sedov explosion
                       ! Ensure that no more that 25% of the gas content is removed
                       mload=min(f_w*mSN(iSN),0.25d0*d*vol_loc) ! f_w is mass loading factor for gas from cell
                       ! At this point mload is either 25% of the cells mass or f_w * mSN
                       mloadSN(iSN)=mSN(iSN)+mload ! mloadSN is the mass of the SN ejecta AND the gas carried along...
                       ! Should mloadSN be used to degrade the PGF immediately? Right here? Talk to Evan... 
!!$                       print *,"ave_SN: a = ",aexp
!!$                       print *,"ave_SN: z = ",1.0/aexp - 1.0
!!$                       print *,"ave_SN: cell index (SN): ",ind_blast(iSN)
!!$                       print *,"ave_SN: mSN (mass of ejecta only) - ",mSN(iSN)
!!$                       print *,"ave_SN: mloadSN (total mass of gas carried out) - ",mloadSN(iSN)
                       ! Update gas mass and metal content in the cell
                       if(metal)then
                          Zload=uold(ind_blast(iSN),imetal)/d ! Zload is blast cell's <Z>
                          ! mload*Zload                   - mass of metals in cell gas carried by SN
                          ! yield*(1.0-zSN(iSN))*mSN(iSN) - new metals created in the explosion
                          ! zSN(iSN)*mSN(iSN)             - mass of metals that existed in the star already
                          ZloadSN(iSN) = ( mload*Zload + yield*(1.0-zSN(iSN))*mSN(iSN) + zSN(iSN)*mSN(iSN)) / mloadSN(iSN)
                          ! ZloadSN is now the mass FRACTION of the ejecta that is metals
                          uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal)-Zload*mload/vol_loc
                          ! The above accounts for the metals in ZloadSN LEAVING the cell...
!!$                          print *,"ave_SN: total ejecta mass (SN mass + gas): ", mloadSN(iSN)
!!$                          print *,"ave_SN: total cell gas mass lost: ", mload
!!$                          print *,"ave_SN: ZloadSN : ",ZloadSN(iSN)
                          ! Rick Sarmento
                          ! However, I'm reducing the primordial Z of the gas based on a mass loading factor
                          ! and blast since some of the primordial metals are going to be carried into
                          ! neghboring cells. 
                          ! Primordial Z is the mass of metals created by POP III stars. 
                          if (prist_gas_fraction) then
                             ! Pristine fraction: iprist - the pristine fraction cell density
                             ! Mass (density) is leaving the cell. We need to adjust the cell's iprist
                             ! density based on the mass of stuff leaving the cell. This ensures the
                             ! pristine fraction for the cell doesn't go UP, but remains the same.
                             cpgf = uold(ind_blast(iSN),iprist)/d ! Recover the pure pristine gas fraction
                             if (cpgf > 1.0d0) cpgf = 1.0d0       ! Handle any numerical rounding issues
!!$                             print *,"ave_SN: cell starting iprist: ", uold(ind_blast(iSN),iprist)
!!$                             print *,"ave_SN: cell starting PGF: ", cpgf
                             uold(ind_blast(iSN),iprist)=cpgf * (d - mload/vol_loc) ! iprist with new cell rho

                             ! Primordial Metals
                             ! Update iprimordZ due to primordial metals lost
                             Z3load=uold(ind_blast(iSN),iprimordz)/d ! Z3load is the cell's primordial Z fraction
                             ! Compute Z3loadSN that we'll use to inject primoridal Z into other cells.
                             ! mload*Z3load            - mass of primordial metals in cell gas carried by SN
                             ! yield*mSN(iSN)*PSN(iSn) - mass of metals created by primordial stars
                             ! Compute the mass of the yield in metals from the pristine stars ONLY
                             Z3loadSN(iSN)=( mload*Z3load + yield*mSN(iSN)*PSN(iSN)) / mloadSN(iSN)

                             ! TEST ONLY
                             ! By making Z3loadSN the same as Z3load we can check the primordial Z
                             ! directly against metallicity (imetal) ... they should track exactly... 
!!$               Z3loadSN(iSN)=( mload*Z3load + yield*(1.0-zSN(iSN))*mSN(iSN) + zSN(iSN)*mSN(iSN)) / mloadSN(iSN)

                             ! Remove the mass-density of primordial-Z gas from this cell
                             uold(ind_blast(iSN),iprimordz)=uold(ind_blast(iSN),iprimordz)-(Z3load*mload/vol_loc)
                          endif
                          if (turbulent_velocity) then
                             ! Need to make sure v_t doesn't go down because of rho leaving the cell...
                             ctv = uold(ind_blast(iSN),iturbvel)/d
!!$                             print *,"ave_SN: orig v_t (km/s)", ctv * scale_v/1.0d5
                             uold(ind_blast(iSN),iturbvel) = ctv * (d - mload/vol_loc) ! iprist with new cell rho
                          end if
                       endif
!!                       print *,"ave_SN: cell rho start:", d
                       d=uold(ind_blast(iSN),1)-mload/vol_loc ! Very important! Updates the cell density - RS
!!                       print *,"ave_SN: cell rho end  :", d, ind_blast(iSN)
                       uold(ind_blast(iSN),1)=d
                       uold(ind_blast(iSN),2)=d*u
                       uold(ind_blast(iSN),3)=d*v
                       uold(ind_blast(iSN),4)=d*w
                       uold(ind_blast(iSN),5)=eint+0.5d0*d*(u*u+v*v+w*w)

                       vloadSN(iSN,1)=(mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(iSN)
                       vloadSN(iSN,2)=(mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(iSN)
                       vloadSN(iSN,3)=(mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  !################################################################
#ifndef WITHOUTMPI
  vol_gas_mpi=0d0; dq_mpi=0d0; u2Blast_mpi=0d0; mloadSN_mpi=0d0; ZloadSN_mpi=0d0; Z3loadSN_mpi=0d0; vloadSN_mpi=0d0
  ! Put the nSN size arrays into nSN_tot size arrays to synchronize processors
  do iSN=1,nSN
     ind_SN=iSN_myid(iSN)
     vol_gas_mpi(ind_SN)=vol_gas(iSN)
     mloadSN_mpi(ind_SN)=mloadSN(iSN)
     ZloadSN_mpi(ind_SN)=ZloadSN(iSN)
     Z3loadSN_mpi(ind_SN)=Z3loadSN(iSN) ! - RS
     vloadSN_mpi(ind_SN,1)=vloadSN(iSN,1)
     vloadSN_mpi(ind_SN,2)=vloadSN(iSN,2)
     vloadSN_mpi(ind_SN,3)=vloadSN(iSN,3)
     dq_mpi     (ind_SN,1)=dq     (iSN,1)
     dq_mpi     (ind_SN,2)=dq     (iSN,2)
     dq_mpi     (ind_SN,3)=dq     (iSN,3)
     u2Blast_mpi(ind_SN,1)=u2Blast(iSN,1)
     u2Blast_mpi(ind_SN,2)=u2Blast(iSN,2)
     u2Blast_mpi(ind_SN,3)=u2Blast(iSN,3)
  enddo
  call MPI_ALLREDUCE(vol_gas_mpi,vol_gas_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mloadSN_mpi,mloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZloadSN_mpi,ZloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Z3loadSN_mpi,Z3loadSN_all,nSN_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info) ! - RS
  call MPI_ALLREDUCE(vloadSN_mpi,vloadSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq_mpi     ,dq_all     ,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast_mpi,u2Blast_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_mpi=vol_gas_all
  mloadSN_mpi=mloadSN_all
  ZloadSN_mpi=ZloadSN_all
  Z3loadSN_mpi=Z3loadSN_all ! - RS
  vloadSN_mpi=vloadSN_all
  dq_mpi     =dq_all
  u2Blast_mpi=u2Blast_all
  ! Put the nSN_tot size arrays into nSN size arrays
  do iSN=1,nSN
     ind_SN=iSN_myid(iSN)
     vol_gas(iSN)=vol_gas_mpi(ind_SN)
     mloadSN(iSN)=mloadSN_mpi(ind_SN)
     ZloadSN(iSN)=ZloadSN_mpi(ind_SN)
     Z3loadSN(iSN)=Z3loadSN_mpi(ind_SN) ! - RS
     vloadSN(iSN,1)=vloadSN_mpi(ind_SN,1)
     vloadSN(iSN,2)=vloadSN_mpi(ind_SN,2)
     vloadSN(iSN,3)=vloadSN_mpi(ind_SN,3)
     dq     (iSN,1)=dq_mpi     (ind_SN,1)
     dq     (iSN,2)=dq_mpi     (ind_SN,2)
     dq     (iSN,3)=dq_mpi     (ind_SN,3)
     u2Blast(iSN,1)=u2Blast_mpi(ind_SN,1)
     u2Blast(iSN,2)=u2Blast_mpi(ind_SN,2)
     u2Blast(iSN,3)=u2Blast_mpi(ind_SN,3)
  enddo
#endif
  !################################################################
  ! At this point ekBlast is vol_loc
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
        if (ekBlast(iSN).le.1d-5) ekBlast(iSN)=0d0
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,Z3loadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! 
  !------------------------------------------------------------------------
  integer::ilevel,j,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
!#ifdef SOLVERhydro
!  integer ::imetal=6
!#endif
!#ifdef SOLVERmhd
!  integer ::imetal=9
!#endif
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,d_gas,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax, ctv
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  ! RS - Adding the ave Z3loadSN mass loading for primordialy generated Z
  real(dp),dimension(1:nSN)::mSN,p_gas,vol_gas,uSedov,ekBlast,mloadSN,ZloadSN,Z3loadSN
  real(dp),dimension(1:nSN,1:3)::xSN,dq,vloadSN
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l ! rmax = rbubble converted to cm converted to code units
  rmax2=rmax*rmax
  
  ! Ejecta specific energy (accounting for dilution)
  ! RS - Need a new array passed into Sedov to account for differing
  ! energy for the different types of SN... 
  ESN=(1d51/(10d0*2d33))/scale_v**2

  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then       ! vol_gas is the total vol of gas in all cells in the blast radius
        d_gas=mSN(iSN)/vol_gas(iSN) ! d_gas is fraction of SN eject mass to ALL cells in the blast radius
        if(ekBlast(iSN)==0d0)then
           p_gas (iSN)=d_gas*ESN    ! RS - here's where we'd use a new ESN array - specific to SN type
           uSedov(iSN)=0d0
        else
           p_gas (iSN)=(1d0-f_ek)*d_gas*ESN
           uSedov(iSN)=sqrt(f_ek*2.0*mSN(iSN)*ESN/ekBlast(iSN)/mloadSN(iSN))
        endif
     else
        p_gas(iSN)=mSN(iSN)*ESN/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then ! If we're within the max radius of the blast/bubble... - RS
                       ! vol_gas(iSN) is the total gas volume within the blast radius
                       ! mloadSN is mass of SN and gas carried along from cell
                       d_gas=mloadSN(iSN)/vol_gas(iSN) ! d_gas is density fraction of mloadSN going into a single cell
!!$                       print *,"Sedov_blast: a = ",aexp
!!$                       print *,"Sedov_blast: z = ",1.0/aexp - 1.0
!!$                       print *,"Sedov_blast: cell index ",ind_cell(i)
!!$                       print *,"Sedov_blast: cell starting rho ",uold(ind_cell(i),1)
!!$                       print *,"Sedov_blast: rho incoming (d_gas) ",d_gas

                       ! So we're dividing up the mloadSN equally into the cell's within the blast radius
                       ! Compute the density and the metal density of the cell

                       if (turbulent_velocity) then
                          ctv = uold(ind_cell(i),iturbvel)/uold(ind_cell(i),1) ! current v_t
!!$                          print *,"Sedov_blast: v_t (km/s) before d inc ", ctv * scale_v/1.0d5
                       end if

                       ! ------------------
                       ! NEW GAS ADDED HERE
                       ! ------------------
                       uold(ind_cell(i),1) = uold(ind_cell(i),1)+d_gas 

                       if (turbulent_velocity) then
                          uold(ind_cell(i),iturbvel) = ctv * uold(ind_cell(i),1)
!!$                          print *,"Sedov_blast: v_t (km/s) after inc fixed ", uold(ind_cell(i),iturbvel)/uold(ind_cell(i),1) * scale_v/1.0d5
                       end if

                       ! Update metallicity
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_gas*ZloadSN(iSN)

                       if(prist_gas_fraction)then
                          ! Pristine fraction
                          ! Since the cell's density goes up due to the injection of mass,
                          ! some of which is new Z, the pristine fraction (iprist/rho) will 
                          ! go down since rho has increased. Remember, iprist is a mass density.
                          ! So, nothing to do for iprist
!!$                          print *,"Sedov_blast: updated Z ",uold(ind_cell(i),imetal)/uold(ind_cell(i),1)
!!$                          print *,"Sedov_blast: updated P ",uold(ind_cell(i),iprist)/uold(ind_cell(i),1)

!!$                          uold(ind_cell(i),iprist) = uold(ind_cell(i),iprist) - d_gas*(1.0d0-ZloadSN(iSN))

                          ! Update the primordial Z mass density in the cell due to injected
                          ! metals from POP III stars. (Z3loadSN is computed above...)
                          uold(ind_cell(i),iprimordz)=uold(ind_cell(i),iprimordz)+d_gas*Z3loadSN(iSN)
                       endif
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vloadSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vloadSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vloadSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        d_gas=mloadSN(iSN)/ekBlast(iSN)
        u=vloadSN(iSN,1)
        v=vloadSN(iSN,2)
        w=vloadSN(iSN,3)
        if(indSN(iSN)>0)then
           print *,"Sedov_blast: NO gas in SN blast radius!, rho ratio old/new ",uold(indSN(iSN),1)/(uold(indSN(iSN),1)+d_gas)
           if (turbulent_velocity) then
              ctv = uold(indSN(iSN),iturbvel)/uold(indSN(iSN),1)
           end if
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_gas*ZloadSN(iSN)
           ! Rick Sarmento - TODO
           ! We're looping over the SNe... adding the mass loading to the cells
           ! THIS IS ONLY DONE if the cell's gas volume is 0.0 ... 
           if(prist_gas_fraction) then
              ! Pristine fraction
              ! Since the cell's density goes up due to the injection of mass,
              ! some of which is new Z, the pristine fraction (iprist/rho) will 
              ! go down since rho has increased. Remember, iprist is a mass density.
              ! So, nothing to do for iprist
              ! print *,"Sedov_blast: NOG d_gas:", d_gas
              ! print *,"Sedov_blast: NOG PGF:", uold(indSN(iSN),iprist)/uold(indSN(iSN),1)
              ! print *,"Sedov_blast: NOG Z/(1-PGF) new:", (uold(indSN(iSN),imetal))/(uold(indSN(iSN),1)-uold(indSN(iSN),iprist))
              ! Primordial Z
              ! Update the primordial Z mass density in the cell due to injected
              ! metals from POP III stars. (Z3loadSN is computed above...)
              uold(indSN(iSN),iprimordz)=uold(indSN(iSN),iprimordz)+d_gas*Z3loadSN(iSN) 
!!$              print *,"SB, vol_gas = 0, d_gas    =",d_gas
!!$              print *,"SB, vol_gas = 0, Z3loadSN =",Z3loadSN(iSN)
!!$              print *,"SB, vol_gas = 0, new iprimordz =",uold(indSN(iSN),iprimordz)
!!$              print *,"SB, vol_gas = 0, imetal        =",uold(indSN(iSN),imetal)
           endif
           if (turbulent_velocity) then
              uold(indSN(iSN),iturbvel) = ctv * uold(indSN(iSN),1)              
           end if
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getSNonmyid(iSN_myid,nSN_myid,xSN,nSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid
  integer::nSN_myid,ii,nSN
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,iSN,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(dp),dimension(1:nSN,1:3)::xSN
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp),dimension(1:1)::order_min ! R. Sarmento - hilbert3d expects a double prec array for order min.
  !  real(qdp)::dkey,order_min,oneqdp=1.0
  real(qdp)::dkey,oneqdp=1.0
  
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
  do iSN=1,nSN

     cpu_read=.false.
     ! Compute boundaries for the SN cube of influence
     xxmin=(xSN(iSN,1)-drSN)/scale ; xxmax=(xSN(iSN,1)+drSN)/scale
     yymin=(xSN(iSN,2)-drSN)/scale ; yymax=(xSN(iSN,2)+drSN)/scale
     zzmin=(xSN(iSN,3)-drSN)/scale ; zzmax=(xSN(iSN,3)+drSN)/scale

     if(TRIM(ordering).eq.'hilbert')then
        
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif
        
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min(1)=0.0d0
           endif
           bounding_min(i)=(order_min(1))*dkey
           bounding_max(i)=(order_min(1)+oneqdp)*dkey
        end do
        
        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if
     
     ! Create the index array for SN in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           iSN_myid(ii)=iSN
        endif
     enddo

  enddo
  
  ! Number of SN in processor myid
  nSN_myid=ii

end subroutine getSNonmyid
!################################################################
!################################################################
!################################################################
!################################################################
