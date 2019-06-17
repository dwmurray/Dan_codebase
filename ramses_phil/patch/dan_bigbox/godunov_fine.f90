!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip
  real(dp)::d,u,v,w,e

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0.0
        end do
        do i=1,active(ilevel)%ngrid
           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           if(ndim>1)v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           if(ndim>2)w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           e=uold(active(ilevel)%igrid(i)+iskip,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(active(ilevel)%igrid(i)+iskip,ndim+2+irad)
           end do
#endif          
           enew(active(ilevel)%igrid(i)+iskip)=e
        end do
     end if
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
     if(pressure_fix)then
        do i=1,reception(icpu,ilevel)%ngrid
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end do
     end if
  end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew 
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx,fact,d_old

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Add gravity source terms to unew
  if(poisson)then
     call add_gravity_source_terms(ilevel)
  end if

  ! Add non conservative pdV terms to unew 
  ! for thermal and/or non-thermal energies
  if(pressure_fix.OR.nener>0)then
     call add_pdv_source_terms(ilevel)
  endif

  ! Add turbulent stirring terms
  if(stir)then
     call add_stir_source_terms(ilevel)
  end if

  ! Add jet terms
  if(jet_feedback) then
     call add_jet_source_terms(ilevel)
  end if

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        ! Correct total energy if internal energy is too small
        do i=1,active(ilevel)%ngrid
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(ind_cell,2)/d
           if(ndim>1)v=uold(ind_cell,3)/d
           if(ndim>2)w=uold(ind_cell,4)/d
           e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,ndim+2+irad)
           end do
#endif
           e_cons=uold(ind_cell,ndim+2)-e_kin
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,ndim+2)=e_prim+e_kin
           end if
        end do
     end if
  end do

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the gravity source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=unew(ind_cell,2)/d
        if(ndim>1)v=unew(ind_cell,3)/d
        if(ndim>2)w=unew(ind_cell,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+f(ind_cell,1)*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+f(ind_cell,2)*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+f(ind_cell,3)*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
  end do

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_stir_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use stir_parameters
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the stir source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact
  real(dp)::ax,ay,az
  integer::iax,iay,iaz

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add stirsource term at time t with half time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0.0; v=0.0; w=0.0
        ax=0.0; ay=0.0; az=0.0
        iax=nvar-stir_nvar+1; iay=iax+1; iaz=iay+1
        if(ndim>0) then
           u=unew(ind_cell,2)/d
           ax=unew(ind_cell,iax)/d
        end if
        if(ndim>1) then
           v=unew(ind_cell,3)/d
           ay=unew(ind_cell,iay)/d
        end if
        if(ndim>2) then
           w=unew(ind_cell,4)/d
           az=unew(ind_cell,iaz)/d
        end if
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+ax*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+ay*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+az*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
  end do

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_stir_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_jet_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use jet_parameters
  use pm_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the jet source terms
  ! with a full time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old
  integer :: ind_grid,ix,iy,iz,isink
  real(dp),dimension(1:twotondim,1:ndim) :: xc
  real(dp),dimension(1:ndim) :: skip_loc,boxcenter,z_axis,r_rel,xx,nr_rel
  real(dp)::scale,dx_loc,dx_min,r_jet,jetValuePerVol,dx,vol_loc
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_m,scale_Vol
  real(dp)::vjet,mjet
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! update only on the highest level
  if( ilevel .ne. nlevelmax) then
     return
  end if 

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m =scale_d*scale_l**ndim
  scale_Vol = scale_l**ndim

  ! Set position of cell centers relative to grid center

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  z_axis=(/0.0d0,0.0d0,1.0d0/)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  vol_loc=dx_loc**ndim
  dx_min=scale*0.5D0**nlevelmax*aexp
  boxcenter(:)=boxlen*0.5

  ! modify the sink particles at the highest level
  call updateProtostar( f_jet,  scale_m, scale_d, scale_t, dtnew(ilevel))
  ! Add stirsource term at time t with a full time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_grid=active(ilevel)%igrid(i)
        ind_cell=ind_grid+iskip
        xx(:)=xg(ind_grid,:)+xc(ind,:)
        xx(:)=(xx(:)-skip_loc(:))*scale
        do isink = 1, nsink

           r_rel(1:ndim) = xx(1:ndim)-xsink(isink,1:ndim)
        
           r_jet = sqrt(sum(r_rel*r_rel))
           call jetCalculation( r_rel, dx_loc, dx_min, lsink(isink,:), jet_cell_supersample, jetValuePerVol)
           if( jetValuePerVol .le. 0. ) then
              cycle ! go to the next sink
           end if

           nr_rel = r_rel/sqrt(sum(r_rel*r_rel))
              
           d=max(unew(ind_cell,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0) then
              u=unew(ind_cell,2)
           end if
           if(ndim>1) then
              v=unew(ind_cell,3)
           end if
           if(ndim>2) then
              w=unew(ind_cell,4)
           end if
           e_kin=0.5*(u**2+v**2+w**2)/d
           e_prim=unew(ind_cell,ndim+2)-e_kin
           d_old=max(uold(ind_cell,1),smallr)
           
           !vjet = 0.3*v_jet*sqrt(msink(isink)/1e33/scale_m)
           vjet = 0.3*v_jet*sqrt(msink(isink)/(1e33*scale_m)*(7e11*scale_d)/rsink(isink))
           mjet = f_jet*delta_mass(isink)  ! mjet is code units
          
           d = d + jetValuePerVol*scale_Vol*mjet
           unew(ind_cell,1)=d
           if(ndim>0)then
              u=u+jetValuePerVol*scale_Vol*mjet*vjet/scale_v*nr_rel(1)
              unew(ind_cell,2)=u
           endif
           if(ndim>1)then
              v=v+jetValuePerVol*scale_Vol*mjet*vjet/scale_v*nr_rel(2)
              unew(ind_cell,3)=v
           end if
           if(ndim>2)then
              w=w+jetValuePerVol*scale_Vol*mjet*vjet/scale_v*nr_rel(3)
              unew(ind_cell,4)=w
           endif
           e_kin=0.5*(u**2+v**2+w**2)/d
           unew(ind_cell,ndim+2)=e_prim+e_kin
        end do
     end do
  end do
  


111 format('   Entering add_jet_source_terms for level ',i2)

end subroutine add_jet_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine jetCalculation( center, dx, dx_min, lsink, ndivisions, jetValuePerVol)
  use jet_parameters      
  implicit none
  real(dp), dimension(3), intent(in) :: center, lsink
  real(dp), intent(in) :: dx, dx_min
  real(dp), intent(out) :: jetValuePerVol
  integer, intent(in) :: ndivisions
	real(dp) :: tempjetValuePerVol, rin, rout, r
  !------------------------------------------------------------------------
  ! This routine computes the jet deposition per cell depending on orientation and distance
  ! from sink particle -- relies on jetShape.
  !------------------------------------------------------------------------
  integer :: i, j, k
  real(dp), dimension(3):: origin, pos
	
  rin = dx_min*rin_jet
  rout = dx_min*rout_jet
  
  ! veto stuff that doesn't matter
  ! make sure it is between rin and rout
  r = sqrt(sum(center*center))
  
  if( r < rin .or. r > rout) then
     jetValuePerVol = 0.
     return 
  end if
  
  origin(1) = center(1) - dx/2.
  origin(2) = center(2) - dx/2.
  origin(3) = center(3) - dx/2.
  
  !Nested loop to divide a cell into n divisions
  !Getting the centroid value and averaging 
  jetValuePerVol = 0.
  DO i = 1, ndivisions
     pos(1) = origin(1) + i*dx/(ndivisions) - dx/(2*ndivisions)
     do j = 1, ndivisions
        pos(2) = origin(2) + j*dx/(ndivisions) - dx/(2*ndivisions)
        do k = 1, ndivisions
           pos(3) = origin(3) + k*dx/(ndivisions) - dx/(2*ndivisions)
           !call jetShape(pos, lsink, rin_jet*dx, rout_jet*dx, tempjetValuePerVol)
           call jetShape(pos, lsink, rin, rout, tempjetValuePerVol)
           jetValuePerVol = jetValuePerVol + tempjetValuePerVol
        END DO
     END DO
  END DO
  jetValuePerVol = jetValuePerVol/(ndivisions*ndivisions*ndivisions)
  return
end subroutine jetCalculation
!################################################################
!################################################################
!################################################################
!################################################################
subroutine jetShape( pos, lsink, rin, rout, jetValuePerVol) 
  use jet_parameters
  implicit none
  real(dp), dimension(3), intent(in) :: pos, lsink
  real(dp), intent(in) :: rin, rout
  real(dp), intent(out) :: jetValuePerVol
  !------------------------------------------------------------------------
  ! This routine computes the jet shape depending on orientation and distance
  ! from sink particle
  !------------------------------------------------------------------------
  real(dp) :: cos_theta, r, exponent, norm
  real(dp), dimension(3) :: npos, nlsink
  
  ! make sure it is between rin and rout
  r = sqrt(sum(pos*pos))
  
  if( r < rin .or. r > rout) then
     jetValuePerVol = 0.
     return 
  end if
  
  npos = pos/sqrt(sum(pos(1:3)*pos(1:3)))
  if( sum(pos*pos) <= 0d0 .or. sum(lsink*lsink) <= 0d0) then 
     jetValuePerVol = 0.
     return
  endif 

  nlsink = lsink/sqrt(sum(lsink(1:3)*lsink(1:3)))
  
  cos_theta = abs(sum(npos*nlsink)) ! this is a dot product
  
  cos_theta0_jet = cos(theta0_jet)
  
  exponent = (cos_theta - 1.)/(cos_theta0_jet -1.)
  if(exponent > 4.) then
     jetValuePerVol = 0.
     return
  end if
  
  norm = 2.*2.*3.14159*(1.-cos_theta0_jet)*(rout*rout*rout/3. - rin*rin*rin/3.) ! extra 2 is the 2 directions
  jetValuePerVol = exp(-exponent)/norm
  
  return
end subroutine jetShape
!################################################################
!################################################################
!################################################################
!################################################################
subroutine updateProtostar( f_jet, scale_m, scale_d, scale_t, dt)
  use pm_commons
  implicit none

  double precision, intent(IN) :: f_jet, scale_m, scale_d, scale_t, dt
  double precision :: morig, mnow, deltam, r_cgs, deltaR
  double precision :: L_star, R_star ! ZAMS values
  double precision :: beta, deltaT, mdot
  double precision, parameter :: f_acc = 1.0 ! This is a correction factor in Nakano 95.
  double precision, parameter :: a_e = 3./4. ! 3/4. is for n = 3
  double precision, parameter :: Grav = 6.7d-8
  integer :: isink

  deltaT = scale_t*dt
  do isink = 1, nsink
     morig = (msink(isink) - delta_mass(isink))*scale_m ! scaled to cgs
     deltaM = ((1.-f_jet)*delta_mass(isink))*scale_m ! scaled to cgs
     msink(isink) = msink(isink) - f_jet * delta_mass(isink)
     mnow = msink(isink)*scale_m
     call set_ZAMS( mnow, L_star, R_star, beta)

    !write(*,*) "called set_ZAMS", m+deltaM, L_star, R_star, beta
    mdot = deltaM/deltaT
 
    ! Nakano et al. 95 ApJ 450...183N Appendix A equation 27
    if( rsink(isink) <= 0.) then
        rsink(isink) = Grav*a_e*mnow*mdot/L_star*(1. + beta - (1.+f_acc)/a_e*0.5)/scale_d

    else 
        r_cgs = rsink(isink)*scale_d ! scaled to cgs
        deltaR = (2. - (1.+f_acc)/(2.*a_e))*r_cgs*deltaM/morig &
                 - r_cgs*r_cgs/(a_e*Grav*morig*morig)*L_star*deltaT
        rsink(isink) = rsink(isink) + deltaR/scale_d
    end if
 
    rsink(isink) = max( rsink(isink), R_star/scale_d)

 end do

 return

end subroutine updateProtostar
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine initialize_ZAMS
  use jet_parameters
  implicit none
  character(len=256) :: header
  double precision :: m, r, l, age
  integer :: line, ios
  open( unit=15, file="ZAMS_ezer67.txt", status="old")

  ! skip the first line
  read(15,*) header
  line = 0
  ios = 0
  ! read data
  do while( ios .eq. 0)
     read(15, *, iostat=ios) m, r, l, age
     if( ios .eq. 0) then
        line = line + 1
        lgm_zams(line) = log10(m)
        lgl_zams(line) = log10(l)
        lgr_zams(line) = log10(r)
     end if
  end do
  close(unit=15)

  num_zams = line

  ! figure out beta
  do line = 1, num_zams - 1
     beta_zams(line) = (lgl_zams(line+1)-lgl_zams(line))/(lgm_zams(line+1)-lgm_zams(line))
     betaR_zams(line) = (lgr_zams(line+1)-lgr_zams(line))/(lgm_zams(line+1)-lgm_zams(line))
  end do
  beta_zams(num_zams) = beta_zams(num_zams-1)
  betaR_zams(num_zams) = betaR_zams(num_zams-1)
  initialized_zams = .true.

  return
end subroutine initialize_ZAMS
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_ZAMS( M_star, L_star, R_star, beta)
  use jet_parameters
  implicit none
  double precision, intent(in) :: M_star
  double precision, intent(out) :: L_star, R_star, beta
  double precision :: lgm, lgL, lgR
  double precision, parameter :: m_sun = 2d33, L_sun = 4d33, R_sun=7d10
  integer :: i
  if( .not. initialized_zams ) then
     call initialize_zams
  end if

  lgm = log10(M_star/m_sun)
  beta = 1.

  if( lgm < lgm_zams(1)) then
     lgL = lgl_zams(1) + beta_zams(1)*(lgm - lgm_zams(1))
     lgR = lgr_zams(1) + betaR_zams(1)*(lgm - lgm_zams(1))
     beta = beta_zams(1)
  elseif( lgm >= lgm_zams(num_zams)) then
     lgL = lgl_zams(num_zams) + beta_zams(num_zams)*(lgm - lgm_zams(num_zams))
     lgR = lgr_zams(num_zams) + betaR_zams(num_zams)*(lgm - lgm_zams(num_zams))
     beta = beta_zams(num_zams)
  else
     do i = 1, num_zams-1
        if( lgm .ge. lgm_zams(i) .and. lgm < lgm_zams(i+1)) then
           lgL = lgl_zams(i) + beta_zams(i)*(lgm - lgm_zams(i))
           lgR = lgr_zams(i) + betaR_zams(i)*(lgm - lgm_zams(i))
           beta = beta_zams(i)
        end if
     end do
  end if

  L_star = 1d1**lgL*L_sun
  R_star = 1d1**lgR*R_sun

  return
end subroutine set_ZAMS
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell1
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp),dimension(1:nvector),save::divu_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
   
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather all neighboring velocities
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                 dx_g(i,idim) = dx_loc
              else
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                 dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                 dx_d(i,idim)=dx_loc
              else 
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                 dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
        end do
        ! End loop over dimensions
  
        ! Compute divu = Trace G
        divu_loc(1:ngrid)=0.0d0
        do i=1,ngrid
           do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
           enddo
        end do

        ! Update thermal internal energy 
        if(pressure_fix)then
           do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              eold=uold(ind_cell(i),ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_cell(i),ndim+2+irad)
              end do
#endif
              ! Add -pdV term
              enew(ind_cell(i))=enew(ind_cell(i)) &
                   & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
           end do
        end if

#if NENER>0
        do irad=1,nener
           do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),ndim+2+irad)=unew(ind_cell(i),ndim+2+irad) &
                & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),ndim+2+irad)*divu_loc(i)*dtnew(ilevel)
           end do
        end do
#endif

     enddo
     ! End loop over cells
  end do
  ! End loop over grids

  return

  ! This is the old technique based on the "pressure fix" option.

  ! Update thermal internal energy 
  if(pressure_fix)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           ind_cell1=active(ilevel)%igrid(i)+iskip
           ! Compute old thermal energy
           d=max(uold(ind_cell1,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(ind_cell1,2)/d
           if(ndim>1)v=uold(ind_cell1,3)/d
           if(ndim>2)w=uold(ind_cell1,4)/d
           eold=uold(ind_cell1,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              eold=eold-uold(ind_cell1,ndim+2+irad)
           end do
#endif
           ! Add pdV term
           enew(ind_cell1)=enew(ind_cell1) &
                & +(gamma-1.0d0)*eold*divu(ind_cell1) ! Note: here divu=-div.u*dt
        end do
     end do
  end if

#if NENER>0
  do irad=1,nener
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           ind_cell1=active(ilevel)%igrid(i)+iskip
           unew(ind_cell1,ndim+2+irad)=unew(ind_cell1,ndim+2+irad) &
                & +(gamma_rad(irad)-1.0d0)*uold(ind_cell1,ndim+2+irad)*divu(ind_cell1) ! Note: here divu=-div.u*dt
        end do
     end do
  end do
#endif

111 format('   Entering add_pdv_source_terms for level ',i2)

end subroutine add_pdv_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use stir_parameters
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0d0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0d0

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim
  integer::nvarskip

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     
     ! If not, interpolate hydro variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if
        
        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ncache)

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector

        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(stir) then
           do ivar=nvar-stir_nvar+1,nvar
              do i=1,ncache
                 ! convert to primitive and rescale
                 unew(ind_cell(i),ivar)=unew(ind_cell(i),1)*uold(ind_cell(i),ivar)/uold(ind_cell(i),1)
                 !write(*,*) "unew = ", unew(ind_cell(i), ivar), ivar
              end do
           end do
        end if
        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ncache
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ncache
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  !write(*,*) "ax = ", unew(ind_cell(10), 6), unew(ind_cell(10), 1), stir
  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     nvarskip=0
     if( stir) nvarskip=stir_nvar
     do ivar=1,nvar-nvarskip
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     !if(stir) then
     !   do ivar=nvar-stir_nvar+1,nvar
     !      ! convert to primitive and rescale
     !      unew(ind_buffer(i),ivar)=unew(ind_buffer(i),1)*uold(ind_buffer(i),ivar)/uold(ind_buffer(i),1)
     !   end do
     !end if

     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     nvarskip=0
     if(stir) nvarskip=stir_nvar

     do ivar=1,nvar-nvarskip
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     !if(stir) then
     !   do ivar=nvar-stir_nvar+1,nvar
     !      ! convert to primitive and rescale
     !      unew(ind_buffer(i),ivar)=unew(ind_buffer(i),1)*uold(ind_buffer(i),ivar)/uold(ind_buffer(i),1)
     !   end do
     !end if

     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

  end do
  ! End loop over dimensions

end subroutine godfine1
