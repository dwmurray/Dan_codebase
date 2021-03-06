!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  ! DWM 05/2019 Stirring
  integer::i,j,id,iu,iv,iw,ip,iax,iay,iaz
  real(dp)::lambda,k,rho1,p1,v1!,b1,xx,yy,zz,theta!,expz,v2,xp,yp,zp,expp,v3 !Starting to try to remove vars
  real(dp),dimension(1:nvector,3)::acc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! DWM Creating a box.
  real(dp),dimension(1:ndim):: box_center, xFromCenter


  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........


  ! DWM 05/2019 Added in Stirring.
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call stir_acc_field(x,acc)
  
  id=1; iu=2; iv=3; iw=4; ip=5; iax=10; iay=11; iaz=12; ! Hardcoded for 3D!! DWM iax-iaz = 10, 11 ,12 b/c Rick mods
  rho1=d_region(1)
  v1=0d0
  p1=rho1*T2_star/(scale_v*scale_v)
  do i=1,nn
     q(i,id)=rho1

     q(i,iu)= 0d0
     q(i,iv)= 0d0
     q(i,iw)= 0d0

     q(i,ip)=p1
     
     ! add the stir acceleration field
     q(i,iax)=acc(i,1)
     q(i,iay)=acc(i,2)
     q(i,iaz)=acc(i,3)
     !DWM will begin by attempting a polluted cube of size 0.1*boxlen
     !Centered in the middle of the box.
     box_center(:)=0.5*boxlen !boxlen=1.0 (code units) !x is in cm
     xFromCenter(:) = abs(x(i,:)-box_center(:)) !just care about distance from center.
     if( xFromCenter(1) .le. 0.1*boxlen) then
        if( xFromCenter(2) .le. 0.1*boxlen) then
           if( xFromCenter(3) .le. 0.1*boxlen) then
              !make polluted fraction be 0.9 to check.
              q(i,7)=0.1 ! looking for ivar=pristine (ivar=7)
           end if
        end if
     end if
     !DWM TODO Test Zone for a sphere next.
  end do
!  if(verbose)write(*,*) 'Printing after add stir q(1,i):'
!  if(verbose)write(*,*) q(1,1),q(1,2),q(1,3),q(1,4),q(1,5),q(1,6),q(1,7),q(1,8),q(1,9),q(1,10),q(1,11),q(1,12)


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif

#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     ! Rick Sarmento - 18 Dec 2013
     ! Note that we are also tracking the pristine fraction
     ! as normalized by the density of the cell. We'll
     ! have to account for this.
     ! region_condinit init's the q's to 0.0 except for 
     ! iprist, which gets 1.0
     ! DWM 05/19 N.B. region_condinit (init_flow_fine.f90)
     ! loops through nvar, so if you add additional vars beyond RS
     ! modifications they will be set to 0.0
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

!  do ivar = 1,10
!     write(*,*) 'Q1',q(ivar,1),ivar
!     write(*,*) 'Q2',q(ivar,2),ivar
!     write(*,*) 'Q3',q(ivar,3),ivar
!     write(*,*) 'Q4',q(ivar,4),ivar
!     write(*,*) 'Q5',q(ivar,5),ivar
!     write(*,*) 'Q6',q(ivar,6),ivar
!     write(*,*) 'Q7',q(ivar,7),ivar
!     write(*,*) 'Q8',q(ivar,8),ivar
!     write(*,*) 'Q9',q(ivar,9),ivar
!     write(*,*) 'Q10',q(ivar,10),ivar
!     write(*,*) 'Q11',q(ivar,11),ivar
!     write(*,*) 'Q12',q(ivar,12),ivar
!  end do
!
!  do ivar = 1,10
!     write(*,*) 'U1',u(ivar,1),ivar
!     write(*,*) 'U2',u(ivar,2),ivar
!     write(*,*) 'U3',u(ivar,3),ivar
!     write(*,*) 'U4',u(ivar,4),ivar
!     write(*,*) 'U5',u(ivar,5),ivar
!     write(*,*) 'U6',u(ivar,6),ivar
!     write(*,*) 'U7',u(ivar,7),ivar
!     write(*,*) 'U8',u(ivar,8),ivar
!     write(*,*) 'U9',u(ivar,9),ivar
!     write(*,*) 'U10',u(ivar,10),ivar
!     write(*,*) 'U11',u(ivar,11),ivar
!     write(*,*) 'U12',u(ivar,12),ivar
!x  end do
!     write(*,*) 'U5',u(5,5)
end subroutine condinit

! DWM 05/2019 Stir_init_k_space
subroutine stir_initialize_k_space 
  use stir_parameters
  use ifport
  integer::i,j,k
  real(dp),parameter::twopi=6.2823d0
  real(dp)::total,norm,kvec
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !==============================================================
  !  This routine initialized the k vector space for stirring turbulence
  !==============================================================

  do i=1,nstir
     stir_kx(i) = 1D0*twopi*(i-1)/boxlen
     stir_ky(i) = 1D0*twopi*(i-1)/boxlen
     stir_kz(i) = 1D0*twopi*(i-1)/boxlen
  end do

  call srand(stir_seed)

  do i=1,nstir
     do j=1,nstir
        do k=1,nstir
           stir_amp_x(i,j,k) = rand()
           stir_amp_y(i,j,k) = rand()
           stir_amp_z(i,j,k) = rand()
           stir_phi_x(i,j,k) = rand()*twopi
           stir_phi_y(i,j,k) = rand()*twopi
           stir_phi_z(i,j,k) = rand()*twopi

           total = sqrt(stir_amp_x(i,j,k)**2 + stir_amp_y(i,j,k)**2 + stir_amp_z(i,j,k)**2)
           kvec = sqrt(stir_kx(i)**2 + stir_ky(j)**2 + stir_kz(k)**2)
           norm = 1d0*stir_norm*stir_norm/(scale_v*scale_v)*(kvec/(twopi/boxlen))**stir_index  
           if(kvec < twopi/boxlen*stir_kmin .or. kvec > twopi/boxlen*stir_kmax) norm = 0D0
           !if(kvec <= twopi/boxlen*2D0) write(*,*) kvec/(twopi/boxlen), norm

           stir_amp_x(i,j,k) = norm*stir_amp_x(i,j,k)/total
           stir_amp_y(i,j,k) = norm*stir_amp_y(i,j,k)/total
           stir_amp_z(i,j,k) = norm*stir_amp_z(i,j,k)/total

        end do
     end do
  end do

  stir_initialize = .true.
  return
end subroutine stir_initialize_k_space

! DWM 05/2019 subroutine stir_acc_field
subroutine stir_acc_field(x,acc)
  use amr_parameters, only : nvector,ndim
  use stir_parameters
  real(dp),dimension(1:nvector,1:ndim),intent(IN)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:3),intent(OUT):: acc ! acc
  !==============================================================
  !  This routine generates the initial stir field for the turbulence
  !==============================================================
  real(dp),dimension(1:nvector)::skx,ckx,sky,cky,skz,ckz,imtrigterms
  integer::i

  if( .not. stir_initialized) then
     if(verbose)write(*,*)'Initializing Stir k-space'
     call stir_initialize_k_space
  end if
  acc = 0D0
  do i = 1,nstir  
     do j=1,nstir
        do k=1,nstir
           skx(:) = sin( stir_kx(i)*x(:,1) + stir_phi_x(i,j,k))
           ckx(:) = cos( stir_kx(i)*x(:,1) + stir_phi_x(i,j,k))
           sky(:) = sin( stir_ky(j)*x(:,2) + stir_phi_y(i,j,k))
           cky(:) = cos( stir_ky(j)*x(:,2) + stir_phi_y(i,j,k))
           skz(:) = sin( stir_kz(k)*x(:,3) + stir_phi_z(i,j,k))
           ckz(:) = cos( stir_kz(k)*x(:,3) + stir_phi_z(i,j,k))

           imtrigterms = ckx*cky*skz + ckx*sky*ckz + skx*cky*ckz &
                - skx*sky*skz 
           ! create a div free stirring field
           acc(:,1) = acc(:,1) + (stir_ky(j)*stir_amp_z(i,j,k)-stir_kz(k)*stir_amp_y(i,j,k))*imtrigterms
           acc(:,2) = acc(:,2) + (stir_kz(k)*stir_amp_x(i,j,k)-stir_kx(i)*stir_amp_z(i,j,k))*imtrigterms
           acc(:,3) = acc(:,3) + (stir_kx(i)*stir_amp_y(i,j,k)-stir_ky(j)*stir_amp_x(i,j,k))*imtrigterms
           
           !acc(:,1) = acc(:,1) + (stir_kx(i)*stir_amp_x(i,j,k))*skr
        end do
     end do
  end do

  return
end subroutine stir_acc_field
  
