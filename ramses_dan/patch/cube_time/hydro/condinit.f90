!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use stir_parameters
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
  ! DWM 05/19 N.B. region_condinit (init_flow_fine.f90)
  ! loops through nvar, so if you add additional vars beyond RS
  ! modifications they will be set to 0.0
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........


  ! DWM 05/2019 Added in Stirring.
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call stir_acc_field(x,acc)
  
  id=1; iu=2; iv=3; iw=4; ip=5; iax=nvar-stir_nvar+1; iay=iax+1; iaz=iay+1; ! DWM iax-iaz = 10, 11 ,12 for Ricks code
  rho1=d_region(1)
  v1=0d0
  p1=rho1*T2_star/(scale_v*scale_v)
  do i=1,nn
     !if(myid==1)write(*,*) 'Condinit x(i,:): ', x(i,:)
     q(i,id)=rho1

     q(i,iu)= 0d0
     q(i,iv)= 0d0
     q(i,iw)= 0d0

     q(i,ip)=p1
     
     ! add the stir acceleration field
     q(i,iax)=acc(i,1)
     q(i,iay)=acc(i,2)
     q(i,iaz)=acc(i,3)
  end do
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
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif
!     write(*,*) 'U5',u(5,5)
end subroutine condinit
