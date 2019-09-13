subroutine stir_update
  use stir_parameters
  use hydro_parameters
!  use amr_commons
  implicit none
  integer::ilevel ! Refinement level
  integer::nn ! Number of Cells
  integer::ivar, iax, iay, iaz
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position
  real(dp),dimension(1:nvector,3)::acc  
  !###########################################################
  ! This subroutine is called every coarse step
  ! DWM - Aug 2019
  !###########################################################
  ! All of stirring ASSUMES! that it is the last three passive scalars.
  iax=nvar-stir_nvar+1; iay=iax+1; iaz=iay+1

  ! Stir_acc_field does check to see that acc has been init.
  call stir_acc_field(x,acc) !returns acc
! Next need to modify either uold or unew(iax,iay,iaz)

  ! See amr_commons for definition of 'active'
  ! active%igrid is array of ptrs to grids
  ! ------------------------------------------------------------------
  ! The outer loop, we look over all levels... 
  ! For each grid, get the indices of the active grids at that level

  ! Loop over levels
!  do ilevel=levelmin,nlevelmax ! Loop by amr level
!     ncache=active(ilevel)%ngrid ! active%ngrid field is number of grids - each grid can have 8 octs (cells)
!     do ind=1, twotondim ! Loop over the dimension
!        iskip=ncoarse+(ind-1)*ngridmax
!        do i=1,ncache ! Loop over the cells
!           ind_cell=active(ilevel)%igrid(i)+iskip !determine index of this cell
!           !DWM, atm just pulling values to double check.
!           d=max(unew(ind_cell,1),smallr)
!           u=0.0; v=0.0; w=0.0
!           ax=0.0; ay=0.0; az=0.0
!           if(ndim>0) then
!              u=unew(ind_cell,2)/d
!              ax=unew(ind_cell,iax)/d
!           end if
!           if(ndim>1) then
!              v=unew(ind_cell,3)/d
!              ay=unew(ind_cell,iay)/d
!           end if
!           if(ndim>2) then
!              w=unew(ind_cell,4)/d
!              az=unew(ind_cell,iaz)/d
!           end if
!           write(*,*) 'ax: ', ax, 'ay: ', ay, 'az: ', az
!!           if(ndim>0)then
!!              unew(ind_cell,iax)=d*acc(ind_cell,1)
!!           end if
!        end do
!     end do
!  end do
!
  if(verbose) write(*,*) 'Called to stir update'
end subroutine stir_update

!Subroutines setting up and creating the stir accel field.
! DWM 05/2019 Stir_init_k_space
!subroutine stir_initialize_k_space 
!  use stir_parameters
!  use ifport
!  integer::i,j,k
!  real(dp),parameter::twopi=6.2823d0
!  real(dp)::total,norm,kvec
!  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!
!  !==============================================================
!  ! This routine initializes the k vector space for stirring turbulence
!  ! Using a user defined stir_seed
!  !==============================================================
!  ! Conversion factor from user units to cgs units
!  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!
!  do i=1,nstir
!     stir_kx(i) = 1D0*twopi*(i-1)/boxlen
!     stir_ky(i) = 1D0*twopi*(i-1)/boxlen
!     stir_kz(i) = 1D0*twopi*(i-1)/boxlen
!  end do
!
!  call srand(stir_seed)
!
!  do i=1,nstir
!     do j=1,nstir
!        do k=1,nstir
!           stir_amp_x(i,j,k) = rand()
!           stir_amp_y(i,j,k) = rand()
!           stir_amp_z(i,j,k) = rand()
!           stir_phi_x(i,j,k) = rand()*twopi
!           stir_phi_y(i,j,k) = rand()*twopi
!           stir_phi_z(i,j,k) = rand()*twopi
!
!           total = sqrt(stir_amp_x(i,j,k)**2 + stir_amp_y(i,j,k)**2 + stir_amp_z(i,j,k)**2)
!           kvec = sqrt(stir_kx(i)**2 + stir_ky(j)**2 + stir_kz(k)**2)
!           norm = 1d0*stir_norm*stir_norm/(scale_v*scale_v)*(kvec/(twopi/boxlen))**stir_index  
!           if(kvec < twopi/boxlen*stir_kmin .or. kvec > twopi/boxlen*stir_kmax) norm = 0D0
!           !if(kvec <= twopi/boxlen*2D0) write(*,*) kvec/(twopi/boxlen), norm
!
!           stir_amp_x(i,j,k) = norm*stir_amp_x(i,j,k)/total
!           stir_amp_y(i,j,k) = norm*stir_amp_y(i,j,k)/total
!           stir_amp_z(i,j,k) = norm*stir_amp_z(i,j,k)/total
!
!        end do
!     end do
!  end do
!  stir_initialize = .true.
!  return
!end subroutine stir_initialize_k_space
!
subroutine stir_update_k_space(seed_value)
  use amr_parameters, only : nvector,ndim
  use stir_parameters
!  use ifport
  implicit none
  integer::i,j,k
  integer::rand_seed
  real(dp),parameter::twopi=6.2823d0
  real(dp)::total,norm,kvec
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::stir_initialize

  !==============================================================
  ! This routine uses a random seed to init the k vector space for stirring turbulence
  ! It returns
  !==============================================================
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do i=1,nstir
     stir_kx(i) = 1D0*twopi*(i-1)/boxlen
     stir_ky(i) = 1D0*twopi*(i-1)/boxlen
     stir_kz(i) = 1D0*twopi*(i-1)/boxlen
  end do
  write(*,*) 'Current Seed value: ', seed_value
  call srand(seed_value)
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
  write(*,*) 'stir_amp_x(1,1,1)', stir_amp_x(1,1,1)
  return
end subroutine stir_update_k_space


! DWM 05/2019 subroutine stir_acc_field
subroutine stir_acc_field(x,acc)
  use amr_parameters, only : nvector,ndim
  use stir_parameters
  real(dp),dimension(1:nvector,1:ndim),intent(IN)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:3),intent(OUT):: acc ! acc
  !==============================================================
  !  This routine generates the stir field for the turbulence
  !==============================================================
  real(dp),dimension(1:nvector)::skx,ckx,sky,cky,skz,ckz,imtrigterms
  integer::i

  if( .not. stir_initialized) then
     if(verbose)write(*,*)'Initializing Stir k-space'
     call stir_update_k_space(stir_seed)
     !call stir_initialize_k_space
  else
     write(*,*) 'Stir acc field is now updating.'
     rand_seed = rand()
     write(*,*) 'random new seed', rand_seed
     call stir_update_k_space(rand_seed)
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
