!Subroutines setting up and creating the stir accel field.

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
  !DWM My impression is I'll have to turn this flag off.
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
  

