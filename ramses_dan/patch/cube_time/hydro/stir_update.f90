subroutine stir_update
  use stir_parameters
  use hydro_parameters
  use ifport ! use for rand()
  use hydro_commons, ONLY: uold ! We mod the last three vars in uold.
  implicit none
  integer::ilevel, ivar, acc_var
  integer::iax, iay, iaz !accelerations
  integer::i,igrid,ncache,iskip,ngrid !ngrid is number of grids
  integer::ind,idim,ix,iy,iz,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,scale,dx_loc
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position
  real(dp),dimension(1:nvector,3)::acc  
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer:: rand_stir_seed
  rand_stir_seed = 0

  !###########################################################
  ! This subroutine is called by amr/amr_step every coarse step
  ! DWM - Aug 2019
  !###########################################################

  !Check if we need to update at this time.
  if (t<stir_timescale(stir_tout))return

  if(myid==1) write(*,89)t,stir_timescale(stir_tout) !visual check
89 format('t=',1pe12.5,' stir_timescale=',1pe12.5)

  stir_tout = stir_tout + 1 !
  !Condition to turn off stirring without having to restart run.
  !Also used to trigger if we want polluted shapes.
  if(stir_tout > 10) then !Change to user set value, not hardcoded !TODO DWM
     stir = .false.
     if(myid==1) then
        write(*,*) 'System has evolved past user requested time of stirring.'
        write(*,*) 'Stir flag is now false: ', stir
     end if
     ! DWM Implementing Polluted Shapes
     if( .not. pollute_initialized) then
        call dump_all
        if(myid==1) write(*,*) 'Now Implementing the polluted shapes.'
        call polluted_shapes
        call dump_all
     end if
     !End Polluted block
     !With stirring now turned off, we simply return.
     return
  end if

  !Randomize the K space
  if( .not. stir_initialized) then
     if(myid==1) write(*,*)'Initializing Stir k-space'
     call stir_update_k_space(stir_seed) !stir_seed is provided by stir_parameters or user.
  else
     ! rand() returns a value between 0 - 1,
     rand_stir_seed = INT(rand() * 1d6)! Convert to int and make mag same as user supplied.
     if(myid==1) write(*,*) 'Stir Seed, Rand_seed: ', stir_seed, rand_stir_seed
     call stir_update_k_space(rand_stir_seed)
  end if

  ! All of stirring ASSUMES! that it is the last three passive scalars.
  iax=nvar-stir_nvar+1; iay=iax+1; iaz=iay+1
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2) 

  do ilevel=1,nlevelmax !Loop through all levels
     ! Mesh size at level ilevel in coarse cell units
     dx=0.5D0**ilevel

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
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     ncache=active(ilevel)%ngrid

     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 x(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 x(i,idim)=(x(i,idim)-skip_loc(idim))*scale
              end do
           end do
           !Now that we have the cell centres in cgs
           call stir_acc_field(x,acc) !returns acc
           ! Now implement the updated acc variables
           do i=1,ngrid
              do ivar=nvar-stir_nvar+1, nvar ! This loops 10 - 12 for uold.
                 acc_var = 0 !Pick the correct acc index for 
                 if(ivar==iax) acc_var = 1
                 if(ivar==iay) acc_var = 2
                 if(ivar==iaz) acc_var = 3
                 if(acc_var == 0)then
                    write(*,*) 'acc_var is zero, are we stirring?'
                    call clean_stop
                 end if
                 ! Recall uold passes everything around weighted by density
                 uold(ind_cell(i),ivar) = acc(i,acc_var) * uold(ind_cell(i),1) 
              end do
           end do
        end do
     end do
  end do





111 format('   Entering Stir Update for level ',I2)

end subroutine stir_update

subroutine stir_update_k_space(seed_value)
  use amr_parameters, only : nvector,ndim
  use stir_parameters
  use ifport
  implicit none
  integer::i,j,k
  integer::seed_value
  real(dp),parameter::twopi=6.2823d0
  real(dp)::total,norm,kvec
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::stir_initialize

  !==============================================================
  ! This routine uses a random seed to init or update the k vector space for stirring turbulence
  ! It returns stir_amp and stir_k 
  !==============================================================
  if(verbose) write(*,*) 'Current Seed value: ', seed_value
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do i=1,nstir
     stir_kx(i) = 1D0*twopi*(i-1)/boxlen
     stir_ky(i) = 1D0*twopi*(i-1)/boxlen
     stir_kz(i) = 1D0*twopi*(i-1)/boxlen
  end do

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
  stir_initialized = .true.
  return
end subroutine stir_update_k_space

! DWM 05/2019 subroutine stir_acc_field
subroutine stir_acc_field(x,acc)
  use amr_parameters, only : nvector,ndim
  use stir_parameters
  implicit none
  real(dp),dimension(1:nvector,1:ndim),intent(IN)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:3),intent(OUT):: acc ! acc
  !==============================================================
  !  This routine generates the stir field for the turbulence
  !  It also either initializes the field according to a user specified
  !  seed value, or updates the field via a random seed.
  !==============================================================
  real(dp),dimension(1:nvector)::skx,ckx,sky,cky,skz,ckz,imtrigterms
  integer::i,j,k

  if( .not. stir_initialized) then
     if(myid==1) write(*,*)'Initializing Stir k-space field'
     call stir_update_k_space(stir_seed) !stir_seed is provided by stir_parameters or user.
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
