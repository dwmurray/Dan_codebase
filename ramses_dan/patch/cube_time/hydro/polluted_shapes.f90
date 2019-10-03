!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================
subroutine polluted_shapes
  use hydro_parameters
  use stir_parameters
  implicit none
  !=========================================================================
  !This subroutine is for the set up of various shapes of polluted material.
  !=========================================================================
  integer :: ilevel
  integer :: ngrid !number of cells
  integer :: i,ivar,igrid,ncache,iskip,idim,ind
  integer :: ix,iy,iz,nx_loc
  real(dp):: scale,dx
  integer, dimension(1:nvector)::ind_grid,ind_cell
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx ! Cell center position.
  !id=1; iu=2; iv=3; iw=4; ip=5; imetal=6;iPrist =7;v_turb=8;primordz=9;
  !iax=nvar-stir_nvar+1; iay=iax+1; iaz=iay+1; ! DWM iax-iaz = 10, 11 ,12 for Ricks code
  !DWM will begin by attempting a polluted cube of size 0.1*boxlen
  !Centered in the middle of the box.

  do ilevel=1,nlevelmax
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
     !dx_loc=dx*scale
     ncache=active(ilevel)%ngrid

     if(myid==1) write(*,*)'ilevel: ', ilevel, 'ncache: ', ncache, 'nvector: ', nvector
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(myid==1) write(*,*)'Polluted init - ngrid: ', ngrid
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
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                 !if(myid==1)write(*,*) 'xx(i,:): ', xx(i,:)
              end do
           end do
           !DWM where the call to condinit was
           !if(myid==1)write(*,*) 'Have found cell centers for this ilevel: ', ilevel
           call polluted_cube(xx,ngrid)
        end do
     end do
  end do
  !DWM TODO Test Zone for a sphere next.  
  !if(myid==1)write(*,*) 'Finished looping to create the cube.'
  pollute_initialized = .True.
  return
end subroutine polluted_shapes

subroutine polluted_cube(x,nn)
  use amr_commons, ONLY: myid
  use hydro_parameters
  use hydro_commons, ONLY: uold
  !use stir_parameters
  implicit none
  integer :: nn !number of cells
  integer :: i, ivar
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp),dimension(1:ndim):: box_center, xFromCenter

  do i=1,nn
     box_center(:)=0.5*boxlen !boxlen is in cgs.!x is in cm
     !if(myid==1)write(*,*) 'boxlen: ', boxlen, 'box_center(:)', box_center(:)
     !if(myid==1)write(*,*) 'x(i,:): ', x(i,:)!, 'box_center(:)', box_center(:)
     xFromCenter(:) = abs(x(i,:)-box_center(:)) !just care about distance from center.
     !if(myid==1)write(*,*) 'xFromCenter(:): ', xFromCenter(:)
     if( xFromCenter(1) .le. 0.1*boxlen) then
        if( xFromCenter(2) .le. 0.1*boxlen) then
           if( xFromCenter(3) .le. 0.1*boxlen) then
              !make polluted fraction be 0.9 to check.
              !q(i,7)=0.1 ! looking for ivar=pristine (ivar=7)
              do ivar=ndim+3, nvar
                 !if(myid==1)write(*,*)'ivar: ', ivar, 'iprist: ', iprist
                 if (ivar==iprist)then
                    !if(myid==1)write(*,*) 'Pollution setting uold.'
                    uold(1:nn,ivar)=0.1d0 * uold(1:nn, 1) !recall everything is in density units
                 end if
              end do
           end if
        end if
     end if
  end do
  !End of Polluted shaped block.
  return
end subroutine polluted_cube
