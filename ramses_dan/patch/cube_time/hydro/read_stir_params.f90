subroutine read_stir_params(nml_ok)
  use stir_parameters
  implicit none
  logical::nml_ok
  integer:: i !DWM setup stir timescale
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/stir_params/stir_seed,stir_norm,stir_index,& 
        stir_kmin,stir_kmax,stir_delta_t
  rewind(1)
  read(1,NML=stir_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &STIR_PARAMS in parameter file'
  nml_ok=.false.
  call clean_stop
102 rewind(1)
  nml_ok=.true.

  ! ------------ DWM
  ! Compute time step to modify the stirring acc field.
  !---------------
  do i=1,noutput
     stir_timescale(i)=dble(i)*stir_delta_t
  end do
end subroutine read_stir_params

