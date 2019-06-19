subroutine read_stir_params(nml_ok)
  use stir_parameters
  implicit none
  logical::nml_ok
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/stir_params/stir_seed,stir_norm,stir_index,& 
        stir_kmin,stir_kmax
  rewind(1)
  read(1,NML=stir_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &STIR_PARAMS in parameter file'
  nml_ok=.false.
  call clean_stop
102 rewind(1)
  nml_ok=.true.
end subroutine read_stir_params

