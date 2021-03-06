module stir_parameters
  use amr_commons
  implicit none
  logical :: stir_initialized=.false.
  integer :: stir_seed=232342
  integer,parameter :: nstir=8, stir_nvar=3
  real(dp) :: stir_norm=1d0,stir_index=-2d0, stir_kmin=0.1d0, stir_kmax=3d0
  real(dp),dimension(nstir,nstir,nstir)::stir_amp_x,stir_amp_y,stir_amp_z,stir_phi_x, stir_phi_y, stir_phi_z
  real(dp),dimension(nstir)::stir_kx,stir_ky,stir_kz
end module stir_parameters
