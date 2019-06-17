module jet_parameters
  use amr_parameters, only : dp
  implicit none
  real(dp) :: v_jet = 1d7,  theta0_jet=0.3d0, f_jet = 0.3d0
  integer :: rin_jet = 4, rout_jet = 8
  real(dp) :: cos_theta0_jet = 0.955 ! cos theta of 0.3 
  integer :: jet_cell_supersample = 3 
  logical :: jet_refine=.false.
  logical :: test_jet_problem=.false.
  real(dp) :: test_rho=3d-21
  real(dp) :: test_radius=1.5d18
  real(dp) :: test_spin=1.4d-15
  real(dp) :: time_jet=1d3 ! time to turn on jet (after 1e3 years)
  real(dp), dimension(100) :: lgm_zams, lgl_zams, lgr_zams, beta_zams, betaR_zams
  integer :: num_zams
  logical :: initialized_zams = .false.

end module jet_parameters
