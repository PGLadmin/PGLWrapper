!> \file load_pure_and_kij_parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE load_pure_and_kij_parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module ideal_gas_enthalpy

  use PARAMETERS, only: dp, machine_eps, RGAS
  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: cp_coefficients
  implicit none

  public :: enthalpy_ig

contains



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    
subroutine enthalpy_ig ( t, p, x, cpig, hig, sig )

  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x
  real(dp), intent(out)                      :: cpig
  real(dp), intent(out)                      :: hig
  real(dp), intent(out)                      :: sig
  !-----------------------------------------------------------------------------

  call enthalpy_ig_joback ( t, p, x, cpig, hig, sig )
    
end subroutine enthalpy_ig
    
    
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    
subroutine enthalpy_ig_joback ( t, p, x, cpig, hig, sig )

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x
  real(dp), intent(out)                      :: cpig
  real(dp), intent(out)                      :: hig
  real(dp), intent(out)                      :: sig

  !-----------------------------------------------------------------------------
  integer                                    :: i
  real(dp), dimension(ncomp,4)               :: coef_ig
  real(dp)                                   :: t2, t3, t4
  real(dp)                                   :: p0, t0, t02, t03, t04
  real(dp), dimension(ncomp)                 :: cp_ig, h_ig, s_ig
  !-----------------------------------------------------------------------------

  p0 = 1.E5_dp

  t0 = 298.15_dp
  t02 = t0 * t0
  t03 = t02 * t0
  t04 = t02 * t02

  t2 = t * t
  t3 = t2 * t
  t4 = t2 * t2

  do i = 1, ncomp

     coef_ig(i,1) = ( cp_coefficients(i,1) - 37.93_dp )
     coef_ig(i,2) = ( cp_coefficients(i,2) + 0.21_dp )
     coef_ig(i,3) = ( cp_coefficients(i,3) - 3.91E-4_dp )
     coef_ig(i,4) = ( cp_coefficients(i,4) + 2.06E-7_dp )

     cp_ig(i)= coef_ig(i,1) + coef_ig(i,2) * t + coef_ig(i,3) * t2 + coef_ig(i,4) * t3
     h_ig(i) = coef_ig(i,1) *( t - t0 ) + coef_ig(i,2)/2.0_dp *( t2 - t02 )  &
          + coef_ig(i,3)/3.0_dp *( t3 - t03 ) + coef_ig(i,4)/4.0_dp *( t4 - t04 )
     s_ig(i) = coef_ig(i,1) *log( t/t0 ) + coef_ig(i,2) *( t - t0 )  &
          + coef_ig(i,3)/2.0_dp *( t2 - t02 ) + coef_ig(i,4)/3.0_dp *( t3 - t03 ) - RGAS *log( p/p0 )

  end do

  hig = sum( x(1:ncomp ) * h_ig(1:ncomp) )
  cpig= sum( x(1:ncomp ) * cp_ig(1:ncomp) )
  sig = sum( x(1:ncomp) * ( s_ig(1:ncomp) - RGAS * log( x(1:ncomp) ) ), MASK = x > 0.0_dp )

end subroutine enthalpy_ig_joback


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW    

subroutine enthalpy_ig_QSPR ( t, p, x, cpig, hig, sig )

  use pcsaft_pure_and_binary_parameters, only: mseg, sigma, epsilon_k, eps_hb,  &
       kap_hb, quadru_moment

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x
  real(dp), intent(out)                      :: cpig
  real(dp), intent(out)                      :: hig
  real(dp), intent(out)                      :: sig

  !-----------------------------------------------------------------------------
  integer                                    :: i
  real(dp)                                   :: p0
  real(dp)                                   :: t0
  real(dp), dimension(ncomp)                 :: cp_ig, h_ig, s_ig
  real(dp)                                   :: p1, p2, p3, p4, p5, p6
  real(dp)                                   :: c1a, c2a, c3a, c4a, c5a, c6a
  real(dp)                                   :: c1b, c2b, c3b, c4b, c5b, c6b
  real(dp)                                   :: a_coeff, b_coeff
  real(dp)                                   :: t300, t400
  real(dp)                                   :: ICPC300, ICPC400
  real(dp)                                   :: sig3
  !-----------------------------------------------------------------------------

  t0 = 298.0_dp
  p0 = 1.E5_dp

  t300 = 300.0_dp
  t400 = 400.0_dp

  !-----------------------------------------------------------------------------
  ! nAnP
  !-----------------------------------------------------------------------------

  ! ICPC @ 300 K

  c1a = -5763.04893_dp
  c2a = 1232.30607_dp
  c3a = -239.3513996_dp
  c4a = 0.0_dp
  c5a = 0.0_dp
  c6a = -15174.28321_dp

  ! ICPC @ 400 K

  c1b = -8171.26676935062_dp
  c2b = 1498.01217504596_dp
  c3b = -315.515836223387_dp
  c4b = 0.0_dp
  c5b = 0.0_dp
  c6b = -19389.5468655708_dp

  !-----------------------------------------------------------------------------
  ! nAP
  !-----------------------------------------------------------------------------

  !
  !! ICPC @ 300K:
  ! c1a = 5177.19095226181
  ! c2a = 919.565206504576
  ! c3a = -108.829105648889
  ! c4a = 0
  ! c5a = -3.93917830677682
  ! c6a = -13504.5671858292
  !
  !! ICPC @ 400K:
  ! c1b = 10656.1018362315
  ! c2b = 1146.10782703748
  ! c3b = -131.023645998081
  ! c4b = 0
  ! c5b = -9.93789225413177
  ! c6b = -24430.12952497

  !-----------------------------------------------------------------------------
  ! AP
  !-----------------------------------------------------------------------------

  !! ICPC @ 300K:
  ! c1a = 3600.32322462175
  ! c2a = 1006.20461224949
  ! c3a = -151.688378113974
  ! c4a = 7.81876773647109d-07
  ! c5a = 8.01001754473385
  ! c6a = -8959.37140957179

  !! ICPC @ 400K:
  ! c1b = 7248.0697641199
  ! c2b = 1267.44346171358
  ! c3b = -208.738557800023
  ! c4b = 0.000170238690157906
  ! c5b = -6.7841792685616
  ! c6b = -12669.4196622924


  do i = 1, ncomp

     sig3 = sigma(i)**3
     p1 = mseg(i) * epsilon_k(i) / t
     p2 = mseg(i) * sig3
     p3 = mseg(i) * sig3 * epsilon_k(i) / t
     p4 = mseg(i) * sig3*sig3 * kap_hb(i,i) *( exp( eps_hb(i,i,1,2)/t ) - 1.0_dp )
     p5 = mseg(i) * sig3 * quadru_moment(i)
     p6 = 1.0_dp

     ICPC300 = c1a*p1/t300 + c2a*p2 + c3a*p3/t300 + c4a*p4/t300 + c5a*p5 + c6a*p6 
     ICPC300 = ICPC300 / 1000.0_dp

     ICPC400 = c1b*p1/t400 + c2b*p2 + c3b*p3/t400 + c4b*p4/t400 + c5b*p5 + c6b*p6 
     ICPC400 = ICPC400 / 1000.0_dp

     !--------------------------------------------------------------------------
     ! linear approxmation
     !--------------------------------------------------------------------------

     b_coeff = ( icpc400 - icpc300 ) / ( t400 - t300 )
     a_coeff = icpc300 - b_coeff * t300

     cp_ig(i)= a_coeff + b_coeff * t
     h_ig(i) = a_coeff * ( t - t0 ) + 0.5_dp * b_coeff *( t*t - t0*t0 )
     s_ig(i) = a_coeff * log( t/t0 ) + b_coeff *( t - t0 ) - RGAS * log( p/p0 )

  end do

  !-----------------------------------------------------------------------------
  ! Ideal gas properties
  !-----------------------------------------------------------------------------

  hig = sum( x(1:ncomp ) * h_ig(1:ncomp) )
  cpig= sum( x(1:ncomp ) * cp_ig(1:ncomp) )
  sig = sum( x(1:ncomp) * ( s_ig(1:ncomp) - RGAS * log( x(1:ncomp) ) ), MASK = x > 0.0_dp )

end subroutine enthalpy_ig_QSPR



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine enthalpy_ig_corr ( t, p, x, coef_ig, cpig, hig, sig )

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x
  real(dp), dimension(ncomp,5), intent(in)   :: coef_ig
  real(dp), intent(out)                      :: cpig
  real(dp), intent(out)                      :: hig
  real(dp), intent(out)                      :: sig

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  integer                                    :: i
  real(dp)                                   :: p0
  real(dp)                                   :: t0, t02, t03, t04
  real(dp)                                   :: t2, t3, t4
  real(dp), dimension(ncomp)                 :: cp_ig, h_ig, s_ig
  !-----------------------------------------------------------------------------

  p0 = 1.E5_dp

  t0 = 298.15_dp
  t02 = t0 * t0
  t03 = t02 * t0
  t04 = t02 * t02

  t2 = t * t
  t3 = t2 * t
  t4 = t2 * t2

  do i = 1, ncomp

     cp_ig(i)= coef_ig(i,1) + coef_ig(i,2) * t + coef_ig(i,3) * t2 + coef_ig(i,4) * t3
     h_ig(i) = coef_ig(i,1) *( t - t0 ) + coef_ig(i,2)/2.0_dp *( t2 - t02 )  &
             + coef_ig(i,3)/3.0_dp *( t3 - t03 ) + coef_ig(i,4)/4.0_dp *( t4 - t04 )
     s_ig(i) = coef_ig(i,1) * log( t/t0 ) + coef_ig(i,2) *( t - t0 )  &
             + coef_ig(i,3)/2.0_dp *( t2 - t02 ) + coef_ig(i,4)/3.0_dp *( t3 - t03 ) - RGAS *log( p/p0 )

  end do

  hig = sum( x(1:ncomp ) * h_ig(1:ncomp) )
  cpig= sum( x(1:ncomp ) * cp_ig(1:ncomp) )
  sig = sum( x(1:ncomp) * ( s_ig(1:ncomp) - RGAS * log( x(1:ncomp) ) ), MASK = x > 0.0_dp )

end subroutine enthalpy_ig_corr

end module ideal_gas_enthalpy
