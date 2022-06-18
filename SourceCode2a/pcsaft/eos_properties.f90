module properties

  use PARAMETERS, only: dp, NAV, RGAS
  use BASIC_VARIABLES, only: ncomp
  implicit none

  real(dp), public                                  :: pcalc_z_transfer
  real(dp), public                                  :: ztot_rho

  private
  public :: chemical_potential_tp, chemical_potential_trho,  &							  ! These are the names of the subroutines in this module. 
       state_trho, state_tp, state_ps, state_ph, Helmholtz_density,  &
       chemical_potential_tp_derivative_tp, calculate_pressure, pressure,  &
       calculate_density

contains

subroutine Helmholtz_density ( tF, rhoi_in, f_dens )

  use module_eos_derivatives, only: F_density, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), intent(out)                             :: f_dens
  !-----------------------------------------------------------------------------

  real(dp), dimension(ncomp)                        :: xF
  real(dp)                                          :: rhoF
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  rhoF = sum( rhoi_in( 1:ncomp ) )
  if ( rhoF < 1.E-50_dp ) then
     write (*,*) 'F_density: density too low'
     stop
  end if
  xF( 1:ncomp ) = rhoi_in( 1:ncomp ) / rhoF

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! Helmholtz energy density
  !-----------------------------------------------------------------------------

  call F_density ( rhoi_in, f_dens )

end subroutine Helmholtz_density



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pressure
!
! calculate the pressure in unit (Pa) and derivatives to packing
! fraction eta (=z3). The first derivative to pressure is always delivered. The
! second and third derivatives are optional arguments.
!
! this subroutine should preferentially be called from subroutines only within
! this module. From outside, rather use SR "calculate_pressure"
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure ( eta_specified, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  use module_eos_derivatives, only: f_rho, f_rho_rho, f_rho_4, z3t, rho, kT
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_specified
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out), optional            :: pcalc_z
  real(dp), intent(out), optional            :: pcalc_z2
  real(dp), intent(out), optional            :: pcalc_z3
  !-----------------------------------------------------------------------------
  real(dp)                                   :: z_res, ztot, ztot_z
  !-----------------------------------------------------------------------------

  if ( .NOT.present( pcalc_z ) ) then

     call f_rho ( eta_specified, z_res )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT

  else if ( .NOT.present( pcalc_z2 ) ) then

     call f_rho_rho ( eta_specified, z_res, ztot_z )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT
     pcalc_z = ( ztot_z*rho + ztot/z3t ) * kT

     ztot_rho = ztot_z * z3t

  else

     call f_rho_4 ( eta_specified, z_res, ztot_z, pcalc_z2, pcalc_z3 )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT
     pcalc_z = ( ztot_z*rho + ztot/z3t ) * kT

  end if

end subroutine pressure



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE calculate_pressure
!
! calculate the pressure in unit (Pa) and derivatives to packing
! fraction eta (=z3). The first derivative to pressure is always delivered. The
! second and third derivatives are optional arguments.
!
! this subroutine is called from code outside of this module
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine calculate_pressure ( tF, xF, eta_specified, rho_out, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  use module_eos_derivatives, only: rho, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  real(dp), dimension(ncomp), intent(in)     :: xF
  real(dp), intent(in)                       :: eta_specified
  real(dp), intent(out)                      :: rho_out
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out), optional            :: pcalc_z
  real(dp), intent(out), optional            :: pcalc_z2
  real(dp), intent(out), optional            :: pcalc_z3
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  call pressure ( eta_specified, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  rho_out = rho

end subroutine calculate_pressure



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE calculate_density
!
! calculate the density (of a mixture) for defined (T,p,x). A starting value
! for the packing fraction (eta_start) needs to be given. A liquid-like value
! of the starting value leads to a liquid density (eta_out and rho_out), if a
! liquid root exists at the given conditions. For ! vapor-like starting values
! one gets a vapor density (if vapor root exists).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine calculate_density ( tF, pF, xF, eta_start, eta_out, rho_out )

  use module_eos_derivatives, only: rho, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  real(dp), intent(in)                       :: pF
  real(dp), dimension(ncomp), intent(in)     :: xF
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta_out
  real(dp), intent(out)                      :: rho_out

  real(dp)                                   :: eta
  real(dp)                                   :: pcalc, pcalc_z
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rho_out = rho
  eta_out = eta

end subroutine calculate_density



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_tp
!
!> \brief calculation of the chemical potentials to specified (t,p,x)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, p, x(ncomp), eta_start).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_tp ( tF, pF, xF, eta_start, rhoi_out, mu_res, mu, lnphi_nj )

  use module_eos_derivatives, only: kT, rhoi, rho, z3t, rho_independent_quantities,  &
       F_density_rhok, F_density_rhok_rho, F_density_rhok_rhol
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out), optional :: mu
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: lnphi_nj
  !-----------------------------------------------------------------------------

  integer                                           :: i
  real(dp)                                          :: eta
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: ztot
  real(dp), dimension(ncomp)                        :: w_rk
  real(dp), dimension(ncomp)                        :: w_rk_r
  real(dp), dimension(ncomp,ncomp)                  :: w_rkrl
  real(dp), dimension(ncomp,ncomp)                  :: wig_rkrl
  real(dp), dimension(ncomp)                        :: p_rhok
  real(dp), dimension(ncomp)                        :: partial_molar_v
  real(dp)                                          :: p_rho
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, pF=pF, xF=xF )
  
  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )
  pcalc_z_transfer = pcalc_z

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------

  if ( present( lnphi_nj ) ) then

     call F_density_rhok_rho ( w_rk, w_rk_r )

     call F_density_rhok_rhol ( w_rkrl, wig_rkrl )

  else

     call F_density_rhok ( w_rk )

  end if

  if ( ztot > 0.0_dp ) mu_res( 1:ncomp ) = w_rk( 1:ncomp ) - LOG( ztot )      ! ln( phi ) for given ( T,p,x )

  !-----------------------------------------------------------------------------
  ! optional: calculate     mu_i(T,p,x)/kT
  !-----------------------------------------------------------------------------

  if ( present( mu ) ) then

     do i = 1, ncomp
        if ( rhoi(i) > 1.E-200_dp ) then
           mu(i) = w_rk(i) + log( rhoi(i) )
        else
           mu(i) = - 1.E200_dp
        end if
     end do

  end if


  !-----------------------------------------------------------------------------
  ! optional: calculate derivative matrix    n * d( ln(phi_i) ) / d( rho_j )
  !-----------------------------------------------------------------------------

  if ( present( lnphi_nj ) ) then

     p_rhok( 1:ncomp ) = rho * kT * w_rk_r( 1:ncomp ) + kT         ! in [Pa*Angstrom**3]
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     partial_molar_v( 1:ncomp ) = p_rhok( 1:ncomp ) / p_rho / rho  ! in [Angstrom**3]

     do i = 1, ncomp
        lnphi_nj( i,1:ncomp ) = rho * w_rkrl( i,1:ncomp )  &
             - partial_molar_v( i ) *partial_molar_v( 1:ncomp ) *rho**2 *p_rho / kT + 1.0_dp
     end do

  end if

end subroutine chemical_potential_tp


  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_tp_derivative_t
!
!> \brief calculation of the chemical potentials to specified (t,p,x)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, p, x(ncomp), eta_start).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_tp_derivative_tp ( tF, pF, xF, eta_start, rhoi_out, mu_res, lnphi_t, lnphi_p )

  use module_eos_derivatives, only: kT, rhoi, rho, z3t, rho_independent_quantities,  &
       F_density_rhok, F_density_rhok_rho, F_density_rhok_T, f_temp_rho
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out), optional :: lnphi_t
  real(dp), dimension(ncomp), intent(out), optional :: lnphi_p
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: ztot
  real(dp), dimension(ncomp)                        :: w_rk
  real(dp), dimension(ncomp)                        :: w_rk_r
  real(dp), dimension(ncomp)                        :: w_rk_t
  real(dp)                                          :: f_res, f_t, f_tr
  real(dp), dimension(ncomp)                        :: p_rhok
  real(dp), dimension(ncomp)                        :: partial_molar_v
  real(dp)                                          :: p_rho
  real(dp)                                          :: p_t
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, pF=pF, xF=xF )
  
  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------

  call F_density_rhok_rho ( w_rk, w_rk_r )

  mu_res( 1:ncomp ) = w_rk( 1:ncomp ) - LOG( ztot )      ! ln( phi ) for given ( T,p,x )

  p_rhok( 1:ncomp ) = rho * kT * w_rk_r( 1:ncomp ) + kT         ! in [Pa*Angstrom**3]
  p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
  partial_molar_v( 1:ncomp ) = p_rhok( 1:ncomp ) / p_rho / rho  ! in [Angstrom**3]

  if ( present( lnphi_p ) ) then

     lnphi_p( 1:ncomp ) = partial_molar_v( 1:ncomp) / kT -1.0_dp / pcalc ! in [1/Pa]

  else if ( present( lnphi_t ) ) then

     call F_density_rhok_T ( w_rk_t )

     call f_temp_rho ( f_res, f_t, f_tr )

     p_t = rho * rho * kT * f_tr + pcalc / tF

     lnphi_t( 1:ncomp ) = w_rk_t( 1:ncomp )  + 1.0_dp/tF - partial_molar_v( 1:ncomp ) * p_t / kT

  else

     write (*,*) ' use SR chemical_potential_tp_derivative_tp with at least one optional argument'
     stop

  end if

end subroutine chemical_potential_tp_derivative_tp


  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_trho
!
!> \brief calculation of the chemical potentials to specified (t,rhoi)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, rhoi(ncomp) ).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_trho ( tF, rhoi_in, mu_res, mu, w_rkrl, wig_rkrl )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: rhoi, dhs, rho_independent_quantities,  &
       density_terms, Iterate_Association, F_density_rhok, F_density_rhok_rhol
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out)           :: mu
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: w_rkrl
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: wig_rkrl
  !-----------------------------------------------------------------------------

  integer                                           :: i
  real(dp)                                          :: eta
  real(dp), dimension(ncomp)                        :: xF
  real(dp), dimension(ncomp)                        :: w_rk
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, rhoiF=rhoi_in )

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  rhoi( 1:ncomp ) = rhoi_in( 1:ncomp )
  xF( 1:ncomp ) = rhoi( 1:ncomp ) / sum( rhoi( 1:ncomp ) )

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! precalculate some density dependent terms
  !-----------------------------------------------------------------------------

  eta = PI_6 * sum( rhoi( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )

  call density_terms ( eta )

  !-----------------------------------------------------------------------------
  ! iterate fraction of unbonded association sites
  !-----------------------------------------------------------------------------
  call Iterate_Association

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------
  call F_density_rhok ( w_rk )

  mu_res( 1:ncomp ) = w_rk( 1:ncomp )             ! = ln( phi ) for given ( T,rhoi )

  do i = 1, ncomp
     if ( rhoi(i) > 1.E-200_dp ) then
        mu(i) = mu_res(i) + log( rhoi(i) )
     else
        mu(i) = - 1.E200_dp
     end if
  end do

  !-----------------------------------------------------------------------------
  ! if second derivative is required (optional argument)
  !-----------------------------------------------------------------------------

  if ( present( w_rkrl ) .AND. present( wig_rkrl ) ) then

     call F_density_rhok_rhol ( w_rkrl, wig_rkrl )

  end if

end subroutine chemical_potential_trho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_tp
!
!> \brief molar enthalpy, entropy, rho for specified (T,p,x)
!!
!! This subroutine gives p, s, h, rho, cp (optional). The required input (T,p,x)
!! and starting value for density.
!! The residual entropy s_res for variables (T,rho,x)! or (T,rho_k)! is given as
!! an optional output. This quantity is used in entropy scaling.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_tp ( tF, pF, xF, eta_start, rhoi_out, molar_rho, s, h, s_res, cp )

  use PARAMETERS, only: NAV
  use module_eos_derivatives, only: rhoi, kT, rho, z3t, rho_independent_quantities,  &
       density_terms, f_temp, f_temp_rho, f_temp_temp
  use ideal_gas_enthalpy, only: enthalpy_ig 
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: s_res
  real(dp), intent(out), optional                   :: cp
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: ztot
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: f_res, f_t, f_t2, f_tr
  real(dp)                                          :: sres
  real(dp)                                          :: p_rho, p_t
  real(dp)                                          :: cv_residual
  real(dp)                                          :: hig, sig, cpig
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate molar enthalpy
  !-----------------------------------------------------------------------------

  if ( .NOT.present( cp ) ) then

     call f_temp ( f_res, f_t )

  else

     !--------------------------------------------------------------------------
     ! partial derivatives needed for cp
     !--------------------------------------------------------------------------

     call f_temp_rho ( f_res, f_t, f_tr )

     ! ------ calculate pressure derivatives -----------------------------------
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     p_t = rho * rho * kT * f_tr + pcalc / tF                      ! in [Pa/K]

     ! ------ derivative to temperature ----------------------------------------
     call f_temp_temp ( f_t, f_t2 )

     cv_residual = - ( tF * f_t2 + 2.0_dp * f_t ) * RGAS * tF      ! in [J/(mol*K)]
     cp = cv_residual + tF * ( p_t / rho )**2 / p_rho * NAV/1.E30_dp - RGAS

  end if

  !-----------------------------------------------------------------------------
  ! molar residual enthalpy and molar residual entropy
  !-----------------------------------------------------------------------------

  h = - tF * f_t + ( ztot - 1.0_dp )              ! h_res/kT, dimensionless
  h = h * RGAS * tF                               ! h_res in [J/mol]

  sres = - tF * f_t - f_res                       ! s_res/k for given (T,rho,x), dimensionless
  s = sres + log( ztot )                          ! s_res/k = s_res_molar/R for given (T,p,x), dimensionless
  s = s * RGAS                                    ! s_res in [J/(mol*K)]

  if ( present( s_res ) ) s_res = sres

  !-----------------------------------------------------------------------------
  ! add ideal gas contributions to molar enthalpy and molar entropy
  !-----------------------------------------------------------------------------

  call enthalpy_ig ( tF, pcalc, xF, cpig, hig, sig )

  h = h + hig
  s = s + sig
  if ( present( cp ) ) cp = cp + cpig

  !-----------------------------------------------------------------------------
  ! molar density
  !-----------------------------------------------------------------------------

  molar_rho = rho / NAV * 1.E30_dp                                   ! in units [mol/m**3]
   
end subroutine state_tp



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_trho
!
!> \brief molar enthalpy, entropy, rho for specified (T, rhoi(nc))
!!
!! This subroutine gives p, s, h, rho, cp (optional). The required input (T, rhoi(nc)).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_trho ( tF, rhoi_in, pcalc, molar_rho, s, h, s_res, cp )

  use PARAMETERS, only: PI_6, NAV
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs, kT, rho, z3t, rho_independent_quantities,  &
       density_terms, f_temp, f_temp_rho, f_temp_temp
  use ideal_gas_enthalpy, only: enthalpy_ig 
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), intent(out)                             :: pcalc
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: s_res
  real(dp), intent(out), optional                   :: cp
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: ztot
  real(dp)                                          :: pcalc_z
  real(dp)                                          :: f_res, f_t, f_t2, f_tr
  real(dp)                                          :: sres
  real(dp)                                          :: p_rho, p_t
  real(dp)                                          :: cv_residual
  real(dp), dimension(ncomp)                        :: xF
  real(dp)                                          :: hig, sig, cpig
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  xF( 1:ncomp ) = rhoi_in( 1:ncomp ) / sum( rhoi_in(1:ncomp) )
  eta = PI_6 * sum( rhoi_in( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! calculate pressure for given (t,rhoi(:)) or actually for given (t,eta,x)
  !-----------------------------------------------------------------------------

  call pressure ( eta, pcalc, pcalc_z )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate molar enthalpy
  !-----------------------------------------------------------------------------

  if ( .NOT.present( cp ) ) then

     call f_temp ( f_res, f_t )

  else

     !--------------------------------------------------------------------------
     ! partial derivatives needed for cp
     !--------------------------------------------------------------------------

     call f_temp_rho ( f_res, f_t, f_tr )

     ! ------ calculate pressure derivatives -----------------------------------
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     p_t = rho * rho * kT * f_tr + pcalc / tF                      ! in [Pa/K]

     ! ------ derivative to temperature ----------------------------------------
     call f_temp_temp ( f_t, f_t2 )

     cv_residual = - ( tF * f_t2 + 2.0_dp * f_t ) * RGAS * tF      ! in [J/(mol*K)]
     cp = cv_residual + tF * ( p_t / rho )**2 / p_rho * NAV/1.E30_dp - RGAS

  end if

  !-----------------------------------------------------------------------------
  ! molar residual enthalpy and molar residual entropy
  !-----------------------------------------------------------------------------

  h = - tF * f_t + ( ztot - 1.0_dp )              ! h_res/kT, dimensionless
  h = h * RGAS * tF                               ! h_res in [J/mol]

  sres = - tF * f_t - f_res                      ! s_res/k for given (T,rho,x), dimensionless
  s = sres + log( ztot )                         ! s_res/k = s_res_molar/R for given (T,p,x), dimensionless
  s = s * RGAS                                    ! s_res in [J/(mol*K)]

  if ( present( s_res ) ) s_res = sres

  !-----------------------------------------------------------------------------
  ! add ideal gas contributions to molar enthalpy and molar entropy
  !-----------------------------------------------------------------------------

  call enthalpy_ig ( tF, pcalc, xF, cpig, hig, sig )

  h = h + hig
  s = s + sig
  if ( present( cp ) ) cp = cp + cpig

  !-----------------------------------------------------------------------------
  ! molar density
  !-----------------------------------------------------------------------------

  molar_rho = rho / NAV * 1.E30_dp                                   ! in units [mol/m**3]
   
end subroutine state_trho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_ps
!
! input, specifications:     pF, sF, xF(:)
! input, starting values:    eta_start, tcalc
! output:                    tcalc, rho, h,  cp (optional)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_ps ( pF, sF, xF, eta_start, tcalc, rhoi_out, molar_rho, h, cp )

  use BASIC_VARIABLES, only: ncomp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: pF
  real(dp), intent(in)                              :: sF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), intent(inout)                           :: tcalc
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: cp

  !-----------------------------------------------------------------------------
  integer, parameter                                :: maxiter = 100
  integer                                           :: i_iter
  real(dp)                                          :: t_org
  real(dp)                                          :: f_error, dfdt
  real(dp), parameter                               :: f_tolerance = 1.E-10_dp
  real(dp), parameter                               :: delta_t = 1.E-4_dp
  real(dp)                                          :: sF1, sF2
  real(dp)                                          :: t_step
  logical                                           :: converged
  !-----------------------------------------------------------------------------

  i_iter = 0

  converged = .false.

  do while ( .NOT.converged .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     t_org = tcalc
     tcalc = t_org - delta_t
     call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF1, h )

     tcalc = t_org
     if ( present( cp ) ) then
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF2, h, cp=cp )
     else
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF2, h )
     end if

     f_error = ( sF2 - sF )
     dfdt = ( sF2 - sF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  50.0_dp, t_step )
     t_step = max( -50.0_dp, t_step )

     tcalc = tcalc - t_step

     IF ( ABS( f_error ) < f_tolerance ) converged = .true.

  end do

  if ( .NOT.converged ) write (*,*) 'ps_state not converged'

end subroutine state_ps



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_ph
!
! input, specifications:     pF, hF, xF(:)
! input, starting values:    eta_start, tcalc
! output:                    tcalc, rho, s,  cp (optional)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_ph ( pF, hF, xF, eta_start, tcalc, rhoi_out, molar_rho, s, cp )

  use BASIC_VARIABLES, only: ncomp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: pF
  real(dp), intent(in)                              :: hF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), intent(inout)                           :: tcalc
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out), optional                   :: cp

  !-----------------------------------------------------------------------------
  integer, parameter                                :: maxiter = 100
  integer                                           :: i_iter
  real(dp)                                          :: t_org
  real(dp)                                          :: f_error, dfdt
  real(dp), parameter                               :: f_tolerance = 1.E-10_dp
  real(dp), parameter                               :: delta_t = 1.E-4_dp
  real(dp)                                          :: hF1, hF2
  real(dp)                                          :: t_step
  logical                                           :: converged
  !-----------------------------------------------------------------------------

  i_iter = 0

  converged = .false.

  do while ( .NOT.converged .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     t_org = tcalc
     tcalc = t_org - delta_t
     call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF1 )

     tcalc = t_org
     if ( present( cp ) ) then
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF2, cp=cp )
     else
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF2 )
     end if

     f_error = ( hF2 - hF )
     dfdt = ( hF2 - hF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  50.0_dp, t_step )
     t_step = max( -50.0_dp, t_step )

     tcalc = tcalc - t_step

     IF ( ABS( f_error ) < f_tolerance ) converged = .true.

  end do

  if ( .NOT.converged ) write (*,*) 'ph_state not converged'

end subroutine state_ph



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE check_input
!
! verify that all input variables are physically reasonable, otherwise stop
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_input ( tF, pF, xF, rhoiF )

  use PARAMETERS, only: machine_eps
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in), optional                    :: pF
  real(dp), dimension(ncomp), intent(in), optional  :: xF
  real(dp), dimension(ncomp), intent(in), optional  :: rhoiF
  !-----------------------------------------------------------------------------
  integer                                           :: i
  real(dp)                                          :: sum_x
  logical                                           :: immediate_exit
  !-----------------------------------------------------------------------------

  immediate_exit = .false.

  if ( present( pF ) ) then

     !--------------------------------------------------------------------------
     ! consider case where (T,p,x) are specified variables
     !--------------------------------------------------------------------------

     !-----------------------------------------------------------------------------
     ! check for NaN and expected array-size
     !-----------------------------------------------------------------------------
     if ( tF /= tF .or. pF /= pF ) then
        write (*,*) ' input is "not a number" !'
        immediate_exit = .true.
     end if

     if ( tF < 0.0_dp .or. pF < 0.0_dp ) then
        write (*,*) ' input is not physical !'
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of pressure (pF)
     !--------------------------------------------------------------------------
     if ( pF < 1.E-50_dp ) then
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of mole fractions xF
     !--------------------------------------------------------------------------
     if ( size( xF, 1) /= ncomp ) then
        write (*,*) ' array-size of xF is not appropriate !'
        immediate_exit = .true.
     end if

     do i = 1, ncomp
        if ( xF(i) /= xF(i) ) then
           write (*,*) ' composition input x is "not a number" ! Species:',i
           immediate_exit = .true.
        end if
     end do

     if ( MINVAL( xF( 1:ncomp ) ) < 0.0_dp ) then
        write (*,*) ' mole fraction x is negative, of species', MINLOC( xF(1:ncomp) )
        immediate_exit = .true.
     end if

     sum_x = sum( xF( 1:ncomp ) )
     if ( abs( sum_x - 1.0_dp ) > 1.E4_dp * machine_eps ) then
        if ( abs( sum_x ) < machine_eps ) then
           write (*,*) ' mole fracitons are zero'
        else
           write (*,*) ' sum of mole fracitons is not equal one'
           ! xF( 1:ncomp ) = xF( 1:ncomp ) / sum_x    ! normalize mole fractions
        end if
        immediate_exit = .true.
     end if

  else

     !--------------------------------------------------------------------------
     ! consider case where (T,rhoi) are specified variables
     !--------------------------------------------------------------------------

     !--------------------------------------------------------------------------
     ! check for NaN and expected array-size
     !--------------------------------------------------------------------------
     if ( tF /= tF ) then
        write (*,*) ' input is "not a number" !'
        immediate_exit = .true.
     end if

     if ( tF < 0.0_dp ) then
        write (*,*) ' input is not physical !'
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of species densities rhoi
     !--------------------------------------------------------------------------
     if ( size( rhoiF, 1) /= ncomp ) then
        write (*,*) ' array-size of rhoiF is not appropriate !'
        immediate_exit = .true.
     end if

     do i = 1, ncomp
        if ( rhoiF(i) /= rhoiF(i) ) then
           write (*,*) ' species desnity rhoi is "not a number" ! Species:',i
           immediate_exit = .true.
        end if
     end do

     if ( MINVAL( rhoiF( 1:ncomp ) ) < 0.0_dp ) then
        write (*,*) ' species desnity rhoi is negative, of species', MINLOC( rhoiF(1:ncomp) )
        immediate_exit = .true.
     end if

     if ( sum( rhoiF(:) ) < 1.E-50_dp ) then
        write (*,*) ' overall density too low', sum( rhoiF(:) )
        immediate_exit = .true.
     end if

  end if

  !-----------------------------------------------------------------------------
  ! if unphysical input specifications: stop
  !-----------------------------------------------------------------------------
  if ( immediate_exit ) then
     write (*,*) 'check_input: unphysical input specs', tF, pF, xF(:)
     stop
  end if

end subroutine check_input

end module properties

