!> \file pure_par_fit.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module diagrams
!> \brief module for generating phase diagrams or property diagrams
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module module_pure_diagrams

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                         :: outp = 0

  private
  public :: pure_vapor_pressure, pure_boiling_temperature, pure_diagrams

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_diagrams
!> \brief determine pure component diagrams (phase diagrams, isotherms, isobars)
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_diagrams

  !-----------------------------------------------------------------------------
  integer                                    :: calcopt
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! make choice for calculation option
  !-----------------------------------------------------------------------------

  if ( ncomp /= 1 ) then
     write (*,*) 'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE:'
     write (*,*) '    ./input_file/INPUT.INP'
     stop
  end if

  write (*,*) ' '
  write (*,'(a)') ' calculate VLE                      (1)'
  write (*,'(a)') ' calculate isotherm                 (2)'
  write (*,'(a)') ' calculate isobar                   (3)'
  write (*,'(a)') ' calculate critical point           (4)'
  write (*,'(a)') ' calculate second virial coeff.     (5)'
  read  (*,*) calcopt

  !-----------------------------------------------------------------------------
  ! calculation options
  !-----------------------------------------------------------------------------

  if ( calcopt == 1 ) then

     call pure_VLE_diagram

  else if ( calcopt == 2 .or. calcopt == 3 ) then

     call pure_isotherm_or_isobar_diagram ( calcopt )

  else if ( calcopt == 4 ) then

     call pure_critical_point_with_properties

  else if ( calcopt == 5 ) then

     call pure_virial_coefficients

  else
 
     write (*,*) 'pure_diagrams: make proper choice for calculation option'
     stop

  end if

end subroutine pure_diagrams


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_VLE_diagram
!> \brief determine pure component VLE properties
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


subroutine pure_VLE_diagram

  use BASIC_VARIABLES, only: t_input
  use input_output, only: u_in_T, pure_output_coex, unit_t
  use module_critical_point, only: critical_point_pure
  use properties, only: state_trho

  !-----------------------------------------------------------------------------
  integer                                    :: i
  integer                                    :: steps
  real(dp)                                   :: t0
  real(dp)                                   :: tc, pc, rhoc
  real(dp)                                   :: tlv, plv
  real(dp)                                   :: end_t
  real(dp)                                   :: rhoL, rhoV
  real(dp)                                   :: hL, hV, hc
  real(dp)                                   :: sL, sV, sc
  real(dp)                                   :: cpL, cpV
  real(dp)                                   :: molar_rho
  real(dp)                                   :: pcalc
  real(dp), dimension(ncomp)                 :: rhoi
  logical                                    :: converg
  logical                                    :: converg_c
  !-----------------------------------------------------------------------------

  steps = 40

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call pure_output_coex ( tlv, plv, rhoL, rhoV, sL, sV, hL, hV, cpL=cpL, cpV=cpV, header=.true. )

  t0 = t_input

  !-----------------------------------------------------------------------------
  ! calculate the critical point to guide in specifying a T-range
  !-----------------------------------------------------------------------------
  tc = t0 * 3.0_dp        ! starting value for critical point
  call critical_point_pure ( tc, converg_c, pc, rhoc )

  if ( converg_c ) write (*,'(a,G18.9,a)') ' the critical temperature is: Tc= ', tc,' [K]'
  write (*,*) ' '

  !-----------------------------------------------------------------------------
  ! let user define the end-temperature for the VLE diagram
  !-----------------------------------------------------------------------------
  write (*,'(3a)') ' specify the end-temperature, in [',unit_t,']:'
  read (*,*) end_t
  end_t = end_t + u_in_T

  !-----------------------------------------------------------------------------
  ! calc. the VLE diagram in 'steps' as the number of grid points
  !-----------------------------------------------------------------------------
  do i = 0, steps

     tlv = t0 + real( i, KIND=dp ) / real( steps, KIND=dp ) * ( end_t - t0 )

     !--------------------------------------------------------------------------
     ! solve phase equilibrium condition for defined temp. tlc
     !--------------------------------------------------------------------------
     call pure_vapor_pressure ( tlv, converg, plv, rhoL, rhoV )

     if ( converg ) then

        !-----------------------------------------------------------------------
        ! calc. coexisting properties, such as enthalpies etc.
        !-----------------------------------------------------------------------
        rhoi(1) = rhoL
        call state_trho ( tlv, rhoi, pcalc, molar_rho, sL, hL, cp=cpL )
        rhoi(1) = rhoV
        call state_trho ( tlv, rhoi, pcalc, molar_rho, sV, hV, cp=cpV )

        !-----------------------------------------------------------------------
        ! output to file (and to terminal)
        !-----------------------------------------------------------------------
        call pure_output_coex ( tlv, plv, rhoL, rhoV, sL, sV, hL, hV, cpL=cpL, cpV=cpV, header=.false. )

     end if


  end do

  !-----------------------------------------------------------------------------
  ! as a last point, add the critical point to the output file
  !-----------------------------------------------------------------------------
  if ( converg_c ) then

     rhoi(1) = rhoc
     call state_trho ( tc, rhoi, pcalc, molar_rho, sc, hc )

     call pure_output_coex( tc, pc, rhoc, rhoc, sc, sc, hc, hc, header=.false. )

  end if

end subroutine pure_VLE_diagram


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_vapor_pressure
!> \brief determine pure component vapor pressure for given T
!!
!! direct substitution, taking advantage of the fact that chemical potential
!! mu_L is a weak function of p, and that
!! mu_V(p) approx.=  mu_V(p0) + ln( p /p0 ) + integral(vol_V/kT * dp)
!! where mu is chem. pot. divided by kT.
!!
!! The method is refined by estimating the pressure dependence of the vapor
!! phase and of the liquid phase. That improves convergence rate.
!!
!! The subroutine does not require starting values in pressure (p_start). If
!! for the relevant initial value for pressure a liquid and a vapor root are
!! not feasible, then the inflection point of the isotherm is found and taken as
!! a starting point for p.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_vapor_pressure ( tF, converg, plv, rhoL, rhoV, p_start )

  use module_critical_point, only: critical_point_pure
  use properties, only: chemical_potential_tp, pcalc_z_transfer
  use module_eos_derivatives, only: kT, z3

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  logical, intent(out)                       :: converg
  real(dp), intent(out)                      :: plv
  real(dp), intent(out)                      :: rhoL
  real(dp), intent(out)                      :: rhoV
  real(dp), intent(in), optional             :: p_start
  !-----------------------------------------------------------------------------

  real(dp), parameter                        :: tolerance = 1.E-10_dp
  integer, parameter                         :: max_iter = 100

  integer                                    :: i
  integer                                    :: count
  real(dp), dimension(2)                     :: eta_start
  real(dp)                                   :: eta_sta
  real(dp), dimension(ncomp)                 :: xF
  real(dp)                                   :: pF
  real(dp)                                   :: eta
  real(dp)                                   :: pcalc_z, pcalc_z2, pcalc_z3
  real(dp), dimension(ncomp)                 :: rhoi_L, rhoi_V
  real(dp), dimension(ncomp)                 :: mu_res_L, mu_res_V
  real(dp), dimension(ncomp)                 :: mu_L, mu_V
  logical                                    :: converg_inflection

  real(dp)                                   :: p0
  real(dp)                                   :: eta_L, eta_V
  real(dp)                                   :: vV_kt
  real(dp)                                   :: error
  real(dp)                                   :: p_inflex
  real(dp)                                   :: alphaL
  real(dp)                                   :: mu_V_correction
  !-----------------------------------------------------------------------------

  converg = .false.
  converg_inflection = .false.

  !-----------------------------------------------------------------------------
  ! check if input is ok
  !-----------------------------------------------------------------------------
  if ( ncomp /= 1 ) then
     write (*,*) 'use SR pure_vapor_pressure only for pure species'
     stop
  end if
  
  !-----------------------------------------------------------------------------
  ! define xF
  !-----------------------------------------------------------------------------
  xF(1) = 1.0_dp

  if ( outp >= 1 ) write (*,*) ' '
  if ( outp >= 1 ) write (*,*) 'ENTERING AT T = ',tF

  !-----------------------------------------------------------------------------
  ! starting value for pressure, and for liquid and vapor densities
  !-----------------------------------------------------------------------------
  if ( present( p_start ) ) then
     pF = p_start
  else
     pF = 1.0E3_dp               ! default starting value for pressure
  end if

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.0E-5_dp

  !-----------------------------------------------------------------------------
  ! phase equilibrium calculation, by sucessive substitution
  !-----------------------------------------------------------------------------

  count = 0
  error = tolerance + 1.0_dp
 
  do while ( error > tolerance .AND. count < max_iter )

     count = count + 1

     call chemical_potential_tp ( tF, pF, xF, eta_start(1), rhoi_L, mu_res_L, mu=mu_L )
     eta_L = z3
     rhoL = rhoi_L(1)
     !vL_kt = 1.0_dp / ( rhoi_L(1) * kT )
     alphaL = rhoL / eta_L / pcalc_z_transfer

     call chemical_potential_tp ( tF, pF, xF, eta_start(2), rhoi_V, mu_res_V, mu=mu_V )
     eta_V = z3
     rhoV = rhoi_V(1)
     vV_kt = 1.0_dp / ( rhoV * kT ) -  1.0_dp / pF
     !alphaV = rhoi_V(1) / eta_V / pcalc_z_transfer

     !--------------------------------------------------------------------------
     ! if vapor & liquid root are not both feasible (weak value of p, or t > tc)
     !--------------------------------------------------------------------------
     if ( abs( rhoV / rhoL - 1.0_dp ) < 1.E-6_dp ) then

        ! ----- determine inflection point of pressure isothem -----------------
        error = 100.0_dp
        if ( .NOT.converg_inflection ) then
           eta_sta = 0.2_dp
           call pressure_inflection_point( eta_sta, eta, p_inflex, pcalc_z, pcalc_z2,  &
                                           pcalc_z3, converg_inflection )
           if ( converg_inflection .AND. pcalc_z > 0.d0 ) then
              if ( outp >= 1 ) write (*,*) 'certainly above critical point'
              RETURN
           else
              pF = p_inflex
              if ( outp >= 1 ) write (*,*) 'new pressure after inflex', pF
              CYCLE
           end if
        else
           if ( outp >= 1 ) write (*,*) 'this case should not appear'
           pF = p_inflex
        end if

     else

        ! ----- normal procedure -----------------------------------------------
        eta_start(1) = eta_L
        eta_start(2) = eta_V

     end if
     p0 = pF

     !--------------------------------------------------------------------------
     ! next (successive substitution) value for p. Solve eq. self-consistently
     !--------------------------------------------------------------------------
     pF = p0 * exp( mu_L(1) - mu_V(1) )      ! starting value for pF
     if ( (rhoL+alphaL*(pF-p0)) > 1.0E-4_dp ) then
        do i = 1, 10                           ! self-consistency by 10 simple recursive steps
           ! pF = p0 * exp( mu_L(1) + vL_kt*(pF-p0) - mu_V(1) - vV_kt*(pF-p0) )
           mu_V_correction = vV_kt*(pF-p0)
           if ( abs( mu_V_correction ) >  1.0_dp ) mu_V_correction = sign( 1.0_dp, mu_V_correction )
           if ( count <= 2 ) mu_V_correction = 0.0_dp
           pF = p0 * exp( mu_L(1) + 1.0_dp/alphaL/kT*log( (rhoL+alphaL*(pF-p0)) / rhoL ) - mu_V(1) - mu_V_correction )
           ! pF = p0 * exp( mu_L(1) + 1.0_dp/alphaL/kT*log( (rhoL+alphaL*(pF-p0)) / rhoL )  &
           !             - mu_V(1) - 1.0_dp/alphaV/kT*log( (rhoV+alphaV*(pF-p0)) / rhoV ) + log(pF/p0) )
        end do
     end if
     ! pF = ( mu_V(1) - p0/(rhoV*kT) - mu_L(1) + p0/(rhoL*kT) ) / ( 1.0_dp/(rhoL*kT) - 1.0_dp/(rhoV*kT) )

     ! if ( pF > 4.0_dp * p0 ) pF = 5.0_dp * p0
     ! if ( pF < 0.25_dp * p0 ) pF = 0.2_dp * p0
     if ( pF < 1.E-10_dp )   pF = 1.E-10_dp

     !--------------------------------------------------------------------------
     ! error
     !--------------------------------------------------------------------------
     error = abs( pF / p0 - 1.0_dp )

     if ( outp >= 1 ) write (*,'(i2,4G22.14)') count, error, pF

  end do

  !-----------------------------------------------------------------------------
  ! for convergence set plv
  !-----------------------------------------------------------------------------
  if ( error < tolerance ) then
     converg = .true.
     plv = pF
  end if

end subroutine pure_vapor_pressure


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_boiling_temperature
!> \brief determine pure component vapor pressure for given T
!!
!! The subroutine does not require starting values in temperature (t_start). The
!! code first searches a temperature, where a liquid exists and a temperature,
!! where a vapor exists. Then, the Clausius-Clapeyron equation is used to 
!! iterate the boiling temperature, using the SR pure_vapor_pressure.
!!
!! NOT implemented:
!!
!! A good alternative is a direct substitution method, taking advantage of the
!! fact that chemical potential mu in both phases (L and V) can be Taylor
!! expanded to first order with -s as the derivative of mu to T. For a starting
!! value for T in L and a different starting value for T in the V, one finds a
!! next guess for the bubble-point T as a crossing point of the Taylor expanded
!! mu-lines.
!!
!! tF = ( mu_res_L(1) + s_res_L*tFL ) / ( mu_V(1) - log( pF / ( KBOL30 * tF ) ) + s_res_L )
!! (this eq. needs to be solved self-consistently, by recursively evaluating it)
!!
!! note: mu (or mu_V, mu_L) are actually mu/(kT)
!!       mu_res (or mu_res_V mu_res_L) are actually mu_res/(kT)
!!       mu/k  =  T * mu/(kT)  =  mu_res/(kT) + T * log( p / ( KBOL30 * T ) )
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_boiling_temperature ( pF, converg, tlv, rhoL, rhoV, t_start )

  use PARAMETERS, only: KBOL30
  use module_critical_point, only: critical_point_pure
  use properties, only: chemical_potential_tp, state_tp, state_trho
  use module_eos_derivatives, only: z3t
  use bubble_dew_point, only: bubble_point_rachford_rice

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: pF
  logical, intent(out)                       :: converg
  real(dp), intent(out)                      :: tlv
  real(dp), intent(out)                      :: rhoL
  real(dp), intent(out)                      :: rhoV
  real(dp), intent(in), optional             :: t_start
  !-----------------------------------------------------------------------------

  real(dp), parameter                        :: tolerance = 1.E-11_dp
  integer, parameter                         :: max_iter = 100

  integer                                    :: count
  real(dp), dimension(2)                     :: eta_start
  real(dp)                                   :: tF
  real(dp), dimension(ncomp)                 :: xF
  real(dp), dimension(ncomp)                 :: rhoi_L, rhoi_V
  real(dp)                                   :: tFL, tFV
  real(dp)                                   :: hL, hV, sL, sV
  real(dp)                                   :: sres_L, sres_V
  real(dp)                                   :: s_res_L, s_res_V
  real(dp)                                   :: rho_ideal_gas
  logical                                    :: found_L, found_V

  real(dp)                                   :: t0
  real(dp)                                   :: error
  real(dp)                                   :: damp
  real(dp)                                   :: delta_inversT
  real(dp)                                   :: molar_rhoL, molar_rhoV
  real(dp)                                   :: dlnp_d_inversT
  real(dp)                                   :: pcalc, plv
  real(dp)                                   :: tc, pc, rhoc
  logical                                    :: converg_c
  !-----------------------------------------------------------------------------

  converg = .false.
  converg_c = .false.

  !-----------------------------------------------------------------------------
  ! check if input is ok
  !-----------------------------------------------------------------------------
  if ( ncomp /= 1 ) then
     write (*,*) 'use SR pure_vapor_pressure only for pure species'
     stop
  end if
  
  !-----------------------------------------------------------------------------
  ! define xF
  !-----------------------------------------------------------------------------
  xF(1) = 1.0_dp

  if ( outp >= 1 ) write (*,*) ' '
  if ( outp >= 1 ) write (*,*) 'ENTERING AT p = ',pF

  !-----------------------------------------------------------------------------
  ! starting value for T, and for liquid and vapor densities
  !-----------------------------------------------------------------------------
  if ( present( t_start ) ) then
     tF = t_start
  else
     tF = 300.0_dp               ! default starting value for temperature
  end if

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.0E-5_dp

  ! p0 = pF
  ! call bubble_point_rachford_rice ( .true., tF, p0, eta_start, xF, xi2, rhoi_L, rhoi_V, .true., converg )
  ! return


  !-----------------------------------------------------------------------------
  ! improve starting value for T
  !-----------------------------------------------------------------------------
  found_L = .false.
  found_V = .false.
  tFV = tF
  tFL = tF
  count = 0

  do while ( .NOT. ( found_L .AND. found_V ) .AND. count < 8 )

     count = count + 1

     call state_tp ( tF, pF, xF, eta_start(1), rhoi_L, molar_rhoL, sL, hL, s_res=sres_L )
     rho_ideal_gas = pF / ( KBOL30 * tF )

     rho_ideal_gas = pF / ( KBOL30 * tF )
     eta_start(2) = min( rho_ideal_gas * z3t, 0.4_dp )    ! eta for ideal gas
     call state_tp ( tF, pF, xF, eta_start(2), rhoi_V, molar_rhoV, sV, hV, s_res=sres_V )

     if ( rhoi_L(1) > 3.0_dp * rho_ideal_gas ) then
        tFL = tF
        s_res_L = sres_L + log( pF / ( rhoi_L(1) * KBOL30 * tF ) )        ! convert to variables (p,T,x)
        found_L = .true.
     end if
     if ( rhoi_V(1) < 3.0_dp * rho_ideal_gas ) then
        tFV = tF
        s_res_V = sres_V + log( pF / ( rhoi_V(1) * KBOL30 * tF ) )        ! convert to variables (p,T,x)
        found_V = .true.
     end if
     if ( .NOT. found_V .AND. rhoi_V(1) > 3.0_dp * rho_ideal_gas ) then
        tF = tF * 1.5_dp
     end if
     if ( .NOT. found_L .AND. rhoi_L(1) < 3.0_dp * rho_ideal_gas ) then
        tF = tF * 0.6_dp
     end if
     ! write (*,*) 'T_L, T_V', tFL, tFV, found_L, found_V

  end do

  !-----------------------------------------------------------------------------
  ! use Clausius-Clapeyron equation to iterate boiling temperature
  !-----------------------------------------------------------------------------
  error = 1.0_dp
  damp = 1.0_dp

  tF = min( tFL, tFV )
  tlv = tF
  if ( outp >= 2 ) write (*,*) ' entering iteration of boiling temperature with T = ', tF
     count = 0

     do while ( error > tolerance .AND. count < max_iter )

        count = count + 1

        call pure_vapor_pressure ( tF, converg, plv, rhoL, rhoV )
        if ( outp >= 2 ) write (*,*) 'calcul. of vapor pressure for T = ', tF, ' converg = ', converg

        if ( converg ) then

           rhoi_L(1) = rhoL
           call state_trho ( tF, rhoi_L, pcalc, molar_rhoL, sL, hL )
           rhoi_V(1) = rhoV
           call state_trho ( tF, rhoi_V, pcalc, molar_rhoV, sV, hV )

           dlnp_d_inversT = - tF / plv * ( hV - hL ) / ( 1.0_dp/molar_rhoV - 1.0_dp/molar_rhoL )
           delta_inversT = ( log(pF) - log(plv) ) / dlnp_d_inversT
           error = abs( delta_inversT )

           t0 = tF
           tF = 1.0_dp / ( 1.0_dp / t0 + delta_inversT * damp )
           if (outp >= 1 ) write (*,*) 'Clausius-Clapeyron iter.', count, tF, delta_inversT

        else

           tc = tFL * 3.0_dp        ! starting value for critical point
           call critical_point_pure ( tc, converg_c, pc, rhoc )
           if ( .NOT. converg_c ) then
              write (*,*) 'pure_boiling_temperature: critical pt. not converged. Starting T=', tFL * 3.0_dp
              stop
           end if
           if ( outp >= 1 ) write (*,*) 'crit. pressure, tc, pc, pF:', tc, pc, pF
           if ( pF > pc ) then
              exit
           else
              if ( outp >= 1 ) write (*,*) ' difficult case (?): now reducing T'
              tF = 0.995_dp * tc
              damp = 0.9_dp
           end if
              
        end if

        if ( outp >= 3 ) read (*,*)

     end do

  !-----------------------------------------------------------------------------
  ! for convergence set tlv
  !-----------------------------------------------------------------------------
  if ( error < tolerance ) then
     converg = .true.
     tlv = tF
  end if

end subroutine pure_boiling_temperature



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_isotherm_or_isobar_diagram
!> \brief determine pure component properties along an isotherm or isobar
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_isotherm_or_isobar_diagram ( calcopt )

  use BASIC_VARIABLES, only: t_input, p_input
  use input_output, only: u_in_T, u_in_p, pure_output_state, unit_t, unit_p
  use properties, only: state_tp

  !-----------------------------------------------------------------------------
  integer, intent(in)                        :: calcopt

  !-----------------------------------------------------------------------------
  integer                                    :: i
  integer                                    :: steps
  real(dp)                                   :: tF, pF
  real(dp)                                   :: start_t, start_p
  real(dp)                                   :: end_t, end_p
  real(dp)                                   :: rhoL, rhoV
  real(dp)                                   :: hL, hV
  real(dp)                                   :: sL, sV
  real(dp)                                   :: cpL, cpV
  real(dp)                                   :: gL, gV
  real(dp)                                   :: molar_rho
  real(dp)                                   :: eta_start
  real(dp), dimension(ncomp)                 :: rhoi
  real(dp), dimension(ncomp)                 :: xF
  !-----------------------------------------------------------------------------

  steps = 40

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call pure_output_state ( tF, pF, rhoL, sL, hL, cp=cpL, header=.true. )

  !-----------------------------------------------------------------------------
  ! let user define the conditions for the isotherm or isobar
  !-----------------------------------------------------------------------------
  start_p = p_input
  start_t = t_input
  if ( calcopt == 2 ) then
     write (*,'(3a)') ' specify the end-pressure, in [',unit_p,']:'
     read (*,*) end_p
     end_p = end_p * u_in_p
     end_t = t_input
  else
     write (*,'(3a)') ' specify the end-temperature, in [',unit_t,']:'
     read (*,*) end_t
     end_t = end_t + u_in_T
     end_p = p_input
  end if

  !-----------------------------------------------------------------------------
  ! calc. the isotherm or isobar along a grid of 'steps' points
  !-----------------------------------------------------------------------------
  do i = 0, steps

     pF = start_p + ( end_p - start_p ) * real( i, KIND=dp ) / real( steps, KIND=dp )
     tF = start_t + ( end_t - start_t ) * real( i, KIND=dp ) / real( steps, KIND=dp )

     !--------------------------------------------------------------------------
     ! calc. state properties for a trial liquid and trial vapor
     !--------------------------------------------------------------------------
     xF(1) = 1.0_dp
     eta_start = 0.5_dp                        ! starting density for a liquid phase
     call state_tp ( tF, pF, xF, eta_start, rhoi, molar_rho, sL, hL, cp=cpL )
     rhoL = rhoi(1)

     eta_start = 1.E-5_dp                      ! starting density for a vapor phase
     call state_tp ( tF, pF, xF, eta_start, rhoi, molar_rho, sV, hV, cp=cpV )
     rhoV = rhoi(1)

     !--------------------------------------------------------------------------
     ! determine if liquid or vapor trial phase is stable
     !--------------------------------------------------------------------------
     gL = hL - tF * sL
     gV = hV - tF * sV

     !--------------------------------------------------------------------------
     ! output to file (and to terminal)
     !--------------------------------------------------------------------------
     if ( gL < gV ) then
        call pure_output_state ( tF, pF, rhoL, sL, hL, cp=cpL, header=.false. )
     else
        call pure_output_state ( tF, pF, rhoV, sV, hV, cp=cpV, header=.false. )
     end if

  end do

end subroutine pure_isotherm_or_isobar_diagram


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_critical_point_with_properties
!> \brief determine pure component crit. point and properties
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_critical_point_with_properties

  use BASIC_VARIABLES, only: t_input
  use input_output, only: pure_output_state
  use module_critical_point, only: critical_point_pure
  use properties, only: state_trho

  !-----------------------------------------------------------------------------
  real(dp)                                   :: tc, pc, rhoc
  real(dp)                                   :: hc
  real(dp)                                   :: sc
  real(dp)                                   :: molar_rho
  real(dp)                                   :: pcalc
  real(dp), dimension(ncomp)                 :: rhoi
  logical                                    :: converg
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call pure_output_state ( tc, pc, rhoc, sc, hc, header=.true. )

  tc = t_input * 3.0_dp       ! starting value for critical point

  !-----------------------------------------------------------------------------
  ! calculate critial point
  !-----------------------------------------------------------------------------
  call critical_point_pure ( tc, converg, pc, rhoc )

  if ( converg ) then

     !--------------------------------------------------------------------------
     ! calc. state properties at the critical point
     !--------------------------------------------------------------------------
     rhoi(1) = rhoc
     call state_trho ( tc, rhoi, pcalc, molar_rho, sc, hc )

     !--------------------------------------------------------------------------
     ! output to file (and to terminal)
     !--------------------------------------------------------------------------
     call pure_output_state ( tc, pc, rhoc, sc, hc, header=.false. )

     write (*,'(a,G18.9,a)') ' the critical temperature is: Tc= ', tc,' [K]'

  else

     write (*,*) 'calculation of critical point did not converge !'

  end if

end subroutine pure_critical_point_with_properties


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_virial_coefficients
!> \brief determine pure component virical coefficients (T-diagram)
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_virial_coefficients

  use PARAMETERS, only: NAV
  use BASIC_VARIABLES, only: t_input
  use input_output, only: u_in_T, unit_t
  use properties, only: calculate_pressure, ztot_rho

  !-----------------------------------------------------------------------------
  integer                                    :: i
  integer                                    :: steps
  real(dp)                                   :: tF
  real(dp)                                   :: start_t
  real(dp)                                   :: end_t
  real(dp)                                   :: rhoV
  real(dp)                                   :: pcalc, pcalc_z
  real(dp)                                   :: B_virial
  real(dp)                                   :: eta
  real(dp), dimension(ncomp)                 :: xF
  logical                                    :: is_open
  !-----------------------------------------------------------------------------

  steps = 40

  eta = 0.0_dp
  xF(1) = 1.0_dp

  !-----------------------------------------------------------------------------
  ! let user define the T-range
  !-----------------------------------------------------------------------------
  start_t = t_input
  write (*,'(3a)') ' specify the end-temperature, in [',unit_t,']:'
  read (*,*) end_t
  end_t = end_t + u_in_T

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  inquire( unit=40, opened=is_open )
  if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

  write  (*,'(3a)') '   T/',unit_t,'       B/(cm**3/mol)'
  write (40,'(3a)') '   T/',unit_t,'       B/(cm**3/mol)'

  !-----------------------------------------------------------------------------
  ! calc. second virial coeff. along a grid of 'steps' points
  !-----------------------------------------------------------------------------
  do i = 0, steps

     tF = start_t + ( end_t - start_t ) * real( i, KIND=dp ) / real( steps, KIND=dp )

     call calculate_pressure ( tF, xF, eta, rhoV, pcalc, pcalc_z )

     B_virial = ztot_rho / 1.E30_dp * NAV * 1.E6_dp

     write  (*,'(2(2x,f14.8))') tF-u_in_T, B_virial
     write (40,'(2(2x,f14.8))') tF-u_in_T, B_virial

  enddo

  if ( .NOT.is_open ) close (40)

end subroutine pure_virial_coefficients


end module module_pure_diagrams
