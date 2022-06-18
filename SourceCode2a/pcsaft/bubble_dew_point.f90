!> \file starting_value.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module bubble_dew_point

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                              :: outp = 0
  real(dp)                                        :: t_transfer, p_transfer
  real(dp), dimension(2)                          :: eta_start_transfer
  real(dp), dimension(:), allocatable             :: x1_transfer
  real(dp), dimension(:), allocatable             :: rhoi1_transfer, rhoi2_transfer
  logical                                         :: iterate_t_tranfer

  private
  public :: bubble_point_rachford_rice

contains

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine bubble_point_rachford_rice
!
! This subroutine performs a bubble- or dew-point calculation for a defined
! composition (xi1(...)). A bubble point calculation for a given liquid composition
! is done, when eta_start(1) is a liquid density (starting value, needs to be
! possible to converge). A dew point calculation is done, when eta_start(1) (low
! starting value) converges to a vapor density.
!
! If eta_start(1) is for a liquid and eta_start(2) is for a vapor, i.e. if a
! VLE is considered, the starting value for the vapor phase is here generated.
! Otherwise, a suitable starting value for xi2 must be an input !!
! 
! Either, the temperature or the pressure is iterated. The algorith is based on
! the Rachford-Rice equation.
!
! input:   xi1(:),  t or p
!          And starting values: eta_start(1), eta_start(2), p or t
!          if not VLE considered: also provide starting value for xi2 !
!
! output:  p or t, xi2(:),  rhoi1(:), rhoi1(:)
!          if converg=.true. : lnphi(:,:) are also converged values (but no output).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine bubble_point_rachford_rice ( iterate_t, t, p, eta_start, xi1, xi2, rhoi1, rhoi2, get_start_val, converg )

  use PARAMETERS, only: machine_eps, PI_6
  use pcsaft_pure_and_binary_parameters, only: sigma, mseg
  use stability, only: get_eta
  use properties, only: chemical_potential_tp, chemical_potential_tp_derivative_tp
  implicit none

  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  real(dp), intent(inout)                         :: t
  real(dp), intent(inout)                         :: p
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), dimension(ncomp), intent(in)          :: xi1
  real(dp), dimension(ncomp), intent(inout)       :: xi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  logical, intent(in)                             :: get_start_val
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer                                         :: k1, k2
  integer                                         :: nr_inside_steps, nr_outside_steps
  integer, parameter                              :: prog_no = 10
  real(dp), parameter                             :: tol_out = 1.E-11_dp
  real(dp), parameter                             :: tol_in_p = 1.E-10_dp
  real(dp), parameter                             :: tol_in_t = 1.E-8_dp
  real(dp), parameter                             :: damping = 0.6_dp
  real(dp), dimension(2)                          :: eta_start_var
  real(dp)                                        :: tol_in
  real(dp)                                        :: p0, t0
  real(dp)                                        :: d_p, d_t
  real(dp)                                        :: deltap, deltat
  real(dp)                                        :: p_sav, t_sav
  real(dp), dimension(ncomp)                      :: xi2_save
  real(dp), dimension(ncomp)                      :: lnphi_1, lnphi_2
  real(dp), dimension(ncomp)                      :: Ki
  real(dp), dimension(ncomp)                      :: xi_compare
  real(dp), dimension(prog_no)                    :: progress
  real(dp)                                        :: f0, f1, dfdp
  real(dp)                                        :: error1, error2
  real(dp)                                        :: similarity
  real(dp)                                        :: stepsize
  logical                                         :: dew
  logical                                         :: SR_exit
  logical                                         :: trivial
  logical                                         :: no_progress
  logical                                         :: always_damping
  !-----------------------------------------------------------------------------

  t_sav = t
  p_sav = p

  converg = .false.

  trivial = .false.
  always_damping = .false.
  SR_exit = .false.
  no_progress = .false.

  eta_start_var(1:2) = eta_start(1:2)

  !-----------------------------------------------------------------------------
  ! bubble point or dew point
  !-----------------------------------------------------------------------------

  if ( eta_start_var(1) < 0.25_dp .AND. eta_start_var(2) > 0.25_dp ) then
     dew = .true.
  else
     dew = .false.
  end if

  !-----------------------------------------------------------------------------
  ! defining some tolerances and maximum number of iterations
  !-----------------------------------------------------------------------------

  if ( ncomp == 2 ) then
     nr_outside_steps = 100
  else
     nr_outside_steps = 400
  end if
  nr_inside_steps = 5

  if ( iterate_t ) then
     tol_in = tol_in_t
  else
     tol_in = tol_in_p
  end if

  d_p = - 1.0_dp                    ! pressure step in unit [Pa]
  d_t = - 1.E-3_dp                  ! temperature step in unit [K]

  progress(:) = 0.0_dp
  
  !-----------------------------------------------------------------------------
  ! for specific case of VLE (phase 1 as liquid, phase 2 as vapor):
  ! generate starting values for Rachford-Rice iteration
  ! liquid with xi1=xiF(:); vapor xi2 assuming ideal gas (phi_2(i) = 1)
  !-----------------------------------------------------------------------------
  if ( eta_start_var(1) > 0.25_dp .AND. eta_start_var(2) < 0.25_dp .AND. get_start_val ) then

     stepsize = 0.05_dp
     call determine_starting_value_bubble_point ( iterate_t, t, p,  &
                            eta_start_var, stepsize, xi1, xi2, rhoi1, rhoi2, lnphi_1, lnphi_2 )

  end if

  !-----------------------------------------------------------------------------
  ! for specific case of VLE (phase 1 as vapor, phase 2 as liquid):
  ! generate starting values for Rachford-Rice iteration
  ! vapor with xi1=xiF(:); liquid xi2 by calculating chem. potential of trial phase
  !-----------------------------------------------------------------------------
  if ( dew .AND. get_start_val ) then

     xi2(:) = xi1(:)                          ! initial guess for xi2:  xi2 = xi1
     stepsize = 0.02_dp
     call determine_starting_value_dew_point ( iterate_t, t, p,  &
                            eta_start_var, stepsize, xi1, xi2, rhoi1, rhoi2, lnphi_1, lnphi_2 )

  end if

  xi2_save(:) = xi2(:)
  xi_compare(:) = xi2(:)

  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'bbp-RR starting values: t, p, x1(1), x2(1)', t, p, xi1(1), xi2(1)
  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'bbp-RR,   N_out, N_inner, xi1(1),   xi2(1),   p,   t,   error1'

  !-----------------------------------------------------------------------------
  ! start iteration
  !-----------------------------------------------------------------------------

  k1 = 0
  error1 = tol_out + 1.0_dp
  do while ( error1 > tol_out .AND. k1 < nr_outside_steps )

     !--------------------------------------------------------------------------
     ! outer loop iteration: converging the mole fractions
     !--------------------------------------------------------------------------

     k1 = k1 + 1

     k2 = 0
     error2 = tol_in + 1.0_dp
     do while ( error2 > tol_in .AND. k2 < nr_inside_steps )

        !=======================================================================
        ! inner loop iteration: converging bubble/dew point T or p
        !=======================================================================

        k2 = k2 + 1

        !--------------------------------------------------------------------------
        ! perform numerical derivative of chemical potential to T or to p.
        ! Numerical derivative is more efficient than analytic derivative
        !--------------------------------------------------------------------------
        p0 = p
        t0 = t
        if ( iterate_t ) then
           !call chemical_potential_tp_derivative_tp ( t, p, xi1, eta_start_var(1), rhoi1, lnphi_1, lnphi_t=lnphi_t_1 )
           !call chemical_potential_tp_derivative_tp ( t, p, xi2, eta_start_var(2), rhoi2, lnphi_2, lnphi_t=lnphi_t_2 )
           !dfdp = SUM( xi1(1:ncomp) * EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) ) * ( lnphi_t_1(1:ncomp) - lnphi_t_2(1:ncomp) ) )
           t = t0 - d_t
        else
           !call chemical_potential_tp_derivative_tp ( t, p, xi1, eta_start_var(1), rhoi1, lnphi_1, lnphi_p=lnphi_p_1 )
           !call chemical_potential_tp_derivative_tp ( t, p, xi2, eta_start_var(2), rhoi2, lnphi_2, lnphi_p=lnphi_p_2 )
           !dfdp = SUM( xi1(1:ncomp) * EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) ) * ( lnphi_p_1(1:ncomp) - lnphi_p_2(1:ncomp) ) )
           !dfdp = dfdp * p0    ! convert to derivative d( sum(x(:)*Ki(:)) ) / dln(p) instead of d(...) / d(p)
           p = p0 - d_p
        end if

        call chemical_potential_tp ( t, p, xi1, eta_start_var(1), rhoi1, lnphi_1 )
        call chemical_potential_tp ( t, p, xi2, eta_start_var(2), rhoi2, lnphi_2 )

        if ( dew .AND. abs(get_eta(rhoi2)/get_eta(rhoi1)-1.0_dp) < 0.3_dp ) then
           stepsize = 0.005_dp
           call determine_starting_value_dew_point ( iterate_t, t, p,  &
                eta_start_var, stepsize, xi1, xi2, rhoi1, rhoi2, lnphi_1, lnphi_2 )
           if ( iterate_t ) then
              t0 = t + d_t
           else
              p0 = p + d_p
           end if
        end if

        Ki(1:ncomp) = EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

        f0 = SUM( xi1(1:ncomp) * Ki(1:ncomp) ) - 1.0_dp

        eta_start_var(1) = get_eta( rhoi1 )
        eta_start_var(2) = get_eta( rhoi2 )

        if ( ( eta_start_var(1) < 0.13_dp .AND. eta_start_var(2) < 0.13_dp ) .OR.  &
             ( eta_start_var(1) > 0.15_dp .AND. eta_start_var(2) > 0.15_dp ) ) then
           eta_start_var(:) = eta_start(:)
        end if

        p = p0
        t = t0
        call chemical_potential_tp ( t, p, xi1, eta_start_var(1), rhoi1, lnphi_1 )
        call chemical_potential_tp ( t, p, xi2, eta_start_var(2), rhoi2, lnphi_2 )

        Ki(1:ncomp) = EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

        f1 = SUM( xi1(1:ncomp) * Ki(1:ncomp) ) - 1.0_dp

        if ( ABS(1.0_dp-sum(rhoi1(:))/sum(rhoi2(:))) < 1.E-5_dp  &
          .AND. SUM(ABS(Ki(1:ncomp)-1.0_dp)) < 1.E-6_dp ) exit

        if ( iterate_t ) then

           dfdp = ( f1 - f0 ) / d_t
           if ( abs( dfdp ) < machine_eps ) dfdp = 0.0000001_dp

           !--------------------------------------------------------------------
           ! Newton inner loop T-iteration
           !--------------------------------------------------------------------
           deltat = f1 / dfdp
           if ( deltat >  20.0_dp ) deltat =  20.0_dp   ! 20 K as max. T-step
           if ( deltat < -20.0_dp ) deltat = -20.0_dp
           if ( ( t0 - deltat ) < 0.0_dp ) deltat = t0 - 1.E-2_dp
           t = t0 - deltat

           !--------------------------------------------------------------------
           ! error for inner loop T-iteration
           !--------------------------------------------------------------------
           error2 = ABS( t - t0 )
           if ( outp >= 3 ) write (*,'(a,3G20.12)') 'inner loop error', t, t0, error2

        else

           dfdp = ( f1 - f0 ) / ( LOG( p0 ) - LOG( p0 - d_p ) )
           if ( abs( dfdp ) < machine_eps ) dfdp = 0.0000001_dp

           !--------------------------------------------------------------------
           ! Newton inner loop ln(p) iteration
           !--------------------------------------------------------------------
           ! p = EXP( LOG(p0) - f1 / dfdp  )
           deltap = f1 / dfdp
           if ( deltap >  2.0_dp ) deltap =  2.0_dp   ! exp(2) = 7.38 as a max. factor for p-step in non-logarithmic unit
           if ( deltap < -2.0_dp ) deltap = -2.0_dp
           p = LOG( p0 ) - deltap
           p = EXP( p )
           p = min( p, 1.E9_dp )

           !--------------------------------------------------------------------
           ! damping inner loop ln(p) iteration
           !--------------------------------------------------------------------
           if ( ABS( p / p0 - 1.0_dp ) > 0.3_dp .OR. always_damping ) p = damping * p + ( 1.0_dp - damping ) * p0

           !--------------------------------------------------------------------
           ! error inner loop ln(p) iteration
           !--------------------------------------------------------------------
           error2 = ABS( p / p0 - 1.0_dp ) * 0.1_dp

           if ( outp >= 3 ) write (*,'(a,3G20.12)') 'inner loop error', p, p0, error2
           if ( error2 > 10.0_dp ) SR_exit = .true.
           if ( error2 > 10.0_dp ) exit

        end if

     end do

     if ( SR_exit ) exit

     !--------------------------------------------------------------------------
     ! determine x and error, outer loop iteration
     !--------------------------------------------------------------------------
     xi2(1:ncomp) = Ki(1:ncomp)* xi1(1:ncomp) / SUM( xi1(1:ncomp) * Ki(1:ncomp) )
     xi2(1:ncomp) = max( xi2(1:ncomp), 1.E-50_dp )

     if ( sum( xi2(1:ncomp) ) < 0.4_dp ) xi2(1:ncomp) = xi2(1:ncomp) / sum( xi2(1:ncomp) )

     error1 =  SUM( ABS( xi_compare(1:ncomp) - xi2(1:ncomp) ) )

     !--------------------------------------------------------------------------
     ! check for trivial solution
     !--------------------------------------------------------------------------
     similarity = sum( ABS( rhoi2(1:ncomp)/(rhoi1(1:ncomp)+1.E-12_dp) - 1.0_dp )  &
                       * PI_6 *mseg(1:ncomp)*sigma(1:ncomp)**3 )
     ! similarity = sum( ABS( rhoi2(1:ncomp) - rhoi1(1:ncomp) ) * PI_6 *mseg(1:ncomp)*sigma(1:ncomp)**3 )

     if ( similarity < 1.0_dp ) trivial = .true.

     !--------------------------------------------------------------------------
     ! for slow convergence: accelerate convergence. This should soon be modified
     ! to be a full Anderson mixing scheme
     !--------------------------------------------------------------------------

     !--------------------------------------------------------------------------
     ! error outer iteration
     !--------------------------------------------------------------------------
     xi_compare(:) = xi2(:)

     if ( error2 > tol_in ) then
        always_damping = .true.
        error1 = error1 + tol_out * 1000.0_dp
     end if

     !--------------------------------------------------------------------------
     ! for slow convergence: relax the tolerance somewhat
     !--------------------------------------------------------------------------
     if ( k1 >= 50 ) error1 = error1 / 10.0_dp
     if ( k1 >= 90 ) error1 = error1 / 100.0_dp

     if ( outp >= 2 ) write (*,'(a,2i5,5G18.10)') 'bbp-RR',k1, k2, xi1(1), xi2(1), p, t, error1

     !--------------------------------------------------------------------------
     ! monitor progress
     !--------------------------------------------------------------------------
     progress(2:prog_no) = progress(1:prog_no-1)
     progress(1) = error1

     if ( progress(1) >= progress(prog_no) / 1.001_dp .AND. k1 > ( prog_no + 2 ) ) no_progress = .true.

     !--------------------------------------------------------------------------
     ! for slow convergence: switch to direct solution procedure
     !--------------------------------------------------------------------------
     if ( ( k1 == 20 .AND. error1 < 1.E-2 .AND. error1 > 1.E-7 .AND. ncomp <= 2 ) .OR.  &
            ( k1 == 50 .AND. error1 < 1.E-3 ) .OR. k1 == 100 .OR. trivial .OR. no_progress ) then

        if ( trivial ) then
           xi2(:) = xi2_save(:)
           t = t_sav
           p = p_sav
        end if

        if ( outp >= 2 ) write (*,'(a,4G20.12)') 'bbp-RR starting values: t, p, x1(1), x2(1)', t, p, xi1(1), xi2(1)

        call two_phase_equilibrium ( iterate_t, t, p, xi1, eta_start, xi2, rhoi1, rhoi2, converg )

        if ( converg ) then
           error1 = 0.0_dp
           error2 = 0.0_dp
           trivial = .false.
           SR_exit = .false.
           no_progress = .false.
        else
           progress(:) = 10.0_dp
           no_progress = .false.
        end if

     end if

     !--------------------------------------------------------------------------
     ! for no progress or for trivial solution: exit
     !--------------------------------------------------------------------------
     if ( trivial .OR. no_progress ) then
        if ( outp > 0 ) write (*,*) 'bubble-point RR no progress. Conventional calc. failed.', trivial, no_progress
        exit
     end if

  end do


  !-----------------------------------------------------------------------------
  ! if convergence: accept solution
  !-----------------------------------------------------------------------------

  if ( error1 < tol_out .AND. error2 < tol_in ) then

     converg = .true.
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'bubble point Rachford-Rice converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  T, p   = ',t, p
     if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) 'eta 1,2 = ', get_eta( rhoi1 ), get_eta( rhoi2 )
     if ( outp >= 2 ) write (*,*) ' '

  else

     t = t_sav
     p = p_sav
     if ( outp >= 1 ) write (*,*) 'bbp-Rachf.-R.: not converged'

  end if

end subroutine bubble_point_rachford_rice


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine determine_starting_value_bubble_point
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine determine_starting_value_bubble_point ( iterate_t, t, p, eta_start, stepsize,  &
     xi1, xi2, rhoi1, rhoi2, lnphi_1, lnphi_2 )

  use stability, only: get_eta
  use properties, only: chemical_potential_tp, calculate_density
  implicit none

  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  real(dp), intent(inout)                         :: t
  real(dp), intent(inout)                         :: p
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), intent(in)                            :: stepsize
  real(dp), dimension(ncomp), intent(in)          :: xi1
  real(dp), dimension(ncomp), intent(out)         :: xi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  real(dp), dimension(ncomp), intent(out)         :: lnphi_1
  real(dp), dimension(ncomp), intent(out)         :: lnphi_2

  !-----------------------------------------------------------------------------
  real(dp)                                        :: eta1, eta2
  real(dp)                                        :: eta_target
  real(dp), dimension(ncomp)                      :: xi2_ideal_gas
  logical                                         :: towards_vapor
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! chemical potential and density of liquid for (T,p,xi1)-starting value
  !-----------------------------------------------------------------------------
  call chemical_potential_tp ( t, p, xi1, eta_start(1), rhoi1, lnphi_1 )
  eta1 = get_eta( rhoi1 )

  !-----------------------------------------------------------------------------
  ! if liquid does not exist at (T,p,x)-starting values, then modify T or p
  !-----------------------------------------------------------------------------
  eta_target = 0.12_dp           ! minimum acceptable value
  if ( eta1 < eta_target ) then

     towards_vapor = .false.
     call modify_starting_condition_for_VLE ( iterate_t, towards_vapor, stepsize, xi1, eta_start(1), eta_target, t, p )
     call chemical_potential_tp ( t, p, xi1, eta_start(1), rhoi1, lnphi_1 )

  end if

  !-----------------------------------------------------------------------------
  ! determine initial guess for xi2 of vapor, assuming ideal gas (phi_2(i) = 1)
  !-----------------------------------------------------------------------------
  xi2( 1:ncomp ) = xi1( 1:ncomp ) * exp( lnphi_1(1:ncomp) )
  xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )
  xi2_ideal_gas( 1:ncomp ) = xi2( 1:ncomp )

  if ( outp >= 3 ) write (*,'(a,4G20.12)') 'id.Gas starting values: t, p, x1(1), x2(1)', t, p, xi1(1), xi2(1)
  ! call calculate_density ( t, p, xi2, eta_start(2), eta2, rho2 )

  !-----------------------------------------------------------------------------
  ! refine estimate for xi2 (vapor) by considering fug. coefficient of gas phase
  !-----------------------------------------------------------------------------
  call chemical_potential_tp ( t, p, xi2, eta_start(2), rhoi2, lnphi_2 )
  eta2 = get_eta( rhoi2 )
  xi2(:) = xi2_ideal_gas(:) / exp( lnphi_2(1:ncomp) )
  xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )

  !-----------------------------------------------------------------------------
  ! close to critical conditions, repeat the procedure to refine estimate
  !-----------------------------------------------------------------------------
  if ( abs( eta1/eta2 - 1.0_dp ) < 5.0_dp ) then
     call chemical_potential_tp ( t, p, xi2, eta_start(2), rhoi2, lnphi_2 )
     eta2 = get_eta( rhoi2 )
     xi2(:) = xi2_ideal_gas(:) / exp( lnphi_2(1:ncomp) )
     xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )
     call chemical_potential_tp ( t, p, xi2, eta_start(2), rhoi2, lnphi_2 )
     eta2 = get_eta( rhoi2 )
     xi2(:) = xi2_ideal_gas(:) / exp( lnphi_2(1:ncomp) )
     xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )
  end if

  !-----------------------------------------------------------------------------
  ! if vapor does not exist at (T,p,xi2)-values, then modify T or p
  !-----------------------------------------------------------------------------
  eta_target = 0.25_dp           ! maximum acceptable value
  if ( eta2 > eta_target ) then

     towards_vapor = .true.
     call modify_starting_condition_for_VLE ( iterate_t, towards_vapor, stepsize, xi2, eta_start(2), eta_target, t, p )

  end if

  call chemical_potential_tp ( t, p, xi2, eta_start(2), rhoi2, lnphi_2 )
  eta2 = get_eta( rhoi2 )

  !-----------------------------------------------------------------------------
  ! a (not conclusive) indication that a classical VLE-solution will not be found
  !-----------------------------------------------------------------------------
  if ( abs( eta1/eta2 - 1.0_dp ) < 0.8 .AND. outp > 0 ) write (*,*) 'probably a VLE can not be found'

end subroutine determine_starting_value_bubble_point


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine determine_starting_value_dew_point
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine determine_starting_value_dew_point ( iterate_t, t, p, eta_start, stepsize,  &
     xi1, xi2, rhoi1, rhoi2, lnphi_1, lnphi_2 )

  use stability, only: get_eta
  use properties, only: chemical_potential_tp, calculate_density
  implicit none

  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  real(dp), intent(inout)                         :: t
  real(dp), intent(inout)                         :: p
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), intent(in)                            :: stepsize
  real(dp), dimension(ncomp), intent(in)          :: xi1
  real(dp), dimension(ncomp), intent(inout)       :: xi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  real(dp), dimension(ncomp), intent(out)         :: lnphi_1
  real(dp), dimension(ncomp), intent(out)         :: lnphi_2

  !-----------------------------------------------------------------------------
  integer                                         :: count
  real(dp)                                        :: eta1, eta2
  real(dp)                                        :: eta_target
  integer, parameter                              :: max_steps = 15
  logical                                         :: towards_vapor
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! chemical potential and density of vapor for (T,p,xi1)-starting value
  !-----------------------------------------------------------------------------
  call chemical_potential_tp ( t, p, xi1, eta_start(1), rhoi1, lnphi_1 )
  eta1 = get_eta( rhoi1 )
  call chemical_potential_tp ( t, p, xi2, eta_start(1), rhoi2, lnphi_2 )
  eta2 = get_eta( rhoi2 )

     eta_target = 0.12_dp           ! maximum acceptable value
  !-----------------------------------------------------------------------------
  ! if vapor and liquid are not enough different at (T,p,x), modify T or p
  !-----------------------------------------------------------------------------
  count = 0

  do while ( abs( eta2 / eta1 - 1.0_dp ) < 0.3_dp .AND. count < max_steps )

     count = count + 1

     if ( eta1 > eta_target ) then

        !--------------------------------------------------------------------------
        ! find vapor phase (by increasing t or decreasing p)
        !--------------------------------------------------------------------------
        towards_vapor = .true.
        call modify_starting_condition_for_VLE ( iterate_t, towards_vapor, stepsize, xi1, eta_start(1), eta_target, t, p )
        call chemical_potential_tp ( t, p, xi1, eta_start(1), rhoi1, lnphi_1 )
        eta1 = get_eta( rhoi1 )
        if ( outp > 0 ) write (*,*) 'find vapor: t, p', t, p

        xi2(:) = xi1(:) * exp( lnphi_1(1:ncomp) ) / exp( lnphi_2(1:ncomp) )
        xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )
        call chemical_potential_tp ( t, p, xi2, eta_start(1), rhoi2, lnphi_2 )
        eta2 = get_eta( rhoi2 )
        if ( outp > 1 ) write (*,*) 'find vapor: eta1, eta2', eta1, eta2

     else

        !--------------------------------------------------------------------------
        ! find liquid phase (by decreasing t or increasing p)
        !--------------------------------------------------------------------------
        towards_vapor = .false.
        call modify_starting_condition_for_VLE ( iterate_t, towards_vapor, stepsize, xi2, eta_start(2), eta_target, t, p )
        call chemical_potential_tp ( t, p, xi2, eta_start(2), rhoi2, lnphi_2 )
        eta2 = get_eta( rhoi2 )
        if ( outp > 0 ) write (*,*) 'find liquid: t, p', t, p

        xi2(:) = xi1(:) * exp( lnphi_1(1:ncomp) ) / exp( lnphi_2(1:ncomp) )
        xi2( 1:ncomp ) = xi2( 1:ncomp ) / sum( xi2( 1:ncomp ) )

        call chemical_potential_tp ( t, p, xi1, eta_start(1), rhoi1, lnphi_1 )
        eta1 = get_eta( rhoi1 )
        if ( outp > 1 ) write (*,*) 'find liquid: eta1, eta2', eta1, eta2

     end if

  end do

end subroutine determine_starting_value_dew_point


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine modify_starting_condition_for_VLE
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine modify_starting_condition_for_VLE ( iterate_t, towards_vapor, stepsize, xi, eta_start, eta_target, t, p )

  use properties, only: chemical_potential_tp, calculate_density
  implicit none

  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  logical, intent(in)                             :: towards_vapor
  real(dp), intent(in)                            :: stepsize
  real(dp), dimension(ncomp), intent(in)          :: xi
  real(dp), intent(in)                            :: eta_start
  real(dp), intent(in)                            :: eta_target
  real(dp), intent(inout)                         :: t
  real(dp), intent(inout)                         :: p

  !-----------------------------------------------------------------------------
  integer                                         :: count
  integer, parameter                              :: max_steps = 15
  real(dp)                                        :: def_factor_p
  real(dp)                                        :: def_factor_t
  real(dp)                                        :: factor_t, factor_p
  real(dp)                                        :: eta, rho
  logical                                         :: repeat
  logical                                         :: success
  !-----------------------------------------------------------------------------

  repeat = .true.
  success = .false.

  def_factor_t = 1.0_dp + stepsize
  def_factor_p = 1.0_dp + 3.0_dp * stepsize

  count = 0

  !-----------------------------------------------------------------------------
  ! if towards_vapor=.true.  increase T (iterate_t) or decrease p (.NOT.iterate_t)
  ! if towards_vapor=.false. decrease T (iterate_t) or increase p (.NOT.iterate_t)
  !-----------------------------------------------------------------------------
  if ( iterate_t ) then

     factor_p = 1.0_dp

     if ( towards_vapor ) then
        factor_t = def_factor_t
     else
        factor_t = 1.0_dp / def_factor_t
     end if

  else

     if ( towards_vapor ) then
        factor_p = 1.0_dp / def_factor_p
     else
        factor_p = def_factor_p
     end if

     factor_t = 1.0_dp

  end if
  
  !-----------------------------------------------------------------------------
  ! modify T or p to approach vapor, or approach liquid (depending on towards_vapor)
  !-----------------------------------------------------------------------------
  do while ( repeat )

     count = count + 1

     t = t * factor_t
     p = p * factor_p

     !----- get density for modified T or p, see if eta_target is reached -----
     call calculate_density ( t, p, xi, eta_start, eta, rho )

     ! write (*,'(i3,2G15.5,F15.2)') count, eta, t, p

     !----- accelerate modifying T or p ----------------------------------------
     if ( count == 10 ) then
        factor_t = factor_t**3
        factor_p = factor_p**3
     end if

     !----- exit criterion -----------------------------------------------------
     if ( ( .NOT.towards_vapor .AND. eta > eta_target ) .OR.  &
               ( towards_vapor .AND. eta < eta_target ) ) then
        success = .true.
     end if

     if ( count >= max_steps .OR. success )  repeat = .false.

  end do

end subroutine modify_starting_condition_for_VLE



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine two_phase_equilibrium
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine two_phase_equilibrium ( iterate_t, tF, pF, x1, eta_start, x2, rhoi1, rhoi2, converg )

  use Solve_NonLin, only: hbrd
  implicit none
  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  real(dp), intent(inout)                         :: tF
  real(dp), intent(inout)                         :: pF
  real(dp), dimension(ncomp), intent(in)          :: x1
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), dimension(ncomp), intent(inout)       :: x2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
!!$  interface
!!$     subroutine bin_3_phase_objec_fct ( ndim, y, residu, dummy )
!!$       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision       
!!$       integer, intent(in)                        :: ndim
!!$       real(dp), dimension(ndim), intent(in)      :: y
!!$       real(dp), dimension(ndim), intent(out)     :: residu
!!$       integer, intent(in out)                    :: dummy
!!$     end subroutine bin_3_phase_objec_fct
!!$  end interface

  integer                                         :: i
  integer                                         :: info
  integer                                         :: ndim
  real(dp), allocatable, dimension(:)             :: y
  real(dp), allocatable, dimension(:)             :: diag
  real(dp), allocatable, dimension(:)             :: residu
  real(dp)                                        :: tolerance
  real(dp)                                        :: step_size
  real(dp)                                        :: t_save, p_save
  real(dp), dimension(ncomp)                      :: x2_save
  !-----------------------------------------------------------------------------

  info = 1

  t_save = tF
  p_save = pF
  x2_save(:) = x2(:)

  ndim = ncomp + 1
  do i = 1, ncomp
     if ( x1(i) < 1.E-100_dp ) then
        ndim = ndim - 1
        write (*,*) ' two_phase_equilibrium: extend SR for case where xi1(i) can be zero'
        return
     end if
  end do
  
  allocate( y(ndim), diag(ndim), residu(ndim) )
  allocate( x1_transfer(ncomp), rhoi1_transfer(ncomp), rhoi2_transfer(ncomp) )

  tolerance = 1.E-9_dp
  step_size = 1.E-8_dp

  converg = .false.

  !-----------------------------------------------------------------------------
  ! assign starting values for the iterated variables (to vector y)
  !-----------------------------------------------------------------------------

  y(1:ncomp) = log( x2(1:ncomp) )
  if ( iterate_t ) then
     y(ncomp+1) = log( tF )
  else
     y(ncomp+1) = log( pF / 1000.0_dp )
  end if
  
  !-----------------------------------------------------------------------------
  ! solve the phase equilibrium conditions
  !-----------------------------------------------------------------------------

  t_transfer = tF
  p_transfer = pF
  iterate_t_tranfer = iterate_t             ! transfer iterate_t to the objective fct
  x1_transfer(1:ncomp) = x1(1:ncomp)
  eta_start_transfer(1:2) = eta_start(1:2)

  CALL hbrd ( two_phase_objec_fct, ndim, y, residu, step_size, tolerance, info, diag, N_max_itr=30 )


  !-----------------------------------------------------------------------------
  ! determine, whether the solution is acceptable (if yes: convergence = .true.)
  !-----------------------------------------------------------------------------

  if ( sum( abs( residu( 1:ndim ) ) ) < tolerance * 10.0_dp .OR. info == 1 ) converg = .true.

  if ( converg ) then

     if ( iterate_t ) then
        tF = exp( y(ncomp+1) )
     else
        pF = exp( y(ncomp+1) ) * 1000.0_dp
     end if
     x2(1:ncomp) = exp( y(1:ncomp) )
     x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )

     rhoi1(1:ncomp) = rhoi1_transfer(1:ncomp)
     rhoi2(1:ncomp) = rhoi2_transfer(1:ncomp)

     if ( outp >= 1 ) write (*,*) 'two_phase_equilibrium converged', info, sum( abs( residu( 1:ndim ) ) )

  else

     tF = t_save
     pF = p_save
     x2(:) = x2_save(:)
     ! write (*,*) 'error two_phase_equilibrium', info, sum( abs( residu( 1:ndim ) ) )
     if ( outp >= 1 ) write (*,*) 'bubble point conventional two_phase solver not converged', info

  end if

  deallocate( y, diag, residu )
  deallocate( x1_transfer, rhoi1_transfer, rhoi2_transfer )

end subroutine two_phase_equilibrium



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE two_phase_objec_fct
!
! equations solved for 2-phase equilibrium
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine two_phase_objec_fct ( ndim, y, residu, dummy )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: sigma, mseg
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in)           :: y
  real(dp), dimension(ndim), intent(out)          :: residu
  integer, intent(in out)                         :: dummy

  !-----------------------------------------------------------------------------
  real(dp)                                        :: tF
  real(dp)                                        :: pF
  real(dp), dimension(ncomp)                      :: x1
  real(dp), dimension(ncomp)                      :: x2
  real(dp), dimension(ncomp)                      :: rhoi1
  real(dp), dimension(ncomp)                      :: rhoi2
  real(dp), dimension(ncomp)                      :: lnphi_1
  real(dp), dimension(ncomp)                      :: lnphi_2
  real(dp), dimension(2)                          :: eta_start
  real(dp)                                        :: similarity
  logical                                         :: iterate_t
  !-----------------------------------------------------------------------------

  iterate_t = iterate_t_tranfer
  x1(1:ncomp) = x1_transfer(1:ncomp)

  !-----------------------------------------------------------------------------
  ! assign values for the iterated quantities
  !-----------------------------------------------------------------------------

  x2(1:ncomp) = exp( y(1:ncomp) )
  if ( iterate_t ) then
     tF = exp( y(ncomp+1) )
     pF = p_transfer
  else
     pF = exp( y(ncomp+1) ) * 1000.0_dp
     tF = t_transfer
  end if

  !-----------------------------------------------------------------------------
  ! determine first 3 values of residuum-vector
  !-----------------------------------------------------------------------------

  residu(ndim) = sum( x2(1:ncomp) ) - 1.0_dp

  !-----------------------------------------------------------------------------
  ! ensure summation relation
  !-----------------------------------------------------------------------------
  x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! standard starting values for densities
  !-----------------------------------------------------------------------------

  eta_start(1:2) = eta_start_transfer(1:2)

  !-----------------------------------------------------------------------------
  ! calculate chemical potentials
  !-----------------------------------------------------------------------------

  call chemical_potential_tp ( tF, pF, x1, eta_start(1), rhoi1(:), lnphi_1(:) )
  call chemical_potential_tp ( tF, pF, x2, eta_start(2), rhoi2(:), lnphi_2(:) )

  rhoi1_transfer(1:ncomp) = rhoi1(1:ncomp)
  rhoi2_transfer(1:ncomp) = rhoi2(1:ncomp)

  !-----------------------------------------------------------------------------
  ! check for trivial solution
  !-----------------------------------------------------------------------------
  similarity = sum( ABS( rhoi2(1:ncomp)/(rhoi1(1:ncomp)+1.E-12_dp) - 1.0_dp ) * PI_6 *mseg(1:ncomp)*sigma(1:ncomp)**3 )

  if ( similarity < 1.0_dp ) then
     residu(ndim) = residu(ndim) + sign( ( ( 1.0_dp - similarity ) * 10.0_dp ), residu(ndim) )
     if ( outp > 1 ) write (*,*) 'adding penalty function', residu(ndim), ( ( 1.0_dp - similarity ) * 10.0_dp )
  end if

  !-----------------------------------------------------------------------------
  ! determine remaining values of residuum-vector (isofugacity conditions)
  !-----------------------------------------------------------------------------

  residu(1:ncomp) = log( x1(1:ncomp) ) + lnphi_1(1:ncomp) - log( x2(1:ncomp) ) - lnphi_2(1:ncomp)

  if ( outp > 0 .AND. dummy == 1 ) write (*,'(a,4G18.8)') 'two_phase_objec_fct',tF, pF, sum( abs( residu ) ), similarity
  if ( outp > 1 .AND. dummy /= 1 ) write (*,'(a,3G18.8)') 'two_phase_objec_fct',tF, pF, sum( abs( residu ) )

end subroutine two_phase_objec_fct


end module bubble_dew_point
