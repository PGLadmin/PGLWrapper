!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine DENSITY_ITERATION
!
!> \brief density iteration
!!
!! iterates the density until the calculated pressure 'pcalc' is equal to
!! the specified pressure 'p'. \n  A Newton-scheme is used for determining
!! the root to the objective function \f$ f(eta) = (p_{ges} / p ) - 1.0. \f$
!! The starting value for the packing fraction (= dimensionless density) is
!! eta_start.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine density_iteration ( eta_start, p, eta, pcalc, pcalc_z )

  use PARAMETERS, only: dp, machine_eps
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(in)                       :: p
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z

  !-----------------------------------------------------------------------------
  integer, parameter                         :: outp = 0
  integer, parameter                         :: max_iterations = 50
  real(dp), parameter                        :: density_tolerance = 1.E-12

  integer                                    :: i
  integer                                    :: count_num_noise
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: pcalc_z2, pcalc_z3
  real(dp)                                   :: error
  real(dp)                                   :: abs_error
  real(dp)                                   :: dydx
  real(dp)                                   :: delta_eta
  real(dp)                                   :: eta_old
  logical                                    :: write_output
  !-----------------------------------------------------------------------------

  eta_iteration = eta_start

  i = 0
  count_num_noise = 0
  eta_old = eta_start
  delta_eta = 0.0_dp

  write_output = .false.

  !-----------------------------------------------------------------------------
  ! iterate density until p_calc = p
  !-----------------------------------------------------------------------------

  abs_error = density_tolerance + 1.0_dp

  do while ( abs_error > density_tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z )

     error = ( pcalc / p ) - 1.0_dp

     if ( write_output .OR. outp > 1 ) write (*,'(i4,4G25.15)') i, error, error / pcalc_z*p, pcalc, eta_iteration

     !--------------------------------------------------------------------------
     ! correction for instable region
     !--------------------------------------------------------------------------
     if ( pcalc_z < 0.0_dp .AND. i < max_iterations ) then

        if ( .NOT.write_output .AND. outp == 1 ) write (*,'(i4,5G25.15)') i, error, error / pcalc_z*p,  &
                                                                          pcalc, eta_iteration, eta_start
        if ( outp == 1 ) write_output = .true.
        if ( i == 1 .AND. eta_old <= 0.25_dp ) delta_eta = 1.9_dp * eta_old
        if ( i == 1 .AND. eta_old >  0.25_dp ) delta_eta = - 0.08_dp
           
        eta_iteration = eta_old - 0.5_dp * delta_eta

        i = i + 1
        call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

        error = ( pcalc / p ) - 1.0_dp

        if ( write_output .OR. outp > 1 ) write (*,'(i4,5G25.15,a)') i, error, error / pcalc_z*p,  &
                                                                     pcalc, eta_iteration, eta_old, ' BBB'

        if ( eta_iteration > 0.55_dp ) then                               ! artificial density root
           call pressure_spinodal ( eta_iteration, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( eta_iteration > 0.55_dp ) then
              if ( error <  0.0_dp ) return                                  ! no liquid solution possible
              if ( error >= 0.0_dp ) eta_iteration = eta_iteration * 0.98_dp ! no solution found so far
           else
              if ( error >  0.0_dp ) eta_iteration = 0.001_dp                ! no liquid solution possible
              if ( error <= 0.0_dp ) eta_iteration = eta_iteration * 1.1_dp  ! no solution found so far
           end if
        else if ( error > 0.0_dp .AND. pcalc_z2 > 0.0_dp ) then           ! no liquid density found
           call pressure_spinodal ( eta_start, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( error >  0.0_dp ) eta_iteration = 0.001_dp                ! no liquid solution possible
           if ( error <= 0.0_dp ) eta_iteration = eta_iteration * 1.1_dp  ! no solution found so far
        else if ( error < 0.0_dp .AND. pcalc_z2 < 0.0_dp ) then           ! no vapor density found
           call pressure_spinodal ( eta_start, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( error <  0.0_dp ) eta_iteration = 0.5_dp                  ! no vapor solution possible
           if ( error >= 0.0_dp ) eta_iteration = eta_iteration * 0.9_dp  ! no solution found so far
        else
           if ( i /= 2 ) then           ! first approach
              eta_iteration = (eta_iteration + eta_start) / 2.0_dp
              if ( abs( eta_iteration - eta_start) < machine_eps ) eta_iteration = eta_iteration + 0.2_dp
           end if
        end if

        CYCLE

     end if


     !--------------------------------------------------------------------------
     ! Newton iteration of density (eta)
     !--------------------------------------------------------------------------
     eta_old = eta_iteration

     dydx = pcalc_z / p
     delta_eta = error / dydx
     if ( abs( delta_eta ) >  0.05_dp ) delta_eta = sign( 0.05_dp, delta_eta )     

     eta_iteration = eta_iteration - delta_eta

     if (eta_iteration > 0.9_dp)  eta_iteration = 0.6_dp
     if (eta_iteration <= 0.0_dp) eta_iteration = 1.E-15_dp

     !--------------------------------------------------------------------------
     ! convergence criteria:
     ! (1) abs_error < density_tolerance, the normal criterium
     ! (2) if step size 'delta_eta' reaches limit of machine-precision, two or
     ! three times in a row
     !--------------------------------------------------------------------------
     abs_error = abs( error )
     if ( eta_iteration < 0.05_dp ) abs_error = abs_error * 1.E4_dp

     if ( abs( delta_eta ) < 100.0_dp * eta_iteration * machine_eps ) then
        count_num_noise = count_num_noise + 1
        if ( ( count_num_noise >= 2 .AND. abs_error < 100.0_dp * density_tolerance ) .OR.  &
               count_num_noise >= 3 )  abs_error = 0.0_dp         ! exit condition
     end if

  end do

  if ( abs_error > density_tolerance ) then
     write (*,*) 'density iteration failed', abs_error, density_tolerance
     if ( outp > 0 ) write (*,*) 'error step size criterion', abs( delta_eta ), 50.0_dp * eta_iteration * machine_eps
     if ( outp > 0 ) read (*,*)
     ! stop
  end if

  eta = eta_iteration
  if ( write_output .OR. outp > 0 ) write (*,*) 'now leaving density_iteration'
  if ( write_output .OR. outp > 0 ) write (*,*) ' '

end subroutine density_iteration


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pressure_spinodal
!
!> \brief pressure_spinodal
!!
!! iterates the density until the derivative of pressure 'pcalc' to
!! density is equal to zero. A Newton-scheme is used for determining
!! the root to the objective function.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure_spinodal ( eta_start, eta, pcalc, pcalc_z )

  use PARAMETERS, only: dp
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z

  !-----------------------------------------------------------------------------
  real(dp), parameter                        :: tolerance = 1.E-5_dp
  integer, parameter                         :: max_iterations = 100

  integer                                    :: i
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: delta_eta
  real(dp)                                   :: error
  real(dp)                                   :: pcalc_z2, pcalc_z3
  !-----------------------------------------------------------------------------

  i = 0
  eta_iteration = eta_start

  !-----------------------------------------------------------------------------
  ! iterate density until d(p)/d(eta) = 0
  !-----------------------------------------------------------------------------

  error = tolerance + 1.0_dp
  do while ( abs( error ) > tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

     error = pcalc_z

     delta_eta = error / pcalc_z2
     if ( delta_eta >  0.05_dp ) delta_eta = 0.05_dp
     if ( delta_eta < -0.02_dp ) delta_eta = -0.02_dp

     eta_iteration   = eta_iteration - delta_eta
     ! write (*,'(a,i3,3G18.10)') 'iter',i, error, eta_iteration

     if ( eta_iteration > 0.7_dp )  eta_iteration = 0.5_dp
     if ( eta_iteration <= 0.0_dp ) eta_iteration = 1.E-16_dp

  end do

  eta = eta_iteration

end subroutine pressure_spinodal


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pressure_inflection_point
!
!> \brief calculates inflection point of pressure isotherm
!!
!! iterates the density until the second derivative of pressure 'pcalc' to
!! density is equal to zero (constant T,x). A Newton-scheme is used for determining
!! the root to the objective function.
!!
!! A recommended intitial value for eta_start is eta_start = 0.2
!! (at low temperatures, it might be advisable to use eta_start = 0.25 )
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure_inflection_point ( eta_start, eta, pcalc, pcalc_z, pcalc_z2, pcalc_z3, converg )

  use PARAMETERS, only: dp
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z
  real(dp), intent(out)                      :: pcalc_z2
  real(dp), intent(out)                      :: pcalc_z3

  !-----------------------------------------------------------------------------
  real(dp), parameter                        :: tolerance = 1.E-2_dp
  integer, parameter                         :: max_iterations = 100

  integer                                    :: i
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: delta_eta
  real(dp)                                   :: error
  logical                                    :: converg
  !-----------------------------------------------------------------------------

  eta_iteration = eta_start

  i = 0
  converg = .false.

  !-----------------------------------------------------------------------------
  ! iterate density until d^2(p)/d(eta)^2 = 0
  !-----------------------------------------------------------------------------

  error = tolerance + 1.0_dp
  do while ( abs(error) > tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

     error = pcalc_z2

     delta_eta = error / pcalc_z3
     if ( delta_eta >  0.08_dp ) delta_eta = 0.08_dp
     if ( delta_eta < -0.05_dp ) delta_eta = -0.05_dp

     eta_iteration   = eta_iteration - delta_eta
     ! write (*,'(a,i3,3G18.10)') 'iter inflex',i, error, eta_iteration

     if ( eta_iteration > 0.7_dp )  eta_iteration = 0.5_dp
     if ( eta_iteration <= 0.0_dp ) eta_iteration = 0.01_dp

  end do

  eta = eta_iteration
  
  if ( abs( error ) < tolerance ) converg = .true.

  if ( pcalc < 1.E-20_dp ) pcalc = 1.E-10_dp

end subroutine pressure_inflection_point

