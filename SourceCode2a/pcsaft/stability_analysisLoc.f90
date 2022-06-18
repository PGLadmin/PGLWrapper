!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module stability

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                              :: outp = 1
  real(dp), allocatable, dimension(:)             :: rhoi_coexisting

  private
  public :: stability_analysis, get_eta

CONTAINS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine stability_analysis
!> \brief stability_analysis
!!
!!
!! input :   t, p, x_feed(:)
!!           eta_trial    :  starting value for density of trial phase
!!           rhoi_coex(:) :  ( optional ) species densities of a phase that is
!!                           already known to be in phase equilibrium. This
!!                           subroutine ignores this phase as a trial phase
!!
!! output:   stable=.true.:  for stable phase or otherwise for a phase split
!!           rhoi_out(:)  :  species densities of trial phase (minimum
!!                           of the constrained tangent plane distance).
!!           rhoi_feed    :  density vector of feed phase
!!           ln_zi_phi    :  residual chemical potential of feed phase
!!                           mu_res(:)=ln( x_feed(:) * phi_feed(:) )
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_analysis ( t, p, x_feed, eta_start, stable, rhoi_feed, ln_zi_phi, rhoi_out, rhoi_coex )

  use PARAMETERS, only: KBOL30 !, machine_eps
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), intent(in)                            :: eta_start
  logical, intent(out)                            :: stable
  real(dp), dimension(ncomp), intent(out)         :: rhoi_feed
  real(dp), dimension(ncomp), intent(out)         :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_out
  real(dp), dimension(ncomp), intent(in), optional :: rhoi_coex
  !-----------------------------------------------------------------------------

  integer                                         :: trial_steps

  integer                                         :: i_trial

  real(dp)                                        :: tpd
  real(dp), dimension(ncomp)                      :: rhoi_trial
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu_feed
  real(dp)                                        :: delta_rhoi
  real(dp)                                        :: p_kT
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! is a coexisting phase to (t,p,x_feed) already known?  Then allocate array.
  !-----------------------------------------------------------------------------
  if ( present( rhoi_coex ) ) then
     allocate( rhoi_coexisting(ncomp) )
     rhoi_coexisting( 1:ncomp ) = rhoi_coex(1:ncomp)
  end if

  !-----------------------------------------------------------------------------
  ! starting values: one vapor & ncomp liquids, each with one species as dominating
  !-----------------------------------------------------------------------------
  trial_steps = ncomp + 1

  stable = .true.

  call chemical_potential_tp ( t, p, x_feed, eta_start, rhoi_feed, mu_res, mu=mu_feed )

  p_kT = p / ( KBOL30 * t )

  do i_trial = 1, trial_steps

     !--------------------------------------------------------------------------
     ! setting trial-phase mole-fractions (exp(ln_zi_phi): id.gas estimate for vapor)
     !--------------------------------------------------------------------------
     ln_zi_phi( 1:ncomp ) = log( x_feed( 1:ncomp ) ) + mu_res( 1:ncomp )

     call define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_trial )
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase eta ', get_eta( rhoi_trial(:) )
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase xi  ',  &
          rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )

     !--------------------------------------------------------------------------
     ! minimizing the objective fct. Phase split for values of tpd < 0.0
     !--------------------------------------------------------------------------
     call minimize_tpd ( t, mu_feed, p_kT, rhoi_feed, rhoi_trial, tpd )

     !--------------------------------------------------------------------------
     ! check for trivial solution ( where delta_rhoi = 0, approximately )
     !--------------------------------------------------------------------------
     delta_rhoi = maxval( abs( 1.0_dp - rhoi_trial(1:ncomp) / rhoi_feed(1:ncomp) ) )

     if ( tpd < -1.E-8_dp .AND. delta_rhoi > 0.001_dp ) then
        stable = .false.
        rhoi_out( 1:ncomp ) = rhoi_trial( 1:ncomp )
        if ( outp >= 1 ) then
           write (*,'(a,i4,G20.12)') ' unstable: N, tpd  ', i_trial, tpd
           write (*,'(a, 10G20.12)') ' unstable: eta, xi ', get_eta( rhoi_trial(:) ),  &
                 rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )
        end if
        exit
     end if

  end do

  !-----------------------------------------------------------------------------
  ! output to terminal
  !-----------------------------------------------------------------------------
  if ( outp >= 1 .AND. stable ) then
     write (*,'(a,i4,2G20.12)') ' stability analysis: stable phase'
  end if
  if ( outp >= 1 ) write (*,*) ' '

  if ( present( rhoi_coex ) ) deallocate( rhoi_coexisting )

end subroutine stability_analysis


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine define_initial_rhoi
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_guess )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i_trial
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), dimension(ncomp), intent(in)          :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_guess
  !-----------------------------------------------------------------------------

  integer                                         :: i
  integer                                         :: dominant_component
  real(dp)                                        :: x_dominant
  real(dp)                                        :: eta
  real(dp)                                        :: rho
  real(dp)                                        :: factor
  real(dp), dimension(ncomp)                      :: x
  !-----------------------------------------------------------------------------

  if ( i_trial == 1 ) then
     !------- vapor trial phase ------------------------------------------------
     x( 1:ncomp ) = exp( ln_zi_phi( 1:ncomp ) )
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )
     eta = 0.0001_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  else !if ( i_trial == 2 ) then
     !------- liquid trial phase ------------------------------------------------
     dominant_component = i_trial - 1
     x_dominant = 0.99_dp
     x( dominant_component ) = x_dominant
     factor = ( 1.0_dp - x_dominant ) / ( sum( x_feed(1:ncomp) ) - x_feed(dominant_component) )
     do i = 1, ncomp
        if ( i /= dominant_component) then
           x(i) = x_feed(i) * factor
        end if
     end do
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )  ! this line is actually not needed
     eta = 0.4_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  end if


end subroutine define_initial_rhoi


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE minimize_tpd
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine minimize_tpd ( t, mu_feed, p_kT, rhoi_feed, rhoi, tpd )

  use properties, only: chemical_potential_tp, chemical_potential_trho,  &
       Helmholtz_density
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), dimension(ncomp), intent(in)          :: mu_feed
  real(dp), intent(in)                            :: p_kT
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), dimension(ncomp), intent(inout)       :: rhoi
  real(dp), intent(out)                           :: tpd

  !-----------------------------------------------------------------------------
  integer                                         :: i_iter, i
  integer                                         :: critical_comp(1)
  integer, parameter                              :: N_max_iter = 100
  real(dp), parameter                             :: tolerance = 1.E-3_dp
  real(dp)                                        :: scaled_tolerance
  real(dp), dimension(ncomp)                      :: rhoi_old
  real(dp)                                        :: tpd_old
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu
  real(dp)                                        :: error
  real(dp)                                        :: f_dens
  logical                                         :: trivial
  real(dp)                                        :: eta_new

  logical                                         :: repeat
  real(dp)                                        :: damp

  real(dp), dimension(ncomp,ncomp)                :: Hessian
  real(dp), dimension(ncomp,ncomp)                :: Hesse
  real(dp), dimension(ncomp,ncomp)                :: Hesse_ig
  real(dp)                                        :: eta_H
  logical                                         :: adjust_hessian
  integer, parameter                              :: it_max = 200
  integer                                         :: it_num
  integer                                         :: rot_num
  real(dp), dimension(ncomp,ncomp)                :: eigenvec
  real(dp), dimension(ncomp)                      :: eigenval
  real(dp), dimension(ncomp)                      :: in_out_vector
  real(dp)                                        :: DET

  logical                                         :: newton
  !-----------------------------------------------------------------------------

  rhoi_old(1:ncomp) = rhoi(1:ncomp)

  scaled_tolerance = tolerance
  damp = 1.0_dp     ! damping factor of direct substituation scheme
  
  i_iter = 0
  error = tolerance + 1.0_dp
  newton = .false.
  tpd_old = 1.E50_dp

  do while ( error > scaled_tolerance .AND. i_iter < N_max_iter )

     i_iter = i_iter + 1

     if ( .NOT.newton ) then

        !-----------------------------------------------------------------------
        ! case: direct substitution
        !-----------------------------------------------------------------------

        ! call direct_substitution_step ( t, mu_feed, p_kT, rhoi, mu_res, mu, tpd, damp, rhoi_old, tpd_old )
        ! if ( i_iter > 1 ) then
           call chemical_potential_trho ( t, rhoi, mu_res, mu )
        ! else
        !   x(:) = rhoi(:) / sum(rhoi(1:ncomp))
        !   p = p_kT * ( KBOL30 * t )
        !   call chemical_potential_tp ( t, p_kT*( KBOL30 * t ), x, 0.0001_dp, rhoi, mu_res, mu )
        ! end if
        g_function(1:ncomp) = mu(1:ncomp) - mu_feed(1:ncomp)

        !-----------------------------------------------------------------------
        ! constrain step size, if needed
        !-----------------------------------------------------------------------
        critical_comp = MINLOC( exp( mu_feed(1:ncomp) - mu(1:ncomp) ) )
        damp = min( damp, exp( mu_feed(critical_comp(1) ) - mu(critical_comp(1) ) ) / 1.E-2_dp )
        if ( outp >= 2 ) write (*,*) 'minloc', exp( mu_feed(critical_comp(1) ) -mu(critical_comp(1) )), damp        
        critical_comp = MAXLOC( exp( mu_feed(1:ncomp) - mu(1:ncomp) ) )
        damp = min( damp, 1.E2_dp / exp( mu_feed(critical_comp(1) ) - mu(critical_comp(1) ) ) )
        if ( outp >= 2 ) write (*,*) 'maxloc', exp( mu_feed(critical_comp(1) ) -mu(critical_comp(1) )), damp

        repeat = .true.
        do while ( repeat )

           repeat = .false.

           !--------------------------------------------------------------------
           ! direct substitution step
           !--------------------------------------------------------------------
           rhoi(1:ncomp) = damp * exp( mu_feed(1:ncomp) - mu(1:ncomp) ) * rhoi_old(1:ncomp)  &
                + ( 1.0_dp - damp ) * rhoi_old(1:ncomp)

           !--------------------------------------------------------------------
           ! for too high densities, reduce rhoi
           !--------------------------------------------------------------------
           eta_new = get_eta ( rhoi )
           if ( eta_new > 0.6_dp ) rhoi(:) = rhoi(:) * 0.6_dp / eta_new
           ! do i = 1, ncomp
           !    write (*,*) 'rhoi',rhoi_old(i), rhoi(i)
           ! end do

           !--------------------------------------------------------------------
           ! calculate tangent plane distance tpd
           !--------------------------------------------------------------------
           call Helmholtz_density ( t, rhoi, f_dens )

           tpd = f_dens + p_kT - sum( mu_feed(1:ncomp) * rhoi(1:ncomp) )
           tpd = tpd * 1000.0_dp
           !write (70,'(i3,11G25.15)') i_iter, tpd, sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) ),  &
           !     eta_new, rhoi(1:ncomp)

           !--------------------------------------------------------------------
           ! if tpd is not successfully decreased, use stronger damping
           !--------------------------------------------------------------------
           if ( tpd > tpd_old .AND. damp > 0.04_dp ) then
              repeat = .true.
              damp = damp * 0.4_dp
              if ( outp >= 2 )  write (*,*) 'inner loop back', damp
           end if

        end do

        !-----------------------------------------------------------------------
        ! if tpd is successfully decreased
        !-----------------------------------------------------------------------
        if ( tpd < tpd_old + 1.E-5_dp ) then
           if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error DS    ', i_iter, tpd,  &
                sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) ), damp
           rhoi_old(1:ncomp) = rhoi(1:ncomp)
           tpd_old = tpd
           if (damp < 1.0_dp ) damp = min( 1.0_dp, damp * 10.0_dp )
           if ( i_iter >= 15 .OR. ( i_iter >= 10 .AND. eta_new > 0.25_dp ) ) then
              newton = .true.
              if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton'
           end if
        else
           if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton', i_iter, tpd,  &
                sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) )
           rhoi(1:ncomp) = rhoi_old(1:ncomp)
           tpd = tpd_old
           newton = .true.
        end if

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) )

     else

        !-----------------------------------------------------------------------
        ! case Newton step
        !-----------------------------------------------------------------------

        !if ( i_iter == 1 ) then
        !   x(:) = rhoi(:) / sum(rhoi(1:ncomp))
        !   p = p_kT * ( KBOL30 * t )
        !   call chemical_potential_tp ( t, p_kT*( KBOL30 * t ), x, 0.0001_dp, rhoi, mu_res, mu )
        !end if

        call chemical_potential_trho ( t, rhoi, mu_res, mu, Hesse, Hesse_ig )
        g_function(1:ncomp) = ( mu(1:ncomp) - mu_feed(1:ncomp) ) * sqrt( rhoi(1:ncomp) )

        do i = 1, ncomp
           Hesse(i,1:ncomp) = Hesse(i,1:ncomp) * sqrt( rhoi(i)*rhoi(1:ncomp) )
           Hesse_ig(i,1:ncomp) = Hesse_ig(i,1:ncomp) * sqrt( rhoi(i)*rhoi(1:ncomp) ) ! unity matrix
        end do

        !-----------------------------------------------------------------------
        ! use method of Murray, by adding a unity matrix to Hessian, if:
        ! (1) H is not positive definite
        ! (2) step size is too large
        ! (3) objective function (tpd) does not descent
        !-----------------------------------------------------------------------
        adjust_hessian = .true.
        eta_H = 1.0_dp

        do while ( adjust_hessian )

           adjust_hessian = .false.

           Hessian(:,:) = Hesse(:,:) + eta_H * Hesse_ig(:,:)

           call jacobi_eigenvalue ( ncomp, Hessian, it_max, eigenvec, eigenval, it_num, rot_num )
           if ( outp >= 3 ) write (*,*) 'smallest eigenvalue',eigenval(1)

           if ( eigenval(1) < 0.001_dp .AND. eta_H < 20.0_dp ) then
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H I, eigenval=',eigenval(1)
              eta_H = eta_H + 0.5_dp
              adjust_hessian = .true.
              CYCLE      ! cycle, because of Hessian-criterion (1): H not positive definite
           end if

           !--------------------------------------------------------------------
           ! solving  AX = B  with Gauss-Jordan method.
           ! Matrix A='Hessian' and vector B='in_out_vector' are destroyed and
           ! redefined to the inverse of A and solution vector X, respectively.
           !--------------------------------------------------------------------

           in_out_vector = g_function
           call MATINV ( ncomp, 1, Hessian, in_out_vector, DET )
           do i = 1,ncomp
              if ( abs( ( in_out_vector(i) / 2.0_dp )**2 / rhoi(i) ) > 0.99_dp ) then
                 adjust_hessian = .true.
              end if
           end do
           if ( adjust_hessian ) then
              eta_H = eta_H + 0.5_dp
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H II', eta_H
              CYCLE      ! cycle, because of Hessian-criterion (2): too large step-size
           end if

           !--------------------------------------------------------------------
           ! new density vector rhoi
           !--------------------------------------------------------------------
           rhoi(1:ncomp) = ( sqrt(rhoi_old(1:ncomp)) - in_out_vector(1:ncomp) / 2.0_dp )**2
           ! do i = 1, ncomp
           !    write (*,*) 'rhoi',rhoi_old(i), rhoi(i)
           ! end do

           !--------------------------------------------------------------------
           ! for too high densities, reduce rhoi
           !--------------------------------------------------------------------
           eta_new = get_eta( rhoi )
           if ( outp >= 3 )  write (*,*) 'eta_new',eta_new
           if ( eta_new > 0.6_dp ) rhoi(:) = rhoi(:) * 0.6_dp / eta_new

           !--------------------------------------------------------------------
           ! calculate tangent plane distance tpd
           !--------------------------------------------------------------------
           call Helmholtz_density ( t, rhoi, f_dens )

           tpd = f_dens + p_kT - sum( mu_feed(1:ncomp) * rhoi(1:ncomp) )
           tpd = tpd * 1000.0_dp

           if ( tpd > tpd_old .AND. eta_H < 30.0_dp ) then
              eta_H = eta_H + 0.5_dp
              adjust_hessian = .true.    ! cycle, because of Hessian-criterion (3): tpd does not descent
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H III'
           end if

        end do

        if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error Newton', i_iter, tpd, sum( abs( g_function(1:ncomp) ) )

        !-----------------------------------------------------------------------
        ! if tpd is successfully decreased
        !-----------------------------------------------------------------------
        if ( tpd < tpd_old + 1.E-5_dp ) then
           rhoi_old(1:ncomp) = rhoi(1:ncomp)
           tpd_old = tpd
        else
           if ( outp >= 1 ) write (*,*) 'minimize_tpd: no further decrease of objective fct possible', tpd
           rhoi(1:ncomp) = rhoi_old(1:ncomp)
           tpd = tpd_old
        end if

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( g_function(1:ncomp) ) )

     end if

     !--------------------------------------------------------------------------
     ! check if approaching trivial solution or identified coexisting phase, exit
     !--------------------------------------------------------------------------
     call check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )
     if ( trivial ) then
        tpd = 100.0_dp
        exit
     end if

     !--------------------------------------------------------------------------
     ! relax on convergence criterion for clearly unstable trial phase
     !--------------------------------------------------------------------------
     if ( tpd < - 0.1_dp ) scaled_tolerance = tolerance * 10.0_dp
     if ( tpd < - 1.0_dp ) scaled_tolerance = tolerance * 100.0_dp

  end do     

end subroutine minimize_tpd


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )

  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), intent(in)                            :: tpd
  logical, intent(out)                            :: trivial
  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: ratio
  real(dp)                                        :: sum_difference
  !-----------------------------------------------------------------------------

  trivial = .false.

  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards feed (rhoi_feed)
  !-----------------------------------------------------------------------------
  sum_difference = 0.0_dp
  do i = 1, ncomp
     if ( rhoi_feed(i) > 1.E-8_dp ) then
        ratio = rhoi(i) / rhoi_feed(i) - 1.0_dp
     else
        ratio = rhoi(i) - rhoi_feed(i)
     end if
     sum_difference = sum_difference + abs( ratio )
  end do
  sum_difference = sum_difference / real( ncomp, KIND=dp )

  if ( sum_difference < 0.01_dp .AND. tpd < 1.E-3_dp ) then
     trivial = .true.
     if ( outp >= 1 ) write (*,*) 'trivial solution, exit'
  end if

  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards coexising phase (rhoi_coexisting)
  !-----------------------------------------------------------------------------
  if ( allocated( rhoi_coexisting ) .AND. .NOT.trivial ) then

     sum_difference = 0.0_dp
     do i = 1, ncomp
        if ( rhoi_feed(i) > 1.E-8_dp ) then
           ratio = rhoi(i) / rhoi_coexisting(i) - 1.0_dp
        else
           ratio = rhoi(i) - rhoi_coexisting(i)
        end if
        sum_difference = sum_difference + abs( ratio )
     end do
     sum_difference = sum_difference / real( ncomp, KIND=dp )

     if ( sum_difference < 0.03_dp .AND. tpd < 1.E-2_dp ) then
        trivial = .true.
        if ( outp >= 1 ) write (*,*) 'found coexisting phase, exit'
     end if

  end if
     

end subroutine check_trivial_solution



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

function get_eta( rhoi )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi
  real(dp)                                        :: get_eta
  !-----------------------------------------------------------------------------

  get_eta = PI_6 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

end function get_eta


end module stability
