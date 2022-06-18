!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module stability_G

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  use stability, only: get_eta
  implicit none

  integer, parameter                              :: outp = 0
  real(dp), allocatable, dimension(:)             :: rhoi_coexisting

  private
  public :: stability_analysis_G

CONTAINS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine stability_analysis
!> \brief stability_analysis_G
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

subroutine stability_analysis_G ( t, p, x_feed, eta_start, stable, rhoi_feed, ln_zi_phi, rhoi_out, rhoi_coex )

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

  integer                                         :: i
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
     allocate( rhoi_coexisting( ncomp ) )
     rhoi_coexisting( 1:ncomp ) = rhoi_coex( 1:ncomp )
  end if

  !-----------------------------------------------------------------------------
  ! is the composition of one or more species = 0. Then, only here, set to 1.E-200
  !-----------------------------------------------------------------------------
  !do i = 1, ncomp
  !   if ( x_feed(i) <= 1.E-100_dp ) x_feed(i) = 1.E-100_dp
  !end do

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
     do i = 1, ncomp
        if ( x_feed(i) > 1.E-100_dp ) then
           ln_zi_phi(i) = log( x_feed(i) ) + mu_res(i)
        else
           ln_zi_phi(i) = log( 1.E-100_dp ) + mu_res(i)
        end if
     end do

     call define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_trial )
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase eta ', get_eta( rhoi_trial(:) )
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase xi  ',  &
          rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )

     !--------------------------------------------------------------------------
     ! minimizing the objective fct. Phase split for values of tpd < 0.0
     !--------------------------------------------------------------------------
     call minimize_tpd_G ( t, mu_feed, p_kT, rhoi_feed, rhoi_trial, tpd )

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

end subroutine stability_analysis_G


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
! SUBROUTINE minimize_tpd_G
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine minimize_tpd_G ( t, mu_feed, p_kT, rhoi_feed, rhoi, tpd )

  use PARAMETERS, only: KBOL30
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
  integer                                         :: i
  integer                                         :: i_iter
  integer, parameter                              :: N_max_iter = 100
  real(dp), parameter                             :: tolerance = 1.E-6_dp
  real(dp)                                        :: scaled_tolerance
  real(dp), dimension(ncomp)                      :: x_old
  real(dp), dimension(ncomp)                      :: Y_old
  real(dp)                                        :: tpd_old
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu
  real(dp)                                        :: error
  logical                                         :: trivial
  
  real(dp)                                        :: p
  real(dp)                                        :: eta_in
  real(dp), dimension(ncomp)                      :: x
  real(dp), dimension(ncomp)                      :: di, lnphi_feed
  real(dp), dimension(ncomp)                      :: Y
  real(dp), dimension(ncomp)                      :: x_feed

  logical                                         :: newton
  !-----------------------------------------------------------------------------

  scaled_tolerance = tolerance
  
  !-----------------------------------------------------------------------------
  ! fugacity coeff. of feed-phase
  !-----------------------------------------------------------------------------
  do i = 1, ncomp

     if ( rhoi_feed(i) > 1.E-200_dp ) then

        lnphi_feed(i) = mu_feed(i) - log( rhoi_feed(i) ) - log( p_kT/sum( rhoi_feed(:) ) )
        x_feed(i) = rhoi_feed(i) / sum( rhoi_feed(:) )

     else

        lnphi_feed(i) = mu_feed(i) - log( 1.E-200_dp ) - log( p_kT/sum( rhoi_feed(:) ) )
        x_feed(i) = 1.E-200_dp

     end if

  end do

  di(:) = log( x_feed(:) ) + lnphi_feed(:)

  ! JG: this is not a good starting value: it is not used in the successive-substitution routine
  Y(1:ncomp) = x_feed(1:ncomp)
  Y_old(1:ncomp) = x_feed(1:ncomp)

  x(1:ncomp) = rhoi(1:ncomp) / sum( rhoi(1:ncomp) )
  x_old(1:ncomp) = x(1:ncomp)

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

        p = p_kT * ( KBOL30 * t )

        eta_in = get_eta( rhoi )

        call chemical_potential_tp ( t, p, x, eta_in, rhoi, mu_res, mu = mu )

        !-----------------------------------------------------------------------
        ! direct substitution step
        !-----------------------------------------------------------------------
        Y(1:ncomp) = exp( di(1:ncomp) - mu_res(1:ncomp) )

        tpd = 1.0_dp + sum( Y(:) * ( log( Y(:) ) + mu_res(:) - di(:) - 1.0_dp ) )

        x(:) = Y(:) / sum( Y(:) )
        g_function(1:ncomp) = abs( x(1:ncomp) - x_old(1:ncomp) )

        if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error S.Sub.', i_iter, tpd,  &
                sum( abs( g_function(1:ncomp) ) )

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( g_function(1:ncomp) ) )

        x_old( 1:ncomp ) = x(1:ncomp)
        Y_old(:) = Y(:)                   ! this line is needed only for Newton-schme
        tpd_old = tpd

        !-----------------------------------------------------------------------
        ! exit condition towards Newton-schme
        !-----------------------------------------------------------------------
        if ( i_iter >= 15 .AND. error > scaled_tolerance .OR.  &
           ( tpd > tpd_old + 1.E-5_dp .AND. i_iter > 2 ) ) then
                 newton = .true.
                 if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton'
        end if

     else

        call stability_G_Newton_step ( i_iter, t, p_kT, rhoi, di, Y_old, x_old, tpd_old, Y, x, tpd, error )

        !-----------------------------------------------------------------------
        ! store previous values
        !-----------------------------------------------------------------------
        Y_old(1:ncomp) = Y(1:ncomp)
        x_old(1:ncomp) = x(1:ncomp)
        tpd_old = tpd

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
     if ( tpd < - 0.01_dp ) scaled_tolerance = tolerance * 10.0_dp
     if ( tpd < - 0.1_dp ) scaled_tolerance = tolerance * 100.0_dp
     if ( tpd < - 0.1_dp .AND. i_iter > 5) scaled_tolerance = tolerance * 1000.0_dp

  end do     

end subroutine minimize_tpd_G


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE stability_G_Newton_step
!
! perform Newton step in stability analysis. A preceeding direct substitution
! step is required to have initial values for Y-vector.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_G_Newton_step ( i_iter, t, p_kT, rhoi, di, Y_old, x_old, tpd_old, Y, x, tpd, error )

  use PARAMETERS, only: KBOL30
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                         :: i_iter
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p_kT
  real(dp), dimension(ncomp), intent(in out)      :: rhoi
  real(dp), dimension(ncomp), intent(in)          :: di
  real(dp), dimension(ncomp), intent(in)          :: Y_old
  real(dp), dimension(ncomp), intent(in)          :: x_old
  real(dp), intent(in)                            :: tpd_old
  real(dp), dimension(ncomp), intent(out)         :: Y
  real(dp), dimension(ncomp), intent(in out)      :: x
  real(dp), intent(out)                           :: tpd
  real(dp), intent(out)                           :: error

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  
  real(dp)                                        :: p
  real(dp)                                        :: eta_in

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
  !-----------------------------------------------------------------------------

  p = p_kT * ( KBOL30 * t )

  eta_in = get_eta( rhoi )

  call chemical_potential_tp ( t, p, x, eta_in, rhoi, mu_res, lnphi_nj=Hesse )

  g_function(1:ncomp) = ( log( Y(1:ncomp) ) + mu_res(1:ncomp) - di(1:ncomp) ) * sqrt( Y(1:ncomp) )

  Hesse_ig(:,:) = 0.0_dp
  ! Hesse(:,:) = 0.0_dp
  do i = 1, ncomp
     Hesse(i,1:ncomp) = Hesse(i,1:ncomp) * sqrt( Y(i)*Y(1:ncomp) )
     Hesse(i,i) = Hesse(i,i) + log( Y(i) ) + mu_res(i) - di(i)
     Hesse_ig(i,i) = 1.0_dp
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
        if ( abs( ( in_out_vector(i) / 2.0_dp )**2 / Y_old(i) ) > 5.00_dp ) then
           adjust_hessian = .true.
        end if
     end do
     if ( adjust_hessian ) then
        eta_H = eta_H + 0.5_dp
        if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H II', eta_H
        CYCLE      ! cycle, because of Hessian-criterion (2): too large step-size
     end if

     !--------------------------------------------------------------------
     ! new vector Y
     !--------------------------------------------------------------------

     Y(1:ncomp) = ( sqrt(Y_old(1:ncomp)) - in_out_vector(1:ncomp) / 2.0_dp )**2

     tpd = 1.0_dp + sum( Y(:) * ( log( Y(:) ) + mu_res(:) - di(:) - 1.0_dp ) )

     x(:) = Y(:) / sum( Y(:) )

     if ( tpd > tpd_old .AND. eta_H < 30.0_dp .AND. i_iter >= 3 ) then
        eta_H = eta_H + 0.5_dp
        adjust_hessian = .true.    ! cycle, because of Hessian-criterion (3): tpd does not descent
        if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H III', tpd, tpd_old
     end if

     if ( outp >= 1 ) write (*,'(a,i4,4G25.15)') ' N, tpd, error Newton', i_iter, tpd,  &
          sum( abs( g_function(1:ncomp) ) ), sum( abs( x(1:ncomp) - x_old(1:ncomp) ) ), eta_H

  end do

  !-----------------------------------------------------------------------
  ! determine error
  !-----------------------------------------------------------------------

  error = sum( abs( g_function(1:ncomp) ) )


end subroutine stability_G_Newton_step


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


end module stability_G
