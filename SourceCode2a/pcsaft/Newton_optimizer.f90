!> \file Newton_optimizer.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
! subroutine Newton_minimizer
!
! minimizes a function f(x) with respect to vector x(1:n) as degrees of freedom.
! The objective function f(x) is provided by "objective_fct". The gradient g of
! f with respect to x and the second derivative H are provided by
! "f_grad_hessian". Voilations of constraints on x are detected in subroutine
! "violation_of_parameter_bound".
!
! IMPORTANT: 
! A prerequisite for this algorithm is a Hessian with positive valued diagonal
! elements (see matrix "H_diagonal"). The method of Murray is applied, ensuring
! a descent step by adding (eta_H * H_diagonal) to Hessian with eta_H chosen as
! zero for default. We set eta_H > 0 if
!    (1) H is not positive definite
!    (2) step size is too large
!    (3) objective function (tpd) does not descent
! Increasing eta_H leads to increasingly small steps in a direction that increas-
! ingly represents the steepest descent direction.
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module optimizer_JG

  use PARAMETERS, only: dp
  implicit none

  integer, parameter                    :: outp = 0

  private
  public :: Newton_minimizer

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine Newton_minimizer
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine Newton_minimizer ( objective_fct, f_grad_hessian, violation_of_parameter_bound,  &
       n, grad_tolerance, para_tolerance, N_max_iter, x, f, g, converged )

  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: n                 ! number of degrees of freedom of minimization problem, x(1:n)
  real(dp), intent(in)                            :: grad_tolerance    ! convergence tolerance for gradient of obc. fct.
  real(dp), intent(in)                            :: para_tolerance    ! convergence tolerance for parameter step
  integer, intent(in)                             :: N_max_iter        ! maximum number of iterations
  real(dp), dimension(n), intent(in out)          :: x                 ! value of degrees of freedom
  real(dp), intent(out)                           :: f                 ! objective function, f(x) with x(n)
  real(dp), dimension(n), intent(out)             :: g                 ! gradient of objective function
  logical, intent(out)                            :: converged         ! is .true. if one of the convergence tolerances were achieved

  !-----------------------------------------------------------------------------
  interface
     subroutine objective_fct ( n, x, f )                              ! objective function, f(x), with x(n)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(out)                      :: f
     end subroutine objective_fct

     subroutine f_grad_hessian ( n, x, g, H, H_diagonal )             ! gradient g(n) & Hessian H(n,n) of obj. fct. f(x), with x(n)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       !real(dp), intent(out)                      :: f
       real(dp), dimension(n), intent(out)        :: g
       real(dp), dimension(n,n), intent(out)      :: H
       real(dp), dimension(n,n), intent(out)      :: H_diagonal
     end subroutine f_grad_hessian

     subroutine violation_of_parameter_bound ( n, delta_x, x_old, bound_violation )  ! upper or lower bound of x(n) violated or max step size exceeded ?
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: delta_x
       real(dp), dimension(n), intent(in)         :: x_old
       logical, intent(out)                       :: bound_violation
     end subroutine violation_of_parameter_bound
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: count
  integer                                         :: it_num, rot_num
  integer, parameter                              :: it_max = 200
  real(dp)                                        :: error_in_x, error_in_g
  real(dp)                                        :: f_old
  real(dp), dimension(n)                          :: x_old
  real(dp), dimension(n)                          :: delta_x
  real(dp), dimension(n,n)                        :: hessian, H, H_diagonal
  real(dp), dimension(n,n)                        :: eigenvec
  real(dp), dimension(n)                          :: eigenval
  real(dp)                                        :: eta_H
  real(dp)                                        :: delta_eta
  real(dp)                                        :: trace_diagonal
  logical                                         :: adjust_hessian
  logical                                         :: bound_violation
  real(dp)                                        :: DET
  character(LEN=2)                                :: char_len  
  !-----------------------------------------------------------------------------

  count = 0
  converged = .false.

  x_old(1:n) = x(1:n)
  f_old = 1.E50_dp

  do while ( .NOT.converged .AND. count <= N_max_iter )

     count = count + 1

     call f_grad_hessian ( n, x, g, H, H_diagonal )

     !--------------------------------------------------------------------------
     ! use method of Murray, by adding a unity matrix to Hessian, if:
     ! (1) H is not positive definite
     ! (2) step size is too large
     ! (3) objective function (tpd) does not descent
     !--------------------------------------------------------------------------
     adjust_hessian = .true.
     eta_H = 0.0_dp

     do while ( adjust_hessian )

        adjust_hessian = .false.

        Hessian(:,:) = H(:,:) + eta_H * H_diagonal(:,:)

        call jacobi_eigenvalue ( n, Hessian, it_max, eigenvec, eigenval, it_num, rot_num )
        if ( outp >= 3 )  write (*,*) 'lowest eigenvalue',eigenval(1)

        if ( eigenval(1) < 0.0001_dp .AND. eta_H < 20.0_dp ) then
           if ( outp >= 2 ) write (*,*) 'Newton_minimizer: increasing eta_H I, eigenval=',eigenval(1)
           trace_diagonal = sum( H_diagonal(1:n,1:n) ) / real( n, KIND=dp )
           delta_eta = max( 0.001_dp, 1.2_dp*abs( 0.0001_dp-eigenval(1) ) / trace_diagonal )
           if ( eta_H > 0.01 ) delta_eta = 0.02_dp
           if ( eta_H > 0.1 ) delta_eta = 0.05_dp
           if ( eta_H > 1.0 ) delta_eta = 0.5_dp
           eta_H = eta_H + delta_eta
           adjust_hessian = .true.
           CYCLE      ! cycle, because of Hessian-criterion (1): H not positive definite
        end if

        !-----------------------------------------------------------------------
        ! solving  AX = B  with Gauss-Jordan method.
        ! Matrix A='Hessian' and vector B='delta_x' are destroyed and
        ! redefined to the inverse of A and solution vector X, respectively.
        !-----------------------------------------------------------------------

        delta_x(:) = - g(:)
        call MATINV ( n, 1, Hessian, delta_x, DET )

        write (char_len,'(I2)') n
        if ( outp >= 2 ) write (*,'(a,'//char_len//'G18.8)') ' Newton_minimizer: Newton-step for delta_x(:)', delta_x(:)

        call violation_of_parameter_bound ( n, delta_x, x_old, bound_violation )

        if ( bound_violation ) then
           if ( eta_H < 30.0_dp ) then
              delta_eta = 0.01_dp
              if ( eta_H > 0.1 ) delta_eta = 0.05_dp
              if ( eta_H > 1.0 ) delta_eta = 0.5_dp
              eta_H = eta_H + delta_eta
              adjust_hessian = .true.
              if ( outp >= 2 ) write (*,*) 'Newton_minimizer: increasing eta_H II', eta_H
              CYCLE      ! cycle, because of Hessian-criterion (2): too large step-size
           else
              write (*,*) 'do something here: remove bound_violation and proceed for one or two steps'
              stop
           end if
        end if

        !-----------------------------------------------------------------------
        ! new x vector
        !-----------------------------------------------------------------------
        x(:) = x_old(:) + delta_x(:)

        !-----------------------------------------------------------------------
        ! calculate value of function f
        !-----------------------------------------------------------------------
        call objective_fct ( n, x, f )

        if ( f > f_old + 1.E-5_dp .AND. eta_H < 30.0_dp ) then
           delta_eta = 0.025_dp
           if ( eta_H >= 0.1 ) delta_eta = 0.2_dp
           if ( eta_H >= 1.0 ) delta_eta = 0.5_dp
           eta_H = eta_H + delta_eta
           adjust_hessian = .true.      ! cycle, because of Hessian-criterion (3): tpd does not descent
           if ( outp >= 2 ) write (*,*) 'Newton_minimizer: increasing eta_H III'
        end if

     end do

     !--------------------------------------------------------------------------
     ! if f is successfully decreased
     !--------------------------------------------------------------------------
     if ( f < f_old + 1.E-5_dp ) then
        x_old(1:n) = x(1:n)
        f_old = f
     else
        if ( outp >= 1 ) write (*,*) 'Newton_minimizer: no decrease of objective fct possible'
        x(1:n) = x_old(1:n)
        f = f_old
     end if

     !--------------------------------------------------------------------------
     ! determine error
     !--------------------------------------------------------------------------

     error_in_g = sum( abs( g(1:n) ) )

     error_in_x = sum( abs( delta_x(1:n) ) )

     if ( error_in_x < para_tolerance .OR. error_in_g < grad_tolerance ) then
        converged = .true.
     end if

     if ( outp >= 1 ) write (*,'(a,i3,3G25.15)') ' finished Newton step: N, f, errors ',count, f, error_in_g, error_in_x
     if ( outp >= 3 ) read(*,*)


  end do

end subroutine Newton_minimizer

end module optimizer_JG
