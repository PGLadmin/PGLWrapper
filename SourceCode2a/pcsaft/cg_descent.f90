!> \file cg_descent.f90
!! \brief conjugate gradient method
!!
!!           A conjugate gradient method with guaranteed descent
!!                  Version 1.1  (December 10, 2004)
!!                  Version 1.2  (June 4, 2005)
!!                  Version 1.3  (October 6, 2005)
!!                  Version 1.4  (November 14, 2005)
!!
!!                William W. Hager    and   Hongchao Zhang
!!               hager@math.ufl.edu       hzhang@math.ufl.edu
!!                        Department of Mathematics
!!                          University of Florida
!!                      Gainesville, Florida 32611 USA
!!                           352-392-0281 x 244
!!
!!                   Copyright 2004 by William W. Hager
!!
!!     This program is free software; you can redistribute it and/or   \n
!!     modify it under the terms of the GNU General Public License as  \n
!!     published by the Free Software Foundation; either version 2 of  \n
!!     the License, or (at your option) any later version.             \n
!!     This program is distributed in the hope that it will be useful, \n
!!     but WITHOUT ANY WARRANTY; without even the implied warranty of  \n
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   \n
!!     GNU General Public License for more details.                    \n
!!                                                                     \n
!!     You should have received a copy of the GNU General Public       \n
!!     License along with this program; if not, write to the Free      \n
!!     Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, \n
!!     MA  02110-1301  USA                                             \n
!!                                                                     \n
!!               http://www.math.ufl.edu/~hager/papers/CG
!
!          INPUT:
!
!      (double) grad_tol-- StopRule = T: |g|_infty <= max (grad_tol,
!                                StopFac*initial |g|_infty) [default]
!                          StopRule = F: |g|_infty <= grad_tol(1+|f|)
!
!      (double) x       --starting guess (length n)
!
!      (int)    dim     --problem dimension (also denoted n)
!
!               cg_value--name of cost evaluation subroutine
!                        (external in main program, cg_value(f, x, n)
!                         puts value of cost function at x in f
!                         f is double precision scalar and x is
!                         double precision array of length n)
!
!               cg_grad --name gradient evaluation subroutine
!                        (external in main program, cg_grad (g, x, n)
!                         puts gradient at x in g, g and x are
!                         double precision arrays of length n)
!
!      (double) gnorm   --if the parameter Step in cg_descent.parm is
!                         .true., then gnorm contains the initial step
!                         used at iteration 0 in the line search
!
!      (double) d       --direction (work array, length n)
!
!      (double) g       --gradient (work array, length n)
!
!      (double) xtemp   --work array (work array, length n)
!
!      (double) gtemp   --work array (work array, length n)
!
!          OUTPUT:
!
!      (int)    status  -- 0 (convergence tolerance satisfied)
!                          1 (change in func <= feps*|f|)
!                          2 (total iterations exceeded maxit)
!                          3 (slope always negative in line search)
!                          4 (number secant iterations exceed nsecant)
!                          5 (search direction not a descent direction)
!                          6 (line search fails in initial interval)
!                          7 (line search fails during bisection)
!                          8 (line search fails during interval update)
!
!      (double) gnorm   --max abs component of gradient
!
!      (double) f       --function value at solution
!
!      (double) x       --solution (length n)
!
!      (int)    iter    --number of iterations
!
!      (int)    nfunc   --number of function evaluations
!
!      (int)    ngrad   --number of gradient evaluations
!
!      Note: The file cg_descent.parm must be placed in the directory
!            where the code is run

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE cg_minimization
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module cg_minimization

  implicit none
  save

  integer                               :: n, n5, n6, nf, ng
  integer                               :: info                   !< same as status
  integer                               :: nrestart, nexpand      !< range [0, \f$ \infty \f$), number of grow/shrink allowed in bracket
  integer                               :: nsecant                !< range [0, \f$ \infty \f$), maximum number of secant steps
  integer                               :: maxit                  !< maxmimum number of iterations
  DoublePrecision                          :: delta                  !< range (0, .5), used in the Wolfe conditions
  DoublePrecision                          :: sigma                  !< range [delta, 1), used in the Wolfe conditions
  DoublePrecision                          :: eps                    !< range [0, \f$ \infty \f$, used to compute line search perturbation
  DoublePrecision                          :: gamma                  !< range (0,1), determines when to perform bisection step
  DoublePrecision                          :: rho                    !< range (1, \f$ \infty \f$), growth factor when finding initial interval
  DoublePrecision                          :: tol                    !< range (0, \f$ \infty \f$), convergence tolerance !< convergence tolerance
  DoublePrecision                          :: fpert, f0, eta         !< range (0, \f$ \infty \f$), used in lower bound for beta
  DoublePrecision                          :: ck                     !< average of prior |f|
  DoublePrecision                          :: qdecay
  DoublePrecision                          :: wolfe_hi, wolfe_lo, awolfe_hi
  DoublePrecision                          :: quadcutoff             !< perform QuadStep if relative change in f > QuadCutOff
  DoublePrecision                          :: stopfac                !< used in StopRule
  DoublePrecision                          :: awolfefac              !< used to decide when to switch from Wolfe to AWolfe
  DoublePrecision                          :: zero, feps
  DoublePrecision                          :: psi0                   !< range (0, 1), factor used in very initial starting guess
  DoublePrecision                          :: psi1                   !< range (0, 1), factor previous step multiplied by in QuadStep
  DoublePrecision                          :: psi2                   !< range (1, \f$ \infty \f$), factor previous step is multipled by for startup
  logical                               :: pertrule               !< gives the rule used for the perturbation in f \n  .false. (fpert = eps) \n  .true (fpert = eps*Ck)
  logical                               :: quadok                 !< QuadOK = .true. means that a quadratic step was taken
  logical                               :: quadstep               !< .false. (no quadratic step) \n .true. (use quadratic step)
  logical                               :: printlevel             !< .false. (no printout) \n .true. (print intermediate results)
  logical                               :: printfinal             !< .false. (no printout) \n .true. (print messages, final error)
  logical                               :: stoprule               !< .false. (... <= tol*(1+|f|)) \n .true. (max abs grad <= max (tol, StopFac*initial abs grad))
  logical                               :: awolfe                 !< .false. (use standard Wolfe initially) \n .true. (use approximate + standard Wolfe)
  logical                               :: step                   !< .false. (program computing starting step at iteration 0) \n .true. (user provides starting step in gnorm argument of cg_descent
  logical                               :: debug                  !< .false. (no debugging) \n .true. (check that function values do not increase)

  DoublePrecision                              :: f_best
  DoublePrecision, dimension(:), allocatable   :: x_best

  private
  public  :: cg_descent

contains


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !SUBROUTINE cg_descent
  !
  !> \brief cg_descent
  !!
  !! \todo Detailed description.
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine cg_descent (grad_tol, x, dim, cg_value, cg_grad,  &
       STATUS, gnorm, f, iter, nfunc, ngrad, d, g, xtemp, gtemp)

    implicit none


    integer, intent(IN)                      :: dim         !< (= n) dimensionality of problem
    DoublePrecision, intent(IN)                     :: grad_tol    !< used in stopping rule \n stoprule = .false. \f$ (|g|_\infty <= grad_{tol}(1+|f|)) \f$ \n stoprule = .true. \f$ (|g|_\infty <= max (grad_{tol} , StopFac * initial |g|_\infty)) \f$ [default]
    DoublePrecision, intent(IN OUT)                 :: x(dim)      !< initial guess / for output: solution (length n)
    integer, intent(OUT)                     :: STATUS      !<  0 (convergence tolerance satisfied)\n 1 (change in func <= feps*|f|)\n 2 (total iterations exceeded maxit) \n 3 (slope always negative in line search)\n 4 (number secant iterations exceed nsecant)\n 5 (search direction not a descent direction) \n 6 (line search fails in initial interval) \n 7 (line search fails during bisection)\n 8 (line search fails during interval update)
    DoublePrecision, intent(IN OUT)                 :: gnorm       !< max abs component of gradient
    DoublePrecision, intent(IN OUT)                 :: f           !< function value at solution
    integer, intent(OUT)                     :: iter        !< number of iterations
    integer, intent(OUT)                     :: nfunc       !< number of function evaluations
    integer, intent(OUT)                     :: ngrad       !< number of gradient evaluations
    DoublePrecision, intent(IN OUT)                 :: d(dim)      !< direction (work array, length n)
    DoublePrecision, intent(IN OUT)                 :: g(dim)      !< gradient (work array, length n)
    DoublePrecision, intent(IN OUT)                 :: xtemp(dim)  !< work array (length n)
    DoublePrecision, intent(IN OUT)                 :: gtemp(dim)  !< work array (length n)

    interface

       !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
       ! SUBROUTINE cg_value
       !
       !> \brief cost evaluation subroutine
       !!
       !!(external in main program, cg_value(f, x, n) \n
       !!puts value of cost function at x in f        \n
       !!f is double precision scalar and x is double precision array of length n)
       !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
       subroutine cg_value (f, x, n)
         implicit none
         integer, intent(IN)            :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: f
       end subroutine cg_value

       !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
       !SUBROUTINE cg_grad
       !
       !> \brief gradient evaluation subroutine
       !!
       !!(external in main program, cg_grad (g, x, n) \n
       !!puts gradient at x in g, g and x are double precision arrays of length n)
       !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
       subroutine cg_grad (g, x, n)
         implicit none
         integer, intent(IN)            :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: g(:)
       end subroutine cg_grad

    end interface

    integer                               :: i, j, i1, i2, i3, i4
    DoublePrecision                              :: delta2, eta_sq, qk
    DoublePrecision                              :: ftemp, xnorm, gnorm2
    DoublePrecision                              :: dnorm2, denom
    DoublePrecision                              :: t, t1, t2, t3, t4
    DoublePrecision                              :: dphi, dphi0
    DoublePrecision                              :: alpha, talpha
    DoublePrecision                              :: yk, yk2, ykgk, dkyk
    DoublePrecision                              :: beta
    !---------------------------------------------------------------------------

    allocate( x_best( dim ) )
    f_best = 1.E30

    ! initialize the parameters

    call cg_init (grad_tol, dim)
    if ( step ) then
       alpha = gnorm
    end if
    delta2 = 2.0 * delta - 1.0
    eta_sq = eta * eta
    iter = 0
    ck = 0
    qk = 0

    ! initial function and gradient evaluations, initial direction

    call cg_value (f, x, n)
    if ( f < f_best ) x_best(1:n) = x(1:n)
    if ( f < f_best ) f_best = f
    nf = nf + 1
    call cg_grad (g, x, n)
    ng = ng + 1
    f0 = f + f
    gnorm = zero
    xnorm = zero
    gnorm2 = zero
    do i = 1, n5
       xnorm = max(xnorm, abs(x(i)))
       t = g(i)
       d(i) = -t
       gnorm = max(gnorm, abs(t))
       gnorm2 = gnorm2 + t*t
    end do
    do i = n6, n, 5
       xnorm = max(xnorm, abs(x(i)))
       t = g(i)
       gnorm = max(gnorm, abs(t))
       d(i)   = -t
       j = i + 1
       t1 = g(j)
       d(j) = -t1
       gnorm = max(gnorm, abs(t1))
       xnorm = max(xnorm, abs(x(j)))
       j = i + 2
       t2 = g(j)
       d(j) = -t2
       gnorm = max(gnorm, abs(t2))
       xnorm = max(xnorm, abs(x(j)))
       j = i + 3
       t3 = g(j)
       d(j) = -t3
       gnorm = max(gnorm, abs(t3))
       xnorm = max(xnorm, abs(x(j)))
       j = i + 4
       t4 = g(j)
       d(j) = -t4
       gnorm = max(gnorm, abs(t4))
       xnorm = max(xnorm, abs(x(j)))
       gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
    end do

    if ( stoprule ) then
       tol = max(gnorm*stopfac, tol)
    end if

    if ( printlevel ) then
       write (*, 10) iter, f, gnorm, awolfe
10     format ('iter: ', i5, ' f= ', e14.6,  &
            ' gnorm= ', e14.6, ' AWolfe= ', l2)
    end if

    if ( cg_tol(f, gnorm) ) GO TO 100

    dphi0 = -gnorm2
    if ( .not.step ) then
       alpha = psi0*xnorm/gnorm
       if ( xnorm == zero ) then
          if ( f /= zero ) then
             alpha = psi0*abs(f)/gnorm2
          else
             alpha = 1.0
          end if
       end if
    end if

    ! start the conjugate gradient iteration


    !   alpha starts as old step, ends as initial step for next iteration
    !   f is function value for alpha = 0
    !   QuadOK = .true. means that a quadratic step was taken


    do iter = 1, maxit
       quadok = .false.
       alpha = psi2*alpha
       if ( quadstep ) then
          if ( f /= zero ) then
             t = abs((f-f0)/f)
          else
             t = 1.0
          end if
          if ( t > quadcutoff ) then
             talpha = psi1*alpha
             call cg_step(xtemp, x, d, talpha)
             call cg_value(ftemp, xtemp, n)
             if ( ftemp < f_best ) x_best(1:n) = xtemp(1:n)
             if ( ftemp < f_best ) f_best = ftemp
             nf = nf + 1
             if ( ftemp < f ) then
                denom = 2.0D0*(((ftemp-f)/talpha)-dphi0)
                if ( denom > zero ) then
                   quadok = .true.
                   alpha = -dphi0*talpha/denom
                end if
             end if
          end if
       end if
       f0 = f

       if ( printlevel ) then
          write (*, 20) quadok, alpha, f0, dphi0
20        format ('QuadOK:', l2, ' initial a:',  &
               e14.6,' f0:', e14.6, ' dphi', e14.6)
       end if

       ! parameters in Wolfe and approximiate Wolfe conditions, and in update

       qk = qdecay*qk + 1.
       ck = ck + (abs(f) - ck)/qk

       if ( pertrule ) then
          fpert = f + eps*ck
       else
          fpert = f + eps
       end if

       wolfe_hi = delta*dphi0
       wolfe_lo = sigma*dphi0
       awolfe_hi = delta2*dphi0
       if ( awolfe ) then
          call cg_line (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,  &
               cg_value, cg_grad)
       else
          call cg_linew (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,  &
               cg_value, cg_grad)
       end if

       if ( info > 0 ) GO TO 100

       ! Test for convergence to within machine epsilon
       ! (set feps to zero to remove this test)

       if ( -alpha*dphi0 <= feps*abs(f) ) then
          info = 1
          GO TO 100
       end if

       ! compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g,

       if ( mod(iter, nrestart) /= 0 ) then
          gnorm = zero
          dnorm2 = zero
          yk2 = zero
          ykgk = zero
          do i = 1, n5
             x(i) = xtemp(i)
             t = gtemp(i)
             yk = t - g(i)
             yk2 = yk2 + yk**2
             ykgk = ykgk + yk*t
             g(i) = t
             gnorm = max(gnorm, abs(t))
             dnorm2 = dnorm2 + d(i)**2
          end do
          do i = n6, n, 5
             x(i) = xtemp(i)
             t = gtemp(i)
             yk = t - g(i)
             yk2 = yk2 + yk**2
             ykgk = ykgk + yk*t
             i1 = i + 1
             x(i1) = xtemp(i1)
             t1 = gtemp(i1)
             i2 = i + 2
             x(i2) = xtemp(i2)
             t2 = gtemp(i2)
             i3 = i + 3
             x(i3) = xtemp(i3)
             t3 = gtemp(i3)
             i4 = i + 4
             x(i4) = xtemp(i4)
             t4 = gtemp(i4)
             yk2 = yk2 + (t1-g(i1))**2 + (t2-g(i2))**2  &
                  + (t3-g(i3))**2 + (t4-g(i4))**2
             ykgk = ykgk + (t1-g(i1))*t1 + (t2-g(i2))*t2  &
                  + (t3-g(i3))*t3 + (t4-g(i4))*t4
             g(i) = t
             gnorm = max(gnorm, abs(t))
             g(i1) = t1
             gnorm = max(gnorm, abs(t1))
             g(i2) = t2
             gnorm = max(gnorm, abs(t2))
             g(i3) = t3
             gnorm = max(gnorm, abs(t3))
             g(i4) = t4
             gnorm = max(gnorm, abs(t4))
             dnorm2 = dnorm2 + d(i)**2 + d(i1)**2 + d(i2)**2  &
                  + d(i3)**2 + d(i4)**2
          end do
          if ( cg_tol(f, gnorm) ) GO TO 100
          dkyk = dphi - dphi0
          beta = (ykgk - 2.0*dphi*yk2/dkyk)/dkyk

          !   faster: initialize dnorm2 = gnorm2 at start, then
          !             dnorm2 = gnorm2 + beta**2*dnorm2 - 2.0*beta*dphi
          !             gnorm2 = ||g_{k+1}||^2
          !             dnorm2 = ||d_{k+1}||^2
          !             dpi = g_{k+1}' d_k

          beta = max(beta, -1.0/sqrt(min(eta_sq, gnorm2)*dnorm2))

          !     update search direction d = -g + beta*dold

          gnorm2 = zero
          do i = 1, n5
             t = g(i)
             d(i) = -t + beta*d(i)
             gnorm2 = gnorm2 + t*t
          end do
          do i = n6, n, 5
             d(i) = -g(i) + beta*d(i)
             i1 = i + 1
             d(i1) = -g(i1) + beta*d(i1)
             i2 = i + 2
             d(i2) = -g(i2) + beta*d(i2)
             i3 = i + 3
             d(i3) = -g(i3) + beta*d(i3)
             i4 = i + 4
             d(i4) = -g(i4) + beta*d(i4)
             gnorm2 = gnorm2 + g(i)**2 + g(i1)**2 + g(i2)**2  &
                  + g(i3)**2 + g(i4)**2
          end do
          dphi0 = -gnorm2 + beta*dphi

       else

          !     search direction d = -g

          if ( printlevel ) then
             write (*, *) "RESTART CG"
          end if
          gnorm = zero
          gnorm2 = zero
          do i = 1, n5
             x(i) = xtemp(i)
             t = gtemp(i)
             g(i) = t
             d(i) = -t
             gnorm = max(gnorm, abs(t))
             gnorm2 = gnorm2 + t*t
          end do
          do i = n6, n, 5
             x(i) = xtemp(i)
             t = gtemp(i)
             g(i) = t
             d(i) = -t
             gnorm = max(gnorm, abs(t))
             j = i + 1
             x(j) = xtemp(j)
             t1 = gtemp(j)
             g(j) = t1
             d(j) = -t1
             gnorm = max(gnorm, abs(t1))
             j = i + 2
             x(j) = xtemp(j)
             t2 = gtemp(j)
             g(j) = t2
             d(j) = -t2
             gnorm = max(gnorm, abs(t2))
             j = i + 3
             x(j) = xtemp(j)
             t3 = gtemp(j)
             g(j) = t3
             d(j) = -t3
             gnorm = max(gnorm, abs(t3))
             j = i + 4
             x(j) = xtemp(j)
             t4 = gtemp(j)
             g(j) = t4
             d(j) = -t4
             gnorm = max(gnorm, abs(t4))
             gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
          end do
          if ( cg_tol(f, gnorm) ) GO TO 100
          dphi0 = -gnorm2
       end if
       if ( .not.awolfe ) then
          if ( abs(f-f0) < awolfefac*ck ) then
             awolfe = .true.
          end if
       end if

       if ( printlevel ) then
          write (*, 10) iter, f, gnorm, awolfe
       end if

       if ( debug ) then
          if ( f > f0 + 1.e-10*ck ) then
             write (*, 270)
             write (*, 50) f, f0
50           format (' new value:', e30.16, 'old value:', e30.16)
             stop
          end if
       end if

       if ( dphi0 > zero ) then
          info = 5
          GO TO 100
       end if
    end do
    info = 2
100 nfunc = nf
    ngrad = ng
    STATUS = info
    if ( info > 2 ) then
       gnorm = zero
       do i = 1, n
          x(i) = xtemp(i)
          g(i) = gtemp(i)
          gnorm = max(gnorm, abs(g(i)))
       end do
    end if
    if ( printfinal ) then
       write (6, *) 'Termination status:', STATUS
       if ( STATUS == 0 ) then
          write (6, 200)
       else if ( STATUS == 1 ) then
          write (6, 210)
       else if ( STATUS == 2 ) then
          write (6, 220) maxit
          write (6, 300)
          write (6, 400) grad_tol
       else if ( STATUS == 3 ) then
          write (6, 230)
          write (6, 300)
          write (6, 430)
          write (6, 410)
       else if ( STATUS == 4 ) then
          write (6, 240)
          write (6, 300)
          write (6, 400) grad_tol
       else if ( STATUS == 5 ) then
          write (6, 250)
       else if ( STATUS == 6 ) then
          write (6, 260)
          write (6, 300)
          write (6, 400) grad_tol
          write (6, 410)
          write (6, 420)
       else if ( STATUS == 7 ) then
          write (6, 260)
          write (6, 300)
          write (6, 400) grad_tol
       else if ( STATUS == 8 ) then
          write (6, 260)
          write (6, 300)
          write (6, 400) grad_tol
          write (6, 410)
          write (6, 420)
       end if
       write (6, 500) gnorm
       write (6, *) 'function value:', f
       write (6, *) 'cg iterations:', iter
       write (6, *) 'function evaluations:', nfunc
       write (6, *) 'gradient evaluations:', ngrad
    end if

    if ( f_best < f ) then
       f = f_best
       x(1:dim) = x_best(1:dim)
    end if

    deallocate( x_best )

    return

200 format (' Convergence tolerance for gradient satisfied')
210 format (' Terminating since change in function value <= feps*|f|')
220 format (' Total number of iterations exceed max allow:', i10)
230 format (' Slope always negative in line search')
240 format (' Line search fails, too many secant steps')
250 format (' Search direction not a descent direction')
260 format (' Line search fails')
270 format (' Debugger is on, function value does not improve')
300 format (' Possible causes of this error message:')
400 format ('   - your tolerance (grad_tol = ', d11.4,  &
         ') may be too strict')
410 format ('   - your gradient routine has an error')
420 format ('   - parameter epsilon in cg_descent.parm is too small')
430 format ('   - your cost function has an error')
500 format (' absolute largest component of gradient: ', d11.4)

  end subroutine cg_descent




  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !     PARAMETERS:

  !     delta - range (0, .5), used in the Wolfe conditions
  !     sigma - range [delta, 1), used in the Wolfe conditions
  !     eps - range [0, infty), used to compute line search perturbation
  !     gamma - range (0,1), determines when to perform bisection step
  !     rho   - range (1, infty), growth factor when finding initial interval
  !     eta   - range (0, infty), used in lower bound for beta
  !     psi0  - range (0, 1), factor used in very initial starting guess
  !     psi1  - range (0, 1), factor previous step multiplied by in QuadStep
  !     psi2  - range (1, infty), factor previous step is multipled by for startup
  !     QuadCutOff - perform QuadStep if relative change in f > QuadCutOff
  !     StopFac - used in StopRule
  !     AWolfeFac - used to decide when to switch from Wolfe to AWolfe
  !     restart_fac - range (0, infty) restart cg when iter = n*restart
  !     maxit_fac - range (0, infty) terminate in maxit = maxit_fac*n iterations
  !     feps - stop when -alpha*dphi0 (est. change in value) <= feps*|f|
  !            (feps = 0 removes this test, example: feps = eps*1.e-5
  !             where eps is machine epsilon)
  !     tol   - range (0, infty), convergence tolerance
  !     nexpand - range [0, infty), number of grow/shrink allowed in bracket
  !     nsecant - range [0, infty), maximum number of secant steps
  !     PertRule - gives the rule used for the perturbation in f
  !                 F => fpert = eps
  !                 T => fpert = eps*Ck, Ck is an average of prior |f|
  !     QuadStep- .true. (use quadratic step) .false. (no quadratic step)
  !     PrintLevel- .false. (no printout) .true. (print intermediate results)
  !     PrintFinal- .false. (no printout) .true. (print messages, final error)
  !     StopRule - .true. (max abs grad <= max (tol, StopFac*initial abs grad))
  !                .false. (... <= tol*(1+|f|))
  !     AWolfe - .false. (use standard Wolfe initially)
  !            - .true. (use approximate + standard Wolfe)
  !     Step - .false. (program computing starting step at iteration 0)
  !          - .true. (user provides starting step in gnorm argument of cg_descent
  !     debug - .false. (no debugging)
  !           - .true. (check that function values do not increase)
  !     info  - same as status

  !     DEFAULT PARAMETER VALUES:

  !         delta : 0.1
  !         sigma : 0.9
  !         eps : 1.e-6
  !         gamma : 0.66
  !         rho   : 5.0
  !         restart: 1.0
  !         eta   : 0.01
  !         psi0  : 0.01
  !         psi1  : 0.1
  !         psi2  : 2.0
  !         QuadCutOff: 1.E-12
  !         StopFac: 0.0
  !         AWolfeFac: 1.E-3
  !         tol   : grad_tol
  !         nrestart: n (restart_fac = 1)
  !         maxit : 500*n (maxit_fac = 500)
  !         feps : 0.0
  !         Qdecay : 0.7
  !         nexpand: 50
  !         nsecant: 50
  !         PertRule: .true.
  !         QuadStep: .true.
  !         PrintLevel: .false.
  !         PrintFinal: .true.
  !         StopRule: .true.
  !         AWolfe: .false.
  !         Step: .false.
  !         debug: .false.
  !         info  : 0
  !         feps  : 0.0


  !      (double) grad_tol-- used in stopping rule
  !      (int)    dim     --problem dimension (also denoted n)
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !SUBROUTINE cg_init
  !
  !> \brief Initialisierung
  !!
  !! Initialisierung von cg_descent
  !! \todo Detailed description.
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine cg_init (grad_tol, dim)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(IN)                      :: grad_tol     !< used in stopping rule \n stoprule = .false. (|g|_infty <= grad_tol(1+|f|)) \n stoprule = .true. (|g|_infty <= max (grad_tol,StopFac*initial |g|_infty)) [default]
    integer, intent(IN)                   :: dim          !< (= n) dimensionality of problem

    DoublePrecision                                  :: restart_fac  !< range (0, infty) restart cg when iter = n*restart
    DoublePrecision                                  :: maxit_fac    !< range (0, infty) terminate in maxit = maxit_fac*n iterations
    !---------------------------------------------------------------------------

    n = dim
    tol = grad_tol
    open (10, FILE='cg_descent_f.parm')
    read (10, *) delta
    read (10, *) sigma
    read (10, *) eps
    read (10, *) gamma
    read (10, *) rho
    read (10, *) eta
    read (10, *) psi0
    read (10, *) psi1
    read (10, *) psi2
    read (10, *) quadcutoff
    read (10, *) stopfac
    read (10, *) awolfefac
    read (10, *) restart_fac
    read (10, *) maxit_fac
    read (10, *) feps
    read (10, *) qdecay
    read (10, *) nexpand
    read (10, *) nsecant
    read (10, *) pertrule
    read (10, *) quadstep
    read (10, *) printlevel
    read (10, *) printfinal
    read (10, *) stoprule
    read (10, *) awolfe
    read (10, *) step
    read (10, *) debug
    nrestart = n*restart_fac
    maxit = n*maxit_fac
    zero = 0.0
    info = 0
    n5 = mod(n, 5)
    n6 = n5 + 1
    nf = 0
    ng = 0
    close (10)

  end subroutine cg_init



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !FUNCTION cg_wolfe
  !
  !> \brief check whether the Wolfe or the approximate Wolfe conditions are satisfied
  !      (double) alpha   -- stepsize
  !      (double) f       -- function value associated with stepsize alpha
  !      (double) dphi    -- derivative value associated with stepsize alpha
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  logical function cg_wolfe (alpha, f, dphi)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(IN)             :: alpha                 !< stepsize
    DoublePrecision, intent(IN)             :: f                     !< function value associated with stepsize alpha
    DoublePrecision, intent(IN)             :: dphi                  !< derivative value associated with stepsize alpha
    !---------------------------------------------------------------------------

    if ( dphi >= wolfe_lo ) then

       ! test original Wolfe conditions

       if ( f-f0 <= alpha*wolfe_hi ) then
          cg_wolfe = .true.
          if ( printlevel ) then
             write (*, 10) f, f0, alpha*wolfe_hi, dphi
10           format (' wolfe f:', e14.6, ' f0: ',  &
                  e14.6, e14.6, ' dphi:', e14.6)
          end if
          return

          ! test approximate Wolfe conditions

       else if ( awolfe ) then
          if ( (f <= fpert).and.(dphi <= awolfe_hi) ) then
             cg_wolfe = .true.
             if ( printlevel ) then
                write (*, 20) f, fpert, dphi, awolfe_hi
20              format ('f:', e14.6, ' fpert:', e14.6,  &
                     ' dphi: ', e14.6, ' fappx:', e14.6)
             end if
             return
          end if
       end if
    end if
    cg_wolfe = .false.

  end function cg_wolfe



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !FUNCTION cg_tol
  !
  !> \brief check for convergence of the cg iterations
  !      (double) f       -- function value associated with stepsize
  !      (double) gnorm   -- gradient (infinity) norm
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  logical function cg_tol (f, gnorm)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(IN OUT)         :: f                    !< function value associated with stepsize
    DoublePrecision, intent(IN)             :: gnorm                !< gradient (infinity) norm
    !---------------------------------------------------------------------------

    if ( stoprule ) then
       if ( gnorm <= tol ) then
          cg_tol = .true.
          return
       end if
    else
       if ( gnorm <= tol*(1.0 + abs(f)) ) then
          cg_tol = .true.
          return
       end if
    end if
    cg_tol = .false.

  end function cg_tol



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !FUNCTION cg_dot
  !
  !> \brief compute dot product of x and y, vectors of length n \f$ x \cdot y \f$
  !      (double) x       -- first vector
  !      (double) y       -- second vector
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  DoublePrecision function cg_dot (x, y)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(IN)                      :: x(:)        !< first vector
    DoublePrecision, intent(IN)                      :: y(:)        !< second vector
    DoublePrecision                                  :: t

    integer                               :: i
    !---------------------------------------------------------------------------

    t = zero
    do i = 1, n5
       t = t + x(i)*y(i)
    end do
    do i = n6, n, 5
       t = t + x(i)*y(i) + x(i+1)*y(i+1) + x(i+2)*y(i+2)  &
            + x(i+3)*y(i+3) + x(i+4)*y(i+4)
    end do
    cg_dot = t

  end function cg_dot


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE cg_step
  !
  !> \brief  compute \f$ xtemp = x + alpha \cdot d \f$

  !     (double) xtemp   -- output vector
  !     (double) x       -- initial vector
  !     (double) d       -- search direction vector
  !     (double) alpha   -- stepsize along search direction vector
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine cg_step (xtemp, x, d, alpha)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(OUT)            :: xtemp(:)             !< output vector
    DoublePrecision, intent(IN)             :: x(:)                 !< initial vector
    DoublePrecision, intent(IN)             :: d(:)                 !< search direction vector
    DoublePrecision, intent(IN)             :: alpha                !< stepsize along search direction vector

    integer :: i, j
    !---------------------------------------------------------------------------

    do i = 1, n5
       xtemp(i) = x(i) + alpha*d(i)
    end do
    do i = n6, n, 5
       xtemp(i) = x(i) + alpha*d(i)
       j = i + 1
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 2
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 3
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 4
       xtemp(j) = x(j) + alpha*d(j)
    end do
  end subroutine cg_step



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !      (double) alpha   -- stepsize along search direction vector
  !      (double) phi     -- function value for step alpha
  !      (double) dphi    -- function derivative for step alpha
  !      (double) dphi0   -- function derivative at starting point (alpha = 0)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine cg_line (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp,  &
       cg_value, cg_grad)

    !  use define_cg_variables
    implicit none


    DoublePrecision, intent(IN OUT)         :: alpha
    DoublePrecision, intent(IN OUT)         :: phi
    DoublePrecision, intent(OUT)            :: dphi
    DoublePrecision, intent(IN)             :: dphi0
    DoublePrecision, intent(IN OUT)         :: x(:)
    DoublePrecision, intent(IN OUT)         :: xtemp(:)
    DoublePrecision, intent(IN OUT)         :: d(:)
    DoublePrecision, intent(IN OUT)         :: gtemp(:)

    interface

       subroutine cg_value (f, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: f
       end subroutine cg_value

       subroutine cg_grad (g, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: g(:)
       end subroutine cg_grad

    end interface

    DoublePrecision                                  :: a, dphia, b, dphib, c
    DoublePrecision                                  :: a0, da0, b0, db0
    DoublePrecision                                  :: width, fquad

    integer                               :: ngrow, nshrink, iter, flag
    !---------------------------------------------------------------------------

    call cg_step (xtemp, x, d, alpha)
    call cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)

    ! Find initial interval [a,b] such that dphia < 0, dphib >= 0,
    !        and phia <= phi0 + feps*dabs (phi0)

    a = zero
    dphia = dphi0
    ngrow = 0
    nshrink = 0
    do while ( dphi < zero )
       call cg_value (phi, xtemp, n)
       nf = nf + 1

       ! if quadstep in effect and quadratic conditions hold, check wolfe condition

       if ( quadok ) then
          if ( ngrow == 0 ) fquad = min(phi, f0)
          if ( phi <= fquad ) then
             if ( printlevel ) then
                write (*, 10) alpha, phi, fquad
10              format ('alpha:', e14.6, ' phi:', e14.6,  &
                     ' fquad:', e14.6)
             end if
             if ( cg_wolfe(alpha, phi, dphi) ) return
          end if
       end if
       if ( phi <= fpert ) then
          a = alpha
          dphia = dphi
       else

          ! contraction phase

          b = alpha
          do while ( .true. )
             alpha = .5D0*(a+b)
             nshrink = nshrink + 1
             if ( nshrink > nexpand ) then
                info = 6
                return
             end if
             call cg_step(xtemp, x, d, alpha)
             call cg_grad(gtemp, xtemp, n)
             ng = ng + 1
             dphi = cg_dot(gtemp, d)
             if ( dphi >= zero ) GO TO 100
             call cg_value(phi, xtemp, n)
             nf = nf + 1
             if ( printlevel ) then
                write (6, 20) a, b, alpha, phi, dphi
20              format ('contract, a:', e14.6,  &
                     ' b:', e14.6, ' alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
             end if
             if ( quadok .and. (phi <= fquad) ) then
                if ( cg_wolfe(alpha, phi, dphi) ) return
             end if
             if ( phi <= fpert ) then
                a = alpha
                dphia = dphi
             else
                b = alpha
             end if
          end do
       end if

       ! expansion phase

       ngrow = ngrow + 1
       if ( ngrow > nexpand ) then
          info = 3
          return
       end if
       alpha = rho*alpha
       call cg_step(xtemp, x, d, alpha)
       call cg_grad(gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       if ( printlevel ) then
          write (*, 30) a, alpha, phi, dphi
30        format ('expand,   a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       end if
    end do
100 continue
    b = alpha
    dphib = dphi
    if ( quadok ) then
       call cg_value(phi, xtemp, n)
       nf = nf + 1
       if ( ngrow + nshrink == 0 ) fquad = min(phi, f0)
       if ( phi <= fquad ) then
          if ( cg_wolfe(alpha, phi, dphi) ) return
       end if
    end if
    do iter = 1, nsecant
       if ( printlevel ) then
          write (*, 40) a, b, dphia, dphib
40        format ('secant, a:', e14.6, ' b:', e14.6,  &
               ' da:', e14.6, ' db:', e14.6)
       end if
       width = gamma*(b - a)
       if ( -dphia <= dphib ) then
          alpha = a - (a-b)*(dphia/(dphia-dphib))
       else
          alpha = b - (a-b)*(dphib/(dphia-dphib))
       end if
       c = alpha
       a0 = a
       b0 = b
       da0 = dphia
       db0 = dphib
       flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
            dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
       if ( flag > 0 ) then
          return
       else if ( flag == 0 ) then
          if ( c == a ) then
             if ( dphi > da0 ) then
                alpha = c - (c-a0)*(dphi/(dphi-da0))
             else
                alpha = a
             end if
          else
             if ( dphi < db0 ) then
                alpha = c - (c-b0)*(dphi/(dphi-db0))
             else
                alpha = b
             end if
          end if
          if ( (alpha > a) .and. (alpha < b) ) then
             if ( printlevel ) write (*, *) "2nd secant"
             flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
                  dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
             if ( flag > 0 ) return
          end if
       end if

       !    bisection iteration

       if ( (b-a) >= width ) then
          alpha = .5D0*(b+a)
          if ( printlevel ) write (*, *) "bisection"
          flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
               dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
          if ( flag > 0 ) return
       else
          if ( b <= a ) then
             info = 7
             return
          end if
       end if
    end do
    info = 4

  end subroutine cg_line


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE cg_linew
  !
  !> \brief cg_linew
  !!
  !! This routine is identical to cg_line except that the function
  !! psi (a) = phi (a) - phi (0) - a*delta*dphi(0) is miniminized instead of
  !! the function phi
  !      (double) alpha   -- stepsize along search direction vector
  !      (double) phi     -- function value for step alpha
  !      (double) dphi    -- function derivative for step alpha
  !      (double) dphi0   -- function derivative at starting point (alpha = 0)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine cg_linew (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(IN OUT)         :: alpha
    DoublePrecision, intent(IN OUT)             :: phi
    DoublePrecision, intent(OUT)            :: dphi
    DoublePrecision, intent(IN)             :: dphi0
    DoublePrecision, intent(IN OUT)         :: x(:)
    DoublePrecision, intent(IN OUT)         :: xtemp(:)
    DoublePrecision, intent(IN OUT)         :: d(:)
    DoublePrecision, intent(IN OUT)         :: gtemp(:)

    interface

       subroutine cg_value (f, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: f
       end subroutine cg_value

       subroutine cg_grad (g, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: g(:)
       end subroutine cg_grad

    end interface

    DoublePrecision :: a, dpsia, b, dpsib, c,  &
         a0, da0, b0, db0, width, fquad, psi, dpsi

    integer :: ngrow, nshrink, iter, flag
    !---------------------------------------------------------------------------


    call cg_step (xtemp, x, d, alpha)
    call cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)
    dpsi = dphi - wolfe_hi

    ! Find initial interval [a,b] such that dpsia < 0, dpsib >= 0,
    !        and psia <= phi0 + feps*dabs(phi0)

    a = zero
    dpsia = dphi0 - wolfe_hi
    ngrow = 0
    nshrink = 0
    do while ( dpsi < zero )
       call cg_value(phi, xtemp, n)
       psi = phi - alpha*wolfe_hi

       nf = nf + 1

       ! if quadstep in effect and quadratic conditions hold, check wolfe condition

       if ( quadok ) then
          if ( ngrow == 0 ) fquad = min(phi, f0)
          if ( phi <= fquad ) then
             if ( printlevel ) then
                write (*, 10) alpha, phi, fquad
10              format ('alpha:', e14.6, ' phi:', e14.6,  &
                     ' fquad:', e14.6)
             end if
             if ( cg_wolfe(alpha, phi, dphi) ) return
          end if
       end if
       if ( psi <= fpert ) then
          a = alpha
          dpsia = dpsi
       else

          ! contraction phase

          b = alpha
          do while ( .true. )
             alpha = .5D0*(a+b)
             nshrink = nshrink + 1
             if ( nshrink > nexpand ) then
                info = 6
                return
             end if
             call cg_step (xtemp, x, d, alpha)
             call cg_grad (gtemp, xtemp, n)
             ng = ng + 1
             dphi = cg_dot(gtemp, d)
             dpsi = dphi - wolfe_hi
             if ( dpsi >= zero ) GO TO 100
             call cg_value (phi, xtemp, n)
             psi = phi - alpha*wolfe_hi
             nf = nf + 1
             if ( printlevel ) then
                write (6, 20) a, b, alpha, phi, dphi
20              format ('contract, a:', e14.6,  &
                     ' b:', e14.6, ' alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
             end if
             if ( quadok .and. (phi <= fquad) ) then
                if ( cg_wolfe(alpha, phi, dphi) ) return
             end if
             if ( psi <= fpert ) then
                a = alpha
                dpsia = dpsi
             else
                b = alpha
             end if
          end do
       end if

       ! expansion phase

       ngrow = ngrow + 1
       if ( ngrow > nexpand ) then
          info = 3
          return
       end if
       alpha = rho*alpha
       call cg_step (xtemp, x, d, alpha)
       call cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       dpsi = dphi - wolfe_hi
       if ( printlevel ) then
          write (*, 30) a, alpha, phi, dphi
30        format ('expand,   a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
          write (6, *) "expand, alpha:", alpha, "dphi:", dphi
       end if
    end do
100 continue
    b = alpha
    dpsib = dpsi
    if ( quadok ) then
       call cg_value(phi, xtemp, n)
       nf = nf + 1
       if ( ngrow + nshrink == 0 ) fquad = min(phi, f0)
       if ( phi <= fquad ) then
          if ( cg_wolfe(alpha, phi, dphi) ) return
       end if
    end if
    do iter = 1, nsecant
       if ( printlevel ) then
          write (*, 40) a, b, dpsia, dpsib
40        format ('secant, a:', e14.6, ' b:', e14.6,  &
               ' da:', e14.6, ' db:', e14.6)
       end if
       width = gamma*(b - a)
       if ( -dpsia <= dpsib ) then
          alpha = a - (a-b)*(dpsia/(dpsia-dpsib))
       else
          alpha = b - (a-b)*(dpsib/(dpsia-dpsib))
       end if
       c = alpha
       a0 = a
       b0 = b
       da0 = dpsia
       db0 = dpsib
       flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
            phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
       if ( flag > 0 ) then
          return
       else if ( flag == 0 ) then
          if ( c == a ) then
             if ( dpsi > da0 ) then
                alpha = c - (c-a0)*(dpsi/(dpsi-da0))
             else
                alpha = a
             end if
          else
             if ( dpsi < db0 ) then
                alpha = c -(c-b0)*(dpsi/(dpsi-db0))
             else
                alpha = b
             end if
          end if
          if ( (alpha > a) .and. (alpha < b) ) then
             if ( printlevel ) write (*, *) "2nd secant"
             flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
                  phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
             if ( flag > 0 ) return
          end if
       end if

       !    bisection iteration

       if ( (b-a) >= width ) then
          alpha = .5D0*(b+a)
          if ( printlevel ) write (*, *) "bisection"
          flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
               phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
          if ( flag > 0 ) return
       else
          if ( b <= a ) then
             info = 7
             return
          end if
       end if
    end do
    info = 4

  end subroutine cg_linew


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! FUNCTION cg_update
  !
  !> \brief cg_update
  !!
  !!        returns 1 if Wolfe condition is satisfied or too many iterations
  !!        returns  0 if the interval updated successfully
  !!        returns -1 if search done

  !      (double) a       -- left side of bracketting interval
  !      (double) dphia   -- derivative at a
  !      (double) b       -- right side of bracketting interval
  !      (double) dphib   -- derivative at b
  !      (double) alpha   -- trial step (between a and b)
  !      (double) phi     -- function value at alpha (returned)
  !      (double) dphi    -- function derivative at alpha (returned)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  integer function cg_update (a, dphia, b, dphib, alpha, phi,  &
       dphi, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(OUT)            :: a                    !< left side of bracketting interval
    DoublePrecision, intent(OUT)            :: dphia                !< derivative at a
    DoublePrecision, intent(OUT)            :: b                    !< right side of bracketting interval
    DoublePrecision, intent(OUT)            :: dphib                !< derivative at b
    DoublePrecision, intent(IN OUT)         :: alpha                !< trial step (between a and b)
    DoublePrecision, intent(IN OUT)         :: phi                  !< function value at alpha (returned)
    DoublePrecision, intent(OUT)            :: dphi                 !< function derivative at alpha (returned)
    DoublePrecision, intent(IN OUT)         :: x(:)                 !< current iterate
    DoublePrecision, intent(IN OUT)         :: xtemp(:)             !< x + alpha*d
    DoublePrecision, intent(IN OUT)         :: d(:)                 !< current search direction
    DoublePrecision, intent(IN OUT)         :: gtemp(:)             !< gradient at x + alpha*d

    interface

       subroutine cg_value (f, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: f
       end subroutine cg_value

       subroutine cg_grad (g, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: g(:)
       end subroutine cg_grad

    end interface


    integer :: nshrink
    !---------------------------------------------------------------------------

    call cg_step (xtemp, x, d, alpha)
    call cg_value (phi, xtemp, n)
    nf = nf + 1
    call cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)
    if ( printlevel ) then
       write (*, 10) alpha, phi, dphi
10     format ('update alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
    end if
    cg_update = 0
    if ( cg_wolfe(alpha, phi, dphi) ) then
       cg_update = 1
       GO TO 110
    end if
    if ( dphi >= zero ) then
       b = alpha
       dphib = dphi
       GO TO 110
    else
       if ( phi <= fpert ) then
          a = alpha
          dphia = dphi
          GO TO 110
       end if
    end if
    nshrink = 0
    b = alpha
    do while ( .true. )
       alpha = .5D0*(a+b)
       nshrink = nshrink + 1
       if ( nshrink > nexpand ) then
          info = 8
          cg_update = 1
          GO TO 110
       end if
       call cg_step (xtemp, x, d, alpha)
       call cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       call cg_value (phi, xtemp, n)
       nf = nf + 1
       if ( printlevel ) then
          write (6, 20) a, alpha, phi, dphi
20        format ('contract, a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       end if
       if ( cg_wolfe(alpha, phi, dphi) ) then
          cg_update = 1
          GO TO 110
       end if
       if ( dphi >= zero ) then
          b = alpha
          dphib = dphi
          GO TO 100
       end if
       if ( phi <= fpert ) then
          if ( printlevel ) then
             write (6, *) "update a:", alpha, "dphia:", dphi
          end if
          a = alpha
          dphia = dphi
       else
          b = alpha
       end if
    end do
100 continue
    cg_update = -1
110 continue
    if ( printlevel ) then
       write (*, 200) a, b, dphia, dphib, cg_update
200    format ('UP a:', e14.6, ' b:', e14.6,  &
            ' da:', e14.6, ' db:', e14.6, ' up:', i2)
    end if

  end function cg_update


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  This routine is identical to cg_update except that the function
  !  psi (a) = phi(a) - phi(0) - a*delta*dphi(0) is miniminized instead of
  !  the function phi

  ! update returns 1 if Wolfe condition is satisfied or too many iterations
  !        returns  0 if the interval updated successfully
  !        returns -1 if search done

  !      (double) a       -- left side of bracketting interval
  !      (double) dpsia   -- derivative at a
  !      (double) b       -- right side of bracketting interval
  !      (double) dpsib   -- derivative at b
  !      (double) alpha   -- trial step (between a and b)
  !      (double) phi     -- function value at alpha(returned)
  !      (double) dphi    -- derivative of phi at alpha(returned)
  !      (double) dpsi    -- derivative of psi at alpha(returned)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  integer function cg_updatew (a, dpsia, b, dpsib, alpha, phi, dphi,  &
       dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    DoublePrecision, intent(OUT)            :: a
    DoublePrecision, intent(OUT)            :: dpsia
    DoublePrecision, intent(OUT)            :: b
    DoublePrecision, intent(OUT)            :: dpsib
    DoublePrecision, intent(IN OUT)         :: alpha
    DoublePrecision, intent(IN OUT)         :: phi
    DoublePrecision, intent(OUT)            :: dphi
    DoublePrecision, intent(OUT)            :: dpsi
    DoublePrecision, intent(IN OUT)         :: x(:)
    DoublePrecision, intent(IN OUT)         :: xtemp(:)
    DoublePrecision, intent(IN OUT)         :: d(:)
    DoublePrecision, intent(IN OUT)         :: gtemp(:)

    interface

       subroutine cg_value (f, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: f
       end subroutine cg_value

       subroutine cg_grad (g, x, n)
         implicit none
         integer, intent(IN)        :: n
         DoublePrecision, intent(IN)           :: x(:)
         DoublePrecision, intent(IN OUT)       :: g(:)
       end subroutine cg_grad

    end interface

    DoublePrecision :: psi

    integer :: nshrink
    !---------------------------------------------------------------------------

    call cg_step (xtemp, x, d, alpha)
    call cg_value (phi, xtemp, n)
    psi = phi - alpha*wolfe_hi
    nf = nf + 1
    call cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot (gtemp, d)
    dpsi = dphi - wolfe_hi
    if ( printlevel ) then
       write (*, 10) alpha, psi, dpsi
10     format ('update alpha:', e14.6, ' psi:', e14.6, ' dpsi:', e14.6)
    end if
    cg_updatew = 0
    if ( cg_wolfe(alpha, phi, dphi) ) then
       cg_updatew = 1
       GO TO 110
    end if
    if ( dpsi >= zero ) then
       b = alpha
       dpsib = dpsi
       GO TO 110
    else
       if ( psi <= fpert ) then
          a = alpha
          dpsia = dpsi
          GO TO 110
       end if
    end if
    nshrink = 0
    b = alpha
    do while ( .true. )
       alpha = .5D0*(a+b)
       nshrink = nshrink + 1
       if ( nshrink > nexpand ) then
          info = 8
          cg_updatew = 1
          GO TO 110
       end if
       call cg_step (xtemp, x, d, alpha)
       call cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       dpsi = dphi - wolfe_hi
       call cg_value (phi, xtemp, n)
       psi = phi - alpha*wolfe_hi
       nf = nf + 1
       if ( printlevel ) then
          write (6, 20) a, alpha, phi, dphi
20        format ('contract, a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       end if
       if ( cg_wolfe(alpha, phi, dphi) ) then
          cg_updatew = 1
          GO TO 110
       end if
       if ( dpsi >= zero ) then
          b = alpha
          dpsib = dpsi
          GO TO 100
       end if
       if ( psi <= fpert ) then
          if ( printlevel ) then
             write (6, *) "update a:", alpha, "dpsia:", dpsi
          end if
          a = alpha
          dpsia = dpsi
       else
          b = alpha
       end if
    end do
100 continue
    cg_updatew = -1
110 continue
    if ( printlevel ) then
       write (*, 200) a, b, dpsia, dpsib, cg_updatew
200    format ('UP a:', e14.6, ' b:', e14.6,  &
            ' da:', e14.6, ' db:', e14.6, ' up:', i2)
    end if

  end function cg_updatew
  ! Version 1.2 Changes:

  !   1. Fix problem with user specified initial step (overwriting step)
  !   2. Change dphi to dpsi at lines 1228 and 1234 in cg_lineW
  !   3. Add comment about how to compute dnorm2 by an update of previous dnorm2
  !   4. In comment statements for cg_lineW and cg_updateW, insert "delta"
  !      in definition of psi (a)
  !   5. In dimension statements, change "(1)" to "(*)"

  ! Version 1.3 Changes:
  !   1. Remove extraneous write in line 985 (same thing written out twice)
  !   2. Remove the parameter theta from cg_descent.parm and from the code
  !      (we use theta = .5 in the cg_update)

  ! Version 1.4 Change:
  !   1. The variable dpsi needs to be included in the argument list for
  !      subroutine updateW (update of a Wolfe line search)

end module cg_minimization
