!> \file levenberg_marquardt.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!> \brief Levenberg_Marquardt Solver
!!
!! MINPACK routines which are used by both LMDIF & LMDER
!! 25 October 2001:
!!    Changed INTENT of iflag in several places to IN OUT.
!!    Changed INTENT of fvec to IN OUT in user routine FCN.
!!    Removed arguments diag and qtv from LMDIF & LMDER.
!!    Replaced several DO loops with array operations.
!! amiller @ bigpond.net.au
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module Levenberg_Marquardt

  use PARAMETERS, only: dp, machine_eps
  implicit none

  private
  public                                :: lmdif1, lmdif, lmder1, lmder, enorm

contains


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine lmdif1
  !> \brief Minimization of the sum of the squares
  !!
  !!  The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear
  !!  functions in n variables by a modification of the Levenberg-Marquardt
  !!  algorithm.  This is done by using the more general least-squares
  !!  solver lmdif.  The user must provide a subroutine which calculates the
  !!  functions.  The jacobian is then calculated by a forward-difference
  !!  approximation.
  !!
  !!  the subroutine statement is
  !!
  !!    subroutine lmdif1(fcn, m, n, x, fvec, tol, epsfcn, info, iwa)
  !!    subprograms called
  !!
  !!    user-supplied ...... fcn
  !!
  !!    minpack-supplied ... lmdif
  !!
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !!
  !! Code converted using TO_F90 by Alan Miller
  !! Date: 1999-12-11  Time: 00:51:44
  !! N.B. Arguments WA & LWA have been removed.
  !
  !  where
  !
  !    fcn is the name of the user-supplied subroutine which calculates
  !      the functions.  fcn must be declared in an external statement in the
  !      user calling program, and should be written as follows.

  !      subroutine fcn(m, n, x, fvec, iflag)
  !      integer m, n, iflag
  !      REAL x(n), fvec(m)
  !      ----------
  !      calculate the functions at x and return this vector in fvec.
  !      ----------
  !      return
  !      end
  !
  !      the value of iflag should not be changed by fcn unless
  !      the user wants to terminate execution of lmdif1.
  !      In this case set iflag to a negative integer.
  !
  !    m is a positive integer input variable set to the number of functions.
  !
  !    n is a positive integer input variable set to the number of variables.
  !      n must not exceed m.
  !
  !    x is an array of length n.  On input x must contain an initial estimate
  !      of the solution vector.  On output x contains the final estimate of
  !      the solution vector.
  !
  !    fvec is an output array of length m which contains
  !      the functions evaluated at the output x.
  !
  !    tol is a nonnegative input variable.  Termination occurs when the
  !      algorithm estimates either that the relative error in the sum of
  !      squares is at most tol or that the relative error between x and the
  !      solution is at most tol.
  !
  !    factor is a positive input variable used in determining the initial step
  !      bound.  This bound is set to the product of factor and the euclidean
  !      norm of diag*x if nonzero, or else to factor itself.  In most cases
  !      factor should lie in the interval (.1,100.). 100. is a generally
  !      recommended value.
  !
  !    info is an integer output variable.  If the user has terminated execution,
  !      info is set to the (negative) value of iflag.  See description of fcn.
  !      Otherwise, info is set as follows.
  !
  !      info = 0  improper input parameters.
  !
  !      info = 1  algorithm estimates that the relative error
  !                in the sum of squares is at most tol.
  !
  !      info = 2  algorithm estimates that the relative error
  !                between x and the solution is at most tol.
  !
  !      info = 3  conditions for info = 1 and info = 2 both hold.
  !
  !      info = 4  fvec is orthogonal to the columns of the
  !                jacobian to machine precision.
  !
  !      info = 5  number of calls to fcn has reached or exceeded 200*(n+1).
  !
  !      info = 6  tol is too small. no further reduction in
  !                the sum of squares is possible.
  !
  !      info = 7  tol is too small.  No further improvement in
  !                the approximate solution x is possible.
  !
  !    iwa is an integer work array of length n.
  !
  !    wa is a work array of length lwa.
  !
  !    lwa is a positive integer input variable not less than m*n+5*n+m.
  !
  !  subprograms called
  !
  !    user-supplied ...... fcn
  !
  !    minpack-supplied ... lmdif
  !
  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 1999-12-11  Time: 00:51:44
  ! N.B. Arguments WA & LWA have been removed.
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine lmdif1(fcn, m, n, x, fvec, tol, epsfcn, factor, info, iwa)

    integer, intent(IN)                 :: m         !< positive integer input variable, set to number of functions
    integer, intent(IN)                 :: n         !< positive integer input variable, set to number of variables. n must not exceed m
    real(dp), intent(IN OUT)            :: x(:)      !< array of length n. On input x must contain initial estimate of solution.  On output x contains the final solution
    real(dp), intent(OUT)               :: fvec(:)   !< output array, length m which contains the functions evaluated at the output x
    real(dp), intent(IN)                :: tol       !< nonnegative input variable.  Termination for rel. error (sum of squares) < tol or rel. error between x and real solution x0 is < tol
    real(dp), intent(IN OUT)            :: epsfcn
    real(dp), intent(IN)                :: factor
    integer, intent(OUT)                :: info      !< integer output variable.  See details in levenberg_marquardt.f90
    integer, intent(OUT)                :: iwa(:)    !< integer work array of length n

    ! EXTERNAL fcn

    interface
       subroutine fcn(m, n, x, fvec, iflag)
         implicit none
         integer, parameter             :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
         integer, intent(IN)            :: m, n
         real(dp), intent(IN)           :: x(:)
         real(dp), intent(IN OUT)       :: fvec(:)
         integer, intent(IN OUT)        :: iflag
       end subroutine fcn
    end interface

    integer                             :: maxfev, mode, nfev, nprint
    real(dp)                            :: ftol, gtol, xtol, fjac(m,n)
    ! REAL(DP), PARAMETER                 :: factor = 100.0_dp, zero = 0.0_dp
    real(dp), parameter                 :: zero = 0.0_dp


    info = 0

    ! check the input parameters for errors.

    if (n <= 0 .or. m < n .or. tol < zero) return

    ! call lmdif.

    maxfev = 200*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    ! epsfcn = zero
    mode = 1
    nprint = 0
    call lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,   &
         mode, factor, nprint, info, nfev, fjac, iwa)
    if (info == 8) info = 4

  end subroutine lmdif1



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine lmdif
  !> \brief Minimization of the sum of the squares
  !!
  !!  The purpose of lmdif is to minimize the sum of the squares of m nonlinear
  !!  functions in n variables by a modification of the Levenberg-Marquardt
  !!  algorithm.  The user must provide a subroutine which calculates the
  !!  functions.  The jacobian is then calculated by a forward-difference
  !!  approximation.
  !!  the subroutine statement is
  !!    subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
  !!                     diag, mode, factor, nprint, info, nfev, fjac,
  !!                     ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
  !! N.B. 7 of these arguments have been removed in this version.
  !!  subprograms called
  !!    user-supplied ...... fcn
  !!    minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
  !!    fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod \n
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more \n
  !! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.
  !  where

  !    fcn is the name of the user-supplied subroutine which calculates the
  !      functions.  fcn must be declared in an external statement in the user
  !      calling program, and should be written as follows.

  !      subroutine fcn(m, n, x, fvec, iflag)
  !      integer m, n, iflag
  !      REAL x(:), fvec(m)
  !      ----------
  !      calculate the functions at x and return this vector in fvec.
  !      ----------
  !      return
  !      end

  !      the value of iflag should not be changed by fcn unless
  !      the user wants to terminate execution of lmdif.
  !      in this case set iflag to a negative integer.

  !    m is a positive integer input variable set to the number of functions.

  !    n is a positive integer input variable set to the number of variables.
  !      n must not exceed m.

  !    x is an array of length n.  On input x must contain an initial estimate
  !      of the solution vector.  On output x contains the final estimate of the
  !      solution vector.

  !    fvec is an output array of length m which contains
  !      the functions evaluated at the output x.

  !    ftol is a nonnegative input variable.  Termination occurs when both the
  !      actual and predicted relative reductions in the sum of squares are at
  !      most ftol.  Therefore, ftol measures the relative error desired
  !      in the sum of squares.

  !    xtol is a nonnegative input variable.  Termination occurs when the
  !      relative error between two consecutive iterates is at most xtol.
  !      Therefore, xtol measures the relative error desired in the approximate
  !      solution.

  !    gtol is a nonnegative input variable.  Termination occurs when the cosine
  !      of the angle between fvec and any column of the jacobian is at most
  !      gtol in absolute value.  Therefore, gtol measures the orthogonality
  !      desired between the function vector and the columns of the jacobian.

  !    maxfev is a positive integer input variable.  Termination occurs when the
  !      number of calls to fcn is at least maxfev by the end of an iteration.

  !    epsfcn is an input variable used in determining a suitable step length
  !      for the forward-difference approximation.  This approximation assumes
  !      that the relative errors in the functions are of the order of epsfcn.
  !      If epsfcn is less than the machine precision, it is assumed that the
  !      relative errors in the functions are of the order of the machine
  !      precision.

  !    diag is an array of length n.  If mode = 1 (see below), diag is
  !      internally set.  If mode = 2, diag must contain positive entries that
  !      serve as multiplicative scale factors for the variables.

  !    mode is an integer input variable.  If mode = 1, the variables will be
  !      scaled internally.  If mode = 2, the scaling is specified by the input
  !      diag. other values of mode are equivalent to mode = 1.

  !    factor is a positive input variable used in determining the initial step
  !      bound.  This bound is set to the product of factor and the euclidean
  !      norm of diag*x if nonzero, or else to factor itself.  In most cases
  !      factor should lie in the interval (.1,100.). 100. is a generally
  !      recommended value.

  !    nprint is an integer input variable that enables controlled printing of
  !      iterates if it is positive.  In this case, fcn is called with iflag = 0
  !      at the beginning of the first iteration and every nprint iterations
  !      thereafter and immediately prior to return, with x and fvec available
  !      for printing.  If nprint is not positive, no special calls
  !      of fcn with iflag = 0 are made.

  !    info is an integer output variable.  If the user has terminated
  !      execution, info is set to the (negative) value of iflag.
  !      See description of fcn.  Otherwise, info is set as follows.

  !      info = 0  improper input parameters.

  !      info = 1  both actual and predicted relative reductions
  !                in the sum of squares are at most ftol.

  !      info = 2  relative error between two consecutive iterates <= xtol.

  !      info = 3  conditions for info = 1 and info = 2 both hold.

  !      info = 4  the cosine of the angle between fvec and any column of
  !                the Jacobian is at most gtol in absolute value.

  !      info = 5  number of calls to fcn has reached or exceeded maxfev.

  !      info = 6  ftol is too small. no further reduction in
  !                the sum of squares is possible.

  !      info = 7  xtol is too small. no further improvement in
  !                the approximate solution x is possible.

  !      info = 8  gtol is too small. fvec is orthogonal to the
  !                columns of the jacobian to machine precision.

  !    nfev is an integer output variable set to the number of calls to fcn.

  !    fjac is an output m by n array. the upper n by n submatrix
  !      of fjac contains an upper triangular matrix r with
  !      diagonal elements of nonincreasing magnitude such that

  !             t     t           t
  !            p *(jac *jac)*p = r *r,

  !      where p is a permutation matrix and jac is the final calculated
  !      Jacobian.  Column j of p is column ipvt(j) (see below) of the
  !      identity matrix. the lower trapezoidal part of fjac contains
  !      information generated during the computation of r.

  !    ldfjac is a positive integer input variable not less than m
  !      which specifies the leading dimension of the array fjac.

  !    ipvt is an integer output array of length n.  ipvt defines a permutation
  !      matrix p such that jac*p = q*r, where jac is the final calculated
  !      jacobian, q is orthogonal (not stored), and r is upper triangular
  !      with diagonal elements of nonincreasing magnitude.
  !      Column j of p is column ipvt(j) of the identity matrix.

  !    qtf is an output array of length n which contains
  !      the first n elements of the vector (q transpose)*fvec.

  !    wa1, wa2, and wa3 are work arrays of length n.

  !    wa4 is a work array of length m.

  !  subprograms called

  !    user-supplied ...... fcn

  !    minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac

  !    fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
  !
  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !
  ! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
       mode, factor, nprint, info, nfev, fjac, ipvt)

    integer, intent(IN)                 :: m         !< positive integer input variable, set to number of functions
    integer, intent(IN)                 :: n         !< positive integer input variable, set to number of variables. n must not exceed m
    real(dp), intent(IN OUT)            :: x(:)      !< array of length n. On input x must contain initial estimate of solution.  On output x contains the final solution
    real(dp), intent(OUT)               :: fvec(:)   !< output array, length m which contains the functions evaluated at the output x
    real(dp), intent(IN)                :: ftol      !< nonnegative input variable.  Termination when actual and predicted relative reductions in the sum of squares are < ftol.  Therefore, ftol measures the relative error desired in the sum of squares
    real(dp), intent(IN)                :: xtol      !< nonnegative input variable.  Termination when relative error between two iterates is < xtol. Therefore, xtol measures relative error desired in the approximate solution.
    real(dp), intent(IN OUT)            :: gtol      !< nonnegative input variable.  Termination occurs when the cosine  of the angle between fvec and any column of the jacobian is at most  gtol in absolute value.  Therefore, gtol measures the orthogonality  desired between the function vector and the columns of the jacobian.
    integer, intent(IN OUT)             :: maxfev    !< positive integer input variable.  Termination occurs when the  number of calls to fcn is at least maxfev by the end of an iteration.
    real(dp), intent(IN OUT)            :: epsfcn    !< input variable used in determining a suitable step length for the forward-difference approximation.  This approximation assumes  that the relative errors in the functions are of the order of epsfcn.  If epsfcn is less than the machine precision, it is assumed that the  relative errors in the functions are of the order of the machine  precision.
    integer, intent(IN)                 :: mode      !< integer input variable.  If mode = 1, the variables will be scaled internally.  If mode = 2, the scaling is specified by the input diag. other values of mode are equivalent to mode = 1.
    real(dp), intent(IN)                :: factor    !< positive input variable used in determining the initial step bound.  This bound is set to the product of factor and the euclidean norm of diag*x if nonzero, or else to factor itself.  In most cases  factor should lie in the interval (.1,100.). 100. is a generally recommended value.
    integer, intent(IN)                 :: nprint    !< integer input variable that enables controlled printing of iterates if it is positive.  In this case, fcn is called with iflag = 0 at the beginning of the first iteration and every nprint iterations  thereafter and immediately prior to return, with x and fvec available for printing.  If nprint is not positive, no special calls of fcn with iflag = 0 are made.
    integer, intent(OUT)                :: info      !< integer output variable.  See details in levenberg_marquardt.f90
    integer, intent(OUT)                :: nfev      !< integer output variable set to the number of calls to fcn.
    real(dp), intent(OUT)               :: fjac(:,:) !< output m by n array. the upper n by n submatrix of fjac contains an upper triangular matrix r withdiagonal elements of nonincreasing magnitude such that \f$ p^t *(jac^t *jac)*p = r^t *r \f$ where p is a permutation matrix and jac is the final calculated Jacobian.  Column j of p is column ipvt(j) (see below) of the identity matrix. the lower trapezoidal part of fjac contains information generated during the computation of r.
    integer, intent(OUT)                :: ipvt(:)   !< integer output array of length n.  ipvt defines a permutation matrix p such that jac*p = q*r, where jac is the final calculated jacobian, q is orthogonal (not stored), and r is upper triangular  with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.

    ! EXTERNAL fcn

    interface
       subroutine fcn(m, n, x, fvec, iflag)
         implicit none
         integer, parameter             :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
         integer, intent(IN)            :: m, n
         real(dp), intent(IN)           :: x(:)
         real(dp), intent(IN OUT)       :: fvec(:)
         integer, intent(IN OUT)        :: iflag
       end subroutine fcn
    end interface

    integer                             :: i, iflag, iter, j, l
    real(dp)                            :: actred, delta, dirder
    real(dp)                            :: epsmch, fnorm, fnorm1, gnorm
    real(dp)                            :: par, pnorm, prered, ratio
    real(dp)                            :: sum, temp, temp1, temp2, xnorm
    real(dp)                            :: diag(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m)
    real(dp), parameter                 :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp
    real(dp), parameter                 :: p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp
    real(dp), parameter                 :: zero = 0.0_dp

    ! epsmch is the machine precision.

    epsmch = epsilon(zero)

    info = 0
    iflag = 0
    nfev = 0

    ! check the input parameters for errors.

    if (n <= 0 .or. m < n .or. ftol < zero .or. xtol < zero .or. gtol < zero  &
         .or. maxfev <= 0 .or. factor <= zero) GO TO 300
    if (mode /= 2) GO TO 20
    do  j = 1, n
       if (diag(j) <= zero) GO TO 300
    end do

    ! evaluate the function at the starting point and calculate its norm.

20  iflag = 1
    call fcn(m, n, x, fvec, iflag)
    nfev = 1
    if (iflag < 0) GO TO 300
    fnorm = enorm(m, fvec)

    ! initialize levenberg-marquardt parameter and iteration counter.

    par = zero
    iter = 1

    ! beginning of the outer loop.

    ! calculate the jacobian matrix.

30  iflag = 2
    call fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
    nfev = nfev + n
    if (iflag < 0) GO TO 300

    ! If requested, call fcn to enable printing of iterates.

    if (nprint <= 0) GO TO 40
    iflag = 0
    if (mod(iter-1,nprint) == 0) call fcn(m, n, x, fvec, iflag)
    if (iflag < 0) GO TO 300

    ! Compute the qr factorization of the jacobian.

40  call qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

    ! On the first iteration and if mode is 1, scale according
    ! to the norms of the columns of the initial jacobian.

    if (iter /= 1) GO TO 80
    if (mode == 2) GO TO 60
    do  j = 1, n
       diag(j) = wa2(j)
       if ( abs( wa2(j) ) < machine_eps ) diag(j) = one
    end do

    ! On the first iteration, calculate the norm of the scaled x
    ! and initialize the step bound delta.

60  wa3(1:n) = diag(1:n)*x(1:n)
    xnorm = enorm(n, wa3)
    delta = factor*xnorm
    if ( abs( delta ) < machine_eps ) delta = factor

    ! Form (q transpose)*fvec and store the first n components in qtf.

80  wa4(1:m) = fvec(1:m)
    do  j = 1, n
       if ( abs( fjac(j,j) ) < machine_eps ) GO TO 120
       sum = DOT_product( fjac(j:m,j), wa4(j:m) )
       temp = -sum/fjac(j,j)
       do  i = j, m
          wa4(i) = wa4(i) + fjac(i,j)*temp
       end do
120    fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
    end do

    ! compute the norm of the scaled gradient.

    gnorm = zero
    if ( abs( fnorm ) < machine_eps ) GO TO 170
    do  j = 1, n
       l = ipvt(j)
       if ( abs( wa2(l) ) < machine_eps ) cycle
       sum = zero
       do  i = 1, j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
       end do
       gnorm = max(gnorm, abs(sum/wa2(l)))
    end do

    ! test for convergence of the gradient norm.

170 if (gnorm <= gtol) info = 4
    if (info /= 0) GO TO 300

    ! rescale if necessary.

    if (mode == 2) GO TO 200
    do  j = 1, n
       diag(j) = max(diag(j), wa2(j))
    end do

    ! beginning of the inner loop.

    !   determine the Levenberg-Marquardt parameter.

200 call lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

    !   store the direction p and x + p. calculate the norm of p.

    do  j = 1, n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    end do
    pnorm = enorm(n, wa3)

    !   on the first iteration, adjust the initial step bound.

    if (iter == 1) delta = min(delta, pnorm)

    !   evaluate the function at x + p and calculate its norm.

    iflag = 1
    call fcn(m, n, wa2, wa4, iflag)
    nfev = nfev + 1
    if (iflag < 0) GO TO 300
    fnorm1 = enorm(m, wa4)

    !   compute the scaled actual reduction.

    actred = -one
    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

    !   Compute the scaled predicted reduction and
    !   the scaled directional derivative.

    do  j = 1, n
       wa3(j) = zero
       l = ipvt(j)
       temp = wa1(l)
       do  i = 1, j
          wa3(i) = wa3(i) + fjac(i,j)*temp
       end do
    end do
    temp1 = enorm(n,wa3)/fnorm
    temp2 = (sqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)

    !   compute the ratio of the actual to the predicted reduction.

    ratio = zero
    if ( abs( prered ) > machine_eps ) ratio = actred/prered

    !   update the step bound.

    if (ratio <= p25) then
       if (actred >= zero) temp = p5
       if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
       if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
       delta = temp*min(delta,pnorm/p1)
       par = par/temp
    else
       if ( abs( par ) > machine_eps .and. ratio < p75 ) GO TO 260
       delta = pnorm/p5
       par = p5*par
    end if

    !   test for successful iteration.

260 if (ratio < p0001) GO TO 290

    !   successful iteration. update x, fvec, and their norms.

    do  j = 1, n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
    end do
    fvec(1:m) = wa4(1:m)
    xnorm = enorm(n, wa2)
    fnorm = fnorm1
    iter = iter + 1

    !   tests for convergence.

290 if (abs(actred) <= ftol .and. prered <= ftol .and. p5*ratio <= one) info = 1
    if (delta <= xtol*xnorm) info = 2
    if (abs(actred) <= ftol .and. prered <= ftol  &
         .and. p5*ratio <= one .and. info == 2) info = 3
    if (info /= 0) GO TO 300

    !   tests for termination and stringent tolerances.

    if (nfev >= maxfev) info = 5
    if (abs(actred) <= epsmch .and. prered <= epsmch  &
         .and. p5*ratio <= one) info = 6
    if (delta <= epsmch*xnorm) info = 7
    if (gnorm <= epsmch) info = 8
    if (info /= 0) GO TO 300

    !   end of the inner loop. repeat if iteration unsuccessful.

    if (ratio < p0001) GO TO 200

    ! end of the outer loop.

    GO TO 30

    ! termination, either normal or user imposed.

300 if (iflag < 0) info = iflag
    iflag = 0
    if (nprint > 0) call fcn(m, n, x, fvec, iflag)

  end subroutine lmdif



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine lmder1
  !> \brief minimize the sum of the squares
  !!
  !!  The purpose of lmder1 is to minimize the sum of the squares of
  !!  m nonlinear functions in n variables by a modification of the
  !!  levenberg-marquardt algorithm.  This is done by using the more
  !!  general least-squares solver lmder.  The user must provide a
  !!  subroutine which calculates the functions and the jacobian.
  !!  the subroutine statement is \n
  !!    subroutine lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)\n
  !!  subprograms called \n
  !!
  !!    user-supplied ...... fcn\n
  !!
  !!    minpack-supplied ... lmder \n
  !!
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !!
  !!  N.B. Arguments LDFJAC, WA & LWA have been removed.
  !  where

  !    fcn is the name of the user-supplied subroutine which
  !      calculates the functions and the jacobian.  fcn must
  !      be declared in an interface statement in the user
  !      calling program, and should be written as follows.

  !      subroutine fcn(m, n, x, fvec, fjac, iflag)
  !      integer   :: m, n, ldfjac, iflag
  !      REAL :: x(:), fvec(:), fjac(:,:)
  !      ----------
  !      if iflag = 1 calculate the functions at x and
  !      return this vector in fvec. do not alter fjac.
  !      if iflag = 2 calculate the jacobian at x and
  !      return this matrix in fjac. do not alter fvec.
  !      ----------
  !      return
  !      end

  !      the value of iflag should not be changed by fcn unless
  !      the user wants to terminate execution of lmder1.
  !      in this case set iflag to a negative integer.

  !    m is a positive integer input variable set to the number of functions.

  !    n is a positive integer input variable set to the number
  !      of variables.  n must not exceed m.

  !    x is an array of length n. on input x must contain
  !      an initial estimate of the solution vector. on output x
  !      contains the final estimate of the solution vector.

  !    fvec is an output array of length m which contains
  !      the functions evaluated at the output x.

  !    fjac is an output m by n array. the upper n by n submatrix
  !      of fjac contains an upper triangular matrix r with
  !      diagonal elements of nonincreasing magnitude such that

  !             t     t           t
  !            p *(jac *jac)*p = r *r,

  !      where p is a permutation matrix and jac is the final calculated
  !      Jacobian.  Column j of p is column ipvt(j) (see below) of the
  !      identity matrix.  The lower trapezoidal part of fjac contains
  !      information generated during the computation of r.

  !    ldfjac is a positive integer input variable not less than m
  !      which specifies the leading dimension of the array fjac.

  !    tol is a nonnegative input variable. termination occurs
  !      when the algorithm estimates either that the relative
  !      error in the sum of squares is at most tol or that
  !      the relative error between x and the solution is at most tol.

  !    info is an integer output variable.  If the user has terminated
  !      execution, info is set to the (negative) value of iflag.
  !      See description of fcn.  Otherwise, info is set as follows.

  !      info = 0  improper input parameters.

  !      info = 1  algorithm estimates that the relative error
  !                in the sum of squares is at most tol.

  !      info = 2  algorithm estimates that the relative error
  !                between x and the solution is at most tol.

  !      info = 3  conditions for info = 1 and info = 2 both hold.

  !      info = 4  fvec is orthogonal to the columns of the
  !                jacobian to machine precision.

  !      info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).

  !      info = 6  tol is too small.  No further reduction in
  !                the sum of squares is possible.

  !      info = 7  tol is too small.  No further improvement in
  !                the approximate solution x is possible.

  !    ipvt is an integer output array of length n. ipvt
  !      defines a permutation matrix p such that jac*p = q*r,
  !      where jac is the final calculated jacobian, q is
  !      orthogonal (not stored), and r is upper triangular
  !      with diagonal elements of nonincreasing magnitude.
  !      column j of p is column ipvt(j) of the identity matrix.

  !    wa is a work array of length lwa.

  !    lwa is a positive integer input variable not less than 5*n+m.

  !  subprograms called

  !    user-supplied ...... fcn

  !    minpack-supplied ... lmder

  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !
  ! N.B. Arguments LDFJAC, WA & LWA have been removed.
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)

    integer, intent(IN)                 :: m          !< is a positive integer input variable set to the number of functions.
    integer, intent(IN)                 :: n        !< is a positive integer input variable set to the number of variables.  n must not exceed m.
    real(dp), intent(IN OUT)            :: x(:)        !< is an array of length n.  On input x must contain an initial estimate of the solution vector.  On output x contains the final estimate of the  solution vector.
    real(dp), intent(OUT)               :: fvec(:)        !< s an output array of length m which contains  the functions evaluated at the output x.
    real(dp), intent(IN OUT)            :: fjac(:,:)  !< is an output m by n array. the upper n by n submatrix of fjac contains an upper triangular matrix r withdiagonal elements of nonincreasing magnitude such that \f$ p^t *(jac^t *jac)*p = r^t *r \f$ where p is a permutation matrix and jac is the final calculated Jacobian.  Column j of p is column ipvt(j) (see below) of the identity matrix. the lower trapezoidal part of fjac contains information generated during the computation of r.
    real(dp), intent(IN)                :: tol        !< is a nonnegative input variable.  Termination occurs when the algorithm estimates either that the relative error in the sum of squares is at most tol or that the relative error between x and the  solution is at most tol.
    integer, intent(OUT)                :: info  !< is an integer output variable.  If the user has terminated execution,      info is set to the (negative) value of iflag.  See description of fcn. Otherwise, info is set as follows. \n   info = 0  improper input parameters. \n   info = 1  algorithm estimates that the relative error in the sum of squares is at most tol. \n   info = 2  algorithm estimates that the relative error between x and the solution is at most tol. \n   info = 3  conditions for info = 1 and info = 2 both hold. \n   info = 4  fvec is orthogonal to the columns of the  jacobian to machine precision. \n   info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1). \n   info = 6  tol is too small. no further reduction in the sum of squares is possible. \n   info = 7  tol is too small.  No further improvement in the approximate solution x is possible.
    integer, intent(IN OUT)             :: ipvt(:)!< is an integer output array of length n.  ipvt defines a permutation matrix p such that jac*p = q*r, where jac is the final calculated jacobian, q is orthogonal (not stored), and r is upper triangular  with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.


    ! EXTERNAL fcn

    interface
       subroutine fcn(m, n, x, fvec, fjac, iflag)
         implicit none
         integer, parameter             :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
         integer, intent(IN)            :: m, n
         real(dp), intent(IN)           :: x(:)
         real(dp), intent(IN OUT)       :: fvec(:)
         real(dp), intent(OUT)          :: fjac(:,:)
         integer, intent(IN OUT)        :: iflag
       end subroutine fcn
    end interface

    integer                             :: maxfev, mode, nfev, njev, nprint
    real(dp)                            :: ftol, gtol, xtol
    real(dp), parameter                 :: factor = 100.0_dp, zero = 0.0_dp

    info = 0

    ! check the input parameters for errors.

    if ( n <= 0 .or. m < n .or. tol < zero ) GO TO 10

    ! call lmder.

    maxfev = 100*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    mode = 1
    nprint = 0
    call lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,  &
         mode, factor, nprint, info, nfev, njev, ipvt)
    if (info == 8) info = 4

10  return

    ! last card of subroutine lmder1.

  end subroutine lmder1



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine lmder
  !> \brief minimize the sum of the squares
  !!
  !!  the purpose of lmder is to minimize the sum of the squares of
  !!  m nonlinear functions in n variables by a modification of
  !!  the levenberg-marquardt algorithm. the user must provide a
  !!  subroutine which calculates the functions and the jacobian.
  !!
  !!  the subroutine statement is \n
  !!
  !!    subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
  !!                     maxfev,diag,mode,factor,nprint,info,nfev,
  !!                     njev,ipvt,qtf,wa1,wa2,wa3,wa4)\n
  !!  subprograms called\n
  !!    user-supplied ...... fcn \n
  !!    minpack-supplied ... dpmpar,enorm,lmpar,qrfac \n
  !!    fortran-supplied ... ABS,MAX,MIN,SQRT,mod \n
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !  where

  !    fcn is the name of the user-supplied subroutine which
  !      calculates the functions and the jacobian. fcn must
  !      be declared in an external statement in the user
  !      calling program, and should be written as follows.

  !      subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
  !      integer m,n,ldfjac,iflag
  !      REAL x(:),fvec(m),fjac(ldfjac,n)
  !      ----------
  !      if iflag = 1 calculate the functions at x and
  !      return this vector in fvec. do not alter fjac.
  !      if iflag = 2 calculate the jacobian at x and
  !      return this matrix in fjac.  Do not alter fvec.
  !      ----------
  !      return
  !      end

  !      the value of iflag should not be changed by fcn unless
  !      the user wants to terminate execution of lmder.
  !      in this case set iflag to a negative integer.

  !    m is a positive integer input variable set to the number
  !      of functions.

  !    n is a positive integer input variable set to the number
  !      of variables. n must not exceed m.

  !    x is an array of length n. on input x must contain
  !      an initial estimate of the solution vector. on output x
  !      contains the final estimate of the solution vector.

  !    fvec is an output array of length m which contains
  !      the functions evaluated at the output x.

  !    fjac is an output m by n array. the upper n by n submatrix
  !      of fjac contains an upper triangular matrix r with
  !      diagonal elements of nonincreasing magnitude such that

  !             t     t           t
  !            p *(jac *jac)*p = r *r

  !      where p is a permutation matrix and jac is the final calculated
  !      jacobian.  Column j of p is column ipvt(j) (see below) of the
  !      identity matrix.  The lower trapezoidal part of fjac contains
  !      information generated during the computation of r.

  !    ldfjac is a positive integer input variable not less than m
  !      which specifies the leading dimension of the array fjac.

  !    ftol is a nonnegative input variable.  Termination occurs when both
  !      the actual and predicted relative reductions in the sum of squares
  !      are at most ftol.   Therefore, ftol measures the relative error
  !      desired in the sum of squares.

  !    xtol is a nonnegative input variable. termination
  !      occurs when the relative error between two consecutive
  !      iterates is at most xtol. therefore, xtol measures the
  !      relative error desired in the approximate solution.

  !    gtol is a nonnegative input variable.  Termination occurs when the
  !      cosine of the angle between fvec and any column of the jacobian is
  !      at most gtol in absolute value.  Therefore, gtol measures the
  !      orthogonality desired between the function vector and the columns
  !      of the jacobian.

  !    maxfev is a positive integer input variable.  Termination occurs when
  !      the number of calls to fcn with iflag = 1 has reached maxfev.

  !    diag is an array of length n.  If mode = 1 (see below), diag is
  !      internally set.  If mode = 2, diag must contain positive entries
  !      that serve as multiplicative scale factors for the variables.

  !    mode is an integer input variable.  if mode = 1, the
  !      variables will be scaled internally.  if mode = 2,
  !      the scaling is specified by the input diag.  other
  !      values of mode are equivalent to mode = 1.

  !    factor is a positive input variable used in determining the
  !      initial step bound. this bound is set to the product of
  !      factor and the euclidean norm of diag*x if nonzero, or else
  !      to factor itself. in most cases factor should lie in the
  !      interval (.1,100.).100. is a generally recommended value.

  !    nprint is an integer input variable that enables controlled printing
  !      of iterates if it is positive.  In this case, fcn is called with
  !      iflag = 0 at the beginning of the first iteration and every nprint
  !      iterations thereafter and immediately prior to return, with x, fvec,
  !      and fjac available for printing.  fvec and fjac should not be
  !      altered.  If nprint is not positive, no special calls of fcn with
  !      iflag = 0 are made.

  !    info is an integer output variable.  If the user has terminated
  !      execution, info is set to the (negative) value of iflag.
  !      See description of fcn.  Otherwise, info is set as follows.

  !      info = 0  improper input parameters.

  !      info = 1  both actual and predicted relative reductions
  !                in the sum of squares are at most ftol.

  !      info = 2  relative error between two consecutive iterates
  !                is at most xtol.

  !      info = 3  conditions for info = 1 and info = 2 both hold.

  !      info = 4  the cosine of the angle between fvec and any column of
  !                the jacobian is at most gtol in absolute value.

  !      info = 5  number of calls to fcn with iflag = 1 has reached maxfev.

  !      info = 6  ftol is too small.  No further reduction in
  !                the sum of squares is possible.

  !      info = 7  xtol is too small.  No further improvement in
  !                the approximate solution x is possible.

  !      info = 8  gtol is too small.  fvec is orthogonal to the
  !                columns of the jacobian to machine precision.

  !    nfev is an integer output variable set to the number of
  !      calls to fcn with iflag = 1.

  !    njev is an integer output variable set to the number of
  !      calls to fcn with iflag = 2.

  !    ipvt is an integer output array of length n.  ipvt
  !      defines a permutation matrix p such that jac*p = q*r,
  !      where jac is the final calculated jacobian, q is
  !      orthogonal (not stored), and r is upper triangular
  !      with diagonal elements of nonincreasing magnitude.
  !      column j of p is column ipvt(j) of the identity matrix.

  !    qtf is an output array of length n which contains
  !      the first n elements of the vector (q transpose)*fvec.

  !    wa1, wa2, and wa3 are work arrays of length n.

  !    wa4 is a work array of length m.

  !  subprograms called

  !    user-supplied ...... fcn

  !    minpack-supplied ... dpmpar,enorm,lmpar,qrfac

  !    fortran-supplied ... ABS,MAX,MIN,SQRT,mod

  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more

  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev, &
       mode, factor, nprint, info, nfev, njev, ipvt)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:50

    ! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.

    integer, intent(IN)                 :: m        !< is a positive integer input variable set to the number of functions.
    integer, intent(IN)                 :: n        !< is a positive integer input variable set to the number of variables.  n must not exceed m.
    real(dp), intent(IN OUT)            :: x(:)        !< is an array of length n.  On input x must contain an initial estimate of the solution vector.  On output x contains the final estimate of the  solution vector.
    real(dp), intent(OUT)               :: fvec(m)        !< is an output array of length m which contains  the functions evaluated at the output x.
    real(dp), intent(OUT)               :: fjac(:,:)   !< is an output m by n array. the upper n by n submatrix of fjac contains an upper triangular matrix r withdiagonal elements of nonincreasing magnitude such that \f$ p^t *(jac^t *jac)*p = r^t *r \f$ where p is a permutation matrix and jac is the final calculated Jacobian.  Column j of p is column ipvt(j) (see below) of the identity matrix. the lower trapezoidal part of fjac contains information generated during the computation of r.
    real(dp), intent(IN)                :: ftol        !< is a nonnegative input variable.  Termination occurs when both the actual and predicted relative reductions in the sum of squares are at most ftol.  Therefore, ftol measures the relative error desired in the sum of squares.
    real(dp), intent(IN)                :: xtol        !< is a nonnegative input variable.  Termination occurs when the relative error between two consecutive iterates is at most xtol. Therefore, xtol measures the relative error desired in the approximate solution.
    real(dp), intent(IN OUT)            :: gtol        !< is a nonnegative input variable.  Termination occurs when the cosine  of the angle between fvec and any column of the jacobian is at most  gtol in absolute value.  Therefore, gtol measures the orthogonality  desired between the function vector and the columns of the jacobian.
    integer, intent(IN OUT)             :: maxfev!< is a positive integer input variable.  Termination occurs when the  number of calls to fcn is at least maxfev by the end of an iteration.
    integer, intent(IN)                 :: mode        !< is an integer input variable.  If mode = 1, the variables will be scaled internally.  If mode = 2, the scaling is specified by the input diag. other values of mode are equivalent to mode = 1.
    real(dp), intent(IN)                :: factor        !< is a positive input variable used in determining the initial step bound.  This bound is set to the product of factor and the euclidean norm of diag*x if nonzero, or else to factor itself.  In most cases  factor should lie in the interval (.1,100.). 100. is a generally recommended value.
    integer, intent(IN)                 :: nprint!< is an integer input variable that enables controlled printing of iterates if it is positive.  In this case, fcn is called with iflag = 0 at the beginning of the first iteration and every nprint iterations  thereafter and immediately prior to return, with x and fvec available for printing.  If nprint is not positive, no special calls of fcn with iflag = 0 are made.
    integer, intent(OUT)                :: info        !< is an integer output variable.  If the user has terminated  execution, info is set to the (negative) value of iflag. See description of fcn.  Otherwise, info is set as follows. \n info = 0  improper input parameters. \n info = 1  both actual and predicted relative reductions in the sum of squares are at most ftol. \n info = 2  relative error between two consecutive iterates <= xtol.    \n info = 3  conditions for info = 1 and info = 2 both hold. \n info = 4  the cosine of the angle between fvec and any column of  the Jacobian is at most gtol in absolute value. \n info = 5  number of calls to fcn has reached or exceeded maxfev. \n   info = 6  ftol is too small. no further reduction in  the sum of squares is possible. \n    info = 7  xtol is too small. no further improvement in the approximate solution x is possible. \n info = 8  gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.
    integer, intent(OUT)                :: nfev  !< is an integer output variable set to the number of calls to fcn with iflag = 1.
    integer, intent(OUT)                :: njev         !< is an integer output variable set to the number of calls to fcn with iflag = 2.
    integer, intent(OUT)                :: ipvt(:)!< is an integer output array of length n.  ipvt defines a permutation matrix p such that jac*p = q*r, where jac is the final calculated jacobian, q is orthogonal (not stored), and r is upper triangular  with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.

    interface
       subroutine fcn(m, n, x, fvec, fjac, iflag)
         implicit none
         integer, parameter             :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
         integer, intent(IN)            :: m, n
         real(dp), intent(IN)           :: x(:)
         real(dp), intent(IN OUT)       :: fvec(:)
         real(dp), intent(OUT)          :: fjac(:,:)
         integer, intent(IN OUT)        :: iflag
       end subroutine fcn
    end interface


    integer                             :: i, iflag, iter, j, l
    real(dp)                            :: actred, delta, dirder, epsmch
    real(dp)                            :: fnorm, fnorm1, gnorm
    real(dp)                            :: par, pnorm, prered, ratio
    real(dp)                            :: sum, temp, temp1, temp2, xnorm
    real(dp)                            :: diag(n), qtf(n)
    real(dp)                            :: wa1(n), wa2(n), wa3(n), wa4(m)
    real(dp), parameter                 :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp
    real(dp), parameter                 :: p25 = 0.25_dp, p75 = 0.75_dp
    real(dp), parameter                 :: p0001 = 0.0001_dp, zero = 0.0_dp

    ! epsmch is the machine precision.

    epsmch = epsilon(zero)

    info = 0
    iflag = 0
    nfev = 0
    njev = 0

    ! check the input parameters for errors.

    if (n <= 0 .or. m < n .or. ftol < zero .or. xtol < zero .or. gtol < zero  &
         .or. maxfev <= 0 .or. factor <= zero) GO TO 300
    if (mode /= 2) GO TO 20
    do  j = 1, n
       if (diag(j) <= zero) GO TO 300
    end do

    ! evaluate the function at the starting point and calculate its norm.

20  iflag = 1
    call fcn(m, n, x, fvec, fjac, iflag)
    nfev = 1
    if (iflag < 0) GO TO 300
    fnorm = enorm(m, fvec)

    ! initialize levenberg-marquardt parameter and iteration counter.

    par = zero
    iter = 1

    ! beginning of the outer loop.

    ! calculate the jacobian matrix.

30  iflag = 2
    call fcn(m, n, x, fvec, fjac, iflag)
    njev = njev + 1
    if (iflag < 0) GO TO 300

    ! if requested, call fcn to enable printing of iterates.

    if (nprint <= 0) GO TO 40
    iflag = 0
    if (mod(iter-1,nprint) == 0) call fcn(m, n, x, fvec, fjac, iflag)
    if (iflag < 0) GO TO 300

    ! compute the qr factorization of the jacobian.

40  call qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

    ! on the first iteration and if mode is 1, scale according
    ! to the norms of the columns of the initial jacobian.

    if (iter /= 1) GO TO 80
    if (mode == 2) GO TO 60
    do  j = 1, n
       diag(j) = wa2(j)
       if ( abs( wa2(j) ) < machine_eps ) diag(j) = one
    end do

    ! on the first iteration, calculate the norm of the scaled x
    ! and initialize the step bound delta.

60  wa3(1:n) = diag(1:n)*x(1:n)
    xnorm = enorm(n,wa3)
    delta = factor*xnorm
    if ( abs( delta ) < machine_eps ) delta = factor

    ! form (q transpose)*fvec and store the first n components in qtf.

80  wa4(1:m) = fvec(1:m)
    do  j = 1, n
       if ( abs( fjac(j,j) ) < machine_eps ) GO TO 120
       sum = DOT_product( fjac(j:m,j), wa4(j:m) )
       temp = -sum/fjac(j,j)
       do  i = j, m
          wa4(i) = wa4(i) + fjac(i,j)*temp
       end do
120    fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
    end do

    ! compute the norm of the scaled gradient.

    gnorm = zero
    if ( abs( fnorm ) < machine_eps ) GO TO 170
    do  j = 1, n
       l = ipvt(j)
       if ( abs( wa2(l) ) < machine_eps ) cycle
       sum = zero
       do  i = 1, j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
       end do
       gnorm = max(gnorm,abs(sum/wa2(l)))
    end do

    ! test for convergence of the gradient norm.

170 if (gnorm <= gtol) info = 4
    if (info /= 0) GO TO 300

    ! rescale if necessary.

    if (mode == 2) GO TO 200
    do  j = 1, n
       diag(j) = max(diag(j), wa2(j))
    end do

    ! beginning of the inner loop.

    !   determine the levenberg-marquardt parameter.

200 call lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

    !   store the direction p and x + p. calculate the norm of p.

    do  j = 1, n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    end do
    pnorm = enorm(n, wa3)

    !   on the first iteration, adjust the initial step bound.

    if (iter == 1) delta = min(delta,pnorm)

    !   evaluate the function at x + p and calculate its norm.

    iflag = 1
    call fcn(m, n, wa2, wa4, fjac, iflag)
    nfev = nfev + 1
    if (iflag < 0) GO TO 300
    fnorm1 = enorm(m, wa4)

    !   compute the scaled actual reduction.

    actred = -one
    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

    !   compute the scaled predicted reduction and
    !   the scaled directional derivative.

    do  j = 1, n
       wa3(j) = zero
       l = ipvt(j)
       temp = wa1(l)
       wa3(1:j) = wa3(1:j) + fjac(1:j,j)*temp
    end do
    temp1 = enorm(n,wa3)/fnorm
    temp2 = (sqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)

    !   compute the ratio of the actual to the predicted reduction.

    ratio = zero
    if ( abs( prered ) > machine_eps ) ratio = actred/prered

    !   update the step bound.

    if (ratio <= p25) then
       if (actred >= zero) temp = p5
       if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
       if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
       delta = temp*min(delta, pnorm/p1)
       par = par/temp
    else
       if ( abs( par ) > machine_eps .and. ratio < p75) GO TO 260
       delta = pnorm/p5
       par = p5*par
    end if

    !   test for successful iteration.

260 if (ratio < p0001) GO TO 290

    !   successful iteration. update x, fvec, and their norms.

    do  j = 1, n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
    end do
    fvec(1:m) = wa4(1:m)
    xnorm = enorm(n,wa2)
    fnorm = fnorm1
    iter = iter + 1

    !   tests for convergence.

290 if (abs(actred) <= ftol .and. prered <= ftol .and. p5*ratio <= one) info = 1
    if (delta <= xtol*xnorm) info = 2
    if (abs(actred) <= ftol .and. prered <= ftol  &
         .and. p5*ratio <= one .and. info == 2) info = 3
    if (info /= 0) GO TO 300

    !   tests for termination and stringent tolerances.

    if (nfev >= maxfev) info = 5
    if (abs(actred) <= epsmch .and. prered <= epsmch  &
         .and. p5*ratio <= one) info = 6
    if (delta <= epsmch*xnorm) info = 7
    if (gnorm <= epsmch) info = 8
    if (info /= 0) GO TO 300

    !   end of the inner loop. repeat if iteration unsuccessful.

    if (ratio < p0001) GO TO 200

    ! end of the outer loop.

    GO TO 30

    ! termination, either normal or user imposed.

300 if (iflag < 0) info = iflag
    iflag = 0
    if (nprint > 0) call fcn(m, n, x, fvec, fjac, iflag)

  end subroutine lmder


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine lmpar
  !> \brief determine a value for the parameter par to solve the system
  !!
  !! Given an m by n matrix a, an n by n nonsingular diagonal matrix d,
  !!  an m-vector b, and a positive number delta, the problem is to determine a
  !!  value for the parameter par such that if x solves the system \n
  !!     \f$    a*x = b \f$,    \f$ sqrt(par)*d*x = 0 \f$, \n
  !!  in the least squares sense, and dxnorm is the euclidean
  !!  norm of d*x, then either par is zero and
  !!        (dxnorm-delta) <= 0.1*delta ,
  !!  or par is positive and
  !!        abs(dxnorm-delta) <= 0.1*delta .
  !!  This subroutine completes the solution of the problem if it is provided
  !!  with the necessary information from the r factorization, with column
  !!  qpivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
  !!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
  !!  elements of nonincreasing magnitude, then lmpar expects the full upper
  !!  triangle of r, the permutation matrix p, and the first n components of
  !!  (q transpose)*b.
  !!  On output lmpar also provides an upper triangular matrix s such that \n
  !!   \f$     p^t *(a^t *a + par*d*d)*p = s^t *s .\f$ \n
  !!  s is employed within lmpar and may be of separate interest.
  !!  Only a few iterations are generally needed for convergence of the algorithm.
  !!  If, however, the limit of 10 iterations is reached, then the output par
  !!  will contain the best value obtained so far. \n
  !!  the subroutine statement is \n
  !!    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2) \n
  !!  subprograms called \n
  !!     minpack-supplied ... dpmpar,enorm,qrsolv \n
  !!    fortran-supplied ... ABS,MAX,MIN,SQRT \n
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more \n
  !! Code converted using TO_F90 by Alan Miller
  !! Date: 1999-12-09  Time: 12:46:12
  !! N.B. Arguments LDR, WA1 & WA2 have been removed.
  !  where
  !    n is a positive integer input variable set to the order of r.
  !    r is an n by n array. on input the full upper triangle
  !      must contain the full upper triangle of the matrix r.
  !      On output the full upper triangle is unaltered, and the
  !      strict lower triangle contains the strict upper triangle
  !      (transposed) of the upper triangular matrix s.
  !    ldr is a positive integer input variable not less than n
  !      which specifies the leading dimension of the array r.
  !    ipvt is an integer input array of length n which defines the
  !      permutation matrix p such that a*p = q*r. column j of p
  !      is column ipvt(j) of the identity matrix.
  !    diag is an input array of length n which must contain the
  !      diagonal elements of the matrix d.
  !    qtb is an input array of length n which must contain the first
  !      n elements of the vector (q transpose)*b.
  !    delta is a positive input variable which specifies an upper
  !      bound on the euclidean norm of d*x.
  !    par is a nonnegative variable. on input par contains an
  !      initial estimate of the levenberg-marquardt parameter.
  !      on output par contains the final estimate.
  !    x is an output array of length n which contains the least
  !      squares solution of the system a*x = b, sqrt(par)*d*x = 0,
  !      for the output par.
  !    sdiag is an output array of length n which contains the
  !      diagonal elements of the upper triangular matrix s.
  !    wa1 and wa2 are work arrays of length n.
  !  subprograms called
  !    minpack-supplied ... dpmpar,enorm,qrsolv
  !    fortran-supplied ... ABS,MAX,MIN,SQRT
  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more \n
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)

    integer, intent(IN)                 :: n
    real(dp), intent(IN OUT)            :: r(:,:)
    integer, intent(IN)                 :: ipvt(:)
    real(dp), intent(IN)                :: diag(:)
    real(dp), intent(IN)                :: qtb(:)
    real(dp), intent(IN)                :: delta
    real(dp), intent(OUT)               :: par
    real(dp), intent(OUT)               :: x(:)
    real(dp), intent(OUT)               :: sdiag(:)

    integer                             :: iter, j, jm1, jp1, k, l, nsing
    real(dp)                            :: dxnorm, dwarf, fp, gnorm
    real(dp)                            :: parc, parl, paru, sum, temp
    real(dp)                            :: wa1(n), wa2(n)
    real(dp)                            :: p001_times_paru
    real(dp), parameter                 :: p1 = 0.1_dp, p001 = 0.001_dp, zero = 0.0_dp

    ! dwarf is the smallest positive magnitude.

    dwarf = tiny(zero)

    ! compute and store in x the gauss-newton direction. if the
    ! jacobian is rank-deficient, obtain a least squares solution.

    nsing = n
    do  j = 1, n
       wa1(j) = qtb(j)
       if ( abs( r(j,j) ) < machine_eps .and. nsing == n) nsing = j - 1
       if (nsing < n) wa1(j) = zero
    end do

    do  k = 1, nsing
       j = nsing - k + 1
       wa1(j) = wa1(j)/r(j,j)
       temp = wa1(j)
       jm1 = j - 1
       wa1(1:jm1) = wa1(1:jm1) - r(1:jm1,j)*temp
    end do

    do  j = 1, n
       l = ipvt(j)
       x(l) = wa1(j)
    end do

    ! initialize the iteration counter.
    ! evaluate the function at the origin, and test
    ! for acceptance of the gauss-newton direction.

    iter = 0
    wa2(1:n) = diag(1:n)*x(1:n)
    dxnorm = enorm(n, wa2)
    fp = dxnorm - delta
    if (fp <= p1*delta) GO TO 220

    ! if the jacobian is not rank deficient, the newton
    ! step provides a lower bound, parl, for the zero of
    ! the function.  Otherwise set this bound to zero.

    parl = zero
    if (nsing < n) GO TO 120
    do  j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    end do
    do  j = 1, n
       sum = DOT_product( r(1:j-1,j), wa1(1:j-1) )
       wa1(j) = (wa1(j) - sum)/r(j,j)
    end do
    temp = enorm(n,wa1)
    parl = ((fp/delta)/temp)/temp

    ! calculate an upper bound, paru, for the zero of the function.

120 do  j = 1, n
       sum = DOT_product( r(1:j,j), qtb(1:j) )
       l = ipvt(j)
       wa1(j) = sum/diag(l)
    end do
    gnorm = enorm(n,wa1)
    paru = gnorm/delta
    if ( abs( paru ) < machine_eps ) paru = dwarf/min(delta,p1)

    ! if the input par lies outside of the interval (parl,paru),
    ! set par to the closer endpoint.

    par = max(par,parl)
    par = min(par,paru)
    if ( abs( par ) < machine_eps ) par = gnorm/dxnorm

    ! beginning of an iteration.

150 iter = iter + 1

    ! evaluate the function at the current value of par.

    if ( LOG10(p001)+LOG10(paru) < -307.0_dp ) then
       p001_times_paru = 0.0_dp
    else
       p001_times_paru = p001 * paru
    end if
    if ( abs( par ) < machine_eps ) par = max( dwarf, p001_times_paru )
    temp = sqrt(par)
    wa1(1:n) = temp*diag(1:n)
    call qrsolv(n, r, ipvt, wa1, qtb, x, sdiag)
    wa2(1:n) = diag(1:n)*x(1:n)
    dxnorm = enorm(n, wa2)
    temp = fp
    fp = dxnorm - delta

    ! if the function is small enough, accept the current value
    ! of par. also test for the exceptional cases where parl
    ! is zero or the number of iterations has reached 10.

    if (abs(fp) <= p1*delta .or. abs( parl ) < machine_eps  .and. fp <= temp  &
         .and. temp < zero .or. iter == 10) GO TO 220

    ! compute the newton correction.

    do  j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    end do
    do  j = 1, n
       wa1(j) = wa1(j)/sdiag(j)
       temp = wa1(j)
       jp1 = j + 1
       wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
    end do
    temp = enorm(n,wa1)
    parc = ((fp/delta)/temp)/temp

    ! depending on the sign of the function, update parl or paru.

    if (fp > zero) parl = max(parl,par)
    if (fp < zero) paru = min(paru,par)

    ! compute an improved estimate for par.

    par = max(parl, par+parc)

    ! end of an iteration.

    GO TO 150

    ! termination.

220 if (iter == 0) par = zero

  end subroutine lmpar


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    !  subroutine qrfac
    !>  computes a qr factorization
    !!
    !!  This subroutine uses Householder transformations with column pivoting
    !!  (optional) to compute a qr factorization of the m by n matrix a.
    !!  That is, qrfac determines an orthogonal matrix q, a permutation matrix p,
    !!  and an upper trapezoidal matrix r with diagonal elements of nonincreasing
    !!  magnitude, such that a*p = q*r.  The householder transformation for
    !!   column k, k = 1,2,...,min(m,n), is of the form \n \f$  i - (1/u(k))*u*u^t \f$ \n
    !!  where u has zeros in the first k-1 positions.  The form of this
    !!  transformation and the method of pivoting first appeared in the
    !!  corresponding linpack subroutine.\n
    !!  the subroutine statement is \n
    !!    subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa) \n
    !! N.B. 3 of these arguments have been omitted in this version. \n
    !!  subprograms called \n
    !!    minpack-supplied ... dpmpar,enorm \n
    !!    fortran-supplied ... MAX,SQRT,MIN \n
    !!  argonne national laboratory. minpack project. march 1980.
    !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
    !  where

    !    m is a positive integer input variable set to the number of rows of a.

    !    n is a positive integer input variable set to the number of columns of a.

    !    a is an m by n array.  On input a contains the matrix for
    !      which the qr factorization is to be computed.  On output
    !      the strict upper trapezoidal part of a contains the strict
    !      upper trapezoidal part of r, and the lower trapezoidal
    !      part of a contains a factored form of q (the non-trivial
    !      elements of the u vectors described above).

    !    lda is a positive integer input variable not less than m
    !      which specifies the leading dimension of the array a.

    !    pivot is a logical input variable.  If pivot is set true,
    !      then column pivoting is enforced.  If pivot is set false,
    !      then no column pivoting is done.

    !    ipvt is an integer output array of length lipvt.  ipvt
    !      defines the permutation matrix p such that a*p = q*r.
    !      Column j of p is column ipvt(j) of the identity matrix.
    !      If pivot is false, ipvt is not referenced.

    !    lipvt is a positive integer input variable.  If pivot is false,
    !      then lipvt may be as small as 1.  If pivot is true, then
    !      lipvt must be at least n.

    !    rdiag is an output array of length n which contains the
    !      diagonal elements of r.

    !    acnorm is an output array of length n which contains the norms of the
    !      corresponding columns of the input matrix a.
    !      If this information is not needed, then acnorm can coincide with rdiag.

    !    wa is a work array of length n.  If pivot is false, then wa
    !      can coincide with rdiag.

    !  subprograms called

    !    minpack-supplied ... dpmpar,enorm

    !    fortran-supplied ... MAX,SQRT,MIN

    !  argonne national laboratory. minpack project. march 1980.
    !  burton s. garbow, kenneth e. hillstrom, jorge j. more

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)

    integer, intent(IN)                 :: m
    integer, intent(IN)                 :: n
    real(dp), intent(IN OUT)            :: a(:,:)
    logical, intent(IN)                 :: pivot
    integer, intent(OUT)                :: ipvt(:)
    real(dp), intent(OUT)               :: rdiag(:)
    real(dp), intent(OUT)               :: acnorm(:)

    integer                             :: i, j, jp1, k, kmax, minmn
    real(dp)                            :: ajnorm, epsmch, sum, temp, wa(n)
    real(dp), parameter                 :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

    ! epsmch is the machine precision.

    epsmch = epsilon(zero)

    ! compute the initial column norms and initialize several arrays.

    do  j = 1, n
       acnorm(j) = enorm(m,a(1:,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       if (pivot) ipvt(j) = j
    end do

    ! Reduce a to r with Householder transformations.

    minmn = min(m,n)
    do  j = 1, minmn
       if (.not.pivot) GO TO 40

       ! Bring the column of largest norm into the pivot position.

       kmax = j
       do  k = j, n
          if (rdiag(k) > rdiag(kmax)) kmax = k
       end do
       if (kmax == j) GO TO 40
       do  i = 1, m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
       end do
       rdiag(kmax) = rdiag(j)
       wa(kmax) = wa(j)
       k = ipvt(j)
       ipvt(j) = ipvt(kmax)
       ipvt(kmax) = k

       ! Compute the Householder transformation to reduce the
       ! j-th column of a to a multiple of the j-th unit vector.

40     ajnorm = enorm(m-j+1, a(j:,j))
       if ( abs( ajnorm ) < machine_eps ) cycle
       if (a(j,j) < zero) ajnorm = -ajnorm
       a(j:m,j) = a(j:m,j)/ajnorm
       a(j,j) = a(j,j) + one

       ! Apply the transformation to the remaining columns and update the norms.

       jp1 = j + 1
       do  k = jp1, n
          sum = DOT_product( a(j:m,j), a(j:m,k) )
          temp = sum/a(j,j)
          a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
          if (.not.pivot .or. abs( rdiag(k) ) < machine_eps ) cycle
          temp = a(j,k)/rdiag(k)
          rdiag(k) = rdiag(k)*sqrt(max(zero, one-temp**2))
          if (p05*(rdiag(k)/wa(k))**2 > epsmch) cycle
          rdiag(k) = enorm(m-j, a(jp1:,k))
          wa(k) = rdiag(k)
       end do
       rdiag(j) = -ajnorm
    end do

  end subroutine qrfac


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  subroutine qrsolv
  !
  !!  Given an m x n matrix a, an n x n diagonal matrix d, and an m-vector b,
  !!  the problem is to determine an x which solves the system
  !!        a*x = b ,     d*x = 0 ,
  !!  in the least squares sense.
  !!
  !!  This subroutine completes the solution of the problem if it is provided
  !!  with the necessary information from the qr factorization, with column
  !!  pivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
  !!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
  !!  elements of nonincreasing magnitude, then qrsolv expects the full upper
  !!  triangle of r, the permutation matrix p, and the first n components of
  !!  (q transpose)*b.  The system a*x = b, d*x = 0, is then equivalent to
  !!               t       t
  !!        r*z = q *b ,  p *d*p*z = 0 ,
  !!  where x = p*z. if this system does not have full rank,
  !!  then a least squares solution is obtained.  On output qrsolv
  !!  also provides an upper triangular matrix s such that
  !!         t   t               t
  !!        p *(a *a + d*d)*p = s *s .
  !!  s is computed within qrsolv and may be of separate interest.
  !!  the subroutine statement is
  !!    subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
  !! N.B. Arguments LDR and WA have been removed in this version.
  !!  where
  !!    n is a positive integer input variable set to the order of r.
  !!    r is an n by n array.  On input the full upper triangle must contain
  !!      the full upper triangle of the matrix r.
  !!      On output the full upper triangle is unaltered, and the strict lower
  !!      triangle contains the strict upper triangle (transposed) of the
  !!      upper triangular matrix s.
  !!    ldr is a positive integer input variable not less than n
  !!      which specifies the leading dimension of the array r.
  !!    ipvt is an integer input array of length n which defines the
  !!      permutation matrix p such that a*p = q*r.  Column j of p
  !!      is column ipvt(j) of the identity matrix.
  !!    diag is an input array of length n which must contain the
  !!      diagonal elements of the matrix d.
  !!    qtb is an input array of length n which must contain the first
  !!      n elements of the vector (q transpose)*b.
  !!    x is an output array of length n which contains the least
  !!      squares solution of the system a*x = b, d*x = 0.
  !!    sdiag is an output array of length n which contains the
  !!      diagonal elements of the upper triangular matrix s.
  !!    wa is a work array of length n.
  !!
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine qrsolv(n, r, ipvt, diag, qtb, x, sdiag)

    integer, intent(IN)                 :: n
    real(dp), intent(IN OUT)            :: r(:,:)
    integer, intent(IN)                 :: ipvt(:)
    real(dp), intent(IN)                :: diag(:)
    real(dp), intent(IN)                :: qtb(:)
    real(dp), intent(OUT)               :: x(:)
    real(dp), intent(OUT)               :: sdiag(:)

    integer                             :: i, j, k, kp1, l, nsing
    real(dp)                            :: COS, cotan, qtbpj, SIN, sum, TAN, temp, wa(n)
    real(dp), parameter                 :: p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

    ! Copy r and (q transpose)*b to preserve input and initialize s.
    ! In particular, save the diagonal elements of r in x.

    do  j = 1, n
       r(j:n,j) = r(j,j:n)
       x(j) = r(j,j)
       wa(j) = qtb(j)
    end do

    ! Eliminate the diagonal matrix d using a givens rotation.

    do  j = 1, n

       ! Prepare the row of d to be eliminated, locating the
       ! diagonal element using p from the qr factorization.

       l = ipvt(j)
       if ( abs( diag(l) ) < machine_eps ) cycle
       sdiag(j:n) = zero
       sdiag(j) = diag(l)

       ! The transformations to eliminate the row of d modify only a single
       ! element of (q transpose)*b beyond the first n, which is initially zero.

       qtbpj = zero
       do  k = j, n

          ! Determine a givens rotation which eliminates the
          ! appropriate element in the current row of d.

          if ( abs( sdiag(k) ) < machine_eps ) cycle
          if (abs(r(k,k)) < abs(sdiag(k))) then
             cotan = r(k,k)/sdiag(k)
             SIN = p5/sqrt(p25 + p25*cotan**2)
             COS = SIN*cotan
          else
             TAN = sdiag(k)/r(k,k)
             COS = p5/sqrt(p25 + p25*TAN**2)
             SIN = COS*TAN
          end if

          ! Compute the modified diagonal element of r and
          ! the modified element of ((q transpose)*b,0).

          r(k,k) = COS*r(k,k) + SIN*sdiag(k)
          temp = COS*wa(k) + SIN*qtbpj
          qtbpj = -SIN*wa(k) + COS*qtbpj
          wa(k) = temp

          ! Accumulate the tranformation in the row of s.

          kp1 = k + 1
          do  i = kp1, n
             temp = COS*r(i,k) + SIN*sdiag(i)
             sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
             r(i,k) = temp
          end do
       end do

       ! Store the diagonal element of s and restore
       ! the corresponding diagonal element of r.

       sdiag(j) = r(j,j)
       r(j,j) = x(j)
    end do

    ! Solve the triangular system for z.  If the system is singular,
    ! then obtain a least squares solution.

    nsing = n
    do  j = 1, n
       if ( abs( sdiag(j) ) < machine_eps .and. nsing == n) nsing = j - 1
       if (nsing < n) wa(j) = zero
    end do

    do  k = 1, nsing
       j = nsing - k + 1
       sum = DOT_product( r(j+1:nsing,j), wa(j+1:nsing) )
       wa(j) = (wa(j) - sum)/sdiag(j)
    end do

    ! Permute the components of z back to components of x.

    do  j = 1, n
       l = ipvt(j)
       x(l) = wa(j)
    end do

  end subroutine qrsolv



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  function enorm
  !> \brief Calculation of the euclidean norm
  !!
  !!  given an n-vector x, this function calculates the euclidean norm of x.
  !!  the euclidean norm is computed by accumulating the sum of squares in
  !!  three different sums.  The sums of squares for the small and large
  !!  components are scaled so that no overflows occur.  Non-destructive
  !!  underflows are permitted.  Underflows and overflows do not occur in the
  !!  computation of the unscaled sum of squares for the intermediate
  !!  components.  The definitions of small, intermediate and large components
  !!  depend on two constants, rdwarf and rgiant.  The main restrictions on
  !!  these constants are that rdwarf**2 not underflow and rgiant**2 not
  !!  overflow.  The constants given here are suitable for every known computer. \n
  !!  the function statement is \n
  !!    REAL (dp) function enorm(n,x) \n
  !!  subprograms called \n
  !!    fortran-supplied ... ABS,SQRT\n
  !!  argonne national laboratory. minpack project. march 1980.
  !!  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !  where
  !    n is a positive integer input variable.
  !    x is an input array of length n.
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  function enorm(n,x) result(fn_val)

    integer, intent(IN)                 :: n           !< is a positive integer input variable.
    real(dp), intent(IN)                :: x(:)        !< is an input array of length n.
    real(dp)                            :: fn_val

    integer                             :: i
    real(dp)                            :: agiant, floatn
    real(dp)                            :: s1, s2, s3
    real(dp)                            :: xabs, x1max, x3max
    real(dp), parameter                 :: one = 1.0_dp, zero = 0.0_dp
    real(dp), parameter                 :: rdwarf = 3.834E-20_dp, rgiant = 1.304E+19_dp

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = n
    agiant = rgiant/floatn
    do  i = 1, n
       xabs = abs(x(i))
       if (xabs > rdwarf .and. xabs < agiant) GO TO 70
       if (xabs <= rdwarf) GO TO 30

       ! sum for large components.

       if (xabs <= x1max) GO TO 10
       s1 = one + s1*(x1max/xabs)**2
       x1max = xabs
       GO TO 20

10     s1 = s1 + (xabs/x1max)**2

20     GO TO 60

       ! sum for small components.

30     if (xabs <= x3max) GO TO 40
       s3 = one + s3*(x3max/xabs)**2
       x3max = xabs
       GO TO 60

40     if ( abs( xabs ) > machine_eps ) s3 = s3 + (xabs/x3max)**2

60     cycle

       ! sum for intermediate components.

70     s2 = s2 + xabs**2
    end do

    ! calculation of norm.

    if ( abs( s1 ) > machine_eps ) then
       fn_val = x1max*sqrt(s1 + (s2/x1max)/x1max)
    else

       if ( abs( s2 ) < machine_eps ) then
          fn_val = x3max*sqrt(s3)
       else
          if (s2 >= x3max) then
             fn_val = sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
          else
             fn_val = sqrt(x3max*((s2/x3max) + (x3max*s3)))
          end if
       end if
    end if

  end function enorm


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! subroutine fdjac2
  !
  ! this subroutine computes a forward-difference approximation
  ! to the m by n jacobian matrix associated with a specified
  ! problem of m functions in n variables.
  ! fcn   user-supplied subroutine which calculates the functions.
  !       fcn must be declared in an external statement in the user
  !       calling program, and should be written as follows.
  !           subroutine fcn(m,n,x,fvec,iflag)
  !       the value of iflag should not be changed by fcn unless
  !       the user wants to terminate execution of fdjac2.
  !       in this case set iflag to a negative integer.
  ! m     positive integer input variable set to the number of functions.
  ! n     is a positive integer input variable set to the number of variables.
  !       n must not exceed m.
  ! x     input array of length n.
  ! fvec  input array of length m which must contain the functions evaluated at x.
  ! fjac  output m by n array which contains the
  !       approximation to the jacobian matrix evaluated at x.
  ! ldfjac positive integer input variable not less than m
  !       which specifies the leading dimension of the array fjac.
  ! iflag integer variable which can be used to terminate
  !       the execution of fdjac2.  see description of fcn.
  ! epsfcn input variable used in determining a suitable step length
  !       for the forward-difference approximation.  This approximation assumes
  !       that the relative errors in the functions are of the order of epsfcn.
  !       If epsfcn is less than the machine precision, it is assumed that the
  !       relative errors in the functions are of the order of the machine
  !       precision.
  ! wa    work array of length m.
  !
  !  argonne national laboratory. minpack project. march 1980.
  !  burton s. garbow, kenneth e. hillstrom, jorge j. more
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)

    integer, intent(IN)                 :: m
    integer, intent(IN)                 :: n
    real(dp), intent(IN OUT)            :: x(n)
    real(dp), intent(IN)                :: fvec(m)
    real(dp), intent(OUT)               :: fjac(:,:)    ! fjac(ldfjac,n)
    integer, intent(IN OUT)             :: iflag
    real(dp), intent(IN)                :: epsfcn

    interface
       subroutine fcn(m, n, x, fvec, iflag)
         implicit none
         integer, parameter             :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
         integer, intent(IN)            :: m, n
         real(dp), intent(IN)           :: x(:)
         real(dp), intent(IN OUT)       :: fvec(:)
         integer, intent(IN OUT)        :: iflag
       end subroutine fcn
    end interface

    integer                             :: j
    real(dp)                            :: eps, epsmch, h, temp, wa(m)
    real(dp), parameter                 :: zero = 0.0_dp

    ! epsmch is the machine precision.

    epsmch = epsilon(zero)

    eps = sqrt(max(epsfcn, epsmch))

    do  j = 1, n
       temp = x(j)
       h = eps*abs(temp)
       if ( abs( h ) < machine_eps ) h = eps
       x(j) = temp + h
       call fcn(m, n, x, wa, iflag)
       if (iflag < 0) exit
       x(j) = temp
       fjac(1:m,j) = (wa(1:m) - fvec(1:m))/h
    end do

  end subroutine fdjac2


end module Levenberg_Marquardt
