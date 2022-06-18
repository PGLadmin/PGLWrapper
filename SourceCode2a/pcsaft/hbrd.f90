!> \file hbrd.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE Solve_NonLin
!> \brief Solve_NonLin
!!
!! \todo Detailed description.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module Solve_NonLin

  ! Corrections to FUNCTION Enorm - 28 November 2003

  use PARAMETERS, only: dp, machine_eps
  implicit none

  private
  public  :: hbrd, hybrd


contains


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE HBRD
  !> \brief  HBRD solve N nonlinear eqs. in N variables by modified POWELL HYBRID METHOD.
  !!
  !! THE PURPOSE OF HBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR! FUNCTIONS
  !! IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.  THIS IS DONE
  !! BY USING THE MORE GENERAL NONLINEAR EQUATION SOLVER HYBRD. USER MUST PROVIDE
  !! A SR WHICH CALCULATES THE FUNCTIONS.
  !! THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!   SUBROUTINE HBRD(N, X, FVEC, EPSFCN, TOL, INFO, WA, LWA)
  !!
  !! SUBPROGRAMS CALLED
  !!   USER-SUPPLIED ...... FCN
  !!   MINPACK-SUPPLIED ... HYBRD
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !!
  !! Reference:
  !! Powell, M.J.D. 'A hybrid method for nonlinear equations' in Numerical Methods
  !!      for Nonlinear Algebraic Equations', P.Rabinowitz (editor), Gordon and
  !!      Breach, London 1970.
  !!
  ! WHERE FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
  !      THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT
  !      IN THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
  !           SUBROUTINE FCN(N, X, FVEC, IFLAG)
  !      THE VALUE OF IFLAG NOT BE CHANGED BY FCN UNLESS
  !      THE USER WANTS TO TERMINATE THE EXECUTION OF HBRD.
  !      IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
  !
  ! N    POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !      OF FUNCTIONS AND VARIABLES.
  !
  ! X    ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN AN INITIAL
  !      ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE
  !      FINAL ESTIMATE OF THE SOLUTION VECTOR.
  ! FVEC OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
  !      THE FUNCTIONS EVALUATED AT THE OUTPUT X.
  ! EPSFCN: INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
  !      FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION ASSUMES
  !      THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF EPSFCN.
  !      IF EPSFCN IS LESS THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE
  !      RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
  !      PRECISION.
  !
  ! TOL  NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
  !      ALGORITHM ESTIMATES THAT THE RELATIVE ERROR BETWEEN X AND THE SOLUTION
  !      IS AT MOST TOL.
  ! INFO INTEGER OUTPUT VARIABLE.  IF THE USER HAS TERMINATED
  !      EXECUTION, INFO IS SET TO THE (NEGATIVE) VALUE OF IFLAG.
  !      SEE DESCRIPTION OF FCN.  OTHERWISE, INFO IS SET AS FOLLOWS.
  !      INFO = 0   IMPROPER INPUT PARAMETERS.
  !      INFO = 1   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
  !                 BETWEEN X AND THE SOLUTION IS AT MOST TOL.
  !      INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED 200*(N+1).
  !      INFO = 3   TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
  !                 THE APPROXIMATE SOLUTION X IS POSSIBLE.
  !      INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS.
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine hbrd(fcn, n, x, fvec, epsfcn, tol, info, diag, N_max_itr )

    integer, intent(IN)               :: n          !< POSITIVE INTEGER INPUT VARIABLE SET TO NUMBER OF FUNCTIONS AND VARIABLES
    real(dp), intent(IN OUT)          :: x(n)       !< ON INPUT X MUST CONTAINS INITIAL VALUES.\n  ON OUTPUT X CONTAINS SOLUTION VECTOR
    real(dp), intent(IN OUT)          :: fvec(n)    !< OUTPUT ARRAY OF LENGTH N WHICH CONTAIN THE FUNCTIONS EVALUATED AT THE OUTPUT X
    real(dp), intent(IN)              :: epsfcn     !< INPUT VARIABLE, USED FOR STEP-LENGTH OF FORWARD-DIFFERENCE SCHEME
    real(dp), intent(IN)              :: tol        !< INPUT VARIABLE.  TERMINATION OCCURS WHEN RELATIVE ERROR BETWEEN X AND THE SOLUTION IS < TOL
    integer, intent(OUT)              :: info       !< INTEGER OUTPUT VARIABLE.  SEE FILE HBRD.F90 FOR DETAILS.
    real(dp), intent(OUT)             :: diag(n)
    integer, intent(IN), optional     :: N_max_itr

    ! EXTERNAL fcn
    interface
       subroutine fcn(N, X, FVEC, IFLAG)
         implicit none
         integer, parameter           :: dp = selected_real_kind(15, 307)  !< double precision
         integer, intent(IN)          :: n
         real(dp), intent(IN)         :: x(n)
         real(dp), intent(OUT)        :: fvec(n)
         integer, intent(IN OUT)      :: iflag
       end subroutine fcn
    end interface

    integer                           :: maxfev, ml, mode, mu, nfev, nprint
    real(dp)                          :: xtol
    real(dp), parameter               :: factor = 1.0_dp, zero = 0.0_dp

    info = 0

    ! CHECK THE INPUT PARAMETERS FOR ERRORS.
    if (n > 0 .and. epsfcn >= zero .and. tol >= zero) then

       ! CALL HYBRD.

       if ( present( N_max_itr ) ) then
          maxfev = N_max_itr
       else
          maxfev = 200*(n + 1)
       end if

       xtol = tol
       ml = n - 1
       mu = n - 1
       mode = 2
       nprint = 0
       call hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
            factor, nprint, info, nfev)
       if (info == 5) info = 4

    end if

  end subroutine hbrd



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE HYBRD
  !> \brief THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
  !!
  !! THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
  !! FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
  !! THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
  !! THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!
  !!   SUBROUTINE HYBRD(FCN, N, X, FVEC, XTOL, MAXFEV, ML, MU, EPSFCN,
  !!                    DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, FJAC,
  !!                    LDFJAC, R, LR, QTF, WA1, WA2, WA3, WA4)
  !! SUBPROGRAMS CALLED
  !!
  !!   USER-SUPPLIED ...... FCN
  !!
  !!   MINPACK-SUPPLIED ... DOGLEG,SPMPAR,ENORM,FDJAC1,
  !!                        QFORM,QRFAC,R1MPYQ,R1UPDT
  !!
  !!   FORTRAN-SUPPLIED ... ABS,MAX,MIN,MIN,MOD
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !
  ! WHERE
  !
  !   FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
  !     THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
  !     THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
  !          SUBROUTINE FCN(N, X, FVEC, IFLAG)
  !     THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
  !     THE USER WANTS TO TERMINATE EXECUTION OF HYBRD.
  !     IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !     OF FUNCTIONS AND VARIABLES.
  !
  !   X IS AN ARRAY OF LENGTH N.  ON INPUT X MUST CONTAIN AN INITIAL
  !     ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE FINAL
  !     ESTIMATE OF THE SOLUTION VECTOR.
  !
  !   FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
  !     THE FUNCTIONS EVALUATED AT THE OUTPUT X.
  !
  !   XTOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
  !     RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES IS AT MOST XTOL.
  !
  !   MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE.  TERMINATION OCCURS WHEN
  !     THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN
  !     ITERATION.
  !
  !   ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE
  !     NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
  !     IF THE JACOBIAN IS NOT BANDED, SET ML TO AT LEAST N - 1.
  !
  !   MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE NUMBER
  !     OF SUPERDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
  !     IF THE JACOBIAN IS NOT BANDED, SET MU TO AT LEAST N - 1.
  !
  !   EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
  !     FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION
  !     ASSUMES THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER
  !     OF EPSFCN. IF EPSFCN IS LESS THAN THE MACHINE PRECISION,
  !     IT IS ASSUMED THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE
  !     ORDER OF THE MACHINE PRECISION.
  !
  !   DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE BELOW),
  !     DIAG IS INTERNALLY SET.  IF MODE = 2, DIAG MUST CONTAIN POSITIVE
  !     ENTRIES THAT SERVE AS MULTIPLICATIVE SCALE FACTORS FOR THE
  !     VARIABLES.
  !
  !   MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE VARIABLES WILL BE
  !     SCALED INTERNALLY.  IF MODE = 2, THE SCALING IS SPECIFIED BY THE
  !     INPUT DIAG.  OTHER VALUES OF MODE ARE EQUIVALENT TO MODE = 1.
  !
  !   FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
  !     INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
  !     FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
  !     TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
  !     INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.
  !
  !   NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
  !     PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
  !     FCN IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
  !     ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
  !     IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
  !     FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
  !     OF FCN WITH IFLAG = 0 ARE MADE.
  !
  !   INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
  !     TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
  !     VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
  !     INFO IS SET AS FOLLOWS.
  !
  !     INFO = 0   IMPROPER INPUT PARAMETERS.
  !
  !     INFO = 1   RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
  !                IS AT MOST XTOL.
  !
  !     INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED MAXFEV.
  !
  !     INFO = 3   XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
  !                THE APPROXIMATE SOLUTION X IS POSSIBLE.
  !
  !     INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS, AS
  !                MEASURED BY THE IMPROVEMENT FROM THE LAST
  !                FIVE JACOBIAN EVALUATIONS.
  !
  !     INFO = 5   ITERATION IS NOT MAKING GOOD PROGRESS, AS MEASURED BY
  !                THE IMPROVEMENT FROM THE LAST TEN ITERATIONS.
  !
  !   NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.
  !
  !   FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE ORTHOGONAL MATRIX Q
  !     PRODUCED BY THE QR FACTORIZATION OF THE FINAL APPROXIMATE JACOBIAN.
  !
  !   LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
  !     WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
  !
  !   R IS AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE
  !     UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
  !     OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.
  !
  !   LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN (N*(N+1))/2.
  !
  !   QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
  !     THE VECTOR (Q TRANSPOSE)*FVEC.
  !
  !   WA1, WA2, WA3, AND WA4 ARE WORK ARRAYS OF LENGTH N.
  !
  ! SUBPROGRAMS CALLED
  !
  !   USER-SUPPLIED ...... FCN
  !
  !   MINPACK-SUPPLIED ... DOGLEG,SPMPAR,ENORM,FDJAC1,
  !                        QFORM,QRFAC,R1MPYQ,R1UPDT
  !
  !   FORTRAN-SUPPLIED ... ABS,MAX,MIN,MIN,MOD
  !
  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
       factor, nprint, info, nfev)

    integer, intent(IN)               :: n          !< POSITIVE INTEGER INPUT VARIABLE SET TO NUMBER OF FUNCTIONS AND VARIABLES
    real(dp), intent(IN OUT)          :: x(n)       !< ON INPUT X MUST CONTAINS INITIAL VALUES.\n  ON OUTPUT X CONTAINS SOLUTION VECTOR
    real(dp), intent(IN OUT)          :: fvec(n)    !< OUTPUT ARRAY OF LENGTH N WHICH CONTAIN THE FUNCTIONS EVALUATED AT THE OUTPUT X
    real(dp), intent(IN)              :: xtol       !< INPUT VARIABLE. TERMINATION OCCURS WHEN RELATIVE ERROR BETWEEN 2 ITERATES IS < XTOL
    integer, intent(IN OUT)           :: maxfev     !< INPUT VARIABLE. TERMINATION OCCURS WHEN THE NUMBER OF CALLS TO FCN IS > MAXFEV
    integer, intent(IN OUT)           :: ml         !< INPUT VARIABLE, SPECIFIES NUMBER OF SUBDIAGONALS WITHIN BAND OF JACOBIAN. IF JACOBIAN ISN'T BANDED, SET ML>= N-1
    integer, intent(IN)               :: mu         !< INPUT VARIABLE, SPECIFIES NUMBER OF SUPERDIAGONALS WITHIN BAND OF JACOBIAN. IF JACOBIAN ISN'T BANDED, SET MU>= N-1
    real(dp), intent(IN)              :: epsfcn     !< INPUT VARIABLE, USED FOR STEP-LENGTH OF FORWARD-DIFFERENCE SCHEME
    real(dp), intent(OUT)             :: diag(n)    !< IF MODE = 1, DIAG IS INTERNALLY SET. IF MODE = 2, DIAG MUST CONTAIN POSITIVE ENTRIES AS DESCRIBED IN HBRD.F90
    integer, intent(IN)               :: mode       !< IF MODE = 1, THE VARIABLES WILL BE SCALED INTERNALLY. IF MODE = 2, THE SCALING IS SPECIFIED BY INPUT DIAG
    real(dp), intent(IN)              :: factor     !< INPUT VARIABLE USED FOR INITIAL STEP BOUND. SEE HBRD.F90 FOR DETAILS. VALUE OF 100 IS A RECOMMENDED VALUE
    integer, intent(IN OUT)           :: nprint     !< INPUT VARIABLE FOR PRINTING. SEE HBRD.F90 FOR DETAILS
    integer, intent(OUT)              :: info       !< INTEGER OUTPUT VARIABLE.  SEE FILE HBRD.F90 FOR DETAILS
    integer, intent(OUT)              :: nfev       !< OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN

    ! EXTERNAL fcn
    interface
       subroutine FCN(N, X, FVEC, IFLAG)
         implicit none
         integer, parameter           :: dp = selected_real_kind(15, 307)  !< double precision
         integer, intent(IN)          :: n
         real(dp), intent(IN)         :: x(n)
         real(dp), intent(OUT)        :: fvec(n)
         integer, intent(IN OUT)      :: iflag
       end subroutine FCN
    end interface

    integer                           :: i, iflag, iter, j
    integer                           :: jm1, l, lr, msum
    integer                           :: ncfail, ncsuc, nslow1, nslow2
    integer                           :: iwa(1)
    logical                           :: jeval, sing
    real(dp)                          :: actred, delta, epsmch
    real(dp)                          :: fnorm, fnorm1, pnorm
    real(dp)                          :: prered, ratio, sum, temp, xnorm
    real(dp), parameter               :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp
    real(dp), parameter               :: p001 = 0.001_dp, p0001 = 0.0001_dp, zero = 0.0_dp

    ! The following were workspace arguments
    real(dp)                          :: fjac(n,n), r(n*(n+1)/2), qtf(n)
    real(dp)                          :: wa1(n), wa2(n), wa3(n), wa4(n)

    ! EPSMCH IS THE MACHINE PRECISION.

    epsmch = epsilon( 1.0_dp )

    info = 0
    iflag = 0
    nfev = 0
    lr = n*(n+1)/2

    ! CHECK THE INPUT PARAMETERS FOR ERRORS.

    if (n > 0 .and. xtol >= zero .and. maxfev > 0 .and. ml >= 0 .and. mu >=  &
         0 .and. factor > zero ) then
       if (mode == 2) then
          diag(1:n) = one
       end if

       ! EVALUATE THE FUNCTION AT THE STARTING POINT AND CALCULATE ITS NORM.

       iflag = 1
       call fcn(n, x, fvec, iflag)
       nfev = 1
       if (iflag >= 0) then
          fnorm = enorm(n, fvec)

          ! DETERMINE THE NUMBER OF CALLS TO FCN NEEDED TO COMPUTE THE JACOBIAN MATRIX.

          msum = min(ml+mu+1,n)

          ! INITIALIZE ITERATION COUNTER AND MONITORS.

          iter = 1
          ncsuc = 0
          ncfail = 0
          nslow1 = 0
          nslow2 = 0

          ! BEGINNING OF THE OUTER LOOP.

20        jeval = .true.

          ! CALCULATE THE JACOBIAN MATRIX.

          iflag = 2
          call fdjac1(fcn, n, x, fvec, fjac, n, iflag, ml, mu, epsfcn, wa1, wa2)
          nfev = nfev + msum
          if (iflag >= 0) then

             ! COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.

             call qrfac(n, n, fjac, n, .false., iwa, 1, wa1, wa2, wa3)

             ! ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
             ! TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.

             if (iter == 1) then
                if (mode /= 2) then
                   do  j = 1, n
                      diag(j) = wa2(j)
                      if ( abs( wa2(j) ) < machine_eps ) diag(j) = one
                   end do
                end if

                ! ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
                ! AND INITIALIZE THE STEP BOUND DELTA.

                wa3(1:n) = diag(1:n) * x(1:n)
                xnorm = enorm(n, wa3)
                delta = factor * xnorm
                if ( abs( delta ) < machine_eps ) delta = factor
             end if

             ! FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.

             qtf(1:n) = fvec(1:n)
             do  j = 1, n
                if ( abs( fjac(j,j) ) > machine_eps ) then
                   sum = zero
                   do  i = j, n
                      sum = sum + fjac(i,j) * qtf(i)
                   end do
                   temp = -sum / fjac(j,j)
                   do  i = j, n
                      qtf(i) = qtf(i) + fjac(i,j) * temp
                   end do
                end if
             end do

             ! COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.

             sing = .false.
             do  j = 1, n
                l = j
                jm1 = j - 1
                if (jm1 >= 1) then
                   do  i = 1, jm1
                      r(l) = fjac(i,j)
                      l = l + n - i
                   end do
                end if
                r(l) = wa1(j)
                if ( abs( wa1(j) ) < machine_eps ) sing = .true.
             end do

             ! ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.

             call qform(n, n, fjac, n, wa1)

             ! RESCALE IF NECESSARY.

             if (mode /= 2) then
                do  j = 1, n
                   diag(j) = max(diag(j), wa2(j))
                end do
             end if

             ! BEGINNING OF THE INNER LOOP.

             ! IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.

120          if (nprint > 0) then
                iflag = 0
                if (mod(iter-1, nprint) == 0) call fcn(n, x, fvec, iflag)
                if (iflag < 0) GO TO 190
             end if

             ! DETERMINE THE DIRECTION P.

             call dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3)

             ! STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.

             do  j = 1, n
                wa1(j) = -wa1(j)
                wa2(j) = x(j) + wa1(j)
                wa3(j) = diag(j) * wa1(j)
             end do
             pnorm = enorm(n, wa3)

             ! ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.

             if (iter == 1) delta = min(delta, pnorm)

             ! EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.

             iflag = 1
             call fcn(n, wa2, wa4, iflag)
             nfev = nfev + 1
             if (iflag >= 0) then
                fnorm1 = enorm(n, wa4)

                ! COMPUTE THE SCALED ACTUAL REDUCTION.

                actred = -one
                if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm) ** 2

                ! COMPUTE THE SCALED PREDICTED REDUCTION.

                l = 1
                do  i = 1, n
                   sum = zero
                   do  j = i, n
                      sum = sum + r(l) * wa1(j)
                      l = l + 1
                   end do
                   wa3(i) = qtf(i) + sum
                end do
                temp = enorm(n, wa3)
                prered = zero
                if (temp < fnorm) prered = one - (temp/fnorm) ** 2

                ! COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED REDUCTION.

                ratio = zero
                if (prered > zero) ratio = actred / prered

                ! UPDATE THE STEP BOUND.

                if (ratio < p1) then
                   ncsuc = 0
                   ncfail = ncfail + 1
                   delta = p5 * delta
                else
                   ncfail = 0
                   ncsuc = ncsuc + 1
                   if (ratio >= p5 .or. ncsuc > 1) delta = max(delta,pnorm/p5)
                   if (abs(ratio-one) <= p1) delta = pnorm / p5
                end if

                ! TEST FOR SUCCESSFUL ITERATION.

                if (ratio >= p0001) then

                   ! SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.

                   do  j = 1, n
                      x(j) = wa2(j)
                      wa2(j) = diag(j) * x(j)
                      fvec(j) = wa4(j)
                   end do
                   xnorm = enorm(n, wa2)
                   fnorm = fnorm1
                   iter = iter + 1
                end if

                ! DETERMINE THE PROGRESS OF THE ITERATION.

                nslow1 = nslow1 + 1
                if (actred >= p001) nslow1 = 0
                if (jeval) nslow2 = nslow2 + 1
                if (actred >= p1) nslow2 = 0

                ! TEST FOR CONVERGENCE.

                if (delta <= xtol*xnorm .or. abs( fnorm ) < machine_eps ) info = 1
                if (info == 0) then

                   ! TESTS FOR TERMINATION AND STRINGENT TOLERANCES.

                   if (nfev >= maxfev) info = 2
                   if (p1*max(p1*delta, pnorm) <= epsmch*xnorm) info = 3
                   if (nslow2 == 5) info = 4
                   if (nslow1 == 10) info = 5
                   if (info == 0) then

                      ! CRITERION FOR RECALCULATING JACOBIAN APPROXIMATION
                      ! BY FORWARD DIFFERENCES.

                      if (ncfail /= 2) then

                         ! CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
                         ! AND UPDATE QTF IF NECESSARY.

                         do  j = 1, n
                            sum = zero
                            do  i = 1, n
                               sum = sum + fjac(i,j) * wa4(i)
                            end do
                            wa2(j) = (sum-wa3(j)) / pnorm
                            wa1(j) = diag(j) * ((diag(j)*wa1(j))/pnorm)
                            if (ratio >= p0001) qtf(j) = sum
                         end do

                         ! COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.

                         call r1updt(n, n, r, lr, wa1, wa2, wa3, sing)
                         call r1mpyq(n, n, fjac, n, wa2, wa3)
                         call r1mpyq(1, n, qtf, 1, wa2, wa3)

                         ! END OF THE INNER LOOP.

                         jeval = .false.
                         GO TO 120
                      end if

                      ! END OF THE OUTER LOOP.

                      GO TO 20
                   end if
                end if
             end if
          end if
       end if
    end if

    ! TERMINATION, EITHER NORMAL OR USER IMPOSED.

190 if (iflag < 0) info = iflag
    iflag = 0
    if (nprint > 0) call fcn(n, x, fvec, iflag)

  end subroutine hybrd


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE DOGLEG
  !> \brief THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
  !!  !! GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
  !! MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA, THE
  !! PROBLEM IS TO DETERMINE THE CONVEX COMBINATION X OF THE
  !! GAUSS-NEWTON AND SCALED GRADIENT DIRECTIONS THAT MINIMIZES
  !! (A*X - B) IN THE LEAST SQUARES SENSE, SUBJECT TO THE
  !! RESTRICTION THAT THE EUCLIDEAN NORM OF D*X BE AT MOST DELTA.
  !!
  !! THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
  !! IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
  !! QR FACTORIZATION OF A. THAT IS, IF A = Q*R, WHERE Q HAS
  !! ORTHOGONAL COLUMNS AND R IS AN UPPER TRIANGULAR MATRIX,
  !! THEN DOGLEG EXPECTS THE FULL UPPER TRIANGLE OF R AND
  !! THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!
  !!   SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
  !!   SUBPROGRAMS CALLED
  !!
  !!   MINPACK-SUPPLIED ... SPMPAR,ENORM
  !!
  !!   FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !
  ! WHERE
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
  !
  !   R IS AN INPUT ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
  !     TRIANGULAR MATRIX R STORED BY ROWS.
  !
  !   LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
  !     (N*(N+1))/2.
  !
  !   DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
  !     DIAGONAL ELEMENTS OF THE MATRIX D.
  !
  !   QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
  !     N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
  !
  !   DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
  !     BOUND ON THE EUCLIDEAN NORM OF D*X.
  !
  !   X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
  !     CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION AND THE
  !     SCALED GRADIENT DIRECTION.
  !
  !   WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.
  !
  ! SUBPROGRAMS CALLED
  !
  !   MINPACK-SUPPLIED ... SPMPAR,ENORM
  !
  !   FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
  !
  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine dogleg(n, r, lr, diag, qtb, delta, x, wa1, wa2)

    integer, intent(IN)           :: n       !< INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
    integer, intent(IN)           :: lr      !< INTEGER INPUT VARIABLE NOT LESS THAN (N*(N+1))/2.
    real(dp), intent(IN)          :: r(lr)   !< INPUT ARRAY OF LENGTH LR, CONTAINS UPPER TRIANGULAR MATRIX R
    real(dp), intent(IN)          :: diag(n) !< INPUT ARRAY OF LENGTH N, CONTAINS DIAGONAL ELEMENTS OF MATRIX D.
    real(dp), intent(IN)          :: qtb(n)  !< INPUT ARRAY OF LENGTH N, CONTAINS FIRST N ELEMENTS OF VECTOR (Q TRANSPOSE)*B
    real(dp), intent(IN)          :: delta   !< INPUT VARIABLE, SPECIFIES AN UPPER BOUND ON EUCLIDEAN NORM OF D*X
    real(dp), intent(IN OUT)      :: x(n)    !< OUTPUT ARRAY OF LENGTH N, CONTAINS DESIRED DIRECTION
    real(dp), intent(OUT)         :: wa1(n)  !< WORK ARRAY OF LENGTH N.
    real(dp), intent(OUT)         :: wa2(n)  !< WORK ARRAY OF LENGTH N.

    integer                       :: i, j, jj, jp1, k, l
    real(dp)                      :: alpha, bnorm, epsmch
    real(dp)                      :: gnorm, qnorm, sgnorm, sum, temp

    ! EPSMCH IS THE MACHINE PRECISION.

    epsmch = epsilon(1.0_dp)

    ! FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.

    jj = (n*(n+1)) / 2 + 1
    do  k = 1, n
       j = n - k + 1
       jp1 = j + 1
       jj = jj - k
       l = jj + 1
       sum = 0.0_dp
       if (n >= jp1) then
          do  i = jp1, n
             sum = sum + r(l) * x(i)
             l = l + 1
          end do
       end if
       temp = r(jj)
       if ( abs( temp ) < machine_eps ) then
          l = j
          do  i = 1, j
             temp = max(temp,abs(r(l)))
             l = l + n - i
          end do
          temp = epsmch * temp
          if ( abs( temp ) < machine_eps ) temp = epsmch
       end if
       x(j) = (qtb(j)-sum) / temp
    end do

    ! TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.

    do  j = 1, n
       wa1(j) = 0.0_dp
       wa2(j) = diag(j) * x(j)
    end do
    qnorm = enorm(n, wa2)
    if (qnorm > delta) then

       ! THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
       ! NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.

       l = 1
       do  j = 1, n
          temp = qtb(j)
          do  i = j, n
             wa1(i) = wa1(i) + r(l) * temp
             l = l + 1
          end do
          wa1(j) = wa1(j) / diag(j)
       end do

       ! CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
       ! THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.

       gnorm = enorm(n, wa1)
       sgnorm = 0.0_dp
       alpha = delta / qnorm
       if ( abs( gnorm ) > machine_eps ) then

          ! CALCULATE THE POINT ALONG THE SCALED GRADIENT
          ! AT WHICH THE QUADRATIC IS MINIMIZED.

          do  j = 1, n
             wa1(j) = (wa1(j)/gnorm) / diag(j)
          end do
          l = 1
          do  j = 1, n
             sum = 0.0_dp
             do  i = j, n
                sum = sum + r(l) * wa1(i)
                l = l + 1
             end do
             wa2(j) = sum
          end do
          temp = enorm(n, wa2)
          sgnorm = (gnorm/temp) / temp

          ! TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.

          alpha = 0.0_dp
          if (sgnorm < delta) then

             ! THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
             ! FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
             ! AT WHICH THE QUADRATIC IS MINIMIZED.

             bnorm = enorm(n, qtb)
             temp = (bnorm/gnorm) * (bnorm/qnorm) * (sgnorm/delta)
             temp = temp - (delta/qnorm) * (sgnorm/delta) ** 2 + sqrt((  &
                  temp-(delta/qnorm))**2+(1.0_dp-(delta/qnorm)**2)*(1.0_dp-( sgnorm/delta)**2))
             alpha = ((delta/qnorm)*(1.0_dp-(sgnorm/delta)**2)) / temp
          end if
       end if

       ! FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
       ! DIRECTION AND THE SCALED GRADIENT DIRECTION.

       temp = (1.0_dp-alpha) * min(sgnorm,delta)
       do  j = 1, n
          x(j) = temp * wa1(j) + alpha * x(j)
       end do
    end if

  end subroutine dogleg


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE FDJAC1
  !> \brief THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
  !!
  !! THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION TO THE N BY N
  !! JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED PROBLEM OF N FUNCTIONS IN N
  !! VARIABLES.  IF THE JACOBIAN HAS A BANDED FORM, THEN FUNCTION EVALUATIONS
  !! ARE SAVED BY ONLY APPROXIMATING THE NONZERO TERMS.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!
  !!   SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,WA2)
  !!
  !! SUBPROGRAMS CALLED
  !!
  !!   MINPACK-SUPPLIED ... SPMPAR
  !!
  !!   FORTRAN-SUPPLIED ... ABS,MAX,SQRT
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !
  ! WHERE
  !
  !   FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
  !     THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
  !     THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
  !
  !     SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !     INTEGER N,IFLAG
  !     REAL X(N),FVEC(N)
  !     ----------
  !     CALCULATE THE FUNCTIONS AT X AND
  !     RETURN THIS VECTOR IN FVEC.
  !     ----------
  !     END
  !
  !     THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
  !     THE USER WANTS TO TERMINATE EXECUTION OF FDJAC1.
  !     IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !     OF FUNCTIONS AND VARIABLES.
  !
  !   X IS AN INPUT ARRAY OF LENGTH N.
  !
  !   FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
  !     FUNCTIONS EVALUATED AT X.
  !
  !   FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
  !     APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
  !
  !   LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
  !     WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
  !
  !   IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
  !     THE EXECUTION OF FDJAC1.  SEE DESCRIPTION OF FCN.
  !
  !   ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
  !     THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
  !     JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
  !     ML TO AT LEAST N - 1.
  !
  !   EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
  !     STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
  !     APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
  !     FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
  !     THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
  !     ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE PRECISION.
  !
  !   MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
  !     THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
  !     JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
  !     MU TO AT LEAST N - 1.
  !
  !   WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.  IF ML + MU + 1 IS AT
  !     LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, AND WA2 IS
  !     NOT REFERENCED.
  !
  ! SUBPROGRAMS CALLED
  !
  !   MINPACK-SUPPLIED ... SPMPAR
  !
  !   FORTRAN-SUPPLIED ... ABS,MAX,SQRT
  !
  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn,   &
       wa1, wa2)

    integer, intent(IN)           :: n       !< INPUT VARIABLE SET TO THE NUMBER OF FUNCTIONS AND VARIABLES
    real(dp), intent(IN OUT)      :: x(n)    !< INPUT ARRAY OF LENGTH N
    real(dp), intent(IN)          :: fvec(n) !< INPUT ARRAY OF LENGTH N, CONTAINS THE FUNCTIONS EVALUATED AT X
    integer, intent(IN)           :: ldfjac  !< INPUT, NOT LESS THAN N WHICH SPECIFIES LEADING DIMENSION OF ARRAY FJAC
    real(dp), intent(OUT)         :: fjac(ldfjac,n) !< OUTPUT N x N ARRAY, CONTAINS APPROXIMATION TO JACOBIAN AT X
    integer, intent(IN OUT)       :: iflag   !< CAN BE USED TO TERMINATE THE EXECUTION OF FDJAC1.  SEE DESCRIPTION OF FCN.
    integer, intent(IN)           :: ml      !< INPUT, SPECIFIES NUMBER OF SUBDIAGONALS. IF JACOBIAN ISN'T BANDED, SET ML >= N - 1
    integer, intent(IN)           :: mu      !< INPUT, SPECIFIES NUMBER OF SUPERDIAGONALS. IF JACOBIAN ISN'T BANDED, SET MU >= N - 1
    real(dp), intent(IN)          :: epsfcn  !< INPUT, USED FOR SUITABLE STEP LENGTH OF FORWARD-DIFFERENCE SCHEME
    real(dp), intent(IN OUT)      :: wa1(n)  !< ARE WORK ARRAY OF LENGTH N. IF ML + MU + 1 >= N, JACOBIAN IS CONSIDERED DENSE
    real(dp), intent(OUT)         :: wa2(n)  !< ARE WORK ARRAY OF LENGTH N. IF ML + MU + 1 >= N, JACOBIAN IS CONSIDERED DENSE

    ! EXTERNAL fcn
    interface
       subroutine FCN(N, X, FVEC, IFLAG)
         implicit none
         integer, parameter           :: dp = selected_real_kind(15, 307)  !< double precision
         integer, intent(IN)          :: n
         real(dp), intent(IN)         :: x(n)
         real(dp), intent(OUT)        :: fvec(n)
         integer, intent(IN OUT)      :: iflag
       end subroutine FCN
    end interface

    integer                       :: i, j, k, msum
    real(dp)                      :: eps, epsmch, h, temp
    real(dp), parameter           :: zero = 0.0_dp

    ! EPSMCH IS THE MACHINE PRECISION.

    epsmch = epsilon(1.0_dp)

    eps = sqrt(max(epsfcn, epsmch))
    msum = ml + mu + 1
    if (msum >= n) then

       ! COMPUTATION OF DENSE APPROXIMATE JACOBIAN.

       do  j = 1, n
          temp = x(j)
          h = eps * abs(temp)
          if ( abs( h ) < machine_eps ) h = eps
          x(j) = temp + h
          call fcn(n, x, wa1, iflag)
          if (iflag < 0) exit
          x(j) = temp
          do  i = 1, n
             fjac(i,j) = (wa1(i)-fvec(i)) / h
          end do
       end do
    else

       ! COMPUTATION OF BANDED APPROXIMATE JACOBIAN.

       do  k = 1, msum
          do  j = k, n, msum
             wa2(j) = x(j)
             h = eps * abs(wa2(j))
             if ( abs( h ) < machine_eps ) h = eps
             x(j) = wa2(j) + h
          end do
          call fcn(n, x, wa1, iflag)
          if (iflag < 0) exit
          do  j = k, n, msum
             x(j) = wa2(j)
             h = eps * abs(wa2(j))
             if ( abs( h ) < machine_eps ) h = eps
             do  i = 1, n
                fjac(i,j) = zero
                if (i >= j-mu .and. i <= j+ml) fjac(i,j) = (wa1(i)-fvec(i)) / h
             end do
          end do
       end do
    end if

  end subroutine fdjac1


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE QFORM
  !> \brief QFORM
  !!
  !! THIS SUBROUTINE PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF AN M BY N
  !! MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX Q FROM ITS FACTORED FORM.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!
  !!   SUBROUTINE QFORM(M,N,Q,LDQ,WA)
  !!
  !! SUBPROGRAMS CALLED
  !!
  !!   FORTRAN-SUPPLIED ... MIN
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !
  ! WHERE
  !
  !   M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !     OF ROWS OF A AND THE ORDER OF Q.
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.
  !
  !   Q IS AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
  !     THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
  !     ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.
  !
  !   LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
  !     WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.
  !
  !   WA IS A WORK ARRAY OF LENGTH M.
  !
  ! SUBPROGRAMS CALLED
  !
  !   FORTRAN-SUPPLIED ... MIN
  !
  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine qform(m, n, q, ldq, wa)

    integer, intent(IN)           :: m        !< INPUT, SET TO NUMBER OF ROWS OF A AND ORDER OF Q
    integer, intent(IN)           :: n        !< INPUT, SET TO NUMBER OF COLUMNS OF A
    integer, intent(IN)           :: ldq      !< INPUT, NOT LESS THAN M, SPECIFIES LEADING DIMENSION OF ARRAY Q
    real(dp), intent(OUT)         :: q(ldq,m) !< AN M x M ARRAY. SEE HBRD.F90 FOR DETAILS.
    real(dp), intent(OUT)         :: wa(m)    !< IS A WORK ARRAY OF LENGTH M.


    integer                       :: i, j, jm1, k, l, minmn, np1
    real(dp)                      :: sum, temp
    real(dp), parameter           :: one = 1.0_dp, zero = 0.0_dp

    ! ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.

    minmn = min(m,n)
    if (minmn >= 2) then
       do  j = 2, minmn
          jm1 = j - 1
          do  i = 1, jm1
             q(i,j) = zero
          end do
       end do
    end if

    ! INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.

    np1 = n + 1
    if (m >= np1) then
       do  j = np1, m
          do  i = 1, m
             q(i,j) = zero
          end do
          q(j,j) = one
       end do
    end if

    ! ACCUMULATE Q FROM ITS FACTORED FORM.

    do  l = 1, minmn
       k = minmn - l + 1
       do  i = k, m
          wa(i) = q(i,k)
          q(i,k) = zero
       end do
       q(k,k) = one
       if ( abs( wa(k) ) > machine_eps ) then
          do  j = k, m
             sum = zero
             do  i = k, m
                sum = sum + q(i,j) * wa(i)
             end do
             temp = sum / wa(k)
             do  i = k, m
                q(i,j) = q(i,j) - temp * wa(i)
             end do
          end do
       end if
    end do

  end subroutine qform


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE QRFAC
  !> \brief COMPUTE A QR FACTORIZATION
  !!
  !! THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN PIVOTING
  !! (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE M BY N MATRIX A.
  !! THAT IS, QRFAC DETERMINES AN ORTHOGONAL MATRIX Q, A PERMUTATION MATRIX P,
  !! AND AN UPPER TRAPEZOIDAL MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING
  !! MAGNITUDE, SUCH THAT A*P = Q*R.  THE HOUSEHOLDER TRANSFORMATION FOR
  !! COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM
  !!
  !!                       T
  !!       I - (1/U(K))*U*U
  !!
  !! WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS.  THE FORM OF THIS
  !! TRANSFORMATION AND THE METHOD OF PIVOTING FIRST APPEARED IN THE
  !! CORRESPONDING LINPACK SUBROUTINE.
  !!
  !! THE SUBROUTINE STATEMENT IS
  !!
  !!   SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
  !!
  !! SUBPROGRAMS CALLED
  !!
  !!   MINPACK-SUPPLIED ... SPMPAR,ENORM
  !!
  !!   FORTRAN-SUPPLIED ... MAX,SQRT,MIN
  !!
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !
  ! WHERE
  !
  !   M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !     OF COLUMNS OF A.
  !
  !   A IS AN M BY N ARRAY.  ON INPUT A CONTAINS THE MATRIX FOR WHICH THE
  !     QR FACTORIZATION IS TO BE COMPUTED.  ON OUTPUT THE STRICT UPPER
  !     TRAPEZOIDAL PART OF A CONTAINS THE STRICT UPPER TRAPEZOIDAL PART OF R,
  !     AND THE LOWER TRAPEZOIDAL PART OF A CONTAINS A FACTORED FORM OF Q
  !     (THE NON-TRIVIAL ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).
  !
  !   LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
  !     WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
  !
  !   PIVOT IS A LOGICAL INPUT VARIABLE.  IF PIVOT IS SET TRUE,
  !     THEN COLUMN PIVOTING IS ENFORCED.  IF PIVOT IS SET FALSE,
  !     THEN NO COLUMN PIVOTING IS DONE.
  !
  !   IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT.  IPVT DEFINES THE
  !     PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
  !     COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
  !     IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.
  !
  !   LIPVT IS A POSITIVE INTEGER INPUT VARIABLE.  IF PIVOT IS FALSE,
  !     THEN LIPVT MAY BE AS SMALL AS 1.  IF PIVOT IS TRUE, THEN
  !     LIPVT MUST BE AT LEAST N.
  !
  !   RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
  !     DIAGONAL ELEMENTS OF R.
  !
  !   ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NORMS OF
  !     THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
  !     IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE WITH RDIAG.
  !
  !   WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
  !     CAN COINCIDE WITH RDIAG.
  !
  ! SUBPROGRAMS CALLED
  !
  !   MINPACK-SUPPLIED ... SPMPAR,ENORM
  !
  !   FORTRAN-SUPPLIED ... MAX,SQRT,MIN
  !
  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

    integer, intent(IN)            :: m           !< INPUT, SET TO THE NUMBER OF ROWS OF A
    integer, intent(IN)            :: n           !< INPUT, SET TO THE NUMBER OF COLUMNS OF A
    integer, intent(IN)            :: lda         !< INPUT, NOT LESS THAN M, SPECIFIES LEADING DIMENSION OF ARRAY A.
    real(dp), intent(IN OUT)       :: a(lda,n)    !< AN M x N ARRAY. SEE HBRD.F90 FOR DETAILS
    logical, intent(IN)            :: pivot       !< LOGICAL INPUT, IF PIVOT=TRUE, THEN COLUMN PIVOTING
    integer, intent(IN)            :: lipvt       !< INPUT, IF PIVOT=FALSE, LIPVT MAY BE AS SMALL AS 1. IF PIVOT=TRUE, LIPVT>=N
    integer, intent(OUT)           :: ipvt(lipvt) !< OUTPUT ARRAY(LIPVT). SEE HBRD.F90 FOR DETAILS
    real(dp), intent(OUT)          :: rdiag(n)    !< OUTPUT ARRAY(N), CONTAINS DIAGONAL ELEMENTS OF R.
    real(dp), intent(OUT)          :: acnorm(n)   !< OUTPUT ARRAY(N), CONTAINS NORMS OF COLUMNS OF INPUT MATRIX A
    real(dp), intent(OUT)          :: wa(n)       !< IS A WORK ARRAY(N). IF PIVOT=FALSE, WA CAN COINCIDE WITH RDIAG.

    integer                        :: i, j, jp1, k, kmax, minmn
    real(dp)                       :: ajnorm, epsmch, sum, temp
    real(dp), parameter            :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

    ! EPSMCH IS THE MACHINE PRECISION.

    epsmch = epsilon(1.0_dp)

    ! COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.

    do  j = 1, n
       acnorm(j) = enorm(m, a(1:,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       if (pivot) ipvt(j) = j
    end do

    ! REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.

    minmn = min(m,n)
    do  j = 1, minmn
       if (pivot) then

          ! BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.

          kmax = j
          do  k = j, n
             if (rdiag(k) > rdiag(kmax)) kmax = k
          end do
          if (kmax /= j) then
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
          end if
       end if

       ! COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
       ! J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.

       ajnorm = enorm(m-j+1, a(j:,j))
       if ( abs( ajnorm ) > machine_eps ) then
          if (a(j,j) < zero) ajnorm = -ajnorm
          do  i = j, m
             a(i,j) = a(i,j) / ajnorm
          end do
          a(j,j) = a(j,j) + one

          ! APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS AND UPDATE THE NORMS.

          jp1 = j + 1
          if (n >= jp1) then
             do  k = jp1, n
                sum = zero
                do  i = j, m
                   sum = sum + a(i,j) * a(i,k)
                end do
                temp = sum / a(j,j)
                do  i = j, m
                   a(i,k) = a(i,k) - temp * a(i,j)
                end do
                if ( .not.( .not.pivot .or. abs( rdiag(k) ) < machine_eps ) ) then
                   temp = a(j,k) / rdiag(k)
                   rdiag(k) = rdiag(k) * sqrt(max(zero,one-temp**2))
                   if (p05*(rdiag(k)/wa(k))**2 <= epsmch) then
                      rdiag(k) = enorm(m-j, a(jp1:,k))
                      wa(k) = rdiag(k)
                   end if
                end if
             end do
          end if
       end if
       rdiag(j) = -ajnorm
    end do

  end subroutine qrfac


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE R1MPYQ

  ! GIVEN AN M BY N MATRIX A, THIS SUBROUTINE COMPUTES A*Q WHERE
  ! Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

  !       GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

  ! AND GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
  ! ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
  ! Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
  ! GV, GW ROTATIONS IS SUPPLIED.

  ! THE SUBROUTINE STATEMENT IS

  !   SUBROUTINE R1MPYQ(M, N, A, LDA, V, W)

  ! WHERE

  !   M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

  !   A IS AN M BY N ARRAY.  ON INPUT A MUST CONTAIN THE MATRIX TO BE
  !     POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q DESCRIBED ABOVE.
  !     ON OUTPUT A*Q HAS REPLACED A.

  !   LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
  !     WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

  !   V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE INFORMATION
  !     NECESSARY TO RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

  !   W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE INFORMATION
  !     NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

  ! SUBROUTINES CALLED

  !   FORTRAN-SUPPLIED ... ABS, SQRT

  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine r1mpyq(m, n, a, lda, v, w)

    integer, intent(IN)            :: m
    integer, intent(IN)            :: n
    integer, intent(IN)            :: lda
    real(dp), intent(IN OUT)       :: a(lda,n)
    real(dp), intent(IN)           :: v(n)
    real(dp), intent(IN)           :: w(n)


    integer                        :: i, j, nmj, nm1
    real(dp)                       :: COS_r, SIN_r, temp
    real(dp), parameter            :: one = 1.0_dp

    ! APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.

    nm1 = n - 1
    if (nm1 >= 1) then
       do  nmj = 1, nm1
          j = n - nmj
          if (abs(v(j)) > one) then
             COS_r = one / v(j)
             SIN_r = sqrt(one-COS_r**2)
          else
             SIN_r = v(j)
             COS_r = sqrt(one-SIN_r**2)
          end if
          do  i = 1, m
             temp = COS_r * a(i,j) - SIN_r * a(i,n)
             a(i,n) = SIN_r * a(i,j) + COS_r * a(i,n)
             a(i,j) = temp
          end do
       end do

       ! APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.

       do  j = 1, nm1
          if (abs(w(j)) > one) then
             COS_r = one / w(j)
             SIN_r = sqrt(one-COS_r**2)
          else
             SIN_r = w(j)
             COS_r = sqrt(one-SIN_r**2)
          end if
          do  i = 1, m
             temp = COS_r * a(i,j) + SIN_r * a(i,n)
             a(i,n) = -SIN_r * a(i,j) + COS_r * a(i,n)
             a(i,j) = temp
          end do
       end do
    end if

  end subroutine r1mpyq



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! SUBROUTINE R1UPDT
  !> \brief R1UPDT
  !!
  !! GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
  !! AND AN N-VECTOR V, THE PROBLEM IS TO DETERMINE AN
  !! ORTHOGONAL MATRIX Q SUCH THAT
  !!               T
  !!       (S + U*V )*Q
  !! IS AGAIN LOWER TRAPEZOIDAL.
  !! THIS SUBROUTINE DETERMINES Q AS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS
  !!       GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
  !! WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
  !! WHICH ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
  !! Q ITSELF IS NOT ACCUMULATED, RATHER THE INFORMATION TO RECOVER THE GV,
  !! GW ROTATIONS IS RETURNED.
  !! THE SUBROUTINE STATEMENT IS
  !!   SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)
  !! SUBPROGRAMS CALLED
  !!   MINPACK-SUPPLIED ... SPMPAR
  !!   FORTRAN-SUPPLIED ... ABS,SQRT
  !! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  !! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE, JOHN L. NAZARETH

  ! WHERE

  !   M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF S.

  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !     OF COLUMNS OF S.  N MUST NOT EXCEED M.

  !   S IS AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
  !     TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
  !     THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.

  !   LS IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
  !     (N*(2*M-N+1))/2.

  !   U IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE VECTOR U.

  !   V IS AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR V.
  !     ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
  !     RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

  !   W IS AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
  !     NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

  !   SING IS A LOGICAL OUTPUT VARIABLE.  SING IS SET TRUE IF ANY OF THE
  !     DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO.  OTHERWISE SING IS
  !     SET FALSE.

  ! SUBPROGRAMS CALLED

  !   MINPACK-SUPPLIED ... SPMPAR

  !   FORTRAN-SUPPLIED ... ABS,SQRT

  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE, JOHN L. NAZARETH
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine r1updt(m, n, s, ls, u, v, w, sing)

    integer, intent(IN)            :: m       !< INPUT, SET TO NUMBER OF ROWS OF S
    integer, intent(IN)            :: n       !< INPUT, SET TO NUMBER OF COLUMNS OF S. N MUST <= M
    integer, intent(IN)            :: ls      !< INPUT VARIABLE NOT LESS THAN (N*(2*M-N+1))/2
    real(dp), intent(IN OUT)       :: s(ls)   !< ARRAY(LS). SEE HBRD.F90 FOR DETAILS
    real(dp), intent(IN)           :: u(m)    !< INPUT ARRAY(M), MUST CONTAIN THE VECTOR U
    real(dp), intent(IN OUT)       :: v(n)    !< IS ARRAY(N). SEE HBRD.F90 FOR DETAILS
    real(dp), intent(OUT)          :: w(m)    !< OUTPUT ARRAY(M). SEE HBRD.F90 FOR DETAILS
    logical, intent(OUT)           :: sing    !< LOGICAL OUTPUT, SING=TRUE IF ANY DIAGONAL ELEMENTS OF OUTPUT S ARE ZERO

    integer                        :: i, j, jj, l, nmj, nm1
    real(dp)                       :: COS, cotan, SIN, TAN, tau, temp
    real(dp)                       :: giant
    real(dp), parameter            :: one = 1.0_dp, p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

    ! GIANT IS THE LARGEST MAGNITUDE.

    giant = huge(1.0_dp)

    ! INITIALIZE THE DIAGONAL ELEMENT POINTER.

    jj = (n*(2*m-n+1)) / 2 - (m-n)

    ! MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.

    l = jj
    do  i = n, m
       w(i) = s(l)
       l = l + 1
    end do

    ! ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
    ! IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.

    nm1 = n - 1
    if (nm1 >= 1) then
       do  nmj = 1, nm1
          j = n - nmj
          jj = jj - (m-j+1)
          w(j) = zero
          if ( abs( v(j) ) > machine_eps ) then

             ! DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE J-TH ELEMENT OF V.

             if (abs(v(n)) < abs(v(j))) then
                cotan = v(n) / v(j)
                SIN = p5 / sqrt(p25+p25*cotan**2)
                COS = SIN * cotan
                tau = one
                if (abs(COS)*giant > one) tau = one / COS
             else
                TAN = v(j) / v(n)
                COS = p5 / sqrt(p25+p25*TAN**2)
                SIN = COS * TAN
                tau = SIN
             end if

             ! APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
             ! NECESSARY TO RECOVER THE GIVENS ROTATION.

             v(n) = SIN * v(j) + COS * v(n)
             v(j) = tau

             ! APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.

             l = jj
             do  i = j, m
                temp = COS * s(l) - SIN * w(i)
                w(i) = SIN * s(l) + COS * w(i)
                s(l) = temp
                l = l + 1
             end do
          end if
       end do
    end if

    ! ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.

    do  i = 1, m
       w(i) = w(i) + v(n) * u(i)
    end do

    ! ELIMINATE THE SPIKE.

    sing = .false.
    if (nm1 >= 1) then
       do  j = 1, nm1
          if ( abs( w(j) ) > machine_eps ) then

             ! DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
             ! J-TH ELEMENT OF THE SPIKE.

             if (abs(s(jj)) < abs(w(j))) then
                cotan = s(jj) / w(j)
                SIN = p5 / sqrt(p25 + p25*cotan**2)
                COS = SIN * cotan
                tau = one
                if (abs(COS)*giant > one) tau = one / COS
             else
                TAN = w(j) / s(jj)
                COS = p5 / sqrt(p25+p25*TAN**2)
                SIN = COS * TAN
                tau = SIN
             end if

             ! APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.

             l = jj
             do  i = j, m
                temp = COS * s(l) + SIN * w(i)
                w(i) = -SIN * s(l) + COS * w(i)
                s(l) = temp
                l = l + 1
             end do

             ! STORE THE INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION.

             w(j) = tau
          end if

          ! TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.

          if ( abs( s(jj) ) < machine_eps ) sing = .true.
          jj = jj + (m-j+1)
       end do
    end if

    ! MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.

    l = jj
    do  i = n, m
       s(l) = w(i)
       l = l + 1
    end do
    if ( abs( s(jj) ) < machine_eps ) sing = .true.

  end subroutine r1updt


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! FUNCTION ENORM

  ! GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE EUCLIDEAN NORM OF X.

  ! THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF SQUARES IN THREE
  ! DIFFERENT SUMS.  THE SUMS OF SQUARES FOR THE SMALL AND LARGE COMPONENTS
  ! ARE SCALED SO THAT NO OVERFLOWS OCCUR.  NON-DESTRUCTIVE UNDERFLOWS ARE
  ! PERMITTED.  UNDERFLOWS AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
  ! SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
  ! THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS DEPEND ON
  ! TWO CONSTANTS, RDWARF AND RGIANT.  THE MAIN RESTRICTIONS ON THESE CONSTANTS
  ! ARE THAT RDWARF**2 NOT UNDERFLOW AND RGIANT**2 NOT OVERFLOW.
  ! THE CONSTANTS GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.

  ! THE FUNCTION STATEMENT IS

  !   REAL FUNCTION ENORM(N, X)

  ! WHERE

  !   N IS A POSITIVE INTEGER INPUT VARIABLE.

  !   X IS AN INPUT ARRAY OF LENGTH N.

  ! SUBPROGRAMS CALLED

  !   FORTRAN-SUPPLIED ... ABS,SQRT

  ! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
  ! BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  function enorm(n, x) result(fn_val)

    integer, intent(IN)                 :: n
    real(dp), intent(IN)                :: x(n)
    real(dp)                            :: fn_val

    integer                             :: i
    real(dp)                            :: agiant, floatn
    real(dp)                            :: s1, s2, s3
    real(dp)                            :: xabs, x1max, x3max
    real(dp), parameter                 :: rdwarf = 1.0E-100_dp, rgiant = 1.0E+100_dp

    s1 = 0.0_dp
    s2 = 0.0_dp
    s3 = 0.0_dp
    x1max = 0.0_dp
    x3max = 0.0_dp
    floatn = n
    agiant = rgiant / floatn
    do  i = 1, n
       xabs = abs(x(i))
       if (xabs <= rdwarf .or. xabs >= agiant) then
          if (xabs > rdwarf) then

             ! SUM FOR LARGE COMPONENTS.

             if (xabs > x1max) then
                s1 = 1.0_dp + s1 * (x1max/xabs) ** 2
                x1max = xabs
             else
                s1 = s1 + (xabs/x1max) ** 2
             end if
          else

             ! SUM FOR SMALL COMPONENTS.

             if (xabs > x3max) then
                s3 = 1.0_dp + s3 * (x3max/xabs) ** 2
                x3max = xabs
             else
                if ( abs( xabs ) > machine_eps ) s3 = s3 + (xabs/x3max) ** 2
             end if
          end if
       else

          ! SUM FOR INTERMEDIATE COMPONENTS.

          s2 = s2 + xabs ** 2
       end if
    end do

    ! CALCULATION OF NORM.

    if ( abs( s1 ) > machine_eps ) then
       fn_val = x1max * sqrt(s1 + (s2/x1max)/x1max)
    else
       if ( abs( s2 ) > machine_eps ) then
          if ( s2 >= x3max ) then
             fn_val = sqrt( s2 * ( 1.0_dp + (x3max/s2)*(x3max*s3) ) )
          else
             fn_val = sqrt(x3max*((s2/x3max) + (x3max*s3)))
          end if
       else
          fn_val = x3max * sqrt(s3)
       end if
    end if

  end function enorm

end module Solve_NonLin
