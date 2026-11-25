MODULE OLS
  IMPLICIT NONE

  ! Generic interface to overload the name "polyFit()" The required version is inferred from the number of arguments.
  INTERFACE polyFit
	MODULE PROCEDURE polyFit_basic			!(x, y, n, coeff)
	MODULE PROCEDURE polyFit_withStats		!(x, y, n, coeff, stdErr)
	MODULE PROCEDURE polyFit_specIntercept	!(x, y, n, coeff, stdErr, yint): explicit y-intercept is set by user.
  END INTERFACE

CONTAINS

  subroutine ols_fit(X, y, beta, stderr, r2)
    ! Inputs
    real(8), intent(in)  :: X(:,:)   ! design matrix (n × p)
    real(8), intent(in)  :: y(:)     ! response vector (n)
    ! Outputs
    real(8), intent(out) :: beta(:)  ! coefficients (p)
    real(8), intent(out) :: stderr(:)! standard errors (p)
    real(8), intent(out) :: r2       ! coefficient of determination

    ! Locals
    integer :: n, p
    real(8) :: XtX(size(X,2), size(X,2)), Xty(size(X,2))
    real(8) :: yhat(size(y)), resid(size(y))
    real(8) :: sigma2
    real(8) :: cov(size(X,2), size(X,2))

    n = size(X,1)
    p = size(X,2)

    ! Normal equations
    XtX = matmul(transpose(X), X)
    Xty = matmul(transpose(X), y)

    ! Coefficients
    beta = matmul(inv(XtX), Xty)

    ! Fitted values and residuals
    yhat = matmul(X, beta)
    resid = y - yhat

    ! Error variance estimate
    sigma2 = sum(resid**2) / real(n - p,8)

    ! Variance–covariance matrix of coefficients
    cov = sigma2 * inv(XtX)

    ! Standard errors
    stderr = sqrt(diagonal(cov))

    ! R²
    r2 = 1.0d0 - sum(resid**2) / sum((y - sum(y)/n)**2)
  end subroutine ols_fit

  ! Simple matrix inversion (Gauss–Jordan)
  function inv(A) result(Ainv)
    real(8), intent(in) :: A(:,:)
    real(8) :: Ainv(size(A,1), size(A,2))
    real(8) :: aug(size(A,1), 2*size(A,2))
    integer :: i, j, n

    n = size(A,1)
    aug(:,1:n) = A
    aug(:,n+1:2*n) = 0.0d0
    do i=1,n
       aug(i,n+i) = 1.0d0
    end do

    ! Gauss–Jordan elimination
    do i=1,n
       aug(i,:) = aug(i,:) / aug(i,i)
       do j=1,n
          if (j /= i) aug(j,:) = aug(j,:) - aug(j,i)*aug(i,:)
       end do
    end do
    Ainv = aug(:,n+1:2*n)
  end function inv

  ! Extract diagonal of a square matrix
  function diagonal(A) result(d)
    real(8), intent(in) :: A(:,:)
    real(8) :: d(size(A,1))
    integer :: i
    do i=1,size(A,1)
       d(i) = A(i,i)
    end do
  end function diagonal

	! Version 1: coefficients only, intercept regressed
	SUBROUTINE polyFit_basic(x, y, n, coeff)
	  REAL*8, INTENT(IN) :: x(:), y(:)
	  INTEGER, INTENT(IN) :: n
	  REAL*8, INTENT(OUT) :: coeff(:)
	  REAL*8 :: Xmat(SIZE(x), n+1)
	  REAL*8 :: DUMMY_STDERR(n+1), DUMMY_R2
	  INTEGER :: i

	  if(n > 6)then
		write(*,*)'polyFit_basic:Sorry. n=6 is the maximum.'
		return
	  endif
	  ! Build design matrix with intercept
	  DO i = 1, SIZE(x)
		Xmat(i,1) = 1.0D0
		Xmat(i,2) = x(i)
		IF (n > 1) Xmat(i,3) = x(i)**2
		IF (n > 2) Xmat(i,4) = x(i)**3
		IF (n > 3) Xmat(i,5) = x(i)**4
		IF (n > 4) Xmat(i,6) = x(i)**5
		IF (n > 5) Xmat(i,7) = x(i)**6
	  END DO

	  CALL OLS_fit(Xmat, y, coeff, DUMMY_STDERR, DUMMY_R2)
	END SUBROUTINE polyFit_basic


	! Version 2: coefficients + stdErr, intercept regressed
	SUBROUTINE polyFit_withStats(x, y, n, coeff, stdErr)
	  REAL*8, INTENT(IN) :: x(:), y(:)
	  INTEGER, INTENT(IN) :: n
	  REAL*8, INTENT(OUT) :: coeff(:), stdErr(:)
	  REAL*8 :: Xmat(SIZE(x), n+1)
	  REAL*8 :: R2
	  INTEGER :: i

	  if(n > 6)then
		write(*,*)'polyFit_withStats:Sorry. n=6 is the maximum.'
		return
	  endif
	  ! Build design matrix with intercept
	  DO i = 1, SIZE(x)
		Xmat(i,1) = 1.0D0
		Xmat(i,2) = x(i)
		IF (n > 1) Xmat(i,3) = x(i)**2
		IF (n > 2) Xmat(i,4) = x(i)**3
		IF (n > 3) Xmat(i,5) = x(i)**4
		IF (n > 4) Xmat(i,6) = x(i)**5
		IF (n > 5) Xmat(i,7) = x(i)**6
	  END DO

	  CALL OLS_fit(Xmat, y, coeff, stdErr, R2)
	END SUBROUTINE polyFit_withStats


	! Version 3: coefficients + stdErr + explicit intercept
	SUBROUTINE polyFit_specIntercept(x, y, n, coeff, stdErr, yint)
	  REAL*8, INTENT(IN) :: x(:), y(:)
	  INTEGER, INTENT(IN) :: n
	  REAL*8, INTENT(OUT) :: coeff(:), stdErr(:)
	  REAL*8, INTENT(IN)  :: yint
	  REAL*8 :: Xmat(SIZE(x), n)
	  REAL*8 :: beta(n), stderr_terms(n), R2   !beta and stdErr are reduced	to n because the intercept is not a variable.
	  INTEGER :: i

	  if(n > 6)then
		write(*,*)'polyFit_specIntercept:Sorry. n=6 is the maximum.'
		return
	  endif
	  ! Build design matrix without intercept
	  DO i = 1, SIZE(x)
		Xmat(i,1) = x(i)
		IF (n > 1) Xmat(i,2) = x(i)**2
		IF (n > 2) Xmat(i,3) = x(i)**3
		IF (n > 3) Xmat(i,4) = x(i)**4
		IF (n > 4) Xmat(i,5) = x(i)**5
		IF (n > 5) Xmat(i,6) = x(i)**6
	  END DO

	  ! Run regression on y - yint
	  CALL OLS_fit(Xmat, y - yint, beta, stderr_terms, R2)

	  ! Assemble coefficient vector
	  coeff(1) = yint
	  coeff(2:) = beta(:)

	  ! Assemble stdErr vector
	  stdErr(1) = 0.0d0            ! fixed intercept ? no uncertainty
	  stdErr(2:) = stderr_terms(:) ! polynomial terms
	END SUBROUTINE polyFit_specIntercept

END MODULE OLS
