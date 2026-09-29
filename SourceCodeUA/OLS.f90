MODULE OLS
! Small ordinary-least-squares helper required by FuEsd26a.
! Coefficients are returned in ascending powers: c(0)+c(1)x+...
IMPLICIT NONE
CONTAINS

SUBROUTINE polyFit(x, y, degree, coeff, stdErr, fixedIntercept)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: degree
    DOUBLE PRECISION, INTENT(IN) :: x(:), y(:)
    DOUBLE PRECISION, INTENT(OUT) :: coeff(0:), stdErr(0:)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: fixedIntercept
    DOUBLE PRECISION, ALLOCATABLE :: normal(:,:), rhs(:), solution(:), basis(:), inverseColumn(:)
    DOUBLE PRECISION :: fitted, residual, sse, sigma2
    INTEGER :: n, firstPower, nUnknown, i, j, k, info, dof

    n = SIZE(x)
    coeff = 0.D0
    stdErr = 0.D0
    IF (SIZE(y) /= n .OR. degree < 0 .OR. SIZE(coeff) < degree+1 .OR. SIZE(stdErr) < degree+1) RETURN

    firstPower = 0
    IF (PRESENT(fixedIntercept)) THEN
        coeff(0) = fixedIntercept
        firstPower = 1
    END IF
    nUnknown = degree - firstPower + 1
    IF (nUnknown <= 0) RETURN

    ALLOCATE(normal(nUnknown,nUnknown), rhs(nUnknown), solution(nUnknown), basis(nUnknown), inverseColumn(nUnknown))
    normal = 0.D0
    rhs = 0.D0
    DO i = 1, n
        DO j = 1, nUnknown
            basis(j) = x(i)**(firstPower+j-1)
        END DO
        DO j = 1, nUnknown
            rhs(j) = rhs(j) + basis(j)*(y(i)-coeff(0))
            DO k = 1, nUnknown
                normal(j,k) = normal(j,k) + basis(j)*basis(k)
            END DO
        END DO
    END DO

    CALL solveLinear(normal, rhs, solution, info)
    IF (info /= 0) THEN
        coeff = 0.D0
        stdErr = HUGE(1.D0)
        DEALLOCATE(normal, rhs, solution, basis, inverseColumn)
        RETURN
    END IF
    DO j = 1, nUnknown
        coeff(firstPower+j-1) = solution(j)
    END DO

    sse = 0.D0
    DO i = 1, n
        fitted = 0.D0
        DO j = 0, degree
            fitted = fitted + coeff(j)*x(i)**j
        END DO
        residual = y(i)-fitted
        sse = sse + residual*residual
    END DO
    dof = MAX(1, n-nUnknown)
    sigma2 = sse/DBLE(dof)

    DO j = 1, nUnknown
        rhs = 0.D0
        rhs(j) = 1.D0
        CALL solveLinear(normal, rhs, inverseColumn, info)
        IF (info == 0) stdErr(firstPower+j-1) = SQRT(MAX(0.D0,sigma2*inverseColumn(j)))
    END DO
    DEALLOCATE(normal, rhs, solution, basis, inverseColumn)
END SUBROUTINE polyFit

SUBROUTINE solveLinear(matrix, vector, solution, info)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: matrix(:,:), vector(:)
    DOUBLE PRECISION, INTENT(OUT) :: solution(:)
    INTEGER, INTENT(OUT) :: info
    DOUBLE PRECISION, ALLOCATABLE :: a(:,:), b(:), rowBuffer(:)
    DOUBLE PRECISION :: factor, pivotValue, valueBuffer
    INTEGER :: n, i, j, pivot

    n = SIZE(vector)
    info = 0
    ALLOCATE(a(n,n), b(n), rowBuffer(n))
    a = matrix
    b = vector
    DO i = 1, n
        pivot = i
        DO j = i+1, n
            IF (ABS(a(j,i)) > ABS(a(pivot,i))) pivot = j
        END DO
        IF (ABS(a(pivot,i)) <= EPSILON(1.D0)*MAX(1.D0,MAXVAL(ABS(a)))) THEN
            info = i
            solution = 0.D0
            DEALLOCATE(a,b,rowBuffer)
            RETURN
        END IF
        IF (pivot /= i) THEN
            rowBuffer = a(i,:)
            a(i,:) = a(pivot,:)
            a(pivot,:) = rowBuffer
            valueBuffer = b(i)
            b(i) = b(pivot)
            b(pivot) = valueBuffer
        END IF
        pivotValue = a(i,i)
        DO j = i+1, n
            factor = a(j,i)/pivotValue
            a(j,i:n) = a(j,i:n)-factor*a(i,i:n)
            b(j) = b(j)-factor*b(i)
        END DO
    END DO
    DO i = n, 1, -1
        solution(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),solution(i+1:n)))/a(i,i)
    END DO
    DEALLOCATE(a,b,rowBuffer)
END SUBROUTINE solveLinear

END MODULE OLS
