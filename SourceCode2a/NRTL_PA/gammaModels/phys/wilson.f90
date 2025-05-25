! This program calculates activity coefficients using the wilson equation.
! Each call to wilson evaluates one composition.
! x = vector of n mole fractions for an n-component mixture, length n.
! Lambda = Lambda values, square matrix, n x n. Diagonals are UNITY.
! Lambda is calculated before this routine is used. In that way, the
! user can include whatever parameter temperature dependence is desired
! without modifying this routine.

SUBROUTINE wilson(gamma,n,xin,Lambda)
IMPLICIT NONE
INTENT(IN)::n,Lambda,xin
INTENT(OUT)::gamma

! declare variables
INTEGER n
REAL*8 Lambda(n,n),Y(n,n)
REAL*8 term1(n),inverseterm1(n),term2(n),gamma(n),x(n),xin(n)


! make sure x is not identically zero
! ctl - Looking at the formula, I don't think this is necessary
! x = xin+1e-37
x = xin

Y = spread(x,1,n)

term1 = matmul(x,transpose(Lambda))
inverseterm1 = 1D0/term1

term2 = matmul(inverseterm1,(transpose(Y)*Lambda))

!Note: the term 'gamma' is actually ln(gamma).
gamma = 1D0 - log(term1) - term2


END SUBROUTINE wilson
