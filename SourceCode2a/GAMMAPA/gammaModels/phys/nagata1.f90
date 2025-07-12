SUBROUTINE nagata1(gamma, n, x, rpara, a)
IMPLICIT NONE
INTENT(IN)::x, rpara, a
INTENT(OUT)::gamma

! This program calculates activity coefficients gamma based on
! Nagata, 1997, Fluid Phase Equilibria, 227-247, Equation (20)
! n = number of components, integer scalar.
! x = vector of n mole fractions for an n-component mixture, shape n.
! rpara = vector of molecular geometric-size parameter of pure components.
! a = adjustable parameters, square matrix with zeros on diagonal, shape n,n.
! T = temperature in Kelvin.

! declear variables
INTEGER n;
REAL*8 r_mix;
REAL*8 x(n),rpara(n);
REAL*8 term2(n),phi(n),term1(n),lnGamma(n),gamma(n);
REAL*8 a(n,n),tau(n,n);

! mixture molecular geometric-size parameter
r_mix = dot_product(x,rpara);
! segment fraction
phi = x*rpara/r_mix;
tau = dexp(a);

term1 = matmul(transpose(tau),phi);
term2 = matmul(tau,x/term1);
gamma = -dlog(term1) + rpara * (1D0-term2) / r_mix;

END SUBROUTINE nagata1
