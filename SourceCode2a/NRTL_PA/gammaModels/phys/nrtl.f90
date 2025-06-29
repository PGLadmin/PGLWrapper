subroutine nrtl(gamma, n, x, tau, alpha)
implicit none
intent(out):: gamma
intent(in):: n, x,tau,alpha
integer n
real*8 gamma(n), x(n), tau(n,n), alpha(n,n), G(n,n), Y(n,n)
real*8 term1(n), term2(n), inverseterm2(n), squareinverseterm2(n), part1(n), part2(n), part3(n), loggamma(n)

! This program calculates activity coefficients using the NRTL equation.
! Each call to nrtl evaluates one composition.
! x = vector of n mole fractions for an n-component mixture, length n.
! tau = tau values, square matrix, n x n.
! alpha = alpha values, square matrix, n x n.
!   The alpha matrix is symmetrical with zeros on the diagonals.


!this will create a matrix with all the rows equal to the vector of
!compositions
! Y = kron(x,ones(nComp,1));
Y = spread(x,1,n)

! G=exp(-alpha.*tau);
G=exp(-alpha*tau);

! term1=x*(tau.*G);
! term2=x*G;
! inverseterm2 = (1./term2);
! squareinverseterm2=inverseterm2.^2;
term1=matmul(x,(tau*G))
term2=matmul(x,G)
inverseterm2 = (1D0/term2)
squareinverseterm2=inverseterm2**2


! part1=(term1./term2);
part1=(term1/term2)
! part2=inverseterm2*((Y.*tau.*G)');
part2= matmul(inverseterm2,transpose(Y*tau*G));
! part3=(term1.*squareinverseterm2)*((Y.*G)');
part3 = matmul(term1*squareinverseterm2,transpose(Y*G));

gamma = part1 + part2 - part3;

return

end subroutine nrtl
