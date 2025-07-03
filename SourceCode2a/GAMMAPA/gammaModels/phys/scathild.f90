subroutine scathild(gamma, x, V, n, Aij, R, T)
implicit none
intent(out):: gamma
intent(in):: x, n, V, Aij , R, T
integer n
real*8 gamma(n), x(n), V(n), Vmix, phi(1,n), R, T, Aij(n,n)
real*8 term1(1,1), term2(1,n), RTlnGamma(1,n)

!This function calculates scatchard-hildebrand residual term

Vmix = dot_product(x,V)
phi = (reshape(x,shape(phi))*reshape(V,shape(phi)))/Vmix
term1 = 0.5D0 * matmul(phi,matmul(Aij,transpose(phi)));
term2 = transpose(matmul(Aij,transpose(phi)));
RTlnGamma = reshape(V,shape(phi))*(term2-term1(1,1));
! gamma is actually ln(gamma) for aspens
gamma = reshape((RTlnGamma/R/T),shape(gamma));
return

end subroutine scathild
