subroutine stavgugg(gamma,x,r,q,n)
implicit none
intent(out):: gamma
intent(in):: x, n, r, q
integer n
real*8 gamma(n), x(n), r(n), q(n), phioxi(n), sumxr, sumxq, thetaoxi(n)

!This function calculates Staverman-Guggenheim term

sumxr = dot_product(x,r)
sumxq = dot_product(x,q)

!Note: x cancels in the numerator and denominator of phi/theta so here we
!calculate the terms phi and theta over x to avoid a division by zero
phioxi = r/sumxr
thetaoxi = q/sumxq

gamma = -5D0*q*(log(phioxi/thetaoxi)+ 1D0- phioxi/thetaoxi);

return

end subroutine stavgugg

