subroutine flory(gamma,x,V,n)
implicit none
intent(out):: gamma
intent(in):: x, n, V
integer n
real*8 gamma(n), x(n), V(n), Vmix, phi(n)

!This function calculates Flory's term

Vmix = dot_product(x,V)
gamma = log(V/Vmix)+ 1D0- V/Vmix

return

end subroutine flory
