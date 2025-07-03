subroutine correqn1262(gamma, x, V, n, bvol)
implicit none
intent(out):: gamma
intent(in):: x, n, V
integer n
real*8 gamma(n), x(n), V(n), bvol(n), Vmix, bmix, denom(n), num(n)

!This function calculates the correction to Flory's combinatorial term as calculated in Bala lab book pg. 119

Vmix = dot_product(x,V) !cc/mol
bmix = dot_product(bvol,x)

denom = 1D0-bvol/V
num = 1D0-bmix/Vmix

gamma = -log(num/denom)+(1D0/num)*(bvol/Vmix - bmix*V/Vmix**2)

return
end subroutine correqn1262
