subroutine VLU(V,IDI,VEQN,DNLDIP,T)

use CONSTANTS, only: outfile

implicit none

intent(out):: V
intent(in)::IDI, T, VEQN, DNLDIP

! V		molar volume [cm3/mol]
! IDI	component idx value (order in component list) [-]
! T		temperature [K]
! VEQN	DIPPR equation [-]
! DNLDIP DIPPR equation constants
!   DNLDIP(6) is lower T limit
!	DNLDIP(7) is upper temperature limit
! Outside the limits, the values are adjust to the limit values.


integer I,IDI,VEQN, LTC, IPROG
real*8 V, T, TC, T1, DNLDIP(7), TAU
T1 = T ! Assign value for calculation that is not returned if changed
if (T.gt.DNLDIP(7)) T1 = DNLDIP(7)
if (T.lt.DNLDIP(6)) T1 = DNLDIP(6)
Select Case (VEQN)
	Case (105)
		V = 1D3* (DNLDIP(2)**(1D0 + (1D0-T1/DNLDIP(3))**DNLDIP(4)))/DNLDIP(1) ! cm3/mol
	Case (116)
		TC = DNLDIP(7) ! critical temperature K for this case
		tau = 1D0-T1/TC
		V = 1D3/(DNLDIP(1) + DNLDIP(2)*tau**(3.5D-1) + DNLDIP(3)*tau**(2D0/3D0) &
		    + DNLDIP(4)*tau + DNLDIP(5)*tau**(4D0/3D0)) ! cm3/mol
	Case (100)
		V = 1D3/(DNLDIP(1) + DNLDIP(2)*T1 + DNLDIP(3)*T1**2 + DNLDIP(4)*T1**3 + DNLDIP(5)*T1**5)
	Case DEFAULT
		WRITE(outfile,99) IDI
 99		FORMAT('Density is not programmed for component',I4)
		WRITE(outfile,100) VEQN
 100	FORMAT ('Equation ',I4, ' is not programmed. Contact the .dll authors,')
End Select
end subroutine vlu
