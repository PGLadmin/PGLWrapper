subroutine VLU(V,IDI,VEQN,DNLDIP,T)
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

#include "dms_plex.cmn"
#include "ppexec_user.cmn" ! for imiss, user_lmsg
#include "dms_errout.cmn" ! for writing errors
#include "dms_maxwrt.cmn" ! for writing to terminal

        REAL*8 B(1)
        EQUIVALENCE (B(1), IB(1))
INTEGER DMS_IFCMNC, DMS_IRRCHK
integer I,IDI,VEQN, LTC, IPROG
real*8 V, T, TC, T1, DNLDIP(7), TAU
T1 = T ! Assign value for calculation that is not returned if changed
if (T.gt.DNLDIP(7)) T1 = DNLDIP(7)
if (T.lt.DNLDIP(6)) T1 = DNLDIP(6)
Select Case (VEQN)
	Case (105)
		V = 1D3* (DNLDIP(2)**(1D0 + (1D0-T1/DNLDIP(3))**DNLDIP(4)))/DNLDIP(1) ! cm3/mol
	Case (116)
		LTC = DMS_IFCMNC('TC')
		TC = B(LTC + IDI) ! critical temperature K
		tau = 1D0-T1/TC
		V = 1D3/(DNLDIP(1) + DNLDIP(2)*tau**(3.5D-1) + DNLDIP(3)*tau**(2D0/3D0) &
		    + DNLDIP(4)*tau + DNLDIP(5)*tau**(4D0/3D0)) ! cm3/mol
	Case (100)
		V = 1D3/(DNLDIP(1) + DNLDIP(2)*T1 + DNLDIP(3)*T1**2 + DNLDIP(4)*T1**3 + DNLDIP(5)*T1**5)
	Case DEFAULT
	    ! need to call IRRCHK to set error severity.
		! Arguments: routine name, severity, user err code, diagnostic level, unused, unused,
		! indicate physical property calculation, use normal header.
		! Setting as an error that will terminate.
		IPROG = 'VLU     '
		I = DMS_IRRCHK(IPROG,0,1,USER_LMSG,USER_IUMISS,0,1,2)
		! Value returned by IRRCHK is not used since all errors of this case are fatal.
		WRITE(ERROUT_IEROUT(1),99) IDI
 99		FORMAT('Density is not programmed for component',I4)
		WRITE(ERROUT_IEROUT(2),100) VEQN
 100	FORMAT ('Equation ',I4, ' is not programmed. Contact the .dll authors,')
        ERROUT_IEROUT(3) = 'or enter DNLDIP constants and set THRSWT(2) to 100, 105, or 116.'
		MAXWRT_MAXBUF(1:3) = ERROUT_IEROUT(1:3)
		CALL DMS_WRTTRM(3)
        CALL DMS_ERRPRT(3)
End Select
end subroutine vlu
