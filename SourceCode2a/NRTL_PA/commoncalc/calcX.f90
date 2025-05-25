SUBROUTINE calcX(Y, ns, rho, xhost, noccur, del)
! Calculate fraction of sites that are unbonded
!
! Output:
!	Y		(ns)
!
! Input:
!	ns		(scalar)	number of sites
!	rho		(scalar)	mixture molar density [mol/cc]
!	xhost	(ns)		mole fraction of hosts of sites [-]
!	noccur	(ns)		number of occurance of a site in its host [-]
!	del		(ns x ns)	equilibrium constant [cc/mol]
!
IMPLICIT NONE
INTENT(IN)	::	ns, rho, xhost, noccur, del
INTENT(OUT)	::	Y

! comment the next four lines if not compiling for Aspen,
! as well as the calls to DMS_WRTTRM, DMS_ERRPRT,
! lines with global_ldiag, DMS_IRRCHK
#include "ppexec_user.cmn" ! for imiss
#include "dms_global.cmn" ! for diagnostic level
#include "dms_maxwrt.cmn" ! for writing to control panel
#include "dms_errout.cmn" ! for writing errors

INTEGER i, j, ns, nmaxc, nmaxr, iprt, DMS_IRRCHK, iprog(2)
INTEGER, DIMENSION(ns)	::	noccur
REAL*8, DIMENSION(ns, ns)	::	del
REAL*8, DIMENSION(ns)	::	vec, Y, Y_old, xhost
REAL*8 Y_diff, rho
CHARACTER*80 NTEXT, FMT
LOGICAL, PARAMETER	::	debug = .False. ! .True. for debugging mode, set .False. to mute

CALL debug_print_start

nmaxc = 6 ! print for up to six sites
vec = rho * xhost * noccur
Y = 2D0 / (1D0 + DSQRT(1D0 + 4D0 * MATMUL(vec, TRANSPOSE(del))))

i = 0
Y_diff = 1D0
DO WHILE (Y_diff > 1D-5)
	IF (i > 1000) THEN
        ! MAXWRT_MAXBUF is size (20,80)
        nmaxc = 6
        if (ns < 6) nmaxc = ns
        WRITE(MAXWRT_MAXBUF(1), 1000)
1000    FORMAT("X FAILED TO CONVERGE IN 1000 ITERATIONS. Up to 6 columns. X=")
        write(FMT,'("(",I0,"(G12.4))")') nmaxc ! make format string for ns columns
        WRITE(NTEXT, FMT) Y(1:nmaxc)
        MAXWRT_MAXBUF(2) = NTEXT
        WRITE(MAXWRT_MAXBUF(3), 1010) ns
1010    FORMAT("N SITES = ", I4)
        MAXWRT_MAXBUF(4)= "OCCURENCES"
        WRITE(MAXWRT_MAXBUF(5), 1020) noccur
1020    FORMAT(I3)
        MAXWRT_MAXBUF(6) = "RHO (CC/MOL)"
        write(FMT,'("(",I0,"(G12.4))")') nmaxc ! make format string for ns columns
        WRITE(NTEXT, FMT) rho
        MAXWRT_MAXBUF(7) = NTEXT
        MAXWRT_MAXBUF(8) = "HOST MOL FRACTION"
        WRITE(NTEXT, FMT) xhost(1:nmaxc)
        MAXWRT_MAXBUF(9) = NTEXT
        MAXWRT_MAXBUF(10) = "DELTA MATRIX (first 10 rows, up to 6 columns)"
        nmaxr = 10 ! at most 10 more rows
        if (ns < nmaxr) nmaxr = ns
        DO j = 1, nmaxr
            WRITE(NTEXT, FMT) del(j,1:nmaxc)
            MAXWRT_MAXBUF(10+j) = NTEXT
        END DO
        CALL DMS_WRTTRM(nmaxr+10)
        ! this is a severe error
		iprog = 'calcX   '
        if(DMS_IRRCHK(iprog, 0, 1, USER_LMSG, USER_IUMISS, 0, 1, 2) .ne. 0) then
			! copy MAXWRT buffer
			ERROUT_IEROUT(1:nmaxr+10) = MAXWRT_MAXBUF(1:nmaxr+10)
			CALL DMS_ERRPRT(nmaxr+10)
			Y = -1D0 ! not converged
			EXIT
		END IF
	END IF
	i = i + 1
	Y_old = Y;
	! If Y_succ is the successive substitution value, and half the step is used,
	! Y_new = Y + 0.5(Y_succ - Y) = 0.5*(Y + Y_succ)
	Y = 0.5D0 * (Y + 1D0 / (1D0 + matmul(vec * Y, transpose(del))))
	Y_diff = maxval(abs(Y-Y_old)/Y_old);
	CALL debug_print_loop
END DO

CALL debug_print_end

CONTAINS
	SUBROUTINE debug_print_start
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'CALL calcX',REPEAT('-',40)
			PRINT *, '>>>>Input:'
			PRINT *, 'ns=', ns
			PRINT *, 'rho=', rho
			PRINT *, 'xhost=', xhost
			PRINT *, 'noccur=', noccur
			PRINT *, 'del='
			DO i = 1, ns
				PRINT *, del(i,:)
			END DO
			PRINT *, '>>>>end of Input'; PRINT *, ''
		END IF
	END SUBROUTINE debug_print_start

	SUBROUTINE debug_print_loop
		IF (debug) THEN
			PRINT *, 'i='; PRINT *, i
			PRINT *, 'Y_diff='; PRINT *, Y_diff
			PRINT *, 'Y='; PRINT *, Y
		END IF
	END SUBROUTINE debug_print_loop

	SUBROUTINE debug_print_end
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'ENDOF calcX',REPEAT('-',40)
		ENDIF
	END SUBROUTINE debug_print_end
END SUBROUTINE calcX
