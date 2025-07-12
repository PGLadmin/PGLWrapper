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

INTEGER i, j, ns, nmaxc, nmaxr
INTEGER, DIMENSION(ns)	::	noccur
REAL*8, DIMENSION(ns, ns)	::	del
REAL*8, DIMENSION(ns)	::	vec, Y, Y_old, xhost
REAL*8 Y_diff, rho
CHARACTER*80 NTEXT, FMT
LOGICAL	::	debug = .False. ! .True. for debugging mode, set .False. to mute

CALL debug_print_start

nmaxc = 6 ! print for up to six sites
vec = rho * xhost * noccur
Y = 2D0 / (1D0 + DSQRT(1D0 + 4D0 * MATMUL(vec, TRANSPOSE(del))))

i = 0
Y_diff = 1D0
DO WHILE (Y_diff > 1D-5)
	IF (i > 1000) THEN
		debug = .True.
		Print *, 'CalcX did not converge in 1000 iterations'
		Call debug_print_start
		Call debug_print_loop
		Call debug_print_end
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
