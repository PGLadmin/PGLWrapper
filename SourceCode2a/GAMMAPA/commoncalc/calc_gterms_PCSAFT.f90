SUBROUTINE calc_gterms_PCSAFT(g, rhogdgdrho, ngdgdnk, n, T, rho, x, sigma, m, epsok)
!
! Output:
!   g           (n x n) radial distribution function [-]
!   rhogdgdrho  (rho/g)(dg/drho)(n x n) [-]
!   ngdgdnk     (n/g)(dg/dnk)(n x n x n) [-]
!
! Input:
!   n           (scalar)  number of components
!   T           (scalar)  temperature [K]
!   rho         (scalar)  mixture molar density [mol/cc]
!   x           (n)      mole fraction of component
!   sigma       (n)      segment diameter [A]
!   m           (n)      number of segments per chain
!   epsok       (n)      depth of pair potential over boltzman constant [K]
!
! Parameters:
!   cf          (scalar)  conversion factor for mole density to number density [cc/mol/A^3]
!   xi2         (scalar)  abbreviation defined in Gross2001 eqn9 [A^-1]
!   xi3         (scalar)  abbreviation defined in Gross2001 eqn9 [-]
!   nxi2nk      (n)      n * molar derivative for xi2 [A^-1]
!   nxi3nk      (n)      n * molar derivative for xi3 [-]
!   d           (n)      temperature dependent segment diameter [A]
!   dmat        (n x n) [A]
!   rhodgdrho   (n x n) [-]
!   ndgdnk      (n x n x n) [-]

IMPLICIT NONE
INTENT(IN)      ::      n, T, rho, x, sigma, m, epsok
INTENT(OUT)     ::      g, rhogdgdrho, ngdgdnk
INTEGER n
LOGICAL, Parameter	::	 debug = .false. ! .True. for debugging mode, set .False. to mute
REAL*8  T, rho, cf, pi, xi2, xi3
REAL*8, DIMENSION(n)           ::      x, sigma, m, epsok, d, nxi2nk, nxi3nk
REAL*8, DIMENSION(n, n)       ::      dmat, rhodgdrho, g, rhogdgdrho
REAL*8, DIMENSION(n, n, n)   ::      ndgdnk, ngdgdnk

cf = 0.60221408D0
pi = 3.14159265D0

CALL debug_print_start

d = sigma * (1D0 - 1.2D-1 * DEXP(-3D0 * epsok / T))
nxi2nk = calc_nxink(2)
nxi3nk = calc_nxink(3)
xi2 = DOT_PRODUCT(x, nxi2nk)
xi3 = DOT_PRODUCT(x, nxi3nk)
CALL calc_dmat()
CALL calc_rhodgdrho()
CALL calc_ndgdnk()
CALL calc_g()
rhogdgdrho = rhodgdrho / g
ngdgdnk = ndgdnk / SPREAD(g, 3, n) ! add a third dimension to g by replicating (n x n) g in the third dimension

call debug_print_end

CONTAINS
        FUNCTION calc_nxink(l)
                INTEGER l
                REAL*8, DIMENSION(n)        ::      dl, calc_nxink
                dl = d ** l
                calc_nxink = pi / 6D0 * rho * cf * m * dl
        END FUNCTION calc_nxink

        SUBROUTINE calc_dmat()
                REAL*8, DIMENSION(n, n)       ::      d_spread
                d_spread = SPREAD(d, 1, n) ! replicate row1 n times
                dmat = d_spread * TRANSPOSE(d_spread) / (d_spread + TRANSPOSE(d_spread))
        END SUBROUTINE calc_dmat

        SUBROUTINE calc_rhodgdrho()
                REAL*8 c1, c2, c3
                c1 = xi3 / (1D0 - xi3)**2
                c2 = 3D0 * xi2 / (1D0 - xi3)**2 + 6D0 * xi2 * xi3 / (1D0 - xi3)**3
                c3 = 4D0 * xi2**2 / (1D0 - xi3)**3 + 6D0 * xi2**2 * xi3 / (1D0 - xi3)**4
                rhodgdrho = c1 + c2 * dmat + c3 * dmat**2
        END SUBROUTINE calc_rhodgdrho

        SUBROUTINE calc_ndgdnk()
                REAL*8, DIMENSION(n)   ::      c1, c2, c3
                REAL*8, DIMENSION(n, n, n)       :: c1_spread, c2_spread, c3_spread
                REAL*8, DIMENSION(n, n, n)   :: dmat_spread
                c1 = nxi3nk / (1D0 - xi3)**2 ! vector of deriviatives
                c2 = 3D0 * nxi2nk / (1D0 - xi3)**2 + 6D0 * xi2 * nxi3nk / (1D0 -xi3)**3
                c3 = 4D0 * xi2 * nxi2nk / (1D0 - xi3)**3 + 6D0 * xi2**2 * nxi3nk / (1D0 - xi3)**4
                c1_spread = SPREAD(SPREAD(c1,1,n),1,n) ! Result is an n x n) x n where
                c2_spread = SPREAD(SPREAD(c2,1,n),1,n) ! each (n x n) is a copy of c1(k), c2(k), c3(k).
                c3_spread = SPREAD(SPREAD(c3,1,n),1,n) ! The first to indices are dgij and the third is dnk.
                dmat_spread = spread(dmat, 3, n)
                ndgdnk = c1_spread + c2_spread * dmat_spread + c3_spread * dmat_spread**2
        END SUBROUTINE calc_ndgdnk

        SUBROUTINE calc_g()
                REAL*8 c1, c2, c3
                c1 = 1D0 / (1D0 - xi3)
                c2 = 3D0 * xi2 / (1D0 - xi3)**2
                c3 = 2D0 * xi2**2 / (1D0 - xi3)**3
                g = c1 + c2 * dmat + c3 * dmat**2
        END SUBROUTINE calc_g
        ! n, T, rho, x, sigma, m, epsok
        SUBROUTINE debug_print_start
		INTEGER i
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'CALL calc_gterms_PCSAFT',REPEAT('-',40)
			PRINT *, '>>>>Input:'
			PRINT *, 'n='; PRINT *, n
			PRINT *, 'T='; PRINT *, T
			PRINT *, 'rho='; PRINT *, rho
			PRINT *, 'x='; PRINT *, x
			PRINT *, 'sigma='; PRINT *, sigma
			PRINT *, 'm='; PRINT *, m
			PRINT *, 'epsok='; PRINT *, epsok
			PRINT *, '>>>>end of Input'; PRINT *, ''
		END IF
	END SUBROUTINE debug_print_start

	SUBROUTINE debug_print_end
		INTEGER i
		IF (debug) THEN
			PRINT *, 'g='
			DO i = 1, n
				PRINT *, g(i,:)
			END DO
			PRINT *, 'rhogdgdrho='
			DO i = 1, n
				PRINT *, rhogdgdrho(i,:)
			END DO
			PRINT *, 'ngdgdnk='
			DO i = 1, n
				PRINT *, ngdgdnk(i,:,:)
			END DO
			PRINT *, REPEAT('-',40),'ENDOF calc_gammaw',REPEAT('-',40)
		END IF
	END SUBROUTINE debug_print_end
END SUBROUTINE calc_gterms_PCSAFT
