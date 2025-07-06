SUBROUTINE gcalc(del, gterm1, gterm2, gcalcFlag, n, ns, KADplus, &
		comp, site, T, rho_mix, bvol, rho_pure, sigma, m, epsok)
!
!
! Note: no input checkings are done, need to be done before calling the subroutine
USE sitenspecies, only: siteinfo, species ! for legacy reasons
IMPLICIT NONE

INTENT(OUT)	::	del, gterm1, gterm2
INTENT(IN)	::	gcalcFlag, n, ns, KADplus, comp, site, T, &
			rho_mix, bvol, rho_pure, sigma, m, epsok
!OPTIONAL	::	bvol, rho_pure, sigma, m, epsok ! inputs that depend on gcalcFlag

LOGICAL, Parameter	::	 debug = .False. ! .True. for debugging mode, set .False. to mute
INTEGER gcalcFlag, n, ns, i, k
REAL*8, DIMENSION(ns, ns, n)	::	gterm2
REAL*8, DIMENSION(ns, ns)	::	KADplus, del, gterm1
REAL*8, DIMENSION(n)		::	bvol, rho_pure, sigma, m, epsok
REAL*8 T, rho_mix
TYPE(species), DIMENSION(n)	::	comp
TYPE(siteinfo), DIMENSION(ns)	::	site

CALL debug_print_start

IF (gcalcFlag .EQ. 0) THEN ! VdW
	CALL gcalc_generic(1D0)
ELSEIF (gcalcFlag .EQ. 1) THEN ! ESD
	CALL gcalc_generic(1.9D0)
ELSEIF (gcalcFlag .EQ. 2) THEN ! CPA
	CALL gcalc_generic(1.9D0 / 4D0)
ELSEIF (gcalcFlag .EQ. 3) THEN ! PCSAFT
	CALL gcalc_PCSAFT
END IF

CALL debug_print_end

! gterm1 = (1 + (rho/g)(dg/d(rho))
! gterm2 = (n/g)(dg/dn_k)

CONTAINS
	SUBROUTINE gcalc_generic(c)
		INTENT(IN)	::	c
		REAL*8 c, g, eta_mix
		eta_mix = DOT_PRODUCT(bvol, comp%x) * rho_mix
		g = 1D0 / (1D0 - c * eta_mix)
		del = KADplus * g ! (ns, ns)
		gterm1 = g ! scalar
		DO k = 1, n
		    ! (ns, ns, n) where (ns x ns) is a copy for each k
			gterm2(:,:,k) = c * bvol(k) * rho_mix / (1D0 - c * eta_mix)
		END DO
	END SUBROUTINE gcalc_generic

	SUBROUTINE gcalc_PCSAFT
		INTEGER l
		REAL*8, DIMENSION(n, n, n)	::	ngdgdnk
		REAL*8, DIMENSION(n, n)		::	g, rhogdgdrho, epsmat, dmat
		REAL*8, DIMENSION(ns, ns) :: dmatsite

		CALL calc_gterms_PCSAFT(g, rhogdgdrho, ngdgdnk, n, T, rho_mix, comp%x, sigma, m, epsok)
		! transform from component matrix to site matrix and then calculate
		CALL comp2site(del, g, n, ns, comp, site) ! del is a temporary placeholder for gij matrix
		CALL comp2site(gterm1, rhogdgdrho, n, ns, comp, site)
		DO k = 1, n
			CALL comp2site(gterm2(:,:,k), ngdgdnk(:,:,k), n, ns, comp, site)
		END DO

		! Calculate dij matrix
		! todo - need to incorporate physical kij
		epsmat(1,:) = epsok ! set first row
		epsmat = spread(epsmat(1,:),1,n) ! copy first row
		epsmat = dsqrt(epsmat*transpose(epsmat)) ! matrix of eps using kij = 0
		dmat(1,:) = sigma ! fill with sigma values
		dmat = spread(dmat(1,:),1,n)
		dmat = (dmat + transpose(dmat))/2D0 ! matrix of average
		dmat = dmat * (1D0 - 1.2D-1 * DEXP(-3D0 * epsmat / T)) ! matrix of dij
		CALL comp2site(dmatsite, dmat, n, ns, comp, site) ! del is a temporary placeholder for gij matrix
        ! pcsaft Delta_ij = N_A*d^3*KADplus*g_ij
		del = 6.0221408D-1*dmatsite**3 * KADplus * del ! del on RHS is gij matrix
		gterm1 = gterm1 + 1D0
	END SUBROUTINE gcalc_PCSAFT

	SUBROUTINE comp2site(sitemat, compmat, n, ns, comp, site)
	! Transform matrix by component to matrix by site, all sites on same
	! component will have same values
		INTENT(OUT)	::	sitemat
		INTENT(IN)	::	compmat, n, ns, comp, site

		INTEGER n, ns, i, j
		REAL*8, DIMENSION(ns, ns)	::	sitemat
		REAL*8, DIMENSION(n, n)		::	compmat
		TYPE(siteinfo), DIMENSION(ns)	::	site
		TYPE(species), DIMENSION(n)	::	comp

		sitemat = 0D0
		DO i = 1, ns
			DO j = 1, ns
				sitemat(i,j) = compmat(site(i)%host, site(j)%host)
			END DO
		END DO
	END SUBROUTINE comp2site

	SUBROUTINE debug_print_start
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'CALL gcalc',REPEAT('-',40)
			PRINT *, '>>>>Input:'
			PRINT *, 'gcalcFlag='; PRINT *, gcalcFlag
			PRINT *, 'n='; PRINT *, n
			PRINT *, 'ns='; PRINT *, ns
			PRINT *, 'T='; PRINT *, T
			PRINT *, 'rho_mix='; PRINT *, rho_mix
			PRINT *, 'rho_pure='; PRINT *, rho_pure
			PRINT *, 'bvol='; PRINT *, bvol
			PRINT *, 'sigma='; PRINT *, sigma
			PRINT *, 'm='; PRINT *, m
			PRINT *, 'epsok='; PRINT *, epsok
			PRINT *, 'KADplus'
			DO i = 1, ns
				PRINT *, KADplus(i,:)
			END DO
			PRINT *, '>>>>end of Input'; PRINT *, ''
		END IF
	END SUBROUTINE debug_print_start

	SUBROUTINE debug_print_end
		IF (debug) THEN
			PRINT *, 'del='
			DO i = 1, ns
				PRINT *, del(i,:)
			END DO
			PRINT *, 'gterm1='
			DO i = 1, ns
				PRINT *, gterm1(i,:)
			END DO
			PRINT *, 'gterm2='
			DO k = 1, n
				PRINT *, k
				DO i = 1, ns
					PRINT *, gterm2(i,:,k)
				END DO
			END DO
			PRINT *, REPEAT('-',40),'ENDOF gcalc',REPEAT('-',40)
		ENDIF
	END SUBROUTINE debug_print_end
END SUBROUTINE gcalc
