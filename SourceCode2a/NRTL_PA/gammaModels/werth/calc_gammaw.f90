SUBROUTINE calc_gammaw(gammaw, kcalc, kop1, gcalcFlag, n, ns, comp, site, T, &
			rho_mix, rho_pure, KAD, epsADok, bvol, sigma, m, epsok)
!Description here...
! Output:
!	- gammaw	(n)		natural log of gamma from association contribution
!   - Y         (ns)     fraction free
! Input:
!	- gcalcFlag 	(scalar)	select method for calculating gterms
!						0. VdW
!						1. ESD
!						2. CPA
!						3. PCSAFT
!   - n		    (scalar)	number of components present
!	- ns 		(scalar)	numer of sites present
!	- T		    (scalar)	temperature [K]
!	- rho_mix	(scalar)	density of mixture [mol/cc]
!	- rho_pure	(n)		density of pure component [mol/cc]
!	- comp		(n)		component info array of sites present
!	- bvol		(n)		covol [??]
!	- sigma		(n)		(gcalcFlag=3)
!	- m		    (n)		(gcalcFlag=3)
!	- epsok		(n)		(gcalcFlag=3)
!	- site		(ns)		site info array for sites present
!	- KAD		(ns,ns)		[??]
!	- epsADok	(ns,ns)		[K]
!
! Parameter:
!	- KADplus	(ns,ns)		[??]

USE sitenspecies
IMPLICIT NONE

INTENT(OUT)	::	gammaw
INTENT(IN)	::	 kcalc, kop1,gcalcFlag, n, ns, comp, site, T, rho_mix, rho_pure, &
			KAD, epsADok, bvol, sigma, m, epsok

LOGICAL, PARAMETER	::	debug = .False. ! .True. for debugging mode, set .False. to mute
INTEGER  kcalc, kop1,gcalcFlag, n, ns
REAL*8 T, rho_mix
REAL*8, DIMENSION(n)		::	gammaw, rho_pure, dAdnk, dAdnk_pure
REAL*8, DIMENSION(ns, ns)	::	Y,KAD, epsADok, KADplus
TYPE(species), DIMENSION(n)	::	comp
TYPE(siteinfo), DIMENSION(ns)	::	site
!REAL*8, OPTIONAL, DIMENSION(n)	::	bvol, sigma, m, epsok
REAL*8, DIMENSION(n)	::	bvol, sigma, m, epsok

CALL debug_print_start

KADplus = KAD * (DEXP(epsADok / T) - 1D0)
CALL calc_dAdnk(dAdnk, kcalc, kop1, gcalcFlag, n, ns, comp, site, T, rho_mix, &
			rho_pure, KADplus, bvol, sigma, m, epsok)
!
CALL calc_dAdnk_pure
gammaw = dAdnk - dAdnk_pure

CALL debug_print_end

CONTAINS
	SUBROUTINE calc_dAdnk_pure
		! Loop through components and call calc_dAdnk to calculate dAdnk for
		! self associating component.
		INTEGER i, j, k, ns_pure

		TYPE(species), DIMENSION(1)	::	comp_pure
		TYPE(siteinfo), DIMENSION(:), ALLOCATABLE	::	site_pure
		REAL*8, DIMENSION(:,:), ALLOCATABLE	::	KADplus_pure

		dAdnk_pure = 0D0 ! default value
		DO k = 1,n
			IF (comp(k)%selfassoc) THEN
				ns_pure = comp(k)%nsite
				! create variables for pure component
				ALLOCATE(site_pure(ns_pure))
				ALLOCATE(KADplus_pure(ns_pure, ns_pure))

				site_pure = site(comp(k)%sid(1:comp(k)%nsite)) ! copy all site info for component k
				site_pure%host = 1 ! There is only one dummy component structure
				site_pure%xhost = 1D0 ! set composition to purity
				comp_pure(1) = comp(k) ! copy all of component information
				comp_pure(1)%x = 1D0 ! overwrite composition to purity
				! Becasue the KAD matrix will be full, all sites will
				! be used and filled in the order. Thus set the sid to
				! successive integers.
				DO i = 1, ns_pure
					comp_pure(1)%sid(i) = i
				END DO
                ! Fill KADplus matrix with only the sites on the host k
				DO i = 1, ns_pure
					DO j = 1, ns_pure
						KADplus_pure(i,j) = KADplus(comp(k)%sid(i), comp(k)%sid(j))
					END DO
				END DO

				CALL calc_dAdnk(dAdnk_pure(k:k), kcalc, kop1, gcalcFlag, 1, ns_pure, &
					comp_pure, site_pure, T, rho_pure(k), rho_pure(k:k), &
					KADplus_pure, bvol(k:k), sigma(k:k), m(k:k), epsok(k:k))
				DEALLOCATE(site_pure, KADplus_pure)
			END IF
		END DO
	END SUBROUTINE calc_dAdnk_pure

	SUBROUTINE debug_print_start
		INTEGER i
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'CALL calc_gammaw',REPEAT('-',40)
			PRINT *, '>>>>Input:'
			PRINT *, 'gcalcFlag='; PRINT *, gcalcFlag
			PRINT *, 'n='; PRINT *, n
			PRINT *, 'ns='; PRINT *, ns
			PRINT *, 'T='; PRINT *, T
            PRINT *, 'x='; PRINT *, comp%x
			PRINT *, 'rho_mix='; PRINT *, rho_mix
			PRINT *, 'rho_pure='; PRINT *, rho_pure
			PRINT *, 'KAD='
			DO i = 1, ns
				PRINT *, KAD(i,:)
			END DO
			PRINT *, 'epsADok='
			DO i = 1, ns
				PRINT *, epsADok(i,:)
			END DO
			PRINT *, 'bvol='; PRINT *, bvol
			PRINT *, 'sigma='; PRINT *, sigma
			PRINT *, 'm='; PRINT *, m
			PRINT *, 'epsok='; PRINT *, epsok
			PRINT *, '>>>>end of Input'; PRINT *, ''
		END IF
	END SUBROUTINE debug_print_start

	SUBROUTINE debug_print_end
		INTEGER i
		IF (debug) THEN
			PRINT *, 'KADplus='
			DO i = 1, ns
				PRINT *, KADplus(i,:)
			END DO
			PRINT *, 'dAdnk='; PRINT *, dAdnk
			PRINT *, 'dAdnk_pure='; PRINT *, dAdnk_pure
            PRINT *, 'GAMMAW='; PRINT *, gammaw
			PRINT *, REPEAT('-',40),'ENDOF calc_gammaw',REPEAT('-',40)
		END IF
	END SUBROUTINE debug_print_end
END SUBROUTINE calc_gammaw
