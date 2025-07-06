SUBROUTINE calc_dAdnk(dAdnk, kcalc, kop1, gcalcFlag, n, ns, comp, site, T, rho_mix, &
			rho_pure, KADplus, bvol, sigma, m, epsok)

use CONSTANTS
USE sitenspecies, only: siteinfo, species !pass other variables for legacy reasons

IMPLICIT NONE
INTENT(OUT)	::	dAdnk
INTENT(IN)	::	gcalcFlag, n, ns, comp, site, T, rho_mix, rho_pure, &
			KADplus, bvol, sigma, m, epsok

LOGICAL, PARAMETER	::	debug = .False. ! .True. for debugging mode, set .False. to mute
INTEGER gcalcFlag, n, ns, i, j, k, kop1, kcalc
INTEGER, DIMENSION(:), ALLOCATABLE	::	sidk
REAL*8, DIMENSION(ns, ns, n)	::	gterm2
REAL*8, DIMENSION(ns, ns)	::	KADplus, mat, del, gterm1
REAL*8, DIMENSION(ns)		::	vec, Y
REAL*8, DIMENSION(n)		::	dAdnk, rho_pure
REAL*8, DIMENSION(n)	::	bvol, sigma, m, epsok
REAL*8	T, rho_mix, quad_sum_term, sum_term
CHARACTER(LEN=11) FMT
TYPE(species), DIMENSION(n)	::	comp
TYPE(siteinfo), DIMENSION(ns)	::	site

CALL gcalc(del, gterm1, gterm2, gcalcFlag, n, ns, KADplus, &
	comp, site, T, rho_mix, bvol, rho_pure, sigma, m, epsok)
CALL calcX(Y, ns, rho_mix, site%xhost, site%noccur, del)

if((kcalc.eq.1).and.(kop1.gt.1)) then
    if(n.eq.1) then
        write(outfile,*)'Pure component-------------'
    else
        write(outfile,*)'Site hosts (for function call)-----------------'
        write(FMT,'("(",I0,"(I5))")') ns ! make format string for ns columns
        write(outfile,FMT) site(:)%host
    endif
	write(outfile,*) 'Site name'
	write(FMT,'("(",I0,"(A8,2X))")') ns ! make format string for ns columns
    write(outfile,FMT) site(:)%name
    write(outfile,*)'Site fractions, in order of sites'
	write(FMT,'("(",I0,"(F12.5))")') ns ! make format string for ns columns
    write(outfile,FMT) Y
    write(outfile,*)'del matrix by row'
	write(FMT,'("(",I0,"(G12.5))")') ns ! make format string for ns columns
    write(outfile,FMT) ((del(i,j), j=1,ns), i=1,ns)
endif

CALL debug_print_start

vec = Y * site%xhost * site%noccur
dAdnk = 0D0
DO k = 1, n
    sum_term = 0D0
    if (comp(k)%nsite.gt.0)then
        ! sidk is a vector for component k of the site ids on host k
	    ALLOCATE(sidk(comp(k)%nsite))
	    sidk = comp(k)%sid(1:comp(k)%nsite)
	    ! sum term is only for sites on host k
	    sum_term = DOT_PRODUCT(DLOG(Y(sidk)), &
			    site(sidk)%noccur)
	    DEALLOCATE(sidk)
    end if !comp
	mat = del * (rho_mix / rho_pure(k) * gterm1 - gterm2(:,:,k))
	quad_sum_term = DOT_PRODUCT(vec, MATMUL(vec, mat))
	dAdnk(k) = 0.5D0 * rho_mix * quad_sum_term + sum_term
	CALL debug_print_loop
END DO

CALL debug_print_end

CONTAINS
	SUBROUTINE debug_print_start
		INTEGER i
		IF (debug) THEN
			PRINT *, REPEAT('-',40),'CALL calc_dAdnk',REPEAT('-',40)
			PRINT *, '>>>>Input:'
			PRINT *, 'gcalcFlag='; PRINT *, gcalcFlag
			PRINT *, 'n='; PRINT *, n
			PRINT *, 'ns='; PRINT *, ns
			PRINT *, 'T='; PRINT *, T
			PRINT *, 'rho_mix='; PRINT *, rho_mix
			PRINT *, 'rho_pure='; PRINT *, rho_pure
			PRINT *, 'KADplus='
			DO i = 1, ns
				PRINT *, KADplus(i,:)
			END DO
			PRINT *, 'Delta='
			DO i = 1, ns
				PRINT *, Del(i,:)
			END DO
			PRINT *, 'bvol='; PRINT *, bvol
			PRINT *, 'sigma='; PRINT *, sigma
			PRINT *, 'm='; PRINT *, m
			PRINT *, 'epsok='; PRINT *, epsok
			PRINT *, '>>>>end of Input'; PRINT *, ''
			PRINT *, 'site X ='; PRINT *, Y
			PRINT *, 'gterm1 ='; PRINT *, gterm1
			PRINT *, 'gterm2 ='; PRINT *, gterm2
		END IF
	END SUBROUTINE debug_print_start

	SUBROUTINE debug_print_loop
		IF (debug) THEN
			PRINT *, 'k='; PRINT *, k
			PRINT *, 'quad_sum_term='; PRINT *, quad_sum_term
			PRINT *, 'sum_term='; PRINT *, sum_term
		END IF
	END SUBROUTINE debug_print_loop

	SUBROUTINE debug_print_end
		IF (debug) THEN
			PRINT *, 'dAdnk='; PRINT *, dAdnk
			PRINT *, REPEAT('-',40),'ENDOF calc_dAdnk',REPEAT('-',40)
		END IF
	END SUBROUTINE debug_print_end
END SUBROUTINE calc_dAdnk
