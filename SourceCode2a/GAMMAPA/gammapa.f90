! gammapa.f90
!
!*********************************************************
!
! gammapa.f90 - main entry point for application
!
!*********************************************************
! Option codes kop() set by user input on the methods panel for GMUSR
! Debug options
! kop(1) == 0 no debug printing
! kop(1) == 1 for debug printing of only gammas and dll version
! kop(1) == 2 adds printing of X, rho, deltas, etc.
! kop(1) == 3 adds printing of delta matrix indexes, and gamma output to gammas.txt
! rdf equation to use. The calculation uses packing factor based on bvol for CPA, ESD, VDW
! kop(2) == 0 VdW
! kop(2) == 1 ESD
! kop(2) == 2 CPA (ESD with packing factor divided by 4)
! kop(2) == 3 PCSAFT
! Liquid volume
! kop(3) == 0 use standard liquid volume
! kop(3) == 1 use temperature-dependent volume
! residual model
! kop(4) == 0 residual term: NRTL. Parameters should be entered as GMNRTL.
! kop(4) == 1 residual term: Wilson. Parameters should be entered as GMUW.
! kop(4) == 2 residual term: Scatchard Hildebrand. Parameters should be entered as GMUSH.
! kop(4) == 3 residual term: Nagata1. Parameters should be entered as GMUNAG.
! combinatorial correction (ignored when kop(4) = 1 (Wilson), kop(6) == 2 (no comb term))
! kop(5) == 0 combinatorial correction term: vdW with bvdw (calculated by uniquac r*15.17)
! kop(5) == 1 combinatorial correction term: None
! kop(5) == 2 combinatorial correction term: Staverman-Guggenheim
! kop(5) == 3 combinatorial correction term: vdW with the value of bvol
! combinatorial term (ignored for kop(4) == 1 (Wilson))
! kop(6) == 0 combinatorial term: Flory with normal volume
! kop(6) == 1 combinatorial term: Flory term with modified volume (v^2/3)
! kop(6) == 2 combinatorial term: None

PROGRAM GAMMAPA

use constants
! constants includes the outfile id which should not be used for another file
use sitenspecies
use VOLUMES
USE PHYS_PARMS

IMPLICIT NONE

! Please update the version when making more than minor fixes
! Append RTPT or TPT1 so that debug code shows the option used for compiled code.
character(len=10) :: ver = '0.1TPT1'

! todo - track sub and Henry's law components
! INTEGER NSUB,  NSUP

INTEGER KOP(10), KCALC
real*8 T, P
REAL*8, dimension(:), allocatable :: gamma, dgamma

integer i,j,k

! for physical and combinatorial models
! aparam, bparam for NRTL, Nagata, Wilson, tau for NRTL, Aij for SH

print *, 'Welcome to the GAMMAPA program.'
print *, ' '
print *, 'Output will be written to <project>/Output/GAMMAPAout.txt'
print *, ' '
print *, 'You will be prompted for the case-insensitive file names of the two input files.'
print *, 'The files should reside in the folder <project>/Input/GAMMAPA.'
print *, ' '


call loadsites(KOP)

! nc is now known
IF(.NOT.ALLOCATED(x))allocate(x(nc))
IF(.NOT.ALLOCATED(gamma))allocate(gamma(nc))
IF(.NOT.ALLOCATED(dgamma))allocate(dgamma(nc))
IF(.NOT.ALLOCATED(aparam))allocate(aparam(nc,nc))
IF(.NOT.ALLOCATED(bparam))allocate(bparam(nc,nc))
IF(.NOT.ALLOCATED(alpha))allocate(alpha(nc,nc))
IF(.NOT.ALLOCATED(Aij))allocate(Aij(nc,nc))
DO i = 1, nc
    x(i) = 0D0
    DO j = 1, nc
        aparam(i,j) = 0D0
        bparam(i,j) = 0D0
        alpha(i,j) = 0D0
        Aij(i,j) = 0D0
    END DO
END DO

! Prepare for output
OPEN(outfile, file="..\..\Output\GAMMAPAout.txt")
OPEN(debugfile, file="..\..\Output\GAMMAPAdebug.txt")


if(kop(1).gt.0) then
    WRITE(debugfile,"(120('-'))")
	WRITE(debugfile,'(A,A10)') 'Print debug info. GAMMAPA version ',ver
	WRITE(debugfile,'(A,10I3)') 'Option codes ', (kop(i),i=1,10)
    WRITE(debugfile,900) NC
    WRITE(debugfile,*) '*** Species Identifiers'
    do i=1,nc
        WRITE (debugfile, '(3X, I4, 3X, I4, 3X, A20)') i, comp(i)%id, comp(i)%name
    enddo
endif
900  FORMAT('n Species ', I3)

! load physical parameters
select case (kop(4))
    case (0)
        call loadnrtl()
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMNRTL Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMNRTL Physical bij'
            write(debugfile,'(6G12.5)') bparam
            write(debugfile,*) ' *** GMNRTL Physical alpha'
            write(debugfile,'(6G12.5)') alpha
        endif

    case (1) ! wilson
!        call loadwilson(aparam, bparam) !not yet written
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMUwilson Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMUwilson Physical bij'
            write(debugfile,'(6G12.5)') bparam
        endif

    case (2) !scatchard hildebrand (eqn Elliott/Lira 12.24,12.25)
!        call loadsh(aparam, bparam)  !not yet written
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMUscathild Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMUscathild Physical bij'
            write(debugfile,'(6G12.5)') bparam
        endif

	case (3) ! Nagata1
!        call loadnagata(aparam, bparam)  !not yet written
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMUnagata Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMUnagata Physical bij'
            write(debugfile,'(6G12.5)') bparam
            write(debugfile,*) ' *** GMUnagata Physical aij + bij/T'
            write(debugfile,'(6G12.5)') aparam + bparam/T
        endif

end select

IF (kop(2) .EQ. 3) THEN
	if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) then
        do i = 1, nc
            write(debugfile,'(I4, 3G12.5)') (i), PCSAFT_m(i), PCSAFT_sigma(i), PCSAFT_epsok(i)
        enddo
	endif
ENDIF



! Kcalc
! 1 - calculate only gammas
! 2 - calculate only gamma derivative, d(ln gamma)/dT
! 3 - calculate both gamma and gamma derivative

Kcalc = 3
T = 298.15D0 ! K
P = 1D0 ! bar
write(outfile,'(A, I3 )') 'Kcalc ', Kcalc
write(outfile, '(A, 2F10.2)') 'T(K) P(bar) ', T, P
! set composition of interest
! example loop for a binary.
! recall x is shared in sitenspecies
write(outfile,'(A)') 'x, lngamma, gamma, hex'
do i = 1, 11
    x(1) = dble(i-1)/10D0
    x(2) = 1-x(1)
    call gammacalc(kop, Kcalc, T, P, gamma, dgamma)
    write(outfile, '(8F15.6)') X(:), gamma(:), dexp(gamma(:)), -R*T**2*(dot_product(x,dgamma))
enddo

close(outfile)
close(debugfile)

print *, 'Program has ended normally. The output is in Output/GAMMAPAout.txt.'
print *, 'Press any key to close this window.'
pause

if(allocated(gamma)) deallocate(gamma)
if(allocated(dgamma)) deallocate(dgamma)
if(allocated(aparam)) deallocate(aparam)
if(allocated(bparam)) deallocate(bparam)
if(allocated(alpha)) deallocate(alpha)
if(allocated(Aij)) deallocate(Aij)
if(allocated(site)) deallocate(site)
if(allocated(comp)) deallocate(comp)
if(allocated(x)) deallocate(x)
if(allocated(KAD)) deallocate(KAD)
if(allocated(eps)) deallocate(eps)
if(allocated(PCSAFT_sigma)) deallocate(PCSAFT_sigma)
if(allocated(PCSAFT_m)) deallocate(PCSAFT_m)
if(allocated(PCSAFT_epsok)) deallocate(PCSAFT_epsok)
if(allocated(veq)) deallocate(veq)
if(allocated(bvol)) deallocate(bvol)
if(allocated(r_uniquac)) deallocate(r_uniquac)
if(allocated(q_uniquac)) deallocate(q_uniquac)
if(allocated(vstd)) deallocate(vstd)
if(allocated(vparms)) deallocate(vparms)

END PROGRAM GAMMAPA

!****************************************************
SUBROUTINE FINDINDEX(index,compid)
! Find the component index for a given component id

use sitenspecies, only: nc, comp
integer index, compid

do i = 1, nc
 if(comp(i)%id .EQ. compid) EXIT
END DO
index = i
END SUBROUTINE FINDINDEX
!*********************************************************
SUBROUTINE GAMMACALC(kop, Kcalc, T, P, gamma, dgamma)

use CONSTANTS
use sitenspecies
use VOLUMES
use PHYS_PARMS

implicit none

integer, INTENT(IN) :: Kcalc
real*8, INTENT(IN) :: T, P
real*8, dimension(nc), INTENT(OUT) :: gamma, dgamma

INTEGER i, j, k, l, kop(10)
! The intent of KOP is input, but the user for certain models inconsistent values
! can be specified, so the values can be overwritten below.
REAL*8 delT, rhomix, rhomixu, rhomixd
real*8, dimension(nc) :: gammares, dgammares, dgammaresd, gammacomb, &
 dgammacomb, dgammacombd, gammacombcorr, dgammacombcorr, dgammacombcorrd
real*8, dimension(nc) :: gammaw, dgammaw, dgammawd
real*8, dimension(nc) :: bvdw, VU, VD, V, rho, rhou, rhod
real*8, dimension(nc, nc) :: tau, Lambda, Vratio, VratioU, VratioD
! VrationU = volume ratio at T+delT/2
! VrationD = volume ratio at T-deltT/2

delT = 0.1D0

! Set default component molar volumes and size-related covolumes

do i = 1, nc
    V(i) = vstd(i) ! standard volume V in cm3/mol, overwrite later for T-dependent
	VU(i) = V(i) ! upper value used for T-derivative, same as V for standard volume
	VD(i) = V(i) ! down value used for T-derivative, same as V for standard volume
    bvdw(i) = r_uniquac(i)*15.17D0 !vdW covolume parameter
enddo

! these are the ln(gamma) and d(ln(gamma)/dT)
gamma = (/ (0D0, k=1,nc) /)
dgamma = (/ (0.D0, k=1,nc) /)  ! this is d(ln gamma)/dT
gammares = (/ (0D0, k=1,nc) /)
dgammares = (/ (0D0, k=1,nc) /)
gammacomb = (/ (0D0, k=1,nc) /)
dgammacomb = (/ (0D0, k=1,nc) /)
gammacombcorr = (/ (0D0, k=1,nc) /)
dgammacombcorr = (/ (0D0, k=1,nc) /)
gammaw = (/ (0D0, k=1,nc) /)
dgammaw = (/ (0D0, k=1,nc) /)
comp%x = (/ (x(k), k=1,nc) /)

! set site mole fractions
DO k = 1, nsites
    site(k)%xhost = comp(site(k)%host).x
END DO

if((kop(1).gt.1)) THEN
  write(debugfile,*) 'Sites Present: site present in call, site in file, id, name, host, host x'
  do i = 1,nsites
   write(debugfile,'(I3,2X,I4,2X,A20,2X,I3,2X,F10.5)') i, site(i)%id, site(i)%name, site(i)%host, site(i)%xhost
  enddo !i
 write(debugfile,*)'Site1, Site2, index1, index2, Keps, eps, kad'
 do l = 1, nsites
    do k = 1, nsites
        write(debugfile,'(2A20, 2I4, 3G12.5)') site(l)%name,site(k)%name,k,l,kad(k,l)*(dexp(eps(k,l)/T)-1D0), eps(k,l), kad(k,l)
    enddo !k
 enddo !l
endif


! ************** calculate the physical contribution *********************************************

if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) write(debugfile,*) 'Component id, volume and covolume'

! KCALC=1, Calculate property only
! KCALC=2, Calculate temperature derivative of property only
! KCALC=3, Calculate property and temperature derivative.

! overwrite molar volume if kop(3) .gt. 0
if (kop(3) .gt. 0) then
    if (KCALC .gt. 1) then ! calculate values needed for T derivative
        do i = 1,nc
            CALL VLU(VU(i),i,VEQ(i),vparms(:,i),T + delT/2D0)
            CALL VLU(VD(i),i,VEQ(i),vparms(:,i),T - delT/2D0)
        enddo
    endif !kcalc > 1
	if (MOD(KCALC,2).gt.0) then	! calculate molar volume
        do i = 1,nc
            CALL VLU(V(i),i,VEQ(i),vparms(:,i),T)
        enddo
    endif !kcalc is odd
endif

if (kop(1).gt.1) then
! write components present in simulation, molar volume to be used, and user covolume.
    if (mod(kcalc,2).gt.0) then
        write(debugfile,*) ' volume and covolume'
        do i = 1, nc
            write(debugfile,'(I4, 2G12.5)') i, V(i), bvol(i)
        enddo
    endif
    if (kcalc.gt.1) then
        write(debugfile,*) ' vols at T+-deltT/2 for finite diff'
        do i = 1, nc
            write(debugfile,'(I4, 2G12.5)') i, VU(i), VD(i)
        enddo
    endif
endif

select case (kop(4))
    case (0)
		if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call nrtl(dgammares, nc, x, aparam+bparam/(T + delT/2D0), alpha)
            call nrtl(dgammaresd, nc, x, aparam+bparam/(T - delT/2D0), alpha)
            dgammares = (dgammares-dgammaresd)/delT
		endif ! kcalc > 1
		if(MOD(kcalc,2) .gt. 0) then ! kcalc is odd
            tau = aparam+bparam/T
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMNRTL Physical tau'
                write(debugfile,'(6G12.5)') tau
            endif
            call nrtl(gammares, nc, x, tau, alpha)
        endif ! kcalc is odd

    case (1)
        do i=1,nc
            do j = 1,nc
                if (i.EQ.j) then
                    Vratio(i,j)= 1D0
					VratioU(i,j)=1D0
					VratioD(i,j)=1D0
                else
                    Vratio(i,j)= V(j)/V(i)
					! calculate upper and lower values for T-derivative
                    VratioU(i,j) = VU(j)/VU(i)
                    VratioD(i,j) = VD(j)/VD(i)
                endif
            enddo ! j
        enddo ! i
		if(kcalc .gt. 1) then ! todo implement wilson
		    ! calculate upper and lower values for T-derivative
            ! call wilson(dgammares, nc, x, VratioU*dexp(aparam+bparam/(T+delT/2D0)))
            ! call wilson(dgammaresd, nc, x, VratioD*dexp(aparam+bparam/(T-delT/2D0)))
            dgammares=(dgammares-dgammaresd)/delT
        endif ! kcalc > 1
        if( MOD(kcalc,2) .gt. 0) then ! kcalc is odd
            Lambda = Vratio*dexp(aparam+bparam/T)
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUwilson Physical Lambda'
                write(debugfile,'(6G12.5)') Lambda
            endif
            ! call wilson(gammares, n, x, Lambda)
        endif ! kcalc is odd

    case (2) !scatchard hildebrand (eqn Elliott/Lira 12.24,12.25) todo implement SH
        if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            ! call scathild(dgammares, x, VU, nc, aparam*(T + delT/2D0)+bparam,R, T + delT/2D0)
            ! call scathild(dgammaresd, x, VD, nc, aparam*(T - delT/2D0)+bparam,R, T - delT/2D0)
            dgammares = (dgammares - dgammaresd)/delT
        endif ! kcalc < 1
        if (MOD(kcalc,2) .gt. 0) then !kcalc is odd
            Aij = aparam*T+bparam
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUscathild Physical AIJ'
                write(debugfile,'(6G12.5)') AIJ
            endif
            ! call scathild(gammares, x, V, nc, Aij, R, T)
        endif ! kcalc is odd

	case (3) ! Nagata1 todo implement Nagata
        if (kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            ! call nagata1(dgammares, nc, x, VU, aparam + bparam/(T + delT/2D0))
            ! call nagata1(dgammaresd, nc, x, VD, aparam + bparam/(T - delT/2D0))
            dgammares = (dgammares - dgammaresd)/delT
			WRITE(debugfile,*) 'Calculating temperature derivative of residual'
        endif
        if (MOD(kcalc,2).gt.0) then !kcalc is odd
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUnagata Physical aij + bij/T'
                write(debugfile,'(6G12.5)') aparam + bparam/T
            endif
            ! call nagata1(gammares, n, x, r_uniquac, aparam)
            ! call nagata1(gammares, n, x, V, aparam + bparam/T)
        endif

end select
! ********************************** end of gamma phyiscal ********************

! ********************************** Wertheim contribution ***********************

if (aspmx.eq.1) then ! only execute if association exists

	if (mod(kcalc,2) .gt. 0) then ! odd kcalc, calculate properties
		rho = 1D0/V ! molar density mol/cc
		rhomix = 1D0/(dot_product(x,V)) ! mol/cc

		CALL calc_gammaw(gammaw, kcalc, kop(1), kop(2), nc, nsites, comp, site, &
							   T, rhomix, rho, KAD, eps, bvol, &
							   PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
	endif

	if (kcalc .gt. 1) then ! need temperature derivative
		rhoU = 1D0/VU ! molar density mol/cc
		rhomixU = 1D0/(dot_product(x,VU)) ! mol/cc
		rhoD = 1D0/VD ! molar density mol/cc
		rhomixD = 1D0/(dot_product(x,VD)) ! mol/cc

		CALL calc_gammaw(dgammaw,  kcalc, kop(1),kop(2), nc, nsites, comp, site, &
							   (T + delT / 2D0), rhomixU, rhoU, KAD, eps, &
							   bvol, PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
		CALL calc_gammaw(dgammawd, kcalc, kop(1), kop(2), nc, nsites, comp, site, &
							   (T - delT / 2D0), rhomixD, rhoD, KAD, eps, &
							   bvol, PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
		dgammaw = (dgammaw - dgammawd) / delT
	endif
endif !aspmx eq 1

! ------------ end of section if association occurs in mixture-----------

! ************************ end of Werthiem contribution ********************************

! combinatorial term always uses volume, not r
if (kop(4) .ne. 1) then  ! if not Wilson
    if(kcalc.gt.1) then
	if (kop(6) .eq. 0) then
		call flory(dgammacomb, x, VU, nc)
		call flory(dgammacombd, x, VD, nc)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	elseif (kop(6) .eq. 1) then
		call flory(dgammacomb, x, VU**(2D0/3D0), nc)
		call flory(dgammacombd, x, VD**(2D0/3D0), nc)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	endif
    endif !kcalc > 1
    if(mod(kcalc,2).gt.0) then
	if (kop(6) .eq. 0) then
		call flory(gammacomb, x, V, nc)
        if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Method: Flory using molar volume'
	elseif (kop(6) .eq. 1) then
		call flory(gammacomb, x, V**(2D0/3D0), nc)
        if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Method: Flory using molar volume^(2/3)'
	endif
	endif
endif

! **********Combinatorial Correction
if ((kop(4) .eq. 1) .or. (kop(6) .eq. 2)) then
    kop(5) = 1; ! force combinatorial correction off if using Wilson, or if combinatorial term is off
    ! note the gammacombcorr is set to zero at the beginning of the routine
    if(kop(1).gt.1) WRITE(debugfile,*) '*** Your selection of Comb Correction May Be Inconsistent with Other Settings'
    if(kop(1).gt.1) WRITE(debugfile,*) '*** Combinatorial Correction Method: None'
endif
if (kop(5) .ne. 1) then
    if(kcalc.gt.1) then
        select case(kop(5))
        case (0)
            call correqn1262(dgammacombcorr, x, VU, nc, bvdw)
            call correqn1262(dgammacombcorrd, x, VD, nc, bvdw)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        case (2)
            dgammacombcorr = 0
        case (3)
	        call correqn1262(dgammacombcorr, x, VU, nc, bvol)
            call correqn1262(dgammacombcorrd, x, VD, nc, bvol)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        end select
    endif
    if(mod(kcalc,2).gt.0) then
        select case(kop(5))
        case (0)
            call correqn1262(gammacombcorr, x, V, nc, bvdw)
            if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: vdW using bvdw'
        case (2)
            call stavgugg(gammacombcorr, x, r_uniquac,q_uniquac, nc)
            if (kop(1) .gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: Staverman-Guggenheim usig r and q'
        case (3)
	        call correqn1262(gammacombcorr, x, V, nc, bvol)
            if (kop(1) .gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: vdW using COVOL'
        end select
    endif
end if

! add logs of physical and wertheim contributions
gamma = gammares + gammaw + gammacomb + gammacombcorr
dgamma = dgammares + dgammaw + dgammacomb + dgammacombcorr

if(kop(1) .gt. 0) then ! write to history and/or .csv files
	  WRITE(debugfile,*) '*** GAMMA RESULTS'
	  WRITE(debugfile,900) KCALC, nc, T, P
 900  FORMAT(' KCALC, N, T, P =', I3, 1X, I3, F10.3, 2X, E13.6)
	  if(mod(kcalc,2).gt.0) then !kcalc is odd, write gammas to history
		WRITE(debugfile,1010)
		DO 400 I = 1, nc! write gamma contributions to history
		 WRITE(debugfile,1000) i, comp(i)%name, x(i), gammares(i),gammacomb(i),gammacombcorr(i), gammaw(i), GAMMA(I), dgamma(i)
 400   CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
	   if(kop(1) .gt. 2) then
		! write to vol.csv
		open(unit=12, file="vol.csv", status='unknown', action='write', position='append')
		WRITE(12,1002) T, ',', 1D0/rhomix, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", V(i), i=1,nc)
		close(12)
		! write to gammas.csv
		open(unit=10, file="gammas.csv", status='unknown', action='write', position='append')
		WRITE(10,1001) T, ',', P, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", gammares(i),",",gammacomb(i),",",gammacombcorr(i), ",",gammaw(i), ",",GAMMA(I), ",",dgamma(i), i=1,nc)
		close(10)
	   endif !kop(1).gt.2
	  endif !mod(kcalc,2)
	  if(kcalc.gt.1) then !write derivative
		WRITE(debugfile,1011)
		DO 410 I = 1, nc
		 WRITE(debugfile,1000) i, comp(i)%name, X(I), dgammares(i),dgammacomb(i),dgammacombcorr(i), dgammaw(i), GAMMA(I), dgamma(i)
 410    CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
		if(kop(1).gt.2) then
		! write to HXS.csv
		open(unit=66, file="HXS.csv", status='unknown', action='write', position='append')
		WRITE(66,1001) T, ",", P, (', ', i, comp(i)%name, ",", X(I), ",", -dgammares(i)*R*(T**2), ",", -dgammacomb(i)*R*(T**2), ",", -dgammacombcorr(i)*R*(T**2), ",", -dgammaw(i)*R*(T**2), ",", -dgamma(i)*R*(T**2), ',', dgamma(i), i=1,nc)
		close(66)
		endif !if kcalc.gt.1
		endif !kop(1).gt.2

		WRITE(debugfile,*) 'End of function call'
	  WRITE(debugfile,"(120('*'))")
endif !kop(1).gt.0

!
!
!
!     Format statements
!
 1000   FORMAT (I3,1X, A8, 2X, 7F15.6)
 1001	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, A20, 7(A1,F15.6))) ! format for gamma, hex csv files, up to 100 components
 1002	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, A20, 2(A1,F15.6))) ! format for vol csv files, up to 100 components
 1010   FORMAT (' I     NAME          X         lnGAMMAres    lnGAMMAcomb  lnGAMMAcombcorr       lnGAMMAw       lnGAMMA      dlnGAMMA')
 1011   FORMAT (' I     NAME          X         dlnGAMMAres   dlnGAMMAcomb dlnGAMMAcombcorr     dlnGAMMAw       lnGAMMA      dlnGAMMA')
!
return
end subroutine GAMMACALC
