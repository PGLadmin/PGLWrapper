! Aspen simulation indexes refer to the order components or groups are entered into the GUI forms.
! Not all components or groups are present in a function call.
! The n components present in a function call are identified by the aspen simuluation indexes loaded
! into idx. If simulation components 2, 4, 5 are present in a function call, n = 3, idx = [2,4,5].
! Similarly, not all simulation groups will be present in a function call. User UNIFAC groups
! are used to load group information. Only user groups (id > 4500) are used for Wertheim sites.
! The number of sites present in a function call, nsitesp, is found by looping over n,
! and finding the user groups that are present in those components.

! When switching between public and local version, search for kop(4).
! For public version remove comment status where it is forced to zero.

SUBROUTINE GMU (T,       P,      X,      N,      IDX, &
     &                IRW,     IIW,    KCALC,  KOP,    NDS, &
     &                KDIAG,   GAMMA,  DGAMMA, KER,    NSUB, &
     &                NSUP,    IDXSUB, IDXSUP, GAMSC,  DGAMSC)
use findsites
use sitenspecies
IMPLICIT NONE

! The code is debugged using the procedure in Aspen Solution ID 122317.
! We are currently compiling with Intel Fortran 10.1 and Visual Studio 2005.
!
!     DECLARE ARGUMENTS
!
! ver is the software version
! To find version number of dll, turn debugging printing on and search history file for 'GMU version'
! 0.2 adding debug printing of intermediate calculations
! 0.5 fix bug in finding packed site index for pure, improved debugging printing for kop(1) = 2
! 0.7 add code to permit interactive debugging running regression
! 0.8 reinsert calculation of V as a function of T, kop(3) = 1
! 0.9 adding combinatorial term (Flory+vdW correction)
! 1.0 adding staverman-guggenheim correction option and binary (non-generalized) scatchard-hildebrand
! 1.1 adding Nagata1 and modified Flory2-3rd. Adding option for vdW combinatorial correction using COVOL.
! 2.0 adding temperature derivative
! 2.01 fix label on gamma
! 2.1 adding gammas.txt output and CPA parameter option, printing gammas instead of ln(gamma)s.
!     for CPA parameter input and combinatorial correction. CPA input will now use b instead of b/4.
! 2.2 fix hmix for wilson
! 2.3 output HXS.txt when in verbose debug mode
! 2.4 changing GMUC to GMUNRTL and GMUNT to GMUNAG and a couple of pointer names
! 2.5 modify to use external subroutine loadsites in external module
! 2.51 GMUNRTL is too long, changing to GMNRTL
! 2.60 add PCSAFT calculation, move gamma werth calculation to external subroutine
! 2.6.2 add user volume calculation
! 2.6.3 fix loadsites refference to idx and added debug code to calc_dAdnk.f90
! 2.6.4 print PCSAFT parameters
! 2.6.5 version for distribution. Remove the factor of 2 on symmetrical association constants
! 2.6.6 add capability for KADN and EPSN for CPA, ESD, VDW
! 2.6.7 modify loadsites to save Ksave, Esave
! revise 'ver' when making changes to code.

! character(len=10) :: ver = '2.6.7TPT1'
! Append RTPT or TPT1 so that debug code shows the option used for compiled code.
character(len=10) :: ver = '2.6.7TPT1'

INTEGER n, NSUB,  NSUP
!
#include "ppexec_user.cmn"
#include "dms_maxwrt.cmn" ! for writing to control panel

      EQUIVALENCE (RMISS, USER_RUMISS)
      EQUIVALENCE (IMISS, USER_IUMISS)
!
#include "pputl_ppglob.cmn"
#include "dms_global.cmn"   !this includes the history file id.
!
!
#include "dms_ppctbl.cmn"
#include "dms_ppwork.cmn"
#include "dms_ipwork.cmn"
!
#include "dms_plex.cmn"
        REAL*8 B(1)
        EQUIVALENCE (B(1), IB(1))

#include "dms_errout.cmn"
#include "dms_ncomp.cmn"


      INTEGER IDX(N),KOP(10),      IDXSUB(NSUB), &
             IDXSUP(NSUP), IRW,   IIW,   KCALC, NDS, &
             KDIAG, KER
      REAL*8 X(N),  xsav(n), GAMMA(N),     DGAMMA(N), &
            GAMSC(NSUB,NSUP),    DGAMSC(NSUB,NSUP), &
            T,     P
	  real*8 delT; ! increment for numerical derivatives.
      REAL*8 V(n), VU(n), VD(n), PHI(n), H(n),ENTR(n),G(n),DPHI(n),DH(n),DS(n),DG(n),DV(n),RMISS
	  ! VU = V at T + delT/2
	  ! VD = V at T - delT/2
      INTEGER NBOPST(6), NAME(2),IMISS
      ! next three lines used for PPMON call
      INTEGER, PARAMETER:: KV = 1 !This is a calculation code that tells the subroutine to calculate only the property and not its temperatue derivative
      INTEGER, PARAMETER:: KPHI = 0,KH=0,KS=0,KG=0 !This is a calculation code that tells the subroutine to NOT calculate the property nor its temperatue derivative
      INTEGER, PARAMETER:: KBASE = 1 !Thermodynamic function base code.elements in their standard states at 298.15 K.
      INTEGER DMS_IFCMNC, DMS_LOCATI, LGAMPAR, LCOV, LUFGRP, LUBMUW, LUGMUW, LNDS, LNDSW, LNDSCOV, NSITE, LUGRP, LPCSFTM, LPCSFTV, LPCSFTU

!     DECLARE LOCAL VARIABLES
INTEGER lgrp, lgrpe, JGRP, kcalcs
! aspmx =1 if association or solvation is present in mixture (or pure components)
! nsites is the number of sites in the aspen simulation file, nsitesp is the number of sites present in the function call based on components.
integer nsites, nsitesp
integer i,j,k,l
real*8 bmix, eta, etamix,etamixU, etamixD, R, Tsav
character(len=8) name_grp1, name_grp2
integer name1, name2
integer iset, idatij, nel, nelem, loff, LBGV, VEQN(n)
! to store plex offsets
integer LBPROC, NBG, LBG, NG, loffk, loffe, lvlstd, LIDSCC, LURVDW,LUNIQQ, LGMUSH, LTHRSWT,LDNLDIP
real*8 gammaw(n), dgammaw(n), dgammawd(n), rho(n), rhoU(n), rhoD(n), bvol(n), bvdw(n), UNIQQ(n), UNIQR(n), DNLDIP(n,7)
real*8 rhomix, rhomixU, rhomixD, gs, gsU, gsD, gspure, gspureU, gspureD
real*8 hmix, hmixU, hmixD, hpure, hpureU, hpureD, term1, term1U, term1D, gterm, gtermU, gtermD !gs is rdf at contact (sigma)
real*8, dimension(n) :: PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok ! PCSAFT parameters
real*8, dimension(:,:), allocatable :: KAD, eps
real*8, dimension(:), allocatable :: Ymix, YmixU, YmixD, Ypure, YpureU, YpureD
integer, dimension(:), allocatable :: nspure, ids_pack
character*256 buffer
! for physical model
! added Lamda for wilson
real*8 gammares(n),dgammares(n), dgammaresd(n), gammacomb(n), dgammacomb(n), dgammacombd(n), gammacombcorr(n), dgammacombcorr(n),dgammacombcorrd(n)
real*8 aparam(N,N), bparam(N,N), alpha(n,n), tau(n,n), Aij(n,n), Lambda(n,n),Vratio(n,n), VratioU(n,n), VratioD(n,n)
! VrationU = volume ratio at T+delT/2
! VrationD = volume ratio at T-deltT/2
! for debugging
! integer nametemp(80)
! real*8 temp(50)
INTEGER IENTER
DATA IENTER /0/

type (siteinfo), dimension(:), allocatable :: site, sitep ! dimensioned later to be the number of sites, sites present
type (species) :: compp(n) ! components present

! Uncomment this section to force pause to attach debugger when using regression or aspen properties.
! To run a properties or regression, view the case in the properties environment. Then click 'input'.
! Use 'Save as..' and save as a .inp file, e.g. 'myfile.inp'. Run the file as an aspen run e.g.
! aspen myfile /dlopt="dlopt.txt"
! When command line debugging is complete or when you need to run from GUI agin, do not forget
! to comment these lines and recompile, or the run will freeze waiting for input in the background,
! and you won't see the prompt from the GUI.
!  IF (ienter .EQ. 0) THEN
!  WRITE(*,'(A)') 'Enter a number'
!  READ (*,*) ienter
!  ienter = 1
!  ENDIF

delT = 0.1D0

! these are the ln(gamma) and d(ln(gamma)/dT)
gamma = (/ (0D0, k=1,n) /)
dgamma = (/ (0.D0, k=1,n) /)  ! this is d(ln gamma)/dT
gammares = (/ (0D0, k=1,n) /)
dgammares = (/ (0D0, k=1,n) /)
gammacomb = (/ (0D0, k=1,n) /)
dgammacomb = (/ (0D0, k=1,n) /)
gammacombcorr = (/ (0D0, k=1,n) /)
dgammacombcorr = (/ (0D0, k=1,n) /)
gammaw = (/ (0D0, k=1,n) /)
dgammaw = (/ (0D0, k=1,n) /)
R=PPGLOB_RGAS/1000D0 !universal gas constant converted to J/mol-K

if (n .eq. 1) return
!     RETRIEVE THE GMU PARAMETER AREA OFFSETS
lvlstd = DMS_IFCMNC('VLSTD') ! standard liquid volume
LURVDW = DMS_IFCMNC('GMUQR') ! UNIFAC R group
LUNIQQ = DMS_IFCMNC('GMUQQ') ! UNIFAC Q group
LCOV = DMS_IFCMNC('COVOL') ! covolume
LUFGRP = DMS_IFCMNC('UFGRP') ! UNIFAC group used for site information
LIDSCC = DMS_IFCMNC('IDSCC') !component name
LDNLDIP = DMS_IFCMNC('DNLDIP') ! DIPPR density constants
LTHRSWT = DMS_IFCMNC('THRSWT') ! Thermo switch
IF (PPCTBL_NBGV .GT. 0) LBGV = DMS_LOCATI(PPCTBL_NBGV)

if(kop(1).gt.0) then
    WRITE(GLOBAL_NH,"(120('-'))")
	WRITE(GLOBAL_NH,'(A57 ,A10)') 'GMU - Print debug info. Start of func call. GMU version ',ver
	WRITE(GLOBAL_NH,'(A20,10I3)') 'GMUSR Option codes ', (kop(i),i=1,10)
    WRITE(GLOBAL_NH,900) KCALC, N, T, P
    WRITE(global_nh,*) ' *** Aspen Component Identifiers'
    do i=1,n
        j = LIDSCC + 2*(idx(i)-1)
        WRITE (global_nh, '(3X, I4, 3X, 2A4)') idx(i),(IB(j+k), k=1,2)
    enddo
endif

IF (kop(2) .EQ. 3) THEN
	LPCSFTM = DMS_IFCMNC('PCSFTM') ! PCSAFT m
	LPCSFTV = DMS_IFCMNC('PCSFTV') ! PCSAFT sigma
	LPCSFTU = DMS_IFCMNC('PCSFTU') ! PCSAFT epsok

	if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) write(global_nh,*) 'Component id, m, sigma, epsok'

	DO i = 1, N
		PCSAFT_m(i) = B(LPCSFTM + idx(i))
		PCSAFT_sigma(i) = B(LPCSFTV + idx(i))
		PCSAFT_epsok(i) = B(LPCSFTU + idx(i))
	ENDDO

	if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) then
        do i = 1, n
            write(global_nh,'(I4, 3G12.5)') idx(i), PCSAFT_m(i), PCSAFT_sigma(i), PCSAFT_epsok(i)
        enddo
	endif
ENDIF

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

! NDS is passed from the main function call and may not apply to Wertheim groups.
! NEL = 3 ! specifically for NTRL
! LNDS = LGAMPAR+NEL*NCOMP_NCC*NCOMP_NCC*(NDS-1)
! LNDSW = LGAMPAR
! LNDSCOV = LCOV

!todo
! COVOL is needed only for CPA, ESD, VDW g(sigma) and/or combinatorial correction. Code can be cleaned
! so that is not required when running PC-SAFT.

! ************** calculate the physical contribution *********************************************
! For public version force kop(4)=0
! kop(4) = 0

! temp = B(lvlstd+1:lvlstd+1+N)*1000 ! V in cm3/mol - for debugging
if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) write(global_nh,*) 'Component id, volume and covolume'

! Find component molar volumes and size-related covolumes

do i = 1, n
    V(i) = B(lvlstd+idx(i))*1000D0 ! standard volume V in cm3/mol, overwrite later for T-dependent
	VU(i) = V(i) ! upper value used for T-derivative, same as V for standard volume
	VD(i) = V(i) ! down value used for T-derivative, same as V for standard volume
    UNIQR(i)= B(LURVDW+idx(i)) ! UNIQUAC value for volume
    UNIQQ(i)= B(LUNIQQ+idx(i)) ! UNIQUAC value for surface area
    bvol(i) = B(LCOV+idx(i)) ! user covolume (cm3/mol)
    bvdw(i) = B(LURVDW+idx(i))*15.17D0 !vdW covolume parameter
enddo

! KCALC=1, Calculate property only
! KCALC=2, Calculate temperature derivative of property only
! KCALC=3, Calculate property and temperature derivative.

! overwrite molar volume if kop(3) .gt. 0
if (kop(3) .gt. 0) then
    xsav = x  ! calling the pure property calculator clobbers composition, save composition
    kcalcs = kcalc ! calling the pure property calculator clobbers kcalc, save setting
    do i = 1,n
        VEQN(i) = B(LTHRSWT + 8*idx(i) - 6) ! look up density equation
        DNLDIP(i,:) = B(LDNLDIP + 7*idx(i) - 6:LDNLDIP + 7*idx(i))
    enddo
    if (KCALC .gt. 1) then ! calculate values needed for T derivative
        do i = 1,n
            CALL VLU(VU(i),idx(i),VEQN(i),DNLDIP(i,:),T + delT/2D0)
            CALL VLU(VD(i),idx(i),VEQN(i),DNLDIP(i,:),T - delT/2D0)
        enddo
    endif !kcalc > 1
	if (MOD(KCALC,2).gt.0) then	! calculate molar volume
        do i = 1,n
            CALL VLU(V(i),idx(i),VEQN(i),DNLDIP(i,:),T)
        enddo
    endif !kcalc is odd
endif

if (kop(1).gt.1) then
! write components present in simulation, molar volume to be used, and user covolume.
    if (mod(kcalc,2).gt.0) then
        write(global_nh,*) ' volume and covolume'
        do i = 1, n
            write(global_nh,'(I4, 2G12.5)') idx(i), V(i), bvol(i)
        enddo
    endif
    if (kcalc.gt.1) then
        write(global_nh,*) ' vols at T+-deltT/2 for finite diff'
        do i = 1, n
            write(global_nh,'(I4, 2G12.5)') idx(i), VU(i), VD(i)
        enddo
    endif
endif

select case (kop(4))
    case (0)
        LGAMPAR=DMS_IFCMNC('GMNRTL')  ! NRTL binary parameters
        NEL = 3 ! specifically for NTRL
        LNDS = LGAMPAR+NEL*NCOMP_NCC*NCOMP_NCC*(NDS-1)
        LNDSCOV = LCOV

        ! when tau = aij + bij/T, if aij is first element, IEL = 1, for bij, IEL = 2.
        do i=1,n
            do j = 1,n
                if (i.EQ.j) then
                    aparam (i,j)= 0D0
                    bparam(i,j)= 0D0
                    alpha(i,j) = 0D0
                else
                    L = LNDS+NEL*NCOMP_NCC*(IDX(J) -1)+NEL*(IDX(I)-1)
                    K = L+1 ! get element number 1, see pg 162 of v8.8 user guide
                    aparam (i,j)= B(K)
                    bparam(i,j)= B(K+1)
                    alpha(i,j) = B(K+2)
                endif
            enddo ! j
        enddo ! i
		if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call nrtl(dgammares, n, x, aparam+bparam/(T + delT/2D0), alpha)
            call nrtl(dgammaresd, n, x, aparam+bparam/(T - delT/2D0), alpha)
            dgammares = (dgammares-dgammaresd)/delT
		endif ! kcalc > 1
		if(MOD(kcalc,2) .gt. 0) then ! kcalc is odd
            tau = aparam+bparam/T
            if(kop(1).gt.1) then
                write(global_nh,*) ' *** GMNRTL Physical aij'
                write(global_nh,'(6G12.5)') aparam
                write(global_nh,*) ' *** GMNRTL Physical bij'
                write(global_nh,'(6G12.5)') bparam
                write(global_nh,*) ' *** GMNRTL Physical alpha'
                write(global_nh,'(6G12.5)') alpha
                write(global_nh,*) ' *** GMNRTL Physical tau'
                write(global_nh,'(6G12.5)') tau
            endif
            call nrtl(gammares, n, x, tau, alpha)
        endif ! kcalc is odd

    case (1)
        LGAMPAR=DMS_IFCMNC('GMUW')  ! wilson binary parameters
        NEL = 2 ! specifically for wilson
        LNDS = LGAMPAR+NEL*NCOMP_NCC*NCOMP_NCC*(NDS-1)
        LNDSCOV = LCOV

        ! when Lambda = aij + bij/T, if aij is first element, IEL = 1, for bij, IEL = 2.
        do i=1,n
            do j = 1,n
                if (i.EQ.j) then
                    aparam (i,j)= 0D0
                    bparam(i,j)= 0D0
                    Vratio(i,j)= 1D0
					VratioU(i,j)=1D0
					VratioD(i,j)=1D0
                else
                    L = LNDS+NEL*NCOMP_NCC*(IDX(J) -1)+NEL*(IDX(I)-1)
                    K = L+1 ! get element number 1, see pg 162 of v8.8 user guide
                    aparam (i,j)= B(K)
                    bparam(i,j)= B(K+1)
                    Vratio(i,j)= V(j)/V(i)
					! calculate upper and lower values for T-derivative
                    VratioU(i,j) = VU(j)/VU(i)
                    VratioD(i,j) = VD(j)/VD(i)
                endif
            enddo ! j
        enddo ! i
		if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call wilson(dgammares, n, x, VratioU*dexp(aparam+bparam/(T+delT/2D0)))
            call wilson(dgammaresd, n, x, VratioD*dexp(aparam+bparam/(T-delT/2D0)))
            dgammares=(dgammares-dgammaresd)/delT
        endif ! kcalc > 1
        if( MOD(kcalc,2) .gt. 0) then ! kcalc is odd
        Lambda = Vratio*dexp(aparam+bparam/T)
        if(kop(1).gt.1) then
            write(global_nh,*) ' *** GMUwilson Physical aij'
            write(global_nh,'(6G12.5)') aparam
            write(global_nh,*) ' *** GMUwilson Physical bij'
            write(global_nh,'(6G12.5)') bparam
            write(global_nh,*) ' *** GMUwilson Physical Lambda'
            write(global_nh,'(6G12.5)') Lambda
        endif
        call wilson(gammares, n, x, Lambda)
        endif ! kcalc is odd

    case (2) !scatchard hildebrand (eqn Elliott/Lira 12.24,12.25)
        LGAMPAR=DMS_IFCMNC('GMUSH')  ! scatchard-hildebrand binary parameters
        NEL = 2 ! specifically for scatchard-hildebrand
        LNDS = LGAMPAR+NEL*NCOMP_NCC*NCOMP_NCC*(NDS-1)
        LNDSCOV = LCOV

        do i=1,n
            do j = 1,n
                if (i.EQ.j) then
                    aparam (i,j)= 0D0
                    bparam(i,j)= 0D0
                else
                    L = LNDS+NEL*NCOMP_NCC*(IDX(J) -1)+NEL*(IDX(I)-1)
                    K = L+1 ! get element number 1, see pg 162 of v8.8 user guide
                    aparam (i,j)= B(K)
                    bparam(i,j)= B(K+1)
                endif
            enddo ! j
        enddo ! i
        if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call scathild(dgammares, x, VU, n, aparam*(T + delT/2D0)+bparam,R, T + delT/2D0)
            call scathild(dgammaresd, x, VD, n, aparam*(T - delT/2D0)+bparam,R, T - delT/2D0)
            dgammares = (dgammares - dgammaresd)/delT
        endif ! kcalc < 1
        if (MOD(kcalc,2) .gt. 0) then !kcalc is odd
            Aij = aparam*T+bparam
            if(kop(1).gt.1) then
                write(global_nh,*) ' *** GMUscathild Physical aij'
                write(global_nh,'(6G12.5)') aparam
                write(global_nh,*) ' *** GMUscathild Physical bij'
                write(global_nh,'(6G12.5)') bparam
                write(global_nh,*) ' *** GMUscathild Physical AIJ'
                write(global_nh,'(6G12.5)') AIJ
            endif
            call scathild(gammares, x, V, n, Aij,R, T)
        endif ! kcalc is odd

	case (3) ! Nagata1
        LGAMPAR=DMS_IFCMNC('GMUNAG')  ! Nagata1 parameters
        NEL = 2
        LNDS = LGAMPAR+NEL*NCOMP_NCC*NCOMP_NCC*(NDS-1)
        LNDSCOV = LCOV

        do i=1,n
            do j = 1,n
                if (i.EQ.j) then
                    aparam (i,j)= 0D0
                    bparam(i,j)= 0D0
                else
                    L = LNDS+NEL*NCOMP_NCC*(IDX(J) -1)+NEL*(IDX(I)-1)
                    K = L+1 ! get element number 1, see pg 162 of v8.8 user guide
                    aparam (i,j)= B(K)
                    bparam(i,j)= B(K+1)
               endif
            enddo ! j
        enddo ! i

        if (kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call nagata1(dgammares, n, x, VU, aparam + bparam/(T + delT/2D0))
            call nagata1(dgammaresd, n, x, VD, aparam + bparam/(T - delT/2D0))
            dgammares = (dgammares - dgammaresd)/delT
			WRITE(GLOBAL_NH,*) 'Calculating temperature derivative of residual'
        endif
        if (MOD(kcalc,2).gt.0) then !kcalc is odd
            if(kop(1).gt.1) then
                write(global_nh,*) ' *** GMUnagata Physical aij'
                write(global_nh,'(6G12.5)') aparam
                write(global_nh,*) ' *** GMUnagata Physical bij'
                write(global_nh,'(6G12.5)') bparam
                write(global_nh,*) ' *** GMUnagata Physical aij + bij/T'
                write(global_nh,'(6G12.5)') aparam + bparam/T
            endif
            ! call nagata1(gammares, n, x, UNIQR, aparam)
            call nagata1(gammares, n, x, V, aparam + bparam/T)
        endif

end select
! ********************************** end of gamma phyiscal ********************

! ********************************** Wertheim contribution ***********************
call loadsites(T,n,x,idx,kop, nsitesp, KAD, eps, compp, sitep, aspmx)

if (aspmx.eq.1) then ! only execute if association exists
	if (mod(kcalc,2) .gt. 0) then ! odd kcalc, calculate properties
		rho = 1D0/V ! molar density mol/cc
		rhomix = 1D0/(dot_product(x,V)) ! mol/cc

		CALL calc_gammaw(gammaw, kcalc, kop(1), kop(2), n, nsitesp, compp, sitep, &
							   T, rhomix, rho, KAD, eps, bvol, &
							   PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
	endif

	if (kcalc .gt. 1) then ! need temperature derivative
		rhoU = 1D0/VU ! molar density mol/cc
		rhomixU = 1D0/(dot_product(x,VU)) ! mol/cc
		rhoD = 1D0/VD ! molar density mol/cc
		rhomixD = 1D0/(dot_product(x,VD)) ! mol/cc

		CALL calc_gammaw(dgammaw,  kcalc, kop(1),kop(2), n, nsitesp, compp, sitep, &
							   (T + delT / 2D0), rhomixU, rhoU, KAD, eps, &
							   bvol, PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
		CALL calc_gammaw(dgammawd, kcalc, kop(1), kop(2), n, nsitesp, compp, sitep, &
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
		call flory(dgammacomb, x, VU, n)
		call flory(dgammacombd, x, VD, n)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	elseif (kop(6) .eq. 1) then
		call flory(dgammacomb, x, VU**(2D0/3D0), n)
		call flory(dgammacombd, x, VD**(2D0/3D0), n)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	endif
    endif !kcalc > 1
    if(mod(kcalc,2).gt.0) then
	if (kop(6) .eq. 0) then
		call flory(gammacomb, x, V, n)
        if (kop(1).gt. 1) WRITE(global_nh,*) '*** Combinatorial Method: Flory using molar volume'
	elseif (kop(6) .eq. 1) then
		call flory(gammacomb, x, V**(2D0/3D0), n)
        if (kop(1).gt. 1) WRITE(global_nh,*) '*** Combinatorial Method: Flory using molar volume^(2/3)'
	endif
	endif
endif

! **********Combinatorial Correction
if ((kop(4) .eq. 1) .or. (kop(6) .eq. 2)) then
    kop(5) = 1; ! force combinatorial correction off if using Wilson, or if combinatorial term is off
    ! note the gammacombcorr is set to zero at the beginning of the routine
    if(kop(1).gt.1) WRITE(global_nh,*) '*** Combinatorial Correction Method: None'
endif
if (kop(5) .ne. 1) then
    if(kcalc.gt.1) then
        select case(kop(5))
        case (0)
            call correqn1262(dgammacombcorr, x, VU, n, bvdw)
            call correqn1262(dgammacombcorrd, x, VD, n, bvdw)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        case (2)
            dgammacombcorr = 0
        case (3)
	        call correqn1262(dgammacombcorr, x, VU, n, bvol)
            call correqn1262(dgammacombcorrd, x, VD, n, bvol)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        end select
    endif
    if(mod(kcalc,2).gt.0) then
        select case(kop(5))
        case (0)
            call correqn1262(gammacombcorr, x, V, n, bvdw)
            if (kop(1).gt. 1) WRITE(global_nh,*) '*** Combinatorial Correction Method: vdW using bvdw'
        case (2)
            call stavgugg(gammacombcorr, x, UNIQR,UNIQQ, n)
            if (kop(1) .gt. 1) WRITE(global_nh,*) '*** Combinatorial Correction Method: Staverman-Guggenheim usig r and q'
        case (3)
	        call correqn1262(gammacombcorr, x, V, n, bvol)
            if (kop(1) .gt. 1) WRITE(global_nh,*) '*** Combinatorial Correction Method: vdW using COVOL'
        end select
    endif
end if

! add logs of physical and wertheim contributions
gamma = gammares + gammaw + gammacomb + gammacombcorr
dgamma = dgammares + dgammaw + dgammacomb + dgammacombcorr

if(kop(1) .gt. 0) then ! write to history and/or .csv files
	  WRITE(GLOBAL_NH,*) '*** GAMMA RESULTS'
	  WRITE(GLOBAL_NH,900) KCALC, N, T, P
 900  FORMAT(' KCALC, N, T, P =', I3, 1X, I3, F10.3, 2X, E13.6)
	  if(mod(kcalc,2).gt.0) then !kcalc is odd, write gammas to history
		WRITE(GLOBAL_NH,1010)
		DO 400 I = 1, N ! write gamma contributions to history
		 j = LIDSCC + 2*(idx(i)-1)
		 WRITE(GLOBAL_NH,1000) IDX(I), (IB(j+k), k=1,2), X(I), gammares(i),gammacomb(i),gammacombcorr(i), gammaw(i), GAMMA(I), dgamma(i)
 400   CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
	   if(kop(1) .gt. 2) then
		! write to vol.csv
		open(unit=12, file="vol.csv", status='unknown', action='write', position='append')
		WRITE(12,1002) T, ',', 1D0/rhomix, (', ', IDX(I), ',  ', (IB(LIDSCC + 2*(idx(i)-1)+k), k=1,2), ",", X(I), ",", V(i), i=1,N)
		close(12)
		! write to gammas.csv
		open(unit=10, file="gammas.csv", status='unknown', action='write', position='append')
		WRITE(10,1001) T, ',', P, (', ', IDX(I), ',  ', (IB(LIDSCC + 2*(idx(i)-1)+k), k=1,2), ",", X(I), ",", gammares(i),",",gammacomb(i),",",gammacombcorr(i), ",",gammaw(i), ",",GAMMA(I), ",",dgamma(i), i=1,N)
		close(10)
	   endif !kop(1).gt.2
	  endif !mod(kcalc,2)
	  if(kcalc.gt.1) then !write derivative
		WRITE(GLOBAL_NH,1011)
		DO 410 I = 1, N
		 j = LIDSCC + 2*(idx(i)-1)
		 WRITE(GLOBAL_NH,1000) IDX(I), (IB(j+k), k=1,2), X(I), dgammares(i),dgammacomb(i),dgammacombcorr(i), dgammaw(i), GAMMA(I), dgamma(i)
 410    CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
		if(kop(1).gt.2) then
		! write to HXS.csv
		open(unit=66, file="HXS.csv", status='unknown', action='write', position='append')
		WRITE(66,1001) T, ",", P, (', ', IDX(I), ",  ", (IB(LIDSCC + 2*(idx(i)-1)+k), k=1,2), ",", X(I), ",", -dgammares(i)*R*(T**2), ",", -dgammacomb(i)*R*(T**2), ",", -dgammacombcorr(i)*R*(T**2), ",", -dgammaw(i)*R*(T**2), ",", -dgamma(i)*R*(T**2), ',', dgamma(i), i=1,N)
		close(66)
		endif !if kcalc.gt.1
		endif !kop(1).gt.2

		WRITE(GLOBAL_NH,*) 'End of function call'
	  WRITE(GLOBAL_NH,"(120('*'))")
endif !kop(1).gt.0

!
!
!
!     Format statements
!
 1000   FORMAT (I3,1X, 2A4, 7F15.6)
 1001	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, 2A4, 7(A1,F15.6))) ! format for gamma, hex csv files, up to 100 components
 1002	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, 2A4, 2(A1,F15.6))) ! format for vol csv files, up to 100 components
 1010   FORMAT (' I     NAME          X         lnGAMMAres    lnGAMMAcomb  lnGAMMAcombcorr       lnGAMMAw       lnGAMMA      dlnGAMMA')
 1011   FORMAT (' I     NAME          X         dlnGAMMAres   dlnGAMMAcomb dlnGAMMAcombcorr     dlnGAMMAw       lnGAMMA      dlnGAMMA')
!
return
end subroutine GMU
