!*****************************************************
!
!  The file is divided into the following sections.
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!  To quickly find a section search for '* Section X' where X is the section.
!
!*********************************************************
! Option codes kop() set by the association parameter input file
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
!
! In addition, see Kcalc below to determine if temperature derivative is calculated.
! 	1 - calculate only gammas
! 	2 - calculate only gamma derivative, d(ln gamma)/dT
!	3 - calculate both gamma and gamma derivative
! 	The temperature derivative is determined by finite differences in GAMMAPA.F90.
!

! **************************************************
!
!  ** Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
! *******************************************************

MODULE GPA_CONSTANTS
	REAL*8, PARAMETER :: avoNum = 602.214076d0, kB = 0.01380649D0, R = avoNum*kB
	INTEGER, PARAMETER :: outfile = 2000 ! unit number for regular output
	INTEGER, PARAMETER :: debugfile = 2001 ! unit for debug file output
	! Units 10, 11, 12 are used for gamma.csv, hxs.csv, vol.csv respectively, within code
	CHARACTER(4), PARAMETER :: missing = '1D35'
END MODULE GPA_CONSTANTS

MODULE GPA_SITENSPECIES
	! See loadsites for descriptions
    TYPE siteinfo
        INTEGER	:: host = 0	    ! component index for site host. Currently each site on only one host
        CHARACTER(20) :: name = ''    !site name
        INTEGER	:: id = 0	! site id value (e.g. 2 or 502)
        INTEGER	:: noccur = 0	! for use when each site on only one host
        INTEGER	:: packed = 0	! packed id
        REAL*8	:: xhost = 0D0	! mole fraction of host
    END TYPE siteinfo

    TYPE species
        CHARACTER(50):: name
        INTEGER :: id
        LOGICAL	:: selfassoc = .FALSE.	! False by default, True if self associating
        LOGICAL	:: dnr = .FALSE.	! False by default, True if electron donors on species
        LOGICAL	:: acpt = .FALSE.	! False by default, True if electron acceptors on species
        INTEGER	:: nsite = 0		! number of sites on species
        INTEGER	:: sid(5) = (/ 0,0,0,0,0 /)	! site (defined above) indexes to point to the sites on the species, vector size limits the # sites on a species
        REAL*8	:: x = 0D0 		! mole fraction
    END TYPE species

    type (siteinfo), dimension(:), allocatable :: site ! dimensioned later to be the number of sites present
    type (species), dimension(:), allocatable :: comp ! components present
    real*8, dimension(:), allocatable :: x
    real*8, dimension(:,:), allocatable :: KAD, eps
    real*8, dimension(:), allocatable :: PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok ! PCSAFT
    integer :: nc, nsites, aspmx
    ! aspmx =1 if association or solvation is present in mixture (or pure components)

END MODULE GPA_SITENSPECIES

MODULE GPA_VOLUMES
	integer, allocatable :: veq(:)
	real*8, dimension(:), allocatable :: bvol, r_uniquac, q_uniquac, vstd
	real*8, dimension(:,:), allocatable :: vparms
END MODULE GPA_VOLUMES

MODULE GPA_PHYS_PARMS
	real*8, dimension(:,:), allocatable :: aparam, bparam, alpha, Aij
END MODULE GPA_PHYS_PARMS

! **************************************************
!
!  Section 1 - Modules
!  ** Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
!****************************************************
! Subroutine FINDINDEX
!
! Find the component index for a given component id in comp structure
!
!***************************************************
SUBROUTINE FINDINDEX(index,compid)

use GPA_SITENSPECIES, only: nc, comp
integer, intent(in) :: compid
integer, intent(out) :: index

do i = 1, nc
 if(comp(i)%id .EQ. compid) EXIT
END DO
index = i
END SUBROUTINE FINDINDEX
!********************************************

! Section 2 continued

!*********************************************
!
! Subroutine PARSE_LINE
!
! Used for reading tab separated input files one line at at time.
!
!*******************************************************
subroutine parse_line(start_index,data_line,data_field)
! This routine parses lines that are tab delimited.
! The data_line is the complete line.
! The start_index is where in data_line to start searching.
! The data_field is string before the next tab.
implicit none
character*500, intent(IN) :: data_line
character*100, intent(OUT) :: data_field
integer, intent(INOUT) :: start_index
integer :: delim_index
CHARACTER(4) :: missing = '1D35'

delim_index = SCAN(data_line(start_index:), achar(9))
IF (delim_index .EQ. 0) THEN
data_field = data_line(start_index:)
ELSE IF (delim_index .EQ. 1) THEN ! EMPTY FIELD
data_field = missing ! INSERT value for missing
ELSE
data_field = data_line(start_index : start_index + delim_index - 2) ! Adjust for 1-based indexing of SCAN
END IF
start_index = start_index + delim_index

end subroutine parse_line
!*********************************************************

! Section 2 continued

!*********************************************************
! Subroutine VLU
!
! Calculate molar volume given equation number and coefficients
!
!***********************************************************

subroutine VLU(V,IDI,VEQN,DNLDIP,T)

use GPA_CONSTANTS, only: debugfile

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
		WRITE(debugfile,99) IDI
 99		FORMAT('Density is not programmed for component',I4)
		WRITE(debugfile,100) VEQN
 100	FORMAT ('Equation ',I4, ' is not programmed. Contact the .dll authors,')
End Select
end subroutine vlu

! **************************************************
!
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  ** Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
! *******************************************************

! gammapa_main.f90
!
! main entry point for GAMMA_PA application
! gammapa_main loads parameters and shows example calls to gamma_pa
! where the calculations are executed.
!
!*********************************************************
!

PROGRAM GAMMAPA_MAIN

use GPA_CONSTANTS
! constants includes the outfile id which should not be used for another file
use GPA_SITENSPECIES
use GPA_VOLUMES
USE GPA_PHYS_PARMS

IMPLICIT NONE

! Please update the version when making more than minor fixes
! Append RTPT or TPT1 so that debug code shows the option used for compiled code.
character(len=10) :: ver = '0.2TPT1'

! todo - track sub and Henry's law components
! INTEGER NSUB,  NSUP

INTEGER KOP(10), KCALC, ioErr
real*8 T, P
REAL*8, dimension(:), allocatable :: gamma, dgamma

integer i,j,k

! for physical and combinatorial models
! aparam, bparam for NRTL, Nagata, Wilson, tau for NRTL, Aij for SH

print *, 'Welcome to the GAMMAPA program.'
print *, ' '
print *, 'Output will be written to <project>/Output/GAMMAPAout.txt'
print *, ' '
! Prepare for output
OPEN(outfile, ioStat=ioErr, file="Output\GAMMAPAout.txt")
if (ioErr.ne.0) then
	print *, 'Could not open OUTPUT/GAMMAPAout.txt. Error Code = ', ioErr
else
	print *, 'File Output/GAMMAPAout.txt opened successfully.'
end if
OPEN(debugfile, ioStat=ioErr, file="Output\GAMMAPAdebug.txt")
if (ioErr.ne.0) then
	print *, 'Could not open OUTPUT/GAMMAPAdebug.txt. Error Code = ', ioErr
else
	print *, 'File Output/GAMMAPAdebug.txt opened successfully.'
end if
print *, ' '
print *, 'You will be prompted for the file names of the two input files.'
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
        call loadwilson()
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMUwilson Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMUwilson Physical bij'
            write(debugfile,'(6G12.5)') bparam
        endif

    case (2) !scatchard hildebrand (eqn Elliott/Lira 12.24,12.25)
        call loadsh()
        if(kop(1).gt.1) then
            write(debugfile,*) ' *** GMUscathild Physical aij'
            write(debugfile,'(6G12.5)') aparam
            write(debugfile,*) ' *** GMUscathild Physical bij'
            write(debugfile,'(6G12.5)') bparam
        endif

	case (3) ! Nagata1
        call loadnagata()
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

print *, 'Enter the Temperature(K) and Pressure (bar)'
print *, 'Use D0 on the end to read as double precision, e.g. 298.15D0'
read(*,*) T, P

print *, 'Enter the Kcalc variable to indicate the calculations desired.'
print *, '1 - calculate only gammas'
print *, '2 - calculate only gamma derivative, d(ln gamma)/dT'
print *, '3 - calculate both gamma and gamma derivative'
print *, 'Note that Kcalc > 1 does a lot of calcs with T derivatives,'
print *, 'so some debugging output is supressed relative to Kcalc = 1.'
read(*,*) Kcalc
print *, 'You entered:', Kcalc

write(outfile,'(A, I3 )') 'Kcalc ', Kcalc
write(outfile, '(A, 2F10.3)') 'T(K) P(bar) ', T, P

! set composition of interest
! recall x is shared in gpa_sitenspecies

!******for meoh-cyclhex-assoc.txt and meoh-cyclhex-nrtl.txt ******
! example loop for a binary.
! set T = 298.15D0 and P = 1D0
!write(outfile,'(A)') 'x, lngamma, gamma, hex'
!do i = 1, 11
!    x(1) = dble(i-1)/10D0
!    x(2) = 1D0-x(1)
!    call gamma_pa(kop, Kcalc, T, P, gamma, dgamma)
!    write(outfile, '(8F15.6)') X(:), gamma(:), dexp(gamma(:)), -R*T**2*(dot_product(x,dgamma))
!enddo
! ****** end meoh-cyclhex ******************

!****** for casestudy1-assoc.txt and casestudy1-nrtl.txt *****
! set T = 298.15K for this composition
! x = (/ 0.33D0, 0.33D0, 0.34D0 /)
! set T = 347.125 (nrtl) T = 343.358 (nrtla) for this composition
!x = (/ 0.165D0, 0.165D0, 0.67D0 /)
! call gamma_pa(kop, Kcalc, T, P, gamma, dgamma)
! write(outfile, '(10F15.6)') X(:), gamma(:), dexp(gamma(:)), -R*T**2*(dot_product(x,dgamma))
! set T = 298.15K for these loops
!write(outfile,'(A)') 'x, lngamma, gamma, hex'
!do i = 1, 11
! use one line or the other
!    x = (/ 0D0, dble(i-1)/10D0, 1D0-dble(i-1)/10D0 /)
!   x = (/ dble(i-1)/10D0, 0D0, 1D0-dble(i-1)/10D0 /)
!   call gamma_pa(kop, Kcalc, T, P, gamma, dgamma)
!   write(outfile, '(8F15.6)') X(:), gamma(:), dexp(gamma(:)), -R*T**2*(dot_product(x,dgamma))
!enddo
!************ end casestudy1 ******************

!***********CaseStudy2i,v,ix***************
!*****CaseStudy2iCPA-assoc.txt, CaseStudy2iCPA-nrtl.txt*****
!*****CaseStudy2iPCS-assoc.txt, CaseStudy2iPCS-nrtl.txt*****
 x = (/ 0.05D0, 0.0D0, 0.25D0, 0.2D0, 0.05D0, 0.019D0, 0.431D0 /) ! Case 2v
 x = (/ 0.25D0, 0.1D0, 0.05D0, 0.1D0, 0.15D0, 0.2D0, 0.15D0 /) ! Case 2ix
 x = (/ 0.05D0, 0.075D0, 0.3D0, 0.2D0, 0.1D0, 0.01D0, 0.265D0 /) ! Case 2i
 call gamma_pa(kop, Kcalc, T, P, gamma, dgamma)
 write(outfile,'(A)') 'Gammas are not calculated when Kcalc = 2'
 write(outfile, '(A5,7F15.6,/,A5,7F15.6,/,A5,7F15.6,/)') 'x', X(:), 'ln(g)', gamma(:), 'g', dexp(gamma(:))
 write(outfile,'(A)') 'Hxs is not calculated if Kcalc = 1'
 write(outfile,'(A10, F15.6)') 'Hxs(J/mol)', -R*T**2*(dot_product(x,dgamma))
!*********** end CaseStudy2i,v,ix***************

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

END PROGRAM GAMMAPA_MAIN
!****************************************************************

! Section 3 continued

!****************************************************************************
!
!  Subrouting LOADSITES
!        read component and site parameters
!		 -'site' for each site holds the index to the site name and id, the species host
!          index and host mole fraction, the number of occurences on the host, the 'packed id'
!		   holds the ids for the species present, which is the same for this application.
!		 - 'comp' for each species holds the mole fraction, the number of sites on the species,
!			the site ids for the sites, a flag to indicate if the species self associates.
!		 - 'aspmx' - 0 if a mixture has no self- or cross-association, 1 if association is present.
!
!****************************************************************************

SUBROUTINE LOADSITES(KOP)

use GPA_CONSTANTS, only: outfile
use GPA_SITENSPECIES
use GPA_VOLUMES

implicit none

! passed variables
integer :: start_index, KOP(10)
! local variables
integer ::  i, j, k, ioErr, hostid, nsiteparm, id1, id2, dnrflagmx, acptflagmx
character*500 :: data_line
character*100 :: data_field, filesites
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
    print *, 'Enter a the name of the tab-delimited sites input file:'
	read(*,'(A)') filesites
	print *, 'You entered:', trim(filesites)
	filesites = 'Input\GAMMAPA\'//filesites
	INQUIRE(FILE=trim(filesites), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'Pure/Assoc file: '//trim(filesites)
OPEN(UNIT=1001, ioStat = ioErr, FILE=trim(filesites))
if (ioErr.ne.0) then
	print *, 'Could not open sites file. Error Code = ', ioErr
else
	print *, 'Sites file opened successfully.'
end if
READ(UNIT=1001, END=106, FMT='(A)') data_line ! title
READ(UNIT=1001, END=106, FMT='(A)') data_line ! KOP
start_index = 1
do i=1,10
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) kop(i)
end do
READ(UNIT=1001, END=106, FMT='(A)') data_line ! title with nc
start_index = 1
call parse_line(start_index,data_line,data_field)
read(data_field,*) nc
if(.not.allocated(comp)) allocate(comp(nc)) ! allocate storage for components
if(.not.allocated(bvol)) allocate(bvol(nc))
if(.not.allocated(r_uniquac)) allocate(r_uniquac(nc))
if(.not.allocated(q_uniquac)) allocate(q_uniquac(nc))
if(.not.allocated(vstd)) allocate(vstd(nc))
if(.not.allocated(PCSAFT_m)) allocate(PCSAFT_m(nc))
if(.not.allocated(PCSAFT_sigma)) allocate(PCSAFT_sigma(nc))
if(.not.allocated(PCSAFT_epsok)) allocate(PCSAFT_epsok(nc))
if(.not.allocated(veq)) allocate(veq(nc))
if(.not.allocated(vparms)) allocate(vparms(7,nc))
READ(UNIT=1001, END=106, FMT='(A)') data_line ! title

! read component information
do i = 1, nc
    READ(UNIT=1001, END=106, FMT='(A)') data_line
    start_index = 1
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) comp(i)%id
    call parse_line(start_index,data_line,data_field)
    comp(i)%name = data_field
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vstd(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) bvol(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) r_uniquac(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) q_uniquac(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) PCSAFT_m(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) PCSAFT_sigma(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) PCSAFT_epsok(i)
    call parse_line(start_index,data_line,data_field) ! source not used
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) veq(i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(1,i) ! store in columns for ease in using debugger
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(2,i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(3,i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(4,i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(5,i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(6,i)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) vparms(7,i)
end do ! nc

! Read site information
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header for sites
start_index = 1
call parse_line(start_index,data_line,data_field)
read(data_field,*) nsites
if(.not.allocated(site)) allocate(site(nsites))
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header line

do i = 1, nsites
    READ(UNIT=1001, END=106, FMT='(A)') data_line
    start_index = 1
    ! The %packed index provides a method to read in more sites that are present in
    ! a function call. In that case the index should be set for the function call.
    ! For this implementation, all sites are always present so the packed id = site index
    site(i)%packed = i
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) site(i)%id
    call parse_line(start_index,data_line,data_field)
    site(i)%name = data_field
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) hostid
        do j = 1, nc ! determine host index
        if(hostid .EQ. comp(j)%id) exit
        end do
    site(i)%host = j
    if (site(i)%id .LT. 500) comp(j).dnr = .TRUE.
    if ((site(i)%id .GT. 499) .AND. (site(i)%id .LT. 1000)) comp(j).acpt = .TRUE.
    if (site(i)%id .GT. 999) THEN
        comp(j).dnr = .TRUE.
        comp(j).acpt = .TRUE.
    end if
    do k = 1, 5 ! add site to component list, find if other sites exist
        if(comp(j)%sid(k) .EQ. 0) exit
    end do
    comp(j)%sid(k) = i
    comp(j)%nsite = k
    call parse_line(start_index,data_line,data_field) ! host name not used
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) site(i)%noccur
end do ! site assigments

! Read site parameters

READ(UNIT=1001, END=106, FMT='(A)') data_line ! header for site params
start_index = 1
call parse_line(start_index,data_line,data_field)
read(data_field,*) nsiteparm

if(.not.allocated(KAD)) allocate(KAD(nsites,nsites))
if(.not.allocated(eps)) allocate(eps(nsites,nsites))
KAD = reshape( (/ (0D0, k=1, nsites*nsites) /), [nsites,nsites])
eps = reshape( (/ (0D0, k=1, nsites*nsites) /), [nsites,nsites])

READ(UNIT=1001, END=106, FMT='(A)') data_line ! header for site parms

DO i = 1, nsiteparm
    READ(UNIT=1001, END=106, FMT='(A)') data_line
    start_index = 1
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) id1
    call parse_line(start_index,data_line,data_field) ! name not used
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) id2
        DO j = 1, nsites
            if(site(j)%id .EQ. id1) EXIT ! found site1 index
        END DO
        DO k = 1, nsites
            if(site(k)%id .EQ. id2) EXIT ! found site2 index
        END DO
    call parse_line(start_index,data_line,data_field) ! name not used
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) KAD(j,k)
    KAD(k,j) = KAD(j,k)
    call parse_line(start_index,data_line,data_field)
    read(data_field,*) eps(j,k)
eps(k,j) = eps(j,k)
END DO ! nsiteparm

! Determine if components self-associate
! Flags for mixture, nonzero if that type of site is found in mixuture
dnrflagmx = 0
acptflagmx = 0
aspmx = 0 ! =1 if Association or Solvation is Present in mixture

do i = 1, nc
    ! program stucture assumes that each site appears in only on species
    if (comp(i)%dnr) dnrflagmx = 1
    if (comp(i)%acpt) acptflagmx = 1
    if( (comp(i)%dnr .AND. comp(i)%acpt) ) THEN
        comp(i)%selfassoc = 1
    END IF
enddo

IF ((dnrflagmx + acptflagmx) .GT. 1) aspmx = 1

106 CLOSE(1001)

END SUBROUTINE LOADSITES
!*************************************************

! Section 3 continued

!**********************************************
!
!  Subroutine loadnrtl - loads NRTL parameters
!
!**********************************************

SUBROUTINE LOADNRTL()

use GPA_CONSTANTS, only: outfile
use GPA_SITENSPECIES, only: comp, nc
USE GPA_PHYS_PARMS, ONLY: aparam, bparam, alpha

integer :: start_index, npair, index1, index2, compid1, compid2

! local variables
integer :: i, j, id1, id2, ioErr
character*500 :: data_line
character*100 :: data_field
character*100 :: filenrtl
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
	print *, 'Enter a the name of the tab-delimited NRTL parameter file:'
	read(*,'(A)') filenrtl
	print *, 'You entered:', trim(filenrtl)
	filenrtl = 'Input\GAMMAPA\'//filenrtl
	INQUIRE(FILE=trim(filenrtl), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'NRTL file: '//trim(filenrtl)
OPEN(UNIT=1001, ioStat = ioErr, FILE=trim(filenrtl))
if (ioErr.ne.0) then
	print *, 'Could not open NRTL file. Error Code = ', ioErr
else
	print *, 'NRTL file opened successfully.'
end if
start_index = 1
READ(UNIT=1001, END=106, FMT='(A)') data_line
call parse_line(start_index,data_line,data_field)
read(data_field,*) npair
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header

DO i = 1, npair
	start_index = 1
	READ(UNIT=1001, END=106, FMT='(A)') data_line
	call parse_line(start_index,data_line,data_field) !pairid
	call parse_line(start_index,data_line,data_field) !id1
	read(data_field,*) compid1
	call parse_line(start_index,data_line,data_field) !name
	call parse_line(start_index,data_line,data_field) !id2
	read(data_field,*) compid2
	call parse_line(start_index,data_line,data_field) !name
	! find index for components
	call findindex(index1,compid1)
	call findindex(index2,compid2)
	call parse_line(start_index,data_line,data_field) !a12
	read(data_field,*) aparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !b12
	read(data_field,*) bparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !a21
	read(data_field,*) aparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !b21
	read(data_field,*) bparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !alpha
	read(data_field,*) alpha(index1,index2)
	alpha(index2,index1) = alpha(index1,index2)
END DO

106 CLOSE(1001)

END SUBROUTINE LOADNRTL
!****************************************************

! Section 3 continued

!**********************************************
!
!  Subroutine LOADWILSON - loads WILSON parameters
!
!**********************************************

SUBROUTINE LOADWILSON()

use GPA_CONSTANTS, only: outfile
use GPA_SITENSPECIES, only: comp, nc
USE GPA_PHYS_PARMS, ONLY: aparam, bparam

integer :: start_index, npair, index1, index2, compid1, compid2

! local variables
integer :: i, j, id1, id2, ioErr
character*500 :: data_line
character*100 :: data_field
character*100 :: filewilson
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
	print *, 'Enter a the name of the tab-delimited Wilson parameter file:'
	read(*,'(A)') filewilson
	print *, 'You entered:', trim(filewilson)
	filewilson = 'Input\GAMMAPA\'//filewilson
	INQUIRE(FILE=trim(filewilson), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'Wilson file: '//trim(filewilson)
OPEN(UNIT=1001, ioStat = ioErr, FILE=trim(filewilson))
if (ioErr.ne.0) then
	print *, 'Could not open WILSON file. Error Code = ', ioErr
else
	print *, 'WILSON file opened successfully.'
end if

start_index = 1
READ(UNIT=1001, END=106, FMT='(A)') data_line
call parse_line(start_index,data_line,data_field)
read(data_field,*) npair
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header

DO i = 1, npair
	start_index = 1
	READ(UNIT=1001, END=106, FMT='(A)') data_line
	call parse_line(start_index,data_line,data_field) !pairid
	call parse_line(start_index,data_line,data_field) !id1
	read(data_field,*) compid1
	call parse_line(start_index,data_line,data_field) !name
	call parse_line(start_index,data_line,data_field) !id2
	read(data_field,*) compid2
	call parse_line(start_index,data_line,data_field) !name
	! find index for components
	call findindex(index1,compid1)
	call findindex(index2,compid2)
	call parse_line(start_index,data_line,data_field) !a12
	read(data_field,*) aparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !b12
	read(data_field,*) bparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !a21
	read(data_field,*) aparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !b21
	read(data_field,*) bparam(index2,index1)
END DO

106 CLOSE(1001)

END SUBROUTINE LOADWILSON
!*********************************************

! Section 3 continued

!**********************************************
!
!  Subroutine LOADSH - loads Scatchard-Hildebrand parameters
!
!**********************************************

SUBROUTINE LOADSH()

use GPA_CONSTANTS, only: outfile
use GPA_SITENSPECIES, only: comp, nc
USE GPA_PHYS_PARMS, ONLY: aparam, bparam

integer :: start_index, npair, index1, index2, compid1, compid2

! local variables
integer :: i, j, id1, id2, ioErr
character*500 :: data_line
character*100 :: data_field
character*100 :: filesh
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
	print *, 'Enter the name of the tab-delimited Scatchard-Hildebrand parameter file:'
	read(*,'(A)') filesh
	print *, 'You entered:', trim(filesh)
	filesh = 'Input\GAMMAPA\'//filesh
	INQUIRE(FILE=trim(filesh), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'Scatchard-Hildebrand file: '//trim(filesh)
OPEN(UNIT=1001, ioStat = ioErr, FILE=trim(filesh))
if (ioErr.ne.0) then
	print *, 'Could not open Scatchard-Hildebrand file. Error Code = ', ioErr
else
	print *, 'Scatchard-Hildebrand file opened successfully.'
end if
start_index = 1
READ(UNIT=1001, END=106, FMT='(A)') data_line
call parse_line(start_index,data_line,data_field)
read(data_field,*) npair
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header

DO i = 1, npair
	start_index = 1
	READ(UNIT=1001, END=106, FMT='(A)') data_line
	call parse_line(start_index,data_line,data_field) !pairid
	call parse_line(start_index,data_line,data_field) !id1
	read(data_field,*) compid1
	call parse_line(start_index,data_line,data_field) !name
	call parse_line(start_index,data_line,data_field) !id2
	read(data_field,*) compid2
	call parse_line(start_index,data_line,data_field) !name
	! find index for components
	call findindex(index1,compid1)
	call findindex(index2,compid2)
	call parse_line(start_index,data_line,data_field) !a12
	read(data_field,*) aparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !b12
	read(data_field,*) bparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !a21
	read(data_field,*) aparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !b21
	read(data_field,*) bparam(index2,index1)
END DO

106 CLOSE(1001)

END SUBROUTINE LOADSH
!************************************************

! Section 3 continued

!**********************************************
!
!  Subroutine LOADNAGATA - loads Nagata parameters
!
!**********************************************

SUBROUTINE LOADNAGATA()

use GPA_CONSTANTS, only: outfile
use GPA_SITENSPECIES, only: comp, nc
USE GPA_PHYS_PARMS, ONLY: aparam, bparam

integer :: start_index, npair, index1, index2, compid1, compid2

! local variables
integer :: i, j, id1, id2, ioErr
character*500 :: data_line
character*100 :: data_field
character*100 :: filenagata
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
	print *, 'Enter a the name of the tab-delimited Nagata parameter file:'
	read(*,'(A)') filenagata
	print *, 'You entered:', trim(filenagata)
	filenagata = 'Input\GAMMAPA\'//filenagata
	INQUIRE(FILE=trim(filenagata), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'Nagata file: '//trim(filenagata)
OPEN(UNIT=1001, ioStat = ioErr, FILE=trim(filenagata))
if (ioErr.ne.0) then
	print *, 'Could not open Nagata file. Error Code = ', ioErr
else
	print *, 'Nagata file opened successfully.'
end if
start_index = 1
READ(UNIT=1001, END=106, FMT='(A)') data_line
call parse_line(start_index,data_line,data_field)
read(data_field,*) npair
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header

DO i = 1, npair
	start_index = 1
	READ(UNIT=1001, END=106, FMT='(A)') data_line
	call parse_line(start_index,data_line,data_field) !pairid
	call parse_line(start_index,data_line,data_field) !id1
	read(data_field,*) compid1
	call parse_line(start_index,data_line,data_field) !name
	call parse_line(start_index,data_line,data_field) !id2
	read(data_field,*) compid2
	call parse_line(start_index,data_line,data_field) !name
	! find index for components
	call findindex(index1,compid1)
	call findindex(index2,compid2)
	call parse_line(start_index,data_line,data_field) !a12
	read(data_field,*) aparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !b12
	read(data_field,*) bparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !a21
	read(data_field,*) aparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !b21
	read(data_field,*) bparam(index2,index1)
END DO

106 CLOSE(1001)

END SUBROUTINE LOADNAGATA

! **************************************************
!
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  ** Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
! *****************************************************

!*********************************************************
!
! SUBROUTINE GAMMA_PA
!
! Implement calls to calculate Wertheim, Physical, Combinatorial
! calculations and combine the results for Gammas and temperature
! derivative of gamma
!
!**********************************************************
SUBROUTINE GAMMA_PA(kop, Kcalc, T, P, gamma, dgamma)

use GPA_CONSTANTS
use GPA_SITENSPECIES
use GPA_VOLUMES
use GPA_PHYS_PARMS

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
        write(debugfile,'(2A20, 2I4, 3G12.5)') site(l)%name,site(k)%name,l,k,kad(l,k)*(dexp(eps(l,k)/T)-1D0), eps(l,k), kad(l,k)
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
		open(unit=12, file="Output\vol.csv", status='unknown', action='write', position='append')
		WRITE(12,1002) T, ',', 1D0/rhomix, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", V(i), i=1,nc)
		close(12)
		! write to gammas.csv
		open(unit=10, file="Output\gammas.csv", status='unknown', action='write', position='append')
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
		open(unit=11, file="Output\HXS.csv", status='unknown', action='write', position='append')
		WRITE(11,1001) T, ",", P, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", -dgammares(i)*R*(T**2), ",", -dgammacomb(i)*R*(T**2), ",", -dgammacombcorr(i)*R*(T**2), ",", -dgammaw(i)*R*(T**2), ",", -dgamma(i)*R*(T**2), ',', dgamma(i), i=1,nc)
		close(11)
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
end subroutine GAMMA_PA

! **************************************************
!
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  ** Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
! *****************************************************


!***************************************************
!
! Subroutine CALC_GAMMAW
!
! Calculate the wertheim contribution to activity coefficients
!   Calls for X calculation
!   Calls Helmholtz derivative calculations
!
!***************************************************************

SUBROUTINE calc_gammaw(gammaw, kcalc, kop1, gcalcFlag, n, ns, comp, site, T, &
			rho_mix, rho_pure, KAD, epsADok, bvol, sigma, m, epsok)
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
!	- bvol		(n)		covol [cm3/mol]
!	- sigma		(n)		(gcalcFlag=3)
!	- m		    (n)		(gcalcFlag=3)
!	- epsok		(n)		(gcalcFlag=3)
!	- site		(ns)		site info array for sites present
!	- KAD		(ns,ns)		[cm3/mol] for CPA or ESD, dimensionless for PCSAFT
!	- epsADok	(ns,ns)		[K]
!
! Parameter:
!	- KADplus	(ns,ns)		Same units as KAD, KADplus = KAD(exp(epsok/T)-1)

USE GPA_SITENSPECIES, only: siteinfo, species !pass other variables for legacy reasons
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
!**********************************************************

! Section 5 continued

!*********************************************************
!
!   Subroutine CALC_DADNK
!
!   Calculates the derivative of Helmholtz wrt n_k
!
!********************************************************

SUBROUTINE calc_dAdnk(dAdnk, kcalc, kop1, gcalcFlag, n, ns, comp, site, T, rho_mix, &
			rho_pure, KADplus, bvol, sigma, m, epsok)

use GPA_CONSTANTS
USE GPA_SITENSPECIES, only: siteinfo, species !pass other variables for legacy reasons

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
        write(debugfile,*)'Pure component-------------'
    else
        write(debugfile,*)'Site hosts (for function call)-----------------'
        write(FMT,'("(",I0,"(I5))")') ns ! make format string for ns columns
        write(debugfile,FMT) site(:)%host
        write(debugfile,*)'Site mole fraction (for function call)-----------------'
        write(FMT,'("(",I0,"(F12.5))")') ns ! make format string for ns columns
        write(debugfile,FMT) site(:)%xhost
    endif
	write(debugfile,*) 'Site name'
	write(FMT,'("(",I0,"(A,2X))")') ns ! make format string for ns columns
    write(debugfile,FMT) site(:)%name
    write(debugfile,*)'Site fractions, in order of sites'
	write(FMT,'("(",I0,"(F12.5))")') ns ! make format string for ns columns
    write(debugfile,FMT) Y
    write(debugfile,*)'del matrix by row'
	write(FMT,'("(",I0,"(G12.5))")') ns ! make format string for ns columns
    write(debugfile,FMT) ((del(i,j), j=1,ns), i=1,ns)
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
!**************************************************

! Section 5 continued

!****************************************************
!
! Subroutine GCALC
!
!   Calculates the g(d) term
!       Calls for PCSAFT calculation if necessary
!       For PCSAFT, each site is mapped to the component rdf
!
!********************************************************

SUBROUTINE gcalc(del, gterm1, gterm2, gcalcFlag, n, ns, KADplus, &
		comp, site, T, rho_mix, bvol, rho_pure, sigma, m, epsok)
!
!
! Note: no input checkings are done, need to be done before calling the subroutine
USE GPA_SITENSPECIES, only: siteinfo, species ! for legacy reasons
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
!****************************************************

! Section 5 continued

!*****************************************************
!
! Subroutine calc_gterms_PCSAFT
!
!   The rdf derivatives for each component for PCSAFT are calculated
!
!****************************************************

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
!*************************************************

! Section 5 continued

!************************************************
!
! Subroutine CALCX
!
!   Determines the fraction unbonded for each site
!
!***************************************************

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

! **************************************************
!
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  ** Section 6 - Physical Models
!  Section 7 - Combinatorial Models and Combinatorial Corrections
!
!*****************************************************
!
!   Subroutine NRTL
!       Requires tau and alpha as matrices
!
!*********************************************************

subroutine nrtl(gamma, n, x, tau, alpha)
implicit none
intent(out):: gamma
intent(in):: n, x,tau,alpha
integer n
real*8 gamma(n), x(n), tau(n,n), alpha(n,n), G(n,n), Y(n,n)
real*8 term1(n), term2(n), inverseterm2(n), squareinverseterm2(n), part1(n), part2(n), part3(n), loggamma(n)

! This program calculates activity coefficients using the NRTL equation.
! Each call to nrtl evaluates one composition.
! x = vector of n mole fractions for an n-component mixture, length n.
! tau = tau values, square matrix, n x n.
! alpha = alpha values, square matrix, n x n.
!   The alpha matrix is symmetrical with zeros on the diagonals.

!This will create a matrix with all the rows equal to the vector of
!compositions
Y = spread(x,1,n)

G=exp(-alpha*tau)

term1=matmul(x,(tau*G))
term2=matmul(x,G)
inverseterm2 = (1D0/term2)
squareinverseterm2=inverseterm2**2

part1=(term1/term2)
part2= matmul(inverseterm2,transpose(Y*tau*G));
part3 = matmul(term1*squareinverseterm2,transpose(Y*G));

gamma = part1 + part2 - part3;

return

end subroutine nrtl

!**************************************
!
!	Subroutine WILSON
!		Wilson physical method
!
!***************************************

! This program calculates activity coefficients using the wilson equation.
! Each call to wilson evaluates one composition.
! x = vector of n mole fractions for an n-component mixture, length n.
! Lambda = Lambda values, square matrix, n x n. Diagonals are UNITY.
! Lambda is calculated before this routine is used. In that way, the
! user can include whatever parameter temperature dependence is desired
! without modifying this routine.

SUBROUTINE wilson(gamma,n,xin,Lambda)
IMPLICIT NONE
INTENT(IN)::n,Lambda,xin
INTENT(OUT)::gamma

! declare variables
INTEGER n
REAL*8 Lambda(n,n),Y(n,n)
REAL*8 term1(n),inverseterm1(n),term2(n),gamma(n),x(n),xin(n)

! make sure x is not identically zero
! ctl - Looking at the formula, I don't think this is necessary
! x = xin+1e-37
x = xin

Y = spread(x,1,n)

term1 = matmul(x,transpose(Lambda))
inverseterm1 = 1D0/term1

term2 = matmul(inverseterm1,(transpose(Y)*Lambda))

!Note: the term 'gamma' is actually ln(gamma).
gamma = 1D0 - log(term1) - term2

END SUBROUTINE wilson
!***************************************************

! Section 6 continued

!*************************************************
!
!	Subroutine SCATHILD
!		Calculate Schatchard-Hildebrand Physical Model
!
!************************************************
subroutine scathild(gamma, x, V, n, Aij, R, T)
implicit none
intent(out):: gamma
intent(in):: x, n, V, Aij , R, T
integer n
real*8 gamma(n), x(n), V(n), Vmix, phi(1,n), R, T, Aij(n,n)
real*8 term1(1,1), term2(1,n), RTlnGamma(1,n)

!This function calculates scatchard-hildebrand residual term

Vmix = dot_product(x,V)
phi = (reshape(x,shape(phi))*reshape(V,shape(phi)))/Vmix
term1 = 0.5D0 * matmul(phi,matmul(Aij,transpose(phi)));
term2 = transpose(matmul(Aij,transpose(phi)));
RTlnGamma = reshape(V,shape(phi))*(term2-term1(1,1));
! gamma is actually ln(gamma) for aspens
gamma = reshape((RTlnGamma/R/T),shape(gamma));
return

end subroutine scathild
!***************************************************

! Section 6 continued

!***************************************************
!
!	Subroutine NAGATA1
!		Unnamed model of Nagata, which is really a
!		modification of Wilson
!
!***************************************************
SUBROUTINE nagata1(gamma, n, x, rpara, a)
IMPLICIT NONE
INTENT(IN)::x, rpara, a
INTENT(OUT)::gamma

! This program calculates activity coefficients gamma based on
! Nagata, 1997, Fluid Phase Equilibria, 227-247, Equation (20)
! n = number of components, integer scalar.
! x = vector of n mole fractions for an n-component mixture, shape n.
! rpara = vector of molecular geometric-size parameter of pure components.
! a = adjustable parameters, square matrix with zeros on diagonal, shape n,n.
! T = temperature in Kelvin.

! declear variables
INTEGER n;
REAL*8 r_mix;
REAL*8 x(n),rpara(n);
REAL*8 term2(n),phi(n),term1(n),lnGamma(n),gamma(n);
REAL*8 a(n,n),tau(n,n);

! mixture molecular geometric-size parameter
r_mix = dot_product(x,rpara);
! segment fraction
phi = x*rpara/r_mix;
tau = dexp(a);

term1 = matmul(transpose(tau),phi);
term2 = matmul(tau,x/term1);
gamma = -dlog(term1) + rpara * (1D0-term2) / r_mix;

END SUBROUTINE nagata1

! **************************************************
!
!  Section 1 - Modules
!  Section 2 - Utility and property routines
!  Section 3 - Main routine and load parameters
!  Section 4 - Gamma_PA routine that calls for Wertheim and Physical Models and combines them
!  Section 5 - Wertheim contribution Code and subroutines
!  Section 6 - Physical Models
!  ** Section 7 - Combinatorial Models and Combinatorial Corrections
!
!*****************************************************
!
!   Subroutine FLORY
!       Flory contribution
!       The 2/3 and 3/4 modification shoud be done before
!       calling the routine.
!
!**********************************************************
subroutine flory(gamma,x,V,n)
implicit none
intent(out):: gamma
intent(in):: x, n, V
integer n
real*8 gamma(n), x(n), V(n), Vmix, phi(n)

!This function calculates Flory's term

Vmix = dot_product(x,V)
gamma = log(V/Vmix)+ 1D0- V/Vmix

return

end subroutine flory
!***************************************************

! Section 7 continued

!****************************************************
!
!   Subroutine STAVGUGG
!       Staverman/Guggenheim correction to Flory
!
!************************************************
subroutine stavgugg(gamma,x,r,q,n)
implicit none
intent(out):: gamma
intent(in):: x, n, r, q
integer n
real*8 gamma(n), x(n), r(n), q(n), phioxi(n), sumxr, sumxq, thetaoxi(n)

!This function calculates Staverman-Guggenheim term

sumxr = dot_product(x,r)
sumxq = dot_product(x,q)

!Note: x cancels in the numerator and denominator of phi/theta so here we
!calculate the terms phi and theta over x to avoid a division by zero
phioxi = r/sumxr
thetaoxi = q/sumxq

gamma = -5D0*q*(log(phioxi/thetaoxi)+ 1D0- phioxi/thetaoxi);

return

end subroutine stavgugg
!***********************************************************

! Section 7 continued

!************************************************************
!
!   Subroutine CORREQN1262
!       A vdw correction to entropy using the covolume parameters
!
!***************************************************************

subroutine correqn1262(gamma, x, V, n, bvol)
implicit none
intent(out):: gamma
intent(in):: x, n, V
integer n
real*8 gamma(n), x(n), V(n), bvol(n), Vmix, bmix, denom(n), num(n)

!This function calculates the correction to Flory's combinatorial term as calculated in Bala lab book pg. 119

Vmix = dot_product(x,V) !cc/mol
bmix = dot_product(bvol,x)

denom = 1D0-bvol/V
num = 1D0-bmix/Vmix

gamma = -log(num/denom)+(1D0/num)*(bvol/Vmix - bmix*V/Vmix**2)

return
end subroutine correqn1262
