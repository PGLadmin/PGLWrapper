!  loadsites.f90
!
!****************************************************************************
!
!  SUBROUTINES: loadsites - read component and site parameters
!                  fills all
!              parse_line - support utility for parsing each line via tab chars
!
!****************************************************************************

SUBROUTINE LOADSITES(KOP)

use sitenspecies
use VOLUMES

implicit none

! passed variables
integer :: start_index, KOP(10)
! local variables
integer ::  i, j, k, hostid, nsiteparm, id1, id2, dnrflagmx, acptflagmx
character*500 :: data_line
character*100 :: data_field

OPEN(UNIT=1001, FILE="..\..\Input\GAMMAPA\MEOH-ETOH-CYCLHEX.txt")
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
