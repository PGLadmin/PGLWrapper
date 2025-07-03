module findsites
    use sitenspecies
    integer, save :: aspmx
    type (species), allocatable, save  :: compa(:) ! components present in simulation file
contains
subroutine loadsites(T,n,x,idx,kop, nsitesp, KAD, eps, compp, sitep, aspmx)
implicit none

!************Declare variables
integer, intent(IN):: n, idx(n),kop(10) ! n: number of components present in the function call
real*8, intent(IN):: x(n) ! n: number of components present in the function call
integer, intent(OUT):: nsitesp, aspmx ! nsitesp is packed indexes in function call
real*8, dimension(:,:), allocatable, intent(OUT) :: KAD, eps
integer, dimension(:), allocatable, save :: ids_pack, ids_call ! for sites present in call, file ids and call host index
type (species), intent(OUT)  :: compp(n) ! components present in function call
type (siteinfo), dimension(:), allocatable, save :: site(:) ! sites in simulation file
type (siteinfo), dimension(:), allocatable, intent(OUT) :: sitep(:) ! dimensioned later to be the number of sites present
real*8 RMISS
integer IMISS, ISET, NELEM, NBG, LBG, NG, I, J, K, L, M, IDATIJ, JGRP, LBGV, LOFF, LOFFK, LUGMK, LOFFE, LUGME, LUFGRP, LGRP, LGRPE, nc
integer, save :: nsites, acptflag, acptflagmx, cflag, cflagmx, dnrflag, dnrflagmx, countcall
integer name1, name2, name(2)
real*8 T, term1
real*8, dimension(:,:), allocatable, save :: Ksave, Esave ! copy of Plex values of KAD and eps
character(len=8) name_grp1, name_grp2

#include "ppexec_user.cmn"
#include "dms_maxwrt.cmn" ! for writing to control panel

      EQUIVALENCE (RMISS, USER_RUMISS)
      EQUIVALENCE (IMISS, USER_IUMISS)
!
#include "pputl_ppglob.cmn"
#include "dms_global.cmn"
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

INTEGER DMS_LOCATI, LBPROC, DMS_IFCMNC, LIDSCC

!************End of declarations
data countcall /0/
! count function calls to identify initial and subsequent calls
countcall = countcall + 1
nc = NCOMP_NCC

LUFGRP = DMS_IFCMNC('UFGRP') ! UNIFAC group used for site information
IF (kop(2) .NE. 3) THEN
LUGMK = DMS_IFCMNC('KAD') ! Wertheim binary parameter K
LUGME = DMS_IFCMNC('EPS') ! Wertheim binary parameter eps
ELSE
LUGMK = DMS_IFCMNC('KADPCS') ! PCSAFT Wertheim binary parameter K
LUGME = DMS_IFCMNC('EPSPCS') ! PCSAFT Wertheim binary parameter eps
ENDIF
LIDSCC = DMS_IFCMNC('IDSCC') !component name
IF (PPCTBL_NBGV .GT. 0) LBGV = DMS_LOCATI(PPCTBL_NBGV)

! ********************************** Wertheim contribution ***********************
! Count number of groups in function call
! All subroutine vectors use conventional components only, where N is the number of
! conventional components. IDX is used to get the index from the simulation

if(countcall.eq.1) THEN ! load sites and parameters the first time

! find all user groups in simulation file (some may not be present in subroutine call
LBPROC = DMS_LOCATI(2)
NBG = IB(LBPROC + 131) ! listing of groups in simulation file
IF (NBG .NE. 0) THEN
   LBG = DMS_LOCATI(NBG)
   NG = IB(LBG) ! aspen number of groups
   ! each group is stored with a repeat unit of 6
   ! first element - Group ID - first word (4 letters) of group name
   ! second element - Group ID - second word (4 letters) of group name
   ! third element - blank
   ! fouth element - Group number (e.g. 1010, 4500, 4700) will be >= 4500 for user group
   ! fifth element - 0
   ! sixth element - group sequence number (order entered in user form) for user-defined groups
   ! Group IDs and group sequence numbers are set to zero for built-in groups
   DO J = 1, NG ! loop through groups to find user groups used for sites
    LOFF = LBG + 1 + 6*(J-1)
    k = IB(LOFF + 5) ! check group sequence number, this will increment for user groups, but is zero for built-in groups, user groups are listed first
    if (k .eq. 0) exit ! done with user groups when 0 is found
   ENDDO
ENDIF

nsites = J - 1 ! number of UNIFAC user groups (used for sites) in simulation file
if(.not.allocated(site)) allocate(site(nsites)) ! allocate storage for all sites in simulation file
if(.not.allocated(ids_pack)) allocate(ids_pack(nsites))
if(.not.allocated(ids_call)) allocate(ids_call(nsites))
if(.not.allocated(Ksave)) allocate(Ksave(nsites,nsites)) ! for storage of values between function calls
if(.not.allocated(Esave)) allocate(Esave(nsites,nsites)) ! for storage of values between function calls
! note that KAD and eps are (nsitesp,nsitesp) while the saved values are (nsites,nsites)

! store id values for user groups, separate loop because nsites is not known
do j = 1, nsites
    LOFF = LBG + 1 + 6*(J-1)
    site(j)%idaspen = IB(LOFF+3) ! store sites IDs in the order entered in the Group form
    write(site(j)%name,'(2A4)') IB(LOFF), IB(LOFF+1)
enddo

ISET = 1 ! Currently written for one parameter set
NELEM = 1 ! One parameter element in each set for each site

! load del matrix with information for all sites in file. Uncomment conditional statements for serious debugging memory pointers
! if ((kop(1).gt.2)) write(global_nh,'site index, site index, Aspen site index, Aspen site index, memory index') ! for MSU debugging
do l = 1, nsites
    do k = 1, nsites
        IDATIJ = NELEM*(k+nsites*(l-1)+nsites*nsites*(ISET-1)-1)
        loffk = lugmk + idatij
        loffe = lugme + idatij
!        if ((kop(1).gt.2)) write(global_nh,'(2x,i3,2x,i3,2x,i3,2x,i3,2x,i3)') k,l,i,j,idatij+1  ! for MSU debugging
        Ksave(k,l) =  B(loffk+1) ! stores the 1st element of the parameter for the i-j pair.
        if (Ksave(k,l) .eq. RMISS) Ksave(k,l) = 0D0 ! for parameters not entered by user, set value to zero
!        if ( k .eq. l ) Ksave(k,l) = 2D0 * Ksave(k,l) ! This was incorrect, CTL 11/24/2021
        Esave(k,l) = B(loffe+1)
        if (Esave(k,l) .eq. RMISS) Esave(k,l) = 0D0
    enddo ! j
enddo ! l

! user-entered vector has only half the matrix, now fill zeros
do j = 1, nsites
    do i = j, nsites
    if (Ksave(i,j) .eq. 0.D0) Ksave(i,j) = Ksave(j,i)
    if (Ksave(j,i) .eq. 0.D0) Ksave(j,i) = Ksave(i,j)
    if (Esave(i,j) .eq. 0.D0) Esave(i,j) = Esave(j,i)
    if (Esave(j,i) .eq. 0.D0) Esave(j,i) = Esave(i,j)
    enddo ! i
enddo ! j

if(.not.allocated(compa)) allocate(compa(nc)) ! allocate memory for file component information

do i = 1, nc  ! loop over conventional components in simulation file to find site hosts.
! program stucture assumes that each site appears in only on species
! flags to track if component is self associating, nonzero if site type if found in component
     dnrflag = 0
     acptflag = 0
     cflag = 0
	 LGRP = (i-1)*24 + 1 ! find index of first group for component - up to 24 group entries for each molecule
	 LGRPE = LGRP + 24 ! ending index for group information for the component
!     temp = B(lufgrp+25:(lufgrp+50+25))   ! for debugging
	do j = LGRP, LGRPE, 2 ! stored as group number and number of occurrences, increment by two
	   ! identify hosts of sites
	   term1 = B(LUFGRP + j)
	   if (term1 .eq. RMISS) term1 = 0.D0
	   JGRP = term1 ! aspen site number stored as real, converted to integer.
	   ! aspen site number is zero (user groups) or RMISS (built-in groups) at end of component structure.
	   ! RMISS direct conversion gives a negative integer but also throws an error in aspen. Thus convert to zero first.
	   if (JGRP .lt. 1) then
		! Before exiting, check to see if species self-associates, only important if more thanone
		! one site is present (loop must be done more than once).
		if ( (cflag + dnrflag + acptflag) .gt. 1 ) compa(i)%selfassoc = 1
	     EXIT
	    end if
	    if (JGRP .LT. 4500) cycle ! built-in group rather than user group
	    ! find the site matches the site number with the index number, k
        do k = 1, nsites
            if (jgrp .eq. site(k)%idaspen) exit
        end do
		! comp is a structure of type 'species', size NCOMP_NCC for the components present in
		! the simulation file.
		! compa(i)%nsite - number of sites on the component
		! compa(i)%sid - vector of side ids on the component
	   compa(i)%nsite = compa(i)%nsite + 1 ! increment number of sites found on host
	   compa(i)%sid(compa(i)%nsite) = k ! store site id
	   ! site is a structure of type 'siteinfo' size nsites (all sites in simulation file)
	   !site(k)%idaspen = jgrp ! store aspen simulation file group id
       site(k)%host = i ! host component simulation index
       site(k)%noccur = B(LUFGRP+j+1) ! number of occurences of the site
!      site(k)%xhost = x(i) ! mole fraction of the host is not used. This field is set below
!	   site(k)%nhost = site(k)%nhost + 1    ! this will permit sites to appear on more than one host in the future
!	   site(k)%host(site(k)%nhost)=i        ! this will permit sites to appear on more than one host in the future
!	   site(k)%noccur(site(k)%nhost) = B(LUFGRP+j+1) ! this will permit sites to appear on more than one host in the future
!		Determine if group is a C-type or donor or acceptor.
		if (JGRP.GT. 4949) THEN ! bivalent site present
			cflag = 2
		else if( JGRP .LT. 4750) THEN ! donor present
			dnrflag = 1
			compa(i)%dnr = 1
		else ! acceptor present
			acptflag = 1
			compa(i)%acpt = 1
		end if
    enddo
enddo

endif ! countcall = 1

if ((kop(1).gt.1)) then
    write(global_nh,*) 'User Sites in Aspen file (subset may be used), count, name, idaspen, file host'
    do j=1,nsites
    write(global_nh,'(I3,2X,A8,2X,I4,2X,I4)') j, site(j)%name, site(j)%idaspen, site(j)%host ! Write all group names in simulation file
    enddo
end if

! now find the sites used only in this function call
! flags for mixture, nonzero if that type of site is found in mixuture
dnrflagmx = 0
acptflagmx = 0
cflagmx = 0
aspmx = 0 ! =1 if Association or Solvation is Present in mixture

do i = 1, n  ! loop over conventional components present to find site hosts.
! program stucture assumes that each site appears in only on species
     m = idx(i) ! look up only components used in this function call, aspen component table index
     compp(i) = compa(m) ! copy species information
     compp(i)%x = x(i)
		if (compp(i)%selfassoc) THEN
			cflagmx = 2
		else if(compp(i)%dnr) THEN
			dnrflagmx = 1
		else
			acptflagmx = 1
		end if
enddo

if((kop(1).gt.1)) write(global_nh,*) 'Site in call, file site, siteid, name, file hostid, call hostid, hostx'

! count sites present
! todo - counting in the loop over n would work only if the sites appear on only one species
nsitesp = 0
do i = 1,nsites
	! ids_pack is a vector of length nsites that points to the groups present in the function call.
	! If nsites = 5 and nsitesp = 3 and the present sites are 1, 4, 5
	! after this loop, ids_pack will be [1,4,5,0,0], ids_call will be the call index of the host
    ids_pack(i) = 0 ! set to zero if site is not present in function call
    ids_call(i) = 0
    do j = 1,n
    if(site(i)%host .eq. idx(j)) then ! site is present in simulation
        nsitesp = nsitesp+1
        ids_pack(nsitesp) = i ! point to the packed site index if site is present in function call
        ids_call(nsitesp) = j
        site(i)%packed = nsitesp ! stored the packed index for later cross reference
        site(i)%xhost = x(j)
        if((kop(1).gt.1)) write(global_nh,'(I3,2X,I3,2X,I4,2X,A8,2(2X,I3),2X,F10.5)')nsitesp, i, site(i)%idaspen, site(i)%name, site(i)%host, j, site(i)%xhost
    endif
    enddo !j
enddo !i

! reset the compp%sid() to be the packed id identifier
do i=1,n
    do j=1,compp(i)%nsite
        compp(i)%sid(j)=site(compp(i)%sid(j))%packed
    enddo !j
enddo !i

! ---------------decide if association occurs at all in mixture------------------
! if there are not acceptors + donors, acceptor + C, or donors + C, then there is no wertheim component
! and all association calculations are skipped, and the wertheim contribution is zero.

if( (cflagmx + dnrflagmx + acptflagmx) .gt. 1 ) then

aspmx = 1

! now build del matrix only for the sites that are present in the function call
allocate(KAD(nsitesp,nsitesp),eps(nsitesp,nsitesp),sitep(nsitesp))

! load del matrix with information for only the sites present. Uncomment conditional statements for serious debugging memory pointers
! if ((kop(1).gt.2)) write(global_nh,'site index, site index, Aspen site index, Aspen site index, memory index') ! for MSU debugging
do l = 1, nsitesp
    j = ids_pack(l) ! index for site that is present
    !transfer information for sites present
    sitep(l) = site(j)
    sitep(l)%host = ids_call(l) ! host index for function call
    do k = l, nsitesp
        i = ids_pack(k) !index for site that is present
        KAD(k,l) = Ksave(i,j) ! stores the 1st element of the parameter for the i-j pair.
        KAD(l,k) = KAD(k,l)
        eps(k,l) = Esave(i,j)
        eps(l,k) = eps(k,l)
    enddo ! j
enddo ! l

if((kop(1).gt.1)) then
! extra debugging information
write(global_nh,*)'Site1, Site2, index1, index2, Keps, eps, kad'
do l = 1, nsitesp
    ! for debugging examine name_grp1
    j = ids_pack(l) ! index for site that is present
    do k = 1, nsitesp
        ! for debugging examine name_grp2
        i = ids_pack(k) !index for site that is present
        write(global_nh,'(2A10,2I4, 3G12.5)') site(i)%name,site(j)%name,k,l,kad(k,l)*(dexp(eps(k,l)/T)-1D0), eps(k,l), kad(k,l)
    enddo !k
enddo !l
endif

endif ! if asp

END SUBROUTINE LOADSITES
end module findsites
