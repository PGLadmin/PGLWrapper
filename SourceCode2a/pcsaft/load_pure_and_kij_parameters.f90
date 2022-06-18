!> \file load_pure_and_kij_parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE pcsaft_pure_and_binary_parameters
!!
!! Description, see below.
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module pcsaft_pure_and_binary_parameters

  use PARAMETERS, only: dp, machine_eps, str30
  use BASIC_VARIABLES, only: ncomp
  implicit none

  private
  save

  integer, parameter                                   :: str10 = 10

  !-----------------------------------------------------------------------------
  ! component name and molecular mass
  !-----------------------------------------------------------------------------
  character(LEN=str30), public, allocatable, dimension(:)  :: compna     ! names of substances
  character(LEN=str10), public, allocatable, dimension(:)  :: CAS_name   ! CAS-number of substances (character: 7732-18-5)
  integer, public, allocatable, dimension(:)           :: CAS_No         ! CAS-number of substances (integer: 7732185)
  real(dp), public, allocatable, dimension(:)          :: mm             ! molecular mass

  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------
  real(dp), public, allocatable, dimension(:)          :: mseg           ! number of segments
  real(dp), public, allocatable, dimension(:)          :: sigma          ! segment size parameter
  real(dp), public, allocatable, dimension(:)          :: epsilon_k      ! segment energy parameter
  real(dp), public, allocatable, dimension(:)          :: dipole_moment  ! dipole moment
  real(dp), public, allocatable, dimension(:)          :: quadru_moment  ! quadrupole moment
  
  integer, public, allocatable, dimension(:)           :: nhb_typ        ! number different types(!) of assoc. sites
  real(dp), public, allocatable, dimension(:,:)        :: nhb_no         ! number sites per type of assoc. sites
  real(dp), public, allocatable, dimension(:,:,:,:)    :: eps_hb         ! association energy parameter
  real(dp), public, allocatable, dimension(:,:)        :: kap_hb         ! association volume parameter

  character(LEN=20), public, allocatable, dimension(:) :: parameter_set  ! flag identifying a pure comp. param. set
  character(LEN=2), public, allocatable, dimension(:)  :: assoc_scheme   ! association in notation of Huang & Radosz

  !-----------------------------------------------------------------------------
  ! ideal gas heat capacity of pure substances
  !-----------------------------------------------------------------------------
  real(dp), public, allocatable, dimension(:,:)        :: cp_coefficients

  !-----------------------------------------------------------------------------
  ! binary parameters of fluid theory
  !-----------------------------------------------------------------------------
  real(dp), public, allocatable, dimension(:,:)        :: kij            ! regular binary param., for cross-energy param. epsilon
  real(dp), public, allocatable, dimension(:,:)        :: lij            ! asymmetric parameter, for cross-energy param. epsilon
  real(dp), public, allocatable, dimension(:,:)        :: kij_assoc      ! binary parameter, for cross-association energy epsilon_AB
  real(dp), public, allocatable, dimension(:,:)        :: kij_t          ! T-dependent correction to binary parameter kij

  !-----------------------------------------------------------------------------
  ! flags specifying active eos-contributions: association term, dipole-term, ...
  !-----------------------------------------------------------------------------

  logical, public                                      :: assoc
  logical, public                                      :: dipole
  logical, public                                      :: qudpole
  logical, public                                      :: dipole_quad
  logical, public                                      :: lij_correction
  logical, public                                      :: kij_T_dependent

  public :: load_pure_and_binary_parameters, file_open

contains



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE load_pure_and_kij_parameters
!> \brief  pure component parameters and kij parameters.
!!
!! This subroutine assigns pure component parameters and binary parameters.
!!
!! INPUT
!! -------
!! compna_input :      component name\n
!!
!! OUTPUT
!! -------
!! compna(i) :        component name\n
!! mm(i) :            molec. mass\n
!!
!! mseg(i) :          number of segments, [/] \n
!! sigma(i) :         segment diameter parameter, [Angstrom]\n
!! epsilon_k(i) :     segment energy param. epsilon/k, [K]\n
!! dipole_moment(i) : dipolar moment, [Ddebye]\n
!! quadru_moment(i) : quadrupolar moment, [Ddebye*Angstrom]\n
!! nhb_typ(i) :       number of types of association sites (integer), [/]\n
!! nhb_no(i,k) :      number of association sites of type k on molec. i (real), [/] \n
!! kap_hb(i,i) :      effective width of assoc. potential (angle-averg.), [/]\n
!! eps_hb(i,i,k,l) :  depth of association potential betw. site k and l, [K]\n
!!
!! assoc :            logical, .true. if an associating substance is considered
!! dipole :           logical, .true. if a dipolar substance is considered
!! qudpole :          logical, .true. if a quadrupolar substance is considered
!! dipole_quad :      logical, .true. if ( dipole .AND. qudpole )
!!
!! kij(i,j)           binary correction parameter acting on dispersion term, [/]\n
!! lij(i,j)           asymmetric binary correction para. (Tang and Gross, 2010), [/]\n
!! kij_assoc(i,j)     correction para. for association energy. [/]\n
!!                    It is not necessarily close to zero, because combining rules
!!                    for hydrogen-bonding interactions can be weak approximations.
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine load_pure_and_binary_parameters

  use BASIC_VARIABLES, only: compna_input

  !-----------------------------------------------------------------------------
  integer                                    :: i, j
  !-----------------------------------------------------------------------------

  compna( 1:ncomp ) = compna_input( 1:ncomp )

  call pcsaft_pure_parameters

  call pcsaft_binary_parameters

  !-----------------------------------------------------------------------------
  ! flags for association and polar terms. And flag if the second, asymmetric
  ! binary interaction parameter lij is used
  !-----------------------------------------------------------------------------
  assoc   = .false.
  qudpole = .false.
  dipole  = .false.
  do i = 1, ncomp
     if ( nhb_typ(i) /= 0 ) assoc = .true.
     if ( abs( quadru_moment(i) ) > machine_eps )   qudpole = .true.
     if ( abs( dipole_moment(i) ) > machine_eps )   dipole  = .true.
     do j = 1, ncomp
        if ( abs( lij(i,j) ) > machine_eps ) lij_correction = .true.
     end do
  end do
  if ( dipole .AND. qudpole )   dipole_quad = .true.

end subroutine load_pure_and_binary_parameters

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pcsaft_pure_parameters
!> \brief pcsaft paramaters
!!
!! pure component parameters from file
!! (as described in SUBROUTINE para_input)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pcsaft_pure_parameters

  use PARAMETERS, only: str30, nsite
  implicit none

  !-----------------------------------------------------------------------------
  integer                                    :: i
  character(LEN=str30)                       :: r_compna
  character(LEN=20)                          :: r_parameter_set
  real(dp)                                   :: r_mm
  real(dp)                                   :: r_mseg
  real(dp)                                   :: r_sigma
  real(dp)                                   :: r_epsilon_k
  real(dp)                                   :: r_dipole_moment
  real(dp)                                   :: r_quadru_moment
  integer                                    :: r_nhb_typ
  character(LEN=2)                           :: r_assoc_scheme
  real(dp), dimension(nsite)                 :: r_nhb_no
  real(dp)                                   :: r_kap_hb
  real(dp), dimension(nsite,nsite)           :: r_eps_hb
  character(LEN=str10)                       :: r_CAS_name

  character (LEN=300)                        :: entire_line
  character (LEN=100)                        :: rest_of_line
  integer                                    :: pos
  integer                                    :: nhb_no_var1, nhb_no_var2
  integer                                    :: read_info
  logical                                    :: parameter_assigned
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! all parameters are already initialized as zero (SR initialize_eos_quantities )
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! find pure-component parameters for each component
  !-----------------------------------------------------------------------------
  do  i = 1, ncomp

     !--------------------------------------------------------------------------
     ! Open parameter-input file and irgnore first line
     !--------------------------------------------------------------------------
     call file_open ( 53, './parameter_database/pcsaft_pure_parameters.txt' )

     parameter_assigned = .false.
     read_info = 0

     !--------------------------------------------------------------------------
     ! cycle line by line through pcsaft_pure_pararameters.txt
     !--------------------------------------------------------------------------
     do while ( read_info == 0 .AND. .NOT.parameter_assigned )

        !-----------------------------------------------------------------------
        ! CAS-Number, name, MM[g/mol], m[-], sigma[A], eps/k[K], kappaAB[-],
        ! epsAB/k[K], DM[Debye], Q[DA], NumberOfAssociationSites
        !-----------------------------------------------------------------------
        read ( 53, '(a)', IOSTAT = read_info) entire_line
        if ( read_info == 0 ) read (entire_line, *) r_CAS_name, r_compna,  &
             r_parameter_set, r_mm, r_mseg, r_sigma, r_epsilon_k,  &
             r_dipole_moment, r_quadru_moment, rest_of_line

        !-----------------------------------------------------------------------
        ! if substance found
        !-----------------------------------------------------------------------
        if ( trim(adjustl(compna(i))) == trim(adjustl(r_compna)) ) then

           rest_of_line = adjustl( rest_of_line )
           read( rest_of_line, '(i2,a)' ) r_nhb_typ, rest_of_line
           
           !--------------------------------------------------------------------
           ! assign component name and CAS number (as string and as integer)
           !--------------------------------------------------------------------
           compna(i) = r_compna

           CAS_name(i) = trim(adjustl( r_CAS_name ) )

           pos = scan( trim( r_CAS_name ), "-", BACK= .true. )
           if ( pos > 0 ) r_CAS_name = r_CAS_name(1:pos-1)//r_CAS_name(pos+1:str10)
           pos = scan( trim( r_CAS_name ), "-", BACK= .true. )
           if ( pos > 0 ) r_CAS_name = r_CAS_name(1:pos-1)//r_CAS_name(pos+1:str10)
           read( r_CAS_name, '(i9)' ) CAS_No(i)

           !--------------------------------------------------------------------
           ! assign pure component parameters
           !--------------------------------------------------------------------
           parameter_set(i) = r_parameter_set
           mm(i) = r_mm
           mseg(i) = r_mseg
           sigma(i) = r_sigma
           epsilon_k(i) = r_epsilon_k
           dipole_moment(i) = r_dipole_moment
           quadru_moment(i) = r_quadru_moment
           nhb_typ(i) = r_nhb_typ

           assoc_scheme(i) = ''
           nhb_no(i,:) = 0.0_dp
           kap_hb(i,i) = 0.0_dp
           eps_hb(i,i,:,:) = 0.0_dp

           !if ( r_nhb_typ == 0 ) write (*,'(a,x,a,x,6G14.7,i3)') compna(i),  &
           !     parameter_set(i), mm(i), mseg(i), sigma(i), epsilon_k(i),  &
           !     dipole_moment(i), quadru_moment(i), nhb_typ(i)

           !--------------------------------------------------------------------
           ! associating parameters, note kappaAB /= 0.0 defines an associating species
           !--------------------------------------------------------------------
           if ( r_nhb_typ > 0 ) then

              read (entire_line, *) r_CAS_name, r_compna, r_parameter_set, r_mm, r_mseg,  &
              r_sigma, r_epsilon_k, r_dipole_moment, r_quadru_moment, r_nhb_typ,  &
              r_assoc_scheme, rest_of_line

              !-----------------------------------------------------------------
              ! association scheme
              !-----------------------------------------------------------------
              if ( r_assoc_scheme == '2B' .or. r_assoc_scheme == '3B' .or. r_assoc_scheme == '4C' ) then

                 read (entire_line, *) r_CAS_name, r_compna, r_parameter_set, r_mm, r_mseg,  &
                 r_sigma, r_epsilon_k, r_dipole_moment, r_quadru_moment, r_nhb_typ,  &
                 r_assoc_scheme, nhb_no_var1, nhb_no_var2, r_kap_hb, r_eps_hb(1,2)

                 r_nhb_no(1) = real( nhb_no_var1, KIND=dp )
                 r_nhb_no(2) = real( nhb_no_var2, KIND=dp )

                 nhb_no(i,1) = r_nhb_no(1)
                 nhb_no(i,2) = r_nhb_no(2)
                 kap_hb(i,i) = r_kap_hb
                 eps_hb(i,i,1,2) = r_eps_hb(1,2)
                 eps_hb(i,i,2,1) = r_eps_hb(1,2)
                 eps_hb(i,i,1,1) = 0.0_dp
                 eps_hb(i,i,2,2) = 0.0_dp
                 assoc_scheme(i) = r_assoc_scheme

                 !write (*,'(a,x,a,x,6G14.7,i3,a,a,a,4G14.7)') compna(i), parameter_set(i),  &
                 !mm(i), mseg(i), sigma(i), epsilon_k(i), dipole_moment(i), quadru_moment(i),  &
                 !nhb_typ(i),' ',assoc_scheme(i),' ', nhb_no(i,1:2), kap_hb(i,i), eps_hb(i,i,1,2)

              else if ( r_assoc_scheme == '1A' ) then

                 read (entire_line, *) r_CAS_name, r_compna, r_parameter_set, r_mm, r_mseg,  &
                 r_sigma, r_epsilon_k, r_dipole_moment, r_quadru_moment, r_nhb_typ,  &
                 r_assoc_scheme, nhb_no_var1, r_kap_hb, r_eps_hb(1,2)

                 r_nhb_no(1) = real( nhb_no_var1, KIND=dp )

                 nhb_no(i,1) = r_nhb_no(1)
                 kap_hb(i,i) = r_kap_hb
                 eps_hb(i,i,1,1) = r_eps_hb(1,2)
                 assoc_scheme(i) = r_assoc_scheme

              end if

           end if

           !--------------------------------------------------------------------
           ! flag for exit of loop
           !--------------------------------------------------------------------
           parameter_assigned = .true.

        end if

        !-----------------------------------------------------------------------
        ! error message if EOF and component not found
        !-----------------------------------------------------------------------
        if ( read_info /= 0 .AND. .NOT.parameter_assigned )  then

           write (*,*) ' no PC-SAFT pure species parameters found for ', i,' ',trim(adjustl(compna(i)))
           stop

        end if

     end do

     close (53)

  end do

end subroutine pcsaft_pure_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pcsaft_par_kij
!> \brief  kij parameters
!!
!! draw binary interaction parameters from file
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pcsaft_binary_parameters

  implicit none

  !-----------------------------------------------------------------------------
  integer                                         :: i, j
  integer                                         :: read_info         ! Status-Flag fuer Einlesevorgang
  real(dp)                                        :: kij_var, lij_var, kij_assoc_var
  character(LEN=25)                               :: species_1, species_2
  character(LEN=10)                               :: CAS_No1, CAS_No2
  character(LEN=10)                               :: par_set_1, par_set_2
  character(LEN=150)                              :: entireline
  logical                                         :: parameter_assigned
  logical                                         :: verbose
  !-----------------------------------------------------------------------------

  verbose = .false.
  
  kij(:,:) = 0.0_dp
  lij(:,:) = 0.0_dp
  kij_assoc(:,:) = 0.0_dp

  do i = 1, ncomp

     do j = (i+1), ncomp

        !-----------------------------------------------------------------------
        ! Open parameter-input file and irgnore first line
        !-----------------------------------------------------------------------
        call file_open ( 54, './parameter_database/pcsaft_binary_parameters.txt' )

        parameter_assigned = .false.
        read_info = 0

        do while ( read_info == 0 .AND. .NOT.parameter_assigned )

           read( 54, '(A)', IOSTAT = read_info ) entireline
           read( entireline, * ) CAS_No1, CAS_No2, species_1, species_2,  &
                par_set_1, par_set_2, kij_var, lij_var, kij_assoc_var

           !--------------------------------------------------------------------
           ! Is a pair of species from file equal to the (ij)-pair?
           ! Two separate if-statements, because of non-symmetric
           ! lij(i,j) = -lij(j.i), so that the sequence of indices matters
           !--------------------------------------------------------------------
           if ( trim( adjustl( species_1 )) == trim( adjustl( compna(i) )) .AND. &
                trim( adjustl( species_2 )) == trim( adjustl( compna(j) )) ) then

              !if ( trim(adjustl( par_set_1 )) == trim( adjustl( parameter_set(i) )) .AND.  &
              !     trim(adjustl( par_set_2 )) == trim( adjustl( parameter_set(j) )) ) then

                 kij( i,j ) = kij_var
                 kij( j,i ) = kij_var
                 lij( i,j ) = lij_var
                 lij( j,i ) = - lij_var
                 kij_assoc( i,j ) = kij_assoc_var
                 kij_assoc( j,i ) = kij_assoc_var

                 parameter_assigned = .true.

              !end if

           else if ( trim( adjustl( species_2 )) == trim( adjustl( compna(i) )) .AND. &
                     trim( adjustl( species_1 )) == trim( adjustl( compna(j) )) ) then

              !if ( trim(adjustl( par_set_2 )) == trim( adjustl( parameter_set(i) )) .AND.  &
              !     trim(adjustl( par_set_1 )) == trim( adjustl( parameter_set(j) )) ) then

                 kij( i,j ) = kij_var
                 kij( j,i ) = kij_var
                 lij( i,j ) = - lij_var
                 lij( j,i ) = lij_var
                 kij_assoc( i,j ) = kij_assoc_var
                 kij_assoc( j,i ) = kij_assoc_var

                 parameter_assigned = .true.

              !end if

           end if

           if ( parameter_assigned .AND. verbose ) write (*,'(a,2i4,2x,a,x,a,2x,3G12.5)')  &
                 ' assign kij parameter  :', i, j, trim( adjustl( compna(i) )),  &
                 trim( adjustl( compna(j) )), kij(i,j), lij(i,j), kij_assoc(i,j)

           !--------------------------------------------------------------------
           ! if end-of-file is reached without loading binary parameters
           !--------------------------------------------------------------------
           if ( read_info /= 0 .AND. .NOT.parameter_assigned .AND. verbose )  then

              write (*,'(a,2i4,x,a,x,a,x,3G12.5)') ' no kij parameter found:', i, j,  &
                 trim(adjustl(compna(i))), trim(adjustl(compna(j))), kij(i,j), lij(i,j), kij_assoc(i,j)

           end if

        end do

        close (54)

     end do

  end do


end subroutine pcsaft_binary_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine file_open
!!
!! This subroutine opens files for reading. \n
!! Beforehand, it checks whether this file is available.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine file_open ( file_number, filename )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: file_number
  character (LEN=*), intent(in)                   :: filename
  !-----------------------------------------------------------------------------

  logical                                         :: filefound
  !-----------------------------------------------------------------------------

  inquire (FILE=filename, EXIST = filefound)
  if (filefound) then
     open (file_number, FILE = filename)
  else
     write (*,*) ' FOLLOWING FILE CAN NOT BE OPENED', filename
     stop
  end if

end subroutine file_open

end module pcsaft_pure_and_binary_parameters

