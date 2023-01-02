    
!> \file mod_basic_variables.f90
!! \brief Several modules with all basic variables.
!!
!! \todo Detailed file description.

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module PARAMETERS
!> \brief Module with important global constants and parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module PARAMETERS	! JRE: This is an unfortunate choice of module name because many programs may have their own "parameters."  These "PARAMETERS" are fairly universal, so let's just make due and hope to avoid collisions.

  implicit none
  save

  integer, parameter                            :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
  real(dp), parameter                           :: machine_eps = epsilon( 1.0_dp )

  integer, parameter                            :: nsite = 5                 !< Max. number of assoc. sites

  real(dp), parameter                           :: PI = 3.141592653589793_dp !< Def. of Pi
  real(dp), parameter                           :: PI_6 = PI / 6.0_dp        !< Def. of Pi/6
  real(dp), parameter                           :: NAv  = 6.022140857E23_dp  !< Def. of Avogadro's number
  real(dp), parameter                           :: KBOL = 1.38064852E-23_dp  !< Def. of Boltzmann's const. in J/K
  real(dp), parameter                           :: RGAS = KBOL * NAv         !< Def. of universal gas const.
  real(dp), parameter                           :: KBOL30 = KBOL * 1.E30_dp  !< Def. of Boltzmann's const. in Pa/(K*Angstrom**3)

  integer, parameter                            :: str30 = 30


End Module PARAMETERS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module BASIC_VARIABLES
!> \brief This module contains basic global variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module BASIC_VARIABLES

  use PARAMETERS, only: dp, nsite, str30
  implicit none
  save

  !-----------------------------------------------------------------------------
  ! basic quantities defining the mixture
  !-----------------------------------------------------------------------------
  integer                                         :: ncomp         ! number of substances

  real(dp)                                        :: t_input       ! temperature in K
  real(dp)                                        :: p_input       ! pressure in Pa
  real(dp), allocatable, dimension(:)             :: x_input
  character(LEN=str30), allocatable, dimension(:) :: compna_input  ! names of substances

End Module BASIC_VARIABLES

subroutine consoleIO  !JRE simplification

  use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
  implicit DoublePrecision(a-h,o-z)
  character*1 answer
  ncomp=2
  if( .not.allocated(compna_input) )allocate( compna_input(ncomp) )
  if( .not.allocated(x_input) )allocate( x_input(ncomp) )
  compna_input(1)='carbon-dioxide_a'
  compna_input(2)='ethanol_n'
  write(*,*)' Current components are: ',(  TRIM( compna_input(i) ),i=1,ncomp  )
  print*,'Do you want to change the components? (y/n)'
  read*,answer
  if(answer=='y'.or. answer=='Y')then
      print*,'Enter nComp'
      read*,nComp
      do i=1,nComp
          print*,'Enter comp name for component',i
          read*,compna_input(i)
      enddo
  endif
  print*,'Enter T(K),P(MPa)'
  read*,t_input,p_input
  p_input=p_input*1.D6 !Internal units for Gross are Pa.
  if(ncomp==1)then
    x_input(1)=1
    return
  elseif(ncomp==2)then
      print*,'Enter x1'
      read*,x_input(1)
      x_input(2)=1-x_input(1)
  else
      print*,'Enter x1-xn'
      read*,( x_input(i), i=1,ncomp )
  endif
  
  print*,'consoleIO Done!'

end subroutine consoleIO

    
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
  real(dp), public, allocatable, dimension(:)          :: sigma          ! segment size parameter (Angstrom)
  real(dp), public, allocatable, dimension(:)          :: epsilon_k      ! segment energy parameter (K)
  real(dp), public, allocatable, dimension(:)          :: dipole_moment  ! dipole moment (Debye)
  real(dp), public, allocatable, dimension(:)          :: quadru_moment  ! quadrupole moment (Debye-Angstrom)
  
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
  use GlobConst, only: LOUD, bVolCc_mol, isPcSaft, iEosOpt !JRE LOUD=.FALSE. to silence I/O feedback.

  !-----------------------------------------------------------------------------
  integer                                    :: i, j
  LOGICAL									LOUDER	  !jre to help debugging
  !-----------------------------------------------------------------------------
  LOUDER=LOUD	 !jre to help debugging
  !LOUDER=.TRUE.	 !jre to help debugging
  compna( 1:ncomp ) = compna_input( 1:ncomp )

  call pcsaft_pure_parameters
  if(LOUDER)print*,'load_pure...: Pure OK. Trying binary.'

  call pcsaft_binary_parameters
  if(LOUDER)print*,'load_pure...: Binary OK. Setting HB and polar flags.'
  if(LOUDER)print*,'comp     1      2 from load_pure_and_binary_parameters'
  do i=1,ncomp
	IF(LOUDER)write(*,'(i3,11f8.4)')i,(kij(i,j),j=1,ncomp)
  enddo

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
     call file_open ( 53, './Input/PcSaft_database/pcsaft_pure_parameters.txt' )   !JRE

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
        if ( read_info == 0 ) read (entire_line, *, IOSTAT = read_info) r_CAS_name, r_compna,  &	   ! JRE20210524, added ", read_info)"
             r_parameter_set, r_mm, r_mseg, r_sigma, r_epsilon_k,  &
             r_dipole_moment, r_quadru_moment, rest_of_line
			 if(read_info /= 0)then
			   write (*,*) ' pcsaft_pure_parameters: bad line in pcsaft_pure_parameters.txt. Offender is:'
			   write (*,*) TRIM(entire_line)
			   pause 'Possible offense: comma(s) in compound name.' 
               stop
			 endif


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

     end do	 ! while ( read_info == 0 .AND. .NOT.parameter_assigned )

     close (53)

  end do ! i=1,nComp

end subroutine pcsaft_pure_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pcsaft_par_kij
!> \brief  kij parameters
!!
!! draw binary interaction parameters from file
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pcsaft_binary_parameters

  use GlobConst, only: LOUD, bVolCc_mol, isPcSaft, iEosOpt, ID !JRE LOUD=.FALSE. to silence I/O feedback.
  implicit none

  !-----------------------------------------------------------------------------
  integer                                         :: i, j, line, id1, id2, idCas1,idCas2
  integer                                         :: read_info,readLineErr         ! Status-Flag fuer Einlesevorgang
  real(dp)                                        :: kij_var, lij_var, kij_assoc_var
  character(LEN=25)                               :: species_1, species_2
  !character(LEN=10)                               :: CAS_No1, CAS_No2
  character(LEN=10)                               :: par_set_1, par_set_2
  character(LEN=150)                              :: entireline
  logical                                         :: parameter_assigned
  logical                                         :: verbose
  logical                                         :: LOUDER
  !-----------------------------------------------------------------------------

  verbose = .false.
  LOUDER=LOUD
  LOUDER=.FALSE.
  !LOUDER=.TRUE.
  verbose = LOUDER
  
  kij(:,:) = 0.0_dp
  lij(:,:) = 0.0_dp
  kij_assoc(:,:) = 0.0_dp

  do i = 1, ncomp

     do j = (i+1), ncomp

        !-----------------------------------------------------------------------
        ! Open parameter-input file and irgnore first line
        !-----------------------------------------------------------------------
		if(LOUDER)print*,'pcsaft_binary... opening ./Input/PcSaft_database/pcsaft_binary_parametersJre.txt'  !JRE
        call file_open ( 54, './Input/PcSaft_database/pcsaft_binary_parametersJre.txt' )				   !JRE

        parameter_assigned = .false.
        read_info = 0
		line=0
        do while ( read_info == 0 .AND. .NOT.parameter_assigned )
		   line=line+1
           read( 54, '(A)', IOSTAT = read_info ) entireline
			!if(LOUDER)print*,'pcsaft_binary... entireline=',TRIM(entireline)
           read( entireline, *,IOSTAT=readLineErr ) id1,id2,idCas1, idCas2, kij_var,lij_var,kij_assoc_var,par_set_1,par_set_2
           !JRE read( entireline, *,IOSTAT=readLineErr ) CAS_No1, CAS_No2, species_1, species_2,  &
           !     par_set_1, par_set_2, kij_var, lij_var, kij_assoc_var
			if(readLineErr /= 0)then
				if(LOUDER)print*,'pcsaft_binary... readLineErr for line=',line,' a ', TRIM(entireline)
				if(LOUDER)pause 'Check which line is causing the problem. Maybe you can fix it.'
				cycle
			else
				if(LOUDER)write(*,'(a,2i5,f8.4)')'id1,id2,kij_var',id1,id2,kij_var
			endif 

           !--------------------------------------------------------------------
           ! Is a pair of species from file equal to the (ij)-pair?
           ! Two separate if-statements, because of non-symmetric
           ! lij(i,j) = -lij(j.i), so that the sequence of indices matters
           !--------------------------------------------------------------------
           !JRE if ( trim( adjustl( species_1 )) == trim( adjustl( compna(i) )) .AND. &
           !JRE     trim( adjustl( species_2 )) == trim( adjustl( compna(j) )) ) then
		   if( id1==MAX(ID(1),ID(2)) .and. id2==MIN(ID(1),ID(2)) )then  !JRE Wrapper is based on DIPPR id's.

              !if ( trim(adjustl( par_set_1 )) == trim( adjustl( parameter_set(i) )) .AND.  &
              !     trim(adjustl( par_set_2 )) == trim( adjustl( parameter_set(j) )) ) then

                 kij( i,j ) = kij_var
                 kij( j,i ) = kij_var
                 lij( i,j ) = lij_var
                 lij( j,i ) = - lij_var
                 kij_assoc( i,j ) = kij_assoc_var
                 kij_assoc( j,i ) = kij_assoc_var

                 parameter_assigned = .true.
				 !pause 'iGotit !'

              !end if

           else if ( trim( adjustl( species_2 )) == trim( adjustl( compna(i) )) .AND. &
                     trim( adjustl( species_1 )) == trim( adjustl( compna(j) )) ) then

              !if ( trim(adjustl( par_set_2 )) == trim( adjustl( parameter_set(i) )) .AND.  &
              !     trim(adjustl( par_set_1 )) == trim( adjustl( parameter_set(j) )) ) then
                 !kij( i,j ) = kij_var	 !JRE: disabled this to ensure it doesn't override above
                 !kij( j,i ) = kij_var
                 lij( i,j ) = - lij_var
                 lij( j,i ) = lij_var
                 kij_assoc( i,j ) = kij_assoc_var
                 kij_assoc( j,i ) = kij_assoc_var

                 parameter_assigned = .true.

              !end if

           end if

           if ( parameter_assigned .AND. LOUDER ) write (*,'(a,2i4,2x,a,x,a,2x,3G12.5)')  &
                 ' assign kij parameter  :', i, j, trim( adjustl( compna(i) )),  &
                 trim( adjustl( compna(j) )), kij(i,j), lij(i,j), kij_assoc(i,j)
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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module EOS_CONSTANTS
!> \brief  Module contains EOS-constants
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module eos_constants

  use PARAMETERS, only: dp
  implicit none
  save

  private
  public :: load_eos_constants, load_dd_constants, load_qq_constants, load_dq_constants

  real(dp), public, dimension( 0:6, 3 )             :: ap, bp

  real(dp), public, allocatable, dimension(:,:,:)   :: qqp2, qqp4, ddp2, ddp4, dqp2, dqp4
  real(dp), public, allocatable, dimension(:,:,:,:) :: qqp3, ddp3, dqp3

contains


!> \file eos_const.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine eos_const
!
!> \brief provides the constants of the PC-SAFT EOS
!
! This subroutine provides the constants of the PC-SAFT EOS.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_eos_constants !( ap, bp )

  use pcsaft_pure_and_binary_parameters, only: dipole, qudpole, dipole_quad
  !-----------------------------------------------------------------------------

  call load_pc_saft_constants

  if ( dipole ) call load_dd_constants !( ddp2, ddp3, ddp4 )

  if ( qudpole ) call load_qq_constants !( qqp2, qqp3, qqp4 )

  if ( dipole_quad ) call load_dq_constants !( dqp2, dqp3, dqp4 )
  
end subroutine load_eos_constants



!> \file eos_const.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine eos_const
!
!> \brief provides the constants of the PC-SAFT EOS
!
! This subroutine provides the constants of the PC-SAFT EOS.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_pc_saft_constants

  !-----------------------------------------------------------------------------
  ! constants of the dispersion term
  !-----------------------------------------------------------------------------

  ap(0,1) =  0.91056314451539_dp
  ap(0,2) = -0.30840169182720_dp
  ap(0,3) = -0.09061483509767_dp
  ap(1,1) =  0.63612814494991_dp
  ap(1,2) =  0.18605311591713_dp
  ap(1,3) =  0.45278428063920_dp
  ap(2,1) =  2.68613478913903_dp
  ap(2,2) = -2.50300472586548_dp
  ap(2,3) =  0.59627007280101_dp
  ap(3,1) = -26.5473624914884_dp
  ap(3,2) =  21.4197936296668_dp
  ap(3,3) = -1.72418291311787_dp
  ap(4,1) =  97.7592087835073_dp
  ap(4,2) = -65.2558853303492_dp
  ap(4,3) = -4.13021125311661_dp
  ap(5,1) = -159.591540865600_dp
  ap(5,2) =  83.3186804808856_dp
  ap(5,3) =  13.7766318697211_dp
  ap(6,1) =  91.2977740839123_dp
  ap(6,2) = -33.7469229297323_dp
  ap(6,3) = -8.67284703679646_dp

  bp(0,1) =  0.72409469413165_dp
  bp(0,2) = -0.57554980753450_dp
  bp(0,3) =  0.09768831158356_dp
  bp(1,1) =  2.23827918609380_dp
  bp(1,2) =  0.69950955214436_dp
  bp(1,3) = -0.25575749816100_dp
  bp(2,1) = -4.00258494846342_dp
  bp(2,2) =  3.89256733895307_dp
  bp(2,3) = -9.15585615297321_dp
  bp(3,1) = -21.0035768148465_dp
  bp(3,2) = -17.2154716477721_dp
  bp(3,3) =  20.6420759743972_dp
  bp(4,1) =  26.8556413626615_dp
  bp(4,2) =  192.672264465249_dp
  bp(4,3) = -38.8044300520628_dp
  bp(5,1) =  206.551338406619_dp
  bp(5,2) = -161.826461648765_dp
  bp(5,3) =  93.6267740770146_dp
  bp(6,1) = -355.602356122079_dp
  bp(6,2) = -165.207693455561_dp
  bp(6,3) = -29.6669055851473_dp

end subroutine load_pc_saft_constants



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine dq_const
!
!> \brief provides the constants of the dipole-quadrupole term
!
! This subr. provides the constants of the dipole-quadrupole term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_dq_constants

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mdq_i, mdq_j, mdq_k
  real(dp)                                   :: mf1, mf2, msegij
  !-----------------------------------------------------------------------------

  do i = 1, ncomp

     mdq_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mdq_j = min( mseg(j), 2.0_dp )

        msegij = ( mdq_i * mdq_j )**0.5
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        dqp2(i,j,0) =  0.697094963_dp + mf1*(-0.673459279_dp) + mf2*0.670340770_dp
        dqp2(i,j,1) = -0.633554144_dp + mf1*(-1.425899106_dp) + mf2*(-4.338471826_dp)
        dqp2(i,j,2) =  2.945509028_dp + mf1 * 4.19441392_dp   + mf2*7.234168360_dp
        dqp2(i,j,3) = -1.467027314_dp + mf1 * 1.0266216_dp
        dqp2(i,j,4) = 0.0_dp

        dqp4(i,j,0) = -0.484038322_dp + mf1 * 0.67651011_dp   + mf2*(-1.167560146_dp)
        dqp4(i,j,1) =  1.970405465_dp + mf1*(-3.013867512_dp) + mf2*2.13488432_dp
        dqp4(i,j,2) = -2.118572671_dp + mf1 * 0.46742656_dp
        dqp4(i,j,3) = 0.0_dp
        dqp4(i,j,4) = 0.0_dp


        do k = 1, ncomp
           mdq_k = min( mseg(k), 2.0_dp )
           msegij = ( mdq_i * mdq_j * mdq_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = ( msegij - 2.0_dp ) / msegij
           dqp3(i,j,k,0) = 0.795009692_dp + mf1*(-2.099579397_dp)
           dqp3(i,j,k,1) = 3.386863396_dp + mf1*(-5.941376392_dp)
           dqp3(i,j,k,2) = 0.475106328_dp + mf1*(-0.178820384_dp)
           dqp3(i,j,k,3) = 0.0_dp
           dqp3(i,j,k,4) = 0.0_dp
        end do

     end do

  end do

end subroutine load_dq_constants


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine dd_const
!
!> \brief provides the constants of the dipole-term
! This subroutine provides the constants of the dipole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_dd_constants

  use PARAMETERS, only: PI
  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mdd_i, mdd_j, mdd_k
  real(dp)                                   :: mf1, mf2, msegij, sin2t
  !-----------------------------------------------------------------------------

  sin2t = SIN( 0.0_dp * PI / 180.0_dp )
  sin2t = sin2t*sin2t

  do i = 1, ncomp

     mdd_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mdd_j = min( mseg(j), 2.0_dp )

        msegij = ( mdd_i * mdd_j )**0.5_dp
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        ddp2(i,j,0) =  0.30435038064_dp + mf1*(0.95346405973_dp+0.201436_dp*sin2t)  &
                      + mf2*(-1.16100802773_dp-1.74114_dp*sin2t)
        ddp2(i,j,1) = -0.13585877707_dp + mf1*(-1.83963831920_dp+1.31649_dp*sin2t)  &
                      + mf2*4.52586067320_dp
        ddp2(i,j,2) =  1.44933285154_dp + mf1 * 2.01311801180_dp  + mf2*0.97512223853_dp
        ddp2(i,j,3) =  0.35569769252_dp + mf1*(-7.37249576667_dp) + mf2*(-12.2810377713_dp)
        ddp2(i,j,4) = -2.06533084541_dp + mf1 * 8.23741345333_dp  + mf2*5.93975747420_dp

        ddp4(i,j,0) =  0.21879385627_dp + mf1*(-0.58731641193_dp) + mf2*3.48695755800_dp
        ddp4(i,j,1) = -1.18964307357_dp + mf1 * 1.24891317047_dp  + mf2*(-14.9159739347_dp)
        ddp4(i,j,2) =  1.16268885692_dp + mf1*(-0.50852797392_dp) + mf2*15.3720218600_dp
        ddp4(i,j,3) =  0.0_dp
        ddp4(i,j,4) =  0.0_dp

        do k = 1, ncomp

           mdd_k = min( mseg(k), 2.0_dp )

           msegij = ( mdd_i * mdd_j * mdd_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = mf1 * ( msegij - 2.0_dp ) / msegij
           ddp3(i,j,k,0) = -0.06467735252_dp + mf1*(-0.95208758351_dp+0.28503_dp*sin2t)  &
                + mf2*(-0.62609792333_dp+2.2195*sin2t)
           ddp3(i,j,k,1) =  0.19758818347_dp + mf1 * 2.99242575222_dp  + mf2*1.29246858189_dp
           ddp3(i,j,k,2) = -0.80875619458_dp + mf1*(-2.38026356489_dp) + mf2*1.65427830900_dp
           ddp3(i,j,k,3) =  0.69028490492_dp + mf1*(-0.27012609786_dp) + mf2*(-3.43967436378_dp)
           ddp3(i,j,k,4) =  0.0_dp

        end do

     end do

  end do

end subroutine load_dd_constants


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine qq_constants
!> \brief  provides the constants of the quadrupole-term
! This subroutine provides the constants of the quadrupole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_qq_constants

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mqq_i, mqq_j, mqq_k
  real(dp)                                   :: mf1, mf2, msegij
  !-----------------------------------------------------------------------------

  do i = 1, ncomp

     mqq_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mqq_j = min( mseg(j), 2.0_dp )

        msegij = ( mqq_i * mqq_j )**0.5_dp
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        qqp2(i,j,0) =  1.237830788_dp + mf1 * 1.285410878_dp  + mf2*1.794295401_dp
        qqp2(i,j,1) =  2.435503144_dp + mf1*(-11.46561451_dp) + mf2*0.769510293_dp
        qqp2(i,j,2) =  1.633090469_dp + mf1 *22.08689285_dp   + mf2*7.264792255_dp
        qqp2(i,j,3) = -1.611815241_dp + mf1 * 7.46913832_dp   + mf2*94.48669892_dp
        qqp2(i,j,4) =  6.977118504_dp + mf1*(-17.19777208_dp) + mf2*(-77.1484579_dp)

        qqp4(i,j,0) =  0.454271755_dp + mf1*(-0.813734006_dp) + mf2*6.868267516_dp
        qqp4(i,j,1) = -4.501626435_dp + mf1 * 10.06402986_dp  + mf2*(-5.173223765_dp)
        qqp4(i,j,2) =  3.585886783_dp + mf1*(-10.87663092_dp) + mf2*(-17.2402066_dp)
        qqp4(i,j,3) =  0.0_dp
        qqp4(i,j,4) =  0.0_dp

        do k = 1, ncomp

           mqq_k = min( mseg(k), 2.0_dp )

           msegij = ( mqq_i * mqq_j * mqq_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = mf1 * ( msegij - 2.0_dp ) / msegij
           qqp3(i,j,k,0) = -0.500043713_dp + mf1 * 2.000209381_dp + mf2*3.135827145_dp
           qqp3(i,j,k,1) =  6.531869153_dp + mf1*(-6.78386584_dp) + mf2*7.247588801_dp
           qqp3(i,j,k,2) = -16.01477983_dp + mf1 * 20.38324603_dp + mf2*3.075947834_dp
           qqp3(i,j,k,3) =  14.42597018_dp + mf1*(-10.89598394_dp)
           qqp3(i,j,k,4) =  0.0_dp

        end do

     end do

  end do

end subroutine load_qq_constants

end module eos_constants

module module_eos_derivatives

  use PARAMETERS, only: dp, machine_eps, PI_6, KBOL30, nsite
  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: assoc, dipole, qudpole, dipole_quad,  &
       mseg, nhb_typ, nhb_no
  implicit none

  private

  !-----------------------------------------------------------------------------
  ! temperature and composition
  !-----------------------------------------------------------------------------
  real(dp)                                   :: t
  real(dp), allocatable,dimension(:)         :: x

  !-----------------------------------------------------------------------------
  ! various densities
  !-----------------------------------------------------------------------------
  real(dp), public                           :: rho         ! defined public

  real(dp)                                   :: z0
  real(dp)                                   :: z1
  real(dp)                                   :: z2
  real(dp), public                           :: z3          ! defined public
  real(dp)                                   :: z22
  real(dp)                                   :: z23
  real(dp)                                   :: z32
  real(dp)                                   :: z33
  real(dp)                                   :: z0t
  real(dp)                                   :: z1t
  real(dp)                                   :: z2t
  real(dp), public                           :: z3t         ! defined public

  real(dp), allocatable, dimension(:), public  :: rhoi      ! defined public

  real(dp)                                   :: ome, ome2, ome3, ome4, ome5

  real(dp), dimension(0:6)                   :: z3_m_m0
  real(dp), dimension(0:6)                   :: z3_m_m1
  real(dp), dimension(0:6)                   :: z3_m_m2

  !-----------------------------------------------------------------------------
  ! derivatives of densities, e.g. z3_rk = d(zeta3)/d(rho_k)
  !-----------------------------------------------------------------------------
  real(dp), allocatable, dimension(:)        :: z0_rk
  real(dp), allocatable, dimension(:)        :: z1_rk
  real(dp), allocatable, dimension(:)        :: z2_rk
  real(dp), allocatable, dimension(:)        :: z3_rk

  !-----------------------------------------------------------------------------
  ! abbreviation
  !-----------------------------------------------------------------------------
  real(dp), public                           :: kT

  !-----------------------------------------------------------------------------
  ! pure and mixture quantities of fluid theory
  !-----------------------------------------------------------------------------
  real(dp), allocatable, dimension(:), public :: dhs        ! defined public
  real(dp), allocatable, dimension(:,:)      :: uij
  real(dp), allocatable, dimension(:,:)      :: uij_t
  real(dp), allocatable, dimension(:,:)      :: sig_ij
  real(dp), allocatable, dimension(:,:)      :: sig3_ij
  real(dp), allocatable, dimension(:,:)      :: eps_ij
  real(dp), allocatable, dimension(:,:)      :: dij_ab

  !-----------------------------------------------------------------------------
  ! radial distribution fct at contact. Needed for chain term and association term
  !-----------------------------------------------------------------------------
  real(dp), allocatable, dimension(:,:)      :: gij
  real(dp)                                   :: gij_t1, gij_t2, gij_t3

  !-----------------------------------------------------------------------------
  ! parameter for hard-sphere model (BMCSL-extension to metastable region, Paricaud 2015)
  !-----------------------------------------------------------------------------
  logical                                    :: mod_BMCSL = .false.
  real(dp), parameter                        :: eta_jam = 0.655_dp
  real(dp), parameter                        :: jam_fact = ( 1.0_dp - eta_jam ) / eta_jam

  !-----------------------------------------------------------------------------
  ! quantities needed for dispersion term (van der Waals attraction)
  !-----------------------------------------------------------------------------
  real(dp)                                   :: m_mean
  real(dp), dimension(0:6)                   :: apar
  real(dp), dimension(0:6)                   :: bpar
  real(dp), dimension(0:6)                   :: ap_rk_aux
  real(dp), dimension(0:6)                   :: bp_rk_aux

  real(dp)                                   :: order1
  real(dp)                                   :: order2
  real(dp), allocatable, dimension(:)        :: ord1_lij_aux
  real(dp), allocatable, dimension(:)        :: ord1_lij
  real(dp), allocatable, dimension(:,:)      :: m_sig_eps

  real(dp)                                   :: I1
  real(dp)                                   :: I1_z
  real(dp)                                   :: I1_z2
  real(dp)                                   :: I2
  real(dp)                                   :: I2_z
  real(dp)                                   :: I2_z2

  real(dp)                                   :: ome_t, ome_t_2
  real(dp)                                   :: c1_con, c2_con, c3_con

  !-----------------------------------------------------------------------------
  ! quantities needed for association term
  !-----------------------------------------------------------------------------
  real(dp), allocatable, dimension(:,:)      :: mx

  real(dp), allocatable, dimension(:,:,:,:)  :: ass_d
  real(dp), allocatable, dimension(:,:,:,:)  :: ass_d_dt
  real(dp), allocatable, dimension(:,:,:,:)  :: ass_d_dt2

  real(dp), parameter                        :: hb_damp = 0.7_dp
  real(dp), parameter                        :: hb_damp_inv = 1.0_dp - hb_damp
  real(dp), parameter                        :: hb_tol = 1.E-10_dp

  !-----------------------------------------------------------------------------
  ! quantities needed for polar terms (dipole-dipole, quad.-quad., dipole-quad.)
  !-----------------------------------------------------------------------------
  real(dp), allocatable, dimension(:,:)      :: Idd2
  real(dp), allocatable, dimension(:,:)      :: Idd2_z
  real(dp), allocatable, dimension(:,:)      :: Idd2_z2
  real(dp), allocatable, dimension(:,:)      :: Iqq2
  real(dp), allocatable, dimension(:,:)      :: Iqq2_z
  real(dp), allocatable, dimension(:,:)      :: Iqq2_z2
  real(dp), allocatable, dimension(:,:)      :: Idq2
  real(dp), allocatable, dimension(:,:)      :: Idq2_z
  real(dp), allocatable, dimension(:,:)      :: Idq2_z2

  real(dp), allocatable, dimension(:,:,:)    :: Idd3
  real(dp), allocatable, dimension(:,:,:)    :: Idd3_z
  real(dp), allocatable, dimension(:,:,:)    :: Idd3_z2
  real(dp), allocatable, dimension(:,:,:)    :: Iqq3
  real(dp), allocatable, dimension(:,:,:)    :: Iqq3_z
  real(dp), allocatable, dimension(:,:,:)    :: Iqq3_z2
  real(dp), allocatable, dimension(:,:,:)    :: Idq3
  real(dp), allocatable, dimension(:,:,:)    :: Idq3_z
  real(dp), allocatable, dimension(:,:,:)    :: Idq3_z2

  real(dp), allocatable, dimension(:,:)      :: pi_dd
  real(dp), allocatable, dimension(:,:,:)    :: psi_dd
  real(dp), allocatable, dimension(:,:)      :: pi_qq
  real(dp), allocatable, dimension(:,:,:)    :: psi_qq
  real(dp), allocatable, dimension(:,:)      :: pi_dq
  real(dp), allocatable, dimension(:,:,:)    :: psi_dq

  real(dp), allocatable, dimension(:,:)      :: pi_dd_no_T
  real(dp), allocatable, dimension(:,:,:)    :: psi_dd_no_T
  real(dp), allocatable, dimension(:,:)      :: pi_qq_no_T
  real(dp), allocatable, dimension(:,:,:)    :: psi_qq_no_T
  real(dp), allocatable, dimension(:,:)      :: pi_dq_no_T
  real(dp), allocatable, dimension(:,:,:)    :: psi_dq_no_T

  real(dp), allocatable, dimension(:)        :: my_factor
  real(dp), allocatable, dimension(:)        :: Q_factor
  real(dp), allocatable, dimension(:)        :: my_fac_dq
  real(dp), allocatable, dimension(:)        :: Q_fac_dq

  real(dp), allocatable, dimension(:)        :: my_factor_no_T
  real(dp), allocatable, dimension(:)        :: Q_factor_no_T
  real(dp), allocatable, dimension(:)        :: my_fac_dq_no_T
  real(dp), allocatable, dimension(:)        :: Q_fac_dq_no_T

  real(dp)                                   :: wdd2, wdd3
  real(dp)                                   :: wdd2_z_term, wdd3_z_term
  real(dp)                                   :: wdd2_z2_term, wdd3_z2_term
  real(dp)                                   :: wqq2, wqq3
  real(dp)                                   :: wqq2_z_term, wqq3_z_term
  real(dp)                                   :: wqq2_z2_term, wqq3_z2_term
  real(dp)                                   :: wdq2, wdq3
  real(dp)                                   :: wdq2_z_term, wdq3_z_term
  real(dp)                                   :: wdq2_z2_term, wdq3_z2_term

  real(dp)                                   :: rdd2, rdd3
  real(dp)                                   :: rdd2_z_term, rdd3_z_term
  real(dp)                                   :: rdd2_z2_term, rdd3_z2_term
  real(dp)                                   :: rqq2, rqq3
  real(dp)                                   :: rqq2_z_term, rqq3_z_term
  real(dp)                                   :: rqq2_z2_term, rqq3_z2_term
  real(dp)                                   :: rdq2, rdq3
  real(dp)                                   :: rdq2_z_term, rdq3_z_term
  real(dp)                                   :: rdq2_z2_term, rdq3_z2_term

  !=============================================================================
  ! public expressions
  !=============================================================================

  ! ...(see above)...                        :: rho
  ! ...(see above)...                        :: rhoi
  ! ...(see above)...                        :: kT
  ! ...(see above)...                        :: dhs
  ! ...(see above)...                        :: z3t
  ! ...(see above)...                        :: z3

  public :: rho_independent_quantities, density_terms, F_density, f_rho, f_rho_rho,  &
       f_rho_4, F_density_rhok, F_density_rhok_T, F_density_rhok_rho,  &
       F_density_rhok_rhol, f_temp, f_temp_temp, f_temp_rho, allocate_eos_quantities,  &
       deallocate_eos_quantities, initialize_eos_quantities, Iterate_Association,  &
       cross_association_parameters, aggregate_polar_parameters


contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE density_terms
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine density_terms ( eta )

  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta
  !-----------------------------------------------------------------------------
  integer                                    :: i, j, l, m
  real(dp)                                   :: p_factor
  real(dp)                                   :: gij_t2_rho, gij_t3_rho
  real(dp)                                   :: rho_2, rho_3
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! abbreviations
  !-----------------------------------------------------------------------------
  rho = eta / z3t
  z0 = z0t * rho
  z1 = z1t * rho
  z2 = z2t * rho
  z3 = z3t * rho

  ! m_mean = z0t / PI_6
  ome = 1.0_dp - z3
  ome2 = ome * ome
  ome3 = ome2 * ome
  ome4 = ome2 * ome2
  ome5 = ome4 * ome

  z22 = z2 * z2
  z23 = z2 * z22
  z32 = z3 * z3
  z33 = z3 * z32

  rhoi( 1:ncomp ) = x(1:ncomp ) * rho

  z3_m_m0(:) = 0.0_dp
  z3_m_m1(:) = 0.0_dp
  z3_m_m2(:) = 0.0_dp
  if ( z3 > 1.E-40_dp ) then
     do m = 0, 6
        z3_m_m0(m) = z3**m
        if ( m > 0 ) z3_m_m1(m) = REAL( m, KIND=dp ) * z3**(m-1)
        if ( m > 1 ) z3_m_m2(m) = REAL( m*(m-1), KIND=dp ) * z3**(m-2)
     end do
  else
     z3_m_m0(0) = 1.0_dp
     z3_m_m1(1) = 1.0_dp
     z3_m_m2(2) = 2.0_dp
  end if

  !-----------------------------------------------------------------------------
  ! expressions needed for chain term (for i==j) and for association term (for i/=j)
  !-----------------------------------------------------------------------------
  gij_t1 = 1.0_dp/ome
  gij_t2 = 3.0_dp*z2t/ome2
  gij_t3 = 2.0_dp*z2t*z2/ome3
  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3
  do i = 1, ncomp
     gij(i,i) = gij_t1 + dij_ab(i,i)*( gij_t2_rho + dij_ab(i,i) *gij_t3_rho )
  end do

  !-----------------------------------------------------------------------------
  ! expressions needed for dispersion term
  !-----------------------------------------------------------------------------
  I1 = sum( apar(0:6) * z3_m_m0(0:6) )
  I2 = sum( bpar(0:6) * z3_m_m0(0:6) )
  I1_z = sum( apar(1:6) * z3_m_m1(1:6) )
  I2_z = sum( bpar(1:6) * z3_m_m1(1:6) )
  I1_z2 = sum( apar(2:6) * z3_m_m2(2:6) )
  I2_z2 = sum( bpar(2:6) * z3_m_m2(2:6) )

  ome_t = 1.0_dp / (ome*(2.0_dp-z3))
  ome_t_2 = ome_t * ome_t
  c1_con= 1.0_dp/ (  1.0_dp + m_mean*(8.0_dp*z3-2.0_dp*z32 )/ome4   &
       + (1.0_dp - m_mean)*(20.0_dp*z3-27.0_dp*z32 +12.0_dp*z33 -2.0_dp*z32*z32 ) * ome_t_2   )
  c2_con= - c1_con*c1_con *(m_mean*(-4.0_dp*z32 +20.0_dp*z3+8.0_dp)/ome5  &
       + (1.0_dp - m_mean) *(2.0_dp*z33 +12.0_dp*z32 -48.0_dp*z3+40.0_dp) * ome_t_2*ome_t  )
  c3_con= 2.0_dp * c2_con*c2_con/c1_con - c1_con*c1_con  &
       *( m_mean*(-12.0_dp*z32 +72.0_dp*z3+60.0_dp) / ome**6   &
       + (1.0_dp - m_mean)  &
       *(-6.0_dp*z32*z32 -48.0_dp*z33 +288.0_dp*z32- 480.0_dp*z3+264.0_dp) * ome_t_2*ome_t_2  )


  !-----------------------------------------------------------------------------
  ! expressions needed association term
  !-----------------------------------------------------------------------------
  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij(i,j) = gij_t1 + dij_ab(i,j)*( gij_t2_rho + dij_ab(i,j) *gij_t3_rho )
           gij(j,i) = gij(i,j)
        end do
     end do

  end if
  

  if ( dipole .OR. qudpole .OR. dipole_quad ) then
     rho_2 = rho * rho
     rho_3 = rho_2 * rho
  end if

  !-----------------------------------------------------------------------------
  ! expressions needed dipole-dipole term
  !-----------------------------------------------------------------------------
  rdd2 = 0.0_dp
  rdd3 = 0.0_dp
  rdd2_z_term = 0.0_dp
  rdd3_z_term = 0.0_dp
  rdd2_z2_term = 0.0_dp
  rdd3_z2_term = 0.0_dp
  if ( dipole ) then

     do i = 1, ncomp
        if ( abs( my_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j) ) > machine_eps ) then
                 Idd2(i,j)    = sum( ( ddp2(i,j,0:4) + eps_ij(i,j) * ddp4(i,j,0:4) ) * z3_m_m0(0:4) )
                 Idd2_z(i,j)  = sum( ( ddp2(i,j,1:4) + eps_ij(i,j) * ddp4(i,j,1:4) ) * z3_m_m1(1:4) )
                 Idd2_z2(i,j) = sum( ( ddp2(i,j,2:4) + eps_ij(i,j) * ddp4(i,j,2:4) ) * z3_m_m2(2:4) )
                 do l = 1, ncomp
                    Idd3(i,j,l)    = sum( ddp3(i,j,l,0:4) * z3_m_m0(0:4) )
                    Idd3_z(i,j,l)  = sum( ddp3(i,j,l,1:4) * z3_m_m1(1:4) )
                    Idd3_z2(i,j,l) = sum( ddp3(i,j,l,2:4) * z3_m_m2(2:4) )
                 end do
              end if
           end do
        end if
     end do

     do i = 1, ncomp
        if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i)*x(j)*pi_dd(i,j)
                 rdd2 = rdd2 + p_factor * Idd2(i,j)
                 rdd2_z_term = rdd2_z_term + p_factor * Idd2_z(i,j)
                 rdd2_z2_term = rdd2_z2_term + p_factor * Idd2_z2(i,j)
                 do l = 1, ncomp
                    p_factor = x(i)*x(j)*x(l)*psi_dd(i,j,l)
                    rdd3 = rdd3 + p_factor * Idd3(i,j,l)
                    rdd3_z_term = rdd3_z_term + p_factor * Idd3_z(i,j,l)
                    rdd3_z2_term = rdd3_z2_term + p_factor * Idd3_z2(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     wdd2 = rdd2 * rho_2
     wdd3 = rdd3 * rho_3
     wdd2_z_term = rdd2_z_term * rho_2
     wdd3_z_term = rdd3_z_term * rho_3
     wdd2_z2_term = rdd2_z2_term * rho_2
     wdd3_z2_term = rdd3_z2_term * rho_3

  end if

  !-----------------------------------------------------------------------------
  ! expressions needed quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  rqq2 = 0.0_dp
  rqq3 = 0.0_dp
  rqq2_z_term = 0.0_dp
  rqq3_z_term = 0.0_dp
  rqq2_z2_term = 0.0_dp
  rqq3_z2_term = 0.0_dp
  if ( qudpole ) then

     do i = 1, ncomp
        if ( abs( Q_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j) ) > machine_eps ) then
                 Iqq2(i,j)    = sum( ( qqp2(i,j,0:4) + eps_ij(i,j) * qqp4(i,j,0:4) ) * z3_m_m0(0:4) )
                 Iqq2_z(i,j)  = sum( ( qqp2(i,j,1:4) + eps_ij(i,j) * qqp4(i,j,1:4) ) * z3_m_m1(1:4) )
                 Iqq2_z2(i,j) = sum( ( qqp2(i,j,2:4) + eps_ij(i,j) * qqp4(i,j,2:4) ) * z3_m_m2(2:4) )
                 do l = 1, ncomp
                    Iqq3(i,j,l)    = sum( qqp3(i,j,l,0:4) * z3_m_m0(0:4) )
                    Iqq3_z(i,j,l)  = sum( qqp3(i,j,l,1:4) * z3_m_m1(1:4) )
                    Iqq3_z2(i,j,l) = sum( qqp3(i,j,l,2:4) * z3_m_m2(2:4) )
                 end do
              end if
           end do
        end if
     end do

     do i = 1, ncomp
        if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i)*x(j)*pi_qq(i,j)
                 rqq2 = rqq2 + p_factor * Iqq2(i,j)
                 rqq2_z_term = rqq2_z_term + p_factor * Iqq2_z(i,j)
                 rqq2_z2_term = rqq2_z2_term + p_factor * Iqq2_z2(i,j)
                 do l = 1, ncomp
                    p_factor = x(i)*x(j)*x(l)*psi_qq(i,j,l)
                    rqq3 = rqq3 + p_factor * Iqq3(i,j,l)
                    rqq3_z_term = rqq3_z_term + p_factor * Iqq3_z(i,j,l)
                    rqq3_z2_term = rqq3_z2_term + p_factor * Iqq3_z2(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     wqq2 = rqq2 * rho_2
     wqq3 = rqq3 * rho_3
     wqq2_z_term = rqq2_z_term * rho_2
     wqq3_z_term = rqq3_z_term * rho_3
     wqq2_z2_term = rqq2_z2_term * rho_2
     wqq3_z2_term = rqq3_z2_term * rho_3

  end if
  
  !-----------------------------------------------------------------------------
  ! expressions needed dipole-quadrupole term
  !-----------------------------------------------------------------------------
  rdq2 = 0.0_dp
  rdq3 = 0.0_dp
  rdq2_z_term = 0.0_dp
  rdq3_z_term = 0.0_dp
  rdq2_z2_term = 0.0_dp
  rdq3_z2_term = 0.0_dp
  if ( dipole_quad ) then

     do i = 1, ncomp
        do j = 1, ncomp
           Idq2(i,j)    = sum( ( dqp2(i,j,0:4) + eps_ij(i,j) * dqp4(i,j,0:4) ) * z3_m_m0(0:4) )
           Idq2_z(i,j)  = sum( ( dqp2(i,j,1:4) + eps_ij(i,j) * dqp4(i,j,1:4) ) * z3_m_m1(1:4) )
           Idq2_z2(i,j) = sum( ( dqp2(i,j,2:4) + eps_ij(i,j) * dqp4(i,j,2:4) ) * z3_m_m2(2:4) )
           do l = 1, ncomp
              Idq3(i,j,l)    = sum( dqp3(i,j,l,0:4) * z3_m_m0(0:4) )
              Idq3_z(i,j,l)  = sum( dqp3(i,j,l,1:4) * z3_m_m1(1:4) )
              Idq3_z2(i,j,l) = sum( dqp3(i,j,l,2:4) * z3_m_m2(2:4) )
           end do
        end do
     end do

     do i = 1, ncomp
        if ( abs( my_fac_dq(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_fac_dq(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i) * x(j) * pi_dq(i,j)
                 rdq2 = rdq2 + p_factor * Idq2(i,j)
                 rdq2_z_term = rdq2_z_term + p_factor * Idq2_z(i,j)
                 rdq2_z2_term = rdq2_z2_term + p_factor * Idq2_z2(i,j)
                 do l = 1, ncomp
                    p_factor = x(i) * x(j) * x(l) * psi_dq(i,j,l)
                    rdq3 = rdq3 + p_factor * Idq3(i,j,l)
                    rdq3_z_term = rdq3_z_term + p_factor * Idq3_z(i,j,l)
                    rdq3_z2_term = rdq3_z2_term + p_factor * Idq3_z2(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     wdq2 = rdq2 * rho_2
     wdq3 = rdq3 * rho_3
     wdq2_z_term = rdq2_z_term * rho_2
     wdq3_z_term = rdq3_z_term * rho_3
     wdq2_z2_term = rdq2_z2_term * rho_2
     wdq3_z2_term = rdq3_z2_term * rho_3

  end if

end subroutine density_terms


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine rho_independent_quantities
!
!> \brief parameters of EOS
!!
!! calculates density-independent parameters of the equation of state.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine rho_independent_quantities ( tF, xF )

  use PARAMETERS, only: PI
  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k,  &
       nhb_typ, eps_hb, kap_hb, kij, lij, kij_t, lij_correction, kij_T_dependent
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  real(dp), intent(in), dimension(ncomp)     :: xF

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, ii, jj, k, m
  real(dp)                                   :: term1, term2, one_third
  real(dp)                                   :: kij_par
  !-----------------------------------------------------------------------------

  t = tF
  x( 1:ncomp ) = xF( 1:ncomp )

  kT = KBOL30 * t                        ! in units [Pa*Angstrom**3]

  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------
  do  i = 1, ncomp
     dhs(i) = sigma(i) * ( 1.0_dp - 0.12_dp * EXP( -3.0_dp * epsilon_k(i) / t ) )
  end do

  !-----------------------------------------------------------------------------
  ! combination rules
  !-----------------------------------------------------------------------------
  do i = 1, ncomp
     do j = 1, ncomp
        sig_ij(i,j) = 0.5_dp * ( sigma(i) + sigma(j) )
        sig3_ij(i,j) = ( sig_ij(i,j) )**3
        kij_par = kij(i,j)
        if ( kij_T_dependent ) then
           kij_par = kij_par + kij_t(i,j) * ( t - 300.0_dp ) / 100.0_dp
        end if
        uij(i,j) = ( 1.0_dp - kij_par ) * ( epsilon_k(i) * epsilon_k(j) )**0.5
        uij_t(i,j) = uij(i,j) / t
        eps_ij(i,j) = ( epsilon_k(i) * epsilon_k(j) )**0.5 / t
     end do
  end do


  !-----------------------------------------------------------------------------
  ! abbreviations
  !-----------------------------------------------------------------------------
  z0t = PI_6 * SUM( x(1:ncomp) * mseg(1:ncomp) )
  z1t = PI_6 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp) )
  z2t = PI_6 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**2 )
  z3t = PI_6 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  do  i = 1, ncomp                       ! derivatives of z0,..z3 to density rho_k
     z0_rk(i) = PI_6 * mseg(i)
     z1_rk(i) = z0_rk(i) * dhs(i)
     z2_rk(i) = z1_rk(i) * dhs(i)
     z3_rk(i) = z2_rk(i) * dhs(i)
  end do
  
  m_mean = z0t / PI_6                    ! mean segment number

  do i = 1, ncomp
     do j = 1, ncomp
        dij_ab(i,j) = dhs(i)*dhs(j) / ( dhs(i) + dhs(j) )
     end do
  end do

  !-----------------------------------------------------------------------------
  ! dispersion term parameters for chain molecules
  !-----------------------------------------------------------------------------
  term1 = 1.0_dp - 1.0_dp / m_mean
  term2 = term1 * ( 1.0_dp - 2.0_dp / m_mean )
  do m = 0, 6
     apar(m) = ap(m,1) + term1 * ap(m,2) + term2 * ap(m,3)
     bpar(m) = bp(m,1) + term1 * bp(m,2) + term2 * bp(m,3)
  end do

  !------ auxiliary quantity needed for d(apar) / d(rho_k)
  ap_rk_aux(0:6) = ( ap(0:6,2) + (3.0_dp -4.0_dp/m_mean) *ap(0:6,3) ) / m_mean**2
  bp_rk_aux(0:6) = ( bp(0:6,2) + (3.0_dp -4.0_dp/m_mean) *bp(0:6,3) ) / m_mean**2  

  !-----------------------------------------------------------------------------
  ! van der Waals mixing rules for perturbation terms
  !-----------------------------------------------------------------------------
  order1 = 0.0_dp
  order2 = 0.0_dp
  do i = 1, ncomp
     do j = 1, ncomp
        term1 = x(i)*x(j)* mseg(i)*mseg(j)*sig3_ij(i,j) * uij_t(i,j)
        order1 = order1 + term1
        order2 = order2 + term1 * uij_t(i,j)
     end do
  end do
  order1 = - 2.0_dp * PI * order1
  order2 = - PI * order2

  if ( lij_correction ) then  ! lij is non-symmetric: caution with sequence of indices

     one_third = 1.0_dp / 3.0_dp
     do i = 1, ncomp
        ord1_lij_aux(i) = 0.0_dp   
        do j = 1, ncomp
           if ( abs( lij(i,j) ) > machine_eps ) then
              m_sig_eps(i,j) = mseg(j)*sig_ij(i,j) *sign( abs( uij_t(i,j)*lij(i,j) )**one_third, lij(i,j) )
              ord1_lij_aux(i) = ord1_lij_aux(i) + x(j)*m_sig_eps(i,j)
           end if
        end do
        order1 = order1 - 2.0_dp * PI * x(i)*mseg(i) * ord1_lij_aux(i)**3
     end do

     do k = 1, ncomp
        ord1_lij(k) = mseg(k) * ord1_lij_aux(k)**3
        do i = 1, ncomp
           if ( abs( ord1_lij_aux(i) ) > machine_eps ) then
              ord1_lij(k) = ord1_lij(k) + x(i)*mseg(i) * ord1_lij_aux(i)**2  &
                                          * ( 3.0_dp*m_sig_eps(i,k) - 2.0_dp*ord1_lij_aux(i) )
           end if
        end do
        ord1_lij(k) = - 2.0_dp *PI * ord1_lij(k)
     end do
     
  end if  ! lij correction

  !-----------------------------------------------------------------------------
  ! constants and parameters of polar terms
  !-----------------------------------------------------------------------------

  ! ------ dipole-dipole term --------------------------------------------------
  if ( dipole ) then

     my_factor(:) = my_factor_no_T(:) / t

     pi_dd(:,:) = pi_dd_no_T(:,:) / t**2
     psi_dd(:,:,:) = psi_dd_no_T(:,:,:) / t**3

  end if  ! dipole
  
  ! ------ quadrupole-quadrupole term ------------------------------------------
  if ( qudpole ) then

     Q_factor(:) = Q_factor_no_T(:) / t

     pi_qq(:,:) = pi_qq_no_T(:,:) / t**2
     psi_qq(:,:,:) = psi_qq_no_T(:,:,:) / t**3

  end if  !qudpole
  
  ! ------ dipole-quadrupole term ----------------------------------------------
  if ( dipole_quad ) then

     my_fac_dq(:) = my_fac_dq_no_T(:) / t

     Q_fac_dq(:) = Q_fac_dq_no_T(:) / t

     pi_dq(:,:) = pi_dq_no_T(:,:) / t**2
     psi_dq(:,:,:) = psi_dq_no_T(:,:,:) / t**3

  end if



  !-----------------------------------------------------------------------------
  ! TPT-1-association parameters
  !-----------------------------------------------------------------------------
  if ( assoc ) then

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 term1 = kap_hb(i,j) *sig3_ij(i,j)
                 ass_d(i,j,ii,jj) = term1 * ( EXP(eps_hb(i,j,ii,jj)/t) - 1.0_dp )
                 term2 = term1 *eps_hb(i,j,ii,jj)/t/t * EXP( eps_hb(i,j,ii,jj)/t )
                 ass_d_dt(i,j,ii,jj) = - term2
                 ass_d_dt2(i,j,ii,jj) = - term2 * (-2.0_dp/t - eps_hb(i,j,ii,jj)/t/t)
              end do
           end do
        end do
     end do

  end if

end subroutine rho_independent_quantities

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine cross_association_parameters
!
!> \brief define binary association parameters
!!
!! these parameters are independent of T, rho, x
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine cross_association_parameters

  use pcsaft_pure_and_binary_parameters, only: sigma, nhb_typ, eps_hb, kap_hb,  &
       assoc_scheme, kij_assoc
  implicit none

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, ii, jj
  real(dp)                                   :: sig3
  !-----------------------------------------------------------------------------

  if ( assoc ) then

     do i = 1, ncomp

        do j = 1, ncomp

           if ( i /= j .AND. (nhb_typ(i) /= 0 .AND. nhb_typ(j) /= 0) ) then

              sig3 = ( 0.5_dp * ( sigma(i) + sigma(j) ) )**3

              kap_hb(i,j) = ( kap_hb(i,i) * kap_hb(j,j) )**0.5_dp  &
                   *( (sigma(i)*sigma(j))**3 )**0.5_dp / sig3

              !-----------------------------------------------------------------
              ! association scheme: both substances either 2B, 3B or 4C
              !-----------------------------------------------------------------
              if ( ( assoc_scheme(i) == '2B' .OR. assoc_scheme(i) == '3B'  &
                                             .OR. assoc_scheme(i) == '4C' ) .AND. &
                   ( assoc_scheme(j) == '2B' .OR. assoc_scheme(j) == '3B'  &
                                             .OR. assoc_scheme(j) == '4C' ) ) then

                 !JRE forall( ii = 1:nhb_typ(i), jj = 1:nhb_typ(j), ii /= jj )
                do jj=1,nhb_typ(j)
                    do ii=jj+1,nhb_typ(i)
                         eps_hb(i,j,ii,jj) = 0.5_dp * ( eps_hb(i,i,ii,jj) + eps_hb(j,j,jj,ii) )*( 1.0_dp - kij_assoc(i,j) )
                         eps_hb(i,j,jj,ii) =eps_hb(i,j,ii,jj) !JRE
                    enddo
                enddo
                 !JRE end forall
              end if

              !-----------------------------------------------------------------
              ! association scheme: one species 1A, the other either 2B, 3B or 4C
              !-----------------------------------------------------------------
              if ( assoc_scheme(i) == '1A' .AND. &
                 ( assoc_scheme(j) == '2B' .OR. assoc_scheme(j) == '3B'  &
                                           .OR. assoc_scheme(j) == '4C' ) ) then

                 do jj = 1, nhb_typ(j)

                    eps_hb(i,j,1,jj) = 0.5_dp * ( eps_hb(i,i,1,1) + eps_hb(j,j,jj,(3-jj)) )  &
                         * ( 1.0_dp - kij_assoc(i,j) )
                    eps_hb(j,i,jj,1) = eps_hb(i,j,1,jj)

                 end do

              end if

              !-----------------------------------------------------------------
              ! association scheme: two species 1A
              !-----------------------------------------------------------------
              if ( assoc_scheme(i) == '1A' .AND. assoc_scheme(j) == '1A' ) then

                 eps_hb(i,j,1,1) = 0.5_dp * ( eps_hb(i,i,1,1) + eps_hb(j,j,1,1) )  &
                                   * ( 1.0_dp - kij_assoc(i,j) )

              end if

           end if

        end do

     end do

  end if

end subroutine cross_association_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine aggregate_polar_parameters
!> \brief combine some parameters of polar terms (independent of T, rho, x)
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine aggregate_polar_parameters

  use PARAMETERS, only: PI, KBOL
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k, dipole_moment,  &
       quadru_moment
  implicit none

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, l
  real(dp)                                   :: factor2, factor3
  real(dp)                                   :: my_square, Q_square
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! combination rules, these quantities are actually determined in
  ! SR rho_independent_quantities, but are already needed here
  !-----------------------------------------------------------------------------
  do i = 1, ncomp

     do j = 1, ncomp

        sig_ij(i,j) = 0.5_dp * ( sigma(i) + sigma(j) )
        sig3_ij(i,j) = ( sig_ij(i,j) )**3

     end do

  end do

  !-----------------------------------------------------------------------------
  ! dipole-dipole term
  !-----------------------------------------------------------------------------
  if ( dipole ) then

     do i = 1, ncomp
        my_square = dipole_moment(i)**2 *1.E-19_dp / ( epsilon_k(i)*KBOL*mseg(i)*sig3_ij(i,i) )
        my_factor_no_T(i) = my_square * epsilon_k(i) * sig3_ij(i,i)
     end do

     factor2 = - PI
     factor3 = - 4.0_dp/3.0_dp*PI*PI

     do i = 1, ncomp
        do j = 1, ncomp
           pi_dd_no_T(i,j) = factor2 * my_factor_no_T(i) * my_factor_no_T(j) / sig3_ij(i,j)
           do l = 1, ncomp
              psi_dd_no_T(i,j,l) = factor3 / sig_ij(i,j) / sig_ij(i,l) / sig_ij(j,l)  &
                   * my_factor_no_T(i) * my_factor_no_T(j) * my_factor_no_T(l)

           end do
        end do
     end do

  end if

  !-----------------------------------------------------------------------------
  ! quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  if ( qudpole ) then

     do i = 1, ncomp
        Q_square = quadru_moment(i)**2 *1.E-19_dp / ( epsilon_k(i)*kbol*mseg(i)*sig_ij(i,i)**5 )
        Q_factor_no_T(i) = Q_square * epsilon_k(i) * sig_ij(i,i)**5
     end do

     factor2 = -9.0_dp/16.0_dp*PI
     factor3 =  9.0_dp/16.0_dp*PI**2

     do i = 1, ncomp
        do j = 1, ncomp
           pi_qq_no_T(i,j) = factor2 * Q_factor_no_T(i) * Q_factor_no_T(j) / sig_ij(i,j)**7
           do l = 1, ncomp
              psi_qq_no_T(i,j,l) = factor3 * Q_factor_no_T(i) * Q_factor_no_T(j) * Q_factor_no_T(l) &
                   /sig3_ij(i,j)/sig3_ij(i,l)/sig3_ij(j,l)
           end do
        end do
     end do

  end if

  !-----------------------------------------------------------------------------
  ! dipole-quadrupole term
  !-----------------------------------------------------------------------------
  if ( dipole_quad ) then

     DO i=1,ncomp
        my_square = dipole_moment(i)**2 *1.E-19_dp / ( epsilon_k(i)*KBOL*mseg(i)*sig3_ij(i,i) )
        my_fac_dq_no_T(i) = my_square * epsilon_k(i) * sig_ij(i,i)**4     ! defined differently than for dipole-dipole term
        Q_square = (quadru_moment(i))**2 *1.E-19_dp / (epsilon_k(i)*kbol*mseg(i)*sig_ij(i,i)**5 )
        Q_fac_dq_no_T(i) = Q_square * epsilon_k(i) * sig_ij(i,i)**4       ! defined differently than for quad.-quad. term
     END DO

     factor2 = -9.0_dp/4.0_dp*PI
     factor3 =  PI**2

     ! caution: the matrix pi_dq and psi_dq are non-symmetric
     do i = 1, ncomp
        do j = 1, ncomp
           pi_dq_no_T(i,j) = factor2 * my_fac_dq_no_T(i) * Q_fac_dq_no_T(j) / sig_ij(i,j)**5
           do l = 1, ncomp
              psi_dq_no_T(i,j,l) = factor3 / (sig_ij(i,j)*sig_ij(i,l)*sig_ij(j,l))**2  &
                   * ( my_fac_dq_no_T(i)*Q_fac_dq_no_T(j)*my_fac_dq_no_T(l)  &
                   + my_fac_dq_no_T(i)*Q_fac_dq_no_T(j)*Q_fac_dq_no_T(l)*1.1937350_dp )

           end do
        end do
     end do

  end if

end subroutine aggregate_polar_parameters



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine f_rho
!
! calculates the derivative of a = A/(NkT) to density rho.  The actual output
! argument is the compressibility factor z_res = rho* d(a)/d(rho)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_rho ( eta, z_res )

  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: eta
  real(dp), intent(out)                           :: z_res
  !-----------------------------------------------------------------------------
  integer                                         :: i, j, ii, jj
  integer                                         :: ass_cnt, max_eval
  real(dp)                                        :: factor_i, factor_ij
  real(dp)                                        :: exp_term
  real(dp)                                        :: one_minus_phi
  real(dp)                                        :: phi_term

  real(dp)                                        :: zhs
  real(dp)                                        :: zhc
  real(dp), dimension(ncomp,ncomp)                :: gij_z
  real(dp)                                        :: gij_z_t1, gij_z_t2, gij_z_t3

  real(dp)                                        :: zdsp
  real(dp)                                        :: edI1dz, edI2dz

  real(dp)                                        :: zhb
  real(dp), dimension(ncomp,ncomp,nsite,nsite)    :: delta, dq_dz
  real(dp), dimension(ncomp,nsite)                :: mx_itr
  real(dp)                                        :: err_sum, sum0, tol

  real(dp)                                        :: rho_2, rho_fac, f3_to_2, diff_factor
  real(dp)                                        :: zdd
  real(dp)                                        :: fdd2, fdd2_z, fdd3, fdd3_z, fdd_z

  real(dp)                                        :: zqq
  real(dp)                                        :: fqq2, fqq2_z, fqq3, fqq3_z, fqq_z

  real(dp)                                        :: zdq
  real(dp)                                        :: fdq2, fdq2_z, fdq3, fdq3_z, fdq_z
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! density-expressions used in many partial derivatives of Helmholtz energy
  !-----------------------------------------------------------------------------
  call density_terms ( eta )

  !-----------------------------------------------------------------------------
  ! p : hard sphere contribution
  !-----------------------------------------------------------------------------
  ! fhs = ( 3.0_dp*z1*z2/ome + z2**3 /z3/ome2 + (z2**3 /z3/z3-z0)*LOG(ome) ) / rhoPI_6
  if ( .NOT. mod_BMCSL ) then
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t*(3.0_dp-z3)/ome3 ) / PI_6
  else
     exp_term = z3 * jam_fact * exp( 50.0_dp * ( z3 - eta_jam ) )
     one_minus_phi = 1.0_dp - z3 - exp_term
     phi_term = 50.0_dp * z3 * exp_term
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z2*z2t*z2t/z3t/ome/one_minus_phi  &
          * ( (1.0_dp+phi_term)/one_minus_phi - one_minus_phi + z3/ome ) ) / PI_6
  end if

  !-----------------------------------------------------------------------------
  ! p : chain term
  !-----------------------------------------------------------------------------
  gij_z_t1 = 1.0_dp / ome2
  gij_z_t2 = gij_t2 * ( 1.0_dp/z3t + 2.0_dp*rho/ome )
  gij_z_t3 = gij_t3 * (2.0_dp+z3) / z3t / ome
  do i = 1, ncomp
     gij_z(i,i) = gij_z_t1 + dij_ab(i,i)*( gij_z_t2 + dij_ab(i,i) *gij_z_t3 )
  end do

  zhc = 0.0_dp
  do i= 1, ncomp
     zhc = zhc + x(i) * (1.0_dp-mseg(i)) * z3 / gij(i,i)* gij_z(i,i)
  end do

  !-----------------------------------------------------------------------------
  ! p : PC-SAFT dispersion contribution
  !     note: edI1dz is equal to d(z3*I1)/d(z3), analogous for edI2dz
  !-----------------------------------------------------------------------------
  edI1dz = z3 * I1_z + I1
  edI2dz = z3 * I2_z + I2

  zdsp = rho * edI1dz * order1  &
       + order2 * m_mean * rho * (c2_con*I2*z3 + c1_con*edI2dz)


  !-----------------------------------------------------------------------------
  ! p: TPT-1-association accord. to Chapman et al.
  !-----------------------------------------------------------------------------
  zhb = 0.0_dp

  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij_z(i,j) = gij_z_t1 + dij_ab(i,j)*( gij_z_t2 + dij_ab(i,j) *gij_z_t3 )
           gij_z(j,i) = gij_z(i,j)
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 delta(i,j,ii,jj) = gij(i,j)    * ass_d(i,j,ii,jj)
                 dq_dz(i,j,ii,jj) = gij_z(i,j) * ass_d(i,j,ii,jj)
              end do
           end do
        end do
     end do

     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-4_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-5_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx(i,ii) = 1.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum0 = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    sum0 = sum0 + x(j)*nhb_no(j,jj)* mx(j,jj)* delta(i,j,ii,jj)
                 end do
              end do
              mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + sum0 * rho )
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS( mx_itr(i,ii) - mx(i,ii) )
              mx(i,ii) = mx_itr(i,ii) * hb_damp + mx(i,ii) * hb_damp_inv
           end do
        end do

        if ( ass_cnt == max_eval .AND. err_sum > SQRT(tol) ) then
           write (*,'(a,2G15.7)') 'F_rho: Max_eval violated (mx) Err_Sum= ', err_sum, tol
           exit
        end if

     end do

     ! --- calculate the hydrogen-bonding contribution -------------------------
     !     zhb = rho* d(F_hb/NkT) / d(rho) = rho* d(F_hb/NkT) / d(z3) * z3t
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_i = 0.5_dp * rhoi(i) * nhb_no(i,ii) * z3t
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 factor_ij = factor_i * rhoi(j) * nhb_no(j,jj)
                 zhb = zhb - factor_ij * mx(i,ii)*mx(j,jj) * dq_dz(i,j,ii,jj)
              end do
           end do
        end do
        zhb = zhb - 0.5_dp * x(i) * sum( nhb_no( i, 1:nhb_typ(i) )* (1.0_dp -mx( i, 1:nhb_typ(i) )) )
     end do

  end if


  !-----------------------------------------------------------------------------
  ! p: polar terms
  !-----------------------------------------------------------------------------
  zdd = 0.0_dp
  zqq = 0.0_dp
  zdq = 0.0_dp

  if ( dipole .OR. qudpole .OR. dipole_quad ) then
     rho_2 = rho * rho
     rho_fac = 2.0_dp * rho / z3t
  end if

  ! ------ dipole-dipole term --------------------------------------------------
  if ( dipole ) then

     if ( abs( rdd2 ) > 1.E-50_dp ) then

        f3_to_2 = rdd3 / rdd2 * rho
        diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

        fdd2 = rdd2 * rho
        fdd3 = rdd3 * rho_2
        fdd2_z = rdd2/z3t + rdd2_z_term * rho
        fdd3_z = rho_fac * rdd3 + rdd3_z_term * rho_2

        !fdd_z = fdd2* (fdd2*fdd2_z - 2.0_dp*fdd3*fdd2_z + fdd2*fdd3_z) / (fdd2-fdd3)**2
        fdd_z = ( fdd2_z + ( fdd3_z - f3_to_2*fdd2_z ) * diff_factor ) * diff_factor

        zdd   = fdd_z * z3

     end if

  end if
  
  ! ------ quadrupole-quadrupole term ------------------------------------------
  if ( qudpole ) then

     if ( abs( rqq2 ) > 1.E-50_dp ) then

        f3_to_2 = rqq3 / rqq2 * rho
        diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

        fqq2 = rqq2 * rho
        fqq3 = rqq3 * rho_2
        fqq2_z = rqq2/z3t + rqq2_z_term * rho
        fqq3_z = rho_fac * rqq3 + rqq3_z_term * rho_2

        fqq_z = ( fqq2_z + ( fqq3_z - f3_to_2*fqq2_z ) * diff_factor ) * diff_factor

        zqq   = fqq_z * z3

     end if

  end if
  
  ! ------ dipole-quadrupole term ----------------------------------------------
  if ( dipole_quad ) then

     if ( abs( rdq2 ) > 1.E-50_dp ) then

        f3_to_2 = rdq3 / rdq2 * rho
        diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

        fdq2 = rdq2 * rho
        fdq3 = rdq3 * rho_2
        fdq2_z = rdq2/z3t + rdq2_z_term * rho
        fdq3_z = rho_fac * rdq3 + rdq3_z_term * rho_2

        fdq_z = ( fdq2_z + ( fdq3_z - f3_to_2*fdq2_z ) * diff_factor ) * diff_factor

        zdq   = fdq_z * z3

     end if

  end if


  !-----------------------------------------------------------------------------
  ! compressibility factor z and total p
  ! as well as derivatives d(z)/d(eta) and d(p)/d(eta) with unit [Pa]
  !-----------------------------------------------------------------------------
  z_res = zhs + zhc + zdsp + zhb + zdd + zqq + zdq

  ! f_r = z_res / rho

  ! pcalc = ( z_res + 1.0_dp ) * rho * kT

end subroutine f_rho


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine f_rho_rho
!
! calculates the first and second derivative of a = A/(NkT) to density rho.
! The actual outputargument is the residual compressibility factor
! z_res = rho* d(a)/d(rho).
! The second derivative of a = A/(NkT), is written as the first derivative of  z
! to density (to packing fraction eta), with
! ztot_z = ( rho*d2(a)/d(rho)**2 + d(a)/d(rho) ) * rho/eta
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_rho_rho ( eta, z_res, ztot_z )

  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: eta
  real(dp), intent(out)                           :: z_res
  real(dp), intent(out)                           :: ztot_z
  !-----------------------------------------------------------------------------

  integer                                         :: i, j, ii, jj
  integer                                         :: ass_cnt,max_eval
  real(dp)                                        :: abbrev
  real(dp)                                        :: factor_i, factor_ij
  real(dp)                                        :: exp_term
  real(dp)                                        :: one_minus_phi
  real(dp)                                        :: phi_term, phi_r, z3_50

  real(dp)                                        :: zhs, zhs_z
  real(dp)                                        :: zhc, zhc_z
  real(dp), dimension(ncomp,ncomp)                :: gij_z, gij_z2
  real(dp)                                        :: gij_z_t1, gij_z_t2, gij_z_t3
  real(dp)                                        :: gij_z2_t1, gij_z2_t2, gij_z2_t3

  real(dp)                                        :: zdsp, zdsp_z
  real(dp)                                        :: edI1dz, edI2dz, edI1d2, edI2d2

  real(dp)                                        :: zhb, zhb_z
  real(dp), dimension(ncomp,nsite)                :: mx_itr, mx_z, mx_z_itr
  real(dp)                                        :: err_sum, sum0, sum1, tol

  real(dp)                                        :: zdd, zdd_z
  real(dp)                                        :: zqq, zqq_z
  real(dp)                                        :: zdq, zdq_z
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! density-expressions used in many partial derivatives of Helmholtz energy
  !-----------------------------------------------------------------------------
  call density_terms ( eta )

  !-----------------------------------------------------------------------------
  ! p : hard sphere contribution
  !-----------------------------------------------------------------------------
  ! fhs= ( 3.0_dp*z1*z2/ome + z2**3 /z3/ome2 + (z2**3 /z3/z3-z0)*LOG(ome) ) / rhoPI_6
  if ( .NOT. mod_BMCSL ) then
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t*(3.0_dp-z3)/ome3 ) / PI_6
     zhs_z = ( z0t + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome + 6.0_dp*z2*z2t**2 /z3t/ome2 ) /ome2 / PI_6
  else
     z3_50 = 50.0_dp * z3
     exp_term = z3 * jam_fact * exp( z3_50 - 50.0_dp * eta_jam  )
     one_minus_phi = 1.0_dp - z3 - exp_term
     phi_term = z3 + ( 1.0_dp + z3_50 ) * exp_term
     phi_r = phi_term + ( 2.0_dp + z3_50 ) * z3_50 * exp_term
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z2*z2t*z2t/z3t/ome/one_minus_phi  &
          * ( phi_term/one_minus_phi - one_minus_phi + 1.0_dp/ome ) ) / PI_6
     !zhs = rho*( -(z2t**3/z3t**2-z0t)*z3t/ome + 3.0_dp*z1t*z2t/ome2 + z2t**3/z3t/ome/one_minus_phi  &
     !     * ( 1.0_dp + z3/ome + phi_term/one_minus_phi ) ) / PI_6
     zhs_z = ( z0t/ome2 + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome3  &
          + z2t**3/z3t**2/ome/one_minus_phi * ( (1.0_dp/ome + phi_term/one_minus_phi)**2  &
          - one_minus_phi/ome + z3/ome2 + (phi_term/one_minus_phi)**2 + phi_r/one_minus_phi ) ) / PI_6
  end if

  !-----------------------------------------------------------------------------
  ! p : chain term
  !-----------------------------------------------------------------------------
  gij_z_t1 = 1.0_dp / ome2
  gij_z_t2 = gij_t2 * ( 1.0_dp/z3t + 2.0_dp*rho/ome )
  gij_z_t3 = gij_t3 * (2.0_dp+z3) / z3t / ome
  gij_z2_t1 = 2.0_dp / ome3
  gij_z2_t2 = 6.0_dp*z2t/z3t/ome4 *(2.0_dp+z3)
  gij_z2_t3 = 4.0_dp*(z2t/z3t)**2 /ome5  *(1.0_dp+4.0_dp*z3+z32)
  do i = 1, ncomp
     gij_z(i,i) = gij_z_t1 + dij_ab(i,i)*( gij_z_t2 + dij_ab(i,i) *gij_z_t3 )
     gij_z2(i,i) = gij_z2_t1 + dij_ab(i,i)*( gij_z2_t2 + dij_ab(i,i) *gij_z2_t3 )
  end do

  zhc = 0.0_dp
  zhc_z = 0.0_dp
  do i= 1, ncomp
     abbrev = x(i) * (1.0_dp-mseg(i)) / gij(i,i)
     zhc = zhc + abbrev * z3 * gij_z(i,i)
     zhc_z = zhc_z + abbrev * ( gij_z(i,i) * ( 1.0_dp - z3 * gij_z(i,i) /gij(i,i) )   &
          + z3 * gij_z2(i,i) )
  end do

  !-----------------------------------------------------------------------------
  ! p : PC-SAFT dispersion contribution
  !     note: edI1dz is equal to d(z3*I1)/d(z3), analogous for edI2dz
  !-----------------------------------------------------------------------------
  edI1dz = z3 * I1_z + I1
  edI2dz = z3 * I2_z + I2

  edI1d2 = z3 * I1_z2 + 2.0_dp * I1_z
  edI2d2 = z3 * I2_z2 + 2.0_dp * I2_z

  abbrev = edI1dz*order1 + order2 * m_mean * (c2_con*I2*z3 + c1_con*edI2dz)
  zdsp = rho * abbrev
  zdsp_z = abbrev/z3t + rho* ( edI1d2*order1  &
       + order2*m_mean *(c3_con*I2*z3 + 2.0_dp*c2_con*edI2dz + c1_con*edI2d2) )
     
  !-----------------------------------------------------------------------------
  ! p: TPT-1-association accord. to Chapman et al.
  !-----------------------------------------------------------------------------
  zhb   = 0.0_dp
  zhb_z = 0.0_dp

  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij_z(i,j) = gij_z_t1 + dij_ab(i,j)*( gij_z_t2 + dij_ab(i,j) *gij_z_t3 )
           gij_z2(i,j) = gij_z2_t1 + dij_ab(i,j)*( gij_z2_t2 + dij_ab(i,j) *gij_z2_t3 )
           gij_z(j,i) = gij_z(i,j)
           gij_z2(j,i) = gij_z2(i,j)
        end do
     end do

     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-4_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-5_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx(i,ii) = 1.0_dp
           mx_z(i,ii) = 0.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum0 = 0.0_dp
              sum1 = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_ij = x(j) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
                    sum0 = sum0 + factor_ij * mx(j,jj) * gij(i,j)
                    sum1 = sum1 + factor_ij * ( mx(j,jj) * gij_z(i,j) + mx_z(j,jj)* gij(i,j) )
                 end do
              end do
              mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + sum0*rho )
              mx_z_itr(i,ii) = -( mx_itr(i,ii)*mx_itr(i,ii) )* ( sum0/z3t + sum1*rho )

           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS(mx_itr(i,ii) - mx(i,ii))  &
                   + ABS(mx_z_itr(i,ii) - mx_z(i,ii))
              mx(i,ii)   = mx_itr(i,ii)*hb_damp   + mx(i,ii) * hb_damp_inv
              mx_z(i,ii) = mx_z_itr(i,ii)*hb_damp + mx_z(i,ii) * hb_damp_inv
           end do
        end do

        if ( ass_cnt == max_eval .AND. err_sum > SQRT(tol) ) then
           write (*,'(a,2G15.7)') 'F_rho_rho: Max_eval violated (mx) Err_Sum= ', err_sum, tol
           exit
        end if

     end do

     ! --- calculate the hydrogen-bonding contribution -------------------------
     !     zhb = rho* d(F_hb/NkT) / d(rho) = rho* d(F_hb/NkT) / d(z3) * z3t
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_i = 0.5_dp * rhoi(i) * nhb_no(i,ii) * z3t
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 factor_ij = factor_i * x(j) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
                 zhb = zhb - factor_ij * mx(i,ii)*mx(j,jj) * rho * gij_z(i,j)
                 zhb_z = zhb_z - factor_ij  &
                      * ( mx(i,ii)*mx(j,jj)* ( 2.0_dp*gij_z(i,j)/z3t + rho * gij_z2(i,j) )  &
                         + ( mx_z(i,ii)*mx(j,jj)+mx(i,ii)*mx_z(j,jj) ) * rho * gij_z(i,j) )
              end do
           end do
        end do
        zhb = zhb - 0.5_dp * x(i) * sum( nhb_no( i, 1:nhb_typ(i) )* (1.0_dp -mx( i, 1:nhb_typ(i) )) )
        zhb_z = zhb_z + 0.5_dp * x(i) * sum( nhb_no( i, 1:nhb_typ(i) )* mx_z( i, 1:nhb_typ(i) ) )
     end do

  end if

  !-----------------------------------------------------------------------------
  ! p: polar terms
  !-----------------------------------------------------------------------------
  zdd = 0.0_dp
  zdd_z = 0.0_dp
  zqq = 0.0_dp
  zqq_z = 0.0_dp
  zdq = 0.0_dp
  zdq_z = 0.0_dp
  if ( dipole ) call f_dd_rho_rho ( zdd, zdd_z )
  if ( qudpole ) call f_qq_rho_rho ( zqq, zqq_z )
  if ( dipole_quad ) call f_dq_rho_rho ( zdq, zdq_z )

  !-----------------------------------------------------------------------------
  ! compressibility factor z and total p
  ! as well as derivatives d(z)/d(z3) and d(p)/d(z3) with unit [Pa]
  !-----------------------------------------------------------------------------
  z_res = zhs + zhc + zdsp + zhb + zdd + zqq + zdq	  !JRE: pause here for sample calcs

  ztot_z = zhs_z + zhc_z + zdsp_z + zhb_z + zdd_z + zqq_z + zdq_z

  ! f_r = z_res / rho
  ! f_r2 = ( ztot_z*z3t - f_r ) / rho

  ! ztot = 1.0_dp + z_res
  ! pges = ztot * rho * kT
  ! pgesdz = ( ztot_z*rho + ztot/z3t ) * kT

end subroutine f_rho_rho

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine f_rho_4
!
!> \brief first, second, third and fourth derivative of Helmholtz energy to density.
!!
! Calculates the first and second derivative of a = A/(NkT) to density rho.
! The actual outputargument is the residual compressibility factor
! z_res = rho* d(a)/d(rho).
! The second derivative of a = A/(NkT), is written as the first derivative of  z
! to density (to packing fraction eta), with
! ztot_z = ( rho*d2(a)/d(rho)**2 + d(a)/d(rho) ) * rho/eta
!
! The third and fourth derivative are given as derivative of pressure (in unit [Pa])
! to packing fraction eta ( = z3).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_rho_4 ( eta, z_res, ztot_z, pcalc_z2, pcalc_z3 )

  !-----------------------------------------------------------------------------
  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta
  real(dp), intent(out)                      :: z_res
  real(dp), intent(out)                      :: ztot_z
  real(dp), intent(out)                      :: pcalc_z2
  real(dp), intent(out)                      :: pcalc_z3

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, ii, jj, m
  integer                                    :: ass_cnt, max_eval

  real(dp)                                   :: zhs, zhs_z, zhs_z2, zhs_z3
  real(dp)                                   :: abbrev
  real(dp)                                   :: z2to3
  real(dp)                                   :: exp_term
  real(dp)                                   :: one_minus_phi
  real(dp)                                   :: phi_term, z3_50, phi, zhs1_z2, zhs2_z2
  real(dp)                                   :: chi1, chi2, chi2_aux, chi1_z, phi_z, chi2_z

  real(dp)                                   :: zhc, zhc_z, zhc_z2, zhc_z3
  real(dp), dimension(ncomp,ncomp)           :: gij_z, gij_z2, gij_z3, gij_z4
  real(dp)                                   :: gij_z_t1, gij_z_t2, gij_z_t3
  real(dp)                                   :: gij_z2_t1, gij_z2_t2, gij_z2_t3
  real(dp)                                   :: gij_z3_t1, gij_z3_t2, gij_z3_t3
  real(dp)                                   :: gij_z4_t1, gij_z4_t2, gij_z4_t3

  real(dp)                                   :: zdsp, zdsp_z, zdsp_z2, zdsp_z3
  real(dp)                                   :: c4_con, c5_con
  real(dp)                                   :: edI1dz, edI2dz, edI1d2, edI2d2
  real(dp)                                   :: edI1d3, edI2d3, edI1d4, edI2d4

  real(dp)                                   :: zhb, zhb_z, zhb_z2, zhb_z3
  real(dp), dimension(ncomp,ncomp,nsite,nsite)  :: delta, dq_dz, dq_d2, dq_d3, dq_d4
  real(dp), dimension(ncomp,nsite)           :: mx_itr, mx_z, mx_z_itr, mx_z2, mx_z2_itr
  real(dp), dimension(ncomp,nsite)           :: mx_z3, mx_z3_itr, mx_z4, mx_z4_itr
  real(dp)                                   :: err_sum, tol
  real(dp)                                   :: sum0, sum1, sum2, sum3, sum4
  real(dp)                                   :: sum_d1, sum_d2, sum_d3, sum_d4

  real(dp)                                   :: zdd, zdd_z, zdd_z2, zdd_z3
  real(dp)                                   :: zqq, zqq_z, zqq_z2, zqq_z3
  real(dp)                                   :: zdq, zdq_z, zdq_z2, zdq_z3
  real(dp)                                   :: ztot_z2, ztot_z3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! density-expressions used in many partial derivatives of Helmholtz energy
  !-----------------------------------------------------------------------------
  call density_terms ( eta )

  z2to3 = z2t / z3t
  !-----------------------------------------------------------------------------
  ! hard sphere contribution
  !-----------------------------------------------------------------------------
  ! fhs= ( 3.0_dp*z1*z2/ome + z2**3 /z3/ome2 + (z2**3 /z3/z3-z0)*LOG(ome) ) / rhoPI_6
  if ( .NOT. mod_BMCSL ) then
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t*(3.0_dp-z3)/ome3 ) / PI_6
     zhs_z = ( z0t + 3.0_dp*z1t*z2to3*(1.0_dp+z3)/ome + 6.0_dp*z2*z2t*z2to3/ome2 ) /ome2 / PI_6
     zhs_z2 = m_mean*(  2.0_dp/ome3  + 6.0_dp*z1t*z2to3/z0t*(2.0_dp+z3)/ome4   &
       + 6.0_dp*z2t* z2to3**2 /z0t *(1.0_dp+3.0_dp*z3)/ome5  )
     zhs_z3 = m_mean*(  6.0_dp/ome4  + 18.0_dp*z1t*z2to3/z0t*(3.0_dp+z3)/ome5   &
       + 24.0_dp*z2t*z2to3**2 /z0t*(2.0_dp+3.0_dp*z3)/ome**6  )
  else
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi = z3 + z3 * exp_term
     one_minus_phi = 1.0_dp - phi
     phi_term = z3 + ( 1.0_dp + z3_50 ) * z3 * exp_term
     chi1 = phi_term / one_minus_phi
     chi2_aux = phi_term + ( 2.0_dp + z3_50 ) * z3_50 * z3 * exp_term
     chi2 = chi2_aux / one_minus_phi
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t/z3/ome/one_minus_phi  &
          * ( phi_term/one_minus_phi - one_minus_phi + 1.0_dp/ome ) ) / PI_6
     zhs_z = ( z0t/ome2 + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome3  &
          + z2t**3/z3t**2/ome/one_minus_phi * ( phi/ome + 2.0_dp*chi1/ome  &
          + 2.0_dp*z3/ome2 + 2.0_dp*chi1**2 + chi2 ) ) / PI_6

     phi_z = 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term        ! phi_z = phi_term/z3
     chi1_z = ( 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term ) / one_minus_phi + chi1*phi_z/one_minus_phi     
     chi2_z = ( 1.0_dp + ( 1.0_dp + 7.0_dp*z3_50 + 6.0_dp*z3_50*z3_50 + z3_50**3) * exp_term  ) / one_minus_phi +  &
          chi2_aux / one_minus_phi**2 * phi_z
     zhs_z2 = zhs_z/ome + ( z0t/ome3 + 3.0_dp*z1t*z2t/z3t*(3.0_dp+z3)/ome4  &
          + z2t**3/z3t**2/ome/one_minus_phi  &
          * ( phi_z/one_minus_phi *(phi/ome + 2.0_dp*(chi1/ome+z3/ome2+chi1*chi1) + chi2 )  &
          + phi_z/ome + phi/ome2 + 2.0_dp*(chi1_z/ome+(chi1+1.0_dp)/ome2)  &
          + 4.0_dp*z3/ome3 + 4.0_dp*chi1*chi1_z + chi2_z ) ) / PI_6
     !zhs_z2 = ( 2.0_dp*z0t/ome3 + 6.0_dp*z1t*z2t/z3t*(2.0_dp+z3)/ome4  &
     !     + z2t**3/z3t**2/ome/one_minus_phi  &
     !     * ( (1.0_dp/ome + phi_z/one_minus_phi) *(phi/ome + 2.0_dp*(chi1/ome+z3/ome2+chi1*chi1) + chi2 )  &
     !     + phi_z/ome + phi/ome2 + 2.0_dp*(chi1_z/ome+(chi1+1.0_dp)/ome2)  &
     !     + 4.0_dp*z3/ome3 + 4.0_dp*chi1*chi1_z + chi2_z ) ) / PI_6
     !write (*,*) zhs_z2
     ! *(phi/ome+2.0_dp*chi1/ome+2.0_dp*z3/ome2+2.0_dp*chi1*chi1+ chi2)  &
     !     +chi1/z3/ome+phi/ome2+2.0_dp*chi2/z3/ome+2.0_dp*chi1/ome2+2.0_dp/ome2+4.0_dp*z3/ome3+2.0_dp/z3*chi1*chi2 + chi3/z3 )

     rho = rho - 1.E-7_dp
     z0 = z0t * rho
     z1 = z1t * rho
     z2 = z2t * rho
     z3 = z3t * rho
     ome = 1.0_dp - z3
     ome2 = ome * ome
     ome3 = ome2 * ome
     ome4 = ome2 * ome2
     ome5 = ome4 * ome
     z22 = z2 * z2
     z23 = z2 * z22
     z32 = z3 * z3
     z33 = z3 * z32
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi = z3 + z3 * exp_term
     one_minus_phi = 1.0_dp - phi
     phi_term = z3 + ( 1.0_dp + z3_50 ) * z3 * exp_term
     chi1 = phi_term / one_minus_phi
     chi2_aux = phi_term + ( 2.0_dp + z3_50 ) * z3_50 * z3 * exp_term
     chi2 = chi2_aux / one_minus_phi
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t/z3/ome/one_minus_phi  &
          * ( phi_term/one_minus_phi - one_minus_phi + 1.0_dp/ome ) ) / PI_6
     zhs_z = ( z0t/ome2 + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome3  &
          + z2t**3/z3t**2/ome/one_minus_phi * ( phi/ome + 2.0_dp*chi1/ome  &
          + 2.0_dp*z3/ome2 + 2.0_dp*chi1**2 + chi2 ) ) / PI_6

     phi_z = 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term        ! phi_z = phi_term/z3
     chi1_z = ( 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term ) / one_minus_phi + chi1*phi_z/one_minus_phi     
     chi2_z = ( 1.0_dp + ( 1.0_dp + 7.0_dp*z3_50 + 6.0_dp*z3_50*z3_50 + z3_50**3) * exp_term  ) / one_minus_phi +  &
          chi2_aux / one_minus_phi**2 * phi_z
     zhs_z2 = zhs_z/ome + ( z0t/ome3 + 3.0_dp*z1t*z2t/z3t*(3.0_dp+z3)/ome4  &
          + z2t**3/z3t**2/ome/one_minus_phi  &
          * ( phi_z/one_minus_phi *(phi/ome + 2.0_dp*(chi1/ome+z3/ome2+chi1*chi1) + chi2 )  &
          + phi_z/ome + phi/ome2 + 2.0_dp*(chi1_z/ome+(chi1+1.0_dp)/ome2)  &
          + 4.0_dp*z3/ome3 + 4.0_dp*chi1*chi1_z + chi2_z ) ) / PI_6
     zhs1_z2 = zhs_z2

     rho = rho + 2.E-7_dp
     z0 = z0t * rho
     z1 = z1t * rho
     z2 = z2t * rho
     z3 = z3t * rho
     ome = 1.0_dp - z3
     ome2 = ome * ome
     ome3 = ome2 * ome
     ome4 = ome2 * ome2
     ome5 = ome4 * ome
     z22 = z2 * z2
     z23 = z2 * z22
     z32 = z3 * z3
     z33 = z3 * z32
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi = z3 + z3 * exp_term
     one_minus_phi = 1.0_dp - phi
     phi_term = z3 + ( 1.0_dp + z3_50 ) * z3 * exp_term
     chi1 = phi_term / one_minus_phi
     chi2_aux = phi_term + ( 2.0_dp + z3_50 ) * z3_50 * z3 * exp_term
     chi2 = chi2_aux / one_minus_phi
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t/z3/ome/one_minus_phi  &
          * ( phi_term/one_minus_phi - one_minus_phi + 1.0_dp/ome ) ) / PI_6
     zhs_z = ( z0t/ome2 + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome3  &
          + z2t**3/z3t**2/ome/one_minus_phi * ( phi/ome + 2.0_dp*chi1/ome  &
          + 2.0_dp*z3/ome2 + 2.0_dp*chi1**2 + chi2 ) ) / PI_6

     phi_z = 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term        ! phi_z = phi_term/z3
     chi1_z = ( 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term ) / one_minus_phi + chi1*phi_z/one_minus_phi     
     chi2_z = ( 1.0_dp + ( 1.0_dp + 7.0_dp*z3_50 + 6.0_dp*z3_50*z3_50 + z3_50**3) * exp_term  ) / one_minus_phi +  &
          chi2_aux / one_minus_phi**2 * phi_z
     zhs_z2 = zhs_z/ome + ( z0t/ome3 + 3.0_dp*z1t*z2t/z3t*(3.0_dp+z3)/ome4  &
          + z2t**3/z3t**2/ome/one_minus_phi  &
          * ( phi_z/one_minus_phi *(phi/ome + 2.0_dp*(chi1/ome+z3/ome2+chi1*chi1) + chi2 )  &
          + phi_z/ome + phi/ome2 + 2.0_dp*(chi1_z/ome+(chi1+1.0_dp)/ome2)  &
          + 4.0_dp*z3/ome3 + 4.0_dp*chi1*chi1_z + chi2_z ) ) / PI_6
     zhs2_z2 = zhs_z2

     rho = rho - 1.E-7_dp
     z0 = z0t * rho
     z1 = z1t * rho
     z2 = z2t * rho
     z3 = z3t * rho
     ome = 1.0_dp - z3
     ome2 = ome * ome
     ome3 = ome2 * ome
     ome4 = ome2 * ome2
     ome5 = ome4 * ome
     z22 = z2 * z2
     z23 = z2 * z22
     z32 = z3 * z3
     z33 = z3 * z32
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi = z3 + z3 * exp_term
     one_minus_phi = 1.0_dp - phi
     phi_term = z3 + ( 1.0_dp + z3_50 ) * z3 * exp_term
     chi1 = phi_term / one_minus_phi
     chi2_aux = phi_term + ( 2.0_dp + z3_50 ) * z3_50 * z3 * exp_term
     chi2 = chi2_aux / one_minus_phi
     zhs = ( z3t*z0/ome + 3.0_dp*z1t*z2/ome2 + z22*z2t/z3/ome/one_minus_phi  &
          * ( phi_term/one_minus_phi - one_minus_phi + 1.0_dp/ome ) ) / PI_6
     zhs_z = ( z0t/ome2 + 3.0_dp*z1t*z2t/z3t*(1.0_dp+z3)/ome3  &
          + z2t**3/z3t**2/ome/one_minus_phi * ( phi/ome + 2.0_dp*chi1/ome  &
          + 2.0_dp*z3/ome2 + 2.0_dp*chi1**2 + chi2 ) ) / PI_6

     phi_z = 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term        ! phi_z = phi_term/z3
     chi1_z = ( 1.0_dp + ( 1.0_dp + z3_50 ) * exp_term ) / one_minus_phi + chi1*phi_z/one_minus_phi     
     chi2_z = ( 1.0_dp + ( 1.0_dp + 7.0_dp*z3_50 + 6.0_dp*z3_50*z3_50 + z3_50**3) * exp_term  ) / one_minus_phi +  &
          chi2_aux / one_minus_phi**2 * phi_z
     zhs_z2 = zhs_z/ome + ( z0t/ome3 + 3.0_dp*z1t*z2t/z3t*(3.0_dp+z3)/ome4  &
          + z2t**3/z3t**2/ome/one_minus_phi  &
          * ( phi_z/one_minus_phi *(phi/ome + 2.0_dp*(chi1/ome+z3/ome2+chi1*chi1) + chi2 )  &
          + phi_z/ome + phi/ome2 + 2.0_dp*(chi1_z/ome+(chi1+1.0_dp)/ome2)  &
          + 4.0_dp*z3/ome3 + 4.0_dp*chi1*chi1_z + chi2_z ) ) / PI_6
     zhs_z3 = ( zhs2_z2 - zhs1_z2 ) / 2.E-7_dp / z3t
  end if
  

  !-----------------------------------------------------------------------------
  ! chain term
  !-----------------------------------------------------------------------------
  gij_z_t1 = 1.0_dp / ome2
  gij_z_t2 = gij_t2 * ( 1.0_dp/z3t + 2.0_dp*rho/ome )
  gij_z_t3 = gij_t3 * (2.0_dp+z3) / z3t / ome
  gij_z2_t1 = 2.0_dp / ome3
  gij_z2_t2 = 6.0_dp*z2to3/ome4 *(2.0_dp+z3)
  gij_z2_t3 = 4.0_dp*z2to3**2/ome5  *(1.0_dp+4.0_dp*z3+z32)
  gij_z3_t1 = 6.0_dp / ome4
  gij_z3_t2 = 18.0_dp*z2to3/ome5 *(3.0_dp+z3)
  gij_z3_t3 = 12.0_dp*(z2to3/ome3)**2  *(3.0_dp+6.0_dp*z3+z32)
  gij_z4_t1 = 24.0_dp / ome5
  gij_z4_t2 = 72.0_dp*z2to3/ome3**2 *(4.0_dp+z3)
  gij_z4_t3 = 48.0_dp*z2to3**2 /ome**7  *(6.0_dp+8.0_dp*z3+z3*z3)
  do i = 1, ncomp
     gij_z(i,i)  = gij_z_t1  + dij_ab(i,i)*( gij_z_t2  + dij_ab(i,i) *gij_z_t3  )
     gij_z2(i,i) = gij_z2_t1 + dij_ab(i,i)*( gij_z2_t2 + dij_ab(i,i) *gij_z2_t3 )
     gij_z3(i,i) = gij_z3_t1 + dij_ab(i,i)*( gij_z3_t2 + dij_ab(i,i) *gij_z3_t3 )
     gij_z4(i,i) = gij_z4_t1 + dij_ab(i,i)*( gij_z4_t2 + dij_ab(i,i) *gij_z4_t3 )
  end do

  zhc = 0.0_dp
  zhc_z = 0.0_dp
  zhc_z2 = 0.0_dp
  zhc_z3 = 0.0_dp
  do i= 1, ncomp
     abbrev = x(i) * (1.0_dp-mseg(i)) * z3
     zhc = zhc + abbrev / gij(i,i)* gij_z(i,i)
     zhc_z = zhc_z + abbrev *( - (gij_z(i,i)/gij(i,i))**2   &
          + gij_z(i,i)/gij(i,i)/z3 + gij_z2(i,i)/gij(i,i) )
     zhc_z2 = zhc_z2 + x(i)*(1.0_dp-mseg(i))  &
          *( 2.0_dp*z3*(gij_z(i,i)/gij(i,i))**3   &
          -2.0_dp*(gij_z(i,i)/gij(i,i))**2   &
          -3.0_dp*z3/gij(i,i)**2 *gij_z(i,i)*gij_z2(i,i)  &
          +2.0_dp/gij(i,i)*gij_z2(i,i) +z3/gij(i,i)*gij_z3(i,i) )
     zhc_z3 = zhc_z3 + x(i)*(1.0_dp-mseg(i)) *( 6.0_dp*(gij_z(i,i)/gij(i,i))**3   &
          -6.0_dp*z3*(gij_z(i,i)/gij(i,i))**4   &
          +12.0_dp*z3/gij(i,i)**3 *gij_z(i,i)**2 *gij_z2(i,i)  &
          -9.0_dp/gij(i,i)**2 *gij_z(i,i)*gij_z2(i,i) +3.0_dp/gij(i,i)*gij_z3(i,i)  &
          -3.0_dp*z3*(gij_z2(i,i)/gij(i,i))**2   &
          -4.0_dp*z3/gij(i,i)**2 *gij_z(i,i)*gij_z3(i,i)  &
          +z3/gij(i,i)*gij_z4(i,i) )
  end do


  !-----------------------------------------------------------------------------
  ! p : PC-SAFT dispersion contribution
  !     note: edI1dz is equal to d(eta*I1)/d(eta), analogous for edI2dz
  !-----------------------------------------------------------------------------
  edI1dz = z3 * I1_z + I1
  edI2dz = z3 * I2_z + I2
  edI1d2 = z3 * I1_z2 + 2.0_dp * I1_z
  edI2d2 = z3 * I2_z2 + 2.0_dp * I2_z

  edI1d3 = 0.0_dp
  edI2d3 = 0.0_dp
  edI1d4 = 0.0_dp
  edI2d4 = 0.0_dp
  do  m = 2, 6
     edI1d3 = edI1d3 + apar(m) * real( (m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
     edI2d3 = edI2d3 + bpar(m) * real( (m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
     edI1d4 = edI1d4 + apar(m) * real( (m+1)*m*(m-1)*(m-2), KIND=dp ) * z3**(m-3)
     edI2d4 = edI2d4 + bpar(m) * real( (m+1)*m*(m-1)*(m-2), KIND=dp ) * z3**(m-3)
  end do

  c4_con= 6.0_dp*c2_con*c3_con/c1_con -6.0_dp*c2_con**3 /c1_con**2   &
       - c1_con*c1_con  &
       *( m_mean*(-48.0_dp*z32 +336.0_dp*z3+432.0_dp)/ome**7   &
       + (1.0_dp - m_mean)  &
       *(24.0_dp*z3**5 +240.0_dp*z3**4 -1920.0_dp*z33   &
       +4800.0_dp*z32 -5280.0_dp*z3+2208.0_dp) /(ome*(2.0_dp-z3))**5  )
  c5_con= 6.0_dp*c3_con**2 /c1_con - 36.0_dp*c2_con**2 /c1_con**2 *c3_con  &
       + 8.0_dp*c2_con/c1_con*c4_con+24.0_dp*c2_con**4 /c1_con**3   &
       - c1_con*c1_con  &
       *( m_mean*(-240.0_dp*z32 +1920.0_dp*z3+3360.0_dp)/ome**8   &
       + (1.0_dp - m_mean)  &
       *(-120.0_dp*z3**6 -1440.0_dp*z3**5 +14400.0_dp*z3**4   &
       -48000.0_dp*z33 +79200.0_dp*z32  -66240.0_dp*z3+22560.0_dp)  &
       /(ome*(2.0_dp-z3))**6  )

  zdsp = rho*edI1dz*order1  &
          + rho*order2*m_mean*(c2_con*I2*z3 + c1_con*edI2dz)
  zdsp_z = zdsp/z3 + rho*edI1d2*order1  &
          + rho*order2*m_mean*(c3_con*I2*z3 + 2.0_dp*c2_con*edI2dz + c1_con*edI2d2)
  zdsp_z2 = -2.0_dp*zdsp/z3/z3 + 2.0_dp*zdsp_z/z3  &
          + rho*edI1d3*order1 + rho*order2*m_mean*( c4_con*I2*z3  &
          + 3.0_dp*c3_con*edI2dz + 3.0_dp*c2_con*edI2d2 + c1_con*edI2d3 )
  zdsp_z3 = 6.0_dp*zdsp/z3**3  - 6.0_dp*zdsp_z/z3/z3  &
          + 3.0_dp*zdsp_z2/z3 + rho*edI1d4*order1  &
          + rho*order2*m_mean*( c5_con*I2*z3  &
          + 4.0_dp*c4_con*edI2dz + 6.0_dp*c3_con*edI2d2  &
          + 4.0_dp*c2_con*edI2d3 + c1_con*edI2d4 )



  !-----------------------------------------------------------------------------
  ! p: TPT-1-association accord. to Chapman et al.
  !-----------------------------------------------------------------------------
  zhb = 0.0_dp
  zhb_z = 0.0_dp
  zhb_z2 = 0.0_dp
  zhb_z3 = 0.0_dp
  assoc = .false.
  do i = 1,ncomp
     if ( nhb_typ(i) /= 0 ) assoc = .true.
  end do
  if (assoc) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij_z(i,j) = gij_z_t1 + dij_ab(i,j)*( gij_z_t2 + dij_ab(i,j) *gij_z_t3 )
           gij_z2(i,j) = gij_z2_t1 + dij_ab(i,j)*( gij_z2_t2 + dij_ab(i,j) *gij_z2_t3 )
           gij_z3(i,j) = gij_z3_t1 + dij_ab(i,j)*( gij_z3_t2 + dij_ab(i,j) *gij_z3_t3 )
           gij_z4(i,j) = gij_z4_t1 + dij_ab(i,j)*( gij_z4_t2 + dij_ab(i,j) *gij_z4_t3 )
           gij_z(j,i) = gij_z(i,j)
           gij_z2(j,i) = gij_z2(i,j)
           gij_z3(j,i) = gij_z3(i,j)
           gij_z4(j,i) = gij_z4(i,j)
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 delta(i,j,ii,jj) = gij(i,j)    * ass_d(i,j,ii,jj)
                 dq_dz(i,j,ii,jj) = gij_z(i,j)  * ass_d(i,j,ii,jj)
                 dq_d2(i,j,ii,jj) = gij_z2(i,j) * ass_d(i,j,ii,jj)
                 dq_d3(i,j,ii,jj) = gij_z3(i,j) * ass_d(i,j,ii,jj)
                 dq_d4(i,j,ii,jj) = gij_z4(i,j) * ass_d(i,j,ii,jj)
              end do
           end do
        end do
     end do

     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-5_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx(i,ii) = 1.0_dp
           mx_z(i,ii) = 0.0_dp
           mx_z2(i,ii) = 0.0_dp
           mx_z3(i,ii) = 0.0_dp
           mx_z4(i,ii) = 0.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .and. ass_cnt <= max_eval )
        ass_cnt = ass_cnt + 1
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum0 = 0.0_dp
              sum1 = 0.0_dp
              sum2 = 0.0_dp
              sum3 = 0.0_dp
              sum4 = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    sum0 =sum0 +x(j)*nhb_no(j,jj)*     mx(j,jj)* delta(i,j,ii,jj)
                    sum1 =sum1 +x(j)*nhb_no(j,jj)*(    mx(j,jj)* dq_dz(i,j,ii,jj)  &
                         +      mx_z(j,jj)* delta(i,j,ii,jj))
                    sum2 =sum2 +x(j)*nhb_no(j,jj)*(    mx(j,jj)* dq_d2(i,j,ii,jj)  &
                         + 2.0_dp*mx_z(j,jj)* dq_dz(i,j,ii,jj)  &
                         +      mx_z2(j,jj)* delta(i,j,ii,jj))
                    sum3 =sum3 +x(j)*nhb_no(j,jj)*(    mx(j,jj)* dq_d3(i,j,ii,jj)  &
                         + 3.0_dp*mx_z(j,jj)* dq_d2(i,j,ii,jj)  &
                         + 3.0_dp*mx_z2(j,jj)* dq_dz(i,j,ii,jj)  &
                         +      mx_z3(j,jj)* delta(i,j,ii,jj))
                    sum4 =sum4 + x(j)*nhb_no(j,jj)*(   mx(j,jj)* dq_d4(i,j,ii,jj)  &
                         + 4.0_dp*mx_z(j,jj)* dq_d3(i,j,ii,jj)  &
                         + 6.0_dp*mx_z2(j,jj)* dq_d2(i,j,ii,jj)  &
                         + 4.0_dp*mx_z3(j,jj)* dq_dz(i,j,ii,jj)  &
                         +      mx_z4(j,jj)* delta(i,j,ii,jj))
                 end do
              end do
              mx_itr(i,ii)= 1.0_dp / (1.0_dp + sum0 * rho)
              mx_z_itr(i,ii)= -(mx_itr(i,ii)*mx_itr(i,ii))* (sum0/z3t +sum1*rho)
              mx_z2_itr(i,ii)= + 2.0_dp/mx_itr(i,ii)*mx_z_itr(i,ii)*mx_z_itr(i,ii)  &
                   - (mx_itr(i,ii)*mx_itr(i,ii)) * (2.0_dp/z3t*sum1 + rho*sum2)
              mx_z3_itr(i,ii)= - 6.0_dp/mx_itr(i,ii)**2 *mx_z_itr(i,ii)**3   &
                   + 6.0_dp/mx_itr(i,ii)*mx_z_itr(i,ii)*mx_z2_itr(i,ii) - mx_itr(i,ii)*mx_itr(i,ii)  &
                   * (3.0_dp/z3t*sum2 + rho*sum3)
              mx_z4_itr(i,ii)= 24.0_dp/mx_itr(i,ii)**3 *mx_z_itr(i,ii)**4   &
                   -36.0_dp/mx_itr(i,ii)**2 *mx_z_itr(i,ii)**2 *mx_z2_itr(i,ii)  &
                   +6.0_dp/mx_itr(i,ii)*mx_z2_itr(i,ii)**2   &
                   +8.0_dp/mx_itr(i,ii)*mx_z_itr(i,ii)*mx_z3_itr(i,ii) - mx_itr(i,ii)**2   &
                   *(4.0_dp/z3t*sum3 + rho*sum4)
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + abs(mx_itr(i,ii) - mx(i,ii))  &
                   + abs(mx_z_itr(i,ii) - mx_z(i,ii)) + abs(mx_z2_itr(i,ii) - mx_z2(i,ii))
              mx(i,ii)     = mx_itr(i,ii)*hb_damp +     mx(i,ii) * hb_damp_inv
              mx_z(i,ii) = mx_z_itr(i,ii)*hb_damp + mx_z(i,ii) * hb_damp_inv
              mx_z2(i,ii) = mx_z2_itr(i,ii)*hb_damp + mx_z2(i,ii) * hb_damp_inv
              mx_z3(i,ii) = mx_z3_itr(i,ii)*hb_damp + mx_z3(i,ii) * hb_damp_inv
              mx_z4(i,ii) = mx_z4_itr(i,ii)*hb_damp + mx_z4(i,ii) * hb_damp_inv
           end do
        end do
     end do

     if ( ass_cnt >= max_eval .and. err_sum > sqrt(tol) ) then
        write (*,'(a,2G15.7)') 'f_rho_4: Max_eval violated (mx) Err_Sum= ',err_sum,tol
        ! stop
     end if


     ! --- calculate the hydrogen-bonding contribution -------------------------
     do i = 1, ncomp
        sum_d1 = 0.0_dp
        sum_d2 = 0.0_dp
        sum_d3 = 0.0_dp
        sum_d4 = 0.0_dp
        do ii = 1, nhb_typ(i)
           sum_d1 = sum_d1 +nhb_no(i,ii)* mx_z(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)
           sum_d2 = sum_d2 +nhb_no(i,ii)*(mx_z2(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)  &
                -(mx_z(i,ii)/mx(i,ii))**2 )
           sum_d3 = sum_d3 +nhb_no(i,ii)*(mx_z3(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)  &
                -3.0_dp/mx(i,ii)**2 *mx_z(i,ii)*mx_z2(i,ii) + 2.0_dp*(mx_z(i,ii)/mx(i,ii))**3 )
           sum_d4 = sum_d4 +nhb_no(i,ii)*(mx_z4(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)  &
                -4.0_dp/mx(i,ii)**2 *mx_z(i,ii)*mx_z3(i,ii)  &
                + 12.0_dp/mx(i,ii)**3 *mx_z(i,ii)**2 *mx_z2(i,ii)  &
                - 3.0_dp/mx(i,ii)**2 *mx_z2(i,ii)**2  - 6.0_dp*(mx_z(i,ii)/mx(i,ii))**4 )
        end do
        zhb = zhb + x(i) * z3 * sum_d1
        zhb_z = zhb_z + x(i) * z3 * sum_d2
        zhb_z2 = zhb_z2 + x(i) * z3 * sum_d3
        zhb_z3 = zhb_z3 + x(i) * z3 * sum_d4
     end do
     zhb_z = zhb_z + zhb/z3
     zhb_z2 = zhb_z2 + 2.0_dp/z3*zhb_z-2.0_dp/z32 *zhb
     zhb_z3 = zhb_z3 - 6.0_dp/z32 *zhb_z + 3.0_dp/z3*zhb_z2 + 6.0_dp/z33 * zhb
  end if


  !-----------------------------------------------------------------------------
  ! p: polar terms
  !-----------------------------------------------------------------------------
  call fdd_rho_4 ( zdd, zdd_z, zdd_z2, zdd_z3 )
  call fqq_rho_4 ( zqq, zqq_z, zqq_z2, zqq_z3 )
  call fdq_rho_4 ( zdq, zdq_z, zdq_z2, zdq_z3 )


  !-----------------------------------------------------------------------------
  ! compressibility factor z and total p
  ! as well as derivatives d(z)/d(z3) and d(p)/d(z3) with unit [Pa]
  !-----------------------------------------------------------------------------
  z_res = zhs + zhc + zdsp + zhb + zdd + zqq + zdq
  ztot_z = zhs_z + zhc_z + zdsp_z + zhb_z + zdd_z + zqq_z + zdq_z
  ztot_z2 = zhs_z2 + zhc_z2 + zdsp_z2 + zhb_z2 + zdd_z2 +zqq_z2 + zdq_z2
  ztot_z3 = zhs_z3 + zhc_z3 + zdsp_z3 + zhb_z3 + zdd_z3 +zqq_z3 + zdq_z3

  ! f_r = z_res / rho
  ! f_r2 = ( ztot_z*z3t - f_r ) / rho

  ! pcalc = ztot * rho * kT
  ! pcalc_z = ( ztot_z + ztot/z3 ) * rho * kT
  pcalc_z2 = ( ztot_z2 * rho + 2.0_dp *ztot_z / z3t ) * kT
  pcalc_z3 = ( ztot_z3 * rho + 3.0_dp *ztot_z2 / z3t ) * kT

end subroutine f_rho_4


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine Iterate_Association
!
! calculates the Helmholtz energy density f = A/(VkT) = A/(NkT)*rho
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine Iterate_Association

  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  integer                                         :: i, j, ii, jj
  integer                                         :: ass_cnt, max_eval

  real(dp), dimension(ncomp,nsite)                :: mx_itr
  real(dp)                                        :: err_sum, sum0, tol
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! w: TPT-1-association accord. to Chapman et al.
  !-----------------------------------------------------------------------------
  if ( assoc ) then

     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-4_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-5_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     mx(:,:) = 1.0_dp

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum0 = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    sum0 = sum0 + rhoi(j)*nhb_no(j,jj)* mx(j,jj)* gij(i,j)*ass_d(i,j,ii,jj)
                 end do
              end do
              mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + sum0 )
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS( mx_itr(i,ii) - mx(i,ii) )
              mx(i,ii) = mx_itr(i,ii) * hb_damp + mx(i,ii) * hb_damp_inv
           end do
        end do

        if ( ass_cnt == max_eval ) then
           write (*,'(a,2G15.7)') 'Iterate_Association: Max_eval violated (mx) Err_Sum= ', err_sum, tol
           exit
        end if

     end do

  end if

end subroutine Iterate_Association


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine F_density
!
! calculates the Helmholtz energy density f = A/(VkT) = A/(NkT)*rho
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_density ( rhoi_in, f_dens )

  use EOS_CONSTANTS
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi_in
  real(dp), intent(out)                           :: f_dens
  !-----------------------------------------------------------------------------
  integer                                         :: i, j, ii, jj, l, m
  integer                                         :: ass_cnt, max_eval

  real(dp)                                        :: wig
  real(dp)                                        :: whs
  real(dp)                                        :: z3_50, phi, one_minus_phi

  real(dp)                                        :: whc
  real(dp)                                        :: gij_t2_rho, gij_t3_rho

  real(dp)                                        :: wdsp

  real(dp)                                        :: whb
  real(dp), dimension(ncomp,nsite)                :: mx_itr
  real(dp)                                        :: err_sum, sum0, tol

  real(dp)                                        :: wdd
  real(dp)                                        :: wdd2, wdd3

  real(dp)                                        :: wqq
  real(dp)                                        :: wqq2, wqq3

  real(dp)                                        :: wdq
  real(dp)                                        :: wdq2, wdq3
  !-----------------------------------------------------------------------------

  rhoi( 1:ncomp ) = rhoi_in(1:ncomp )

  !-----------------------------------------------------------------------------
  ! density-expressions used in many partial derivatives of Helmholtz energy
  !-----------------------------------------------------------------------------
  rho = sum ( rhoi( 1:ncomp ) )
  z0 = z0t * rho
  z1 = z1t * rho
  z2 = z2t * rho
  z3 = z3t * rho

  ome = 1.0_dp - z3
  ome2 = ome * ome
  ome3 = ome2 * ome
  ome4 = ome2 * ome2
  ome5 = ome4 * ome

  z22 = z2 * z2
  z23 = z2 * z22
  z32 = z3 * z3
  z33 = z3 * z32

  do m = 0, 6
     z3_m_m0(m) = z3**m
  end do



  !-----------------------------------------------------------------------------
  ! w : ideal gas contribution
  !-----------------------------------------------------------------------------
  wig = sum( rhoi(1:ncomp) * ( log( rhoi(1:ncomp) ) - 1.0_dp ) )


  !-----------------------------------------------------------------------------
  ! w : hard sphere contribution
  !-----------------------------------------------------------------------------
  if ( .NOT. mod_BMCSL ) then
     whs = ( 3.0_dp*z1*z2/ome + z22*z2t/z3t/ome2 + (z2*(z2t/z3t)**2-z0)*LOG(ome) ) / PI_6
  else
     z3_50 = 50.0_dp * z3
     phi = z3 + z3 * jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     one_minus_phi = 1.0_dp - phi
     whs = ( 3.0_dp*z1*z2/ome + z22*z2t/z3t/ome/one_minus_phi + (z2*(z2t/z3t)**2-z0)*LOG(ome) ) / PI_6
  end if


  !-----------------------------------------------------------------------------
  ! w : chain term
  !-----------------------------------------------------------------------------
  gij_t1 = 1.0_dp/ome
  gij_t2 = 3.0_dp*z2t/ome2
  gij_t3 = 2.0_dp*z2t*z2/ome3
  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3
  do i = 1, ncomp
     gij(i,i) = gij_t1 + dij_ab(i,i)*( gij_t2_rho + dij_ab(i,i) *gij_t3_rho )
  end do

  whc = 0.0_dp
  do i = 1, ncomp
     whc = whc + rhoi(i) * (1.0_dp-mseg(i)) * LOG( gij(i,i) )
  end do

  !-----------------------------------------------------------------------------
  ! w : PC-SAFT dispersion contribution
  !-----------------------------------------------------------------------------
  I1 = sum( apar(0:6) * z3_m_m0(0:6) )
  I2 = sum( bpar(0:6) * z3_m_m0(0:6) )

  ome_t = 1.0_dp / (ome*(2.0_dp-z3))
  ome_t_2 = ome_t * ome_t
  c1_con= 1.0_dp/ (  1.0_dp + m_mean*(8.0_dp*z3-2.0_dp*z32 )/ome4   &
       + (1.0_dp - m_mean)*(20.0_dp*z3-27.0_dp*z32 +12.0_dp*z33 -2.0_dp*z32*z32 ) * ome_t_2   )

  wdsp = rho**2 * ( I1 * order1 + m_mean * c1_con * I2 * order2 )


  !-----------------------------------------------------------------------------
  ! w: TPT-1-association accord. to Chapman et al.
  !-----------------------------------------------------------------------------
  whb = 0.0_dp

  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij(i,j) = gij_t1 + dij_ab(i,j)*( gij_t2_rho + dij_ab(i,j) *gij_t3_rho )
           gij(j,i) = gij(i,j)
        end do
     end do

     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-4_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-5_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx(i,ii) = 1.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum0 = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    sum0 = sum0 + rhoi(j)*nhb_no(j,jj)* mx(j,jj)* gij(i,j)*ass_d(i,j,ii,jj)
                 end do
              end do
              mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + sum0 )
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS( mx_itr(i,ii) - mx(i,ii) )
              mx(i,ii) = mx_itr(i,ii) * hb_damp + mx(i,ii) * hb_damp_inv
              if ( mx(i,ii) < 0.0_dp ) mx(i,ii) = 1.E-20_dp
           end do
        end do

        if ( ass_cnt == max_eval ) then
           write (*,'(a,2G15.7)') 'F_density: Max_eval violated (mx) Err_Sum= ', err_sum, tol
           exit
        end if

     end do

     ! --- calculate the hydrogen-bonding contribution -------------------------
     whb = 0.0_dp
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           whb = whb + rhoi(i) * nhb_no(i,ii) * ( 0.5_dp - 0.5_dp * mx(i,ii) + LOG( mx(i,ii) ) )
        end do
     end do

  end if


  !-----------------------------------------------------------------------------
  ! w: polar terms
  !-----------------------------------------------------------------------------
  wdd = 0.0_dp
  wqq = 0.0_dp
  wdq = 0.0_dp

  ! ------ dipole-dipole term --------------------------------------------------
  if ( dipole ) then

     do i = 1, ncomp
        if ( abs( my_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j) ) > machine_eps ) then
                 Idd2(i,j) = sum( ( ddp2(i,j,0:4) + eps_ij(i,j) * ddp4(i,j,0:4) ) * z3_m_m0(0:4) )
                 do l = 1, ncomp
                    Idd3(i,j,l) = sum( ddp3(i,j,l,0:4) * z3_m_m0(0:4) )
                 end do
              end if
           end do
        end if
     end do

     wdd2 = 0.0_dp
     wdd3 = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_factor(i)*rhoi(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j)*rhoi(j) ) > machine_eps ) then
                 wdd2 = wdd2 + rhoi(i)*rhoi(j)*pi_dd(i,j) * Idd2(i,j)
                 do l = 1, ncomp
                    wdd3 = wdd3 + rhoi(i)*rhoi(j)*rhoi(l)*psi_dd(i,j,l) * Idd3(i,j,l)
                 end do
              end if
           end do
        end if
     end do

     if ( abs( wdd2 ) > 1.E-50_dp ) then

        wdd = wdd2 *wdd2 / ( wdd2 -wdd3 )

     end if

  end if
  
  ! ------ quadrupole-quadrupole term ------------------------------------------
  if ( qudpole ) then

     do i = 1, ncomp
        if ( abs( Q_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j) ) > machine_eps ) then
                 Iqq2(i,j) = sum( ( qqp2(i,j,0:4) + eps_ij(i,j) * qqp4(i,j,0:4) ) * z3_m_m0(0:4) )
                 do l = 1, ncomp
                    Iqq3(i,j,l) = sum( qqp3(i,j,l,0:4) * z3_m_m0(0:4) )
                 end do
              end if
           end do
        end if
     end do

     wqq2 = 0.0_dp
     wqq3 = 0.0_dp
     do i = 1, ncomp
        if ( abs( Q_factor(i)*rhoi(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j)*rhoi(j) ) > machine_eps ) then
                 wqq2 = wqq2 + rhoi(i)*rhoi(j)*pi_qq(i,j) * Iqq2(i,j)
                 do l = 1, ncomp
                    wqq3 = wqq3 + rhoi(i)*rhoi(j)*rhoi(l)*psi_qq(i,j,l) * Iqq3(i,j,l)
                 end do
              end if
           end do
        end if
     end do

     if ( abs( wqq2 ) > 1.E-50_dp ) then

        wqq = wqq2 *wqq2 / ( wqq2 -wqq3 )

     end if

  end if
  
  ! ------ dipole-quadrupole term ----------------------------------------------
  if ( dipole_quad ) then

     do i = 1, ncomp
        do j = 1, ncomp
           Idq2(i,j) = sum( ( dqp2(i,j,0:4) + eps_ij(i,j) * dqp4(i,j,0:4) ) * z3_m_m0(0:4) )
           do l = 1, ncomp
              Idq3(i,j,l) = sum( dqp3(i,j,l,0:4) * z3_m_m0(0:4) )
           end do
        end do
     end do

     wdq2 = 0.0_dp
     wdq3 = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_fac_dq(i)*rhoi(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_fac_dq(j)*rhoi(j) ) > machine_eps ) then
                 wdq2 = wdq2 + rhoi(i) * rhoi(j) * pi_dq(i,j) * Idq2(i,j)
                 do l = 1, ncomp
                    wdq3 = wdq3 + rhoi(i) * rhoi(j) * rhoi(l) * psi_dq(i,j,l) * Idq3(i,j,l)
                 end do
              end if
           end do
        end if
     end do


     if ( abs( wdq2 ) > 1.E-50_dp ) then

        wdq = wdq2 *wdq2 / ( wdq2 -wdq3 )

     end if

  end if


  !-----------------------------------------------------------------------------
  ! Helmholtz energy density f_dens = A/(VkT) = A/(NkT)*rho
  !-----------------------------------------------------------------------------
  f_dens = wig + whs + whc + wdsp + whb + wdd + wqq + wdq

end subroutine F_density


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_density_rhok
!
! calculate derivatives of w = A/(VkT) = a*rho to component density: d(w)/d(rho_k)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_density_rhok ( myres )

  use PARAMETERS, only: PI
  USE EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: lij_correction
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(out)         :: myres
  !-----------------------------------------------------------------------------
  integer                                         :: i, j, k, ii, jj

  real(dp)                                        :: whs_rk
  real(dp)                                        :: z2to3, z2to3_2, z2to3_3
  real(dp)                                        :: hs_term1, hs_term2, hs_term3
  real(dp)                                        :: z3_50, exp_term
  real(dp)                                        :: phi_0, phi_1
  real(dp)                                        :: one_minus_phi
  real(dp)                                        :: lamda, log_ome

  real(dp)                                        :: whc_rk
  real(dp)                                        :: gij_rk_t1, gij_rk_t2, gij_rk_t3
  real(dp), dimension(ncomp,ncomp)                :: gij_rk

  real(dp)                                        :: m_rk
  real(dp)                                        :: wdsp_rk
  real(dp)                                        :: I1_rk, I2_rk
  real(dp)                                        :: I1_m, I2_m
  real(dp)                                        :: abbrev
  real(dp)                                        :: ord1_rk, ord2_rk
  real(dp)                                        :: c1_rk
  real(dp)                                        :: c_factor

  real(dp)                                        :: whb_rk
  real(dp)                                        :: factor_i, factor_ij

  real(dp)                                        :: wdd_rk
  real(dp)                                        :: wdd2_rk, wdd3_rk
  real(dp)                                        :: wqq_rk
  real(dp)                                        :: wqq2_rk, wqq3_rk
  real(dp)                                        :: wdq_rk
  real(dp)                                        :: wdq2_rk, wdq3_rk
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! quantities independent of k. Index k cycles for the derivative  d(w)/d(rho_k)
  !-----------------------------------------------------------------------------

  I1_m = sum(  ap_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I1 to mean segment number m_mean
  I2_m = sum(  bp_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I2 to mean segment number m_mean
  c_factor = (8.0_dp*z3-2.0_dp*z32)/ome4  &
       - (-2.0_dp*z3**4 +12.0_dp*z33 -27.0_dp*z32+20.0_dp*z3) *ome_t_2
  c_factor = c1_con * c1_con * c_factor

  z2to3 = z2t / z3t
  z2to3_2 = z2to3 * z2to3
  z2to3_3 = z2to3_2 * z2to3
  log_ome = LOG(ome)
  if ( .NOT. mod_BMCSL ) then
     hs_term1 = 3.0_dp * z2/ome
     hs_term2 = 3.0_dp * ( z1/ome + z2*z2to3/ome2 + z2to3_2*log_ome )
     hs_term3 = 3.0_dp*z1*z2/ome2 + z2*z2to3_2 *(3.0_dp*z3-1.0_dp)/ome3  &
          + (z0-z2*z2to3_2)/ome - 2.0_dp*z2to3_3*log_ome
  else
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi_0 = 1.0_dp + exp_term
     phi_1 = phi_0 + z3_50 * exp_term
     one_minus_phi = 1.0_dp - phi_0 * z3
     lamda = z2*z2to3/ome/one_minus_phi
     hs_term1 = 3.0_dp * z2/ome
     hs_term2 = 3.0_dp * ( z1/ome + lamda + z2to3_2*log_ome )
     hs_term3 = 3.0_dp*z1*z2/ome2 + (z0-z2*z2to3_2)/ome + lamda*(z2/ome - z2to3) - 2.0_dp*z2to3_3*log_ome  &
                + z2*lamda/one_minus_phi * phi_1
  end if

  !-----------------------------------------------------------------------------
  ! calculate derivative d(w)/d(rho_k) for each substance k
  !-----------------------------------------------------------------------------
  do  k = 1, ncomp

     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k) : hard sphere contribution
     !--------------------------------------------------------------------------
     whs_rk = ( -z0_rk(k)*log_ome + z1_rk(k)*hs_term1 + z2_rk(k)*hs_term2 + z3_rk(k)*hs_term3 ) / PI_6

     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k) : chain term
     !--------------------------------------------------------------------------
     gij_rk_t1 = z3_rk(k)/ome2
     gij_rk_t2 = 3.0_dp*(z2_rk(k)+2.0_dp*z2*z3_rk(k)/ome)/ome2
     gij_rk_t3 = z2/ome3  *(4.0_dp*z2_rk(k)+6.0_dp*z2*z3_rk(k)/ome)
     do i = 1, ncomp
        gij_rk(i,i) = gij_rk_t1 + dij_ab(i,i)*( gij_rk_t2 + dij_ab(i,i) *gij_rk_t3 )
     end do

     whc_rk = 0.0_dp
     do i = 1, ncomp
        whc_rk = whc_rk + x(i)*rho * (1.0_dp-mseg(i)) / gij(i,i) * gij_rk(i,i)
     end do
     whc_rk = whc_rk + ( 1.0_dp-mseg(k)) * LOG( gij(k,k) )


     !--------------------------------------------------------------------------
     ! PC-SAFT:  d(f)/d(rho_k) : dispersion contribution
     !--------------------------------------------------------------------------

     ! --- derivative d(m_mean)/d(rho_k) ---------------------------------------
     m_rk = ( mseg(k) - m_mean ) ! / rho                 ! multiplied by rho, as m_rk=rho*d(m_mean)/d(rho_k)

     I1_rk = rho * z3_rk(k) * I1_z + m_rk * I1_m         ! multiplied by rho, as I1_rk=rho*d(I1)/d(rho_k)
     I2_rk = rho * z3_rk(k) * I2_z + m_rk * I2_m         ! multiplied by rho, as I2_rk=rho*d(I2)/d(rho_k)

     abbrev = - 2.0_dp *PI *mseg(k)*rho
     ord1_rk = 2.0_dp *abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k) )
     ord2_rk = abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k)**2 )

     if ( lij_correction ) ord1_rk = ord1_rk + rho * ord1_lij(k)

     c1_rk = rho * z3_rk(k) * c2_con - m_rk * c_factor   ! multiplied by rho, as c1_rk=rho*d(c1)/d(rho_k)

     wdsp_rk = order1*rho*I1_rk + ord1_rk*I1  &
          +  c1_con*m_mean * ( order2*rho*I2_rk + ord2_rk*I2 )  &
          +  ( c1_con*m_rk + c1_rk*m_mean ) * order2*rho*I2


     !--------------------------------------------------------------------------
     ! TPT-1-association according to Chapman et al.
     !--------------------------------------------------------------------------
     whb_rk = 0.0_dp

     if ( assoc ) then

        do i = 1, ncomp
           do j = (i+1), ncomp
              gij_rk(i,j) = gij_rk_t1 + dij_ab(i,j)*( gij_rk_t2 + dij_ab(i,j) *gij_rk_t3 )
              gij_rk(j,i) = gij_rk(i,j)
           end do
        end do

        whb_rk = sum ( nhb_no( k, 1:nhb_typ(k) ) * LOG( mx( k, 1:nhb_typ(k) ) ) )
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              factor_i = rho * x(i) * mx(i,ii) * nhb_no(i,ii)
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_ij = factor_i * rho * x(j) * mx(j,jj) * nhb_no(j,jj)
                    whb_rk = whb_rk - factor_ij / 2.0_dp * gij_rk(i,j) *ass_d(i,j,ii,jj)
                 end do
              end do
           end do
        end do

     end if


     !--------------------------------------------------------------------------
     ! polar terms
     !--------------------------------------------------------------------------

     ! ------ dipole-dipole term -----------------------------------------------
     wdd_rk = 0.0_dp
     if ( dipole ) then

        if ( abs( wdd2 ) > 1.E-50_dp ) then

           wdd2_rk = wdd2_z_term * z3_rk(k)
           wdd3_rk = wdd3_z_term * z3_rk(k)
           if ( abs( my_factor(k) ) > machine_eps ) then
              do i = 1, ncomp
                 if ( abs( my_factor(i)*rhoi(i) ) > machine_eps ) then
                    wdd2_rk = wdd2_rk + rhoi(i) * 2.0_dp * pi_dd(i,k) * Idd2(i,k)
                    do j = 1, ncomp
                       wdd3_rk = wdd3_rk + 3.0_dp *rhoi(i) *rhoi(j) *psi_dd(i,j,k) *Idd3(i,j,k)
                    end do
                 end if
              end do
           end if

           wdd_rk = wdd2 *( wdd2*wdd2_rk - 2.0_dp*wdd3*wdd2_rk+wdd2*wdd3_rk ) / (wdd2-wdd3)**2
           ! wdd_rk = wdd2 * ( wdd2_rk + ( wdd2*wdd3_rk - wdd3*wdd2_rk ) / (wdd2-wdd3) ) / (wdd2-wdd3)

        end if

     end if

     ! ------ quadrupole-quadrupole term ---------------------------------------
     wqq_rk = 0.0_dp
     if ( qudpole ) then

        if ( abs( wqq2 ) > 1.E-50_dp ) then

           wqq2_rk = wqq2_z_term * z3_rk(k)
           wqq3_rk = wqq3_z_term * z3_rk(k)
           if ( abs( Q_factor(k) ) > machine_eps ) then
              do i = 1, ncomp
                 if ( abs( Q_factor(i)*rhoi(i) ) > machine_eps ) then
                    wqq2_rk = wqq2_rk + rhoi(i) * 2.0_dp * pi_qq(i,k) * Iqq2(i,k)
                    do j = 1, ncomp
                       wqq3_rk = wqq3_rk + 3.0_dp *rhoi(i) *rhoi(j) *psi_qq(i,j,k) *Iqq3(i,j,k)
                    end do
                 end if
              end do
           end if

           wqq_rk = wqq2 *( wqq2*wqq2_rk - 2.0_dp*wqq3*wqq2_rk+wqq2*wqq3_rk ) / (wqq2-wqq3)**2
           ! wqq_rk = wqq2 * ( wqq2_rk + ( wqq2*wqq3_rk - wqq3*wqq2_rk ) / (wqq2-wqq3) ) / (wqq2-wqq3)

        end if

     end if

     ! ------ dipole-quadrupole term -------------------------------------------
     wdq_rk = 0.0_dp
     if ( dipole_quad ) then

        if ( abs( wdq2 ) > 1.E-50_dp ) then

           wdq2_rk = wdq2_z_term * z3_rk(k)
           wdq3_rk = wdq3_z_term * z3_rk(k)
           do i = 1, ncomp
              wdq2_rk = wdq2_rk + rhoi(i) * ( pi_dq(i,k) + pi_dq(k,i) ) * Idq2(i,k)
              do j = 1, ncomp
                 wdq3_rk = wdq3_rk + rhoi(i) * rhoi(j)  &
                      * ( psi_dq(i,j,k)+psi_dq(i,k,j)+psi_dq(k,i,j) ) * Idq3(i,j,k)
              end do
           end do

           wdq_rk = wdq2 *( wdq2*wdq2_rk - 2.0_dp*wdq3*wdq2_rk+wdq2*wdq3_rk ) / (wdq2-wdq3)**2
           ! wdq_rk = wdq2 * ( wdq2_rk + ( wdq2*wdq3_rk - wdq3*wdq2_rk ) / (wdq2-wdq3) ) / (wdq2-wdq3)

        end if

     end if


     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k) : summation of all contributions
     !--------------------------------------------------------------------------
     myres(k) = whs_rk + whc_rk + wdsp_rk + whb_rk + wdd_rk + wqq_rk + wdq_rk !JRE: pause here for sample calcs.

  end do

end subroutine F_density_rhok



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_density_rhok
!
! calculate the second derivatives of w = A/(VkT) = a*rho to
! mole fraction x and to T: ( d2(w)/d(rho_k)d(T) )
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_density_rhok_T ( w_rk_t )

  use PARAMETERS, only: PI
  USE EOS_CONSTANTS, only: ddp4, qqp4, dqp4
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k, lij_correction

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(out)         :: w_rk_t
  !-----------------------------------------------------------------------------
  integer                                         :: i, j, ii, jj, k, kk
  real(dp)                                        :: rhoPI_6
  real(dp)                                        :: z1_t, z2_t, z3_t
  real(dp)                                        :: abbrev
  real(dp)                                        :: zfr1t, zfr2t, zfr3t, zfr_omet
  real(dp)                                        :: zfr1k, zfr2k, zfr3k, zfr_omek
  real(dp)                                        :: z1_rk_t, z2_rk_t, z3_rk_t
  real(dp)                                        :: zfr1kt, zfr2kt, zfr3kt, zfr_ome_kt
  real(dp), dimension(ncomp)                      :: dhs_t

  real(dp)                                        :: m_rk
  real(dp)                                        :: whs_rk_T_t1, whs_rk_T_t2, whs_rk_T_t3
  real(dp)                                        :: whs_rk_t

  real(dp)                                        :: whc_rk, whc_rk_t
  real(dp)                                        :: gij_t2_rho, gij_t3_rho
  real(dp)                                        :: gij_t_t1, gij_t_t2, gij_t_t3
  real(dp)                                        :: gij_rk_t1, gij_rk_t2, gij_rk_t3
  real(dp)                                        :: gij_rk_t_t1, gij_rk_t_t2, gij_rk_t_t3
  real(dp), dimension(ncomp,ncomp)                :: dijdt
  real(dp), dimension(ncomp,ncomp)                :: gij_t
  real(dp), dimension(ncomp,ncomp)                :: gij_rk
  real(dp), dimension(ncomp,ncomp)                :: gij_rk_t

  real(dp)                                        :: wdsp_rk_t
  real(dp)                                        :: I1_m, I2_m, I1_z_m, I2_z_m
  real(dp)                                        :: I1_t, I2_t
  real(dp)                                        :: I1_rk, I2_rk
  real(dp)                                        :: I1_rk_t, I2_rk_t
  real(dp)                                        :: ord1_rk, ord2_rk
  real(dp)                                        :: c1_rk, c1_t, c1_rk_t
  real(dp)                                        :: c_factor, c_factor_t

  integer                                         :: ass_cnt, max_eval
  real(dp), dimension(ncomp,ncomp,nsite,nsite)    :: delta, delta_t, delta_rk, delta_rk_t
  real(dp), dimension(ncomp,nsite)                :: mx_t, mx_rk
  real(dp), dimension(ncomp,nsite)                :: mx_itr, mx_itr_t, mx_itr_rk
  real(dp)                                        :: tol
  real(dp)                                        :: sum0, sum_t, sum_rk, err_sum
  real(dp)                                        :: ass_s2
  real(dp)                                        :: fhb, fhb_t
  real(dp)                                        :: whb_rk, fhb_rk
  real(dp)                                        :: whb_rk_t, fhb_rk_t
  real(dp)                                        :: factor_i, factor_j

  real(dp)                                        :: wdd_rk, wdd_rk_t
  real(dp)                                        :: wdd2_t, wdd3_t
  real(dp)                                        :: wdd2_rk, wdd3_rk
  real(dp)                                        :: wdd2_rk_t, wdd3_rk_t
  real(dp)                                        :: p_factor
  real(dp)                                        :: wdd2B_term, wdd2B_z_term
  real(dp)                                        :: Idd2B_z_ij
  real(dp), dimension(ncomp,ncomp)                :: Idd2B
  real(dp)                                        :: diff_f, diff_factor

  real(dp)                                        :: wqq_rk, wqq_rk_t
  real(dp)                                        :: wqq2_t, wqq3_t
  real(dp)                                        :: wqq2_rk, wqq3_rk
  real(dp)                                        :: wqq2_rk_t, wqq3_rk_t
  real(dp)                                        :: wqq2B_term, wqq2B_z_term
  real(dp)                                        :: Iqq2B_z_ij
  real(dp), dimension(ncomp,ncomp)                :: Iqq2B

  real(dp)                                        :: wdq_rk, wdq_rk_t
  real(dp)                                        :: wdq2_t, wdq3_t
  real(dp)                                        :: wdq2_rk, wdq3_rk
  real(dp)                                        :: wdq2_rk_t, wdq3_rk_t
  real(dp)                                        :: wdq2B_term, wdq2B_z_term
  real(dp)                                        :: Idq2B_z_ij
  real(dp), dimension(ncomp,ncomp)                :: Idq2B
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! derivative of some auxilliary properties to temperature
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     dhs_t(i) = sigma(i) *(-3.0_dp*epsilon_k(i)/t/t)*0.12_dp*EXP(-3.0_dp*epsilon_k(i)/t)
  end do

  rhoPI_6 = PI_6 * rho

  z1_t = 0.0_dp
  z2_t = 0.0_dp
  z3_t = 0.0_dp
  do i = 1, ncomp
     abbrev = x(i) * mseg(i) * dhs_t(i)
     z1_t = z1_t + abbrev
     z2_t = z2_t + abbrev * 2.0_dp * dhs(i)
     z3_t = z3_t + abbrev * 3.0_dp * dhs(i) * dhs(i)
  end do
  z1_t  = rhoPI_6 * z1_t
  z2_t  = rhoPI_6 * z2_t
  z3_t  = rhoPI_6 * z3_t

  zfr1t = z1_t / z1
  zfr2t = z2_t / z2
  zfr3t = z3_t / z3
  zfr_omet = z3_t / ome

  !-----------------------------------------------------------------------------
  ! quantities independent of k. Index k cycles for the derivative  d(w)/d(rho_k)
  !-----------------------------------------------------------------------------
  !gij_t1 = 1.0_dp/ome
  !gij_t2 = 3.0_dp*z2t/ome2
  !gij_t3 = 2.0_dp*z2t*z2/ome3
  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3
  gij_t_t1 = z3_t/ome2
  gij_t_t2 = 3.0_dp*(z2_t+2.0_dp*z2*z3_t/ome)/ome2
  gij_t_t3 = z2/ome3  *(4.0_dp*z2_t+6.0_dp*z2*z3_t/ome)
  do i = 1, ncomp
     do j = i, ncomp
        dijdt(i,j) = (dhs_t(i)*dhs(j) + dhs(i)*dhs_t(j)) / (dhs(i)+dhs(j))  &
                  - dhs(i)*dhs(j)/(dhs(i)+dhs(j))**2  *(dhs_t(i)+dhs_t(j))
        !gij(i,j) = gij_t1 + dij_ab(i,j)*( gij_t2_rho + dij_ab(i,j) *gij_t3_rho )
        gij_t(i,j) = gij_t_t1 + dij_ab(i,j)*( gij_t_t2 + dij_ab(i,j) *gij_t_t3 ) &
                  + dijdt(i,j) * ( gij_t2_rho + 2.0_dp*dij_ab(i,j) *gij_t3_rho )
        !gij(j,i) = gij(i,j)
        gij_t(j,i) = gij_t(i,j)
        dijdt(j,i) = dijdt(i,j)
     end do
  end do

  c1_t = c2_con*z3_t
  c_factor = (8.0_dp*z3-2.0_dp*z32)/ome4  &
       - (-2.0_dp*z3**4 +12.0_dp*z33 -27.0_dp*z32+20.0_dp*z3) *ome_t_2
  c_factor = c1_con * c1_con * c_factor
  c_factor_t = (8.0_dp-4.0_dp*z3)/ome4 + 8.0_dp*(4.0_dp*z3-z32)/ome5  &
       - (-8.0_dp*z33 +36.0_dp*z32 -54.0_dp*z3+20.0_dp) *ome_t_2  &
       + (-2.0_dp*z3**4 +12.0_dp*z33 -27.0_dp*z32+20.0_dp*z3) *ome_t_2*ome_t*4.0_dp*(z3-1.5_dp)
  c_factor_t = c_factor_t * z3_t * c1_con * c1_con


  I1_m = sum(  ap_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I1 to mean segment number m_mean
  I2_m = sum(  bp_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I2 to mean segment number m_mean
  I1_z_m = sum( ap_rk_aux(1:6) * z3_m_m1(1:6) ) ! derivative d(I1)/d(z3) to m_mean
  I2_z_m = sum( bp_rk_aux(1:6) * z3_m_m1(1:6) ) ! derivative d(I2)/d(z3) to m_mean

  ! ------ dipole-dipole term --------------------------------------------------
  if ( dipole ) then

     Idd2B(:,:) = 0.0_dp
     wdd2B_term = 0.0_dp
     wdd2B_z_term = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j) ) > machine_eps ) then
                 Idd2B(i,j) = eps_ij(i,j) * sum( ddp4(i,j,0:4) * z3_m_m0(0:4) )
                 Idd2B_z_ij = eps_ij(i,j) * sum( ddp4(i,j,1:4) * z3_m_m1(1:4) )
                 p_factor = rhoi(i) * rhoi(j) * pi_dd(i,j)
                 wdd2B_term = wdd2B_term + p_factor * Idd2B(i,j)
                 wdd2B_z_term = wdd2B_z_term + p_factor * Idd2B_z_ij
              end if
           end do
        end if
     end do
     wdd2_t = wdd2_z_term * z3_t - 2.0_dp * wdd2 / t - wdd2B_term / t
     wdd3_t = wdd3_z_term * z3_t - 3.0_dp * wdd3 / t

  end if

  ! ------ quadrupole-quadrupole term ------------------------------------------

  if ( qudpole ) then

     Iqq2B(:,:) = 0.0_dp
     wqq2B_term = 0.0_dp
     wqq2B_z_term = 0.0_dp
     do i = 1, ncomp
        if ( abs( Q_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j) ) > machine_eps ) then
                 Iqq2B(i,j) = eps_ij(i,j) * sum( qqp4(i,j,0:4) * z3_m_m0(0:4) )
                 Iqq2B_z_ij = eps_ij(i,j) * sum( qqp4(i,j,1:4) * z3_m_m1(1:4) )
                 p_factor = rhoi(i) * rhoi(j) * pi_qq(i,j)
                 wqq2B_term = wqq2B_term + p_factor * Iqq2B(i,j)
                 wqq2B_z_term = wqq2B_z_term + p_factor * Iqq2B_z_ij
              end if
           end do
        end if
     end do
     wqq2_t = wqq2_z_term * z3_t - 2.0_dp * wqq2 / t - wqq2B_term / t
     wqq3_t = wqq3_z_term * z3_t - 3.0_dp * wqq3 / t

  end if

  
  ! ------ dipole-quadrupole term ----------------------------------------------

  if ( dipole_quad ) then

     Idq2B(:,:) = 0.0_dp
     wdq2B_term = 0.0_dp
     wdq2B_z_term = 0.0_dp
     do i = 1, ncomp
        do j = 1, ncomp
           Idq2B(i,j) = eps_ij(i,j) * sum( dqp4(i,j,0:4) * z3_m_m0(0:4) )
           Idq2B_z_ij = eps_ij(i,j) * sum( dqp4(i,j,1:4) * z3_m_m1(1:4) )
           p_factor = rhoi(i) * rhoi(j) * pi_dq(i,j)
           wdq2B_term = wdq2B_term + p_factor * Idq2B(i,j)
           wdq2B_z_term = wdq2B_z_term + p_factor * Idq2B_z_ij
        end do
     end do
     wdq2_t = wdq2_z_term * z3_t - 2.0_dp * wdq2 / t - wdq2B_term / t
     wdq3_t = wdq3_z_term * z3_t - 3.0_dp * wdq3 / t

  end if

  if ( assoc ) then

     tol = hb_tol
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-3_dp

     max_eval = 200
     mx_t(:,:) = 0.0_dp

  end if

     
  !-----------------------------------------------------------------------------
  ! calculate derivative d(w)/d(rho_k) for each substance k
  !-----------------------------------------------------------------------------
  do  k = 1, ncomp

     zfr1k = z1_rk(k) / z1
     zfr2k = z2_rk(k) / z2
     zfr3k = z3_rk(k) / z3
     zfr_omek = z3_rk(k) / ome

     z1_rk_t  = PI_6 * mseg(k) * dhs_t(k)
     z2_rk_t  = 2.0_dp * dhs(k) * z1_rk_t
     z3_rk_t  = 1.5_dp * dhs(k) * z2_rk_t

     ! --- derivative d(m_mean)/d(rho_k) ---------------------------------------
     m_rk = ( mseg(k) - m_mean ) / rho

     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : hard sphere contribution
     !--------------------------------------------------------------------------
     !whs = (  3.0_dp*z1*z2/ome + z2**3 /z3/ome2 + (z2**3 /z3/z3-z0)*LOG(ome)  ) / PI_6
     !whs_rk = (  3.0_dp*(z1_rk(k)*z2+z1*z2_rk(k))/ome + 3.0_dp*z1*z2*z3_rk(k)/ome2  &
     !        + 3.0_dp*z2*z2*z2_rk(k)/z3/ome2 + z2**3 *z3_rk(k)*(3.0_dp*z3-1.0_dp)/z3/z3/ome3   &
     !        + ((3.0_dp*z2*z2*z2_rk(k)*z3-2.0_dp*z2**3 *z3_rk(k))/z3**3 -z0_rk(k)) *LOG(ome)  &
     !        + (z0-z2**3 /z3/z3)*z3_rk(k)/ome  ) / PI_6
     !whs_rk_t = (  3.0_dp*(z1_rk(k)_t*z2+z1_rk(k)*z2_t+z1_t*z2_rk(k)+z1*z2_rk_t)/ome  &
     !        + 3.0_dp*(z1_rk(k)*z2+z1*z2_rk(k))/ome2*zfr_ome
     !        + 3.0_dp*z1*z2*z3_rk(k)/ome2*(zfr1+zfr2+z3_rk_t/z3_rk(k)+2.0_dp*zfr_ome)  &
     !        + 3.0_dp*z2*z2*z2_rk(k)/z3/ome2*(2.0_dp*zfr2+z2_rk_t/z2_rk(k)+2.0_dp*zfr_ome)
     !        + z2**3 *z3_rk(k)*(3.0_dp*z3-1.0_dp)/z3/z3/ome3 &
     !            *(3.0_dp*zfr3+z3_rk_t/z3_rk(k)+z3_t/(z3-1.0_dp/3.0_dp)-2.0_dp*z3_t/z3+3.0_dp*zfr_ome)   &
     !   + ((3.0_dp*z2*z2*z2_rk(k)*z3-2.0_dp*z2**3 *z3_rk)/z3**3 - z0_rk(k)) *LOG(ome)  &
     !        + (z0-z2**3 /z3/z3)*z3_rk(k)/ome  ) / PI_6

     zfr1kt = z1_rk_t / z1 - zfr1k*zfr1t
     zfr2kt = z2_rk_t / z2 - zfr2k*zfr2t
     zfr3kt = z3_rk_t / z3 - zfr3k*zfr3t
     zfr_ome_kt = z3_rk_t / ome + zfr_omek * zfr_omet
     whs_rk_T_t1 = 3.0_dp*z1*z2/ome * ( (zfr1k +zfr2k +zfr_omek) *(zfr1t +zfr2t +zfr_omet)  &
          + zfr1kt + zfr2kt + zfr_ome_kt )
     whs_rk_T_t2 = z23 /z3/ome2 * ( (3.0_dp*zfr2k -zfr3k +2.0_dp*zfr_omek)  &
          * (3.0_dp*zfr2t -zfr3t +2.0_dp*zfr_omet) + 3.0_dp*zfr2kt - zfr3kt + 2.0_dp*zfr_ome_kt )
     whs_rk_T_t3 = z23/z32 *( ( 3.0_dp*zfr2k-2.0_dp*zfr3k) * (3.0_dp*zfr2t-2.0_dp*zfr3t)  &
          +  3.0_dp*zfr2kt-2.0_dp*zfr3kt ) * LOG(ome)  &
          - (z23/z32 *(3.0_dp*zfr2k-2.0_dp*zfr3k) -z0_rk(k))*zfr_omet  &
          - (z23/z32 *(3.0_dp*zfr2t-2.0_dp*zfr3t) )*zfr_omek - (z23/z32-z0)*zfr_ome_kt
     whs_rk_T = ( whs_rk_T_t1 + whs_rk_T_t2 + whs_rk_T_t3 )  / PI_6

             
     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : chain term
     !--------------------------------------------------------------------------
     gij_rk_t1 = z3_rk(k)/ome2
     gij_rk_t2 = gij_t2*rho*(zfr2k+2.0_dp*zfr_omek)
     gij_rk_t3 = gij_t3*rho*(2.0_dp*zfr2k+3.0_dp*zfr_omek)
     gij_rk_t_t1 = zfr_ome_kt / ome + z3_rk(k)*z3_t/ome3
     gij_rk_t_t2 = gij_t2*rho*( (zfr2k+2.0_dp*zfr_omek)*(zfr2t+2.0_dp*zfr_omet)  &
                   + zfr2kt+2.0_dp*zfr_ome_kt )
     gij_rk_t_t3 = gij_t3*rho*( (2.0_dp*zfr2k+3.0_dp*zfr_omek)*(2.0_dp*zfr2t+3.0_dp*zfr_omet)  &
                   + 2.0_dp*zfr2kt+3.0_dp*zfr_ome_kt )
     do i = 1, ncomp
        gij_rk(i,i) = gij_rk_t1 + dij_ab(i,i)*( gij_rk_t2 + dij_ab(i,i) *gij_rk_t3 )
        gij_rk_t(i,i) = gij_rk_t_t1 + dij_ab(i,i)*( gij_rk_t_t2 + dij_ab(i,i) *gij_rk_t_t3 ) &
                                    + dijdt(i,i)* ( gij_rk_t2 + 2.0_dp*dij_ab(i,i) *gij_rk_t3 )
     end do

     whc_rk = 0.0_dp
     whc_rk_t = 0.0_dp
     do i = 1, ncomp
        whc_rk = whc_rk + x(i)*rho * (1.0_dp-mseg(i)) * gij_rk(i,i) / gij(i,i)
        whc_rk_t = whc_rk_t + x(i)*rho * (1.0_dp-mseg(i))  &
                   * ( gij_rk_t(i,i) - gij_rk(i,i)*gij_t(i,i)/gij(i,i) ) / gij(i,i)
     end do
     whc_rk = whc_rk + ( 1.0_dp-mseg(k)) * LOG( gij(k,k) )
     whc_rk_t = whc_rk_t + ( 1.0_dp-mseg(k)) / gij(k,k) * gij_t(k,k)


     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : dispersion contribution
     !--------------------------------------------------------------------------

     I1_t = I1_z * z3_t
     I2_t = I2_z * z3_t
     I1_rk = z3_rk(k) * I1_z + m_rk * I1_m
     I2_rk = z3_rk(k) * I2_z + m_rk * I2_m
     I1_rk_t = I1_z * z3_rk_t + ( I1_z2 * z3_rk(k) + m_rk * I1_z_m ) * z3_t
     I2_rk_t = I2_z * z3_rk_t + ( I2_z2 * z3_rk(k) + m_rk * I2_z_m ) * z3_t

     abbrev = - 2.0_dp *PI *mseg(k)*rho
     ord1_rk = 2.0_dp *abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k) )
     ord2_rk = abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k)**2 )

     if ( lij_correction ) ord1_rk = ord1_rk + rho * ord1_lij(k)

     c1_rk = z3_rk(k) * c2_con - m_rk * c_factor
     c1_rk_t = z3_rk(k) * z3_t * c3_con + z3_rk_t * c2_con  &
          - m_rk * ( 2.0_dp * c1_t /c1_con * c_factor + c_factor_t )

     !wdsp_rk = order1*rho*rho*I1_rk + ord1_rk*I1  &
     !     + c1_con*m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
     !     + ( c1_con*m_rk + c1_rk*m_mean ) * order2*rho*rho*I2
     wdsp_rk_t = order1*rho*rho*(I1_rk_t-I1_rk/t) + ord1_rk*(I1_t-I1/t)  &
          + c1_con*m_mean * ( order2*rho*rho*(I2_rk_t-2.0_dp*I2_rk/t)  &
                             + ord2_rk*(I2_t-2.0_dp*I2/t) )  &
          + c1_t *m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
          + order2*rho*rho*( ( c1_t*m_rk + c1_rk_t*m_mean ) * I2  &
                                   +( c1_con*m_rk + c1_rk*m_mean ) * (I2_t-2.0_dp*I2/t) )


     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : TPT-1-association according to Chapman et al.
     !--------------------------------------------------------------------------

     fhb = 0.0_dp
     fhb_t = 0.0_dp
     fhb_rk = 0.0_dp
     fhb_rk_t = 0.0_dp

     whb_rk = 0.0_dp
     whb_rk_t = 0.0_dp
     
     if ( assoc ) then

        ass_s2  = 0.0_dp
        do kk = 1, nhb_typ(k)
           ass_s2  = ass_s2  + nhb_no(k,kk) * LOG(mx(k,kk))
        end do

        do i = 1, ncomp
           do j = (i+1), ncomp
              gij_rk(i,j) = gij_rk_t1 + dij_ab(i,j)*( gij_rk_t2 + dij_ab(i,j) *gij_rk_t3 )
              gij_rk_t(i,j) = gij_rk_t_t1 + dij_ab(i,j)*( gij_rk_t_t2 + dij_ab(i,j) *gij_rk_t_t3 ) &
                   + dijdt(i,j)* ( gij_rk_t2 + 2.0_dp*dij_ab(i,j) *gij_rk_t3 )
              gij_rk(j,i) = gij_rk(i,j)
              gij_rk_t(j,i) = gij_rk_t(i,j)
           end do
        end do

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    delta(i,j,ii,jj) = gij(i,j)*ass_d(i,j,ii,jj)
                    delta_t(i,j,ii,jj) = gij_t(i,j)*ass_d(i,j,ii,jj) + gij(i,j)*ass_d_dt(i,j,ii,jj)
                    delta_rk(i,j,ii,jj) = gij_rk(i,j) * ass_d(i,j,ii,jj)
                    delta_rk_t(i,j,ii,jj) = gij_rk_t(i,j) * ass_d(i,j,ii,jj)  &
                         + gij_rk(i,j) * ass_d_dt(i,j,ii,jj)
                 end do
              end do
           end do
        end do


        !--- initialize mxdt(i,j) -------------------------------------------------
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              mx_rk(i,ii) = 0.0_dp
           end do
        end do

        !--- iterate over all components and all sites ----------------------------
        err_sum = 10.0_dp * tol
        ass_cnt = 0

        do while ( err_sum > tol .AND. ass_cnt <= max_eval )

           do i = 1, ncomp
              do ii = 1, nhb_typ(i)
                 sum0 = 0.0_dp
                 sum_t = 0.0_dp
                 sum_rk = 0.0_dp
                 do j = 1, ncomp
                    do jj = 1, nhb_typ(j)
                       sum0  = sum0  + rhoi(j)*nhb_no(j,jj)*  mx(j,jj) *delta(i,j,ii,jj)
                       sum_t = sum_t + rhoi(j)*nhb_no(j,jj)*( mx(j,jj) *delta_t(i,j,ii,jj)  &
                            + mx_t(j,jj)*delta(i,j,ii,jj) )
                       sum_rk = sum_rk + rhoi(j)*nhb_no(j,jj)*( mx(j,jj)*delta_rk(i,j,ii,jj)  &
                            + mx_rk(j,jj)*delta(i,j,ii,jj) )
                    end do
                 end do
                 do kk = 1, nhb_typ(k)
                    sum_rk = sum_rk + nhb_no(k,kk) * mx(k,kk)*delta(i,k,ii,kk)
                 end do
                 mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + sum0 )
                 mx_itr_t(i,ii) = - mx_itr(i,ii)**2 * sum_t
                 mx_itr_rk(i,ii) = - mx_itr(i,ii)**2 * sum_rk
              end do
           end do

           err_sum = 0.0_dp
           do i = 1, ncomp
              do ii = 1, nhb_typ(i)
                 err_sum = err_sum + ABS(mx_itr(i,ii) - mx(i,ii))  &
                      + ABS(mx_itr_t(i,ii) - mx_t(i,ii)) + ABS(mx_itr_rk(i,ii) - mx_rk(i,ii))
                 mx(i,ii) = mx_itr(i,ii) * hb_damp + mx(i,ii) * hb_damp_inv
                 mx_t(i,ii) = mx_itr_t(i,ii) * hb_damp + mx_t(i,ii) * hb_damp_inv
                 mx_rk(i,ii) = mx_itr_rk(i,ii) * hb_damp + mx_rk(i,ii) * hb_damp_inv
              end do
           end do

           ass_cnt = ass_cnt + 1

           if ( ass_cnt == max_eval ) then
              write (6,*) 'F_density_rhok_T: max_eval violated err_sum = ',err_sum,tol
              exit
           end if

        end do

        whb_rk = sum( nhb_no( k, 1:nhb_typ(k) ) * LOG( mx( k, 1:nhb_typ(k) ) ) )
        whb_rk_t = sum( nhb_no( k,1:nhb_typ(k) ) *mx_t( k,1:nhb_typ(k) ) /mx( k,1:nhb_typ(k) ) )
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              fhb = fhb + x(i)* nhb_no(i,ii)* ( 0.5_dp * ( 1.0_dp - mx(i,ii) ) + LOG(mx(i,ii)) )
              fhb_t = fhb_t  +x(i)*nhb_no(i,ii) * mx_t(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)
              factor_i = 0.5_dp * rhoi(i) * nhb_no(i,ii) * mx(i,ii)
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_j = rhoi(j) * nhb_no(j,jj) * mx(j,jj)
                    whb_rk = whb_rk - factor_i * factor_j * delta_rk(i,j,ii,jj)
                    whb_rk_t = whb_rk_t - factor_i*factor_j * delta_rk_t(i,j,ii,jj) &
                         -  ( 0.5_dp*rhoi(i) *nhb_no(i,ii) *mx_t(i,ii) *factor_j  &
                         + factor_i *rhoi(j) *nhb_no(j,jj) *mx_t(j,jj) ) *delta_rk(i,j,ii,jj)
                 end do
              end do
           end do
        end do
        fhb_rk = ( whb_rk - fhb ) / rho
        fhb_rk_t = ( whb_rk_t - fhb_t ) / rho

     end if

     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : dipole-dipole contribution
     !--------------------------------------------------------------------------

     wdd_rk_t = 0.0_dp

     if ( dipole ) then

        if ( abs( wdd2 ) > 1.E-50_dp ) then

           wdd2_rk = wdd2_z_term * z3_rk(k)
           wdd3_rk = wdd3_z_term * z3_rk(k)
           wdd2_rk_t = ( wdd2_z2_term * z3_t - wdd2B_z_term / t ) * z3_rk(k) + wdd2_z_term * z3_rk_t
           wdd3_rk_t = wdd3_z2_term * z3_t * z3_rk(k) + wdd3_z_term * z3_rk_t
           if ( abs( my_factor(k) ) > machine_eps ) then
              do i = 1, ncomp
                 if ( abs( my_factor(i)*rhoi(i) ) > machine_eps ) then
                    p_factor = 2.0_dp * rhoi(i)* pi_dd(i,k)
                    wdd2_rk = wdd2_rk + p_factor * Idd2(i,k)
                    wdd2_rk_t = wdd2_rk_t + p_factor * ( Idd2_z(i,k) * z3_t - Idd2B(i,k)/t )
                    do j = 1, ncomp
                       p_factor = 3.0_dp * rhoi(i) * rhoi(j) * psi_dd(i,j,k)
                       wdd3_rk = wdd3_rk + p_factor * Idd3(i,j,k)
                       wdd3_rk_t = wdd3_rk_t + p_factor * Idd3_z(i,j,k) * z3_t
                    end do
                 end if
              end do
           end if
           wdd2_rk_t = wdd2_rk_t - 2.0_dp / t * wdd2_rk
           wdd3_rk_t = wdd3_rk_t - 3.0_dp / t * wdd3_rk
           diff_f = wdd2 - wdd3
           diff_factor = wdd2 / diff_f / diff_f

           wdd_rk = wdd2 * ( wdd2_rk - 2.0_dp*wdd3/wdd2*wdd2_rk + wdd3_rk ) * diff_factor

           wdd_rk_t = ( 2.0_dp*wdd2_t*( wdd2_rk + wdd3_rk ) + wdd2*( wdd2_rk_t + wdd3_rk_t )  &
                - 2.0_dp*( wdd2_rk*wdd2_t*wdd3/wdd2 + wdd3_t*wdd2_rk + wdd3*wdd2_rk_t ) ) * diff_factor  &
                + 2.0_dp * wdd_rk * ( wdd3_t - wdd2_t ) / diff_f

        end if


     end if

     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : quadrupole-quadrupole contribution
     !--------------------------------------------------------------------------

     wqq_rk_t = 0.0_dp

     if ( qudpole ) then

        if ( abs( wqq2 ) > 1.E-50_dp ) then

           wqq2_rk = wqq2_z_term * z3_rk(k)
           wqq3_rk = wqq3_z_term * z3_rk(k)
           wqq2_rk_t = ( wqq2_z2_term * z3_t - wqq2B_z_term / t ) * z3_rk(k) + wqq2_z_term * z3_rk_t
           wqq3_rk_t = wqq3_z2_term * z3_t * z3_rk(k) + wqq3_z_term * z3_rk_t
           if ( abs( Q_factor(k) ) > machine_eps ) then
              do i = 1, ncomp
                 if ( abs( Q_factor(i)*rhoi(i) ) > machine_eps ) then
                    p_factor = 2.0_dp * rhoi(i)* pi_qq(i,k)
                    wqq2_rk = wqq2_rk + p_factor * Iqq2(i,k)
                    wqq2_rk_t = wqq2_rk_t + p_factor * ( Iqq2_z(i,k) * z3_t - Iqq2B(i,k)/t )
                    do j = 1, ncomp
                       p_factor = 3.0_dp * rhoi(i) * rhoi(j) * psi_qq(i,j,k)
                       wqq3_rk = wqq3_rk + p_factor * Iqq3(i,j,k)
                       wqq3_rk_t = wqq3_rk_t + p_factor * Iqq3_z(i,j,k) * z3_t
                    end do
                 end if
              end do
           end if
           wqq2_rk_t = wqq2_rk_t - 2.0_dp / t * wqq2_rk
           wqq3_rk_t = wqq3_rk_t - 3.0_dp / t * wqq3_rk
           diff_f = wqq2 - wqq3
           diff_factor = wqq2 / diff_f / diff_f

           wqq_rk = wqq2 * ( wqq2_rk - 2.0_dp*wqq3/wqq2*wqq2_rk + wqq3_rk ) * diff_factor

           wqq_rk_t = ( 2.0_dp*wqq2_t*( wqq2_rk + wqq3_rk ) + wqq2*( wqq2_rk_t + wqq3_rk_t )  &
                - 2.0_dp*( wqq2_rk*wqq2_t*wqq3/wqq2 + wqq3_t*wqq2_rk + wqq3*wqq2_rk_t ) ) * diff_factor  &
                + 2.0_dp * wqq_rk * ( wqq3_t - wqq2_t ) / diff_f

        end if

     end if


     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : quadrupole-quadrupole contribution
     !--------------------------------------------------------------------------

     wdq_rk_t = 0.0_dp

     if ( dipole_quad ) then

        if ( abs( wdq2 ) > 1.E-50_dp ) then

           wdq2_rk = wdq2_z_term * z3_rk(k)
           wdq3_rk = wdq3_z_term * z3_rk(k)
           wdq2_rk_t = ( wdq2_z2_term * z3_t - wdq2B_z_term / t ) * z3_rk(k) + wdq2_z_term * z3_rk_t
           wdq3_rk_t = wdq3_z2_term * z3_t * z3_rk(k) + wdq3_z_term * z3_rk_t
           if ( abs( my_fac_dq(k) ) > machine_eps .OR. abs( Q_fac_dq(k) ) > machine_eps ) then
              do i = 1, ncomp
                 if ( abs( my_fac_dq(i)*rhoi(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*rhoi(i) ) > machine_eps ) then
                    p_factor = rhoi(i)*( pi_dq(i,k) + pi_dq(k,i) )
                    wdq2_rk = wdq2_rk + p_factor * Idq2(i,k)
                    wdq2_rk_t = wdq2_rk_t + p_factor * ( Idq2_z(i,k) * z3_t - Idq2B(i,k)/t )
                    do j = 1, ncomp
                       p_factor =  rhoi(i) * rhoi(j) * ( psi_dq(i,j,k)+psi_dq(i,k,j)+psi_dq(k,i,j) )
                       wdq3_rk = wdq3_rk + p_factor * Idq3(i,j,k)
                       wdq3_rk_t = wdq3_rk_t + p_factor * Idq3_z(i,j,k) * z3_t
                    end do
                 end if
              end do
           end if
           wdq2_rk_t = wdq2_rk_t - 2.0_dp / t * wdq2_rk
           wdq3_rk_t = wdq3_rk_t - 3.0_dp / t * wdq3_rk
           diff_f = wdq2 - wdq3
           diff_factor = wdq2 / diff_f / diff_f

           wdq_rk = wdq2 * ( wdq2_rk - 2.0_dp*wdq3/wdq2*wdq2_rk + wdq3_rk ) * diff_factor

           wdq_rk_t = ( 2.0_dp*wdq2_t*( wdq2_rk + wdq3_rk ) + wdq2*( wdq2_rk_t + wdq3_rk_t )  &
                - 2.0_dp*( wdq2_rk*wdq2_t*wdq3/wdq2 + wdq3_t*wdq2_rk + wdq3*wdq2_rk_t ) ) * diff_factor  &
                + 2.0_dp * wdq_rk * ( wdq3_t - wdq2_t ) / diff_f

        end if

     end if


     !--------------------------------------------------------------------------
     ! d2(f) / d(rho_k)d(T) : summation of all contributions
     !--------------------------------------------------------------------------
     w_rk_t(k) = whs_rk_t + whc_rk_t + wdsp_rk_t + whb_rk_t + wdd_rk_t + wqq_rk_t + wdq_rk_t

  end do

end subroutine F_density_rhok_T



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_density_rhok
!
! calculate the second derivatives of w = A/(VkT) = a*rho to densities rho_k
! and to density rho ( d2(w) / d(rho_k)d(rho) )
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_density_rhok_rho ( w_rk, w_rk_r )

  use PARAMETERS, only: PI
  USE EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: lij_correction
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(out)         :: w_rk
  real(dp), dimension(ncomp), intent(out)         :: w_rk_r
  !-----------------------------------------------------------------------------
  integer                                         :: i, j, k, ii, jj, kk
  real(dp)                                        :: z1tfr, z2tfr, z3tfr

  real(dp)                                        :: m_rk, m_rk_r
  real(dp)                                        :: z1fr_rk, z2fr_rk, z3fr_rk
  
  !real(dp)                                        :: whs
  real(dp)                                        :: whs_rk, whs_rk_r
  real(dp)                                        :: whs1, whs2, whs3
  real(dp)                                        :: whs1_rk, whs2_rk, whs3_rk
  real(dp)                                        :: whs1_rk_r, whs2_rk_r, whs3_rk_r

  real(dp)                                        :: gij_r_t1, gij_r_t2, gij_r_t3
  real(dp)                                        :: gij_rk_t1, gij_rk_t2, gij_rk_t3
  real(dp)                                        :: gij_rk_r_t1, gij_rk_r_t2, gij_rk_r_t3
  real(dp), dimension(ncomp,ncomp)                :: gij_rk, gij_r, gij_rk_r
  real(dp)                                        :: whc_rk, whc_rk_r

  real(dp)                                        :: wdsp_rk, wdsp_r, wdsp_rk_r
  real(dp)                                        :: I1_rk, I2_rk, I1_rk_r, I2_rk_r
  real(dp)                                        :: abbrev
  real(dp)                                        :: ord1_rk, ord2_rk
  real(dp)                                        :: c1_rk, c1_rk_r
  real(dp)                                        :: I1_r, I2_r
  real(dp)                                        :: I1_m, I2_m
  real(dp)                                        :: I1_r_m, I2_r_m
  real(dp)                                        :: c_factor, c_factor_r

  integer                                         :: ass_cnt, max_eval
  real(dp)                                        :: whb, whb_r, whb_rk, whb_rk_r
  real(dp)                                        :: factor_i, factor_ij, factor_ij_mx
  real(dp), dimension(ncomp,nsite)                :: mx_r, mx_rk
  real(dp), dimension(ncomp,nsite)                :: mx_r_itr, mx_rk_itr
  real(dp)                                        :: sum_r, sum_rk
  real(dp)                                        :: err_sum, tol

  real(dp), dimension(ncomp)                      :: wdd_rk_r, wqq_rk_r, wdq_rk_r
  real(dp), dimension(ncomp)                      :: wdd_rk, wqq_rk, wdq_rk
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! some abbreviations
  !-----------------------------------------------------------------------------
  whs1 = 3.0_dp*z1*z2/ome
  whs2 = z23 / z3 / ome2
  whs3 = z23 / z32

  z1tfr = z1t / z1
  z2tfr = z2t / z2
  z3tfr = z3t / ome
     
  !-----------------------------------------------------------------------------
  ! quantities independent of k. Index k cycles for the derivative  d(w)/d(rho_k)
  !-----------------------------------------------------------------------------
  gij_r_t1 = z3t / ome2
  gij_r_t2 = gij_t2*rho * ( z2tfr + 2.0_dp*z3tfr )
  gij_r_t3 = gij_t3*rho * ( 2.0_dp*z2tfr + 3.0_dp*z3tfr )
  do i = 1, ncomp
     gij_r(i,i) = gij_r_t1 + dij_ab(i,i)*( gij_r_t2 + dij_ab(i,i) *gij_r_t3 )
  end do

  ! --- auxilliary for derivatives of apar, bpar to rho_k ----------------------
  ! ap_rk_aux(0:6) = ( ap(0:6,2) + (3.0_dp -4.0_dp/m_mean) *ap(0:6,3) ) / m_mean**2
  ! bp_rk_aux(0:6) = ( bp(0:6,2) + (3.0_dp -4.0_dp/m_mean) *bp(0:6,3) ) / m_mean**2

  I1_r = I1_z * z3t
  I2_r = I2_z * z3t
  I1_m = sum( ap_rk_aux(0:6) * z3_m_m0(0:6) )          ! derivative of I1 to mean segment number m_mean
  I2_m = sum( bp_rk_aux(0:6) * z3_m_m0(0:6) )          ! derivative of I2 to mean segment number m_mean
  I1_r_m = z3t * sum( ap_rk_aux(1:6) * z3_m_m1(1:6) )  ! derivative d(I1)/d(rho) to m_mean
  I2_r_m = z3t * sum( bp_rk_aux(1:6) * z3_m_m1(1:6) )  ! derivative d(I2)/d(rho) to m_mean

  c_factor = (8.0_dp*z3-2.0_dp*z32)/ome4  &
       - (-2.0_dp*z3**4 +12.0_dp*z33 -27.0_dp*z32+20.0_dp*z3) *ome_t_2
  c_factor = c_factor * c1_con * c1_con
  c_factor_r = (-4.0_dp*z32+20.0_dp*z3+8.0_dp)/ome5  &
       - (2.0_dp*z33 +12.0_dp*z32-48.0_dp*z3+40.0_dp) *ome_t_2*ome_t
  c_factor_r = c_factor_r * z3t * c1_con * c1_con
  
  ! wdsp = rho *rho *I1 *order1 + rho *rho *c1_con *m_mean *I2 *order2
  wdsp_r = rho *order1*(2.0_dp*I1 + rho*I1_r)  &
           + m_mean *order2 *rho *rho *c1_con *I2 *(2.0_dp/rho + c2_con*z3t/c1_con + I2_r/I2)

  if ( assoc ) then
     ! --- constants for iteration ---------------------------------------------
     tol = hb_tol
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-4_dp
     max_eval = 1000

     ! --- initialize mx(i,j) --------------------------------------------------
     mx_r(:,:) = 0.0_dp

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij_r(i,j) = gij_r_t1 + dij_ab(i,j)*( gij_r_t2 + dij_ab(i,j) *gij_r_t3 )
           gij_r(j,i) = gij_r(i,j)
        end do
     end do
  end if


  !-----------------------------------------------------------------------------
  ! calculate derivative d(w)/d(rho_k) for each substance k
  !-----------------------------------------------------------------------------
  do  k = 1, ncomp

     !--------------------------------------------------------------------------
     ! more abbreviations
     !--------------------------------------------------------------------------
     ! z0_rk = PI_6 * mseg(k)
     ! z1_rk = z0_rk * dhs(k)
     ! z2_rk = z1_rk * dhs(k)
     ! z3_rk = z2_rk * dhs(k)

     z1fr_rk = z1_rk(k) / z1
     z2fr_rk = z2_rk(k) / z2
     z3fr_rk = z3_rk(k) / ome

     ! --- derivative d(m_mean)/d(rho_k) ---------------------------------------
     m_rk = ( mseg(k) - m_mean ) / rho
     m_rk_r = - m_rk / rho

     !--------------------------------------------------------------------------
     ! d(w)/d(rho_k) : hard sphere contribution
     !--------------------------------------------------------------------------

     whs1_rk = whs1 * ( z1fr_rk + z2fr_rk + z3fr_rk )
     whs2_rk = whs2 * ( 3.0_dp*z2fr_rk - z3_rk(k)/z3 + 2.0_dp*z3fr_rk )
     whs3_rk = whs3 * ( 3.0_dp*z2fr_rk - 2.0_dp*z3_rk(k)/z3 )

     ! whs= (  whs1 + whs2 + ( whs3-z0 )*LOG(ome)  ) / PI_6

     whs_rk = (  whs1_rk + whs2_rk + ( whs3_rk-z0_rk(k) )*LOG(ome)  &
                    - (whs3-z0) * z3fr_rk  ) / PI_6

     whs1_rk_r = whs1_rk * (z1tfr + z2tfr + z3tfr)  &
                + whs1 * (-z1fr_rk*z1tfr - z2fr_rk*z2tfr + z3fr_rk*z3tfr )
     whs2_rk_r = whs2_rk * ( 3.0_dp*z2tfr - z3t/z3 + 2.0_dp*z3tfr )  &
                + whs2 * ( -3.0_dp *z2fr_rk*z2tfr + z3_rk(k)/z32 *z3t +2.0_dp*z3fr_rk*z3tfr )
     whs3_rk_r = whs3_rk*(3.0_dp*z2tfr-2.0_dp*z3t/z3)  &
                + whs3 * (-3.0_dp*z2fr_rk*z2tfr +2.0_dp *z3_rk(k)/z32 *z3t)
     whs_rk_r = (  whs1_rk_r + whs2_rk_r  + whs3_rk_r * LOG(ome) &
                - ( whs3_rk-z0_rk(k) )*z3tfr - (whs3*(3.0_dp*z2tfr-2.0_dp*z3t/z3)-z0t)*z3fr_rk  &
                - (whs3-z0)*z3fr_rk*z3tfr ) / PI_6

     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k) : chain term
     !--------------------------------------------------------------------------
     gij_rk_t1 = z3_rk(k) / ome2
     gij_rk_t2 = gij_t2*rho * ( z2fr_rk + 2.0_dp*z3fr_rk )
     gij_rk_t3 = gij_t3*rho * ( 2.0_dp*z2fr_rk + 3.0_dp*z3fr_rk )
     gij_rk_r_t1 = z3_rk(k) * 2.0_dp*z3t/ome3
     gij_rk_r_t2 = gij_rk_t2 *( z2tfr +2.0_dp*z3tfr )  &
                   + gij_t2*rho *( -z2fr_rk*z2tfr +2.0_dp*z3fr_rk*z3tfr )
     gij_rk_r_t3 = gij_rk_t3 *( 2.0_dp*z2tfr + 3.0_dp*z3tfr )  &
                   + gij_t3*rho *( -2.0_dp*z2fr_rk*z2tfr + 3.0_dp*z3fr_rk*z3tfr )
     do i = 1, ncomp
        gij_rk(i,i)  = gij_rk_t1  + dij_ab(i,i) * ( gij_rk_t2  + dij_ab(i,i) * gij_rk_t3 )
        gij_rk_r(i,i)= gij_rk_r_t1 +dij_ab(i,i) * ( gij_rk_r_t2 +dij_ab(i,i) * gij_rk_r_t3 )
     end do

     whc_rk = 0.0_dp
     whc_rk_r = 0.0_dp
     do i = 1, ncomp
        whc_rk = whc_rk + rhoi(i) * (1.0_dp-mseg(i)) / gij(i,i) * gij_rk(i,i)
        whc_rk_r = whc_rk_r + rhoi(i) * (1.0_dp-mseg(i)) / gij(i,i)  &
                             * ( gij_rk_r(i,i) - gij_rk(i,i)*gij_r(i,i)/gij(i,i) )
     end do
     whc_rk_r = whc_rk_r + whc_rk / rho
     whc_rk = whc_rk + ( 1.0_dp-mseg(k)) * LOG( gij(k,k) )
     whc_rk_r = whc_rk_r + ( 1.0_dp-mseg(k)) * gij_r(k,k)/gij(k,k)


     !--------------------------------------------------------------------------
     ! PC-SAFT:  d(f)/d(rho_k) : dispersion contribution
     !--------------------------------------------------------------------------

     I1_rk = z3_rk(k) * I1_z + m_rk * I1_m
     I2_rk = z3_rk(k) * I2_z + m_rk * I2_m

     I1_rk_r = z3t * z3_rk(k) * I1_z2 + m_rk * I1_r_m + m_rk_r * I1_m
     I2_rk_r = z3t * z3_rk(k) * I2_z2 + m_rk * I2_r_m + m_rk_r * I2_m

     abbrev = - 2.0_dp *PI *mseg(k)*rho
     ord1_rk = 2.0_dp *abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k) )
     ord2_rk = abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k)**2 )

     if ( lij_correction ) ord1_rk = ord1_rk + rho * ord1_lij(k)

     c1_rk = z3_rk(k) * c2_con - m_rk * c_factor
     c1_rk_r = z3_rk(k) * z3t * c3_con  - m_rk_r * c_factor   &
               - m_rk * ( 2.0_dp * z3t * c2_con / c1_con * c_factor + c_factor_r )

    ! wdsp = rho *rho *I1 *order1 + rho *rho *c1_con *m_mean *I2 *order2
     wdsp_rk = order1*rho*rho*I1_rk + ord1_rk*I1  &
          + c1_con*m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
          + ( c1_con*m_rk + c1_rk*m_mean ) * order2*rho*rho*I2
     wdsp_rk_r = order1*rho*(2.0_dp*I1_rk +rho*I1_rk_r) + ord1_rk*I1_r + ord1_rk/rho*I1  &
          + c2_con*z3t*m_mean *( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
          + c1_con*m_mean *( order2*rho*( 2.0_dp*I2_rk + rho*I2_rk_r )  &
          + ord2_rk*I2_r + ord2_rk/rho*I2 )  &
          + ( c1_con*m_rk_r + c2_con*z3t*m_rk + c1_rk_r*m_mean ) * order2*rho*rho*I2  &
          + ( c1_con*m_rk + c1_rk*m_mean ) * order2*rho*( 2.0_dp*I2 + rho*I2_r)

     !--------------------------------------------------------------------------
     ! TPT-1-association according to Chapman et al.
     !--------------------------------------------------------------------------
     whb = 0.0_dp
     whb_r = 0.0_dp
     whb_rk = 0.0_dp
     whb_rk_r = 0.0_dp

     if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           gij_rk(i,j) = gij_rk_t1  + dij_ab(i,j) * ( gij_rk_t2 + dij_ab(i,j) * gij_rk_t3 )
           gij_rk_r(i,j) = gij_rk_r_t1 + dij_ab(i,j) * ( gij_rk_r_t2 + dij_ab(i,j)* gij_rk_r_t3 )
           gij_rk(j,i) = gij_rk(i,j)
           gij_rk_r(j,i) = gij_rk_r(i,j)
        end do
     end do

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx_rk(i,ii) = 0.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum = tol + 1.0_dp
     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum_r = 0.0_dp
              sum_rk = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_ij = rhoi(j)*nhb_no(j,jj) * ass_d(i,j,ii,jj)
                    sum_r = sum_r + factor_ij * ( mx_r(j,jj)*gij(i,j)  &
                         + mx(j,jj) * ( gij_r(i,j) + gij(i,j) / rho ) )
                    sum_rk = sum_rk + factor_ij * ( mx(j,jj)*gij_rk(i,j)  &
                         + mx_rk(j,jj)*gij(i,j) )
                 end do
              end do
              do kk = 1, nhb_typ(k)
                 sum_rk = sum_rk + nhb_no(k,kk) *mx(k,kk) *gij(i,k) *ass_d(i,k,ii,kk)
              end do
              mx_r_itr(i,ii) = - mx(i,ii)**2 * sum_r
              mx_rk_itr(i,ii) = - mx(i,ii)**2 * sum_rk
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS(mx_r_itr(i,ii) - mx_r(i,ii))  &
                   + ABS(mx_rk_itr(i,ii) - mx_rk(i,ii))
              mx_r(i,ii) = mx_r_itr(i,ii)*hb_damp + mx_r(i,ii) * hb_damp_inv
              mx_rk(i,ii) = mx_rk_itr(i,ii)*hb_damp + mx_rk(i,ii) * hb_damp_inv
           end do
        end do

        if ( ass_cnt == max_eval .AND. err_sum > SQRT(tol) ) then
           WRITE (*,'(a,2G15.7)') 'F_density_rhok_rho: Max_eval violated (mx) Err_Sum= ',err_sum,tol
           exit
        end if

     end do


     ! --- calculate the hydrogen-bonding contribution -------------------------
     whb_rk = sum( nhb_no( k, 1:nhb_typ(k) ) * LOG( mx( k, 1:nhb_typ(k) ) ) )
     whb_rk_r = sum( nhb_no( k,1:nhb_typ(k) ) *mx_r( k,1:nhb_typ(k) ) /mx( k,1:nhb_typ(k) ) )
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_i = 0.5_dp * rhoi(i) * nhb_no(i,ii)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 factor_ij = factor_i * rhoi(j) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
                 factor_ij_mx =  mx(i,ii) * mx(j,jj) * factor_ij
                 whb_rk = whb_rk - factor_ij_mx * gij_rk(i,j)
                 whb_rk_r = whb_rk_r - factor_ij_mx * ( gij_rk_r(i,j) + 2.0_dp/rho*gij_rk(i,j) ) &
                     - factor_ij * ( mx(i,ii)*mx_r(j,jj) + mx_r(i,ii)*mx(j,jj) ) * gij_rk(i,j)
              end do
           end do
        end do
     end do

     end if

     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k) : summation of contributions. Polar contributions are added below
     !--------------------------------------------------------------------------

     w_rk(k) = whs_rk + whc_rk + wdsp_rk + whb_rk
     w_rk_r(k) = whs_rk_r + whc_rk_r + wdsp_rk_r + whb_rk_r

  end do

  !-----------------------------------------------------------------------------
  ! polar terms
  !-----------------------------------------------------------------------------
  wdd_rk(:) = 0.0_dp
  wqq_rk(:) = 0.0_dp
  wdq_rk(:) = 0.0_dp
  wdd_rk_r(:) = 0.0_dp
  wqq_rk_r(:) = 0.0_dp
  wdq_rk_r(:) = 0.0_dp
  if ( dipole ) call F_dd_density_rhok_rho( wdd_rk_r, wdd_rk )
  if ( qudpole ) call F_qq_density_rhok_rho( wqq_rk_r, wqq_rk )
  if ( dipole_quad ) call F_dq_density_rhok_rho( wdq_rk_r, wdq_rk )


  !-----------------------------------------------------------------------------
  ! d(f)/d(rho_k) : summation of all contributions
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

     w_rk(k) = w_rk(k) + wdd_rk(k) + wqq_rk(k) + wdq_rk(k)
     w_rk_r(k) = w_rk_r(k) + wdd_rk_r(k) + wqq_rk_r(k) + wdq_rk_r(k)

  end do

end subroutine F_density_rhok_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_density_rhok_rhol
!
!> \brief calculates second eos-derivatives to component-density
!!
!! This subroutine gives second eos-derivatives to component-density
!!       dd( F/VkT ) / d(rhoi)d(rhoj)
!! The variables are (T, rhoi), corrsponding to the NVT ensemble
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_density_rhok_rhol ( w_rkrl, wig_rkrl )

  use PARAMETERS, only: PI
  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: lij_correction
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp,ncomp), intent(out) :: w_rkrl
  real(dp), dimension(ncomp,ncomp), intent(out) :: wig_rkrl

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k, l
  integer                                    :: n_dim, m_dim
  integer                                    :: ii, jj, ll, iii, jjj, kk, lll
  real(dp)                                   :: whs_rkrl, whc_rkrl, wdsp_rkrl
  real(dp), dimension(ncomp)                 :: m_rk
  real(dp)                                   :: m_rkrl
  real(dp)                                   :: prod_z3z3_rkrl, prod_z2z3_rkrl
  real(dp)                                   :: abbrev

  real(dp)                                   :: hs_factor1, hs_factor2, hs_factor3, hs_factor4

  real(dp)                                   :: factor1, factor2, factor3
  real(dp)                                   :: z2fr_rk, z3fr_rk
  real(dp)                                   :: gij_rk_t1, gij_rk_t2, gij_rk_t3
  real(dp)                                   :: gij_rkrl_t1, gij_rkrl_t2, gij_rkrl_t3
  real(dp), allocatable, dimension(:,:,:)    :: gij_rk
  real(dp), allocatable, dimension(:,:,:,:)  :: gij_rkrl

  real(dp), allocatable, dimension(:,:)      :: q_XX, q_Xr, q_Xr_transpose, q_rkrl
  real(dp)                                   :: determinant

  real(dp)                                   :: I1_m, I2_m
  real(dp)                                   :: I1_z_m, I2_z_m
  real(dp), dimension(ncomp)                 :: I1_rk, I2_rk
  real(dp)                                   :: I1_rkrl, I2_rkrl
  real(dp), dimension(ncomp)                 :: ord1_rk, ord2_rk
  real(dp)                                   :: ord1_rkrl, ord2_rkrl
  real(dp)                                   :: ord1_lij_rkrl
  real(dp)                                   :: c2_dm
  real(dp), dimension(ncomp)                 :: c1_rk
  real(dp)                                   :: c1_rkrl
  real(dp)                                   :: chi_dm, chi_dmdeta


  integer                                    :: ass_cnt, max_eval  
  real(dp), dimension(ncomp,ncomp,nsite)     :: mx_rk
  real(dp), dimension(ncomp,nsite)           :: mx_rk_itr
  real(dp), dimension(ncomp)                 :: sum_rk
  real(dp)                                   :: k_sum
  real(dp), dimension(ncomp)                 :: err_sum
  real(dp)                                   :: tol
  real(dp)                                   :: factor_hb, factor_ij
  real(dp)                                   :: factor_ij_A, factor_ij_B
  real(dp), allocatable, dimension(:,:)      :: whb_rkrl

  real(dp), dimension(ncomp,ncomp)           :: wpolar_rkrl
  real(dp), allocatable, dimension(:,:)      :: wdd_rkrl, wqq_rkrl, wdq_rkrl

  logical                                    :: assoc_by_inversion
  !-----------------------------------------------------------------------------

  assoc_by_inversion = .true.
  
  allocate( gij_rk(ncomp,ncomp,ncomp) )
  allocate( gij_rkrl(ncomp,ncomp,ncomp,ncomp) )

  !-----------------------------------------------------------------------------
  ! quantities independent of k and l.
  !-----------------------------------------------------------------------------

  I1_m = sum(  ap_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I1 to mean segment number m_mean
  I2_m = sum(  bp_rk_aux(0:6) * z3_m_m0(0:6) )  ! derivative of I2 to mean segment number m_mean
  I1_z_m = sum( ap_rk_aux(1:6) * z3_m_m1(1:6) ) ! derivative d(I1)/d(z3) to m_mean
  I2_z_m = sum( bp_rk_aux(1:6) * z3_m_m1(1:6) ) ! derivative d(I2)/d(z3) to m_mean

  chi_dm = (8.0_dp*z3-2.0_dp*z32)/ome4  &
       - (-2.0_dp*z3**4 +12.0_dp*z33 -27.0_dp*z32+20.0_dp*z3) * ome_t_2
  chi_dmdeta = (-4.0_dp*z32+20.0_dp*z3+8.0_dp)/ome5  &
       - (2.0_dp*z33 +12.0_dp*z32-48.0_dp*z3+40.0_dp) *ome_t_2 *ome_t
  c2_dm = - 2.0_dp *c2_con *c1_con *chi_dm - c1_con *c1_con *chi_dmdeta

  hs_factor1 = 3.0_dp*z2/ome2
  hs_factor2 = 6.0_dp*z1*z2/ome3 + 3.0_dp*z23/z32/ome3  &
       + (3.0_dp*z3-1.0_dp)*(5.0_dp*z3-2.0_dp)*z23/z33/ome4  &
       + (z0-z23/z3/z3)/ome2 + 4.0_dp*z2/z3*z22/z32/ome + 6.0_dp/z3*z23/z33*LOG(ome)
  hs_factor3 = 6.0_dp*z2/z3/ome2 + 6.0_dp*z2/z32*LOG(ome)
  hs_factor4 = 3.0_dp * ( z1/ome2  &
       + z22/z32 * ( (3.0_dp*z3-1.0_dp)/ome3 - 2.0_dp/z3*LOG(ome) - 1.0_dp/ome ) )

  if ( assoc .AND. .NOT.assoc_by_inversion ) then

     ! ------ initialize mx(i,j) -----------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx_rk(:,i,ii) = 0.0_dp
        end do
     end do

     ! ------ tolerance for iterating hydrogen bonding monomer-fraction --------
     tol = hb_tol
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-4_dp
     max_eval = 1000

  end if
   
  !=============================================================================
  ! calculate some derivatives to rho_k ( d(f)/d(rho_k) )
  !=============================================================================
  DO  k = 1, ncomp

     z2fr_rk = z2_rk(k) / z2
     z3fr_rk = z3_rk(k) / ome

     gij_rk_t1 = z3_rk(k) / ome2
     gij_rk_t2 = gij_t2*rho * ( z2fr_rk + 2.0_dp*z3fr_rk )
     gij_rk_t3 = gij_t3*rho * ( 2.0_dp*z2fr_rk + 3.0_dp*z3fr_rk )
     do i = 1, ncomp
        gij_rk(k,i,i) = gij_rk_t1  + dij_ab(i,i) * ( gij_rk_t2  + dij_ab(i,i) * gij_rk_t3 )
        if ( assoc ) then
           do j = (i+1), ncomp
              gij_rk(k,i,j) = gij_rk_t1  + dij_ab(i,j) * ( gij_rk_t2  + dij_ab(i,j) * gij_rk_t3 )
              gij_rk(k,j,i) = gij_rk(k,i,j)
           end do
        end if
     end do

     m_rk(k) = ( mseg(k) - m_mean ) / rho

     I1_rk(k) = z3_rk(k) * I1_z + m_rk(k) * I1_m
     I2_rk(k) = z3_rk(k) * I2_z + m_rk(k) * I2_m

     abbrev = - 2.0_dp *PI *mseg(k)*rho
     ord1_rk(k) = 2.0_dp *abbrev *sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k) )
     ord2_rk(k) = abbrev * sum( x(1:ncomp)*mseg(1:ncomp)*sig3_ij(1:ncomp,k)*uij_t(1:ncomp,k)**2 )

     if ( lij_correction ) then
        ord1_rk(k) = ord1_rk(k) + rho * ord1_lij(k)
     end if

     c1_rk(k)= c2_con*z3_rk(k) - c1_con*c1_con * chi_dm * m_rk(k)

  end do

  !=============================================================================
  ! calculate the derivative of f to rho_k and rho_l ( dd(f)/d(rho_k)d(rho_l) )
  !=============================================================================

  DO  k = 1, ncomp

  DO  l = k, ncomp

     !--------------------------------------------------------------------------
     ! dd(f)/d(rho_k)d(rho_l) : hard sphere contribution
     !--------------------------------------------------------------------------
     prod_z3z3_rkrl = z3_rk(k)*z3_rk(l)
     prod_z2z3_rkrl = z2_rk(k)*z3_rk(l) + z2_rk(l)*z3_rk(k)
     
     gij_rkrl_t1 = 2.0_dp * prod_z3z3_rkrl / ome3
     gij_rkrl_t2 = 6.0_dp*( prod_z2z3_rkrl +3.0_dp*z2/ome*prod_z3z3_rkrl ) / ome3
     gij_rkrl_t3 = ( 4.0_dp*z2_rk(k)*z2_rk(l) + 12.0_dp*z2/ome *prod_z2z3_rkrl  &
                + 24.0_dp*z22/ome2*prod_z3z3_rkrl ) / ome3
     do i = 1, ncomp
        gij_rkrl(k,l,i,i) = gij_rkrl_t1 + dij_ab(i,i) * ( gij_rkrl_t2 + dij_ab(i,i) * gij_rkrl_t3 )
        if ( assoc ) then
           gij_rkrl(l,k,i,i) = gij_rkrl(k,l,i,i)
           do j = (i+1), ncomp
              gij_rkrl(k,l,i,j) = gij_rkrl_t1 + dij_ab(i,j) *( gij_rkrl_t2 + dij_ab(i,j) *gij_rkrl_t3 )
              gij_rkrl(k,l,j,i) = gij_rkrl(k,l,i,j)
              gij_rkrl(l,k,i,j) = gij_rkrl(k,l,i,j)
              gij_rkrl(l,k,j,i) = gij_rkrl(k,l,i,j)
           end do
        end if
     end do


     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k)d(rho_l) : hard sphere contribution
     !--------------------------------------------------------------------------
     whs_rkrl = ( 3.0_dp/ome *(z1_rk(k)*z2_rk(l)+z1_rk(l)*z2_rk(k))  &
               + hs_factor1 * (z1_rk(k)*z3_rk(l)+z1_rk(l)*z3_rk(k))  &
               + hs_factor2 * prod_z3z3_rkrl  &
               + hs_factor3 * z2_rk(k)*z2_rk(l)  &
               + hs_factor4 * prod_z2z3_rkrl + (z0_rk(k)*z3_rk(l)+z0_rk(l)*z3_rk(k)) /ome  )  /  PI_6


     !--------------------------------------------------------------------------
     ! d(f)/d(rho_k)d(rho_l) : chain term
     !--------------------------------------------------------------------------

     whc_rkrl = 0.0_dp
     DO i = 1, ncomp
        whc_rkrl = whc_rkrl + rhoi(i) * (1.0_dp-mseg(i)) / gij(i,i)  &
                                       * ( gij_rkrl(k,l,i,i) - gij_rk(k,i,i)*gij_rk(l,i,i)/gij(i,i) )
     END DO
     whc_rkrl = whc_rkrl + (1.0_dp-mseg(k)) / gij(k,k) * gij_rk(l,k,k) + (1.0_dp-mseg(l)) / gij(l,l) * gij_rk(k,l,l)


     !--------------------------------------------------------------------------
     ! PC-SAFT:  d(f)/d(rho_k)d(rho_l) : dispersion contribution
     !--------------------------------------------------------------------------

     m_rkrl = ( 2.0_dp*m_mean - mseg(k) - mseg(l) ) / rho/rho

     abbrev = 2.0_dp / m_mean * m_rk(k) * m_rk(l)
     factor1 = m_rkrl - abbrev
     factor2 = abbrev * 2.0_dp / m_mean**3
     factor3 = z3_rk(k)*m_rk(l) + z3_rk(l)*m_rk(k)
     I1_rkrl = factor1 *I1_m + factor2 *sum( ap(0:6,3)*z3_m_m0(0:6) )  &
          + factor3 *I1_z_m + prod_z3z3_rkrl *I1_z2
     I2_rkrl = factor1 *I2_m + factor2 *sum( bp(0:6,3)*z3_m_m0(0:6) )  &
          + factor3 *I2_z_m + prod_z3z3_rkrl *I2_z2

     ord1_rkrl = - 4.0_dp * PI * mseg(k)*mseg(l)*sig3_ij(k,l) * uij_t(k,l)
     ord2_rkrl =  0.5_dp * ord1_rkrl * uij_t(k,l)

     if ( lij_correction ) then
        ord1_lij_rkrl = - 2.0_dp *( mseg(k)*ord1_lij_aux(k)**3 + mseg(l)*ord1_lij_aux(l)**3 ) &
                        + 3.0_dp *( ord1_lij_aux(k)**2 * mseg(k)* m_sig_eps(k,l) &
                                  + ord1_lij_aux(l)**2 * mseg(l)* m_sig_eps(l,k) )
        do i = 1, ncomp
           ord1_lij_rkrl = ord1_lij_rkrl + 6.0_dp*x(i)*mseg(i) * ord1_lij_aux(i)  &
                             *( m_sig_eps(i,k)*m_sig_eps(i,l) + ord1_lij_aux(i)**2  &
                                - ord1_lij_aux(i)*(m_sig_eps(i,k)+m_sig_eps(i,l)) )
        end do        
        ord1_rkrl = ord1_rkrl - 2.0_dp * PI * ord1_lij_rkrl
     end if

     c1_rkrl = c3_con*prod_z3z3_rkrl + c2_dm*m_rk(l)*z3_rk(k)  &
          - c1_con* ( 2.0_dp*c1_rk(l)*chi_dm*m_rk(k)  &
          + c1_con*chi_dmdeta*z3_rk(l)*m_rk(k)+ c1_con*chi_dm*m_rkrl)

     wdsp_rkrl = ord1_rkrl*I1 +ord1_rk(k)*I1_rk(l)+ord1_rk(l)*I1_rk(k)+order1*rho*rho*I1_rkrl  &
          + ( m_rkrl*c1_con + m_rk(k)*c1_rk(l)+m_rk(l)*c1_rk(k)+ m_mean*c1_rkrl ) * order2*rho*rho*I2  &
          + (m_rk(k)*c1_con+m_mean*c1_rk(k)) * ( ord2_rk(l)*I2 + order2*rho*rho*I2_rk(l) )  &
          + (m_rk(l)*c1_con+m_mean*c1_rk(l)) * ( ord2_rk(k)*I2 + order2*rho*rho*I2_rk(k) )  &
          + m_mean*c1_con * ( ord2_rkrl*I2 +ord2_rk(k)*I2_rk(l)+ord2_rk(l)*I2_rk(k)+order2*rho*rho*I2_rkrl)


     !--------------------------------------------------------------------------
     ! second derivative of Helmholtz energy to rho_k and rho_l
     !--------------------------------------------------------------------------

     w_rkrl(k,l) = whs_rkrl + whc_rkrl + wdsp_rkrl
     w_rkrl(l,k) = w_rkrl(k,l)

  end do
  end do

  !-----------------------------------------------------------------------------
  ! d(f)/d(rho_k)d(rho_l) : TPT-1-association
  ! (solution by inversion according to Michelsen & Hendriks)
  !-----------------------------------------------------------------------------
  IF ( assoc .AND. assoc_by_inversion ) THEN

     n_dim = 0
     DO i = 1, ncomp
        DO ii = 1, nhb_typ(i)
           n_dim = n_dim + 1
        END DO
     END DO

     m_dim = ncomp

     allocate( q_XX( n_dim, n_dim ) )
     allocate( q_Xr( n_dim, m_dim ) )
     allocate( q_Xr_transpose( m_dim, n_dim ) )
     allocate( q_rkrl( m_dim, m_dim ) )
     allocate( whb_rkrl( m_dim, m_dim ) )

     iii = 0
     DO i = 1, ncomp
        DO ii = 1, nhb_typ(i)
           iii = iii + 1
           jjj = 0
           DO j = 1, ncomp
              DO jj = 1, nhb_typ(j)
                 jjj = jjj + 1
                 q_XX(iii,jjj) = - x(i) * x(j) * gij(i,j) *ass_d(i,j,ii,jj)
                 if ( iii == jjj ) q_XX(iii,jjj) = q_XX(iii,jjj) -x(i)*nhb_no(i,ii)/rho /( mx(i,ii)*nhb_no(i,ii) )**2
              END DO
           END DO
        END DO
     END DO

     lll = 0
     do l = 1, ncomp
        DO ll = 1, nhb_typ(l)
           lll = lll + 1
           q_Xr(lll,1:ncomp) = 0.0_dp
           DO i = 1, ncomp
              do ii = 1, nhb_typ(i)
                 factor_ij = x(l)*x(i)*mx(i,ii)*nhb_no(i,ii) *ass_d(i,l,ii,ll)
                 q_Xr(lll,1:ncomp) = q_Xr(lll,1:ncomp) - factor_ij *gij_rk(1:ncomp,i,l)
              end do
           END DO
           do k = 1, ncomp
           do kk = 1, nhb_typ(k)
              q_Xr(lll,k) = q_Xr(lll,k) - x(l)/rho* mx(k,kk)*nhb_no(k,kk) *gij(k,l) *ass_d(k,l,kk,ll)
           end do
           q_Xr_transpose(k,lll) = q_Xr(lll,k)
           end do
        END DO
     end do

     q_rkrl(:,:) = 0.0_dp
     DO i = 1, ncomp
        DO ii = 1, nhb_typ(i)
           do l = 1, ncomp
           DO ll = 1, nhb_typ(l)
              factor_ij = rhoi(i) * mx(i,ii)*nhb_no(i,ii)  &
                                      * mx(l,ll)*nhb_no(l,ll) * ass_d(i,l,ii,ll)
              q_rkrl(1:ncomp,l) = q_rkrl(1:ncomp,l) - factor_ij * gij_rk(1:ncomp,i,l)
           END DO
           end do
           do k = 1, ncomp
           DO kk = 1, nhb_typ(k)
              factor_ij = rhoi(i) * mx(i,ii)*nhb_no(i,ii)  &
                                      * mx(k,kk)*nhb_no(k,kk) *ass_d(i,k,ii,kk)
              q_rkrl(k,1:ncomp) = q_rkrl(k,1:ncomp) - factor_ij * gij_rk(1:ncomp,i,k)
           END DO
           end do
        END DO
     END DO

     do k = 1, ncomp
        do l = 1, ncomp
           do kk = 1, nhb_typ(k)
              do ll = 1, nhb_typ(l)
                 q_rkrl(k,l) = q_rkrl(k,l) - nhb_no(k,kk)*mx(k,kk)*nhb_no(l,ll)*mx(l,ll)*gij(k,l)*ass_d(k,l,kk,ll)
              end do
           end do
        end do
     end do
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_hb = 0.5_dp * rhoi(i) * nhb_no(i,ii) * mx(i,ii)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 factor_ij = factor_hb * rhoi(j) * mx(j,jj)*nhb_no(j,jj) * ass_d(i,j,ii,jj)
                 !JRE forall ( k = 1:ncomp,  l = 1:ncomp )
                 do l=1,ncomp
                     do k=1,ncomp
                         q_rkrl(k,l) = q_rkrl(k,l) - factor_ij * gij_rkrl(k,l,i,j)
                     enddo
                 enddo
                 !JRE end forall
              end do
           end do
        end do
     end do


     !whb_rkrl = q_rkrl - q_Xr_transpose * inv( q_XX ) * q_Xr
     call MATINVg( n_dim, m_dim, q_XX, q_Xr, determinant )  ! output q_Xr := inv( q_XX ) * q_Xr

     whb_rkrl = rho * rho * MATMUL( q_Xr_transpose, q_Xr )
     whb_rkrl = q_rkrl - whb_rkrl

!!$     !------------------------------------
!!$     deallocate( q_XX, q_Xr, q_Xr_transpose, q_rkrl, whb_rkrl )
!!$     end do
!!$     CALL cpu_time(t2)
!!$     !------------------------------------

     do k = 1, ncomp
        do l = 1, ncomp
           w_rkrl(k,l) = w_rkrl(k,l) + whb_rkrl(k,l)
        end do
     end do

     deallocate( q_XX, q_Xr, q_Xr_transpose, q_rkrl, whb_rkrl )

  end if

  !-----------------------------------------------------------------------------
  ! d(f)/d(rho_k)d(rho_l) : TPT-1-association
  ! (solution by classical iteration)
  !-----------------------------------------------------------------------------
  if ( assoc .AND. .NOT.assoc_by_inversion ) then

     ! --- initialize mx(i,j) --------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx_rk(:,i,ii) = 0.0_dp
        end do
     end do

     ! --- iterate over all components and all sites ---------------------------
     ass_cnt = 0
     err_sum(:) = tol + 1.0_dp
     do while ( sum( err_sum(1:ncomp) ) > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)

              sum_rk(1:ncomp) = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_hb = rhoi(j) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
                    sum_rk(1:ncomp) = sum_rk(1:ncomp) + factor_hb  &
                         * ( mx(j,jj) * gij_rk(1:ncomp,i,j) + gij(i,j) * mx_rk(1:ncomp,j,jj) )
                 end do
              end do
              do k = 1, ncomp
                 k_sum = 0.0_dp
                 do kk = 1, nhb_typ(k)
                    k_sum = k_sum + nhb_no(k,kk) * mx(k,kk) * gij(i,k) * ass_d(i,k,ii,kk)
                 end do
                 mx_rk_itr(i,ii) = -( mx(i,ii)*mx(i,ii) ) * ( k_sum + sum_rk(k) )
                 err_sum(k) = ABS(mx_rk_itr(i,ii) - mx_rk(k,i,ii))
                 mx_rk(k,i,ii) = mx_rk_itr(i,ii)*hb_damp + mx_rk(k,i,ii) * hb_damp_inv
              end do
              ! write (*,'(i4,4G24.14)') ass_cnt, err_sum(1:ncomp),mx_rk(1:ncomp,i,ii)

           end do
        end do

        if ( ass_cnt == max_eval .AND. sum( err_sum(1:ncomp) ) > SQRT(tol) ) then
           WRITE (*,'(a,2G15.7)') 'F_density_rhok_rhol: Max_eval violated (mx) Err_Sum= ', err_sum, tol
           exit
        end if

     end do


     ! --- calculate the hydrogen-bonding contribution -------------------------
!!$     do k = 1, ncomp
!!$        whb_rk(k) = sum ( nhb_no( k, 1:nhb_typ(k) ) * LOG( mx( k, 1:nhb_typ(k) ) ) )
!!$     end do
!!$     
!!$     do i = 1, ncomp
!!$        do ii = 1, nhb_typ(i)
!!$           factor_hb = 0.5_dp * rhoi(i) * mx(i,ii) * nhb_no(i,ii)
!!$           do j = 1, ncomp
!!$              do jj = 1, nhb_typ(j)
!!$                 factor_ij = factor_hb * rhoi(j) * mx(j,jj) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
!!$                 whb_rk(1:ncomp) = whb_rk(1:ncomp) - factor_ij * gij_rk(1:ncomp,i,j)
!!$              end do
!!$           end do
!!$        end do
!!$     end do

     allocate( whb_rkrl( ncomp, ncomp ) )

     do k = 1, ncomp
        whb_rkrl(k,1:ncomp) = 0.0_dp
        do kk = 1, nhb_typ(k)
           if ( abs( mx(k,kk) ) > 1.E-20_dp .and. nhb_no(k,kk) > 1.E-20_dp ) then
              whb_rkrl(k,1:ncomp) = whb_rkrl(k,1:ncomp) + nhb_no(k,kk) * mx_rk(1:ncomp,k,kk) / mx(k,kk)
           end if
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_hb = 0.5_dp * rhoi(i) * nhb_no(i,ii)
           do l = 1, ncomp
              factor_ij = 2.0_dp * factor_hb * mx(i,ii)  &
                   * sum( nhb_no(l,1:nhb_typ(l))*mx(l,1:nhb_typ(l))*ass_d(i,l,ii,1:nhb_typ(l)) )
              whb_rkrl(1:ncomp,l) = whb_rkrl(1:ncomp,l) - factor_ij * gij_rk(1:ncomp,i,l)
           end do
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 factor_ij = factor_hb * rhoi(j) * nhb_no(j,jj) * ass_d(i,j,ii,jj)
                 factor_ij_A = factor_ij * mx(i,ii) * mx(j,jj)
                 do l = 1, ncomp
                 factor_ij_B = factor_ij * (mx_rk(l,i,ii)*mx(j,jj)+mx(i,ii)*mx_rk(l,j,jj))
                 whb_rkrl(1:ncomp,l) = whb_rkrl(1:ncomp,l) - factor_ij_A *gij_rkrl(1:ncomp,l,i,j)  &
                                      - factor_ij_B *gij_rk(1:ncomp,i,j)
                 end do
              end do
           end do
        end do
     end do

     do k = 1, ncomp
        do l = 1, ncomp
           w_rkrl(k,l) = w_rkrl(k,l) + whb_rkrl(k,l)
        end do
     end do

     deallocate( whb_rkrl )

  end if

  
  !-----------------------------------------------------------------------------
  ! d(f)/d(rho_k)d(rho_l) : polar terms
  !-----------------------------------------------------------------------------

  wpolar_rkrl(:,:) = 0.0_dp

  if ( dipole ) then
     allocate( wdd_rkrl( ncomp, ncomp ) )
     call F_dd_density_rhok_rhol( wdd_rkrl )
     wpolar_rkrl(:,:) = wdd_rkrl(:,:)
     deallocate( wdd_rkrl )
  end if
  
  if ( qudpole ) then
     allocate( wqq_rkrl( ncomp, ncomp ) )
     call F_qq_density_rhok_rhol( wqq_rkrl )
     wpolar_rkrl(:,:) = wpolar_rkrl(:,:) + wqq_rkrl(:,:)
     deallocate( wqq_rkrl )
  end if
  
  if ( dipole_quad ) then
     allocate( wdq_rkrl( ncomp, ncomp ) )
     call F_dq_density_rhok_rhol( wdq_rkrl )
     wpolar_rkrl(:,:) = wpolar_rkrl(:,:) + wdq_rkrl(:,:)
     deallocate( wdq_rkrl )
  end if

  w_rkrl(:,:) = w_rkrl(:,:) + wpolar_rkrl(:,:)

  !-----------------------------------------------------------------------------
  ! d(f)/d(rho_k)d(rho_l) : ideal gas term
  !-----------------------------------------------------------------------------
  wig_rkrl(:,:) = 0.0_dp
  do k = 1, ncomp
     if ( rhoi(k) > 1.E-200_dp ) then
        wig_rkrl(k,k) = 1.0_dp / rhoi(k)
     else
        wig_rkrl(k,k) = 1.E200_dp
     end if
  end do

  deallocate( gij_rk )
  deallocate( gij_rkrl )

end subroutine F_density_rhok_rhol


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_temp_temp ( f_t, f_t2 )

  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: f_t
  real(dp), intent(out)                      :: f_t2

  !-----------------------------------------------------------------------------
  real(dp)                                   :: abbrev
  real(dp)                                   :: zfr1, zfr2, zfr3, zfr_ome
  real(dp)                                   :: zfr1d2, zfr2d2, zfr3d2, zfr_omed2

  integer                                    :: i, ii, j, jj

  real(dp)                                   :: z1_t0, z2_t0, z3_t0
  real(dp)                                   :: z1_t20, z2_t20, z3_t20
  real(dp)                                   :: z2to3
  real(dp)                                   :: z3_t
  real(dp)                                   :: z3_t2

  real(dp)                                   :: fhs_t, fhs_t2
  real(dp)                                   :: hs_term1, hs_term2
  real(dp)                                   :: z_term
  real(dp)                                   :: z3_50
  real(dp)                                   :: exp_term
  real(dp)                                   :: phi_0, phi_1, phi_t2
  real(dp)                                   :: one_minus_phi

  real(dp)                                   :: fhc_t, fhc_t2
  real(dp), dimension(ncomp)                 :: dhs_t, dhs_t2
  real(dp)                                   :: gij_t2_rho, gij_t3_rho
  real(dp)                                   :: gij_t_t1, gij_t_t2, gij_t_t3
  real(dp)                                   :: gij_t2_t1, gij_t2_t2, gij_t2_t3
  real(dp), dimension(ncomp,ncomp)           :: gij_dt, gij_dt2

  real(dp)                                   :: fdsp_t, fdsp_t2
  real(dp)                                   :: I1_t, I2_t, I1_t2, I2_t2
  real(dp)                                   :: c1_t, c1_t2

  real(dp)                                   :: fhb_t, fhb_t2
  integer                                    :: ass_cnt, max_eval
  real(dp)                                   :: factor_hb, factor_i, factor_ij
  real(dp)                                   :: d_sum, d_sum_t, d_sum_t2
  real(dp)                                   :: d_prod, d_prod_t, d_prod_t2
  real(dp)                                   :: dij, dij_t, dij_t2
  real(dp), dimension(ncomp,ncomp,nsite,nsite)   :: delta, del_t
  real(dp)                                   :: del_t2
  real(dp), dimension(ncomp,nsite)           :: mx_t, mx_t_itr
  real(dp)                                   :: tol, err_sum
  real(dp)                                   :: sum_t

  real(dp)                                   :: fdd_t, fdd_t2
  real(dp)                                   :: p_factor
  real(dp)                                   :: rho_2
  real(dp)                                   :: Idd2B_ij, Idd2B_z_ij, rdd2B_term, rdd2B_z_term
  real(dp)                                   :: fdd2, fdd2_t, fdd2_t2, fdd3, fdd3_t, fdd3_t2

  real(dp)                                   :: fqq_t, fqq_t2
  real(dp)                                   :: Iqq2B_ij, Iqq2B_z_ij, rqq2B_term, rqq2B_z_term
  real(dp)                                   :: fqq2, fqq2_t, fqq2_t2, fqq3, fqq3_t, fqq3_t2

  real(dp)                                   :: fdq_t, fdq_t2
  real(dp)                                   :: Idq2B_ij, Idq2B_z_ij, rdq2B_term, rdq2B_z_term
  real(dp)                                   :: fdq2, fdq2_t, fdq2_t2, fdq3, fdq3_t, fdq3_t2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! derivative of some auxilliary properties to temperature
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     dhs_t(i) = sigma(i) *(-3.0_dp*epsilon_k(i)/t/t)*0.12_dp*EXP(-3.0_dp*epsilon_k(i)/t)
     dhs_t2(i) = dhs_t(i)*3.0_dp*epsilon_k(i)/t/t  &
          + 6.0_dp*sigma(i)*epsilon_k(i)/t**3  *0.12_dp*EXP(-3.0_dp*epsilon_k(i)/t)
  end do

  ! derivative of z1 to T divided by rho
  z1_t0 = 0.0_dp
  z2_t0 = 0.0_dp
  z3_t0 = 0.0_dp
  do i = 1, ncomp
     abbrev = x(i) * mseg(i) * dhs_t(i)
     z1_t0 = z1_t0 + abbrev
     z2_t0 = z2_t0 + abbrev * 2.0_dp * dhs(i)
     z3_t0 = z3_t0 + abbrev * 3.0_dp * dhs(i) * dhs(i)
  end do
  z1_t0 = PI_6 * z1_t0
  z2_t0 = PI_6 * z2_t0
  z3_t0 = PI_6 * z3_t0
  z3_t = z3_t0 * rho


  z1_t20 = 0.0_dp
  z2_t20 = 0.0_dp
  z3_t20 = 0.0_dp
  do i = 1, ncomp
     abbrev = x(i) * mseg(i)
     z1_t20 = z1_t20 + abbrev*dhs_t2(i)
     z2_t20 = z2_t20 + abbrev*2.0_dp *( dhs_t(i)*dhs_t(i)+dhs(i)*dhs_t2(i) )
     z3_t20 = z3_t20 + abbrev*3.0_dp *( 2.0_dp*dhs(i)*dhs_t(i)*dhs_t(i)  &
                                    + dhs(i)*dhs(i)*dhs_t2(i) )
  end do
  z1_t20 = PI_6 * z1_t20
  z2_t20 = PI_6 * z2_t20
  z3_t20 = PI_6 * z3_t20
  z3_t2 = z3_t20 * rho

  zfr1 = z1_t0 / z1t
  zfr2 = z2_t0 / z2t
  zfr3 = z3_t0 / z3t
  zfr_ome = z3_t0 / ome

  zfr1d2 = z1_t20 / z1t - zfr1*zfr1
  zfr2d2 = z2_t20 / z2t - zfr2*zfr2
  zfr3d2 = z3_t20 / z3t - zfr3*zfr3
  zfr_omed2 = z3_t20 / ome + z3_t*z3_t0 / ome2

  z2to3 = z2t / z3t

  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT hard spheres to temp. (fhs_t)
  !-----------------------------------------------------------------------------
  if ( .NOT. mod_BMCSL ) then

     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome2
     z_term = z2t * z2to3 * z2to3
     ! fhs= (  hs_term1 + hs_term2 + (z_term-z0t)*LOG(ome)  ) / PI_6
     fhs_t = (  hs_term1* ( zfr1 + zfr2 + zfr_ome*rho)  &
       + hs_term2 * ( 3.0_dp*zfr2 + (3.0_dp*z3-1.0_dp)/z3t * zfr_ome )  &
       + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *LOG(ome)  &
       + (z0t-z_term)*zfr_ome*rho  )  / PI_6

     fhs_t2 = (  hs_term1 *( ( zfr1 + zfr2 + zfr_ome*rho)**2 + zfr1d2 + zfr2d2 + zfr_omed2*rho ) &
       + hs_term2 *( ( 3.0_dp*zfr2 + (3.0_dp*z3-1.0_dp)/z3t * zfr_ome )**2  &
                                 + 3.0_dp*zfr2d2 + (3.0_dp*rho-1.0_dp/z3t)*zfr_omed2 + zfr3/z3t*zfr_ome ) &
       + z_term *LOG(ome) *( ( 3.0_dp*zfr2 - 2.0_dp*zfr3 )**2 + 3.0_dp*zfr2d2 - 2.0_dp*zfr3d2 ) &
       - 2.0_dp*z_term*zfr_ome*rho*( 3.0_dp*zfr2 - 2.0_dp*zfr3)  &
       + (z0t-z_term)*zfr_omed2*rho  )  / PI_6

  else

     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi_0 = 1.0_dp + exp_term
     phi_1 = phi_0 + z3_50 * exp_term
     phi_t2 = z3_t2 * phi_1 + z3_t*z3_t * ( 2.0_dp + z3_50 ) * 50.0_dp * exp_term
     one_minus_phi = 1.0_dp - phi_0 * z3
     z_term = z2t * z2to3 * z2to3
     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome/one_minus_phi
     ! fhs = (  hs_term1 + hs_term2 + (z_term-z0t)*LOG(ome)  ) / PI_6
     fhs_t = (  hs_term1 * ( zfr1 + zfr2 + zfr_ome*rho )  &
       + hs_term2 * ( 3.0_dp*zfr2 + (2.0_dp*z3-1.0_dp)/z3t * zfr_ome + phi_1 *z3_t /one_minus_phi )  &
       + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *LOG(ome)  &
       + (z0t-z_term)*zfr_ome*rho  )  / PI_6
     fhs_t2 = (  hs_term1 *( ( zfr1 + zfr2 + zfr_ome*rho)**2 + zfr1d2 + zfr2d2 + zfr_omed2*rho ) &
       + hs_term2 *( ( 3.0_dp*zfr2 + (2.0_dp*z3-1.0_dp)/z3t * zfr_ome + phi_1 *z3_t /one_minus_phi )**2  &
       + 3.0_dp*zfr2d2 + (2.0_dp*rho-1.0_dp/z3t)*zfr_omed2 + zfr3/z3t*zfr_ome  &
       + (phi_1 *z3_t /one_minus_phi)**2 + phi_t2 /one_minus_phi) &
       + z_term *LOG(ome) *( ( 3.0_dp*zfr2 - 2.0_dp*zfr3 )**2 + 3.0_dp*zfr2d2 - 2.0_dp*zfr3d2 ) &
       - 2.0_dp*z_term*zfr_ome*rho*( 3.0_dp*zfr2 - 2.0_dp*zfr3)  &
       + (z0t-z_term)*zfr_omed2*rho  )  / PI_6

  end if

  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT of chain term to T (fhc_t)
  !-----------------------------------------------------------------------------

  fhc_t  = 0.0_dp
  fhc_t2 = 0.0_dp

  gij_t_t1 = z3_t/ome2
  gij_t_t2 = gij_t2*rho * ( zfr2 + 2.0_dp*zfr_ome*rho )
  gij_t_t3 = gij_t3*rho * ( 2.0_dp*zfr2 + 3.0_dp*zfr_ome*rho )

  gij_t2_t1 = ( z3_t2 + 2.0_dp*z3_t*zfr_ome*rho ) / ome2
  gij_t2_t2 = gij_t_t2 * ( zfr2 + 2.0_dp*zfr_ome*rho ) + gij_t2*rho * ( zfr2d2 + 2.0_dp*zfr_omed2*rho )
  gij_t2_t3 = gij_t_t3 * ( 2.0_dp*zfr2 + 3.0_dp*zfr_ome*rho ) + gij_t3*rho * ( 2.0_dp*zfr2d2 + 3.0_dp*zfr_omed2*rho )

  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3
  do i = 1, ncomp
     dij = 0.5_dp * dhs(i)
     dij_t = 0.5_dp * dhs_t(i)
     dij_t2 = 0.5_dp * dhs_t2(i)
     gij_dt(i,i) = gij_t_t1 + dij_t *gij_t2_rho  &
          + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
     gij_dt2(i,i) = gij_t2_t1 + dij_t2 *gij_t2_rho   &
          + 2.0_dp * dij_t *( gij_t_t2 + dij_t *gij_t3_rho + 2.0_dp*dij *gij_t_t3) &
          + dij *( gij_t2_t2 + 2.0_dp*dij_t2 *gij_t3_rho + dij *gij_t2_t3 )
  end do

  do i = 1, ncomp
     fhc_t = fhc_t + x(i) * (1.0_dp-mseg(i)) * gij_dt(i,i) / gij(i,i)
     fhc_t2= fhc_t2+ x(i) * (1.0_dp-mseg(i))  &
          * (gij_dt2(i,i) / gij(i,i) - (gij_dt(i,i)/gij(i,i))**2 )
  end do


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT dispersion term to T (fdsp_t)
  !-----------------------------------------------------------------------------

  I1_t = I1_z * z3_t
  I2_t = I2_z * z3_t
  I1_t2 = I1_z2 * z3_t * z3_t + I1_z * z3_t2
  I2_t2 = I2_z2 * z3_t * z3_t + I2_z * z3_t2

  c1_t = c2_con*z3_t
  c1_t2 = c3_con*z3_t*z3_t + c2_con*z3_t2

  fdsp_t  = rho * (I1_t-I1/t) * order1  &
       + rho * m_mean * (c1_t*I2+c1_con*I2_t-2.0_dp*c1_con*I2/t) * order2

  fdsp_t2 = rho * (I1_t2-2.0_dp*I1_t/t+2.0_dp*I1/t/t) * order1  &
       + rho * m_mean * order2 * ( c1_t2*I2 +2.0_dp*c1_t*I2_t -4.0_dp*c1_t*I2/t  &
       + 6.0_dp*c1_con*I2/t/t -4.0_dp*c1_con*I2_t/t +c1_con*I2_t2 )


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT association term to T (fhb_t)
  !-----------------------------------------------------------------------------

  fhb_t  = 0.0_dp
  fhb_t2 = 0.0_dp

  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           d_sum = dhs(i) + dhs(j)
           d_sum_t = dhs_t(i) + dhs_t(j)
           d_sum_t2 = dhs_t2(i) + dhs_t2(j)
           d_prod = dhs(i)*dhs(j)
           d_prod_t = dhs_t(i)*dhs(j) + dhs(i)*dhs_t(j)
           d_prod_t2 = dhs_t2(i)*dhs(j) + 2.0_dp*dhs_t(i)*dhs_t(j) + dhs(i)*dhs_t2(j)
           dij = d_prod / d_sum
           dij_t = ( d_prod_t - dij * d_sum_t ) / d_sum
           dij_t2 = - dij_t*d_sum_t /d_sum + ( d_prod_t2 -dij_t*d_sum_t -dij*d_sum_t2 ) /d_sum
           gij_dt(i,j) = gij_t_t1 + dij_t *gij_t2_rho  &
                + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
           gij_dt2(i,j) = gij_t2_t1 + dij_t2 *gij_t2_rho   &
                + 2.0_dp * dij_t *( gij_t_t2 + dij_t *gij_t3_rho + 2.0_dp*dij *gij_t_t3) &
                + dij *( gij_t2_t2 + 2.0_dp*dij_t2 *gij_t3_rho + dij *gij_t2_t3 )
           gij_dt(j,i)  = gij_dt(i,j)
           gij_dt2(j,i)  = gij_dt2(i,j)
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 delta(i,j,ii,jj) = gij(i,j)*ass_d(i,j,ii,jj)
                 del_t(i,j,ii,jj) = gij_dt(i,j)*ass_d(i,j,ii,jj) + gij(i,j)*ass_d_dt(i,j,ii,jj)
              end do
           end do
        end do
     end do


     !--- constants for iteration ----------------------------------------------
     tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-3_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-4_dp
     max_eval = 200

     !--- initialize mx_t(i,j) -------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx_t(i,ii) = 0.0_dp
        end do
     end do


     !--- iterate over all components and all sites ----------------------------
     err_sum = 10.0_dp * tol
     ass_cnt = 0

     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              sum_t = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_hb = rhoi(j)*nhb_no(j,jj)
                    sum_t = sum_t + factor_hb*( mx(j,jj) *del_t(i,j,ii,jj)  &
                         + mx_t(j,jj)*delta(i,j,ii,jj) )
                 end do
              end do
              mx_t_itr(i,ii)= - mx(i,ii)**2 * sum_t
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS( mx_t_itr(i,ii) - mx_t(i,ii) )
              mx_t(i,ii) = mx_t_itr(i,ii) * hb_damp + mx_t(i,ii) * hb_damp_inv
           end do
        end do

        ass_cnt = ass_cnt + 1

        if ( ass_cnt == max_eval ) then
           write (6,*) 'f_temp_temp: max_eval violated err_sum = ',err_sum,tol
           exit
        end if

     end do

     fhb_t = 0.0_dp
     fhb_t2 = 0.0_dp
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           ! fhb = fhb + x(i)* nhb_no(i,ii)* ( 0.5_dp * ( 1.0_dp - mx(i,ii) ) + LOG(mx(i,ii)) )
           factor_i = 0.5_dp * rhoi(i)*nhb_no(i,ii)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 del_t2 = gij_dt2(i,j)*ass_d(i,j,ii,jj)  &
                      + 2.0_dp*gij_dt(i,j)*ass_d_dt(i,j,ii,jj) + gij(i,j)*ass_d_dt2(i,j,ii,jj)
                 factor_ij = factor_i * x(j)*nhb_no(j,jj)
                 fhb_t = fhb_t - factor_ij *mx(i,ii)*mx(j,jj) * del_t(i,j,ii,jj)
                 fhb_t2 = fhb_t2 - factor_ij *( mx(i,ii)*mx(j,jj) * del_t2  &
                      + (mx_t(i,ii)*mx(j,jj)+mx(i,ii)*mx_t(j,jj)) * del_t(i,j,ii,jj) )
              end do
           end do
        end do
     end do

  end if



  if ( dipole .OR. qudpole .OR. dipole_quad ) then
     rho_2 = rho * rho
  end if

  !-----------------------------------------------------------------------------
  ! derivatives of f/kT of dipole-dipole term to temp. (fdd_t)
  !-----------------------------------------------------------------------------
  fdd_t  = 0.0_dp
  fdd_t2 = 0.0_dp

  if ( dipole ) then

     if ( abs( wdd2 ) > 1.E-50_dp ) then

        rdd2B_term = 0.0_dp
        rdd2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
              do j = 1, ncomp
                 if ( abs( my_factor(j)*x(j) ) > machine_eps ) then
                    Idd2B_ij = eps_ij(i,j) * sum( ddp4(i,j,0:4) * z3_m_m0(0:4) )
                    Idd2B_z_ij = eps_ij(i,j) * sum( ddp4(i,j,1:4) * z3_m_m1(1:4) )
                    p_factor = x(i) * x(j) * pi_dd(i,j)
                    rdd2B_term = rdd2B_term + p_factor * Idd2B_ij
                    rdd2B_z_term = rdd2B_z_term + p_factor * Idd2B_z_ij
                 end if
              end do
           end if
        end do

        fdd2 = rdd2 * rho
        fdd3 = rdd3 * rho_2
        fdd2_t = ( rdd2_z_term * z3_t - 2.0_dp * rdd2 / t - rdd2B_term / t ) * rho
        fdd3_t = ( rdd3_z_term * z3_t - 3.0_dp * rdd3 / t ) * rho_2
        fdd2_t2 = ( 6.0_dp/t * rdd2/t + rdd2_z2_term*z3_t*z3_t  &
             + rdd2_z_term *(z3_t2 -4.0_dp/t *z3_t)  &
             + 6.0_dp*rdd2B_term/t**2 - 2.0_dp*rdd2B_z_term *z3_t /t ) * rho
        fdd3_t2 = ( 12.0_dp/t *rdd3/t + rdd3_z2_term*z3_t*z3_t + rdd3_z_term *(z3_t2 -6.0_dp/t *z3_t) ) *rho_2

        fdd_t = fdd2* (fdd2*fdd2_t - 2.0_dp*fdd3*fdd2_t+fdd2*fdd3_t) / (fdd2-fdd3)**2 
        fdd_t2 = ( 2.0_dp*fdd2*fdd2_t*fdd2_t +fdd2*fdd2*fdd2_t2  &
             - 2.0_dp*fdd2_t**2 *fdd3  -2.0_dp*fdd2*fdd2_t2*fdd3 +fdd2*fdd2*fdd3_t2 )  &
             / (fdd2-fdd3)**2  + fdd_t * 2.0_dp*(fdd3_t-fdd2_t)/(fdd2-fdd3)

     end if

  end if


  !-----------------------------------------------------------------------------
  ! derivatives f/kT of quadrupole-quadrup. term to T  (fqq_t)
  !-----------------------------------------------------------------------------

  fqq_t  = 0.0_dp
  fqq_t2 = 0.0_dp

  if ( qudpole ) then

     if ( abs( wqq2 ) > 1.E-50_dp ) then

        rqq2B_term = 0.0_dp
        rqq2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
              do j = 1, ncomp
                 if ( abs( Q_factor(j)*x(j) ) > machine_eps ) then
                    Iqq2B_ij = eps_ij(i,j) * sum( qqp4(i,j,0:4) * z3_m_m0(0:4) )
                    Iqq2B_z_ij = eps_ij(i,j) * sum( qqp4(i,j,1:4) * z3_m_m1(1:4) )
                    p_factor = x(i) * x(j) * pi_qq(i,j)
                    rqq2B_term = rqq2B_term + p_factor * Iqq2B_ij
                    rqq2B_z_term = rqq2B_z_term + p_factor * Iqq2B_z_ij
                 end if
              end do
           end if
        end do

        fqq2 = rqq2 * rho
        fqq3 = rqq3 * rho_2
        fqq2_t = ( rqq2_z_term * z3_t - 2.0_dp * rqq2 / t - rqq2B_term / t ) * rho
        fqq3_t = ( rqq3_z_term * z3_t - 3.0_dp * rqq3 / t ) * rho_2
        fqq2_t2 = ( 6.0_dp/t * rqq2/t + rqq2_z2_term*z3_t*z3_t  &
             + rqq2_z_term *(z3_t2 -4.0_dp/t *z3_t)  &
             + 6.0_dp*rqq2B_term/t**2 - 2.0_dp*rqq2B_z_term *z3_t /t ) * rho
        fqq3_t2 = ( 12.0_dp/t *rqq3/t + rqq3_z2_term*z3_t*z3_t + rqq3_z_term *(z3_t2 -6.0_dp/t *z3_t) ) *rho_2

        fqq_t = fqq2* (fqq2*fqq2_t - 2.0_dp*fqq3*fqq2_t+fqq2*fqq3_t) / (fqq2-fqq3)**2 
        fqq_t2 = ( 2.0_dp*fqq2*fqq2_t*fqq2_t +fqq2*fqq2*fqq2_t2  &
             - 2.0_dp*fqq2_t**2 *fqq3  -2.0_dp*fqq2*fqq2_t2*fqq3 +fqq2*fqq2*fqq3_t2 )  &
             / (fqq2-fqq3)**2  + fqq_t * 2.0_dp*(fqq3_t-fqq2_t)/(fqq2-fqq3)

     end if

  end if


  !-----------------------------------------------------------------------------
  ! derivatives f/kT of dipole-quadruppole term to T  (fdq_t)
  !-----------------------------------------------------------------------------

  fdq_t = 0.0_dp
  fdq_t2= 0.0_dp

  if ( dipole_quad ) then

     if ( abs( wdq2 ) > 1.E-50_dp ) then

        rdq2B_term = 0.0_dp
        rdq2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*x(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*x(i) ) > machine_eps ) then
              do j = 1, ncomp
                 Idq2B_ij = eps_ij(i,j) * sum( dqp4(i,j,0:4) * z3_m_m0(0:4) )
                 Idq2B_z_ij = eps_ij(i,j) * sum( dqp4(i,j,1:4) * z3_m_m1(1:4) )
                 p_factor = x(i) * x(j) * pi_dq(i,j)
                 rdq2B_term = rdq2B_term + p_factor * Idq2B_ij
                 rdq2B_z_term = rdq2B_z_term + p_factor * Idq2B_z_ij
              end do
           end if
        end do

        fdq2 = rdq2 * rho
        fdq3 = rdq3 * rho_2
        fdq2_t = ( rdq2_z_term * z3_t - 2.0_dp * rdq2 / t - rdq2B_term / t ) * rho
        fdq3_t = ( rdq3_z_term * z3_t - 3.0_dp * rdq3 / t ) * rho_2
        fdq2_t2 = ( 6.0_dp/t * rdq2/t + rdq2_z2_term*z3_t*z3_t  &
             + rdq2_z_term *(z3_t2 -4.0_dp/t *z3_t)  &
             + 6.0_dp*rdq2B_term/t**2 - 2.0_dp*rdq2B_z_term *z3_t /t ) * rho
        fdq3_t2 = ( 12.0_dp/t *rdq3/t + rdq3_z2_term*z3_t*z3_t + rdq3_z_term *(z3_t2 -6.0_dp/t *z3_t) ) *rho_2

        fdq_t = fdq2* (fdq2*fdq2_t - 2.0_dp*fdq3*fdq2_t+fdq2*fdq3_t) / (fdq2-fdq3)**2 
        fdq_t2 = ( 2.0_dp*fdq2*fdq2_t*fdq2_t +fdq2*fdq2*fdq2_t2  &
             - 2.0_dp*fdq2_t**2 *fdq3  -2.0_dp*fdq2*fdq2_t2*fdq3 +fdq2*fdq2*fdq3_t2 )  &
             / (fdq2-fdq3)**2  + fdq_t * 2.0_dp*(fdq3_t-fdq2_t)/(fdq2-fdq3)

     end if

  end if


  !-----------------------------------------------------------------------------
  ! total derivative of fres/kT to temperature
  !-----------------------------------------------------------------------------

  f_t = fhs_t + fhc_t + fdsp_t + fhb_t + fdd_t + fqq_t + fdq_t

  !-----------------------------------------------------------------------------
  ! second derivative of fres/kT to T
  !-----------------------------------------------------------------------------

  f_t2 = fhs_t2 +fhc_t2 +fdsp_t2 +fhb_t2 +fdd_t2 +fqq_t2 +fdq_t2

end subroutine f_temp_temp

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_temp_rho ( f_res, f_t, f_tr )

  USE EOS_CONSTANTS
  use GlobConst, only: LOUD
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: f_res
  real(dp), intent(out)                      :: f_t
  real(dp), intent(out)                      :: f_tr
  !-----------------------------------------------------------------------------

  real(dp)                                   :: abbrev
  real(dp)                                   :: zfr1, zfr2, zfr3, zfr_ome
  real(dp)                                   :: z2to3


  integer                                    :: i, ii, j, jj
  real(dp), dimension(ncomp)                 :: dhs_t
  real(dp)                                   :: z1_t, z2_t, z3_t
  real(dp), dimension(ncomp,ncomp)           :: gij_tr

  real(dp)                                   :: fhs, fhs_t, fhs_tr
  real(dp)                                   :: hs_term1, hs_term2
  real(dp)                                   :: z_term
  real(dp)                                   :: z3_50
  real(dp)                                   :: exp_term
  real(dp)                                   :: phi_0, phi_1, phi_2
  real(dp)                                   :: one_minus_phi

  real(dp)                                   :: fhc, fhc_t, fhc_tr
  real(dp), dimension(ncomp,ncomp)           :: gij_r

  real(dp)                                   :: fdsp, fdsp_t, fdsp_tr
  real(dp)                                   :: I1_r, I2_r
  real(dp)                                   :: I1_t, I2_t, I1_t_r, I2_t_r
  real(dp)                                   :: c1_t, c1_tr

  real(dp)                                   :: fhb, fhb_t, fhb_tr
  integer                                    :: ass_cnt, max_eval
  real(dp)                                   :: dij, dij_t
  real(dp)                                   :: gij_t2_rho, gij_t3_rho
  real(dp)                                   :: gij_r_t1, gij_r_t2, gij_r_t3
  real(dp)                                   :: gij_t_t1, gij_t_t2, gij_t_t3
  real(dp)                                   :: gij_tr_t1, gij_tr_t2, gij_tr_t3
  real(dp), dimension(ncomp,ncomp)           :: gij_t
  real(dp), dimension(ncomp,ncomp,nsite,nsite)  :: delta, del_t, del_r, del_tr
  real(dp), dimension(ncomp,nsite)           :: mx_t, mx_r, mx_tr
  real(dp), dimension(ncomp,nsite)           :: mx_r_itr, mx_t_itr, mx_tr_itr
  real(dp)                                   :: factor_hb
  real(dp)                                   :: suma, sum_r, sum_t, sum_tr
  real(dp)                                   :: tol, err_sum

  real(dp)                                   :: fdd, fdd_t, fdd_tr
  real(dp)                                   :: fdd2, fdd2_t, fdd2_r, fdd2_tr
  real(dp)                                   :: fdd3, fdd3_t, fdd3_r, fdd3_tr
  real(dp)                                   :: p_factor
  real(dp)                                   :: rho_2
  real(dp)                                   :: diff_f
  real(dp)                                   :: Idd2B_ij, Idd2B_z_ij, rdd2B_term, rdd2B_z_term

  real(dp)                                   :: fqq, fqq_t, fqq_tr
  real(dp)                                   :: fqq2, fqq2_t, fqq2_r, fqq2_tr
  real(dp)                                   :: fqq3, fqq3_t, fqq3_r, fqq3_tr
  real(dp)                                   :: Iqq2B_ij, Iqq2B_z_ij, rqq2B_term, rqq2B_z_term

  real(dp)                                   :: fdq, fdq_t, fdq_tr
  real(dp)                                   :: fdq2, fdq2_t, fdq2_r, fdq2_tr
  real(dp)                                   :: fdq3, fdq3_t, fdq3_r, fdq3_tr
  real(dp)                                   :: Idq2B_ij, Idq2B_z_ij, rdq2B_term, rdq2B_z_term
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! derivative of some auxilliary properties to temperature
  !-----------------------------------------------------------------------------
  do i = 1, ncomp
     dhs_t(i) = sigma(i) *(-3.0_dp*epsilon_k(i)/t/t)*0.12_dp*EXP(-3.0_dp*epsilon_k(i)/t)
  end do

  z1_t = 0.0_dp
  z2_t = 0.0_dp
  z3_t = 0.0_dp
  do i = 1, ncomp
     abbrev = x(i) * mseg(i) * dhs_t(i)
     z1_t = z1_t + abbrev
     z2_t = z2_t + abbrev * 2.0_dp * dhs(i)
     z3_t = z3_t + abbrev * 3.0_dp * dhs(i) * dhs(i)
  end do
  z1_t  = PI_6 * z1_t
  z2_t  = PI_6 * z2_t
  z3_t  = PI_6 * z3_t

  zfr1 = z1_t / z1t
  zfr2 = z2_t / z2t
  zfr3 = z3_t / z3t
  zfr_ome = z3_t / ome
  
  z2to3 = z2t / z3t

  !-----------------------------------------------------------------------------
  ! derivative of f/kT hard spheres to temp. and to rho (fhs_tr)
  !-----------------------------------------------------------------------------
  if ( .NOT. mod_BMCSL ) then
     z_term = z2t * z2to3 * z2to3
     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome2
     fhs = (  hs_term1 + hs_term2 + (z_term-z0t)*LOG(ome)  ) / PI_6
     ! fhs   = ( 3.0_dp*z1t*z2/ome + z2*z2t*z2to3/ome2 + (z_term-z0t)*LOG(ome)  ) / PI_6
     fhs_t = (  hs_term1 * ( zfr1 + zfr2 + zfr_ome*rho )  &
          + hs_term2 * ( 3.0_dp*zfr2 + (3.0_dp*z3-1.0_dp)/z3t * zfr_ome )  &
          + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *LOG(ome)  &
          + (z0t-z_term)*zfr_ome*rho  )  / PI_6
     ! fhs_r = ( 3.0_dp*z1t*z2t/ome2 + z2 *z2t*z2t *(3.0_dp-z3)/ome3 + z0t*z3t/ome)  / PI_6

     fhs_tr = (3.0_dp*z2t*ome2*z1_t + 3.0_dp*(z2*z2t *(3.0_dp - z3) + z1t*ome)*ome*z2_t &
          + (2.0_dp*z22*z2t*(4.0_dp - z3) + 6.0_dp*z1t*z2*ome + z0t*ome2)*z3_t) / ome4 / PI_6
  else
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi_0 = 1.0_dp + exp_term
     phi_1 = phi_0 + z3_50 * exp_term
     phi_2 = phi_1 + ( 2.0_dp + z3_50 ) * z3_50 * exp_term
     one_minus_phi = 1.0_dp - phi_0 * z3
     z_term = z2t * z2to3 * z2to3
     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome/one_minus_phi
     fhs = (  hs_term1 + hs_term2 + (z_term-z0t)*LOG(ome)  ) / PI_6
     fhs_t = (  hs_term1 * ( zfr1 + zfr2 + zfr_ome*rho )  &
       + hs_term2 * ( 3.0_dp*zfr2 + (2.0_dp*z3-1.0_dp)/z3t * zfr_ome + phi_1 *z3_t*rho /one_minus_phi )  &
       + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *LOG(ome)  &
       + (z0t-z_term)*zfr_ome*rho  )  / PI_6
     ! fhs_r = ( 3.0_dp*z1t*z2t/ome2 + z2t**3/z3t/ome/one_minus_phi*(1.0_dp+z3/ome+rho*phi_r/one_minus_phi)  &
     !          - (z2t*z2to3*z2to3-z0t)*z3t/ome)  / PI_6
     fhs_tr = (3.0_dp*(z2t*z1_t+z2_t*z1t)/ome2 + 6.0_dp*z1t*z2t/ome3*z3_t*rho &
          + z2t**3/z3t/ome/one_minus_phi*( (1.0_dp+z3/ome + z3*phi_1/one_minus_phi)*(3.0_dp*zfr2 -zfr3  &
                        + zfr_ome*rho + phi_1 *z3_t*rho /one_minus_phi) &
          + zfr_ome*rho + z3_t*rho*z3/ome2 + z3_t*rho*phi_2/one_minus_phi + z3_t*rho*z3* (phi_1/one_minus_phi)**2) &
          - (z2t*z2to3*z2to3-z0t)*z3t/ome*(zfr3+zfr_ome*rho) - z2t*z2to3*z2to3*z3t/ome*(3.0_dp*zfr2-2.0_dp*zfr3) ) / PI_6
  end if

  !-----------------------------------------------------------------------------
  ! derivative of f/kT chain term to temp. and to rho (fhc_tr)
  !-----------------------------------------------------------------------------
  !gij_t1 = 1.0_dp/ome
  !gij_t2 = 3.0_dp*z2t/ome2
  !gij_t3 = 2.0_dp*z2t*z2/ome3

  fhc = 0.0_dp
  fhc_t = 0.0_dp
  fhc_tr = 0.0_dp

  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3

  gij_r_t1 = z3t/ome2
  gij_r_t2 = gij_t2 * ( 1.0_dp + 2.0_dp*z3/ome )
  gij_r_t3 = gij_t3 * ( 2.0_dp + 3.0_dp*z3/ome )

  gij_t_t1 = z3_t*rho/ome2
  gij_t_t2 = gij_t2*rho * ( zfr2 + 2.0_dp*zfr_ome*rho )
  gij_t_t3 = gij_t3*rho * ( 2.0_dp*zfr2 + 3.0_dp*zfr_ome*rho )

  gij_tr_t1 = z3_t/ome2 * ( 2.0_dp*z3/ome + 1.0_dp )
  abbrev = zfr_ome * ( 1.0_dp + z3/ome )
  gij_tr_t2 = gij_r_t2 * ( zfr2 + 2.0_dp*zfr_ome*rho ) + 2.0_dp * gij_t2*rho * abbrev
  gij_tr_t3 = gij_r_t3 * ( 2.0_dp*zfr2 + 3.0_dp*zfr_ome*rho ) + 3.0_dp * gij_t3*rho * abbrev
  do i = 1, ncomp
     dij = 0.5_dp * dhs(i)
     dij_t = 0.5_dp * dhs_t(i)
     gij_r(i,i) = gij_r_t1 + dij *( gij_r_t2 + dij * gij_r_t3 )
     gij_t(i,i) = gij_t_t1 + dij_t *gij_t2_rho  &
          + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
     gij_tr(i,i) = gij_tr_t1 + dij_t * gij_r_t2  &
          + dij *( gij_tr_t2 + 2.0_dp *dij_t *gij_r_t3 + dij *gij_tr_t3 )
  end do

  do i = 1, ncomp
     fhc = fhc + x(i) * (1.0_dp-mseg(i)) * LOG( gij(i,i) )
     fhc_t = fhc_t + x(i) * (1.0_dp-mseg(i)) * gij_t(i,i) / gij(i,i)
     fhc_tr = fhc_tr + x(i) * (1.0_dp-mseg(i))   &
                * ( gij_tr(i,i) / gij(i,i) - gij_t(i,i) / gij(i,i)**2 *gij_r(i,i) )
  end do

  
  !-----------------------------------------------------------------------------
  ! derivative of f/kT dispersion term to temp. and to rho (fdsp_tr)
  !-----------------------------------------------------------------------------

  I1_t = I1_z * z3_t * rho
  I2_t = I2_z * z3_t * rho
  I1_r = I1_z * z3t
  I2_r = I2_z * z3t
  I1_t_r = z3_t * ( I1_z + I1_z2 * z3 )
  I2_t_r = z3_t * ( I2_z + I2_z2 * z3 )

  c1_t = c2_con * z3_t * rho
  c1_tr = c3_con*z3_t*z3 + c2_con*z3_t

  fdsp = rho * I1 * order1 + rho * m_mean * c1_con * I2 * order2
  fdsp_t = rho * (I1_t-I1/t) * order1  &
       + rho * m_mean * (c1_t*I2+c1_con*I2_t-2.0_dp*c1_con*I2/t) * order2
  fdsp_tr = fdsp_t / rho + rho * order1 * ( I1_t_r-I1_r/t )  &
       + rho*m_mean * order2 * ( c1_tr*I2 + c1_t*I2_r + c2_con*z3t*I2_t &
       + c1_con*I2_t_r - 2.0_dp*c2_con*z3t*I2/t - 2.0_dp*c1_con*I2_r/t )

  
  !-----------------------------------------------------------------------------
  ! derivative of f/kT association term to temp. and to rho (fhb_tr)
  !-----------------------------------------------------------------------------

  fhb = 0.0_dp
  fhb_t = 0.0_dp
  fhb_tr = 0.0_dp

  if (assoc) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           dij = dhs(i)*dhs(j) / (dhs(i)+dhs(j))
           !dij_t =(dhs_t(i)*dhs(j) + dhs(i)*dhs_t(j)) / (dhs(i)+dhs(j))  &
           !     - dhs(i)*dhs(j)/(dhs(i)+dhs(j))**2  *(dhs_t(i)+dhs_t(j))
           dij_t = dij * ( dhs_t(i)/dhs(i) +dhs_t(j)/dhs(j)  &
                - (dhs_t(i)+dhs_t(j)) /(dhs(i)+dhs(j)) )
           gij_r(i,j) = gij_r_t1 + dij *( gij_r_t2 + dij * gij_r_t3 )
           gij_t(i,j) = gij_t_t1 + dij_t * gij_t2_rho  &
                + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
           gij_tr(i,j) = gij_tr_t1 + dij_t * gij_r_t2  &
                + dij *( gij_tr_t2 + 2.0_dp *dij_t *gij_r_t3 + dij *gij_tr_t3 )
           gij_r(j,i) = gij_r(i,j)
           gij_t(j,i) = gij_t(i,j)
           gij_tr(j,i) = gij_tr(i,j)
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 delta(i,j,ii,jj)=gij(i,j)*ass_d(i,j,ii,jj)
                 del_r(i,j,ii,jj) = gij_r(i,j)*ass_d(i,j,ii,jj)
                 del_t(i,j,ii,jj) = gij_t(i,j)*ass_d(i,j,ii,jj) + gij(i,j)*ass_d_dt(i,j,ii,jj)
                 del_tr(i,j,ii,jj) = gij_tr(i,j)*ass_d(i,j,ii,jj)  &
                                      + gij_r(i,j)*ass_d_dt(i,j,ii,jj)
              end do
           end do
        end do
     end do


     !--- constants for iteration ----------------------------------------------
     tol = hb_tol
     if ( z3 < 0.2_dp  ) tol = hb_tol * 1.E-1_dp
     if ( z3 < 0.01_dp ) tol = hb_tol * 1.E-2_dp
     if ( z3 < 1.E-6_dp) tol = hb_tol * 1.E-3_dp
     max_eval = 200

     !--- initialize mx_t(i,j) -------------------------------------------------
     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           mx_t(i,ii) = 0.0_dp
           mx_r(i,ii) = 0.0_dp
           mx_tr(i,ii) = 0.0_dp
        end do
     end do

     !--- iterate over all components and all sites ----------------------------
     err_sum = 10.0_dp * tol
     ass_cnt = 0

     do while ( err_sum > tol .AND. ass_cnt <= max_eval )

        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              suma = 0.0_dp
              sum_r = 0.0_dp
              sum_t = 0.0_dp
              sum_tr = 0.0_dp
              do j = 1, ncomp
                 do jj = 1, nhb_typ(j)
                    factor_hb = rhoi(j) * nhb_no(j,jj)
                    suma = suma + factor_hb * mx(j,jj) * delta(i,j,ii,jj)
                    sum_r = sum_r + factor_hb * ( mx(j,jj) *del_r(i,j,ii,jj) &
                         + mx_r(j,jj) *delta(i,j,ii,jj) )
                    sum_t = sum_t + factor_hb  &
                             * ( mx(j,jj) *del_t(i,j,ii,jj) + mx_t(j,jj) *delta(i,j,ii,jj) )
                    sum_tr = sum_tr + factor_hb * ( mx_r(j,jj)*del_t(i,j,ii,jj)  &
                         + mx(j,jj)*del_tr(i,j,ii,jj) + mx_t(j,jj)*del_r(i,j,ii,jj)  &
                         + mx_tr(j,jj)*delta(i,j,ii,jj) )
                 end do
              end do
              !mx_itr(i,ii) = 1.0_dp / ( 1.0_dp + suma )
              mx_r_itr(i,ii) = -( mx(i,ii)*mx(i,ii) )* ( suma/rho + sum_r )
              mx_t_itr(i,ii) = - mx(i,ii)**2 * sum_t
              mx_tr_itr(i,ii) = + 2.0_dp*mx(i,ii)**3 * sum_t * ( suma/rho + sum_r )  &
                               - mx(i,ii)**2 *( sum_tr + sum_t/rho )
           end do
        end do

        err_sum = 0.0_dp
        do i = 1, ncomp
           do ii = 1, nhb_typ(i)
              err_sum = err_sum + ABS(mx_r_itr(i,ii) - mx_r(i,ii))  &
                   + ABS(mx_t_itr(i,ii) - mx_t(i,ii)) + ABS(mx_tr_itr(i,ii) - mx_tr(i,ii))
              !mx(i,ii) = mx_itr(i,ii) * hb_damp + mx(i,ii) * hb_damp_inv
              mx_r(i,ii) = mx_r_itr(i,ii) * hb_damp + mx_r(i,ii) * hb_damp_inv
              mx_t(i,ii) = mx_t_itr(i,ii) * hb_damp + mx_t(i,ii) * hb_damp_inv
              mx_tr(i,ii) = mx_tr_itr(i,ii)*hb_damp + mx_tr(i,ii) * hb_damp_inv
           end do
        end do

        ass_cnt = ass_cnt + 1

        if ( ass_cnt == max_eval ) then
           if(LOUD)write (6,*) 'f_temp_rho: max_eval violated for ass_crit. err_sum,tol = ',err_sum,tol
           exit
        end if

     end do	! while( errsum > tol ...  ; err_sum = err_sum + ABS(mx_r_itr(i,ii) - mx_r(i,ii)) + ABS(mx_t_itr(i,ii) - mx_t(i,ii)) + ABS(mx_tr_itr(i,ii) - mx_tr(i,ii))


     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           fhb = fhb + x(i)*nhb_no(i,ii) *( 0.5_dp * ( 1.0_dp - mx(i,ii) ) + LOG(mx(i,ii)) )
           fhb_t = fhb_t + x(i)*nhb_no(i,ii) *mx_t(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp)
           fhb_tr = fhb_tr + x(i)*nhb_no(i,ii) *(mx_tr(i,ii)*(1.0_dp/mx(i,ii)-0.5_dp) &
                                                 - mx_r(i,ii)*mx_t(i,ii)/mx(i,ii)**2)
        end do
     end do

  end if

  
  if ( dipole .OR. qudpole .OR. dipole_quad ) then
     rho_2 = rho * rho
  end if

  !-----------------------------------------------------------------------------
  ! derivative of f/kT dipole-dipole term to temp. and to rho (fdd_tr)
  !-----------------------------------------------------------------------------

  fdd = 0.0_dp
  fdd_t = 0.0_dp
  fdd_tr = 0.0_dp
  if ( dipole ) then

     if ( abs( wdd2 ) > 1.E-50_dp ) then

        fdd2 = rdd2 * rho
        fdd3 = rdd3 * rho_2
        fdd2_r = rdd2 + rdd2_z_term * z3
        fdd3_r = ( 2.0_dp * rdd3 + rdd3_z_term * z3 ) * rho

        rdd2B_term = 0.0_dp
        rdd2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
              Idd2B_ij = eps_ij(i,i) * sum( ddp4(i,i,0:4) * z3_m_m0(0:4) )
              Idd2B_z_ij = eps_ij(i,i) * sum( ddp4(i,i,1:4) * z3_m_m1(1:4) )
              p_factor = x(i) * x(i) * pi_dd(i,i)
              rdd2B_term = rdd2B_term + p_factor * Idd2B_ij
              rdd2B_z_term = rdd2B_z_term + p_factor * Idd2B_z_ij
              do j = i+1, ncomp
                 if ( abs( my_factor(j)*x(j) ) > machine_eps ) then
                    Idd2B_ij = eps_ij(i,j) * sum( ddp4(i,j,0:4) * z3_m_m0(0:4) )
                    Idd2B_z_ij = eps_ij(i,j) * sum( ddp4(i,j,1:4) * z3_m_m1(1:4) )
                    p_factor = 2.0_dp * x(i) * x(j) * pi_dd(i,j)
                    rdd2B_term = rdd2B_term + p_factor * Idd2B_ij
                    rdd2B_z_term = rdd2B_z_term + p_factor * Idd2B_z_ij
                 end if
              end do
           end if
        end do
        fdd2_t = ( rdd2_z_term * rho * z3_t - 2.0_dp * rdd2 / t - rdd2B_term / t ) * rho
        fdd3_t = ( rdd3_z_term * rho * z3_t - 3.0_dp * rdd3 / t ) * rho_2
        fdd2_tr = ( rdd2_z2_term * z3 + 2.0_dp * rdd2_z_term ) * z3_t * rho  &
             - rdd2B_term /t - z3 /t *rdd2B_z_term - 2.0_dp / t *fdd2_r
        fdd3_tr = ( rdd3_z2_term * z3 + 3.0_dp * rdd3_z_term ) * z3_t * rho_2 - 3.0_dp / t * fdd3_r
        diff_f = fdd2 - fdd3

        fdd = fdd2 * fdd2 / diff_f
        fdd_t = fdd *( fdd2_t -2.0_dp*fdd3/fdd2*fdd2_t +fdd3_t ) / diff_f
        fdd_tr = 2.0_dp * ( fdd2_r*fdd2_t/diff_f + fdd2 * ( fdd2_r*fdd3_t - fdd3_r*fdd2_t  &
             - fdd3*fdd2_tr + 0.5_dp*fdd2*(fdd2_tr+fdd3_tr) ) / diff_f**2 + fdd_t * (fdd3_r-fdd2_r) / diff_f )

     end if

  end if

 
  !-----------------------------------------------------------------------------
  ! derivative of f/kT quadrupole-quadrupole term to temp. and to rho (fqq_tr)
  !-----------------------------------------------------------------------------

  fqq = 0.0_dp
  fqq_t = 0.0_dp
  fqq_tr = 0.0_dp

  if ( qudpole ) then

     if ( abs( wqq2 ) > 1.E-50_dp ) then

        fqq2 = rqq2 * rho
        fqq3 = rqq3 * rho_2
        fqq2_r = rqq2 + rqq2_z_term * z3
        fqq3_r = ( 2.0_dp * rqq3 + rqq3_z_term * z3 ) * rho

        rqq2B_term = 0.0_dp
        rqq2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
              Iqq2B_ij = eps_ij(i,i) * sum( qqp4(i,i,0:4) * z3_m_m0(0:4) )
              Iqq2B_z_ij = eps_ij(i,i) * sum( qqp4(i,i,1:4) * z3_m_m1(1:4) )
              p_factor = x(i) * x(i) * pi_qq(i,i)
              rqq2B_term = rqq2B_term + p_factor * Iqq2B_ij
              rqq2B_z_term = rqq2B_z_term + p_factor * Iqq2B_z_ij
              do j = i+1, ncomp
                 if ( abs( Q_factor(j)*x(j) ) > machine_eps ) then
                    Iqq2B_ij = eps_ij(i,j) * sum( qqp4(i,j,0:4) * z3_m_m0(0:4) )
                    Iqq2B_z_ij = eps_ij(i,j) * sum( qqp4(i,j,1:4) * z3_m_m1(1:4) )
                    p_factor = 2.0_dp * x(i) * x(j) * pi_qq(i,j)
                    rqq2B_term = rqq2B_term + p_factor * Iqq2B_ij
                    rqq2B_z_term = rqq2B_z_term + p_factor * Iqq2B_z_ij
                 end if
              end do
           end if
        end do
        fqq2_t = ( rqq2_z_term * rho * z3_t - 2.0_dp * rqq2 / t - rqq2B_term / t ) * rho
        fqq3_t = ( rqq3_z_term * rho * z3_t - 3.0_dp * rqq3 / t ) * rho_2
        fqq2_tr = ( rqq2_z2_term * z3 + 2.0_dp * rqq2_z_term ) * z3_t * rho  &
             - rqq2B_term /t - z3 /t *rqq2B_z_term - 2.0_dp / t *fqq2_r
        fqq3_tr = ( rqq3_z2_term * z3 + 3.0_dp * rqq3_z_term ) * z3_t * rho_2 - 3.0_dp / t * fqq3_r
        diff_f = fqq2 - fqq3

        fqq = fqq2 * fqq2 / diff_f
        fqq_t = fqq *( fqq2_t -2.0_dp*fqq3/fqq2*fqq2_t +fqq3_t ) / diff_f
        fqq_tr = 2.0_dp * ( fqq2_r*fqq2_t/diff_f + fqq2 * ( fqq2_r*fqq3_t - fqq3_r*fqq2_t  &
             - fqq3*fqq2_tr + 0.5_dp*fqq2*(fqq2_tr+fqq3_tr) ) / diff_f**2 + fqq_t * (fqq3_r-fqq2_r) / diff_f )

     end if

  end if


  !-----------------------------------------------------------------------------
  ! derivative of f/kT dipole-quadrupole term to temp. and to rho (fdq_tr)
  !-----------------------------------------------------------------------------

  fdq = 0.0_dp
  fdq_t = 0.0_dp
  fdq_tr= 0.0_dp
  if ( dipole_quad ) then

     if ( abs( wdq2 ) > 1.E-50_dp ) then

        fdq2 = rdq2 * rho
        fdq3 = rdq3 * rho_2
        fdq2_r = rdq2 + rdq2_z_term * z3
        fdq3_r = ( 2.0_dp * rdq3 + rdq3_z_term * z3 ) * rho

        rdq2B_term = 0.0_dp
        rdq2B_z_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*x(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*x(i) ) > machine_eps ) then
              Idq2B_ij = eps_ij(i,i) * sum( dqp4(i,i,0:4) * z3_m_m0(0:4) )
              Idq2B_z_ij = eps_ij(i,i) * sum( dqp4(i,i,1:4) * z3_m_m1(1:4) )
              p_factor = x(i) * x(i) * pi_dq(i,i)
              rdq2B_term = rdq2B_term + p_factor * Idq2B_ij
              rdq2B_z_term = rdq2B_z_term + p_factor * Idq2B_z_ij
              do j = i+1, ncomp
                 if ( abs( my_fac_dq(j)*x(j) ) > machine_eps .OR. abs( Q_fac_dq(j)*x(j) ) > machine_eps ) then
                    Idq2B_ij = eps_ij(i,j) * sum( dqp4(i,j,0:4) * z3_m_m0(0:4) )
                    Idq2B_z_ij = eps_ij(i,j) * sum( dqp4(i,j,1:4) * z3_m_m1(1:4) )
                    p_factor = x(i) * x(j) * ( pi_dq(i,j) + pi_dq(j,i) )
                    rdq2B_term = rdq2B_term + p_factor * Idq2B_ij
                    rdq2B_z_term = rdq2B_z_term + p_factor * Idq2B_z_ij
                 end if
              end do
           end if
        end do
        fdq2_t = ( rdq2_z_term * rho * z3_t - 2.0_dp * rdq2 / t - rdq2B_term / t ) * rho
        fdq3_t = ( rdq3_z_term * rho * z3_t - 3.0_dp * rdq3 / t ) * rho_2
        fdq2_tr = ( rdq2_z2_term * z3 + 2.0_dp * rdq2_z_term ) * z3_t * rho  &
             - rdq2B_term /t - z3 /t *rdq2B_z_term - 2.0_dp / t *fdq2_r
        fdq3_tr = ( rdq3_z2_term * z3 + 3.0_dp * rdq3_z_term ) * z3_t * rho_2 - 3.0_dp / t * fdq3_r
        diff_f = fdq2 - fdq3

        fdq = fdq2 * fdq2 / diff_f
        fdq_t = fdq *( fdq2_t -2.0_dp*fdq3/fdq2*fdq2_t +fdq3_t ) / diff_f
        fdq_tr = 2.0_dp * ( fdq2_r*fdq2_t/diff_f + fdq2 * ( fdq2_r*fdq3_t - fdq3_r*fdq2_t  &
             - fdq3*fdq2_tr + 0.5_dp*fdq2*(fdq2_tr+fdq3_tr) ) / diff_f**2 + fdq_t * (fdq3_r-fdq2_r) / diff_f )

     end if

  end if


  !-----------------------------------------------------------------------------
  ! residual Helmholtz energy f_res/kT and derivative of f_res/kT to temperature
  !-----------------------------------------------------------------------------

  f_res = fhs + fhc + fdsp + fhb + fdd + fqq + fdq

  f_t = fhs_t + fhc_t + fdsp_t + fhb_t + fdd_t + fqq_t + fdq_t

  !-----------------------------------------------------------------------------
  ! second derivative of fres/kT to temperature (T)  and to density (rho)
  !-----------------------------------------------------------------------------

  f_tr = fhs_tr +fhc_tr +fdsp_tr +fhb_tr +fdd_tr +fqq_tr +fdq_tr
  ! write (*,*) 'f_dtdr ana.',f_tr

end subroutine f_temp_rho

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_temp ( f_res, f_t )

  USE EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: sigma, epsilon_k
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: f_res
  real(dp), intent(out)                      :: f_t

  !-----------------------------------------------------------------------------
  real(dp)                                   :: abbrev
  real(dp)                                   :: zfr1, zfr2, zfr3, zfr_ome
  real(dp)                                   :: z2to3, z_term

  integer                                    :: i, ii, j, jj

  real(dp), dimension(ncomp)                 :: dhs_t
  real(dp)                                   :: z1_t, z2_t, z3_t

  real(dp)                                   :: fhs, fhs_t
  real(dp)                                   :: hs_term1, hs_term2
  real(dp)                                   :: log_ome
  real(dp)                                   :: z3_50
  real(dp)                                   :: exp_term
  real(dp)                                   :: phi_0, phi_1
  real(dp)                                   :: one_minus_phi

  real(dp)                                   :: fhc, fhc_t
  real(dp)                                   :: dij, dij_t
  real(dp)                                   :: gij_t2_rho, gij_t3_rho
  real(dp)                                   :: gij_t_t1, gij_t_t2, gij_t_t3
  real(dp), dimension(ncomp,ncomp)           :: gij_t

  real(dp)                                   :: fdsp, fdsp_t
  real(dp)                                   :: I1_t, I2_t
  real(dp)                                   :: c1_t

  real(dp)                                   :: fhb, fhb_t
  real(dp)                                   :: del_t
  real(dp)                                   :: factor_hb
  
  real(dp)                                   :: fdd, fdd_t
  real(dp)                                   :: p_factor
  real(dp)                                   :: diff_f
  real(dp)                                   :: rho_2
  real(dp)                                   :: fdd2, fdd3, fdd2_t, fdd3_t
  real(dp)                                   :: Idd2B_ij, rdd2B_term

  real(dp)                                   :: fqq, fqq_t
  real(dp)                                   :: fqq2, fqq2_t, fqq3, fqq3_t
  real(dp)                                   :: Iqq2B_ij, rqq2B_term

  real(dp)                                   :: fdq, fdq_t
  real(dp)                                   :: fdq2, fdq2_t, fdq3, fdq3_t
  real(dp)                                   :: Idq2B_ij, rdq2B_term

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! derivative of some auxilliary properties to temperature
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     dhs_t(i) = sigma(i) *(-3.0_dp*epsilon_k(i)/t/t)*0.12_dp*EXP(-3.0_dp*epsilon_k(i)/t)
  end do

  z1_t = 0.0_dp
  z2_t = 0.0_dp
  z3_t = 0.0_dp
  do i = 1, ncomp
     abbrev = x(i) * mseg(i) * dhs_t(i)
     z1_t = z1_t + abbrev
     z2_t = z2_t + abbrev * 2.0_dp * dhs(i)
     z3_t = z3_t + abbrev * 3.0_dp * dhs(i) * dhs(i)
  end do
  z1_t  = PI_6 * z1_t
  z2_t  = PI_6 * z2_t
  z3_t  = PI_6 * z3_t

  zfr1 = z1_t / z1t
  zfr2 = z2_t / z2t
  zfr3 = z3_t / z3t
  zfr_ome = z3_t / ome


  !-----------------------------------------------------------------------------
  ! 1st derivative of f/kT hard spheres to temp. (fhs_t)
  !-----------------------------------------------------------------------------
  z2to3 = z2t / z3t
  z_term = z2t * z2to3 * z2to3
  log_ome = LOG(ome)
  if ( .NOT. mod_BMCSL ) then
     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome2
     fhs = (  hs_term1 + hs_term2 + (z_term-z0t)*log_ome  ) / PI_6
     fhs_t = (  hs_term1 * ( zfr1 + zfr2 + zfr_ome*rho )  &
       + hs_term2 * ( 3.0_dp*zfr2 + (3.0_dp*z3-1.0_dp)/z3t * zfr_ome )  &
       + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *log_ome  &
       + (z0t-z_term)*zfr_ome*rho  )  / PI_6
  else
     z3_50 = 50.0_dp * z3
     exp_term = jam_fact * exp( z3_50 - 50.0_dp * eta_jam )
     phi_0 = 1.0_dp + exp_term
     phi_1 = phi_0 + z3_50 * exp_term
     one_minus_phi = 1.0_dp - phi_0 * z3
     hs_term1 = 3.0_dp*z1t*z2/ome
     hs_term2 = z2*z2t*z2to3/ome/one_minus_phi
     fhs = (  hs_term1 + hs_term2 + (z_term-z0t)*log_ome  ) / PI_6
     fhs_t = (  hs_term1 * ( zfr1 + zfr2 + zfr_ome*rho )  &
       + hs_term2 * ( 3.0_dp*zfr2 + (2.0_dp*z3-1.0_dp)/z3t * zfr_ome + phi_1 *z3_t*rho /one_minus_phi )  &
       + ( 3.0_dp*zfr2 - 2.0_dp*zfr3) * z_term *log_ome  &
       + (z0t-z_term)*zfr_ome*rho  )  / PI_6
  end if
  !-----------------------------------------------------------------------------
  ! 1st derivative of f/kT of chain term to T (fhc_t)
  !-----------------------------------------------------------------------------

  fhc = 0.0_dp
  fhc_t = 0.0_dp

  gij_t2_rho = rho * gij_t2
  gij_t3_rho = rho * gij_t3

  gij_t_t1 = z3_t*rho/ome2
  gij_t_t2 = gij_t2_rho * ( zfr2 + 2.0_dp*zfr_ome*rho )
  gij_t_t3 = gij_t3_rho * ( 2.0_dp*zfr2 + 3.0_dp*zfr_ome*rho )
  do i = 1, ncomp
     dij = 0.5_dp * dhs(i)
     dij_t = 0.5_dp * dhs_t(i)
     gij_t(i,i) = gij_t_t1 + dij_t *gij_t2_rho  &
          + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
  end do

  do i = 1, ncomp
     fhc = fhc + x(i) * (1.0_dp-mseg(i)) * LOG( gij(i,i) )
     fhc_t = fhc_t + x(i) * (1.0_dp-mseg(i)) * gij_t(i,i) / gij(i,i)
  end do


  !-----------------------------------------------------------------------------
  ! 1st derivative of f/kT dispersion term to T (fdsp_t)
  !-----------------------------------------------------------------------------

  I1_t = I1_z * z3_t * rho
  I2_t = I2_z * z3_t * rho

  c1_t = c2_con * z3_t * rho

  fdsp = rho * I1 * order1 + rho * m_mean * c1_con * I2 * order2
  fdsp_t  = rho * (I1_t-I1/t) * order1  &
       + rho * m_mean * (c1_t*I2+c1_con*I2_t-2.0_dp*c1_con*I2/t) * order2


  !-----------------------------------------------------------------------------
  ! 1st derivative of f/kT association term to T (fhb_t)
  !-----------------------------------------------------------------------------

  fhb = 0.0_dp
  fhb_t = 0.0_dp

  if ( assoc ) then

     do i = 1, ncomp
        do j = (i+1), ncomp
           dij = dhs(i)*dhs(j) / (dhs(i)+dhs(j))
           dij_t = dij * ( dhs_t(i)/dhs(i) +dhs_t(j)/dhs(j)  &
                - (dhs_t(i)+dhs_t(j)) /(dhs(i)+dhs(j)) )
           gij_t(i,j) = gij_t_t1 + dij_t * gij_t2_rho  &
                + dij *( gij_t_t2 + 2.0_dp *dij_t *gij_t3_rho + dij *gij_t_t3 )
           gij_t(j,i) = gij_t(i,j)
        end do
     end do

     do i = 1, ncomp
        do ii = 1, nhb_typ(i)
           factor_hb = 0.5_dp * rhoi(i)*nhb_no(i,ii)*mx(i,ii)
           fhb = fhb + x(i)*nhb_no(i,ii) *( 0.5_dp *( 1.0_dp - mx(i,ii) ) + LOG(mx(i,ii)) )
           do j = 1, ncomp
              do jj = 1, nhb_typ(j)
                 del_t = gij_t(i,j)*ass_d(i,j,ii,jj) + gij(i,j)*ass_d_dt(i,j,ii,jj)
                 fhb_t = fhb_t - factor_hb * x(j)*nhb_no(j,jj)*mx(j,jj) * del_t
              end do
           end do
        end do
     end do

  end if


  if ( dipole .OR. qudpole .OR. dipole_quad ) then
     rho_2 = rho * rho
  end if

  !-----------------------------------------------------------------------------
  ! derivative of f/kT dipole-dipole term to temp. (fdd_t)
  !-----------------------------------------------------------------------------

  fdd = 0.0_dp
  fdd_t = 0.0_dp

  if ( dipole ) then

     if ( abs( wdd2 ) > 1.E-50_dp ) then

        fdd2 = rdd2 * rho
        fdd3 = rdd3 * rho_2

        rdd2B_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
              Idd2B_ij = eps_ij(i,i) * sum( ddp4(i,i,0:4) * z3_m_m0(0:4) )
              p_factor = x(i) * x(i) * pi_dd(i,i)
              rdd2B_term = rdd2B_term + p_factor * Idd2B_ij
              do j = i+1, ncomp
                 if ( abs( my_factor(j)*x(j) ) > machine_eps ) then
                    Idd2B_ij = eps_ij(i,j) * sum( ddp4(i,j,0:4) * z3_m_m0(0:4) )
                    p_factor = 2.0_dp * x(i) * x(j) * pi_dd(i,j)
                    rdd2B_term = rdd2B_term + p_factor * Idd2B_ij
                 end if
              end do
           end if
        end do
        fdd2_t = ( rdd2_z_term * rho * z3_t - 2.0_dp * rdd2 / t - rdd2B_term / t ) * rho
        fdd3_t = ( rdd3_z_term * rho * z3_t - 3.0_dp * rdd3 / t ) * rho_2
        diff_f = fdd2 - fdd3

        fdd = fdd2 * fdd2 / diff_f
        fdd_t = fdd * ( fdd2_t - 2.0_dp*fdd3/fdd2*fdd2_t + fdd3_t ) / diff_f

     end if

  end if


  !-----------------------------------------------------------------------------
  ! derivative f/kT of quadrupole-quadrup. term to T  (fqq_t)
  !-----------------------------------------------------------------------------

  fqq = 0.0_dp
  fqq_t = 0.0_dp

  if ( qudpole ) then

     if ( abs( wqq2 ) > 1.E-50_dp ) then

        fqq2 = rqq2 * rho
        fqq3 = rqq3 * rho_2

        rqq2B_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
              Iqq2B_ij = eps_ij(i,i) * sum( qqp4(i,i,0:4) * z3_m_m0(0:4) )
              p_factor = x(i) * x(i) * pi_qq(i,i)
              rqq2B_term = rqq2B_term + p_factor * Iqq2B_ij
              do j = i+1, ncomp
                 if ( abs( Q_factor(j)*x(j) ) > machine_eps ) then
                    Iqq2B_ij = eps_ij(i,j) * sum( qqp4(i,j,0:4) * z3_m_m0(0:4) )
                    p_factor = 2.0_dp * x(i) * x(j) * pi_qq(i,j)
                    rqq2B_term = rqq2B_term + p_factor * Iqq2B_ij
                 end if
              end do
           end if
        end do
        fqq2_t = ( rqq2_z_term * rho * z3_t - 2.0_dp * rqq2 / t - rqq2B_term / t ) * rho
        fqq3_t = ( rqq3_z_term * rho * z3_t - 3.0_dp * rqq3 / t ) * rho_2
        diff_f = fqq2 - fqq3

        fqq = fqq2 * fqq2 / diff_f
        fqq_t = fqq * ( fqq2_t - 2.0_dp*fqq3/fqq2*fqq2_t + fqq3_t ) / diff_f

     end if

  end if


  !-----------------------------------------------------------------------------
  ! derivative f/kT of dipole-quadruppole term to T  (fdq_t)
  !-----------------------------------------------------------------------------

  fdq = 0.0_dp
  fdq_t = 0.0_dp

  if ( dipole_quad ) then

     if ( abs( wdq2 ) > 1.E-50_dp ) then

        fdq2 = rdq2 * rho
        fdq3 = rdq3 * rho_2

        rdq2B_term = 0.0_dp
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*x(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*x(i) ) > machine_eps ) then
              Idq2B_ij = eps_ij(i,i) * sum( dqp4(i,i,0:4) * z3_m_m0(0:4) )
              p_factor = x(i) * x(i) * pi_dq(i,i)
              rdq2B_term = rdq2B_term + p_factor * Idq2B_ij
              do j = i+1, ncomp
                 if ( abs( my_fac_dq(j)*x(j) ) > machine_eps .OR. abs( Q_fac_dq(j)*x(j) ) > machine_eps ) then
                    Idq2B_ij = eps_ij(i,j) * sum( dqp4(i,j,0:4) * z3_m_m0(0:4) )
                    p_factor = x(i) * x(j) * ( pi_dq(i,j) + pi_dq(j,i) )
                    rdq2B_term = rdq2B_term + p_factor * Idq2B_ij
                 end if
              end do
           end if
        end do
        fdq2_t = ( rdq2_z_term * rho * z3_t - 2.0_dp * rdq2 / t - rdq2B_term / t ) * rho
        fdq3_t = ( rdq3_z_term * rho * z3_t - 3.0_dp * rdq3 / t ) * rho_2
        diff_f = fdq2 - fdq3

        fdq = fdq2 * fdq2 / diff_f
        fdq_t = fdq * ( fdq2_t - 2.0_dp*fdq3/fdq2*fdq2_t + fdq3_t ) / diff_f

     end if

  end if


  !-----------------------------------------------------------------------------
  ! residual Helmholtz energy f_res/kT and derivative of f_res/kT to temperature
  !-----------------------------------------------------------------------------

  f_res = fhs + fhc + fdsp + fhb + fdd + fqq + fdq
  f_t = fhs_t + fhc_t + fdsp_t + fhb_t + fdd_t + fqq_t + fdq_t

end subroutine f_temp

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE f_dd_rho_rho
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_dd_rho_rho ( zdd, zdd_z )

  !-----------------------------------------------------------------------------
  real(dp), intent(IN OUT)                   :: zdd, zdd_z
  !-----------------------------------------------------------------------------
  real(dp)                                   :: fdd_z, fdd_z2
  real(dp)                                   :: rho_2, z3_fac, f3_to_2
  real(dp)                                   :: diff_factor, abbrev
  real(dp)                                   :: fdd2, fdd2_z, fdd2_z2
  real(dp)                                   :: fdd3, fdd3_z, fdd3_z2
  
  !-----------------------------------------------------------------------------

  zdd = 0.0_dp
  zdd_z = 0.0_dp

  if ( abs( rdd2 ) > 1.E-50_dp ) then

     f3_to_2 = rdd3 / rdd2
     f3_to_2 = f3_to_2 * rho
     diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fdd2 = rdd2 * rho
     fdd3 = rdd3 * rho_2
     fdd2_z = rdd2/z3t + rdd2_z_term * rho
     fdd3_z = z3_fac * rho * rdd3 + rdd3_z_term * rho_2
     fdd2_z2 = z3_fac * rdd2_z_term + rdd2_z2_term *rho
     fdd3_z2 = z3_fac /z3t  *rdd3 + 4.0_dp *rho *rdd3_z_term / z3t + rdd3_z2_term *rho_2

     ! fdd_z = fdd2 * (fdd2*fdd2_z - 2.0_dp*fdd3*fdd2_z + fdd2*fdd3_z) / (fdd2-fdd3)**2
     fdd_z = ( fdd2_z + ( fdd3_z - f3_to_2*fdd2_z ) * diff_factor ) * diff_factor

     if ( abs( fdd2 ) > 1.E-30_dp) then
        abbrev = ( fdd2_z**2 + fdd_z *(fdd3_z-fdd2_z) ) / fdd2
     else
        abbrev = 0.0_dp
     end if
     fdd_z2 = ( 2.0_dp* abbrev + ( ( fdd2_z2 + fdd3_z2 ) - 2.0_dp*f3_to_2*fdd2_z2 ) *diff_factor ) * diff_factor

     zdd = fdd_z * z3
     zdd_z = fdd_z2 * z3 + fdd_z

  end if

end subroutine f_dd_rho_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE f_qq_rho_rho
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_qq_rho_rho ( zqq, zqq_z )

  !-----------------------------------------------------------------------------
  real(dp), intent(IN OUT)                   :: zqq, zqq_z
  !-----------------------------------------------------------------------------
  real(dp)                                   :: fqq_z, fqq_z2
  real(dp)                                   :: rho_2, z3_fac, f3_to_2
  real(dp)                                   :: diff_factor, abbrev
  real(dp)                                   :: fqq2, fqq2_z, fqq2_z2
  real(dp)                                   :: fqq3, fqq3_z, fqq3_z2
  !-----------------------------------------------------------------------------

  zqq = 0.0_dp
  zqq_z = 0.0_dp

  if ( abs( rqq2 ) > 1.E-50_dp ) then

     f3_to_2 = rqq3 / rqq2
     f3_to_2 = f3_to_2 * rho
     diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fqq2 = rqq2 * rho
     fqq3 = rqq3 * rho_2
     fqq2_z = rqq2/z3t + rqq2_z_term * rho
     fqq3_z = z3_fac * rho * rqq3 + rqq3_z_term * rho_2
     fqq2_z2 = z3_fac * rqq2_z_term + rqq2_z2_term *rho
     fqq3_z2 = z3_fac /z3t  *rqq3 + 4.0_dp *rho *rqq3_z_term / z3t + rqq3_z2_term *rho_2

     fqq_z = ( fqq2_z + ( fqq3_z - f3_to_2*fqq2_z ) * diff_factor ) * diff_factor

     if ( abs( fqq2 ) > 1.E-30_dp) then
        abbrev = ( fqq2_z**2 + fqq_z *(fqq3_z-fqq2_z) ) / fqq2
     else
        abbrev = 0.0_dp
     end if
     fqq_z2 = ( 2.0_dp* abbrev + ( ( fqq2_z2 + fqq3_z2 ) - 2.0_dp*f3_to_2*fqq2_z2 ) *diff_factor ) * diff_factor

     zqq = fqq_z * z3
     zqq_z = fqq_z2 * z3 + fqq_z

  end if

end subroutine f_qq_rho_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE f_dq_rho_rho
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_dq_rho_rho ( zdq, zdq_z )

  !-----------------------------------------------------------------------------
  real(dp), intent(IN OUT)                   :: zdq, zdq_z
  !-----------------------------------------------------------------------------
  real(dp)                                   :: fdq_z, fdq_z2
  real(dp)                                   :: rho_2, z3_fac, f3_to_2
  real(dp)                                   :: diff_factor, abbrev
  real(dp)                                   :: fdq2, fdq2_z, fdq2_z2
  real(dp)                                   :: fdq3, fdq3_z, fdq3_z2
  !-----------------------------------------------------------------------------

  zdq = 0.0_dp
  zdq_z = 0.0_dp

  if ( abs( rdq2 ) > 1.E-50_dp ) then

     f3_to_2 = rdq3 / rdq2
     f3_to_2 = f3_to_2 * rho
     diff_factor = 1.0_dp / ( 1.0_dp - f3_to_2 )

     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fdq2 = rdq2 * rho
     fdq3 = rdq3 * rho_2
     fdq2_z = rdq2/z3t + rdq2_z_term * rho
     fdq3_z = z3_fac * rho * rdq3 + rdq3_z_term * rho_2
     fdq2_z2 = z3_fac * rdq2_z_term + rdq2_z2_term *rho
     fdq3_z2 = z3_fac /z3t  *rdq3 + 4.0_dp *rho *rdq3_z_term / z3t + rdq3_z2_term *rho_2

     fdq_z = ( fdq2_z + ( fdq3_z - f3_to_2*fdq2_z ) * diff_factor ) * diff_factor

     if ( abs( fdq2 ) > 1.E-30_dp) then
        abbrev = ( fdq2_z**2 + fdq_z *(fdq3_z-fdq2_z) ) / fdq2
     else
        abbrev = 0.0_dp
     end if
     fdq_z2 = ( 2.0_dp* abbrev + ( ( fdq2_z2 + fdq3_z2 ) - 2.0_dp*f3_to_2*fdq2_z2 ) *diff_factor ) * diff_factor

     zdq = fdq_z * z3
     zdq_z = fdq_z2 * z3 + fdq_z

  end if

end subroutine f_dq_rho_rho


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_dd_density_rhok_rho( wdd_rk_r, wdd_rk )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in out)  :: wdd_rK_r
  real(dp), dimension(ncomp), intent(in out)  :: wdd_rK

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k

  real(dp)                                   :: rho_2
  real(dp)                                   :: factorA, factorB
  real(dp)                                   :: factor1, factor2
  real(dp)                                   :: wdd2_r, wdd3_r
  real(dp)                                   :: wdd2_rk_r, wdd3_rk_r
  real(dp)                                   :: rdd2_rk, rdd3_rk
  real(dp)                                   :: wdd2_rk, wdd3_rk
  !-----------------------------------------------------------------------------

  wdd_rk(:) = 0.0_dp
  wdd_rk_r(:) = 0.0_dp

  if ( abs( wdd2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  rho_2 = rho*rho
  factorA = wdd2 / (wdd2-wdd3)
  factorA = factorA * factorA
  factorB = factorA * ( 1.0_dp - 2.0_dp * wdd3 / wdd2 )

  factor1 = 2.0_dp / ( wdd2 - wdd3 )
  factor2 = 2.0_dp * factorA / wdd2

  !-----------------------------------------------------------------------------
  ! first derivatives to density rho
  !-----------------------------------------------------------------------------

  wdd2_r = ( 2.0_dp * rdd2 + rdd2_z_term * z3 ) * rho
  wdd3_r = ( 3.0_dp * rdd3 + rdd3_z_term * z3 ) * rho_2


  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k. rdd2_rk and rdd3_rk are defined
  ! as rdd2_rk = rdd2_rk_(usually) * rho = d(wdd)/d(rho_k) / rho.
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

     rdd2_rk = rdd2_z_term * rho * z3_rk(k)
     rdd3_rk = rdd3_z_term * rho * z3_rk(k)
     if ( abs( my_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
              rdd2_rk = rdd2_rk + 2.0_dp *x(i) *pi_dd(i,k) *Idd2(i,k)
              rdd3_rk = rdd3_rk + 3.0_dp *x(i)*x(i) *psi_dd(i,i,k) *Idd3(i,i,k)
              do j = i+1, ncomp
                 rdd3_rk = rdd3_rk + 6.0_dp *x(i)*x(j) *psi_dd(i,j,k) *Idd3(i,j,k)
              end do
           end if
        end do
     end if

     wdd2_rk = rdd2_rk * rho
     wdd3_rk = rdd3_rk * rho_2

     ! wdd_rk = ( wdd2*wdd2_rk - 2.0_dp*wdd3*wdd2_rk+wdd2*wdd3_rk ) * wdd2 / (wdd2-wdd3)**2
     wdd_rk(k) = factorB * wdd2_rk + factorA * wdd3_rk

     !-----------------------------------------------------------------------------
     ! second derivatives to rho_k and to rho
     !-----------------------------------------------------------------------------

     wdd2_rk_r = rdd2_rk + rdd2_z_term*rho *z3_rk(k) + wdd2_z2_term *z3t *z3_rk(k)
     wdd3_rk_r = ( 2.0_dp *rdd3_rk + rdd3_z_term *rho *z3_rk(k) ) *rho + wdd3_z2_term *z3t *z3_rk(k)
     if ( abs( my_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_factor(i)*rhoi(i) ) > machine_eps ) then
              wdd2_rk_r = wdd2_rk_r + 2.0_dp *rhoi(i) *pi_dd(i,k) *Idd2_z(i,k) *z3t
              wdd3_rk_r = wdd3_rk_r + 3.0_dp *rhoi(i)*rhoi(i) *psi_dd(i,i,k) *Idd3_z(i,i,k) *z3t
              do j = i+1, ncomp
                 wdd3_rk_r = wdd3_rk_r + 6.0_dp *rhoi(i)*rhoi(j) *psi_dd(i,j,k) *Idd3_z(i,j,k) *z3t
              end do
           end if
        end do
     end if


     wdd_rk_r(k) = factor1 * wdd2_r*wdd2_rk + factorB * wdd2_rk_r + factorA * wdd3_rk_r  &
                 + factor2 * ( wdd2_r*wdd3_rk - wdd3_r*wdd2_rk )  &
                 + factor1 * wdd_rk(k) * ( wdd3_r - wdd2_r )

  end do

  end if

end subroutine F_dd_density_rhok_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_qq_density_rhok_rho( wqq_rk_r, wqq_rk )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in out)  :: wqq_rK_r
  real(dp), dimension(ncomp), intent(in out)  :: wqq_rK

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k

  real(dp)                                   :: rho_2
  real(dp)                                   :: factorA, factorB
  real(dp)                                   :: factor1, factor2
  real(dp)                                   :: wqq2_r, wqq3_r
  real(dp)                                   :: wqq2_rk_r, wqq3_rk_r
  real(dp)                                   :: rqq2_rk, rqq3_rk
  real(dp)                                   :: wqq2_rk, wqq3_rk
  !-----------------------------------------------------------------------------

  wqq_rk(:) = 0.0_dp
  wqq_rk_r(:) = 0.0_dp

  if ( abs( wqq2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  rho_2 = rho*rho
  factorA = wqq2 / (wqq2-wqq3)
  factorA = factorA * factorA
  factorB = factorA * ( 1.0_dp - 2.0_dp * wqq3 / wqq2 )

  factor1 = 2.0_dp / ( wqq2 - wqq3 )
  factor2 = 2.0_dp * factorA / wqq2

  !-----------------------------------------------------------------------------
  ! first derivatives to density rho
  !-----------------------------------------------------------------------------

  wqq2_r = ( 2.0_dp * rqq2 + rqq2_z_term * z3 ) * rho
  wqq3_r = ( 3.0_dp * rqq3 + rqq3_z_term * z3 ) * rho_2


  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k. rqq2_rk and rqq3_rk are defined
  ! as rqq2_rk = rqq2_rk_(usually) * rho = d(wqq)/d(rho_k) / rho.
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

     rqq2_rk = rqq2_z_term * rho * z3_rk(k)
     rqq3_rk = rqq3_z_term * rho * z3_rk(k)
     if ( abs( Q_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
              rqq2_rk = rqq2_rk + 2.0_dp *x(i) *pi_qq(i,k) *Iqq2(i,k)
              rqq3_rk = rqq3_rk + 3.0_dp *x(i)*x(i) *psi_qq(i,i,k) *Iqq3(i,i,k)
              do j = i+1, ncomp
                 rqq3_rk = rqq3_rk + 6.0_dp *x(i)*x(j) *psi_qq(i,j,k) *Iqq3(i,j,k)
              end do
           end if
        end do
     end if

     wqq2_rk = rqq2_rk * rho
     wqq3_rk = rqq3_rk * rho_2

     ! wqq_rk = ( wqq2*wqq2_rk - 2.0_dp*wqq3*wqq2_rk+wqq2*wqq3_rk ) * wqq2 / (wqq2-wqq3)**2
     wqq_rk(k) = factorB * wqq2_rk + factorA * wqq3_rk

     !-----------------------------------------------------------------------------
     ! second derivatives to rho_k and to rho
     !-----------------------------------------------------------------------------

     wqq2_rk_r = rqq2_rk + rqq2_z_term*rho *z3_rk(k) + wqq2_z2_term *z3t *z3_rk(k)
     wqq3_rk_r = ( 2.0_dp *rqq3_rk + rqq3_z_term *rho *z3_rk(k) ) *rho + wqq3_z2_term *z3t *z3_rk(k)
     if ( abs( Q_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( Q_factor(i)*rhoi(i) ) > machine_eps ) then
              wqq2_rk_r = wqq2_rk_r + 2.0_dp *rhoi(i) *pi_qq(i,k) *Iqq2_z(i,k) *z3t
              wqq3_rk_r = wqq3_rk_r + 3.0_dp *rhoi(i)*rhoi(i) *psi_qq(i,i,k) *Iqq3_z(i,i,k) *z3t
              do j = i+1, ncomp
                 wqq3_rk_r = wqq3_rk_r + 6.0_dp *rhoi(i)*rhoi(j) *psi_qq(i,j,k) *Iqq3_z(i,j,k) *z3t
              end do
           end if
        end do
     end if


     wqq_rk_r(k) = factor1 * wqq2_r*wqq2_rk + factorB * wqq2_rk_r + factorA * wqq3_rk_r  &
                 + factor2 * ( wqq2_r*wqq3_rk - wqq3_r*wqq2_rk )  &
                 + factor1 * wqq_rk(k) * ( wqq3_r - wqq2_r )

  end do

  end if

end subroutine F_qq_density_rhok_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_dq_density_rhok_rho( wdq_rk_r, wdq_rk )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in out)  :: wdq_rK_r
  real(dp), dimension(ncomp), intent(in out)  :: wdq_rK

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k

  real(dp)                                   :: rho_2
  real(dp)                                   :: factorA, factorB
  real(dp)                                   :: factor1, factor2
  real(dp)                                   :: wdq2_r, wdq3_r
  real(dp)                                   :: wdq2_rk_r, wdq3_rk_r
  real(dp)                                   :: rdq2_rk, rdq3_rk
  real(dp)                                   :: wdq2_rk, wdq3_rk
  !-----------------------------------------------------------------------------

  wdq_rk(:) = 0.0_dp
  wdq_rk_r(:) = 0.0_dp

  if ( abs( wdq2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  rho_2 = rho*rho
  factorA = wdq2 / (wdq2-wdq3)
  factorA = factorA * factorA
  factorB = factorA * ( 1.0_dp - 2.0_dp * wdq3 / wdq2 )

  factor1 = 2.0_dp / ( wdq2 - wdq3 )
  factor2 = 2.0_dp * factorA / wdq2

  !-----------------------------------------------------------------------------
  ! first derivatives to density rho
  !-----------------------------------------------------------------------------

  wdq2_r = ( 2.0_dp * rdq2 + rdq2_z_term * z3 ) * rho
  wdq3_r = ( 3.0_dp * rdq3 + rdq3_z_term * z3 ) * rho_2


  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k. rdq2_rk and rdq3_rk are defined
  ! as rdq2_rk = rdq2_rk_(usually) * rho = d(wdq)/d(rho_k) / rho.
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

     rdq2_rk = rdq2_z_term * rho * z3_rk(k)
     rdq3_rk = rdq3_z_term * rho * z3_rk(k)
     if ( abs( my_fac_dq(k) ) > machine_eps .OR. abs( Q_fac_dq(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*x(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*x(i) ) > machine_eps ) then
              rdq2_rk = rdq2_rk + x(i) *( pi_dq(i,k) + pi_dq(k,i) ) *Idq2(i,k)
              do j = 1, ncomp
                 rdq3_rk = rdq3_rk + x(i)*x(j)  &
                      *( psi_dq(i,j,k) +psi_dq(i,k,j) +psi_dq(k,i,j) ) *Idq3(i,j,k)
              end do
           end if
        end do
     end if

     wdq2_rk = rdq2_rk * rho
     wdq3_rk = rdq3_rk * rho_2

     ! wdq_rk = ( wdq2*wdq2_rk - 2.0_dp*wdq3*wdq2_rk+wdq2*wdq3_rk ) * wdq2 / (wdq2-wdq3)**2
     wdq_rk(k) = factorB * wdq2_rk + factorA * wdq3_rk

     !-----------------------------------------------------------------------------
     ! second derivatives to rho_k and to rho
     !-----------------------------------------------------------------------------

     wdq2_rk_r = rdq2_rk + rdq2_z_term*rho *z3_rk(k) + wdq2_z2_term *z3t *z3_rk(k)
     wdq3_rk_r = ( 2.0_dp *rdq3_rk + rdq3_z_term *rho *z3_rk(k) ) *rho + wdq3_z2_term *z3t *z3_rk(k)
     if ( abs( my_fac_dq(k) ) > machine_eps .OR. abs( Q_fac_dq(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*x(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*x(i) ) > machine_eps ) then
              wdq2_rk_r = wdq2_rk_r + rhoi(i) *( pi_dq(i,k)+pi_dq(k,i) ) *Idq2_z(i,k) *z3t
              do j = 1, ncomp
                 wdq3_rk_r = wdq3_rk_r + rhoi(i)*rhoi(j)  &
                      *( psi_dq(i,j,k) +psi_dq(i,k,j) +psi_dq(k,i,j) ) *Idq3_z(i,j,k) *z3t
              end do
           end if
        end do
     end if


     wdq_rk_r(k) = factor1 * wdq2_r*wdq2_rk + factorB * wdq2_rk_r + factorA * wdq3_rk_r  &
                 + factor2 * ( wdq2_r*wdq3_rk - wdq3_r*wdq2_rk )  &
                 + factor1 * wdq_rk(k) * ( wdq3_r - wdq2_r )

  end do

  end if

end subroutine F_dq_density_rhok_rho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_dd_density_rhok_rhol ( wdd_rkrl )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp,ncomp), INTENT(IN OUT)  :: wdd_rkrl

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k, l

  real(dp)                                   :: abbrev1, abbrev2
  real(dp)                                   :: factorA
  real(dp)                                   :: factor1, factor2, factor3, factor4
  real(dp), dimension(ncomp)                 :: wdd2_z_k
  real(dp), dimension(ncomp)                 :: wdd3_z_k
  real(dp), dimension(ncomp)                 :: wdd_rk
  real(dp), dimension(ncomp)                 :: wdd2_rk
  real(dp), dimension(ncomp)                 :: wdd3_rk
  real(dp)                                   :: wdd2_rkrl, wdd3_rkrl
  !-----------------------------------------------------------------------------

  if ( abs( wdd2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  factor1 = 2.0_dp / ( wdd2 - wdd3 )
  factor2 = 0.5_dp * factor1 * factor1 * wdd2
  factor3 = 0.5_dp * factor2 * ( wdd2 - 2.0_dp * wdd3 )
  factor4 = 0.5_dp * factor2 * wdd2

  factorA = wdd2 / ( wdd2 - wdd3 )
  factorA = factorA * factorA

  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k
  !-----------------------------------------------------------------------------
  do k = 1, ncomp

     wdd2_rk(k) = 0.0_dp
     wdd3_rk(k) = 0.0_dp
     wdd2_z_k(k) = 0.0_dp
     wdd3_z_k(k) = 0.0_dp
     if ( abs( my_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_factor(i)*rhoi(i) ) > machine_eps ) then
              abbrev1 = rhoi(i) * pi_dd(i,k)
              wdd2_rk(k) = wdd2_rk(k) + abbrev1 * Idd2(i,k)
              wdd2_z_k(k) = wdd2_z_k(k) + abbrev1 * Idd2_z(i,k)
              abbrev2 = 0.5_dp * rhoi(i) * rhoi(i) * psi_dd(i,i,k)
              wdd3_rk(k) = wdd3_rk(k) + abbrev2 * Idd3(i,i,k)
              wdd3_z_k(k) = wdd3_z_k(k) + abbrev2 * Idd3_z(i,i,k)
              do j = i+1, ncomp
                 abbrev2 = rhoi(i) * rhoi(j) * psi_dd(i,j,k)
                 wdd3_rk(k) = wdd3_rk(k) + abbrev2 * Idd3(i,j,k)
                 wdd3_z_k(k) = wdd3_z_k(k) + abbrev2 * Idd3_z(i,j,k)
              end do
           end if
        end do
     end if
     wdd2_rk(k) = 2.0_dp * wdd2_rk(k) + wdd2_z_term * z3_rk(k)
     wdd3_rk(k) = 6.0_dp * wdd3_rk(k) + wdd3_z_term * z3_rk(k)
     wdd2_z_k(k) = 2.0_dp * wdd2_z_k(k)
     wdd3_z_k(k) = 6.0_dp * wdd3_z_k(k)

     wdd_rk(k) = factor3 * wdd2_rk(k) + factorA * wdd3_rk(k)

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives to rho_k and to rho_l
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = k, ncomp

     wdd2_rkrl = 2.0_dp * pi_dd(k,l) * Idd2(k,l) + wdd2_z_k(k)*z3_rk(l) + wdd2_z_k(l)*z3_rk(k)  &
                 + wdd2_z2_term * z3_rk(k)*z3_rk(l)
     wdd3_rkrl = wdd3_z_k(k)*z3_rk(l) + wdd3_z_k(l)*z3_rk(k) + wdd3_z2_term * z3_rk(k)*z3_rk(l) &
                 + 6.0_dp * sum( rhoi(1:ncomp) * psi_dd(1:ncomp,k,l) * Idd3(1:ncomp,k,l) )

     wdd_rkrl(k,l) = factor1 * wdd2_rk(l)*wdd2_rk(k)  &
                 + factor2 * (wdd3_rk(k)*wdd2_rk(l)-wdd3_rk(l)*wdd2_rk(k))  &
                 + factor3 * wdd2_rkrl + factor4 * wdd3_rkrl  &
                 + factor1 * wdd_rk(k) * ( wdd3_rk(l) - wdd2_rk(l) )
     wdd_rkrl(l,k) = wdd_rkrl(k,l)

  end do
  end do

  end if

end subroutine F_dd_density_rhok_rhol




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_qq_density_rhok_rhol ( wqq_rkrl )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp,ncomp), INTENT(IN OUT)  :: wqq_rkrl

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k, l

  real(dp)                                   :: abbrev1, abbrev2
  real(dp)                                   :: factorA
  real(dp)                                   :: factor1, factor2, factor3, factor4
  real(dp), dimension(ncomp)                 :: wqq2_z_k
  real(dp), dimension(ncomp)                 :: wqq3_z_k
  real(dp), dimension(ncomp)                 :: wqq_rk
  real(dp), dimension(ncomp)                 :: wqq2_rk
  real(dp), dimension(ncomp)                 :: wqq3_rk
  real(dp)                                   :: wqq2_rkrl, wqq3_rkrl
  !-----------------------------------------------------------------------------

  if ( abs( wqq2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  factor1 = 2.0_dp / ( wqq2 - wqq3 )
  factor2 = 0.5_dp * factor1 * factor1 * wqq2
  factor3 = 0.5_dp * factor2 * ( wqq2 - 2.0_dp * wqq3 )
  factor4 = 0.5_dp * factor2 * wqq2

  factorA = wqq2 / ( wqq2 - wqq3 )
  factorA = factorA * factorA

  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k
  !-----------------------------------------------------------------------------
  do k = 1, ncomp

     wqq2_rk(k) = 0.0_dp
     wqq3_rk(k) = 0.0_dp
     wqq2_z_k(k) = 0.0_dp
     wqq3_z_k(k) = 0.0_dp
     if ( abs( Q_factor(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( Q_factor(i)*rhoi(i) ) > machine_eps ) then
              abbrev1 = rhoi(i) * pi_qq(i,k)
              wqq2_rk(k) = wqq2_rk(k) + abbrev1 * Iqq2(i,k)
              wqq2_z_k(k) = wqq2_z_k(k) + abbrev1 * Iqq2_z(i,k)
              abbrev2 = 0.5_dp * rhoi(i) * rhoi(i) * psi_qq(i,i,k)
              wqq3_rk(k) = wqq3_rk(k) + abbrev2 * Iqq3(i,i,k)
              wqq3_z_k(k) = wqq3_z_k(k) + abbrev2 * Iqq3_z(i,i,k)
              do j = i+1, ncomp
                 abbrev2 = rhoi(i) * rhoi(j) * psi_qq(i,j,k)
                 wqq3_rk(k) = wqq3_rk(k) + abbrev2 * Iqq3(i,j,k)
                 wqq3_z_k(k) = wqq3_z_k(k) + abbrev2 * Iqq3_z(i,j,k)
              end do
           end if
        end do
     end if
     wqq2_rk(k) = 2.0_dp * wqq2_rk(k) + wqq2_z_term * z3_rk(k)
     wqq3_rk(k) = 6.0_dp * wqq3_rk(k) + wqq3_z_term * z3_rk(k)
     wqq2_z_k(k) = 2.0_dp * wqq2_z_k(k)
     wqq3_z_k(k) = 6.0_dp * wqq3_z_k(k)

     wqq_rk(k) = factor3 * wqq2_rk(k) + factorA * wqq3_rk(k)

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives to rho_k and to rho_l
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = k, ncomp

     wqq2_rkrl = 2.0_dp * pi_qq(k,l) * Iqq2(k,l) + wqq2_z_k(k)*z3_rk(l) + wqq2_z_k(l)*z3_rk(k)  &
                 + wqq2_z2_term * z3_rk(k)*z3_rk(l)
     wqq3_rkrl = wqq3_z_k(k)*z3_rk(l) + wqq3_z_k(l)*z3_rk(k) + wqq3_z2_term * z3_rk(k)*z3_rk(l) &
                 + 6.0_dp * sum( rhoi(1:ncomp) * psi_qq(1:ncomp,k,l) * Iqq3(1:ncomp,k,l) )

     wqq_rkrl(k,l) = factor1 * wqq2_rk(l)*wqq2_rk(k)  &
                 + factor2 * (wqq3_rk(k)*wqq2_rk(l)-wqq3_rk(l)*wqq2_rk(k))  &
                 + factor3 * wqq2_rkrl + factor4 * wqq3_rkrl  &
                 + factor1 * wqq_rk(k) * ( wqq3_rk(l) - wqq2_rk(l) )
     wqq_rkrl(l,k) = wqq_rkrl(k,l)

  end do
  end do

  end if

end subroutine F_qq_density_rhok_rhol




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine F_dq_density_rhok_rhol ( wdq_rkrl )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp,ncomp), INTENT(IN OUT)  :: wdq_rkrl

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k, l

  real(dp)                                   :: abbrev1, abbrev2
  real(dp)                                   :: factorA
  real(dp)                                   :: factor1, factor2, factor3, factor4
  real(dp)                                   :: dq_sum

  real(dp), dimension(ncomp)                 :: wdq2_z_k
  real(dp), dimension(ncomp)                 :: wdq3_z_k
  real(dp), dimension(ncomp)                 :: wdq_rk
  real(dp), dimension(ncomp)                 :: wdq2_rk
  real(dp), dimension(ncomp)                 :: wdq3_rk
  real(dp)                                   :: wdq2_rkrl, wdq3_rkrl
  !-----------------------------------------------------------------------------

  if ( abs( wdq2 ) > 1.E-50_dp ) then

  !-----------------------------------------------------------------------------
  ! some quantities independent of components
  !-----------------------------------------------------------------------------
  factor1 = 2.0_dp / ( wdq2 - wdq3 )
  factor2 = 0.5_dp * factor1 * factor1 * wdq2
  factor3 = 0.5_dp * factor2 * ( wdq2 - 2.0_dp * wdq3 )
  factor4 = 0.5_dp * factor2 * wdq2

  factorA = wdq2 / ( wdq2 - wdq3 )
  factorA = factorA * factorA

  !-----------------------------------------------------------------------------
  ! first derivatives to component-density rho_k
  !-----------------------------------------------------------------------------
  do k = 1, ncomp

     wdq2_rk(k) = wdq2_z_term * z3_rk(k)
     wdq3_rk(k) = wdq3_z_term * z3_rk(k)
     wdq2_z_k(k) = 0.0_dp
     wdq3_z_k(k) = 0.0_dp
     if ( abs( my_fac_dq(k) ) > machine_eps .OR. abs( Q_fac_dq(k) ) > machine_eps ) then
        do i = 1, ncomp
           if ( abs( my_fac_dq(i)*rhoi(i) ) > machine_eps .OR. abs( Q_fac_dq(i)*rhoi(i) ) > machine_eps ) then
              abbrev1 = rhoi(i) * ( pi_dq(i,k) + pi_dq(k,i) )
              wdq2_rk(k) = wdq2_rk(k) + abbrev1 * Idq2(i,k)
              wdq2_z_k(k) = wdq2_z_k(k) + abbrev1 * Idq2_z(i,k)
              do j = 1, ncomp
                 abbrev2 = rhoi(i) * rhoi(j) * ( psi_dq(i,j,k) + psi_dq(i,k,j) + psi_dq(k,i,j) )
                 wdq3_rk(k) = wdq3_rk(k) + abbrev2 * Idq3(i,j,k)
                 wdq3_z_k(k) = wdq3_z_k(k) + abbrev2 * Idq3_z(i,j,k)
              end do
           end if
        end do
     end if

     wdq_rk(k) = factor3 * wdq2_rk(k) + factorA * wdq3_rk(k)

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives to rho_k and to rho_l
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = k, ncomp

     dq_sum = 0.0_dp
     do i = 1, ncomp
        dq_sum = dq_sum + rhoi(i) *( psi_dq(i,k,l)+psi_dq(i,l,k)+psi_dq(k,i,l)  &
                                    +psi_dq(l,i,k)+psi_dq(k,l,i)+psi_dq(l,k,i) ) * Idq3(i,k,l)
     end do

     wdq2_rkrl = (pi_dq(k,l)+pi_dq(l,k)) * Idq2(k,l) + wdq2_z_k(k)*z3_rk(l) + wdq2_z_k(l)*z3_rk(k)  &
                 + wdq2_z2_term * z3_rk(k)*z3_rk(l)
     wdq3_rkrl = wdq3_z_k(k)*z3_rk(l) + wdq3_z_k(l)*z3_rk(k) + wdq3_z2_term *z3_rk(k)*z3_rk(l) + dq_sum

     wdq_rkrl(k,l) = factor1 * wdq2_rk(l)*wdq2_rk(k)  &
                 + factor2 * (wdq3_rk(k)*wdq2_rk(l)-wdq3_rk(l)*wdq2_rk(k))  &
                 + factor3 * wdq2_rkrl + factor4 * wdq3_rkrl  &
                 + factor1 * wdq_rk(k) * ( wdq3_rk(l) - wdq2_rk(l) )
     wdq_rkrl(l,k) = wdq_rkrl(k,l)

  end do
  end do

  end if

end subroutine F_dq_density_rhok_rhol


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fdd_rho_4
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fdd_rho_4 ( zdd, zdd_z, zdd_z2, zdd_z3 )

  use EOS_CONSTANTS, only: ddp2, ddp3, ddp4

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: zdd
  real(dp), intent(out)                      :: zdd_z
  real(dp), intent(out)                      :: zdd_z2
  real(dp), intent(out)                      :: zdd_z3
  !-----------------------------------------------------------------------------
  integer                                    :: i, j, l, m
  real(dp)                                   :: fdd_z, fdd_z2, fdd_z3, fdd_z4
  real(dp)                                   :: fdd2, fdd2_z, fdd2_z2, fdd2_z3, fdd2_z4
  real(dp)                                   :: fdd3, fdd3_z, fdd3_z2, fdd3_z3, fdd3_z4
  real(dp), dimension(ncomp,ncomp)           :: etaIdd2_z3
  real(dp), dimension(ncomp,ncomp)           :: etaIdd2_z4
  real(dp)                                   :: p_factor
  real(dp)                                   :: rho_2, diff_f, z3_fac
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Idd3_z3
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Idd3_z4
  real(dp), dimension(2:4)                   :: z3_m_t3, z3_m_t4
  real(dp), dimension(1:4)                   :: z3_m_tsq3, z3_m_tsq4
  !-----------------------------------------------------------------------------


  zdd = 0.0_dp
  zdd_z = 0.0_dp
  zdd_z2 = 0.0_dp
  zdd_z3 = 0.0_dp

  if ( abs( rdd2 ) > 1.E-50_dp ) then

     do m = 2, 4
        z3_m_t3(m) = REAL( (m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
        z3_m_t4(m) = REAL( (m+1)*m*(m-1)*(m-2), KIND=dp ) * z3**(m-3)
     end do
     do m = 1, 4
        z3_m_tsq3(m) = REAL( (m+2)*(m+1)*m, KIND=dp ) * z3**(m-1)
        z3_m_tsq4(m) = REAL( (m+2)*(m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
     end do

     etaIdd2_z3(:,:) = 0.0_dp
     etaIdd2_z4(:,:) = 0.0_dp
     eta2_Idd3_z3(:,:,:) = 0.0_dp
     eta2_Idd3_z4(:,:,:) = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j) ) > machine_eps ) then
                 etaIdd2_z3(i,j) = sum( ( ddp2(i,j,2:4) + eps_ij(i,j) * ddp4(i,j,2:4) ) * z3_m_t3(2:4) )
                 etaIdd2_z4(i,j) = sum( ( ddp2(i,j,3:4) + eps_ij(i,j) * ddp4(i,j,3:4) ) * z3_m_t4(3:4) )
                 do l = 1, ncomp
                    eta2_Idd3_z3(i,j,l) = sum( ddp3(i,j,l,1:4) * z3_m_tsq3(1:4) )
                    eta2_Idd3_z4(i,j,l) = sum( ddp3(i,j,l,2:4) * z3_m_tsq4(2:4) )
                 end do
              end if
           end do
        end if
     end do


     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fdd2 = rdd2 * rho
     fdd3 = rdd3 * rho_2
     fdd2_z = rdd2/z3t + rdd2_z_term * rho
     fdd3_z = z3_fac * rho * rdd3 + rdd3_z_term * rho_2
     fdd2_z2 = z3_fac * rdd2_z_term + rdd2_z2_term *rho
     fdd3_z2 = z3_fac /z3t  *rdd3 + 4.0_dp *rho *rdd3_z_term / z3t + rdd3_z2_term *rho_2
     diff_f = fdd2 - fdd3

     fdd2_z3 = 0.0_dp
     fdd2_z4 = 0.0_dp
     fdd3_z3 = 0.0_dp
     fdd3_z4 = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_factor(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( my_factor(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i) * x(j) * pi_dd(i,j)
                 fdd2_z3 = fdd2_z3 + p_factor * etaIdd2_z3(i,j)
                 fdd2_z4 = fdd2_z4 + p_factor * etaIdd2_z4(i,j)
                 do l = 1, ncomp
                    p_factor = x(i) * x(j) * x(l) * psi_dd(i,j,l)
                    fdd3_z3 = fdd3_z3 + p_factor * eta2_Idd3_z3(i,j,l)
                    fdd3_z4 = fdd3_z4 + p_factor * eta2_Idd3_z4(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     fdd2_z3 = fdd2_z3 / z3t
     fdd2_z4 = fdd2_z4 / z3t
     fdd3_z3 = fdd3_z3 / z3t**2
     fdd3_z4 = fdd3_z4 / z3t**2

     fdd_z = fdd2 * ( fdd2_z + ( fdd2*fdd3_z - fdd3*fdd2_z ) / diff_f ) / diff_f
     fdd_z2 = ( 2.0_dp*( fdd2_z**2 + fdd_z *(fdd3_z-fdd2_z) )  &
          + fdd2 *( fdd2 *( fdd2_z2 + fdd3_z2 ) - 2.0_dp*fdd3*fdd2_z2 ) / diff_f ) / diff_f
     fdd_z3=(  (6.0_dp*fdd2_z*fdd2_z2 + fdd2*fdd2_z3 - 2.0_dp*fdd2_z3*fdd3 -2.0_dp*fdd2_z2*fdd3_z &
          + 2.0_dp*fdd2_z*fdd3_z2 + fdd2*fdd3_z3 ) / ( 1.0_dp-fdd3/fdd2 )  &
          + 2.0_dp*fdd2_z*( fdd2_z**2 - 3.0_dp*fdd2_z2*fdd3 - fdd2_z *fdd3_z) / diff_f ) / diff_f  &
          + 2.0_dp * ( 2.0_dp*fdd_z2*(fdd3_z-fdd2_z)  &
          +     fdd_z*(fdd3_z2-fdd2_z2) - fdd_z/diff_f*(fdd3_z-fdd2_z)**2 ) / diff_f
     fdd_z4=( 12.0_dp*fdd2_z**2 *fdd2_z2+6.0_dp*fdd2*fdd2_z2**2  &
          +8.0_dp*fdd2*fdd2_z*fdd2_z3+fdd2*fdd2*fdd2_z4-6.0_dp*fdd2_z2**2 *fdd3  &
          -12.0_dp*fdd2_z*fdd2_z2*fdd3_z -8.0_dp*fdd2_z*fdd2_z3*fdd3  &
          -2.0_dp*fdd2*fdd2_z4*fdd3-4.0_dp*fdd2*fdd2_z3*fdd3_z  &
          +4.0_dp*fdd2*fdd2_z*fdd3_z3+fdd2**2 *fdd3_z4 ) /diff_f**2  &
          + 6.0_dp/diff_f* ( fdd_z3*(fdd3_z-fdd2_z)  &
          -fdd_z2/diff_f*(fdd3_z-fdd2_z)**2  &
          - fdd_z/diff_f*(fdd3_z-fdd2_z)*(fdd3_z2-fdd2_z2)  &
          + fdd_z2*(fdd3_z2-fdd2_z2) +1.0_dp/3.0_dp*fdd_z*(fdd3_z3-fdd2_z3) )

     zdd = fdd_z * z3
     zdd_z = fdd_z2 * z3 + fdd_z
     zdd_z2 = fdd_z3 * z3 + 2.0_dp * fdd_z2
     zdd_z3 = fdd_z4 * z3 + 3.0_dp * fdd_z3
     
  end if

end subroutine fdd_rho_4



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fqq_rho_4
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fqq_rho_4 ( zqq, zqq_z, zqq_z2, zqq_z3 )

  use EOS_CONSTANTS, only: qqp2, qqp3, qqp4

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: zqq
  real(dp), intent(out)                      :: zqq_z
  real(dp), intent(out)                      :: zqq_z2
  real(dp), intent(out)                      :: zqq_z3
  !-----------------------------------------------------------------------------
  integer                                    :: i, j, l, m
  real(dp)                                   :: fqq_z, fqq_z2, fqq_z3, fqq_z4
  real(dp)                                   :: fqq2, fqq2_z, fqq2_z2, fqq2_z3, fqq2_z4
  real(dp)                                   :: fqq3, fqq3_z, fqq3_z2, fqq3_z3, fqq3_z4
  real(dp), dimension(ncomp,ncomp)           :: etaIqq2_z3
  real(dp), dimension(ncomp,ncomp)           :: etaIqq2_z4
  real(dp)                                   :: p_factor
  real(dp)                                   :: rho_2, diff_f, z3_fac
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Iqq3_z3
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Iqq3_z4
  real(dp), dimension(2:4)                   :: z3_m_t3, z3_m_t4
  real(dp), dimension(1:4)                   :: z3_m_tsq3, z3_m_tsq4
  !-----------------------------------------------------------------------------


  zqq = 0.0_dp
  zqq_z = 0.0_dp
  zqq_z2 = 0.0_dp
  zqq_z3 = 0.0_dp

  if ( abs( rqq2 ) > 1.E-50_dp ) then

     do m = 2, 4
        z3_m_t3(m) = REAL( (m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
        z3_m_t4(m) = REAL( (m+1)*m*(m-1)*(m-2), KIND=dp ) * z3**(m-3)
     end do
     do m = 1, 4
        z3_m_tsq3(m) = REAL( (m+2)*(m+1)*m, KIND=dp ) * z3**(m-1)
        z3_m_tsq4(m) = REAL( (m+2)*(m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
     end do

     etaIqq2_z3(:,:) = 0.0_dp
     etaIqq2_z4(:,:) = 0.0_dp
     eta2_Iqq3_z3(:,:,:) = 0.0_dp
     eta2_Iqq3_z4(:,:,:) = 0.0_dp
     do i = 1, ncomp
        if ( abs( Q_factor(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j) ) > machine_eps ) then
                 etaIqq2_z3(i,j) = sum( ( qqp2(i,j,2:4) + eps_ij(i,j) * qqp4(i,j,2:4) ) * z3_m_t3(2:4) )
                 etaIqq2_z4(i,j) = sum( ( qqp2(i,j,3:4) + eps_ij(i,j) * qqp4(i,j,3:4) ) * z3_m_t4(3:4) )
                 do l = 1, ncomp
                    eta2_Iqq3_z3(i,j,l) = sum( qqp3(i,j,l,1:4) * z3_m_tsq3(1:4) )
                    eta2_Iqq3_z4(i,j,l) = sum( qqp3(i,j,l,2:4) * z3_m_tsq4(2:4) )
                 end do
              end if
           end do
        end if
     end do

     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fqq2 = rqq2 * rho
     fqq3 = rqq3 * rho_2
     fqq2_z = rqq2/z3t + rqq2_z_term * rho
     fqq3_z = z3_fac * rho * rqq3 + rqq3_z_term * rho_2
     fqq2_z2 = z3_fac * rqq2_z_term + rqq2_z2_term *rho
     fqq3_z2 = z3_fac /z3t  *rqq3 + 4.0_dp *rho *rqq3_z_term / z3t + rqq3_z2_term *rho_2
     diff_f = fqq2 - fqq3

     fqq2_z3 = 0.0_dp
     fqq2_z4 = 0.0_dp
     fqq3_z3 = 0.0_dp
     fqq3_z4 = 0.0_dp
     do i = 1, ncomp
        if ( abs( Q_factor(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_factor(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i) * x(j) * pi_qq(i,j)
                 fqq2_z3 = fqq2_z3 + p_factor * etaIqq2_z3(i,j)
                 fqq2_z4 = fqq2_z4 + p_factor * etaIqq2_z4(i,j)
                 do l = 1, ncomp
                    p_factor = x(i) * x(j) * x(l) * psi_qq(i,j,l)
                    fqq3_z3 = fqq3_z3 + p_factor * eta2_Iqq3_z3(i,j,l)
                    fqq3_z4 = fqq3_z4 + p_factor * eta2_Iqq3_z4(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     fqq2_z3 = fqq2_z3 / z3t
     fqq2_z4 = fqq2_z4 / z3t
     fqq3_z3 = fqq3_z3 / z3t**2
     fqq3_z4 = fqq3_z4 / z3t**2

     fqq_z = fqq2 * ( fqq2_z + ( fqq2*fqq3_z - fqq3*fqq2_z ) / diff_f ) / diff_f
     fqq_z2 = ( 2.0_dp*( fqq2_z**2 + fqq_z *(fqq3_z-fqq2_z) )  &
          + fqq2 *( fqq2 *( fqq2_z2 + fqq3_z2 ) - 2.0_dp*fqq3*fqq2_z2 ) / diff_f ) / diff_f
     fqq_z3=(2.0_dp*fqq2_z**3  +6.0_dp*fqq2*fqq2_z*fqq2_z2+fqq2*fqq2*fqq2_z3  &
          -6.0_dp*fqq2_z*fqq2_z2*fqq3-2.0_dp*fqq2_z**2 *fqq3_z  &
          -2.0_dp*fqq2*fqq2_z3*fqq3 -2.0_dp*fqq2*fqq2_z2*fqq3_z  &
          +2.0_dp*fqq2*fqq2_z*fqq3_z2+fqq2*fqq2*fqq3_z3) / diff_f**2  &
          + 2.0_dp/diff_f* ( 2.0_dp*fqq_z2*(fqq3_z-fqq2_z)  &
          +     fqq_z*(fqq3_z2-fqq2_z2)  &
          -     fqq_z/diff_f*(fqq3_z-fqq2_z)**2 )
     fqq_z4=( 12.0_dp*fqq2_z**2 *fqq2_z2+6.0_dp*fqq2*fqq2_z2**2  &
          +8.0_dp*fqq2*fqq2_z*fqq2_z3+fqq2*fqq2*fqq2_z4-6.0_dp*fqq2_z2**2 *fqq3  &
          -12.0_dp*fqq2_z*fqq2_z2*fqq3_z -8.0_dp*fqq2_z*fqq2_z3*fqq3  &
          -2.0_dp*fqq2*fqq2_z4*fqq3-4.0_dp*fqq2*fqq2_z3*fqq3_z  &
          +4.0_dp*fqq2*fqq2_z*fqq3_z3+fqq2**2 *fqq3_z4 ) /diff_f**2  &
          + 6.0_dp/diff_f* ( fqq_z3*(fqq3_z-fqq2_z)  &
          -fqq_z2/diff_f*(fqq3_z-fqq2_z)**2  &
          - fqq_z/diff_f*(fqq3_z-fqq2_z)*(fqq3_z2-fqq2_z2)  &
          + fqq_z2*(fqq3_z2-fqq2_z2) +1.0_dp/3.0_dp*fqq_z*(fqq3_z3-fqq2_z3) )

     zqq = fqq_z * z3
     zqq_z = fqq_z2 * z3 + fqq_z
     zqq_z2 = fqq_z3 * z3 + 2.0_dp * fqq_z2
     zqq_z3 = fqq_z4 * z3 + 3.0_dp * fqq_z3

  end if

end subroutine fqq_rho_4



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fdq_rho_4
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fdq_rho_4 ( zdq, zdq_z, zdq_z2, zdq_z3 )

  use EOS_CONSTANTS, only: dqp2, dqp3, dqp4

  !-----------------------------------------------------------------------------
  real(dp), intent(out)                      :: zdq
  real(dp), intent(out)                      :: zdq_z
  real(dp), intent(out)                      :: zdq_z2
  real(dp), intent(out)                      :: zdq_z3
  !-----------------------------------------------------------------------------
  integer                                    :: i, j, l, m
  real(dp)                                   :: fdq_z, fdq_z2, fdq_z3, fdq_z4
  real(dp)                                   :: fdq2, fdq2_z, fdq2_z2, fdq2_z3, fdq2_z4
  real(dp)                                   :: fdq3, fdq3_z, fdq3_z2, fdq3_z3, fdq3_z4
  real(dp), dimension(ncomp,ncomp)           :: etaIdq2_z3
  real(dp), dimension(ncomp,ncomp)           :: etaIdq2_z4
  real(dp)                                   :: p_factor
  real(dp)                                   :: rho_2, diff_f, z3_fac
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Idq3_z3
  real(dp), dimension(ncomp,ncomp,ncomp)     :: eta2_Idq3_z4
  real(dp), dimension(2:4)                   :: z3_m_t3, z3_m_t4
  real(dp), dimension(1:4)                   :: z3_m_tsq3, z3_m_tsq4
  !-----------------------------------------------------------------------------


  zdq = 0.0_dp
  zdq_z = 0.0_dp
  zdq_z2 = 0.0_dp
  zdq_z3 = 0.0_dp

  if ( abs( rdq2 ) > 1.E-50_dp ) then

     do m = 2, 4
        z3_m_t3(m) = REAL( (m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
        z3_m_t4(m) = REAL( (m+1)*m*(m-1)*(m-2), KIND=dp ) * z3**(m-3)
     end do
     do m = 1, 4
        z3_m_tsq3(m) = REAL( (m+2)*(m+1)*m, KIND=dp ) * z3**(m-1)
        z3_m_tsq4(m) = REAL( (m+2)*(m+1)*m*(m-1), KIND=dp ) * z3**(m-2)
     end do

     etaIdq2_z3(:,:) = 0.0_dp
     etaIdq2_z4(:,:) = 0.0_dp
     eta2_Idq3_z3(:,:,:) = 0.0_dp
     eta2_Idq3_z4(:,:,:) = 0.0_dp
     do i = 1, ncomp
        do j = 1, ncomp
           etaIdq2_z3(i,j) = sum( ( dqp2(i,j,2:4) + eps_ij(i,j) * dqp4(i,j,2:4) ) * z3_m_t3(2:4) )
           etaIdq2_z4(i,j) = sum( ( dqp2(i,j,3:4) + eps_ij(i,j) * dqp4(i,j,3:4) ) * z3_m_t4(3:4) )
           do l = 1, ncomp
              eta2_Idq3_z3(i,j,l) = sum( dqp3(i,j,l,1:4) * z3_m_tsq3(1:4) )
              eta2_Idq3_z4(i,j,l) = sum( dqp3(i,j,l,2:4) * z3_m_tsq4(2:4) )
           end do
        end do
     end do

     rho_2 = rho*rho
     z3_fac = 2.0_dp / z3t

     fdq2 = rdq2 * rho
     fdq3 = rdq3 * rho_2
     fdq2_z = rdq2/z3t + rdq2_z_term * rho
     fdq3_z = z3_fac * rho * rdq3 + rdq3_z_term * rho_2
     fdq2_z2 = z3_fac * rdq2_z_term + rdq2_z2_term *rho
     fdq3_z2 = z3_fac /z3t  *rdq3 + 4.0_dp *rho *rdq3_z_term / z3t + rdq3_z2_term *rho_2
     diff_f = fdq2 - fdq3

     fdq2_z3 = 0.0_dp
     fdq2_z4 = 0.0_dp
     fdq3_z3 = 0.0_dp
     fdq3_z4 = 0.0_dp
     do i = 1, ncomp
        if ( abs( my_fac_dq(i)*x(i) ) > machine_eps ) then
           do j = 1, ncomp
              if ( abs( Q_fac_dq(j)*x(j) ) > machine_eps ) then
                 p_factor = x(i) * x(j) * pi_dq(i,j)
                 fdq2_z3 = fdq2_z3 + p_factor * etaIdq2_z3(i,j)
                 fdq2_z4 = fdq2_z4 + p_factor * etaIdq2_z4(i,j)
                 do l = 1, ncomp
                    p_factor = x(i) * x(j) * x(l) * psi_dq(i,j,l)
                    fdq3_z3 = fdq3_z3 + p_factor * eta2_Idq3_z3(i,j,l)
                    fdq3_z4 = fdq3_z4 + p_factor * eta2_Idq3_z4(i,j,l)
                 end do
              end if
           end do
        end if
     end do
     fdq2_z3 = fdq2_z3 / z3t
     fdq2_z4 = fdq2_z4 / z3t
     fdq3_z3 = fdq3_z3 / z3t**2
     fdq3_z4 = fdq3_z4 / z3t**2

     fdq_z = fdq2 * ( fdq2_z + ( fdq2*fdq3_z - fdq3*fdq2_z ) / diff_f ) / diff_f
     fdq_z2 = ( 2.0_dp*( fdq2_z**2 + fdq_z *(fdq3_z-fdq2_z) )  &
          + fdq2 *( fdq2 *( fdq2_z2 + fdq3_z2 ) - 2.0_dp*fdq3*fdq2_z2 ) / diff_f ) / diff_f
     fdq_z3=(2.0_dp*fdq2_z**3  +6.0_dp*fdq2*fdq2_z*fdq2_z2+fdq2*fdq2*fdq2_z3  &
          -6.0_dp*fdq2_z*fdq2_z2*fdq3-2.0_dp*fdq2_z**2 *fdq3_z  &
          -2.0_dp*fdq2*fdq2_z3*fdq3 -2.0_dp*fdq2*fdq2_z2*fdq3_z  &
          +2.0_dp*fdq2*fdq2_z*fdq3_z2+fdq2*fdq2*fdq3_z3) / diff_f**2  &
          + 2.0_dp/diff_f* ( 2.0_dp*fdq_z2*(fdq3_z-fdq2_z)  &
          +     fdq_z*(fdq3_z2-fdq2_z2)  &
          -     fdq_z/diff_f*(fdq3_z-fdq2_z)**2 )
     fdq_z4=( 12.0_dp*fdq2_z**2 *fdq2_z2+6.0_dp*fdq2*fdq2_z2**2  &
          +8.0_dp*fdq2*fdq2_z*fdq2_z3+fdq2*fdq2*fdq2_z4-6.0_dp*fdq2_z2**2 *fdq3  &
          -12.0_dp*fdq2_z*fdq2_z2*fdq3_z -8.0_dp*fdq2_z*fdq2_z3*fdq3  &
          -2.0_dp*fdq2*fdq2_z4*fdq3-4.0_dp*fdq2*fdq2_z3*fdq3_z  &
          +4.0_dp*fdq2*fdq2_z*fdq3_z3+fdq2**2 *fdq3_z4 ) /diff_f**2  &
          + 6.0_dp/diff_f* ( fdq_z3*(fdq3_z-fdq2_z)  &
          -fdq_z2/diff_f*(fdq3_z-fdq2_z)**2  &
          - fdq_z/diff_f*(fdq3_z-fdq2_z)*(fdq3_z2-fdq2_z2)  &
          + fdq_z2*(fdq3_z2-fdq2_z2) +1.0_dp/3.0_dp*fdq_z*(fdq3_z3-fdq2_z3) )

     zdq = fdq_z * z3
     zdq_z = fdq_z2 * z3 + fdq_z
     zdq_z2 = fdq_z3 * z3 + 2.0_dp * fdq_z2
     zdq_z3 = fdq_z4 * z3 + 3.0_dp * fdq_z3

  end if

end subroutine fdq_rho_4



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine allocate_eos_quantities
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine allocate_eos_quantities

  use PARAMETERS, only: nsite
  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: compna, CAS_name, CAS_No, mm,  &
       sigma, epsilon_k, dipole_moment, quadru_moment, nhb_typ, nhb_no, eps_hb,  &
       kap_hb, parameter_set, assoc_scheme, kij, lij, kij_assoc, kij_t, cp_coefficients
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! species-names and composition of considered mixture
  !-----------------------------------------------------------------------------

  if(allocated(compna) ) return
  allocate( compna(ncomp) )
  allocate ( CAS_name( ncomp ) )
  allocate ( CAS_No( ncomp ) )

  allocate ( x( ncomp ) )

  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------

  allocate ( mm( ncomp ) )

  allocate ( mseg( ncomp ) )
  allocate ( sigma( ncomp ) )
  allocate ( epsilon_k( ncomp ) )
  allocate ( dipole_moment( ncomp ) )
  allocate ( quadru_moment( ncomp ) )
  
  allocate ( nhb_typ( ncomp ) )
  allocate ( nhb_no( ncomp, nsite ) )
  allocate ( eps_hb( ncomp, ncomp, nsite, nsite ) )
  allocate ( kap_hb( ncomp, ncomp ) )

  allocate ( parameter_set( ncomp ) )
  allocate ( assoc_scheme( ncomp ) )

  allocate ( cp_coefficients( ncomp, 4 ) )

  !-----------------------------------------------------------------------------
  ! PCP-SAFT eos polar constants
  !-----------------------------------------------------------------------------

  allocate ( qqp2( ncomp, ncomp, 0:4 ) )
  allocate ( qqp4( ncomp, ncomp, 0:4 ) )
  allocate ( ddp2( ncomp, ncomp, 0:4 ) )
  allocate ( ddp4( ncomp, ncomp, 0:4 ) )
  allocate ( dqp2( ncomp, ncomp, 0:4 ) )
  allocate ( dqp4( ncomp, ncomp, 0:4 ) )
  allocate ( qqp3( ncomp, ncomp, ncomp, 0:4 ) )
  allocate ( ddp3( ncomp, ncomp, ncomp, 0:4 ) )
  allocate ( dqp3( ncomp, ncomp, ncomp, 0:4 ) )

  !-----------------------------------------------------------------------------
  ! eos pure component and mixture parameters or parameter combinations
  !-----------------------------------------------------------------------------

  allocate ( dhs( ncomp ) )
  allocate ( uij( ncomp, ncomp ) )
  allocate ( uij_t( ncomp, ncomp ) )
  allocate ( sig_ij( ncomp, ncomp ) )
  allocate ( sig3_ij( ncomp, ncomp ) )
  allocate ( eps_ij( ncomp, ncomp ) )

  allocate ( kij( ncomp, ncomp ) )
  allocate ( lij( ncomp, ncomp ) )
  allocate ( kij_assoc( ncomp, ncomp ) )
  allocate ( kij_t( ncomp, ncomp ) )

  allocate ( ord1_lij_aux( ncomp ) )
  allocate ( ord1_lij( ncomp ) )
  allocate ( m_sig_eps( ncomp, ncomp ) )

  allocate ( ass_d( ncomp, ncomp, nsite, nsite ) )
  allocate ( ass_d_dt( ncomp, ncomp, nsite, nsite ) )
  allocate ( ass_d_dt2( ncomp, ncomp, nsite, nsite ) )
  allocate ( dij_ab( ncomp, ncomp ) )

  allocate ( pi_dd( ncomp, ncomp ) )
  allocate ( psi_dd( ncomp, ncomp, ncomp ) )
  allocate ( pi_qq( ncomp, ncomp ) )
  allocate ( psi_qq( ncomp, ncomp, ncomp ) )
  allocate ( pi_dq( ncomp, ncomp ) )
  allocate ( psi_dq( ncomp, ncomp, ncomp ) )

  allocate ( pi_dd_no_T( ncomp, ncomp ) )
  allocate ( psi_dd_no_T( ncomp, ncomp, ncomp ) )
  allocate ( pi_qq_no_T( ncomp, ncomp ) )
  allocate ( psi_qq_no_T( ncomp, ncomp, ncomp ) )
  allocate ( pi_dq_no_T( ncomp, ncomp ) )
  allocate ( psi_dq_no_T( ncomp, ncomp, ncomp ) )

  allocate ( my_factor( ncomp ) )
  allocate ( Q_factor( ncomp ) )
  allocate ( my_fac_dq( ncomp ) )
  allocate ( Q_fac_dq( ncomp ) )

  allocate ( my_factor_no_T( ncomp ) )
  allocate ( Q_factor_no_T( ncomp ) )
  allocate ( my_fac_dq_no_T( ncomp ) )
  allocate ( Q_fac_dq_no_T( ncomp ) )

  !-----------------------------------------------------------------------------
  ! density and composition dependent quantities
  !-----------------------------------------------------------------------------

  allocate ( gij( ncomp, ncomp ) )
  allocate ( mx( ncomp, nsite ) )

  allocate ( rhoi( ncomp ) )

  allocate ( z0_rk( ncomp ) )
  allocate ( z1_rk( ncomp ) )
  allocate ( z2_rk( ncomp ) )
  allocate ( z3_rk( ncomp ) )

  allocate ( Idd2( ncomp, ncomp ) )
  allocate ( Idd2_z( ncomp, ncomp ) )
  allocate ( Idd2_z2( ncomp, ncomp ) )
  allocate ( Iqq2( ncomp, ncomp ) )
  allocate ( Iqq2_z( ncomp, ncomp ) )
  allocate ( Iqq2_z2( ncomp, ncomp ) )
  allocate ( Idq2( ncomp, ncomp ) )
  allocate ( Idq2_z( ncomp, ncomp ) )
  allocate ( Idq2_z2( ncomp, ncomp ) )

  allocate ( Idd3( ncomp, ncomp, ncomp ) )
  allocate ( Idd3_z( ncomp, ncomp, ncomp ) )
  allocate ( Idd3_z2( ncomp, ncomp, ncomp ) )
  allocate ( Iqq3( ncomp, ncomp, ncomp ) )
  allocate ( Iqq3_z( ncomp, ncomp, ncomp ) )
  allocate ( Iqq3_z2( ncomp, ncomp, ncomp ) )
  allocate ( Idq3( ncomp, ncomp, ncomp ) )
  allocate ( Idq3_z( ncomp, ncomp, ncomp ) )
  allocate ( Idq3_z2( ncomp, ncomp, ncomp ) )

end subroutine allocate_eos_quantities


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine deallocate_eos_quantities
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deallocate_eos_quantities

  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: compna, CAS_name, CAS_No, mm,  &
       sigma, epsilon_k, dipole_moment, quadru_moment, nhb_typ, nhb_no, eps_hb,  &
       kap_hb, parameter_set, assoc_scheme, kij, lij, kij_assoc, kij_t, cp_coefficients
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! species-names and composition of considered mixture
  !-----------------------------------------------------------------------------

  deallocate ( compna )
  deallocate ( CAS_name )
  deallocate ( CAS_No )

  deallocate ( x )

  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------

  deallocate ( mm )

  deallocate ( mseg )
  deallocate ( sigma )
  deallocate ( epsilon_k )
  deallocate ( dipole_moment )
  deallocate ( quadru_moment )
  
  deallocate ( nhb_typ )
  deallocate ( nhb_no )
  deallocate ( eps_hb )
  deallocate ( kap_hb )

  deallocate ( parameter_set )
  deallocate ( assoc_scheme )

  deallocate ( cp_coefficients )

  !-----------------------------------------------------------------------------
  ! PCP-SAFT eos polar constants
  !-----------------------------------------------------------------------------

  deallocate ( qqp2 )
  deallocate ( qqp4 )
  deallocate ( ddp2 )
  deallocate ( ddp4 )
  deallocate ( dqp2 )
  deallocate ( dqp4 )
  deallocate ( qqp3 )
  deallocate ( ddp3 )
  deallocate ( dqp3 )

  !-----------------------------------------------------------------------------
  ! eos pure component and mixture parameters or parameter combinations
  !-----------------------------------------------------------------------------

  deallocate ( dhs )
  deallocate ( uij )
  deallocate ( uij_t )
  deallocate ( sig_ij )
  deallocate ( sig3_ij )
  deallocate ( eps_ij )

  deallocate ( kij )
  deallocate ( lij )
  deallocate ( kij_assoc )
  deallocate ( kij_t )

  deallocate ( ord1_lij_aux )
  deallocate ( ord1_lij )
  deallocate ( m_sig_eps )

  deallocate ( ass_d )
  deallocate ( ass_d_dt )
  deallocate ( ass_d_dt2 )
  deallocate ( dij_ab )

  deallocate ( pi_dd )
  deallocate ( psi_dd )
  deallocate ( pi_qq )
  deallocate ( psi_qq )
  deallocate ( pi_dq )
  deallocate ( psi_dq )

  deallocate ( pi_dd_no_T )
  deallocate ( psi_dd_no_T )
  deallocate ( pi_qq_no_T )
  deallocate ( psi_qq_no_T )
  deallocate ( pi_dq_no_T )
  deallocate ( psi_dq_no_T )

  deallocate ( my_factor )
  deallocate ( Q_factor )
  deallocate ( my_fac_dq )
  deallocate ( Q_fac_dq )

  deallocate ( my_factor_no_T )
  deallocate ( Q_factor_no_T )
  deallocate ( my_fac_dq_no_T )
  deallocate ( Q_fac_dq_no_T )
  
  !-----------------------------------------------------------------------------
  ! density and composition dependent quantities
  !-----------------------------------------------------------------------------

  deallocate ( gij )
  deallocate ( mx )

  deallocate ( rhoi )

  deallocate ( z0_rk )
  deallocate ( z1_rk )
  deallocate ( z2_rk )
  deallocate ( z3_rk )

  deallocate ( Idd2 )
  deallocate ( Idd2_z )
  deallocate ( Idd2_z2 )
  deallocate ( Iqq2 )
  deallocate ( Iqq2_z )
  deallocate ( Iqq2_z2 )
  deallocate ( Idq2 )
  deallocate ( Idq2_z )
  deallocate ( Idq2_z2 )

  deallocate ( Idd3 )
  deallocate ( Idd3_z )
  deallocate ( Idd3_z2 )
  deallocate ( Iqq3 )
  deallocate ( Iqq3_z )
  deallocate ( Iqq3_z2 )
  deallocate ( Idq3 )
  deallocate ( Idq3_z )
  deallocate ( Idq3_z2 )

end subroutine deallocate_eos_quantities

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine initialize_eos_quantities
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine initialize_eos_quantities

  use EOS_CONSTANTS
  use pcsaft_pure_and_binary_parameters, only: compna, mm, sigma, epsilon_k, dipole_moment,  &
       quadru_moment, nhb_typ, nhb_no, eps_hb, kap_hb, parameter_set, assoc_scheme,  &
       kij, lij, kij_assoc, kij_t, kij_T_dependent, lij_correction, cp_coefficients

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! species-names and composition of considered mixture
  !-----------------------------------------------------------------------------

  compna(:) = ''
  x(:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------

  mm(:) = 0.0_dp

  mseg(:) = 0.0_dp
  sigma(:) = 0.0_dp
  epsilon_k(:) = 0.0_dp
  dipole_moment(:) = 0.0_dp
  quadru_moment(:) = 0.0_dp
  
  nhb_typ(:) = 0
  nhb_no(:,:) = 0.0_dp
  eps_hb(:,:,:,:) = 0.0_dp
  kap_hb(:,:) = 0.0_dp

  parameter_set(:) = ''
  assoc_scheme(:) = ''

  cp_coefficients(:,:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! PCP-SAFT eos polar constants
  !-----------------------------------------------------------------------------

  qqp2(:,:,:) = 0.0_dp
  qqp4(:,:,:) = 0.0_dp
  ddp2(:,:,:) = 0.0_dp
  ddp4(:,:,:) = 0.0_dp
  dqp2(:,:,:) = 0.0_dp
  dqp4(:,:,:) = 0.0_dp
  qqp3(:,:,:,:) = 0.0_dp
  ddp3(:,:,:,:) = 0.0_dp
  dqp3(:,:,:,:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! eos pure component and mixture parameters or parameter combinations
  !-----------------------------------------------------------------------------

  dhs(:) = 0.0_dp
  uij(:,:) = 0.0_dp
  uij_t(:,:) = 0.0_dp
  sig_ij(:,:) = 0.0_dp
  sig3_ij(:,:) = 0.0_dp
  eps_ij(:,:) = 0.0_dp

  kij(:,:) = 0.0_dp
  lij(:,:) = 0.0_dp
  kij_assoc(:,:) = 0.0_dp
  kij_t(:,:) = 0.0_dp

  ord1_lij_aux(:) = 0.0_dp
  ord1_lij(:) = 0.0_dp
  m_sig_eps(:,:) = 0.0_dp

  ass_d(:,:,:,:) = 0.0_dp
  ass_d_dt(:,:,:,:) = 0.0_dp
  ass_d_dt2(:,:,:,:) = 0.0_dp
  dij_ab(:,:) = 0.0_dp

  assoc = .false.
  qudpole = .false.
  dipole = .false.
  dipole_quad = .false.

  lij_correction = .false.
  kij_T_dependent = .false.

  pi_dd(:,:) = 0.0_dp
  psi_dd(:,:,:) = 0.0_dp
  pi_qq(:,:) = 0.0_dp
  psi_qq(:,:,:) = 0.0_dp
  pi_dq(:,:) = 0.0_dp
  psi_dq(:,:,:) = 0.0_dp

  pi_dd_no_T(:,:) = 0.0_dp
  psi_dd_no_T(:,:,:) = 0.0_dp
  pi_qq_no_T(:,:) = 0.0_dp
  psi_qq_no_T(:,:,:) = 0.0_dp
  pi_dq_no_T(:,:) = 0.0_dp
  psi_dq_no_T(:,:,:) = 0.0_dp

  my_factor(:) = 0.0_dp
  Q_factor(:) = 0.0_dp
  my_fac_dq(:) = 0.0_dp
  Q_fac_dq(:) = 0.0_dp

  my_factor_no_T(:) = 0.0_dp
  Q_factor_no_T(:) = 0.0_dp
  my_fac_dq_no_T(:) = 0.0_dp
  Q_fac_dq_no_T(:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! density and composition dependent quantities
  !-----------------------------------------------------------------------------

  gij(:,:) = 0.0_dp
  mx(:,:) = 0.0_dp

  rhoi(:) = 0.0_dp

  z0_rk(:) = 0.0_dp
  z1_rk(:) = 0.0_dp
  z2_rk(:) = 0.0_dp
  z3_rk(:) = 0.0_dp

  Idd2(:,:) = 0.0_dp
  Idd2_z(:,:) = 0.0_dp
  Idd2_z2(:,:) = 0.0_dp
  Iqq2(:,:) = 0.0_dp
  Iqq2_z(:,:) = 0.0_dp
  Iqq2_z2(:,:) = 0.0_dp
  Idq2(:,:) = 0.0_dp
  Idq2_z(:,:) = 0.0_dp
  Idq2_z2(:,:) = 0.0_dp

  Idd3(:,:,:) = 0.0_dp
  Idd3_z(:,:,:) = 0.0_dp
  Idd3_z2(:,:,:) = 0.0_dp
  Iqq3(:,:,:) = 0.0_dp
  Iqq3_z(:,:,:) = 0.0_dp
  Iqq3_z2(:,:,:) = 0.0_dp
  Idq3(:,:,:) = 0.0_dp
  Idq3_z(:,:,:) = 0.0_dp
  Idq3_z2(:,:,:) = 0.0_dp

end subroutine initialize_eos_quantities

end module module_eos_derivatives

module properties

  use PARAMETERS, only: dp, NAV, RGAS
  use BASIC_VARIABLES, only: ncomp
  implicit none

  real(dp), public                                  :: pcalc_z_transfer
  real(dp), public                                  :: ztot_rho

  private
  public :: chemical_potential_tp, chemical_potential_trho,  &
       state_trho, state_tp, state_ps, state_ph, Helmholtz_density,  &
       chemical_potential_tp_derivative_tp, calculate_pressure, pressure,  &
       state_trhoJre,calculate_density

contains

subroutine Helmholtz_density ( tF, rhoi_in, f_dens )

  use module_eos_derivatives, only: F_density, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), intent(out)                             :: f_dens
  !-----------------------------------------------------------------------------

  real(dp), dimension(ncomp)                        :: xF
  real(dp)                                          :: rhoF
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  rhoF = sum( rhoi_in( 1:ncomp ) )
  if ( rhoF < 1.E-50_dp ) then
     write (*,*) 'F_density: density too low'
     stop
  end if
  xF( 1:ncomp ) = rhoi_in( 1:ncomp ) / rhoF

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! Helmholtz energy density
  !-----------------------------------------------------------------------------

  call F_density ( rhoi_in, f_dens )

end subroutine Helmholtz_density



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pressure
!
! calculate the pressure in unit (Pa) and derivatives to packing
! fraction eta (=z3). The first derivative to pressure is always delivered. The
! second and third derivatives are optional arguments.
!
! this subroutine should preferentially be called from subroutines only within
! this module. From outside, rather use SR "calculate_pressure"
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure ( eta_specified, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  use module_eos_derivatives, only: f_rho, f_rho_rho, f_rho_4, z3t, rho, kT
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_specified
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out), optional            :: pcalc_z
  real(dp), intent(out), optional            :: pcalc_z2
  real(dp), intent(out), optional            :: pcalc_z3
  !-----------------------------------------------------------------------------
  real(dp)                                   :: z_res, ztot, ztot_z
  !-----------------------------------------------------------------------------

  if ( .NOT.present( pcalc_z ) ) then

     call f_rho ( eta_specified, z_res )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT

  else if ( .NOT.present( pcalc_z2 ) ) then

     call f_rho_rho ( eta_specified, z_res, ztot_z )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT
     pcalc_z = ( ztot_z*rho + ztot/z3t ) * kT

     ztot_rho = ztot_z * z3t

  else

     call f_rho_4 ( eta_specified, z_res, ztot_z, pcalc_z2, pcalc_z3 )
     ztot = z_res + 1.0_dp
     pcalc = ztot * rho * kT
     pcalc_z = ( ztot_z*rho + ztot/z3t ) * kT

  end if

end subroutine pressure



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE calculate_pressure
!
! calculate the pressure in unit (Pa) and derivatives to packing
! fraction eta (=z3). The first derivative to pressure is always delivered. The
! second and third derivatives are optional arguments.
!
! this subroutine is called from code outside of this module
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine calculate_pressure ( tF, xF, eta_specified, rho_out, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  use module_eos_derivatives, only: rho, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  real(dp), dimension(ncomp), intent(in)     :: xF
  real(dp), intent(in)                       :: eta_specified
  real(dp), intent(out)                      :: rho_out
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out), optional            :: pcalc_z
  real(dp), intent(out), optional            :: pcalc_z2
  real(dp), intent(out), optional            :: pcalc_z3
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  call pressure ( eta_specified, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

  rho_out = rho

end subroutine calculate_pressure



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE calculate_density
!
! calculate the density (of a mixture) for defined (T,p,x). A starting value
! for the packing fraction (eta_start) needs to be given. A liquid-like value
! of the starting value leads to a liquid density (eta_out and rho_out), if a
! liquid root exists at the given conditions. For ! vapor-like starting values
! one gets a vapor density (if vapor root exists).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine calculate_density ( tF, pF, xF, eta_start, eta_out, rho_out )

  use module_eos_derivatives, only: rho, rho_independent_quantities
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: tF
  real(dp), intent(in)                       :: pF
  real(dp), dimension(ncomp), intent(in)     :: xF
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta_out
  real(dp), intent(out)                      :: rho_out

  real(dp)                                   :: eta
  real(dp)                                   :: pcalc, pcalc_z
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rho_out = rho
  eta_out = eta

end subroutine calculate_density



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_tp
!
!> \brief calculation of the chemical potentials to specified (t,p,x)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, p, x(ncomp), eta_start).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_tp ( tF, pF, xF, eta_start, rhoi_out, mu_res, mu, lnphi_nj )

  use module_eos_derivatives, only: kT, rhoi, rho, z3t, rho_independent_quantities,  &
       F_density_rhok, F_density_rhok_rho, F_density_rhok_rhol
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out), optional :: mu
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: lnphi_nj
  !-----------------------------------------------------------------------------

  integer                                           :: i
  real(dp)                                          :: eta
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: ztot
  real(dp), dimension(ncomp)                        :: w_rk
  real(dp), dimension(ncomp)                        :: w_rk_r
  real(dp), dimension(ncomp,ncomp)                  :: w_rkrl
  real(dp), dimension(ncomp,ncomp)                  :: wig_rkrl
  real(dp), dimension(ncomp)                        :: p_rhok
  real(dp), dimension(ncomp)                        :: partial_molar_v
  real(dp)                                          :: p_rho
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, pF=pF, xF=xF )
  
  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )	 !JRE pause here and enter to get sample calcs
  pcalc_z_transfer = pcalc_z

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------

  if ( present( lnphi_nj ) ) then

     call F_density_rhok_rho ( w_rk, w_rk_r )

     call F_density_rhok_rhol ( w_rkrl, wig_rkrl )

  else

     call F_density_rhok ( w_rk ) !JRE: w_rk = myres inside this function.

  end if

  if ( ztot > 0.0_dp ) mu_res( 1:ncomp ) = w_rk( 1:ncomp ) - LOG( ztot )      ! ln( phi ) for given ( T,p,x )

  !-----------------------------------------------------------------------------
  ! optional: calculate     mu_i(T,p,x)/kT = muRes+log(rho(i))
  !-----------------------------------------------------------------------------

  if ( present( mu ) ) then

     do i = 1, ncomp
        if ( rhoi(i) > 1.E-200_dp ) then
           mu(i) = w_rk(i) + log( rhoi(i) )
        else
           mu(i) = - 1.E200_dp
        end if
     end do

  end if


  !-----------------------------------------------------------------------------
  ! optional: calculate derivative matrix    n * d( ln(phi_i) ) / d( rho_j )
  !-----------------------------------------------------------------------------

  if ( present( lnphi_nj ) ) then

     p_rhok( 1:ncomp ) = rho * kT * w_rk_r( 1:ncomp ) + kT         ! in [Pa*Angstrom**3]
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     partial_molar_v( 1:ncomp ) = p_rhok( 1:ncomp ) / p_rho / rho  ! in [Angstrom**3]

     do i = 1, ncomp
        lnphi_nj( i,1:ncomp ) = rho * w_rkrl( i,1:ncomp )  &
             - partial_molar_v( i ) *partial_molar_v( 1:ncomp ) *rho**2 *p_rho / kT + 1.0_dp
     end do

  end if

end subroutine chemical_potential_tp


  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_tp_derivative_t
!
!> \brief calculation of the chemical potentials to specified (t,p,x)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, p, x(ncomp), eta_start).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_tp_derivative_tp ( tF, pF, xF, eta_start, rhoi_out, mu_res, lnphi_t, lnphi_p )

  use module_eos_derivatives, only: kT, rhoi, rho, z3t, rho_independent_quantities,  &
       F_density_rhok, F_density_rhok_rho, F_density_rhok_T, f_temp_rho
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out), optional :: lnphi_t
  real(dp), dimension(ncomp), intent(out), optional :: lnphi_p
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: ztot
  real(dp), dimension(ncomp)                        :: w_rk
  real(dp), dimension(ncomp)                        :: w_rk_r
  real(dp), dimension(ncomp)                        :: w_rk_t
  real(dp)                                          :: f_res, f_t, f_tr
  real(dp), dimension(ncomp)                        :: p_rhok
  real(dp), dimension(ncomp)                        :: partial_molar_v
  real(dp)                                          :: p_rho
  real(dp)                                          :: p_t
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, pF=pF, xF=xF )
  
  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------

  call F_density_rhok_rho ( w_rk, w_rk_r )

  mu_res( 1:ncomp ) = w_rk( 1:ncomp ) - LOG( ztot )      ! ln( phi ) for given ( T,p,x )

  p_rhok( 1:ncomp ) = rho * kT * w_rk_r( 1:ncomp ) + kT         ! in [Pa*Angstrom**3]
  p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
  partial_molar_v( 1:ncomp ) = p_rhok( 1:ncomp ) / p_rho / rho  ! in [Angstrom**3]

  if ( present( lnphi_p ) ) then

     lnphi_p( 1:ncomp ) = partial_molar_v( 1:ncomp) / kT -1.0_dp / pcalc ! in [1/Pa]

  else if ( present( lnphi_t ) ) then

     call F_density_rhok_T ( w_rk_t )

     call f_temp_rho ( f_res, f_t, f_tr )

     p_t = rho * rho * kT * f_tr + pcalc / tF

     lnphi_t( 1:ncomp ) = w_rk_t( 1:ncomp )  + 1.0_dp/tF - partial_molar_v( 1:ncomp ) * p_t / kT

  else

     write (*,*) ' use SR chemical_potential_tp_derivative_tp with at least one optional argument'
     stop

  end if

end subroutine chemical_potential_tp_derivative_tp


  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine chemical_potential_trho
!
!> \brief calculation of the chemical potentials to specified (t,rhoi)
!!
!! This subroutine gives the residual chemical potential and the chemical
!! potential:
!!     \f$ mu_i^{res}(T,p,x)/kT = ln( phi_i ) \f$
!! The required input (T, rhoi(ncomp) ).
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine chemical_potential_trho ( tF, rhoi_in, mu_res, mu, w_rkrl, wig_rkrl )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: rhoi, dhs, rho_independent_quantities,  &
       density_terms, Iterate_Association, F_density_rhok, F_density_rhok_rhol
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), dimension(ncomp), intent(out)           :: mu_res
  real(dp), dimension(ncomp), intent(out)           :: mu
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: w_rkrl
  real(dp), dimension(ncomp,ncomp), intent(out), optional :: wig_rkrl
  !-----------------------------------------------------------------------------

  integer                                           :: i
  real(dp)                                          :: eta
  real(dp), dimension(ncomp)                        :: xF
  real(dp), dimension(ncomp)                        :: w_rk
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! verify reasonable input
  !-----------------------------------------------------------------------------

  call check_input ( tF, rhoiF=rhoi_in )

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  rhoi( 1:ncomp ) = rhoi_in( 1:ncomp )
  xF( 1:ncomp ) = rhoi( 1:ncomp ) / sum( rhoi( 1:ncomp ) )

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! precalculate some density dependent terms
  !-----------------------------------------------------------------------------

  eta = PI_6 * sum( rhoi( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )

  call density_terms ( eta )

  !-----------------------------------------------------------------------------
  ! iterate fraction of unbonded association sites
  !-----------------------------------------------------------------------------
  call Iterate_Association

  !-----------------------------------------------------------------------------
  ! calculate     mu_i^res(T,p,x)/kT = ln( phi_i )
  !-----------------------------------------------------------------------------
  call F_density_rhok ( w_rk )

  mu_res( 1:ncomp ) = w_rk( 1:ncomp )             ! = ln( phi ) for given ( T,rhoi )

  do i = 1, ncomp
     if ( rhoi(i) > 1.E-200_dp ) then
        mu(i) = mu_res(i) + log( rhoi(i) )
     else
        mu(i) = - 1.E200_dp
     end if
  end do

  !-----------------------------------------------------------------------------
  ! if second derivative is required (optional argument)
  !-----------------------------------------------------------------------------

  if ( present( w_rkrl ) .AND. present( wig_rkrl ) ) then

     call F_density_rhok_rhol ( w_rkrl, wig_rkrl )

  end if

end subroutine chemical_potential_trho



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_tp
!
!> \brief molar enthalpy, entropy, rho for specified (T,p,x)
!!
!! This subroutine gives p, s, h, rho, cp (optional). The required input (T,p,x)
!! and starting value for density.
!! The residual entropy s_res for variables (T,rho,x)! or (T,rho_k)! is given as
!! an optional output. This quantity is used in entropy scaling.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_tp ( tF, pF, xF, eta_start, rhoi_out, molar_rho, s, h, s_res, cp )

  use PARAMETERS, only: NAV
  use module_eos_derivatives, only: rhoi, kT, rho, z3t, rho_independent_quantities,  &
       density_terms, f_temp, f_temp_rho, f_temp_temp
  !JRE use ideal_gas_enthalpy, only: enthalpy_ig 							 !JRE disabled because it would require adding another module.
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in)                              :: pF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: s_res
  real(dp), intent(out), optional                   :: cp
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: ztot
  real(dp)                                          :: pcalc, pcalc_z
  real(dp)                                          :: f_res, f_t, f_t2, f_tr
  real(dp)                                          :: sres
  real(dp)                                          :: p_rho, p_t
  real(dp)                                          :: cv_residual
  real(dp)                                          :: cmprsblty
  !real(dp)                                          :: hig, sig, cpig						!JRE disabled because it would require adding another module.
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call density_iteration ( eta_start, pF, eta, pcalc, pcalc_z )

  rhoi_out( 1:ncomp ) = rhoi( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate molar enthalpy
  !-----------------------------------------------------------------------------

  if ( .NOT.present( cp ) ) then

     call f_temp ( f_res, f_t )

  else

     !--------------------------------------------------------------------------
     ! partial derivatives needed for cp
     !--------------------------------------------------------------------------

     call f_temp_rho ( f_res, f_t, f_tr )

     ! ------ calculate pressure derivatives -----------------------------------
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     cmprsblty=p_rho/(RGAS*tF)
     p_t = rho * rho * kT * f_tr + pcalc / tF                      ! in [Pa/K]

     ! ------ derivative to temperature ----------------------------------------
     call f_temp_temp ( f_t, f_t2 )

     cv_residual = - ( tF * f_t2 + 2.0_dp * f_t ) * RGAS * tF      ! in [J/(mol*K)]
     cp = cv_residual + tF * ( p_t / rho )**2 / p_rho * NAV/1.E30_dp - RGAS

  end if

  !-----------------------------------------------------------------------------
  ! molar residual enthalpy and molar residual entropy
  !-----------------------------------------------------------------------------

  h = - tF * f_t + ( ztot - 1.0_dp )              ! h_res/kT, dimensionless
  h = h * RGAS * tF                               ! h_res in [J/mol]

  sres = - tF * f_t - f_res                       ! s_res/k for given (T,rho,x), dimensionless
  s = sres + log( ztot )                          ! s_res/k = s_res_molar/R for given (T,p,x), dimensionless
  s = s * RGAS                                    ! s_res in [J/(mol*K)]

  if ( present( s_res ) ) s_res = sres

  !-----------------------------------------------------------------------------
  ! add ideal gas contributions to molar enthalpy and molar entropy
  !-----------------------------------------------------------------------------

  !call enthalpy_ig ( tF, pcalc, xF, cpig, hig, sig )				!JRE disabled because it would require adding another module.
  !h = h + hig
  !s = s + sig
  !if ( present( cp ) ) cp = cp + cpig

  !-----------------------------------------------------------------------------
  ! molar density
  !-----------------------------------------------------------------------------

  molar_rho = rho / NAV * 1.E30_dp                                   ! in units [mol/m**3]
   
end subroutine state_tp



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_trho
!
!> \brief molar enthalpy, entropy, rho for specified (T, rhoi(nc))
!!
!! This subroutine gives p, s, h, rho, cp (optional). The required input (T, rhoi(nc)).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_trho ( tF, rhoi_in, pcalc, molar_rho, s, h, s_res, cp )

  use PARAMETERS, only: PI_6, NAV
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs, kT, rho, z3t, rho_independent_quantities,  &
       density_terms, f_temp, f_temp_rho, f_temp_temp
  !use ideal_gas_enthalpy, only: enthalpy_ig 	  !JRE disabled because it would require adding another module.
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), intent(out)                             :: pcalc
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: s_res
  real(dp), intent(out), optional                   :: cp
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: ztot
  real(dp)                                          :: pcalc_z
  real(dp)                                          :: f_res, f_t, f_t2, f_tr
  real(dp)                                          :: sres
  real(dp)                                          :: p_rho, p_t
  real(dp)                                          :: cv_residual
  real(dp), dimension(ncomp)                        :: xF
  !real(dp)                                          :: hig, sig, cpig	 !JRE disabled because it would require adding another module.
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  xF( 1:ncomp ) = rhoi_in( 1:ncomp ) / sum( rhoi_in(1:ncomp) )
  eta = PI_6 * sum( rhoi_in( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! calculate pressure for given (t,rhoi(:)) or actually for given (t,eta,x)
  !-----------------------------------------------------------------------------

  call pressure ( eta, pcalc, pcalc_z )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate molar enthalpy
  !-----------------------------------------------------------------------------

  if ( .NOT.present( cp ) ) then

     call f_temp ( f_res, f_t )

  else

     !--------------------------------------------------------------------------
     ! partial derivatives needed for cp
     !--------------------------------------------------------------------------

     call f_temp_rho ( f_res, f_t, f_tr )

     ! ------ calculate pressure derivatives -----------------------------------
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
     p_t = rho * rho * kT * f_tr + pcalc / tF                      ! in [Pa/K]

     ! ------ derivative to temperature ----------------------------------------
     call f_temp_temp ( f_t, f_t2 )

     cv_residual = - ( tF * f_t2 + 2.0_dp * f_t ) * RGAS * tF      ! in [J/(mol*K)]
     cp = cv_residual + tF * ( p_t / rho )**2 / p_rho * NAV/1.E30_dp - RGAS

  end if

  !-----------------------------------------------------------------------------
  ! molar residual enthalpy and molar residual entropy
  !-----------------------------------------------------------------------------

  h = - tF * f_t + ( ztot - 1.0_dp )              ! h_res/kT, dimensionless
  h = h * RGAS * tF                               ! h_res in [J/mol]

  sres = - tF * f_t - f_res                      ! s_res/k for given (T,rho,x), dimensionless
  s = sres + log( ztot )                         ! s_res/k = s_res_molar/R for given (T,p,x), dimensionless
  s = s * RGAS                                    ! s_res in [J/(mol*K)]

  if ( present( s_res ) ) s_res = sres

  !-----------------------------------------------------------------------------
  ! add ideal gas contributions to molar enthalpy and molar entropy
  !-----------------------------------------------------------------------------

  !call enthalpy_ig ( tF, pcalc, xF, cpig, hig, sig )	  !JRE disabled because it would require adding another module.
  !h = h + hig														  !JRE disabled because it would require adding another module.
  !s = s + sig														  !JRE disabled because it would require adding another module.
  !if ( present( cp ) ) cp = cp + cpig						   !JRE disabled because it would require adding another module.

  !-----------------------------------------------------------------------------
  ! molar density
  !-----------------------------------------------------------------------------

  molar_rho = rho / NAV * 1.E30_dp                                   ! in units [mol/m**3]
   
end subroutine state_trho

subroutine state_trhoJre ( tF, rhoi_in, pcalc, molar_rho, eta, ztot, sRes_R, hRes_RT, cmprsblty, cvRes_R, cpRes_R )

  use PARAMETERS, only: PI_6, NAV
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs, kT, rho, z3t, rho_independent_quantities,  &
       density_terms, f_temp, f_temp_rho, f_temp_temp
  !JRE use ideal_gas_enthalpy, only: enthalpy_ig 
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), dimension(ncomp), intent(in)            :: rhoi_in
  real(dp), intent(out)                             :: pcalc
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: sRes_R, hRes_RT, cmprsblty, cvRes_R, cpRes_R
  !-----------------------------------------------------------------------------

  real(dp)                                          :: eta
  real(dp)                                          :: ztot
  real(dp)                                          :: pcalc_z
  real(dp)                                          :: f_res, f_t, f_t2, f_tr
  real(dp)                                          :: p_rho, p_t
  real(dp), dimension(ncomp)                        :: xF
  !real(dp)                                          :: sres
  !real(dp)                                          :: hig, sig, cpig
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define variables used in EOS
  !-----------------------------------------------------------------------------

  xF( 1:ncomp ) = rhoi_in( 1:ncomp ) / sum( rhoi_in(1:ncomp) )
  eta = PI_6 * sum( rhoi_in( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, xF )

  !-----------------------------------------------------------------------------
  ! calculate pressure for given (t,rhoi(:)) or actually for given (t,eta,x)
  !-----------------------------------------------------------------------------

  call pressure ( eta, pcalc, pcalc_z )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  ztot = pcalc / ( kT * rho )

  !-----------------------------------------------------------------------------
  ! calculate molar enthalpy
  !-----------------------------------------------------------------------------

   call f_temp ( f_res, f_t )

     !--------------------------------------------------------------------------
     ! partial derivatives needed for cp
     !--------------------------------------------------------------------------

     call f_temp_rho ( f_res, f_t, f_tr )

  !-----------------------------------------------------------------------------
  ! molar density
  !-----------------------------------------------------------------------------

     molar_rho = rho / NAV * 1.E30_dp                                   ! in units [mol/m**3]

     ! ------ calculate pressure derivatives -----------------------------------
     p_rho = pcalc_z * z3t                                         ! in [Pa*Angstrom**3]
   
     cmprsblty=p_rho/(RGAS*tF)/(1.D30/NAV) !convert rho to molar_rho. Pa to MPa and mol/m^3 to mol/cc cancel.
     
     p_t = rho * rho * kT * f_tr + pcalc / tF                      ! in [Pa/K]

     ! ------ derivative to temperature ----------------------------------------
     call f_temp_temp ( f_t, f_t2 )

     cvRes_R = - ( tF * f_t2 + 2.0_dp * f_t ) * tF      ! in [J/(mol*K)]
     cpRes_R = cvRes_R + tF * ( p_t / rho )**2 / p_rho * NAV/1.E30_dp/RGAS - 1

  !-----------------------------------------------------------------------------
  ! molar residual enthalpy and molar residual entropy
  !-----------------------------------------------------------------------------

  hRes_RT = - tF * f_t + ( ztot - 1.0_dp )              ! h_res/kT, dimensionless
  !h = h * RGAS * tF                               ! h_res in [J/mol]

  sres_R = - tF * f_t - f_res                      ! s_res/k for given (T,rho,x), dimensionless
  !s = sres + log( ztot )                         ! s_res/k = s_res_molar/R for given (T,p,x), dimensionless
  !s = s * RGAS                                    ! s_res in [J/(mol*K)]


end subroutine state_trhoJre !JRE see state_trho above for the original, including implementation of Hig,Sig,Cpig.  



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_ps
!
! input, specifications:     pF, sF, xF(:)
! input, starting values:    eta_start, tcalc
! output:                    tcalc, rho, h,  cp (optional)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_ps ( pF, sF, xF, eta_start, tcalc, rhoi_out, molar_rho, h, cp )

  use BASIC_VARIABLES, only: ncomp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: pF
  real(dp), intent(in)                              :: sF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), intent(inout)                           :: tcalc
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: h
  real(dp), intent(out), optional                   :: cp

  !-----------------------------------------------------------------------------
  integer, parameter                                :: maxiter = 100
  integer                                           :: i_iter
  real(dp)                                          :: t_org
  real(dp)                                          :: f_error, dfdt
  real(dp), parameter                               :: f_tolerance = 1.E-10_dp
  real(dp), parameter                               :: delta_t = 1.E-4_dp
  real(dp)                                          :: sF1, sF2
  real(dp)                                          :: t_step
  logical                                           :: converged
  !-----------------------------------------------------------------------------

  i_iter = 0

  converged = .false.

  do while ( .NOT.converged .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     t_org = tcalc
     tcalc = t_org - delta_t
     call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF1, h )

     tcalc = t_org
     if ( present( cp ) ) then
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF2, h, cp=cp )
     else
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, sF2, h )
     end if

     f_error = ( sF2 - sF )
     dfdt = ( sF2 - sF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  50.0_dp, t_step )
     t_step = max( -50.0_dp, t_step )

     tcalc = tcalc - t_step

     IF ( ABS( f_error ) < f_tolerance ) converged = .true.

  end do

  if ( .NOT.converged ) write (*,*) 'ps_state not converged'

end subroutine state_ps



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_ph
!
! input, specifications:     pF, hF, xF(:)
! input, starting values:    eta_start, tcalc
! output:                    tcalc, rho, s,  cp (optional)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_ph ( pF, hF, xF, eta_start, tcalc, rhoi_out, molar_rho, s, cp )

  use BASIC_VARIABLES, only: ncomp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: pF
  real(dp), intent(in)                              :: hF
  real(dp), dimension(ncomp), intent(in)            :: xF
  real(dp), intent(in)                              :: eta_start
  real(dp), intent(inout)                           :: tcalc
  real(dp), dimension(ncomp), intent(out)           :: rhoi_out
  real(dp), intent(out)                             :: molar_rho
  real(dp), intent(out)                             :: s
  real(dp), intent(out), optional                   :: cp

  !-----------------------------------------------------------------------------
  integer, parameter                                :: maxiter = 100
  integer                                           :: i_iter
  real(dp)                                          :: t_org
  real(dp)                                          :: f_error, dfdt
  real(dp), parameter                               :: f_tolerance = 1.E-10_dp
  real(dp), parameter                               :: delta_t = 1.E-4_dp
  real(dp)                                          :: hF1, hF2
  real(dp)                                          :: t_step
  logical                                           :: converged
  !-----------------------------------------------------------------------------

  i_iter = 0

  converged = .false.

  do while ( .NOT.converged .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     t_org = tcalc
     tcalc = t_org - delta_t
     call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF1 )

     tcalc = t_org
     if ( present( cp ) ) then
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF2, cp=cp )
     else
        call state_tp ( tcalc, pF, xF, eta_start, rhoi_out, molar_rho, s, hF2 )
     end if

     f_error = ( hF2 - hF )
     dfdt = ( hF2 - hF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  50.0_dp, t_step )
     t_step = max( -50.0_dp, t_step )

     tcalc = tcalc - t_step

     IF ( ABS( f_error ) < f_tolerance ) converged = .true.

  end do

  if ( .NOT.converged ) write (*,*) 'ph_state not converged'

end subroutine state_ph



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE check_input
!
! verify that all input variables are physically reasonable, otherwise stop
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_input ( tF, pF, xF, rhoiF )

  use PARAMETERS, only: machine_eps
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                              :: tF
  real(dp), intent(in), optional                    :: pF
  real(dp), dimension(ncomp), intent(in), optional  :: xF
  real(dp), dimension(ncomp), intent(in), optional  :: rhoiF
  !-----------------------------------------------------------------------------
  integer                                           :: i
  real(dp)                                          :: sum_x
  logical                                           :: immediate_exit
  !-----------------------------------------------------------------------------

  immediate_exit = .false.

  if ( present( pF ) ) then

     !--------------------------------------------------------------------------
     ! consider case where (T,p,x) are specified variables
     !--------------------------------------------------------------------------

     !-----------------------------------------------------------------------------
     ! check for NaN and expected array-size
     !-----------------------------------------------------------------------------
     if ( tF /= tF .or. pF /= pF ) then
        write (*,*) ' input is "not a number" !'
        immediate_exit = .true.
     end if

     if ( tF < 0.0_dp .or. pF < 0.0_dp ) then
        write (*,*) ' input is not physical !'
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of pressure (pF)
     !--------------------------------------------------------------------------
     if ( pF < 1.E-50_dp ) then
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of mole fractions xF
     !--------------------------------------------------------------------------
     if ( size( xF, 1) /= ncomp ) then
        write (*,*) ' array-size of xF is not appropriate !'
        immediate_exit = .true.
     end if

     do i = 1, ncomp
        if ( xF(i) /= xF(i) ) then
           write (*,*) ' composition input x is "not a number" ! Species:',i
           immediate_exit = .true.
        end if
     end do

     if ( MINVAL( xF( 1:ncomp ) ) < 0.0_dp ) then
        write (*,*) ' mole fraction x is negative, of species', MINLOC( xF(1:ncomp) )
        immediate_exit = .true.
     end if

     sum_x = sum( xF( 1:ncomp ) )
     if ( abs( sum_x - 1.0_dp ) > 1.E4_dp * machine_eps ) then
        if ( abs( sum_x ) < machine_eps ) then
           write (*,*) ' mole fracitons are zero'
        else
           write (*,*) ' sum of mole fracitons is not equal one'
           ! xF( 1:ncomp ) = xF( 1:ncomp ) / sum_x    ! normalize mole fractions
        end if
        immediate_exit = .true.
     end if

  else

     !--------------------------------------------------------------------------
     ! consider case where (T,rhoi) are specified variables
     !--------------------------------------------------------------------------

     !--------------------------------------------------------------------------
     ! check for NaN and expected array-size
     !--------------------------------------------------------------------------
     if ( tF /= tF ) then
        write (*,*) ' input is "not a number" !'
        immediate_exit = .true.
     end if

     if ( tF < 0.0_dp ) then
        write (*,*) ' input is not physical !'
        immediate_exit = .true.
     end if

     !--------------------------------------------------------------------------
     ! verify proper specification of species densities rhoi
     !--------------------------------------------------------------------------
     if ( size( rhoiF, 1) /= ncomp ) then
        write (*,*) ' array-size of rhoiF is not appropriate !'
        immediate_exit = .true.
     end if

     do i = 1, ncomp
        if ( rhoiF(i) /= rhoiF(i) ) then
           write (*,*) ' species desnity rhoi is "not a number" ! Species:',i
           immediate_exit = .true.
        end if
     end do

     if ( MINVAL( rhoiF( 1:ncomp ) ) < 0.0_dp ) then
        write (*,*) ' species desnity rhoi is negative, of species', MINLOC( rhoiF(1:ncomp) )
        immediate_exit = .true.
     end if

     if ( sum( rhoiF(:) ) < 1.E-50_dp ) then
        write (*,*) ' overall density too low', sum( rhoiF(:) )
        immediate_exit = .true.
     end if

  end if

  !-----------------------------------------------------------------------------
  ! if unphysical input specifications: stop
  !-----------------------------------------------------------------------------
  if ( immediate_exit ) then
     write (*,*) 'check_input: unphysical input specs', tF, pF, xF(:)
     stop
  end if

end subroutine check_input

end module properties

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine DENSITY_ITERATION
!
!> \brief density iteration
!!
!! iterates the density until the calculated pressure 'pcalc' is equal to
!! the specified pressure 'p'. \n  A Newton-scheme is used for determining
!! the root to the objective function \f$ f(eta) = (p_{ges} / p ) - 1.0. \f$
!! The starting value for the packing fraction (= dimensionless density) is
!! eta_start.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine density_iteration ( eta_start, p, eta, pcalc, pcalc_z )

  use PARAMETERS, only: dp, machine_eps
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(in)                       :: p
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z

  !-----------------------------------------------------------------------------
  integer, parameter                         :: outp = 0
  integer, parameter                         :: max_iterations = 50
  real(dp), parameter                        :: density_tolerance = 1.E-12

  integer                                    :: i
  integer                                    :: count_num_noise
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: pcalc_z2, pcalc_z3
  real(dp)                                   :: error
  real(dp)                                   :: abs_error
  real(dp)                                   :: dydx
  real(dp)                                   :: delta_eta
  real(dp)                                   :: eta_old
  logical                                    :: write_output
  !-----------------------------------------------------------------------------

  eta_iteration = eta_start

  i = 0
  count_num_noise = 0
  eta_old = eta_start
  delta_eta = 0.0_dp

  write_output = .false.

  !-----------------------------------------------------------------------------
  ! iterate density until p_calc = p
  !-----------------------------------------------------------------------------

  abs_error = density_tolerance + 1.0_dp

  do while ( abs_error > density_tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z )

     error = ( pcalc / p ) - 1.0_dp

     if ( write_output .OR. outp > 1 ) write (*,'(i4,4G25.15)') i, error, error / pcalc_z*p, pcalc, eta_iteration

     !--------------------------------------------------------------------------
     ! correction for instable region
     !--------------------------------------------------------------------------
     if ( pcalc_z < 0.0_dp .AND. i < max_iterations ) then

        if ( .NOT.write_output .AND. outp == 1 ) write (*,'(i4,5G25.15)') i, error, error / pcalc_z*p,  &
                                                                          pcalc, eta_iteration, eta_start
        if ( outp == 1 ) write_output = .true.
        if ( i == 1 .AND. eta_old <= 0.25_dp ) delta_eta = 1.9_dp * eta_old
        if ( i == 1 .AND. eta_old >  0.25_dp ) delta_eta = - 0.08_dp
           
        eta_iteration = eta_old - 0.5_dp * delta_eta

        i = i + 1
        call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

        error = ( pcalc / p ) - 1.0_dp

        if ( write_output .OR. outp > 1 ) write (*,'(i4,5G25.15,a)') i, error, error / pcalc_z*p,  &
                                                                     pcalc, eta_iteration, eta_old, ' BBB'

        if ( eta_iteration > 0.55_dp ) then                               ! artificial density root
           call pressure_spinodal ( eta_iteration, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( eta_iteration > 0.55_dp ) then
              if ( error <  0.0_dp ) return                                  ! no liquid solution possible
              if ( error >= 0.0_dp ) eta_iteration = eta_iteration * 0.98_dp ! no solution found so far
           else
              if ( error >  0.0_dp ) eta_iteration = 0.001_dp                ! no liquid solution possible
              if ( error <= 0.0_dp ) eta_iteration = eta_iteration * 1.1_dp  ! no solution found so far
           end if
        else if ( error > 0.0_dp .AND. pcalc_z2 > 0.0_dp ) then           ! no liquid density found
           call pressure_spinodal ( eta_start, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( error >  0.0_dp ) eta_iteration = 0.001_dp                ! no liquid solution possible
           if ( error <= 0.0_dp ) eta_iteration = eta_iteration * 1.1_dp  ! no solution found so far
        else if ( error < 0.0_dp .AND. pcalc_z2 < 0.0_dp ) then           ! no vapor density found
           call pressure_spinodal ( eta_start, eta_iteration, pcalc, pcalc_z )
           error  = ( pcalc / p ) - 1.0_dp
           if ( error <  0.0_dp ) eta_iteration = 0.5_dp                  ! no vapor solution possible
           if ( error >= 0.0_dp ) eta_iteration = eta_iteration * 0.9_dp  ! no solution found so far
        else
           if ( i /= 2 ) then           ! first approach
              eta_iteration = (eta_iteration + eta_start) / 2.0_dp
              if ( abs( eta_iteration - eta_start) < machine_eps ) eta_iteration = eta_iteration + 0.2_dp
           end if
        end if

        CYCLE

     end if


     !--------------------------------------------------------------------------
     ! Newton iteration of density (eta)
     !--------------------------------------------------------------------------
     eta_old = eta_iteration

     dydx = pcalc_z / p
     delta_eta = error / dydx
     if ( abs( delta_eta ) >  0.05_dp ) delta_eta = sign( 0.05_dp, delta_eta )     

     eta_iteration = eta_iteration - delta_eta

     if (eta_iteration > 0.9_dp)  eta_iteration = 0.6_dp
     if (eta_iteration <= 0.0_dp) eta_iteration = 1.E-15_dp

     !--------------------------------------------------------------------------
     ! convergence criteria:
     ! (1) abs_error < density_tolerance, the normal criterium
     ! (2) if step size 'delta_eta' reaches limit of machine-precision, two or
     ! three times in a row
     !--------------------------------------------------------------------------
     abs_error = abs( error )
     if ( eta_iteration < 0.05_dp ) abs_error = abs_error * 1.E4_dp

     if ( abs( delta_eta ) < 100.0_dp * eta_iteration * machine_eps ) then
        count_num_noise = count_num_noise + 1
        if ( ( count_num_noise >= 2 .AND. abs_error < 100.0_dp * density_tolerance ) .OR.  &
               count_num_noise >= 3 )  abs_error = 0.0_dp         ! exit condition
     end if

  end do

  if ( abs_error > density_tolerance ) then
     write (*,*) 'density iteration failed', abs_error, density_tolerance
     if ( outp > 0 ) write (*,*) 'error step size criterion', abs( delta_eta ), 50.0_dp * eta_iteration * machine_eps
     if ( outp > 0 ) read (*,*)
     ! stop
  end if

  eta = eta_iteration
  if ( write_output .OR. outp > 0 ) write (*,*) 'now leaving density_iteration'
  if ( write_output .OR. outp > 0 ) write (*,*) ' '

  call pressure ( eta, pcalc, pcalc_z ) !JRE: pause here for sample calcs. call one last time just to check.

end subroutine density_iteration


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pressure_spinodal
!
!> \brief pressure_spinodal
!!
!! iterates the density until the derivative of pressure 'pcalc' to
!! density is equal to zero. A Newton-scheme is used for determining
!! the root to the objective function.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure_spinodal ( eta_start, eta, pcalc, pcalc_z )

  use PARAMETERS, only: dp
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z

  !-----------------------------------------------------------------------------
  real(dp), parameter                        :: tolerance = 1.E-5_dp
  integer, parameter                         :: max_iterations = 100

  integer                                    :: i
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: delta_eta
  real(dp)                                   :: error
  real(dp)                                   :: pcalc_z2, pcalc_z3
  !-----------------------------------------------------------------------------

  i = 0
  eta_iteration = eta_start

  !-----------------------------------------------------------------------------
  ! iterate density until d(p)/d(eta) = 0
  !-----------------------------------------------------------------------------

  error = tolerance + 1.0_dp
  do while ( abs( error ) > tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

     error = pcalc_z

     delta_eta = error / pcalc_z2
     if ( delta_eta >  0.05_dp ) delta_eta = 0.05_dp
     if ( delta_eta < -0.02_dp ) delta_eta = -0.02_dp

     eta_iteration   = eta_iteration - delta_eta
     ! write (*,'(a,i3,3G18.10)') 'iter',i, error, eta_iteration

     if ( eta_iteration > 0.7_dp )  eta_iteration = 0.5_dp
     if ( eta_iteration <= 0.0_dp ) eta_iteration = 1.E-16_dp

  end do

  eta = eta_iteration

end subroutine pressure_spinodal


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pressure_inflection_point
!
!> \brief calculates inflection point of pressure isotherm
!!
!! iterates the density until the second derivative of pressure 'pcalc' to
!! density is equal to zero (constant T,x). A Newton-scheme is used for determining
!! the root to the objective function.
!!
!! A recommended intitial value for eta_start is eta_start = 0.2
!! (at low temperatures, it might be advisable to use eta_start = 0.25 )
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pressure_inflection_point ( eta_start, eta, pcalc, pcalc_z, pcalc_z2, pcalc_z3, converg )

  use PARAMETERS, only: dp
  use properties, only: pressure
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: eta_start
  real(dp), intent(out)                      :: eta
  real(dp), intent(out)                      :: pcalc
  real(dp), intent(out)                      :: pcalc_z
  real(dp), intent(out)                      :: pcalc_z2
  real(dp), intent(out)                      :: pcalc_z3

  !-----------------------------------------------------------------------------
  real(dp), parameter                        :: tolerance = 1.E-2_dp
  integer, parameter                         :: max_iterations = 100

  integer                                    :: i
  real(dp)                                   :: eta_iteration
  real(dp)                                   :: delta_eta
  real(dp)                                   :: error
  logical                                    :: converg
  !-----------------------------------------------------------------------------

  eta_iteration = eta_start

  i = 0
  converg = .false.

  !-----------------------------------------------------------------------------
  ! iterate density until d^2(p)/d(eta)^2 = 0
  !-----------------------------------------------------------------------------

  error = tolerance + 1.0_dp
  do while ( abs(error) > tolerance .AND. i < max_iterations )

     i = i + 1

     call pressure ( eta_iteration, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

     error = pcalc_z2

     delta_eta = error / pcalc_z3
     if ( delta_eta >  0.08_dp ) delta_eta = 0.08_dp
     if ( delta_eta < -0.05_dp ) delta_eta = -0.05_dp

     eta_iteration   = eta_iteration - delta_eta
     ! write (*,'(a,i3,3G18.10)') 'iter inflex',i, error, eta_iteration

     if ( eta_iteration > 0.7_dp )  eta_iteration = 0.5_dp
     if ( eta_iteration <= 0.0_dp ) eta_iteration = 0.01_dp

  end do

  eta = eta_iteration
  
  if ( abs( error ) < tolerance ) converg = .true.

  if ( pcalc < 1.E-20_dp ) pcalc = 1.E-10_dp

end subroutine pressure_inflection_point
    
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module stability

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                              :: outp = 1
  real(dp), allocatable, dimension(:)             :: rhoi_coexisting

  private
  public :: stability_analysis, get_eta

CONTAINS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine stability_analysis
!> \brief stability_analysis
!!
!!
!! input :   t, p, x_feed(:)
!!           eta_trial    :  starting value for density of trial phase
!!           rhoi_coex(:) :  ( optional ) species densities of a phase that is
!!                           already known to be in phase equilibrium. This
!!                           subroutine ignores this phase as a trial phase
!!
!! output:   stable=.true.:  for stable phase or otherwise for a phase split
!!           rhoi_out(:)  :  species densities of trial phase (minimum
!!                           of the constrained tangent plane distance).
!!           rhoi_feed    :  density vector of feed phase
!!           ln_zi_phi    :  residual chemical potential of feed phase
!!                           mu_res(:)=ln( x_feed(:) * phi_feed(:) )
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_analysis ( t, p, x_feed, eta_start, stable, rhoi_feed, ln_zi_phi, rhoi_out, rhoi_coex )

  use PARAMETERS, only: KBOL30 !, machine_eps
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), intent(in)                            :: eta_start
  logical, intent(out)                            :: stable
  real(dp), dimension(ncomp), intent(out)         :: rhoi_feed
  real(dp), dimension(ncomp), intent(out)         :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_out
  real(dp), dimension(ncomp), intent(in), optional :: rhoi_coex
  !-----------------------------------------------------------------------------

  integer                                         :: trial_steps

  integer                                         :: i_trial

  real(dp)                                        :: tpd
  real(dp), dimension(ncomp)                      :: rhoi_trial
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu_feed
  real(dp)                                        :: delta_rhoi
  real(dp)                                        :: p_kT
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! is a coexisting phase to (t,p,x_feed) already known?  Then allocate array.
  !-----------------------------------------------------------------------------
  if ( present( rhoi_coex ) ) then
     allocate( rhoi_coexisting(ncomp) )
     rhoi_coexisting( 1:ncomp ) = rhoi_coex(1:ncomp)
  end if

  !-----------------------------------------------------------------------------
  ! starting values: one vapor & ncomp liquids, each with one species as dominating
  !-----------------------------------------------------------------------------
  trial_steps = ncomp + 1

  stable = .true.

  call chemical_potential_tp ( t, p, x_feed, eta_start, rhoi_feed, mu_res, mu=mu_feed )

  p_kT = p / ( KBOL30 * t )

  do i_trial = 1, trial_steps

     !--------------------------------------------------------------------------
     ! setting trial-phase mole-fractions (exp(ln_zi_phi): id.gas estimate for vapor)
     !--------------------------------------------------------------------------
     ln_zi_phi( 1:ncomp ) = log( x_feed( 1:ncomp ) ) + mu_res( 1:ncomp )

     call define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_trial )
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase eta ', get_eta( rhoi_trial(:) )
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase xi  ',  &
          rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )

     !--------------------------------------------------------------------------
     ! minimizing the objective fct. Phase split for values of tpd < 0.0
     !--------------------------------------------------------------------------
     call minimize_tpd ( t, mu_feed, p_kT, rhoi_feed, rhoi_trial, tpd )

     !--------------------------------------------------------------------------
     ! check for trivial solution ( where delta_rhoi = 0, approximately )
     !--------------------------------------------------------------------------
     delta_rhoi = maxval( abs( 1.0_dp - rhoi_trial(1:ncomp) / rhoi_feed(1:ncomp) ) )

     if ( tpd < -1.E-8_dp .AND. delta_rhoi > 0.001_dp ) then
        stable = .false.
        rhoi_out( 1:ncomp ) = rhoi_trial( 1:ncomp )
        if ( outp >= 1 ) then
           write (*,'(a,i4,G20.12)') ' unstable: N, tpd  ', i_trial, tpd
           write (*,'(a, 10G20.12)') ' unstable: eta, xi ', get_eta( rhoi_trial(:) ),  &
                 rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )
        end if
        exit
     end if

  end do

  !-----------------------------------------------------------------------------
  ! output to terminal
  !-----------------------------------------------------------------------------
  if ( outp >= 1 .AND. stable ) then
     write (*,'(a,i4,2G20.12)') ' stability analysis: stable phase'
  end if
  if ( outp >= 1 ) write (*,*) ' '

  if ( present( rhoi_coex ) ) deallocate( rhoi_coexisting )

end subroutine stability_analysis


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine define_initial_rhoi
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_guess )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i_trial
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), dimension(ncomp), intent(in)          :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_guess
  !-----------------------------------------------------------------------------

  integer                                         :: i
  integer                                         :: dominant_component
  real(dp)                                        :: x_dominant
  real(dp)                                        :: eta
  real(dp)                                        :: rho
  real(dp)                                        :: factor
  real(dp), dimension(ncomp)                      :: x
  !-----------------------------------------------------------------------------

  if ( i_trial == 1 ) then
     !------- vapor trial phase ------------------------------------------------
     x( 1:ncomp ) = exp( ln_zi_phi( 1:ncomp ) )
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )
     eta = 0.0001_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  else !if ( i_trial == 2 ) then
     !------- liquid trial phase ------------------------------------------------
     dominant_component = i_trial - 1
     x_dominant = 0.99_dp
     x( dominant_component ) = x_dominant
     factor = ( 1.0_dp - x_dominant ) / ( sum( x_feed(1:ncomp) ) - x_feed(dominant_component) )
     do i = 1, ncomp
        if ( i /= dominant_component) then
           x(i) = x_feed(i) * factor
        end if
     end do
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )  ! this line is actually not needed
     eta = 0.4_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  end if


end subroutine define_initial_rhoi


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE minimize_tpd
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine minimize_tpd ( t, mu_feed, p_kT, rhoi_feed, rhoi, tpd )

  use properties, only: chemical_potential_tp, chemical_potential_trho,  &
       Helmholtz_density
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), dimension(ncomp), intent(in)          :: mu_feed
  real(dp), intent(in)                            :: p_kT
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), dimension(ncomp), intent(inout)       :: rhoi
  real(dp), intent(out)                           :: tpd

  !-----------------------------------------------------------------------------
  integer                                         :: i_iter, i , iLoc(1)
  integer                                         :: critical_comp
  integer, parameter                              :: N_max_iter = 100
  real(dp), parameter                             :: tolerance = 1.E-3_dp
  real(dp)                                        :: scaled_tolerance
  real(dp), dimension(ncomp)                      :: rhoi_old
  real(dp)                                        :: tpd_old
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu
  real(dp)                                        :: error
  real(dp)                                        :: f_dens
  logical                                         :: trivial
  real(dp)                                        :: eta_new

  logical                                         :: repeat
  real(dp)                                        :: damp

  real(dp), dimension(ncomp,ncomp)                :: Hessian
  real(dp), dimension(ncomp,ncomp)                :: Hesse
  real(dp), dimension(ncomp,ncomp)                :: Hesse_ig
  real(dp)                                        :: eta_H
  logical                                         :: adjust_hessian
  integer, parameter                              :: it_max = 200
  integer                                         :: it_num
  integer                                         :: rot_num
  real(dp), dimension(ncomp,ncomp)                :: eigenvec
  real(dp), dimension(ncomp)                      :: eigenval
  real(dp), dimension(ncomp)                      :: in_out_vector
  real(dp)                                        :: DET

  logical                                         :: newton
  !-----------------------------------------------------------------------------

  rhoi_old(1:ncomp) = rhoi(1:ncomp)

  scaled_tolerance = tolerance
  damp = 1.0_dp     ! damping factor of direct substituation scheme
  
  i_iter = 0
  error = tolerance + 1.0_dp
  newton = .false.
  tpd_old = 1.E50_dp

  do while ( error > scaled_tolerance .AND. i_iter < N_max_iter )

     i_iter = i_iter + 1

     if ( .NOT.newton ) then

        !-----------------------------------------------------------------------
        ! case: direct substitution
        !-----------------------------------------------------------------------

        ! call direct_substitution_step ( t, mu_feed, p_kT, rhoi, mu_res, mu, tpd, damp, rhoi_old, tpd_old )
        ! if ( i_iter > 1 ) then
           call chemical_potential_trho ( t, rhoi, mu_res, mu )
        ! else
        !   x(:) = rhoi(:) / sum(rhoi(1:ncomp))
        !   p = p_kT * ( KBOL30 * t )
        !   call chemical_potential_tp ( t, p_kT*( KBOL30 * t ), x, 0.0001_dp, rhoi, mu_res, mu )
        ! end if
        g_function(1:ncomp) = mu(1:ncomp) - mu_feed(1:ncomp)

        !-----------------------------------------------------------------------
        ! constrain step size, if needed
        !-----------------------------------------------------------------------
        iLoc =  MINLOC( exp( mu_feed(1:ncomp) - mu(1:ncomp) ) )
        critical_comp = iLoc(1)
		damp = min( damp, exp( mu_feed(critical_comp) - mu(critical_comp) ) / 1.E-2_dp )
        if ( outp >= 2 ) write (*,*) 'minloc', exp( mu_feed(critical_comp) -mu(critical_comp)), damp        
        iLoc = MAXLOC( exp( mu_feed(1:ncomp) - mu(1:ncomp) ) )
        critical_comp = iLoc(1)
        damp = min( damp, 1.E2_dp / exp( mu_feed(critical_comp) - mu(critical_comp) ) )
        if ( outp >= 2 ) write (*,*) 'maxloc', exp( mu_feed(critical_comp) -mu(critical_comp)), damp

        repeat = .true.
        do while ( repeat )

           repeat = .false.

           !--------------------------------------------------------------------
           ! direct substitution step
           !--------------------------------------------------------------------
           rhoi(1:ncomp) = damp * exp( mu_feed(1:ncomp) - mu(1:ncomp) ) * rhoi_old(1:ncomp)  &
                + ( 1.0_dp - damp ) * rhoi_old(1:ncomp)

           !--------------------------------------------------------------------
           ! for too high densities, reduce rhoi
           !--------------------------------------------------------------------
           eta_new = get_eta ( rhoi )
           if ( eta_new > 0.6_dp ) rhoi(:) = rhoi(:) * 0.6_dp / eta_new
           ! do i = 1, ncomp
           !    write (*,*) 'rhoi',rhoi_old(i), rhoi(i)
           ! end do

           !--------------------------------------------------------------------
           ! calculate tangent plane distance tpd
           !--------------------------------------------------------------------
           call Helmholtz_density ( t, rhoi, f_dens )

           tpd = f_dens + p_kT - sum( mu_feed(1:ncomp) * rhoi(1:ncomp) )
           tpd = tpd * 1000.0_dp
           !write (70,'(i3,11G25.15)') i_iter, tpd, sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) ),  &
           !     eta_new, rhoi(1:ncomp)

           !--------------------------------------------------------------------
           ! if tpd is not successfully decreased, use stronger damping
           !--------------------------------------------------------------------
           if ( tpd > tpd_old .AND. damp > 0.04_dp ) then
              repeat = .true.
              damp = damp * 0.4_dp
              if ( outp >= 2 )  write (*,*) 'inner loop back', damp
           end if

        end do

        !-----------------------------------------------------------------------
        ! if tpd is successfully decreased
        !-----------------------------------------------------------------------
        if ( tpd < tpd_old + 1.E-5_dp ) then
           if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error DS    ', i_iter, tpd,  &
                sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) ), damp
           rhoi_old(1:ncomp) = rhoi(1:ncomp)
           tpd_old = tpd
           if (damp < 1.0_dp ) damp = min( 1.0_dp, damp * 10.0_dp )
           if ( i_iter >= 15 .OR. ( i_iter >= 10 .AND. eta_new > 0.25_dp ) ) then
              newton = .true.
              if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton'
           end if
        else
           if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton', i_iter, tpd,  &
                sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) )
           rhoi(1:ncomp) = rhoi_old(1:ncomp)
           tpd = tpd_old
           newton = .true.
        end if

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( sqrt(rhoi(1:ncomp))*g_function(1:ncomp) ) )

     else

        !-----------------------------------------------------------------------
        ! case Newton step
        !-----------------------------------------------------------------------

        !if ( i_iter == 1 ) then
        !   x(:) = rhoi(:) / sum(rhoi(1:ncomp))
        !   p = p_kT * ( KBOL30 * t )
        !   call chemical_potential_tp ( t, p_kT*( KBOL30 * t ), x, 0.0001_dp, rhoi, mu_res, mu )
        !end if

        call chemical_potential_trho ( t, rhoi, mu_res, mu, Hesse, Hesse_ig )
        g_function(1:ncomp) = ( mu(1:ncomp) - mu_feed(1:ncomp) ) * sqrt( rhoi(1:ncomp) )

        do i = 1, ncomp
           Hesse(i,1:ncomp) = Hesse(i,1:ncomp) * sqrt( rhoi(i)*rhoi(1:ncomp) )
           Hesse_ig(i,1:ncomp) = Hesse_ig(i,1:ncomp) * sqrt( rhoi(i)*rhoi(1:ncomp) ) ! unity matrix
        end do

        !-----------------------------------------------------------------------
        ! use method of Murray, by adding a unity matrix to Hessian, if:
        ! (1) H is not positive definite
        ! (2) step size is too large
        ! (3) objective function (tpd) does not descent
        !-----------------------------------------------------------------------
        adjust_hessian = .true.
        eta_H = 1.0_dp

        do while ( adjust_hessian )

           adjust_hessian = .false.

           Hessian(:,:) = Hesse(:,:) + eta_H * Hesse_ig(:,:)

           call jacobi_eigenvalue ( ncomp, Hessian, it_max, eigenvec, eigenval, it_num, rot_num )
           if ( outp >= 3 ) write (*,*) 'smallest eigenvalue',eigenval(1)

           if ( eigenval(1) < 0.001_dp .AND. eta_H < 20.0_dp ) then
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H I, eigenval=',eigenval(1)
              eta_H = eta_H + 0.5_dp
              adjust_hessian = .true.
              CYCLE      ! cycle, because of Hessian-criterion (1): H not positive definite
           end if

           !--------------------------------------------------------------------
           ! solving  AX = B  with Gauss-Jordan method.
           ! Matrix A='Hessian' and vector B='in_out_vector' are destroyed and
           ! redefined to the inverse of A and solution vector X, respectively.
           !--------------------------------------------------------------------

           in_out_vector = g_function
           call MATINVg ( ncomp, 1, Hessian, in_out_vector, DET )
           do i = 1,ncomp
              if ( abs( ( in_out_vector(i) / 2.0_dp )**2 / rhoi(i) ) > 0.99_dp ) then
                 adjust_hessian = .true.
              end if
           end do
           if ( adjust_hessian ) then
              eta_H = eta_H + 0.5_dp
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H II', eta_H
              CYCLE      ! cycle, because of Hessian-criterion (2): too large step-size
           end if

           !--------------------------------------------------------------------
           ! new density vector rhoi
           !--------------------------------------------------------------------
           rhoi(1:ncomp) = ( sqrt(rhoi_old(1:ncomp)) - in_out_vector(1:ncomp) / 2.0_dp )**2
           ! do i = 1, ncomp
           !    write (*,*) 'rhoi',rhoi_old(i), rhoi(i)
           ! end do

           !--------------------------------------------------------------------
           ! for too high densities, reduce rhoi
           !--------------------------------------------------------------------
           eta_new = get_eta( rhoi )
           if ( outp >= 3 )  write (*,*) 'eta_new',eta_new
           if ( eta_new > 0.6_dp ) rhoi(:) = rhoi(:) * 0.6_dp / eta_new

           !--------------------------------------------------------------------
           ! calculate tangent plane distance tpd
           !--------------------------------------------------------------------
           call Helmholtz_density ( t, rhoi, f_dens )

           tpd = f_dens + p_kT - sum( mu_feed(1:ncomp) * rhoi(1:ncomp) )
           tpd = tpd * 1000.0_dp

           if ( tpd > tpd_old .AND. eta_H < 30.0_dp ) then
              eta_H = eta_H + 0.5_dp
              adjust_hessian = .true.    ! cycle, because of Hessian-criterion (3): tpd does not descent
              if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H III'
           end if

        end do

        if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error Newton', i_iter, tpd, sum( abs( g_function(1:ncomp) ) )

        !-----------------------------------------------------------------------
        ! if tpd is successfully decreased
        !-----------------------------------------------------------------------
        if ( tpd < tpd_old + 1.E-5_dp ) then
           rhoi_old(1:ncomp) = rhoi(1:ncomp)
           tpd_old = tpd
        else
           if ( outp >= 1 ) write (*,*) 'minimize_tpd: no further decrease of objective fct possible', tpd
           rhoi(1:ncomp) = rhoi_old(1:ncomp)
           tpd = tpd_old
        end if

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( g_function(1:ncomp) ) )

     end if

     !--------------------------------------------------------------------------
     ! check if approaching trivial solution or identified coexisting phase, exit
     !--------------------------------------------------------------------------
     call check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )
     if ( trivial ) then
        tpd = 100.0_dp
        exit
     end if

     !--------------------------------------------------------------------------
     ! relax on convergence criterion for clearly unstable trial phase
     !--------------------------------------------------------------------------
     if ( tpd < - 0.1_dp ) scaled_tolerance = tolerance * 10.0_dp
     if ( tpd < - 1.0_dp ) scaled_tolerance = tolerance * 100.0_dp

  end do     

end subroutine minimize_tpd


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )

  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), intent(in)                            :: tpd
  logical, intent(out)                            :: trivial
  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: ratio
  real(dp)                                        :: sum_difference
  !-----------------------------------------------------------------------------

  trivial = .false.

  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards feed (rhoi_feed)
  !-----------------------------------------------------------------------------
  sum_difference = 0.0_dp
  do i = 1, ncomp
     if ( rhoi_feed(i) > 1.E-8_dp ) then
        ratio = rhoi(i) / rhoi_feed(i) - 1.0_dp
     else
        ratio = rhoi(i) - rhoi_feed(i)
     end if
     sum_difference = sum_difference + abs( ratio )
  end do
  sum_difference = sum_difference / real( ncomp, KIND=dp )

  if ( sum_difference < 0.01_dp .AND. tpd < 1.E-3_dp ) then
     trivial = .true.
     if ( outp >= 1 ) write (*,*) 'trivial solution, exit'
  end if

  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards coexising phase (rhoi_coexisting)
  !-----------------------------------------------------------------------------
  if ( allocated( rhoi_coexisting ) .AND. .NOT.trivial ) then

     sum_difference = 0.0_dp
     do i = 1, ncomp
        if ( rhoi_feed(i) > 1.E-8_dp ) then
           ratio = rhoi(i) / rhoi_coexisting(i) - 1.0_dp
        else
           ratio = rhoi(i) - rhoi_coexisting(i)
        end if
        sum_difference = sum_difference + abs( ratio )
     end do
     sum_difference = sum_difference / real( ncomp, KIND=dp )

     if ( sum_difference < 0.03_dp .AND. tpd < 1.E-2_dp ) then
        trivial = .true.
        if ( outp >= 1 ) write (*,*) 'found coexisting phase, exit'
     end if

  end if
     

end subroutine check_trivial_solution



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

function get_eta( rhoi )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi
  real(dp)                                        :: get_eta
  !-----------------------------------------------------------------------------

  get_eta = PI_6 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

end function get_eta

end module stability


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module stability_G

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  use stability, only: get_eta
  implicit none

  integer, parameter                              :: outp = 0
  real(dp), allocatable, dimension(:)             :: rhoi_coexisting

  private
  public :: stability_analysis_G

CONTAINS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine stability_analysis
!> \brief stability_analysis_G
!!
!!
!! input :   t, p, x_feed(:)
!!           eta_trial    :  starting value for density of trial phase
!!           rhoi_coex(:) :  ( optional ) species densities of a phase that is
!!                           already known to be in phase equilibrium. This
!!                           subroutine ignores this phase as a trial phase
!!
!! output:   stable=.true.:  for stable phase or otherwise for a phase split
!!           rhoi_out(:)  :  species densities of trial phase (minimum
!!                           of the constrained tangent plane distance).
!!           rhoi_feed    :  density vector of feed phase
!!           ln_zi_phi    :  residual chemical potential of feed phase
!!                           mu_res(:)=ln( x_feed(:) * phi_feed(:) )
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_analysis_G ( t, p, x_feed, eta_start, stable, rhoi_feed, ln_zi_phi, rhoi_out, rhoi_coex )

  use PARAMETERS, only: KBOL30 !, machine_eps
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), intent(in)                            :: eta_start
  logical, intent(out)                            :: stable
  real(dp), dimension(ncomp), intent(out)         :: rhoi_feed
  real(dp), dimension(ncomp), intent(out)         :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_out
  real(dp), dimension(ncomp), intent(in), optional :: rhoi_coex
  !-----------------------------------------------------------------------------

  integer                                         :: trial_steps

  integer                                         :: i_trial

  integer                                         :: i
  real(dp)                                        :: tpd
  real(dp), dimension(ncomp)                      :: rhoi_trial
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu_feed
  real(dp)                                        :: delta_rhoi
  real(dp)                                        :: p_kT
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! is a coexisting phase to (t,p,x_feed) already known?  Then allocate array.
  !-----------------------------------------------------------------------------
  if ( present( rhoi_coex ) ) then
     allocate( rhoi_coexisting( ncomp ) )
     rhoi_coexisting( 1:ncomp ) = rhoi_coex( 1:ncomp )
  end if

  !-----------------------------------------------------------------------------
  ! is the composition of one or more species = 0. Then, only here, set to 1.E-200
  !-----------------------------------------------------------------------------
  !do i = 1, ncomp
  !   if ( x_feed(i) <= 1.E-100_dp ) x_feed(i) = 1.E-100_dp
  !end do

  !-----------------------------------------------------------------------------
  ! starting values: one vapor & ncomp liquids, each with one species as dominating
  !-----------------------------------------------------------------------------
  trial_steps = ncomp + 1

  stable = .true.

  call chemical_potential_tp ( t, p, x_feed, eta_start, rhoi_feed, mu_res, mu=mu_feed )

  p_kT = p / ( KBOL30 * t )

  do i_trial = 1, trial_steps

     !--------------------------------------------------------------------------
     ! setting trial-phase mole-fractions (exp(ln_zi_phi): id.gas estimate for vapor)
     !--------------------------------------------------------------------------
     do i = 1, ncomp
        if ( x_feed(i) > 1.E-100_dp ) then
           ln_zi_phi(i) = log( x_feed(i) ) + mu_res(i)
        else
           ln_zi_phi(i) = log( 1.E-100_dp ) + mu_res(i)
        end if
     end do

     call define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_trial )
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase eta ', get_eta( rhoi_trial(:) )
     if ( outp >= 1 ) write (*,'(a,10G20.12)') ' trialphase xi  ',  &
          rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )

     !--------------------------------------------------------------------------
     ! minimizing the objective fct. Phase split for values of tpd < 0.0
     !--------------------------------------------------------------------------
     call minimize_tpd_G ( t, mu_feed, p_kT, rhoi_feed, rhoi_trial, tpd )

     !--------------------------------------------------------------------------
     ! check for trivial solution ( where delta_rhoi = 0, approximately )
     !--------------------------------------------------------------------------
     delta_rhoi = maxval( abs( 1.0_dp - rhoi_trial(1:ncomp) / rhoi_feed(1:ncomp) ) )

     if ( tpd < -1.E-8_dp .AND. delta_rhoi > 0.001_dp ) then
        stable = .false.
        rhoi_out( 1:ncomp ) = rhoi_trial( 1:ncomp )
        if ( outp >= 1 ) then
           write (*,'(a,i4,G20.12)') ' unstable: N, tpd  ', i_trial, tpd
           write (*,'(a, 10G20.12)') ' unstable: eta, xi ', get_eta( rhoi_trial(:) ),  &
                 rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )
        end if
        exit
     end if

  end do

  !-----------------------------------------------------------------------------
  ! output to terminal
  !-----------------------------------------------------------------------------
  if ( outp >= 1 .AND. stable ) then
     write (*,'(a,i4,2G20.12)') ' stability analysis: stable phase'
  end if
  if ( outp >= 1 ) write (*,*) ' '

  if ( present( rhoi_coex ) ) deallocate( rhoi_coexisting )

end subroutine stability_analysis_G


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine define_initial_rhoi
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine define_initial_rhoi ( i_trial, x_feed, ln_zi_phi, rhoi_guess )

  use PARAMETERS, only: PI_6
  use pcsaft_pure_and_binary_parameters, only: mseg
  use module_eos_derivatives, only: dhs
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i_trial
  real(dp), dimension(ncomp), intent(in)          :: x_feed
  real(dp), dimension(ncomp), intent(in)          :: ln_zi_phi
  real(dp), dimension(ncomp), intent(out)         :: rhoi_guess
  !-----------------------------------------------------------------------------

  integer                                         :: i
  integer                                         :: dominant_component
  real(dp)                                        :: x_dominant
  real(dp)                                        :: eta
  real(dp)                                        :: rho
  real(dp)                                        :: factor
  real(dp), dimension(ncomp)                      :: x
  !-----------------------------------------------------------------------------

  if ( i_trial == 1 ) then
     !------- vapor trial phase ------------------------------------------------
     x( 1:ncomp ) = exp( ln_zi_phi( 1:ncomp ) )
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )
     eta = 0.0001_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  else !if ( i_trial == 2 ) then
     !------- liquid trial phase ------------------------------------------------
     dominant_component = i_trial - 1
     x_dominant = 0.99_dp
     x( dominant_component ) = x_dominant
     factor = ( 1.0_dp - x_dominant ) / ( sum( x_feed(1:ncomp) ) - x_feed(dominant_component) )
     do i = 1, ncomp
        if ( i /= dominant_component) then
           x(i) = x_feed(i) * factor
        end if
     end do
     x( 1:ncomp ) = x( 1:ncomp ) / sum( x( 1:ncomp ) )  ! this line is actually not needed
     eta = 0.4_dp
     rho = eta / ( PI_6 * sum( x( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 ) )
     rhoi_guess( 1:ncomp ) = x( 1:ncomp ) * rho
  end if


end subroutine define_initial_rhoi


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE minimize_tpd_G
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine minimize_tpd_G ( t, mu_feed, p_kT, rhoi_feed, rhoi, tpd )

  use PARAMETERS, only: KBOL30
  use properties, only: chemical_potential_tp, chemical_potential_trho,  &
       Helmholtz_density
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), dimension(ncomp), intent(in)          :: mu_feed
  real(dp), intent(in)                            :: p_kT
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), dimension(ncomp), intent(inout)       :: rhoi
  real(dp), intent(out)                           :: tpd

  !-----------------------------------------------------------------------------
  integer                                         :: i
  integer                                         :: i_iter
  integer, parameter                              :: N_max_iter = 100
  real(dp), parameter                             :: tolerance = 1.E-6_dp
  real(dp)                                        :: scaled_tolerance
  real(dp), dimension(ncomp)                      :: x_old
  real(dp), dimension(ncomp)                      :: Y_old
  real(dp)                                        :: tpd_old
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp), dimension(ncomp)                      :: mu
  real(dp)                                        :: error
  logical                                         :: trivial
  
  real(dp)                                        :: p
  real(dp)                                        :: eta_in
  real(dp), dimension(ncomp)                      :: x
  real(dp), dimension(ncomp)                      :: di, lnphi_feed
  real(dp), dimension(ncomp)                      :: Y
  real(dp), dimension(ncomp)                      :: x_feed

  logical                                         :: newton
  !-----------------------------------------------------------------------------

  scaled_tolerance = tolerance
  
  !-----------------------------------------------------------------------------
  ! fugacity coeff. of feed-phase
  !-----------------------------------------------------------------------------
  do i = 1, ncomp

     if ( rhoi_feed(i) > 1.E-200_dp ) then

        lnphi_feed(i) = mu_feed(i) - log( rhoi_feed(i) ) - log( p_kT/sum( rhoi_feed(:) ) )
        x_feed(i) = rhoi_feed(i) / sum( rhoi_feed(:) )

     else

        lnphi_feed(i) = mu_feed(i) - log( 1.E-200_dp ) - log( p_kT/sum( rhoi_feed(:) ) )
        x_feed(i) = 1.E-200_dp

     end if

  end do

  di(:) = log( x_feed(:) ) + lnphi_feed(:)

  ! JG: this is not a good starting value: it is not used in the successive-substitution routine
  Y(1:ncomp) = x_feed(1:ncomp)
  Y_old(1:ncomp) = x_feed(1:ncomp)

  x(1:ncomp) = rhoi(1:ncomp) / sum( rhoi(1:ncomp) )
  x_old(1:ncomp) = x(1:ncomp)

  i_iter = 0
  error = tolerance + 1.0_dp
  newton = .false.
  tpd_old = 1.E50_dp

  
  do while ( error > scaled_tolerance .AND. i_iter < N_max_iter )

     i_iter = i_iter + 1

     if ( .NOT.newton ) then

        !-----------------------------------------------------------------------
        ! case: direct substitution
        !-----------------------------------------------------------------------

        p = p_kT * ( KBOL30 * t )

        eta_in = get_eta( rhoi )

        call chemical_potential_tp ( t, p, x, eta_in, rhoi, mu_res, mu = mu )

        !-----------------------------------------------------------------------
        ! direct substitution step
        !-----------------------------------------------------------------------
        Y(1:ncomp) = exp( di(1:ncomp) - mu_res(1:ncomp) )

        tpd = 1.0_dp + sum( Y(:) * ( log( Y(:) ) + mu_res(:) - di(:) - 1.0_dp ) )

        x(:) = Y(:) / sum( Y(:) )
        g_function(1:ncomp) = abs( x(1:ncomp) - x_old(1:ncomp) )

        if ( outp >= 1 ) write (*,'(a,i4,3G25.15)') ' N, tpd, error S.Sub.', i_iter, tpd,  &
                sum( abs( g_function(1:ncomp) ) )

        !-----------------------------------------------------------------------
        ! determine error
        !-----------------------------------------------------------------------

        error = sum( abs( g_function(1:ncomp) ) )

        x_old( 1:ncomp ) = x(1:ncomp)
        Y_old(:) = Y(:)                   ! this line is needed only for Newton-schme
        tpd_old = tpd

        !-----------------------------------------------------------------------
        ! exit condition towards Newton-schme
        !-----------------------------------------------------------------------
        if ( i_iter >= 15 .AND. error > scaled_tolerance .OR.  &
           ( tpd > tpd_old + 1.E-5_dp .AND. i_iter > 2 ) ) then
                 newton = .true.
                 if ( outp >= 2 ) write (*,'(a,i4,3G25.15)') 'exiting towards Newton'
        end if

     else

        call stability_G_Newton_step ( i_iter, t, p_kT, rhoi, di, Y_old, x_old, tpd_old, Y, x, tpd, error )

        !-----------------------------------------------------------------------
        ! store previous values
        !-----------------------------------------------------------------------
        Y_old(1:ncomp) = Y(1:ncomp)
        x_old(1:ncomp) = x(1:ncomp)
        tpd_old = tpd

     end if

     !--------------------------------------------------------------------------
     ! check if approaching trivial solution or identified coexisting phase, exit
     !--------------------------------------------------------------------------
     call check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )
     if ( trivial ) then
        tpd = 100.0_dp
        exit
     end if

     !--------------------------------------------------------------------------
     ! relax on convergence criterion for clearly unstable trial phase
     !--------------------------------------------------------------------------
     if ( tpd < - 0.01_dp ) scaled_tolerance = tolerance * 10.0_dp
     if ( tpd < - 0.1_dp ) scaled_tolerance = tolerance * 100.0_dp
     if ( tpd < - 0.1_dp .AND. i_iter > 5) scaled_tolerance = tolerance * 1000.0_dp

  end do     

end subroutine minimize_tpd_G


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE stability_G_Newton_step
!
! perform Newton step in stability analysis. A preceeding direct substitution
! step is required to have initial values for Y-vector.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_G_Newton_step ( i_iter, t, p_kT, rhoi, di, Y_old, x_old, tpd_old, Y, x, tpd, error )

  use PARAMETERS, only: KBOL30
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                         :: i_iter
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p_kT
  real(dp), dimension(ncomp), intent(in out)      :: rhoi
  real(dp), dimension(ncomp), intent(in)          :: di
  real(dp), dimension(ncomp), intent(in)          :: Y_old
  real(dp), dimension(ncomp), intent(in)          :: x_old
  real(dp), intent(in)                            :: tpd_old
  real(dp), dimension(ncomp), intent(out)         :: Y
  real(dp), dimension(ncomp), intent(in out)      :: x
  real(dp), intent(out)                           :: tpd
  real(dp), intent(out)                           :: error

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp), dimension(ncomp)                      :: g_function
  real(dp), dimension(ncomp)                      :: mu_res
  
  real(dp)                                        :: p
  real(dp)                                        :: eta_in

  real(dp), dimension(ncomp,ncomp)                :: Hessian
  real(dp), dimension(ncomp,ncomp)                :: Hesse
  real(dp), dimension(ncomp,ncomp)                :: Hesse_ig
  real(dp)                                        :: eta_H
  logical                                         :: adjust_hessian
  integer, parameter                              :: it_max = 200
  integer                                         :: it_num
  integer                                         :: rot_num
  real(dp), dimension(ncomp,ncomp)                :: eigenvec
  real(dp), dimension(ncomp)                      :: eigenval
  real(dp), dimension(ncomp)                      :: in_out_vector
  real(dp)                                        :: DET
  !-----------------------------------------------------------------------------

  p = p_kT * ( KBOL30 * t )

  eta_in = get_eta( rhoi )

  call chemical_potential_tp ( t, p, x, eta_in, rhoi, mu_res, lnphi_nj=Hesse )

  g_function(1:ncomp) = ( log( Y(1:ncomp) ) + mu_res(1:ncomp) - di(1:ncomp) ) * sqrt( Y(1:ncomp) )

  Hesse_ig(:,:) = 0.0_dp
  ! Hesse(:,:) = 0.0_dp
  do i = 1, ncomp
     Hesse(i,1:ncomp) = Hesse(i,1:ncomp) * sqrt( Y(i)*Y(1:ncomp) )
     Hesse(i,i) = Hesse(i,i) + log( Y(i) ) + mu_res(i) - di(i)
     Hesse_ig(i,i) = 1.0_dp
  end do

  !-----------------------------------------------------------------------
  ! use method of Murray, by adding a unity matrix to Hessian, if:
  ! (1) H is not positive definite
  ! (2) step size is too large
  ! (3) objective function (tpd) does not descent
  !-----------------------------------------------------------------------
  adjust_hessian = .true.
  eta_H = 1.0_dp

  do while ( adjust_hessian )

     adjust_hessian = .false.

     Hessian(:,:) = Hesse(:,:) + eta_H * Hesse_ig(:,:)

     call jacobi_eigenvalue ( ncomp, Hessian, it_max, eigenvec, eigenval, it_num, rot_num )
     if ( outp >= 3 ) write (*,*) 'smallest eigenvalue',eigenval(1)

     if ( eigenval(1) < 0.001_dp .AND. eta_H < 20.0_dp ) then
        if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H I, eigenval=',eigenval(1)
        eta_H = eta_H + 0.5_dp
        adjust_hessian = .true.
        CYCLE      ! cycle, because of Hessian-criterion (1): H not positive definite
     end if

     !--------------------------------------------------------------------
     ! solving  AX = B  with Gauss-Jordan method.
     ! Matrix A='Hessian' and vector B='in_out_vector' are destroyed and
     ! redefined to the inverse of A and solution vector X, respectively.
     !--------------------------------------------------------------------

     in_out_vector = g_function
     call MATINVg ( ncomp, 1, Hessian, in_out_vector, DET )
     do i = 1,ncomp
        if ( abs( ( in_out_vector(i) / 2.0_dp )**2 / Y_old(i) ) > 5.00_dp ) then
           adjust_hessian = .true.
        end if
     end do
     if ( adjust_hessian ) then
        eta_H = eta_H + 0.5_dp
        if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H II', eta_H
        CYCLE      ! cycle, because of Hessian-criterion (2): too large step-size
     end if

     !--------------------------------------------------------------------
     ! new vector Y
     !--------------------------------------------------------------------

     Y(1:ncomp) = ( sqrt(Y_old(1:ncomp)) - in_out_vector(1:ncomp) / 2.0_dp )**2

     tpd = 1.0_dp + sum( Y(:) * ( log( Y(:) ) + mu_res(:) - di(:) - 1.0_dp ) )

     x(:) = Y(:) / sum( Y(:) )

     if ( tpd > tpd_old .AND. eta_H < 30.0_dp .AND. i_iter >= 3 ) then
        eta_H = eta_H + 0.5_dp
        adjust_hessian = .true.    ! cycle, because of Hessian-criterion (3): tpd does not descent
        if ( outp >= 2 ) write (*,*) 'stability analysis: increasing eta_H III', tpd, tpd_old
     end if

     if ( outp >= 1 ) write (*,'(a,i4,4G25.15)') ' N, tpd, error Newton', i_iter, tpd,  &
          sum( abs( g_function(1:ncomp) ) ), sum( abs( x(1:ncomp) - x_old(1:ncomp) ) ), eta_H

  end do

  !-----------------------------------------------------------------------
  ! determine error
  !-----------------------------------------------------------------------
  error = sum( abs( g_function(1:ncomp) ) )
end subroutine stability_G_Newton_step


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine check_trivial_solution ( rhoi, rhoi_feed, tpd, trivial )
  implicit none
  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi
  real(dp), dimension(ncomp), intent(in)          :: rhoi_feed
  real(dp), intent(in)                            :: tpd
  logical, intent(out)                            :: trivial
  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: ratio
  real(dp)                                        :: sum_difference
  !-----------------------------------------------------------------------------
  trivial = .false.
  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards feed (rhoi_feed)
  !-----------------------------------------------------------------------------
  sum_difference = 0.0_dp
  do i = 1, ncomp
     if ( rhoi_feed(i) > 1.E-8_dp ) then
        ratio = rhoi(i) / rhoi_feed(i) - 1.0_dp
     else
        ratio = rhoi(i) - rhoi_feed(i)
     end if
     sum_difference = sum_difference + abs( ratio )
  end do
  sum_difference = sum_difference / real( ncomp, KIND=dp )
  if ( sum_difference < 0.01_dp .AND. tpd < 1.E-3_dp ) then
     trivial = .true.
     if ( outp >= 1 ) write (*,*) 'trivial solution, exit'
  end if
  !-----------------------------------------------------------------------------
  ! is trial phase (rhoi) converging towards coexising phase (rhoi_coexisting)
  !-----------------------------------------------------------------------------
  if ( allocated( rhoi_coexisting ) .AND. .NOT.trivial ) then
     sum_difference = 0.0_dp
     do i = 1, ncomp
        if ( rhoi_feed(i) > 1.E-8_dp ) then
           ratio = rhoi(i) / rhoi_coexisting(i) - 1.0_dp
        else
           ratio = rhoi(i) - rhoi_coexisting(i)
        end if
        sum_difference = sum_difference + abs( ratio )
     end do
     sum_difference = sum_difference / real( ncomp, KIND=dp )

     if ( sum_difference < 0.03_dp .AND. tpd < 1.E-2_dp ) then
        trivial = .true.
        if ( outp >= 1 ) write (*,*) 'found coexisting phase, exit'
     end if

  end if
end subroutine check_trivial_solution
end module stability_G
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!JRE CalculateSinglePhase was written by Vladimir. It includes stability analysis in the assessment of model properties. 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine CalculateSinglePhase (prp_id, res, ierr)
  use PARAMETERS, only: dp, KBOL, machine_eps
  use properties, only: state_trho, state_tp
  use stability, only: get_eta
  use stability_G, only: stability_analysis_G
  use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
implicit none
  real(dp) res  !result of calculation
  integer prp_id, ierr
  real(dp)                                        :: tF, pF
  real(dp), dimension(ncomp)                      :: xF
  real(dp), dimension(ncomp)                      :: rhoiL, rhoiV
  real(dp)                                        :: sL, sV
  real(dp)                                        :: hL, hV
  real(dp)                                        :: cpL, cpV
  real(dp)                                        :: eta_L, eta_V
  real(dp)                                        :: eta_start
  real(dp), dimension(ncomp)                      :: ln_zi_phi
  real(dp), dimension(ncomp)                      :: rhoi_trial
  real(dp)                                        :: molar_rho
  logical                                         :: stable_L, stable_V
    real(dp) s_resL, s_resG
  tF = t_input
  pF = p_input
  xF( 1 ) = x_input( 1 )
  !-----------------------------------------------------------------------------
  ! for liquid trial phase: test stability, then determine state variables
  !-----------------------------------------------------------------------------
  eta_start = 0.45_dp                        ! starting density for a liquid phase
  call stability_analysis_G ( tF, pF, xF, eta_start, stable_L, rhoiL, ln_zi_phi, rhoi_trial )
  eta_L = get_eta( rhoiL(:) )
  if ( stable_L ) call state_tp ( tF, pF, xF, eta_L, rhoiL, molar_rho, sL, hL, s_res=s_resL, cp=cpL )
  !-----------------------------------------------------------------------------
  ! for vapor trial phase: test stability, then determine state variables
  !-----------------------------------------------------------------------------
  eta_start = 1.E-5_dp                      ! starting density for a vapor phase
  call stability_analysis_G ( tF, pF, xF, eta_start, stable_V, rhoiV, ln_zi_phi, rhoi_trial )
  eta_V = get_eta( rhoiV(:) )
  if ( stable_V ) call state_tp ( tF, pF, xF, eta_V, rhoiV, molar_rho, sV, hV, s_res=s_resG, cp=cpV )
  !-----------------------------------------------------------------------------
  ! write output
  !-----------------------------------------------------------------------------
      res = 0
      if ( stable_L ) then
        if(prp_id==1) res=rhoiL(1)
        if(prp_id==-3) res=sL
        if(prp_id==-1) res=hL
        if(prp_id==12) res=cpL
      else if ( stable_V ) then
        if(prp_id==2) res=rhoiV(1)
        if(prp_id==-4) res=sV
        if(prp_id==-2) res=hV
        if(prp_id==14) res=cpV
      else
	    ierr = 1    !failed
      end if
      if (res.eq.0) ierr = 1
end subroutine CalculateSinglePhase
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!JRE These routines provide the interface from Gross's PcSaft to UaWrapper and CalcEos. They are at the end because the modules must come first.
!       
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine GetPcSaft(nComps,casrn,iErr)
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	!use input_output, only: read_problem_definition
	use GlobConst, only: LOUD, bVolCc_mol, isPcSaft, iEosOpt, SetNewEos !JRE LOUD=.FALSE. to silence I/O feedback.
	use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
	use PARAMETERS, only: dp, nsite
	use eos_constants
	use module_eos_derivatives, only: allocate_eos_quantities, initialize_eos_quantities,  &
       cross_association_parameters, aggregate_polar_parameters
	use pcsaft_pure_and_binary_parameters, only: load_pure_and_binary_parameters
    !use PARAMETERS, only: dp, KBOL, machine_eps
    implicit none
    integer nComps, casrn(nComps),iErr,iStat,i,L
	CHARACTER*128 masterDir
    character(16) aString
    logical :: parameter_assigned
    integer :: read_info
    character(LEN=16)    :: r_CAS_name   !actually, CASRN
    character(LEN=64)    :: alias, r_compna   !compound alias
    character (LEN=300)                        :: entire_line
    character (LEN=100)                        :: rest_of_line
    character(LEN=20)                          :: r_parameter_set
    real(dp)                                   :: r_mm
    real(dp)                                   :: r_mseg
    real(dp)                                   :: r_sigma
    real(dp)                                   :: r_epsilon_k
    real(dp)                                   :: r_dipole_moment
    real(dp)                                   :: r_quadru_moment
    !integer                                    :: r_nhb_typ
    !character(LEN=2)                           :: r_assoc_scheme
    !real(dp), dimension(nsite)                 :: r_nhb_no
    !real(dp)                                   :: r_kap_hb
    !real(dp), dimension(nsite,nsite)           :: r_eps_hb
	CHARACTER($MAXPATH)  CURDIR,parmFile,readFile !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	integer line
	!-----------------------------------------------------------------------------
	! read definition of problem
	!-----------------------------------------------------------------------------
	!call read_problem_definition - already done
	!  Get current directory
    if(LOUD)print*,'GetPcSaft: nComps,casrn=',nComps,casrn
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	masterDir=TRIM(curDir)
	parmFile=TRIM(masterDir)//'\Input\PcSaft_database\pcsaft_pure_parameters.txt' 
	if(LOUD)print*,'parmFile=',TRIM(parmFile)
	isPcSaft=.FALSE.
    iErr=0
	if(iEosOpt.ne.10 .and. iEosOpt.ne.15 .and. iEosOpt.ne.16)iErr=1
	if(LOUD.and.iErr==1)print*,'GetPcSaft: wrong iEosOpt=',iEosOpt
	if(iErr==1)return
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isTPT or other similar.
	isPcSaft=.TRUE.
	if(iEosOpt==10)readFile=TRIM(masterDir)//'\input\PcSaft_database\pcsaft_pure_parametersGross.txt' 
	if(iEosOpt==15)readFile=TRIM(masterDir)//'\input\PcSaft_database\pcsaft_pure_parametersGc.txt' 
	if(iEosOpt==16)readFile=TRIM(masterDir)//'\input\PcSaft_database\pcsaft_pure_parametersGcTb.txt'
	open(551,file=readFile)
	open(661,file=parmFile)
	read_info=0 
	do while(read_info==0)
		read (551, '(a)', IOSTAT = read_info) entire_line
		write(661,'(a)') TRIM(entire_line)
	enddo
	close(661)
	close(551)
	ncomp=nComps
	if ( .not. allocated(compna_input) ) allocate( compna_input( ncomp ) )  !names
	!if (.not. allocated(x_input)) allocate( x_input( ncomp ) )       ! mole fractions
	do i=1,nComps
        !create casrn string and find compound alias
        write(aString, *) casrn(i)
        aString=ADJUSTL((aString))
		!print*,'aString=',aString
        l=LEN_TRIM(aString)
        aString = aString(1:l-3) // '-' // aString(l-2:l-1) // '-' // aString(l:l)
        open (53, FILE = parmFile )
        parameter_assigned = .false.
        read_info = 0
		if(LOUD)print*,'iComp,CasNo=',i,'  ',TRIM(aString)
		line=1
        do while ( read_info == 0 .AND. .NOT.parameter_assigned )
            read ( 53, '(a)', IOSTAT = read_info) entire_line
			!if(LOUD)print*,TRIM(entire_line)
			if ( read_info /= 0 .and. LOUD) print*,'Error reading entire line=',line
            if ( read_info == 0 ) read (entire_line, *, ioStat=read_info ) r_CAS_name, r_compna ,  &
                    r_parameter_set , r_mm , r_mseg, r_sigma, r_epsilon_k ,  &
                    r_dipole_moment, r_quadru_moment, rest_of_line  !These are read but not stored here.
			!print*,r_mm,trim(r_parameter_set),' ',trim(entire_line)
			if ( read_info /= 0 .and. LOUD) print*,'Error reading internal line=',line
			line=line+1
			if(read_info > 0)read_info=0   ! ioStat < 0 if end of file or end of record is found.
            if ( trim(adjustl(aString)) == trim(adjustl(r_CAS_name)) ) then
                    alias = r_compna
                    parameter_assigned = .true.
					if(LOUD)print*,'Got it! iComp,compna=',i,'  ',TRIM(alias)
					if(LOUD)write(*,*)r_mm,r_mseg,r_sigma,r_epsilon_k, r_dipole_moment,r_quadru_moment, TRIM(rest_of_line)
					bVolCc_mol(i)=r_mseg*(r_sigma/10)**3*3.14159265D0/6*602.22 !bVol just needs to be close for eta initialization etc.
            end if
        end do
        close (53)
        if (.NOT.parameter_assigned) then
            iErr=iErr+1  !count up the number of compounds not found.
        else 
            compna_input(i) = alias   !Only compna_input is stored from this part. 
        endif
    enddo !i=1,nComps
    if( iErr.ne.0 )return
	!-----------------------------------------------------------------------------
	! allocate and initialize array quantities
	!-----------------------------------------------------------------------------
	call allocate_eos_quantities
	call initialize_eos_quantities
	if(LOUD)print*,'GetPcSaft: allocated and initialized.'
	!-----------------------------------------------------------------------------
	! read pure component and mixture parameters of PC-SAFT EOS
	!-----------------------------------------------------------------------------
	call load_pure_and_binary_parameters  !compna_input is used to load values
	if(LOUD)print*,'GetPcSaft: loaded.'
	call cross_association_parameters
	if(LOUD)print*,'GetPcSaft: cross_assoc ok.'
	call aggregate_polar_parameters
	if(LOUD)print*,'GetPcSaft: aggregated polar.'
	!-----------------------------------------------------------------------------
	! EOS constants
	!-----------------------------------------------------------------------------
	call load_eos_constants
	if(LOUD)print*,'GetPcSaft: Loaded EOS. Done.'
end subroutine GetPcSaft
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!JRE The Fugi routines provide the primary interface..
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 subroutine FugiPcSaftVtot(nComps,tKelvin,vTotCc,gmol,pMPa,zFactor,chemPoRes,iErr)
    use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
    use PARAMETERS, only: dp, nsite, NAV , RGAS
    use GlobConst, only: LOUD, uRes_RT, sRes_R, aRes_RT, hRes_RT, cpRes_R, cvRes_R, cmprsblty
	use properties !JRE add. Gives access to the chemical_potential___ and state___  subroutines.  
	!use stability, only: get_eta  !JRE add
	use eos_constants
	use module_eos_derivatives, only: kT,rho,allocate_eos_quantities, initialize_eos_quantities,  &
       cross_association_parameters, aggregate_polar_parameters
	use pcsaft_pure_and_binary_parameters, only: load_pure_and_binary_parameters
	!-----------------------------------------------------------------------------
	!JRE add
    implicit none
    integer nComps,iErr !,LIQ
    DoublePrecision tKelvin,pMPa,zFactor
    DoublePrecision gmol(nComps),chemPoRes(nComps)
    !real(dp) tlv, plv, rhoL, rhoV
    !real(dp) res
	integer i
	!real(dp)                                        :: sRes_R
	!real(dp)                                        :: hRes_RT
	!real(dp)                                        :: cpRes_R,cvRes_R,cmprsblty
	real(dp)                                        :: eta
	!real(dp)                                        :: eta_start
	real(dp), dimension(ncomp)                      :: fugacity ![=] MPa
	real(dp), dimension(ncomp)                      :: mu,mu_res ![=] J/mol
	real(dp), dimension(ncomp)                      :: rhoiJre ![=] molecules/Angst^3
	!real(dp), dimension(ncomp,ncomp)                :: w_rkrl, wig_rkrl  !JRE These are optional in chemical_potential_trho and not needed for Fugi at this time.
	real(dp)                                        :: molar_rho,totMol,vTotCc
	real(dp)                                        :: pcalc
	logical LOUDER
	LOUDER=.FALSE.
	!LOUDER=LOUD
    if(.not.allocated(x_input))allocate( x_input(ncomp) )
    ncomp=nComps
    t_input=tKelvin
    !p_input=pMPa*1.D6
    totMol=sum( gmol(1:ncomp) )
    molar_rho = totMol/vTotCc*1.D6 ! [=] gmol/m^3
    do i=1,ncomp
        x_input(i) = gmol(i)/totMol
        rhoiJre(i)=gmol(i)/totMol*molar_rho*(NAV/1.D30)
    enddo
    !call initialize_unit xxx must call GetPcSaft before Fugi.
    iErr=0
	if(LOUDER)print*,'FugiPcSaftVtot: t,p,x1=',t_input,( x_input(1:ncomp) )
	if(LOUDER)print*,'rho:',rhoiJre( 1:ncomp )   !number densities of each component
	if(LOUDER)print*,'calling chempo'
	call chemical_potential_trho ( t_input, rhoiJre, mu_res, mu) !, w_rkrl, wig_rkrl) !function overload: last two are optional.
	!eta = get_eta( rhoiJre )
	!rGasJre=RGAS FYI: RGAS=8.31445986144D0 by Gross. close enough.
	if(LOUDER)print*,'calling state_trho'
	!call state_trho ( t_input, rhoiJre, pcalc, molar_rho, sRes_R, hRes_RT, cmprsblty, cvRes_R )
	!subroutine state_trho ( tF, rhoi_in, pcalc, molar_rho, s, h, s_res, cp )
	call state_trhoJre ( t_input, rhoiJre, pcalc, molar_rho, eta, zFactor, sRes_R, hRes_RT, cmprsblty, cvRes_R, cpRes_R )
	!subroutine state_trhoJre ( tF, rhoi_in, pcalc, molar_rho, eta, ztot, sRes_R, hRes_RT, cmprsblty, cvRes_R, cpRes_R )
	!if(ABS(pcalc-p_input)/p_input > 1.D-3) iErr=3
	if(zFactor > 0)chemPoRes(1:nComp)=mu_res(1:ncomp)-LOG(zFactor)
	if(LOUDER)print*,'mu_res:',chemPoRes( 1:ncomp )
	pMPa=pcalc/1.d6
	fugacity(1:ncomp)=x_input(1:ncomp)*pMPa*exp( chemPoRes(1:ncomp) )
	if(LOUDER)print*,'fugacity(MPa):',fugacity( 1:ncomp )
	if(LOUDER)print*,'pcalc,eta,Z',pcalc,eta,zFactor
	if(LOUDER)print*,'rho_mol,S,H:',molar_rho, sRes_R,hRes_RT
	if(LOUDER)print*,'cmprsblty,CvRes/R,CpRes/R:',cmprsblty, cvRes_R, cpRes_R
	uRes_RT = hRes_RT-(zFactor-1)
	aRes_RT = uRes_RT - sRes_R
	if(LOUDER)print*,'FugiPcSaftVtot Done!'
	!pause 'check output'
	!JRE end
	return
end subroutine FugiPcSaftVtot

subroutine FugiPcSaft(nComps,tKelvin,pMPa,xFrac,LIQ,zFactor,chemPoRes,iErr)
    use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
    use PARAMETERS, only: dp, nsite, NAV, RGAS
    use GlobConst, only: LOUD, uRes_RT, sRes_R, aRes_RT, hRes_RT, cpRes_R, cvRes_R, cmprsblty,etaPass
	use properties !JRE add
	use stability, only: get_eta  !JRE add
	use eos_constants
	use module_eos_derivatives, only: kT,rho,allocate_eos_quantities, initialize_eos_quantities,  &
       cross_association_parameters, aggregate_polar_parameters
	use pcsaft_pure_and_binary_parameters, only: load_pure_and_binary_parameters
	!-----------------------------------------------------------------------------
	!JRE add
	implicit none
	integer nComps,iErr,LIQ
	DoublePrecision tKelvin,pMPa,zFactor !,rGasJre
	DoublePrecision xFrac(nComps),chemPoRes(nComps)
	!integer i
	!real(dp)                                        :: sRes_R	 !these are declared in GlobConst
	!real(dp)                                        :: hRes_RT
	!real(dp)                                        :: cpRes_R,cvRes_R,cmprsblty
	real(dp)                                        :: eta
	real(dp)                                        :: eta_start
	real(dp), dimension(ncomp)			:: mu_res, rhoiJre, fugacity
	!real(dp), dimension(ncomp)                      :: mu				  !JRE These are optional in chemical_potential_tp and not needed for Fugi at this time.
	!real(dp), dimension(ncomp,ncomp)                :: lnphi_nj		 !JRE they can be activated by tacking them on to the end of the call to chemical_potential_tp.
	real(dp)                                        :: molar_rho
	real(dp)                                        :: pcalc
	logical LOUDER  !JRE This can be set to LOUD in order to facilitate debugging if desired.
    LOUDER=.FALSE.
    LOUDER=LOUD
    ncomp=nComps
    t_input=tKelvin
    p_input=pMPa*1.D6
	!if(LIQ==3)LIQ=1	 !No distinction for PcSaft for now.
    if(.not.allocated(x_input))allocate( x_input(ncomp) )
    x_input(1:ncomp) = xFrac(1:ncomp)/sum( xFrac(1:ncomp) )
    !call initialize_unit xxx must call GetPcSaft before Fugi.
    iErr=0
	!print*,'Calling consoleIO'
	!call consoleIO !JRE get T,P, ncomp,compnames.
	!JRE end
	if(LOUDER)print*,'FugiPcSaft: t,p=',t_input, p_input
	if(LOUDER)print*,'FugiPcSaft: x=',x_input(1:ncomp)
	eta_start=1.d-5
	if( LIQ==1 .or. LIQ==3)eta_start=0.45D0
	!print*,'calling chempo'
	call chemical_potential_tp ( t_input, p_input, x_input, eta_start, rhoiJre, mu_res) !, mu, lnphi_nj ) !function overload, last two are optional.
	chemPoRes(1:ncomp)=mu_res(1:ncomp)
	if(LOUDER)print*,'mu_res:',chemPoRes( 1:ncomp )
	fugacity(1:ncomp)=x_input(1:ncomp)*pMPa*exp( chemPoRes(1:ncomp) )
	if(LOUDER)print*,'fugacity(MPa):',fugacity( 1:ncomp )
	!print*,'calling state_trho'
	call state_trhoJre ( t_input, rhoiJre, pcalc, molar_rho, eta, zFactor, sRes_R, hRes_RT, cmprsblty, cvRes_R, cpRes_R )
	if(ABS(pcalc-p_input)/p_input > 1.D-3) iErr=13
	if(LOUDER)print*,'pcalc,eta,Z',pcalc,eta,zFactor
	if(LOUDER)print*,'rho_mol,S,H:',molar_rho, sRes_R,hRes_RT
	if(LOUDER)print*,'cmprsblty,CvRes/R,CpRes/R:',cmprsblty, cvRes_R, cpRes_R
	etaPass=eta !Needed for PsatEar.
    uRes_RT = hRes_RT-(zFactor-1)
    aRes_RT = uRes_RT - sRes_R
    if(LOUDER)print*,'FugiPcSaft Done!'
	!pause 'check output'
	!JRE end
	return
end subroutine FugiPcSaft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine QueryParPurePcSaft(iComp,iParm,value,iErr)
	USE pcsaft_pure_and_binary_parameters      !Just for PcSaft
	IMPLICIT NONE
	DoublePrecision value
	Integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	! The following were omitted JRE 20210421
	!character(LEN=20), public, allocatable, dimension(:) :: parameter_set  ! flag identifying a pure comp. param. set
	!character(LEN=2), public, allocatable, dimension(:)  :: assoc_scheme   ! association in notation of Huang & Radosz
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

	iErr=0
	if(iParm==1)then
		value=mseg(iComp)
	elseif(iParm==2)then
		value=sigma(iComp)
	elseif(iParm==3)then
		value=epsilon_k(iComp)
	elseif(iParm==4)then
		value=kap_hb(iComp,iComp)
	elseif(iParm==5)then
		value=eps_hb(iComp,iComp,1,2) !This assumes only one bonding site type per molecule. OK for 20210421.
	elseif(iParm==6)then
		value=dipole_moment(iComp)
	elseif(iParm==7)then
		value=quadru_moment(iComp)
	elseif(iParm==8)then
		value=nhb_typ(iComp)
	elseif(iParm==9)then
		value=nhb_no(iComp,1) !This assumes only one bonding site type per molecule. OK for 20210421.
	else
		iErr=1
	endif
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine SetParPurePcSaft(iComp,iParm,value,iErr)
	USE pcsaft_pure_and_binary_parameters      !Just for PcSaft
	IMPLICIT NONE
	DoublePrecision value
	Integer iComp,iParm,iErr
  !-----------------------------------------------------------------------------
  ! pure component parameters
  !-----------------------------------------------------------------------------
  ! The following were omitted JRE 20210421
  !character(LEN=20), public, allocatable, dimension(:) :: parameter_set  ! flag identifying a pure comp. param. set
  !character(LEN=2), public, allocatable, dimension(:)  :: assoc_scheme   ! association in notation of Huang & Radosz
	iErr=0
	if(iParm==1)then
		mseg(iComp)=value
	elseif(iParm==2)then
		sigma(iComp)=value
	elseif(iParm==3)then
		epsilon_k(iComp)=value
	elseif(iParm==4)then
		kap_hb(iComp,iComp)=value
	elseif(iParm==5)then
		eps_hb(iComp,iComp,1,1)=value	!This assumes only one bonding site type per molecule. OK for 20210421.
	elseif(iParm==6)then
		dipole_moment(iComp)=value
	elseif(iParm==7)then
		quadru_moment(iComp)=value
	elseif(iParm==8)then
		nhb_typ(iComp)=NINT(value)  ! roundup to ensure integer is right
	elseif(iParm==9)then
		nhb_no(iComp,1)=value
	else
		iErr=1
	endif
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine QueryParMixPcSaft(iParm,value,iErr)
	USE pcsaft_pure_and_binary_parameters      !Just for PcSaft.
	DoublePrecision value
	!! kij(i,j)          binary correction parameter acting on dispersion term, [/]\n
	!! lij(i,j)          asymmetric binary correction para. (Tang and Gross, 2010), [/]\n
	!! kij_assoc(i,j)    correction para. for association energy. [/]\n
	iErr=0
	if(iParm==1)then
		value=Kij(1,2)
	elseif(iParm==2)then
		value=Lij(1,2)
	elseif(iParm==3)then
		value=Kij_Assoc(1,2)
	else
		iErr=1
	endif
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine SetParMixPcSaft(iParm,value,iErr) ! called from UaEosTools.
	USE pcsaft_pure_and_binary_parameters      !Just for PcSaft.
	DoublePrecision value
	!! kij(i,j)          binary correction parameter acting on dispersion term, [/]\n
	!! lij(i,j)          asymmetric binary correction para. (Tang and Gross, 2010), [/]\n
	!! kij_assoc(i,j)    correction para. for association energy. [/]\n
	!! These values are initially read in: subroutine pcsaft_binary_parameters
	iErr=0
	if(iParm==1)then
		Kij(1,2)=value
		Kij(2,1)=value
	elseif(iParm==2)then
		Lij(1,2)=value
		Lij(2,1)=value
	elseif(iParm==3)then
		Kij_Assoc(1,2)=value
		Kij_Assoc(2,1)=value
	else
		iErr=1
	endif
	return
	end

