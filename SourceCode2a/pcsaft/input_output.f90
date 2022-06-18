!> \file input_output.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module input_output
!> \brief module for writing properties to output-files
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module input_output

  use PARAMETERS, only: dp, NAV
  implicit none

  character (LEN=4), public                  :: unit_p
  character (LEN=1), public                  :: unit_t

  real(dp), public                           :: u_in_T
  real(dp), public                           :: u_in_p

  private
  public :: read_problem_definition, pure_output_state, pure_output_coex,  &
       mixture_output_state, mixture_output_bin_coex, get_mass_density, output_flash

contains

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE read_problem_definition
!> \brief  This subroutine draws user input data from the input file
!!
!! This subroutine draws user input data from the file
!!     INPUT.INP      which is located:    ./INPUT_FILE/INPUT.INP
!!
!! SUMMERY OF OUTPUT VARIABLES:
!! t_input        temperature
!! p_input        pressure
!! ncomp          number of components
!! x              array of feed molefractions. For binary mixtures
!!                most calculation options don't require a definition
!!                of feed molefractions. Rather, values are generated
!!                automatically if all feed molefractions are set to
!!                zero in the INPUT-file.
!! compna_input   name (string) and molec.mass of pure component
!! unit_t,unit_p  units of T and p for input and output
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine read_problem_definition

  use PARAMETERS, only: dp, KBOL, machine_eps
  use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input, compna_input
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  integer                                    :: i_ncomp, ioErr
  integer                                    :: read_info
  character (LEN=50)                         :: filename
  character (LEN=30), dimension(100)         :: reading1
  real(dp), dimension(100)                   :: reading2
  real(dp)                                   :: sumfeed
  logical                                    :: filefound
  !-----------------------------------------------------------------------------

  filename = '..\input_file\INPUT.txt'
  filename = 'INPUT.txt'
	inquire (FILE=filename, EXIST = filefound)
!  if ( filefound ) then
     open ( 32, FILE = filename, ioStat=ioErr )
!  else
	if( ioErr .ne. 0)print*,'after opening: ioStat=',ioErr
    ! write (*,*) ' input_output: FOLLOWING FILE CANNOT BE OPENED: ', filename
!     stop
!  endif

  unit_t = 'K'
  unit_p = 'MPA'
  read (32,*,ioStat=ioErr) t_input,  p_input
  !read (32,*,ioStat=ioErr) t_input, unit_t, p_input, unit_p
  !read ( 32,'(a50)',ioStat=ioErr ) filename
  if(ioErr.ne.0)print*,'failed to read t,p'
  !print*,'dumstring=',filename
  !read ( 32,*,ioStat=ioErr) t_input, unit_t,  p_input, unit_p
  if(ioErr.ne.0)print*,' input_output: trouble reading T,P line'
  print*,'T,P=',t_input,p_input

  i_ncomp = 0
  read_info = 0
  do while ( read_info == 0 )
     i_ncomp = i_ncomp + 1
     read ( 32, *, IOSTAT=read_info ) reading1( i_ncomp ), reading2( i_ncomp )
	  !if(read_info .ne.0)print*,' input_output: trouble reading i_ncomp line. i_ncomp=',i_ncomp
  end do

  CLOSE (32)

  ncomp = i_ncomp - 1                        ! number of species in problem

  allocate( compna_input( ncomp ) )

  compna_input( 1:ncomp ) = reading1( 1:ncomp )    ! names of all species of the problem
  print*,'ncomp,comps=',ncomp,(compna_input(i_ncomp),i_ncomp=1,ncomp)

  allocate( x_input( ncomp ) )                     ! mole fractions of all species of the problem
  x_input( 1:ncomp ) = reading2( 1:ncomp )
  sumfeed = sum( x_input( 1:ncomp ) )

  if ( abs( sumfeed ) > machine_eps .and. abs( sumfeed - 1.0_dp ) > machine_eps ) then
     x_input(1:ncomp) = x_input(1:ncomp) / sumfeed
  end if

  !-----------------------------------------------------------------------------
  ! convert units to [p]=bar and [T]=K, for various possible units in input file
  !-----------------------------------------------------------------------------
  if ( unit_t == 'C' ) then
     u_in_t = 273.15_dp
  else
     u_in_t = 0.0_dp
  end if

  if ( unit_p == 'bar' ) then
     u_in_p = 1.E5_dp
  else if ( unit_p == 'mbar' ) then
     u_in_p = 1.E2_dp
  else if ( unit_p == 'MPa' ) then
     u_in_p = 1.E6_dp
  else if ( unit_p == 'kPa' ) then
     u_in_p = 1.E3_dp
  else if ( unit_p == '*' ) then
     u_in_p = 1.E30_dp * KBOL
  else
     u_in_p = 1.0_dp
  end if

  t_input = t_input + u_in_t
  p_input = p_input * u_in_p

end subroutine read_problem_definition


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_output_state
!> \brief writing properties (t, p, mass_density, s, h, ... ) to output-file
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_output_state ( t, p, rho, s, h, cp, header )

  use BASIC_VARIABLES, only: ncomp
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t, p, rho, s, h
  real(dp), intent(in), optional             :: cp
  logical, optional                          :: header
  !-----------------------------------------------------------------------------
  real(dp)                                   :: mass_density
  real(dp), dimension(ncomp)                 :: rhoi
  real(dp), dimension(ncomp)                 :: mass_frac
  logical                                    :: header_flag
  logical                                    :: is_open
  !-----------------------------------------------------------------------------

  if ( present( header ) ) then
     header_flag = header
  else
     header_flag = .false.
  end if


  if ( header ) then

     !--------------------------------------------------------------------------
     ! inital preparation: open output file and write header
     !--------------------------------------------------------------------------

     inquire( unit=40, opened=is_open )
     if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

     if ( present( cp ) ) then
        write (40,'(6a)') 'T/',unit_t,'  p/',unit_p,'  density/(kg/m**3)  density/(mol/m**3)  ',  &
                                'h/(J/mol)  s/(J/mol/K)  cpres/(J/mol/K)'
     else
        write (40,'(6a)') 'T/',unit_t,'  p/',unit_p,'  density/(kg/m**3)  density/(mol/m**3)  ',  &
                                'h/(J/mol)  s/(J/mol/K)'
     end if

  else

     !--------------------------------------------------------------------------
     ! output to terminal
     !--------------------------------------------------------------------------
     write (*,'(t2,a,f7.2,3a,f9.4,2a)') ' T =',t-u_in_T,' ',unit_t,'   p =',p/u_in_p,' ',unit_p

     !--------------------------------------------------------------------------
     ! output to file
     !--------------------------------------------------------------------------

     rhoi(:) = rho
     call get_mass_density ( rhoi, mass_density, mass_frac )

     if ( present( cp ) ) then
        write (40,'(7(2x,f18.8))') t-u_in_t, p/u_in_p, mass_density, rho*1.E30/NAV, h, s, cp
     else
        write (40,'(6(2x,f18.8))') t-u_in_t, p/u_in_p, mass_density, rho*1.E30/NAV, h, s
     end if

  end if

end subroutine pure_output_state


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pure_output_coex
!> \brief writing liquid-vapor properties of pure substance to output-file
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_output_coex ( t, p, rhoL, rhoV, sL, sV, hL, hV, cpL, cpV, header )

  use BASIC_VARIABLES, only: ncomp
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t, p, rhoL, rhoV, sL, sV, hL, hV
  real(dp), intent(in), optional             :: cpL
  real(dp), intent(in), optional             :: cpV
  logical, optional                          :: header
  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                 :: rhoi
  real(dp)                                   :: mass_density_L, mass_density_V
  real(dp), dimension(ncomp)                 :: mass_frac
  logical                                    :: header_flag
  logical                                    :: is_open
  !-----------------------------------------------------------------------------

  if ( present( header ) ) then
     header_flag = header
  else
     header_flag = .false.
  end if


  if ( header ) then

     !--------------------------------------------------------------------------
     ! inital preparation: open output file and write header
     !--------------------------------------------------------------------------

     inquire( unit=40, opened=is_open )
     if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

     if ( present( cpL ) .OR. present( cpV ) ) then
        write (40,'(7a)') 'T/',unit_t,'  p/',unit_p,'  density_L/(kg/m**3)  density_V/(kg/m**3)  ',  &
             'density_L/(mol/m**3)  density_V/(mol/m**3)  sL/(J/mol)  sV/(J/mol)  ', &
             'hL/(J/mol/K)  hV/(J/mol/K)  cpL/(J/mol/K)  cpV/(J/mol/K)'
     else
        write (40,'(7a)') 'T/',unit_t,'  p/',unit_p,'  density_L/(kg/m**3)  density_V/(kg/m**3)  ',  &
             'density_L/(mol/m**3)  density_V/(mol/m**3)  sL/(J/mol)  sV/(J/mol)  ', &
             'hL/(J/mol/K)  hV/(J/mol/K)'
     end if

  else

     !--------------------------------------------------------------------------
     ! output to terminal
     !--------------------------------------------------------------------------
     write (*,'(t2,a,f7.2,3a,f9.4,2a)') ' T =',t-u_in_T,' ',unit_t,'   p =',p/u_in_p,' ',unit_p

     !--------------------------------------------------------------------------
     ! output to file
     !--------------------------------------------------------------------------

     rhoi(:) = rhoL
     call get_mass_density ( rhoi, mass_density_L, mass_frac )
     rhoi(:) = rhoV
     call get_mass_density ( rhoi, mass_density_V, mass_frac )

     if ( present( cpL ) .OR. present( cpV ) ) then
        write (40,'(12(2x,f18.8))') t-u_in_t, p/u_in_p, mass_density_L, mass_density_V,  &
             rhoL*1.E30/NAV, rhoV*1.E30/NAV, sL, sv, hL, hV, cpL, cpV
     else
        write (40,'(10(2x,f18.8))') t-u_in_t, p/u_in_p, mass_density_L, mass_density_V,  &
             rhoL*1.E30/NAV, rhoV*1.E30/NAV, sL, sv, hL, hV
     end if

  end if

end subroutine pure_output_coex


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine mixture_output_state
!> \brief writing properties (t, p, mass_density, s, h, ... ) to output-file
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine mixture_output_state ( t, p, x, rhoi, s, h, cp, header )

  use BASIC_VARIABLES, only: ncomp
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x
  real(dp), dimension(ncomp), intent(in)     :: rhoi
  real(dp), intent(in)                       :: s
  real(dp), intent(in)                       :: h
  real(dp), intent(in), optional             :: cp
  logical, optional                          :: header
  !-----------------------------------------------------------------------------
  real(dp)                                   :: mass_density
  real(dp)                                   :: rho
  real(dp), dimension(ncomp)                 :: mass_frac
  logical                                    :: header_flag
  logical                                    :: is_open
  character(LEN=3), allocatable, dimension(:):: x_lable
  character(LEN=3)                          :: char_len
  !-----------------------------------------------------------------------------

  if ( present( header ) ) then
     header_flag = header
  else
     header_flag = .false.
  end if


  if ( header ) then

     !--------------------------------------------------------------------------
     ! inital preparation: open output file and write header
     !--------------------------------------------------------------------------

     inquire( unit=40, opened=is_open )
     if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

     write (char_len,'(I3)') ncomp + 6
     allocate( x_lable( ncomp ) )

     x_lable( 1:ncomp ) = '  x'
     if ( present( cp ) ) then
        write (40,'('//char_len//'a)') 'T/',unit_t,'  p/',unit_p, x_lable(1:ncomp),'  density/(kg/m**3)  ',  &
                                'density/(mol/m**3)  h/(J/mol)  s/(J/mol/K)  cpres/(J/mol/K)'
     else
        write (40,'('//char_len//'a)') 'T/',unit_t,'  p/',unit_p, x_lable(1:ncomp),'  density/(kg/m**3)  ',  &
                                'density/(mol/m**3)  h/(J/mol)  s/(J/mol/K)'
     end if

     deallocate( x_lable )

  else

     !--------------------------------------------------------------------------
     ! output to terminal
     !--------------------------------------------------------------------------
     write (*,'(t2,a,f7.2,3a,f9.4,2a)') ' T =',t-u_in_T,' ',unit_t,'   p =',p/u_in_p,' ',unit_p

     !--------------------------------------------------------------------------
     ! output to file
     !--------------------------------------------------------------------------

     call get_mass_density ( rhoi, mass_density, mass_frac )

     rho = sum( rhoi( 1:ncomp ) )

     if ( present( cp ) ) then
        write (char_len,'(I3)') ncomp + 7
        write (40,'('//char_len//'(2x,f18.8))') t-u_in_t, p/u_in_p, x(1:ncomp), mass_density, rho*1.E30/NAV, h, s, cp
     else
        write (char_len,'(I3)') ncomp + 6
        write (40,'('//char_len//'(2x,f18.8))') t-u_in_t, p/u_in_p, x(1:ncomp), mass_density, rho*1.E30/NAV, h, s
     end if

  end if

end subroutine mixture_output_state


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine mixture_output_bin_coex
!> \brief writing properties (t, p, mass_density, s, h, ... ) to output-file
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine mixture_output_bin_coex ( t, p, x1, x2, rhoi1, rhoi2, s1, s2, h1, h2, cp1, cp2, header )

  use BASIC_VARIABLES, only: ncomp
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x1
  real(dp), dimension(ncomp), intent(in)     :: x2
  real(dp), dimension(ncomp), intent(in)     :: rhoi1
  real(dp), dimension(ncomp), intent(in)     :: rhoi2
  real(dp), intent(in)                       :: s1
  real(dp), intent(in)                       :: s2
  real(dp), intent(in)                       :: h1
  real(dp), intent(in)                       :: h2
  real(dp), intent(in), optional             :: cp1
  real(dp), intent(in), optional             :: cp2
  logical, optional                          :: header
  !-----------------------------------------------------------------------------
  real(dp)                                   :: mass_density1, mass_density2
  real(dp), dimension(ncomp)                 :: mass_frac1, mass_frac2
  logical                                    :: header_flag
  logical                                    :: is_open
  !-----------------------------------------------------------------------------

  if ( present( header ) ) then
     header_flag = header
  else
     header_flag = .false.
  end if


  if ( header ) then

     !--------------------------------------------------------------------------
     ! inital preparation: open output file and write header
     !--------------------------------------------------------------------------

     inquire( unit=40, opened=is_open )
     if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

     if ( present( cp1 ) .OR. present( cp2 ) ) then
        write (40,'(6a)') ' T/',unit_t,'  p/',unit_p, ' x1L  x1V  densityL/(kg/m**3)  ',  &
                'densityV/(kg/m**3)  hL/(J/mol)  hV/(J/mol)  sL/(J/mol/K)  sV/(J/mol/K)  cpresL/(J/mol/K)  cpresV/(J/mol/K)'
     else
        write (40,'(6a)') ' T/',unit_t,'  p/',unit_p, ' x1L  x1V  densityL/(kg/m**3)  ',  &
                'densityV/(kg/m**3)  hL/(J/mol)  hV/(J/mol)  sL/(J/mol/K)  sV/(J/mol/K)'
     end if

  else

     !--------------------------------------------------------------------------
     ! output to terminal
     !--------------------------------------------------------------------------
     write (*,'(t2,a,f7.2,3a,f9.4,2a,a,2f10.6)') ' T =',t-u_in_T,' ',unit_t,'   p =',p/u_in_p,' ',unit_p,' x1_L x1_V',x1(1),x2(1)

     !--------------------------------------------------------------------------
     ! output to file
     !--------------------------------------------------------------------------

     call get_mass_density ( rhoi1, mass_density1, mass_frac1 )
     call get_mass_density ( rhoi2, mass_density2, mass_frac2 )

     !rho1 = 1.E30 / NAV * sum( rhoi1( 1:ncomp ) )       ! molar density
     !rho2 = 1.E30 / NAV * sum( rhoi2( 1:ncomp ) )       ! molar density

     if ( present( cp1 ) .OR. present( cp2 ) ) then
        write (40,'(14(2x,f18.8))') t-u_in_t, p/u_in_p, x1(1), x2(1), mass_density1, mass_density2, h1, h2, s1, s2, cp1, cp2
     else
        write (40,'(12(2x,f18.8))') t-u_in_t, p/u_in_p, x1(1), x2(1), mass_density1, mass_density2, h1, h2, s1, s2
     end if

  end if

end subroutine mixture_output_bin_coex


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine output_flash
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine output_flash ( t, p, x1, x2, x3, rhoi1, rhoi2, rhoi3, s1, s2, s3,  &
     h1, h2, h3, cp1, cp2, cp3, converg_3 )

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: compna
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                       :: t
  real(dp), intent(in)                       :: p
  real(dp), dimension(ncomp), intent(in)     :: x1
  real(dp), dimension(ncomp), intent(in)     :: x2
  real(dp), dimension(ncomp), intent(in)     :: x3
  real(dp), dimension(ncomp), intent(in)     :: rhoi1
  real(dp), dimension(ncomp), intent(in)     :: rhoi2
  real(dp), dimension(ncomp), intent(in)     :: rhoi3
  real(dp), intent(in)                       :: s1, s2, s3
  real(dp), intent(in)                       :: h1, h2, h3
  real(dp), intent(in), optional             :: cp1, cp2, cp3
  logical, intent(in)                        :: converg_3
  !-----------------------------------------------------------------------------
  integer                                    :: i
  real(dp)                                   :: mass_density1, mass_density2, mass_density3
  real(dp), dimension(ncomp)                 :: mass_frac1, mass_frac2, mass_frac3
  logical                                    :: is_open
  !-----------------------------------------------------------------------------

  call get_mass_density ( rhoi1, mass_density1, mass_frac1 )
  call get_mass_density ( rhoi2, mass_density2, mass_frac2 )
  if ( converg_3 ) call get_mass_density ( rhoi3, mass_density3, mass_frac3 )

  !-----------------------------------------------------------------------------
  ! open output file and write header
  !-----------------------------------------------------------------------------

  inquire( unit=40, opened=is_open )
  if ( .NOT.is_open ) open ( 40, FILE = './output_file/output.dat' )

  write (40,'(4a)') ' T/',unit_t,'  p/',unit_p
  write (40,'(2(2x,f18.8))') t-u_in_t, p/u_in_p

  if ( .NOT.converg_3 ) then
     write (40,'(6a)') ' substance  x1  x2  substance-name'
     do i = 1, ncomp
        write (40,'(2x,i5, 2(2x,f18.8), 2a)') i, x1(i), x2(i), '  ', compna(i)
     end do
     write (40,'(a)') ' density1/(kg/m**3)   density2/(kg/m**3)'
     write (40,'(2(2x,f18.8))') mass_density1, mass_density2
     write (40,'(a)') ' h1/(J/mol)  h2/(J/mol)'
     write (40,'(2(2x,f18.8))') h1, h2
     write (40,'(a)') ' s1/(J/mol/K)  s2/(J/mol/K)'
     write (40,'(2(2x,f18.8))') s1, s2
     write (40,'(a)') ' cpres1/(J/mol/K)  cpres2/(J/mol/K)'
     write (40,'(2(2x,f18.8))') cp1, cp2
  else
     write (40,'(6a)') ' substance  x1  x2  x3  substance-name'
     do i = 1, ncomp
        write (40,'(2x,i5, 2(2x,f18.8), a)') i, x1(i), x2(i), x3(i), compna(i)
     end do
     write (40,'(a)') ' density1/(kg/m**3)   density2/(kg/m**3)   density3/(kg/m**3)'
     write (40,'(2(2x,f18.8))') mass_density1, mass_density2, mass_density3
     write (40,'(a)') ' h1/(J/mol)  h2/(J/mol)  h3/(J/mol)'
     write (40,'(2(2x,f18.8))') h1, h2, h3
     write (40,'(a)') ' s1/(J/mol/K)  s2/(J/mol/K)  s3/(J/mol/K)'
     write (40,'(2(2x,f18.8))') s1, s2, s3
     write (40,'(a)') ' cpres1/(J/mol/K)  cpres2/(J/mol/K)  cpres3/(J/mol/K)'
     write (40,'(2(2x,f18.8))') cp1, cp2, cp3
  end if

  close (40)

  !-----------------------------------------------------------------------------
  ! output to terminal
  !-----------------------------------------------------------------------------
  write (*,'(a)') ' flash condition:'
  write (*,'(t2,a,f7.2,3a,f9.4,2a)') ' T =',t-u_in_T,' ',unit_t,'   p =',p/u_in_p,' ',unit_p
  write (*,'(a)') ' '
  write (*,'(a)') ' output is written to file: ./output_file/output.dat'

end subroutine output_flash


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine get_mass_density
!> \brief for vector rhoi of a phase, calc. mass-based density and mass fraction
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine get_mass_density ( rhoi, mass_density, mass_frac )

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mm
  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)     :: rhoi
  real(dp), intent(out)                      :: mass_density
  real(dp), dimension(ncomp), intent(out)    :: mass_frac
  !-----------------------------------------------------------------------------

  mass_density = sum( rhoi(1:ncomp) * mm(1:ncomp) )

  mass_frac(1:ncomp) = rhoi(1:ncomp) * mm(1:ncomp) / mass_density

  mass_density = mass_density * 1.E27 /NAV   ! convert units to  kg/m**3

end subroutine get_mass_density

end module input_output
