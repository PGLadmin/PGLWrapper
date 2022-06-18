!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! program PC_SAFT
!
! 1) build a problem (i.e. define components and T and p)
! 2) perform property calculation
! 3) deallocate quantities
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

program PC_SAFT
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ		!JRE
	USE PORTLIB																		  !JRE

  !-----------------------------------------------------------------------------
  use system_definition, only: build_problem
  use module_eos_derivatives, only: deallocate_eos_quantities
  IMPLICIT NONE
  !-----------------------------------------------------------------------------
	  integer(2) iStat
	  doublePrecision tKelvin
	  integer ioErr
	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	!iChange=VERIFY('C:\MYPROJEX\CalcEos',CURDIR)
	!iChange2=VERIFY('C:\MYPROJEX\CALCEOS',CURDIR)
	!if(iChange2.eq.0)iChange=0
	!if(iChange.eq.0 .or. iChange.gt.26)DEBUG=.TRUE.
	!outFile=TRIM(masterDir)//'\output.txt' ! // is the concatenation operator
	!OPEN(UNIT=52,FILE=outFile,RECL=8192)
	print*,'Main: curDir=',TRIM(curDir)
    open(UNIT=51,FILE='input.txt',ioStat=ioErr)
	if(ioErr)print*,'Main: failed to open input.txt'
	read(51,*,ioStat=ioErr)tKelvin
	if(ioErr)print*,'Main: failed to read first line of input.txt'
	if(ioErr==0)print*,'T(K)',tKelvin
	close(51)

  call build_problem

  call choose_calculation_option

  call deallocate_eos_quantities

end program PC_SAFT



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine choose_calculation_option
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine choose_calculation_option

  use BASIC_VARIABLES, only: ncomp, t_input, p_input, x_input
  use pcsaft_pure_and_binary_parameters, only: compna
  use module_mixture_diagrams, only: state_property_calculation
  use flash, only: flash_with_output
  !JRE 20200630 Sorry, this feature does not compile properly.
  !use pure_fit_parameters, only: fitting_pure_component_parameters_driver
  use kij_fitting, only: kij_lij_fitting
  implicit none

  !-----------------------------------------------------------------------------
  integer                                         :: i
  integer                                         :: option
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! write system info to terminal
  !-----------------------------------------------------------------------------
  write (*,'(a)') ' '
  write (*,'(a)') ' This is a program for calculating physical properties'
  write (*,'(a)') ' and phase equilibria using the PC-SAFT equation of state'
  write (*,'(a)') ' '
  write (*,'(a)') ' The program was written by'
  write (*,'(a)') ' '
  write (*,'(a)') ' Joachim Gross'
  write (*,'(a)') ' University of Stuttgart'
  write (*,'(a)') ' Institute for Thermodynamics'
  write (*,'(a)') ' and Thermal Process Engineering'
  write (*,'(a)') ' gross@itt.uni-stuttgart.de'
  write (*,'(a)') ' '
  write (*,'(a)') ' '
  write (*,'(a)') ' The input file (./input_file/input.dat) currently'
  write (*,'(a)') ' specifies the following substances:'
  write (*,'(a)') ' -----------------------------------'
  do i = 1, ncomp
     write (*,'(a,i2,2a)') ' component ',i,': ',compna(i)
  end do
  write (*,'(a)') ' '
  write (*,'(a)') ' Choose the calculation option: '
  write (*,'(a)') ' -----------------------------------'
  write (*,'(a)') ' property diagram or phase diagram            (1)'
  if ( ncomp > 1) write (*,'(a)') ' flash calculation                            (2)'
  write (*,'(a)') ' calculate properties at given state point    (3)'
  write (*,'(a)') ' adjusting pure component parameters          (4)'
  write (*,'(a)') ' adjusting binary parameters kij              (5)'
  write (*,'(a)') ' '
  read (*,*) option

  select case ( option )
  case(1)
     call select_prop_or_phase_diagram
  case(2)
     call flash_with_output ( t_input, p_input, x_input )
  case(3)
     call state_property_calculation
  case(4)
	print*,'Main: Sorry, this feature does not compile properly. JRE 20200630'
     !call fitting_pure_component_parameters_driver
  case(5)
     call kij_lij_fitting
  case DEFAULT
     write(*,*) ' '
     write(*,*) ' Your choice is not valid!'
     write(*,*) ' '
     stop
  end select

end subroutine choose_calculation_option


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine select_phase_diagram
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine select_prop_or_phase_diagram

  use BASIC_VARIABLES, only: ncomp
  use module_pure_diagrams, only: pure_diagrams
  use module_mixture_diagrams, only: binary_phase_diagram, binary_critical_points_diagram
  implicit none

  integer                               :: option
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! offer various diagram-options, depending on the number of species
  !-----------------------------------------------------------------------------
  if ( ncomp == 1 ) then

     call pure_diagrams

  else if ( ncomp == 2 .OR. ncomp == 3 ) then

     write (*,'(a)') ' select type of property diagram or phase diagram'
     if ( ncomp == 2 ) then
        write (*,'(a)') ' binary phase equilibrium diagram             (1)'
        !write (*,'(a)') ' property-x diagrams (at (T,p)=const.)        (2)'
        !write (*,'(a)') ' binary VLLE equilibrium                      (3)'
        write (*,'(a)') ' binary critical points diagram               (4)'
        !write (*,'(a)') ' binary Joule-Thomson inversion points        (5)'
        !write (*,'(a)') ' p-T diagram (at x=const.)                    (6)'
        read (*,*) option
        select case (option)
        case(1)
           call binary_phase_diagram
        !case(2)
        !   call property_x_diagram
        !case(3)
        !   call binary_VLLE
        case(4)
           call binary_critical_points_diagram
        !case(5)
        !   call Joule_Thomson_inversion_points
        !case(6)
        !   call p_t_diagram
        case DEFAULT
           write(*,*) ' '
           write(*,*) ' Your choice is not valid!'
           stop
        end select
     else
        write (*,'(a)') ' ternary diagram                              (1)'
        write (*,'(a)') ' ternary 3-phase points                       (2)'
        write (*,'(a)') ' p-T diagram (at x=const.)                    (3)'
        read (*,*) option
        select case (option)
        !case(1) call ternary_phase_diagram
        !case(2) call ternary_3_phase_equilibrium
        !case(3) call p_t_diagram
        case DEFAULT
           write(*,*) ' '
           write(*,*) ' Your choice is not valid!'
           stop
        end select
     end if

  else

     write (*,'(a)') ' you are calculating a p-T diagram (at x=const.)'
     !call p_t_diagram

  end if

end subroutine select_prop_or_phase_diagram
