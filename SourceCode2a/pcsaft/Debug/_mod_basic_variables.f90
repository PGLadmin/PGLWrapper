!> \file mod_basic_variables.f90
!! \brief Several modules with all basic variables.
!!
!! \todo Detailed file description.

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module PARAMETERS
!> \brief Module with important global constants and parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module PARAMETERS

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
