module system_definition

  implicit none

  private
  public :: build_problem

contains

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine build_system
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine build_problem

  use input_output, only: read_problem_definition
  use eos_constants
  use module_eos_derivatives, only: allocate_eos_quantities, initialize_eos_quantities,  &
       cross_association_parameters, aggregate_polar_parameters
  use pcsaft_pure_and_binary_parameters, only: load_pure_and_binary_parameters
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! read definition of problem
  !-----------------------------------------------------------------------------
  call read_problem_definition

  !-----------------------------------------------------------------------------
  ! allocate and initialize array quantities
  !-----------------------------------------------------------------------------
  call allocate_eos_quantities
  call initialize_eos_quantities

  !-----------------------------------------------------------------------------
  ! read pure component and mixture parameters of PC-SAFT EOS
  !-----------------------------------------------------------------------------
  call load_pure_and_binary_parameters
  call cross_association_parameters
  call aggregate_polar_parameters

  !-----------------------------------------------------------------------------
  ! EOS constants
  !-----------------------------------------------------------------------------
  call load_eos_constants

end subroutine build_problem


end module system_definition
