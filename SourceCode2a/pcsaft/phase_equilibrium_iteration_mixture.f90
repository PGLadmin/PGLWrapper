!> \file pure_par_fit.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module diagrams
!> \brief module for generating phase diagrams or property diagrams
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module module_iterate_pe

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                              :: outp = 0

  ! integer, parameter                              :: str3 = 3
  ! character(LEN=str3)                             :: specType
  ! real(dp)                                        :: specValue
  real(dp), dimension(2)                          :: eta_start_transfer
  real(dp), dimension(:), allocatable             :: rhoi1_transfer
  real(dp), dimension(:), allocatable             :: rhoi2_transfer

  private
  public :: iterate_phase_equilibrium

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine binary_phase_diagram
!> \brief determine binary T-x or p-x phase diagram
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine iterate_phase_equilibrium ( t, p, eta_start, xi1, xi2, rhoi1, rhoi2, converg )  ! ( spec_type, spec_value, t, p, eta_start, xi1, xi2, converg )

  use BASIC_VARIABLES, only: t_input, p_input
  ! use properties, only: state_trho, state_tp
  ! use stability, only: get_eta
  use stability_G, only: stability_analysis_G
  ! use bubble_dew_point, only: bubble_point_rachford_rice
  ! use input_output, only: mixture_output_bin_coex
  ! use flash, only: Gibbs_flash
  use Solve_NonLin, only: hbrd
  !-----------------------------------------------------------------------------
  ! character(LEN=str3), intent(in)                 :: spec_type
  ! real(dp), intent(in)                            :: spec_value
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), dimension(ncomp), intent(inout)       :: xi1
  real(dp), dimension(ncomp), intent(inout)       :: xi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer                                         :: info
  integer                                         :: ndim
  real(dp), allocatable, dimension(:)             :: y
  real(dp), allocatable, dimension(:)             :: diag
  real(dp), allocatable, dimension(:)             :: residu
  real(dp)                                        :: tolerance
  real(dp)                                        :: step_size
  real(dp), dimension(ncomp)                      :: xi1_initial, xi2_initial
  logical                                         :: trivial
  !-----------------------------------------------------------------------------

  info = 1

  ndim = 2 * ncomp
  allocate( y(ndim), diag(ndim), residu(ndim) )

  allocate( rhoi1_transfer(ncomp), rhoi2_transfer(ncomp) )

  tolerance = 1.E-7_dp
  step_size = 1.E-9_dp

  converg = .false.

  xi1_initial(:) = xi1(:)
  xi2_initial(:) = xi2(:)
  !-----------------------------------------------------------------------------
  ! assign starting values for the iterated variables (to vector y)
  !-----------------------------------------------------------------------------

  y(1:ncomp) = log( xi1(1:ncomp) )
  y(ncomp+1:2*ncomp) = log( xi2(1:ncomp) )
  ! y(2*ncomp+1:3*ncomp) = lnK(1:ncomp)

  eta_start_transfer(1:2) = eta_start(1:2)
  t_input = t
  p_input = p

  !-----------------------------------------------------------------------------
  ! solve the phase equilibrium conditions
  !-----------------------------------------------------------------------------

  CALL hbrd (phase_equil_objec_fct, ndim, y, residu, step_size, tolerance, info, diag)

  !-----------------------------------------------------------------------------
  ! determine, whether the solution is acceptable (if yes: convergence = .true.)
  !-----------------------------------------------------------------------------
  call check_trivial_solution ( rhoi1_transfer, rhoi2_transfer, trivial )

  if ( sum( abs( residu( 1:ndim ) ) ) < tolerance .AND. .NOT.trivial ) converg = .true.

  if ( converg ) then
     if ( outp >= 1 ) write (*,*) 'equilibrium converged',info, sum( abs( residu( 1:ndim ) ) )
     xi1(1:ncomp) = exp( y(1:ncomp) )
     xi2(1:ncomp) = exp( y(ncomp+1:2*ncomp) )
     rhoi1(:) = rhoi1_transfer(:)
     rhoi2(:) = rhoi2_transfer(:)
  else
     if ( outp >= 1 ) write (*,*) 'iterate_phase_equilibrium: no convergence',info, sum( abs( residu( 1:ndim ) ) )
     xi1(:) = xi1_initial(:)
     xi2(:) = xi2_initial(:)
     rhoi1(:) = 0.0_dp
     rhoi2(:) = 0.0_dp
  end if

  deallocate( y, diag, residu )
  deallocate( rhoi1_transfer, rhoi2_transfer )

end subroutine iterate_phase_equilibrium



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE check_trivial_solution
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_trivial_solution ( rhoi1, rhoi2, trivial )

  use stability, only: get_eta
  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi1
  real(dp), dimension(ncomp), intent(in)          :: rhoi2
  logical, intent(out)                            :: trivial
  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: eta_measure
  real(dp), dimension(ncomp)                      :: diff_vector
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     diff_vector(i) = abs( rhoi1(i) - rhoi2(i) )
  end do

  eta_measure = get_eta( diff_vector )

  if ( eta_measure < 1.E-5 ) then
     trivial = .true.
  else
     trivial = .false.
  end if

end subroutine check_trivial_solution


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE bin_3_phase_objec_fct
!
! equations solved for binary 3-phase equilibrium
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine phase_equil_objec_fct ( ndim, y, residu, dummy )

  use BASIC_VARIABLES, only: t_input, p_input
  use properties, only: chemical_potential_tp

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in)           :: y
  real(dp), dimension(ndim), intent(out)          :: residu
  integer, intent(in out)                         :: dummy

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: tF
  real(dp)                                        :: pF
  real(dp), dimension(ncomp)                      :: xi1
  real(dp), dimension(ncomp)                      :: xi2
  ! real(dp), dimension(ncomp)                      :: lnK
  real(dp), dimension(ncomp)                      :: rhoi1
  real(dp), dimension(ncomp)                      :: rhoi2
  real(dp), dimension(ncomp)                      :: lnphi_1
  real(dp), dimension(ncomp)                      :: lnphi_2
  real(dp), dimension(2)                          :: eta_start
  ! character(LEN=2)                                :: comp_no
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! assign values for the iterated quantities
  !-----------------------------------------------------------------------------

  if ( maxval( y(1:ncomp) ) < 200.0_dp ) then
     do i = 1, ncomp
        xi1(i) = 0.0_dp
        if ( abs( y(i) ) < 200.0_dp ) then
           xi1(i) = exp( y(i) )
        end if
        xi2(i) = 0.0_dp
        if ( abs( y(ncomp+i) ) < 200.0_dp ) then
           xi2(i) = exp( y(ncomp+i) )
        end if
     end do
  else
     xi1(:) = 0.0_dp
     xi1( maxloc( y(1:ncomp) ) ) = 1.0_dp
  end if

  if ( maxval( y(ncomp+1:2*ncomp) ) < 200.0_dp ) then
     do i = 1, ncomp
        xi2(i) = 0.0_dp
        if ( abs( y(ncomp+i) ) < 200.0_dp ) then
           xi2(i) = exp( y(ncomp+i) )
        end if
     end do
  else
     xi2(:) = 0.0_dp
     xi2( maxloc( y(ncomp+1:2*ncomp) ) ) = 1.0_dp
  end if

  ! lnK(1:ncomp) = y(2*ncomp+1:3*ncomp)

  tF = t_input
  pF = p_input
  
  !-----------------------------------------------------------------------------
  ! determine first 3 values of residuum-vector
  !-----------------------------------------------------------------------------

  residu(ndim-1) = sum( xi1(1:ncomp) ) - 1.0_dp
  residu(ndim  ) = sum( xi2(1:ncomp) ) - 1.0_dp

  !-----------------------------------------------------------------------------
  ! ensure summation relation
  !-----------------------------------------------------------------------------
  xi1(1:ncomp) = xi1(1:ncomp) / sum( xi1(1:ncomp) )
  xi2(1:ncomp) = xi2(1:ncomp) / sum( xi2(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! starting values for densities
  !-----------------------------------------------------------------------------

  eta_start(1:2) = eta_start_transfer(1:2)

  !-----------------------------------------------------------------------------
  ! calculate chemical potentials
  !-----------------------------------------------------------------------------

  call chemical_potential_tp ( tF, pF, xi1, eta_start(1), rhoi1(:), lnphi_1(:) )
  call chemical_potential_tp ( tF, pF, xi2, eta_start(2), rhoi2(:), lnphi_2(:) )
  rhoi1_transfer(:) = rhoi1(:)
  rhoi2_transfer(:) = rhoi2(:)

  !-----------------------------------------------------------------------------
  ! determine remaining values of residuum-vector (isofugacity conditions)
  !-----------------------------------------------------------------------------

  residu(1:ncomp) = y(1:ncomp) + lnphi_1(1:ncomp) - y(ncomp+1:2*ncomp) - lnphi_2(1:ncomp)
  ! residu(1:ncomp) = lnK(1:ncomp) + lnphi_1(1:ncomp) - lnphi_2(1:ncomp)
  ! residu(ncomp+1:2*ncomp) = xi2(1:ncomp) - exp( lnK(1:ncomp) ) * xi1(1:ncomp)

  ! if ( dummy == 1 ) write (*,*) residu
  ! if ( dummy == 1 ) read (*,*)

end subroutine phase_equil_objec_fct

end module module_iterate_pe
