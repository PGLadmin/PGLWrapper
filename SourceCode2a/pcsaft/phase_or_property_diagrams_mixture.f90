!> \file pure_par_fit.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module diagrams
!> \brief module for generating phase diagrams or property diagrams
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module module_mixture_diagrams

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  logical                                         :: iterate_t_tranfer

  private
  public :: state_property_calculation, binary_phase_diagram, binary_critical_points_diagram

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine binary_phase_diagram
!> \brief determine binary T-x or p-x phase diagram
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine binary_phase_diagram

  use BASIC_VARIABLES, only: t_input, p_input
  use properties, only: state_trho, state_tp
  use stability, only: get_eta
  use stability_G, only: stability_analysis_G
  use bubble_dew_point, only: bubble_point_rachford_rice
  use input_output, only: mixture_output_bin_coex
  use flash, only: Gibbs_flash
  use module_critical_point, only: critical_point_mixtures_binary_defined_t_or_p
  implicit none
  !-----------------------------------------------------------------------------
  integer                                         :: i
  integer                                         :: i_max
  integer                                         :: steps
  integer                                         :: iso_x
  real(dp)                                        :: tF, pF
  real(dp), dimension(ncomp)                      :: xF
  real(dp), dimension(ncomp)                      :: xV
  real(dp), dimension(ncomp)                      :: x_L1
  real(dp), dimension(ncomp)                      :: x_L2
  real(dp), dimension(2)                          :: eta_start
  real(dp), dimension(3)                          :: eta_start3
  real(dp), dimension(ncomp)                      :: rhoiL, rhoiV, rhoiL2
  real(dp), dimension(ncomp)                      :: rhoi_trial
  real(dp), dimension(ncomp)                      :: ln_zi_phi
  real(dp), dimension(ncomp)                      :: xc
  real(dp)                                        :: sV, sL, sL2
  real(dp)                                        :: hV, hL, hL2
  real(dp)                                        :: cpV, cpL, cpL2
  real(dp)                                        :: molar_rhoL, molar_rhoV
  real(dp)                                        :: pcalc
  real(dp)                                        :: x_max
  logical                                         :: iterate_t
  logical                                         :: converg
  logical                                         :: stable
  logical                                         :: get_start_val
  !-----------------------------------------------------------------------------

  steps = 80

  tF = t_input
  pF = p_input

  eta_start(1) = 0.4_dp
  eta_start(2) = 0.1E-4_dp

  i_max = 0
  x_max = 0.0_dp

  get_start_val = .true.

  !-----------------------------------------------------------------------------
  ! user input: calculate isothermal or isobaric diagram
  !-----------------------------------------------------------------------------
  write (*,*) ' '
  write (*,*) ' choose your calculation option for binary phase-diagram:'
  write (*,*) ' isotherm (1) or isobar (2)'
  read (*,*) iso_x
  if ( iso_x == 1 ) then
     iterate_t = .false.
  else
     iterate_t = .true.
  end if

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call mixture_output_bin_coex ( tF, pF, xF, xV, rhoiL, rhoiV, sL, sV, hL, hV, cp1=cpL, cp2=cpV, header=.true. )
  print*,'phase_or...: done mixture_output...'

  !-----------------------------------------------------------------------------
  ! perform bubble point calculation for given x (and vary x along a grid)
  !-----------------------------------------------------------------------------
  do i = 0, steps

     xF(1) = real( i, KIND=dp ) / real( steps, KIND=dp )
     if ( x_max > 0.0_dp ) xF(1) = x_max + real( i-i_max, KIND=dp ) / real( steps-i_max, KIND=dp ) * ( 1.0_dp - x_max )
     xF(2) = 1.0_dp - xF(1)

     if ( xF(1) < x_max ) CYCLE

     call bubble_point_rachford_rice ( iterate_t, tF, pF, eta_start, xF, xV, rhoiL, rhoiV, get_start_val, converg )

     if ( converg ) then

        get_start_val = .false.

        !-----------------------------------------------------------------------
        ! check stability of liquid phase
        !-----------------------------------------------------------------------
        if ( i /= 0 .AND. i /= steps ) then

           call stability_analysis_G ( tF, pF, xF, eta_start(1), stable, rhoiL, ln_zi_phi, rhoi_trial )

           if ( .NOT.stable ) then

              !x_L1(1:ncomp) = xF(1:ncomp)                                   ! initial estimate
              !x_L2(1:ncomp) = rhoi_trial(1:ncomp) / sum( rhoi_trial(1:ncomp) )  ! initial estimate from stability analysis
              !eta_start2(1) = get_eta( rhoiL )
              !eta_start2(2) = get_eta( rhoi_trial )
              !phase_frac2 = 0.001_dp

              !call Gibbs_flash ( tF, pF, xF, ln_zi_phi, eta_start2, phase_frac2, x_L1, x_L2, rhoiL, rhoiL2, converg_GF )

              eta_start3(1:2) = eta_start(1)
              eta_start3(3) = eta_start(2)
              rhoiL2 = rhoi_trial
              call binary_three_phase_equilibrium ( iterate_t, tF, pF, eta_start3, rhoiL, rhoiL2,  &
                                                    rhoiV, x_L1, x_L2, xV, converg )

              if ( converg ) then

                 x_max = max( x_L1(1), x_L2(1), xV(1) )
                 i_max = i

                 !--------------------------------------------------------------
                 ! output to file (and to terminal)
                 !--------------------------------------------------------------
                 call state_tp ( tF, pF, x_L1, eta_start3(1), rhoiL, molar_rhoL, sL, hL, cp=cpL )
                 call state_tp ( tF, pF, x_L2, eta_start3(2), rhoiL2, molar_rhoL, sL2, hL2, cp=cpL2 )
                 call state_tp ( tF, pF, xV, eta_start3(3), rhoiV, molar_rhoV, sV, hV, cp=cpV )
                 call mixture_output_bin_coex ( tF, pF, x_L1, xV, rhoiL, rhoiV, sL, sV, hL, hV,  &
                                                cp1=cpL, cp2=cpV, header=.false. )
                 call mixture_output_bin_coex ( tF, pF, x_L2, xV, rhoiL2, rhoiV, sL2, sV, hL2, hV,  &
                                                cp2=cpL2, cp1=cpV, header=.false. )
                 CYCLE

              end if

           end if

        end if

        !-----------------------------------------------------------------------
        ! calculate state point properties for liquid and vapor phase
        !-----------------------------------------------------------------------
        call state_trho ( tF, rhoiL, pcalc, molar_rhoL, sL, hL, cp=cpL )
        call state_trho ( tF, rhoiV, pcalc, molar_rhoV, sV, hV, cp=cpV )

        !-----------------------------------------------------------------------
        ! output to file (and to terminal)
        !-----------------------------------------------------------------------
        call mixture_output_bin_coex ( tF, pF, xF, xV, rhoiL, rhoiV, sL, sV, hL, hV,  &
                                       cp1=cpL, cp2=cpV, header=.false. )

     else

        get_start_val = .true.

        if ( i /= 0 .AND. i /= steps ) then

           ! call stability_analysis_G ( tF, pF, xF, eta_start(1), stable, rhoiL, ln_zi_phi, rhoi_trial )

           ! if ( .NOT. stable ) then
           !    write (*,*) 'bubble_point_rachford_rice should converge, but did not. Analyse this case'
           ! end if

           !---------------------------------------------------------------------
           ! calculate tF, pF, for given xF
           !---------------------------------------------------------------------
           xc = xF
           call critical_point_mixtures_binary_defined_t_or_p ( iterate_t, tF, pF, xc, rhoiL, converg )

           if ( converg ) then
              call state_trho ( tF, rhoiL, pcalc, molar_rhoL, sL, hL, cp=cpL )
              rhoiV = rhoiL
              call state_trho ( tF, rhoiV, pcalc, molar_rhoV, sV, hV, cp=cpV )
              write (*,*) 'calculated mixture critical point'
              call mixture_output_bin_coex ( tF, pF, xc, xc, rhoiL, rhoiV, sL, sV, hL, hV,  &
                                           cp1=cpL, cp2=cpV, header=.false. )
              exit
           end if

        end if

     end if

  end do
 
end subroutine binary_phase_diagram


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine three_phase_equilibrium
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine binary_three_phase_equilibrium ( iterate_t, tF, pF, eta_start3, rhoi1, rhoi2, rhoi3, x1, x2, x3, converg )

  use Solve_NonLin, only: hbrd
  implicit none
  !-----------------------------------------------------------------------------
  logical, intent(in)                             :: iterate_t
  real(dp), intent(inout)                         :: tF
  real(dp), intent(inout)                         :: pF
  real(dp), dimension(3), intent(in)              :: eta_start3
  real(dp), dimension(ncomp), intent(in)          :: rhoi1
  real(dp), dimension(ncomp), intent(in)          :: rhoi2
  real(dp), dimension(ncomp), intent(in)          :: rhoi3
  real(dp), dimension(ncomp), intent(out)         :: x1
  real(dp), dimension(ncomp), intent(out)         :: x2
  real(dp), dimension(ncomp), intent(out)         :: x3
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
!!$  interface
!!$     subroutine bin_3_phase_objec_fct ( ndim, y, residu, dummy )
!!$       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision       
!!$       integer, intent(in)                        :: ndim
!!$       real(dp), dimension(ndim), intent(in)      :: y
!!$       real(dp), dimension(ndim), intent(out)     :: residu
!!$       integer, intent(in out)                    :: dummy
!!$     end subroutine bin_3_phase_objec_fct
!!$  end interface

  integer                                         :: info
  integer                                         :: ndim
  real(dp), allocatable, dimension(:)             :: y
  real(dp), allocatable, dimension(:)             :: diag
  real(dp), allocatable, dimension(:)             :: residu
  real(dp)                                        :: tolerance
  real(dp)                                        :: step_size
  !-----------------------------------------------------------------------------

  info = 1

  ndim = 2*ncomp + 3   ! = 7
  allocate( y(ndim), diag(ndim), residu(ndim) )

  tolerance = 1.E-7_dp
  step_size = 1.E-9_dp

  converg = .false.

  !-----------------------------------------------------------------------------
  ! assign starting values for the iterated variables (to vector y)
  !-----------------------------------------------------------------------------

  y(1:ncomp) = log( rhoi1(1:ncomp) / sum( rhoi1(1:ncomp) ) )
  y(ncomp+1:2*ncomp) = log( rhoi2(1:ncomp) / sum( rhoi2(1:ncomp) ) )
  y(2*ncomp+1:3*ncomp) = log( rhoi3(1:ncomp) / sum( rhoi3(1:ncomp) ) )
  if ( iterate_t ) then
     y(3*ncomp+1) = log(tF)
     !p = pF
  else
     y(3*ncomp+1) = log(pF)
     !t = tF
  end if
  
  !-----------------------------------------------------------------------------
  ! solve the phase equilibrium conditions
  !-----------------------------------------------------------------------------

  iterate_t_tranfer = iterate_t             ! transfer iterate_t to the objective fct

  CALL hbrd (bin_3_phase_objec_fct, ndim, y, residu, step_size, tolerance, info, diag)


  !-----------------------------------------------------------------------------
  ! determine, whether the solution is acceptable (if yes: convergence = .true.)
  !-----------------------------------------------------------------------------

  if ( sum( abs( residu( 1:ndim ) ) ) < tolerance ) converg = .true.

  if ( converg ) then

     if ( iterate_t ) then
        tF = exp( y(3*ncomp+1) )
     else
        pF = exp( y(3*ncomp+1) )
     end if
     write (*,*) '3-phase-equilibrium converged',info, sum( abs( residu( 1:ndim ) ) )
     x1(1:ncomp) = exp( y(1:ncomp) )
     x2(1:ncomp) = exp( y(ncomp+1:2*ncomp) )
     x3(1:ncomp) = exp( y(2*ncomp+1:3*ncomp) )
     !------ normalize to reduce numerical noise -------------------------------
     x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
     x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )
     x3(1:ncomp) = x3(1:ncomp) / sum( x3(1:ncomp) )

  end if

  if ( .NOT.converg ) write (*,*) 'error 3-phase-equilibrium',info, sum( abs( residu( 1:ndim ) ) )

  deallocate( y, diag, residu )

end subroutine binary_three_phase_equilibrium


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE bin_3_phase_objec_fct
!
! equations solved for binary 3-phase equilibrium
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine bin_3_phase_objec_fct ( ndim, y, residu, dummy )

  use BASIC_VARIABLES, only: t_input, p_input
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in)           :: y
  real(dp), dimension(ndim), intent(out)          :: residu
  integer, intent(in out)                         :: dummy

  !-----------------------------------------------------------------------------
  real(dp)                                        :: tF
  real(dp)                                        :: pF
  real(dp), dimension(ncomp)                      :: x1
  real(dp), dimension(ncomp)                      :: x2
  real(dp), dimension(ncomp)                      :: x3
  real(dp), dimension(ncomp)                      :: rhoi1
  real(dp), dimension(ncomp)                      :: rhoi2
  real(dp), dimension(ncomp)                      :: rhoi3
  real(dp), dimension(ncomp)                      :: lnphi_1
  real(dp), dimension(ncomp)                      :: lnphi_2
  real(dp), dimension(ncomp)                      :: lnphi_3
  real(dp), dimension(3)                          :: eta_start3
  logical                                         :: iterate_t
  !-----------------------------------------------------------------------------

  iterate_t = iterate_t_tranfer

  !-----------------------------------------------------------------------------
  ! assign values for the iterated quantities
  !-----------------------------------------------------------------------------

  x1(1:ncomp) = exp( y(1:ncomp) )
  x2(1:ncomp) = exp( y(ncomp+1:2*ncomp) )
  x3(1:ncomp) = exp( y(2*ncomp+1:3*ncomp) )
  if ( iterate_t ) then
     tF = exp( y(3*ncomp+1) )
     pF = p_input
  else
     pF = exp( y(3*ncomp+1) )
     tF = t_input
  end if

  !-----------------------------------------------------------------------------
  ! determine first 3 values of residuum-vector
  !-----------------------------------------------------------------------------

  residu(1) = sum( x1(1:ncomp) ) - 1.0_dp
  residu(2) = sum( x2(1:ncomp) ) - 1.0_dp
  residu(3) = sum( x3(1:ncomp) ) - 1.0_dp

  !-----------------------------------------------------------------------------
  ! ensure summation relation
  !-----------------------------------------------------------------------------
  x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
  x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )
  x3(1:ncomp) = x3(1:ncomp) / sum( x3(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! standard starting values for densities
  !-----------------------------------------------------------------------------

  eta_start3(1) = 0.4_dp
  eta_start3(2) = 0.4_dp
  eta_start3(3) = 1.E-5_dp

  !-----------------------------------------------------------------------------
  ! calculate chemical potentials
  !-----------------------------------------------------------------------------

  call chemical_potential_tp ( tF, pF, x1, eta_start3(1), rhoi1(:), lnphi_1(:) )
  call chemical_potential_tp ( tF, pF, x2, eta_start3(2), rhoi2(:), lnphi_2(:) )
  call chemical_potential_tp ( tF, pF, x3, eta_start3(3), rhoi3(:), lnphi_3(:) )

  !-----------------------------------------------------------------------------
  ! determine remaining values of residuum-vector (isofugacity conditions)
  !-----------------------------------------------------------------------------

  residu(3+1:3+ncomp) = log( x1(1:ncomp) ) + lnphi_1(1:ncomp) - log( x3(1:ncomp) ) - lnphi_3(1:ncomp)
  residu(3+ncomp+1:3+2*ncomp) = log( x2(1:ncomp) ) + lnphi_2(1:ncomp) - log( x3(1:ncomp) ) - lnphi_3(1:ncomp)
  ! residu(3+1:3+ncomp) = y(1:ncomp) + lnphi_1(1:ncomp) - y(2*ncomp+1:3*ncomp) - lnphi_3(1:ncomp)
  ! residu(3+ncomp+1:3+2*ncomp) = y(ncomp+1:2*ncomp) + lnphi_2(1:ncomp) - y(2*ncomp+1:3*ncomp) - lnphi_3(1:ncomp)
  ! if ( dummy == 1 ) write (*,*) pF, sum( abs( residu ) )
  ! if ( dummy == 1 ) read (*,*)

end subroutine bin_3_phase_objec_fct



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine state_property_calculation
!> \brief determine state properties of mixtures to given input specifications
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine state_property_calculation

  use PARAMETERS, only: machine_eps
  use BASIC_VARIABLES, only: t_input, p_input, x_input
  use properties, only: state_tp
  use stability, only: get_eta
  use stability_G, only: stability_analysis_G
  use input_output, only: mixture_output_state
  implicit none
  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------

  tF = t_input
  pF = p_input
  xF( 1:ncomp ) = x_input( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call mixture_output_state ( tF, pF, xF, rhoiL, sL, hL, cp=cpL, header=.true. )

  !-----------------------------------------------------------------------------
  ! for liquid trial phase: test stability, then determine state variables
  !-----------------------------------------------------------------------------
  eta_start = 0.45_dp                        ! starting density for a liquid phase
  call stability_analysis_G ( tF, pF, xF, eta_start, stable_L, rhoiL, ln_zi_phi, rhoi_trial )
  eta_L = get_eta( rhoiL(:) )

  if ( get_eta( rhoiL(:) ) < 0.15_dp ) then

     write (*,*) 'LIQUID: A high-density liquid does not exist at given (T,p,x).'
     write (*,*) 'LIQUID: The fluid density (packing fraction, eta) is: ', eta_L
     write (*,*) ' '
     if ( .NOT.stable_L ) write (*,*) 'LIQUID: trial phase with liquid-density as starting value is NOT stable !!!'

  else

     if ( .NOT.stable_L ) then
        write (*,*) 'LIQUID: liquid phase is NOT stable !!!'
     else
        write (*,*) 'LIQUID: liquid phase is stable.'
     end if
     write (*,*) ' '

  end if

  if ( stable_L ) call state_tp ( tF, pF, xF, eta_L, rhoiL, molar_rho, sL, hL, cp=cpL )

  !-----------------------------------------------------------------------------
  ! for vapor trial phase: test stability, then determine state variables
  !-----------------------------------------------------------------------------
  eta_start = 1.E-5_dp                      ! starting density for a vapor phase
  call stability_analysis_G ( tF, pF, xF, eta_start, stable_V, rhoiV, ln_zi_phi, rhoi_trial )
  eta_V = get_eta( rhoiV(:) )

  if ( get_eta( rhoiV(:) ) > 0.2_dp ) then

     write (*,*) 'VAPOR:  A low-density vapor does not exist at given (T,p,x).'
     write (*,*) 'VAPOR:  The fluid density (packing fraction, eta) is: ', eta_V
     write (*,*) ' '
     if ( .NOT.stable_V ) write (*,*) 'VAPOR:  trial phase with vapor-density as starting value is NOT stable !!!'

  else

     if ( .NOT.stable_V ) then
        write (*,*) 'VAPOR:  vapor phase is NOT stable !!!'
     else
        write (*,*) 'VAPOR:  vapor phase is stable.'
     end if
     write (*,*) ' '

  end if

  if ( stable_V ) call state_tp ( tF, pF, xF, eta_V, rhoiV, molar_rho, sV, hV, cp=cpV )

  
  !-----------------------------------------------------------------------------
  ! write output
  !-----------------------------------------------------------------------------
  if ( stable_L .AND. stable_V ) then

     !--------------------------------------------------------------------------
     ! this case: only one solution is found for two initial values of eta_start
     !--------------------------------------------------------------------------
     if ( abs( eta_V/eta_L - 1.0_dp ) > 1.E3_dp * machine_eps ) then
        write (*,*) 'state_property_calculation: unexpected case, stop'
        stop
     else
        call mixture_output_state ( tF, pF, xF, rhoiL, sL, hL, cp=cpL, header=.false. )
     end if
     
     !gL = hL - tF * sL
     !gV = hV - tF * sV

  else if ( stable_L ) then

     !--------------------------------------------------------------------------
     ! liquid is stable: output to file (and to terminal)
     !--------------------------------------------------------------------------
     call mixture_output_state ( tF, pF, xF, rhoiL, sL, hL, cp=cpL, header=.false. )

  else if ( stable_V ) then

     !--------------------------------------------------------------------------
     ! vapor is stable: output to file (and to terminal)
     !--------------------------------------------------------------------------
     call mixture_output_state ( tF, pF, xF, rhoiV, sV, hV, cp=cpV, header=.false. )

  else

     write (*,*) ' '
     write (*,*) 'FINAL:  the system is not not stable at given (T,p,x). Perform a flash calculation.'

  end if

end subroutine state_property_calculation

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine binary_critical_points_diagram
!> \brief determine critical points of mixtures
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine binary_critical_points_diagram

  use BASIC_VARIABLES, only: t_input, p_input, x_input
  use properties, only: state_trho
  use input_output, only: mixture_output_state
  use module_critical_point, only: critical_point_mixtures
  implicit none

  !-----------------------------------------------------------------------------
  real(dp)                                        :: tF, pF
  real(dp), dimension(ncomp)                      :: xF
  real(dp), dimension(ncomp)                      :: rhoi
  real(dp)                                        :: sc, hc
  real(dp)                                        :: molar_rho
  real(dp)                                        :: pcalc
  logical                                         :: converg
  !-----------------------------------------------------------------------------

  tF = t_input
  pF = p_input
  xF( 1:ncomp ) = x_input( 1:ncomp )

  !-----------------------------------------------------------------------------
  ! open the output file and write the header of output file
  !-----------------------------------------------------------------------------
  call mixture_output_state ( tF, pF, xF, rhoi, sc, hc, header=.true. )

  !-----------------------------------------------------------------------------
  ! for liquid trial phase: test stability, then determine state variables
  !-----------------------------------------------------------------------------
  call critical_point_mixtures ( xF, tF, pF, rhoi, converg )

  if ( converg ) call state_trho ( tF, rhoi, pcalc, molar_rho, sc, hc )

  
  !-----------------------------------------------------------------------------
  ! write output
  !-----------------------------------------------------------------------------
  if ( converg ) then

     call mixture_output_state ( tF, pF, xF, rhoi, sc, hc, header=.false. )

  else

     write (*,*) 'critical point calculation did not converg for T_initial=', t_input

  end if

end subroutine binary_critical_points_diagram

end module module_mixture_diagrams
