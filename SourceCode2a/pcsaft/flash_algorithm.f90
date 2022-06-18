!> \file starting_value.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module flash

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                    :: outp = 0

  !-----------------------------------------------------------------------------
  ! variables that are not passed through the argument list to a Newton solver.
  !-----------------------------------------------------------------------------
  real(dp)                              :: t_flash
  real(dp)                              :: p_flash
  real(dp), allocatable, dimension(:)   :: zi_flash
  real(dp), allocatable, dimension(:)   :: ln_zi_phi_flash
  real(dp), dimension(2)                :: dense_start_flash

  private
  public :: flash_driver, G_flash_grads, G_tpd_objective, G_flash_parameter_bound,  &
       pt_flash, ps_flash, ph_flash, Gibbs_flash, rachford_rice, flash_with_output

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine flash_with_output
!> \brief Flash calculation for specified t, p, composition xiF(:).
!!
!! This subroutine performs a flash calculation for a mixtue with defined t, p,
!! and composition xiF. Output is written to a file: ./output_file/output.dat.
!!
!! input:    t, p, xiF(:)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine flash_with_output ( t, p, xiF )

  use input_output, only: output_flash
  use properties, only: state_trho
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: rhoi1
  real(dp), dimension(ncomp)                      :: rhoi2
  real(dp), dimension(ncomp)                      :: rhoi3
  real(dp), dimension(ncomp)                      :: x1
  real(dp), dimension(ncomp)                      :: x2
  real(dp), dimension(ncomp)                      :: x3
  real(dp)                                        :: h1, h2, h3
  real(dp)                                        :: s1, s2, s3
  real(dp)                                        :: cp1, cp2, cp3
  real(dp)                                        :: v
  real(dp)                                        :: pcalc
  logical                                         :: converg
  logical                                         :: converg_3
  !-----------------------------------------------------------------------------

  call flash_driver ( t, p, xiF, rhoi1, rhoi2, rhoi3, converg, converg_3 )

  if ( converg ) then

     x1(:) = rhoi1(:) / sum( rhoi1(1:ncomp) )
     x2(:) = rhoi2(:) / sum( rhoi2(1:ncomp) )
     call state_trho ( t, rhoi1, pcalc, v, s1, h1, cp=cp1 )
     call state_trho ( t, rhoi2, pcalc, v, s2, h2, cp=cp2 )

     if ( converg_3 ) then
        x3(:) = rhoi3(:) / sum( rhoi3(1:ncomp) )
        call state_trho ( t, rhoi3, pcalc, v, s3, h3, cp=cp3 )
     end if

     call output_flash ( t, p, x1, x2, x3, rhoi1, rhoi2, rhoi3, s1, s2, s3,  &
          h1, h2, h3, cp1, cp2, cp3, converg_3 )

  else

     write (*,*) 'the fluid is (probably) stable at these conditions'

  end if

end subroutine flash_with_output



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine flash_driver
!> \brief Flash calculation for specified t, p, composition xiF(:).
!!
!! This subroutine performs a flash calculation for a mixtue with defined t, p,
!! and composition xiF. A starting value for a liquid density is assumed !!
!! A stability analysis is performed for identifying a possible third phase
!! or another stable second phase. And a subsequent (two or three) phase equi-
!! librium calculation follows.
!!
!! input:    t, p, xiF(:)
!! output:   rhoi1, rhoi2, rhoi3 as converged component densities, if solution is found
!!           converg is .true. if a solution is found
!!           converg_3 is .true. if 3-phase solution is found
!!           rhoi3 is meaningful only if (converg_3 = .true.)
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine flash_driver ( t, p, xiF, rhoi1, rhoi2, rhoi3, converg, converg_3 )

  use stability, only: stability_analysis, get_eta
  use stability_G, only: stability_analysis_G
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi3
  logical, intent(out)                            :: converg
  logical, intent(out)                            :: converg_3

  !-----------------------------------------------------------------------------
  integer                                         :: nphase
  real(dp), dimension(ncomp)                      :: rhoi_out
  logical                                         :: stable
  real(dp)                                        :: phase_frac2
  real(dp)                                        :: etaF
  real(dp)                                        :: dense_start1
  real(dp), dimension(2)                          :: dense_start2
  real(dp), dimension(3)                          :: dense_start3
  real(dp), dimension(ncomp)                      :: ln_zi_phi
  real(dp), dimension(ncomp)                      :: zi
  real(dp), dimension(ncomp)                      :: x1, x2
  real(dp), dimension(ncomp,3)                    :: xi
  real(dp), dimension(ncomp)                      :: rhoiA
  real(dp), dimension(ncomp)                      :: rhoiB
  real(dp), dimension(ncomp)                      :: rhoiF
  real(dp), dimension(3)                          :: rho3
  logical                                         :: slow_iteration
  logical                                         :: converg_GF
  !-----------------------------------------------------------------------------

  converg_3 = .false.
  rhoi3(:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! start with a Rachford-Rice iteration that searches for VLE
  !-----------------------------------------------------------------------------

  call rachford_rice ( t, p, xiF, rhoi1, rhoi2, phase_frac2, slow_iteration, converg )

  !-----------------------------------------------------------------------------
  ! determine for feed:  ln_zi_phi(:) = log( rhoiF(:) ) + lnphi_F(:)  of feed
  ! that quantity is needed by both: stability analysis and alternative flash SRs
  !-----------------------------------------------------------------------------
  call determine_ln_zi_phi ( t, p, xiF, ln_zi_phi, etaF )

  !-----------------------------------------------------------------------------
  ! if Rachford-Rice did not converge but phase split will occur, proceed.
  !-----------------------------------------------------------------------------

  if ( .NOT.converg .AND. slow_iteration ) then
     dense_start2(1) = 0.4_dp
     dense_start2(2) = 1.E-5_dp
     x2(1:ncomp) = rhoi2(1:ncomp) / sum( rhoi2(1:ncomp) )
     call Gibbs_flash ( t, p, xiF, ln_zi_phi, dense_start2, phase_frac2, x1, x2, rhoi1, rhoi2, converg )
  end if

  !-----------------------------------------------------------------------------
  ! if flash calc. converged, perform stability analysis for liquid phase. If flash
  ! calc. did not converge, perform stability analysis for feed composition
  !-----------------------------------------------------------------------------

  if ( converg ) then

     x1(:) = rhoi1(:) / sum( rhoi1(1:ncomp) )
     dense_start1 = get_eta( rhoi1 )  ! assuming phase 1 from RR-algorithm is a liquid phase
     call stability_analysis_G ( t, p, x1, dense_start1, stable, rhoiF, ln_zi_phi, rhoi_out, rhoi2 )
     ! rhoiF is identical to rhoi1
     zi(1:ncomp) = x1(1:ncomp)

  else

     dense_start1 = etaF
     ! call stability_analysis ( t, p, xiF, dense_start1, stable, rhoiF, ln_zi_phi, rhoi_out )
     call stability_analysis_G ( t, p, xiF, dense_start1, stable, rhoiF, ln_zi_phi, rhoi_out )
     zi(1:ncomp) = xiF(1:ncomp)

  end if
  !-----------------------------------------------------------------------------
  ! note: zi(:) corresponds to rhoiF(:), according to zi(:) = rhoiF(:) / sum( rhoiF(:) )
  !-----------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------
  ! if either feed condition is ustable, or coexisting liquid (obtained from RR)
  ! is unstable, determine other coexisting phase
  !-----------------------------------------------------------------------------
  if ( .NOT.stable ) then

     x1(1:ncomp) = zi(1:ncomp)                                   ! initial estimate
     x2(1:ncomp) = rhoi_out(1:ncomp) / sum( rhoi_out(1:ncomp) )  ! initial estimate from stability analysis
     dense_start2(1) = get_eta( rhoiF )
     dense_start2(2) = get_eta( rhoi_out )
     phase_frac2 = 0.01_dp
     call Gibbs_flash ( t, p, zi, ln_zi_phi, dense_start2, phase_frac2, x1, x2, rhoiA, rhoiB, converg_GF )

     if ( converg_GF .AND. converg ) then

        dense_start3(1) = get_eta( rhoiA )
        dense_start3(2) = get_eta( rhoiB )
        dense_start3(3) = get_eta( rhoi2 )
        xi(:,1) = x1(:)
        xi(:,2) = x2(:)
        xi(:,3) = rhoi2(:) / sum( rhoi2(:) )

        nphase = 3
        call three_phase_flash ( nphase, t, p, xiF, zi, dense_start3, xi, rho3, converg_3 )

        if ( converg_3 ) then
           rhoi1(:) = rho3(1) * xi(:,1)
           rhoi2(:) = rho3(2) * xi(:,2)
           rhoi3(:) = rho3(3) * xi(:,3)
           converg = .true.
           write (*,*) 'three phase equilibrium converged'
        else
           rhoi1(:) = rhoiA(:)
           rhoi2(:) = rhoiB(:)
           x1(:) = rhoiA(:) / sum( rhoiA(1:ncomp) )
           x2(:) = rhoiB(:) / sum( rhoiB(1:ncomp) )
           converg = .true.
        end if

     else if ( converg_GF .AND. .NOT.converg ) then
        rhoi1(:) = rhoiA(:)
        rhoi2(:) = rhoiB(:)
        x1(:) = rhoiA(:) / sum( rhoiA(1:ncomp) )
        x2(:) = rhoiB(:) / sum( rhoiB(1:ncomp) )
        converg = .true.
     end if


  end if

end subroutine flash_driver

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine pt_flash
!> \brief Flash calculation for specified t, p, composition xiF(:).
!!
!! This subroutine performs a flash calculation for a mixtue with defined t, p,
!! and composition xiF.
!! A stability analysis for a possible third phase or another stable second
!! phase is not performed
!!
!! input:    t, p, xiF(:)
!! output:   converg
!!           rhoi1, rhoi2 as converged component densities, if solution is found
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pt_flash ( t, p, xiF, rhoi1, rhoi2, phase_frac2, converg )

  use stability, only: get_eta
  use stability_G, only: stability_analysis_G
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  real(dp), intent(out)                           :: phase_frac2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: rhoi_out
  logical                                         :: stable
  real(dp)                                        :: etaF
  real(dp)                                        :: dense_start1
  real(dp), dimension(2)                          :: dense_start2
  real(dp), dimension(ncomp)                      :: ln_zi_phi
  real(dp), dimension(ncomp)                      :: zi
  real(dp), dimension(ncomp)                      :: x1, x2
  real(dp), dimension(ncomp)                      :: rhoiA
  real(dp), dimension(ncomp)                      :: rhoiB
  real(dp), dimension(ncomp)                      :: rhoiF
  logical                                         :: converg_GF
  logical                                         :: slow_iteration
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! start with a Rachford-Rice iteration that searches for VLE
  !-----------------------------------------------------------------------------

  call rachford_rice ( t, p, xiF, rhoi1, rhoi2, phase_frac2, slow_iteration, converg )

  !-----------------------------------------------------------------------------
  ! determine for feed:  ln_zi_phi(:) = log( rhoiF(:) ) + lnphi_F(:)  of feed
  ! that quantity is needed by both: stability analysis and alternative flash SRs
  !-----------------------------------------------------------------------------
  call determine_ln_zi_phi ( t, p, xiF, ln_zi_phi, etaF )

  !-----------------------------------------------------------------------------
  ! if Rachford-Rice did not converge but phase split will occur, proceed.
  !-----------------------------------------------------------------------------

  if ( .NOT.converg .AND. slow_iteration ) then
     dense_start2(1) = 0.4_dp
     dense_start2(2) = 1.E-5_dp
     x2(1:ncomp) = rhoi2(1:ncomp) / sum( rhoi2(1:ncomp) )
     call Gibbs_flash ( t, p, xiF, ln_zi_phi, dense_start2, phase_frac2, x1, x2, rhoi1, rhoi2, converg )
  end if

  !-----------------------------------------------------------------------------
  ! if RR converged, perform stability analysis for liquid phase. If RR did not
  ! converge, perform stability analysis for feed composition
  !-----------------------------------------------------------------------------

  if ( .NOT.converg ) then

     dense_start1 = etaF
     call stability_analysis_G ( t, p, xiF, dense_start1, stable, rhoiF, ln_zi_phi, rhoi_out )
     zi(1:ncomp) = xiF(1:ncomp)
     !-----------------------------------------------------------------------------
     ! note: zi(:) corresponds to rhoiF(:), according to zi(:) = rhoiF(:) / sum( rhoiF(:) )
     !-----------------------------------------------------------------------------

  
     !-----------------------------------------------------------------------------
     ! if either feed condition is ustable, or coexisting liquid (obtained from RR)
     ! is unstable, determine other coexisting phase
     !-----------------------------------------------------------------------------
     if ( .NOT.stable ) then

        x1(1:ncomp) = zi(1:ncomp)                                   ! initial estimate
        x2(1:ncomp) = rhoi_out(1:ncomp) / sum( rhoi_out(1:ncomp) )  ! initial estimate from stability analysis
        dense_start2(1) = get_eta( rhoiF )
        dense_start2(2) = get_eta( rhoi_out )
        phase_frac2 = 0.001_dp
        call Gibbs_flash ( t, p, zi, ln_zi_phi, dense_start2, phase_frac2, x1, x2, rhoiA, rhoiB, converg_GF )

        if ( converg_GF ) then
           rhoi1(:) = rhoiA(:)
           rhoi2(:) = rhoiB(:)
           x1(:) = rhoiA(:) / sum( rhoiA(1:ncomp) )
           x2(:) = rhoiB(:) / sum( rhoiB(1:ncomp) )
           converg = .true.
        end if

     end if

  end if

end subroutine pt_flash


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine Gibbs_flash
!> \brief Gibbs_flash
!!
!! algorithm described in book of Michelsen & Mollerup. Start with successive
!! substitution. Only if convergence is slow, switch to Newton scheme with full
!! Hessian.
!!
!! Input variables:  t, p, zi(:) (feed composition)
!!                   ln_zi_phi(:) (must be calculated prior to this subroutine)
!!                   dense_start(2) (starting values density (eta) of both phases)
!!                   x2(:) (starting value for composition of phase 2)
!!                   phase_frac2 (starting value mole-based fraction of phase 2)
!!
!! Output variables: converg
!!                   x1(:), x2(:), rhoi1(:), rhoi2(:)
!!                   (coexisting properties if converg=.true.)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine Gibbs_flash ( t, p, zi, ln_zi_phi, dense_start, phase_frac2, x1, x2, rhoi1, rhoi2, converg )

  use optimizer_JG, only: Newton_minimizer
  use stability, only: get_eta
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: zi
  real(dp), dimension(ncomp), intent(in)          :: ln_zi_phi
  real(dp), dimension(2), intent(in)              :: dense_start
  real(dp), intent(inout)                         :: phase_frac2
  real(dp), dimension(ncomp), intent(inout)       :: x1
  real(dp), dimension(ncomp), intent(inout)       :: x2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer                                         :: N_max_iter
  real(dp), dimension(ncomp)                      :: vi
  real(dp)                                        :: f_tpd
  real(dp), dimension(ncomp)                      :: grad
  real(dp)                                        :: grad_tolerance
  real(dp)                                        :: para_tolerance
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! define quantities needed within subroutines called by a Newton optimizer
  ! composition vector zi_flash is also needed for Gibbs_SS_minimizer (needed in
  ! SR G_flash_parameter_bound)
  !-----------------------------------------------------------------------------
  allocate( zi_flash( ncomp ) )
  allocate( ln_zi_phi_flash( ncomp ) )

  zi_flash( 1:ncomp ) = zi( 1:ncomp )
  t_flash = t
  p_flash = p
  ln_zi_phi_flash(1:ncomp) = ln_zi_phi(1:ncomp)
  dense_start_flash(1:2) = dense_start(1:2)
  
  !-----------------------------------------------------------------------------
  ! define degrees of freedom of flash optimization
  !-----------------------------------------------------------------------------
  vi(1:ncomp) = x2(1:ncomp) * phase_frac2

  !-----------------------------------------------------------------------------
  ! Successive substitution
  !-----------------------------------------------------------------------------
  grad_tolerance = 1.E-7_dp
  N_max_iter = 20
  zi_flash(1:ncomp) = zi(1:ncomp)   ! composition vector zi needed in the SR G_flash_parameter_bound

  call Gibbs_SS_minimizer ( t, p, zi, ln_zi_phi, dense_start, vi,  &
                            grad_tolerance, N_max_iter, f_tpd, grad, rhoi1, rhoi2, converg )

  !-----------------------------------------------------------------------------
  ! if unsuccessful, perform Newton iteration with fully determined Hessian
  !-----------------------------------------------------------------------------
  if ( .NOT.converg ) then

     !--------------------------------------------------------------------------
     ! T, p and composition vector zi, etc. are needed in the SR G_flash_grads
     !--------------------------------------------------------------------------
     if ( outp >= 1 ) write (*,*) 'Gibbs_SS_min not converged'
     ! ln_zi_phi_flash(1:ncomp) = ln_zi_phi(1:ncomp)   ! is this line needed, I don't think so (see above)
     ! dense_start_flash(1:2) = dense_start(1:2)       ! is this line needed, I don't think so (see above)

     N_max_iter = 100
     para_tolerance = 1.E-12

     !--------------------------------------------------------------------------
     ! Newton minimization
     !--------------------------------------------------------------------------
     call Newton_minimizer ( G_tpd_objective, G_flash_grads, G_flash_parameter_bound,  &
        ncomp, grad_tolerance, para_tolerance, N_max_iter, vi, f_tpd, grad, converg )

  end if

  !-----------------------------------------------------------------------------
  ! if flash converged, convert parameters back to compositions of both phases
  !-----------------------------------------------------------------------------
  if ( converg ) then

     x1(1:ncomp) = zi(1:ncomp) - vi(1:ncomp)
     x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )

     phase_frac2 = sum( vi(1:ncomp) )
     x2(1:ncomp) = vi(1:ncomp) / phase_frac2

     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'Gibbs flash converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  x p1   = ',x1( 1:ncomp )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',x2( 1:ncomp )
     if ( outp >= 2 ) write (*,*) 'eta 1,2 = ',get_eta( rhoi1 ), get_eta( rhoi2 )
     if ( outp >= 2 ) write (*,*) ' '

  end if

  deallocate( zi_flash, ln_zi_phi_flash )

end subroutine Gibbs_flash


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine Gibbs_SS_minimizer
!> \brief Gibbs_SS_minimizer
!!
!! Gibbs flash using successive substitution.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine Gibbs_SS_minimizer ( t, p, zi, ln_zi_phi, dense_start, vi, grad_tolerance, N_max_iter,  &
                                f, g, rhoi1, rhoi2, converged )

  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: zi
  real(dp), dimension(ncomp), intent(in)          :: ln_zi_phi
  real(dp), dimension(2), intent(in)              :: dense_start
  real(dp), dimension(ncomp), intent(in out)      :: vi
  real(dp), intent(in)                            :: grad_tolerance
  integer, intent(in)                             :: N_max_iter
  real(dp)                                        :: f
  real(dp), dimension(ncomp), intent(out)         :: g
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  logical, intent(out)                            :: converged

  !-----------------------------------------------------------------------------
  integer                                         :: count
  real(dp), dimension(ncomp)                      :: x1
  real(dp), dimension(ncomp)                      :: x2
  real(dp), dimension(ncomp)                      :: si
  real(dp), dimension(ncomp)                      :: vi_old
  real(dp), dimension(ncomp)                      :: delta_vi
  real(dp), dimension(ncomp)                      :: gv
  real(dp), dimension(ncomp)                      :: gl
  real(dp), dimension(ncomp)                      :: lnphi
  real(dp)                                        :: tpd1, tpd2
  real(dp)                                        :: f_old
  real(dp)                                        :: phase_frac2
  real(dp)                                        :: damp
  real(dp)                                        :: term
  real(dp)                                        :: error_in_g
  logical                                         :: bound_violation
  logical                                         :: adjust_damping
  !-----------------------------------------------------------------------------

  if ( outp >= 2 ) write (*,*) ' '
  if ( outp >= 2 ) write (*,*) 'entering Gibbs_SS_minimizer'

  !-----------------------------------------------------------------------------
  ! initialize iteration
  !-----------------------------------------------------------------------------
  count = 0
  converged = .false.

  damp = 1.0_dp

  vi_old(1:ncomp) = vi(1:ncomp)
  f_old = 1.E50_dp

  do while ( .NOT.converged .AND. count < N_max_iter )

     count = count + 1

     !--------------------------------------------------------------------------
     ! chemical potentials phase 2
     !--------------------------------------------------------------------------
     phase_frac2 = sum( vi(1:ncomp) )
     x2(1:ncomp) = vi(1:ncomp) / phase_frac2
     call chemical_potential_tp ( t, p, x2, dense_start(2), rhoi2(:), lnphi(:) )

     gv(1:ncomp) = log( x2(1:ncomp) ) + lnphi(1:ncomp)
     tpd2 = sum( x2(1:ncomp) * ( gv(1:ncomp) - ln_zi_phi(1:ncomp) ) )

     !--------------------------------------------------------------------------
     ! chemical potentials phase 1
     !--------------------------------------------------------------------------
     x1(1:ncomp) = zi(1:ncomp) - vi(1:ncomp)
     x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
     call chemical_potential_tp ( t, p, x1, dense_start(1), rhoi1(:), lnphi(:) )

     gl(1:ncomp) = log( x1(1:ncomp) ) + lnphi(1:ncomp)
     tpd1 = sum( x1(1:ncomp) * ( gl(1:ncomp) - ln_zi_phi(1:ncomp) ) )

     !--------------------------------------------------------------------------
     ! gradient of tangent plane fct ( isofugacity condition, iterated to zero )
     !--------------------------------------------------------------------------
     g(1:ncomp) = gv(1:ncomp) - gl(1:ncomp) 
     if ( outp >= 3 ) write (*,*) 'error in g:', g(1:ncomp)

     !--------------------------------------------------------------------------
     ! calculate value of Gibbs tangent plane function f
     !--------------------------------------------------------------------------
     f = phase_frac2 * tpd2 + ( 1.0_dp - phase_frac2 ) * tpd1

     !--------------------------------------------------------------------------
     ! determine si
     !--------------------------------------------------------------------------
     si(1:ncomp) = x1(1:ncomp) * x2(1:ncomp) / zi(1:ncomp)

     !--------------------------------------------------------------------------
     ! vector delta_vi(:) of successive substitution step
     !--------------------------------------------------------------------------
     term = sum( si(1:ncomp) * g(1:ncomp) ) / ( sum( si(1:ncomp) ) - 1.0_dp )
     delta_vi(:) = phase_frac2 * (1.0_dp-phase_frac2) * si(1:ncomp) * ( term - g(1:ncomp) )

     adjust_damping = .true.

     do while ( adjust_damping )

        !-----------------------------------------------------------------------
        ! are bounds of parameters violated?
        !-----------------------------------------------------------------------
        call G_flash_parameter_bound  ( ncomp, damp*delta_vi, vi_old, bound_violation )

        if ( bound_violation ) then
           if ( damp > 0.2_dp ) then
              damp = damp * 0.5_dp
              if ( outp >= 2 ) write (*,*) 'damping the step-size'
           else
              if ( outp > 0 ) write (*,*) 'Gibbs_SS_minimizer: do something here:'
              if ( outp > 0 ) write (*,*) 'remove bound_violation and proceed for one or two steps'
              exit
           end if
        else
           adjust_damping = .false.
        end if

     end do

     if ( bound_violation ) exit

     !--------------------------------------------------------------------------
     ! new vi vector
     !--------------------------------------------------------------------------
     vi(1:ncomp) = vi_old(1:ncomp) + damp * delta_vi(1:ncomp)

        
     !--------------------------------------------------------------------------
     ! if f is successfully decreased
     !--------------------------------------------------------------------------
     if ( f < f_old + 1.E-3_dp ) then
        vi_old(1:ncomp) = vi(1:ncomp)
        f_old = f
        damp = min( damp * 2.0_dp, 1.0_dp )
     else
        if ( outp >= 1 ) write (*,*) 'Gibbs_SS_minimizer: no decrease of obj. fct possible', f, f_old
        if ( outp >= 1 ) write (*,*) 'damping needed?????'
        vi(1:ncomp) = vi_old(1:ncomp)
        f = f_old
     end if

     !--------------------------------------------------------------------------
     ! determine error
     !--------------------------------------------------------------------------

     error_in_g = sum( abs( g(1:ncomp) ) )

     if ( error_in_g < grad_tolerance ) then
        converged = .true.
     end if

     if ( outp >= 2 ) write (*,'(a,i3,2G25.15)') 'Gibbs_SS_minimizer: finished SS-step: N, f, errors ',count, f, error_in_g

  end do

end subroutine Gibbs_SS_Minimizer


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE G_tpd_objective
!> \brief G_tpd_objective
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine G_tpd_objective ( n, vi, f )

  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real(dp), dimension(n), intent(in)     :: vi
  real(dp), intent(out)                  :: f

  !-----------------------------------------------------------------------------
  real(dp)                               :: t
  real(dp)                               :: p
  real(dp), dimension(ncomp)             :: x1
  real(dp), dimension(ncomp)             :: x2
  real(dp), dimension(ncomp)             :: rhoi1, rhoi2
  real(dp), dimension(ncomp)             :: lnphi_1, lnphi_2
  real(dp), dimension(2)                 :: dense_start
  real(dp)                               :: tpd1, tpd2
  real(dp)                               :: phase_frac2
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! temperature and pressure and feed composition
  !-----------------------------------------------------------------------------
  t = t_flash
  p = p_flash
  dense_start(1:2) = dense_start_flash(1:2)

  !-----------------------------------------------------------------------------
  ! rewrite mol numbers of phase 2 (vi) to mole fraction of both phases
  !-----------------------------------------------------------------------------
  phase_frac2 = sum( vi(1:ncomp) )
  x2(1:ncomp) = vi(1:ncomp) / phase_frac2

  x1(1:ncomp) = zi_flash(1:ncomp) - vi(1:ncomp)
  x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! Helmholtz energy densities (of first and second trial phase)
  !-----------------------------------------------------------------------------
  call chemical_potential_tp ( t, p, x2, dense_start(2), rhoi2(:), lnphi_2(:) )
  tpd2 = sum( x2(1:ncomp) * ( log( x2(1:ncomp) ) + lnphi_2(1:ncomp) - ln_zi_phi_flash(1:ncomp) ) )

  call chemical_potential_tp ( t, p, x1, dense_start(1), rhoi1(:), lnphi_1(:) )
  tpd1 = sum( x1(1:ncomp) * ( log( x1(1:ncomp) ) + lnphi_1(1:ncomp) - ln_zi_phi_flash(1:ncomp) ) )

  !-----------------------------------------------------------------------------
  ! calculate value of Gibbs tangent plane function f
  !-----------------------------------------------------------------------------
  f = phase_frac2 * tpd2 + ( 1.0_dp - phase_frac2 ) * tpd1

end subroutine G_tpd_objective



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE G_flash_grads
!> \brief G_flash_grads
!!
!! the diagonal matrix is not entirely diagonal. Is that an advantage? Other-
!! wise, that should be changed.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine G_flash_grads ( n, vi, grad, hessian, diagonal )

  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real(dp), dimension(n), intent(in)     :: vi
  real(dp), dimension(n), intent(out)    :: grad
  real(dp), dimension(n,n), intent(out)  :: hessian
  real(dp), dimension(n,n), intent(out)  :: diagonal

  !-----------------------------------------------------------------------------
  integer                                :: i, j
  real(dp)                               :: t
  real(dp)                               :: p
  real(dp)                               :: phase_frac2
  real(dp), dimension(2)                 :: dense_start
  real(dp), dimension(ncomp)             :: x1
  real(dp), dimension(ncomp)             :: x2
  real(dp), dimension(ncomp)             :: rhoi_out
  real(dp), dimension(ncomp)             :: lnphi_1, lnphi_2
  real(dp), dimension(ncomp,ncomp)       :: lnphi_1_nj, lnphi_2_nj
  !-----------------------------------------------------------------------------

  if ( n /= ncomp ) then
     write (*,*) 'error in G_flash_grads'
     stop
  end if
  
  !-----------------------------------------------------------------------------
  ! temperature and pressure and feed composition
  !-----------------------------------------------------------------------------
  t = t_flash
  p = p_flash
  dense_start(1:2) = dense_start_flash(1:2)

  !-----------------------------------------------------------------------------
  ! rewrite mol numbers of phase 2 (vi) to mole fraction of both phases
  !-----------------------------------------------------------------------------
  phase_frac2 = sum( vi(1:ncomp) )
  x2(1:ncomp) = vi(1:ncomp) / phase_frac2

  x1(1:ncomp) = zi_flash(1:ncomp) - vi(1:ncomp)
  x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! first derivative and second derivatives of Helmholtz energy densities
  ! (i.e. chemical potentials and Hessian matrix) of first and second trial phase.
  ! All quantities are calculated at constant T, p.
  !-----------------------------------------------------------------------------
  call chemical_potential_tp ( t, p, x2, dense_start(2), rhoi_out, lnphi_2, lnphi_nj=lnphi_2_nj )
  !tpd2 = sum( x2(1:ncomp) * ( log( x2(1:ncomp) ) + lnphi_2(1:ncomp) - ln_zi_phi(1:ncomp) ) )

  call chemical_potential_tp ( t, p, x1, dense_start(1), rhoi_out, lnphi_1, lnphi_nj=lnphi_1_nj )
  !tpd1 = sum( x1(1:ncomp) * ( log( x1(1:ncomp) ) + lnphi_1(1:ncomp) - ln_zi_phi(1:ncomp) ) )

  !-----------------------------------------------------------------------------
  ! calculate value of Gibbs tangent plane function f (objective function)
  !-----------------------------------------------------------------------------
  !f = phase_frac2 * tpd2 + ( 1.0_dp - phase_frac2 ) * tpd1

  !-----------------------------------------------------------------------------
  ! derivatives of objective function
  !-----------------------------------------------------------------------------

  grad(1:ncomp) = log( x2(1:ncomp) ) + lnphi_2(1:ncomp) - log( x1(1:ncomp) ) - lnphi_1(1:ncomp)

  !-----------------------------------------------------------------------------
  ! second derivatives of Helmholtz energy densities (Hessian matrix)
  ! of first and second trial phase
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     do j = 1, ncomp
        hessian(i,j) = ( 1.0_dp - phase_frac2 ) * lnphi_2_nj(i,j) + phase_frac2 * lnphi_1_nj(i,j) - 1.0_dp
     end do
  end do
  hessian(:,:) = hessian(:,:) / ( phase_frac2 * ( 1.0_dp - phase_frac2 ) )

  diagonal(:,:) = 0.0_dp
  do i = 1, ncomp
     diagonal(i,i) = diagonal(i,i) + zi_flash(i) / x1(i) / x2(i)
  end do
  diagonal(:,:) = diagonal(:,:) / ( phase_frac2 * ( 1.0_dp - phase_frac2 ) )

  hessian(:,:) = hessian(:,:) + diagonal(:,:)

end subroutine G_flash_grads



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE G_flash_parameter_bound
!> \brief G_flash_parameter_bound
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
subroutine G_flash_parameter_bound ( n, delta_vi, vi_old, bound_violation )

  use stability, only: get_eta
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real(dp), dimension(n), intent(in)     :: delta_vi
  real(dp), dimension(n), intent(in)     :: vi_old
  logical, intent(out)                   :: bound_violation

  !-----------------------------------------------------------------------------
  integer                                :: i
  real(dp), dimension(ncomp)             :: vi
  real(dp), dimension(ncomp)             :: li
  real(dp), dimension(ncomp)             :: li_old
  real(dp)                               :: max_step_1, max_step_2
  !-----------------------------------------------------------------------------

  bound_violation = .false.

  vi(:) = vi_old(:) + delta_vi(:)
  li(:) = zi_flash(:) - vi(:)
  ! write (*,*) 'vi=',vi
  ! write (*,*) 'li=',li
  ! write (*,*) 'zi=',zi_flash

  !-----------------------------------------------------------------------------
  ! assess parameter-bounds
  !-----------------------------------------------------------------------------
  do i = 1, ncomp
     ! if ( vi(i) < 0.0_dp ) write (*,*) ' xi(i) phase 2 below zero, i=',i
     ! if ( li(i) < 0.0_dp ) write (*,*) ' xi(i) phase 1 below zero, i=',i
     if ( vi(i) < 0.0_dp ) bound_violation = .true.
     if ( li(i) < 0.0_dp ) bound_violation = .true.
  end do

  !-----------------------------------------------------------------------------
  ! assess maximum step size, start by calculating old component density values
  !-----------------------------------------------------------------------------

  li_old(:) = zi_flash(:) - vi_old(:)
  
  max_step_1 = maxval( abs( vi(:) / sum( vi(:) ) - vi_old(:) / sum( vi_old(:) ) ) )
  max_step_2 = maxval( abs( li(:) / sum( li(:) ) - li_old(:) / sum( li_old(:) ) ) )

  if ( max_step_1 > 0.2_dp .OR. max_step_2 > 0.2_dp ) bound_violation = .true.

end subroutine G_flash_parameter_bound



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine rachford_rice
!> \brief rachford_rice
!!
!! This subroutine performs a p-T flash for a defined composition (xiF) in a
!! feed-trial phase by iterating the Rachford-Rice equation.
!!
!! input:   t, p, xiF(:)
!!          optional input: x1s(:) and x2s(:) as starting values for mole fractions
!!
!! output:  rhoi1(:), rhoi2(:)
!!          xi(1,:), xi(2,:)
!!          values of: dense(:), lnphi(:,:) are also converged values, if converg=.true.
!!          phase_frac2 (mole-based fraction of phase 2)
!!          slow_iteration is .true. if convergence not reached but solution seems possible
!!          converg is .true. if a solution is found
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine rachford_rice ( t, p, xiF, rhoi1, rhoi2, phase_frac2, slow_iteration, converg, x1s, x2s )

  use PARAMETERS, only: machine_eps
  use stability, only: get_eta
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(out)         :: rhoi1
  real(dp), dimension(ncomp), intent(out)         :: rhoi2
  real(dp), intent(out)                           :: phase_frac2
  logical, intent(out)                            :: converg
  logical, intent(out)                            :: slow_iteration
  real(dp), dimension(ncomp),intent(in), optional :: x1s
  real(dp), dimension(ncomp),intent(in), optional :: x2s

  !-----------------------------------------------------------------------------
  integer                                         :: k1, k2
  integer                                         :: nr_inside_steps, nr_outside_steps
  real(dp), parameter                             :: tol_out = 1.E-10_dp
  real(dp), parameter                             :: tol_in = 1.E-13_dp
  real(dp)                                        :: phase_frac2_0, phase_frac2_00
  real(dp), dimension(ncomp)                      :: lnphi_1, lnphi_2
  real(dp), dimension(ncomp)                      :: x1, x2
  real(dp), dimension(ncomp)                      :: x1_old, x2_old
  real(dp), dimension(ncomp)                      :: Ki
  real(dp)                                        :: dense_start1, dense_start2
  real(dp)                                        :: f1, d_ph, dfdph
  real(dp)                                        :: error1, error2
  real(dp)                                        :: phase_frac2_min, phase_frac2_max
  real(dp)                                        :: g_0, g_1
  logical                                         :: RR_vle

  integer, parameter                              :: N_acc_cycle = 8
  real(dp), dimension(ncomp)                      :: lnK_old1, lnK_old2
  real(dp), dimension(ncomp)                      :: del1, del2
  real(dp)                                        :: lambda
  real(dp), dimension(2)                          :: error1_hist
  real(dp), parameter                             :: damping = 0.5_dp
  !-----------------------------------------------------------------------------

  converg = .false.
  slow_iteration = .false.

  nr_outside_steps = 80
  nr_inside_steps = 10

  error1_hist(:) = 1.0_dp

!!$  phase_frac2 = 0.5_dp
!!$  phase_frac2_0 = phase_frac2
!!$  d_ph = 0.00001_dp
!!$  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'RR, N_out, N_inner, x1(1), xiF(1), x2(1), phase_frac2, error1'

  dense_start1 = 0.4_dp       ! Index 1 is for liquid density (here: packing fraction eta)
  dense_start2 = 1.E-5_dp     ! Index 2 is for vapour density (here: packing fraction eta)

  !-----------------------------------------------------------------------------
  ! starting values for Ki-values:
  ! liquid with x1=xiF(:); vapor x2 assuming ideal gas (phi_2(i) = 1)
  !-----------------------------------------------------------------------------
  if ( present( x1s ) .AND. present( x2s ) ) then
     
     x1( 1:ncomp ) = x1s( 1:ncomp )
     x2( 1:ncomp ) = x2s( 1:ncomp )

     call chemical_potential_tp ( t, p, x1, dense_start1, rhoi1, lnphi_1 )

     call chemical_potential_tp ( t, p, x2, dense_start2, rhoi2, lnphi_2 )

     Ki(1:ncomp) = exp( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

     !--------------------------------------------------------------------------
     ! test for RR-conditions
     !--------------------------------------------------------------------------
     g_0 = sum( Ki(1:ncomp) * xiF(1:ncomp) ) - 1.0_dp
     g_1 = 1.0_dp - sum( xiF(1:ncomp) / Ki(1:ncomp) )

     if ( g_0 < 0.0_dp .OR. g_1 > 0.0_dp ) then
        call correct_vapor_or_liquid_phase ( t, p, xiF, g_0, g_1, x1, x2, lnphi_1, lnphi_2, Ki )
     end if

     phase_frac2 = 0.5_dp

  else

     call check_rachford_rice_VLE ( t, p, xiF, x1, x2, Ki, phase_frac2, RR_vle )

  end if

!!$  !x2( 1:ncomp ) = xiF( 1:ncomp )
!!$  x2( 1:ncomp ) = x1( 1:ncomp ) * exp( lnphi_1(1:ncomp) )
!!$  x2( 1:ncomp ) = x2( 1:ncomp ) / sum( x2( 1:ncomp ) )
!!$  call chemical_potential_tp ( t, p, x2, dense_start2, rhoi2, lnphi_2 )
!!$  write (*,*) 'sum( rhoi2(:) )',sum( rhoi2(:) )
!!$
!!$!  if ( abs( get_eta( rhoi2 ) / get_eta( rhoi1 ) - 1.0_dp ) < machine_eps*1.E4_dp ) then
!!$!     x2( 1:ncomp ) = x1( 1:ncomp ) * exp( lnphi_1(1:ncomp) )
!!$!     x2( 1:ncomp ) = x2( 1:ncomp ) / sum( x2( 1:ncomp ) )
!!$!     call chemical_potential_tp ( t, p, x2, dense_start2, rhoi2, lnphi_2 )
!!$!  end if
!!$
!!$  !-----------------------------------------------------------------------------
!!$  ! initial estimate for composition, based on Ki and starting value for phase_frac2
!!$  !-----------------------------------------------------------------------------
!!$
!!$  Ki(1:ncomp) = EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )
!!$  x1(1:ncomp) = xiF(1:ncomp) / ( 1.0_dp - phase_frac2 + phase_frac2 *Ki(1:ncomp) )
!!$  x2(1:ncomp) = Ki(1:ncomp)* x1(1:ncomp)
!!$  x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
!!$  x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )

  phase_frac2_0 = phase_frac2
  d_ph = 0.00001_dp
  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'RR, N_out, N_inner, x1(1), xiF(1), x2(1), phase_frac2, error1'

  !-----------------------------------------------------------------------------
  ! start iteration
  !-----------------------------------------------------------------------------

  x1_old(:) = x1(:)
  x2_old(:) = x2(:)
  lnK_old1(:) = 0.0_dp
  ! Ki(:) = 1.0_dp

  error2 = 0.0_dp

  k1 = 0
  !i_acc = 0
  !i_cycle = 0
  !lnK_old(:) = 0.0_dp
  error1 = tol_out + 1.0_dp

  do while ( error1 > tol_out .AND. k1 < nr_outside_steps )

     !--------------------------------------------------------------------------
     ! outer loop (converging the mole fractions)
     !--------------------------------------------------------------------------

     k1 = k1 + 1
     !i_cycle = i_cycle + 1

     ! if ( k1 > 1 ) then
        call chemical_potential_tp ( t, p, x1, dense_start1, rhoi1, lnphi_1 )
        call chemical_potential_tp ( t, p, x2, dense_start2, rhoi2, lnphi_2 )
     ! end if
     ! use converged value of eta as starting value for next iteration
     dense_start1 = get_eta( rhoi1 )
     dense_start2 = get_eta( rhoi2 )

     lnK_old2(1:ncomp) = lnK_old1(1:ncomp)
     lnK_old1(1:ncomp) = log( Ki(1:ncomp) )
     Ki(1:ncomp) = EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

     del2( 1:ncomp ) = lnK_old1(1:ncomp) - lnK_old2(1:ncomp)
     del1( 1:ncomp ) = lnphi_1( 1:ncomp ) - lnphi_2( 1:ncomp ) - lnK_old1(1:ncomp)

     !--------------------------------------------------------------------------
     ! accelerate convergence
     !--------------------------------------------------------------------------
     if ( mod( k1, N_acc_cycle ) == 0 .AND. k1 > 3 * N_acc_cycle .AND.  &
                      phase_frac2 > 1.E-8_dp .AND. phase_frac2 < 1.0_dp-1.E-8_dp ) then
        lambda = sum( xiF(1:ncomp) * del1( 1:ncomp )**2 )
        lambda = lambda / sum( xiF(1:ncomp) * del1( 1:ncomp ) * del2( 1:ncomp ) )
        !lambda = sum( del1( 1:ncomp )**2 )
        !lambda = lambda / sum( del1(1:ncomp) * del2(1:ncomp) )
        lambda = min( lambda, 0.85_dp )
        Ki(1:ncomp) = exp( log( Ki(1:ncomp) ) + del1(1:ncomp) * lambda/(1.0_dp-lambda) )
        error1_hist(:) = 1.0_dp
        if ( outp >= 2 ) write (*,*) 'accelerating RR',lambda
     end if

     !--------------------------------------------------------------------------
     ! check for trivial solution
     !--------------------------------------------------------------------------
     if ( ABS(1.0_dp-sum(rhoi1(:))/sum(rhoi2(:))) < 1.E-5_dp  &
          .AND. SUM(ABS(Ki(1:ncomp)-1.0_dp)) < 1.E-6_dp ) then

        call chemical_potential_tp ( t, p, x1, 0.4_dp, rhoi1, lnphi_1 )
        call chemical_potential_tp ( t, p, x2, 1.E-5_dp, rhoi2, lnphi_2 )
        dense_start1 = get_eta( rhoi1 )
        dense_start2 = get_eta( rhoi2 )
        Ki(1:ncomp) = EXP( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )
        if ( ABS(1.0_dp-sum(rhoi1(:))/sum(rhoi2(:))) < 1.E-5_dp  &
             .AND. SUM(ABS(Ki(1:ncomp)-1.0_dp)) < 1.E-6_dp ) then
           if ( outp >= 1 ) write (*,*) 'RR running into trivial solution, exit.'
           exit
        end if

     end if

     !--------------------------------------------------------------------------
     ! can a solution be expected with these iteration values?
     !--------------------------------------------------------------------------
     g_0 = sum( xiF(:) * Ki(:) ) - 1.0_dp
     g_1 = 1.0_dp - sum( xiF(:) / Ki(:) )
     !if ( g_0  > 0.0_dp .OR. g_1 < 0.0_dp ) then
     !   call correct_vapor_or_liquid_phase ( t, p, xiF, g_0, g_1, x1, x2, lnphi_1, lnphi_2, Ki )
     !end if

     phase_frac2_min = max( 0.0_dp, maxval( ( Ki(:)*xiF(:) - 1.0_dp ) / ( Ki(:) - 1.0_dp ), MASK = Ki > 1.0_dp ) )
     phase_frac2_max = min( 1.0_dp, minval( ( xiF(:) - 1.0_dp ) / ( Ki(:) - 1.0_dp ), MASK = Ki < 1.0_dp ) )
     if ( k1 <= 2 ) phase_frac2 = 0.5_dp * ( phase_frac2_min + phase_frac2_max )

     k2 = 0
     phase_frac2_00 = phase_frac2
     error2 = tol_in + 1.0_dp
     do while ( error2 > tol_in .AND. k2 < nr_inside_steps )

        !-----------------------------------------------------------------------
        ! inner loop (converging the fraction of one phase)
        !-----------------------------------------------------------------------

        k2 = k2 + 1

        phase_frac2_0 = phase_frac2
        f1 = SUM( xiF(1:ncomp)*(Ki(1:ncomp)-1.0_dp) / ( 1.0_dp - phase_frac2 + phase_frac2*Ki(1:ncomp) ) )

        dfdph = - SUM( xiF(1:ncomp) / ( 1.0_dp/(Ki(1:ncomp)-1.0_dp) + phase_frac2 )**2 )
        if ( abs( dfdph ) < machine_eps ) dfdph = 0.0000001_dp

        if ( f1 < 0.0_dp ) phase_frac2_max = phase_frac2
        if ( f1 > 0.0_dp ) phase_frac2_min = phase_frac2

        phase_frac2 = phase_frac2_0 - f1 / dfdph
        
        if ( phase_frac2 < phase_frac2_min .OR. phase_frac2 > phase_frac2_max ) then
           phase_frac2 = 0.5_dp * ( phase_frac2_min + phase_frac2_max )
        end if

        error2 = ABS( phase_frac2 - phase_frac2_0 ) + abs( f1 / dfdph )
        if ( outp >= 3 ) write (*,'(a,5G25.15)') 'phase_fraction, error', phase_frac2, f1/ dfdph, error2

     end do

     !--------------------------------------------------------------------------
     ! for slow convergence and likely numerical noise, apply damping
     !--------------------------------------------------------------------------
     if ( sum( error1_hist(1:2) ) < tol_out * 100.0_dp ) then
        phase_frac2 = damping * phase_frac2 + ( 1.0_dp - damping ) * phase_frac2_0
     end if

     !--------------------------------------------------------------------------
     ! determine x and error of outer loop
     !--------------------------------------------------------------------------
     x1(1:ncomp) = xiF(1:ncomp) / ( 1.0_dp - phase_frac2 + phase_frac2 *Ki(1:ncomp) )
     x2(1:ncomp) = Ki(1:ncomp)* x1(1:ncomp)

     ! if ( sum( x1(1:ncomp) ) < 0.2_dp ) x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
     ! if ( sum( x2(1:ncomp) ) < 0.2_dp ) x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )
     x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
     x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )

     error1 =  SUM( ABS( x1_old(1:ncomp) - x1(1:ncomp) ) ) + SUM( ABS( x2_old(1:ncomp) - x2(1:ncomp) ) )

     if ( sum( error1_hist(1:2) ) < tol_out * 80.0_dp .AND. k1 > 15 ) then
        if ( outp > 2 )  write (*,*) 'now relaxing tolerance', error1, error1 * 0.05_dp
        error1 = error1 * 0.05_dp
     end if

     x1_old(:) = x1(:)
     x2_old(:) = x2(:)

     error1_hist(2) = error1_hist(1)
     error1_hist(1) = error1
 
     if ( outp >= 2 ) write (*,'(a,2i5,5G20.12)') 'RR', k1, k2, x1(1), xiF(1), x2(1), phase_frac2, error1

     if ( k1 == 100 .AND. outp > 0 ) write (*,*) 'Rachford-Rice: outside loop slow convergence'

  end do


  !-----------------------------------------------------------------------------
  ! if convergence: accept solution
  !-----------------------------------------------------------------------------

  if ( error1 < tol_out .AND. error2 < tol_in ) then

     converg = .true.
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'Rachford-Rice converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) 'eta 1,2 = ',get_eta( rhoi1 ), get_eta( rhoi2 )
     if ( outp >= 2 ) write (*,*) ' '

  else

     if ( outp >= 1 ) write (*,*) 'Rachford-Rice NOT converged'

  end if

  if ( .NOT.converg .AND. k1 == nr_outside_steps ) then

     slow_iteration = .true.
     if ( outp >= 2 ) write (*,*) 'slow_iteration = .true.'

  end if

end subroutine rachford_rice


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine determine_ln_zi_phi
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine determine_ln_zi_phi ( t, p, xiF, ln_zi_phi, etaF )

  use stability, only: get_eta
  use properties, only: chemical_potential_tp, state_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(out)         :: ln_zi_phi
  real(dp), intent(out)                           :: etaF

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(ncomp)                      :: lnphi_1, lnphi_2
  real(dp), dimension(ncomp)                      :: mu_1, mu_2
  real(dp), dimension(ncomp)                      :: lnphi_F
  real(dp)                                        :: dense_start
  real(dp)                                        :: gL, gV
  !-----------------------------------------------------------------------------

  dense_start = 0.45_dp                        ! starting density for a liquid phase
  call chemical_potential_tp ( t, p, xiF, dense_start, rhoi1, lnphi_1, mu=mu_1 )
  gL = sum( mu_1(:) )

  dense_start = 1.E-6_dp                       ! starting density for a vapor phase
  call chemical_potential_tp ( t, p, xiF, dense_start, rhoi2, lnphi_2, mu=mu_2 )
  gV = sum( mu_2(:) )

  !-----------------------------------------------------------------------------
  ! determine condition with lower Gibbs-energy (condition with higher value is instable)
  !-----------------------------------------------------------------------------
  if ( gL < gV ) then
     lnphi_F(:) = lnphi_1(:)
     etaF = get_eta( rhoi1(:) )
  else
     lnphi_F(:) = lnphi_2(:)
     etaF = get_eta( rhoi2(:) )
  end if

  do i = 1, ncomp
     if ( xiF(i) > 1.E-100_dp ) then
        ln_zi_phi(i) = log( xiF(i) ) + lnphi_F(i)
     else
        ln_zi_phi(i) = log( 1.E-100_dp ) + lnphi_F(i)
     end if
  end do

end subroutine determine_ln_zi_phi


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine check_rachford_rice_VLE
!> \brief check_rachford_rice_VLE
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_rachford_rice_VLE ( t, p, xiF, x1, x2, Ki, phase_frac2, RR_vle )

  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(out)         :: x1
  real(dp), dimension(ncomp), intent(out)         :: x2
  real(dp), dimension(ncomp), intent(out)         :: Ki
  real(dp), intent(out)                           :: phase_frac2
  logical, intent(out)                            :: RR_vle

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(ncomp)                      :: lnphi_1, lnphi_2
  real(dp)                                        :: g_0, g_1
  real(dp)                                        :: dense_start
  !-----------------------------------------------------------------------------

  dense_start = 0.4_dp
  call chemical_potential_tp ( t, p, xiF, dense_start, rhoi1, lnphi_1 )

  dense_start = 1.E-6_dp
  call chemical_potential_tp ( t, p, xiF, dense_start, rhoi2, lnphi_2 )

  !-----------------------------------------------------------------------------
  ! if a vapor and a liquid phase can not be calculated for the feed composition
  !-----------------------------------------------------------------------------
  if ( abs( sum( rhoi1(1:ncomp) ) / sum( rhoi2(1:ncomp) ) - 1.0_dp ) < 1.E-5_dp ) then

     !--------------------------------------------------------------------------
     ! assess vapor with composition from ideal gas assumption (phi=1 for vapor)
     !--------------------------------------------------------------------------
     x2( 1:ncomp ) = xiF( 1:ncomp ) * exp( lnphi_1( 1:ncomp ) )
     x2( 1:ncomp ) = x2( 1:ncomp ) / sum( x2( 1:ncomp ) )
     dense_start = 1.E-6_dp
     call chemical_potential_tp ( t, p, x2, dense_start, rhoi2, lnphi_2 )

  else

     if ( outp >= 2 ) write (*,*) ' feed can exist in both, vapor or liquid density, at given t,p'

  end if

  Ki(1:ncomp) = exp( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

  !--------------------------------------------------------------------------
  ! test for RR-conditions
  !--------------------------------------------------------------------------
  g_0 = sum( Ki(1:ncomp) * xiF(1:ncomp) ) - 1.0_dp
  g_1 = 1.0_dp - sum( xiF(1:ncomp) / Ki(1:ncomp) )

  if ( g_0 < 0.0_dp .OR. g_1 > 0.0_dp ) then
     call correct_vapor_or_liquid_phase ( t, p, xiF, g_0, g_1, x1, x2, lnphi_1, lnphi_2, Ki )
  end if

  if ( g_0 > 0.0_dp .AND. g_1 < 0.0_dp ) then
     RR_vle = .true.
  else
     RR_vle = .false.
  end if

  !-----------------------------------------------------------------------------
  ! initial estimate for composition, based on Ki and starting value for phase_frac2
  !-----------------------------------------------------------------------------

  phase_frac2 = 0.5_dp

  x1(1:ncomp) = xiF(1:ncomp) / ( 1.0_dp - phase_frac2 + phase_frac2 *Ki(1:ncomp) )
  x2(1:ncomp) = Ki(1:ncomp)* x1(1:ncomp)
  x1(1:ncomp) = x1(1:ncomp) / sum( x1(1:ncomp) )
  x2(1:ncomp) = x2(1:ncomp) / sum( x2(1:ncomp) )

end subroutine check_rachford_rice_VLE


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine correct_vapor_or_liquid_phase
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine correct_vapor_or_liquid_phase ( t, p, xiF, g_0, g_1, x1, x2, lnphi_1, lnphi_2, Ki )

  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), intent(inout)                         :: g_0
  real(dp), intent(inout)                         :: g_1
  real(dp), dimension(ncomp), intent(inout)       :: x1
  real(dp), dimension(ncomp), intent(inout)       :: x2
  real(dp), dimension(ncomp), intent(inout)       :: lnphi_1
  real(dp), dimension(ncomp), intent(inout)       :: lnphi_2
  real(dp), dimension(ncomp), intent(out)         :: Ki

  !-----------------------------------------------------------------------------
  integer                                         :: count
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp)                                        :: dense_start
  !-----------------------------------------------------------------------------

  count = 0

  do while ( ( g_0 < 0.0_dp .OR. g_1 > 0.0_dp ) .AND. count < 20 )

     count = count + 1

     if ( g_0 < 0.0_dp) then
        x1(:) = xiF(:)
        x2(:) = Ki(:) * xiF(:) / sum( Ki(:) * xiF(:) )
        ! phase_frac2 = 1.0E-5_dp
        dense_start = 1.E-6_dp
        call chemical_potential_tp ( t, p, x2, dense_start, rhoi2, lnphi_2 )
        if ( outp > 1 ) write (*,*) 'case: g(0)<0,  x2=', x2(:)
     end if

     if ( g_1 > 0.0_dp) then
        x2(:) = xiF(:)
        x1(:) =  xiF(:) / Ki(:) / sum( xiF(:) / Ki(:) )
        ! phase_frac2 = 1.0_dp - 1.0E-5_dp
        dense_start = 0.4_dp
        call chemical_potential_tp ( t, p, x1, dense_start, rhoi1, lnphi_1 )
        if ( outp > 1 ) write (*,*) 'case: g(1)>0,  x1=', x1(:)
     end if

     Ki(1:ncomp) = exp( lnphi_1(1:ncomp) - lnphi_2(1:ncomp) )

     g_0 = sum( Ki(1:ncomp) * xiF(1:ncomp) ) - 1.0_dp
     g_1 = 1.0_dp - sum( xiF(1:ncomp) / Ki(1:ncomp) )

     if ( g_0 < 0.0_dp .AND. sum( abs( x2(:) - Ki(:) * xiF(:) / sum( Ki(:) * xiF(:) ) ) ) < 1.E-7_dp ) exit
     if ( g_1 > 0.0_dp .AND. sum( abs( x1(:) -  xiF(:) / Ki(:) / sum( xiF(:) / Ki(:) ) ) ) < 1.E-7_dp ) exit

  end do

end subroutine correct_vapor_or_liquid_phase


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine three_phase_flash
!> \brief three_phase_flash
!!
!! algorithm described in book of Michelsen & Mollerup. Start with successive
!! substitution. Only if convergence is slow, switch to Newton scheme with full
!! Hessian.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


subroutine three_phase_flash ( nphase, t, p, xiF, zi, dense_start, xi, rho, converg )

  use stability, only: get_eta
  use properties, only: chemical_potential_tp
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: nphase
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), dimension(ncomp), intent(in)          :: zi
  real(dp), dimension(nphase), intent(in)         :: dense_start
  real(dp), dimension(ncomp,nphase), intent(inout) :: xi
  real(dp), dimension(nphase), intent(out)        :: rho
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer                                         :: i,iLoc(1) !JRE added iLoc(1) for MASK error
  integer                                         :: j, k
  integer                                         :: ref
  integer                                         :: k1, k2
  integer                                         :: nr_outside_steps
  integer                                         :: nr_inside_steps
  integer                                         :: i_min
  real(dp), dimension(ncomp,nphase)               :: xi_old
  real(dp), dimension(ncomp,nphase)               :: rhoi
  real(dp), dimension(ncomp,nphase)               :: lnphi
  real(dp), dimension(ncomp,nphase)               :: Ki
  real(dp)                                        :: error1, error2
  real(dp), dimension(nphase)                     :: phase_frac
  real(dp), dimension(nphase)                     :: phase_frac_old
  real(dp), dimension(nphase)                     :: delta_phase_frac
  real(dp)                                        :: increment

  real(dp)                                        :: q
  real(dp), dimension(ncomp)                      :: E
  real(dp), dimension(nphase)                     :: g
  real(dp), dimension(nphase,nphase)              :: H
  real(dp)                                        :: damp
  real(dp)                                        :: DET
  real(dp)                                        :: min_phase_frac
  real(dp), dimension(ncomp)                      :: psi

  real(dp), parameter                             :: tol_out = 1.E-10_dp
  real(dp), parameter                             :: tol_in = 1.E-13_dp
  !-----------------------------------------------------------------------------

  converg = .false.

  nr_outside_steps = 200
  nr_inside_steps = 10

  !-----------------------------------------------------------------------------
  ! initial value for fraction of all three phases
  !-----------------------------------------------------------------------------
  if ( nphase /= 3 ) then
     write (*,*) ' this part of the code is not general, but limited to 3 phases. Extend.'
     stop
  end if

  iLoc = MAXLOC( abs( xi(1:ncomp,3) - xi(1:ncomp,1) ) )
  i=iLoc(1)
  phase_frac(3) = ( xiF(i) - zi(i) ) / ( xi(i,3) - zi(i) )

  iLoc = MAXLOC( abs( xi(1:ncomp,2) - xi(1:ncomp,1) ) )
  i=iLoc(1)
  phase_frac(2) = ( zi(i) - xi(i,1) ) / ( xi(i,2) - xi(i,1) ) * ( 1.0_dp - phase_frac(3) )

  phase_frac(1) = 1.0_dp - phase_frac(2) - phase_frac(3)

  
  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'RR, N_out, N_inner, x1(1), xiF(1), x2(1), phase_frac, error1'

  !-----------------------------------------------------------------------------
  ! start iteration
  !-----------------------------------------------------------------------------

  xi_old( :, 1:nphase ) = xi( :, 1:nphase )

  error2 = 0.0_dp
  !comp_balance_violation = .false.

  k1 = 0
  error1 = tol_out + 1.0_dp
  do while ( error1 > tol_out .AND. k1 < nr_outside_steps )

     !--------------------------------------------------------------------------
     ! outer loop (converging the mole fractions)
     !--------------------------------------------------------------------------

     k1 = k1 + 1

     !--------------------------------------------------------------------------
     ! identify phase with highest molar phase fraction
     !--------------------------------------------------------------------------
     iLoc = MAXLOC( phase_frac(1:nphase) )
	 ref=iLoc(1)

     do j = 1, nphase
        call chemical_potential_tp ( t, p, xi(:,j), dense_start(j), rhoi(:,j), lnphi(:,j) )
     end do
     ! use converged value of eta as starting value for next iteration
     !dense_start(1) = get_eta( rhoi(1) )
     !dense_start(2) = get_eta( rhoi(2) )
     !dense_start(3) = get_eta( rhoi(3) )

     do j = 1, nphase
        Ki(1:ncomp,j) = EXP( lnphi(1:ncomp,ref) - lnphi(1:ncomp,j) )
     end do

     !--------------------------------------------------------------------------
     ! check for trivial solution
     !--------------------------------------------------------------------------
     !if ( ABS(1.0_dp-sum(rhoi1(:))/sum(rhoi2(:))) < 1.E-5_dp  &
     !     .AND. SUM(ABS(Ki(1:ncomp)-1.0_dp)) < 1.E-6_dp ) exit

     k2 = 0
     increment = 1.E-5_dp   ! increase diagonal elements for notorious cases, see book Michelsen
     error2 = tol_in + 1.0_dp
     do while ( error2 > tol_in .AND. k2 < nr_inside_steps )

        !-----------------------------------------------------------------------
        ! inner loop (converging the fraction of one phase)
        !-----------------------------------------------------------------------

        k2 = k2 + 1

        phase_frac_old(:) = phase_frac(:)

        !JRE forall( i = 1:ncomp ) E(i) = sum( phase_frac(1:nphase) / exp( lnphi(i,1:nphase) ) )
        Do i = 1,ncomp 
			E(i) = sum( phase_frac(1:nphase) / exp( lnphi(i,1:nphase) ) )
		enddo

        q = sum( phase_frac(:) ) - sum( xiF(:) * log( E(:) ) )

        do j = 1, nphase
           g(j) = 1.0_dp - sum( xiF(1:ncomp) / ( E(1:ncomp) * exp( lnphi(1:ncomp,j) ) ) )
           do k = 1, nphase
              H(j,k) = sum( xiF(1:ncomp) / ( E(1:ncomp) * E(1:ncomp) * exp( lnphi(1:ncomp,j) ) * exp( lnphi(1:ncomp,k) ) ) )
              if ( j == k ) H(j,k) = H(j,k) + increment
           end do
        end do

        delta_phase_frac(:) = - g(:)
        call MATINV ( nphase, 1, H, delta_phase_frac, DET )

        damp = 1.0_dp

        min_phase_frac = MINVAL( phase_frac_old(1:nphase) + delta_phase_frac(1:nphase) )
        if ( min_phase_frac < 0.0_dp ) then
           iLoc = MINLOC( phase_frac_old(1:nphase) + delta_phase_frac(1:nphase) )
		   i_min=iLoc(1)
           if ( abs( phase_frac_old(i_min) ) > 1.E-20 ) then
              damp = - phase_frac_old(i_min) * ( 1.0_dp - 1.E-10_dp ) / delta_phase_frac(i_min)
           else
              damp = 1.0E-4_dp
           end if
           increment = increment * 10.0_dp
           if ( outp >= 1 ) write (*,*) 'now reducing step size. Test this functionality !',damp
        end if

        !-----------------------------------------------------------------------
        ! new x vector
        !-----------------------------------------------------------------------
        phase_frac(1:nphase) = phase_frac_old(1:nphase) + damp * delta_phase_frac(1:nphase)

        error2 = sum( abs( phase_frac(1:nphase) - phase_frac_old(1:nphase) ) )
        if ( outp >= 3 ) write (*,'(a,4G20.8)') 'phase_fraction, error', phase_frac(1:3), error2

     end do

     !--------------------------------------------------------------------------
     ! constrain the phase fraction to ( 0 < phase_frac < 1 )
     !--------------------------------------------------------------------------
     !comp_balance_violation = .false.
     !if ( phase_frac > 1.0_dp ) comp_balance_violation = .true.
     !if ( phase_frac < 0.0_dp ) comp_balance_violation = .true.
     if ( k2 >= nr_inside_steps .AND. outp >= 2) write (*,*) 'three_phase_flash: Newton-loop not converged'

     !--------------------------------------------------------------------------
     ! determine x and error of outer loop
     !--------------------------------------------------------------------------
     do i = 1, ncomp
        psi(i) = 1.0_dp + sum( phase_frac(1:nphase) * ( Ki(i,1:nphase) - 1.0_dp ) )
     end do

     xi( 1:ncomp, ref ) = xiF(1:ncomp) / psi(1:ncomp)
     !JRE: forall( j = 1:nphase, j /= ref )
	 do j=1,nphase
        xi( 1:ncomp, j ) = xiF(1:ncomp) * Ki( 1:ncomp, j ) / psi(1:ncomp)
     enddo !JRE forall

     !JRE: forall( j = 1:nphase )
	 do j=1,nphase
        xi( 1:ncomp, j ) = xi( 1:ncomp, j ) / sum( xi( 1:ncomp, j ) )
     enddo !JRE forall

     error1 = 0.0_dp
     do j = 1, nphase
        error1 = error1 + sum( abs( xi_old(1:ncomp,j) - xi(1:ncomp,j) ) )
     end do
     xi_old(:,:) = xi(:,:)

     if ( outp >= 2 ) write (*,'(a,2i5,7G20.12)') 'RR', k1, k2, xi(1,1:nphase), phase_frac(1:3), error1

  end do


  !-----------------------------------------------------------------------------
  ! if convergence: accept solution
  !-----------------------------------------------------------------------------

  if ( error1 < tol_out .AND. error2 < tol_in ) then

     converg = .true.
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'three_phase_flash converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi( 1:ncomp,1 ) / sum( rhoi( 1:ncomp,1 ) )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi( 1:ncomp,2 ) / sum( rhoi( 1:ncomp,2 ) )
     if ( outp >= 2 ) write (*,*) '  x p3   = ',rhoi( 1:ncomp,3 ) / sum( rhoi( 1:ncomp,3 ) )
     if ( outp >= 2 ) write (*,*) 'eta 1,2,3 = ',get_eta( rhoi(:,1) ), get_eta( rhoi(:,2) ), get_eta( rhoi(:,3) )
     if ( outp >= 2 ) write (*,*) ' '
     do j = 1, nphase
        rho(j) = sum( abs( rhoi(1:ncomp,j) ) )
     end do

  end if

  if ( .NOT.converg .AND. outp >= 1) write (*,*) 'three_phase_flash: not converged'

end subroutine three_phase_flash


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine ps_flash
!> \brief Flash calculation for specified p, s, and composition xiF(:).
!!
!! This subroutine performs a flash calculation for a mixtue with defined p, s,
!! and composition xiF.
!! A stability analysis for a possible third phase or another stable second
!! phase is not performed
!!
!! input:    p, s, xiF(:)
!!           starting value for temperature T
!! output:   converg, T, phase_frac2 (mole-based fraction of phase 2)
!!           rhoi1, rhoi2 as converged component densities, if solution is found
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine ps_flash ( p, sF, xiF, tcalc, h, rhoiL, rhoiV, phase_frac2, converg )

  use properties, only: state_trho
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: p
  real(dp), intent(in)                            :: sF
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), intent(inout)                         :: tcalc
  real(dp), intent(out)                           :: h
  real(dp), dimension(ncomp), intent(out)         :: rhoiL
  real(dp), dimension(ncomp), intent(out)         :: rhoiV
  real(dp), intent(out)                           :: phase_frac2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer, parameter                              :: maxiter = 100
  integer                                         :: i_iter
  real(dp)                                        :: temp
  real(dp)                                        :: t_org
  real(dp)                                        :: t_step
  real(dp)                                        :: f_error
  real(dp)                                        :: dfdt
  real(dp), parameter                             :: f_tolerance = 1.E-10_dp
  real(dp), parameter                             :: delta_t = 1.E-4_dp
  real(dp)                                        :: pcalc
  real(dp)                                        :: v
  real(dp)                                        :: sV, sL
  real(dp)                                        :: hV, hL
  real(dp)                                        :: sF1, sF2
  logical                                         :: pt_converg
  !-----------------------------------------------------------------------------

  converg = .false.

  i_iter = 0
  temp = tcalc

  do while ( .NOT.converg .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     !--------------------------------------------------------------------------
     ! first evaluation of pT-flash, for difference scheme (Newton iteration)
     !--------------------------------------------------------------------------

     t_org = temp
     temp = t_org - delta_t
     call pt_flash ( temp, p, xiF, rhoiL, rhoiV, phase_frac2, pt_converg )

     if ( pt_converg ) then

        call state_trho ( temp, rhoiL, pcalc, v, sL, hL )
        call state_trho ( temp, rhoiV, pcalc, v, sV, hV )

        sF1 = phase_frac2 * sV + ( 1.0_dp - phase_frac2 ) * sL

     else

        exit

     end if

     !--------------------------------------------------------------------------
     ! second evaluation of pT-flash
     !--------------------------------------------------------------------------

     temp = t_org
     call pt_flash ( temp, p, xiF, rhoiL, rhoiV, phase_frac2, pt_converg )

     if ( pt_converg ) then

        call state_trho ( temp, rhoiL, pcalc, v, sL, hL )
        call state_trho ( temp, rhoiV, pcalc, v, sV, hV )

        sF2 = phase_frac2 * sV + ( 1.0_dp - phase_frac2 ) * sL

     else

        exit

     end if

     !--------------------------------------------------------------------------
     ! Newton step for iterating the temperature
     !--------------------------------------------------------------------------
     f_error = ( sF2 - sF )
     dfdt = ( sF2 - sF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  20.0_dp, t_step )
     t_step = max( -20.0_dp, t_step )

     temp = temp - t_step

     if ( ABS( f_error ) < f_tolerance ) then
        tcalc = temp
        h = phase_frac2 * hV + ( 1.0_dp - phase_frac2 ) * hL
        converg = .true.
     end if

  end do

  if ( i_iter >= maxiter ) write (*,*) 'ps_flash not converged'

end subroutine ps_flash



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine ph_flash
!> \brief Flash calculation for specified p, h, and composition xiF(:).
!!
!! This subroutine performs a flash calculation for a mixtue with defined p, h,
!! and composition xiF.
!! A stability analysis for a possible third phase or another stable second
!! phase is not performed
!!
!! input:    p, h, xiF(:)
!!           starting value for temperature T
!! output:   converg, T, phase_frac2 (vapor fraction)
!!           rhoi1, rhoi2 as converged component densities, if solution is found
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine ph_flash ( p, hF, xiF, tcalc, s, rhoiL, rhoiV, phase_frac2, converg )

  use properties, only: state_trho
  implicit none

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: p
  real(dp), intent(in)                            :: hF
  real(dp), dimension(ncomp), intent(in)          :: xiF
  real(dp), intent(inout)                         :: tcalc
  real(dp), intent(out)                           :: s
  real(dp), dimension(ncomp), intent(out)         :: rhoiL
  real(dp), dimension(ncomp), intent(out)         :: rhoiV
  real(dp), intent(out)                           :: phase_frac2
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  integer, parameter                              :: maxiter = 100
  integer                                         :: i_iter
  real(dp)                                        :: temp
  real(dp)                                        :: t_org
  real(dp)                                        :: t_step
  real(dp)                                        :: f_error
  real(dp)                                        :: dfdt
  real(dp), parameter                             :: f_tolerance = 1.E-9_dp
  real(dp), parameter                             :: t_tolerance = 1.E-10_dp
  real(dp), parameter                             :: delta_t = 1.E-3_dp
  real(dp)                                        :: pcalc
  real(dp)                                        :: v
  real(dp)                                        :: sV, sL
  real(dp)                                        :: hV, hL
  real(dp)                                        :: hF1, hF2
  logical                                         :: pt_converg
  !-----------------------------------------------------------------------------

  converg = .false.

  i_iter = 0
  temp = tcalc

  do while ( .NOT.converg .AND. i_iter <= maxiter )

     i_iter = i_iter + 1

     !--------------------------------------------------------------------------
     ! first evaluation of pT-flash, for difference scheme (Newton iteration)
     !--------------------------------------------------------------------------

     t_org = temp
     temp = t_org - delta_t
     call pt_flash ( temp, p, xiF, rhoiL, rhoiV, phase_frac2, pt_converg )

     if ( pt_converg ) then

        call state_trho ( temp, rhoiL, pcalc, v, sL, hL )
        call state_trho ( temp, rhoiV, pcalc, v, sV, hV )

        hF1 = phase_frac2 * hV + ( 1.0_dp - phase_frac2 ) * hL

     else

        exit

     end if

     !--------------------------------------------------------------------------
     ! second evaluation of pT-flash
     !--------------------------------------------------------------------------

     temp = t_org
     call pt_flash ( temp, p, xiF, rhoiL, rhoiV, phase_frac2, pt_converg )

     if ( pt_converg ) then

        call state_trho ( temp, rhoiL, pcalc, v, sL, hL )
        call state_trho ( temp, rhoiV, pcalc, v, sV, hV )

        hF2 = phase_frac2 * hV + ( 1.0_dp - phase_frac2 ) * hL

     else

        exit

     end if

     !--------------------------------------------------------------------------
     ! Newton step for iterating the temperature
     !--------------------------------------------------------------------------
     f_error = ( hF2 - hF )
     dfdt = ( hF2 - hF1 ) / delta_t

     t_step = f_error / dfdt
     t_step = min(  20.0_dp, t_step )
     t_step = max( -20.0_dp, t_step )

     temp = temp - t_step

     if ( ABS( f_error ) < f_tolerance .OR. abs( t_step ) < t_tolerance ) then
        tcalc = temp
        s = phase_frac2 * sV + ( 1.0_dp - phase_frac2 ) * sL
        converg = .true.
     end if

  end do

  if ( i_iter >= maxiter ) write (*,*) 'ph_flash not converged'

end subroutine ph_flash


end module flash
