!> \file critical_point.f90
!! \brief critical point of pure substances or of mixtures
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module module_critical_point
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module module_critical_point

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  private
  public :: critical_point_mixtures, critical_point_pure,  &
            critical_point_mixtures_binary_defined_t_or_p

CONTAINS

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine critical_point_mixtures
!> \brief Heidemann-Khalil-Michelsen scheme for critical point of mixtures
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine critical_point_mixtures ( zi, t, p, rhoi, converg )

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: zi
  real(dp), intent(inout)                         :: t
  real(dp), intent(out)                           :: p
  real(dp), dimension(ncomp), intent(out)         :: rhoi
  logical, intent(out)                            :: converg
  !-----------------------------------------------------------------------------

  integer                                         :: i
  integer                                         :: i_iter
  integer, parameter                              :: N_max_iter = 50
  real(dp), parameter                             :: tolerance = 1.E-5_dp
  real(dp), parameter                             :: d_t = 0.001_dp
  real(dp), parameter                             :: d_eta = 0.0001_dp
  real(dp), dimension(2)                          :: fct
  real(dp)                                        :: t_pert, eta_pert
  real(dp), dimension(2)                          :: fct_t, fct_e
  real(dp), dimension(2,2)                        :: H, inv_H
  real(dp)                                        :: det
  real(dp)                                        :: error
  real(dp), dimension(2)                          :: delta
  real(dp)                                        :: eta
  !-----------------------------------------------------------------------------

  converg = .false.

  eta = 0.14_dp

  i_iter = 0
  error = tolerance + 1.0_dp

  do while ( error > tolerance .AND. i_iter < N_max_iter )

     i_iter = i_iter + 1

     !--------------------------------------------------------------------------
     ! calculate residuals and derivative of residuals by numerical differentiation
     !--------------------------------------------------------------------------
     t_pert = t + d_t
     call critical_point_mixtures_objective ( t_pert, eta, zi, p, rhoi, fct_t )

     eta_pert = eta + d_eta
     call critical_point_mixtures_objective ( t, eta_pert, zi, p, rhoi, fct_e )

     call critical_point_mixtures_objective ( t, eta, zi, p, rhoi, fct )

     !--------------------------------------------------------------------------
     ! derivative of residuals
     !--------------------------------------------------------------------------

     H(1,1) = ( fct(1) - fct_t(1) ) / d_t
     H(2,1) = ( fct(1) - fct_e(1) ) / d_eta
     H(1,2) = ( fct(2) - fct_t(2) ) / d_t
     H(2,2) = ( fct(2) - fct_e(2) ) / d_eta

     det = H(1,1) * H(2,2) - H(2,1) * H(1,2)

     inv_H(1,1) = H(2,2)
     inv_H(1,2) = - H(1,2)
     inv_H(2,1) = - H(2,1)
     inv_H(2,2) = H(1,1)
     inv_H = inv_H / det

     do i = 1, 2
        delta(i) = sum( inv_H(1:2,i) * fct(1:2) )
     end do

     !--------------------------------------------------------------------------
     ! limit maximum step size
     !--------------------------------------------------------------------------
     if ( abs( delta(1) ) > 0.25_dp * t .OR. abs( delta(2) ) > 0.015_dp ) delta(:) = delta(:) * 0.5_dp

     delta(1) = min(   0.3_dp * t, delta(1) )
     delta(1) = max( - 0.3_dp * t, delta(1) )
     delta(2) = min(   0.01_dp, delta(2) )
     delta(2) = max( - 0.01_dp, delta(2) )

     !--------------------------------------------------------------------------
     ! Newton step for critical point Tc and eta_c
     !--------------------------------------------------------------------------
     t = t + delta(1)
     eta = eta + delta(2)

     eta = max( eta, 0.0001_dp )

     error = sum( abs( fct(1:2) ) )

     if ( error < tolerance ) converg = .true.

     ! write (*,*) 'step, t, eta  ', delta(1), delta(2)
     ! write (*,*) 'objective fct ', fct
     ! write (*,*) 'values, t, eta', t, eta
     ! write (*,*) 'error ',error, i_iter
     ! read (*,*)

  end do

end subroutine critical_point_mixtures



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine critical_point_mixtures_binary_defined_t_or_p
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine critical_point_mixtures_binary_defined_t_or_p ( iterate_t, t, p, zi, rhoi, converg )

   !-----------------------------------------------------------------------------
   logical, intent(in)                             :: iterate_t
   real(dp), intent(inout)                         :: t
   real(dp), intent(inout)                         :: p
   real(dp), dimension(ncomp), intent(inout)       :: zi
   real(dp), dimension(ncomp), intent(out)         :: rhoi
   logical, intent(out)                            :: converg
   !-----------------------------------------------------------------------------
 
   integer                                         :: i_iter
   integer, parameter                              :: N_max_iter = 50
   real(dp), parameter                             :: tolerance = 1.E-7_dp
   real(dp)                                        :: d_z
   real(dp)                                        :: t1, t2
   real(dp)                                        :: p1, p2
   real(dp), dimension(ncomp)                      :: zi_aux
   real(dp), dimension(ncomp)                      :: rhoi_aux
   real(dp)                                        :: error, rel_error
   real(dp)                                        :: dt_dz
    !-----------------------------------------------------------------------------
 
   converg = .false.

   d_z = 0.001_dp
   i_iter = 0
   rel_error = tolerance + 1.0_dp
 
   do while ( rel_error > tolerance .AND. i_iter < N_max_iter )
 
      i_iter = i_iter + 1
 
      !--------------------------------------------------------------------------
      ! calculate residuals and derivative of residuals by numerical differentiation
      !--------------------------------------------------------------------------
      t1 = t
      p1 = p
      call critical_point_mixtures ( zi, t1, p1, rhoi, converg )
 
      t2 = t1
      p2 = p1
      if ( ( zi(1) + d_z ) >= 1.0_dp ) d_z = - d_z
      zi_aux(1) = zi(1) + d_z
      zi_aux(2) = 1.0_dp - zi_aux(1)
      call critical_point_mixtures ( zi_aux, t2, p2, rhoi_aux, converg )
 
      !--------------------------------------------------------------------------
      ! Newton step
      !--------------------------------------------------------------------------
 
      if ( .NOT.iterate_t ) then
         error = t1 - t
         rel_error = abs( error ) / t
         dt_dz = ( t2 - t1 ) / d_z
      else
         error = p1 - p
         rel_error = abs( error ) / p
         dt_dz = ( p2 - p1 ) / d_z
      end if

      zi(1) = zi(1) - error / dt_dz
      zi(2) = 1.0_dp - zi(1)

      if ( rel_error < tolerance .OR. abs( error / dt_dz ) < 1.E-8 ) then
         converg = .true.
         t = t1
         p = p1
      end if
 
      ! write (*,*) 'i, rel_error ', i_iter, rel_error, error / dt_dz
      ! write (*,*) 'zi_old(1), zi(1)', zi(1) + error / dt_dz, zi(1)
      ! read (*,*)
 
   end do

   if ( .NOT.converg ) write (*,*) 'critical_point_mixtures_binary_defined_t_or_p: calc. not converged'
 
 end subroutine critical_point_mixtures_binary_defined_t_or_p
 
 
 
 !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE critical_point_mixtures_objective
!> \brief critical_point_mixtures_objective
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine critical_point_mixtures_objective ( tF, eta, zF, pF, rhoi, obj_fct )

  use properties, only: chemical_potential_trho, pressure
  use module_eos_derivatives, only: z3t, rho_independent_quantities

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: tF
  real(dp), intent(in)                            :: eta
  real(dp), dimension(ncomp), intent(in)          :: zF
  real(dp), intent(out)                           :: pF
  real(dp), dimension(ncomp), intent(out)         :: rhoi
  real(dp), dimension(2), intent(out)             :: obj_fct

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: pcalc_z
  real(dp), dimension(ncomp)                      :: lnphi
  real(dp), dimension(ncomp,ncomp)                :: w_rkrl, wig_rkrl
  real(dp), dimension(ncomp,ncomp)                :: qij
  real(dp), dimension(ncomp)                      :: mu
  real(dp), dimension(ncomp)                      :: mu_res
  real(dp)                                        :: lamda
  real(dp), dimension(ncomp)                      :: eigenvector
  
  integer, parameter                              :: it_max = 200
  integer                                         :: it_num
  integer                                         :: rot_num
  real(dp), dimension(ncomp,ncomp)                :: eigenvec
  real(dp), dimension(ncomp)                      :: eigenval

  real(dp)                                        :: epsilon
  real(dp), dimension(ncomp)                      :: rhoi_eps
  real(dp)                                        :: d3F_ds3
  real(dp)                                        :: dF_ds_1, dF_ds_2
  real(dp), dimension(ncomp)                      :: mu_eps
  real(dp), dimension(ncomp)                      :: g_function
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! obtain parameters and density independent expressions
  !-----------------------------------------------------------------------------

  call rho_independent_quantities ( tF, zF )

  !-----------------------------------------------------------------------------
  ! density iteration: (pTx)-ensemble
  !-----------------------------------------------------------------------------

  call pressure ( eta, pF, pcalc_z )

  rhoi(1:ncomp) = eta / z3t * zF(1:ncomp)


  !-----------------------------------------------------------------------------
  ! derivative of d2(F) / (d(rhoi)*d(rhoj)) = w_rkrl + wig_rkrl,  with F as the tpd
  ! The lowest eigenvalue of that matrix is zero at a critical point
  ! (and at any spinodal point)
  !-----------------------------------------------------------------------------

  call chemical_potential_trho ( tF, rhoi, lnphi, mu, w_rkrl=w_rkrl, wig_rkrl=wig_rkrl )
  do i = 1, ncomp
     qij(i,1:ncomp) = ( w_rkrl(i,1:ncomp) + wig_rkrl(i,1:ncomp) ) * sqrt( rhoi(i) * rhoi(1:ncomp) )
  end do

  call jacobi_eigenvalue ( ncomp, qij, it_max, eigenvec, eigenval, it_num, rot_num )

  lamda = eigenval(1)
  eigenvector(:) = eigenvec(:,1)

  !-----------------------------------------------------------------------------
  ! numerical derivative of d3(F)/d(s)**3 where F is the tpd
  ! derivative is calculated as d3(F)/d(s)**3 = ( d(F)/d(s)!(s=-eps) + d(F)/d(s)!(s=+eps) ) / eps**2
  ! d(F)/d(s)!(s) is calculated as d(F)/d(s)!(s) = sum( d(F)/d(rhoi) *sqrt(rhoi) * eigenvector(:) )
  ! with direction vector d(rhoi(:)) = sqrt(rhoi) * d( s*eigenvector(:) ) = sqrt(rhoi)*eigenvector(:) * d(s)
  ! note sqrt(rhoi) is a scalar
  !-----------------------------------------------------------------------------

  epsilon = - 1.E-4_dp
  rhoi_eps(1:ncomp) = rhoi(1:ncomp) + sqrt( rhoi(1:ncomp) ) * eigenvector(1:ncomp) * epsilon
  call chemical_potential_trho ( tF, rhoi_eps, mu_res, mu_eps )
  g_function(1:ncomp) = ( mu_eps(1:ncomp) - mu(1:ncomp) ) * sqrt( rhoi(1:ncomp) )
  dF_ds_1 = sum( eigenvector(1:ncomp) * g_function(1:ncomp) )

  epsilon = + 1.E-4_dp
  rhoi_eps(1:ncomp) = rhoi(1:ncomp) + sqrt( rhoi(1:ncomp) ) * eigenvector(1:ncomp) * epsilon
  call chemical_potential_trho ( tF, rhoi_eps, mu_res, mu_eps )
  g_function(1:ncomp) = ( mu_eps(1:ncomp) - mu(1:ncomp) ) * sqrt( rhoi(1:ncomp) )
  dF_ds_2 = sum( eigenvector(1:ncomp) * g_function(1:ncomp) )

  d3F_ds3 = ( dF_ds_1 + dF_ds_2 ) / epsilon**2

  !-----------------------------------------------------------------------------
  ! objective function for critical point
  !-----------------------------------------------------------------------------

  obj_fct(1) = lamda * 10.0_dp
  obj_fct(2) = d3F_ds3

end subroutine critical_point_mixtures_objective


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine critical_point_pure
!> \brief determine critical point of pure substance
!!
!! This subroutine determines the critical point of a pure substance.
!! The density and temperature of the pure component are iterated
!! until the derivative of (dP/d_rho) and the second derivative
!! (dP2/d2_rho) are zero.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine critical_point_pure ( tc, converg, pc, rhoc )

  use properties, only: pressure
  use module_eos_derivatives, only: z3t, rho_independent_quantities
  !-----------------------------------------------------------------------------
  real(dp), intent(inout)                         :: tc
  logical, intent(out)                            :: converg
  real(dp), intent(out)                           :: pc
  real(dp), intent(out)                           :: rhoc
  !-----------------------------------------------------------------------------

  integer, parameter                              :: max_eta_count = 200
  real(dp), parameter                             :: t_tol = 0.005_dp  ! tolerance derivative d(p)/d(t)
  real(dp), parameter                             :: eta_tol = 0.2_dp  ! tolerance derivative d2(p)/d(t)d(t)
  real(dp), parameter                             :: tdelt = 0.1_dp

  integer                                         :: count
  integer                                         :: count_t
  real(dp)                                        :: tF
  real(dp), dimension(ncomp)                      :: xF
  real(dp)                                        :: pdz1
  real(dp)                                        :: eta_step
  real(dp)                                        :: t_step
  real(dp)                                        :: eta
  real(dp)                                        :: pcalc
  real(dp)                                        :: pcalc_z
  real(dp)                                        :: pcalc_z2
  real(dp)                                        :: pcalc_z3
  !-----------------------------------------------------------------------------

  eta = 0.2_dp                                    ! starting value density (i.e. packing fraction)
  tF = tc

  xF(1) = 1.0_dp

  count_t = 0
  pcalc_z = t_tol + 1.0_dp

  do while ( ABS( pcalc_z ) > t_tol .AND. count_t < 30 )

     tF = tF - tdelt
     call rho_independent_quantities ( tF, xF )

     count_t = count_t + 1

     count = 0
     pcalc_z2 = eta_tol + 1.0_dp

     do while ( ABS(pcalc_z2) > eta_tol .AND. count < max_eta_count )

        count = count + 1
        call pressure ( eta, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

        eta_step = pcalc_z2 / pcalc_z3
        if ( abs( eta_step ) > 0.1_dp ) eta_step = eta_step / abs( eta_step ) * 0.1_dp
        eta = eta - eta_step
        IF ( eta < 0.0_dp ) eta = 0.001_dp
        IF ( eta > 0.7_dp ) eta = 0.5_dp
        ! write(*,'(a,i4,4(E16.7))') 'c_1', count, eta, pcalc, pcalc_z2, pcalc_z3

     end do
     pdz1 = pcalc_z



     tF = tF + tdelt
     call rho_independent_quantities ( tF, xF )

     count = 0
     pcalc_z2 = eta_tol + 1.0_dp

     do while ( abs( pcalc_z2 ) > eta_tol .and. count < max_eta_count )

        count = count + 1
        call pressure ( eta, pcalc, pcalc_z, pcalc_z2, pcalc_z3 )

        eta_step = pcalc_z2 / pcalc_z3
        if ( abs( eta_step ) > 0.1_dp ) eta_step = eta_step / abs( eta_step ) * 0.1_dp
        eta = eta - eta_step
        if ( eta < 0.0_dp ) eta = 0.001_dp
        if ( eta > 0.7_dp ) eta = 0.5_dp

     end do

     t_step = pcalc_z / ( ( pcalc_z-pdz1 ) / tdelt )
     if ( ABS( t_step ) >= 0.05_dp * tF ) t_step = t_step / abs( t_step ) * ( 0.05_dp * tF )
     tF = tF - t_step

     if ( tF < 0.0_dp ) tF = 600.0_dp

  end do

  if ( ABS( pcalc_z ) < t_tol .AND. ABS( pcalc_z2 ) < eta_tol ) then
     converg = .true.
     tc = tF
     pc = pcalc
     rhoc = eta / z3t
     ! zc = p /( kbol30 * t * rho )
  else
     converg = .false.
     tc = 600.0_dp
     pc = 1000.0_dp
     rhoc = 0.05_dp
     if ( ABS( pcalc_z ) > t_tol )  write (*,*) 'critical_point_pure: T-iteration not converged'
     if ( ABS( pcalc_z2 ) > eta_tol ) write (*,*) 'critical_point_pure: density-iteration not converged'
  end if
     

end subroutine critical_point_pure


end module module_critical_point

