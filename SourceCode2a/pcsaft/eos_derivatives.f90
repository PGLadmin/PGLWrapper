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

  real(dp)                                   :: z0
  real(dp)                                   :: z1
  real(dp)                                   :: z2
  real(dp)                                   :: z22
  real(dp)                                   :: z23
  real(dp)                                   :: z32
  real(dp)                                   :: z33
  real(dp)                                   :: z0t
  real(dp)                                   :: z1t
  real(dp)                                   :: z2t

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

  !-----------------------------------------------------------------------------
  ! pure and mixture quantities of fluid theory
  !-----------------------------------------------------------------------------
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

  public :: rho_independent_quantities, aggregate_polar_parameters, density_terms, F_density, f_rho, f_rho_rho,  &
       f_rho_4, F_density_rhok, F_density_rhok_T, F_density_rhok_rho,  &
       F_density_rhok_rhol, f_temp, f_temp_temp, f_temp_rho, allocate_eos_quantities,  &
       deallocate_eos_quantities, initialize_eos_quantities, Iterate_Association,  &
       cross_association_parameters
  real(dp), public                           :: z3          ! defined public
  real(dp), public                           :: rho         ! defined public
  real(dp), public                           :: z3t         ! defined public
  real(dp), allocatable, dimension(:), public  :: rhoi      ! defined public
  real(dp), public                           :: kT
  real(dp), allocatable, dimension(:), public :: dhs        ! defined public


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
     
  end if

  !-----------------------------------------------------------------------------
  ! constants and parameters of polar terms
  !-----------------------------------------------------------------------------

  ! ------ dipole-dipole term --------------------------------------------------
  if ( dipole ) then

     my_factor(:) = my_factor_no_T(:) / t

     pi_dd(:,:) = pi_dd_no_T(:,:) / t**2
     psi_dd(:,:,:) = psi_dd_no_T(:,:,:) / t**3

  end if
  
  ! ------ quadrupole-quadrupole term ------------------------------------------
  if ( qudpole ) then

     Q_factor(:) = Q_factor_no_T(:) / t

     pi_qq(:,:) = pi_qq_no_T(:,:) / t**2
     psi_qq(:,:,:) = psi_qq_no_T(:,:,:) / t**3

  end if
  
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
				 !JRE rewrite deprecated forall as simple do.
                 !forall( ii = 1:nhb_typ(i), jj = 1:nhb_typ(j), ii /= jj )
                 !   eps_hb(i,j,ii,jj) = 0.5_dp * ( eps_hb(i,i,ii,jj) + eps_hb(j,j,jj,ii) )  * ( 1.0_dp - kij_assoc(i,j) )
                 !end forall
				 do jj=1,nhb_typ(j)
					do ii=1,nhb_typ(i)
						eps_hb(i,j,ii,jj) = 0.5_dp * ( eps_hb(i,i,ii,jj) + eps_hb(j,j,jj,ii) )  * ( 1.0_dp - kij_assoc(i,j) )
					 enddo
				 enddo

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
  z_res = zhs + zhc + zdsp + zhb + zdd + zqq + zdq

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
     myres(k) = whs_rk + whc_rk + wdsp_rk + whb_rk + wdd_rk + wqq_rk + wdq_rk

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
				 !JRE  replaced deprecated forall with simple do.
                 !forall ( k = 1:ncomp,  l = 1:ncomp )
                 !end forall
				 DO l=1,ncomp
					DO k=1,ncomp
						q_rkrl(k,l) = q_rkrl(k,l) - factor_ij * gij_rkrl(k,l,i,j)
					enddo ! k
				enddo ! l
				!END JRE		
              end do
           end do
        end do
     end do


     !whb_rkrl = q_rkrl - q_Xr_transpose * inv( q_XX ) * q_Xr
     call MATINV( n_dim, m_dim, q_XX, q_Xr, determinant )  ! output q_Xr := inv( q_XX ) * q_Xr

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
           write (6,*) 'f_temp_rho: max_eval violated err_sum = ',err_sum,tol
           exit
        end if

     end do

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

  allocate ( compna( ncomp ) )
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
