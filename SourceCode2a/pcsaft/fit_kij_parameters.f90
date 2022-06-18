!> \file kij-fitting.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE kij_fitting
!!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module kij_fitting

  use PARAMETERS, only: dp, machine_eps
  use BASIC_VARIABLES, only: ncomp
  implicit none
  save

  private
  integer, parameter                              :: outp = 0
  integer                                         :: n_exp
  integer                                         :: adjust_parameter_option
  integer                                         :: objective_type
  real(dp), allocatable, dimension(:)             :: x_exp, y_exp, t_exp, p_exp
  real(dp), allocatable, dimension(:,:)           :: deviation
  real(dp), parameter                             :: lij_exponent = 3.0_dp
  real(dp), allocatable, dimension(:,:)           :: dev_previous
  real(dp), allocatable, dimension(:,:)           :: par_previous
  logical, parameter                              :: get_start_val = .true.

  public :: kij_lij_fitting

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE KIJ_FIT_LIJ
!> \brief routine for fitting kij (and lij, if needed) to exp. binary data
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine kij_lij_fitting

  !use optimizer_2D
  use BASIC_VARIABLES, only: compna_input
  use eos_constants
  use module_eos_derivatives, only: allocate_eos_quantities, initialize_eos_quantities,  &
       cross_association_parameters, aggregate_polar_parameters, deallocate_eos_quantities
  use pcsaft_pure_and_binary_parameters, only: load_pure_and_binary_parameters, kij, lij,  &
       kij_assoc, kij_t, lij_correction, kij_T_dependent, file_open, assoc
  use input_output, only: u_in_p
  use minimizer, only: cg_optimizer, bfgs, steepest_descent

  !-----------------------------------------------------------------------------
  integer                                         :: index_read_loop
  integer                                         :: i, nadjust, PRIN
  integer                                         :: read_info
  real(dp), allocatable, dimension(:)             :: optpars
  real(dp)                                        :: fmin, t0, h0, MACHEP
  real(dp)                                        :: dummy2, dummy3, dummy4, dummy5
  real(dp)                                        :: rms1, rms2, rms3, aad_m1, aad_m2, aad_m3
  character                                       :: fitin*50, fitin2*71, fitout*80, dummy1*30

  ! integer                                       :: STATUS,  iter, nfunc, ngrad
  ! real(dp)                                      :: gnorm
  ! real(dp), allocatable, dimension(:)           :: d, g, xtemp, gtemp
  real(dp), allocatable, dimension(:)             :: g
  integer                                         :: maxiter

  integer                                         :: iflag
  real(dp), parameter                             :: xtol = 1.E-4_dp
  real(dp), parameter                             :: gtol = 1.E-6_dp
  real(dp), parameter                             :: ftol = 1.E-5_dp
  integer, parameter                              :: iprint = 2
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! forget what was read from INPUT-file, because the input is given by user
  !-----------------------------------------------------------------------------
  call deallocate_eos_quantities
  deallocate( compna_input )

  !-----------------------------------------------------------------------------
  ! define some quantities
  !-----------------------------------------------------------------------------
  u_in_P = 1.0_dp

  ncomp = 2
  allocate( compna_input( ncomp ) )

  write (*,*) ' CHOOSE THE INPUT FILE (./input_file/kij-dat/...) :'
  write (*,*) ' SPECIFY FULL NAME:'
  read (*,*) fitin
  ! fitin = 'h2o-co2.dat'
  ! fitin = 'h2o-acetic-acid.dat'
  ! fitin = '1-propanol-water.dat'
  ! fitin = 'h2o-nh3.dat'
  ! fitin = 'butanone-hexane.dat'
  ! fitin = 'butanone-hexane_TAMie015.dat'
  ! fitin = 'butanone-hexane_TAMie0107.dat'
  ! fitin = 'butanone-hexane_TAMie0000.dat'
  ! fitin = 'co2-heptane.dat'
  ! fitin = 'methane-tetracontane.dat'
  ! fitin = 'n2-acetone.dat'
  ! fitin = 'meoh-octane-lle.dat'
  ! fitin = 'ethygly-hexane-lle.dat'
  ! fitin = 'co2-benzene_all.dat'
  ! fitin = '2-propanol-benzene_318K.dat'
  ! fitin = 'water-nitrogen.dat'
  fitin2 = './input_file/kij-dat/'//fitin

  !-----------------------------------------------------------------------------
  ! read experimental data (phase equilibrium data)
  !-----------------------------------------------------------------------------
  do index_read_loop = 1, 2
     call file_open( 87, fitin2 )

     i = 0
     read (87,*) compna_input(1)
     read (87,*) compna_input(2)

     read_info = 0
     do while ( read_info == 0 )
        read ( 87, *, IOSTAT=read_info ) dummy1, dummy2, dummy3, dummy4, dummy5
        !------ pure components are excluded, because kij not sensitive --------
        if ( read_info == 0 .AND.  &
             .NOT.( abs(dummy2-1.0_dp) < machine_eps .or. abs(dummy3-1.0_dp) < machine_eps ) .AND.  &
             .NOT.( abs(dummy2) < machine_eps .and. abs(dummy3) < machine_eps ) ) then
           i = i + 1
           if ( index_read_loop == 2 ) then
              x_exp(i) = dummy2
              y_exp(i) = dummy3
              p_exp(i) = dummy4 * 1.E5_dp
              t_exp(i) = dummy5 + 273.15_dp
           end if
        end if
     end do
     close (87)
     n_exp = i
     if ( index_read_loop == 1 ) then
        allocate ( x_exp(n_exp), y_exp(n_exp), t_exp(n_exp), p_exp(n_exp), deviation(n_exp,3) )
        deviation(:,:) = 0.0_dp
     end if
  end do

  !-----------------------------------------------------------------------------
  ! select type of objective function
  !-----------------------------------------------------------------------------
  write (*,*) ' SPECIFY THE OBJECTIVE FUNCTION :'
  write (*,*) '      1: regular'
  write (*,*) '      2: (p/p_exp-1) for liquid, only'
  write (*,*) '      3: distance criterion'
  write (*,*) '      4: chem. pot. criterion (only for complete (T,p,x,y)-data sets)'
  read (*,*) objective_type

  !-----------------------------------------------------------------------------
  ! define adjustable parameters
  !-----------------------------------------------------------------------------
  write (*,*) ' CHOOSE THE BINARY PARAMETERS TO BE ADJUSTED :'
  write (*,*) '      1:    kij'
  write (*,*) '      2:    kij AND lij'
  write (*,*) '      3:    kij AND kij_assoc'
  write (*,*) '      4:    kij AND lij AND kij_assoc'
  write (*,*) '      5:    kij AND lij asymmetric to lji'
  write (*,*) '      6:    kij = a + b * (T/K-300)/100'
  write (*,*) '      7:    kij = a + b * (T/K-300)/100 AND kij_assoc'
  write (*,*) '      8:    kij_assoc'
  read (*,*) adjust_parameter_option

  !-----------------------------------------------------------------------------
  ! set the number of adjustable parameters
  !-----------------------------------------------------------------------------
  if ( adjust_parameter_option == 1 ) then
     nadjust = 1
  else if ( adjust_parameter_option == 2 ) then
     nadjust = 2
  else if ( adjust_parameter_option == 3 ) then
     nadjust = 2
  else if ( adjust_parameter_option == 4 ) then
     nadjust = 3
  else if ( adjust_parameter_option == 5 ) then
     nadjust = 3
  else if ( adjust_parameter_option == 6 ) then
     nadjust = 2
   else if ( adjust_parameter_option == 7 ) then
      nadjust = 3
   else if ( adjust_parameter_option == 8 ) then
      nadjust = 1
    else
     write (*,*) 'make correct choice. stop.'
     stop
  end if

  !-----------------------------------------------------------------------------
  ! allocate ad initialize parameter-arrays and data-arrays
  !-----------------------------------------------------------------------------
  allocate ( optpars(nadjust) )
  allocate ( dev_previous(n_exp, nadjust+1) )
  allocate ( par_previous(nadjust, nadjust+1) )
  dev_previous(:,:) = 0.0_dp
  par_previous(:,:) = 0.0_dp

  fitout = fitin
  fitout = './output_file/kij_pcsaft/'//trim( adjustl( fitout ) )
  open ( 88, FILE = fitout )

  !-----------------------------------------------------------------------------
  ! allocate and initialize eos-array quantities
  !-----------------------------------------------------------------------------
  call allocate_eos_quantities
  call initialize_eos_quantities

  !-----------------------------------------------------------------------------
  ! read pure component and mixture parameters of PC-SAFT EOS
  !-----------------------------------------------------------------------------
  write (*,*) compna_input(1), '  ', compna_input(2)
  call load_pure_and_binary_parameters
  call cross_association_parameters
  call aggregate_polar_parameters

  !-----------------------------------------------------------------------------
  ! EOS constants
  !-----------------------------------------------------------------------------
  call load_eos_constants


  !-----------------------------------------------------------------------------
  ! starting value for kij
  !-----------------------------------------------------------------------------
  write (*,'(a,F12.5)') ' current value of kij =     ', kij(1,2)
  write (*,'(a,F12.5)') ' current value of lij =     ', lij(1,2)
  if ( assoc ) write (*,'(a,F12.5)') ' current value of kij_assoc=', kij_assoc(1,2)

  write (*,'(a)') ' do you want to give another starting value for kij? ( 1 : yes, 0 : no)'
  read (*,*) i

  if ( i == 1 ) then
     write (*,*) 'give starting value (or fixed value) for kij'
     read (*,*) kij(1,2)
     kij(2,1) = kij(1,2)
     write (*,*) 'give starting value (or fixed value) for lij'
     read (*,*) lij(1,2)
     lij(2,1) = - lij(1,2)
     if ( abs( lij(1,2) ) > machine_eps ) lij_correction = .true.
     if ( assoc ) then
        write (*,*) 'give starting value (or fixed value) for kij_assoc'
        read (*,*) kij_assoc(1,2)
        kij_assoc(2,1) = kij_assoc(1,2)
     end if
  end if

  if ( ( adjust_parameter_option == 2 .OR. adjust_parameter_option == 4 .OR.  &
       adjust_parameter_option == 5 ) .AND. abs( lij(1,2) ) < 1.E-5_dp ) then
     write (*,*) 'a non-zero starting value for lij has to be used. Try positve, then negative.'
     stop
  end if

  !-----------------------------------------------------------------------------
  ! numerical constants for optimizer
  !-----------------------------------------------------------------------------
  t0 = 1.E-4_dp
  h0 = 0.01_dp
  PRIN = 0
  MACHEP = 1.E-15_dp

  !-----------------------------------------------------------------------------
  ! write starting values to array of parameters for optimization
  !-----------------------------------------------------------------------------
  optpars(1) = kij(1,2)  + 1.0_dp
  if ( adjust_parameter_option == 2 ) then
     optpars(2) = sign( abs(lij(1,2))**(1.0_dp/lij_exponent), lij(1,2) ) + 1.0_dp
     lij_correction = .true.
  else if ( adjust_parameter_option == 3 ) then
     optpars(2) = kij_assoc(1,2) + 1.0_dp
  else if ( adjust_parameter_option == 4 ) then
     optpars(2) = sign( abs(lij(1,2))**(1.0_dp/lij_exponent), lij(1,2) ) + 1.0_dp
     lij_correction = .true.
     optpars(3) = kij_assoc(1,2) + 1.0_dp
  else if ( adjust_parameter_option == 5 ) then
     lij_correction = .true.
     optpars(2) = sign( abs(lij(1,2))**(1.0_dp/lij_exponent), lij(1,2) ) + 1.0_dp
     optpars(3) = sign( abs(lij(1,2))**(1.0_dp/lij_exponent), - lij(1,2) ) + 1.0_dp
  else if ( adjust_parameter_option == 6 ) then
     kij_T_dependent = .true.
     optpars(2) = 1.0_dp
   else if ( adjust_parameter_option == 7 ) then
      kij_T_dependent = .true.
      optpars(2) = 1.0_dp
      optpars(3) = kij_assoc(1,2) + 1.0_dp
   else if ( adjust_parameter_option == 8 ) then
      optpars(1) = kij_assoc(1,2) + 1.0_dp   ! overwrite first line above
    end if


  !-----------------------------------------------------------------------------
  ! call optimizer, minimizing objective function (deviation of calc. to exp data)
  !-----------------------------------------------------------------------------
  if ( nadjust == 1 ) then

     allocate ( g(nadjust) )
     call Newton_Opt_2D ( kij_objective_fct, optpars, nadjust, 1.E-8_dp, 1.E-8_dp, g, fmin)
     deallocate ( g )

  else
     
     maxiter = 1
     call steepest_descent ( nadjust,optpars,fmin,xtol,gtol,ftol,maxiter,iprint,iflag,kij_objective_fct ) !, grad )

     ! call bfgs_optimize_binary_par ( nadjust, optpars, fmin )
     maxiter = 7
     call bfgs ( nadjust,optpars,fmin,xtol,gtol,ftol,maxiter,iprint,iflag,kij_objective_fct ) !,kij_grad)
     ! write (*,*) ' now doing CG, compare convergence characteristics !!'
     ! call cg_optimizer ( nadjust,optpars,fmin,xtol,gtol,ftol,maxiter,iprint,iflag,kij_objective_fct ) !,kij_grad)

  end if

  
  !-----------------------------------------------------------------------------
  ! assign array of optimized parameters back to binary parameters kij, ...
  !-----------------------------------------------------------------------------
  if ( adjust_parameter_option /= 8 ) then
     kij(1,2) = optpars(1) - 1.0_dp
  end if
  if ( adjust_parameter_option == 2 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
  else if ( adjust_parameter_option == 3 ) then
     kij_assoc(1,2) = optpars(2) - 1.0_dp
  else if ( adjust_parameter_option == 4 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
     kij_assoc(1,2) = optpars(3) - 1.0_dp
  else if ( adjust_parameter_option == 5 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
     lij(2,1) = ( optpars(3) - 1.0_dp )**lij_exponent
  else if ( adjust_parameter_option == 6 ) then
     kij_t(1,2) = optpars(2) - 1.0_dp
  else if ( adjust_parameter_option == 7 ) then
     kij_t(1,2) = optpars(2) - 1.0_dp
     kij_assoc(1,2) = optpars(3) - 1.0_dp
  else if ( adjust_parameter_option == 8 ) then
      kij_assoc(1,2) = optpars(1) - 1.0_dp
  end if


  !-----------------------------------------------------------------------------
  ! output results
  !-----------------------------------------------------------------------------
  write (*,*) ' Result  kij=',kij(1,2),'  AAD=', fmin / real( n_exp, KIND=dp )
  write(88,*) ' Result  kij=',kij(1,2),'  AAD=', fmin / real( n_exp, KIND=dp )
  if ( adjust_parameter_option == 6 .OR. adjust_parameter_option == 7 ) then
     write (*,*) ' Result  beta =', kij_t(1,2),' with kij(T) = kij + beta * ( T/K - 300 ) / 100'
     write(88,*) ' Result  beta =', kij_t(1,2),' with kij(T) = kij + beta * ( T/K - 300 ) / 100'
  end if
  write (*,*) ' Result  lij=', lij(1,2)
  write(88,*) ' Result  lij=', lij(1,2)
  if ( assoc ) write (*,*) ' Result  kij_assoc=', kij_assoc(1,2)
  if ( assoc ) write(88,*) ' Result  kij_assoc=', kij_assoc(1,2)


  rms1 = 0.0_dp
  rms2 = 0.0_dp
  aad_m1 = 0.0_dp
  aad_m2 = 0.0_dp
  write (88,*) 'experimental data'
  do i = 1, n_exp
     write (88,'(i4,4(2x,G13.6))') i, x_exp(i), y_exp(i), p_exp(i)/1.E5_dp, t_exp(i)-273.15_dp
  end do
  write (88,*) '  '
  write (88,*) 'calculated values'
  do i = 1, n_exp
     write (88,'(i4,4(2x,G13.6))') i, x_exp(i) +deviation(i,1), y_exp(i) +deviation(i,2), p_exp(i)/1.E5_dp, t_exp(i)-273.15_dp
  end do
  write (88,*) '  '
  do i = 1, n_exp
     write (*,'(3(a,G20.10))')  ' DEV% x_ph_1=',100.0_dp*deviation(i,1),' DEV% x_ph_2=',100.0_dp*deviation(i,2),  &
                                ' DEV% p=',100.0_dp*deviation(i,3)
     write (88,'(3(a,G20.10))') ' DEVIATION% x_phase1=',100.0_dp*deviation(i,1), &
                                ' DEVIATION% x_phase2=',100.0_dp*deviation(i,2),' DEVIATION% p=',100.0_dp*deviation(i,2)
     rms1 = rms1 + deviation(i,1)**2
     rms2 = rms2 + deviation(i,2)**2
     rms3 = rms3 + deviation(i,3)**2
     aad_m1 = aad_m1 + abs( deviation(i,1) )
     aad_m2 = aad_m2 + abs( deviation(i,2) )
     aad_m3 = aad_m3 + abs( deviation(i,3) )
  end do
  rms1 = 100.0_dp * ( rms1 / real( n_exp, KIND=dp ) )**0.5_dp
  rms2 = 100.0_dp * ( rms2 / real( n_exp, KIND=dp ) )**0.5_dp
  rms3 = 100.0_dp * ( rms3 / real( n_exp, KIND=dp ) )**0.5_dp
  aad_m1 = 100.0_dp * aad_m1 / real( n_exp, KIND=dp )
  aad_m2 = 100.0_dp * aad_m2 / real( n_exp, KIND=dp )
  aad_m3 = 100.0_dp * aad_m3 / real( n_exp, KIND=dp )

  write(*,*)' '
  write(*,*)' AAD% x_phase_1 =', aad_m1
  write(*,*)' AAD% x_phase_2 =', aad_m2
  write(*,*)' AAD%         p =', aad_m3
  write(*,*)' '
  write(*,*)' RMS% x_phase_1 =', rms1
  write(*,*)' RMS% x_phase_2 =', rms2
  write(*,*)' RMS%         p =', rms3

  write(88,*)' '
  write(88,*)' AAD% x_phase_1 =', aad_m1
  write(88,*)' AAD% x_phase_2 =', aad_m2
  write(88,*)' AAD%         p =', aad_m3
  write(88,*)' '
  write(88,*)' RMS% x_phase_1 =', rms1
  write(88,*)' RMS% x_phase_2 =', rms2
  write(88,*)' RMS%         p =', rms3

  write (*,*) ' '
  write(*,*)' Result  kij=', kij(1,2)
  write(*,*)'         lij=', lij(1,2)
  if ( assoc ) write(*,*)'         kij_assoc=', kij_assoc(1,2)
  write(88,*)' Result  kij=', kij(1,2)
  write(88,*)'         lij=', lij(1,2)
  if ( assoc ) write(88,*)'         kij_assoc=', kij_assoc(1,2)
  write (*,*) ' '
  write (88,*) ' '

  !-----------------------------------------------------------------------------
  ! close file and deallocate
  !-----------------------------------------------------------------------------
  close (88)
  deallocate ( x_exp, y_exp, t_exp, p_exp, deviation )
  deallocate ( optpars )
  deallocate ( dev_previous, par_previous )
  ! DEALLOCATE( d, g, xtemp, gtemp )

end subroutine kij_lij_fitting




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE kij_objective_fct
!> \brief kij_objective_fct
!!
!! \todo Detailed description.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine kij_objective_fct ( nadjust, optpars, fmin )

  use pcsaft_pure_and_binary_parameters, only: kij, lij, kij_assoc, assoc, kij_t
  use module_eos_derivatives, only: cross_association_parameters

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: nadjust
  real(dp), dimension(nadjust), intent(in)        :: optpars
  real(dp), intent(in out)                        :: fmin

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: dev_x_i, dev_y_i, weigh_x_p, weigh_x_y
  real(dp), allocatable, dimension(:)             :: dev
  real(dp)                                        :: p_trial
  real(dp)                                        :: error
  real(dp)                                        :: distance
  real(dp)                                        :: distanceL, distanceV
  real(dp)                                        :: delta_dev, delta_par
  real(dp)                                        :: penalty
  real(dp)                                        :: kij_bound
  real(dp), dimension(2)                          :: exp_vec
  real(dp), dimension(ncomp)                      :: residual
  logical                                         :: converg
  logical                                         :: converg_bbp
  logical                                         :: convergL, convergV
  character (LEN=1)                               :: report_success_step
  character (LEN=30)                              :: report_character
  logical                                         :: lle_flag
  !-----------------------------------------------------------------------------

  kij_bound = 0.4_dp

  penalty = 0.0_dp

  !-----------------------------------------------------------------------------
  ! assign array of optimized parameters back to binary parameters kij, lij, ...
  ! For unreasonable (high or low) values of kij, lij, .. define a penalty function
  !-----------------------------------------------------------------------------
  if ( adjust_parameter_option /= 8 ) then
     kij(1,2) = optpars(1) - 1.0_dp
     if ( abs( kij(1,2) ) > kij_bound ) then
        penalty = abs( kij(1,2) ) - kij_bound
        kij(1,2) = sign( kij_bound, kij(1,2) )
     end if
     kij(2,1) = kij(1,2)
     if ( assoc ) call cross_association_parameters
  end if
  if ( adjust_parameter_option == 2 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
     if ( abs( lij(1,2) ) > 0.1_dp ) then
        penalty = penalty + abs( lij(1,2) ) - 0.1_dp
        lij(1,2) = sign( 0.1_dp, lij(1,2) )
     end if
     lij(2,1) = - lij(1,2)
     if ( assoc ) call cross_association_parameters
  else if ( adjust_parameter_option == 3 ) then
     kij_assoc(1,2) = optpars(2) - 1.0_dp
     if ( abs( kij_assoc(1,2) ) > 0.9_dp ) then
        penalty = penalty + abs( kij_assoc(1,2) ) - 0.9_dp
        kij_assoc(1,2) = sign( 0.9_dp, kij_assoc(1,2) )
     end if
     kij_assoc(2,1) = kij_assoc(1,2)
     call cross_association_parameters
  else if ( adjust_parameter_option == 4 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
     if ( abs( lij(1,2) ) > 0.1_dp ) then
        penalty = penalty + abs( lij(1,2) ) - 0.1_dp
        lij(1,2) = sign( 0.1_dp, lij(1,2) )
     end if
     lij(2,1) = - lij(1,2)
     kij_assoc(1,2) = optpars(3) - 1.0_dp
     if ( abs( kij_assoc(1,2) ) > 0.9_dp ) then
        penalty = penalty + abs( kij_assoc(1,2) ) - 0.9_dp
        kij_assoc(1,2) = sign( 0.9_dp, kij_assoc(1,2) )
     end if
     kij_assoc(2,1) = kij_assoc(1,2)
     call cross_association_parameters
  else if ( adjust_parameter_option == 5 ) then
     lij(1,2) = ( optpars(2) - 1.0_dp )**lij_exponent
     if ( abs( lij(1,2) ) > 0.1_dp ) then
        penalty = penalty + abs( lij(1,2) ) - 0.1_dp
        lij(1,2) = sign( 0.1_dp, lij(1,2) )
     end if
     lij(2,1) = ( optpars(3) - 1.0_dp )**lij_exponent
     if ( abs( lij(2,1) ) > 0.1_dp ) then
        penalty = penalty + abs( lij(2,1) ) - 0.1_dp
        lij(2,1) = sign( 0.1_dp, lij(2,1) )
     end if
     call cross_association_parameters
  else if ( adjust_parameter_option == 6 ) then
     kij_t(1,2) = ( optpars(2) - 1.0_dp )
     kij_t(2,1) = kij_t(1,2)
     call cross_association_parameters
  else if ( adjust_parameter_option == 7 ) then
     kij_t(1,2) = ( optpars(2) - 1.0_dp )
     kij_t(2,1) = kij_t(1,2)
     kij_assoc(1,2) = optpars(3) - 1.0_dp
     if ( abs( kij_assoc(1,2) ) > 0.9_dp ) then
        penalty = penalty + abs( kij_assoc(1,2) ) - 0.9_dp
        kij_assoc(1,2) = sign( 0.9_dp, kij_assoc(1,2) )
     end if
     kij_assoc(2,1) = kij_assoc(1,2)
     call cross_association_parameters
  else if ( adjust_parameter_option == 8 ) then
      kij_assoc(1,2) = optpars(1) - 1.0_dp
      if ( abs( kij_assoc(1,2) ) > 0.9_dp ) then
         penalty = penalty + abs( kij_assoc(1,2) ) - 0.9_dp
         kij_assoc(1,2) = sign( 0.9_dp, kij_assoc(1,2) )
      end if
      kij_assoc(2,1) = kij_assoc(1,2)
      call cross_association_parameters
  end if
  ! write (*,'(a,4(f14.8))') ' binary par entering OF ====', optpars(1:nadjust) - 1.0_dp

  !-----------------------------------------------------------------------------
  ! For unreasonable values of kij, lij.., set objectiv fct = penalty fct and exit
  !-----------------------------------------------------------------------------

  if ( penalty > machine_eps ) then
     fmin = 100000.0_dp + 100000.0_dp * penalty
     write (*,*) 'penalty function because the bound were exceeded', penalty, fmin
     return
  end if

  allocate( dev( n_exp ) )

  !-----------------------------------------------------------------------------
  ! loop for experimental mixture data
  !-----------------------------------------------------------------------------

  do i = 1, n_exp

     report_character = ' '
     dev(i) = 0.0_dp

     !--------------------------------------------------------------------------
     if ( objective_type == 1 ) then
     !--------------------------------------------------------------------------

        weigh_x_p = 0.75_dp   ! weight of errors in x to errors in p
        weigh_x_y = 0.8_dp    ! weight of errors in x to errors in y

        call deviation_x_and_y ( i, dev_x_i, dev_y_i, converg, report_success_step, lle_flag )

        deviation(i,1) = dev_x_i
        deviation(i,2) = dev_y_i

        if ( .NOT. lle_flag ) then

           call deviation_bubble_point_p ( i, p_trial, converg_bbp )
           deviation(i,3) = p_trial / p_exp(i) - 1.0_dp

        else

           deviation(i,3) = 0.0_dp
           converg_bbp = .false.

        end if

        !-----------------------------------------------------------------------
        ! calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------
        if ( converg .AND. converg_bbp ) then

           dev(i) = ( weigh_x_p / abs(dev_x_i) + (1.0_dp -weigh_x_p) * abs( p_trial / p_exp(i) - 1.0_dp )**(-1) )**(-2)
           ! dev(i) = ( weigh_x_p / dev_x_i**2 + (1.0_dp -weigh_x_p) * ( t_trial / t_exp(i) - 1.0_dp )**(-2) )**(-1)
           dev(i) = weigh_x_y * dev(i) + ( 1.0_dp - weigh_x_y ) * dev_y_i**2

        else if ( converg_bbp ) then

           dev(i) = 5.0_dp * ( p_trial/p_exp(i) - 1.0_dp )**2
           report_character = 'only delta(p) converged'

        else if ( converg ) then

           dev(i) = weigh_x_y * dev_x_i**2 + ( 1.0_dp - weigh_x_y ) * dev_y_i**2
           report_character = 'only delta(x,y) converged'

        else

           write (*,'(a,i4,a,2G18.11)') ' pt. no. ',i,' not converged, x,y=', x_exp(i), y_exp(i)

        end if

     !--------------------------------------------------------------------------
     else if ( objective_type == 2 ) then
     !--------------------------------------------------------------------------

        call deviation_bubble_point_p ( i, p_trial, converg_bbp )

        deviation(i,3) = p_trial / p_exp(i) - 1.0_dp

        !-----------------------------------------------------------------------
        ! calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------
        if ( converg_bbp ) then
           dev(i) = 5.0_dp * ( p_trial/p_exp(i) - 1.0_dp )**2
        else
           write (*,'(a,i4,a,2G18.11)') ' pt. no. ',i,' not converged, x,y=', x_exp(i), y_exp(i)
        end if

     !--------------------------------------------------------------------------
     else if ( objective_type == 3 ) then
     !--------------------------------------------------------------------------

        call deviation_L_normal_to_envelope ( i, distanceL, exp_vec, convergL )
        if ( convergL ) then
           distance = distanceL
           ! deviation(i,1) = exp_vec(1)
           ! deviation(i,3) = exp_vec(2)
           else
           distance = 0.0_dp
           write (*,'(a,i4,a,2G18.11)') ' L-pt. no. ',i,' not converged, x,y=', x_exp(i), y_exp(i)
        end if

!!$        call deviation_V_normal_to_envelope ( i, distanceV, convergV )
!!$        if ( convergV ) then
!!$           distance = distance + distanceV
!!$        else
!!$           write (*,'(a,i4,a,2G18.11)') ' V-pt. no. ',i,' not converged, x,y=', x_exp(i), y_exp(i)
!!$           ! distance = deviation(i,3)
!!$           ! dev(i) = distance**2
!!$        end if

        deviation(i,3) = distance

        !-----------------------------------------------------------------------
        ! calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------
        ! if ( convergL .AND. convergV ) then
        if ( convergL ) then
           dev(i) = distance**2
        end if

     !--------------------------------------------------------------------------
     else if ( objective_type == 4 ) then
     !--------------------------------------------------------------------------

        call deviation_fugacity ( i, residual )

        deviation(i,3) = sum ( residual(:) )

        !-----------------------------------------------------------------------
        ! calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------
        dev(i) = sum( ( residual(:) )**2 )

     !--------------------------------------------------------------------------
     end if
     !--------------------------------------------------------------------------

     !--------------------------------------------------------------------------
     ! point-wise: write output, if applicable
     !--------------------------------------------------------------------------
     if ( outp >= 1) write (*,'(a,i3,G20.12,2x,a,x,a)') 'deviation',i,dev(i), report_success_step, report_character

  end do

  !-----------------------------------------------------------------------------
  ! write output
  !-----------------------------------------------------------------------------
  error = sum( dev(1:n_exp) ) / real( n_exp, kind=dp ) * 1000.0_dp
  fmin = error

  if ( adjust_parameter_option == 2 ) then
     write (*,'(a,f14.8,a,2(f14.8))') 'deviation =', error,'     kij =', optpars(1) -1.0_dp, lij(1,2)
  else if ( adjust_parameter_option == 4 ) then
     write (*,'(a,2f14.8,a,3(f14.8))') 'deviation =', error, sum(abs(deviation(:,1))+abs(deviation(:,2))),  &
          '     kij =', optpars(1)-1.0_dp, lij(1,2), optpars(3)-1.0_dp
  else if ( adjust_parameter_option == 5 ) then
     write (*,'(a,2f14.8,a,3(f14.8))') 'deviation =', error, sum(abs(deviation(:,1))+abs(deviation(:,2))),  &
          '     kij =', optpars(1)-1.0_dp, lij(1,2), lij(2,1)
  else
     write (*,'(a,f14.8,a,3(f14.8))') 'deviation =', error,'     kij =', optpars(1:nadjust) - 1.0_dp
  end if

  !-----------------------------------------------------------------------------
  ! ignore point that did not converge or did not converge previously, although
  ! parameters are very close
  !-----------------------------------------------------------------------------
  do i = 1, nadjust
     dev_previous(:,i+1) = dev_previous(:,i)
     par_previous(:,i+1) = par_previous(:,i)
  end do
  dev_previous(:,1) = dev(:)
  par_previous(:,1) = optpars(:)

  delta_par = sum( abs( optpars(1:nadjust) - par_previous(1:nadjust,2) ) )

  do i = 1, n_exp
     delta_dev = ( dev(i) - dev_previous(i,2) )  /  ( sum( dev(1:n_exp) ) / real( n_exp, kind=dp ) )
     if ( delta_dev > 0.2_dp .AND. delta_par < 1.E-4 ) then
        write (*,'(a,i4,2f14.8)') 'correcting unexpected jump in objective fct.',i, delta_dev, delta_par
        dev_previous(i,1) = dev_previous(i,2)
        error = sum( dev_previous(1:n_exp,1) ) / real( n_exp, kind=dp ) * 1000.0_dp
        fmin = error
        write (*,'(a,f14.8,a)') 'deviation =', error,' revised'
     end if
     ! delta_par = optpars(i) - par_previous(i,2)
     ! write (*,*) i,delta_dev
     ! if (  .AND.  ) then
     ! dev_previous(:,i+1) = dev_previous(:,i)
     ! par_previous(:,i+1) = par_previous(:,i)
  end do

  deallocate( dev )

end subroutine kij_objective_fct



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE deviation_L_normal_to_envelope
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deviation_L_normal_to_envelope ( i, distance, exp_vec, converg_tot )

  use bubble_dew_point, only: bubble_point_rachford_rice

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i
  real(dp), intent(out)                           :: distance
  real(dp), dimension(2), intent(out)             :: exp_vec
  logical, intent(out)                            :: converg_tot

  !-----------------------------------------------------------------------------
  integer                                         :: it
  real(dp), dimension(ncomp)                      :: xiF
  real(dp), dimension(ncomp)                      :: xiB
  real(dp), dimension(ncomp)                      :: xi2_trial
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(2)                          :: eta_start
  real(dp), dimension(2)                          :: line_vec
  real(dp), dimension(2)                          :: shift_vec
  real(dp)                                        :: error
  real(dp)                                        :: p_trial, t_trial
  real(dp)                                        :: p_trial2
  real(dp)                                        :: delta_x
  real(dp)                                        :: damping
  real(dp), parameter                             :: delta_x_abs = 1.E-4_dp
  real(dp), parameter                             :: tolerance = 1.E-8_dp
  logical                                         :: converg
  logical                                         :: iterate_t
  !-----------------------------------------------------------------------------

  iterate_t = .false.

  p_trial = p_exp(i)          ! starting value
  t_trial = t_exp(i)
  xiF(2) = x_exp(i)
  xiF(1) = 1.0_dp - xiF(2)
  xi2_trial(:) = xiF(:)

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.E-5_dp

  if ( x_exp(i) < 0.5_dp) then
     delta_x = delta_x_abs
  else
     delta_x = - delta_x_abs
  end if

  shift_vec(1) = 0.0_dp
  distance = 0.0_dp

  !-----------------------------------------------------------------------------
  ! find xi, where distance of equilibrium line to exp point is minial
  !-----------------------------------------------------------------------------
  it = 0

  error = tolerance + 1.0_dp
  do while ( error > tolerance .AND. it < 60 )

     it = it + 1

     damping = 1.0_dp
     if ( it <= 2 ) damping = 0.75_dp
     if ( it > 8 .AND. abs( shift_vec(1) ) < 1.E-5_dp ) damping = 0.5_dp
     if ( it > 25 ) damping = 0.25_dp
     xiF(2) = xiF(2) + shift_vec(1) * damping
     xiF(1) = 1.0_dp - xiF(2)
     if ( xiF(2) > 1.0_dp ) write (*,*) 'extend deviation_L_normal_to_envelope. xiF(2) > 1.0_dp. stop'
     if ( xiF(2) < 0.0_dp ) write (*,*) 'extend deviation_L_normal_to_envelope. xiF(2) < 0.0_dp. stop'
     xi2_trial(:) = xiF(:)
     call bubble_point_rachford_rice ( iterate_t, t_trial, p_trial, eta_start, xiF, xi2_trial,  &
          rhoi1, rhoi2, get_start_val, converg )
     if ( .NOT. converg ) exit

     if ( xiF(2) > 1.0_dp-delta_x ) write (*,*) 'extend deviation_L_normal_to_envelope. stop'
     if ( xiF(2) > 1.0_dp-delta_x ) delta_x = - delta_x
     !if ( xiF(2) > 1.0_dp-delta_x ) stop
     xiB(2) = xiF(2) + delta_x
     xiB(1) = 1.0_dp - xiB(2)
     p_trial2 = p_trial
     call bubble_point_rachford_rice ( iterate_t, t_trial, p_trial2, eta_start, xiB, xi2_trial,  &
          rhoi1, rhoi2, get_start_val, converg )
     if ( .NOT. converg ) exit

     ! if ( .NOT. converg ) call VLE_for_given_tp ( t_trial, p_trial, eta_start, xi1, xi2, rhoi1, rhoi2, converg, report_success_step )

     line_vec(1:2) = (/ delta_x, (p_trial2-p_trial)/p_exp(i) /)
     line_vec(1:2) = line_vec(1:2) / sqrt( dot_product(line_vec(1:2),line_vec(1:2)) )
     exp_vec(1:2) = (/ (x_exp(i)-xiF(2)), (p_exp(i)-p_trial)/p_exp(i) /)
     shift_vec(:) = line_vec(1:2) * dot_product( line_vec, exp_vec )

     distance = sqrt( dot_product( exp_vec, exp_vec ) )

     if ( shift_vec(1) > 1.0_dp - xiF(2) ) then
        shift_vec(1) = 1.0_dp - xiF(2)
        if ( outp > 1 .AND. it > 2 ) write (*,*) 'criterion A: point i no real contribution to OF',i
     end if
     if ( shift_vec(1) < - xiF(2) ) then
        shift_vec(1) = - xiF(2)
        if ( outp > 1 .AND. it > 2 )  write (*,*) 'criterion B: point i no real contribution to OF',i
     end if

     error = abs( shift_vec(1) )
     if ( outp > 1 ) write (*,'(a,i4,3G18.8)') 'length-shift, shift-x2L-dir., distance',i,  &
          sqrt( dot_product( shift_vec, shift_vec ) ), shift_vec(1), distance

  end do

  !-----------------------------------------------------------------------------
  ! converged result for point i ?
  ! (accept errors lower than 100-fold tolerance as solutions)
  !-----------------------------------------------------------------------------
  if ( error > tolerance * 100.0_dp ) then
     if ( error < tolerance * 1.E4_dp ) write (*,*) 'point',i,'converged to ',error, tolerance
     if ( outp > 2 ) read (*,*)
     distance = 0.0_dp
     converg_tot = .false.
  else
     converg_tot = .true.
     ! write (*,'(a,4G18.8)') 'tpx', t_trial, p_trial2, xiB(1), xi2_trial(1)
  end if

end subroutine deviation_L_normal_to_envelope

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE deviation_V_normal_to_envelope
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deviation_V_normal_to_envelope ( i, distance, converg_tot )

  use bubble_dew_point, only: bubble_point_rachford_rice

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i
  real(dp), intent(out)                           :: distance
  logical, intent(out)                            :: converg_tot

  !-----------------------------------------------------------------------------
  integer                                         :: it
  real(dp), dimension(ncomp)                      :: xiF
  real(dp), dimension(ncomp)                      :: xiB
  real(dp), dimension(ncomp)                      :: xi2_trial
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(2)                          :: eta_start
  real(dp), dimension(2)                          :: line_vec
  real(dp), dimension(2)                          :: exp_vec
  real(dp), dimension(2)                          :: shift_vec
  real(dp)                                        :: error
  real(dp)                                        :: p_trial, t_trial
  real(dp)                                        :: p_trial2
  real(dp)                                        :: delta_x
  real(dp)                                        :: damping
  real(dp), parameter                             :: delta_x_abs = 1.E-4_dp
  real(dp), parameter                             :: tolerance = 1.E-8_dp
  logical                                         :: converg
  logical                                         :: iterate_t
  !-----------------------------------------------------------------------------

  iterate_t = .false.

  p_trial = p_exp(i)          ! starting value
  t_trial = t_exp(i)
  xiF(2) = y_exp(i)
  xiF(1) = 1.0_dp - xiF(2)
  xi2_trial(:) = xiF(:)

  eta_start(1) = 1.E-5_dp
  eta_start(2) = 0.4_dp

  if ( y_exp(i) < 0.5_dp) then
     delta_x = delta_x_abs
  else
     delta_x = - delta_x_abs
  end if

  shift_vec(1) = 0.0_dp
  distance = 0.0_dp

  !-----------------------------------------------------------------------------
  ! find xi, where distance of equilibrium line to exp point is minial
  !-----------------------------------------------------------------------------
  it = 0

  error = tolerance + 1.0_dp
  do while ( error > tolerance .AND. it < 60 )

     it = it + 1

     damping = 1.0_dp
     if ( it <= 2 ) damping = 0.25_dp
     if ( it > 8 .AND. abs( shift_vec(1) ) < 1.E-5_dp ) damping = 0.5_dp
     if ( it > 25 ) damping = 0.25_dp
     xiF(2) = xiF(2) + shift_vec(1) * damping
     xiF(1) = 1.0_dp - xiF(2)
     if ( xiF(2) > 1.0_dp ) write (*,*) 'extend deviation_V_normal_to_envelope. xiF(2) > 1.0_dp. stop'
     if ( xiF(2) < 0.0_dp ) write (*,*) 'extend deviation_V_normal_to_envelope. xiF(2) < 0.0_dp. stop'
     ! xi2_trial(:) = xiF(:)
     xi2_trial(2) = x_exp(i)
     xi2_trial(1) = 1.0_dp - xi2_trial(2)
     call bubble_point_rachford_rice ( iterate_t, t_trial, p_trial, eta_start, xiF, xi2_trial,  &
          rhoi1, rhoi2, .false., converg )
     if ( .NOT. converg ) exit

     if ( xiF(2) > 1.0_dp-delta_x ) write (*,*) 'extend deviation_V_normal_to_envelope. stop'
     if ( xiF(2) > 1.0_dp-delta_x ) stop
     xiB(2) = xiF(2) + delta_x
     xiB(1) = 1.0_dp - xiB(2)
     p_trial2 = p_trial
     xi2_trial(2) = x_exp(i)
     xi2_trial(1) = 1.0_dp - xi2_trial(2)
     call bubble_point_rachford_rice ( iterate_t, t_trial, p_trial2, eta_start, xiB, xi2_trial,  &
          rhoi1, rhoi2, .false., converg )
     if ( .NOT. converg ) exit

     ! if ( .NOT. converg ) call VLE_for_given_tp ( t_trial, p_trial, eta_start, xi1, xi2, rhoi1, rhoi2, converg, report_success_step )

     line_vec(1:2) = (/ delta_x, (p_trial2-p_trial)/p_exp(i) /)
     line_vec(1:2) = line_vec(1:2) / sqrt( dot_product(line_vec(1:2),line_vec(1:2)) )
     exp_vec(1:2) = (/ (y_exp(i)-xiF(2)), (p_exp(i)-p_trial)/p_exp(i) /)
     shift_vec(:) = line_vec(1:2) * dot_product( line_vec, exp_vec )

     distance = sqrt( dot_product( exp_vec, exp_vec ) )

     if ( shift_vec(1) > 1.0_dp - xiF(2) ) then
        shift_vec(1) = 1.0_dp - xiF(2)
        if ( outp > 1  .AND. it > 2 ) write (*,*) 'criterion A: point i no real contribution to OF',i
     end if        
     if ( shift_vec(1) < - xiF(2) ) then
        shift_vec(1) = - xiF(2)
        if ( outp > 1 .AND. it > 2 )  write (*,*) 'criterion B: point i no real contribution to OF',i
     end if

     error = abs( shift_vec(1) )
     if ( outp > 1 ) write (*,'(a,i4,3G18.8)') 'length-shift, shift-x2V-dir., distance',i,  &
          sqrt( dot_product( shift_vec, shift_vec ) ), shift_vec(1), distance

  end do

  !-----------------------------------------------------------------------------
  ! converged result for point i ?
  ! (accept errors lower than 1.E4-fold tolerance as solutions)
  !-----------------------------------------------------------------------------
  if ( error > tolerance * 1.E4_dp ) then
     if ( error < tolerance * 1.E5_dp ) write (*,*) 'point',i,'converged to ',error, tolerance
     if ( outp > 2 ) read (*,*)
     distance = 0.0_dp
     converg_tot = .false.
  else
     converg_tot = .true.
  end if

end subroutine deviation_V_normal_to_envelope


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE deviation_bubble_point_p
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deviation_bubble_point_p ( i, p_trial, converg )

  use bubble_dew_point, only: bubble_point_rachford_rice

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i
  real(dp), intent(out)                           :: p_trial
  logical, intent(out)                            :: converg

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: xiF
  real(dp), dimension(ncomp)                      :: xi2_trial
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(2)                          :: eta_start
  real(dp)                                        :: t_trial
  logical                                         :: iterate_t
  !-----------------------------------------------------------------------------

  iterate_t = .false.

  p_trial = p_exp(i)          ! starting value
  t_trial = t_exp(i)
  xiF(2) = x_exp(i)
  xiF(1) = 1.0_dp - xiF(2)
  xi2_trial(:) = xiF(:)

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.E-5_dp

  call bubble_point_rachford_rice ( iterate_t, t_trial, p_trial, eta_start, xiF, xi2_trial,  &
       rhoi1, rhoi2, get_start_val, converg )

  if ( .not. converg ) p_trial = p_exp(i)

end subroutine deviation_bubble_point_p


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE deviation_x_and_y
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deviation_x_and_y ( i, dev_x, dev_y, converg, report_success_step, lle_flag )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i
  real(dp), intent(out)                           :: dev_x
  real(dp), intent(out)                           :: dev_y
  logical, intent(out)                            :: converg
  character (LEN=1), intent(out)                  :: report_success_step
  logical, intent(out)                            :: lle_flag

  !-----------------------------------------------------------------------------
  real(dp)                                        :: t, p
  real(dp), dimension(ncomp)                      :: xi1, xi2
  real(dp), dimension(ncomp)                      :: xi1_calc, xi2_calc
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(2)                          :: eta_start
  real(dp)                                        :: x_sta, y_sta
  real(dp)                                        :: xi_save
  real(dp)                                        :: eta1, eta2
  !-----------------------------------------------------------------------------

  converg = .false.
  lle_flag = .false.

  dev_x = 0.0_dp
  dev_y = 0.0_dp

  report_success_step = ' '

  t = t_exp(i)
  p = p_exp(i)

  !--------------------------------------------------------------------------
  ! starting values for composition x and density eta
  !--------------------------------------------------------------------------

  x_sta = x_exp(i)
  y_sta = y_exp(i)
  if ( abs(x_sta) < machine_eps ) x_sta = 0.001_dp          ! these are only starting values
  if ( abs(x_sta-1.0_dp) < machine_eps ) x_sta = 0.999_dp   ! these are only starting values
  if ( abs(y_sta) < machine_eps ) y_sta = 0.001_dp          ! these are only starting values
  if ( abs(y_sta-1.0_dp) < machine_eps ) y_sta = 0.999_dp   ! these are only starting values

  xi1(2) = x_sta
  xi2(2) = y_sta
  xi1(1) = 1.0_dp - x_sta
  xi2(1) = 1.0_dp - y_sta

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.E-5_dp

  call VLE_for_given_tp ( t, p, eta_start, xi1, xi2, rhoi1, rhoi2, converg, report_success_step )

  if ( converg ) then

     call occupy_target ( rhoi1, rhoi2, xi1_calc, xi2_calc, eta1, eta2 )

     if ( abs( eta1 / eta2 - 1.0_dp ) < 0.25_dp ) then

        lle_flag = .true.
        !-----------------------------------------------------------------------
        ! converged to LLE, first check, whether the phase index has to be exchanged
        !-----------------------------------------------------------------------
        if ( abs(x_exp(i)) > machine_eps .and. abs(x_exp(i)-1.0_dp) > machine_eps .and.  &
             abs(xi1_calc(2)-x_exp(i)) > abs(xi2_calc(2)-x_exp(i)) ) then
           write (*,'(a,i3,6G20.12)') 'swap1', i, xi1_calc(2), x_exp(i), xi2_calc(2), y_exp(i), t_exp(i), p_exp(i)
           xi_save = xi1_calc(2)
           xi1_calc(2) = xi2_calc(2)
           xi2_calc(2) = xi_save
        end if
        if ( abs(y_exp(i)) > machine_eps .and. abs(y_exp(i)-1.0_dp) > machine_eps .and.  &
             abs(xi2_calc(2)-y_exp(i)) > abs(xi1_calc(2)-y_exp(i)) ) then
           write (*,'(a,i3,6G20.12)') 'swap2', i, xi1_calc(2), x_exp(i), xi2_calc(2), y_exp(i), t_exp(i), p_exp(i)
           xi_save = xi1_calc(2)
           xi1_calc(2) = xi2_calc(2)
           xi2_calc(2) = xi_save
        end if

     end if

     !--------------------------------------------------------------------------
     ! converged results, calculate the x-deviations and y-deviations
     !--------------------------------------------------------------------------
     if ( abs(x_exp(i)) > machine_eps .and. abs(x_exp(i)-1.0_dp) > machine_eps ) dev_x = xi1_calc(2) - x_exp(i)
     if ( abs(y_exp(i)) > machine_eps .and. abs(y_exp(i)-1.0_dp) > machine_eps ) dev_y = xi2_calc(2) - y_exp(i)

  else

     write (*,*) 'phase equilibrium for defined T,p not converged, point', i

  end if

end subroutine deviation_x_and_y



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE deviation_fugacity
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine deviation_fugacity ( i, residual )

  use properties, only: chemical_potential_tp

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: i
  real(dp), intent(out), dimension(ncomp)         :: residual

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: xiL
  real(dp), dimension(ncomp)                      :: xiV
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(ncomp)                      :: lnphi_1, mu1
  real(dp), dimension(ncomp)                      :: lnphi_2, mu2
  real(dp), dimension(2)                          :: eta_start
  real(dp)                                        :: tF, pF
  !-----------------------------------------------------------------------------

  pF = p_exp(i)
  tF = t_exp(i)
  xiL(2) = x_exp(i)
  xiL(1) = 1.0_dp - xiL(2)
  xiV(2) = y_exp(i)
  xiV(1) = 1.0_dp - xiV(2)

  eta_start(1) = 0.4_dp
  eta_start(2) = 1.E-5_dp

  call chemical_potential_tp ( tF, pF, xiL, eta_start(1), rhoi1, lnphi_1, mu1 )
  call chemical_potential_tp ( tF, pF, xiV, eta_start(2), rhoi2, lnphi_2, mu2 )

  residual(1:ncomp) = log( xiL(1:ncomp) ) + lnphi_1(1:ncomp) - log( xiV(1:ncomp) ) - lnphi_2(1:ncomp)
  ! residual(1:ncomp) = xiL(1:ncomp) * exp( lnphi_1(1:ncomp) ) - ( xiV(1:ncomp) * exp( lnphi_2(1:ncomp) ) )
  ! residual(1:ncomp) = xiL(1:ncomp) * exp( lnphi_1(1:ncomp) ) / ( xiV(1:ncomp) * exp( lnphi_2(1:ncomp) ) ) - 1.0_dp

end subroutine deviation_fugacity


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE VLE_for_given_tp
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine VLE_for_given_tp ( tF, pF, eta_start, xi1, xi2, rhoi1, rhoi2, converg, report_success_step )

  use flash, only: flash_driver, rachford_rice
  use module_iterate_pe, only: iterate_phase_equilibrium

  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: tF
  real(dp), intent(in)                            :: pF
  real(dp), dimension(2), intent(in)              :: eta_start
  real(dp), dimension(ncomp), intent(inout)       :: xi1, xi2
  real(dp), dimension(ncomp), intent(out)         :: rhoi1, rhoi2
  logical, intent(out)                            :: converg
  character (LEN=1), intent(out)                  :: report_success_step

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp)                      :: xiF
  real(dp), dimension(ncomp)                      :: xi
  real(dp), dimension(ncomp)                      :: rhoi3
  real(dp)                                        :: phase_frac2
  logical                                         :: slow_iteration
  logical                                         :: likely_phase_split
  logical                                         :: converg_3
  !-----------------------------------------------------------------------------

  converg = .false.

  report_success_step = ' '

  !--------------------------------------------------------------------------
  ! step A for solving the pT-flash problem: Rachford-Rice
  !--------------------------------------------------------------------------

  xiF(1) = 0.5_dp * ( xi1(1) + xi2(1) )
  xiF(2) = 1.0_dp - xiF(1)

  call rachford_rice ( tF, pF, xiF, rhoi1, rhoi2, phase_frac2, slow_iteration, converg, x1s=xi1, x2s=xi2 )
  if ( converg ) report_success_step = 'A'


  !--------------------------------------------------------------------------
  ! step B for solving the pT-flash problem: iterating compositions of both phases directly
  !--------------------------------------------------------------------------

  if ( .NOT.converg ) then

     call iterate_phase_equilibrium ( tF, pF, eta_start, xi1, xi2, rhoi1, rhoi2, converg )
     !if ( converg ) eta1 = get_eta( rhoi1 )
     !if ( converg ) eta2 = get_eta( rhoi2 )
     if ( converg ) report_success_step = 'B'

  end if

  !--------------------------------------------------------------------------
  ! step C for solving the pT-flash problem: VLE-scan
  !--------------------------------------------------------------------------

  if ( .NOT.converg ) then

     call binary_x_scan_for_vle_or_LLE ( tF, pF, likely_phase_split, xi )
     if ( likely_phase_split ) then
        call flash_driver ( tF, pF, xi, rhoi1, rhoi2, rhoi3, converg, converg_3 )
        if ( converg ) report_success_step = 'C'
        if ( converg_3 ) write (*,*) '3-phase equilibrium found'
     end if

  end if

  !--------------------------------------------------------------------------
  ! step D for solving the pT-flash problem: x-scan and flash
  !--------------------------------------------------------------------------

!!$     if ( .NOT.converg ) then
!!$
!!$        outp = 0                    ! output to terminal
!!$        call scan_compositions ( converg )
!!$        if ( converg == 1 ) then
!!$           val_conv = val_init
!!$           call restore_converged
!!$           report_success_step = 'D'
!!$        end if
!!$
!!$     end if
  ! write (*,*) report_success_step
  ! read (*,*)

end subroutine VLE_for_given_tp


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE occupy_target
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine occupy_target ( rhoi1, rhoi2, xi1_calc, xi2_calc, eta1, eta2 )

  use stability, only: get_eta

  !-----------------------------------------------------------------------------
  real(dp), dimension(ncomp), intent(in)          :: rhoi1
  real(dp), dimension(ncomp), intent(in)          :: rhoi2
  real(dp), dimension(ncomp), intent(out)         :: xi1_calc
  real(dp), dimension(ncomp), intent(out)         :: xi2_calc
  real(dp), intent(out)                           :: eta1
  real(dp), intent(out)                           :: eta2
  !-----------------------------------------------------------------------------
  
  xi1_calc(:) = rhoi1(:) / sum( rhoi1(1:ncomp) )
  xi2_calc(:) = rhoi2(:) / sum( rhoi2(1:ncomp) )
  eta1 = get_eta( rhoi1 )
  eta2 = get_eta( rhoi2 )

end subroutine occupy_target


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE binary_x_scan_for_vle_or_LLE
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine binary_x_scan_for_vle_or_LLE ( t, p, likely_phase_split, xi )

  use stability, only: get_eta
  use properties, only: chemical_potential_tp
  !-----------------------------------------------------------------------------
  real(dp), intent(in)                            :: t
  real(dp), intent(in)                            :: p
  logical, intent(out)                            :: likely_phase_split
  real(dp), dimension(ncomp), intent(out)         :: xi

  !-----------------------------------------------------------------------------
  integer, parameter                              :: steps = 120

  integer                                         :: i, j, k
  integer, dimension(0:steps)                     :: phases
  real(dp), dimension(ncomp)                      :: lnphi_1, lnphi_2
  real(dp), dimension(ncomp)                      :: mu1, mu2
  real(dp), dimension(ncomp)                      :: rhoi1, rhoi2
  real(dp), dimension(2)                          :: eta_start
  real(dp)                                        :: rho1, rho2
  real(dp)                                        :: gibbs1, gibbs2
  real(dp), dimension(0:steps)                    :: vlemin, llemin, xval, start_xv, start_xl
  real(dp)                                        :: dg_dx2
  logical                                         :: lle_check  ! currently not any more used
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! initialize quantities
  !-----------------------------------------------------------------------------
  j = 0
  k = 0

  start_xv(:) = 0.0_dp
  start_xl(:) = 0.0_dp

  eta_start(1) = 0.45_dp
  eta_start(2) = 1.E-6_dp

  !-----------------------------------------------------------------------------
  ! scan composition x in equidistant grid of 'steps' number of steps
  !-----------------------------------------------------------------------------
  do i = 0, steps

     xi(1) = 1.0_dp - REAL( i, KIND=dp ) / REAL( steps, KIND=dp )
     if ( xi(1) <= 1.E-50_dp ) xi(1) = 1.E-50_dp
     xi(2) = 1.0_dp - xi(1)

     call chemical_potential_tp ( t, p, xi, eta_start(1), rhoi1, lnphi_1, mu1 )
     call chemical_potential_tp ( t, p, xi, eta_start(2), rhoi2, lnphi_2, mu2 )
     rho1 = sum( rhoi1 )
     rho2 = sum( rhoi2 )
     gibbs1 = sum( xi(1:ncomp) * mu1(1:ncomp) )
     gibbs2 = sum( xi(1:ncomp) * mu2(1:ncomp) )

     xval(i) = xi(1)
     llemin(i) = gibbs1

     !--------------------------------------------------------------------------
     ! can a low-density and a high-density be found, for given (T,p,x)?
     !--------------------------------------------------------------------------
     if ( abs( 1.0_dp - rho1 / rho2 ) > 0.0001_dp ) then
        vlemin(i) = gibbs1 - gibbs2               ! divided by kT
        phases(i) = 2
     else
        phases(i) = 1
     end if

     !--------------------------------------------------------------------------
     ! does (gibbs1-gibbs2) change sign? That's a composition where phase split occurs.
     !--------------------------------------------------------------------------
     if ( i > 0 .AND. phases(i) == 2 ) then

        if ( phases(i-1) == 2 .AND. ABS( vlemin(i) + vlemin(i-1) ) <  &
                                    ABS( vlemin(i) ) + ABS( vlemin(i-1) ) ) then
           j = j + 1
           start_xv(j) = xval(i-1) - ( xval(i) - xval(i-1) ) / ( vlemin(i) - vlemin(i-1) ) * vlemin(i-1)
        end if

     end if

  end do

  !-----------------------------------------------------------------------------
  ! look for indication of a LLE
  !-----------------------------------------------------------------------------
  do i = 1, steps - 1

     dg_dx2 = ( llemin(i-1) - 2.0_dp * llemin(i) + llemin(i+1) ) / ( xval(i) - xval(i-1) )**2
     if ( dg_dx2 < 0.0_dp ) then
        k = k + 1
        start_xl(k) = xval(i)
     end if

  end do

  !-----------------------------------------------------------------------------
  ! define composition vector xi, where phase split is likely. Otherwise set
  ! likely_phase_split = .false.
  !-----------------------------------------------------------------------------
  likely_phase_split = .true.

  if ( abs(start_xl(1)) < machine_eps .AND. abs(start_xv(1)) > machine_eps ) then
     xi(1) = start_xv(1)
     xi(2) = 1.0_dp - xi(1)
     lle_check = .false.
     if ( outp >= 2 ) write (*,*) ' VLE is likely', xi(1), xi(2)
  else if ( abs(start_xl(1)) > machine_eps .AND. abs(start_xv(1)) < machine_eps ) then
     xi(1) = start_xl(1)
     xi(2) = 1.0_dp - xi(1)
     if ( outp >= 2 ) write (*,*) ' LLE is likely', xi(1), xi(2)
     lle_check = .true.
  else if ( abs(start_xl(1)) > machine_eps .AND. abs(start_xv(1)) > machine_eps ) then
     xi(1) = start_xv(1)
     xi(2) = 1.0_dp - xi(1)
     if ( outp >= 2 ) write(*,*) ' starting with VLE and check for LLE'
     lle_check = .true.
  else
     if ( outp >= 1 ) write (*,*) ' no promising starting values found in SR binary_x_scan_for_vle_or_LLE'
     likely_phase_split = .false.
  end if

end subroutine binary_x_scan_for_vle_or_LLE

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine Newton_Opt_2D ( f_val, x, n, min_grad, eps, g, f )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: n
  real(dp), dimension(n), intent(in out)          :: x
  real(dp), intent(in out)                        :: f
  real(dp), dimension(n), intent(in out)          :: g
  real(dp), intent(in)                            :: min_grad
  real(dp), intent(in)                            :: eps
  !-----------------------------------------------------------------------------
  interface
     subroutine f_val (n, x, f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine f_val
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: i, count, max_c
  integer                                         :: line_steps
  real(dp)                                        :: ftrial, det
  real(dp)                                        :: error_in_x
  real(dp), dimension(n)                          :: x0
  real(dp), dimension(n,n)                        :: H, invers_H
  real(dp)                                        :: ft, ft_min, length_scale
  real(dp), dimension(n)                          :: xt, direction
  logical                                         :: line_search
  !-----------------------------------------------------------------------------

  if ( n > 2 ) write (*,*) 'Newton_Opt_2D is so far only defined for n<=2'
  if ( n > 2 ) stop

  max_c = 50                    ! maximum number of iterations
  count = 0
  error_in_x = 10.0_dp * eps
  line_steps = 5

  do while ( error_in_x > eps .and. count <= max_c )

     count = count + 1

     x0 = x

     call f_grad_hessian ( f_val, n, x, f, g, H )

     if ( n == 1 ) then
        det = H(1,1)
        invers_H(1,1) = 1.0_dp / H(1,1)
        ! write (*,*) 'g',g
        ! write (*,*) 'H',H
     else if ( n == 2 ) then
        det = H(1,1) * H(2,2) - H(1,2) * H(2,1)
        if ( sum( abs( g(1:n) ) ) < min_grad .or. abs(det) < 1.E-20_dp ) exit

        invers_H(1,1) =   H(2,2)  / det
        invers_H(2,2) =   H(1,1)  / det
        invers_H(1,2) = - H(2,1)  / det
        invers_H(2,1) = - H(1,2)  / det
     else
        ! call jacobi_eigenvalue ( n, H, it_max, eigenvec, eigenval, it_num, rot_num )
        ! if ( outp >= 3 )  write (*,*) 'lowest eigenvalue',eigenval(1)

        ! if ( eigenval(1) < 0.0001_dp .AND. eta_H < 20.0_dp ) then
        !-----------------------------------------------------------------------
        ! solving  AX = B  with Gauss-Jordan method.
        ! Matrix A='Hessian' and vector B='delta_x' are destroyed and
        ! redefined to the inverse of A and solution vector X, respectively.
        !-----------------------------------------------------------------------
        ! delta_x(:) = - g(:)
        ! call MATINV ( n, 1, H, delta_x, DET )
     end if

     if ( det > 0.0_dp ) then
        do i = 1, n
           x(i) = x0(i) - sum( invers_H(i,1:n) * g(1:n) )
        end do
        line_search = .false.
        if ( outp >= 2 ) write (*,*) 'making regular Newton step x_new(1:2), x_old(1:2)', x, x0
     else
        ! if ( outp >= 1 )
        write (*,*) 'Newton_Opt_2D: problem concave, line search',f
        ft_min = f
        length_scale = 0.1_dp * dot_product( x0, x0 )
        direction(1:n) = g(1:n) / dot_product( g, g )
        xt(1:n) = x0(1:n) - length_scale * direction(1:n)
        do i = 1, 2 * line_steps
           xt(1:n) = xt(1:n) + length_scale * direction(1:n) / real( line_steps, KIND=dp )
           call f_val ( n, xt, ft )
           if ( ft < ft_min ) then
              x(1:n) = xt(1:n)
              ft_min = ft
           end if
        end do
        if ( outp >= 1 ) write (*,*) 'finished line search', x(1:n)
        line_search = .true.
     end if

     call f_val ( n, x, ftrial )

     if ( ftrial > f * ( 1.0_dp + 1.E-5_dp ) .and. .not. line_search ) then
        ! if ( outp >= 1 )
        write (*,*) 'reduced step size'
        do i = 1, n
           x(i) = x0(i) - 0.5_dp * sum( invers_H(i,1:n) * g(1:n) )
        end do
        ! call f_val ( n, x, ftrial )
        ! if ( ftrial > f ) then
        !   do i = 1,n
        !      x(i) = x0(i) - 0.25 * SUM( invers_H(i,1:n)*g(1:n) )
        !   end do
        ! end if
     end if

     error_in_x = sum( abs( x(1:n) - x0(1:n) ) )

     if (outp > 1) write (*,'(a,4(G15.7))') ' finished Newton step: f,err,par ',f,sum( abs( g(1:n) ) ), x(1:n)

     ! write (*,*) 'f= ',fmin
     ! write (*,*) 'par',exp( x(1) ), exp( x(2) )
     ! write (*,'(a,4(G15.7))') ' I= ',invers_H(1,1),invers_H(1,2),invers_H(2,1),invers_H(2,2)
     ! write (*,'(a,4(G15.7))') ' H= ',hessian(1,1),hessian(1,2),hessian(2,1),hessian(2,2)

  end do


end subroutine Newton_Opt_2D



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_grad_hessian ( f_val, n, optpara, fmin, gtrans, hessian )

  integer, intent(in)                             :: n
  real(dp), dimension(n), intent(in)              :: optpara
  real(dp), intent(in out)                        :: fmin
  real(dp), dimension(n), intent(in out)          :: gtrans
  real(dp), dimension(n,n), intent(in out)        :: hessian
  !-----------------------------------------------------------------------------
  interface
     subroutine f_val (n, x, f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine f_val
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: i, j
  real(dp)                                        :: delta
  real(dp), dimension(n)                          :: optpara_mod
  real(dp), dimension(n)                          :: g
  real(dp), dimension(n,n)                        :: gi
  real(dp), dimension(n,n)                        :: gi_left
  !-----------------------------------------------------------------------------

  delta = 1.0E-4_dp

  optpara_mod = optpara

  do i = 1, n

     optpara_mod = optpara
     optpara_mod(i) = optpara(i) * ( 1.0_dp + delta )

     call f_grad ( f_val, n, optpara_mod, g )
     gi(i,1:n) = g(1:n)

     optpara_mod(i) = optpara(i) * ( 1.0_dp - delta )

     call f_grad ( f_val, n, optpara_mod, g )
     gi_left(i,1:n) = g(1:n)

  enddo

  call f_grad ( f_val, n, optpara, g )


  do i = 1, n
     do j = 1, n

        hessian(i,j) = ( gi(i,j) - gi_left(i,j) ) / ( 2.0_dp * optpara(i) * delta )
        ! hessian(j,i) = hessian(i,j)
        ! write (*,*) 'hessian',i,j,hessian(i,j)

     enddo
  enddo

  gtrans = g
  call f_val ( n, optpara, fmin )

end subroutine f_grad_hessian


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE f_grad 
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine f_grad ( f_val, n, optpara, g )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: n
  real(dp), dimension(n), intent(in)              :: optpara
  real(dp), dimension(n), intent(in out)          :: g
  !-----------------------------------------------------------------------------
  interface
     subroutine f_val (n, x, f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine f_val
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: delta, fmin
  real(dp)                                        :: f
  real(dp), dimension(n)                          :: optpara_mod, fi
  !-----------------------------------------------------------------------------


  delta = 1.E-5_dp

  optpara_mod = optpara

  do i = 1, n

     optpara_mod = optpara
     optpara_mod(i) = optpara(i) * ( 1.0_dp + delta )

     call f_val ( n, optpara_mod, fmin )
     fi(i) = fmin

  enddo

  call f_val ( n, optpara, fmin )
  f = fmin


  do i = 1, n

     g(i) = ( fi(i) - f ) / ( optpara(i) * delta )
     ! write (*,*) 'g difference scheme', g
     ! write (*,*) fi(i), f, optpara(i)*delta

  enddo

end subroutine f_grad


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!                SIMPLE DRIVER FOR L-BFGS-B (version 2.3)
!     --------------------------------------------------------------

!        L-BFGS-B is a code for solving large nonlinear optimization
!        problems with simple bounds on the variables.

!        The code can also be used for unconstrained problems and is
!        as efficient for these problems as the earlier limited memory
!        code L-BFGS.

!        This is the simplest driver in the package. It uses all the
!        default settings of the code.

!     References:

!        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!        memory algorithm for bound constrained optimization'',
!        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

!        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!        Subroutines for Large Scale Bound Constrained Optimization''
!        Tech. Report, NAM-11, EECS Department, Northwestern University, 1994.

!        (Postscript files of these papers are available via anonymous
!        ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

!                              *  *  *

!        NEOS, November 1994. (Latest revision April 1997.)
!        Optimization Technology Center.
!        Argonne National Laboratory and Northwestern University.
!        Written by
!                           Ciyou Zhu
!        in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

!     NOTE: The user should adapt the subroutine 'timer' if 'etime' is
!           not available on the system.  An example for system
!           AIX Version 3.2 is available at the end of this driver.
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine bfgs_optimize_binary_par ( n, x, f )

  use L_BFGS_B

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: n
  real(dp), dimension(n), intent(in out)          :: x
  real(dp), intent(in out)                        :: f
  ! real(dp), dimension(n), intent(out)             :: g

  !-----------------------------------------------------------------------------
  ! integer, parameter :: nmax  = 30           ! dimension of the largest problem to be solved.

  !     Declare the variables needed by the code.
  !     A description of all these variables is given at the end of the driver.

  character (LEN=60) :: task
  integer            :: m, iprint, nbd(n)
  real(dp)           :: factr, pgtol, l(n), u(n), g(n)

  !     Declare a few additional variables for this sample problem.

  ! real (dp) :: t1, t2        JG: commented line
  integer   :: i
  !-----------------------------------------------------------------------------

  !     We wish to have output at every iteration.

  iprint = 1

  !     We specify the tolerances in the stopping criteria.

  factr  = 1.0E+7_dp
  pgtol  = 1.0E-5_dp

  !     We specify the dimension n of the sample problem and the number
  !     m of limited memory corrections stored.  (n and m should not
  !     exceed the limits nmax.)

  ! n = 25
  m =  5

  !     We now provide nbd which defines the bounds on the variables:
  !     l   specifies the lower bounds,
  !     u   specifies the upper bounds.

  !     First set bounds on the odd-numbered variables.

  do i = 1, n
     nbd(i) = 2
  end do
  l(1)   = 0.6_dp
  u(1)   = 1.5_dp

  l(2)   = 0.6_dp
  u(2)   = 1.4_dp

  l(3)   = 0.5_dp
  u(3)   = 1.5_dp

  !     We now define the starting point.

  ! x(1:n) = 3.0_dp

  !     We now write the heading of the output.

  write (6,16)
16 format(/ '     Solving sample problem.'  &
       / '      (f = 0.0 at the optimal solution.)' /)

  !     We start the iteration by initializing task.

  task = 'START'

  !     ------- The beginning of the loop ----------

  !     This is the call to the L-BFGS-B code.

111 call mainlb (n, m, x, l, u, nbd, f, g, factr, pgtol, task, iprint)

  if (task(1:2) == 'FG') then

     !        The minimization routine has returned to request the
     !        function f and gradient g values at the current x.

     !        Compute function value f for the sample problem.

     call kij_objective_fct ( n, x, f )
     !f = .25_dp*(x(1) - 1._dp)**2
     !do i = 2, n
     !   f = f + (x(i) - x(i-1)**2)**2
     !end do
     !f = 4._dp*f

     !        Compute gradient g for the sample problem.

     call f_grad ( kij_objective_fct, n, x, g )
     !t1   = x(2) - x(1)**2
     !g(1) = 2._dp*(x(1) - 1._dp) - 16._dp*x(1)*t1
     !do i = 2, n-1
     !   t2   = t1
     !   t1   = x(i+1) - x(i)**2
     !   g(i) = 8._dp*t2 - 16._dp*x(i)*t1
     !end do
     !g(n) = 8._dp*t1

     !        Go back to the minimization routine.
     GO TO 111

  else if (task(1:5) == 'NEW_X') then

     !        The minimization routine has returned with a new iterate,
     !        and we have opted to continue the iteration.

     GO TO 111

  else

     !        We terminate execution when task is neither FG nor NEW_X.
     !        We print the information contained in the string task
     !        if the default output is not used and the execution is
     !        not stopped intentionally by the user.

     if (iprint <= -1 .and. task(1:4) /= 'STOP') write(6,*) task

  end if

end subroutine bfgs_optimize_binary_par


end module kij_fitting
