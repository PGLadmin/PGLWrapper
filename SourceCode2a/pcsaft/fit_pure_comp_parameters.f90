!> \file pure_par_fit.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE pure_fit_parameters
!> \brief module-variables for pure component parameter regression
!!
!! Module contains parameters and variables needed for pure component parameter regression.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module pure_fit_parameters

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none
  save

  !-----------------------------------------------------------------------------
  integer, parameter                              :: n_types_of_data = 3  ! number considered properties: rho(T,p), p^sat(T), B(T)
  integer                                         :: n_dense              ! number of selected density data points, rho(T,p)
  integer                                         :: n_lv                 ! number of selected vap. press data points, p^sat(T)
  integer                                         :: n_virial             ! number of selected virial coeff. values, B(T)
  real(dp), dimension(n_types_of_data)            :: type_weight          ! weight for considered properties within objective function

  real(dp)                                        :: mm_pure              ! molecular mass of considered species

  real(dp), dimension(:), allocatable             :: v_calc
  real(dp), dimension(:), allocatable             :: plv_calc
  real(dp), dimension(:), allocatable             :: B_calc

  real(dp), dimension(:), allocatable             :: plv, tlv
  real(dp), dimension(:), allocatable             :: pliq, tliq, vliq
  real(dp), dimension(:), allocatable             :: b_vir, t_vir
  real(dp)                                        :: v_crit

  integer, dimension(:), allocatable              :: i_par                ! index-vector pointing at adjusted pure comp. parameters
  character (LEN=25), dimension(:), allocatable   :: aa_txt

  integer                                         :: outp = 0

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fitting_pure_component_parameters_driver
!> \brief driver for fitting of pure component parameters
!
!! Driver for the fitting of pure component parameters.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fitting_pure_component_parameters_driver

  character (LEN=40)                              :: pure_fit_file
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! read input-file
  !-----------------------------------------------------------------------------

  write (*,*) '  SPECIFY INPUT FILE: ( ./input_file/pure_comp/ )'
  read (*,*) pure_fit_file
  ! pure_fit_file = 'nh3.dat'             ! JG
  ! pure_fit_file = 'acetic_acid.dat'     ! JG

  !-----------------------------------------------------------------------------
  ! optimize pure component parameters
  !-----------------------------------------------------------------------------

  call fitting_pure_component_parameters ( pure_fit_file )

end subroutine fitting_pure_component_parameters_driver
  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fitting_pure_component_parameters
!
!! SR for the fitting pure component parameters. This subroutine
!! uses LMDIF1 (Levenberg-Marquardt MINPACK routine) in order to minimze
!! the deviation of calculation results to experimental values. The
!! deviations are calculated in pure_component_objective_fct and stored in 'fvec'.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fitting_pure_component_parameters ( pure_fit_file )

  use Levenberg_Marquardt
  use module_eos_derivatives, only: allocate_eos_quantities, initialize_eos_quantities,  &
       deallocate_eos_quantities

  !-----------------------------------------------------------------------------
  character (LEN=40), intent(in)                  :: pure_fit_file

  !-----------------------------------------------------------------------------
  integer                                         :: info
  integer                                         :: lm_n_par               ! number of adjustable parameters
  integer                                         :: lm_m_dat               ! number of data points
  integer, allocatable                            :: ipvt(:)
  real(dp), allocatable                           :: para(:)                ! degrees of freedom. Adjusted parameters
  real(dp), allocatable                           :: fvec(:)                ! vector of residuals (to be squared and inimized)
  logical                                         :: flag_exit

  ! --- numerical constants for Levenberg-Marquardt algorithm ------------------
  real(dp), parameter                             :: tol = 1.0E-8_dp
  real(dp)                                        :: epsfcn
  real(dp), parameter                             :: lm_factor = 0.01_dp     ! maximum initial step size (parameter iteration)
  !-----------------------------------------------------------------------------

  epsfcn = 1.0E-6_dp**2  ! sqrt of relat. step size (finite differences)
  flag_exit = .false.

  !-----------------------------------------------------------------------------
  ! forget what was read from INPUT-file, newly allocate & initialize array quantities
  !-----------------------------------------------------------------------------
  call deallocate_eos_quantities
  ncomp = 1
  call allocate_eos_quantities
  call initialize_eos_quantities

  !-----------------------------------------------------------------------------
  ! open output-file
  !-----------------------------------------------------------------------------
  open (23,FILE='./output_file/pure_comp/'//pure_fit_file)

  !-----------------------------------------------------------------------------
  ! read experimental data
  !-----------------------------------------------------------------------------
  call read_pure_component_data ( pure_fit_file, flag_exit, lm_m_dat )
  if ( flag_exit ) return

  !-----------------------------------------------------------------------------
  ! read EOS-parameters
  !-----------------------------------------------------------------------------
  call define_adjusted_pure_parameters ( lm_n_par, para )

  allocate( fvec(lm_m_dat) )
  allocate( ipvt(lm_n_par) )


  !-----------------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  !-----------------------------------------------------------------------------

  call lmdif1 (pure_component_objective_fct,lm_m_dat,lm_n_par,para,fvec,tol,epsfcn,lm_factor,info,ipvt)


  call output_pure_component_parameters ( pure_fit_file, info, lm_n_par, para )

  deallocate( para, i_par, aa_txt )
  deallocate( fvec, ipvt )
  if ( n_dense > 0 ) deallocate( v_calc, pliq, tliq, vliq )
  if ( n_lv > 0 ) deallocate( plv_calc, plv, tlv )
  if ( n_virial > 0 ) deallocate( B_calc, b_vir, t_vir )

  close (23)

end subroutine fitting_pure_component_parameters



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pure_component_objective_fct
!> \brief objective function for pure component parameter fitting
!!
!! Routine for pure component parameter fitting. The residual between
!! calculated values and experimental data is calculated and stored in
!! vector 'fvec'. The vector 'fvec' serves as the objective function for
!! the parameter.
!! PURE_RESIDUAL is called by FITTING via the Levenberg-Marquardt routine.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine pure_component_objective_fct ( lm_m_dat, lm_n_par, para, fvec, iflag )

  use PARAMETERS, only: NAV, machine_eps
  use eos_constants, only: load_dd_constants, load_qq_constants, load_dq_constants
  use module_eos_derivatives, only: z3t
  use properties, only: calculate_density, calculate_pressure, ztot_rho
  use pcsaft_pure_and_binary_parameters, only: dipole, qudpole, dipole_quad,  &
       mseg, sigma, epsilon_k, kap_hb, eps_hb, dipole_moment, quadru_moment,  &
       assoc, assoc_scheme, nhb_typ, nhb_no
  use input_output, only: get_mass_density
  use module_pure_diagrams, only: pure_vapor_pressure
  use module_critical_point, only: critical_point_pure
  use module_eos_derivatives, only: aggregate_polar_parameters

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: lm_m_dat
  integer, intent(in)                             :: lm_n_par
  real(dp), dimension(:), intent(in)              :: para
  real(dp), dimension(:), intent(in out)          :: fvec
  integer, intent(in out)                         :: iflag

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: weight, actual_weight
  real(dp)                                        :: eta_start
  real(dp)                                        :: eta, rho, rho_V
  real(dp)                                        :: pcalc, pcalc_z
  real(dp)                                        :: mass_density
  real(dp), dimension(ncomp)                      :: xF
  real(dp), dimension(ncomp)                      :: rhoi, mass_frac
  real(dp)                                        :: tc, pc, rhoc
  real(dp)                                        :: penalty
  logical                                         :: converg
  character (LEN=4)                               :: char_len
  !-----------------------------------------------------------------------------

  fvec(:) = 0.0_dp

  xF(1) = 1.0_dp

  if ( minval( para(:) ) < 0.0_dp ) then
     fvec(:) = 1.0_dp + abs( minval( para(:) ) ) * 1000.0_dp
     ! fvec(:) = fvec(:) * 2.0_dp ! penalty function for negative parameters
     write (*,*) 'warning negative pure component parameters encountered',para !(minloc(para(:)))
     return
  end if

  tc = 0.0_dp

  !-----------------------------------------------------------------------------
  ! assign pure component parameters from para-vector
  !-----------------------------------------------------------------------------
  do i = 1, lm_n_par

     if ( i_par(i) == 1 ) mseg(1) = para(i)
     if ( i_par(i) == 2 ) sigma(1) = para(i)
     if ( i_par(i) == 3 ) epsilon_k(1) = para(i)
     if ( i_par(i) == 4 ) kap_hb(1,1) = exp( - para(i) )
     if ( i_par(i) == 5 ) then
        if ( assoc_scheme(1) == '1A' ) then
           eps_hb(1,1,1,1) = para(i)
        else
           eps_hb(1,1,1,2) = para(i)
           eps_hb(1,1,2,1) = para(i)
        end if
     end if
     if ( i_par(i) == 6 ) dipole_moment(1) = para(i)
     if ( i_par(i) == 7 ) quadru_moment(1) = para(i)

  end do

  !-----------------------------------------------------------------------------
  ! load eos constants
  !-----------------------------------------------------------------------------
  if ( dipole_moment(1) > 0.0_dp .AND. .NOT.dipole ) then
     write (*,*) 'ERROR1',dipole, dipole_moment(1)
     stop
  end if
  if ( dipole_moment(1) < machine_eps .AND. dipole ) then
     write (*,*) 'ERROR2',dipole, dipole_moment(1)
     stop
  end if

  if ( dipole ) call load_dd_constants
  if ( qudpole ) call load_qq_constants
  if ( dipole_quad ) call load_dq_constants

  if ( dipole .OR. qudpole ) call aggregate_polar_parameters

  if ( outp > 0 ) write (*,*) mseg(1), sigma(1), epsilon_k(1)
  if ( assoc ) then
     if ( outp > 0 ) write (*,*) assoc_scheme(1), nhb_typ(1), nhb_no(1,1:nhb_typ(1))
     if ( outp > 0 ) write (*,*) kap_hb(1,1), eps_hb(1,1,1,1), eps_hb(1,1,1,2)
  end if

  !-----------------------------------------------------------------------------
  ! check parameter bounds
  !-----------------------------------------------------------------------------
  penalty = 0.0_dp
  if ( sigma(1) < 1.8_dp ) penalty = penalty + abs( sigma(1) - 1.8_dp )
  if ( sigma(1) > 7.0_dp ) penalty = penalty + abs( sigma(1) - 7.0_dp )
  if ( epsilon_k(1) > 600.0_dp ) penalty = penalty + abs( epsilon_k(1) - 600.0_dp ) / 100.0_dp
  if ( kap_hb(1,1) > 1.0_dp ) penalty = penalty + abs( kap_hb(1,1) - 1.0_dp ) * 10.0_dp
  if ( eps_hb(1,1,1,1) > 10.E3_dp ) penalty = penalty + abs( eps_hb(1,1,1,1) - 10.E3_dp ) / 10.E3_dp
  if ( eps_hb(1,1,1,2) > 10.E3_dp ) penalty = penalty + abs( eps_hb(1,1,1,1) - 10.E3_dp ) / 10.E3_dp
  if ( dipole_moment(1) > 8.0_dp ) penalty = penalty + abs( dipole_moment(1) - 8.0_dp )
  if ( quadru_moment(1) > 10.0_dp ) penalty = penalty + abs( quadru_moment(1) - 10.0_dp )

  if ( penalty > machine_eps ) then
     fvec(:) = 100.0_dp + penalty**2 * 100.0_dp 
     write (*,*) 'exceeding bounds of pure component parameters',para(:)
     return
  end if

  !-----------------------------------------------------------------------------
  ! liquid densit data
  !-----------------------------------------------------------------------------
  do i = 1, n_dense

     if ( vliq(i) < v_crit ) then
        eta_start = 0.5_dp
     else
        eta_start = 0.0001_dp
     end if

     if ( pliq(i) > 0.0_dp ) then

        !-----------------------------------------------------------------------
        ! case 1: calculate rho to given (p,T)
        !-----------------------------------------------------------------------
        call calculate_density ( tliq(i), pliq(i), xF, eta_start, eta, rho )

     else

        !-----------------------------------------------------------------------
        ! case 2, when pliq(i) = 0, indicating rho_L at coexistence
        !-----------------------------------------------------------------------
        call pure_vapor_pressure ( tliq(i), converg, pcalc, rho, rho_V )
        eta = rho * z3t

        if ( .NOT.converg ) then
           call calculate_density ( tliq(i), 1.E6_dp, xF, eta_start, eta, rho )
           if ( eta < 0.1_dp ) call calculate_density ( tliq(i), 1.E7_dp, xF, eta_start, eta, rho )
           write (*,*) 'vapor press. calc. for rho_L not converged',tliq(i),eta_start,eta
        end if

     end if

     if ( eta > 0.55_dp ) write (*,*)' high density ', i, eta

     rhoi(1) = rho
     call get_mass_density ( rhoi, mass_density, mass_frac )
     v_calc(i) = 1.0_dp / mass_density

     !--------------------------------------------------------------------------
     ! entry to deviation-vector (element of objective funciton)
     !--------------------------------------------------------------------------
     weight = 0.5_dp
     actual_weight = weight * sqrt( type_weight(1) /  ( vliq(i) * v_calc(i) ) )

     fvec(i) = actual_weight * ( v_calc(i) - vliq(i) )
     if ( outp > 0 )  write (*,*) 'dense ',i, ( v_calc(i) - vliq(i) ) / vliq(i)*100.0_dp

  end do



  !-----------------------------------------------------------------------------
  ! vapor pressure data
  !-----------------------------------------------------------------------------
  do  i = 1, n_lv

     call pure_vapor_pressure ( tlv(i), converg, plv_calc(i), rho, rho_V ) !, plv(i) )

     !--------------------------------------------------------------------------
     ! entry to deviation-vector (element of objective funciton)
     !--------------------------------------------------------------------------
     if ( converg ) then

        weight = 0.2_dp
        actual_weight = weight * sqrt( type_weight(2) ) / ( plv(i) * plv_calc(i) )**0.425_dp

        fvec( i + n_dense ) = actual_weight * abs( plv_calc(i) - plv(i) )
        if ( outp > 0 )  write (*,*) 'psat  ',i, ( plv_calc(i) - plv(i) ) / plv(i)*100.0_dp,tlv(i)

     else

        if ( tc < machine_eps ) then     ! has tc not already been calculated (?), then ...
           tc = tlv(i)
           call critical_point_pure ( tc, converg, pc, rhoc )
        end if
        fvec( i + n_dense ) = 5.0_dp * ( tlv(i) - tc )
        ! write (*,*) 'exp.boiling temp. above current crit. pt.', i, tlv(i), tc

     end if

  end do


  !-----------------------------------------------------------------------------
  ! 2nd virial coefficients
  !-----------------------------------------------------------------------------
  do  i = 1, n_virial

     call calculate_pressure ( t_vir(i), xF, 0.0_dp, rho_V, pcalc, pcalc_z )

     B_calc(i) = ztot_rho / 1.E30_dp * NAV * 1.E6_dp        ! in unit (cm**3/mol)

     weight = ( epsilon_k(1) / t_vir(i) )**3
     actual_weight = weight * type_weight(3) / sqrt( abs( b_vir(i) * B_calc(i) ) )
     fvec( i + n_dense + n_lv ) = actual_weight * ( B_calc(i) - b_vir(i) )
     ! write (*,*) 'B_vir ',i, ( B_calc(i) - b_vir(i) ) / b_vir(i) * 100.0_dp

  end do

  if ( ( n_dense + n_lv + n_virial ) /= lm_m_dat ) write (*,*) 'array length not matching!'

  !-----------------------------------------------------------------------------
  ! ensure mseg >= 1.0
  !-----------------------------------------------------------------------------
  if ( mseg(1) < 1.0_dp ) then
     fvec(:) = fvec(:) + fvec(:) * ( 1.0_dp - mseg(1) )**2 * 100.0_dp
     if ( iflag == 1 ) write (*,*) 'correcting for m<1', mseg(1)
  end if
  
  !-----------------------------------------------------------------------------
  ! output iteration progress to screen
  !-----------------------------------------------------------------------------
  write (char_len,'(I3)') lm_n_par + 1
  if (iflag == 1) write (*,'(a,'//char_len//'(G15.7))') ' parameters, error',( para(i), i=1,lm_n_par ), sum(abs(fvec(:)))
  if ( outp > 0 ) read (*,*)

end subroutine pure_component_objective_fct



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE read_pure_component_data
!>  \brief Read experimental pure component data
!!
!! Read experimental pure component data for the regression of pure
!! component parameters. This file also writes the experimental data to
!! the output file.
!! So far, the experimental data is stored in vectors of fixed length
!! (see module 'pure_fit_parameters'). At some point this should be
!! changed - allocate the vector length here.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine read_pure_component_data ( pure_fit_file, flag_exit, lm_m_dat )

  use pcsaft_pure_and_binary_parameters, only: mm

  !-----------------------------------------------------------------------------
  character (LEN=40), intent(in)                  :: pure_fit_file
  integer, intent(out)                            :: lm_m_dat
  logical, intent(inout)                          :: flag_exit

  !-----------------------------------------------------------------------------
  integer                                         :: i, ii, k
  integer                                         :: alldata
  integer                                         :: select_data
  integer                                         :: data_no
  integer, dimension(n_types_of_data)             :: datatype_available
  integer, dimension(n_types_of_data)             :: n_low, n_high
  character (LEN=80)                              :: info_line
  logical                                         :: filefound
  logical                                         :: repeat
  !-----------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------
  ! open input-file (22) for reading
  !-----------------------------------------------------------------------------
  inquire ( FILE='./input_file/pure_comp/'//pure_fit_file, EXIST = filefound )
  if ( filefound ) then
     open (22, FILE = './input_file/pure_comp/'//pure_fit_file )
  else
     write (*,*) ' FILE CANNOT BE OPENED ', './input_file/pure_comp/'//trim(adjustl(pure_fit_file))
     flag_exit = .true.
     return
  end if

  n_dense = 0
  n_lv = 0
  n_virial = 0

  datatype_available(:) = 0

  !-----------------------------------------------------------------------------
  ! reading info-lines (input-file) and writing them (output-file)
  !-----------------------------------------------------------------------------

  do  i = 1, 15
     read (22,'(A80)') info_line
     write (23,'(a)') info_line
  end do

  ! --- reading molecular mass and an estimate of the crit. density ------------
  read (22,'(A80)') info_line
  write (23,'(A)') info_line
  read (22,*) mm_pure, v_crit
  write (23, '(2(2x,G13.6))') mm_pure, v_crit
  mm(1) = mm_pure                        ! needed for calculating mass-densities

  do i= 1, 14
     read (22,'(A80)') info_line
     write (23,'(A)') info_line
  end do

  ! --- reading the number of experimental data points -------------------------
  read (22,'(A80)') info_line
  read (22,*) ( datatype_available(i), i = 1, n_types_of_data )
  write (*,*) ' '

  !-----------------------------------------------------------------------------
  ! selection of data sets by user
  !-----------------------------------------------------------------------------

  write (*,*) ' '
  write (*,*) ' number of data-sets:'
  write (*,'(T2,A40,I3)') 'PVT - data                :', datatype_available(1)
  write (*,'(T2,A40,I3)') 'vapor pressure data       :', datatype_available(2)
  if ( datatype_available(3) > 0 ) then
     write (*,'(T2,A40,I3)') '2nd virial coefficients   :', datatype_available(3)
  end if
  write (*,*) ' '

  ! --- option to select all data sets -----------------------------------------
  write (*,*) ' select data set !'
  write (*,*) '      consider all available data:     1'
  write (*,*) '      select data-sets individually:   0'
  read (*,*) alldata

  do i = 1, n_types_of_data

     if ( alldata == 1 ) then

        n_high(i) = datatype_available(i)
        if ( n_high(i) > 0 ) n_low(i) = 1

     else

        n_low(i) = 0
        n_high(i) = 0
        if ( i == 1) write (*,*) ' select data-types !'
        if ( i == 1) write (*,*) ' '
        if (datatype_available(i) > 0) then
           if (i == 1) write (*,*) ' PVT-data ?  (0/1) '
           if (i == 2) write (*,*) ' vapor pressure-data ?  (0/1) '
           if (i == 3) write (*,*) ' virial coeff. B(T)-data ?  (0/1) '
           read (*,*) select_data
           if ( select_data == 1 ) then
              if (i == 1) write (*,*) ' CHOOSE FROM PVT-DATA:'
              if (i == 2) write (*,*) ' CHOOSE FROM VAPOR PRESSURE DATA:'
              if (i == 3) write (*,*) ' CHOOSE FROM 2nd VIRIAL COEFFICIENT DATA:'
              repeat = .true.
              do while ( repeat )
                 write (*,*) ' Specify the data-sets to be considered'
                 write (*,*) ' Lower number of data-set (e.g.  1):'
                 read (*,*) n_low(i)
                 write (*,'(a,I3,a)') ' Upper number of data-set (e.g.',datatype_available(i),'):'
                 read (*,*) n_high(i)
                 ! --- check the input --------------------------------------------------
                 if ( (n_low(i) > n_high(i)) .OR. (n_low(i) < 1) .OR. (n_high(i) > datatype_available(i)) ) then
                    write (*,*) ' Erroneous input! The lower data-set number must be smaller'
                    write (*,*) ' than the upper value.'
                    write (*,*) ' The upper data-set number must not be greater than ',datatype_available(i)
                    write (*,*) ' Lower number of data-set (e.g.  1):'
                    read (*,*) n_low(i)
                    write (*,'(a,I3,a)') ' Upper number of data-set (e.g.',datatype_available(i),'):'
                    read (*,*) n_high(i)
                 else
                    repeat = .false.
                 end if
              end do
           end if
        end if

     end if

     write (*,*) ' '

     !--------------------------------------------------------------------------
     ! reading data sets from input file
     !--------------------------------------------------------------------------

     ! --- two headers ---------------------------------------------------------
     if ( datatype_available(i) /= 0 ) read (22,'(A80)') info_line
     if ( n_high(i) > 0 )  write (23,'(a)') info_line
     if ( datatype_available(i) /= 0 ) read (22,'(A80)') info_line
     if ( n_high(i) > 0 )  write (23,'(a)') info_line

     ! --- reading data sets if a certain data type i was chosen ---------------
     write (*,*) ' '
     if ( n_high(i) > 0 ) then

        ! --- reading data -----------------------------------------------------
        do k = 1, n_low(i) - 1
           read (22,*)
        end do

        if ( i == 1 .AND. n_high(i) > 0 ) then
           n_dense = n_high(i) - ( n_low(i) - 1 )
           allocate( pliq(n_dense), tliq(n_dense), vliq(n_dense), v_calc(n_dense) )
           write (*,*) i, n_low(i), n_high(i)
           do k = n_low(i), n_high(i)
              ii = k - ( n_low(i) - 1 )
              read (22,*) data_no, pliq(ii), tliq(ii), vliq(ii)
              write (*,'(i3,3(2x,G12.5))') data_no, pliq(ii), tliq(ii), vliq(ii)
              write (23,'(i3,3(2x,G12.5))') data_no, pliq(ii), tliq(ii), vliq(ii)
           end do
        else if ( i == 2 .AND. n_high(i) > 0 ) then
           n_lv = n_high(i) - ( n_low(i) - 1 )
           allocate( plv(n_lv), tlv(n_lv), plv_calc(n_lv) )
           do k = n_low(i), n_high(i)
              ii = k - ( n_low(i) - 1 )
              read (22,*) data_no, plv(ii), tlv(ii)
              write (*,'(i3,2(2x,G12.5))') data_no, plv(ii), tlv(ii)
              write (23,'(i3,2(2x,G12.5))') data_no, plv(ii), tlv(ii)
           end do
        else if ( i == 3 .AND. n_high(i) > 0 ) then
           n_virial = n_high(i) - ( n_low(i) - 1 )
           allocate( b_vir(n_virial), t_vir(n_virial), B_calc(n_virial) )
           do k = n_low(i), n_high(i)
              ii = k - ( n_low(i) - 1 )
              read (22,*) data_no, b_vir(ii), t_vir(ii)
              write (*,'(i3,2(2x,G12.5))') data_no, b_vir(ii), t_vir(ii)
              write (23,'(i3,2(2x,G12.5))') data_no, b_vir(ii), t_vir(ii)
           end do
        end if

        do k = n_high(i)+1, datatype_available(i)
           read (22,*)
        end do

     else

        ! --- the data type was not selected for the fitting -------------
        if (datatype_available(i) > 0) then
           do  k = 1, datatype_available(i)
              read (22,*)
           end do
        end if

     end if

  end do

  !-----------------------------------------------------------------------------
  ! total number of data sets considered in parameter fitting
  !-----------------------------------------------------------------------------

  lm_m_dat = n_dense + n_lv + n_virial          ! excluding crit. point

  close (22)

end subroutine read_pure_component_data



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE define_adjusted_pure_parameters
!> \brief Subroutine for reading in starting values of the pure component parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine define_adjusted_pure_parameters ( lm_n_par, para )

  use PARAMETERS, only: machine_eps
  use pcsaft_pure_and_binary_parameters, only: assoc, dipole, qudpole, dipole_quad,  &
       mseg, sigma, epsilon_k, kap_hb, eps_hb, dipole_moment, quadru_moment, assoc_scheme,  &
       nhb_typ, nhb_no

  !-----------------------------------------------------------------------------
  integer, intent(out)                            :: lm_n_par
  real(dp), allocatable, intent(out)              :: para(:)

  !-----------------------------------------------------------------------------
  integer                                         :: fitopt, position
  integer                                         :: i, j
  integer                                         :: i_assoc_scheme
  integer                                         :: default_starting_values
  logical                                         :: i_fit
  logical                                         :: assoc_scheme_defined
  !-----------------------------------------------------------------------------

  assoc = .false.
  dipole = .false.
  qudpole = .false.
  dipole_quad = .false.

  nhb_typ(:) = 0
  nhb_no(:,:) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! selecting adjustable parameters
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! assign pure component parameters from para-vector
  !-----------------------------------------------------------------------------
  write (*,*) ' '
  write (*,*) '*********** PARAMETER SPECIFICATION ***********'
  write (*,*) ' '

  write (*,*) ' You may choose from two default fitting options'
  write (*,*) ' or individually specify all parameters.'
  write (*,*) ' (Dipolar and quadrupolar moments are set to zero'
  write (*,*) ' for both default options; for dipolar or'
  write (*,*) ' quadrupolar compounds choose option 3)'
  write (*,*) '  '
  write (*,*) ' Choose one the following options'
  write (*,*) ' Default 1: Fitting 3 parameters'
  write (*,*) '            * segment number           m        [ ]'
  write (*,*) '            * segment diameter         sigma    [A]'
  write (*,*) '            * segment energy parameter eps/k    [K]'
  write (*,*) ' Default 2: Fitting 5 parameters for associating'
  write (*,*) '            compounds with 2 association sites'
  write (*,*) '            * segment number           m        [ ]'
  write (*,*) '            * segment diameter         sigma    [A]'
  write (*,*) '            * segment energy parameter eps/k    [K]'
  write (*,*) '            * assoc.-volume parameter  kappa_AB [ ]'
  write (*,*) '            * assoc.-energy parameter  eps_AB/k [K]'
  write (*,*) ' Option  3: Individual specification of all'
  write (*,*) '            parameters to be fitted.'
  write (*,*) ' Choose 1, 2, or 3.'
  read (*,*) fitopt
  ! fitopt = 3  ! JG

  if ( fitopt == 1 ) then
     lm_n_par = 3
     allocate( para( lm_n_par ), i_par( lm_n_par ), aa_txt( lm_n_par ) )
     aa_txt(1) = ' Segment number m       :'
     aa_txt(2) = ' Segment size par. sigma:'
     aa_txt(3) = ' Segment energy eps/k   :'
     i_par(1) = 1
     i_par(2) = 2
     i_par(3) = 3
  else if ( fitopt == 2 ) then
     lm_n_par = 5
     allocate( para( lm_n_par ), i_par( lm_n_par ), aa_txt( lm_n_par ) )
     aa_txt(1) = ' Segment number m       :'
     aa_txt(2) = ' Segment size par. sigma:'
     aa_txt(3) = ' Segment energy eps/k   :'
     aa_txt(4) = ' Assoc.-volume kappa_AB :'
     aa_txt(5) = ' Assoc.-energy eps_AB/k :'
     i_par(1) = 1
     i_par(2) = 2
     i_par(3) = 3
     i_par(4) = 4
     i_par(5) = 5
     assoc_scheme(1) = '2B'
     nhb_typ(1) = 2
     nhb_no(1,1:2) = 1.0_dp
     assoc = .true.
  else
     write (*,*) ' The following pure component parameters can be fitted'
     write (*,*) ' (1) segment number                  m        [ ]'
     write (*,*) ' (2) segment diameter                sigma    [A]'
     write (*,*) ' (3) segment energy parameter        eps/k    [K]'
     write (*,*) ' (4) assoc.-volume parameter         kappa_AB [ ]'
     write (*,*) ' (5) assoc.-energy parameter         eps_AB/k [K]'
     write (*,*) ' (6) dipolar moment                  my       [D]'
     write (*,*) ' (7) quadrupolar moment              Q        [D]'
     write (*,*) ' '
     write (*,*) ' Note, that the number of association sites has'
     write (*,*) ' to be user-specified'
     write (*,*) ' '
     write (*,*) ' Choose the total number of fitting-parameters'
     read (*,*) lm_n_par
     ! lm_n_par = 4  ! JG
     allocate( para( lm_n_par ), i_par( lm_n_par ), aa_txt( lm_n_par ) )
     write (*,*) ' Choose the position-number of fitting-param.'
     write (*,*) ' (see above list of parameter 1 to 7)'
     do i = 1, lm_n_par
        write (*,'(a,i1,a,i2)') '  position-number for param. ',i,' of ',lm_n_par
        read (*,*) position
        if ( position == 1 ) aa_txt(i) = ' Segment number m       :'
        if ( position == 2 ) aa_txt(i) = ' Segment size par. sigma:'
        if ( position == 3 ) aa_txt(i) = ' Segment energy eps/k   :'
        if ( position == 4 .or. position == 5 ) assoc = .true.
        if ( position == 4 ) aa_txt(i) = ' Assoc.-volume kappa_AB :'
        if ( position == 5 ) aa_txt(i) = ' Assoc.-energy eps_AB/k :'
        if ( position == 6 ) dipole = .true.
        if ( position == 6 ) aa_txt(i) = ' Dipolar moment         :'
        if ( position == 7 ) qudpole = .true.
        if ( position == 7 ) aa_txt(i) = ' Quadrupolar moment     :'
        i_par(i) = position
     end do
     if ( dipole .AND. qudpole ) dipole_quad = .true.
  end if

  !-----------------------------------------------------------------------------
  ! specifying starting values of adjustable parameters
  !-----------------------------------------------------------------------------

  write (*,*) ' '
  write (*,*) ' Input of starting values for fitting-parameters:'
  write (*,*) ' '

  write (*,*) ' GIVE STARTING VALUES FOR PARAMETERS  '
  write (*,*) ' '
  write (*,*) ' choose 0:   for manually specifying starting values'
  write (*,*) ' choose 1:   for default starting values'
  write (*,*) '             (sigma=3.7A, eps/k=250.K, m=0.04*M, kappa_hb=0.03, eps_hb/k=1800.K'
  write (*,*) '             kappa_hb=0.03, eps_hb/k=1800.K, dipole_moment=2.D, Q_moment=2.DA)'
  read (*,*) default_starting_values
  ! default_starting_values = 0   ! JG

  do i = 1, lm_n_par
     write (*,'(2a,i2,a,i2,a)') aa_txt(i),' (param.',i,' of',lm_n_par,')'
     if ( default_starting_values == 1 ) then
        if ( i_par(i) == 1 ) para(i) = 0.04_dp * mm_pure
        if ( i_par(i) == 2 ) para(i) = 3.4_dp
        if ( i_par(i) == 3 ) para(i) = 270.0_dp
        if ( i_par(i) == 4 ) para(i) = 0.03_dp
        if ( i_par(i) == 5 ) para(i) = 1800.0_dp
        if ( i_par(i) == 6 ) para(i) = 2.0_dp
        if ( i_par(i) == 7 ) para(i) = 2.0_dp
     else
        read (*,*) para(i)
        ! if ( i == 1 ) para(i) = 2.0_dp        ! JG
        ! if ( i == 2 ) para(i) = 3.4_dp        ! JG
        ! if ( i == 3 ) para(i) = 280.0_dp      ! JG
        ! if ( i == 4 ) para(i) = 0.03_dp       ! JG
        ! if ( i == 5 ) para(i) = 3500.0_dp     ! JG
     end if

     if ( i_par(i) == 4 ) para(i) = - log( para(i) )
  end do
  write (*,*) ' '

  !-----------------------------------------------------------------------------
  ! for case, where param. are selected individually, determine rest
  !-----------------------------------------------------------------------------
  assoc_scheme_defined = .false.

  if ( fitopt >= 3 ) then
     write (*,*) ' '
     write (*,*) ' The remaining parameters (not-fitted) now'
     write (*,*) ' have to be defined:'
     write (*,*) ' '
     do i = 1, 7
        i_fit = .false.
        do j = 1, lm_n_par
           if ( i_par(j) == i ) i_fit = .true.
        end do
        if ( .NOT.i_fit ) then
           if ( i == 1 ) then
              !write (*,*) ' Give segm.No.m / Molar Mass :'
              write (*,*) ' Give number of segments, m :'
              read (*,*) mseg(1)
           else if ( i == 2 ) then
              write (*,*) ' Give segment size parameter, sigma / [A]  :'
              read (*,*) sigma(1)
           else if ( i == 3 ) then
              write (*,*) ' Give segment energy parameter eps/k / [K]:'
              read (*,*) epsilon_k(1)
           else if ( i == 4 .OR. i == 5 ) then
              if ( .NOT.assoc_scheme_defined ) then
                 assoc_scheme_defined = .true.
                 write (*,*) ' Give association scheme:'
                 write (*,*) ' choose 0:   for non-associating compounds'
                 write (*,*) ' choose 1:   for 1 association site,  scheme 1A'
                 write (*,*) ' choose 2:   for 2 association sites, scheme 2B'
                 write (*,*) ' choose 3:   for 3 association sites, scheme 3B'
                 write (*,*) ' choose 4:   for 4 association sites, scheme 4C'
                 read (*,*) i_assoc_scheme
                 ! i_assoc_scheme = 1   ! JG
                 if ( i_assoc_scheme == 0 ) then
                    assoc = .false.
                 else
                    assoc = .true.
                    if ( i_assoc_scheme == 1 ) then
                       assoc_scheme(1) = '1A'
                       nhb_typ(1) = 1
                       nhb_no(1,1) = 1.0_dp
                    else if ( i_assoc_scheme == 2 ) then
                       assoc_scheme(1) = '2B'
                       nhb_typ(1) = 2
                       nhb_no(1,1:2) = 1.0_dp
                    else if ( i_assoc_scheme == 3 ) then
                       assoc_scheme(1) = '3B'
                       nhb_typ(1) = 2
                       nhb_no(1,1) = 1.0_dp
                       nhb_no(1,2) = 2.0_dp
                    else if ( i_assoc_scheme == 4 ) then
                       assoc_scheme(1) = '4C'
                       nhb_typ(1) = 2
                       nhb_no(1,1:2) = 2.0_dp
                    end if
                 end if
              end if
              if ( i == 4 .AND. assoc ) then
                 write (*,*) ' Give association volume kappa_AB :'
                 read (*,*) kap_hb(1,1)
              end if
              if ( i == 5 .AND. assoc ) then
                 write (*,*) ' Give assoc.-energy eps_AB/k :'
                 if ( i_assoc_scheme == 1 ) then
                    read (*,*) eps_hb(1,1,1,1)
                 else
                    read (*,*) eps_hb(1,1,1,2)
                    eps_hb(1,1,2,1) = eps_hb(1,1,1,2)
                 end if
              end if
           else if ( i == 6 ) then
              write (*,*) ' Give dipolar moment         :'
              read (*,*) dipole_moment(1)
              if ( abs( dipole_moment(1) ) > machine_eps ) dipole = .true.
           else if ( i == 7 ) then
              write (*,*) ' Give quadrupolar moment     :'
              read (*,*) quadru_moment(1)
              ! quadru_moment(1) = 0.0_dp    ! JG
              if ( abs( quadru_moment(1) ) > machine_eps ) qudpole = .true.
           end if
        else
           if ( .NOT.assoc_scheme_defined ) then
              assoc_scheme_defined = .true.
              write (*,*) ' Give association scheme:'
              write (*,*) ' choose 0:   for non-associating compounds'
              write (*,*) ' choose 1:   for 1 association site,  scheme 1A'
              write (*,*) ' choose 2:   for 2 association sites, scheme 2B'
              write (*,*) ' choose 3:   for 3 association sites, scheme 3B'
              write (*,*) ' choose 4:   for 4 association sites, scheme 4C'
              read (*,*) i_assoc_scheme
              if ( i_assoc_scheme == 0 ) then
                 write (*,*) 'that is an unexpected selection. stop'
                 assoc = .false.
              else
                 if ( i_assoc_scheme == 1 ) then
                    assoc_scheme(1) = '1A'
                    nhb_typ(1) = 1
                    nhb_no(1,1) = 1.0_dp
                 else if ( i_assoc_scheme == 2 ) then
                    assoc_scheme(1) = '2B'
                    nhb_typ(1) = 2
                    nhb_no(1,1:2) = 1.0_dp
                 else if ( i_assoc_scheme == 3 ) then
                    assoc_scheme(1) = '3B'
                    nhb_typ(1) = 2
                    nhb_no(1,1) = 1.0_dp
                    nhb_no(1,2) = 2.0_dp
                 else if ( i_assoc_scheme == 4 ) then
                    assoc_scheme(1) = '4C'
                    nhb_typ(1) = 2
                    nhb_no(1,1:2) = 2.0_dp
                 end if
              end if
           end if
        end if
     end do
  end if

  if ( dipole .AND. qudpole ) dipole_quad = .true.

  !-----------------------------------------------------------------------------
  ! weighting for the different data types (currently all =1.0_dp)
  !-----------------------------------------------------------------------------

  write (*,*) ' '
  write (*,*) 'Input of weighting factors for the different data-sets'
  write (*,*) ' '

  if ( n_dense > 0 ) then
     write (*,*) 'weights for density data ?'
     type_weight(1) = 1.0_dp
     write (*,*) type_weight(1)
     ! read (*,*) type_weight(1)
  end if
  if ( n_lv > 0 ) then
     write (*,*) 'weights for vapor pressure data ?'
     type_weight(2) = 1.0_dp
     write (*,*) type_weight(2)
     ! read (*,*) type_weight(2)
  end if
  if ( n_virial > 0 ) then
     write (*,*) 'weights for 2nd virial coefficients ?'
     type_weight(3) = 0.01_dp
     write (*,*) type_weight(3)
     ! read (*,*) type_weight(3)
  end if

  ! --- summarizing weighting of the different data types ----------------------
  write (*,*) 'the following weights were specified:'
  write (*,*) ' '
  if ( n_dense > 0 )  write (*,*) ' density data            : ', type_weight(1)
  if ( n_lv > 0 )     write (*,*) ' vapor pressure data     : ', type_weight(2)
  if ( n_virial > 0 ) write (*,*) ' 2nd virial coefficients : ', type_weight(3)

  !-----------------------------------------------------------------------------
  ! summarizing the starting values and the scaling of parameters
  !-----------------------------------------------------------------------------

  write (*,*) ' '
  write (23,*) ' '
  write (23,*) ' '
  write (23,'(a)') ' PARAMETER-TYPE           PARAM.-START   SCALING         para'
  write (23,'(a)') ' --------------------------------------------------------------'
  do i=1,lm_n_par
     write (23,'(a,t26,G13.6)') aa_txt(i), para(i)
  end do
  write (23,*) ' '
  write (23,'(a)') ' WEIGHTS'
  write (23,'(a)') ' -------'
  if ( n_dense > 1 )  write (23,'(A,T26,F5.2)')' density data           :', type_weight(1)
  if ( n_lv > 1 )     write (23,'(A,T26,F5.2)')' vapor pressure data    :', type_weight(2)
  if ( n_virial > 1 ) write (23,'(A,T26,F5.2)')' 2nd virial coeff.      :', type_weight(3)
  write (23,*) ' '

end subroutine define_adjusted_pure_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pure_output
!> \brief output result of pure component parameter regression.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine output_pure_component_parameters ( pure_fit_file, info, lm_n_par, para )

  use pcsaft_pure_and_binary_parameters, only: assoc, dipole, qudpole,  &
       mseg, sigma, epsilon_k, kap_hb, eps_hb, dipole_moment, quadru_moment, assoc_scheme,  &
       nhb_typ, nhb_no

  !-----------------------------------------------------------------------------
  character (LEN=40), intent(in)                  :: pure_fit_file
  integer, intent(in)                             :: info
  integer, intent(in)                             :: lm_n_par
  real(dp), dimension(:), intent(in)              :: para

  !-----------------------------------------------------------------------------
  integer                                         :: i, j
  integer                                         :: ps
  real(dp)                                        :: devi, devimax, devisum, deviav, rms
  real(dp)                                        :: p_deviav, p_devimax, r_deviav, r_devimax
  real(dp)                                        :: assoc_eps
  character (LEN=4)                               :: char_len
  character (LEN=140)                             :: string1
  character (LEN=9)                               :: string2
  character (LEN=50)                              :: pure_fit_file_in
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! algorithm output: optimization status
  !-----------------------------------------------------------------------------

  write (*,*)  ' '
  write (23,*) ' '

  select case (info)
  case (:-1)
     write (23,*) 'SOLVER STATUS: Users FCN returned INFO = ', -info
     write (*, *) 'SOLVER STATUS: Users FCN returned INFO = ', -info
  case (0)
     write (23,*) 'SOLVER STATUS: Improper values for input parameters'
     write (*, *) 'SOLVER STATUS: Improper values for input parameters'
  case (1:3)
     write (23,*) 'SOLVER STATUS: Convergence'
     write (*, *) 'SOLVER STATUS: Convergence'
  case (4)
     write (23,*) 'SOLVER STATUS: Residuals orthogonal to the Jacobian'
     write (*, *) 'SOLVER STATUS: Residuals orthogonal to the Jacobian'
     write (*, *) 'There may be an error in FCN'
  case (5)
     write (23,*) 'SOLVER STATUS: Too many calls to FCN'
     write (*, *) 'SOLVER STATUS: Too many calls to FCN'
     write (*, *) 'Either slow convergence, or an error in FCN'
  case (6:7)
     write (23,*) 'SOLVER STATUS: TOL was set too small'
     write (*, *) 'SOLVER STATUS: TOL was set too small'
  case DEFAULT
     write (23,*) 'SOLVER STATUS: INFO =', info, ' ???'
     write (*, *) 'SOLVER STATUS: INFO =', info, ' ???'
  end select

  write (*,*)  '-------------'
  write (23,*) '-------------'

  !-----------------------------------------------------------------------------
  ! writing resulting parameters to screen and to output file
  !-----------------------------------------------------------------------------

  write (*,*) ' '
  write (*,*) '-----------------------------------------------'
  do i = 1, lm_n_par
     if ( i_par(i) == 4 ) then
        write (*,'(a,2x,f12.6)') aa_txt(i), exp( - para(i) )
     else
        write (*,'(a,2x,f12.6)') aa_txt(i), para(i)
     end if
  end do
  !WRITE (*,'(a,3(2x,f15.5))') ' STD ',(stdval(i), i=1,3)
  write (*,*) '-----------------------------------------------'

  write (23,*) ' '
  do i = 1, lm_n_par
     write (23,'(T2,A,I2,T26,G16.9)') 'PARAMETER ', i, para(i)
  end do


  !-----------------------------------------------------------------------------
  ! writing resulting parameters in the format of para_input.f90
  !-----------------------------------------------------------------------------
  if ( assoc ) then
     if ( assoc_scheme(1) == '1A' ) then
        assoc_eps = eps_hb(1,1,1,1)
     else
        assoc_eps = eps_hb(1,1,1,2)
     end if
  end if     

  write (23, *) ' '
  write (23, *)'-----------------------------------------------'
  write (23, *)'       mm(i)       = ', mm_pure
  do i = 1, 7
     ps = 0
     do j = 1, lm_n_par
        if ( i_par(j) == i ) ps = j
     end do

     if ( ps >= 1 ) then
        if ( i == 1 ) mseg(1) = para(ps)
        if ( i == 2 ) sigma(1) = para(ps)
        if ( i == 3 ) epsilon_k(1) = para(ps)
        if ( i == 4 ) kap_hb(1,1) = exp( - para(ps) )
        if ( i == 5 ) assoc_eps = para(ps)
        if ( i == 6 ) dipole_moment(1) = para(ps)
        if ( i == 7 ) quadru_moment(1) = para(ps)
     end if

  end do

  write (23, *)'       parame(i,1) = ', mseg(1)
  write (23, *)'       parame(i,2) = ', sigma(1)
  write (23, *)'       parame(i,3) = ', epsilon_k(1)

  if ( assoc ) then       

     write (23, *)'       nhb_typ(i)     = ', nhb_typ(1)
     if ( assoc_scheme(1) == '1A' ) then
        write (23, *)'       nhb_no(i,1)    = ', nhb_no(1,1)
        write (23, *)'       kap_hb(i,i)    = ', kap_hb(1,1)
        write (23, *)'       eps_hb(i,i,1,1)= ', assoc_eps
     else
        write (23, *)'       nhb_no(i,1)    = ', nhb_no(1,1)
        write (23, *)'       nhb_no(i,2)    = ', nhb_no(1,2)
        write (23, *)'       kap_hb(i,i)    = ', kap_hb(1,1)
        write (23, *)'       eps_hb(i,i,1,2)= ', assoc_eps
        write (23, *)'       eps_hb(i,i,2,1)= ', assoc_eps
        write (23, *)'       eps_hb(i,i,1,1)=  0.0'
        write (23, *)'       eps_hb(i,i,2,2)=  0.0'
     end if

  end if

  if ( dipole ) then
     write (23, *)'       parame(i,6) = ', dipole_moment(1)
  end if

  if ( qudpole ) then
     write (23, *)'       parame(i,8) = ', quadru_moment(1)
  end if

  write (23, *) '-----------------------------------------------'
  write (*,*) ' '
  write (23, *)


  !-----------------------------------------------------------------------------
  ! comparing calculated values to experimental data
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! 1. density data
  !-----------------------------------------------------------------------------

  if ( n_dense > 0 ) then
     write (23,*) '------------- fluid density data --------------'
     write (23,'(/A)') 'DATNR   AD (%)     T/K    v_exp/(m3/kg)   v_calc/(m3/kg) '
     devisum = 0.0_dp
     devimax = 0.0_dp
     devi = 0.0_dp
     rms  = 0.0_dp
     ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation -
     do i = 1, n_dense
        devi = abs( v_calc(i) - vliq(i) ) / sqrt( vliq(i)*v_calc(i) ) * 100.0_dp
        rms = rms + ( ( v_calc(i) - vliq(i) ) / vliq(i) )**2
        if ( devi > devimax ) devimax = devi
        devisum = devisum + devi
        write (23,'(i3,2x,G10.3,3(2x,G12.5))') i, devi, tliq(i), vliq(i), v_calc(i)
     end do
     deviav = devisum / real( n_dense, KIND=dp )
     rms = sqrt( rms / real( n_dense, KIND=dp ) )*100.0_dp
     write (23,*) ' '
     write (23,'(A)') 'Deviation of calculated to exp. volume data'
     write (23,*) ' '
     write (23,'(A, t35, F7.3, A)') 'Max. Deviation:',devimax, ' %'
     write (23,'(A, t35, F7.3, A)') 'Average Deviation: ',deviav, ' %'
     write (23,'(A, t35, F7.3, A)') 'RMS: ',rms, ' %'
     write (23,*) ' '
     r_devimax = devimax
     r_deviav  = deviav
  end if

  !-----------------------------------------------------------------------------
  ! 2. vapor pressure data
  !-----------------------------------------------------------------------------

  if ( n_lv > 0 ) then
     write (23,*) '------------- vapor pressure data -------------'
     write (23,'(/A)') 'DATNR  AD (%)     T/K        PEXP/PA         PCAL/PA'
     devisum = 0.0_dp
     devimax = 0.0_dp
     devi = 0.0_dp
     rms  = 0.0_dp
     do i = 1, n_lv
        ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation
        devi = abs( plv_calc(i) - plv(i) )/ sqrt( plv(i) * plv_calc(i) ) * 100.0_dp
        rms = rms + ( (plv_calc(i) - plv(i) ) / plv(i) )**2
        if ( devi > devimax ) devimax = devi
        devisum = devisum + devi
        write (23,'(i3,2x,G10.3,2x,G12.5,2(2x,G14.7))') i, devi, tlv(i), plv(i), plv_calc(i)
     end do
     deviav = devisum / real( n_lv, KIND=dp )
     rms = sqrt( rms / real( n_lv, KIND=dp ) ) * 100.0_dp
     write (23,*) ' '
     write (23,'(A)') 'Deviation of calculated to exp. vapor pressure data'
     write (23,*) ' '
     write (23, '(3(A, t35, F7.3, A/))') 'Max. Deviation:',devimax, ' %',  &
          'Average Deviation: ',deviav, ' %', 'RMS: ',rms, ' %'
     write (23,*) ' '
     p_devimax = devimax
     p_deviav  = deviav
  end if

  !-----------------------------------------------------------------------------
  ! 3.: 2nd virial coefficient data
  !-----------------------------------------------------------------------------

  if ( n_virial > 0 ) then
     write (23,*) '---------------2nd virial coeff. --------------'
     write (23,'(/A,A)') 'DATNR  AD (%)   T/K     B_exp     B_calc'
     devisum = 0.0_dp
     devimax = 0.0_dp
     devi = 0.0_dp
     rms  = 0.0_dp
     ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation -
     do  i = 1, n_virial
        devi = abs( B_calc(i) - b_vir(i) ) / (sqrt(abs(b_vir(i)))*sqrt(abs(B_calc(i))))*100.0_dp
        rms = rms + ((B_calc(i)-b_vir(i))/b_vir(i))**2
        if (devi > devimax) devimax = devi
        devisum = devisum + devi
        write (23,'(i3,2x,4(2x,G12.5))') i, devi, t_vir(i), b_vir(i), B_calc(i)
     end do
     deviav = devisum / real( n_virial, KIND=dp )
     rms = sqrt( rms / real( n_virial, KIND=dp ) ) * 100.0_dp
     write (23,*) ' '
     write (23,'(A)') 'Deviation of calculated to exp. 2nd virial coeff. data'
     write (23,*) ' '
     write (23, '(3(A, t35, F6.2, A/))') 'Max. Deviation:',devimax, ' %',  &
          'Average Deviation: ',deviav, ' %', 'RMS: ',rms, ' %'
     write (23,*) ' '
  end if

  write (23, *)' '
  pure_fit_file_in ='./input_file/pure_comp/'//trim( adjustl( pure_fit_file ) )
  write (string1,'(a,6G15.8)') pure_fit_file_in, mm_pure, mseg(1), sigma(1), epsilon_k(1), &
          dipole_moment(1), quadru_moment(1)
  if ( assoc ) then
     write (string2,'(i4,2a)') nhb_typ(1), '   ', assoc_scheme(1)
     if ( assoc_scheme(1) == '1A' ) then
        write (23,'(2a,i4,2G15.8)') string1, string2, nint(nhb_no(1,1)), kap_hb(1,1), eps_hb(1,1,1,1)
     else
        write (23,'(2a,2i4,2G15.8)') string1, string2, nint(nhb_no(1,1:2)), kap_hb(1,1), eps_hb(1,1,1,2)
     end if
  else
     write (23,'(a)') string1
  end if
  write (char_len,'(I3)') lm_n_par+5
  write (23,'(a,'//char_len//'G15.8)') pure_fit_file_in, mm_pure,  &
       ( para(i), i= 1, lm_n_par ),r_deviav,r_devimax,p_deviav,p_devimax
  write (23, *)' '

end subroutine output_pure_component_parameters

end module pure_fit_parameters
