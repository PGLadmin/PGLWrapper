      module load_bdd_mod
            
      implicit none
            
      real*8, parameter :: null_real = 8888.d0
      integer, parameter :: null_int = 8888
      character(4), parameter :: null_char = '8888'
            
      logical :: zero_corr, n_corr
            
      !> Module permettant de charger la dippr en mémoire
            
      character(50), parameter :: bdd_filename = "dippr.csv"
      integer, parameter :: nb_ids = 2
      integer, parameter :: nb_cstvalues = 7
      integer, parameter :: nb_tpty = 15
      integer, parameter :: nb_coeff = 7
            
      integer, parameter :: RK_idx = 1
      integer, parameter :: PR_idx = 2
                
      !   indices modif APM
      integer, parameter :: mw_idx = 1    !masse molaire
      integer, parameter :: tc_idx = 2    !Temperature critique
      integer, parameter :: pc_idx = 3    !Pression critique
      integer, parameter :: acen_idx = 4  !Facteur acentrique
      integer, parameter :: vc_idx = 5    !Volume critique
      integer, parameter :: zc_idx = 6    !Fac. comp. ctritique
      integer, parameter :: zr_idx = 7    !Z Rackett
    
      type :: bdd_ids
          character(200) :: val = ''
      end type
    
      type :: bdd_cstvalues
          real*8 :: val = null_real
          real*8 :: err = null_real
      end type
    
      type :: bdd_eq
          integer :: eqnid = null_int
          real*8 :: coeff(nb_coeff) = null_real
          real*8 :: tc = null_real
          character(500) :: note = null_char
        contains
          procedure :: y => bdd_y
      end type
            
      type :: bdd_molecules
          character(500) :: note
          type(bdd_ids) :: id(nb_ids)
          type(bdd_cstvalues) :: cstvalue(nb_cstvalues)
          type(bdd_eq), allocatable :: eq(:)
          real*8 :: parTwu(2,3)
          real*8 :: trans(2)
          logical :: valable
      end type
    
      type :: bdd_files
          integer :: nb_row
          type(bdd_molecules), allocatable :: molecule(:)
      contains
          procedure :: read_bdd_file
      end type
    
      type(bdd_files) :: bdd
    
      contains
            
        subroutine read_bdd_file(self)
            
                    class(bdd_files) :: self
                    integer :: chemid = null_int, chemid_max
                    integer :: ioerr
                    integer :: i, j, k ,l
                    integer, allocatable :: chemid_col(:)
                    integer, allocatable :: chemid_uniq(:)
                    type(bdd_molecules) :: null_molecule
            
                    write(*,'(A$)')'* Loading DIPPR'
            
                    open(42, file = trim(bdd_filename))
            
                    read(42,*) ! le séparateur du fichier *.csv
                    read(42,*) ! la ligne de titres des colonnes
            
                    ! recherche du chemid max dans le fichier et comptage du nombre de lignes
                    chemid_max = 0
                    do
                        read(42,*, iostat = ioerr)chemid
                        if (ioerr /= 0) then
                            exit
                        else
                            chemid_max = max(chemid, chemid_max)
                            self%nb_row = self%nb_row + 1
                        end if
                    end do
            
                    allocate(self%molecule(chemid_max))
            
                    ! stockage de la colonne contennant les chemid
                    allocate(chemid_col(self%nb_row))
            
                    rewind(42)
                    read(42,*)
                    read(42,*)
            
                    do i = 1, self%nb_row
                        read(42,*)chemid_col(i)
                    end do
            
                    ! lecture de la bdd complète
                    rewind(42)
                    read(42,*)
                    read(42,*)
            
                    do i = 1, self%nb_row
            
                        read(42,*, iostat = ioerr)chemid, & !ChemID DIPPR
                            & (self%molecule(chemid)%id(j),j = 1,nb_ids), & !Nom complet, shortname
            
                        & (self%molecule(chemid)%cstvalue(j)%val, &
                            & j = 1, nb_cstvalues), &
            
                        &  (self%molecule(chemid)%parTwu(RK_idx,j), j = 1,3), &
                        &  self%molecule(chemid)%trans(RK_idx), &
                       
                        &  (self%molecule(chemid)%parTwu(PR_idx,j), j = 1,3), &
                        &  self%molecule(chemid)%trans(PR_idx), &            
                            
                            & self%molecule(chemid)%eq(1)%eqnid, &
                            & (self%molecule(chemid)%eq(1)%coeff(j), j = 1, nb_coeff)
            
                        !if (chemid_uniq(i) /= chemid) stop '-> Probleme dans la lecture de la DIPPR !'
            
                    end do
            
                    close(42)
            
                    write(*,*)'-> OK'
            
                end subroutine
            
                !------
            
                function bdd_y(self, t) result(y)
            
                    class(bdd_eq) :: self
                    real*8 :: t, tc, tr, tt
                    real*8 :: y
                    real*8 :: a, b, c, d, e, f, g
            
                    a = self%coeff(1)
                    b = self%coeff(2)
                    c = self%coeff(3)
                    d = self%coeff(4)
                    e = self%coeff(5)
                    f = self%coeff(6)
                    g = self%coeff(7)
            
                    tc = self%tc
            
                    selectcase(self%eqnid)
                        case(100)
                            y = a + b * t + c * t**2 + d * t**3 + e * t**4
                        case(101)
                            y = exp(a + b / t + c * log(t) + d * t**e)
                        case(102)
                            y = a * t**b / (1.0_dp + c / t + d / t**2)
                        case(103)
                            y = a + b * exp(-c / t**d)
                        case(104)
                            y = a + b / t + c / t**3 + d / t**6 + e / t**9
                        case(105)
                            y = a / (b**(1.0_dp + (1.0_dp - t / c)**d))
                        case(106)
                            tr = t / tc
                            y = a * (1.0_dp - tr)**(b + c * tr + d * tr**2 + e * tr**3)
                            !write(*,*)tr
                        case(107)
                            y = a + b * (c / t / (sinh(c / t)))**2 + d * (e / t / (cosh(e / t)))**2
                        case(114)
                            tt = 1.0_dp - t / tc
                            y = a**2 / tt + b - 2.0_dp * a * c * tt - a * d * tt**2 - c**2 * tt**3 / 3.0_dp &
                                & - c * d * tt**4 / 2.0_dp - d**2 * tt**5 / 5.0_dp
                        case(115)
                            y = a + b / t + c * log(t) + d * t**2 + e / t**2
                        case(116)
                            tr = t / tc
                            y = a + b * (1.0_dp - tr)**0.35_dp + c * (1.0_dp - tr)**(2.0_dp/3.0_dp) &
                                & + d * (1.0_dp - tr) + e * (1.0_dp - tr)**(4.0_dp/3.0_dp)
                        case(119)
                            tt = 1.0_dp - t / tc
                            y = a + b * tt**(1.0_dp/3.0_dp) + c * tt**(2.0_dp/3.0_dp) + d * tt**(5.0_dp/3.0_dp) &
                                & + e * tt**(16.0_dp/3.0_dp) + f * tt**(43.0_dp/3.0_dp) + g * tt**(110.0_dp/3.0_dp)
                        case(123)
                            tt = 1.0_dp - t / tc
                            y = a * (1.0_dp + b * tt**(1.0_dp/3.0_dp) + c * tt**(2.0_dp/3.0_dp) + d * tt)
                        case(124)
                            tt = 1.0_dp - t / tc
                            y = a + b / tt + c * tt + d * tt**2 + e * tt**3
                        case(127)
                            y = a + b * (c/t)**2 * exp(c/t) / (exp(c/t) - 1)**2 &
                                & + d * (e/t)**2 * exp(e/t) / (exp(e/t) - 1)**2 &
                                & + f * (g/t)**2 * exp(g/t) / (exp(g/t) - 1)**2
                        case default
                            y = huge(1.0_dp)
                    end select
            
                end function
            
                !------
            
                subroutine load_bdd()
            
                    ! lecture du fichier contenant toute la dippr
                    call bdd%read_bdd_file()
            
                end subroutine
            
                !------
            
                function nan() result(x)
            
                    real*8 :: x
            
                    x = 0.0_dp
                    x = 0.0_dp / x
            
                end function
            
            end module
            
            !------
