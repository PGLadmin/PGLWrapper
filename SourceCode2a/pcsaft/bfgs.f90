!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
! This file contains the LBFGS algorithm and supporting routines
!  
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module moduel_lbfgs

  use PARAMETERS, only: dp
  implicit none

  public :: LBFGS, BFGS_optimizer

contains

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
! N       is an INTEGER variable that must be set by the user to the
!         number of variables. It is not altered by the routine.
!         Restriction: N>0.
! 
! M       is an INTEGER variable that must be set by the user to
!         the number of corrections used in the BFGS update. It
!         is not altered by the routine. Values of M less than 3 are
!         not recommended; large values of M will result in excessive
!         computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
! X       is a DOUBLE PRECISION array of length N. On initial entry
!         it must be set by the user to the values of the initial
!         estimate of the solution vector. On exit with IFLAG=0, it
!         contains the values of the variables at the best point
!         found (usually a solution).
! 
! F       is a DOUBLE PRECISION variable. Before initial entry and on
!         a re-entry with IFLAG=1, it must be set by the user to
!         contain the value of the function F at the point X.
! 
! G       is a DOUBLE PRECISION array of length N. Before initial
!         entry and on a re-entry with IFLAG=1, it must be set by
!         the user to contain the components of the gradient G at
!         the point X.
! 
! DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
!         user  wishes to provide the diagonal matrix Hk0 at each
!         iteration. Otherwise it should be set to .FALSE., in which
!         case  LBFGS will use a default value described below. If
!         DIAGCO is set to .TRUE. the routine will return at each
!         iteration of the algorithm with IFLAG=2, and the diagonal
!          matrix Hk0  must be provided in the array DIAG.
! 
! 
! DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
!         then on initial entry or on re-entry with IFLAG=2, DIAG
!         it must be set by the user to contain the values of the 
!         diagonal matrix Hk0.  Restriction: all elements of DIAG
!         must be positive.
! 
! IPRINT  is an INTEGER array of length two which must be set by the
!         user.
! 
!         IPRINT(1) specifies the frequency of the output:
!         IPRINT(1) < 0 : no output is generated,
!         IPRINT(1) = 0 : output only at first and last iteration,
!         IPRINT(1) > 0 : output every IPRINT(1) iterations.
! 
!         IPRINT(2) specifies the type of output generated:
!         IPRINT(2) = 0 : iteration count, number of function 
!                         evaluations, function value, norm of the
!                            gradient, and steplength,
!         IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
!                         variables and  gradient vector at the
!                         initial point,
!         IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
!                         variables,
!         IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
! 
! 
! EPS     is a positive DOUBLE PRECISION variable that must be set by
!         the user, and determines the accuracy with which the solution
!         is to be found. The subroutine terminates when
!
!                     ||G|| < EPS max(1,||X||),
!
!         where ||.|| denotes the Euclidean norm.
! 
! XTOL    is a  positive DOUBLE PRECISION variable that must be set by
!         the user to an estimate of the machine precision (e.g.
!         10**(-16) on a SUN station 3/60). The line search routine will
!         terminate if the relative width of the interval of uncertainty
!         is less than XTOL.
! 
! W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
!         workspace for LBFGS. This array must not be altered by the
!         user.
! 
! IFLAG   is an INTEGER variable that must be set to 0 on initial entry
!         to the subroutine. A return with IFLAG<0 indicates an error,
!         and IFLAG=0 indicates that the routine has terminated without
!         detecting errors. On a return with IFLAG=1, the user must
!         evaluate the function F and gradient G. On a return with
!         IFLAG=2, the user must provide the diagonal matrix Hk0.
! 
!         The following negative values of IFLAG, detecting an error,
!         are possible:
! 
!          IFLAG=-1  The line search routine MCSRCH failed. The
!                    parameter INFO provides more detailed information
!                    (see also the documentation of MCSRCH):
!
!                   INFO = 0  IMPROPER INPUT PARAMETERS.
!
!                   INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
!                             UNCERTAINTY IS AT MOST XTOL.
!
!                   INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
!                             REQUIRED AT THE PRESENT ITERATION.
!
!                   INFO = 4  THE STEP IS TOO SMALL.
!
!                   INFO = 5  THE STEP IS TOO LARGE.
!
!                   INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
!                             THERE MAY NOT BE A STEP WHICH SATISFIES
!                             THE SUFFICIENT DECREASE AND CURVATURE
!                             CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!
! 
!          IFLAG=-2  The i-th diagonal element of the diagonal inverse
!                    Hessian approximation, given in DIAG, is not
!                    positive.
!       
!          IFLAG=-3  Improper input parameters for LBFGS (N or M are
!                    not positive).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine BFGS_optimizer( N, X, G, F, EPS, XTOL, GTOL, STPMIN, STPMAX, IFLAG )

  !-----------------------------------------------------------------------------
  ! Variables that need to be initialized
  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: N
  real(dp), dimension(N)          :: X, G
  real(dp)          :: F
  real(dp), intent(in)                            :: EPS, XTOL
  real(dp), intent(in)                            :: GTOL
  real(dp), intent(in)                            :: STPMIN
  real(dp), intent(in)                            :: STPMAX
  integer, intent(out)                            :: IFLAG

  !-----------------------------------------------------------------------------
  integer, dimension(2)                           :: IPRINT
  integer                                         :: M
  integer                                         :: MP, LP
  logical                                         :: DIAGCO
  real(dp), allocatable, dimension(:)             :: DIAG
  real(dp), allocatable, dimension(:)             :: W
  real(dp)                                        :: GNORM, STP1, FTOL, STP
  real(dp)                                        :: YS, YY, SQ, YR
  real(dp)                                        :: BETA, XNORM

  integer                                         :: ITER, NFUN, POINT
  integer                                         :: ISPT, IYPT, MAXFEV, INFO
  integer                                         :: BOUND, NPT, CP, I, NFEV
  integer                                         :: INMC, IYCN, ISCN
  logical                                         :: FINISH

  integer                                         :: INFOC
  logical                                         :: BRACKT, STAGE1
  real(dp)                                        :: DG, DGM, DGINIT, DGTEST
  real(dp)                                        :: DGX, DGXM, DGY, DGYM
  real(dp)                                        :: FINIT, FTEST1
  real(dp)                                        :: FM, FX, FXM, FY, FYM, STX, ST
  real(dp)                                        :: STY,STMIN,STMAX,WIDTH,WIDTH1
  !-----------------------------------------------------------------------------

  M = 7

  LP = 6
  MP = 6
  DIAGCO = .false.

  IFLAG = 0

  IPRINT(1) = 1
  IPRINT(2) = 1

  allocate( DIAG(N) )
  allocate( W(N*(2*M+1)+2*M) )

  DIAG(:) = 0.0_dp
  W(:) = 0.0_dp

  call LBFGS(N,M,IPRINT,X,G,DIAG,W,F,EPS,XTOL,DIAGCO,MP,LP,GTOL,STPMIN,STPMAX,&
  GNORM,STP1,FTOL,STP,YS,YY,SQ,YR,BETA,XNORM,IFLAG,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,&
  INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,FINISH,INFOC,BRACKT,STAGE1,DG,DGM,DGINIT,&
  DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY,STMIN,STMAX,&
  WIDTH,WIDTH1)

  deallocate( diag, W )

end subroutine BFGS_optimizer

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine LBFGS(N,M,IPRINT,X,G,DIAG,W,F,EPS,XTOL,DIAGCO,MP,LP,GTOL,STPMIN,STPMAX,&
     GNORM,STP1,FTOL,STP,YS,YY,SQ,YR,BETA,XNORM,IFLAG,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,&
     INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,FINISH,INFOC,BRACKT,STAGE1,DG,DGM,DGINIT,&
     DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY,STMIN,STMAX,&
     WIDTH,WIDTH1)

  ! Variables that need to be initialized
  integer :: N,M,IPRINT(2)
  real(dp):: X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
  real(dp):: F,EPS,XTOL
  logical ::DIAGCO
  integer ::MP,LP
  real(dp):: GTOL,STPMIN,STPMAX

  ! Variables that should not be initialized
  real(dp):: GNORM,STP1,FTOL,STP,YS,YY,SQ,YR,BETA,XNORM
  integer :: IFLAG,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV
  integer :: INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
  logical :: FINISH

  !ss Variables for linear earch alg - noit initialized
  integer ::INFOC
  logical:: BRACKT,STAGE1
  real(dp)::DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY,STMIN,STMAX,WIDTH,WIDTH1
  !
  ! LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
  !                          JORGE NOCEDAL
  !                        *** July 1990 ***
  ! 
  !     This subroutine solves the unconstrained minimization problem
  ! 
  !                      min F(x),    x= (x1,x2,...,xN),
  !
  !      using the limited memory BFGS method. The routine is especially
  !      effective on problems involving a large number of variables. In
  !      a typical iteration of this method an approximation Hk to the
  !      inverse of the Hessian is obtained by applying M BFGS updates to
  !      a diagonal matrix Hk0, using information from the previous M steps.
  !      The user specifies the number M, which determines the amount of
  !      storage required by the routine. The user may also provide the
  !      diagonal matrices Hk0 if not satisfied with the default choice.
  !      The algorithm is described in "On the limited memory BFGS method
  !      for large scale optimization", by D. Liu and J. Nocedal,
  !      Mathematical Programming B 45 (1989) 503-528.
  ! 
  !      The user is required to calculate the function value F and its
  !      gradient G. In order to allow the user complete control over
  !      these computations, reverse  communication is used. The routine
  !      must be called repeatedly under the control of the parameter
  !      IFLAG. 
  !
  !      The steplength is determined at each iteration by means of the
  !      line search routine MCVSRCH, which is a slight modification of
  !      the routine CSRCH written by More and Thuente.
  ! 
  !      The calling statement is 
  ! 
  !          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
  ! 
  !      where
  ! 
  !     N       is an INTEGER variable that must be set by the user to the
  !             number of variables. It is not altered by the routine.
  !             Restriction: N>0.
  ! 
  !     M       is an INTEGER variable that must be set by the user to
  !             the number of corrections used in the BFGS update. It
  !             is not altered by the routine. Values of M less than 3 are
  !             not recommended; large values of M will result in excessive
  !             computing time. 3<= M <=7 is recommended. Restriction: M>0.
  ! 
  !     X       is a DOUBLE PRECISION array of length N. On initial entry
  !             it must be set by the user to the values of the initial
  !             estimate of the solution vector. On exit with IFLAG=0, it
  !             contains the values of the variables at the best point
  !             found (usually a solution).
  ! 
  !     F       is a DOUBLE PRECISION variable. Before initial entry and on
  !             a re-entry with IFLAG=1, it must be set by the user to
  !             contain the value of the function F at the point X.
  ! 
  !     G       is a DOUBLE PRECISION array of length N. Before initial
  !             entry and on a re-entry with IFLAG=1, it must be set by
  !             the user to contain the components of the gradient G at
  !             the point X.
  ! 
  !     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
  !             user  wishes to provide the diagonal matrix Hk0 at each
  !             iteration. Otherwise it should be set to .FALSE., in which
  !             case  LBFGS will use a default value described below. If
  !             DIAGCO is set to .TRUE. the routine will return at each
  !             iteration of the algorithm with IFLAG=2, and the diagonal
  !              matrix Hk0  must be provided in the array DIAG.
  ! 
  ! 
  !     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
  !             then on initial entry or on re-entry with IFLAG=2, DIAG
  !             it must be set by the user to contain the values of the 
  !             diagonal matrix Hk0.  Restriction: all elements of DIAG
  !             must be positive.
  ! 
  !     IPRINT  is an INTEGER array of length two which must be set by the
  !             user.
  ! 
  !             IPRINT(1) specifies the frequency of the output:
  !             IPRINT(1) < 0 : no output is generated,
  !             IPRINT(1) = 0 : output only at first and last iteration,
  !             IPRINT(1) > 0 : output every IPRINT(1) iterations.
  ! 
  !             IPRINT(2) specifies the type of output generated:
  !             IPRINT(2) = 0 : iteration count, number of function 
  !                             evaluations, function value, norm of the
  !                                gradient, and steplength,
  !             IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
  !                             variables and  gradient vector at the
  !                             initial point,
  !             IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
  !                             variables,
  !             IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
  ! 
  ! 
  !     EPS     is a positive DOUBLE PRECISION variable that must be set by
  !             the user, and determines the accuracy with which the solution
  !             is to be found. The subroutine terminates when
  !
  !                         ||G|| < EPS max(1,||X||),
  !
  !             where ||.|| denotes the Euclidean norm.
  ! 
  !     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
  !             the user to an estimate of the machine precision (e.g.
  !             10**(-16) on a SUN station 3/60). The line search routine will
  !             terminate if the relative width of the interval of uncertainty
  !             is less than XTOL.
  ! 
  !     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
  !             workspace for LBFGS. This array must not be altered by the
  !             user.
  ! 
  !     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
  !             to the subroutine. A return with IFLAG<0 indicates an error,
  !             and IFLAG=0 indicates that the routine has terminated without
  !             detecting errors. On a return with IFLAG=1, the user must
  !             evaluate the function F and gradient G. On a return with
  !             IFLAG=2, the user must provide the diagonal matrix Hk0.
  ! 
  !             The following negative values of IFLAG, detecting an error,
  !             are possible:
  ! 
  !              IFLAG=-1  The line search routine MCSRCH failed. The
  !                        parameter INFO provides more detailed information
  !                        (see also the documentation of MCSRCH):
  !
  !                       INFO = 0  IMPROPER INPUT PARAMETERS.
  !
  !                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
  !                                 UNCERTAINTY IS AT MOST XTOL.
  !
  !                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
  !                                 REQUIRED AT THE PRESENT ITERATION.
  !
  !                       INFO = 4  THE STEP IS TOO SMALL.
  !
  !                       INFO = 5  THE STEP IS TOO LARGE.
  !
  !                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
  !                                 THERE MAY NOT BE A STEP WHICH SATISFIES
  !                                 THE SUFFICIENT DECREASE AND CURVATURE
  !                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
  !
  ! 
  !              IFLAG=-2  The i-th diagonal element of the diagonal inverse
  !                        Hessian approximation, given in DIAG, is not
  !                        positive.
  !           
  !              IFLAG=-3  Improper input parameters for LBFGS (N or M are
  !                        not positive).
  !
  !    ON THE DRIVER:
  !
  !    The program that calls LBFGS must contain the declaration:
  !
  !                          EXTERNAL LB2
  !
  !    LB2 is a BLOCK DATA that defines the default values of several
  !    parameters described in the COMMON section. 
  ! 
  !    COMMON:
  ! 
  !    The subroutine contains one common area, which the user may wish to
  !    reference:
  ! 
  !     MP  is an INTEGER variable with default value 6. It is used as the
  !     unit number for the printing of the monitoring information
  !     controlled by IPRINT.
  ! 
  !     LP  is an INTEGER variable with default value 6. It is used as the
  !     unit number for the printing of error messages. This printing
  !     may be suppressed by setting LP to a non-positive value.
  ! 
  !     GTOL is a DOUBLE PRECISION variable with default value 0.9, which
  !     controls the accuracy of the line search routine MCSRCH. If the
  !     function and gradient evaluations are inexpensive with respect
  !     to the cost of the iteration (which is sometimes the case when
  !     solving very large problems) it may be advantageous to set GTOL
  !     to a small value. A typical small value is 0.1.  Restriction:
  !     GTOL should be greater than 1.E-04.
  ! 
  !     STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
  !     specify lower and uper bounds for the step in the line search.
  !     Their default values are 1.E-20 and 1.E+20, respectively. These
  !     values need not be modified unless the exponents are too large
  !     for the machine being used, or unless the problem is extremely
  !     badly scaled (in which case the exponents should be increased).
  !
  !  MACHINE DEPENDENCIES
  !
  !     The only variables that are machine-dependent are XTOL,
  !     STPMIN and STPMAX.
  !
  !  GENERAL INFORMATION
  ! 
  !    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
  ! 
  !    Input/Output  :  No input; diagnostic messages on unit MP and
  !                     error messages on unit LP.
  ! 
  ! 
  !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  real(dp):: DDOT
  real(dp)::ONE,ZERO

  ONE = 1.0_dp
  ZERO = 0.0_dp

  if (IFLAG == 0) then

     ITER= 0

     if ( N <= 0 .or. M <= 0 ) then
        write (0,*) 'N or M is less than 0. Must be positive'
        return
     end if

     if ( GTOL <= 1.E-04_dp ) then
        if ( LP > 0 ) then
           write (0,*) ' GTOL IS LESS THAN OR EQUAL TO 1.E-04'
           write (0,*) ' IT HAS BEEN RESET TO 9.E-01'
        end if
        GTOL = 9.E-01_dp
     end if

     NFUN = 1
     POINT = 0
     FINISH = .false.

     if ( DIAGCO ) then
        do I = 1, N
           if ( DIAG(I) <= ZERO ) then
              IFLAG = -2
              if ( LP > 0 ) then
                 write (0,*) 'IFLAG=-2, The ',N,'th diagonal element of the '
                 write (0,*) 'inverse Hessian approximation is not positive'
              end if
              return
           end if
        end do   
     else
        do I = 1, N
           DIAG(I)= 1.0_dp
        end do
     end if

     !!    THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
     !!      ---------------------------------------
     !!    THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
     !!      OTHER TEMPORARY INFORMATION.
     !!      LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
     !!      LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
     !!      IN THE FORMULA THAT COMPUTES H*G.
     !!      LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
     !!      STEPS.
     !!      LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
     !!      GRADIENT DIFFERENCES.
     !!  
     !!      THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
     !!      CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
     !!  
     ISPT= N+2*M
     IYPT= ISPT+N*M     
     do I=1,N
        W(ISPT+I)= -G(I)*DIAG(I)
     enddo
     GNORM= DSQRT(DDOTjre(N,G,1,G,1))
     STP1= ONE/GNORM

     !! PARAMETERS FOR LINE SEARCH ROUTINE
     FTOL= 1.0E-4_dp
     MAXFEV= 20

     if (IPRINT(1) >= 0) then
        call LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH,MP,LP,GTOL,STPMIN,STPMAX)
     end if
  end if


  !!  MAIN ITERATION LOOP
  mainwhile: do while (FINISH)

     if (IFLAG  ==  0) then
        ITER= ITER+1
        INFO=0
        BOUND=ITER-1

        if (ITER == 1) then

           if (ITER  >  M) BOUND=M

           YS= DDOTjre(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)

           if (.not.DIAGCO) then
              YY= DDOTjre(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
              do  I=1,N
                 DIAG(I)= YS/YY
              end do
           else
              IFLAG=2
              return
           end if

        end if

     end if

     if (IFLAG  ==  0 .or. IFLAG  ==  2) then
        !! For IFLAG = 0 or 2 - must recalculate Hk0
        !!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
        !!     "Updating quasi-Newton matrices with limited storage",
        !!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.

        if (DIAGCO) then
           do  I=1,N
              if (DIAG(I) <= ZERO) then
                 IFLAG=-2
                 if (LP >  0) then
                    write (0,*) 'IFLAG=-2, The ',N,'th diagonal element of the '
                    write (0,*) 'inverse Hessian approximation is not positive'
                 end if
                 return
              end if
           end do
        end if

        CP= POINT
        if (POINT == 0) CP=M
        W(N+CP)= ONE/YS
        do  I=1,N
           W(I)= -G(I)
        end do
        CP= POINT

        do I= 1,BOUND
           CP=CP-1
           if (CP ==  -1)CP=M-1
           SQ= DDOTjre(N,W(ISPT+CP*N+1),1,W,1)
           INMC=N+M+CP+1
           IYCN=IYPT+CP*N
           W(INMC)= W(N+CP+1)*SQ
           call DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
        end do

        do I=1,N
           W(I)=DIAG(I)*W(I)
        end do

        do I=1,BOUND
           YR= DDOTjre(N,W(IYPT+CP*N+1),1,W,1)
           BETA= W(N+CP+1)*YR
           INMC=N+M+CP+1
           BETA= W(INMC)-BETA
           ISCN=ISPT+CP*N
           call DAXPY(N,BETA,W(ISCN+1),1,W,1)
           CP=CP+1
           if (CP == M)CP=0
        end do

        !! STORE THE NEW SEARCH DIRECTION
        do  I=1,N
           W(ISPT+POINT*N+I)= W(I)
        end do

        !! OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
        !! BY USING THE LINE SEARCH ROUTINE MCSRCH

        NFEV=0
        STP=ONE
        if (ITER == 1) STP=STP1
        do I=1,N
           W(I)=G(I)
        end do

     end if

     !! . . DO for IFLAG=0,1 or 2.
     !!
     call MCSRCH (N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,XTOL,MAXFEV,INFO,&
          NFEV,DIAG,LP,GTOL,STPMIN,STPMAX,INFOC,BRACKT,STAGE1,&
          DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,&
          FY,FYM,STX,STY,STMIN,STMAX,WIDTH,WIDTH1)

     if (INFO  ==  -1) then
        IFLAG=1
        return
     end if

     if (INFO  /=  1) then
        IFLAG = -1
        if (LP > 0) then
           write (0,*) 'IFLAG= -1:  LINE SEARCH FAILED. SEE'
           write (0,*) 'DOCUMENTATION OF ROUTINE MCSRCH:  ERROR RETURN'
           write (0,*) 'OF LINE SEARCH: INFO= ',INFO
           write (0,*) 'POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT'
           write (0,*) 'OR INCORRECT TOLERANCES'
        end if
        return
     end if

     NFUN= NFUN + NFEV

     !! . . COMPUTE THE NEW STEP AND GRADIENT CHANGE 
     NPT=POINT*N

     do  I=1,N
        W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
        W(IYPT+NPT+I)= G(I)-W(I)
     end do

     POINT=POINT+1

     if (POINT == M) POINT=0

     !!  TERMINATION TEST
     GNORM= DSQRT(DDOTjre(N,G,1,G,1))
     XNORM= DSQRT(DDOTjre(N,X,1,X,1))
     XNORM= DMAX1(1.0_dp,XNORM)
     if (GNORM/XNORM  <=  EPS) FINISH=.true.

     if (IPRINT(1) >= 0) then
        call LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH,MP,LP,GTOL,STPMIN,STPMAX)
     end if

     if (FINISH) then
        IFLAG=0
        return
     end if

  end do mainwhile

end subroutine LBFGS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine LB1 (IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH,MP,LP,GTOL,STPMIN,STPMAX)
  !! THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
  !! AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
  integer :: IPRINT(2),ITER,NFUN,LP,MP,N,M, I
  real(dp) :: X(N),G(N),F,GNORM,STP,GTOL,STPMIN,STPMAX
  logical ::FINISH

  if (ITER == 0) then
     write (0,*) '*************************************************'
     write (0,*) 'N= ',N,'Number of Corrections=',ITER
     write (0,*) 'INITIAL VALUES:  F',F,'  GNORM= ',GNORM

     if (IPRINT(2) >= 1) then
        write (0,*) 'VECTOR X= '
        do I=1,N
           write (0,*) X(I)
        end do
        write (0,*) ' GRADIENT VECTOR G= '
        do I=1,N
           write (0,*) G(I)
        end do
     end if
     write (0,*) '*************************************************'
     write (0,*) '   I   NFUN',NFUN,'FUNC',F,'GNORM',GNORM,'STEPLENGTH',STP
  else
     if ((IPRINT(1) == 0).and.(ITER /= 1.and..not.FINISH)) return

     if (IPRINT(1) /= 0)then
        if (mod(ITER-1,IPRINT(1)) == 0.or.FINISH)then
           if (IPRINT(2) > 1.and.ITER > 1) then
              write (0,*) 'I NFUN',NFUN,'FUNC',F,'GNORM',GNORM,'STEPLENGTH',STP
              write (0,*) '   I   NFUN',NFUN,'FUNC',F,'GNORM',GNORM,'STEPLENGTH',STP
           end if
        else
           return
        end if
     else
        if ( IPRINT(2) > 1.and.FINISH)  then
           write (0,*) 'I NFUN',NFUN,'FUNC',F,'GNORM',GNORM,'STEPLENGTH',STP
           write (0,*) '   I   NFUN',NFUN,'FUNC',F,'GNORM',GNORM,'STEPLENGTH',STP
        end if
     end if

     if (IPRINT(2) == 2.or.IPRINT(2) == 3) then
        if (FINISH)then
           write (0,*) ' FINAL POINT X= '
        else
           write (0,*) 'VECTOR X= '
        end if

        do I=1,N
           write (0,*) X(I)
        end do

        if (IPRINT(2) == 3) then
           write (0,*) ' GRADIENT VECTOR G= '
           do I=1,N
              write (0,*) G(I)
           end do
        end if

     end if
     if ( FINISH ) write (*,*) 'THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS IFLAG = 0'

  end if

end subroutine LB1



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine daxpy(n,da,dx,incx,dy,incy)
  !! constant times a vector plus a vector.
  !! uses unrolled loops for increments equal to one.
  !! jack dongarra, linpack, 3/11/78.
  real(dp) ::dx(1),dy(1),da
  integer :: i,incx,incy,ix,iy,m,mp1,n

  if (n <= 0) return

  if (da  ==  0.0_dp) return

  if (incx  /=  1.or. incy  /=  1) then 
     !! code for unequal increments or equal increment not equal to 1
     ix = 1
     iy = 1
     if (incx < 0)ix = (-n+1)*incx + 1
     if (incy < 0)iy = (-n+1)*incy + 1
     do i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
     end do
     return
  end if

  m = mod(n,4)

  if ( m  /=  0 ) then 
     do i = 1,m
        dy(i) = dy(i) + da*dx(i)
     end do
     if ( n  <  4 ) return
  end if

  mp1 = m + 1

  do  i = mp1,n,4
     dy(i) = dy(i) + da*dx(i)
     dy(i + 1) = dy(i + 1) + da*dx(i + 1)
     dy(i + 2) = dy(i + 2) + da*dx(i + 2)
     dy(i + 3) = dy(i + 3) + da*dx(i + 3)
  end do

end subroutine daxpy



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
  !!   forms the dot product of two vectors.
  !!   uses unrolled loops for increments equal to one.
  !!   jack dongarra, linpack, 3/11/78.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

DoublePrecision function ddot( n,dx,incx,dy,incy )	!for some reason, this cannot link. possible conflict with tom778.f90. copied as ddotJre
  integer i,incx,incy,ix,iy,m,mp1,n
  DoublePrecision dx(n),dy(n),dtemp
  dtemp = 0.0_dp
  ddot = dtemp

  if (n <= 0) return

  if (incx .ne. 1.or. incy .ne. 1) then  ! /= means ".ne."
     !!   code for unequal increments or equal increments not equal to 1
     ix = 1
     iy = 1
     if (incx < 0)ix = (-n+1)*incx + 1
     if (incy < 0)iy = (-n+1)*incy + 1
     do  i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
     end do
     ddot = dtemp
     return
  end if

  m = mod(n,5)

  if ( m  .ne.  0 ) then
     do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
     end do
     if ( n  <  5 ) then 
        ddot = dtemp
        return
     end if
  end if

  mp1 = m + 1
  do i = mp1,n,5
     dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +  &
             dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do

  ddot = dtemp       

return
end function ddot

DoublePrecision Function DDOTjre(N,dx,incx,dy,incy)	!written by JRE 20200630. cf. https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-fortran/top/blas-and-sparse-blas-routines/blas-routines/blas-level-1-routines-and-functions/dot.html
integer i,j,incx,incy,ix,iy,m,mp1,n
DoublePrecision dx(N), dy(N), dtemp

  dtemp = 0.0_dp
  ddotjre = dtemp

  if (n <= 0) return

  if (incx .ne. 1.or. incy .ne. 1) then  ! /= means ".ne."
     !!   code for unequal increments or equal increments not equal to 1
     ix = 1
     iy = 1
     if (incx < 0)ix = (-n+1)*incx + 1
     if (incy < 0)iy = (-n+1)*incy + 1
     do  i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
     enddo
     ddotjre = dtemp
     return
  end if

  m = mod(n,5)

  if ( m  .ne.  0 ) then
     do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
     end do
     if ( n  <  5 ) then 
        ddotjre = dtemp
        return
     end if
  end if

  mp1 = m + 1
  do i = mp1,n,5
     dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +  &
             dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do

  ddotjre = dtemp       

return
end function DDOTjre



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!   LINE SEARCH ROUTINE MCSRCH
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine MCSRCH ( N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA,&
     LP,GTOL,STPMIN,STPMAX,INFOC,BRACKT,STAGE1,DG,DGM,DGINIT,DGTEST,&
     DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY,STMIN,&
     STMAX,WIDTH,WIDTH1)

  ! Variables that are init in calling routine
  integer :: N,MAXFEV,INFO,NFEV
  real(dp):: F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
  real(dp):: X(N),G(N),S(N),WA(N)
  integer                                         :: LP

  ! Variables that are set internally but saved
  integer :: INFOC
  logical :: BRACKT,STAGE1
  real(dp) :: DG,DGM,DGINIT,DGTEST,DGX,DGXM,STMIN,STMAX,WIDTH,WIDTH1
  real(dp) :: DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY

  !SUBROUTINE MCSRCH
  !
  !  A slight modification of the subroutine CSRCH of More' and Thuente.
  !    The changes are to allow reverse communication, and do not affect
  !    the performance of the routine. 
  !
  !    THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
  !    A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
  !    AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
  !    UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
  !    UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
  !    MINIMIZER OF THE MODIFIED FUNCTION
  !
  !      F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
  !
  !   IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
  !   HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
  !   THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
  !   CONTAINS A MINIMIZER OF F(X+STP*S).
  !
  !   THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
  !   THE SUFFICIENT DECREASE CONDITION
  !
  !   F(X+STP*S)  <=  F(X) + FTOL*STP*(GRADF(X)'S),
  !
  !   AND THE CURVATURE CONDITION
  !
  !  ABS(GRADF(X+STP*S)'S))  <=  GTOL*ABS(GRADF(X)'S).
  !
  !  IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
  !  IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
  !  BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
  !  CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
  !  ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
  !  SATISFIES THE SUFFICIENT DECREASE CONDITION.
  !
  !  THE SUBROUTINE STATEMENT IS
  !
  !  SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
  !     WHERE
  !
  !   N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !   OF VARIABLES.
  !
  !   X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
  !   BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
  !   X + STP*S.
  !
  !   F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
  !   AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
  !
  !   G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
  !   GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
  !   OF F AT X + STP*S.
  !
  !   S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
  !   SEARCH DIRECTION.
  !
  !   STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
  !   INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
  !   STP CONTAINS THE FINAL ESTIMATE.
  !
  !  FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
  !  communication implementation GTOL is defined in a COMMON
  !  statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
  !  CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
  !  SATISFIED.
  !
  !  XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
  !  WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
  !  IS AT MOST XTOL.
  !
  !  STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
  !  SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
  !  communication implementatin they are defined in a COMMON
  !  statement).
  !
  !  MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
  !  OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
  !  MAXFEV BY THE END OF AN ITERATION.
  !
  !  INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
  !
  !  INFO = 0  IMPROPER INPUT PARAMETERS.
  !
  !  INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
  !
  !  INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
  !  DIRECTIONAL DERIVATIVE CONDITION HOLD.
  !
  !  INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
  !  IS AT MOST XTOL.
  !
  !  INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
  !
  !  INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
  !
  !  INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
  !
  !  INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
  !  THERE MAY NOT BE A STEP WHICH SATISFIES THE
  !  SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
  !  TOLERANCES MAY BE TOO SMALL.
  !
  !  NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
  !  CALLS TO FCN.
  !
  !  WA IS A WORK ARRAY OF LENGTH N.
  !
  !  SUBPROGRAMS CALLED
  !
  !  MCSTEP
  !
  !  FORTRAN-SUPPLIED...ABS,MAX,MIN
  !
  !  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
  !  JORGE J. MORE, DAVID J. THUENTE

  integer :: J
  real(dp):: P5,P66,XTRAPF,ZERO
  P5=0.5D0
  P66=0.66D0
  XTRAPF=4.0_dp
  ZERO=0.0_dp

  if (INFO  /= -1 ) then 

     INFOC = 1

     !! CHECK THE INPUT PARAMETERS FOR ERRORS.
     if (N <= 0.or.STP <= ZERO.or.FTOL < ZERO.or.&
          GTOL < ZERO.or.XTOL < ZERO.or.STPMIN < ZERO &
          .or.STPMAX < STPMIN.or.MAXFEV  <=  0) return

     !! COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
     !! AND CHECK THAT S IS A DESCENT DIRECTION.
     DGINIT = ZERO
     do  J = 1, N
        DGINIT = DGINIT + G(J)*S(J)
     end do

     if (DGINIT  >=  ZERO) then
        if (LP  >=  0) then
           write (0,*) '  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION'
        end if
        return
     end if

     !! INITIALIZE LOCAL VARIABLES.

     BRACKT = .false.
     STAGE1 = .true.
     NFEV = 0
     FINIT = F
     DGTEST = FTOL*DGINIT
     WIDTH = STPMAX - STPMIN
     WIDTH1 = WIDTH/P5
     do  J = 1, N
        WA(J) = X(J)
     end do

     !!  THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
     !!  FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
     !!  THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
     !!  FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
     !!  THE INTERVAL OF UNCERTAINTY.
     !!  THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
     !!  FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
     STX = ZERO
     FX = FINIT
     DGX = DGINIT
     STY = ZERO
     FY = FINIT
     DGY = DGINIT

  end if

  !! START OF ITERATION.
  mainiter: do

     if (INFO  /=  -1) then

        !! SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
        !! TO THE PRESENT INTERVAL OF UNCERTAINTY.

        if (BRACKT) then
           STMIN = min(STX,STY)
           STMAX = max(STX,STY)
        else
           STMIN = STX
           STMAX = STP + XTRAPF*(STP - STX)
        end if

        !! FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
        STP = max(STP,STPMIN)
        STP = min(STP,STPMAX)

        !! IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
        !! STP BE THE LOWEST POINT OBTAINED SO FAR.
        if ( ( BRACKT .and. ( STP <= STMIN .or. STP >= STMAX ) ) &
             .or. NFEV >= MAXFEV-1 .or. INFOC == 0  &
             .or. ( BRACKT .and. STMAX-STMIN <= XTOL*STMAX ) ) STP = STX

        !! EVALUATE THE FUNCTION AND GRADIENT AT STP
        !! AND COMPUTE THE DIRECTIONAL DERIVATIVE.
        !! We return to main program to obtain F and G.
        do J = 1, N
           X(J) = WA(J) + STP*S(J)
        end do

        INFO=-1
        return

     end if

     INFO=0
     NFEV = NFEV + 1
     DG = ZERO

     do J = 1, N
        DG = DG + G(J)*S(J)
     end do

     FTEST1 = FINIT + STP*DGTEST

     !!   TEST FOR CONVERGENCE.
     if ((BRACKT .and. (STP  <=  STMIN .or. STP  >=  STMAX)) .or. INFOC  ==  0) INFO = 6
     if (STP  ==  STPMAX .and. F  <=  FTEST1 .and. DG  <=  DGTEST) INFO = 5
     if (STP  ==  STPMIN .and.(F  >  FTEST1 .or. DG  >=  DGTEST)) INFO = 4
     if (NFEV  >=  MAXFEV) INFO = 3
     if (BRACKT .and. STMAX-STMIN  <=  XTOL*STMAX) INFO = 2
     if (F  <=  FTEST1 .and. abs(DG)  <=  GTOL*(-DGINIT)) INFO = 1

     !! CHECK FOR TERMINATION.
     if (INFO  /=  0) return

     !! IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
     !! FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
     if (STAGE1 .and. F  <=  FTEST1 .and. DG  >=  min(FTOL,GTOL)*DGINIT) STAGE1 = .false.

     !!  A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
     !!  WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
     !!  FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
     !!  DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
     !!  OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
     if (STAGE1 .and. F  <=  FX .and. F  >  FTEST1) then

        !! DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
        FM = F - STP*DGTEST
        FXM = FX - STX*DGTEST
        FYM = FY - STY*DGTEST
        DGM = DG - DGTEST
        DGXM = DGX - DGTEST
        DGYM = DGY - DGTEST

        !! CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
        !! AND TO COMPUTE THE NEW STEP.
        call MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,STMAX,INFOC)

        !! RESET THE FUNCTION AND GRADIENT VALUES FOR F
        FX = FXM + STX*DGTEST
        FY = FYM + STY*DGTEST
        DGX = DGXM + DGTEST
        DGY = DGYM + DGTEST
     else

        !! CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
        !! AND TO COMPUTE THE NEW STEP.
        call MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,BRACKT,STMIN,STMAX,INFOC)
     end if

     !! FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE INTERVAL OF UNCERTAINTY.
     if (BRACKT) then
        if (abs(STY-STX)  >=  P66*WIDTH1) STP=STX+P5*(STY-STX)
        WIDTH1 = WIDTH
        WIDTH = abs(STY-STX)
     end if

  end do mainiter
end subroutine MCSRCH



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP_step,BRACKT,STPMIN,STPMAX,INFO)

  integer :: INFO
  real(dp) ::STX,FX,DX,STY,FY,DY,STP,FP,DP_step,STPMIN,STPMAX
  logical ::BRACKT,BOUND
  !  SUBROUTINE MCSTEP
  !
  !   THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
  !   A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
  !   A MINIMIZER OF THE FUNCTION.
  !
  !   THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
  !   VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
  !   ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
  !   DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
  !   MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
  !   WITH ENDPOINTS STX AND STY.
  !
  !     THE SUBROUTINE STATEMENT IS
  !
  !       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP_step,BRACKT,
  !                        STPMIN,STPMAX,INFO)
  !
  !    WHERE
  !
  !   STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
  !   THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
  !   SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
  !   OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
  !   SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
  !
  !   STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
  !   THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
  !   THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
  !   UPDATED APPROPRIATELY.
  !
  !   STP, FP, AND DP_step ARE VARIABLES WHICH SPECIFY THE STEP,
  !   THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
  !   IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
  !   BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
  !
  !   BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
  !   HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
  !   THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
  !   IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
  !
  !   STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
  !   AND UPPER BOUNDS FOR THE STEP.
  !
  !   INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
  !   IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
  !   ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
  !   INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
  !
  !   SUBPROGRAMS CALLED
  !
  !   FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
  !
  !   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
  !   JORGE J. MORE, DAVID J. THUENTE

  real(dp) :: GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
  INFO = 0

  !!  CHECK THE INPUT PARAMETERS FOR ERRORS.
  if ((BRACKT.and.(STP <= min(STX,STY).or.STP >= max(STX,STY))).or. &
       DX*(STP-STX) >= 0..or.STPMAX < STPMIN) return

  !! DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
  SGND = DP_step*(DX/abs(DX))

  if (FP  >  FX) then
     !!  FIRST CASE. A HIGHER FUNCTION VALUE.
     !!  THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
     !!  TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
     !!  ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
     INFO = 1
     BOUND = .true.
     THETA = 3.*(FX - FP)/(STP - STX) + DX + DP_step
     S = max(abs(THETA),abs(DX),abs(DP_step))
     GAMMA = S*sqrt((THETA/S)**2 - (DX/S)*(DP_step/S))
     if (STP  <  STX) GAMMA = -GAMMA
     P = (GAMMA - DX) + THETA
     Q = ((GAMMA - DX) + GAMMA) + DP_step
     R = P/Q
     STPC = STX + R*(STP - STX)
     STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
     if (abs(STPC-STX)  <  abs(STPQ-STX)) then
        STPF = STPC
     else
        STPF = STPC + (STPQ - STPC)/2
     end if
     BRACKT = .true.

  else if (SGND  <  0.0) then
     !  SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
     !  OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
     !  STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
     !  THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
     INFO = 2
     BOUND = .false.
     THETA = 3*(FX - FP)/(STP - STX) + DX + DP_step
     S = max(abs(THETA),abs(DX),abs(DP_step))
     GAMMA = S*sqrt((THETA/S)**2 - (DX/S)*(DP_step/S))
     if (STP  >  STX) GAMMA = -GAMMA
     P = (GAMMA - DP_step) + THETA
     Q = ((GAMMA - DP_step) + GAMMA) + DX
     R = P/Q
     STPC = STP + R*(STX - STP)
     STPQ = STP + (DP_step/(DP_step-DX))*(STX - STP)
     if (abs(STPC-STP)  >  abs(STPQ-STP)) then
        STPF = STPC
     else
        STPF = STPQ
     end if
     BRACKT = .true.

  else if (abs(DP_step)  <  abs(DX)) then
     ! THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
     ! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
     ! THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
     ! IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
     ! IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
     ! EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
     ! COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
     ! CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN
     INFO = 3
     BOUND = .true.
     THETA = 3.*(FX - FP)/(STP - STX) + DX + DP_step
     S = max(abs(THETA),abs(DX),abs(DP_step))

     ! THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
     ! TO INFINITY IN THE DIRECTION OF THE STEP.
     GAMMA = S*sqrt(max(0.0_dp,(THETA/S)**2 - (DX/S)*(DP_step/S)))
     if (STP  >  STX) GAMMA = -GAMMA
     P = (GAMMA - DP_step) + THETA
     Q = (GAMMA + (DX - DP_step)) + GAMMA
     R = P/Q
     if (R  <  0.0 .and. GAMMA  /=  0.0) then
        STPC = STP + R*(STX - STP)
     else if (STP  >  STX) then
        STPC = STPMAX
     else
        STPC = STPMIN
     end if
     STPQ = STP + (DP_step/(DP_step-DX))*(STX - STP)
     if (BRACKT) then
        if (abs(STP-STPC)  <  abs(STP-STPQ)) then
           STPF = STPC
        else
           STPF = STPQ
        end if
     else
        if (abs(STP-STPC)  >  abs(STP-STPQ)) then
           STPF = STPC
        else
           STPF = STPQ
        end if
     end if

  else
     ! FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
     ! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
     ! NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
     ! IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
     INFO = 4
     BOUND = .false.
     if (BRACKT) then
        THETA = 3*(FP - FY)/(STY - STP) + DY + DP_step
        S = max(abs(THETA),abs(DY),abs(DP_step))
        GAMMA = S*sqrt((THETA/S)**2 - (DY/S)*(DP_step/S))
        if (STP  >  STY) GAMMA = -GAMMA
        P = (GAMMA - DP_step) + THETA
        Q = ((GAMMA - DP_step) + GAMMA) + DY
        R = P/Q
        STPC = STP + R*(STY - STP)
        STPF = STPC
     else if (STP  >  STX) then
        STPF = STPMAX
     else
        STPF = STPMIN
     end if
  end if

  !  UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
  !  DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
  if (FP  >  FX) then
     STY = STP
     FY = FP
     DY = DP_step
  else
     if (SGND  <  0.0) then
        STY = STX
        FY = FX
        DY = DX
     end if
     STX = STP
     FX = FP
     DX = DP_step
  end if

  ! COMPUTE THE NEW STEP AND SAFEGUARD IT.
  STPF = min(STPMAX,STPF)
  STPF = max(STPMIN,STPF)
  STP = STPF
  if (BRACKT .and. BOUND) then
     if (STY  >  STX) then
        STP = min(STX+0.66*(STY-STX),STP)
     else
        STP = max(STX+0.66*(STY-STX),STP)
     end if
  end if
end subroutine MCSTEP



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine dump(N,M,IPRINT,X,G,DIAG,W,F,EPS,XTOL,DIAGCO,MP,LP,GTOL,&
     STPMIN,STPMAX,GNORM,STP1,FTOL,STP,YS,YY,SQ,YR,BETA,XNORM,IFLAG,ITER,&
     NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,&
     FINISH,INFOC,BRACKT,STAGE1,DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,&
     FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY,STMIN,STMAX,WIDTH,WIDTH1)

  ! Variables that need to be initialized
  integer ::N,M,IPRINT(2)
  real(dp):: X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
  real(dp):: F,EPS,XTOL,GTOL,STMIN,STMAX
  logical ::DIAGCO
  integer:: MP,LP

  ! Variables that should not be initialized
  real(dp):: GNORM,STP1,FTOL,STP,YS,YY,SQ,YR,BETA,XNORM
  integer ::IFLAG,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO
  integer :: BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
  logical ::FINISH

  ! Variables for linear earch alg - noit initialized
  integer ::INFOC
  logical ::BRACKT,STAGE1
  real(dp):: DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM
  real(dp):: FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY
  real(dp):: STPMIN,STPMAX,WIDTH,WIDTH1

  write (0,*)'***********x=   ************'
  ! write (0,*) x
  write (0,*)'***********g=   ************'
  ! write (0,*) g
  write (0,*)'***********dag= ************'
  ! write (0,*) diag
  write (0,*)'***********w=   ************'
  ! write (0,*) w
  write (0,*)'***********************'
  write (0,*) 'F',F
  write (0,*) 'EPS',EPS
  write (0,*) 'XTOL',XTOL
  write (0,*) 'DIAGCO',DIAGCO
  write (0,*) 'MP',MP
  write (0,*) 'LP',LP
  write (0,*) 'GTOL',GTOL
  write (0,*) 'STPMIN',STPMIN
  write (0,*) 'STPMAX',STPMAX
  write (0,*) 'GNORM',GNORM
  write (0,*) 'STP1',STP1
  write (0,*) 'FTOL',FTOL
  write (0,*) 'STP',STP
  write (0,*) 'YS',YS
  write (0,*) 'YY',YY
  write (0,*) 'SQ',SQ
  write (0,*) 'YR',YR
  write (0,*) 'BETA',BETA
  write (0,*) 'XNORM',XNORM
  write (0,*) 'IFLAG',IFLAG
  write (0,*) 'ITER',ITER
  write (0,*) 'NFUN',NFUN
  write (0,*) 'POINT',POINT
  write (0,*) 'ISPT',ISPT
  write (0,*) 'IYPT',IYPT
  write (0,*) 'MAXFEV',MAXFEV
  write (0,*) 'INFO',INFO
  write (0,*) 'BOUND',BOUND
  write (0,*) 'NPT',NPT
  write (0,*) 'CP',CP
  write (0,*) 'I',I
  write (0,*) 'NFEV',NFEV
  write (0,*) 'INMC',INMC
  write (0,*) 'IYCN',IYCN
  write (0,*) 'ISCN',ISCN
  write (0,*) 'FINISH',FINISH
  write (0,*) 'INFOC',INFOC
  write (0,*) 'BRACKT',BRACKT
  write (0,*) 'STAGE1',STAGE1
  write (0,*) 'DG',DG
  write (0,*) 'DGM',DGM
  write (0,*) 'DGINIT',DGINIT
  write (0,*) 'DGTEST',DGTEST
  write (0,*) 'DGX',DGX
  write (0,*) 'DGXM',DGXM
  write (0,*) 'DGY',DGY
  write (0,*) 'DGYM',DGYM
  write (0,*) 'FINIT',FINIT
  write (0,*) 'FTEST1',FTEST1
  write (0,*) 'FM',FM
  write (0,*) 'FX',FX
  write (0,*) 'FXM',FXM
  write (0,*) 'FY',FY
  write (0,*) 'FYM',FYM
  write (0,*) 'STX',STX
  write (0,*) 'STY',STY
  write (0,*) 'STPMIN',STMIN
  write (0,*) 'STPMAX',STMAX
  write (0,*) 'WIDTH',WIDTH
  write (0,*) 'WIDTH1',WIDTH1

end subroutine dump

end module moduel_lbfgs
