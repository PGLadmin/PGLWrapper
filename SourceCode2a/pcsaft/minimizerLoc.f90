
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module minimizer
!> \brief parameters and variables for  phase stability
!!
!! This module contains parameters and variables for a phase stability
!! analysis and flash calculations.
!! \todo variables
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module minimizer

  use PARAMETERS, only: dp
  use BASIC_VARIABLES, only: ncomp
  implicit none

  integer, parameter                              :: outp = 0

  private
  public :: cg_optimizer, bfgs, steepest_descent

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine cg
!
! conjugate gradient minimization
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine cg_optimizer ( ndim, x, f, xtol, gtol, ftol, maxiter, iprint, iflag, func ) !,grad )
  
   use PARAMETERS, only: dp
  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in out)       :: x
  real(dp), intent(in out)                        :: f
  real(dp), intent(in)                            :: xtol
  real(dp), intent(in)                            :: gtol
  real(dp), intent(in)                            :: ftol
  integer, intent(in)                             :: maxiter
  integer, intent(in)                             :: iprint
  integer, intent(in out)                         :: iflag
  !-----------------------------------------------------------------------------
  interface
     subroutine func (n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
!!$     subroutine grad (func,n,x,g)
!!$       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
!!$       integer, intent(in)                        :: n
!!$       real(dp), intent(in)                       :: x(n)
!!$       real(dp)                                   :: g(n)
!!$     end subroutine grad
  end interface
  !-----------------------------------------------------------------------------

  integer                                         :: iter
  real(dp)                                        :: alpha,fp,gnorm,gnormp
  real(dp), save, allocatable                     :: g(:),u(:)
  !-----------------------------------------------------------------------------

  if ( .not. allocated(g) ) allocate(g(ndim),u(ndim))

  iter = 0

  ! call func ( ndim, x, f )
  call grad ( func, ndim, x, f, g )

  gnorm = dot_product( g, g )
  ! g(1:ndim) = g(1:ndim)/sqrt(gnorm)
  gnorm = gnorm / ndim
  if ( iprint == 1 ) then
     write(6,'(a,i8,10es15.7)') ' iter,f,gnorm =', iter, f, gnorm
  else if ( iprint >= 2 ) then
     write(6,'(a,i8,10es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
  end if
  u(1:ndim) = - g(1:ndim)

  do iter = 1, maxiter

     fp = f

     !--------------------------------------------------------------------------
     ! line minimization
     !--------------------------------------------------------------------------
     call quad_interpolate (ndim,x,u,f,xtol,ftol,alpha,iprint,iflag,func)

     !--------------------------------------------------------------------------
     ! if quad interpolation failed, perform golden section
     !--------------------------------------------------------------------------
     if ( iflag/100 /= 0 ) then
        iflag = iflag -(iflag/100)*100
        call golden_section(ndim,x,u,f,xtol,ftol,alpha,iprint,iflag,func)
     end if
     if ( iflag/100 /= 0 ) return
     x(1:ndim) = x(1:ndim) +alpha*u(1:ndim)
     ! call func ( ndim, x, f )
     call grad ( func, ndim, x, f, g )

     !--------------------------------------------------------------------------
     ! store previous gnorm
     !--------------------------------------------------------------------------
     gnormp = gnorm
     gnorm = dot_product( g, g )
     ! g(1:ndim) = g(1:ndim) / sqrt(gnorm)
     gnorm = gnorm / ndim
     u(1:ndim) = -g(1:ndim) + gnorm / gnormp * u(1:ndim)
     if ( iprint == 1 ) then
        write(6,'(a,i8,10es15.7)') ' iter,f,gnorm =', iter, f, gnorm
     else if ( iprint >= 2 ) then
        write(6,'(a,i8,10es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
     end if

     !--------------------------------------------------------------------------
     ! check convergence 
     !--------------------------------------------------------------------------
!!$    if ( abs(alpha) < xtol ) then
!!$      print *,'>>> CG converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol =',alpha,xtol
!!$      iflag = iflag +1 
!!$      return
     if ( gnorm < gtol ) then
        print *,'>>> CG converged wrt gtol'
        write(6,'(a,2es15.7)') '   gnorm, gtol =', gnorm, gtol
        iflag = iflag + 2
        return
     else if ( abs( f - fp ) / abs(fp) < ftol ) then
        print *,'>>> CG converged wrt ftol'
        write(6,'(a,2es15.7)') '   f-fp,ftol =', abs( f - fp ) / abs(fp), ftol
        iflag = iflag + 3
        return
     end if
  end do

  print *,'*** maxiter exceeded ***'
  print *,'steepest_descent done.'
  iflag = iflag + 10

end subroutine cg_optimizer


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine bfgs
!
! Broyden-Fletcher-Goldfarb-Shanno type of Quasi-Newton method.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine bfgs ( ndim, x0, f, xtol, gtol, ftol, maxiter, iprint, iflag, func ) !,grad)

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in out)       :: x0
  real(dp), intent(in out)                        :: f
  real(dp), intent(in)                            :: xtol
  real(dp), intent(in)                            :: gtol
  real(dp), intent(in)                            :: ftol
  integer, intent(in)                             :: maxiter
  integer, intent(in)                             :: iprint
  integer, intent(in out)                         :: iflag
  !-----------------------------------------------------------------------------
  interface
     subroutine func(n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
!!$     subroutine grad(func,n,x,g)
!!$       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
!!$       integer, intent(in)                        :: n
!!$       real(dp), intent(in)                       :: x(n)
!!$       real(dp)                                   :: g(n)
!!$     end subroutine grad
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: i, j, iter
  real(dp), save, allocatable                     :: gg(:,:), aa(:,:), cc(:,:)
  real(dp), save, allocatable                     :: x(:), u(:), v(:), y(:), gp(:), ggy(:), ygg(:), g(:)
  real(dp)                                        :: tmp1, tmp2
  real(dp)                                        :: b
  real(dp)                                        :: svy,svyi
  real(dp)                                        :: fp
  real(dp)                                        :: f_out
  real(dp)                                        :: alpha
  real(dp)                                        :: gnorm
  !-----------------------------------------------------------------------------

  if ( .not.allocated(gg) ) then
     allocate( gg(ndim,ndim), aa(ndim,ndim), cc(ndim,ndim))
     allocate( x(ndim), u(ndim), v(ndim), y(ndim), g(ndim), gp(ndim), ggy(ndim), ygg(ndim))
  end if

  !-----------------------------------------------------------------------------
  ! initial G = I
  !-----------------------------------------------------------------------------
  gg(1:ndim,1:ndim) = 0.0_dp
  do i = 1 ,ndim
     gg(i,i) = 1.0_dp
  end do

  ! call func ( ndim, x0, f )
  call grad ( func, ndim, x0, f, g )
  gnorm = dot_product( g, g )
  ! g(1:ndim)= g(1:ndim) / sqrt(gnorm)
  x(1:ndim) = x0(1:ndim)

  iter = 0
  gnorm = gnorm / ndim
  if ( iprint == 1 ) then
     write (6,'(a,i8,100es15.7)') ' iter, f, gnorm =', iter, f, gnorm
  else if ( iprint >= 2 ) then
     write (6,'(a,i8,100es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
  end if

  !-----------------------------------------------------------------------------
  ! iteration
  !-----------------------------------------------------------------------------
  do iter = 1, maxiter

     u(1:ndim) = 0.0_dp
     do i = 1, ndim
        u(1:ndim) = u(1:ndim) - gg(1:ndim,i) * g(i)
     end do

     !--------------------------------------------------------------------------
     ! store previous func and grad values
     !--------------------------------------------------------------------------
     fp = f
     gp(1:ndim) = g(1:ndim)

     !--------------------------------------------------------------------------
     ! line minimization
     !--------------------------------------------------------------------------
     call quad_interpolate ( ndim, x, u, f, xtol, ftol, alpha, iprint, iflag, func )

     !--------------------------------------------------------------------------
     ! if quad interpolation failed, perform golden section
     !--------------------------------------------------------------------------
     if ( iflag/100 /= 0 ) then
        iflag = iflag - ( iflag / 100 ) * 100
        call golden_section ( ndim, x, u, f, xtol, ftol, alpha, iprint, iflag, func )
     end if
     if ( iflag/100 /= 0 ) then
        x0(1:ndim) = x(1:ndim)
        return
     end if

     x(1:ndim) = x(1:ndim) + alpha * u(1:ndim)
     call grad ( func, ndim, x, f_out, g )
     gnorm = dot_product( g, g )
     ! g(1:ndim) = g(1:ndim)/sqrt(gnorm)
     gnorm = gnorm / ndim
     if ( iprint == 1 ) then
        write(6,'(a,i8,100es15.7)') ' iter,f,gnorm =', iter, f, gnorm
     else if ( iprint >= 2 ) then
        write(6,'(a,i8,100es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
     end if

     !--------------------------------------------------------------------------
     ! check convergence
     !--------------------------------------------------------------------------
!!$    if ( abs(alpha) < xtol ) then
!!$      print *,'>>> BFGS converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol =',alpha,xtol
!!$      x0(1:ndim) = x(1:ndim)
!!$      iflag = iflag +1
!!$      return
     if ( gnorm < gtol ) then
        print *,'>>> BFGS converged wrt gtol'
        write(6,'(a,2es15.7)') '   gnorm, gtol =', gnorm, gtol
        x0(1:ndim) = x(1:ndim)
        iflag = iflag +2
        return
     else if ( abs(f-fp)/abs(fp) < ftol ) then
        print *,'>>> BFGS converged wrt ftol'
        write(6,'(a,2es15.7)') '   f-fp,ftol =', abs( f - fp ) / abs(fp), ftol
        x0(1:ndim) = x(1:ndim)
        iflag = iflag + 3
        return
     end if

     v(1:ndim) = alpha * u(1:ndim)
     y(1:ndim) = g(1:ndim) - gp(1:ndim)

     !--------------------------------------------------------------------------
     ! update G matrix, gg, according to BFGS
     !--------------------------------------------------------------------------
     svy = dot_product( v, y )
     svyi = 1.0_dp / svy
     do i = 1, ndim
        tmp1 = 0.0_dp
        tmp2 = 0.0_dp
        do j =1,ndim
           aa(j,i) = v(j) * v(i) * svyi
           tmp1 = tmp1 + gg(i,j) * y(j)
           tmp2 = tmp2 + y(j) * gg(j,i)
        end do
        ggy(i) = tmp1
        ygg(i) = tmp2
     end do
     b = 1.0_dp
     do i = 1, ndim
        b = b + y(i) * ggy(i) * svyi
     end do
     aa(1:ndim,1:ndim) = aa(1:ndim,1:ndim) * b
     cc(1:ndim,1:ndim) = 0.0_dp
     do j =1,ndim
        do i =1,ndim
           cc(i,j) = cc(i,j) + ( v(i)*ygg(j) + ggy(i)*v(j) ) * svyi
        end do
     end do
     gg(1:ndim,1:ndim) = gg(1:ndim,1:ndim) + aa(1:ndim,1:ndim) - cc(1:ndim,1:ndim)
  end do

  print *,'maxiter exceeded in bfgs'
  iflag = iflag + 10
  x0(1:ndim) = x(1:ndim)

end subroutine bfgs



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine get_bracket ( ndim, x0, d, a, b, c, fa, fb, fc, iprint, iflag, func )

  !-----------------------------------------------------------------------------
  integer, intent(in):: ndim,iprint
  integer, intent(in out):: iflag
  real(dp), intent(in):: x0(ndim),d(ndim)
  real(dp), intent(in out):: a,b,fa,fb
  real(dp), intent(out):: c,fc
  !-----------------------------------------------------------------------------
  interface
     subroutine func(n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
  end interface

  !-----------------------------------------------------------------------------
  real(dp), parameter:: RATIO = 1.61803398875_dp
  real(dp), parameter:: RATIOI= 1.0_dp/RATIO
  real(dp), parameter:: MAXITER = 50
  ! real(dp):: dum,r,q,u,ulim,fu
  integer:: iter
  !-----------------------------------------------------------------------------

  call func (ndim,x0+a*d,fa)
  call func (ndim,x0+b*d,fb)

  iter = 0

10 continue

  iter = iter +1

  if ( iter > MAXITER ) then
     print *,'[Error] iter > MAXITER in get_bracket'
     print *,'  Search direction may not be a descent direction.'
     iflag = iflag + 1000
     return
  end if

  if ( abs(b-a) < 1.E-12_dp) then
     print *,'[Error] a and b is too close in get_bracket'
     print *,'  Search direction may not be a descent direction.'
     iflag = iflag + 2000
     return
  end if

  if ( fa < fb ) then
     c = a + RATIOI * ( b - a )
     call func ( ndim, x0+c*d, fc )
     call exchange ( c, b )
     call exchange ( fc, fb )
     if ( iprint == 3 ) then
        write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc =',a,b,c,fa,fb,fc
     end if
     goto 10
  else
     c = a + RATIO * ( b - a )
     call func ( ndim, x0+c*d, fc )
     if ( fb > fc ) then
        b = a + RATIO * ( c - a )
        call func ( ndim, x0+b*d, fb )
        call exchange ( b, c )
        call exchange ( fb, fc )
        if ( iprint == 3 ) then
           write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc =',a,b,c,fa,fb,fc
        end if
        goto 10
     end if
!!$      write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc =',a,b,c,fa,fb,fc
  end if

end subroutine get_bracket



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine exchange ( a, b )

  real(dp), intent(in out)                        :: a, b
  real(dp)                                        :: tmp
  tmp = a
  a = b
  b = tmp

end subroutine exchange


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine quad_interpolate ( ndim, x0, g, f, xtol, ftol, a, iprint, iflag, func )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in)           :: x0
  real(dp), dimension(ndim), intent(in)           :: g
  real(dp), intent(out)                           :: f
  real(dp), intent(in)                            :: xtol
  real(dp), intent(in)                            :: ftol
  real(dp), intent(out)                           :: a
  integer, intent(in)                             :: iprint
  integer, intent(in out)                         :: iflag

  !-----------------------------------------------------------------------------
  interface
     subroutine func(n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
  end interface
  !-----------------------------------------------------------------------------

  real(dp), parameter                             :: STP0    = 1.E-2_dp
  real(dp), parameter                             :: STPMAX  = 1.E+1_dp
  real(dp), parameter                             :: TINY    = 1.E-15_dp
  integer, parameter                              :: MAXITER = 20

  integer                                         :: iter, imin, iminJre(1), imax, ix
  real(dp)                                        :: r, q, fmin, fmax, dmin, dmax, d, xmin
  real(dp), save, allocatable                     :: xi(:), fi(:)
!!$  real(dp), dimension(ndim)                       :: g_JG
!!$  real(dp), dimension(ndim)                       :: max_step_size
!!$  real(dp)                                        :: factor
  !-----------------------------------------------------------------------------

  if ( .not. allocated(xi) ) allocate(xi(4),fi(4))

  !-----------------------------------------------------------------------------
  ! JG: limit maximum step size. Test this functionality
  !-----------------------------------------------------------------------------
!!$  max_step_size(:) = 5.0_dp
!!$  factor = maxval( abs(g(:)) / max_step_size(:) )
!!$
!!$  ! g(:) = sign( min( abs(g(:)), max_step_size(:) ), g(:) )
!!$  print *,' factor =', factor
!!$  if ( factor < 1.0_dp ) factor = 1.0_dp
!!$  g_JG(:) = g(:) / factor
!!$
!!$  print *,' g =',g(:)
!!$  read (*,*)


  xi(1) = 0.0_dp
  xi(2) = xi(1) + STP0
  call get_bracket ( ndim,x0,g,xi(1),xi(2),xi(3),fi(1),fi(2),fi(3),iprint,iflag,func )

  if ( iflag / 1000 /= 0 ) return

!!$    fi(1) = func(ndim,x0+xi(1)*g)
!!$    fi(2) = func(ndim,x0+xi(2)*g)
!!$    if ( fi(1) > fi(2) ) then
!!$      xi(3) = xi(1) +2*STP0
!!$      fi(3) = func(ndim,x0+xi(3)*g)
!!$    else
!!$      xi(3) = xi(1) -STP0
!!$      fi(3) = func(ndim,x0+xi(3)*g)
!!$    end if

  iter = 0

10 continue

  iter = iter + 1

  if ( iter > MAXITER ) then
     print *,' [Error] iter > MAXITER in quad_interpolate !!!'
     print *,'   iter,MAXITER = ', iter, MAXITER
     iflag = iflag + 100
     return
  end if

  !-----------------------------------------------------------------------------
  ! step 3; compute turning point
  !-----------------------------------------------------------------------------
  r = ( xi(2)-xi(1) ) * ( fi(2)-fi(3) )
  q = ( xi(2)-xi(3) ) * ( fi(2)-fi(1) )
  xi(4) = xi(2) - ( (xi(2)-xi(3) ) * q -(xi(2)-xi(1)) * r) / (2.0_dp*sign(max(abs(q-r),TINY),q-r))
  call func ( ndim, x0+xi(4)*g, fi(4) )
!!$    write(6,'(a,2(2x,4f11.2))') ' xi,fi =',xi(1:4),fi(1:4)

  !-----------------------------------------------------------------------------
  ! step4
  !-----------------------------------------------------------------------------
  fmin = min( fi(1), fi(2), fi(3) )
  fmax = max( fi(1), fi(2), fi(3) )
  dmin = min( abs(xi(4)-xi(1)), abs(xi(4)-xi(2)), abs(xi(4)-xi(3)) )

  if ( fi(4) < fmin .and. dmin > STPMAX ) then
     imax = 0
     dmax = 0.0_dp
     do ix = 1, 3
        d = abs( xi(4) - xi(ix) )
        if ( dmax < d ) then
           dmax = d
           imax = ix
        end if
     end do

     !--------------------------------------------------------------------------
     ! eliminate max function point and set point at STPMAX
     !--------------------------------------------------------------------------
     imax = 0
     fmax = - 1.E+30_dp
     do ix  = 1, 3
        if ( fmax < fi(ix) ) then
           fmax = fi(ix)
           imax = ix
        end if
     end do
     do ix = imax+1, 3
        xi(ix-1) = xi(ix)
        fi(ix-1) = fi(ix)
     end do
     if ( fi(2) > fi(1) ) then
        xi(3) = xi(1) + STPMAX
     else
        xi(3) = xi(2) + STPMAX
     end if
     call func ( ndim, x0+xi(3)*g, fi(3) )
     goto 10

  else if ( fi(4) > fmax ) then ! fi(4) is maximum

     !imin = minloc(  abs( xi(4) - xi(1:3) ), 1  )	!JRE fails data type for argument MASK
     iminJre = minloc(  abs( xi(4) - xi(1:3) )  )
     dmin = abs( xi(4) - xi(iminJre(1) ) )
     xmin = xi( iminJre(1) )

     !--------------------------------------------------------------------------
     ! eliminate nearest point to xi(4) and add a new point
     !--------------------------------------------------------------------------
     do ix = iminJre(1)+1, 3
        xi(ix-1) = xi(ix)
        fi(ix-1) = fi(ix)
     end do
!!$      if ( fi(2) > fi(1) ) then
!!$        xi(3) = xi(1) +STPMAX
!!$      else
!!$        xi(3) = xi(2) +STPMAX
!!$      end if
     xi(3) = ( xmin + xi(4) ) * 0.5_dp
     call func ( ndim, x0+xi(3)*g, fi(3) )
     goto 10
  end if

  !-----------------------------------------------------------------------------
  ! step 5: check convergence
  !-----------------------------------------------------------------------------
  if ( dmin < xtol ) then
     imin = 0
     fmin = 1.E+30_dp
     do ix = 1, 4
        if ( fmin > fi(ix) ) then
           imin = ix
           fmin = fi(ix)
        end if
     end do
     f = fi(imin)
     a = xi(imin)
     return
  end if

  !-----------------------------------------------------------------------------
  ! step 6: discard point of highest f value and replace it by xi(4)
  !-----------------------------------------------------------------------------
  imax = 0
  fmax = - 1.E+30_dp
  do ix = 1, 3
     if ( fmax < fi(ix) ) then
        imax = ix
        fmax = fi(ix)
     end if
  end do
  xi(imax) = xi(4)
  fi(imax) = fi(4)
  goto 10

end subroutine quad_interpolate


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine golden_section ( ndim,x0,g,f,xtol,ftol,alpha,iprint,iflag,func )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in)           :: x0
  real(dp), dimension(ndim), intent(in)           :: g
  real(dp), intent(out)                           :: f
  real(dp), intent(in)                            :: xtol
  real(dp), intent(in)                            :: ftol
  real(dp), intent(out)                           :: alpha
  integer, intent(in)                             :: iprint
  integer, intent(in out)                         :: iflag

  !-----------------------------------------------------------------------------
  interface
     subroutine func(n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
  end interface
  !-----------------------------------------------------------------------------

  real(dp), parameter:: STP0 = 1.E-1_dp
  real(dp), parameter:: GR   = 0.61803398875_dp
  real(dp), parameter:: GR2  = 1.0_dp -GR
  integer, parameter:: MAXITER = 100

  integer:: iter
  real(dp):: a,b1,b2,c,fa,fb1,fb2,fc,xl
  !-----------------------------------------------------------------------------

  a = 0.0_dp
  b1 = STP0
  call get_bracket ( ndim,x0,g,a,b1,c,fa,fb1,fc,iprint,iflag,func )
  if ( iflag / 1000 /= 0 ) return
  xl = ( c - a )
  b1 = a +GR2*xl
  b2 = a +GR *xl
  call func (ndim,x0+b1*g,fb1)
  call func (ndim,x0+b2*g,fb2)

  iter = 0
10 continue
  iter = iter + 1
  if ( iter > MAXITER ) then
     print *,'[Error] iter > MAXITER in golden_section.'
     print *,'  iter,MAXITER = ', iter, MAXITER
     iflag = iflag + 100
     return
!!$      stop
  end if
!!$    write(6,'(a,2(2x,4es11.3))') ' a,b1,b2,c,fa,fb1,fb2,fc =',a,b1,b2,c,fa,fb1,fb2,fc
  if ( fb1 > fb2 ) then
     a = b1
     fa = fb1
     b1 = b2
     fb1 = fb2
     xl = ( c - a )
     b2 = a + GR * xl
     call func (ndim,x0+b2*g,fb2)
  else
     c = b2
     fc = fb2
     b2 = b1
     fb2 = fb1
     xl = ( c - a )
     b1 = a + GR2 * xl
     call func (ndim,x0+b1*g,fb1)
  end if
!!$    print *,' xl,c,a,xtol =',xl,c,a,xtol
  if ( xl < xtol ) then
     if ( fb1 > fb2 ) then
        alpha = b2
        f = fb2
     else
        alpha = b1
        f = fb1
     end if
     return
  end if
  goto 10

end subroutine golden_section

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE grad
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine grad ( func, n, optpara, f, g )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: n
  real(dp), dimension(n), intent(in)              :: optpara
  real(dp), intent(in out)                        :: f
  real(dp), dimension(n), intent(in out)          :: g
  !-----------------------------------------------------------------------------
  interface
     subroutine func (n, x, f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
  end interface

  !-----------------------------------------------------------------------------
  integer                                         :: i
  real(dp)                                        :: delta, fmin
  real(dp), dimension(n)                          :: optpara_mod, fi
  !-----------------------------------------------------------------------------


  delta = 1.E-5_dp

  optpara_mod = optpara

  do i = 1, n

     optpara_mod = optpara
     optpara_mod(i) = optpara(i) * ( 1.0_dp + delta )

     call func ( n, optpara_mod, fmin )
     fi(i) = fmin

  enddo

  call func ( n, optpara, fmin )
  f = fmin


  do i = 1, n

     g(i) = ( fi(i) - f ) / ( optpara(i) * delta )
     ! write (*,*) fi(i), f, optpara(i)*delta

  enddo

end subroutine grad


!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$subroutine test_linmin
!!$
!!$  real(dp)                                        :: f
!!$  real(dp)                                        :: x0(10)
!!$  ! real(dp)                                        :: x0(10),g(10),a,b,c,fa,fb,fc
!!$  real(dp), external                              :: f1,f2
!!$  real(dp), parameter                             :: xtol = 1.E-3_dp
!!$  real(dp), parameter                             :: gtol = 1.E-5_dp
!!$  real(dp), parameter                             :: ftol = 1.E-4_dp
!!$  integer, parameter                              :: iprint = 2
!!$  integer:: iflag
!!$
!!$  x0(1) = -2.0_dp
!!$  x0(2) =  1.0_dp
!!$  call cg (2,x0,f,xtol,gtol,ftol,2000,iprint,iflag,rosenbrock,drosen)
!!$  ! call bfgs (2,x0,f,xtol,gtol,ftol,2000,iprint,iflag,rosenbrock,drosen)
!!$  print *,'iflag,f =',iflag,f
!!$end subroutine test_linmin
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function f1 (n,x)
!!$
!!$  integer, intent(in)                             :: n
!!$  real(dp), intent(in)                            :: x(n)
!!$  real(dp)                                        :: f1
!!$
!!$  f1 = x(1)*x(1) + 2.0_dp * exp(-x(1))
!!$
!!$end function f1
!!$
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function f2 (n,x)
!!$
!!$  integer, intent(in)                             :: n
!!$  real(dp), intent(in)                            :: x(n)
!!$  real(dp)                                        :: f2
!!$
!!$  f2 = - x(1) * cos(x(1))
!!$
!!$end function f2
!!$
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function rosenbrock (n,x)
!!$
!!$  integer, intent(in)                             :: n
!!$  real(dp), intent(in)                            :: x(n)
!!$  real(dp)                                        :: rosenbrock
!!$
!!$  rosenbrock = 100.0_dp*(x(2)-x(1)*x(1))**2 +(1.0_dp-x(1))*(1.0_dp-x(1))
!!$
!!$end function rosenbrock
!!$
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function drosen (n,x)
!!$
!!$  integer, intent(in)                             :: n
!!$  real(dp), intent(in)                            :: x(n)
!!$  real(dp)                                        :: drosen(n)
!!$
!!$  drosen(1) = -400.0_dp *x(1) *(x(2)-x(1)*x(1)) -2.0_dp*(1.0_dp-x(1))
!!$  drosen(2) = 200.0_dp *(x(2)-x(1)*x(1))
!!$
!!$end function drosen
!!$
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function f3 (n,x)
!!$
!!$  integer, intent(in)                             :: n
!!$  real(dp), intent(in)                            :: x(n)
!!$  real(dp)                                        :: f3
!!$  f3 = (x(1)+2*x(2)+7.0_dp)**2 +(2*x(1)+x(2)-5)**2
!!$
!!$end function f3
!!$
!!$
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$!
!!$!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!$
!!$function df3 (n,x)
!!$
!!$  integer, intent(in):: n
!!$  real(dp), intent(in):: x(n)
!!$  real(dp):: df3(n)
!!$  df3(1) = 2.0_dp*(x(1)+2*x(2)+7.0_dp) +4.0_dp*(2*x(1)+x(2)-5)
!!$  df3(2) = 4.0_dp*(x(1)+2*x(2)+7.0_dp) +2.0_dp*(2*x(1)+x(2)-5)
!!$
!!$end function df3



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine steepest_descent ( ndim, x, f, xtol, gtol, ftol, maxiter, iprint, iflag, func ) !, grad )

  !-----------------------------------------------------------------------------
  integer, intent(in)                             :: ndim
  real(dp), dimension(ndim), intent(in out)       :: x
  real(dp), intent(in out)                        :: f
  ! real(dp), intent(in)                            :: step_reduction    ! reduce step size. Use step_reduction=1.0_dp for no reduction. Otherwise positive value, lower than unity.
  real(dp), intent(in)                            :: xtol
  real(dp), intent(in)                            :: gtol
  real(dp), intent(in)                            :: ftol
  integer, intent(in)                             :: maxiter
  integer, intent(in)                             :: iprint
  integer, intent(in out)                         :: iflag
  !-----------------------------------------------------------------------------
  interface
     subroutine func (n,x,f)
       integer, parameter                         :: dp = selected_real_kind(15, 307)  !< double precision ( quad: (33, 4931) )
       integer, intent(in)                        :: n
       real(dp), dimension(n), intent(in)         :: x
       real(dp), intent(in out)                   :: f
     end subroutine func
  end interface
  !-----------------------------------------------------------------------------

  integer                                         :: iter
  real(dp)                                        :: alpha,fp,gnorm
  real(dp), save, allocatable                     :: g(:)
  !-----------------------------------------------------------------------------

  if ( .not. allocated(g) ) allocate(g(ndim))

  iter = 0

  ! call func ( ndim, x, f )
  call grad ( func, ndim, x, f, g )

  gnorm = dot_product( g, g )
  g(1:ndim) = - g(1:ndim) / sqrt(gnorm)
  gnorm = gnorm / ndim

  if ( iprint == 1 ) then
     write(6,'(a,i8,10es15.7)') ' iter, f, gnorm =', iter, f, gnorm
  else if ( iprint >= 2 ) then
     write(6,'(a,i8,10es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
  end if

  do iter = 1, maxiter

     fp = f

     !--------------------------------------------------------------------------
     ! line minimization
     !--------------------------------------------------------------------------
     call quad_interpolate ( ndim, x, g, f, xtol, ftol, alpha, iprint, iflag, func )

     !--------------------------------------------------------------------------
     ! if quad interpolation failed, perform golden section
     !--------------------------------------------------------------------------
     if ( iflag / 100 /= 0 ) then
        iflag = iflag - ( iflag / 100 ) * 100
        call golden_section ( ndim, x, g, f, xtol, ftol, alpha, iprint, iflag, func )
     end if

     if ( iflag / 100 /= 0 ) return

     x(1:ndim) = x(1:ndim) + alpha * g(1:ndim)

     ! call func( ndim, x, f )
     call grad( func, ndim, x, f, g )

     gnorm = dot_product( g, g )
     g(1:ndim) = - g(1:ndim) / sqrt(gnorm)
     gnorm = gnorm / ndim

     if ( iprint == 1 ) then
        write(6,'(a,i8,10es15.7)') ' iter, f, gnorm =', iter, f, gnorm
     else if ( iprint >= 2 ) then
        write(6,'(a,i8,10es15.7)') ' iter, x, f, gnorm =', iter, x(1:ndim), f, gnorm
     end if

     !--------------------------------------------------------------------------
     ! check convergence
     !--------------------------------------------------------------------------
!!$    if ( abs(alpha) < xtol ) then
!!$      print *,'>>> SD converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol =',alpha,xtol
!!$      iflag = iflag +1
!!$      return
     if ( gnorm < gtol ) then
        print *,'>>> SD converged wrt gtol'
        write(6,'(a,2es15.7)') '   gnorm, gtol =', gnorm, gtol
        iflag = iflag + 2
        return
     else if ( abs( f - fp ) / abs(fp) < ftol ) then
        print *,'>>> Sd converged wrt ftol'
        write(6,'(a,2es15.7)') '   f-fp, ftol =', abs( f - fp ) / abs(fp),ftol
        iflag = iflag + 3
        return
     end if
  end do

  print *,'maxiter exceeded in steepest_descent'
  print *,'steepest_descent done.'
  iflag = iflag + 10

end subroutine steepest_descent

end module minimizer
