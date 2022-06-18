!
!     MSIMSLC.F90 - Declare common MS IMSL MATH/LIBRARY and 
!                   MS IMSL STAT/LIBRARY routines
!
!     The modules defined by this file (all named msimslc? or msimsl?c)
!     are intended for use by msimslmd, msimslms, msimslsd, and msimslss
!     only: you should NOT use them directly in your programs.  
!     Instead, use the other modules named above, or use the 
!     module msimsl.
!
!     Copyright (c) 1994 by Visual Numerics, Inc.  All Rights Reserved.
!

module msimslc		! MS IMSL common routines

!
!     Level 1 BLAS
!

      interface
        subroutine iset (n, ia, ix, incx)
          integer    n, ia, incx, ix(*)
        end subroutine
      end interface

      interface
        subroutine icopy (n, ix, incx, iy, incy)
          integer    n, incx, incy, ix(*), iy(*)
        end subroutine
      end interface

      interface
        subroutine iadd (n, ia, ix, incx)
          integer    n, ia, incx, ix(*)
        end subroutine
      end interface

      interface
        subroutine isub (n, ia, ix, incx)
          integer    n, ia, incx, ix(*)
        end subroutine
      end interface

      interface
        subroutine iswap (n, ix, incx, iy, incy)
          integer    n, incx, incy, ix(*), iy(*)
        end subroutine
      end interface

      interface
        double precision function dsdot (n, sx, incx, sy,incy)
          integer    n, incx, incy
          real       sx(*), sy(*)
        end function
      end interface

      interface
        integer function isum (n, ix, incx)
          integer    n, incx, ix(*)
        end function
      end interface

      interface
        integer function iimin (n, ix, incx)
          integer    n, incx, ix(*)
        end function
      end interface

      interface
        integer function iimax (n, ix, incx)
          integer    n, incx, ix(*)
        end function
      end interface

!
!     Chapter 10:  Utilities
!

      interface
        subroutine wrirn (title, nrmat, ncmat, mat, ldmat,itring)
          integer    nrmat, ncmat, ldmat, itring, mat(ldmat,*)
          character  title*(*)
        end subroutine
      end interface

      interface
        subroutine wrirl (title, nrmat, ncmat, mat, ldmat,itring, &
			  fmt, rlabel, clabel)
          integer    nrmat, ncmat, ldmat, itring, mat(ldmat,*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*)
        end subroutine
      end interface

      interface
        subroutine wropt (iopt, iset, iscope)
          integer    iopt, iset, iscope
        end subroutine
      end interface

      interface
        subroutine pgopt (iopt, ipage)
          integer    iopt, ipage
        end subroutine
      end interface

      interface
        subroutine svign (n, ia, ib)
          integer    n, ia(*), ib(*)
        end subroutine
      end interface

      interface
        subroutine svigp (n, ia, ib, iperm)
          integer    n, ia(*), ib(*), iperm(*)
        end subroutine
      end interface

      interface
        subroutine isrch (n, ivalue, ix, incx, index)
          integer    n, ivalue, incx, index, ix(*)
        end subroutine
      end interface

      interface
        subroutine ssrch (n, string, chx, incx, index)
          integer    n, incx, index
          character  string*(*), chx(*)*(*)
        end subroutine
      end interface

      interface
        character *1 function achar (i)
          integer    i
        end function
      end interface

      interface
        integer function iachar (ch)
          character  ch
        end function
      end interface

      interface
        integer function icase (ch)
          character  ch
        end function
      end interface

      interface
        integer function iicsr (str1, str2)
          character  str1*(*), str2*(*)
        end function
      end interface

      interface
        integer function iidex (chrstr, key)
          character  chrstr*(*), key*(*)
        end function
      end interface

      interface
        subroutine cvtsi (string, number)
          integer    number
          character  string*(*)
        end subroutine
      end interface

      interface
        real function cpsec ()
        end function
      end interface

      interface
        subroutine timdy (ihour, minute, isec)
          integer    ihour, minute, isec
        end subroutine
      end interface

      interface
        subroutine tdate (iday, month, iyear)
          integer    iday, month, iyear
        end subroutine
      end interface

      interface
        integer function ndays (iday, imonth, iyear)
          integer    iday, imonth, iyear
        end function
      end interface

      interface
        subroutine ndyin (ndays, iday, month, iyear)
          integer    ndays, iday, month, iyear
        end subroutine
      end interface

      interface
        integer function idywk (iday, imonth, iyear)
          integer    iday, imonth, iyear
        end function
      end interface

      interface
        subroutine rnget (iseed)
          integer    iseed
        end subroutine
      end interface

      interface
        subroutine rnset (iseed)
          integer    iseed
        end subroutine
      end interface

      interface
        subroutine rnopt (iopt)
          integer    iopt
        end subroutine
      end interface

!
!     Reference Material
!

      interface
        subroutine erset (iersvr, ipact, isact)
          integer    iersvr, ipact, isact
        end subroutine
      end interface

      interface
        integer function iercd ()
        end function
      end interface

      interface
        integer function n1rty (iopt)
          integer    iopt
        end function
      end interface

      interface
        subroutine iwkcin (nelmts, nalc)
          integer    nelmts, nalc
        end subroutine
      end interface

      interface
        subroutine iwkin (nsu)
          integer    nsu
        end subroutine
      end interface

      interface
        integer function imach (n)
          integer    n
        end function
      end interface

      interface
        subroutine umach (n, nunit)
          integer    n, nunit
        end subroutine
      end interface

end module msimslc

module msimslcs		! MS IMSL common single-precision routines

      interface
        subroutine lslrg (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          real       a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine l2lrg (n, a, lda, b, ipath, x, fac, ipvt,wk)
          integer    n, lda, ipath, ipvt(*)
          real       a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine lftrg (n, a, lda, fac, ldfac, ipvt)
          integer    n, lda, ldfac, ipvt(*)
          real       a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine l2trg (n, a, lda, fac, ldfac, ipvt, scale)
          integer    n, lda, ldfac, ipvt(*)
          real       a(lda,*), fac(ldfac,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine lfsrg (n, fac, ldfac, ipvt, b, ipath, x)
          integer    n, ldfac, ipath, ipvt(*)
          real       fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine lfdrg (n, fac, ldfac, ipvt, det1, det2)
          integer    n, ldfac, ipvt(*)
          real       det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine linrg (n, a, lda, ainv, ldainv)
          integer    n, lda, ldainv
          real       a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine l2nrg (n, a, lda, ainv, ldainv, wk, iwk)
          integer    n, lda, ldainv, iwk(*)
          real       a(lda,*), ainv(ldainv,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine lsqrr (nra, nca, a, lda, b, tol, x, res,kbasis)
          integer    nra, nca, lda, kbasis
          real       tol, a(*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine l2qrr (nra, nca, a, lda, b, tol, x, res,kbasis, qr, &
			  qraux, ipvt, work)
          integer    nra, nca, lda, kbasis, ipvt(*)
          real       tol, a(lda,*), b(*), x(*), res(*), qr(*), qraux(*),   &
                     work(*)
        end subroutine
      end interface

      interface
        subroutine lsbrr (nra, nca, a, lda, b, tol, x, res,kbasis)
          integer    nra, nca, lda, kbasis
          real       tol, a(lda,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine l2brr (nra, nca, a, lda, b, tol, x, res,kbasis, qr,  &
			  brrux, ipvt, work)
          integer    nra, nca, lda, kbasis, ipvt(*)
          real       tol, a(lda,*), b(*), x(*), res(*), qr(*), &
		     brrux(*), work(*)
        end subroutine
      end interface

      interface
        subroutine lsvrr (nra, nca, a, lda, ipath, tol,irank, s, u, &
		          ldu, v, ldv)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          real       tol, a(lda,*), s(*), u(ldu,*), v(ldv,*)
        end subroutine
      end interface

      interface
        subroutine l2vrr (nra, nca, a, lda, ipath, tol,irank, s, u, ldu, &
			  v, ldv, wka, wk)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          real       tol, a(lda,*), s(*), u(*), v(*), wka(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine lsgrr (nra, nca, a, lda, tol, irank,ginva, ldginv)
          integer    nra, nca, lda, irank, ldginv
          real       tol, a(lda,*), ginva(ldginv,*)
        end subroutine
      end interface

      interface
        subroutine l2grr (nra, nca, a, lda, tol, irank,ginva, ldginv, wka, wk)
          integer    nra, nca, lda, irank, ldginv
          real       tol, a(lda,*), ginva(ldginv,*), wka(*), wk(*)
        end subroutine
      end interface

!

      interface
        subroutine evcsf (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          real       a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine e5csf (n, a, lda, eval, evec, ldevec, wk,iwk)
          integer    n, lda, ldevec, iwk(n)
          real       a(lda,*), eval(*), evec(ldevec,*), wk(n,3)
        end subroutine
      end interface

      interface
        real function episf (n, neval, a, lda, eval, evec,ldevec)
          integer    n, neval, lda, ldevec
          real       a(lda,*), eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        real function e2isf (n, neval, a, lda, eval, evec,ldevec, work)
          integer    n, neval, lda, ldevec
          real       a(lda,*), eval(*), evec(ldevec,*), work(*)
        end function
      end interface

      interface
        subroutine gvlsp (n, a, lda, b, ldb, eval)
          integer    n, lda, ldb
          real       a(lda,*), b(ldb,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine g3lsp (n, a, lda, b, ldb, eval, iwk, wk,s)
          integer    n, lda, ldb, iwk(*)
          real       a(lda,*), b(ldb,*), eval(*), wk(n,*), s(n+1,*)
        end subroutine
      end interface

!

      interface
        subroutine rline (nobs, xdata, ydata, b0, b1, stat)
          integer    nobs
          real       b0, b1, xdata(*), ydata(*), stat(12)
        end subroutine
      end interface

      interface
        subroutine rcurv (nobs, xdata, ydata, ndeg, b,sspoly, stat)
          integer    nobs, ndeg
          real       xdata(*), ydata(*), b(*), sspoly(*), stat(10)
        end subroutine
      end interface

      interface
        subroutine r2urv (nobs, x, y, ndeg, b, sspoly, stat,wk, iwk)
          integer    nobs, ndeg, iwk(*)
          real       x(*), y(*), b(*), sspoly(*), stat(10), wk(*)
        end subroutine
      end interface

!

      interface
        subroutine qdags (f, a, b, errabs, errrel, result,errest)
          real       f, a, b, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine q2ags (f, a, b, errabs, errrel, result, errest, maxsub, &
			  neval, nsubin, alist, blist, rlist,elist, iord)
          integer    maxsub, neval, nsubin, iord(*)
          real       f, a, b, errabs, errrel, result, errest, alist(*),    &
                     blist(*), rlist(*), elist(*)
          external   f
        end subroutine
      end interface

!

      interface
        subroutine fftrf (n, seq, coef)
          integer    n
          real       seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine f2trf (n, seq, coef, wfftr)
          integer    n
          real       seq(*), coef(*), wfftr(*)
        end subroutine
      end interface

      interface
        subroutine fftrb (n, coef, seq)
          integer    n
          real       coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine f2trb (n, coef, seq, wfftr)
          integer    n
          real       coef(*), seq(*), wfftr(*)
        end subroutine
      end interface

      interface
        subroutine fftri (n, wfftr)
          integer    n
          real       wfftr(*)
        end subroutine
      end interface

!

      interface
        subroutine uminf (fcn, n, xguess, xscale, fscale, iparam, &
			  rparam, x, fvalue)
          integer    n, iparam(*)
          real       fscale, fvalue, xguess(*), xscale(*), rparam(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine u2inf (fcn, n, xguess, xscale, fscale,iparam, rparam, &
		          x, fvalue, wk)
          integer    n, iparam(*)
          real       fscale, fvalue, xguess(*), xscale(*), rparam(*),    &
                     x(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine u4inf (iparam, rparam)
          integer    iparam(*)
          real       rparam(*)
        end subroutine
      end interface

!

      interface
        subroutine trnrr (nra, nca, a, lda, nrb, ncb, b, ldb)
          integer    nra, nca, lda, nrb, ncb, ldb
          real       a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine mxtxf (nra, nca, a, lda, nb, b, ldb)
          integer    nra, nca, lda, nb, ldb
          real       a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine mxtyf (nra, nca, a, lda, nrb, ncb, b, ldb,nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          real       a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine mrrrr (nra, nca, a, lda, nrb, ncb, b, ldb,nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          real       a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        real function blinf (nra, nca, a, lda, x, y)
          integer    nra, nca, lda
          real       a(*), x(*), y(*)
        end function
      end interface

!
!     Level 1 BLAS
!

      interface
        subroutine sset (n, sa, sx, incx)
          integer    n, incx
          real       sa, sx(*)
        end subroutine
      end interface

      interface
        subroutine cset (n, ca, cx, incx)
          integer    n, incx
          complex    ca, cx(*)
        end subroutine
      end interface

      interface
        subroutine scopy (n, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sx(*), sy(*)
        end subroutine
      end interface

      interface
        subroutine ccopy (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end subroutine
      end interface

      interface
        subroutine sscal (n, sa, sx, incx)
          integer    n, incx
          real       sa, sx(*)
        end subroutine
      end interface

      interface
        subroutine cscal (n, ca, cx, incx)
          integer    n, incx
          complex    ca, cx(*)
        end subroutine
      end interface

      interface
        subroutine csscal (n, sa, cx, incx)
          integer    n, incx
          real       sa
          complex    cx(*)
        end subroutine
      end interface

      interface
        subroutine svcal (n, sa, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sa, sx(*), sy(*)
        end subroutine
      end interface

      interface
        subroutine cvcal (n, ca, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    ca, cx(*), cy(*)
        end subroutine
      end interface

      interface
        subroutine csvcal (n, sa, cx, incx, cy, incy)
          integer    n, incx, incy
          real       sa
          complex    cx(*), cy(*)
        end subroutine
      end interface

      interface
        subroutine sadd (n, sa, sx, incx)
          integer    n, incx
          real       sa, sx(*)
        end subroutine
      end interface

      interface
        subroutine cadd (n, ca, cx, incx)
          integer    n, incx
          complex    ca, cx(*)
        end subroutine
      end interface

      interface
        subroutine ssub (n, sa, sx, incx)
          integer    n, incx
          real       sa, sx(*)
        end subroutine
      end interface

      interface
        subroutine csub (n, ca, cx, incx)
          integer    n, incx
          complex    ca, cx(*)
        end subroutine
      end interface

      interface
        subroutine saxpy (n, sa, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sa, sx(*), sy(*)
        end subroutine
      end interface

      interface
        subroutine caxpy (n, ca, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    ca, cx(*), cy(*)
        end subroutine
      end interface

      interface
        subroutine sswap (n, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sx(*), sy(*)
        end subroutine
      end interface

      interface
        subroutine cswap (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end subroutine
      end interface

      interface
        real function sdot (n, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sx(*), sy(*)
        end function
      end interface

      interface
        complex function cdotu (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end function
      end interface

      interface
        complex function cdotc (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end function
      end interface

      interface
        complex function czdotc (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end function
      end interface

      interface
        complex function czdotu (n, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    cx(*), cy(*)
        end function
      end interface

      interface
        real function sdsdot (n, sb, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sb, sx(*), sy(*)
        end function
      end interface

      interface
        complex function czcdot (n, ca, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    ca, cx(*), cy(*)
        end function
      end interface

      interface
        complex function czudot (n, ca, cx, incx, cy, incy)
          integer    n, incx, incy
          complex    ca, cx(*), cy(*)
        end function
      end interface

      interface
        real function sddoti (n, sa, acc, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sa, sx(*), sy(*)
          double precision acc(2)
        end function
      end interface

      interface
        real function sddota (n, sa, acc, sx, incx, sy, incy)
          integer    n, incx, incy
          real       sa, sx(*), sy(*)
          double precision acc(2)
        end function
      end interface

      interface
        complex function czdoti (n, ca, acc, cx, incx, cy,incy)
          integer    n, incx, incy
          double precision acc(4)
          complex    ca, cx(*), cy(*)
        end function
      end interface

      interface
        complex function czdota (n, ca, acc, cx, incx, cy,incy)
          integer    n, incx, incy
          double precision acc(4)
          complex    ca, cx(*), cy(*)
        end function
      end interface

      interface
        subroutine shprod (n, sx, incx, sy, incy, sz, incz)
          integer    n, incx, incy, incz
          real       sx(*), sy(*), sz(*)
        end subroutine
      end interface

      interface
        real function sxyz (n, sx, incx, sy, incy, sz, incz)
          integer    n, incx, incy, incz
          real       sx(*), sy(*), sz(*)
        end function
      end interface

      interface
        real function ssum (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        real function sasum (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        real function scasum (n, cx, incx)
          integer    n, incx
          complex    cx(*)
        end function
      end interface

      interface
        real function snrm2 (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        real function scnrm2 (n, cx, incx)
          integer    n, incx
          complex    cx(*)
        end function
      end interface

      interface
        real function sprdct (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        integer function ismin (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        integer function ismax (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        integer function isamin (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        integer function icamin (n, cx, incx)
          integer    n, incx
          complex    cx(*)
        end function
      end interface

      interface
        integer function isamax (n, sx, incx)
          integer    n, incx
          real       sx(*)
        end function
      end interface

      interface
        integer function icamax (n, cx, incx)
          integer    n, incx
          complex    cx(*)
        end function
      end interface

      interface
        subroutine srotg (sa, sb, sc, ss)
          real       sa, sb, sc, ss
        end subroutine
      end interface

      interface
        subroutine srot (n, sx, incx, sy, incy, c, s)
          integer    n, incx, incy
          real       c, s, sx(*), sy(*)
        end subroutine
      end interface

      interface
        subroutine csrot (n, cx, incx, cy, incy, c, s)
          integer    n, incx, incy
          real       c, s
          complex    cx(*), cy(*)
        end subroutine
      end interface

      interface
        subroutine srotmg (sd1, sd2, sx1, sy1, sparam)
          real       sd1, sd2, sx1, sy1, sparam(5)
        end subroutine
      end interface

      interface
        subroutine srotm (n, sx, incx, sy, incy, sparam)
          integer    n, incx, incy
          real       sx(*), sy(*), sparam(5)
        end subroutine
      end interface

      interface
        subroutine csrotm (n, cx, incx, cy, incy, sparam)
          integer    n, incx, incy
          real       sparam(5)
          complex    cx(*), cy(*)
        end subroutine
      end interface

!
!     Level 2 BLAS
!

      interface
        subroutine sgemv (trans, m, n, alpha, a, lda, x,incx, beta, y, incy)
          integer    m, n, lda, incx, incy
          real       alpha, beta, a(*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine cgemv (trans, m, n, alpha, a, lda, x,incx, beta, y, incy)
          integer    m, n, lda, incx, incy
          complex    alpha, beta, a(*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine sgbmv (trans, m, n, nlca, nuca, alpha, a, lda, x, &
			  incx, beta, y, incy)
          integer    m, n, nlca, nuca, lda, incx, incy
          real       alpha, beta, a(lda,*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine cgbmv (trans, m, n, nlca, nuca, alpha, a,lda, x, &
			  incx, beta, y, incy)
          integer    m, n, nlca, nuca, lda, incx, incy
          complex    alpha, beta, a(lda,*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine chemv (uplo, n, alpha, a, lda, x, incx,beta, y, incy)
          integer    n, lda, incx, incy
          complex    alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine chbmv (uplo, n, ncoda, alpha, a, lda, x,incx, beta, y, incy)
          integer    n, ncoda, lda, incx, incy
          complex    alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine ssymv (uplo, n, alpha, a, lda, x, incx,beta, y, incy)
          integer    n, lda, incx, incy
          real       alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine ssbmv (uplo, n, ncoda, alpha, a, lda, x,incx, beta, y, incy)
          integer    ncoda, n, lda, incx, incy
          real       alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine strmv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          real       a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctrmv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          complex    a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine stbmv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          real       a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctbmv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          complex    a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine strsv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          real       a(*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctrsv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          complex    a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine stbsv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          real       a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctbsv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          complex    a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine sger (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          real       alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine cgeru (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          complex    alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine cgerc (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          complex    alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine cher (uplo, n, alpha, x, incx, a, lda)
          integer    n, incx, lda
          real       alpha
          complex    x(*), a(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine cher2 (uplo, n, alpha, x, incx, y, incy,a, lda)
          integer    n, incx, incy, lda
          complex    alpha, x(*), y(*), a(lda,*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine ssyr (uplo, n, alpha, x, incx, a, lda)
          integer    n, incx, lda
          real       alpha, x(*), a(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine ssyr2 (uplo, n, alpha, x, incx, y, incy,a, lda)
          integer    n, incx, incy, lda
          real       alpha, x(*), y(*), a(lda,*)
          character  uplo*1
        end subroutine
      end interface

!
!     Level 3 BLAS
!

      interface
        subroutine sgemm (transa, transb, m, n, k, alpha, a,lda, b, &
			  ldb, beta, c, ldc)
          integer    m, n, k, lda, ldb, ldc
          real       alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  transa*1, transb*1
        end subroutine
      end interface

      interface
        subroutine cgemm (transa, transb, m, n, k, alpha, a,lda, b, &
			  ldb, beta, c, ldc)
          integer    m, n, k, lda, ldb, ldc
          complex    alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  transa*1, transb*1
        end subroutine
      end interface

      interface
        subroutine ssymm (side, uplo, m, n, alpha, a, lda, b, &
			  ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          real       alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine csymm (side, uplo, m, n, alpha, a, lda, b, &
			  ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          complex    alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine chemm (side, uplo, m, n, alpha, a, lda, b, &
			  ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          complex    alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine ssyrk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          real       alpha, beta, a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine csyrk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          complex    alpha, beta, a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine cherk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          real       alpha, beta
          complex    a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine ssyr2k (uplo, trans, n, k, alpha, a, lda,b, &
			   ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          real       alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine csyr2k (uplo, trans, n, k, alpha, a, lda,b, &
			   ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          complex    alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine cher2k (uplo, trans, n, k, alpha, a, lda,b, &
			   ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          real       beta
          complex    alpha, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine strmm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          real       alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctrmm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          complex    alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine strsm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          real       alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ctrsm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          complex    alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

!

      interface
        subroutine wrrrn (title, nra, nca, a, lda, itring)
          integer    nra, nca, lda, itring
          real       a(*)
          character  title*(*)
        end subroutine
      end interface

      interface
        subroutine wrrrl (title, nra, nca, a, lda, itring,fmt, rlabel, clabel)
          integer    nra, nca, lda, itring
          real       a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*)
        end subroutine
      end interface

      interface
        subroutine w2rrl (title, nra, nca, a, lda, itring,fmt, &
			  rlabel, clabel, chwk)
          integer    nra, nca, lda, itring
          real       a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*),  &
                     chwk(*)*10
        end subroutine
      end interface

      interface
        subroutine permu (n, x, ipermu, ipath, xpermu)
          integer    n, ipath, ipermu(*)
          real       x(*), xpermu(*)
        end subroutine
      end interface

      interface
        subroutine perma (nra, nca, a, lda, ipermu, ipath,aper, ldaper)
          integer    nra, nca, lda, ipath, ldaper, ipermu(*)
          real       a(lda,*), aper(ldaper,*)
        end subroutine
      end interface

      interface
        subroutine p2rma (nra, nca, a, lda, ipermu, ipath,aper, ldaper, work)
          integer    nra, nca, lda, ipath, ldaper, ipermu(*)
          real       a(lda,*), aper(ldaper,*), work(*)
        end subroutine
      end interface

      interface
        subroutine svrgn (n, ra, rb)
          integer    n
          real       ra(*), rb(*)
        end subroutine
      end interface

      interface
        subroutine svrgp (n, ra, rb, iperm)
          integer    n, iperm(*)
          real       ra(*), rb(*)
        end subroutine
      end interface

      interface
        subroutine srch (n, value, x, incx, index)
          integer    n, incx, index
          real       value, x(*)
        end subroutine
      end interface

      interface
        real function rnunf ()
        end function
      end interface

      interface
        subroutine rnun (nr, r)
          integer    nr
          real       r(*)
        end subroutine
      end interface


      interface
        subroutine plotp (ndata, nfun, x, a, lda, inc, range, &
			  symbol, xtitle, ytitle, title)
          integer    ndata, nfun, lda, inc
          real       x(*), a(lda,*), range(4)
          character  symbol*(*), xtitle*(*), ytitle*(*), title*(*)
        end subroutine
      end interface

!

      interface
        real function binom (n, m)
          integer    n, m
        end function
      end interface

      interface
        real function gamma (x)
          real       x
        end function
      end interface

      interface
        subroutine r9gaml (xmin, xmax)
          real       xmin, xmax
        end subroutine
      end interface

      interface
        real function beta (a, b)
          real       a, b
        end function
      end interface

!
!
!     Special Functions Chapter 11:  Probability Distribution Functions
!                                    and Inverses
!

      interface
        real function bindf (k, n, p)
          integer    k, n
          real       p
        end function
      end interface

      interface
        real function binpr (k, n, p)
          integer    k, n
          real       p
        end function
      end interface

      interface
        real function hypdf (k, n, m, l)
          integer    k, n, m, l
        end function
      end interface

      interface
        real function hyppr (k, n, m, l)
          integer    k, n, m, l
        end function
      end interface

      interface
        real function poidf (k, theta)
          integer    k
          real       theta
        end function
      end interface

      interface
        real function poipr (k, theta)
          integer    k
          real       theta
        end function
      end interface

      interface
        real function aks1df (nobs, d)
          integer    nobs
          real       d
        end function
      end interface

      interface
        real function ak21df (nobs, d, wk)
          integer    nobs
          real       d, wk(*)
        end function
      end interface

      interface
        real function aks2df (nobsx, nobsy, d)
          integer    nobsx, nobsy
          real       d
        end function
      end interface

      interface
        real function ak22df (nobsx, nobsy, d, wk)
          integer    nobsx, nobsy
          real       d, wk(*)
        end function
      end interface

      interface
        real function anordf (x)
          real       x
        end function
      end interface

      interface
        real function anorin (p)
          real       p
        end function
      end interface

      interface
        real function betdf (x, pin, qin)
          real       x, pin, qin
        end function
      end interface

      interface
        real function betin (p, pin, qin)
          real       p, pin, qin
        end function
      end interface

      interface
        real function bnrdf (x, y, rho)
          real       x, y, rho
        end function
      end interface

      interface
        real function chidf (chsq, df)
          real       chsq, df
        end function
      end interface

      interface
        real function chiin (p, df)
          real       p, df
        end function
      end interface

      interface
        real function csndf (chsq, df, alam)
          real       chsq, df, alam
        end function
      end interface

      interface
        real function fdf (f, dfn, dfd)
          real       f, dfn, dfd
        end function
      end interface

      interface
        real function fin (p, dfn, dfd)
          real       p, dfn, dfd
        end function
      end interface

      interface
        real function gamdf (x, a)
          real       x, a
        end function
      end interface

      interface
        real function tdf (t, df)
          real       t, df
        end function
      end interface

      interface
        real function tin (p, df)
          real       p, df
        end function
      end interface

      interface
        real function tndf (t, idf, delta)
          integer    idf
          real       t, delta
        end function
      end interface

      interface
        real function gcdf (x0, iopt, m, x, f)
          integer    iopt, m
          real       x0, x(*), f(*)
        end function
      end interface

      interface
        real function g4df (x0, iopt, m, x, f, wk, iwk)
          integer    iopt, m, iwk(*)
          real       x0, x(*), f(*), wk(*)
        end function
      end interface

      interface
        real function gcin (p, iopt, m, x, f)
          integer    iopt, m
          real       p, x(*), f(*)
        end function
      end interface

      interface
        real function g3in (p, iopt, m, x, f, wk, iwk)
          integer    iopt, m, iwk(*)
          real       p, x(*), f(*), wk(*)
        end function
      end interface

      interface
        real function amach (n)
          integer    n
        end function
      end interface

      interface
        logical function ifnan(x)
          real       x
        end function
      end interface

      interface
        real function g2df (x0, iopt, m, x, f, wk)
          integer    iopt, m
          real       x0, x(*), f(*), wk(*)
        end function
      end interface

      interface
        real function g3df (xx0, iopt, m, xx, f, c)
          integer    iopt, m
          real       xx0, xx(*), f(*), c(*)
        end function
      end interface

      interface
        real function g2in (p, iopt, m, x, f, wk)
          integer    iopt, m
          real       p, x(*), f(*), wk(*)
        end function
      end interface

      interface
        subroutine shouap (n, sy, k, j, sa, nrh, nch, sh,ldh)
          integer    k, j, nch, ldh
          real       sa, sy(*), sh(ldh,*)
        end subroutine
      end interface

      interface
        subroutine shoutr (n, x, k, j, v, beta)
          integer    n, k, j
          real       beta, x(*), v(*)
        end subroutine
      end interface

end module msimslcs

module msimslcd		! MS IMSL common double-precision routines

!
!     Chapter 1:  Linear Systems
!

      interface
        subroutine dlslrg (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lrg (n, a, lda, b, ipath, x, fac, ipvt,wk)
          integer    n, lda, ipath, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlftrg (n, a, lda, fac, ldfac, ipvt)
          integer    n, lda, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2trg (n, a, lda, fac, ldfac, ipvt,scale)
          integer    n, lda, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine dlfsrg (n, fac, ldfac, ipvt, b, ipath, x)
          integer    n, ldfac, ipath, ipvt(*)
          double precision fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfdrg (n, fac, ldfac, ipvt, det1, det2)
          integer    n, ldfac, ipvt(*)
          double precision det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlinrg (n, a, lda, ainv, ldainv)
          integer    n, lda, ldainv
          double precision a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine dl2nrg (n, a, lda, ainv, ldainv, wk, iwk)
          integer    n, lda, ldainv, iwk(*)
          double precision a(lda,*), ainv(ldainv,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlsqrr (nra, nca, a, lda, b, tol, x, res,kbasis)
          integer    nra, nca, lda, kbasis
          double precision tol, a(*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dl2qrr (nra, nca, a, lda, b, tol, x, res, &
			   kbasis, qr, qraux, ipvt, work)
          integer    nra, nca, lda, kbasis, ipvt(*)
          double precision tol, a(lda,*), b(*), x(*), res(*), qr(*),   &
                 qraux(*), work(*)
        end subroutine
      end interface

      interface
        subroutine dlsbrr (nra, nca, a, lda, b, tol, x, res,kbasis)
          integer    nra, nca, lda, kbasis
          double precision tol, a(lda,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dl2brr (nra, nca, a, lda, b, tol, x, res, &
			  kbasis, qr, brrux, ipvt, work)
          integer    nra, nca, lda, kbasis, ipvt(*)
          double precision tol, a(lda,*), b(*), x(*), res(*), qr(*),  &
               brrux(*), work(*)
        end subroutine
      end interface

      interface
        subroutine dlsvrr (nra, nca, a, lda, ipath, tol, &
			   irank, s, u, ldu, v, ldv)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          double precision tol, a(lda,*), s(*), u(ldu,*), v(ldv,*)
        end subroutine
      end interface

      interface
        subroutine dl2vrr (nra, nca, a, lda, ipath, tol, &
			   irank, s, u, ldu, v, ldv, wka, wk)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          double precision tol, a(lda,*), s(*), u(*), v(*), wka(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlsgrr (nra, nca, a, lda, tol, irank,ginva, ldginv)
          integer    nra, nca, lda, irank, ldginv
          double precision tol, a(lda,*), ginva(ldginv,*)
        end subroutine
      end interface

      interface
        subroutine dl2grr (nra, nca, a, lda, tol, irank,ginva, ldginv, wka, wk)
          integer    nra, nca, lda, irank, ldginv
          double precision tol, a(lda,*), ginva(ldginv,*), wka(*), wk(*)
        end subroutine
      end interface

!
!     Chapter 2:  Eigensystem Analysis
!

      interface
        subroutine devcsf (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de5csf (n, a, lda, eval, evec, ldevec, wk,iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision a(lda,*), eval(*), evec(ldevec,*), wk(n,3)
        end subroutine
      end interface

      interface
        double precision function depisf (n, neval, a, lda,eval, evec, ldevec)
          integer    n, neval, lda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function de2isf (n, neval, a, lda, &
			 		  eval, evec, ldevec, work)
          integer    n, neval, lda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*), work(*)
        end function
      end interface

      interface
        subroutine dgvlsp (n, a, lda, b, ldb, eval)
          integer    n, lda, ldb
          double precision a(lda,*), b(ldb,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine dg3lsp (n, a, lda, b, ldb, eval, iwk, wk,s)
          integer    n, lda, ldb, iwk(*)
          double precision a(lda,*), b(ldb,*), eval(*), wk(n,*), s(n+1,*)
        end subroutine
      end interface

!
!     Chapter 3:  Interpolation and Approximation
!

      interface
        subroutine drline (nobs, xdata, ydata, b0, b1, stat)
          integer    nobs
          double precision b0, b1, xdata(*), ydata(*), stat(12)
        end subroutine
      end interface

      interface
        subroutine drcurv (nobs, xdata, ydata, ndeg, b,sspoly, stat)
          integer    nobs, ndeg
          double precision xdata(*), ydata(*), b(*), sspoly(*), stat(10)
        end subroutine
      end interface

      interface
        subroutine dr2urv (nobs, x, y, ndeg, b, sspoly, stat,wk, iwk)
          integer    nobs, ndeg, iwk(*)
          double precision x(*), y(*), b(*), sspoly(*), stat(10), wk(*)
        end subroutine
      end interface

!
!     Chapter 4:  Integration and Differentiation
!

      interface
        subroutine dqdags (f, a, b, errabs, errrel, result,errest)
          double precision f, a, b, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2ags (f, a, b, errabs, errrel, result, errest, maxsub, &
			   neval, nsubin, alist, blist, rlist,elist, iord)
          integer    maxsub, neval, nsubin, iord(*)
          double precision f, a, b, errabs, errrel, result, errest,     &
                alist(*), blist(*), rlist(*), elist(*)
          external   f
        end subroutine
      end interface

!
!     Chapter 6:  Transforms
!

      interface
        subroutine dfftrf (n, seq, coef)
          integer    n
          double precision seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine df2trf (n, seq, coef, wfftr)
          integer    n
          double precision seq(*), coef(*), wfftr(*)
        end subroutine
      end interface

      interface
        subroutine dfftrb (n, coef, seq)
          integer    n
          double precision coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine df2trb (n, coef, seq, wfftr)
          integer    n
          double precision coef(*), seq(*), wfftr(*)
        end subroutine
      end interface

      interface
        subroutine dfftri (n, wfftr)
          integer    n
          double precision wfftr(*)
        end subroutine
      end interface

!
!     Chapter 8:  Optimization
!

      interface
        subroutine duminf (fcn, n, xguess, xscale, fscale,iparam, &
			   rparam, x, fvalue)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),   &
                 rparam(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du2inf (fcn, n, xguess, xscale, fscale,iparam, &
			   rparam, x, fvalue, wk)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),  &
                 rparam(*), x(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du4inf (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

!
!     Chapter 9:  Basic Matrix/Vector Operations
!

      interface
        subroutine dtrnrr (nra, nca, a, lda, nrb, ncb, b,ldb)
          integer    nra, nca, lda, nrb, ncb, ldb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dmxtxf (nra, nca, a, lda, nb, b, ldb)
          integer    nra, nca, lda, nb, ldb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dmxtyf (nra, nca, a, lda, nrb, ncb, b,ldb, nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          double precision a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dmrrrr (nra, nca, a, lda, nrb, ncb, b,ldb, nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          double precision a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        double precision function dblinf (nra, nca, a, lda,x, y)
          integer    nra, nca, lda
          double precision a(*), x(*), y(*)
        end function
      end interface

!
!     Level 1 BLAS
!

      interface
        subroutine dset (n, da, dx, incx)
          integer    n, incx
          double precision da, dx(*)
        end subroutine
      end interface

      interface
        subroutine zset (n, za, zx, incx)
          integer    n, incx
          complex    *16 za, zx(*)
        end subroutine
      end interface

      interface
        subroutine dcopy (n, dx, incx, dy, incy)
          integer    n, incx, incy
          double precision dx(*), dy(*)
        end subroutine
      end interface

      interface
        subroutine zcopy (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end subroutine
      end interface

      interface
        subroutine dscal (n, da, dx, incx)
          integer    n, incx
          double precision da, dx(*)
        end subroutine
      end interface

      interface
        subroutine zscal (n, za, zx, incx)
          integer    n, incx
          complex    *16 za, zx(*)
        end subroutine
      end interface

      interface
        subroutine zdscal (n, da, zx, incx)
          integer    n, incx
          double precision da
          complex    *16 zx(*)
        end subroutine
      end interface

      interface
        subroutine dvcal (n, da, dx, incx, dy, incy)
          integer    n, incx, incy
          double precision da, dx(*), dy(*)
        end subroutine
      end interface

      interface
        subroutine zvcal (n, za, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 za, zx(*), zy(*)
        end subroutine
      end interface

      interface
        subroutine zdvcal (n, da, zx, incx, zy, incy)
          integer    n, incx, incy
          double precision da
          complex    *16 zx(*), zy(*)
        end subroutine
      end interface

      interface
        subroutine dadd (n, da, dx, incx)
          integer    n, incx
          double precision da, dx(*)
        end subroutine
      end interface

      interface
        subroutine zadd (n, za, zx, incx)
          integer    n, incx
          complex    *16 za, zx(*)
        end subroutine
      end interface

      interface
        subroutine dsub (n, da, dx, incx)
          integer    n, incx
          double precision da, dx(*)
        end subroutine
      end interface

      interface
        subroutine zsub (n, za, zx, incx)
          integer    n, incx
          complex    *16 za, zx(*)
        end subroutine
      end interface

      interface
        subroutine daxpy (n, da, dx, incx, dy, incy)
          integer    n, incx, incy
          double precision da, dx(*), dy(*)
        end subroutine
      end interface

      interface
        subroutine zaxpy (n, za, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 za, zx(*), zy(*)
        end subroutine
      end interface

      interface
        subroutine dswap (n, dx, incx, dy, incy)
          integer    n, incx, incy
          double precision dx(*), dy(*)
        end subroutine
      end interface

      interface
        subroutine zswap (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end subroutine
      end interface

      interface
        double precision function ddot (n, dx, incx, dy,incy)
          integer    n, incx, incy
          double precision dx(*), dy(*)
        end function
      end interface

      interface
        complex *16function zdotu (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end function
      end interface

      interface
        complex *16function zdotc (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end function
      end interface

      interface
        complex *16function zqdotc (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end function
      end interface

      interface
        complex *16function zqdotu (n, zx, incx, zy, incy)
          integer    n, incx, incy
          complex    *16 zx(*), zy(*)
        end function
      end interface

      interface
        double precision function dqddot (n, db, dx, incx,dy, incy)
          integer    n, incx, incy
          double precision db, dx(*), dy(*)
        end function
      end interface

      interface
        complex *16function zqcdot (n, za, zx, incx, zy,incy)
          integer    n, incx, incy
          complex    *16 za, zx(*), zy(*)
        end function
      end interface

      interface
        complex *16function zqudot (n, za, zx, incx, zy,incy)
          integer    n, incx, incy
          complex    *16 za, zx(*), zy(*)
        end function
      end interface

      interface
        double precision function dqdoti (n, da, acc, dx,incx, dy, incy)
          integer    n, incx, incy
          double precision da, acc(2), dx(*), dy(*)
        end function
      end interface

      interface
        double precision function dqdota (n, da, acc, dx,incx, dy, incy)
          integer    n, incx, incy
          double precision da, acc(2), dx(*), dy(*)
        end function
      end interface

      interface
        complex *16function zqdoti (n, za, acc, zx, incx, zy,incy)
          integer    n, incx, incy
          double precision acc(4)
          complex    *16 za, zx(*), zy(*)
        end function
      end interface

      interface
        complex *16function zqdota (n, za, acc, zx, incx, zy,incy)
          integer    n, incx, incy
          double precision acc(4)
          complex    *16 za, zx(*), zy(*)
        end function
      end interface

      interface
        subroutine dhprod (n, dx, incx, dy, incy, dz, incz)
          integer    n, incx, incy, incz
          double precision dx(*), dy(*), dz(*)
        end subroutine
      end interface

      interface
        double precision function dxyz (n, dx, incx, dy,incy, dz, incz)
          integer    n, incx, incy, incz
          double precision dx(*), dy(*), dz(*)
        end function
      end interface

      interface
        double precision function dsum (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        double precision function dasum (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        double precision function dzasum (n, zx, incx)
          integer    n, incx
          complex    *16 zx(*)
        end function
      end interface

      interface
        double precision function dnrm2 (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        double precision function dznrm2 (n, zx, incx)
          integer    n, incx
          complex    *16 zx(*)
        end function
      end interface

      interface
        double precision function dprdct (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        integer function idmin (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        integer function idmax (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        integer function idamin (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        integer function izamin (n, zx, incx)
          integer    n, incx
          complex    *16 zx(*)
        end function
      end interface

      interface
        integer function idamax (n, dx, incx)
          integer    n, incx
          double precision dx(*)
        end function
      end interface

      interface
        integer function izamax (n, zx, incx)
          integer    n, incx
          complex    *16 zx(*)
        end function
      end interface

      interface
        subroutine drotg (da, db, dc, ds)
          double precision da, db, dc, ds
        end subroutine
      end interface

      interface
        subroutine drot (n, dx, incx, dy, incy, c, s)
          integer    n, incx, incy
          double precision c, s, dx(*), dy(*)
        end subroutine
      end interface

      interface
        subroutine zdrot (n, zx, incx, zy, incy, c, s)
          integer    n, incx, incy
          double precision c, s
          complex    *16 zx(*), zy(*)
        end subroutine
      end interface

      interface
        subroutine drotmg (dd1, dd2, dx1, dy1, dparam)
          double precision dd1, dd2, dx1, dy1, dparam(5)
        end subroutine
      end interface

      interface
        subroutine drotm (n, dx, incx, dy, incy, dparam)
          integer    n, incx, incy
          double precision dx(*), dy(*), dparam(5)
        end subroutine
      end interface

      interface
        subroutine zdrotm (n, zx, incx, zy, incy, dparam)
          integer    n, incx, incy
          double precision dparam(5)
          complex    *16 zx(*), zy(*)
        end subroutine
      end interface

!
!     Level 2 BLAS
!

      interface
        subroutine dgemv (trans, m, n, alpha, a, lda, x,incx, beta, y, incy)
          integer    m, n, lda, incx, incy
          double precision alpha, beta, a(*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine zgemv (trans, m, n, alpha, a, lda, x,incx, beta, y, incy)
          integer    m, n, lda, incx, incy
          complex    *16 alpha, beta, a(*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine dgbmv (trans, m, n, nlca, nuca, alpha, a, &
			  lda, x, incx, beta, y, incy)
          integer    m, n, nlca, nuca, lda, incx, incy
          double precision alpha, beta, a(lda,*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine zgbmv (trans, m, n, nlca, nuca, alpha, a, &
			  lda, x, incx, beta, y, incy)
          integer    m, n, nlca, nuca, lda, incx, incy
          complex    *16 alpha, beta, a(lda,*), x(*), y(*)
          character  trans*1
        end subroutine
      end interface

      interface
        subroutine zhemv (uplo, n, alpha, a, lda, x, incx,beta, y, incy)
          integer    n, lda, incx, incy
          complex    *16 alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine zhbmv (uplo, n, ncoda, alpha, a, lda, x,incx, beta, y, incy)
          integer    n, ncoda, lda, incx, incy
          complex    *16 alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine dsymv (uplo, n, alpha, a, lda, x, incx,beta, y, incy)
          integer    n, lda, incx, incy
          double precision alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine dsbmv (uplo, n, ncoda, alpha, a, lda, x,incx, beta, y, incy)
          integer    n, ncoda, lda, incx, incy
          double precision alpha, beta, a(lda,*), x(*), y(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine dtrmv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          double precision a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztrmv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          complex    *16 a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine dtbmv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          double precision a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztbmv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          complex    *16 a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine dtrsv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          double precision a(*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztrsv (uplo, trans, diag, n, a, lda, x,incx)
          integer    n, lda, incx
          complex    *16 a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine dtbsv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          double precision a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztbsv (uplo, trans, diag, n, ncoda, a,lda, x, incx)
          integer    n, ncoda, lda, incx
          complex    *16 a(lda,*), x(*)
          character  uplo*1, trans*1, diag*1
        end subroutine
      end interface

      interface
        subroutine dger (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          double precision alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine zgeru (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          complex    *16 alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine zgerc (m, n, alpha, x, incx, y, incy, a,lda)
          integer    m, n, incx, incy, lda
          complex    *16 alpha, x(*), y(*), a(*)
        end subroutine
      end interface

      interface
        subroutine zher (uplo, n, alpha, x, incx, a, lda)
          integer    n, incx, lda
          double precision alpha
          complex    *16 x(*), a(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine zher2 (uplo, n, alpha, x, incx, y, incy,a, lda)
          integer    n, incx, incy, lda
          complex    *16 alpha, x(*), y(*), a(lda,*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine dsyr (uplo, n, alpha, x, incx, a, lda)
          integer    n, incx, lda
          double precision alpha, x(*), a(*)
          character  uplo*1
        end subroutine
      end interface

      interface
        subroutine dsyr2 (uplo, n, alpha, x, incx, y, incy,a, lda)
          integer    n, incx, incy, lda
          double precision alpha, x(*), y(*), a(lda,*)
          character  uplo*1
        end subroutine
      end interface

!
!     Level 3 BLAS
!

      interface
        subroutine dgemm (transa, transb, m, n, k, alpha, a, &
 			  lda, b, ldb, beta, c, ldc)
          integer    m, n, k, lda, ldb, ldc
          double precision alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  transa*1, transb*1
        end subroutine
      end interface

      interface
        subroutine zgemm (transa, transb, m, n, k, alpha, a, &
			  lda, b, ldb, beta, c, ldc)
          integer    m, n, k, lda, ldb, ldc
          complex    *16 alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  transa*1, transb*1
        end subroutine
      end interface

      interface
        subroutine dsymm (side, uplo, m, n, alpha, a, lda, b,ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          double precision alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine zsymm (side, uplo, m, n, alpha, a, lda, b,ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          complex    *16 alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine zhemm (side, uplo, m, n, alpha, a, lda, b,ldb, beta, c, ldc)
          integer    m, n, lda, ldb, ldc
          complex    *16 alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  side*1, uplo*1
        end subroutine
      end interface

      interface
        subroutine dsyrk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          double precision alpha, beta, a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine zsyrk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          complex    *16 alpha, beta, a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine zherk (uplo, trans, n, k, alpha, a, lda,beta, c, ldc)
          integer    n, k, lda, ldc
          double precision alpha, beta
          complex    *16 a(lda,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine dsyr2k (uplo, trans, n, k, alpha, a, lda, &
			   b, ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          double precision alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine zsyr2k (uplo, trans, n, k, alpha, a, lda, &
			   b, ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          complex    *16 alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine zher2k (uplo, trans, n, k, alpha, a, lda, &
			   b, ldb, beta, c, ldc)
          integer    n, k, lda, ldb, ldc
          double precision beta
          complex    *16 alpha, a(lda,*), b(ldb,*), c(ldc,*)
          character  uplo*1, trans*1
        end subroutine
      end interface

      interface
        subroutine dtrmm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          double precision alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztrmm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          complex    *16 alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine dtrsm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          double precision alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

      interface
        subroutine ztrsm (side, uplo, transa, diag, m, n,alpha, a, lda, b, ldb)
          integer    m, n, lda, ldb
          complex    *16 alpha, a(lda,*), b(ldb,*)
          character  side*1, uplo*1, transa*1, diag*1
        end subroutine
      end interface

!
!     Chapter 10:  Utilities
!

      interface
        subroutine dwrrrn (title, nra, nca, a, lda, itring)
          integer    nra, nca, lda, itring
          double precision a(*)
          character  title*(*)
        end subroutine
      end interface

      interface
        subroutine dwrrrl (title, nra, nca, a, lda, itring,fmt, rlabel, clabel)
          integer    nra, nca, lda, itring
          double precision a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*)
        end subroutine
      end interface

      interface
        subroutine dw2rrl (title, nra, nca, a, lda, itring, &
			   fmt, rlabel, clabel, chwk)
          integer    nra, nca, lda, itring
          double precision a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*),  &
                     chwk(*)*10
        end subroutine
      end interface

      interface
        subroutine dpermu (n, x, ipermu, ipath, xpermu)
          integer    n, ipath, ipermu(*)
          double precision x(*), xpermu(*)
        end subroutine
      end interface

      interface
        subroutine dperma (nra, nca, a, lda, ipermu, ipath,aper, ldaper)
          integer    nra, nca, lda, ipath, ldaper, ipermu(*)
          double precision a(lda,*), aper(ldaper,*)
        end subroutine
      end interface

      interface
        subroutine dp2rma (nra, nca, a, lda, ipermu, ipath,aper, ldaper, work)
          integer    nra, nca, lda, ipath, ldaper, ipermu(*)
          double precision a(lda,*), aper(ldaper,*), work(*)
        end subroutine
      end interface

      interface
        subroutine dsvrgn (n, ra, rb)
          integer    n
          double precision ra(*), rb(*)
        end subroutine
      end interface

      interface
        subroutine dsvrgp (n, ra, rb, iperm)
          integer    n, iperm(*)
          double precision ra(*), rb(*)
        end subroutine
      end interface

      interface
        subroutine dsrch (n, value, x, incx, index)
          integer    n, incx, index
          double precision value, x(*)
        end subroutine
      end interface

      interface
        double precision function drnunf ()
        end function
      end interface

      interface
        subroutine drnun (nr, r)
          integer    nr
          double precision r(*)
        end subroutine
      end interface

      interface
        subroutine dplotp (ndata, nfun, x, a, lda, inc, &
			   range, symbol, xtitle, ytitle, title)
          integer    ndata, nfun, lda, inc
          double precision x(*), a(lda,*), range(4)
          character  symbol*(*), xtitle*(*), ytitle*(*), title*(*)
        end subroutine
      end interface

!
!     Special Functions Chapter 4:  Gamma Function and Related Functions
!

      interface
        double precision function dbinom (n, m)
          integer    n, m
        end function
      end interface

      interface
        double precision function dgamma (x)
          double precision x
        end function
      end interface

      interface
        subroutine d9gaml (xmin, xmax)
          double precision xmin, xmax
        end subroutine
      end interface

      interface
        double precision function dbeta (a, b)
          double precision a, b
        end function
      end interface

!
!     Special Functions Chapter 11:  Probability Distribution Functions
!                                    and Inverses
!

      interface
        double precision function dbindf (k, n, p)
          integer    k, n
          double precision p
        end function
      end interface

      interface
        double precision function dbinpr (k, n, p)
          integer    k, n
          double precision p
        end function
      end interface

      interface
        double precision function dhypdf (k, n, m, l)
          integer    k, n, m, l
        end function
      end interface

      interface
        double precision function dhyppr (k, n, m, l)
          integer    k, n, m, l
        end function
      end interface

      interface
        double precision function dpoidf (k, theta)
          integer    k
          double precision theta
        end function
      end interface

      interface
        double precision function dpoipr (k, theta)
          integer    k
          double precision theta
        end function
      end interface

      interface
        double precision function dks1df (nobs, d)
          integer    nobs
          double precision d
        end function
      end interface

      interface
        double precision function dk21df (nobs, d, wk)
          integer    nobs
          double precision d, wk(*)
        end function
      end interface

      interface
        double precision function dks2df (nobsx, nobsy, d)
          integer    nobsx, nobsy
          double precision d
        end function
      end interface

      interface
        double precision function dk22df (nobsx, nobsy, d,wk)
          integer    nobsx, nobsy
          double precision d, wk(*)
        end function
      end interface

      interface
        double precision function dnordf (x)
          double precision x
        end function
      end interface

      interface
        double precision function dnorin (p)
          double precision p
        end function
      end interface

      interface
        double precision function dbetdf (x, pin, qin)
          double precision x, pin, qin
        end function
      end interface

      interface
        double precision function dbetin (p, pin, qin)
          double precision p, pin, qin
        end function
      end interface

      interface
        double precision function dbnrdf (x, y, rho)
          double precision x, y, rho
        end function
      end interface

      interface
        double precision function dchidf (chsq, df)
          double precision chsq, df
        end function
      end interface

      interface
        double precision function dchiin (p, df)
          double precision p, df
        end function
      end interface

      interface
        double precision function dcsndf (chsq, df, alam)
          double precision chsq, df, alam
        end function
      end interface

      interface
        double precision function dfdf (f, dfn, dfd)
          double precision f, dfn, dfd
        end function
      end interface

      interface
        double precision function dfin (p, dfn, dfd)
          double precision p, dfn, dfd
        end function
      end interface

      interface
        double precision function dgamdf (x, a)
          double precision x, a
        end function
      end interface

      interface
        double precision function dtdf (t, df)
          double precision t, df
        end function
      end interface

      interface
        double precision function dtin (p, df)
          double precision p, df
        end function
      end interface

      interface
        double precision function dtndf (t, idf, delta)
          integer    idf
          double precision t, delta
        end function
      end interface

      interface
        double precision function dgcdf (x0, iopt, m, x, f)
          integer    iopt, m
          double precision x0, x(*), f(*)
        end function
      end interface

      interface
        double precision function dg4df (x0, iopt, m, x, f,wk, iwk)
          integer    iopt, m, iwk(*)
          double precision x0, x(*), f(*), wk(*)
        end function
      end interface

      interface
        double precision function dgcin (p, iopt, m, x, f)
          integer    iopt, m
          double precision p, x(*), f(*)
        end function
      end interface

      interface
        double precision function dg3in (p, iopt, m, x, f,wk, iwk)
          integer    iopt, m, iwk(*)
          double precision p, x(*), f(*), wk(*)
        end function
      end interface

!
!     Reference Material
!

      interface
        double precision function dmach (n)
          integer    n
        end function
      end interface

      interface
        logical function difnan(x)
          double precision x
        end function
      end interface

!
!     Deprecated Routines
!

      interface
        double precision function dg2df (x0, iopt, m, x, f,wk)
          integer    iopt, m
          double precision x0, x(*), f(*), wk(*)
        end function
      end interface

      interface
        double precision function dg3df (xx0, iopt, m, xx, f,c)
          integer    iopt, m
          double precision xx0, xx(*), f(*), c(*)
        end function
      end interface

      interface
        double precision function dg2in (p, iopt, m, x, f,wk)
          integer    iopt, m
          double precision p, x(*), f(*), wk(*)
        end function
      end interface

      interface
        subroutine dhouap (n, dy, k, j, da, nrh, nch, dh,ldh)
          integer    k, j, nch, ldh
          double precision da, dy(*), dh(ldh,*)
        end subroutine
      end interface

      interface
        subroutine dhoutr (n, x, k, j, v, beta)
          integer    n, k, j
          double precision beta, x(*), v(*)
        end subroutine
      end interface

end module msimslcd

module msimslmc		! MS IMSL common math routines

!
!     Chapter 9:  Basic Matrix/Vector Operations
!

      interface
        subroutine dqadd (a, acc)
          double precision a, acc(2)
        end subroutine
      end interface

      interface
        subroutine dqini (da, dacc)
          double precision da, dacc(2)
        end subroutine
      end interface

      interface
        subroutine dqmul (da, db, dacc)
          double precision da, db, dacc(2)
        end subroutine
      end interface

      interface
        subroutine dqsto (dacc, da)
          double precision da, dacc(2)
        end subroutine
      end interface

      interface
        subroutine zqadd (za, zacc)
          double precision zacc(4)
          complex    *16 za
        end subroutine
      end interface

      interface
        subroutine zqini (za, zacc)
          double precision zacc(4)
          complex    *16 za
        end subroutine
      end interface

      interface
        subroutine zqmul (za, zb, zacc)
          double precision zacc(4)
          complex    *16 za, zb
        end subroutine
      end interface

      interface
        subroutine zqsto (zacc, za)
          double precision zacc(4)
          complex    *16 za
        end subroutine
      end interface

!
!     Chapter 10:  Utilities
!

      interface
        subroutine svibn (n, ia, ib)
          integer    n, ia(*), ib(*)
        end subroutine
      end interface

      interface
        subroutine svibp (n, ia, ib, iperm)
          integer    n, ia(*), ib(*), iperm(*)
        end subroutine
      end interface

      interface
        character *(*)function verml (iselct)
          integer    iselct
        end function
      end interface

      interface
        subroutine iumag (prodnm, ichp, iact, numopt, iopts,ivals)
          integer    ichp, iact, numopt, iopts(*), ivals(*)
          character  prodnm*4
        end subroutine
      end interface

      interface
        subroutine sumag (prodnm, ichp, iact, numopt, iopts,svals)
          integer    ichp, iact, numopt, iopts(*)
          real       svals(*)
          character  prodnm*4
        end subroutine
      end interface

      interface
        subroutine dumag (prodnm, ichp, iact, numopt, iopts,dvals)
          integer    ichp, iact, numopt, iopts(*)
          double precision dvals(*)
          character  prodnm*4
        end subroutine
      end interface

      interface
        subroutine prime (n, npf, ipf, iexp, ipw)
          integer    n, npf, ipf(13), iexp(13), ipw(13)
        end subroutine
      end interface

!
!     Chapter 10:  Utilities
!

      interface
        complex function cexprl (z)
          complex    z
        end function
      end interface

      interface
        complex function clnrel (z)
          complex    z
        end function
      end interface

!
!     Special Functions Chapter 4:  Gamma Function and Related Functions
!

      interface
        complex function cgamma (z)
          complex    z
        end function
      end interface

      interface
        complex function cgamr (z)
          complex    z
        end function
      end interface

      interface
        complex function clngam (zin)
          complex    zin
        end function
      end interface

      interface
        complex function cpsi (zin)
          complex    zin
        end function
      end interface

      interface
        complex function cbeta (a, b)
          complex    a, b
        end function
      end interface

      interface
        complex function clbeta (a, b)
          complex    a, b
        end function
      end interface

!
!     Deprecated Routines
!

      interface
        subroutine czadd (a, acc)
          double precision acc(4)
          complex    a
        end subroutine
      end interface

      interface
        subroutine czini (z, acc)
          double precision acc(4)
          complex    z
        end subroutine
      end interface

      interface
        subroutine czmul (a, b, acc)
          double precision acc(4)
          complex    a, b
        end subroutine
      end interface

      interface
        subroutine czsto (acc, z)
          double precision acc(4)
          complex    z
        end subroutine
      end interface

      interface
        subroutine sdadd (a, acc)
          real       a
          double precision acc(2)
        end subroutine
      end interface

      interface
        subroutine sdini (s, acc)
          real       s
          double precision acc(2)
        end subroutine
      end interface

      interface
        subroutine sdmul (a, b, acc)
          real       a, b
          double precision acc(2)
        end subroutine
      end interface

      interface
        subroutine sdsto (acc, s)
          real       s
          double precision acc(2)
        end subroutine
      end interface

end module msimslmc

module msimslsc		! MS IMSL common stat routines

!
!     Chapter 11:  Cluster Analysis
!

      interface
        subroutine cnumb (node, iclson, icrson, k, iclus,nclus)
          integer    node, k, iclson(*), icrson(*), iclus(*), nclus(*)
        end subroutine
      end interface

      interface
        subroutine c2umb (node, iclson, icrson, k, iclus,nclus, ipt)
          integer    node, iclson(*), icrson(*), k, iclus(*), nclus(*), ipt(*)
        end subroutine
      end interface

!
!     Chapter 18:  Random Number Generation
!

      interface
        subroutine rnopg (iopt)
          integer    iopt
        end subroutine
      end interface

      interface
        subroutine rnsef (iarray)
          integer    iarray(1565)
        end subroutine
      end interface

      interface
        subroutine rngef (iarray)
          integer    iarray(1565)
        end subroutine
      end interface

      interface
        subroutine rnisd (iseed1, iseed2)
          integer    iseed1, iseed2
        end subroutine
      end interface

      interface
        subroutine rnbin (nr, n, p, ir)
          integer    nr, n, ir(*)
          real       p
        end subroutine
      end interface

      interface
        subroutine rngeo (nr, p, ir)
          integer    nr, ir(*)
          real       p
        end subroutine
      end interface

      interface
        subroutine rnhyp (nr, n, m, l, ir)
          integer    nr, n, m, l, ir(*)
        end subroutine
      end interface

      interface
        subroutine rnlgr (nr, a, ir)
          integer    nr, ir(*)
          real       a
        end subroutine
      end interface

      interface
        subroutine rnnbn (nr, rk, p, ir)
          integer    nr, ir(*)
          real       rk, p
        end subroutine
      end interface

      interface
        subroutine rnpoi (nr, theta, ir)
          integer    nr, ir(*)
          real       theta
        end subroutine
      end interface

      interface
        subroutine rnund (nr, k, ir)
          integer    nr, k, ir(*)
        end subroutine
      end interface

      interface
        subroutine rnmtn (nr, n, k, p, ir, ldir)
          integer    nr, n, k, ldir, ir(ldir,*)
          real       p(*)
        end subroutine
      end interface

      interface
        subroutine rntab (ido, nrow, ncol, nrtot, nctot,itab, lditab)
          integer    ido, nrow, ncol, lditab, nrtot(*), nctot(*), itab(*)
        end subroutine
      end interface

      interface
        subroutine r2tab (ido, nrow, ncol, nrtot, nctot, &
			  itab, lditab, iopt, irsum, iwk, wk)
          integer    ido, nrow, ncol, lditab, iopt, irsum, nrtot(*),  &
                nctot(*), itab(lditab,*), iwk(*)
          real       wk(*)
        end subroutine
      end interface

      interface
        subroutine rnper (k, iper)
          integer    k, iper(*)
        end subroutine
      end interface

      interface
        subroutine rnsri (nsamp, npop, index)
          integer    nsamp, npop, index(*)
        end subroutine
      end interface

!
!     Chapter 19:  Utilities
!

      interface
        character *(*)function versl (iselct)
          integer    iselct
        end function
      end interface

end module msimslsc
