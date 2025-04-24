      SUBROUTINE MINSSQ(FUNCS,X,Y,NPTS,NPARMS,A,stdErr,CHISQ)
C      PROGRAM minssq
C     driver for routine mrqmin
C  THE USER MUST SUPPLY FUNCS(NDATA,NPARMS,X,Y,DEVIAT)
C    WHERE ALL DEVIATIONS ARE EVALUATED AT ONCE.  
C  PROGRAMS REQUIRED:  FUNCS,DMRQMIN,FDMRQCOF,FDJAC,DGAUSSJ,DCOVSRT
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(MAXA=20,maxdata=1111)
C  NOTE:  MAXA SET HERE MUST BE IDENTICAL TO THE ORIGINAL DIMENSION OF 
C    COVAR IN THE CALLING ROUTINE.      
      DIMENSION x(*),y(*),ia(maxa),ipvt2(maxa),wa(maxdata)
     *,a(NPARMS),covar(MAXA,MAXA),alpha(MAXA,MAXA)
	1,stdErr(nParms),sig(maxdata)
C  ALPHA HERE IS JUST A WORKSPACE     
      EXTERNAL FUNCS
      mfit=NPARMS
      do 10 i=1,mfit
        ia(i)=1
10    continue

	SIG(2)=1.D-11
      IF(ABS( SIG(2) ).LT.1.D-5)THEN
        DO 15 I=1,NPTS
          SIG(I)=1
          !SIG(I)=1.D-5
15      CONTINUE
      END IF
                
      alamda=-1
      call mrqmin(x,y,sig,NPTS,a,ia,NPARMS,covar,alpha,
     *       MAXA,chisq,FUNCS,alamda)
      k=1
      itst=0
      ALAMDA=.001
20    CONTINUE
        k=k+1
        if(k.gt.111)write(*,*)'too many minmrqs'
        if(k.gt.111)return
        ochisq=chisq
        call mrqmin(x,y,sig,NPTS,a,ia,NPARMS,covar,alpha,
     *       MAXA,chisq,FUNCS,alamda)
      if ( (abs(ochisq-chisq)/chisq).gt.1.D-4 ) go to 20

      alamda=0.0d0
      call mrqmin(x,y,sig,NPTS,a,ia,NPARMS,covar,alpha,
     *       MAXA,chisq,FUNCS,alamda)
C        write(*,*) 'Uncertainties:'
C        write(*,'(1x,6f8.4/)') (sqrt(covar(i,i)),i=1,NPARMS)
      call MATINV(covar,MAXA,NPARMS,IPVT2,WA,INFO)

	do iParm=1,nParms
		do jParm=iParm,nParms
			covar(iParm,jParm)=sqrt( chisq*ABS(covar(iParm,jParm)) )
			covar(jParm,iParm)=covar(iParm,jParm)
		enddo
		stdErr(iParm)=sqrt( covar(ipvt2(iParm),ipvt2(iParm)) )
	enddo
      RETURN
      END
      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      DOUBLE PRECISION alamda,chisq,a(ma),alpha(nca,nca),
     * covar(nca,nca),sig(ndata),x(*),y(*)
      PARAMETER (MMAX=20)
	EXTERNAL funcs
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,m,mfit
      DOUBLE PRECISION ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.d0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.01d0
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.d0+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.d0)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1d0*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.d0*alamda
        chisq=ochisq
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j[)1y.
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs)
C  THIS ROUTINE USES FINITE DIFFERENCES TO ESTIMATE ALPHA AND BETA. 
C    IT IS VERY GENERAL IN THAT ANY SET OF DEVIATION FUNCTIONS CAN 
C    BE MINIMIZED BY SUM OF SQUARES, INCLUDING THOSE WHERE THE VALUE OF
C    THE DEVIATIONS MUST BE EVALUATED FOR ALL "DATA" AT THE SAME TIME.
C    AN EXAMPLE WHERE THIS IS USEFUL IS WHEN THE FUNCTION VALUES ARE
C    FOURIER TRANSFORMS OF THE ORIGINAL DATA.  
C  SINCE THE ROUTINE USES FINITE DIFFERENCES, THE USER MUST ONLY SUPPLY
C    A ROUTINE FUNCS(NDATA,NPARMS,X,Y,DEVIAT).THE COMPUTER DOES THE REST.
C                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MAXDAT=200,MAXA=20)
      EXTERNAL funcs
C  MAXDAT IS THE MAXIMUM NUMBER OF DATA POINTS ALLOWABLE.  IT IS LOCAL
C    TO THIS SUBROUTINE.
C  MAXA=NALP IS THE MAXIMUM NUMBER OF ADJUSTABLE PARAMATERS
C    NALP IS USED FOR DIMENSIONING PURPOSES AND MUST BE EXACTLY
C    AS SPECIFIED IN THE ORIGINAL DIMENSION STATEMENT FOR ALPHA
C
      DIMENSION IA(*),A(*),X(*),Y(*),SIG(*)
     & ,FJAC(MAXDAT,MAXA),ALPHA(NALP,NALP),BETA(NALP),DEVIAT(MAXDAT)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      NPARMS=MFIT
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.d0
12      continue
        beta(j)=0.d0
13    continue
      CALL FUNCS(NDATA,NPARMS,a,X,Y,DEVIAT) 
C      CALL FUNCS(NDATA,NPARMS,a,DEVIAT) 
      CALL FDJAC(FUNCS,NDATA,NPARMS,A,X,Y,DEVIAT,FJAC,MAXDAT)
      chisq=0.d0
      do 16 i=1,ndata
        sig2i=1.d0/(sig(i)*sig(i))
        CHISQ=CHISQ+DEVIAT(I)*DEVIAT(I)*SIG2i
        DO 15 J=1,NPARMS
          BETA(J)=BETA(J)-DEVIAT(I)*FJAC(I,J)*SIG2I
          DO 14 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+FJAC(I,J)*FJAC(I,K)*SIG2I
14        CONTINUE
15      CONTINUE
16    CONTINUE
C      WRITE(*,*)'IN FDMRQCOF,CHISQ,C',CHISQ,(A(I),I=1,4)

C  FILL IN THE SYMMETRIC SIDE.

      do 18 j=2,MFIT
        DO 17 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
17      CONTINUE
18    CONTINUE          
      RETURN
      END            
      SUBROUTINE FDJAC(FUNCS,M,N,X,XDATA,YDATA,FVEC,FJAC,LDFJAC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(MAXWA=200)
      DIMENSION X(N),XDATA(M),YDATA(M),FVEC(M),FJAC(LDFJAC,N),WA(MAXWA)
	external funcs
C     **********
C
C     THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
C     TO THE M BY N JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED
C     PROBLEM OF M FUNCTIONS IN N VARIABLES.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE FDJAC(FCN,M,N,X,XDATA,YDATA,FVEC,FJAC,LDFJAC)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES THE FUNCTIONS. FCN MUST BE DECLARED
C         IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(M,N,X,XDATA,YDATA,FVEC)
C         INTEGER M,N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(M),XDATA(M),YDATA(M)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF FUNCTIONS.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES. N MUST NOT EXCEED M.
C
C       X IS AN INPUT ARRAY OF LENGTH N GIVING THE PARAMETER VALUES.
C
C       FVEC IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE
C         FUNCTIONS EVALUATED AT X.
C
C       FJAC IS AN OUTPUT M BY N ARRAY WHICH CONTAINS THE
C         APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       WA IS A WORK ARRAY OF LENGTH M.
C
C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ...... FCN
C
C     **********
      DATA ZERO /0.0D0/
      EPS = 1.D-7
      DO 20 J = 1, N
         TEMP = X(J)
         H = EPS*DABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
c      write(*,*)'n,m,h,xfx',n,m,h,x
         CALL FUNCS(M,N,X,XDATA,YDATA,WA)
C         CALL FUNCS(M,N,X,WA)
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
c            write(*,*)'i,j,fjac(i,j)',i,j,fjac(i,j)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE FDJAC.
C
      END
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) pause 'singular matrix in gaussj'
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j[)1y.
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      DOUBLE PRECISION covar(npc,npc)
      INTEGER i,j,k
      DOUBLE PRECISION swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.d0
          covar(j,i)=0.d0
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j[)1y.
      SUBROUTINE MATINV(A,NROWA,N,IPVT,WORK,INFO)
	!In-place matrix inverse.  A becomes inverse of A on return.
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION A(NROWA,N),IPVT(N),WORK(*),DET(2)
      CALL DGEFA(A,NROWA,N,IPVT,INFO)
      CALL DGEDI(A,NROWA,N,IPVT,DET,WORK,01)
      RETURN
      END
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),DET(2),WORK(1)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
C  10240 BYTES REQUIRED FOR MEMBER DGEDI
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C  3760 BYTES REQUIRED FOR MEMBER DAXPY
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
C  3280 BYTES REQUIRED FOR MEMBER DSCAL
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
C  4400 BYTES REQUIRED FOR MEMBER DSWAP
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
C  8240 BYTES REQUIRED FOR MEMBER DGEFA
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
C  2960 BYTES REQUIRED FOR MEMBER IDAMAX

