!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KIJOPT(NC)
C
C  PURPOSE:  BP CALCULATIONS TO DETERMINE THE OPTIMAL
C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
C  PROGRAMMED BY:  JRE 2/96
C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
C  REFERENCE:NUMERICAL RECIPES
C
      IMPLICIT DOUBLE PRECISION(A-H,K,O-Z)
      CHARACTER*44 FN
      PARAMETER (NMX=55)
      PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
      
      DIMENSION TC(NMX),PC(NMX),ACEN(NMX),ID(NMX)
      DIMENSION deviate(222),parm(1)
	common/ppData/TC,PC,ACEN,ID
	common/tpxyData/tDat(111),PDAT(111),XDAT(111),npts
      COMMON/eta/etaL,etaV,ZL,ZV
	dimension KIJ(NMX,NMX),KTIJ(NMX,NMX)
	&		,HIJ(NMX,NMX),HTIJ(NMX,NMX)
	&   ,xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)
      COMMON/BIPs/KIJ,KTIJ,HIJ,HTIJ,xsTau,xsTauT,xsAlpha

      WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (PATH\FN.FT)'
      WRITE(*,*)'FIRST LINE =NPTS, 2ND LINE = TK, PMPA, X,Y, ...'
      READ(*,'(A)')FN
      OPEN(51,FILE=FN)
      READ(51,*)NPTS
      DO 5 I=1,NPTS
        READ(51,*)tDat(I),PDAT(I),XDAT(I)
        WRITE(*,*)tDat(I),PDAT(I),XDAT(I)
5	CONTINUE
	MAXIT=1111
      INIT=1
      WRITE(*,*)'THIS ROUTINE ASSUMES A BINARY SYSTEM FROM MAIN PROGRAM'
      WRITE(6,*)'ENTER RANGE FOR KIJ TO OPTIMIZE, (KIJLO KIJHI)  '
      READ(5,*)KIJ0,KIJ3                                          
      KIJ1=KIJ0+CGOLD*(KIJ3-KIJ0)
      KIJ2=KIJ0+RGOLD*(KIJ3-KIJ0)
	nParms=1
      L=0
	parm(1)=KIJ0
      call KIJVAL(nc,nParms,parm,deviate,PAAD0,L)
	parm(1)=KIJ1
      call KIJVAL(nc,nParms,parm,deviate,PAAD1,L)
	parm(1)=KIJ2
      call KIJVAL(nc,nParms,parm,deviate,PAAD2,L)
	parm(1)=KIJ3
      call KIJVAL(nc,nParms,parm,deviate,PAAD3,L)
      TOL=0.002      
1000  IF(ABS(KIJ3-KIJ0).GT.TOL*(ABS(KIJ1)+ABS(KIJ2)))THEN
        IF(PAAD2.LT.PAAD1)THEN
          KIJ0=KIJ1
          KIJ1=KIJ2
          KIJ2=RGOLD*KIJ1+CGOLD*KIJ3
          PAAD0=PAAD1
          PAAD1=PAAD2
			parm(1)=KIJ2
			call KIJVAL(nc,nParms,parm,deviate,PAAD2,L)
        ELSE
          KIJ3=KIJ2
          KIJ2=KIJ1
          KIJ1=RGOLD*KIJ2+CGOLD*KIJ0
          PAAD3=PAAD2
          PAAD2=PAAD1
			parm(1)=KIJ1
			call KIJVAL(nc,nParms,parm,deviate,PAAD1,L)
        ENDIF
      GOTO 1000
      ENDIF
      IF(PAAD1.LT.PAAD2)THEN
        PAADP=PAAD1
        KIJB=KIJ1
      ELSE
        PAADP=PAAD2
        KIJB=KIJ2
      ENDIF    
	OPEN(67,FILE='KIJOUT.txt')
      WRITE(*,*)'KIJOPT, PAADP'
      WRITE(*,606)KIJB,PAADP
      WRITE(67,*)'KIJOPT, PAADP'
      WRITE(67,606)KIJB,PAADP
	write(*,600)
	write(67,600)
600   format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )

     	LAST=1
	parm(1)=KIJb
      call KIJVAL(nc,nParms,parm,deviate,PAAD1,Last)

606   FORMAT(1X,F8.4,1X,F10.2)      
	IF(LAST.EQ.1)CLOSE(67)
      close(51)
	pause 'These results are tabulated in KIJOut.txt.'
	RETURN
	END
