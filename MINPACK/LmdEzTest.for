C
C  DRIVER FOR LMDIF EXAMPLE.
C  For FCN, THE ANSWERS ARE: 0.0824,1.1330,2.3437 => FNORM=.0906
C
c  Example. Polynomial Regression
c	xData/0.05,0.08,0.12,0.15,0.17,0.2,0.23,0.25,0.28,0.3,0.33,0.35,0.38,0.4,0.42,0.43,0.45/
c	yData/0.5084,0.8424,1.3272,1.7183,1.9867,2.4122,2.8478,3.1429,3.5972,3.8954,4.3444,4.6342,5.0583,5.3285,5.5840,5.7083,5.9486/
c	LMDifEZ gives with unweighted deviates:
c	9.411	13.657	5.638	-38.289 	0.058	0.065	0.595	1.317
c	Poly0 Coeffs:											StdErr:
c	9.4099	13.6615	5.6290	-38.2874		0.06  	0.07  	0.60  	1.32
c	Excel gives with unweighted deviates:
c	Poly0 Coeffs:											StdErr:
c	9.4099	13.6615	5.6290	-38.2874		0.0545	0.6231	2.2380	2.5313
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(nData=17,nParms=4,maxParms=20,maxData=111)
c      dimension iPvt(nParms),workArray(LWA)
      DOUBLE PRECISION parm(nParms)
	1                 ,sig(maxData),x(nData),y(nData),fvec(maxData)
	doublePrecision stdErr(maxParms)
      EXTERNAL DevalSsq,FCN
	data x/0.05,0.08,0.12,0.15,0.17,0.2,0.23,0.25,0.28,0.3,0.33,0.35
	1     ,0.38,0.4,0.42,0.43,0.45/
	data y/0.5084,0.8424,1.3272,1.7183,1.9867,2.4122,2.8478,3.1429
	1  ,3.5972,3.8954,4.3444,4.6342,5.0583,5.3285,5.5840,5.7083,5.9486/
	common/data/xDat(maxData),yDat(maxData)
	do iDat=1,nData	!this is necessary b/c the data are initialized in a data statement w/o block data.
		xDat(iDat)=x(iDat)
		yDat(iDat)=y(iDat)
	enddo
C
C  THE FOLLOWING STARTING VALUES PROVIDE A ROUGH FIT.
C
      parm(1)=1.0
      parm(2)=1.0
      parm(3)=1.0
      parm(4)=1.0
	nDataFcn=15	!For FCN example
      nParmsFcn=3
C  
C  THE FOLLOWING STARTING VALUES PROVIDE A ROUGH FIT.
C

C
C  SET tol TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C  UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED, 
C  THIS IS THE RECOMMENDED SETTING.
C
      tol=1e-4
	factor=0.001D0
	do idat=1,nData
		sig(idat)=1
	enddo
C
c      call minssq(Deval,nData,nParms,parm,factor,deviate,tol,
c	1								iErrCode,stdErr)
c      call LmDifEz(FCN,nDataFcn,nParmsFcn,parm,factor,FVEC,TOL,
c	1								iErrCode,stdErr)
      call LmDifEz(DevalSsq,nData,nParms,parm,factor,FVEC,TOL,
	1								iErrCode,stdErr)
c      call MINSSQ(DevalSsq,X,Y,nData,nParms,parm
c	1            ,stdErr,CHISQ)
      !fNorm=ENORM(nData,deviate)
      !WRITE(*,1000) fNorm,iErrCode
	!WRITE(*,'(8F10.3)')(parm(J),J=1,nParms)
	do iParm=1,nParms
		WRITE(*,'(8F10.3)')parm(iParm),stdErr(iParm)
	enddo
      STOP
1000  FORMAT(5X,31H FINAL L2 NORM OF THE RESIDUALS,D15.7 //
     &       5X,15H EXIT PARAMETER,16X,I10 //
     &       5X,27H FINAL APPROXIMATE SOLUTION // 5X,3D15.7)
C
C  LAST CARD OF DRIVER FOR LMDIF1 EXAMPLE.
      END
      
      SUBROUTINE FCN(M,N,X,FVEC,IFLAG) !don't forget to make FCN external in main
      INTEGER M,N,IFLAG
      DOUBLE PRECISION X(N),FVEC(M)
C
C  SUBROUTINE FCN FOR LMDIF1 EXAMPLE.
C  THE ANSWERS ARE: 0.0824,1.1330,2.3437 => FNORM=.0906
C
      INTEGER I
      DOUBLE PRECISION TMP1,TMP2,TMP3
      DOUBLE PRECISION Y(15)
      DATA Y/1.4D-1,1.8D-1,2.2D-1,2.5D-1,2.9D-1,3.2D-1,3.5D-1,3.9D-1,
     &       3.7D-1,5.8D-1,7.3D-1,9.6D-1,1.34D0,2.1D0,4.39D0/
      DO 10 I=1,15
        TMP1=I
        TMP2=16-I
        TMP3=TMP1
        IF(I.GT.8)TMP3=TMP2
        FVEC(I)=Y(I)-(X(1)+TMP1/(X(2)*TMP2+X(3)*TMP3))
10    CONTINUE
      RETURN
C
C  LAST CARD OF SUBROUTINE FCN.
C
      END



      SUBROUTINE DevalSsq(nData,nParms,parm,deviate,iFlag)
	implicit doublePrecision(A-H,O-Z)
	parameter(maxData=111,maxParms=20)
      dimension parm(nParms),deviate(nData) !,xDat(nParms),yDat(nData)
	common/data/xDat(maxData),yDat(maxData)
C
C  SUBROUTINE FCN FOR LMDIF1 EXAMPLE.
C  THE ANSWERS ARE: 0.0824,1.1330,2.3437 => FNORM=.0906
C
c	DATA X/.05,.08,.1,.12,.15,.17,.2,.23,.25,.28,.3,.33,.35,.38,.4/
c	DATA Y/0.5073,.8433,1.08,1.3285,1.7209,1.992,2.412,2.84815,
c	1 3.145,3.5921,3.89094,4.34246,4.6326,5.06,5.3309475/
C-23.491x3 + 20.358x2 + 8.9679x
C      DATA Y/1.4D-1,1.8D-1,2.2D-1,2.5D-1,2.9D-1,3.2D-1,3.5D-1,3.9D-1,
C     &       3.7D-1,5.8D-1,7.3D-1,9.6D-1,1.34D0,2.1D0,4.39D0/
c	LMDifEZ gives with unweighted deviates:
c	9.411	13.657	5.638	-38.289 	0.058	0.065	0.595	1.317
c	Poly0 Coeffs:						StdErr:
c	9.4099	13.6615	5.6290	-38.2874	0.06  	0.07  	0.60  	1.32
c	Excel gives with unweighted deviates:
c	Poly0 Coeffs:						StdErr:
c	9.4099	13.6615	5.6290	-38.2874	0.0545	0.6231	2.2380	2.5313
	!dummy=iflag
	rmsErr=0
      DO 10 I=1,nData
        TMP1=I
        TMP2=16-I
        TMP3=TMP1
        IF(I.GT.8)TMP3=TMP2
C        deviate(I)=Y(I)-(parm(1)+TMP1/(parm(2)*TMP2+parm(3)*TMP3))
		poly=parm(1)*xDAT(i)
		if(nParms.gt.1)then
			do iParm=2,nParms
				poly=poly+parm(iParm)*xDAT(i)**(iParm)
			enddo
		endif
        deviate(I)=YDAT(I)-poly
	  rmsErr=rmsErr+deviate(i)*deviate(i)
10    CONTINUE
	rmsErr=sqrt(rmsErr/nData)
	write(*,'(e11.4,4(f11.5))')rmsErr,(parm(i),i=1,nParms)
      RETURN
C
C  LAST CARD OF SUBROUTINE FCN.
C
      END
