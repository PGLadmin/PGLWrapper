
      Subroutine GetPRWS(NC,iErrCode)
C  
C  PURPOSE:  LOOKS UP THE PRWS PARAMETERS AND STORES THEM IN COMMON
C            BINARY INTERACTION PARAMETERS ARE STORED IN /BIPs/ AND
C            STRYJEK-VERA PARAMETERS ARE STORED IN /PRSVWS/.
C
C  INPUT
C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
C  OUTPUT
	USE GlobConst
	USE VpDb ! for VpCoeffs
      IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
      PARAMETER(ndb=350,numSV=20)
	character bipFile*234
	integer GetBIPs
c	dimension IDDippr(numSV),svKap1(numSV)
      COMMON /PRSVWS/ svKappa1(NMX)
C           Stryjek and Vera, can j chem eng, 64:323 (1986).

	write(*,*)'ID  NAME       TCK   PCMPa      w     '
	do 8 i=1,nc
	  write(*,606) ID(i),NAME(i),TC(i),PC(i),ACEN(i)
8	continue
606	format(i4,1x,a11,f6.1,f6.3,1x,f6.3,f8.2,f6.1,f8.1,i4,f8.6,f8.4)
	GetPRSVWS=0
	IF(DEBUG)THEN
		bipFile='c:\spead\CalcEos\input\ParmsPrws.txt'
	ELSE
		bipFile=TRIM(masterDir)//'\input\ParmsPrws.txt' ! // is the concatenation operator
	ENDIF
	open(55,file=bipFile)
	do iComp=1,NC
		svKappa1(iComp)=0
c		idArray(iComp)=idDippr((iComp))
		read(55,*)numKappa
		do iKappa=1,numKappa
			read(55,*)idKappa,svKap
			if( idKappa.eq.ID(iComp) )svKappa1(iComp)=svKap
		enddo
		rewind(55)
	enddo
	close(55)
c  note:  bips are passed back through common/BIPs/
	IF(DEBUG)THEN
		bipFile='c:\spead\CalcEos\input\BipWongSand.txt'
	ELSE
		bipFile=TRIM(masterDir)//'\input\BipWongSand.txt' ! // is the concatenation operator
	ENDIF
      IERRCODE=GetBIPs(bipFile,ID,NC)

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc PROGRAMED BY AV 06/26/06cccccccccccccccccccccccccccccccccccccccccccccc
cccccPURPOSE: CALCULATION OF THERMODYNAMIC PROPERTIES WITH NRTLccccccccc

	SUBROUTINE GetNRTL (nComps,idComp,iErrCode)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=1555)
	character bipFile*234,inFile*234
	integer GetBIPs
	DIMENSION idComp(nComps)
	dimension solParmD(ndb),vLiqD(ndb),IDA(ndb)
	!common/NRTL/vLiq(nmx)
	write(*,*)'ID  NAME       TCK   PCMPa      w     '
	iErrCode=0
c  note:  bips are passed back through common/BIPs/
	IF(DEBUG)then 
		OPEN(50,FILE='c:\SPEAD\CalcEos\input\ParmsCrit.dta',
	1	FORM='BINARY')
      ELSE 
		inFile=TRIM(masterDir)//'\input\ParmsCrit.dta' ! // is the concatenation operator
		OPEN(50,FILE=inFile,FORM='BINARY')
	ENDIF

	READ(50,ERR=861)NDECK 
	DO I=1,NDECK
	READ(50,ERR=861)IDA(I),tcd,pcd,zcd,acend,rmwd,solParmD(i),vLiqD(i)
	enddo   
	CLOSE(50)

	IF(DEBUG)then 
		OPEN(50,FILE='c:\SPEAD\CalcEos\input\parmsCrAdd.txt')
      ELSE 
		inFile=TRIM(masterDir)//'\input\parmsCrAdd.txt' ! // is the concatenation operator
		OPEN(50,FILE=inFile)
	ENDIF
	read(50,*,ERR=863)NDECKADD
	nTot=nDeck+nDeckAdd
	DO I=nDeck+1,nTot
		READ(50,*,ERR=863)IDA(I),tcd,pcd,zcd,acend,rmwd,solParmD(i),
	1	vLiqD(i)
	enddo   
	CLOSE(50)

	DO J=1,nComps
		DO I=1,NDECK
			IF(IDA(I).EQ.idComp(J))THEN
			  vLiq(j)=vLiqD(i)
			ENDIF
		enddo
	enddo
		
	IF(DEBUG)THEN
		bipFile='c:\spead\CalcEos\input\BipWongSand.txt'
	ELSE
		bipFile=TRIM(masterDir)//'\input\BipWongSand.txt' ! // is the concatenation operator
	ENDIF
      IERRCODE=GetBIPs(bipFile,idComp,nComps)

      RETURN
861	continue
	write(*,*)'GetNRTL error - error reading ParmsCrit.txt. Path?'
	write(*,*)'nDeck,iCompo',NDECK,i
	if(LOUD)pause
	return                      
863	continue
	write(*,*)'NRTL error - error reading ParmsCrAdd.txt. Path?'
	write(*,*)'nDeckAdd,iCompo',nDeckAdd,i
	if( (i-1).gt.0)write(*,*)'last good parms',ida(i),solParmD(i),
	1vLiqD(i)
	if(LOUD)pause
	return                      
	END

      

	!PROGRAMED BY AV 06/22
	!PURPOSE: CALCULATION OF VAPOR PRESSURES AND FUGACITY COEFFICIENTS BY NRTL ACTIVITY COEFFICIENT MODEL

	subroutine FuNRTL(tK,pMpa,moFrac,nComps,LIQ,FUGC,zFactor,IER)
c  use NRTL style of EXPression
	USE GlobConst
	USE BIPs
	USE VpDb ! for VpCoeffs
	implicit doubleprecision(A-H,K,O-Z)
	double precision moFrac(Nmx)
	DIMENSION IER(12)
	integer kComp
	dimension fugc(nComps),xsGamma(nComps)
	DIMENSION omega(Nmx,Nmx),sumDenom(Nmx),sumNumer(Nmx)
	!COMMON/NRTL/vLiq(nmx) 
	!DIMENSION vpCoeffs(NMX,5)
	!dimension pSat(NMX)
  	do iErr=1,6
		ier(iErr)=0
	enddo
	!DEBUG=.TRUE.
	vLiqMix=0
	do iComp = 1,nComps
		vLiqMix=vLiqMix+moFrac(iComp)*vLiq(iComp)
	enddo
	zFactor=1
	fugc(iComp)=1
	if(LIQ.eq.1)then
		zFactor=pMpa*vLiqMix/(rGas*tK)

		sumLog=0.0
		do 10 jComp=1,nComps
			sumDenom(jComp)=0.0
			sumNumer(jComp)=0.0
			do 9 kComp=1,nComps
				omega(kComp,jComp)=EXP(-xsAlpha(kComp,jComp)
	1		  		                    *xsTau(kComp,jComp)/tK)
				sumDenom(jComp)=sumDenom(jComp)+
	1			moFrac(kComp)*omega(kComp,jComp)
				sumNumer(jComp)=sumNumer(jComp)+moFrac(kComp)
	1			 *omega(kComp,jComp)*xsTau(kComp,jComp)/tK
9			continue
			sumLog=sumLog+moFrac(jComp)*sumDenom(jComp)
10		continue

		xsFreeEn =0
		do 20 kComp=1,nComps
			bigSumK=0.d0
			do 15 jComp=1,nComps
				bigSumK=bigSumK
	1			 +moFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)
     2		*( xsTau(kComp,jComp)/tK-sumNumer(jComp)/sumDenom(jComp) )
15			continue
			xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
			xsFreeEn =xsFreeEn+xsLnGamma*moFrac(kComp)
			fugc(kComp)=xsLnGamma
			xsGamma(kComp)=exp(xsLnGamma)
20		continue
		do kComp=1,nComps
			!pSat=pc(kComp)*10**( 7*(1+acen(kComp))/3*(1-tc(kComp)/
	1		!tK) )
			!if(pSat.lt.1e-33)pSat=1E-33 !logArg in fugc() cannot be zero
  			!IF(kComp.eq.1)pSat=27491.07E-6
			!if(kComp.eq.nComps)pSat=15145.42E-6 

			!if(tDat(1).eq.tDat(2).AND.XDAT(NPTS).eq.1.and.kComp.eq.1)then
				!pSat=PDAT(NPTS)
			!elseif(tDat(1).eq.tDat(2).AND.XDAT(1).eq.0.and.kComp.eq.2)then
			!	pSat=PDAT(1)
			!else
						!pSat=pc(kComp)*10**( 7*(1+acen(kComp))/3*
	1					!(1-tc(kComp)/tK) )
	!THE ABOVE ESTIMATION OF vp WILL CAUSE NRTL TO GOVE BAD AAD IN vp. THE FOLLOWING EQN FROM DIPPR DATABASE 
	!IS RECOMMENDED
				        CALL GetVp(nComps,ID,iErrCode)
			            pSat=exp(vpCoeffs(kComp,1)+vpCoeffs(kComp,2)/tK
	1					+vpCoeffs(kComp,3)*DLOG(tK)
     1                    +vpCoeffs(kComp,4)*tK**vpCoeffs(kComp,5))
     1                    /1000000
					!endif
				!enddo
			!endif		
			fugc(kComp)=fugC(kComp)+LOG(pSat/pMPa)
26		  enddo
	   endif
	   return
	   end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FugiPRWS(T,P,X,NC,LIQ,FUGC,Z,IER)
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C      PURPOSE --- THIS ROUTINE CALCULATES THE FUGACITIES AND
C                  USING THE PENG-ROBINSON EQUATION WITH STRYJEK-VERA
C                  CORRECTION OF THE PR ALPHA AND WONG-SANDLER MIXING
C                  RULE.  WHERE POSSIBLE, IT FOLLOWS THE CODE FROM NIST14
C                  FROM THE SUMMER OF 99.
c
C      CODED BY:	JRE 7/99 INC. WS MRULE
C		ref Wong and Sandler, aicheJ, 38:677 (1992).
C
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
	USE GlobConst
	USE BIPs
      IMPLICIT DOUBLEPRECISION (A-H,K,O-Z)
	DIMENSION IER(12)
      DIMENSION bigQSum(Nmx),xsChemPo(Nmx)
	dimension X(Nmx),ZR(3),AZ(3),A(Nmx),B(Nmx),fugc(Nmx)
C
	dimension b2KIJ(NMX,NMX),b2KTIJ(NMX,NMX)
c      COMMON/PRSVWS/ svKappa1(NMX)
      COMMON/PRSVWS/ svKappa1(Nmx)
c	data initEos/0/
c  cf. stryjek and vera, can j chem eng., 64:323 (1986). 
      DATA  kappa0,kappa1,kappa2,kappa3
	1  /0.378893d0,1.4897153,-0.17131848d0,0.0196554d0/
	b2KIJ=KIJ
	b2KTIJ=KTIJ
	do iErr=1,6
		ier(iErr)=0
	enddo
C
	if(initEos.eq.0)then
	  initEos=1
C
C                        SET UP PRS PARAMETERS
	   aOmega=0.45723553d0
	   bOmega=0.07779607d0
	   root2=SQRT(2.d0)
	   con1=  1+root2
	   con2=-(1-root2)
	   con3=2*root2
	   RDUM=RGAS
c	   rGas=R
	   bigC = 1/root2*LOG(root2-1)  !WS eq A4 !
      ENDIF
      DO 030 I = 1, NC
		PR = P / PC(I)
		TR = T / TC(I)
	svKappa0=kappa0+acen(I)*( kappa1+acen(I)*(kappa2+acen(I)*kappa3) )
	    kappa=svKappa0+svKappa1(i)*(1+SQRT(TR))*(0.7d0-TR) 
		ALPHA = 1.D0+kappa*(1.D0-SQRT(TR))
c  note:  A(i) = SQRT(a*p/(RT)^2)
c			B(i) = bP/(RT)
		A(I) = ALPHA * SQRT(AOMEGA*PR) / TR
		B(I) = BOMEGA * PR / TR
030   CONTINUE
C                        GIVEN COMPOSITIONS OF LIQUID AND VAPOR,
C                        CALCULATE THE FUGACITIES AND K-VALUES
	call xsPRWS(xsChemPo,xsFreeEn,X,T,NC)
	AMX = 0.0
	BMX = 0.0
	aObMx = 0.0
	bigQ=0.0
C                        CALCULATE MIXTURE EOS PARAMETERS
	DO 080 I = 1, NC
		aObI = A(I)*A(I)/B(I)
		b2i=( 1-aObI )*B(I)
		aObMx = aObMx + X(I) * aObI
		bigQSum(I) = 0.0
		DO 060 J = 1, NC
			aObJ = A(J)*A(J)/B(J)
			b2j=( 1-aObJ )*B(J)
c  eq 11
			b2ij=(b2i+b2j)/2*(1 -(b2Kij(I,J)+b2KTij(I,J)/T) )
			bigQSum(I) = bigQSum(I) + X(J)*b2ij
			bigQ=bigQ+X(J)*X(I)*b2ij
  060		CONTINUE
  080	CONTINUE
c  eq A9 => bigD = (a/bRT)mix
c  eq A8 => bigQ = (b2)mix *P/RT  the *P/RT is different but keeps dimensionless.
	bigD=aObMx+xsFreeEn/bigC
	BMX=bigQ/(1-bigD)
c  noting A/B = a/bRT = bigD, then
	AMX=bigD*BMX
C                        SET UP AND SOLVE THE CUBIC EQUATION
C
	AZ(1) = (-AMX+BMX+BMX*BMX)*BMX
	AZ(2) = AMX-BMX*(3.0*BMX+2.0)
	AZ(3) = BMX-1.0
	CALL CUBIC(AZ,ZR)
	Z = ZR(2)
	if(liq.eq.1)Z=ZR(1)
C                        CALCULATE FUGACITIES
C
	QUEST=Z-BMX
	IF(QUEST.LT.0)IER(1)=1
	TMP1 = (Z-1.0) / BMX
	TMP2 = -LOG(MAX(Z-BMX,1.0D-20))
	TMP3 = -AMX * LOG((Z+CON1*BMX)/(Z-CON2*BMX)) / (CON3*BMX)
	DO 100 I = 1, NC
c		original pr values:
c			dnbdni=B(i)
c			dn2adni=2*ASUM(I)
c  note:  (a/bRT)i=Ai*Ai/Bi
	  dnDdni=A(I)*A(I)/B(I)+xsChemPo(I)/bigC
	  dnbdni=2*bigQSum(I)/(1-bigD)-bigQ*(1-dnDdni)/(1-bigD)**2
	  dn2adni=bigD*dnbdni+BMX*dnDdni
		fugc(I)=dnbdni*TMP1+TMP2+TMP3*(dn2adni/AMX-dnbdni/BMX)
100	CONTINUE
      RETURN
      END

	subroutine xsPRWS(fugXS,xsFreeEn,moFrac,tK,nComps)
c  compute the excess mixture value by the PR Wong-Sandler mixture rule
c  use NRTL style of EXPression
	USE BIPs
	implicit doubleprecision(A-H,K,O-Z)
	double precision moFrac(Nmx)
	integer kComp
	dimension fugXS(nComps),xsGamma(nComps)
	DIMENSION omega(Nmx,Nmx),sumDenom(Nmx),sumNumer(Nmx)
      common/PRSVWS/svKappa1(Nmx)

c  begin the computation of the xsGamma.  could use NRTL or anything else here.
c  ref.  ReID et al. 1987, PGL p 274,256
c  sumDenom is the sumK(OmegaKI*xK), 
c  OmegaKI = Gki = EXP(-alphaKI*wsTau/T)
	sumLog=0.0
	do 10 jComp=1,nComps
	  sumDenom(jComp)=0.0
	  sumNumer(jComp)=0.0
	  do 9 kComp=1,nComps
			omega(kComp,jComp)=EXP(-xsAlpha(kComp,jComp)
	1		  		                    *xsTau(kComp,jComp)/tK)
			sumDenom(jComp)=sumDenom(jComp)+
	1     moFrac(kComp)*omega(kComp,jComp)
			sumNumer(jComp)=sumNumer(jComp)+moFrac(kComp)
	1     *omega(kComp,jComp)*xsTau(kComp,jComp)/tK
9	  continue
	  sumLog=sumLog+moFrac(jComp)*sumDenom(jComp)
10	continue

	xsFreeEn =0
	do 20 kComp=1,nComps
		bigSumK=0.d0
	  do 15 jComp=1,nComps
	    bigSumK=bigSumK
	1	  +moFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)
     2    *( xsTau(kComp,jComp)/tK-sumNumer(jComp)/sumDenom(jComp) )
15		continue
		xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
	  xsFreeEn =xsFreeEn+xsLnGamma*moFrac(kComp)
		fugXS(kComp)=xsLnGamma
		xsGamma(kComp)=exp(xsLnGamma)
20	continue

	return
	end
C
