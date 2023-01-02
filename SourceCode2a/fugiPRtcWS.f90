
      Subroutine GetPRtcWS(NC,iErrCode)
!C  
!C  PURPOSE:  LOOKS UP THE PRWS PARAMETERS AND STORES THEM IN COMMON
!C            BINARY INTERACTION PARAMETERS ARE STORED IN /BIPs/ AND
!C            STRYJEK-VERA PARAMETERS ARE STORED IN /PRSVWS/.
!C
!C  INPUT
!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
!C  OUTPUT
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=350,numSV=20)
	character bipFile*234
	integer GetBIPs
!c	dimension IDDippr(numSV),svKap1(numSV)
	common/ParmsPrTc/zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX),OMA,OMB
!C           This version of PRWS uses Jaubert's correlation for alpha. 

	write(*,*)'ID  NAME       TCK   PCMPa      w     '
	do i=1,nc
	  write(*,606) ID(i),NAME(i),TC(i),PC(i),ACEN(i)
	enddo
606	format(i4,1x,a11,f6.1,f6.3,1x,f6.3,f8.2,f6.1,f8.1,i4,f8.6,f8.4)
	GetPRSVWS=0
	iErrCode=0
	Call GetPRtc(NC,iErrGet) !use the existing code to loade common parms

	if(iErrGet.ne.0)then
		iErrCode=1
		if(LOUD)print*,'GetPRtcWs: failed to load pure parms.'
		return
	endif
!c  note:  bips are passed back through common/BIPs/
	IF(DEBUG)THEN
		bipFile='c:\spead\CalcEos\input\BipPRtcWS.txt'
	ELSE
		bipFile=TRIM(masterDir)//'\input\BipPRtcWS.txt' ! // is the concatenation operator
	ENDIF
	iErrGet=GetBIPs(bipFile,ID,NC)
	if(iErrGet.ne.0)then
		iErrCode=2
		if(LOUD)print*,'GetPRtcWs: failed to load BIPs.'
		return
	endif

	RETURN
	END

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FugiPRtcWS( T,P,gMol,NC,LIQ,FUGC,Z,IER)
!C
!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!C
!C      PURPOSE --- THIS ROUTINE CALCULATES THE FUGACITIES AND
!C                  USING THE PENG-ROBINSON EQUATION WITH STRYJEK-VERA
!C                  CORRECTION OF THE PR ALPHA AND WONG-SANDLER MIXING
!!C                  RULE.  WHERE POSSIBLE, IT FOLLOWS THE CODE FROM NIST14
!C                  FROM THE SUMMER OF 99.
!c
!C      CODED BY:	JRE 7/99 INC. WS MRULE
!C		ref Wong and Sandler, aicheJ, 38:677 (1992).
!C
!C
!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION ( A-H,K,O-Z)
	DIMENSION IER(*)
	DIMENSION bigQSum(Nmx),xsChemPo(Nmx)
	dimension X(Nmx),ZR(3),AZ(3),A(Nmx),B(Nmx),fugc(Nmx)
	dimension ALA(NMX,NMX),TdLAL_dT(NMX),T2d2LAL_dT2(NMX),gMol(NMX),bVol(NMX)
!C
	dimension b2KIJ(NMX,NMX),b2KTIJ(NMX,NMX)
!c      COMMON/PRSVWS/ svKappa1(NMX)
	common/ParmsPrTc/zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX),OMA,OMB
	data initCall/1/
	b2KIJ=KIJ
	b2KTIJ=KTIJ
	do iErr=1,6
		ier(iErr)=0
	enddo
!C                        SET UP PRS PARAMETERS
	tAbs=T
	pMPa=P
	if(initCall==1)then
		initEos=1
		root2=SQRT(2.d0)
		con1=  1+root2
		con2=-(1-root2)
		con3=2*root2
		bigC = 1/root2*LOG(root2-1)  !WS eq A4 !
	ENDIF
	totMol=sum(gMol)
	DO iComp = 1, NC
		aCrit = OMA*RGAS*RGAS*TC(iComp)*TC(iComp)/PC(iComp)
		X(i)=gMol(i)/totMol
		!sTmp = 0.37464 + 1.54226*ACEN(iComp) - 0.26993*ACEN(iComp)*ACEN(iComp)
		!DLALDT(iComp)= -sTmp*SQRT(Tr/ ALPHA) ! dAlp/dT = 2*sqrt(alp)*S*(-0.5/sqrt(Tr)) = -sqrt(alp)*S/sqrt(Tr); (T/alp)*dAlp/dT= -S*sqrt(Tr/Alp)	 EL2ed Eq.7.18,8.35
		Tr = tAbs / Tc(iComp)
		Pow = alphaM(iComp)*alphaN(iComp)
		TrPow = Tr**Pow 
		!ALPHA = TrPow/ Tr**alphaN(iComp) *DEXP( alphaL(iComp)*(1-TrPow) )	 ! = Tr^[N*(M-1)]*exp[ L*(1-Tr^MN) ] Eq. 4
		rLnAlpha=(Pow-alphaN(iComp))*DLOG(Tr) + alphaL(iComp)*(1-TrPow)	 ! -> -L*TrPow as Tr -> inf
		if( (rLnAlpha) < -33)rLnAlpha= -33.d0 !avoid exponential underflow
		ALPHA = EXP(rLnAlpha)  ! -> 0
		!ln(ALPHA) = pow*ln(Tr) - N*ln(Tr) + L*(1-Tr^pow)
		! dLn(ALPHA)/dTr = pow/Tr - N/Tr - L*pow*Tr^(pow-1) 
		! d2Ln(alpha)/dTr^2 = -pow/Tr^2 + N/Tr^2 -L*pow*(pow-1)*Tr^(pow-2)
		!(Tr/Alp)*dAlp/dTr = Tr*dLn(alpah)//dTr = Pow - N - L*pow*Tr^Pow   -> -L*pow*Tr^Pow at Tr->inf
		! T^2*d2(ln(alpha)/dT^2 = -pow + N -L*pow*(pow-1)*Tr^pow

		TdLAL_dT(iComp)= Pow - alphaN(iComp) - alphaL(iComp)*pow*TrPow	 ! Ures/RT= -Ares/RT+(a/bRT)*TdLAL_dT*ln(...) = -Ares/RT + Ares/RT*TdLAL_dT
		! Td/dT[ TdLALdT ] = T*d/dT[ Pow - alphaN(iComp) - pow*alphaL(iComp)*TrPow ]
		T2d2LAL_dT2(iComp) = -pow+alphaN(iComp) -pow*(pow-1)*alphaL(iComp)*TrPow
		ALA(iComp,iComp) = aCrit*ALPHA
		!bVol(iComp) = OMB*RGAS*TC(iComp)/PC(iComp) - cVolCc_mol(iComp)	
		bVol(iComp) = bVolCc_mol(iComp)	
		bMix=bMix+x(iComp)*bVol(iComp)
		FUGC(iComp)=0 !initial to zero for isZiter=1.
		!c  note: A(i) = SQRT(a*p/(RT)^2)
		!c			B(i) = bP/(RT)
		Pr= pMPa/ Pc(iComp)
		A(I) =  SQRT( ALPHA*OMA*Pr ) / Tr
		B(I) = OMB * Pr / Tr
	enddo

030   CONTINUE
!C                        GIVEN COMPOSITIONS OF LIQUID AND VAPOR,
!C                        CALCULATE THE FUGACITIES AND K-VALUES
	call xsPRtcWS(xsChemPo,xsFreeEn,X,T,NC)
	AMX = 0.0
	BMX = 0.0
	aObMx = 0.0
	bigQ=0.0
!C                        CALCULATE MIXTURE EOS PARAMETERS
	DO i = 1, NC
		aObi = A(i)*A(i)/B(i)
		b2i=( 1-aObi )*B(i)
		aObMx = aObMx + X(i) * aObi
		bigQSum(i) = 0.0
		DO J = 1, NC
			aObJ = A(J)*A(J)/B(J)
			b2j=( 1-aObJ )*B(J)
			!c  eq 11
			b2ij=(b2i+b2j)/2*(1 -(b2Kij(I,J)+b2KTij(I,J)/T) )
			bigQSum(I) = bigQSum(I) + X(J)*b2ij
			bigQ=bigQ+X(J)*X(I)*b2ij
		enddo ! j
	enddo ! i
!c  eq A9 => bigD = (a/bRT)mix
!c  eq A8 => bigQ = (b2)mix *P/RT  the *P/RT is different but keeps dimensionless.
	bigD=aObMx+xsFreeEn/bigC
	BMX=bigQ/(1-bigD)
!c  noting A/B = a/bRT = bigD, then
	AMX=bigD*BMX
!C                        SET UP AND SOLVE THE CUBIC EQUATION
!C
	AZ(1) = (-AMX+BMX+BMX*BMX)*BMX
	AZ(2) = AMX-BMX*(3*BMX+2)
	AZ(3) = BMX-1
	CALL CUBIC(AZ,ZR)
	Z = ZR(2)
	if(liq==1 .or. liq==3)Z=ZR(1)
	if(Z > 0)etaPass=BMX/Z
!C                        CALCULATE FUGACITIES
!C
	QUEST=Z-BMX
	IF(QUEST.LT.0)IER(1)=1
	TMP1 = (Z-1.0) / BMX
	TMP2 = -LOG(MAX(Z-BMX,1.0D-20))
	TMP3 = -AMX * LOG((Z+CON1*BMX)/(Z-CON2*BMX)) / (CON3*BMX)
	DO 100 I = 1, NC
!c		original pr values:
!c			dnbdni=B(i)
!c			dn2adni=2*ASUM(I)
!c  note:  (a/bRT)i=Ai*Ai/Bi
	  dnDdni=A(I)*A(I)/B(I)+xsChemPo(I)/bigC
	  dnbdni=2*bigQSum(I)/(1-bigD)-bigQ*(1-dnDdni)/(1-bigD)**2
	  dn2adni=bigD*dnbdni+BMX*dnDdni
		fugc(I)=dnbdni*TMP1+TMP2+TMP3*(dn2adni/AMX-dnbdni/BMX)
100	CONTINUE
    RETURN
    END

	subroutine xsPRtcWS( fugXS,xsFreeEn,moFrac,tK,nComps)
!c  compute the excess mixture value by the PR Wong-Sandler mixture rule
!c  use NRTL style of EXPression
	USE BIPs
	implicit doubleprecision(A-H,K,O-Z)
	double precision moFrac(Nmx)
	integer kComp
	dimension fugXS(nComps),xsGamma(nComps)
	DIMENSION omega(Nmx,Nmx),sumDenom(Nmx),sumNumer(Nmx)

!  begin the computation of the xsGamma.  could use NRTL or anything else here.
!c  ref.  ReID et al. 1987, PGL p 274,256
!c  sumDenom is the sumK(OmegaKI*xK), 
!c  OmegaKI = Gki = EXP(-alphaKI*wsTau/T)
	sumLog=0.0
	do jComp=1,nComps
		sumDenom(jComp)=0.0
		sumNumer(jComp)=0.0
		do kComp=1,nComps
			omega(kComp,jComp)=EXP(-xsAlpha(kComp,jComp)*xsTau(kComp,jComp)/tK)
			sumDenom(jComp)=sumDenom(jComp)+moFrac(kComp)*omega(kComp,jComp)
			sumNumer(jComp)=sumNumer(jComp)+moFrac(kComp)*omega(kComp,jComp)*xsTau(kComp,jComp)/tK
		enddo
		sumLog=sumLog+moFrac(jComp)*sumDenom(jComp)
	enddo

	xsFreeEn =0
	do kComp=1,nComps
		bigSumK=0.d0
		do jComp=1,nComps
			bigSumK=bigSumK+moFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)*( xsTau(kComp,jComp)/tK-sumNumer(jComp)/sumDenom(jComp) )
		enddo ! jComp
		xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
		xsFreeEn =xsFreeEn+xsLnGamma*moFrac(kComp)
		fugXS(kComp)=xsLnGamma
		xsGamma(kComp)=exp(xsLnGamma)
	enddo ! kComp
	return
	end
!C
!C

