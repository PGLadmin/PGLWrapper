!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PrTcParms
	USE GlobConst
	USE BIPs
	DoublePrecision zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX)
	DoublePrecision TCj(NMX),PCj(NMX),acenj(NMX)
	DoublePrecision OMA,OMB	! these are not adjustable parms but it seems silly to recompute them all the time.
END MODULE PrTcParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetPrTc(nComps,iErrCode)
	!  
	!  PURPOSE:  LOOKS UP THE PrTc PARAMETERS AND STORES THEM IN COMMON
	!  Reference:   Jaubert et al., JCED,63, 3980-3988 (2018)
	!  Background: TC stands for "translated consistent." it is a volume-translated PR EOS that uses the Twu alpha function with constraints to make the parameters "consistent."
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	!    commons: ppData, BIPs through GetBips
	!  Programmed by:  JRE 07/00
	USE PrTcParms ! GlobConst(bVol)+alphaL-M+cvol
	implicit doublePrecision( A-H,K,O-Z)
	character*99 bipFile,inFile,dumString
	integer GetBIPs
	LOGICAL LOUDER
	LOUDER=LOUD
	LOUDER=.TRUE.
	NC=nComps
	!common/ParmsPrTc/zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX),OMA,OMB
    iErrCode=0
	etaMax=1-zeroTol
	! note:  bips are passed back through module /BIPs/
		inFile=TRIM(PGLinputDir)//'\ParmsPrTcJaubert.txt' ! // is the concatenation operator
		bipFile=TRIM(PGLinputDir)//'\BipPRtc.txt' ! // is the concatenation operator
	if(LOUD)write(*,*)'GetPrTc: inFile= ',TRIM(inFile)
	open(40,file=inFile,ioStat=ioErr)
	if(ioErr.and.LOUDER)then
		print*,'GetPrTc: PGLinputDir= ',TRIM(PGLinputDir) 
		print*,'GetPrTc: Error opening: ', TRIM(inFile)
		if(LOUDER)pause 'Check that the parms file is in the right dir.'
	endif

	THIRD=1.D0/3
	sqrt2=DSQRT(2.D0)
	sqrt8=2*sqrt2
	etac = 1+(4-sqrt8)**third+(4+sqrt8)**third	 !Jaubert, Eq.7
	etac = 1/etac
	OMA = (40*etac+8)/(49-37*etac)
	OMB = etac/(etac+3)						  
	if(LOUDER)write(*,*)'    ID          TwuL      TwuM      TwuN       cVolCc_mol    Tmin(K) '
	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup and only once if nComps=1.
	do iComp=1,nComps
		iGotIt=0
		rewind(UNIT=40,ioStat=ioErr)
		read(40,*,ioStat=ioErr)nDeck
		!print*,'nDeck=',nDeck
		jComp=0 
		DO while (jComp < nDeck .and. iGotIt==0) !rewind and loop through hard drive till you find the comp of interest.
			jComp=jComp+1
			read(40,'(A99)',ioStat=ioErr)dumString
			!Print*,TRIM(dumString)
			if(ioErr)print*,'GetPrTc: i,j,ioErr,String: ',iComp,jComp,ioErr,TRIM(dumString)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!  NOTE!  Jaubert replaces Tc,Pc,acen with his values!!!      !!!!!!!!!!!!!!!!!
			read(dumString,*,ioStat=ioErr)idBase,Tcj(iComp),PcBar,acenj(iComp),alphaL(iComp),alphaM(iComp),alphaN(iComp),cVolCc_mol(iComp),zRa(iComp),TminK(iComp)
			!if( TminK(iComp) < 50.and.LOUD) print*,TRIM(dumString)
			!if( TminK(iComp) < 50.and.LOUD)print*,'GetPrTc: Warning 50>Tmin=',TminK(iComp)
			if(cVolCc_mol(iComp)==86)cVolCc_mol(iComp)=0
			!cVolCc_mol(iComp)=0 !for debugging, but it shouldn't matter for VLE.
			Pcj(iComp)=PcBar !/10 !JRE changed the database to MPa JRE 20200504
			!read(40,*,ERR=861)idBase,idCc,(zRefDb(iCoeff),iCoeff=1,3),(a1Db(iCoeff),iCoeff=1,nTptCoeffs),(a2Db(iCoeff),iCoeff=1,nTptCoeffs),vMolecDb,tKmin(iComp),nTypes(iComp),(idType(iComp,iType),iType=1,nTypes(iComp)),(nFg(iComp,iType),iType=1,nTypes(iComp))
			IF(idBase==id(iComp))THEN
				iGotIt=1  !this will kick us to the next component
				bVolCC_mol(iComp)=OMB*rGas*TCj(iComp)/PCj(iComp)-cVolCc_mol(iComp)
				if(LOUDER)write(*,'(i7,5f13.4)')ID(iComp),alphaL(iComp),alphaM(iComp),alphaN(iComp),cVolCc_mol(iComp),TminK(iComp)
				if(LOUDER)write(*,'(a,5f13.4)')' Jaubert Tc,Pc,acen,bVol= ', Tcj(iComp),Pcj(iComp),acenj(iComp),bVolCC_mol(iComp)
				if(bVolCC_mol(iComp) < 1)then
					iErrCode=8686
					if(LOUDER)pause 'GetPrTc: bVol < 1??? Could happen if Pc[=]bar or cVol too big.'
				endif 	  
				exit !quit searching if found
			ENDIF !idBase==id(iComp)
		enddo !while(iGotIt.eq.0)
		if(iGotIt==0)then
			iErrCode=iErrCode*10+iComp
			if(LOUDER)write(*,*)'GetPrTc: no data for id#:',id(iComp),' in file=',TRIM(inFile)
		endif
	enddo !all iComps


	IERRGET=GetBIPs(bipFile,ID,NC)
	if(LOUDER)then
		print*,'From GetPrTc:'
		write(*,'(a,11i7)')' kij ',(ID(j),j=1,NC)
		do i=1,NC
			do j=1,NC
				if(kij(i,j) > 0.9)kij(i,j)=0.9d0
				if(kij(i,j) <  -1)kij(i,j)= -1
			enddo
			write(*,'(i5,11f7.4)')ID(j),(kij(i,j),j=1,NC)
		enddo
	endif ! LOUDER, print bips.

	iErrCode=iErrCode*10 !+iErrGet

	RETURN
END	! Subroutine GetPrTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QueryParPurePrTc(iComp,iParm,value,iErr)
	USE PrTcParms      !Just for ESD
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	iErr=0
	if(iParm==1)then
		value=alphaL(iComp)
	elseif(iParm==2)then
		value=alphaM(iComp)
	elseif(iParm==3)then
		value=alphaN(iComp)
	elseif(iParm==4)then
		value=cVolCC_mol(iComp)
	elseif(iParm==5)then
		value=bVolCC_mol(iComp)
	else
		iErr=1
	endif
	return
end	!Subroutine QueryParPurePrTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPurePrTc(iComp,iParm,value,iErr)
	USE PrTcParms      !Just for PrTc
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	iErr=0
	if(iParm==1)then
		alphaL(iComp)=value
	elseif(iParm==2)then
		alphaM(iComp)=value
	elseif(iParm==3)then
		alphaN(iComp)=value
	elseif(iParm==4)then
		cVolCC_mol(iComp)=value
	elseif(iParm==5)then
		bVolCC_mol(iComp)=value
	else
		iErr=1
	endif
	return
end	!Subroutine SetParPurePrTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!$ FuPrTcVtot
	!
	!   REVISION DATE:  7/2020, implemented Jaubert's "tcPR" modifications starting from original PR(1976) code. 
	!   REVISION DATE:  OCTOBER 31, 1985
	!										April 6, 2001, convert to f90
	!   PROGRAMMED BY:  J.R. ELLIOTT, JR. (JAN. 1983)
	!
	!   PURPOSE:  CALCULATE THE FUGACITY COEFFICIENTS AND COMPRESSIBILITY
	!      FACTOR OF EITHER A GAS OR LIQUID ACCORDING TO THE
	!      PENG-ROBINSON EQN OF STATE.
	!
	!   ARGUMENTS:
	!
	!     INPUT:
	!        TC       VECTOR CRITICAL TEMPERATURES OF THE COMPONENTS
	!        PC       VECTOR CRITICAL PRESSURES OF THE COMPONENTS
	!        ACEN     VECTOR ACENTRIC FACTORS OF THE COMPONENTS
	!        RGAS     GAS CONSTANT ( EG. 8.31434 CC-MPA/(GMOL-K) ) IN PHASE LIQ
	!        tAbs     ABSOLUTE TEMPERATURE
	!        pAbs     ABSOLUTE PRESSURE
	!        xFrac    VECTOR MOLE FRACTIONS OF COMPONENTS IN PHASE LIQ
	!        NC       NUMBER OF COMPONENTS
	!        LIQ      PARAMETER SPECIFYING DESIRED PHASE TO BE CONSIDERED
	!                 0 FOR VAPOR PHASE
	!                 1 FOR LIQUID PHASE
	!
	!     OUTPUT:
	!        FUGC     log FUGACITY COEFFICIENT OF COMPONENTS
	!        zFactor  COMPRESSIBILITY FACTOR OF MIXTURE
	!        IER      VECTOR ERROR PARAMETERS
	!          IER(1) = 1 IF ONE OF IER(4)-IER(6) ARE NOT ZERO
	!          IER(2) = 1 IF LIQUID ROOT WAS FOUND BUT WAS NOT REAL
	!          IER(3) = 1 IF VAPOR  ROOT WAS FOUND BUT WAS NOT REAL
	!          IER(4) = 1 NEGATIVE LOG CALCULATED
	!          IER(5) = 1 THE LOG OF THE FUGACITY COEFFICIENT CALCULATED
	!                     TO BE LARGE ENOUGH TO INDICATE OVERFLOW.
	!          IER(6) = 1 SRKNR DID NOT CONVERGE
	!
	!   NOTE:           UNITS OF ALL THE INPUTS SHOULD BE
	!                   CONSISTENT WITH UNITS OF RGAS.  EXCEPT
	!                   FOR THIS, THE USER MAY CHOOSE HIS OWN UNITS.
	!
	!   REQD. ROUTINES:
	!         THE USER MUST SUPPLY A ROUTINE FOR CALCULATING BINARY
	!         INTERACTION COEFFICIENTS OF THE FORM:
	!         FUNCTION ACTION(T,I,J)
	!         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	!         ACTION=
	!         RETURN
	!         END
	!         WHERE T IS THE TEMPERATURE
	!               I,J ARE THE SUBSCRIPTS OF THE BINARY PAIR
	!         ACTION MAY BE SET TO ZERO FOR HYDROCARBON-HYDROCARBON MIXTURES
	!         FOR OTHER MIXTURES, ACTION SHOULD BE ESTIMATED AS
	!         RECOMMENDED IN PROCEDURE 8D1.1.
	!
	!   ROUTINES CALLED:
	!         ACTION,SRKNR
	!
	!   SUBPROGRAM RESTRICTIONS:
	!         AS WRITTEN, THE MAXIMUM NUMBER OF COMPONENTS
	!         THAT CAN BE CONSIDERED IS TEN.  THIS RESTRICTION CAN BE
	!         OVERCOME BY CHANGING THE DIMENSION STATEMENTS.
	!
	!   GENERAL COMMENTS:
	!         BRIEFLY, THE ALGORITHM CALLS FOR THE CALCULATION
	!         OF THE EXTREMUM(Z**3-Z**2+Q*Z-R) WHEN EITHER A LIQUID
	!         OR VAPOR ROOT IS DETERMINED NOT TO BE REAL, EVEN THOUGH
	!         THE CALCULATION CALLS FOR A ROOT.  IF SUCH ROOTS
	!         ARE FOUND, IER SIGNIFIES THIS.  WHEN IT IS DETERMINED
	!         THAT THE DESIRED ROOT EXISTS, IT IS CALCULATED USING
	!         NEWTON-RAPHSON ITERATION IN THE ROUTINE "SRKNR".
	!
	!   REFERENCES:
	!		  J.R. Elliott, C.T. Lira, Introductory Chemical Engineering Thermodynamics, p 334 (1999). 2ed p262 (2012).
	!		
	!         PENG, D-Y, ROBINSON, D.B. IEC Fun,15: 59–64
	!
	!         SOAVE, G., CHEM. ENG. SCI., 27:1197 (1972).
	!
	!         Cismondi: J. Chem. Eng. Data 2019, 64, 2093-2109
	!         Jaubert:FPE, 429 (2016) 301-312
	!
	!****************************************************************
	SUBROUTINE FuPrTcVtot(isZiter,tKelvin,vTotCc,gMol,NC,FUGC,zFactor,aDep,uDep,iErrZ)
	USE PrTcParms ! GlobConst(bVol)+alphaL-M+cvol
	USE BIPs
	Implicit DoublePrecision(A-H,K,O-Z)
	DIMENSION FUGC(NC),xFrac(NC),gMol(NC) !,IER(*)
	DoublePrecision ALA(NMX,NMX),bVol(NMX),TDLAL_DT(NMX),T2d2LAL_dT2(NMX),NdC_b_dni
	!common/ParmsPrTc/zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX),OMA,OMB
	COMMON/DEPFUN/dU_NKT,dA_NKT,dS_NK,dH_NKT
	COMMON/eta/etaL,etaV,zFactorL,zFactorV
	data initCall/1/
	!  prBIPs are passed in from GetPrBIPs()
	iErrZ=0
	!if(zFactor < BIGB)iErrZ=13
	!IF( (zFactor+(1+sqrt2)*BIGB) < 0) iErrZ=14
	!IF(dChemPo > 33) iErrZ=15
	!if(iErrZ=6)'FuPrTcVtot: 0~cmprsblty=',cmprsblty
	!if(iErrTmin.ne.0)iErrZ=5
	!iErrZ=17  if dChemPo > 33
	sqrt2=DSQRT(2.D0)
	sqrt8=2*sqrt2
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	if(LOUD.and.initCall)print*,'FuPrTcVtot:OMA,OMB',OMA,OMB
	totMol=sum(gMol(1:NC))
	if( vTotCc < 1E-11 .and. LOUD)print*,'FuPrTcVtot: 0 ~ vTotCc=',vTotCc
	rhoMol_cc=totMol/vTotCc
	bMix=0
	iErrTmin=0
	TcVolatile=1234
	do iComp = 1,NC
		if(PCj(iComp)==0)then
			iErrZ=11
			if(LOUD)write(*,'(a,F7.2,5F8.3)')' FuPrTcVtot: nonsense input. Tc,Pc,L,M,N,c: ',TCj(iComp),PCj(iComp),alphaL(iComp),alphaM(iComp),alphaN(iComp),cVolCc_mol(iComp),TminK(iComp)
			cycle
		endif
		if(Tcj(iComp) < TcVolatile)then
			TcVolatile=Tcj(iComp)
			iVolatile=iComp
		endif
		aCrit = OMA*RGAS*RGAS*TCj(iComp)*TCj(iComp)/PCj(iComp)
		!sTmp = 0.37464 + 1.54226*ACEN(iComp) - 0.26993*ACEN(iComp)*ACEN(iComp)
		!DLALDT(iComp)= -sTmp*SQRT(Tr/ ALPHA) ! dAlp/dT = 2*sqrt(alp)*S*(-0.5/sqrt(Tr)) = -sqrt(alp)*S/sqrt(Tr); (T/alp)*dAlp/dT= -S*sqrt(Tr/Alp)	 EL2ed Eq.7.18,8.35
		Tr = tKelvin / Tcj(iComp)
		Pow = alphaM(iComp)*alphaN(iComp)
		TrPow = Tr**Pow 
		!ALPHA = TrPow/ Tr**alphaN(iComp) *DEXP( alphaL(iComp)*(1-TrPow) )	 ! = Tr^[N*(M-1)]*exp[ L*(1-Tr^MN) ] Eq. 4
		rLnAlpha=(Pow-alphaN(iComp))*DLOG(Tr) + alphaL(iComp)*(1-TrPow)	 ! -> -L*TrPow as Tr -> inf
		if( (rLnAlpha) < -33)rLnAlpha= -33.d0 !avoid exponential underflow
		ALPHA = EXP(rLnAlpha)  ! -> 0
		!ALPHA = (  1 + ( 0.37464D0 + acen(iComp)*(1.54226D0 - 0.26993D0*ACEN(iComp)) )*( 1 - DSQRT( Tr) )  )**2	 ! for debugging
		!ln(ALPHA) = pow*ln(Tr) - N*ln(Tr) + L*(1-Tr^pow)
		! dLn(ALPHA)/dTr = pow/Tr - N/Tr - L*pow*Tr^(pow-1) 
		! d2Ln(alpha)/dTr^2 = -pow/Tr^2 + N/Tr^2 -L*pow*(pow-1)*Tr^(pow-2)
		!(Tr/Alp)*dAlp/dTr = Tr*dLn(alpah)//dTr = Pow - N - L*pow*Tr^Pow   -> -L*pow*Tr^Pow at Tr->inf
		! T^2*d2(ln(alpha)/dT^2 = -pow + N -L*pow*(pow-1)*Tr^pow

		TdLAL_dT(iComp)= Pow - alphaN(iComp) - alphaL(iComp)*pow*TrPow	 ! Ures/RT= -Ares/RT+(a/bRT)*TdLAL_dT*ln(...) = -Ares/RT + Ares/RT*TdLAL_dT
		!              = alphaN_i*[alphaM_i-1-alphaL_i*alphaM_i*Tr^(alphaM*alphaN)]
		! Td/dT[ TdLALdT ] = T*d/dT[ Pow - alphaN(iComp) - pow*alphaL(iComp)*TrPow ]
		T2d2LAL_dT2(iComp) = -pow+alphaN(iComp) -pow*(pow-1)*alphaL(iComp)*TrPow
		ALA(iComp,iComp) = aCrit*ALPHA
		!bVol(iComp) = OMB*RGAS*TC(iComp)/PC(iComp) - cVolCc_mol(iComp)	
		bVol(iComp) = bVolCc_mol(iComp)	
		xFrac(iComp)=gMol(iComp)/totMol
		bMix=bMix+xFrac(iComp)*bVol(iComp)
		FUGC(iComp)=0 !initial to zero for isZiter=1.
	enddo
	if( tKelvin < TminK(iVolatile) )iErrTmin=1
	if(iErrTmin.ne.0)then
		if(LOUD.and.initCall)write(*,'(a,5f8.2)' )' FuVtot: T,Tmin(i)',tKelvin,( TminK(i),i=1,NC)
		iErrZ=5
	endif
	if(iErrZ>10)return
	eta=bMix*rhoMol_cc
	if(LOUD.and.initCall)print*,'FuPrTcVtot: eta,bMix',eta,bMix
	if(LOUD.and.initCall)write(*,'(a,3(1PE11.4))')' FuPrTcVtot:TdLAL_dT,T2d2LAL_dT2',TdLAL_dT(1),T2d2LAL_dT2(1)
	if(ABS(totMol-1) > 1e-5)then
		if(LOUD)print*, 'FuPrTcVtot warning: totMol= ',totMol
		!iErrZ=2
		!goto 861
	endif

	aMix = 0
	TdaMixDt=0
	T2d2aMixDt2=0
	cMix = 0
	DO iComp = 1,NC
		DO jComp = 1,NC
			ALA(iComp,jComp) = SQRT( ALA(iComp,iComp)*ALA( jComp,jComp) )*(1-KIJ(iComp,jComp))
			aMix = aMix + ALA(iComp,jComp)*xFrac(iComp)*xFrac( jComp)
			TdALA_dT= ALA(iComp,jComp)*( TdLAL_dT(iComp)+TdLAL_dT( jComp) )/2	 ! ~ [aMix]
			TdTdALA_dT_dT =     TdALA_dT      * ( TdLAL_dT(iComp) + TdLAL_dT( jComp) )/2 &	   ! ~ [aMix* dAmix/aMix] ~ [aMix]
			                  +  ALA(iComp,jComp)*( T2d2LAL_dT2(iComp)+T2d2LAL_dT2( jComp) )/2 ! where d2LALdT2= d/dT[ DLALDT ]		 ! ~ [aMix]
			TdaMixDt = TdaMixDt + TdALA_dT*xFrac(iComp)*xFrac( jComp)
			T2d2aMixDt2 = T2d2aMixDt2 + TdTdALA_dT_dT*xFrac(iComp)*xFrac( jComp)
		enddo
		cMix = cMix + cVolCc_mol(iComp)*xFrac(iComp) ! dCmix_dNi=c(i)
	enddo
	!TdaMixDt = TdaMixDt*aMix
 	!T2d2aMixDt2 = T2d2aMixDt2*aMix
	!BIGMES=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!dU_NKT= BIGMES*(-1+daMixDt/aMix)	!EL2ed Eq.8.35
	if(aMix < zeroTol .and. LOUD)print*,'FuPrTcVtot:0~aMix,aCrit,alpha,Tr=',aMix,aCrit,alpha,Tr
	if(bMix < zeroTol .and. LOUD)print*,'FuPrTcVtot:bMix~0=',bMix

	eta=bMix*rhoMol_Cc
	crho=cMix*rhoMol_Cc
	!print*,'FuPrTcVtot: eta,bMix,cMix=',eta,bMix !,cMix


!	Z=1/(1-brho) - (A/Bsqrt8)*brho/[(1+crho)*(1+brho+2crho)+(brho+crho)*(1-brho)]
!	(1+crho)*(1+brho+2crho)+(brho+crho)*(1-brho)=1+brho+2crho+crho+brho*crho+2crho^2
!                                                                          +brho+crho-brho^2-brho*crho
!                                                                          = 1+2*brho-brho^2 +crho*(4-brho+2crho)
!                                                                         = 1 + eta*bq + eta^2*cq; 
!                                                                         bq = 2+4c/b; cq = -[1-2(c/b)^2]; -qq=bq^2-4cq
!                                                                        Aatt =(1/sqrt(-qq))*{ ln[ (2cq*eta+bq-sqrt(-qq))/(2cq*eta+bq+sqrt(-qq) ] - ln[ (bq-sqrt(-qq))/(bq+sqrt(-qq) ] } EL2ed Eq. B.34
!  Alt: RKPR(Cismondi)=>Aatt=(a/bRT)*(1/(d1-d2))*y/[ (1+d1y)(1+d2y) ] ; 1+di*y = 0 => di = -1/yi; yi = [ -bq +/- sqrt(bq^2 - 4cq) ]/2; 
!        denom = [ (1+d1y)(1+d2y) ] = 1 + (d1+d2) y + d1d2 y^2 => bq = (d1+d2) = C+1; -cq = -d1*d2 = 1-2(c/b)^2 = C    ; d2 = cq/d1 => d1 + cq/d1 = bq => d1^2 + cq -bqd1 = 0 => d1 = ( bq +/- sqrt(bq^2-4*cq) )/2	= (bq +/- sqrtNqq)/2
	c_b = cMix/bMix
	bq=2+4*c_b
	cq= 2*(c_b)*(c_b)-1
	! qq = -(bq^2-4cq) ; sqrtNqq=sqrt(-qq) = (d1-d2)=(1+c/b)*sqrt8
	sqrtNqq=(1+c_b)*sqrt8 ! x = c/b & -qq = (bq*bq-4*cq) = (2+4x)^2-4(2x^2-1)= 4+16x+16x^2-8x^2+4 = 8+16x+8x^2 = 8*(1+x)^2  
	voidFrac=1-eta
	if( voidFrac.le. 0 .and. LOUD)print*,'FuPrTcVtot:0~ voidFrac = ',voidFrac
	dArep_dEta=1/ voidFrac
	zRep=eta * dArep_dEta
	d1=1+sqrt2+c_b*(2+sqrt2) !=( bq+sqrtNqq )/2=[(2+4x)+(1+x)sqrt8]/2 = 1+2x+(1+x)*sqrt2 = 1+sqrt2+x*(2+sqrt2) ! = 2cq/(bq-sqrtNqq)
	d2=1-sqrt2+c_b*(2-sqrt2) !=( bq-sqrtNqq )/2=[(2+4x)-(1+x)sqrt8]/2 = 1+2x-(1+x)*sqrt2 = 1-sqrt2+x*(2-sqrt2) ! = 2cq/(bq+sqrtNqq)
	!  2	  +c_b*2          = d1+d2 = 2*(1+x)
	denom = 1+eta*bq + cq*eta*eta ! (1+d1*eta)*(1+d2*eta) = (1+(d2+d1)*eta+d1*d2*eta^2); 
	dAatt_dEta= -aMix/(bMix*RGAS*tKelvin)/denom	 ! Jaubert Eq.7
	zAtt=eta*dAatt_dEta
	zFactor = 1+zRep+zAtt 
	pMPa=zFactor*rhoMol_Cc*RGAS*tKelvin

	!print*,'qq,sqrt(-qq)',qq,sqrtNqq
	!print*,'b,c/b',bMix,c_b
	!print*,'bq,cq',bq,cq
	!print*,'zRep,zAtt',zRep,zAtt
	!print*,'d1+,d2+',d1Plus,d2Plus
	!print*,'d1-,d2-',d1Minus,d2Minus
	BIGA = aMix*pMPa/(RGAS*RGAS*tKelvin*tKelvin)
	BIGB = bMix*pMPa/(RGAS*tKelvin)
	if(LOUD.and.initCall)print*,'FuPrTcVtot: zFactor,pMPa=',zFactor,pMPa
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
	!
	!BIGMES=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!argLog1 = (2*cq*eta+bq-sqrtNqq)/(2*cq*eta+bq+sqrtNqq)	! argLog1 may be < 0, but then so is argLog2
	!argLog2 = (         bq-sqrtNqq     )/(         bq+sqrtNqq     )
	argLog1 = 1+d1*eta
	argLog2 = 1+d2*eta

	argLog = argLog1/argLog2
	if(LOUD.and.initCall)print*,'FuPrTcVtot: argLog1,2=',argLog1,argLog2
	if(loud .and. argLog < zeroTol)print*, 'FuPrTcVtot: argLog 1or2 < 0',argLog1,argLog2
	aResAtt= -BIGA/BIGB/sqrtNqq*DLOG( argLog )
	!print*,'a/bRT,alpha', BIGA/BIGB,alpha
	!print*,'zAtt,aAtt',zAtt,aResAtt
	aAttRKPR = -BIGA/BIGB/(d1-d2)*DLOG( (1+d1*eta)/(1+d2*eta) )
	zAttRKPR = -BIGA/BIGB*eta/( (1+d1*eta)*(1+d2*eta) )
	denomRKPR= (1+d1*eta)*(1+d2*eta)
	!print*,'RKPR: zAtt,aAtt',zAttRKPR, aAttRKPR
	!print*,'denom: tcPR,RKPR', denom,denomRKPR
	!pause 'check values'
	!BIGMESOld=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!print*,'BigMes,BigMesOld:',BigMes,BigMesOld
	aRes_RT= -LOG(1-eta) + aResAtt
	!uDep/RT = beta*d(A/RT)/dBeta = -T*d(A/RT)/dT = -T* [ -A/RT^2 + (1/RT)*dA/dT ] = A/RT - (dA/dT)/R = A/RT*[ 1-(T/alpha)*dAlpha/dT ]  
	uRes_RT= aResAtt*(1-TdaMixDt/aMix)	! Cv = (dU/dT) = U/T + T*d(U/T)/dT; U/RT = A/RT*(1-TdaMixDt/aMix) => Td(U/RT)/dT = T(dA/RT)/dT -(A/RT)*Td/dT[TdaMixDt/aMix] = -U/RT -A/RT*[ Td(TdaMixDt)/dT -(TdaMixDt/aMix)^2] 
	!uDep_R = aAtt*T*(1-TdLnAlpha/dT) ~ (a/bR)
	!d(uDep_R)/dT = 
	! Td(Ures/RT)/dT = bigmes*[ TdTdaMix_dT_dT/aMix - TdaMixDt^2/aMix^2 ] 
	CvRes_R=  aResAtt*( -T2d2aMixDt2 - TdaMixDt*TdaMixDt/aMix )/aMix
	CvRes_R=  aResAtt*( -TdLAL_dT(1)**2 - T2d2LAL_dT2(1) )	 !todo: work on checking the generaliztion to mixtures. Above is slightly off. 
	!T*dZ/dT = T*dZatt/dT = eta*T*d/dT[dAatt/dEta]	; dAatt/dEta= -(aMix/bRT)/denom
	! T*d/dT[dAatt/dEta]= -T*dLnaMix/dT*(a/bRT)/denom -dAatt_dEta = (TdaMixDt/aMix)*dAatt_dEta - dAatt_dEta	= dAatt_dEta*(TdaMixDt/aMix - 1); NOTE: the -dAatt_dEta comes from the 1/T in a/bRT.
	dA_NKT=aRes_RT
	dU_NKT=uRes_RT	 !for passing through GlobConst
	TdZ_dT= zAtt*(1-TdaMixDt/aMix) !since aAtt = f1(T)*f2(rho), the derivative wrt T in dZ/dT is analogous to derivative wrt T in dAatt/dT. 
	dH_NKT=dU_NKT+zFactor-1
	dS_NK=dU_NKT-dA_NKT !+LOG(zFactor)	  !Don't subtract lnZ here. it may cause trouble if Z < 0. 
	hRes_RT=dH_NKT
	sRes_R  =dS_NK

	aDep=dA_NKT	!to pass as calling argument as well as common
	uDep=dU_NKT

	!zRep=eta/(1-eta) = eta*dArep/dEta
	!d2ARep/dEta2 = 1/(1-eta)^2
	!dZrep/dEta = dArep/dEta+eta*d2Arep/dEta2 
	!eta*dZrep/dEta = zRep + eta^2*d2Arep/dEta2 = eta/(1-eta) + eta^2/(1-eta)^2 = eta/(1-eta)*[ 1+eta/(1-eta) ] = eta/(1-eta)^2
	!zAtt= -eta*aMix/(bMix*RGAS*tKelvin)/( 1+eta*bq + cq*eta*eta ) = eta*dAatt_dEta	 ! Jaubert Eq.7
	!dAatt_dEta= -aMix/(bMix*RGAS*tKelvin)/denom = -(a/bRT)/denom	 ! Jaubert Eq.7
	!d2Aatt_dEta2 = +(a/bRT)/( 1+eta*bq + cq*eta*eta )^2 * (bq + 2*cq*eta ) 
	!d2Aatt_dEta2 = -dAatt_dEta/( 1+eta*bq + cq*eta*eta ) * (bq + 2*cq*eta ) 
	d2Aatt_dEta2 = -dAatt_dEta/denom * (bq + 2*cq*eta )  
	!dZAtt/dEta = dAatt_dEta + eta*d2Aatt_dEta2 => eta*dZAtt/dEta = zAtt + eta^2*d2Aatt_dEta2
	! Z = P/rhoRT => dZ/dRho = (P/RT)(-1/rho^2) + 1/(rhoRT)*dP/dRho; rho*dZ/dRho = (P/RT)(-1/rho) + 1/(RT)*dP/dRho; 
	!cmprsblty=(dP/dRho)T*(1/RT) = Z+rho*dZ/dRho 
	dZ_dEta = 1/( voidFrac*voidFrac ) + dAatt_dEta + eta*d2Aatt_dEta2
	cmprsblty=zFactor+eta*dZ_dEta
	!cmprsblty=cmprsblty*1.804 ! dunno why, but works for CH4.
	!PGL6edEq.6.24=> CpRes_R = CvRes_R-1-[Z+TdZ/dT]^2/[Z+rho*dZ/dRho]
	if( ABS(cmprsblty) > 1.D-11)then
		!if(LOUD)print*,'(dP/dT)/rhoR=',(zFactor+TdZ_dT)
		CpRes_R = CvRes_R-1 + (zFactor+TdZ_dT)*(zFactor+TdZ_dT)/cmprsblty
	else
		if(LOUD)print*,'FuPrTcVtot: 0~cmprsblty=',cmprsblty
		iErrZ=6
		CpRes_R=86.8686D0
	endif

	if(isZiter.ne.0)goto 861  ! It's good to know cmprsblty even if isZiter==1. e.g. slope of spinodal, critPure, ...

	if(zFactor < BIGB)then
		iErrZ=13
		goto 861
	endif
	if(LOUD.and.initCall)print*,'aRes,uRes',aRes_RT,uRes_RT
	DO iComp = 1,NC
		SUMXA = 0
		DO jComp=1,NC
			SUMXA = SUMXA + xFrac( jComp)*ALA(iComp,jComp)
		enddo

		IF( (zFactor+(1+sqrt2)*BIGB) < 0 .or. zFactor < BIGB) THEN
			!	AVOID CALCULATION OF NEGATIVE LOGARITHMS BUT INDICATE ERROR.
			ChemPoRes = 33
			iErrZ=14
		ELSE
			NdC_b_dni= ( cVolCC_mol(iComp) - cMix*bVolCC_mol(iComp)/bMix )/bMix	!NOTE: declared DoublePrecision.
			!BIGB=b*P/RT  ; c_b = cMix/bMix
			ChemPoRes = bVol(iComp)/bMix*(zFactor-1) - DLOG(zFactor-BIGB) & 
			!+ aResAtt*( 2*SUMXA/aMix - bVol(ICOMP)/bMix ) + zFactor*(cMix-cVolCC_mol(iComp))*rhoMol_cc	! Privat(2016a) Eq 14.
			+ aResAtt*( 2*SUMXA/aMix - bVol(ICOMP)/bMix - NdC_b_dni/(1+c_b) ) &
			+ zAtt*(1-eta)*NdC_b_dni/(1+c_b)  
		END IF
		IF(ChemPoRes > 33) THEN
			!	AVOID EXPONENT OVERFLOW. DON'T INDICATE ERROR. IT MAY BE N-HEXANE AT P>1GPa AND 298K.
			iErrZ=7
			ChemPoRes = 33
		endif
		FUGC(ICOMP) = (ChemPoRes)
	enddo
	if(LOUD.and.initCall)print*,'FUGC(1-NC)',(fugc(i),i=1,nc)
!	if(LIQ.eq.0)then
!	  etaV=BIGB/zFactor
!	  zFactorV=zFactor
!	ELSE
!	  etaL=BIGB/zFactor
!	  zFactorL=zFactor
!	ENDIF

861	continue
	initCall=0
	RETURN
	END

	!SUBROUTINE FugiPrTc(tKelvin,pMPa,gMol,NC,LIQ,FUGC,zFactor,IER)
	Subroutine FugiPrTc( tKelvin,pMPa,gmol,NC,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,ier )	  ! JRE 2020
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!$ FUGI
	!
	!   REVISION DATE:  OCTOBER 31, 1985
	!										April 6, 2001, convert to f90
	!   PROGRAMMED BY:  J.R. ELLIOTT, JR. (JAN. 1983)
	!
	!   PURPOSE:  CALCULATE THE FUGACITY COEFFICIENTS AND COMPRESSIBILITY
	!      FACTOR OF EITHER A GAS OR LIQUID ACCORDING TO THE
	!      PENG-ROBINSON EQN OF STATE.
	!
	!   ARGUMENTS:
	!
	!     INPUT:
	!        TC       VECTOR CRITICAL TEMPERATURES OF THE COMPONENTS
	!        PC       VECTOR CRITICAL PRESSURES OF THE COMPONENTS
	!        ACEN     VECTOR ACENTRIC FACTORS OF THE COMPONENTS
	!        RGAS     GAS CONSTANT ( EG. 8.31434 CC-MPA/(GMOL-K) ) IN PHASE LIQ
	!        tAbs     ABSOLUTE TEMPERATURE
	!        pAbs     ABSOLUTE PRESSURE
	!        xFrac    VECTOR MOLE FRACTIONS OF COMPONENTS IN PHASE LIQ
	!        NC       NUMBER OF COMPONENTS
	!        LIQ      PARAMETER SPECIFYING DESIRED PHASE TO BE CONSIDERED
	!                 0 FOR VAPOR PHASE
	!                 1 FOR LIQUID PHASE
	!
	!     OUTPUT:
	!        FUGC     log FUGACITY COEFFICIENT OF COMPONENTS
	!        zFactor  COMPRESSIBILITY FACTOR OF MIXTURE
	!        IER      VECTOR ERROR PARAMETERS
	!          IER(1) = 1 IF ONE OF IER(4)-IER(6) ARE NOT ZERO
	!          IER(2) = 1 IF LIQUID ROOT WAS FOUND BUT WAS NOT REAL
	!          IER(3) = 1 IF VAPOR  ROOT WAS FOUND BUT WAS NOT REAL
	!          IER(4) = 1 NEGATIVE LOG CALCULATED
	!          IER(5) = 1 THE LOG OF THE FUGACITY COEFFICIENT CALCULATED
	!                     TO BE LARGE ENOUGH TO INDICATE OVERFLOW.
	!          IER(6) = 1 SRKNR DID NOT CONVERGE
	!
	!   NOTE:           UNITS OF ALL THE INPUTS SHOULD BE
	!                   CONSISTENT WITH UNITS OF RGAS.  EXCEPT
	!                   FOR THIS, THE USER MAY CHOOSE HIS OWN UNITS.
	!
	!   REQD. ROUTINES:
	!         THE USER MUST SUPPLY A ROUTINE FOR CALCULATING BINARY
	!         INTERACTION COEFFICIENTS OF THE FORM:
	!         FUNCTION ACTION(T,I,J)
	!         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	!         ACTION=
	!         RETURN
	!         END
	!         WHERE T IS THE TEMPERATURE
	!               I,J ARE THE SUBSCRIPTS OF THE BINARY PAIR
	!         ACTION MAY BE SET TO ZERO FOR HYDROCARBON-HYDROCARBON MIXTURES
	!         FOR OTHER MIXTURES, ACTION SHOULD BE ESTIMATED AS
	!         RECOMMENDED IN PROCEDURE 8D1.1.
	!
	!   ROUTINES CALLED:
	!         ACTION,SRKNR
	!
	!   SUBPROGRAM RESTRICTIONS:
	!         AS WRITTEN, THE MAXIMUM NUMBER OF COMPONENTS
	!         THAT CAN BE CONSIDERED IS TEN.  THIS RESTRICTION CAN BE
	!         OVERCOME BY CHANGING THE DIMENSION STATEMENTS.
	!
	!   GENERAL COMMENTS:
	!         BRIEFLY, THE ALGORITHM CALLS FOR THE CALCULATION
	!         OF THE EXTREMUM(Z**3-Z**2+Q*Z-R) WHEN EITHER A LIQUID
	!         OR VAPOR ROOT IS DETERMINED NOT TO BE REAL, EVEN THOUGH
	!         THE CALCULATION CALLS FOR A ROOT.  IF SUCH ROOTS
	!         ARE FOUND, IER SIGNIFIES THIS.  WHEN IT IS DETERMINED
	!         THAT THE DESIRED ROOT EXISTS, IT IS CALCULATED USING
	!         NEWTON-RAPHSON ITERATION IN THE ROUTINE "SRKNR".
	!
	!   REFERENCES:
	!		      J.R. Elliott, C.T. Lira, Introductory Chemical Engineering Thermodynamics, p 334 (1999).
	!		
	!         GUNDERSEN, T., COMPUTERS IN CHEM. ENG., 6:245 (1982).
	!
	!         SOAVE, G., CHEM. ENG. SCI., 27:1197 (1972).
	!
	!         PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!****************************************************************
	USE PrTcParms ! GlobConst(bVol)+alphaL-M+cvol
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	DIMENSION FUGC(NC),xFrac(NC),gMol(NC),IER(12)
	!DIMENSION bVol(NMX),DLALDT(NMX), ALA(NMX,NMX),
	COMMON/DEPFUN/dU_NKT,dA_NKT,dS_NK,dH_NKT
	COMMON/eta/etaL,etaV,zFactorL,zFactorV
	!common/ParmsPrTc/zRa(NMX),cVolCc_mol(NMX),alphaL(NMX),alphaM(NMX),alphaN(NMX),TminK(NMX),OMA,OMB
	data initCall/1/
	!  prBIPs are passed in from GetPrBIPs()
	IER=0 !sets all array values to zero.
    iErr=0 !error tracker internal for FUGI()

	!zeroTol=1D-11 !zeroTol is now global
	if(tKelvin < zeroTol .and.LOUD)pause 'FuPrTc: Nonsense on input. 0~tKelvin'
	totMoles=sum(gMol) 
	xFrac=gMol/totMoles
	if( ABS(totMoles-1) > zeroTol .and. LOUD)print*,'FuPrTc: ??? totMoles,x(1)=',totMoles,xFrac(1) 
	bMix=SumProduct(NC,xFrac,bVolCc_mol)
	if(LOUD.and.initCall)print*,'FuPrTc: bMix=',bMix
	aRes_RT=86.8686	  ! initialize to avoid NAN on error return
	uRes_RT=86.8686
	cvRes_R=86.8686
	cpRes_R=86.8686
	cmprsblty=86.8686

	isZiter=1

	pb_RT = pMPa*bMix/(rGas*tKelvin)
	eta=pb_RT/1.001d0  	!super high pressures can generate eta>1 because Z>>1. Check eta in if() below.
	if(LIQ==1 .or. LIQ==3 .or. eta>etaMax)eta=etaMax/1.001d0 !organize so "new" eta for 1st iteration will be exactly ideal gas value, to improve precision when P~1E-11.
	rhoMol_Cc=eta/bMix
	if(LOUD.and.initCall)write(*,'( a,3(1PE11.4) )')' FuPrTc: 1st FuVtot. rhoMol_cc,eta=', rhoMol_cc,eta 
	call FuPrTcVtot(isZiter,tKelvin,1/rhoMol_Cc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
	!pause 'check initial values'
	if(iErrZ > 10)then
		iErr=13
		if(LOUD)then
            write(*,*)'FuPrTc: Initial FuVtot returned with error code:',iErrZ
		    pause
        end if
		goto 861
	elseif(iErrZ.ne.0)then
		iErr=iErrZ
	endif
	etaOld=eta
	errOld=pb_rt-eta*zFactor
	eta=eta/1.001d0
	!ETA=rhoMol_Cc*bVolMix	!calculate here before entering loop, then at the end with each new rhoMol_Cc
	change=1234
	itMax=55
	NITER=0
	pMax= -1.d11 !for crude goldenz
	pMin=1.d11
	isZiter=1
	do while(ABS(change) > 1.e-10 .and. iErr < 10)   ! This gives ~10 sig figs on rho for liq and above initialization gives exact ideal gas value when P->1E-11.
		NITER=NITER+1
		if(niter > itMax)iErr=16
		if(niter > itMax .and. LOUD)print*,'FugiPrTc: iterations exceeded. eta=',eta
		if(niter > itMax)exit
		!ETA=rhoMol_Cc*bVolMix
		if( (1-eta) < 1.d-11)then
			if(LOUD)write(*,*) 'FuPrTc: eta > 1. nIter=',nIter
			if(LOUD)pause
		endif
		rhoMol_cc=eta/bMix
		vTotCc=totMoles/rhoMol_cc
		if(LOUD.and.initCall)print*,'FuPrTc:',NITER,'th FuVtot. eta=', eta 
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
		if(iErrZ>10)then
			iErr=13
			if(LOUD)write(*,*)'FuPrTc: FuVtot returned with error code:',iErrZ
			!call BeepMsg(errMsg(iErr))
			cycle
		elseif(iErrZ.ne.0)then
			iErr=iErrZ
		endif
		error=pb_rt-eta*zFactor
		if(eta*zFactor > pMax .and. eta < 0.15)then
			pMax=eta*zFactor
			etaAtPmax=eta	 !for crude goldenz. ~randomly searches for lo eta max in p
		endif
		if(eta*zFactor < pMin .and. eta > 0.15)then
			pMin=eta*zFactor
			etaAtPmin=eta	 !for crude goldenz. ~randomly searches for hi eta min in p
		endif
		change=error/(error-errOld)*(eta-etaOld)
		if(LOUD.and.initCall)write(*,'(a,2F8.4,2E11.4)')' eta,Z,error,change',eta,zFactor,error,change
		!pause 'check error'
		etaOld=eta
		errOld=error
        if(zFactor < 0)change=1 !force another iteration if Z < 0. This happens when P ~ 1E-12.
		if(ABS(change/eta) > 0.1)change=eta*DSIGN(0.1d0,change)	!step limiting. 0.3>0.2 to give max looseness
		!if(ABS(change/zFactor) > 1)change=zFactor*DSIGN(0.5d0,change)	!step limiting. Stop zFactor from going negative.
		if(LOUD.and.initCall)print*,'Fugi:eta,change',eta,change
		eta=eta-change
		rhoMol_cc=eta/bMix
		if(eta < 0 .or. eta > etaMax)then !NOTE: this should not happen given step limiting
			if(niter < (itMax-1))niter=itMax-1 !restrict tries with new guess so we can crash if necessary.
			if(LOUD)write(*,'(a,f8.4)')' Warning in FuTpt: next guess is eta=',eta
			eta=0
			if(liq==1)eta=etaMax
			if(LOUD)write(*,*) 'Restarting iteration with eta=',eta
			if(LOUD)pause
		endif
	enddo  !iteration on eta to find P=P(input)
	!print*,'Hello! iteration concluded. iErr,eta=',iErr,eta
	etaPass=eta											! Accurate rho is essential for PsatEar at low T (e.g. propane).	
	if(iErr > 10)then ! use crude goldenz to keep iterations going. 
		!print*,'Hello! iErr.ne.0 . LIQ,etaAtPmax,etaAtPmin=',LIQ,etaAtPmax,etaAtPmin
		if(liq==0 .or. LIQ==2)eta=etaAtPmax  !crude goldenz
		if(liq==1 .or. LIQ==3)eta=etaAtPmin  !crude goldenz
		if(eta < zeroTol)eta=zeroTol
		!print*,'Hello! iErr.ne.0. eta,bMix=',eta,bMix
		rhoMol_cc=eta/bMix
		vTotCc=totMoles/rhoMol_cc
		if(LOUD)write(*,'(a,i3,1x,f7.4,a,f10.5)')' FuPrTc: iErr,goldenEta=',iErr,eta
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
		if(LOUD)write(*,'(a,f10.5)')' FuPrTc:  using goldenZ=',zFactor
		iErr=10
	endif
	IF (eta < 0)then
		iErr=12
		!call BeepMsg(errMsg(iErr))
	endif
	IF (zFactor < 0 )then
		iErr=14
		!IF(LOUD)call BeepMsg(errMsg(iErr))
	endif
	!  ITERATION ON RHO HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.
	rhoMol_cc=eta/bMix
	vTotCc=totMoles/rhoMol_cc
	IF (LIQ==1 .or. LIQ==3) THEN
	  ETAL=ETA
	  ZL=zFactor  
	ELSE
	  ETAV=ETA
	  ZV=zFactor  
	ENDIF
	dU_RT=uDep
	dH_RT=uDep+zFactor-1
	dA_RT=aDep
	aRes_RT=aDep
	uRes_RT=uDep
	hRes_RT=dH_RT
    if(LIQ > 1)goto 861
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS and derivative props (cmprsblty,CvRes_R,CpRes_R passed by GlobConst) 
	!
	!print*,'aRes,uDep',aDep,uDep
	Sres_R=uDep-aDep					  ! A = U - TS => Sres_R = Ures/RT - Ares/RT
	if(Loud.and.initCall)print*,'Sres,hRes',Sres_R,hRes_RT
	isZiter=0
	call FuPrTcVtot(isZiter,tKelvin,1/rhoMol_Cc,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrZ)
	if(LOUD.and.initCall)print*,'aRes,uDep',aRes,uRes
	if(iErrZ>10)then
		iErr=13
	elseif(iErrZ.ne.0)then
		iErr=iErrZ
	endif  	
	zFactor=pMPa/(rhoMol_cc*rGas*tKelvin) ! Seemingly unnecessary, this should improve precision when computing rho from return. Z = 1+zRef+zAtt+zAssoc is subject to roundoff when Z->1E-9, but Z=P/(rhoRT)=>rho=P/(ZRT) with precision.
	etaPass=eta											! Accurate rho is essential for PsatEar at low T (e.g. propane).	
861	continue
	if(iErrZ.and.LOUD)print*,'FugiPrTc: FuVtot last call. iErrZ=',iErrZ
	initCall=0
	IER(1)=iErr
	RETURN
	END

	SUBROUTINE ZITERDum(zFactor,BIGB,A1,A0,IER)	 !ZITER() is defined in FugiPR. 
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!
	!  PURPOSE    - CALCULATE zFactor FROM NEWTON-RAPHSON ITERATION
	!
	!      IER    - 200 DIDNT CONVERGE IN 25 ITERATIONS
	!
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit doublePrecision(A-H,O-Z)
	IER = 0
	do kount=1,25
		F = -A0 + zFactor*(  A1+zFactor*( -(1-BIGB)+zFactor )  )
		DF = 3*zFactor*zFactor - (1-BIGB)*2*zFactor + A1
		ZN = zFactor - F/DF
		ERR = DABS((ZN - zFactor)/ZN)
		zFactor = ZN
		IF(ERR.lt. 1.D-9) GO TO 86
	enddo
	IER = 200 !only to reach here is if do loop has exceeded iterations
86 continue
	RETURN
	END

