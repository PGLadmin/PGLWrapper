!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PrTcParms
	USE GlobConst
	USE BIPs
	USE PREosParms !OMA,OMB, ...
	DoublePrecision zRa(nmx),cVolCc_mol(nmx),alphaL(nmx),alphaM(nmx),alphaN(nmx) !,TminK(nmx)
	DoublePrecision Tcj(nmx),Pcj(nmx),acenj(nmx)
END MODULE PrTcParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetPrTc(nComps,iErrCode)
	!  
	!  PURPOSE:  LOOKS UP THE PrTc PARAMETERS AND STORES THEM IN USEd PREosParms
	!  Reference:   Jaubert et al., JCED,63, 3980-3988 (2018)
	!  Background: TC stands for "translated consistent." a volume-translated PR EOS usin Twu alpha function 
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	!    USEd: GlobConst, BIPs 
	!  Programmed by:  JRE 07/00
	USE PrTcParms ! GlobConst(bVol)+alphaL-M+cvol & PREosParms
	Implicit DoublePrecision( A-H,K,O-Z)
	character*99 bipFile,inFile,dumString
	integer GetBIPs
	LOGICAL LOUDER
	LOUDER=LOUD
	LOUDER=.TRUE.
	NC=nComps
    iErrCode=0
	etaMax=1-zeroTol
	! note:  bips are passed back through module /BIPs/

	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)

		inFile=TRIM(PGLinputDir)//'\ParmsPrTcJaubert.txt' ! // is the concatenation operator
		bipFile=TRIM(PGLinputDir)//'\BipPRtc.txt' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,*)'GetPrTc: inFile= ',TRIM(inFile)
	open(40,file=inFile,ioStat=ioErr)
	if(ioErr.and.LOUDER)then
		write(dumpUnit,*)'GetPrTc: PGLinputDir= ',TRIM(PGLinputDir) 
		write(dumpUnit,*)'GetPrTc: Error opening: ', TRIM(inFile)
		if(LOUDER)write(dumpUnit,*)'Check that the parms file is in the right dir.'
	endif

	if(LOUDER)write(dumpUnit,*)'    ID          TwuL      TwuM      TwuN       cVolCc_mol    Tmin(K) '
	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup and only once if nComps=1.
	do iComp=1,nComps
		iGotIt=0
		rewind(UNIT=40,ioStat=ioErr)
		read(40,*,ioStat=ioErr)nDeck
		!write(dumpUnit,*)'nDeck=',nDeck
		jComp=0 
		DO while (jComp < nDeck .and. iGotIt==0) !rewind and loop through hard drive till you find the comp of interest.
			jComp=jComp+1
			read(40,'(A99)',ioStat=ioErr)dumString
			!write(dumpUnit,*)TRIM(dumString)
			if(ioErr)write(dumpUnit,*)'GetPrTc: i,j,ioErr,String: ',iComp,jComp,ioErr,TRIM(dumString)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!  NOTE!  Jaubert replaces Tc,Pc,acen with his values!!!      !!!!!!!!!!!!!!!!!
			read(dumString,*,ioStat=ioErr)idBase,Tcj(iComp),PcBar,acenj(iComp),alphaL(iComp),alphaM(iComp),alphaN(iComp) &
			                                                                           ,cVolCc_mol(iComp),zRa(iComp),tKmin(iComp)
			!if( TminK(iComp) < 50.and.LOUD) write(dumpUnit,*)TRIM(dumString)
			!if( TminK(iComp) < 50.and.LOUD)write(dumpUnit,*)'GetPrTc: Warning 50>Tmin=',TminK(iComp)
			if(cVolCc_mol(iComp)==86)cVolCc_mol(iComp)=0
			!cVolCc_mol(iComp)=0 !for debugging, but it shouldn't matter for VLE.
			Pcj(iComp)=PcBar !/10 !JRE changed the database to MPa JRE 20200504
			IF(idBase==id(iComp))THEN
				iGotIt=1  !this will kick us to the next component
				bVolCc_mol(iComp)=OMB*Rgas*TCj(iComp)/PCj(iComp)-cVolCc_mol(iComp)
				if(LOUDER)write(dumpUnit,'(i7,5f13.4)')ID(iComp),alphaL(iComp),alphaM(iComp),alphaN(iComp),cVolCc_mol(iComp) &
				                                                                                                     ,tKmin(iComp)
				if(LOUDER)write(dumpUnit,'(a,5f13.4)')' Jaubert Tc,Pc,acen,bVol= ', Tcj(iComp),Pcj(iComp),acenj(iComp), &
				                                                                                                 bVolCc_mol(iComp)
				if(bVolCc_mol(iComp) < 1)then
					iErrCode=8686
					if(LOUDER)write(dumpUnit,*)'GetPrTc: bVol < 1??? Could happen if Pc[=]bar or cVol too big.'
				endif 	  
				exit !quit searching if found
			ENDIF !idBase==id(iComp)
		enddo !while(iGotIt.eq.0)
		if(iGotIt==0)then
			iErrCode=iErrCode*10+iComp
			if(LOUDER)write(dumpUnit,*)'GetPrTc: no data for id#:',id(iComp),' in file=',TRIM(inFile)
		endif
	enddo !all iComps


	IERRGET=GetBIPs(bipFile,ID,NC)
	if(LOUDER)then
		write(dumpUnit,*)'From GetPrTc:'
		write(dumpUnit,'(a,11i7)')' kij ',(ID(j),j=1,NC)
		do i=1,NC
			do j=1,NC
				if(KIJ(i,j) > 0.9)KIJ(i,j)=0.9d0
				if(KIJ(i,j) <  -1)KIJ(i,j)= -1
			enddo
			write(dumpUnit,'(i5,11f7.4)')ID(j),(KIJ(i,j),j=1,NC)
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
		value=cVolCc_mol(iComp)
	elseif(iParm==5)then
		value=bVolCc_mol(iComp)
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
		cVolCc_mol(iComp)=value
	elseif(iParm==5)then
		bVolCc_mol(iComp)=value
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
	!        Rgas     GAS CONSTANT ( EG. 8.31434 CC-MPA/(GMOL-K) ) IN PHASE LIQ
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
	!        iErr = 1-9 warning
	!				12 IF LIQUID ROOT WAS FOUND BUT WAS NOT REAL
	!				13 IF VAPOR  ROOT WAS FOUND BUT WAS NOT REAL
	!				14 NEGATIVE LOG CALCULATED
	!				15 THE LOG OF THE FUGACITY COEFFICIENT CALCULATED TO BE LARGE ENOUGH TO INDICATE OVERFLOW.
	!				16 NR DID NOT CONVERGE
	!
	!   NOTE:           UNITS OF ALL THE INPUTS SHOULD BE
	!                   CONSISTENT WITH UNITS OF Rgas.  EXCEPT
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
	SUBROUTINE FuPrTcVtot(isZiter,tKelvin,vTotCc,gMol,NC,FUGC,zFactor,aRes,uRes,iErrZ)
	USE PrTcParms ! GlobConst(bVol)+alphaL-M+cvol
	USE BIPs
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision FUGC(NC),xFrac(NC),gMol(NC) 
	DoublePrecision ALA(nmx,nmx),TdLAL_dT(nmx),T2d2LAL_dT2(nmx),NdC_b_dni
    DoublePrecision argLog,tKelvin,eta,zFactor !,bVol(nmx)
	data initCall/1/
	!  prBIPs are passed in from GetPrBIPs()
	iErrZ=0
	!if(zFactor < BIGB)iErrZ=13
	!IF( (zFactor+(1+sqrt2)*BIGB) < 0) iErrZ=14
	!IF(dChemPo > 33) iErrZ=15
	!if(iErrZ=6)'FuPrTcVtot: 0~cmprsblty=',cmprsblty
	!if(iErrTmin.ne.0)iErrZ=5
	!iErrZ=17  if dChemPo > 33
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTcVtot:OMA,OMB',OMA,OMB
	totMol=sum(gMol(1:NC))
	if( vTotCc < 1E-11 .and. LOUD)write(dumpUnit,*)'FuPrTcVtot: 0 ~ vTotCc=',vTotCc
	rhoMol_cc=totMol/vTotCc
	bMix=0
	iErrTmin=0
	TcVolatile=1234
	do iComp = 1,NC
		if(PCj(iComp)==0)then
			iErrZ=11
			if(LOUD)write(dumpUnit,'(a,F7.2,5F8.3)')' FuPrTcVtot: nonsense input. Tc,Pc,L,M,N,c: ',TCj(iComp),PCj(iComp), &
			                                              alphaL(iComp),alphaM(iComp),alphaN(iComp),cVolCc_mol(iComp),tKmin(iComp)
			cycle
		endif
		if(Tcj(iComp) < TcVolatile)then
			TcVolatile=Tcj(iComp)
			iVolatile=iComp
		endif
		aCrit = OMA*Rgas*Rgas*TCj(iComp)*TCj(iComp)/PCj(iComp)
		!sTmp = 0.37464 + 1.54226*ACEN(iComp) - 0.26993*ACEN(iComp)*ACEN(iComp)
		!DLALDT(iComp)= -sTmp*SQRT(Tr/ ALPHA) ! dAlp/dT = 2*sqrt(alp)*S*(-0.5/sqrt(Tr)) 
		!                                        = -sqrt(alp)*S/sqrt(Tr); (T/alp)*dAlp/dT= -S*sqrt(Tr/Alp)    EL2ed Eq.7.18,8.35
		Tr = tKelvin / Tcj(iComp)
		Pow = alphaM(iComp)*alphaN(iComp)
		TrPow = Tr**Pow 
		!ALPHA = TrPow/ Tr**alphaN(iComp) *DEXP( alphaL(iComp)*(1-TrPow) )	 ! = Tr^[N*(M-1)]*exp[ L*(1-Tr^MN) ] Eq. 4
		rLnAlpha=(Pow-alphaN(iComp))*DLOG(Tr) + alphaL(iComp)*(1-TrPow)	 ! -> -L*TrPow as Tr -> inf
		if( (rLnAlpha) < -33)rLnAlpha= -33.d0 !avoid exponential underflow
		ALPHA = EXP(rLnAlpha)  ! -> 0
		!ALPHA = (  1 + ( 0.37464D0 + acen(iComp)*(1.54226D0 - 0.26993D0*ACEN(iComp)) )*( 1 - DSQRT( Tr) )  )**2 ! for debugging
		!ln(ALPHA) = pow*ln(Tr) - N*ln(Tr) + L*(1-Tr^pow)
		! dLn(ALPHA)/dTr = pow/Tr - N/Tr - L*pow*Tr^(pow-1) 
		! d2Ln(alpha)/dTr^2 = -pow/Tr^2 + N/Tr^2 -L*pow*(pow-1)*Tr^(pow-2)
		!(Tr/Alp)*dAlp/dTr = Tr*dLn(alpah)//dTr = Pow - N - L*pow*Tr^Pow   -> -L*pow*Tr^Pow at Tr->inf
		! T^2*d2(ln(alpha)/dT^2 = -pow + N -L*pow*(pow-1)*Tr^pow

		TdLAL_dT(iComp)= Pow - alphaN(iComp) - alphaL(iComp)*pow*TrPow 
		! Ures/RT= -Ares/RT+(a/bRT)*TdLAL_dT*ln(...) = -Ares/RT + Ares/RT*TdLAL_dT
		!              = alphaN_i*[alphaM_i-1-alphaL_i*alphaM_i*Tr^(alphaM*alphaN)]
		! Td/dT[ TdLALdT ] = T*d/dT[ Pow - alphaN(iComp) - pow*alphaL(iComp)*TrPow ]
		T2d2LAL_dT2(iComp) = -pow+alphaN(iComp) -pow*(pow-1)*alphaL(iComp)*TrPow
		ALA(iComp,iComp) = aCrit*ALPHA
		!bVolCc_mol(iComp) = OMB*Rgas*TC(iComp)/Pc(iComp) - cVolCc_mol(iComp)	
		xFrac(iComp)=gMol(iComp)/totMol
		bMix=bMix+xFrac(iComp)*bVolCc_mol(iComp)
		FUGC(iComp)=0 !initial to zero for isZiter=1.
	enddo
	if( tKelvin < tKmin(iVolatile) )iErrTmin=1
	if(iErrTmin.ne.0)then
		if(LOUD.and.initCall)write(dumpUnit,'(a,5f8.2)' )' FuVtot: T,Tmin(i)',tKelvin,( tKmin(i),i=1,NC)
		iErrZ=5
	endif
	if(iErrZ>10)return
	eta=bMix*rhoMol_cc
	if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTcVtot: eta,bMix',eta,bMix
	if(LOUD.and.initCall)write(dumpUnit,'(a,3(1PE11.4))')' FuPrTcVtot:TdLAL_dT,T2d2LAL_dT2',TdLAL_dT(1),T2d2LAL_dT2(1)
	if(ABS(totMol-1) > 1e-5)then
		if(LOUD)write(dumpUnit,*) 'FuPrTcVtot warning: totMol= ',totMol
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
			                  +  ALA(iComp,jComp)*( T2d2LAL_dT2(iComp)+T2d2LAL_dT2( jComp) )/2 ! where d2LALdT2= d/dT[ DLALDT ]
			TdaMixDt = TdaMixDt + TdALA_dT*xFrac(iComp)*xFrac( jComp)
			T2d2aMixDt2 = T2d2aMixDt2 + TdTdALA_dT_dT*xFrac(iComp)*xFrac( jComp)
		enddo
		cMix = cMix + cVolCc_mol(iComp)*xFrac(iComp) ! dCmix_dNi=c(i)
	enddo
	!TdaMixDt = TdaMixDt*aMix
 	!T2d2aMixDt2 = T2d2aMixDt2*aMix
	!BIGMES=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!dU_NKT= BIGMES*(-1+daMixDt/aMix)	!EL2ed Eq.8.35
	if(aMix < zeroTol .and. LOUD)write(dumpUnit,*)'FuPrTcVtot:0~aMix,aCrit,alpha,Tr=',aMix,aCrit,alpha,Tr
	if(bMix < zeroTol .and. LOUD)write(dumpUnit,*)'FuPrTcVtot:bMix~0=',bMix

	eta=bMix*rhoMol_Cc
	crho=cMix*rhoMol_Cc
	!write(dumpUnit,*)'FuPrTcVtot: eta,bMix,cMix=',eta,bMix !,cMix


!	Z=1/(1-brho) - (A/Bsqrt8)*brho/[(1+crho)*(1+brho+2crho)+(brho+crho)*(1-brho)]
!	(1+crho)*(1+brho+2crho)+(brho+crho)*(1-brho)=1+brho+2crho+crho+brho*crho+2crho^2
!                                                                          +brho+crho-brho^2-brho*crho
!                                                                          = 1+2*brho-brho^2 +crho*(4-brho+2crho)
!                                                                         = 1 + eta*bq + eta^2*cq; 
!                                                                         bq = 2+4c/b; cq = -[1-2(c/b)^2]; -qq=bq^2-4cq
!   Aatt =(1/sqrt(-qq))*{ ln[ (2cq*eta+bq-sqrt(-qq))/(2cq*eta+bq+sqrt(-qq) ] - ln[ (bq-sqrt(-qq))/(bq+sqrt(-qq) ] } EL2ed Eq. B.34
!  Alt: RKPR(Cismondi)=>Aatt=(a/bRT)*(1/(d1-d2))*y/[ (1+d1y)(1+d2y) ] ; 
!       1+di*y = 0 => di = -1/yi; yi = [ -bq +/- sqrt(bq^2 - 4cq) ]/2; 
!        denom = [ (1+d1y)(1+d2y) ] = 1 + (d1+d2) y + d1d2 y^2 => bq = (d1+d2) = C+1; -cq = -d1*d2 = 1-2(c/b)^2 = C    ; 
!           d2 = cq/d1 => d1 + cq/d1 = bq => d1^2 + cq -bqd1 = 0 => d1 = ( bq +/- sqrt(bq^2-4*cq) )/2	= (bq +/- sqrtNqq)/2
	c_b = cMix/bMix
	bq=2+4*c_b
	cq= 2*(c_b)*(c_b)-1
	! qq = -(bq^2-4cq) ; sqrtNqq=sqrt(-qq) = (d1-d2)=(1+c/b)*sqrt8
	sqrtNqq=(1+c_b)*sqrt8 ! x = c/b & -qq = (bq*bq-4*cq) = (2+4x)^2-4(2x^2-1)= 4+16x+16x^2-8x^2+4 = 8+16x+8x^2 = 8*(1+x)^2  
	voidFrac=1-eta
	if( voidFrac.le. 0 .and. LOUD)write(dumpUnit,*)'FuPrTcVtot:0~ voidFrac = ',voidFrac
	dArep_dEta=1/ voidFrac
	zRep=eta * dArep_dEta
	d1=1+sqrt2+c_b*(2+sqrt2) !=( bq+sqrtNqq )/2=[(2+4x)+(1+x)sqrt8]/2 = 1+2x+(1+x)*sqrt2 = 1+sqrt2+x*(2+sqrt2) ! = 2cq/(bq-sqrtNqq)
	d2=1-sqrt2+c_b*(2-sqrt2) !=( bq-sqrtNqq )/2=[(2+4x)-(1+x)sqrt8]/2 = 1+2x-(1+x)*sqrt2 = 1-sqrt2+x*(2-sqrt2) ! = 2cq/(bq+sqrtNqq)
	!  2	  +c_b*2          = d1+d2 = 2*(1+x)
	denom = 1+eta*bq + cq*eta*eta ! (1+d1*eta)*(1+d2*eta) = (1+(d2+d1)*eta+d1*d2*eta^2); 
	dAatt_dEta= -aMix/(bMix*Rgas*tKelvin)/denom	 ! Jaubert Eq.7
	zAtt=eta*dAatt_dEta
	zFactor = 1+zRep+zAtt 
	pMPa=zFactor*rhoMol_Cc*Rgas*tKelvin

	!write(dumpUnit,*)'qq,sqrt(-qq)',qq,sqrtNqq
	!write(dumpUnit,*)'b,c/b',bMix,c_b
	!write(dumpUnit,*)'bq,cq',bq,cq
	!write(dumpUnit,*)'zRep,zAtt',zRep,zAtt
	!write(dumpUnit,*)'d1+,d2+',d1Plus,d2Plus
	!write(dumpUnit,*)'d1-,d2-',d1Minus,d2Minus
	BIGA = aMix*pMPa/(Rgas*Rgas*tKelvin*tKelvin)
	BIGB = bMix*pMPa/(Rgas*tKelvin)
	if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTcVtot: zFactor,pMPa=',zFactor,pMPa
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
	!
	!BIGMES=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!argLog1 = (2*cq*eta+bq-sqrtNqq)/(2*cq*eta+bq+sqrtNqq)	! argLog1 may be < 0, but then so is argLog2
	!argLog2 = (         bq-sqrtNqq     )/(         bq+sqrtNqq     )
	argLog1 = 1+d1*eta
	argLog2 = 1+d2*eta

	argLog = argLog1/argLog2
	if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTcVtot: argLog1,2=',argLog1,argLog2
	if(loud .and. argLog < zeroTol)write(dumpUnit,*) 'FuPrTcVtot: argLog 1or2 < 0',argLog1,argLog2
	aResAtt= -BIGA/BIGB/sqrtNqq*DLOG( argLog )
	!write(dumpUnit,*)'a/bRT,alpha', BIGA/BIGB,alpha
	!write(dumpUnit,*)'zAtt,aAtt',zAtt,aResAtt
	aAttRKPR = -BIGA/BIGB/(d1-d2)*DLOG( (1+d1*eta)/(1+d2*eta) )
	zAttRKPR = -BIGA/BIGB*eta/( (1+d1*eta)*(1+d2*eta) )
	denomRKPR= (1+d1*eta)*(1+d2*eta)
	!write(dumpUnit,*)'RKPR: zAtt,aAtt',zAttRKPR, aAttRKPR
	!write(dumpUnit,*)'denom: tcPR,RKPR', denom,denomRKPR
	!write(dumpUnit,*)'check values'
	!BIGMESOld=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!write(dumpUnit,*)'BigMes,BigMesOld:',BigMes,BigMesOld
	aRes= -LOG(1-eta) + aResAtt
	!uDep/RT = beta*d(A/RT)/dBeta = -T*d(A/RT)/dT = -T* [-A/RT^2 + (1/RT)*dA/dT ] = A/RT - (dA/dT)/R = A/RT*[1-(T/alpha)*dAlpha/dT]
	uRes= aResAtt*(1-TdaMixDt/aMix)	! Cv = (dU/dT) = U/T + T*d(U/T)/dT; U/RT = A/RT*(1-TdaMixDt/aMix) => Td(U/RT)/dT = T(dA/RT)/dT -(A/RT)*Td/dT[TdaMixDt/aMix] = -U/RT -A/RT*[ Td(TdaMixDt)/dT -(TdaMixDt/aMix)^2] 
	!uDep_R = aAtt*T*(1-TdLnAlpha/dT) ~ (a/bR)
	!d(uDep_R)/dT = 
	! Td(Ures/RT)/dT = bigmes*[ TdTdaMix_dT_dT/aMix - TdaMixDt^2/aMix^2 ] 
	CvRes_R=  aResAtt*( -T2d2aMixDt2 - TdaMixDt*TdaMixDt/aMix )/aMix
	CvRes_R=  aResAtt*( -TdLAL_dT(1)**2 - T2d2LAL_dT2(1) ) !todo: check generaliztion to mixtures. Above is slightly off(?) 
	!T*dZ/dT = T*dZatt/dT = eta*T*d/dT[dAatt/dEta]	; dAatt/dEta= -(aMix/bRT)/denom
	! T*d/dT[dAatt/dEta]= -T*dLnaMix/dT*(a/bRT)/denom -dAatt_dEta 
	!                                               = (TdaMixDt/aMix)*dAatt_dEta - dAatt_dEta	= dAatt_dEta*(TdaMixDt/aMix - 1)
	!                                                                 NOTE: the -dAatt_dEta comes from the 1/T in a/bRT.

	!zRep=eta/(1-eta) = eta*dArep/dEta
	!d2ARep/dEta2 = 1/(1-eta)^2
	!dZrep/dEta = dArep/dEta+eta*d2Arep/dEta2 
	!eta*dZrep/dEta = zRep + eta^2*d2Arep/dEta2 = eta/(1-eta) + eta^2/(1-eta)^2 = eta/(1-eta)*[ 1+eta/(1-eta) ] = eta/(1-eta)^2
	!zAtt= -eta*aMix/(bMix*Rgas*tKelvin)/( 1+eta*bq + cq*eta*eta ) = eta*dAatt_dEta	 ! Jaubert Eq.7
	!dAatt_dEta= -aMix/(bMix*Rgas*tKelvin)/denom = -(a/bRT)/denom	 ! Jaubert Eq.7
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
		!if(LOUD)write(dumpUnit,*)'(dP/dT)/rhoR=',(zFactor+TdZ_dT)
		! Z=1/(1-eta)-a/bRT*F(eta)=> T*dZ_dT= +a/bRT*F(eta)-F(eta)/bRT*T*da/dT = ZATT*(-1+TdaMixDt/a)
		TdZ_dT = ZATT*(-1+TdaMixDt/aMix)
		CpRes_R = CvRes_R-1 + (zFactor+TdZ_dT)*(zFactor+TdZ_dT)/cmprsblty
	else
		if(LOUD)write(dumpUnit,*)'FuPrTcVtot: 0~cmprsblty=',cmprsblty
		iErrZ=6
		CpRes_R=86.8686D0
	endif

	if(isZiter.ne.0)goto 861  ! It's good to know cmprsblty even if isZiter==1. e.g. slope of spinodal, critPure, ...

	if(zFactor < BIGB)then
		iErrZ=13
		goto 861
	endif
	if(LOUD.and.initCall)write(dumpUnit,*)'aRes,uRes',aRes,uRes
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
			NdC_b_dni= ( cVolCc_mol(iComp) - cMix*bVolCc_mol(iComp)/bMix )/bMix	!NOTE: declared DoublePrecision.
			!BIGB=b*P/RT  ; c_b = cMix/bMix
			ChemPoRes = bVolCc_mol(iComp)/bMix*(zFactor-1) - DLOG(1-eta) & 
			!+ aResAtt*( 2*SUMXA/aMix - bVolCc_mol(ICOMP)/bMix ) + zFactor*(cMix-cVolCc_mol(iComp))*rhoMol_cc ! Privat(2016a) Eq 14
			+ aResAtt*( 2*SUMXA/aMix - bVolCc_mol(ICOMP)/bMix - NdC_b_dni/(1+c_b) ) &
			+ zAtt*(1-eta)*NdC_b_dni/(1+c_b)  
		END IF
		IF(ChemPoRes > 33) THEN
			!	AVOID EXPONENT OVERFLOW. DON'T INDICATE ERROR. IT MAY BE N-HEXANE AT P>1GPa AND 298K.
			iErrZ=7	 ! Warning because exp(33)~ infinity may be a good enough approximation.
			ChemPoRes = 33
		endif
		FUGC(ICOMP) = (ChemPoRes)
	enddo
	if(LOUD.and.initCall)write(dumpUnit,*)'FUGC(1-NC)',(fugc(i),i=1,nc)
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
	END	!Subroutine FuPrTcVtot()

	!SUBROUTINE FugiPrTc(tKelvin,pMPa,gMol,NC,LIQ,FUGC,zFactor,IER)
	Subroutine FugiPrTc( tKelvin,pMPa,gmol,NC,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )	  ! JRE 2020
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
	!        Rgas     GAS CONSTANT ( EG. 8.31434 CC-MPA/(GMOL-K) ) IN PHASE LIQ
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
	!        iErr = 1-9 warning
	!				12 IF LIQUID ROOT WAS FOUND BUT WAS NOT REAL
	!				13 IF VAPOR  ROOT WAS FOUND BUT WAS NOT REAL
	!				14 NEGATIVE LOG CALCULATED
	!				15 THE LOG OF THE FUGACITY COEFFICIENT CALCULATED TO BE LARGE ENOUGH TO INDICATE OVERFLOW.
	!				16 NR DID NOT CONVERGE
	!
	!   NOTE:           UNITS OF ALL THE INPUTS SHOULD BE
	!                   CONSISTENT WITH UNITS OF Rgas.  EXCEPT
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
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision FUGC(NC),xFrac(NC),gMol(NC)
    DoublePrecision tKelvin,zFactor 
	!DIMENSION bVol(nmx),DLALDT(nmx), ALA(nmx,nmx),
	data initCall/1/
	!  prBIPs are passed in from GetPrBIPs()
	IER=0 !sets all array values to zero.
    iErr=0 !error tracker internal for FUGI()

	!zeroTol=1D-11 !zeroTol is now global
	if(tKelvin < zeroTol .and.LOUD)write(dumpUnit,*)'FuPrTc: Nonsense on input. 0~tKelvin'
	totMoles=sum( gMol(1:NC) ) 
	xFrac(1:NC)=gMol(1:NC)/totMoles
	if( ABS(totMoles-1) > zeroTol .and. LOUD)write(dumpUnit,*)'FuPrTc: ??? totMoles,x(1)=',totMoles,xFrac(1) 
	bMix=SUM( xFrac(1:NC)*bVolCc_mol(1:NC) )
	if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTc: bMix=',bMix
	aRes=86.8686	  ! initialize to avoid NAN on error return
	uRes=86.8686
	cvRes_R=86.8686
	cpRes_R=86.8686
	cmprsblty=86.8686

	isZiter=1

	pb_RT = pMPa*bMix/(Rgas*tKelvin)
	eta=pb_RT/1.001d0  	!super high pressures can generate eta>1 because Z>>1. Check eta in if() below.
	if(LIQ==1 .or. LIQ==3 .or. eta>etaMax)eta=etaMax/1.001d0 !organize to improve precision when P~1E-11.
	rhoMol_Cc=eta/bMix
	if(LOUD.and.initCall)write(dumpUnit,'( a,3(1PE11.4) )')' FuPrTc: 1st FuVtot. rhoMol_cc,eta=', rhoMol_cc,eta 
	call FuPrTcVtot(isZiter,tKelvin,1/rhoMol_Cc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
	!write(dumpUnit,*)'check initial values'
	if(iErrZ > 10)then
		iErr=13
		if(LOUD)then
            write(dumpUnit,*)'FuPrTc: Initial FuVtot returned with error code:',iErrZ
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
	pMax= -1234 !for crude goldenz
	pMin=  1234
	isZiter=1
	do while(ABS(change) > 1.D-10 .and. iErr < 10) ! This gives ~10 sig figs on rho for liq and above 
	!                                                                    initialization gives exact ideal gas value when P->1E-11.
		NITER=NITER+1
		if(niter > itMax)iErr=16
		if(niter > itMax .and. LOUD)write(dumpUnit,*)'FugiPrTc: iterations exceeded. eta=',eta
		if(niter > itMax)exit
		!ETA=rhoMol_Cc*bVolMix
		if( (1-eta) < 1.d-11)then
			if(LOUD)write(dumpUnit,*) 'FuPrTc: eta > 1. nIter=',nIter
			if(LOUD)pause
		endif
		rhoMol_cc=eta/bMix
		vTotCc=totMoles/rhoMol_cc
		if(LOUD.and.initCall)write(dumpUnit,*)'FuPrTc:',NITER,'th FuVtot. eta=', eta 
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
		if(iErrZ>10)then
			iErr=13
			if(LOUD)write(dumpUnit,*)'FuPrTc: FuVtot returned with error code:',iErrZ
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
		if(LOUD.and.initCall)write(dumpUnit,'(a,2F8.4,2E11.4)')' eta,Z,error,change',eta,zFactor,error,change
		!write(dumpUnit,*)'check error'
		etaOld=eta
		errOld=error
        if(zFactor < 0)change=1 !force another iteration if Z < 0. This happens when P ~ 1E-12.
		if(ABS(change/eta) > 0.1)change=eta*DSIGN(0.1d0,change)	!step limiting. 0.3>0.2 to give max looseness
		!if(ABS(change/zFactor) > 1)change=zFactor*DSIGN(0.5d0,change)	!step limiting. Stop zFactor from going negative.
		if(LOUD.and.initCall)write(dumpUnit,*)'Fugi:eta,change',eta,change
		eta=eta-change
		rhoMol_cc=eta/bMix
		if(eta < 0 .or. eta > etaMax)then !NOTE: this should not happen given step limiting
			if(niter < (itMax-1))niter=itMax-1 !restrict tries with new guess so we can crash if necessary.
			if(LOUD)write(dumpUnit,'(a,f8.4)')' Warning in FuTpt: next guess is eta=',eta
			eta=0
			if(liq==1)eta=etaMax
			if(LOUD)write(dumpUnit,*) 'Restarting iteration with eta=',eta
			if(LOUD)pause
		endif
	enddo  !iteration on eta to find P=P(input)
	!write(dumpUnit,*)'Hello! iteration concluded. iErr,eta=',iErr,eta
	etaPass=eta											! Accurate rho is essential for PsatEar at low T (e.g. propane).	
	if(iErr > 10)then ! use crude goldenz to keep iterations going. 
		!write(dumpUnit,*)'Hello! iErr.ne.0 . LIQ,etaAtPmax,etaAtPmin=',LIQ,etaAtPmax,etaAtPmin
		if(liq==0 .or. LIQ==2)eta=etaAtPmax  !crude goldenz
		if(liq==1 .or. LIQ==3)eta=etaAtPmin  !crude goldenz
		if(eta < zeroTol)eta=zeroTol
		!write(dumpUnit,*)'Hello! iErr.ne.0. eta,bMix=',eta,bMix
		rhoMol_cc=eta/bMix
		vTotCc=totMoles/rhoMol_cc
		if(LOUD)write(dumpUnit,'(a,i3,1x,f7.4,a,f10.5)')' FuPrTc: iErr,goldenEta=',iErr,eta
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,xFrac,NC,FUGC,zFactor,aDep,uDep,iErrZ)
		if(LOUD)write(dumpUnit,'(a,f10.5)')' FuPrTc:  using goldenZ=',zFactor
		iErr=10
	endif
	IF (eta < 0)then
		iErr=12
		goto 861
		!call BeepMsg(errMsg(iErr))
	endif
	!  ITERATION ON RHO HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.
	rhoMol_cc=eta/bMix
	vTotCc=totMoles/rhoMol_cc
    if(LIQ > 1)goto 861
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS and derivative props 
	!                       (cmprsblty,CvRes_R,CpRes_R passed by GlobConst) 
	!
	isZiter=0
	call FuPrTcVtot(isZiter,tKelvin,1/rhoMol_Cc,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrZ)
	IF (zFactor < 0 )then
		iErr=14
		goto 861
	else
		FUGC(1:NC)=FUGC(1:NC)-DLOG(zFactor)	! convert from T,V to T,P. 
	endif
	if(LOUD)write(dumpUnit,*)'aRes,uRes',aRes,uRes
	if(iErrZ>10)then
		iErr=13
		goto 861
	elseif(iErrZ.ne.0)then
		iErr=iErrZ
	endif  	
861	continue
	if(iErr > 0 .and.LOUD)write(dumpUnit,*)'FugiPrTc: FuVtot last call. iErr,iErrZ=',iErr,iErrZ
	initCall=0
	RETURN
	END

	SUBROUTINE ZITERDum(zFactor,BIGB,A1,A0,iErr)	 !ZITER() is defined in FugiPR. 
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!
	!  PURPOSE    - CALCULATE zFactor FROM NEWTON-RAPHSON ITERATION
	!
	!      IER    - 200 DIDNT CONVERGE IN 25 ITERATIONS
	!
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit doublePrecision(A-H,O-Z)
	iErr = 0
	do kount=1,25
		F = -A0 + zFactor*(  A1+zFactor*( -(1-BIGB)+zFactor )  )
		DF = 3*zFactor*zFactor - (1-BIGB)*2*zFactor + A1
		ZN = zFactor - F/DF
		ERR = DABS((ZN - zFactor)/ZN)
		zFactor = ZN
		IF(ERR.lt. 1.D-9) GO TO 86
	enddo
	iErr = 200 !only to reach here is if do loop has exceeded iterations
86 continue
	RETURN
	END

