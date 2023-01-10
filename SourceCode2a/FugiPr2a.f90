
	subroutine GetPR(NC,iErrGet)
	!  
	!  PURPOSE:  LOOKS UP THE PR PARAMETERS AND STORES THEM IN COMMON
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	!    commons: ppData, BIPs through GetBips
	!  Programmed by:  JRE 07/00
	USE GlobConst
	implicit doublePrecision(A-H,K,O-Z)
	character bipFile*88
	integer GetBIPs
	iErrGet=0
	OMB = 0.07779607d0						  
	if(LOUD)write(*,*)'ID  NAME       TCK   PCMPa      w     '
	do i=1,nc
	  if(LOUD)write(*,'(i4,1x,a11,f6.1,f6.3,1x,f6.3)')ID(i),NAME(i),TC(i),PC(i),ACEN(i)
	  bVolCC_mol(i)=OMB*8.31434*TC(i)/PC(i)
	enddo

	! note:  bips are passed back through common/BIPs/
		bipFile=TRIM(PGLinputDir)//'\BipPengRob.txt' ! // is the concatenation operator
	IERRCODE=GetBIPs(bipFile,ID,NC)

	RETURN
	END

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!$ FuPrVtot
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
	!		  J.R. Elliott, C.T. Lira, Introductory Chemical Engineering Thermodynamics, p 334 (1999). 2ed p262 (2012).
	!		
	!         PENG, D-Y, ROBINSON, D.B. IEC Fun,15: 59–64
	!
	!         GUNDERSEN, T., COMPUTERS IN CHEM. ENG., 6:245 (1982).
	!
	!         SOAVE, G., CHEM. ENG. SCI., 27:1197 (1972).
	!
	!         PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!****************************************************************
	SUBROUTINE FuPrVtot(isZiter,tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
	USE GlobConst
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	DIMENSION FUGC(NC),xFrac(NC),IER(12)
	DIMENSION ALA(NMX,NMX),bVol(NMX),DLALDT(NMX)
	COMMON/DEPFUN/dU_NKT,dA_NKT,dS_NK,dH_NKT
	COMMON/eta/etaL,etaV,zFactorL,zFactorV
	!  prBIPs are passed in from GetPrBIPs()
	DATA initKALL/0/
		!OMA = 0.45723553d0
		!OMB = 0.07779607d0
	IF(initKALL.EQ.0)then
		initKALL=1
		THIRD=1.D0/3
		IDDum=ID(nc)
		sqrt2=DSQRT(2.D0)
		sqrt8=2*sqrt2
		etac = 1+(4-sqrt8)**third+(4+sqrt8)**third	 !Jaubert, Eq.7
		etac = 1/etac
		OMA = (40*etac+8)/(49-37*etac)
		OMB = etac/(etac+3)						  
	endif

	ier = 0 ! vector init
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	totMol=0
	do iComp = 1,NC
		aCrit = OMA*RGAS*RGAS*TC(iComp)*TC(iComp)/PC(iComp)
		sTmp = 0.37464D0 + 1.54226D0*ACEN(iComp) - 0.26993D0*ACEN(iComp)*ACEN(iComp)
		Tr=tAbs/TC(iComp)
		ALPHA = (  1 + sTmp*( 1 - DSQRT( Tr) )  )**2
		DLALDT(iComp)= -sTmp*SQRT( Tr / ALPHA)	   !EL2ed Eq. 7.18
		ALA(iComp,iComp) = aCrit*ALPHA
		bVol(iComp) = OMB*RGAS*TC(iComp)/PC(iComp)
		totMol=totMol+xFrac(iComp)
	enddo
	if(ABS(totMol-1) > 1e-5)then
		if(LOUD)pause 'FuPrVtot only works for xFrac not mole numbers. Sorry.'
		goto 861
	endif

	aMix = 0
	daMixDt=0
	bMix = 0
	DO iComp = 1,NC
		DO jComp = 1,NC
			ALA(iComp,jComp) = SQRT(ALA(iComp,iComp)*ALA(jComp,jComp))*(1-KIJ(iComp,jComp))
			aMix = aMix + ALA(iComp,jComp)*xFrac(iComp)*xFrac(jComp)
			daMixDt = daMixDt + ALA(iComp,jComp)*xFrac(iComp)*xFrac(jComp)*(DLALDT(iComp)+DLALDT(jComp))/2
		enddo
		bMix = bMix + bVol(iComp)*xFrac(iComp)
	enddo

	eta=bMix*rhoMol_Cc


	!	qCoeff = BIGA - BIGB - BIGB*BIGB	!for SRK
	!	rCoeff = BIGA*BIGB
	zFactor=( 1/(1-eta)-aMix/(bMix*RGAS*tAbs)*eta/(1+2*eta-eta*eta) )!EL2ed 7.15
	pAbs=zFactor*rhoMol_Cc*RGAS*tAbs
	if(isZiter)return
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
	!
	BIGA = aMix*pAbs/(RGAS*RGAS*tAbs*tAbs)
	BIGB = bMix*pAbs/(RGAS*tAbs)
	qCoeff = BIGA - 2*BIGB - 3*BIGB*BIGB
	rCoeff = BIGB*( BIGA-BIGB*(1+BIGB) )

	if(zFactor < BIGB)then
		ier(4)=1
		goto 861
	endif
	FATT= -BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	DO iComp = 1,NC
		SUMXA = 0
		DO jComp=1,NC
			SUMXA = SUMXA + xFrac(jComp)*ALA(iComp,jComp)
		enddo
		
		IF( (zFactor-BIGB) < 0 .OR. (zFactor+(1+sqrt2)*BIGB) < 0) THEN
			!	AVOID CALCULATION OF NEGATIVE LOGARITHMS BUT INDICATE ERROR.
			dChemPo = 0
			IER(4) = 1
		ELSE
			dChemPo = bVol(iComp)/bMix*(zFactor-1) - DLOG(zFactor-BIGB) + FATT*(2*SUMXA/aMix - bVol(ICOMP)/bMix)
		END IF
		IF(dChemPo .GT.  33) THEN
			!	AVOID EXPONENT OVERFLOW BUT INDICATE ERROR.
			IER(5) = 1
			dChemPo = 0
		endif
		FUGC(ICOMP) = (dChemPo)
	enddo
	IF(IER(4).NE.0.OR.IER(5).NE.0.OR.IER(6).NE.0)IER(1)=1
	dU_NKT= FATT*(1-daMixDt/aMix)
	dH_NKT=dU_NKT+zFactor-1
	dA_NKT= -LOG(1-BIGB/zFactor)+FATT
	dS_NK=dU_NKT-dA_NKT !+LOG(zFactor)
	aDep=dA_NKT	!to pass as calling argument as well as common
	uDep=dU_NKT
	sRes_R = dS_NK

	if(LIQ==0 .or. LIQ==2)then
	  etaV=BIGB/zFactor
	  zFactorV=zFactor
	ELSE
	  etaL=BIGB/zFactor
	  zFactorL=zFactor
	ENDIF

861	continue
	RETURN
	END

	SUBROUTINE FugiPR(tAbs,pAbs,xFrac,NC,LIQ,FUGC,zFactor,IER)
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
	USE GlobConst
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	DIMENSION FUGC(NC),xFrac(NC),IER(12)
	DIMENSION ALA(NMX,NMX),bVol(NMX),DLALDT(NMX)
	COMMON/DEPFUN/dU_NKT,dA_NKT,dS_NK,dH_NKT
	COMMON/eta/etaL,etaV,zFactorL,zFactorV
	!  prBIPs are passed in from GetPrBIPs()
	DATA initKALL/0/
	IF(initKALL.EQ.0)then
		initKALL=1
		THIRD=1.D0/3
		IDDum=ID(nc)
		sqrt2=DSQRT(2.D0)
		sqrt8=2.D0*sqrt2
		etac = 1+(4-sqrt8)**third+(4+sqrt8)**third	 !Jaubert, Eq.7
		etac = 1/etac
		OMA = (40*etac+8)/(49-37*etac)
		OMB = etac/(etac+3)						  
	endif
	IER = 0 ! vector init
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	do iComp = 1,NC
		aCrit = OMA*RGAS*RGAS*TC(iComp)*TC(iComp)/PC(iComp)
		sTmp = 0.37464 + 1.54226*ACEN(iComp) - 0.26993*ACEN(iComp)*ACEN(iComp)
		Tr = tAbs/TC(iComp)
		ALPHA = (  1 + sTmp*( 1-DSQRT( Tr) )  )**2
		DLALDT(iComp)= -sTmp*SQRT( Tr/ ALPHA) ! dAlp/dT = 2*sqrt(alp)*S*(-0.5/sqrt(Tr)) = -sqrt(alp)*S/sqrt(Tr); (T/alp)*dAlp/dT= -S*sqrt(Tr/Alp)	 EL2ed Eq.7.18,8.35
		ALA(iComp,iComp) = aCrit*ALPHA
		bVol(iComp) = OMB*RGAS*TC(iComp)/PC(iComp)
	enddo

	aMix = 0
	daMixDt=0
	bMix = 0
	DO iComp = 1,NC
		DO jComp = 1,NC
			ALA(iComp,jComp) = SQRT(ALA(iComp,iComp)*ALA(jComp,jComp))*(1-KIJ(iComp,jComp))
			aMix = aMix + ALA(iComp,jComp)*xFrac(iComp)*xFrac(jComp)
			daMixDt = daMixDt + ALA(iComp,jComp)*xFrac(iComp)*xFrac(jComp)*(DLALDT(iComp)+DLALDT( jComp))/2
		enddo
		bMix = bMix + bVol(iComp)*xFrac(iComp)
	enddo
 	!FATT= -BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	!dU_NKT= FATT*(1-daMixDt/aMix)	!EL2ed Eq.8.35

	BIGA = aMix*pAbs/(RGAS*RGAS*tAbs*tAbs)
	BIGB = bMix*pAbs/(RGAS*tAbs)

	!	qCoeff = BIGA - BIGB - BIGB*BIGB	!for SRK
	!	rCoeff = BIGA*BIGB
	qCoeff = BIGA - 2*BIGB - 3*BIGB*BIGB
	rCoeff = BIGB*(BIGA-BIGB*(1+BIGB))
	IF(MOD(LIQ,2) == 0)THEN !LIQ=2 for spinodal vapor calc
		zFactor = 3
		CALL ZITER(zFactor,BIGB,qCoeff,rCoeff,IERF)
		IF(IERF > 9)IER(6)=1
	ELSE  !LIQ=3 for spinodal LIQ calc
		zFactor=0
		CALL ZITER(zFactor,BIGB,qCoeff,rCoeff,IERF)
		IF(IERF > 9)IER(6)=1
	END IF

	if(MOD(LIQ,2)==0)then
	  etaV=BIGB/zFactor
	  zFactorV=zFactor
	ELSE
	  etaL=BIGB/zFactor
	  zFactorL=zFactor
	ENDIF
	!if(LIQ.gt.1)return !for spinodal, don't need fugc's


!
!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
!
	if(zFactor < BIGB)then
		ier(4)=1
		goto 861
	endif
	if(zFactor > 0)etaPass=BIGB/zFactor	! PsatEar needs eta directly.

	BIGMES=BIGA/BIGB/sqrt8*DLOG( (zFactor+(1+sqrt2)*BIGB)/(zFactor+(1-sqrt2)*BIGB) )
	dU_NKT= BIGMES*(-1+daMixDt/aMix)	!EL2ed Eq. 8.35, 8.37
	dH_NKT=dU_NKT+zFactor-1
	dA_NKT= -LOG(1-BIGB/zFactor)-BIGMES
	aRes_RT= dA_NKT
	uRes_RT= dU_NKT
	hRes_RT= dH_NKT
	dS_NK=dU_NKT-dA_NKT !+LOG(zFactor)

	DO iComp = 1,NC
		SUMXA = 0
		DO jComp=1,NC
			SUMXA = SUMXA + xFrac(jComp)*ALA(iComp,jComp)
		enddo
		
		IF( (zFactor-BIGB).LT. 0 .OR. (zFactor+(1+sqrt2)*BIGB).LT.0) THEN
			!	AVOID CALCULATION OF NEGATIVE LOGARITHMS BUT INDICATE ERROR.
			dChemPo = 0
			IER(4) = 1
		ELSE
			dChemPo = bVol(iComp)/bMix*(zFactor-1) - DLOG(zFactor-BIGB)- BIGMES*(2*SUMXA/aMix - bVol(ICOMP)/bMix)
		END IF
		IF(dChemPo .GT.  33) THEN
			!	AVOID EXPONENT OVERFLOW BUT INDICATE ERROR.
			IER(5) = 1
			dChemPo = 0
		endif
		FUGC(ICOMP) = (dChemPo)
	enddo
	IF(IER(4).NE.0.OR.IER(5).NE.0.OR.IER(6).NE.0)IER(1)=1

861	continue
	RETURN
	END

	SUBROUTINE ZITER(zFactor,BIGB,A1,A0,IER)
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

