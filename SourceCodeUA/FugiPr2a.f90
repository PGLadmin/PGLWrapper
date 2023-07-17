
MODULE PREosParms
	USE GlobConst, only:nmx
	Parameter( THIRD=1.D0/3,sqrt2=1.414213562373095d0,sqrt8=2*sqrt2 )
	!etac = 1/( 1+(4-sqrt8)**third+(4+sqrt8)**third )
	Parameter( etac=0.253076586541599d0,OMA = (40*etac+8)/(49-37*etac),OMB = etac/(etac+3)  )
	DoublePrecision svKappa1(nmx)! For the PRsvWS model. Stryjek and Vera, can j chem eng, 64:323 (1986)
END MODULE PREosParms
subroutine GetPR(NC,iErrGet)
	!  
	!  PURPOSE:  LOOKS UP THE PR PARAMETERS AND STORES THEM IN USEd PREosParms
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	!    USEd PREosParms, GlobConst, BIPs 
	!  Programmed by:  JRE 07/00
	USE GlobConst
	USE PREosParms
	Implicit NONE !doublePrecision(A-H,K,O-Z)
	character*88 bipFile
	integer GetBIPs,iErrGet,NC,nComps
	iErrGet=0
	etaMax=1-zeroTol
	bVolCc_mol(1:NC)=OMB*Rgas*Tc(1:NC)/Pc(1:NC)
	! note:  bips are passed back through USEd BIPs/

	nComps=NC
	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)

	bipFile=TRIM(PGLinputDir)//'\BipPengRob.txt' ! // is the concatenation operator
	iErrGet=GetBIPs(bipFile,ID,NC)
	tKmin(1:NC)=0.4d0*Tc(1:NC)

	RETURN
END	!subroutine GetPR(NC,iErrGet)

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
	!        Rgas     GAS CONSTANT ( EG. 8.31434 CC-MPA/(GMOL-K) ) IN PHASE LIQ
	!        tKelvin  ABSOLUTE TEMPERATURE
	!        PMPa     ABSOLUTE PRESSURE
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
	!         GUNDERSEN, T., COMPUTERS IN CHEM. ENG., 6:245 (1982).
	!
	!         SOAVE, G., CHEM. ENG. SCI., 27:1197 (1972).
	!
	!         PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!****************************************************************
	SUBROUTINE FuPrVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,aRes,uRes,iErr)
	USE GlobConst
	USE PREosParms
	USE BIPs
	Implicit DoublePrecision(a-h,k,o-z)
	DoublePrecision FUGC(NC),gmol(NC),eta,Tr !,zFactor,qCoeff,rCoeff
	DoublePrecision xFrac(nmx),ALA(nmx,nmx),DLALDT(nmx)
	!DoublePrecision, STATIC:: THIRD,sqrt2,sqrt8,etac,OMA,OMB ! STATIC keeps these variables in memory so they can be reused with recomputing them at every call.
	!  prBIPs are passed in from GetPrBIPs()

	iErr = 0 ! vector init
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	totMol=SUM( gmol(1:NC) )
	rhoMol_cc=totMol/vTotCc
	do iComp = 1,NC
		aCrit = OMA*Rgas*Rgas*Tc(iComp)*Tc(iComp)/Pc(iComp)
		sTmp = 0.37464D0 + 1.54226D0*ACEN(iComp) - 0.26993D0*ACEN(iComp)*ACEN(iComp)
		Tr=tKelvin/Tc(iComp)
		ALPHA = (  1 + sTmp*( 1 - DSQRT( Tr) )  )*(  1 + sTmp*( 1 - DSQRT( Tr) )  )
		DLALDT(iComp)= -sTmp*SQRT( Tr / ALPHA)	   !EL2ed Eq. 7.18
		ALA(iComp,iComp) = aCrit*ALPHA
		xFrac(iComp)=gmol(iComp)/totMol
	enddo
	if(ABS(totMol-1) > 1e-5)then
		if(LOUD)write(dumpUnit,*)'FuPrVtot only works for xFrac not mole numbers. Sorry.'
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
		bMix = bMix + bVolCc_mol(iComp)*xFrac(iComp)
	enddo

	eta=bMix*rhoMol_Cc


	!	qCoeff = BIGA - BIGB - BIGB*BIGB	!for SRK
	!	rCoeff = BIGA*BIGB
	zFactor=( 1/(1-eta)-aMix/(bMix*Rgas*tKelvin)*eta/(1+2*eta-eta*eta) )!EL2ed 7.15
	PMPa=zFactor*rhoMol_Cc*Rgas*tKelvin
	if(isZiter==1)return
	!
	!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
	!
	BIGA = aMix*PMPa/(Rgas*Rgas*tKelvin*tKelvin)
	BIGB = bMix*PMPa/(Rgas*tKelvin)
	qCoeff = BIGA - 2*BIGB - 3*BIGB*BIGB
	rCoeff = BIGB*( BIGA-BIGB*(1+BIGB) )

	if(zFactor < BIGB)then
		iErr=14
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
			iErr=14
		ELSE
			dChemPo = bVolCc_mol(iComp)/bMix*(zFactor-1) - DLOG(zFactor-BIGB) + FATT*(2*SUMXA/aMix - bVolCc_mol(ICOMP)/bMix)
		END IF
		IF(dChemPo > 33) THEN
			!	AVOID EXPONENT OVERFLOW BUT INDICATE ERROR.
			iErr=15
			dChemPo = 0
		endif
		FUGC(ICOMP) = (dChemPo)
	enddo
	aRes= -DLOG(1-eta)+FATT
	uRes=FATT*(1-daMixDt/aMix)

861	continue
	RETURN
	END

	Subroutine FugiPR( tKelvin,pMPa,xFrac,NC,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )
	!SUBROUTINE FugiPR(tKelvin,pMPa,xFrac,NC,LIQ,FUGC,zFactor,IER)
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
	!        tKelvin     ABSOLUTE TEMPERATURE
	!        pMPa     ABSOLUTE PRESSURE
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
	USE GlobConst
	USE PREosParms
	USE BIPs
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision FUGC(NC),xFrac(NC) !
	DoublePrecision ALA(nmx,nmx),bVol(nmx),DLALDT(nmx)
    DoublePrecision zFactor,rCoeff,qCoeff
	LOGICAL bEven
	!  prBIPs are passed in from GetPrBIPs()
	!  Note: In most cases, it's better to call FuVtot to iterate on Z, but for PREOS it really doesn't matter. 
	iErr = 0 ! 
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	do iComp = 1,NC
		aCrit = OMA*Rgas*Rgas*Tc(iComp)*Tc(iComp)/Pc(iComp)
		sTmp = 0.37464 + 1.54226*ACEN(iComp) - 0.26993*ACEN(iComp)*ACEN(iComp)
		Tr = tKelvin/Tc(iComp)
		ALPHA = (  1 + sTmp*( 1-DSQRT( Tr) )  )**2
		DLALDT(iComp)= -sTmp*SQRT( Tr/ ALPHA) ! dAlp/dT = 2*sqrt(alp)*S*(-0.5/sqrt(Tr)) = -sqrt(alp)*S/sqrt(Tr); (T/alp)*dAlp/dT= -S*sqrt(Tr/Alp)	 EL2ed Eq.7.18,8.35
		ALA(iComp,iComp) = aCrit*ALPHA
		bVol(iComp) = OMB*Rgas*Tc(iComp)/Pc(iComp)
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

	BIGA = aMix*pMPa/(Rgas*Rgas*tKelvin*tKelvin)
	BIGB = bMix*pMPa/(Rgas*tKelvin)

	!	qCoeff = BIGA - BIGB - BIGB*BIGB	!for SRK
	!	rCoeff = BIGA*BIGB
	qCoeff = BIGA - 2*BIGB - 3*BIGB*BIGB
	rCoeff = BIGB*(BIGA-BIGB*(1+BIGB))
	IF( bEven(LIQ) )THEN !LIQ=2 for spinodal vapor calc
		zFactor = 3
		CALL ZITER(zFactor,BIGB,qCoeff,rCoeff,IERZ)
		IF(IERZ > 9)iErr=16
	ELSE  !LIQ=3 for spinodal LIQ calc
		zFactor=0
		CALL ZITER(zFactor,BIGB,qCoeff,rCoeff,IERZ)
		IF(IERZ > 9)iErr=16
	END IF

	if( bEven(LIQ) )then
	  etaV=BIGB/zFactor
	  zFactorV=zFactor
	ELSE
	  etaL=BIGB/zFactor
	  zFactorL=zFactor
	ENDIF
!
!  CALCULATE FUGACITY COEFFICIENTS OF INDIVIDUAL COMPONENTS
!
	if(zFactor < BIGB)then
		iErr=14
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
	aRes= dA_NKT
	uRes= dU_NKT
	rhoMol_cc=pMPa/(zFactor*Rgas*tKelvin)

	DO iComp = 1,NC
		SUMXA = 0
		DO jComp=1,NC
			SUMXA = SUMXA + xFrac(jComp)*ALA(iComp,jComp)
		enddo
		
		IF( (zFactor-BIGB) < zeroTol .OR. (zFactor+(1+sqrt2)*BIGB) < zeroTol) THEN
			!	AVOID CALCULATION OF NEGATIVE LOGARITHMS BUT INDICATE ERROR.
			dChemPo = 0
			iErr = 14
		ELSE
			dChemPo = bVol(iComp)/bMix*(zFactor-1) - DLOG(zFactor-BIGB)- BIGMES*(2*SUMXA/aMix - bVol(ICOMP)/bMix)
		END IF
		IF(dChemPo > 33) THEN
			!	AVOID EXPONENT OVERFLOW BUT INDICATE ERROR.
			iErr = 15
			dChemPo = 0
		endif
		FUGC(ICOMP) = (dChemPo)
	enddo

861	continue
	RETURN
	END

	SUBROUTINE ZITER(zFactor,BIGB,A1,A0,iErr)
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!
	!  PURPOSE    - CALCULATE zFactor FROM NEWTON-RAPHSON ITERATION
	!
	!      IER    - 200 DIDNT CONVERGE IN 25 ITERATIONS
	!
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit NONE !DoublePrecision(A-H,O-Z)
	DoublePrecision zFactor,BIGB,A1,A0, F,DF,ZN,ERR
	Integer iErr,kount
	iErr = 0
	do kount=1,25
		F = -A0 + zFactor*(  A1+zFactor*( -(1-BIGB)+zFactor )  )
		DF = 3*zFactor*zFactor - (1-BIGB)*2*zFactor + A1
		ZN = zFactor - F/DF
		ERR = DABS((ZN - zFactor)/ZN)
		zFactor = ZN
		IF(ERR < 1.D-9) GO TO 86
	enddo
	iErr = 200 !only to reach here is if do loop has exceeded iterations
86	continue
	RETURN
	END

