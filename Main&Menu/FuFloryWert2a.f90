	!  
	!  PROGRAMMED BY:  JRE 2/02
	!  REVISION DATE:  
	! 
	!  LOOKS UP THE TPT PROPERTIES AND RETURNS them in commons
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
    SUBROUTINE GetFloryWert(nComps,idComp,iErrCode)
    USE GlobConst
	USE Assoc !FloryWert calls Wertheim()
	USE BIPs  !BIPs for physical interactions.
    IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555)
	character bipFile*88
	character inFile*251
	!DIMENSION XC(nComps)
    DIMENSION idComp(*)
	dimension IDA(ndb),bondVolA(ndb),DHA(ndb),DHD(ndb),NDSA(ndb),NASA(ndb),NDegreeA(ndb)
	!common/FloryWert/solParm(nmx),vLiq(nmx),eHbKcalMol(NMX),bondVolNm3FH(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	!solParmD(ndb),vLiqD(ndb),eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX) from ParmsFloryWert.txt

!	COMMON/BIPs/KIJ,KTIJ,aBipAD,aBipDA,xsTau,xsTauT,xsAlpha	
!		COMMON/BIPs/KIJ(NMX,NMX),KTIJ(NMX,NMX),HIJ(NMX,NMX),HTIJ(NMX,NMX),xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)


	iErrCode=0

	nC=nComps

	! NOTE: ParmsCrit is read elsewhere, including solParm and vLiq.
	IF(DEBUG)then 
	  OPEN(50,FILE='c:\SPEAD\CalcEos\input\ParmsFloryWert.TXT')
    ELSE 
		inFile=TRIM(masterDir)//'\input\ParmsFloryWert.txt' ! // is the concatenation operator
		OPEN(50,FILE=inFile)
	ENDIF

	READ(50,*,ERR=862)NdeckFloryWert 
	DO I=1,NdeckFloryWert
		!nDegree(i,1)=1
		READ(50,*,ERR=862)IDA(I),bondVolA(I),nDegreeA(I),NASA(I),NDSA(I),DHA(I),DHD(i)
		nTypes(i)=1	! only one type per molecule in this model.
		idType(i,1)=IDA(i)
	enddo   
	CLOSE(50)

	DO J=1,nComps
		iGotIt=0
		NDegree(J,1)=0
		NDonors(J,1)=0
		nAcceptors(J,1)=0
		DO I=1,NdeckFloryWert
			IF(IDA(I).EQ.idComp(J))THEN
				iGotIt=1
				bVolCC_mol(i)=vLiq(i)*0.4d0 ! let 0.4 be the standard liquid eta.
				nTypes(j)=1	! only one type per molecule in this model.
				idType(j,1)=IDA(i)
				bondVolNm3(J,1)=bondVolA(I)
				eAcceptorKcal_Mol(J,1)=DHA(I)
				eDonorKcal_Mol(J,1)=DHD(I)
				NDegree(J,1)=NDegreeA(I)
				nAcceptors(J,1)=NASA(I)
				nDonors(J,1)=NDSA(I)
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			write(*,*)'GetFloryWert Warning: Parms not found in ParmsFloryWert.txt for compId=',idComp(j)
			if(LOUD)pause     'Check that hbonding is not applicable or define in ".txt file'
		endif
	enddo

	write(*,*)'ID    Nd   bondVolNm3    eHbKcalMol'
	do i=1,nComps
	  write(*,601)idComp(i),NDegree(i,1),bondVolNm3(i,1),eAcceptorKcal_Mol(i,1),eDonorKcal_Mol(i,1)
	enddo
601	format(i6,1x,i4,6E12.4)
606	format(i6,1x,a11,f6.1,f6.3,1x,f7.2,f8.2,f6.1,f8.1,i4,f8.6,f8.4)


	! note:  bips are passed back through common/BIPs/
	IF(DEBUG)THEN
		bipFile='c:\SPEAD\CalcEos\input\BipFloryWert.txt'
	ELSE
		bipFile=TRIM(masterDir)//'\input\BipFloryWert.txt' ! // is the concatenation operator
	ENDIF
	!IERRCODE=GetBIPs(bipFile,idComp,nComps)

	RETURN
862	continue
	write(*,*)'GetFlory error - error reading ParmsFloryWert.txt. Path?'
	write(*,*)'nDeck,iCompo',NdeckFloryWert,j
	if(LOUD)pause
	return                      
	END

	SUBROUTINE FuFloryWert(tKelvin,pMpa,xFrac,nComps,LIQ,FUGC,zFactor,iErr)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!$ FUGI
	!
	!   REVISION DATE:  OCTOBER 31, 1995
	!					April 6, 2001, convert to f90
	!   ORIGINAL PROGRAM BY:  J.R. ELLIOTT, JR. (JAN. 1983)
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
	!             Franzen, Elliott, and Kyu, Macro 28:5147(95).
	!		
	!Example1.  oPvp+pvb at 406K,0.1MPa	Franzen, Elliott, and Kyu, Macro 28:5147(95).
	! id	  xi	  mw	 solp	   vL	eHb(kcal)	KcStar	Nd	NAS	NDS	XA	XD	lnPhiAssoc	lnPhiPhys	gamma   
	!23041	0.3255	39953	10.35	37500	7.948	1.64E-6	377	1	0	0.87	1.00  -40.517	41.615	2.998	oPvPyridineFEK40
	!24089	0.6745	87991	 9.14	64900	7.948	3.31E-6	473	0	1	1.00	0.95   -5.587	 5.600	1.013	pvButyralFEK  88	!
	!dXA1/dN1=0.087; dXA1/dN2=-0.042;  dXD2/dN1=-0.115; dXD2/dN2=0.0558;  aAssoc=-16.96
	!FYI: LCST at ~404.3, xLo=0.87, xUp=0.92
	!Example2.  pPvp+pEVOH at 333K,0.1MPa	Keskin, Elliott, IECR 42:6334(03).
	!id		x	  mw	solp	vL	KcStar	eHb(kcal)	Nd	NAS	NDS	XA  	XD  	const	sumdxn	lnPhiAssoc	gamma   
	!26068	0.5	67879	8.6	66637	2.76E-05	4.00	474	1	0	0.9598	1.0000	-9.918	17.552	7.634	pBuMeAcry068
	!27064	0.5	63310	12	60296	2.33E-05	4.00	878	1	1	0.9630	0.9413	-44.247	-17.552	-61.799	pEt56VOh44_063
	!dXAdN	 0.039885405	 0.036887261	dXDdN	0.00000	 0.01502
	!		-0.039885402	-0.036887259			0.00000	-0.01502
	!FYI: LCST ~296K.  Amazingly sensitive to temperature.  Remarkable asymmetry in gamma's.  Gibbs-Duhem OK?
	!****************************************************************
	USE GlobConst !, only: avoNum,zeroTol,bVolCC_mol,etaMax,etaPure,half,ID,isTPT,isESD,iEosOpt,LOUD,DEBUG,Rgas,MasterDir,Tc,Pc,acen
	USE Assoc !XA,XD,XC passed.
	USE VpDb ! for VpCoeffs
	USE BIPs
	implicit doublePrecision(A-H,O-Z)

	character bipFile*88,bipHbFile*88
!	character*123 ErrMsg(11)

	!doublePrecision KIJ,KTIJ
	DIMENSION FUGC(nComps),xFrac(nComps) !,dFUGASSOC_dT(NMX),dFUGASSOC_dRHO(NMX)
	DIMENSION volFrac(NMX),fugAssoc(NMX) !,vMolecNm3(NMX) !,KIJ(NMX,NMX),KTIJ(NMX,NMX)
	!COMMON/BIPs/KIJ,KTIJ,HIJ,HTIJ,xsTau,xsTauT,xsAlpha
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!COMMON/eta/etaL,etaV,zFactorL,zFactorV
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	!  prBIPs are passed in from GetPrBIPs()

!	COMMON/BIPs/KIJ(NMX,NMX),KTIJ(NMX,NMX),HIJ(NMX,NMX),HTIJ(NMX,NMX),xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)

	!eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX) from ParmsFloryWert.txt
	!common/FloryWert/solParm(nmx),vLiq(nmx),eHbKcalMol(NMX),bondVolNm3FH(NMX),ND(NMX),NDS(NMX),NAS(NMX)

	!for wertheim1c: common/Assoc/eHbKcal_mol(nmx,maxTypes),bondVolNm3(nmx,maxTypes),nDegree(nmx,maxTypes),nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),idType(nmx,maxTypes),localType(maxTypesGlobal),idLocalType(maxTypes),nTypes(NMX),nTypesTot
	Bip(iComp,jComp)=kij(iComp,jComp)+kTij(iComp,jComp)/tKelvin
	
	!  COMPUTE THE MOLECULAR PARAMETERS A AND B AND THEIR CROSS COEFFS
	vLiqMix=0
	do iComp = 1,nComps
		vLiqMix=vLiqMix+xFrac(iComp)*vLiq(iComp)
	enddo

	!  note:  bips are passed back through common/BIPs/
	IF(DEBUG)then 
		bipFile='c:\spead\CalcEos\input\BipSpead.txt'
	ELSE 
		bipFile=TRIM(masterDir)//'\input\BipSpead.txt' ! // is the concatenation operator
    ENDIF
	!iErrCode=GetBIPs(bipFile,ID,nC)
	!Note: no need to check for switching because kij and ktij are symmetric
    if(LOUD)then
	    if(iErrCode.ne.0)pause 'GetTpt Error: GetBIPs returned error.'
    end if

	!load aBipAD matrix
	IF(DEBUG)then 
		bipHbFile='c:\Spead\CalcEos\Input\BipAD.txt'
	ELSE 
		bipHbFile=TRIM(masterDir)//'\input\BipAD.txt' ! // is the concatenation operator
	ENDIF
	call GetAssocBips(bipHbFile,aBipAD,ierABip) !in WertheimVv.f90. idLocalType,nTypesTot in common/assoc 

	!load aBipDA matrix
	IF(DEBUG)then 
		BipHbFile='c:\Spead\CalcEos\Input\BipDA.txt'
	ELSE 
		BipHbFile=TRIM(masterDir)//'\input\BipDA.txt' ! // is the concatenation operator
	ENDIF
	call GetAssocBips(bipHbFile,aBipDA,ierABip) !in WertheimVv.f90

	if(iEosOpt.eq.8)CALL GetVp(nComps,ID,iErrVp)	
 
  	solParmAvg=0
	bipAvg=0
	DO iComp = 1,nComps
		volFrac(iComp)=xFrac(iComp)*vLiq(iComp)/vLiqMix
		solParmAvg = solParmAvg + solParm(iComp)*volFrac(iComp)
		do jComp=1,nComps
			!bip=kij(iComp,jComp)+kTij(iComp,jComp)/tKelvin
			vFracJ=xFrac(jComp)*vLiq(jComp)/vLiqMix 
			bipAvg=bipAvg+volFrac(iComp)*solparm(iComp)*vFracJ*solParm(jComp)*Bip(iComp,jComp)
		enddo
	enddo
			
	iErr=0

	zFactor=1
	fugc(iComp)=1
	etaV=0
	if(LIQ.eq.1)then
		etaL=0.4
		zFactor=pMpa*vLiqMix/(rGas*tKelvin)
		do kComp=1,nComps
			delDelta=(solParm(kComp)-solParmAvg) ! cf Sandler 2ed p.341
			bipPart=0
			do iComp=1,nComps
				bipTmp=kij(iComp,kComp)+kTij(iComp,kComp)/tKelvin
				bipPart=bipPart+volFrac(iComp)*solParm(iComp)*bipTmp
			enddo
			actCoeff=vLiq(kComp)*( delDelta*delDelta+2*solParm(kComp)*bipPart-bipAvg )/(RgasCal*tKelvin)
			fugc(kComp)= actCoeff !actCoeff is really LnGam  
			!vMolecNm3(kComp)=vLiq(kComp)*etaL/602.22
		enddo
		rhoMol_cc=1/vLiqMix
		isZiter=0 ! V is assumed in FloryWert, so no need to iterate on Z.
		!vMolecNm3(1:nComps)=bVolCC_mol(1:nComps)/avoNum
		call Wertheim(isZiter,etaL,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,fugAssoc,iErrCode) !nDs,nAs,ND,
		call WertheimFugc(xFrac,nComps,etaL,fugAssoc,h_nMichelsen,iErrF)	!NDS,NAS,ND
		if(iErrF > 0)iErr=11

		TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
		PcEos(1:nComps)=Pc(1:nComps)
		ZcEos(1:nComps)=Zc(1:nComps)

		write(*,'(a,1x,2(a),a,a,a)')' Compo','  xFrac','   vFrac  ','     LnGamFH','     LnPhiAssoc','     gamma  '
		do kComp=1,nComps
			pSat=pc(kComp)*10**( 7*(1+acen(kComp))/3*(1-tc(kComp)/tKelvin) )
			if(pSat.lt.1e-33)pSat=1E-33 !logArg in fugc() cannot be zero
			rLnGamma=fugc(kComp)+fugAssoc(kComp)
			if(rLnGamma.gt.709)rLnGamma=709 !exp overflow for larger
			gamma=exp(rLnGamma)
			write(*,'(i3,1x,2(f10.7),f13.7,f13.7,e11.4)')kComp,xFrac(kComp),volFrac(kComp),fugc(kComp),fugAssoc(kComp),gamma
			fugc(kComp)=LOG(pSat/pMPa) + fugc(kComp) + fugAssoc(kComp)	!log of product is sum of logs 
		enddo
	endif

	RETURN
	END

