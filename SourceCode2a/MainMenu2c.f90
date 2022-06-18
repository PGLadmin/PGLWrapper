	!MAIN PROGRAM FOR CALLING EOS SUBROUTINES. Includes GetCrit, GetBips, BipIo.
	!	SETS UP THE CRITS, EOSPARMS, bips, DEBUG status, FUGI(iEosOpt)
	!	Echoes user IO to Output.txt, and reports error checks. 
	!	INITIALIZATION AND CALLING SEQUENCE FOR VLE, LLE, VLLE SUBROUTINES.
	!reqd routines:
	!	bubpl.for, bubtl.for, dewtv.for, flashsub.for, FuEsdMy.for FuEsdXs2.for, FugiPr.f90, FugiPrws.for, RegPure.f90
	!	LmDifEzCov2.for, Mintools.for
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE PORTLIB
	USE GlobConst
	USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	CHARACTER ANSWER*1
	CHARACTER calcType*2
	CHARACTER errMsgPas*77
	CHARACTER outFile*251

	DIMENSION BIPTMP(NMX)
	DIMENSION ZFEED(NMX),S(NMX)
    INTEGER ier(11),idCas(NMX)
	
	COMMON/CrOpt/calcType
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/FEED/ZFEED
	COMMON/KVALUES/S
	!common/FloryWert/vLiq(nmx)

	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	iChange =VERIFY('c:\msdev\Projects\calceos',CURDIR)
	iChange2=VERIFY('C:\MSDEV\Projects\CalcEos',CURDIR)
	if(iChange2==0)iChange=0
	if(iChange==0 .or. iChange > 26)DEBUG=.TRUE.
	DEBUG=.FALSE.   !setting to true mostly directs input from c:\spead...	 Logic above automatically updates to TRUE if CURDIR=c:\msdev\Projects\calceos
	!DEBUG=.TRUE.    ! Logic not working JRE 20200407
	LOUD = .FALSE.
	!LOUD = .TRUE.
	masterDir=TRIM(curDir)
	outFile=TRIM(masterDir)//'\output\output.txt' ! // is the concatenation operator
	if(LOUD)print*,'MasterDir=',TRIM(masterDir)
	OPEN(UNIT=52,FILE=outFile)
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! AFG 2011: Option 9 is added
	write(*,*)'Enter EOS option: 1=PR,2=SymEsd,3=PRWS,4=Esd,5=SPEAD,6=FloryWert,7=NRTL,'
	write(*,*)'    8=SpeadGamma, 9=SPEAD11, 10=PcSaft, 11=tcPR, 12=GCESD, 13=TransSPEAD'
	!NewEos: Add here for listing new iEosOpt.
	read(*,*)iEosOpt
	write(52,*)'Enter EOS option: 1=PR,2=SymEsd,3=PRWS,4=Esd,5=SPEAD,6=FloryWert,7=NRTL,8=SpeadGamma, 9=SPEAD11, 10=PcSaft, 11=tcPR, 12=GCESD, 13=TransSPEAD'
	write(52,*)iEosOpt
	if(iEosOpt==9) then
		Write(*,*)' Please note that SPEAD11 is a correlation for n-alkanes > C7'
		Write(*,*)' So, using it for other components will result in nonsense results. '
		Write(*,*)' Also, its applicability to mixtures has not been tested, yet.'
		Write(*,*)' AFG'
	endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!NewEos: Add here for listing new iEosOpt.
	NC=-1
	do while(NC.LT.0 .or. nc.gt.55)
		WRITE(6,*)'ENTER NUMBER OF COMPONENTS (0 TO QUIT)   '
		READ(5,*)NC
		WRITE(52,*)' NUMBER OF COMPONENTS (0 TO QUIT)   '
		write(52,*)NC
		IF(NC.EQ.0)goto 86
		IF(NC.LT.0.or.nc.gt.55)then
			write(*,*)'PLEASE TYPE -1 < NC < 55'
		ENDIF
	enddo
	INITIAL=0

	WRITE(6,*)'ENTER ID NOS. OF THE COMPONENTS OR ZEROS FOR LIST  '
	READ(5,*)(ID(I),I=1,NC)
	WRITE(52,*)' ID NOS. OF THE COMPONENTS OR ZEROS FOR LIST  '
	write(52,*)(ID(I),I=1,NC)
	call IdCasLookup(NC,idCas,iErrCas,errMsgPas)
	if(LOUD.and.iErrCas)pause 'Main: idCAS lookup failure. '
	CALL GETCRIT(NC,iErrCrit)
	if(iErrCrit.ne.0)then
		write(*,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		if(LOUD)pause
		goto 86
	endif
	iErrGet=0 
	if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	if(iEosOpt.eq.2)CALL GetEsdCas(NC,idCas,iErrGet)	 !Results placed in USE EsdParms
	if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet) 
	if(iEosOpt.eq.4)CALL GetEsdCas(NC,idCas,iErrGet)
	if(iEosOpt.eq.5)CALL GetTpt(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	if(iEosOpt.eq.8)CALL GetTpt(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.9)CALL GetTpt(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	!if(iEosOpt.eq.9)CALL GetRG(NC,ID,iErrGet)		!AFG 2011 : Reads White's RG method parameters, disabled 20191121 by JRE. unused.
	if(iEosOpt.eq.10)CALL GetPcSaft(NC,ID,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters
	if(iEosOpt.eq.11)CALL GetPrTc(NC,iErrGet) 	!JRE 2019 : Reads Gross's PcSaft parameters
	!NewEos: Add here for initializing parms.
	if(iErrGet.gt.0)then
		if(LOUD)write(*,*)'Error in Main: failed to get required dbase props.'
		if(LOUD)pause
		goto 86
	endif
	if(LOUD)WRITE(*,*)'THE KIJ MATRIX AT 298K IS ESTIMATED AS'
	if(LOUD)WRITE(*,'(8X,9I8)')(ID(I),I=1,NC)
	DO I=2,NC
		DO J=1,I-1
			BIPTMP(J)=KIJ(I,J)+KTIJ(I,J)/298
		ENDDO
		if(LOUD)WRITE(6,601)I,ID(I),(BIPTMP(J),J=1,I-1)
	enddo
	if(LOUD)WRITE(*,*)'THE HIJ MATRIX AT 298K IS ESTIMATED AS'
	if(LOUD)WRITE(*,'(8X,9I8)')(ID(I),I=1,NC)
	DO I=2,NC
		DO J=1,I-1
			BIPTMP(J)=HIJ(I,J)+HTIJ(I,J)/298
		ENDDO
		if(LOUD)WRITE(6,601)I,ID(I),(BIPTMP(J),J=1,I-1)
	enddo
	ANSWER='N'
	WRITE(*,*)'DO YOU WANT TO MODIFY THE BIP MATRIX? (Y/N)'
	if(NC.gt.1)READ(*,'(A1)')ANSWER
	WRITE(52,*)'DO YOU WANT TO MODIFY THE BIP MATRIX? (Y/N)'
	write(52,'(A2)')ANSWER
	IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')call BIPIO(NC) !note: results passed through common/bips

	iQuit=0
	do while(iQuit.ne.1) 
		WRITE(6,*)'ENTER TYPE OF PHASE EQUILIBRIUM CALCULATION'
		WRITE(6,*)'BD FOR 1-comp Bifurcation (Multiroot) Diagram   '
		WRITE(6,*)'BI FOR LLE Binodal   '
		WRITE(6,*)'BF FOR BUBBLE PRESSURE AT FILE CONDITIONS'
		WRITE(6,*)'BP FOR BUBBLE POINT PRESSURE'
		WRITE(6,*)'BT FOR BUBBLE POINT TEMPERATURE'
		WRITE(6,*)'DT FOR DEW POINT TEMPERATURE'
		WRITE(6,*)'FU FOR FUGACITY CALCULATION Given T, P'
		WRITE(6,*) 'FV FOR FUGACITY CALCULATION GIVEN T,V'
		WRITE(6,*)'CR FOR MIXTURE CRITICAL POINT CALCULATION-can be used for Options 4 & 5'		  !Added by AGF, Dec. 2009
		WRITE(6,*)'CP FOR PURE COMPONENT CRITICAL POINT CALCULATION-can be used for Options 4 & 5'	  !Added by AGF, Dec. 2009
		WRITE(6,*)'GE FOR Excess Gibbs and Enthalpy CALCULATION '
		WRITE(6,*)'GX FOR Excess Gibbs at T,eta '
		WRITE(6,*)'IC FOR isochore plot on single component system'
		WRITE(6,*)'IT FOR isotherm plot on single component system'	  !Modified by AGF, Dec. 2009
		WRITE(6,*)'KD FOR Kij regression on multi-system Database'
		WRITE(6,*)'KI FOR KIJ ITERATION ON BINARY AT SINGLE POINT'
		WRITE(6,*)'KN FOR NEW BIP INPUT'
		WRITE(6,*)'KO FOR SINGLE KIJ ITERATION FOR FILE OF EXPTL DATA'
		WRITE(6,*)'LL FOR LLE   '
		WRITE(6,*)'OB FOR OPTIMUM BIPS WITH MULTIPLE BIPS'
		WRITE(6,*)'PE FOR PHASE ENVELOPE (P,T given X,Y)'
		WRITE(6,*)'PS FOR POLYMER-SOLVENT PARTITIONING ESTIMATION '
		WRITE(6,*)'PT FOR Pure component Property Table '
		WRITE(6,*)'PX FOR P,X,Y given T'
		WRITE(6,*)'RP FOR REGRESSION OF PURE COMPONENT ESD PARAMETERS'
		WRITE(6,*)'RS FOR Kii regression on temperature-vapor pressure Database for SPEAD EOS'
		WRITE(6,*)'RX FOR REGRESSION OF XA,XD (e.g. from FTIR data)'
		WRITE(6,*)'SC FOR SUPERCRITICAL FLUID SOLID SOLUBILITY'
		WRITE(6,*)'TP FOR 3-PHASE CALCULATION  '
		WRITE(6,*)'TX FOR T,X,Y given P'
		WRITE(6,*)'VL FOR VLE FLASH   '
		WRITE(6,*)'VP FOR vapor pressure   '
		WRITE(6,*)'FC FOR FREEZING CURVE CALC 4 CRYSTAL solid' 
		WRITE(6,*)'FO FREEZCURVE KIJ optimization'	
		WRITE(6,*)'KS FOR Tpt State Interpolation'
		WRITE(6,*)'EA FOR Equal Area rule calculation (vapor pressure)'
		WRITE(6,*)'FW FOR FOR FITTING THE WHITE RG METHOD PARAMETER'
		WRITE(6,*)'VC FOR VIRIAL COEFFICIENTS FOR PURE COMPONENTS, EOS OPTIONS: 4 & 5'  !AFG 2011
		WRITE(6,*)'QT TO QUIT'


		WRITE(6,*)'Note: previous output files may be deleted after your selection. '
		READ(5,'(A2)')calcType
		write(52,*)' The requested type of calculation is:'
		write(52,'(A2)')calcType

		!Added by AGF
		!AUG 10
		IF(calcType.EQ.'CR'.OR.calcType.EQ.'cr') THEN
			isZiter=3		   !Will calculate first derivatives of desired variables in respect to N
		ELSEIF(calcType.EQ.'BP'.OR.calcType.EQ.'bp'.OR.calcType.EQ.'CP'.OR.calcType.EQ.'cp') THEN
			isZiter=4		   !Will calculate first derivatives of desired variables in respect to T,RHO and P
		ELSE
			isZiter=0
		ENDIF 
		!---------------------------------

		IF(calcType.ne.'QT'.and.calcType.ne.'qt')then
			outFile=TRIM(masterDir)//'\output\TPXY.txt'
			open(61,file=outFile)
			write(61,*)' ' !open and clear old contents
			close(61)
		endif

		IF(calcType.EQ.'BD'.OR.calcType.EQ.'bd')CALL BifurcationDiagram(NC)
		IF(calcType.EQ.'BI'.OR.calcType.EQ.'bi')CALL BINODAL(NC)
		IF(calcType.EQ.'BF'.OR.calcType.EQ.'bf')CALL BpFromFile(NC)
		IF(calcType.EQ.'BP'.OR.calcType.EQ.'bp')CALL BPITER(NC)	   !In the future we can add isZiter to other options if it is necessary, AUG 10
		IF(calcType.EQ.'BT'.OR.calcType.EQ.'bt')CALL BTITER(NC)
		IF(calcType.EQ.'DT'.OR.calcType.EQ.'dt')CALL DTITER(NC)
		IF(calcType.EQ.'FC'.OR.calcType.EQ.'fc')CALL FREEZINGCURVECALC(NC)
		IF(calcType.EQ.'FO'.OR.calcType.EQ.'fo')CALL FreezOptKij(NC)
		IF(calcType.EQ.'FU'.OR.calcType.EQ.'fu')CALL FITER(NC)
		IF(calcType.EQ.'CR'.OR.calcType.EQ.'cr')CALL CRITER(NC)	                !Added by AGF, Dec. 2009	AUG 10
		IF(calcType.EQ.'CP'.OR.calcType.EQ.'cp')CALL CPITER(NC) 	            	!Added by AGF, Dec. 2009	AUG 10
   		IF(calcType.EQ.'GE'.OR.calcType.EQ.'ge')CALL GibbsXS(NC)
		IF(calcType.EQ.'GX'.OR.calcType.EQ.'gx')CALL GibbsXsEta(NC)
		IF(calcType.EQ.'IC'.OR.calcType.EQ.'ic')CALL Isochore(NC,iErrCode,errMsgPas)
		IF(calcType.EQ.'IT'.OR.calcType.EQ.'it')CALL Isotherm(NC)				!Modified by AGF, Dec. 2009
		IF(calcType.EQ.'KD'.OR.calcType.EQ.'kd')CALL KIJDB(NC)
		IF(calcType.EQ.'KI'.OR.calcType.EQ.'ki')CALL KIJITR(NC)
		IF(calcType.EQ.'KN'.OR.calcType.EQ.'kn')call BipIo(NC)
		IF(calcType.EQ.'KO'.OR.calcType.EQ.'ko')CALL KIJOPT(NC)
		IF(calcType.EQ.'LL'.OR.calcType.EQ.'ll')CALL LLITER(NC)
		IF(calcType.EQ.'OB'.OR.calcType.EQ.'ob')CALL OptiBIPs(NC)
		IF(calcType.EQ.'PE'.OR.calcType.EQ.'pe')CALL PhasEnv(NC)
		IF(calcType.EQ.'PS'.OR.calcType.EQ.'ps')CALL PSITER(NC,iErrCode)
		IF(calcType.EQ.'PT'.OR.calcType.EQ.'pt')CALL PropTable(NC,ier)
		IF(calcType.EQ.'PX'.OR.calcType.EQ.'px')CALL PXYT(NC)
		IF(calcType.EQ.'RP'.OR.calcType.EQ.'rp')CALL RegPureIo(NC)
		IF(calcType == 'RS'.OR.calcType == 'rs')CALL RegSpeadIo(NC)
		IF(calcType.EQ.'RX'.OR.calcType.EQ.'rx')CALL RegXaXdIo(NC)
		IF(calcType.EQ.'SC'.OR.calcType.EQ.'sc')CALL SCFITER(NC)
		IF(calcType.EQ.'TP'.OR.calcType.EQ.'tp')CALL TEITER(NC)
		IF(calcType.EQ.'TX'.OR.calcType.EQ.'tx')CALL TXYP(NC)
		IF(calcType.EQ.'VL'.OR.calcType.EQ.'vl')CALL VLITER(NC)
		IF(calcType.EQ.'VP'.OR.calcType.EQ.'vp')CALL VpIter(NC) !,gmol,tKelvin,Psat,vLCc,vVCc)
		IF(calcType.EQ.'FW'.OR.calcType.EQ.'fw')CALL FitWhite(NC)				  	!AFG 2011
		IF(calcType.EQ.'VC'.OR.calcType.EQ.'vc')CALL VirialCoeff(NC)				!AFG 2011
		IF(calcType=='KS'.OR.calcType=='ks')CALL RegA0A1A2TptIo(NC)
		!IF(calcType.EQ.'IT'.OR.calcType.EQ.'it')CALL Isotherm(NC,ier,errMsnPas)
		IF(calcType.EQ.'QT'.OR.calcType.EQ.'qt')iQuit=1
	enddo  !while(iQuit.ne.1)

601	FORMAT(I3,I5,9(1X,F7.4))
86	continue
	close(52)
	stop
	END

	!***********************************************************************

	Subroutine BipIo(NC) !note: results passed through common/bips
	USE GlobConst
	USE BIPs
	implicit DoublePrecision(a-h,k,o-z)
	character answer*1
	DO J=2,NC
		DO I=1,J-1
			WRITE(*,*)'ENTER KIJ,KTIJ FOR THE PAIR',I,J,'    '
			READ(*,*)KIJ(I,J),KTIJ(I,J)
			WRITE(52,*)' KIJ,KTIJ FOR THE PAIR',I,J,'    '
			write(52,*)KIJ(I,J),KTIJ(I,J)
			KIJ(J,I)=KIJ(I,J)
			KTIJ(J,I)=KTIJ(I,J)
			if(iEosOpt.EQ.3 .or. iEosOpt.EQ.7)THEN
				WRITE(*,602)i,j,j,i
				READ(*,*)xsTau(I,J),xsTau(J,I)
				WRITE(52,602)i,j,j,i
				write(52,*)xsTau(I,J),xsTau(J,I)

				xsAlpha(I,J)=0.3d0
				if(ABS(xsTau(I,J)).ge.0.1.or.ABS(xsTau(J,I)).ge.0.1)then
					write(*,*)'Accept xsNrtlAlpha=0.3? (y/n) '
					read(*,'(a)')answer
					write(52,*)'Accept xsNrtlAlpha=0.3? (y/n) '
					write(52,'(a)')answer
					if(answer.ne.'y'.and.answer.ne.'Y')then
						write(*,*)'Enter xsNrtlAlpha'
						read(*,*)xsAlpha(I,J)
						write(52,*)'Enter xsNrtlAlpha'
						write(52,*)xsAlpha(I,J)
					endif
				endif
				xsAlpha(J,I)=xsAlpha(I,J)
			endif !iEosOpt.eq.3
		enddo
		WRITE(6,601)J,ID(J),(KIJ(J,I),I=1,J-1)
	enddo !J=2,NC

	if(iEosOpt.ge.4.and. ItsEven(iEosOpt))then
		DO J=2,NC
			DO I=1,J-1
				WRITE(*,*)'ENTER HIJ,HTIJ FOR THE PAIR',I,J,'    '
				READ(*,*)HIJ(I,J),HTIJ(I,J)
				WRITE(52,*)'ENTER HIJ,HTIJ FOR THE PAIR',I,J,'    '
				write(52,*)HIJ(I,J),HTIJ(I,J)
				HIJ(J,I)=HIJ(I,J)
				HTIJ(J,I)=HTIJ(I,J)
			enddo
			WRITE(6,601)J,ID(J),(HIJ(J,I),I=1,J-1)
		enddo
	endif !iEosOpt.ge.4

601	FORMAT(I3,I5,9(1X,F7.4))
602	format(' ENTER xsTau(',i2,',',i2,'),xsTau(',i2,',',i2,') OR ZEROES TO NULLIFY XS EFFECT')
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer Function ItsEven(iArg)
    integer iArg
	ItsEven=1-MOD(iArg,2)
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ErrCheck(ier,iflag)
	DIMENSION ier(*)
      IFLAG=0
	DO I=1,11
		IF (ier(I).NE.0)IFLAG=1
	enddo
	IF (IFLAG.EQ.1)WRITE(6,*)'ERROR !!- CHECK ANSWERS CAREFULLY'
	IF(ier(2).EQ.1)WRITE(*,*)'ERROR ALL UPPER PHASE'
	IF(ier(3).EQ.1)WRITE(*,*)'ERROR ALL LOWER PHASE'
	IF(ier(4).EQ.4)WRITE(*,*)'ERROR IN FUGI-NEG LOG CALCD'
	IF(ier(4).EQ.5)WRITE(*,*)'ERROR IN FUGI-FUGACITY OVERFLOWS'
	IF(ier(4).EQ.6)WRITE(*,*)'ERROR IN FUGI-Z ITERATION NO CNVRG'
	IF(ier(5).EQ.1)WRITE(*,*)'ERROR VAPOR AND LIQUID ROOTS CLOSE'
	IF(ier(6).EQ.1)WRITE(*,*)'ERROR VLE ITERATION NO CNVRG'
	IF(ier(7).EQ.1)WRITE(*,*)'ERROR VLE ITERATION FAILED TO IMPROVE'
	IF(ier(8).EQ.1)WRITE(*,*)'ERROR IN TC,PC,OR X,Y'
	IF(ier(9).EQ.1)WRITE(*,*)'ERROR P SPECIFIED < 0'
	IF(ier(10).EQ.1)WRITE(*,*)'ERROR T SPECIFIED IS UNREASONABLE'
	IF(ier(11).EQ.1)WRITE(*,*)'ERROR MORE THAN 10 COMPONENTS OR ITMAX<1'
	IF(ier(12).EQ.1)WRITE(*,*)'FUGACITIES NOT EQUAL ON LAST ITER'
	if(iflag.ne.0)write(*,*)'Error on last iteration for next pt.'
	if(iflag.ne.0)ier(1)=1
	RETURN
	END
	    
	subroutine BeepMsg(msg)
	USE GlobConst
	character msg
	if(LOUD) then
	    write(*,*)TRIM(msg)
        pause
    end if
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdDipprLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst
    integer ier
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas
	character inFile*251
	dimension idCas(NC)
	dimension idCasDb(maxDb),idDb(maxDb)
	inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT'
	if(DEBUG)inFile='c:\spead\idTrcDipCas.TXT'
	OPEN(50,FILE=inFile)
	!OPEN(50,FILE='c:\spead\idTrcDipCas.TXT')
!	IF(DEBUG)THEN
!		OPEN(50,FILE='c:\spead\idTrcDipCas.TXT')
!	ELSE
!		inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT' ! // is the concatenation operator
!		OPEN(50,FILE=inFile)
!	ENDIF
	errMsg(0)='No Problem'
	errMsg(1)='IdCasLookup Error: at least one id not found'
	ier=0
	iGotIt=0
	read(50,*)nDeck
	!print*,'IdLookup: nDeck=',nDeck
	do i=1,nDeck
		read(50,*,ioStat=ioErr)idcc,idDb(i),idCasDb(i)
		if(ioErr)print*,'i,idDippr=',i,idDb(i)
		if( idCas(1)==idCasDb(i) )then
			iGotIt=1
			id(1)=idDb(i)
		endif
	enddo
	if(iGotIt==0)then
		ier=1
		if(LOUD)pause 'IdDipprLookup: idDippr not found for comp1'
		goto 861
	endif
	if(NC==1)goto 861
	do iComp=2,NC
		iGotIt=0
		do i=1,nDeck
			if(idCasDb(i).eq.idCas(iComp))then
				id(iComp)=id(i)
				iGotit=1
				cycle
			endif
		enddo
		ier=1 !only way here is if cycle failed => iGotIt==0.
		write(*,*)'Did not find idCas(i).i,id=',iComp,id(i)
		goto 861
	enddo
861	errMsgPas=errMsg(ier)
	close(50)
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!**********************************************

	!C      AARON'S SUBROUTINE
	!C      $GETCRIT.FOR   VERSION 1.0
	SUBROUTINE GETCRIT(NC,iErrCode)
	!C  
	!C  PROGRAMMED BY:  A.S. PUHALA - AUG 93
	!C  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	!C  REVISION DATE:  2/02 - FOR binary data file
	!C 
	!C  LOOKS UP THE CRITICAL PROPERTIES AND RETURNS JUST TC,PC,ACEN,NAME
	!C
	!C  INPUT
	!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!C  OUTPUT
	!C    TC - CRITICAL TEMPERATURE
	!C    PC - CRITICAL PRESSURE
	!C    ACEN - ACENTRIC FACTOR
	!C    NAME - COMPONENT NAME
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	CHARACTER NAMED*30,tCode*4,pCode*4,vCode*4,form*12
	character readText*132,dumText*77 
	PARAMETER(ndb=1555)
	character infile*251
	DIMENSION IDNUM(ndb),TCD(ndb),PCD(ndb),ACEND(ndb),NAMED(ndb)
	DIMENSION ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb)
	common/FloryWert/solParm(nmx),vLiq(nmx),vMolec(NMX),&
		eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	iErrCode=0
	inFile=TRIM(masterDir)//'\input\ParmsCrit.dta' ! // is the concatenation operator
	IF(DEBUG)inFile='c:\SPEAD\CalcEos\input\ParmsCrit.dta' 
	if(LOUD)write(*,*)'CritFile=',TRIM(inFile)
	OPEN(40,FILE=inFile,FORM='BINARY')
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
	READ(40,ERR=861)NDECK1
    if(LOUD)then
	    if(ndeck1.gt.ndb)pause 'GetCrit: more data in file than allocated'
    end if
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
		READ (40,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
!		write(61,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
!			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
!	close(61)
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1).EQ.0)THEN
		DO I=1,NDECK1,3
			WRITE(*,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
            if(LOUD)then
			    IF(((I-1)/45)*45.EQ.(I-1))PAUSE
            end if
		enddo
		RETURN
	END IF
	CLOSE(40)
    if(LOUD)then
	    if((nDeck1).gt.ndb)pause 'GetCrit: too much data in file'
    end if

	IF(DEBUG)THEN
		OPEN(40,FILE='c:\spead\calceos\input\ParmsCrAdd.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	ENDIF
	READ(40,*,ERR=862)dumText

	nDeck2=0
	DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ(40,'(A132)',ERR=862,END=200)readText
		nDeck2=nDeck2+1
		I=nDeck1+nDeck2
		read(readText,*,ERR=862)IDNUM(I),TCD(I),PCD(I),ZCD(I) &
			,ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
		cycle
200		EXIT !terminate loop
	enddo
	close(40)
	nDeck=nDeck1+nDeck2
    if(nDeck > ndb .and. LOUD)pause 'GetCrit: too much data in ParmsCrAdd file'
	if(LOUD)Write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	DO iComp=1,NC
		iGotIt=0
		DO J=1,NDECK
			IF(IDNUM(J).EQ.ID(iComp)) THEN
				NAME(iComp)=NAMED(J)
				TC(iComp)=TCD(J)
				PC(iComp)=PCD(J)
				ID(iComp)=IDNUM(J)
				ACEN(iComp)=ACEND(J)
				ZC(iComp)=ZCD(J)
				rMwPlus(iComp)=rMwD(j)
				solParm(iComp)=solParmD(j)
				vLiq(iComp)=vLiqD(j)
				iGotIt=1
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			if(LOUD)write(*,*)'Error in GetCrit: not found for iComp=',iComp
			if(LOUD)pause
		endif
		if(LOUD)write(*,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp),TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	if(LOUD)write(*,*)'GetCrit error - error reading ParmsCrit.dta. Path? Debug?'
	if(LOUD)write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	if(LOUD)pause
	return                      
862	continue
	if(LOUD)write(*,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	if(LOUD)write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	if(LOUD)pause
	return                      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetVp(NC,IdDum,iErrCode,VpCoeffs)
	!	PROGRAMMED BY: AV 06/22/06
	!THIS SUBROUTINE GIVES VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	!INPUT:
	!NC - NUMBER OF COMPONENTS
	!	   ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!OUTPUT:
	!VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	USE GlobConst  !Passes ID, LOUD, NMX,masterDir, and DEBUG values, so don't redeclare them here.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1975,maxPts=333)
	character infile*251 !masterDir*234,
	DIMENSION IDNUM(NDB),rMINTD(NDB) ,VALMIND(NDB) ,rMAXTD(NDB),VALMAXD(NDB),AVGDEVD(NDB),NUMCOEFFD(NDB) ,vpCoeffsd(NDB,5)
	DIMENSION IdDum(NMX),rMINT(NMX) ,VALMIN(NMX) ,rMAXT(NMX),VALMAX(NMX),AVGDEV(NMX),NUMCOEFF(NMX) ,vpCoeffs(NMX,5)  
	common PDAT(maxPts),XDAT(maxPts),TDAT(maxPts),ntps
	iErrCode=0
	IF(DEBUG)then 
		OPEN(40,FILE='c:\SPEAD\CalcEos\input\CoeffsVp2a.txt') !dta',FORM='BINARY')
    ELSE 
		inFile=TRIM(masterDir)//'\input\CoeffsVp2a.txt' !dta' ! // is the concatenation operator
		OPEN(40,FILE=inFile,ioStat=ioErr) !,FORM='BINARY')
		if(ioErr.and.LOUD)pause 'GetVp: a copy of CoeffsVp2a must be in the ..\CalcEos\input dir'
	ENDIF
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	if(IDdum(1).ne.ID(1)) then
		if (LOUD)pause 'GetVp: ID passed .ne. ID global. Is that intended???'
	end if
	I=0
	ndeck1=1974
	READ(40,*,ERR=861)NDECK1
    if(LOUD)then
	    if(ndeck1.gt.ndb)pause 'GetVp: more data in file than allocated'
    end if
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
		READ(40,*,ERR=861)IDNUM(I),rMINTD(I) ,VALMIND(I) ,rMAXTD(I),VALMAXD(I),AVGDEVD(I),NUMCOEFFD(I) ,(vpCoeffsd(i,icoeff),icoeff=1,5)  
	enddo
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1).EQ.0)THEN
		DO I=1,NDECK1,3
			!WRITE(*,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
            if(LOUD)then
			    IF(((I-1)/45)*45.EQ.(I-1))PAUSE
            end if
		enddo
		RETURN
	END IF
	CLOSE(40)

	!IF(DEBUG)THEN
		!OPEN(41,FILE='c:\spead\calceos\input\CoeffsVp.TXT')
	!ELSE
		!inFile=TRIM(masterDir)//'\input\CoeffsVp.TXT' ! // is the concatenation operator
	!	OPEN(41,FILE=inFile)
	!ENDIF
	!READ(41,*,ERR=862)dumText

	!if((nDeck).gt.ndb)pause 'GetCrit: too much data in file'
	!nDeck2=0
	!DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		!READ(40,'(A132)',ERR=862,END=200)readText
		!nDeck2=nDeck2+1
		!I=nDeck1+nDeck2
		!read(readText,*,ERR=862)IDNUM(I),rMINTD(I) ,VALMIND(I) ,rMAXTD(I),VALMAXD(I),AVGDEVD(I),NUMCOEFFD(I) ,(vpCoeffsd(icoeff),icoeff=1,numcoeffd)  
		!indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		!read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
	!	cycle
!200		EXIT !terminate loop
!	enddo
!	close(41)
	nDeck=nDeck1 !+nDeck2

	!Write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	DO iComp=1,NC
		iGotIt=0
		DO J=1,NDECK
			IF(IDNUM(J).EQ.ID(iComp)) THEN
				rMINT(iComp)=rMINTD(J)
				VALMIN(iComp)=VALMIND(J)
				rMAXT(iComp)=rMAXTD(J)
				ID(iComp)=IDNUM(J)
				VALMAX(iComp)=VALMAXD(J)
				AVGDEV(iComp)=AVGDEVD(J)
				NUMCOEFF(iComp)=NUMCOEFFD(j)
				do k=1,5
					vpCoeffs(iComp,k)=vpCoeffsd(j,k)
				enddo
				iGotIt=1
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			write(*,*)'Error in GetCrit: not found for iComp=',iComp
			if(LOUD)pause
		endif
		!write(*,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
		!,TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
		
	enddo
	!do i=1,npts
		!T=TDAT(i)
		!do j=1,iComp
			!pSatExp(iComp)=exp(vpCoeffs(iComp,1)+vpCoeffs(iComp,2)/T+vpCoeffs(iComp,3)*DLOG(T)+vpCoeffs(iComp,4)*T**vpCoeffs(iComp,5))/1000000
		!enddo
	!enddo

	RETURN
861	continue
	write(*,*)'GetVp error - error reading CoeffsVp.dta. Path? Debug?'
	write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,I
	iErrCode=1
	if(LOUD)pause
	return                      
!862	continue
!	write(*,*)'GetVp error - error reading CoeffsVp.txt. Path?'
!	write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,I
!	iErrCode=1
!	pause
!	return                      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	integer function GetBIPs(bipFile,idComp,NC)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION (A-H,K,O-Z)
	PARAMETER (listPool=1000)
	logical switched
	character bipFile*88, dumString*88
	!integer GetBIPs
	dimension idBinarY(listPool),idComp(NC),alphaDB(listPool),KIJDB(listPool),KTIJDB(listPool),&
	wsTAUij (listPool),wsTAUji (listPool),wsTauTij(listPool),wsTauTji(listPool)
	common/ksvall/ks0ij,ks1ij
	dimension ks0ij(nmx,nmx),ks1ij(nmx,nmx)
	!data initial/0/
	initial=0 !just assume that the database needs to be reloaded if being called

	GetBIPs=0	
	if(initial.eq.0)then
		initial=1
		open(55,file=bipFile,ERR=861)
		read(55,'(a88)',ERR=861)dumString
		item=1
		do while(item.ge.0)
			read(55,*,ERR=861,end=100) idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),&
				wsTAUji(item),wsTauTij(item),wsTauTji(item),alphaDB(item)
			item=item+1	!increment for subsequent loop if any.
			cycle
100			exit
		enddo
		close(55)
	endif
	numBips=item
	if(LOUD)write(*,*)'# of bips in database =',numBips

	do iComp=1,NC
		do jComp=1,NC
			idBin=10000*idComp(iComp)+idComp(jComp)
			switched=.FALSE.
			if(idComp(iComp).lt.idComp(jComp))then
				switched=.TRUE.
				idBin=10000*idComp(jComp)+idComp(iComp)
			endif

			itemFound=0
			do item=1,numBips
				if(idBinarY(item).eq.idBin)itemFound=item
			enddo

			xsAlpha(iComp,jComp)=0.3d0
			KIJ (iComp,jComp)=0
			KTIJ(iComp,jComp)=0
			KS0IJ(iComp,jComp)=0
			KS1IJ(iComp,jComp)=0
			xsTau (iComp,jComp)=0
			xsTau (jComp,iComp)=0
			xsTauT(iComp,jComp)=0
			xsTauT(jComp,iComp)=0

			if(itemFound.ne.0)then
				xsAlpha(iComp,jComp)=alphaDB(itemFound)
				KIJ(iComp,jComp)=KIJDB(itemFound)
				KTIJ(iComp,jComp)=KTIJDB(itemFound)

				xsTau (iComp,jComp)=wsTAUij(itemFound)
				xsTau (jComp,iComp)=wsTAUji(itemFound)
				xsTauT(iComp,jComp)=wsTauTij(itemFound)
				xsTauT(jComp,iComp)=wsTauTji(itemFound)
				if(switched)then
					xsTau (jComp,iComp)=wsTAUij(itemFound)
					xsTau (iComp,jComp)=wsTAUji(itemFound)
					xsTauT(jComp,iComp)=wsTauTij(itemFound)
					xsTauT(iComp,jComp)=wsTauTji(itemFound)
				endif
			endif
			xsAlpha(jComp,iComp)=xsAlpha(iComp,jComp)
			KIJ(jComp,iComp)=KIJ(iComp,jComp)
			KTIJ(jComp,iComp)=KTIJ(iComp,jComp)
		enddo
	enddo

	return
861	continue
	if(LOUD)write(*,'(a33,a22)')'GetBip error - error reading ',bipFile
	if(LOUD)write(*,*)'numBips,item',numBips,item
	if(LOUD)write(*,*)'Data read as:'
	if(LOUD)write(*,'(i4,7F8.4)')idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),wsTAUji(item),wsTauTij(item),wsTauTji(item),alphaDB(item)
	GetBIPs= -1
	return                      
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdCcLookup(NC,idcc,ier,errMsgPas)
    USE GlobConst
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas
	character inFile*251
	dimension idcc(NC)
	dimension idCcDb(maxDb),idDb(maxDb)
	IF(DEBUG)THEN
		OPEN(50,FILE='c:\spead\calceos\input\idTrcDipCas.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT' ! // is the concatenation operator
		OPEN(50,FILE=inFile)
	ENDIF
	errMsg(0)='No Problem'
	errMsg(1)='IdCcLookup Error: at least one id not found'
	ier=0
	i=0
	id(1)=0
	do i=1,maxDb
		read(50,*,END=101)idCcDb(i),idDb(i)
		if(idCcDb(i).eq.idcc(1))id(1)=idDb(i)
		cycle
101		continue
		numDb=i-1
		exit !terminate do loop
	enddo
	if(id(1).eq.0)then
		ier=1
		write(*,*)'Did not find idcc(1)=',idcc(1)
		errMsgPas=errMsg(ier)
		return
	endif
	do iComp=2,NC
		id(iComp)=0
		iGotIt=0
		i=0
		do while(iGotIt.eq.0.and.i.lt.numDb)
			i=i+1
			if(idCcDb(i).eq.idcc(iComp))then
				id(iComp)=idDb(i)
				iGotit=1
			endif
		enddo
		if(id(iComp).eq.0)then		
			ier=1
			write(*,*)'Did not find idcc(i).i,idcc=',iComp,idcc(1)
			errMsgPas=errMsg(ier)
			return
		endif
	enddo
	close(50)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdCasLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas
	character inFile*251
	dimension idCas(NC)
	dimension idCasDb(maxDb),idDb(maxDb)
	inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT'
	if(DEBUG)inFile='c:\spead\calceos\input\idTrcDipCas.TXT'
	OPEN(50,FILE=inFile)
!	IF(DEBUG)THEN
!		OPEN(50,FILE='c:\spead\idTrcDipCas.TXT')
!	ELSE
!		inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT' ! // is the concatenation operator
!		OPEN(50,FILE=inFile)
!	ENDIF
	errMsg(0)='No Problem'
	errMsg(1)='IdCasLookup Error: at least one id not found'
	ier=0
	iGotIt=0
	read(50,*)nDeck
	!print*,'IdLookup: nDeck=',nDeck
	do i=1,nDeck
		read(50,*,ioStat=ioErr)idcc,idDb(i),idCasDb(i)
		if(ioErr .and. LOUD)print*,'i,idDippr=',i,idDb(i)
		if( id(1)==idDb(i) )then
			iGotIt=1
			idCas(1)=idCasDb(i)
		endif
	enddo
	if(iGotIt==0)then
		ier=1
		if(LOUD)pause 'IdCasLookup: idCas not found for comp1'
		goto 861
	endif
	if(NC==1)goto 861
	do iComp=2,NC
		if(LOUD)print*,'Looking up',iComp,'th compound. id=', id(iComp),' nDeck= ',nDeck
		notFound=1
		do i=1,nDeck
			if(idDb(i)==id(iComp))then
				idCas(iComp)=idCasDb(i)
				if(LOUD)print*,'iGotIt!,idCas=',idCasDb(i)
				notFound=0
				exit
			endif
		enddo
		if(notFound)ier=1 !only way here is if cycle failed => iGotIt==0.
		if(LOUD .and. ier)write(*,*)'Did not find idCas(i).i,id=',iComp,id(i)
		goto 861
	enddo
861	errMsgPas=errMsg(ier)
	close(50)
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine PrintErrFlash(iErrCode)
	USE GlobConst
	write(*,*)'PrintErrFlash:iErrCode=',iErrCode
	IF(iErrCode.EQ.1)WRITE(*,*)'LLeFL:Initial guess REQUIRED.'
	IF(iErrCode.EQ.2)WRITE(*,*)'Flash:ERROR ALL UPPER PHASE'
	IF(iErrCode.EQ.3)WRITE(*,*)'Flash:ERROR ALL LOWER PHASE'
	IF(iErrCode.EQ.4)WRITE(*,*)'Flash:LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.5)WRITE(*,*)'Flash:VAPOR  FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.6)WRITE(*,*)'Flash:ITMAX EXCEEDED		 '
	IF(iErrCode.EQ.7)WRITE(*,*)' Flash:xFrac(i) < 0 for some i'
	IF(iErrCode.EQ.11)WRITE(*,*)' zV~zL - Flash: warning'
	!IF(iErrCode.EQ.21)	!NOTE: not relevant for polymers because vLiq >> 1 and p >> pSat
	!& WRITE(*,*)' zL > 0.33 - Flash:warning '
	IF(iErrCode.EQ.22)WRITE(*,*)' zV < 0.33 - Flash:warning'
	IF(iErrCode.EQ.23)WRITE(*,*)' Flash warning:gEx < 0 on at least one LLE iteration.'
	!C
	if(LOUD)pause 'PrintErrFlash: Check errors'
	return
	end


	!**********************************************
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!C      AARON'S SUBROUTINE
	!C      $GETCRIT.FOR   VERSION 1.4
	SUBROUTINE GETCRITCAS(NC,iErrCode)	!GetCritCas assumes ID(GlobConst)=IdCas

	!C  
	!C  PROGRAMMED BY:  A.S. PUHALA - AUG 93
	!C  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	!C  REVISION DATE:  2/02 - FOR binary data file
	!C 
	!C  LOOKS UP THE CRITICAL PROPERTIES AND RETURNS JUST TC,PC,ACEN,NAME
	!C
	!C  INPUT
	!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!C  OUTPUT
	!C    TC - CRITICAL TEMPERATURE
	!C    PC - CRITICAL PRESSURE
	!C    ACEN - ACENTRIC FACTOR
	!C    NAME - COMPONENT NAME
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	CHARACTER NAMED*30,tCode*4,pCode*4,vCode*4,form*12
	character readText*132,dumText*77 
	PARAMETER(ndb=1555)
	character infile*251
	DIMENSION IDNUM(ndb),TCD(ndb),PCD(ndb),ACEND(ndb),NAMED(ndb),iCasd(ndb)
	DIMENSION ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb)
	common/FloryWert/solParm(nmx),vLiq(nmx),vMolec(NMX),&
		eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	iErrCode=0
	inFile=TRIM(masterDir)//'\ParmsCrit.txt' ! // is the concatenation operator
	IF(DEBUG)inFile='c:\SPEAD\CalcEos\input\ParmsCrit.txt' 
	print*,'CritFile=',TRIM(inFile)
	!OPEN(40,FILE=inFile)
	OPEN(40,FILE='ParmsCrit.txt')
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
	READ(40,*,ERR=861)NDECK1
	print*,'nDeck1=',nDeck1
	print*,'ID(1)=',ID(1)
    if(LOUD)then
	    if(ndeck1.gt.ndb)pause 'GetCrit: more data in file than allocated'
    end if
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
		READ (40,*,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCasd(i),tCode,pCode,vCode,form,NAMED(I)
!		write(61,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
!			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
!	close(61)
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1).EQ.0)THEN
		DO I=1,NDECK1,3
			WRITE(*,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
            if(LOUD)then
			    IF(((I-1)/45)*45.EQ.(I-1))PAUSE
            end if
		enddo
		RETURN
	ENDIF
	CLOSE(40)
    if(LOUD)then
	    if((nDeck).gt.ndb)pause 'GetCrit: too much data in file'
    end if

	IF(DEBUG)THEN
		OPEN(40,FILE='c:\spead\calceos\input\ParmsCrAdd.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	ENDIF
	READ(40,*,ERR=862)dumText

	nDeck2=0
	DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ(40,'(A132)',ERR=862,END=200)readText
		nDeck2=nDeck2+1
		I=nDeck1+nDeck2
		read(readText,*,ERR=862)IDNUM(I),TCD(I),PCD(I),ZCD(I) &
			,ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCasd(I)
		indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
		cycle
200		EXIT !terminate loop
	enddo
	close(40)
	nDeck=nDeck1+nDeck2
    if(LOUD)then
	    if((nDeck).gt.ndb)pause 'GetCrit: too much data in ParmsCrAdd file'
    end if
	Write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	DO iComp=1,NC
		iGotIt=0
		DO J=1,NDECK
			IF(ICasd(J).EQ.ID(iComp)) THEN
				NAME(iComp)=NAMED(J)
				TC(iComp)=TCD(J)
				PC(iComp)=PCD(J)
				ID(iComp)=IDNUM(J)
				ACEN(iComp)=ACEND(J)
				ZC(iComp)=ZCD(J)
				rMwPlus(iComp)=rMwD(j)
				solParm(iComp)=solParmD(j)
				vLiq(iComp)=vLiqD(j)
				iGotIt=1
				exit
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			write(*,*)'Error in GetCrit: not found for iComp=',iComp
			if(LOUD)pause
		endif
		write(*,'(1x,i11,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
		,TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	write(*,*)'GetCrit error - error reading ParmsCrit.dta. Path? Debug?'
	write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	if(LOUD)pause
	return                      
862	continue
	write(*,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	if(LOUD)pause
	return                      
	END

