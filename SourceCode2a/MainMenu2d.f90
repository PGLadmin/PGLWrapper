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
	parameter (nEOSs=17)

	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	CHARACTER*1 ANSWER
	CHARACTER*2 calcType
	CHARACTER*77 errMsgPas
	CHARACTER*251 outFile
	CHARACTER*12 EosName(nEOSs)

	DIMENSION BIPTMP(NMX)
	DIMENSION ZFEED(NMX),S(NMX)
    INTEGER ier(12) !,idCas(NMX)
	
	COMMON/CrOpt/calcType
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/FEED/ZFEED
	COMMON/KVALUES/S
	!common/FloryWert/vLiq(nmx)
	data EosName/'PR','Esd96','PRWS','Esd','SPEAD','FloryWert','NRTL','SpeadGamma','SPEAD11','PcSaft','tcPRq','GCESD','GcEsdTb','TransSPEAD','GcPcSaft','GcPcSaft(Tb)','tcPR-GE(W)'/

	LOUD =.FALSE.
	!LOUD = .TRUE.
	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	iChange =VERIFY('c:\msdev\Projects\calceos',CURDIR)
	iChange2=VERIFY('C:\MSDEV\Projects\CalcEos',CURDIR)
	if(iChange2==0)iChange=0
	if(iChange==0 .or. iChange > 26)DEBUG=.TRUE.
	DEBUG=.FALSE.   !setting to true mostly directs input from c:\spead...	 Logic above automatically updates to TRUE if CURDIR=c:\msdev\Projects\calceos
	DEBUG=.TRUE.    ! Logic not working JRE 20200407

	masterDir=TRIM(curDir)
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(61,file=outFile)
	write(61,*)' ' ! clear out old contents.
	close(61)
	outFile=TRIM(masterDir)//'\output\output.txt' ! // is the concatenation operator

	if(LOUD)print*,'MasterDir=',TRIM(masterDir)
	OPEN(UNIT=52,FILE=outFile)
	call LoadCritParmsDb(iErr)
    if(iErr/=0)pause 'Main: LoadCritParmsDb failed'
    call GetVpDb(iErr)
    if(iErr/=0)pause 'Main: GetVpDb failed'
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! AFG 2011: Option 9 is added
	write(*,*)'Enter EOS option: '
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=1,7)
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=8,13)
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=14,nEOSs)
	read(*,*)iEosOpt
	write(52,*)'EOS=',iEosOpt,' ',EosName(iEosOpt)
	if(iEosOpt==9) then
		Write(*,*)' Please note that SPEAD11 is a correlation for n-alkanes > C7'
		Write(*,*)' So, using it for other components will result in nonsense results. '
		Write(*,*)' Also, its applicability to mixtures has not been tested, yet.'
		Write(*,*)' AFG'
	endif
	if(iEosOpt==17)print*,'Note: BIPs for tcPR-GE(W) are stored in ...EtcFiles\Param_Wilson.txt'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!NewEos: Add here for listing new iEosOpt.
	NC= -1
	do while(NC < 0 .or. NC > 55)
		WRITE(6,*)'ENTER NUMBER OF COMPONENTS (0 TO QUIT)   '
		READ(5,*)NC
		WRITE(52,*)' NUMBER OF COMPONENTS (0 TO QUIT)   '
		write(52,*)NC
		IF(NC==0)goto 86
		IF(NC < 0 .or. NC > 55)then
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
	if(iEosOpt.eq.2)CALL GetEsdCas(NC,idCas,iErrGet)	 !Results placed in USE EsdParms				!Diky model 12
	if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet) 
	if(iEosOpt.eq.4)CALL GetEsdCas(NC,idCas,iErrGet)
	if(iEosOpt.eq.5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms			  !Diky model 6
	if(iEosOpt.eq.6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	if(iEosOpt.eq.8)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.9)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	!if(iEosOpt.eq.9)CALL GetRG(NC,ID,iErrGet)		!AFG 2011 : Reads White's RG method parameters, disabled 20191121 by JRE. unused.
	if(iEosOpt.eq.10)CALL GetPcSaft(NC,idCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters   !Diky model 21
	if(iEosOpt.eq.11)CALL GetPrTc(NC,iErrGet) 	!JRE 2019 : Reads Jaubert's pure parameters				!Diky model 22
	if(iEosOpt.eq.12)CALL GetEsdCas(NC,idCas,iErrGet)	 !ParmsEsdEmami in USE EsdParms			  !Diky model 23
	if(iEosOpt.eq.13)CALL GetEsdCas(NC,idCas,iErrGet)	 !Results placed in USE EsdParms					 !Diky model 17
	if(iEosOpt.eq.14)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!ParmsTptTransferable placed in common:      			  !Diky model 18 
	if(iEosOpt.eq.15)CALL GetPcSaft(NC,idCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters   !Diky model 21
	if(iEosOpt.eq.16)CALL GetPcSaft(NC,idCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters   !Diky model 2_
	if(iEosOpt.eq.17)CALL GetPrLorraine(NC,iErrGet)		!JRE 20210510 : Reads Jaubert's tcPR parameters including BIPs for GE mixing rule.   !Diky model 
	!NewEos: Add here for initializing parms.
	if(iErrGet > 0)then
		if(LOUD)pause  'Error in Main: failed to get required dbase props.'
		if(isESD)then
			print*, 'Main: iErrGetEos=',iErrGet
			pause 'Main: check error.'
			if(NC > 1)goto 86
			calcType='RP'
			goto 10 
		elseif(isTPT)then
			print*, 'Main: ',TRIM(errMsgPas)
			pause 'Main: check error.'
		else
			print*, 'Main: iErrGetEos=',iErrGet
			pause 'Main: check error.'
			goto 86
		endif
	endif
	call QueryParMix(1,Bip12,iErrBip) ! Query what was stored during the Get__ function.
	write(*,*)' ID     Tc     Pc     w    bVol    Name'
	do i=1,nc
		write(*,'(i5,f8.2,f8.4,f7.4,f8.2,1x,a33)') ID(i),Tc(i),Pc(I),acen(i),bVolCc_mol(i),TRIM(Name(i))
	enddo    
	IF(DEBUG)write(*,*)'Main: Kij(1,2) = ',Bip12
	if(DEBUG)WRITE(*,*)'THE KIJ MATRIX AT 298K IS ESTIMATED AS'
	if(DEBUG)WRITE(*,'(8X,9I8)')(ID(I),I=1,NC)
	DO I=2,NC
		DO J=1,I-1
			BIPTMP(J)=KIJ(I,J)+KTIJ(I,J)/298
		ENDDO
		if(DEBUG)WRITE(6,601)I,ID(I),(BIPTMP(J),J=1,I-1)
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
		WRITE(6,*)'AC FOR INFINITE DILUTION ACTIVITY COEFFICIENTS   '
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
		WRITE(6,*)'IB FOR isoBar plot on single component system'
		WRITE(6,*)'IC FOR isochore plot on single component system'
		WRITE(6,*)'IT FOR isotherm plot on single component system'	  !Modified by AGF, Dec. 2009
		WRITE(6,*)'KD FOR Kij regression on multi-system Database'
		WRITE(6,*)'KI FOR Kij ITERATION ON BINARY AT SINGLE POINT'
		WRITE(6,*)'KL FOR Kij regression on multi-system Database'
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
10    continue	! transfer in.. e.g. if compd is missing from dbase for ESD, then RP.
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

		IF(calcType.EQ.'AC'.OR.calcType.EQ.'ac')CALL IDACs(NC)
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
		IF(calcType.EQ.'IB'.OR.calcType.EQ.'ib')CALL Isobar(NC,iErrCode,errMsgPas)
		IF(calcType.EQ.'IC'.OR.calcType.EQ.'ic')CALL Isochore(NC,iErrCode,errMsgPas)
		IF(calcType.EQ.'IT'.OR.calcType.EQ.'it')CALL Isotherm(NC)				!Modified by AGF, Dec. 2009
		IF(calcType.EQ.'KD'.OR.calcType.EQ.'kd')CALL KIJDB(NC)
		IF(calcType.EQ.'KI'.OR.calcType.EQ.'ki')CALL KIJITR(NC)
		IF(calcType.EQ.'KL'.OR.calcType.EQ.'kl')CALL KIJDBLLE(NC)
		IF(calcType.EQ.'KS'.OR.calcType.EQ.'ks')CALL KIJDBSLE(NC)
		IF(calcType.EQ.'KN'.OR.calcType.EQ.'kn')call BipIo(NC)
		IF(calcType.EQ.'KO'.OR.calcType.EQ.'ko')CALL KIJOPT(NC)
		IF(calcType.EQ.'LL'.OR.calcType.EQ.'ll')CALL LLITER(NC)
		IF(calcType.EQ.'OB'.OR.calcType.EQ.'ob')CALL OptiBIPs(NC)
		IF(calcType.EQ.'PE'.OR.calcType.EQ.'pe')CALL PhasEnv(NC)
		IF(calcType.EQ.'PS'.OR.calcType.EQ.'ps')CALL PSITER(NC,iErrCode)
		IF(calcType.EQ.'PT'.OR.calcType.EQ.'pt')CALL PropTable(NC,ier)
		IF(calcType.EQ.'PX'.OR.calcType.EQ.'px')CALL PXYT(NC)
		IF(calcType.EQ.'RP'.OR.calcType.EQ.'rp')CALL RegPureEsd2(NC)
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


