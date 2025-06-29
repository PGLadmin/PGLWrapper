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
	USE Assoc !to display assoc parms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	CHARACTER*1 ANSWER
	CHARACTER*2 calcType
	CHARACTER*77 errMsgPas
	CHARACTER*251 outFile
    LOGICAL RELEASE
	!CHARACTER*15 EosName(nEOSs) ! moved to GlobConst

	DIMENSION BIPTMP(NMX)
	DIMENSION ZFEED(NMX),S(NMX)
    !INTEGER ier(12) !,idCas(NMX)
	Integer QueryNparMix	
	COMMON/CrOpt/calcType
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/FEED/ZFEED
	COMMON/KVALUES/S
	!common/FloryWert/vLiq(nmx)
	!cf. GlobConst  1     2       3       4          5          6         7           8              9        10       11      12       13          14         15           16            17        18		19       20
	!data EosName/'PR','ESD96','PRWS','ESD-MEM2','SPEADMD','Flory-MEM2','NRTL','SpeadGamma-MEM2','SPEAD11','PcSaft','tcPRq','GCESD','GcEsdTb','TransSPEAD','GcPcSaft','GcPcSaft(Tb)','tcPR-GE(W)','ESD2','LsgMem2','SptPcSaft'/
	
	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
    !print*,'testing'
	iChange =VERIFY('c:\msdev\Projects\calceos',CURDIR)
	iChange2=VERIFY('C:\MSDEV\Projects\CalcEos',CURDIR)
	if(iChange2==0)iChange=0
	if(iChange==0 .or. iChange > 26)DEBUG=.TRUE.
	DEBUG=.FALSE.   !setting to true mostly directs input from c:\spead...	 Logic above automatically updates to TRUE if CURDIR=c:\msdev\Projects\calceos
	!DEBUG=.TRUE.    ! Logic not working JRE 20200407
    !print*,'DEBUG=',DEBUG
														   
	masterDir=TRIM(curDir)
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
    !print*,'KijOutFile=',TRIM(outFile)
	open(655,file=outFile)
	write(655,*)' ' ! clear out old contents.
	close(655)
	outFile=TRIM(masterDir)//'\output\output.txt' ! // is the concatenation operator
    !print*,'outFile=',TRIM(outFile)
	LOUD =.FALSE.
	!LOUD = .TRUE.
    RELEASE=.FALSE.
    !RELEASE=.TRUE.
	PGLInputDir='c:\PGLWrapper\input'
    IF(RELEASE)PGLInputDir=TRIM(masterDir)//'\input'
	dumpUnit=6
	if(dumpUnit/=6)then
		dumpUnit=686
		outFile=TRIM(masterDir)//'\output\JreDebugDLL.txt'
		open(dumpUnit,file=outFile)
    endif
    !print*,'dumpUnit=',dumpUnit
	write(dumpUnit,*)'MasterDir=',TRIM(masterDir)
	write(dumpUnit,*)'PGLInputDir=',TRIM(PGLInputDir)
    !print*,'Opening outFile=',TRIM(outFile)
	OPEN(UNIT=52,FILE=outFile)
    iErr=0
	call LoadCritParmsDb(iErr)
    !print*,'Main: from LoadCritParmsDb... iErr=',iErr
	if(iErr /=0)then
		write(*,*)'Main: LoadCritParmsDb() -> iErr=',iErr
		pause 'Did you add a new compound and forget to update nCritSet?'
		stop
	endif
    if(iErr/=0)pause 'Main: LoadCritParmsDb failed'
    call GetVpDb(iErr)
    if(iErr/=0)pause 'Main: GetVpDb failed'
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! AFG 2011: Option 9 is added
	write(*,*)'Enter EOS option: '
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=1,7)
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=8,13)
	write(  *,'( 1x,11(i2,a,a,1x) )'  )( i,'=',TRIM(EosName(i)),i=14,nModels)
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
	Call PGLStartup(NC,iEosOpt,1,iErrStart)	! idOpt=1 means ID(dippr) is USEd.
	if(iErrStart.ne.0)then
		write(*,*)'Main: ErrorCode from PGLStartup = ',iErrStart
		if(LOUD)pause
		goto 86
	endif
	if(bTPT)print*, 'Main: etaMax(TPT)=',etaMax
	write(*,*)' ID     Tc     Pc     w    bVol    Name'
	do i=1,nc
		write(*,'(i5,f8.2,f8.4,f7.4,f8.2,1x,a)') ID(i),Tc(i),Pc(I),acen(i),bVolCc_mol(i),TRIM(Name(i))
	enddo
	if(bESD.or.bTPT.or.iEosOpt==19)then !display hb parms
		write(*,*)' ID  IdType  nDeg  KAD   nDon eDon   nAcc eAcc    '
		do i=1,NC
			do j=1,nTypes(i)
				write(*,622)ID(i),idType(i,j),nDegree(i,j),bondVolNm3(i,j),nDonors(i,j),eDonorKcal_mol(i,j),nAcceptors(i,j),eAcceptorKcal_mol(i,j)
			enddo
		enddo
	endif !display hb parms
	if(bESD)then
		write(*,*)' ID      q     eps/kB    bVol'
		do i=1,NC
			write(*,form601)ID(i),q(i),eokP(i),bVolCc_mol(i)
		enddo
	endif !display hb parms
622 format(1x,3i6,f8.5,i3,f8.4,i3,f8.4) 
	pause 'Main: Check basic info.'   
	IF(DEBUG)write(*,*)'Main: Kij(1,2) = ',Bip12
	if(DEBUG)WRITE(*,*)'THE KIJ MATRIX AT 298K IS ESTIMATED AS'
	if(DEBUG)WRITE(*,'(8X,9I8)')(ID(I),I=1,NC)
	nParms= QueryNparMix()	!NOTE: calling routine must declare "integer QueryNparMix"
	call QueryParMix(1,Bip12,iErrBip) ! Query what was stored during the Get__ function.
	if(nParms > 1)call QueryParMix(1,Bip21,iErrBip) ! Query what was stored during the Get__ function.
	DO I=1,NC
		DO J=1,NC
			BIPTMP(J)=KIJ(I,J)+KTIJ(I,J)/298
			if(nParms > 1)BIPTMP(J)=xsTau(i,j)
		ENDDO
		if(nParms==1.or.nParms==3)then
			WRITE(6,601)I,ID(I),(BIPTMP(J),J=1,NC)
		elseif(nParms>1)then
			WRITE(6,602)I,ID(I),(BIPTMP(J),J=1,NC)
		endif
	enddo
601	FORMAT(I3,I5,9(1X,F7.4))
602	FORMAT(I3,I5,9(1X,F9.1))
!	if(LOUD)WRITE(*,*)'THE HIJ MATRIX AT 298K IS ESTIMATED AS'
!	if(LOUD)WRITE(*,'(8X,9I8)')(ID(I),I=1,NC)
!	DO I=2,NC
!		DO J=1,I-1
!			BIPTMP(J)=HIJ(I,J)+HTIJ(I,J)/298
!		ENDDO
!		if(LOUD)WRITE(6,601)I,ID(I),(BIPTMP(J),J=1,I-1)
!	enddo
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
		WRITE(6,*)'DF FOR Diffusion Coefficient Database'
		WRITE(6,*)'FU FOR FUGACITY CALCULATION Given T, P'
		WRITE(6,*)'FV FOR FUGACITY CALCULATION GIVEN T,V'
		WRITE(6,*)'CR FOR MIXTURE CRITICAL POINT CALCULATION-can be used for Options 4 & 5'		  !Added by AGF, Dec. 2009
		WRITE(6,*)'CP FOR PURE COMPONENT CRITICAL POINT CALCULATION-can be used for Options 4 & 5'	  !Added by AGF, Dec. 2009
		WRITE(6,*)'CD FOR PURE COMPONENT CRITICAL POINTS OF DATABASE'	  !Added by JRE, FEB. 2023
		WRITE(6,*)'GE FOR Excess Gibbs and Enthalpy CALCULATION '
		WRITE(6,*)'GV FOR Excess Gibbs at T,Vtot '
		WRITE(6,*)'GX FOR Excess Gibbs at T,eta '
		WRITE(6,*)'IB FOR isoBar plot on single component system'
		WRITE(6,*)'IC FOR isochore plot on single component system'
		WRITE(6,*)'IT FOR isotherm plot on single component system'	  !Modified by AGF, Dec. 2009
		WRITE(6,*)'HD FOR HE calculation for multi-system VLE Database'
		WRITE(6,*)'KD FOR Kij regression on multi-system VLE Database'
		WRITE(6,*)'KI FOR Kij ITERATION ON BINARY AT SINGLE POINT'
		WRITE(6,*)'KL FOR Kij regression on multi-system LLE Database'
		WRITE(6,*)'KN FOR NEW BIP INPUT'
		WRITE(6,*)'KO FOR SINGLE KIJ ITERATION FOR FILE OF EXPTL DATA'
		WRITE(6,*)'KS FOR Kij regression on multi-system SLE Database'
		WRITE(6,*)'LL FOR LLE   '
		WRITE(6,*)'MS FOR MemSced Optimization.   '
		WRITE(6,*)'OB FOR OPTIMUM BIPS WITH MULTIPLE BIPS'
		WRITE(6,*)'PE FOR PHASE ENVELOPE (P,T given X,Y)'
		WRITE(6,*)'PS FOR POLYMER-SOLVENT PARTITIONING ESTIMATION '
		WRITE(6,*)'PT FOR Pure component Property Table '
		WRITE(6,*)'PX FOR P,X,Y given T'
		WRITE(6,*)'RP FOR REGRESSION OF PURE COMPONENT ESD PARAMETERS'
		WRITE(6,*)'RS FOR Kii regression on temperature-vapor pressure Database for SPEAD EOS'
		WRITE(6,*)'RX FOR REGRESSION OF XA,XD (e.g. from FTIR data)'
		WRITE(6,*)'SC FOR SUPERCRITICAL FLUID SOLID SOLUBILITY'
		WRITE(6,*)'TH FOR ThermoCheck using numerical derivatives of FUGI() outputs  '
		WRITE(6,*)'TP FOR 3-PHASE CALCULATION  '
		WRITE(6,*)'TX FOR T,X,Y given P'
		WRITE(6,*)'VL FOR VLE FLASH   '
		WRITE(6,*)'VP FOR vapor pressure   '
		WRITE(6,*)'FC FOR FREEZING CURVE CALC 4 CRYSTAL solid' 
		WRITE(6,*)'FO FREEZCURVE KIJ optimization'	
		!WRITE(6,*)'EA FOR Equal Area rule calculation (vapor pressure)'
		!WRITE(6,*)'FW FOR FOR FITTING THE WHITE RG METHOD PARAMETER'
		WRITE(6,*)'VC FOR VIRIAL COEFFICIENTS FOR PURE COMPONENTS, EOS OPTIONS: 4 & 5'  !AFG 2011
		WRITE(6,*)'QT TO QUIT'

		WRITE(6,*)'Note: previous output files may be deleted after your selection. '
		READ(5,'(A2)')calcType
10		continue	! transfer in.. e.g. if compd is missing from dbase for ESD, then RP.
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
			open(671,file=outFile)
			write(671,*)' ' !open and clear old contents
			close(671)
		endif

		IF(calcType=='AC'.OR.calcType=='ac')CALL IDACs(NC)
		IF(calcType=='BD'.OR.calcType=='bd')CALL BifurcationDiagram(NC)
		IF(calcType=='BI'.OR.calcType=='bi')CALL BINODAL(NC)
		IF(calcType=='BF'.OR.calcType=='bf')CALL BpFromFile(NC)
		IF(calcType=='BP'.OR.calcType=='bp')CALL BPITER(NC)	   !In the future we can add isZiter to other options if it is necessary, AUG 10
		IF(calcType=='BT'.OR.calcType=='bt')CALL BTITER(NC)
		IF(calcType=='DT'.OR.calcType=='dt')CALL DTITER(NC)
		IF(calcType=='DF'.OR.calcType=='df')CALL DifCoDb
		IF(calcType=='FC'.OR.calcType=='fc')CALL FREEZINGCURVECALC(NC)
		IF(calcType=='FO'.OR.calcType=='fo')CALL FreezOptKij(NC)
		IF(calcType=='FU'.OR.calcType=='fu')CALL FITER(NC)
		IF(calcType=='FV'.OR.calcType=='fv')CALL FVITER(NC)
		IF(calcType=='CR'.OR.calcType=='cr')CALL CRITER(NC)	                !Added by AGF, Dec. 2009	AUG 10
		IF(calcType=='CP'.OR.calcType=='cp')CALL CPITER(NC) 	            	!Added by AGF, Dec. 2009	AUG 10
		IF(calcType=='CD'.OR.calcType=='cd')CALL CritDB(iErrCode) 	            	!Added by AGF, Dec. 2009	AUG 10
   		IF(calcType=='GE'.OR.calcType=='ge')CALL GibbsXS(NC)
   		IF(calcType=='GV'.OR.calcType=='gv')CALL GibbsXsVtot(NC)
		IF(calcType=='GX'.OR.calcType=='gx')CALL GibbsXsEta(NC)
		IF(calcType=='IB'.OR.calcType=='ib')CALL Isobar(NC,iErrCode,errMsgPas)
		IF(calcType=='IC'.OR.calcType=='ic')CALL Isochore(NC,iErrCode,errMsgPas)
		IF(calcType=='IT'.OR.calcType=='it')CALL Isotherm(NC)				!Modified by AGF, Dec. 2009
		IF(calcType=='HD'.OR.calcType=='hd')CALL HEDB(NC)
		IF(calcType=='KI'.OR.calcType=='ki')CALL KIJITR(NC)
		IF(calcType=='KD'.OR.calcType=='kd')CALL KIJDB(NC)
		IF(calcType=='KL'.OR.calcType=='kl')CALL KIJDBLLE(NC)
		IF(calcType=='KS'.OR.calcType=='ks')CALL KIJDBSLE(NC)
		IF(calcType=='KN'.OR.calcType=='kn')call BipIo(NC)
		IF(calcType=='KO'.OR.calcType=='ko')CALL KIJOPT(NC)
		IF(calcType=='LL'.OR.calcType=='ll')CALL LLITER(NC)
		IF(calcType=='MS'.OR.calcType=='ms')CALL MemScedOpt(NC)
		IF(calcType=='OB'.OR.calcType=='ob')CALL OptiBIPs(NC)
		IF(calcType=='PE'.OR.calcType=='pe')CALL PhasEnv(NC)
		IF(calcType=='PS'.OR.calcType=='ps')CALL PSITER(NC,iErrCode)
		IF(calcType=='PT'.OR.calcType=='pt')CALL PropTable(NC,iErrCode)
		IF(calcType=='PX'.OR.calcType=='px')CALL PXYT(NC)
		IF(calcType=='RP'.OR.calcType=='rp')CALL RegPureEsd2(NC)
		IF(calcType=='RS'.OR.calcType=='rs')CALL RegSpeadIo(NC)
		IF(calcType=='RX'.OR.calcType=='rx')CALL RegXaXdIo(NC)
		IF(calcType=='SC'.OR.calcType=='sc')CALL SCFITER(NC)
		IF(calcType=='TH'.OR.calcType=='th')CALL ThermoCheck(NC)
		IF(calcType=='TP'.OR.calcType=='tp')CALL TEITER(NC)
		IF(calcType=='TX'.OR.calcType=='tx')CALL TXYP(NC)
		IF(calcType=='VL'.OR.calcType=='vl')CALL VLITER(NC)
		IF(calcType=='VP'.OR.calcType=='vp')CALL VpIter(NC) !,gmol,tKelvin,Psat,vLCc,vVCc)
		!IF(calcType=='FW'.OR.calcType=='fw')CALL FitWhite(NC)				  	!AFG 2011
		IF(calcType=='VC'.OR.calcType=='vc')CALL VirialCoeff(NC)				!AFG 2011
		IF(calcType=='KS'.OR.calcType=='ks')CALL RegA0A1A2TptIo(NC)
		!IF(calcType=='IT'.OR.calcType=='it')CALL Isotherm(NC,ier,errMsnPas)
		IF(calcType=='QT'.OR.calcType=='qt')iQuit=1
	enddo  !while(iQuit.ne.1)

86	continue
	close(52)
	stop
	END


