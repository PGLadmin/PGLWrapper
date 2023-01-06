	!Subroutine PGLWrapper()
	!MAIN PROGRAM FOR CALLING EOS SUBROUTINES. Includes GetCrit, GetBips, BipIo.
	!	SETS UP THE CRITS, EOSPARMS, bips, DEBUG status, FUGI(iEosOpt)
	!	Echoes user IO to Output.txt, and reports error checks. 
	!	INITIALIZATION AND CALLING SEQUENCE FOR VLE, LLE, VLLE SUBROUTINES.
	!reqd routines:
	!	bubpl.for, bubtl.for, dewtv.for, flashsub.for, FuEsdMy.for FuEsdXs2.for, FugiPr.f90, FugiPrws.for, RegPure.f90
	!	LmDifEzCov2.for, Mintools.for
	!           1     2       3     4      5        6         7        8           9        10       11      12       13          14         15           16            17            18
	! EosName/'PR','Esd96','PRWS','Esd','SPEAD','FloryWert','NRTL','SpeadGamma','SPEAD11','PcSaft','tcPRq','GCESD','GcEsdTb','TransSPEAD','GcPcSaft','GcPcSaft(Tb)','tcPR-GE(W)','ESD-MEM2'/
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE PORTLIB
	USE GlobConst
	USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	!CHARACTER*1 ANSWER
	!CHARACTER*2 calcType
	CHARACTER*77 errMsgPas,readString,Property(22)
	CHARACTER*251 outFile
	CHARACTER*3 aPhase

	!DIMENSION BIPTMP(NMX)
	!DIMENSION ZFEED(NMX),S(NMX)
    DoublePrecision xFrac(NMX),FUGC(NMX),xPure1(NMX),xPure2(NMX) !FUGI requires mole fraction specification because it is written generally for mixtures.
    INTEGER ier(12), localCas(NMX)
	
	COMMON/eta/etaL,etaV,ZL,ZV

	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	!iChange=VERIFY('C:\MYPROJEX\CalcEos',CURDIR)
	!iChange2=VERIFY('C:\MYPROJEX\CALCEOS',CURDIR)
	!if(iChange2.eq.0)iChange=0
	!if(iChange.eq.0 .or. iChange.gt.26)DEBUG=.TRUE.
	masterDir=TRIM(curDir)
	outFile=TRIM(masterDir)//'\output.txt' ! // is the concatenation operator
	DEBUG=.TRUE.
	DEBUG=.FALSE. ! read input files from c:\Projects\...\input dir.
    LOUD = .TRUE.
    LOUD = .FALSE.
	OPEN(UNIT=52,FILE=outFile,RECL=8192)
    open(UNIT=51,FILE='input.txt')
    NC=1 !assume one component for all calculations (as of 1/1/2020 this is all we need).
    xFrac(1)=1  !   "
	read(51,'(a77)')readString
	if(LOUD)print*,'readString: ',readString
    read(readString,*)iEosOpt,iProperty,localCas(1)
	if(iProperty > 10) then
		read(readString,*)iEosOpt,iProperty,localCas(1),localCas(2) !,xFrac(1)
		NC=2
		if(LOUD)print*,'Mixture calculation. id1,id2',localCas(1),localCas(2)
	endif
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (g/cc) given tKelvin
    !iProperty = 3: vapor fluid density (g/cc) given tKelvin, pKPa, iPhase (=0 for vapor)
    !iProperty = 4: liquid fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid)
    !iProperty = 5: critical fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 6: sound speed
    !iProperty = 7: Hres/RT,CpRes_R/R,CvResR,cmprsblty given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)
	!iProperty = 11: Vex (cc/g),Hex(J/mo)
    Property(1)='  vapor pressure (kPa) given tKelvin '
    Property(2)='  saturated liquid density (g/cc) given tKelvin '
    Property(3)='  vapor fluid density (g/cc) given tKelvin, pKPa, iPhase (=0 for vapor) '
    Property(4)='  liquid fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid) '
    Property(5)='  critical fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)'
    Property(6)='  sound speed '
    Property(7)='  Hres/RT,CpRes_R/R,CvResR,cmprsblty given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)	'
	Property(11)='  Vex (cc/g),Hex(J/mo) '
	if(LOUD)print*,TRIM(Property(iProperty))

    !iEosOpt=10 !Temporary for testing

	INITIAL=0
	ierCode=8686 ! Initialize to failure in the iEosOpt write slot to make it easier to terminate if any Get_ functions fail.
	call LoadCritParmsDb(iErrCrit)
	if(iErrCrit > 0 .and. LOUD)print*,'PGLWRapperMain: From LoadCritParmsDb, iErrCrit=',iErrCrit
	call IdDipprLookup(NC,localCas,iErrCas,errMsgPas)
	if(LOUD)print*,'localCas,idDippr=',(localCas(i),id(i), i=1,NC)
	if(iErrCas)then
        if(LOUD)write(6 ,*)' Sorry, must abort.  Not found for CAS number(s)=', (localCas(i),i=1,NC)
        write(52,*)ierCode,' Sorry, must abort.  Not found for CAS number(s)=', (localCas(i),i=1,NC)
        goto 86
    endif
	if(LOUD)print*,'idDippr,class=',ID(1),TRIM(class(1))
	CALL GETCRIT(NC,iErrCrit)
	if(iErrCrit)then
		if(LOUD)write(6 ,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		write(52,*)ierCode,' Error in Main: ErrorCode from GetCrit = ',iErrCrit
		goto 86
	endif
	if(LOUD)print*,'calling GetEOS'
	iErrGet=0 
	if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	if(iEosOpt.eq.2)CALL GetEsd96Cas(NC,localCas,iErrGet)	 !Results placed in USE EsdParms					!Diky model 12
	if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet) 
	if(iEosOpt.eq.5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms					!Diky model 6
	if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	if(iEosOpt.eq.8)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.9)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.10)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters			!Diky model 25
	if(iEosOpt.eq.11)CALL GetPrTc(NC,iErrGet)		!JRE 2019 : Reads Jaubert's parameters						  !Diky model 22
	if(iEosOpt.eq.12)CALL GetEsd96Cas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmami					 !Diky model 23
	if(iEosOpt.eq.13)CALL GetEsd96Cas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmamiTb					 !Diky model 24
	if(iEosOpt.eq.14)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in common: TptParms, HbParms				    !Diky model 18
	if(iEosOpt.eq.15)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Emami's PcSaft parameters			!Diky model 26
	if(iEosOpt.eq.16)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Emami's(Tb) PcSaft parameters			!Diky model 27
	!if(iEosOpt.eq.17)This is the model number for COSMOtherm, if one is required JRE 20200119.		!Diky model ??
	!NewEos: Add here for initializing parms.
	if(iErrGet .ne. 0)then
		if(LOUD)write(* ,*)ierCode,' Error in Main: failed to get required dbase props.'
		write(52,*)ierCode,' Error in Main: failed to get required dbase props.'
		goto 86
	endif
    if(LOUD)write(*,*)'PGLWRapperMain: Success so far... Got EOS Parms for iEosOpt=',iEosOpt
	if(LOUD)pause 'PGLWRapperMain: starting calcs'
	if(iProperty==1)then
		write(52,'(2i4,i11,a)')iEosOpt, 1, localCas(1),' = iEosOpt, iProp, localCas: T, P,expVal,pSatKPa,ierCode '
	elseif(iProperty > 1 .and. iProperty < 6)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,rho(g/cc),ierCode '
	elseif(iProperty>5 .and. iProperty < 11)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode '
	else
		write(52,'(i4,a,5i11)')iEosOpt,' = iEosOpt. localCas()= ',(localCas(i),i=1,NC)
	endif
	ierCode=0

	errDummy=86.86D0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
    if(iProperty==1 .or. iProperty==2)then !calculate psat or rhoSatLiq given T(K)
        notDone=1
        line=0  !read statement
        do while(notDone)
            read(51,*,ioStat=ioErr)tKelvin,expVal
            if(ioErr < 0)then   !this means we reached the end of the file
                exit            !break the loop and return to remaining execution 
            elseif(ioErr > 0)then  
                write(52,*)'Unexpected error (e.g. type mismatch?) reading input line=',line
                cycle
            endif
            line=line+1
            call PsatEar(tKelvin,pMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
!            if(ierCode)then
!                write(52,*)'Unexpected error from Psat calculation. iErrCode,tKelvin,line=',ierCode,tKelvin,line
!                cycle
!            endif
            pKPa=pMPa*1000
            if(ierCode.ne.0)then
                pKPa=errDummy
                if(ierCode > 10)rhoLiq=errDummy
            endif
            if(iProperty==1)write(52,'( 3(1pE11.4,1x),i4)')tKelvin,expVal,pKPa,ierCode
            if(iProperty==2)write(52,'( 3(1pE11.4,1x),i4)')tKelvin,expVal,1000*rhoLiq,ierCode
        enddo
        GOTO 86
    endif !iProperty <= 2.	###########################################################################################################3333
    notDone=1
    line=0  !read statement	debugging
	if( ABS(Zc(1)) > 1E-11)rhoCritG_cc = Pc(1)*rMw(1)/( Zc(1)*rGas*Tc(1) )
	if(LOUD)print*,'PGLWRapper:MW,rhoCrit',rMw(1),rhoCritG_cc
	if(iProperty < 11)then ! calculate density, departure functions, and derivative properties given T,P
		do while(notDone)
			read(51,*,ioStat=ioErr)tKelvin,pKPa,expRho,iPhase
			if(ioErr < 0)then   !this means we reached the end of the file
				exit            !break the loop and return to remaining execution (notDone ignored really)
			elseif(ioErr > 0)then  
				write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
				cycle
			endif
			line=line+1
			!iPhase=1 !set as default
			!if(expRho/1000 < rhoCritG_cc)iPhase=0
			pMPa=pKPa/1000
			if(LOUD)write(*,'(a,f7.2,1x,3f9.5,i4)')'PGLWRapper: calling fugi. T,P,rhoc,expRho,iPhase:',tKelvin,pMPa,rhoCritG_cc,expRho/1000,iPhase
			call FUGI(tKelvin,pMPa,xFrac,NC,iPhase,FUGC,zFactor,ier)
	!        if(ier(1) .or. zFactor <= 0)then
	!            write(52,*)'Unexpected error from Psat calculation. iErrCode,tKelvin,line=',ierCode,tKelvin,line
	!            cycle
	!        endif
			ierCode=0
			if(ier(1) > 10)ierCode=sum(ier)*10+ier(1)
			if(LOUD.and.ierCode.ne.0)print*,'PGLWRapper: ier()=',(ier(i),i=1,11)
			if(ierCode > 10)then
				vpt=errDummy
				rhoG_cc=errDummy
			elseif(zFactor==0) then
				rhoG_cc=errDummy
				ierCode=ierCode*100+1
			else ! All good!
				rhoG_cc=pMPa*rMw(1)/(zFactor*rGas*tKelvin)
			endif
			rhoKg_m3=1000*rhoG_cc
			if(iProperty > 2 .and. iProperty < 6)write(52,'( f8.3,1x,3(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,1000*rhoG_cc,ierCode
			if(iProperty > 2 .and. iProperty < 6 .and.LOUD)write(*,*) ' tKelvin,pMPa,zFactor,expRho,rhoKg_m3,ierCode '
			if(iProperty > 2 .and. iProperty < 6 .and.LOUD)write(*,'( f8.3,1x,4(1pE11.4,1x),i4)')tKelvin,pMPa,zFactor,expRho,rhoKg_m3,ierCode
			if(cmprsblty==0)then
				vpt=errDummy
 				ierCode=ierCode*100+2
			else !cmprsblty good!
				vpt = 1/(cmprsblty*pKPa/zFactor)  !cmprsblty=(dP/dRho)T*(1/RT) => vpt = Z/(P/RT*dP/dRho) = 1/(rho*dP/dRho) = (1/rho)*(dRho/dP)_T
			endif
			! Res props are members of GlobConst, including CpRes and CvRes.
			if(iProperty > 5.and.LOUD)write(*,*) ' tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode '
			if(iProperty > 5.and.LOUD)write(*,'( f8.3,1x,6(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode
			if(iProperty > 5)write(52,'( f8.3,1x,6(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode
			if(LOUD)pause 'check output.'
		enddo !while(notDone)
	endif ! 2 < iProperty < 11 next for mixtures	#########################################################################################################
	tPrevious=0
	pPrevious=0
	xPure1=0
	xPure2=0
	xPure1(1)=1
	xPure2(2)=1
	! above we read: if(iProperty > 10)read(readString,*)iEosOpt,iProperty,localCas(1),localCas(2) !,xFrac(1)
	if(iProperty == 11)then	! compute Vex, Hex, Gex, and lnGamma's	for binary system
		write(52,'(a)')'      T(K)       p(kPa)         x1   phas Vex(cc/mol) rho(mol/cc)  Hex(J/mol)  Gex(J/mol)   ln(gam1)   ln(gam2) '
		do while(notDone)
			read(51,*,ioStat=ioErr)tKelvin,pKPa,iPhase,xFrac(1)
			aPhase='VAP'
			if(iPhase==1)aPhase='LIQ'
			xFrac(2) = 1-xFrac(1)
			if(LOUD .and. xFrac(2) < 0)pause 'PGLWRapper: Error xFrac(1) > 1  ??? '
			if(ioErr < 0)then   !this means we reached the end of the file
				exit            !break the loop and return to remaining execution 
			elseif(ioErr > 0)then  
				write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
				cycle
			endif
			line=line+1
			iErr=0
			pMPa=pKPa/1000
			if( ABS(tKelvin-tPrevious) > 0.1 .or. ABS(pKPa-pPrevious) > 0.1 )then
				if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure1(1)
				if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure2(1)
				call FUGI(tKelvin,pMPa,xPure2,NC,iPhase,FUGC,zFactor,ier)
				if(LOUD)print*,'PGLWRapperMain: Check output from fugi()pure2 call. ier(1) = ',ier(1)
				if(ier(1) > 10)iErr=iErr+10*(ier(1)-10)
				vPure2=rGas*tKelvin*zFactor/pMPa
				chemPoPure2=FUGC(2) ! = ln(fPure2/P)
				hRes2 = hRes_RT*rGas*tKelvin ! from module
				NC1 = 1 ! for pure compound1, we can ignore compound2. This may reduce error indications e.g. when T < Tmin for compound2 but OK for compound1. Too much trouble for compound2.
				call FUGI(tKelvin,pMPa,xPure1,NC,iPhase,FUGC,zFactor,ier)
				if(LOUD)print*,'PGLWRapperMain: Check output from fugi()pure1 call. ier(1) = ',ier(1)
				if(ier(1) > 10)iErr=ier(1)-10
				vPure1=rGas*tKelvin*zFactor/pMPa
				chemPoPure1=FUGC(1) ! = ln(fPure1/P)
				hRes1 = hRes_RT*rGas*tKelvin ! from module
				tPrevious=tKelvin
				pPrevious=pKPa
			endif
			!iPhase=1 !set as default
			!if(expRho/1000 < rhoCritG_cc)iPhase=0
            if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xFrac(1)
			call FUGI(tKelvin,pMPa,xFrac,NC,iPhase,FUGC,zFactor,ier)
            if(LOUD)print*,'PGLWRapperMain: Check output from fugi()mix call. ier(1) = ',ier(1)
			if(ier(1) > 10)iErr=iErr+100*(ier(1)-10)
			!iWarn=0
			!if(iErr==5 .or. iErr==500 .or. iErr==50 .or. iErr==55 .or. iErr==550 .or. iErr==505 .or. iErr==555)iWarn=1
			if( iErr > 10 )then
				write(52,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				if(LOUD)write(*,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				cycle
			endif
																														 
			vMix=rGas*tKelvin*zFactor/pMPa
			rhoMol_cc=1 / vMix
			activity1=FUGC(1)-chemPoPure1 ! = ln(fi/xi*fiPure) 
			activity2=FUGC(2)-chemPoPure2 ! = ln(fi/xi*fiPure)
			hResMix = hRes_RT*rGas*tKelvin ! from module
			vXsCc_mol = vMix - (xFrac(1)*vPure1	+ xFrac(2)*vPure2)
			hXsJ_mol = hResMix - (xFrac(1)*hRes1 + xFrac(2)*hRes2)
			gXsJ_mol =  (xFrac(1)*activity1 + xFrac(2)*activity2)*rGas*tKelvin
			write(52,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,vXsCc_mol,rhoMol_cc,hXsJ_mol,gXsJ_mol,activity1,activity2,iErr
		enddo !	while(notDone)
	elseif(iProperty == 12)then ! end iProperty=11, start iProperty = 12: rhoMol_cc,zFactor,hRes_RT
	!if(iProperty == 11)then	! compute Vex, Hex, Gex, and lnGamma's
		write(52,'(a)')'     T(K)    rho(mol/cc)      x1         pKPa         Z        Hres/RT        iErr  '
		do while(notDone)
			read(51,*,ioStat=ioErr)tKelvin,rhoMol_cc,xFrac(1)
			xFrac(2) = 1-xFrac(1)
			if(LOUD .and. xFrac(2) < 0)pause 'PGLWRapper: Error xFrac(1) > 1  ??? '
			if(ioErr < 0)then   !this means we reached the end of the file
				exit            !break the loop and return to remaining execution 
			elseif(ioErr > 0)then  
				write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
				cycle
			endif
			line=line+1
			iErr=0
			NC=2
			isZiter=0
			call FuVtot(isZiter,tKelvin,1/rhoMol_cc,xFrac,NC,FUGC,zFactor,iErr)
			if(ier(1) > 10)iErr=iErr+100*(ier(1)-10)
			!iWarn=0
			!if(iErr==5 .or. iErr==500 .or. iErr==50 .or. iErr==55 .or. iErr==550 .or. iErr==505 .or. iErr==555)iWarn=1
			if( iErr > 10 )then
				write(52,'(1x,3(E11.4,1x),1x,3E12.4,i8)')tKelvin,rhoMol_cc,xFrac(1),errDummy,errDummy,errDummy,iErr
				if(LOUD)write(*,'(1x,3(E11.4,1x),1x,3E12.4,i8)')tKelvin,pKPa,xFrac(1),errDummy,errDummy,errDummy,iErr
				cycle
			endif
			pKPa = zFactor*rhoMol_cc*rGas*tKelvin *1000
			write(52,'(1x,3(1PE11.4,1x),3(1PE11.4,1x),i8)')tKelvin,pKPa,xFrac(1),rhoMol_cc,zFactor,hRes_RT,iErr
		enddo ! 
	elseif(iProperty == 13)then	! compute Vex, Hex, Gex, and lnGamma's	for ternary system
		write(52,'(a)')'      T(K)       p(kPa)         x1   phas Vex(cc/mol) rho(mol/cc)  Hex(J/mol)  Gex(J/mol)   ln(gam1)   ln(gam2) '
		do while(notDone)
			read(51,*,ioStat=ioErr)xFrac(1),xFrac(2),tKelvin,pKPa
			aPhase='VAP'
			if(iPhase==1)aPhase='LIQ'
			xFrac(3) = 1-xFrac(1)-xFrac(2)
			if(xFrac(3) < 0)ioErr=123 !Use existing error traps in this case.
			if(LOUD .and. xFrac(2) < 0)pause 'PGLWRapper: Error xFrac(1) > 1  ??? '
			if(ioErr < 0)then   !this means we reached the end of the file
				exit            !break the loop and return to remaining execution 
			elseif(ioErr > 0)then  
				write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
				cycle
			endif
			line=line+1
			iErr=0
			pMPa=pKPa/1000
			if( ABS(tKelvin-tPrevious) > 0.1 .or. ABS(pKPa-pPrevious) > 0.1 )then
				if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure1(1)
				if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure2(1)
				call FUGI(tKelvin,pMPa,xPure2,NC,iPhase,FUGC,zFactor,ier)
				if(LOUD)print*,'PGLWRapperMain: Check output from fugi()pure2 call. ier(1) = ',ier(1)
				if(ier(1) > 10)iErr=iErr+10*(ier(1)-10)
				vPure2=rGas*tKelvin*zFactor/pMPa
				chemPoPure2=FUGC(2) ! = ln(fPure2/P)
				hRes2 = hRes_RT*rGas*tKelvin ! from module
				NC1 = 1 ! for pure compound1, we can ignore compound2. This may reduce error indications e.g. when T < Tmin for compound2 but OK for compound1. Too much trouble for compound2.
				call FUGI(tKelvin,pMPa,xPure1,NC,iPhase,FUGC,zFactor,ier)
				if(LOUD)print*,'PGLWRapperMain: Check output from fugi()pure1 call. ier(1) = ',ier(1)
				if(ier(1) > 10)iErr=ier(1)-10
				vPure1=rGas*tKelvin*zFactor/pMPa
				chemPoPure1=FUGC(1) ! = ln(fPure1/P)
				hRes1 = hRes_RT*rGas*tKelvin ! from module
				tPrevious=tKelvin
				pPrevious=pKPa
			endif
			!iPhase=1 !set as default
			!if(expRho/1000 < rhoCritG_cc)iPhase=0
            if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' PGLWRapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xFrac(1)
			call FUGI(tKelvin,pMPa,xFrac,NC,iPhase,FUGC,zFactor,ier)
            if(LOUD)print*,'PGLWRapperMain: Check output from fugi()mix call. ier(1) = ',ier(1)
			if(ier(1) > 10)iErr=iErr+100*(ier(1)-10)
			!iWarn=0
			!if(iErr==5 .or. iErr==500 .or. iErr==50 .or. iErr==55 .or. iErr==550 .or. iErr==505 .or. iErr==555)iWarn=1
			if( iErr > 10 )then
				write(52,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				if(LOUD)write(*,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				cycle
			endif
																														 
			vMix=rGas*tKelvin*zFactor/pMPa
			rhoMol_cc=1 / vMix
			activity1=FUGC(1)-chemPoPure1 ! = ln(fi/xi*fiPure) 
			activity2=FUGC(2)-chemPoPure2 ! = ln(fi/xi*fiPure)
			hResMix = hRes_RT*rGas*tKelvin ! from module
			vXsCc_mol = vMix - (xFrac(1)*vPure1	+ xFrac(2)*vPure2)
			hXsJ_mol = hResMix - (xFrac(1)*hRes1 + xFrac(2)*hRes2)
			gXsJ_mol =  (xFrac(1)*activity1 + xFrac(2)*activity2)*rGas*tKelvin
			write(52,'(1x,3(E11.4,1x),a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,vXsCc_mol,rhoMol_cc,hXsJ_mol,gXsJ_mol,activity1,activity2,iErr
		enddo !	while(notDone)
	endif ! 13 < iProperty 	 (ie. mixtures) ############################################################################################
601	FORMAT(I3,I5,9(1X,F7.4))
86	continue
    if(LOUD)write(*,*)'Success! Results are written to output.txt'
	close(52)
	close(51)
	stop
	END

