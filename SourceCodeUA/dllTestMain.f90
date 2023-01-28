	!MAIN PROGRAM FOR CALLING EOS SUBROUTINES. Includes GetCrit, GetBips, BipIo.
	!	SETS UP THE CRITS, EOSPARMS, bips, DEBUG status, FUGI(iEosOpt)
	!	Echoes user IO to Output.txt, and reports error checks. 
	!	INITIALIZATION AND CALLING SEQUENCE FOR VLE, LLE, VLLE SUBROUTINES.
	!reqd routines:
	!	bubpl.for, bubtl.for, dewtv.for, flashsub.for, FuEsdMy.for FuEsdXs2.for, FugiPr.f90, FugiPrws.for, RegPure.f90
	!	LmDifEzCov2.for, Mintools.for
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	!USE PORTLIB
	!USE GlobConst
	!USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	!CHARACTER*1 ANSWER
	!CHARACTER*2 calcType
	CHARACTER*77 errMsgPas
	CHARACTER*251 outFile

	!DIMENSION BIPTMP(NMX)
	!DIMENSION ZFEED(NMX),S(NMX)
    DoublePrecision xFrac(NMX),FUGC(NMX) !FUGI requires mole fraction specification because it is written generally for mixtures.
    INTEGER ier(12),localCas(NMX)
    DoublePrecision CalcResult
    integer ieos, casrn, prp_id, ierr
    doublePrecision var1, var2
    character*255 tag, value
	
	COMMON/eta/etaL,etaV,ZL,ZV

	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	!iChange=VERIFY('C:\MYPROJEX\CalcEos',CURDIR)
	!iChange2=VERIFY('C:\MYPROJEX\CALCEOS',CURDIR)
	!if(iChange2.eq.0)iChange=0
	!if(iChange.eq.0 .or. iChange.gt.26)DEBUG=.TRUE.
	masterDir=TRIM(curDir)
    tag='LOCATION'
    value='C:\Soft\myPGLwrapper\PGLwrapper'
    masterDir=trim(value)
    PGLinputDir='C:\Soft\myPGLwrapper\PGLwrapper\Input'
    !ierr=SETSTRING(tag, value)
	call LoadCritParmsDb(iErrCrit)
	if(iErrCrit > 0 .and. LOUD)print*,'UaWrapperMain: From LoadCritParmsDb, iErrCrit=',iErrCrit
    call Test1()
    stop
    call Test2()
    stop
    call Test3()
    stop
	outFile=TRIM(masterDir)//'\output.txt' ! // is the concatenation operator
	OPEN(UNIT=52,FILE=outFile,RECL=8192)
    open(UNIT=51,FILE='input.txt')
    NC=1 !assume one component for all calculations (as of 1/1/2020 this is all we need).
    xFrac(1)=1  !   "
    read(51,*)iEosOpt,iProperty,localCas(1)
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (g/cc) given tKelvin
    !iProperty = 3: vapor fluid density (g/cc) given tKelvin, pKPa, iPhase (=0 for vapor)
    !iProperty = 4: liquid fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid)
    !iProperty = 5: critical fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 6: sound speed
    !iProperty = 7: Hres/RT,CpRes_R/R,CvResR,cmprsblty given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)

    !iEosOpt=10 !Temporary for testing

	INITIAL=0
	DEBUG=.TRUE.
	DEBUG=.FALSE.
    LOUD = .TRUE.
    LOUD = .FALSE.
	ierCode=8686 ! Initialize to failure in the iEosOpt write slot to make it easier to terminate if any Get_ functions fail.
    
    ieos=5
    casrn=67641
    prp_id=8
    var1=300.0
    var2=0.0
    ierr=0
    CalcResult = CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
    CalcResult = CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
    stop
    
	call IdDipprLookup(NC,localCas,iErrCas,errMsgPas)
	if(LOUD)print*,'localCas,idDippr=',localCas(1),id(1)
	if(iErrCas)then
        if(LOUD)write(6 ,*)' Sorry, must abort.  Not found for CAS number=', localCas(1)
        write(52,*)ierCode,' Sorry, must abort.  Not found for CAS number=', localCas(1)
        goto 86
    endif
	if(LOUD)print*,'idDippr,class=',ID(1),TRIM(class(1))
	CALL GETCRIT(NC,iErrCrit)
	if(iErrCrit)then
		if(LOUD)write(5 ,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		write(52,*)ierCode,' Error in Main: ErrorCode from GetCrit = ',iErrCrit
		goto 86
	endif
	if(LOUD)print*,'calling GetEOS'
	iErrGet=0 
	if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	if(iEosOpt.eq.2)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParms					!Diky model 12
	if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet) 
	if(iEosOpt.eq.4)CALL GetEsdCas(NC,localCas,iErrGet)
	if(iEosOpt.eq.5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms					!Diky model 6
!	if(iEosOpt.eq.6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	if(iEosOpt.eq.8)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.9)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.10)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters			!Diky model 21
	if(iEosOpt.eq.11)CALL GetPrTc(NC,iErrGet)		!JRE 2019 : Reads Jaubert's parameters						  !Diky model 22
	if(iEosOpt.eq.12)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmami					 !Diky model 23
	if(iEosOpt.eq.13)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmamiTb					 !Diky model 24
	if(iEosOpt.eq.14)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms				    !Diky model 18
	!NewEos: Add here for initializing parms.
	if(iErrGet .ne. 0)then
		if(LOUD)write(* ,*)ierCode,' Error in Main: failed to get required dbase props.'
		write(52,*)ierCode,' Error in Main: failed to get required dbase props.'
		goto 86
	endif
    if(LOUD)write(*,*)iEosOpt,' = iEosOpt'
	if(iProperty==1)then
		write(52,'(2i4,i11,a)')iEosOpt, 1, localCas(1),' = iEosOpt, iProp, localCas: T, P,expVal,pSatKPa,ierCode '
	elseif(iProperty > 1 .and. iProperty < 6)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,rho(g/cc),ierCode '
	elseif(iProperty>5)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode '
	else
		write(52,*)iEosOpt,' = iEosOpt'
	endif
	ierCode=0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
    if(iProperty==1 .or. iProperty==2)then
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
                pKPa=86.8686D0
                rhoLiq=86.8686D0
            endif
            if(iProperty==1)write(52,'( 3(1pE11.4,1x),i4)')tKelvin,expVal,pKPa,ierCode
            if(iProperty==2)write(52,'( 3(1pE11.4,1x),i4)')tKelvin,expVal,1000*rhoLiq,ierCode
        enddo
        GOTO 86
    endif !iProperty <= 2.
    notDone=1
    line=0  !read statement
	if( ABS(Zc(1)) > 1E-11)rhoCritG_cc = Pc(1)*rMw(1)/( Zc(1)*rGas*Tc(1) )
	if(LOUD)print*,'UaWrapper:MW,rhoCrit',rMw(1),rhoCritG_cc
    do while(notDone)
        read(51,*,ioStat=ioErr)tKelvin,pKPa,expRho,iPhase
        if(ioErr < 0)then   !this means we reached the end of the file
            exit            !break the loop and return to remaining execution 
        elseif(ioErr > 0)then  
            write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
            cycle
        endif
        line=line+1
		!iPhase=1 !set as default
		!if(expRho/1000 < rhoCritG_cc)iPhase=0
        pMPa=pKPa/1000
		if(LOUD)write(*,'(a,f7.2,1x,3f9.5,i4)')'UaWrapper: calling fugi. T,P,rhoc,expRho,iPhase:',tKelvin,pMPa,rhoCritG_cc,expRho/1000,iPhase
        call FUGI(tKelvin,pMPa,xFrac,NC,iPhase,FUGC,zFactor,iErrF)
!        if(ier(1) .or. zFactor <= 0)then
!            write(52,*)'Unexpected error from Psat calculation. iErrCode,tKelvin,line=',ierCode,tKelvin,line
!            cycle
!        endif
		ierCode=0
		if(ier(1).ne.0)ierCode=sum(ier)*10+ier(1)
		if(LOUD.and.ierCode.ne.0)print*,'UaWrapper: ier()=',(ier(i),i=1,11)
        if(ierCode.ne.0)then
            vpt=86.8686D0
            rhoG_cc=86.8686D0
        elseif(zFactor==0) then
            rhoG_cc=86.8686D0
			ierCode=ierCode*100+1
       else
            rhoG_cc=pMPa*rMw(1)/(zFactor*rGas*tKelvin)
        endif
		rhoKg_m3=1000*rhoG_cc
        if(iProperty > 2 .and. iProperty < 6)write(52,'( f8.3,1x,3(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,1000*rhoG_cc,ierCode
        if(iProperty > 2 .and. iProperty < 6 .and.LOUD)write(*,*) ' tKelvin,pMPa,zFactor,expRho,rhoKg_m3,ierCode '
        if(iProperty > 2 .and. iProperty < 6 .and.LOUD)write(*,'( f8.3,1x,4(1pE11.4,1x),i4)')tKelvin,pMPa,zFactor,expRho,rhoKg_m3,ierCode
        if(cmprsblty==0)then
            vpt=86.8686D0
 			ierCode=ierCode*100+2
		else
            vpt = 1/(cmprsblty*pKPa/zFactor)  !cmprsblty=(dP/dRho)T*(1/RT) => vpt = Z/(P/RT*dP/dRho) = 1/(rho*dP/dRho) = (1/rho)*(dRho/dP)_T
		endif
        if(iProperty > 5.and.LOUD)write(*,*) ' tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode '
        if(iProperty > 5.and.LOUD)write(*,'( f8.3,1x,6(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode
        if(iProperty > 5)write(52,'( f8.3,1x,6(1pE11.4,1x),i4)')tKelvin,pKPa,expRho,hRes_RT,CpRes_R,CvRes_R,vpt,ierCode
		if(LOUD)pause 'check output.'
    enddo
601	FORMAT(I3,I5,9(1X,F7.4))
86	continue
    if(LOUD)write(*,*)'Success! Results are written to output.txt'
	close(52)
	close(51)
	stop
    END

subroutine Test1
IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
integer ieos, Rn1, Rn2, Rn3, prp_id, ierr
doublePrecision var1, var2, prp, uncert
DoublePrecision CalcResult
    ieos=5
    Rn1=67641
    Rn2=0
    Rn3=0
    prp_id=1
    var1=300.0
    var2=0.0
    ierr=0
    ierr=INITIALIZE_MODEL(ieos, Rn1, Rn2, Rn3)
    call CalculateProperty1local(ieos, Rn1, prp_id, var1, var2, prp, ierr)
    CalcResult = CalculateProperty(ieos, Rn1, prp_id, var1, var2, ierr)
!    prp_id=8
!    ierr=Calculate1(Rn1, ieos, prp_id, var1, var2, prp, uncert)
    write(*,*) prp
	end subroutine

subroutine Test2
integer ieos, casrn1, casrn2, prp_id, ierr
doublePrecision var1, var2, var3, prp
    ieos = 5
    casrn1 = 75150
    casrn2 = 2551624
    var1 = 298.136
    var2 = 1000
    var3 = 0.0013779105351375952
    prp_id = 201
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    var2 = 91911.609110881604
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    var2 = 1000
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    return

    ieos = 2
    casrn1 = 123911
    casrn2 = 127184
    var1 = 298.15
    var2 = 101.3
    var3 = 0.085
    prp_id = 36
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    return

    ieos = 2
    casrn1 = 75525
    casrn2 = 75058
    var1 = 298.136
    var2 = 101
    var3 = 0.9428
    prp_id = 89
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
	end subroutine

subroutine Test3	!testing dll calls
integer i, ieos, casrn1, casrn2, casrn3, prp_id, ierr
doublePrecision var1, var2, var3, prp, p
integer GETNPAR,GETPAR
	i=INITIALIZE_MODEL(2,71363,7732185,0)
	i=GETNPAR(2)
	i=GETPAR(1,p)
	return
    ieos = 5
    casrn1 = 75150
    casrn2 = 2551624
    var1 = 298.136
    var2 = 1000
    var3 = 0.0013779105351375952
    prp_id = 201
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    var2 = 91911.609110881604
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    var2 = 1000
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    return

    ieos = 2
    casrn1 = 123911
    casrn2 = 127184
    var1 = 298.15
    var2 = 101.3
    var3 = 0.085
    prp_id = 36
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
    return

    ieos = 2
    casrn1 = 75525
    casrn2 = 75058
    var1 = 298.136
    var2 = 101
    var3 = 0.9428
    prp_id = 89
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, prp, ierr)
    write(*,*) prp
end subroutine
	