MODULE DLLConst
	!IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!SAVE
	!PUBLIC
	Implicit NONE
	Integer oldRN1,oldRN2,oldRN3,oldEOS
	Character*15 EosName(18) 
	!               1     2       3       4          5          6         7           8              9            10          11      12        13         14          15           16            17        18
	data EosName/'PR','ESD96','PRWS','ESD-MEM2','SPEADMD','Flory-MEM2','NRTL','SpeadGamma-MEM2','SPEAD11','PcSaft(Gross)','tcPRq','GCESD','GCESD(Tb)','TransSPEAD','GcPcSaft','GcPcSaft(Tb)','tcPR-GE(W)','ESD2'/
END MODULE DLLConst

subroutine CalculateProperty1local(ieos, casrn, prp_id, var1, var2, res, ierr)
!	CalculateProperty1 & local is for pure compounds.
	!USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE GlobConst
    USE DllConst
	Implicit NONE
    integer ieos, casrn, prp_id, ierr
    double Precision var1, var2, res
    double Precision xFrac(nmx),FUGC(nmx) !FUGI requires mole fraction specification because it is written generally for mixtures.
    INTEGER localCas(nmx)
    integer NC,IPROPERTY,IPHASE,INITIAL,iErrStart,line,ierCode,iErrF,iErrDerv,isZiter,notDone
    double Precision tKelvin,vTotCc,Z,aRes,uRes,rhoLiq,rhoVap,uSatL,uSatV,chemPot,pKPa,pMPa,zFactor,aRes_RT,uRes_RT,rhoMol_cc,rhoG_cc,hRes_RT
	if (LOUD.and. dumpUnit==6)return
    if(LOUD)write(dumpUnit,*)'CalculateProperty1local: ieos,casrn,prp_id=',ieos,casrn 
    if(LOUD)write(dumpUnit,610)'CalculateProperty1local: var1,var2=',var1,var2
610 format(1x,a,12E12.4)

    NC=1 !assume one component for all calculations (as of 1/1/2020 this is all we need).
    xFrac(1)=1  !   "
    if(iEosOpt/=ieos)then
        if(LOUD)write(dumpUnit,611)'CalculateProperty1local: failed ieos check. iEosOpt=',iEosOpt
        iErr=15
        goto 86
    endif
    iProperty=ABS(prp_id)
    if (prp_id==-2) iProperty=3
    localCas(1)=casrn
    !write(*,*)casrn
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (g/cc) given tKelvin
    !iProperty = 3: fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 4: Hres/RT(41),CpRes/R(42),CvResR(43),cmprsblty(44) given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)
    iPhase=1
    if (prp_id==3) iPhase=0   !vapor or gas
    res=0
    ierr=0
	
	INITIAL=0
    if(casrn.ne.oldRN1.or.iEosOpt.ne.oldEOS)then
		Call PGLWrapperStartup(NC,iEosOpt,localCas,iErrStart) !LoadCritParmsDb,IdDipprLookup,GetCrit
        if(LOUD)write(dumpUnit,*)'CalculateProperty1local: After Start, iErr,ID=',iErrStart,ID(1:NC) 
	    if(iErrStart/=0)then
            ierr=1
            goto 86
        endif
        oldEOS=iEosOpt
        oldRN1=casrn
    endif
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
    if (iProperty==34)then
        tKelvin=var2
        vTotCc = 1000*rMw(1)/var1
        isZiter=0
        call FuVtot(isZiter,tKelvin,vTotCc,xFrac,NC,FUGC,Z,aRes,uRes,iErr)
        res=Z*rGas*tKelvin*var1/rMw(1)
        goto 86
    endif

    if(iProperty==1 .or. iProperty==2 .or. iProperty==13)then
        notDone=1
        line=0  !read statement
!        do while(notDone)
            tKelvin=var1
            line=line+1
            call PsatEar(tKelvin,pMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
            if(LOUD)write(dumpUnit,611)'CalculateProperty1local: After Psat, iErr,P,rhoL,rhoV=',ierCode,pMPa,rhoLiq,rhoVap 
            if(ierCode.NE.0)then
!                write(52,*)'Unexpected error from Psat calculation. iErrCode,tKelvin,line=',ierCode,tKelvin,line
                ierr=3
                goto 86
            endif
            pKPa=pMPa*1000
            if(iProperty==1) res=pKPa
            if(iProperty==2) res=1000*rhoLiq*rMw(1)
            if(iProperty==13) res=1000*rhoVap*rMw(1)
            goto 86
!        enddo
    endif !iProperty <= 2.
611 format(1x,a,i11,12E12.4)
    notDone=1
    line=0  !read statement
!    do while(notDone)
        tKelvin=var1
        pKPa=var2
        line=line+1
        pMPa=pKPa/1000
        !call FUGI(tKelvin,pMPa,xFrac,NC,iPhase,FUGC,zFactor,iErrF)
		CALL FugiTP( tKelvin,pMPa,xFrac,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
        if (iErrF.ne.0) then
            iPhase=1-iPhase
            CALL FugiTP( tKelvin,pMPa,xFrac,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
        endif
        if(LOUD)write(dumpUnit,611)'CalculateProperty1local: After FugiTP, iErr,Z,rho(g/cc)=',iErrF,zFactor,rhoMol_cc*rMw(1) 
        if(iErrF > 0 .or. zFactor <= 0)then
!            write(52,*)'Unexpected error from Psat calculation. iErrCode,tKelvin,line=',ierCode,tKelvin,line
            iErr=4
            goto 86
        endif
        rhoG_cc=rhoMol_cc*rMw(1)
		hRes_RT=uRes_RT+zFactor-1
		iErrDerv=0
		if(iProperty > 41)then
            CALL NUMDERVS(NC,xFrac,tKelvin,rhoMol_cc,zFactor,iErrDerv) !cmprsblty,CpRes_R,CvRes_R
            if(LOUD)write(dumpUnit,611)'CalculateProperty1local: After NUMDERVS, iErr,...=',iErrDerv,cmprsblty,CpRes_R,CvRes_R 
        endif
		if(iErrDerv > 0)then
			iErr=5
			goto 86
		endif
        if (iProperty==3) res=1000*rhoG_cc
        if (iProperty==4) res=1000*rhoG_cc
        !if(iProperty==4)write(52,*)tKelvin,pKPa,hRes_RT,CpRes_R,CvRes_R,cmprsblty  !cmprsblty=(dP/dRho)T*(1/RT)
        if (iProperty==41) res=hRes_RT
        if (iProperty==42) res=CpRes_R
        if (iProperty==43) res=CvRes_R
        if (iProperty==44) res=1/(cmprsblty*pKPa/zFactor)  !cmprsblty=(dP/dRho)T*(1/RT)
        if (iProperty==45) res=CpRes_R
!    enddo
86  if(LOUD)write(dumpUnit,611)'CalculateProperty1local: returning, iErr,result=',iErr,res 
    return
end subroutine CalculateProperty1local

function CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE GlobConst
	IMPLICIT double Precision(A-H,K,O-Z)
    double Precision CalculateProperty
    integer ieos, casrn, prp_id, ierr
    double Precision var1, var2
	!DEC$ ATTRIBUTES DLLEXPORT::CalculateProperty
    !!MS$ ATTRIBUTES DLLEXPORT::CalculateProperty
    call CalculateProperty1local(ieos, casrn, prp_id, var1, var2, res, ierr)
    CalculateProperty=res
    return
end function CalculateProperty

subroutine CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, res, ierr)
!	CalculateProperty2 & local is for binary mixtures.
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE GlobConst
    use DllConst
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    double Precision res
    integer ieos, casrn1, casrn2, prp_id, ierr
    double Precision var1, var2, var3
    double Precision xFrac(nmx),FUGC(nmx),xPure1(nmx),xPure2(nmx) !FUGI requires mole fraction specification because it is written generally for mixtures.
    INTEGER localCas(nmx)
	!CHARACTER*77 errMsgPas
!	COMMON/eta/etaL,etaV,ZL,ZV
	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
    IF (LOUD.and.dumpUnit==6) return
	DEBUG=.FALSE.
	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	!iChange=VERIFY('C:\MYPROJEX\CalcEos',CURDIR)
	!iChange2=VERIFY('C:\MYPROJEX\CALCEOS',CURDIR)
	!if(iChange2.eq.0)iChange=0
	!if(iChange.eq.0 .or. iChange.gt.26)DEBUG=.TRUE.
	NC=2 !no of components
    iEosOpt=ieos
    iProperty=ABS(prp_id)
    localCas(1)=casrn1
    localCas(2)=casrn2
    !write(*,*)casrn
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (g/cc) given tKelvin
    !iProperty = 3: fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 4: Hres/RT,CpRes/R,CvResR,cmprsblty given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)
    iPhase=1
    if (prp_id==201) iPhase=0   !gas
    if (prp_id==202) iPhase=0   !gas
    res=0
    ierr=0
	
	INITIAL=0

    if (oldRN1.ne.casrn1 .or. oldRN2.ne.casrn2 .or. oldEOS.ne.ieos) then
		Call PGLWrapperStartup(NC,iEosOpt,localCas,iErrStart) !LoadCritParmsDb,IdDipprLookup,GetCrit
	    if(iErrStart.NE.0)then
		    ierr=3
		    return
        endif
        oldRN1=casrn1
        oldRN2=casrn2
        oldEOS=ieos
	endif
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
        tPrevious = 0
        pPrevious = 0
	    xPure1(2)=0
	    xPure2(1)=0
	    xPure1(1)=1
	    xPure2(2)=1
		!write(52,'(a)')'   T(K)   p(kPa)     x1  phas Vex(cc/mol) rho(mol/cc)  Hex(J/mol)  Gex(J/mol)   ln(gam1)   ln(gam2) '
		!do while(notDone)
			!read(51,*,ioStat=ioErr)tKelvin,pKPa,iPhase,xFrac(1)
        tKelvin = var1
        pKPa = var2
        xFrac(1) = var3
			aPhase='VAP'
			if(iPhase==0)aPhase='GAS'
			if(iPhase==1)aPhase='LIQ'
			xFrac(2) = 1-xFrac(1)
			!if(xFrac(2) < 0)pause 'UaWrapper: Error xFrac(1) > 1  ??? '
			!if(ioErr < 0)then   !this means we reached the end of the file
			!	exit            !break the loop and return to remaining execution 
			!elseif(ioErr > 0)then  
			!	write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
			!	cycle
			!endif
			!line=line+1
			iErr=0
			pMPa=pKPa/1000
			if( ABS(tKelvin-tPrevious) > 0.1 .or. ABS(pKPa-pPrevious) > 0.1 )then
				CALL FugiTP( tKelvin,pMPa,xPure1,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
				if(iErrF > 0)iErr=iErrF
				hRes_RT=uRes_RT+zFactor-1
				vPure1=rGas*tKelvin*zFactor/pMPa
				chemPoPure1=FUGC(1) ! = ln(fPure1/P)
				hRes1 = hRes_RT*rGas*tKelvin ! from module
				CALL FugiTP( tKelvin,pMPa,xPure2,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
				if(iErrF > 0)iErr=iErrF
				hRes_RT=uRes_RT+zFactor-1
				vPure2=rGas*tKelvin*zFactor/pMPa
				chemPoPure2=FUGC(2) ! = ln(fPure2/P)
				hRes2 = hRes_RT*rGas*tKelvin ! from module
				tPrevious=tKelvin
				pPrevious=pKPa
			endif
			!iPhase=1 !set as default
			!if(expRho/1000 < rhoCritG_cc)iPhase=0
				!if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' UaWrapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xFrac(1)
			CALL FugiTP( tKelvin,pMPa,xFrac,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
			hRes_RT=uRes_RT+zFactor-1
				!if(LOUD)print*,'UaWrapperMain: Check output from fugi()mix call. ier(1) = ',ier(1)
			if(iErrF > 0)iErr=iErr+100*iErrF
			if(iErr > 0)then
				!write(52,'(1x,2f8.3,F7.4,1x,a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				!if(LOUD)write(*,'(1x,2f8.3,F7.4,1x,a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				!cycle
                res = 0
                return
			endif
																														 
			vMix=rGas*tKelvin*zFactor/pMPa
			rhoMol_cc=1/vMix
			activity1=FUGC(1)-chemPoPure1 ! = ln(fi/xi*fiPure) 
			activity2=FUGC(2)-chemPoPure2 ! = ln(fi/xi*fiPure)
			hResMix = hRes_RT*rGas*tKelvin ! from module
			vXsCc_mol = vMix - (xFrac(1)*vPure1	+ xFrac(2)*vPure2)
			hXsJ_mol = hResMix - (xFrac(1)*hRes1 + xFrac(2)*hRes2)
			gXsJ_mol =  (xFrac(1)*activity1 + xFrac(2)*activity2)
			!write(52,'(1x,2f8.3,F7.4,1x,a3,1x,7E12.4)')tKelvin,pKPa,xFrac(1),aPhase,vXsCc_mol,rhoMol_cc,hXsJ_mol,gXsJ_mol,activity1,activity2
            if (iProperty==4) res = exp(activity1)
            if (iProperty==5) res = exp(activity2)
            if (iProperty==81) res = 1000*rhoMol_cc
            if (iProperty==89) res = 0.000001*vXsCc_mol
            if (iProperty==36) res = 0.001*hXsJ_mol
            if (iProperty==201) res = exp(FUGC(1))
            if (iProperty==202) res = exp(FUGC(2))
            if (iProperty==203) res = exp(FUGC(1))
            if (iProperty==204) res = exp(FUGC(2))
		!enddo !	while(notDone)
        return
end subroutine CalculateProperty2local

function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE GlobConst
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    DoublePrecision CalculateProperty2
    integer ieos, casrn1, casrn2, prp_id, ierr
    DoublePrecision var1, var2, var3, res ! t,p,x, result
    !DoublePrecision xFrac(nmx),FUGC(nmx),xPure1(nmx),xPure2(nmx) !FUGI requires mole fraction specification because it is written generally for mixtures.
!	COMMON/eta/etaL,etaV,ZL,ZV
  !DEC$ATTRIBUTES DLLEXPORT::CalculateProperty2
  !!MS$ATTRIBUTES DLLEXPORT::CalculateProperty2
    call CalculateProperty2local(ieos, casrn1, casrn2, prp_id, var1, var2, var3, res, ierr)
    CalculateProperty2=res
    return
end function CalculateProperty2

integer function Activate(hello)
    use DllConst
    use GlobConst
    character(255) hello,local
    !DEC$ATTRIBUTES DLLEXPORT::Activate
    LOUD=.FALSE.
    !LOUD=.TRUE.
    if (LOUD)then
        dumpUnit=686
        open(dumpUnit,file='c:\myprojex\pgldll\pgldll\DebugDLL.txt')
        write(dumpUnit,*)'Activate: DLL has started.'
    endif
    oldRN1=0
    oldRN2=0
    oldRN3=0
    oldEOS=0
    Activate=1
	local=TRIM(hello)
    
    return
end function Activate

integer function QUERYMODEL(no, model_type, level, modelname)
    use GlobConst
    integer no, model_type, level
    character(255) modelname
    !DEC$ATTRIBUTES DLLEXPORT::QUERYMODEL
    QUERYMODEL=0
    if (no.eq.1) then
        model_type=1
        level=2
        QUERYMODEL=1
    endif
    if (no.eq.2) then
        model_type=1
        level=2
        QUERYMODEL=2
    endif
    if (no.eq.3) then
        model_type=1
        level=2
        QUERYMODEL=3
    endif
    if (no.eq.4) then
        model_type=1
        level=2
        QUERYMODEL=4
    endif
    if (no.eq.5) then
        model_type=1
        level=2
        QUERYMODEL=5
    endif
    if (no.eq.6) then
        model_type=1
        level=2
        QUERYMODEL=6
    endif
    if (no.eq.7) then
        model_type=2
        level=2
        QUERYMODEL=7
    endif
    if (no.eq.8) then
        model_type=1
        level=2
        QUERYMODEL=8
    endif
    if (no.eq.9) then
        model_type=1
        level=2
        QUERYMODEL=9
    endif
    if (no.eq.10) then
        model_type=1
        level=2
        QUERYMODEL=10
    endif
    if (no.eq.11) then
        model_type=1
        level=2
        QUERYMODEL=11
    endif
    if (no.eq.12) then
        model_type=1
        level=2
        QUERYMODEL=12
    endif
    if (no.eq.13) then
        model_type=1
        level=2
        QUERYMODEL=13
    endif
    if (no.eq.14) then
        model_type=1
        level=2
        QUERYMODEL=14
    endif
    if (no.eq.15) then
        model_type=1
        level=2
        QUERYMODEL=15
    endif
    if (no.eq.16) then
        model_type=1
        level=2
        QUERYMODEL=16
    endif
    if (no.eq.17) then
        model_type=1
        level=2
        QUERYMODEL=3
    endif
    if (no.eq.18) then
        model_type=3
        level=2
        QUERYMODEL=2
    endif
    if (no.eq.19) then
        model_type=3
        level=2
        QUERYMODEL=3
    endif
    if (no.eq.20) then
        model_type=3
        level=2
        QUERYMODEL=4
    endif
    if (no.eq.21) then
        model_type=3
        level=2
        QUERYMODEL=5
    endif
    if (no.eq.22) then
        model_type=3
        level=2
        QUERYMODEL=6
    endif
    if (no.eq.23) then
        model_type=3
        level=2
        QUERYMODEL=8
    endif
    if (no.eq.24) then
        model_type=3
        level=2
        QUERYMODEL=9
    endif
    if (no.eq.25) then
        model_type=3
        level=2
        QUERYMODEL=10
    endif
    if (no.eq.26) then
        model_type=3
        level=2
        QUERYMODEL=11
    endif
    if (no.eq.27) then
        model_type=3
        level=2
        QUERYMODEL=12
    endif
    if (no.eq.28) then
        model_type=3
        level=2
        QUERYMODEL=13
    endif
    if (no.eq.29) then
        model_type=3
        level=2
        QUERYMODEL=14
    endif
    if (no.eq.30) then
        model_type=3
        level=2
        QUERYMODEL=15
    endif
    if (no.eq.31) then
        model_type=3
        level=2
        QUERYMODEL=16
    endif
    if (QUERYMODEL>0) modelname=EosName(QUERYMODEL)
    return
    end function QUERYMODEL

integer function Calculate2(casrn1, casrn2, modelid, propertyid, t, p, x, res, uncert)
!int *cmpid1, int* cmpid2, int* modelid, int* propertyid, double* t, double* p, double* x, double* res, double* uncert
!double Precision function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    integer modelid, casrn1, casrn2, propertyid, ierr
    double Precision t, p, x, res, uncert
    !DEC$ ATTRIBUTES DLLEXPORT::Calculate2
    !!MS$ ATTRIBUTES DLLEXPORT::Calculate2
    if (x.le.0) x=0.000001
    if (x.ge.1) x=0.999999
    call CalculateProperty2local(modelid, casrn1, casrn2, propertyid, t, p, x, res, ierr)
    uncert=0
    Calculate2=ierr
    return
    end function Calculate2

integer function Calculate3(casrn1, casrn2, casrn3, modelid, propertyid, t, p, x1, x2, res, uncert)
!int *cmpid1, int* cmpid2, int* modelid, int* propertyid, double* t, double* p, double* x, double* res, double* uncert
!double Precision function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    integer modelid, casrn1, casrn2, casrn3, propertyid, ierr
    double Precision t, p, x1, x2, res, uncert
    !DEC$ ATTRIBUTES DLLEXPORT::Calculate3
    call CalculateProperty3local(modelid, casrn1, casrn2, casrn3, propertyid, t, p, x1, x2, res, ierr)
    uncert=0
    Calculate2=ierr
    uncert=0
    Calculate3=ierr
    return
end function Calculate3

integer function Calculate(casrn, modelid, propertyid, t, p, x, res, uncert)
!int *cmpid1, int* cmpid2, int* modelid, int* propertyid, double* t, double* p, double* x, double* res, double* uncert
!double Precision function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    integer modelid, casrn(255), propertyid,localCas !, ierr
    double Precision t, p, x(255), res, uncert
    !DEC$ ATTRIBUTES DLLEXPORT::Calculate
    !!MS$ ATTRIBUTES DLLEXPORT::Calculate
    res=0
    uncert=0
    Calculate=1
	ieos=modelid
	localCas=casrn(1)
	var1=t
	var2=p
	var3=x(1)
	prp_id=propertyid
    return
end function Calculate

integer function Calculate1(casrn1, modelid, propertyid, t, p, res, uncert)
!int *cmpid1, int* cmpid2, int* modelid, int* propertyid, double* t, double* p, double* x, double* res, double* uncert
!double Precision function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
	USE GlobConst
    Implicit NONE
    integer modelid, casrn1, propertyid, localprpid, ierr
    double Precision t, p, res, uncert
    !DEC$ ATTRIBUTES DLLEXPORT::Calculate1
    !!MS$ ATTRIBUTES DLLEXPORT::Calculate1
    localprpid = 0
    if (propertyid.eq.1) localprpid=4   !L
    if (propertyid.eq.2) localprpid=3   !G
    if (propertyid.eq.3) localprpid=2   !L+G
    if (propertyid.eq.4) localprpid=13   !G+L
    if (propertyid.eq.8) localprpid=1   !VP
    if (propertyid.eq.12) localprpid=42   !CP
    if (propertyid.eq.14) localprpid=45   !CP
    if (propertyid.eq.34) localprpid=34   !fluid pressure
    if (propertyid.eq.5) then
        res=Tc(1)
        Calculate1=0
    else if (propertyid.eq.6) then
        res=1000*Pc(1)
        Calculate1=0
    else if (propertyid.eq.7) then
        res=1000*Pc(1)*rMw(1)/(Zc(1)*Rgas*Tc(1))
        Calculate1=0
    else if (localprpid.eq.0) then
        res=0
        Calculate1=1
    else
        call CalculateProperty1local(modelid, casrn1, localprpid, t, p, res, ierr)
        if (propertyid.eq.12)then
            res=res*rGas
        endif
        uncert=0
        Calculate1=ierr
    end if
    return
end function Calculate1

integer function FIND_COMP(name)
    character(255) name
    !DEC$ ATTRIBUTES DLLEXPORT::FIND_COMP
    FIND_COMP=0
    return
end function FIND_COMP

integer function SUPPORTS_COMP(modelid, id1)
    integer modelid, id1
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_COMP
    SUPPORTS_COMP=0
    return
end function SUPPORTS_COMP

integer function SUPPORTS_BIN(modelid, id1, id2)
    integer modelid, id1, id2
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_BIN
    SUPPORTS_BIN=0
    return
end function SUPPORTS_BIN

integer function SUPPORTS_PRP1(modelid, propertyid)
    integer modelid, propertyid
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_PRP1
    SUPPORTS_PRP1=0
    if (propertyid==1) SUPPORTS_PRP1=1
    if (propertyid==2) SUPPORTS_PRP1=1
    if (propertyid==3) SUPPORTS_PRP1=1
    if (propertyid==4) SUPPORTS_PRP1=1
    if (propertyid==5) SUPPORTS_PRP1=1
    if (propertyid==6) SUPPORTS_PRP1=1
    if (propertyid==7) SUPPORTS_PRP1=1
    if (propertyid==8) SUPPORTS_PRP1=1
    if (propertyid==12) SUPPORTS_PRP1=1
    if (propertyid==14) SUPPORTS_PRP1=1
    if (propertyid==34) SUPPORTS_PRP1=1 !fluid pressure
    return
end function SUPPORTS_PRP1

integer function SUPPORTS_PRP2(modelid, propertyid)
    integer modelid, propertyid
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_PRP2
    SUPPORTS_PRP2=0
    if (propertyid==4) SUPPORTS_PRP2=1
    if (propertyid==5) SUPPORTS_PRP2=1
    if (propertyid==81) SUPPORTS_PRP2=1
    if (propertyid==89) SUPPORTS_PRP2=1
    if (propertyid==36) SUPPORTS_PRP2=1
    if (propertyid==201) SUPPORTS_PRP2=1
    if (propertyid==202) SUPPORTS_PRP2=1
    if (propertyid==203) SUPPORTS_PRP2=1
    if (propertyid==204) SUPPORTS_PRP2=1
    return
end function SUPPORTS_PRP2

integer function SUPPORTS_PRP3(modelid, propertyid)
    integer modelid, propertyid
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_PRP3
    SUPPORTS_PRP3=0
    return
end function SUPPORTS_PRP3

integer function SUPPORTS_PRP4(modelid, propertyid)
    integer modelid, propertyid
    !DEC$ ATTRIBUTES DLLEXPORT::SUPPORTS_PRP4
    SUPPORTS_PRP4=0
    return
end function SUPPORTS_PRP4

integer function GETNPAR(modelid)
    use DllConst
    integer modelid
    integer nPar
    integer QueryNparPure
    integer QueryNparMix
    !DEC$ ATTRIBUTES DLLEXPORT::GETNPAR
    iEosOpt=modelid
    if (oldRn1.gt.0.and.oldRn2.eq.0) then
        nPar=QueryNparPure()
    elseif (oldRn1.gt.0) then
        nPar=QueryNparMix()
    else
        nPar=0
        res=1
    end if
    GETNPAR=nPar
    return
end function GETNPAR

integer function GETPAR(n, retvalue)
    use DllConst
    integer n
    double precision retvalue
    integer res
    !DEC$ ATTRIBUTES DLLEXPORT::GETPAR
    if (oldRn1.gt.0.and.oldRn2.eq.0) then
        call QueryParPure(1,n,retvalue,res)
    elseif (oldRn1.gt.0) then
        call QueryParMix(n, retvalue, res)
    else
        retvalue=0
        res=1
    end if
    GETPAR=res
    return
end function GETPAR

integer function SETPAR(n, newvalue)
    use DllConst
    integer n
    double precision newvalue
    integer res
    !DEC$ ATTRIBUTES DLLEXPORT::SETPAR
    if (oldRn1.gt.0.and.oldRn2.eq.0) then
        call SetParPure(1, n, newvalue, res)
    elseif (oldRn1.gt.0) then
        call SetParMix(n, newvalue, res)
    else
        res=1
    end if
    SETPAR=res
    return
end function SETPAR

integer function INITIALIZE_MODEL(iEosLocal, Rn1, Rn2, Rn3)
	USE GlobConst
    use DllConst
	IMPLICIT double Precision(A-H,K,O-Z)
    integer iEosLocal, Rn1, Rn2, Rn3
    !DEC$ ATTRIBUTES DLLEXPORT::INITIALIZE_MODEL
    integer ieos, casrn1, casrn2, casrn3, ierr
    INTEGER localCas(nmx)
    IF (LOUD.and.dumpUnit==6) then
	    INITIALIZE_MODEL=5
	    return
    endif
    if (LOUD)write(dumpUnit,*)'INITIALIZE_MODEL: starting'
    INITIALIZE_MODEL=0
	ieos=iEosLocal
    casrn1=Rn1
    casrn2=Rn2
    casrn3=Rn3
    iEosOpt=ieos
    if (oldEOS.ne.0 .and. iEosOpt.ne.oldEOS) then
	    !INITIALIZE_MODEL=4
	    !return
    endif
    ierr=0
	INITIAL=0
        localCas(1)=casrn1
        localCas(2)=casrn2
        localCas(3)=casrn3
    NC=3
    if (Rn2.eq.0) then
        NC=1 !no of components
    elseif (Rn3.eq.0) then
        NC=2 !no of components
    endif
	iErrLoad=0 
    if (NC==1) then
        if (oldRN1.ne.casrn1 .or. oldEOS.ne.ieos) then
	        call PGLWrapperStartup(NC,iEosLocal,localCas,ierLoad)
	        if(iErrLoad > 0 .and. LOUD)write(dumpUnit,*)'InitializeModel: From LoadCritParmsDb, iErrCrit=',iErrCrit
	        if(iErrLoad.NE.0)then
		        INITIALIZE_MODEL=3
		        return
            endif
            oldRN1=casrn1
            oldRN2=0
            oldRN3=0
            oldEOS=ieos
        endif
    elseif (NC==2) then
        if (oldRN1.ne.casrn1 .or. oldRN2.ne.casrn2 .or. oldEOS.ne.ieos) then
	        call PGLWrapperStartup(NC,iEosLocal,localCas,ierLoad)
	        if(iErrLoad > 0 .and. LOUD)write(dumpUnit,*)'InitializeModel: From LoadCritParmsDb, iErrCrit=',iErrCrit
	        if(iErrLoad.NE.0)then
		        INITIALIZE_MODEL=3
		        return
            endif
            oldRN1=casrn1
            oldRN2=casrn2
            oldRN3=0
            oldEOS=ieos
	    endif
    else
        if (oldRN1.ne.casrn1 .or. oldRN2.ne.casrn2 .or. oldRN3.ne.casrn3 .or. oldEOS.ne.ieos) then
	        call PGLWrapperStartup(NC,iEosLocal,localCas,ierLoad)
	        if(iErrLoad > 0 .and. LOUD)write(dumpUnit,*)'InitializeModel: From LoadCritParmsDb, iErrCrit=',iErrCrit
	        if(iErrLoad.NE.0)then
		        INITIALIZE_MODEL=3
		        return
            endif
            oldRN1=casrn1
            oldRN2=casrn2
            oldRN3=casrn3
            oldEOS=ieos
	    endif
    endif
	if (LOUD)write(dumpUnit,*)'InitializeModel: returning. iErr=',INITIALIZE_MODEL
    return
    end function INITIALIZE_MODEL

subroutine CalculateProperty3local(ieos, casrn1, casrn2, casrn3, prp_id, var1, var2, var3, var4, res, ierr)
	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE GlobConst
    use DllConst
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    double Precision res
    integer ieos, casrn1, casrn2, casrn3, prp_id, ierr
    double Precision var1, var2, var3, var4
    double Precision xFrac(nmx),FUGC(nmx),xPure1(nmx),xPure2(nmx),xPure3(nmx) !FUGI requires mole fraction specification because it is written generally for mixtures.
    INTEGER localCas(nmx)
	CHARACTER*77 errMsgPas
!	COMMON/eta/etaL,etaV,ZL,ZV
	CHARACTER($MAXPATH)  CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
    IF (LOUD.and.dumpUnit==6) return
	DEBUG=.FALSE.
	!  Get current directory
	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	!iChange=VERIFY('C:\MYPROJEX\CalcEos',CURDIR)
	!iChange2=VERIFY('C:\MYPROJEX\CALCEOS',CURDIR)
	!if(iChange2.eq.0)iChange=0
	!if(iChange.eq.0 .or. iChange.gt.26)DEBUG=.TRUE.
	NC=3 !no of components
    iEosOpt=ieos
    iProperty=ABS(prp_id)
    localCas(1)=casrn1
    localCas(2)=casrn2
    localCas(3)=casrn3
    !write(*,*)casrn
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (g/cc) given tKelvin
    !iProperty = 3: fluid density (g/cc) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 4: Hres/RT,CpRes/R,CvResR,cmprsblty given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)
    iPhase=1
    if (prp_id==201) iPhase=0   !gas
    if (prp_id==202) iPhase=0   !gas
    if (prp_id==205) iPhase=0   !gas
    res=0
    ierr=0
	
	INITIAL=0

	call IdDipprLookup(NC,localCas,iErrCas,errMsgPas)
	if(iErrCas)then
        ierr=1
        return
    endif
	CALL GETCRIT(NC,iErrCrit)
	if(iErrCrit)then
		ierr=2
		return
	endif
	iErrGet=0 
    if (oldRN1.ne.casrn1 .or. oldRN2.ne.casrn2 .or. oldRN3.ne.casrn3 .or. oldEOS.ne.ieos) then
	    if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	    if(iEosOpt.eq.2)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParms					!Diky model 12
	    if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet) 
	    if(iEosOpt.eq.4)CALL GetEsdCas(NC,localCas,iErrGet)
	    if(iEosOpt.eq.5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms					!Diky model 6
!	    if(iEosOpt.eq.6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	    if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	    if(iEosOpt.eq.8)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	    if(iEosOpt.eq.9)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms
	    if(iEosOpt.eq.10)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters			!Diky model 25
	    if(iEosOpt.eq.11)CALL GetPrTc(NC,iErrGet)		!JRE 2019 : Reads Jaubert's parameters						  !Diky model 22
	    if(iEosOpt.eq.12)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmami					 !Diky model 23
	    if(iEosOpt.eq.13)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParmsEmamiTb					 !Diky model 24
	    if(iEosOpt.eq.14)CALL GetTpt(NC,ID,iErrGet,errMsgPas)!Results placed in common: TptParms, HbParms				    !Diky model 18
	    if(iEosOpt.eq.15)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters			!Diky model 26
	    if(iEosOpt.eq.16)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters			!Diky model 27
	    !NewEos: Add here for initializing parms.
	    if(iErrGet.NE.0)then
		    ierr=3
		    return
        endif
        oldRN1=casrn1
        oldRN2=casrn2
        oldRN3=casrn3
        oldEOS=ieos
	endif
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
        tPrevious = 0
        pPrevious = 0
	    xPure1(2)=0
	    xPure1(3)=0
	    xPure2(1)=0
	    xPure2(3)=0
	    xPure3(1)=0
	    xPure3(2)=0
	    xPure1(1)=1
	    xPure2(2)=1
	    xPure2(3)=1
		!write(52,'(a)')'   T(K)   p(kPa)     x1  phas Vex(cc/mol) rho(mol/cc)  Hex(J/mol)  Gex(J/mol)   ln(gam1)   ln(gam2) '
		!do while(notDone)
			!read(51,*,ioStat=ioErr)tKelvin,pKPa,iPhase,xFrac(1)
        tKelvin = var1
        pKPa = var2
        xFrac(1) = var3
        xFrac(2) = var4
			aPhase='VAP'
			if(iPhase==0)aPhase='GAS'
			if(iPhase==1)aPhase='LIQ'
			xFrac(3) = 1-xFrac(1)-xFrac(2)
			!if(xFrac(2) < 0)pause 'UaWrapper: Error xFrac(1) > 1  ??? '
			!if(ioErr < 0)then   !this means we reached the end of the file
			!	exit            !break the loop and return to remaining execution 
			!elseif(ioErr > 0)then  
			!	write(52,*)'Unexpected error (e.g. type mismatch?) reading line=',line
			!	cycle
			!endif
			!line=line+1
			iErr=0
			pMPa=pKPa/1000
			if( ABS(tKelvin-tPrevious) > 0.1 .or. ABS(pKPa-pPrevious) > 0.1 )then
				!if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' UaWrapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure1(1)
				CALL FugiTP( tKelvin,pMPa,xPure1,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
				if(iErrF > 0)iErr=iErrF
				hRes_RT=uRes_RT+zFactor-1
				vPure1=rGas*tKelvin*zFactor/pMPa
				chemPoPure1=FUGC(1) ! = ln(fPure1/P)
				hRes1 = hRes_RT*rGas*tKelvin ! from module
				!if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' UaWrapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xPure2(1)
				CALL FugiTP( tKelvin,pMPa,xPure2,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
				if(iErrF > 0)iErr=iErr+10*iErrF
				hRes_RT=uRes_RT+zFactor-1
				vPure2=rGas*tKelvin*zFactor/pMPa
				chemPoPure2=FUGC(2) ! = ln(fPure2/P)
				hRes2 = hRes_RT*rGas*tKelvin ! from module
				tPrevious=tKelvin
				pPrevious=pKPa
			endif
			!iPhase=1 !set as default
			!if(expRho/1000 < rhoCritG_cc)iPhase=0
				!if(LOUD)write(*,'(a,f7.2,1x,f9.5,i4,f9.4)')' UaWrapper: calling fugi. T,P,iPhase,x1:',tKelvin,pMPa,iPhase,xFrac(1)
			CALL FugiTP( tKelvin,pMPa,xFrac,NC,iPhase,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrF )
			if(iErrF > 0)iErr=iErr+10*iErrF
			hRes_RT=uRes_RT+zFactor-1
			if(iErr > 0)then
				!write(52,'(1x,2f8.3,F7.4,1x,a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				!if(LOUD)write(*,'(1x,2f8.3,F7.4,1x,a3,1x,6E12.4,i8)')tKelvin,pKPa,xFrac(1),aPhase,errDummy,errDummy,errDummy,errDummy,errDummy,errDummy,iErr
				!cycle
                res = 0
                return
			endif
																														 
			vMix=rGas*tKelvin*zFactor/pMPa
			rhoMol_cc=1/vMix
			activity1=FUGC(1)-chemPoPure1 ! = ln(fi/xi*fiPure) 
			activity2=FUGC(2)-chemPoPure2 ! = ln(fi/xi*fiPure)
			hResMix = hRes_RT*rGas*tKelvin ! from module
			vXsCc_mol = vMix - (xFrac(1)*vPure1	+ xFrac(2)*vPure2)
			hXsJ_mol = hResMix - (xFrac(1)*hRes1 + xFrac(2)*hRes2)
			gXsJ_mol =  (xFrac(1)*activity1 + xFrac(2)*activity2)
			!write(52,'(1x,2f8.3,F7.4,1x,a3,1x,7E12.4)')tKelvin,pKPa,xFrac(1),aPhase,vXsCc_mol,rhoMol_cc,hXsJ_mol,gXsJ_mol,activity1,activity2
            if (iProperty==4) res = exp(activity1)
            if (iProperty==5) res = exp(activity2)
            if (iProperty==81) res = 1000*rhoMol_cc
            if (iProperty==89) res = 0.000001*vXsCc_mol
            if (iProperty==36) res = 0.001*hXsJ_mol
            if (iProperty==201) res = exp(FUGC(1))
            if (iProperty==202) res = exp(FUGC(2))
            if (iProperty==203) res = exp(FUGC(1))
            if (iProperty==204) res = exp(FUGC(2))
            if (iProperty==205) res = exp(FUGC(3))
            if (iProperty==206) res = exp(FUGC(3))
		!enddo !	while(notDone)
        return
end subroutine CalculateProperty3local

integer function SETSTRING(tag, value)
	USE GlobConst
    character*255 tag, value !, local
    !DEC$ ATTRIBUTES DLLEXPORT::SETSTRING
    if(LOUD)write(dumpUnit,*)'SETSTRING: value,tag',value,' ',tag 
    if (tag(1:8).eq.'LOCATION') then
        do i1=1,255
            if (value(i1:i1).eq.'|'.and. .NOT. LOUD) then !if(LOUD), assume debugging and use hard coded PGLInputDir
                masterDir=value(1:i1-1)
                PGLInputDir=trim(masterDir)//'\input'
                goto 666
            end if
        enddo
666     continue        
        SETSTRING=1
    else
        SETSTRING=0
    endif
    if(LOUD)write(dumpUnit,*)'SetString: returning. SETSTRING=',SETSTRING 
    return
    end function SETSTRING
