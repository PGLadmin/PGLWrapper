!MAIN PROGRAM FOR CALLING PGLDLL.DLL. This is untested code using iso_c_binding
	character(255) errMsg 
    character*255 tag, value,hello !, local
    LOGICAL LOUD

INTERFACE
	integer function InitPGLDLL(hello)
		!DEC$ATTRIBUTES DLLIMPORT :: InitPGLDLL  
		character(255) hello
	end function InitPGLDLL
	integer function SetLoudTrue(hello)
		!DEC$ATTRIBUTES DLLIMPORT :: SetLoudTrue  
		character(255) hello
	end function SetLoudTrue
    integer function iSetDumpUnit(aPlace)
        character(4) aPlace
        !DEC$ ATTRIBUTES DLLIMPORT::iSetDumpUnit
	end function iSetDumpUnit
	integer function INITIALIZE_MODEL(iEosLocal, Rn1, Rn2, Rn3)
		!DEC$ATTRIBUTES DLLIMPORT :: INITIALIZE_MODEL
		integer iEosLocal, Rn1, Rn2, Rn3
	end function INITIALIZE_MODEL
	integer function Calculate1(casrn1, modelid, propertyid, t, p, res, uncert)
		!DEC$ATTRIBUTES DLLIMPORT :: Calculate1 
		integer modelid, casrn1, propertyid, localprpid, ierr
		double Precision t, p, res, uncert
	end function Calculate1
	DoublePrecision function CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
		!DEC$ATTRIBUTES DLLIMPORT :: CalculateProperty 
		integer ieos, casrn, prp_id, ierr
		double Precision var1, var2
	end function CalculateProperty
END INTERFACE
    !tag='LOCATION'
    !value='c:\PGLWrapper|'
    !iErr=SetString(tag,value)
    iErr=iSetLoudTrue(errMsg)
    iDumpUnit=iSetDumpUnit('FILE')
    !if (LOUD)iDumpUnit=iSetDumpUnit('FILE')
	iErr=InitPGLDLL(errMsg)
	if(iErr > 10)then
		write(*,'(a)')' DLLTestMain: InitPGLDLL! ErrMsg=',errMsg
		write(*,*)' DLLTestMain: iErr=',iErr
	else
		write(*,*)'DLLTestMain: InitPGLDLL returned iErr=0. Woohoo!'
        !write(*,*)'DLLTestMain: PGLInputDir=',TRIM(PGLInputDir)
    endif
    
	call Test1
	call Test2
	call Test3

    stop
END !main Program

subroutine Test1
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	integer ieos, Rn1, Rn2, Rn3, prp_id, ierr
	doublePrecision var1, var2, prp, uncert
	DoublePrecision CalcResult
    ieos=2
    Rn1=67641   !acetone
    Rn2=0
    Rn3=0
    prp_id=1    !vapor pressure
    var1=300.0  !tKelvin
    var2=0.0
    ierr=0
    ierr=INITIALIZE_MODEL(ieos, Rn1, Rn2, Rn3)
    prp= CalculateProperty(ieos, Rn1, prp_id, var1, var2, ierr)
    !CalcResult = CalculateProperty(ieos, Rn1, prp_id, var1, var2, ierr)
!    prp_id=8
!    ierr=Calculate1(Rn1, ieos, prp_id, var1, var2, prp, uncert)
    write(*,*) 'Test1: T,prp=',var1,prp
end subroutine

subroutine Test2
!            if (prp_id==4) res = exp(activity1)
!            if (prp_id==5) res = exp(activity2)
!            if (prp_id==6) res = Rgas*(xFrac(1)*(activity1)+xFrac(2)*(activity2))
!            if (prp_id==36) res = hXsJ_mol/1000	!kJ/mol
!            if (prp_id==81) res = 1000*rhoMol_cc	!mol/L
!            if (prp_id==89) res = vXsCc_mol/1D6	!m^3/mol
!            if (prp_id==201) res = exp(FUGC(1))	! vapor
!            if (prp_id==202) res = exp(FUGC(2))	! vapor
!            if (prp_id==203) res = exp(FUGC(1))	! liq
!            if (prp_id==204) res = exp(FUGC(2))	! liq
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	integer ieos, casrn1, casrn2, prp_id, ierr
	doublePrecision var1, var2, var3, prp
    ieos=2
    casrn1 = 75150      !CarbonDisulfide
    casrn2 = 2551624    !SulfurHexafluoride
    var1 = 298.136
    var2 = 1000
    var3 = 0.0013779105351375952
    prp_id = 201        !vapor chemical potential of component 1. 
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: iErr,P,x1,prp,=',iErr,var2,var3,prp
    var2 = 91911.609110881604
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp,iErr=',var2,var3,prp,iErr
    var2 = 1000
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp,iErr=',var2,var3,prp,iErr
    return

    ieos = 2
    casrn1 = 123911
    casrn2 = 127184
    var1 = 298.15
    var2 = 101.3
    var3 = 0.085
    prp_id = 36
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp=',var2,var3,prp
    return

    ieos = 2
    casrn1 = 75525
    casrn2 = 75058
    var1 = 298.136
    var2 = 101
    var3 = 0.9428
    prp_id = 89
    prp= CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    write(*,*) 'Test2: P,x1,prp=',var2,var3,prp
end subroutine

subroutine Test3	!testing dll calls
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    integer i, ieos, casrn1, casrn2, casrn3, prp_id, ierr,nPar
    doublePrecision var1, var2, var3, prp, p
    integer GETNPAR,GETPAR
	i=INITIALIZE_MODEL(2,71363,7732185,0)       ! (iEos,casrn1,casrn2,casrn3)
	nPar=GETNPAR(2)
	i=GETPAR(1,p)
    write(*,*) 'Test3: nPar,par=',nPar,p
	return
    ieos=2
    casrn1 = 123911     !14dioxane
    casrn2 = 127184     !tetrachloroethene
    var1 = 298.15
    var2 = 101.3
    var3 = 0.085
    prp_id = 36         !Hxs
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    write(*,*) 'Test3: P,x1,prp,iErr=',var2,var3,prp,iErr
    return

    ieos = 2
    casrn1 = 75525
    casrn2 = 75058
    var1 = 298.136
    var2 = 101
    var3 = 0.9428
    prp_id = 89
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test3: P,x1,prp=',var2,var3,prp
end subroutine
	