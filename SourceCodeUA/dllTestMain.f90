!MAIN PROGRAM FOR CALLING PGLDLL.DLL. This is untested code using iso_c_binding
    USE GlobConst
    use iso_c_binding
    character(255) errMsg 
    character*255 tag, value,hello !, local
    integer i1, cnt, clen, status, nbsl
    character arg*100,cmd*100

interface
    integer function INITPGLDLL(errMsg)
        !DEC$ ATTRIBUTES DLLIMPORT, ALIAS:"INITPGLDLL" :: InitPGLDLL
        character(255) errMsg
    end function INITPGLDLL
	integer function SetLoudTrue(hello)
    !DEC$ ATTRIBUTES, ALIAS:"SETLOUDTRUE" :: SetLoudTrue
		character(255) hello
	end function SetLoudTrue
    integer function iSetDumpUnit(aPlace)
        character(4) aPlace
        !DEC$ ATTRIBUTES DLLIMPORT::iSetDumpUnit
	end function iSetDumpUnit
    integer(c_int) function ISETMASTERDIR(input) bind(C, name="ISETMASTERDIR")
        use iso_c_binding
        character(kind=c_char), dimension(*) :: input
    end function ISETMASTERDIR
	integer function INITIALIZE_MODEL(iEosLocal, Rn1, Rn2, Rn3)
		!DEC$ATTRIBUTES :: INITIALIZE_MODEL
		integer iEosLocal, Rn1, Rn2, Rn3
	end function INITIALIZE_MODEL
	integer function Calculate1(casrn1, modelid, propertyid, t, p, res, uncert)
		!DEC$ATTRIBUTES :: Calculate1 
		integer modelid, casrn1, propertyid, localprpid, ierr
		double Precision t, p, res, uncert
	end function Calculate1
	DoublePrecision function CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
		!DEC$ATTRIBUTES :: CalculateProperty 
		integer ieos, casrn, prp_id, ierr
		double Precision var1, var2
	end function CalculateProperty
	DoublePrecision function CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
		!!!DEC$ATTRIBUTES :: CalculateProperty2 
		integer ieos, casrn1, casrn2, prp_id, ierr
		double Precision var1, var2, var3
	end function CalculateProperty2
END INTERFACE
    call get_command (cmd, clen, status)
    call get_command_argument (0, cmd, clen, status)
    cnt = command_argument_count ()
    do i1 = 1, cnt
        call get_command_argument (i, arg, clen, status)
    end do
    nbsl=0
    do i1=100,1,-1
        if (cmd(i1:i1).eq.'\') then
            nbsl=nbsl+1
            if (nbsl.eq.2) then
                MasterDir=cmd(1:i1-1)
                PGLInputDir=trim(masterDir)//'\input'
            exit
            end if
        end if
    enddo
    ! Above code is for running from command line. For debugging&testing in VS Studio:
    MasterDir='c:\PGLWrapper'
    PGLInputDir=trim(MasterDir)//'\input'
    iErr=iSetMasterDir(MasterDir) !Copy the location of MasterDir into the DLL GlobConst.       
    !tag='LOCATION'
    !value='c:\PGLWrapper|'
    !iErr=SetString(tag,value)
    iErr=iSetLoudTrue(errMsg) ! LOUD=.TRUE. => debug info 
    iDumpUnit=iSetDumpUnit('FILE')
    !if (LOUD)iDumpUnit=iSetDumpUnit('FILE')
	iErr=INITPGLDLL(errMsg)
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
    pause 'Check the output. It should say: 34.47, 0.762,0.0095,0.762,5.75'
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
    write(*,*) 'Test2: iErr,  T(K)   P(kPa)   x1    prp'
    write(*,'(i4,1x,2f9.2,1x,f8.5,1x,E11.4)') iErr,var1,var2,var3,prp
    var2 = 91911.609110881604
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    !write(*,*) 'Test2: iErr,  T(K)   P(kPa)   x1    prp'
    write(*,'(i4,1x,2f9.2,1x,f8.5,1x,E11.4)') iErr,var1,var2,var3,prp
    var2 = 1000
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    !write(*,*) 'Test2: iErr,  T(K)   P(kPa)   x1    prp'
    write(*,'(i4,1x,2f9.2,1x,f8.5,1x,E11.4)') iErr,var1,var2,var3,prp
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
	