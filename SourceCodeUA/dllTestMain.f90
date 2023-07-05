!MAIN PROGRAM FOR CALLING PGLDLL.DLL. This is untested code using iso_c_binding
!	use iso_c_binding
!    DoublePrecision xFrac(NMX),FUGC(NMX) !FUGI requires mole fraction specification because it is written generally for mixtures.
!    INTEGER ier(12),localCas(NMX)
!    DoublePrecision CalcResult
!    integer ieos, casrn, prp_id, ierr
!    doublePrecision var1, var2
!    character*255 tag, value
!
!	interface
!		function Calculate1() bind(c)	!
!		  integer(c_int) :: Calculate1
!		end function Calculate1
!	end interface
!
!	! Declare the DLL function
!	procedure(Calculate1), bind(c) :: PGLDLL.DLL
!
!	! Load the DLL
!	call c_f_pointer(c_loc(PGLDLL.DLL), Calculate1)
!
!	! Call the DLL function
!	print *, Calculate1()
!	stop
!end !DLLTestMain
	
!MAIN PROGRAM FOR CALLING PGLDLL.DLL. This is untested code using iso_c_binding
    !DoublePrecision xFrac(NMX),FUGC(NMX) !FUGI requires mole fraction specification because it is written generally for mixtures.
    !INTEGER ier(12),localCas(NMX)
    !DoublePrecision CalcResult
    !integer ieos, casrn, prp_id, ierr
    !doublePrecision var1, var2
    !character(255) tag, value
	character(255) errMsg 

INTERFACE
	integer function Activate(hello)
		!MS$ATTRIBUTES DLLIMPORT :: Activate 
		character(255) hello
	end function Activate
	integer function INITIALIZE_MODEL(iEosLocal, Rn1, Rn2, Rn3)
		!MS$ATTRIBUTES DLLIMPORT :: INITIALIZE_MODEL 
		integer iEosLocal, Rn1, Rn2, Rn3
	end function INITIALIZE_MODEL
	integer function Calculate1(casrn1, modelid, propertyid, t, p, res, uncert)
		!MS$ATTRIBUTES DLLIMPORT :: Calculate1 
		integer modelid, casrn1, propertyid, localprpid, ierr
		double Precision t, p, res, uncert
	end function Calculate1
	DoublePrecision function CalculateProperty(ieos, casrn, prp_id, var1, var2, ierr)
		!MS$ATTRIBUTES DLLIMPORT :: CalculateProperty 
		integer ieos, casrn, prp_id, ierr
		double Precision var1, var2
	end function CalculateProperty
end interface
	iErr=Activate(errMsg)
	if(iErr > 0)then
		write(*,'(a)')' DLLTestMain: Activate! ErrMsg=',errMsg
		write(*,*)' DLLTestMain: iErr=',iErr
	else
		write(*,*)'DLLTestMain: Activate returned iErr=0. Woohoo!'
	endif
	call Test1
	!call Test2
	call Test3
	stop
END !main Program

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
    prp= CalculateProperty(ieos, Rn1, prp_id, var1, var2, ierr)
    !CalcResult = CalculateProperty(ieos, Rn1, prp_id, var1, var2, ierr)
!    prp_id=8
!    ierr=Calculate1(Rn1, ieos, prp_id, var1, var2, prp, uncert)
    write(*,*) 'Test1: prp=',prp
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
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp=',var2,var3,prp
    var2 = 91911.609110881604
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp=',var2,var3,prp
    var2 = 1000
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test2: P,x1,prp=',var2,var3,prp
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
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test3: P,x1,prp=',var2,var3,prp
    var2 = 91911.609110881604
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test3: P,x1,prp=',var2,var3,prp
    var2 = 1000
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3,  ierr)
    write(*,*) 'Test3: P,x1,prp=',var2,var3,prp
    return

    ieos = 2
    casrn1 = 123911
    casrn2 = 127184
    var1 = 298.15
    var2 = 101.3
    var3 = 0.085
    prp_id = 36
    prp=CalculateProperty2(ieos, casrn1, casrn2, prp_id, var1, var2, var3, ierr)
    write(*,*) 'Test3: P,x1,prp=',var2,var3,prp
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
	