MODULE DllConst
	!IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!SAVE
	!PUBLIC
	Implicit NONE
	Integer oldRN1,oldRN2,oldRN3,oldEOS
END MODULE DllConst

double Precision function FORTRAN_DLL1(i1, d1)
    integer i1
    double Precision d1

  ! Expose subroutine FORTRAN_DLL1 to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::FORTRAN_DLL1

  ! Variables
    FORTRAN_DLL1 = i1*d1
    return
end function FORTRAN_DLL1
