	!Subroutine MEMTest()  ! Comment this line out to make this the main program
	USE GlobConst
    USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

    DoublePrecision xFrac(NMX)
    dimension rLnPhiAssoc(NMX)
	
    nTypes(1)=4
    nTypes(2)=3
    nDegree(1,4)=1
    nDegree(2,3)=1
    nDonors(1,4)=1
    nDonors(2,3)=1
    nAcceptors(2,3)=1
    eDonorKcal_mol(1,4)=2.3
    eDonorKcal_mol(2,3)=6
    eAcceptorKcal_mol(2,3)=6
	bondVolNm3(1,4)=0.0007d0
	bondVolNm3(2,3)=0.0005d0
	isTpt=.TRUE.
    bVolCC_mol(2)=.0665308*602.22
    bVolCC_mol(1)=.0901502*602.22
    
    isZiter=0
    tKelvin=300
    xFrac=(0.5,0.5)
    nComps=2
    rhoMol_cc=0.01778
    write(*,*)xFrac(1),xFrac(2)
    call  MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr)
    write(*,*)zAssoc,aAssoc,uAssoc
	stop
	END
