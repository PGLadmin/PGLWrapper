	!Subroutine MEMTest()  ! Comment this line out to make this the main program
	USE GlobConst
    USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)

    DoublePrecision xFrac(NMX)
    dimension rLnPhiAssoc(NMX)
	!Example. Acetone(1)+chloroform(2) with speadmd model.
	nTypes(1)=3
	nTypes(2)=4
	nDegree(1,3)=1
	nDegree(2,1)=1
	nDonors(1,3)=1
	nDonors(2,1)=0
	nAcceptors(1,3)=1
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(1,3)=1.0
	eAcceptorKcal_mol(2,1)=1.4
	eDonorKcal_mol(1,3)=1.6
	eDonorKcal_mol(2,1)=0
	bondVolNm3(1,3)=0.0025d0
	bondVolNm3(2,1)=0.0025d0
	isTpt=.TRUE.
	bVolCC_mol(2)=35.96
	bVolCC_mol(1)=35.34
	rhoMol_cc=0.013 ! this is a reasonable estimate at all compositions for speadmd. 
	isZiter=0
	tKelvin=300
	xFrac=(0.5,0.5)
	nComps=2
	write(*,*)xFrac(1),xFrac(2)
	call  MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr)
	write(*,*)zAssoc,aAssoc,uAssoc
	stop
	END
