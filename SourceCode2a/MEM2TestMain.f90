	!Subroutine MEMTest()  ! Comment this line out to make this the main program
	USE GlobConst
    USE Assoc
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    DoublePrecision xFrac(NMX)
    DoublePrecision rLnPhiAssoc(NMX),rLnPhiAssoc1(NMX),rLnPhiAssoc2(NMX),gam(NMX) !,vLiq(NMX)

	Call Example3Setup(iErr)

	isZiter=0
	tKelvin=300
	nComps=2
	xFrac(2)=zeroTol
	xFrac(1)=1-xFrac(2)
	call MEM2(isZiter,tKelvin,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1,uAssoc1,rLnPhiAssoc1,iErr)
	write(*,'(a)')' xFrac(1),GE_RT,rhoMol_cc,gam(1),gam(2),UE_RT,zAssoc,aAssoc,chemPoAssoc(1),chemPoAssoc(2)'
	do i=0,10,1
		xFrac(1)=i/10.d0
		if(xFrac(1) < zeroTol)xFrac(1)=zeroTol
		xFrac(2)=1-xFrac(1)
		if(xFrac(2) < zeroTol)xFrac(2)=zeroTol
		bMix=SUM( xFrac(1:nComps)*bVolCC_mol(1:nComps) )
		vMix=SUM( xFrac(1:nComps)*vLiq(1:nComps) )
		rhoMol_cc=1/vMix
		call MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr)
		if(i==0)then
			zAssoc2=zAssoc
			aAssoc2=aAssoc
			uAssoc2=uAssoc
			rLnPhiAssoc2(1:nComps)=rLnPhiAssoc(1:nComps)
		endif
		gam(1)=exp( rLnPhiAssoc(1)-zAssoc*bVolCC_mol(1)/bMix-(rLnPhiAssoc1(1)-zAssoc1) )
		gam(2)=exp( rLnPhiAssoc(2)-zAssoc*bVolCC_mol(2)/bMix-(rLnPhiAssoc2(2)-zAssoc2) )
		GE_RT=xFrac(1)*LOG(gam(1))+xFrac(2)*LOG(gam(2))
		UE_RT=aAssoc-(xFrac(1)*aAssoc1+xFrac(2)*aAssoc2)
		write(*,'(1x,3f8.4,12f8.2)')xFrac(1),GE_RT,rhoMol_cc,gam(1),gam(2),UE_RT,zAssoc,aAssoc,rLnPhiAssoc(1),rLnPhiAssoc(2)
		if(i==0)gamInf1=gam(1)
		if(i==10)gamInf2=gam(2)
	enddo
	write(*,610)'T(K),IDAC1,IDAC2=',tKelvin,gamInf1,gamInf2
610	format(1x,a,12E12.4)
	pause 'Success! Check results above.'
	stop
	END

	Subroutine Example1Setup(iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	iErr=0
	!Example. Acetone(1)+chloroform(2) with speadmd model.
	isTpt=.TRUE.				!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=3
	nDegree(1,3)=1
	nDonors(1,3)=1
	nAcceptors(1,3)=1
	eAcceptorKcal_mol(1,3)=1.0
	eDonorKcal_mol(1,3)=1.6
	bondVolNm3(1,3)=0.0025d0
	bVolCC_mol(1)=36			!cf. ParmsTpt.txt
	vLiq(1)=74					!cf. ParmsPrTcJaubert.txt
	nTypes(2)=4
	nDegree(2,1)=1
	nDonors(2,1)=0
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=1.4
	eDonorKcal_mol(2,1)=0
	bondVolNm3(2,1)=0.0025d0
	bVolCC_mol(2)=38
	vLiq(2)=81
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example2Setup(iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	iErr=0
	!Example. THF(1)+methanol(2) with speadmd model.
	nTypes(1)=1
	nTypes(2)=1
	nDegree(1,1)=1
	nDegree(2,1)=1
	nDonors(1,1)=1
	nDonors(2,1)=1
	nAcceptors(1,1)=0
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(1,1)=0		!cf.ParmsHB4.txt
	eAcceptorKcal_mol(2,1)=5
	eDonorKcal_mol(1,1)=1.3d0
	eDonorKcal_mol(2,1)=5
	bondVolNm3(1,1)=0.0007d0
	bondVolNm3(2,1)=0.0010d0
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	bVolCC_mol(1)=.0702793d0*avoNum	!cf. ParmsTpt.txt
	bVolCC_mol(2)=.0310936d0*avoNum
	vLiq(1)=59.11d0/0.7214d0		!cf. ParmsPrTcJaubert.txt
	vLiq(2)=32.04d0/0.7896d0		
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example3Setup(iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	iErr=0
	!Example. Dioxane(1)+chloroform(2) with speadmd model. cf. EL2ed Example 19.3, p779. 
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!methanol...
	nDegree(1,1)=2
	nDonors(1,1)=1
	nAcceptors(1,1)=0
	eAcceptorKcal_mol(1,1)=0		!cf.ParmsHB4.txt, adapted from THF
	eDonorKcal_mol(1,1)=1.3d0
	bondVolNm3(1,1)=0.0007d0
	bVolCC_mol(1)=.0748917d0*avoNum	!cf. ParmsTpt.txt
	vLiq(1)=62.11d0/0.7251d0		!cf. ParmsPrTcJaubert.txt
	nTypes(2)=4						!chloroform...
	nDegree(2,1)=1
	nDonors(2,1)=0
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=1.4
	eDonorKcal_mol(2,1)=0
	bondVolNm3(2,1)=0.0025d0
	bVolCC_mol(2)=38
	vLiq(2)=81
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example4Setup(iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	iErr=0
	!Example. Dioxane(1)+methanol(2) with speadmd model. cf. EL2ed Example 19.3, p779. 
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1
	nDegree(1,1)=2
	nDonors(1,1)=1
	nAcceptors(1,1)=0
	eAcceptorKcal_mol(1,1)=0		!cf.ParmsHB4.txt, adapted from THF
	eDonorKcal_mol(1,1)=1.3d0
	bondVolNm3(1,1)=0.0007d0
	bVolCC_mol(1)=.0748917d0*avoNum	!cf. ParmsTpt.txt
	vLiq(1)=62.11d0/0.7251d0		!cf. ParmsPrTcJaubert.txt
	nTypes(2)=1
	nDegree(2,1)=1
	nDonors(2,1)=1
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=5
	eDonorKcal_mol(2,1)=5
	bondVolNm3(2,1)=0.0010d0
	bVolCC_mol(2)=.0310936d0*avoNum
	vLiq(2)=32.04d0/0.7896d0		
	return
	end
