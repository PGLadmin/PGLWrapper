	!Subroutine MEMTest()  ! Comment this line out to make this the main program
	USE GlobConst
    USE Assoc
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    DoublePrecision xFrac(NMX),IDAC(NMX),gamInf(NMX)
    DoublePrecision rLnPhiAssoc(NMX),rLnPhiAssoc1(NMX),rLnPhiAssoc2(NMX),gam(NMX) !,vLiq(NMX)
	!LOUD=.TRUE.
	dumpUnit=6
	isZiter=0

	tKelvin=300
	Call Example3Setup(tKelvin,IDAC,iErr)

	nComps=2
	xFrac(2)=zeroTol
	xFrac(1)=1-xFrac(2)
	tStep=1.d-4
	tPlus =tKelvin*(1+tStep)
	tMinus=tKelvin*(1-tStep)
	call MEM2(isZiter,tPlus,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1p,uAssoc1,rLnPhiAssoc1,iErr)
	call MEM2(isZiter,tMinus,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1,uAssoc1,rLnPhiAssoc1,iErr)
	uAssoc1= -tKelvin*(aAssoc1p-aAssoc)/(tPlus-tMinus)
	write(*,'(a)')' xFrac(1),GE_RT,UE_RT,rhoMol_cc,gam(1),gam(2),zAssoc,uAssoc,chemPoAssoc(1),chemPoAssoc(2)'
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
		GE_RT=aAssoc-(xFrac(1)*aAssoc1+xFrac(2)*aAssoc2)
		UE_RT=uAssoc-(xFrac(1)*uAssoc1+xFrac(2)*uAssoc2)
		write(*,'(1x,4f8.4,12f8.3)')xFrac(1),GE_RT,UE_RT,rhoMol_cc,gam(1),gam(2),zAssoc,uAssoc,rLnPhiAssoc(1),rLnPhiAssoc(2)
		if(i== 0)gamInf(1)=gam(1)
		if(i==10)gamInf(2)=gam(2)
		if(i==4.and.LOUD)pause 'Main: halfway there! Check so far.'
	enddo
	write(*,610)'T(K),IDAC1,IDAC2=',tKelvin,gamInf(1:nComps),IDAC(1:nComps)
610	format(1x,a,12E12.4)
	pause 'Success! Check results above.'
	stop
	END

	Subroutine Example1Setup(tKelvin,IDAC,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDAC(NMX)
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
	IDAC(1:2)=1 				!default initialization
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example2Setup(tKelvin,IDAC,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDAC(NMX)
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
	IDAC(1:2)=1 				!default initialization
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example3Setup(tKelvin,IDAC,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDAC(NMX)
	iErr=0
	!Example. Dioxane(1)+chloroform(2) with speadmd model. cf. EL2ed Example 19.3, p779. 
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!methanol...
	nDegree(1,1)=1
	nDonors(1,1)=2
	nAcceptors(1,1)=0
	eAcceptorKcal_mol(1,1)=0		!cf.ParmsHB4.txt, adapted from THF
	eDonorKcal_mol(1,1)=1.3d0
	bondVolNm3(1,1)=0.0007d0
	bVolCC_mol(1)=.0748917d0*avoNum	!cf. ParmsTpt.txt
	vLiq(1)=62.11d0/0.7251d0		!cf. ParmsPrTcJaubert.txt
	nTypes(2)=1						!chloroform...
	nDegree(2,1)=1
	nDonors(2,1)=0
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=1.4
	eDonorKcal_mol(2,1)=0
	bondVolNm3(2,1)=0.0025d0
	bVolCC_mol(2)=38
	vLiq(2)=81
	Call RdfCalc(rdfContact2,dAlpha,bVolCC_mol(2)/vLiq(2))
	ralph1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	ralph2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralph1inf1,ralph2inf1=',ralph1inf1,ralph2inf1
	rhoDelta=ralph1inf1*ralph2inf1
	IDAC(1)= -nDegree(1,1)*nDonors(1,1)*LOG(1+nDegree(2,1)*nAcceptors(2,1)*rhoDelta)						   !Eq. 64, Only good for acetone+chloroform, dioxane+chloroform, ...
	Call RdfCalc(rdfContact1,dAlpha,bVolCC_mol(1)/vLiq(1))
	ralph1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralph2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(*,610)'ralph1inf2,ralph2inf2=',ralph1inf2,ralph2inf2
	rhoDelta=ralph1inf2*ralph2inf2
	IDAC(2)= -nDegree(2,1)*nAcceptors(2,1)*LOG(1+nDegree(1,1)*nDonors(1,1)*rhoDelta)
	IDAC(1:2)=EXP(IDAC(1:2))
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example4Setup(tKelvin,IDAC,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDAC(NMX)
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
	IDAC(1:2)=1 				!default initialization
610	format(1x,a,12E12.4)
	return
	end
