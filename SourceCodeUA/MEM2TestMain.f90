!	     Here we consider application of the master equation method (MEM) by considering several examples.			
!	MEM has been discussed by Elliott (2022). Here, MEM1 is a special case called by MEM2 when applicable.			
!	The parameters are taken from the SPEADMD model extant 2022. This model is characterized by 			
!	two sets of parameters: ParmsHB4.txt and ParmsTpt.txt. The ParmsHB file characterizes the bonding			
!	energies and volumes of the potential models adapted by Wertheim (1984) in his TPT1 theory.			
!	Wertheim's potential can be envisioned as "blisters" of square-well energy decorating the surface			
!	of a hard sphere. This kind of potential suffices to represent the short-range and sterospecific nature			
!	of hydrogen bonding and similar complexation interactions. The results of TPT1 can be shown to be 			
!	equivalent to chemical reaction theory with specific assumptions about heat capacity and reaction 			
!	equilibrium constants. Hence we refer to this generally as a "chemical" (^chem) contribution.			
!	     The discussion here pertains to direct computation of the chemical contributions to activity coefficients.			
!	Dispersion and polarity interactions would need to be added through implementing other theories.			
!	SPEADMD was chosen as the basis for characterization (as opposed to ESD or PC-SAFT for example)			
!	because its molecular volumes and Helmholtz energy perturbation contributions (A0,A1,A2,…) are			
!	derived from molecular simulations of a well-defined transferable potential model. In the ParmsTpt			
!	file, the results of the tranferable potential model have been refined to match experimental data.			
!	This refinement is achieved by scaling mii and kii in bVol(i)=bVol^TFF*(1-mii) and A1=A1^TFF*(1-kii),			
!	A2=A2^TFF*(1-kii), where bvol^TFF is the molecular volume from simulation of the transferable force			
!	field. The parameters mii and kii are optimized by minimizing deviations from saturation pressure and			
!	liquid density data. Only the scaled values of bVol and A_i parameters are stored in ParmsTpt.txt.			
!	It is presumed that simulations with other transferable potential models (e.g., TraPPE-UA, TAMie-AUA, 			
!	OPLS-AA, ab initio/ML, ...) would yield qualitatively similar values for bVol and {A_i}. Time will tell.			
!	     The examples cover: (1) acetone+chloroform (2) THF+chloroform (3) 14dioxane+chloroform 			
!	(4) 135 trioxane+chloroform (5) methanol+THF (6) acetone+methanol. Example 1 might look unexpected			
!	for three reasons: (1) It uses the complete site type descriptions from SPEADMD. (2) Acetone is represented			
!	as being a weak electron acceptor. The point of the complete site type description is that multiple			
!	site types can be represented in the same molecule with MEM2. For acetone, only one of those site			
!	types has nonzero chemical contributions, but ethanolamine would provide an example of a molecule			
!	with two site types. The weak acceptor of acetone is supported by the shape of its vapor pressure as well			
!	spectroscopic observations (cf. Abraham, 1993). The values in ParmsHB4(2022) are only preliminary.			
!	Because acetone associates in this model, this is an example of case 5.3 in Elliott (2022) for the infinite			
!	dilution activity coefficients (IDACs). This differs from the ESD analysis by Elliott (2022). Note that			
!	the formula for estimating the IDAC differs from the complete calculation because the estimation 			
!	formula makes approximations for FA and FD. 			
!	    Examples 2-4 illustrate the effect of increasing the number of donors on one molecule while holding			
!	the second molecule constant (as chloroform). These are cases of section 5.2. The chemical contributions			
!	of all three examples are assumed the same except for the number of donors. The dioxane and trioxane 			
!	are approximated as having the same molecular volume. As expected, increasing donor number decreases			
!	the IDAC of the donor compound. The chloroform IDAC decreases slightly less than the donor compound.			
!	A single site type is assumed for each molecule in all examples except Example 1.			
!	    Example 5 covers a case similar to Example 1 but THF is an electron donor whereas chlroform is an 			
!	acceptor. This is still covered by case 5.3, but the weak electron acceptor of acetone causes ambiguous			
!	mixing, where the solvation interaction is slightly weaker than the association of acetone.			
!	    Example 6 is a case for which the IDAC cannot be expressed as a simple limit because acetone's			
!	acceptor is weak, but nonzero. It is somewhere between cases 5.1 and 5.3. The IDAC expression in the 			
!	setup routine appies the symmetric association expression (case 5.1).			
!	    We can gain more practical insight by considering a couple of the examples in detail. The analysis			
!	below focuses on Examples 1 and 5. 			
!	Example 1. Acetone(1)+chloroform(2)			
!	Table 1	eAcc(kcal/mol)	eDon(kcal/mol)	
!	acetone	1	1.6	
!	chloroform	1.4	0	
!	From this small table, it should be apparent that the acid/base solvation interaction is stronger than acetone’s association interaction. A simple way to think about these interactions is:			
!	KADij ~ sqrt{ [exp(eAcc(i)/RT)-1]* [exp(eDon(j)/RT)-1] } , ie. the geometric mean.			
!	It’s hard to think in terms of exponentials, however. This suggests tabulating values of:			
!	sqrt[exp(energy/RT)-1], where the “energy” can be either electron acceptor or electron donor.  			
!	Let’s call this a tabulation of acidity(alpha^L) and basicity(beta^L), where alpha^L=sqrt[exp(eAcc/RT)-1]:			
!	The superscript "L" serves as a reminder of the Lewis acid/base perspective applied here.			
!	For acetone+chloroform at		300	K,
!	Table 2	alpha	beta	
!	acetone	2.09	3.69	
!	chloroform	3.08	0	
!	Note that the “A” for acceptor goes with the “alpha” for acidity. If you keep that straight, then the D and beta must go together my process of elimination.			
!	Then KADij=alphai*betaj			
!	Table 3	KADij		
!	acetone	7.7	0	
!	chloroform	11.36	0	
!	In the KADij matrix (where i designates row and j designates column), we see that chloroform has zero “acid(1)*base(2)” and “acid(2)*base(2)” because base(2) is zero.			
!	We also see that the solvation interaction is 48% stronger than the association interaction. That means that these molecules will prefer solvation and the GE should be negative.			
!	These observations are validated by the computed IDACs, 0.81, 0.78, both of which are less than 1. 			
!	Next, let’s consider MeOH(1)+THF(2). From ParmsHB4, we have			
!		T(K) =	300	
!	Table 4	eAcc	eDon	
!	Methanol	5	5	
!	THF	0	1.3	
!		alpha	beta	
!	Methanol	66.25	66.25	
!	THF	0	2.8	
!		KADij		
!	Methanol	4388.67	185.63	
!	THF	0	0	
!	In this case, the association interaction is 24 times (2300%) stronger than the solvation interaction, so we expect that GE > 0.			
!	Once again, the computed IDACs bear out this intuitive analysis: 15.5 and 5.2.  			
	!Subroutine MEMTest()  ! Comment this line out to make this the main program
	USE GlobConst
    USE Assoc
	IMPLICIT DoublePrecision(A-H,K,O-Z)
    DoublePrecision xFrac(NMX),IDACest(NMX),gamInf(NMX)
    DoublePrecision rLnPhiAssoc(NMX),rLnPhiAssoc1(NMX),rLnPhiAssoc2(NMX),gam(NMX) !,vLiq(NMX)
	common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	!LOUD=.TRUE.
	dumpUnit=6														
	isZiter=0	! This means we compute uAssoc.

	tKelvin=300
	Call Example5Setup(tKelvin,IDACest,iErr) !Set parameters in USEd Assoc and GlobConst(bVol,vLiq). cf. MEM2.f90
	nComps=2
	!Check the uAssoc formula against numerical derivative.
	xFrac(2)=.5d0
	xFrac(1)=1-xFrac(2)
	vMix=SUM( xFrac(1:nComps)*vLiq(1:nComps) )	! vLiq USEd in GlobConst
	call MEM2(isZiter,tKelvin,xFrac,nComps,1/vMix,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr)
	write(*,form610)' Main: XA1, XD1, XA2, XD2, FA,FD ',XA(1,1),XD(1,1),XA(2,1),XD(2,1),FA,FD
	step=1.d-4
	tPlus =tKelvin*(1+step)
	tMinus=tKelvin*(1-step)
	call MEM2(isZiter,tPlus ,xFrac,nComps,1/vMix,zAssoc,aPlus ,uAssoc,rLnPhiAssoc,iErr)
	FAPlus=FA
	FDPlus=FD
	write(*,form610)' Main: XA1+,XD1+,XA2+,XD2+',XA(1,1),XD(1,1),XA(2,1),XD(2,1)
	XA1Plus=XA(1,1)
	XD1Plus=XD(1,1)
	XD2Plus=XD(2,1)
	call MEM2(isZiter,tMinus,xFrac,nComps,1/vMix,zAssoc,aMinus,uAssoc,rLnPhiAssoc,iErr)
	FAMinus=FA
	FDMinus=FD
	write(*,form610)' Main: XA1-,XD1-,XA2-,XD2-',XA(1,1),XD(1,1),XA(2,1),XD(2,1)
	XA1Minus=XA(1,1)
	XD1Minus=XD(1,1)
	XD2Minus=XD(2,1)
	uCalc= -tKelvin*(aPlus-aMinus)/(tPlus-tMinus)
	uCalc= -tKelvin*(aPlus-aMinus)/(tPlus-tMinus)
	uCalc= -tKelvin*(aPlus-aMinus)/(tPlus-tMinus)
	write(*,form610)' Main: uAssoc,uCalc',uAssoc,uCalc
	dXA1= -tKelvin*(XA1Plus-XA1Minus)/(tPlus-tMinus)
	dXD1= -tKelvin*(XD1Plus-XD1Minus)/(tPlus-tMinus)
	dXD2= -tKelvin*(XD2Plus-XD2Minus)/(tPlus-tMinus)
	write(*,form610)' Main: dXA1,dXD1,dXD2',dXA1,dXD1,dXD2
	pause 'Main: Check calculation'
	!Start by evaluating at ~ pure x1. Evaluation at ~ pure x2 occurs as part of first step in loop.
	xInf=1.d-5	! When I made this smaller, Ex 6 had a problem with uAssoc at x1=xInf. JRE 20230206.
	xFrac(2)=xInf
	xFrac(1)=1-xFrac(2)
	call MEM2(isZiter,tKelvin,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1,uAssoc1,rLnPhiAssoc1,iErr)
	if(LOUD)write(dumpUnit,610)'Main: z1,chemPo1=',zAssoc1,rLnPhiAssoc1(1:nComps)
	!tStep=1.d-4
	!tPlus =tKelvin*(1+tStep)
	!tMinus=tKelvin*(1-tStep)
	!call MEM2(isZiter,tPlus,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1p,uAssoc1,rLnPhiAssoc1,iErr)
	!call MEM2(isZiter,tMinus,xFrac,nComps,1/vLiq(1),zAssoc1,aAssoc1,uAssoc1,rLnPhiAssoc1,iErr)
	!uAssoc1= -tKelvin*(aAssoc1p-aAssoc1)/(tPlus-tMinus)
	write(*,'(a)')' xFrac(1),GE_RT,UE_RT,rhoMol_cc,gam(1),gam(2),zAssoc,uAssoc,chemPoAssoc(1),chemPoAssoc(2)'
	do i=0,10,1
		xFrac(1)=i/10.d0
		if(xFrac(1) < xInf)xFrac(1)=xInf
		xFrac(2)=1-xFrac(1)
		if(xFrac(2) < xInf)xFrac(2)=xInf
		bMix=SUM( xFrac(1:nComps)*bVolCC_mol(1:nComps) )
		vMix=SUM( xFrac(1:nComps)*vLiq(1:nComps) )
		rhoMol_cc=1/vMix
		call MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr)
		if(i==0)then
			zAssoc2=zAssoc
			aAssoc2=aAssoc
			uAssoc2=uAssoc
			rLnPhiAssoc2(1:nComps)=rLnPhiAssoc(1:nComps)
			if(LOUD)write(dumpUnit,610)'Main: z2,chemPo2=',zAssoc2,rLnPhiAssoc2(1:nComps)
		endif
		gam(1)=exp( rLnPhiAssoc(1)-zAssoc*bVolCC_mol(1)/bMix-(rLnPhiAssoc1(1)-zAssoc1) )
		gam(2)=exp( rLnPhiAssoc(2)-zAssoc*bVolCC_mol(2)/bMix-(rLnPhiAssoc2(2)-zAssoc2) )
		GE_RT=aAssoc-(xFrac(1)*aAssoc1+xFrac(2)*aAssoc2)
		UE_RT=uAssoc-(xFrac(1)*uAssoc1+xFrac(2)*uAssoc2)
		write(*,'(1x,4f8.4,12f8.3)')xFrac(1),GE_RT,UE_RT,rhoMol_cc,gam(1),gam(2),zAssoc,uAssoc,rLnPhiAssoc(1),rLnPhiAssoc(2)
		if(i== 0)gamInf(1)=gam(1)
		if(i==10)gamInf(2)=gam(2)
		if(i==0.and.LOUD)pause 'Main: halfway there! Check so far.'
	enddo
	write(*,610)'T(K),IDAC1,IDAC2=',tKelvin,gamInf(1:nComps),IDACest(1:nComps)
610	format(1x,a,12E12.4)
	pause 'Success! Check results above.'
	stop
	END

	Subroutine Example1Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	!Example 1. Acetone(1)+chloroform(2) with speadmd model.
	write(dumpUnit,*)'Example 1. Acetone(1)+chloroform(2)'
	isTpt=.TRUE.				!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=3					!acetone
	nDegree(1,3)=1
	nDonors(1,3)=1
	nAcceptors(1,3)=1
	eDonorKcal_mol(1,3)=1.6
	eAcceptorKcal_mol(1,3)=1.0
	bondVolNm3(1,3)=0.0025d0
	bVolCC_mol(1)=36			!cf. ParmsTpt.txt
	vLiq(1)=74					!cf. ParmsPrTcJaubert.txt
	nTypes(2)=4					!chloroform
	nDegree(2,1)=1
	nDonors(2,1)=0
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=1.4
	eDonorKcal_mol(2,1)=0
	bondVolNm3(2,1)=0.0025d0
	bVolCC_mol(2)=38
	vLiq(2)=81

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex1 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,3)*avoNum*( EXP(eAcceptorKcal_mol(1,3)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,3)*avoNum*( EXP(   eDonorKcal_mol(1,3)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1									!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/sqrt(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	FDinf1=FD0inf1/sqrt(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	sumLnXPure2(1)= -nDegree(1,3)*(  nDonors(1,3)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,3)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex1: FA0,FD0,FA,FD,SumLnXiInf1=',FA0inf1,FD0inf1,FAinf1,FDinf1,sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,3)*avoNum*( EXP(eAcceptorKcal_mol(1,3)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,3)*avoNum*( EXP(   eDonorKcal_mol(1,3)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,3)*nDonors(1,3)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,3)*nAcceptors(1,3)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/sqrt(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/sqrt(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,3)*(  nDonors(1,3)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,3)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: FA0,FD0,FA,FD,SumLnXiInf2=',FA0,FD0,FAinf2,FDinf2,sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	IDACest(1:2)=EXP(IDACest(1:2))
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
	if(LOUD)pause 'Ex1 setup done. Check IDACs.'
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example2Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	!Example 2. THF(1)+chloroform(2) with speadmd model.
	write(dumpUnit,*)'Example 2. THF(1)+chloroform(2)'
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!THF
	nDegree(1,1)=1
	nDonors(1,1)=1
	nAcceptors(1,1)=0
	eAcceptorKcal_mol(1,1)=0		!cf.ParmsHB4.txt
	eDonorKcal_mol(1,1)=1.3d0
	bondVolNm3(1,1)=0.0007d0
	bVolCC_mol(1)=.0702793d0*avoNum	!cf. ParmsTpt.txt
	vLiq(1)=59.11d0/0.7214d0		!cf. ParmsPrTcJaubert.txt
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
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralph1inf1,ralph2inf1=',ralphD1inf1,ralphA2inf1
	rhoDelta=ralphD1inf1*ralphA2inf1
	IDACest(1)= -nDegree(1,1)*nDonors(1,1)*LOG(1+nDegree(2,1)*nAcceptors(2,1)*rhoDelta)						   !Eq. 64, Only good for acetone+chloroform, dioxane+chloroform, ...
	Call RdfCalc(rdfContact1,dAlpha,bVolCC_mol(1)/vLiq(1))
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(*,610)'ralph1inf2,ralph2inf2=',ralphD1inf2,ralphA2inf2
	rhoDelta=ralphD1inf2*ralphA2inf2
	IDACest(2)= -nDegree(2,1)*nAcceptors(2,1)*LOG(1+nDegree(1,1)*nDonors(1,1)*rhoDelta)
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest1=',IDACest(1:2)

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex1 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1									!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/sqrt(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	FDinf1=FD0inf1/sqrt(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	sumLnXPure2(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,1)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf1=',sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,1)*nDonors(1,1)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,1)*nAcceptors(1,1)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/sqrt(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/sqrt(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,1)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf2=',sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
	IDACest(1:2)=EXP(IDACest(1:2))
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example3Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	!Example 3. Dioxane(1)+chloroform(2) with speadmd model. cf. EL2ed Example 19.3, p779. 
	write(dumpUnit,*)'Example 3. Dioxane(1)+chloroform(2)'
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!14dioxane...
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

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex1 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1									!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/sqrt(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	FDinf1=FD0inf1/sqrt(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	sumLnXPure2(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,1)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf1=',sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,1)*nDonors(1,1)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,1)*nAcceptors(1,1)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/sqrt(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/sqrt(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,1)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf2=',sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
	IDACest(1:2)=EXP(IDACest(1:2))
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example4Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	!Example 4. 135trioxane(1)+chloroform(2) with speadmd model.
	write(dumpUnit,*)'Example 4. 135trioxane(1)+chloroform(2)'
	nTypes(1)=1						!135trioxane...
	nDegree(1,1)=1
	nAcceptors(1,1)=0
	nDonors(1,1)=3
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

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex1 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1*nDegree(2,1)*nDonors(2,1)									!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/sqrt(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	FDinf1=FD0inf1/sqrt(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	sumLnXPure2(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,1)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf1=',sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,1)*nDonors(1,1)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,1)*nAcceptors(1,1)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/SQRT(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/SQRT(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,1)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: SumLnXiInf2=',sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	IDACest(1:2)=EXP(IDACest(1:2))
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example5Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	!Example 5. methanol(1)+THF(2) with speadmd model.  
	write(dumpUnit,*)'Example 5. methanol(1)+THF(2)'
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!methanol
	nDegree(1,1)=1
	nDonors(1,1)=1
	nAcceptors(1,1)=1
	eAcceptorKcal_mol(1,1)=5
	eDonorKcal_mol(1,1)=5
	bondVolNm3(1,1)=0.0010d0
	bVolCC_mol(1)=.0310936d0*avoNum
	vLiq(1)=32.04d0/0.7896d0		
	nTypes(2)=1						!THF
	nDegree(2,1)=1
	nDonors(2,1)=1
	nAcceptors(2,1)=0
	eAcceptorKcal_mol(2,1)=0		!cf.ParmsHB4.txt
	eDonorKcal_mol(2,1)=1.3d0
	bondVolNm3(2,1)=0.0007d0
	bVolCC_mol(2)=.0702793d0*avoNum	!cf. ParmsTpt.txt
	vLiq(2)=59.11d0/0.7214d0		!cf. ParmsPrTcJaubert.txt

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex5 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1*nDegree(2,1)*nDonors(2,1)		!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/SQRT(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	FDinf1=FD0inf1/SQRT(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	sumLnXPure2(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,1)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex5: FA0,FD0,FA,FD,SumLnXiInf1=',FA0inf1,FD0inf1,FAinf1,FDinf1,sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1yinf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,1)*nDonors(1,1)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,1)*nAcceptors(1,1)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/SQRT(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/SQRT(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,1)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: FA0,FD0,FA,FD,SumLnXiInf2=',FA0,FD0,FAinf2,FDinf2,sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	IDACest(1:2)=EXP(IDACest(1:2))
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
	if(LOUD)pause 'Ex5Setup: done. returning. check IDACs'
610	format(1x,a,12E12.4)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Example6Setup(tKelvin,IDACest,iErr)
	USE GlobConst
	USE Assoc
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision IDACest(NMX),sumLnXPure1(2),sumLnxPure2(2)
	iErr=0
	!Example 6. acetone(1)+methanol(2) with speadmd model.  
	write(dumpUnit,*)'Example 6. acetone(1)+methanol(2)'
	isTpt=.TRUE.					!for rdfCalc, isTPT uses Carnahan-Starling instead of vdw or ESD. 
	nTypes(1)=1						!acetone
	nDegree(1,1)=1
	nDonors(1,1)=1
	nAcceptors(1,1)=1
	eAcceptorKcal_mol(1,1)=1.0
	eDonorKcal_mol(1,1)=1.6
	bondVolNm3(1,1)=0.0025d0
	bVolCC_mol(1)=36				!cf. ParmsTpt.txt
	vLiq(1)=74						!cf. ParmsPrTcJaubert.txt
	nTypes(2)=1						!methanol
	nDegree(2,1)=1
	nDonors(2,1)=1
	nAcceptors(2,1)=1
	eAcceptorKcal_mol(2,1)=5
	eDonorKcal_mol(2,1)=5
	bondVolNm3(2,1)=0.0010d0
	bVolCC_mol(2)=.0310936d0*avoNum
	vLiq(2)=32.04d0/0.7896d0		

	Call RdfCalc(rdfContact2,dAlphaInf1,bVolCC_mol(2)/vLiq(2))
	if(LOUD)write(dumpUnit,610)'Ex5 at inf1: eta,rdfContact,dAlpha=',bVolCC_mol(2)/vLiq(2),rdfContact2,dAlphaInf1
	ralphA1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on acetone
	ralphD1inf1=SQRT(  rdfContact2*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eDonor on acetone
	ralphA2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  ) !eAcceptor on chloroform
	ralphD2inf1=SQRT(  rdfContact2*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(2)  )
	if(LOUD)write(*,610)'ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1=',ralphA1inf1,ralphD1inf1,ralphA2inf1,ralphD2inf1
	FA0inf1=0+1*ralphD2inf1*nDegree(2,1)*nDonors(2,1)		!Eq. 18	  =0 for Ex2.
	FD0inf1=0+1*ralphA2inf1*nDegree(2,1)*nAcceptors(2,1)	!Eq. 19
	FAinf1=FA0inf1/SQRT(1+ralphD1inf1*FD0inf1)						! ~Eq. 61,62	
	FDinf1=FD0inf1/SQRT(1+ralphA1inf1*FA0inf1)						! ~Eq. 61,62
	if(FAinf1>1)FAinf1=1
	if(FDinf1>1)FDinf1=1
	sumLnXPure2(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf1*FDinf1) )+nAcceptors(1,1)*LOG(1+ralphA1inf1*FAinf1)  )
	sumLnXPure2(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf1*FDinf1) )+nAcceptors(2,1)*LOG(1+ralphA2inf1*FAinf1)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex5: FA0,FD0,FA,FD,SumLnXiInf1=',FA0inf1,FD0inf1,FAinf1,FDinf1,sumLnXPure2(1:2)
	Call RdfCalc(rdfContact1,dAlphaInf2,bVolCC_mol(1)/vLiq(1))
	ralphA1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(eAcceptorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD1inf2=SQRT(  rdfContact1*bondVolNm3(1,1)*avoNum*( EXP(   eDonorKcal_mol(1,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphA2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(eAcceptorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	ralphD2inf2=SQRT(  rdfContact1*bondVolNm3(2,1)*avoNum*( EXP(   eDonorKcal_mol(2,1)*1000/(RgasCal*tKelvin))-1 )/vLiq(1)  )
	if(LOUD)write(dumpUnit,610)'ralphA1inf2,ralphD1yinf2,ralphA2inf2,ralphD2inf2=',ralphA1inf2,ralphD1inf2,ralphA2inf2,ralphD2inf2
	FA0=nDegree(1,1)*nDonors(1,1)   *ralphD1inf2+0	!Eq. 18
	FD0=nDegree(1,1)*nAcceptors(1,1)*ralphA1inf2+0	!Eq. 19
	FAinf2=FA0/SQRT(1+ralphD1inf2*FD0)						! ~Eq. 61,62
	FDinf2=FD0/SQRT(1+ralphA1inf2*FA0)						! ~Eq. 61,62
	if(FAinf2>1)FAinf2=1
	if(FDinf2>1)FDinf2=1
	sumLnXPure1(1)= -nDegree(1,1)*(  nDonors(1,1)*LOG( (1+ralphD1inf2*FDinf2) )+nAcceptors(1,1)*LOG(1+ralphA1inf2*FAinf2)  )
	sumLnXPure1(2)= -nDegree(2,1)*(  nDonors(2,1)*LOG( (1+ralphD2inf2*FDinf2) )+nAcceptors(2,1)*LOG(1+ralphA2inf2*FAinf2)  )	!Eq. 65, adapted for THF replacing chloroform. 
	if(LOUD)write(dumpUnit,610)'Ex2: FA0,FD0,FA,FD,SumLnXiInf2=',FA0,FD0,FAinf2,FDinf2,sumLnXPure1(1:2)
	IDACest(1)=sumLnXPure2(1)-sumLnXPure1(1)  -FAinf2*FDinf2+FAinf1*FDinf1*bVolCC_mol(1)/bVolCC_mol(2) 	
	IDACest(2)=sumLnXPure1(2)-sumLnXPure2(2)  -FAinf1*FDinf1+FAinf2*FDinf2*bVolCC_mol(2)/bVolCC_mol(1) 		
	IDACest(1:2)=EXP(IDACest(1:2))
	if(LOUD)write(dumpUnit,610)'Ex2: IDACest3=',IDACest(1:2)
	if(LOUD)pause 'Ex 6 setup done. Check IDACs.'
610	format(1x,a,12E12.4)
	return
	end
