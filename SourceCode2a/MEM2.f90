!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Assoc  ! This module is site-based (similar to Group Contribution (GC) basis). Sums are over sites per molecule then over molecules.
	USE GlobConst, only:nmx,isTPT,isESD,iEosOpt,LOUD,dumpUnit ! nmx is the maximum number of compounds, typically 55. 
	implicit NONE
	integer maxTypes,nsx,maxTypesGlobal,localPool	
	PARAMETER (maxTypes=44,nsx=maxTypes,maxTypesGlobal=999) 
    parameter(localPool=9999) !this must be long enough to cover all SpeadMd site types (e.g. 1401=methanol hydroxy)
	!maxTypes is the max # of site types for all molecules. maxTypes > sum^NC(count(Types))	!nsx = maxTypes (dunno why redundant), nbx=Max Bonds, 
	DoublePrecision eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes)
	DoublePrecision eHbKcal_mol(nmx,maxTypes) ! eHbKcal_mol=(eDonorKcal_mol+eAcceptorKcal_mol)/2
	DoublePrecision bondVolNm3(nmx,maxTypes)
	Integer idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	Integer nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),nDegree(nmx,maxTypes)
	Integer localType(localPool),idLocalType(maxTypes)	!these are to accumulate site lists in Wertheim so the aBipAd,aBipDa arrays don't get too large. cf. AlphaSp
	DoublePrecision aBipAD(maxTypes,maxTypes),aBipDA(maxTypes,maxTypes) !association bips
	DoublePrecision XA(NMX,maxTypes),XD(NMX,maxTypes),XC(NMX,maxTypes),cvAssoc
	DoublePrecision rLnGamRep(nmx),rLnGamAtt(nmx),rLnGamAssoc(nmx) ! These replace common/FugiParts
	LOGICAL LouderWert
	 
	!localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!idLocalType points back to localType for double-linking. E.g. idLocalType(1)=101.
contains
	Subroutine RdfCalc(rdfContact,dAlpha,eta)
	USE GlobConst
	DoublePrecision denom,eta,denom2,void,void2,void4,dLng,rdfContact,d2Lng,d2g,dAlpha,dg_dEta,dAlph_dEta
	Integer iRdfOpt,iErrCode
	! Input:
	! eta = packing fraction
	! iEosOpt, isESD, isTPT cf. GlobConst
	! Output:
	! rdfContact= radial distribution function at contact
	! dAlpha=dLn(alpha)/dLn(rho)
	! Miscellaneous:
	! alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/dEta) = 1+(eta/rdf)*(dRdf/dEta)
	! dLng = (eta/rdf)*(dRdf/dEta)     
	! iRdfOpt	- option for characterizing the radial distribution funciton at contact.	
	!			= 0, not specified => error
	!			= 1, ESD form
	!			= 2, Carnahan-Starling form
	!			= 3, activity coefficient model form (rdf=universal constant=g(eta=0.4)
	!			= 4, ESD form, geometric association rule.
	iRdfOpt=0					!ESD non-geometric form
	if(isESD)then
		iRdfOpt=4	!ESD geometric form
	elseif(isTpt)then
		iRdfOpt=2	!CS form
	elseif(iEosOpt==6)then
		iRdfOpt=3	!activity model
	endif

	if(iRdfOpt.eq.0)then
		iErrCode=11
		if(LOUD)write(dumpUnit,*) 'Error in RdfCalc: iRdfOpt not specified'
		if(LOUD)write(dumpUnit,*) '1=ESD form, 2=CS form, 3=activity coeff form, 4=ESD+geometric association'
		return
	endif
	!alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/dEta) = 1+(eta/rdf)*(dRdf/dEta)
	!dLng = (eta/rdf)*(dRdf/dEta) = dLng/dLnEta
	denom=1.d0-1.9d0*eta 
	denom2=denom*denom
	void=1.d0-eta
	void2=void*void
	void4=void2*void2    

	dLng=1.9d0*eta/denom
	if(iRdfOpt==2)dLng=eta*( 3.d0/void-1.d0/(2.d0-eta) )
	if(iRdfOpt==3)dLng= -1.d0
	rdfContact=1/denom
	if(iRdfOpt==2)rdfContact=(1.d0-eta/2.d0)/void/void2
	if(iRdfOpt==3)rdfContact=1.d0
	d2Lng=1.9d0*eta/denom2
	d2g=d2Lng-dLng
	dAlpha=1+dLng
	dg_deta=1.9d0/denom2
	dAlph_deta=dg_deta
	if(iRdfOpt==2)dg_deta=(2.5d0-eta)/void4
  	if(iRdfOpt==2)dAlph_deta=3.d0/void2-2.d0/((2.d0-eta)*(2.d0-eta))
	return
	END	Subroutine RdfCalc
END MODULE Assoc

	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!  ELLIOTT'S SUBROUTINE 
	!  20221202 Initial adaptation of MEM2 for speadmd. Started as accelerator for initial guesses of SS method. Found SS no longer necessary. 90 minutes->63 seconds for regression of Jaubert database.
	!  20230114 Unified basis for ESD96 and ESD-MEM2.
	SUBROUTINE MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr )
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID,PGLinputDir,dumpUnit,RgasCal
	USE Assoc ! XA,XD,XC,NAcceptors(NMX),NDonors(NMX),NDegree(NMX),...
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (FA,FD) AND properties (zAssoc, lnPhiAssoc,...) given T,rho,x
	!  Reference: Elliott, IECR, 61(42):15724 (2022). doi: 10.1021/acs.iecr.2c02723
	!  INPUT:
	!  isZiter = 1 if only zAssoc,aAssoc are required (for zFactor iterations), 0 for lnPhiAssoc,aAssoc,uAssoc, -1 if derivatives are required, integer
	!  tKelvin = T(K), real*8
	!  xFrac   = component mole fractions, real*8
	!  nComps  = # of components, integer
	!  rhoMol_cc = molar density in gmol/cm3, real*8
	!  OUTPUT:
	!  zAssoc  = association effect on compressibility factor, real*8 
	!  aAssoc  = association effect on residual Helmholtz energy, real*8 
	!  uAssoc  = association effect on residual internal energy, real*8 
	!  rLnPhiAssoc  = association effect on log of fugacity coefficient, real*8
	!  iErr    = error code, integer (see errMsg vector below for code definitions).
	!  Incidentals that may be of interest: 
	!  FA,FD   = THE CHARACTERISTIC ASSOCIATION = (1/XAi-1)/RALPHi, Eqs.7-8
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei*rdfContact, Eq. 5
	!  
	Implicit NONE !DoublePrecision(A-H,K,O-Z)
	DoublePrecision xFrac(NMX),RALPHA(NMX,maxTypes),RALPHD(NMX,maxTypes),rLnPhiAssoc(NMX) !,KVE(NMX)
	DoublePrecision tKelvin,rhoMol_cc,zAssoc,aAssoc,uAssoc 
	DoublePrecision eta,picard,Ftol,bVolMix,rdfContact,dAlpha,ralphAmax,ralphDmax,avgNAS,avgNDS,FA0,FD0,FC0_ralphC
	DoublePrecision FA,betadFD_dBeta,FD,betadFA_dBeta,FC,betadFC_dBeta,epsCC,ralphC,sqArg,FAold,error,errOld,avgFo,delFo,delFA,delFD,sqArgA,sqArgD,sumA,sumD,errA,errD,FDold,FA1,XCtemp,hAD,hDA,hCC,bepsA,bepsD,bepsC,sumLnXi,dLnAlpha_dLnBeta
	DoublePrecision, STATIC:: xOld(NMX),ralphAmean,ralphDmean,etaOld,rdfOld ! static is equivalent to the "SAVE" attribute
	Integer, STATIC:: IDold(NMX)  								            ! static is equivalent to the "SAVE" attribute
	Integer isZiter,nComps,iErr, iComp,isDonor,isAcceptor,moreDonors,iErrMEM1,iType,i,j,itMax,nIter
	LOGICAL LOUDER,isAcid
	Character*123 errMsg(22)
	data etaOld,rdfOld,ralphAmean,ralphDmean/0.3583,3.0,1.0,1.0/ !This makes these values static for Intel Compiler without STATIC. We need to reuse them on next call.
	iErr=0
	errMsg( 1)=' MEM2: ERROR - NO CNVRG. ' ! This is a warning because future calls (like iterating on Z) may converge and remove the issue.
	errMsg( 2)=' MEM2: Failed bond site balance.' ! This is a warning because some (crude) approximations might not care about balance. Some users may choose to elevate this to a critical error.
	errMsg(11)=' MEM2: rdfContact<1'
	errMsg(12)=' MEM2: sqArg for ralphMean update.'
	errMsg(13)=' MEM2: sqArg(A or D)'
	errMsg(14)=' MEM2: XA,XD, or XC < 0' 
	errMsg(15)=' MEM2: etaOld < 0?' 
	errMsg(16)=' MEM2: Sorry. Derivatives beyond zAssoc and uAssoc have not been implemented yet.'
	if(isZiter < 0)then
		iErr=16
		if(Louder)write(dumpUnit,601)TRIM(errMsg(iErr))
		return
	endif 
601	format(1x,a,8E12.4)
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LouderWert
	!LOUDER = .TRUE. ! for local debugging.
	picard=0.83D0
	Ftol=1.D-6
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	Call RdfCalc(rdfContact,dAlpha,eta)	!dAlpha = dLnAlpha/dLnRho = 1+dLn(g)_dLn(eta)
	if(rdfContact < 1)then
		iErr=11
		if(Louder)write(dumpUnit,'(a,3f9.4)')' MEM2: eta,rdfContact(<1)=',eta,rdfContact
		if(Louder)write(dumpUnit,601)TRIM(errMsg(iErr))
		return
	endif
	zAssoc=0
	aAssoc=0
	uAssoc=0
	rLnPhiAssoc(1:nComps)=0
    iErr=0  
610	format(1x,2i5,i12,i8,i11,i8,2E14.4)
	if(LOUDER)write(dumpUnit,601)' MEM2: T,rho,eta,g^c',tKelvin,rhoMol_cc,eta,rdfContact
	!XA=1  !				   ! initializing the entire array is slow.
	!XD=1
	ralphAmax=0
	ralphDmax=0
	isAcid=.FALSE.
	do iComp=1,nComps
		if(nTypes(iComp)==0 .and.LOUDER)write(dumpUnit,*)'MEM2: nTypes=0 for iComp=',iComp
		do iType=1,nTypes(iComp)
			XA(iComp,iType)=1  ! initializing the entire array is slow.
			XD(iComp,iType)=1
			XC(iComp,iType)=1
			if(idType(iComp,iType)==1603)then
				isAcid=.TRUE.
				epsCC=3*( eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(iComp,iType) )/2
				!ralphC is not subscripted because 1603 is the only acid type and always the same.
				ralphC=SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(epsCC/tKelvin/RgasCal*1000)-1.d0 )  ) ! Note the "3*" in the exponent, which makes the C bond "special."
				if(LOUDER)write(dumpUnit,'(a,F10.6,3E12.4)')' MEM2: epsCC,bondVol,ralphC',epsCC,bondVolNm3(iComp,iType),ralphC
			endif
			isAcceptor=0
			if(nAcceptors(iComp,iType) > 0)isAcceptor=1
			!sqArg=rho*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eAcceptorKcal_Mol(iComp,iType)/tKelvin/RgasCal*1000)-1.d0 )
			!if(sqArg < 0)write(dumpUnit,'(a,2i3,f8.2,f8.4)')' MEM2: iComp,iType,T',iComp,iType,tKelvin,eAcceptorKcal_Mol(iComp,iType)
			ralphA(iComp,iType)=isAcceptor*SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eAcceptorKcal_Mol(iComp,iType)/tKelvin/RgasCal*1000)-1.d0 )  )
			if(ralphA(iComp,iType) > ralphAmax)ralphAmax=ralphA(iComp,iType)
			isDonor=0
			if(nDonors(iComp,iType) > 0)isDonor=1
			!sqArg=rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eDonorKcal_Mol(iComp,iType)   /tKelvin/RgasCal*1000)-1.d0 )
			!if(sqArg < 0)write(dumpUnit,'(a,2i3,f8.2,f8.4)')' MEM2: iComp,iType,T,eDonor',iComp,iType,tKelvin,eDonorKcal_Mol(iComp,iType)
			ralphD(iComp,iType)=isDonor*   SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eDonorKcal_Mol(iComp,iType)   /tKelvin/RgasCal*1000)-1.d0 )  )
			if(ralphD(iComp,iType) > ralphDmax)ralphDmax=ralphD(iComp,iType)
			if(LOUDER)write(dumpUnit,'(a,2i3,2F10.4)')' MEM2:iComp,iType,ralphA,ralphD',iComp,iType,ralphA(iComp,iType),ralphD(iComp,iType) 
		enddo
	enddo
	FA0=0
	FD0=0
	FC0_ralphC=0
	avgNAS=0
	avgNDS=0
	do i=1,nComps
		do j=1,nTypes(i)
			FA0=FA0+xFrac(i)*nDegree(i,j)*nDonors(i,j)   *ralphD(i,j)	
			FD0=FD0+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*ralphA(i,j)
			if(idType(i,j)==1603)FC0_ralphC=FC0_ralphC+xFrac(i)*nDegree(i,j) !only one C allowed on type 1603, so nCarboxy(i,j)=1 always.
			avgNDS=avgNDS+xFrac(i)*nDegree(i,j)*nDonors(i,j)	
			avgNAS=avgNAS+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)
		enddo	
	enddo
	moreDonors=0
	if(avgNAS < avgNDS)moreDonors=1
	if(FA0==0 .or. FD0==0)then
		if(LOUDER)write(dumpUnit,601)' MEM2: No assoc. FA0,FD0=',FA0,FD0
		return ! no need to calculate if either no donors or no acceptors
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if( ABS(FA0-FD0) < Ftol)then ! All we need is MEM1 if FA0=FD0
		Call MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,FA,rLnPhiAssoc,iErrMEM1 )!,rLnPhiAssoc,ier)
		if(iErrMEM1 > 0)then
			iErr=12
			if(LOUDER)write(dumpUnit,*)'MEM2: error from MEM1=',iErrMEM1
		endif
		if(.not.isAcid)return
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(LOUDER)write(dumpUnit,'(a,2E12.4,1x,8i6,1x)')' MEM2: x1,x1old,ID,IDold ',xFrac(1),xOld(1),(ID(i),i=1,nComps),(IDold(i),i=1,nComps)
	if(LOUDER)write(dumpUnit,601)' MEM2: etaOld,rdfOld ',etaOld,rdfOld
	if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 ) > Ftol/10.or.etaOld.le.0  )then	!Only use default estimate if compounds or composition have changed.
	!if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. sumxmy2(nComps,xFrac,xOld) > Ftol/10  )then	!Only use default estimate if compounds or composition have changed.
		if(LOUDER)write(dumpUnit,'(a,i6,1x,6E12.4)')' MEM2:New IDs or X ',SUM( ID(1:nComps)-IDold(1:nComps) ),SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 )
		ralphAmean=ralphAmax
		ralphDmean=ralphDmax
		if(moreDonors)then
			ralphDmean=ralphDmean*avgNAS/avgNDS
		else
			ralphAmean=ralphAmean*avgNDS/avgNAS
		endif
	!if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 ) > Ftol/10  )then	!Only use default estimate if compounds or composition have changed.
		if(LOUDER)write(dumpUnit,601)' MEM2: fresh ralphAmean,ralphDmean= ',ralphAmean,ralphDmean
	else !if( SUM(xFrac(1:nComps)) < 0)then
		if(etaOld > 0)then ! adapt old values.to accelerate Z iterations
			sqArg=eta/etaOld*rdfContact/rdfOld
			if(sqArg < zeroTol)then
				iErr=12
				if(LOUDER)write(dumpUnit,'(a,4f10.4)')' MEM2: sqArg for ralphMean update. eta,etaOld,rdfContact,rdfOld',eta,etaOld,rdfContact,rdfOld
			endif
			ralphAmean=ralphAmean*SQRT(sqArg)
			ralphDmean=ralphDmean*SQRT(sqArg)
		if(LOUDER)write(dumpUnit,601)' MEM2: reused&updated ralphAmean,ralphDmean= ',ralphAmean,ralphDmean
		else
			if(LOUDER)write(dumpUnit,601)' MEM2: etaOld < 0???',etaOld
			iErr=15
			return
		endif
	endif
	FAold =0 ! low guess
	errOld= -FA0 ! ie. when guessing FA=FD=0, sumD=FA0. So, FA(guess)-sumD = 0 - FA0 
	FA    =1 ! Take FA0 and FD0 as initial guesses in secant iteration..
	FD    =1 ! 
	NITER=0
	avgFo=(ralphDmean*FD0+ralphAmean*FA0)/2
	if(Louder)write(dumpUnit,'(a,f8.2,11f8.5)')' MEM2:T,eta',tKelvin,eta !,(etaPure(i),i=1,nComps)
	error=1234
	ITMAX=44
	if(LOUDER)write(dumpUnit,'(a,i3,2f8.2,2f10.4,4E12.4)')' iter,ralphAmean,ralphDmean,FA,FD,FA0,FD0  ',nIter,ralphAmean,ralphDmean,FA,FD,FA0,FD0
	do while(ABS(error)>Ftol.and.NITER<ITMAX)
		NITER=NITER+1
		!if(NITER > 33)picard=0.5D0
		delFo=(ralphDmean*FD0-ralphAmean*FA0)
		delFA=1+delFo
		delFD=1-delFo
		sqArgA=delFA*delFA+4*ralphAmean*FA0
		sqArgD=delFD*delFD+4*ralphDmean*FD0
		if(sqArgA < zeroTol .or. sqArgD < zeroTol)then
			if(LOUDER)write(dumpUnit,'(a,i3,4E12.4)')' MEM2: sqArg(A or D). iter,ralphAmean,FA0,ralphDmean,FD0',nIter,ralphAmean,FA0,ralphDmean,FD0
			iErr=13
			return
		endif 
		FA = 2*FA0/( delFA+DSQRT(sqArgA) )
		FD = 2*FD0/( delFD+DSQRT(sqArgD) )	! 
		sumA=0
		sumD=0
		do i=1,nComps !If ralphMeans are correct, then FA,FD are correct and we should get sumA=FD and sumD=FA.
			do j=1,nTypes(i)
				sumA=sumA+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*ralphA(i,j)/( 1+FA*ralphA(i,j) ) !=FD=FD0/(1+ralphAmean*FA)
				sumD=sumD+xFrac(i)*nDegree(i,j)*nDonors(i,j)   *ralphD(i,j)/( 1+FD*ralphD(i,j) ) !=FA=FA0/(1+ralphDmean*FD) 
			enddo	
		enddo
		errA=FA-sumD
		errD=FD-sumA
		error=MAX(ABS(errA),ABS(errD))
		FAold=FA
		FDold=FD
		!Compute new ralphMeans.
		if(LOUDER)write(dumpUnit,'(a,i3,2f8.2,2f10.4,4E12.4)')' iter,ralphAmean,ralphDmean,FA,FD,errA,errD',nIter,ralphAmean,ralphDmean,FA,FD,errA,errD
		ralphDmean= picard*(-1+FA0/sumD)/(FD+1D-9)+(1-picard)*ralphDmean !Using new sumD to compute new ralphMeans.
		if(ralphDmean < 0)ralphDmean=(-1+FA0/sumD)/(FD+1D-9)
		if(ralphDmean < 0)ralphDmean=zeroTol
		ralphAmean= picard*(-1+FD0/sumA)/(FA+1D-9)+(1-picard)*ralphAmean !add 1D-9 to avoid possible zero divide if FD->0
		if(ralphAmean < 0)ralphAmean=(-1+FD0/sumA)/(FA+1D-9)
		if(ralphAmean < 0)ralphAmean=zeroTol
		if(NITER==1)FA1=FA
	enddo !while abs(error) > 1.D-9. 
	if(NITER>ITMAX-1)then
		iErr=1	!warning level
		if(LOUDER)write(dumpUnit,*)'ERROR - NO CNVRG. NITER,FA1,FA=',NITER,FA1,FA
		!GOTO 86
	ENDIF
	if(NITER > itMax-1 .and. LOUDER)write(dumpUnit,'(a,i3,2i5,2F7.4,3f8.2,2f10.4,4E12.4)')' iter,id(1),id(2),eta,xFrac(1),tKelvin,ralphAmean,ralphDmean,FA,FD,errA,errD',nIter,id(1),id(2),xFrac(1),tKelvin,ralphAmean,ralphDmean,FA,FD,errA,errD
	!  fAssoc ITERATION HAS CONCLUDED
	if(isAcid)then
		FC=2*FC0_ralphC*ralphC/( 1+SQRT(1+4*FC0_ralphC*ralphC*ralphC) )
		XCtemp=1/(1+FC*ralphC)
		if(XCtemp < zeroTol)then
			iErr=14
			if(LOUDER)write(dumpUnit,'(a,2f12.4,2F10.4)')' MEM2:ralphC,FC0,FC,XC',ralphC,FC0_ralphC*ralphC,FC,XCtemp 
			if(LOUDER)write(dumpUnit,'(a,E12.4)')' MEM2: XC=',XCtemp
		endif
	endif
	if(FA < zeroTol .or.FD < zeroTol .or.FC < 0 )then
		if(LOUDER)write(dumpUnit,601)' MEM2: FB < 0? FA,FD,FC=',FA,FD,FC
		if(LOUDER)write(dumpUnit,601)' MEM2: ralphMeans=',ralphAmean,ralphDmean
		if(LOUDER)write(dumpUnit,*) 'MEM2: Check this out!'
	endif
	hAD=0
	hDA=0
	hCC=0
	do i=1,nComps  ! XA and XD USE Assoc.
		sumLnXi=0
		do j=1,nTypes(i)
			XA(i,j)= 1/(1+ralphA(i,j)*FA)
			XD(i,j)= 1/(1+ralphD(i,j)*FD)
			hAD=hAD+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*(1-XA(i,j))
			hDA=hDA+xFrac(i)*nDegree(i,j) * nDonors(i,j) *(1-XD(i,j))
			XC(i,j)= 1
			if(idType(i,j)==1603)then
				XC(i,j)=XCtemp
				hCC=hCC+xFrac(i)*nDegree(i,j) * (1-XC(i,j))
			endif
		enddo
		if(LOUDER)write(dumpUnit,'(a,i3,8F10.7)')' MEM2: i,XA(i,j)=',i,(XA(i,j),j=1,nTypes(i))
		if(LOUDER)write(dumpUnit,'(a,i3,8F10.7)')' MEM2: i,XD(i,j)=',i,(XD(i,j),j=1,nTypes(i))
	enddo

	zAssoc= -dAlpha*(hAD+hDA+hCC)/2
	if(LOUDER)write(dumpUnit,'(a,8f10.4)')' MEM2:FA,FD,zAssoc,dAlpha,h^M,zAcid',FA,FD,zAssoc,dAlpha,hAD+hDA+hCC,-dAlpha*hCC/2
	xOld(1:nComps)=xFrac(1:nComps)
	IDold(1:nComps)=ID(1:nComps)
	etaOld=eta
	rdfOld=rdfContact
	if(isZiter==1)return 
	!aAssoc= 0 !already initialized at the top.
	betadFA_dBeta=0
	betadFD_dBeta=0
	betadFC_dBeta=0
	do i=1,nComps  ! XA and XD USE Assoc.
		sumLnXi=0.d0
		do j=1,nTypes(i)
			if( XA(i,j) < zeroTol .or.XD(i,j) < zeroTol .or.XC(i,j) < zeroTol )then
				if(LOUDER)write(dumpUnit,'(a,2i3,3E12.4)')' MEM2: XBij < 0? i,j,XA,D,C=',i,j,XA(i,j),XD(i,j),XC(i,j)
				iErr=14
				return
			endif
			sumLnXi=sumLnXi+nDegree(i,j)*(  nAcceptors(i,j)*DLOG( XA(i,j) )+nDonors(i,j)*DLOG( XD(i,j) )+ DLOG( XC(i,j) )  )	! 
			bepsA=(eAcceptorKcal_mol(i,j)/RgasCal*1000/tKelvin)
			dLnAlpha_dLnBeta=DEXP(bepsA)	! Eqs. 44,45
			if(bepsA > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsA)/(dLnAlpha_dLnBeta-1)
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralphA(i,j)*dLnAlpha_dLnBeta
			bepsD=(eDonorKcal_mol(i,j)/RgasCal*1000/tKelvin)
			dLnAlpha_dLnBeta=DEXP(bepsD)
			if(bepsD > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsD)/(dLnAlpha_dLnBeta-1)
			betadFD_dBeta=betadFD_dBeta+xFrac(i)*nDegree(i,j) * nDonors(i,j) *XD(i,j)*ralphD(i,j)*dLnAlpha_dLnBeta
			if(idType(i,j)==1603)then
				bepsC=3*(eDonorKcal_mol(i,j)+eDonorKcal_mol(i,j))/(2*RgasCal/1000*tKelvin)
				dLnAlpha_dLnBeta=DEXP(bepsC)
				if(bepsC > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsC)/(dLnAlpha_dLnBeta-1)
				betadFC_dBeta=betadFC_dBeta+xFrac(i) *nDegree(i,j) *XC(i,j)*ralphC*dLnAlpha_dLnBeta
			endif
		enddo
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA+hCC)/2*(dAlpha-1)*bVolCC_mol(i)/bVolMix !*( 1+(dAlpha-1)*rhoMol_cc*bVolCC_mol(1:nComps) )	! Eq. 42
	enddo
	aAssoc=  aAssoc+(hAD+hDA+hCC)/2  ! Eqs. 1 & 40
	uAssoc= -( FA*betadFD_dBeta+FD*betadFA_dBeta+FC*betadFC_dBeta )/2 ! Eq. 43
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' MEM2: aAssoc,fugAssoc= ',aAssoc,(rLnPhiAssoc(i),i=1,nComps)
	if(ABS(hAD-hDA) > Ftol)then
		iErr=2 !Warning
		if(LOUDER)write(dumpUnit,*)'MEM2: Failed bond site balance.'
	endif
	if(LOUDER)then
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid)write(dumpUnit,'(a,3E12.4)')' MEM2: ralphC,XC= ',ralphC,XCtemp
	endif
	return
	end	!subroutine MEM2()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!  ELLIOTT'S SUBROUTINE 
	!  20221202 Initial adaptation MEM2 paper: IECR,61(42):15725. Start as accelerator for initial guesses of SS method.
	SUBROUTINE MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc,rLnPhiAssoc,iErr )!,aAssoc,uAssoc,rLnPhiAssoc,ier)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID,dumpUnit,RgasCal
	USE Assoc
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (fAssoc) AND zAssoc given rho,VM
	Implicit NONE !DoublePrecision(A-H,K,O-Z)
	!PARAMETER(NMX=55)
	DoublePrecision xFrac(NMX),RALPH(NMX,maxTypes),rLnPhiAssoc(NMX) !,KVE(NMX)
	DoublePrecision tKelvin,rhoMol_cc,zAssoc,aAssoc,uAssoc 
	DoublePrecision eta,bVolMix,rdfContact,dAlpha,FA0 !,picard,ralphAmax !,Ftol,ralphDmax,avgNAS,avgNDS !,ralphDmean,FD0,FC0_ralphC
	DoublePrecision betadFA_dBeta,sqArg,errOld
	DoublePrecision sumA,hAD,hDA,bepsA,sumLnXi,dLnAlpha_dLnBeta,Fold,fAssoc1,fAssoc,ERR,change,bigYhb
	DoublePrecision, STATIC:: xOld(NMX),ralphMean,etaOld,rdfOld ! static is equivalent to the "SAVE" attribute
	Integer, STATIC:: IDold(NMX)  								            ! static is equivalent to the "SAVE" attribute
	Integer isZiter,nComps,iErr, iComp,iType,i,j,itMax,nIter,nSitesTot,iComplex,iErrCode !,iErrMEM,moreDonors,isDonor,isAcceptor
	LOGICAL LOUDER,isAcid
	!DoublePrecision ralphMatch(NMX)
    !Integer ier(12)	 !NAS(NMX),NDS(NMX),ND(NMX),
	!  NDi     = THE DEGREE OF POLYMERIZATION OF COMPO i
	!  fAssoc  = THE CHARACTERISTIC ASSOCIATION = 1/RALPHi(1/XAi-1)
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei/(1-1.9eta)
	!  ier     = 4  eta > 0.53
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LouderWert
	!LOUDER = .TRUE. ! for local debugging.
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	Call RdfCalc(rdfContact,dAlpha,eta)
	if(rdfContact < 1)then
		if(Louder)write(dumpUnit,'(a,3f9.4)')' MEM1: Returning zero. eta,rdfContact(<1)=',eta,rdfContact
		if(Louder)write(dumpUnit,*)
		return
	endif
	zAssoc=0
	aAssoc=0
	uAssoc=0
    iErr=0  
	!Note: 1st approximation assumes equal donors and acceptors on each site.
	nSitesTot=0
	isAcid=0
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)		
			nSitesTot=nSitesTot+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)+nDonors(iComp,iType) ) !ndHb(iType) should be zero if iType does not hBond.
			if(idType(iComp,iType)==1603)isAcid=1
			iComplex=nDonors(iComp,iType)*nAcceptors(iComp,iType)
			eHbKcal_mol(iComp,iType)=( eDonorKcal_mol(iComp,iType)+eAcceptorKcal_mol(iComp,iType) )/2.d0 
			if(iComplex.ne.0)iComplex=1  !in general, nDs or nAs could be 2, 3,... (e.g. water).  But complexation is either true or false.
			bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/RgasCal*1000)-1) !the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.

			SQARG=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
			IF(SQARG < 0)THEN
				if(Louder)write(dumpUnit,*) 'neg alpha in wertheim. ETA=',ETA
				iErrCode=3
				if(Louder)write(dumpUnit,*)
				RETURN
			ENDIF
			ralph(iComp,iType)=DSQRT( SQARG )
		enddo
	enddo

	!BEGIN SECANT ITERATION TO DETERMINE fAssoc	assuming symmetric case.
	fAssoc=0            
	FA0=0  ! cf. IECR 61(42):15724 Eq. 18
	ralphMean=0
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)      
			iComplex=nAcceptors(iComp,iType)*nDonors(iComp,iType)
			if(iComplex.ne.0)then !only apply symmetric rule when acceptors AND donors on a site.
				FA0=FA0+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType) !/(1+fAssoc*ralph(iComp,iType))
				if(LOUDER)write(dumpUnit,'(a,2i3,8E12.4)')' MEM1: iComp,iType,ralph',iComp,iType,ralph(iComp,iType)
				if(ralph(iComp,iType) > ralphMean)ralphMean=ralph(iComp,iType) ! Infinity norm for initial estimate.
			endif
		enddo
	enddo
	errOld=fAssoc-FA0
	fOld=fAssoc
	IF(ABS(errOld) < 1.D-9.and.LOUDER)write(dumpUnit,601)' MEM1: No assoc? FA0=',FA0	!don't iterate if errOld < 1D-9
	IF(ABS(errOld) > 1.D-9)then	!don't iterate if errOld < 1D-9
		fAssoc1=2*FA0/( 1+DSQRT(1+4*ralphMean*FA0) )  ! Eq.23 where FD0=FA0. This is exact for binary like benzene+methanol.
		fAssoc=fAssoc1 
		NITER=0
		ERR=1111
		itMax=123
		do while(ABS(ERR) > 1.D-9.AND.NITER < itMax) ! Don't use ERR/fAssoc because:(1)fAssoc~0 sometimes (2) max(fAssoc) ~ 1 which isn't "large."
			NITER=NITER+1                  
			SUMA=0
			DO iComp=1,nComps
				do iType=1,nTypes(iComp)      
					iComplex=nAcceptors(iComp,iType)*nDonors(iComp,iType)
					if(iComplex > 0)then !only apply symmetric rule when acceptors AND donors on a site.
						SUMA=SUMA+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType)/(1+fAssoc*ralph(iComp,iType))
					endif
				enddo
			enddo
			ERR=fAssoc-SUMA
			!checkSum=FA0/(1+fAssoc*ralphMean)
			if( ABS((ERR-errOld)/errOld) < zeroTol .and. Louder)write(dumpUnit,601)'MEM1: err~errOld=',err,errOld
			CHANGE=ERR/(ERR-errOld)*(fAssoc-fOld)
			if(LOUDER)write(dumpUnit,601)' MEM1: FA,SUMA',fAssoc,SUMA 
			fOld=fAssoc
			errOld=ERR
			fAssoc=fAssoc-CHANGE
			!ralphMean=(-1+FA0/sumA)/(fAssoc+1D-9)
		ENDDO
		IF(NITER > itMax-1)THEN
			if(LouderWert)write(dumpUnit,*)'ERROR - NO CNVRG ON fAssoc'
			iErrCode=4
			if(Louder)write(dumpUnit,*)
			return
		ENDIF
	endif

	!fAssoc ITERATION HAS CONCLUDED
	if(fAssoc < zeroTol)then
		if(LOUDER)write(dumpUnit,601)' MEM1: 0 > FA=',fAssoc
	endif
601 format(1x,a,8E12.4)

	!Elliott'96 Eq.35 (mod for CS or ESD rdf)
	ZASSOC= -fAssoc*fAssoc*dAlpha
	if(LOUDER)write(dumpUnit,601)'MEM1: eta,fAssoc,Zassoc=',eta,fAssoc,Zassoc
	if(isZiter==1)return
	!XA=1 !set to unity means no association !Setting the entire array to 1 is slow.	          
	!XD=1 !set to unity means no association	          
	!XC=1 !set to unity means no association	          
	hAD=0
	do i=1,nComps  ! XA and XD USE Assoc.
		do j=1,nTypes(i)
			XA(i,j)= 1/(1+ralph(i,j)*fAssoc)
			XD(i,j)=XA(i,j)
			XC(i,j)=1
			hAD=hAD+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*(1-XA(i,j))
		enddo
		if(LOUDER)write(dumpUnit,'(a,i3,4F10.7)')' MEM1: i,XA(i,j)=',i,(XA(i,j),j=1,nTypes(i))
	enddo

	if(LOUDER)write(dumpUnit,'(a,8f10.4)')' MEM1:fAssoc,zAssoc,dAlpha,h^M,sumLnXi',fAssoc,zAssoc,dAlpha,hAD+hDA,sumLnXi
	hDA=hAD
	zAssoc= -dAlpha*(hAD+hDA)/2 ! check same?
	!aAssoc= 0 	!already initialized
	do i=1,nComps  ! XA and XD USE Assoc.
		sumLnXi=0.d0
		do j=1,nTypes(i)
			if( XA(i,j) < zeroTol)then
				if(LOUDER)write(dumpUnit,'(a,2i3,8E12.4)')' MEM1: XAij < 0? i,j,XA=',i,j,XA(i,j)
				iErr=14
				return
			endif
			sumLnXi=sumLnXi+nDegree(i,j)*(  nAcceptors(i,j)*DLOG( XA(i,j) )*2  )	! ln(XD)=ln(XA) so *2.
		enddo
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA)/2*(dAlpha-1)*bVolCC_mol(i)/bVolMix !*( 1+(dAlpha-1)*rhoMol_cc*bVolCC_mol(1:nComps) )	! Eq. 42
		!write(dumpUnit,'(a,i3,6E12.4)')' MEM1: iComp,sumLnXi,lnPhi=',i,sumLnXi,rLnPhiAssoc(i)
	enddo
	aAssoc=  aAssoc+(hAD+hDA)/2  ! Eqs. 1 & 40
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' MEM1:         fugAssocBefore=',(rLnPhiAssoc(i),i=1,nComps)
	
	betadFA_dBeta=0	!cf E'96(19)
	do i=1,nComps !cf E'96(19)	 TODO: rewrite this in terms of ralph.
		do j=1,nTypes(iComp)
			bepsA=(eAcceptorKcal_mol(i,j)/RgasCal*1000/tKelvin)
			dLnAlpha_dLnBeta=DEXP(bepsA)	! Eqs. 44,45
			if(bepsA > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsA)/(dLnAlpha_dLnBeta-1)
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralph(i,j)*dLnAlpha_dLnBeta
		enddo
	enddo
	! Note: betadFA_dBeta here is actually betadFA_dBeta/(-half)
	uAssoc= -fAssoc*betadFA_dBeta   !cf E'22(Eq43): FA=FD=> -(FA*betadFD_dBeta+FD*betadFA_dBeta)/2 = -FA*betadFA_dBeta
	xOld(1:nComps)=xFrac(1:nComps)
	IDold(1:nComps)=ID(1:nComps)
	etaOld=eta
	rdfOld=rdfContact
	if(LOUDER)write(dumpUnit,'(a,I4,2F10.6,2E12.4)')' MEM1:F1,fAssoc,zAssoc,rmsErr,nIter=',NITER,fAssoc1,fAssoc,zAssoc
	return
	end	!subroutine MEM1()

