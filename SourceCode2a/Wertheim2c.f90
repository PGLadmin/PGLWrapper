!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Assoc  ! This module is site-based (similar to Group Contribution (GC) basis). Sums are over sites per molecule then over molecules.
	USE GlobConst, only:nmx,isTPT,isESD,iEosOpt,LOUD ! nmx is the maximum number of compounds, typically 55. 
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
	LOGICAL LouderWert 
	!localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!idLocalType points back to localType for double-linking. E.g. idLocalType(1)=101.
contains
	Subroutine RdfCalc(rdfContact,dAlpha,eta)
	USE GlobConst
	DoublePrecision denom,eta,denom2,void,void2,void4,dLng,rdfContact,d2Lng,d2g,dAlpha,dg_dEta,dAlph_dEta
	Integer iRdfOpt,iErrCode
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
		if(LOUD)write(*,*) 'Error in RdfCalc: iRdfOpt not specified'
		if(LOUD)pause '1=ESD form, 2=CS form, 3=activity coeff form, 4=ESD+geometric association'
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
	!  20221202 Initial adaptation MEM2 for speadmd: IECR,61(42):15725. Started as accelerator for initial guesses of SS method. Found SS no longer necessary. 90 minutes->63 seconds for regression of Jaubert database.
	SUBROUTINE MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr )!,ier)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	!USE SpeadParms
	USE Assoc ! XA,XD,XC,NAcceptors(NMX),NDonors(NMX),NDegree(NMX),...
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (FA,FD) AND properties (zAssoc, lnPhiAssoc,...) given T,rho,x
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision xFrac(NMX),RALPHA(NMX,maxTypes),RALPHD(NMX,maxTypes),xOld(NMX),rLnPhiAssoc(NMX) !,KVE(NMX)
	Integer IDold(NMX)
	!DoublePrecision ralphMatch(NMX)
    
	LOGICAL LOUDER
	Character*77 errMsg(22)
	!  Reference: Elliott, IECR, 61:15724 (2022).
	!  INPUT:
	!  isZiter = 1 if only zAssoc is required (for zFactor iterations), 0 for lnPhiAssoc,aAssoc,uAssoc, -1 if derivatives are required, integer
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
	iErr=0
	errMsg( 1)=' ERROR - NO CNVRG. '
	errMsg( 2)=' MEM2: Failed bond site balance.'
	errMsg(11)=' MEM2: rdfContact<1'
	errMsg(12)=' MEM2: sqArg for ralphMean update.'
	errMsg(13)=' MEM2: sqArg(A or D)'
	errMsg(14)=' MEM2: XA,XD, or XC < 0' 
	errMsg(15)=' MEM2: etaOld < 0?' 
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LouderWert
	!LOUDER = .TRUE. ! for local debugging.
	picard=0.83D0
	Ftol=1.D-6
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	voidFrac=1-eta
	voidFrac3=voidFrac*voidFrac*voidFrac
	!rdfContact=(1-eta/2)/voidFrac3
	!dAlpha=1+eta*( 3/voidFrac-1/(2-eta) )  !dAlpha=(rho/alpha)*dAlpha/dRho
	Call RdfCalc(rdfContact,dAlpha,eta)
	if(rdfContact < 1)then
		iErr=11
		if(Louder)write(*,'(a,3f9.4)')' MEM2: eta,rdfContact(<1)=',eta,rdfContact
		if(Louder)pause
		return
	endif
	zAssoc=0
	aAssoc=0
	uAssoc=0
	rLnPhiAssoc(1:nComps)=0
    iErr=0  
	if(LOUDER)write(*,601)' MEM2: T,rho,eta,g^c',tKelvin,rhoMol_cc,eta,rdfContact
	!XA=1  !				   ! initializing the entire array is slow.
	!XD=1
	ralphAmax=0
	ralphDmax=0
	isAcid=0
	do iComp=1,nComps
		if(nTypes(iComp)==0 .and.LOUDER)print*,'MEM2: nTypes=0 for iComp=',iComp
		do iType=1,nTypes(iComp)
			XA(iComp,iType)=1  ! initializing the entire array is slow.
			XD(iComp,iType)=1
			XC(iComp,iType)=1
			if(idType(iComp,iType)==1603)then
				isAcid=1
				epsCC=3*( eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(iComp,iType) )/2
				!ralphC is not subscripted because 1603 is the only acid type and always the same.
				ralphC=SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(epsCC/tKelvin/1.987D-3)-1.d0 )  ) ! Note the "3*" in the exponent, which makes the C bond "special."
				if(LOUDER)write(*,'(a,F10.6,3E12.4)')' MEM2: epsCC,bondVol,ralphC',epsCC,bondVolNm3(iComp,iType),ralphC
			endif
			isAcceptor=0
			if(nAcceptors(iComp,iType) > 0)isAcceptor=1
			!sqArg=rho*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eAcceptorKcal_Mol(iComp,iType)/tKelvin/1.987D-3)-1.d0 )
			!if(sqArg < 0)write(*,'(a,2i3,f8.2,f8.4)')' MEM2: iComp,iType,T',iComp,iType,tKelvin,eAcceptorKcal_Mol(iComp,iType)
			ralphA(iComp,iType)=isAcceptor*SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eAcceptorKcal_Mol(iComp,iType)/tKelvin/1.987D-3)-1.d0 )  )
			if(ralphA(iComp,iType) > ralphAmax)ralphAmax=ralphA(iComp,iType)
			isDonor=0
			if(nDonors(iComp,iType) > 0)isDonor=1
			!sqArg=rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eDonorKcal_Mol(iComp,iType)   /tKelvin/1.987D-3)-1.d0 )
			!if(sqArg < 0)write(*,'(a,2i3,f8.2,f8.4)')' MEM2: iComp,iType,T,eDonor',iComp,iType,tKelvin,eDonorKcal_Mol(iComp,iType)
			ralphD(iComp,iType)=isDonor*   SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(eDonorKcal_Mol(iComp,iType)   /tKelvin/1.987D-3)-1.d0 )  )
			if(ralphD(iComp,iType) > ralphDmax)ralphDmax=ralphD(iComp,iType)
			if(LOUDER)write(*,'(a,2i3,2F10.4)')' MEM2:iComp,iType,ralphA,ralphD',iComp,iType,ralphA(iComp,iType),ralphD(iComp,iType) 
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
	if(FA0==0 .or. FD0==0)return ! no need to calculate if either no donors or no acceptors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if( ABS(FA0-FD0) < Ftol)then ! All we need is MEM1 if FA0=FD0
		Call MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,FA,rLnPhiAssoc,iErrMEM1 )!,rLnPhiAssoc,ier)
		if(iErrMEM1 > 0)iErr=12
		if(LOUDER)print*,'MEM2: error from MEM1=',iErrMEM1
		if(isAcid==0)return
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 ) > Ftol/10  )then	!Only use default estimate if compounds or composition have changed.
	!if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. sumxmy2(nComps,xFrac,xOld) > Ftol/10  )then	!Only use default estimate if compounds or composition have changed.
		if(LOUDER)write(*,'(a,i4,6E12.4)')' MEM2:New IDs or X',SUM( ID(1:nComps)-IDold(1:nComps) ),SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 )
		ralphAmean=ralphAmax
		ralphDmean=ralphDmax
		if(moreDonors)then
			ralphDmean=ralphDmean*avgNAS/avgNDS
		else
			ralphAmean=ralphAmean*avgNDS/avgNAS
		endif
	else
		if(etaOld.ge.0)then ! adapt old values.to accelerate Z iterations
			sqArg=eta/etaOld*rdfContact/rdfOld
			if(sqArg < zeroTol)then
				iErr=12
				if(LOUDER)write(*,'(a,4f10.4)')' MEM2: sqArg for ralphMean update. eta,etaOld,rdfContact,rdfOld',eta,etaOld,rdfContact,rdfOld
			endif
			ralphAmean=ralphAmean*SQRT(sqArg)
			ralphDmean=ralphDmean*SQRT(sqArg)
		else
			if(LOUDER)write(*,601)' MEM2: etaOld < 0???',etaOld
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
	if(Louder)write(*,'(a,f8.2,11f8.5)')' MEM2:T,eta',tKelvin,eta !,(etaPure(i),i=1,nComps)
	error=1234
	ITMAX=44
	if(LOUDER)write(*,'(a,i3,2f8.2,2f10.4,4E12.4)')' iter,ralphAmean,ralphDmean,FA,FD,FA0,FD0  ',nIter,ralphAmean,ralphDmean,FA,FD,FA0,FD0
	do while(ABS(error)>Ftol.and.NITER<ITMAX)
		NITER=NITER+1
		!if(NITER > 33)picard=0.5D0
		delFo=(ralphDmean*FD0-ralphAmean*FA0)
		delFA=1+delFo
		delFD=1-delFo
		sqArgA=delFA*delFA+4*ralphAmean*FA0
		sqArgD=delFD*delFD+4*ralphDmean*FD0
		if(sqArgA < zeroTol .or. sqArgD < zeroTol)then
			if(LOUDER)write(*,'(a,i3,4E12.4)')' MEM2: sqArg(A or D). iter,ralphAmean,FA0,ralphDmean,FD0',nIter,ralphAmean,FA0,ralphDmean,FD0
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
		if(LOUDER)write(*,'(a,i3,2f8.2,2f10.4,4E12.4)')' iter,ralphAmean,ralphDmean,FA,FD,errA,errD',nIter,ralphAmean,ralphDmean,FA,FD,errA,errD
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
		if(LOUDER)write(*,*)'ERROR - NO CNVRG. NITER,FA1,FA=',NITER,FA1,FA
		!GOTO 86
	ENDIF
	if(NITER > itMax-1 .and. LOUDER)write(*,'(a,i3,2i5,2F7.4,3f8.2,2f10.4,4E12.4)')' iter,id(1),id(2),eta,xFrac(1),tKelvin,ralphAmean,ralphDmean,FA,FD,errA,errD',nIter,id(1),id(2),xFrac(1),tKelvin,ralphAmean,ralphDmean,FA,FD,errA,errD
	!  fAssoc ITERATION HAS CONCLUDED
	if(isAcid==1)then
		FC=2*FC0_ralphC*ralphC/( 1+SQRT(1+4*FC0_ralphC*ralphC*ralphC) )
		XCtemp=1/(1+FC*ralphC)
		if(XCtemp < zeroTol)then
			iErr=14
			if(LOUDER)write(*,'(a,2f12.4,2F10.4)')' MEM2:ralphC,FC0,FC,XC',ralphC,FC0_ralphC*ralphC,FC,XCtemp 
			if(LOUDER)write(*,'(a,E12.4)')' MEM2: XC=',XCtemp
		endif
	endif
	if(FA < zeroTol .or.FD < zeroTol .or.FC < 0 )then
		if(LOUDER)write(*,601)' MEM2: FB < 0? FA,FD,FC=',FA,FD,FC
		if(LOUDER)write(*,601)' MEM2: ralphMeans=',ralphAmean,ralphDmean
		if(LOUDER)pause 'MEM2: Check this out!'
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
		if(LOUDER)write(*,'(a,i3,8F10.7)')' MEM2: i,XA(i,j)=',i,(XA(i,j),j=1,nTypes(i))
		if(LOUDER)write(*,'(a,i3,8F10.7)')' MEM2: i,XD(i,j)=',i,(XD(i,j),j=1,nTypes(i))
	enddo
601	format(1x,a,8E12.4)

	zAssoc= -dAlpha*(hAD+hDA+hCC)/2
	if(LOUDER)write(*,'(a,8f10.4)')' MEM2:FA,FD,zAssoc,dAlpha,h^M,zAcid',FA,FD,zAssoc,dAlpha,hAD+hDA+hCC,-dAlpha*hCC/2
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
				if(LOUDER)write(*,'(a,2i3,3E12.4)')' MEM2: XBij < 0? i,j,XA,D,C=',i,j,XA(i,j),XD(i,j),XC(i,j)
				iErr=14
				return
			endif
			sumLnXi=sumLnXi+nDegree(i,j)*(  nAcceptors(i,j)*DLOG( XA(i,j) )+nDonors(i,j)*DLOG( XD(i,j) )+ DLOG( XC(i,j) )  )	! 
			bepsA=(eAcceptorKcal_mol(i,j)/1.987D-3/tKelvin)
			beta_alphadAlpha_dBeta=DEXP(bepsA)	! Eqs. 44,45
			if(bepsA > 1.D-4)beta_alphadAlpha_dBeta=beta_alphadAlpha_dBeta*(bepsA)/(beta_alphadAlpha_dBeta-1)
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralphA(i,j)*beta_alphadAlpha_dBeta
			bepsD=(eDonorKcal_mol(i,j)/1.987D-3/tKelvin)
			beta_alphadAlpha_dBeta=DEXP(bepsD)
			if(bepsD > 1.D-4)beta_alphadAlpha_dBeta=beta_alphadAlpha_dBeta*(bepsD)/(beta_alphadAlpha_dBeta-1)
			betadFD_dBeta=betadFD_dBeta+xFrac(i)*nDegree(i,j) * nDonors(i,j) *XD(i,j)*ralphD(i,j)*beta_alphadAlpha_dBeta
			if(idType(i,j)==1603)then
				bepsC=3*(eDonorKcal_mol(i,j)+eDonorKcal_mol(i,j))/(2*1.987D-3*tKelvin)
				beta_alphadAlpha_dBeta=DEXP(bepsC)
				if(bepsC > 1.D-4)beta_alphadAlpha_dBeta=beta_alphadAlpha_dBeta*(bepsC)/(beta_alphadAlpha_dBeta-1)
				betadFC_dBeta=betadFC_dBeta+xFrac(i) *nDegree(i,j) *XC(i,j)*ralphC*beta_alphadAlpha_dBeta
			endif
		enddo
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA+hCC)/2*(dAlpha-1)*bVolCC_mol(i)/bVolMix !*( 1+(dAlpha-1)*rhoMol_cc*bVolCC_mol(1:nComps) )	! Eq. 42
	enddo
	aAssoc=  aAssoc+(hAD+hDA+hCC)/2  ! Eqs. 1 & 40
	uAssoc= -( FA*betadFD_dBeta+FD*betadFA_dBeta+FC*betadFC_dBeta )/2 ! Eq. 43
	if(LOUDER)write(*,'(a,6E12.4)')' MEM2: aAssoc,fugAssoc= ',aAssoc,(rLnPhiAssoc(i),i=1,nComps)
	if(ABS(hAD-hDA) > Ftol)then
		iErr=2 !Warning
		if(LOUDER)print*,'MEM2: Failed bond site balance.'
	endif
	if(LOUDER)then
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid==1)write(*,'(a,3E12.4)')' MEM2: ralphC,XC= ',ralphC,XCtemp
	endif
	return
	end	!subroutine MEM2()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!  ELLIOTT'S SUBROUTINE 
	!  20221202 Initial adaptation MEM2 paper: IECR,61(42):15725. Start as accelerator for initial guesses of SS method.
	SUBROUTINE MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc,rLnPhiAssoc,iErr )!,aAssoc,uAssoc,rLnPhiAssoc,ier)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	USE Assoc
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (fAssoc) AND zAssoc given rho,VM
	Implicit DoublePrecision(A-H,K,O-Z)
	!PARAMETER(NMX=55)
	DoublePrecision xFrac(NMX),RALPH(NMX,maxTypes),xOld(NMX),rLnPhiAssoc(NMX) !,KVE(NMX)
	Integer IDold(NMX)
	!DoublePrecision ralphMatch(NMX)
    !Integer ier(12)	 !NAS(NMX),NDS(NMX),ND(NMX),
	LOGICAL LOUDER 
	!  NDi     = THE DEGREE OF POLYMERIZATION OF COMPO i
	!  fAssoc  = THE CHARACTERISTIC ASSOCIATION = 1/RALPHi(1/XAi-1)
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei/(1-1.9eta)
	!  ier     = 4  eta > 0.53
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LouderWert
	!LOUDER = .TRUE. ! for local debugging.
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	voidFrac=1-eta
	voidFrac3=voidFrac*voidFrac*voidFrac
	Call RdfCalc(rdfContact,dAlpha,eta)
	if(rdfContact < 1)then
		if(Louder)write(*,'(a,3f9.4)')' MEM1: Returning zero. eta,rdfContact(<1)=',eta,rdfContact
		if(Louder)pause
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
			bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/1.987D-3)-1) !the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.

			SQARG=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
			IF(SQARG < 0)THEN
				if(Louder)write(*,*) 'neg alpha in wertheim. ETA=',ETA
				iErrCode=3
				if(Louder)pause
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
				if(LOUDER)write(*,'(a,2i3,8E12.4)')' MEM1: iComp,iType,ralph',iComp,iType,ralph(iComp,iType)
				if(ralph(iComp,iType) > ralphMean)ralphMean=ralph(iComp,iType) ! Infinity norm for initial estimate.
			endif
		enddo
	enddo
	errOld=fAssoc-FA0
	fOld=fAssoc
	IF(ABS(errOld) < 1.D-9.and.LOUDER)write(*,601)' MEM1: No assoc? FA0=',FA0	!don't iterate if errOld < 1D-9
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
			if( ABS((ERR-errOld)/errOld) < zeroTol .and. Louder)write(*,601)'MEM1: err~errOld=',err,errOld
			CHANGE=ERR/(ERR-errOld)*(fAssoc-fOld)
			if(LOUDER)write(*,601)' MEM1: FA,SUMA',fAssoc,SUMA 
			fOld=fAssoc
			errOld=ERR
			fAssoc=fAssoc-CHANGE
			!ralphMean=(-1+FA0/sumA)/(fAssoc+1D-9)
		ENDDO
		IF(NITER > itMax-1)THEN
			if(LouderWert)write(*,*)'ERROR - NO CNVRG ON fAssoc'
			iErrCode=4
			if(Louder)pause
			return
		ENDIF
	endif

	!fAssoc ITERATION HAS CONCLUDED
	if(fAssoc < zeroTol)then
		if(LOUDER)write(*,601)' MEM1: 0 > FA=',fAssoc
	endif
601 format(1x,a,8E12.4)

	!Elliott'96 Eq.35 (mod for CS or ESD rdf)
	ZASSOC= -fAssoc*fAssoc*dAlpha
	if(LOUDER)write(*,601)'MEM1: eta,fAssoc,Zassoc=',eta,fAssoc,Zassoc
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
		if(LOUDER)write(*,'(a,i3,4F10.7)')' MEM1: i,XA(i,j)=',i,(XA(i,j),j=1,nTypes(i))
	enddo

	if(LOUDER)write(*,'(a,8f10.4)')' MEM1:fAssoc,zAssoc,dAlpha,h^M,sumLnXi',fAssoc,zAssoc,dAlpha,hAD+hDA,sumLnXi
	hDA=hAD
	zAssoc= -dAlpha*(hAD+hDA)/2 ! check same?
	!aAssoc= 0 	!already initialized
	do i=1,nComps  ! XA and XD USE Assoc.
		sumLnXi=0.d0
		do j=1,nTypes(i)
			if( XA(i,j) < zeroTol)then
				if(LOUDER)write(*,'(a,2i3,8E12.4)')' MEM1: XAij < 0? i,j,XA=',i,j,XA(i,j)
				iErr=14
				return
			endif
			sumLnXi=sumLnXi+nDegree(i,j)*(  nAcceptors(i,j)*DLOG( XA(i,j) )*2  )	! ln(XD)=ln(XA) so *2.
		enddo
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA)/2*(dAlpha-1)*bVolCC_mol(i)/bVolMix !*( 1+(dAlpha-1)*rhoMol_cc*bVolCC_mol(1:nComps) )	! Eq. 42
		!write(*,'(a,i3,6E12.4)')' MEM1: iComp,sumLnXi,lnPhi=',i,sumLnXi,rLnPhiAssoc(i)
	enddo
	aAssoc=  aAssoc+(hAD+hDA)/2  ! Eqs. 1 & 40
	if(LOUDER)write(*,'(a,6E12.4)')' MEM1:         fugAssocBefore=',(rLnPhiAssoc(i),i=1,nComps)
	
	betadFA_dBeta=0	!cf E'96(19)
	do i=1,nComps !cf E'96(19)	 TODO: rewrite this in terms of ralph.
		do j=1,nTypes(iComp)
			bepsA=(eAcceptorKcal_mol(i,j)/1.987D-3/tKelvin)
			beta_alphadAlpha_dBeta=DEXP(bepsA)	! Eqs. 44,45
			if(bepsA > 1.D-4)beta_alphadAlpha_dBeta=beta_alphadAlpha_dBeta*(bepsA)/(beta_alphadAlpha_dBeta-1)
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralph(i,j)*beta_alphadAlpha_dBeta
		enddo
	enddo
	! Note: betadFA_dBeta here is actually betadFA_dBeta/(-half)
	uAssoc= -fAssoc*betadFA_dBeta   !cf E'22(Eq43): FA=FD=> -(FA*betadFD_dBeta+FD*betadFA_dBeta)/2 = -FA*betadFA_dBeta
	xOld(1:nComps)=xFrac(1:nComps)
	IDold(1:nComps)=ID(1:nComps)
	etaOld=eta
	rdfOld=rdfContact
	if(LOUDER)write(*,'(a,I4,2F10.6,2E12.4)')' MEM1:F1,fAssoc,zAssoc,rmsErr,nIter=',NITER,fAssoc1,fAssoc,zAssoc
	return
	end	!subroutine MEM1()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
	! THE WERTHEIM ROUTINES OF SPEAD AND ESD TOGETHER COMBINED BY AGF, Oct. 2009
	subroutine Wertheim( isZiter,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErrCode)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID ,etaMax
	USE Assoc !XA,XD,XC
	!USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	!PARAMETER(NMX=55,maxTypes=44,maxTypesGlobal=999)
	!USE GlobConst
	character*77 errMsg(11)
    logical HyBond,bAssoc(NMX)
	DoublePrecision ralph(NMX,maxTypes),xFrac(NMX),rLnPhiAssoc(NMX),vMolecNm3(NMX) !bondVolNm3(NMX),
	DoublePrecision xaOld(nmx,maxTypes),xdOld(nmx,maxTypes) !,alphAD(nmx,maxTypes),alphDA(nmx,maxTypes)
	!dimension FAsolv(nmx,maxTypes),FDsolv(nmx,maxTypes),deltaAlphaA(nmx,maxTypes),deltaAlphaD(nmx,maxTypes) ! solvation factors: e.g. FAsolv = sum(xj*Ndj*XDj*alphADij) = (-1+1/XAi) ~ sqrt(alphaii)*Fassoc
	data initCall/1/
	LouderWert=LOUD
	!LouderWert=.TRUE.

	!COMMON/ALPHA/alphAD(NMX,NMX),alphDA(NMX,NMX)
	! NDi		= THE DEGREE OF POLYMERIZATION OF COMPO i
	! fAssoc	= THE CHARACTERISTIC ASSOCIATION = 1/ralphi(1/XAi-1)
	! ralphi	= ROOT OF ALPHA WHERE ALPHAi=RHO*VXi*KADi*Ei/(1-1.9ETA)
	! iErrCode
	!	0		= no error
	!	1		= no convergence on ss iteration.
	errMsg(11)= 'Wertheim Error: this code has not been written to handle asymmetric bonding types in a single molecule.'
	errMsg(1)= 'Wertheim Error: Nonsense in input parameters.'
	vMolecNM3(1:nComps)=bVolCC_mol(1:nComps)/avoNum


	!E.g.1 For water, cf. Nezbeda, I; Pavlicek, J.; FPE, 116:530-536 (1996).
	!sigma=0.2949, Kad=0.003677, eokHb=1293.3, bVol=avoNum*pi*sigma^3/6=8.0868, tK=300, eta=0.5
	! Yhb=73.515, ralph=7.6651, fAssoc=1.3505, XA=0.088094, aAssoc=-7.8936, dAlpha=3.6667, zAssoc=-6.6873
	!E.g.2 For water, 
	!sigma=0.306, Kad=0.0018, eokKcal_Mol=3.2, bVol=avoNum*pi*sigma^3/6=9.03, tK=583.15, eta=0.5
	! Yhb=14.826, ralph=2.3102, fAssoc=1.2142, XA=0.2628, aAssoc=-3.8711, dAlpha=3.6667, zAssoc=-5.406
    iErrCode=1 !XsolJre failed to converge
    iErrCode=0 !normal return is 0
!	do iType=1,nComps
!		if(nComps.eq.1.and.nDs(iType).ne.nAs(iType))then
!			iErrCode=1 !solving the general problem of asymmetric fractional bonding of multiple 
!			!types is equivalent to solving fractional bonding in mixtures.  
!			!todo: rewrite this code to use the general wertheim code for mixtures.  
!			!only then can we apply multiple bonding types to pure fluids.
!			call BeepMsg(errMsg(iErrCode))
!		endif
!	enddo
	!e.g. ethanol: rho=.68, tK=453 zAssoc=, aAssoc=
	if(LOUD.and.initCall)print*,'Wertheim: starting.'
	sumx=SUM( xFrac(1:nComps) )
	if( ABS(sumx-1) > zeroTol .and. sumx > zeroTol)xFrac(1:nComps)=xFrac(1:nComps)/sumx 
	bVolMix=SUM( xFrac(1:nComps)*bVolCc_mol(1:nComps) )
	!sumx=0
	nSitesTot=0
	nAccTot=0
	nDonTot=0
	do iComp=1,nComps
		!bVolCc_mol(iComp)=vMolecNm3(iComp)*avoNum
		!sumx=sumx+xFrac(iComp)
		if(bVolCc_mol(iComp) < 1)iErrCode=1
		!bVolMix=bVolMix+xFrac(iComp)*vMolecNm3(iComp)*avoNum
        nDonComp=0
        nAccComp=0
		do iType=1,nTypes(iComp)		
			nSitesTot=nSitesTot+nDegree(iComp,iType)*( nAcceptors(iComp,iType)+nDonors(iComp,iType) ) !ndHb(iType) should be zero if iType does not hBond.
			nAccComp=nAccComp+nDegree(iComp,iType)*nAcceptors(iComp,iType)
			nDonComp=nDonComp+nDegree(iComp,iType)*nDonors(iComp,iType)
			XA(iComp,iType)=1
			XD(iComp,iType)=1
			XC(iComp,iType)=1
			if(LOUD.and.initCall)print*,'Wertheim: iComp,iType,#Acc,#Don', iComp,iType,nAcceptors(iComp,iType),nDonors(iComp,iType)
		enddo
        bAssoc(iComp)=.false. ! then check if each component associates.
		if(	nAccComp*nDonComp > 0 )bAssoc(iComp)=.true.
		nAccTot=nAccTot+nAccComp
		nDonTot=nDonTot+nDonComp
	enddo
    zAssoc=0 !initialize to zero as base assumption.
    aAssoc=0
    uAssoc=0
	HyBond=.false. ! then check if any prospect for HyBonding, association or solvation 
	if(nDonTot*nAccTot > 0)HyBond=.true.
	if(.not.HyBond)then
		if(LouderWert.and.initCall)print*,'Wertheim: No HyBond. nDonTot,nAccTot=',nDonTot,nAccTot
		if(LouderWert.and.initCall)pause 'Wertheim: no HyBond. Dont come back!'
		initCall=0
		return !no need for wertheim in this case.
	endif
	if(LouderWert)print*,'nAcceptors(iComp1)=',( nAcceptors(1,iType),iType=1,nTypes(1) )
	if(LouderWert)print*,'nDonors(iComp1)=',( nDonors(1,iType),iType=1,nTypes(1) )
	if(initCall.and.LouderWert.and.HyBond)pause 'FuTptVtot: HyBonding identified in this fluid. Correct?'
	!pause 'Wertheim: .not.HyBond, but still in Wertheim???'
	if(tKelvin < 2)iErrCode=1 ! < 2 is too low even for He.
	if(eta > etaMax .or. eta < 0)iErrCode=1
	if(bVolMix.le.0)iErrCode=1
	if(sumx.le.0)iErrCode=1
	if(iErrCode==1)then
		if(LouderWert)write(*,'(a,f10.2)')' Wertheim: bVolMixCc_mol= ',bVolMix
		if(LouderWert)write(*,'(a,1PE11.4,a,0Pf7.2,a,f7.4)')' eta=',eta,' tKelvin=',tKelvin,' sumx=',sumx
		if(LouderWert)write(*,'(a,11F8.4)')' vMolecNm3=',(vMolecNm3(iComp),iComp=1,nComps)
		!call BeepMsg(errMsg(iErrCode))
		if(LouderWert)write(*,*)errMsg(iErrCode)
		if(LouderWert)pause
	endif
	if(iErrCode.ne.0)return !sorry for bad return
	
	rhoMol_cc=eta/bVolMix
	!pause 'Calling MEM1'
	if(LouderWert)call MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc,rLnPhiAssoc,iErrMEM1 )! XA,XD USE Assoc
	if(LouderWert)write(*,'(a,8f9.4)')' Wertheim-MEM1:eta,fAssoc,zAssoc,aAssoc,uAssoc',eta,fAssoc,zAssoc,aAssoc,uAssoc
	!pause 'Results from MEM1'
	!if(LouderWert.and.initCall)write(*,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
	if(LouderWert.and.initCall)pause 'Wertheim: MEM1 estimate. Calling MEM2'
	call MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr )!,rLnPhiAssoc,ier)
	if(LouderWert)write(*,'(a,8f9.4)')' Wertheim-MEM2:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
	!if(LouderWert)write(*,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
	!if(LouderWert)pause 'Wertheim: MEM2 estimate. Calling XsolJre'
    iDoSolvation=0
    if(isESD)iDoSolvation=0
    !if(isAcid > 0)iDoSolvation=1	 ! FYI: rmsErr on XsolJre2a is 1E-6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111111111
	if(iDoSolvation > 0)then  !checkup: 1E-111  ! Refine from ESD96 if necessary
		nSitesTot=0
		one=1
		half=one/2 !defined in GlobConst
		Call RdfCalc(rdfContact,dAlpha,eta)
		if(rdfContact < 1)then
			if(Loud)write(*,'(a,3f9.4)')' Wertheim: eta,rdfContact(<1)=',eta,rdfContact
			if(Loud)pause
			return
		endif

		isAcid=0
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)		
				if(idType(iComp,iType)==1603)isAcid=1
				nSitesTot=nSitesTot+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)+nDonors(iComp,iType) ) !ndHb(iType) should be zero if iType does not hBond.
				iComplex=nDonors(iComp,iType)*nAcceptors(iComp,iType)
				IF( isTpt ) THEN
					eHbKcal_mol(iComp,iType)=( eDonorKcal_mol(iComp,iType)+eAcceptorKcal_mol(iComp,iType) )/2.d0 
				ENDIF		
				if(iComplex.ne.0)iComplex=1  !in general, nDs or nAs could be 2, 3,... (e.g. water).  But complexation is either true or false.
				bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/1.987D-3)-1) !the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.

				SQARG=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
				IF(SQARG < 0)THEN
					if(LouderWert)write(*,*) 'neg alpha in wertheim. ETA=',ETA
					iErrCode=3
					if(LouderWert)pause
					RETURN
				ENDIF
				ralph(iComp,iType)=DSQRT( SQARG )
			enddo
		enddo
		!XA=1 !set to unity means no association	          
		!XD=1 !set to unity means no association	          
		!XC=1 !set to unity means no association	          
		if(LouderWert.and.initCall) write(*,'(a,5F10.6,3F10.1)')'Wertheim:solvation. X(1),X(2)=',xFrac(1),xFrac(2)
		!pause 'Results from MEM1'
		!print*,'Wert:uAssoc,CvAssoc',uAssoc,CvAssoc
		!print*,'uAssoc,b_xaDxa_db',uAssoc,b_xaDxa_db

		!Note: extension of uAssoc to asymmetric association below.

		!return	 !endpoint for geometric association rule, 
		!but I moved the if-then lower to catch cases when geometric rule is followed though not specified.(eg. methanol+nHexane).
		!Geometric initial guess is ready. Prepare to check if Old result is better.
	!
	!	!NOT Computing alphDA,alphAD because 4 subscripts are necessary for general form & tha's too many.
	!	DO I=1,nComps
	!		DO J=1,nComps
	!			dHadij=( eHbKcal_Mol(I)+eHbKcal_Mol(J) )/2.0*(1-aBip(I,J))
	!			kVi=bondVolNm3(i)*avoNum
	!			kVj=bondVolNm3(j)*avoNum
	!			yHBij=EXP(dHadij/tKelvin/1.987D-3)-1
	!			KVEij=SQRT(kVi*kVj)*yHBij
	!			iComplexAD=nAs(I)*nDs(J)
	!			if(iComplexAD.ne.0)iComplexAD=1
	!			iComplexDA=nDs(I)*nAs(J)
	!			if(iComplexDA.ne.0)iComplexDA=1
	!			alphAD(I,J)=iComplexAD*KVEij*ETA/bVolMix*rdfContact
	!			alphDA(I,J)=iComplexDA*KVEij*ETA/bVolMix*rdfContact
	!		enddo                                                  
	!	enddo

		!Check if Old guesses are better. During Z iteration, Old result may be close. 20201224: How can that be when we initialize to 1 in line 61?  Wiping out "old" stuff for Wertheim2d
		!Also, compute FAsolv, FDsolv, and rmsErrSolv, without SRCR.
		if(isAcid)call XsolJre2a(xFrac,rhoMol_cc,tKelvin,nComps,ierXsol)	 ! XA,XD,XC passed through module Assoc.
		!if(LouderWert)write(*,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
        if(ierXsol.ne.0)then
            iErrCode=1
			if(LouderWert)write(*,'(a,15i3)')' Wertheim: ier from Xsol=',ierXsol
            return !sorry.
        endif
		!TooDone: DIOXANE+CHLOROFORM SYSTEM WORKS LIKE EL2ed P799 but it overestimates the solvation.  Maybe because difunctionality of dioxane(?). 

		!C     SOLVING FOR AN INITIAL DXD GUESS USING ELLIOTT 96 EQUATION 19
		!C     NOTE THAT ALPHAij=SQRT(ALPHAi*ALPHAj) IS STILL APPLIED HERE.
		if(ierXsol)iErrCode=5 !Set error, but finish calculations anyway.

		h_nMichelsen=0
		hCC=0           !M&H Eq 11.
		aAssoc=0				 !M&H Eq 13.
		uAssoc=0
		do iComp=1,nComps		 
			do iType=1,nTypes(iComp)
				if(XA(iComp,iType) < zeroTol .and. LouderWert)pause 'Error in Wertheim: XA<1E-11'
				if(XD(iComp,iType) < zeroTol .and. LouderWert)pause 'Error in Wertheim: XD<1E-11'
				h_nMichelsen=h_nMichelsen+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1-XA(iComp,iType))+nDonors(iComp,iType)*(1-XD(iComp,iType)) )          !M&H Eq 11.
				hCC=hCC+xFrac(iComp)*nDegree(iComp,iType)*( (1-XC(iComp,iType)) )          !M&H Eq 11.
				sumDon=nDegree(iComp,iType)*nDonors(iComp,iType)*( DLOG(XD(iComp,iType))+(1-XD(iComp,iType))/2 )
				sumAcc=nDegree(iComp,iType)*nAcceptors(iComp,iType)*( DLOG(XA(iComp,iType))+(1-XA(iComp,iType))/2 )
				sumCrb=nDegree(iComp,iType)*( DLOG(XC(iComp,iType))+(1-XC(iComp,iType))/2 )
				aAssoc=aAssoc+xFrac(iComp)*(sumAcc+sumDon+sumCrb)
				xaOld(iComp,iType)=xa(iComp,iType)
				xdOld(iComp,iType)=xd(iComp,iType)
			enddo
		enddo
		zAcid = -half*(hCC)*dAlpha	!M&H Eq 10&13.
		zAssoc= -half*(h_nMichelsen)*dAlpha+zAcid	!M&H Eq 10&13.
		if(LOUDERwert)write(*,'(a,8E12.4)')' Wertheim:hAD+hDA,hCC,fAssoc,zAssoc,zAcid',h_nMichelsen,hCC,fAssoc,zAssoc,zAcid

		!if(LouderWert)write(*,'(a,8f9.4)')' From Xsol:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
		uAssocAD=0	 !cf ML Michelsen and EM Hendriks (M&H), Physical Properties from Association Models, FPE 180:165-174(2001)
		uAssocDA=0	 !use analogy to Passoc in vicinity of Eqs 6-11.
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)
				SUMDONi=nDegree(iComp,iType)*nDonors(iComp,iType)*XD(iComp,iType)
				sumAcci=nDegree(iComp,iType)*nAcceptors(iComp,iType)*XA(iComp,iType)
				DO jComp=1,nComps
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rhoMol_cc,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						SUMDONj=nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)
						sumAccj=nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)
						uAssocAD=uAssocAD+xFrac(iComp)*xFrac(jComp)*sumAcci*sumDonj*dAlphADij	 ! dAlpha = beps*d(alpha/dbeps), M&H(01)
						uAssocDA=uAssocDA+xFrac(iComp)*xFrac(jComp)*sumDoni*sumAccj*dAlphDAij
					enddo
				enddo
			enddo
		enddo
		uAssoc= -half*(uAssocAD+uAssocDA)	!this appears to be faulty. TODO: fix. JRE 20191008			
	endif
	if(LouderWert)write(*,'(a,8f9.4)')' Wertheim-Xsol:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
	if(LouderWert)pause 'Wertheim: After Calling XsolJre'
	initCall=0
	RETURN
	END	 ! subroutine Wertheim()
	      
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!	SUBROUTINE DESIGNED BY MARTY YUZWA 3/97                              
	!	Revised JRE 6/02 - use successive substitution/wegstein, switch to f90, 
	!		include CS or ESD rdf.
	!	Purpose: Solve for ln(PhiAssoc) of all components even when NDS.ne.NAS
	!   Input:
	!		XA,XD,X,NDS,NAS,ND,vMolecNm3,nComps,ETA,fAssoc from Wertheim call.
	!	Output:
	!		FUGASSOC	
	!																		 
	!	DXAN=DERIVATIVE OF XA WRT N(K) MULTIPLIED BY N  	                 
	!	iRdfOpt	- option for characterizing the radial distribution funciton at contact.	
	!			= 0, not specified => error
	!			= 1, ESD form
	!			= 2, Carnahan-Starling form
	!			= 3, activity coefficient model form (rdf=universal constant=g(eta=0.4)
	!			= 4, ESD form, geometric association rule.
	!                                                                          
	!     EQUATION USED FOR FUGASSOC IS FOUND IN YUZWA 97 EQUATION #12         
	!	References:
	!		Michelsen and Hendriks, Fluid Phase Equilibria 180 (2001) 165-174
	!
	!	ALL DERIVATIVES ARE ADDED BY AFG, Oct. 2009
	!	ALSO AFG COMBINED THE WERTHEIM ROUTINES OF SPEAD AND ESD INTO ONE SINGLE SUBROUTINE
	!
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
	SUBROUTINE WertheimFugc(xFrac,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	USE Assoc
	!USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	LOGICAL LOUDER
	DIMENSION xFrac(NMX),FUGASSOC(NMX) !,vMolecNm3(NMX),NAS(NMX),NDS(NMX),ND(NMX),ralph(NMX),XAmol(NMX),XDmol(NMX)
	DIMENSION dFUGASSOC_dRHO(NMX),dFUGASSOC_dT(NMX) !,FUGASSOC_num(NMX)
	DIMENSION dXAdRHO(NMX,maxTypes),dXDdRHO(NMX,maxTypes),dXCdRHO(NMX,maxTypes),dXAdT(NMX,maxTypes),dXDdT(NMX,maxTypes),dXCdT(NMX,maxTypes)
	COMMON/dalpha_dx/dalphad_dT,dalphda_dT,dalphcc_dT,dalphad_dRHO,dalphda_dRHO,dalphcc_dRHO
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	!COMMON/num/FUGASSOC_num

	!ETAp=ETA !this is necessary for old method

	!C     SOLVING FOR THE DERIVATIVE OF g (rdf) WRT ETA
	!C     THIS IS NOT EXACTLY EQUAL TO THE DERIVATIVE
	!C     IT IS EQUAL TO eta/g TIMES THE DERIVATIVE OF g WRT eta

	LOUDER=LouderWert
	!LOUDER=.TRUE.

	needDerivs=0
	if(needDerivs==0)return	! Let's forget about derivatives etc. for now and stick to VLE calculations.

	bVolMix=SUM(xFrac(1:nComps)*bVolCC_mol(1:nComps))
	Call RdfCalc(rdfContact,dAlpha,eta)
	rho=eta/bVolMix
	one=1.d0
	half=one/2.d0	 !defined in GlobConst 9/24/19

	call RdfCalc(rdfContact,dAlpha,eta)	!dAlpha=eta/alpha*(dAlpha/dEta)=1+eta/g*(dg/dEta)
	h_nMichelsen=0.d0           !M&H Eq 11.
	do iComp=1,nComps		 
		do iType=1,nTypes(iComp)
			h_nMichelsen=h_nMichelsen+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1.d0-XA(iComp,iType))+nDonors(iComp,iType)*(1.d0-XD(iComp,iType))+(1.d0-XC(iComp,iType)) )          !M&H Eq 11.
		enddo
	enddo
	if(LOUDER)write(*,'(a,6E12.4)')' WertheimFugc: fugAssocBefore=',(fugAssoc(i),i=1,nComps)
	aAssoc=0
	isAcid=0
	do iComp=1,nComps		 
		sumLnXi=0.d0
		do iType=1,nTypes(iComp)
			if(XA(iComp,iType) < zeroTol .and. LouderWert)pause 'Error in WertheimFugc: XA<zeroTol'
			if(XD(iComp,iType) < zeroTol .and. LouderWert)pause 'Error in WertheimFugc: XD<zeroTol'
			if(XC(iComp,iType) < zeroTol .and. LouderWert)pause 'Error in WertheimFugc: XC<zeroTol'
			if(XC(iComp,iType) < 1-zeroTol)then
				XCtemp=XC(iComp,iType)
				isAcid=1
			endif
			sumLnXi=sumLnXi+nDegree(iComp,iType)*( nAcceptors(iComp,iType)*LOG(XA(iComp,iType))+nDonors(iComp,iType)*LOG(XD(iComp,iType))+LOG(XC(iComp,iType)) )          !M&H Eq 13
		enddo
		aAssoc=aAssoc+xFrac(iComp)*sumLnXi
		!if(LOUDER)write(*,'(a,i3,4F10.7)')' Fugc: i,XA(i,j)=',iComp,(XA(iComp,j),j=1,nTypes(iComp))
		!if(LOUDER)write(*,'(a,i3,4F10.7)')' Fugc: i,XD(i,j)=',iComp,(XD(iComp,j),j=1,nTypes(iComp))
		!if(LOUDER)write(*,'(a,i3,4F10.7)')' Fugc: i,XC(i,j)=',iComp,(XC(iComp,j),j=1,nTypes(iComp))
		fugAssoc(iComp)=sumLnXi-half*h_nMichelsen*(dAlpha-1)*bVolCC_mol(iComp)/bVolMix  !Eq 13
	enddo
	if(LOUDER)then
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid==1)write(*,'(a,3E12.4)')' Fugc: XC= ',XCtemp
	endif
	if(LOUDER)write(*,'(a,6E12.4)')' WertheimFugc: hMich,aAssoc=',h_nMichelsen,aAssoc
	if(LOUDER)write(*,'(a,6E12.4)')' WertheimFugc: fugAssocAfter =',(fugAssoc(i),i=1,nComps)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Derivatives are added by AFG, Oct. 2009
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	isZiter=1 ! we don't need T and rho derivatives to compute SUM's
	SUMA=0.d0
	SUMD=0.d0
	SUMC=0.d0
	DO J=1,nComps
		DO jType=1,nTypes(J)
			call AlphaSp( isZiter,J,jType,J,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)			
			SUMA=SUMA+xFrac(J)*nAcceptors(J,jType)*nDegree(J,jType)*alphADij*(XA(J,jType))**2
			SUMD=SUMD+xFrac(J)*nDonors(J,jType)*nDegree(J,jType)*alphDAij*(XD(J,jType))**2
			SUMC=SUMC+xFrac(J)*nDegree(J,jType)*alphCCij*(XC(J,jType))**2
		enddo
	ENDDO
	isZiter=0 ! we need T derivatives for below
 	DO I=1,nComps						  !Elliott'96 Eq.37
		DO iType=1,nTypes(I)
			dXAdT(I,iType)=0.d0
			dXDdT(I,iType)=0.d0
			dXCdT(I,iType)=0.d0
 			DO K=1,nComps 
				DO kType=1,nTypes(K)
 					! First Derivative of X in respect to T while N & V are constant
					call AlphaSp(isZiter,I,iType,K,kType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
					dXAdT(I,iType)=dXAdT(I,iType)-xFrac(K)*nDonors(K,kType)*nDegree(K,kType)*XD(K,kType)*alphADij*dalphad_dT
					dXDdT(I,iType)=dXDdT(I,iType)-xFrac(K)*nAcceptors(K,kType)*nDegree(K,kType)*XA(K,kType)*alphDAij*dalphda_dT
					dXCdT(I,iType)=dXCdT(I,iType)-xFrac(K)*nDegree(K,kType)*XC(K,kType)*alphCCij*dalphcc_dT
				ENDDO
			ENDDO 
			! First Derivative of X in respect to RHO while N & T are constant
			dXAdRHO(I,iType)=-XA(I,iType)*( XD(I,iType)*(1.d0/XA(I,iType)-1.d0)*(1.d0+dLng) )/(1.d0+SUMD)/rho
			dXDdRHO(I,iType)=-XD(I,iType)*( XA(I,iType)*(1.d0/XD(I,iType)-1.d0)*(1.d0+dLng) )/(1.d0+SUMA)/rho
			dXCdRHO(I,iType)=-XC(I,iType)*( XC(I,iType)*(1.d0/XC(I,iType)-1.d0)*(1.d0+dLng) )/(1.d0+SUMC)/rho
		ENDDO
	ENDDO
	DO I=1,nComps
		DO iType=1,nTypes(I)						  !Elliott'96 Eq.37
			dXAdT(I,iType)=dXAdT(I,iType)*XA(I,iType)*XD(I,iType)/(1.d0+SUMD)
			dXDdT(I,iType)=dXDdT(I,iType)*XD(I,iType)*XA(I,iType)/(1.d0+SUMA)
			dXCdT(I,iType)=dXCdT(I,iType)*XC(I,iType)*XC(I,iType)/(1.d0+SUMC)
		ENDDO
	ENDDO
	dh_dT=0.d0
	dh_dRHO=0.d0
	h_nMichelsen=0.d0           !M&H Eq 11.
	do iComp=1,nComps		 
		do iType=1,nTypes(iComp)
			h_nMichelsen=h_nMichelsen+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1.d0-XA(iComp,iType))+nDonors(iComp,iType)*(1.d0-XD(iComp,iType))+(1.d0-XC(iComp,iType)) )          !M&H Eq 11.
			dh_dT=dh_dT-xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*dXAdT(iComp,iType)+nDonors(iComp,iType)*dXDdT(iComp,iType) ) 
 			!first derivative at constant T and N
			dh_dRHO=dh_dRHO-xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*dXAdRHO(iComp,iType)+nDonors(iComp,iType)*dXDdRHO(iComp,iType)+dXCdRHO(iComp,iType) )
		enddo
	enddo
	do iComp=1,nComps		 
		dsumLnXi_dT=0.d0
		dsumLnXi_dRHO=0.d0				 !M&H Eq 13.
		do iType=1,nTypes(iComp)
			! for M&H Eq.13, note that sum(Ai) is for sites on ith comp only.
			dsumLnXi_dT=dsumLnXi_dT+nDegree(iComp,iType)*( nAcceptors(iComp,iType)*dXAdT(iComp,iType)/XA(iComp,iType)+nDonors(iComp,iType)*dXDdT(iComp,iType)/XD(iComp,iType)+dXCdT(iComp,iType)/XC(iComp,iType))
			dsumLnXi_dRHO=dsumLnXi_dRHO+nDegree(iComp,iType)*(nAcceptors(iComp,iType)*dXAdRHO(iComp,iType)/XA(iComp,iType)+nDonors(iComp,iType)*dXDdRHO(iComp,iType)/XD(iComp,iType)+dXCdRHO(iComp,iType)/XC(iComp,iType))
		enddo
		do iType=1,nTypes(iComp)
			!first derivative in respect to T while N & V are constant
			dFUGASSOC_dT(iComp)=dsumLnXi_dT-half*(dalpha-1.d0)*bVolCC_mol(iComp)/bVolMix*dh_dT
			!first derivative in respect to RHO while N & T are constant
			dFUGASSOC_dRHO(iComp)=dsumLnXi_dRHO-half*((dalpha-1.d0)*bVolCC_mol(iComp)/bVolMix*dh_dRHO+h_nMichelsen*dAlph_deta*bVolCC_mol(iComp)/avoNum)
		enddo	
	enddo

	return 
	end           

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	SUBROUTINE DESIGNED BY JRE 9/05 using type-type basis instead of compd-compd basis                            
!C                                                                       
!C	THIS SUBROUTINE SOLVES FOR THE XD AND XA    
!C	IT USES Successive substitution for the off-diagonals and solves analytically for the XA(i) or XD(i)    
!C	1-XA(1,t1)=XA(1,t1)*X1D*fracAD+XA(1,t1)*sumDon => 1/X1A=(1+sumDon+X1D*fracAD) ; fracIJ, sumAcc, sumDon def below
!C	sumDon=sumDoni+sumDonj; 
!C	sumDoni=xFrac1*sum(tj.ne.ti)[XD(1,tj)*Nd(1,tj)+nDs(1,tj)*alphaAD(1,ti;i,tj)]
!C	sumDonj=Sum{xFracj*Sum(tj=1,nTypesj)[XD(j,tj)*Nd(j,tj)+nDs(j,tj)*alphaAD(i,ti;j,tj)]
!C	1-X1D=X1A*X1D*fracDA+X1D*sumAcc => 1/X1D=(1+sumAcc+X1A*fracDA)
!C  1/X1A=[1+sumDon+fracAD/(1+sumAcc+X1A*fracDA)] = F1A definition (X1A gives poss zero in denom)
!C  F1A=[1+sumDon+fracAD/(1+sumAcc+fracDA/F1A)]
!C  F1A*(1+sumAcc+fracDA/F1A)=[(1+sumDon)*(1+sumAcc+fracDA/F1A)+fracAD]
!C  F1A*(1+sumAcc)+fracDA=(1+sumDon)*(1+sumAcc+fracDA/F1A)+fracAD
!C  F1A^2*(1+sumAcc)+F1A*fracDA=(1+sumDon)*(1+sumAcc)*F1A+(1+sumDon)*fracDA+fracAD*F1A]
!C  F1A^2*(1+sumAcc)-F1A*[-fracDA+(1+sumDon)*(1+sumAcc)+fracAD]-(1+sumDon)*fracDA=0
!C  quadB= (1+sumDon)*(1+sumAcc)+fracAD-fracDA
!C  4ac  =-(1+sumDon)*(1+sumAcc)*fracDA*4
!C  X1A=1/F1A=2*(1+sumAcc)/[quadB+SQRT(quadB^2-4ac)]
!C                                                                          
!C	XA,XD are initial guesses on input and answers on output
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJre2a(xFrac,rho,tKelvin,nComps,ierCode)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	USE Assoc !XA,XD,XC,aBipAD,aBipDA
	IMPLICIT doublePrecision(A-H,K,O-Z)
	PARAMETER (iWegFreq=4)
	DIMENSION xFrac(NMX)
	DIMENSION fa(NMX,maxTypes),fd(NMX,maxTypes),faOld(NMX,maxTypes),fdOld(NMX,maxTypes)
	DIMENSION xaOld(NMX,maxTypes),xdOld(NMX,maxTypes),alphCCavg(maxTypes)
	DoublePrecision sumDoni(maxTypes),sumDonj(maxTypes),sumAcci(maxTypes),sumAccj(maxTypes)
	LOGICAL LOUDER
	data initCall/1/
	WegFreq(iter)=ABS(MOD(FLOAT(iter),FLOAT(iWegFreq)))	! WegFreq=0 if iter is multiple of iWegFreq
	LOUDER=LOUD
	LOUDER=LouderWert ! LouderWert is USEd from Assoc
	!LOUDER=.TRUE.
!	isEven(iter)= ( iter/2.D0 - iter/2 )*2.01
	ierCode=0
	wegParm=0
	rmsErr=111 !avoid termination of doWhile before doing anything.
	iter=0
	itMax=988
	isAcid=0
	do while(rmsErr > 1E-6.and.iter < itMax) ! Note: XA,XD,XC guesses are USEd from Assoc.
		iter=iter+1
		rmsErr=0
        rmsFerr=0
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)
				if(idType(iComp,iType)==1603)isAcid=1
				SUMDONi(iType)=0 !We must store these so we can update all types simultaneously on a consistent basis.
				SUMACCi(iType)=0
				do jType=1,nTypes(iComp)
					if(jType.ne.iType)then
						call AlphaSp(isZiter,iComp,iType,iComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						sumDoni(iType)=sumDoni(iType)+xFrac(iComp)*nDegree(iComp,jType)*nDonors(iComp,jType)*XD(iComp,jType)*ALPHADij
						sumAcci(iType)=sumAcci(iType)+xFrac(iComp)*nDegree(iComp,jType)*nAcceptors(iComp,jType)*XA(iComp,jType)*ALPHDAij
					endif
				enddo
				SUMDONj(iType)=0 ! SUMDONj =/= SUMDONi.
				SUMACCj(iType)=0
				alphCCavg(iType)=0
				!C	sumDon=sumDoni+sumDonj; 
				!C	sumDoni=xFrac1*sum(tj.ne.ti)[XD(1,tj)*Nd(1,tj)*nDs(1,tj)*alphaAD(1,ti;i,tj)]
				!C	sumDonj=Sum{xFracj*Sum(tj=1,nTypesj)[XD(j,tj)*Nd(j,tj)*nDs(j,tj)*alphaAD(i,ti;j,tj)]
				DO jComp=1,nComps
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						if(idType(iComp,iType)==1603 .and. idType(jComp,jType)==1603)then
							alphCCavg(iType)=alphCCavg(iType)+xFrac(jComp)*nDegree(jComp,jType)*alphCCij
						endif
					enddo
					if(jComp==2)ralphC=SQRT(  MaxVal( alphCCavg(1:nTypes(jComp)) )  )
					if(jComp.eq.iComp)cycle	!here we treat types on different comps as different (if no intra effect, they actually should have the same XA,XD)
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						SUMDONj(iType)=SUMDONj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*ALPHADij
						if(SUMDONj(iType)<0)then
							if(LouderWert)write(*,*)'iType,sumdonj(iType)',iType,sumdonj(iType)
							if(LouderWert)pause 'XsolJre:sumdon<0?'
						endif
						SUMACCj(iType)=SUMACCj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*ALPHDAij
						if(SUMACCj(iType)<0)then
							if(LouderWert)write(*,*)'iType,sumaccj(iType)',iType,sumaccj(iType)
							if(LouderWert)pause 'XsolJre:sumacc<0?'
						endif
					enddo      
				enddo
			enddo !terminate the loop over all types such that sumDon's and sumAcc's are all done at once before starting next stage
			do iType=1,nTypes(iComp)
				sumDon=sumDoni(iType)+sumDonj(iType)
				sumAcc=sumAcci(iType)+sumAccj(iType)
				call AlphaSp(isZiter,iComp,iType,iComp,iType,tKelvin,rho,bVolMix,alphADii,alphDAii,alphCCii,dAlphADii,dAlphDAii,dAlphCCii)
				fracAD=xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*alphADii
				fracDA=xFrac(iComp)*nDegree(iComp,iType)*nDonors(iComp,iType)*alphDAii
				fa(iComp,iType)=( XA(iComp,iType)-1/( 1+sumdon+fracAD*XD(iComp,iType) ) )/XA(iComp,iType) ! divide by XA so relative error is tracked
				fd(iComp,iType)=( XD(iComp,iType)-1/( 1+sumacc+fracDA*XA(iComp,iType) ) )/XD(iComp,iType)
                rmsFerr=rmsFerr+fa(iComp,iType)*fa(iComp,iType)+fd(iComp,iType)*fd(iComp,iType)
				wegParmAcc= 0.5D0
				wegParmDon= 0.5D0
				if(WegFreq(iter).lt.0.01)then !check for iter hit.
					dFa=fa(iComp,iType)-faOld(iComp,iType)  ! for many siteTypes, this will be zero (converged), so skip weg.
					dFd=fd(iComp,iType)-fdOld(iComp,iType)
					if(ABS(dFa).gt.1e-11)then ! .and. isEven(iter)==1)then
						changeAcc=fa(iComp,iType)/( dFa ) !must be inside if, or Old may not exist on first iteration.
						wegParmAcc= -changeAcc
                        if(wegParmAcc > 1) wegParmAcc=0.9D0
                        if(wegParmAcc < 0) wegParmAcc=0.1D0
					elseif(ABS(dFd).gt.1e-11)then  
						changeDon=fd(iComp,iType)/( dFd ) !must be inside if, or Old may not exist on first iteration.
						wegParmDon= -(changeDon) !/2
                        if(wegParmDon > 1) wegParmDon=0.9D0
                        if(wegParmDon < 0) wegParmDon=0.1D0
					endif
				endif	  
				xaOld(iComp,iType)=xa(iComp,iType)
				xdOld(iComp,iType)=xd(iComp,iType)
				faOld(iComp,iType)=fa(iComp,iType)
				fdOld(iComp,iType)=fd(iComp,iType)
				!solving the quadratic between XAi and XDi to get next guess
				quadB=(1+sumDon)*(1+sumAcc)+fracAD-fracDA
				sqArg=quadB*quadB+4*fracDA*(1+sumAcc)*(1+sumDon)
				if(sqArg<0)then
					if(LouderWert)write(*,*)'quadB,fracDA,sumAcc,sumDon',quadB,fracDA,sumAcc,sumDon
					if(LouderWert)pause 'XsolJre: sqArg < 0.'
				endif
				if( (1+sumAcc) < 0) sumAcc= -0.99 !pause 'XsolJre: 1+sumAcc < 0'
				if(ABS( quadB+SQRT(sqArg) ) < 1e-11)then
					xaNew=SQRT(1+sumAcc)/( (1+sumDon)*SQRT(1+sumAcc)+fracAD )
				else
					xaNew=2*(1+sumAcc)/( quadB+SQRT(sqArg) )
				endif
				if(xaNew<0)xaNew=1e-11
				XA(iComp,iType)=XA(iComp,iType)*wegParmAcc+(1-wegParmAcc)*xaNew
				if(XA(iComp,iType)<0)XA(iComp,iType)=1e-11	!note: WegParm>1 sometimes. CanNe cause overstep.
				if(XA(iComp,iType)>1)XA(iComp,iType)=1
				!XA(iComp)=XA(iComp)-(1-wegParm)*fa(iComp) !same as above
				!XD(iComp)=XD(iComp)-(1-wegParm)*fd(iComp)
				xdNew=1/( 1+sumAcc+fracDA*xaNew )  !using freshly computed XAi here makes this the analytical solution.
				if(xdNew<0)xdNew=1e-11
				XD(iComp,iType)=XD(iComp,iType)*wegParmDon+(1-wegParmDon)*xdNew
				if(XD(iComp,iType)<0)XD(iComp,iType)=1e-11	!note: WegParm>1 sometimes. Can cause overstep.
				if(XD(iComp,iType)>1)XD(iComp,iType)=1   	!note: WegParm<0 sometimes. Can cause overstep.
				!XD(iComp)=xdNew
				rmsErr=rmsErr+( xa(iComp,iType)-xaOld(iComp,iType) )**2
				rmsErr=rmsErr+( xd(iComp,iType)-xdOld(iComp,iType) )**2
				!Note: XCi=XCj b/c we assume alphCC is universal.  So we can factor out alphCCavg.
				!xCC=( -1+SQRT(1+4*alphCCavg(iType)) )/( 2*alphCCavg(iType) ) !divide by zero when type is not CC
				!xCC=( (1+4*alphCCavg(iType)-1) )/( 2*alphCCavg(iType)*(1+SQRT(1+4*alphCCavg(iType)) ) !divide by zero when type is not CC
				if( (1+4*alphCCavg(iType)) < 0) pause 'XsolJre: 1+4*alphCCavg(iType) < 0'
				xCC=2/( 1+SQRT(1+4*alphCCavg(iType)) )
				XC(iComp,iType)=xCC
			enddo  ! iType
		enddo  !  iComp
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsFerr=sqrt( rmsFerr/(2*nComps) )
		rmsOld=rmsErr
		if(LOUDER)write(*,'(a,i3,2E12.4)')' XsolJre2a: iter,rmsErr,rmsFerr=',iter,rmsErr,rmsFerr
		!if(LOUD.and.initCall)write(*,'(a,7F10.6,3f7.3)')' XAs&XDs,...',XA(1,2),XA(1,3),XA(2,2),XD(1,2),XD(2,2),rmsErr,wegParm,float(iter) !,deltaAlphaA(1,3) 
	enddo ! while(rmsErr > tol)
	if(LOUDER)then
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(*,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid==1)write(*,'(a,3E12.4)')' Xsol: alphCCavg,XC= ',ralphC*ralphC,xCC
		write(*,'(a,3E12.4)')' rmsErr,wegParm,iter',rmsErr,wegParm,float(iter) !,deltaAlphaA(1,3)
	endif 
	if(iter > itMax)then
		ierCode=2
		if(Louder)pause 'Error in XsolJre: iter>itMax'
	else
!		if(LOUD)write(*,*)'XsolJre converged. iter=',iter
!		if(LOUD)write(*,*)'xa1,xa2',XA(1),XA(2)
!		if(LOUD)write(*,*)'xd1,xd2',XD(1),XD(2)
	endif
	initCall=0
	RETURN
	END	 ! XsolJre2a

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJre1a(xFrac,rho,tKelvin,nComps,ierCode)	! Simple successive substitution
	USE Assoc !XA,XD,XC,aBipAD,aBipDA
	IMPLICIT doublePrecision(A-H,K,O-Z)
	PARAMETER (iWegFreq=99)
	DIMENSION xFrac(NMX)
	DIMENSION fa(NMX,maxTypes),fd(NMX,maxTypes),alphCCavg(maxTypes)
	DIMENSION faOld(NMX,maxTypes),fdOld(NMX,maxTypes) ,xaOld(NMX,maxTypes),xdOld(NMX,maxTypes)
	DoublePrecision sumDoni(maxTypes),sumDonj(maxTypes),sumAcci(maxTypes),sumAccj(maxTypes)
	WegFreq(iter)=ABS(MOD(FLOAT(iter),FLOAT(iWegFreq)))	! WegFreq=0 if iter is multiple of iWegFreq
	sumDoni = 0
	sumAcci = 0
	ierCode=0
	wegParm=0
	rmsErr=111 !avoid termination of doWhile before doing anything.
	iter=0
	itMax=99
	do while(rmsErr > 1E-6.and.iter < itMax)
		iter=iter+1
		rmsErr=0
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)
				SUMDONj(iType)=0 ! SUMDONj =/= SUMDONi.
				SUMACCj(iType)=0
				alphCCavg(iType)=0
				!C	sumDon=sumDoni+sumDonj; 
				!C	sumDoni=xFrac1*sum(tj.ne.ti)[XD(1,tj)*Nd(1,tj)*nDs(1,tj)*alphaAD(1,ti;i,tj)]
				!C	sumDonj=Sum{xFracj*Sum(tj=1,nTypesj)[XD(j,tj)*Nd(j,tj)*nDs(j,tj)*alphaAD(i,ti;j,tj)]
				DO jComp=1,nComps
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						if(idType(iComp,iType)==1603 .and. idType(jComp,jType)==1603)then
							alphCCavg(iType)=alphCCavg(iType)+xFrac(jComp)*nDegree(jComp,jType)*alphCCij
						endif
					enddo
					!if(jComp.eq.iComp)cycle	! sum i and j no different for simple iteration
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						SUMDONj(iType)=SUMDONj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*alphADij
						if(SUMDONj(iType)<0)then
							if(LouderWert)write(*,*)'iType,sumdonj(iType)',iType,sumdonj(iType)
							if(LouderWert)pause 'sumdon<0?'
						endif
						SUMACCj(iType)=SUMACCj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*alphDAij
						if(SUMACCj(iType)<0)then
							if(LouderWert)write(*,*)'iType,sumaccj(iType)',iType,sumaccj(iType)
							if(LouderWert)pause 'sumacc<0?'
						endif
					enddo ! jType
				enddo ! jComp
				!XD(iComp,iType)=1/( 1+sumAccj(iType) ) ! Acceptors are more likely to be in excess, so fix the donors first.     
			enddo !terminate the loop over all types such that sumDon's and sumAcc's are all done at once before starting next stage
			do iType=1,nTypes(iComp)
				sumDon=sumDoni(iType)+sumDonj(iType)
				sumAcc=sumAcci(iType)+sumAccj(iType)
				call AlphaSp(isZiter,iComp,iType,iComp,iType,tKelvin,rho,bVolMix,alphADii,alphDAii,alphCCii,dAlphADii,dAlphDAii,dAlphCCii)
				fa(iComp,iType)=( 1/XA(iComp,iType) - (1+sumdon) )
				fd(iComp,iType)=( 1/XD(iComp,iType) - (1+sumacc) )
				wegParm= 0.0
				if(WegFreq(iter).lt.0.01)then !check for iter hit.
					dFa=fa(iComp,iType)-faOld(iComp,iType)
					dFd=fd(iComp,iType)-fdOld(iComp,iType)
					if(ABS(dFa).gt.1e-11)then
						changeAcc=fa(iComp,iType)/( dFa ) !must be inside if, or Old may not exist on first iteration.
						wegParm= -changeAcc
					endif
					if(ABS(dFd).gt.1e-11)then
						changeDon=fd(iComp,iType)/( dFd ) !must be inside if, or Old may not exist on first iteration.
						wegParm= -changeDon
					endif
					wegParm= -0.1D0
				endif	  
				xaOld(iComp,iType)=xa(iComp,iType)
				xdOld(iComp,iType)=xd(iComp,iType)
				faOld(iComp,iType)=fa(iComp,iType)
				fdOld(iComp,iType)=fd(iComp,iType)
				xaNew=1/(1+sumDon)
				if(xaNew<0)xaNew=1e-11
				XA(iComp,iType)=XA(iComp,iType)*(wegParm)+(1-wegParm)*xaNew
				if(XA(iComp,iType)<0)XA(iComp,iType)=1e-11	!note: WegParm>1 sometimes. Can cause overstep.
				if(XA(iComp,iType)>1)XA(iComp,iType)=1
				!XA(iComp)=XA(iComp)-(1-wegParm)*fa(iComp) !same as above
				!XD(iComp)=XD(iComp)-(1-wegParm)*fd(iComp)
				xdNew=1/( 1+sumAcc )  !using freshly computed XAi here makes this the analytical solution.
				if(xdNew<0)xdNew=1e-11
				XD(iComp,iType)=XD(iComp,iType)*(wegParm)+(1-wegParm)*xdNew
				if(XD(iComp,iType)<0)XD(iComp,iType)=1e-11	!note: WegParm>1 sometimes. Can cause overstep.
				if(XD(iComp,iType)>1)XD(iComp,iType)=1   	!note: WegParm<0 sometimes. Can cause overstep.
				!XD(iComp)=xdNew
				!rmsErr=rmsErr+fa(iComp,iType)*fa(iComp,iType)
				!rmsErr=rmsErr+fd(iComp,iType)*fd(iComp,iType)
				rmsErr=rmsErr+(xa(iComp,iType)-xaOld(iComp,iType))**2
				rmsErr=rmsErr+(xd(iComp,iType)-xdOld(iComp,iType))**2
				!Note: XCi=XCj b/c we assume alphCC is universal.  So we can factor out alphCCavg.
				!xCC=( -1+SQRT(1+4*alphCCavg(iType)) )/( 2*alphCCavg(iType) ) !divide by zero when type is not CC
				!xCC=( (1+4*alphCCavg(iType)-1) )/( 2*alphCCavg(iType)*(1+SQRT(1+4*alphCCavg(iType)) ) !divide by zero when type is not CC
				if( (1+4*alphCCavg(iType)) < 0) pause 'Wertheim: 1+4*alphCCavg(iType)2nd < 0'
				xCC=2/( 1+SQRT(1+4*alphCCavg(iType)) )
				XC(iComp,iType)=xCC
			enddo
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsOld=rmsErr
		write(*,'(a,6F10.6,3F10.1)')' XAs&XDs:',XA(1,2),XA(1,3),XA(2,2),XD(1,2),XD(2,2),rmsErr,float(iter), wegParm !,deltaAlphaA(1,3) 
	enddo ! while(rmsErr > tol)
	if(iter > itMax)then
		ierCode=2
		if(LouderWert)pause 'Error in XsolJre: iter>itMax'
	else
!		if(LOUD)write(*,*)'XsolJre converged. iter=',iter
!		if(LOUD)write(*,*)'xa1,xa2',XA(1),XA(2)
!		if(LOUD)write(*,*)'xd1,xd2',XD(1),XD(2)
	endif
	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Subroutine AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	USE Assoc !aBipAD,aBipDA
	USE BIPs  !H((I,J) for ESD
	implicit doublePrecision(a-h,k,o-z)
	!COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	!COMMON/dalpha_dx/dalphad_dT,dalphda_dT,dalphcc_dT,dalphad_dRHO,dalphda_dRHO,dalphcc_dRHO
	!note: because aBip is on a type-type basis, we don't need to quadruply subscript it.
	indexi=localType( idType(iComp,iType) )	!reference to local type list so aBip matrix is small.
	indexj=localType( idType(jComp,jType) )	!reference to local type list so aBip matrix is small.
	YDi=EXP( eDonorKcal_Mol(iComp,iType)/(tKelvin*1.987D-3) ) - 1.d0
	YDj=EXP( eDonorKcal_Mol(jComp,jType)/(tKelvin*1.987D-3) ) - 1.d0
	YAi=EXP( eAcceptorKcal_Mol(iComp,iType)/(tKelvin*1.987D-3) ) - 1.d0
	YAj=EXP( eAcceptorKcal_Mol(jComp,jType)/(tKelvin*1.987D-3) ) - 1.d0
	dHdaij=( (eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(jComp,jType))*aBipDA(indexi,indexj) ) 
	dHadij=( (eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(jComp,jType))*aBipAD(indexi,indexj) ) 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	eta=rho*bVolMix
	call RdfCalc(eta,rdfContact,dAlpha)
	kVi=bondVolNm3(iComp,iType)*avoNum								   
	kVj=bondVolNm3(jComp,jType)*avoNum
	rootRhogKi=SQRT( rho*rdfContact*kVi )
	rootRhogKj=SQRT( rho*rdfContact*kVj )
	ralphAi=rootRhogKi*SQRT(  YAi  )
	ralphDi=rootRhogKi*SQRT(  YDi  )
	ralphAj=rootRhogKj*SQRT(  YAj  )
	ralphDj=rootRhogKj*SQRT(  YDj  )

	yHBadij=EXP(dHadij/tKelvin/1.987D-3)-1.d0
	yHBdaij=EXP(dHdaij/tKelvin/1.987D-3)-1.d0
	KVEadij=SQRT(kVi*kVj)*yHBadij
	KVEdaij=SQRT(kVi*kVj)*yHBdaij
	dKVEadij=SQRT(kVi*kVj)*(dHadij/tKelvin/1.987D-3)*(yHBadij+1)   ! beta*dKVE/dBeta
	dKVEdaij=SQRT(kVi*kVj)*(dHdaij/tKelvin/1.987D-3)*(yHBdaij+1)
	iComplexAD=nAcceptors(iComp,iType)*nDonors(jComp,jType)
	if(iComplexAD > 0)iComplexAD=1
	iComplexDA=nDonors(iComp,iType)*nAcceptors(jComp,jType)
	if(iComplexDA > 0)iComplexDA=1
    !if(iComp==jComp .and. iType/=jType)iComplexAD=0  !ignore intramolecular association for now. !! This would exclude e.g. ethoxyGlycol ether from solvating with the hydroxy in pure fluid. 
    !if(iComp==jComp .and. iType/=jType)iComplexDA=0  !ignore intramolecular association for now.
	!alphADij=iComplexAD*KVEadij*rho*rdfContact
	!alphDAij=iComplexDA*KVEdaij*rho*rdfContact
	alphADij=iComplexAD*( ralphAi*ralphDj+aBipAD(indexi,indexj)*KVEadij*rho*rdfContact )
	alphDAij=iComplexAD*( ralphDi*ralphAj+aBipDA(indexi,indexj)*KVEdaij*rho*rdfContact )
	if(alphADij<0 .or. alphDAij<0)then
		if(LouderWert)write(*,*)'alphADij,alphDAij',alphADij,alphDAij
		if(LouderWert)pause 'AlphaSp: -ve alpha? Thats weird.'
	endif
	if(isZiter)return
	!ralph=(rootRhogK*Y)^0.5
	!dLnRalph=(beta/ralph)*0.5/ralph*dAlpha/dBeta=	0.5*dLnAlpha/dLnBeta 
	dLnAlphAi=eAcceptorKcal_mol(iComp,iType)/(1.987D-3*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphAi > 1.D-4)dLnAlphAi=dLnAlphAi*(YAi+1)/YAi
	dLnAlphDi=eDonorKcal_mol(iComp,iType)/(1.987D-3*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphDi > 1.D-4)dLnAlphDi=dLnAlphDi*(YDi+1)/YDi
	dLnAlphAj=eAcceptorKcal_mol(jComp,jType)/(1.987D-3*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphAj > 1.D-4)dLnAlphAj=dLnAlphAj*(YAj+1)/YAj
	dLnAlphDj=eDonorKcal_mol(jComp,jType)/(1.987D-3*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphDj > 1.D-4)dLnAlphDj=dLnAlphDj*(YDj+1)/YDj

	dAlphADij=iComplexAD*0.5*(ralphDj*ralphAi*dLnAlphAi+ralphAi*ralphDj*dLnAlphDj)+dKVEadij*rho*rdfContact
	dAlphDAij=iComplexDA*0.5*(ralphAj*ralphDi*dLnAlphDi+ralphDi*ralphAj*dLnAlphAj)+dKVEdaij*rho*rdfContact
	alphCCij=0
	if(idType(iComp,iType).eq.1603 .and. idType(jComp,jType).eq.1603)then 
		dHCCij=3*( eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(iComp,iType) )/2 
		!dHccij=3.00*eHbKcal_Mol(iComp,iType) !we can still vary AD and CC association ~independently by changing bondVolNm3 in ParmsHbVv.txt b/c bondVolCC=0.01=universal constant
		yHBccij=EXP(dHccij/tKelvin/1.987D-3)-1
		KVEccij=bondVolNm3(iComp,iType)*avoNum*yHBccij ! value for CC bonding volume taken from AD bonding volume.
		dKVEccij=bondVolNm3(iComp,iType)*avoNum*(dHccij/tKelvin/1.987D-3)*(yHBccij+1) !universal value for CC bonding volume
		alphCCij=KVEccij*rho*rdfContact
		dAlphCCij=dKVEccij*rho*rdfContact
	endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2009
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! First derivatives in respect to T while N & V are constant
	dyHBadij_dT= -(dHadij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHadij/tKelvin/1.987D-3)
	dyHBdaij_dT= -(dHdaij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHdaij/tKelvin/1.987D-3)
	dyHBccij_dT= -(dHccij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHccij/tKelvin/1.987D-3)

	dKVEadij_dT=SQRT(kVi*kVj)*dyHBadij_dT
	dKVEdaij_dT=SQRT(kVi*kVj)*dyHBdaij_dT
	dKVEccij_dT=SQRT(kVi*kVj)*dyHBccij_dT

	if (KVEadij > zeroTol)	then
		dalphad_dT=dKVEadij_dT/KVEadij   !Actually this is equal to T/alpha*dalpha_dT
	else 
		dalphad_dT=0
	endif
	if (KVEdaij > zeroTol)	then
		dalphda_dT=dKVEdaij_dT/KVEdaij   !Actually this is equal to T/alpha*dalpha_dT
	else 
		dalphda_dT=0
	endif
	if (KVEccij > zeroTol)	then
		dalphcc_dT=dKVEccij_dT/KVEccij   !Actually this is equal to T/alpha*dalpha_dT
	else 
		dalphcc_dT=0
	endif
										 !g=(1-1.9eta)^(-1) => dg/dEta= 1.9/(1-1.9eta)^2
	! First derivatives in respect to RHO while N & T are constant
	if (kVj > zeroTol)	then
		dg_dEta=(dAlpha-1)*rdfContact/eta ! = (dLng/dLnEta)*g/eta
		if(rho < zeroTol .or. rdfContact < zeroTol .and. LouderWert)write(*,'(a,2E12.4)')' AlphaSp: rho,rdfContact=',rho,rdfContact
		dalphad_dRHO=(rdfContact+rho*dg_deta*bVolMix)/(rho*rdfContact)   !Actually this is equal to RHO/alpha*dalpha_dRHO
		dalphda_dRHO=dalphad_dRHO
		dalphcc_dRHO=dalphad_dRHO
	else 
		dalphad_dRHO=0.d0
		dalphda_dRHO=0.d0
		dalphcc_dRHO=0.d0
	endif

	!if(LouderWert)write(*,'(a,4i3,2E12.4)')' AlphaSp: i,k,j,l,alphAD,alphDA=',iComp,iType,jComp,jType,alphADij,alphDAij

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Subroutine GetAssocBips(bipHbFile,aBip,ier)  !idLocalType,nTypesTot are in common/assoc
	USE GlobConst, only: avoNum,zeroTol,bVolCC_mol,ID
	USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000) ! listPool could be 10000x10000, but very few BIP corrections in reality. 
	character bipHbFile*88,dumString*123
	logical switched,LOUDER
	!character*123 ErrMsg(11)
	dimension idBinarY(listPool)
	doublePrecision KIJDB(listPool),KTIJDB(listPool) !,KIJ(NMX,NMX),KTIJ(NMX,NMX)
	dimension aBip(MaxTypes,MaxTypes)!,aBipDA(NMX,NMX) !,xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)
	!common/HbParms/dHkcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	LOUDER=LOUD
	!LOUDER=.TRUE.
	ier=0
	OPEN(40,FILE=BipHbFile)
	read(40,'(a)',ERR=861)dumString !1st line is a dummy header
	item=1
	do while(item > 0)
		read(40,*,ERR=861,end=200) idDon,idAcc,KIJDB(item),KTIJDB(item)
		idBinarY(item)=10000*idDon+idAcc
		item=item+1	!increment for subsequent loop if any.
		cycle
200		exit
	enddo
	nHbBips=item-1
	close(40)
	
	do iType=1,nTypesTot
		do jType=1,nTypesTot
			aBip(iType,jType)=0 
			switched=.FALSE. !disable switching b/c ADij is asymmetric. ie. the acceptor on i(eg.O=< ) and donor on j(eg.OH) may be different from the acceptor on j(OH) and donor on i(O=,has no donor).
			!if(idLocalType(iType).lt.idLocalType(jType))switched=.TRUE.
			idBin=10000*idLocalType(iType)+idLocalType(jType)
			if(switched)idBin=10000*idLocalType(jType)+idLocalType(iType)
			do item=1,nHbBips
				if(idBinarY(item)==idBin)then
					aBip(iType,jType)=KIJDB(item) 
					IF(LOUDER)WRITE(*,'(a,2i5,f8.3)')' GetAssocBIPs: FOUND - idi,idj,BipIJ ',idLocalType(iType),idLocalType(jType),KijDB(item)
					exit !found it so terminate the item loop
				endif
			enddo
			!Omit below b/c aBipAD is symmetric (???)
			!if(switched)then
			!	HTIJ(i,j)=xsTau(i,j)
			!endif
		enddo
	enddo
	if(LOUDER)THEN
		write(*,'(7x,11i7)')(idLocalType(iType),iType=1,nTypesTot)
		DO jType=1,nTypesTot
			write(*,'(i7,11f7.3)')idLocalType(jType),(aBip(jType,iType),iType=1,nTypesTot)
		enddo
	endif
	return
861	continue
	if(LOUD)write(*,*)'GetAssocBips: error opening file=',bipHbFile
	ier=1
	return 
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C

