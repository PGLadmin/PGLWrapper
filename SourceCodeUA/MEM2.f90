
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!  ELLIOTT'S SUBROUTINE 
	!  20221202 Initial adaptation of MEM2 for speadmd. 
	!  20230114 Unified basis for ESD96 and ESD-MEM2.
	SUBROUTINE MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr )
	USE GlobConst
	USE Assoc ! XA,XD,XC,NAcceptors(nmx),NDonors(nmx),NDegree(nmx),...
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (FA,FD) AND properties (zAssoc, lnPhiAssoc,...) given T,rho,x
	!  Reference: Elliott, IECR, 61(42):15724 (2022). doi: 10.1021/acs.iecr.2c02723
	!  INPUT:
	!  isZiter = 1 if only zAssoc,aAssoc are required (for zFactor iterations), 0 for lnPhiAssoc,uAssoc, -1 if Cv etc are required.
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
	DoublePrecision xFrac(nmx),ralphA(nmx,maxTypes),ralphD(nmx,maxTypes),rLnPhiAssoc(nmx) !,KVE(nmx)
	DoublePrecision tKelvin,rhoMol_cc,zAssoc,aAssoc,uAssoc,Fwertheim,xDiff !,aAssocPas,uAssocPas,zAssocPas
	DoublePrecision eta,picard,Ftol,bVolMix,rdfContact,dAlpha,ralphAmax,ralphDmax,avgNAS,avgNDS,FA0,FD0,FC0_ralphC,sqArg
	DoublePrecision FA,betadFD_dBeta,FD,betadFA_dBeta,FC,betadFC_dBeta,epsCC,ralphC,FAold,avgFo,delFo,delFA,delFD,error,errOld
	DoublePrecision sqArgA,sqArgD,sumA,sumD,errA,errD,FDold,FA1,XCtemp,hAD,hDA,hCC,bepsA,bepsD,bepsC,sumLnXi,dLnAlpha_dLnBeta
	DoublePrecision, SAVE:: xOld(nmx),ralphAmean,ralphDmean,etaOld,rdfOld ! "SAVE" is equivalent to the static attribute
	Integer, SAVE:: IDold(nmx)  								            ! static is equivalent to the "SAVE" attribute
	Integer isZiter,nComps,iErr, iComp,isDonor,isAcceptor,iErrMEM1,iType,i,j,itMax,nIter
	LOGICAL LOUDER,isAcid,moreDonors
	Character*123 errMsg(22)
	data etaOld,rdfOld,ralphAmean,ralphDmean/0.3583,3.0,1.0,1.0/ !This makes these values static for Intel Compiler without STATIC.
	!common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	iErr=0
	errMsg( 1)=' MEM2: ERROR - NO CNVRG. ' ! Warning because future calls (e.g., iterating Z) may converge and remove the issue.
	errMsg( 2)=' MEM2: Failed bond site balance.' ! Warning because some approximations might not require balance. 
	errMsg(11)=' MEM2: rdfContact<1'
	errMsg(12)=' MEM2: sqArg for ralphMean update.'
	errMsg(13)=' MEM2: sqArg(A or D)'
	errMsg(14)=' MEM2: XA,XD, or XC < 0' 
	errMsg(15)=' MEM2: etaOld < 0?' 
	errMsg(16)=' MEM2: Sorry. Derivatives beyond zAssoc and uAssoc have not been implemented yet.'
	if(isZiter < 0)then
		iErr=16
		if(LOUDER)write(dumpUnit,601)TRIM(errMsg(iErr))
		return
	endif 
601	format(1x,a,8E12.4)
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LOUDERWert
	!LOUDER= .TRUE. ! for local debugging.
	FC=0 !initialize
	FA=0
	FD=0
    zAssoc=0
    aAssoc=0
    uAssoc=0
	rLnPhiAssoc(1:nComps)=0
    iErr=0  
	picard=0.83D0
	Ftol=1.D-7
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	Call RdfCalc(rdfContact,dAlpha,eta)	!dAlpha = dLnAlpha/dLnRho = 1+dLn(g)_dLn(eta)
	if(rdfContact < 1)then
		iErr=11
		if(LOUDER)write(dumpUnit,'(a,3f9.4)')' MEM2: eta,rdfContact(<1)=',eta,rdfContact
		if(LOUDER)write(dumpUnit,601)TRIM(errMsg(iErr))
		return
	endif
610	format(1x,2i5,i12,i8,i11,i8,2E14.4)
	if(LOUDER)write(dumpUnit,601)' MEM2: T,rho,eta,g^c',tKelvin,rhoMol_cc,eta,rdfContact
	!XA=1   				   ! initializing the entire array is slow.
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
				epsCC=3*( eAcceptorKcal_mol(iComp,iType)+eDonorKcal_mol(iComp,iType) )/2
				!ralphC is not subscripted because 1603 is the only acid type and always the same.
				ralphC=SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*( EXP(epsCC/tKelvin/RgasCal*1000)-1.d0 )  ) 
				! Note the "3*" in the exponent, which makes the C bond "special."
				if(LOUDER)write(dumpUnit,'(a,F10.6,3E12.4)')' MEM2: epsCC,bondVol,ralphC',epsCC,bondVolNm3(iComp,iType),ralphC
			endif
			isAcceptor=0
			if(nAcceptors(iComp,iType) > 0)isAcceptor=1
			Fwertheim=EXP(eAcceptorKcal_mol(iComp,iType)/tKelvin/RgasCal*1000)-1.d0
			!sqArg=rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*Fwertheim
			ralphA(iComp,iType)=isAcceptor*SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*Fwertheim  )
			if(ralphA(iComp,iType) > ralphAmax)ralphAmax=ralphA(iComp,iType)
			isDonor=0
			if(nDonors(iComp,iType) > 0)isDonor=1
			Fwertheim=EXP(eDonorKcal_mol(iComp,iType)   /tKelvin/RgasCal*1000)-1.d0
			!sqArg=rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*Fwertheim
			!if(sqArg < 0)write(dumpUnit,form612)' MEM2: iComp,iType,T,eDonor',iComp,iType,tKelvin,eDonorKcal_Mol(iComp,iType)
			ralphD(iComp,iType)=isDonor*   SQRT(  rhoMol_cc*rdfContact*bondVolNm3(iComp,iType)*avoNum*Fwertheim  )
			if(ralphD(iComp,iType) > ralphDmax)ralphDmax=ralphD(iComp,iType)
			if(LOUDER)write(dumpUnit,form612)' MEM2:iComp,iType,ralphA,ralphD',iComp,iType,ralphA(iComp,iType),ralphD(iComp,iType)
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
			if(idType(i,j)==1603)FC0_ralphC=FC0_ralphC+xFrac(i)*nDegree(i,j) !one C allowed on type 1603 => nCarboxy(i,j)=1 always.
			avgNDS=avgNDS+xFrac(i)*nDegree(i,j)*nDonors(i,j)	
			avgNAS=avgNAS+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)
		enddo	
	enddo
	moreDonors=.FALSE.
	if(avgNAS < avgNDS)moreDonors=.TRUE.
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
	if(LOUDER)write(dumpUnit,'(a,2E12.4,1x,8i6,1x)')' MEM2: x1,x1old,ID,IDold ',xFrac(1),xOld(1),(ID(i),i=1,nComps),IDold(1:nComps)
	if(LOUDER)write(dumpUnit,601)' MEM2: etaOld,rdfOld ',etaOld,rdfOld
	xDiff=SUM( (xFrac(1:nComps)-xOld(1:nComps))**2 )
	if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. xDiff > Ftol/10.or.etaOld.le.0  )then	
	!Only use default estimate if compounds or composition have changed.
		if(LOUDER)write(dumpUnit,form611)' MEM2:New IDs or X ',SUM( ID(1:nComps)-IDold(1:nComps) ),xDiff 
		ralphAmean=ralphAmax
		ralphDmean=ralphDmax
		if(moreDonors)then
			ralphDmean=ralphDmean*avgNAS/avgNDS
		else
			ralphAmean=ralphAmean*avgNDS/avgNAS
		endif
	!if(  SUM( ID(1:nComps)-IDold(1:nComps) )/=0 .or. xDiff ) > Ftol/10  )then	!Use default estimate if compds or x's have changed
		if(LOUDER)write(dumpUnit,601)' MEM2: fresh ralphAmean,ralphDmean= ',ralphAmean,ralphDmean
	else !if( SUM(xFrac(1:nComps)) < 0)then
		if(etaOld > 0)then ! adapt old values.to accelerate Z iterations
			sqArg=eta/etaOld*rdfContact/rdfOld
			if(sqArg < zeroTol)then
				iErr=12
				if(LOUDER)write(dumpUnit,form610)' MEM2: sqArg err ralphMean. eta,etaOld,rdf,rdfOld',eta,etaOld,rdfContact,rdfOld
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
	nIter=0
	avgFo=(ralphDmean*FD0+ralphAmean*FA0)/2
	if(LOUDER)write(dumpUnit,'(a,f8.2,11f8.5)')' MEM2:T,eta',tKelvin,eta !,(etaPure(i),i=1,nComps)
	error=1234
	itMax=88
	if(LOUDER)write(dumpUnit,form611)' iter,ralphAmean,ralphDmean,FA,FD,FA0,FD0  ',nIter,ralphAmean,ralphDmean,FA,FD,FA0,FD0
	do while(ABS(error)>Ftol.and.nIter<itMax)
		nIter=nIter+1
		!if(nIter > 33)picard=0.5D0
		delFo=(ralphDmean*FD0-ralphAmean*FA0)
		delFA=1+delFo
		delFD=1-delFo
		sqArgA=delFA*delFA+4*ralphAmean*FA0
		sqArgD=delFD*delFD+4*ralphDmean*FD0
		if(sqArgA < zeroTol .or. sqArgD < zeroTol)then
			if(LOUDER)write(dumpUnit,form611)' MEM2: sqArg(A or D). iter,ralphA,FA0,ralphD,FD0',nIter,ralphAmean,FA0,ralphDmean,FD0
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
		if(LOUDER)write(dumpUnit,form611)' iter,ralph,ralphD,FA,FD,errA,errD',nIter,ralphAmean,ralphDmean,FA,FD,errA,errD
		ralphDmean= picard*(-1+FA0/sumD)/(FD+1D-9)+(1-picard)*ralphDmean !Using new sumD to compute new ralphMeans.
		if(ralphDmean < 0)ralphDmean=(-1+FA0/sumD)/(FD+1D-9)
		if(ralphDmean < 0)ralphDmean=zeroTol
		ralphAmean= picard*(-1+FD0/sumA)/(FA+1D-9)+(1-picard)*ralphAmean !add 1D-9 to avoid possible zero divide if FD->0
		if(ralphAmean < 0)ralphAmean=(-1+FD0/sumA)/(FA+1D-9)
		if(ralphAmean < 0)ralphAmean=zeroTol
		if(nIter==1)FA1=FA
	enddo !while abs(error) > 1.D-9. 
	if(nIter>itMax-1)then
		iErr=1	!warning level
		if(LOUDER)write(dumpUnit,*)'ERROR - NO CNVRG. nIter,FA1,FA=',nIter,FA1,FA
		!GOTO 86
	ENDIF
	if(nIter > itMax-1 .and. LOUDER)write(dumpUnit,form613)' iter,id1,id2,eta,x1,tKelvin,ralphA,ralphD,FA,FD,errA,errD',&
	                                                       nIter,ID(1),ID(2),xFrac(1),tKelvin,ralphAmean,ralphDmean,FA,FD,errA,errD
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
			betadFD_dBeta=betadFD_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralphA(i,j)*0.5d0*dLnAlpha_dLnBeta
			bepsD=(eDonorKcal_mol(i,j)/RgasCal*1000/tKelvin)
			dLnAlpha_dLnBeta=DEXP(bepsD)
			if(bepsD > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsD)/(dLnAlpha_dLnBeta-1)
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j) * nDonors(i,j) *XD(i,j)*ralphD(i,j)*0.5d0*dLnAlpha_dLnBeta  
			! the 0.5 comes from dLnRalph = 0.5*dLnAlpha
			if(idType(i,j)==1603)then
				bepsC=3*(eDonorKcal_mol(i,j)+eDonorKcal_mol(i,j))/(2*RgasCal/1000*tKelvin)
				dLnAlpha_dLnBeta=DEXP(bepsC)
				if(bepsC > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsC)/(dLnAlpha_dLnBeta-1)
				betadFC_dBeta=betadFC_dBeta+xFrac(i) *nDegree(i,j) *XC(i,j)*ralphC*dLnAlpha_dLnBeta
			endif
		enddo
		if(LOUDER)write(dumpUnit,611)'i,sumLnXi=',i,sumLnXi
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA+hCC)/2*(dAlpha-1)*bVolCc_mol(i)/bVolMix ! Eq. 42
	enddo
611 format(1x,a,i11,12E12.4)
	aAssoc=  aAssoc+(hAD+hDA+hCC)/2  ! Eqs. 1 & 40. Validation that dSUM(xi*SumLnXi)/dBeta=0 is given at end of paper.
	!aAssocPas=aAssoc
	if(isZiter==1)return 
	!(hAD+hDA)/2=hAD = FA*FD by Eq 40.
	uAssoc= -( FA*betadFD_dBeta+FD*betadFA_dBeta+FC*betadFC_dBeta/2 ) ! Eqs. 43,44
	!uAssocPas=uAssoc
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' MEM2: aAssoc,fugAssoc= ',aAssoc,(rLnPhiAssoc(i),i=1,nComps)
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' MEM2: BdFA/dB,BdFD/dB=',betadFA_dBeta,betadFD_dBeta
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
	SUBROUTINE MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc,rLnPhiAssoc,iErr )
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID,dumpUnit,RgasCal,form612
	USE Assoc
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (fAssoc) AND zAssoc given rho,VM
	Implicit NONE !DoublePrecision(A-H,K,O-Z)
	!PARAMETER(nmx=55)
	DoublePrecision xFrac(nmx),ralph(nmx,maxTypes),rLnPhiAssoc(nmx) !,KVE(nmx)
	DoublePrecision tKelvin,rhoMol_cc,zAssoc,aAssoc,uAssoc  ,FA,FD,betadFD_dBeta,prefix  !,aAssocPas,uAssocPas,zAssocPas 
	DoublePrecision eta,bVolMix,rdfContact,dAlpha,FA0 !,picard,ralphAmax !,Ftol,ralphDmax,avgNAS,avgNDS !,ralphDmean,FD0,FC0_ralphC
	DoublePrecision betadFA_dBeta,sqArg,errOld
	DoublePrecision sumA,hAD,hDA,bepsA,sumLnXi,dLnAlpha_dLnBeta,fOld,fAssoc1,fAssoc,ERR,change,bigYhb
	DoublePrecision, STATIC:: xOld(nmx),ralphMean,etaOld,rdfOld ! static is equivalent to the "SAVE" attribute
	Integer, STATIC:: IDold(nmx)  								            ! static is equivalent to the "SAVE" attribute
	Integer isZiter,nComps,iErr, iComp,iType,i,j,itMax,nIter,nSitesTot,iComplex,iErrCode !,iErrMEM,moreDonors,isDonor,isAcceptor
	LOGICAL LOUDER,isAcid
	!common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	!  NDi     = THE DEGREE OF POLYMERIZATION OF COMPO i
	!  fAssoc  = THE CHARACTERISTIC ASSOCIATION = 1/RALPHi(1/XAi-1)
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei/(1-1.9eta)
	!  ier     = 4  eta > 0.53
	LOUDER=LOUD	! from GlobConst
	!LOUDER=LOUDERWert
	!LOUDER= .TRUE. ! for local debugging.
	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
	eta=rhoMol_cc*bVolMix
	Call RdfCalc(rdfContact,dAlpha,eta)
	if(rdfContact < 1)then
		if(LOUDER)write(dumpUnit,610)' MEM1: Returning zero. eta,rdfContact(<1)=',eta,rdfContact
		if(LOUDER)write(dumpUnit,*)
		return
	endif
	zAssoc=0
	aAssoc=0
	uAssoc=0
	rLnPhiAssoc(1:nComps)=0
    iErr=0  
	!Note: 1st approximation assumes equal donors and acceptors on each site.
	nSitesTot=0
	isAcid=.FALSE.
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)		
			nSitesTot=nSitesTot+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)+nDonors(iComp,iType) ) 
			!ndHb(iType) should be zero if iType does not hBond.
			if(idType(iComp,iType)==1603)isAcid=.TRUE.
			iComplex=nDonors(iComp,iType)*nAcceptors(iComp,iType)
			eHbKcal_mol(iComp,iType)=( eDonorKcal_mol(iComp,iType)+eAcceptorKcal_mol(iComp,iType) )/2.d0 
			if(iComplex.ne.0)iComplex=1  !nDs or nAs could be 2, 3,... (e.g. water).  But complexation is either true or false.
			bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/RgasCal*1000)-1) 
			! iComplex permits eHb to be non-zero, but still have zero complexation. 
			!eg. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0. So complexation can occur between nAs1*nDs2, but not nDs1*nAs2
			sqArg=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
			IF(sqArg < 0)THEN
				if(LOUDER)write(dumpUnit,*) 'neg alpha in wertheim. eta=',eta
				iErrCode=3
				if(LOUDER)write(dumpUnit,*)
				RETURN
			ENDIF
			ralph(iComp,iType)=DSQRT( sqArg )
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
				FA0=FA0+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType) !
				if(LOUDER)write(dumpUnit,612)' MEM1: iComp,iType,ralph',iComp,iType,ralph(iComp,iType)
				if(ralph(iComp,iType) > ralphMean)ralphMean=ralph(iComp,iType) ! Infinity norm for initial estimate.
			endif
		enddo
	enddo
	errOld=fAssoc-FA0
	fOld=fAssoc
	SUMA=FA0
	IF(ABS(errOld) < 1.D-9.and.LOUDER)write(dumpUnit,610)' MEM1: No assoc? FA0=',FA0	!don't iterate if errOld < 1D-9
	IF(ABS(errOld) > 1.D-9)then	!don't iterate if errOld < 1D-9
		fAssoc1=2*FA0/( 1+DSQRT(1+4*ralphMean*FA0) )  ! Eq.23 where FD0=FA0. This is exact for binary like benzene+methanol.
		fAssoc=fAssoc1 
		nIter=0
		ERR=1111
		itMax=123
		do while(ABS(ERR) > 1.D-9.AND.nIter < itMax) ! Don't use ERR/fAssoc because:(1)fAssoc->0 (2) max(fAssoc) ~ 1.
			nIter=nIter+1                  
			sumA=0 ! = SUM( xi*NDi*NAi*ralphi/(1+ralphi*FA) ) = FA = FA0 if no association
			DO iComp=1,nComps
				do iType=1,nTypes(iComp)      
					iComplex=nAcceptors(iComp,iType)*nDonors(iComp,iType)
					if(iComplex > 0)then !only apply symmetric rule when acceptors AND donors on a site.
						prefix=xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)
						sumA=sumA+prefix*ralph(iComp,iType)/(1+fAssoc*ralph(iComp,iType))
					endif
				enddo
			enddo
			ERR=fAssoc-sumA
			!checkSum=FA0/(1+fAssoc*ralphMean)
			if( ABS((ERR-errOld)/errOld) < zeroTol .and. LOUDER)write(dumpUnit,610)'MEM1: err~errOld=',ERR,errOld
			change=ERR/(ERR-errOld)*(fAssoc-fOld)
			if(LOUDER)write(dumpUnit,610)' MEM1: FA,sumA',fAssoc,sumA 
			fOld=fAssoc
			errOld=ERR
			fAssoc=fAssoc-change
			!ralphMean=(-1+FA0/sumA)/(fAssoc+1D-9)
		ENDDO
		IF(nIter > itMax-1)THEN
			if(LOUDERWert)write(dumpUnit,*)'ERROR - NO CNVRG ON fAssoc'
			iErrCode=4
			if(LOUDER)write(dumpUnit,*)
			return
		ENDIF
	endif

	!fAssoc ITERATION HAS CONCLUDED
	if(fAssoc < zeroTol)then
		if(LOUDER)write(dumpUnit,610)' MEM1: 0 > FA=',fAssoc
	endif
610 format(1x,a,8E12.4)
611 format(1x,a,i8,8E12.4)
612 format(1x,a,2i8,8E12.4)

	!Elliott'96 Eq.35 (mod for CS or ESD rdf)
	zAssoc= -fAssoc*fAssoc*dAlpha
	if(LOUDER)write(dumpUnit,610)'MEM1: eta,fAssoc,zAssoc=',eta,fAssoc,zAssoc
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
		if(LOUDER)write(dumpUnit,611)' MEM1: i,XA(i,j)=',i,(XA(i,j),j=1,nTypes(i))
	enddo

	if(LOUDER)write(dumpUnit,610)' MEM1:fAssoc,zAssoc,dAlpha,h^M,sumLnXi',fAssoc,zAssoc,dAlpha,hAD+hDA,sumLnXi
	hDA=hAD
	zAssoc= -dAlpha*(hAD+hDA)/2 ! check same?
	!aAssoc= 0 	!already initialized
	do i=1,nComps  ! XA and XD USE Assoc.
		sumLnXi=0.d0
		do j=1,nTypes(i)
			if( XA(i,j) < zeroTol)then
				if(LOUDER)write(dumpUnit,612)' MEM1: XAij < 0? i,j,XA=',i,j,XA(i,j)
				iErr=14
				return
			endif
			sumLnXi=sumLnXi+nDegree(i,j)*(  nAcceptors(i,j)*DLOG( XA(i,j) )*2  )	! ln(XD)=ln(XA) so *2.
		enddo
		aAssoc=aAssoc+xFrac(i)*sumLnXi	! ln(XD)=ln(XA) so *2.
		rLnPhiAssoc(i)=sumLnXi-(hAD+hDA)/2*(dAlpha-1)*bVolCc_mol(i)/bVolMix 	! Eq. 42
		!write(dumpUnit,611)' MEM1: iComp,sumLnXi,lnPhi=',i,sumLnXi,rLnPhiAssoc(i)
	enddo
	aAssoc=  aAssoc+(hAD+hDA)/2  ! Eqs. 1 & 40
	!aAssocPas=aAssoc
	if(LOUDER)write(dumpUnit,610)' MEM1:         fugAssocBefore=',(rLnPhiAssoc(i),i=1,nComps)
	if(LOUDER)write(dumpUnit,*)' MEM1: nTypes()=',nTypes(1:nComps)
	
	betadFA_dBeta=0	!cf E'96(19)
	do i=1,nComps !cf E'96(19)	 
		do j=1,nTypes(i)
			bepsA=(eAcceptorKcal_mol(i,j)/RgasCal*1000/tKelvin)
			dLnAlpha_dLnBeta=DEXP(bepsA)	! Eqs. 44,45
			if(bepsA > 1.D-4)dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsA)/(dLnAlpha_dLnBeta-1)
			if(LOUDER)write(dumpUnit,form612)' MEM1:i,j,eA,bepsA,BdLnAlph_dLnB',i,j,eAcceptorKcal_mol(i,j),bepsA,dLnAlpha_dLnBeta
			betadFA_dBeta=betadFA_dBeta+xFrac(i)*nDegree(i,j)*nAcceptors(i,j)*XA(i,j)*ralph(i,j)*0.5d0*dLnAlpha_dLnBeta	  
			! 0.5 because dLnRalph=0.5*dLnAlpha
		enddo
	enddo
	! Note: betadFA_dBeta here is actually betadFA_dBeta/(-half)
	uAssoc= -betadFA_dBeta*fAssoc*2   !cf E'22(Eq43): FA=FD=> -(FA*betadFD_dBeta+FD*betadFA_dBeta) = -2*FA*betadFA_dBeta
	!uAssocPas=uAssoc
	FA = fAssoc ! for common/MEM2parts/...
	FD = fAssoc ! for common/MEM2parts/...
	betadFD_dBeta=betadFA_dBeta
	xOld(1:nComps)=xFrac(1:nComps)
	IDold(1:nComps)=ID(1:nComps)
	etaOld=eta
	rdfOld=rdfContact
	if(LOUDER)write(dumpUnit,611)' MEM1:nIter,F1,fAssoc,zAssoc,rmsErr=',nIter,fAssoc1,fAssoc,zAssoc
	return
	end	!subroutine MEM1()

