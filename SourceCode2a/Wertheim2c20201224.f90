	! THE WERTHEIM ROUTINES OF SPEAD AND ESD TOGETHER COMBINED BY AGF, Oct. 2009

	subroutine Wertheim( vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode)
	USE Assoc !XA,XD,XC
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!PARAMETER(NMX=55,maxTypes=44,maxTypesGlobal=999)
	!USE GlobConst
	character*77 errMsg(11)
    logical HyBond,bAssoc(NMX)
	dimension vMolecNm3(NMX),xFrac(NMX),ralph(NMX,maxTypes) !bondVolNm3(NMX),
	dimension xaOld(nmx,maxTypes),xdOld(nmx,maxTypes)
	data initCall/1/

	!COMMON/ALPHA/alphAD(NMX,NMX),alphDA(NMX,NMX)
	! NDi		= THE DEGREE OF POLYMERIZATION OF COMPO i
	! fAssoc	= THE CHARACTERISTIC ASSOCIATION = 1/ralphi(1/XAi-1)
	! ralphi	= ROOT OF ALPHA WHERE ALPHAi=RHO*VXi*KADi*Ei/(1-1.9ETA)
	! iErrCode
	!	0		= no error
	!	1		= no convergence on ss iteration.
	errMsg(11)= 'Wertheim Error: this code has not been written to handle asymmetric bonding types in a single molecule.'
	errMsg(1)= 'Wertheim Error: Nonsense in input parameters.'


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
	bVolMix=0
	sumx=0
	nSitesTot=0
	nAccTot=0
	nDonTot=0
	do iComp=1,nComps
		bVolCc_mol(iComp)=vMolecNm3(iComp)*avoNum
		sumx=sumx+xFrac(iComp)
		if(bVolCc_mol(iComp).lt.1)iErrCode=1
		bVolMix=bVolMix+xFrac(iComp)*vMolecNm3(iComp)*avoNum
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
		if(LOUD.and.initCall)print*,'Wertheim: No HyBond. nDonTot,nAccTot=',nDonTot,nAccTot
		if(LOUD.and.initCall)pause 'Wertheim: no HyBond. Dont come back!'
		initCall=0
		return !no need for wertheim in this case.
	endif
	if(initCall.and.LOUD)print*,'nAcceptors(iComp1)=',( nAcceptors(1,iType),iType=1,nTypes(1) )
	if(initCall.and.LOUD)print*,'nDonors(iComp1)=',( nDonors(1,iType),iType=1,nTypes(1) )
	if(initCall.and.LOUD.and.HyBond)pause 'FuTptVtot: HyBonding identified in this fluid. Correct?'
	!pause 'Wertheim: .not.HyBond, but still in Wertheim???'
	if(tKelvin < 2)iErrCode=1 ! < 2 is too low even for He.
	if(eta > etaMax .or. eta < 0)iErrCode=1
	if(bVolMix.le.0)iErrCode=1
	if(sumx.le.0)iErrCode=1
	if(iErrCode.eq.1)then
		if(LOUD)write(*,*)' eta=',eta,' tKelvin=',tKelvin,'sumx=',sumx
		if(LOUD)write(*,*)' vMolecNm3=',(vMolecNm3(iComp),iComp=1,nComps)
		!call BeepMsg(errMsg(iErrCode))
		if(LOUD)write(*,*)errMsg(iErrCode)
		if(LOUD)pause
	endif
	if(iErrCode.ne.0)return !sorry for bad return
	
	rhoMol_cc=eta/bVolMix

	one=1
	!half=one/2 !defined in GlobConst
	Call RdfCalc(rdfContact,dAlpha,eta)
	if(rdfContact.le.1)then
		if(LOUD)write(*,'(a,3f9.4)')' Wertheim: eta,rdfContact=',eta,rdfContact
		if(LOUD)pause
		return
	endif

	
	!c1=zRefCoeff(1)+3   !todo: check vs. saft form.
	!effNumSegments=(c1-2.5)/1.5 !Note: xFrac ignores the presence of sites that do not hbond,
	!but tangent sphere segmental approach makes it clear that these sites dilute the tendency to hbond.
	!NonTangent approach of DmdTpt means we need a correlation for the effective number of tangent sphere sites.
	!Q:Why not sufficient to know # of bonding sites per volume?  ie. the dilution is by volume, regardless of vaccuum or segments in that volume.

	!Note: 1st approximation assumes equal donors and acceptors on each site.
	nSitesTot=0
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)		
			nSitesTot=nSitesTot+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)+nDonors(iComp,iType) ) !ndHb(iType) should be zero if iType does not hBond.
			iComplex=nDonors(iComp,iType)*nAcceptors(iComp,iType)
			IF( isTpt ) THEN
				eHbKcal_mol(iComp,iType)=( eDonorKcal_mol(iComp,iType)+eAcceptorKcal_mol(iComp,iType) )/2.d0 
			ENDIF		
			if(iComplex.ne.0)iComplex=1  !in general, nDs or nAs could be 2, 3,... (e.g. water).  But complexation is either true or false.
			bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/1.987e-3)-1) !the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.

			if(rdfContact.le.1)then
				if(LOUD)write(*,'(a,3f9.4)')' Wertheim: eta,rdfContact=',eta,rdfContact
				if(LOUD)pause
				iErrCode=2
				return
			endif
			SQARG=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
			IF(SQARG.LT.0)THEN
				if(LOUD)write(*,*) 'neg alpha in wertheim. ETA=',ETA
				iErrCode=3
				if(LOUD)pause
				RETURN
			ENDIF
			ralph(iComp,iType)=DSQRT( SQARG )
		enddo
	enddo

	!BEGIN SECANT ITERATION TO DETERMINE fAssoc	assuming symmetric case.
	fAssoc=0            
	SUM=0
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)      
			iComplex=nAcceptors(iComp,iType)*nDonors(iComp,iType)
			if(iComplex.ne.0)then !only apply symmetric rule when acceptors AND donors on a site.
				SUM=SUM+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType)/(1+fAssoc*ralph(iComp,iType))
				if(LOUD.and.initCall)print*,'Wertheim: iComp,iType,ralph',iComp,iType,ralph(iComp,iType)
			endif
		enddo
	enddo
	errOld=fAssoc-SUM
	fOld=fAssoc
	IF(ABS(errOld).gt.1.D-9)then
		fAssoc=1 
		NITER=0
		ERR=1111
		do while(ABS(ERR/fAssoc).GT.1.D-9.AND.NITER.LT.113)
			NITER=NITER+1                  
			SUM=0
			DO iComp=1,nComps
				do iType=1,nTypes(iComp)      
					iComplex=nAcceptors(iComp,iType)*nDonors(iComp,iType)
					if(iComplex.ne.0)then !only apply symmetric rule when acceptors AND donors on a site.
						SUM=SUM+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType)/(1+fAssoc*ralph(iComp,iType))
					endif
				enddo
			enddo
			ERR=fAssoc-SUM
			CHANGE=ERR/(ERR-errOld)*(fAssoc-fOld)
			fOld=fAssoc
			errOld=ERR
			fAssoc=fAssoc-CHANGE
		ENDDO
		IF(NITER.GE.111)THEN
			if(LOUD)write(*,*)'ERROR - NO CNVRG ON fAssoc'
			iErrCode=4
			if(LOUD)pause
			return
		ENDIF
	endif

	!fAssoc ITERATION HAS CONCLUDED

	!Elliott'96 Eq.35 (mod for CS or ESD rdf)
	ZASSOC=-fAssoc*fAssoc*dAlpha
	if(LOUD.and.initCall)print*,'Wertheim: fAssoc,zAssoc=',fAssoc,zAssoc
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)
			XA(iComp,iType)=1/(fAssoc*ralph(iComp,iType)+1)
			XD(iComp,iType)=XA(iComp,iType)
		enddo
	ENDDO
	aAssoc=0
	do iComp=1,nComps
		!aAssoc=aAssoc+xFrac(iComp)*nD(iComp)*( 2*LOG(XA(iComp))+(1-XA(iComp)) )
		do iType=1,nTypes(iComp)  !We assume here that XA=XD.
			aAssoc=aAssoc+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*( 2*LOG(XA(iComp,iType))+(1-XA(iComp,iType)) )
		enddo
	enddo
	
	sumUp=0	!cf E'96(19)
	sumDn=1	!cf E'96(19)
	sumUa=0	!cf E'96(19)
	!sumZup=0
	do iComp=1,nComps !cf E'96(19)
		do iType=1,nTypes(iComp)
			alphaTmp=ralph(iComp,iType)*ralph(iComp,iType)
			iComplex=nDonors(iComp,iType)*nAcceptors(iComp,iType)
			if(iComplex.ne.0) iComplex=1 !in general, nDs or nAs could be 2,3... (e.g.water).
			! But complexation is either true or false.
			bepsAD=eHbKcal_mol(iComp,iType)/tKelvin/1.987e-3
			bigYhb=iComplex*(EXP(bepsAD)-1)!the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.
			beps_Yhb=1  !anticipate divide by zero problem if Yhb=0
			if(bigYhb.gt.1e-4) beps_Yhb=bepsAD/bigYhb
			b_alpDalp_Db=beps_Yhb*(bigYhb+1)
			! 1-XA = alphaX^2 => -bepsDXa_dbeps = alpha*2XA*bepsDXa_dbeps + alphaX^2*b_alpDalp_Db 
			!                      => b_xaDXa_db= -(1/XA-1)*b_alpDalp_Db/(1+2XA)
			!For pure with Nd=1, b_xaDxa_Db= -XA*alpha*(beta/alpha)*(dAlpha/dBeps)/(1+2XA*alpha)=-XA*alpha/(1+2XA*alpha)*beps_Yhb*(bigYhb+1)
			!For mix, (1/XAi-1)=sum(Ndi*XDi*alphai) => -(beta/XA^2)dXA/dBeta = sum(Nd*XD*alpha*[beta/XD*dXD/dBeta+beta/alpha*dAlpha/dBeta)]
			! E'96(19)=> (beta/XA)dXA/dBeta [1+sum(Nd*alpha*XA^2)] = -XA*sum[Nd*XAj*alphaj*(beta/alpha*dAlpha/dBeta)]
			b_xaDxa_Db=-XA(iComp,iType)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*alphaTmp*b_alpDalp_Db/(1+2*XA(iComp,iType)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*alphaTmp)
			!uAssoc=uAssoc+xFrac(iComp)*nD(iComp)*b_xaDxaDb*(2-xA(iComp))
			!sumZup=sumZup-xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType)*XA(iComp,iType)*dAlpha
			sumUp=sumUp-xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*ralph(iComp,iType)*XA(iComp,iType)*b_alpDalp_Db
			sumDn=sumDn+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*alphaTmp*XA(iComp,iType)**2
			sumUa=sumUa+xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*XA(iComp,iType)*ralph(iComp,iType)*(2-XA(iComp,iType))
		enddo
	enddo
	!zAssocCk=sumZup*sumUa/sumDn   !cf E'96(19)
	uAssoc=sumUp*sumUa/sumDn   !cf E'96(Eq19)
	!CvAssoc below is for NC=1 only. JRE 20191004
	!uAssoc = -Nd*(1-XA)*bepsY1_Yad.
	iTypeAssoc=1
	do i=1,nTypes(1)
		if(ralph(1,i) > 1D-5)iTypeAssoc=i
	enddo
	XA1=XA(1,iTypeAssoc)
	Nd1=nDegree(1,iTypeAssoc)*nAcceptors(1,iTypeAssoc)
	!if(iEosOpt==5)XA1=XA(1,2) !Spead counts donors as separate types from acceptors. 
	bepsY1_Yad=1
	if(bepsAD > 1D-5)bepsY1_Yad=bepsAD*exp(bepsAD)/( exp(bepsAD)-1 ) !
	uAssocCk= -(1-XA1)*bepsY1_Yad*Nd1
	!b*dUassoc/dBeta = bepsY1_Yad*bDXa_dBeta -(1-XA)*bepsAD*{ 1*(Y+1)/Y+beps*(Y+1)/Y -[beps(Y+1)/Y^2]*(Y+1) }
	!b*dUassoc/dBeta = -bepsY1_Yad*(1-XA)*b_alpDalp_Db/(1+2XA) -(1-XA)*bepsAD*{ 1*(Y+1)/Y+beps*(Y+1)/Y -[beps(Y+1)/Y^2]*(Y+1) }
	CvAssoc= uAssoc + Nd1*bepsY1_Yad**2*XA1*(1-XA1)/(2-XA1) + Nd1*(1-XA1)*bepsY1_Yad*(1+bepsAD-bepsY1_Yad)
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

	!Check if Old guesses are better. During Z iteration, Old result may be close.
	ssqErr=0
	ssqerrOld=0
	rho=eta/bVolMix !for alphaSp()
	isAcid=0
	DO iComp=1,nComps
		do iType=1,nTypes(iComp)
			if(idType(iComp,iType).eq.1603)isAcid=1
			XC(iComp,iType)=1 !set to unity means no association	          
			sumAd=0
			sumAdOld=0
			sumDa=0
			sumDaOld=0
			DO jComp=1,nComps
				do jType=1,nTypes(jComp)
					call AlphaSp(iComp,iType,jComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
					sumAd=sumAd+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*alphADij
					sumAdOld=sumAdOld+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*xdOld(jComp,jType)*alphADij
					sumDa=sumDa+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*alphDAij
					sumDaOld=sumDaOld+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*xaOld(jComp,jType)*alphDAij
				enddo
			enddo
			ssqErr=ssqErr+( XA(iComp,iType)-1/(1+sumAd) )**2+( XD(iComp,iType)-1/(1+sumDa) )**2		 ! e.g. 1/XAi = 1+sum( xj*Ndj*XDj*alphADij ) => (1/XAi-1)/sqrt(alpha) = 1+sum(xj*Ndj*XDj*sqrtAlpha  => OK ESD96
			ssqerrOld=ssqerrOld+( xaOld(iComp,iType)-1/(1+sumAdOld) )**2+( xdOld(iComp,iType)-1/(1+sumDaOld) )**2
			if(ssqerrOld > ssqErr)then
				if(xaOld(iComp,iType) > 1e-33)XA(iComp,iType)=xaOld(iComp,iType)
				if(xdOld(iComp,iType) > 1e-33)XD(iComp,iType)=xdOld(iComp,iType)
			endif
		enddo
	enddo
	rmsErr=SQRT( ssqErr/(2*nComps) )
	if(isAcid)rmsErr=1
    iDoSolvation=1
    if(iEosOpt==2)iDoSolvation=0
    if(rmsErr < 1.D-10)iDoSolvation=0
    if(nDonTot .ne. nAccTot .and. nDonTot > 1)iDoSolvation=0  !sorry, this is just crashing.  20201221 e.g. diethyleneGlycol+ethyleneGlycol
    if(nDonTot .ne. nAccTot .and. nAccTot > 1)iDoSolvation=0
    !print*,'nDonTot,nAccTot=',nDonTot,nAccTot
    !pause 'check don acc'
  	!return !Sorry. Below disagrees slightly with above for nComps=1. Needs work. JRE 20191008.
	if(iDoSolvation==1)then  !checkup: 1E-111  ! Refine from ESD96 if necessary.  
		!if(LOUD)pause 'Wertheim: Calling XsolJre'
		call XsolJre2a(xFrac,rho,tKelvin,rdfContact,nComps,ierXsol)	 !/Assoc/ passed through common
        if(ierXsol.ne.0)then
            iErrCode=1
            return !sorry.
        endif
		!TODO: Check DIOXANE+CHLOROFORM SYSTEM WORKS LIKE E&L P541

		!C     SOLVING FOR AN INITIAL DXD GUESS USING ELLIOTT 96 EQUATION 19
		!C     NOTE THAT ALPHAij=SQRT(ALPHAi*ALPHAj) IS STILL APPLIED HERE.
		if(ierXsol)iErrCode=5 !Set error, but finish calculations anyway.

		h_nMichelsen=0           !M&H Eq 11.
		aAssoc=0				 !M&H Eq 13.
		do iComp=1,nComps		 
			do iType=1,nTypes(iComp)
				if(XA(iComp,iType) < 1.e-11 .and. LOUD)pause 'Error in WertheimFugc: XA<1E-11'
				if(XD(iComp,iType) < 1.e-11 .and. LOUD)pause 'Error in WertheimFugc: XD<1E-11'
				h_nMichelsen=h_nMichelsen+xFrac(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1-XA(iComp,iType))+nDonors(iComp,iType)*(1-XD(iComp,iType))+(1-XC(iComp,iType)) )          !M&H Eq 11.
				sumDon=nDegree(iComp,iType)*nDonors(iComp,iType)*( DLOG(XD(iComp,iType))+(1-XD(iComp,iType))/2 )
				sumAcc=nDegree(iComp,iType)*nAcceptors(iComp,iType)*( DLOG(XA(iComp,iType))+(1-XA(iComp,iType))/2 )
				sumCrb=nDegree(iComp,iType)*( DLOG(XC(iComp,iType))+(1-XC(iComp,iType))/2 )
				aAssoc=aAssoc+xFrac(iComp)*(sumAcc+sumDon+sumCrb)
				xaOld(iComp,iType)=xa(iComp,iType)
				xdOld(iComp,iType)=xd(iComp,iType)
			enddo
		enddo
		zAssoc=-half*h_nMichelsen*dAlpha	!M&H Eq 10&13.
		!print*,'Wertheim:h_nMichelsen,fAssoc',h_nMichelsen,fAssoc

		uAssocAD=0	 !cf ML Michelsen and EM Hendriks, Physical Properties from Association Models, FPE 180:165-174(2001)
		uAssocDA=0	 !use analogy to Passoc in vicinity of Eqs 6-11.
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)
				SUMDONi=nDegree(iComp,iType)*nDonors(iComp,iType)*XD(iComp,iType)
				sumAcci=nDegree(iComp,iType)*nAcceptors(iComp,iType)*XA(iComp,iType)
				DO jComp=1,nComps
					do jType=1,nTypes(jComp)
						SUMDONj=nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)
						sumAccj=nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)
						indexi=localType( idType(iComp,iType) )	!reference to local type list so aBip matrix is small.
						indexj=localType( idType(jComp,jType) )	!reference to local type list so aBip matrix is small.
						IF (iEosOpt.EQ.5.or.iEosOpt.EQ.8.or.iEosOpt.eq.9) THEN
							dHadij=( eAcceptorKcal_mol(iComp,iType)+eDonorKcal_mol(jComp,jType) )/2.0*(1-aBipAD(indexi,indexj))	!UpdatedbyEM
							dHdaij=( eDonorKcal_mol(iComp,iType)+eAcceptorKcal_mol(jComp,jType) )/2.0*(1-aBipAD(indexi,indexj)) !UpdatedbyEM
						ELSEIF (iEosOpt.EQ.4) THEN
							dHadij=( eHbKcal_mol(iComp,iType)+eHbKcal_mol(jComp,jType) )/2.0*(1-HIJ(iComp,jComp))
						ENDIF
						yHBadij=EXP(dHadij/tKelvin/1.987D-3)-1
						yHBdaij=EXP(dHdaij/tKelvin/1.987D-3)-1
						beps_Yhb=1  !anticipate divide by zero problem if Yhb=0
						if(yHBadij.gt.1e-4) beps_Yhb=dHadij/1.987e-3/(tKelvin*yHBadij)
						call AlphaSp(iComp,iType,jComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
						dAlphaAD= alphADij*beps_Yhb*(yHBadij+1)  !dAlphAD
						beps_Yhb=1  !anticipate divide by zero problem if Yhb=0
						if(yHBdaij.gt.1e-4) beps_Yhb=dHdaij/1.987e-3/(tKelvin*yHBdaij)
						dAlphaDA= alphDAij*beps_Yhb*(yHBdaij+1)  !dAlphDA 
						uAssocAD=uAssocAD+xFrac(iComp)*xFrac(jComp)*sumAcci*sumDonj*dAlphaAD
						uAssocDA=uAssocDA+xFrac(iComp)*xFrac(jComp)*sumDoni*sumAccj*dAlphaDA
						!KVEij=SQRT(kVi*kVj)*yHBij
						!iComplexAD=nAs(I)*nDs(J)
						!if(iComplexAD.ne.0)iComplexAD=1
						!iComplexDA=nDs(I)*nAs(J)
						!alphAD(I,J)=iComplexAD*KVEij*ETA/bVolMix*rdfContact
						!alphDA(I,J)=iComplexDA*KVEij*ETA/bVolMix*rdfContact
					enddo
				enddo
			enddo
		enddo
		!uAssoc= -half*(uAssocAD+uAssocDA)	!this appears to be faulty. TODO: fix. JRE 20191008			
	endif

	initCall=0
	RETURN
	END
      
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

	SUBROUTINE WertheimFugc(X,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
	USE Assoc
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(NMX),vMolecNm3(NMX),FUGASSOC(NMX) !,NAS(NMX),NDS(NMX),ND(NMX),ralph(NMX),XAmol(NMX),XDmol(NMX)
	DIMENSION dFUGASSOC_dRHO(NMX),dFUGASSOC_dT(NMX),FUGASSOC_num(NMX)
	DIMENSION dXAdRHO(NMX,maxTypes),dXDdRHO(NMX,maxTypes),dXCdRHO(NMX,maxTypes),dXAdT(NMX,maxTypes),dXDdT(NMX,maxTypes),dXCdT(NMX,maxTypes)
	COMMON/dalpha_dx/dalphad_dT,dalphda_dT,dalphcc_dT,dalphad_dRHO,dalphda_dRHO,dalphcc_dRHO
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	COMMON/num/FUGASSOC_num

	!ETAp=ETA !this is necessary for old method

	!C     SOLVING FOR THE DERIVATIVE OF g (rdf) WRT ETA
	!C     THIS IS NOT EXACTLY EQUAL TO THE DERIVATIVE
	!C     IT IS EQUAL TO eta/g TIMES THE DERIVATIVE OF g WRT eta
	bVolMix=0.d0
	do iComp=1,nComps
		bVolMix=bVolMix+X(iComp)*vMolecNm3(iComp)*avoNum
	enddo
	Call RdfCalc(rdfContact,dAlpha,eta)
	rho=eta/bVolMix
	SUMA=0.d0
	SUMD=0.d0
	SUMC=0.d0
	DO J=1,nComps
		DO jType=1,nTypes(J)
			call AlphaSp(J,jType,J,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)			
			SUMA=SUMA+X(J)*nAcceptors(J,jType)*nDegree(J,jType)*alphADij*(XA(J,jType))**2
			SUMD=SUMD+X(J)*nDonors(J,jType)*nDegree(J,jType)*alphDAij*(XD(J,jType))**2
			SUMC=SUMC+X(J)*nDegree(J,jType)*alphCCij*(XD(J,jType))**2
		enddo
	ENDDO
	one=1.d0
	!half=one/2.d0	 !defined in GlobConst 9/24/19
	vMix=0.d0
	do iComp=1,nComps
		vMix=vMix+x(iComp)*vMolecNm3(iComp)
	enddo

	call RdfCalc(rdfContact,dAlpha,eta)	!dAlpha=eta/alpha*(dAlpha/dEta)=1+eta/g*(dg/dEta)
	h_nMichelsen=0.d0           !M&H Eq 11.
	do iComp=1,nComps		 
		do iType=1,nTypes(iComp)
			h_nMichelsen=h_nMichelsen+x(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1.d0-XA(iComp,iType))+nDonors(iComp,iType)*(1.d0-XD(iComp,iType))+(1.d0-XC(iComp,iType)) )          !M&H Eq 11.
		enddo
	enddo

	do iComp=1,nComps		 
		sumLnXi=0.d0
		do iType=1,nTypes(iComp)
			if(XA(iComp,iType) < 1.e-11 .and. LOUD)pause 'Error in WertheimFugc: XA<1E-11'
			if(XD(iComp,iType) < 1.e-11 .and. LOUD)pause 'Error in WertheimFugc: XD<1E-11'
			if(XC(iComp,iType) < 1.e-11 .and. LOUD)pause 'Error in WertheimFugc: XC<1E-11'
			sumLnXi=sumLnXi+nDegree(iComp,iType)*( nAcceptors(iComp,iType)*LOG(XA(iComp,iType))+nDonors(iComp,iType)*LOG(XD(iComp,iType))+LOG(XC(iComp,iType)) )          !M&H Eq 13
		enddo
		fugAssoc(iComp)=sumLnXi-half*h_nMichelsen*(dAlpha-1)*vMolecNm3(iComp)/vMix  !Eq 13
	enddo

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Derivatives are added by AFG, Oct. 2009
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	DO I=1,nComps						  !Elliott'96 Eq.37
		DO iType=1,nTypes(I)
			dXAdT(I,iType)=0.d0
			dXDdT(I,iType)=0.d0
			dXCdT(I,iType)=0.d0
 			DO K=1,nComps 
				DO kType=1,nTypes(K)
 					! First Derivative of X in respect to T while N & V are constant
					call AlphaSp(I,iType,K,kType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
					dXAdT(I,iType)=dXAdT(I,iType)-X(K)*nDonors(K,kType)*nDegree(K,kType)*XD(K,kType)*alphADij*dalphad_dT
					dXDdT(I,iType)=dXDdT(I,iType)-X(K)*nAcceptors(K,kType)*nDegree(K,kType)*XA(K,kType)*alphDAij*dalphda_dT
					dXCdT(I,iType)=dXCdT(I,iType)-X(K)*nDegree(K,kType)*XC(K,kType)*alphCCij*dalphcc_dT
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
			h_nMichelsen=h_nMichelsen+x(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*(1.d0-XA(iComp,iType))+nDonors(iComp,iType)*(1.d0-XD(iComp,iType))+(1.d0-XC(iComp,iType)) )          !M&H Eq 11.
			dh_dT=dh_dT-x(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*dXAdT(iComp,iType)+nDonors(iComp,iType)*dXDdT(iComp,iType) ) 
 			!first derivative at constant T and N
			dh_dRHO=dh_dRHO-x(iComp)*nDegree(iComp,iType)*( nAcceptors(iComp,iType)*dXAdRHO(iComp,iType)+nDonors(iComp,iType)*dXDdRHO(iComp,iType)+dXCdRHO(iComp,iType) )
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
			dFUGASSOC_dT(iComp)=dsumLnXi_dT-half*(dalpha-1.d0)*vMolecNm3(iComp)/vMix*dh_dT
			!first derivative in respect to RHO while N & T are constant
			dFUGASSOC_dRHO(iComp)=dsumLnXi_dRHO-half*((dalpha-1.d0)*vMolecNm3(iComp)/vMix*dh_dRHO+h_nMichelsen*dAlph_deta*vMolecNm3(iComp))
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
	SUBROUTINE XsolJre2a(xFrac,rho,tKelvin,rdfContact,nComps,ierCode)
	USE Assoc !XA,XD,XC,aBipAD,aBipDA
	IMPLICIT doublePrecision(A-H,K,O-Z)
	PARAMETER (iWegFreq=3)
	DIMENSION xFrac(NMX)
	DIMENSION fa(NMX,maxTypes),fd(NMX,maxTypes),faOld(NMX,maxTypes),fdOld(NMX,maxTypes)
	DIMENSION xaOld(NMX,maxTypes),xdOld(NMX,maxTypes),alphCCavg(maxTypes)
	DoublePrecision sumDoni(maxTypes),sumDonj(maxTypes),sumAcci(maxTypes),sumAccj(maxTypes)
	WegFreq(iter)=ABS(MOD(FLOAT(iter),FLOAT(iWegFreq)))	! WegFreq=0 if iter is multiple of iWegFreq
	ierCode=0
	wegParm=0
	rmsErr=111 !avoid termination of doWhile before doing anything.
	iter=0
	do while(rmsErr.gt.1E-10.and.iter.lt.1113)
		iter=iter+1
		rmsErr=0
		DO iComp=1,nComps
			do iType=1,nTypes(iComp)
				SUMDONi(iType)=0 !We must store these so we can update all types simultaneously on a consistent basis.
				SUMACCi(iType)=0
				do jType=1,nTypes(iComp)
					if(jType.ne.iType)then
						call AlphaSp(iComp,iType,iComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
						sumDoni(iType)=sumDoni(iType)+xFrac(iComp)*nDegree(iComp,jType)*nDonors(iComp,jType)*XD(iComp,jType)*ALPHADij
						sumAcci(iType)=sumAcci(iType)+xFrac(iComp)*nDegree(iComp,jType)*nAcceptors(iComp,jType)*XA(iComp,jType)*ALPHDAij
					endif
				enddo
				SUMDONj(iType)=0
				SUMACCj(iType)=0
				alphCCavg(iType)=0
				!C	sumDon=sumDoni+sumDonj; 
				!C	sumDoni=xFrac1*sum(tj.ne.ti)[XD(1,tj)*Nd(1,tj)*nDs(1,tj)*alphaAD(1,ti;i,tj)]
				!C	sumDonj=Sum{xFracj*Sum(tj=1,nTypesj)[XD(j,tj)*Nd(j,tj)*nDs(j,tj)*alphaAD(i,ti;j,tj)]
				DO jComp=1,nComps
					do jType=1,nTypes(jComp)
						call AlphaSp(iComp,iType,jComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
						if(idType(iComp,iType).eq.1603 .and. idType(jComp,jType).eq.1603)then
							alphCCavg(iType)=alphCCavg(iType)+xFrac(jComp)*nDegree(jComp,jType)*alphCCij
						endif
					enddo
					if(jComp.eq.iComp)cycle	!here we treat types on different comps as different (if no intra effect, they actually should have the same XA,XD)
					do jType=1,nTypes(jComp)
						call AlphaSp(iComp,iType,jComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
						SUMDONj(iType)=SUMDONj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*ALPHADij
						if(SUMDONj(iType)<0)then
							if(LOUD)write(*,*)'iType,sumdonj(iType)',iType,sumdonj(iType)
							if(LOUD)pause 'sumdon<0?'
						endif
						SUMACCj(iType)=SUMACCj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*ALPHDAij
						if(SUMACCj(iType)<0)then
							if(LOUD)write(*,*)'iType,sumaccj(iType)',iType,sumaccj(iType)
							if(LOUD)pause 'sumacc<0?'
						endif
					enddo      
				enddo
			enddo !terminate the loop over all types such that sumDon's and sumAcc's are all done at once before starting next stage
			do iType=1,nTypes(iComp)
				sumDon=sumDoni(iType)+sumDonj(iType)
				sumAcc=sumAcci(iType)+sumAccj(iType)
				call AlphaSp(iComp,iType,iComp,iType,tKelvin,rho,rdfContact,bVolMix,alphADii,alphDAii,alphCCii)
				fracAD=xFrac(iComp)*nDegree(iComp,iType)*nAcceptors(iComp,iType)*alphADii
				fracDA=xFrac(iComp)*nDegree(iComp,iType)*nDonors(iComp,iType)*alphDAii
				fa(iComp,iType)=( XA(iComp,iType)-1/( 1+sumdon+fracAD*XD(iComp,iType) ) )/XA(iComp,iType)
				fd(iComp,iType)=( XD(iComp,iType)-1/( 1+sumacc+fracDA*XA(iComp,iType) ) )/XD(iComp,iType)
				!solving the quadratic between XAi and XDi to get next guess
				quadB=(1+sumDon)*(1+sumAcc)+fracAD-fracDA
				sqArg=quadB*quadB+4*fracDA*(1+sumAcc)*(1+sumDon)
				if(sqArg<0)then
					if(LOUD)write(*,*)'quadB,fracDA,sumAcc,sumDon',quadB,fracDA,sumAcc,sumDon
					if(LOUD)pause 'XsolJre: sqArg < 0.'
				endif
				wegParm= 0.0
				if(WegFreq(iter).lt.0.01)then !check for iter hit.
					dFa=fa(iComp,iType)-faOld(iComp,iType)
					dFd=fd(iComp,iType)-fdOld(iComp,iType)
					if(ABS(dFa).gt.1e-11)then
						changeAcc=fa(iComp,iType)/( dFa ) !must be inside if, or Old may not exist on first iteration.
						wegParm=changeAcc
					endif
					if(ABS(dFd).gt.1e-11)then
						changeDon=fd(iComp,iType)/( dFd ) !must be inside if, or Old may not exist on first iteration.
						wegParm=changeDon
					endif
				endif	  
				xaOld(iComp,iType)=xa(iComp,iType)
				xdOld(iComp,iType)=xd(iComp,iType)
				faOld(iComp,iType)=fa(iComp,iType)
				fdOld(iComp,iType)=fd(iComp,iType)
				if(ABS( quadB+SQRT(sqArg) ) < 1e-11)then
					xaNew=SQRT(1+sumAcc)/( (1+sumDon)*SQRT(1+sumAcc)+fracAD )
				else
					xaNew=2*(1+sumAcc)/( quadB+SQRT(sqArg) )
				endif
				if(xaNew<0)xaNew=1e-11
				XA(iComp,iType)=XA(iComp,iType)*wegParm+(1-wegParm)*xaNew
				if(XA(iComp,iType)<0)XA(iComp,iType)=1e-11	!note: WegParm>1 sometimes. Can cause overstep.
				if(XA(iComp,iType)>1)XA(iComp,iType)=1
				!XA(iComp)=XA(iComp)-(1-wegParm)*fa(iComp) !same as above
				!XD(iComp)=XD(iComp)-(1-wegParm)*fd(iComp)
				xdNew=1/( 1+sumAcc+fracDA*xaNew )  !using freshly computed XAi here makes this the analytical solution.
				if(xdNew<0)xdNew=1e-11
				XD(iComp,iType)=XD(iComp,iType)*wegParm+(1-wegParm)*xdNew
				if(XD(iComp,iType)<0)XD(iComp,iType)=1e-11	!note: WegParm>1 sometimes. Can cause overstep.
				if(XD(iComp,iType)>1)XD(iComp,iType)=1   	!note: WegParm<0 sometimes. Can cause overstep.
				!XD(iComp)=xdNew
				rmsErr=rmsErr+fa(iComp,iType)*fa(iComp,iType)
				rmsErr=rmsErr+fd(iComp,iType)*fd(iComp,iType)
				!Note: XCi=XCj b/c we assume alphCC is universal.  So we can factor out alphCCavg.
				!xCC=( -1+SQRT(1+4*alphCCavg(iType)) )/( 2*alphCCavg(iType) ) !divide by zero when type is not CC
				!xCC=( (1+4*alphCCavg(iType)-1) )/( 2*alphCCavg(iType)*(1+SQRT(1+4*alphCCavg(iType)) ) !divide by zero when type is not CC
				xCC=2/( 1+SQRT(1+4*alphCCavg(iType)) )
				XC(iComp,iType)=xCC
			enddo
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsOld=rmsErr
	enddo
	if(iter.gt.1111)then
		ierCode=2
		if(LOUD)pause 'Error in XsolJre: iter>1111'
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
	Subroutine AlphaSp(iComp,iType,jComp,jType,tKelvin,rho,rdfContact,bVolMix,alphADij,alphDAij,alphCCij)
	USE Assoc !aBipAD,aBipDA
	USE BIPs  !H((I,J) for ESD
	implicit doublePrecision(a-h,k,o-z)
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	COMMON/dalpha_dx/dalphad_dT,dalphda_dT,dalphcc_dT,dalphad_dRHO,dalphda_dRHO,dalphcc_dRHO

	!note: because aBip is on a type-type basis, we don't need to quadruply subscript it.

!	if(
    indexi=localType( idType(iComp,iType) )	!reference to local type list so aBip matrix is small.
	indexj=localType( idType(jComp,jType) )	!reference to local type list so aBip matrix is small.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Added by AFG
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	IF (iEosOpt==5 .or. iEosOpt==8 .or. iEosOpt==9) THEN
		dHadij=( eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(jComp,jType) )/2.0*( 1-aBipAD(indexi,indexj) )
		dHdaij=( eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(jComp,jType) )/2.0*( 1-aBipDA(indexi,indexj) )
	ELSE
		dHadij=( eHbKcal_Mol(iComp,iType)+eHbKcal_Mol(jComp,jType) )/2.d0*( 1-HIJ(indexi,indexj) )
		dHdaij=( eHbKcal_Mol(iComp,iType)+eHbKcal_Mol(jComp,jType) )/2.d0*( 1-HIJ(indexi,indexj) )
	ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	kVi=bondVolNm3(iComp,iType)*avoNum								   
	kVj=bondVolNm3(jComp,jType)*avoNum
	yHBadij=EXP(dHadij/tKelvin/1.987D-3)-1.d0
	yHBdaij=EXP(dHdaij/tKelvin/1.987D-3)-1.d0
	KVEadij=SQRT(kVi*kVj)*yHBadij
	KVEdaij=SQRT(kVi*kVj)*yHBdaij
	iComplexAD=nAcceptors(iComp,iType)*nDonors(jComp,jType)
	if(iComplexAD.ne.0)iComplexAD=1
	iComplexDA=nDonors(iComp,iType)*nAcceptors(jComp,jType)
	if(iComplexDA.ne.0)iComplexDA=1
	alphADij=iComplexAD*KVEadij*rho*rdfContact
	alphDAij=iComplexDA*KVEdaij*rho*rdfContact
	if(alphADij<0 .or. alphDAij<0)then
		if(LOUD)write(*,*)'alphADij,alphDAij',alphADij,alphDAij
		if(LOUD)pause 'AlphaSp: -ve alpha? Thats weird.'
	endif
	alphCCij=0
	if(idType(iComp,iType).eq.1603 .and. idType(jComp,jType).eq.1603)then 
		dHccij=3.00*eHbKcal_Mol(iComp,iType) !we can still vary AD and CC association ~independently by changing bondVolNm3 in ParmsHbVv.txt b/c bondVolCC=0.01=universal constant
		yHBccij=EXP(dHccij/tKelvin/1.987D-3)-1
		kVEccij=bondVolNm3(iComp,iType)*avoNum*yHBccij !universal value for CC bonding volume
		alphCCij=KVEccij*rho*rdfContact
	endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2009
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! First derivatives in respect to T while N & V are constant
	dyHBadij_dT=-(dHadij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHadij/tKelvin/1.987D-3)
	dyHBdaij_dT=-(dHdaij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHdaij/tKelvin/1.987D-3)
	dyHBccij_dT=-(dHccij/1.987D-3)/(tKelvin*tKelvin)*EXP(dHccij/tKelvin/1.987D-3)

	dKVEadij_dT=SQRT(kVi*kVj)*dyHBadij_dT
	dKVEdaij_dT=SQRT(kVi*kVj)*dyHBdaij_dT
	dKVEccij_dT=SQRT(kVi*kVj)*dyHBccij_dT

	if (KVEadij.ne.0)	then
		dalphad_dT=dKVEadij_dT/KVEadij   !Actually this is equal to N/alpha*dalpha_dT
	else 
		dalphad_dT=0
	endif
	if (KVEdaij.ne.0)	then
		dalphda_dT=dKVEdaij_dT/KVEdaij   !Actually this is equal to N/alpha*dalpha_dT
	else 
		dalphda_dT=0
	endif
	if (KVEccij.ne.0)	then
		dalphcc_dT=dKVEccij_dT/KVEccij   !Actually this is equal to N/alpha*dalpha_dT
	else 
		dalphcc_dT=0
	endif

	! First derivatives in respect to RHO while N & T are constant
	if (kVj.ne.0)	then
		dalphad_dRHO=(rdfContact+rho*dg_deta*bVolMix)/(rho*rdfContact)   !Actually this is equal to RHO/alpha*dalpha_dRHO
		dalphda_dRHO=dalphad_dRHO
		dalphcc_dRHO=dalphad_dRHO
	else 
		dalphad_dRHO=0.d0
		dalphda_dRHO=0.d0
		dalphcc_dRHO=0.d0
	endif

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Subroutine GetAssocBips(bipHbFile,aBip,ier)  !idLocalType,nTypesTot are in common/assoc
	USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000)
	character bipHbFile*88,dumString*123
	logical switched
	!character*123 ErrMsg(11)
	dimension idBinarY(listPool)
	doublePrecision KIJDB(listPool),KTIJDB(listPool) !,KIJ(NMX,NMX),KTIJ(NMX,NMX)
	dimension aBip(MaxTypes,MaxTypes)!,aBipDA(NMX,NMX) !,xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)
	!common/HbParms/dHkcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	ier=0
	OPEN(40,FILE=BipHbFile)
	read(40,'(a)',ERR=861)dumString !1st line is a dummy header
	item=1
	do while(item.ge.0)
		read(40,*,ERR=861,end=200) id1,id2,KIJDB(item),KTIJDB(item)
		idBinarY(item)=10000*id1+id2
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
				if(idBinarY(item).eq.idBin)then
					aBip(iType,jType)=KIJDB(item) 
					exit !found it so terminate the item loop
				endif
			enddo
			!Omit below b/c aBipAD is symmetric (???)
			!if(switched)then
			!	HTIJ(i,j)=xsTau(i,j)
			!endif
		enddo
	enddo
	return
861	continue
	if(LOUD)write(*,*)'GetAssocBips: error opening file=',bipHbFile
	ier=1
	return 
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine RdfCalc(rdfContact,dAlpha,eta)
	USE GlobConst
	implicit doublePrecision(a-h,o-z)

	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta

	! alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/dEta) = 1+(eta/rdf)*(dRdf/dEta)
	! dLng = (eta/rdf)*(dRdf/dEta)     
	! iRdfOpt	- option for characterizing the radial distribution funciton at contact.	
	!			= 0, not specified => error
	!			= 1, ESD form
	!			= 2, Carnahan-Starling form
	!			= 3, activity coefficient model form (rdf=universal constant=g(eta=0.4)
	!			= 4, ESD form, geometric association rule.
	iRdfOpt=0					!ESD non-geometric form
	if(iEosOpt==4)then
		iRdfOpt=3
	elseif(isESD)then
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
	!dLng = (eta/rdf)*(dRdf/dEta)
	denom=1.d0-1.9d0*eta 
	denom2=denom*denom
	void=1.d0-eta
	void2=void*void
	void4=void2*void2    

	dLng=1.9d0*eta/denom
	if(iRdfOpt.eq.2)dLng=eta*( 3.d0/void-1.d0/(2.d0-eta) )
	if(iRdfOpt.eq.3)dLng=-1.d0
	rdfContact=1/(1-1.9*eta)
	if(iRdfOpt.eq.2)rdfContact=(1.d0-eta/2.d0)/void/void2
	if(iRdfOpt.eq.3)rdfContact=1.d0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Added by AFG
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	d2Lng=1.9d0*eta/denom2
	d2g=d2Lng-dLng
	dAlpha=1+dLng
	dg_deta=1.9d0/denom2
	dAlph_deta=dg_deta
	if(iRdfOpt.eq.2)dg_deta=(5.d0/2.d0-eta)/void4
  	if(iRdfOpt.eq.2)dAlph_deta=3.d0/void2-2.d0/((2.d0-eta)*(2.d0-eta))
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	return
	end

