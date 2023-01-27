!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
	! THE WERTHEIM ROUTINES OF SPEAD AND ESD TOGETHER COMBINED BY AGF, Oct. 2009
	subroutine Wertheim( isZiter,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErrCode)
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID ,etaMax,dumpUnit,RgasCal
	USE Assoc !XA,XD,XC
	!USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	!PARAMETER(nmx=55,maxTypes=44,maxTypesGlobal=999)
	!USE GlobConst
	character*77 errMsg(11)
    logical HyBond,bAssoc(nmx)
	DoublePrecision ralph(nmx,maxTypes),xFrac(nmx),rLnPhiAssoc(nmx),vMolecNm3(nmx) !bondVolNm3(nmx),
	DoublePrecision xaOld(nmx,maxTypes),xdOld(nmx,maxTypes) !,alphAD(nmx,maxTypes),alphDA(nmx,maxTypes)
    DoublePrecision rdfContact,dAlpha,eta
    DoublePrecision tKelvin,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij
    !DoublePrecision alphADii,alphDAii,alphCCii,dAlphADii,dAlphDAii,dAlphCCii
    DoublePrecision rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc	!dimension FAsolv(nmx,maxTypes),FDsolv(nmx,maxTypes),deltaAlphaA(nmx,maxTypes),deltaAlphaD(nmx,maxTypes) ! solvation factors: e.g. FAsolv = sum(xj*Ndj*XDj*alphADij) = (-1+1/XAi) ~ sqrt(alphaii)*Fassoc
	data initCall/1/
	LouderWert=LOUD
	!LouderWert=.TRUE.

	! NDi		= THE DEGREE OF POLYMERIZATION OF COMPO i
	! fAssoc	= THE CHARACTERISTIC ASSOCIATION = 1/ralphi(1/XAi-1)
	! ralphi	= ROOT OF ALPHA WHERE ALPHAi=RHO*VXi*KADi*Ei/(1-1.9ETA)
	! iErrCode
	!	0		= no error
	!	1		= no convergence on ss iteration.
	errMsg(11)= 'Wertheim Error: this code has not been written to handle asymmetric bonding types in a single molecule.'
	errMsg(1)= 'Wertheim Error: Nonsense in input parameters.'
	vMolecNm3(1:nComps)=bVolCc_mol(1:nComps)/avoNum


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
	if(LOUD.and.initCall)write(dumpUnit,*)'Wertheim: starting.'
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
			if(LOUD.and.initCall)write(dumpUnit,*)'Wertheim: iComp,iType,#Acc,#Don', iComp,iType,nAcceptors(iComp,iType),nDonors(iComp,iType)
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
		if(LouderWert.and.initCall)write(dumpUnit,*)'Wertheim: No HyBond. nDonTot,nAccTot=',nDonTot,nAccTot
		if(LouderWert.and.initCall)write(dumpUnit,*) 'Wertheim: no HyBond. Dont come back!'
		initCall=0
		return !no need for wertheim in this case.
	endif
	if(LouderWert)write(dumpUnit,*)'nAcceptors(iComp1)=',( nAcceptors(1,iType),iType=1,nTypes(1) )
	if(LouderWert)write(dumpUnit,*)'nDonors(iComp1)=',( nDonors(1,iType),iType=1,nTypes(1) )
	if(initCall.and.LouderWert.and.HyBond)write(dumpUnit,*) 'FuTptVtot: HyBonding identified in this fluid. Correct?'
	!write(dumpUnit,*) 'Wertheim: .not.HyBond, but still in Wertheim???'
	if(tKelvin < 2)iErrCode=1 ! < 2 is too low even for He.
	if(eta > etaMax .or. eta < 0)iErrCode=1
	if(bVolMix.le.0)iErrCode=1
	if(sumx.le.0)iErrCode=1
	if(iErrCode==1)then
		if(LouderWert)write(dumpUnit,'(a,f10.2)')' Wertheim: bVolMixCc_mol= ',bVolMix
		if(LouderWert)write(dumpUnit,'(a,1PE11.4,a,0Pf7.2,a,f7.4)')' eta=',eta,' tKelvin=',tKelvin,' sumx=',sumx
		if(LouderWert)write(dumpUnit,'(a,11F8.4)')' vMolecNm3=',(vMolecNm3(iComp),iComp=1,nComps)
		!call BeepMsg(errMsg(iErrCode))
		if(LouderWert)write(dumpUnit,*)errMsg(iErrCode)
		if(LouderWert)write(dumpUnit,*)
	endif
	if(iErrCode.ne.0)return !sorry for bad return
	
	rhoMol_cc=eta/bVolMix
	!write(dumpUnit,*) 'Calling MEM1'
	if(LouderWert)call MEM1(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fAssoc,rLnPhiAssoc,iErrMEM1 )! XA,XD USE Assoc
	if(LouderWert)write(dumpUnit,'(a,8f9.4)')' Wertheim-MEM1:eta,fAssoc,zAssoc,aAssoc,uAssoc',eta,fAssoc,zAssoc,aAssoc,uAssoc
	!write(dumpUnit,*) 'Results from MEM1'
	!if(LouderWert.and.initCall)write(dumpUnit,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
	if(LouderWert.and.initCall)write(dumpUnit,*) 'Wertheim: MEM1 estimate. Calling MEM2'
	call MEM2(isZiter,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,iErr )!,rLnPhiAssoc,ier)
	if(LouderWert)write(dumpUnit,'(a,8f9.4)')' Wertheim-MEM2:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
	!if(LouderWert)write(dumpUnit,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
	!if(LouderWert)write(dumpUnit,*) 'Wertheim: MEM2 estimate. Calling XsolJre'
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
			if(Loud)write(dumpUnit,'(a,3f9.4)')' Wertheim: eta,rdfContact(<1)=',eta,rdfContact
			if(Loud)write(dumpUnit,*)
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
				bigYhb=iComplex*(EXP(eHbKcal_mol(iComp,iType)/tKelvin/RgasCal*1000)-1) !the iComplex factor permits eHb to be non-zero, but still have zero complexation.  e.g. acetone+chloroform: nDs1=0, nAs1=1, nDs2=1, nAs2=0.  So complexation can occur between nAs1*nDs2, but not nDs1*nAs2.

				SQARG=bondVolNm3(iComp,iType)*avoNum*rhoMol_cc*rdfContact*bigYhb
				IF(SQARG < 0)THEN
					if(LouderWert)write(dumpUnit,*) 'neg alpha in wertheim. ETA=',ETA
					iErrCode=3
					if(LouderWert)write(dumpUnit,*)
					RETURN
				ENDIF
				ralph(iComp,iType)=DSQRT( SQARG )
			enddo
		enddo
		!XA=1 !set to unity means no association	          
		!XD=1 !set to unity means no association	          
		!XC=1 !set to unity means no association	          
		if(LouderWert.and.initCall) write(dumpUnit,'(a,5F10.6,3F10.1)')'Wertheim:solvation. X(1),X(2)=',xFrac(1),xFrac(2)
		!write(dumpUnit,*) 'Results from MEM1'
		!write(dumpUnit,'Wert:uAssoc,CvAssoc',uAssoc,CvAssoc
		!write(dumpUnit,'uAssoc,b_xaDxa_db',uAssoc,b_xaDxa_db

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
	!			yHBij=EXP(dHadij/tKelvin/RgasCal*1000)-1
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
		!if(LouderWert)write(dumpUnit,'(a,  5F10.6,3F10.1)')' XA(1,4),XA(2,3),XD(2,3):',XA(1,4),XA(2,3),XD(2,3) ,deltaAlphaA(1,4),deltaAlphaD(2,3),deltaAlphaD(2,3) 
        if(ierXsol.ne.0)then
            iErrCode=1
			if(LouderWert)write(dumpUnit,'(a,15i3)')' Wertheim: ier from Xsol=',ierXsol
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
				if(XA(iComp,iType) < zeroTol .and. LouderWert)write(dumpUnit,*) 'Error in Wertheim: XA<1E-11'
				if(XD(iComp,iType) < zeroTol .and. LouderWert)write(dumpUnit,*) 'Error in Wertheim: XD<1E-11'
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
		if(LOUDERwert)write(dumpUnit,'(a,8E12.4)')' Wertheim:hAD+hDA,hCC,fAssoc,zAssoc,zAcid',h_nMichelsen,hCC,fAssoc,zAssoc,zAcid

		!if(LouderWert)write(dumpUnit,'(a,8f9.4)')' From Xsol:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
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
	if(LouderWert)write(dumpUnit,'(a,8f9.4)')' Wertheim-Xsol:eta,zAssoc,aAssoc,uAssoc',eta,zAssoc,aAssoc,uAssoc
	if(LouderWert)write(dumpUnit,*) 'Wertheim: After Calling XsolJre'
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
	SUBROUTINE WertheimFugc(xFrac,nComps,eta,fugAssoc,h_nMichelsen,iErr)
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID,dumpUnit
	USE Assoc
	!USE BIPs
	IMPLICIT DoublePrecision(A-H,O-Z)
	LOGICAL LOUDER
	DoublePrecision xFrac(nmx),fugAssoc(nmx) !,vMolecNm3(nmx),NAS(nmx),NDS(nmx),ND(nmx),ralph(nmx),XAmol(nmx),XDmol(nmx)
    DoublePrecision rdfContact,dAlpha,eta
	!DIMENSION dfugAssoc_dRHO(nmx),dfugAssoc_dT(nmx) !,fugAssoc_num(nmx)
	!DIMENSION dXAdRHO(nmx,maxTypes),dXDdRHO(nmx,maxTypes),dXCdRHO(nmx,maxTypes),dXAdT(nmx,maxTypes),dXDdT(nmx,maxTypes),dXCdT(nmx,maxTypes)

	!ETAp=ETA !this is necessary for old method

	!C     SOLVING FOR THE DERIVATIVE OF g (rdf) WRT ETA
	!C     THIS IS NOT EXACTLY EQUAL TO THE DERIVATIVE
	!C     IT IS EQUAL TO eta/g TIMES THE DERIVATIVE OF g WRT eta

	LOUDER=LouderWert
	!LOUDER=.TRUE.
	iErr=0

	needDerivs=0
	if(needDerivs==0)return	! Let's forget about derivatives etc. for now and stick to VLE calculations.

	bVolMix=SUM(xFrac(1:nComps)*bVolCc_mol(1:nComps))
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
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' WertheimFugc: fugAssocBefore=',(fugAssoc(i),i=1,nComps)
	aAssoc=0
	isAcid=0
	do iComp=1,nComps		 
		sumLnXi=0.d0
		do iType=1,nTypes(iComp)
			if(XA(iComp,iType) < zeroTol .and. LouderWert)write(dumpUnit,*) 'Error in WertheimFugc: XA<zeroTol'
			if(XD(iComp,iType) < zeroTol .and. LouderWert)write(dumpUnit,*) 'Error in WertheimFugc: XD<zeroTol'
			if(XC(iComp,iType) < zeroTol .and. LouderWert)write(dumpUnit,*) 'Error in WertheimFugc: XC<zeroTol'
			if(XC(iComp,iType) < 1-zeroTol)then
				XCtemp=XC(iComp,iType)
				isAcid=1
			endif
			sumLnXi=sumLnXi+nDegree(iComp,iType)*( nAcceptors(iComp,iType)*LOG(XA(iComp,iType))+nDonors(iComp,iType)*LOG(XD(iComp,iType))+LOG(XC(iComp,iType)) )          !M&H Eq 13
		enddo
		aAssoc=aAssoc+xFrac(iComp)*sumLnXi
		!if(LOUDER)write(dumpUnit,'(a,i3,4F10.7)')' Fugc: i,XA(i,j)=',iComp,(XA(iComp,j),j=1,nTypes(iComp))
		!if(LOUDER)write(dumpUnit,'(a,i3,4F10.7)')' Fugc: i,XD(i,j)=',iComp,(XD(iComp,j),j=1,nTypes(iComp))
		!if(LOUDER)write(dumpUnit,'(a,i3,4F10.7)')' Fugc: i,XC(i,j)=',iComp,(XC(iComp,j),j=1,nTypes(iComp))
		fugAssoc(iComp)=sumLnXi-half*h_nMichelsen*(dAlpha-1)*bVolCc_mol(iComp)/bVolMix  !Eq 13
	enddo
	if(LOUDER)then
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid==1)write(dumpUnit,'(a,3E12.4)')' Fugc: XC= ',XCtemp
	endif
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' WertheimFugc: hMich,aAssoc=',h_nMichelsen,aAssoc
	if(LOUDER)write(dumpUnit,'(a,6E12.4)')' WertheimFugc: fugAssocAfter =',(fugAssoc(i),i=1,nComps)

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
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID,dumpUnit
	USE Assoc !XA,XD,XC,aBipAD,aBipDA
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	PARAMETER (iWegFreq=4)
	DoublePrecision xFrac(nmx)
	DoublePrecision fa(nmx,maxTypes),fd(nmx,maxTypes),faOld(nmx,maxTypes),fdOld(nmx,maxTypes)
	DoublePrecision xaOld(nmx,maxTypes),xdOld(nmx,maxTypes),alphCCavg(maxTypes)
	DoublePrecision sumDoni(maxTypes),sumDonj(maxTypes),sumAcci(maxTypes),sumAccj(maxTypes)
    DoublePrecision tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij
    DoublePrecision alphADii,alphDAii,alphCCii,dAlphADii,dAlphDAii,dAlphCCii
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
				sumDoni(iType)=0 !We must store these so we can update all types simultaneously on a consistent basis.
				sumAcci(iType)=0
				do jType=1,nTypes(iComp)
					if(jType.ne.iType)then
						call AlphaSp(isZiter,iComp,iType,iComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						sumDoni(iType)=sumDoni(iType)+xFrac(iComp)*nDegree(iComp,jType)*nDonors(iComp,jType)*XD(iComp,jType)*ALPHADij
						sumAcci(iType)=sumAcci(iType)+xFrac(iComp)*nDegree(iComp,jType)*nAcceptors(iComp,jType)*XA(iComp,jType)*ALPHDAij
					endif
				enddo
				sumDonj(iType)=0 ! SUMDONj =/= SUMDONi.
				sumAccj(iType)=0
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
                    sqArg=MaxVal( alphCCavg(1:nTypes(jComp)) )
					if(jComp==2)ralphC=SQRT(  sqArg  )
					if(jComp.eq.iComp)cycle	!here we treat types on different comps as different (if no intra effect, they actually should have the same XA,XD)
					do jType=1,nTypes(jComp)
						call AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
						sumDonj(iType)=sumDonj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*ALPHADij
						if(sumDonj(iType)<0)then
							if(LouderWert)write(dumpUnit,*)'iType,sumDonj(iType)',iType,sumDonj(iType)
							if(LouderWert)write(dumpUnit,*) 'XsolJre:sumdon<0?'
						endif
						sumAccj(iType)=sumAccj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*ALPHDAij
						if(sumAccj(iType)<0)then
							if(LouderWert)write(dumpUnit,*)'iType,sumAccj(iType)',iType,sumAccj(iType)
							if(LouderWert)write(dumpUnit,*) 'XsolJre:sumacc<0?'
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
					if(LouderWert)write(dumpUnit,*)'quadB,fracDA,sumAcc,sumDon',quadB,fracDA,sumAcc,sumDon
					if(LouderWert)write(dumpUnit,*) 'XsolJre: sqArg < 0.'
				endif
				if( (1+sumAcc) < 0) sumAcc= -0.99 !write(dumpUnit,*) 'XsolJre: 1+sumAcc < 0'
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
				if( (1+4*alphCCavg(iType)) < 0) write(dumpUnit,*) 'XsolJre: 1+4*alphCCavg(iType) < 0'
				xCC=2/( 1+SQRT(1+4*alphCCavg(iType)) )
				XC(iComp,iType)=xCC
			enddo  ! iType
		enddo  !  iComp
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsFerr=sqrt( rmsFerr/(2*nComps) )
		rmsOld=rmsErr
		if(LOUDER)write(dumpUnit,'(a,i3,2E12.4)')' XsolJre2a: iter,rmsErr,rmsFerr=',iter,rmsErr,rmsFerr
		!if(LOUD.and.initCall)write(dumpUnit,'(a,7F10.6,3f7.3)')' XAs&XDs,...',XA(1,2),XA(1,3),XA(2,2),XD(1,2),XD(2,2),rmsErr,wegParm,float(iter) !,deltaAlphaA(1,3) 
	enddo ! while(rmsErr > tol)
	if(LOUDER)then
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs1,...',(XA(1,i),i=1,4),(XD(1,i),i=1,4)
		write(dumpUnit,'(a,8F9.6,3f7.3)')' XAs&XDs2,...',(XA(2,i),i=1,4),(XD(2,i),i=1,4)
		if(isAcid==1)write(dumpUnit,'(a,3E12.4)')' Xsol: alphCCavg,XC= ',ralphC*ralphC,xCC
		write(dumpUnit,'(a,3E12.4)')' rmsErr,wegParm,iter',rmsErr,wegParm,float(iter) !,deltaAlphaA(1,3)
	endif 
	if(iter > itMax)then
		ierCode=2
		if(Louder)write(dumpUnit,*) 'Error in XsolJre: iter>itMax'
	else
!		if(LOUD)write(dumpUnit,*)'XsolJre converged. iter=',iter
!		if(LOUD)write(dumpUnit,*)'xa1,xa2',XA(1),XA(2)
!		if(LOUD)write(dumpUnit,*)'xd1,xd2',XD(1),XD(2)
	endif
	initCall=0
	RETURN
	END	 ! XsolJre2a

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJre1a(xFrac,rho,tKelvin,nComps,ierCode)	! Simple successive substitution
	USE Assoc !XA,XD,XC,aBipAD,aBipDA
	IMPLICIT doublePrecision(A-H,K,O-Z)
	PARAMETER (iWegFreq=99)
	DIMENSION xFrac(nmx)
	DIMENSION fa(nmx,maxTypes),fd(nmx,maxTypes),alphCCavg(maxTypes)
	DIMENSION faOld(nmx,maxTypes),fdOld(nmx,maxTypes) ,xaOld(nmx,maxTypes),xdOld(nmx,maxTypes)
	DoublePrecision sumDoni(maxTypes),sumDonj(maxTypes),sumAcci(maxTypes),sumAccj(maxTypes)
    DoublePrecision tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij
    DoublePrecision alphADii,alphDAii,alphCCii,dAlphADii,dAlphDAii,dAlphCCii
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
				sumDonj(iType)=0 ! SUMDONj =/= SUMDONi.
				sumAccj(iType)=0
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
						sumDonj(iType)=sumDonj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nDonors(jComp,jType)*XD(jComp,jType)*alphADij
						if(sumDonj(iType)<0)then
							if(LouderWert)write(dumpUnit,*)'iType,sumDonj(iType)',iType,sumDonj(iType)
							if(LouderWert)write(dumpUnit,*) 'sumdon<0?'
						endif
						sumAccj(iType)=sumAccj(iType)+xFrac(jComp)*nDegree(jComp,jType)*nAcceptors(jComp,jType)*XA(jComp,jType)*alphDAij
						if(sumAccj(iType)<0)then
							if(LouderWert)write(dumpUnit,*)'iType,sumAccj(iType)',iType,sumAccj(iType)
							if(LouderWert)write(dumpUnit,*) 'sumacc<0?'
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
				if( (1+4*alphCCavg(iType)) < 0) write(dumpUnit,*) 'Wertheim: 1+4*alphCCavg(iType)2nd < 0'
				xCC=2/( 1+SQRT(1+4*alphCCavg(iType)) )
				XC(iComp,iType)=xCC
			enddo
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsOld=rmsErr
		write(dumpUnit,'(a,6F10.6,3F10.1)')' XAs&XDs:',XA(1,2),XA(1,3),XA(2,2),XD(1,2),XD(2,2),rmsErr,float(iter), wegParm !,deltaAlphaA(1,3) 
	enddo ! while(rmsErr > tol)
	if(iter > itMax)then
		ierCode=2
		if(LouderWert)write(dumpUnit,*) 'Error in XsolJre: iter>itMax'
	else
!		if(LOUD)write(dumpUnit,*)'XsolJre converged. iter=',iter
!		if(LOUD)write(dumpUnit,*)'xa1,xa2',XA(1),XA(2)
!		if(LOUD)write(dumpUnit,*)'xd1,xd2',XD(1),XD(2)
	endif
	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Subroutine AlphaSp(isZiter,iComp,iType,jComp,jType,tKelvin,rho,bVolMix,alphADij,alphDAij,alphCCij,dAlphADij,dAlphDAij,dAlphCCij)
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID,dumpUnit,RgasCal
	USE Assoc !aBipAD,aBipDA
	USE BIPs  !H((I,J) for ESD
	implicit doublePrecision(a-h,k,o-z)
	!note: because aBip is on a type-type basis, we don't need to quadruply subscript it.
	indexi=localType( idType(iComp,iType) )	!reference to local type list so aBip matrix is small.
	indexj=localType( idType(jComp,jType) )	!reference to local type list so aBip matrix is small.
	YDi=EXP( eDonorKcal_Mol(iComp,iType)/(tKelvin*RgasCal/1000) ) - 1.d0
	YDj=EXP( eDonorKcal_Mol(jComp,jType)/(tKelvin*RgasCal/1000) ) - 1.d0
	YAi=EXP( eAcceptorKcal_Mol(iComp,iType)/(tKelvin*RgasCal/1000) ) - 1.d0
	YAj=EXP( eAcceptorKcal_Mol(jComp,jType)/(tKelvin*RgasCal/1000) ) - 1.d0
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

	yHBadij=EXP(dHadij/tKelvin/RgasCal/1000)-1.d0
	yHBdaij=EXP(dHdaij/tKelvin/RgasCal/1000)-1.d0
	KVEadij=SQRT(kVi*kVj)*yHBadij
	KVEdaij=SQRT(kVi*kVj)*yHBdaij
	dKVEadij=SQRT(kVi*kVj)*(dHadij/tKelvin/RgasCal/1000)*(yHBadij+1)   ! beta*dKVE/dBeta
	dKVEdaij=SQRT(kVi*kVj)*(dHdaij/tKelvin/RgasCal/1000)*(yHBdaij+1)
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
		if(LouderWert)write(dumpUnit,*)'alphADij,alphDAij',alphADij,alphDAij
		if(LouderWert)write(dumpUnit,*) 'AlphaSp: -ve alpha? Thats weird.'
	endif
	if(isZiter)return
	!ralph=(rootRhogK*Y)^0.5
	!dLnRalph=(beta/ralph)*0.5/ralph*dAlpha/dBeta=	0.5*dLnAlpha/dLnBeta 
	dLnAlphAi=eAcceptorKcal_mol(iComp,iType)/(RgasCal/1000*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphAi > 1.D-4)dLnAlphAi=dLnAlphAi*(YAi+1)/YAi
	dLnAlphDi=eDonorKcal_mol(iComp,iType)/(RgasCal/1000*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphDi > 1.D-4)dLnAlphDi=dLnAlphDi*(YDi+1)/YDi
	dLnAlphAj=eAcceptorKcal_mol(jComp,jType)/(RgasCal/1000*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphAj > 1.D-4)dLnAlphAj=dLnAlphAj*(YAj+1)/YAj
	dLnAlphDj=eDonorKcal_mol(jComp,jType)/(RgasCal/1000*tKelvin)	! = dLnAlph__/dLnBeta
	if(	dLnAlphDj > 1.D-4)dLnAlphDj=dLnAlphDj*(YDj+1)/YDj

	dAlphADij=iComplexAD*0.5*(ralphDj*ralphAi*dLnAlphAi+ralphAi*ralphDj*dLnAlphDj)+dKVEadij*rho*rdfContact
	dAlphDAij=iComplexDA*0.5*(ralphAj*ralphDi*dLnAlphDi+ralphDi*ralphAj*dLnAlphAj)+dKVEdaij*rho*rdfContact
	alphCCij=0
	if(idType(iComp,iType).eq.1603 .and. idType(jComp,jType).eq.1603)then 
		dHCCij=3*( eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(iComp,iType) )/2 
		!dHccij=3.00*eHbKcal_Mol(iComp,iType) !we can still vary AD and CC association ~independently by changing bondVolNm3 in ParmsHbVv.txt b/c bondVolCC=0.01=universal constant
		yHBccij=EXP(dHccij/tKelvin/RgasCal*1000)-1
		KVEccij=bondVolNm3(iComp,iType)*avoNum*yHBccij ! value for CC bonding volume taken from AD bonding volume.
		dKVEccij=bondVolNm3(iComp,iType)*avoNum*(dHccij/tKelvin/RgasCal*1000)*(yHBccij+1) !universal value for CC bonding volume
		alphCCij=KVEccij*rho*rdfContact
		dAlphCCij=dKVEccij*rho*rdfContact
	endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2009
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	! First derivatives in respect to T while N & V are constant
	dyHBadij_dT= -(dHadij/RgasCal*1000)/(tKelvin*tKelvin)*EXP(dHadij/tKelvin/RgasCal*1000)
	dyHBdaij_dT= -(dHdaij/RgasCal*1000)/(tKelvin*tKelvin)*EXP(dHdaij/tKelvin/RgasCal*1000)
	dyHBccij_dT= -(dHccij/RgasCal*1000)/(tKelvin*tKelvin)*EXP(dHccij/tKelvin/RgasCal*1000)

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
		if(rho < zeroTol .or. rdfContact < zeroTol .and. LouderWert)write(dumpUnit,'(a,2E12.4)')' AlphaSp: rho,rdfContact=',rho,rdfContact
		dalphad_dRHO=(rdfContact+rho*dg_deta*bVolMix)/(rho*rdfContact)   !Actually this is equal to RHO/alpha*dalpha_dRHO
		dalphda_dRHO=dalphad_dRHO
		dalphcc_dRHO=dalphad_dRHO
	else 
		dalphad_dRHO=0.d0
		dalphda_dRHO=0.d0
		dalphcc_dRHO=0.d0
	endif

	!if(LouderWert)write(dumpUnit,'(a,4i3,2E12.4)')' AlphaSp: i,k,j,l,alphAD,alphDA=',iComp,iType,jComp,jType,alphADij,alphDAij

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Subroutine GetAssocBips(bipHbFile,aBip,ier)  !idLocalType,nTypesTot are in USEd Assoc
	USE GlobConst, only: avoNum,zeroTol,bVolCc_mol,ID,dumpUnit
	USE Assoc
	IMPLICIT DoublePrecision(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000) ! listPool could be 10000x10000, but very few BIP corrections in reality. 
	character bipHbFile*88,dumString*123
	logical switched,LOUDER
	Integer idBinarY(listPool)
	DoublePrecision KIJDB(listPool),KTIJDB(listPool) !,KIJ(nmx,nmx),KTIJ(nmx,nmx)
	DoublePrecision aBip(maxTypes,maxTypes)!,aBipDA(nmx,nmx) !,xsTau(nmx,nmx),xsTauT(nmx,nmx),xsAlpha(nmx,nmx)
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
					IF(LOUDER)write(dumpUnit,'(a,2i5,f8.3)')' GetAssocBIPs: FOUND - idi,idj,BipIJ ',idLocalType(iType),idLocalType(jType),KijDB(item)
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
		write(dumpUnit,'(7x,11i7)')(idLocalType(iType),iType=1,nTypesTot)
		DO jType=1,nTypesTot
			write(dumpUnit,'(i7,11f7.3)')idLocalType(jType),(aBip(jType,iType),iType=1,nTypesTot)
		enddo
	endif
	return
861	continue
	if(LOUD)write(dumpUnit,*)'GetAssocBips: error opening file=',bipHbFile
	ier=1
	return 
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C

