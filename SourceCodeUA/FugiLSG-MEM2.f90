!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc PURPOSE: CALCULATION OF THERMODYNAMIC PROPERTIES WITH NRTLccccccccc
	SUBROUTINE GetLsgMem2(nComps,idComp,iErrCode)
	USE GlobConst
	USE SpeadParms
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=2555)
	Character*345 bipFile,ParmsTptFile,ParmsHbFile,dumString
	Integer GetBIPs !,LSGindex(9999)
	Integer idComp(nComps)
	LOGICAL LOUDER
	!Data initCall/1/
	LOUDER=LOUD
	iErrCode=0
	ParmsTptFile=TRIM(PGLinputDir)//'\ParmsTpt.txt' ! // is the concatenation operator
	if(iEosOpt==14)ParmsTptFile=TRIM(PGLinputDir)//'\ParmsTptTransPGL6ed.txt' ! // is the concatenation operator
	if(LOUDER)write(dumpUnit,*)'ParmsTptFile=',TRIM(ParmsTptFile)
	OPEN(40,FILE=ParmsTptFile)

	READ(40,*,ioStat=ioErr)nDeck	 !check that there is something to read.
	if(ioErr.ne.0)write(dumpUnit,*)'GetLsgMem2: error opening ParmsTptFile'
	if(nDeck < 1 .and. LOUDER)write(dumpUnit,*) 'nDeck<1 in TptParms.txt.  Please check.'

	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup.
	!The basic idea is to search the hard drive until the descriptors of each component are found, then quit and rewind for next.
	!For BIPs, a linked list is used so the relevant BIPs can be accessed through a short vector that includes only nonzero BIPs
	nTypesTot=0
	do iComp=1,nComps
		iGotIt=0
		rewind(UNIT=40)
		read(40,*,ioStat=ioErr)nDeck,nTptCoeffs
		if(LOUDER)write(dumpUnit,*)'GetLsgMem2: nDeck,nTptCoeffs=',nDeck,nTptCoeffs
		jComp=0 
		DO while(jComp < nDeck.and.iGotIt==0) !rewind and loop through hard drive till you find the comp of interest.
			jComp=jComp+1
			read(40,'(a345)',ioStat=ioErr)dumString
			if(ioErr.ne.0)then
				if(LOUDER)write(dumpUnit,*)'GetLsgMem2: Error reading dumString. line=',jComp
				iErrCode=12
				return
			endif
			read(dumString,*,ioStat=ioErr)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3), &
			                             (a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),&
				vMolecNm3(iComp),tKmin(iComp),rMwDum,nTypes(iComp),(nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				bVolCc_mol(iComp)=vMolecNm3(iComp)*AvoNum
			if(ioErr.ne.0)then
				if(LOUDER)write(dumpUnit,*)'GetLsgMem2: Error reading dumString. line,dumString=',jComp,TRIM(dumString)
				iErrCode=13
				return
			endif
			IF(idBase==idComp(iComp))THEN
				iGotIt=1  !this will kick us to the next component

				if(LOUDER)write(dumpUnit,'(a,5f13.4)')' Mw,bVol(cc/mol),vEff(nm3)',rMwDum,bVolCc_mol(iComp),vMolecNm3(iComp)
				if(LOUDER)write(dumpUnit,'(a,5e13.5)')' zRefCof',(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
				if(LOUDER)write(dumpUnit,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUDER)write(dumpUnit,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUDER)write(dumpUnit,'(a,i3,9(i3,i5))')' nTypes,Nd,Id',nTypes(iComp), & 
				                                                  (nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				if(LOUDER)write(dumpUnit,'(a,i3,11i5)')' nTypes,#Fgi',nTypes(iComp),(nDegree(iComp,iType),iType=1,nTypes(iComp))
				if(nTypesTot==0)then
					nTypesTot=1
					localType( idType(1,1) )= 1 !this must happen 1st, next step will work fine even if redundant 1st time thru.
					idLocalType(1)=idType(1,1)	!these pairs point back and forth at each other
				endif
				do iType=1,nTypes(iComp) !
					isNewType=1
					do iCheck=1,nTypesTot !check if this site in this molecule is new.
						if(idType(iComp,iType).eq.idLocalType(iCheck))isNewType=0
					enddo
					if(isNewType)then
						nTypesTot=nTypesTot+1
						localType( idType(iComp,iType) )=nTypesTot !Note: nTypesTot is incrementing during this part of the code.
						idLocalType(nTypesTot)=idType(iComp,iType) !these pairs point back and forth at each other
					endif
				enddo
				etaFactor=bVolRef(iComp,Tc(iComp))/bVolCc_mol(iComp) ! tKelvin=Tc(i) here. JRE 20200406
				!etaFactor=1
				exit !quit searching if found
			ENDIF !idBase==idComp
		enddo !while(iGotIt.eq.0)
		if(iGotIt == 0)then
			iErrCode=iErrCode*10+iComp
			if(LOUDER)write(dumpUnit,*)' GetTpt: ParmsTpt.txt has no data for component #:',idComp(iComp)
			if(LOUDER)write(dumpUnit,*)
		endif
	enddo
	if(nTypesTot > maxTypes)then
		iErrCode=14
		if(LOUDER)write(dumpUnit,*)'GetLsgMem2: nTypesTot>maxTypes'
		return
	endif

	CLOSE(40)
        
	if(iErrCode.ne.0)return
		
!c  note:  bips are passed back through USEd BIPs/
	 	!phys part done.  do chem part.
	ParmsHbFile=TRIM(PGLinputDir)//'\ParmsHb4.txt' ! // is the concatenation operator
    inHbFile=40
	if(LOUDER)write(dumpUnit,*)'ParmsHbFile=',TRIM(ParmsHbFile)
	OPEN(inHbFile,FILE=ParmsHbFile)

	READ(inHbFile,*,ioStat=ioErr)nDeck	 !check that there is something to read.
	if(nDeck < 1 .and. LOUDER)write(dumpUnit,*) 'nDeck<1 in ParmsHb4.txt.  Please check.'

	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup.
	DO iComp=1,nComps	!Note: assume all types hbond. It just means summing over some zeroes.
		!nHbTypes(iComp)=nTypes(iComp) 
		bondSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType)*(nAs+nDs) ) 
		                !so no increment for nAs=nDs=0.
		numSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType) )
		do iType=1,nTypes(iComp)
			numSitesTot=numSitesTot+nDegree(iComp,iType)
			isHbType=0
			eDonorKcal_mol(iComp,iType)=0        
			eAcceptorKcal_mol(iComp,iType)=0
			!eHbKcal_mol(iComp,iType)=0  
			bondVolNm3(iComp,iType)=0 !note: bondVol is same for same type on diff components, but keeping separate for each 
			!                  component allows the potential of adding a comp specific accessibility factor at a later date.
			nDonors(iComp,iType)=0
			nAcceptors(iComp,iType)=0
			!dHkcalMol(iComp)=0
			!bondVolNm3Esd(iComp)=0
			!ND(iComp)=0
			!NDS(iComp)=0
			!NAS(iComp)=0

			iGotIt=0  !we may have valid idTypes with no hbonding. These will be set to zero hbonding.
			rewind(UNIT=inHbFile)
			read(inHbFile,*,ioStat=ioErr)nDeck
			if(ioErr.ne.0)then
				iErrCode=18
				return
			endif
			DO jType=1,nDeck
				read(inHbFile,'(a222)')dumString
				if(LOUDER)write(dumpUnit,*)TRIM(dumString)
				read(dumString,*,ioStat=ioErr)mainType,iSubType,ndsBase,nasBase,bondVolDb,bondRateDb,dHDonorKcalDb, &
				                                                                                                  dHAcceptorKcalDb
				if(ioErr.ne.0)then
					iErrCode=18
					return
				endif
				idTypeDb=mainType*100+iSubType
				IF(idTypeDb==idType(iComp,iType))THEN
					iGotIt=1
					eDonorKcal_mol(iComp,iType)=dHDonorKcalDb  
					eAcceptorKcal_mol(iComp,iType)=dHAcceptorKcalDb 
					!eHbKcal_mol(iComp,iType)=dHkcalDb 
					!nDegree(iComp,iType)=nFg(iComp,iType)
					nDonors(iComp,iType)=ndsBase
					nAcceptors(iComp,iType)=nasBase
					bondSitesTot=bondSitesTot+nDegree(iComp,iType)*(ndsBase+nasBase)
					bondVolNm3(iComp,iType)=bondVolDb
					bondRate(iComp,iType)=bondRateDb
					!dHkcalMol(iComp)=(dHDonorKcalDb+dHAcceptorKcalDb)/2.d0
					!bondVolNm3Esd(iComp)=bondVolDb
					!ND(iComp)=nDegree(iComp,iType)
					!NDS(iComp)=ndsBase
					!NAS(iComp)=nasBase
					exit !terminate the search for this type because you found it.
				ENDIF
			enddo !loop over ParmsHb4.txt
			!if(iGotIt.eq.0)iErrCode=iErrCode*10+iComp
		enddo !loop over iType
		iBondExp=2
		!bondVolNm3Esd(iComp)=bondVolNm3Esd(iComp)*( 1+bondRate(iComp,1)*(bondSitesTot/numSitesTot)**iBondExp )
		!above is for a single bSite/molec.  Below is more general, but not implemented yet.
		do iHbType=1,nTypes(iComp)
			screenFactor=( 1+bondRate(iComp,iHbType)*( nDegree(iComp,iHbType)/DFLOAT(numSitesTot) )**iBondExp )
			bondVolNm3(iComp,iHbType)=bondVolNm3(iComp,iHbType)*screenFactor
            if(bondVolNm3(iComp,iHbType) < 0)bondVolNm3(iComp,iHbType)=2D-8 !this only happens for formic acid. JRE 20200303.
		enddo
		if(LOUDER)write(dumpUnit,*)' SiteType:',( idType(iComp,i),i=1,nTypes(iComp))
		if(LOUDER)write(dumpUnit,'(a,11i4)')' nDegree   :',( nDegree(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(dumpUnit,'(a,11i4)')' nAcceptors:',( nAcceptors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(dumpUnit,'(a,11i4)')' nDonors:   ',( nDonors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(dumpUnit,'(a,9E11.4)')' bondVolNm3:',( bondVolNm3(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(dumpUnit,'(a,11f7.1)')' Don energy:',( eDonorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(dumpUnit,'(a,11f7.1)')' Acc energy:',( eAcceptorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
	enddo !iComp. All compounds have their parameters defined.
	CLOSE(inHbFile)

	bipFile=TRIM(PGLinputDir)//'\BipLsgMem2.txt' ! // is the concatenation operator
    iErrCode=GetBIPs(bipFile,idComp,nComps)	 

	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)
	if(LOUDER)then
		write(dumpUnit,*)'       aij'
		do i=1,nComps
			write(dumpUnit,form600)(xsTau(i,j),j=1,nComps)
		enddo
	endif

	return                      
	END
	!PROGRAMED BY AV 06/22
	!PURPOSE: CALCULATION OF VAPOR PRESSURES AND FUGACITY COEFFICIENTS BY NRTL ACTIVITY COEFFICIENT MODEL

	Subroutine FuLsgMEM2(tKelvin,pMpa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr)
!c  use NRTL style of EXPression
	USE GlobConst ! includes vLiq when you call CritParms first.
	USE BIPs
	USE VpDb ! for VpCoeffs
	Implicit DoublePrecision(A-H,K,O-Z)
	Parameter(CoordNum=6)
	DoublePrecision xFrac(Nmx)
	Integer kComp
	DoublePrecision fugc(nmx),xsGamma(nmx),volFrac(nmx),vSat(nmx),pSat(nmx)
	DoublePrecision omega(Nmx,Nmx),sumDenom(Nmx),xInf(nmx),GpureAssoc(nmx),rLnGamAssoc(nmx) !,sumNumer(Nmx)
	Logical LOUDER
	Common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),Zrep,Zatt,Zassoc,aRep,aAtt,aAssoc
	iErr=0

	NC=nComps
	LOUDER=LOUD
	!LOUDER=.TRUE.
	zFactor=1
	fugc(1:NC)=0
	aRes=0
	uRes=0
	zRep=0	 ! so, e.g., Zrep*bVol(i)/bMix=0
	zAtt=0
	aRep=0
	aAtt=0
	if(LIQ==0)return ! ideal gas for vapor means... That's all folks!
	expo=2.d0/7.d0
	vLiqMix=0
	if(LOUDER)write(dumpUnit,form610)' FuLsgMem2: T,P,x1,Tc1,Tc2=',tKelvin,pMPa,xFrac(1),Tc(1:NC)
	do iComp = 1,nComps
		phi=0
		if( tKelvin < Tc(iComp) )phi=(1-tKelvin/Tc(iComp))**expo	!Zc = PcVc/RTc => Zc*RTc/Pc = Vc => Assume vSat=Vc if T>Tc.
		vSat(iComp)=Rgas*Tc(iComp)/Pc(iComp)*( 0.29056+0.08775*acen(iComp) )**(1+phi) !PGL6ed Eq. 5.2-2.
		if(Tc(iComp) > 300)then
			phi=phi - (1-298/Tc(iComp))**expo  !PGL6ed Eq. 5.2-6, aka Yamada-Gunn adaptation of Rackett eq.
			vSat(iComp)=vLiq(iComp)*( 0.29056+0.08775*acen(iComp) )**phi
		endif
		vLiqMix=vLiqMix+xFrac(iComp)*vSat(iComp)
	enddo
	rhoMol_cc=1/vLiqMix
	bMix=SUM( xFrac(1:NC)*bVolCC_mol(1:NC) )
	volFrac(1:NC)=xFrac(1:NC)*bVolCC_mol(1:NC)/bMix
	eta=rhoMol_cc*bMix
	zFactor=pMpa*vLiqMix/(rGas*tKelvin)
	xInf(1:NC)=zeroTol
	do iComp = 1,nComps
		xInf(iComp)=1-(NC-1)*zeroTol
		rhoPure=1/vSat(iComp)
		Call MEM2(0,tKelvin,xInf,nComps,rhoPure,zAssoc,aAssoc,uAssoc,fugAssoc,iErrF )
		GpureAssoc(iComp)=fugAssoc(iComp)-zAssoc
		xInf(iComp)=zeroTol
	enddo
	if(LOUDER)write(dumpUnit,form610)' FuLsgMem2: GpureAssoc=',GpureAssoc(1:NC)

	sumLog=0
	do jComp=1,nComps
		sumDenom(jComp)=0
		!sumNumer(jComp)=0
		do kComp=1,nComps
			omega(kComp,jComp)=EXP( -xsTau(kComp,jComp)/tKelvin )
			sumDenom(jComp)=sumDenom(jComp)+volFrac(kComp)*omega(kComp,jComp)  !denom1=phi1*Om11+phi2*Om2,1
			!sumNumer(jComp)=sumNumer(jComp)+xFrac(kComp)*omega(kComp,jComp)*xsTau(kComp,jComp)/tKelvin
		enddo
		sumLog=sumLog+xFrac(jComp)*sumDenom(jComp)
	enddo

	Call MEM2(0,tKelvin,xFrac,nComps,rhoMol_cc,zAssoc,aAssoc,uAssoc,fugAssoc,iErrF )
	CALL GetVp(nComps,ID,iErrVp)
	if(iErrVp.ne.0)iErrCode=16
	xsFreeEn =0
	GmixIdSoln=0
	GeFlory=0
	do kComp=1,nComps
		bigSumK=0.d0
		do jComp=1,nComps
			bigSumK=bigSumK+volFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)  ! = phi1*Om11/denom1+phi2*Om2,1/denom2
		enddo
		GeFlory =GeFlory+LOG(volFrac(kComp)/xFrac(kComp))*xFrac(kComp)
		GmixIdSoln =GmixIdSoln+LOG(xFrac(kComp))*xFrac(kComp)
		qFactor=2*eta*CoordNum/2*bVolCC_mol(kComp)/(0.01484*avoNum)
		fugAtt(kComp)=qFactor*( 1-LOG(sumDenom(kComp))-bigSumK )	! assuming qi~bVol(i)/bVol(water)
		fugRep(kComp)=LOG(volFrac(kComp)/xFrac(kComp))+(1-volFrac(kComp)/xFrac(kComp)) ! FHA contributoin
		rLnGamAssoc(kComp)=fugAssoc(kComp)-zAssoc*bVolCc_mol(kComp)/bMix-GpureAssoc(kComp) 
		fugc(kComp)=( fugRep(kComp)+fugAtt(kComp) )+rLnGamAssoc(kComp)
		if(fugc(kComp) < -25) fugc(kComp)= -25		! Avoid gamma < zeroTol. Only happens for bad guesses of xsTau(i,j)
		xsGamma(kComp)=exp( fugc(kComp) )
		xsFreeEn =xsFreeEn+fugc(kComp)*xFrac(kComp)	! = SUM(xi*lnGami)
		if(tKelvin < Tc(kComp))then
			pSat(kComp)=exp(vpCoeffs(kComp,1)+vpCoeffs(kComp,2)/tKelvin+vpCoeffs(kComp,3) &
			                                                  *DLOG(tKelvin)+vpCoeffs(kComp,4)*tKelvin**vpCoeffs(kComp,5))/1000000
		else
			pSat(kComp)=Pc(kComp)*10**( 7*(1+acen(kComp))/3*(1-tKelvin/Tc(kComp)) - 3*exp(-3*Tc(kComp)/tKelvin) )
		endif
		fugc(kComp)=fugc(kComp)+LOG(pSat(kComp)/pMPa)	! = ln(gami*PiSat/P)
		if(LOUDER)write(dumpUnit,form610)' FuLsgMem2: Psat=',Psat
	enddo
	pTot = SUM( xFrac(1:NC)*xsGamma(1:NC)*pSat(1:NC) )	! = SUM( xi*P*phi(i) ) where phi(i)=exp(fugc(i))=gami*piSat/P
	if(LOUDER)write(dumpUnit,612)xFrac(1),eta, vLiqMix, (fugc(i),i=1,NC),(fugRep(i),i=1,NC),(fugAtt(i),i=1,NC),rLnGamAssoc(1:NC)
612	format(2(f7.4,1x),(f7.2,1x),12(f7.3,1x))    
	aRes=GmixIdSoln+xsFreeEn !assuming Vxs=0
	if(LOUD)then
		write(dumpUnit,form610)' FuLsgMem2: x1,vFrac1,eta,bMix,vMix',xFrac(1),volFrac(1),eta,bMix,vLiqMix
		write(dumpUnit,*)'       aij'
		do i=1,nComps
			write(dumpUnit,form600)(xsTau(i,j),j=1,nComps)
		enddo
		write(dumpUnit,form610)' FuLsgMem2: sumDenom=',sumDenom(1:NC)
		write(dumpUnit,form610)' FuLsgMem2: ln(gamAtt)=',fugAtt(1:NC)
		write(dumpUnit,form610)' FuLsgMem2: fugAssoc=',fugAssoc(1:NC)
		write(dumpUnit,form610)' FuLsgMem2: xsFreeEn,gam()=',xsFreeEn,xsGamma(1:NC)
	endif
	return
	end
