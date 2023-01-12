!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SpeadParms
	USE GlobConst, only:NMX,ID,Tc,CvRes_R,etaMax,Rgas,zeroTol,etaPure,bVolCC_mol,LOUD !need most of GlobConst here  
	USE Assoc !need all of assoc so USE SpeadParms automatically includes Assoc
	DoublePrecision zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX) !,rMwPlus(NMX)
	DoublePrecision ks0,ks1	      ! Parameters to estimate Sxs of Vahid and Elliott (2014)
	!Parameter(ks0= -0.04d0,ks1=0) ! best I could do based on avg of the entire system AV 11/29/08
	Integer nTptCoeffs
	LOGICAL bGaussEx
	DoublePrecision GexG0(NMX),GexG1(NMX),GexEta0(NMX) ! Gaussian Extrapolation parameters
	LOGICAL isHeNeH2(NMX) ! flag to indicate ultracryo compounds.
	!nnx = Total number of united-atom sites in a molecule.
END MODULE SpeadParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetTptCas(nComps,idCasPas,iErrCode,errMsgPas)
	!  Purpose:	load up the parms for tpt, including HbParms. We only load hbonding parameters for types that actually hbond.  This makes the 
	!			summations and arrays more compact in the Wertheim functions.
	!  
	!  PROGRAMMED BY:  JRE 2/02
	!  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	! 
	!  LOOKS UP THE TPT PROPERTIES AND RETURNS them in commons
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	USE GlobConst,  only: PGLinputDir,idCas,ID
	USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000)
	character*77 errMsgPas
	!	integer GetBIPs
	integer idComp(nComps),idCasPas(nComps) !localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	idCas(1:nComps)=idCasPas(1:nComps)
	call IdDipprLookup(nComps,idCasPas,ier,errMsgPas) ! ID is USEs GlobConst
	idComp(1:nComps)=ID(1:nComps)
	call GetTpt(nComps,idComp,iErrCode,errMsgPas)

	return
	end



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetTpt(nComps,idComp,ierrCode,errMsgPas)
	!  Purpose:	load up the parms for tpt, including HbParms. We only load hbonding parameters for types that actually hbond.  This makes the 
	!			summations and arrays more compact in the Wertheim functions.
	!  
	!  PROGRAMMED BY:  JRE 2/02
	!  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	! 
	!  LOOKS UP THE TPT PROPERTIES AND RETURNS them in commons
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	USE GlobConst
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	USE VpDb	   ! For iEosOpt=8, SpeadGamma
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000)
	integer GetBIPs
	integer idComp(nComps) !,ier(12) !localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!integer localType(maxTypesGlobal),idLocalType(maxTypes)  !idLocalType points back to localType for double-linking
	doublePrecision bondRate(nmx,maxTypes) !local to this subroutine.
	!DoublePrecision epskAdKelvin(maxTypes,maxTypes),epskDaKelvin(maxTypes,maxTypes) !local to this subroutine, just for info.
	DoublePrecision gmol(NMX),chemPo(NMX)
	!dimension zRefDb(5),a1Db(5),a2Db(5)
	character*88 bipFile*88,bipHbFile*88,ParmsHbFile*88,ParmsTptFile*88
	character*333 dumString*333
	character*77 ErrMsg(0:11),errMsgPas
	character*234 outFile
	logical LOUDER
	!integer nFg(nmx,maxTypes) !not in common because we transcribe these into HbParms and they do not need to be passed otherwise.
	!common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!common/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	!COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
    !common/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!common/siteParm/nFg	   !This is needed in GetRG
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	LOUDER=LOUD
	!LOUDER=.TRUE. !COMMENT OUT FOR RELEASE
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isESD or isPcSaft.
	isTPT=.TRUE. ! in GlobConst simplifies calls in FUGI and FuVtot
    etaMax = 0.99D0 ! if eta > etaMax, errors will be declared
	bGaussEx=.FALSE.
	if(iEosOpt==20)bGaussEx=.TRUE.
	iErrCode=0
	ErrMsg(0)='GetTpt: Success!'
	ErrMsg(1)='GetTpt: error reading ParmsTpt.txt. Path? Debug?'
	ErrMsg(2)='GetTpt: error reading ParmsHb4.txt. Path?'
	ErrMsg(3)='GetTpt: error reading ParmsTpt.txt at line='
	ErrMsg(4)='GetTpt: nTypesTot > maxTypes'
	ErrMsg(5)='GetTpt: a1 > 0 for at least one component when 0 < eta < 0.85'
	ErrMsg(6)='GetTpt: a2 > 0 for at least one component when 0 < eta < 0.85'
		ParmsTptFile=TRIM(PGLinputDir)//'\ParmsTpt.txt' ! // is the concatenation operator
		if(iEosOpt==14)ParmsTptFile=TRIM(PGLinputDir)//'\ParmsTptTransPGL6ed.txt' ! // is the concatenation operator
	if(LOUDER)write(*,*)'ParmsTptFile=',TRIM(ParmsTptFile)
	OPEN(40,FILE=ParmsTptFile,ERR=861)

	READ(40,*,ERR=861)nDeck	 !check that there is something to read.
    if(nDeck < 1 .and. LOUDER)pause 'nDeck<1 in TptParms.txt.  Please check.'

	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup.
    !The basic idea is to search the hard drive until the descriptors of each component are found, then quit and rewind for the next.
    !For BIPs, a linked list is created so the relevant BIPs can be accessed through a short vector that includes just the nonzero BIPs(?).
	nTypesTot=0
	do iComp=1,nComps
		iGotIt=0
		rewind(UNIT=40,ERR=861)
		read(40,*,ERR=861)nDeck,nTptCoeffs
		jComp=0 
		DO while(jComp.lt.nDeck.and.iGotIt.eq.0) !rewind and loop through hard drive till you find the comp of interest.
			jComp=jComp+1
			read(40,'(A333)',ioStat=ioErr,END=111)dumString
			if(ioErr.ne.0)goto 863
			read(dumString,*,ioStat=ioErr)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),&
				vMolecNm3(iComp),tKmin(iComp),rMwDum,nTypes(iComp),(nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				bVolCC_mol(iComp)=vMolecNm3(iComp)*AvoNum
			if(ioErr.ne.0)goto 863
111			continue  !transfer from read(40,...,END=111)dumString
			!read(40,*,ERR=861)idBase,idCc,(zRefDb(iCoeff),iCoeff=1,3),(a1Db(iCoeff),iCoeff=1,nTptCoeffs),(a2Db(iCoeff),iCoeff=1,nTptCoeffs),vMolecDb,tKmin(iComp),nTypes(iComp),(idType(iComp,iType),iType=1,nTypes(iComp)),(nFg(iComp,iType),iType=1,nTypes(iComp))
			IF(idBase.EQ.idComp(iComp))THEN
				iGotIt=1  !this will kick us to the next component

				if(LOUDER)write(*,'(a,5f13.4)')' Mw,bVol(cc/mol),vEff(nm3)',rMwDum,bVolCc_mol(iComp),vMolecNm3(iComp)
				if(LOUDER)write(*,'(a,5e13.5)')' zRefCof',(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
				if(LOUDER)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUDER)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUDER)write(*,'(a,i3,9(i3,i5))')' nTypes,Nd,Id',nTypes(iComp),(nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				if(LOUDER)write(*,'(a,i3,<nTypes(iComp)>i5)')' nTypes,#Fgi',nTypes(iComp),(nDegree(iComp,iType),iType=1,nTypes(iComp))
				if(nTypesTot.eq.0)then
					nTypesTot=1
					localType( idType(1,1) )= 1 !this must happen 1st, the next step will work fine even if redundant 1st time thru.
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
                etaFactor=bVolRef(iComp,Tc(iComp))/bVolCC_mol(iComp) ! tKelvin=Tc(i) here. JRE 20200406
				!etaFactor=1
				if(LOUDER)write(*,*)' bVol(Tc)/bVolCC_mol =  ' ,etaFactor
                if(LOUDER)write(*,*)'  eta     zRef     A1/100    A2/10000   rho(g/cc)   zRefCs'
				etaStd=0.4d0
				increment=5
                do i=5,85,increment
                    eta=i/1.D2
                    etaRef=eta*etaFactor
					rhoG_cc=eta*rMw(iComp)/bVolCC_mol(iComp)
					ZrefCs = 4*eta*(1-eta/2)/(1-eta)**3
	                call TptTerms(isZiter,nComps,iComp,eta,etaStd,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
                    if(LOUDER)write(*,'(f7.4,6f10.5)')eta,z0i,a1i/100,a2i/(100*100),rhoG_cc,zRefCs
                    if(a1i > 0)iErrCode=5
                    if(a2i > 0)iErrCode=6
					if(iErrCode > 0 .and. eta < etaMax)etaMax=eta-increment/1.D2 !We might still get useful results under reasonable conditions.
                enddo
                if(iErrCode .and. LOUDER)write(*,*) 'GetTpt: check Ais. A1 or A2 > 0 ? iErrCode=',iErrCode
				if(iErrCode==5 .or. iErrCode==6)then
					iErrCode=0 !convert to warning, while reducing etaMax to accommodate.
					if(LOUDER)print*,'etaMax reduced to:',etaMax
				endif
				exit !quit searching if found
			ENDIF !idBase==idComp
		enddo !while(iGotIt.eq.0)
		if(iGotIt == 0)then
			iErrCode=iErrCode*10+iComp
			if(LOUDER)write(*,*)' GetTpt: ParmsTpt.txt has no data for component #:',idComp(iComp)
			if(LOUDER)pause
		endif
	enddo
	if(nTypesTot > maxTypes)then
		iErrCode=4
		if(LOUDER)write(*,*)ErrMsg(4)
		if(LOUDER)pause
	endif
	CLOSE(40)
        
	if(iErrCode.ne.0)return

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG, 2011
!C This function replaces the SPEAD parameters with those of SPEAD11 correlation
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	if (iEosOpt.Eq.9) then
		CALL swap(nComps,nDegree)
	endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
 	!ref part done.  do att part.
				
		ParmsHbFile=TRIM(PGLinputDir)//'\ParmsHb4.txt' ! // is the concatenation operator
    inHbFile=40
	if(LOUDER)write(*,*)'ParmsHbFile=',TRIM(ParmsHbFile)
	OPEN(inHbFile,FILE=ParmsHbFile)

	READ(inHbFile,*,ERR=862)nDeck	 !check that there is something to read.
	if(nDeck < 1 .and. LOUDER)pause 'nDeck<1 in ParmsHb4.txt.  Please check.'

	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup.
	DO iComp=1,nComps	!Note: assume all types hbond. It just means summing over some zeroes.
		!nHbTypes(iComp)=nTypes(iComp) 
		bondSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType)*(nAs+nDs) ) !so no increment for nAs=nDs=0.
		numSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType) )
		do iType=1,nTypes(iComp)
			numSitesTot=numSitesTot+nDegree(iComp,iType)
			isHbType=0
			eDonorKcal_mol(iComp,iType)=0        
			eAcceptorKcal_mol(iComp,iType)=0
			!eHbKcal_mol(iComp,iType)=0  
			bondVolNm3(iComp,iType)=0 !note: bondVol is same for same type on diff components, but keeping separate for each component allows the potential of adding a comp specific accessibility factor at a later date.
			nDonors(iComp,iType)=0
			nAcceptors(iComp,iType)=0
			!dHkcalMol(iComp)=0
			!bondVolNm3Esd(iComp)=0
			!ND(iComp)=0
			!NDS(iComp)=0
			!NAS(iComp)=0

			iGotIt=0  !we may have valid idTypes with no hbonding. These will be set to zero hbonding.
			rewind(UNIT=40)
			read(inHbFile,*,ERR=862)nDeck
			DO jType=1,nDeck
				read(inHbFile,*,ERR=862)mainType,iSubType,ndsBase,nasBase,bondVolDb,bondRateDb,dHDonorKcalDb,dHAcceptorKcalDb
				idTypeDb=mainType*100+iSubType
				IF(idTypeDb.EQ.idType(iComp,iType))THEN
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
			bondVolNm3(iComp,iHbType)=bondVolNm3(iComp,iHbType)*( 1+bondRate(iComp,iHbType)*( nDegree(iComp,iHbType)/DFLOAT(numSitesTot) )**iBondExp )
            if(bondVolNm3(iComp,iHbType) < 0)bondVolNm3(iComp,iHbType)=2D-8 !this only happens for formic acid. JRE 20200303.
		enddo
		if(LOUDER)write(*,*)' SiteType:',( idType(iComp,i),i=1,nTypes(iComp))
		if(LOUDER)write(*,'(a,11i4)')' nDegree   :',( nDegree(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(*,'(a,11i4)')' nAcceptors:',( nAcceptors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(*,'(a,11i4)')' nDonors:   ',( nDonors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(*,'(a,9E11.4)')' bondVolNm3:',( bondVolNm3(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(*,'(a,11f7.1)')' Don energy:',( eDonorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUDER)write(*,'(a,11f7.1)')' Acc energy:',( eAcceptorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
		gmol=0 ! Q&D estimate of eta where P~0
		gmol(iComp)=1
		isHeNeH2(iComp)=.FALSE.
		if(ID(iComp)==913.or.ID(iComp)==923.or.ID(iComp)==919.or.ID(iComp)==902.or.ID(iComp)==925)isHeNeH2(iComp)=.TRUE.
		if(tKmin(iComp) < 50 .and. .not.isHeNeH2(iComp) )tKmin(iComp)=0.4*Tc(iComp) 
		vTotCc=bVolCc_mol(iComp)/0.5D0
		if(LOUDER)print*,'Calling Fu: Tmin,vTot=',tKmin(iComp),vTotCc
		isZiter=1 
		call FuTptVtot(isZiter,zFactorLo,aRes,uRes,chemPo,vTotCc,tKmin(iComp),gmol,nComps,iErrFu)
		if(iErrFu > 0 .and. LOUDER)pause 'GetTpt: error from FuTptVtot at eta=0.5' 
		vTotCc=bVolCc_mol(iComp)/0.7D0
		if(LOUDER)print*,'Calling Fu: Tmin,vTot=',tKmin(iComp),vTotCc
		call FuTptVtot(isZiter,zFactorHi,aRes,uRes,chemPo,vTotCc,tKmin(iComp),gmol,nComps,iErrFu)
		if(iErrFu > 0 .and. LOUDER)pause 'GetTpt: error from FuTptVtot at eta=0.7'
		if(iErrFu > 0)etaPure=0.5d0
		!Interpolate on V, assuming that V is more linear in ZsatLiq than rho is. 
		!v=vLo+(Ztarget-ZLo)*dV/dZ = vLo-ZLo*(vHi-vLo)/(Zhi-ZLo)
		!vTot0=bVolCc_mol(iComp)*(  1/0.5d0-zFactorLo*(1/0.7d0-1/0.5d0)/(zFactorHi-zFactorLo)  )
		!etaPure=bVolCc_mol/vTot0
		etaPure(iComp)=1/(  1/0.5d0-zFactorLo*(1/0.7d0-1/0.5d0)/(zFactorHi-zFactorLo)  )
		if(LOUDER)write(*,'(a,f10.5)')' etaPure(iComp)~',etaPure(iComp) 
        if(LOUDER)pause 'Check HB parms&etaPure.'
	enddo !iComp. All compounds have their parameters defined.
	CLOSE(inHbFile)

	nC=nComps

!  note:  bips USE BIPs
		bipFile=TRIM(PGLinputDir)//'\BipSpead.txt' ! // is the concatenation operator
	iErrCode=GetBIPs(bipFile,ID,nC)
	!Note: no need to check for switching because kij and ktij are symmetric
	if(iErrCode > 10)then
		iErr=11 ! 
	elseif(iErrCode.ne.0 .and. LOUDER)then
		pause 'GetTpt Error: GetBIPs returned error.'
	endif
    if(LOUDER)then
	    if(nC > 1)pause 'GetTpt: you need to reincorporate GetBips for mixtures.'
		if(iErrCode > 10) print*,'GetTpt: BIPs missing for ',iErrCode-10,' binary combinations.'
		print*,'GetTpt: bipFile=',TRIM(bipFile)
		if(iErrCode.ne.0) print*,'GetTpt: bipFile read error=',iErrCode
    end if

	!load aBipAD matrix
		bipHbFile=TRIM(PGLinputDir)//'\BipDA.txt' ! // is the concatenation operator
	!if(LOUDER)print*,'GetTpt: Only BipDA is needed.'
	call GetAssocBips(bipHbFile,aBipDA,ierABip) !in WertheimVv.f90. idLocalType,nTypesTot USE Assoc 
	do i=1,nTypesTot
		do j=1,nTypesTot
			aBipAD(j,i)=aBipDA(i,j)
		enddo
	enddo

	!load aBipDA matrix

	if(LOUDER)THEN ! DISPLAY RELEVANT kijAD
		write(*,*)'bipHbFile=',TRIM(bipHbFile)
		write(*,*)' Solvation BIPs'
		write(*,'(5x,11i7)')(idLocalType(iType),iType=1,nTypesTot)
		do jType=1,nTypesTot
			write(*,'(i5,11f7.3)')idLocalType(jType),( aBipAD(jType,iType),iType=1,nTypesTot )
		enddo
		do iComp=1,nComps
			do jComp=1,nComps
				do iType=1,nTypes(iComp)
					do jType=1,nTypes(jComp)
						indexi=localType( idType(iComp,iType) )	!reference to local type list so aBip matrix is small.
						indexj=localType( idType(jComp,jType) )	!reference to local type list so aBip matrix is small.
						!iComplexDA=0
						!if(nDonors(iComp,iType)*nAcceptors(jComp,jType) > 0)iComplexDA=1
						!iComplexAD=0
						!if(nDonors(jComp,jType)*nAcceptors(iComp,iType) > 0)iComplexAD=1										 !if aBip < 0 then energy is made more positive (stronger hbonding energy).
						!epskDaKelvin(indexi,indexj)=iComplexDA*( eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(jComp,jType) )/2-( eDonorKcal_Mol(iComp,iType)+eAcceptorKcal_Mol(jComp,jType) )/2*(aBipDA(indexi,indexj) )
						!epskAdKelvin(indexi,indexj)=iComplexAD*( eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(jComp,jType) )/2-( eAcceptorKcal_Mol(iComp,iType)+eDonorKcal_Mol(jComp,jType) )/2*(aBipAD(indexi,indexj) )
					enddo
				enddo
			enddo
		enddo
!		write(*,*)' AiDj Energies (kCal/mol)'
!		write(*,'(a5,11i7)')' j\i ',(idLocalType(iType),iType=1,nTypesTot)
!		do jType=1,nTypesTot
!			write(*,'(i5,11f7.2)')idLocalType(jType),( epskAdKelvin(jType,iType),iType=1,nTypesTot )
!		enddo
!		write(*,*)' DiAj Energies (kCal/mol)'
!		write(*,'(a5,11i7)')' j\i ',(idLocalType(iType),iType=1,nTypesTot)
!		do jType=1,nTypesTot
!			write(*,'(i5,11f7.2)')idLocalType(jType),( epskDAKelvin(jType,iType),iType=1,nTypesTot )
!		enddo

	ENDIF ! DISPLAY RELEVANT kijAD

	if(iEosOpt.eq.8)CALL GetVp(nComps,ID,iErrVp)	
	outFile=TRIM(PGLinputDir)//'\CheckSpeadmdReading.txt'
	open(619,file=outFile)
	write(619,'(a)')' iComp,jType,       ID,   nDegree,nAcceptors,nDonors,eAcceptorKcal_mol,eDonorKcal_mol	'
	do i=1,nComps
		do j=1,nTypes(i)
			write(619,601)i,j,id(i),nDegree(i,j),nAcceptors(i,j),nDonors(i,j),eAcceptorKcal_mol(i,j),eDonorKcal_mol(i,j)
		enddo
	enddo ! i=1,nComps
	close(619)
601	format(1x,2i5,i12,i8,i11,i8,2E14.4)
	errMsgPas=Trim( ErrMsg(iErrCode) )
	RETURN
861	continue
	iErrCode=1
	if(LOUDER)write(*,*)Trim( ErrMsg(iErrCode) )
	if(LOUDER)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUDER)pause
	errMsgPas=Trim( ErrMsg(iErrCode) )
	return                      
862	continue
	iErrCode=2
	if(LOUD)write(*,*)Trim( ErrMsg(2) )
	if(LOUD)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUD)pause
	errMsgPas=Trim( ErrMsg(iErrCode) )
	return                      
863	continue  ! Error reading line of ParmsTpt. Try reading parts of dumString sequentially.
	iErrCode=3
	if(LOUDER)write(*,*)Trim( ErrMsg(3) ),jComp
	if(LOUDER)write(*,*)trim( dumString )
	read(dumString,*)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
	if(LOUDER)write(*,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
	if(LOUDER)pause 'ids and zRef read ok. Checking A1.'
	read(dumString,*)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUDER)write(*,*)'A1',(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUDER)pause 'A1 read ok. Checking A2.'
	read(dumString,*)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUDER)write(*,*)'A2',(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUDER)pause 'A2 read ok. Checking V,Tmin,Mw'
	read(dumString,*)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),vMolecNm3(iComp),tKmin(iComp),rMwDum
	if(LOUDER)write(*,*)'V,Tmin,Mw',vMolecNm3(iComp),tKmin(iComp),rMw(iComp)
	if(LOUDER)pause 'VTM read ok. Checking nTypes, nFg, idType'
	read(dumString,*)idBase,idCc,idCasTmp,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),vMolecNm3(iComp),tKmin(iComp),rMwDum,nTypes(iComp),(nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
	if(LOUDER)write(*,*)'Types',nTypes(iComp),(nDegree(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
	if(LOUDER)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUDER)pause	'Check the comp#. Thats all folks!'
	errMsgPas=Trim( ErrMsg(iErrCode) )

	return                      
	END	! Subroutine GetTpt()

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG, 2010																C
!C This function replaces the SPEAD parameters with those of SPEAD11 correlation	C
!C Please note that it has been optimized for n-alkanes								C
!C Also, its applicability to mixtures has not been tested, yet.					C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE swap(nComps,nFg)
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DOUBLEPRECISION mShape(NMX)
	DIMENSION nCarbons(NMX),nFg(NMX,maxTypes)
	!common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!common/FileSwap/isFileSwap,isSw
	data initialCall/1/
	DO iComp=1,nComps
		nCarbons(iComp)=0
		DO iType=1,nTypes(iComp)
			nCarbons(iComp)=nCarbons(iComp)+nFg(iComp,iType)
		ENDDO
	ENDDO
	z1=0.d0
	DO iComp=1,nComps	
		mShape(iComp)=1.d0-0.14d0*(nCarbons(iComp)-1)+0.28268d0*(nCarbons(iComp)-1)*(nCarbons(iComp)-1)/nCarbons(iComp)
		dm1=(mShape(iComp)-1.d0)
		dm=dm1/mShape(iComp)
		dm2=dm*dm 
		dm4=dm2*dm2

		zRefCoeff(iComp,1)=mShape(iComp)*(4.d0-2.d0*dm*(1.d0+z1*dm))-3.d0
		zRefCoeff(iComp,2)=mShape(iComp)*(-2.d0+4.d0*dm*(1.d0+z1*dm))+3.d0
		zRefCoeff(iComp,3)=-2.d0*mShape(iComp)*dm*(1.d0+z1*dm)-1.d0

		a1Coeff(iComp,1)=(-5.d0+(mShape(iComp)-1.d0)*167.d0-dm2*mShape(iComp)*174.d0)*100.d0
		a1Coeff(iComp,2)=(-27.d0-(mShape(iComp)-1.d0)*184.d0+dm2*mShape(iComp)*164.d0)*100.d0
		a1Coeff(iComp,3)=(-29.d0-(mShape(iComp)-1.d0)*304.d0+dm2*mShape(iComp)*194.d0)*100.d0
		a1Coeff(iComp,4)=0.d0

		a2Coeff(iComp,1)=(-27.d0-(mShape(iComp)-1.d0)*119.d0-(mShape(iComp)-1.d0)*dm4*81.d0)*10000.d0  
		a2Coeff(iComp,2)=(195.d0+(mShape(iComp)-1.d0)*275.d0+(mShape(iComp)-1.d0)*dm4*593.d0)*10000.d0  
		a2Coeff(iComp,3)=(-407.d0-(mShape(iComp)-1.d0)*150.d0-(mShape(iComp)-1.d0)*dm4*1358.d0)*10000.d0 
		a2Coeff(iComp,4)=(239.d0-(mShape(iComp)-1.d0)*8.d0+(mShape(iComp)-1.d0)*dm4*850.d0)*10000.d0	
	ENDDO

	!added by JRE20161010 for Sai Venkata Ramana to help with PrismTpt
	isSW=0
	isFileSwap=0
	write(*,*)'Enter 1 to swap coeffs from ...input\PrismTptAcoeffs.txt or 0 to use AFG Spead11'
	read(*,*)isFileSwap
	write(*,*)'Enter 1 if this is a SW potential or 0 otherwise'
	if(isFileSwap)read(*,*)isSW

	if(initialCall .and. isSW)then
		!initialCall=0
		write(*,*)'Enter eps_kB'
		read(*,*)eps_kB
	endif
	if(isSW .or. isFileSwap)then
		open(51,file='c:/spead/calceos/input/PrismTptAcoeffs.txt')
		read(51,*)nDeck
		if(nDeck < 1)return
		do i=1,nDeck
			read(51,*)idComp,(zRefCoeff(i,j),j=1,3),(a1Coeff(i,j),j=1,3),(a2Coeff(i,j),j=1,4)
			do j=1,4
				if(isSw)a1Coeff(i,j)=a1Coeff(i,j)*eps_kB
				if(isSw)a2Coeff(i,j)=a2Coeff(i,j)*eps_kB*eps_kB
			enddo
		enddo
		close(51)
	endif !isSW
	if(initialCall)then
		initialCall=0
		iComp=1
		if(LOUD)write(*,'(a,4F10.0)')' zRefCoeff:',(zRefCoeff(1,j),j=1,4)
		if(LOUD)write(*,'(a,4F10.0)')' a1Coeff:',(a1Coeff(1,j),j=1,4)
		if(LOUD)write(*,'(a,4F10.0)')' a2Coeff:',(a2Coeff(1,j),j=1,4)
		if(LOUD)write(*,*)'  eta      A1        A2'
		do i=1,16
			eta=i*.05d0
			!copied from TptTerms()
			eta2=eta*eta
			eta3=eta*eta2
			eta4=eta*eta3
			ePart1=exp(-2*eta3)
			ePart2=1.d0/(0.2d0+eta) 
			denom=(1.d0+50.d0*eta3)
			A1=a1Coeff(iComp,1)*eta+a1Coeff(iComp,2)*eta*ePart1+a1Coeff(iComp,3)*eta3*ePart2		
			A2=(a2Coeff(iComp,1)*eta+a2Coeff(iComp,2)*eta2+a2Coeff(iComp,3)*eta3+a2Coeff(iComp,4)*eta4)/denom
			if(LOUD)write(*,'(f7.4,2F10.0)')eta,A1,A2
		enddo
	endif
	RETURN
	END	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QueryParPureSpead(iComp,iParm,value,iErr)
	USE SpeadParms      !Just for Spead
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr,i
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	! eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes)
	! bondVolNm3(nmx,maxTypes),nDegree(nmx,maxTypes)
	! idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	! nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes)
	! zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX)
	iErr=0
	do i=1,3
		if(iParm==i)value=zRefCoeff(iComp,i)
		if(iParm==i+3)value=a1Coeff(iComp,i)
		if(iParm==i+6)value=a2Coeff(iComp,i)
	enddo
	if(iParm==10)value=a2Coeff(iComp,4)
	if(iParm==11)value=vMolecNm3(iComp)
	if(iParm > 11)iErr=1
	return
end	!Subroutine QueryParPureSpead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPureSpead(iComp,iParm,value,iErr)
	USE SpeadParms      !Just for Spead
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr,i
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	! eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes)
	! bondVolNm3(nmx,maxTypes),nDegree(nmx,maxTypes)
	! idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	! nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes)
	! zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX)
	iErr=0
	do i=1,3
		if(iParm==i)zRefCoeff(iComp,i)=value
		if(iParm==i+3)a1Coeff(iComp,i)=value
		if(iParm==i+6)a2Coeff(iComp,i)=value
	enddo
	if(iParm==10)a2Coeff(iComp,4)=value
	if(iParm==11)vMolecNm3(iComp)=value
	if(iParm > 11)iErr=1
	return
end	!Subroutine QueryParPureSpead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																												  C
!C Modified by AFG 2011 : isZiter is redefined																							  C
!C																												  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	SUBROUTINE FuTpt( tKelvin,pMPa,gmol,nComps,LIQ,chemPo,rhoMol_cc,zFactor,aRes,uRes,iErr )	  !AUG 10
	!SUBROUTINE FuTpt(tKelvin,pMPa,gmol,nComps,LIQ,chemPo,zFactor,ier)
	USE GlobConst  !aRes_RT,uRes_RT etc
	USE SpeadParms !GlobConst(some)+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	LOGICAL LOUDER
	character*111 errMsg(22)
	integer ier(12)	!assocFlag,
	DIMENSION xFrac(NMX),gmol(NMX),chemPo(NMX)
	DIMENSION zRefMix(5),bRef(NMX)!,aMixQuadPart(NMX),a1Mix(5),a2Mix(5)
	!COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	!COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	!COMMON/DEPFUN/dU_RT,dA_RT,DSONK,dH_RT
	!COMMON/ETA/ETAL,ETAV,ZL,ZV
	data initCall/1/
	!COMMON/eta/etaL,etaV,zL,zV
	!COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	!       C$FUGI.f90
	!  LATEST REVISION : 2/02 jre
	!  ROUTINES CALLED:(DNEQNF(HYBRSUBS), MINPACK)               
	!  NAS IS THE NUMBER OF ACCEPTOR SITES
	!  NDS IS THE NUMBER OF DONOR SITES
	!  ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!  KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!  Z IS PV/NoKT HERE  
	!
	!Example 1. nHexane(11) at 300K and BP
	!    z1	z2	z3	a11	a12	a13	a14	a21	a22	a23	a24	vEffNm3	tabDelimited(paste to xls)
	! 1.814934	4.290586	-5.552626	-4016.6704	-6\078.570772	-16772.1984	29014.80941	-805614.5418	3021284.632	-7246589.144	5343634.813	0.098133	202.92	86.1774	2	2	101	4	201	nHexane	
	! eta=0.55: a1= -4183.4, a2= -5257.9; z1= -3638.0; z2= 14673.0; aAtt= -14.00; zAtt= -11.96; uAtt= -14.06;
	! zRep=25.033; zFactor= 14.07; aDep= -5.029; uDep= -14.06 (zAssoc=uAssoc=aAssoc=0)
	! PMPa(init)=0.023956,converge to; eta=0.45397; zFactor= 0.00125; lnPhi=-0.217860
	! Converge to PMPa(finl)= 0.0194 ;etaL=0.45397; zLiquid= 0.00101; etaVap=0.00464; zVap=0.99188
	!Example 2. Acetamide(2853) at 300K and BP
	!    z1	z2	z3	a11	a12	a13	a14	a21	a22	a23	a24	vEffNm3	tabDelimited(paste to xls)
	! 1.277100	2.791500	-3.989700	-5867.236	-14078.092	11230.159	2185.834	-1807958.7	9935591.3	-24396602.4	7055745.1	0.052351	361	59.0680	4	1	115	1	905	1	1605	1	1903	AcetAmide
	! HbMn Sb	nDs	nAs	bondVNm3	bVolSlo	eHbKcal_mol 	
	! 19 3	1	1	0.000800	0	4.5	NH2- in primary amide
	! eta=0.55: a1= -5417.2, a2= -29991.6; z1= -5338.9; z2= 19487.8; aAtt= -18.39; zAtt= -17.58; uAtt= -18.72;
	!  bigYhb=1898; ralph=11.27=sumInit; converge-> fAssoc=0.9566; XA=0.0849
	! zRep=19.6644; zFactor= -0.8387; aDep= -15.107; uDep= -25.635; (zAssoc= -3.9233; uAssoc= -6.9117; aAssoc= -4.017)
	! PMPa(init)=5.25E-5,converge to; eta=0.55622; zFactor= 1.205E-6; lnPhi=-2.494 
	! Converge to PMPa(finl)= 4.377E-6 ;etaL=0.55622; zLiquid= 0.000000; etaVap=0.000000; zVap=0.999996
	!Example 3. Ethanol+water@353K, 0.5x, Pguess=0.08392
	! eta=.55	rhoMol_cc=.0302, zAssoc=-5.67, aAssoc=-5.33							
	! Kij=0	bigYhb	eHb	ralph	   XA		zRep=	16.77	A1	-2236
	! EtOH	1912	5.3	13.96	5.86E-02	zAtt=	-6.128	A2	-2082
	! H2O	94.797	3.2	4.986	1.48E-01	zFactor= 5.974	aDep	-5.38
	! Cnvg eta=0.4588, zAssoc= -4.088, aAssoc=-4.46, zLiq=0.00113								
	!   zRep	zAtt			   XA	     fugAssoc			
	! 8.8001	-5.714			8.20E-02	-9.259210773			
	! 			             	2.00E-01	-7.830941658			
	iErr=0						 
	errMsg(5)=' FuTpt Warning: T < Tmin for at least one component.'
	errMsg(11)=' FuTpt Error: Nonsense on input.'
	errMsg(12)=' FuTpt Error: eta converged to -ve value.'
	errMsg(13)=' FuTpt Error: FuTptVtot returned with error.'
	errMsg(14)=' FuTpt Error: zFactor < 0 on final iteration.'
	errMsg(15)=' FuTpt Error: eta > 1 during iteration.'
	errMsg(16)=' FuTpt Error: z iteration did not converge. Returning crude goldenz.'
	errMsg(17)=' FuTpt Error: etaPure ~0 for at least one compound.'
	LOUDER=LOUD
	!LOUDER=.TRUE.

	ier=0 !vector init
	zRefMix=0      
	nTptCoeffs=4+1;
	totMoles=1 !we may assume that totMoles=1, since sum(xFrac)=1
    isTtooLow=0
    TcVolatile=123456
	totMoles=SUM(gmol(1:nComps))
	if(totMoles < zeroTol)then
		iErr=11
		return
	endif
	if(MINVAL(etaPure(1:nComps)) < zeroTol)then
		iErr=12
		return
	endif
	bVolMix=0
	bRefMix=0
	Vmix=0
	if(initCall.and.LOUDER)print*,'FuTpt: T(K) = ',tKelvin
    DO i=1,nComps
		xFrac(i)=gmol(i)/totMoles
		Vmix=Vmix+xFrac(i)*bVolCc_mol(i)/etaPure(i) !etaPure is evaluated at Tmin for each compd.
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        chemPo(i)=86.8686D0 !initialize to error indicator in case of bad return
        if(Tc(i) < TcVolatile)then
            iVolatile=i
            TcVolatile=Tc(i)
        endif 
    enddo
    if(tKelvin < tKmin(iVolatile))then	   ! Don't round Tmin down here or NUMDERVS may fail later.
        isTtooLow=1
        iErr=5
        if(LOUDER.and.initCall)write(*,'(a,i3,2f7.1,a)')' FuTpt: iVolatile,T(K),TKmin= ',iVolatile,tKelvin,tKmin(iVolatile),TRIM( errMsg(iErr) )
    endif
    if(isTtooLow==1)then
        iErr=5
        if(LOUDER)write(*,*)TRIM(errMsg(iErr))
        !goto 86		! declare warning if T < Tmin, but compute properties anyway. e.g. LDEN probably OK although Psat might be nonsense.
    endif
    hRes_RT=86.8686D0 !initialize to error indicator in case of bad return
    uRes_RT=86.8686D0 !initialize to error indicator in case of bad return
    Ares_RT=86.8686D0 !initialize to error indicator in case of bad return
    Sres_R =86.8686D0 !initialize to error indicator in case of bad return
    zFactor=86.8686D0 !initialize to error indicator in case of bad return
    
	
	!  INITIATE SECANT ITERATION ON eta

	IF (rGas.Le.0 .or. tKelvin.le.0)then
		iErr=11
		if(LOUDER)write(*,*)' FuTpt: Nonsense. rGas=',rGas,' tKelvin=',tKelvin
		!call BeepMsg(errMsg(iErr))
		goto 86
	endif
	Pb_RT=pMPa*bVolMix/RGAS/tKelvin
	P_RT=pMPa/RGAS/tKelvin
	pb_rt=P_RT*bVolMix

	!  GUESS FOR eta
	eta=pb_rt/1.5d0 ! keep initial guess on the low side to avoid metastable region. 
	etaHi=1.2D0*bVolMix/Vmix+0.24d0*( exp(-PMPa/100) ) !increase the initial guess for high pressure. etaMax might be less than 1 for some compounds.
	if(etaHi < 0.6d0)etaHi=0.6d0
	if(etaHi > etaMax)etaHi=etaMax/1.1d0
	IF(LIQ==1.or.LIQ==3.or.eta > etaMax)eta=etaHi
	etaRef=eta*bRefMix/bVolMix
	IF (eta > etaMax .or. etaRef > etaMax .or. eta < 0)then
		iErr=11
		if(LOUDER)write(*,*)' pMPa=',pMPa,' eta=',eta,' etaRef=',etaRef
		!call BeepMsg(errMsg(iErr))
		goto 86
	endif
	rhoMol_cc=eta/bVolMix
	vTotCc=totMoles/rhoMol_cc
	isZiter=1
	if(initCall.and.LOUDER)write(*,'(a,6E12.4)')' FuTpt: first call to FuVtot. T(K),eta=',tKelvin,eta
	call FuTptVtot(isZiter,zFactor,aDep,uDep,chemPo,vTotCc,tKelvin,gmol,nComps,iErrZ)
	if(iErrZ > 10)then
		iErr=13
		if(LOUDER)write(*,*)'FuTpt: FuTptVtot returned with error code:',iErrZ
		!call BeepMsg(errMsg(iErr))
		goto 86
	endif
	etaOld=eta
	errOld=pb_rt-eta*zFactor
	eta=etaOld/1.1D0 !set this "new" value to ideal gas value so if change < 1E-9 on first iteration, we get exactly the ideal gas result.
    !if(LIQ==1.or.LIQ==3)eta=0.8 !Take a big step for liquids to bracket answer sooner. Saves about 4/13 iterations
	!ETA=rhoMol_Cc*bVolMix	!calculate here before entering loop, then at the end with each new rhoMol_Cc
	change=1234
	itMax=55
	NITER=0
	pMax=-1.d11 !for crude goldenz
	pMin=1.d11
	isZiter=1
    iterGold=0
	if(LOUDER.and.initCall)write(*,'(a,i4,6E12.4)')' FuTpt init: liq,eta', liq,eta
100 continue
	do while(ABS(change) > 1.e-9 .and. iErr < 11)	! using change as criterion, instead of change/eta, means ~9 sig figs on liquid density, while ideal gas is exactly returned for vapor density, even when P-> 1E-11.
		NITER=NITER+1
		if(niter > itMax)iErr=16
		if(iErr==16)then
			!call BeepMsg(errMsg(iErr))
			cycle
		endif
		!ETA=rhoMol_Cc*bVolMix
		if( (1-eta) < 1.d-11)then
			if(LOUDER)write(*,*) ' FuTpt: eta > 1. nIter=',nIter
			do iComp=1,nComps
				if(LOUDER)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
				if(LOUDER)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			enddo
			if(LOUDER)pause
		endif
		rhoMol_cc=eta/bVolMix
		vTotCc=totMoles/rhoMol_cc
		call FuTptVtot(isZiter,zFactor,aDep,uDep,chemPo,vTotCc,tKelvin,gmol,nComps,iErrZ)
		if(iErrZ > 10)then
			iErr=13
			if(LOUDER)write(*,*)'FuTpt: FuTptVtot returned with error code:',iErrZ
			!call BeepMsg(errMsg(iErr))
			cycle
		endif
		error=pb_rt-eta*zFactor
		if(LOUDER)write(*,'(a,i3,F9.5,2(E12.4))')' FuTpt: After FuVtot. iter,eta,Z,error=',NITER,eta,zFactor,error
		if(eta*zFactor > pMax .and. eta < 0.15)then
			pMax=eta*zFactor
			etaAtPmax=eta	 !for crude goldenz. ~randomly searches for lo eta max in p
		endif
		if(eta*zFactor < pMin .and. eta > 0.15)then
			pMin=eta*zFactor
			etaAtPmin=eta	 !for crude goldenz. ~randomly searches for hi eta min in p
		endif
		change=error/(error-errOld)*(eta-etaOld)
		!print*,'FuTpt: eta,change=',eta,change
		etaOld=eta
		errOld=error
        if(zFactor < 0)change=1 !force another iteration if Z < 0. This happens when P ~ 1E-12.
		if(ABS(change/eta) > 0.3)change=eta*DSIGN(0.2d0,change)	!step limiting. 0.3>0.2 to give max looseness
		!if(ABS(change/zFactor) > 1)change=zFactor*DSIGN(0.5d0,change)	!step limiting. Stop zFactor from going negative.
		eta=eta-change
		rhoMol_cc=eta/bVolMix
		if(eta < 0 .or. eta > etaHi)then !NOTE: this should not happen given step limiting
			if(niter < (itMax-1))niter=itMax-1 !restrict tries with new guess so we can crash if necessary.
			if(LOUDER)write(*,'(a,f8.4)')' Warning in FuTpt: next guess is eta=',eta
			eta=0
			if(liq==1.or.liq==3)eta=etaHi
			if(LOUDER)write(*,*) 'Restarting iteration with eta=',eta
			if(LOUDER)pause
		endif
	enddo  !iteration on eta to find P=P(input)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Whoohoo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(iErr > 10)then
		if(LOUDER)print*,'FuTpt: Convergence failed, attempting crude estimate.'
        iterGold=iterGold+1
		if(liq==0 .or. LIQ==2)eta=etaAtPmax  !crude goldenz for vapor
		if(liq==1.or.LIQ==3)then
			eta=etaAtPmin*0.99  !crude goldenz for liquid
			if(Pmin < zeroTol)then
				if(Pmax > zeroTol)then
					eta=Pmin+(Pmax-0)*(etaAtPmax-etaAtPmin)/(Pmax-Pmin)+0.1d0
					if(eta > etaMax)eta=etaMax
				else
					iErr=21
					goto 86
				endif
			endif

			rhoMol_cc=eta/bVolMix
			vTotCc=totMoles/rhoMol_cc
			call FuTptVtot(isZiter,zFactor,aDep,uDep,chemPo,vTotCc,tKelvin,gmol,nComps,iErrZ)
		    if(eta < 0 .or. eta > etaMax)eta=etaHi ! liq==1 or 3 in this if.
			if(LOUDER)write(*,'(a,f8.4,1x,a,i3,1x,a,f7.4,a,1PE11.4,a,1PE11.4)')' FuTpt: P=',pMPa,' LIQ=',LIQ,' etaMin=',etaAtPmin,' pMin=',pMin,' pMax=',pMax
            if(iterGold < 2)goto 100
		endif
	endif
	IF (eta < 0)then
		iErr=12
        if(LOUDER)write(*,*)'FuTpt: eta < 0 after converged.'
		!call BeepMsg(errMsg(iErr))
	endif
	etaPass=eta  ! etaPass in GlobConst
	IF (LIQ==1 .or. LIQ==3) THEN
	  ETAL=ETA
	  ZL=zFactor  
	ELSE
	  ETAV=ETA
	  ZV=zFactor  
	ENDIF
	dU_RT=uDep
	dH_RT=uDep+zFactor-1
	dA_RT=aDep
	if(LOUDER)write(*,'(a,6E12.4)')' FuTpt: getting chemPo(T,P). T,eta,Z=',tKelvin,eta,zFactor
	initCall=0
    if(LIQ > 1)return
	!  ITERATION ON RHO HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.
	!call one more time with isZiter=0, also to facilitate debugging
	isZiter=0 ! No more iterations on Z!!!!!!!!!!!!!!!!!!
	rhoMol_cc=eta/bVolMix
	vTotCc=totMoles/rhoMol_cc
	call FuTptVtot(isZiter,zFactor,aRes,uRes,chemPo,vTotCc,tKelvin,gmol,nComps,iErrZ)
	if(iErrZ > 9)iErr=iErr+iErrZ*10
	IF (zFactor < zeroTol)then
		iErr=14
		if(LOUDER)write(*,'(a,8E12.4)')' FuTpt: eta,Z=',eta,zFactor
        if(LOUDER)pause 'FuTpt: zFactor < 0 after converged.'
        goto 86
		!IF(LOUDER)call BeepMsg(errMsg(iErr))
	endif
	if(iErr > 9)goto 86
    Ares_RT=aRes
    uRes_RT=uREs
    hRes_RT=uRes+zFactor-1
    Sres_R =uRes-aRes
	chemPo(1:nComps)=chemPo(1:nComps)-DLOG(zFactor) ! transform from TV to TP
	if(LOUDER)write(*,'(a,6E12.4)')' FuTpt: FuVtot=>chemPo=',(chemPo(i),i=1,nComps)
	if(LOUDER)write(*,'(a,6E12.4)')' FuTpt: ChemPoCalc done. aRes,uRes=',aRes,uRes
	if(zFactor > zeroTol)ChemPoPure=aRes+zFactor-1-DLOG(zFactor)
	if(nComps==1 .and. LOUDER)write(*,'(a,6E12.4)')' FuTpt: fugc(1),Gdep1',chemPo(1),ChemPoPure
	initCall=0
	RETURN

86	if(LOUDER)WRITE(6,'(a,i4)')' FuTpt: error= ',iErr
	IER(1)=iErr
	RETURN
867 continue
	if(LOUDER)write(*,*)'Exiting FuTpt: Please correct xFrac'
	iErr=17
	ier(1)=17
	ier(7)=1
	return
	END	 ! Subroutine FuTpt()

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CODED BY AV 06/26/06. cf. DannerGess/SPEADCI paper. IECR08
! PURPOSE: DEVELOPING SpeadGamma MODEL. cf.SPEADCI paper IECR, 47:7955 (08).
	SUBROUTINE FuSpeadGamma(tKelvin,pMPa,xFrac,nComps,LIQ,chemPo,zFactor,ier)
	USE SpeadParms
	USE VpDb
	USE Assoc
	IMPLICIT DoublePrecision(A-H,O-Z)
	DIMENSION xFrac(NMX),chemPo(NMX),fugcPure(NMX),fugcL(nmx),rLnGam(nmx),x(nmx),fugc(nmx),ier(12)
	iErr=0

	do index=1,12
		ier(index)=0
	enddo      
	!iEosOpt=5
   
		xInf=1e-11
	do iPure=1,nComps
		sum=0
		do iComp=1,nComps
			x(iComp)=xInf
			sum=sum+x(iComp)
		enddo
		x(iPure)=1-sum
		nc=2
 		call FuTpt( tKelvin,pMPa,xFrac,nComps,1,fugcL,rhoMol_cc,zL,aRes,uRes,iErr )	  
		!call FuTpt(tKelvin,pMPa,x,nc,1,fugcL,zL,ier)
		fugcPure(iPure)=fugcL(iPure)
		!hPure(iPure)=DHONKT
	enddo
	!call FuTpt(tKelvin,pMPa,xFrac,nComps,1,fugcL,zL,ier)
 	call FuTpt( tKelvin,pMPa,xFrac,nComps,1,fugcL,rhoMol_cc,zL,aRes,uRes,iErr )	  
	!CALL GetVp(nComps,ID,iErrVp)	


	gIs=0
	gE=0
	gMix=0
	!hIs=0
	do iComp=1,nComps
		gIs=gIs+x(iComp)*LOG(x(iComp))
		!hIs=hIs+x(iComp)*hPure(iComp)
		rLnGam(iComp)=fugcL(iComp)-fugcPure(iComp) 
		gE=gE+xFrac(iComp)*rLnGam(iComp)
		gMix=gIs+gE
	enddo
	!hE=DHONKT-hIs

	zFactor=1
	fugc(iComp)=1
	etaV=0
	if(LIQ.eq.1)then
		etaL=0.5
		zFactor=zL !why not use the spead Z?
		do kComp=1,nComps
			pSatExp=exp(vpCoeffs(kComp,1)+vpCoeffs(kComp,2)/tKelvin+vpCoeffs(kComp,3)*DLOG(tKelvin)+vpCoeffs(kComp,4)*tKelvin**vpCoeffs(kComp,5))/1000000 !todo: fill in the vp eqn...
			if(pSat.lt.1e-33)pSat=1E-33 !logArg in fugc() cannot be zero
			gamma=exp(rLnGam(kComp))
			!if(LOUD)write(*,'(i3,1x,2(f7.4),f10.3,f10.3,e11.4)')kComp,xFrac(kComp),volFrac(kComp),fugc(kComp),fugAssoc(kComp),gamma
			chemPo(kComp)=LOG(pSatExp/pMPa) + rLnGam(kComp)	!log of product is sum of logs 
		enddo
	endif
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine TptTerms(isZiter,nComps,iComp,eta,etaRef,a0i,a1i,a2i,z0i,z1i,z2i,iErr)
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	character*77 errMsg(11)
	!common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!DIMENSION dEta_dN(NMX),dEtaRef_dN(NMX)
	!COMMON/dervZ/dz0i_dN(NMX),dz1i_dN(NMX),dz2i_dN(NMX)
	!COMMON/derva/da1i_dN(NMX),da2i_dN(NMX),da0i_dN(NMX)
	!COMMON/dervZ2/dz0i_deta,dz1i_deta,dz2i_deta,d2z0i_deta2,d2z1i_deta2,d2z2i_deta2
	!COMMON/dervZ3/d3z0i_deta3,d3z1i_deta3,d3z2i_deta3
	!COMMON/derva2/da1i_deta,da2i_deta,da0i_deta,d2a1i_deta2,d3a1i_deta3,d2a2i_deta2

	iErr=0
	NC=nComps
	errMsg(1)='TptTerms Error: etaRef out of range.'
	errMsg(2)='TptTerms Error: eta out of range.'
	if(etaRef > etaMax .or. etaRef < 0)then
		iErr=1
		if(LOUD)write(*,*)'In TptTerms: eta=',eta
		return
	endif
	if(eta > etaMax .or. eta < 0)then
		iErr=2
		if(LOUD)write(*,*)'In TptTerms: eta=',eta
		return
	endif
	c1=zRefCoeff(iComp,1)+3
	c2=zRefCoeff(iComp,2)-3
	c3=zRefCoeff(iComp,3)+1
	void=1.d0-etaRef
	void2=void*void
	void3=void*void2
	if(void.le.0 .and. LOUD)write(*,*)'TptTerms: LOG error for a0i. iComp,etaRef,void=',iComp,etaRef,void
	if(void.le.0 .and. LOUD)write(*,*)'TptTerms: LOG error for a0i. iComp,etaRef,void=',iComp,etaRef,void
	a0i=0.5d0*(c1+c2+c3)/void2-(c2+2*c3)/void-c3*LOG(void)-0.5d0*(c1-c2-3*c3)
	dA0iDeta=(c1+etaRef*(c2+etaRef*c3))/void3  ! =zRef/eta	 
	etaDA0iDeta=etaRef*dA0iDeta  ! =zRef	 

	eta2=eta*eta
	eta3=eta2*eta
	eta4=eta2*eta2
	eta5=eta4*eta	
	const=2.d0
	ePart1=exp(-const*eta3)
	dePart1=-const*3.d0*eta2*ePart1
	ePart2=1.d0/(0.2d0+eta) 
	ePart2Sq=ePart2*ePart2
	ePart2Qu=ePart2Sq*ePart2Sq
	dePart2=-1.d0*ePart2Sq
	denom=(1.d0+50.d0*eta3)
	a1iCrude=a1Coeff(iComp,1)*eta+a1Coeff(iComp,2)*eta*ePart1+a1Coeff(iComp,3)*eta3*ePart2+a1Coeff(iComp,4)		!a1Coeff(iComp,2)*eta2+
	da1i_detaCrude=a1Coeff(iComp,1)+a1Coeff(iComp,2)*(ePart1+eta*dePart1)+a1Coeff(iComp,3)*(3.d0*eta2*ePart2+eta3*dePart2)	  !+2.d0*a1Coeff(iComp,2)*eta+

	a2iCrude=(a2Coeff(iComp,1)*eta+a2Coeff(iComp,2)*eta2+a2Coeff(iComp,3)*eta3+a2Coeff(iComp,4)*eta4)/denom
	da2i_detaCrude=(a2Coeff(iComp,1)+2.d0*a2Coeff(iComp,2)*eta+3.d0*a2Coeff(iComp,3)*eta2+4.d0*a2Coeff(iComp,4)*eta3-150.d0*eta2*a2iCrude)/denom

	bipKij=KIJ(iComp,iComp)+KTIJ(iComp,iComp)*eta
	dbipKij_deta=KTIJ(iComp,iComp)
	a1i=a1iCrude*(1.d0-bipKij)	
	a2i=a2iCrude*(1.d0-bipKij)
	da1i_deta=da1i_detaCrude*(1.d0-bipKij)+a1iCrude*dbipKij_deta
	da2i_deta=da2i_detaCrude*(1.d0-bipKij)+a2iCrude*dbipKij_deta

	z0i=etaDA0iDeta
	z1i=eta*da1i_deta
	z2i=eta*da2i_deta

	if(a1i > 0 .or. a2i > 0)then
		if(LOUD)write(*,'(a,i4,E11.4,3e12.4)')' TptTerms: a1 or a2 > 0.  iComp, eta, a1i, a2i=',iComp,eta,a1i,a2i
		!if(LOUD)pause 'Fatal error.  Sorry.  Check your ParmsTpt.'
		!a2i = -1D-11
        iErr=2
        return
	endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2010																							  C
!C Derivatives of ref. and perturbation terms in respect to eta and N											  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	if (isZiter==0) then
		dz0i_deta=(c1+2.d0*c2*etaRef+3.d0*c3*etaRef*etaRef)/void3+3.d0*z0i/void
		d2z0i_deta2=(2.d0*c2+6.d0*c3*etaRef)/void3+6.d0*dz0i_deta/void-6.d0*z0i/void2
		d3z0i_deta3=6.d0*(c3+z0i)/void3-18.d0*dz0i_deta/void2+9.d0*d2z0i_deta2/void
		d2ePart1=-const*6.d0*eta*ePart1-const*3.d0*eta2*dePart1
		d3ePart1=-const*6.d0*ePart1-const*12.d0*eta*dePart1-const*3.d0*eta2*d2ePart1
		d4ePart1=-18.d0*const*dePart1-18.d0*const*eta*d2ePart1-const*3.d0*eta2*d3ePart1
		d2ePart2=2.d0*ePart2Sq*ePart2
		d3ePart2=-6.d0*ePart2Qu
		d4ePart2=24.d0*ePart2Qu*ePart2
		
		d2a1i_deta2Crude=a1Coeff(iComp,2)*(2.d0*dePart1+eta*d2ePart1)+a1Coeff(iComp,3)*(6.d0*eta*ePart2+6.d0*eta2*dePart2+eta3*d2ePart2)
		d3a1i_deta3Crude=a1Coeff(iComp,2)*(3.d0*d2ePart1+eta*d3ePart1)+a1Coeff(iComp,3)*(6.d0*ePart2+18.d0*eta*dePart2+9.d0*eta2*d2ePart2+eta3*d3ePart2)
		d4a1i_deta4Crude=a1Coeff(iComp,2)*(4.d0*d3ePart1+eta*d4ePart1)+a1Coeff(iComp,3)*(24.d0*dePart2+36.d0*eta*d2ePart2+12.d0*eta2*d3ePart2+eta3*d4ePart2)
		d2a1i_deta2=d2a1i_deta2Crude*(1.d0-bipKij)+da1i_detaCrude*dbipKij_deta
		d3a1i_deta3=d3a1i_deta3Crude*(1.d0-bipKij)+d2a1i_deta2Crude*dbipKij_deta
		d4a1i_deta4=d4a1i_deta4Crude*(1.d0-bipKij)+d3a1i_deta3Crude*dbipKij_deta

		dz1i_deta=da1i_deta+d2a1i_deta2*eta
		d2z1i_deta2=d2a1i_deta2+d2a1i_deta2+d3a1i_deta3*eta
		d3z1i_deta3=d3a1i_deta3+d3a1i_deta3+d3a1i_deta3+d4a1i_deta4*eta
		
		d2a2i_deta2Crude=(2.d0*a2Coeff(iComp,2)+6.d0*a2Coeff(iComp,3)*eta+12.d0*a2Coeff(iComp,4)*eta2-300.d0*eta*a2iCrude-300.d0*eta2*da2i_detaCrude)/denom
		d3a2i_deta3Crude=(6.d0*a2Coeff(iComp,3)+24.d0*a2Coeff(iComp,4)*eta-900.d0*eta*da2i_detaCrude-450.d0*eta2*d2a2i_deta2Crude-300.d0*a2iCrude)/denom
		d4a2i_deta4Crude=(24.d0*a2Coeff(iComp,4)-1200.d0*da2i_detaCrude-600.d0*eta2*d3a2i_deta3Crude-1800.d0*eta*d2a2i_deta2Crude)/denom
		d2a2i_deta2=d2a2i_deta2Crude*(1.d0-bipKij)+da2i_detaCrude*dbipKij_deta
		d3a2i_deta3=d3a2i_deta3Crude*(1.d0-bipKij)+d2a2i_deta2Crude*dbipKij_deta
		d4a2i_deta4=d4a2i_deta4Crude*(1.d0-bipKij)+d3a2i_deta3Crude*dbipKij_deta
		
		dz2i_deta=da2i_deta+d2a2i_deta2*eta
		d2z2i_deta2=d2a2i_deta2+d2a2i_deta2+d3a2i_deta3*eta
		d3z2i_deta3=d3a2i_deta3+d3a2i_deta3+d3a2i_deta3+d4a2i_deta4*eta

	endif

	return
	end
	
 !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 ! Compute the derivatives of the third order GaussEx taking Gex parms, eta and T, as input.
 ! phi=sum(iComp)sum(iType) [bVol0(iComp,iType)*sum(jComp,jType)epsij*ratioA3ij]	where ratioA3=A3/A2.
 ! Since 3*A3/A2 = (g0)*exp(-g1*[(eta-eta0)^2-eta0^2]), then phi is already a Gaussian. For Ahmad, gFun=exp(phi/tKelvin)
	Subroutine GaussFun3(tKelvin,eta,Gf,dGf_dEta,d2Gf_dEta2,cpGf,uDepGf,iErr)
	USE SpeadParms, only:GexG0,GexG1,GexEta0
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!common/GaussFun/G0,G1,eta0
	G0=GexG0(1)
	G1=GexG1(1)
	eta0=GexEta0(1)
	iErr=0
	phi=G0*EXP(  -G1*( (eta-eta0)*(eta-eta0)-eta0*eta0 )  )
	dPhi_dEta=phi*( -2*G1*(eta-eta0) )
	d2Phi_dEta2=dPhi_dEta*( -2*G1*(eta-eta0) )-2*G1*phi
	gFun = phi/tKelvin
	dgFun_dEta = dphi_dEta/tKelvin
	d2gFun_dEta2 = d2phi_dEta2/tKelvin
	gFun2 = gFun*gFun
	gFun3 = gFun2*gFun
	gFun4 = gFun3*gFun
	gFun8 = gFun4*gFun4
	dgFun2_deta = 2.d0*gFun*dgFun_deta
	d2gFun2_deta2 = 2.d0*dgFun_deta*dgFun_deta + 2.d0*gFun*d2gFun_deta2
	dgFun4_deta = 4.d0*gFun3*dgFun_deta
	d2gFun4_deta2 = 4.d0*3.d0*gFun2*dgFun_deta*dgFun_deta + 4.d0*gFun3*d2gFun_deta2
  	expg = exp(gFun)
	dexpg_deta = dgFun_deta*expg
	d2expg_deta2 = d2gFun_deta2*expg+dgFun_deta*dexpg_deta

!	cpFun
	cpPhi=eta*dPhi_dEta	 !JRE: I think this is true, but I'm not 100% sure of Ahmad's definition.
	cpgFun = cpPhi/tKelvin
	cpgFun2 = 2.d0*gFun*cpgFun
	cpdgFun4 = 4.d0*gFun3*cpgFun
	cpexpg = cpgFun*expg
	if (gFun < 1.d-5) then
		Gf = 1.d0 + gFun/3
		dGf_dEta = dGfun_dEta/3
		d2Gf_dEta2 = d2Gfun_dEta2/3
		uDepGf=gFun/3
	else
		!Gf = 1+2.d0/gFun2*(expg - 1 - gFun - gFun2/2)
		Gf = 1.d0 + 2.d0*expg/gFun2 - 2.d0/gFun2 - 2.d0/gFun - 1.d0
		dGf_dEta = ((2.d0*dexpg_deta*gFun2 - dgFun2_deta*2.d0*expg)/gFun4 + 2.d0*dgFun2_deta/gFun4 + 2.d0*dgFun_deta/gFun2)
		term1 = ( (2.d0*d2expg_deta2*gFun2 - 2*expg*d2gFun2_deta2)*gFun4 - dgFun4_deta*(2.d0*dexpg_deta*gFun2 - dgFun2_deta*2.d0*expg) )/gFun8
		d2Gf_dEta2 = (term1 + (2.d0*d2gFun2_deta2*gFun4-dgFun4_deta*2.d0*dgFun2_deta)/gFun8 + (2.d0*d2gFun_deta2*gFun2 - dgFun2_deta*2.d0*dgFun_deta)/gFun4)
		cpGf = ((2.d0*cpexpg*gFun2 - cpgFun2*2.d0*expg)/gFun4 + 2.d0*cpgFun2/gFun4 + 2.d0*cpgFun/gFun2)
		!uDepGf	= beta*dGf_dBeta=4/gfun3*(beta*dgFun_dBeta)*(expg-1-gFun-gFun2/2)+2/gFun2*(beta*dgFun_dBeta)*(expg-1-gfun)
		!beta*dgFun_dBeta = gFun => uDepGf=4/gfun2*(expg-1-gFun-gFun2/2)+2/gFun*(expg-1-gFun)
		uDepGf=2/gfun*( 2*(expg-1-gFun-gFun2/2)/gFun+(expg-1-gFun) )
		!=> expg -2*expg/gFun +2/gFun +1 = -2*(expg-1-G-G^2/2)/G + (expg-1-G) 
		!=> expg -2*expg/gFun +2/gFun +1 = -2*expg/G+2/G+2+G + expg-1-G 
		!=> expg  +2/gFun +1 =  2/G+2+G + expg-1-G 	= OK!
		!uDepGf=2/gfun*( -2*(expg-1-G-G^2/2)/G + (expg-1-G) )
		!uDepGf=2/gfun*( -2*(G^3/6)/G + (G^2/2) )
		!uDepGf=2/gfun*( -2*(G^2/6) + (G^2/2) )
		!uDepGf=2*gFun*( -2/6 + 3/6 )
	endif
	return
	end

	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine MixRule(isZiter,gmol,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErr) 
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*77 errMsg(22)
	!doublePrecision KIJ,KTIJ
	!dimension cRef(NMX,3),cRefMix(3)
	DIMENSION xFrac(NMX),aMixQuadPart(NMX),gmol(NMX)
	DIMENSION solParEntro(NMX),solParEnerg(NMX) !,zRefMix(5), !a1Mix(5),a2Mix(5),
	!DoublePrecision ks0ij(NMX,NMX),ks1ij(NMX,NMX) ! USE BIPs
	DIMENSION bRef(NMX)
	!DIMENSION dSPEntro_deta(NMX),dSPEnerg_deta(NMX),dEta_dN(NMX),dEtaRef_dN(NMX),dbVolMix_dN(NMX),dbRefMix_dN(NMX)

	!COMMON/VCoeff/B2,B3

	!COMMON/z_aComps/zRefij,aRefij,dKijDeta,z1ij,a1ij,z2ij,a2ij
	!COMMON/dervZ/dz0i_dN(NMX),dz1i_dN(NMX),dz2i_dN(NMX)
	!COMMON/dervZ2/dz0i_deta,dz1i_deta,dz2i_deta,d2z0i_deta2,d2z1i_deta2,d2z2i_deta2
	!COMMON/dervZ3/d3z0i_deta3,d3z1i_deta3,d3z2i_deta3
	!COMMON/derva/da1i_dN(NMX),da2i_dN(NMX),da0i_dN(NMX)
	!COMMON/derva2/da1i_deta,da2i_deta,da0i_deta,d2a1i_deta2,d3a1i_deta3,d2a2i_deta2
	!COMMON/DervMix/dz0Mix_dN(NMX),da0Mix_dN(NMX),dz1Mix_dN(NMX),da1Mix_dN(NMX),dz2Mix_dN(NMX),da2Mix_dN(NMX),daMixQuadPart_dN(NMX,NMX),daMixQuadPart_deta(NMX)
	!COMMON/DervMix2/dz0Mix_deta,da0Mix_deta,dz1Mix_deta,da1Mix_deta,dz2Mix_deta,da2Mix_deta,d2z0Mix_deta2,d2z1Mix_deta2,d2z2Mix_deta2
	!COMMON/DervMix3/d2a1Mix_deta2,d3a1Mix_deta3,d2a2Mix_deta2
	!COMMON/DervZMix3/d3z0Mix_deta3,d3z1Mix_deta3,d3z2Mix_deta3
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!common/XsProps/a0Xs
	!data initial/1/
	! AGF. It should be 0 not -0.044  
	data init/1/		! best I could do based on avg of the entire system AV 11/29/08
							! Note that polymer solutions are very sensitive to the BIPs and these values should be optimized via KD option
	iErr=0
	errMsg(11)=' MixRule Error: a1i or a2i > 0'
	errMsg(12)=' MixRule Error: Wertheim returned error.'
	errMsg(13)=' MixRule Error: sqArg = a1i*a1j or a2i*a2j < 0'
	errMsg(14)=' MixRule Error: rMwMix .le. 0.  Check input.'
	errMsg(15)=' MixRule Error: eta out of range.'
	errMsg(5)=' MixRule Warn: T < TKmin.'
	errMsg(17)=' MixRule Error: Error returned from TptTerms().'

	totMoles=SUM(gmol)
	bVolMix=0
	bRefMix=0
    TcVolatile=123456
    DO i=1,nComps
        xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        !chemPo(iComp)=86 !initialize to error indicator in case of bad return !No chemPo here !!!
        if(Tc(i) < TcVolatile)then
            iVolatile=i
            TcVolatile=Tc(i)
        endif 
   enddo
   if(tKelvin < tKmin(iVolatile)*0.998D0)then	   !round TKmin down a little for possible lost precision
        iErr=5
        if(LOUD.and.init)write(*,'(a,i3,2f7.1,a)')' FuVtot: iVolatile,T(K),TKmin= ',iVolatile,tKelvin,tKmin(iVolatile),TRIM( errMsg(iErr) )
    endif
	if(bVolMix==0 .and. LOUD)pause 'MixRule error: bVolMix=0'
    etaRef=eta*bRefMix/bVolMix
    rhoMol_cc=eta/bVolMix
	!if(LOUD)write(*,*)'init,ks0=',init,ks0
		ks0= -0.04*18*nComps/ SUM(bRef) !ks0= -0.04*18*Nc/sum(bRef) is a crude rule so that kijS0->0 for polymer blends
		ks1=0 ! see data statement
		!if(LOUD)write(*,*)'Enter estimate for ks0'
		!read(*,*)ks0
		!if(LOUD)write(*,*)'kSo=',ks0
		do iComp=1,nComps
			do jComp=1,nComps
				ks0ij(iComp,jComp)=0
				ks1ij(iComp,jComp)=0
				if(iComp.ne.jComp)ks0ij(iComp,jComp)=ks0  !only use ks0 when i.ne.j !!!
				if(iComp.ne.jComp)ks1ij(iComp,jComp)=ks1
			enddo
		enddo
		etaStd=0.4
        etaRefStd=etaStd*bRefMix/bVolMix
		!print*,'MixRule: Calling TptTerms for SolPar calcs.'
		do iComp=1,nComps
			call TptTerms(isZiter,nComps,iComp,etaStd,etaRefStd,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
			solParEntro(iComp)=SQRT( a0i*1.987d0*298.d0*etaRefStd/bRef(iComp) )
			solParEnerg(iComp)=SQRT(-a1i*1.987d0*etaStd/bVolCC_mol(iComp) )
			!if(LOUD)write(*,'(a,i4,2(f8.2))')' i,delSi,delUi',iComp,solParEntro(iComp),solParEnerg(iComp)
		enddo

	if(eta > etaMax .or. eta < 0)then
		iErr=15
		if(LOUD)write(*,*)'MixRule: eta=',eta
		!if(LOUD)call BeepMsg( TRIM(errMsg(iErr)) )
		return
	endif
	if(etaRef > etaMax .or. etaRef < 0)then
		iErr=15
		if(LOUD)write(*,*)'MixRule: etaRef=',etaRef
		!if(LOUD)call BeepMsg( TRIM(errMsg(iErr)) )
		return
	endif

	if(nTptCoeffs.ne.4)nTptCoeffs=4;

	a0Mix=0
	z0Mix=0
	a1Mix=0
	a2Mix=0
	z1Mix=0
	z2Mix=0
	a0KijMix=0 ! new for kij=kij(eta), maybe not necessary.
	a1KijMix=0 ! new for kij=kij(eta)
	a2KijMix=0 ! new for kij=kij(eta)
	do iComp=1,nComps
		bVoli=bVolCC_mol(iComp)
		aMixQuadPart(iComp)=0.d0
		!print*,'MixRule: Calling TptTerms for Compo ai zi calcs. iComp=',iComp
		call TptTerms(isZiter,nComps,iComp,eta,etaRef,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
		if(iErrTpt.ne.0)then
			if(LOUD)print*, 'MixRule: iErrTpt=',iErrTpt
			if(LOUD)pause 'MixRule: error from TptTerms'
			iErr=17
			return
		endif
		if(a0i < 0 .or. a1i > 0 .or. a2i > 0)then
			if(LOUD)write(*,*)'MixRule: error in TptTerms:a0i<0 or a1i or a2i > 0. a0i,a1i,a2i'
			if(LOUD)write(*,*)'   eta         a0i        a1i         a2i'
			if(LOUD)write(*,'(1x,E11.4,3e13.5)')eta,a0i,a1i,a2i
			if(LOUD)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)pause 'Check your input.'
			iErr=11
			exit !terminate the do loop, Note: "cycle" would go to enddo and keep looping.
		endif
		if (nComps==1) then
            bipKii=KIJ(1,1)+KTIJ(1,1)*eta
            !bipKii=KIJ(1,1)-KTIJ(1,1)*TC(1)/tKelvin
			a0Mix=a0i !*(1-bipKii)
			a1Mix=a1i*(1-bipKii)
			a2Mix=a2i*(1-bipKii)
			z0Mix=z0i !*(1-bipKii)
			z1Mix=z1i*(1-bipKii)
			z2Mix=z2i*(1-bipKii)
			aMixQuadPart=a0i+(a1i+a2i/tKelvin)/tKelvin
			goto 100 !skip the mixture calculation and go to derivative calcs.
		!Analytical second and third virial coefficients
			!B2=((zRefCoeff(1,1)+3))+(a1Coeff(1,1)+a1Coeff(1,2))/tKelvin+(a2Coeff(1,1))/(tKelvin*tKelvin) 
			!B3=(3.d0*(zRefCoeff(1,1)+3)+(zRefCoeff(1,2)-3)+(2.d0*a2Coeff(1,2))/(tKelvin*tKelvin)) 	!+(2.d0*a1Coeff(1,2))/tKelvin
			!if(LOUD)write(*,*)B2,' ',B3 
		endif
		if(iErr > 10)then	   ! allow to continue for iErr=6(T< Tmin)
			if(LOUD)pause
			return
		endif
		chiConst=1/1.987d0/298.d0
		do jComp=1,nComps
			bVolj=bVolCC_mol(jComp)
			!print*,'MixRule: Calling TptTerms for jComp calcs. jComp=',jComp
			call TptTerms(isZiter,nComps,jComp,eta,etaRef,a0j,a1j,a2j,z0j,z1j,z2j,iErrTpt)
		   	vij=(bRef(iComp)+bRef(jComp))/2
			chiS=vij*(solParEntro(iComp)-solParEntro(jComp))**2/(1.987*298)
			!ks0= -0.04*18*nComps/ SUM(bRef) !ks0= -0.04*18*Nc/sum(bRef) is a crude rule so that kijS0->0 for polymer blends
			ksij=( ks0ij(iComp,jComp)+etaRef*ks1ij(iComp,jComp) )*chiS
			! ksij =  -0.04*18*nComps/ SUM(bRef) *chiS=  -0.04*18*nComps/ SUM(bRef) *vij(deliS-deljS)^2/298R =  
			!solParEntro(iComp)=SQRT( a0i*1.987d0*298.d0*etaStd/bRef(iComp) )
			! => ksij =  -0.04*18*nComps/ SUM(bRef) *vij(a0i/bi-a0j/bj)^2*298R*0.4/298R =  -0.04*0.4*(bi+bj)*(a0i/bi-a0j/bj)^2*18*nComps/ SUM(bRef) 
			!NOTE: genlzd value for ksij is -0.04, but only when i.ne.j!!!!!!!!!!!!!!!!!!!!!!
			!ksij=0 ! this is what Neil assumes.
			!_ASSERT( (a0i*a0j).gt.0 )
			aRefij=SQRT( bRef(iComp)*bRef(jComp)*a0i*a0j )/bRefMix  !NOTE: genlzd value for ksij is -0.04
			a0Mix=a0Mix+xFrac(iComp)*xFrac(jComp)*aRefij*(1.d0-ksij) !moved 1-ksij to here so dKijDeta works right for z0Mix JRE 2012-07
			if((ABS(a0i)<1e-22 .or. ABS(a0j)<1e-22) .and. LOUD)pause 'MixRule:a0i or a0j ~ 0. Watch for divide error.' 
			zRefij=SQRT( bRef(iComp)*bRef(jComp)/(a0i*a0j) ) /2.d0
			zRefij=zRefij*( z0i*a0j+a0i*z0j )/bRefMix*(1-ksij)	! check: a0j/sqrt(a0ij) and sqrt(bi*bj)/bMix cancel so units are zRef. 
			z0Mix=z0Mix+xFrac(iComp)*xFrac(jComp)*( zRefij+kETAij(iComp,jComp)*aRefij )
			aMixQuadPart(iComp)=aMixQuadPart(iComp)+xFrac(jComp)*aRefij*(1-ksij)
			!a0KijMix=a0KijMix+xFrac(iComp)*xFrac(jComp)*aRefij*dKijDeta 	!this is a new term b/c kij=kij(eta)
			!Ref part done. Do att part.
			!	bipKij=KU0IJ(iComp,jComp)+KU1IJ(iComp,jComp)*eta
			!Inserted by AV 8/5/07 to get rid of the SpeadGamma model
            bipKij=KIJ(iComp,jComp)+KTIJ(iComp,jComp)*eta !NOTE: nComps==1 omits this part by goto 100
			if(iComp==jComp.and.nComps==2)bipKij=0
			if(a1j > 0 .or. a2j > 0)then
				if(LOUD)write(*,*)'Mixrule: a1j or a2j .ge. 0. a1j,a2j',a1j,a2j
				if(LOUD)pause 'Check your parmstpt.'
				iErr=11
				exit
			endif
			a1ijSq=( bVoli*bVolj * a1i*a1j )
			if(a1ijSq < 0)then
				if(LOUD)pause 'Mixrule Error:sqArg < 0 for a1ijSq'
				iErr=13
				exit
			else
				a1ij= -SQRT(a1ijSq)*(1-bipKij)/bVolMix
			endif
							
			!a1ij=0.5*(a1i+a1j)*(1-bipKij)		!jre082304
			a1Mix=a1Mix+xFrac(iComp)*xFrac(jComp)*a1ij
			!z1ij=eta*da1ij/dEta
			!z1ij=(1-bipKij)*( etaDA1jDeta+etaDA1iDeta )/2								!jre082304
			if( (a1i/a1j) < 0 .and. LOUD)pause 'MixRule: a1i/a1j < 0'
			z1ij=(1-bipKij)*( z1j*SQRT(a1i/a1j)+z1i*SQRT(a1j/a1i) )/2.d0
			z1ij=z1ij*SQRT( bVoli*bVolj )/bVolMix !this should be same result as form for z0, but simpler connection to old molfrac basis. 
			z1Mix=z1Mix+xFrac(iComp)*xFrac(jComp)*(z1ij+kETAij(iComp,jComp)*a1ij)
			a1KijMix=a1KijMix+xFrac(iComp)*xFrac(jComp)*a1ij*kETAij(iComp,jComp) 	!this is a new term b/c kij=kij(eta)
			a2ijSq=( bVoli*bVolj * a2i*a2j )
			if(a2ijSq < 0)then
				if(LOUD)pause 'MixRule Error:sqArg < 0 for a2ijSq'
				iErr=13
				exit
			else
				a2ij=-SQRT(a2ijSq)*(1-bipKij)/bVolMix
			endif
			!a2ij=(a2i+a2j)/2*(1-bipKij)
			a2Mix=a2Mix+xFrac(iComp)*xFrac(jComp)*a2ij
			z2iPart =0
			if(a2i < 0)	z2iPart=z2i*SQRT(a2j/a2i)
			z2jPart =0
			if(a2j < 0)	z2jPart=z2j*SQRT(a2i/a2j)
			z2ij=(1-bipKij)*( z2jPart+z2iPart )/2  !jre072304  THIS IS WHY A2J CANNOT BE ZERO
			z2ij=z2ij*SQRT( bVoli*bVolj )/bVolMix !this should be same result as form for z0, but simpler connection to old molfrac basis. 
			z2Mix=z2Mix+xFrac(iComp)*xFrac(jComp)*(z2ij+kETAij(iComp,jComp)*a2ij)
			a2KijMix=a2KijMix+xFrac(iComp)*xFrac(jComp)*a2ij*kETAij(iComp,jComp) 	!this is a new term b/c kij=kij(eta)
   			aMixQuadPart(iComp)=aMixQuadPart(iComp)+xFrac(jComp)*(a1ij+a2ij/tKelvin)/tKelvin
		enddo  !jComp
	enddo !iComp
	!if(LOUD.and.init)print*,'MixRule: Calling all done except derivative calcs. isZiter = ',isZiter

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2010
!C Deleted by JRE: See FuTpt2d20221229.f90 for old code to compute derivatives of phys terms. Deprecated in favor of NumDervs for all.																							  C
!C Derivatives of ref. and perturbation terms in respect to eta and N											  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
100	if (isZiter== -1) then !compute derivatives analytically. I gave up on this...(?) JRE
		continue
	endif

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine AxsCalc(xFrac,eta,nComps,tKelvin,a0Mix,iErr) 
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*77 errMsg(11)
	!doublePrecision KIJ,KTIJ
	DIMENSION xFrac(NMX),aMixQuadPart(NMX),bRef(NMX)
	DIMENSION solParEntro(NMX),solParEnerg(NMX) !,zRefMix(5), !a1Mix(5),a2Mix(5),
	!dimension cRef(NMX,3),cRefMix(3)
	!dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	!common/ksvall/ks0ij,ks1ij
	!common/XsProps/a0XssS
	!common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!data initial/1/
	!common/ksvals/ks0(nmx,nmx),ks1(nmx,nmx)
	!data ks0,ks1/0, 0/		! best I could do based on avg of the entire system AV 11/29/08
	!vMolecNm3(1)=.0260
	!vMolecNm3(2)=.15734						! Note that polymer solutions are very sensitive to the BIPs and these values should be optimized via KD option
    do i=1,nComps
        bRef(i)=bVolRef(i,tKelvin)
    enddo
    bVolMix=SumProduct(nComps,xFrac,bVolCC_mol)
    bRefMix=SumProduct(nComps,xFrac,bRef)
  	etaStd=0.4
    etaRefStd=etaStd*bRefMix/bVolMix
	do iComp=1,nComps
		call TptTerms(isZiter,nComps,iComp,etaStd,etaRef,a0i,a1i,a2i,etaDA0iDeta,etaDA1iDeta,etaDA2iDeta,iErrTpt)
		!if(iComp==1)then
		!	A01=a0i
		!else
		!	A02=a0i
		!endif
		!if(iComp==1)a0i=3.1286
		!if(iComp==2)a0i=5.227
		solParEntro(iComp)=SQRT( a0i*1.987*298*etaStd/bRef(iComp) )
		solParEnerg(iComp)=SQRT(-a1i*1.987*etaStd/bVolCC_mol(i) )
		if(LOUD)write(*,'(a,i4,2(f8.2))')' i,delSi,delUi',iComp,solParEntro(iComp),solParEnerg(iComp)
	enddo
!	do iComp=1,nComps
!	call TptTerms(iComp,eta,a0i,a1i,a2i,etaDA0iDeta,etaDA1iDeta,etaDA2iDeta)
!		if(iComp==1)then
!			A01=a0i
!		else
!			A02=a0i
!		endif
!	enddo


	iErr=0
	errMsg(1)=' AxsCalc Error: a1i or a2i > 0'
	errMsg(2)=' AxsCalc Error: Wertheim returned error.'
	errMsg(3)=' AxsCalc Error: sqArg = a1i*a1j or a2i*a2j < 0'
	errMsg(4)=' AxsCalc Error: rMwMix .le. 0.  Check input.'
	errMsg(5)=' AxsCalc Error: eta out of range.'
	if(eta > 0.85 .or. eta < 0)then
		iErr=5
		if(LOUD)write(*,*)'AxsCalc: eta=',eta
		!call BeepMsg(errMsg(iErr))
		return
	endif

	if(nTptCoeffs.ne.4)nTptCoeffs=4;
	bVolMix=0  !must compute rMwMix,bVolMix first to get eta for A1,A2 mixrule.
	do iComp=1,nComps
		bVolMix=bVolMix+xFrac(iComp)*vMolecNm3(iComp)*602.22
	enddo
	if(bVolMix==0 .and. LOUD)pause 'AxsCalc: bVolMix=0'
	rhoMol_cc=eta/bVolMix 

	a0Mix=0
	z0Mix=0
	a1Mix=0
	a2Mix=0
	z1Mix=0
	z2Mix=0
	a0KijMix=0 ! new for kij=kij(eta), maybe not necessary.
	a1KijMix=0 ! new for kij=kij(eta)
	a2KijMix=0 ! new for kij=kij(eta)
	do iComp=1,nComps
		bVoli=bVolCC_mol(iComp)
		aMixQuadPart(iComp)=0
		call TptTerms(isZiter,nComps,iComp,eta,etaRef,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
		if(a0i.le.0 .or. a1i.ge.0 .or. a2i.ge.0)then
			if(LOUD)write(*,*)'FuTptVtot error:a0i<0 or a1i or a2i .ge. 0. a0i,a1i,a2i'
			if(LOUD)write(*,*)'   eta         a0i        a1i         a2i'
			if(LOUD)write(*,'(1x,f8.4,3e13.5)')eta,a0i,a1i,a2i
			if(LOUD)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)pause 'Check your input.'
			iErr=1
			exit !terminate the do loop, Note: "cycle" would go to enddo and keep looping.
		endif
		if(iErr.ne.0)then
			if(LOUD)pause
			return
		endif
		do jComp=1,nComps
			bVolj=bVolCC_mol(jComp)
			call TptTerms(isZiter,nComps,jComp,eta,etaRef,a0j,a1j,a2j,z0j,z1j,z2j,iErrTpt)
			vij=( bRef(iComp)+bRef(jComp) )/2
			chiS=vij/1.987/298*(solParEntro(iComp)-solParEntro(jComp))**2
			ksij=(ks0ij(iComp,jComp)+eta*ks1ij(iComp,jComp))*chiS
			!NOTE: genlzd value for ksij is -0.04, but only when i.ne.j!!!!!!!!!!!!!!!!!!!!!!
			!ksij=0 ! this is what Neil assumes.
			!_ASSERT( (a0i*a0j).gt.0 )
			aRefij=SQRT( bRef(iComp)*bRef(jComp)*a0i*a0j )*(1-ksij)/bRefMix  !NOTE: genlzd value for ksij is -0.04
			a0Mix=a0Mix+xFrac(iComp)*xFrac(jComp)*aRefij
		enddo  !jComp
	enddo !iComp
			!zRefij=SQRT( bVoli*bVolj/(a0i*a0j) ) /2
			!zRefij=zRefij*( z0i*a0j+a0i*z0j )/bVolMix*(1-ksij)	! check: a0j/sqrt(a0ij) and sqrt(bi*bj)/bMix cancel so units are zRef. 
			!dKijDeta=(-ks1ij(iComp,jComp)*eta)/(1-ksij) !really this is eta*dKijDeta/(1-ksij)
			!z0Mix=z0Mix+xFrac(iComp)*xFrac(jComp)*( zRefij+dKijDeta*aRefij )
			!aMixQuadPart(iComp)=aMixQuadPart(iComp)+xFrac(jComp)*aRefij
			!a0KijMix=a0KijMix+xFrac(iComp)*xFrac(jComp)*aRefij*dKijDeta 	!this is a new term b/c kij=kij(eta)
			!Ref part done. Do att part.
		!	bipKij=KU0IJ(iComp,jComp)+KU1IJ(iComp,jComp)*eta
		!Inserted by AV 8/5/07 to get rid of the SpeadGamma model
			!if(nComps==1)then
			!	bipKij=KII(iComp,jComp)+KTII(iComp,jComp)*eta
			!else
			!	bipKij=KU0IJ(iComp,jComp)+KU1IJ(iComp,jComp)*eta
			!endif
			!if(iComp==jComp.and.nComps==2)bipKij=0
			!if(a1j.ge.0 .or. a2j.ge.0)then
			!	if(LOUD)write(*,*)'FuTptVtot error:a1j or a2j .ge. 0. a1j,a2j',a1j,a2j
			!	if(LOUD)pause 'Check your input.'
			!	iErr=1
			!	exit
			!endif
			!a1ijSq=( bVoli*bVolj * a1i*a1j )
			!if(a1ijSq.lt.0)then
			!	if(LOUD)pause 'ZCalcTpt Error:sqArg < 0 for a1ijSq'
			!	iErr=3
			!	exit
			!else
			!	a1ij=-SQRT(a1ijSq)*(1-bipKij)/bVolMix
			!endif
			!a1ij=0.5*(a1i+a1j)*(1-bipKij)		!jre082304
			!a1Mix=a1Mix+xFrac(iComp)*xFrac(jComp)*a1ij
			!z1ij=eta*da1ij/dEta
			!z1ij=(1-bipKij)*( etaDA1jDeta+etaDA1iDeta )/2								!jre082304
			!z1ij=(1-bipKij)*( z1j*SQRT(a1i/a1j)+z1i*SQRT(a1j/a1i) )/2
			!z1ij=z1ij*SQRT( bVoli*bVolj )/bVolMix !this should be same result as form for z0, but simpler connection to old molfrac basis. 
			!dKijDeta=(-KU1IJ(iComp,jComp)*eta)/(1-bipKij) !really this is dKijDeta/(1-bipKij)
			!z1Mix=z1Mix+xFrac(iComp)*xFrac(jComp)*(z1ij+dKijDeta*a1ij)
			!a1KijMix=a1KijMix+xFrac(iComp)*xFrac(jComp)*a1ij*dKijDeta 	!this is a new term b/c kij=kij(eta)
			!a2ijSq=( bVoli*bVolj * a2i*a2j )
			!if(a2ijSq.lt.0)then
			!	if(LOUD)pause 'ZCalcTpt Error:sqArg < 0 for a2ijSq'
			!	iErr=3
			!	exit
			!else
			!	a2ij=-SQRT(a2ijSq)*(1-bipKij)/bVolMix
			!endif
			!a2ij=(a2i+a2j)/2*(1-bipKij)
			!a2Mix=a2Mix+xFrac(iComp)*xFrac(jComp)*a2ij
			!z2ij=(1-bipKij)*( z2j*SQRT(a2i/a2j)+z2i*SQRT(a2j/a2i) )/2  !jre072304
			!z2ij=z2ij*SQRT( bVoli*bVolj )/bVolMix !this should be same result as form for z0, but simpler connection to old molfrac basis. 
			!dKijDeta=(bipKij-KU0IJ(iComp,jComp))/(1-bipKij) !same as a1
			!z2Mix=z2Mix+xFrac(iComp)*xFrac(jComp)*(z2ij+dKijDeta*a2ij)
			!a2KijMix=a2KijMix+xFrac(iComp)*xFrac(jComp)*a2ij*dKijDeta 	!this is a new term b/c kij=kij(eta)

			!aMixQuadPart(iComp)=aMixQuadPart(iComp)+xFrac(jComp)*(a1ij+a2ij/tKelvin)/tKelvin
		!enddo  !jComp
	!enddo !iComp
	!a0Xsss=a0Mix-(xFrac(1)*A01+xFrac(2)*A02)
	!if(LOUD)write(*,*)'eta,Ais',eta,a0Mix,a1Mix,a2Mix
	!aRefKij=a0KijMix
	!aAttKij=(a1KijMix+a2KijMix/tKelvin)/tKelvin
	return
	end
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Written by AGF, Nov. 2009																					  C
!C	Having T,V and nMols, it calculates the Z factor and all derivatives which are needed in critical point,		  C
!C	flash and bubble point calculations.																		  C
!	Example 1.  nC7+benzene at 458K,0.9MPa, Hij=0, Kij=0, eta=0.435481025
!	id  xi     b	   Z0   Z1	     Z2 	NAS NDS	XA 	XD	     lnPhi
!	17  0.5	 69.05	11.084	-4592	42776	0	0	1	1      	-2.5705
!	501 0.5	 44.13	10.414	-3948	20788	0	0	1	1   	-1.9403
!	bMix	A0Mix	A1Mix	A2Mix	Z0Mix	Z1Mix	Z2Mix	zRep	zAtt	zFactor	
!	56.59	4.858	-3330	-11343	10.66	-4249	32165	10.656	11.617	0.0393
!	Example 2.  EtOH+water at 353K,0.1MPa, x1=0.5, Hij=0, Kij= -0.02, eta=0.8 (1st iteration)
!	id     xi  	   b	   Z0     Z1	Z2  	NAS NDS	bigYhb	  XA 	  XD	lnPhiAssoc
!	1102  0.5	27.24	173.6	-3722	9873	001	001	1912	0.016	0.016 	-33.03
!	1921  0.5	 8.94	240.0	-3645	9435	2	2	94.8	0.043	0.043   -20.70
!	bMix	A0Mix	A1Mix	A2Mix	Z0Mix	Z1Mix	Z2Mix	zRep	zAtt	zAssoc	zFactor	
!	18.09	27.8	-3301	-2937	184.9	-3477	9129	184.9	-9.78	-17.9	158.25
!	Example 2b.  EtOH+water at 353K,0.1MPa, x1=0.5, Hij=0, Kij= -0.02, eta=0.4654 (converged)
!	id  	xi    b	     Z0   	 Z1		Z2 		NAS NDS	bigYhb	  XA 	  XD	lnPhi
!	1102  0.5	27.24	????	????	????	001	001	????	0.??	0.?? 	  0.24375
!	1921  0.5	 8.94	????	????	????	2	2	????	0.??	0.??   	 -0.45329
!	zFactor = 0.0013; 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine FuTptVtot(isZiter,zFactor,aRes,uRes,ChemPo,vTotCc,tKelvin,gmol,nComps,iErr)
	! isZiter = 1 (omit aDep,uDep,lnPhi,derivatives), 0 (omit derivatives), -1 (compute derivatives and all) 
	USE GlobConst, only: aRes_RT,uRes_RT,HRes_RT,Sres_R
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs !NOTE: USE Assoc is not here because no assoc is computed.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!PARAMETER(etaStep=2.d0/1000.d0,maxRhoStep=1.d0/etaStep)
	character*77 errMsg(22)
	!doublePrecision moleStep,mShape(NMX)
	LOGICAL LOUDER
	!integer assocFlag
	DIMENSION xFrac(NMX),aMixQuadPart(NMX),FUGASSOC(NMX),bRef(NMX),gmol(NMX),chemPo(NMX) !,gmol_old(NMX)

	DIMENSION dFUGASSOC_dT(NMX),dFUGASSOC_dRHO(NMX) !,dh_dN(NMX),dFUGASSOC_dN(NMX,NMX)	  
	!DIMENSION dh_dN_num(NMX),dFUGASSOC_dN_num(NMX,NMX),fugAssocLoop(NMX)
	!DIMENSION dFUGREP_dN(NMX,NMX),dFUGATT_dN(NMX,NMX),dFugCof_dN(NMX,NMX),dFUGQUAD_dN(NMX,NMX),dzFactor_dN(NMX) 
	!DIMENSION dzAtt_dN(NMX),daAtt_dN(NMX),dzAssoc_dN_num(NMX),dbVolMix_dN(NMX),dbRefMix_dN(NMX) 

	!COMMON/z_aComps/zRefij,aRefij,dKijDeta,z1ij,a1ij,z2ij,a2ij
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	!COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
	!common/Derv1RG/dPRG_dV,d2PRG_dV2
	!COMMON/DervMix/dz0Mix_dN(NMX),da0Mix_dN(NMX),dz1Mix_dN(NMX),da1Mix_dN(NMX),dz2Mix_dN(NMX),da2Mix_dN(NMX),daMixQuadPart_dN(NMX,NMX),daMixQuadPart_deta(NMX)
	!COMMON/DervMix2/dz0Mix_deta,da0Mix_deta,dz1Mix_deta,da1Mix_deta,dz2Mix_deta,da2Mix_deta,d2z0Mix_deta2,d2z1Mix_deta2,d2z2Mix_deta2
	!COMMON/DervMix3/d2a1Mix_deta2,d3a1Mix_deta3,d2a2Mix_deta2
	!COMMON/DervZMix3/d3z0Mix_deta3,d3z1Mix_deta3,d3z2Mix_deta3
	!COMMON/ETA2/ETA
	!COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	!COMMON/dfug2/dFUGASSOC_dN_num,dh_dN_num
	!COMMON/dFug/h_nMichelsen,dh_dT,dh_dN,dFugassoc_dT,dFugassoc_dN,dFUGASSOC_dRHO,dh_dRHO
	!COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	!COMMON/A/a0Mix,a1Mix,a2Mix,aAtt,a_d,zRep
	!common/aRG/zNC,aDepNC
	!common/rg/iFlagRG,iFlagFit,mShape,cutOff(NMX),Phi(NMX),Zi(NMX)
	data initCall/1/
	LOUDER=LOUD
	!LOUDER=.TRUE.

	iErr=0
	errMsg(11)=' FuTptVtot Error: Nonsense input. e.g. eta > 1'
	!errMsg(4)=' FuTptVtot Error: a1i or a2i > 0'
	errMsg(12)=' FuTptVtot Error: Wertheim returned error.'
	errMsg(13)=' FuTptVtot Error: Error from MixRule. Terminal, sorry.'
	errMsg(5)=' FuTptVtot warning: T < tKmin for at least one component.'
	totMoles=SUM(gmol)
	bVolMix=0
	bRefMix=0
    TcVolatile=123456
	isAcid=0
    DO i=1,nComps
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        if(Tc(i) < TcVolatile)then
            iVolatile=i
            TcVolatile=Tc(i)
        endif
		do j=1,nTypes(i)
			if(idType(i,j)==1603)isAcid=1
		enddo 
    enddo
   if(tKelvin < tKmin(iVolatile)*0.998D0)then	   !round TKmin down a little for possible lost precision
        iErr=5
        if(LOUDER.and.initCall)write(*,'(a,i3,2f7.1,a)')' FuVtot: iVolatile,T(K),TKmin= ',iVolatile,tKelvin,tKmin(iVolatile),TRIM( errMsg(iErr) )
    endif
    hRes_RT=86.8686D0 !initialize to error indicator in case of bad return
    uRes_RT=86.8686D0 !initialize to error indicator in case of bad return
    Ares_RT=86.8686D0 !initialize to error indicator in case of bad return
    Sres_R =86.8686D0 !initialize to error indicator in case of bad return
    zFactor=86.8686D0 !initialize to error indicator in case of bad return
    if(iErr > 10)return	  ! let iErr==5 be a warning. e.g. liquid density may be ok although Psat would be bogus.  
    !bVolPure=bVolMix
	rhoMol_cc=totMoles/vTotCc
	eta=rhoMol_cc*bVolMix
	eta2=eta*eta
	eta3=eta2*eta
	eta4=eta3*eta
	if( eta > etaMax)then
		if(LOUDER)write(*,*) 'FuTptVtot: eta > etaMax. eta,etaMax=',eta,etaMax
		if(LOUDER)pause
        iErr=11
        return
	endif

	!call XsAs(nComps,iErrXs)  !for debugging
	!print*,'FuTptVtot: calling MixRule'
	call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
	if(iErrMix > 10)then
		iErr=13
        if(LOUDER)write(*,*)'FuTptVtot: iErrMix=',iErrMix
		!call BeepMsg(errMsg(iErr))
		return
	endif 
	zRep=z0Mix
	aTvRef=a0Mix
	aAtt=( a1Mix+   a2Mix/tKelvin )/tKelvin 	 ! => a3 = a2*g0*exp{-g1*[(eta-eta0)^2-eta0^2]}/3; g0[=]K
	zAtt=( z1Mix+   z2Mix/tKelvin )/tKelvin	 ! 
	uAtt=( a1Mix+ 2*a2Mix/tKelvin )/tKelvin
	if(bGaussEx)then 
		Gf=1
		dGf_dEta=0
		d2Gf_dEta2=0
		uDepGf=0
		!if(nComps==1)call GaussFun3(tKelvin,eta,Gf,dGf_dEta,d2Gf_dEta2,cpGf,uDepGf,iErr)
		!denom=(1.d0-eta)*(1.d0-eta)*(1.d0-eta)
		!Adep=A0+A1/T+A2/T^2*Gf=> Z-1 = Z0 + Z1/T + (eta*dA2/dEta*Gf + A2*eta*dGf_dEta)/T^2
		!& dZ_dEta = dZ_dEta0 + dZ_dEta1/T + (dA2/dEta*Gf + A2*dGf_dEta)/T^2+ (eta*d2A2/dEta2*Gf +2*eta*dA2/dEta*dGf_dEta + A2*eta*d2Gf_dEta2)/T^2
		!P = Z*rho*R*T => dP_dRho = Z*R*T + rho*R*T*dZ_dRho = R*T*(Z+eta*dZ_dEta)
		!dRho/rho = -dV/V => V*dP_dV = -rho*dP_dRho => dP_dV = -rho^2*dP_dRho  
		z2Gex=z2Mix*Gf+a2Mix*eta*dGf_dEta
		aAtt=( a1Mix+ a2Mix*Gf/tKelvin )/tKelvin 	 ! => a3 = a2*g0*exp{-g1*[(eta-eta0)^2-eta0^2]}/3; g0[=]K
		zAtt=( z1Mix+ z2Gex/tKelvin )/tKelvin	 ! 
		uAtt=( a1Mix+ 2.d0*a2Mix/tKelvin*Gf+a2Mix/tKelvin*uDepGf )/tKelvin 
	endif
	!a_d=a2Mix/a1Mix/TC(1)
	
	!C We define a flag for association term here to avoid calling Werthiem routine when it is not necessary...
	!print*,'FuTptVtot: calling Wertheim'
    zAssoc=0
    aAssoc=0
    uAssoc=0
    if(iErr < 11)call Wertheim(isZiter,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,fugAssoc,iErrCode)
!	if( ABS(zAssoc) < 1e-11 )print*,'zAssoc ~ 0???'		 !for debugging
	IF(iErrCode)then
		!iErr=iErrCode
		if(LOUDER)pause 'FuTptVtot: Wertheim gave error.'
		iErr=12
		return
	endif
	iDoSolvation=0 ! all of aBipDA is zero as of 20221226, so no need to do solvation. 
	if ( ABS(zAssoc) > zeroTol .and. isZiter==0) then
		if(iDoSolvation==1)call WertheimFugc(xFrac,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
	endif

	if(iErr > 10)then
		if(LOUDER)write(*,*)errMsg(iErr)
		return
	endif
	
	!if(LOUDER)print*,'uAtt,uAssoc',uAtt,uAssoc,uDep
	zFactor=(1.d0+zRep+zAtt+zAssoc) !+zNC
	PMpa=(zFactor)*rGas*rhoMol_cc*tKelvin
	if(isZiter > 0)return 
	aRes=aTvRef+aAtt+aAssoc !+aDepNC
	uRes=uAtt+uAssoc
	
	aRes_RT=aRes
	uRes_RT=uRes
	sRes_R =-(aRes-uRes) 
	hRes_RT=uRes+zFactor-1.d0
	DO iComp=1,nComps
		bVolk=bVolCC_mol(iComp) 
		bRefk=bVolRef(iComp,tKelvin) !Function bVolRef()...
		FUGREP=bRefk/bRefMix*(z0Mix-a0Mix)
		!FUGREP=bVolk/bVolMix*(z0Mix-a0Mix)
		FUGATT=bVolk/bVolMix*(zAtt - aAtt)
		!aAttQuadPart=aAttQuadPart+etaPow*aMixQuadPart(kComp,iCoeff)
		FUGQUAD=2.d0*aMixQuadPart(iComp) !we lump all the quadparts together in the current version.  
		FUGBON=fugAssoc(iComp)
		chemPo(iComp)=FUGREP+FUGATT+FUGQUAD+fugAssoc(iComp) !+(zNC-aDepNC)*bVolk/bVolMix+2.d0*aDepNC 
		if(LOUDER)write(*,'(a,5E12.4)')' FuVtot:fugrep,fugatt,fugquad,fugbon,Z',fugrep,fugatt,fugquad,fugbon,zFactor
	enddo

	if(isZiter < 0)then
		if(LOUDER.and.initCall)print*,'FuTptVtot: T,rho,z=',tKelvin,rhoMol_cc,zFactor
        if(LOUDER.and.zFactor < zeroTol) then
            pause 'FuTptVtot: zFactor < 0?'
            continue
        endif
		call NUMDERVS(nComps,xFrac,tKelvin,rhoMol_cc,zFactor,cmprsblty,iErrNum)
		if(LOUDER.and.initCall)print*,'FuTptVtot: cmprsblty,CvRes_R=',cmprsblty,CvRes_R
		if(iErrNum > 1.and.LOUDER)then
			print*,'FuTptVtot: NUMDERVS failed. No derivatives for you!'
			iErr=16
		endif
	endif
	initCall=0
	return
	end	 ! subroutine FuTptVtot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine XsAs(nComps,iErr)
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	character*77 errMsg(11)
	!DIMENSION XA(NMX,maxTypes),XD(NMX,maxTypes),XC(NMX,maxTypes)
	dimension a0Pure(2),a1Pure(2),a2Pure(2),xFrac(NMX),aMixQuadPart(NMX)
	!common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	errMsg(1)='XsAs Error: This routine permits only 2 components.'
	iErr=0
	if(nComps.ne.2)then
		iErr=1
		if(LOUD)write(*,*)errMsg(1)
		if(LOUD)pause
		return
	endif
	if(LOUD)write(*,*)'Input eta for xs Ai calculations'
	read(*,*)eta
	tKelvin=300 !tKelvin should not matter for the Ai values, they are only density dependent.
	do iComp=1,2
		xFrac(1)=0
		xFrac(2)=0
		xFrac(iComp)=1
		call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Pure(iComp),a1Pure(iComp),a2Pure(iComp),z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
	enddo
	do iFrac=0,10,1
		xFrac(1)=float(iFrac)/10.0
		xFrac(2)=1-xFrac(1)
		call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
		a0xs=a0Mix-( xFrac(1)*a0Pure(1)+xFrac(2)*a0Pure(2) )
		a1xs=a1Mix-( xFrac(1)*a1Pure(1)+xFrac(2)*a1Pure(2) )
		a2xs=a2Mix-( xFrac(1)*a2Pure(1)+xFrac(2)*a2Pure(2) )
		if(LOUD)write(*,'(2f8.4,2f10.3,2E11.4)')a0Mix,a0Xs,a1Mix,a1Xs,a2Mix,a2Xs
	enddo
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DoublePrecision Function bVolRef(i,tKelvin)
	USE GlobConst !, only:ID 
	USE SpeadParms, only:isHeNeH2
    Implicit DoublePrecision(A-H,O-Z)
	!Purpose: Enable softening of the reference contribution using the E&D FPE'86 correlation.
	DoublePrecision nEff ! simplifies the equations while accurately conveying the name.
    dEhs_sigmaCube=1
    if(isHeNeH2(i))then  !isHeNeH2(i) is part of GlobConst
        nEff=9	 ! H2(902),D2(925),Ne(919)
		if(id(i)==913.or.id(i)==923)nEff = 15	  !He3,He4
		p2 = (.0093d0*nEff-.0592d0)*nEff !E&D(1986)
		p2 = ( (nEff-6)/6 )**2 /exp(  gammln( (nEff-3)/nEff )  )**(2*nEff/3) ! ESK(2019)
		!print*,'bVolRef: p2 = ', p2
		!pause
		p1 = (nEff-1)  
		p1 = 1.13*nEff-3.28 !ESK(2019) 
		eps_kB=2	!Helium
		if(id(i)==902.or.id(i)==925.or.id(i)==919)eps_kB=15	 ! H2,D2,Ne
        Tstar=tKelvin/eps_kB !eps/kB for softness correction.
		rMin_sigmaCube=(nEff/6)**(3/(nEff-6)) 
		E0=1  !E&D
		E0=0 !ESK     
        dEhs_sigmaCube=( p2*Tstar*Tstar+p1*Tstar+1 )**( -3.D0/(2*nEff+E0) )*rMin_sigmaCube ! JRE FPE, 86. rMin_Sigma^3=(n/6)^( 3/(n-6) ) = 1.5 if n=9.
    endif
    bVolRef=bVolCC_mol(i)*dEhs_sigmaCube
    return
    end
      FUNCTION gammln(xx)  ! returns the natural log of the gamma function.
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5/
	  data stp/2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine ChemPoPhysNum(chemPo,eta,tKelvin,gmol,nComps,iErr)
	!Purpose: Eval the chemPo of phys part by numerical derivative to check analytical.
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE SpeadParms !GlobConst+Assoc(XA,XD,XC)+AiCoeffs etc.
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!DIMENSION XA(NMX),XD(NMX),chemPoAssoc(NMX)!,aAttk(NMX)
	dimension xFracPert(NMX),chemPo(NMX),xFrac(NMX),gmol(NMX),bRef(NMX)
	DIMENSION aMixQuadPart(NMX) !,zRefMix(5)!,a1Mix(5),a2Mix(5)
	!common/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	!COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	!common/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT

	iErr=0
	totMoles=SUM(gmol) ! SUM() is an intrinsic function
	if(ABS(totMols) < 0.001 .and. LOUD)pause 'FuTpt: totMols ~0? Will try renormalizing, but you should check input.'
	bVolMix=0
	bRefMix=0
    DO i=1,nComps
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
	ENDDO
	!Basis: 1mole to start. Central diff.
	redStep=1e-4
	do k=1,nComps
		do i=1,nComps
			dn=0
			if(i.eq.k)dn=redStep
			xFracPert(i)=xFrac(i)*(1-dn)/(1-redStep)
		enddo
		etaMinus=eta*(1-redStep*bVolCC_mol(k)/bVolMix)
		call MixRule(isZiter,xFracPert,tKelvin,etaMinus,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
		if(iErrMix > 10 .and. LOUD)pause 'ChemPoPhysNum: fatal error from MixRule.'
		aMinus=a0Mix+(a1Mix+a2Mix/tKelvin)/tKelvin
		zMinus=1+z0Mix+(z1Mix+z2Mix/tKelvin)/tKelvin
		do i=1,nComps
			dn=0
			if(i.eq.k)dn=redStep
			xFracPert(i)=xFrac(i)*(1+dn)/(1+redStep)
		enddo
		etaPlus=eta*(1+redStep*bVolCC_mol(k)/bVolMix)
		call MixRule(isZiter,xFracPert,tKelvin,etaPlus,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
		if(iErrMix > 10 .and. LOUD)pause 'ChemPoPhysNum: fatal error from MixRule.'
		aPlus=a0Mix+(a1Mix+a2Mix/tKelvin)/tKelvin
		zPlus=1+z0Mix+(z1Mix+z2Mix/tKelvin)/tKelvin
		zMid=(zPlus+zMinus)/2
		chemPo(k)=( aPlus*(1+redStep)-aMinus*(1-redStep) )/2/redStep - DLOG(zMid)
	enddo

	return
	end
	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Programmed by AFG 2011.																						  C
!C This part till the end of the file is not used anywhere in the CalcEos.										  C
!C I used these routines to compare PC-SAFT with SPEAD11 														  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE PCSAFT(tKelvin,eta,zFactor)
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION etaList(21)
	DOUBLEPRECISION mShape
	nEta=21
 	!tKelvin=469.70d0
	mShape=2.002 !29.7282 !2.002  !7.9849 !4.6627d0 !2.6896d0			  
	sigma=3.6184 !4.2787d0 !3.6184d0 !3.8384d0  !3.7729d0
	eok=208.11 !293.608d0 !208.11 !257.75d0 !231.2d0			 
	!PC=2.140d6  !Pa
	!TC=613.70d0 !K
	etaList(1)=0.00001
	DO I=2,nEta
		etaList(I)=etaList(I-1)+0.03
	ENDDO
	!DO I=1,nEta   !AUG, this routine just works for options 4 and 5, we need to add lines for other options as well
		!eta=etaList(I)
		call zPcSaftCalc(tKelvin,eta,mShape,sigma,eok,zFactor,aRes,a_d)
		!if(LOUD)write(*,*)eta,'  ',a_d
	!ENDDO
	return 
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  	SUBROUTINE zPcSaftCalc(tKelvin,eta,mShape,sigma,eok,zFactor,aRes,a_d)
	USE GlobConst ! for LOUD
	IMPLICIT DOUBLEPRECISION(a-h,o-z)
	!INTEGER NITER 
	DOUBLEPRECISION mShape,mReverse,mdReverse,md2Reverse,md3Reverse
	DOUBLEPRECISION ap(7,3),bp(7,3),am(7),bm(7)
	DOUBLEPRECISION I1,I2
	dEhs=sigma !*(1-0.12*EXP(-3*eok/tKelvin))
	md=mShape*dEhs
	md2=mShape*dEhs*dEhs
	mReverse=1/mShape
	mdReverse=mReverse/dEhs
	md2Reverse=mdReverse/dEhs
	md3Reverse=md2Reverse/dEhs
	kBoltzman= 8.31447/6.02214d23
	rho=6/pi*eta*md3Reverse
	Xitha0=pi/6*rho*mShape
	Xitha1=Xitha0*dEhs
	Xitha2=Xitha1*dEhs

	ap(1,1)=  0.91056314451539d0 !from p1234 of ref1
	ap(1,2)= -0.30840169182720d0
	ap(1,3)= -0.09061483509767d0
	ap(2,1)=  0.63612814494991d0
	ap(2,2)=  0.18605311591713d0
	ap(2,3)=  0.45278428063920d0
	ap(3,1)=  2.68613478913903d0
	ap(3,2)= -2.50300472586548d0
	ap(3,3)=  0.59627007280101d0
	ap(4,1)= -26.5473624914884d0
	ap(4,2)=  21.4197936296668d0
	ap(4,3)= -1.72418291311787d0
	ap(5,1)=  97.7592087835073d0
	ap(5,2)= -65.2558853303492d0
	ap(5,3)= -4.13021125311661d0
	ap(6,1)= -159.591540865600d0
	ap(6,2)=  83.3186804808856d0
	ap(6,3)=  13.7766318697211d0
	ap(7,1)=  91.2977740839123d0
	ap(7,2)= -33.7469229297323d0
	ap(7,3)= -8.67284703679646d0

	bp(1,1)=  0.72409469413165d0
	bp(1,2)= -0.57554980753450d0
	bp(1,3)=  0.09768831158356d0
	bp(2,1)=  1.11913959304690d0  *2.d0
	bp(2,2)=  0.34975477607218d0  *2.d0
	bp(2,3)= -0.12787874908050d0  *2.d0
	bp(3,1)= -1.33419498282114d0  *3.d0
	bp(3,2)=  1.29752244631769d0  *3.d0
	bp(3,3)= -3.05195205099107d0  *3.d0
	bp(4,1)= -5.25089420371162d0  *4.d0
	bp(4,2)= -4.30386791194303d0  *4.d0
	bp(4,3)=  5.16051899359931d0  *4.d0
	bp(5,1)=  5.37112827253230d0  *5.d0
	bp(5,2)=  38.5344528930499d0  *5.d0
	bp(5,3)= -7.76088601041257d0  *5.d0
	bp(6,1)=  34.4252230677698d0  *6.d0
	bp(6,2)= -26.9710769414608d0  *6.d0
	bp(6,3)=  15.6044623461691d0  *6.d0
	bp(7,1)= -50.8003365888685d0  *7.d0
	bp(7,2)= -23.6010990650801d0  *7.d0
	bp(7,3)= -4.23812936930675d0  *7.d0

	am(1)=ap(1,1) + (mShape-1)/mShape*ap(1,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(1,3)
	am(2)=ap(2,1) + (mShape-1)/mShape*ap(2,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(2,3)
	am(3)=ap(3,1)+ (mShape-1)/mShape*ap(3,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(3,3)
	am(4)=ap(4,1) + (mShape-1)/mShape*ap(4,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(4,3)
	am(5)=ap(5,1) + (mShape-1)/mShape*ap(5,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(5,3)
	am(6)=ap(6,1) + (mShape-1)/mShape*ap(6,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(6,3)
	am(7)=ap(7,1) + (mShape-1)/mShape*ap(7,2) + (mShape-1)/mShape*(mShape-2)/mShape*ap(7,3)

	bm(1)=bp(1,1) + (mShape-1)/mShape*bp(1,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(1,3)
	bm(2)=bp(2,1) + (mShape-1)/mShape*bp(2,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(2,3)
	bm(3)=bp(3,1) + (mShape-1)/mShape*bp(3,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(3,3)
	bm(4)=bp(4,1)+ (mShape-1)/mShape*bp(4,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(4,3)
	bm(5)=bp(5,1) + (mShape-1)/mShape*bp(5,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(5,3)
	bm(6)=bp(6,1) + (mShape-1)/mShape*bp(6,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(6,3)
	bm(7)=bp(7,1) + (mShape-1)/mShape*bp(7,2) + (mShape-1)/mShape*(mShape-2)/mShape*bp(7,3)      
	order1=mShape*mShape*sigma**3*eok/tKelvin
	order2=mShape*mShape*sigma**3*(eok/tKelvin)**2	  

!***************************
!		eta Vapor
!		b2Virial=2.5+1.5*mShape-12*md3Reverse*am(1)*order1-6/dEhs3*bm(1)*order2
!		eta=(-1+ SQRT(1+(4*pMpa*pi/6*md2*dEhs*b2Virial)/8.314/tKelvin))/2/b2Virial
!***************************

	voidFrac=1.d0-eta
	voidFrac2=voidFrac*voidFrac
	voidFrac3=voidFrac*voidFrac2
	voidFrac4=voidFrac*voidFrac3
	voidFrac5=voidFrac*voidFrac4
	voidFrac6=voidFrac*voidFrac5

	I1=am(1)
	I2=bm(1)
	dEtaI1_dEta=am(1)
	dEtaI2_dEta=bm(1)
	etaProd=eta
	do iCoeff=2,7
		I1=I1+etaProd*am(iCoeff)
		I2=I2+etaProd*bm(iCoeff)
		dEtaI1_dEta=dEtaI1_dEta + iCoeff*am(iCoeff)*etaProd 
		dEtaI2_dEta=dEtaI2_dEta + iCoeff*bm(iCoeff)*etaProd 
		etaProd=etaProd*eta
	enddo

	denom=voidFrac*(2-eta)
	denom2=denom*denom
	denom3=denom*denom2
	IF(voidFrac.le.0 .and. LOUD)PAUSE 'voidFrac<0'
	!zHs=eta/(1-eta) + 3*Xitha1*Xitha2/Xitha0/(1-eta)**2+(3*Xitha2**3-eta*Xitha2**3)/Xitha0/(1-eta)**3
	zHs=4.d0*eta*(1.d0-eta/2.d0)/voidfrac3	  !it is ZHS-1
	aHs=1/Xitha0* (3*Xitha1*Xitha2/voidFrac+Xitha2**3/eta/(1-eta)**2+(Xitha2**3/(eta*eta)-Xitha0)*DLOG(voidFrac))
	!gHs=1/(1-eta)+(dEhs/2)*3*Xitha2/(1-eta)**2+(dEhs/2)**2*2*Xitha2**2/voidFrac3
	gHs=(1.d0-eta/2.d0)/voidfrac3
	rhogHs_rho=eta*(2.5d0-eta)/voidfrac4
	!rhogHs_rho=eta/(1-eta)**2+dEhs/2*(3*Xitha2/(1-eta)**2+6*Xitha2*eta/(1-eta)**3)+(dEhs/2)**2*(4*Xitha2**2/(1-eta)**3+6*Xitha2**2*eta/(1-eta)**4)
	zHc=mShape*(zHs) - (mShape-1.d0)*(rhogHs_rho/gHs)
	IF(gHs.le.0 .and. LOUD)PAUSE 'gHs<0'
	aHc=mShape*aHs - (mShape-1)*DLOG(gHs)

	C1inv=1+mShape*(8*eta-2*eta*eta)/voidFrac4 + &
		(1-mShape)*eta*(  20+eta*(-27+eta*(12-2*eta) )  )/denom2
	C1=1/C1inv
	C2=-(C1*C1)* (mShape*(-4*eta**2+20*eta+8)/voidFrac5+(1-mShape)*&
	 (2*eta**3+12*eta**2-48*eta+40)/denom3   )
	zDisp=-2*pi*rho*dEtaI1_dEta*order1-pi*rho*mShape*&
	 (C1*dEtaI2_dEta+C2*eta*I2)*order2
	aDisp=-2*pi*rho*I1*order1-pi*rho*mshape*C1*I2*order2	  ! mShape*mShape*sigma**3*(eok)**2

	!a_d=(zHc)/80*(1.d0-eta)**3 ! pi*rho*mshape*C1*I2*order2/(2*pi*rho*I1*order1)
	!a_d=-12*mShape*eta*I1*eok/5
	a_d=(-6*eta*mShape*mShape*C1*I2*(eok)**2)/80
	zFactor=1+zHc+zDisp
	aRes=aHc+aDisp !+aAssoc

	return 
	END 

