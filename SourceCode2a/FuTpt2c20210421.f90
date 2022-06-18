	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetTptCas(nComps,idCas,iErrCode)
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
	USE Assoc !GlobConst+XA,XD,XC
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000)
	character errMsgPas*77
	!	integer GetBIPs
	integer idComp(nComps),idCas(nComps) !localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	call IdDipprLookup(nComps,idCas,ier,errMsgPas)
	call GetTpt(nComps,idComp,iErrCode)
	return
	end



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetTpt(nComps,idComp,ierrCode)
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
	USE Assoc !GlobConst+XA,XD,XC
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(ndb=1555,listPool=1000)
	integer GetBIPs
	integer idComp(nComps) !,ier(11) !localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!integer localType(maxTypesGlobal),idLocalType(maxTypes)  !idLocalType points back to localType for double-linking
	doublePrecision bondRate(nmx,maxTypes) !local to this subroutine.
	!dimension zRefDb(5),a1Db(5),a2Db(5)
	character bipFile*88,bipHbFile*88,ParmsHbFile*88,ParmsTptFile*88,dumString*333
	character*123 ErrMsg(11)
	dimension nFg(nmx,maxTypes) !not in common because we transcribe these into HbParms and they do not need to be passed otherwise.
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	common/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
    common/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	common/siteParm/nFg	   !This is needed in GetRG
	common/ppVpCoeffs/vpCoeffs(NMX,5)
	isTPT=.TRUE. ! in GlobConst simplifies calls in FUGI and FuVtot
	etaMax=1-zeroTol
	iErrCode=0
	ErrMsg(1)='GetTpt error - error reading ParmsTpt.txt. Path? Debug?'
	ErrMsg(2)='GetTpt error - error reading ParmsHb4.txt. Path?'
	ErrMsg(3)='GetTpt error - error reading ParmsTpt.txt at line='
	ErrMsg(4)='GetTpt error - nTypesTot > maxTypes'
	ErrMsg(5)='GetTpt error - a1 > 0 for at least one component when 0 < eta < 0.75'
	ErrMsg(6)='GetTpt error - a2 > 0 for at least one component when 0 < eta < 0.75'
	!if(LOUD)print*,'Debug=', DEBUG
	IF(DEBUG)then 
		ParmsTptFile='c:\SPEAD\CalcEos\input\ParmsTpt.txt'	 !use this form internal dev since spead\ParmsTpt is most reliable.
		if(iEosOpt==14)ParmsTptFile='c:\SPEAD\CalcEos\input\ParmsTptTransPGL6ed.txt'	 !use this form internal dev since spead\ParmsTpt is most reliable.
	ELSE 
		ParmsTptFile=TRIM(masterDir)//'\input\ParmsTpt.txt' ! // is the concatenation operator
		if(iEosOpt==14)ParmsTptFile=TRIM(masterDir)//'\input\ParmsTptTransPGL6ed.txt' ! // is the concatenation operator
		!ParmsTptFile='G:\My Drive\MYPROJEX\UaWrapper\UaWrapper'//'\input\ParmsTpt.txt' ! // is the concatenation operator
		!if(iEosOpt==13)ParmsTptFile='G:\My Drive\MYPROJEX\UaWrapper\UaWrapper'//'\input\ParmsTptTransPGL6ed.txt' ! // is the concatenation operator
		!ParmsTptFile='ParmsTpt.txt' !use this when preparing a release version.
	ENDIF

	if(LOUD)write(*,*)'ParmsTptFile=',TRIM(ParmsTptFile)
	OPEN(40,FILE=ParmsTptFile,ERR=861)

	READ(40,*,ERR=861)nDeck	 !check that there is something to read.
    if(nDeck < 1 .and. LOUD)pause 'nDeck<1 in TptParms.txt.  Please check.'

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
			read(40,'(A333)',ERR=863,END=111)dumString
			read(dumString,*,ERR=863)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),&
				vMolecNm3(iComp),tKmin(iComp),rMw(iComp),nTypes(iComp),(nFg(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				bVolCC_mol(iComp)=vMolecNm3(iComp)*AvoNum
111			continue  !transfer from read(40,...,END=111)dumString
			!read(40,*,ERR=861)idBase,idCc,(zRefDb(iCoeff),iCoeff=1,3),(a1Db(iCoeff),iCoeff=1,nTptCoeffs),(a2Db(iCoeff),iCoeff=1,nTptCoeffs),vMolecDb,tKmin(iComp),nTypes(iComp),(idType(iComp,iType),iType=1,nTypes(iComp)),(nFg(iComp,iType),iType=1,nTypes(iComp))
			IF(idBase.EQ.idComp(iComp))THEN
				iGotIt=1  !this will kick us to the next component

				if(LOUD)write(*,'(a,5f13.4)')' Mw,bVol(cc/mol),vEff(nm3)',rMw(iComp),bVolCc_mol(iComp),vMolecNm3(iComp)
				if(LOUD)write(*,'(a,5e13.5)')' zRefCof',(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
				if(LOUD)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUD)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
				if(LOUD)write(*,'(a,i3,9(i3,i5))')' nTypes,Nd,Id',nTypes(iComp),(nFg(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
				if(LOUD)write(*,'(a,i3,<nTypes(iComp)>i5)')' nTypes,#Fgi',nTypes(iComp),(nFg(iComp,iType),iType=1,nTypes(iComp))
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
				if(LOUD)write(*,*)' bVol(Tc)/bVolCC_mol =  ' ,etaFactor
                if(LOUD)write(*,*)'  eta     zRef     A1/100    A2/10000   rho(g/cc)   zRefCs'
                do i=5,85,5
                    eta=i/1.D2
                    etaRef=eta*etaFactor
					rhoG_cc=eta*rMw(iComp)/bVolCC_mol(iComp)
					ZrefCs = 4*eta*(1-eta/2)/(1-eta)**3
	                call TptTerms(isZiter,nComps,iComp,eta,etaRef,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
                    if(LOUD)write(*,'(f7.4,6f10.5)')eta,z0i,a1i/100,a2i/(100*100),rhoG_cc,zRefCs
                    if(a1i > 0)iErrCode=5
                    if(a2i > 0)iErrCode=6
                enddo
                if(iErrCode .and. LOUD)pause 'GetTpt: check Ais. A1 or A2 > 0 ?'
				exit !quit searching if found
			ENDIF !idBase==idComp
		enddo !while(iGotIt.eq.0)
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			if(LOUD)write(*,*)'Error in GetTpt: ParmsTpt.txt has no data for component #:',idComp(iComp)
		endif
	enddo
	if(nTypesTot.gt.maxTypes)then
		iErrCode=4
		if(LOUD)write(*,*)ErrMsg(4)
		if(LOUD)pause
	endif
	CLOSE(40)
        
	if(iErrCode.ne.0)return

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG, 2011
!C This function replaces the SPEAD parameters with those of SPEAD11 correlation
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	if (iEosOpt.Eq.9) then
		CALL swap(nComps,nFg)
	endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
 	!ref part done.  do att part.
				
	IF(DEBUG)then 
		!ParmsHbFile='c:\Spead\CalcEos\Input\ParmsHb3.txt'
		ParmsHbFile='c:\Spead\CalcEos\Input\ParmsHb4.txt'
	ELSE 
		!ParmsHbFile=TRIM(masterDir)//'\input\ParmsHb3.txt' ! // is the concatenation operator
		ParmsHbFile=TRIM(masterDir)//'\input\ParmsHb4.txt' ! // is the concatenation operator
	ENDIF
    inHbFile=40
	if(LOUD)write(*,*)'ParmsHbFile=',TRIM(ParmsHbFile)
	OPEN(inHbFile,FILE=ParmsHbFile)

	READ(inHbFile,*,ERR=862)nDeck	 !check that there is something to read.
	if(nDeck < 1 .and. LOUD)pause 'nDeck<1 in ParmsHb4.txt.  Please check.'

	!we do not store the entire database then operate on it because that would take a lot of space.
	!instead, we rewind and re-read it from the hard drive multiple times.  this happens only at startup.
	DO iComp=1,nComps	!Note: assume all types hbond. It just means summing over some zeroes.
		!nHbTypes(iComp)=nTypes(iComp) 
		bondSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType)*(nAs+nDs) ) !so no increment for nAs=nDs=0.
		numSitesTot=0	!this total is for just this molecule. = sum( iType, nFg(iComp,iType) )
		do iType=1,nTypes(iComp)
			numSitesTot=numSitesTot+nFg(iComp,iType)
			isHbType=0
			eDonorKcal_mol(iComp,iType)=0        
			eAcceptorKcal_mol(iComp,iType)=0
			eHbKcal_mol(iComp,iType)=0  
			bondVolNm3(iComp,iType)=0 !note: bondVol is same for same type on diff components, but keeping separate for each component allows the potential of adding a comp specific accessibility factor at a later date.
			nDegree(iComp,iType)=0
			nDonors(iComp,iType)=0
			nAcceptors(iComp,iType)=0
			dHkcalMol(iComp)=0
			bondVolNm3Esd(iComp)=0
			ND(iComp)=0
			NDS(iComp)=0
			NAS(iComp)=0

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
					nDegree(iComp,iType)=nFg(iComp,iType)
					nDonors(iComp,iType)=ndsBase
					nAcceptors(iComp,iType)=nasBase
					bondSitesTot=bondSitesTot+nFg(iComp,iType)*(ndsBase+nasBase)
					bondVolNm3(iComp,iType)=bondVolDb
					bondRate(iComp,iType)=bondRateDb
					dHkcalMol(iComp)=(dHDonorKcalDb+dHAcceptorKcalDb)/2.d0
					bondVolNm3Esd(iComp)=bondVolDb
					ND(iComp)=nFg(iComp,iType)
					NDS(iComp)=ndsBase
					NAS(iComp)=nasBase
					exit !terminate the search for this type because you found it.
				ENDIF
			enddo !loop over ParmsHb4.txt
			!if(iGotIt.eq.0)iErrCode=iErrCode*10+iComp
		enddo !loop over iType
		iBondExp=2
		bondVolNm3Esd(iComp)=bondVolNm3Esd(iComp)*( 1+bondRate(iComp,1)*(bondSitesTot/numSitesTot)**iBondExp )
		!above is for a single bSite/molec.  Below is more general, but not implemented yet.
		do iHbType=1,nTypes(iComp)
			bondVolNm3(iComp,iHbType)=bondVolNm3(iComp,iHbType)*( 1+bondRate(iComp,iHbType)*( nFg(iComp,iHbType)/DFLOAT(numSitesTot) )**iBondExp )
            if(bondVolNm3(iComp,iHbType) < 0)bondVolNm3(iComp,iHbType)=2D-8 !this only happens for formic acid. JRE 20200303.
		enddo
		if(LOUD)write(*,*)' SiteType:',( idType(iComp,i),i=1,nTypes(iComp))
		if(LOUD)write(*,'(a,11i4)')' nAcceptors:',( nAcceptors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUD)write(*,'(a,11i4)')' nDonors:   ',( nDonors(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUD)write(*,'(a,9E11.4)')' bondVolNm3:',( bondVolNm3(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUD)write(*,'(a,11f7.1)')' Don energy:',( eDonorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
		if(LOUD)write(*,'(a,11f7.1)')' Acc energy:',( eAcceptorKcal_mol(iComp,iType),iType=1,nTypes(iComp) )
        if(LOUD)pause 'Check HB parms.'
	enddo !iComp
	CLOSE(inHbFile)

	nC=nComps

!  note:  bips are passed back through common/BIPs/
	IF(DEBUG)then 
		bipFile='c:\spead\CalcEos\input\BipSpead.txt'
	ELSE 
		bipFile=TRIM(masterDir)//'\input\BipSpead.txt' ! // is the concatenation operator
	ENDIF
	iErrCode=GetBIPs(bipFile,ID,nC)
	!Note: no need to check for switching because kij and ktij are symmetric
	if(iErrCode.ne.0 .and. LOUD)pause 'GetTpt Error: GetBIPs returned error.'

	!load aBipAD matrix
	IF(DEBUG)then 
		bipHbFile='c:\Spead\CalcEos\Input\BipAD.txt'
	ELSE 
		bipHbFile=TRIM(masterDir)//'\input\BipAD.txt' ! // is the concatenation operator
	ENDIF
	call GetAssocBips(bipHbFile,aBipAD,ierABip) !in WertheimVv.f90. idLocalType,nTypesTot in common/assoc 

	!load aBipDA matrix
	IF(DEBUG)then 
		BipHbFile='c:\Spead\CalcEos\Input\BipDA.txt'
	ELSE 
		BipHbFile=TRIM(masterDir)//'\input\BipDA.txt' ! // is the concatenation operator
	ENDIF
	call GetAssocBips(bipHbFile,aBipDA,ierABip) !in WertheimVv.f90

	if(iEosOpt.eq.8)CALL GetVp(nComps,ID,iErrVp,vpCoeffs)	

	RETURN
861	continue
	iErrCode=1
	if(LOUD)write(*,*)ErrMsg(iErrCode)
	if(LOUD)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUD)pause
	return                      
862	continue
	iErrCode=2
	if(LOUD)write(*,*)ErrMsg(2)
	if(LOUD)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUD)pause
	return                      
863	continue
	iErrCode=3
	if(LOUD)write(*,*)ErrMsg(3),jComp
	if(LOUD)write(*,*)dumString
	read(dumString,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
	if(LOUD)write(*,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3)
	if(LOUD)pause 'ids and zRef read ok. Checking A1.'
	read(dumString,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUD)write(*,*)'A1',(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUD)pause 'A1 read ok. Checking A2.'
	read(dumString,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUD)write(*,*)'A2',(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs)
	if(LOUD)pause 'A2 read ok. Checking V,Tmin,Mw'
	read(dumString,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),vMolecNm3(iComp),tKmin(iComp),rMw(iComp)
	if(LOUD)write(*,*)'V,Tmin,Mw',vMolecNm3(iComp),tKmin(iComp),rMw(iComp)
	if(LOUD)pause 'VTM read ok. Checking nTypes, nFg, idType'
	read(dumString,*)idBase,idCc,idCAS,(zRefCoeff(iComp,iCoeff),iCoeff=1,3),(a1Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),(a2Coeff(iComp,iCoeff),iCoeff=1,nTptCoeffs),vMolecNm3(iComp),tKmin(iComp),rMw(iComp),nTypes(iComp),(nFg(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
	if(LOUD)write(*,*)'Types',nTypes(iComp),(nFg(iComp,iType),idType(iComp,iType),iType=1,nTypes(iComp))
	if(LOUD)write(*,*)'nDeck,iCompo',NDECK,jComp
	if(LOUD)pause
	return                      
END	! GetTpt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG, 2010																C
!C This function replaces the SPEAD parameters with those of SPEAD11 correlation	C
!C Please note that it has been optimized for n-alkanes								C
!C Also, its applicability to mixtures has not been tested, yet.					C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE swap(nComps,nFg)
	USE Assoc !parameters + nTypes to compute nCarbons
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DOUBLEPRECISION mShape(NMX)
	DIMENSION nCarbons(NMX),nFg(NMX,maxTypes)
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	common/FileSwap/isFileSwap,isSw
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
	write(*,*)'Enter 1 to swap coeffs from ...calceos\input\PrismTptAcoeffs.txt or 0 to use AFG Spead11'
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
END	! Subroutine SWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QueryParPureSpead(iComp,iParm,value,iErr)
	USE EsdParms      !Just for ESD
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	!DoublePrecision EOKP(NMX),KCSTAR(NMX),DH(NMX),C(NMX),Q(NMX),VX(NMX)
	!Integer         ND(NMX),NDS(NMX),NAS(NMX)
	iErr=0
	if(iParm==1)then
		value=c(iComp)
	elseif(iParm==2)then
		value=vx(iComp)
	elseif(iParm==3)then
		value=eokp(iComp)
	elseif(iParm==4)then
		value=KcStar(iComp)
	elseif(iParm==5)then
		value=DH(iComp)
	elseif(iParm==6)then
		value=ND(iComp)
	elseif(iParm==7)then
		value=NAS(iComp)
	elseif(iParm==8)then
		value=NDS(iComp)
	else
		iErr=1
	endif
	return
end	!Subroutine QueryParPureSpead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPureSpead(iComp,iParm,value,iErr)
	USE EsdParms      !Just for ESD
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	!DoublePrecision EOKP(NMX),KCSTAR(NMX),DH(NMX),C(NMX),Q(NMX),VX(NMX)
	!Integer         ND(NMX),NDS(NMX),NAS(NMX)
	iErr=0
	if(iParm==1)then
		c(iComp)=value
		q(iComp)=1+(value-1)*1.9076D0
	elseif(iParm==2)then
		vx(iComp)=value
	elseif(iParm==3)then
		eokp(iComp)=value
	elseif(iParm==4)then
		KcStar(iComp)=value
	elseif(iParm==5)then
		DH(iComp)=value
	elseif(iParm==6)then
		ND(iComp)=value
	elseif(iParm==7)then
		NAS(iComp)=value  ! Not used for ESD96
	elseif(iParm==8)then
		NDS(iComp)=value  ! Not used for ESD96
	else
		iErr=1
	endif
	return
end	! SetParPureSpead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																												  C
!C Modified by AFG 2011 : isZiter is redefined																							  C
!C																												  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuTpt(tKelvin,pMPa,gmol,nComps,LIQ,chemPo,zFactor,ier)
	USE Assoc !XA,XB,XC, + GlobConst{Tc,Pc,...RGAS,AvoNum,...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	character*111 errMsg(22)
	integer assocFlag,ier(11)
	DIMENSION xFrac(NMX),gmol(NMX),chemPo(NMX)
	DIMENSION zRefMix(5),bRef(NMX)!,aMixQuadPart(NMX),a1Mix(5),a2Mix(5)
	COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	COMMON/DEPFUN/dU_RT,dA_RT,DSONK,dH_RT
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

	do index=1,11
		ier(index)=0
	enddo      
	nTptCoeffs=4+1;
	do iCoeff=1,nTptCoeffs
		zRefMix(iCoeff)=0.0
	enddo
	totMoles=1 !we may assume that totMoles=1, since sum(xFrac)=1
	totMoles=SUM(gmol)
	bVolMix=0
	bRefMix=0
	if(initCall.and.LOUD)print*,'FuTpt: T(K) = ',tKelvin
    isTtooLow=0
	TcVolatile=12345
    DO i=1,nComps
		if(Tc(i) < TcVolatile)then
			TcVolatile=Tc(i)
			iVolatile=i
		endif
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        chemPo(i)=86.8686D0 !initialize to error indicator in case of bad return
    enddo
	if( tKelvin < tKmin(iVolatile) )isTtooLow=1
    if(isTtooLow)then
		iErr=5
		if(LOUD.and.initCall)write(*,*)'FuTpt:iVolatile,TKmin= ',iVolatile,TKmin(iVolatile),TRIM( errMsg(iErr) )
		if(LOUD)write(*,*)TRIM(errMsg(iErr))
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
		if(LOUD)write(*,*)' rGas=',rGas,' tKelvin=',tKelvin
		call BeepMsg(errMsg(iErr))
		goto 86
	endif
	Pb_RT=pMPa*bVolMix/RGAS/tKelvin
	P_RT=pMPa/RGAS/tKelvin
	pb_rt=P_RT*bVolMix

	!  GUESS FOR eta
	eta=pb_rt/1.001D0 ! Make "old" value imprecise, so "new" value (below) will be ideal gas result at low rho. => precision when P -> 1E-11.
	IF(LIQ==1.or.LIQ==3.or.eta>0.75)eta=0.75D0
	etaRef=eta*bRefMix/bVolMix
	IF (eta > 0.9999 .or. etaRef > 0.9999 .or. eta < 0)then
		iErr=11
		if(LOUD)write(*,*)' pMPa=',pMPa,' eta=',eta,' etaRef=',etaRef
		!call BeepMsg(errMsg(iErr))
		goto 86
	endif
	rhoMol_cc=eta/bVolMix
	vTotCc=totMoles/rhoMol_cc
	isZiter=1
	if(initCall.and.LOUD)print*,'FuTpt: first call to FuVtot. T(K),eta=',tKelvin,eta
	call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErrZ)
	if(iErrZ > 10)then
		iErr=13
		if(LOUD)write(*,*)'FuTpt: FuTptVtot returned with error code:',iErrZ
		!call BeepMsg(errMsg(iErr))
		goto 86
	endif
	etaOld=eta
	errOld=pb_rt-eta*zFactor
	eta=etaOld*1.001D0 !set this "new" value to ideal gas value so if change < 1E-9 on first iteration, we get exactly the ideal gas result.
    !if(LIQ==1.or.LIQ==3)eta=0.8 !Take a big step for liquids to bracket answer sooner. Saves about 4/13 iterations
	!ETA=rhoMol_Cc*bVolMix	!calculate here before entering loop, then at the end with each new rhoMol_Cc
	change=1234
	itMax=55
	NITER=0
	pMax=-1.d11 !for crude goldenz
	pMin=1.d11
	isZiter=1
    iterGold=0
	if(LOUD.and.initCall)print*,'FuTpt init: liq,eta', liq,eta
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
			if(LOUD)write(*,*) 'Error in FuTpt: eta > 1. nIter=',nIter
			do iComp=1,nComps
				if(LOUD)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
				if(LOUD)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			enddo
			if(LOUD)pause
		endif
		rhoMol_cc=eta/bVolMix
		vTotCc=totMoles/rhoMol_cc
		if(initCall.and.LOUD)print*,'FuTpt: iterative call to FuVtot. iter,eta=',NITER,eta
		call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErrZ)
		if(iErrZ > 10)then
			iErr=13
			if(LOUD)write(*,*)'FuTpt: FuTptVtot returned with error code:',iErrZ
			!call BeepMsg(errMsg(iErr))
			cycle
		endif
		error=pb_rt-eta*zFactor
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
		if(eta < 0 .or. eta > 0.9999)then !NOTE: this should not happen given step limiting
			if(niter < (itMax-1))niter=itMax-1 !restrict tries with new guess so we can crash if necessary.
			if(LOUD)write(*,'(a,f8.4)')' Warning in FuTpt: next guess is eta=',eta
			eta=0
			if(liq==1)eta=0.98
			if(LOUD)write(*,*) 'Restarting iteration with eta=',eta
			if(LOUD)pause
		endif
	enddo  !iteration on eta to find P=P(input)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Whoohoo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(iErr.ne.0)then
        iterGold=iterGold+1
		if(liq==0 .or. LIQ==2)eta=etaAtPmax  !crude goldenz
		if(liq==1.or.LIQ==3)then
			eta=etaAtPmin*1.01  !crude goldenz
			rhoMol_cc=eta/bVolMix
			vTotCc=totMoles/rhoMol_cc
			call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErrZ)
			if(LOUD)write(*,'(a,i3,1x,a,f7.4,a,f7.4,a,f7.4)')' FuTpt: LIQ=',LIQ,' etaMin=',etaAtPmin,' pMin=',pMin,' pMax=',pMax
            if(iterGold < 2)goto 100
		endif
	endif
	IF (eta < 0)then
		iErr=12
		!call BeepMsg(errMsg(iErr))
	endif
	IF (zFactor < 0)then
		iErr=14
		IF(LOUD)call BeepMsg(errMsg(iErr))
	endif
	!call one more time to facilitate debugging
	rhoMol_cc=eta/bVolMix
	vTotCc=totMoles/rhoMol_cc
	zFactor=pMPa/(rhoMol_cc*rGas*tKelvin) ! Seemingly unnecessary, this should improve precision when computing rho from return. Z = 1+zRef+zAtt+zAssoc is subject to roundoff when Z->1E-9, but Z=P/(rhoRT)=>rho=P/(ZRT) with precision.
	IF (LIQ==1 .or. LIQ==3) THEN
	  ETAL=ETA
	  etaPass=eta  ! etaPass in GlobConst
	  ZL=zFactor  
	ELSE
	  ETAV=ETA
	  etaPass=eta  ! etaPass in GlobConst
	  ZV=zFactor  
	ENDIF
	dU_RT=uDep
	dH_RT=uDep+zFactor-1
	dA_RT=aDep
	if(LOUD.and.initCall)print*,'FuTpt: done. T(K),eta=',tKelvin,eta
	initCall=0
    if(LIQ > 1)return
	isZiter=0
	if(LOUD.and.initCall)print*,'FuTpt: Calculating fugacity. LIQ = ', LIQ
	call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErrZ)
	if(iErr.eq.2 .or. iErr.eq.4)goto 86
   	
	!  ITERATION ON RHO HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.

	!call ChemPoPhysNum(chemPo,eta,tKelvin,xFrac,nComps,iErrCp)
	if(LIQ < 2)call ChemPoCalcTpt(chemPo,eta,tKelvin,xFrac,nComps,iErrCp)
	ChemPoPure=aDep+zFactor-1-DLOG(zFactor)
	if(nComps==1 .and. LOUD .and. initCall)print*,'fugc(1),Gdep1',chemPo(1),ChemPoPure
	
	RETURN

86	if(LOUD)WRITE(6,'(a)')' ERROR IN FuTpt.  '
	IER(1)=iErr
	RETURN
867 continue
	if(LOUD)write(*,*)'Exiting FuTpt: Please correct xFrac'
	ier(1)=17
	ier(7)=1
	return
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CODED BY AV 06/26/06. cf. DannerGess/SPEADCI paper. IECR08
! PURPOSE: DEVELOPING SpeadGamma MODEL. cf.SPEADCI paper IECR, 47:7955 (08).
	SUBROUTINE FuSpeadGamma(tKelvin,pMPa,xFrac,nComps,LIQ,chemPo,zFactor,ier)
	USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION xFrac(NMX),chemPo(NMX),fugcPure(NMX),fugcL(nmx),rLnGam(nmx),x(nmx),fugc(nmx),ier(11)
	common/ppVpCoeffs/vpCoeffs(NMX,5)
	common/eta/etaL,etaV,zL,zV
	!common/EosOpt/masterDir,iEosOpt,initEos,DEBUG
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
		call FuTpt(tKelvin,pMPa,x,nc,1,fugcL,zL,ier)
		fugcPure(iPure)=fugcL(iPure)
		!hPure(iPure)=DHONKT
	enddo
	call FuTpt(tKelvin,pMPa,xFrac,nComps,1,fugcL,zL,ier)

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
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	character*77 errMsg(11)
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!DIMENSION dEta_dN(NMX),dEtaRef_dN(NMX)
	!COMMON/dervZ/dz0i_dN(NMX),dz1i_dN(NMX),dz2i_dN(NMX)
	!COMMON/derva/da1i_dN(NMX),da2i_dN(NMX),da0i_dN(NMX)
	COMMON/dervZ2/dz0i_deta,dz1i_deta,dz2i_deta,d2z0i_deta2,d2z1i_deta2,d2z2i_deta2
	COMMON/dervZ3/d3z0i_deta3,d3z1i_deta3,d3z2i_deta3
	COMMON/derva2/da1i_deta,da2i_deta,da0i_deta,d2a1i_deta2,d3a1i_deta3,d2a2i_deta2

	iErr=0
	NC=nComps
	errMsg(1)='TptTerms Error: etaRef out of range.'
	if(etaRef > 0.9999 .or. etaRef < 0)then
		iErr=1
		if(LOUD)write(*,*)'In TptTerms: eta=',eta
		!call BeepMsg(errMsg(iErr))
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
	dA0iDeta=(c1+etaRef*(c2+etaRef*c3))/void3  ! =zRef	 
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
		if(LOUD)write(*,'(a,i4,E11.4,3e12.4)')' Error in TptTerms: a1 or a2 > 0.  iComp, eta, a1i, a2i=',iComp,eta,a1i,a2i
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
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	common/GaussFun/G0,G1,eta0
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
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*77 errMsg(22)
	!doublePrecision KIJ,KTIJ
	!dimension cRef(NMX,3),cRefMix(3)
	DIMENSION xFrac(NMX),aMixQuadPart(NMX),gmol(NMX)
	DIMENSION solParEntro(NMX),solParEnerg(NMX) !,zRefMix(5), !a1Mix(5),a2Mix(5),
	DIMENSION ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	DIMENSION KII(NMX,NMX),KTII(NMX,NMX)
	DIMENSION bRef(NMX)
	DIMENSION dSPEntro_deta(NMX),dSPEnerg_deta(NMX),dEta_dN(NMX),dEtaRef_dN(NMX),dbVolMix_dN(NMX),dbRefMix_dN(NMX)

	COMMON/z_aComps/zRefij,aRefij,dKijDeta,z1ij,a1ij,z2ij,a2ij
	COMMON/SIPs/KII,KTII
	!COMMON/dervZ/dz0i_dN(NMX),dz1i_dN(NMX),dz2i_dN(NMX)
	COMMON/dervZ2/dz0i_deta,dz1i_deta,dz2i_deta,d2z0i_deta2,d2z1i_deta2,d2z2i_deta2
	COMMON/dervZ3/d3z0i_deta3,d3z1i_deta3,d3z2i_deta3
	!COMMON/derva/da1i_dN(NMX),da2i_dN(NMX),da0i_dN(NMX)
	COMMON/derva2/da1i_deta,da2i_deta,da0i_deta,d2a1i_deta2,d3a1i_deta3,d2a2i_deta2
	COMMON/VCoeff/B2,B3
	COMMON/DervMix/dz0Mix_dN(NMX),da0Mix_dN(NMX),dz1Mix_dN(NMX),da1Mix_dN(NMX),dz2Mix_dN(NMX),da2Mix_dN(NMX),daMixQuadPart_dN(NMX,NMX),daMixQuadPart_deta(NMX)
	COMMON/DervMix2/dz0Mix_deta,da0Mix_deta,dz1Mix_deta,da1Mix_deta,dz2Mix_deta,da2Mix_deta,d2z0Mix_deta2,d2z1Mix_deta2,d2z2Mix_deta2
	COMMON/DervMix3/d2a1Mix_deta2,d3a1Mix_deta3,d2a2Mix_deta2
	COMMON/DervZMix3/d3z0Mix_deta3,d3z1Mix_deta3,d3z2Mix_deta3
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!common/XsProps/a0Xs
	!data initial/1/
	! AGF. It should be 0 not -0.044  
	data init,ks0,ks1/1, -0.04, 0/		! best I could do based on avg of the entire system AV 11/29/08
							! Note that polymer solutions are very sensitive to the BIPs and these values should be optimized via KD option
	iErr=0
	errMsg(11)=' MixRule Error: a1i or a2i > 0'
	errMsg(12)=' MixRule Error: Wertheim returned error.'
	errMsg(13)=' MixRule Error: sqArg = a1i*a1j or a2i*a2j < 0'
	errMsg(14)=' MixRule Error: rMwMix .le. 0.  Check input.'
	errMsg(15)=' MixRule Error: eta out of range.'
	errMsg(5)=' MixRule Error: T < TKmin.'
	errMsg(17)=' MixRule Error: Error returned from TptTerms().'

	totMoles=SUM(gmol)
	bVolMix=0
	bRefMix=0
    DO i=1,nComps
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        if(tKelvin < tKmin(i)*0.998)then   ! round down to enable NUMDERVS
            iErr=5
            if(LOUD)write(*,*)'iComp= ',i,TRIM( errMsg(iErr) )
        endif
        !chemPo(iComp)=86 !initialize to error indicator in case of bad return !No chemPo here !!!
    enddo
	if(bVolMix==0 .and. LOUD)pause 'MixRule error: bVolMix=0'
    etaRef=eta*bRefMix/bVolMix
    rhoMol_cc=eta/bVolMix
	!if(LOUD)write(*,*)'init,ks0=',init,ks0
		init=1  !Sorry. These values aren't being stored for some reason. Better add to Tpt Module?
	if(init.ne.0)then
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
	endif

	if(eta > 0.9999 .or. eta < 0)then
		iErr=15
		if(LOUD)write(*,*)'MixRule: eta=',eta
		!if(LOUD)call BeepMsg( TRIM(errMsg(iErr)) )
		return
	endif
	if(etaRef > 0.9999 .or. etaRef < 0)then
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
			if(LOUD)write(*,*)'MixRule error in TptTerms:a0i<0 or a1i or a2i > 0. a0i,a1i,a2i'
			if(LOUD)write(*,*)'   eta         a0i        a1i         a2i'
			if(LOUD)write(*,'(1x,E11.4,3e13.5)')eta,a0i,a1i,a2i
			if(LOUD)write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)pause 'Check your input.'
			iErr=11
			exit !terminate the do loop, Note: "cycle" would go to enddo and keep looping.
		endif
		if (nComps==1) then
            bipKii=KII(1,1)+KTII(1,1)*eta
            !bipKii=KII(1,1)-KTII(1,1)*TC(1)/tKelvin
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
			if(ABS(a0i)<1e-22 .or. ABS(a0j)<1e-22 .and. LOUD)pause 'Mixrule:a0i or a0j ~ 0. Watch for divide error.' 
			zRefij=SQRT( bRef(iComp)*bRef(jComp)/(a0i*a0j) ) /2.d0
			zRefij=zRefij*( z0i*a0j+a0i*z0j )/bRefMix*(1-ksij)	! check: a0j/sqrt(a0ij) and sqrt(bi*bj)/bMix cancel so units are zRef. 
			z0Mix=z0Mix+xFrac(iComp)*xFrac(jComp)*( zRefij+dKijDeta*aRefij )
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
			z1Mix=z1Mix+xFrac(iComp)*xFrac(jComp)*(z1ij+dKijDeta*a1ij)
			a1KijMix=a1KijMix+xFrac(iComp)*xFrac(jComp)*a1ij*dKijDeta 	!this is a new term b/c kij=kij(eta)
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
			z2Mix=z2Mix+xFrac(iComp)*xFrac(jComp)*(z2ij+dKijDeta*a2ij)
			a2KijMix=a2KijMix+xFrac(iComp)*xFrac(jComp)*a2ij*dKijDeta 	!this is a new term b/c kij=kij(eta)
   			aMixQuadPart(iComp)=aMixQuadPart(iComp)+xFrac(jComp)*(a1ij+a2ij/tKelvin)/tKelvin
		enddo  !jComp
	enddo !iComp
	!if(LOUD.and.init)print*,'MixRule: Calling all done except derivative calcs. isZiter = ',isZiter

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2010																							  C
!C Derivatives of ref. and perturbation terms in respect to eta and N											  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
100	if (isZiter==0) then
		DO I=1,nComps
			deta_dN(I)=eta*bVolCC_mol(I)/bVolMix
			dbVolMix_dN(I)=bVolCC_mol(I)-bVolMix
			detaRef_dN(I) =etaRef*bRef(I)/bRefMix
			dbRefMix_dN(I)=bRef(I)-bRefMix
		ENDDO
		do iComp=1,nComps
			call TptTerms(isZiter,nComps,iComp,etaStd,etaRefStd,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
			solParEntro(iComp)=SQRT( a0i*1.987d0*298.d0*etaStd/bRef(iComp) )
			solParEnerg(iComp)=SQRT(-a1i*1.987d0*etaStd/bVolCC_mol(iComp) )
			dSPEntro_deta(iComp)=0.5d0*1.987d0*298.d0*etaRefStd/bRef(iComp)*da0i_deta/solParEntro(iComp)
			dSPEnerg_deta(iComp)=-0.5d0*1.987d0*etaStd/bVolCC_mol(iComp)*da1i_deta/solParEnerg(iComp)
		enddo

		rhoMol_cc=eta/bVolMix 
		a0KijMix=0 ! new for kij=kij(eta), maybe not necessary.
		a1KijMix=0 ! new for kij=kij(eta)
		a2KijMix=0 ! new for kij=kij(eta)
		dz0Mix_deta=0.d0
		dz1Mix_deta=0.d0
		dz2Mix_deta=0.d0
		da0Mix_deta=0.d0
		da1Mix_deta=0.d0
		da2Mix_deta=0.d0
		do iComp=1,nComps
			bVoli=bVolCC_mol(iComp)
			daMixQuadPart_deta(iComp)=0.d0
			call TptTerms(isZiter,nComps,iComp,eta,etaRef,a0i,a1i,a2i,z0i,z1i,z2i,iErrTpt)
			if (nComps.eq.1) then
				dz0Mix_deta=dz0i_deta !dz0i_dEta comes from common/derv...
				dz1Mix_deta=dz1i_deta
				dz2Mix_deta=dz2i_deta
				d2z0Mix_deta2=d2z0i_deta2
				d2z1Mix_deta2=d2z1i_deta2
				d2z2Mix_deta2=d2z2i_deta2
				d3z0Mix_deta3=d3z0i_deta3
				d3z1Mix_deta3=d3z1i_deta3
				d3z2Mix_deta3=d3z2i_deta3

  				da1Mix_deta=da1i_deta
				d2a1Mix_deta2=d2a1i_deta2
				d3a1Mix_deta3=d3a1i_deta3 
				
				da2Mix_deta=da2i_deta
				d2a2Mix_deta2=d2a2i_deta2
				return 
			endif
			dz0i_old_deta=dz0i_deta
			da0i_old_deta=da0i_deta
			dz1i_old_deta=dz1i_deta
			da1i_old_deta=da1i_deta
			dz2i_old_deta=dz2i_deta
			da2i_old_deta=da2i_deta
			do jComp=1,nComps
				bVolj=bVolCC_mol(jComp)
				call TptTerms(isZiter,nComps,jComp,eta,etaRef,a0j,a1j,a2j,z0j,z1j,z2j,iErrTpt)
				dz0j_deta=dz0i_deta
				dz0i_deta=dz0i_old_deta
				da0j_deta=da0i_deta
				da0i_deta=da0i_old_deta

				dz1j_deta=dz1i_deta
				dz1i_deta=dz1i_old_deta
				da1j_deta=da1i_deta
				da1i_deta=da1i_old_deta

				dz2j_deta=dz2i_deta
				dz2i_deta=dz2i_old_deta
				da2j_deta=da2i_deta
				da2i_deta=da2i_old_deta

		   		vij=(bRef(iComp)+bRef(jComp))/2.d0
				chiS=vij*chiConst*(solParEntro(iComp)-solParEntro(jComp))**2
				ksij=(ks0ij(iComp,jComp)+eta*ks1ij(iComp,jComp))*chiS
				dchiS_deta=vij*chiConst*(dSPEntro_deta(iComp)-dSPEntro_deta(jComp))**2
				dksij_deta=ks1ij(iComp,jComp)*chiS+eta*ks1ij(iComp,jComp)*dchiS_deta

				aRefij=SQRT( bRef(iComp)*bRef(jComp)*a0i*a0j )*(1.d0-ksij)/bVolMix  !NOTE: genlzd value for ksij is -0.04
				zRefij=SQRT( bRef(iComp)*bRef(jComp)/(a0i*a0j) ) /2.d0
				daRefij_deta=(bRef(iComp)*bRef(jComp)*(da0i_deta*a0j+da0j_deta*a0i))/(2.d0*SQRT( bRef(iComp)*bRef(jComp)*a0i*a0j ))*(1.d0-ksij)/bVolMix+(-dksij_deta/bVolMix)*SQRT( bRef(iComp)*bRef(jComp)*a0i*a0j )
				da0Mix_deta=da0Mix_deta+xFrac(iComp)*xFrac(jComp)*daRefij_deta
				dzRef=(1/2.d0)*(-bVoli*bVolj*(da0i_deta*a0j+da0j_deta*a0i)/((a0i*a0j)*(a0i*a0j)))/(2.d0*SQRT( bVoli*bVolj/(a0i*a0j)))
				zRefij=zRefij*( z0i*a0j+a0i*z0j )/bVolMix*(1-ksij)	! check: a0j/sqrt(a0ij) and sqrt(bi*bj)/bMix cancel so units are zRef. 
				dzRefij_deta=dzRef*( z0i*a0j+a0i*z0j )/bVolMix*(1-ksij)+(SQRT( bRef(iComp)*bRef(jComp)/(a0i*a0j) ) /2.d0)*( ( z0i*a0j+a0i*z0j )*(-dksij_deta/bRefMix)+((1-ksij)/bRefMix)*((dz0i_deta*a0j+da0j_deta*z0i)+(da0i_deta*z0j+dz0j_deta*a0i)) )
				dKijDeta=(-ks1ij(iComp,jComp)*eta)/(1.d0-ksij) !really this is eta*dKijDeta/(1-ksij)
				d2KijDeta2=(-ks1ij(iComp,jComp)*(1.d0-ksij)+dksij_deta*(-ks1ij(iComp,jComp)*eta))/(1.d0-ksij)**2
				dz0Mix_deta=dz0Mix_deta+xFrac(iComp)*xFrac(jComp)*( dzRefij_deta+d2KijDeta2*aRefij+dKijDeta*daRefij_deta )
				if(nComps==1)then
					dbipKij_deta=KTII(iComp,jComp)
				else
					dbipKij_deta=KTIJ(iComp,jComp)
				endif
				if(iComp==jComp.and.nComps==2)bipKij=0

				a1ijSq=( bVoli*bVolj * a1i*a1j )
				a1ij=-SQRT(a1ijSq)*(1-bipKij)/bVolMix
				!if(LOUD)print*,'MixRule: i,j,a1ij,a2ij = ',iComp,jComp,a1ij,a2ij
				da1ij_deta=-(bVoli*bVolj*(da1i_deta*a1j+da1j_deta*a1i))/(2.d0*SQRT(a1ijSq))*(1-bipKij)/bVolMix-SQRT(a1ijSq)*(-dbipKij_deta/bVolMix)
				da1Mix_deta=da1Mix_deta+xFrac(iComp)*xFrac(jComp)*da1ij_deta
				z1ij=(1-bipKij)*( z1j*SQRT(a1i/a1j)+z1i*SQRT(a1j/a1i) )/2.d0
				dz1=(-dbipKij_deta/2.d0)*( z1j*SQRT(a1i/a1j)+z1i*SQRT(a1j/a1i) )+((1-bipKij)/2.d0)*(dz1j_deta*SQRT(a1i/a1j)+z1j*(da1i_deta*a1j-da1j_deta*a1i)/(a1j**2)/(2.d0*SQRT(a1i/a1j))+dz1i_deta*SQRT(a1j/a1i)+z1i*(da1j_deta*a1i-da1i_deta*a1j)/(a1i**2)/(2.d0*SQRT(a1j/a1i)))
				z1ij=z1ij*SQRT( bVoli*bVolj )/bVolMix  
				dz1ij_deta=dz1*SQRT( bVoli*bVolj )/bVolMix
				dKijDeta=(-KTIJ(iComp,jComp)*eta)/(1-bipKij) !really this is dKijDeta/(1-bipKij)
				d2KijDeta2=(-KTIJ(iComp,jComp)*(1-bipKij)+dbipKij_deta*(-KTIJ(iComp,jComp)*eta))/((1-bipKij)*(1-bipKij))
				dz1Mix_deta=dz1Mix_deta+xFrac(iComp)*xFrac(jComp)*(dz1ij_deta+d2KijDeta2*a1ij+dKijDeta*da1ij_deta)
				a1KijMix=a1KijMix+xFrac(iComp)*xFrac(jComp)*a1ij*dKijDeta 	!this is a new term b/c kij=kij(eta)
				!if(LOUD)print*,'MixRule: i,j,dz1ij_deta = ',iComp,jComp,dz1ij_deta
				
				a2ijSq=( bVoli*bVolj * a2i*a2j )
				a2ij=-SQRT(a2ijSq)*(1-bipKij)/bVolMix
				if(a2ij .ge. 0)cycle
				da2ij_deta=-(bVoli*bVolj*(da2i_deta*a2j+da2j_deta*a2i))/(2.d0*SQRT(a2ijSq))*(1.d0-bipKij)/bVolMix-SQRT(a2ijSq)*(-dbipKij_deta/bVolMix)
				da2Mix_deta=da2Mix_deta+xFrac(iComp)*xFrac(jComp)*da2ij_deta
				!if(LOUD)print*,'MixRule: i,j,a1ij,a2ij = ',iComp,jComp,a1ij,a2ij
				z2ij=(1-bipKij)*( z2j*SQRT(a2i/a2j)+z2i*SQRT(a2j/a2i) )/2  !jre072304
				dz2=(-dbipKij_deta/2.d0)*( z2j*SQRT(a2i/a2j)+z2i*SQRT(a2j/a2i) )+((1.d0-bipKij)/2.d0)*(dz2j_deta*SQRT(a2i/a2j)+z2j*(da2i_deta*a2j-da2j_deta*a2i)/(a2j**2)/(2.d0*SQRT(a2i/a2j))+dz2i_deta*SQRT(a2j/a2i)+z2i*(da2j_deta*a2i-da2i_deta*a2j)/(a2i**2)/(2.d0*SQRT(a2j/a2i)))
				z2ij=z2ij*SQRT( bVoli*bVolj )/bVolMix  
				dz2ij_deta=dz2*SQRT( bVoli*bVolj )/bVolMix
				dz2Mix_deta=dz2Mix_deta+xFrac(iComp)*xFrac(jComp)*(dz2ij_deta+d2KijDeta2*a2ij+dKijDeta*da2ij_deta)
				a2KijMix=a2KijMix+xFrac(iComp)*xFrac(jComp)*a2ij*dKijDeta 	!this is a new term b/c kij=kij(eta)
				daMixQuadPart_deta(iComp)=daMixQuadPart_deta(iComp)+xFrac(jComp)*(daRefij_deta+(da1ij_deta+da2ij_deta/tKelvin)/tKelvin)
			enddo  !jComp
		enddo !iComp

		DO I=1,nComps
			dz0Mix_dN(I)=dz0Mix_deta*detaRef_dN(I)
			da0Mix_dN(I)=da0Mix_deta*detaRef_dN(I)
			dz1Mix_dN(I)=dz1Mix_deta*deta_dN(I)
			da1Mix_dN(I)=da1Mix_deta*deta_dN(I)
			dz2Mix_dN(I)=dz2Mix_deta*deta_dN(I)
			da2Mix_dN(I)=da2Mix_deta*deta_dN(I)
			DO m=1,nComps
				daMixQuadPart_dN(I,m)=daMixQuadPart_deta(I)*deta_dN(m)  !I think this is slightly off for He and H2 mixes. JRE 20200406
			ENDDO
		ENDDO
	endif

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine AxsCalc(xFrac,eta,nComps,tKelvin,a0Mix,iErr) 
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*77 errMsg(11)
	!doublePrecision KIJ,KTIJ
	DIMENSION xFrac(NMX),aMixQuadPart(NMX),bRef(NMX)
	DIMENSION solParEntro(NMX),solParEnerg(NMX) !,zRefMix(5), !a1Mix(5),a2Mix(5),
	!dimension cRef(NMX,3),cRefMix(3)
	dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	common/ksvall/ks0ij,ks1ij
	DIMENSION KII(NMX,NMX),KTII(NMX,NMX)
	COMMON/SIPs/KII,KTII
	!common/XsProps/a0XssS
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
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
		call BeepMsg(errMsg(iErr))
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
	subroutine FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErr) 
	USE Assoc !GlobConst+XA,XD,XC
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs !NOTE: USE Assoc is not here because no assoc is computed.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	PARAMETER(etaStep=2.d0/1000.d0,maxRhoStep=1.d0/etaStep)
	character*77 errMsg(22)
	doublePrecision moleStep,mShape(NMX)
	integer assocFlag
	DIMENSION xFrac(NMX),gmol_old(NMX),aMixQuadPart(NMX),gmol(NMX),bRef(NMX)

	DIMENSION dh_dN(NMX),dFUGASSOC_dT(NMX),dFUGASSOC_dN(NMX,NMX),dFUGASSOC_dRHO(NMX),FUGASSOC(NMX)	  
	DIMENSION dh_dN_num(NMX),dFUGASSOC_dN_num(NMX,NMX),fugAssocLoop(NMX)
	DIMENSION dFUGREP_dN(NMX,NMX),dFUGATT_dN(NMX,NMX),dFugCof_dN(NMX,NMX),dFUGQUAD_dN(NMX,NMX),dzFactor_dN(NMX) 
	DIMENSION dzAtt_dN(NMX),daAtt_dN(NMX),dzAssoc_dN_num(NMX),dbVolMix_dN(NMX),dbRefMix_dN(NMX) 

	COMMON/z_aComps/zRefij,aRefij,dKijDeta,z1ij,a1ij,z2ij,a2ij
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
	!common/Derv1RG/dPRG_dV,d2PRG_dV2
	COMMON/DervMix/dz0Mix_dN(NMX),da0Mix_dN(NMX),dz1Mix_dN(NMX),da1Mix_dN(NMX),dz2Mix_dN(NMX),da2Mix_dN(NMX),daMixQuadPart_dN(NMX,NMX),daMixQuadPart_deta(NMX)
	COMMON/DervMix2/dz0Mix_deta,da0Mix_deta,dz1Mix_deta,da1Mix_deta,dz2Mix_deta,da2Mix_deta,d2z0Mix_deta2,d2z1Mix_deta2,d2z2Mix_deta2
	COMMON/DervMix3/d2a1Mix_deta2,d3a1Mix_deta3,d2a2Mix_deta2
	COMMON/DervZMix3/d3z0Mix_deta3,d3z1Mix_deta3,d3z2Mix_deta3
	COMMON/ETA2/ETA
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	COMMON/dfug2/dFUGASSOC_dN_num,dh_dN_num
	COMMON/dFug/h_nMichelsen,dh_dT,dh_dN,dFugassoc_dT,dFugassoc_dN,dFUGASSOC_dRHO,dh_dRHO
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	COMMON/A/a0Mix,a1Mix,a2Mix,aAtt,a_d,zRep
	common/aRG/zNC,aDepNC
	common/rg/iFlagRG,iFlagFit,mShape,cutOff(NMX),Phi(NMX),Zi(NMX)
	data initCall/1/

	iErr=0
	errMsg(5)=' FuTptVtot Warning: T < tKmin for at least one component.'
	errMsg(11)=' FuTptVtot Error: Nonsense input. e.g. eta > 1'
	!errMsg(4)=' FuTptVtot Error: a1i or a2i > 0'
	errMsg(12)=' FuTptVtot Error: Wertheim returned error.'
	errMsg(13)=' FuTptVtot Error: Error from MixRule. Terminal, sorry.'
	if( (1-eta) < 1.e-11)then
		if(LOUD)write(*,*) 'Error in FuTptVtot: eta > 1'
		if(LOUD)pause
        iErr=11
        return
	endif

	totMoles=SUM(gmol)
	bVolMix=0
	bRefMix=0
    isTtooLow=0
	TcVolatile=12345
    DO i=1,nComps
		if(Tc(i) < TcVolatile)then
			TcVolatile=Tc(i)
			iVolatile=i
		endif
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
        !chemPo(iComp)=86 !initialize to error indicator in case of bad return !No chemPo here !!!
    enddo
    if(tKelvin < tKmin(iVolatile)*0.998D0)then	   !round TKmin down a little for possible lost precision
        iErr=5
        if(LOUD.and.initCall)write(*,'(a,i3,2f7.1,a)')' FuVtot: iVolatile,T(K),TKmin= ',iVolatile,tKelvin,tKmin(iVolatile),TRIM( errMsg(iErr) )
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
	!call XsAs(nComps,iErrXs)  !for debugging
	!print*,'FuTptVtot: calling MixRule'
	call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
	if(iErrMix > 10)then
		iErr=13
		call BeepMsg(errMsg(iErr))
		return
	endif 
	Gf=1
	dGf_dEta=0
	d2Gf_dEta2=0
	uDepGf=0
	!if(nComps==1)call GaussFun3(tKelvin,eta,Gf,dGf_dEta,d2Gf_dEta2,cpGf,uDepGf,iErr)
	!denom=(1.d0-eta)*(1.d0-eta)*(1.d0-eta)
	zRep=z0Mix
	aTvRef=a0Mix
	!Adep=A0+A1/T+A2/T^2*Gf=> Z-1 = Z0 + Z1/T + (eta*dA2/dEta*Gf + A2*eta*dGf_dEta)/T^2
	!& dZ_dEta = dZ_dEta0 + dZ_dEta1/T + (dA2/dEta*Gf + A2*dGf_dEta)/T^2+ (eta*d2A2/dEta2*Gf +2*eta*dA2/dEta*dGf_dEta + A2*eta*d2Gf_dEta2)/T^2
	!P = Z*rho*R*T => dP_dRho = Z*R*T + rho*R*T*dZ_dRho = R*T*(Z+eta*dZ_dEta)
	!dRho/rho = -dV/V => V*dP_dV = -rho*dP_dRho => dP_dV = -rho^2*dP_dRho  
	z2Gex=z2Mix*Gf+a2Mix*eta*dGf_dEta
	aAtt=( a1Mix+ a2Mix*Gf/tKelvin)/tKelvin 	 ! => a3 = a2*g0*exp{-g1*[(eta-eta0)^2-eta0^2]}/3; g0[=]K
	zAtt=( z1Mix+ z2Gex/tKelvin)/tKelvin	 ! 
	uAtt=( a1Mix+ 2.d0*a2Mix/tKelvin*Gf+a2Mix/tKelvin*uDepGf )/tKelvin 
	!a_d=a2Mix/a1Mix/TC(1)
	
	!C We define a flag for association term here to avoid calling Werthiem routine when it is not necessary...
	!print*,'FuTptVtot: calling Wertheim'
	if(iErr < 10)call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode)
!	if( ABS(zAssoc) < 1e-11 )print*,'zAssoc ~ 0???'		 !for debugging
	IF(iErrCode)then
		!iErr=iErrCode
		if(LOUD)pause 'FuTptVtot: Wertheim gave error.'
		iErr=12
		return
	endif
	if ( ABS(zAssoc)>1D-33 ) then
		if(isZiter==0)call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
		assocFlag=1
	else
		assocFlag=0
	endif

	if(iErr > 10)then
		if(LOUD)write(*,*)errMsg(iErr)
		return
	endif
! RG method disabled. Prefer GaussEx method. JRE 20191010 
!	if (iFlagRG.EQ.1) then
!		call HelmRG(nComps,gmol,tKelvin,eta,bVolMix,aDepNC,zNC,PMpaNC,dP_detaNC,d2P_deta2NC,d3P_deta3NC)
!	else 																			
		aDepNC=0.d0
		zNC=0.d0
		PMpaNC=0.d0
		dP_detaNC=0.d0
		d2P_deta2NC=0.d0
		d3P_deta3NC=0.d0
!	endif
	
	aDep=aTvRef+aAtt+aAssoc+aDepNC
	uDep=uAtt+uAssoc
	!if(LOUD)print*,'uAtt,uAssoc',uAtt,uAssoc,uDep
	zFactor=(1.d0+zRep+zAtt+zAssoc)+zNC
	PMpa=(zFactor)*rGas*rhoMol_cc*tKelvin 
	DAONKT=aDep
	DUONKT=uDep
	dSoNk=-(aDep-uDep) 
	DHONKT=uDep+zFactor-1.d0
	
	aRes_RT=aDep
	uRes_RT=uDep
	sRes_R =-(aDep-uDep) 
	hRes_RT=uDep+zFactor-1.d0
	if(isZiter==0)then
		if(LOUD.and.initCall)print*,'FuTptVtot: T,rho,z=',tKelvin,rhoMol_cc,zFactor
		call NUMDERVS(nComps,xFrac,tKelvin,rhoMol_cc,zFactor,cmprsblty,iErrNum)
		if(LOUD.and.initCall)print*,'FuTptVtot: cmprsblty,CvRes_R=',cmprsblty,CvRes_R
		if(iErrNum > 1.and.LOUD)then
			print*,'FuTptVtot: NUMDERVS failed. No derivatives for you!'
			iErr=16
		endif
		return	!all or nothing on numdervs vs analytical.
	endif
	!all below is obsoleted by NUMDERVS. JRE 20200425
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Added by AFG 2010
!C Derivatives of ref. and perturbation terms in respect to eta 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	dzRep_deta=dz0Mix_deta
	d2zRep_deta2=d2z0Mix_deta2
	d3zRep_deta3=d3z0Mix_deta3
	aTvRef=a0Mix

	!aAtt=( a1Mix+ a2Mix*Gf/tKelvin)/tKelvin 	 !
	!zAtt=( z1Mix+ z2Gex/tKelvin)/tKelvin	 
	!z2Gex=z2Mix*Gf+a2Mix*eta*dGf_dEta
	!=>dZ2Gex_dEta=dZ2Mix_dEta*Gf+2*z2Mix*dGf_dEta+a2Mix*dGf_dEta+a2Mix*eta*d2Gf_dEta2
	dz2Gex_dEta=dZ2Mix_dEta*Gf+2*z2Mix*dGf_dEta+a2Mix*dGf_dEta+a2Mix*eta*d2Gf_dEta2
	dzAtt_deta=(dz1Mix_deta+dz2Gex_deta/tKelvin)/tKelvin	 !
	d2zAtt_deta2=(d2z1Mix_deta2+d2z2Mix_deta2/tKelvin)/tKelvin 	!WARNING: d2z and d3z do not include Gex
	d3zAtt_deta3=(d3z1Mix_deta3+d3z2Mix_deta3/tKelvin)/tKelvin	!WARNING: d2z and d3z do not include Gex
	deta_dV=-eta/vTotCc
	d2eta_dV2=(-deta_dV*vTotCc+eta)/(vTotCc*vTotCc)
	d3eta_dV3=-3.d0*(d2eta_dV2)/vTotCc
	const1=rGas*tKelvin/bVolMix
	
	dzAssoc_deta=-0.5d0*(dh_drho*(1.d0/bVolMix)*(1.d0+dLng)+dAlph_deta*h_nMichelsen)
	d2zAssoc_deta2=0.d0
	d3zAssoc_deta3=0.d0
	dzFactor_deta=dzRep_deta+dzAtt_deta+dzAssoc_deta 
	d2zFactor_deta2=d2zRep_deta2+d2zAtt_deta2+d2zAssoc_deta2 
	d3zFactor_deta3=d3zRep_deta3+d3zAtt_deta3+d3zAssoc_deta3
	dzFactor_dV=dzFactor_deta*deta_dV
	d2zFactor_dV2=d2zFactor_deta2*deta_dV**2+dzFactor_deta*d2eta_dV2
	d3zFactor_dV3=d3zFactor_deta3*deta_dV**3+d2zFactor_deta2*2.d0*deta_dV*d2eta_dV2
	d3zFactor_dV3=d3zFactor_dV3+d2zFactor_deta2*deta_dV*d2eta_dV2+d3eta_dV3*dzFactor_deta
	cmprsblty=zFactor+eta*dzFactor_deta
	!if(LOUD)print*,'cmprsblty=',cmprsblty
	! First derivatives of P,ZREP,ZATT and ZASSOC in respect to T according to the ESD EOS while V and N are constant
	betaDZREP_dBeta=0.d0
	!zAtt=eta*dA_dEta=Z1/T+Z2/T^2> betaDZAtt_dBeta = Z1/T + 2Z2/T^2 
	!zAtt=( z1Mix+ z2Gex/tKelvin)/tKelvin	 ! 
	betaDZatt_dBeta= (z1Mix+2*z2Mix/tKelvin)/tKelvin !omitting Gex contribution for now. JRE 20191004 
	dZASSOC_dT=-half*(dLng+1)*dh_dT !dh_dT is computed in WertheimFugc.
	betaDZassoc_dBeta= -tKelvin*dZASSOC_dT
	betaDZ_dBeta= betaDZrep_dBeta+betaDZatt_dBeta+betaDZassoc_dBeta 

	dP_deta=const1*(dzFactor_deta*eta+(zFactor-zNC))+dP_detaNC 
	d2P_deta2=const1*(d2zFactor_deta2*eta+2.d0*dzFactor_deta)+d2P_deta2NC 
	d3P_deta3=const1*(d3zFactor_deta3*eta+3.d0*d2zFactor_deta2)+d3P_deta3NC  

	dP_dV=deta_dV*dP_deta 
	d2P_dV2=d2P_deta2*deta_dV**2+dP_deta*d2eta_dV2 
	d3P_dV3=d3P_deta3*deta_dV**3+d2P_deta2*2.d0*deta_dV*d2eta_dV2
	d3P_dV3=d3P_dV3+d2P_deta2*d2eta_dV2*deta_dV+d3eta_dV3*dP_deta 
	!print*,'FuTptVtot: done with bulk derivatives. Starting numerical assoc derivatives.'
    
	if (isZiter==0) then !Determine assoc derivatives by numerical differentiation
	    P_old=PMpa
	    Z_old=zFactor
	    a_old=aDep
	    u_old=uDep
		if (assocFlag==1) then
			two=2.d0
			!half=1.d0/two	!half is defined in GlobConst
			moleStep=10000.d0
			rhoMol_cc_old=rhoMol_cc
			totMoles_old=totMoles
			bVolMix_old=bVolMix
			eta_old=eta
			DO I=1,nComps  !compute dFugAssoc(j)/dNi and similar numerically. For CritMix().
				if(gmol(i)< 1d-11)gmol(i)=1d-11 !avoid zero divide error for inf dilution calculations
				gmol_old(I)=gmol(I)
				gmol(I)=gmol(I)+gmol(I)/moleStep !do we need to store (calling argument) gmol before altering it? JRE 20200406
				totMoles=SUM(gmol)
                xFrac=gmol/totMoles
				bVolMix=SumProduct(nComps,xFrac,bVolCC_mol)
				rhoMol_cc=totMoles/vTotCc
				eta=rhoMol_cc*bVolMix
				call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
				call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode)
				call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
				DO J=1,nComps
					fugAssocLoop(J)=FUGASSOC(J)
				ENDDO
				zAssocLoop=zAssoc
				gmol(I)=gmol_old(I)-gmol_old(I)/moleStep
				totMoles=SUM(gmol)
                xFrac=gmol/totMoles
				bVolMix=SumProduct(nComps,xFrac,bVolCC_mol)
				bRefMix=SumProduct(nComps,xFrac,bRef)
				rhoMol_cc=totMoles/vTotCc
				eta=rhoMol_cc*bVolMix
				etaRef=rhoMol_cc*bRefMix
				call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
				call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode)
				call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
				gmol(I)=gmol_old(I)
				totmoles=totmoles_old
				dzAssoc_dN_num(I)= (zAssocLoop-zAssoc)/(two*gmol(I)/moleStep)
				DO J=1,nComps
					dFUGASSOC_dN_num(J,I)=((fugAssocLoop(J)-FUGASSOC(J))/(two*gmol(I)/moleStep))*totmoles
				ENDDO
 			ENDDO
			rhoMol_cc=rhoMol_cc_old
			eta=eta_old
			totmoles=totmoles_old
			bVolMix=bVolMix_old
			PMpa=P_old
			zFactor=Z_old
			do ii=1,nComps
				xFrac(ii)=gmol(ii)/totMoles
			enddo
			call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
			call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode)
			call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
		endif
		DO I=1,nComps
			!First derivative in respect to Ni at constant T and V
			dbVolMix_dN(I)=bVolCC_mol(I)-bVolMix
			dbRefMix_dN(I)=bRef(I)-bVolMix
			dzAtt_dN(I)=( dz1Mix_dN(I)+ dz2Mix_dN(I)/tKelvin )/tKelvin
			daAtt_dN(I)=( da1Mix_dN(I)+ da2Mix_dN(I)/tKelvin )/tKelvin
			dzFactor_dN(I)=dz0Mix_dN(I)+dzAtt_dN(I)+dzAssoc_dN_num(I)
			dP_dN(I)=rhoMol_cc*rGas*tKelvin*(zFactor+dzFactor_dN(I))
		ENDDO
		DO kComp=1,nComps
			!First derivative in respect to Ni at constant T and V
			DO I=1,nComps	
				dFUGREP_dN(kComp,I)=bRef(kComp)/bRefMix*(dz0Mix_dN(I)-da0Mix_dN(I))-bRef(kComp)*(z0Mix-a0Mix)*dbRefMix_dN(I)/(bRefMix*bRefMix)
				dFUGATT_dN(kComp,I)=bVolCC_mol(kComp)/bVolMix*(dzAtt_dN(I)-daAtt_dN(I))-bVolCC_mol(kComp)*(zAtt - aAtt)*dbVolMix_dN(I)/(bVolMix*bVolMix)
				dFUGQUAD_dN(kComp,I)=2.d0*daMixQuadPart_dN(kComp,I)
				dFugCof_dN(kComp,I)=dFUGREP_dN(kComp,I)+dFUGATT_dN(kComp,I)+dFUGQUAD_dN(kComp,I)+dFugAssoc_dN_num(kComp,I)-dzFactor_dN(I)/zFactor
				IF (kComp.eq.I) THEN
					!Which is actually dLnf(i)/dN(k)
					dFUG_dN(kComp,I)=(dFugCof_dN(kComp,I)+dP_dN(I)/P_old+(1.d0-xFrac(kComp))/(xFrac(kComp)+1d-33))! +1d-33 to avoid zero divide at infinite dilution.
				ELSE
					dFUG_dN(kComp,I)=(dFugCof_dN(kComp,I)+dP_dN(I)/P_old-1.d0)!/totMoles
				ENDIF
			ENDDO
		ENDDO
	endif !isZiter===0
	
!	!print*,'FuTptVtot: Derivatives done. Getting pure comp derivative props.'
!	cvRes_R=0
!	if(nComps==1)then	!added by JRE 2018. I'm not confident of formulas for NC>1 as of 20191004.
!		!uAssoc,CvAssoc computed in Wertheim() and passed in module Assoc.
!		!if(LOUD)print*,'uAssoc,CvAssoc=',uAssoc,CvAssoc
!		!cvAtt = uAtt - beta*d(Ures_RT)/dBeta
!		!uAtt = A1/T+2*A2/T^2 => beta*dUAtt_dBeta= A1/T+4*A2/T^2
!!		cvAtt = uAtt - (a1Mix+4*a2Mix/tKelvin)/tKelvin !=4*A2/T^2
!		!if(LOUD)print*,'uAtt,uAssoc',uAtt,uAssoc,uRes_RT
!		cvRes_R = CvAtt + CvAssoc
!		!if(LOUD)print*,'FuTptVtot: uAtt,CvRes_R',uAtt,cvRes_R
!		!if(LOUD)print*,'FuTptVtot: TdZ_dT',-betaDZ_dBeta
!		TdP_dT=zFactor-betaDZ_dBeta
!		!if(LOUD)print*,'dZASSOC_dT,dh_dT ',dZASSOC_dT,dh_dT
!		if(ABS(cmprsblty) < 1e-11 .and. LOUD)pause 'FuTptVtot: cmprsblty ~ 0'
!		CpRes_R = cvRes_R-1+TdP_dT**2/cmprsblty 
!	endif
 
	initCall=0
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine XsAs(nComps,iErr)
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	character*77 errMsg(11)
	!DIMENSION XA(NMX,maxTypes),XD(NMX,maxTypes),XC(NMX,maxTypes)
	dimension a0Pure(2),a1Pure(2),a2Pure(2),xFrac(NMX),aMixQuadPart(NMX)
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine ChemPoCalcTpt(chemPo,eta,tKelvin,gmol,nComps,iErr)
	!need bVolMix to get rho from eta because alpha(i) ~ rho*bondVol(i) ~ eta*bondVol(i)/bVolMix
	!need vMolecNm3 in Wertheim	just because Marty wrote the code badly.  
	USE Assoc !GlobConst+XA,XD,XC
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DoublePrecision chemPoAssoc(NMX),chemPo(NMX),xFrac(NMX),gmol(NMX)!,aAttk(NMX)
	!DIMENSION bVolCC_mol(NMX)!,xFrac_old(NMX),dXA_dx(NMX)
	DIMENSION aMixQuadPart(NMX),bRef(NMX) !,zRefMix(5)!,a1Mix(5),a2Mix(5)
	DIMENSION dFUGASSOC_dT(NMX),dFUGASSOC_dRHO(NMX)!,FUGASSOC(NMX),dh_dN(NMX),dFUGASSOC_dN(NMX,NMX)

	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	common/aRG/zNC,aDepNC

	iErr=0
	!iErr=10 => SUMX .ne. 1
	!iErr=11 => -ve zFactor calculated
	totMoles=SUM(gmol) ! SUM() is an intrinsic function
	if(ABS(totMoles) < 0.001 .and. LOUD)pause 'FuTpt: totMoles ~0? Will try renormalizing, but you should check input.'
	bVolMix=0
	bRefMix=0
    DO i=1,nComps
		xFrac(i)=gmol(i)/totMoles
        bVolMix=bVolMix+xFrac(i)*bVolCC_mol(i)
		bRef(I)=bVolRef(i,tKelvin) !Function bVolRef()
        bRefMix=bRefMix+xFrac(i)*bRef(i)
	ENDDO

	call MixRule(isZiter,xFrac,tKelvin,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,z0Mix,z1Mix,z2Mix,aMixQuadPart,iErrMix)
	if(iErrMix > 10 .and. LOUD)pause 'ChemPoCalctpt: fatal error from MixRule.'
	aAtt=(a1Mix+a2Mix/tKelvin)/tKelvin  
	zAtt=(z1Mix+z2Mix/tKelvin)/tKelvin
	Gf=1
	dGf_dEta=0
	if(nComps==1)call GaussFun3(tKelvin,eta,Gf,dGf_dEta,d2Gf_dEta2,cpGf,uDepGf,iErr)
	aAtt=( a1Mix+ a2Mix/tKelvin*Gf)/tKelvin 	 !
	zAtt=( z1Mix+ z2Mix/tKelvin*Gf+a2Mix/tKelvin*eta*dGf_dEta)/tKelvin	 ! 
	uAtt=(a1Mix+ 2.d0*a2Mix/tKelvin*Gf+a2Mix/tKelvin*uDepGf )/tKelvin 
	rhoMol_Cc=eta/bVolMix
	call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode) !NDS,NAS,ND,
	IF(iErrCode.ne.0)then
		iErr=iErrCode
		if(LOUD)pause 'Error in FuTpt: Wertheim gave error.'
		return
	endif
! The last term is the nonClassical part of Z
	zFactor=(1+z0Mix+zAtt+zAssoc)+zNC
	
	if(zFactor.lt.1e-33)then
		if(LOUD)write(*,*)'zFactor=',zFactor 
		if(LOUD)pause 'ChemPoCalcTpt: zFactor < 1E-55. LOG(Z)?'
		iErr=11
		return
	endif
	call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,eta,chemPoAssoc,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)! ,NDS,NAS,ND
	DO kComp=1,nComps
		bVolk=bVolCC_mol(kComp) 
		bRefk=bVolRef(kComp,tKelvin) !Function bVolRef()...
		FUGREP=bRefk/bRefMix*(z0Mix-a0Mix)
		FUGATT=bVolk/bVolMix*(zAtt - aAtt)
		!aAttQuadPart=aAttQuadPart+etaPow*aMixQuadPart(kComp,iCoeff)
		FUGQUAD=2.d0*aMixQuadPart(kComp) !we lump all the quadparts together in the current version.  
		FUGBON=chemPoAssoc(kComp)
		
! The last two terms comprise the nonClassical part of Fugacity
		chemPo(kComp)=FUGREP+FUGATT+FUGQUAD+FUGBON-DLOG(zFactor) !+(zNC-aDepNC)*bVolk/bVolMix+2.d0*aDepNC 
		!print*,'fugrep,fugatt,fugquad',fugrep,fugatt,fugquad
	enddo
	!print*,'a0mix,aAtt,z0mix,zAtt',a0mix,aAtt,z0mix,zAtt
	!print*,'gRes,ChemPo',(a0mix+aAtt+z0mix+zAtt),(fugrep+fugAtt+fugQuad)
	!pause 'check parts.'
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DoublePrecision Function bVolRef(i,tKelvin)
	USE GlobConst !id()
    Implicit DoublePrecision(A-H,O-Z)
	!Purpose: Enable softening of the reference contribution using the E&D FPE'86 correlation.
	DoublePrecision nEff ! simplifies the equations while accurately conveying the name.
    dEhs_sigmaCube=1
    if(id(i)==913 .or. id(i)==902)then  !id(i) is part of GlobConst
        nEff=9
		if(id(i)==913)nEff = 15
		p2 = (.0093d0*nEff-.0592d0)*nEff !E&D(1986)
		p2 = ( (nEff-6)/6 )**2 /exp(  gammln( (nEff-3)/nEff )  )**(2*nEff/3) ! ESK(2019)
		!print*,'bVolRef: p2 = ', p2
		!pause
		p1 = (nEff-1)  
		p1 = 1.13*nEff-3.28 !ESK(2019) 
		eps_kB=2	!Helium
		if(id(i)==902)eps_kB=15
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
	USE GlobConst !RGAS,...NMX, ...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!DIMENSION XA(NMX),XD(NMX),chemPoAssoc(NMX)!,aAttk(NMX)
	dimension xFracPert(NMX),chemPo(NMX),xFrac(NMX),gmol(NMX),bRef(NMX)
	DIMENSION aMixQuadPart(NMX) !,zRefMix(5)!,a1Mix(5),a2Mix(5)
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	common/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	common/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT

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

