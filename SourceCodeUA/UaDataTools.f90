
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine PGLStartup(NC,iEosLocal,idOpt,ierCode)
	!Purpose: CALL EOS Get___ ROUTINES, also includes LoadCritParmsDb(if needed),GetCrit, Get(EOS), GetBips, BipIo, 
	!	SETS UP THE CRITS, EOSPARMS, bips, DEBUG status, FUGI(iEosOpt)
	!	Echoes user IO to Output.txt, and reports error checks. 
	!	INITIALIZATION AND CALLING SEQUENCE FOR VLE, LLE, VLLE SUBROUTINES.
	!reqd routines:
	!	bubpl.for, bubtl.for, dewtv.for, flashsub.for, FuEsdMy.for FuEsdXs2.for, FugiPr.f90, FugiPrws.for, RegPure.f90
	!	LmDifEzCov2.for, Mintools.for
	!               1     2       3       4          5          6         7           8              9            10          11
	USE GlobConst
	USE CritParmsDb
	USE BIPs
	USE EsdParms
	Implicit DoublePrecision(A-H,K,O-Z)
	CHARACTER*77 errMsgPas !,readString,Property(22)
	!CHARACTER*251 dumpFile
    Integer  localCas(NMX),localID(NMX),iEosLocal,ierCode,NC 
	if(idOpt==2)then
		localCas(1:NC)=idCas(1:NC)	! idCas USEd from GlobConst
	elseif(idOpt==1)then
		localID(1:NC)=ID(1:NC)	! idCas USEd from GlobConst
		call IdCasLookup(NC,localID,iErrLook,errMsgPas) ! idCas USEd in GlobConst to set.
		if(iErrLook > 0)then
			ierCode=11
			if(LOUD)write(dumpUnit,*)'PGLStartup: from idDipprLookup-'//TRIM(errMsgPas)
			return
		endif
		localCas(1:NC)=idCas(1:NC)	! idCas USEd from GlobConst
	else
		ierCode=12
		return
	endif
	Call PGLWrapperStartup(NC,iEosLocal,localCas,ierCode)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine PGLWrapperStartup(NC,iEosLocal,localCas,ierCode)
	!Purpose: CALL EOS Get___ ROUTINES, also includes LoadCritParmsDb(if needed),GetCrit, Get(EOS), GetBips, BipIo, 
	!	SETS UP THE CRITS, EOSPARMS, bips, DEBUG status, FUGI(iEosOpt)
	!	Echoes user IO to Output.txt, and reports error checks. 
	!	INITIALIZATION AND CALLING SEQUENCE FOR VLE, LLE, VLLE SUBROUTINES.
	!reqd routines:
	!	bubpl.for, bubtl.for, dewtv.for, flashsub.for, FuEsdMy.for FuEsdXs2.for, FugiPr.f90, FugiPrws.for, RegPure.f90
	!	LmDifEzCov2.for, Mintools.for
	!               1     2       3       4          5          6         7           8              9            10          11
	USE GlobConst
	USE CritParmsDb
	USE BIPs
	USE EsdParms
	Implicit DoublePrecision(A-H,K,O-Z)
	CHARACTER*77 errMsgPas !,readString,Property(22)
	!CHARACTER*251 dumpFile
    Integer  localCas(NMX),iEosLocal,ierCode,NC 
	
	idCas(1:NC)=localCas(1:NC)
    if(LOUD)write(dumpUnit,*)'PGLWrapperStartup: starting.NC,iEosOpt,idCas=',NC,iEosOpt ,idCas(1:NC)

	INITIAL=0
	iEosOpt=iEosLocal
	ierCode=0 ! Initialize to failure in the iEosOpt write slot to make it easier to terminate if any Get_ functions fail.
	iErrCrit=0
	if(nDeckDb.ne.nCritSet)call LoadCritParmsDb(iErrCrit) !nDeckDb.ne.nCritSet indicates that LoadCrit has not been called.
	if(iErrCrit > 0)then
		if(LOUD)write(dumpUnit,*)'PGLWRapperStartup: Error from LoadCritParmsDb, iErrCrit=',iErrCrit
		ierCode=11
		goto 86
	endif
	call IdDipprLookup(NC,localCas,iErrCas,errMsgPas) ! ID USEd in GlobConst to set.
	if(LOUD)write(dumpUnit,*)'PGLWRapperStartup:localCas,idDippr=',(localCas(i),id(i), i=1,NC)
	if(iErrCas/=0)then
        if(LOUD)write(dumpUnit,*)' Sorry, must abort.  Not found for CAS number(s)=', (localCas(i),i=1,NC)
        write(dumpUnit,*)ierCode,' PGLWRapperStartup:Sorry, must abort.  Not found for CAS number(s)=', (localCas(i),i=1,NC)
		ierCode=12
        goto 86
    endif
	if(LOUD)write(dumpUnit,*)'idDippr,class=',ID(1),TRIM(class(1))
	CALL GETCRIT(NC,iErrCrit)
	if(iErrCrit/=0)then
		if(LOUD)write(dumpUnit,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		write(dumpUnit,*)ierCode,' Error in Main: ErrorCode from GetCrit = ',iErrCrit
		ierCode=13
		goto 86
	endif
	if(LOUD)write(dumpUnit,*)'calling GetEOS'
	iErrGet=0
	if(iEosOpt==1)CALL GetPR(NC,iErrGet)
	if(iEosOpt==2)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USEd EsdParms					!Diky model 12
	if(iEosOpt==3)CALL GetPRWS(NC,iErrGet) 
	if(iEosOpt==4)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USEd EsdParms					!Diky model 12
	if(iEosOpt==5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)   !Results placed in USEd: SpeadParms, Assoc			!Diky model 6
	if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
	if(iEosOpt==8)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in USEd: SpeadParms, Assoc
	if(iEosOpt==9)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in USEd: SpeadParms, Assoc
	if(iEosOpt==10)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Gross's PcSaft parameters		!Diky model 25
	if(iEosOpt==11)CALL GetPrTc(NC,iErrGet)		!JRE 2019 : Reads Jaubert's parameters						!Diky model 22
	if(iEosOpt==12)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParms, Assoc				!Diky model 23
	if(iEosOpt==13)CALL GetEsdCas(NC,localCas,iErrGet)	 !Results placed in USE EsdParms, Assoc				!Diky model 24
	if(iEosOpt==14)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!Results placed in USEd: SpeadParms, Assoc			!Diky model 18
	if(iEosOpt==15)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Emami's PcSaft parameters		!Diky model 26
	if(iEosOpt==16)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Emami's(Tb) PcSaft parameters	!Diky model 27
	if(iEosOpt==20)CALL GetPcSaft(NC,localCas,iErrGet)		!JRE 2019 : Reads Emami's(Tb) PcSaft parameters	!Diky model 27
	!if(iEosOpt.eq.17)This is the model number for COSMOtherm, if one is required JRE 20200119.				!Diky model ??
	!NewEos: Add here for initializing parms.
	if(iErrGet .ne. 0)then
		if(LOUD)write(dumpUnit,*)ierCode,' Error in Main: failed to get required dbase props.'
		if(LOUD)write(dumpUnit,*)iErrGet,' PGLWrapperStartup:Get(EOS) failed. iErrGet=',iErrGet
		ierCode=iErrGet*10
		goto 86
	endif
    if(LOUD)write(dumpUnit,*)'PGLWrapperStartup: Success so far... Got EOS Parms for iEosOpt=',iEosOpt
86	return
	end	   !PGLWrapperStartup



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine ReadNext(inUnit,nexString3,maxLines,maxVars,bHeader,header77,nVars,var,nLines,iErr)
	! Purpose: Read uncounted data with "next" character string between data sets.
	! Notes:   inFile for inUnit must be opened before calling.
	!	       var(maxLines,maxVars) must be declared in calling routine with exactly these dimensions.
	!		   datasets may be read with or without header per dataset.
	! 
	Implicit DoublePrecision(a-h,o-z)
	DoublePrecision	var(maxLines,maxVars)
	Character*77 header77,dumString !,inFile77
	Character*3 nexString3,nexCheck
	LOGICAL bHeader
	ierRead=0
	lineCount=0
	do while(ierRead==0)
		READ(inUnit,'(a77)',ioStat=ioErr)dumString
		if(ioErr == -1)then ! end of file.
			iErr= -1
			if(lineCount==0)iErr= -2
			return
		endif
		read(dumString,'(a3)')nexCheck
		if( TRIM(nexCheck)==TRIM(nexString3) )then
			nLines=lineCount-1
			return
			!cycle	! don't cycle here. The data for this system should be loaded now. So do regression, then cycle.
 		elseif(lineCount==0.and.bHeader)then
			header77=TRIM(dumString)
			lineCount=1
		else
			read(dumString,*,ioStat=ioErr)(var(lineCount,iVar),iVar=1,nVars)
			lineCount=lineCount+1
		endif
	enddo ! cycle back to get more data. goto 101 when "***Next***"
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine GetTmHfus(id,Tm,Hfus,iErr)
	USE GlobConst, Only:PGLinputDir,LOUD,dumpUnit 
	! Purpose: Read Tm(K) and Hfus(J/mol) from input file.
	Implicit DoublePrecision(a-h,o-z)
	Character*77 HfusFile,FN !,inFile77
    FN="Hfus.txt"
	HfusFile=TRIM(PGLinputDir)//'\'//TRIM(FN)
	if(LOUD)write(dumpUnit,*)'HfusFile=',TRIM(HfusFile)
	inUnitHfus=51
	open(inUnitHfus,file=HfusFile)
	read(inUnitHfus,*)nHfus
	notFound=1
	iErr=0
	Tm=0   ! wipe out any previously stored values to avoid confusion if error occurs.
	Hfus=0
	do i=1,nHfus
		read(inUnitHfus,*,ioStat=ioErr)idCas,idDippr,TmDb,HfusTde,HfusDip
		if(ioErr/=0)cycle
		if(id==idDippr)then
			notFound=0
			Tm=TmDb
			Hfus=HfusTde*1000 ! convert to J/mol so Rgas=8.314 when computing SLE.
			if(idCas==120127.or.idCas==92524)then
				write(dumpUnit,*)'idCas,ID,Hfus=',idCas,idDippr,HfusTde,HfusDip,Hfus
				!write(dumpUnit,*) 'GetTmHfus: check HfusTde,HfusDip,Hfus again'
			endif
			exit ! leave if done.
		endif
	enddo
	if(notFound==1)iErr=11
	if(ioErr /= 0)iErr=12
	close(inUnitHfus)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine BeepMsg(msg)
	USE GlobConst
	character*222 msg
	if(LOUD)write(dumpUnit,*)TRIM(msg) !-- causes crash???
!	write(dumpUnit,*)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdDipprLookup(NC,idCas,iErr,errMsgPas)
    USE GlobConst, only:ID,nmx,LOUD,dumpUnit,nCritSet
	USE CritParmsDb, only:idCasDb,IDnum,CrIndex,nDeckDb
	Implicit NONE
    integer idCas(nmx),iErr,NC,iComp,iGotIt,I
	character*77 errMsg(0:22),errMsgPas
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	errMsg(0)='No Problem'
	errMsg(11)='IdDipprLookup Error: at least one id not found'
	errMsg(12)='IdDipprLookup: You must call LoadCritParmsDb first. Goodbye!'
	iErr=0
	!write(dumpUnit,*)'IdLookup: nDeckDb=',nDeckDb
	if(nDeckDb.ne.nCritSet)then
		iErr=12
		if(LOUDER)write(dumpUnit,*)TRIM(errMsg(iErr))
		goto 861
	endif

	do iComp=1,NC
		iGotIt=0
		do i=1,nDeckDb
			if(idCasDb(i)==idCas(iComp))then
				ID(iComp)=IDnum(i)
				!class(iComp)=classdb(i)	!class() is not available. Call idCasLookup() if you need class().
				iGotit=1
				exit ! breaks the inner loop only?
			endif
		enddo
		if(iGotIt==0)iErr=11 !only way here is if cycle failed => iGotIt==0.
		if(LOUD.and.iErr>0)write(dumpUnit,*)'Did not find idCas(i).i,id=',iComp,id(i)
		if(iErr > 0)goto 861
	enddo
861	errMsgPas=TRIM(errMsg(iErr))
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE LoadCritParmsDb(iErrCode)
	!C  
	!C  PROGRAMMED BY:  JRE 12/21
	!C  Loads THE CRITICAL PROPERTIES (Mw, etc) into GlobConst. 
    !C      Includes CrIndex(idDippr)=line where idDippr was found (linked list)
	!C      This should be a faster way of loading properties, e.g. when running VLE evaluations for a large db.
	!C  INPUT
	!C    idOpt 1 if ID is set and need to lookup idCas, 2 if idCas is set and need to lookup ID.
	!C    (ID	 VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS USEd from GlobConst if idOpt=1)
	!C    (idCas VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS USEd from GlobConst if idOpt=2)
	!C  OUTPUT
	!C    TC		CRITICAL TEMPERATURE(K)
	!C    PC		CRITICAL PRESSURE (MPa)
	!C    ACEN		ACENTRIC FACTOR
	!C    Class		COMPONENT class {norml,polar,gases,assoc,Asso+...}
	!C    CrIndex	e.g., Tc(icomp)=TcD(CrIndex(ID(iComp)))
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet
	USE CritParmsDb 
	IMPLICIT DoublePrecision(A-H,O-Z)
    Integer ndbLocal
	!CHARACTER*4 tCode,pCode,vCode
	!character*132 readText,dumText
	!character*12  form
	character*251 inFile,dumString
	iErrCode=0
    ndbLocal=3000
	CrIndex=ndb ! vector initialize to ndb. if CrIndex(idDippr)==ndb, compd was not found in ParmsCrit.txt.
	inFile=TRIM(PGLinputDir)//'\ParmsPrTcJaubert.txt' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: CritFile=',TRIM(inFile)
!	OPEN(40,FILE=inFile,FORM='BINARY')
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
	READ(40,'(a251)',ERR=861)dumString
	READ(dumString,*,ERR=861)NDECK1
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
!		READ (40,ERR=861)IDnum(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
		READ (40,'(a222)',ioStat=ioErr)dumString
!		READ (dumString,*,ioStat=ioErr)IDnum(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
!			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,idCasDb(I) !,tCode,pCode,vCode,form,NAMED(I)
!		READ (dumString,'(a127,3a4,a12,a30)')readText,tCode,pCode,vCode,form,NAMED(I)
		READ (dumString,*)IDnum(I),TCD(I),PcTemp,ACEND(I),TbD(i),TwuL,TwuM,TwuN,cVt,ZCD(I),Tmin,idCasDb(I), &
		                                                              solParmD(i),rhoG_cc,rMwD(i),classDb(i),formDb(i),NAMED(i)
!1	190.56	4.599	0.0115	0.1473	0.9075	1.8243	-3.5604	0.2894	85	 74828	11.62	0.4224	16.04
		PCD(I)=PcTemp !/10
		if(rhoG_cc < zeroTol)rhoG_cc=1
		vLiqD(i)=rMwD(i)/rhoG_cc
		CrIndex(IDnum(i))=i
		if(ioErr.ne.0.and.LOUD)write(dumpUnit,*)'LoadCritParmsDb: error reading ParmsCrit.txt. line=',i
		if(i==1.and.LOUD)then
			write(dumpUnit,602)IDnum(I),idCasDb(I),TCD(I),PCD(I),ZCD(I),ACEND(I),solParmD(I),vLiqD(i),rMwD(i) !&
			write(dumpUnit,*)'LoadCrit: classDb,formDb,nameD=',classDb(i),formDb(i),NAMED(i)
		endif
		!	,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo

100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.3,i10,3a4,1x,a12,a33)
602	format(i5,i11,11f10.3,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor      &
	                                                                                          ! Cas#  tC  pC  vC  FORM        Name 
	CLOSE(40)
	!if((nDeck1).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in file'
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: So far so good! ParmsCrit.txt is loaded. Skipping ParmsCrAdd.'
	nDeckDb=NDECK1
	if(nDeckDb.ne.nCritSet.and.LOUD)write(dumpUnit,*) 'LoadCritParmsDb: nCritSet is inconsistent.nDeck,nCrit=',nDeckDb,nCritSet
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: Success! DB is loaded.'
	!if((nDeck).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in ParmsCrAdd file'
	!write(dumpUnit,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	RETURN
861	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=11
	!write(dumpUnit,*)
	close(40)
	return                      
862	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=12
	close(40)
	!write(dumpUnit,*)
	return                      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GETCRIT(NC,iErrCode)	! ID() from USEd GlobConst
	!C  
	!C  PROGRAMMED BY:  A.S. PUHALA - AUG 93
	!C  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	!C  REVISION DATE:  2/02 - FOR binary data file
	!C 
	!C  LOOKS UP THE CRITICAL PROPERTIES AND RETURNS JUST TC,PC,ACEN,NAME
	!C
	!C  INPUT
	!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!C  OUTPUT
	!C    TC - CRITICAL TEMPERATURE
	!C    PC - CRITICAL PRESSURE
	!C    ACEN - ACENTRIC FACTOR
	!C    NAME - COMPONENT NAME
	USE GlobConst
	USE CritParmsDb
	Implicit DoublePrecision(A-H,O-Z)
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	!	eHbKcalMol(nmx),bondVolNm3(nmx),ND(nmx),NDS(nmx),NAS(nmx)
	if(nDeckDb.ne.nCritSet)Call LoadCritParmsDb(iErrLoadCrit)
	if(iErrLoadCrit > 0)then
		if(LOUDER)write(dumpUnit,*)'GetCrit: error from LoadCritParmsDb. Thats all folks!'
		return
	endif
	iErrCode=0
	DO iComp=1,NC
		J=CrIndex(ID(iComp)) !USEd from CritParmsDb
		if(J==ndb)then		 !This is the indicator that ID(iComp) is not found in Db. cf. LoadCritParmsDb
			iGotIt=0
		else
			NAME(iComp)=NAMED(J)
			Tc(iComp)=TCD(J)
			Tb(iComp)=TbD(J)
			Pc(iComp)=PCD(J)
			ID(iComp)=IDnum(J)
			ACEN(iComp)=ACEND(J)
			ZC(iComp)=ZCD(J)
			rMw(iComp)=rMwD(j)
			solParm(iComp)=solParmD(j)
			vLiq(iComp)=vLiqD(j)
			class(iComp)=classDb(j)
			idCas(iComp)=idCasDb(j)
			iGotIt=1
		endif
		if(iGotIt==0)then
			iErrCode=iErrCode*100+iComp*10
			!write(dumpUnit,*)'Error in GetCrit: not found for iComp=',iComp
			!write(dumpUnit,*)
		endif
		if(LOUDER)write(dumpUnit,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp),Tc(iComp),Pc(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
	END	!GetCrit()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetVpDb(iErrCode)
	!	PROGRAMMED BY: AV 06/22/06
	!THIS SUBROUTINE GIVES VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	!INPUT:
	!NC - NUMBER OF COMPONENTS
	!	   ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!OUTPUT:
	!VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
    USE GlobConst, only: PGLinputDir,LOUD,dumpUnit
    USE VpDb ! Stores all the coefficients.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    character*255 inFile
	iErrCode=0
    inFile=TRIM(PGLinputDir)//'\CoeffsVp2a.TXT'
    OPEN(662,FILE=inFile,ioStat=ioErr)
    if(ioErr/=0.and.LOUD)write(dumpUnit,*) 'GetVpDb: error opening CoeffsVp2a.txt'
    !open(662,file='junk.txt')
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	!ndeck1=1974
	READ(662,*,ioStat=ioErr)NDECK1 ! ioStat= -1,endofFile; -2,endOfLine
    if(ioErr/=0)write(dumpUnit,*)'ioErr,nDeck1,inFile',ioErr,nDeck1,TRIM(inFile)
    if(ioErr/=0)write(dumpUnit,*) 'GetVpDb: Error reading NDECK1'
    
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetVp: more data in file than allocated'
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
		READ(662,*,ioStat=ioErr)IDnum(I),rMINTD(I) ,VALMIND(I) ,rMAXTD(I),VALMAXD(I),AVGDEVD(I),NUMCOEFFD(I), &
		                                                                                      (vpCoeffsd(i,iCoeff),iCoeff=1,5)  
	    if(ioErr/=0.and.LOUD)write(dumpUnit,*) 'inFile=',TRIM(inFile)
	    if(ioErr/=0.and.LOUD)write(dumpUnit,*) 'ioErr,line,id,vpCoeffs=',ioErr,I,IDnum(I),(vpCoeffsd(i,iCoeff),iCoeff=1,5)
	    if(ioErr/=0.and.LOUD)write(dumpUnit,*) 'GetVp: error reading CoeffsVp2a.txt'
		indexVpDb(IDnum(I))=I ! vpCoeffs(iComp,iCoeff)=vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)
    enddo
    if(LOUD)write(dumpUnit,*)'GetVpDb: Success! USE VpDb for vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)'
	CLOSE(662)
    return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetVp(NC,ID,iErrCode) !returns VpCoeffs(NC,5) in USEd VpDb
	!	PROGRAMMED BY: AV 06/22/06
	!THIS SUBROUTINE GIVES VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	!INPUT:
	!NC - NUMBER OF COMPONENTS
	!	   ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!OUTPUT:
	!VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	USE GlobConst, only:DEBUG,nmx,dumpUnit
	USE VpDb
	LOGICAL notFound
	!DoublePrecision vpCoeffs(nmx,5)
	Integer ID(nmx)   
	iErrCode=0

	DO iComp=1,NC
		notFound=.TRUE.
		if(indexVpDb(ID(iComp)) > 0)then
			do iCoeff=1,5
				vpCoeffs(iComp,iCoeff)=vpCoeffsd(indexVpDb(ID(iComp)),iCoeff)
			enddo
			notFound=.FALSE.
		endif
		if(notFound)then
			iErrCode=iErrCode*10+iComp
			!write(dumpUnit,*)'Error in GetCrit: not found for iComp=',iComp
			!write(dumpUnit,*)
		endif
		!write(dumpUnit,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
		!,Tc(iComp),Pc(iComp),acen(iComp),ZC(iComp)
	enddo

	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	integer function GetBIPs(bipFile,idComp,NC)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION (A-H,K,O-Z)
	PARAMETER (listPool=1000)
	logical switched
	character bipFile*88, dumString*88
	!integer GetBIPs
	dimension idBinarY(listPool),idComp(NC),alphaDB(listPool),KIJDB(listPool),KTIJDB(listPool),&
	wsTAUij (listPool),wsTAUji (listPool),wsTauTij(listPool),wsTauTji(listPool)
	!dimension ks0ij(nmx,nmx),ks1ij(nmx,nmx)
	!data initial/0/
	initial=0 !just assume that the database needs to be reloaded if being called

	GetBIPs=0	
	if(initial.eq.0)then
		initial=1
		open(55,file=bipFile,ERR=861)
		read(55,'(a88)',ERR=861)dumString
		item=1
		do while(item.ge.0)
			read(55,*,ERR=861,end=100) idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),&
				wsTAUji(item),wsTauTij(item),wsTauTji(item),alphaDB(item)
			item=item+1	!increment for subsequent loop if any.
			cycle
100			exit
		enddo
		close(55)
	endif
	numBips=item
	!write(dumpUnit,*)'# of bips in database =',numBips

	do iComp=1,NC
		do jComp=iComp+1,NC
			idBin=10000*idComp(iComp)+idComp(jComp)
			switched=.FALSE.
			if(idComp(iComp).lt.idComp(jComp))then
				switched=.TRUE.
				idBin=10000*idComp(jComp)+idComp(iComp)
			endif

			itemFound=0
			do item=1,numBips
				if(idBinarY(item)==idBin)then
					itemFound=item
					exit ! stops searching after finding the first instance.
				endif
			enddo

			xsAlpha(iComp,jComp)=0.3d0
			KIJ (iComp,jComp)=0
			KTIJ(iComp,jComp)=0
			xsTau (iComp,jComp)=0
			xsTau (jComp,iComp)=0
			xsTauT(iComp,jComp)=0
			xsTauT(jComp,iComp)=0
			kETAij(iComp,jComp)=0	! For speadmd.
			KS0IJ(iComp,jComp)=0
			KS1IJ(iComp,jComp)=0

			if(itemFound > 0)then
				xsAlpha(iComp,jComp)=alphaDB(itemFound)
				KIJ(iComp,jComp)=KIJDB(itemFound)
				KTIJ(iComp,jComp)=KTIJDB(itemFound)

				xsTau (iComp,jComp)=wsTAUij(itemFound)
				xsTau (jComp,iComp)=wsTAUji(itemFound)
				xsTauT(iComp,jComp)=wsTauTij(itemFound)
				xsTauT(jComp,iComp)=wsTauTji(itemFound)
				if(switched)then
					xsTau (jComp,iComp)=wsTAUij(itemFound)
					xsTau (iComp,jComp)=wsTAUji(itemFound)
					xsTauT(jComp,iComp)=wsTauTij(itemFound)
					xsTauT(iComp,jComp)=wsTauTji(itemFound)
				endif
				xsAlpha(jComp,iComp)=xsAlpha(iComp,jComp)
				KIJ(jComp,iComp)=KIJ(iComp,jComp)
				KTIJ(jComp,iComp)=KTIJ(iComp,jComp)
			endif
		enddo
	enddo

	return
861	continue
	if(LOUD)write(dumpUnit,'(a33,a22)')'GetBip error - error reading ',bipFile
	if(LOUD)write(dumpUnit,*)'numBips,item',numBips,item
	if(LOUD)write(dumpUnit,*)'Data read as:'
	if(LOUD)write(dumpUnit,'(i4,7F8.4)')idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),wsTAUji(item), &
	                                                                                  wsTauTij(item),wsTauTji(item),alphaDB(item)
	GetBIPs= -1
	return                      
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdCcLookup(NC,idcc,ier,errMsgPas)
    USE GlobConst
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas
	character inFile*251
	dimension idcc(NC)
	dimension idCcDb(maxDb),idDb(maxDb)
		inFile=TRIM(PGLinputDir)//'\idTrcDipCas.TXT' ! // is the concatenation operator
		OPEN(50,FILE=inFile)
	errMsg(0)='No Problem'
	errMsg(1)='IdCcLookup Error: at least one id not found'
	ier=0
	i=0
	ID(1)=0
	do i=1,maxDb
		read(50,*,END=101)idCcDb(i),idDb(i)
		if(idCcDb(i).eq.idcc(1))ID(1)=idDb(i)
		cycle
101		continue
		numDb=i-1
		exit !terminate do loop
	enddo
	if(id(1).eq.0)then
		ier=1
		if(LOUD)write(dumpUnit,*)'GetCrit:Did not find idcc(1)=',idcc(1)
		errMsgPas=errMsg(ier)
		return
	endif
	do iComp=2,NC
		ID(iComp)=0
		iGotIt=0
		i=0
		do while(iGotIt.eq.0.and.i.lt.numDb)
			i=i+1
			if(idCcDb(i).eq.idcc(iComp))then
				ID(iComp)=idDb(i)
				iGotit=1
			endif
		enddo
		if(id(iComp).eq.0)then		
			ier=1
			if(LOUD)write(dumpUnit,*)'Did not find idcc(i).i,idcc=',iComp,idcc(1)
			errMsgPas=errMsg(ier)
			return
		endif
	enddo
	close(50)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdCasLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst, only:ID,idTrc,nmx,class,DEBUG,LOUD,PGLinputDir,dumpUnit,name
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas,dumString
	character*251 inFile
	character*28 NameDb(maxDb)
	!character*5 classdb(maxDb)
	integer idCas(nmx)
	integer idCasDb(maxDb),idDb(maxDb),idTrcDb(maxDb)
    !OPEN(50,FILE='c:\spead\idTrcDipCas.TXT')
		inFile=TRIM(PGLinputDir)//'\idTrcDipCas.TXT' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,*)'idCasLookup: inFile=',TRIM(inFile)
	OPEN(50,FILE=inFile)
	errMsg(0)='No Problem'
	errMsg(1)='IdCasLookup Error: at least one id not found'
	ier=0
	iGotIt=0
	read(50,*)nDeck
	!write(dumpUnit,*)'IdLookup: nDeck=',nDeck
	do i=1,nDeck
		read(50,'(a77)')dumString
		read(dumString,'(2i6,i12,1x,a28,a5)',ioStat=ioErr)idTrcDb(i),idDb(i),idCasDb(i),NameDb(i) !,classdb(i)  
		if(ioErr/=0.and.LOUD)write(dumpUnit,*)'idCas ioErr. i,idDippr=',i,idDb(i)
		if( id(1)==idDb(i) )then
			iGotIt=1
			idCas(1)=idCasDb(i)
			idTrc(1)=idTrcDb(i)
			!class(1)=classdb(i)
			name(1)=TRIM(NameDb(i))
			if(LOUD)write(dumpUnit,*)TRIM(dumString)
			if(LOUD)write(dumpUnit,*)'IdCasLookup:id1,name,class1:',id(1),TRIM(name(1)) !,class(1)
		endif
	enddo
	if(iGotIt==0)then
		ier=1
		!write(dumpUnit,*) 'IdCasLookup: idCas not found for comp1'
		goto 861
	endif
	if(NC==1)goto 861
	do iComp=2,NC
		iGotIt=0
		do i=1,nDeck
			if(idDb(i).eq.id(iComp))then
				idCas(iComp)=idCasDb(i)
				idTrc(iComp)=idTrcDb(i)
				!class(iComp)=classdb(i)
				name(iComp)=TRIM(NameDb(i))
				if(LOUD)write(dumpUnit,*)'IdCasLookup:idi,idCasi:',id(iComp),idCas(iComp) !,class(i)
				iGotit=1
				exit
			endif
		enddo
		if(iGotIt==0)then 
			ier=1 
			if(LOUD)write(dumpUnit,*)'IdCasLookup:Did not find idCas(i).i,id=',iComp,id(i)
			goto 861
		endif
	enddo
861	errMsgPas=errMsg(ier)
	close(50)
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine PrintErrFlash(iErrCode)
	USE GlobConst
	if(LOUD)write(dumpUnit,*)'PrintErrFlash:iErrCode=',iErrCode
	IF(iErrCode.EQ.1)write(dumpUnit,*)'LLeFL:Initial guess REQUIRED.'
	IF(iErrCode.EQ.2)write(dumpUnit,*)'Flash:ERROR ALL UPPER PHASE'
	IF(iErrCode.EQ.3)write(dumpUnit,*)'Flash:ERROR ALL LOWER PHASE'
	IF(iErrCode.EQ.4)write(dumpUnit,*)'Flash:LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.5)write(dumpUnit,*)'Flash:VAPOR  FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.6)write(dumpUnit,*)'Flash:ITMAX EXCEEDED		 '
	IF(iErrCode.EQ.7)write(dumpUnit,*)' Flash:xFrac(i) < 0 for some i'
	IF(iErrCode.EQ.11)write(dumpUnit,*)' zV~zL - Flash: warning'
	!IF(iErrCode.EQ.21)	!NOTE: not relevant for polymers because vLiq >> 1 and p >> pSat
	!& write(dumpUnit,*)' zL > 0.33 - Flash:warning '
	IF(iErrCode.EQ.22)write(dumpUnit,*)' zV < 0.33 - Flash:warning'
	IF(iErrCode.EQ.23)write(dumpUnit,*)' Flash warning:gEx < 0 on at least one LLE iteration.'
	!C
	!write(dumpUnit,*) 'PrintErrFlash: Check errors'
	return
	end


SUBROUTINE GetCritCas(NC,idCas,iErrCode)	!GetCritCas assumes ID(GlobConst)=IdCas
	!**********************************************
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!C      AARON'S SUBROUTINE
	!C      $GETCRIT.FOR   VERSION 1.4
	!C  
	!C  PROGRAMMED BY:  A.S. PUHALA - AUG 93
	!C  REVISION DATE:  2/93 - FOR ESD COMPATIBILITY JRE
	!C  REVISION DATE:  2/02 - FOR binary data file
	!C                  1/23 - To use CritParmsDb and GetCrit()
	!C 
	!C  LOOKS UP THE CRITICAL PROPERTIES AND RETURNS JUST TC,PC,ACEN,NAME
	!C
	!C  INPUT
	!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!C  OUTPUT
	!C    TC - CRITICAL TEMPERATURE
	!C    PC - CRITICAL PRESSURE
	!C    ACEN - ACENTRIC FACTOR
	!C    NAME - COMPONENT NAME
	USE GlobConst, ONLY:LOUD,dumpUnit,nmx 
	USE CritParmsDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	character*77 errMsgPas
	integer idCas(nmx),NC,iErrCode,iErrCrit
	iErrCode=0
	Call IdDipprLookup(NC,idCas,iErrCode,errMsgPas) ! ID passed through USEd GlobConst
	if(iErrCode.ne.0)then 
		if(LOUD)write(dumpUnit,*)'GetCritCas: iErr(IdDipprLookup)=',iErrCode,TRIM(errMsgPas)
		return
	endif
	Call GetCrit(NC,iErrCrit)
	if(iErrCrit.ne.0)then 
		if(LOUD)write(dumpUnit,*)'GetCritCas: iErr(GetCrit)=',iErrCrit
		iErrCode=iErrCrit
		return
	endif
	return                      
END ! subroutine GetCritCas()

      DOUBLE PRECISION FUNCTION AvgAbsRelDev(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N),Y(N)
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      AvgAbsRelDev=0
      DO I=1,N
          if( ABS(Y(I)) < 1D-33 )then
              AvgAbsRelDev=86.8686
              return
          else
              AvgAbsRelDev=DABS( X(I)-Y(I)/Y(I) )
          endif
      ENDDO
      AvgAbsRelDev=AvgAbsRelDev/N
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION BiasRelDev(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N),Y(N)
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      BiasRelDev=0
      DO I=1,N
          if( ABS(Y(I)) < 1D-33 )then
              BiasRelDev=86.8686
              return
          else
              BiasRelDev=( X(I)-Y(I)/Y(I) )
          endif
      ENDDO
      BiasRelDev=BiasRelDev/N
      RETURN
      END
      
     
      DOUBLE PRECISION FUNCTION AvgAbsDev(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N),Y(N)
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      AvgAbsDev=0
      DO I=1,N
          AvgAbsDev=DABS( X(I)-Y(I) )
      ENDDO
      AvgAbsDev=AvgAbsDev/N
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION SumProduct(N,X,Y)
	USE GlobConst !for LOUD
      INTEGER N
      DOUBLE PRECISION X(N),Y(N)
	!zeroTol=1.D-11 !zeroTol is now global
      SumProduct=0
      do i=1,N
          SumProduct=SumProduct+X(i)*Y(i)
      enddo
	if( ABS(SumProduct) < zeroTol .and. LOUD)then
		write(dumpUnit,*)'SumProduct: N,x1,y1',N,X(1),Y(1)
	endif
      return
      end
      
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      INTEGER N
      DOUBLE PRECISION X(N),SqrtArg,SumSq
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      SqrtArg=SumSq(N,X)
      if(SqrtArg .lt. 0.D0)then
          ENORM=86.8686
      else
          ENORM=SQRT(SqrtArg)
      endif
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION RmsJre(N,X)
      INTEGER N
      DOUBLE PRECISION X(N),SqrtArg,SumSq
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      SqrtArg=SumSq(N,X)
      if(SqrtArg < 0 .or. N < 1)then
          RmsJre=86.8686
      else
          RmsJre=DSQRT(SqrtArg/N)
      endif
      
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION SumXmY2(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N),Y(N),DIF(N),SumSq
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      DO I=1,N
          DIF(I)=X(I)-Y(I)
      ENDDO
      SumXmY2=SumSq(N,DIF)
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION RmsDev(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N),Y(N),SqrtArg,SumSq
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
      DO I=1,N
          X(I)=X(I)-Y(I)
      ENDDO
      SqrtArg=SumSq(N,X)/N
      if(SqrtArg < 0)then
          RmsDev=86.8686
      else
          RmsDev=DSQRT(SqrtArg)
      endif
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION SumSq(N,X)
	  USE GlobConst !LOUD
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(N)
!C     **********
!C
!C     FUNCTION ENORM
!C
!C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!C     EUCLIDEAN NORM OF X.
!C
!C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
!C     SQUARES IN THREE DIFFERENT SUMS. THE SUMS OF SQUARES FOR THE
!C     SMALL AND LARGE COMPONENTS ARE SCALED SO THAT NO OVERFLOWS
!C     OCCUR. NON-DESTRUCTIVE UNDERFLOWS ARE PERMITTED. UNDERFLOWS
!C     AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!C     SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!C     THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS
!C     DEPEND ON TWO CONSTANTS, RDWARF AND RGIANT. THE MAIN
!C     RESTRICTIONS ON THESE CONSTANTS ARE THAT RDWARF**2 NOT
!C     UNDERFLOW AND RGIANT**2 NOT OVERFLOW. THE CONSTANTS
!C     GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.
!C
!C     THE FUNCTION STATEMENT IS
!C
!C       DOUBLE PRECISION FUNCTION ENORM(N,X)
!C
!C     WHERE
!C
!C       N IS A POSITIVE INTEGER INPUT VARIABLE.
!C
!C       X IS AN INPUT ARRAY OF LENGTH N.
!C
!C     SUBPROGRAMS CALLED
!C
!C       FORTRAN-SUPPLIED ... DABS,DSQRT
!C
!C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!C
!C     **********
      INTEGER I
      DOUBLE PRECISION AGIANT,FLOATN,ONE,RDWARF,RGIANT,S1,S2,S3,XABS,X1MAX,X3MAX,ZERO !,SQARG
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = DABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
!C
!C              SUM FOR LARGE COMPONENTS.
!C
               IF (XABS .LE. X1MAX) GO TO 10
				if((XABS.EQ.0).AND.LOUD)write(dumpUnit,*) 'ENORM:XABS=0'
                  S1 = ONE + S1*(X1MAX/XABS)**2
				!if(X1MAX.EQ.0)write(dumpUnit,*) 'ENORM:X1MAX=0'
                  X1MAX = XABS
                  GO TO 20
10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
20          CONTINUE
               GO TO 60
30       CONTINUE
!C
!C              SUM FOR SMALL COMPONENTS.
!C
               IF (XABS .LE. X3MAX) GO TO 40
				if((XABS.EQ.0).AND.LOUD)write(dumpUnit,*) 'SumSq:XABS=0'
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
40          CONTINUE
				!if(X3MAX.EQ.0)write(dumpUnit,*) 'ENORM:X3MAX=0'
                  IF ( X3MAX.NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
50          CONTINUE
60       CONTINUE
            GO TO 80
70    CONTINUE
!C
!C           SUM FOR INTERMEDIATE COMPONENTS.
!C
            S2 = S2 + XABS**2
80    CONTINUE
90    CONTINUE
!C
!C     CALCULATION OF NORM.
!C
!      IF (S1 .EQ. ZERO) GO TO 100
	IF (S1 .NE. ZERO)THEN
         SumSq = X1MAX*X1MAX*(S1+(S2/X1MAX)/X1MAX)
!         GO TO 130
!  100 CONTINUE
      ELSE IF (S2 .NE. ZERO) THEN
            IF (S2 .GE. X3MAX)THEN
			SumSq=S2*(ONE+(X3MAX/S2)*(X3MAX*S3))
            ELSE
			SumSq=X3MAX*((S2/X3MAX)+(X3MAX*S3))
            ENDIF
		  !GO TO 120
  !110    CONTINUE
	ELSE
            SumSq = X3MAX*X3MAX*S3
	ENDIF
  !120    CONTINUE
  !130 CONTINUE
      RETURN
!C
!C     LAST CARD OF FUNCTION SumSq.
!C
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer Function ItsEven(iArg) ! returns 0 for odd and 1 for even.
    integer iArg
	ItsEven=1-MOD(iArg,2)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	LOGICAL Function bEven(iArg) ! returns 0 for odd and 1 for even.
    integer iArg,iTest
	bEven=.FALSE.
	iTest=(iArg/2.d0-iArg/2)*2.001d0
	if(iTest==1)bEven=.TRUE.
	return
	end

