!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GlobConst
	!SAVE
	!PUBLIC
	Implicit NONE
	DoublePrecision pi,twoPi,fourPi,AvoNum,Rgas,RgasCal,kB,half,zeroTol
	Integer NMX 
	PARAMETER (nmx=55,pi=3.14159265359796d0,twoPi=2.d0*pi, fourPi = 4.d0*pi, half=0.5d0)
	PARAMETER (AvoNum=602.214076d0,kB=0.01380649D0,Rgas=AvoNum*kB,RgasCal=Rgas/4.184d0,zeroTol=1.D-12)!kB[=]MPa.nm3/K) AvoNum[=]cm3/(nm3*mol), 
	!          https://www.nist.gov/si-redefinition 
    !nmx is the max allowed number of Compounds
	!integer :: idComp(nmx),nsTypes(nmx),IDs(nsx),IDsBase(nmx,nsx),siteNum(nmx,maxTypes)
	!integer :: nComps, nsTypesTot,iTPT,iFlagFF,nNormGrid

	DoublePrecision :: TC(nmx), PC(nmx), ACEN(nmx), ZC(nmx), rMw(nmx), bVolCC_mol(NMX),solParm(NMX),vLiq(NMX)
	DoublePrecision uRes_RT, sRes_R, aRes_RT, hRes_RT, cpRes_R, cvRes_R, cmprsblty !cmprsblty=(dP/dRho)T*(1/RT)	= Z+rho*dZ/dRho
	character*234 masterDir,PGLinputDir
	character*30 NAME(NMX)
	character*5 class(NMX) ! Allowed: norml,heavy,polar,assoc,Asso+,gases,siloa,salty,ormet,metal,inorg (cf. Ch06CompoundsList.xls, PGL6edClasses.xls)
    Logical LOUD   !LOUD=.TRUE. means writing debug info to the screen.
	LOGICAL DEBUG, isESD, isTPT, isPcSaft 
	integer ID(nmx), idCas(nmx), idTrc(nmx), iEosOpt, initEos, dumpUnit 
	DoublePrecision etaPass !TODO:Consider if rho needs to be a calling argument of any FUGI(). Otherwise, precision in rho is lost when Z->0 by passing zFactor then computing rho, esp for PREOS (incl Jaubert version).
    DoublePrecision etaMax  !each EOS has a max value for eta, e.g. PR,TPT: etaMax=1-zeroTol. This must be set in the Get_ function for the EOS
	DoublePrecision etaPure(NMX) !store eta for each compound at P=0 and Tmin(K). 
	Character*15 EosName(18) 
	!               1     2       3       4          5          6         7           8              9            10          11      12        13         14          15           16            17        18
	data EosName/'PR','ESD96','PRWS','ESD-MEM2','SPEADMD','Flory-MEM2','NRTL','SpeadGamma-MEM2','SPEAD11','PcSaft(Gross)','tcPRq','GCESD','GCESD(Tb)','TransSPEAD','GcPcSaft','GcPcSaft(Tb)','tcPR-GE(W)','ESD2'/
	!LOUD = .TRUE.		  !!!!!!!!!!!!!!! YOU CAN'T SET VARIABLES IN A MODULE, ONLY PARAMETERS AND DATA !!!!!!!!!!!!!
	!LOUD = .FALSE.
contains
	integer function SetNewEos(newEosOpt)
	!returns 0 if no error.
	!USE GlobConst ! eliminated because GlobConst contains SetNewEos.
	implicit none
	integer newEosOpt
	SetNewEos=0
	isESD=.FALSE.
	isTPT=.FALSE.
	isPcSaft=.FALSE.
	iEosOpt=newEosOpt
	return
	end function SetNewEos
END MODULE GlobConst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CritParmsDb
	Integer ndb
	Parameter (ndb=3000)
	character*30 NAMED(ndb)	!LoadCrit() loads ParmsCrit database.
	Integer IDNUM(ndb),CrIndex(9999),idCasDb(ndb),nDeckDb ! e.g. TCD(CrIndex(2)) returns Tc of ethane. 
	DoublePrecision TCD(ndb),PCD(ndb),ACEND(ndb),ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb) ! LoadCrit uses CrIndex to facilitate lookup. TCD(ndb)=8686. CrIndex()=ndb initially.
END MODULE CritParmsDb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE BIPs
	USE GlobConst, only: NMX 
	integer nConverged,maxPts,nPtsBipDat
	parameter(maxPts=1777) ! this defines the max allowed # of experimental data points in a single binary system.
	DoublePrecision KIJ(NMX,NMX),KTIJ(NMX,NMX),kETAij(NMX,NMX) !usual dispersive BIPs & k^eta_ij
	DoublePrecision KS0IJ(NMX,NMX),KS1IJ(NMX,NMX)              !entropic BIPs & k^eta_ij
	DoublePrecision HIJ(NMX,NMX),HTIJ(NMX,NMX) !molecular hBonding BIPs for ESD. (Spead aBipAd,aBipDa are site based.)
	DoublePrecision Lij(NMX,NMX) !covolume adjustment.  bVolMix=sum(sum(xi*xj*bij)); bij=(1-Lij)*(bi+bj)/2
	DoublePrecision xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)	!this is for the PRWS/xsNRTL mixing rule.
	DoublePrecision tDat(maxPts),PDAT(maxPts),XDAT(maxPts),YDAT(maxPts) !,deviate(maxPts) !can't include Deviate() cuz it's an argument for LmDif.
	DoublePrecision pDatMin,pDatMax
	Integer iDat(maxPts)  ! sometimes need to indicate whether the data are for comp1 or comp2. e.g. SLE.
END MODULE BIPs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE VpDb
	USE GlobConst, only:NMX,dumpUnit
	IMPLICIT NONE !DoublePrecision(A-H,O-Z)
	Integer nVpDb
	PARAMETER(nVpDb=2475)
	Integer IDNUM(nVpDb),NUMCOEFFD(nVpDb) ,indexVpDb(9999)
	DoublePrecision rMINTD(nVpDb) ,VALMIND(nVpDb) ,rMAXTD(nVpDb),VALMAXD(nVpDb),AVGDEVD(nVpDb),vpCoeffsd(nVpDb,5)
	DoublePrecision vpCoeffs(NMX,5)
END MODULE VpDb

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
	USE GlobConst, Only:DEBUG,PGLInputDir,LOUD,dumpUnit 
	! Purpose: Read Tm(K) and Hfus(J/mol) from input file.
	Implicit DoublePrecision(a-h,o-z)
	Character*77 HfusFile,FN !,inFile77
    FN="Hfus.txt"
	HfusFile=TRIM(PGLInputDir)//'\'//TRIM(FN)
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
		if(ioErr)cycle
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
	if(notFound==1)iErr=1
	if(ioErr /= 0)iErr=2
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
	subroutine IdDipprLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst, only:ID,nmx,LOUD,dumpUnit
	USE CritParmsDb, only:idCasDb,idNum,CrIndex,nDeckDb
    integer ier,	idCas(nmx)
	character*77 errMsg(0:11),errMsgPas
	errMsg(0)='No Problem'
	errMsg(1)='IdDipprLookup Error: at least one id not found'
	ier=0
	!write(dumpUnit,*)'IdLookup: nDeck=',nDeck
	do iComp=1,NC
		iGotIt=0
		do i=1,nDeckDb
			if(idCasDb(i)==idCas(iComp))then
				id(iComp)=idNum(i)
				!class(iComp)=classdb(i)	!class() is not available. Call idCasLookup() if you need class().
				iGotit=1
				exit ! breaks the inner loop only?
			endif
		enddo
		if(iGotIt==0)ier=1 !only way here is if cycle failed => iGotIt==0.
		if(LOUD.and.ier>0)write(dumpUnit,*)'Did not find idCas(i).i,id=',iComp,id(i)
		if(ier > 0)goto 861
	enddo
861	errMsgPas=errMsg(ier)
	close(50)
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
	!C    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!C  OUTPUT
	!C    TC - CRITICAL TEMPERATURE
	!C    PC - CRITICAL PRESSURE
	!C    ACEN - ACENTRIC FACTOR
	!C    NAME - COMPONENT NAME
	USE GlobConst
	USE CritParmsDb 
	IMPLICIT DoublePrecision(A-H,O-Z)
	CHARACTER*4 tCode,pCode,vCode
	character readText*132,dumText*77,form*12 
	character*251 infile,dumString
	iErrCode=0
	CrIndex=ndb ! vector initialize to ndb. if CrIndex(idDippr)==ndb, compd was not found in ParmsCrit.txt.
	inFile=TRIM(PGLInputDir)//'\ParmsCrit.txt' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: CritFile=',TRIM(inFile)
!	OPEN(40,FILE=inFile,FORM='BINARY')
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
!	READ(40,ERR=861)NDECK1
	READ(40,*,ERR=861)NDECK1
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
!		READ (40,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
		READ (40,'(a188)',ioStat=ioErr)dumString
		READ (dumString,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,idCasDb(I),tCode,pCode,vCode,form,NAMED(I)
		CrIndex(IDNUM(i))=i
		if(ioErr.and.LOUD)write(dumpUnit,*)'GetCrit: error reading ParmsCrit.txt. line=',i
		if(i==1.and.LOUD)write(dumpUnit,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) !&
		!	,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo

100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.3,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	CLOSE(40)
	!if((nDeck1).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in file'
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: So far so good! ParmsCrit.txt is loaded. Trying ParmsCrAdd.'
		inFile=TRIM(PGLinputDir)//'\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	!ENDIF
	READ(40,*,ERR=862)dumText

	nDeck2=0
	DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ(40,'(A132)',ERR=862,END=200)readText
		nDeck2=nDeck2+1
		I=nDeck1+nDeck2
		read(readText,*,ERR=862)IDNUM(I),TCD(I),PCD(I),ZCD(I) &
			,ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		if(idNum(i) > 9999)then
			nDeck2=nDeck2-1 ! Take one step back!
			cycle ! LoadCritParmsDb does not permit this because CrIndex is limited to 9999.
		endif
		CrIndex(IDNUM(i))=i
		indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
		cycle
200		EXIT !terminate loop
	enddo
	close(40)
	nDeckDb=nDeck1+nDeck2
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: Success! DB is loaded.'
	!if((nDeck).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in ParmsCrAdd file'
	!write(dumpUnit,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	RETURN
861	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
862	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GETCRIT(NC,iErrCode)
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
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	CHARACTER NAMED*30,tCode*4,pCode*4,vCode*4,form*12
	character readText*132,dumText*77 
	PARAMETER(ndb=1555)
	character*251 infile,dumString
	DIMENSION IDNUM(ndb),TCD(ndb),PCD(ndb),ACEND(ndb),NAMED(ndb)
	DIMENSION ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb)
	!common/FloryWert/solParm(nmx),vLiq(nmx),vMolec(NMX),&
	!	eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	iErrCode=0
	inFile=TRIM(PGLinputDir)//'\ParmsCrit.txt' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,'(2a)')' CritFile=',TRIM(inFile)
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
!	READ(40,ERR=861)NDECK1
	READ(40,*,ERR=861)NDECK1
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
!		READ (40,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
		READ (40,'(a188)',ioStat=ioErr)dumString
		READ (dumString,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ (dumString,'(a126,3a4,a12,a30)',ioStat=ioErr)readText,tCode,pCode,vCode,form,NAMED(I)
		if(ioErr.and.LOUD)write(dumpUnit,*)'GetCrit: error reading ParmsCrit.txt. line=',i
		if(i==1.and.LOUD)write(dumpUnit,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) !&
		!	,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo

100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1)==0)THEN
		DO I=1,NDECK1,3
			IF(LOUD)write(dumpUnit,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
			!IF(((I-1)/45)*45.EQ.(I-1))write(dumpUnit,*)
		enddo
		RETURN
	END IF
	CLOSE(40)
	!if((nDeck1).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in file'

	inFile=TRIM(PGLinputDir)//'\ParmsCrAdd.TXT' ! // is the concatenation operator
	OPEN(40,FILE=inFile)
	READ(40,*,ERR=862)dumText

	nDeck2=0
	DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ(40,'(A132)',ERR=862,END=200)readText
		nDeck2=nDeck2+1
		I=nDeck1+nDeck2
		read(readText,*,ERR=862)IDNUM(I),TCD(I),PCD(I),ZCD(I) &
			,ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
		cycle
200		EXIT !terminate loop
	enddo
	close(40)
	nDeck=nDeck1+nDeck2
	!if((nDeck).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in ParmsCrAdd file'
	!write(dumpUnit,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	DO iComp=1,NC
		iGotIt=0
		DO J=1,NDECK
			IF(IDNUM(J).EQ.ID(iComp)) THEN
				NAME(iComp)=NAMED(J)
				TC(iComp)=TCD(J)
				PC(iComp)=PCD(J)
				ID(iComp)=IDNUM(J)
				ACEN(iComp)=ACEND(J)
				ZC(iComp)=ZCD(J)
				rMw(iComp)=rMwD(j)
				solParm(iComp)=solParmD(j)
				vLiq(iComp)=vLiqD(j)
				iGotIt=1
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			!write(dumpUnit,*)'Error in GetCrit: not found for iComp=',iComp
			!write(dumpUnit,*)
		endif
		if(LOUD)write(dumpUnit,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp),TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
862	continue
	!write(dumpUnit,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    if(ioErr.and.LOUD)write(dumpUnit,*) 'GetVpDb: error opening CoeffsVp2a.txt'
    !open(662,file='junk.txt')
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	!ndeck1=1974
	READ(662,*,ioStat=ioErr)NDECK1 ! ioStat= -1,endofFile; -2,endOfLine
    if(ioErr)write(dumpUnit,*)'ioErr,nDeck1,inFile',ioErr,nDeck1,TRIM(inFile)
    if(ioErr)write(dumpUnit,*) 'GetVpDb: Error reading NDECK1'
    
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetVp: more data in file than allocated'
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
		READ(662,*,ioStat=ioErr)IDNUM(I),rMINTD(I) ,VALMIND(I) ,rMAXTD(I),VALMAXD(I),AVGDEVD(I),NUMCOEFFD(I) ,(vpCoeffsd(i,iCoeff),iCoeff=1,5)  
	    if(ioErr.and.LOUD)write(dumpUnit,*) 'inFile=',TRIM(inFile)
	    if(ioErr.and.LOUD)write(dumpUnit,*) 'ioErr,line,id,vpCoeffs=',ioErr,I,IDNUM(I),(vpCoeffsd(i,iCoeff),iCoeff=1,5)
	    if(ioErr.and.LOUD)write(dumpUnit,*) 'GetVp: error reading CoeffsVp2a.txt'
		indexVpDb(IDNUM(I))=I ! vpCoeffs(iComp,iCoeff)=vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)
    enddo
    if(LOUD)write(dumpUnit,*)'GetVpDb: Success! USE VpDb for vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)'
	CLOSE(662)
    return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE GetVp(NC,ID,iErrCode) !returns VpCoeffs(NC,5)
	!	PROGRAMMED BY: AV 06/22/06
	!THIS SUBROUTINE GIVES VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	!INPUT:
	!NC - NUMBER OF COMPONENTS
	!	   ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!OUTPUT:
	!VAPOR PRESSURE COEFFICIENTS FROM DIPPR DATABASE
	USE GlobConst, only:DEBUG,NMX,dumpUnit
	USE VpDb
	LOGICAL notFound
	!DoublePrecision vpCoeffs(NMX,5)
	Integer ID(NMX)   
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
		!,TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
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
	!common/ksvall/ks0ij,ks1ij
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
			KS0IJ(iComp,jComp)=0
			KS1IJ(iComp,jComp)=0
			xsTau (iComp,jComp)=0
			xsTau (jComp,iComp)=0
			xsTauT(iComp,jComp)=0
			xsTauT(jComp,iComp)=0

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
	if(LOUD)write(dumpUnit,'(i4,7F8.4)')idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),wsTAUji(item),wsTauTij(item),wsTauTji(item),alphaDB(item)
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
	id(1)=0
	do i=1,maxDb
		read(50,*,END=101)idCcDb(i),idDb(i)
		if(idCcDb(i).eq.idcc(1))id(1)=idDb(i)
		cycle
101		continue
		numDb=i-1
		exit !terminate do loop
	enddo
	if(id(1).eq.0)then
		ier=1
		if(LOUD)write(dumpUnit,*)'Did not find idcc(1)=',idcc(1)
		errMsgPas=errMsg(ier)
		return
	endif
	do iComp=2,NC
		id(iComp)=0
		iGotIt=0
		i=0
		do while(iGotIt.eq.0.and.i.lt.numDb)
			i=i+1
			if(idCcDb(i).eq.idcc(iComp))then
				id(iComp)=idDb(i)
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
    USE GlobConst, only:ID,idTrc,nmx,class,DEBUG,LOUD,PGLinputDir,dumpUnit
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas,dumString
	character*251 inFile
	character*28 dumName
	character*5 classdb(maxDb)
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
		read(dumString,'(2i6,i12,1x,a28,a5)',ioStat=ioErr)idTrcDb(i),idDb(i),idCasDb(i),dumName,classdb(i)  
		if(ioErr.and.LOUD)write(dumpUnit,*)'idCas ioErr. i,idDippr=',i,idDb(i)
		if( id(1)==idDb(i) )then
			iGotIt=1
			idCas(1)=idCasDb(i)
			idTrc(1)=idTrcDb(i)
			class(1)=classdb(i)
			if(LOUD)write(dumpUnit,*)TRIM(dumString)
			if(LOUD)write(dumpUnit,*)'IdCasLookup:id1,name,class1:',id(1),TRIM(dumName),class(1)
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
				class(iComp)=classdb(i)
				if(LOUD)write(dumpUnit,*)'IdCasLookup:idi,classi:',id(i),class(i)
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


	!**********************************************
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!C      AARON'S SUBROUTINE
	!C      $GETCRIT.FOR   VERSION 1.4
	SUBROUTINE GETCRITCAS(NC,iErrCode)	!GetCritCas assumes ID(GlobConst)=IdCas

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
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	CHARACTER NAMED*30,tCode*4,pCode*4,vCode*4,form*12
	character*132 readText,dumText 
	PARAMETER(ndb=1555)
	character infile*251
	DIMENSION IDNUM(ndb),TCD(ndb),PCD(ndb),ACEND(ndb),NAMED(ndb),iCasd(ndb)
	DIMENSION ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb) !,vLiq(NMX)
	!common/FloryWert/solParm(nmx),vLiq(nmx),vMolec(NMX),&
	!	eHbKcalMol(NMX),bondVolNm3(NMX),ND(NMX),NDS(NMX),NAS(NMX)
	iErrCode=0
	inFile=TRIM(PGLinputDir)//'\ParmsCrit.txt' ! // is the concatenation operator
	!write(dumpUnit,*)'CritFile=',TRIM(inFile)
	!OPEN(40,FILE=inFile)
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
	READ(40,*,ERR=861)NDECK1
	write(dumpUnit,*)'nDeck1=',nDeck1
	write(dumpUnit,*)'ID(1)=',ID(1)
	!if(ndeck1.gt.ndb)write(dumpUnit,*) 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)write(dumpUnit,*)
		READ (40,'(a188)',ioStat=ioErr)readText
		READ (readText,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
!		write(61,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
!			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1).EQ.0)THEN
		DO I=1,NDECK1,3
			if(LOUD)write(dumpUnit,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
			!IF(((I-1)/45)*45.EQ.(I-1))write(dumpUnit,*)
		enddo
		RETURN
	ENDIF
	CLOSE(40)
	!if((nDeck).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in file'

		inFile=TRIM(PGLinputDir)//'\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	READ(40,*,ERR=862)readText	! clear the header

	nDeck2=0
	DO while(nDeck2.ge.0) !use while loop to avoid adding/incrementing compd counter at beginning
		!Note: formatted read of tab delimited text did not work(?) JRE 042806
		!ADVANCE feature did not work either b/c it can't accept unformatted read (list directed i/o). 
		!READ(dumString,101,ERR=862,ADVANCE='NO')IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ(40,'(A132)',ERR=862,END=200)readText
		nDeck2=nDeck2+1
		I=nDeck1+nDeck2
		read(readText,*,ERR=862)IDNUM(I),TCD(I),PCD(I),ZCD(I) &
			,ACEND(I),rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCasd(I)
		indChars=INDEX(readText,'P 5',BACK=.TRUE.)
		read(readText,'(a<indChars+3>,a12,a25)')dumText,form,NAMED(I)
		cycle
200		EXIT !terminate loop
	enddo
	close(40)
	nDeck=nDeck1+nDeck2
	!if((nDeck).gt.ndb)write(dumpUnit,*) 'GetCrit: too much data in ParmsCrAdd file'
	if(LOUD)write(dumpUnit,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	DO iComp=1,NC
		iGotIt=0
		DO J=1,NDECK
			IF(ICasd(J).EQ.ID(iComp)) THEN
				NAME(iComp)=NAMED(J)
				TC(iComp)=TCD(J)
				PC(iComp)=PCD(J)
				ID(iComp)=IDNUM(J)
				ACEN(iComp)=ACEND(J)
				ZC(iComp)=ZCD(J)
				rMw(iComp)=rMwD(j)
				solParm(iComp)=solParmD(j)
				vLiq(iComp)=vLiqD(j)
				iGotIt=1
				exit
			ENDIF
		enddo
		if(iGotIt.eq.0)then
			iErrCode=iErrCode*10+iComp
			!write(dumpUnit,*)'Error in GetCrit: not found for iComp=',iComp
			!write(dumpUnit,*)
		endif
		if(LOUD)write(dumpUnit,'(1x,i11,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
		,TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	if(LOUD)write(dumpUnit,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
862	continue
	if(LOUD)write(dumpUnit,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	if(LOUD)write(dumpUnit,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!write(dumpUnit,*)
	return                      
END ! subroutine GetCrit

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer Function ItsEven(iArg) ! returns 0 for odd and 1 for even.
    integer iArg
	ItsEven=1-MOD(iArg,2)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!int QueryNparPure();    //number of adjustable parameters for the current model and mixture
	integer Function QueryNparPure()
	USE GlobConst, only: iEosOpt !, class
	parameter(nModels=16)
!   cf. GlobConst for iEosOpt's as listed below
!    1     2      3          4        5        6         7         8          9         10      11     12       13         14          15         16             17        18
!	'PR','Esd96','PRWS','Esd-MEM2','SPEAD','FloryWert','NRTL','SpeadGamma','SPEAD11','PcSaft','tcPR','GCESD','GcEsdTb','TffSPEAD','GcPcSaft','GcPcSaft(Tb)','PRLorraine','ESD2'/
	integer nParInert(nModels),nParAssoc(nModels),nParPolar(nModels) ! the # of pure compound parameters for the ith iEosOpt
	!               1 2 3 4 5  6 7 8  9 10 11 12 13 14 15 16
	data nParInert /0,3,0,3,11,0,0,11,0, 3, 4, 3,11, 0, 3, 3/	! e.g. m, sigma, eps/kB. for PcSaft & ESD
	data nParAssoc /0,2,0,2, 0,0,0, 0,0, 2, 0, 2, 0, 0, 2, 2/	! GC&Tff EOS's have zero adj par's even if assoc. JRE 20210421
	data nParPolar /0,0,0,0, 0,0,0, 0,0, 2, 0, 0, 0, 0, 0, 0/	! GC&Tff EOS's have zero adj par's even if assoc. JRE 20210421
	! PR & ESD's use Tc,Pc,w. AC models have no pure pars. 
	! SPEAD: 11=A0coeff(3),A1coeff(3),A2Coeff(4),vEffNm3. Does not allow adjustment of Assoc parameters because these must be transferable. 
	! tcPR has alphaN&M, Gc&Tff models have zero adj parameters 
	QueryNparPure=nParInert(iEosOpt)+nParPolar(iEosOpt)+nParAssoc(iEosOpt) 
	!if(class=='assoc')QueryNparPure=QueryNparPure+nParAssoc(iEosOpt)		! m, sigma, eps/kB, Kad, epsHB
	return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!int QueryNparMix();    //number of adjustable parameters for the current model and mixture
	Integer Function QueryNparMix()	!NOTE: calling routine must declare "integer QueryNparMix"
	USE GlobConst, only: iEosOpt
	!USE BIPs      !ESD, SPEADMD,
	!integer QueryNparMix 
	QueryNparMix=1 ! NparMix=1 for ESD, TPT,
	if(iEosOpt==10)then		!PcSaft(Gross)
		!! kij(i,j)          binary correction parameter acting on dispersion term, [/]\n
		!! lij(i,j)          asymmetric binary correction para. (Tang and Gross, 2010), [/]\n
		!! kij_assoc(i,j)    correction para. for association energy. [/]\n
		QueryNparMix=1		!JRE 20210531: just optimize Kij for now. We can try more parameters later.
	elseif(iEosOpt==11)then	!tcPRq: kij (for now, could include betaij later). JRE 20210531
		QueryNparMix=1
	elseif(iEosOpt== 3)then !PR76+WSmixing
		QueryNparMix=3 ! kij,tau12,tau21
	elseif(iEosOpt==17)then		! PrLorraine
		QueryNparMix=2	 !  age0(1,2) and age0(2,1)
	endif
	return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QueryParPure(iComp,iParm,value,iErr)
	USE GlobConst      ! iEosOpt,
	USE ESDParms
	DoublePrecision value
	Integer iComp,iParm,iErr
	iErr=0
	if(isEsd)then
		Call QueryParPureEsd(iComp,iParm,value,iErr) 
	elseif(isTpt)then
		Call QueryParPureSpead(iComp,iParm,value,iErr) 
	elseif(isPcSaft)then
		Call QueryParPurePcSaft(iComp,iParm,value,iErr) 
	elseif(iEosOpt==11)then
		Call QueryParPurePrTc(iComp,iParm,value,iErr)
	else
		iErr=1 
	endif
	return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPure(iComp,iParm,value,iErr)
	USE GlobConst      ! iEosOpt,
	USE ESDParms
	DoublePrecision value
	Integer iComp,iParm,iErr
	iErr=0
	if(isEsd)then
		Call SetParPureEsd(iComp,iParm,value,iErr) 
	elseif(isTpt)then
		Call SetParPureSpead(iComp,iParm,value,iErr) 
	elseif(isPcSaft)then
		Call SetParPurePcSaft(iComp,iParm,value,iErr) 
	elseif(iEosOpt==11)then
		Call SetParPurePrTc(iComp,iParm,value,iErr) 
	endif
	return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine QueryParMix(iParm,value,iErr)
	USE GlobConst
	USE BIPs      !ESD, SPEADMD,
	DoublePrecision value
	iErr=0
	if(iParm > 1)then
		if(isTpt)iErr=1
		if(isEsd)iErr=1
		if(iErr)return
	endif
	if(iEosOpt==10)then
		call QueryParMixPcSaft(iParm,value,iErr) ! PcSaft uses the name Kij() in Module pcsaft_pure_and_binary_parameters, so it needs to be separate from Module BIPs. 
		return
	endif
	if(iEosOpt==17)then
		call QueryParMixPrLorraine(iParm,value,iErr) ! PcSaft uses the name Kij() in Module pcsaft_pure_and_binary_parameters, so it needs to be separate from Module BIPs. 
		return
	endif
	if(iParm==1)then
		if(iEosOpt==7)then
			value=xsTau(1,2)
		else
			value=Kij(1,2)
		endif
		return
	elseif(iParm==2)then
		if(iEosOpt==7)then		!NRTL
		elseif(iEosOpt==11)then	!tcPR, kij and betaij
			value=Lij(1,2)
		elseif(iEosOpt== 3)then !PR76+WSmixing
			value=xsTau(1,2)  ! NOTE: tau12 =/= tau21
		endif
	elseif(iParm==3)then
		if(iEosOpt==3)then		!NRTL
			value=xsTau(1,2)  ! NOTE: tau12 =/= tau21
		else
			iErr=2
		endif
	endif
	return
	end
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParMix(iParm,value,iErr)
!    1     2      3      4      5        6         7         8          9         10      11     12       13         14          15           16
!	'PR','Esd96','PRWS','Esd','SPEAD','FloryWert','NRTL','SpeadGamma','SPEAD11','PcSaft','tcPR','GCESD','GcEsdTb','TffSPEAD','GcPcSaft','GcPcSaft(Tb)'/
!   Diky   12                   6                                                 25      22     23       24         18          26			   2
	USE BIPs      !ESD, SPEADMD, PR, PRtc,
	USE GlobConst 
	DoublePrecision value
	data initCall/1/
	if(initCall)then
		initCall=0
		Kij=0
		KTij=0
		Hij=0
		HTij=0
		xsTau=0
		xsTauT=0
		xsAlpha=0
		if(iEosOpt==2)then
			if(ID(1)==1921 .or. ID(2)==1921)Kij(1,2)=0.2d0
			Kij(2,1)=Kij(1,2)
		endif
	endif
	iErr=0
	if(iParm > 1)then
		if(isTpt)iErr=1
		if(isEsd)iErr=1
		if(iErr)return
	endif
	if(iEosOpt==10)then
		call SetParMixPcSaft(iParm,value,iErr) ! PcSaft uses the name Kij() in Module pcsaft_pure_and_binary_parameters, so it needs to be separate from Module BIPs. 
		return
	endif
	if(iEosOpt==17)then
		call SetParMixPrLorraine(iParm,value,iErr) ! PcSaft uses the name Kij() in Module pcsaft_pure_and_binary_parameters, so it needs to be separate from Module BIPs. 
		return
	endif
	if(iParm==1)then
		if(iEosOpt==7)then
			xsTau(1,2)=value
		else
			Kij(1,2)=value
			Kij(2,1)=value
		endif
		return
	elseif(iParm==2)then
		if(iEosOpt==7)then		!NRTL
		elseif(iEosOpt==11)then	!tcPR, kij and betaij
			Lij(1,2)=value
			Lij(2,1)=value
		elseif(iEosOpt== 3)then !PR76+WSmixing
			xsTau(1,2)=value  ! NOTE: tau12 =/= tau21
		endif
	elseif(iParm==3)then
		if(iEosOpt==3)then		!NRTL
			xsTau(1,2)=value  ! NOTE: tau12 =/= tau21
		else
			iErr=2
		endif
	endif
	return
end	! SetParMix

