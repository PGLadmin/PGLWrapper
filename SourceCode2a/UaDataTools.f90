
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	USE GlobConst, only:NMX
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
	USE GlobConst, Only:DEBUG,masterDir,LOUD 
	! Purpose: Read Tm(K) and Hfus(J/mol) from input file.
	Implicit DoublePrecision(a-h,o-z)
	Character*77 HfusFile,FN !,inFile77
    FN="Hfus.txt"
	HfusFile=TRIM(masterDir)//'\input\'//TRIM(FN)
	if(DEBUG)HfusFile='c:\spead\calceos\input\'//TRIM(FN)
	if(LOUD)print*,'HfusFile=',TRIM(HfusFile)
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
				write(*,*)'idCas,ID,Hfus=',idCas,idDippr,HfusTde,HfusDip,Hfus
				!pause 'GetTmHfus: check HfusTde,HfusDip,Hfus again'
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
	if(LOUD)write(*,*)TRIM(msg) !-- causes crash???
!	pause
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdDipprLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst, only:ID,nmx,LOUD
	USE CritParmsDb, only:idCasDb,idNum,CrIndex,nDeckDb
    integer ier,	idCas(nmx)
	character*77 errMsg(0:11),errMsgPas
	errMsg(0)='No Problem'
	errMsg(1)='IdDipprLookup Error: at least one id not found'
	ier=0
	!print*,'IdLookup: nDeck=',nDeck
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
		if(LOUD.and.ier>0)write(*,*)'Did not find idCas(i).i,id=',iComp,id(i)
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
!	inFile=TRIM(masterDir)//'\input\ParmsCrit.dta' ! // is the concatenation operator
	inFile=TRIM(masterDir)//'\input\ParmsCrit.txt' ! // is the concatenation operator
	IF(DEBUG)inFile='c:\SPEAD\CalcEos\input\ParmsCrit.txt' 
	if(LOUD)print*,'LoadCritParmsDb: CritFile=',TRIM(inFile)
!	OPEN(40,FILE=inFile,FORM='BINARY')
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
!	READ(40,ERR=861)NDECK1
	READ(40,*,ERR=861)NDECK1
	!if(ndeck1.gt.ndb)pause 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
!		READ (40,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
		READ (40,'(a188)',ioStat=ioErr)dumString
		READ (dumString,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,idCasDb(I),tCode,pCode,vCode,form,NAMED(I)
		CrIndex(IDNUM(i))=i
		if(ioErr.and.LOUD)print*,'GetCrit: error reading ParmsCrit.txt. line=',i
		if(i==1.and.LOUD)write(*,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) !&
		!	,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
!	close(61)
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.3,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	CLOSE(40)
	!if((nDeck1).gt.ndb)pause 'GetCrit: too much data in file'
	if(LOUD)print*,'LoadCritParmsDb: So far so good! ParmsCrit.txt is loaded. Trying ParmsCrAdd.'
	IF(DEBUG)THEN
		OPEN(40,FILE='c:\spead\calceos\input\ParmsCrAdd.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	ENDIF
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
	if(LOUD)print*,'LoadCritParmsDb: Success! DB is loaded.'
	!if((nDeck).gt.ndb)pause 'GetCrit: too much data in ParmsCrAdd file'
	!Write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

	RETURN
861	continue
	!write(*,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
	return                      
862	continue
	!write(*,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	!write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
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
!	inFile=TRIM(masterDir)//'\input\ParmsCrit.dta' ! // is the concatenation operator
	inFile=TRIM(masterDir)//'\input\ParmsCrit.txt' ! // is the concatenation operator
	IF(DEBUG)inFile='c:\SPEAD\CalcEos\input\ParmsCrit.txt' 
	if(LOUD)write(*,'(2a)')' CritFile=',TRIM(inFile)
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
!	READ(40,ERR=861)NDECK1
	READ(40,*,ERR=861)NDECK1
	!if(ndeck1.gt.ndb)pause 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
!		READ (40,ERR=861)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
		READ (40,'(a188)',ioStat=ioErr)dumString
		READ (dumString,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas
		READ (dumString,'(a126,3a4,a12,a30)',ioStat=ioErr)readText,tCode,pCode,vCode,form,NAMED(I)
		if(ioErr.and.LOUD)print*,'GetCrit: error reading ParmsCrit.txt. line=',i
		if(i==1.and.LOUD)write(*,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) !&
		!	,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
!	close(61)
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1)==0)THEN
		DO I=1,NDECK1,3
			IF(LOUD)WRITE(*,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
			!IF(((I-1)/45)*45.EQ.(I-1))PAUSE
		enddo
		RETURN
	END IF
	CLOSE(40)
	!if((nDeck1).gt.ndb)pause 'GetCrit: too much data in file'

	IF(DEBUG)THEN
		OPEN(40,FILE='c:\spead\calceos\input\ParmsCrAdd.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	ENDIF
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
	!if((nDeck).gt.ndb)pause 'GetCrit: too much data in ParmsCrAdd file'
	!Write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

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
			!write(*,*)'Error in GetCrit: not found for iComp=',iComp
			!pause
		endif
		if(LOUD)write(*,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp),TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	!write(*,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
	return                      
862	continue
	!write(*,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	!write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
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
    USE GlobConst, only: masterDir,DEBUG,LOUD
    USE VpDb ! Stores all the coefficients.
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    character*255 inFile
	iErrCode=0
    inFile=TRIM(masterDir)//'\input\CoeffsVp2a.TXT'
	IF(DEBUG)INFILE='c:\spead\calceos\input\CoeffsVp2a.TXT'
    OPEN(662,FILE=inFile,ioStat=ioErr)
    if(ioErr.and.LOUD)pause 'GetVpDb: error opening CoeffsVp2a.txt'
    !open(662,file='junk.txt')
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	!ndeck1=1974
	READ(662,*,ioStat=ioErr)NDECK1 ! ioStat= -1,endofFile; -2,endOfLine
    if(ioErr)print*,'ioErr,nDeck1,inFile',ioErr,nDeck1,TRIM(inFile)
    if(ioErr)pause 'GetVpDb: Error reading NDECK1'
    
	!if(ndeck1.gt.ndb)pause 'GetVp: more data in file than allocated'
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
		READ(662,*,ioStat=ioErr)IDNUM(I),rMINTD(I) ,VALMIND(I) ,rMAXTD(I),VALMAXD(I),AVGDEVD(I),NUMCOEFFD(I) ,(vpCoeffsd(i,iCoeff),iCoeff=1,5)  
	    if(ioErr.and.LOUD)print*, 'inFile=',TRIM(inFile)
	    if(ioErr.and.LOUD)print*, 'ioErr,line,id,vpCoeffs=',ioErr,I,IDNUM(I),(vpCoeffsd(i,iCoeff),iCoeff=1,5)
	    if(ioErr.and.LOUD)pause 'GetVp: error reading CoeffsVp2a.txt'
		indexVpDb(IDNUM(I))=I ! vpCoeffs(iComp,iCoeff)=vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)
    enddo
    if(LOUD)print*,'GetVpDb: Success! USE VpDb for vpCoeffsd(indexVpDb(idDippr(iComp),iCoeff)'
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
	USE GlobConst, only:DEBUG,NMX
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
			!write(*,*)'Error in GetCrit: not found for iComp=',iComp
			!pause
		endif
		!write(*,'(1x,i5,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
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
	!write(*,*)'# of bips in database =',numBips

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
	if(LOUD)write(*,'(a33,a22)')'GetBip error - error reading ',bipFile
	if(LOUD)write(*,*)'numBips,item',numBips,item
	if(LOUD)write(*,*)'Data read as:'
	if(LOUD)write(*,'(i4,7F8.4)')idBinarY(item),KIJDB(item),KTIJDB(item),wsTAUij(item),wsTAUji(item),wsTauTij(item),wsTauTji(item),alphaDB(item)
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
	IF(DEBUG)THEN
		OPEN(50,FILE='c:\spead\calceos\input\idTrcDipCas.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT' ! // is the concatenation operator
		OPEN(50,FILE=inFile)
	ENDIF
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
		if(LOUD)write(*,*)'Did not find idcc(1)=',idcc(1)
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
			if(LOUD)write(*,*)'Did not find idcc(i).i,idcc=',iComp,idcc(1)
			errMsgPas=errMsg(ier)
			return
		endif
	enddo
	close(50)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine IdCasLookup(NC,idCas,ier,errMsgPas)
    USE GlobConst, only:ID,idTrc,nmx,class,DEBUG,LOUD,masterDir
	parameter(maxDb=3000)
	character*77 errMsg(0:11),errMsgPas,dumString
	character*251 inFile
	character*28 dumName
	character*5 classdb(maxDb)
	integer idCas(nmx)
	integer idCasDb(maxDb),idDb(maxDb),idTrcDb(maxDb)
    !OPEN(50,FILE='c:\spead\idTrcDipCas.TXT')
	IF(DEBUG)THEN
		inFile='c:\spead\calceos\input\idTrcDipCas.TXT'
	ELSE
		inFile=TRIM(masterDir)//'\input\idTrcDipCas.TXT' ! // is the concatenation operator
	ENDIF
	if(LOUD)print*,'idCasLookup: inFile=',TRIM(inFile)
	OPEN(50,FILE=inFile)
	errMsg(0)='No Problem'
	errMsg(1)='IdCasLookup Error: at least one id not found'
	ier=0
	iGotIt=0
	read(50,*)nDeck
	!print*,'IdLookup: nDeck=',nDeck
	do i=1,nDeck
		read(50,'(a77)')dumString
		read(dumString,'(2i6,i12,1x,a28,a5)',ioStat=ioErr)idTrcDb(i),idDb(i),idCasDb(i),dumName,classdb(i)  
		if(ioErr.and.LOUD)print*,'idCas ioErr. i,idDippr=',i,idDb(i)
		if( id(1)==idDb(i) )then
			iGotIt=1
			idCas(1)=idCasDb(i)
			idTrc(1)=idTrcDb(i)
			class(1)=classdb(i)
			if(LOUD)print*,TRIM(dumString)
			if(LOUD)print*,'IdCasLookup:id1,name,class1:',id(1),TRIM(dumName),class(1)
		endif
	enddo
	if(iGotIt==0)then
		ier=1
		!pause 'IdCasLookup: idCas not found for comp1'
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
				if(LOUD)print*,'IdCasLookup:idi,classi:',id(i),class(i)
				iGotit=1
				exit
			endif
		enddo
		if(iGotIt==0)then 
			ier=1 
			if(LOUD)write(*,*)'IdCasLookup:Did not find idCas(i).i,id=',iComp,id(i)
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
	if(LOUD)write(*,*)'PrintErrFlash:iErrCode=',iErrCode
	IF(iErrCode.EQ.1)WRITE(*,*)'LLeFL:Initial guess REQUIRED.'
	IF(iErrCode.EQ.2)WRITE(*,*)'Flash:ERROR ALL UPPER PHASE'
	IF(iErrCode.EQ.3)WRITE(*,*)'Flash:ERROR ALL LOWER PHASE'
	IF(iErrCode.EQ.4)WRITE(*,*)'Flash:LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.5)WRITE(*,*)'Flash:VAPOR  FUGACITY COEFFICIENT CALCULATION HAD ERROR'
	IF(iErrCode.EQ.6)WRITE(*,*)'Flash:ITMAX EXCEEDED		 '
	IF(iErrCode.EQ.7)WRITE(*,*)' Flash:xFrac(i) < 0 for some i'
	IF(iErrCode.EQ.11)WRITE(*,*)' zV~zL - Flash: warning'
	!IF(iErrCode.EQ.21)	!NOTE: not relevant for polymers because vLiq >> 1 and p >> pSat
	!& WRITE(*,*)' zL > 0.33 - Flash:warning '
	IF(iErrCode.EQ.22)WRITE(*,*)' zV < 0.33 - Flash:warning'
	IF(iErrCode.EQ.23)WRITE(*,*)' Flash warning:gEx < 0 on at least one LLE iteration.'
	!C
	!pause 'PrintErrFlash: Check errors'
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
	inFile=TRIM(masterDir)//'\input\ParmsCrit.txt' ! // is the concatenation operator
	IF(DEBUG)inFile='c:\SPEAD\CalcEos\input\ParmsCrit.txt' 
	!print*,'CritFile=',TRIM(inFile)
	!OPEN(40,FILE=inFile)
	OPEN(40,FILE=inFile)
	!C	open(61,FILE='ParmsCrit.dta',FORM='BINARY')
	I=0
	READ(40,*,ERR=861)NDECK1
	print*,'nDeck1=',nDeck1
	print*,'ID(1)=',ID(1)
	!if(ndeck1.gt.ndb)pause 'GetCrit: more data in file than allocated'
!	open(61,file='c:\spead\CalcEos\ParmsCrit.txt')
	DO I=1,NDECK1
		!NOTE: Can NOT read dumString here b/c unformatted read from dumString is not allowed.
		!if(i.eq.691)pause
		READ (40,'(a188)',ioStat=ioErr)readText
		READ (readText,*,ioStat=ioErr)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
!		write(61,102)IDNUM(I),TCD(I),PCD(I),ZCD(I),ACEND(I) &
!			,rMwD(i),solParmD(i),vLiqD(i),tBoil,tMelt,hFor,gFor,iCas,tCode,pCode,vCode,form,NAMED(I)
	enddo
!	close(61)
100	FORMAT(I5,2X,A20,2X,F7.2,2X,F8.4,2X,F7.4,2X,F7.4)
101	format(i5,11f10.0,i10,3a4,1x,a12,a33)
102	format(i5,11f10.0,i10,3a4,1x,a12,a33)
599	FORMAT(3(1X,i5,1X,A20))
	!Tc          PcMpa     Zc       acen      MW      solParm     vL        NBP        MP  hfor(kJ/mol)  gfor       Cas#  tC  pC  vC  FORM        Name 
	IF(ID(1).EQ.0)THEN
		DO I=1,NDECK1,3
			if(LOUD)write(*,599)IDNUM(I),NAMED(I),IDNUM(I+1),NAMED(I+1),IDNUM(I+2),NAMED(I+2)
			!IF(((I-1)/45)*45.EQ.(I-1))PAUSE
		enddo
		RETURN
	ENDIF
	CLOSE(40)
	!if((nDeck).gt.ndb)pause 'GetCrit: too much data in file'

	IF(DEBUG)THEN
		OPEN(40,FILE='c:\spead\calceos\input\ParmsCrAdd.TXT')
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsCrAdd.TXT' ! // is the concatenation operator
		OPEN(40,FILE=inFile)
	ENDIF
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
	!if((nDeck).gt.ndb)pause 'GetCrit: too much data in ParmsCrAdd file'
	if(LOUD)write(*,*)' ID    Name                  Tc(K)   Pc(MPa)    acen      Zc'

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
			!write(*,*)'Error in GetCrit: not found for iComp=',iComp
			!pause
		endif
		if(LOUD)write(*,'(1x,i11,1x,a,1x,f7.0,3(1x,f8.2))')id(iComp),NAME(iComp) &
		,TC(iComp),PC(iComp),acen(iComp),ZC(iComp)
	enddo
	RETURN
861	continue
	if(LOUD)write(*,*)'GetCrit error - error reading ParmsCrit.txt. Path? Debug?'
	!write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
	return                      
862	continue
	if(LOUD)write(*,*)'GetCrit error - error reading ParmsCrAdd.txt. Path?'
	if(LOUD)write(*,*)'nDeckCrit,nDeckCrAdd,iCompo',NDECK1,NDECK2,I
	iErrCode=1
	!pause
	return                      
END ! subroutine GetCrit
