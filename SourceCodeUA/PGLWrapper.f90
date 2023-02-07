	USE MSFLIB !For FILE$CURDRIVE AND GETDRIVEDIRQQ
	USE PORTLIB
	USE GlobConst
	USE CritParmsDb
	USE BIPs
	USE EsdParms
	Implicit DoublePrecision(A-H,K,O-Z)

	CHARACTER($MAXPATH) CURDIR !TO DETERMINE WHERE TO LOOK FOR PARM FILES ETC.
	CHARACTER*77 readString,Property(22) !,errMsgPas
	CHARACTER*251 outFile ,dumpFile
	!CHARACTER*3 aPhase
    DoublePrecision xFrac(NMX) !,FUGC(NMX),xPure1(NMX),xPure2(NMX) !FUGI requires mole fraction specification because it is written generally for mixtures.
    Integer  localCas(NMX),prp_id

	CURDIR = FILE$CURDRIVE
	iStat = GETDRIVEDIRQQ(CURDIR)
	masterDir=TRIM(curDir)
	PGLinputDir='c:\PGLWrapper\input'
	DEBUG=.TRUE.
	DEBUG=.FALSE. ! read input files from c:\Projects\...\input dir.
    LOUD = .TRUE.
    LOUD = .FALSE.
	CheckDLL=.TRUE.
	CheckDLL=.FALSE.
	dumpUnit=6
	if(CheckDLL)then
		LOUD=.TRUE. ! Force this, just in case forgotten above.
		dumpUnit=686
		dumpFile=TRIM(masterDir)//'\DebugDLL.txt'
		OPEN(dumpUnit,file=dumpFile)
		write(dumpUnit,*)'Debug results for PGLWrapper'	! clears old contents. Future writes will APPEND.
		!close(dumpUnit) ! Don't close dumpUnit until the end of main program so we can keep appending easily.
	endif
	outFile=TRIM(masterDir)//'\output.txt' ! // is the concatenation operator
	OPEN(52,FILE=outFile)
    open(51,FILE='input.txt')
    NC=1 !assume one component for all calculations (as of 1/1/2020 this is all we need).
    xFrac(1)=1  !   "
	read(51,'(a77)')readString
	if(LOUD)write(dumpUnit,*)'readString: ',TRIM(readString)
    read(readString,*)iEosOpt,iProperty,localCas(1)
	if(iProperty > 10) then
		read(readString,*)iEosOpt,iProperty,localCas(1),localCas(2) !,xFrac(1)
		NC=2
		if(LOUD)write(dumpUnit,*)'Mixture calculation. id1,id2',localCas(1),localCas(2)
	endif
	if(iProperty==1)then
		write(52,'(2i4,i11,1x,a,1x,a,1x,a)')iEosOpt, 1, localCas(1),TRIM(EosName(iEosOpt)),TRIM(Property(iProperty)),' = iEosOpt, iProp, idCas: T, P(expVal),pSatKPa,ierCode '
	elseif(iProperty > 1 .and. iProperty < 6)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,rho(g/cc),ierCode ',TRIM(EosName(iEosOpt)),' ',TRIM(Property(iProperty))
	elseif(iProperty>5 .and. iProperty < 11)then
		write(52,*)iEosOpt,' = iEosOpt: T, P,expVal,(uRes+zFactor-1),CpRes_R,CvRes_R,vpt,ierCode ',TRIM(EosName(iEosOpt)),' ',TRIM(Property(iProperty))
	else
		write(52,'(i4,a,5i11,1x,a,1x,a,1x,a)')iEosOpt,' = iEosOpt. localCas()= ',(localCas(i),i=1,NC),TRIM(EosName(iEosOpt)),TRIM(Property(iProperty))
	endif
	errDummy=86.86D0
	prp_id=iProperty				!iPhase=1 (liquid)
    !iProperty = 1: vapor pressure (kPa) given tKelvin
    !iProperty = 2: saturated liquid density (kg/m3) given tKelvin
    !iProperty = 13: saturated vapor density (kg/m3) given tKelvin
    !iProperty = 3: vapor density (kg/m3) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 4: liquid density (kg/m3) given tKelvin, pKPa, iPhase (=1 for liquid, 0 for vapor)
    !iProperty = 4_: Hres/RT(41),CpRes/R(42),CvResR(43),cmprsblty(44) given tKelvin, pKPa  !cmprsblty=(dP/dRho)T*(1/RT)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   READY TO CALCULATE   !!!!!!!!!!!!!!!!!!!!!!!!!!
    notDone=1
    line=0  !read statement	debugging
	if(NC==1)then
		do while(notDone>0)
			if(iProperty < 3)then
				read(51,*,ioStat=ioErr)var1	  !tKelvin
				var2=0
			else
				read(51,*,ioStat=ioErr)var1,var2 !tKelvin,pMPa
			endif
			tKelvin=var1		!This appears to be universally true.
			if(ioErr < 0)then   !this means we reached the end of the file
				notDone=0
				exit            !break the loop and return to remaining execution 
			elseif(ioErr > 0)then  
				write(52,*)'Unexpected ioErr (e.g. type mismatch?) reading input line,T=',line,tKelvin
				cycle
			endif
			line=line+1
			result = CalculateProperty(ieosOpt, localCas(1), prp_id, var1, var2, ierCode)
			write(52,'( 3(1pE11.4,1x),i4)')var1,var2,result,ierCode
		enddo
	else
		write(52,*)'PGLWrapper: Sorry, mixture calculations are not supported by PGLWrapperMain now. Try PGLDLL.'
	endif

601	FORMAT(I3,I5,9(1X,F7.4))
86	continue
    if(LOUD)write(dumpUnit,*)'Success! Results are written to output.txt'
	if(LOUD.and.CheckDLL)close(dumpUnit)
	close(52)
	close(51)
	stop
	END

