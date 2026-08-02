!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE EsdParms
	USE GlobConst, only:nmx
	DoublePrecision eokP(nmx),KCSTAR(nmx),DH(nmx),c(nmx),q(nmx),vx(nmx)
	DoublePrecision mShape(nmx),KadNm3(nmx),epsA_kB(nmx),epsD_kB(nmx) 
	DoublePrecision ESD2_B0,ESD2_k0,ESD2_B1,ESD2_alpha2,ESD2_K10(nmx),ESD2_K11(nmx),ESD2_K12,ESD2_qCorr(0:2),ESD2_zCorr(0:2)!ESD2 coefficients
	Integer         ND(nmx),NDS(nmx),NAS(nmx)
	LOGICAL         isMEM2
END MODULE EsdParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE EsdMEM2ParmsDb	! create a linked list for ESDMEM2 to expedite lookup in MemScedOpt calcs.
	USE CritParmsDb, only:ndb		! ndb,..
	USE EsdParms		! need ESD2 parameters at least.
	Parameter(nMemSced=118,nTau=nMemSced-12,nBeta=37,nAlpha=12,nAssoc=18)	!for MemSced
	!NOTE: ndb(ESD)=ndb(CritParms) because we must first compute "Exact" values so all possible values will be tabulated. .
	DoublePrecision eokPdb(ndb),vxDb(ndb),mShapeDb(ndb),KadNm3Db(ndb),epsA_kBdb(ndb),epsD_kBdb(ndb),ZcEsdDb(ndb),tauDb(ndb)
	DoublePrecision eTotAssoc(nAssoc) ,tau(nTau),k11Db(ndb),k10Db(ndb) ! for ESD2 													
	Integer, SAVE::         IndexEsd(99999),NDdb(ndb),NDSdb(ndb),NASdb(ndb)   ! e.g. eokP(i)=eokPdb( IndexEsd(ID(i)) )
	LOGICAL, SAVE::         isReadEsd	
	Integer idTau(nTau),idBeta(nBeta),idAlpha(nAlpha),idAssoc(nAssoc)	! for MemSced.
!	data idBeta/ 501, 502, 504, 507, 510, 518, 840, 1071, 1080, 1090, 1093, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1114,& 
!	1132, 1180, 1181, 1183, 1301, 1312, 1313, 1314, 1315, 1359, 1391, 1402, 1403, 1404, 1405, 1421, 1446, 1457, 1461, 1479,&
!	1571, 1680, 1706, 1748, 1772, 1773, 1781, 1782, 1790, 1791, 1792, 1844, 1845, 1876, 2375, 2391, 2796, 2852, 2856, 2861, 2900,&
!	6854, 9855, 9858 /
!	data idAlpha/1511, 1521, 1522,1523, 1527, 1541, 1586, 1681, 1760, 1761, 1762, 1763, 1886, 1938/
!	data idTau/3, 5, 7, 8, 11, 12, 13, 14, 15, 17, 19, 21, 23, 27, 35, 41, 43, 46, 55, 56, 64, 66, 68, 104, 105, 138, 140, 159,&
!	 160, 209, 216, 250, 501, 502, 504, 507, 510, 518, 840, 908,1071, 1080, 1090, 1093, 1101, 1102, 1103, 1104, 1105, 1106, 1107,&
!	 1108, 1109, 1114, 1132, 1180, 1181, 1183, 1252, 1301, 1312, 1313,1314, 1315, 1359, 1391, 1402, 1403, 1404, 1405, 1421,&
!	 1446, 1457, 1461, 1479, 1501, 1511, 1521, 1522, 1523, 1527, 1541, 1571, 1586, 1645,1680, 1681, 1682, 1692, 1706, 1748,&
!	 1760, 1761, 1762, 1763, 1772, 1773, 1781, 1782, 1790, 1791, 1792, 1844, 1845, 1876, 1886, 1921,1938, 2375, 2391, 2796,&
!	 2852, 2856, 2861, 2900, 6854, 9855, 9858/
	data idTau/3, 7, 8, 11, 12, 13, 14, 15, 17, 21, 23, 27, 35, 41, 43, 46, 56, 68, 104, 105, 138, 140, 159,&
	 160, 209, 216, 250, 501, 502, 504, 507, 510, 518, 840,1071, 1080, 1090, 1093, 1101, 1102, 1103, 1104, 1105, 1106, 1107,&
	 1108, 1109, 1114, 1132, 1180, 1181, 1183, 1301, 1312, 1313,1314, 1315, 1359, 1391, 1402, 1403, 1405, 1421,&
	 1446, 1457, 1461, 1479, 1501, 1511, 1521, 1522, 1523, 1527, 1541, 1571, 1586, 1645,1680, 1681, 1682, 1692, 1748,&
	 1760, 1761, 1762, 1763, 1772, 1773, 1781, 1782, 1790, 1791, 1792, 1845, 1876, 1886,1938, 2375, 2391, 2796,&
	 2852, 2856, 2861, 2900, 6854, 9855/
	data idBeta/ 501, 1080, 1090, 1093, 1101, 1102, 1104, 1105, 1108, 1109,& 
	1132, 1180, 1181, 1183, 1301, 1313, 1314, 1315, 1402, 1403, 1405, 1446, 1461,&
	1571, 1772, 1773, 1782, 1790, 1791, 1792, 1845, 1876, 2375, 2391, 2796, 2856, 2861/
	data idAlpha/1511, 1521, 1522,1523, 1527, 1541, 1681, 1760, 1761, 1762, 1763, 1886/
	data idAssoc/ 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1114, 1132, 1180, 1181, 1183, 1792, 2852, 2861, 6854 /
	data eTotAssoc/ 2620, 2620, 2620, 2620, 2620, 2620, 2620, 2620, 2620, 2620, 2620, 2620, 1250, 1400, 1100, 1100, 575, 1150/
	Contains
	SUBROUTINE LoadEsdDb(iErrCode)
	!C  PROGRAMMED BY:  JRE 9/23
	!C  Loads THE EsdMem2 PROPERTIES (mShape,eokP,Vx) into MODULE EsdMem2ParmsDb. 
    !C      Includes IndexEsd(idDippr)= "line where idDippr was found" (linked list)
	!C      This should be a faster way of loading properties, e.g. when running VLE evaluations for a large db.
	!C  INPUT
	!C    idOpt 1 if ID is set and need to lookup idCas, 2 if idCas is set and need to lookup ID.
	!C    (ID	 VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS USEd from GlobConst if idOpt=1)
	!C    (idCas VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS USEd from GlobConst if idOpt=2)
	!C  OUTPUT
	!C    mShapeDb	
	!C    eokPdb	
	!C    VxDb		
	!C    IndexEsd	e.g., mShape(icomp)=mShapeDb(IndexEsd(ID(iComp)))
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet,ID,iEosOpt,ToUpper
	USE CritParmsDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	character*251 inFile,dumString
	Character*5 upperClass
    LOGICAL LOUDER,bEsd1
    LOUDER=LOUD
    !LOUDER=.TRUE.
	!EQUIVALENCE(IndexEsd(1),CrIndex(1))   !This doesn't work. CrIndex is in module, similar to "common" prohibition. JRE20230911
	bEsd1=.TRUE.
	if(iEosOpt==23)bEsd1=.FALSE.
	iErrCode=0
	if(.NOT.isReadCrit)then
		iErrCode=13
		return
	endif
	if(bEsd1)then
		do iComp=1,nCritSet,1	! Initialize to Exact parameter value.
			iTemp=CrIndex( IDnum(iComp) )
			IndexEsd( IDnum(iComp) )=iTemp ! All compds in CritParmsDb included in EsdParmsDb, 
			upperClass=ToUpper( classDb(iTemp) )
			!i=CrIndex( IDnum(iComp) )
			i=iComp	! i=line in ParmsCrit. line=[1,nCritSet] Here, we effectively append columns of ESD parms to the columns of ParmsCrit for all.
			if(upperClass(1:2)/='AS')then
				if(iEosOpt <19)call ExactEsd1(IDnum(i),VxDb(i),c1,mShapeDb(i),eokPdb(i),ZcEsdDb(i),iErr)
				if(iEosOpt==23)call ExactEsd2(IDnum(i),VxDb(i),c1,mShapeDb(i),eokPdb(i),ZcEsdDb(i),iErr)
			endif
			if(iEosOpt==18)tauDb(i)=0
		enddo
	endif

	inFile=TRIM(PGLinputDir)//'\ParmsEsdMem2.txt'	! replace Exact values when available.
	if(iEosOpt==18)inFile=TRIM(PGLinputDir)//'\ParmsMemSced.txt'	! replace Exact values when available.
	if(iEosOpt==23)inFile=TRIM(PGLinputDir)//'\ParmsEsd2Mem2.txt'	! replace Exact values when available.
	if(LOUD)write(dumpUnit,*)'LoadEsdParmsDb: File=',TRIM(inFile)
	OPEN(40,FILE=inFile)
	READ(40,'(a251)',ioStat=ioErr)dumString
	if(ioErr/=0)write(dumpUnit,*) 'LoadEsdDb: Failed to load ',TRIM(inFile) 
	READ(dumString,*,ioStat=ioErr)NDECK1
	if(ioErr/=0)pause 'LoadEsdDb: Failed to read nDeck '
	DO iEsd=1,NDECK1
		READ (40,'(a222)',ioStat=ioErr)dumString
		READ (dumString,*,ioStat=ioErr)IdEsd
		i=CrIndex( IdEsd )	!Every comp in critParms has EsdParms assigned.
		IndexEsd(IdEsd)=i
        if(i < 1 .or. i > nCritSet)then
            IF(LOUDER)write(dumpUnit,*)' ExactEsd: missing from ParmsPrTcJaubert, idEsd=',idEsd
            iErrCode=14		!probably, idEsd is not included in ParmsPrTcJaubert
            return
        endif
		if(iEosOpt==18)then
			READ (dumString,*,ioStat=ioErr)IdEsd,mShapeDb(I),eokPdb(i),VxDb(I),tauDb(i),KadNm3Db(i),&
		                                                 epsD_kBdb(i),epsA_kBdb(i),NDdb(i),NDSdb(i),NASdb(i)
		elseif(iEosOpt==23)then
			READ (dumString,*,ioStat=ioErr)IdEsd,mShapeDb(I),eokPdb(i),VxDb(I),k10Db(i),k11Db(i),KadNm3Db(i),&
		                                                 epsD_kBdb(i),epsA_kBdb(i),NDdb(i),NDSdb(i),NASdb(i)
		else
			READ (dumString,*,ioStat=ioErr)IdEsd,mShapeDb(I),eokPdb(i),VxDb(I),KadNm3Db(i),&
		                                                 epsD_kBdb(i),epsA_kBdb(i),NDdb(i),NDSdb(i),NASdb(i)
		endif
		if(ioErr.ne.0.and.LOUD)write(dumpUnit,*)'LoadEsdParmsDb: error reading ParmsEsdMEM2.txt. line=',iEsd
		if(ioErr.ne.0)goto 862
	enddo
	CLOSE(40)
	if(LOUD)write(dumpUnit,*)'LoadCritParmsDb: So far so good! ParmsCrit.txt is loaded. Skipping ParmsCrAdd.'
	if(LOUD)write(dumpUnit,*)'LoadEsdParmsDb: Success! DB is loaded.'
	isReadEsd=.TRUE.
	RETURN	!end LoadEsdDb
861	continue
	iErrCode=11
	close(40)
	return                      
862	continue
	iErrCode=12
	close(40)
	return                      
	END SUBROUTINE LoadEsdDb
END MODULE EsdMem2ParmsDb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BuildEsd2Corr(iErrCode)
	!C  PROGRAMMED BY:  JRE 2026
	!C  Builds THE ParmsEsd2Mem2. 
    !C      Includes IndexEsd(idDippr)= "line where idDippr was found" (linked list)
	!C      This is a faster way of loading properties, e.g., when running VLE evaluations for a large db.
	!C  INPUT
	!C    IF IFLAG = -1, write the avgDev for each compound.
	!C  OUTPUT
	!C    mShapeDb	
	!C    eokPdb	
	!C    VxDb		
	!C    IndexEsd	e.g., mShape(icomp)=mShapeDb(IndexEsd(ID(iComp)))
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet,ID,iEosOpt,nmx,form7103
	USE CritParmsDb
	!USE EsdMem2ParmsDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	parameter(nParms=3)
	!Character*251 inFile,dumString
	DoublePrecision	pDev(9999),ESD2parms(nParms),stdErr(nParms) !,chemPot(nmx)
	!Integer idTemp(nmx)
    LOGICAL LOUDER !,bEsd1
	EXTERNAL DevalPsatCorr
    LOUDER=LOUD
    !LOUDER=.TRUE.
	ESD2parms(1)= 1.33  	!qCorr(0)
	ESD2parms(2)= 6.33  	!qCorr(1)
	ESD2parms(3)= 2.33		!qCorr(2)
	ESD2parms(1)=19.5		!These values give about 4.5%dev. The values starting with 16.0 give ~4.0%dev.
	ESD2parms(3)=0.056d0
	ESD2parms(4)=0.078d0
	ESD2parms(2)=0.78d0
	ESD2parms(1)=16.0d0 	!B1 =17.11 for 2ss
	ESD2parms(2)= 1.77d0	!k11=0.114
	ESD2parms(3)= 0.12d0	!alpha2=0.064
	ESD2parms(4)= 0.18d0	!k10=0
	ESD2parms(1)= 60 		!These values give about 3.7%dev. The values starting with 16.0 give ~4.0%dev.
	ESD2parms(2)=  2.d0
	ESD2parms(3)=  -0.7d0
	ESD2parms(4)=  0.25d0

	open(61,file='c:\temp\ESD2out.txt')
	iFlag= -11	!signals to build the DB with no printing
	Call DevalPsatCorr(nPts,nParms,ESD2parms,Pdev,iFlag)
	pause 'BuildEsd2: check avgDev.'
	factor=1.d-1
	tol=1.d-5
	Call LmDifEz(DevalPsatCorr,nPts,nParms,ESD2parms,factor,Pdev,tol,iErrCode,stdErr)
	write(*,*)'iErrCode,factor,nPtsTot=',iErrCode,factor,nPts
	write(*,form7103)' ESD2parms:',ESD2parms
	write(*,form7103)' StdErr:   ',stdErr
	write(61,form7103)' ESD2parms:',ESD2parms
	Call DevalPsatCorr(nPts,nParms,ESD2parms,Pdev, -1)	!iFlag= -1 signals last call, so write(61...
	close(61)
	return
	end	SUBROUTINE BuildEsd2Corr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE DevalPsatCorr(nPts,nParms,ESD2parms,Pdev,iFlag)
	!C  PROGRAMMED BY:  JRE 2026
	!C  Builds THE ParmsEsd2Mem2. 
    !C      Includes IndexEsd(idDippr)= "line where idDippr was found" (linked list)
	!C      This should be a faster way of loading properties, e.g. when running VLE evaluations for a large db.
	!C         IF IFLAG = 1 CALCULATE THE FUNCTIONS AT X AND
	!C         RETURN THIS VECTOR IN FVEC. DO NOT ALTER FJAC.
	!C         IF IFLAG = 2 CALCULATE THE JACOBIAN AT X AND
	!C         RETURN THIS MATRIX IN FJAC. DO NOT ALTER FVEC.
	!C         IF IFLAG = -1, write the avgDev for each compound.
	!C  OUTPUT
	!C    Pdev - %AAD in Pvp for entire db. 	
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet,ID,iEosOpt,nmx,TcEos,PcEos,Tc,Pc,bVolCC_mol,ZcEos,form7123
	USE CritParmsDb
	USE EsdParms
	USE VpDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	!Character*251 inFile !,dumString
	DoublePrecision	pDev(9999),pDevComp(22),chemPot(nmx),ESD2parms(4)!,mShapeLocal,eokLocal,bVolLocal
	Dimension PvpDippr(9999),TvpDippr(9999),IDvpDippr(9999)
    LOGICAL LOUDER !,bEsd1
	Character*255 errMsg(22)
	errMsg(11)='DevalPsatCorr: BuildEsd2Corr returned error.'
	!External DevalPsatDb
    LOUDER=LOUD
	nComps=1	!required for BuildEsd2
	iComp =1	!required for BuildEsd2
    !LOUDER=.TRUE.
	if(ABS(iFlag)==11)then	!Build the DB
		nPtsTot=0
		do i=1,nVpDb
			!read(51,*)idTemp(1),TcLocal,PcLocal,acenLocal,TbLocal,TwuL,TwuM,TwuN,cCC_mol,zRa,TminLocal
			!read(51,*)idTemp,rMin_T,rMin_Val,rMax_T,rMax_Val,Avg_Dev,Num_Coeffs,vpA,vpB,vpC,vpD,vpE' !	Max_Dev	Max_Dev_T	Value
			!call MatchCritPt(TcLocal,PcLocal,acenLocal,mShapeLocal,eokLocal,bVolLocal)
			IDi=IdVp(i)
			if(IDi > 898)exit
			TminLocal=rMinTd(i)
			TmaxLocal=rMaxTd(i)
			if( CrIndex(IDi)==0 )cycle	!Compd is not in Jaubert DB.
			Tci=TcD( CrIndex(IDi) )
			do iTr=95,40,-5
				Tkelvin=iTr*Tci/100
				if(Tkelvin < TminLocal)exit
				if(Tkelvin > TmaxLocal)cycle
				PVPi=EXP(vpCoeffsd(i,1)+vpCoeffsd(i,2)/Tkelvin+vpCoeffsd(i,3)*LOG(Tkelvin)+vpCoeffsd(i,4)*Tkelvin**vpCoeffsd(i,5))/1.D6
				if(PVPi < 1.d-4)cycle
				nPtsTot=nPtsTot+1
				PvpDippr(nPtsTot)=PVPi
				IDvpDippr(nPtsTot)=IDi
				TvpDippr(nPtsTot)=Tkelvin
			enddo
		enddo
		write(*,*)'nPtsTot=',nPtsTot
		nPts=nPtsTot	
	endif
	!ESD2_qCorr(0)=ESD2parms(1)
	!ESD2_qCorr(1)=ESD2parms(2)
	!ESD2_qCorr(2)=ESD2parms(3)
	ESD2_B1 =ESD2parms(1)
	ESD2_K11(1)=ESD2parms(2)
	ESD2_K10(1)=ESD2parms(3)
	if(nParms>3)ESD2_alpha2=ESD2parms(4)
	!ESD2_alpha2=ESD2parms(3)
	Call BuildESD2qzCorr(ESD2_qCorr,ESD2_zCorr,iErr)
	if(iErr>9)then
		Pdev(1:nptsTot)=100
		return
	endif
	avgDev=0
	nCompsTot=0
	nPtsComp=0
	IDold=0
	iErrCount=0
	do i=1,nPtsTot
		!read(51,*)idTemp(1),TcLocal,PcLocal,acenLocal,TbLocal,TwuL,TwuM,TwuN,cCC_mol,zRa,TminLocal
		!read(51,*)idTemp,rMin_T,rMin_Val,rMax_T,rMax_Val,Avg_Dev,Num_Coeffs,vpA,vpB,vpC,vpD,vpE' !	Max_Dev	Max_Dev_T	Value
		!call MatchCritPt(TcLocal,PcLocal,acenLocal,mShapeLocal,eokLocal,bVolLocal)
		IDi=IDvpDippr(i)
		if(IDi/=IDold)then	!write results for previous comp and initialize for new comp.
			if(nPtsComp>0)avgDevComp=avgDevComp/nPtsComp
			if(nPtsComp>0.and.iFlag < 0)write(*,'(a,i5,1x,f7.2,i5)')' id(i),avgDev(i):',IDold,avgDevComp,nPtsComp
			if(nPtsComp>0.and.iFlag < 0)write(61,'(i5,1x,f7.2,i5,22f9.2)')IDold,avgDevComp,nPtsComp,pDevComp(1:nPtsComp)
			IDold=IDi
			nPtsComp=0
			avgDevComp=0
			nCompsTot=nCompsTot+1
		endif
		Tkelvin=TvpDippr(i)
		Tci = TcD( CrIndex(IDi) )
		Pci = PcD( CrIndex(IDi) )
		Wci=acenD( CrIndex(IDi) )
		Tr=Tkelvin/Tci
		Tc(1)=Tci
		Pc(1)=Pci
		TcEos(1)=Tci
		PcEos(1)=Pci
		Call ExactEsd2(IDi,vx(1),c(1),q(1),eokP(1),ZcEsd,iErr)	!vx,c,q,eokP USEd in ESDPARMS.
		if(iErr>9)then
			Pdev(i)=100
			cycle
		endif
		bVolCc_mol(1)=vx(1)
		mShape(1)=q(1)
		ZcEos(1)=ZcEsd
		Call PsatEar(Tkelvin,PsatMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
		!PsatExpt=PvpMPa(nComps,iComp,Tkelvin,iErrVp)
		nPtsComp=nPtsComp+1
		if(ierCode<10)then
			Pdev(i)=(PsatMPa-PvpDippr(i))/PvpDippr(i)*100
			if(IDi==IDold.or.nPtsComp==1)PdevComp(nPtsComp)=Pdev(i)
		else
			Pdev(i)=100
			iErrCount=iErrCount+1
		endif
		avgDev=avgDev+ABS( Pdev(i) )
		avgDevComp=avgDevComp+ABS( Pdev(i) )
	enddo
	avgDev=avgDev/nPtsTot
	write(*,form7123)'nComps,iErrCount,avgDev,parms:',nCompsTot,iErrCount,avgDev,ESD2parms(1:nParms),ESD2_qCorr(0:2)
	return
end SUBROUTINE DevalPsatCorr

SUBROUTINE BuildEsd2Db(iErrCode)
	!C  PROGRAMMED BY:  JRE 2026
	!C  Builds THE ParmsEsd2Mem2. 
    !C      Includes IndexEsd(idDippr)= "line where idDippr was found" (linked list)
	!C      This should be a faster way of loading properties, e.g. when running VLE evaluations for a large db.
	!C  INPUT
	!C    IF IFLAG = -1, write the avgDev for each compound.
	!C  OUTPUT
	!C    mShapeDb	
	!C    eokPdb	
	!C    VxDb		
	!C    IndexEsd	e.g., mShape(icomp)=mShapeDb(IndexEsd(ID(iComp)))
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet,ID,iEosOpt,nmx,form7123,ToUpper	!Must avoid Zc.
	USE CritParmsDb
	USE VpDb
	!USE EsdMem2ParmsDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	parameter(nParms=2)
	!Character*251 inFile,dumString
	DoublePrecision	pDev(9999),ESD2parms(nParms),stdErr(nParms) !,chemPot(nmx)
	!Integer idTemp(nmx)
    LOGICAL LOUDER !,bEsd1
	Character*5 classLocal,upperClass	!,ToUpper
	EXTERNAL DevalPsatComp
    LOUDER=LOUD
    !LOUDER=.TRUE.

	factor=1.d-1
	tol=1.d-0
	open(61,file='c:\temp\ESD2out.txt')
	do i=1,nVpDb
		IDi=idVp(i)
		classLocal=( classDb(CrIndex(IDi)) )
		upperClass=ToUpper( classDb(CrIndex(IDi)) )
		if ( upperClass(1:2) == 'AS' ) then	!try associating compounds.
			!Call LookupKAB(IDi)	!output parms USEd in EsdMEM2ParmsDb
			cycle	!Use RegDbEsd (RD option) for associating compounds.	
		endif 
		ESD2parms(1)=   2.d0	!k11=  2.0 for generalized correlation
		ESD2parms(2)=  -0.7d0	!k10= -0.7 for generalized correlation
		if(nParms>2)ESD2parms(3)=   60.d0	!k10= -0.7 for generalized correlation
		iFlag=IDi+10000
		Call DevalPsatComp(nPts,nParms,ESD2parms,Pdev,iFlag)	!First call sends IDi and initialization signal.
		if(iFlag > 10 .or. nPts < 3)cycle		!Signals error. 
		!Call LmDifEz(DevalPsatComp,nPts,nParms,ESD2parms,factor,Pdev,tol,iErrCode,stdErr)
		Call DevalPsatComp(nPts,nParms,ESD2parms,Pdev, -1)
		PAAD=SUM( ABS(Pdev(1:nPts)) )/nPts
		!write(61,form7123)' ID,nPts,ESD2parms:',IDi,nPts,ESD2parms,PAAD  **** Written in Deval... ****
	enddo
	close(61)
	return
END	SUBROUTINE BuildEsd2Db
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DevalPsatComp(nPts,nParms,ESD2parms,Pdev,iFlag)
	!C  PROGRAMMED BY:  JRE 2026
	!C  Builds THE ParmsEsd2Mem2. 
    !C      Includes IndexEsd(idDippr)= "line where idDippr was found" (linked list)
	!C      This should be a faster way of loading properties, e.g. when running VLE evaluations for a large db.
	!C         IF IFLAG = 1 CALCULATE THE FUNCTIONS AT X AND
	!C         RETURN THIS VECTOR IN FVEC. DO NOT ALTER FJAC.
	!C         IF IFLAG = 2 CALCULATE THE JACOBIAN AT X AND
	!C         RETURN THIS MATRIX IN FJAC. DO NOT ALTER FVEC.
	!C         IF IFLAG = -1, write the avgDev for each compound.
	!C  OUTPUT
	!C    Pdev - %AAD in Pvp for entire db. 	
	USE GlobConst, ONLY:LOUD,dumpUnit,zeroTol,PGLinputDir,nCritSet,ID,iEosOpt,nmx,TcEos,PcEos,Tc,Pc,bVolCC_mol,ZcEos,form7123
	USE CritParmsDb
	USE EsdParms
	USE VpDb
	IMPLICIT DoublePrecision(A-H,O-Z)
	!Character*251 inFile !,dumString
	DoublePrecision	pDev(22),chemPot(nmx),ESD2parms(nParms)!,mShapeLocal,eokLocal,bVolLocal
	Dimension Pvp(22),Tvp(22) !,IDvpDippr(9999)
    LOGICAL LOUDER !,bEsd1
	!Character*255 errMsg(22)
	!errMsg(11)='DevalPsatDb: BuildEsd2Corr returned error.'
	!External DevalPsatDb
    LOUDER=LOUD
	nComps=1	!required for BuildEsd2
	iComp =1	!required for BuildEsd2
    !LOUDER=.TRUE.
	if(iFlag > 10000)then
		IDi=iFlag-10000	!First call stores IDi and PvpDippr.
		Call PureVpTable(IDi,nPts,Tvp,Pvp,TCi,PCi,iErrVp)
		if(iErrVp > 9)then
			iFlag=11
			return
		endif
		Tc(1)=TCi
		Pc(1)=PCi
		TcEos(1)=TCi
		PcEos(1)=PCi
		iFlag=0
	endif	!Compd initialization complete
	ESD2_K11(1)=ESD2parms(1)
	ESD2_K10(1)=ESD2parms(2)
	if(nParms > 2)ESD2_B1=ESD2parms(3)
	Call ExactEsd2(IDi,vx(1),c(1),q(1),eokP(1),ZcEsd,iErr)	!vx,c,q,eokP USEd in ESDPARMS.
	if(iErr>9)then
		Pdev(1:nPts)=100
		iFlag= 13
		return
	endif
	bVolCc_mol(1)=vx(1)
	mShape(1)=q(1)
	ZcEos(1)=ZcEsd
	avgDevComp=0
	iErrCount=0
	do i=1,nPts
		Tkelvin=Tvp(i)
		Tr=Tkelvin/Tci
		Call PsatEar(Tkelvin,PsatMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
		if(ierCode<10)then
			Pdev(i)=(PsatMPa-Pvp(i))/Pvp(i)*100
		else
			Pdev(i)=100
			iErrCount=iErrCount+1
		endif
		avgDevComp=avgDevComp+ABS( Pdev(i) )
	enddo
	avgDevComp=avgDevComp/nPts
	if(nParms>2)pDev(nPts)=(ESD2_B1-60)/60*100	!try to keep B1 close to 60
	if(nPts>0.and.iFlag < 0)then
		write(* ,form7123)' id(i),nPts,avgDev:',IDi,nPts,avgDevComp
		write(61,form7123)' id(i),nPts,avgDev:',IDi,nPts,avgDevComp,q(1),vx(1),eokP(1),ESD2parms
	endif
	!write(*,form7123)'nComps,iErrCount,avgDev,parms:',nCompsTot,iErrCount,avgDev,ESD2parms(1:nParms)!,ESD2_qCorr(0:2)
	return
end SUBROUTINE DevalPsatComp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine GetEsdCas(NC,idCasPas,iErr) !ID is USEd in GlobConst
	!
	!  PURPOSE:  LOOKS UP THE ESD PARAMETERS AND STORES THEM IN USEd EsdParms
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT(to EsdParms)   qShape,eokP,VX,KadNm3,epsA_kB,epsD_kB,ND,NDS,NAS, 
	USE GlobConst   !is implied by USE ASSOC. GlobConst includes ID,Tc,Pc,Zc,acen,...
	USE EsdParms	!Note: EsdParms reads the parameters  from disk 
	USE EsdMem2ParmsDb
	USE Assoc ! For eAcceptor,eDonor,... 
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(listPool=1000)
	Character*222 bipFile,inFile,dumString !,dumString,ParmsTptFile*50 !,bipHbFile*50
	Integer iGotIt(NC),idCasPas(NC),idCasa(ndb),GetBIPs,ierCompExact(NC)  ! ndb USEd from EsdMem2ParmsDb(CritParmsDb) 
	DoublePrecision QA(ndb),KCSTA(ndb),eDonEpsK(ndb),eAccEpsK(ndb),bVolA(ndb),eokA(ndb),DHA(ndb)
	Integer IDA(ndb),NDSA(ndb),NASA(ndb),NDA(ndb)
	!doublePrecision bondRate(nmx,maxTypes)
	!iErr=1 !Parms missing and iErrExact.ne.0 for at least one component
	LOGICAL LOUDER
	LOUDER=LOUD
	ESD2_qCorr=0	!Initialize to zero to indicate not set.
	ESD2_zCorr=0
	ESD2_B0=4.6d0
	ESD2_k0=1.74d0
	ESD2_B1=19.5
	ESD2_alpha2=0.056d0
	ESD2_k10(1)=0.078d0
	ESD2_k11(1)=0.78d0
	ESD2_B1 =16.0d0		!B1 =17.11 for 2ss
	ESD2_k11(1)= 1.27d0	!k11=0.114
	ESD2_alpha2= 0.12d0	!alpha2=0.064
	ESD2_k10(1)= 0.18d0	!k10=0
	ESD2_B1 =60.0d0			!B1 =64 for 2ss
	ESD2_k11(1)= 2.d0			!k11=1.65
	ESD2_k10(1)= -0.7d0		!k10= -0.378
	ESD2_alpha2= 0.25d0		!alpha2=0.25
	ESD2_qCorr(0)= 1.1355d0		!YYE
	ESD2_qCorr(1)= 4.7120d0		!YYE
	ESD2_qCorr(2)= 0.9693d0		!YYE
	ESD2_zCorr(0)= third		!YYE
	ESD2_zCorr(1)= 0.0808d0		!YYE
	ESD2_zCorr(2)=-0.0407d0		!YYE
	!LOUDER=.TRUE.
	idCas(1:NC)=idCasPas(1:NC) ! workaround after promoting idCas to GlobConst 
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of bTPT or bPcSaft.
	bESD=.TRUE. ! USEd in GlobConst, simplifies calls in FuVtot or FUGI
    bNeedFullWertheim=.FALSE. ! USEd in Assoc, ESD applies MEM2&MEM1 by convention. 
	if(LOUDER)write(dumpUnit,*)' GetEsdCas: idCas()=',idCas(1:NC) 
	if(LOUDER)write(dumpUnit,*)' GetEsdCas: ID()=',ID(1:NC) 
	if(LOUDER)write(dumpUnit,610)' GetEsdCas: Tc()=',Tc(1:NC)
610 format(1x,a,12E12.4)
	k0=1.9d0
	B0=4
	if(iEosOpt==23 .and. ESD2_k0 > 0)then
		k0=ESD2_k0
		B0=ESD2_B0
	endif
	qFactor=B0/(B0-k0)
	etaMax=1/k0-zeroTol
	isMEM2=.FALSE.
	if(iEosOpt==4.or.iEosOpt==18 .or. iEosOpt==23)isMEM2=.TRUE.
	
	nComps=NC
	TcEos(1:nComps)=Tc(1:nComps)	!This EOS is consistent with experimental values for Tc,Pc.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)	!This is a placeholder

	if(isMEM2)then ! The MEM2 option uses the approach of loading the entire database (~1800 compounds) at outset. 
		if(.NOT.isReadEsd)Call LoadEsdDb(iErr)		! LoadEsdDb puts ESD parms into memory for all nCritSet compds.
		if(iErr > 10 .and. LOUD)write(dumpUnit,*)' GetEsdCas: iErr(LoadEsdDb)=',iErr
		if(iErr > 10)return
		do J=1,NC
			i=IndexEsd( ID(J) ) 
			q(J)=mShapeDb(I)
			c(J)=1+(q(J)-1)/qFactor ! 4/(4-1.9)=1.90476
			eokP(J)=eokPdb(I)
			vx(J)=VxDb(I)
			bVolCc_mol(J)=vx(J)		! need bVol generally for every EOS. Included in GlobConst
			KadNm3(j)=KadNm3Db(i)	! KcStar[=] nm^3 since 2021.
			nTypes(j)=1				! ESD is always 1.
			nDegree(j,1)=NdDb(i)
			nAcceptors(j,1)=nAsDb(i)
			nDonors(j,1)=nDsDb(i)
			eAcceptorKcal_mol(j,1)=epsA_kBdb(I)/1000*(RgasCal)	! cf. Table 6.1 of PGL6ed
			eDonorKcal_mol(j,1)=epsD_kBdb(I)/1000*(RgasCal)
			bondVolNm3(j,1)=KadNm3Db(i) 
        	tKmin(J)=0.4d0*Tc(J) ! this is the general rule for ESD.
			ZcEos(J)=ZcEsdDb(i)
			tau(j)=tauDb(i)
			ESD2_K10(j)=k10Db(i)
			ESD2_K11(j)=k11Db(i)
			ND(j)=NdDb(i)
			NAs(j)=nAsDb(i)
			NDs(j)=nDsDb(i)
		enddo
	else ! read data from hard drive
		inFile=TRIM(PGLinputDir)//'\ParmsEsd96.TXT'
		if(iEosOpt==12)inFile=TRIM(PGLinputDir)//'\ParmsEsdEmami.txt' ! // is the concatenation operator
		if(iEosOpt==13)inFile=TRIM(PGLinputDir)//'\ParmsEsdEmamiTb.txt' ! // is the concatenation operator
		OPEN(31,FILE=inFile)
		if(LOUDER)write(dumpUnit,*)'GetEsdCas:inFile=',TRIM(inFile)
		if(LOUDER)write(dumpUnit,*) 'Check the ESD parms file location.'
		READ(31,*)nDeck
		if(LOUDER)write(dumpUnit,602)' GetEsdCas: nDeck=',nDeck
		do i=1,nDeck
			READ(31,'(a222)')dumString
			!if(isMEM2)then
			!	READ(dumString,*,ioStat=ioErr)IDA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I)&
			!	                                                                            ,NDA(i),NDSA(I),NASA(I),idCasa(i)
				!if(LOUDER)write(dumpUnit,602)'From inFile:',IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I)
			!else ! iEosOpt=2,12,13 all use ESD96.
				READ(dumString,*,ioStat=ioErr)IDA(I),CAi,QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I) ,idCasa(i)
				!READ(dumString,*,ioStat=ioErr)IDA(I),QA(I) ,eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I),NDA(I),NASA(I),idCasa(i)
				if(NASA(I).ne.NDSA(I))NASA(I)=0	 ! ESD96 requires that NAS=NDS for all compounds.
				NDSA(I)=NASA(I)
				eAccEpsK(i)=DHA(I)*1000/RgasCal
				eDonEpsK(i)=eAccEpsK(i)
				!if(LOUDER)write(dumpUnit,603)'GetEsdCas: inFile~',i,IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I)
			!endif
			!write(dumpUnit,*),*)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I)  ,idCasa(i)
			if(ioErr/=0 .and. LOUDER)write(dumpUnit,'(a,a)')' GetESDCas: error reading ',TRIM(inFile),' line=',TRIM(dumString)
			if(  ( idCasa(i)==id(1) .or. idCasa(i)==id(2) ) .and. LOUDER  )write(dumpUnit,*)'Found in ParmsEsd idCas=',idCasa(i) 
		enddo !i=1,NC
		CLOSE(31)
		if(LOUDER)write(dumpUnit,*)'nDeck,id(nDeck)=',nDeck,ida(nDeck)
		if(LOUDER)write(dumpUnit,*)'nDeck,idCas()=',nDeck,(idCas(i),i=1,NC)

		!  Begin by computing corr states values.  these will be replaced if in dbase
		ierCompExact=0
		if(iEosOpt > 4)then
			ierCompExact=11 !Declare error because iEosOpt > 4 means using only GC parameters from one of ParmsEsdEmami__
		else
			if(Tc(1) < 4)call GetCritCas(NC,idCas,iErrCrit) !GetCritCas assumes ID(GlobConst)=IdCas
			call ExactEsd(NC,vx,c,q,eokP,iErrExact,ierCompExact) !iErrExact = 100+iComp if compd is assoc. Check ParmsEsd before fail.
			!mShape(1:nmx)=Q(1:nmx)
			if(iErrExact>0)then
				if(LOUDER)write(dumpUnit,*)'GetESDWarning: iErrExact=',iErrExact,' ierComp='	,ierCompExact(1:NC)
			endif
		endif

		nTypes(1:NC)=1	 !all esd versions use nTypes=1. Multifunctional molecules require SPEADMD.
		DO J=1,NC
			bVolCc_mol(j)=vx(j) !Copy ExactEsd value first. Replaced below if in dbase.
			ND(J)=0
			NDS(J)=0
			NAS(J)=0
			nDegree(J,1)=0
			nDonors(J,1)=0
			nAcceptors(J,1)=0
			eAcceptorKcal_mol(j,1)=0
			eDonorKcal_mol(j,1)=0
			bondVolNm3(j,1)=0
			iGotIt(J)=0
			DO I=1,NDECK
				IF(IDCASA(I).EQ.IDCas(J))THEN
					iGotIt(J)=1
					eokP(J)=eokA(I)
					!ND(J)=NDA(I)
					q(J)=QA(I)
					c(J)=1+(q(J)-1)*(4-1.9d0)/4 ! 4/(4-1.9)=1.90476
					vx(J)=bVolA(I)
					bVolCc_mol(J)=vx(J) !need bVol generally for every EOS. Included in GlobConst
					KCSTAR(J)=KCSTA(I) ! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
					epsA_kB(J)=eAccEpsK(I)
					epsD_kB(J)=eDonEpsK(I)
					ND(J)=NDA(i)
					DH(J)=DHA(I)*1000/RgasCal/ Tc( j)	! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
        			tKmin(J)=0.4d0*Tc(J) ! this is the general rule for ESD.
					NAS(J)=NASA(I)
					NDS(J)=NDSA(I)
					iComplex=0
					if(NAS(j)*NDS(j) > 0)iComplex=1
					if(iComplex==0 .and. iEosOpt==2)KCSTAR(j)=0 ! Disable solvation for iEosOpt==2 unless self-assoc.
					KadNm3(j)=KcStar(j) !KcStar[=] nm^3 since 2021.
					nDegree(j,1)=NDA(i)
					nAcceptors(j,1)=NASA(i)
					nDonors(j,1)=NDSA(i)
					eAcceptorKcal_mol(j,1)=eAccEpsK(I)/1000*(RgasCal)	! cf. Table 6.1 of PGL6ed
					eDonorKcal_mol(j,1)=eDonEpsK(I)/1000*(RgasCal)
					bondVolNm3(j,1)=KcSta(i) !*bVolCc_mol(j)/avoNum
					if(LOUDER)write(dumpUnit,602)'GetEsdCAS:iGotIt! id,bVol,eAcc=',ID(j),vx(J),eAccEpsK(I)
					if(louder)write(dumpUnit,*)'GetEsdCas:nTypes,nDeg,nAcc,nDon=',nTypes(j),nDegree(j,1),nAcceptors(j,1),nDonors(j,1)
					exit !exits this do loop, not the outer one.
				ENDIF
			enddo
			if(iGotIt(J)==0)then
				if(ierCompExact( j).ne.0 )then
					iErr=11 !Parms missing and iErrExact.ne.0 for at least one component 
					if(LOUDER)write(dumpUnit,*)'GetEsdCas:Parms missing and iErrExact.ne.0 for component = ', j
					goto 861
				else
					if(LOUD)write(dumpUnit,*)'Corr. States used for EsdParms of component=',j
				endif
			end if
		enddo
		if(iErr/=0)then
			if(LOUDER)write(dumpUnit,*) 'GetEsdCas: error for at least one compound'
			continue
		endif
        
		if(LOUDER)then
			write(dumpUnit,*)'  ID     NAME       mESD    eok      bVol    Nd   KADnm3   eDon   eAcc(kcal/mol)'
			do i=1,NC
				write(dumpUnit,606)IDCas(i),NAME(i),q(i),eokP(i),vx(i),NDS(i),NAS(i),KadNm3(i),&
																 eDonorKcal_mol(i,1),eAcceptorKcal_mol(i,1) !/1000
			enddo
        endif
    endif ! MEM2 or not MEM2
61	FORMAT(I5,2(F8.4,1X),2(F9.3,1X),I3,1X,E11.4,1X,F8.4)
601 format(1x,a,8e12.4)
602 format(1x,a,i11,8e12.4)
603 format(1x,a,2i11,8e12.4)
606	format(i9,1x,a11,f9.3,f8.2,f8.2,2i3,1x,f8.6,2f8.0)

	!note:  bips are passed back through USEd BIPs
711	bipFile=TRIM(PGLinputDir)//'\BipEsd96.txt' ! // is the concatenation operator
	if(isMEM2)bipFile=TRIM(PGLinputDir)//'\BipEsdMEM2.txt' ! // is the concatenation operator
	if(iEosOpt==18)bipFile=TRIM(PGLinputDir)//'\BipMemSced.txt' ! // is the concatenation operator
	if(iEosOpt==23)bipFile=TRIM(PGLinputDir)//'\BipEsd2MEM2.txt' ! // is the concatenation operator
	if(LOUD)write(dumpUnit,*)'GetEsdCas: bipFile=',TRIM(bipFile)
	if(NC > 1)iErrCode=GetBIPs(bipFile,ID,NC) !not necessary for pure fluids
	if(iErrCode > 10)iErr=11 ! 
    if(LOUDER)then
		write(dumpUnit,*)'GetEsdCas: bipFile=',TRIM(bipFile)
		write(dumpUnit,*)'     ',(id(j),j=1,NC)
		do i=1,NC
			write(dumpUnit,'(i5,11f8.4)')id(i),(Kij(i,j),j=1,NC)
		enddo
	    write(dumpUnit,*) 'GetEsdCas: check BIPs.'
		if(iErrCode > 0) write(dumpUnit,*)'GetEsdCas: BIPs missing for ',iErrCode-10,' binary combinations.'
    end if
	if(LOUD)write(dumpUnit,*) 'GetEsdCas: done. Returning'
	RETURN
	
861	continue
	!trap file reading errors
	if(LOUDER)write(dumpUnit,*)'GetEsd Error - error reading EsdParms.txt'
	if(LOUDER)write(dumpUnit,*)'nDeck,iCompo',NDECK,I

	!	do i=1,NC
	!		ID(i)=idStore(i)
	!	enddo

	return
end	!GetEsdCas
!------------------------------------------------------------------------------------
subroutine ExactEsd(NC,vx,c,q,eokP,iErr,ierComp)
	USE GlobConst
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD2 parameters for hydrocarbons based on the exact solution
	!  for Zc, Bc, and Yc vs. qShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999), also wikipedia.
	!parameter (nmx=55)
	DoublePrecision k1
	DoublePrecision vx(nmx),c(nmx),q(nmx),eokP(nmx)
    Integer ierComp(NC) !,ID(nmx)
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	ierComp=0
	iErr=0				  
	do i=1,NC
		isAssoc=0
		if(TRIM(class(i))=='assoc' .or. TRIM(class(i))=="Asso+")isAssoc=1
		if(isAssoc>0)then
			ierComp(i)=1
			iErr=100+i
			if(LOUDER)write(dumpUnit,*)'ExactESD: no parms for ID,class=',ID(i),TRIM(class(i))
			if(LOUDER)write(dumpUnit,*) 'ExactESD:check ID.'
			cycle
		endif	
		isHelium=0
		if( id(i)==913 .or. ID(i)==7440597)isHelium=1
		if(isHelium==1)then
			iErr=3
			if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for helium ID=',ID(i)
			ierComp(i)=3
			cycle
		endif	
		isH2=0
		if( id(i)==902 .or. ID(i)==133740 .or. id(i)==925 .or. id(i)==7782390)isH2=1
		if(isH2==1)then
			iErr=2
			if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for H2 or D2 ID=',ID(i)
			ierComp(i)=2
			cycle
        endif	

		k1=1.7745d0
		k2=1.0617d0
        Zm=9.5d0
        cqFactor=4/(4-1.9d0)
        Wci=ACEN(i)
		cShape=1+3.535*Wci+	0.533*Wci*Wci
		qShape=1+(cShape-1)*cqFactor
		!RooTCinv=1/SQRT(cShape)
		!ZcTmp=1.d0/3.d0+RooTCinv*(.0384+RooTCinv*(-.062+RooTCinv*(0.0723-.0577*RooTCinv)))
        RootQinv=1/sqrt(qShape)
        ZcTmp=( 1+RootQinv*(.1451d0+RootQinv*(-.2046d0+RootQinv*(0.0323d0-0*RootQinv))) )/3
 		atemp=Zm*qShape*1.9d0+4*cShape*k1-k1*1.9d0
		quadB=k1*1.9d0*ZcTmp+3*atemp
		sqArg=quadB*quadB+4*atemp*(4*cShape-1.9d0)*(Zm*qShape-k1)/ZcTmp

		Bc=ZcTmp*ZcTmp*(-quadB+SQRT(sqArg))/(2*atemp*(4*cShape-1.9d0))
		Yc=ZcTmp*ZcTmp*ZcTmp/(atemp*Bc*Bc)
		rlnY1=LOG(Yc+k2)
		c(i) = cShape
		q(i) = qShape
		eokP(i)=TC(i)*rlnY1
		bVolCc_mol(i)=Rgas*Tc(i)/Pc(i)*Bc
		vx(i)=bVolCc_mol(i)
	enddo
	return
end	!subroutine ExactEsd
!------------------------------------------------------------------------------------  
subroutine ExactEsd1(ID1,vx1,c1,q1,eokP1,ZcEsd,iErr)
! ExactEsd1 gets the exact solution for EsdParms for a single compound by referencing ID1 through the linked list of MEM2 style.
	USE GlobConst
	USE CritParmsDb
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD parameters for non-associating compounds based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
	!parameter (nmx=55)
	DoublePrecision k1
    character*5 tempClass !,ToUpper
	!DoublePrecision vx(nmx),c(nmx),q(nmx),eokP(nmx)
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	ierComp=0
	iErr=0				  
	isAssoc=0
    tempClass=ToUpper(classDb(CrIndex(ID1)))
	if(tempClass(1:2)=='AS')isAssoc=1
	if(isAssoc>0)then
		iErr=11
		if(LOUDER)write(dumpUnit,*)'ExactESD: TcPcw parms for ID,class=',ID1,tempClass
		if(LOUDER)write(dumpUnit,*) 'ExactESD:check ID. if Class isAssoc, '
	endif	
	isHelium=0
	if( id1==913 .or. ID1==7440597)isHelium=1
	if(isHelium==1)then
		iErr=13
		if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for helium ID=',ID1
		return
	endif	
	isH2=0
	if( id1==902 .or. ID1==133740 .or. id1==925 .or. id1==7782390)isH2=1
	if(isH2==1)then
		iErr=12
		if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for H2 or D2 ID=',ID1
		return
	endif	

	k1=1.7745d0
	k2=1.0617d0
    Zm=9.5d0
    cqFactor=4/(4-1.9d0)
	Wci=ACEND(CrIndex(ID1))
	cShape=1+3.535d0*Wci+0.533d0*Wci*Wci
	if(cShape < zeroTol)then
		if(LOUDER)write(dumpUnit,form611) ' ExactEsd1: ID,cShape =',ID1,cShape
		iErr=13
		cShape=1
		return
	endif
	RooTCinv=1/SQRT(cShape)
	qShape=1+cqFactor*(cShape-1)
	!ZcTmp=1.d0/3.d0+RooTCinv*(.0384+RooTCinv*(-.062+RooTCinv*(0.0723-.0577*RooTCinv)))
	if(qShape < 0)then
		if(LOUD)write(dumpUnit,*)'ExactESD1: q < 0??? ID,q=',ID1,qShape
		iErr=14
		return
	endif
    RootQinv=1/sqrt(qShape)
    ZcTmp=( 1+RootQinv*(.1451d0+RootQinv*(-.2046d0+RootQinv*(0.0323d0-0*RootQinv))) )/3
 	atemp=Zm*qShape*1.9d0+4*cShape*k1-k1*1.9d0
	quadB=k1*1.9d0*ZcTmp+3*atemp
	sqArg=quadB*quadB+4*atemp*(4*cShape-1.9d0)*(Zm*qShape-k1)/ZcTmp
	if(sqArg < 0)then
		if(LOUD)write(dumpUnit,*)'ExactESD1: sqArg(Bc) < 0??? ID,q,sqArg=',ID1,qShape,sqArg
		iErr=15
		return
	endif
	Bc=ZcTmp*ZcTmp*(-quadB+SQRT(sqArg))/(2*atemp*(4*cShape-1.9d0))
	Yc=ZcTmp*ZcTmp*ZcTmp/(atemp*Bc*Bc)
	rlnY1=LOG(Yc+k2)
	c1 = cShape
	q1 = qShape
	eokP1=TcD(CrIndex(ID1))*rlnY1
	Vx1=Rgas*TcD(CrIndex(ID1))/PcD(CrIndex(ID1))*Bc
	ZcEsd=ZcTmp
	return
end	!subroutine ExactEsd1
!------------------------------------------------------------------------------------  
subroutine ExactEsd2(ID1,vx1,c1,q1,eokP1,ZcEsd,iErr)
! ExactEsd2 gets the exact solution for Esd2Parms for a single compound by referencing ID1 through the linked list of MEM2 style.
	USE GlobConst !, ONLY:
	USE CritParmsDb
	USE EsdParms
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD parameters for non-associating compounds based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
	!parameter (nmx=55)
	character*255 errMsg(22)
	!DoublePrecision k1
    character*5 tempClass !,ToUpper
	!DoublePrecision vx(nmx),c(nmx),q(nmx),eokP(nmx)
	LOGICAL LOUDER
	LOUDER=LOUD
	LOUDER=.TRUE.
	errMsg(11)='ExactESD2: NA for assoc'
	errMsg(12)='ExactESD2: NA for He. Use Esd2Parms.txt'
	errMsg(13)='ExactESD2: NA for H2. Use Esd2Parms.txt'
	errMsg(14)='ExactESD2: NA for cShape(w) < 0'
	errMsg(15)='ExactESD2: NA for qShape(w) < 0'
	errMsg(16)='ExactCritESD: sqArg(-pCub/3) < 0???	'
	errMsg(17)='ExactCritESD: sqArg(3*cubQ/cubM/cubP) < 0??? '
	errMsg(18)='ExactCritESD: ? < 0??? '
	errMsg(19)='ExactCritESD: bepsc < 0??? '
	ierComp=0
	iErr=0				  
	isAssoc=0
    tempClass=ToUpper(classDb(CrIndex(ID1)))
	if(tempClass(1:2)=='AS')isAssoc=1
	if(isAssoc>0)then
		iErr=11
		if(LOUDER)write(dumpUnit,*)'ExactESD: TcPcw parms for ID,class=',ID1,tempClass
		if(LOUDER)write(dumpUnit,*) 'ExactESD:check ID. if Class isAssoc, '
		return
	endif	
	isHelium=0
	if( id1==913 .or. ID1==7440597)isHelium=1
	if(isHelium==1)then
		iErr=13
		if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for helium ID=',ID1
		errMsgPass=errMsg(iErr)
		return
	endif	
	isH2=0
	if( id1==902 .or. ID1==133740 .or. id1==925 .or. id1==7782390)isH2=1
	if(isH2==1)then
		iErr=12
		if(LOUDER)write(dumpUnit,*)'ExactESD2: Parms not available for H2 or D2 ID=',ID1
		errMsgPass=errMsg(iErr)
		return
	endif	

	k10=ESD2_k10(1)
	k11=ESD2_k11(1)
	a2 =ESD2_alpha2
    Zm =ESD2_B1
    cqFactor=ESD2_B0/(ESD2_B0-ESD2_k0)
	Wci=ACEND(CrIndex(ID1))
	if(iEosOpt==23 .and. ESD2_qCorr(0)>0)then
		qShape=ESD2_qCorr(0)+Wci*( ESD2_qCorr(1)+ESD2_qCorr(2)*Wci )
		cShape=1+(qShape-1)/cqFactor
	else
		cShape=1+3.535d0*Wci+0.533d0*Wci*Wci
		qShape=1+cqFactor*(cShape-1)
	endif
	RooTCinv=1/SQRT(cShape)
	if(cShape < zeroTol .or. qShape < zeroTol)then
		if(LOUDER)write(dumpUnit,form611) ' ExactEsd2: ID,cShape =',ID1,cShape
		iErr=14
		cShape=1
		errMsgPass=errMsg(iErr)
		return
	endif
	!ZcTmp=1.d0/3.d0+RooTCinv*(.0384+RooTCinv*(-.062+RooTCinv*(0.0723-.0577*RooTCinv)))
	CALL ExactCritESD(qShape,bepsc,Fc,alphac,Bc,etac,ZcEsd,iErr)
	if(iErr > 0)then
		print*,errMsgPass
		errMsgPass=errMsg(iErr)
		return
	endif
	Fc=k10+k11*bepsc		!ESD2a
	!bepsc=(Fc-k10)/k11		!ESD2a
	!Fc=k10+k11*(EXP(alpha2*bepsc)-1)	!ESD2b
	!Y=(Fc-k10)/k11						!ESD2b
	!bepsc=4*LOG( 1+(Fc-k10)/k11 )   	!ESD2b
	eokP1=TcD(CrIndex(ID1))*bepsc
	Vx1=Rgas*TcD(CrIndex(ID1))/PcD(CrIndex(ID1))*Bc
	q1=qShape
	c1=cShape
	if(iErr > 0)errMsgPass=errMsg(iErr)
	return
end	!subroutine ExactEsd2
!------------------------------------------------------------------------------------  
Subroutine ExactCritESD(qShape,bepsc,Fc,alphac,Bc,etac,Zc,iErr)
	USE GlobConst, ONLY: pi,LOUD,dumpUnit,errMsgPass,iEosOpt,zeroTol
	USE ESDParms
	Implicit DoublePrecision(A-H,K,O-Z)
	character*255 errMsg(22)
	LOGICAL ESD2a
	ESD2a=.FALSE.
	errMsg(16)='ExactCritESD: sqArg(-pCub/3) < 0???	 '
	errMsg(17)='ExactCritESD: sqArg(3*cubQ/cubM/cubP) < 0???	'
	errMsg(18)='ExactCritESD: iteration on Zc failed to converge'
	errMsg(19)='ExactCritESD: argLog < 0'
	iErr=0
	B0 =ESD2_B0
	k0 =ESD2_k0
	k10=ESD2_k10(1)
	k11=ESD2_k11(1)
	alpha2 =ESD2_alpha2
    B1 =ESD2_B1
    cqFactor=B0/(B0-k0)
	cShape=1+(qShape-1)/cqFactor
	C1=B0*cShape
    RootQinv=1/SQRT(qShape)
	if(iEosOpt==23 .and. ESD2_zCorr(0)>0)then
		Zc=ESD2_zCorr(0)+RootQinv*( ESD2_zCorr(1)+ESD2_zCorr(2)*RootQinv )
	else
		Zc=( 1+RootQinv*(.1451d0+RootQinv*(-.2046d0+RootQinv*(0.0323d0-0*RootQinv))) )/3
	endif
	itMax=11
	do iter=1,itMax
		a2= 3*k0*(2-1/Zc )-C1*(1-6*Zc+9*Zc*Zc)/Zc**3
		a1= 3*k0*k0*( 4- 4/Zc+1/Zc**2 )
		a0=k0*k0*k0*( 8-12/Zc+6/Zc**2-1/Zc**3 )+C1*k0*k0*(1+Zc*(-6+9*Zc))/Zc**3
		cubP=(3*a1-a2*a2)/3
		cubQ=(2*a2**3-9*a2*a1+27*a0)/27
		if(-cubP/3 < 0)then
			if(LOUD)write(dumpUnit,*)'ExactCritESD: sqArg(-cubP/3) < 0??? '
			iErr=16
			return
		endif
		cubM=2*SQRT(-cubP/3)
		cos3theta=(3*cubQ/cubM/cubP)
		theta=ACOS( cos3theta )/3
		cubX=cubM*COS(theta+2*pi/3)
		Fc=cubX-a2/3
		Bc=(1-3*Zc)/(Fc-k0)
		etac=Bc/Zc
		!C2c=C1+( (Fc-k0)**2*(2-3*Zc)+(Fc*Fc+Fc*k0+k0*k0)*(1-3*Zc)**2 )/( 3*(1-3*Zc)*(Fc-k0) )
		if(ESD2a)then
			bepsc=(Fc-k10)/k11				!ESD2a
			alphac=bepsc*(1+alpha2*bepsc)	!ESD2a
		else
			argLog=1+(Fc-k10)/k11
			if(argLog < zeroTol)then
				if(LOUD)write(dumpUnit,*)'ExactCritESD: argLog(1+(Fc-k10)/k11) < 0??? k10,k11=',k10,k11
				iErr=19
				return
			endif
			bepsc=4*LOG( argLog )		!ESD2b
			alphac= EXP( bepsc/4 )-1	!ESD2b
		endif
		C2c=alphac*qShape*B1
		ZcDev=1+C1*etac/(1-k0*etac)-C2c*etac/(1+Fc*etac)-Zc
		if(iter==1)then
			DevOld=ZcDev
			ZcOld=Zc
			Zc=0.35d0
			cycle
		endif
		CHNG=ZcDev/(ZcDev-DevOld)*(Zc-ZcOld)
		ZcOld=Zc
		DevOld=ZcDev
		Zc=Zc-CHNG
		if(ABS(ZcDev)<1.d-3)exit
	enddo
	if(iter>=itMax)iErr=17
	if(iErr > 0)errMsgPass=errMsg(iErr)
	iCheckAcen=1
	if(iCheckAcen==1)then
		beps7=bepsc/0.7d0
		Y7=EXP( beps7/4 )-1						!ESD2b
		F7=ESD2_k10(1)+ESD2_k11(1)*Y7					!ESD2b
		alpha7=Y7								!ESD2b
		C1=ESD2_B0*cShape
		C2=ESD2_B1*alpha7*qShape
		Zliq=0.02d0/qShape	!approximate value. goes to zero quickly for large q. 
		Bvir=C1-C2			!Z = 1+C1*eta/(1-k0*eta)-B1*alpha/(1+F*eta)~0
		!0=(1+F*eta)*(1-k0*eta)+C1*eta*(1+Feta)-C2*eta*(1-k0*eta)=1-Zliq+eta*(F-k0+C1-C2)+eta^2*(C1*F-k0*F+C2*k0)
		if( (F7-k0+C1-C2)**2-4*(C1*F7-k0*F7+C2*k0)*(1-Zliq) < 0)then
			iErr=11
			errMsgPass='BuildEsd2Corr: sqArg(Zliq) < 0?'
			return
		endif
		eta7=(  -(F7-k0+C1-C2) + SQRT( (F7-k0+C1-C2)**2-4*(C1*F7-k0*F7+C2*k0)*(1-Zliq) )  )/( 2*(C1*F7-k0*F7+C2*k0) )
		AresLiq= -C1/k0*LOG(1-k0*eta7)-C2/F7*LOG(1+F7*eta7)
		etaV=0.005
		do i=1,11
			etaVnew=eta7*exp(AresLiq+Zliq-1 -2*Bvir*etaV)	!JRE: a few iterations should suffice
			if(ABS(etaVnew-etaV)/etaV < 1.D-4)exit
			etaV=etaVnew
		enddo
		B7=etaV*(1+Bvir*etaV)	!B=bP/RT=> B7/Bc = (P/0.7)/(Pc/1)=>Pr=B7/B3*0.7
		acenEsd= -LOG10(0.7d0*B7/Bc)-1
	endif
	return
end Subroutine ExactCritESD
!------------------------------------------------------------------------------------  
subroutine BuildEsd2qzCorr(qCorr,zCorr,iErr)
! ExactEsd2 gets the exact solution for Esd2Parms for a single compound by referencing ID1 through the linked list of MEM2 style.
	USE GlobConst, ONLY: errMsgPass,third,LOUD !JRE: cant use the whole thing because "Zc" collides.
	!USE CritParmsDb
	USE EsdParms
	USE OLS	!OLS_fit(Xmat, y, coeff, stdErr, r2)
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD parameters for non-associating compounds based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
	!parameter (nmx=55)
	!character*255 errMsg(22)
	DoublePrecision mSegs(10),sqrtMsegs(10),acenEsd(10),ZcEsd(10),qCorr(0:2),zCorr(0:2),stdErr(0:2)
	LOGICAL LOUDER
	!Character*10 shortForm
	!data shortForm/'(a,13f7.3)'/
	LOUDER=LOUD
	LOUDER=.FALSE.
	LOUDER=.TRUE.
	iErr=0
	qCorr=0
	k0=ESD2_k0
	cqFactor=ESD2_B0/(ESD2_B0-ESD2_k0)
	do iq=1,10
		qShape=iq
		mSegs(iq)=qShape
		sqrtMsegs(iq)=1/SQRT(qShape)
		cShape=1+(qShape-1)/cqFactor
		Call ExactCritESD(qShape,bepsc,Fc,alphac,Bc,etac,ZcEsd(iq),iErr)
		if(iErr>9)then
			return
		endif
		beps7=bepsc/0.7d0
		F7=ESD2_k10(1)+ESD2_k11(1)*beps7			!ESD2a
		alpha7=beps7*(1+ESD2_alpha2*beps7)	!ESD2a
		beps7=bepsc/0.7d0
		Y7=EXP( beps7/4 )-1						!ESD2b
		F7=ESD2_k10(1)+ESD2_k11(1)*Y7					!ESD2b
		alpha7=Y7								!ESD2b
		C1=ESD2_B0*cShape
		C2=ESD2_B1*alpha7*qShape
		Zliq=0.02d0/qShape	!approximate value. goes to zero quickly for large q. 
		Bvir=C1-C2			!Z = 1+C1*eta/(1-k0*eta)-B1*alpha/(1+F*eta)~0
		!0=(1+F*eta)*(1-k0*eta)+C1*eta*(1+Feta)-C2*eta*(1-k0*eta)=1-Zliq+eta*(F-k0+C1-C2)+eta^2*(C1*F-k0*F+C2*k0)
		if( (F7-k0+C1-C2)**2-4*(C1*F7-k0*F7+C2*k0)*(1-Zliq) < 0)then
			iErr=11
			errMsgPass='BuildEsd2Corr: sqArg(Zliq) < 0?'
			return
		endif
		eta7=(  -(F7-k0+C1-C2) + SQRT( (F7-k0+C1-C2)**2-4*(C1*F7-k0*F7+C2*k0)*(1-Zliq) )  )/( 2*(C1*F7-k0*F7+C2*k0) )
		AresLiq= -C1/k0*LOG(1-k0*eta7)-C2/F7*LOG(1+F7*eta7)
		etaV=0.005
		do i=1,11
			etaVnew=eta7*exp(AresLiq+Zliq-1 -2*Bvir*etaV)	!JRE: a few iterations should suffice
			if(ABS(etaVnew-etaV)/etaV < 1.D-4)exit
			etaV=etaVnew
		enddo
		B7=etaV*(1+Bvir*etaV)	!B=bP/RT=> B7/Bc = (P/0.7)/(Pc/1)=>Pr=B7/B3*0.7
		acenEsd(iq)= -LOG10(0.7d0*B7/Bc)-1
	enddo
	Call polyFit(acenEsd, mSegs, 2, qCorr, stdErr)
	yint=third
	Call polyFit(sqrtMsegs, ZcEsd, 2, zCorr, stdErr, yint)
	if(LOUDER)write(*,'(a,11f9.4)')' q =',mSegs(1:10)
	if(LOUDER)write(*,'(a,11f9.4)')' w =',acenEsd(1:10)
	if(LOUDER)write(*,'(a,11f9.4)')' Zc=',ZcEsd(1:10)
	if(LOUDER)write(*,'(a,11f9.4)')' qi =',qCorr(0:2)
	if(LOUDER)write(*,'(a,11f9.4)')' zi =',zCorr(0:2)
	return
end subroutine BuildEsd2qzCorr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	FugiESD
	!   LATEST REVISION : 
	!	9/94 jre
	!	1/96 (swiTChed to chempot, ADDED POLYETHYLENE  jre)
	!	7/96 PS, PPO, PEO, PIB (ram natarajan)
	!	1/97 PS, PPO, PEO, PIB (made consistent, jre)
	!	7/06 jre ->f90, sample calcs, check <k1Yb> rule from 91, cf[2].
	!	Literature:
	!	[1] ESD, IECR, 29:1476 (1990) Note: <Yb>=<qYb>/<q> was superseded in ref[2]apx.
	!	[2] S&E, IECR, 30:524  (1991) Note:	Many typos in the apx here make it worthless. W1,W2 approach superseded in ref[3]
	!	[3] S&E, IECR, 31:2783 (1992) Note: <Yb> mix rule not clarified here, but ref[2]apx clarifies.
	!	[4] P&E, IECR, 32:3174 (1993) Note: <Yb> mix rule is wrong here.  Copied from 1990, not from the program.
	!	[5] JRE, IECR, 35:1624 (1996) Note: a typo of -1 was omitted then canceled. pdf clarifies.
	!	[6] E&N, IECR, 41:1043 (2002) 
	!	Example 1.  nC7+benzene at 458.1K,0.9MPa, Hij=0, Kij=0,
	!	id  xi  	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA 	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt
	!	17  0.6439	 47.8	2.30	280.7	0.000	0.000	0	0	1	1       0   	6.263848455	-9.701784181
	!	501 0.3561	 29.5	1.77	336.5	0.000	0.000	0	0	1	1   	0   	4.266198462	-7.308887938
	!	bMix	cshapemix	cbMix	qYbMix	k1YbMix	etaLiq	zRep	zAtt	zAssoc	
	!	41.273	2.11096846	87.128	108.587	61.8619	0.21885	3.16348	-4.11891	0	
	!	Example 2a.  MeOH+EtOH at 393.15K,0.62MPa, Hij=0, Kij=0.008,Ref[2]system with 1992 parameters and hbonding.
	!	id   xi  	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt
	!	1101 0.942	 20.359	1.1202	326.06	5.266	0.0226	1	1	.2232 .2232	-4.207667435	5.971886911	-6.246698829
	!	1102 0.058	 23.574	1.5655	270.14	4.985	0.0283	1	1	.2222 .2222	-4.407776333	7.447791846	-7.880790386
	!	bMix  	cShapeMix	cbMix 	qYbMix	k1YbMix	etaLiq 	zRep	zAtt	zAssoc
	!	20.545	1.1460274	23.5457	31.5621	44.1112	0.32134	3.78232	-2.7751	-1.9951	
	!	Example 2b.  ~MeOH+~EtOH at 393.15K,0.1MPa, Hij=0, Kij=0.002,~Ref[2]corrected for typos.
	!	id  xi  	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt
	!	1101 0.942	 10.414	2.349	197.01	5.266	0.0226	1	1	.2934 .2934	-3.028929333 7.054903125	-8.020987388
	!	1102 0.058	 15.778	2.500	206.75	4.985	0.0283	1	1	.2652 .2652	-3.528006835 9.328508527	-10.59426183
	!	bMix  	cShapeMix	cbMix 	qYbMix	k1YbMix	etaLiq 	zRep	zAtt
	!	10.728	2.357758	25.2961	22.7629	11.2779	0.24008	4.1634	-3.85983
	!	Example 3.  EtOH+H2O at 363K,0.1MPa, Hij=0, Kij=0.0323 ->etaLiq=0.3441,zLiq=0.001588,k1Yb=40.04375
	!	id  	xi	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt
	!	1102	0.5	23.574	1.1565	270.14	4.985	0.0283	1	1	0.144	0.144	-6.223	10.8003		-10.4607
	!	1921	0.5	 9.411	1.0053	427.25	5.143	0.10	1	1	0.113	0.113	-5.301	5.16093		-6.30716
	!	Example 4.  MeOH+Benzene at 331.08K,0.1MPa, Hij=0, Kij=0.0182 
	!	id  	xi	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD  	lnPhiRep	lnPhiAtt	lnPhiAssoc	
	!	1101 0.5	 20.4	1.12	326.1	5.17	0.0226	1	1	0.1878	0.1878	6.175   	-8.0857 	-3.87025
	!	501  0.5	 29.5	1.77	336.5	0.000	0.000	0	0	1	     1  	9.248   	-14.543     -0.7625
	!	bMix	cbMix	qYbMix	k1YbMix	etaLiq	zRep	zAtt	zAssoc	sqrt(alpha1)	fAssoc	kbe(1)
	!	24.95	36.04	75.85	73.77	0.3228	4.823	-4.7697	-1.0501		6.787   	0.6373	1377
	!   Note: F=x1*ralph1/(1+F*ralph1) => F=2x1*ralph1/(1+sqrt(1+4*ralph1*ralph1*x1)); (1/X-1)=ralph*F => X=1/(1+ralph*F)
	!	Example 5.  EtOH+nC7 at 343.17K,0.09633MPa, Hij=0, Kij=0.0317 
	!	id  	xi	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD  	lnPhiRep	lnPhiAtt	lnPhiAssoc	
	!	1102 0.5671	23.574	1.1565	270.14	4.985	0.0283	1	1	0.2454	0.2454	6.863214679	-9.30332785	-3.234494246
	!	17   0.4329	 47.8	2.30	280.7	0.000	0.000	0	0	1	1       0   12.27391746	-17.3688260	-0.860055613
	!	bMix	cshapemix	cbMix	qYbMix	k1YbMix	etaLiq	zRep	zAtt	zAssoc	sqrt(alpha1)	fAssoc	kbe(1)
	!	34.044	1.883551	64.1246	107.349	71.1151	0.30995	5.68055	-5.6358	-1.0410		4.701303	0.65418	998.0017

	SUBROUTINE FugiESD(tKelvin,pMPa,xFrac,NC,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr)
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	  + GlobConst{Rgas,Tc,Pc,...}
	USE GlobConst !, only:dumpUnit
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision xFrac(nmx),FUGC(nmx) !,chemPoAssoc(nmx)
    DoublePrecision tKelvin,pMPa,rhoMol_cc,zFactor,aRes,uRes,checkKij12
	LOGICAL LOUDER
	!  ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!  eokP IS THE DISPERSE ATTRACTION OVER BOLTZ k FOR PURE SPECIES
	!  KCSTAR IS THE BONDING VOLUME IN NM^3 
	!  DH IS THE BONDING ENERGY /RTC 
	!  C,Q,bVol ARE THE PURE COMPONENT EOS PARAMETERS
	!  KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!  zFactor IS PV/NoKT HERE  
	!  ier = 11 - AT LEAST ONE ERROR
	!        12 - NOT USED
	!        13 - NOT USED
	!        14 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!        15 - rho IS -VE
	!        16 - TOO MANY Z ITERATIONS
	!        17 - eta > etaMax or eta < 0
	!		111 - goldenZ instead of real Z.
	!
	DATA INITIAL/1/  !Zm=9.5 since 1996. It is 9.5 in the EL text, and Ref6. Not mentioned in ref[2,3,4,5]
	LOUDER=LOUD
	!LOUDER=.TRUE.
	if(LOUDER)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUDER)write(dumpUnit,*)'FugiEsd: Kij(1,2)= ',checkKij12
	iErr=0
	!NOTE: iErrTmin is checked in FuVtot
	sumx=SUM( xFrac(1:NC) )
	if(ABS(sumx-1) > 1e-8)then
		if(LOUDER)write(dumpUnit,*) 'FugiEsd: sumx .ne. 1'
	endif
	!INITIATE SECANT ITERATION ON rho
	bMix=SUM( xFrac(1:NC)*bVolCc_mol(1:NC) )
	Pb_RT=pMPa*bMix/(Rgas*tKelvin)
	!GUESS FOR rho
	eta=Pb_RT/1.05D0  !NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
	IF(LIQ==1 .or. LIQ==3 .or. eta>etaMax)eta=etaMax/1.15d0
	rho=eta/bMix
	if(eta > etaMax .and. LOUDER)write(dumpUnit,*)'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
	isZiter=1 ! FUGC calculations are skipped for isZiter=1
	Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	IF(iErr > 10)GOTO 86 ! let iErr=iErrF in FugiTP
	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor

	eta=etaOld/1.15D0
	IF (eta < 0 .and. LOUD) write(dumpUnit,31)LIQ
	rho=eta/bMix
	if(initial==1.and.LOUD)write(dumpUnit,*)'FugiEsd: initial eta,err',etaOld,errOld
	itMax=77
	errBesteta=1234
	do nIter=1,itMax
		Call FuEsdVtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		IF(iErr > 10)EXIT
		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial==1.and.LOUDER)write(dumpUnit,'(a,2e11.4,3f10.5)')'FugiEsd eta,Z', eta,zFactor
		if(initial==1.and.LOUDER)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd eta,CHNG,niter',eta,CHNG,niter 
		etaOld=eta
		ERROLD=ERR
		!  LIMIT THE CHANGE IN Density for liquid..
		IF(liq==1.and.DABS(CHNG/etaOld) > 0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==3.and.DABS(CHNG/etaOld) > 0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		!  Low eta must move from zero, so you can't limit its % change
		IF(liq==0.and.DABS(CHNG) > 0.02d0)CHNG=DSIGN(0.02D0,CHNG)
		IF(liq==2.and.DABS(CHNG) > 0.02d0)CHNG=DSIGN(0.02D0,CHNG)
 		eta=eta-CHNG
		if(ABS(CHNG) < errBesteta)then
			etaBest=eta
			errBesteta=ABS(CHNG)
		endif
		if(eta < 0 .or. eta > 1/1.9)eta=etaOld-DSIGN(0.1D0,CHNG)*etaOld
		IF(DABS(CHNG) < 1.D-9 .and. eta > 0)EXIT  ! Don't use CHNG/eta. Converge quickly to ideal gas if P->0, ~9 sigfigs if liquid
	enddo !nIter=1,itMax
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Iteration Concluded    !!!!!!!!!!!!!!!!!!!!!!!!!
	if(eta < 0 .or. eta > 1.9)then
		if(LOUDER)write(dumpUnit,*) 'FugiESD: eta < 0 or > 1.9 on final iteration. eta=',eta
		iErr=17
		goto 86
	endif
	!One last call to get FUGC.
	if(initial==1.and.LOUD)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd2 cnvrgd: eta,CHNG,niter',eta,CHNG,niter
	etaPass=eta
	rho=eta/bMix
	rhoMol_cc=rho 
	IF (rho < 0)THEN
        iErr=15
		if(LOUDER)write(dumpUnit,31)LIQ
		GOTO 86
	ENDIF
!  ITERATION ON rho HAS CONCLUDED.  DEPARTURES PASSED THROUGH GlobConst.
	if(eta > 0.43.and.LOUDER)write(dumpUnit,'(a,i3,F8.2,2F8.5)')' FugiESD: converged. nIter,eta,T,x1=',nIter,tKelvin,xFrac(1),eta
	if(ABS(eta-rho*bMix) > 1E-11 .and. LOUDER)write(dumpUnit,*) 'FugiESD: eta.ne.(rho*bMix)?'
	if(pMPa==0 .and. LOUDER)write(dumpUnit,*)' FugiEsd: P=0? LIQ,P=',LIQ,pMPa
	!zFactor=P/(rho*Rgas*T)  ! add this to improve precision when computing rho from Z after return.
	zStore=zFactor   
	isZiter=0
	Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErrF)
	if(ABS( (zFactor-zStore)/zStore ) > 1.D-4.and.LOUDER)write(dumpUnit,*) 'FugiESD: zFactor changed on last call???'
	if(zFactor < zeroTol)then
		iErr=11
		if(LOUDER)write(dumpUnit,*)'FugiEsd: converged Z <= 0. eta,Z=',eta,zFactor
		goto 86
	endif
	if(iErrF > 0.or.nIter > itMax-1 .or. eta < 0)then ! if iErr still > 0 on last iteration, then declare error.
		iErr=iErrF
		eta=etaBest
	endif
	FUGC(1:NC)=FUGC(1:NC)-DLOG(zFactor)	 !Must subtract ln(Z) when given Vtot as independent variable.
	if(LOUDER)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')' FugiEsd: eta,CHNG,niter,FUGC',eta,CHNG,niter,(FUGC(i),i=1,NC) 
	initial=0
	return
86	if(LOUDER)write(dumpUnit,*)' ERROR IN FugiEsd.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(NITER.GE.ITMAX)THEN
		if(LOUDER)write(dumpUnit,*)'TOO MANY Z ITERATIONS'
        iErr=16
	END IF
	IF(iErr > 10 .and.LOUDER) write(dumpUnit,*)'ERROR IN FuEsdVtot'
	if(iErr < 10)iErr=11
	initial=0
	RETURN
	END	!Subroutine FugiESD()
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	Written originally by JRE, Oct. 2019																				C
!C	Given T,V and gmol(), this routine calculates rho, zFactor, Ares_RT, Ures_RT, lnPhi 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,Ares,Ures,iErr)
	! Input:
	! isZiter:	1 if Z iteration, 0 if need fugc,Ares&Ures
	! tKelvin:	T(K)
	! vTotCc:	total volume(cm3)
	! gmol():	mole number vector
	! NC:		number of components
	! Output:
	! FUGC():	fugacity coefficients
	! zFactor:	Compressibility factor (PV/nRT)
	! aRes:		Dimensionless residual Helmholtz energy, a(T,v)/RT.
	! uRes:		Dimensionless residual Helmholtz energy, u(T,v)/RT.
	! iErr:		Warning:1=>MEM_() failed to converge; 5=T<Tmin(EOS); 
	!			Severe: 11=input nonsense; 12=xFrac(i) < 0; 13=severe error from MEM_();
	! References:
	! Elliott, J.R., S.J.Suresh, M.D.Donohue, Ind. Eng. Chem. Res., 29:1476 (1990). doi: 10.1021/ie00103a057   
	USE GlobConst !Tc,Pc,...
	USE FugiParts
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...	 
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS, B0,k0,... for ESD2
	USE BIPs
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(nmx),xFrac(nmx),FUGC(nmx),k10,k11,k2,B0,B1 !,KCSTARp(nmx),KVE(nmx)
	DoublePrecision YQVIJ(nmx,nmx),YQVI(nmx),Y(nmx,nmx),EOK(nmx,nmx)
	DoublePrecision CVI(nmx),CVIJ(nmx,nmx),QV(nmx,nmx)
    DoublePrecision voidFrac,tKelvin,zFactor !,zAssoc,aAssoc,uAssoc,rho,fAssoc
	Integer initCall
	DoublePrecision k1(nmx),bLij(nmx,nmx) 
	LOGICAL LOUDER,bEsd2
    Character*133 errMsg(0:22)
	!common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	DATA B0,k0,k10,k11,k2,B1,expE,initCall/4,1.9D0,0,1.7745D0,1.0617D0,9.5D0,1,1/
    errMsg(0)='Success!'
    errMsg(11)='FuEsdVtot: nonsense T(K),totMoles,vTotCc as input.'
    errMsg(12)='FuEsdVtot: Xi<0 for some i.'
    errMsg(5)='FuEsdVtot: T(K) < Tmin(all i).'
    errMsg(11)='FuEsdVtot: nonsense T(K),totMoles,vTotCc as input.'
    errMsg(11)='FuEsdVtot: nonsense T(K),totMoles,vTotCc as input.'
    errMsg(11)='FuEsdVtot: nonsense T(K),totMoles,vTotCc as input.'
    

	LOUDER=LOUD
	!LOUDER=.TRUE.
	stepSize=1.D-4
	iErr=0 !1=warning from AlpSolEz2, 11=input nonsense, 12=xFrac<0, 13=critical error from AlpSol, 14=voidFrac<0..
	totMoles=sum( gmol(1:NC) )
	xFrac(1:NC)=gmol(1:NC)/totMoles

	if( tKelvin	   < zeroTol .or. totMoles < zeroTol .or. vTotCc < zeroTol)then
		if(LOUD)write(dumpUnit,*)'FuEsdVtot: nonsense T(K),totMoles,vTotCc=',tKelvin,totMoles,vTotCc
		iErr=11
	endif
	iErrTmin=0
	TminTot=.01
	TrVolatile=0.1
	DO I=1,NC
		rLogPr=DLOG10( 0.0001/Pc(i) )	! ESD not for Psat<0.0001 MPa. In mixes, ensure most volatile comp has Psat > 0.0001 MPa
		aScvp=7*(1+acen(i))/3	 !SCVP: log10(Pr)=7(1+w)/3*(1-1/Tr) = a*(1-x) 
		xt1= 1-rLogPr/aScvp  ! x = 1/Tr at Psat=0.0001 MPa, first approximation of x = x1 + dx
		xt = xt1 -0.178*acen(i)*acen(i)*(1-xt1)*(1-xt1)  !Crude empirical correlation. cf. PGL6Samples.xlsx(nC19oh).
		if( xt > 2.222)xt = 2.2222	!1/2.2222 = 0.45. If 
		if( xt < 1 .and. LOUDER)write(dumpUnit,*) 'FugiEsd: TrMin > Tc???'
		TrMin = 1/xt ! = min( Tr@Psat=0.0001 or 0.45 )
		if(LOUDER)write(dumpUnit,601)' xt1,xt,TrMin',xt1,xt,TrMin
		if( tKelvin/ Tc(i) < TrMin .and. NC==1)iErrTmin=iErrTmin+1
		if( tKelvin/Tc(i) > TrVolatile) TrVolatile=tKelvin/Tc(i)  ! The largest TrVolatile is the Tr of the compd with lowest Tc. 
		if( Tc(i)*TrMin > TminTot) TminTot=Tc(i)*TrMin	 ! The largest Tmin is the weakest link. 
		if( Tc(i)*TrMin > TminTot .and. LOUDER) write(dumpUnit,*)'i,Tmin(i): ', i,Tc(i)*TrMin
		IF(xFrac(I) < 0 .and. LOUDER)write(dumpUnit,*) 'FuEsdVtot: ERROR - Xi<0'
	enddo 
	if(TrVolatile < 0.4d0)iErrTmin =2 ! it's only a problem if the most volatile compound has Tr < 0.45 or Psat < 0.0001.
	if(iErrTmin > 0) then
		iErr=5 ! warning level because functions like Vxs or Hxs might be insensitive to this issue.
		if(LOUDER)write(dumpUnit,*)'FuEsdVtot: T(K) < Tmin(all i)',tKelvin,TminTot
		!if(LOUDER) write(dumpUnit,*) 'FugiEsd: at least one compound has Tr < TrMin'
	endif
	bMix=0
	do i=1,nc
		xFrac(i)=gmol(i)/totMoles
		IF(xFrac(i) < 0)then
			if(LOUDER)write(dumpUnit,*) 'FuEsdVtotERROR IN FuEsdVtot, Xi<0'
			iErr=12
		endif
		bMix=bMix+xFrac(i)*bVolCc_mol(i)
	enddo
	if(iErr > 10)return
	rho=totMoles/ vTotCc
	eta=rho*bMix 
	if(LOUDER.and.initCall==1)write(dumpUnit,601)' FuEsdVtot: T,x1,bMix,eta=',tKelvin,xFrac(1),bMix,eta
	bEsd2=.FALSE.
	if(iEosOpt==23)then
		bEsd2=.TRUE.
		if(initCall==1)then
			!initCall=0
			B0 =ESD2_B0
			k0 =ESD2_k0
			B1 =ESD2_B1
			k10=ESD2_k10(1)
			k11=ESD2_k11(1)
			expE=0.25d0
		endif
	endif  

601 format(1x,a,8E12.4)
	YQVM=0.d0	!JRE26: This is alpha in ESD2
	VM=0.d0
	CVM=0.d0
	Cmix=0.d0
	K1YVM=0		!JRE26: This is bigF in ESD2
	iType=1	  
	DO I=1,NC
		DO J=1,NC
			bLij(i,j)=0	!JRE26: putting this as a placeholder for Sxs control at some time in future.
			kIJbip=KIJ(I,J)+KTIJ(I,J)/tKelvin 
			EOK(I,J)=DSQRT(eokP(I)*eokP(J))*(1.d0-kIJbip)
			bepsij=EOK(I,J)/tKelvin
			Y(I,J)=DEXP(expE*bepsij)-K2 
			if(bEsd2)Y(I,J)=bepsij*(1+ESD2_alpha2*bepsij)		!ESD2a
			if(bEsd2)Y(I,J)=EXP(expE*bepsij)-1					!ESD2b
			QV(I,J) = (Q(I)*bVolCc_mol(J) + Q(J)*bVolCc_mol(I)) / 2.d0
			YQVIJ(I,J)=QV(I,J)*Y(I,J)
			CVIJ(I,J) = (C(I)*bVolCc_mol(J) + C(J)*bVolCc_mol(I)) / 2.d0*(1-bLij(i,j)) 
			! e.g. (x1*c1+x2*c2)*(x1*b1+x2*b2) = x1^2*c1*b1+x1*x2*(c1*b2+c2*b1)+x2^2*b2^2
			YQVM=YQVM+YQVIJ(I,J)*xFrac(I)*xFrac(J)	   
			CVM = CVM + CVIJ(I,J)*xFrac(I)*xFrac(J)	   !note: above means <c>=sum(xi*ci) and <b>=sum(xj*bj) 
		enddo
		k1(I)=K10+k11*Y(I,I)						!JRE26: so...k1(i) supersedes bigF in ESD2 notation. 
		!if(bEsd2)k1(I)=K10+k11*EOK(i,i)/tKelvin					!k1(i) takes the place of bigF. Comment for ESD2b
		Cmix=Cmix+xFrac(i)*C(i)
		VM=VM+xFrac(I)*bVolCc_mol(I)
		K1YVM=K1YVM+xFrac(I)*K1(I)*bVolCc_mol(I) !1991 form, overwritten if applying 1990 form
	enddo
	Qmix=1+(Cmix-1)*B0/(B0-k0)
	if( ABS(vx(1) - bVolCc_mol(1)) > zeroTol .and. LOUD ) write(dumpUnit,601)' FuEsdVtot: VX.ne.bVol=',vx(1),bVolCc_mol(1)
	if(LOUD.and.k1yvm < zeroTol)write(dumpUnit,*)'FuEsdVtot: 0~k1yvm=',k1yvm 
	if(isMEM2)then
		CALL MEM2(isZiter,tKelvin,xFrac,NC,rho,zAssoc,aAssoc,uAssoc,fugAssoc,iErrMEM )!,ier)
	else
		CALL MEM1(isZiter,tKelvin,xFrac,NC,rho,zAssoc,aAssoc,uAssoc,fAssoc,fugAssoc,iErrMEM )!,aAssoc,uAssoc,rLnPhiAssoc,ier)
	endif
	if(LOUDER)write(dumpUnit,601)' FuEsdVtot: rho,zAssoc=',rho,zAssoc
	if(iErrMEM > 0 .and. LOUDER)write(dumpUnit,601)' FuEsdVtot: iErrMEM > 0. rho,zAssoc=',rho,zAssoc
	if(iErrMEM==1)iErr=1
	if(iErrMEM > 10)iErr=13
	voidFrac=1-k0*eta
	denom=voidFrac
	zRep=  B0*Cmix*eta/denom
	zAtt= -B1*YQVM*rho/(1+K1YVM*rho)
	zFactor=(1+zRep+zAtt+zAssoc)
	pMPa=zFactor*Rgas*rho*tKelvin
    if(voidFrac < 0)then
	    IF(LOUD)write(dumpUnit,*) 'FuEsdVtot:Error! (1-k0*eta) IS -VE. eta,rho=',eta,rho
		iErr=14
    endif
	if(LOUDER)write(dumpUnit,form610)' FuEsdVtot: done with isZiter=1. zAssoc,zFactor=',zAssoc,zFactor
	if(iErr > 10)return


	!aRep= -4.d0/1.9D0*DLOG(voidFrac)*CVM/VM
	aRep= -B0/k0*DLOG(voidFrac)*Cmix
	aAtt= -B1*YQVM/K1YVM*DLOG(1+K1YVM*rho)
	aRes=aRep+aAtt+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes(T,V). Important for EAR.
	if(isZiter==1)return ! don't need the rest if isZiter.
 	DO I=1,NC
	   YQVI(I)=0.D0
	   CVI(I)=0.D0
	enddo     
	BdYb_dB=0
	BdYbq_dB=0
	DO I=1,NC
		!ralph(I)=SQRT(alphAD(I,I))	 ! ralph is computed in alpsolEz.
		!k1(i)=k10+k11*Yii for ESD => dK1_dBeta=k11*(Yii+K2)
		!k1(i)=k10+k11*bepsii for ESD2 => dK1_dBeta=k11
		BdK1_dB=EOK(I,I)/tKelvin*(Y(I,I)+K2)
		if(bEsd2)BdK1_dB=EOK(I,I)/tKelvin*k11
		BdYb_dB=BdYb_dB+xFrac(I)*vx(I)*bdK1_dB  ! = Beta*d<k1Yb>/dBeta     
		DO J=1,NC
			BdY_dB=EOK(I,J)/tKelvin*(Y(I,J)+K2)
			if(bEsd2)BdY_dB=EOK(I,J)/tKelvin*(1+2*esd2_alpha2*EOK(I,J)/tKelvin)
			BdYbq_dB=BdYbq_dB+xFrac(J)*xFrac(i)*QV(i,j)*BdY_dB ! = Beta*d<Ybq>/dBeta
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*xFrac(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*xFrac(J)
 		enddo
	enddo
	UATT= -B1*YQVM*rho/(1+K1YVM*rho)*BdYb_dB/K1YVM + aAtt*(BdYbq_dB/YQVM-BdYb_dB/K1YVM) !FYI:Don't forget the aAtt*... term!

	!     CALLING ASSYMU IF NUMBER OF DONOR SITES AND ACCEPTOR SITES ARE NOT EQUAL 
	!      IF(IFLAG.EQ.2)CALL ASSYMU(ALPHAD,ALPHDA,XA,XD,X,ND,NC,QIJ,uAssoc)

	!if(LOUDER)write(dumpUnit,601)' FuEsdVtot: zAssoc,aAssoc,uAssoc=',zAssoc,aAssoc,uAssoc

	if (isZiter==0) then
		!call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,eta,fugAssoc,h_nMichelsen,dfugAssoc_dT,dfugAssoc_drho,dh_dT,dh_drho)
        if( zFactor < zeroTol)then  ! Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
		    if(LOUDER)write(dumpUnit,*) 'FuEsdVtot: zFactor.le.0 when isZiter=0. zAssoc,Z=',zAssoc,zFactor
			iErr=3	  ! warning level because another call might produce Z > 0.
			goto 86
        endif
		if(LOUDER)write(dumpUnit,'(a,f10.5)')' i,lnGamRep,lnGamAtt,lnGamBon.'
		DO I=1,NC
			!FUGREP(I)=FREP*( 2.d0*C(I)/Cmix-vx(I)/VM ) + zRep*vx(I)/VM
			fugRep(i)=aRep*( C(I)/Cmix ) + zRep*vx(I)/VM ! For pure i, FugRepi= -4ci/1.9*ln(1-1.9eta) + 4ci*eta/(1-1.9eta) 
			fugAtt(i)=aAtt*( 2*YQVI(I)/YQVM-K1(I)*vx(I)/K1YVM )+zAtt*K1(I)*vx(I)/K1YVM !91-pres form,
			!fugAssoc(i)=ND(i)*2*DLOG(XA(i,1)) + zAssoc*1.9D0*vx(i)*rho !JRE'96 Eq.43.
			FUGC(I)=FUGREP(I)+FUGATT(I)+fugAssoc(I)  ! -DLOG(Z)  Don't subtract ln(Z) when given Vtot as independent variable.
			rLnGamRep(i)=FUGREP(i)-zRep*vx(i)/VM  ! cf. Bala and Lira (2016), Eqs A6-A14. to correct from constant volume to P.
			rLnGamAtt(i)=FUGATT(i)-zAtt*vx(i)/VM
			rLnGamAssoc(i)=fugAssoc(i)-zAssoc*vx(i)/VM
			IF(LOUDER)write(dumpUnit,'(i3,f7.4,9f10.4)')i,xFrac(i),rLnGamRep(i),rLnGamAtt(i),rLnGamAssoc(i) !,ralpha(i),ralphd(i)
		ENDDO
	endif
	!I've lost faith in uAssoc for MEM2. Differentiate numerically.
	if(bTPT)then  ! disables because isESD =/= isTPT. Change to isMEM2 if you want to enable. 
		Tplus =tKelvin*(1+stepSize)
		Tminus=tKelvin*(1-stepSize)
		CALL MEM2(isZiter,Tplus ,xFrac,NC,rho,Zdum,aPlus ,uAssoc,CVI,iErrMEM )! reusing CVI here to avoid replacing fugAssoc
		CALL MEM2(isZiter,Tminus,xFrac,NC,rho,Zdum,aMinus,uAssoc,CVI,iErrMEM )!,ier)
		uAssoc= -tKelvin*(aPlus-aMinus)/(Tplus-Tminus)
	endif
	uRes=UATT+uAssoc

	aAssocPas=aAssoc
	uAssocPas=uAssoc

	if(LOUDER)write(dumpUnit,*) 'FuEsdVtot: Check results before returning.'
86	return
	end !subroutine FuEsdVtot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QueryParPureEsd(iComp,iParm,value,iErr)
	USE EsdParms      !Just for ESD
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	!DoublePrecision eokP(nmx),KCSTAR(nmx),DH(nmx),C(nmx),Q(nmx),vx(nmx)
	!Integer         ND(nmx),NDS(nmx),NAS(nmx)
	iErr=0
	if(iParm==1)then
		value=c(iComp)
	elseif(iParm==2)then
		value=vx(iComp)
	elseif(iParm==3)then
		value=eokP(iComp)
	elseif(iParm==4)then
		value=KCSTAR(iComp)
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
end	Subroutine QueryParPureEsd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPureEsd(iComp,iParm,value,iErr)
	USE EsdParms      !Just for ESD
	IMPLICIT NONE
	DoublePrecision value
	integer iComp,iParm,iErr
	!-----------------------------------------------------------------------------
	! pure component parameters
	!-----------------------------------------------------------------------------
	!DoublePrecision eokP(nmx),KCSTAR(nmx),DH(nmx),C(nmx),Q(nmx),vx(nmx)
	!Integer         ND(nmx),NDS(nmx),NAS(nmx)
	iErr=0
	if(iParm==1)then
		c(iComp)=value
		q(iComp)=1+(value-1)*1.9076D0
	elseif(iParm==2)then
		vx(iComp)=value
	elseif(iParm==3)then
		eokP(iComp)=value
	elseif(iParm==4)then
		KCSTAR(iComp)=value
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
end	Subroutine SetParPureEsd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

