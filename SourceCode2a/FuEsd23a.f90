!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine GetEsdCas(nC,idCasPas,iErr) !ID is passed through GlobConst
	!
	!  PURPOSE:  LOOKS UP THE ESD PARAMETERS AND STORES THEM IN COMMON EsdParms
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT(to EsdParms)   qShape,eokP,VX,KadNm3,epsA_kB,epsD_kB,ND,NDS,NAS, 
	USE GlobConst !is implied by USE ASSOC. GlobConst includes ID,Tc,Pc,Zc,acen,...
	USE EsdParms
	USE Assoc ! For eAcceptor,eDonor,... 
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=3000,listPool=1000)
	Character*222 bipFile,inFile,dumString !,dumString,ParmsTptFile*50 !,bipHbFile*50
	Integer iGotIt(NC),idCasPas(NC),idCasa(ndb),GetBIPs,ierCompExact(NC)  !,idStore(NC) 
	DoublePrecision QA(ndb),KCSTA(ndb),eDonEpsK(ndb),eAccEpsK(ndb),bVolA(ndb),eokA(ndb),DHA(ndb) 
	Integer IDA(ndb),NDSA(ndb),NASA(ndb),NDA(ndb)
	!doublePrecision bondRate(nmx,maxTypes)
	!iErr=1 !Parms missing and iErrExact.ne.0 for at least one component
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	idCas(1:NC)=idCasPas(1:NC) ! workaround after promoting idCas to GlobConst 
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isTPT or isPcSaft.
	isESD=.TRUE. ! in GlobConst, simplifies calls in FuVtot or FUGI 
	etaMax=1/1.9D0-zeroTol
	isMEM2=.FALSE.
	if(iEosOpt==4)isMEM2=.TRUE.
	IF(DEBUG)then 
		inFILE='c:\Spead\CalcEos\input\ParmsEsd96.TXT'
		if(iEosOpt==12)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmami.TXT'
		if(iEosOpt==13)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmamiTb.TXT'
		if(isMEM2)inFILE='c:\Spead\CalcEos\input\ParmsEsdMEM2.TXT'
	ELSE
		inFILE=TRIM(masterDir)//'\input\ParmsEsd96.TXT'
		if(iEosOpt==12)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmami.txt' ! // is the concatenation operator
		if(iEosOpt==13)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmamiTb.txt' ! // is the concatenation operator
		if(isMEM2)inFile=TRIM(masterDir)//'\input\ParmsEsdMEM2.txt' ! // is the concatenation operator
	ENDIF
	OPEN(31,FILE=inFile)
	if(LOUDER)print*,'GetEsd:inFile=',TRIM(inFile)
	if(LOUDER)pause 'Check the ESD parms file location.'


	READ(31,*)nDeck
	if(LOUDER)write(*,602)' nDeck=',nDeck
	do i=1,nDeck
		READ(31,'(a222)')dumString
		if(isMEM2)then
			READ(dumString,*,ioStat=ioErr)IDA(I),QA(I) ,eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I),NDA(i),NDSA(I),NASA(I),idCasa(i)
			if(LOUDER)write(*,602)'From inFile:',IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I)
		else ! iEosOpt=2,11,13 all use ESD96.
			READ(dumString,*,ioStat=ioErr)IDA(I),CAi,QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I) ,idCasa(i)
			!READ(dumString,*,ioStat=ioErr)IDA(I),QA(I) ,eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I),NDA(I),NASA(I),idCasa(i)
			if(LOUDER)write(*,602)'From inFile:',IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I)
			if(NASA(I).ne.NDSA(I))NASA(I)=0	 ! ESD96 requires that NAS=NDS for all compounds.
			NDSA(I)=NASA(I)
			eAccEpsK(i)=DHA(I)*1000*RgasCal
			eDonEpsK(i)=eAccEpsK(i)
			NDSA(I)=NASA(I)
		endif
		!write(*,*)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I)  ,idCasa(i)
		if(ioErr .and. LOUDER)write(*,'(a,a)')' GetESD: error reading ',TRIM(inFile),' line=',TRIM(dumString)
		if(  ( idCasa(i)==id(1) .or. idCasa(i)==id(2) ) .and. LOUDER  )print*,'Found in ParmsEsd idCas=',idCasa(i) 
	enddo
	NDI=0
	I=0  !overrides headerless reading. Set I=1 to re-activate.
	DO while(I > 0) !reads to end of file w/o reading nDeck
		READ(31,*,ERR=861,END=50)IDA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I),NDSA(I),NASA(I),idCasa(i)
		i=i+1
		cycle
50		exit !terminate do loop
	enddo
	!nDeck=I   
61	FORMAT(I5,2(F8.4,1X),2(F9.3,1X),I3,1X,E11.4,1X,F8.4)
	CLOSE(31)
	if(LOUDER)print*,'nDeck,id(nDeck)=',nDeck,ida(nDeck)
	if(LOUDER)print*,'nDeck,idCas()=',nDeck,(idCas(i),i=1,NC)

    !  Begin by computing corr states values.  these will be replaced if in dbase
    ierCompExact=0
	if(iEosOpt > 4)then
		ierCompExact=1 ! Set error for all comps. Declare error because iEosOpt > 4 means using only GC parameters from one of ParmsEsdEmami__
	else
		if(Tc(1) < 4)call GetCritCas(nC,iErrCrit) !GetCritCas assumes ID(GlobConst)=IdCas
		call ExactEsd(nC,VX,C,Q,eokP,iErrExact,ierCompExact) !iErrExact = 100+iComp if compd is assoc or Asso+. Wait to see if Parms are in ParmsEsd before failing.
		!mShape(1:NMX)=Q(1:NMX)
		if(LOUDER.and.iErrExact)print*,'GetESDWarning: iErrExact=',iErrExact,' Checking database for iComp='	,( ierCompExact(iComp),iComp=1,NC)
	endif

	nTypes(1:NC)=1	 !all esd versions use nTypes=1. Multifunctional molecules require SPEADMD.
	DO J=1,nC
        bVolCC_mol(j)=Vx(j) !Copy ExactEsd value first. Replaced below if in dbase.
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
				Q(J)=QA(I)
				C(J)=1+(Q(J)-1)*(4-1.9d0)/4 ! 4/(4-1.9)=1.90476
				Vx(J)=bVolA(I)
				bVolCC_mol(J)=Vx(J) !need bVol generally for every EOS. Included in GlobConst
				KCSTAR(J)=KCSTA(I) ! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				epsA_kB(J)=eAccEpsK(I)
				epsD_kB(J)=eDonEpsK(I)
				ND(J)=NDA(i)
				DH(J)=DHA(I)*1000/RgasCal/ Tc( j)			! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				NAS(J)=NASA(I)
				NDS(J)=NDSA(I)
				iComplex=0
				if(NAS(j)*NDS(j) > 0)iComplex=1
				if(iComplex==0 .and. iEosOpt==2)KCSTAR(j)=0 !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				KadNm3(j)=KcStar(j) !KcStar[=] nm^3 since 2021.
				nDegree(j,1)=NDA(i)
				nAcceptors(j,1)=NASA(i)
				nDonors(j,1)=NDSA(i)
				eAcceptorKcal_mol(j,1)=eAccEpsK(I)/1000/(RgasCal)	! cf. Table 6.1 of PGL6ed
				eDonorKcal_mol(j,1)=eDonEpsK(I)/1000/(RgasCal)
				bondVolNm3(j,1)=KcSta(i) !*bVolCC_mol(j)/avoNum
				!iComplex=0
				!if(NAS(j)*NDS(j) > 0)iComplex=1
				!if(iComplex==0 .and. iEosOpt==2)KCSTAR(j)=0 !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				!if(iComplex==0 .and. iEosOpt==2 .and. LOUDER)pause 'iComplex=0 for assoc???' !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				if(LOUDER)write(*,602)'GetEsdCAS:iGotIt! id,bVol,eAcc=',ID(j),Vx(J),eAccEpsK(I)
				if(louder)print*,'GetEsdCas: nTypes,nDegree,nAcc,nDon=',nTypes(j),nDegree(j,1),nAcceptors(j,1),nDonors(j,1)
				exit !exits this do loop, not the outer one.
			ENDIF
		enddo
		if(iGotIt(J)==0)then
			if(ierCompExact( j).ne.0 )then
				iErr=11 !Parms missing and iErrExact.ne.0 for at least one component 
				if(LOUDER)write(*,*)'GetESD:Parms missing and iErrExact.ne.0 for component = ', j
				goto 861
			else
				if(LOUD)write(*,*)'Corr. States used for EsdParms of component=',j
			endif
        end if
    enddo
601 format(1x,a,8e12.4)
602 format(1x,a,i11,8e12.4)
	if(iErr)then
		if(LOUDER)pause 'GetEsd2: error for at least one compound'
		continue
	endif
        
	if(LOUDER)then
        write(*,*)'  ID     NAME       mESD    eok      bVol    Nd   KADnm3   eHB/k(kK)'
	    do i=1,nC
		    write(*,606)IDCas(i),NAME(i),q(i),eokP(i),Vx(i),NDS(i),NAS(i),KadNm3(i),epsD_kB(i),epsA_kB(i) !/1000
        enddo
    end if
606	format(i9,1x,a11,f9.3,f8.2,f8.2,2i3,1x,f8.6,2f8.0)

	!note:  bips are passed back through common/BIPs/
	IF(DEBUG)then 
		bipFile='c:\SPEAD\CalcEos\input\BipEsd96.txt'
		if(isMEM2)bipFile='c:\SPEAD\CalcEos\input\BipEsdMEM2.txt'
	ELSE 
		bipFile=TRIM(masterDir)//'\input\BipEsd96.txt' ! // is the concatenation operator
		if(isMEM2)bipFile=TRIM(masterDir)//'\input\BipEsdMEM2.txt' ! // is the concatenation operator
	ENDIF
	if(nC > 1)iErrCode=GetBIPs(bipFile,ID,nC) !not necessary for pure fluids
	if(iErrCode > 10)iErr=11 ! 
    if(LOUDER)then
		print*,'GetEsd2Cas: bipFile=',TRIM(bipFile)
		print*,'     ',(id(j),j=1,NC)
		do i=1,NC
			write(*,'(i5,11f8.4)')id(i),(Kij(i,j),j=1,NC)
		enddo
	    pause 'GetESDCas: check BIPs.'
		if(iErrCode > 0) print*,'GetEsdCas: BIPs missing for ',iErrCode-10,' binary combinations.'
    end if
	RETURN
	
861	continue
	!trap file reading errors
	if(LOUDER)write(*,*)'GetEsd Error - error reading EsdParms.txt'
	if(LOUDER)write(*,*)'nDeck,iCompo',NDECK,I

	!	do i=1,NC
	!		ID(i)=idStore(i)
	!	enddo

	return
end	!GetEsdCas
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------  
subroutine ExactEsd(nC,VX,c,q,eokP,iErr,ierComp)
	USE GlobConst
	implicit doubleprecision(A-H,K,O-Z)
	!  compute ESD2 parameters for hydrocarbons based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
	!parameter (NMX=55)
	DoublePrecision k1
	dimension VX(nmx),c(nmx),Q(nmx),eokP(nmx),ierComp(NC) !,ID(nmx)
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
	ierComp=0
	iErr=0				  
	do i=1,nC
		isAssoc=0
		if(TRIM(class(i))=='assoc' .or. TRIM(class(i))=="Asso+")isAssoc=1
		if(isAssoc)then
			ierComp(i)=1
			iErr=100+i
			if(LOUDER)print*,'ExactESD2: Parms not available for assoc ID=',ID(i)
			cycle
		endif	
		isHelium=0
		if( id(i)==913 .or. ID(i)==7440597)isHelium=1
		if(isHelium)then
			iErr=3
			if(LOUDER)print*,'ExactESD: Parms not available for helium ID=',ID(i)
			ierComp(i)=3
			cycle
		endif	
		isH2=0
		if( id(i)==902 .or. ID(i)==133740 .or. id(i)==925 .or. id(i)==7782390)isH2=1
		if(isH2)then
			iErr=2
			if(LOUDER)print*,'ExactESD: Parms not available for H2 or D2 ID=',ID(i)
			ierComp(i)=2
			cycle
		endif	

		Wci=ACEN(i)
		cShape=1+3.535*Wci+	0.533*Wci*Wci
		RooTCinv=1/SQRT(cShape)
		k1=1.7745
		qShape=1+1.90476*(cShape-1)
		ZcTmp=1.d0/3.d0+RooTCinv*(.0384+RooTCinv*(-.062+RooTCinv*(0.0723-.0577*RooTCinv)))
		atemp=9.5*qShape*1.9+4*cShape*k1-k1*1.9
		quadB=k1*1.9*ZcTmp+3*atemp
		sqArg=quadB*quadB+4*atemp*(4*cShape-1.9)*(9.5*qShape-k1)/ZcTmp

		Bc=ZcTmp*ZcTmp*(-quadB+SQRT(sqArg))/(2*atemp*(4*cShape-1.9))
		Yc=ZcTmp*ZcTmp*ZcTmp/(atemp*Bc*Bc)
		rlnY1=LOG(Yc+1.0617)
		c(i) = cShape
		Q(i) = qShape
		eokP(i)=TC(i)*rlnY1
		bVolCc_mol(i)=8.314*TC(i)/PC(i)*Bc
		VX(i)=bVolCc_mol(i)
	enddo
	return
end	!subroutine ExactEsd
!------------------------------------------------------------------------------ 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!------------------------------------------------------------------------------ 
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
	!	Example 2b.  ~MeOH+~EtOH at 393.15K,0.1MPa, Hij=0, Kij=0.002,~Ref[2]corrected.(The typos in <qYb> were beyond help and the hBonding has completely changed along with b,c,e.)
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
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	  + GlobConst{rGas,Tc,Pc,...}
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision xFrac(NMX),FUGC(NMX) !,chemPoAssoc(nmx)
	Integer ier(12) 
	LOGICAL LOUDER
	!  ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!  eokP IS THE DISPERSE ATTRACTION OVER BOLTZ k FOR PURE SPECIES
	!  KCSTAR IS THE BONDING VOLUME IN NM^3 
	!  DH IS THE BONDING ENERGY /RTC 
	!  C,Q,bVol ARE THE PURE COMPONENT EOS PARAMETERS
	!  KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!  zFactor IS PV/NoKT HERE  
	!  ier = 1 - AT LEAST ONE ERROR
	!        2 - NOT USED
	!        3 - NOT USED
	!        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!        5 - rho IS -VE
	!        6 - TOO MANY Z ITERATIONS
	!        7 - eta > 0.53
	!		11 - goldenZ instead of real Z.
	!
	!COMMON/ETA/ETAL,ETAV,ZL,ZV
	DATA INITIAL/1/  !Zm=9.5 since 1996, at least.  It is 9.5 in the CHEMCAD documentation, the EL text, and 2002.  It is not mentioned in ref[2,3,4,5,6]
	LOUDER=LOUD
	!LOUDER=.TRUE.
	if(LOUDER)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUDER)print*,'FugiEsd: Kij(1,2)= ',checkKij12
	ier(1:6)=0 !vector initialization
	!NOTE: iErrTmin is checked in FuVtot
	sumx=SUM( xFrac(1:NC) )
	if(ABS(sumx-1) > 1e-8)then
		if(LOUDER)pause 'FugiEsd: sumx .ne. 1'
	endif
	!INITIATE SECANT ITERATION ON rho
	bMix=SUM( xFrac(1:NC)*bVolCC_mol(1:NC) )
	Pb_RT=pMPa*bMix/(Rgas*tKelvin)
	!GUESS FOR rho
	eta=Pb_RT/1.05D0  !NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
	IF(LIQ==1 .or. LIQ==3 .or. eta>etaMax)eta=etaMax/1.15d0
	rho=eta/bMix
	if(eta > 1/1.9 .and. LOUDER)print*,'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
	isZiter=1 ! FUGC calculations are skipped for isZiter=1
	Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	IF(iErr > 10)GOTO 86
	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor

	eta=etaOld/1.15D0
	IF (eta < 0 .and. LOUD) WRITE(6,31)LIQ
	rho=eta/bMix
	if(initial.and.LOUD)print*,'FugiEsd: initial eta,err',etaOld,errOld
	itMax=77
	errBestEta=1234
	do nIter=1,itMax
		Call FuEsdVtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		IF(iErr > 10)EXIT
		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial.and.LOUDER)write(*,'(a,2e11.4,3f10.5)')'FugiEsd eta,Z', eta,zFactor
		if(initial.and.LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd eta,CHNG,niter',eta,CHNG,niter 
		etaOld=eta
		ERROLD=ERR
		!  LIMIT THE CHANGE IN Density for liquid..
		IF(liq==1.and.DABS(CHNG/etaOld) > 0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==3.and.DABS(CHNG/etaOld) > 0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		!  Low eta must move from zero, so you can't limit its % change
		IF(liq==0.and.DABS(CHNG) > 0.02d0)CHNG=DSIGN(0.02D0,CHNG)
		IF(liq==2.and.DABS(CHNG) > 0.02d0)CHNG=DSIGN(0.02D0,CHNG)
 		eta=eta-CHNG
		if(ABS(CHNG) < errBestEta)then
			etaBest=eta
			errBestEta=ABS(CHNG)
		endif
		if(eta < 0 .or. eta > 1/1.9)eta=etaOld-DSIGN(0.1D0,CHNG)*etaOld
		IF(DABS(CHNG) < 1.D-9 .and. eta > 0)EXIT  ! Don't use CHNG/eta here. Converge quickly to ideal gas if P->0, still ~9 sigfigs if liquid.
	enddo !nIter=1,itMax
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Iteration Concluded    !!!!!!!!!!!!!!!!!!!!!!!!!
	if(eta < 0 .or. eta > 1.9)then
		if(LOUDER)pause 'eta < 0 or > 1.9 on final iteration in FugiESD.'
		continue
	endif
	!One last call to get FUGC.
	if(initial.and.LOUD)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd2 cnvrgd: eta,CHNG,niter',eta,CHNG,niter
	etaPass=eta
	rho=eta/bMix
	rhoMol_cc=rho 
	IF (rho < 0)THEN
		ier(5)=1
        ier(1)=15
		if(LOUDER)WRITE(6,31)LIQ
		GOTO 86
	ENDIF
!  ITERATION ON rho HAS CONCLUDED.  DEPARTURES PASSED THROUGH GlobConst.
	if(eta > 0.43)write(*,'(a,i3,F8.2,2F8.5)')' FugiESD: converged. nIter,eta,T,x1=',nIter,tKelvin,xFrac(1),eta
	if(ABS(eta-rho*bMix) > 1E-11 .and. LOUDER)pause 'eta.ne.(rho*bMix)?'
	if(pMPa==0 .and. LOUDER)print*,' FugiEsd: P=0? LIQ,P=',LIQ,pMPa
	!zFactor=P/(rho*rGas*T)  ! add this to improve precision when computing rho from Z after return.   
	isZiter=0
	Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	if(zFactor.le.0)then
		ier(1)=11
		if(LOUDER)print*,'FugiEsd: converged Z <= 0. eta,Z=',eta,zFactor
		goto 86
	endif
	if(iErr > 0.or.nIter > itMax-1 .or. eta < 0)then ! if iErr still > 0 on last iteration, then declare error.
		ier(4)=iErr
		eta=etaBest
		ier(11)=1
		ier(1)=11
	endif
	FUGC(1:NC)=FUGC(1:NC)-DLOG(zFactor)	 !Must subtract ln(Z) when given Vtot as independent variable.
	if(LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FugiEsd: eta,CHNG,niter,FUGC',eta,CHNG,niter,(FUGC(i),i=1,NC) 
	initial=0
	return
86	if(LOUDER)WRITE(6,*)' ERROR IN FugiEsd.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(NITER.GE.ITMAX)THEN
		if(LOUDER)write(*,*)'TOO MANY Z ITERATIONS'
		ier(6)=1
        ier(1)=16
	END IF
	IF(ier(4)==1.and.LOUDER) WRITE(*,*)'ERROR IN FuEsdVtot'
	if(ier(1) < 10)ier(1)=11
	initial=0
	RETURN
	END	!Subroutine FugiESD()
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	Written originally by JRE, Oct. 2019																				C
!C	Given T,V and gmol(), this routine calculates rho, Zfactor, Ares_RT, Ures_RT, lnPhi 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,Ares,Ures,iErr)
	USE GlobConst
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(NMX),xFrac(NMX),FUGC(NMX),k10,k2,Zm !,KCSTARp(NMX),KVE(NMX)
	!integer iera(12) !IER(12),
	DoublePrecision YQVIJ(NMX,NMX),YQVI(NMX),Y(NMX,NMX),EOK(NMX,NMX)
	DoublePrecision CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX)
	Integer initCall
	DoublePrecision k1(NMX) !,RALPHA(NMX),RALPHD(NMX)
	!DoublePrecision FUGASSOC(NMX),FUGREP(NMX),FUGATT(NMX)			  
	!double precision moleStep
	!DoublePrecision dEOK_dT(NMX,NMX),dY_dT(NMX,NMX),dYQVI_dT(NMX),dYQVIJ_dT(NMX,NMX)
	!DoublePrecision dCVM_dN(NMX),dYQVM_dN(NMX) ,dVM_dN(NMX),dK1YVM_dN(NMX) !,dCVI_dN(NMX,NMX)
	!DIMENSION dZREP_dN(NMX),dZATT_dN(NMX),dZASSOC_dN(NMX),dZ_dN(NMX),dP_dN(NMX),dFREP_dN(NMX),dFATT_dN(NMX),dYQVI_dN(NMX,NMX)
	!DIMENSION dFUGC1_dP(NMX),dFUGC1_dT(NMX),dFUGC2_dP(NMX),dFUGC2_dT(NMX),dFUGASSOC_dRHO(NMX)
	!DIMENSION dFUGREP_dN(NMX,NMX),dFUGATT_dN(NMX,NMX),dFUGC1_dN(NMX,NMX),dFUGC2_dN(NMX,NMX),dFUG_dN(NMX,NMX)
	!DIMENSION dh_dN(NMX),dFUGASSOC_dT(NMX),dFUGASSOC_dN(NMX,NMX)
	!DIMENSION dh_dN_num(NMX),dFUGASSOC_dN_num(NMX,NMX),fugassocLoop(NMX),gmol_old(NMX)
	!DoublePrecision vMolecNm3(NMX) !for Wertheim
	LOGICAL LOUDER
	
!	COMMON/ETA/ETAL,ETAV,ZL,ZV
!	COMMON/ETA2/ETA
!	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
!	COMMON/ALPHA/ALPHAD(NMX,NMX), ALPHDA(NMX,NMX)
!	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(nmx)
!	COMMON/Derv1/dbVolMix_dN(NMX),DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
!!	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
!	COMMON/dFug/h_nMichelsen,dh_dT,dh_dN,dFugassoc_dT,dFugassoc_dN,dFUGASSOC_dRHO,dh_dRHO!,FUGASSOC_num
!	COMMON/dfug2/dFUGASSOC_dN_num,dh_dN_num
!	COMMON/num/FUGASSOC_num
!	COMMON/fugCR/PMpa,dFUG_dN,dP_dN
!	COMMON/Helm_derv/Ur,UrV,FREP,FATT,FassocAfg
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),ralph(nmx),Zrep,Zatt,Zassoc,Fassoc
	DATA K10,K2,ZM,initCall/1.7745D0,1.0617D0,9.5D0,1/

	LOUDER=LOUD
	!LOUDER=.TRUE.
	  
	iErr=0 !1=warning from AlpSolEz2, 11=input nonsense, 12=xFrac<0, 13=critical error from AlpSol, 14=voidFrac<0..
	totMoles=sum( gmol(1:NC) )
	xFrac(1:NC)=gmol(1:NC)/totMoles

	if( tKelvin	   < zeroTol .or. totMoles < zeroTol .or. vTotCc < zeroTol)then
		if(LOUD)print*,'FuEsdVtot: nonsense T(K),totMoles,vTotCc=',tKelvin,totMoles,vTotCc
		iErr=11
	endif
	iErrTmin=0
	TminTot=.01
	TrVolatile=0.1
	DO I=1,NC
		rLogPr=DLOG10( 0.0001/Pc(i) )	! ESD not recommended below 0.0001 MPa for pure fluids. For mixes, ensure most volatile comp has Psat > 0.0001 MPa. Think about polymers.
		aScvp=7*(1+acen(i))/3	 !SCVP: log10(Pr)=7(1+w)/3*(1-1/Tr) = a*(1-x) 
		xt1= 1-rLogPr/aScvp  ! x = 1/Tr at Psat=0.0001 MPa, first approximation of x = x1 + dx
		xt = xt1 -0.178*acen(i)*acen(i)*(1-xt1)*(1-xt1)  ! This crude empirical correlation was developed for nonadecanol.  cf. PGL6Samples.xlsx(nC19oh).
		if( xt > 2.222)xt = 2.2222	!1/2.2222 = 0.45. If 
		if( xt < 1 .and. LOUDER)pause 'FugiEsd: TrMin > Tc???'
		TrMin = 1/xt ! = min( Tr@Psat=0.0001 or 0.45 )
		if(LOUDER)write(*,601)' xt1,xt,TrMin',xt1,xt,TrMin
		if( tKelvin/ Tc(i) < TrMin .and. NC==1)iErrTmin=iErrTmin+1
		if( tKelvin/Tc(i) > TrVolatile) TrVolatile=tKelvin/Tc(i)  ! The largest TrVolatile is the Tr of the compd with lowest Tc. 
		if( Tc(i)*TrMin > TminTot) TminTot=Tc(i)*TrMin	 ! The largest Tmin is the weakest link. 
		if( Tc(i)*TrMin > TminTot .and. LOUDER) print*,'i,Tmin(i): ', i,Tc(i)*TrMin
		IF(xFrac(I) < 0 .and. LOUDER)PAUSE 'FugiEsd: ERROR - Xi<0'
	enddo 
	if(TrVolatile < 0.4d0)iErrTmin =2 ! it's only a problem if the most volatile compound has Tr < 0.45 or Psat < 0.0001.
	if(iErrTmin > 0) then
		iErr=5 ! warning level because functions like Vxs or Hxs might be insensitive to this issue.
		if(LOUDER)print*,'FuEsdVtot: T(K) < Tmin(all i)',tKelvin,TminTot
		!if(LOUDER) pause 'FugiEsd: at least one compound has Tr < TrMin'
	endif
	bMix=0
	do i=1,nc
		xFrac(i)=gmol(i)/totMoles
		IF(xFrac(i) < 0)then
			if(LOUDER)PAUSE 'ERROR IN FuEsdVtot, Xi<0'
			iErr=12
		endif
		bMix=bMix+xFrac(i)*bVolCC_mol(i)
	enddo
	if(iErr > 10)return
	rho=totMoles/ vTotCc
	eta=rho*bMix 
	if(LOUDER.and.initCall)write(*,601)' FuEsdVtot: T,x1,bMix,eta=',tKelvin,xFrac(1),bMix,eta
601 format(1x,a,8E12.4)
	YQVM=0.d0
	VM=0.d0
	CVM=0.d0
	Cmix=0.d0
	K1YVM=0
	iType=1	  
	DO I=1,NC
		K1(I)=K10
		Cmix=Cmix+xFrac(i)*C(i)
		DO J=1,NC
			kIJbip=KIJ(I,J)+KTIJ(I,J)/tKelvin 
			EOK(I,J)=DSQRT(EOKP(I)*EOKP(J))*(1.d0-kIJbip)
			Y(I,J)=DEXP(EOK(I,J)/tKelvin)-K2
			QV(I,J) = (Q(I)*VX(J) + Q(J)*VX(I)) / 2.d0
			YQVIJ(I,J)=QV(I,J)*Y(I,J)
			CVIJ(I,J) = (C(I)*VX(J) + C(J)*VX(I)) / 2.d0 ! e.g. (x1*c1+x2*c2)*(x1*b1+x2*b2) = x1^2*c1*b1+x1*x2*(c1*b2+c2*b1)+x2^2*b2^2
			YQVM=YQVM+YQVIJ(I,J)*xFrac(I)*xFrac(J)	   
			CVM = CVM + CVIJ(I,J)*xFrac(I)*xFrac(J)	   !note: above means <c>=sum(xi*ci) and <b>=sum(xj*bj) 
		enddo
		!dHkcalMol(I)=DH(I)*TC(I)*1.987D-3	!JRE 12/00
		!KVE(I)=KCSTAR(I)*( DEXP(dHkcalMol(I)/tKelvin/1.987D-3)-1.d0 )  *avoNum  !KVE[=] cc/mol, KcStar [=] nm3.
		!	if(initCall.and.LOUD)print*,'FuESDVtot: KcStar=',KcStar(i)
		!if(initCall.and.LOUD)pause 'Right?'
		!KCSTARp(I)=KCSTAR(I)
		VM=VM+xFrac(I)*VX(I)
		K1YVM=K1YVM+xFrac(I)*K1(I)*Y(I,I)*VX(I) !1991 form, overwritten if applying 1990 form
		!vMolecNm3(i)=VX(I)/avoNum
		!bondVolNm3Esd(I)=KCSTAR(I) !*vMolecNm3(I)
		!eHbKcal_mol(I,iType)=dHkcalMol(I)
		!bondVolNm3(I,iType)=bondVolNm3Esd(I)
	enddo
	if( ABS(VX(1) - bVolCC_mol(1)) > zeroTol .and. LOUD ) write(*,601)' FuEsdVtot: VX.ne.bVol=',VX(1),bVolCC_mol(1)
	if(LOUD.and.k1yvm < zeroTol)print*,'FuEsdVtot: 0~k1yvm=',k1yvm 
	eta=rho*vm
	if(isMEM2)then
		CALL MEM2(isZiter,tKelvin,xFrac,NC,rho,zAssoc,aAssoc,uAssoc,FUGASSOC,iErrMEM )!,ier)
	else
		CALL MEM1(isZiter,tKelvin,xFrac,NC,rho,zAssoc,aAssoc,uAssoc,fAssoc,FUGASSOC,iErrMEM )!,aAssoc,uAssoc,rLnPhiAssoc,ier)
	endif
!	CALL AlpSolEz2(isZiter,tKelvin,xFrac,Nc,VM,rho,zAssoc,aAssoc,uAssoc,FUGASSOC,iera)
!	IF(iera(1) > 0)then
!       if(iera(1)>10)iErr=13 ! This is a warning level for iErr when AlpSol did not converge.
!		if(iera(1)==1)iErr=1
!   endif
	if(LOUDER)write(*,601)' FuEsdVtot: rho,zAssoc=',rho,zAssoc
	if(iErrMEM > 0 .and. LOUDER)write(*,601)' FuEsdVtot: iErrMEM > 0. rho,zAssoc=',rho,zAssoc
	if(iErrMEM==1)iErr=1
	if(iErrMEM > 10)iErr=13
	voidFrac=1-1.9D0*ETA
	denom=voidFrac
	ZREP= 4.d0*Cmix*eta/denom
	ZATT= -ZM*YQVM*RHO/(1+K1YVM*RHO)
	Z=(1+ZREP+ZATT+ZASSOC)
	PMpa=Z*rGas*RHO*tKelvin
    if(voidFrac < 0)then
	    IF(LOUD)print*, 'FuEsdVtot:Error! (1-1.9*ETA) IS -VE. eta,rho=',eta,rho
		iErr=14
    endif
	if(iErr > 10)return



	if(isZiter==1)return ! don't need the rest if isZiter.



 	DO I=1,NC
	   YQVI(I)=0.D0
	   CVI(I)=0.D0
	enddo     
	UATT=0.d0
	DO I=1,NC
		!ralph(I)=SQRT(alphAD(I,I))	 ! ralph is computed in alpsolEz.
		UATT=UATT+xFrac(I)*VX(I)*EOK(I,I)/tKelvin*(Y(I,I)+K2)*K1(I)       
		DO J=1,NC
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*xFrac(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*xFrac(J)
 		enddo
	enddo
	FATT= -ZM*YQVM/K1YVM*DLOG(1.d0+K1YVM*RHO)
	FREP= -4.d0/1.9D0*DLOG(voidFrac)*CVM/VM
	FREP= -4.d0/1.9D0*DLOG(voidFrac)*Cmix
	UATT= -ZM*YQVM*RHO/(1+K1YVM*RHO)*UATT/K1YVM	!uAtt=...uAtt/k1YVM where uAtt defined in 3rd line of above do loop.

	!     CALLING ASSYMU IF NUMBER OF DONOR SITES AND ACCEPTOR SITES ARE NOT EQUAL 
	!      IF(IFLAG.EQ.2)CALL ASSYMU(ALPHAD,ALPHDA,XA,XD,X,ND,NC,QIJ,UASSOC)

	DUONKT=UATT+uASSOC
	DAONKT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	DSONK =DUONKT-DAONKT
	DHONKT=DUONKT+Z-1.d0
	uRes_RT=UATT+UASSOC
	aRes_RT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	!print*,'aRes_RT=',aRes_RT
	sRes_R =UATT+UASSOC-aRes_RT
	hRes_RT=UATT+UASSOC+Z-1
	Ares=DAONKT
	Ures=DUONKT
	if(LOUDER)write(*,601)' FuEsdVtot: zAssoc,aAssoc,uAssoc=',zAssoc,aAssoc,uAssoc

	if (isZiter==0) then
		!call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
        if( Z < zeroTol)then  ! Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
		    if(LOUDER)pause 'FuEsdVtot: Z.le.0 for fugc calculation.'
			iErr=3	  ! warning level because another call might produce Z > 0.
			goto 86
        endif
		if(LOUDER)write(*,'(a,f10.5)')' i,lnGamRep,lnGamAtt,lnGamBon.'
		DO I=1,NC
			FUGREP(I)=FREP*( 2.d0*C(I)/Cmix-VX(I)/VM ) + ZREP*VX(I)/VM
			FUGREP(I)=FREP*( C(I)/Cmix ) + ZREP*VX(I)/VM ! For pure i, FugRepi= -4ci/1.9*ln(1-1.9eta) + 4ci*eta/(1-1.9eta) 
			FUGATT(I)=FATT*( 2*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM )+ZATT*K1(I)*Y(I,I)*VX(I)/K1YVM !91-pres form, complete w/o dNk1YbNk, cf. EL99 p559 & S&E91apx
			!FUGASSOC(i)=ND(i)*2*DLOG(XA(i,1)) + zAssoc*1.9D0*VX(i)*rho !JRE'96 Eq.43.
			FUGC(I)=FUGREP(I)+FUGATT(I)+FUGASSOC(I)  ! -DLOG(Z)  Don't subtract ln(Z) when given Vtot as independent variable.
			rLnGamRep=FUGREP(i)-Zrep*Vx(i)/VM  ! cf. Bala and Lira (2016), Eqs A6-A14. to correct from constant volume to P.
			rLnGamAtt=FUGATT(i)-Zatt*Vx(i)/VM
			rLnGamBon=FUGASSOC(i)-Zassoc*Vx(i)/VM
			IF(LOUDER)WRITE(*,'(i3,f7.4,9f10.4)')i,xFrac(i),rLnGamRep,rLnGamAtt,rLnGamBon !,ralpha(i),ralphd(i)
		ENDDO
	endif
	if(LOUDER)pause 'FuEsdVtot: Check results before returning.'
86	return
	end !subroutine FuEsdVtot
