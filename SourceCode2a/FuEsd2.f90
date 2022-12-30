!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MODULE EsdParms: Included in ESD96
!	USE GlobConst
!	DoublePrecision EOKP(NMX),KCSTAR(NMX),DH(NMX),C(NMX),Q(NMX),VX(NMX)
!	Integer         ND(NMX),NDS(NMX),NAS(NMX)
!END MODULE EsdParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine GetEsd2Cas(nC,idCasPas,iErr) !ID is passed through GlobConst
	!
	!  PURPOSE:  LOOKS UP THE ESD PARAMETERS AND STORES THEM IN COMMON EsdParms
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT(to EsdParms)   qShape,eokP,VX,KadNm3,epsA_kB,epsD_kB,ND,NDS,NAS, 
	USE GlobConst !is implied by USE ASSOC. GlobConst includes ID,Tc,Pc,Zc,acen,...
	!USE Assoc (not needed for self-contained ESD2 w/ grand master eq.) !nTypes(=1 for ESD), nDegree, nDonors, nAcceptors, bondRate(?), includes GlobConst{Tc,Pc,...}
	USE EsdParms 
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=3000,listPool=1000)
	character*222 bipFile,inFile,dumString !,dumString,ParmsTptFile*50 !,bipHbFile*50
	integer iGotIt(NC),idCasPas(NC),idCasa(ndb),GetBIPs,ierCompExact(NC)  !,idStore(NC) 
	doublePrecision kcSta
	DIMENSION IDA(ndb),QA(ndb),KCSTA(ndb),eDonEpsK(ndb),eAccEpsK(ndb),NDSA(ndb),NASA(ndb),bVolA(ndb),eokA(ndb) !,NDA(ndb)
	!doublePrecision bondRate(nmx,maxTypes)
	!iErr=1 !Parms missing and iErrExact.ne.0 for at least one component
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
    do i=1,NC
        idCas(i)=idCasPas(i)
    enddo
	!idCas(1:NC)=idCasPas(1:NC) ! workaround after promoting idCas to GlobConst 
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isTPT or isPcSaft.
	isESD=.TRUE. ! in GlobConst, simplifies calls in FuVtot or FUGI 
	etaMax=1/1.9D0-zeroTol
!	nComps=1
!	do i=1,nComps
!		idStore(i)=id(i)
!		id(i)=idCas(i)
!	enddo
	IF(DEBUG)then 
		inFILE='c:\Spead\CalcEos\input\ParmsEsd2.TXT'
		!if(iEosOpt==12)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmami.TXT'
		!if(iEosOpt==13)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmamiTb.TXT'
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsEsd2.txt' ! // is the concatenation operator
		!if(iEosOpt==12)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmami.txt' ! // is the concatenation operator
		!if(iEosOpt==13)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmamiTb.txt' ! // is the concatenation operator
	ENDIF
	OPEN(31,FILE=inFile)
	if(LOUDER)print*,'GetEsd2:inFile=',TRIM(inFile)
	if(LOUDER)pause 'Check the ESD2 parms file location.'

	!READ(31,'(a88)',ERR=861)dumString 	!read title line. Ignore nDeck if present.
	READ(31,*)nDeck
	do i=1,nDeck
		READ(31,'(a222)')dumString
		READ(dumString,*,ioStat=ioErr)IDA(I),QA(I) ,eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I),NDSA(I),NASA(I),idCasa(i)
		!write(*,*)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I)  ,idCasa(i)
		if(ioErr .and. LOUDER)print*,'GetESD2: error reading ParmsEsd_.txt. line=',TRIM(dumString)
		if(  ( idCasa(i)==id(1) .or. idCasa(i)==id(2) ) .and. LOUDER  )print*,'Found in ParmsEsd idCas=',idCasa(i) 
	enddo
	NDI=1
	I= -1  !overrides headerless reading. Set I=0 to re-activate.
	DO while(I.ge.0) !reads to end of file w/o reading nDeck
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
	if(iEosOpt == 4)then
		ierCompExact=1 ! Set error for all comps. Declare error because iEosOpt > 4 means using only GC parameters from one of ParmsEsdEmami__
	else
		if(Tc(1) < 4)call GetCritCas(nC,iErrCrit) !GetCritCas assumes ID(GlobConst)=IdCas
		call ExactEsd2(nC,VX,C,Q,eokP,iErrExact,ierCompExact) !iErrExact = 100+iComp if compd is assoc or Asso+. Wait to see if Parms are in ParmsEsd before failing.
		!mShape(1:NMX)=Q(1:NMX)
		if(LOUDER.and.iErrExact)print*,'GetESDWarning: iErrExact=',iErrExact,' Checking database for iComp='	,( ierCompExact(iComp),iComp=1,NC)
	endif

	DO J=1,nC
        bVolCC_mol(j)=Vx(j) !Copy ExactEsd value first. Replaced below if in dbase.
		ND(J)=0
		NDS(J)=0
		NAS(J)=0
		iGotIt(J)=0
		DO I=1,NDECK
			IF(IDCASA(I).EQ.IDCas(J))THEN
				iGotIt(J)=1
				eokP(J)=eokA(I)
				!ND(J)=NDA(I)
				Q(J)=QA(I)
				C(J)=1+(Q(J)-1)/1.90476
				Vx(J)=bVolA(I)
				bVolCC_mol(J)=Vx(J) !need bVol generally for every EOS. Included in GlobConst
				KCSTAR(J)=KCSTA(I) ! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				!DH(J)=DHA(I)*1000/1.987/ Tc( j)			! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				epsA_kB(J)=eAccEpsK(I)
				epsD_kB(J)=eDonEpsK(I)
				NAS(J)=NASA(I)
				NDS(J)=NDSA(I)
				KadNm3(j)=KcStar(j) !KcStar[=] nm^3 since 2021.
				iComplex=0
				!if(NAS(j)*NDS(j) > 0)iComplex=1
				!if(iComplex==0 .and. iEosOpt==2)KCSTAR(j)=0 !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				!if(iComplex==0 .and. iEosOpt==2 .and. LOUDER)pause 'iComplex=0 for assoc???' !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				if(LOUDER)print*,'iGotIt! idCas,bVol=',idCas(j),Vx(J)
				exit !exits this do loop, not the outer one.
			ENDIF
		enddo
		if(iGotIt(J)==0)then
			if(ierCompExact( j).ne.0 )then
				iErr=1 !Parms missing and iErrExact.ne.0 for at least one component 
				if(LOUDER)write(*,*)'GetESD:Parms missing and iErrExact.ne.0 for component = ', j
				goto 861
			else
				if(LOUD)write(*,*)'Corr. States used for EsdParms of component=',j
			endif
        end if
    enddo
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
		bipFile='c:\SPEAD\CalcEos\input\BipEsd2.txt'
	ELSE 
		bipFile=TRIM(masterDir)//'\input\BipEsd2.txt' ! // is the concatenation operator
	ENDIF
	if(nC > 1)iErrCode=GetBIPs(bipFile,ID,nC) !not necessary for pure fluids
	if(iErrCode > 10)iErr=11 ! 
    if(LOUDER)then
		print*,'GetEsd2Cas: bipFile=',TRIM(bipFile)
		print*,'     ',(id(j),j=1,NC)
		do i=1,NC
			write(*,'(i5,11f8.4)')id(i),(Kij(i,j),j=1,NC)
		enddo
	    pause 'GetESD2Cas: check BIPs.'
		if(iErrCode > 0) print*,'GetESD2Cas: BIPs missing for ',iErrCode-10,' binary combinations.'
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
subroutine ExactEsd2(nC,VX,c,q,eokP,iErr,ierComp)
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
subroutine QueryParPureEsd2(iComp,iParm,value,iErr)
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
end	!Subroutine QueryParPureEsd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPureEsd2(iComp,iParm,value,iErr)
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
end	! SetParPureEsd
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

	SUBROUTINE FugiESD2(tKelvin,pMPa,xFrac,NC,LIQ,FUGC,zFactor,ier)
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
	ier=0 !vector initialization

	iErrTmin=0
	TminTot=.01
	TrVolatile=0.1
	if(LOUDER)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUDER)print*,'FugiEsd2: Kij(1,2)= ',checkKij12
	DO I=1,NC
		rLogPr=DLOG10( 0.0001/Pc(i) )	! ESD not recommended below 0.0001 MPa for pure fluids. For mixes, ensure most volatile comp has Psat > 0.0001 MPa. Think about polymers.
		aScvp=7*(1+acen(i))/3	 !SCVP: log10(Pr)=7(1+w)/3*(1-1/Tr) = a*(1-x) 
		xt1= 1-rLogPr/aScvp  ! x = 1/Tr at Psat=0.0001 MPa, first approximation of x = x1 + dx
		xt = xt1 -0.178*acen(i)*acen(i)*(1-xt1)*(1-xt1)  ! This crude empirical correlation was developed for nonadecanol.  cf. PGL6Samples.xlsx(nC19oh).
		if( xt > 2.222)xt = 2.2222	!1/2.2222 = 0.45. If 
		if( xt < 1 .and. LOUD)pause 'FugiEsd: TrMin > Tc???'
		TrMin = 1/xt ! = min( Tr@Psat=0.0001 or 0.45 )
		if(initial .and. LOUD)print*,'xt1,xt,TrMin',xt1,xt,TrMin
		if( tKelvin/ Tc(i) < TrMin .and. NC==1)iErrTmin=iErrTmin+1
		if( tKelvin/Tc(i) > TrVolatile) TrVolatile=tKelvin/Tc(i)  ! The largest TrVolatile is the Tr of the compd with lowest Tc. 
		if( Tc(i)*TrMin > TminTot) TminTot=Tc(i)*TrMin	 ! The largest Tmin is the weakest link. 
		if( Tc(i)*TrMin > TminTot .and. LOUD) print*,'i,Tmin(i): ', i,Tc(i)*TrMin
		IF(xFrac(I) < 0 .and. LOUD)PAUSE 'FugiEsd: ERROR - Xi<0'
	enddo 
	if(TrVolatile < 0.44d0)iErrTmin =2 ! it's only a problem if the most volatile compound has Tr < 0.45 or Psat < 0.0001.
	if(iErrTmin > 0) then
		ier(1)=5
		ier(2)=iErrTmin
		if(LOUD)print*,'FugiEsd: T(K) < Tmin(all i)',tKelvin,TminTot
		!if(LOUD) pause 'FugiEsd: at least one compound has Tr < TrMin'
		!return  !! make this a warning for Vex,Hex etc.  
	endif
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
	Call FuEsd2Vtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	IF(iErr > 10)GOTO 86
	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor

	eta=etaOld/1.15D0
	IF (eta < 0 .and. LOUD) WRITE(6,31)LIQ
	rho=eta/bMix
	if(initial.and.LOUD)print*,'FugiEsd2: initial eta,err',etaOld,errOld
	itMax=77
	errBestEta=1234
	do nIter=1,itMax
		Call FuEsd2Vtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		IF(iErr > 10)EXIT
		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial.and.LOUDER)write(*,'(a,2e11.4,3f10.5)')'FugiEsd2 eta,Z', eta,zFactor
		if(initial.and.LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd2 eta,CHNG,niter',eta,CHNG,niter 
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
	IF (rho < 0)THEN
		ier(5)=1
        ier(1)=15
		if(LOUDER)WRITE(6,31)LIQ
		GOTO 86
	ENDIF
!  ITERATION ON rho HAS CONCLUDED.  DEPARTURES PASSED THROUGH GlobConst.
	if(eta > 0.43)write(*,'(a,i3,F8.2,2F8.5)')' FugiESD2: converged. nIter,eta,T,x1=',nIter,tKelvin,xFrac(1),eta
	if(ABS(eta-rho*bMix) > 1E-11 .and. LOUDER)pause 'eta.ne.(rho*bMix)?'
	if(pMPa==0 .and. LOUDER)print*,' FugiEsd2: P=0? LIQ,P=',LIQ,pMPa
	!zFactor=P/(rho*rGas*T)  ! add this to improve precision when computing rho from Z after return.   
	if(zFactor.le.0)then
		ier(1)=11
		if(LOUDER)print*,'FugiEsd2: converged Z <= 0. eta,Z=',eta,zFactor
	endif
	isZiter=0
	rho=eta/bMix
	Call FuEsd2Vtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	if(iErr > 0.or.nIter > itMax-1 .or. eta < 0)then ! if iErr still > 0 on last iteration, then declare error.
		ier(4)=iErr
		eta=etaBest
		ier(11)=1
		ier(1)=11
	endif
	FUGC(1:NC)=FUGC(1:NC)-DLOG(zFactor)	 !Must subtract ln(Z) when given Vtot as independent variable.
	if(LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd2: eta,CHNG,niter,FUGC',eta,CHNG,niter,(FUGC(i),i=1,NC) 
	initial=0
	return
86	if(LOUDER)WRITE(6,*)' ERROR IN FuEsd96.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(NITER.GE.ITMAX)THEN
		if(LOUDER)write(*,*)'TOO MANY Z ITERATIONS'
		ier(6)=1
        ier(1)=16
	END IF
	IF(ier(4)==1.and.LOUDER) WRITE(*,*)'ERROR IN FuEsd2Vtot'
	if(ier(1) < 10)ier(1)=11
	initial=0
	RETURN
	END	!FugiESD()
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FugiESD2a(T,P,X,nC,LIQ,FUGC,zFactor,ier)
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	  + GlobConst{rGas,Tc,Pc,...}
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	integer kComp,QueryNParPure
	DoublePrecision X(NMX),FUGC(NMX),chemPoAssoc(nmx)
	!DoublePrecision mShapeMix
	integer ier(12),iera(12)
	DoublePrecision YQVIJ(NMX,NMX),KVE(NMX),YQVI(NMX),Y(NMX,NMX),eok(NMX,NMX)
	DoublePrecision CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX),dXsYbN(NMX)	!,XA(NMX),XD(NMX)
	DoublePrecision k1(NMX) !,RALPHA(NMX),RALPHD(NMX) !,rLnPhiAssoc(NMX)
	LOGICAL LOUDER
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),ralph(nmx),Zrep,Zatt,Zassoc,Fassoc
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/eta2/eta
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!      COMMON/xsGibbs/xsBIP(NMX,NMX),xsNrtlAl(NMX,NMX)
	!  ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!  eokP IS THE DISPERSE ATTRACTION OVER BOLTZ k FOR PURE SPECIES
	!  KCSTAR IS THE DIMENSIONLESS BONDING VOLUME FOR PURE 
	!  DH IS THE BONDING ENERGY /RTC 
	!  C,Q,bVol ARE THE PURE COMPONENT EOS PARAMETERS
	!  KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!  Z IS PV/NoKT HERE  
	!  ier = 1 - AT LEAST ONE ERROR
	!        2 - NOT USED
	!        3 - NOT USED
	!        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!        5 - rho IS -VE
	!        6 - TOO MANY Z ITERATIONS
	!        7 - eta > 0.53
	!		11 - goldenZ instead of real Z.
	!
	DATA k10,K2,zM,INITIAL/1.7745,1.0617,9.5,1/  !Zm=9.5 since 1996, at least.  It is 9.5 in the CHEMCAD documentation, the EL text, and 2002.  It is not mentioned in ref[2,3,4,5,6]
	LOUDER=LOUD
	!LOUDER=.TRUE.

	do i=1,11
		ier(i)=0
	enddo
	if(initial==1.and.LOUD)then
		initial=0
		nParPure=QueryNParPure()
		print*,'nParPure=',nParPure
		pause 'FugiEsd: check nParPure'
	endif
	!	if(initial.eq.0)then
	!	  write(*,*)'enter xsScale '
	!	  read(*,*)xsScale
	!	  write(*,*)'enter xsBIP(1,2),xsBIP(2,1) '
	!	  read(*,*)xsBIP(1,2),xsBIP(2,1)
	!	endif
	qYbMix=0.0
	bMix=0.0
	cbMix=0
	k1YbMix=0
	sumx=0
	cShapeMix=0
	iErrTmin=0
	sixth=1.d0/6.d0
	TminTot=.01
	TrVolatile=0.1
	if(LOUD)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUD)print*,'FugiEsd: Kij(1,2)= ',checkKij12
	DO I=1,NC
		sumx=sumx+x(i)
		rLogPr=DLOG10( 0.0001/Pc(i) )	! ESD not recommended below 0.0001 MPa for pure fluids. For mixes, ensure most volatile comp has Psat > 0.0001 MPa. Think about polymers.
		aScvp=7*(1+acen(i))/3	 !SCVP: log10(Pr)=7(1+w)/3*(1-1/Tr) = a*(1-x) 
		xt1= 1-rLogPr/aScvp  ! x = 1/Tr at Psat=0.0001 MPa, first approximation of x = x1 + dx
		xt = xt1 -0.178*acen(i)**2*(1-xt1)*(1-xt1)  ! This crude empirical correlation was developed for nonadecanol.  cf. PGL6Samples.xlsx(nC19oh).
		if( xt > 2.222)xt = 2.2222	!1/2.2222 = 0.45. If 
		if( xt < 1 .and. LOUD)pause 'FugiEsd: TrMin > Tc???'
		TrMin = 1/xt ! = min( Tr@Psat=0.0001 or 0.45 )
		if(initial .and. LOUD)print*,'xt1,xt,TrMin',xt1,xt,TrMin
		if( T/ TC(i) < TrMin .and. NC==1)iErrTmin=iErrTmin+1
		if( T/Tc(i) > TrVolatile) TrVolatile=T/Tc(i)  ! The largest TrVolatile is the Tr of the compd with lowest Tc. 
		if( Tc(i)*TrMin > TminTot) TminTot=Tc(i)*TrMin	 ! The largest Tmin is the weakest link. 
		if( Tc(i)*TrMin > TminTot .and. LOUD) print*,'i,Tmin(i): ', i,Tc(i)*TrMin
		IF(X(I) < 0 .and. LOUD)PAUSE 'FugiEsd: ERROR - Xi<0'
		k1(I)=k10
		DO J=1,nC
			kIJtemp=(KIJ(I,J)+KTIJ(I,J)/T) 
			eok(I,J)=DSQRT(eokP(I)*eokP(J))*( 1-kIJtemp )
			Y(I,J)=DEXP(eok(I,J)/T)-K2
			QV(I,J) = (Q(I)*Vx(J) + Q(J)*Vx(I)) / 2.0
			YQVIJ(I,J)=QV(I,J)*Y(I,J)
			CVIJ(I,J) = (C(I)*Vx(J) + C(J)*Vx(I))/2.0
			qYbMix=qYbMix+YQVIJ(I,J)*X(I)*X(J)
			cbMix = cbMix + CVIJ(I,J)*X(I)*X(J)								   
		enddo
		KVE(I)=KadNm3(I)*avoNum*( DEXP(epsA_kB(I)/T)-1 )  !this KVE is in cc/mol
		!if(initial.and.LOUD)print*,'FugiESD: KcStar=',KcStar(i)
		!if(initial.and.LOUD)pause 'Right?'
		bMix=bMix+X(I)*Vx(I)
		cShapeMix=cShapeMix+X(I)*C(I)
		k1YbMix=k1YbMix+X(I)*k1(I)*Y(I,I)*Vx(I)  !91 form, overwritten if applying 90 form
		dXsYbN(I) = 0
	enddo 
	if(TrVolatile < 0.44d0)iErrTmin =2 ! it's only a problem if the most volatile compound has Tr < 0.45 or Psat < 0.0001.
	if(iErrTmin > 0) then
		ier(1)=5
		ier(2)=iErrTmin
		if(LOUD)print*,'FugiEsd: T(K) < Tmin(all i)',T,TminTot
		!if(LOUD) pause 'FugiEsd: at least one compound has Tr < TrMin'
		!return  !! make this a warning for Vex,Hex etc.  
	endif
	if(ABS(sumx-1) > 1e-8)then
		if(LOUD)pause 'FugiEsd: sumx .ne. 1'
	endif
	! Mixing Rule T derivatives
	TdYQVM_dT=0.d0
	TdK1YVM_dT=0.d0
	DO I=1,NC 
		DO J=1,NC
	  		dEOK_dTij=DSQRT(EOKP(I)*EOKP(J))*(KTIJ(I,J)/(T*T)) 
			dY_dTij=(dEOK_dTij*T-EOK(I,J))*DEXP(EOK(I,J)/T)/(T*T)
			TdYQVM_dT=TdYQVM_dT+X(I)*X(J)*QV(I,J)*dY_dTij
			if(i==j) TdK1YVM_dT=TdK1YVM_dT+X(I)*K1(I)*VX(I)*dY_dTij
		enddo
    enddo

	!call xsMixRule(dXsYbN,xsYb,dXsYbB,bVol,X,T,nC)
	!k1YbMix=k1YbMix+k10*xsYb			!91 form ignores this.

	!qMix=1+1.90476*(cShapeMix-1) 		!90 form, not used otherwise
	!k1YbMix=k10*qYbMix/qMix			!90 form, comment out to eliminate, cf. EL99 p559 & S&E91apx
	!pause 'In FugiESD. Using 1990 form'

	    if(k1YbMix==0 .and. LOUD)pause 'FugiESD: k1YbMix=0' 
      
	!INITIATE SECANT ITERATION ON rho
	Pb_RT=P*bMix/(RGAS*T)
	!GUESS FOR rho
	eta=Pb_RT/1.05D0  !NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
	IF(LIQ==1 .or. LIQ==3 .or. eta>0.5)eta=.5D0
	rho=eta/bMix
	if(eta > 1/1.9 .and. LOUD)print*,'FugiEsd:etaInit > 0.53. P,T=',P,T 

	!ALPSOL WILL USE THIS VALUE OF rho TO CALCULATE A NEW VALUE OF zAssoc
	isZiter=1
	CALL AlpSolEz2(isZiter,T,X,Nc,bMix,rho,zAssoc,aAssoc,uAssoc,chemPoAssoc,iera)
	IF(iera(4)==1)GOTO 86
	voidFrac=(1-1.9D0*eta)
	zRep=4*cbMix*rho/voidFrac
	zAtt=-zM*qYbMix*rho/(1+k1YbMix*rho)
	zFactor=(1+zRep+zAtt+zAssoc)

	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor
	eta=etaOld*1.05D0
	IF (eta < 0 .and. LOUD) WRITE(6,31)LIQ
	rho=eta/bMix
	if(initial.and.LOUD)print*,'FugiEsd: initial eta,err',etaOld,errOld
	itMax=77
	nIter=0
	do nIter=1,itMax
		ier(2)=0
		ier(3)=0
		CALL AlpSolEz2(isZiter,T,X,Nc,bMix,rho,zAssoc,aAssoc,uAssoc,chemPoAssoc,iera)
		IF(iera(4)==1)EXIT
		voidFrac=(1-1.9D0*bMix*rho)
		zRep=4*cbMix*rho/voidFrac
		zAtt=-zM*qYbMix*rho/(1+k1YbMix*rho)
		zFactor=(1+zRep+zAtt+zAssoc)

		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial.and.LOUD)write(*,'(a,2e11.4,3f10.5)')' rho,Z,zAssoc,fAssoc', rho,zFactor,zAssoc,fAssoc
		if(initial.and.LOUD)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' eta,CHNG,niter,ralph',eta,CHNG,niter,(ralph(i),i=1,NC)
		etaOld=eta
		ERROLD=ERR
		!  LIMIT THE CHANGE IN Density for liquid..
		!  Low eta must move from zero, so you can't limit its % change
		IF(liq==1.and.DABS(CHNG/etaOld).GT.0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==0.and.DABS(CHNG).GT.0.02d0)CHNG=DSIGN(0.02D0,CHNG)
		IF(liq==3.and.DABS(CHNG/etaOld).GT.0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==2.and.DABS(CHNG).GT.0.02d0)CHNG=DSIGN(0.02D0,CHNG)
 		eta=eta-CHNG
		if(eta < 0 .or. eta > 1/1.9)eta=etaOld-DSIGN(0.1D0,CHNG)*etaOld
		rho=eta/bMix
		IF(DABS(CHNG) < 1.D-9 .and. eta > 0)EXIT  ! Don't use CHNG/eta here. Converge quickly to ideal gas if P->0, still ~9 sigfigs if liquid.
	enddo !nIter=1,itMax
	if(eta < 0 .or. eta > 1.9)then
		if(LOUD)pause 'eta < 0 or > 1.9 on final iteration in FugiESD.'
		continue
	endif

	if(iera(4).ne.0.or.nIter.ge.itMax.or.eta < 0)then
		ier(4)=iera(4)
		if(LOUD)print*, 'FugiEsd: calling GoldenZ. nIter, eta',nIter,eta
		call GoldenZesd(X,nC,KVE,ND,bMix,NAS,NDS,cbMix,qYbMix,k1YbMix,Pb_RT,LIQ,zFactor,rho,zAssoc,XA,RALPH,fAssoc,ier)
		eta=rho*bMix
		ier(11)=1
		ier(1)=10
	endif
	!One last call of AlpSol to get derivative properties. Should not affect answer.
	isZiter=0
	CALL AlpSolEz2(isZiter,T,X,Nc,bMix,rho,zAssoc,aAssoc,uAssoc,chemPoAssoc,iera)
	IF (rho < 0)THEN
		ier(5)=1
        ier(1)=15
		if(LOUD)WRITE(6,31)LIQ
		GOTO 86
	ENDIF
!  ITERATION ON rho HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.

	if(ABS(eta-rho*bMix) > 1E-11 .and. LOUD)pause 'eta.ne.(rho*bMix)?'
	DO I=1,nC
	   YQVI(I)=0.D0
	   CVI(I)=0.D0
	enddo
	if(P==0 .and. LOUD)print*,'FugiEsd: P=0? LIQ,P=',LIQ,P
	zFactor=P/(rho*rGas*T)  ! add this to improve precision when computing rho from Z after return.   
	if(zFactor.le.0 .and. LOUD)print*,'FugiEsd: converged Z <= 0. eta,Z=',eta,zFactor

	!For uAssoc, cf M&H,FPE,180:165,2001. Eq. 6 can be written
	!Q({XAi},{alphai}) = sum(ni*[lnXAi+(1-XAi)+lnXDi+(1-XDi)]) - h({XAi},{alphai})/2 = aAssoc/RT
	!dQ/dBeps=sum(dQ/dXAi*dXAi/dBeps)+sum(dQ/dAlphai*dAlphai/dBeps)
	!By the stationary property, everything cancels except the derivatives of h wrt {alphai}, so
	!uAssoc/RT= -(1/2)sum{ xixjNdiNdjXAiXDj*beta*[d(rho*gKADij*Yij/dBeta] }
	!         = -(1/2)sum{ xixjNdiNdjXAiXDj*rho*gKADij*bepsADij*(Yij+1) }
	!         = -(1/2)sum{ xixjNdiNdjXAiXDj*alphaij*bepsADij*(Yij+1)/Yij }
	!uAssoc/RT= -(1/2)sum{ xixjNdiNdjXAiXDj*beta*[sqrt(alphai)*dSqrt(alphaj)/dBeta+sqrt(alphai)*dSqrt(alphaj)/dBeta] }
	!		  = -(1/2)sum{ xixjNdiNdjXAiXDj*beta*[sqrt(alphai)*(1/2)alphaj^(-0.5)*epsADi*(Yi+1)+sqrt(alphai)*dSqrt(alphaj)/dBeta] }
	!         = -(1/2)sum{ xixjNdiNdjXAiXDj*alphaij*beta*[(dSqrt(alphaj)/dBeta)/sqrt(alphaj)+(dSqrt(alphai)/dBeta)/sqrt(alphai)] }
	![beta*dSqrt(alphai)/dBeta]/sqrt(alphai)= [(+1/2)/sqrt(alphai)*rho*gKad*bepsADi*(Yi+1)]/sqrt(alphai)=(+1/2)*bepsADi*(Yi+1)/Yi
	!uAssoc/RT= -(1/2)sum{ xi*xjNdjXDjNdiXAi*alphaij* [bepsADi*(Yi+1)/Yi+bepsADj*(Yj+1)/Yj]/2 } ; 
	!         = -(1/2)sum{ xi*xjNdjXDjNdiXAi*sqrt(alphai)*sqrt(alphaj)*[bepsADi*(Yi+1)/Yi+bepsADj*(Yj+1)/Yj]/2 } ; 
	!         = -(1/4)sum{ xiNdiXAi*sqrt(alphai)*xjNdjXDj*sqrt(alphaj)*[bepsADi*(Yi+1)/Yi+bepsADj*(Yj+1)/Yj] } ; 
	!         = -(1/4)sum{ xiNdiXAi*sqrt(alphai)*xjNdjXDj*sqrt(alphaj)*[bepsADi*(Yi+1)/Yi]+sum{ xiNdiXAi*sqrt(alphai)*xjNdjXDj*sqrt(alphaj)*[bepsADj*(Yj+1)/Yj] } ; 2 in 2/4 is by symmetry.
	!         = -(1/4)sum{ xiNdiXAi*sqrt(alphai)*[bepsADi*(Yi+1)/Yi]*sum[xjNdjXDj*sqrt(alphaj)]+sum{ xiNdiXAi*sqrt(alphai)*sum(xjNdjXDj*sqrt(alphaj)*[bepsADj*(Yj+1)/Yj]) } ; 2 in 2/4 is by symmetry.
	!         = -(2/4)sum{ xiNdiXAi*sqrt(alphai)*bepsADi*(Yi+1)/Yi*sum(xjNdjXDjSqrt(alphaj) } ; 2 in 2/4 is by symmetry.
	!         = -(F/2)sum{ xiNdiXAi*Sqrt(alphai)*bepsADi*(Yi+1)/Yi } = (2F/2)*beta*dF/dBeta; beta*dF/dBeta=(+1/2)sum{ xiNdiXAi*Sqrt(alphai)*bepsADi*(Yi+1)/Yi }
	!CvAssoc/R = d(UAssoc/R)/dT. T*d(UAssoc/RT)/dT = (T/T)*d(UAssoc/R)/dT-T*(UAssoc/R)/T^2=>CvAssoc=uAssoc/RT+T*d(UAssoc/RT)/dT
	!CvAssoc=uAssoc/RT-beps*d(UAssoc/RT)/dBeps=uAssoc/RT-(1/2)(beps*dF/dBeps)^2-F/2*(dF/dBeps+beps^2*d2F/dBeps^2)
	!dF/dBeta = (1/2)*sum{ xiNdiXAiSqrt(rho*gKad)*epsADi*(Yi+1)/sqrt(Yi) }
	!d2F/dBeta2 = (1/2)*sum{ xiNdiXAiSqrt(rho*gKadi)*epsADi*[epsADi*(Yi+1)/sqrt(Yi)-(1/2)(Yi+1)/Yi^1.5*epsADi*(Yi+1) }
	!d2F/dBeta2 = (1/2)*sum{ xiNdiXAiSqrt(alphai)*epsADi^2*((Yi+1)/Yi)*[1-(1+Yi)/(2Yi)] }
	!beta^2d2F/dBeta2 = (1/2)*sum{ xiNdiXAiSqrt(alphai)*bepsADi^2*((Yi+1)/Yi)*[1-(1+Yi)/(2Yi)] }
	!beta^2d2F/dBeta2 = (1/2)*sum{ xiNdiXAiSqrt(alphai)*(bepsADi*(Yi+1)/Yi)*[bepsADi-bepsADi*(1+Yi)/(2Yi)] }
	!  d(...)B is beta*d(...)/dbeta.  d2(...)BB is d2(...)/dbeta^2, eTC...
	dYbB=0
	dYqbB=0
	DO I=1,nC
		do jComp=1,nC
			qbij=( Q(I)*Vx(jComp)+Q(jComp)*Vx(I) )/2
			kIJtemp=(KIJ(I,jComp)+KTIJ(I,jComp)/T) 
			eokIJroot=SQRT(eok(I,I)*eok(jComp,jComp))
			dKijB=KTIJ(I,jComp)/T
			dYqbB=dYqbB+X(I)*X(jComp)*qbij*(Y(I,jComp)+k2)*eokIJroot/T*( 1-kIJtemp - dKijB )
		enddo
        DO J=1,nC
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*X(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*X(J)
		enddo
	enddo

	dYbB = dYbB ! + dXsYbB !for xs form
	uAtt = -zM*qYbMix/(k1YbMix)*LOG(1+rho*k1YbMix)*( dYqbB/qYbMix - k10*dYbB/k1YbMix ) 
	uAtt =uAtt-zM*qYbMix/(k1YbMix)*k10*rho*dYbB/(1+rho*k1YbMix)

	eta=rho*bMix
	IF (LIQ.EQ.1) THEN
	  etaL=eta
	  ZL=zFactor  
	ELSE
	  etaV=eta
	  ZV=zFactor  
	ENDIF

	!call acTCoeff(fugXS,bVol,X,T,rho,nC)

	voidFrac=1-1.9D0*eta
	IF (voidFrac.LT.0)THEN
		if(LOUD)WRITE(6,*)'WARNING! (1-1.9*eta) IS -VE IN FUGI. eta=',eta
		ier(7)=1
        ier(1)=17
		GOTO 86
	ENDIF

	!etaYmix=qYbMix/qMix*rho					   !90 form, not used otherwise
	FREP= -4.D0/1.9D0*LOG(voidFrac)*cbMix/bMix
	FATT= -zM*qYbMix/k1YbMix*LOG(1.D0+k1YbMix*rho)
	uRes_RT=uAtt+UASSOC
	aRes_RT=FREP+FATT+aASSOC !-LOG(zFactor)
	DHONKT=uRes_RT+zFactor-1
	DUONKT=uAtt+UASSOC
	DAONKT=FREP+FATT+aASSOC !-LOG(zFactor)
	DSONK =DUONKT-DAONKT
	DGONKT=DAONKT+zFactor-1  !=sum[xi*fugc(i)]. note that bVol(i)/bMix -> 1 for pure i, so these terms sum to zRep+zAtt=Z-1.
	gAssoc=aAssoc+zAssoc
	if(initial.and.LOUD)write(*,'(a, I3,4F9.5)') ' FugiEsd2: eta done. LIQ,Z,zAssoc=',LIQ,zFactor,zAssoc !,XA(1)
	if(initial.and.LOUD)write(*,'(a, 4F9.5)') ' FugiEsd2: aResRep,aResAtt,aAssoc,aRes_RT:',frep,fatt,aAssoc,aRes_RT
	if(initial > 0)initial=initial+1
	if(initial > 3)initial=0										!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    initial   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(liq > 1)return !for EAR method. Don't evaluate LOG(zFactor) in case zFactor < 0.
	if(zFactor.le.0 .and. LOUD)pause 'FugiEsd: LIQ<2 but Z <= 0???'
	IF(LOUDER)WRITE(*,'(a,f10.5)')' kComp,lnGamRep,lnGamAtt,lnGamAssoc'
	DO kComp=1,nC
		FUGREP(kComp)=FREP*( 2*CVI(kComp)/cbMix-Vx(kComp)/bMix )+zRep*Vx(kComp)/bMix
		!dNk1YbNk=k1(kComp)*etaYmix*(2*YQVI(kComp)/qYbMix-Q(kComp)/qMix) !90 form, overwritten by next line for xs option
		dNk1YbNk=k1(kComp)*( Y(kComp,kComp)*Vx(kComp)+dXsYbN(kComp) ) !91 form,comment out to eliminate, cf. EL99 p559 & S&E91apx
		FUGATT(kComp)=zAtt*dNk1YbNk/k1YbMix+FATT*( 2*YQVI(kComp)/qYbMix-dNk1YbNk/k1YbMix )
		!FUGASN=-fAssoc*fAssoc*Vx(kComp)/bMix*1.9D0*eta/voidFrac
		!FUGAssoc(kComp)=2*ND(kComp)*LOG( XA(kComp) )+FUGASN
		!chemPoAssoc(kComp)=FUGAssoc(kComp)
		IF (zFactor<0 .and. LOUD)WRITE(6,*)'WARNING! Z NEGATIVE IN FUGI!'
		FUGC(kComp)=FUGREP(kComp)+FUGATT(kComp)+chemPoAssoc(kComp)-LOG(zFactor)
		!JRE: Pause here to check sample results.
		!rLnGamRep=FUGREP-Zrep*Vx(kComp)/bMix  ! cf. Bala and Lira (2016), Eqs A6-A14.
		!rLnGamAtt=FUGATT-Zatt*Vx(kComp)/bMix
		!rLnGamBon=FUGBON-Zassoc*Vx(kComp)/bMix
		!IF(LOUDER)WRITE(*,'(i3,5E12.4)')kComp,rLnGamRep,rLnGamAtt,rLnGamBon,ralph(kComp)
	enddo
	aRes_RT = DAONKT
	uRes_RT = DUONKT
	hRes_RT = DHONKT
	sRes_R  = DSONK											   
	gDep_RT = DGONKT
	cmprsblty=0
	cvRes_R = 0
	cpRes_R = 0
    if(LOUD .and. initial > 0)print*,'FugiEsd: returning. fugc=',(FUGC(i),i=1,NC)
	if(Nc>1)return	!I'm not confident in these derivatives for multicomp as of 9/21/19.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Derivative Props!!!!!!!!!!!!!!!!!!!!!!
	!cmprsblty=[(dP/dRho)/RT] by definition here.
	!Zassoc=-F^2/(1-1.9eta)=>eta*dZassoc/dEta= -2F/(1-1.9eta)*eta*dF/dEta-F^2*1.9eta/(1-1.9eta)^2
	!F = Nd*XA*ralph => eta*dF/dEta = F*[(eta/XA)*dXA/dEta+(eta/ralph)*dRalph/dEta]
	!1-XA = alpha*XA^2 => -eta*dXA/dEta = alpha*2XA*eta*dXA/dEta + XA^2*eta*dAlpha/dEta
	! -alphaXA^2*(eta/alpha)*dAlpha/dEta = eta*dXA/dEta*(1+2alpha*XA) = -(1-XA)*(eta/alpha)*dAlpha/dEta
	! (eta/XA)*dXA/dEta = -(1-XA)/(X+2alpha*XA^2)*(eta/alpha)*dAlpha/dEta= -(1-XA)/(2-XA)*(eta/alpha)*dAlpha/dEta
	! eta*dF/dEta = F * (eta/alpha)*(dAlpha/dEta)*[1/2-(1-XA)/(2-XA)]
	! alpha=eta*Kad'*Y/(1-1.9eta) =>(eta/alpha)*(dAlpha/dEta)=1/(1-1.9eta)
	! eta*dF/dEta = F * [1/2-(1-XA)/(2-XA)]/(1-1.9eta)
	! eta*dZassoc/dEta = -2F/(1-1.9eta)*{F * [1/2-(1-XA)/(2-XA)]/(1-1.9eta)} -F^2*1.9eta/(1-1.9eta)^2
	! eta*dZassoc/dEta = -F^2/(1-1.9eta)*{2*[1/2-(1-XA)/(2-XA)]/(1-1.9eta)} +zAssoc*1.9eta/(1-1.9eta)
	! eta*dZassoc/dEta = zAssoc*[1-2*(1-XA)/(2-XA)]/(1-1.9eta)} +zAssoc*1.9eta/(1-1.9eta)
	! 2*[1/2-(1-XA)/(2-XA)] = 2[(2-XA)-2(1-XA)]/(2(2-XA))= XA/(2-XA)
	rhoDZASSOC_dRho=0 !zAssoc/voidFrac*( 1.9*eta+XA(1)/(2-XA(1)) ) 
	cmprsblty=zFactor+zRep/voidFrac+zAtt/(1+k1Ybmix*rho)+rhoDZassoc_dRho
	beps=eok(1,1)/T
	!Aassoc=2*lnXA+(1-XA)=>uAssoc=(2-XA)(beps/XA)*dXA/dBeps;1-XA=alpha*XA^2=> -dXA/dBeps=2*alpha*XA*dXA/dBeps+XA^2*dAlpha/dBeps
	!(1+2alpha*X)(beps*dXA/dBeps) = -alphaXA^2(beps/alpha)*dAlpha/dBeps = -(1-XA)*beps(Y+1)/Y=(X+2alphaXA^2)(beps/XA)(dXA/dBeps)=(2-XA)(beps/XA)(dXA/dBeps)
	!(beps/XA)(dXA/dBeps)= -[(1-XA)/(2-XA)]*beps(Y+1)/Y => uAssoc= -(2-XA)[(1-XA)/(2-XA)]*beps(Y+1)/Y = -(1-XA)*beps(Y+1)/Y
	!uAssoc= -(1-XA(1))*bepsY1_Yad	 !multicomp formula above gives same result.
	!bepsDxa = -XA(1)*XA(1)*alpha*bepsDYad_Yad/SQRT(1+4*alpha)
	!CvAssoc/R = d(UAssoc/R)/dT. T*d(UAssoc/RT)/dT = (T/T)*d(UAssoc/R)/dT-T*(UAssoc/R)/T^2=>CvAssoc=uAssoc/RT+T*d(UAssoc/RT)/dT
	!CvAssoc=uAssoc/RT-beps*d(UAssoc/RT)/dBeps
	!beps*d(UAssoc/RT)/dBeps = -bepsY1_Yad*(-beps*dXA/dBeps) -(1-XA)*beps*(Y1/Yad+bepsY1/Yad-bepsY1^2/Y^2)
	!beps*d(UAssoc/RT)/dBeps = +bepsY1_Yad*(beps*dXA/dBeps) -(1-XA)*beps*(Y1/Yad)(1+beps-bepsY1/Y)
	!beps*d(UAssoc/RT)/dBeps = +bepsY1_Yad*[-XA*(1-XA)/(2-XA)]*bepsY1_Yad -(1-XA)*beps*(Y1/Yad)(1+beps-bepsY1/Y)
	!beps*d(UAssoc/RT)/dBeps = -bepsY1_Yad*[XA*(1-XA)/(2-XA)]*bepsY1_Yad -(1-XA)*beps*(Y1/Yad)(1+beps-bepsY1/Y)
	!beps*d(UAssoc/RT)/dBeps = -bepsY1_Yad^2*[XA*(1-XA)/(2-XA)] -(1-XA)*beps*(Y1/Yad)(1+beps-bepsY1/Y)
	bepsAD=DH(1)/T*TC(1)
	alpha=ralph(1)*ralph(1)
	bepsY1_Yad=1
	if(bepsAD > 1D-5)bepsY1_Yad=bepsAD*exp(bepsAD)/( exp(bepsAD)-1 ) !
	!CvAssoc= uAssoc + bepsY1_Yad**2*XA(1)*(1-XA(1))/(2-XA(1)) + (1-XA(1))*bepsY1_Yad*(1+bepsAD-bepsY1_Yad)
	!print*,'uAssoc,CvAssoc',uAssoc,CvAssoc
	cvRes_R = uAtt*( 1.7745*eta*beps*exp(beps)/(1+k1YbMix*rho) - beps ) !+ CvAssoc
	!Z=P/rho*R*T => dZ/dT = (dP/dT)/rho*R*T- P/(rho*RT^2)=> (dP/dT) = rho*R*(Z+T*dZ/dT)
	!V*dP/dV = -rho*dP/dRho => dP/dV = -rho^2*dP/dRho . 
	!cpRes_R =  CvRes_R -1 - (T/R)*(dP/dT)^2/(dP/dV) = CvRes_R-1 + (rho*R)^2*(Z+T*dZ/dT)^2*(T/R)/[rho^2*dP/dRho]
	!CpRes_R =  CvRes_R-1 + R*T*(Z+T*dZ/dT)^2/[dP/dRho] =  CvRes_R-1 + (Z-beps*dZ/dBeps)^2/[(dP/dRho)/RT]
    if(LOUD)then
	    if(exp(beps) < 1.0617)pause 'FuEsd96: beps < ln(1.0617)'
    end if
	bepsDY_Y=beps*exp(beps)/( exp(beps)-1.0617 ) !add 1D-8 to avoid zero divide.
	bepsDZatt_dBeps = zAtt*bepsDY_Y/(1+k1YbMix*rho)
		!zAssoc= -F^2/(1-1.9eta)=> beps*dZAssoc/dBeps = -2F*(beps*dF/dBeps)/denom
		!F = NdiXAiRalphi)=> beps*dF/dBeps = Ndi*ralphi*XA[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}= F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		!beps*dZAssoc/dBeps = -2F/denom *  F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]} = 2Zassoc*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		!bxDxa_db= -(1-XA(1))/(2-XA(1))*bepsY1_Yad
		bepsDZassoc_dBeps=0 !2*zAssoc*(bxDxa_db+bepsY1_Yad/2)
	!bepsDZassoc_dBeps = 2*zAssoc*bepsY1_Yad*( 0.5d0+fAssoc*SQRT(alpha)/SQRT(1+4*alpha) )
	TdZATT_dT=-ZM*RHO*( TdYQVM_dT*(1.d0+k1Ybmix*RHO)-qYbmix*TdK1YVM_dT*RHO )/( (1.d0+k1Ybmix*RHO)*(1.d0+k1Ybmix*RHO) )
	!print*,'FuEsd96: bepsDZatt_dBeps,TdZatt_dT',bepsDZatt_dBeps !,TdZatt_dT
	CpRes_R = cvRes_R-1 + (zFactor-bepsDZatt_dBeps-bepsDZassoc_dBeps)**2/cmprsblty 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Derivative Props!!!!!!!!!!!!!!!!!!!!!!

	initial=0
	RETURN
86	WRITE(6,*)' ERROR IN FuEsd96.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(NITER.GE.ITMAX)THEN
		if(LOUD)write(*,*)'TOO MANY Z ITERATIONS'
		ier(6)=1
        ier(1)=16
	END IF
	IF(ier(4)==1.and.LOUD) WRITE(*,*)'ERROR IN ALPSOL'
	if(ier(1) < 10)ier(1)=11
	initial=0
	RETURN
	END	!FugiESD()
	          
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!  ELLIOTT'S SUBROUTINE 
	!  20220702 Initial adaptation of grand master eq.
	SUBROUTINE AlpSolEz2(isZiter,tKelvin,X,Nc,VM,rho,zAssoc,aAssoc,uAssoc,rLnPhiAssoc,ier)
	USE GlobConst
	USE EsdParms
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (fAssoc) AND zAssoc given rho,VM
	Implicit DoublePrecision(A-H,K,O-Z)
	!PARAMETER(NMX=55)
	DoublePrecision X(NMX),RALPHA(NMX),RALPHD(NMX),XA(NMX),XD(NMX),rLnPhiAssoc(NMX)
	DoublePrecision KVEA(NMX),KVED(NMX) !,ralphTest(NMX)	!store and reuse values until ID's change.
	dimension xOld(NMX),IDold(NMX) !,KVE(NMX)
	!DoublePrecision ralphMatch(NMX)
    Integer ier(12)	 !NAS(NMX),NDS(NMX),ND(NMX),
	LOGICAL LOUDER
	!data xOld/0,0,0,0,0,0,0,0,0,0,0/ 
	!  NDi     = THE DEGREE OF POLYMERIZATION OF COMPO i
	!  fAssoc  = THE CHARACTERISTIC ASSOCIATION = 1/RALPHi(1/XAi-1)
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei/(1-1.9eta)
	!  ier(1)     = 1   Warning: FA iteration failed to converge.
	!  ier(1)     = 11  ralph calculation failed.
	!  ier(1)     = 14  eta > 0.53
	LOUDER=LOUD	! from GlobConst
	!LOUDER = .TRUE. ! for local debugging.
	picard=0.83D0

	eta=rho*VM
	voidFrac=1-1.9d0*eta
	rdfContact=1/voidFrac
	zAssoc=0
	aAssoc=0
	uAssoc=0
	rLnPhiAssoc(1:Nc)=0
	ier=0 !vector initialization.
	if(LOUDER)print*,'eta,rho,x1=',eta,rho,x(1)

	if(  SUM( ABS(ID(1:NC)-IDold(1:NC)) )>0  .or. ABS(tKelvin-tOld) > zeroTol  )then ! use stored values of KVE unless IDs change.
		DO I=1,NC           
			KVEA(I)=KadNm3(I)*AvoNum*( EXP(epsA_kB(I)/tKelvin)-1 )
			KVED(I)=KadNm3(I)*AvoNum*( EXP(epsD_kB(I)/tKelvin)-1 )
		ENDDO
	endif

	do i=1,NC
		RALPHA(I)=DSQRT( rho*rdfContact*KVEA(i) )
		RALPHD(I)=DSQRT( rho*rdfContact*KVED(i) )
		if(LOUDER)write(*,'(a,i3,2F12.4)')' i,ralphA,ralphD=',i,ralphA(i),ralphD(i)
	enddo
	FA0=SUM( x(1:Nc)*NDS(1:Nc)*ralphD(1:Nc) )                       
	FD0=SUM( x(1:Nc)*NAS(1:Nc)*ralphA(1:Nc) )
	if(  sumxmy2(NC,X,xOld) > 1.D-5 .or. SUM( ABS(ID(1:NC)-IDold(1:NC)) )>0  )then ! keep old guesses for ralphMeans unless	composition changes.
		avgNAS=SUM(x(1:NC)*NAS(1:NC))
 		avgNDS=SUM(x(1:NC)*NDS(1:NC))
		moreDonors=1
		if(avgNAS > avgNDS)moreDonors=0
		ralphAmean=( MAXVAL(ralphA) )  ! MAXVAL is perfect if only one compound with A sites.
		ralphDmean=( MAXVAL(ralphD) )
		!ralphTest(1:NC) = ralphA(1:NC)-ralphAmean
		if(avgNAS==0 .or. avgNDS==0)return ! that's all folks!
		if(FD0/avgNAS /= ralphAmean)ralphAmean=FD0/avgNAS ! if only one compound with A sites, then avg=MAX. 
		if(FA0/avgNDS /= ralphDmean)ralphDmean=FA0/avgNDS ! if only one compound with D sites, then avg=MAX. 
		if(moreDonors)then
			ralphDmean=ralphDmean*avgNAS/avgNDS
		else
			ralphAmean=ralphAmean*avgNDS/avgNAS
		endif
	else
		if(etaOld>0)then ! adapt old values.to accelerate Z iterations.
			ralphAmean=ralphAmean*SQRT(eta/etaOld*rdfContact/rdfOld)
			ralphDmean=ralphDmean*SQRT(eta/etaOld*rdfContact/rdfOld)
		endif
	endif
	if(FA0==0 .or. FD0==0)return ! no need to calculate solvation if either no donors or no acceptors
	if(LOUDER)then
		do i=1,NC
			print*,'i,NAS,NDS',i,NAS(i),NDS(i)
		enddo
	endif

	FDold = 0       ! low guess
	errOld= -FD0	! sumA= FD0 when FD=0.
	FAold = 0
	errOld= -FA0
		                       
	NITER=0
	ITMAX=44
	avgFo=(ralphDmean*FD0+ralphAmean*FA0)/2
	delFo=(ralphDmean*FD0-ralphAmean*FA0)
	delFA=1+delFo
	FA = 2*FA0/( delFA+DSQRT(delFA*delFA+4*ralphAmean*FA0) )
	delFD=1-delFo
	FD = 2*FD0/( delFD+DSQRT(delFD*delFD+4*ralphDmean*FD0) )	! How good is this quadratic estimate vs. secant?
	if(ABS(delFo/avgFo)	< 0.1D0)picard=1.1 ! accelerate for symmetric solvation.
	error=1234 !must initialize before while.
	if(LOUDER)write(*,'(a,i3,2f8.3,2f12.6,4E12.4)')' ralphAmean,ralphDmean,FA,FD,FA0,FD0',nIter,ralphAmean,ralphDmean,FA,FD,FA0,FD0
	Ftol=1.D-7
	do while(ABS(error)>Ftol.and.ABS(FA-FAold)>Ftol.and.NITER<ITMAX)
		NITER=NITER+1
		!Solve quadratic
		!delFA=1+(ralphDmean*FD0-ralphAmean*FA0)
		!FA = 2*FA0/( delFA+DSQRT(delFA*delFA+4*ralphAmean*FA0) )
		delFo=(ralphDmean*FD0-ralphAmean*FA0) !FYI:I tried setting FD=FA when delFo<Ftol, but it made things worse...?
		delFD=1-delFo
		FD = 2*FD0/( delFD+DSQRT(delFD*delFD+4*ralphDmean*FD0) )
		delFA=1+delFo
		FA = 2*FA0/( delFA+DSQRT(delFA*delFA+4*ralphAmean*FA0) ) ! This quadratic estimate is much better than secant.
		sumA=0
		do j=1,Nc
			sumA=sumA+x(j)*NAS(j)*ralphA(j)/( 1+FA*ralphA(j) ) ! using best current estimate of FA=sumD(FD)
		enddo
		errD=(FD-sumA) ! If FD is right, then ralphMeans are right and FD=sumA. 
		sumD=0
		do j=1,Nc
			sumD=sumD+x(j)*NDS(j)*ralphD(j)/( 1+FD*ralphD(j) )
		enddo
		errA=FA-sumD
		error=MAX(ABS(errA),ABS(errD))
		!Prepare for next iteration.
		if(LOUDER)write(*,'(a,i3,2f8.3,2F12.6,4E12.4)')' ralphAmean,ralphDmean,FA,FD,erA,erD',nIter,ralphAmean,ralphDmean,FA,FD,errA,errD
		ralphAold=ralphAmean !the point of this loop is to compute the new values and compare to old.
		ralphDold=ralphDmean
		ralphDmean= picard*(-1+FA0/sumD)/(FD+1D-9)+(1-picard)*ralphDmean !Using new sumD to compute new ralphMeans.
		ralphAmean= picard*(-1+FD0/sumA)/(FA+1D-9)+(1-picard)*ralphAmean !add 1D-9 to avoid possible zero divide if FD->0
		if(NITER==1)FA1=FA
		if(NITER==1)FD1=FD
		!set default in case error~errOld crashes. This uses the new ralphMeans so should be ~improvement.
		!delFD=1-(ralphDmean*FD0-ralphAmean*FA0)
		!FD = 2*FD0/( delFD+DSQRT(delFD*delFD+4*ralphDmean*FD0) )	
		!Store current values as old values before updating FA.
	enddo !while
	if(NITER>ITMAX-1)then
		if(ABS(FA-FAold) > Ftol)ier(1)=1 ! Leave this as a warning.
		if(LOUDER)write(*,*)'ERROR - NO CNVRG. NITER,FA1,FA=',NITER,FA1,FA
		!GOTO 86
	ENDIF
	if(LOUDER)write(*,'(a,6E12.4,2i6)')' AlpsolEz2: errA,errD',errA,errD
	IF(NITER > itMax-1.and.eta<0.43)write(*,'(a,3i5,1x,3F8.4,2E12.4)')' AlpsolEz2: NITER,T,x1,eta,erA,erD=',NITER,ID(1),ID(2),tKelvin,x(1),eta,errA,errD
	!write(*,'(a,3i5,1x,3F8.4,2E12.4)')' AlpsolEz2: NITER,T,x1,eta,erA,erD=',NITER,ID(1),ID(2),tKelvin,x(1),eta,errA,errD
	!  fAssoc ITERATION HAS CONCLUDED
	tOld=tKelvin
	xOld(1:NC)=x(1:NC)
	IDold(1:NC)=ID(1:NC)
	etaOld=eta
	rdfOld=rdfContact ! from computation of ralphA(),ralphD()
	zAssoc= -FA*FD/voidFrac
	if(isZiter==1)return
	uFA= 0
	uFD= 0
	aFA=0
	aFD=0
	Abonds=0
	Dbonds=0
	DO I=1,nC
		XA(I)=1/(FA*RALPHA(I)+1)
		XD(I)=1/(FD*RALPHD(I)+1)
		rLnPhiAssoc(i)=NAS(I)*LOG( XA(I) )+NDS(I)*LOG( XD(I) )+zAssoc*1.9d0*bVolCC_mol(i)*rho !ESD convention is one bonding type per molecule.	Use SPEADMD for more.
		Abonds=Abonds+x(i)*NAS(I)*(1-XA(I))
		Dbonds=Dbonds+x(i)*NDS(I)*(1-XD(I))
		aFA=aFA+x(i)*NAS(i)*(  DLOG( XA(i) )+( 1-XA(i) )/2  )
		aFD=aFD+x(i)*NDS(i)*(  DLOG( XD(i) )+( 1-XD(i) )/2  )
		bepsD=epsD_kB(i)/tKelvin
		eBepsD=bepsD*EXP( bepsD )
		if(bepsD > 1.d-5)eBepsD=bepsD*eBepsD/(eBepsD-1) ! Avoid divide by zero when bepsD->0. 
		bepsA=epsA_kB(i)/tKelvin
		eBepsA=bepsA*EXP( bepsA )
		if(bepsA > 1.d-5)eBepsA=bepsA*eBepsA/(eBepsA-1) ! Avoid divide by zero when bepsA->0. 
		uFA=uFA+x(i)*NDS(i)*XD(i)*ralphD(i)*eBepsD
		uFD=uFD+x(i)*NAS(i)*XA(i)*ralphA(i)*eBepsA
	enddo
	aAssoc=(aFA+aFD)
	uAssoc=(FD*uFA+FA*uFD)/2
	if(ABS(Abonds-Dbonds) > 1.d-5 .and. eta<0.43)then
		write(*,'(a,2E12.4,F8.2,2F8.5)')' Abonds,Dbonds,T,x1,eta',Abonds,Dbonds,tKelvin,x(1),eta
		pause 'Check bond site balance.'
	endif

	if(LOUDER)print*,'eta,FA,FD,Zassoc,ralph(),lnPhiAssoc'
	if(LOUDER)write(*,'(8f10.4)')eta,FA,FD,zAssoc,(ralphA(i),i=1,NC),(rLnPhiAssoc(i),i=1,NC)
86	CONTINUE
	RETURN
	END	 !AlpSolEz2()

!**************************************************************************
	SUBROUTINE GoldenZesd2(X,nC,KVE,ND,bMix,NAS,NDS,CbMix,YqbMix,k1YbMix,&
                       PoRT,LIQ,zFactor,rho,zAssoc,XA,RALPH,fAssoc,ier)
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION X(*),XA(*),KVE(*),RALPH(*),ier(*)
	dimension ND(*),NAS(*),NDS(*)
	data rgold,cgold,zM/0.61803399,0.38196602,9.5/
	etaLo = 1.e-5
	etaHi = 0.1
	if(liq.eq.1)then
		etaLo=0.15;
		etaHi=0.5;
	endif

	rhoLo = etaLo/bMix;
	rhoHi = etaHi/bMix;

	rho1=rhoLo+cgold*(rhoHi-rhoLo);
	rho2=rhoLo+rgold*(rhoHi-rhoLo);

	rho=eta/bMix
	CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rhoLo,zAssoc,XA,RALPH,fAssoc,ier)
	zRep=4*CbMix*rhoLo/(1-1.9D0*bMix*rhoLo)
	zAtt=-zM*YqbMix*rhoLo/(1+k1YbMix*rhoLo)
	zFactor=(1+zRep+zAtt+zAssoc)
	errLo = ABS(PoRT-rhoLo*zFactor);

	CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rho1,zAssoc,XA,RALPH,fAssoc,ier)
	zRep = 4 * cbMix * rho1 / (1 - 1.90*bMix*rho1);
	zAtt = -zM * YqbMix * rho1 / (1 + k1YbMix * rho1);
	zFactor = (1. + zRep + zAtt + zAssoc);
	err1 = ABS(PoRT-rho1*zFactor);

	CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rho2,zAssoc,XA,RALPH,fAssoc,ier)
	zRep = 4 * cbMix * rho2 / (1 - 1.90*bMix*rho2);
	zAtt = -zM * YqbMix * rho2 / (1 + k1YbMix * rho2);
	zFactor = (1 + zRep + zAtt + zAssoc);
	err2 = ABS(PoRT-rho2*zFactor);

	CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rhoHi,zAssoc,XA,RALPH,fAssoc,ier)
	zRep = 4 * cbMix * rhoHi / (1 - 1.90*bMix*rhoHi);
	zAtt = -zM * YqbMix * rhoHi / (1 + k1YbMix * rhoHi);
	zFactor = (1. + zRep + zAtt + zAssoc);
	errHi = ABS(PoRT-rhoHi*zFactor);

	do iter=1,111

		if(err2 < err1)then
			rhoLo = rho1;
			rho1  = rho2;
			rho2  = rgold*rho1+cgold*rhoHi;
			errLo = err1;
			err1  = err2;

			CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rho2,zAssoc,XA,RALPH,fAssoc,ier)
			zRep = 4 * cbMix * rho2 / (1 - 1.90*bMix*rho2)
			zAtt = -zM * YqbMix * rho2 / (1 + k1YbMix * rho2)
			zFactor = (1 + zRep + zAtt + zAssoc)
			err2 = ABS(PoRT-rho2*zFactor)

		else

			rhoHi = rho2;
			rho2  = rho1;
			rho1  = rgold*rho2+cgold*rhoLo;
			errHi = err2;
			err2  = err1;

			CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rho1,zAssoc,XA,RALPH,fAssoc,ier)
			zRep = 4 * cbMix * rho1 / (1 - 1.90*bMix*rho1);
			zAtt = -zM * YqbMix * rho1 / (1 + k1YbMix * rho1);
			zFactor = (1. + zRep + zAtt + zAssoc);
			err1 = ABS(PoRT-rho1*zFactor);

		endif
	enddo !iter=1,111

!//that's the best we can do for now folks. get props at best guess and return.
	rho = (rho1+rho2) * 0.5;
	zFactor = PoRT/rho;
!	err1=ABS(PoRT-rho*zFactor);
	eta = bMix*rho;
	return
	end	 !GoldenZb
    
	SUBROUTINE FugiEtaESD2(tKelvin,eta,X,Nc,LIQ,FUGC,zFactor,ier)
	!Purpose: get fugacity coeff and misc for activity coeff/GibbsXs at const eta
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	  + GlobConst{rGas,Tc,Pc,...}
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	dimension X(NMX),FUGC(NMX),ier(12) !,chemPoAssoc(nmx)
	!dimension YQVIJ(NMX,NMX),KVE(NMX),YQVI(NMX),Y(NMX,NMX),eok(NMX,NMX)
	!dimension CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX),XA(NMX),dXsYbN(NMX)
	!Dimension RALPH(NMX),k1(NMX)
	!COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!      COMMON/xsGibbs/xsBIP(NMX,NMX),xsNrtlAl(NMX,NMX)
	!  ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!  eokP IS THE DISPERSE ATTRACTION OVER BOLTZ k FOR PURE SPECIES
	!  KCSTAR IS THE DIMENSIONLESS BONDING VOLUME FOR PURE 
	!  DH IS THE BONDING ENERGY /RTC 
	!  C,Q,bVol ARE THE PURE COMPONENT EOS PARAMETERS
	!  KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!  Z IS PV/NoKT HERE  
	!  ier = 1 - AT LEAST ONE ERROR
	!        2 - NOT USED
	!        3 - NOT USED
	!        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!        5 - rho IS -VE
	!        6 - TOO MANY Z ITERATIONS
	!        7 - eta > 0.53
	!		11 - goldenZ instead of real Z.
	!
	!DATA k10,K2,zM,INITIAL/1.7745,1.0617,9.5,0/  !Zm=9.5 since 1996, at least.  It is 9.5 in the CHEMCAD documentation, the EL text, and 2002.  It is not mentioned in ref[2,3,4,5,6]
	DATA INITIAL/1/
	do i=1,12
		ier(i)=0
	enddo
	INITIAL=1
	if(LIQ==1.and.eta < 0.1 .and. LOUD)print*,'FugiEtaEsd2: eta < 0.1 for LIQ=1?'
	sumx=SUM(x(1:Nc))
	if(ABS(sumx-1) > 1e-8)then
		if(LOUD)pause 'FugiEsd2: sumx .ne. 1'
	endif
	bMix=SUM( bVolCC_mol(1:Nc)*x(1:Nc) )
	vTotCc=bMix/eta	 ! assume 1 mol basis

	CALL FuEsd2Vtot(isZiter,tKelvin,vTotCc,x,Nc,FUGC,zFactor,Ares,Ures,iErr)
	if(iErr > 0)ier(1)=1
	RETURN
86	WRITE(6,*)' ERROR IN FuEsd96.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(ier(4)==1 .and. LOUD) WRITE(*,*)'ERROR IN ALPSOL'
	if(ier(1) < 10)ier(1)=11
	RETURN
	END	  !FugiEtaEsd()
	          
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	Written by AFG, Oct. 2009																				C
!C	With Having T,V and nMols, this routine calculates the Z factor and all derivatives needed in 			C
!C	critical point,	flash and bubble point calculations.													C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuEsd2Vtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,Ares,Ures,iErr)
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(NMX),xFrac(NMX),FUGC(NMX),KCSTARp(NMX)
	integer iera(12) !IER(12),
	DoublePrecision YQVIJ(NMX,NMX),KVE(NMX),YQVI(NMX),Y(NMX,NMX),EOK(NMX,NMX)
	DoublePrecision CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX)
	DoublePrecision k1(NMX) !,RALPHA(NMX),RALPHD(NMX)
	!DoublePrecision FUGASSOC(NMX),FUGREP(NMX),FUGATT(NMX)			  
	!double precision moleStep
	DoublePrecision dEOK_dT(NMX,NMX),dY_dT(NMX,NMX),dYQVI_dT(NMX),dYQVIJ_dT(NMX,NMX)
	DoublePrecision dCVM_dN(NMX),dYQVM_dN(NMX) ,dVM_dN(NMX),dK1YVM_dN(NMX) !,dCVI_dN(NMX,NMX)
	!DIMENSION dZREP_dN(NMX),dZATT_dN(NMX),dZASSOC_dN(NMX),dZ_dN(NMX),dP_dN(NMX),dFREP_dN(NMX),dFATT_dN(NMX),dYQVI_dN(NMX,NMX)
	!DIMENSION dFUGC1_dP(NMX),dFUGC1_dT(NMX),dFUGC2_dP(NMX),dFUGC2_dT(NMX),dFUGASSOC_dRHO(NMX)
	!DIMENSION dFUGREP_dN(NMX,NMX),dFUGATT_dN(NMX,NMX),dFUGC1_dN(NMX,NMX),dFUGC2_dN(NMX,NMX),dFUG_dN(NMX,NMX)
	!DIMENSION dh_dN(NMX),dFUGASSOC_dT(NMX),dFUGASSOC_dN(NMX,NMX)
	!DIMENSION dh_dN_num(NMX),dFUGASSOC_dN_num(NMX,NMX),fugassocLoop(NMX),gmol_old(NMX)
	!DoublePrecision vMolecNm3(NMX) !for Wertheim
	LOGICAL LOUDER
	
	COMMON/ETA/ETAL,ETAV,ZL,ZV
	COMMON/ETA2/ETA
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	COMMON/ALPHA/ALPHAD(NMX,NMX), ALPHDA(NMX,NMX)
	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(nmx)
	COMMON/Derv1/dbVolMix_dN(NMX),DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,assocFlag
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta
	COMMON/dFug/h_nMichelsen,dh_dT,dh_dN,dFugassoc_dT,dFugassoc_dN,dFUGASSOC_dRHO,dh_dRHO!,FUGASSOC_num
	COMMON/dfug2/dFUGASSOC_dN_num,dh_dN_num
	COMMON/num/FUGASSOC_num
	COMMON/fugCR/PMpa,dFUG_dN,dP_dN
	COMMON/Helm_derv/Ur,UrV,FREP,FATT,FassocAfg
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),ralph(nmx),Zrep,Zatt,Zassoc,Fassoc
	DATA K10,K2,ZM,initCall/1.7745D0,1.0617D0,9.5D0,1/

	LOUDER=LOUD
	!LOUDER=.TRUE.
	  
	iErr=0 !1=warning from AlpSolEz2, 11=input nonsense, 12=xFrac<0, 13=critical error from AlpSol, 14=voidFrac<0..
	totMoles=sum(gmol)
	if( tKelvin	   < zeroTol .or. totMoles < zeroTol .or. vTotCc < zeroTol)then
		if(LOUD)print*,'FuEsdVtot: nonsense T(K),totMoles,vTotCc=',tKelvin,totMoles,vTotCc
		iErr=11
	endif
	bMix=0
	do i=1,nc
		xFrac(i)=gmol(i)/totMoles
		IF(xFrac(i) < 0)then
			if(LOUD)PAUSE 'ERROR IN FuEsdVtot, Xi<0'
			iErr=12
		endif
		bMix=bMix+xFrac(i)*bVolCC_mol(i)
	enddo
	if(iErr > 0)return
	rho=totMoles/ vTotCc
	eta=rho*bMix 
	if(LOUD.and.initCall)print*,'FuEsdVtot: bMix,eta=',bMix,eta
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
		dHkcalMol(I)=DH(I)*TC(I)*1.987D-3	!JRE 12/00
		KVE(I)=KCSTAR(I)*( DEXP(dHkcalMol(I)/tKelvin/1.987D-3)-1.d0 )  *avoNum  !KVE[=] cc/mol, KcStar [=] nm3.
		!	if(initCall.and.LOUD)print*,'FuESDVtot: KcStar=',KcStar(i)
		!if(initCall.and.LOUD)pause 'Right?'
		KCSTARp(I)=KCSTAR(I)
		VM=VM+xFrac(I)*VX(I)
		K1YVM=K1YVM+xFrac(I)*K1(I)*Y(I,I)*VX(I) !1991 form, overwritten if applying 1990 form
		!vMolecNm3(i)=VX(I)/avoNum
		bondVolNm3Esd(I)=KCSTAR(I) !*vMolecNm3(I)
		eHbKcal_mol(I,iType)=dHkcalMol(I)
		bondVolNm3(I,iType)=bondVolNm3Esd(I)
	enddo
	if( ABS(VX(1) - bVolCC_mol(1)) > zeroTol .and. LOUD ) print*,'FuEsdVtot: VX.ne.bVol=',VX(1),bVolCC_mol(1)
	if(LOUD.and.k1yvm < zeroTol)print*,'FuEsdVtot: 0~k1yvm=',k1yvm 
	eta=rho*vm
	CALL AlpSolEz2(isZiter,tKelvin,xFrac,Nc,VM,rho,zAssoc,aAssoc,uAssoc,FUGASSOC,iera)
	IF(iera(1) > 0)then
        if(iera(1)>10)iErr=13 ! This is a warning level for iErr when AlpSol did not converge.
		if(iera(1)==1)iErr=1
    endif
	voidFrac=1-1.9D0*ETA
	denom=voidFrac
	ZREP= 4.d0*CVM*RHO/denom  ! old form.
	ZREP= 4.d0*Cmix*eta/denom
	ZATT= -ZM*YQVM*RHO/(1+K1YVM*RHO)
	Z=(1+ZREP+ZATT+ZASSOC)
	PMpa=Z*rGas*RHO*tKelvin
    if(voidFrac < 0)then
	    IF(LOUD)print*, 'FuEsd2Vtot:Error! (1-1.9*ETA) IS -VE. eta,rho=',eta,rho
		iErr=14
    end if
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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Added by AFG 2010
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!JRE: compute these derivatives for thermal props, regardless of isZiter. They are returned as part of GlobalConst
	! Mixing Rule derivatives
	DO I=1,NC 
		dCVM_dN(I)=0.d0
		dYQVM_dN(I)=0.d0
		DO J=1,NC
	  		dEOK_dT(I,J)=DSQRT(EOKP(I)*EOKP(J))*(KTIJ(I,J)/(tKelvin*tKelvin)) 
			dY_dT(I,J)=(dEOK_dT(I,J)*tKelvin-EOK(I,J))*DEXP(EOK(I,J)/tKelvin)/(tKelvin*tKelvin)
			dCVM_dN(I)=dCVM_dN(I)+2.d0*xFrac(J)*CVIJ(I,J)
			dYQVM_dN(I)=dYQVM_dN(I)+2.d0*xFrac(J)*YQVIJ(I,J)
		enddo
	enddo
	DO I=1,NC
		dCVM_dN(I)=dCVM_dN(I)-2.d0*CVM
		dVM_dN(I)=VX(I)-VM
		dK1YVM_dN(I)=K1(I)*Y(I,I)*VX(I)-K1YVM
		dYQVM_dN(I)=dYQVM_dN(I)-2.d0*YQVM
	ENDDO
	dYQVM_dT=0.d0
	dK1YVM_dT=0.d0
	DO I=1,NC
		DO J=1,NC
			dYQVIJ_dT(I,J)=QV(I,J)*dY_dT(I,J)
			dYQVI_dT(I)=dYQVI_dT(I)+xFrac(J)*dYQVIJ_dT(I,J)
			dYQVM_dT=dYQVM_dT+xFrac(I)*xFrac(J)*QV(I,J)*dY_dT(I,J)
		ENDDO
		dK1YVM_dT=dK1YVM_dT+xFrac(I)*K1(I)*VX(I)*dY_dT(I,I)
	ENDDO
	!print*,'dK1YVM_dT,dYQVM_dT',dK1YVM_dT,dYQVM_dT

	! First derivatives of P,ZREP,ZATT and ZASSOC in respect to RHO according to the ESD EOS while T and N are constant
	voidFrac=(1.d0-1.9d0*VM*RHO)
	dZREP_dRHO=4.d0*CVM/(voidFrac*voidFrac) !=zRep/(1-1.9eta)/rho
	dZATT_dRHO=-ZM*YQVM/((1.d0+K1YVM*RHO)*(1.d0+K1YVM*RHO))		  !=zAtt/(1+k1Yeta)/rho
	dZASSOC_dRHO=-half*(h_nMichelsen*(dAlph_deta*VM)+dh_dRHO*(1.d0+dLng)) !(eta/g)*dg/dEta=eta*(1-1.9eta)*1.9/(1-1.9eta)^2=1.9eta/(1-1.9eta)
	dZASSOC_dRHO=-half*(2*fAssoc**2*(dAlph_deta*VM)+dh_dRHO*(1.d0+dLng)) !(eta/g)*dg/dEta=eta*(1-1.9eta)*1.9/(1-1.9eta)^2=1.9eta/(1-1.9eta)
	rhoDZAssoc_dRho2=rho*dZASSOC_dRHO
	rhoDZASSOC_dRho=-(-zAssoc*1.9*eta+XA(1,1)*(1-XA(1,1))/(2-XA(1,1))/voidFrac)/voidFrac !(eta/g)*dg/dEta=eta*(1-1.9eta)*1.9/(1-1.9eta)^2=1.9eta/(1-1.9eta)
	!print*,'rhoDZAssoc_dRho 1&2',rhoDZAssoc_dRho,rhoDZAssoc_dRho2
	dP_dRHO=rGas*tKelvin*Z+rGas*tKelvin*(RHO*dZREP_dRHO+RHO*dZATT_dRHO+rhoDZASSOC_dRho)!=RT*(Z+rho*dZ_dRho)
	!rho*dh_dRho= -(dXAdRho+dXDdRho)
	!rhoDh_dRho=rho*dh_dRHO
	!rhoDh_dRho2= 2*XA(1,1)*(1-XA(1,1))/(2-XA(1,1))/voidFrac
	!print*,'rhoDh_dRho 1&2',rhoDh_dRho,rhoDh_dRho2
	!rhoDZassoc_dRho3= 2*Zassoc*( half-(1-XA(1,1))/(2-XA(1,1)) )/voidFrac
	!print*,'h_nMichelsen/2,f^2',h_nMichelsen/2,fAssoc**2


	dP_dV=dP_dRHO*RHO*(-1.d0/ vTotCc)
	IF (dP_dRHO.EQ.0.d0) THEN
		dP_dRHO=1E-18
	ENDIF
	dRho_dP=1/dP_dRHO
	cmprsblty2=dP_dRHO/rGas/tKelvin
	cmprsblty=Z+zRep/denom+zAtt/(1+K1YVM*RHO) + rhoDZassoc_dRho !cf. EsdDerivatives.jnt
	!print*,'cmprsblty&2=',cmprsblty,cmprsblty2
	dZ_dP=(dZREP_dRHO+dZATT_dRHO+dZASSOC_dRHO)*dRho_dP

	! First derivatives of P,ZREP,ZATT and ZASSOC in respect to T according to the ESD EOS while V and N are constant
	dZREP_dT=0.d0
	dZATT_dT=-ZM*RHO*( dYQVM_dT*(1.d0+K1YVM*RHO)-RHO*YQVM*dK1YVM_dT )/( (1.d0+K1YVM*RHO)*(1.d0+K1YVM*RHO) )
	dZASSOC_dT=-half*(1.9d0*RHO*VM/denom+1)*dh_dT !dh_dT is computed in WertheimFugc.
	dP_dT=rGas*RHO*Z+rGas*tKelvin*RHO*(dZREP_dT+dZATT_dT+dZASSOC_dT)
	dT_dP=1/dP_dT
	dZ_dT=(dZREP_dT+dZATT_dT+dZASSOC_dT)
	cvRes_R=0
	if(NC==1)then	!added by JRE 2018. I'm not confident of formulas for NC>1 as of 20191004.
		beps=eok(1,1)/tKelvin
		bepsAD=DH(1)/tKelvin*TC(1)
		alpha=ralph(1)*ralph(1)
		bepsY1_Yad=1
		if(bepsAD > 1D-5)bepsY1_Yad=bepsAD*exp(bepsAD)/( exp(bepsAD)-1 ) !
		CvAssoc= uAssoc + bepsY1_Yad**2*XA(1,1)*(1-XA(1,1))/(2-XA(1,1)) + (1-XA(1,1))*bepsY1_Yad*(1+bepsAD-bepsY1_Yad)
		!print*,'uAssoc,CvAssoc=',uAssoc,CvAssoc
		cvRes_R = uAtt*( 1.7745*eta*beps*exp(beps)/(1+k1YVM*rho) - beps ) + CvAssoc
		!print*,'FuEsdVtot: TdZ_dT',tKelvin*dZ_dT
		!zAssoc= -F^2/(1-1.9eta)=> beps*dZAssoc/dBeps = -2F*(beps*dF/dBeps)/denom
		!F = NdiXAiRalphi)=> beps*dF/dBeps = Ndi*ralphi*XA[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}= F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		!beps*dZAssoc/dBeps = -2F/denom *  F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]} = 2Zassoc*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		bxDxa_db= -(1-XA(1,1))/(2-XA(1,1))*bepsY1_Yad
		bepsDZassoc_dBeps=2*zAssoc*(bxDxa_db+bepsY1_Yad/2)
		!write(666,'(a,5f10.4)')'beps,rho,zAssoc,rhoDZassoc_dRho',beps,zAssoc,bepsDZassoc_dBeps

		bepsDZatt_dBeps= -tKelvin*dZatt_dT !see derivative formulas above.
		TdP_dT1=Z+tKelvin*dZ_dT
		TdZassoc_dT=tKelvin*dZASSOC_dT
		TdP_dT=Z-bepsDZatt_dBeps-bepsDZassoc_dBeps
		!print*,'dZASSOC_dT,dh_dT ',dZASSOC_dT,dh_dT
		if(ABS(cmprsblty) < 1D-11)then
			if(LOUD.and.initCall)Print*,'FuEsdVtot warning: cmprsblty ~ 0'
			cpRes_R = 86.8686D0
		else
			cpRes_R = cvRes_R-1+TdP_dT**2/cmprsblty 
        end if
	endif
	if (isZiter==0) then
		!call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
        if( Z .le. -8686)then  ! Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
		    if(LOUDER)pause 'FuEsdVtot: Z.le.0 for fugc calculation.'
			iErr=3
			goto 86
        end if
		if(LOUDER)write(*,'(a,f10.5)')' i,lnGamRep,lnGamAtt,lnGamBon,ralph.  F=',Fassoc
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
86	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine RegPureEsd2(NC)!note:tc,pc,acen,id & esd parms passed in common. 
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (4/2020) to use CritPure for PGL package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)	 ! max # of EOS parameters to estimate: c,eok,b,Kad,epsAD
	DoublePrecision	parm(3),error(3),stderr(3)
	integer iErr2(NC) ! iErrExactEsd list for each component.
	LOGICAL LOUDER
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDev2
	LOUDER=LOUD
	LOUDER=.TRUE.
	if(NC > 1)then
		pause 'RegpureEsd2 only works for Nc=1. Sorry.'
		return
	endif
	print*,'Initial guess from ParmsEsd? Enter 1 for yes or 0 to use MW guides as guess'
	read(*,*)iAns
	if(iAns==0)then
		c(1)=1+rmw(1)/150
		eokP(1)=359+.66/rmw(1)
		Vx(1)=17*rmw(1)/50
		parm(1)=c(1)
		if(ND(1)==0)call ExactEsd2(nC,VX,c,q,eokP,iErr,iErr2)
		if(ND(1)==1)call RegPureDev3b(3,1,parm,error,iflag) !replaces eokP and Vx
	endif
	if(LOUDER)write(*,*)'RegPure1:  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!  general ***********************************
	q(1)=1+1.90476*(c(1)-1)
	nParms=3
	PARM(1)=c(1)
	PARM(2)=eokP(1)
	PARM(3)=Vx(1)
	call RegPureDev2(nData,nParms,parm,error,iflag)
	if(LOUDER)write(*,'(a,3f11.3,i3)')' C,Eok,vX',C(1),eokP(1),VX(1),Nd(1)
	if(LOUDER)write(*,*)'initial devA,devT,devP=',(error(i),i=1,3)
	aceCalc=acen(1)+error(1)
	TcCalc=Tc(1)+error(2)
	PcCalc=Pc(1)+error(3)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
	iflag=0
	iterOpt=0
	rguess=1234
	if(iteropt==0)then
		do while(rguess > 0)
			write(*,*)'enter C (-1 to skip to regression using last guess)'
			read(*,*)rguess,eokP(1),Vx(1)
			if(rguess.gt.0)then
				parm(1)=rguess
				if(ND(1)==1)call RegPureDev3b(3,1,parm,error,iflag) ! use old method as guess for eokP and Vx, given c
				parm(2)=eokP(1)
				parm(3)=vx(1)
				call RegPureDev2(nData,nParms,parm,error,iflag)
				if(LOUDER)write(*,'(a,1x,f9.4,f7.2,f9.4)')' initial devA,devT,devP=',(error(i),i=1,3)
				aceCalc=acen(1)+error(1)
				TcCalc=Tc(1)+error(2)
				PcCalc=Pc(1)+error(3)
				if(LOUDER)write(*,'(a,1x,f8.4,2f9.4,i3)')' c,eok,Vx=',(parm(i),i=1,3),Nd(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
			endif
		enddo
	endif
	q(1)=1+1.90476*(c(1)-1)
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!	pause
	factor=1.d-4
	tol=1.d-4
	nParms=3
	nData=3
	call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info > 4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	iflag=1
	call RegPureDev2(nData,nParms0,parm,error,iflag)
	qShape=1+1.90476*( c(1) -1 )
	rKadNm3=KcStar(1)
	dHkCal_mol=dH(1)*1.987*Tc(1)/1000
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kadNm3    dHkCal/mol nAs   nDs        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	if(LOUDER)write(*,*)'final devA,devT,devP=',(error(i),i=1,3)
	pause 'check errors'
	return
	end
	!******************************************************************************************************************************************
	SUBROUTINE RegPureDev2(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	dimension PARM(nParms),ERROR(nData),x(NMX) !,fugc(NMX) !,y(NMX),ierBp(12)
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	q(1)=1+1.90476*(c(1)-1)
	eokp(1)=Parm(2)
	VX(1)=parm(3)
	bVolCc_mol(1)=VX(1)
	isZiter=1 !use numerical derivatives
	toll=1.d-5
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUD)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	call CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,acen_pure,iErrCode)
	error(1)=acen_pure-acen(1)
	error(2)=TC_Pure-Tc(1)
	error(3)=PC_Pure-Pc(1)
	return
	end
!******************************************************************************************************************************************

	subroutine RegPureEsd2c(dHkJ_molP,rKadStarP,nParms0,iterOpt,iDenOptP,idToReg,nComps)!note:tc,pc,acen,id & esd parms passed in common
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (3/01) to put option in PrEsd package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDev2
	dimension PARM(NPMAX),ERROR(55),stdErr(NMX)
	common/ PUREVP /VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ CONSTK /KIJ(NMX,NMX),INITIAL
	common/ AVERR  /AAPERR,RHOERR,objfun
	common/ CONST  /rKadStar, dHkJ_mol, iDenOpt
	rKadStar=rKadStarP
	dHkJ_mol=dHkJ_molP
	iDenOpt = iDenOptP
	!if(rGas.lt.1.987)rGas=8.314
	do iComp=1,nComps
		if(id(iComp).eq.idToReg)iCompReg=iComp
	enddo

	nC=1 ! within RegPure, we must limit to only 1 compo.
	if(iCompReg.ne.1)then !the comp for regression must be in 1 position so calls to fugi etc sum only over that component
		tcStore=tc(1)
		pcStore=pc(1)
		acenStore=acen(1)
		eokPStore=eokP(1)
		kcStarStore=kcStar(1)
		dHStore=dH(1)
		cShapeStore=c(1)
		QStore=Q(1)
		vxStore=VX(1)
		ndStore=nd(1)
		tc(1)=tc(iCompReg)
		pc(1)=pc(iCompReg)
		acen(1)=acen(iCompReg)
	endif

	c(1)=1+acen(1)*2
	c(1)=1+rmw(1)/150
	u0ok=359+.66/rmw(1)
	v00=17*rmw(1)/50
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/8.314/Tc(1)
	!  general ***********************************
	PARM(1)=c(1)
	PARM(2)=v00
	PARM(3)=359+.66/rmw(1)
	PARM(4)=kcStar(1)
	PARM(5)=dH(1)
	nData=16 !this gives lots of points between tMax and tMin
	nParms=1 !this will force crit match estimate of bVol,eokP
	call RegPureDev2(nData,nParms,parm,error,iflag)
	if(LOUD)write(*,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	if(LOUD)write(*,*)'initial aaperr,sgerr',aaperr,rhoerr
	!outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	!open(66,file=outFile)
	!66 should be open from RegPureIo.
	write(66,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	write(66,*)'initial aaperr,sgerr',aaperr,rhoerr
	objo=objfun
	aaperro=aaperr
	if(iteropt.eq.0)then
		do while(rguess.gt.0)
			write(*,*)'enter C,vx,eok (-1 to skip to regression using last guess)'
			read(*,*)rguess,vguess,u0ok
			if(rguess.gt.0)then
				parm(1)=rguess
				parm(2)=vguess
				parm(3)=u0ok
				r=rguess
				v00=vguess
				call RegPureDev2(nData,nParms0,parm,error,iflag)
				write(*,*)'aaperr,sgerr',aaperr,rhoerr
			endif
		enddo
	else
		!	if(nParms.ge.2)then	!search for starting v00
		!		v00=v00/1.1
		!		parm(2)=v00
		!		call DevPure(nData,nParms,parm,error,iflag)
		!	    write(*,*)'aaperr,sgerr',aaperr,rhoerr
		!	    if(aaperr.lt.aaperro)then
		!	      aaperro=aaperr
		!		  goto 3000
		!	    endif
		!		parm(2)=parm(2)*1.1
		!	end if
		nParms=1 !get starting value by optimizing cShape only.
		NPO=ND(1)
		FACTOR=1.d-4
		TOL=1.D-4
		call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
		parm(2)=VX(1)
		if(nParms0.ge.3)then
			nParms=2
			FACTOR=1.d-4
			TOL=1.D-4
			aaperrB4=aaperr
			if(LOUD)write(*,*)'trying nParms=2'
			call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
			if(aaperr.gt.aaperrB4)then
				if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
				if(LOUD)pause 'Warning at nParms=2: aaperr.gt.aaperrB4'
			endif
		endif
	endif
	NPO=ND(1)
	nParms=nParms0
	qShape=1+1.90476*(C(1)-1)
	dHStar=dHkJ_mol*1000/8.314/tc(1)
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	aaperrB4=aaperr												   
	if(LOUD)write(*,*)'one last call to LmDif, aaperr ',aaperr
	!	pause
	factor=1.d-4
	tol=1.d-4
	parm(3)=eokP(1)
	!	nParms=3
	call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info.gt.4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	if(aaperr.gt.aaperrB4)then
		if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
		if(LOUD)pause 'Warning: aaperr.gt.aaperrB4'
	endif

	!evaluate actual vp error
	iFlag=86
	call RegPureDev2(nData,nParms,parm,error,iflag)
	qShape=1+1.90476*(C(1)-1)
	dHout=dh(1)*Tc(1)*8.314
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	WRITE(66,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	WRITE(66,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
601	FORMAT(I5,1X,2(F9.3,1X),2(F9.3,1X),F9.5,2I4,2(F7.2,1X),f7.4,1x,f7.2,1x,f7.4)
	if(LOUD)PAUSE 'check final results'
86	continue
	if(iCompReg.ne.1)then
		tc(iCompReg)=tc(1)
		pc(iCompReg)=pc(1)
		acen(iCompReg)=acen(1)
		eokP(iCompReg)=eokP(1)
		kcStar(iCompReg)=kcStar(1)
		dH(iCompReg)=dH(1)
		C(iCompReg)=C(1)
		Q(iCompReg)=Q(1)
		VX(iCompReg)=VX(1)
		nd(iCompReg)=nd(1)
		tc(1)=tcStore
		pc(1)=pcStore
		acen(1)=acenStore
		eokP(1)=eokPStore
		kcStar(1)=kcStarStore
		dH(1)=dHStore
		C(1)=cShapeStore
		Q(1)=QStore
		VX(1)=vxStore
		nd(1)=ndStore
	endif
	return
	end

	SUBROUTINE RegPureDev2c(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	INTEGER LIQ
	dimension ier(12)
	dimension PARM(*),ERROR(nData),fugc(NMX),x(NMX),y(NMX),ierBp(12)
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP !we need accurate VP for regression
	common/ CONST / rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	IF(TMAXVP/TC(1).LT.0.8)WRITE(66,*)'Warning TMAX<0.8Tc is TOO LOW'
	IF(TMINVP/TC(1).GT.0.6)WRITE(66,*)'Warning TMIN>0.6Tc is TOO HIGH'
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/Tc(1)/8.314

	!TR7=0.7*TC(1)
	!PR7=EXP( VPA+VPC*LOG(TR7)+VPB/TR7+VPD*TR7**VPE)/PC(1)
	!ACEN(1)=-LOG10(PR7)-1
	!???We must know ACEN from input
	IF(nParms.GE.2)VX(1)    =PARM(2)
	IF(nParms.GE.3)eokP(1)  =PARM(3)
	IF(nParms.GE.4)kcStar(1)=PARM(4)
	IF(nParms.GE.5)dH(1)    =PARM(5)
	Q(1)=1+1.90476*(C(1)-1)
	if(dh(1).gt.111)then
		if(LOUD)write(*,*)'dH is too large.  Spec kJ not J.'
		if(LOUD)pause
		return
	endif
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	if(nParms.ge.3)goto 2000
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev2 during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms.eq.1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms.le.2)eokP(1)=LOG(1.0617+Yc)*Tc(1)
2000 continue		
    !if(LOUD)write(*,'(a,5f10.3)')' c,Vx,eps/k',C(1),VX(1),eokP(1)
    !pause 'check values'
    if(nData==1)then !use the acentric factor to estimate the vp error
        pSat7 = Pc(1)*10**( 7*(1+acen(1))/3*(1-1/0.7) )
        call Fugi(TK,P,X,NC,1,FUGC,ZL,IER)
        fugcL=fugc(1)
        call Fugi(TK,pSat7,X,NC,0,FUGC,ZV,IER)
        ERROR(1)=( EXP( FUGCL-fugc(1) )-1 )*100
        return
    endif 
	tmin=TminVP*1.1
	if(tmin.lt. 0.45*Tc(1))tmin=0.45*Tc(1)
	tMax=tMaxVp*0.9
	if(tMax.gt. 700)tMax=700
	DELT=(TMAX-TMIN)/(nData-1)
	PERR=0
	IPTS=0
	iTempMax=nData-1
	DO IT=1,iTempMax,1
!		TK=TMIN+DELT*(IT-1)
		TK=TMIN+DELT*(IT)
		TR=TK/TC(1)
		P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE )
		if(iFlag.eq.86)then
			pCalc=P
			init=1
			itMax=111
			calL BUBPL(tK,x,nC,init,pCalc,itMax,y,ierBp)
			perri=(pCalc-P)/P*100
			if(LOUD)write(*,'(f7.2,4f10.5)')tK,pCalc,P,perri
			IPTS=IPTS+1
		else
			LIQ=1
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
			fugcL=fugc(1)
			LIQ=0
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZV,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZV,IER)
			fugcV=fugc(1)
			IPTS=IPTS+1
			perri=( EXP( FUGCL-FUGCV )-1 )*100
		endif
		ERROR(IPTS)=perri
		if(kcStar(1).lt.0)error(ipts)=error(ipts)*10000
		if(C(1).lt.1)error(ipts)=error(ipts)*10000
		if(abs(zl-zv).lt.1.d-2)error(ipts)=error(ipts-1)*5
		if(ier(1).ne.0)error(ipts)=error(ipts-1)*10
!		if(abs(u0ok-u0okest)/u0ok.gt.0.03)error(ipts)=error(ipts)*1000
		PERR=PERR+ABS(perri)
		objfun=objfun+error(ipts)*error(ipts)
	enddo

	AAPERR=PERR/(nData-1)

	TK=298
	!P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE ) !
	P=PC(1)*10**( 7/3*(1+ACEN(1))*(1-TC(1)/TK) ) !We don't need super accurate P to estimate rhoLiq.
	if(P.lt.0.1)P=0.1
	LIQ=1
	call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
	!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
	RHOL=RMW(1)*P/8.314/TK/ZL
	IPTS=IPTS+1
	RHOERR=(SG-RHOL)/SG*100
	if(iDenOpt.eq.1)ERROR(nData)=iDenOpt*RHOERR
	if(LOUD)write(*,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	WRITE(66,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	!pause 'check deviations'
	!ERROR(nData)=0
	RETURN
	END

!***********************************************************************
	SUBROUTINE RegPureDev3b(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!INTEGER LIQ
	dimension ier(12)
	dimension PARM(*),ERROR(nData) ,fugc(NMX),x(NMX) !,y(NMX),ierBp(12)
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	Q(1)=1+1.90476*(C(1)-1)
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev3 during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms==1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms < 3)eokP(1)=LOG(1.0617+Yc)*Tc(1)
	parm(2)=eokP(1)
	parm(3)=VX(1)
	pSat7=Pc(1)*10**( -1-acen(1) )
	T=Tc(1)*0.7D0
	call FugiESD2(T,pSat7,X,nC,1,FUGC,zFactor,ier)
	chemPoLiq=FUGC(1)
	call FugiESD2(T,pSat7,X,nC,0,FUGC,zFactor,ier)
	chemPoVap=FUGC(1)
	error(1)=chemPoLiq-chemPoVap
	if(iFlag==1.and.LOUD)print*,'RegPureDevEsd3: c,chemPoDev',C(1),error(1) !mostly to avoid iFlag warning.
	RETURN
	END

!***********************************************************************
      SUBROUTINE RegPureIo2(NC)
	!C
	!C  PURPOSE:  Regress Pure component parameters of ESD model.
	!C  PROGRAMMED BY:  JRE 3/01
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER fileNameVp*55
	character outFile*251
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ CONST /  rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	print*,'Minimize vp error based on estd rKadStar and dHkJ_mol '  
	print*,'constrained to match Tc. Pc is ~matched implicitly by vp' 
	print*,'or matched explicitly when nParms=1.'
	print*,'best setting is nParms=2 and iteropt=0 => flex to vx.'   
	print*,'kcStar comes from kcStar~kadStar/cShape so kadStar=.025'
	print*,'Data from c:\thermo\esd\VpCcEsdCorr.txt if available.'
	print*,'Enter dHkcal/mol,rKadStar(eg 0.025),iDenOpt(0=ignore SG)'
	read(*,*)dHkJ_mol,rKadStar,iDenOpt
	dHkJ_mol=dHkJ_mol*4.184
	!rKadStar=0.025
	!iDenOpt=0		!0=ignore density in fit, 1=weight sg as one datum
	!initial=0

	WRITE(*,*)'ENTER nParms,iteropt(0=AutoReg, 1=guess then reg)'
	READ(*,*)nParms0,iteropt
	idToReg=id(1)
	ndTemp=1
	iCompToReg=1
	if(nc.ne.1)then
		write(*,*)'Enter dippr id # and Nd of component for regression '
		read(*,*)idToReg,ndTemp
		do iComp=1,nc
			if(id(iComp).eq.idToReg)iCompToReg=iComp
		enddo
		nd(iCompToReg)=ndTemp
	endif
	fileNameVp='c:\thermo\esd\VpCcEsdCorr.txt'
	!write(*,*)'Enter filename for vpdata'
	!read(*,'(a)')fileNameVp
	iFound=0
	open(55,file=fileNameVp,ERR=86)
	read(55,*,END=86,ERR=86)nDbaseVp !error here means no file so skip cc search.
	write(*,*)'Enter chemcad id # of component for regression '
	read(*,*)idToRegCC
	do iVp=1,nDbaseVp
		READ(55,*,END=86)idDb,rmwDb,sgDb,vpADb,vpBDb,vpCDb,vpDDb,vpEDb,tMinVpDb,tMaxVpDb
		if(idDb.eq.idToRegCC)then
			iFound=1
			rmwDum=rmwDb
			sg=sgDb
			vpA=vpADb
			vpB=vpBDb
			vpC=vpCDb
			vpD=vpDDb
			vpE=vpEDb
			tMinVp=tMinVpDb
			tMaxVp=tMaxVpDb
		endif
	enddo
	CLOSE(55)
86	continue !if the dbase doesn't exist, then we must type the data in.
	if(iFound.eq.0)then
		write(*,*)'Compd not found in std dbase.'
		write(*,*)'Enter Mw,SpeGrav,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp'
		read(*,*)rMwDum,sg,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp
	endif
	
	vpA=vpA-LOG(1.E6)!convert dippr vp from Pa to MPa.

	TRLO=tMinVp/TC(iCompToReg)
	TRHI=tMaxVp/TC(iCompToReg)
	outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	open(66,file=outFile)
	IF(TRLO.GT.0.6.OR.TRHI.LT.0.8)THEN
	  if(LOUD)write(* ,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	  WRITE(66,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	END IF

	!  we need a big function here because of switching comp pos., and 
	!initial guesses to converge.
	call RegPureEsd2c(dHkJ_mol,rKadStar,nParms0,iterOpt,iDenOpt,id(iCompToReg),nC)!note:tc,pc,acen,id & esd parms passed in common
	close(66)
	if(LOUD)pause 'Results in RegPureOut.txt'
	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE KsVAL2(nc,nParms,parm,tKelvin,deviate,kAADk,kBias,kErrMax,LAST)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DOUBLE PRECISION parm(nParms),deviate(maxPts)
	parameter (nEta=5)
	common/SimData/ xFracc(23),vEffNm3(23),a0xss(20,20),A0Mix(20,20),A1Mix(20,20),A2Mix(20,20),Nps
	!common/XsProps/a0XsSs
	!dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	!common/ksvall/ks0ij,ks1ij
	DIMENSION XXFRAC(NMX)
            
	607	FORMAT(1X,F8.4,1X,F10.4,1X,F10.8,1X,F10.8,1X,F10.4,1X,i5)      

	!KS0IJ(1,2)=0
	!KS1IJ(1,2)=0
	!KS0IJ(2,1)=KS0IJ(1,2)
	!KS1IJ(2,1)=KS1IJ(1,2)
	!KS0IJ(1,1)=0
	!KS1IJ(1,1)=0

		if(last.lt.0)then !only kijOpt sets last to -1, so that is the only one that effects this.
			write(*,*)'Enter optimization choice: 1=ks0, 2=ks0,ks1'
			read(*,*)iOpt
		endif
		if (nParms==1)then
			KS0IJ(1,2)=parm(1)
		else
			KS0IJ(1,2)=parm(1)
			KS1IJ(1,2)=parm(2)
		endif
	MAXIT=1111
	KAADK=0
	ssqErr=0
	kBias=0
	eta=0.d0
	!tKelvin=273
	ikk=1
	nTot=nEta*(NPS-2)
	!Do ikk=1,nTot
	nFracs=NPS-2
		Do iEta=1,nEta
			eta=eta+0.1d0

			DO iData=1,nFracs
			!if(iData==(NPS-1))then
			!	xxFrac(1)=1
			!	xxFrac(2)=0
			!elseif(iData==NPS)then
			!	xxFrac(1)=0
			!	xxFrac(2)=1
			!else
				xxFrac(1)=xFracc(iData)
				xxFrac(2)=1-xxFrac(1)
			!endif
				IFLAG=0
				ITMAX=MAXIT
				INIT=1
				nComps=nc  
				CALL AxsCalc(xxFrac,eta,nComps,tKelvin,a0MixCalc,iErr)
				!a0XsCalc=a0MixCalc-(xxFrac(1)*A0Mix((NPS-1),iEta)+xxFrac(2)*A0Mix(NPS,iEta)) 
			!a0Xs=a0Mix-(xxFrac(1)*

				deviate(ikk) =(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0XsCalc-a0xss(iData,iEta))/a0xss(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))**2! !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100	 !
				if(ABS(deviate(ikk)).gt.ABS(kErrMax))kErrMax=deviate(ikk)
				
				!penalize the deviate passed to lmdif, but only compile paadp for converged pts
				KAADK = KAADK+ABS(deviate(ikk))
				ssqErr= ssqErr+deviate(ikk)*deviate(ikk)
				kBias=kBias+deviate(ikk)
		!	if(ABS(deviate(iData)).gt.ABS(pErrMax))pErrMax=deviate(iData)
			
				!WRITE(67,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)

				ikk=ikk+1
			enddo
		enddo
	!enddo

	rmsErr=SQRT(ssqErr/(nEta*(NPS-2)))
	KAADK=KAADK/(nEta*(NPS-2))
	kBias=kBias/(nEta*(NPS-2))
	if(LOUD)write(*,'(7f10.3)')(parm(i),i=1,nparms),KAADK,rmsErr

	RETURN
	END


