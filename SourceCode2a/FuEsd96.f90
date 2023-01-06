!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine GetEsd96Cas(nC,idCasPas,iErr) !ID is passed through GlobConst
	!
	!  PURPOSE:  LOOKS UP THE ESD PARAMETERS AND STORES THEM IN COMMON EsdParms
	!
	!  INPUT
	!    ID - VECTOR OF COMPONENT ID'S INPUT FOR COMPUTATIONS
	!  OUTPUT
	USE GlobConst !is implied by USE ASSOC. GlobConst includes ID,Tc,Pc,Zc,acen,...
	USE Assoc ! ESD has it's own Parms. !nTypes(=1 for ESD), nDegree, nDonors, nAcceptors, bondRate(?), includes GlobConst{Tc,Pc,...}
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=3000,listPool=1000)
	character*222 bipFile,inFile,dumString !,dumString,ParmsTptFile*50 !,bipHbFile*50
	integer iGotIt(NC),idCasPas(NC),idCasa(ndb),GetBIPs,ierCompExact(NC)  !,idStore(NC) 
	doublePrecision KCSTA(ndb),DHA(ndb),bVolA(ndb),eokA(ndb),cA(ndb),qA(ndb)
	integer IDA(ndb),NDSA(ndb),NASA(ndb),NDA(ndb)
	doublePrecision bondRate(nmx,maxTypes)
	!iErr=1 !Parms missing and iErrExact.ne.0 for at least one component
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER=.TRUE.
    !do i=1,NC
    !    idCas(i)=idCasPas(i)
    !enddo
	idCas(1:NC)=idCasPas(1:NC) ! workaround after promoting idCas to GlobConst 
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isTPT or isPcSaft.
	isESD=.TRUE. ! in GlobConst, simplifies calls in FuVtot or FUGI 
	etaMax=1/1.9D0-zeroTol
	IF(DEBUG)then  !DEBUG is part of GlobConst
		inFILE='c:\Spead\CalcEos\input\ParmsEsd96.TXT'
		if(iEosOpt==12)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmami.TXT'
		if(iEosOpt==13)inFILE='c:\Spead\CalcEos\input\ParmsEsdEmamiTb.TXT'
	ELSE
		inFile=TRIM(masterDir)//'\input\ParmsEsd96.txt' ! ESD96(MEM1) with optimized compound parameters for associating compounds
		if(iEosOpt==12)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmami.txt' ! // ESD96 with group contribution parameter estimates
		if(iEosOpt==13)inFILE=TRIM(masterDir)//'\input\ParmsEsdEmamiTb.txt' ! // ESD96 with group contribution parameter estimates except EOK is varied to match Tb760
	ENDIF
	OPEN(31,FILE=inFile)
	if(LOUDER)print*,'GetEsd:inFile=',TRIM(inFile)
	if(LOUDER)pause 'Check the ESD parms file location.'

	!READ(31,'(a88)',ERR=861)dumString 	!read title line. Ignore nDeck if present.
	READ(31,*)nDeck
	do i=1,nDeck
		READ(31,'(a222)')dumString
		READ(dumString,*,ioStat=ioErr)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I) ,idCasa(i)
		!write(*,*)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I)  ,idCasa(i)
		if(ioErr /=0 .and. LOUDER)print*,'GetESD: error reading ParmsEsd_.txt. line=',TRIM(dumString)
		if(  ( idCasa(i)==id(1) .or. idCasa(i)==id(2) ) .and. LOUDER  )print*,'Found in ParmsEsd idCas=',idCasa(i) 
	enddo
	NDI=1
	I= -1  !overrides headerless reading.
	DO while(I.ge.0) !reads to end of file w/o reading nDeck
		READ(31,*,ERR=861,END=50)IDA(I),CA(I),QA(I),eokA(I),bVolA(I)&
			,NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I),idCasa(i)
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
		ierCompExact=1 ! Set error for all comps. Declare error because iEosOpt=12,13 means using only GC parameters from one of ParmsEsdEmami__
	else
		if(Tc(1) < 4)call GetCritCas(nC,iErrCrit) !GetCritCas assumes ID(GlobConst)=IdCas
		call ExactEsd96(nC,VX,C,Q,eokP,iErrExact,ierCompExact) !iErrExact = 100+iComp if compd is assoc or Asso+. Wait to see if Parms are in ParmsEsd before failing.
		if(LOUDER.and.iErrExact/=0)print*,'GetESDWarning: iErrExact=',iErrExact,' Checking database for iComp='	,( ierCompExact(iComp),iComp=1,NC)
	endif
	nTypes(1:NC)=1	 !all esd versions use nTypes=1. Multifunctional molecules require SPEADMD. 

	iType=1 	 ! Only one siteType per molecule in ESD for now. Can generalize more later if desired.
	DO J=1,NC
        bVolCC_mol(j)=Vx(j) !Copy ExactEsd value first. Replaced below if in dbase.
		ND(J)=0
		NDS(J)=0
		NAS(J)=0
		iGotIt(J)=0
		DO I=1,NDECK
			IF(IDCASA(I).EQ.IDCas(J))THEN
				iGotIt(J)=1
				eokP(J)=eokA(I)
				C(J)=CA(I)
				ND(J)=NDA(I)
				Q(J)=QA(I)
				Vx(J)=bVolA(I)
				bVolCC_mol(J)=Vx(J) !need bVol generally for every EOS. Included in GlobConst
				KCSTAR(J)=KCSTA(I) ! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				DH(J)=DHA(I)*1000/1.987/ Tc( j)			! new format for ParmsEsd used bondVolNm3, and epsHbKcal/mol.  JRE 20200428
				NAS(J)=NASA(I)
				NDS(J)=NDSA(I)
				iComplex=0
				if(NAS(j)*NDS(j) > 0)iComplex=1
				if(iComplex==0 .and. iEosOpt==2)KCSTAR(j)=0 !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				if(iComplex==0 .and. iEosOpt==2 .and. LOUD)pause 'iComplex=0 for assoc???' !ParmsEsd may contain parameters for compounds that solvate but do not associate (for iEosOpt==4). These need to be killed for iEosOpt==2.
				if(LOUDER)print*,'iGotIt! idCas,bVol=',idCas(j),Vx(J)
				!exit !exits this do loop, not the outer one.
			ENDIF
		enddo
		nDegree(J,iType)=ND(J)
		nDonors(J,iType)=NDS(J)
		nAcceptors(J,iType)=NAS(J)
		bondRate(J,iType)=0.d0	   
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
		if(LOUDER)pause 'GetEsd: error for at least one compound'
		continue
	endif
    
    nTypesTot=0
    do iComp=1,NC !tabulate the aBip matrix indices.
        nTypes(iComp)=1             !For ESD, assume just one site type for each molecule
        idType(iComp,1)=iComp       !    and each molecule is a separate type.
        !nFg(iComp,1)=1             !Same reasoning.
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
        enddo !iType
    enddo !iComp
    
	if(LOUDER)then
        write(*,*)'  ID     NAME       cESD    eok      bVol    Nd     KAD   eHB/k(K)'
	    do i=1,NC
		    write(*,606)IDCas(i),NAME(i),C(i),eokP(i),Vx(i),ND(i),KCSTAR(i),DH(i)*1.987*Tc(i)/1000
        enddo
    end if
606	format(i9,1x,a11,f9.3,f8.2,f8.2,i3,1x,f8.6,f8.4)

	!note:  bips are passed back through common/BIPs/
	IF(DEBUG)then 
		bipFile='c:\SPEAD\CalcEos\input\BipEsd96.txt'
	ELSE 
		bipFile=TRIM(masterDir)//'\input\BipEsd96.txt' ! // is the concatenation operator
	ENDIF
	if(NC > 1)iErrCode=GetBIPs(bipFile,ID,NC) !not necessary for pure fluids
	if(iErrCode > 10)iErr=11 ! 
    if(LOUDER)then
		print*,'GetEsdCas: bipFile=',TRIM(bipFile)
		print*,'     ',(id(j),j=1,NC)
		do i=1,NC
			write(*,'(i5,11f8.4)')id(i),(Kij(i,j),j=1,NC)
		enddo
	    pause 'GetESD: check BIPs.'
		if(iErrCode > 0) print*,'GetESD96: BIPs missing for ',iErrCode-10,' binary combinations.'
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
subroutine ExactEsd96(nC,VX,c,q,eokP,iErr,ierComp)
	USE GlobConst
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD parameters for hydrocarbons based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
	!parameter (NMX=55)
	DoublePrecision k1
	DoublePrecision VX(nmx),c(nmx),Q(nmx),eokP(nmx)
	integer ierComp(NC) !,ID(nmx)
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
			if(LOUDER)print*,'ExactESD: Parms not available for assoc ID=',ID(i)
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
		cShape=1+3.535d0*Wci+0.533d0*Wci*Wci
		RootCinv=1/SQRT(cShape)
		rootC=SQRT(cShape)
		k1=1.7745d0
		qShape=1+1.90476d0*(cShape-1)
		ZcTmp=(1+0.115d0/rootC-0.186d0/cShape+0.217d0/cShape/rootC-0.173d0/cShape**2)/3	!EL2ed Eq. 19.125
		!Alt: ZcTmp=(1+0.1388/SQRT(q)-0.171/q-0.01874/q/SQRT(q)+0.02361/q**2)/3
		atemp=9.5d0*qShape*1.9d0+4*cShape*k1-k1*1.9d0
		quadB=k1*1.9d0*ZcTmp+3*atemp
		sqArg=quadB*quadB+4*atemp*(4*cShape-1.9d0)*(9.5d0*qShape-k1)/ZcTmp

		Bc=ZcTmp*ZcTmp*(-quadB+SQRT(sqArg))/(2*atemp*(4*cShape-1.9d0))
		Yc=ZcTmp*ZcTmp*ZcTmp/(atemp*Bc*Bc)
		rlnY1=LOG(Yc+1.0617d0)
		c(i) = cShape
		Q(i) = qShape
		eokP(i)=TC(i)*rlnY1
		bVolCc_mol(i)=Rgas*TC(i)/PC(i)*Bc
		VX(i)=bVolCc_mol(i)
	enddo
	return
end	!subroutine ExactEsd
!------------------------------------------------------------------------------ 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------ 
	!	FugiESD
	!   LATEST REVISION : 
	!	9/94 jre implemented ESD96 method
	!	1/96 (switched to chempot, ADDED POLYETHYLENE  jre)
	!	7/96 PS, PPO, PEO, PIB (ram natarajan)
	!	1/97 PS, PPO, PEO, PIB (made consistent, jre)
	!	7/06 jre ->f90, sample calcs, check <k1Yb> rule from 91, cf[2].
	!   11/22 jre: implemented module structure, Trange error checking, and calling FuVtot instead of redundancy in FugiESD
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

	SUBROUTINE FugiESD96(tKelvin,pMPa,xFrac,NC,LIQ,FUGC,zFactor,ier)
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	  + GlobConst{rGas,Tc,Pc,...}
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision xFrac(NMX),FUGC(NMX) !,chemPoAssoc(nmx)
	Integer ier(12) !,iera(12)
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
	ier=0

	if(LOUDER)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUDER)print*,'FugiEsd: Kij(1,2)= ',checkKij12
	
	! NOTE: iErrTmin is checked in FuVtot

	sumx=SUM( xFrac(1:NC) )
	if(ABS(sumx-1) > 1e-8)then
		if(LOUDER)pause 'FugiEsd: sumx .ne. 1'
	endif
	!INITIATE SECANT ITERATION ON rho
	bMix=SUM( xFrac(1:NC)*bVolCC_mol(1:NC) )
	Pb_RT=pMPa*bMix/(Rgas*tKelvin)
	!GUESS FOR rho
	eta=Pb_RT/1.05D0  !NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
	IF(LIQ==1 .or. LIQ==3 .or. eta>etaMax)eta=etaMax-0.05d0
	rho=eta/bMix
	if(eta > 1/1.9 .and. LOUDER)print*,'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
	isZiter=1 ! FUGC calculations are skipped for isZiter=1
	Call FuEsd96Vtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	IF(iErr > 0)GOTO 86
	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor

	eta=etaOld/1.05D0
	IF (eta < 0 .and. LOUD) WRITE(6,31)LIQ
	rho=eta/bMix
	if(initial.and.LOUD)print*,'FugiEsd: initial eta,err',etaOld,errOld
	itMax=77
	errBestEta=1234
	do nIter=1,itMax
		Call FuEsd96Vtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		IF(iErr > 0)EXIT
		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial.and.LOUDER)write(*,'(a,2e11.4,3f10.5)')'FuEsd96 eta,Z', eta,zFactor
		if(initial.and.LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')'FuEsd96 eta,CHNG,niter',eta,CHNG,niter 
		etaOld=eta
		ERROLD=ERR
		!  LIMIT THE CHANGE IN Density for liquid..
		!  Low eta must move from zero, so you can't limit its % change
		IF(liq==1.and.DABS(CHNG/etaOld).GT.0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==0.and.DABS(CHNG).GT.0.02d0)CHNG=DSIGN(0.02D0,CHNG)
		IF(liq==3.and.DABS(CHNG/etaOld).GT.0.1D0)CHNG=DSIGN(0.1D0,CHNG)*etaOld
		IF(liq==2.and.DABS(CHNG).GT.0.02d0)CHNG=DSIGN(0.02D0,CHNG)
 		eta=eta-CHNG
		if(ABS(CHNG) < errBestEta)then
			etaBest=eta
			errBestEta=ABS(CHNG)
		endif
		if(eta < 0 .or. eta > 1/1.9)eta=etaOld-DSIGN(0.1D0,CHNG)*etaOld
		IF(DABS(CHNG) < 1.D-9 .and. eta > 0)EXIT  ! Don't use CHNG/eta here. Converge quickly to ideal gas if P->0, still ~9 sigfigs if liquid.
	enddo !nIter=1,itMax
	if(eta < 0 .or. eta > 1.9)then
		if(LOUDER)pause 'eta < 0 or > 1.9 on final iteration in FugiESD.'
		continue
	endif

	if(iErr.ne.0.or.nIter.ge.itMax.or.eta < 0)then
		ier(4)=iErr
		eta=etaBest
		Call FuEsd96Vtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		ier(11)=1
		ier(1)=10
	endif
	!One last call to get FUGC.
	if(initial.and.LOUD)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd96 cnvrgd: eta,CHNG,niter',eta,CHNG,niter
	etaPass=eta
	rho=eta/bMix 
	IF (rho < 0)THEN
		ier(5)=1
        ier(1)=15
		if(LOUDER)WRITE(6,31)LIQ
		GOTO 86
	ENDIF
!  ITERATION ON rho HAS CONCLUDED.  DEPARTURES PASSED THROUGH GlobConst.
	if(ABS(eta-rho*bMix) > 1E-11 .and. LOUDER)pause 'eta.ne.(rho*bMix)?'
	if(pMPa==0 .and. LOUDER)print*,' FugiEsd: P=0? LIQ,P=',LIQ,pMPa
	!zFactor=P/(rho*rGas*T)  ! add this to improve precision when computing rho from Z after return.   
	if(zFactor.le.0 .and. LOUDER)print*,'FugiEsd: converged Z <= 0. eta,Z=',eta,zFactor
	isZiter=0
	rho=eta/bMix
	Call FuEsd96Vtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	FUGC(1:NC)=FUGC(1:NC)-DLOG(zFactor)	 !Must subtract ln(Z) when given Vtot as independent variable.
	if(LOUDER)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd96: eta,CHNG,niter,FUGC',eta,CHNG,niter,(FUGC(i),i=1,NC) 
	initial=0
	return
86	if(LOUDER)WRITE(6,*)' ERROR IN FuEsd96.  '
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! rho -VE IN FUGI')
	IF(NITER.GE.ITMAX)THEN
		if(LOUDER)write(*,*)'TOO MANY Z ITERATIONS'
		ier(6)=1
        ier(1)=16
	END IF
	IF(ier(4)==1.and.LOUDER) WRITE(*,*)'ERROR IN FuEsd96Vtot'
	if(ier(1) < 10)ier(1)=11
	initial=0
	RETURN
	END	!FugiESD()
	!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!  ELLIOTT'S SUBROUTINE 
	!  REVISED  9/94 FOR POLYSEGMENTED MOLECULES LIKE GLUTARAL
	!  REVISED 10/94 FOR ANY # OF ASSOCIATING SPECIES
	SUBROUTINE AlpSolEz(X,NC,KVE,ND,VM,NAS,NDS,rho,zAssoc,XAPas,RALPH,fAssoc,ier)
	USE GlobConst
	USE ASSOC !GlobConst+XA+XD+XC
	!  PURPOSE:  COMPUTE THE EXTENT OF ASSOCIATION (fAssoc) AND zAssoc given rho,VM
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!PARAMETER(NMX=55)
	DOUBLEPRECISION X(NMX),KVE(NMX),RALPH(NMX),rLnPhiAssoc(NMX),XAPas(NMX,maxTypes)
    Integer ND(NMX),NAS(NMX),NDS(NMX),ier(12)
	LOGICAL LOUDER 
	!  NDi     = THE DEGREE OF POLYMERIZATION OF COMPO i
	!  fAssoc  = THE CHARACTERISTIC ASSOCIATION = 1/RALPHi(1/XAi-1)
	!  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei/(1-1.9eta)
	!  ier     = 4  eta > 0.53
	LOUDER=LOUD	! from GlobConst
	!LOUDER = .TRUE. ! for local debugging.

	eta=rho*VM
      
	DO I=1,6
		ier(I)=0
	ENDDO                                  

	DO I=1,nC           
		SQARG=KVE(I)/VM*eta/(1-1.9D0*eta)
		IF(SQARG < 0)THEN
			if(LOUD)write(*,*)'SQARG < 0 IN ALPSOL '
			ier(4)=1
			RETURN
		ENDIF
		RALPH(I)=DSQRT( SQARG )
	ENDDO

	!BEGIN SECANT ITERATION TO DETERMINE fAssoc
	fAssoc=0            
	FA0=0    ! cf. Elliott IECR 61:15724 Eq.18
	!FA1=0
	ralphMean=0                    
	DO I=1,nC
		iComplex=nAs(I)*nDs(I)	    !Esd96 means that all must bond or none
		if(iComplex.ne.0)iComplex=1 ! and ralph takes care of solvation
		FA0=FA0+X(I)*ND(I)*RALPH(I)*iComplex
		!FA1=FA1+X(I)*ND(I)*RALPH(I)*iComplex/(1+ralph(i)) ! RHS of Eq.18 with FD=1
		if(ralph(i) > ralphMean)ralphMean=ralph(i) ! Use infinity norm to estimate ralphMean initially.
	enddo
	ERROLD=fAssoc-FA0
	IF(ABS(ERROLD).LE.1.D-8)GOTO 65 !no need to iterate if sum=0 (no hbonding)
	fOld=fAssoc
	! FA0/(1+ralphMean*1) ~ FA1 => ralphMean ~ FA0/FA1-1, rearranged form of Eq.18
	!if(FA1 > 0)ralphMean=FA0/FA1-1
	fAssoc=2*FA0/( 1+DSQRT(1+4*ralphMean*FA0) )  ! Eq.23 where FD0=FA0.	This is exact for binary like benzene+methanol.
	NITER=0
55	NITER=NITER+1                  
		IF(NITER.GT.111)THEN
			if(LOUD)write(*,*)'ERROR - NO CNVRG ON fAssoc'
			ier(4)=1
			GOTO 86
		ENDIF
		SUM=0
		DO I=1,nC      
			iComplex=nAs(I)*nDs(I)	    !SymEsd means that all A&D must bond or none
			if(iComplex.ne.0)iComplex=1 ! and ralph takes care of solvation
			SUM=SUM+X(I)*ND(I)*RALPH(I)*iComplex/(1+fAssoc*RALPH(I))
		enddo !For pure F=ralphNd/(1+F*ralph)=> F+ralph*F^2-ralphNd=0 =>F=[-1+sqrt(1+4*Nd*ralph^2)]/(2*ralph)=(1-1+4*alphNd)/[2*ralph*(1+sqrt(1+4*ralph^2))]=2*Nd*ralph/[1+sqrt(1+4*ralph^2)]
		ERR=fAssoc-SUM
		CHANGE=ERR/(ERR-ERROLD)*(fAssoc-FOLD)
		FOLD=fAssoc
		ERROLD=ERR
		fAssoc=fAssoc-CHANGE
	IF(ABS(ERR/fAssoc).GT.1.D-8)GOTO 55
65	CONTINUE

	!  fAssoc ITERATION HAS CONCLUDED

	zAssoc= -fAssoc*fAssoc/(1-1.9D0*eta)
	DO I=1,nC
		XAPas(I,1)=1/(fAssoc*RALPH(I)+1)
		rLnPhiAssoc(i)=2*ND(I)*LOG( XAPas(I,1) )+zAssoc*1.9d0*bVolCC_mol(i)*rho
	enddo

	if(LOUDER)print*,'eta,Fassoc,Zassoc,ralph(),lnPhiAssoc'
	if(LOUDER)write(*,'(8f10.4)')eta,fAssoc,zAssoc,(ralph(i),i=1,NC),(rLnPhiAssoc(i),i=1,NC)
86	CONTINUE
	RETURN
	END	 !AlpSolEz()

!**************************************************************************
	SUBROUTINE GoldenZesd(X,nC,KVE,ND,bMix,NAS,NDS,CbMix,YqbMix,k1YbMix,&
                       PoRT,LIQ,zFactor,rho,zAssoc,XA,RALPH,fAssoc,ier)
	USE GlobConst, only:NMX
	USE Assoc, only:maxTypes
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DoublePrecision X(NMX),XA(NMX,maxTypes),KVE(NMX),RALPH(NMX)
	Integer ND(NMX),NAS(NMX),NDS(NMX),ier(*)
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
    
	          
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	Written by AFG, Oct. 2009																				C
!C	With Having T,V and nMols, this routine calculates the Z factor and all derivatives needed in 			C
!C	critical point,	flash and bubble point calculations.													C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuEsd96Vtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,Ares,Ures,iErr)
	USE GlobConst
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(NMX),xFrac(NMX),FUGC(NMX),KCSTARp(NMX) !,XA(NMX)
	integer IER(12),iera(12) !,kComp
	DoublePrecision YQVIJ(NMX,NMX),KVE(NMX),YQVI(NMX),Y(NMX,NMX),EOK(NMX,NMX)
	DoublePrecision CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX)
	DoublePrecision k1(NMX) !,RALPH(NMX)
	DoublePrecision dEOK_dT(NMX,NMX),dY_dT(NMX,NMX),dYQVI_dT(NMX),dYQVIJ_dT(NMX,NMX)
	DoublePrecision dCVM_dN(NMX),dYQVM_dN(NMX) ,dVM_dN(NMX),dK1YVM_dN(NMX) !,dCVI_dN(NMX,NMX)
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
	  
	iErr=0
	!zeroTol=1D-11  !zeroTol is now global
	totMoles=sum(gmol)
	xFrac(1:NC)=gmol(1:NC)/totMoles
	if( tKelvin	   < zeroTol .or. totMoles < zeroTol .or. vTotCc < zeroTol)then
		if(LOUD)print*,'FuEsdVtot: nonsense T(K),totMoles,vTotCc=',tKelvin,totMoles,vTotCc
		iErr=1
	endif
	bMix=0
	iErrTmin=0
	TminTot=.01
	TrVolatile=0.1
	if(LOUDER)call QueryParMix(1,checkKij12,iErrBip)
	if(LOUDER)print*,'FugiEsd: Kij(1,2)= ',checkKij12
	DO I=1,NC
		rLogPr=DLOG10( 0.0001/Pc(i) )	! ESD not recommended below 0.0001 MPa for pure fluids. For mixes, ensure most volatile comp has Psat > 0.0001 MPa. Think about polymers.
		aScvp=7*(1+acen(i))/3	 !SCVP: log10(Pr)=7(1+w)/3*(1-1/Tr) = a*(1-x) 
		xt1= 1-rLogPr/aScvp  ! x = 1/Tr at Psat=0.0001 MPa, first approximation of x = x1 + dx
		xt = xt1 -0.178*acen(i)*acen(i)*(1-xt1)*(1-xt1)  ! This crude empirical correlation was developed for nonadecanol.  cf. PGL6Samples.xlsx(nC19oh).
		if( xt > 2.222)xt = 2.2222	!1/2.2222 = 0.45. If 
		if( xt < 1 .and. LOUD)pause 'FugiEsd: TrMin > Tc???'
		TrMin = 1/xt ! = min( Tr@Psat=0.0001 or 0.45 )
		if(LOUDER)print*,'xt1,xt,TrMin',xt1,xt,TrMin
		if( tKelvin/ Tc(i) < TrMin .and. NC==1)iErrTmin=iErrTmin+1
		if( tKelvin/Tc(i) > TrVolatile) TrVolatile=tKelvin/Tc(i)  ! The largest TrVolatile is the Tr of the compd with lowest Tc. 
		if( Tc(i)*TrMin > TminTot) TminTot=Tc(i)*TrMin	 ! The largest Tmin is the weakest link. 
		if( Tc(i)*TrMin > TminTot .and. LOUD) print*,'i,Tmin(i): ', i,Tc(i)*TrMin
		IF(xFrac(I) < 0 .and. LOUDER)PAUSE 'FugiEsd: ERROR - Xi<0'
	enddo 
	if(TrVolatile < 0.44d0)iErrTmin =2 ! it's only a problem if the most volatile compound has Tr < 0.45 or Psat < 0.0001.
	if(iErrTmin > 0) then
		ier(1)=5
		ier(2)=iErrTmin
		if(LOUDER)print*,'FugiEsd: T(K) < Tmin(all i)',tKelvin,TminTot
		!if(LOUD) pause 'FugiEsd: at least one compound has Tr < TrMin'
		!return  !! make this a warning for Vex,Hex etc.
	endif  
	rho=totMoles/ vTotCc
	eta=rho*bMix 
	if(LOUDER)print*,'FuEsdVtot: bMix,eta=',bMix,eta
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
		!eHbKcal_mol(I,iType)=dHkcalMol(I)
		!bondVolNm3(I,iType)=bondVolNm3Esd(I)
	enddo
	if( ABS(VX(1) - bVolCC_mol(1)) > zeroTol .and. LOUD ) print*,'FuEsdVtot: VX.ne.bVol=',VX(1),bVolCC_mol(1)
	!endif	!this was not working when computing vp
	if(LOUD.and.k1yvm < zeroTol)print*,'FuEsdVtot: 0~k1yvm=',k1yvm 
	eta=rho*vm
	!call WERTHEIM(vMolecNm3,eta,tKelvin,xFrac,NC,zAssoc,aAssoc,uAssoc,iErrCode)
	call AlpSolEz(xFrac,NC,KVE,ND,VM,NAS,NDS,rho,zAssoc,XA,RALPH,fAssoc,iera) ! passes alphAD,alphDA through common.
	FassocAfg=Fassoc ! I need Fassoc for the GX option. JRE 20210723
	IF(iera(1).ne.0)then
        ier(1)=14
        GOTO 86
    endif
	voidFrac=1-1.9D0*ETA
	denom=voidFrac
	ZREP= 4.d0*CVM*RHO/denom  ! old form.
	ZREP= 4.d0*Cmix*eta/denom
	ZATT= -ZM*YQVM*RHO/(1.d0+K1YVM*RHO)
	zFactor=(1.d0+ZREP+ZATT+ZASSOC)
	PMpa=zFactor*rGas*RHO*tKelvin
    if(LOUD)then
	    IF (voidFrac < 0 .and. LOUD)print*, 'FuEsd96Vtot:Error! (1-1.9*ETA) IS -VE. eta,rho=',eta,rho
    end if
	!this is not a problem for FuEsdVtot. IF (Z.LT.0.0)pause 'Error! Z NEGATIVE IN FuEsdVtot!'

	!     CALLING ASSYMU IF NUMBER OF DONOR SITES AND ACCEPTOR SITES ARE NOT EQUAL 
	!      IF(IFLAG.EQ.2)CALL ASSYMU(ALPHAD,ALPHDA,XA,XD,X,ND,NC,QIJ,UASSOC)

	IF (IER(2).NE.0 .and. LOUD)WRITE(6,*)'ERROR IN LIQUID PHASE IN ALPSOL'
	IF (IER(3).NE.0 .and. LOUD)WRITE(6,*)'ERROR IN VAPOR PHASE IN ALPSOL'


 	DO I=1,NC
	   YQVI(I)=0.D0
	   CVI(I)=0.D0
	enddo     
	aAssoc=0.d0
	UNUMER=0.d0
	UDENOM=1.d0
	UATT=0.d0
	DO I=1,NC
		!ralph(I)=SQRT(alphAD(I,I))	 ! ralph is computed in alpsolEz.
		aAssoc=aAssoc+xFrac(I)*ND(I)*( 2.d0*DLOG(XA(I,1))+1.d0-XA(I,1) )
		!aAssoc=aAssoc+xFrac(I)*ND(I)*( 2.d0*DLOG(XA(I))+1.d0-XA(I) )
		UATT=UATT+xFrac(I)*VX(I)*EOK(I,I)/tKelvin*(Y(I,I)+K2)*K1(I)       
		UFACTI=xFrac(I)*ND(I)*RALPH(I)*XA(I,1)*(2.d0-XA(I,1))
		!UFACTI=xFrac(I)*ND(I)*RALPH(I)*XA(I)*(2.d0-XA(I))
		UDENOM=UDENOM+xFrac(I)*ND(I)*(RALPH(I)*XA(I,1))**2    
		!UDENOM=UDENOM+xFrac(I)*ND(I)*(RALPH(I)*XA(I))**2    
		bepsADi=DH(I)/tKelvin*TC(I)
		EBETHI=1
		IF(bepsADi > 1.D-5)EBETHI=(EXP(bepsADi)-1.d0)/bepsADi
		DO J=1,NC
			bepsADj=DH(J)/tKelvin*TC(J)
			EBETHJ=1 !anticipate limit for small bepsAD                         
			IF(bepsADj > 1.D-5)EBETHJ=(EXP(bepsADj)-1.d0)/bepsADj
			QIJ=(EXP(DH(J)/tKelvin*TC(J))/EBETHJ+EXP(DH(I)/tKelvin*TC(I))/EBETHI)/2.d0
			!ralph(J)=SQRT(alphAD(J,J))
			UNUMER=UNUMER+UFACTI*xFrac(J)*ND(J)*RALPH(J)*XA(J,1)*QIJ
			!UNUMER=UNUMER+UFACTI*xFrac(J)*ND(J)*RALPH(J)*XA(J)*QIJ
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*xFrac(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*xFrac(J)
 		enddo
	enddo
	FATT= -ZM*YQVM/K1YVM*DLOG(1.d0+K1YVM*RHO)
	FREP= -4.d0/1.9D0*DLOG(voidFrac)*CVM/VM
	FREP= -4.d0/1.9D0*DLOG(voidFrac)*Cmix
	UATT= -ZM*YQVM*RHO/(1+K1YVM*RHO)*UATT/K1YVM	!uAtt=...uAtt/k1YVM where uAtt defined in 3rd line of above do loop.
	!uAssoc= -sum[xiNdiRalphi*XAi*(2-XAi)*xjNdj*ralphj*XAj*(bepsY1_Yadi+bepsY1_Yadj)/2]/[1+sum(xiNdiAlphai*XA^2)]
	UASSOC=-UNUMER/UDENOM	 !uAssoc is no longer computed by Wertheim.

	!     CALLING ASSYMU IF NUMBER OF DONOR SITES AND ACCEPTOR SITES ARE NOT EQUAL 
	!      IF(IFLAG.EQ.2)CALL ASSYMU(ALPHAD,ALPHDA,XA,XD,X,ND,NC,QIJ,UASSOC)

	DUONKT=UATT+UASSOC
	DAONKT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	DSONK =DUONKT-DAONKT
	DHONKT=DUONKT+zFactor-1.d0
	uRes_RT=UATT+UASSOC
	aRes_RT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	sRes_R =UATT+UASSOC-aRes_RT
	hRes_RT=UATT+UASSOC+zFactor-1
	Ares=aRes_RT
	Ures=uRes_RT
	if (isZiter==0) then
        if( zFactor .le. -8686)then  ! Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
		    if(LOUDER)pause 'FuEsdVtot: zFactor.le.0 for fugc calculation.'
			iErr=3
			goto 86
        end if
		if(LOUDER)write(*,'(a,f10.5)')' i,lnGamRep,lnGamAtt,lnGamAsn,ralph.  F=',Fassoc
		DO I=1,NC
			!FUGREP(I)=FREP*( 2.d0*CVI(I)/CVM-VX(I)/VM ) + ZREP*VX(I)/VM
			FUGREP(I)=FREP*( C(I)/Cmix ) + ZREP*VX(I)/VM ! For pure i, FugRepi= -4ci/1.9*ln(1-1.9eta) + 4ci*eta/(1-1.9eta) 
			FUGATT(I)=FATT*( 2*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM )+ZATT*K1(I)*Y(I,I)*VX(I)/K1YVM !91-pres form, complete w/o dNk1YbNk, cf. EL99 p559 & S&E91apx
			FUGASSOC(i)=ND(i)*2*DLOG(XA(i,1)) + zAssoc*1.9D0*VX(i)*rho !JRE'96 Eq.43.
			!FUGASSOC(i)=ND(i)*2*DLOG(XA(i)) + zAssoc*1.9D0*VX(i)*rho !JRE'96 Eq.43.
			FUGC(I)=FUGREP(I)+FUGATT(I)+FUGASSOC(I)  ! -DLOG(Z)  !Don't subtract ln(Z) when given Vtot as independent variable.
			rLnGamRep=FUGREP(i)-Zrep*Vx(i)/VM  ! cf. Bala and Lira (2016), Eqs A6-A14. to correct from constant volume to P.
			rLnGamAtt=FUGATT(i)-Zatt*Vx(i)/VM
			rLnGamAsn=FUGASSOC(i)-Zassoc*Vx(i)/VM
			IF(LOUDER)WRITE(*,'(i3,f7.4,9f10.4)')i,xFrac(i),rLnGamRep,rLnGamAtt,rLnGamAsn,ralph(i)
		ENDDO
	endif
	if(LOUD)write(*,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd96vtot: eta,CHNG,niter,FUGC',eta,0.0,1,(FUGC(i),i=1,NC) 
	if(isZiter < 2)return
	!print*,'aRes_RT=',aRes_RT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Added by AFG 2010
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!JRE: compute these derivatives for thermal props, regardless if isZiter. The are returned as part of GlobalConst
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
	!rhoDZASSOC_dRho=-(-zAssoc*1.9*eta+XA(1)*(1-XA(1))/(2-XA(1))/voidFrac)/voidFrac !(eta/g)*dg/dEta=eta*(1-1.9eta)*1.9/(1-1.9eta)^2=1.9eta/(1-1.9eta)
	!print*,'rhoDZAssoc_dRho 1&2',rhoDZAssoc_dRho,rhoDZAssoc_dRho2
	dP_dRHO=rGas*tKelvin*zFactor+rGas*tKelvin*(RHO*dZREP_dRHO+RHO*dZATT_dRHO+rhoDZASSOC_dRho)!=RT*(Z+rho*dZ_dRho)
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
	cmprsblty=zFactor+zRep/denom+zAtt/(1+K1YVM*RHO) + rhoDZassoc_dRho !cf. EsdDerivatives.jnt
	!print*,'cmprsblty&2=',cmprsblty,cmprsblty2
	dZ_dP=(dZREP_dRHO+dZATT_dRHO+dZASSOC_dRHO)*dRho_dP

	! First derivatives of P,ZREP,ZATT and ZASSOC in respect to T according to the ESD EOS while V and N are constant
	dZREP_dT=0.d0
	dZATT_dT=-ZM*RHO*( dYQVM_dT*(1.d0+K1YVM*RHO)-RHO*YQVM*dK1YVM_dT )/( (1.d0+K1YVM*RHO)*(1.d0+K1YVM*RHO) )
	dZASSOC_dT=-half*(1.9d0*RHO*VM/denom+1)*dh_dT !dh_dT is computed in WertheimFugc.
	dP_dT=rGas*RHO*zFactor+rGas*tKelvin*RHO*(dZREP_dT+dZATT_dT+dZASSOC_dT)
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
		!CvAssoc= uAssoc + bepsY1_Yad**2*XA(1)*(1-XA(1))/(2-XA(1)) + (1-XA(1))*bepsY1_Yad*(1+bepsAD-bepsY1_Yad)
		!print*,'uAssoc,CvAssoc=',uAssoc,CvAssoc
		cvRes_R = uAtt*( 1.7745*eta*beps*exp(beps)/(1+k1YVM*rho) - beps ) + CvAssoc
		!print*,'FuEsdVtot: TdZ_dT',tKelvin*dZ_dT
		!zAssoc= -F^2/(1-1.9eta)=> beps*dZAssoc/dBeps = -2F*(beps*dF/dBeps)/denom
		!F = NdiXAiRalphi)=> beps*dF/dBeps = Ndi*ralphi*XA[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}= F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		!beps*dZAssoc/dBeps = -2F/denom *  F*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]} = 2Zassoc*[(beps/XA)(dXA/dBeps)+(beps/2alpha)*(dAlpha/dBeps)]}
		bxDxa_db= -(1-XA(1,1))/(2-XA(1,1))*bepsY1_Yad
		!bxDxa_db= -(1-XA(1))/(2-XA(1))*bepsY1_Yad
		bepsDZassoc_dBeps=2*zAssoc*(bxDxa_db+bepsY1_Yad/2)
		!write(666,'(a,5f10.4)')'beps,rho,zAssoc,rhoDZassoc_dRho',beps,zAssoc,bepsDZassoc_dBeps

		bepsDZatt_dBeps= -tKelvin*dZatt_dT !see derivative formulas above.
		TdP_dT1=zFactor+tKelvin*dZ_dT
		TdZassoc_dT=tKelvin*dZASSOC_dT
		TdP_dT=zFactor-bepsDZatt_dBeps-bepsDZassoc_dBeps
		!print*,'dZASSOC_dT,dh_dT ',dZASSOC_dT,dh_dT
		if(ABS(cmprsblty) < 1D-11)then
			if(LOUD.and.initCall)Print*,'FuEsdVtot warning: cmprsblty ~ 0'
			cpRes_R = 86.8686D0
		else
			cpRes_R = cvRes_R-1+TdP_dT**2/cmprsblty 
        end if
	endif
86	return
	end	!FuEsd96Vtot

