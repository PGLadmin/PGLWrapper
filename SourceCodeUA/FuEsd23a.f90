!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE EsdParms
	USE GlobConst
	DoublePrecision eokP(nmx),KCSTAR(nmx),DH(nmx),c(nmx),q(nmx),vx(nmx)
	DoublePrecision mShape(nmx),KadNm3(nmx),epsA_kB(nmx),epsD_kB(nmx) !for ESD2
	Integer         ND(nmx),NDS(nmx),NAS(nmx)
	LOGICAL         isMEM2
END MODULE EsdParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine GetEsdCas(NC,idCasPas,iErr) !ID is passed through GlobConst
	!
	!  PURPOSE:  LOOKS UP THE ESD PARAMETERS AND STORES THEM IN USEd EsdParms
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
	if(LOUDER)write(dumpUnit,*)' GetEsdCas: idCas()=',idCas(1:NC) 
	if(LOUDER)write(dumpUnit,*)' GetEsdCas: ID()=',ID(1:NC) 
	if(LOUDER)write(dumpUnit,610)' GetEsdCas: Tc()=',Tc(1:NC)
610 format(1x,a,12E12.4)
	etaMax=1/1.9D0-zeroTol
	isMEM2=.FALSE.
	if(iEosOpt==4)isMEM2=.TRUE.
	
	nComps=NC
	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)

	inFile=TRIM(PGLinputDir)//'\ParmsEsd96.TXT'
	if(iEosOpt==12)inFile=TRIM(PGLinputDir)//'\ParmsEsdEmami.txt' ! // is the concatenation operator
	if(iEosOpt==13)inFile=TRIM(PGLinputDir)//'\ParmsEsdEmamiTb.txt' ! // is the concatenation operator
	if(isMEM2)inFile=TRIM(PGLinputDir)//'\ParmsEsdMEM2.txt' ! // is the concatenation operator
	OPEN(31,FILE=inFile)
	if(LOUDER)write(dumpUnit,*)'GetEsdCas:inFile=',TRIM(inFile)
	if(LOUDER)write(dumpUnit,*) 'Check the ESD parms file location.'


	READ(31,*)nDeck
	if(LOUDER)write(dumpUnit,602)' GetEsdCas: nDeck=',nDeck
	do i=1,nDeck
		READ(31,'(a222)')dumString
		if(isMEM2)then
			READ(dumString,*,ioStat=ioErr)IDA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I)&
			                                                                            ,NDA(i),NDSA(I),NASA(I),idCasa(i)
			!if(LOUDER)write(dumpUnit,602)'From inFile:',IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eDonEpsK(I),eAccEpsK(I)
		else ! iEosOpt=2,12,13 all use ESD96.
			READ(dumString,*,ioStat=ioErr)IDA(I),CAi,QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I) ,idCasa(i)
			!READ(dumString,*,ioStat=ioErr)IDA(I),QA(I) ,eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I),NDA(I),NASA(I),idCasa(i)
			if(NASA(I).ne.NDSA(I))NASA(I)=0	 ! ESD96 requires that NAS=NDS for all compounds.
			NDSA(I)=NASA(I)
			eAccEpsK(i)=DHA(I)*1000/RgasCal
			eDonEpsK(i)=eAccEpsK(i)
			!if(LOUDER)write(dumpUnit,603)'GetEsdCas: inFile~',i,IDCASA(I),QA(I),eokA(I),bVolA(I),KCSTA(I),eAccEpsK(I)
		endif
		!write(dumpUnit,*),*)IDA(I),CA(I),QA(I) ,eokA(I),bVolA(I),NDA(I),KCSTA(I),DHA(I),NASA(I),NDSA(I)  ,idCasa(i)
		if(ioErr .and. LOUDER)write(dumpUnit,'(a,a)')' GetESDCas: error reading ',TRIM(inFile),' line=',TRIM(dumString)
		if(  ( idCasa(i)==id(1) .or. idCasa(i)==id(2) ) .and. LOUDER  )write(dumpUnit,*)'Found in ParmsEsd idCas=',idCasa(i) 
	enddo !i=1,NC
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
601 format(1x,a,8e12.4)
602 format(1x,a,i11,8e12.4)
603 format(1x,a,2i11,8e12.4)
	if(iErr)then
		if(LOUDER)write(dumpUnit,*) 'GetEsdCas: error for at least one compound'
		continue
	endif
        
	if(LOUDER)then
        write(dumpUnit,*)'  ID     NAME       mESD    eok      bVol    Nd   KADnm3   eDon   eAcc(kcal/mol)'
	    do i=1,NC
		    write(dumpUnit,606)IDCas(i),NAME(i),q(i),eokP(i),vx(i),NDS(i),NAS(i),KadNm3(i),&
			                                                 eDonorKcal_mol(i,1),eAcceptorKcal_mol(i,1) !/1000
        enddo
    end if
606	format(i9,1x,a11,f9.3,f8.2,f8.2,2i3,1x,f8.6,2f8.0)

	!note:  bips are passed back through USEd BIPs
	bipFile=TRIM(PGLinputDir)//'\BipEsd96.txt' ! // is the concatenation operator
	if(isMEM2)bipFile=TRIM(PGLinputDir)//'\BipEsdMEM2.txt' ! // is the concatenation operator
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
!------------------------------------------------------------------------------------  
subroutine ExactEsd(NC,vx,c,q,eokP,iErr,ierComp)
	USE GlobConst
	Implicit DoublePrecision(A-H,K,O-Z)
	!  compute ESD2 parameters for hydrocarbons based on the more exact solution
	!  for Zc, Bc, and Yc vs. cShape
	!  Ref:  Elliott and Lira, Introductory Chemical Engineering Thermo, p564 (1999).
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
		if(isHelium)then
			iErr=3
			if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for helium ID=',ID(i)
			ierComp(i)=3
			cycle
		endif	
		isH2=0
		if( id(i)==902 .or. ID(i)==133740 .or. id(i)==925 .or. id(i)==7782390)isH2=1
		if(isH2)then
			iErr=2
			if(LOUDER)write(dumpUnit,*)'ExactESD: Parms not available for H2 or D2 ID=',ID(i)
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
		q(i) = qShape
		eokP(i)=TC(i)*rlnY1
		bVolCc_mol(i)=8.314*Tc(i)/Pc(i)*Bc
		vx(i)=bVolCc_mol(i)
	enddo
	return
end	!subroutine ExactEsd
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
	if(eta > 1/1.9 .and. LOUDER)write(dumpUnit,*)'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
	isZiter=1 ! FUGC calculations are skipped for isZiter=1
	Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
	IF(iErr > 10)GOTO 86 ! let iErr=iErrF in FugiTP
	etaOld=eta
	ERROLD=Pb_RT-eta*zFactor

	eta=etaOld/1.15D0
	IF (eta < 0 .and. LOUD) write(dumpUnit,31)LIQ
	rho=eta/bMix
	if(initial.and.LOUD)write(dumpUnit,*)'FugiEsd: initial eta,err',etaOld,errOld
	itMax=77
	errBesteta=1234
	do nIter=1,itMax
		Call FuEsdVtot(isZiter,tKelvin,bMix/eta,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
		IF(iErr > 10)EXIT
		ERR=Pb_RT-eta*zFactor
		CHNG=ERR/(ERR-ERROLD)*(eta-etaOld)
		if(initial.and.LOUDER)write(dumpUnit,'(a,2e11.4,3f10.5)')'FugiEsd eta,Z', eta,zFactor
		if(initial.and.LOUDER)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd eta,CHNG,niter',eta,CHNG,niter 
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
	if(initial.and.LOUD)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')' FuEsd2 cnvrgd: eta,CHNG,niter',eta,CHNG,niter
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
	USE GlobConst
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...	 and FugiParts(rLnGam...)
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	Implicit DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(nmx),xFrac(nmx),FUGC(nmx),k10,k2,Zm !,KCSTARp(nmx),KVE(nmx)
	DoublePrecision YQVIJ(nmx,nmx),YQVI(nmx),Y(nmx,nmx),EOK(nmx,nmx)
	DoublePrecision CVI(nmx),CVIJ(nmx,nmx),QV(nmx,nmx)
    DoublePrecision voidFrac,tKelvin,zAssoc,aAssoc,uAssoc,rho,fAssoc,zFactor
	Integer initCall
	DoublePrecision k1(nmx) 
	LOGICAL LOUDER
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),Zrep,Zatt,Zassoc,aRep,aAtt,aAssoc,uAtt,uAssoc
	!common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	DATA K10,K2,ZM,initCall/1.7745D0,1.0617D0,9.5D0,1/

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
	if(LOUDER.and.initCall)write(dumpUnit,601)' FuEsdVtot: T,x1,bMix,eta=',tKelvin,xFrac(1),bMix,eta
601 format(1x,a,8E12.4)
	YQVM=0.d0
	VM=0.d0
	CVM=0.d0
	Cmix=0.d0
	K1YVM=0
	iType=1	  
	DO I=1,NC
		DO J=1,NC
			kIJbip=KIJ(I,J)+KTIJ(I,J)/tKelvin 
			EOK(I,J)=DSQRT(eokP(I)*eokP(J))*(1.d0-kIJbip)
			Y(I,J)=DEXP(EOK(I,J)/tKelvin)-K2
			QV(I,J) = (Q(I)*bVolCc_mol(J) + Q(J)*bVolCc_mol(I)) / 2.d0
			YQVIJ(I,J)=QV(I,J)*Y(I,J)
			CVIJ(I,J) = (C(I)*bVolCc_mol(J) + C(J)*bVolCc_mol(I)) / 2.d0 
			! e.g. (x1*c1+x2*c2)*(x1*b1+x2*b2) = x1^2*c1*b1+x1*x2*(c1*b2+c2*b1)+x2^2*b2^2
			YQVM=YQVM+YQVIJ(I,J)*xFrac(I)*xFrac(J)	   
			CVM = CVM + CVIJ(I,J)*xFrac(I)*xFrac(J)	   !note: above means <c>=sum(xi*ci) and <b>=sum(xj*bj) 
		enddo
		k1(I)=K10
		Cmix=Cmix+xFrac(i)*C(i)
		VM=VM+xFrac(I)*bVolCc_mol(I)
		K1YVM=K1YVM+xFrac(I)*K1(I)*Y(I,I)*bVolCc_mol(I) !1991 form, overwritten if applying 1990 form
	enddo
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
	voidFrac=1-1.9D0*eta
	denom=voidFrac
	zRep= 4.d0*Cmix*eta/denom
	zAtt= -ZM*YQVM*rho/(1+K1YVM*rho)
	zFactor=(1+zRep+zAtt+zAssoc)
	pMPa=zFactor*Rgas*rho*tKelvin
    if(voidFrac < 0)then
	    IF(LOUD)write(dumpUnit,*) 'FuEsdVtot:Error! (1-1.9*eta) IS -VE. eta,rho=',eta,rho
		iErr=14
    endif
	if(LOUDER)write(dumpUnit,form610)' FuEsdVtot: done with isZiter=1. zAssoc,zFactor=',zAssoc,zFactor
	if(iErr > 10)return


	!aRep= -4.d0/1.9D0*DLOG(voidFrac)*CVM/VM
	aRep= -4.d0/1.9D0*DLOG(voidFrac)*Cmix
	aAtt= -ZM*YQVM/K1YVM*DLOG(1.d0+K1YVM*rho)
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
		BdYb_dB=BdYb_dB+xFrac(I)*vx(I)*EOK(I,I)/tKelvin*(Y(I,I)+K2)*K1(I)  ! = Beta*d<k1Yb>/dBeta     
		DO J=1,NC
			BdYbq_dB=BdYbq_dB+xFrac(J)*xFrac(i)*QV(i,j)*EOK(i,j)/tKelvin*(Y(i,j)+K2) ! = Beta*d<Ybq>/dBeta
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*xFrac(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*xFrac(J)
 		enddo
	enddo
	UATT= -ZM*YQVM*rho/(1+K1YVM*rho)*BdYb_dB/K1YVM + aAtt*(BdYbq_dB/YQVM-BdYb_dB/K1YVM) !FYI:I had omitted the 2nd term previously

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
			fugAtt(i)=aAtt*( 2*YQVI(I)/YQVM-K1(I)*Y(I,I)*vx(I)/K1YVM )+zAtt*K1(I)*Y(I,I)*vx(I)/K1YVM !91-pres form,
			!fugAssoc(i)=ND(i)*2*DLOG(XA(i,1)) + zAssoc*1.9D0*vx(i)*rho !JRE'96 Eq.43.
			FUGC(I)=FUGREP(I)+FUGATT(I)+fugAssoc(I)  ! -DLOG(Z)  Don't subtract ln(Z) when given Vtot as independent variable.
			rLnGamRep(i)=FUGREP(i)-zRep*vx(i)/VM  ! cf. Bala and Lira (2016), Eqs A6-A14. to correct from constant volume to P.
			rLnGamAtt(i)=FUGATT(i)-zAtt*vx(i)/VM
			rLnGamAssoc(i)=fugAssoc(i)-zAssoc*vx(i)/VM
			IF(LOUDER)write(dumpUnit,'(i3,f7.4,9f10.4)')i,xFrac(i),rLnGamRep(i),rLnGamAtt(i),rLnGamAssoc(i) !,ralpha(i),ralphd(i)
		ENDDO
	endif
	!I've lost faith in uAssoc for MEM2. Differentiate numerically.
	if(isTPT)then  ! disables because isESD =/= isTPT. Change to isMEM2 if you want to enable. 
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

