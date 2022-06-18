	SUBROUTINE FuEsdTp(tKelvin,P,xFrac,NC,LIQ,FUGC,Z,ier)
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS	 + GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),IER(12) !,dHkcalMol(NMX),KCSTARp(NMX)
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	COMMON/CONSTK/KdumIJ(NMX,NMX),INITIAL 
	COMMON/ALPHA/ALPHAD(NMX,NMX), ALPHDA(NMX,NMX)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/eta2/eta
	DATA INITIAL/0/  

	!DATA K10,K2,ZM,INITIAL/1.7745,1.0617,9.5,0/  
	!	C$FUGI.FOR
	!   LATEST REVISION : 9/94 jre
	!                   : 1/96 (switched to chempot, ADDED POLYETHYLENE  jre)
	!                   : 7/96 PS, PPO, PEO, PIB (ram natarajan)
	!                   : 1/97 PS, PPO, PEO, PIB (made consistent, jre)
	!                   : 3/97 Added ASSYM (Marty Yuzwa) ROUTINE REQUIRES
	!					:	   THE USE OF A NEW MAINESD, DNEQNF, AND ESDPARMS.DAT
	!	NAS IS THE NUMBER OF ACCEPTER SITES
	!	NDS IS THE NUMBER OF DONOR SITES
	!	ND IS THE DEGREE OF POLYMERIZATION OF EACH SPECIES
	!	EOKP IS THE DISPERSE ATTRACTION OVER BOLTZ k FOR PURE SPECIES
	!	KCSTAR IS THE DIMENSIONLESS BONDING VOLUME FOR PURE 
	!	DH IS THE REDUCED BONDING ENERGY
	!	C,Q,VX ARE THE PURE COMPONENT EOS PARAMETERS
	!	KIJ IS THE BINARY INTERACTION COEFFICIENT 
	!	Z IS PV/NoKT HERE  
	!	IER = 1 - AT LEAST ONE ERROR
	!	      2 - TOO MANY ALPHA ITERATIONS
	!	      3 - 
	!	      4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!	      5 - 
	!	      6 - TOO MANY Z ITERATIONS
	!	Literature:
	!	[1] ESD, IECR, 29:1476 (1990) Note: <Yb>=<qYb>/<q> was superseded in ref[2]apx.
	!	[2] S&E, IECR, 30:524  (1991) Note:	Many typos in the apx here make it worthless. W1,W2 approach superseded in ref[3]
	!	[3] S&E, IECR, 31:2783 (1992) Note: <Yb> mix rule not clarified here, but ref[2]apx clarifies.
	!	[4] P&E, IECR, 32:3174 (1993) Note: <Yb> mix rule is wrong here.  Copied from 1990, not from the program.
	!	[5] JRE, IECR, 35:1624 (1996) Note: a typo of -1 was omitted then canceled. pdf clarifies.
	!	[6] E&N, IECR, 41:1043 (2002) 
	!	Example 1.  nC7+benzene at 458.1K,0.9MPa, Hij=0, Kij=0,
	!	id  xi  	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA 	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt		fugc
	!	17  0.6439	 47.8	2.30	280.7	0.000	0.000	0	0	1	1   	-   	6.245970544	-9.676977412   -0.3271
	!	501 0.3561	 29.5	1.77	336.5	0.000	0.000	0	0	1	1   	-     	4.254184969	-7.290251086	0.0682
	!	vm  	cShapeMix	cvm 	qYbMix	k1YbMix	etaLiq 	zRep	zAtt	   zFactor
	!	41.27	2.11096846	87.128	108.587	61.8619	.21845	3.15351	-4.10887   0.045
	!	Example 2.  EtOH+H2O at 363K,0.1MPa, Hij=0, Kij=0.0323 ->etaLiq=0.3441,zLiq=0.001588,k1Yb=40.04375
	!	id  xi	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD	lnPhiAssoc	lnPhiRep	lnPhiAtt		fugc
	!	1102	0.5	23.574	1.1565	270.14	4.985	0.0283	1	1	0.144	0.144	-6.223	10.8003		-10.4607	0.56175
	!	1921	0.5	 9.411	1.0053	427.25	5.143	0.10	1	1	0.113	0.113	-5.301	5.16093		-6.30716   -0.00227
	!	T&E=>  eta=0.3441; Z=0.00158; aAssoc=-3.245
	!          dXAidNj				dXDidNj
	!	 -0.2429	-0.1383 	-0.2429	-0.1383
	!	 -0.1970	-0.1121 	-0.1970	-0.1121
	!	Example 3.  Acetone(nonAssoc)+chloroform at 338K,0.101325MPa, 
	!	T&E=>  eta=0.305,aAssoc=-0.203
	!	id	xi	  c		  eokp	   b	Nd	  KcStar	eHb/RTc		NAS	NDS	XA	  XD	lnPhiAssoc
	!	1051	0.36 2.1001	274.703	26.965	1	1.00E-01	0.51    	 1	 1	0.63 0.98	 -0.680
	!	1521	0.64 1.7949	318.313	25.643	1	0.0143    	6.550   	 0	 1	 1   0.80	 -0.405
	!          dXAidNj				dXDidNj
	!        -0.19	-0.65		-0.073	-0.622
	!	      0		  0			-0.006	-0.050
	!	bMix	cshapemix	cbMix	qYbMix	k1YbMix	etaLiq	zRep	zAtt	zAssoc	sqrt(alpha1)	fAssoc	kbe(1)
	!	34.044	1.883551	64.1246	107.349	71.1151	0.30995	5.68055	-5.6358	-1.0410		4.701303	0.65418	998.0017
	!	Example 4.  MeOH+Benzene at 331.08K,0.1MPa, Hij=0, Kij=0.0182 
	!	id  	xi	   b	   c   	eokp	eHb(kcal)KcStar	NAS NDS	XA  	XD  	lnPhiRep	lnPhiAtt	lnPhiAssoc	
	!	1101 0.5	 20.4	1.12	326.1	5.17	0.0226	1	1	0.1878	0.1878	6.175   	-8.0857 	-3.87025
	!	501  0.5	 29.5	1.77	336.5	0.000	0.000	0	0	1	     1  	9.248   	-14.543     -0.7625
	!	bMix	cbMix	qYbMix	k1YbMix	etaLiq	zRep	zAtt	zAssoc	sqrt(alpha1)	fAssoc	kbe(1)
	!	24.95	36.04	75.85	73.77	0.3228	4.823	-4.7697	-1.0501		6.787   	0.6373	1377

	do I=1,6
		IER(I)=0
	enddo	
	!  ON INITIAL CALL, CHECK FOR POLYMERS THAT NEED MW SPECS TO GET C,VX,ETC	
	!IF(INITIAL.EQ.0)CALL POLYCBE(NC,ID,C,Q,VX,EOKP)
	INITIAL=1
	totMoles=1 !we may assume that totMoles=1, since sum(xFrac)=1
	bMix=0
	TCt=0
	DO I=1,NC
        if(LOUD)then
		    IF(xFrac(I).LT.0)PAUSE 'ERROR IN FUGI, Xi<0'
        end if
		gmol(i)=xFrac(i) !we may assume that totMoles=1
		bMix=bMix+xFrac(i)*vx(i)
    enddo

    if(LOUD)then
	    if(P < 1E-14)pause 'FuEsdTp: PMPa ~ 0???'
	    if(tKelvin < 1)pause 'FuEsdTp: T(K) < 1???'
    end if
	PORT=P/RGAS/tKelvin
	PVMORT=PORT*bMix

	!     GUESS FOR RHO
	ETA=PVMORT
	IF(LIQ.EQ.1.or.liq.eq.3)ETA=.52
	RHO=ETA/bMix
	vTotCc=totMoles/rho

	!  ALPSOL WILL USE THIS VALUE OF RHO TO CALCULATE A NEW VALUE OF RHO
	isZiter=0
	call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)

	vOld=vTotCc
	ERROLD=PORT-Z*totMoles/vTotCc
	vTotCc=vTotCc*1.05
	IF (vTotCc.LT.0) WRITE(6,31)NC
	NITER=0
	relErr=111
	iErrCode=0
	do while(relErr.gt.1E-9.and.iErrCode.eq.0)
		NITER=NITER+1
		if(NITER.ge.113)iErrCode=6
		isZiter=1
10		call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot) !dFuDp,dFuDv,dFuDn,IER)
		err=PORT-Z*totMoles/vTotCc
        if(LOUD)then
		    if(err.eq.errOld)pause 'FuEsdTp: err=errOld during Z iteration'
        end if
		CHNG=err/(err-errOld)*(vTotCc-vOld)
		if(ABS(CHNG/vTotCc).gt.0.3)CHNG=vTotCc*DSIGN(0.2d0,CHNG)
		vOld=vTotCc
		errOld=err
		vTotCc=vTotCc-CHNG
		if(vTotCc.le.0) then
			vTotCc=-0.95d0*vTotCc
			goto 10
		endif !pause 'FuEsdTp: vTotCc.le.0 during vTotCc iteration'
		relErr=ABS(chng/vTotCc)
	enddo
	!if(iErrCode.ne.0)then
		!c		call GoldenZ2(X,nC,KVE,ND,VM,CVM,qYbMix,k1YbMix,PORT,LIQ,
		!c     &                         Z,rho,zAssoc,XA,RALPH,fAssoc,ier)
		!c		if(LOUD)write(*,*)'liq,GoldenZ',liq,z
		!ier(11)=1
		!ier(1)=1
	!endif
	if(ier(2).ne.0.or.ier(3).ne.0)ier(1)=1
	if(ier(1).ne.0)goto 86

	!     if(LOUD)write(*,*)'ETA,Z',ETA,Z
	!      if(LOUD)write(*,*)'ZREP,ZATT',ZREP,ZATT
	IF (vTotCc.LT.0) WRITE(6,31)LIQ

	!  ITERATION ON RHO HAS CONCLUDED.  GET DEPARTURES AND FUGACITY COEFFS.
	isZiter=0
	call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
	eta=bMix/vTotCc
	if (LIQ.EQ.1.or.LIQ.eq.3) then
		etaL=eta
		zL=Z
	else
		etaV=eta
		zV=Z
	endif

	RETURN
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! RHO (1) -VE IN FUGI')
33	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! RHO (2) -VE IN FUGI')
86	continue
	IF(IER(4)==1.and.LOUD) WRITE(*,*)'ERROR IN ALPSOL'
	IER(1)=1
	iErr=1
686	FORMAT(' ERROR IN FuEsdWert.  ')
	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C	Written by AFG, Oct. 2009																				C
!C	With Having T,V and nMols, this routine calculates the Z factor and all derivatives needed in 			C
!C	critical point,	flash and bubble point calculations.													C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErr)
	USE Assoc !includes GlobConst {Tc,Pc,...} + XA,XD,XC...
	USE EsdParms ! eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),KCSTARp(NMX),IER(12)
	DIMENSION YQVIJ(NMX,NMX),KVE(NMX),YQVI(NMX),Y(NMX,NMX),EOK(NMX,NMX)
	DIMENSION CVI(NMX),CVIJ(NMX,NMX),QV(NMX,NMX)
	DIMENSION RALPH(NMX),k1(NMX)
	DIMENSION FUGASSOC(NMX)
	double precision moleStep
	DIMENSION dEOK_dT(NMX,NMX),dY_dT(NMX,NMX),dYQVI_dT(NMX),dYQVIJ_dT(NMX,NMX)
	DIMENSION dCVM_dN(NMX),dYQVM_dN(NMX),dVM_dN(NMX),dK1YVM_dN(NMX),dCVI_dN(NMX,NMX)
	DIMENSION dZREP_dN(NMX),dZATT_dN(NMX),dZASSOC_dN(NMX),dZ_dN(NMX),dP_dN(NMX),dFREP_dN(NMX),dFATT_dN(NMX),dYQVI_dN(NMX,NMX)
	DIMENSION dFUGC1_dP(NMX),dFUGC1_dT(NMX),dFUGC2_dP(NMX),dFUGC2_dT(NMX),dFUGASSOC_dRHO(NMX)
	DIMENSION dFUGREP_dN(NMX,NMX),dFUGATT_dN(NMX,NMX),dFUGC1_dN(NMX,NMX),dFUGC2_dN(NMX,NMX),dFUG_dN(NMX,NMX)
	DIMENSION dh_dN(NMX),dFUGASSOC_dT(NMX),dFUGASSOC_dN(NMX,NMX),FUGREP(NMX),FUGATT(NMX)
	DIMENSION dh_dN_num(NMX),dFUGASSOC_dN_num(NMX,NMX),fugassocLoop(NMX),gmol_old(NMX)
	DIMENSION vMolecNm3(NMX) !for Wertheim
	
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
	COMMON/Helm_derv/Ur,UrV,FREP,FATT,FASSOC
	DATA K10,K2,ZM,initCall/1.7745D0,1.0617D0,9.5D0,1/
	  
	iErr=0
	zero=1D-11
	totMoles=sum(gmol)
	if( tKelvin	   < zero .or. totMoles < zero .or. vTotCc < zero)then
		if(LOUD)print*,'FuEsdVtot: nonsense T(K),totMoles,vTotCc=',tKelvin,totMoles,vTotCc
		iErr=1
	endif
	bMix=0
	do i=1,nc
		xFrac(i)=gmol(i)/totMoles
		IF(xFrac(i) < 0)then
			if(LOUD)PAUSE 'ERROR IN FuEsdVtot, Xi<0'
			iErr=1
		endif
		bMix=bMix+xFrac(i)*bVolCC_mol(i)
	enddo
	if(iErr)return
	rho=totMoles/ vTotCc
	eta=rho*bMix 
	if(LOUD.and.initCall)print*,'FuEsdVtot: bMix,eta=',bMix,eta
	YQVM=0.d0
	VM=0.d0
	CVM=0.d0
	K1YVM=0
	iType=1	  
	DO I=1,NC
		K1(I)=K10
		DO J=1,NC
			kIJbip=KIJ(I,J)+KTIJ(I,J)/tKelvin 
			EOK(I,J)=DSQRT(EOKP(I)*EOKP(J))*(1.d0-kIJbip)
			Y(I,J)=DEXP(EOK(I,J)/tKelvin)-K2
			QV(I,J) = (Q(I)*VX(J) + Q(J)*VX(I)) / 2.d0
			YQVIJ(I,J)=QV(I,J)*Y(I,J)
			CVIJ(I,J) = (C(I)*VX(J) + C(J)*VX(I)) / 2.d0
			YQVM=YQVM+YQVIJ(I,J)*xFrac(I)*xFrac(J)	   
			CVM = CVM + CVIJ(I,J)*xFrac(I)*xFrac(J)	   !note: above means <c>=sum(xi*ci) and <v>=sum(xi*vi)
		enddo
		dHkcalMol(I)=DH(I)*TC(I)*1.987D-3	!JRE 12/00
		KVE(I)=KCSTAR(I)*( DEXP(dHkcalMol(I)/tKelvin/1.987D-3)-1.d0 )  *avoNum  !KVE[=] cc/mol
		!	if(initCall.and.LOUD)print*,'FuESDVtot: KcStar=',KcStar(i)
		!if(initCall.and.LOUD)pause 'Right?'
		KCSTARp(I)=KCSTAR(I)
		VM=VM+xFrac(I)*VX(I)
		K1YVM=K1YVM+xFrac(I)*K1(I)*Y(I,I)*VX(I) !1991 form, overwritten if applying 1990 form
		vMolecNm3(i)=VX(I)/avoNum
		bondVolNm3Esd(I)=KCSTAR(I) !*vMolecNm3(I)
		eHbKcal_mol(I,iType)=dHkcalMol(I)
		bondVolNm3(I,iType)=bondVolNm3Esd(I)
	enddo
	if( ABS(VX(1) - bVolCC_mol(1)) > zero .and. LOUD ) print*,'FuEsdVtot: VX.ne.bVol=',VX(1),bVolCC_mol(1)
	!endif	!this was not working when computing vp
	if(LOUD.and.k1yvm < zero)print*,'FuEsdVtot: 0~k1yvm=',k1yvm 
	eta=rho*vm
	call WERTHEIM(vMolecNm3,eta,tKelvin,xFrac,NC,zAssoc,aAssoc,uAssoc,iErrCode)

	IF(iErrCode.ne.0)GOTO 86
	voidFrac=1-1.9D0*ETA
	denom=voidFrac
	ZREP= 4.d0*CVM*RHO/denom
	ZATT= -ZM*YQVM*RHO/(1.d0+K1YVM*RHO)
	Z=(1.d0+ZREP+ZATT+ZASSOC)
	PMpa=Z*rGas*RHO*tKelvin
    if(LOUD)then
	    IF (voidFrac < 0 .and. LOUD)print*, 'FuEsdVtot:Error! (1-1.9*ETA) IS -VE. eta,rho=',eta,rho
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
		ralph(I)=SQRT(alphAD(I,I))
		aAssoc=aAssoc+xFrac(I)*ND(I)*( 2.d0*DLOG(XA(I,1))+1.d0-XA(I,1) )
		UATT=UATT+xFrac(I)*VX(I)*EOK(I,I)/tKelvin*(Y(I,I)+K2)*K1(I)       
		UFACTI=xFrac(I)*ND(I)*RALPH(I)*XA(I,1)*(2.d0-XA(I,1))
		UDENOM=UDENOM+xFrac(I)*ND(I)*(RALPH(I)*XA(I,1))**2    
		bepsADi=DH(I)/tKelvin*TC(I)
		EBETHI=1
		IF(bepsADi.GT.1.D-5)EBETHI=(EXP(bepsADi)-1.d0)/bepsADi
		DO J=1,NC
			bepsADj=DH(J)/tKelvin*TC(J)
			EBETHJ=1 !anticipate limit for small bepsAD                         
			IF(bepsADi.GT.1.D-5)EBETHJ=(EXP(bepsADj)-1.d0)/bepsADj
			QIJ=(EXP(DH(J)/tKelvin*TC(J))/EBETHJ+EXP(DH(I)/tKelvin*TC(I))/EBETHI)/2.d0
			ralph(J)=SQRT(alphAD(J,J))
			UNUMER=UNUMER+UFACTI*xFrac(J)*ND(J)*RALPH(J)*XA(J,1)*QIJ
			YQVI(I)=YQVI(I)+YQVIJ(I,J)*xFrac(J)
			CVI(I)=CVI(I) + CVIJ(I,J)*xFrac(J)
 		enddo
	enddo
	FATT=-ZM*YQVM/K1YVM*DLOG(1.d0+K1YVM*RHO)
	FREP=-4.d0/1.9D0*DLOG(voidFrac)*CVM/VM
	UATT=-ZM*YQVM*RHO/(1+K1YVM*RHO)*UATT/K1YVM	!uAtt=...uAtt/k1YVM where uAtt defined in 3rd line of above do loop.
	!uAssoc= -sum[xiNdiRalphi*XAi*(2-XAi)*xjNdj*ralphj*XAj*(bepsY1_Yadi+bepsY1_Yadj)/2]/[1+sum(xiNdiAlphai*XA^2)]
	!UASSOC=-UNUMER/UDENOM	 !uAssoc is computed by Wertheim.

	!     CALLING ASSYMU IF NUMBER OF DONOR SITES AND ACCEPTOR SITES ARE NOT EQUAL 
	!      IF(IFLAG.EQ.2)CALL ASSYMU(ALPHAD,ALPHDA,XA,XD,X,ND,NC,QIJ,UASSOC)

	DUONKT=UATT+UASSOC
	DAONKT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	DSONK =DUONKT-DAONKT
	DHONKT=DUONKT+Z-1.d0
	uRes_RT=UATT+UASSOC
	aRes_RT=FREP+FATT+aAssoc !-DLOG(Z) !don't subtract log(z) for aRes)T,V. Important for EAR.
	!print*,'aRes_RT=',aRes_RT
	sRes_R =UATT+UASSOC-aRes_RT
	hRes_RT=UATT+UASSOC+Z-1

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
	fAssoc= sqrt(-zAssoc*voidFrac)		!Zassoc = -F^2/voidFrac
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


	dP_dV=dP_dRHO*RHO*(-1.d0/vTotCc)
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
	if (isZiter.eq.0) then
		call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
        if(LOUD)then
		    if(Z.le.0)pause 'FuEsdVtot: Z.le.0 for fugc calculation.'
        end if
		DO I=1,NC
			FUGREP(I)=FREP*( 2.d0*CVI(I)/CVM-VX(I)/VM ) + ZREP*VX(I)/VM
			FUGATT(I)=ZATT*K1(I)*Y(I,I)*VX(I)/K1YVM+FATT*( 2.d0*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM ) !91-pres form, complete w/o dNk1YbNk, cf. EL99 p559 & S&E91apx
			FUGC(I)=FUGREP(I)+FUGATT(I)+FUGASSOC(I)-DLOG(Z)
		ENDDO

		two=2.d0
		moleStep=100000.d0
		rho_old=rho
		totMoles_old=totMoles
		vm_old=vm
		eta_old=eta
		DO I=1,NC
			gmol_old(I)=gmol(I)
			gmol(I)=gmol(I)+gmol(I)/moleStep
			totMoles=0.d0
			vm=0.d0
			do ii=1,nc
				totMoles=totMoles+gmol(ii)  
			enddo
			do ii=1,nc
				xFrac(ii)=gmol(ii)/totMoles
				vm=vm+xFrac(ii)*VX(ii)
			enddo
			rho=totMoles/vTotCc
			eta=rho*vm
			call WERTHEIM(vMolecNm3,eta,tKelvin,xFrac,NC,zAssoc,aAssoc,uAssoc,iErrCode)
			call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
			DO J=1,NC
				fugassocLoop(J)=FUGASSOC(J)
			ENDDO
			hLoop=h_nMichelsen
			gmol(I)=gmol_old(I)-gmol_old(I)/moleStep
			totMoles=0.d0
			vm=0.d0
			do ii=1,nc
				totMoles=totMoles+gmol(ii)  
			enddo
			do ii=1,nc
				xFrac(ii)=gmol(ii)/totMoles
				vm=vm+xFrac(ii)*VX(ii)
			enddo
			rho=totMoles/vTotCc
			eta=rho*vm
			call WERTHEIM(vMolecNm3,eta,tKelvin,xFrac,NC,zAssoc,aAssoc,uAssoc,iErrCode)
			call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
			gmol(I)=gmol_old(I)
			dh_dN_num(I)=(hLoop-h_nMichelsen)/(two*gmol(I)/moleStep)
			DO J=1,NC
				dFUGASSOC_dN_num(J,I)=((fugassocLoop(J)-FUGASSOC(J))/(two*gmol(I)/moleStep))*totmoles
 			ENDDO
 		ENDDO
		rho=rho_old
		eta=eta_old
		totmoles=totmoles_old
		vm=vm_old
		do ii=1,nc
			xFrac(ii)=gmol(ii)/totMoles
		enddo
 		call WERTHEIM(vMolecNm3,eta,tKelvin,xFrac,NC,zAssoc,aAssoc,uAssoc,iErrCode)
		call WertheimFugc(xFrac,vMolecNm3,tKelvin,NC,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)

		one=1.d0
		!half=one/2.d0	!defined in GlobConst 9/24/19
	! First derivatives of P,ZREP,ZATT and ZASSOC in respect to N according to the ESD EOS while V and T are constant
		DO I=1,NC
			dZREP_dN(I)=4.d0*RHO*(dCVM_dN(I)-1.9d0*VM*dCVM_dN(I)*RHO+CVM+1.9d0*dVM_dN(I)*CVM*RHO)/(denom*denom)
			dZATT_dN(I)=-ZM*RHO*(dYQVM_dN(I)+dYQVM_dN(I)*RHO*K1YVM+YQVM-YQVM*RHO*dK1YVM_dN(I))/((1.d0+K1YVM*RHO)*(1.d0+K1YVM*RHO))
			dZASSOC_dN(I)=-half*(dh_dN_num(I)*(dLng+1.d0)+h_nMichelsen*(d2Lng*VX(I)/VM))
			dZ_dN(I)=dZREP_dN(I)+dZATT_dN(I)+dZASSOC_dN(I)
			dP_dN(I)=RHO*rGas*tKelvin*(Z+dZ_dN(I))
		ENDDO
	! first derivative in respect to N while T and V are constant
		DO I=1,NC
			dYQVI_dT(I)=0.d0     
			DO J=1,NC
				dYQVI_dT(I)=dYQVI_dT(I)+QV(I,J)*xFrac(J)*dY_dT(I,J)
				dCVI_dN(I,J)=CVIJ(I,J)-CVI(I)
				dYQVI_dN(I,J)=YQVIJ(I,J)-YQVI(I)
			ENDDO
		ENDDO
	! First Derivative in respect to RHO while N and T are constant
		dFREP_dRHO=-4.d0*CVM/((1-1.9d0*RHO*VM)*1.9d0)
		dFATT_dRHO=-ZM*YQVM/(1.d0+K1YVM*RHO)

	! First Derivative in respect to T while N and V are constant
		dFATT_dT=-ZM*( (dYQVM_dT*K1YVM-dK1YVM_dT*YQVM)/(K1YVM*K1YVM)*DLOG(1.d0+K1YVM*RHO)+YQVM*RHO*dK1YVM_dT/(K1YVM*(1.d0+K1YVM*RHO)) )

	! First Derivative in respect to N while V and T are constant
		DO I=1,NC
			dFREP_dN(I)=4.d0*VX(I)*RHO*(CVM/VM)/(1.d0-1.9d0*RHO*VM)-4*(DLOG(1.d0-1.9d0*RHO*VM)/1.9d0)*(dCVM_dN(I)*VM-dVM_dN(I)*CVM)/(VM*VM)
			dFATT_dN(I)=-ZM*( ( (dYQVM_dN(I)*K1YVM-dK1YVM_dN(I)*YQVM)/(K1YVM*K1YVM) )*DLOG(1+K1YVM*RHO)+(dK1YVM_dN(I)*RHO+RHO*K1YVM)/(1+K1YVM*RHO)*(YQVM/K1YVM) )
		ENDDO

		DO I=1,NC
		! First derivatives in respect to P while T and N are constant
			dFUGREP_dP=(dFREP_dRHO*(2.d0*CVI(I)/CVM-VX(I)/VM )+dZREP_dRHO*VX(I)/VM)*dRHO_dP
			dFUGATT_dP=(dZATT_dRHO*K1(I)*Y(I,I)*VX(I)/K1YVM+dFATT_dRHO*( 2.d0*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM ))*dRHO_dP
			dFUGBON_dP=dFUGASSOC_dRHO(I)*dRho_dP
			dFUGC1_dP(I)=((dP_dN(I)/dP_dRHO)*vTotCc*Z-1.d0)*(1.d0/PMpa)	  !called FP in Michelsen's book, P.30

		! First derivatives in respect to T while V and N are constant
			dFUGREP_dT=0.d0
			SecTerm=dFATT_dT*( 2.d0*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM )+FATT*( (2.d0*dYQVI_dT(I)*YQVM-2.d0*dYQVM_dT*YQVI(I))/(YQVM*YQVM)-K1(I)*VX(I)*(dY_dT(I,I)*K1YVM-dK1YVM_dT*Y(I,I))/(K1YVM*K1YVM) )
			dFUGATT_dT=K1(I)*VX(I)*( ((dZATT_dT*Y(I,I)+dY_dT(I,I)*ZATT)*K1YVM-dK1YVM_dT*ZATT*Y(I,I))/(K1YVM*K1YVM) )+SecTerm
			dFUGC1_dT(I)=dFUGREP_dT+dFUGATT_dT+dFUGASSOC_dT(I)-dZ_dT/Z

		! First derivative in respect to T while P and N are constant
			dFUGC2_dT(I)=dFUGC1_dT(I)-dFUGC1_dP(I)*dP_dT				!called FT in Michelsen's book, P.30

		! First derivative in respect to P while V and N are constant
			dFUGC2_dP(I)=dFUGC1_dP(I)+dFUGC2_dT(I)*dT_dP

		! First derivative in respect to N while V and T are constant
			DO J=1,NC
				dFUGREP_dN(I,J)=dFREP_dN(J)*(2.d0*CVI(I)/CVM-VX(I)/VM)+FREP*(2.d0*(dCVI_dN(I,J)*CVM-dCVM_dN(J)*CVI(I))/(CVM*CVM)+dVM_dN(J)*VX(I)/(VM*VM))+dZREP_dN(J)*VX(I)/VM-ZREP*(VX(I)*dVM_dN(J))/(VM*VM)
				dFUGATT_dN(I,J)=K1(I)*Y(I,I)*VX(I)*(dZATT_dN(J)*K1YVM-dK1YVM_dN(J)*ZATT)/(K1YVM*K1YVM)+dFATT_dN(J)*( 2.d0*YQVI(I)/YQVM-K1(I)*Y(I,I)*VX(I)/K1YVM )+FATT*(2.d0*(dYQVI_dN(I,J)*YQVM-YQVI(I)*dYQVM_dN(J))/(YQVM*YQVM)+K1(I)*Y(I,I)*VX(I)*dK1YVM_dN(J)/(K1YVM*K1YVM))
				dFUGC1_dN(I,J)=dFUGREP_dN(I,J)+dFUGATT_dN(I,J)+dFUGASSOC_dN_num(I,J)-dZ_dN(J)/Z

				IF (I.eq.J) THEN
					dFUG_dN(I,J)=(dFUGC1_dN(I,J)+dP_dN(J)/PMpa+(1.d0-xFrac(I))/xFrac(I))!/totMoles
				ELSE
					dFUG_dN(I,J)=(dFUGC1_dN(I,J)+dP_dN(J)/PMpa-1.d0)!/totMoles
				ENDIF
			ENDDO

		! First derivative in respect to N while P and T are constant
			DO J=1,NC
				dFUGC2_dN(I,J)=(dFUGC1_dN(I,J)-dFUGC1_dP(I)*dP_dN(J))!/totMoles     !called Fj in Michelsen's book, P.30
			ENDDO
		ENDDO
	endif !(isZiter==0)
	!if(LOUD)write(*,*)'x1,aAssoc',x(1),assoc
31	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! RHO (1) -VE IN FUGI')
33	FORMAT(1X,'LIQ=',1X,I1,2X,',','WARNING! RHO (2) -VE IN FUGI')
86	continue
	IF(IER(4)==1 .and. LOUD) WRITE(*,*)'ERROR IN ALPSOL'
	IER(1)=0
	if(IER(4).ne.0)IER(1)=1
	iErr=ier(1)
	initCall=0
686	FORMAT(' ERROR IN FuEsdWert.  ')
	RETURN
	END
          
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     SUBROUTINE DESIGNED BY JRE 6/02                            
!C                                                                       
!C     THIS SUBROUTINE SOLVES FOR THE MATRIX OF EQUATIONS WITH XD AND XA    
!C     IT USES Successive substitution for the off-diagonals and solves analytically for the XA(i) or XD(i)    
!C                                                                          
!C                                                                          
!C     XA,XD are initial guesses on input and answers on output
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJreAnal(ND,xFrac,XA,XD,nComps)
	USE GlobConst ! for LOUD
	IMPLICIT doublePrecision(A-H,K,O-Z)
	DIMENSION xFrac(NMX),XA(NMX),XD(NMX),ND(NMX)
	COMMON/ALPHA/ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	ierCode=0
	picard=1
	rmsErr=111 !avoid termination of doWhile before doing anything.
	iter=0
	do while(rmsErr.gt.1E-8.and.iter.lt.113)
		iter=iter+1
		rmsErr=0
		DO iComp=1,nComps
			SUMDON=0
			SUMACC=0
			DO J=1,nComps
				if(iComp.ne.J)then
					SUMDON=SUMDON+xFrac(J)*ND(J)*XD(J)*ALPHAD(iComp,J)
					SUMACC=SUMACC+xFrac(J)*ND(J)*XA(J)*ALPHDA(iComp,J)
				endif      
			enddo
			fracDA=xFrac(iComp)*Nd(iComp)*alphDA(iComp,iComp)
			fracAD=xFrac(iComp)*Nd(iComp)*alphAD(iComp,iComp)
			fa=XA(iComp)-1/( 1+sumdon+fracAD*XD(iComp) )
			fd=XD(iComp)-1/( 1+sumacc+fracDA*XA(iComp) )
			!solving the quadratic between XAi and XDi to get next guess
			quadB=(1+sumDon)*(1+sumAcc)+fracAD-fracDA
			sqArg=quadB*quadB+4*fracDA*(1+sumAcc)*(1+sumDon)
			if(sqArg.lt.0)then
				ierCode=1
				if(LOUD)pause 'Error in XsolJre: sqArg < 0'
			endif
			XA(iComp)=2*(1+sumAcc)/( quadB+SQRT(sqArg) )
			xaNew=1/( 1+sumDon+fracAD*XD(iComp) )
			!XA(iComp)=XA(iComp)*(1-picard)+picard*xaNew
			!XD(iComp)=XD(iComp)*(1-picard)+picard/(1+sumacc)
			XD(iComp)=1/( 1+sumAcc+fracDA*XA(iComp)	)  !using freshly computed XAi here makes this the analytical solution.
			rmsErr=rmsErr+fa*fa
			rmsErr=rmsErr+fd*fd
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		rmsOld=rmsErr
	enddo
	if(iter.gt.111)then
		ierCode=2
		if(LOUD)pause 'Error in XsolJre: iter>111'
	endif

	RETURN
	END      
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     SUBROUTINE DESIGNED BY JRE 6/02                            
!C                                                                       
!C     THIS SUBROUTINE SOLVES FOR THE MATRIX OF EQUATIONS WITH XD AND XA    
!C     IT USES Successive substitution    
!C                                                                          
!C                                                                          
!C     XA,XD are initial guesses on input and answers on output
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJreSimple(ND,xFrac,XA,XD,nComps)
	USE GlobConst ! for LOUD
	IMPLICIT doublePrecision(A-H,K,O-Z)
	!DIMENSION ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	!DIMENSION X(NMX),ND(NMX),VX(NMX)
	DIMENSION xFrac(NMX),XA(NMX),XD(NMX),ND(NMX)
	COMMON/ALPHA/ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	ierCode=0
	wegParm=0.25
	rmsErr=111 !avoid termination of do while before doing anything.
	iter=0
	do while(rmsErr.gt.1E-8.and.iter.lt.113)
		iter=iter+1
		rmsErr=0
		DO iComp=1,nComps
			SUMDON=0
			SUMACC=0
			DO J=1,nComps
				SUMDON=SUMDON+xFrac(J)*ND(J)*XD(J)*ALPHAD(iComp,J)      
				SUMACC=SUMACC+xFrac(J)*ND(J)*XA(J)*ALPHDA(iComp,J)      
			enddo
			fa=(XA(iComp)-1/(1+sumdon))
			fd=(XD(iComp)-1/(1+sumacc))
			XA(iComp)=XA(iComp)*wegParm+(1-wegParm)/(1+sumdon)
			XD(iComp)=XD(iComp)*wegParm+(1-wegParm)/(1+sumacc)
			rmsErr=rmsErr+fa*fa
			rmsErr=rmsErr+fd*fd
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
	enddo
	if(iter.gt.111)then
		ierCode=1
		if(LOUD)pause 'Error in XsolJre: iter>111'
	endif

	RETURN
	END      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     SUBROUTINE DESIGNED BY JRE 6/02                            
!C                                                                       
!C     THIS SUBROUTINE SOLVES FOR THE MATRIX OF EQUATIONS WITH XD AND XA    
!C     IT USES Successive substitution    
!C                                                                          
!C                                                                          
!C     XA,XD are initial guesses on input and answers on output
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJreWeg(ND,xFrac,XA,XD,nComps)
	USE GlobConst ! for LOUD
	IMPLICIT doublePrecision(A-H,K,O-Z)
	!DIMENSION ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	!DIMENSION X(NMX),ND(NMX),VX(NMX)
	DIMENSION xFrac(NMX),XA(NMX),XD(NMX),ND(NMX)
	DIMENSION fa(NMX),fd(NMX),faOld(NMX),fdOld(NMX)
	DIMENSION xaOld(NMX),xdOld(NMX)
	COMMON/ALPHA/ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	ierCode=0
	wegParm=0
	iWegFreq=4
	rmsErr=111 !avoid termination of do while before doing anything.
	iter=0
	do while(rmsErr.gt.1E-8.and.iter.lt.113)
		iter=iter+1
		rmsErr=0
		DO I=1,nComps
			SUMDON=0
			SUMACC=0
			DO J=1,nComps
				SUMDON=SUMDON+xFrac(J)*ND(J)*XD(J)*ALPHAD(I,J)      
				SUMACC=SUMACC+xFrac(J)*ND(J)*XA(J)*ALPHDA(I,J)      
			enddo
			fa(i)=(XA(I)-1/(1+sumdon))
			fd(i)=(XD(I)-1/(1+sumacc))
			rmsErr=rmsErr+fa(i)*fa(i)
			rmsErr=rmsErr+fd(i)*fd(i)
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		do iComp=1,nComps
			changeAcc=fa(iComp)/( fa(iComp)-faOld(iComp) ) &
				*( XA(iComp)-xaOld(iComp) )								!crude wegstein
			changeDon=fd(iComp)/( fd(iComp)-fdOld(iComp) ) &
				*( XD(iComp)-xdOld(iComp) )								!crude wegstein
			wegParm=0
			if(ABS(MOD(float(iter),float(iWegFreq))).lt.0.01)&			!check for hit on iWegFreq
				wegParm=changeAcc/( XA(iComp)-xaOld(iComp) )
			xaOld(iComp)=xa(iComp)
			xdOld(iComp)=xd(iComp)
			faOld(iComp)=fa(iComp)
			fdOld(iComp)=fd(iComp)
			!XA(iComp)=XA(iComp)*wegParm+(1-wegParm)*(XA(iComp)-fa(iComp))
			XA(iComp)=XA(iComp)-(1-wegParm)*fa(iComp) !same as above
			XD(iComp)=XD(iComp)-(1-wegParm)*fd(iComp)
			!XA(iComp)=XA(iComp)-changeAcc*1.3
			!XD(iComp)=XD(iComp)-changeDon*1.3
		enddo
		rmsOld=rmsErr

	enddo
	if(iter.gt.111)then
		ierCode=1
		if(LOUD)pause 'Error in XsolJre: iter>111'
	endif

	RETURN
	END      
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     SUBROUTINE DESIGNED BY JRE 6/02                            
!C                                                                       
!C     THIS SUBROUTINE SOLVES FOR THE MATRIX OF EQUATIONS WITH XD AND XA    
!C     IT USES Successive substitution    
!C                                                                          
!C                                                                          
!C     XA,XD are initial guesses on input and answers on output
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	SUBROUTINE XsolJreSec(ND,eta,xFrac,XA,XD,nComps)
	USE GlobConst ! for LOUD
	IMPLICIT doublePrecision(A-H,K,O-Z)
	!DIMENSION ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	!DIMENSION X(NMX),ND(NMX),VX(NMX)
	DIMENSION xFrac(NMX),XA(NMX),XD(NMX),ND(NMX)
	DIMENSION fa(NMX),fd(NMX),faOld(NMX),fdOld(NMX)
	DIMENSION xaOld(NMX),xdOld(NMX)
	COMMON/ALPHA/ALPHAD(NMX,NMX),ALPHDA(NMX,NMX)
	data etaOld/0.0/
	ierCode=0
	if(ABS(eta-etaOld).gt.0.05)then
		do iComp=1,nComps
			xaOld(iComp)=xa(iComp)*0.95
			xdOld(iComp)=xd(iComp)*0.95
			faOld(iComp)=1
			fdOld(iComp)=1
			rmsOld=1
		enddo
	endif
	wegParm=1
	rmsErr=111 !avoid termination of do while before doing anything.
	iter=0
	do while(rmsErr.gt.1E-8.and.iter.lt.113)
		iter=iter+1
		rmsErr=0
		DO I=1,nComps
			SUMDON=0
			SUMACC=0
			DO J=1,nComps
				SUMDON=SUMDON+xFrac(J)*ND(J)*XD(J)*ALPHAD(I,J)      
				SUMACC=SUMACC+xFrac(J)*ND(J)*XA(J)*ALPHDA(I,J)      
			enddo
			fa(i)=(XA(I)-1/(1+sumdon))
			fd(i)=(XD(I)-1/(1+sumacc))
			rmsErr=rmsErr+fa(i)*fa(i)
			rmsErr=rmsErr+fd(i)*fd(i)
		enddo
		rmsErr=sqrt( rmsErr/(2*nComps) )
		do iComp=1,nComps
			changeAcc=fa(iComp)/( fa(iComp)-faOld(iComp) ) &
				*( XA(iComp)-xaOld(iComp) )								!crude wegstein
			changeDon=fd(iComp)/( fd(iComp)-fdOld(iComp) ) &
				*( XD(iComp)-xdOld(iComp) )								!crude wegstein
			xaOld(iComp)=xa(iComp)
			xdOld(iComp)=xd(iComp)
			faOld(iComp)=fa(iComp)
			fdOld(iComp)=fd(iComp)
			XA(iComp)=XA(iComp)*(1-wegParm)+wegParm/(1+sumdon)
			XD(iComp)=XD(iComp)*(1-wegParm)+wegParm/(1+sumacc)
			XA(iComp)=XA(iComp)-changeAcc*1.3
			XD(iComp)=XD(iComp)-changeDon*1.3
		enddo
		rmsOld=rmsErr

	enddo
	if(iter.gt.111)then
		ierCode=1
		if(LOUD)pause 'Error in XsolJre: iter>111'
	endif
	etaOld=eta

	RETURN
	END      
      
	SUBROUTINE GoldenZ2(X,nC,KVE,ND,NAS,NDS,bMix,CbMix,YqbMix,k1YbMix, &
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

	do 50 iter=1,111

	if(err2 < err1)then
		rhoLo = rho1;
		rho1  = rho2;
		rho2  = rgold*rho1+cgold*rhoHi;
		errLo = err1;
		err1  = err2;

		CALL AlpSolEz(X,nC,KVE,ND,bMix,NAS,NDS,rho2,zAssoc,XA,RALPH,fAssoc,ier)
		zRep = 4 * cbMix * rho2 / (1 - 1.90*bMix*rho2);
		zAtt = -zM * YqbMix * rho2 / (1 + k1YbMix * rho2);
		zFactor = (1 + zRep + zAtt + zAssoc);
		err2 = ABS(PoRT-rho2*zFactor);

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

50	continue

	!c//that's the best we can do for now folks. get props at best guess and return.
	rho = (rho1+rho2) * 0.5;
	zFactor = PoRT/rho;
	eta = bMix*rho;
	return
	end

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine xsMixRule(dXsYbN,xsYb,dXsYbB,bVol,moFrac,tK,nComps)
	!c  compute the excess mixture value for the bY parameter in denom of zAtt
	!c  use NRTL style of expression
	USE GlobConst !for NMX
	USE BIPs
	implicit doubleprecision(A-H,K,O-Z)
	double precision moFrac(nComps)
	integer jComp,kComp
	dimension bVol(nComps),dXsYbN(nComps)
	DIMENSION omega(NMX,NMX),sumDenom(NMX),sumNumer(NMX)

	bVolMix=0.0
	do jComp=1,nComps
		bVolMix=bVolMix+moFrac(jComp)*bVol(jComp)
	enddo
	!  begin the computation of the xsY.   use NRTL here.
	!  this xsY will be multiplied by bVolMix to get xsYb at the bottom.
	!  ref.  ReID et al. 1987, PGL p 274,256
	!  bigTauJI = xsTauJI/T + xsTauTJI/T^2
	!  sumNumer is the sumJ(bigTauJI*OmegaJI*xK), 
	!  sumDenom is the sumJ(OmegaJI*xK), 
	!  OmegaKI = Gki = EXP(-alphaKI*bigTauKI)
	xsY =0.0
	dXsYB =0.0
	do iComp=1,nComps
		sumDenom(iComp)=0.0
		sumNumer(iComp)=0.0
		do jComp=1,nComps
			bigTauJI=xsTau(jComp,iComp)/tK+xsTauT(jComp,iComp)/(tK*tK)
			omega(jComp,iComp)=EXP( -xsAlpha(jComp,iComp)*bigTauJI )
			sumDenom(iComp)=sumDenom(iComp)+moFrac(jComp)*omega(jComp,iComp)
			sumNumer(iComp)=sumNumer(iComp)+moFrac(jComp)*omega(jComp,iComp)*bigTauJI
		enddo
		!  compute the excess GE
		xsY=xsY+moFrac(iComp)*sumNumer(iComp)/sumDenom(iComp)
		!  compute the excess enthalpy.  this formula is NOT given by reID et al.
		!  d(...)B is beta*d(...)/dbeta.  d2(...)BB is beta^2*d2(...)/dbeta^2, eTC...
		dSumNumB=0.d0
		dSumDenB=0.d0
		do jComp=1,nComps
			bigTauJI=xsTau(jComp,iComp)/tK+xsTauT(jComp,iComp)/(tK*tK)
			dBigTauB=xsTau(jComp,iComp)/tK+2*xsTauT(jComp,iComp)/(tK*tK)
			dOmegaB=omega(jComp,iComp)*(-xsAlpha(jComp,iComp)*dBigTauB)
			dSumNumB=dSumNumB+moFrac(jComp)*( omega(jComp,iComp)*dBigTauB+bigTauJI*dOmegaB )
			dSumDenB=dSumDenB+moFrac(jComp)*dOmegaB
		enddo
		dXsYB=dXsYB+moFrac(iComp)/sumDenom(iComp)*( dSumNumB-dSumDenB*sumNumer(iComp)/sumDenom(iComp) )
	enddo
	!  get the partials and check the xsY by summing the partials
	xsYck=0
	do kComp=1,nComps
		bigSumK=0.d0
		do jComp=1,nComps
			bigTauKJ=xsTau(kComp,jComp)/tK+xsTauT(kComp,jComp)/(tK*tK)
			bigSumK=bigSumK+moFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)&
			*( bigTauKJ-sumNumer(jComp)/sumDenom(jComp) )
		enddo
		xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
		!  temporary misnomer here. the "b" part is factored in below.
		!  this is really d(nYE)/dnk
		dXsYbN(kComp)=xsLnGamma
		xsYck=xsYck+moFrac(kComp)*dXsYbN(kComp)
    enddo
    if(LOUD)then
	    if(ABS(xsYck-xsY).gt.1.d-4)pause 'partials dont add up in xsMix..'
    end if

	!till now, everything was operating on xsY.  what we want is xsYb.
	xsYb=xsY*bVolMix
	dXsYbB=dXsYB*bVolMix
	do kComp=1,nComps
		dXsYbN(kComp)=dXsYbN(kComp)*bVolMix+bVol(kComp)*xsY
	ENDDO
	return
	end
