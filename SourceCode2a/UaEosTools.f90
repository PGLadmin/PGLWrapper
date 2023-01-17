!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE PsatEar(tK,pMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
	!COMPUTE VAPOR PRESSURE GIVEN tK AND INITIAL GUESS (IE.pMPa)
	!Ref: Eubank et al, IECR 31:942 (1992)
	!AU:  JRE Aug2010
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,..., etaPass, ...
	!USE Assoc ! temporary for checking assoc parms
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!character tabChar*1
	DIMENSION fugc(NMX),xFrac(NMX) 
	DIMENSION IER(12)
	LOGICAL LOUDER
	CHARACTER*77 errMsg(0:11) !,errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data initCall/1/
	errMsg(1)='No spinodal max/min'
	errMsg(2)='Psat iteration did not converge'
	errMsg(3)='Liquid FUGI call failed on last iteration'
	errMsg(4)='Vapor  FUGI call failed on last iteration'
	errMsg(5)='zVap=zLiq on last iteration'
    errMsg(6)='PsatEar: FUGI returned T < Tmin error'
    errMsg(7)='PsatEar: Calculated Psat < 0.0001 MPa error'
	LOUDER=LOUD
	LOUDER=.TRUE.
	LOUDER=.FALSE.
	if(LOUDER.and.initCall)write(dumpUnit,*)'EAR method'
	tMin=Tc(1)*0.25
	xFrac(1)=1
	NC=1
	ierCode=0 
	zCritEff=ZC(1)
	if(iEosOpt==1 .or. iEosOpt==11)zCritEff=0.3074
	if(iEosOpt==2 .or. iEosOpt==4)zCritEff=0.34  !ESD
	rhoCrit=Pc(1)/(zCritEff*Rgas*Tc(1))
	rhoVap=rhoCrit*1.000 !-ve value on input means use default value for etaHi
    if(LOUDER)then 
        if(DEBUG) then
	        OPEN(67,file='c:\spead\calceos\output\isotherm.txt') !bonus: write the isotherm from iterations on spinodal
        else
	        OPEN(67,file='spinodal.txt')
        endif
	    write(67,*)' T(K)    rho(mol/cc)   P(MPA)        Z     DA_NKT    DU_NKT    eta'
    endif
	Tr=tK/Tc(1)
	if(Tr > 0.85)then
		if(LOUDER)write(dumpUnit,*)'Psat: calling spinodal for vapor.'
		call Spinodal(NC,xFrac,tK,0,rhoVap,zVap,aDepVap,dU_NkT,iErrVap)
		if(iErrVap==3)then !first definite call to fuVtot suffices to check for T < Tmin.
			ierCode=6
			goto 86
		endif
		rhoLiq=rhoVap !-ve value on input means use default value for etaLo
		if(LOUDER)write(dumpUnit,*)'Psat: calling spinodal for Liquid.'
		call Spinodal(NC,xFrac,tK,1,rhoLiq,zLiq,aDepLiq,dU_NkT,iErrLiq)
		if(LOUDER)close(67)
		pMax=Rgas*tK*(rhoVap*zVap) 
		pMin=Rgas*tK*(rhoLiq*zLiq)
		if(LOUDER.and.pMax < 1e-11)write(dumpUnit,*)'Psat: 0~pMax=',pMax
		relDiffP=ABS(pMax-pMin)/pMax
		if(iErrVap.ne.0 .or. iErrLiq.ne.0 .or. relDiffP < 1D-3)then
			ierCode=1
			if(LOUDER)write(dumpUnit,'(a,f8.1,a)')' PsatEar: Spinodal error at T = ',tK,' T > Tc?'
			!if(LOUD)write(dumpUnit,*) 
			goto 86
		endif
		!if(pMin < 0)pMin=0
		etaLiq=rhoLiq*bVolCC_mol(1)
		etaVap=rhoVap*bVolCC_mol(1)
		pOld=(pMin+pMax)/2 !when approaching critical, this works well.
		if(pOld < zerotol)pOld=Pc(1)/2 ! OK for Tr > 0.85(?)
		if(LOUDER)write(dumpUnit,'(a,4f11.4)')' PsatEar: after spinodal: pMin,pMax,pInit= ',pMin,pMax,pOld
		if(LOUDER)write(dumpUnit,'(a,6E12.4)')' PsatEar: spinodal=>etaLiq,etaVap=',etaLiq,etaVap
		if( (rhoVap/rhoLiq) < 0 .and. LOUDER)then
			ier(1)=18
			write(dumpUnit,*) 'PsatEar: error after spinodals: rhoVap/rhoLiq < 0'
			goto 86
		endif
		pEubank=Rgas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995), assuming 1.2x shifts on V&L relative to spinodals.
		!pOld=pEubank/10  !JRE: I find Eubank overestimates P sometimes. It even gets higher than pMax
	else ! If Tr < 0.85, compute pSatItic estimate  FPE,501:112236(19), JCP,110:3043(99)
		pOld=1D-1 !If p < 0, then the vapor density goes negative and log(rhoLiq/rhoVap) is indeterminate.
        !NOTE: Don't use pOld=0 when pOld > 0. You might get an error because vdw loop doesn't cross zero!
		if(LOUDER)write(dumpUnit,*)'PsatEar:calling liquid fugi for ITIC initial guess. ID,Tc=',ID(1),Tc(1)
		!if(LOUDER)write(dumpUnit,601)'PsatEar:eAcc,eDon',(eAcceptorKcal_mol(1,j),j=1,nTypes(1)),(eDonorKcal_mol(1,j),j=1,nTypes(1))
		!call FUGI(tK,pOld,xFrac,NC,1,FUGC,zLiq,ier) !LIQ=1=>liquid root with fugc calculation so we get Ares.
		CALL FugiTP( tK,Pold,xFrac,NC,1,rhoLiq,zLiq,aResLiq,FUGC,uSatL,iErrF )
		if(iErrF==2)then
			if(LOUDER)write(dumpUnit,*)'PsatEar:T < Tmin? Sorry.'
			ierCode=16
			return
		elseif(zLiq < zeroTol)then
			if(LOUDER)write(dumpUnit,601)' PsatEar: Initial call to fugi. 0 > zLiq=',zLiq
			ierCode=17
			goto 86
		elseif(iErrF > 0)then
			if(LOUDER)write(dumpUnit,'(a,i12)')' PsatEar: Error from ITIC call to fugi. ier=',(ier(i),i=1,7)
			ierCode=19
			return
		endif
		!AresLiq=Ares_RT						!AresVap  +  Zvap-1 -ln(Zvap) =AresLiq+ZLiq-1-ln(ZLiq)							 
		!rhoLiq=pOld/(zLiq*Rgas*tK)		!=>	!B2*rhoVap+B2*rhoVap+ln(rhoVap/rhoLiq) =AresLiq+0-1  
		!rhoLiq=etaPass/bVolCC_mol(1) ! crazy precision is required for acids (ZL < 1D-7) or when etaLiq > 0.98 (and Z->0). e.g. pentanoic acid (1258)
		etaLiq=rhoLiq*bVolCC_mol(1)
		if(etaLiq < 0.15D0 .and. LOUDER)write(dumpUnit,*)'PsatEar: etaLiq < 0.15? pOld,etaLiq,tK=',pOld,rhoLiq*bVolCC_mol(1),tK
        rhoVap=rhoLiq*EXP( AresLiq-1 ) !pSatItic, FPE, 501:112236(19), Eq 19. Ares_RT from GlobConst
		pSatItic=rhoVap*Rgas*tK
		etaVap=rhoVap*bVolCC_mol(1)
		!call FuVtot(1,tK,pSatItic,xFrac,NC,0,FUGC,zVap,ier)
		if(LOUDER)write(dumpUnit,601)' PsatEar: Calling init FuVtot. T,etaLiq,etaVap,aResLiq=',tK,rhoLiq*bVolCC_mol(1),rhoVap*bVolCC_mol(1),aResLiq    
		call FuVtot(1,tK,1/rhoVap,xFrac,NC,FUGC,zVap,aRes,uRes,iErrF)		!isZiter=1=>vapor Z with no fugc calculation. zVap far from zero.
		if(iErrF > 9)then
			if(LOUDER)write(dumpUnit,'(a,i4)')' PsatEar: Initial Vapor FuVtot Error=',iErrF
			ierCode=19
			return
		endif ! 
		!B2cc_mol*rhoVap=(zVap-1)
		rhoVap=rhoLiq*exp( aResLiq-1-2*(zVap-1) )! => rhoVap=rhoLiq*EXP(AresLiq-1-2*B2*rhoVap).
		pSatItic=rhoVap*Rgas*tK*zVap
		if(pSatItic > zeroTol)then
			pTest=pSatItic
			zLiq=pSatItic/(rhoLiq*Rgas*tK) ! if Psat > 0 then zLiq > 0
			FugcLiq=AresLiq+zLiq-1-LOG(zLiq)
			if(LOUDER)write(dumpUnit,601)' PsatEar: Init-zVap,zLiq,Psat,FugcVap,FugcLiq=',zVap,zLiq,pSatItic,Fugc(1),FugcLiq
		else
			if(LOUDER)write(dumpUnit,601)' PsatEar: 0 > pSatItic = ?',pSatItic
			pTest=zeroTol*10
		endif
		if(LOUDER)write(dumpUnit,'(a,2f8.4,e11.4)')' PsatEar: pSatItic,eta,rhoLiqG/cc',pTest,rhoLiq*bVolCC_mol(1),rhoLiq*rMw(1)
	    !pEuba2=Rgas*tK/(1/rhoVap-1/rhoLiq)*( Ares_RT- 0 +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
        !NOTE: When 1/rhoVap >> 1/rhoLiq, pEuba2=Rgas*tK*rhoVap*( Ares_RT-ln(rhoVap/rhoLiq) )=pTest*( Ares_RT - (Ares_RT-1) )
        pOld=pTest
	endif ! Tr > 0.85
601 format(1x,a,6E12.4)
	etaLiq=rhoLiq*bVolCC_mol(1)
	if(pOld < zeroTol)then
		pOld=zeroTol
		if(LOUDER)write(dumpUnit,*)' PsatEar: zeroTol > pInit = ?',pOld
	endif
	!if(LOUD .and. ABS(etaL-etaLiq)/etaLiq > 0.0001 .and. LOUD)write(dumpUnit,*) 'Psat: warning lost precision in rhoLiq'
	pBest=1234
	fBest=1234
	if(LOUDER.and.initCall)write(dumpUnit,*)'PsatEar:Initial pSat,etaL=',pOld,etaLiq
	if(LOUDER.and.initCall)write(dumpUnit,*) 'Check initial guess for Psat.'
	itMax=33
	do iter=1,itMax !iterate on pOld according to Eubank criterion
		!write(dumpUnit,*)'calling fugi for liquid iteration.'
		!call FUGI(tK,pOld,xFrac,NC,1,FUGC,zLiq,ier) !LIQ=3=>liquid root with no fugc calculation
		CALL FugiTP( tK,Pold,xFrac,NC,1,rhoLiq,zLiq,aResLiq,FUGC,uSatL,iErrF )
		if(iErrF > 0)ierCode=3 !declare error but don't stop. if future iterations give valid fugi(), then no worries
		!chemPotL=fugc(1)
		!aDepLiq=aRes_RT
		!rhoLiq=pOld/(zLiq*Rgas*tK)
		!rhoLiq=etaPass/bVolCC_mol(1) ! crazy precision is required when etaLiq > 0.98 (and Z->0). e.g. pentanoic acid (1258)
		!write(dumpUnit,*)'calling fugi for liquid iteration.'
		!call Fugi(tK,pOld,xFrac,NC,0,FUGC,zVap,ier) !LIQ=2=>vapor root with no fugc calculation
		CALL FugiTP( tK,Pold,xFrac,NC,0,rhoVap,zVap,aResVap,FUGC,uSatV,iErrF )
		if(iErrF > 0)ierCode=4 !declare warning error but don't stop. if future iterations give valid fugi(), then no worries.
		!chemPotV=fugc(1)
		!aDepVap=aRes_RT
		!rhoVap=etaPass/bVolCC_mol(1)
		if( (rhoVap/rhoLiq) < 0 )then
			write(dumpUnit,601)' PsatEar: rho < 0? rhoVap,Liq=',rhoVap,rhoLiq
			If(LOUDER)write(dumpUnit,*) 'Psat error : rhoVap/rhoLiq < 0'
			ierCode=19
			goto 86
		endif
		if(ABS(zVap-zLiq) < 1e-3)then
			if(LOUDER)write(dumpUnit,*)'rho,zLiq,zVap',rhoLiq,zVap,zLiq
			if(LOUDER)write(dumpUnit,*) 'PsatEar: zVap=zLiq'
			ierCode=15
			goto 86
		endif
        if(pOld < 1E-8 .and. iter > 2)then
			etaLiq=rhoLiq*bVolCC_mol(1)
			if(etaLiq < 0.15D0 .and. LOUDER)write(dumpUnit,601)'Psat: etaLiq < 0.15? pOld,zLiq,tK,etaLiq=',pOld,zLiq,tK,etaLiq
            rhoVap=rhoLiq*EXP( aResLiq+zLiq-1-aResVap+zVap-1 ) !pSatItic, FPE, 501:112236(19), Eq 19.  aDepVap might be significant for carbo acids.
            pTest=rhoVap*Rgas*tK*zVap
			ierCode=1 ! assign warning in this case. 
			exit
        endif
		!rLogRhoRat=DLOG(rhoVap/rhoLiq)	!for debugging
		pTest=pOld*( aResLiq-aResVap-DLOG(rhoVap/rhoLiq) )/(zVap-zLiq)	!EAR method of Eubank and Hall (1995)
		!fErr= aDepLiq+zLiq-1 -(aDepVap+zVap-1) -DLOG(zLiq/zVap)
		fErr= aResLiq+zLiq-1 -(aResVap+zVap-1) -DLOG(rhoVap/rhoLiq) 
		if( ABS( fErr) < fBest)then	! for very low pSat, like propane.
			fBest=ABS( fErr)
			pBest=pTest
		endif
		change=pTest-pOld
		if(LOUDER)write(dumpUnit,'(a,i3,8E12.4)')' PsatEar:iter,pOld,pTest',iter,pOld,pTest
		if(LOUDER)write(dumpUnit,601)'PsatEar:aResLiq,aResVap,zLiq,zVap',aResLiq,aResVap,zLiq,zVap
		tol=1D-6
		if(pTest > 0)then
			pOld=pTest
		else
			pOld=pOld/10
        endif
		!pOld =Rgas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
		if(ABS(change/pOld) < tol)EXIT   
	enddo
	pMPa=pOld
	if(ABS(change/pMPa) > tol .or. iter.ge.itMax)then
		ierCode=2
		if(LOUDER)write(dumpUnit,'(a,i3,8E12.4)')' PsatEar:tolerance not met. iter,T,etaL,etaV,relDev',itMax,tK,rhoVap*bVolCc_mol(1),rhoLiq*bVolCc_mol(1),change/pMPa
		pMPa=pBest
		goto 86
	endif
	!iter=0                                  
	!check that Fugi did not give warning on last call
86	continue
	if(LOUDER.and.initCall)write(dumpUnit,*) 'Psat: P,rhoVap,rhoLiq',pMPa,rhoVap,rhoLiq
	if(LOUDER)write(dumpUnit,'(a,3f9.5)') ' Psat: P,etaVap,etaLiq',pMPa,rhoVap*bVolCc_mol(1),rhoLiq*bVolCc_mol(1)
	isZiter=0
	!if(ierCode==0) call FuVtot(isZiter,tK,1/rhoLiq,xFrac,NC,FUGC,ZLiq,iErrFu) !one last call to fugi for chemPo and debugging. We know the density so use it.
	fugc(1)=86.86D0
	if( zLiq > 0)then
		fugc(1) = aResLiq+zLiq-1 -DLOG(zLiq)
	else
		if(LOUDER)write(dumpUnit,*) 'Psat: error zLiq < 0 after converging Psat.'
	endif
	if(LOUDER.and.initCall)write(dumpUnit,*)'Psat: cmprsblty,CvRes_R=',cmprsblty,CvRes_R
	!if(ierCode==0)call Fugi(tK,pMPa,xFrac,NC,1,FUGC,zLiq,ier) !one last call to fugi for chemPo and debugging.
	chemPot=FUGC(1) !since chemPotL=chemPotV at pSat       
	uSatL=DUONKT
	if( (isTPT .or. isESD) .and. pMPa < 1.D-4 .and. LOUDER)write(dumpUnit,*)'Psat: 1.D-4 > pMpa = ',pMPa
	if( (isTPT .or. isESD) .and. pMPa < 1.D-4)ierCode=7
	if(iEosOpt==22 .and. pMPa < 1.D-2)ierCode=7
	!Convert densities to g/cc.
    if(LOUDER)write(dumpUnit,'(f10.2,3f14.10,i2)')tK,pMPa,rhoVap*rMw(1),rhoLiq*rMw(1),ierCode
	if(LOUDER)write(dumpUnit,*) 'check final Psat'
	initCall=0
	RETURN
	END	! subroutine PsatEar()

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
	SUBROUTINE Spinodal(NC,xFrac,tKelvin,LIQ,rho,zFactor,dA_NkT,dU_NkT,iErr) 
	!Compute the points where dP/dRho=0 appear at given T.
	!Method: Golden Section search
	!Ref: NumRep 
	!Ref: Eubank and Hall, AICHEJ, 41:924 (1995)
	!AU: JRE, Aug 2010
	!Notes:
	!  Assume gmol(i)=xFrac(i) => vTotCc = 1/rho
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION xFrac(*)
	dimension FUGC(NMX) !,vMolec(NMX)
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data rgold,cgold,initCall/0.61803399,0.38196602,1/
	!open(67,file='spinodal.txt') this is opened in the calling routine. that way we can compile V and L together.
	isZiter=2 !derivative dP_drho not necessary (although this would be easier if we used them)
	iErr=0
	!if(eta==etaHi)iErr=1 !this means there is no max from vap search b/c T>Tc.
	!if(eta==etaLo)iErr=2 !this means there is no min from liq search b/c T>Tc.
    !iErr=3 means T < Tmin for at least one component
	bMix=0.d0
	do i=1,NC
		bMix=bMix+xFrac(i)*bVolCC_mol(i)
	enddo
	if(bMix < 1E-11 .and. LOUD)write(dumpUnit,*)'Spinodal: nonsense bMix=',bMix
	if(Tc(1) < 1E-11 .and. LOUD)write(dumpUnit,*)'Spinodal: nonsense Tc(1)=',Tc(1)

	Tr=tKelvin/Tc(1)
	!if(isTPT)write(dumpUnit,*)'Spinodal: etaMax=',etaMax
	if(rho<0)then !this is the signal to use default estimates, but really using rhoCrit works much better in the critical region.
		etaLo = 1e-11
		etaHi = 0.17
		if(LIQ==1)then
			etaLo=0.1;
			etaHi=etaMax;	!must be less than 0.526 for ESD.
			!if(isESD)etaHi=0.52
		endif
		if(iEosOpt==1 .or. iEosOpt==11)then	 ! PR flavors
			etaLo=etaLo*2
			!etaHi=etaHi*2
			!if(etaHi > 0.99)etaHi=etaMax
		endif
	else
		etaCrit=rho*bMix
		etaLo=etaCrit*0.0001
		etaHi=etaCrit !if LIQ=0 and rho>0, then use the given value as the guess for rhoVapHi
		if(LIQ==1)then
			etaLo=etaCrit
			etaHi=etaMax
			!etaHi=0.77
			!if(isESD)etaHi=0.52
			!if(iEosOpt==1.or.iEosOpt==11)etaHi=0.99
		endif
	endif
	minimax= -1	!routine was written to minimize, therefore take -ve to maximize.
	if(LIQ)minimax=1
	if(LOUD.and.initCall)write(dumpUnit,*)'Spinodal: LIQ,etaLo,etaHi',LIQ,etaLo,etaHi
	rhoLo = etaLo/bMix;
	rhoHi = etaHi/bMix;
	rho1=rhoLo+cgold*(rhoHi-rhoLo);
	rho2=rhoLo+rgold*(rhoHi-rhoLo);

    if(LOUD .and. rhoLo < 1E-11)write(dumpUnit,*)'Spinodal: nonsense rhoLo =',rhoLo    
	call FuVtot(isZiter,tKelvin,1/rhoLo,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrFu)
	if(iErrFu)then
        if(LOUD)write(dumpUnit,*) 'Spinodal: error in initial calls to FuVtot'
        iErr=3
        goto 86
    endif
	pLo = minimax*rhoLo*zFactor*Rgas*tKelvin
	if(LOUD.and.initCall)write(dumpUnit,601)tKelvin,1/rhoLo,minimax*pLo,zFactor,DUONKT,DAONKT,rhoLo*bMix
	if(LOUD)write(67,601)tKelvin,1/rhoLo,minimax*pLo,zFactor,DUONKT,DAONKT,rhoLo*bMix
	call FuVtot(isZiter,tKelvin,1/rho1,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrFu)
	p1 = minimax*rho1*zFactor*Rgas*tKelvin
	if(LOUD)write(67,601)tKelvin,1/rho1,minimax*p1,zFactor,DUONKT,DAONKT,rho1*bMix
	call FuVtot(isZiter,tKelvin,1/rho2,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrFu)
	p2 = minimax*rho2*zFactor*Rgas*tKelvin
	if(LOUD)write(67,601)tKelvin,1/rho2,minimax*p2,zFactor,DUONKT,DAONKT,rho2*bMix
	call FuVtot(isZiter,tKelvin,1/rhoHi,xFrac,NC,FUGC,zFactor,aRes,uRes,iErrFu)
	pHi = minimax*rhoHi*zFactor*Rgas*tKelvin
	if(LOUD)write(67,601)tKelvin,1/rhoHi,minimax*pHi,zFactor,DUONKT,DAONKT,rhoHi*bMix
	!NOTE: iErrFu=1 for liquid because -ve Z values will cause error in fugacity calculation.
	itMax = 66-33*(14/rMw(1)) !Increase the precision for long chains cuz...
    if(itMax < 33)itMax=33  !H2 and He give problems for above formula
	etaOld=etaHi
101	continue
	do iter=1,itMax !1/2^33=1E-10; (2/3)^33=1E-6
		if(p2 < p1)then
			rhoLo = rho1;
			rho1  = rho2;
			rho2  = rgold*rho1+cgold*rhoHi;
			pLo = p1;
			p1  = p2;
			if(rho2 < 1D-33 .and. LOUD)write(dumpUnit,*)'Spinodal: 0~rho2=',rho2
			call FuVtot(isZiter,tKelvin,1/rho2,xFrac,NC,FUGC,zFactor,aRes,uRes,iErr)
			p2 = minimax*rho2*zFactor*Rgas*tKelvin
			if(LOUD)write(67,601)tKelvin,1/rho2,minimax*p2,zFactor,DUONKT,DAONKT,rho2*bMix
		else
			rhoHi = rho2;
			rho2  = rho1;
			rho1  = rgold*rho2+cgold*rhoLo;
			pHi = p2;
			p2  = p1;
			if(rho2 < 1D-33 .and. LOUD)write(dumpUnit,*)'Spinodal: 0~rho1=',rho1
			call FuVtot(isZiter,tKelvin,1/rho1,xFrac,NC,FUGC,zFactor,aRes,uRes,iErr)
			p1 = minimax*rho1*zFactor*Rgas*tKelvin
			if(LOUD)write(67,601)tKelvin,1/rho1,minimax*p1,zFactor,DUONKT,DAONKT,rho1*bMix
		endif
		eta = bMix*(rho1+rho2)/2
		if(ABS(eta-etaOld) < 1D-6)exit !close enough
		etaOld=eta
		if(LOUD.and.initCall)write(dumpUnit,'(a,i3,1x,2E12.4,i3)')' Spinodal: iter,eta,P,iErr',iter,eta,minimax*(p1+p2)/2,iErr
		!write(dumpUnit,*) 
	enddo
	!c//that's the best we can do for now folks. get props at best guess and return.
	rho = (rho1+rho2) /2
	zFactor = minimax*(p1+p2)/(2*rho*Rgas*tKelvin)
	if(eta==etaHi)iErr=1 !this means there is no max from vap search b/c T>Tc.
	if(eta==etaLo)iErr=2 !this means there is no min from liq search b/c T>Tc.
	dA_NkT=DAONKT
	dU_NkT=DUONKT
	aRes_RT=dA_NKT
	uRes_RT=dU_NKT
601 FORMAT(f7.2,1x,2E12.4,1x,f8.4,1x,3F10.4)
86  continue
	initCall=0
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!C	PURPOSE: COMPUTE THE COMPRESSIBILITY,rhoD2PdRho2_RT, CvRes,& CpRes BY NUMERICAL DIFFERENTIATION
	!C	Input:
	!C  xFrac:			mole fraction vector
	!C	tKelvin:		T(K)
	!C	rhoMol_cc:		density (mol/cm3)
	!C	zFactor:		compressibility factor (used for central part of 2nd derivative)
	!C	Output:			Values are passed through GlobConst  
	!C  cmprsblty:		(1/RT)(dP/drho)_T dimensionless measure of compressibilty
	!C	rhoD2PdRho2_RT: (rho/RT)(d2P/drho^2)_T, dimensionless derivative of cmprsblty
	!C	CvRes_R:		(Cv-Cv^ig)/R, dimensionless measure of residual isochoric heat capacity.   
	!C	CpRes_R:		(Cp-Cp^ig)/R, dimensionless measure of residual isobaric heat capacity.   
	!C
	!C	CODED BY:   JRE, 4/25/2020
	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	SUBROUTINE NUMDERVS(NC,xFrac,tKelvin,rhoMol_cc,zFactor,iErr)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
	dimension xFrac(NC),FUGC(NC) ! ,bVol(NC)
	data initCall/1/
	iErr=0
	step=1.D-4
	bMix=SumProduct(NC,xFrac,bVolCC_mol)
	eta=bMix*rhoMol_Cc
	if(rhoMol_cc < zeroTol .or. tKelvin < zeroTol .or. eta < zeroTol .or. zFactor < zeroTol .or. eta > etaMax)then
		if(LOUD)write(dumpUnit,*)'NUMDERVS: nonsense input. NC,T,rho,Z,eta,b,x1,bMix=',NC,tKelvin, rhoMol_cc,zFactor,eta, bVolCC_mol(1),xFrac(1),bMix
		iErr=11
		return 
	endif
	isZiter=1 !we don't need the fugacity.	
	rhoPlus =rhoMol_cc*(1+step)
	rhoMinus=rhoMol_cc*(1-step)
	CALL FuVtot(isZiter,tKelvin,1/rhoPlus,xFrac,NC,FUGC,zPlus,aRes,uRes,iErrFu)
	if( iErrFu.ne.0 )then
		if(iEosOpt .ne. 5 .or. iErrFu >10)then !ignore warnings for speadmd
			if(LOUD)write(dumpUnit,*)'NUMDERVS: error from FuVtot on first call.'
			iErr=12
			return
		endif
	endif
	CALL FuVtot(isZiter,tKelvin,1/rhoMinus,xFrac,NC,FUGC,zMinus,aRes,uRes,iErrFu)
	rhoDZ_dRho=rhoMol_cc*(zPlus-zMinus)/(rhoPlus-rhoMinus)
	cmprsblty=zFactor+rhoDZ_dRho ! = (dP/dRho)/RT
	if(initCall.and.Loud)then
		!write(dumpUnit,*)'NUMDERVS: Z,rhoDZ_dRho',zFactor,rhoDZ_dRho
		!write(dumpUnit,*)'NUMDERVS: zPlus,zMinus',zPlus,zMinus
		!write(dumpUnit,*)'NUMDERVS: rho+,rho-,cmprsblty',rhoPlus,rhoMinus,cmprsblty
		write(dumpUnit,*)'NUMDERVS: cmprsblty',cmprsblty
	endif
	rho2DZ_dRho2=(zPlus+zMinus-2*zFactor)/(rhoMol_cc*step)**2
	rhoD2PdRho2_RT=2*rhoDZ_dRho+rho2DZ_dRho2
	TPlus =tKelvin*(1+step)
	TMinus=tKelvin*(1-step)
	isZiter=0 ! we need uRes for these
	CALL FuVtot(isZiter,TPlus,1/rhoMol_cc,xFrac,NC,FUGC,zPlus,aRes,uRes,iErrFu)
	uResPlus=uRes*TPlus	!uRes_RT passed through GlobConst
	if(initCall.and.Loud)write(dumpUnit,*)'NUMDERVS: 1-step,Tminus',1-step,Tminus
	CALL FuVtot(isZiter,TMinus,1/rhoMol_cc,xFrac,NC,FUGC,zMinus,aRes,uRes,iErrFu)
	uResMinus=uRes*TMinus
	TdZ_dT=tKelvin*(zPlus-zMinus)/( TPlus-TMinus )
	CvRes_R= (uResPlus-uResMinus)/( TPlus-TMinus )
	if( ABS(cmprsblty) < zeroTol)then
		if(LOUD)write(dumpUnit,*)'NUMDERVS: 0~cmprsblty=',cmprsblty
		CpRes_R=86.8686D0
		iErr=13
		return
	else
		CpRes_R=CvRes_R - 1 + (zFactor+TdZ_dT)**2/cmprsblty
	endif
	if(initCall.and.Loud)then
		!write(dumpUnit,*)'NUMDERVS: Z,rhoDZ_dRho',zFactor,rhoDZ_dRho
		!write(dumpUnit,*)'NUMDERVS: zPlus,zMinus',zPlus,zMinus
		!write(dumpUnit,*)'NUMDERVS: rho+,rho-,cmprsblty',rhoPlus,rhoMinus,cmprsblty
		write(dumpUnit,*)'NUMDERVS: CvRes_R,CpRes_R',CvRes_R,CpRes_R 
	endif
	initCall=0
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!C	PURPOSE: COMPUTE THE COMPRESSIBILITY,rhoD2PdRho2_RT, CvRes,& CpRes BY NUMERICAL DIFFERENTIATION
	!C	Input:
	!C  xFrac:			mole fraction vector
	!C	tKelvin:		T(K)
	!C	rhoMol_cc:		density (mol/cm3)
	!C	zFactor:		compressibility factor (used for central part of 2nd derivative)
	!C	Output:			Values are passed through GlobConst  
	!C  cmprsblty:		(1/RT)(dP/drho)_T dimensionless measure of compressibilty
	!C	rhoD2PdRho2_RT: (rho/RT)(d2P/drho^2)_T, dimensionless derivative of cmprsblty
	!C	CvRes:			(Cv-Cv^ig)/R, dimensionless measure of residual isochoric heat capacity.   
	!C	CpRes:			(Cp-Cp^ig)/R, dimensionless measure of residual isobaric heat capacity.   
	!C
	!C	CODED BY:   JRE, 4/25/2020
	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	SUBROUTINE FUGIDERVS(NC,xFrac,tKelvin,rhoMol_cc,zFactor,iErr)
	USE GlobConst
	!USE PhiDervs
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
	dimension xFrac(NC),FUGC(NC) !,bVol(NC)
	!zeroTol=1D-11  !now global
	iErr=0
	step=1.D-6
	!bVol=bVolCC_mol
	bMix=SumProduct(NC,xFrac,bVolCc_mol)
	eta=bMix*rhoMol_cc
	if(rhoMol_cc < zeroTol .or. tKelvin < zeroTol .or. eta < zeroTol .or. eta > etaMax)then
		if(LOUD)write(dumpUnit,*)'NUMDERVS: nonsense input. T,eta=',tKelvin, eta
		iErr=11
		return 
	endif
	isZiter=1 !We don't need the fugacity or uRes yet. 
	rhoPlus =rhoMol_cc*(1+step)
	rhoMinus=rhoMol_cc*(1-step)
	CALL FuVtot(isZiter,tKelvin,1/rhoPlus,xFrac,NC,FUGC,zPlus,aRes,uRes,iErrFu)
	if( iErrFu.ne.0 )then 
		if(LOUD)write(dumpUnit,*)'NUMDERVS: error from FuVtot on first call.'
		iErr=12
		return
	endif
	CALL FuVtot(isZiter,tKelvin,1/rhoMinus,xFrac,NC,FUGC,zMinus,aRes,uRes,iErrFu)
	rhoDZ_dRho=rhoMol_cc*(zPlus-zMinus)/(rhoPlus-rhoMinus)
	cmprsblty=zFactor+rhoDZ_dRho ! = (dP/dRho)/RT
	rho2DZ_dRho2=(zPlus+zMinus-2*zFactor)/(rhoMol_cc*step)**2
	RHOD2PDRHO2_RT=2*rhoDZ_dRho+rho2DZ_dRho2
	TPlus  =tKelvin*(1+step)
	TMinus=tKelvin*(1-step)
	CALL FuVtot(isZiter,TPlus,1/rhoMol_cc,xFrac,NC,FUGC,zPlus,aRes,uRes,iErrFu)
	uResPlus=uRes_RT*TPlus	!uRes_RT passed through GlobConst
	CALL FuVtot(isZiter,TMinus,1/rhoMol_cc,xFrac,NC,FUGC,zMinus,aRes,uRes,iErrFu)
	uResMinus=uRes_RT*TMinus
	TdZ_dT=tKelvin*(zPlus-zMinus)/( TPlus-TMinus )
	CvRes_R=(uResPlus-uResMinus)/( TPlus-TMinus )
	if(cmprsblty < zeroTol)then
		if(LOUD)write(dumpUnit,*)'FUGIDERVS: 0~cmprsblty=',cmprsblty
	endif
	CpRes_R=CvRes_R - 1 + (zFactor+TdZ_dT)**2/cmprsblty

	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CUBIC( A,Z)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
!C
!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!C
!C      PURPOSE --- THIS ROUTINE FINDS THE REAL ROOTS OF A
!C                  A CUBIC EQUATION: Z**3 + A(3)*Z**2 + A(2)*Z + A(1) = 0
!C                  THE SMALLEST ROOT IS RETURNED IN Z(1) AND THE LARGEST
!C                  IN Z(2)
!C                  THE INTERMEDIATE REAL ROOT (IF THERE IS ONE) IS DISCARDED
!C
!C      CODED BY:   J. F. ELY
!C                  NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
!C                  THERMOPHYSICS DIVISION
!C                  BOULDER, COLORADO  80303
!C
!C
!C
!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!C
      DIMENSION A(3), Z(3)
      LOGICAL initCall
      DATA initCall, PI /.TRUE., 3.14159265358D0 /
!C
!C                        ON FIRST ENTRY INITIALIZE CONSTANTS
!C
	IF (initCall)then
		zeroTol=1D-33
		ENTER = .FALSE.
		HPI = 0.5D0 * PI
		ADD1 = 2.0D0 * PI / 3.0D0
		ADD2 = 2.0D0 * ADD1
	endif
!C
!C                        SEE HOW MANY REAL ROOTS EXIST
	A0 = A(1)
	A1 = A(2)
	A2 = A(3) / 3.0D0
	Q = A1 / 3.0D0 - A2 * A2
	R = (A1 * A2 - A0) / 2.0D0 - A2 * A2 * A2
	TEST = Q * Q * Q + R * R
	IF (TEST < 0.0) then
		!C                        IRREDUCIBLE CASE, THREE REAL ALL DIFFERENT
		!C
		TEST = SQRT( -TEST)
		Q = 2.0D0 * DCBRT(SQRT(R*R+TEST*TEST))
		THETA = HPI
		IF( ABS(R) > zeroTol) THETA = ATAN(TEST/R)
		IF( THETA < 0) THETA = THETA + PI
		THETA = THETA / 3.0D0
		Z(1) = Q * COS(THETA) - A2
		Z(2) = Q * COS(THETA+ADD1) - A2
		Z(3) = Q * COS(THETA+ADD2) - A2
	else
		!C                        EITHER ONE REAL OR THREE REAL WITH TWO
		!C                        IDENTICAL ROOTS
		TEST=SQRT( TEST)
		S1 = DCBRT(R+TEST)
		S2 = DCBRT(R-TEST)
		Z(1) = S1 + S2 - A2
		Z(2) = Z(1)
		Z(3) = Z(1)
		IF ( TEST > 0) RETURN !one real (?)
		!C
		!!C                        THREE REAL WITH TWO IDENTICAL
		!C
		Z(2) = -0.5D0 * (S1+S2) - A2
		Z(3) = Z(2)
	endif
	!c
	!C                        ORDER THE ROOTS AND RETURN
	!c
	ZMIN = MIN(Z(1),Z(2),Z(3))
	Z(2) = MAX(Z(1),Z(2),Z(3))
	Z(1) = ZMIN
	IF (ZMIN < 0.0) Z(1) = Z(2)
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION FUNCTION DCBRT(X)
	IMPLICIT DOUBLEPRECISION ( A-H,O-Z)
	PARAMETER (THRD = 1.0D0/3.0D0)
	SIGN=1.0D0
	IF (X < 0) SIGN= -SIGN
	DCBRT=SIGN*ABS(X)**THRD
	RETURN
	END
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
	!USE ESDParms
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
	!USE ESDParms
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
