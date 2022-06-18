!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BootBP(tK,X,NC,INIT,pMPa,ITMAX,Y,ier,iOpt)
	!
	!  PURPOSE:  bootstrapping BUBBLE PRESSURE EVALUATIONS for difficult convergence.
	!  PROGRAMMED BY:  JRE 6/99
	!
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	double precision moFrac(NMX)
	INTEGER iterBoot,writeOpt,iOpt
	DIMENSION X(NMX),Y(NMX),ier(12)
	dimension tKCalc(NMX),pMPaCalc(NMX),etaLCalc(NMX),etaVCalc(NMX) 
	COMMON/eta/etaL,etaV,ZL,ZV
	writeOpt=0
	ier=0
	MAXIT=1111
	tWanted=tK
	pWanted=pMPa
	TCKay=0
	xMin=1
	xMax=0
	do i=1,nc
		TCKay=TCKay+X(i)*TC(i)
		if(X(i).gt.xMax)xMax=X(i)
		if(X(i).gt.xMax)iMax=i
		if(X(i).lt.xMin)xMin=X(i)
		moFrac(i)=X(i)
	ENDDO
	tK=0.7d0*TCKay-10
	init=0
	iterBoot=0
	delT=5
	pOld=1
	tOld=tK
1000 CONTINUE
	iterBoot=iterBoot+1
	tK=tK+delT					! In this version of bootBp, we vary T. Use BootBPx to vary x1 at const T.

	!Because near the critical point Z will not converge, it is better to stop it well before CR
	!It's a temporary solution for pure components, can be extended to mixtures if it is necessary 
	if (tK > (TC(1)*0.9d0).and.NC==1) then
		ier(1)=2
		RETURN
	endif
	 
	ITMAX=MAXIT
	CALL BUBPL(tK,X,NC,INIT,pMPa,ITMAX,Y,ier)

	!  if bubpl gave an error, maybe it's because our composition is too pure.
	if(ier(1).ne.0.and.xMax.ge.0.988)then
		moFrac(iMax)=.95
		do i=1,nc	
			if(i.ne.iMax)	then
				if(nc.eq.1) moFrac(i)=0.05/(nc)
			else
				moFrac(i)=0.05/(nc-1)
			endif
		ENDDO
		ITMAX=MAXIT
		CALL BUBPL(tK,moFrac,NC,INIT,pMPa,ITMAX,Y,ier)
		if(ier(1).ne.0)then
		write(*,*)'BootBP failed at xMax=0.95.'
			ier(1)=4
			goto 86
		endif

		moFrac(iMax)=.98
		do i=1,nc
			if(i.ne.iMax)moFrac(i)=0.02/(nc-1)
		ENDDO
		init=2
		ITMAX=MAXIT
		CALL BUBPL(tK,moFrac,NC,INIT,pMPa,ITMAX,Y,ier)
		if(ier(1).ne.0)then
			write(*,*)'BootBP failed at xMax=0.98.'
			ier(1)=4
			goto 86
		endif
		moFrac(iMax)=.99
		do i=1,nc
			if(i.ne.iMax)moFrac(i)=0.01/(nc-1)
		ENDDO
		CALL BUBPL(tK,moFrac,NC,INIT,pMPa,ITMAX,Y,ier)
		if(ier(1).ne.0)then
			write(*,*)'BootBP failed at xMax=0.99.'
			ier(1)=4
			goto 86
		endif
	endif

	!  the composition is not too pure.  check for other problems.
	IF(ier(1).ne.0)then
		if(iterBoot.le.1)then
			write(*,*)'BootBP failed on first try'
			ier(1)=2
			GOTO 86
		else if(delT.gt.2.0)then
			pMPa=pOld
			tK=tK-delT
			delT=0.5
			do iErrCount=1,11
				ier(iErrCount)=0
			ENDDO
			goto 1000
		else
			write(*,601)pMPa,tK
			if(iOpt.eq.0)write(61,601)pMPa,tK
			write(*,*)'t,p wanted',tWanted,pWanted
			ier(1)=3
			goto 86
		endif
	endif

	!  we got a convergence!  keep going with init=2 to bootstrap.
	init=2
	tKCalc(iterBoot)=tK
	pMPaCalc(iterBoot)=pMPa
	etaLCalc(iterBoot)=etaL
	etaVCalc(iterBoot)=etaV
	if (writeOpt==0) then	  !old Format
		if(iOpt==0)write(* ,602)tK,pMPa,X(1),Y(1),etaL,etaV,ier(1)
		if(iOpt==0)write(61,602)tK,pMPa,X(1),Y(1),etaL,etaV,ier(1)
	else  !P vs. T 
		if(iOpt==0)write(* ,603)tK,pMPa,etaL,etaV
		if(iOpt==0)write(61,603)tK,pMPa,etaL,etaV
	endif

	if(iOpt < 2)then
		if(tK < tWanted)then
			if(ABS(tWanted-tK).lt.delT.and.delT.gt.2)then
				delT=0.5
			endif
			pOld=pMPa
			tOld=tK
			goto 1000
		endif
		pMPa=pOld
		CALL BUBPL(tWanted,X,NC,INIT,pMPa,ITMAX,Y,ier)
		!call ErrCheck(ier,iFlag)
		if(ier(1).ne.0)then
			write(*,*)'disaster in BootBP on final calc'
		endif
	else if(iOpt.eq.2)then
		if(pMPa.lt.pWanted)then
			if(ABS(pWanted-pMPa).lt.0.05.and.delT.gt.2)then
				delT=0.5
			endif
			pOld=pMPa
			tOld=tK
			goto 1000
		endif
		tWanted=tOld+(tK-tOld)/(pMPa-pOld)*(pWanted-pOld)
		pMPa=pOld
		CALL BUBPL(tWanted,X,NC,INIT,pMPa,ITMAX,Y,ier)
		!call ErrCheck(ier,iFlag)
		tK=tWanted
		if(ier(1).ne.0)then
			write(*,*)'disaster in BootBP on final calc'
		endif
	endif
		
601	FORMAT(' BPMax=',(F9.4,1X),'at tK=',F9.2,' The requested pt is not 2-phase.')
602	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
603	FORMAT(1X,1(F7.2,2X,F10.6,1X),2(F9.5,1X))

86	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BootBPx(tK,X,NC,INIT,pMPa,ITMAX,Y,ier,writeOpt)
	!
	!  PURPOSE:  bootstrapping BUBBLE PRESSURE EVALUATIONS for difficult convergence.
	!  PROGRAMMED BY:  JRE 6/99
	!
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	double precision moFrac(NMX),xFrac(NMX)
	INTEGER ier(12),ierbp(12), writeOpt !,iOpt,iterBoot
	DoublePrecision X(NMX),Y(NMX),Psat(NMX)
	
	!dimension tKCalc(NMX),pMPaCalc(NMX),etaLCalc(NMX),etaVCalc(NMX) 
	COMMON/eta/etaL,etaV,ZL,ZV
	ier=0 ! cant use vector init for assumed size ier(*)
	if(NC .ne. 2)then
		ier(1)=1
		if(writeOpt)print*,'BootBPx: This routine is only for binaries. NC=',NC
		return
	endif
	if( Tc(1) > Tc(2) )then
		ier(1)=2
		if(writeOpt)write(*,'(a,2F9.2)')' BootBPx: Comp1 must be volatile. Tc()=',Tc(1),Tc(2)
	endif
	MAXIT=111
	tWanted=tK
	pWanted=pMPa
	xWanted=x(1) ! We assume that comp1 is the most volatile.
	pLast=0
	TCKay=0
	xMin=1
	xMax=0
	do i=1,nc
		TCKay=TCKay+X(i)*TC(i)
		moFrac(i)=X(i)
	ENDDO
	if ( tK > (TC(2)*0.9d0) ) then ! If comp2 is the heavy comp, and we are so close to its critical, then just give up.
		ier(1)=2
		RETURN
	endif
	CALL BubPSTART(tK,pIdSoln,x,Psat,NC)
	pLast=Psat(2)
	xLast=zeroTol
	y=0
	y(1:2)=x(1:2)*Psat(1:2)/pIdSoln
	sumy=SUM(y)
	y(1:2)=y(1:2)/sumy
	pNow=pLast
	init=2 ! Bootstrap the last guess from this point forward.
	do i=1,99
		xFrac(1)=i/100.d0
		xFrac(2)=1-xFrac(1)
		ITMAX=MAXIT
		CALL BUBPL(tK,xFrac,NC,INIT,pNow,ITMAX,Y,ierbp)	!get vapor pressure of pure heavy comp.
		if(ierbp(1).ne.0)then
			if(writeOpt)print*,'BootBPx: ierbp=',ierbp(1), ' at x1=',xFrac(1)
			print*,'BootBPx: ierbp=',ierbp(1), ' at x1=',xFrac(1)
			ier(1)=4
			goto 86
		endif
		if(xFrac(1) > xWanted)then
			ITMAX=MAXIT
			CALL BUBPL(tK,x,NC,INIT,pNow,ITMAX,Y,ierbp)	!get bp at input x value.
			if(ierbp(1).ne.0)then
				if(writeOpt)print*,'BootBPx: ierbp=',ierbp(1), ' on final check. I give up.'
				print*,'BootBPx: ierbp=',ierbp(1), ' on final check. I give up.'
				ier(1)=5
				goto 86
			endif
			pMPa=pNow
			return ! Done! 
		endif !pNow > pWanted
		pLast=pNow
		xLast=xFrac(1)
	enddo ! x1=0.01,0.99
		
601	FORMAT(' BPMax=',(F9.4,1X),'at tK=',F9.2,' The requested pt is not 2-phase.')
602	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
603	FORMAT(1X,1(F7.2,2X,F10.6,1X),2(F9.5,1X))

86	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BootDT(tK,X,NC,INIT,pMPa,ITMAX,Y,ier,iOpt)
	!
	!  PURPOSE:  bootstrapping Det Temperature EVALUATIONS for difficult convergence.
	!  PROGRAMMED BY:  JRE 6/99
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	double precision moFrac(NMX)
	DIMENSION X(NMX),Y(NMX),ier(12)
	COMMON/eta/etaL,etaV,ZL,ZV

	do i=1,11
		ier(i)=0
	ENDDO
	MAXIT=1111
	tWanted=tK
	pWanted=pMPa
	TCKay=0
	xMin=1
	xMax=0
	do i=1,nc
		TCKay=TCKay+X(i)*TC(i)
		if(X(i).gt.xMax)xMax=X(i)
		if(X(i).gt.xMax)iMax=i
		if(X(i).lt.xMin)xMin=X(i)
		moFrac(i)=X(i)
	ENDDO
	tK=0.7*TCKay-10
	pMPa=1
	init=0
	iterBoot=0
	delP=0.5
	pOld=pMPa
	tOld=tK
1000 CONTINUE
	iterBoot=iterBoot+1
	pMPa=pMPa+delP
	ITMAX=MAXIT
	CALL DEWTV(pMPa,Y,NC,INIT,tK,ITMAX,X,ier)
	!  if DEWTV gave an error, maybe it's because our composition is too pure.
	if(ier(1).ne.0.and.xMax.ge.0.988)then
		moFrac(iMax)=.95
		do i=1,nc
			if(i.ne.iMax)moFrac(i)=0.05/(nc-1)
		ENDDO
		ITMAX=MAXIT
		CALL DEWTV(pMPa,moFrac,NC,INIT,tK,ITMAX,X,ier)
		if(ier(1).ne.0)then
			write(*,*)'BootDT failed at xMax=0.95.'
			ier(1)=4
			goto 86
		endif

		moFrac(iMax)=.98
		do i=1,nc
			if(i.ne.iMax)moFrac(i)=0.02/(nc-1)
		ENDDO
		init=2
		ITMAX=MAXIT
		CALL DEWTV(pMPa,moFrac,NC,INIT,tK,ITMAX,X,ier)
		if(ier(1).ne.0)then
			write(*,*)'BootDT failed at xMax=0.98.'
			ier(1)=4
			goto 86
		endif
		moFrac(iMax)=.99
		do i=1,nc
			if(i.ne.iMax)moFrac(i)=0.01/(nc-1)
		ENDDO
		CALL DEWTV(pMPa,moFrac,NC,INIT,tK,ITMAX,X,ier)
		if(ier(1).ne.0)then
			write(*,*)'BootDT failed at xMax=0.99.'
			ier(1)=4
			goto 86
		endif
	endif

	!  the composition is not too pure.  check for other problems.
	IF(ier(1).ne.0)then
		if(iterBoot.le.1)then
			write(*,*)'BootDT failed on first try'
			ier(1)=2
			GOTO 86
		elseif(delP.gt.0.2)then
			!  maybe we're taking too big of steps.  back up and sneak in.
			pMPa=pOld
			tK=tOld
			delP=0.05
			do iErrCount=1,11
				ier(iErrCount)=0
			ENDDO
		else
			!  you have to give up sometimes
			write(*,601)pMPa,tK
			if(iOpt.eq.0)write(61,601)pMPa,tK
			write(*,*)'t,p wanted',tWanted,pWanted
			ier(1)=3
			goto 86
		endif
	endif
	!  we got a convergence!  keep going with init=2 to bootstrap.
	init=2
	if(iOpt.eq.0)write(* ,602)tK,pMPa,X(1),Y(1),etaL,etaV,ier(1)
	if(iOpt.eq.0)write(61,602)tK,pMPa,X(1),Y(1),etaL,etaV,ier(1)

	if(iOpt.le.1)then
		if(pMPa.lt.pWanted)then
			if(ABS(pWanted-pMPa).lt.delP.and.delP.gt.0.2)then
				delP=0.05
			endif
			pOld=pMPa
			tOld=tK
			goto 1000
		endif
		tK=tOld
		CALL DEWTV(pWanted,Y,NC,INIT,tK,ITMAX,X,ier)
		!call ErrCheck(ier,iFlag)
		if(ier(1).ne.0)then
			write(*,*)'disaster in BootDT on final calc'
		endif
	elseif(iOpt.eq.2)then
		if(tK.lt.tWanted)then
			if(ABS(tWanted-tK).lt.0.5.and.delP.gt.0.2)then
				delP=0.05
			endif
			pOld=pMPa
			tOld=tK
			goto 1000
		endif
		pWanted=pOld+(pMPa-pOld)/(tK-tOld)*(tWanted-tOld)
		tK=tOld
		CALL DEWTV(pWanted,Y,NC,INIT,tK,ITMAX,X,ier)
		!call ErrCheck(ier,iFlag)
		pMPa=pWanted
		if(ier(1).ne.0)then
			write(*,*)'disaster in BootDT on final calc'
		endif
	endif
		
601	FORMAT(' DPMax=',(F9.4,1X),'at tK=',F9.2,' The requested pt is not 2-phase.')
602	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
86	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      SUBROUTINE BUBPL(tKelvin,xFrac,nComps,init,pMpa,itMax,yFrac,ier)
	!   REVISION DATE:   AUG 02 F90
	!   REVISION DATE:   FEB 93 (FOR ESD COMPATIBILITY)
	!   REVISION DATE:   JAN 92 SJS-FOR -VE PRESSURES
	!   REVISION DATE:   SEPTEMBER 5, 1985
	!   PROGRAMMED BY:   J.R. ELLIOTT, JR. (JAN. 1983)
	!
	!   PURPOSE:  CALCULATE BUBBLE POINT PRESSURE OF liqUid BASED
	!      ON TEMPERATURE AND liqUid COMPOSITION.
	!
	!   ARGUMENTS:
	!
	!     INPUT:
	!        tCrit()     VECTOR CRITICAL TEMPERATURES OF THE COMPONENTS
	!        pCrit()     VECTOR CRITICAL PRESSURES OF THE COMPONENTS
	!        acen()   VECTOR acenTRIC FACTORS OF THE COMPONENTS
	!        id()     VECTOR STANDARD id NUMBERS OF THE COMPONENTS
	!        rGas     GAS CONSTANT (EG. 8.31434 CC-MPA/(GMOL-K) )
	!        T        ABSOLUTE TEMPERATURE
	!        xFrac()      VECTOR MOLE FRACTIONS IN THE liqUid PHASE
	!        nComps       NUMBER OF COMPONENTS
	!        init     PARAMETER FOR SPECifICATION OF WHETHER THE initIAL GUESS
	!                 IS PROVIDED BY THE USER OR SHOULD BE CALCULATED.
	!          init = 0 initIAL GUESS FOR P CALCULATED BY PSTART
	!          init = 1 initIAL GUESS FOR P PASSED FROM CALLING ROUTINE
	!
	!     INPUT/OUTPUT:
	!        P        OPTIONAL INPUT initIAL GUESS/OUTPUT CALCULATED
	!                   ABSOLUTE PRESSURE
	!        itMax    INPUT maxIMUM NUMBER OF iterATIONS PERMITTED.
	!                 THE RECOMMENDED VALUE IS 50.
	!                 OUTPUT itMax IS SET TO THE NUMBER OF iterATIONS PERFORMED
	!
	!     OUTPUT:
	!        yFrac()      VECTOR MOLE FRACTIONS IN THE VAPOR PHASE
	!        ier()    VECTOR ERROR PARAMETERS
	!          ier(1)=1     AT LEAST ONE OF ier(2)-ier(11) WAS NOT ZERO
	!          ier(2)=1     liqUid ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST iterATION
	!          ier(3)=1     VAPOR ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST iterATION
	!          ier(4)=4,5,6 TERMINAL ERROR RETURNED FROM FUGI CALCULATION
	!                           THE NUMBER TELLS WHICH COMPONENT OF FUGI'S
	!                           ERROR VECTOR WAS SIGNifICANT
	!                =4     NEGATIVE LOG CALCULATED
	!                =5     LOG OF FUGACITY COEFFICIENT CAUSES OVERFLOW
	!                =6     iterATION ON COMPRESSIBILITY FACTOR DID NOT CONVERGE
	!                =11    gOldenZ used on last iteration instead of real Z.
	!          ier(5)=1     CALCULATIONS DETERMINED VAPOR & liqUID ROOTS EQUAL
	!                           (TRIVIAL SOLUTION)
	!          ier(6)=1     FAILED TO CONVERGE IN itMax LOOPS
	!          ier(7)=1     AN iterATION WAS PERFORMED WITH NO IMPROVEMENT
	!                           IN THE OBJECTIVE FUNCTION
	!          ier(8)=1     AN ELEMENT OF tCrit, pCrit, OR X WAS SPECifIED
	!                          INCORRECTLY.
	!          ier(9)=1     THE T SPECifIED WAS LESS THAN ZERO.
	!          ier(10)=1    AN initIAL GUESS FOR P WAS SPECifIED BUT
	!                          IT WAS UNACCEPTABLE.
	!          ier(11)=1    THE VALUE FOR NC WAS GREATER THAN 10 OR itMax<1.
	!          ier(12)=1    THE FUGACITIES WERE NOT EQUAL ON LAST iterATION.
	!
	!   NOTE:               UNITS OF ALL THE INPUTS SHOULD BE
	!                       CONSISTENT WITH UNITS OF rGas.  EXCEPT
	!                       FOR THIS, THE USER MAY CHOOSE HIS OWN UNITS.
	!
	!   REQD. ROUTINES:
	!          PSTART, FUGI
	!
	!   SUBPROGRAM RESTRICTIONS:
	!          AS WRITTEN, THE maxIMUM NUMBER OF COMPONENTS
	!          THAT CAN BE CONSIDERED IS nmx.
	!
	!   GENERAL COMMENTS:
	!          THIS SUBROUTINE SOLVES THE OBJECTIVE FUNCTION GIVEN
	!          IN THE LiterATURE REFERENCE BY A SECANT iterATION
	!          ON THE PRESSURE VARIABLE.  THE SUBROUTINE CALLED
	!          FOR FUGACITY CALCULATIONS ("FUGI") CONFORMS TO THE
	!          SPECifICATIONS OF THE SOAVE EQUATION OF STATE GIVEN
	!          IN PROCEDURE 8D1.1.
	!
	!   METHOD RELIABILITY:
	!          THE AVERAGE ERRORS QUOTED BELOW ARE EXPECTED WHEN USING
	!          THE CORRELATIONS OF BINARY INTERACTION COEFFICIENT GIVEN
	!          IN CHAP. 8 OF THE API TECHNICAL DATA BOOK based on SRK EOS.
	!          SYSTEM TYPE                 AVERAGE PERCENT ERROR IN P
	!          HYDROCARBON-HYDROCARBON              4.3
	!          HYDROCARBON-HYDROGEN SULFIDE         4.8
	!          HYDROCARBON-NITROGEN                14.0
	!          HYDROCARBON-CARBON MONOXIDE          7.6
	!          HYDROCARBON-CARBON DIOXIDE           7.4
	!
	!   REFERENCES:
	!          ANDERSON, T.F.;  PRAUSNITZ, J.M.;  IND. ENG. CHEM. PROC.
	!          DES. DEV., 19:9-14, (1980).
	!
	!          PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!*******************************************************
	!
	!IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	USE GlobConst
    USE EsdParms
	IMPLICIT NONE
	integer init,itMax,ier,nComps
	doublePrecision tKelvin,xFrac,pMpa,yFrac,bestErr,bestPmpa,pIdSoln,pInit
    DoublePrecision bMix !,value  !for debugging
	integer iComp,liq,kLiq,maxIter,iter,ierFugi,iErrCheck,iTry !,iErrMix
	doublePrecision fugcLiq,fugcVap,vp,zLiq,zVap,sumy,expArg,calck,gObj,gOld,change,pOld,pNew ,etaV,etaL
	DIMENSION xFrac(nmx),yFrac(nmx),ier(12)
	DIMENSION fugcLiq(nmx),fugcVap(nmx),ierFugi(12),VP(nmx)
	COMMON/ETA/ETAL,ETAV,zLiq,zVap
	!
	!  CHECK INPUTS FOR ERRORS.
	!
	ier=0

	DO iComp=1,nComps
		if(Tc(iComp).Le.0)ier(8)=1
		if(Pc(iComp).Le.0)ier(8)=1
		if(xFrac(iComp) < 0 .OR. xFrac(iComp) > 1)ier(8)=1
	enddo
	if(tKelvin < 0)ier(9)=1
	if(init==1 .AND. pMpa < 0)ier(10)=1
	if(nComps > nmx)ier(11)=1
	if(itMax < 1)ier(11)=1
	DO iComp=8,12
		if(ier(iComp).NE.0)ier(1)=1
	enddo
	if(ier(1).NE.0)then
		GOTO 86
	endif

	!  GET A STARTING VALUE OF PRESSURE if NONE IS PROVIDED AND SET
	!  VAPOR FUGACITY COEFFICIENTS TO UNITY.

	CALL BubPSTART(tKelvin,pIdSoln,xFrac,vp,nComps)
	if(init==0)then
		pInit=pIdSoln
	else
		pInit=pMPa
	endif
	pMPa=pInit
    bMix=0
	DO iComp=1,nComps
		if(xFrac(iComp) < 1e-11)xFrac(iComp)=1e-11
		fugcVap(iComp)=0.D0
        bMix=bMix+xFrac(iComp)*bVolCC_mol(iComp)
	enddo

	!  if init=2, COMPUTE VAPOR FUGACITIES USING INPUT VAPOR COMPOSITION.

	if(init.EQ.2)THEN
		liq=0
		CALL Fugi(tKelvin,pMpa,yFrac,nComps,liq,fugcVap,zVap,ierFugi)
		if(ierFugi(1) > 0)then
			ier(1) = 13
			return
		endif
	END if
	kLiq=1
	maxIter=itMax
	iter=0

    !if(LOUD)write(*,'(a,11F10.2)')'BubPL: Vx=',(vx(iComp),iComp=1,nComps)
	!  BEGINNING OF MAIN iterATION LOOP
	!call QueryParMix(1,value,iErrMix) ! for debugging.
	ier(1)= -2 !initIALIZE ier TO iterATE ON pMpa,Y
	bestErr=1E11
	bestPmpa=1E11
	DO WHILE(ier(1) < 0)
		iter=iter+1
		if(iter > maxIter) THEN
			ier(1)=1
			ier(6)=1
            exit ! ends iteration
		ENDif
		itMax=iter

		if(kLiq.NE.0)THEN
			!  CALCULATE liqUid FUGACITIES.  CHECK FOR ERRORS.
			liq=1
			if (pMpa < 0 .AND. LOUD)WRITE(6,*)'BUBPL: P < 0 during iteration.!'
			CALL Fugi(tKelvin,pMpa,xFrac,nComps,liq,fugcLiq,zLiq,ierFugi)
			etaL=etaPass
			if(ierFugi(1).NE.0) THEN
				DO iErrCheck=4,11
					if(ierFugi(iErrCheck)==1)ier(4)=iErrCheck
				ENDDO
				ier(1)=ier(4)
				if(ierFugi(11).ne.1)GOTO 86
				if( ier(1) .ne. 0) goto 86
			END if
		END if
		!  CALCULATE VAPOR COMPOSITION AND ERROR IN OBJECTIVE FUNCTION.
		sumy=0
		DO iComp=1,nComps
			expArg=fugcLiq(iComp)-fugcVap(iComp)
			if(expArg > 222)then
				expArg=222
			endif
			CALCK=EXP(expArg)
			yFrac(iComp)=CALCK*xFrac(iComp)
			sumy=sumy + yFrac(iComp)
		enddo
		gObj=1 - sumy
		if(ABS(gObj) < bestErr)then
			bestErr=ABS(gObj)
			bestPmpa=pMPa
		endif
		if( iter > 1 ) THEN                 !  FOR HIGHER iterATIONS, A SECANT STEP IS TAKEN.
			if(gObj==gOld)THEN			
				ier(7)=1
				GOTO 86
			END if
			change=gObj*(pMpa-pOld)/(gObj-gOld)
			pOld=pMpa
			gOld=gObj
			!  LIMIT THE CHANGE IN Pressure SLIGHTLY.
			if(DABS(change/pOld) > 0.2D0)change=DSIGN(0.2D0,change)*pOld
			pNew=pOld - change
			pMpa=pNew
			if (pMpa < 0) THEN
				pMpa=EXP(pMpa)
				IF(LOUD)WRITE(6,*)'WARNING: P IS -VE IN BUBPL. SHORTENING STEP.'
			ENDif
			!  TEST FOR CONVERGENCE OF P AND Y.  if PRESSURE IS CONVERGED
			!  BUT Y IS NOtKelvin, SKIP RECALCULATION OF liqUid FUGACITIES AT
			!  THE START OF THE NEXT iterATION.
			kLiq=1
			if(DABS(change/pMpa) < 1.D-8 .or. DABS(change) < 1E-12)THEN
				kLiq=0
				ier(1)= -1
				if(DABS(sumy-1) < 1.D-7)ier(1)=0
			ENDif
		ELSE	!  FOR FIRST iterATION, OLD VALUES ARE COMPUTED.
			pOld=pMpa
			gOld=gObj
			!pMpa=pIdSoln*0.99D0	 ! This should give us a good chance to bracket the answer if water is a component.
			pMpa=pInit*0.98D0	 ! pIdSoln causes zVap=zLiq if too high.  Rely on bootstrapping..
		END if

		!  NOT CONVERGED.  GET NEW VAPOR FUGACITIES, CHECK FOR ERRORS.

		DO iComp=1,nComps
			yFrac(iComp)=yFrac(iComp)/sumy
		ENDDO
		liq=0
		iTry=0		  !EM. Check for KO bug
		ierFugi(1)=1 !force call of fugi for vapor at least once.
		DO WHILE(ierFugi(1).NE.0)
			CALl Fugi(tKelvin,pMpa,yFrac,nComps,liq,fugcVap,zVap,ierFugi)
			etaV=etaPass
			if(iTry > 5)exit ! break out of the loop  
			if(ierFugi(11)==1)then
				!  this error probably means that the pressure is too high to give a vapor root.
				pMPa=pMPa/1.1
				do iComp=1,nComps  ! Try enhancing composition of volatile components.
					expArg=fugcLiq(iComp)-fugcVap(iComp)
					if(expArg > 222)then
						expArg=222
					endif
					CALCK=EXP(expArg)
					if(calcK > 1)calcK = calcK *2
					if(calcK < 1)calcK = calcK /2
					yFrac(iComp)=xFrac(iComp)*CALCK
				enddo
				yFrac(1:nComps)=yFrac(1:nComps)/sum( yFrac(1:nComps) )
			else if(ierFugi(1).ne.0.and.iter > 1)then  !by now ierFugi(1) should be zero
				!  this error could occur from a bad guess on p in an iteration greater than 1.
				!  the following step will back up towards the previous (hopefully viable) pressure.
				change=change/2 !EM.?
				pMpa=pOld-change
			endif
			iTry=iTry+1
        ENDDO !WHILE ierFugi(1).NE.0
		if(ierFugi(1).NE.0)THEN
			DO iErrCheck=4,11
				if( ierFugi(iErrCheck) == 1)ier(4)=iErrCheck
			enddo
			if(ierFugi(11).ne.1)GOTO 86	! check 
		ENDif
		if(zVap.LE.0 .AND. LOUD)pause 'Bubp: zVap.LE.0'
		if(DABS((zVap-zLiq)/zVap) < 1.D-4) THEN
			ier(5)=1
		ENDif

	ENDDO !WHILE ier(1).LE.0 iteration on P.
 
	!  END OF MAIN iterATION LOOP

	!  PERFORM FINAL ERROR CHECKS.
 
86	continue
	if(ierFugi(1).ne.0)ier(1)=1
	if(ierFugi(2).EQ.1)ier(2)=1
	if(ierFugi(3).EQ.1)ier(3)=1
	if(ier(6)==1) pMPa=bestPmpa
	if(ierFugi(7).EQ.1)ier(7)=1
	if(ierFugi(11).EQ.1)ier(4)=11 ! e.g. no vapor root found on last P iteration. 
	DO iErrCheck=2,12
		if(ier(iErrCheck).NE.0)ier(1)=1
	enddo
	RETURN
	END
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE BubPSTART(tKelvin,pMpa,xFrac,vp,nComps)
	!  LATEST REVISION - JULY 26, 1985
	!  LATEST REVISION - FEB 13, 1993
	!					 AUG 6, 2002. F90 FORMAtKelvin, IMPLICIT NONE.
	!
	!  PURPOSE -
	!  CALCULATE AN initIAL GUESS FOR ROUTINE BUBPL
	!  FROM A GENERALIZED VAPOR PRESSURE EQUATION.
	!
	!  ARGUMENTS -
	!  Tc   - INPUT VECTOR OF CRITICAL TEMPERATURES
	!  Pc   - INPUT VECTOR OF CRITICAL PRESSURES
	!  acen - INPUT VECTOR OF acenTRIC FACTORS
	!  T    - INPUT TEMPERATURE FOR DESIRED CALCULATION
	!  P    - OUTPUT ESTIMATED PRESSURE
	!  X    - INPUT VECTOR OF MOLE FRACTIONS OF COMPONENTS
	!  nComps   - INPUT NUMBER OF COMPONENTS
	!
	!  GENERAL COMMENTS -
	!  ASSUMING THE LOG OF VAPOR PRESSURE IS LINEAR WITH
	!  RESPECT TO INVERSE TEMPERATURE, THE VAPOR PRESSURE
	!  IS ESTIMATED USING THE CRITICAL PRESSURE AND
	!  acenTRIC FACTOR TO FIX THE SLOPE AND INTERCEPT.
	!
	!  PROGRAMMED BY -
	!  J.R. ELLIOTtKelvin, JR.  (OCT. 1984)
	USE GlobConst
	IMPLICIT NONE
	integer iComp,nComps
	doublePrecision tKelvin,pMpa,xFrac(nmx),VP(nmx)
	doublePrecision TR,sumx
	sumx=SUM(xFrac)
	if(sumx < zeroTol)pause 'BubPStart: sumx < 0?'
	xFrac=xFrac/sumx
	if(nComps < 1 .or. tKelvin < zeroTol .or. ABS( SUM(xFrac)-1)>0.1 )print*,'BubPStart: bad input. NC,TK=',nComps,tKelvin
	pMpa=0
	DO iComp=1,nComps
		TR=tKelvin/Tc(iComp)
		VP(iComp)=10**( 7*(1+acen(iComp))*(1-1/TR)/3 )*Pc(iComp)
		pMpa=pMpa+xFrac(iComp)*vp(iComp)
	enddo
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BUBTL(T,X,NC,INIT,P,ITMAX,Y,IER)
	!
	!   REVISION DATE:   SEPTEMBER 5, 1985
	!   PROGRAMMED BY:   J.R. ELLIOTT, JR. (JULY 1985)
	!
	!   PURPOSE:  CALCULATE BUBBLE POINT TEMPERATURE OF LIQUID BASED
	!      ON PRESSURE AND LIQUID COMPOSITION.
	!
	!   ARGUMENTS:
	!
	!     INPUT:
	!        TC()     VECTOR CRITICAL TEMPERATURES OF THE COMPONENTS
	!        PC()     VECTOR CRITICAL PRESSURES OF THE COMPONENTS
	!        ACEN()   VECTOR ACENTRIC FACTORS OF THE COMPONENTS
	!        ID()     VECTOR STANDARD ID NUMBERS OF THE COMPONENTS
	!        RGAS     GAS CONSTANT (EG. 8.31434 CC-MPA/(GMOL-K) )
	!        P        ABSOLUTE PRESSURE
	!        X()      VECTOR MOLE FRACTIONS IN THE LIQUID PHASE
	!        N!       NUMBER OF COMPONENTS
	!        INIT     PARAMETER FOR SPECIFICATION OF WHETHER THE INITIAL GUESS
	!                 IS PROVIDED BY THE USER OR SHOULD BE CALCULATED.
	!          INIT = 0 INITIAL GUESS FOR T CALCULATED BY TSTART
	!          INIT = 1 INITIAL GUESS FOR T PASSED FROM CALLING ROUTINE
	!
	!     INPUT/OUTPUT:
	!        T        OPTIONAL INPUT INITIAL GUESS/OUTPUT CALCULATED
	!                   ABSOLUTE TEMPERATURE
	!        ITMAX    INPUT MAXIMUM NUMBER OF ITERATIONS PERMITTED.
	!                 THE RECOMMENDED VALUE IS 50.
	!                 OUTPUT ITMAX IS SET TO THE NUMBER OF ITERATIONS PERFORMED
	!
	!     OUTPUT:
	!        Y()      VECTOR MOLE FRACTIONS IN THE VAPOR PHASE
	!        IER()    VECTOR ERROR PARAMETERS
	!          IER(1)=1     AT LEAST ONE OF IER(2)-IER(11) WAS NOT ZERO
	!          IER(2)=1     LIQUID ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST ITERATION
	!          IER(3)=1     VAPOR ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST ITERATION
	!          IER(4)=4,5,6 TERMINAL ERROR RETURNED FROM FUGI CALCULATION
	!                           THE NUMBER TELLS WHICH COMPONENT OF FUGI'S
	!                           ERROR VECTOR WAS SIGNIFICANT
	!                =4     NEGATIVE LOG CALCULATED
	!                =5     LOG OF FUGACITY COEFFICIENT CAUSES OVERFLOW
	!                =6     ITERATION ON COMPRESSIBILITY FACTOR DID NOT CONVERGE
	!          IER(5)=1     CALCULATIONS DETERMINED VAPOR & LIQUID ROOTS EQUAL
	!                           (TRIVIAL SOLUTION)
	!          IER(6)=1     FAILED TO CONVERGE IN ITMAX LOOPS
	!          IER(7)=1     AN ITERATION WAS PERFORMED WITH NO IMPROVEMENT
	!                           IN THE OBJECTIVE FUNCTION
	!          IER(8)=1     AN ELEMENT OF TC, PC, OR X WAS SPECIFIED
	!                          INCORRECTLY.
	!          IER(9)=1     THE P SPECIFIED WAS LESS THAN ZERO.
	!          IER(10)=1    AN INITIAL GUESS FOR T WAS SPECIFIED BUT
	!                          IT WAS UNACCEPTABLE.
	!          IER(11)=1    THE VALUE FOR NC WAS GREATER THAN 10.
	!          IER(12)=1    THE VALUE FOR ITMAX WAS LESS THAN 1.
	!
	!   NOTE:               UNITS OF ALL THE INPUTS SHOULD BE
	!                       CONSISTENT WITH UNITS OF RGAS.  EXCEPT
	!                       FOR THIS, THE USER MAY CHOOSE HIS OWN UNITS.
	!
	!   REQD. ROUTINES:
	!          TSTART, FUGI, ESTACT, SRKNR
	!
	!   SUBPROGRAM RESTRICTIONS:
	!          AS WRITTEN, THE MAXIMUM NUMBER OF COMPONENTS
	!          THAT CAN BE CONSIDERED IS TEN.
	!
	!   GENERAL COMMENTS:
	!          THIS SUBROUTINE SOLVES THE OBJECTIVE FUNCTION GIVEN
	!          IN THE LITERATURE REFERENCE BY A SECANT ITERATION
	!          ON THE TEMPERATURE VARIABLE.  THE SUBROUTINE CALLED
	!          FOR FUGACITY CALCULATIONS ("FUGI") CONFORMS TO THE
	!          SPECIFICATIONS OF THE SOAVE EQUATION OF STATE GIVEN
	!          IN PROCEDURE 8D1.1.
	!
	!   METHOD RELIABILITY:
	!          THE AVERAGE ERRORS QUOTED BELOW ARE EXPECTED WHEN USING
	!          THE CORRELATIONS OF BINARY INTERACTION COEFFICIENT GIVEN
	!          IN CHAP. 8 OF THE API TECHNICAL DATA BOOK.
	!          SYSTEM TYPE                 AVERAGE PERCENT ERROR IN T
	!          HYDROCARBON-HYDROCARBON              4.3
	!          HYDROCARBON-HYDROGEN SULFIDE         4.8
	!          HYDROCARBON-NITROGEN                14.0
	!          HYDROCARBON-CARBON MONOXIDE          7.6
	!          HYDROCARBON-CARBON DIOXIDE           7.4
	!
	!   REFERENCES:
	!          ANDERSON, T.F.;  PRAUSNITZ, J.M.;  IND. ENG. CHEM. PROC.
	!          DES. DEV., 19:9-14, (1980).
	!
	!          PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!******************************************************
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(nmx),Y(nmx),IER(12)
	DIMENSION FUGCL(nmx),FUGCV(nmx),ierf(12),VP(nmx)
	
	! CHECK INPUTS FOR ERRORS.
	DO I=1,12
 		IER(I)=0
	enddo
	DO I=1,NC
		IF(TC(I).LT.0)IER(8)=1
		IF(PC(I).LT.0)IER(8)=1
		IF(X(I) .LT. 0 .OR. X(I) .GT. 1)IER(8)=0
	enddo
	IF(P.LT.0)IER(9)=1
	IF(INIT.EQ.1 .AND. T.LT.0)IER(10)=1
	IF(NC.GT.10)IER(11)=1
	IF(ITMAX.LT.1)IER(12)=1
	DO I=8,12
		IF(IER(I).NE.0)IER(I)=1
		IF(IER(I).NE.0)goto 95
	enddo

	!  GET A STARTING VALUE OF TEMPERATURE IF NONE IS PROVIDED AND SET
	!  VAPOR FUGACITY COEFFICIENTS TO UNITY.
	IF(INIT.EQ.0)CALL BubTSTRT(T,P,X,VP,NC)
	TI=.99/T
	DO  I=1,NC
		FUGCV(I)=0.D0
	enddo
	ZV=0.95
	ZL=0.0003
	KLIQ=1
	MAX=ITMAX
	!C
	!C  BEGINNING OF MAIN ITERATION LOOP
	iCount=0
	DO ITER=1,MAX
		ITMAX=ITER
13		continue !re-entry for bad T guess     
		T=1/TI
		IF(KLIQ.NE.0)THEN
			!C  CALCULATE LIQUID FUGACITIES.  CHECK FOR ERRORS.
			LIQ=1
			CALl Fugi(T,P,X,NC,LIQ,FUGCL,ZL,ierf)
			IF(ierf(1).NE.0.AND.ITER.GT.1)THEN
				DO IE=4,11
					IF(ierf(IE).EQ.1)IER(4)=IE
				enddo
			!c				WRITE(*,*)'ERROR FROM FUGI. ITER,T=',ITER,T
			if(ierf(11).ne.1)GOTO 86
			END IF
		END IF
		!C  CALCULATE VAPOR COMPOSITION AND ERROR IN OBJECTIVE FUNCTION.
		SUMY=0
		DO I=1,NC
			Y(I)=X(I)*EXP(FUGCL(I)-FUGCV(I))
			SUMY=SUMY+Y(I)
		enddo
		G=DLOG(SUMY)
		IF(ITER.EQ.1)THEN
			!C  FOR FIRST ITERATION, OLD VALUES ARE GENERATED
			TIOLD=TI
			GOLD=G
			!CCC     TI=TI/(1.0+G/5.0)
			TI=TI*1.01
			T=1/TI
			IF(ZL.GT.0.33.OR.ierf(1).NE.0)then
				iCount=iCount+1
				!c		  WRITE(*,*)'TI,T',TI,T
				!c			PAUSE 'WARNING - T GUESS TOO HIGH, I WILL TRY TO CORRECT'
				if(iCount.le.5)GOTO 13
			endif
			GOTO 25
		ELSE
			!C  FOR HIGHER ITERATIONS, A SECANT STEP IS TAKEN.
			IF(G.EQ.GOLD)THEN
				IER(7)=1
				GOTO 86
			ENDIF
			CHNG=G*(TI-TIOLD)/(G-GOLD)
		ENDIF
		TIOLD=TI
		GOLD=G
		!C  LIMIT THE CHANGE IN TEMPERATURE SLIGHTLY.
		IF(DABS(CHNG/TI).GT.0.02D0)CHNG=DSIGN(0.02D0,CHNG)*TI
		TINEW=TI - CHNG
		TI=TINEW
		!C  TEST FOR CONVERGENCE OF T AND Y.  IF TEMPERATURE IS CONVERGED
		!C  BUT Y IS NOT, SKIP RECALCULATION OF LIQUID FUGACITIES AT
		!C  THE START OF THE NEXT ITERATION.
		IF(DABS((TIOLD-TI)/TI).LE.1.D-7)THEN
			IF(DABS(G).LE.1.D-6)GOTO 86
			KLIQ=0
		ELSE
			KLIQ=1
		END IF
		!C  NOT CONVERGED.  GET NEW VAPOR FUGACITIES, CHECK FOR ERRORS.
25		continue !skip to point for 1st iteration
		DO I=1,NC
			Y(I)=Y(I)/SUMY
		enddo
		T=1/TI
		LIQ=0
		if(iter.gt.1)CALl Fugi(T,P,Y,NC,LIQ,FUGCV,ZV,ierf)
		IF(DABS((ZV-ZL)/ZV).LT.1.D-4) IER(5)=1
		IF(ierf(1).NE.0)THEN
		DO IE=4,11
			IF(ierf(IE).EQ.1)IER(4)=IE
		ENDDO
		if(ierf(11).ne.1)goto 86
		END IF
	enddo  !iter=1,max
	!C  END OF MAIN ITERATION LOOP
	!C
	!C  ONLY WAY FOR PROGRAM TO REACH THIS NEXT STATEMENT IS FOR
	!C  NUMBER OF ITERATIONS TO EXCEED ITMAX
	!C  THEREFORE INDICATE ERROR.
	IER(6)=1
	!C
	!C  PERFORM FINAL ERROR CHECKS.
86	continue
	IF(ierf(2).EQ.1)IER(2)=1
	IF(ierf(3).EQ.1)IER(3)=1
	IF(ierf(7).EQ.1)IER(7)=1
	IF(IER(5).EQ.1)THEN
		IF(ZV.LT.0.31)IER(3)=1
		IF(ZL.GT.0.31)IER(2)=1
	ENDIF
	DO IE=2,11
		IF(IER(IE).NE.0)IER(1)=1
	enddo
95	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BubTSTRT(T,P,X,VP,NC)
	!
	!  LATEST REVISION - JULY 29, 1985
	!  LATEST REVISION - FEB 13, 1993
	!
	!  PURPOSE -
	!  CALCULATE AN INITIAL GUESS FOR ROUTINES BUBTL
	!  FROM A GENERALIZED VAPOR PRESSURE EQUATION.
	!
	!  ARGUMENTS -
	!  TC()   - INPUT VECTOR OF CRITICAL TEMPERATURES
	!  PC()   - INPUT VECTOR OF CRITICAL PRESSURES
	!  ACEN() - INPUT VECTOR OF ACENTRIC FACTORS
	!  T      - OUTPUT ESTIMATED TEMPERATURE
	!  P      - INPUT PRESSURE FOR DESIRED CALCULATION
	!  VP()   - OUTPUT ESTIMATED VAPOR PRESSURES OF THE COMPONENTS AT
	!           THE GIVEN TEMPERATURE.
	!  X      - INPUT VECTOR OF MOLE FRACTIONS OF COMPONENTS
	!  N!     - INPUT NUMBER OF COMPONENTS
	!
	!  GENERAL COMMENTS -
	!  ASSUMING THE LOG OF VAPOR PRESSURE IS LINEAR WITH
	!  RESPECT TO INVERSE TEMPERATURE, THE VAPOR PRESSURE
	!  IS ESTIMATED USING THE CRITICAL PRESSURE AND
	!  ACENTRI! FACTOR TO FIX THE SLOPE AND INTERCEPT.
	!  BASED ON THE SAME FUNCTIONAL FORM FOR THE MIX, LOG(P)=A-B/T
	!
	!  PROGRAMMED BY -
	!  J.R. ELLIOTT, JR.  (JULY 1985)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(*),VP(*)
	A=0
	B=0
	DO I=1,NC
		PR=P/PC(I)
		TSAT =TC(I)/( 1-3.D0/7.D0*DLOG10(PR)/(1+ACEN(I)) )
		A=A+X(I)*( DLOG10(PC(I))+7/3.D0*(1+ACEN(I)) )
		B=B+X(I)*7/3.D0*(1+ACEN(I))*TC(i)
		T=T+X(I)*TSAT
	enddo
	T=-B/(LOG10(P)-A)
	TI=1/T
	PCALC=0
	DO I=1,NC
		TR=T/TC(I)
		VP(I)=PC(I)*10**(7.D0/3.D0*(1+ACEN(I))*(1-1/TR))
		PCALC=PCALC+X(I)*VP(I)
	enddo
	FOLD=LOG(P/PCALC)
	TIOLD=TI
	TI=1.1*TIOLD
	KOUNT=0
1000 KOUNT=KOUNT+1
     if(LOUD)then
	    IF(KOUNT.GT.25)PAUSE'ERROR IN BTSTRT, NO CONVERGE'
     end if
	IF(KOUNT.GT.25)RETURN
	T=1/TI
	PCALC=0
	DO I=1,NC
		TR=T/TC(I)
		VP(I)=PC(I)*10**(7.D0/3.D0*(1+ACEN(I))*(1-1/TR))
		PCALC=PCALC+X(I)*VP(I)
	enddo
	F=LOG(P/PCALC)
	CHNG=F/(F-FOLD)*(TI-TIOLD)
	TIOLD=TI
	FOLD=F
	!WRITE(*,*)'IN BTSTRT. T,CHNG',T,CHNG
	TI=TI-CHNG
	IF(ABS(CHNG/TI).GT.1.D-4)goto 1000
	T=1/TI
	!WRITE(*,*)'IN BTSTRT. T=',T
	DO I=1,NC
		TR=T/TC(I)
		VP(I)=PC(I)*10**(7.D0/3.D0*(1+ACEN(I))*(1-1/TR))
	enddo
	RETURN
	END
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
	!$ DEWTV
	SUBROUTINE DEWTV(P,Y,NC,INIT,T,ITMAX,X,IER)
	!
	!   REVISION DATE:   SEPTEMBER 5, 1985
	!   PROGRAMMED BY:   J.R. ELLIOTT, JR. (JULY 1985)
	!
	!   PURPOSE:  CALCULATE DEW POINT TEMPERATURE OF VAPOR BASED
	!      ON PRESSURE AND VAPOR COMPOSITION.
	!
	!   ARGUMENTS:
	!
	!     INPUT:
	!        TC()     VECTOR CRITICAL TEMPERATURES OF THE COMPONENTS
	!        PC()     VECTOR CRITICAL PRESSURES OF THE COMPONENTS
	!        ACEN()   VECTOR ACENTRI! FACTORS OF THE COMPONENTS
	!        ID()     VECTOR STANDARD ID NUMBERS OF THE COMPONENTS
	!        RGAS     GAS CONSTANT (EG. 8.31434 CC-MPA/(GMOL-K) )
	!        P        ABSOLUTE PRESSURE
	!        Y()      VECTOR MOLE FRACTIONS IN THE VAPOR PHASE
	!        N!       NUMBER OF COMPONENTS
	!        INIT     PARAMETER FOR SPECIFICATION OF WHETHER THE INITIAL GUESS
	!                 IS PROVIDED BY THE USER OR SHOULD BE CALCULATED.
	!          INIT = 0 INITIAL GUESS FOR T CALCULATED BY TSTART
	!          INIT = 1 INITIAL GUESS FOR T PASSED FROM CALLING ROUTINE
	!
	!     INPUT/OUTPUT:
	!        T        OPTIONAL INPUT INITIAL GUESS/OUTPUT CALCULATED
	!                   ABSOLUTE TEMPERATURE
	!        ITMAX    INPUT MAXIMUM NUMBER OF ITERATIONS PERMITTED.
	!                 THE RECOMMENDED VALUE IS 50.
	!                 OUTPUT ITMAX IS SET TO THE NUMBER OF ITERATIONS PERFORMED
	!
	!     OUTPUT:
	!        X()      VECTOR MOLE FRACTIONS IN THE LIQUID PHASE
	!        IER()    VECTOR ERROR PARAMETERS
	!          IER(1)=1     AT LEAST ONE OF IER(2)-IER(11) WAS NOT ZERO
	!          IER(2)=1     LIQUID ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST ITERATION
	!          IER(3)=1     VAPOR ROOT PASSED FROM FUGI WAS NOT REAL
	!                           ON LAST ITERATION
	!          IER(4)=4,5,6 TERMINAL ERROR RETURNED FROM FUGI CALCULATION
	!                           THE NUMBER TELLS WHICH COMPONENT OF FUGI'S
	!                           ERROR VECTOR WAS SIGNIFICANT
	!                =4     NEGATIVE LOG CALCULATED
	!                =5     LOG OF FUGACITY COEFFICIENT CAUSES OVERFLOW
	!                =6     ITERATION ON COMPRESSIBILITY FACTOR DID NOT CONVERGE
	!          IER(5)=1     CALCULATIONS DETERMINED VAPOR & LIQUID ROOTS EQUAL
	!                           (TRIVIAL SOLUTION)
	!          IER(6)=1     FAILED TO CONVERGE IN ITMAX LOOPS
	!          IER(7)=1     AN ITERATION WAS PERFORMED WITH NO IMPROVEMENT
	!                           IN THE OBJECTIVE FUNCTION
	!          IER(8)=1     AN ELEMENT OF TC, PC, OR X WAS SPECIFIED
	!                          INCORRECTLY.
	!          IER(9)=1     THE P SPECIFIED WAS LESS THAN ZERO.
	!          IER(10)=1    AN INITIAL GUESS FOR T WAS SPECIFIED BUT
	!                          IT WAS UNACCEPTABLE.
	!          IER(11)=1    THE VALUE FOR N! WAS GREATER THAN 10.
	!          IER(12)=1    THE VALUE OF ITMAX WAS LESS THAN 1.
	!
	!   NOTE:               UNITS OF ALL THE INPUTS SHOULD BE
	!                       CONSISTENT WITH UNITS OF RGAS.  EXCEPT
	!                       FOR THIS, THE USER MAY CHOOSE HIS OWN UNITS.
	!
	!   REQD. ROUTINES:
	!          TSTART, FUGI, ESTACT, SRKNR
	!
	!   SUBPROGRAM RESTRICTIONS:
	!          AS WRITTEN, THE MAXIMUM NUMBER OF COMPONENTS
	!          THAT CAN BE CONSIDERED IS TEN.
	!
	!   GENERAL COMMENTS:
	!          THIS SUBROUTINE SOLVES THE OBJECTIVE FUNCTION GIVEN
	!          IN THE LITERATURE REFERENCE BY A SECANT ITERATION
	!          ON THE TEMPERATURE VARIABLE.  THE SUBROUTINE CALLED
	!          FOR FUGACITY CALCULATIONS ("FUGI") CONFORMS TO THE
	!          SPECIFICATIONS OF THE SOAVE EQUATION OF STATE GIVEN
	!          IN PROCEDURE 8D1.1.
	!
	!   METHOD RELIABILITY:
	!          THE AVERAGE ERRORS QUOTED BELOW ARE EXPECTED WHEN USING
	!          THE CORRELATIONS OF BINARY INTERACTION COEFFICIENT GIVEN
	!          IN CHAP. 8 OF THE API TECHNICAL DATA BOOK.
	!          SYSTEM TYPE                 AVERAGE PERCENT ERROR IN T
	!          HYDROCARBON-HYDROCARBON              4.3
	!          HYDROCARBON-HYDROGEN SULFIDE         4.8
	!          HYDROCARBON-NITROGEN                14.0
	!          HYDROCARBON-CARBON MONOXIDE          7.6
	!          HYDROCARBON-CARBON DIOXIDE           7.4
	!
	!   REFERENCES:
	!          ANDERSON, T.F.;  PRAUSNITZ, J.M.;  IND. ENG. CHEM. PROC.
	!          DES. DEV., 19:9-14, (1980).
	!
	!          PROCEDURE 8D1.1 OF TECHNICAL DATA BOOK.
	!
	!*******************************************************
	USE GlobConst !Tc,Pc,acen,id...
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(NC),Y(NC),IER(12)
	DIMENSION FUGCL(NMX),FUGCV(NMX),ierf(12),VP(10)
	COMMON/ETA/ETAL,ETAV,ZL,ZV
	!C
	!C  CHECK INPUTS FOR ERRORS.
	!C
	DO I=1,12
		IER(I)=0
	enddo
	DO I=1,NC
		IF(TC(I).LT.0)IER(8)=1
		IF(PC(I).LT.0)IER(8)=1
		IF(Y(I).LT.0 .OR. Y(I).GT.1)IER(8)=1
		if(Y(i).lt.1e-11)Y(i)=1e-11
	enddo
	IF(P.LT.0)IER(9)=1
	IF(INIT.EQ.1 .AND. T.LT.0)IER(10)=1
	IF(NC.GT.10)IER(11)=1
	IF(ITMAX.LT.1)IER(12)=1
	DO I=8,12
		IF(IER(I).NE.0)IER(1)=1
	enddo
	IF(IER(1).NE.0)GOTO 86
	!
	!C  GET A STARTING VALUE OF TEMPERATURE IF NONE IS PROVIDED AND SET
	!C  LIQUID FUGACITY COEFFICIENTS TO UNITY.
	if(init.le.1)CALL DewTSTRT(TEST,P,Y,VP,NC,X,INIT)
	IF(INIT.EQ.0)T=TEST
	TI=1.01/T
	KLIQ=1
	MAX=ITMAX
	!C
	!C  BEGINNING OF MAIN ITERATION LOOP
	DO ITER=1,MAX
13		continue !re-entry for bad T guess
		T=1.D0/TI
		ITMAX=ITER
		IF(KLIQ.NE.0)THEN
			!  CALCULATE VAPOR FUGACITIES.  CHECK FOR ERRORS.
			LIQ=0
			CALl Fugi(T,P,Y,NC,LIQ,FUGCV,ZV,ierf)
			IF(ierf(1).NE.0)THEN
				DO IE=4,12
				  IF(ierf(IE).EQ.1)IER(4)=IE
				enddo
				if(ierf(11).ne.1)goto 86
			END IF
		END IF
		IF(ITER.EQ.1)THEN
			LIQ=1
			CALl Fugi(T,P,X,NC,LIQ,FUGCL,ZL,ierf)
		END IF
		!C  CALCULATE LIQUID COMPOSITION AND ERROR IN OBJECTIVE FUNCTION.
		SUMX=0
		DO I=1,NC
			X(I)=Y(I)*EXP(FUGCV(I)- FUGCL(I))
			SUMX=SUMX + X(I)
		enddo
		G=DLOG(SUMX)
		IF(ITER.EQ.1)THEN
			!C  FOR FIRST ITERATION, OBTAIN OLD VALUES.
			TIOLD=TI
			GOLD=G
			TI=TI*.99D0
			T=1/TI
			IF(ZV.LT.0.31.OR.ierf(1).NE.0)WRITE(*,*)'TI,T',TI,T
			IF(ZV.LT.0.31.OR.ierf(1).NE.0)then
				if(LOUD)PAUSE 'WARNING - T GUESS TOO LOW,I WILL TRY TO CORRECT'
			endif
			IF(ZV.LT.0.31.OR.ierf(1).NE.0)goto 13
			goto 25
		ELSE
			!  FOR HIGHER ITERATIONS, A SECANT STEP IS TAKEN.
			IF(G.EQ.GOLD)THEN
				IER(7)=1
				goto 86
			END IF
			CHNG=G*(TI-TIOLD)/(G-GOLD)
		END IF
		IF(DABS(CHNG/TI).GT.0.1D0)CHNG=DSIGN(0.10D0,CHNG)*TI
		TIOLD=TI
		GOLD=G
		TINEW=TI - CHNG
		TI=TINEW
		!  TEST FOR CONVERGENCE OF T AND X.  IF TEMPERATURE IS CONVERGED
		!C  BUT X IS NOT, SKIP RECALCULATION OF VAPOR FUGACITIES AT
		!C  THE START OF THE NEXT ITERATION.
		IF(DABS((TIOLD-TI)/TI).LE.1.D-6)THEN
			IF(DABS(G).LE.1.D-5)goto 86
			KLIQ=0
		ELSE
			KLIQ=1
		END IF
		!  NOT CONVERGED.  GET NEW LIQUID FUGACITIES, CHECK FOR ERRORS.
25		Continue !skip point for 1st iteration
		DO I=1,NC
			X(I)=X(I)/SUMX
		enddo
		T=1.D0/TI
		LIQ=1
		CALl Fugi(T,P,X,NC,LIQ,FUGCL,ZL,ierf)
		IF(DABS((ZV-ZL)/ZV).LT.1.D-2) IER(5)=1
		IF(ierf(1).NE.0)THEN
			DO IE=4,12
				IF(ierf(IE).EQ.1)IER(4)=IE
			enddo
			WRITE(*,*)'ERROR IN FUGI. T=',T
			if(ierf(11).ne.1)goto 86
		END IF
	enddo !iter=1,max
	!C  END OF MAIN ITERATION LOOP
	!C
	!C  ONLY WAY FOR PROGRAM TO REACH THIS NEXT STATEMENT IS FOR
	!C  NUMBER OF ITERATIONS TO EXCEED ITMAX
	!C  THEREFORE INDICATE ERROR.
      IER(6)=1
	!C
	!C  PERFORM FINAL ERROR CHECKS.
86	IF(ierf(2).EQ.1)IER(2)=1
	IF(ierf(3).EQ.1)IER(3)=1
	IF(IER(5).EQ.1)THEN
		IF(ZV.LT.0.31)IER(3)=1
		IF(ZL.GT.0.31)IER(2)=1
	ENDIF
	DO IE=2,12
		IF(IER(IE).NE.0)IER(1)=1
	enddo
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE DewTSTRT(T,P,Y,VP,NC,X,INIT)
	!C
	!C  LATEST REVISION - FEB 13, 1993
	!C
	!C  PURPOSE -
	!C  CALCULATE AN INITIAL GUESS FOR ROUTINE DEWTV
	!C   FROM A GENERALIZED VAPOR PRESSURE EQUATION.
	!C
	!C  ARGUMENTS -
	!C  TC   - INPUT VECTOR OF CRITICAL TEMPERATURES
	!C  PC   - INPUT VECTOR OF CRITICAL PRESSURES
	!C  ACEN - INPUT VECTOR OF ACENTRIC FACTORS
	!C  T    - INPUT TEMPERATURE FOR DESIRED CALCULATION
	!C  P    - OUTPUT ESTIMATED PRESSJRE
	!C  Y    - INPUT VECTOR OF MOLE FRACTIONS OF COMPONENTS
	!C  NC   - INPUT NUMBER OF COMPONENTS
	!C
	!C  GENERAL COMMENTS -
	!C  ASSUMING THE LOG OF VAPOR PRESSURE IS LINEAR WITH
	!C  RESPECT TO INVERSE TEMPERATURE, THE VAPOR PRESSURE
	!C  IS ESTIMATED USING THE CRITICAL PRESSURE AND
	!C  ACENTRIC FACTOR TO FIX THE SLOPE AND INTERCEPT.
	!C
	!C  PROGRAMMED BY -
	!C  J.R. ELLIOTT, JR.  (FEB. 1993)
	USE GlobConst !{Tc,Pc,RGAS, ...}
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION Y(*),VP(*),X(*)
	TSTART=0
	DO I=1,NC
		TSAT=TC(I)/( 1-3*LOG10(P/PC(I))/7/(1+ACEN(I)) )
		TSTART=TSTART+Y(I)*TSAT
	enddo
	TOLD=TSTART
	PINV=0
	DO I=1,NC
		TR=TOLD/TC(I)
		VP(I)=PC(I)*10**( 7*(1+ACEN(I))*(1-1/TR)/3)
		PINV=PINV+Y(I)/VP(I)
	enddo
	FOLD=1/P-PINV
	T=TOLD*0.9
	DO ITER=1,25
		PINV=0
		DO I=1,NC
			TR=T/TC(I)
			VP(I)=10**( 7*(1+ACEN(I))*(1-1/TR)/3 )*PC(I)
			PINV=PINV+Y(I)/VP(I)
			X(I)=Y(I)*p/VP(I)
		enddo
		F=1/P-PINV
		CHNG=F*(T-TOLD)/(F-FOLD)
		FOLD=F
		TOLD=T   
		!WRITE(*,*)'IN DTSTRT: T,F',T,F
		T=T-CHNG
		IF(ABS(CHNG/T).LT.1.E-4)RETURN
	enddo
	IF(INIT.EQ.1)RETURN
	WRITE(*,*)'ERROR - IDEAL DEWT ESTIMATE DID NOT CONVERGE '
	WRITE(*,*)'TRY ENTERING YOUR OWN INITIAL GUESS OR GIVE UP'
	WRITE(*,*)'     TC       PC      ACEN        Y'
	DO I=1,NC
		WRITE(*,601)TC(I),PC(I),ACEN(I),Y(I)
	enddo
601	FORMAT(1X,4(F10.4,1X))
	return
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE FlashEar(tKelvin,pMpa,NC,iPhaseUpPas,zFeed,calck,itMax,xFrac,yFrac,vof,iErrCode)
	!C  FLASH ROUTINE FOR VAPOR-LIQUID OR LIQUID-LIQUID	Using Equal Area Rule
	!C	Ref: Y.Tang and G.Stephenson, AIChE J., 2011
	!C
	!C     REQD. ROUTINES:
	!C                      Fugi, 
	!C
	!integer NMX,iComp,nComps,iPhaseUp,itMax,maxIt,iErrCode,ierVof,iter,iFlag,iterVof
	USE GlobConst
	IMPLICIT DoublePrecision(A-H,O-Z)
	integer ierFugi(12) !,ID,iTry,iWarn
	doublePrecision calck(NMX),xFrac(NMX),yFrac(NMX),zFeed(NMX)
	doublePrecision fugcVap(NMX),fugcLiq(NMX)
	!doublePrecision chemPo(NMX),FUGC(NMX),calckNew(NMX),fugcPure(NMX)
	!doublePrecision vof,tKelvin,pMpa,rGas,etaVap,etaLiq,zVap,zLiq,sumx,sumy,devMax,devi
	!doublePrecision difChemPo,TC,PC,ACEN,biggestKvalue,gEx,xInf,picardFactor,devIold
	!doublePrecision rLnGam1,rLnGam2,denom
	!doublePrecision etaPasUp,etaPasLo,zPasUp,zPasLo,zDum1,zDum2
	character*88 errMsg(11)
	common/ETA/etaLiq,etaVap,zDum1,zDum2
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	iErrCode=0
	tK=tKelvin
	iWarn=0
	maxIt=itMax
	iter=0
	iTry=0
	iFlag=1
	xInf=1e-11
	do i=1,2
		xFrac(i)=xInf
	enddo
	zVap=1 !Avoid divide zero error when failing before FUGI(zVap)
	deviOld=1000 !initialize detection for unstable successive substitution
	errMsg(1)='No spinodal max/min'
	errMsg(2)='Practically at a flat spinodal.'
	errMsg(3)='Liquid FUGI call failed on last iteration'
	errMsg(4)='Vapor  FUGI call failed on last iteration'
	write(*,*)'EAR method'
	itMax=44
	call SpinodalMixEz(NC,tKelvin,PMPa,x1FracLo,gTotLo,dG_dxLo,iPhaseLo,x1FracHi,gTotHi,dG_dxHi,iPhaseHi,iErrSpin) 
	write(*,*)'iPhaseLo,iPhaseHi,x1FracLo,dG_dxLo,gTotLo,x1FracHi,dG_dxHi,gTotHi'
	write(*,'(2i3,8F10.3)')iPhaseLo,iPhaseHi,x1FracLo,dG_dxLo,gTotLo,x1FracHi,dG_dxHi,gTotHi
	if(iErrSpin.ne.0)then
		iErrCode=1
		write(*,'(a,f8.1,f9.3,a)')' FlashEar: no minimax at T,P = ',tK,pMPa,' T > Tc?'
		if(LOUD)pause 'check spinodal values'
		goto 86
	endif
	if(ABS(dG_dxHi-dG_dxLo) .lt. 1D-3)then !the difference is too small. It means dG_dx is practically flat so it is the critical point.
		iErrCode=2
		write(*,*)'T,P=',tKelvin,pMPa
		write(*,*)'The spinodal is practically flat. Exiting FlashEar.'
		if(LOUD)pause 'check spinodal values'
		goto 86
	endif
	dG_dx=( gTotHi-gTotLo )/(x1FracHi-x1FracLo)	!EAR method of Eubank and Hall (1995)
	if(dG_dx > MAX(dG_dxHi,dG_dxLo) .or. dgTry < MIN(dG_dxLo,dG_dxHi) )dG_dx=(dG_dxHi+dG_dxLo)/2
	x1SpinLo=x1FracLo !Store this in case there is trouble during iteration and need to restart.
	x1FracLo=x1FracLo - 0.1
	x1FracUp=x1FracHi + 0.1
	if(x1FracLo > x1FracHi)then
		x1FracLo=x1FracHi - 0.1
		x1FracUp=x1FracLo + 0.1
	endif
	if(x1FracLo < 0)x1FracLo=1d-33
	if(x1FracUp > 1)x1FracUp=1-1d-11
	itMax=22
	do iter=1,itMax
		Call EarMixSolver(tK,pMPa,dG_dx,x1FracLo,gTotLo,iPhaseLo,ierEarSolver)
		if(ierEarSolver)then
			iErrCode=3
			write(*,*)'dG_dx,x1FracLo',dG_dx,x1FracLo
			if(LOUD)pause 'FlashEar: EarSolver failed'
			goto 86
		endif 
		Call EarMixSolver(tK,pMPa,dG_dx,x1FracUp,gTotUp,iPhaseUp,ierEarSolver) 
		if(LOUD.and.ierEarSolver)pause 'FlashEar: EarSolver failed'
		dgTest=( gTotUp-gTotLo )/(x1FracUp-x1FracLo)	!EAR method of Eubank and Hall (1995)
		change=dgTest-dG_dx
		dgTry=dG_dx+change
		if(dgTry > MAX(dG_dxHi,dG_dxLo) .or. dgTry < MIN(dG_dxLo,dG_dxHi) )then
			change=SIGN(1.d-4,change)
			dG_dx=MAX(dG_dxHi,dG_dxLo)-change*iter
			if(change < 0)dG_dx=MIN(dG_dxHi,dG_dxLo)+change*iter  !get a fresh start with a little change
			!iter=itMax-1			   !do one last iteration, but this is illegal, so just keep restarting until itMax exceeded.
			x1FracLo=x1SpinLo-0.1
			x1FracUp=x1FracHi+0.1
			if(x1FracLo < 0)x1FracLo=1d-33
			if(x1FracUp > 1)x1FracUp=1-1d-11
			cycle
		endif
		dG_dx=dG_dx+change
		if(ABS(change/dG_dx).lt.1d-5)EXIT
	enddo
	if(ABS(change/dG_dx).gt.1d-5)then
		iErrCode=4
		write(*,'(a,i3,a,2e12.4)')' FlashEar error: exceeded max iterations = ',itMax,' at T,P = ',tK,pMPa
		if(LOUD)pause 'FlashEar: Check spinodal values'
		goto 86
	endif

	itMax=iter

	iPhaseUpSto=iPhaseUpPas !just in case meaningful info there.
	iPhaseUpPas=iPhaseUp+10*iPhaseLo
	iSwap=0
	if(iPhaseUpPas .eq. 1)iSwap=1 !it means they arent the same and x1FracLo is for a vapor, while x1FracUp is for a liquid
	xFrac(1)=x1FracLo
	xFrac(2)=1-xFrac(1)		
	call Fugi(tK,pMPa,xFrac,NC,iPhaseLo,fugcLiq,zLiq,ierFugi) !one last call to fugi for chemPo and debugging.
	uSatL=DUONKT
	yFrac(1)=x1FracUp
	yFrac(2)=1-yFrac(1)		
	call Fugi(tK,pMPa,yFrac,NC,iPhaseUp,fugcVap,zVap,ierFugi) !one last call to fugi for chemPo and debugging.
	uSatV=DUONKT
	if(iSwap)then
		yFrac(1)=x1FracLo
		xFrac(1)=x1FracUp
		xFrac(2)=1-xFrac(1)		
		yFrac(2)=1-yFrac(1)		
	endif
	do i=1,2
		calcK(i)=fugcLiq(i)/fugcVap(i)
		if(iSwap)calcK(i)=fugcVap(i)/fugcLiq(i)
	enddo
	vof=( zFeed(1)-xFrac(1) )/(yFrac(1)-xFrac(1)+1e-33 )
	write(*,*)'FlashEar: iPhase=',iPhaseUpPas

86	continue
	RETURN
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine EarMixSolver(tK,pMPa,dG_dxTarget,x1Frac,gTot,iPhase,ier) 
	USE GlobConst !{Tc,Pc,...,Rgas, ...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION xFrac(NMX)
	dimension FUGCPURE(NMX),fugcVap(NMX),fugcLiq(NMX),ierFugi(12) 
	!dimension vMolec(NMX),gTotStor(100) !,dG_dxStor(100)
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT

	ier=0
	NC=2 !Not allowed if NC.ne.2 (for now)
	tKelvin=tK

	xInf = 1e-11
	do i=1,2
		xFrac(i)=xInf
	enddo
	nComps=NC
	do i=1,2
		xFrac(i)=1
		call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
		call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
		fugcPure(i)=fugcVap(i)
		if( fugcLiq(i) < fugcVap(i) .and. ierFugi(1)==0 )fugcPure(i)=fugcLiq(i) ! min G is most stable.
		xFrac(i)=xInf
	enddo

	xFrac(1)=x1Frac
	xFrac(2)=1-xFrac(1)
	call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
	chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
	gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	ierLiq=ierFugi(1)
	if(zLiq < 0.3 .and. zLiq>0)gTotLiq=gTotLiq-0.001 !give it a nudge to avoid randomness when in the one root region. We don't use it directly anyway.
	call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
	chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
	gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	iPhaseLo=0
	ierCheck=ierFugi(1)
	if(gTotLiq < (gTot-0.0002) .and. ierLiq.eq.0)then !subtract 0.0002 to avoid randomness.
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		iPhaseLo=1
		gTot=gTotLiq
		ierCheck=ierLiq
		zFactor=zLiq
	endif
	if(LOUD.and.ierCheck)pause 'EarMixSolver: Fugi error for best estimate, initializing.'		
	dG_dxOld = ( chemPo1-chemPo2 )
	x1FracOld=x1Frac
	x1Frac=(0.5d0+49*x1Frac)/50	!make the next guess 2% closer to the middle. 
	if(ABS(x1Frac-x1FracOld) < 1d-3)x1Frac=0.499d0 !Just in case guess is 0.5d0.
	itMax=55
	do iter=1,itMax
		xFrac(1)=x1Frac
		xFrac(2)=1-xFrac(1)
		call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		ierLiq=ierFugi(1)
		if(zLiq < 0.3 .and. zLiq>0)gTotLiq=gTotLiq-0.001 !give it a nudge to avoid randomness when in the one root region. We don't use it directly anyway.
		call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zVap,ierFugi)
		chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
		gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		iPhase=0
		ierCheck=ierFugi(1)
		if(gTotLiq < (gTot-0.0002) .and. ierLiq.eq.0)then	  !avoid random switching when in one-root region.  Make it reallly liquid.
			chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
			iPhase=1
			gTot=gTotLiq
			ierCheck=ierLiq
		endif		
		if(ierCheck)then
			write(*,'(a,f7.1,2f9.4)')'T,P=',tKelvin,pMPa
			if(LOUD)pause 'EarMixSolver: Fugi error for best estimate.'
			ier=2
		endif		
		dG_dx = ( chemPo1-chemPo2 )
		fNew=dG_dx - dG_dxTarget
		fOld=dG_dxOld - dG_dxTarget
		change=fNew*(x1Frac-x1FracOld)/(fNew-fOld)
		x1FracOld=x1Frac
		dG_dxOld=dG_dx
		x1Try=x1Frac-change
		if(x1Try < 0)x1Try=x1Frac/2     !step limit to avoid x1 < 0.
		if(x1Try > 1)x1Try=(1+x1Frac)/2 !step limit to avoid x1 > 1.
		x1Frac=x1Try
		if(ABS(x1Frac-x1FracOld)/x1Frac < 1.d-5)exit
	enddo
	if(iter.ge.itMax)then
		write(*,'(a,f7.1,2f9.4)')' T,P,dG=',tKelvin,pMPa,dG_dxTarget
		if(LOUD)pause 'EarMixSolver: could not find dG_dx to match'
		ier=1
	endif
	 
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE SpinodalMixEz(NC,tKelvin,PMPa,x1FracLo,gTotLo,dG_dxLo,iPhaseLo,x1FracHi,gTotHi,dG_dxHi,iPhaseHi,iErr) 
	!Compute the points where d2G/dx1^2=0 appear at given T.
	!Method: Simple linear searches using central differences. Search x1Frac=0->1 for Lo and x1Frac=1->0 for Hi.
	!        This is crude but we just need approx spinodals for EAR and we get both from single calc of gTot(i)
	!Ref: Tang and Stephenson, AIChEJ (2011). Eubank and Hall, AICHEJ, 41:924 (1995)
	!AU: JRE, Aug 2012
	USE GlobConst !{Tc,Pc,...,Rgas, ...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION xFrac(NMX)
	dimension FUGCPURE(NMX),fugcVap(NMX),fugcLiq(NMX),ierFugi(12) !,vMolec(NMX)
	dimension gTotStor(100),dG_dxStor(100),iPhase(100)
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!data rgold,cgold/0.61803399,0.38196602/
	iErr=0
	if(iPhaseHi > 0)iErr=1
	if(NC .ne. 2)iErr=2
    if(LOUD)then
	    if(iPhaseHi > 0)pause 'SpinodalMixEz: This routine only works VLE. Sorry.'
	    if(NC.ne.2)pause 'SpinodalMixEz: This routine only works for NC=2. Sorry.'
		if(iErr > 0)return
    end if
	x1FracLo = 1e-11
	x1FracHi=1
	do i=1,2
		xFrac(i)=x1FracLo
	enddo
	nComps=NC
	do i=1,2
		xFrac(i)=1
		call Fugi(tKelvin,pMPa,xFrac,nComps,iPhaseHi,fugcVap,zFactor,ierFugi)
		call Fugi(tKelvin,pMPa,xFrac,nComps,iPhaseLo,fugcLiq,zFactor,ierFugi)
		fugcPure(i)=fugcVap(i)
		if(fugcLiq(i) < fugcVap(i) .and. ierFugi(1)==0 )fugcPure(i)=fugcLiq(i) ! min G is most stable.
		xFrac(i)=x1FracLo
	enddo
	write(*,*)' x1,mu1,mu2,dG,GV,GL,zLiq,iPhase'
	do i=1,99,1
		xFrac(1)=i/100.d0
		xFrac(2)=1-xFrac(1)
		call Fugi(tKelvin,pMPa,xFrac,nComps,iPhaseLo,fugcLiq,zLiq,ierFugi)
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		if(zLiq < 0.3 .and. zLiq>0)gTotLiq=gTotLiq-0.001 !give it a nudge to avoid randomness when in the one root region. We don't use it directly anyway.
		ierLiq=ierFugi(1)
		call Fugi(tKelvin,pMPa,xFrac,nComps,iPhaseHi,fugcVap,zFactor,ierFugi)
		chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
		gTotVap=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		gTot=gTotVap
		iPhase(i)=0
		ierCheck=ierFugi(1)
		if(gTotLiq < (gTotVap-0.0002) .and. ierLiq.eq.0)then  !Avoid randomness: make sure gTotLiq significantly less than gTotVap after nudge.
			chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
			iPhase(i)=1
			gTot=gTotLiq
			ierCheck=ierLiq
			zFactor=zLiq
		endif
		if(LOUD.and.ierCheck)pause 'SpinodalMixEz: Error from Fugi when using best gTot.'		
		gTotStor(i)=gTot
		dG_dxStor(i)=chemPo1-chemPo2
		write(*,'( f7.3,6(f10.3),i3 )')xFrac(1),chemPo1,chemPo2,dG_dxStor(i),gTotVap,gTotLiq,zFactor,iPhase(i)
	enddo
	minimax=1						 !1 for minimax indicates a natural minimum
	if( dG_dxStor(2) > dG_dxStor(1) )minimax= -1	!Multiply by -1 means maximum -> minimum
	do i=1,98
		if(minimax*dG_dxStor(i+1) > minimax*dG_dxStor(i) )exit
	enddo
	iSpinLo=i	!do loop increments one last time before ending	so check for 99
		
!	d2GminLo=1e11
!	do i=2,98
!		d2G_dx2=100*( dG_dxStor(i+1)-dG_dxStor(i-1) )/2 ! dx=( i+1-(i-1) )/100
!		if( ABS(d2G_dx2) < ABS(d2GminLo) )then
!			d2GminLo=(d2G_dx2)
!			minTypeLo=0 					  !This means the lower spinodal is a maximum
!			if( dG_dxStor(i) < dG_dxStor(2) )minTypeLo=1 !This means the lower spinodal is a minimum
!			iMinLo=i-1 !back it up one to be safe, in case it is a sharp change
!		else
!			exit !if abs(d2G) starts to increase again, then we are leaving the local min 
!		endif
!	enddo
    if(LOUD)then
	    if(iSpinLo.eq.99)pause 'SpinodalMixEz: Failed to find spinodal when increasing x1.'
    end if
	x1FracLo=DFLOAT(iSpinLo)/100
	dG_dxLo=dG_dxStor(iSpinLo)
	gTotLo=gTotStor(iSpinLo)
	iPhaseLo=iPhase(iSpinLo)
	minimax=1						 !1 for minimax indicates a natural minimum
	if( dG_dxStor(99) < dG_dxStor(98) )minimax= -1	!Multiply by -1 means maximum -> minimum
	do i=99,2,-1
		if(minimax*dG_dxStor(i-1) > minimax*dG_dxStor(i) )exit
	enddo
	iSpinHi=i !do loop increments one last time before ending
!	d2GminHi=1e11
!	do i=98,2,-1
!		d2G_dx2=100*( dG_dxStor(i+1)-dG_dxStor(i-1) )/2 ! dx=( i+1-(i-1) )/100
!		if( ABS(d2G_dx2) < ABS(d2GminHi) )then
!			d2GminHi=(d2G_dx2)
!			iMinHi=i+1 !back it up one to be safe, in case it is a sharp change
!			minTypeHi=0 					  !This means the higher spinodal is a maximum
!			if( dG_dxStor(i) < dG_dxStor(98) )minTypeHi=1 !This means the higher spinodal is a minimum
!		else
!			exit
!		endif
!	enddo
	if(iSpinHi.eq.1)then !do loop increments one last time before ending
		if(LOUD)pause 'SpinodalMixEz: Failed to find spinodal when decreasing x1.'
		iErr=2
	endif
	x1FracHi=DFLOAT(iSpinHi)/100
	dG_dxHi=dG_dxStor(iSpinHi)
	gTotHi=gTotStor(iSpinHi)
	iPhaseHi=iPhase(iSpinHi)
	!if(minTypeHi.eq.minTypeLo)pause 'SpinodalMix: minTypeHi=minTypeLo?'
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE SpinodalMix(NC,xFrac,tKelvin,PMPa,MIN,rho,zFactor,gTot,dG_dx,iErr) 
	!Compute the points where dG/dx1=0 appear at given T.
	!Method: Golden Section search
	!Ref: NumRep 
	!Ref: Tang and Stephenson, AIChEJ (2011). Eubank and Hall, AICHEJ, 41:924 (1995)
	!AU: JRE, Aug 2012
	!Notes:
	!  Assume gmol(i)=xFrac(i) => vTotCc = 1/rho
	USE GlobConst !{Tc,Pc,...,Rgas, ...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DIMENSION xFrac(NMX)
	dimension FUGCPURE(NMX),vMolec(NMX),fugcVap(NMX),fugcLiq(NMX),ierFugi(12),gTotStor(100),dG_dxStor(100)
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data rgold,cgold/0.61803399,0.38196602/
    if(LOUD)then
	    if(NC.ne.2)pause 'SpinodalMix: This routine only works for NC=2. Sorry.'
    end if
	isZiter=0 !derivatives necessary for mixes
	LIQ=1     !Use only for LLE b/c dunno how to handle root changes at const P when composition changes.
	iErr=0
	bMix=0.d0
	if (iEosOpt.Eq.5.or.iEosOpt.eq.9) then	   !AFG 2011 , Added EOS Opt. 9
		do i=1,NC
			vMolec(i)=bVolCC_mol(i)
			!xFrac(i)=gmol(i)/totMoles
			bMix=bMix+xFrac(i)*vMolec(i)
		enddo
	else
		do i=1,NC
			bMix=bMix+xFrac(i)*bVolCC_mol(i)
			!rMw(i)=rMwPlus(i)
		enddo
	endif
	!if(rho < 0)etaHi = 0.2
	minimax= -1	!routine was written to minimize, therefore take -ve to maximize.
	if(MIN)minimax=1
	x1FracLo = 1e-11
	x1FracHi=1
	do i=1,2
		xFrac(i)=x1FracLo
	enddo
	nComps=NC
	do i=1,2
		xFrac(i)=1
		call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
		call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
		fugcPure(i)=fugcVap(i)
		if(fugcLiq(i) < fugcVap(i))fugcPure(i)=fugcLiq(i) ! min G is most stable.
		xFrac(i)=x1FracLo
	enddo
	write(*,*)' x1,mu1,mu2,dG,G,iPhase'
	do i=1,99,1
		xFrac(1)=i/100.d0
		xFrac(2)=1-xFrac(1)
		call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zVap,ierFugi)
		chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
		gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
		iPhaseLo=0
		if(gTotLiq < (gTot+0.0005) .and. zLiq<0.3)then
			chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
			iPhaseLo=1
			gTot=gTotLiq
		endif		
		gTotStor(i)=gTot
		dG_dxStor(i)=chemPo1-chemPo2
		write(*,'( f7.3,4(f11.3),i4 )')xFrac(1),chemPo1,chemPo2,dG_dxStor(i),gTotStor(i),iPhaseLo
	enddo
	d2GminLo=1e11
	do i=2,98
		d2G_dx2=100*( dG_dxStor(i+1)-dG_dxStor(i-1) )/2 ! dx=( i+1-(i-1) )/100
		if( ABS(d2G_dx2) < ABS(d2GminLo) )then
			d2GminLo=(d2G_dx2)
			minTypeLo=0 					  !This means the lower spinodal is a maximum
			if( dG_dxStor(i) < dG_dxStor(2) )minTypeLo=1 !This means the lower spinodal is a minimum
			iMinLo=i-1 !back it up one to be safe, in case it is a sharp change
		else
			exit !if abs(d2G) starts to increase again, then we are leaving the local min 
		endif
	enddo
	d2GminHi=1e11
	do i=98,2,-1
		d2G_dx2=100*( dG_dxStor(i+1)-dG_dxStor(i-1) )/2 ! dx=( i+1-(i-1) )/100
		if( ABS(d2G_dx2) < ABS(d2GminHi) )then
			d2GminHi=(d2G_dx2)
			iMinHi=i+1 !back it up one to be safe, in case it is a sharp change
			minTypeHi=0 					  !This means the higher spinodal is a maximum
			if( dG_dxStor(i) < dG_dxStor(98) )minTypeHi=1 !This means the higher spinodal is a minimum
		else
			exit
		endif
    enddo
    if(LOUD)then
	    if(minTypeHi.eq.minTypeLo)pause 'SpinodalMix: minTypeHi=minTypeLo?'
    end if

	!MIN=1 means we are looking for a minimum
	if(MIN)then
		iMin=iMinLo
		if(minTypeHi)iMin=iMinHi
	else
		iMin=iMinLo	
		if(minTypeLo)iMin=iMinHi
	endif

	x1FracLoIni=x1FracLo
	x1FracHiIni=x1FracHi

	x1FracLo=(iMin-10)/100.d0;	!if rho=rhoLoIni, then no min was found and iErr=2
	x1FracHi=(iMin+10)/100.d0; !if rho=rhoHiIni, then no max was found	and iErr=1
	if(x1FracLo < 0)x1FracLo=1d-33
	if(x1FracLo > 1)x1FracHi=1-1d-11

	x1Frac1=x1FracLo+cgold*(x1FracHi-x1FracLo);
	x1Frac2=x1FracLo+rgold*(x1FracHi-x1FracLo);

	xFrac(1)=x1FracLo
	xFrac(2)=1-xFrac(1)
	call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
	chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
	gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
	chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
	gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	iPhaseLo=0
	if(gTotLiq < gTot)then
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		iPhaseLo=1
		gTot=gTotLiq
	endif		
	dG_dxLo = minimax*( chemPo1-chemPo2 )

	xFrac(1)=x1Frac1
	xFrac(2)=1-xFrac(1)
	call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
	chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
	gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
	chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
	gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	iPhase1=0
	if(gTotLiq < gTot)then
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		iPhase1=1
		gTot=gTotLiq
	endif		
	dG_dx1 = minimax*( chemPo1-chemPo2 )

	xFrac(1)=x1Frac2
	xFrac(2)=1-xFrac(1)
	call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
	chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
	gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
	chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
	gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	iPhase2=0
	if(gTotLiq < gTot)then
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		iPhase2=1
		gTot=gTotLiq
	endif		
	dG_dx2 = minimax*( chemPo1-chemPo2 )

	xFrac(1)=x1FracHi-x1FracLo
	xFrac(2)=1-xFrac(1)
	call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
	chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
	gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
	chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
	chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
	gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
	iPhaseHi=0
	if(gTotLiq < gTot)then
		chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
		chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
		iPhaseHi=1
		gTot=gTotLiq
	endif		
	dG_dxHi = minimax*( chemPo1-chemPo2 )

	itMax = 22

	do iter=1,itMax !1/2^33=1E-10; (2/3)^33=1E-6
		if(dG_dx2 < dG_dx1)then
			x1FracLo = x1Frac1;
			x1Frac1  = x1Frac2;
			x1Frac2  = rgold*x1Frac1+cgold*x1FracHi;
			dG_dxLo = dG_dx1;
			dG_dx1  = dG_dx2;
			xFrac(1)=x1Frac2
			xFrac(2)=1-xFrac(1)
			call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
			chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
			gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
			call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
			chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
			gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
			iPhase2=0
			if(gTotLiq < gTot)then
				chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
				chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
				iPhase2=1
				gTot=gTotLiq
			endif		
			dG_dx2 = minimax*( chemPo1-chemPo2 )
		else
			x1FracHi = x1Frac2;
			x1Frac2  = x1Frac1;
			x1Frac1  = rgold*x1Frac2+cgold*x1FracLo;
			dG_dxHi = dG_dx2;
			dG_dx2  = dG_dx1;
			xFrac(1)=x1Frac1
			xFrac(2)=1-xFrac(1)
			call Fugi(tKelvin,pMPa,xFrac,nComps,1,fugcLiq,zFactor,ierFugi)
			chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
			gTotLiq=xFrac(1)*chemPo1+xFrac(2)*chemPo2
			call Fugi(tKelvin,pMPa,xFrac,nComps,0,fugcVap,zFactor,ierFugi)
			chemPo1=fugcVap(1)+LOG(xFrac(1))-fugcPure(1)
			chemPo2=fugcVap(2)+LOG(xFrac(2))-fugcPure(2)
			gTot=xFrac(1)*chemPo1+xFrac(2)*chemPo2
			iPhase1=0
			if(gTotLiq < gTot)then
				chemPo1=fugcLiq(1)+LOG(xFrac(1))-fugcPure(1)
				chemPo2=fugcLiq(2)+LOG(xFrac(2))-fugcPure(2)
				iPhase1=1
				gTot=gTotLiq
			endif		
			dG_dx1 = minimax*( chemPo1-chemPo2 )
		endif
	enddo
	if(x1FracHi.eq.x1FracHiIni)iErr=1 !this means there is no max from vap search b/c T>Tc.
	if(x1FracLo.eq.x1FracLoIni)iErr=2 !this means there is no min from liq search b/c T>Tc.
	!c//that's the best we can do for now folks. get props at best guess and return.
	x1Frac = (x1Frac1+x1Frac2) /2
	!xFrac(2)=1-xFrac(1)	!unnecessary
	dG_dx=minimax*(dG_dx1+dG_dx2)/2	!minimax switches the sign, so need to switch it back
	MIN=iPhase1 ! return the phase indicator through the MIN argument
	!gTot=gTotVap
	!if(MIN)gTot=gTotLiq
	rho=pMPa/(zFactor*rGas*tKelvin) !Z=P/(rho*R*T)

	!dA_NkT=DAONKT
	!dU_NkT=DUONKT
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE FLASH(tKelvin,pMpa,nComps,iPhaseUp,zFeed,calck,itMax,xFrac,yFrac,vof,iErrCode)
	!C  FLASH ROUTINE FOR VAPOR-LIQUID OR LIQUID-LIQUID
	!C
	! iErrCode 1 = ONE OF 2-6 IS SIGNIFICANT
	!C         2 = RACRIC FOUND ALL VAPOR  FOR SOME GUESSES OF K-VALUES
	!C         3 = RACRIC FOUND ALL LIQUID FOR SOME GUESSES OF K-VALUES
	!C         4 = LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         5 = VAPOR  FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         6 = itMax EXCEEDED
	!C         7 = xFrac(iComp) < 0 for some i
	!C         8 = max[K(i)] < 1 for all i => all liquid.
	!          9 = gEx < 0 for some attempted composition and requested LLE.
	!		  11 = zVap~zLiq - warning
	!		  21 = zLiq > 0.33 - warning
	!		  22 = zVap < 0.33 - warning
	!C
	!C     REQD. ROUTINES:
	!C                      Fugi, RACRIC
	!C
	USE GlobConst
	IMPLICIT NONE
	integer iComp,nComps,iPhaseUp,itMax,maxIt,iErrCode,ierVof,iter,iFlag,iterVof
	integer ierFugi(12),iTry,iWarn,ix
	doublePrecision calck(NMX),xFrac(NMX),yFrac(NMX),zFeed(NMX)
	doublePrecision fugcUp(NMX),fugcLiq(NMX),calckNew(NMX),fugcPure(NMX)
	doublePrecision vof,tKelvin,pMpa,etaVap,etaLiq,zVap,zLiq,sumx,sumy,devMax,devi
	doublePrecision difChemPo,biggestKvalue,gEx,xInf,picardFactor,devIold
	doublePrecision rLnGam1,rLnGam2 !,denom
	doublePrecision etaPasUp,etaPasLo,zPasUp,zPasLo,zDum1,zDum2
	character*88 errMsg(11)
	common/ETA/etaLiq,etaVap,zDum1,zDum2
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	iErrCode=0
	iWarn=0
	maxIt=itMax
	iter=0
	iTry=0
	iFlag=1
	xInf=1e-11
	zVap=1 !Avoid divide zero error when failing before FUGI(zVap)
	deviOld=1000 !initialize detection for unstable successive substitution
	picardFactor=.5
	errMsg(1)='FlashSub: nonsense input.'
	errMsg(2)='FlashSub: all lower phase.'
	errMsg(3)='FlashSub: all upper phase.'
	errMsg(4)='FlashSub: fugi error for lower phase.'
	errMsg(5)='FlashSub: fugi error for upper phase.'
	errMsg(6)='FlashSub: itmax exceeded in FlashSub.'
	errMsg(7)='FlashSub: x(i)<0 for flash calc.'
	errMsg(8)='FlashSub: sumx or sumy <= 0'
	errMsg(9)='FlashSub: calck(1)=calck(2) for binary.'
	if(iPhaseUp.eq.1)then !compute gMix to be sure gEx > 0
		do iComp=1,nComps
			xFrac(iComp)=xInf
		enddo
		do iComp=1,nComps
			xFrac(iComp)=1
			CALL Fugi(tKelvin,pMpa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
			fugcPure(iComp)=fugcLiq(iComp)
			xFrac(iComp)=xInf
		enddo
	endif

	if(nComps == 2 .and. iPhaseUp==0)then ! iPhaseUp==0 means VLE
		CALL flashear(tKelvin,pMpa,nComps,iPhaseUp,zFeed,calck,itMax,xFrac,yFrac,vof,iErrCode)
		goto 86
	endif

	if(nComps == 2)then
		do ix=1,21
			xFrac(1)=(ix-1)*0.05d0
			xFrac(2)=1-xFrac(1)
			CALL Fugi(tKelvin,pMpa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
		enddo
		goto 86
	endif

	DO while(iErrCode.eq.0.and.iFlag.ne.0) !iFlag checks convergence of compositions
		iter=iter+1
		! FOR BINARIES, USE ANALYTICAL EXPRESSIONS FOR X(1) & X(2) AS GIVEN. This is handy when you don't know best feed comp to max binodal T.
		! BELOW. FOR MULTICOMPONENTS, CALL RACRIC TO GET VOF AND THEN THE X'S.
		CALL RACRIC(zFeed,calck,vof,ierVof,nComps,iterVof)
		DO iComp=1,nComps
			xFrac(iComp)=zFeed(iComp)/(1+vof*(calck(iComp)-1))
			yFrac(iComp)=calck(iComp)*xFrac(iComp)
			IF (xFrac(iComp).LT.0)then
				iErrCode=7
				goto 86
			endif
		enddo
		if(ierVof.ne.1)iTry=iTry+1
		IF(ierVof.EQ.3.and.iTry.gt.2)then
			iErrCode=2
			GOTO 86
		ELSE IF(ierVof.EQ.2.and.iTry.gt.2)then
			iErrCode=3
			GOTO 86
		END IF

		sumx=0
		sumy=0
		DO iComp=1,nComps
			sumx=sumx+xFrac(iComp)
			sumy=sumy+yFrac(iComp)
		enddo
		IF(sumx.le.0.or.sumy.le.0)then
			iErrCode=8
			if(sumx.le.0)write(*,*)'FlashSub: sumx=',sumx
			if(sumy.le.0)write(*,*)'FlashSub: sumy=',sumy
			goto 86
		endif
		DO iComp=1,nComps
			xFrac(iComp)=xFrac(iComp)/sumx
			yFrac(iComp)=yFrac(iComp)/sumy
		enddo
		CALl Fugi(tKelvin,pMpa,xFrac,nComps,1,fugcLiq,zLiq,ierFugi)
		etaPasLo=etaLiq
		zPasLo=zLiq
		IF(ierFugi(1).NE.0)then
			iErrCode=4
			GOTO 86
		END IF
		CALl Fugi(tKelvin,pMpa,yFrac,nComps,iPhaseUp,fugcUp,zVap,ierFugi)
		if(iPhaseUp.eq.1)then
			etaPasUp=etaLiq
			zPasUp=zVap
			gEx=0
			do iComp=1,nComps
				gEx=gEx+xFrac(iComp)*( fugcLiq(iComp)-fugcPure(iComp) )
			enddo
			if(gEx.le.-1111.and.iter.gt.5)then	  !temporarily disable by setting to -1111.  Some problem. JRE 1/28/04
				rLnGam1=fugcLiq(1)-fugcPure(1) 
				rLnGam2=fugcLiq(2)-fugcPure(2)
				write(*,*)'lnGam1,lnGam2',rLnGam1,rLnGam2 
				write(*,'(a,e11.4,a,7(1x,f9.8))')' gEx=',gEx,' xFrac=',(xFrac(iComp),iComp=1,nComps)
				write(*,*)'Flash warning: LLE request when gEx < 0'
				iWarn=23
				!goto 86
			endif
		endif
		IF(ierFugi(1).NE.0)then
			iErrCode=5
			GOTO 86
		END IF
		iFlag=0
		devMax=0
		DO iComp=1,nComps
			difChemPo=fugcLiq(iComp)-fugcUp(iComp)
			if(difChemPo.gt.708)difChemPo=708 !exp overflow otherwise
			if(tc(iComp).gt.999.and.difChemPo.gt.-11.and.iPhaseUp.eq.0)difChemPo=-11
			calckNew(iComp)=EXP(difChemPo)
			if(calckNew(iComp).LT.1.E-33)calckNew(iComp)=1.0E-33
			DEVI=(calckNew(iComp)-calck(iComp))/calckNew(iComp)
			IF(ABS(DEVI).GT.1.E-8)iFlag=1
			IF(ABS(DEVI).GT.devMax)devMax=ABS(DEVI)
			if(calckNew(iComp).gt.biggestKvalue)biggestKvalue=calckNew(iComp)
		enddo
		biggestKvalue=1e-33
		if(ABS(DEVI).gt.ABS(devIold))then
			if(devi*deviOld.gt.0)then
				picardFactor= -picardFactor/10
			else
				picardFactor= picardFactor/10
			endif
			!if(ABS(devi).lt.1e-4)picardFactor=sign(min(.01,ABS(devi*100)),devi)
		endif
		DO iComp=1,nComps
			calck(iComp)=picardFactor*calckNew(iComp)+(1-picardFactor)*calck(iComp)
			if(calck(iComp).gt.biggestKvalue)biggestKvalue=calck(iComp)
		enddo
		!if(biggestKvalue.le.1)iErrCode=8
		if(iter.gt.maxIt)iErrCode=6
		deviOld=devi
		picardFactor=picardFactor*1.5
		if(picardFactor.gt.1)picardFactor=1
		!if(ABS(DEVI).lt.ABS(devIold))devIold=devi !trap unstable successive substitution
	enddo !While(iErrCode.eq.0)
86  continue
	itMax=iter
    if(LOUD)then
	    if(ABS(zVap).lt.1e-11)pause 'FLASH error: zVap=0.'
    end if
	IF(DABS((zVap-zLiq)/zVap).LT.0.01.and.iPhaseUp.eq.0)iErrCode=11
	!IF (zLiq.GT.0.33)iErrCode=21	 !Note: for polymers zLiq >> 1 because vLiq is so large and P ~ 0.1MPa >> pSat
	if(iErrCode.ne.0)then
		write(*,*)errMsg(iErrCode)
		if(LOUD)pause
	endif
	if(iWarn.ne.0)iErrCode=iErrCode+iWarn
	IF (zVap.LT.0.3.and.iPhaseUp.eq.0 )iErrCode=22
	RETURN
	END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! iErrCode 1 = ONE OF 2-6 IS SIGNIFICANT
	!C         2 = RACRIC FOUND ALL VAPOR  FOR SOME GUESSES OF K-VALUES
	!C         3 = RACRIC FOUND ALL LIQUID FOR SOME GUESSES OF K-VALUES
	!C         4 = LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         5 = UPPER  FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         6 = itMax EXCEEDED
	!C         7 = xFrac(iComp) < 0 for some i
	!C         8 = max[K(i)] < 1 for all i => all liquid.
	!          9 = gEx < 0 for some attempted composition and requested LLE.
	!		  11 = K1~K2 - warning
	!		  21 = zLiq > 0.33 - warning
	!		  22 = zVap < 0.33 - warning
	SUBROUTINE LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,XU,UOF,iErrCode)
	USE GlobConst  !nmx
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DoublePrecision fugcPure(NMX),fugcLo(NMX),fugcUp(NMX),Gmix(21),dG_dx1(21),d2G_dx12(21),xFrac(NMX),actCo(NMX),gamInf(NMX)
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	DIMENSION X(nmx),XU(nmx),ZFEED(nmx),calck(nmx)
	integer ierFugi(12)
	LOGICAL LOUDER
	Character*1 answer
	!COMMON/FEED/ZFEED
	iErrCode=0
	LOUDER=LOUD
	!LOUDER=.TRUE.
	iPhaseUp=1 !liquid-liquid

	if(NC == 2)then	 !compute gMix to be sure gEx > 0
		if(initialize==2)then !check d2G/dx1^2 < 0 somewhere and spinodal,maxGams for initial guess.
			xFrac(1:NC)=zeroTol  ! initialize.
			do iComp=1,NC
				xFrac(iComp)=1
				CALL Fugi(T,P,xFrac,NC,1,fugcLo,zLo,ierFugi)
				fugcPure(iComp)=fugcLo(iComp)
				if(iComp==1)gamInf(2)=fugcLo(2)
				if(iComp==2)gamInf(1)=fugcLo(1)
				xFrac(iComp)=zeroTol
			enddo
			gamInf(1:NC)=exp( gamInf(1:NC)-fugcPure(1:NC) )
			Gmix=0
			gam1max=gamInf(1)
			x1Max1=0
			gam2max=gamInf(2)
			x1Max2=1
			if(LOUDER)print*,'    x1        G/RT        GE/RT        gam1       gam2'
			if(LOUDER)write(*,'(3F10.3,2F10.3)')0.,0.,0.,gamInf(1),1.0
			do ix=1,19
				xFrac(1)=ix*0.05d0
				xFrac(2)=1-xFrac(1)
				CALL Fugi(T,P,xFrac,NC,1,fugcLo,zLo,ierFugi)
				actCo(1:NC)=exp( fugcLo(1:NC)-fugcPure(1:NC) )
				Gmix(ix)=xFrac(1)*DLOG(xFrac(1)*actCo(1))+xFrac(2)*DLOG(xFrac(2)*actCo(2))
				GE_RT=xFrac(1)*DLOG(actCo(1))+xFrac(2)*DLOG(actCo(2))
				if(LOUDER)write(*,'(6F10.3)')xFrac(1),Gmix(ix),GE_RT,actCo(1),actCo(2)
				if(actCo(1) > gam1max)x1Max1=xFrac(1)
				if(actCo(1) > gam1max)gam1max=actCo(1)
				if(actCo(2) > gam2max)x1Max2=xFrac(2)
				if(actCo(2) > gam2max)gam2max=actCo(2)
			enddo
			!if(LOUDER)print*,'gam2max,x1Max2=',gam2max,x1Max2
			if(LOUDER)write(*,'(3F10.6,2F10.3)')1.,0.,0.,1.,gamInf(2)
			d2Gmin=1234
			!if(LOUDER)write(*,'(6F10.6)')x1,Gmix(ix),dG_dx1,d2G_dx12(ix)
			if(LOUDER)print*,'    x1        Gmix       dG_dx1      d2G_dx1^2'
			spinLo=0
			spinHi=1
			do ix=1,19
				x1=ix*0.05d0
				d2G_dx12(ix)=( Gmix(ix-1)-2*Gmix(ix)+Gmix(ix+1) )/0.0025d0
				if(	d2G_dx12(ix)*d2G_dx12(ix-1) < 0 .and. spinLo==0)spinLo=x1
				if(	d2G_dx12(ix)*d2G_dx12(ix-1) < 0 .and. spinLo >0)spinHi=x1
				dG_dx1(ix)=( Gmix(ix-1)-Gmix(ix+1) )/0.1d0
				if(d2G_dx12(ix) < d2Gmin)then
					d2Gmin=d2G_dx12(ix)
					x1min=x1
				endif
				if(LOUDER)write(*,'(6F10.6)')x1,Gmix(ix),dG_dx1(ix),d2G_dx12(ix)
			enddo
			if(LOUDER)print*,'x1min,d2Gmin,x1Max1,gam1Max,x1Max2,gam2Max'
			if(LOUDER)write(*,'(11F10.4)')x1min,d2Gmin,x1Max1,gam1Max,x1Max2,gam2Max
			if(d2Gmin > 0)then
				if(LOUDER)write(*,*)'LLEFL: d2G/dx1^2 > 0 for all x1. Min(d2G)=',d2Gmin
				if(LOUDER)pause 'No LLE for all x1=[0,1]'
				iErrCode=2
				return
			endif
			x(1)=x1Max1
			xU(1)=x1Max2
			answer='y'
			if(LOUDER)write(*,'(a,2f4.2)')' spinLo,spinHi=',spinLo,spinHi
			do while(answer/='n' .and. LOUDER) ! only enable manual looping when LOUDER.
				print*,'Enter guess for x1Lo,x1Up'
				read*,x(1),xU(1)
				x(2)=1-x(1)
				xU(2)=1-xU(1)
				CALL Fugi(T,P,x ,NC,1,fugcLo,zLo,ierFugi)
				if(ierFugi(1) /= 0 .and. LOUDER)iErrCode=4
				actCo1Lo=EXP(fugcLo(1)-fugcPure(1))
				actCo2Lo=EXP(fugcLo(2)-fugcPure(2))
				GmixLo=x(1)*DLOG(x(1)*actCo(1))+x(2)*DLOG(x(2)*actCo(2))
				CALL Fugi(T,P,xU,NC,1,fugcUp,zUp,ierFugi)
				if(ierFugi(1) /= 0 .and. LOUDER)iErrCode=5
				if(iErrCode > 0)exit
				actCo(1)=EXP(fugcUp(1)-fugcPure(1))
				actCo(2)=EXP(fugcUp(2)-fugcPure(2))
				GmixLo=xU(1)*DLOG(xU(1)*actCo(1))+xU(2)*DLOG(xU(2)*actCo(2))
				calcK(1:NC)=EXP(fugcLo(1:NC)-fugcUp(1:NC))
				!if(LOUDER)print*,'Fugi done. K1,K2=',calcK1,calcK2
				if( ABS(calcK(1)-calcK(2)) < zeroTol)pause 'LLEFL: (K1-K2) ~ 0 in denom.'
				x1Lo=(1-calcK(2))/(calcK(1)-calcK(2))
				x1Up=x1Lo*calcK(1)
				if(LOUDER)write(*,'(2f10.6,2E12.4,2F10.6)')x(1),xU(1),calcK(1),calcK(2),x1Lo,x1Up
				print*,'Try manually again? (y/n)'
				read*,answer
			enddo
		elseif(initialize==1)then
				x(1)=(1-calcK(2))/(calcK(1)-calcK(2))
				xU(1)=x(1)*calcK(1)
		endif !initialize==1

		!itMax=21 : itMax is a calling argument
		alpha=0.3  ! enable Picard iteration.
		if(LOUDER)print*,'Enter alpha for picard iteration'
		if(LOUDER)read*,alpha
		x1Lo=x(1)
		!x1Up=spinHi*alpha+(1-alpha)*1
		!x(1)=x1Lo
		x1up=xU(1)
		if(LOUDER)print*,'    x1       x1U      calcK1      calcK2      x1new    x1Upnew'
		do iter=1,itMax
			x(1)=alpha*x1Lo+(1-alpha)*x(1)
			x(2)=1-x(1)
			xU(1)=x1Up*alpha+(1-alpha)*xU(1)
			xU(2)=1-xU(1)
			!if(LOUDER)print*,'Calling Fugi. x1Lo,x1Up=',x(1),xU(1)
			CALL Fugi(T,P,x ,NC,1,fugcLo,zLo,ierFugi)
			if(ierFugi(1) /= 0 .and. LOUDER)iErrCode=4
			CALL Fugi(T,P,xU,NC,1,fugcUp,zUp,ierFugi)
			if(ierFugi(1) /= 0 .and. LOUDER)iErrCode=5
			if(iErrCode > 0)exit
			calcK(1:NC)=EXP(fugcLo(1:NC)-fugcUp(1:NC))
			!if(LOUDER)print*,'Fugi done. K1,K2=',calcK1,calcK2
			if( ABS(calcK(1)-calcK(2)) < zeroTol)then
				if(LOUDER)pause 'LLEFL: (K1-K2) ~ 0 in denom.'
				iErrCode=11
				return
			endif
			x1Lo=(1-calcK(2))/(calcK(1)-calcK(2))
			x1Up=x1Lo*calcK(1)
			if(LOUDER)write(*,'(2f10.6,2E12.4,2F10.6)')x(1),xU(1),calcK(1),calcK(2),x1Lo,x1Up
			if(x1Lo < zeroTol .or. x1Lo > 1)then
				if(LOUDER)pause 'LLEFL: x1Lo out of range.'
				iErrCode=7
				return
			endif
			if(x1Up < zeroTol .or. x1Up > 1)then
				if(LOUDER)pause 'LLEFL: x1Up out of range.'
				iErrCode=7
				return
			endif
			dev=( 0.5d0-ABS(0.5-x(1)) )/( 0.5d0-ABS(0.5-x1Lo) )	! %dev in xSolute where solute is whichever has smallest composition.
			if( ABS(DLOG(dev+zeroTol)) < 0.0008d0)exit ! 3 sig figs on xSolute ~ close enough!
		enddo !iter
		if(iter > itMax)iErrCode=6
		itMax=iter
		zFeed(1:NC)=x(1:NC)/2+xU(1:NC)
		UOF=0.5
	else ! ie. NC > 2.
		IF(Initialize==0)iErrCode=1
		call FLASH(T,P,NC,iPhaseUp,zFeed,CALCK,ITMAX,X,XU,UOF,iErrCode)
	endif
861	continue
	return
	end	 !LLEFLASH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! iErrCode 1 = NC =/= 2
	!C         2 = RACRIC FOUND ALL VAPOR  FOR SOME GUESSES OF K-VALUES
	!C         3 = RACRIC FOUND ALL LIQUID FOR SOME GUESSES OF K-VALUES
	!C         4 = LIQUID FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         5 = UPPER  FUGACITY COEFFICIENT CALCULATION HAD ERROR
	!C         6 = itMax EXCEEDED
	!C         7 = xFrac(iComp) < 0 for some i
	!C         8 = max[K(i)] < 1 for all i => all liquid.
	!          9 = gEx < 0 for some attempted composition and requested LLE.
	!		  11 = K1~K2 - warning
	!		  21 = zLiq > 0.33 - warning
	!		  22 = zVap < 0.33 - warning
	SUBROUTINE SLEFL(T,P,NC,Initialize,ITMAX,iSolute,XidSoln,X,iErrCode)
	USE GlobConst  !nmx
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DoublePrecision fugcPure(NMX),fugcLo(NMX),gamInf(NMX),xFrac(NMX) !,fugcUp(NMX),Gmix(21),dG_dx1(21),d2G_dx12(21),actCo(NMX)
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	DIMENSION X(nmx),XidSoln(nmx) !,ZFEED(nmx),calck(nmx)
	integer ierFugi(12)
	LOGICAL LOUDER
	iErrCode=0
	LOUDER=LOUD
	!LOUDER=.TRUE.

	if(NC /= 2)then	 
		pause 'SLEFL: Sorry, this option only works for NC=2 at this time.'
		iErrCode=1
		return
	endif
	xFrac(1:NC)=zeroTol  ! initialize.
	do iComp=1,NC  !compute FugcPure, gamInf.
		xFrac(iComp)=1
		CALL Fugi(T,P,xFrac,NC,1,fugcLo,zLo,ierFugi)
		fugcPure(iComp)=fugcLo(iComp)
		if(iComp==1)gamInf(2)=fugcLo(2)
		if(iComp==2)gamInf(1)=fugcLo(1)
		xFrac(iComp)=zeroTol
	enddo
	gamInf(1:NC)=exp( gamInf(1:NC)-fugcPure(1:NC) )
	xInf=XidSoln(iSolute)/gamInf(iSolute)
	noSoln=1 ! Check if solution exists.
	if( xInf < 1 )noSoln=0 ! sufficient condition. Not necessary, e.g. if gamInf << 1 and XidSoln(iSolute) < 1
	if( XidSoln(iSolute) < 1)noSoln=0 ! no guarantee, but it's worth a try.
	if(noSoln==1)then
		if(LOUDER)print*,'Xid,gamInf=',XidSoln(iSolute),gamInf(iSolute)
		if(LOUDER)pause 'SLEFL: Sorry, no solution exists. Try lower temperature.'
		iErrCode=2
		return		
	endif
	alpha=0.5  ! enable Picard iteration.
	if(initialize==3)then
		X(iSolute)=xIdSoln(iSolute)*(1-alpha)+(alpha)*xInf
	elseif(initialize==2)then
		X(iSolute)=xIdSoln(iSolute)
	elseif(initialize==1)then
		X(iSolute)=xInf
	endif ! initialize. Otherwise, Xinitial=X(input).
	iSolvent=1
	if(iSolute==1)iSolvent=2
	if(initialize > 0)X(iSolvent)=1-X(iSolute) 

	!itMax=21 : itMax is a calling argument
	if(LOUDER)print*,'Enter alpha for picard iteration'
	if(LOUDER)read*,alpha
	if(LOUDER)write(*,'(a,6E12.4)')'Xid=',XidSoln(iSolute)
	if(LOUDER)print*,'    x(iS)       gamma      x1new    '
	x1Lo=X(iSolute)
	do iter=1,itMax
		x(iSolute)=alpha*x1Lo+(1-alpha)*x(iSolute)
		x(iSolvent)=1-x(iSolute)
		!if(LOUDER)print*,'Calling Fugi. x1Lo,x1Up=',x(1),xU(1)
		CALL Fugi(T,P,x ,NC,1,fugcLo,zLo,ierFugi)
		bVolMix=SUM( x(1:NC)*bVolCC_mol(1:NC) )
		rhoMol_cc=P/(zLo*Rgas*T)
		eta=bVolMix*rhoMol_cc 
		if(eta < 0.15)print*,'SLEFL: eta < 0.5? = ',eta 
		if(ierFugi(1) /= 0 .and. LOUDER)iErrCode=4
		if(iErrCode > 0)exit
		gamma=EXP(fugcLo(iSolute)-fugcPure(iSolute))
		!if(LOUDER)print*,'Fugi done. K1,K2=',calcK1,calcK2
		x1Lo=XidSoln(iSolute)/gamma
		if(LOUDER)write(*,'(f10.6,2E12.4,2F10.6)')x(1),gamma,x1Lo
		if(x1Lo < zeroTol .or. x1Lo > 1)then
			if(LOUDER)pause 'SLEFL: x1Lo out of range.'
			iErrCode=7
			return
		endif
		dev=( x(iSolute) /x1Lo )	! %dev in xSolute where solute is whichever has smallest composition.
		if( ABS(DLOG(dev)) < 0.00008d0)exit ! 4 sig figs on xSolute ~ close enough!
	enddo !iter
	if(iter > itMax)iErrCode=6
	itMax=iter
861	continue
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LLVE(T,zFeed,NC,INIT,P,ITMAX,Y,X,X2,IER,VOF,S)
	USE GlobConst !{Tc,Pc,...,Rgas, ...
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(nmx),Y(nmx),IER(6),X2(nmx),S(nmx),zFeed(nmx)
	ITMAXO=ITMAX
	POLD=P
	CALL PCAL(T,zFeed,NC,INIT,P,ITMAXO,Y,X,X2,IER,VOF,S)
	iter=0
	do while ( ABS((P-POLD)/POLD).GT.1.0E-4.and.iErrCode.eq.0)
		iter=iter+1
		if(iter.gt.25)iErrCode=1
		POLD=P                    
		ITMAXO=ITMAX
		CALL PCAL(T,zFeed,NC,INIT,P,ITMAXO,Y,X,X2,IER,VOF,S)
		if(ier(1).ne.0)iErrCode=2
	enddo
	ITMAX=ITMAXO
	RETURN
	END
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE PCAL(T,zFeed,NC,INIT,P,ITMAX,Y,X,X2,IER,VOF,S)
	USE GlobConst !nmx
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION X(nmx),Y(nmx),IER(6),IERF(12),X2(nmx),CALCK(nmx),S(nmx),zFeed(nmx)
	DO I=1,6
		IER(I)=0
	enddo
	DO I=1,12
		IERF(I)=0
	enddo
	INITF=1
	do iComp=1,NC
		calck(iComp)=s(iComp)
	enddo
	CALL LLEFL(T,P,NC,INITF,CALCK,ITMAX,zFeed,X,X2,VOF,iErrCode)
	ITMAXP=1111
	IF (iErrCode .ne. 0)WRITE(6,*)' LLEFL',iErrCode
	CALL BUBPL(T,X2,NC,INIT,P,ITMAXP,Y,IERF)
	DO I=1,6
		IF (IERF(I).NE.0)IFLAG=2
	enddo
    IF (IFLAG.EQ. 2)WRITE(6,*)' BUBPL',(IERF(I),I=1,6)
	DO J=1,NC
		S(J)=CALCK(J)
	enddo
	ier(1)=iErrCode
	RETURN
	END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE RACRIC(zFeed,kRatio,vof,markd,nComps,iter)
	!IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	implicit none
	!PURPOSE:  SOLVE MULTICOMPONENT FLASH EQUATION BASED ON RACHFORD-RICE ALGORITHM
	!	zFeed	FEED COMPOSITION
	!	kRatio	K-VALUE DISTRIBUTION COEFFICIENTS
	!	vof 	VAPOR/FEED RATIO
	!	MARKD	= 1 (IERVOF=MARKD) MEANS NO ERROR
	!			= 2 all liquid
	!			= 3 all vapor
	!			= 4 iteration limit exceeded
	!	
	integer it,iter,markd,iComp,nComps;
	doublePrecision kRatio(nComps),zFeed(nComps);
	doublePrecision	vof,f0,f1,change,fd,f1d,f2d,c1;
	it=0 ;
	iter=0;
	markd=1;
	f0=0.0;
	f1=0.0;
	do iComp=1,nComps;
		f0=f0+(kRatio(iComp)-1.0)*zfeed(iComp);
		f1=f1+(kRatio(iComp)-1.0)*zfeed(iComp)/(kRatio(iComp)+1E-33);
	enddo;
	if(f0.le.0)then
		vof=0.0;
		markd=2;
		return;
	else if(f1.ge.0)then
		vof=1.0
		markd=3
		return
	endif

	change=1111
	vof=0.5
	do while(markd.eq.1.and.ABS(change).gt.0.001)
		iter=iter+1
		fd=0.0
		f1d=0.0
		f2d=0.0
		do iComp=1,nComps
			c1=(kRatio(iComp)-1)/( (kRatio(iComp)-1)*vof+1 )
			fd=fd+c1*zfeed(iComp)
			f1d=f1d-c1*c1*zfeed(iComp)
			f2d=f2d+2.0*c1*c1*c1*zfeed(iComp)
		enddo
		vof=vof-2.0*fd*f1d/(2.0*f1d*f1d-fd*f2d)
		if (vof.LE.0)then
			vof=0.0001
			it=it+1
			if(it.gt.2)markd=2
		elseif(vof.GE.1.0)then
			vof=0.9999
			it=it+1
			if(it.gt.2)markd=3
		endif
		change=fd/vof/f1d
		if(iter.gt.1111)markd=4
	enddo
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE VLEFL (T,P,NC,INIT,zFeed,CALCK,ITMAX,X,y,VOF,IER)
	USE GlobConst !Tc,Pc, ...
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	DIMENSION CALCK(*),X(*),Y(*),ZFEED(*)
	IF(INIT.EQ.0)THEN
		DO I=1,NC
			VP=PC(I)*EXP( 5.3727*(1+ACEN(I))*(1-TC(I)/T) );
			CALCK(I)=VP/P;
			if(CALCK(I).LT.1.E-33)CALCK(I)=1.0E-33;
		enddo
	ENDIF
	iPhaseUp=0 !vapor-liquid
	call Flash(T,P,NC,iPhaseUp,zFeed,CALCK,ITMAX,X,y,VOF,IER)
	return
	end

