	Subroutine BipIo(NC) !note: results passed through common/bips
	USE GlobConst
	USE BIPs
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	!Integer idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	!DoublePrecision aBipAD(maxTypes,maxTypes),aBipDA(maxTypes,maxTypes) !association bips
	implicit DoublePrecision(a-h,k,o-z)
	character answer*1
	DO J=2,NC
		DO I=1,J-1
			WRITE(*,*)'ENTER KIJ,KTIJ FOR THE PAIR',I,J,'    '
			READ(*,*)KIJ(I,J),KTIJ(I,J)
			WRITE(52,*)' KIJ,KTIJ FOR THE PAIR',I,J,'    '
			write(52,*)KIJ(I,J),KTIJ(I,J)
			KIJ(J,I)=KIJ(I,J)
			KTIJ(J,I)=KTIJ(I,J)
			if(iEosOpt.EQ.3 .or. iEosOpt.EQ.7)THEN
				WRITE(*,602)i,j,j,i
				READ(*,*)xsTau(I,J),xsTau(J,I)
				WRITE(52,602)i,j,j,i
				write(52,*)xsTau(I,J),xsTau(J,I)

				xsAlpha(I,J)=0.3d0
				if(ABS(xsTau(I,J)).ge.0.1.or.ABS(xsTau(J,I)).ge.0.1)then
					write(*,*)'Accept xsNrtlAlpha=0.3? (y/n) '
					read(*,'(a)')answer
					write(52,*)'Accept xsNrtlAlpha=0.3? (y/n) '
					write(52,'(a)')answer
					if(answer.ne.'y'.and.answer.ne.'Y')then
						write(*,*)'Enter xsNrtlAlpha'
						read(*,*)xsAlpha(I,J)
						write(52,*)'Enter xsNrtlAlpha'
						write(52,*)xsAlpha(I,J)
					endif
				endif
				xsAlpha(J,I)=xsAlpha(I,J)
			endif !iEosOpt.eq.3
		enddo
		if(LOUD)write(*,601)J,ID(J),(KIJ(J,I),I=1,J-1)
	enddo !J=2,NC

	if(iEosOpt.ge.4.and.itsEven(iEosOpt))then
		DO J=2,NC
			DO I=1,J-1
				WRITE(*,*)'ENTER HIJ,HTIJ FOR THE PAIR',I,J,'    '
				READ(*,*)HIJ(I,J),HTIJ(I,J)
				WRITE(52,*)'ENTER HIJ,HTIJ FOR THE PAIR',I,J,'    '
				write(52,*)HIJ(I,J),HTIJ(I,J)
				HIJ(J,I)=HIJ(I,J)
				HTIJ(J,I)=HTIJ(I,J)
			enddo
			if(LOUD)write(*,601)J,ID(J),(HIJ(J,I),I=1,J-1)
		enddo
	endif !iEosOpt.ge.4

	if(iEosOpt==5)then
		write(*,'(5x,11i8)')(idLocalType(j),j=1,nTypesTot)
		do i=1,nTypesTot
			write(*,'(i5,11f8.4)')idLocalType(i),(aBipDA(i,j),j=1,nTypesTot)
		enddo
		print*,'Enter i,j,BipDA(i,j)'
		read*,i,j,aBipDA(i,j)
		aBipAD(j,i)=aBipDA(i,j)
		aBipAD(i,j)=aBipDA(i,j)	!keep symmetric for now. JRE 20220108
		aBipDA(j,i)=aBipDA(i,j)
	endif

601	FORMAT(I3,I5,9(1X,F7.4))
602	format(' ENTER xsTau(',i2,',',i2,'),xsTau(',i2,',',i2,') OR ZEROES TO NULLIFY XS EFFECT')
	return
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ErrCheck(ier,iflag)
	DIMENSION ier(12)
      IFLAG=0
	DO I=1,11
		IF (ier(I).NE.0)IFLAG=1
	enddo
	IF (IFLAG.EQ.1)WRITE(*,*)'ERROR !!- CHECK ANSWERS CAREFULLY'
	IF(ier(2).EQ.1)WRITE(*,*)'ERROR ALL UPPER PHASE'
	IF(ier(3).EQ.1)WRITE(*,*)'ERROR ALL LOWER PHASE'
	IF(ier(4).EQ.4)WRITE(*,*)'ERROR IN FUGI-NEG LOG CALCD'
	IF(ier(4).EQ.5)WRITE(*,*)'ERROR IN FUGI-FUGACITY OVERFLOWS'
	IF(ier(4).EQ.6)WRITE(*,*)'ERROR IN FUGI-Z ITERATION NO CNVRG'
	IF(ier(5).EQ.1)WRITE(*,*)'ERROR VAPOR AND LIQUID ROOTS CLOSE'
	IF(ier(6).EQ.1)WRITE(*,*)'ERROR VLE ITERATION NO CNVRG'
	IF(ier(7).EQ.1)WRITE(*,*)'ERROR VLE ITERATION FAILED TO IMPROVE'
	IF(ier(8).EQ.1)WRITE(*,*)'ERROR IN TC,PC,OR X,Y'
	IF(ier(9).EQ.1)WRITE(*,*)'ERROR P SPECIFIED < 0'
	IF(ier(10).EQ.1)WRITE(*,*)'ERROR T SPECIFIED IS UNREASONABLE'
	IF(ier(11).EQ.1)WRITE(*,*)'ERROR MORE THAN 10 COMPONENTS OR ITMAX<1'
	IF(ier(12).EQ.1)WRITE(*,*)'FUGACITIES NOT EQUAL ON LAST ITER'
	if(iflag.ne.0)write(*,*)'Error on last iteration for next pt.'
	if(iflag.ne.0)ier(1)=1
	RETURN
	END
	    
    !***********************************************************************
      
	SUBROUTINE BifurcationDiagram(NC)
	!C
	!C  PURPOSE:  COMPUTE TX1-x2 for LLE GIVEN P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character outFile*251
	!DIMENSION ZFEED(NMX),kPas(NMX)
	DIMENSION X(NMX),FUGC(NMX),ier(12)
	!common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!data nx,x1/13,0.0001,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.9999/
    if(LOUD)then
	    if(nc > 1)pause 'error - this option is only for 1 compo'
    end if
	if(nc > 1)return
601 format(i4,f8.2,e12.4,i4,a)
	x(1) = 1
	write(*,*)'Enter P(MPa) for V vs T plot.'
	read(*,*)PMPa
	outFile=TRIM(masterDir)//'\output\Bifurcation.txt'
	open(61,file=outFile)
	vOld1=0
	vOld2=0
	LIQ=2
	do it=130,1,-1
		iBifurcate=0
		tKelvin=it*TC(1)/100
		CALL FUGI(tKelvin,PMPa,X,NC,LIQ,FUGC,zFactor,ier)
		vNew=rGas*tKelvin*zFactor/PMPa
		if(it==1)then
			vOld2=vNew
			cycle
		elseif(it==2)then
			vOld1=vNew
			cycle
		else
			vNext=2*vOld1-vOld2
			if(ABS(vNew-vNext)/vNew > 0.01)iBifurcate=1
		endif
		if(ier(1).ne.0)cycle
		if(iBifurcate)then
			write(61,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
			write(* ,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
		else
			write(61,601)it,tKelvin,vNew,ier(1)
			write(* ,601)it,tKelvin,vNew,ier(1)
		endif
		vOld2=vOld1
		vOld1=vNew
	enddo
	if(LOUD)write(*,*)'Vapor done. Starting Liquid.'
	vOld1=0
	vOld2=0
	LIQ=3
	do it=1,130 
		iBifurcate=0
		tKelvin=it*TC(1)/100
		CALL FUGI(tKelvin,PMPa,X,NC,LIQ,FUGC,zFactor,ier)
		vNew=rGas*tKelvin*zFactor/PMPa
		if(it==1)then
			vOld2=vNew
			cycle
		elseif(it==2)then
			vOld1=vNew
			cycle
		else
			vNext=2*vOld1-vOld2
			if(ABS(vNew-vNext)/vNew > 0.01)iBifurcate=1
		endif
		if(ier(1).ne.0)cycle
		if(iBifurcate)then
			write(61,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
			write(* ,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
		else
			write(61,601)it,tKelvin,vNew,ier(1)
			write(* ,601)it,tKelvin,vNew,ier(1)
		endif
		vOld2=vOld1
		vOld1=vNew
	enddo
	close(61)
	if(LOUD)pause 'Success! Output in output\Bifurcation.txt'
	return
	end

!***********************************************************************
      
	SUBROUTINE Binodal(NC)
	!C
	!C  PURPOSE:  COMPUTE TX1-x2 for LLE GIVEN P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*251 outFile
	LOGICAL LOUDER
	DIMENSION ZFEED(NMX),kPas(NMX)
	DIMENSION X(NMX),xU(NMX)
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!data nx,x1/13,0.0001,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.9999/
	LOUDER=LOUD
	LOUDER=.TRUE.
    if(LOUDER)then
	    if(nc.gt.2)pause 'error - this option is only for 2 compos'
    end if
	if(nc.gt.2)return
	!if(LOUD)write(*,*)'enter file NAME for output'
	!read(*,'(A1)')outfile
	outFile=TRIM(masterDir)//'\output\Binodal.txt'
	open(61,file=outFile)
600	format('   T(K)',4x,'  P(MPA)',5x,'x1Alpha',4x,'x1Beta',2x,'etaAlpha',3x, 'etaBeta')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),6(F9.5,1X))
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')           
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))      
	WRITE(6,*)'ENTER PRESSURE(MPa) and starting Temperature (K)'
	READ(5,*)P,T
	if(LOUD)write(*,600)
	write(61,600)
	MAXIT=1111
	zero=0

	WRITE(*,*)' '
	WRITE(*,*)'ENTER FEED COMPOSITION      '
	READ(5,*)(ZFEED(I),I=1,NC)
	!c  zFeed is passed through common
	kPas(1)=1111
	kPas(2)=.001
	init=1
	do iTempK=1,1000
		ITMAX=MAXIT
		CALL LLEFL(T,P,NC,INIT,kPas,ITMAX,zFeed,X,xU,UOF,iErrCode)
		if(iErrCode.ne.0)then
			call PrintErrFlash(iErrCode)
		else
			if(LOUDER)write(* ,601)T,P,X(1),xU(1),etaPasLo,etaPasUp
			write(61,'(1X,F7.2,2X,F10.6,1X,4(f13.11,1X))')T,P,X(1),xU(1),etaPasLo,etaPasUp
		endif
		if(iErrCode.ne.0.or.kPas(1).lt.1.05)goto 86

		tIncrease=5
		if(kPas(1).lt.2)tIncrease=1
		if(kPas(1).lt.1.3)tIncrease=0.5
		T=T+tIncrease
	enddo

86	continue
	write(*,*)'These results are tabulated in',TRIM(outFile)
	if(LOUDER)pause 'Success! Check your results.'
	CLOSE(61)   
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE BpFromFile(NC)
!C
!C  PURPOSE:  BP CALCULATIONS FOR A SERIES
!C	OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
!C	WORKS LIKE KI OPTION BUT ON FILE INSTEAD OF SINGLE POINT.
!C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
!C  PROGRAMMED BY:  JRE 2/96
!C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
!C  REFERENCE:NUMERICAL RECIPES
!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER FN*44
	PARAMETER (RGOLD=.61803399,CGOLD=.38196602)

	dimension xFrac(NMX),yFrac(NMX),ier(12)
	dimension idFile(NMX)
	!common/tpxyData/tDat(111),PDAT(111),XDAT(111),npts
	COMMON/eta/etaL,etaV,ZL,ZV

	if(LOUD)write(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (..\VleData\FN.FT)'
	if(LOUD)write(*,*)'FIRST LINE =NPTS,nComps,(id(i),i=1,nComps)'
	if(LOUD)write(*,*)'2ND LINE = TK, PMPA, X1,x2,x3,...'
	READ(*,'(A)')FN
	if(DEBUG)then
		FN='C:\spead\calceos\VleData\'//TRIM(FN)
	else
		FN=TRIM(masterDir)//'\VleData\'//TRIM(FN)
	endif
	OPEN(51,FILE=FN)
	READ(51,*,ERR=861)NPTS,nComps,(idFile(iComp),iComp=1,nComps)
	iCheck=1
	if(nComps.ne.NC)iCheck=0
	do iComp=1,nComps
		if(idFile(iComp).ne.id(iComp))iCheck=0
	enddo
	if(iCheck.eq.0)then
		if(LOUD)write(*,*)'Error 1st line of file: check nComps and id order.'
	endif
	MAXIT=1111
	INIT=1
	iFlag=0
	open(66,file='BpFromFile.txt')
      DO I=1,NPTS
		READ(51,*,ERR=861)tKelvin,pMPa,(xFrac(iComp),iComp=1,nComps)
		ITMAX=MAXIT
		pCalc=pMPa
		calL BUBPL(tKelvin,xFrac,nComps,init,pCalc,itMax,yFrac,ier)
		if(ier(1).ne.0)then
			CALL BootBP(tKelvin,xFrac,nComps,init,pCalc,itMax,yFrac,ier,1)
		endif
		call ErrCheck(ier,iFlag)
		IF(IFLAG.EQ.1)then
			write(66,*)'Error for this point'
			cycle !loop around to next iPt.
		endif
		write(* ,606)tKelvin,pCalc,(xFrac(iComp),iComp=1,nComps-1)
		write(66,606)tKelvin,pCalc,(xFrac(iComp),iComp=1,nComps-1)
	enddo
606   FORMAT(1X,F8.2,1X,F10.7,8(e11.4,1x))      
      close(51)
	close(66)
	if(LOUD)pause 'These results are tabulated in BpFromFile.txt.'
	RETURN
861	continue
	if(LOUD)write(*,*)' Error reading first line of ',TRIM(fn)
	if(LOUD)pause
	return 
862	continue
	if(LOUD)write(*,'(a,i3,a)')' Error reading pt=',i,TRIM(fn)
	if(LOUD)pause
	return 
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE BPITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE BUBBLE PRESSURE EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION X(NMX),Y(NMX),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	!C      if(LOUD)write(*,*)'ENTER DESIRED MAXIMUM FOR NUMBER OF ITERATIONS'
	!C      READ(*,*)MAXIT
	MAXIT=1111
	WRITE(*,*)'INITIALIZATION FROM RAOULTS LAW? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'INITIALIZATION FROM RAOULTS LAW? Y/N?'
	write(52,'(A1)')ANSWER
	INIT=0
	IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
		WRITE(*,*)'ENTER INITIAL GUESS FOR PRESSURE(MPA)'
		READ(*,*)P
		WRITE(52,*)'ENTER INITIAL GUESS FOR PRESSURE(MPA)'
		write(52,*)P
		INIT=1
	ENDIF
	vof=0
	ANSWER='Y'
	DO WHILE(ANSWER.NE.'N'.and.ANSWER.NE.'n')
		WRITE(6,*)'ENTER TEMPERATURE(KELVIN)  '
		READ(5,*)T
		WRITE(52,*)'ENTER TEMPERATURE(KELVIN)  '
		write(52,*)T
		WRITE(6,*)'ENTER LIQD COMPOSITION '
		READ(5,*)(X(J),J=1,NC)
		WRITE(52,*)'ENTER LIQD COMPOSITION '
		write(52,*)(X(J),J=1,NC)
		!renormalize in case mole#'s instead of mole fractions
		sumx=0
		do j=1,nc
			sumx=sumx+x(j)
		enddo
		do j=1,nc
			x(j)=x(j)/sumx
		enddo

		!C      do 6 i=1,1000
		!C      init=0
		IFLAG=0
		ITMAX=MAXIT
		CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)
		!if(ier(1).ne.0)then
			!CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,1)
		!endif

		call ErrCheck(ier,iFlag)
		IF(IFLAG.EQ.1)GOTO 86
		call FlashOut(NC,X,Y,T,P,itmax,vof)
		!C  NOTE:  FOR REPEAT BP CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN
		!C         Routine BEFORE PROCEEDING.
		INIT=1
		WRITE(6,*)'REPEAT BP CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
		WRITE(52,*)'REPEAT BP CALCULATIONS? Y/N?'
		write(52,'(A1)')ANSWER
	ENDDO !WHILE(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GOTO 1000
86	continue
	pause 'Results are in TPXY.txt'
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE BTITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE BUBT EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION X(NMX),Y(NMX),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV

	!C      WRITE(*,*)'ENTER DESIRED MAXIMUM FOR NUMBER OF ITERATIONS'
	!C      READ(*,*)MAXIT
	MAXIT=1111
	WRITE(*,*)'INITIALIZE TEMP FROM VAPOR PRESSURE EQUATION? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'INITIALIZE TEMP FROM VAPOR PRESSURE EQUATION? Y/N?'
	write(52,'(A1)')ANSWER
	INIT=0
	IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
		WRITE(*,*)'ENTER INITIAL GUESS FOR TEMPERATURE(KELVIN) '
		READ(*,*)T
		WRITE(52,*)'ENTER INITIAL GUESS FOR TEMPERATURE(KELVIN) '
		write(52,*)T
		INIT=1
	END IF
	!open(61,file='TPXY.txt')
	!write(61,*)' ' !clear contents for FlashOut to write.
	!close(61)
	vof=0
	ANSWER='Y'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
		WRITE(6,*)'ENTER PRESSURE(MPa)  '
		READ(5,*)P
		WRITE(52,*)'ENTER PRESSURE(MPa)  '
		write(52,*)P

		WRITE(6,*)'ENTER LIQD COMPOSITION '
		READ(5,*)(X(J),J=1,NC)
		WRITE(52,*)'ENTER LIQD COMPOSITION '
		write(52,*)(X(J),J=1,NC)

		IFLAG=0
		ITMAX=MAXIT
		CALL BUBTL(T,X,NC,INIT,P,ITMAX,Y,ier)
		if(ier(1).ne.0)then
		CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
		endif
		call ErrCheck(ier,iFlag)
		call FlashOut(NC,X,Y,T,P,itmax,vof)

		!C  NOTE:  FOR REPEAT BT CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!C         PROGRAM BEFORE PROCEEDING.

		INIT=2

		WRITE(6,*)'REPEAT BT CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
		WRITE(52,*)'REPEAT BT CALCULATIONS? Y/N?'
		write(52,'(A1)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GO TO 1000
	if(LOUD)pause 'Results are in TPXY.txt'
86	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DTITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE DEWT EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION X(NMX),Y(NMX),ier(20)
	!COMMON/ppDataPlus/ZC(NMX),rMw(nmx)
	COMMON/eta/etaL,etaV,ZL,ZV

	!WRITE(*,*)'ENTER DESIRED MAXIMUM FOR NUMBER OF ITERATIONS'
	!READ(*,*)MAXIT
	MAXIT=1111
	WRITE(*,*)'INITIALIZE TEMP FROM VAPOR PRESSURE EQUATION? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'INITIALIZE TEMP FROM VAPOR PRESSURE EQUATION? Y/N?'
	write(52,'(A1)')ANSWER
	INIT=0
	IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
		WRITE(*,*)'ENTER INITIAL GUESS FOR TEMPERATURE(KELVIN) '
		READ(*,*)T
		WRITE(52,*)'ENTER INITIAL GUESS FOR TEMPERATURE(KELVIN) '
		write(52,*)T
		INIT=1
	END IF
	!open(61,file='TPXY.txt')
	!write(61,*)' ' !clear contents for FlashOut to write.
	!close(61)
	vof=1
	ANSWER='Y'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
		WRITE(6,*)'ENTER PRESSURE(MPa)  '
		READ(5,*)P
		WRITE(52,*)'ENTER PRESSURE(MPa)  '
		write(52,*)P

		WRITE(6,*)'ENTER VAPOR COMPOSITION '
		READ(5,*)(Y(J),J=1,NC)
		WRITE(52,*)'ENTER VAPOR COMPOSITION '
		write(52,*)(Y(J),J=1,NC)

		IFLAG=0
		ITMAX=MAXIT
		CALL DEWTV(P,Y,NC,INIT,T,ITMAX,X,ier)
		if(ier(1).ne.0)then
			CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
		endif
		call ErrCheck(ier,iFlag)
		call FlashOut(NC,X,Y,T,P,itmax,vof)

		!C  NOTE:  FOR REPEAT DT CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!C         PROGRAM BEFORE PROCEEDING.
		INIT=1
		T=1.1*T
		WRITE(6,*)'REPEAT DT CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
		WRITE(52,*)'REPEAT DT CALCULATIONS? Y/N?'
		write(52,'(A1)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GO TO 1000
	if(LOUD)pause 'Results are in TPXY.txt'
86	RETURN
	END

!***********************************************************************
	SUBROUTINE FITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE FUGACITY EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	Implicit DoublePrecision(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	character*133 outFile
	LOGICAL LOUDER
	DoublePrecision fugcV(NMX),fugcL(NMX)
	DoublePrecision X(NMX),Y(NMX),VLK(NMX),wFrac(NMX)
	Integer ier(20)
	COMMON/eta/etaL,etaV,zLiqDum,zVapDum
	common/eta2/eta
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data initCall/1/
	LOUDER=LOUD
	LOUDER=.TRUE.
	outFile=TRIM(masterDir)//'\output\PropTable.txt'
	OPEN(6123,file=outFile)
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND PRESSURE(MPa)   '
	READ(5,*)T,P
	WRITE(52,*)'ENTER TEMPERATURE(KELVIN) AND PRESSURE(MPa)   '
	write(52,*)T,P

	IFLAG=0
	write(*,*)' Enter phase type for lower phase (1 for LIQ, 0 or VAP)'
	read(*,*)LIQ
	WRITE(6,*)'ENTER Lower Phase COMPOSITIONS '
	READ(5,*)(X(J),J=1,NC)
	WRITE(52,*)'LIQD COMPOSITIONS '
	write(52,*)(X(J),J=1,NC)
	rMwAvg=0
	do i=1,nc
		rMwAvg=rMwAvg+x(i)*rMw(i)
	enddo
	do i=1,nc
		wFrac(i)=x(i)*rMw(i)/rMwAvg
	enddo
	write(*,*)'Wt Fracs'
	write(*,'(11f8.5)')(wFrac(i),i=1,nc)
	if(initCall)print*,'Fiter: calling Fugi for liq. T(K)=',T
	CALL Fugi(T,P,X,NC,LIQ,fugcL,ZL,ier)	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(initCall)write(*,'(a,f8.2,E12.4)')' Fiter: called Fugi for liq. T(K),ZL=',T,ZL
	etaLo=eta
	rhoLo=P*rMwAvg/(zL*rGas*T)
	iFLAG=0
	DO I=1,6
		IF (ier(I) > 9)IFLAG=1
	enddo
	IF (IFLAG==1 )then
		WRITE(*,'(a,22I4)')'ier',(ier(I),I=1,12)
		WRITE(*,*)'error in lower phase calculation'
		goto 81 
	endif
	DHL=hRes_RT  !these are passed through GlobConst
	DSL=sRes_R
	CvL=CvRes_R
	CpL=CpRes_R
	cmL=cmprsblty

	write(6,*)'Upper phase? Enter 0 for vapor, 1 for liquid.'
	read(*,*)iPhaseUp
	write(52,*)'Upper phase? Enter 0 for vapor, 1 for liquid.'
	write(52,*)iPhaseUp

	WRITE(6,*)'ENTER Upper PHASE COMPOSITIONS    '
	READ(5,*)(Y(I),I=1,NC)
	WRITE(52,*)'Upper PHASE COMPOSITIONS    '
	write(52,*)(Y(I),I=1,NC)
	CALL Fugi(T,P,Y,NC,iPhaseUp,fugcV,ZV,ier)
	rMwAvg=0
	do i=1,nc
		rMwAvg=rMwAvg+y(i)*rMw(i)
	enddo
	etaUp=eta
	rhoUp=P*rMwAvg/(zV*rGas*T)
	if(iPhaseUp.eq.1)etaUp=etaL
	IFLAG=0
	DO I=1,6
		IF (ier(I) > 9)IFLAG=1
	enddo
	IF (IFLAG==1 )then
		WRITE(*,'(a,22I4)')'ier',(ier(I),I=1,12)
		WRITE(*,*)'error in upper phase calculation'
		goto 81 
	endif
	!for Fugi:
	!C  ier = 1 - AT LEAST ONE ERROR
	!C        2 - NOT USED
	!C        3 - NOT USED
	!C        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!C        5 - rho IS -VE
	!C        6 - TOO MANY Z ITERATIONS
	!C        7 - eta > 0.53
	!C		11 - goldenZ instead of real Z.

	IF(ier(2)==1 .and. LOUDER)write(*,*)'ERROR ALL UPPER PHASE'
	IF(ier(3)==1 .and. LOUDER)write(*,*)'ERROR ALL LOWER PHASE'
	!C      PAUSE
	rLogInfinity=707				!larger #'s in xls give #VALUE
	expLogInf=EXP(rLogInfinity)		!useful for trapping overflows in log's and exp's
	DO I=1,NC
		write(*,602)NAME(I),ID(I)
		write(6123,602)NAME(I),ID(I)
		activity=fugcL(I)-fugcV(I)
		if(activity.gt.rLogInfinity)then			!polymers can cause this error.
			pause 'Fiter warning: exp overflow.  Setting kRatio to rLogInfinity.'
			activity=rLogInfinity
		endif
		VLK(I)=EXP(activity)
		write(*,603)
		write(6123,603)
		if(x(i).lt.1/expLogInf)x(i)=1/expLogInf  !fix if x(i)=0 for one component.
		activity=LOG(x(i)*p)+fugcL(i)	!log of the product is the sum of the logs.
		if(activity.gt.rLogInfinity)then			!polymers can cause this error.
			pause 'Fiter warning: exp overflow.  Setting activity to rLogInfinity.'
			activity=rLogInfinity
		endif
		fugLi=EXP(activity)
		if(y(i).lt.1/expLogInf)y(i)=1/expLogInf  !fix if x(i)=0 for one component.
		activity=LOG(y(i)*p)+fugcV(i)	!log of the product is the sum of the logs
		if(activity.gt.rLogInfinity)then
			pause 'Fiter warning: exp overflow.  Setting activity to rLogInfinity.'
			activity=rLogInfinity
		endif
		fugVi=EXP(activity)
		write(*,604)X(I),Y(I),fugLi,fugVi,fugcL(I),fugcV(I)
		write(6123,604)X(I),Y(I),fugLi,fugVi,fugcL(I),fugcV(I)
		write(*,*)'  '
		write(*,605)
		write(*,606)VLK(I),T,P,rhoLo,rhoUp,ZL,ZV
		write(6123,605)
		write(6123,606)VLK(I),T,P,rhoLo,rhoUp,ZL,ZV
		write(*,*)' '
		!if(LOUDER)PAUSE
	enddo
	write(*,*)'DEPARTURE FUNCTIONS:'    
	write(*,*)' (HLo-Hig)/RT SresLo/R CvResLo/R  CpResLo/R  dP/dRho/RT    '
	write(*,610)DHL,DSL,CvL,CpL,cmL
	write(*,*)' (HUp-Hig)/RT SresUp/R CvResUp/R CpResUp/R dP/dRho/RT  '
	write(*,610)hRes_RT,sRes_R,CvRes_R,CpRes_R,cmprsblty  !vapor is latest call, so just print from GlobConst
	!write(6123,*)'DEPARTURE FUNCTIONS:'    
	write(6123,*)'    (H-Hig)/RT   Sres/R   CvRes/R CpRes/Rb dP/dRho/RT    '
	write(6123,'(a,4x,10(F9.4,1X))')' L: ',DHL,DSL,CvL,CpL,cmL
	!write(6123,*)' (HUp-Hig)/RT SresUp/R CvResUp/R  CpResUp/R  cmpUp/R  '
	write(6123,'(a,4x,10(F9.4,1X))')' V: ',hRes_RT,sRes_R,CvRes_R,CpRes_R,cmprsblty  !vapor is latest call, so just print from GlobConst
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(2X,'LIQUID X',2X,'Upper Y',1X,'fugacityLo(MPa)',1X,'fugacityUp(MPa)    lnPhiLo',3X,'lnPhiUp')
605 FORMAT(3X,'K-RATIO',7X,'T(K)',3X,'P(MPa)',2X,'rhoLo',4X,'rhoUp',6X,'ZLo  ',6X,'ZUp')
604	FORMAT(2(3x,F7.4),2(3x,e11.4),2(3x,f12.7))
606	FORMAT(e12.5,3X,F7.2,1X,F7.4,1X,2(F6.4,3X),2(F11.8,1X))
610	FORMAT(4x,10(F9.4,1X))
	initCall=0      
81 WRITE(6,*)'REPEAT FUGACITY CALCULATIONS? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'REPEAT FUGACITY CALCULATIONS? Y/N?'
	write(52,'(A1)')ANSWER
	IF(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GOTO 1000
	close(6123)
	write(*,*)'outFile=',TRIM(outFile)
	pause 'Success! Check results in outFile.'
	RETURN
	END	 !FITER

!***********************************************************************
	SUBROUTINE FVITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE FUGACITY EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2022
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	Implicit DoublePrecision(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	character*133 outFile
	LOGICAL LOUDER
	DoublePrecision fugcL(NMX) !fugcV(NMX),
	DoublePrecision X(NMX),VLK(NMX),wFrac(NMX) !,Y(NMX)
	!Integer ier(20)
	COMMON/eta/etaL,etaV,zLiqDum,zVapDum
	common/eta2/eta
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data initCall/1/
	LOUDER=LOUD
	LOUDER=.TRUE.
	outFile=TRIM(masterDir)//'\output\PropTable.txt'
	OPEN(6123,file=outFile)
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND packing fraction   '
	READ(5,*)T,eta
	WRITE(52,*)'ENTER TEMPERATURE(KELVIN) AND packing fraction   '
	write(52,*)T,eta

	IFLAG=0
	WRITE(6,*)'ENTER COMPOSITIONS '
	READ(5,*)(X(J),J=1,NC)
	WRITE(52,*)' COMPOSITIONS '
	write(52,*)(X(J),J=1,NC)
	rMwAvg=0
	do i=1,nc
		rMwAvg=rMwAvg+x(i)*rMw(i)
	enddo
	do i=1,nc
		wFrac(i)=x(i)*rMw(i)/rMwAvg
	enddo
	write(*,*)'Wt Fracs'
	write(*,'(11f8.5)')(wFrac(i),i=1,nc)
	if(initCall)print*,'FViter: calling Fugi for liq. T(K)=',T
	bMix=SUM( x(1:Nc)*bVolCc_mol(1:Nc) )
	vTotCc=bMix/eta
	isZiter=0
	Call FuVtot(isZiter,T,vTotCc,x,NC,FUGCL,ZL,iErr)
	if(initCall)write(*,'(a,f8.2,E12.4)')' FViter: called Fugi for liq. T(K),ZL=',T,ZL
	etaLo=eta
	rhoLo=1/vTotCc
	P = rhoLo*ZL*Rgas*T
	IF (iErr > 0 )then
		WRITE(*,'(a,22I4)')'iErr=',iErr
		WRITE(*,*)'error in calculation'
		goto 81 
	endif
	DHL=hRes_RT  !these are passed through GlobConst
	DSL=sRes_R
	CvL=CvRes_R
	CpL=CpRes_R
	cmL=cmprsblty

	rLogInfinity=707				!larger #'s in xls give #VALUE
	expLogInf=EXP(rLogInfinity)		!useful for trapping overflows in log's and exp's
	DO I=1,NC
		write(*,602)NAME(I),ID(I)
		write(6123,602)NAME(I),ID(I)
		activity=fugcL(I)-1
		if(activity.gt.rLogInfinity)then			!polymers can cause this error.
			pause 'Fiter warning: exp overflow.  Setting kRatio to rLogInfinity.'
			activity=rLogInfinity
		endif
		VLK(I)=EXP(activity)
		write(*,603)
		write(6123,603)
		if(x(i).lt.1/expLogInf)x(i)=1/expLogInf  !fix if x(i)=0 for one component.
		activity=LOG(x(i)*p)+fugcL(i)	!log of the product is the sum of the logs.
		if(activity.gt.rLogInfinity)then			!polymers can cause this error.
			pause 'Fiter warning: exp overflow.  Setting activity to rLogInfinity.'
			activity=rLogInfinity
		endif
		fugLi=EXP(activity)
		write(*,'(f7.4,E12.4,f12.4,f12.7,E12.5)')X(I),fugLi,fugcL(I),VLK(i)
		write(6123,'(f7.4,E12.4,f12.4,f12.7,E12.5)')X(I),fugLi,fugcL(I),VLK(i)
		!write(6123,604)X(I),Y(I),fugLi,fugVi,fugcL(I),fugcV(I)
		!if(LOUDER)PAUSE
	enddo
	!write(*,*)'DEPARTURE FUNCTIONS:'    
	write(*,*)' (HLo-Hig)/RT SresLo/R CvResLo/R  CpResLo/R  dP/dRho/RT    Z       rho(mol/cc)'
	write(*,610)DHL,DSL,CvL,CpL,cmL,ZL,1/vTotCc
	!write(*,*)' (HUp-Hig)/RT SresUp/R CvResUp/R CpResUp/R dP/dRho/RT  '
	!write(*,610)hRes_RT,sRes_R,CvRes_R,CpRes_R,cmprsblty  !vapor is latest call, so just print from GlobConst
	!write(6123,*)'DEPARTURE FUNCTIONS:'    
	write(6123,*)'    (H-Hig)/RT   Sres/R   CvRes/R CpRes/R dP/dRho/RT    Z       rho(mol/cc)'
	write(6123,'(a,4x,10(F9.4,1X))')' ',DHL,DSL,CvL,CpL,cmL,ZL,1/vTotCc
	!write(6123,*)' (HUp-Hig)/RT SresUp/R CvResUp/R  CpResUp/R  cmpUp/R  '
	!write(6123,'(a,4x,10(F9.4,1X))')' V: ',hRes_RT,sRes_R,CvRes_R,CpRes_R,cmprsblty  !vapor is latest call, so just print from GlobConst
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(2X,'LIQUID X',1X,'fugacityLo(MPa)',1X,'   lnPhiLo',3X,' ')
605 FORMAT(3X,'K-RATIO',7X,'T(K)',3X,'P(MPa)',2X,'rhoLo',4X,'rhoUp',6X,'ZLo  ',6X,'ZUp')
604	FORMAT(2(3x,F7.4),2(3x,e11.4),2(3x,f12.7))
606	FORMAT(e12.5,3X,F7.2,1X,F7.4,1X,2(F6.4,3X),2(F11.8,1X))
610	FORMAT(4x,10(F9.4,1X))
	initCall=0      
81 WRITE(6,*)'REPEAT FUGACITY CALCULATIONS? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'REPEAT FUGACITY CALCULATIONS? Y/N?'
	write(52,'(A1)')ANSWER
	IF(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GOTO 1000
	close(6123)
	write(*,*)'outFile=',TRIM(outFile)
	pause 'Success! Check results in outFile.'
	RETURN
	END	 !FVITER
!***********************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																			   C
!C  PURPOSE:  CRITICAL POINT CALCULATION FOR MIXTURES 						   C
!C  PROGRAMMED BY:  AFG 2011												   C
!C																			   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE CRITER(NC)		   
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER outfile*100
	DIMENSION X(NMX),gmol(NMX)
	COMMON/ETA2/ETA
	COMMON/eta/etaL,etaV,zLiqDum,zVapDum
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	
	WRITE(*,*)'ENTER MOLE FRACTIONS '
	READ(*,*)(x(J),J=1,NC)
	do i=1,nc
		gmol(i)=x(i)
		totMoles=totMoles+gmol(i)  
	enddo
	isZiter=0
	toll=1.d-10
	CALL CP_Mixture(isZiter,toll,gmol,NC,TC_mix,VC_mix,PC_mix,ZC_mix,iErrCode)	  
	if (iErrCode.eq.0) then
	    outFile=TRIM(masterDir)//'\output\CriticalPoint.txt'
	    open(61,file=outFile)
		if (iEosOpt.EQ.5.or.iEosOpt.eq.9) then
			write(61,600)
 			if(LOUD)write(*,*)'Tc=',TC_mix,' etaC=',eta
			if(LOUD)write(*,*)'Pc=',PC_mix,' Zc=',ZC_mix
			write(61,602)TC_mix,PC_mix,ZC_mix,eta
		else
 			write(61,601)
 			if(LOUD)write(*,*)'Tc=',TC_mix,' Pc=',PC_mix
			if(LOUD)write(*,*)'Vc=',VC_mix,' Zc=',ZC_mix
			write(61,602)TC_mix,PC_mix,ZC_mix,VC_mix
		endif
		close(61)
	 	if(LOUD)pause 'Results are saved as output/CriticalPoint.txt'
	endif

600	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, '    EtaC')
601	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, 'VC (cc/mol)')
602 format(2x,f8.3,5x,f6.4,5x,f5.4,5x,f12.5)

	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																			   C
!C  PURPOSE:  CRITICAL POINT CALCULATION FOR PURE COMPOUNDS					   C
!C  PROGRAMMED BY:  AFG 2010												   C
!C																			   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE CPITER(NC)		!AUG 10
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER outfile*100
	COMMON/ETA2/ETA
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT

	isZiter=0
	toll=1.d-4
	!CALL CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,Acen_pure,iErrCode)
	!SUBROUTINE CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,iErrCode)
	CALL CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,Acen_pure,iErrCode)	   
	if (iErrCode.eq.0) then
	    outFile=TRIM(masterDir)//'\output\CriticalPoint.txt'
	    open(61,file=outFile)
		if (iEosOpt.EQ.5.or.iEosOpt.eq.9) then
			write(61,600)
            eta=bVolCC_mol(1)/VC_Pure
 			if(LOUD)write(*,*)'Tc=',TC_Pure,' etaC=',eta
			if(LOUD)write(*,*)'Pc=',PC_pure,' Zc=',ZC_Pure
			rhoc=PC_pure*rMw(1)/(TC_Pure*8.31434*ZC_Pure)
			if(LOUD)write(*,*)'rhoc(g/cc)=',rhoc
			write(61,602)TC_Pure,PC_PURE,ZC_PURE,eta,rhoc
		else
 			write(61,601)
 			if(LOUD)write(*,*)'Tc=',TC_Pure,' Pc=',PC_Pure
			if(LOUD)write(*,*)'Vc=',VC_Pure,' Zc=',ZC_Pure
			write(61,602)TC_Pure,PC_PURE,ZC_PURE,VC_Pure
		endif
		close(61)
	 	if(LOUD)pause 'Results are saved as output/CriticalPoint.txt'
		print*,'Replace Tc,Pc,Zc with CritPure values(0=No)? Improves vp calcs near critical.'
		read*,iAnswer
		if(iAnswer)then
			TC(1)=TC_PURE
			PC(1)=PC_PURE
			ZC(1)=ZC_PURE
		endif
	endif
		 	
600	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, '    EtaC')
601	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, 'VC (cc/mol)')
602 format(2x,f8.3,5x,f6.4,5x,f5.4,5x,f12.5,5x,f12.5)

    RETURN
	END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine FlashOut(nComps,xFrac,yFrac,tKelvin,pMpa,itmax,vof)
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	implicit doublePrecision(a-h,o-z)
	character outFile*251
	LOGICAL LOUDER
	dimension xFrac(nComps),yFrac(nComps)
	COMMON/eta/etaL,etaV,ZL,ZV
	LOUDER=LOUD
	LOUDER=.TRUE.
	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	IF(DEBUG)outFile='C:\SPEAD\CALCEOS\output\TPXY.txt'
	open(61,file=outFile,ACCESS='APPEND')
	if(LOUDER)write(* ,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX
	if(LOUDER)write(* ,'(4(5x,a))')' zLo ',' zUp ','etaLo','etaUp'
	if(LOUDER)write(* ,'(4(f9.6,1x))')zL,zV,etaL,etaV
	if(LOUDER)write(* ,'(a,f8.2,a,f9.4)')' T(K)= ',tKelvin,' P(MPa)=',pMpa
	if(pMpa.lt.1E-3 .or. pMpa > 10E4  .and. LOUDER)write(* ,*)'PMPa=',pMpa
	WRITE(61,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX
	write(61,'(4(5x,a))')' zLo ',' zUp ','etaLo','etaUp'
	write(61,'(4(f9.6,1x))')zL,zV,etaL,etaV
	write(61,'(a,f8.2,a,f9.4)')' T(K)= ',tKelvin,' P(MPa)=',pMpa
	if(pMpa.lt.1E-3 .or. pMpa.gt.10E4)write(61,*)'PMPa=',pMpa
	avgMwLiq=0
	avgMwVap=0
	do i=1,nComps
		avgMwLiq=avgMwLiq+xFrac(i)*rMw(i)
		avgMwVap=avgMwVap+yFrac(i)*rMw(i)
	enddo

	DO I=1,nComps
		if(LOUDER)write(* ,602)NAME(I),ID(I)
		WRITE(61,602)NAME(I),ID(I)
		VLKI=yFrac(I)/(xFrac(I)+1.D-11)
		wFrVap=yFrac(i)*rMw(i)/avgMwVap
		wFrLiq=xFrac(i)*rMw(i)/avgMwLiq
		if(LOUDER)write(* ,603)
		if(xFrac(i).gt.1E-7)then
			if(LOUDER)write(* ,604)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
			WRITE(61,604)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
		else
			if(LOUDER)write(* ,605)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
			WRITE(61,605)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
		endif
		if(LOUDER)write(* ,*)' '					   
		WRITE(61,603)
		WRITE(61,*)' '					   
	enddo
	if(LOUDER)write(* ,182)VOF
	WRITE(61,182)VOF
	close(61)
182	FORMAT(2X,'UpperMoles/Feed = ',F11.5)
602	FORMAT('  COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(3X,'  LIQUID X  ',3X,'Upper Y',3X,' Yi/Xi ',5X,' wFrLo ',6X,'wFrUp')
604	FORMAT(1X,(F11.8,1X),2(E11.4,1X),F11.8,1X,E11.4)
605	FORMAT(1X,(E11.4,1X),2(E11.4,1X),E11.4,1X,E11.4)
	return
	end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE FREEZINGCURVECALC(NC)
	!
	!  PURPOSE:  iterative solubility calculation using activity coefficients for crystalline solids

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*64 Input_File_Name
	character outFile*251
	DIMENSION fugcL(NMX),X(NMX),ier(12)
	DIMENSION XPure(NC)
	!Dimension fugcV(NMX),Y(NMX),VLK(NMX)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!WRITE(*,*)'ENTER delHfusion(J/mol) for solubility equation.'
	!read(*,*)delHfus
	!WRITE(*,*)'ENTER melting Temperature(K) for SOLID'
	!read(*,*)Tmelting
	y1CalcOld=0
	!IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')
	!		call KIJOPTIMUM(NC)
	!CONTINUE

	if(LOUD)write(*,*)'What is the experimental data file name (used for initial guess)>:(?'
	if(LOUD)write(*,*)'The format of the file is:'
	if(LOUD)write(*,*)'First Line:	# of data points, delta heat of fusion(J/mol),Tmelting of the solid(K) '
	if(LOUD)write(*,*)'Following Lines are data points in the format as follows:'
	if(LOUD)write(*,*)'Temp(K)	Pres(MPa)	xSolid xFirstSolvent'
	READ(*,*)Input_File_Name
		outFile=TRIM(masterDir)//'\output\Freezing_Results.txt'
		open(70,file=outFile)
    	!open(70, file="Freezing_Results.txt")
		WRITE(70,*)'T	P	Gamma	x(1)'
	open(68, file=Input_File_Name)	
	read(68,*)nDataPoints, delHfus, Tmelting
	DO Iter_File_Data=1,nDataPoints
		ierScf=0
		READ(68,*)T,P,(X(J),J=1,NC)
		XPure(1)=1
		do i=2,NC
			XPure(i)=0
		enddo		
		CALl Fugi(T,P,XPure,NC,1,fugcL,ZL,ier)
		PhiPureSolid=exp(fugcL(1))
		!     WRITE(6,*)'ENTER Guess for SCF COMPOSITIONS (comp1 = solid)'
		!     READ(5,*)(X(J),J=1,NC)	Disabled when I started reading data from the file
		dev=1111
		iter=0
		do while(ierScf.eq.0 .and.abs(dev).gt.1.d-5)
			iter=iter+1
			if(iter.gt.1111)ierScf=2
			IFLAG=0
			CALl Fugi(T,P,X,NC,1,fugcL,ZL,ier)
			DO I=1,12
				IF (ier(I).NE.0)IFLAG=1
			enddo
			IF (IFLAG.EQ.1)then
				if(LOUD)write(*,'(a,22I4)')'ier',(ier(I),I=1,12)
				if(LOUD)write(*,*)'error in lower phase calculation'
				ierScf=1
			endif
			Gamma=exp(fugcL(1))/PhiPureSolid
			x1Calc=exp((-delHfus)/(rGas*Tmelting)*((Tmelting/T)-1))/Gamma
			dev=x1Calc-x(1)
			if(LOUD)write(*,*)'iter      x(1),       x1Calc         '
			if(LOUD)write(*,'(i4,4e12.4)')iter,x(1),x1Calc,Gamma
			!pause		
			if(x1Calc.ge.0 .and. x1Calc .le. 1)then
				x(1)=x1Calc			
			endif
		enddo !while

		if(ierScf==1  .and. LOUD)write(*,*)'Error returned from FUGI.'
		if(ierScf==2  .and. LOUD)write(*,*)'Error: too many iterations.'

		DHL=DHONKT
		DSL=DSONK
		!for Fugi:
		!C  ier = 1 - AT LEAST ONE ERROR
		!C        2 - NOT USED
		!C        3 - NOT USED
		!C        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
		!C        5 - rho IS -VE
		!C        6 - TOO MANY Z ITERATIONS
		!C        7 - eta > 0.53
		!C		11 - goldenZ instead of real Z.
		IF(ier(2).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL UPPER PHASE'
		IF(ier(3).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL LOWER PHASE'
		if(ierScf.eq.0 .and. LOUD)then
			DO I=1,1
				WRITE(*,602)NAME(I),ID(I)
				WRITE(70,700)T,P,Gamma,X(1)
				WRITE(*,603)
				fugLi=X(i)*p*EXP(fugcL(i))
				WRITE(*,604)X(I),fugLi,fugcL(I)
				WRITE(*,*)'  '
				WRITE(*,605)
				WRITE(*,606)T,P,etaL,ZL
				  
				WRITE(*,*)' '
				!PAUSE
			enddo    
			WRITE(*,*)'DEPARTURE FUNCTIONS: (Hscf-Hig)/NkT  (Sscf-Sig)/Nk'
			WRITE(*,610)DHL,DSL
		endif
	END DO
700 FORMAT(f5.1, 1X, f8.4, 1X,f7.4, 1X,f7.5)
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'  yFrac',3X,'  Fugacity(MPa)',5X,'  lnPhi')
604	FORMAT(1(F12.10),2(5x,e11.4),2(1x,E11.4))
605	FORMAT(3X,'T(K)   ',4X,'P(MPa)',3X,'eta',6X,'zFactor')
606	FORMAT(f7.2,3X,F9.3,1X,F8.5,1X,2(F8.5,3X),2(F9.4,1X))
607	FORMAT(f4.3,' ',f7.5,' ',2(f12.9,' '),4(f12.9,' '))
610	FORMAT(21X,4(F10.4,4X))      
501	FORMAT(A1)
	CLOSE(68)
	CLOSE(70)
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE FreezOptKij(NC)
	!!! For kij optimization in freezing curve calculation
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	DOUBLE PRECISION MinLocal, Next_to_Min
	CHARACTER Input_File_Name*64,outFile*234 
	common XCompErrorFO
	COMMON/eta/etaL,etaV,ZL,ZV

	WRITE(*,*)'What is the experimental data file name (used for initial guess)?'
	WRITE(*,*)'The format of the file is:'
	WRITE(*,*)'First Line:	# of data points, delta heat of fusion(J/mol),Tmelting of the solid(K) '
	WRITE(*,*)'Following Lines are data points in the format as follows:'
	WRITE(*,*)'Temp(K)	Pres(MPa)	xSolid xFirstSolvent'
	READ(*,*)Input_File_Name


	WRITE(6,*)'ENTER RANGE AND KijSTEP FOR KIJ TO OPTIMIZE, (KIJLO KIJHI)  '
	READ(5,*)KIJ0,KIJ3, KIJStep
	outFile=TRIM(masterDir)//'\output\KIJout.txt'
	open(71,file=outFile)
	!OPEN(71,file="KIJOut.txt")
	write(71,*)'Kij	Err'

	!SUBROUTINE FREEZCURVEFromFile(NC, CHARACTER Name_of_File, Error)  
	MinError=100000

	iSteps=(KIJ3-KIJ0)/KIJStep	  
	KIJ_on_this_Loop=KIJ0
	DO IKIJLoop=1,iSteps		
		DO J=2,NC
			DO I=1,J-1			
				KIJ(J,I)=KIJ_on_this_Loop
				KIJ(I, J)=KIJ(J, I)
				KTIJ(J,I)=0
				KTIJ(J,I)=KTIJ(I,J)
			enddo
		enddo
		CALL FREEZCURVEFromFile(NC, Input_File_Name, Error)
		if(XCompErrorFO.lt.MinError)then
			MinError=XCompErrorFO
			Next_to_Min=MinLocal
			MinLocal=KIJ_on_this_Loop
		endif

		write(71,*)KIJ_on_this_Loop, Error
		KIJ_on_this_Loop=KIJ_on_this_Loop + KIJStep
		Error=0.0
	ENDDO
	close(71)
  
	RETURN
	END	 ! FreezOptKij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE FreezCurveFromFile(NC, Input_File_Name, Error)
	!  PURPOSE:  iterative solubility calculation usind activity coefficients for crystalline solids
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	DOUBLE PRECISION Local_SUM_Error
	CHARACTER*64 Input_File_Name
	character outFile*251
	DIMENSION X(NMX),fugcL(NMX),ier(12)
	DIMENSION XPure(NC), XExpt(NC)
	!Dimension fugcV(NMX),Y(NMX),VLK(NMX)
	common XCompErrorFO
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!WRITE(*,*)'ENTER delHfusion(J/mol) for solubility equation.'
	!read(*,*)delHfus
	!WRITE(*,*)'ENTER melting Temperature(K) for SOLID'
	!read(*,*)Tmelting
	y1CalcOld=0
	Local_SUM_Error=0.0
	outFile=TRIM(masterDir)//'\output\Freezing_Results.txt'
	open(70,file=outFile)
	!open(70, file="Freezing_Results.txt")
	WRITE(70,*)'T	P	Gamma	x(1)'
	open(68, file=Input_File_Name)
	read(68,*)nDataPoints, delHfus, Tmelting
	DO Iter_File_Data=1,nDataPoints	
		ierScf=0
		READ(68,*)T,P,(XExpt(J),J=1,NC)
		DO J=1,NC
			X(J)=XExpt(J)
		ENDDO
		XPure(1)=1
		do i=2,NC
			XPure(i)=0
		enddo		
		CALl Fugi(T,P,XPure,NC,1,fugcL,ZL,ier)
		PhiPureSolid=exp(fugcL(1))
		dev=1111
		iter=0
		do while(ierScf.eq.0 .and.abs(dev).gt.1.d-5)
			iter=iter+1
			if(iter.gt.1111)ierScf=2
			IFLAG=0
			CALl Fugi(T,P,X,NC,1,fugcL,ZL,ier)
			DO I=1,12
				IF (ier(I).NE.0)IFLAG=1
			enddo
			IF (IFLAG.EQ.1)then
				if(LOUD)write(*,'(a,22I4)')'ier',(ier(I),I=1,12)
				if(LOUD)write(*,*)'error in lower phase calculation'
				ierScf=1
			endif
			Gamma=exp(fugcL(1))/PhiPureSolid
			x1Calc=exp((-delHfus)/(rGas*Tmelting)*((Tmelting/T)-1))/Gamma
			dev=x1Calc-x(1)
			if(LOUD)write(*,*)'iter      x(1),       x1Calc         '
			if(LOUD)write(*,'(i4,4e12.4)')iter,x(1),x1Calc,Gamma
			!pause		
			if(x1Calc.ge.0 .and. x1Calc .le. 1)then
				x(1)=x1Calc			
			endif
		enddo !while

		if(ierScf.eq.1 .and. LOUD)write(*,*)'Error returned from FUGI.'
		if(ierScf.eq.2 .and. LOUD)write(*,*)'Error: too many iterations.'


		Local_SUM_Error=Local_SUM_Error+abs((x(1) - xExpt(1))/xExpt(1)*100)
		AvgDenom=nDataPoints
		Average=Local_SUM_Error/AvgDenom

	ENDDO
	Error=Average
	CLOSE(68)
	CLOSE(70)
	RETURN
	END !FreezCurveFromFile

!***********************************************************************
	SUBROUTINE GibbsXsEta(NC)
	!C
	!C  PURPOSE:  Gibbs Excess vs. composition at const T,P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER answer*1
	character outFile*251
	DIMENSION X(NMX),gmol(NMX),x1(19)	!,ier(20)
	DIMENSION fugcPure(NMX),fugRepPure(nmx),fugAttPure(nmx),fugAssocPure(nmx),hPure(NMX),sPure(NMX),fugc(NMX)
	DoublePrecision rLnGam(NMX),rLnGamRep(NMX),rLnGamAtt(NMX),rLnGamAssoc(NMX)
	LOGICAL LOUDER
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),ralph(nmx),Zrep,Zatt,Zassoc,Fassoc
	data nx,x1/19,0.0001,.001,.01,.02,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,.98,.99,.999,0.9999/

	LOUDER=LOUD
	LOUDER=.TRUE.

	nComps=NC

	if(nc > 3)then 
		if(LOUDER)pause 'error - this option is only for 2 or 3 compos'
		return
	endif
	!c	write(*,*)'enter file NAME for output'
	!c	read(*,'(A1)')outfile
	outFile=TRIM(masterDir)//'\output\xsGibbs.txt'
	open(61,file=outFile)
!	if(iEosOpt.ne.5.or.iEosOpt.ne.8)then
!		pause 'GX error: this option only for spead model at this time.'
!		return
!	endif
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND eta   '
	READ(5,*)tKelvin,eta
	if(tKelvin < 100)then
		if(LOUDER)write(*,*)'warning T<100K. Sure?(y/n)'
		read(*,*)answer
		if(answer=='n'.or.answer=='N')goto 1000
	endif
	zero=0
	nComps=nc
	if(LOUDER)write(* ,'(a,f7.2,a,f7.3)')' T(K)=',tKelvin,'Eta=',eta
	write(61,*)'  T(K)   Eta'
	write(61,'(f8.2,f8.3)')tKelvin,eta
	tInf=3.0e9
	tStore=tKelvin
	xInf=1e-11
	do iPure=1,nc
		sum=0
		do iComp=1,nc
			x(iComp)=xInf
			sum=sum+x(iComp)
			gmol(iComp)=x(iComp)
		enddo
		x(iPure)=1-sum
		bVolMix=sumproduct(nc,x,bVolCC_mol)
		vTotCc=bVolMix/eta
		!call FuTptVtot(zFactor,aDep,uDep,eta,tKelvin,x,bVolMix,nComps,iErr)
		call FuVtot(isZiter,tStore,vTotCc,x,NC,FUGC,zFactor,iErr)
		if(iErr /= 0)then
			print*,' GibbsXsVsEta: Error returned from pure FuVtot. iPure,Z=',iPure,zFactor
			print*,' Try a larger eta?'
			goto 1000
		endif

		fugcPure(iPure)=FUGC(iPure)-Zfactor-LOG(vTotCc)        ! For the total part, we should also subtract lnZ, but that can give ln(-ve). We really want ln(Vi/V), so we initiate with ln(Vi).. 
		fugRepPure(iPure)=FUGRep(iPure)-Zrep       ! Vk/V = Vk/Vk for pure.	   
		fugAttPure(iPure)=FUGAtt(iPure)-Zatt
		fugAssocPure(iPure)=FUGAssoc(iPure)-Zassoc
		hPure(iPure)=uRes_RT
		call FuVtot(isZiter,tInf,vTotCc,x,NC,FUGC,zFactor,iErr)
		sPure(iPure)=aRes_RT-uRes_RT
	enddo
	!if(nc.eq.2 .and. LOUDER)write(6 ,*)'   x1      Z         uXs      -s0Xs        aXs      gIs         V(cc)','     LnGam'
	if(nc.eq.3 .and. LOUDER)write(6 ,*)'   x1      x2      Z         hE      gMix         gE       gIs        V(cc)       LnGam(i)'
	if(nc.eq.2)write(61,*)'   x1      x2      uXs      -s0Xs        aXs      gIs        Z         V(cc)       LnGam(i)'
	if(nc.eq.3)write(61,*)'   x1      x2      x3      hE      gMix         gE       gIs        Z         V(cc)       LnGam(i)'
	write(* ,*)' xi,PMPa,vTotCc,F,lnGam(1,2),lnGamRep(1,2),lnGamAtt(1,2)',',lnGamAssoc(1,2),ralph(1,2)'
	nX3=nX
	if(nc==2)nX3=1
	do iX3=1,nX3
		do iX1=1,nx
			x(3)=x1(iX3)
			if(nc==2)x(3)=0
			x(1)=x1(iX1)*(1-x(3))
			x(2)=1-(x(1)+x(3))
			sumx=0
			do i=1,nc
                if(LOUDER)then
				    if(x(i).le.0.or.x(i).ge.1)pause 'GibbsXs error in xCalc'
                end if
				sumx=sumx+x(i)
            enddo
            if(LOUDER)then
			    if(ABS(sumx-1).gt.1e-6)pause 'GibbsXs error in xSum'
            end if

			gmol(1:nc)=x(1:nc)
			bVolMix=sumproduct(nc,x,bVolCC_mol)
			vTotCc=bVolMix/eta
			!call FuTptVtot(zFactor,aDep,uDep,eta,tKelvin,x,bVolMix,nComps,iErr)
			call FuVtot(isZiter,tKelvin,vTotCc,x,NC,FUGC,zFactor,iErr)
			!call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErr)
			gIs=0
			gE=0
			gMix=0
			hIs=0
			aIs=0
			sIs=0
			do iComp=1,nc
				gIs=gIs+x(iComp)*LOG(x(iComp))
				hIs=hIs+x(iComp)*hPure(iComp)
				sIs=sIs+x(iComp)*sPure(iComp)
				aIs=aIs+x(iComp)*fugcPure(iComp)
				rLnGamRep(iComp)=fugRep(iComp)-Zrep*bVolCc_mol(iComp)/bVolMix-fugRepPure(iComp)
				rLnGamAtt(iComp)=fugAtt(iComp)-Zatt*bVolCc_mol(iComp)/bVolMix-fugAttPure(iComp)
				rLnGamAssoc(iComp)=fugAssoc(iComp)-Zassoc*bVolCc_mol(iComp)/bVolMix-fugAssocPure(iComp)
				rLnGam(iComp)=fugc(iComp)-Zfactor*bVolCc_mol(iComp)/bVolMix - LOG(vTotCc)-fugcPure(iComp)
				!
				! rLnGam(iComp)=rLnGamRep(iComp)+rLnGamAtt(iComp)+rLnGamAssoc(iComp) + (1-bi/b) -ln(V/Vi) ! gives the same result
			enddo
			hE=uRes_RT-hIs
			aE=aRes_RT-aIs
			call FuVtot(isZiter,tKelvin,vTotCc,x,NC,FUGC,zFactor,iErr)
			!call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErr)
			sDep=aRes_RT-uRes_RT
			sE=sDep-sIs !really this is -ve sXs/R
			!if(LOUDER)write(* ,611)(x(i),i=1,nc-1),zFactor,hE,sE,aE,gIs,vTotCc,(rLnGam(i),i=1,nc)
			PMPa=Zfactor*Rgas*tKelvin/vtotcc
			write(* ,612)x(1),PMPa, vTotCc, Fassoc, (rLnGam(i),i=1,nc),(rLnGamRep(i),i=1,nc),(rLnGamAtt(i),i=1,nc),(rLnGamAssoc(i),i=1,nc),(ralph(i),i=1,nc)
			write(61,611)(x(i),i=1,nc),hE,sE,aE,gIs,zFactor,vTotCc,(rLnGam(i),i=1,nc)
		enddo
	enddo

	CLOSE(61)
	write(*,*)'OutFile=',TRIM(outfile)
	pause 'Check results in OutFile.'   
	RETURN
600	format('   x1  ',3x,'rLnGam1',4x,'rLnGam2',4x,'gE/RT',4x,'hE/RT',6x, 'dGIS/RT' ,6x, 'etaL',6x,'gMix/RT')
601	FORMAT(1X,1(F7.4,1X,F10.3,1X),7(F10.4))
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
605	FORMAT(3X,'K-RATIO',7X,'T(K)',4X,'P(MPa)',2X,'etaL',5X,'etaV',5X,'ZL  ',5X,'ZV')
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))  
611	format(3(f7.4,1x),9(f10.3,1x))    
612	format(f7.4,1x,2(f7.2,1x),12(f7.3,1x))    
	END

!***********************************************************************
	SUBROUTINE GibbsXs(NC)
	!C
	!C  PURPOSE:  Gibbs Excess vs. composition at const T,eta.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	CHARACTER answer*1
	character outFile*251
	Integer ier(20)
	DoublePrecision X(NMX),x1(19),fugcL(NMX)
	DoublePrecision hPure(NMX),fugRepPure(NMX),fugAttPure(NMX),fugAssocPure(NMX),fugcPure(NMX)
	DoublePrecision rLnGam(NMX),rLnGamRep(NMX),rLnGamAtt(NMX),rLnGamAssoc(NMX),Vcc_mol(NMX)
	LOGICAL LOUDER
	COMMON/eta/etaL,etaV,ZL,ZV
	common/FugiParts/fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),ralph(nmx),Zrep,Zatt,Zassoc,Fassoc
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data nx,x1/19,0.0001,.001,.01,.02,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,.98,.99,.999,0.9999/
	LOUDER=LOUD
	LOUDER = .TRUE.
    if(LOUDER)then
	    if(nc.gt.3)pause 'error - this option is only for 2 or 3 compos'
    end if
	if(nc.gt.3)return
	!c	write(*,*)'enter file NAME for output'
	!c	read(*,'(A1)')outfile
	outFile=TRIM(masterDir)//'\output\xsGibbs.txt'
	open(61,file=outFile)
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND PRESSURE(MPa)   '
	READ(5,*)T,P
	if(T < 100)then
		write(*,*)'warning T<100K. Sure?(y/n)'
		read(*,*)answer
		if(answer.eq.'n'.or.answer.eq.'N')goto 1000
	endif
	zero=0
	if(LOUDER)write(* ,'(a,f7.2,a,f7.3)')' T(K)=',T,'P(MPa)=',P
	write(61,*)'  T(K)   P(MPa)'
	write(61,'(f8.2,f8.3)')T,P
	xInf=zeroTol
	do iPure=1,nc
		sum=0
		do iComp=1,nc
			x(iComp)=xInf
			sum=sum+x(iComp)
		enddo
		x(iPure)=1-sum
		call Fugi(t,p,x,nc,1,fugcL,zL,ier)
		!bMix=sumproduct(nc,x,bVolCc_mol) ! no need to calculate this. bMix=bVol(iPure)
		Vcc_mol(iPure)=zL*Rgas*T/P
		fugRepPure(iPure)=fugRep(iPure)-zRep
		fugAttPure(iPure)=fugAtt(iPure)-zAtt
		fugAssocPure(iPure)=fugAssoc(iPure)-zAssoc
		fugcPure(iPure)=fugcL(iPure)
		hPure(iPure)=hRes_RT
	enddo
	if(nc==2 .and. LOUDER)write(6 ,*)'   x1      x2       hE      gMix         gE       gIs        LnGam(i)'
	if(nc==3 .and. LOUDER)write(6 ,*)'   x1      x2      x3      hE      gMix         gE       gIs        LnGam(i)'
	if(nc==2)write(61,*)'   x1      x2       hE      gMix         gE       gIs        LnGam(i)'
	if(nc==3)write(61,*)'   x1      x2      x3      hE      gMix         gE       gIs        LnGam(i)'

	write(* ,'(a,a)')'    xi    eta  vTotCc        F      lnGam(1,2)   lnGamRep(1,2)    lnGamAtt(1,2)','  lnGamAssoc(1,2)    ralph(1,2)'


	nX3=nX
	if(nc.eq.2)nX3=1
	do iX3=1,nX3
		do iX1=1,nx
			x(3)=x1(iX3)
			if(nc.eq.2)x(3)=0
			x(1)=x1(iX1)*(1-x(3))
			x(2)=1-(x(1)+x(3))
			sumx=0
			do i=1,nc
                if(LOUDER)then
				    if(x(i).le.0.or.x(i).ge.1)pause 'GibbsXs error in xCalc'
                end if
				sumx=sumx+x(i)
            enddo
            if(LOUDER)then
			    if(ABS(sumx-1).gt.1e-6)pause 'GibbsXs error in xSum'
            end if
			call Fugi(t,p,x,nc,1,fugcL,zL,ier)
			vTotCc=ZL*Rgas*T/P
			gIs=0
			bMix=SumProduct(nc,bVolCC_mol(1:nc),x(1:nc) )
			eta=bMix/vTotCc
			do iComp=1,nc
				gIs=gIs+x(iComp)*LOG(x(iComp))
				rLnGamRep(iComp)=fugRep(iComp)-Zrep*Vcc_mol(iComp)/vTotCc-fugRepPure(iComp)
				rLnGamAtt(iComp)=fugAtt(iComp)-Zatt*Vcc_mol(iComp)/vTotCc-fugAttPure(iComp)
				rLnGamAssoc(iComp)=fugAssoc(iComp)-Zassoc*Vcc_mol(iComp)/vTotCc-fugAssocPure(iComp)
				!rLnGam(iComp)=fugc(iComp)-Zfactor*bVolCc_mol(iComp)/bVolMix - LOG(vTotCc)-fugcPure(iComp)
				rLnGam(iComp)=fugcL(iComp)-fugcPure(iComp) 
			enddo
			hIs =SumProduct(Nc, x(1:Nc),hPure(1:Nc) )
			gE  =SumProduct(Nc, x(1:Nc),rLnGam(1:Nc) )
			gMix=gIs+gE
			hE=hRes_RT-hIs
			!if(LOUDER)write(* ,611)(x(i),i=1,nc),hE,gMix,gE,gIs,(rLnGam(i),i=1,nc)
			write(61,611)(x(i),i=1,nc),hE,gMix,gE,gIs,(rLnGam(i),i=1,nc)
			write(* ,612)x(1),eta, vTotCc, Fassoc, (rLnGam(i),i=1,nc),(rLnGamRep(i),i=1,nc),(rLnGamAtt(i),i=1,nc),(rLnGamAssoc(i),i=1,nc),(ralph(i),i=1,nc)
		enddo  !iX1
	enddo !iX3

	CLOSE(61)
	write(*,*)'OutFile=',TRIM(outfile)
	pause 'Check results in OutFile.'   
	RETURN
600	format('   x1  ',3x,'rLnGam1',4x,'rLnGam2',4x,'gE/RT',4x,'hE/RT',6x, 'dGIS/RT' ,6x, 'etaL',6x,'gMix/RT')
601	FORMAT(1X,1(F7.4,1X,F10.3,1X),7(F10.4))
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
605	FORMAT(3X,'K-RATIO',7X,'T(K)',4X,'P(MPa)',2X,'etaL',5X,'etaV',5X,'ZL  ',5X,'ZV')
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))  
611	format(3(f7.4,1x),9(f10.3,1x))    
612	format(2(f7.4,1x),(f7.2,1x),12(f7.3,1x))    
	END

!***********************************************************************
      
	SUBROUTINE Isochore(NC,iErrCode,errMsgPas)
	!C
	!C  PURPOSE:  COMPUTE phase envelope GIVEN x,y.
	!C  PROGRAMMED BY:  JRE 6/99, Revised: 3/11	to use FuVtot()
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	parameter(nList=15)
	character outFile*251
    !integer ier(12)
	DIMENSION gMol(NMX),fugc(NMX),tList(nList)  !,ierFugi(12)
	character*77 errMsg(0:11),errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data tList/0.01,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.1/
	iErrCode=0
	errMsg(0)='No Problem in Isochore'
	errMsg(1)='Isochore Error: NC.ne.1'
	errMsg(2)='Isochore Error: itMax exceeded.'
	errMsg(3)='Isochore Error: NC.ne.1'
	if(NC.ne.1)then
		iErrCode=1
		if(LOUD)write(*,*)'Isochore Error: NC=1 only. Restart CalcEos with number of components = 1'
		goto 86
	endif
	gmol(1)=1				     !Arbitrarily setting 1 mole for basis.
	outFile=TRIM(masterDir)//'\output\Isochore.txt'
	open(61,file=outFile)
	!open(666,file='debugData.txt')
	!WRITE(6,*)'ENTER packing fraction   '
	!READ(5,*)eta
	!vTotCc = bVolCc_mol(1)/eta	 !Arbitrarily setting 1 mole for basis.
	!write(6 ,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
	!write(61,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
10	WRITE(6,*)'ENTER density (g/cm3)   '
	READ(5,*)rhoLiqG_cc
	vTotCc = rMw(1)/rhoLiqG_cc	 !Arbitrarily setting 1 mole for basis
	eta=bVolCc_mol(1)/vTotCc
	if(eta > etaMax)print*,'VpIter: etaMax < eta=',eta
	if(eta > etaMax)goto 10
	write(* ,*)'rho(g/cc) =',rhoLiqG_cc,' V(cc/mol) = ',vTotCc
	write(61,*)'rho(g/cc) =',rhoLiqG_cc,' V(cc/mol) = ',vTotCc
	if(LOUD)write(* ,*)'  T(K),      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT,    CvRes/R,     (dP/dRho)/RT'
	write(61,*)'  T(K),      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT,    CvRes/R,     (dP/dRho)/RT'
	iTemp=0
	!do while(iErr.eq.0)
	do iTemp=1,nList
		!iTemp=iTemp+1
		if(iTemp.gt.nList)then
			write(*,*)'Enter value of desired reciprocal reduced temperature.(<0 to stop)'
			read(*,*)trInv
		else
			trInv=tList(iTemp)
		endif
		if(trInv.le.0)exit
		tKelvin=TC(1)/trInv	
		isZiter=1  !avoid ln(-ve Z) during first call.
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,iErrFu)
		if(iErrFu)exit
		pNew = zFactor*8.314*tKelvin/vTotCc
		if(LOUD)write(* ,601)tKelvin,pNew,zFactor,uRes_RT*tKelvin/TC(1),aRes_RT
		write(61,601)tKelvin,pNew,zFactor,uRes_RT*tKelvin/TC(1),aRes_RT
	enddo !new temperatures
	if(LOUD)write(*,*)'  T(K)    PMPa      Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R' !   Sdep'    
	write(61,'(a)')'  T(K)    PMPa      Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R    HRes/RT'    
	print*,'Enter Tstart,Tstop,increment'
	read*,tStart,tStop,tInc
	tKelvin=tStart-tInc
	do while(tKelvin <= tStop)
		tKelvin=tKelvin+tInc
		isZiter=1 !=> derivative props moved separate from ln(Z)
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,iErrFu)
		if(iErrFu)exit
		if(isTPT)call NUMDERVS(NC,gMol,tKelvin,1/vTotCc,zFactor,cmprsblty,iErrNum)	!
!		if(isTPT)print*,'Isotherm: CvRes_R,pas=',CvRes_R
		pNew = zFactor*rGas*tKelvin/vTotCc
		if(LOUD)write(* ,602)tKelvin,pNew,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R
		write(61,602)tKelvin,pNew,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R,hRes_RT
	!DoublePrecision uRes_RT, sRes_R, aRes_RT, hRes_RT, CpRes_R, cvRes_R, cmprsblty !cmprsblty=(dP/dRho)T*(1/RT)
	enddo
601	format(f11.4,5(',',f11.4))
602	FORMAT(2F9.3,1x,1(F7.4),8(F9.4,1X))      
86	continue
	if(LOUD)pause 'Your results are stored in Isochore.txt'
	close(61)
	!close(666)
	errMsgPas=errMsg(iErrCode)
	if(iErrCode.ne.0)then
		if(LOUD)write(*,*)errMsgPas
		if(LOUD)pause
	endif
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE Isobar(NC,iErrCode,errMsgPas)
	!C
	!C  PURPOSE:  COMPUTE phase envelope GIVEN x,y.
	!C  PROGRAMMED BY:  JRE 6/99, Revised: 3/11	to use FuVtot()
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	parameter(nList=15)
	character outFile*251
	LOGICAL LOUDER
    !integer ier(12)
	DIMENSION gMol(NMX),fugc(NMX),tList(nList)  ,ierFugi(12)
	character*77 errMsg(0:11),errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data tList/0.01,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.1/
	data initCall/1/
	LOUDER=LOUD
	LOUDER=.TRUE.
	iErrCode=0
	errMsg(0)='No Problem in IsoBar'
	errMsg(1)='isobar Error: NC.ne.1'
	errMsg(2)='isobar Error: itMax exceeded.'
	errMsg(3)='isobar Error: NC.ne.1'
	if(NC.ne.1)then
		iErrCode=1
		if(LOUDER)write(*,*)'isobar Error: NC=1 only. Restart CalcEos with number of components = 1'
		goto 86
	endif
	gmol(1)=1				     !Arbitrarily setting 1 mole for basis.
	outFile=TRIM(masterDir)//'\output\isobar.txt'
	open(61,file=outFile)
	!open(666,file='debugData.txt')
	!WRITE(6,*)'ENTER packing fraction   '
	!READ(5,*)eta
	!vTotCc = bVolCc_mol(1)/eta	 !Arbitrarily setting 1 mole for basis.
	!write(6 ,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
	!write(61,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
10	WRITE(6,*)'ENTER Pressure (MPa)   '
	READ(5,*)pMPa
	if(LOUDER)write(* ,*)'  T(K),      rhog/cc,     zFactor,   (A-Aig)/RT,  (U-Uig)/RTc,    eta' !   CvRes/R,     (dP/dRho)/RT'
	write(61,'(a155)')'  T(K),     rhog/cc,     zFactor,   (A-Aig)/RT,  (U-Uig)/RTc,      eta'	     
	iTemp=0
	!do while(iErr.eq.0)
	do iTemp=1,nList
		!iTemp=iTemp+1
		if(iTemp.gt.nList)then
			write(*,*)'Enter value of desired reciprocal reduced temperature.(<0 to stop)'
			read(*,*)trInv
		else
			trInv=tList(iTemp)
		endif
		if(trInv.le.0)exit
		tKelvin=TC(1)/trInv	
		isZiter=1  !avoid ln(-ve Z) during first call.
		LIQ=0
		if( trInv > 1)then
			call PsatEar(tKelvin,pSatMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
			if(pMPa > pSatMPa)LIQ=1
		endif
		call FUGI(tKelvin,pMPa,gmol,NC,LIQ,FUGC,zFactor,ierFugi)
		if( ierFugi(1) )exit
		rhoNew = pMPa*rMw(1)/(zFactor*rGas*tKelvin)
		if(LOUDER)write(* ,601)tKelvin,rhoNew,zFactor,aRes_RT,uRes_RT*tKelvin/TC(1),etaPass
		write(61,601)tKelvin,rhoNew,zFactor,aRes_RT,etaPass,uRes_RT*tKelvin/TC(1),etaPass
	enddo !new temperatures
	if(LOUDER)write(*,*)'  T(K)    rhog/cc    Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R' !   Sdep'    
	write(61,'(a)')'  T(K)    rhog/cc    Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes/R    HRes/RT     eta'    
	print*,'Enter Tlo,Thi,increment'
	read*,Tlo,Thi,tInc
	tKelvin=Thi+tInc
	do while(tKelvin >= Tlo)
		tKelvin=tKelvin-tInc
		LIQ=0
		if( tKelvin/ Tc(1) < 1)then
			PsatQnd=Pc(1)*10**( 7*(1+acen(1))/3*(1-Tc(1)/tKelvin) )
			call PsatEar(tKelvin,pSatMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
			if(ierCode .ne. 0)then
				if(ierCode .ne. 7)print*,'Isobar: PsatErrCode=',ierCode	! 7 means you are below PsatMin for the EOS.
				pSatMPa=PsatQnd
			endif
			if(pMPa > pSatMPa)LIQ=1
		endif
		call FUGI(tKelvin,pMPa,gmol,NC,LIQ,FUGC,zFactor,ierFugi)
		if(ierFugi(1))then
			if(LOUDER)print*,'Isobar: T,ierFugi=',tKelvin,ierFugi
			exit
		endif
		if( ABS(zFactor) > 1D-11) rhoMol_cc=pMPa/(zFactor*rGas*tKelvin)
		rhoG_cc = rhoMol_cc*rMw(1)
		eta=etaPass
		call NUMDERVS(NC,gMol,tKelvin,rhoMol_cc,zFactor,cmprsblty,iErrNum)	!It's just way easier and more reliable to get derivative properties numerically.
		if(isTPT.and.initCall)print*,'Isotherm: CvRes_R,pas=',CvRes_R
		if(LOUDER)write(* ,602)tKelvin,rhoNew,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R,eta
		write(61,602)tKelvin,rhoG_cc,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R,hRes_RT,etaPass
		initCall=0
	!DoublePrecision uRes_RT, sRes_R, aRes_RT, hRes_RT, CpRes_R, cvRes_R, cmprsblty !cmprsblty=(dP/dRho)T*(1/RT)
	enddo
601	format(f11.4,6(',',f11.4))
602	FORMAT(F9.3,f9.6,1x,1(F7.4),9(F9.4,1X))      
86	continue
	close(61)
	!close(666)
	IF(LOUDER)PRINT*,'Output file:',TRIM(outFile)
	if(LOUDER)pause 'Success! Check output.'
	errMsgPas=errMsg(iErrCode)
	if(iErrCode.ne.0)then
		if(LOUDER)write(*,*)errMsgPas
		if(LOUDER)pause
	endif
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							 C
!C  PURPOSE:  COMPUTE PRESSURE DENSITY ISOTHERM												 C
!C  MODIFIED BY AFG, Dec. 2009 																 C
!C																							 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE Isotherm(NC)
	USE GlobConst
	USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER (nEta=99)   
	CHARACTER*234 outFile*251
	DIMENSION fugc(NMX),etaList(nEta),gmol(NMX)	  
	CHARACTER*77 errMsg(0:11),errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	COMMON/ETA2/ETA
	COMMON/A/a0Mix,a1Mix,a2Mix,aAtt,a_d,zRep
	COMMON/DervMix2/dz0Mix_deta,da0Mix_deta,dz1Mix_deta,da1Mix_deta,dz2Mix_deta,da2Mix_deta,d2z0Mix_deta2,d2z1Mix_deta2,d2z2Mix_deta2 
	COMMON/DervMix3/d2a1Mix_deta2,d3a1Mix_deta3,d2a2Mix_deta2
	common/VCoeff/B2,B3
	common/TptSaft/A0saft,A1saft,A2saft
	data initCall/1/

	ier=0
	errMsg(0)='No Problem in Isotherm'
	errMsg(1)='Isotherm Error: NC.ne.1'
	errMsg(2)='Isotherm Error: IT option just works with EOS # 4, 5, 9, 10'
	if(NC.ne.1)then
		ier=1
		if(LOUD)write(*,*)'Note: quit and start over with number of components = 1'
		goto 86
	endif
	outFile=TRIM(masterDir)//'\output\Isotherm.txt'
	if(LOUD)write(*,*)'Output in:',TRIM(outFile)
	open(61,file=outFile)
	write(6,*)'ENTER TEMPERATURE '
	read(5,*)tKelvin
	tr=tKelvin/TC(1)
	!tKelvin=TC(1)*tr
	write(61,*)'Temperature =',tKelvin,' Tc = ',TC(1)
	etaList(1)=1D-11
	dEta=etaMax/100
	do I=2,nEta	!nEta is a parameter
		etaList(I)=(I-1)*dEta
	enddo
	gmol(1)=1
	isZiter=1
	if(LOUD)write(* ,600)	!!!!!!!!!!!!!!!!!!!!!!!! Write Header !!!!!!!!!!!!!!!!1
	write(61,600)
	if(bVolCC_mol(1) < 1 .and. LOUD)print*,'Isotherm: 0~bVol=',bVolCc_mol(1)
	if(initCall .and. LOUD)print*,'Isotherm: bVol=',bVolCc_mol(1)
	do i=1,nEta   
		if (I==1) write(61,600)
		eta=etaList(i)
		rhoMol_cc=eta/bVolCC_mol(1)
			if(initCall.and.Loud .and. i < 3)print*,'Isotherm: Calling FuVtot. eta,rhoMol_cc=',eta,rhoMol_cc
		call FuVtot(isZiter,tKelvin,1/rhoMol_cc,gmol,NC,FUGC,Z,iErrFu)
		if(isTPT)call NUMDERVS(NC,gMol,tKelvin,rhoMol_cc,Z,cmprsblty,iErrNum)	!
!		if(isTPT)print*,'Isotherm: cmprsblty,pas=',cmprsblty,cmprsbltyPas
!		if(isTPT)cmprsblty=cmprsbltyPas
		pMPa=Z*rGas*tKelvin*rhoMol_cc
		if(iEosOpt==10)DUONKT=uRes_RT		
		if(iEosOpt==10)DAONKT=aRes_RT		
		if(LOUD)write(* ,601) eta,PMpa,Z,DUONKT,DAONKT,rhoMol_cc,cmprsblty !cmprsblty from GlobConst
		write(61,601) eta,PMpa,Z,DUONKT,DAONKT,rhoMol_cc,cmprsblty
	enddo
	close(61)
	if(LOUD)pause 'Your results are stored in Isotherm.txt'
86	continue
	errMsgPas=errMsg(ier)
	if(ier.ne.0)then
		if(LOUD)write(*,*)errMsgPas
		if(LOUD)pause
	endif
600 format(4x,'eta',7x,'PMPa',9x,'Z',9X,'dU_RT',9x,'dA_RT',9x,'rho(mol/cc)',9x,'cmprsblty')			
601 format(1x,f6.4,6(1x,F11.5),f12.4)			
602 format(1x,f6.4,6(2x,E11.4),e12.4)			
	initCall=0
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							 C
!C  PURPOSE:  compute Virial Coeffs for pure components										 C
!C  PROGRAMMED BY:  AFG,  JAN 2011															 C
!C	AT THIS MOMENT JUST FOR EOS OPTIONS 4 AND 5												 C
!C   																						 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE VirialCoeff(NC)
	USE GlobConst
	USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER (nEta=3,nT=100,etaGrid=1.d-6)   
	CHARACTER*234 outFile*251
	DIMENSION fugc(NMX),etaList(nEta),B2grid(nT),B3grid(nT),deltaEta(nEta),gmol(NMX)!,ierFugi(20)	  
	CHARACTER*150 ErrMsg(4)
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/VCoeff/B2,B3

	ErrMsg(1)='VirialCoeff error - NC>1 , NC must be equal to 1'
	ErrMsg(2)='VirialCoeff error - Please use "VC" option with EOS # 4, 5 & 9'
	if(NC.ne.1)then
		if(LOUD)write(*,*)ErrMsg(1)
		return
	endif
	outFile=TRIM(masterDir)//'\output\VCoeff.txt'
	open(61,file=outFile)
	!WRITE(6,*)'ENTER INITIAL REDUCED TEMPERATURE   '
	!READ(5,*)tr
	Tr = 0.4
	betaRed = 1 / Tr
	etaList(1)=1.d-8	   !There is not a significant difference between 10-6 and 10-7 
	DO I=2,nEta	 ! etaGrid = 1D-6, so we go from ~ 0 to 3E-6 in eta.
		etaList(I)=etaList(I-1)+etaGrid	    
	ENDDO
	gmol(1)=1
	isZiter=1
	deltaBetaRed = 0.1D0
	do J=1,nT
		tKelvin=TC(1)/ betaRed
		if( betaRed < 0.2D0 ) deltaBetaRed = 0.01D0
		if( betaRed < 0.02D0 ) deltaBetaRed = 0.001D0
		if( betaRed < 0.002D0 ) deltaBetaRed = 0.0001D0
		if( betaRed < 0.0001D0 ) exit
		betaRed = betaRed - deltaBetaRed 
		do I=1,nEta   
			eta=etaList(I)
			if (iEosOpt.EQ.4) then
				vTotCc=VX(1)/eta
			elseif (iEosOpt.EQ.5.or.iEosOpt.eq.9) then
				vTotCc=bVolCC_mol(1)/eta
			else
				if(LOUD)write(*,*)ErrMsg(2)
				return
			endif
			CALL FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErr)
			B2grid(I)=(Z-1.d0)/eta
		enddo
		if (J.EQ.1) WRITE(6 ,600)		
		if (J.EQ.1) WRITE(61,600)		
		do k=1,(nEta-1)
			deltaEta(k)=etaList(k+1)-etaList(k)
			B3grid(k)=(B2grid(k+1)-B2grid(k))/deltaEta(k)
		enddo
	!Since Virial Coefficient are defined as eta goes to zero:
		B2=B2grid(1) *bVolCC_mol(1)
		B3=B3grid(1) *bVolCC_mol(1)**2
		WRITE(6  ,601)tKelvin,betaRed,B2,B3	
		WRITE(61,601)tKelvin,betaRed,B2,B3	
	enddo
	close(61)
	if(LOUD)pause 'Your results are stored in VCoeff.txt'

600 format(1x,'T(K)',6x,'1/Tr',4x,'B2(cc/mol)',2x,'B3(cc/mol)^2')			 
!600 format(1x,'1/kT',10x,'B2',10x,'B3')			 
601 format(1x,f8.2,2x,f6.4,4x,f8.2,4x,f10.1)

	return
	end
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
	SUBROUTINE Isotherm_old(NC,ier,errMsgPas)
	!C
	!C  PURPOSE:  COMPUTE phase envelope GIVEN x,y.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER outFile*251
	DIMENSION xFrac(NMX),fugc(NMX),ierFugi(20),etaList(8)
	CHARACTER*77 errMsg(0:11),errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data etaList/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5/
	ier=0
	errMsg(0)='No Problem in Isotherm'
	errMsg(1)='Isotherm Error: NC.ne.1'
	errMsg(2)='Isotherm Error: itMax exceeded.'
	errMsg(3)='Isotherm Error: NC.ne.1'
	if(NC.ne.1)then
		ier=1
		if(LOUD)write(*,*)'Note: quit and start over with number of components = 1'
		goto 86
	endif
	outFile=TRIM(masterDir)//'\output\Isotherm.txt'
	open(61,file=outFile)
	WRITE(6,*)'ENTER temperature   '
	READ(5,*)trInv
	tKelvin=TC(1)/trInv
	write(61,*)'Temperature =',tKelvin,' Tc = ',TC(1)
	!Note: This is wacky since it is easy to compute Z,U given rho,
	!but I don't have a standard interface for ZCalc of all eosOpt's.
	!I only have a standard interface for fugi.
	xFrac(1)=1
	itMax=111
	LIQ=1
	tol=0.001
	pOld=255	!initialize here so we bootstrap
	write(61,*)'  1/Tr,      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT'
	ieta=0
	do while(ier.eq.0)
		pOld=pOld/2  !for boot strapping from hi temp to low, move the next guess down.
		ieta=ieta+1
		if(ieta.gt.8)then
			if(LOUD)write(*,*)'Enter value of desired packing fraction.(<0 to stop)'
			read(*,*)eta
		else
			eta=etaList(ieta)
		endif
		if(trInv.le.0)exit
		tKelvin=TC(1)/trInv
		iter=0
		call Fugi(tKelvin,pOld,xFrac,NC,LIQ,FUGC,zFactor,ierFugi)
		etaOld=etaL
		fOld=etaOld-eta
		pNew=pOld/1.1
		absDev=111
		do while(ier.eq.0 .and. absDev.gt.tol)
			iter=iter+1
			if(iter.gt.itMax)then
				ier=2
				cycle
			endif
			call Fugi(tKelvin,pNew,xFrac,NC,LIQ,FUGC,zFactor,ierFugi)
			fNew=etaL-eta
            if(LOUD)then
			    if(ABS(pNew-pOld).lt.1e-11)pause 'IC error: pNew~pOld'
            end if
			change=fNew/(fNew-fOld)*(pNew-pOld)
			pOld=pNew
			fOld=fNew
			pNew=pNew-change
			absDev=ABS(change/pNew)
		enddo !while(absDev.gt.tol)
		if(LOUD)write(* ,*)'  1/Tr,      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT'
		if(LOUD)write(* ,601)trInv,pNew,zFactor,dUoNkT*tKelvin/TC(1),dAoNkT
		write(61,601)trInv,pNew,zFactor,dUoNkT*tKelvin/TC(1),dAoNkT
	enddo !new temperatures
601	format(f7.4,5(',',f11.4))
86	continue
	if(LOUD)pause 'Your results are stored in Isotherm.txt'
	close(61)
	errMsgPas=errMsg(ier)
	if(ier.ne.0)then
		if(LOUD)write(*,*)errMsgPas
		if(LOUD)pause
	endif
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	SUBROUTINE KIJITR(NC)
	!C
	!C  PURPOSE:  ITERATIVE BP CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) AT ONE 
	!C    DATA POINT.
	!C  PROGRAMMED BY:  JRE 2/95
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION X(NMX),Y(NMX),ier(20),VLK(NMX)
	COMMON/eta/etaL,etaV,ZL,ZV
	MAXIT=1111
	INIT=0
	WRITE(*,*)'THIS ROUTINE ASSUMES A BINARY SYSTEM FROM MAIN PROGRAM'
	WRITE(*,*)'ENTER TEMPERATURE(KELVIN),X(1)  '
	READ(5,*)T,X(1)
	write(52,*)'The number of components is:',NC
	WRITE(52,*)'ENTER TEMPERATURE(KELVIN),X(1)  '
	write(52,*)T,X(1)
	X(2)=1-X(1)

	answer='y'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
		WRITE(*,*)'ENTER KIJ '
		READ(*,*)KIJ(1,2)
		WRITE(52,*)'ENTER KIJ '
		write(52,*)KIJ(1,2)
		KIJ(2,1)=KIJ(1,2)
		IFLAG=0
		ITMAX=MAXIT
		CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)
		DO I=1,6
		IF (ier(I).NE.0)IFLAG=1
		enddo
		IF (IFLAG.EQ.1)WRITE(6,*)'ier',(ier(I),I=1,6)
		IF(ier(2).EQ.1 )WRITE(*,*)'ERROR ALL UPPER PHASE'
		IF(ier(3).EQ.1 )WRITE(*,*)'ERROR ALL LOWER PHASE'

		write(*,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX

		DO I=1,NC
			write(*,602)NAME(I),ID(I)
			VLK(I)=Y(I)/(X(I)+1.D-11)
			write(*,603)
			write(*,604)X(I),Y(I),VLK(I),T,P,ZL,ZV
			WRITE(52,604)X(I),Y(I),VLK(I),T,P,ZL,ZV
			write(*,*)' '
		enddo

		!C  NOTE:  FOR REPEAT KIJ CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!C         PROGRAM BEFORE PROCEEDING.

		INIT=1

		WRITE(6,*)'REPEAT KIJ CALCULATIONS? Y/N?'
		READ(*,'(a)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
	RETURN
602	FORMAT('  COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(3X,'LIQUID X',3X,'VAPOR Y',3X,' Yi/Xi ',6X,'T(K)',8X,'P(MPa)',7X,' ZL ',6X,' ZV ')
604	FORMAT(1X,(F7.4,1X),E11.4,2X,E11.4,1X,2(f10.4,1X),2X,2(F6.4,4X))
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJDB(NC)
	!C
	!C  PURPOSE:  BP CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
	!C  PROGRAMMED BY:  JRE 2/96
	!C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
	!C  REFERENCE:NUMERICAL RECIPES
	!C
	USE PortLib !for time()
	USE GlobConst
	USE BIPs  !includes nConverged,maxPts,tDat,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	character*77 errMsgLookup
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5),temParm(5),TR(NMX)
	Integer idStore(NMX),idcc(2),QueryNParMix !,iErrGet,idCas(NMX)
	LOGICAL LOUDER
	EXTERNAL LMDevFcn
	character*77 inFile,outFile,errMsgPas !,dumString ,name1,name2
	character*1 tabChar
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/FloryWert/ vLiq(nmx)
    data initial/1/
	tabChar = char(9)
	LOUDER=LOUD	!This provides local control if desired.
	!LOUDER=.TRUE. !Uncomment this line for screen writes from this routine.
	!store current values of nc,id.
	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	NC=2
	zero=0
	eightySix= -86
    idOpt=1
	write(*,*)'Enter 1 for dippr ids or 2 for ccids'
	write(*,*)'(samplevle.dat uses ccids)'
!	read(*,*)idOpt
	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in VleData subdir)'
	WRITE(*,*)'FIRST LINE =NPTS,id1,id2 2ND LINE = TK, PMPA, X,Y, ...'
	write(*,*)'Enter 1 for Jaubert, 2 for JaubertAq, 3 for DannerGess, 4 Solvation '
	print*,'Add +10 to skip opt'
	read(*,*)iAns
	iOptimize=1
	if(iAns > 10)then
		iOptimize=0
		iAns=iAns - 10
	elseif(iAns < 0)then
		iOptimize= -1	  ! it means just read the database and rewrite in new format.
		iAns=ABS(iAns)
	endif
    FN="Jaubert.txt"
	if(iAns==2)FN="JaubertAq.txt"
	if(iAns==3)FN="DannerGessVle.txt"															
	if(iAns==4)FN="SolvationSystems.txt"
!	READ(*,'(A)')FN
	inFile=TRIM(masterDir)//'\VleData\'//TRIM(FN)
	if(DEBUG)inFile='c:\spead\calceos\VleData\'//TRIM(FN)
	open(51,file=inFile)

	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFile)
	write(661,*)' ' ! purge old contents of KijOut. KijVal uses APPEND
	if(iOptimize .ge. 0)close(661) ! purge old contents of KijOut. KijVal uses APPEND.

	outFile=TRIM(masterDir)//'\output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\output\KijDb.txt'
	if(initial)then
        initial=0
        open(61,file=outFile)
    else
        open(61,file=outFile,ACCESS='APPEND') ! e.g. JaubertAq is appended to Jaubert if during the same run.
    endif
	write(61,*)'    ID1	  ID2  nConv      Kij       stdErr   %PAADP    rmsErrP    dStdErr'

	if(iEosOpt.eq.4)then
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij'
	elseif(iEosOpt.eq.5.or.iEosOpt.eq.8.or.iEosOpt.eq.9)then	   !AFG 2011 , Added EOS Opt. 9
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij, 4-kij,kTij'
	elseif(iEosOpt.eq.3)then
		write(*,*)'Options:  1-kij, 2-tauIJ,tauJI, 3-kij,tauIJ,tauJI, '
		write(*,*)'          4-kij,kTij,tauIJ,tauJI, 5-kij,kTij,tauIJ,tauJI,alphaIJ'
	elseif(iEosOpt.eq.7)then
		write(*,*)'Options:  2-tauIJ,tauJI, 3-alphaIJ,tauIJ,tauJI'
	else
		write(*,*)'Options:  1-kij, 2-kij,kTij '
    endif
    iRegOpt=1
!	READ(5,*)iRegOpt
	ier=0
    iSystem=0
	iTimeStart=time()
	paadTot=0
	nConvergedTot=0
	do while(ier.eq.0) !loop over data sets in db
		!READ(51,'(4i5,a,a,a77)',ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
		!pause 'KijDb: Reading new mix'
		READ(51,*,ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
        if(ioErr /= 0)then
			timeTot=float(time()-iTimeStart)
			print*,'KijDb: elapsed time(sec)=',timeTot
			if(iSystem > 0)paadTot=paadTot/iSystem
			write(*,'(a,i9,2f10.2)')' Overall nConverged,paad=',nConvergedTot,paadTot
			close(61)
			close(661)
			if(ioErr== -1)then
				pause 'KijDb: end of file.'
			else
				print*,'KijDb: read error. system#,ioErr,File=',iSystem,ioErr,TRIM(inFile)
				pause 'KijDb: check file.'
			endif
			exit ! end of file reached
		endif
        iSystem=iSystem+1
		call IdCasLookup(NC,idCas,ierLookup,errMsgLookup)	! idDippr passed by GlobConst. This also sets the class()
		if(iOptimize < 0)then
			!write(661,'(i4,2i5,2i6,a)')nPtsBipDat,id(1),id(2),idTrc(1),idTrc(2),' ! Line0:nPts,idDip1,idDip2,idTrc1,idTrc2; Line1-nPts:T(K),P(MPa),x1,y1' 	
		endif !iOptimize < 0

		if(idOpt.eq.2)then
			do i=1,NC
				idcc(i)=id(i)
			enddo
		endif
		write(*,*)'KijDb: New mix, NPTS,ID()=',nPtsBipDat,id(1),id(2)
        if(LOUDER)pause 'Starting new binary mixture'
		if(nPtsBipDat.le.0)then
			ier=1 !indicates end of database
			cycle !loop to end will terminate
		endif
!		if(idOpt==2)call IdCcLookup(NC,idcc,ierLookup,errMsgLookup)
!		if(ierLookup.ne.0)then
!			write(61,607)idcc(1),idcc(2),zero,eightySix
!			DO I=1,NPTS	!read the points but do nothing so we can cycle to next.
!				READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
!			enddo
!			cycle !loop to next component but leave ier alone so keep looping.
!		endif
		parm=0
		CALL GETCRIT(NC,iErrCrit)
		if(iErrCrit.ne.0 .and. LOUD)then
			write(*,*)'KijDb Error: ErrorCode from GetCrit = ',iErrCrit
			pause
		endif
		iErrGet=SetNewEos(iEosOpt) !wipe out any instance of other eos.
		iErrGet=1 !should be set to zero by Get__(), so something went wrong if iErrGet>0 after calls.
		if(iEosOpt==1)CALL GetPR(NC,iErrGet)
		if(iEosOpt==2)CALL GetEsdCas(NC,idCas,iErrGet)
		if(iEosOpt==3)CALL GetPRWS(NC,iErrGet)
		if(iEosOpt==5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
		if(iEosOpt==6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
		if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
		if(iEosOpt==10)CALL GetPcSaft(NC,idCas,iErrGet)
		if(iEosOpt==11)CALL GetPrtc(NC,iErrGet)
		if(iEosOpt==17)CALL GetPrLorraine(NC,iErrGet)
		if(iEosOpt==18)CALL GetEsd2Cas(NC,idCas,iErrGet)
		if(iErrGet > 0 .or. iErrCrit.ne.0 .or. ierLookup > 0 .or. iOptimize < 0)then
			write(*,*)'KijDB Error: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(61,607)id(1),id(2),0,zero,eightySix
			DO I=1,nPtsBipDat	!read the points but do nothing so we can cycle to next
				!READ(51,'(a77)')dumString
				READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
				if(iOptimize < 0) write(661,'(2(i6,a),f8.2,a,3(f8.5,a),1x,a)')idTrc(1),tabChar,idTrc(2),tabChar,tDat(I),tabChar,PDAT(I),tabChar,XDAT(I),tabChar,YDAT(I),tabChar !,dumString
			enddo
			cycle !loop to next component but leave ier alone so keep looping.
        endif
        !if(LOUDER)write(*,'(a,11F10.2)')' Vx=',(vx(i),i=1,NC)
		pdatMin=1E11
		pdatMax=0
		DO I=1,nPtsBipDat	!we only get to here if there is no error.  otherwise, cycle takes to enddo while()
			READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I) !	data in module BIPs
			if(pdat(i) < pdatMin)pdatMin=pdat(i)
			if(pdat(i) > pdatMax)pdatMax=pdat(i)
			DO iComp=1,NC
				if(Tc(iComp) < zeroTol)pause 'Main: Tc < 0 for some comp???'
				TR(iComp)=TDAT(I)/TC(iComp)
			ENDDO
		enddo
		nParms=QueryNParMix() ! Fortran requires declaration of integer QueryNParMix above.
		parm(2)=0
		temParm=0
		do i=1,nParms
			call QueryParMix(i,temParm(i),iErrParMix)  
			parm(i)=temParm(i) ! this will take initial guess from database
		enddo
		if(iErrParMix)print*,'KijDb: error returning stored mix parms'
		write(*,'(a,3E12.4)')' KijDb: IniParms=',(temParm(i),i=1,nParms)
		if(iOptimize)then
			MAXIT=111
			INIT=1
			LWA=275	 !C	LWA>5*nParms+NPTS ~ 275
			TOL=0.0001	 !c      TOL=SQRT(DPMPAR(1))
			factor=100 !c  factor controls the magnitude of the first step from the initial guess.
			if(LOUDER)write(*,*)'KIJ,paadp,rmserr'
			if(ABS(temParm(1)) < zeroTol)then	! signal that parameters do not exist for this compound
				bestKij= -0.03
				bestAAD=12345
				do iGrid=1,4
					parm(1)= -0.15d0+(iGrid-1)*0.1d0
					call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,0)	!nConverged passed by BIPs module. 0 means enforce penalty for nonconverged points.
					if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
					write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
					if(PAADP < bestAAD)then
						bestAAD=PAADP
						bestKij=parm(1)
					endif
				enddo
				parm(1)=bestKij
			endif !ParmSrch
			CALL LMDifEZ(LMDevFcn,nPtsBipDat,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
			if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
			if(LOUDER)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		endif ! iOptimize
		if(iEosOpt.le.3.or.iEosOpt.eq.7)then
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ','       ',' xsBIP12',' xsBIP21'
		elseif(iEosOpt.le.5.or. iEosOpt.eq.7.or.iEosOpt.eq.9)then	   !AFG 2011 , Added EOS Opt. 9
			if(LOUDER)write(* ,*)'   KIJ','      ',' HIJ','       ','  kTij'
		else
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ' !,'      ',' xsBIP12',' xsBIP21'
		endif
		!write(* ,'(5f9.4)')(parm(i),i=1,nParms)
		if(LOUDER)write(* ,*)'StdErr'
		if(LOUDER)write(* ,'(5f9.4)')(stdErr(i),i=1,nParms)
		if(LOUDER)pause 'KijDB: Converged. Last call to KijVal.'
		LastCall=1	! Do not apply penalty for unconverged points on LastCall. Just get the deviations for converged points.
		call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,LastCall)	!nConverged passed by BIPs module.
		if(nConverged > 0)rmsErr=SQRT( DFLOAT(nPtsBipDat)/nConverged )*RmsJre(nPtsBipDat,deviate) !adjust for nConverged.nPts needed to sum over all elements.
		!write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==1)write(61,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==2)write(61,608)id(1),id(2),nConverged,(parm(i),i=1,2),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==2)write(6 ,608)id(1),id(2),nConverged,(parm(i),i=1,2),PAADP,rmsErr,DstdErr,name(1),name(2)
		write(6 ,'( 3i5,2F10.2,11(1x,E11.4) )')iSystem,id(1),id(2),PAADP,rmsErr,(parm(i),i=1,nParms)
		paadTot=paadTot+paadP
		nConvergedTot=nConvergedTot+nConverged
	enddo !while(ier.eq.0) loop over all data sets in db

	close(51)
	close(61)
	if(LOUDER)pause 'These results are tabulated in KijDb.txt.'
	!restore stored nc,id
	NC=ncStore
	do i=1,NC
		id(i)=idStore(i)
	enddo
	CALL GETCRIT(NC,iErrCrit)
	call IdCasLookup(NC,idCas,ier,errMsgLookup)	! idDippr comes from GlobConst. This also resets the class(). .
	if(iErrCrit.ne.0)then
		if(LOUDER)write(*,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		if(LOUDER)pause
		stop
	endif 
	if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	if(isESD)CALL GetEsdCas(NC,idCas,iErrGet)
	if(iEosOpt.eq.3)CALL GetPRWS(NC,iErrGet)
	if(isTPT)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
	if(iEosOpt.eq.6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt.eq.7)CALL GetNRTL (NC,ID,iErrGet)
	if(iErrGet.gt.0)then
		if(LOUDER)write(*,*)'KijDB Error: failed to restore required pure props.'
		if(LOUDER)write(*,*)'id1,id2',id(1),id(2)
		if(LOUDER)pause
	endif
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
607	FORMAT(1X,3i6,1X,<nParms>E11.4,1X,<nParms>E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
608	FORMAT(1X,3i6,1X,E11.4,1X,E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
      
	RETURN
	END	   ! KijDb

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJDBLLE(NC)
	!C
	!C  PURPOSE:  FUGI CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED USEing BIPs
	!C  PROGRAMMED BY:  JRE 12/21
	!C  METHOD:   Levenberg-Marquardt 
	!C  REFERENCE:NUMERICAL RECIPES, MinPack
	!C
	USE GlobConst
	USE BIPs  !includes nConverged,maxPts,tDat,...
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	!USE MSIMSL !include IMSL for sort routine: DSROWR
	!USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	parameter(nVarsMax=11)
	CHARACTER*44 FN
	character*77 errMsgLookup
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5),temParm(5),var(maxPts,nVarsMax) !,TR(NMX)
	Integer idStore(NMX),QueryNParMix !,iErrGet,idcc(2),idCas(NMX)
	LOGICAL LOUDER
	EXTERNAL LMDevLLE
	character*77 inFile,outFile,outFilePts,errMsgPas,dumString
	character*3 nexString
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/FloryWert/ vLiq(nmx)
	!common/bipDA/iDA,jDA
	common/DevPas/iOptimize
    data initial/1/
	LOUDER=LOUD	!This provides local control if desired.
	!LOUDER=.TRUE. !Uncomment this line for screen writes from this routine.
	!store current values of nc,id.
	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	NC=2
	zero=0
	eightySix= -86
    idOpt=2
	write(*,*)'Enter 1 for dippr ids or 2 for idCas'
	write(*,*)'(LLeDb.txt uses idCas)'
	iOptimize=1
	print*,'Enter 1 to optimize BIPs or 0 to use stored BIPs.'
	read*,iOptimize
!	read(*,*)idOpt
!	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in VleData subdir)'
!	WRITE(*,*)'FIRST LINE =NPTS,id1,id2 2ND LINE = TK, PMPA, X,Y, ...'
!	write(*,*)'Enter 1 for Jaubert, 2 for JaubertAq, 3 for DannerGess, +10 to skip opt'
!	read(*,*)iAns
!	if(iAns > 10)then
!		iOptimize=0
!		iAns=iAns - 10
!	endif
!	if(iAns==2)FN="JaubertAq.txt"
!	if(iAns==4)FN="JaubertAbb.txt"
!	if(iAns==3)FN="DannerGessVle.txt"
!	READ(*,'(A)')FN
    FN="LLeDbPGL6ed96b.txt"
	inFile=TRIM(masterDir)//'\LLe&SLeData\'//TRIM(FN)
	if(DEBUG)inFile='c:\spead\calceos\LLe&SLeData\'//TRIM(FN)
	inUnit=51
	open(inUnit,file=inFile)

	outFilePts=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFilePts='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFilePts)
		if(iOptimize >1)write(661,*)'   T(K)    ln(K1c)  ln(K1e)  ln(K2c)  ln(K2e)   dev1     dev2 ' ! purge old contents of KijOut. KijVal uses APPEND
		if(iOptimize==0)write(661,*)'   T(K)      x(1)     x1DB       y1     y1DB    dev1     dev2 ' ! purge old contents of KijOut. KijVal uses APPEND
	close(661) ! purge old contents of KijOut. KijVal uses APPEND.

	outFile=TRIM(masterDir)//'\output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\output\KijDb.txt'
	if(initial)then
        initial=0
        open(61,file=outFile)
    else
        open(61,file=outFile,ACCESS='APPEND') ! e.g. JaubertAq is appended to Jaubert if during the same run.
    endif
	write(61,*)'   ID1   ID2  nConv      Kij       stdErr   %PAALDK    rmsErrP    dStdErr'

	if(iEosOpt==4)then
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij'
	elseif(iEosOpt==5.or.iEosOpt==8.or.iEosOpt==9)then	   !AFG 2011 , Added EOS Opt. 9
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij, 4-kij,kTij'
	elseif(iEosOpt==3)then
		write(*,*)'Options:  1-kij, 2-tauIJ,tauJI, 3-kij,tauIJ,tauJI, '
		write(*,*)'          4-kij,kTij,tauIJ,tauJI, 5-kij,kTij,tauIJ,tauJI,alphaIJ'
	elseif(iEosOpt.eq.7)then
		write(*,*)'Options:  2-tauIJ,tauJI, 3-alphaIJ,tauIJ,tauJI'
	else
		write(*,*)'Options:  1-kij, 2-kij,kTij '
    endif
    iRegOpt=1
!	READ(5,*)iRegOpt
	ier=0
    iSystem=1
	nSystems=0 ! only count systems that have coex data.
	nLines=0
	nexString="***"
	nVars=5
	read(inUnit,'(a77)',ioStat=ioErr)dumString ! first line is *** Next *** so just skip it.
	do while(ierRead==0) !loop over data sets in db
		call ReadNext(inUnit,nexString,maxPts,nVarsMax,.TRUE.,dumString,nVars,var,nLines,ierRead)
		print*,'Header= ',TRIM(dumString)
		if(ierRead == -2)exit ! end of file	before reading any lines. ierRead= -1 for end of file at end of dataset.
		read(dumString,*,ioStat=ioErr)i99,i13,idCas(1),idCas(2)
		!sort var by increasing temperature
		iColTK=1
		Call SortTable(var,maxPts,nLines,nVars,iColTK,1) ! 1=yes for LoToHi
		if(ioErr)pause 'KijDbLLE: error reading header of dataset?'
		tMin=1234
		tMax=0
		call IdDipprLookup(NC,idCas,ierLookup,errMsgLookup)	! idDippr goes to GlobConst. This also sets the class()
		if(ierLookup > 0)then
			print*,'KijDbLLE: ID fail!',idCas(1),idCas(2),ID(1),ID(2)
			print*,TRIM(errMsgLookup)
			cycle
		endif
		!print*,'ID=',ID(1),ID(2)
		!call IdCasLookup(NC,idCas,ierLookup,errMsgLookup)	
		!if(LOUDER)pause 'Starting new binary mixture'
		CALL GETCRIT(NC,iErrCrit)
		if(iErrCrit.ne.0)then
			print*,'ID()=',ID(1),ID(2)
			write(*,*)'KijDbLLE Error: iErrCit,Tc() = ',iErrCrit,Tc(1),Tc(2)
			!pause
			cycle
		endif
		iErrGet=SetNewEos(iEosOpt) !wipe out any instance of other eos.
		if(iEosOpt==1)CALL GetPR(NC,iErrGet)
		if(iEosOpt==2)CALL GetEsdCas(NC,idCas,iErrGet)
		if(iEosOpt==3)CALL GetPRWS(NC,iErrGet)
		if(iEosOpt==5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
		if(iEosOpt==6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
		if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
		if(iEosOpt==10)CALL GetPcSaft(NC,idCas,iErrGet)
		if(iEosOpt==11)CALL GetPrtc(NC,iErrGet)
		if(iEosOpt==17)CALL GetPrLorraine(NC,iErrGet)
		if(iErrGet > 0 .or. iErrCrit > 0 .or. ierLookup > 0)then
			write(*,*)'KijDB Error: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(61,607)id(1),id(2),0,zero,eightySix
			cycle
		endif
		iPt=0
		do i=1,nLines ! parse the lines into data
			if(var(i,3).ne.1)then
				print*,'var()=',(var(i,j),j=1,nVars)
				pause 'KijDbLLE: XDAT =/= X1???'
			endif
			if( iOptimize >0 .and. (var(i,4)<0.or.var(i,5)<0) )cycle ! for regression, we require coexisting x,y and subcritical.
			!print*,'KijDbLLE: x1Lo,x1Up=',var(i,4),var(i,5)
			iPt=iPt+1
			XDAT(iPt)=var(i,4)
			YDAT(iPt)=var(i,5)
			TDAT(iPt)=var(i,1)
			PDAT(iPt)=var(i,2)/1000 !convert from kPa to MPa.
			if(TDAT(iPt) > tMax)tMax=TDAT(iPt)
			if(TDAT(iPt) < tMin)tMin=TDAT(iPt)
		enddo
		nPtsBipDat=iPt
		if(nPtsBipDat==0)cycle !omit systems for which no coexistence points available.
		! All DAT read in for this dataset. Time to regress.
		!write(*,*)'KijDbLLE: New mix, NPTS,ID(),TLast=',nPtsBipDat,id(1),id(2),TDAT(nPtsBipDat)
		nSystems=nSystems+1
		write(6 ,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax !,Tc(1),Tc(2)
		!write(61,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,Tc(1),Tc(2)
		!write(661,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,Tc(1),Tc(2)
		if(LOUDER)then
			do i=1,nPtsBipDat
				write(6,'(f8.2,f8.4,2e12.4)')TDAT(i),PDAT(i),XDAT(i),YDAT(i)
			enddo
		endif
					
		write(*,*)'New mixture. iSystem,id(),nPts=',nSystems,id(1),id(2),nPtsBipDat
!		if(nPtsBipDat.le.0)then
!			ier=1 !indicates end of database
!			cycle !loop to end will terminate
!		endif
!		if(idOpt==2)call IdCcLookup(NC,idcc,ierLookup,errMsgLookup)
!		if(ierLookup.ne.0)then
!			write(61,607)idcc(1),idcc(2),zero,eightySix
!			DO I=1,NPTS	!read the points but do nothing so we can cycle to next.
!				READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
!			enddo
!			cycle !loop to next component but leave ier alone so keep looping.
!		endif
		parm=0
		nParms=QueryNParMix() ! Fortran requires declaration of integer QueryNParMix above.
		parm(2)=0
		temParm=0
		do i=1,nParms
			call QueryParMix(i,temParm(i),iErrParMix)  
			parm(i)=temParm(i) ! this will take initial guess from database, if available.
		enddo
		print*,'Stored BIPs=',(parm(k),k=1,nParms)
		if(iErrParMix)print*,'KijDbLLE: error returning stored mix parms'
		write(*,'(a,3E12.4)')' KijDbLLE: IniParms=',(temParm(i),i=1,nParms)
		if(iEosOpt==86)then	! 86 means skip this for now. Change 86 to 5 for evaluating solvation effects.
			nParms=2 ! add one to adjust bipDA
			write(*,'(5x,11i7)')(idLocalType(i),i=1,nTypesTot)
			do i=1,nTypesTot
				write(*,'(i5,11f7.4)')idLocalType(i),(aBipDA(i,j),j=1,nTypesTot)
			enddo
			print*,'Enter i,j for site type optimization'
			read*,iDA,jDA
			parm(2)=aBipDA(iDA,jDA)
		elseif(iEosOpt==4)then
			nParms=2
			parm(2)= -0.5
		endif
		if(iOptimize >0)then
			MAXIT=111
			INIT=1
			LWA=275	 !C	LWA>5*nParms+NPTS ~ 275
			TOL=0.0001	 !c      TOL=SQRT(DPMPAR(1))
			factor= 0.2 !c  factor=100 enables large first step from the initial guess.
			if(iEosOpt==17)factor=0.2
			if(LOUDER)write(*,*)'KIJ,paaLDK,rmserr'
			if(ABS(temParm(1)) < zeroTol)then	! signal that parameters do not exist for this compound
				bestKij= 0.1
				bestP2 = 0.1
				bestAAD=12345
				do iGrid=1,4
					parm(1)= (iGrid-1)*0.05d0-0.125
					if(iEosOpt==17)then
						parm(1)=parm(1)*10000
						parm(2)=parm(1)
					endif
					call KIJDevLLe(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,0)	!nConverged passed by BIPs module. 0 means enforce penalty for nonconverged points.
					if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
					write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAALDK,rmsErr,DstdErr,name(1),name(2)
					if(PAALDK < bestAAD)then
						bestAAD=PAALDK
						bestKij=parm(1)
						if(iEosOpt==17)bestP2=parm(2)
					endif
				enddo
				parm(1)=bestKij
				parm(2)=bestP2
			endif !ParmSrch
			print*,'Initial kij()=',(parm(k),k=1,nParms)
			if(LOUDER)pause 'Ready to regress'
			nCoexValues=2*nPtsBipDat ! Coexisting phases have 2 values each.
			CALL LMDifEZ(LMDevLLE,nCoexValues,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
			if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
			if(LOUDER)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		elseif(ABS(parm(1)) < zeroTol)then
			write(61,*)'KijDbLLE: warning, parm(1)=0'
		endif ! iOptimize
		if(iEosOpt.le.3.or.iEosOpt.eq.7)then
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ','       ',' xsBIP12',' xsBIP21'
		elseif(iEosOpt.le.5.or. iEosOpt.eq.7.or.iEosOpt.eq.9)then	   !AFG 2011 , Added EOS Opt. 9
			if(LOUDER)write(* ,*)'   KIJ','      ',' HIJ','       ','  kTij'
		else
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ' !,'      ',' xsBIP12',' xsBIP21'
		endif
		!write(* ,'(5f9.4)')(parm(i),i=1,nParms)
		if(LOUDER)write(* ,*)'StdErr'
		if(LOUDER)write(* ,'(5f9.4)')(stdErr(i),i=1,nParms)
		!if(LOUDER)pause 'KijDbLLE: Converged. Last call to KijVal.'
		LastCall=1	! Do not apply penalty for unconverged points on LastCall. Just get the deviations for converged points.
		call KIJDevLLE(NC,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LastCall)	!nConverged passed by BIPs module.
		!CALL LMDevLLE(Mdata,Nparms,parm,deviate,IFLAG)	! This function is different from KijDevLLE()
		if(nConverged > 0)rmsErr=SQRT( DFLOAT(nPtsBipDat)/nConverged )*RmsJre(nPtsBipDat,deviate) !adjust for nConverged.nPts needed to sum over all elements.
		!write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==1)write(6 ,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==1)write(61,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(61,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(6 ,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		write(6 ,'( 3i5,2F10.2,11(1x,E11.4) )')nSystems,id(1),id(2),PAALDK,PBIASK,(parm(i),i=1,nParms)
	enddo !while(ier.eq.0) loop over all data sets in db
	write(61,*)'nSystem total = ',nSystems
	write( * ,*)' #systems=',nSystems
	!write(661,*)' #systems=',nSystems
	close(51)
	close(61)
	!close(661)
	write(*,*)'Pointwise details in: ',TRIM(outFilePts)
	write(*,*)'Kij Summary in: ',TRIM(outFile)
	pause 'Success! Check results in output files.'
	!restore stored nc,id
	NC=ncStore
	do i=1,NC
		id(i)=idStore(i)
	enddo
	CALL GETCRIT(NC,iErrCrit)
	call IdCasLookup(NC,idCas,ier,errMsgLookup)	! idDippr comes from GlobConst. This also resets the class(). .
	if(iErrCrit.ne.0)then
		if(LOUDER)write(*,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		if(LOUDER)pause
		stop
	endif 
	if(iErrGet.gt.0)then
		if(LOUDER)write(*,*)'KijDB Error: failed to restore required pure props.'
		if(LOUDER)write(*,*)'id1,id2',id(1),id(2)
		if(LOUDER)pause
	endif
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
607	FORMAT(1X,3i6,1X,<nParms>E11.4,1X,<nParms>E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
608	FORMAT(1X,3i6,1X,f10.5,1X,f10.5,f10.2,f10.2,f10.4,1x,a20,1x,a20)
609	FORMAT(1X,3i6,1X,f10.1,1X,f10.1,f10.2,f10.2,f10.4,1x,a20,1x,a20)
      
	RETURN
	END	   ! KijDbLLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJDBSLE(NC)
	!C
	!C  PURPOSE:  FUGI CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED USEing BIPs
	!C  PROGRAMMED BY:  JRE 12/21
	!C  METHOD:   Levenberg-Marquardt 
	!C  REFERENCE:NUMERICAL RECIPES, MinPack
	!C
	USE GlobConst
	USE BIPs  !includes nConverged,maxPts,tDat,...
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	!USE MSIMSL !include IMSL for sort routine: DSROWR
	!USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	parameter(nVarsMax=11)
	CHARACTER*44 FN
	character*77 errMsgLookup
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5),temParm(5),var(maxPts,nVarsMax) !,TR(NMX)
	Integer idStore(NMX),QueryNParMix !,iErrGet,idcc(2),idCas(NMX)
	LOGICAL LOUDER
	EXTERNAL LMDevSLE
	character*77 inFile,outFile,outFilePts,errMsgPas,dumString
	character*3 nexString
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/FloryWert/ vLiq(nmx)
	!common/bipDA/iDA,jDA
	common/DevPas/iOptimize
    data initial/1/
	LOUDER=LOUD	!This provides local control if desired.
	!LOUDER=.TRUE. !Uncomment this line for screen writes from this routine.
	!store current values of nc,id.
	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	NC=2
	zero=0
	eightySix= -86
    idOpt=2
	write(*,*)'Enter 1 for dippr ids or 2 for idCas'
	write(*,*)'(SLeDb.txt uses idCas)'
	iOptimize=1
	print*,'Enter 2 to optimize BIPs using SLEFL or 1 for Q&D optimization.'
	print*,'Enter 0 to use stored BIPs or -1 to print stored BIPs.'
	read*,iOptimize
!	read(*,*)idOpt
!	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in VleData subdir)'
!	WRITE(*,*)'FIRST LINE =NPTS,id1,id2 2ND LINE = TK, PMPA, X,Y, ...'
!	write(*,*)'Enter 1 for Jaubert, 2 for JaubertAq, 3 for DannerGess, +10 to skip opt'
!	read(*,*)iAns
!	if(iAns > 10)then
!		iOptimize=0
!		iAns=iAns - 10
!	endif
!	if(iAns==2)FN="JaubertAq.txt"
!	if(iAns==4)FN="JaubertAbb.txt"
!	if(iAns==3)FN="DannerGessVle.txt"
!	READ(*,'(A)')FN
    FN="SLeDbPGL6ed20220129.txt"
	inFile=TRIM(masterDir)//'\LLe&SLeData\'//TRIM(FN)
	if(DEBUG)inFile='c:\spead\calceos\LLe&SLeData\'//TRIM(FN)
	print*,'inFile=',TRIM(inFile)
	inUnit=5151
	open(inUnit,file=inFile)

	outFilePts=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFilePts='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFilePts)
	if(iOptimize==1)write(661,*)'   T(K)      x(1)     x1DB       y1DB    dev1     itMax ' ! purge old contents of KijOut. KijVal uses APPEND
	if(iOptimize==1)write(661,*)'   T(K)      y1DB    x(1)     x1DB       dev1     gam1    itMax ' ! purge old contents of KijOut. KijVal uses APPEND
	if(iOptimize==0)write(661,*)'   T(K)      x(1)     x1DB       y1DB    dev1     itMax ' ! purge old contents of KijOut. KijVal uses APPEND
	close(661) ! purge old contents of KijOut. KijVal uses APPEND.

	outFile=TRIM(masterDir)//'\output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\output\KijDb.txt'
	iSumUnit=61
	if(initial)then
        initial=0
        open(iSumUnit,file=outFile)
    else
        open(iSumUnit,file=outFile,ACCESS='APPEND') ! e.g. JaubertAq is appended to Jaubert if during the same run.
    endif
	if(iOptimize >1)write(iSumUnit,*)'   ID1   ID2  nConv      Kij       stdErr   %PAALDK    rmsErrP    dStdErr'
	!if(iOptimize==0)write(iSumUnit,*)'   T(K)      x(1)     x1DB       y1     y1DB    dev1     dev2 ' ! purge old contents of KijOut. KijVal uses APPEND
	if(iOptimize <0)write(iSumUnit,*)' nPts  id1  id2         idCas1         idCas2   Tmin    Tmax     kij  ' ! purge old contents of KijOut. KijVal uses APPEND

	if(iEosOpt==4)then
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij'
	elseif(iEosOpt==5.or.iEosOpt==8.or.iEosOpt==9)then	   !AFG 2011 , Added EOS Opt. 9
		write(*,*)'Options:  1-kij, 2-kij,Hij, 3-kij,Hij,kTij, 4-kij,kTij'
	elseif(iEosOpt==3)then
		write(*,*)'Options:  1-kij, 2-tauIJ,tauJI, 3-kij,tauIJ,tauJI, '
		write(*,*)'          4-kij,kTij,tauIJ,tauJI, 5-kij,kTij,tauIJ,tauJI,alphaIJ'
	elseif(iEosOpt.eq.7)then
		write(*,*)'Options:  2-tauIJ,tauJI, 3-alphaIJ,tauIJ,tauJI'
	else
		write(*,*)'Options:  1-kij, 2-kij,kTij '
    endif
    iRegOpt=1
!	READ(5,*)iRegOpt
	ier=0
    iSystem=1
	nSystems=0 ! only count systems that have coex data.
	nLines=0
	nexString="***"
	nVars=4
	read(inUnit,'(a77)',ioStat=ioErr)dumString ! first line is *** Next *** so just skip it.
	do while(ierRead==0) !loop over data sets in db
		call ReadNext(inUnit,nexString,maxPts,nVarsMax,.TRUE.,dumString,nVars,var,nLines,ierRead)
		print*,'Header= ',TRIM(dumString)
		if(ierRead == -2)exit ! end of file	before reading any lines. ierRead= -1 for end of file at end of dataset.
		read(dumString,*,ioStat=ioErr)i99,i13,idCas(1),idCas(2)
		iColTK=1 ! sort by iDat component instead of temperature 
		Call SortTable(var,maxPts,nLines,nVars,iColTK,1) ! 1=yes for LoToHi
		if(ioErr)pause 'KijDbSLE: error reading header of dataset?'
		tMin=1234
		tMax=0
		call IdDipprLookup(NC,idCas,ierLookup,errMsgLookup)	! idDippr goes to GlobConst. This also sets the class()
		if(ierLookup > 0)then
			print*,'KijDbSLE: ID fail!',idCas(1),idCas(2),ID(1),ID(2)
			print*,TRIM(errMsgLookup)
			cycle
		endif
		!print*,'ID=',ID(1),ID(2)
		!call IdCasLookup(NC,idCas,ierLookup,errMsgLookup)	
		!if(LOUDER)pause 'Starting new binary mixture'
		CALL GETCRIT(NC,iErrCrit)
		if(iErrCrit.ne.0)then
			print*,'ID=',ID(1),ID(2)
			write(*,*)'KijDbSLE Error: iErrCrit,Tc() = ',iErrCrit,Tc(1),Tc(2)
			!pause
			cycle
		endif
		iErrGet=SetNewEos(iEosOpt) !wipe out any instance of other eos.
		if(iEosOpt==1)CALL GetPR(NC,iErrGet)
		if(iEosOpt==2)CALL GetEsdCas(NC,idCas,iErrGet)
		if(iEosOpt==3)CALL GetPRWS(NC,iErrGet)
		if(iEosOpt==5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
		if(iEosOpt==6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
		if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
		if(iEosOpt==10)CALL GetPcSaft(NC,idCas,iErrGet)
		if(iEosOpt==11)CALL GetPrtc(NC,iErrGet)
		if(iEosOpt==17)CALL GetPrLorraine(NC,iErrGet)
		if(iErrGet > 0 .or. iErrCrit > 0 .or. ierLookup > 0)then
			write(*,*)'KijDB Error: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(iSumUnit,'(2i6,3x,2i15,f12.2)')id(1),id(2),idCas(1),idCas(2),eightySix
			cycle
		endif
		iPt=0
		iDatOld=0
		do i=1,nLines ! parse the lines into data
			!if( iOptimize==1 .and. (var(i,4)<0.or.var(i,5)<0) )cycle ! for regression, we require coexisting x,y and subcritical.
			!print*,'KijDbLLE: x1Lo,x1Up=',var(i,4),var(i,5)
			if( var(i,4) > 1-zeroTol)cycle !omit points that are pure
			if( var(i,4) <   zeroTol)cycle !omit points that are pure. Also omits x1<0, which denotes bad data.
			iPt=iPt+1
			XDAT(iPt)=var(i,4)
			!YDAT(iPt)=var(i,5)	! the solid phase is assumed pure.
			TDAT(iPt)=var(i,1)
			PDAT(iPt)=var(i,2)/1000  !convert from kPa to MPa. 
			iDat(iPt)=NINT(var(i,3))!easier for general read as real, but we need integer for indexing.  
			if(TDAT(iPt) > tMax)tMax=TDAT(iPt)
			if(TDAT(iPt) < tMin)tMin=TDAT(iPt)
			CALL GetTmHfus(id(iDat(iPt)),TmKelvin,HfusJ_mol,iErrHfus)
			if( id(1)==1319)then
				print*,'KijDbSle: ID=1319, Hfus,iErrHfus=',HfusJ_mol,iErrHfus
				!pause 'check iErr for 1319'
			endif
			if(iErrHfus.ne.0 .or. TmKelvin < zerotol)then
				print*,'ID(iSolute)=',ID(iDat(iPt))
				iPt=iPt-1 ! disqualify last iPt+1 so nPtsBipDat=0 if nothing.
				write(*,*)'KijDbSLE Error: iErrHfus,Tm,Hfus = ',iErrHfus,TmKelvin,HfusJ_mol
				!pause
				cycle !
			endif
			YDAT(iPt)=exp( -HfusJ_mol/Rgas*(1/TDAT(iPt)-1/TmKelvin) )	 ! use YDAT to pass ideal soly.
		enddo
		nPtsBipDat=iPt
		if(nPtsBipDat==0)cycle !omit systems for which no coexistence points available.
		! All DAT read in for this dataset. Time to regress.
		!write(*,*)'KijDbLLE: New mix, NPTS,ID(),TLast=',nPtsBipDat,id(1),id(2),TDAT(nPtsBipDat)
		nSystems=nSystems+1
					
!		write(*,*)'New mixture. iSystem,id(),nPts=',nSystems,id(1),id(2),nPtsBipDat
!		if(nPtsBipDat.le.0)then
!			ier=1 !indicates end of database
!			cycle !loop to end will terminate
!		endif
!		if(idOpt==2)call IdCcLookup(NC,idcc,ierLookup,errMsgLookup)
!		if(ierLookup.ne.0)then
!			write(61,607)idcc(1),idcc(2),zero,eightySix
!			DO I=1,NPTS	!read the points but do nothing so we can cycle to next.
!				READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
!			enddo
!			cycle !loop to next component but leave ier alone so keep looping.
!		endif
		parm=0
		nParms=QueryNParMix() ! Fortran requires declaration of integer QueryNParMix above.
		parm(2)=0
		temParm=0
		do i=1,nParms
			call QueryParMix(i,temParm(i),iErrParMix)  
			parm(i)=temParm(i) ! this will take initial guess from database, if available.
		enddo
		if(nParms < 2 .and. iOptimize > 2)then
			nParms=2
			parm(1)=parm(1)*10000
			parm(2)=parm(1)	 ! simple initialization to check if initial guess makes a difference.
		endif
		tMin=TDAT(1) ! take advantage of sorted simplicity.
		tMax=TDAT(nPtsBipDat)
		write(6 ,'(3i5,2i15,2f8.2,f8.4)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,parm(1) !,Tc(1),Tc(2)
		!write(iSumUnit,'(3i5,2i15,2f8.2,f8.5)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,parm(1)
		!write(661,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,Tc(1),Tc(2)
		if(iOptimize < 0)cycle ! just printing summary of db and relevant BIPs.
		if(LOUDER)then
			do i=1,nPtsBipDat
				write(6,'(f8.2,f8.4,2e12.4)')TDAT(i),PDAT(i),XDAT(i),YDAT(i)
			enddo
		endif
		print*,'Stored BIPs=',(parm(k),k=1,nParms)
		if(iErrParMix)print*,'KijDbLLE: error returning stored mix parms'
		write(*,'(a,3E12.4)')' KijDbLLE: IniParms=',(temParm(i),i=1,nParms)
		if(iEosOpt==86)then	! 86 means skip this for now. Change 86 to 5 for evaluating solvation effects.
			nParms=2 ! add one to adjust bipDA
			write(*,'(5x,11i7)')(idLocalType(i),i=1,nTypesTot)
			do i=1,nTypesTot
				write(*,'(i5,11f7.4)')idLocalType(i),(aBipDA(i,j),j=1,nTypesTot)
			enddo
			print*,'Enter i,j for site type optimization'
			read*,iDA,jDA
			parm(2)=aBipDA(iDA,jDA)
		elseif(iEosOpt==4)then
			nParms=2
			parm(2)= -0.5
		endif
		if(iOptimize>0)then
			MAXIT=111
			INIT=1
			LWA=275	 !C	LWA>5*nParms+NPTS ~ 275
			TOL=0.0001	 !c      TOL=SQRT(DPMPAR(1))
			factor= 2 !c  factor=100 enables large first step from the initial guess.
			if(nParms > 1)factor=0.2
			if(iOptimize > 2)nParms=2 ! for wilson and NRTL models 
			if(LOUDER)write(*,*)'KIJ,paaLDK,rmserr'
			if(ABS(temParm(1)) < zeroTol .or. iOptimize > 2)then	! signal that parameters do not exist for this compound
				bestKij= 0.1
				bestP2 = 0.1
				bestAAD=12345
				do iGrid=1,4
					parm(1)= (iGrid-1)*0.05d0 -0.025 ! start a little off zero to help Hessian search.
					if(nParms > 1)then
						parm(1)=parm(1)*10000
						parm(2)=parm(1)
					endif
					call KIJDevSLe(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,0)	!nConverged passed by BIPs module. 0 means enforce penalty for nonconverged points.
					if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
					write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAALDK,rmsErr,DstdErr,name(1),name(2)
					if(PAALDK < bestAAD)then
						bestAAD=PAALDK
						bestKij=parm(1)
						if(nParms > 1)bestP2=parm(2)
					endif
				enddo
				parm(1)=bestKij
				parm(2)=bestP2
			endif !ParmSrch
			print*,'Initial kij()=',(parm(k),k=1,nParms)
			if(LOUDER)pause 'Ready to regress'
			nCoexValues=2*nPtsBipDat ! Coexisting phases have 2 values each.
			CALL LMDifEZ(LMDevSLE,nPtsBipDat,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
			if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
			if(LOUDER)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		elseif(ABS(parm(1)) < zeroTol)then
			write(iSumUnit,*)'KijDbSLE: warning, parm(1)=0'
		endif ! iOptimize
		if(iEosOpt.le.3.or.iEosOpt.eq.7)then
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ','       ',' xsBIP12',' xsBIP21'
		elseif(iEosOpt.le.5.or. iEosOpt.eq.7.or.iEosOpt.eq.9)then	   !AFG 2011 , Added EOS Opt. 9
			if(LOUDER)write(* ,*)'   KIJ','      ',' HIJ','       ','  kTij'
		else
			if(LOUDER)write(* ,*)'   KIJ','      ',' KTIJ' !,'      ',' xsBIP12',' xsBIP21'
		endif
		!write(* ,'(5f9.4)')(parm(i),i=1,nParms)
		if(LOUDER)write(* ,*)'StdErr'
		if(LOUDER)write(* ,'(5f9.4)')(stdErr(i),i=1,nParms)
		!if(LOUDER)pause 'KijDbLLE: Converged. Last call to KijVal.'
		LastCall=1	! Do not apply penalty for unconverged points on LastCall. Just get the deviations for converged points.
		call KIJDevSLE(NC,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LastCall)	!nConverged passed by BIPs module.
		!CALL LMDevLLE(Mdata,Nparms,parm,deviate,IFLAG)	! This function is different from KijDevLLE()
		if(nConverged > 0)rmsErr=SQRT( DFLOAT(nPtsBipDat)/nConverged )*RmsJre(nPtsBipDat,deviate) !adjust for nConverged.nPts needed to sum over all elements.
		!write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==1)write(6 ,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==1)write(iSumUnit,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(iSumUnit,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(6 ,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		write(6 ,'( 3i5,2F10.2,11(1x,E11.4) )')nSystems,id(1),id(2),PAALDK,PBIASK,(parm(i),i=1,nParms)
	enddo !while(ier.eq.0) loop over all data sets in db
	write(61,*)'nSystem total = ',nSystems
	write( * ,*)' #systems=',nSystems
	!write(661,*)' #systems=',nSystems
	close(51)
	close(61)
	!close(661)
	write(*,*)'Pointwise details in: ',TRIM(outFilePts)
	write(*,*)'Kij Summary in: ',TRIM(outFile)
	pause 'Success! Check results in output files.'
	!restore stored nc,id
	NC=ncStore
	do i=1,NC
		id(i)=idStore(i)
	enddo
	CALL GETCRIT(NC,iErrCrit)
	call IdCasLookup(NC,idCas,ier,errMsgLookup)	! idDippr comes from GlobConst. This also resets the class(). .
	if(iErrCrit.ne.0)then
		if(LOUDER)write(*,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		if(LOUDER)pause
		stop
	endif 
	if(iErrGet.gt.0)then
		if(LOUDER)write(*,*)'KijDB Error: failed to restore required pure props.'
		if(LOUDER)write(*,*)'id1,id2',id(1),id(2)
		if(LOUDER)pause
	endif
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
607	FORMAT(1X,3i6,1X,<nParms>E11.4,1X,<nParms>E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
608	FORMAT(1X,3i6,1X,f10.5,1X,f10.5,f10.2,f10.2,f10.4,1x,a20,1x,a20)
609	FORMAT(1X,3i6,1X,f10.1,1X,f10.1,f10.2,f10.2,f10.4,1x,a20,1x,a20)
      
	RETURN
	END	   ! KijDbLLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE KIJOPT(NC)
	!C
	!C  PURPOSE:  BP CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
	!C  PROGRAMMED BY:  JRE 2/96
	!C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
	!C  REFERENCE:NUMERICAL RECIPES
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	DIMENSION deviate(maxPts),parm(1)
	character*251 inFile !outFile*251,
	COMMON/eta/etaL,etaV,ZL,ZV

	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (..\calceos\VleData\FN.FT)'
	WRITE(*,*)'FIRST LINE =NPTS, 2ND LINE = TK, PMPA, X, ...'
	FN='Jaubert.txt'
	READ(*,'(A)')FN
	!WRITE(52,*)' NAME OF FILE WITH EXPERIMENTAL DATA (PATH\FN.FT)'
	!write(52,'(1x,A)')FN
	inFile=TRIM(masterDir)//'\VleData\'//TRIM(FN)
	if(DEBUG)inFile='c:\spead\calceos\VleData\'//TRIM(FN)
	write(*,*)'KijOpt: Data dir=',TRIM(inFile)
	open(51,file=inFile)
	READ(51,*)nPtsBipDat
	pDatMin=1E11
	pDatMax=0
	DO I=1,nPtsBipDat
		READ(51,*,ERR=861)tDat(I),PDAT(I),XDAT(I),YDAT(I)
		WRITE(*,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
		if(PDAT(i) > pDatMax)pDatMax=PDAT(i)
		if(PDAT(i) < pDatMin)pDatMin=PDAT(i)
	enddo
	MAXIT=111
	INIT=1
	WRITE(*,*)'THIS ROUTINE ASSUMES A BINARY SYSTEM FROM MAIN PROGRAM'
	WRITE(6,*)'ENTER RANGE TO OPTIMIZE, (LO, HI)  '
	READ(5,*)KIJ0,KIJ3                                          
	!WRITE(52,*)'ENTER RANGE TO OPTIMIZE, (KIJLO KIJHI)  '
	!write(52,*)KIJ0,KIJ3                                          
	KIJ1=KIJ0+CGOLD*(KIJ3-KIJ0)
	KIJ2=KIJ0+RGOLD*(KIJ3-KIJ0)
	nParms=1
	parm(1)=KIJ0
	LastCall= -1 ! -1 means that this is the initial call so we can specially initialize if necessary
	

	call KIJVAL(nc,nParms,parm,deviate,PAAD0,DstdErr,LastCall)
	rms0=RmsJre(nptsBipDat,deviate)
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD0,rms0,(parm(i),i=1,nParms)
	LastCall=0
	parm(1)=KIJ1
	call KIJVAL(nc,nParms,parm,deviate,PAAD1,DstdErr,LastCall)
	rms1=RmsJre(nptsBipDat,deviate)
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD1,rms1,(parm(i),i=1,nParms)
	parm(1)=KIJ2
	call KIJVAL(nc,nParms,parm,deviate,PAAD2,DstdErr,LastCall)
	rms2=RmsJre(nptsBipDat,deviate)
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD2,rms2,(parm(i),i=1,nParms)
	parm(1)=KIJ3
	call KIJVAL(nc,nParms,parm,deviate,PAAD3,DstdErr,LastCall)
	rms3=RmsJre(nptsBipDat,deviate)
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD3,rms3,(parm(i),i=1,nParms)
	TOL=0.0001 ! 4 DECIMAL PLACES ON KIJ IS ENOUGH.
	iter=0      
	DO WHILE( ABS(KIJ3-KIJ0) > TOL ) !*(ABS(KIJ1)+ABS(KIJ2)))
		iter=iter+1
		if(iter > maxit)exit
		!IF(rms2 < rms1)THEN
		IF(PAAD2 < PAAD1)THEN
			KIJ0=KIJ1
			KIJ1=KIJ2
			KIJ2=RGOLD*KIJ1+CGOLD*KIJ3
			PAAD0=PAAD1
			PAAD1=PAAD2
			rms0=rms1
			rms1=rms2
			parm(1)=KIJ2
			call KIJVAL(nc,nParms,parm,deviate,PAAD2,DstdErr,LastCall)
			rms2=RmsJre(nptsBipDat,deviate)
			write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD2,rms2,(parm(i),i=1,nParms)
		ELSE
			KIJ3=KIJ2
			KIJ2=KIJ1
			KIJ1=RGOLD*KIJ2+CGOLD*KIJ0
			PAAD3=PAAD2
			PAAD2=PAAD1
			rms3=rms2
			rms2=rms1
			parm(1)=KIJ1
			call KIJVAL(nc,nParms,parm,deviate,PAAD1,DstdErr,LastCall)
			rms1=RmsJre(nptsBipDat,deviate)
			write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAAD1,rms1,(parm(i),i=1,nParms)
		ENDIF
	ENDDO
	!IF(rms1 < rms2)THEN
	IF(PAAD1 < PAAD2)THEN
		PAADP=PAAD1
		rmse=rms1
		KIJB=KIJ1
	ELSE
		PAADP=PAAD2
		rmse=rms2
		KIJB=KIJ2
	ENDIF    
	!open(67,file=outFile) ! KijVal writes KijOut using append.
	!write(67,*)' ' !clear out old contents. Don't do this if you want to run ko after kd, but retain KijOut of big database. 
	!close(67)
	!WRITE(52,*)'KIJOPT, PAADP'
	!WRITE(52,606)KIJB,PAADP
	write(*,600)

	LastCall=1
	parm(1)=KIJb
	call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,LastCall) ! On LastCall, penalties set to zero and prints to KijOut.txt
	rmse=RmsJre(nPtsBipDat,deviate)
	write(*,*)'KIJOPT, PAADP'
	write(*,606)KIJB,PAADP
	!outFile=TRIM(masterDir)//'\output\KijOut.txt'
	!IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	!open(67,file=outFile, ACCESS='APPEND') !outFile is defined above
	!WRITE(67,*)'KIJOPT, PAADP'
	!WRITE(67,606)KIJB,PAADP
	!write(67,600)
	do i=1,nPtsBipDat
		!dev = 100*ln(calc/expt) => expt*exp( dev/100 ) = calc
		pcalc = pdat(i)*EXP( deviate(i)/100 )
		iErrDev=0
		if(ABS(deviate(i)) < 1.D-8)iErrDev=1
		PAAD = ( pCalc-pdat(i) )/pdat(i)*100
		!write(67,'(1x,f10.3,3f10.5,i5,f10.2)')tdat(i),pcalc,xdat(i),ydat(i),iErrDev,PAAD	
	enddo
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rms,parm()=',PAADP,rmse,(parm(i),i=1,nParms)

	!CLOSE(67)
	close(51)
	!write(*,*) 'KijOut.txt=',TRIM(outFile)
	pause 'Detailed results are tabulated in KijOut.txt.'
	RETURN
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)
861	if(LOUD)write(*,*)'Error reading line=',i,' from data file. Check format.'
	return      
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,LAST)
	USE GlobConst
	USE BIPs ! includes nConverged,Kij(),KTij(), ...
	USE VpDb ! for VpCoeffs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DIMENSION X(NMX),Y(NMX) !,yy(nmx),xx(nmx),
	integer ier(20)
	DOUBLE PRECISION parm(nParms),deviate(maxPts) !,temParm(nParms)
	character*77 errMsg(12),outFile
	DIMENSION GammaCalc(nmx),GammaExp(nmx),pSatExp(nmx)
	COMMON/eta/etaL,etaV,ZL,ZV
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	!dimension vpCoeffs(nmx,5)
	! DATA iOpt/1/		!JRE 20210324 iOpt=2 was previously to optmize Hij instead of Kij. I gave up on that.

	errMsg(1)='KijVal Error: Number of components must equal 2.'
	if(nc.ne.2)then
		if(LOUD)write(*,*)errMsg(1)
		ier=1
		return
	endif
                
607	FORMAT(1X,F8.4,1X,F10.4,1X,F10.8,1X,F10.8,1X,F10.4,1X,i5)      

	hijTmp=HIJ(1,2)
	KIJ(1,2)=0
	KTIJ(1,2)=0
	HIJ(1,2)=0
	HTIJ(1,2)=0
	xsTau(1,2)=0
	xsTau(2,1)=0
	xsAlpha(1,2)=0.3d0
	xsAlpha(2,1)=xsAlpha(1,2)
	do i=1,nParms
		call SetParMix(i,parm(i),iErrSetBip)
	enddo
	!done initializing parm/bip arrays
	!temParm=0
	!do i=1,nParms
	!	call QueryParMix(i,temParm(i),iErrParMix)  !for debugging
	!enddo
	!if(iErrParMix)print*,'KijVal: error returning stored mix parms'
	!write(*,'(a,3E12.4)')' KijVal: Passed Parms=',(temParm(i),i=1,nParms)

	MAXIT=111
	PAADP=0
	ssqErr=0
	!GexCalc=0
	!GexExp=0
	DstdErr=0
	!ITMAX=MAXIT
	!INIT=1
	nConverged=0
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFile,ACCESS='APPEND')
	if(LAST==1)write(661,'(3I7, 2F9.4)')nPtsBipDat,ID(1),ID(2),Kij(1,2),KTij(1,2)
	pMax=5.1d0 !for most systems, anything higher than this just causes more crashes and skews the importance toward scf region
	if(pDatMin > pMax)pMax=pDatMin+0.1 ! for some systems (e.g. nC2+water) the lowest pressure is large. if nConverged=0, LmDif crashes.           
	DO iData=1,nPtsBipDat !	data in module BIPs
		X(1)=XDAT(iData)
		X(2)=1-X(1)
		!Compute self-consistent pSat to enable gamma calculations
		!This needs to be inside the data loop because the temperature may change.
		P=PDAT(iData)
		if(P <= 0)pause 'KijVal: PDAT(iData)=P <= 0???'
		!if(ID(1)==1921 .or. ID(2)==1921)P=pMax ! using a normal guess for P causes BubPL to get stuck in local minimum.
		T=tDat(iData)
		deviate(iData)=0
        if(P > pMax)then
			if(LOUD)print*,'KijVal: Skipping P > Pmax.',P,Pmax
			!No error indicator & no penalty, but nConverged is not incremented and deviate=0 can be "interpreted." 
			cycle  ! Higher pressures led to too many crashes.
		endif
        tMax=0
        do i=1,NC
            if(0.9*TC(i) > tMax)tMax=0.9d0*Tc(i)
        enddo
        if(T > tMax)then ! Near critical or scf temperatures cause crashes and other bad behavior.
			if(LOUD)print*,'KijVal: Skipping T > Tmax.',T,Tmax
			!No error indicator & no penalty, but nConverged is not incremented and deviate=0 can be "interpreted." 
			cycle  ! Higher pressures led to too many crashes.
		endif
		!P1=PDAT(iData)
		IFLAG=0
		ITMAX=MAXIT
		INIT=1    
		if(NC==1)then
	        call PsatEar( T,P,chemPo,rhoLiq,rhoVap,uSatL,uSatV,ierCode)  !using y() as dummy for chemPo()
        else
            CALL BUBPL( T,X,NC,INIT,P,ITMAX,Y,ier)
        endif
		if(ier(1).ne.0 .and. NC > 1)then
			!print*,'KijVal: BubPL failed. iData,Kij=',iData,Kij(1,2)
			iWrite=0
			if(LOUD)iWrite=1
			CALL BootBPx(T,X,NC,INIT,P,ITMAX,Y,ier,iWrite)
		endif

		IF(ier(1).ne.0 .and. LOUD)WRITE(*,*) 'Error in POINT NUMBER = ', iData

		if(ier(1) == 0 .and. P > 0)then
			deviate(iData) = DLOG( P/PDAT(iData) )*100 ! use log deviate to minimize bias for bad fits.
			pDev = (P-PDAT(iData))/PDAT(iData)*100
			nConverged=nConverged+1
			!print*,'got one! nConverged=',nConverged
			PAADP = PAADP+ABS(pDev)
			ssqErr= ssqErr+deviate(iData)*deviate(iData)
		else
			deviate(iData)=300 ! gas+water systems go to nConverged=0 if penalty=200
			if(LAST==1)deviate(iData)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
		endif
		!penalize the deviate passed to lmdif, but only compile paadp for converged pts
        !pSatCalc1=P
		!yy(1)=y(1)
		!yy(2)=1-yy(1)
		!xx(1)=x(1)
		!xx(2)=1-xx(1)
		!IF(XX(1).EQ.0)XX(1)=1E-11
		!IF(XX(2).EQ.0)XX(2)=1E-11
		!compute SEDG (Std Err of Danner Gess)
		!IF(PDAT(iData).GT.2.0)IER(1)=1
		IF(LAST==1 .and. NC > 86)then		 ! Get vapor pressures and activity coefficients. 86=>disabled
			x(2)=1e-11
			x(1)=1-x(2)
			IF(LOUD)pause 'KijVal: How did I get to gamma calcs???'
			!IFLAG=0
			ITMAX=MAXIT
			INIT=0
			if(iEosOpt > 1 )then 
				x(1)=1e-11
				x(2)=1-x(1)
				ITMAX=MAXIT
				INIT=0
				CALL BUBPL(T,X,NC,INIT,pSatCalc2,ITMAX,Y,ier) ! Get EOS vapor pressure
				X(1)=XDAT(iData)
				X(2)=1-X(1)
				if(x(1).eq.0)x(1)=1e-11
				if(x(2).eq.0)x(2)=1e-11
				Y(1)=YDAT(iData)
				Y(2)=1-Y(1)
				do iComp=1,NC
					if (iComp.eq.1) then
						pSatCalc=1
					Else
						pSatCalc=1
					endif
					GammaCalc(iComp)=y(iComp)*P/(x(iComp)*pSatCalc)
					!GexCalc=GexCalc+X(iComp)*DLOG(GammaCalc(iComp))
				
					CALL GetVp(NC,ID,iErrCode)
				
					pSatExp(iComp)=exp(vpCoeffs(iComp,1)+vpCoeffs(iComp,2)/T+vpCoeffs(iComp,3)*DLOG(T)+vpCoeffs(iComp,4)*T**vpCoeffs(iComp,5))/1000000
					if(LOUD)print*,' KijVal(afterGetVp): iComp,T(K),PsatExp',iComp,T,pSatExp(iComp)
					if(iData.eq.1.AND.XDAT(nPtsBipDat).eq.1)pSatExp(1)=PDAT(nPtsBipDat)
					if(iData.eq.nPtsBipDat.and.XDAT(1).eq.0)pSatExp(NC)=PDAT(1)
					if(tDat(1).eq.tDat(2).AND.XDAT(nPtsBipDat).eq.1)pSatExp(1)=PDAT(nPtsBipDat)
					if(tDat(1).eq.tDat(2).AND.XDAT(1).eq.0)pSatExp(2)=PDAT(1)
					gammaExp(iComp)= -86.86D0
					if(y(iComp) < 0)cycle
					if( ABS( x(iComp) ) < zeroTol .and. LOUD)pause 'KijVal: xi=0? ' 
					if( ABS( pSatExp(iComp) ) < zeroTol .and. LOUD)pause 'KijVal: PSatExp=0? ' 
					GammaExp(iComp)=Y(iComp)*PDAT(iData)/(X(iComp)*pSatExp(iComp))
					!GexExp=GexExp+X(iComp)*DLOG(GammaExp(iComp))
				enddo ! iComp=1,NC
				GexCalc=X(1)*DLOG(GammaCalc(1))+X(2)*DLOG(GammaCalc(2))
				if ( abs(GammaExp(1)) > zeroTol .and. abs(GammaExp(2)) > zeroTol)then
					GexExp=X(1)*DLOG(GammaExp(1))+X(2)*DLOG(GammaExp(2))
				else
					GexExp=0
                endif
                if(LOUD)then
				    if(nPtsBipDat.le.1)pause 'KijVal: nPts > 1 is required'
                end if
				DstdErr=DstdErr+(GexExp-GexCalc)*(GexExp-GexCalc)
			endif ! iEosOpt > 1

		endif !last=1 .and. compute activities
		if(LAST==1)WRITE(661,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1) ! write deviations on LastCall
		!WRITE(6 ,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)
	enddo
	close(661)
	if(nPtsBipDat < 1)pause 'KijVal: nPtsBipDat < 1??? '
	DstdErr=SQRT(DstdErr/(nPtsBipDat))
 	!rmsErr=SQRT( ENORM(nConverged,deviate))
	if(nConverged > 0)then
		rmsErr=SQRT(ssqErr/nConverged)
		PAADP=PAADP/nConverged
	else
		rmsErr=8686
		PAADP=8686
	endif
	write(*,'(a,i4,2f10.3,7f12.4)')' nConverged,parms,%AADP,rms%ErrP=',nConverged,PAADP,rmsErr,(parm(i),i=1,nparms)
	RETURN
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJDevLLE(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LAST)
	USE GlobConst
	USE BIPs ! includes nConverged,Kij(),KTij(), ...
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DIMENSION X(NMX),Y(NMX),calcK(NMX),zFeed(NMX) !,yy(nmx),xx(nmx),
	integer ier(20)
	DOUBLE PRECISION parm(nParms),deviate(maxPts),FUGCL(NMX),FUGCU(NMX) !,temParm(nParms)
	character*77 errMsg(12),outFile
	LOGICAL LOUDER
	common/bipDA/iDA,jDA
	common/DevPas/iOptimize
	!DIMENSION GammaCalc(nmx),GammaExp(nmx),pSatExp(nmx)
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	!dimension vpCoeffs(nmx,5)
	! DATA iOpt/1/		!JRE 20210324 iOpt=2 was previously to optmize Hij instead of Kij. I gave up on that.

	errMsg(1)='KijVal Error: Number of components must equal 2.'
	if(nc.ne.2)then
		if(LOUD)write(*,*)errMsg(1)
		ier=1
		return
	endif

	LOUDER=LOUD
	!LOUDER=.TRUE.
                
	hijTmp=HIJ(1,2)
	KIJ(1,2)=0
	KTIJ(1,2)=0
	HIJ(1,2)=0
	HTIJ(1,2)=0
	xsTau(1,2)=0
	xsTau(2,1)=0
	xsAlpha(1,2)=0.3d0
	xsAlpha(2,1)=xsAlpha(1,2)
	nParmSet=nParms
	if(iEosOpt==5.and.nParms > 1)then
		if(nParms==2)nParmSet=nParms-1
		aBipDA(iDA,jDA)=parm(2)	 !iDA,jDA sent by common
		if(parm(2) > 1)aBipDA(iDA,jDA)=1-zeroTol	 !iDA,jDA sent by common
		aBipAD(jDA,iDA)=aBipDA(iDA,jDA)
		aBipDA(jDA,iDA)=aBipDA(iDA,jDA)				 !make symmetric for now
		aBipAD(iDA,jDA)=aBipDA(iDA,jDA)
	endif
	if(iEosOpt==4)then
		if(nParms==2)nParmSet=nParms-1
		Hij(1,2)=parm(2)
		Hij(2,1)=Hij(1,2)
	endif
	! constrain bad guesses
	if(nParms > 1)then
		do i=1,nParms
			if(parm(i) < -700)parm(i)= -700
			if(parm(i) >25000)parm(i)= 25000
		enddo
	else
		if( parm(1) > 0.9)parm(1)=0.9d0
		if( parm(1) <  -1)parm(1)= -1
	endif
	do i=1,nParmSet
		call SetParMix(i,parm(i),iErrSetBip)
	enddo
	print*,'KijDevLLE:parms=',(parm(i),i=1,nParms)

	!done initializing parm/bip arrays
	!temParm=0
	!do i=1,nParms
	!	call QueryParMix(i,temParm(i),iErrParMix)  !for debugging
	!enddo
	!if(iErrParMix)print*,'KijVal: error returning stored mix parms'
	!write(*,'(a,3E12.4)')' KijVal: Passed Parms=',(temParm(i),i=1,nParms)
	MAXIT=111
	PAALDK=0
	PBIASK=0
	ssqErr=0
	!GexCalc=0
	!GexExp=0
	DstdErr=0
	!ITMAX=MAXIT
	!INIT=1
	nConverged=0
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFile,ACCESS='APPEND')
	if(LAST==1)write(661,'(3I5,2i11, 2F9.4)')nPtsBipDat,ID(1),ID(2),idCas(1),idCas(2),Kij(1,2),KTij(1,2)
	if(LOUDER)write( 6 ,'(3I7, 2F9.4)')nPtsBipDat,ID(1),ID(2),Kij(1,2),KTij(1,2)
	if(LAST==1.or.LOUDER)then
		if(iOptimize >0)print*,'tDat(i),LOG(CalcK1),LOG(exptK1),LOG(CalcK2),LOG(exptK2),dev1,dev2' 
		if(iOptimize==0)print*,'tDat(),x(1),xDat(),y(1),yDat(),dev1,dev2,itMax'
	endif 
	pMax=5.1d0 !for most systems, anything higher than this just causes more crashes and skews the importance toward scf region
	if(pDatMin > pMax)pMax=pDatMin+0.1 ! for some systems (e.g. nC2+water) the lowest pressure is large. if nConverged=0, LmDif crashes.           
	rmsErr=0
	DO iData=1,nPtsBipDat !	data in module BIPs
		!This needs to be inside the data loop because the temperature may change.
		P=PDAT(iData)
		if(P < zeroTol)P=0.1 ! P shouldn't matter for LLE.
		if(P <= 0.and. LOUD)print*, 'KijDevLLE: PDAT(iData)=P <= 0???'
		!if(ID(1)==1921 .or. ID(2)==1921)P=pMax ! using a normal guess for P causes BubPL to get stuck in local minimum.
		T=tDat(iData)
		if(iOptimize > 0)then
			X(1)=XDAT(iData)
			X(2)=1-X(1)
			Y(1)=YDAT(iData) ! Y IS THE UPPER PHASE COMPOSITION
			Y(2)=1-Y(1)
			deviate(iData)=0
			LIQ=1
			if(iOptimize==2)then 
				if(xDat(iData) > 0 .and. yDat(iData) > 0 )then	 ! .and. iData < 3
					calcK(1)=yDat(iData)/xDat(iData) ! use experimental value if coexistence available. 
					calcK(2)=(1-yDat(iData))/(1-xDat(iData)) ! for soluble compounds, use experimental value if possible. 
				else ! then yDat must be valid
					cycle ! only use coexistence data when optimizing.
				endif
				itMax=151
				Call LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode)
				if(iErrCode==0)then
					dev1=DLOG(  ( 0.5d0-ABS(0.5d0-x(1)) )/( 0.5d0-ABS(0.5d0-xDat(iData)) )  )*100  !automatically takes the smaller
					dev2=DLOG(  ( 0.5d0-ABS(0.5d0-y(1)) )/( 0.5d0-ABS(0.5d0-yDat(iData)) )  )*100
					nValues=nValues+2
					deviate(iData) = dev1 ! use log deviate to minimize bias for bad fits.
					deviate(iData+nPtsBipDat) = dev2 ! use log deviate to minimize bias for bad fits.
					nConverged=nConverged+2 ! Two dev's in above deviate. Consistent with one-sided evaluations when coexistence unavailable.
					!print*,'got one! nConverged=',nConverged
					PAALDK = PAALDK+ABS(dev1)+ABS(dev2)
					PBIASK = PBIASK+(dev1+dev2)
					ssqErr= ssqErr+dev1*dev1+dev2*dev2
				else
					calcK(1)=1
					calcK(2)=1
					deviate(iData)=150 ! gas+water systems go to nConverged=0 if penalty=200
					deviate(iData+nPtsBipDat)=150 ! gas+water systems go to nConverged=0 if penalty=200
					if(LAST==1)then
						deviate(iData)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
						deviate(iData+nPtsBipDat)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
					endif
				endif !last=1 .and. compute activities
			else !iOptimize=1
				exptK1=Y(1)/X(1)
				exptK2=Y(2)/X(2)
				if(LOUDER)print*,'Calling fugi'
				ierCalc=0
				CALL FUGI( T,P,X,NC,LIQ,FUGCL,ZL,ier ) !FUGC=ln(phi)
				if(ier(1) > 10)ierCalc=1
				CALL FUGI( T,P,Y,NC,LIQ,FUGCU,ZU,ier )
				if(ier(1) > 10)ierCalc=1
				IF(ierCalc.ne.0 .and. LOUDER)WRITE(*,*) 'KijDevLLE: Error in POINT NUMBER = ', iData,ier(1)
				if(ierCalc==0)calcK(1:NC)=exp( fugcL(1:NC)-fugcU(1:NC) ) !FUGC=ln(phi),
				deltaK=calcK(1)-calcK(2)
				if(ABS(deltaK) > zeroTol)then	 
					x1Lo=(1-calcK(2))/deltaK
					if(x1Lo < zeroTol .or. x1Lo > (1-zeroTol) )ierCalc=2
					x1Up=x1Lo*calcK(1)
					if(x1Lo < zeroTol .or. x1Lo > (1-zeroTol) )ierCalc=3
				else
					ierCalc=4
				endif
				if(ierCalc==0)then
					dev1=  DLOG(  calcK(1)/exptK1 )*100  ! ln(K1calc/K1expt)
					dev2=  DLOG(  calcK(2)/exptK2 )*100  ! ln(K1calc/K1expt)
					dev1= DLOG( x1Lo/(1-x1Lo) )-DLOG( x(1)/(1-x(1)) ) ! Vladimir's metric
					dev2= DLOG( x1Up/(1-x1Up) )-DLOG( y(1)/(1-y(1)) ) ! Vladimir's metric
					deviate(iData) = dev1 ! use log deviate to minimize bias for bad fits.
					deviate(iData+nPtsBipDat) = dev2 ! use log deviate to minimize bias for bad fits.
					nConverged=nConverged+2 ! Two dev's in above deviate. Consistent with one-sided evaluations when coexistence unavailable.
					!print*,'got one! nConverged=',nConverged
					PAALDK = PAALDK+ABS(dev1)+ABS(dev2)
					PBIASK = PBIASK+(dev1+dev2)
					ssqErr= ssqErr+deviate(iData)*deviate(iData)
				else
					calcK(1)=1
					calcK(2)=1
					deviate(iData)=150 ! gas+water systems go to nConverged=0 if penalty=200
					deviate(iData+nPtsBipDat)=150 ! gas+water systems go to nConverged=0 if penalty=200
					if(LAST==1)then
						deviate(iData)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
						deviate(iData+nPtsBipDat)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
					endif
				endif !last=1 .and. compute activities
				if(LAST==1)WRITE(661,617)tDat(iData),LOG(calcK(1)),LOG(exptK1),LOG(calcK(2)),LOG(exptK2),dev1,dev2 ! write deviations on LastCall
				if(LAST==1.and.LOUDER)WRITE( 6 ,607)tDat(iData),LOG(calcK(1)),LOG(exptK1),LOG(calcK(2)),LOG(exptK2),dev1,dev2 ! write deviations on LastCall
				!WRITE(6 ,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)
			endif    
		else ! iOptimize < 1 => compute LLE @ T if not optimizing.
			initialize=0  ! use the previous (stored) value of x,y if initialize=0.
			if(iData==1)then !If the first point in dataset, we need to start fresh. Dataset should be sorted from LoT to HiT. 
				initialize=1
				calcK(1)=0.001	 ! x1=(1-K2)/(K1-K2) => x1=small if K1=large.

				if(xDat(iData) > 0 .and. yDat(iData) > 0)then
					calcK(1)=yDat(iData)/xDat(iData) ! use experimental value if coexistence available. 
					calcK(2)=(1-yDat(iData))/(1-xDat(iData)) ! for soluble compounds, use experimental value if possible. 
				elseif(xDat(iData) > 0)then
					if(xDat(iData) < 0.5d0)calcK(1)=(1/xDat(iData)+5)/2 ! for soluble compounds, use experimental value if possible. 
					if(calcK(1)>0.02)calcK(1)=0.02 ! JRE 20220113 added to improve convergence of H2O+isobutyric acid.  Questionable.
					calcK(2)=1/calck(1)	! assume symmetric for initial guess.
				else ! then yDat must be valid
					if(yDat(iData) > 0.5d0)calcK(1)=(1/(1-yDat(iData))+5)/2 ! for soluble compounds, use experimental value if possible. 
					if(calcK(1)>0.02)calcK(1)=0.02 ! JRE 20220113 added to improve convergence of H2O+isobutyric acid.  Questionable.
					calcK(2)=1/calck(1)	! assume symmetric for initial guess.
				endif
			endif
			itMax=151
			Call LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode)
			if(iErrCode > 0)cycle
			nValues=0
			dev1=0
			dev2=0
			if(xDat(iData) > 0)then
				dev1=DLOG( x(1)/xDat(iData) )*100
				if( xDat(iData) > ABS(yDat(iData)) )dev1=DLOG( (1-x(1))/(1-xDat(iData)) )*100
				dev1=DLOG(  ( 0.5d0-ABS(0.5d0-x(1)) )/( 0.5d0-ABS(0.5d0-xDat(iData)) )  )*100  !automatically takes the smaller
				nValues=nValues+1
			endif
			if(yDat(iData) > 0)then
				dev2=DLOG( y(1)/yDat(iData) )*100
				if( xDat(iData) < ABS(yDat(iData)) )dev2=DLOG( (1-y(1))/(1-yDat(iData)) )*100
				dev2=DLOG(  ( 0.5d0-ABS(0.5d0-y(1)) )/( 0.5d0-ABS(0.5d0-yDat(iData)) )  )*100
				nValues=nValues+1
			endif
			deviate(iData) = ( ABS(dev1)+ABS(dev2) )/nValues ! use log deviate to minimize bias for bad fits.
			nConverged=nConverged+nValues
			PAALDK = PAALDK+deviate(iData)
			PBIASK = PBIASK+(dev1+dev2)
			ssqErr= ssqErr+deviate(iData)*deviate(iData)
			if(LAST==1)WRITE(661,607)tDat(iData),x(1),xDat(iData),y(1),yDat(iData),dev1,dev2,itMax ! write deviations on LastCall
			if(LOUDER )WRITE( 6 ,607)tDat(iData),x(1),xDat(iData),y(1),yDat(iData),dev1,dev2,itMax ! write deviations on LastCall
		endif
	enddo
	close(661)
607	FORMAT( 1X,F8.2,4(1X,F9.6),2(1x,f8.2),i4 )      
617	FORMAT( 1X,F8.2,4(1X,F9.4),2(1x,f8.2),i4 )      

	if(nPtsBipDat < 1)pause 'KijVal: nPtsBipDat < 1??? '
	DstdErr=SQRT(DstdErr/(nPtsBipDat))
 	!rmsErr=SQRT( ENORM(nConverged,deviate))
	if(nConverged > 0)then
		rmsErr=SQRT(ssqErr/nConverged)
		PAALDK=PAALDK/nConverged
		PBIASK=PBIASK/nConverged
	else
		rmsErr=8686
		PAALDK=8686
		PBIASK=8686
	endif
	write(*,'(a,i4,2f10.3,7f12.4)')' KijDevLLE:nConverged,%AALDK,rms%Err,parms=',nConverged,PAALDK,rmsErr,(parm(i),i=1,nparms)
	if(LOUDER)pause 'Check K-values'
	RETURN
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE KIJDevSLE(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LAST)
	USE GlobConst
	USE BIPs ! includes nConverged,Kij(),KTij(), ...
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DIMENSION X(NMX),XidSoln(NMX) !,calcK(NMX),Y(NMX),yy(nmx),xx(nmx),
	integer ier(20)
	DOUBLE PRECISION parm(nParms),deviate(maxPts),FUGCPure(NMX),actCo(NMX),FUGCL(NMX),xOld(NMX)!,temParm(nParms)
	character*77 errMsg(12),outFile
	LOGICAL LOUDER
	common/bipDA/iDA,jDA
	common/DevPas/iOptimize
	!DIMENSION GammaCalc(nmx),GammaExp(nmx),pSatExp(nmx)
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	!dimension vpCoeffs(nmx,5)
	! DATA iOpt/1/		!JRE 20210324 iOpt=2 was previously to optmize Hij instead of Kij. I gave up on that.

	errMsg(1)='KijVal Error: Number of components must equal 2.'
	if(nc.ne.2)then
		if(LOUD)write(*,*)errMsg(1)
		ier=1
		return
	endif

	LOUDER=LOUD
	!LOUDER=.TRUE.
                
	hijTmp=HIJ(1,2)
	KIJ(1,2)=0
	KTIJ(1,2)=0
	HIJ(1,2)=0
	HTIJ(1,2)=0
	xsTau(1,2)=0
	xsTau(2,1)=0
	xsAlpha(1,2)=0.3d0
	xsAlpha(2,1)=xsAlpha(1,2)
	nParmSet=nParms
	if(iEosOpt==5.and.nParms > 1)then
		if(nParms==2)nParmSet=nParms-1
		aBipDA(iDA,jDA)=parm(2)	 !iDA,jDA sent by common
		if(parm(2) > 1)aBipDA(iDA,jDA)=1-zeroTol	 !iDA,jDA sent by common
		aBipAD(jDA,iDA)=aBipDA(iDA,jDA)
		aBipDA(jDA,iDA)=aBipDA(iDA,jDA)				 !make symmetric for now
		aBipAD(iDA,jDA)=aBipDA(iDA,jDA)
	endif
	if(iEosOpt==4)then
		if(nParms==2)nParmSet=nParms-1
		Hij(1,2)=parm(2)
		Hij(2,1)=Hij(1,2)
	endif
	if(nParms > 1)then
		do i=1,nParms
			if(parm(i) < -700)parm(i)= -700
			if(parm(i) >25000)parm(i)= 25000
		enddo
	else
		if( parm(1) > 0.9)parm(1)=0.9d0
		if( parm(1) <  -1)parm(1)= -1
	endif
	do i=1,nParmSet
		call SetParMix(i,parm(i),iErrSetBip)
	enddo
	print*,'KijDevSLE: parms=',(parm(i),i=1,nParms)

	!done initializing parm/bip arrays
	!temParm=0
	!do i=1,nParms
	!	call QueryParMix(i,temParm(i),iErrParMix)  !for debugging
	!enddo
	!if(iErrParMix)print*,'KijVal: error returning stored mix parms'
	!write(*,'(a,3E12.4)')' KijVal: Passed Parms=',(temParm(i),i=1,nParms)
	MAXIT=111
	PAALDK=0
	PBIASK=0
	ssqErr=0
	!GexCalc=0
	!GexExp=0
	DstdErr=0
	!ITMAX=MAXIT
	!INIT=1
	nConverged=0
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFile,ACCESS='APPEND')
	if(LAST==1)write(661,'(3I5,2i11, 2F9.4)')nPtsBipDat,ID(1),ID(2),idCas(1),idCas(2),Kij(1,2),KTij(1,2)
	if(LOUDER)write( 6 ,'(3I7, 2F9.4)')nPtsBipDat,ID(1),ID(2),Kij(1,2),KTij(1,2)
	if(LAST==1.or.LOUDER)then
		if(iOptimize >0)print*,'tDat(i),LOG(CalcK1),LOG(exptK1),LOG(CalcK2),LOG(exptK2),dev1,dev2' 
		if(iOptimize==0)print*,'tDat(),x(1),xDat(),y(1),yDat(),dev1,dev2,itMax'
	endif 
	pMax=5.1d0 !for most systems, anything higher than this just causes more crashes and skews the importance toward scf region
	if(pDatMin > pMax)pMax=pDatMin+0.1 ! for some systems (e.g. nC2+water) the lowest pressure is large. if nConverged=0, LmDif crashes.           
	rmsErr=0
	LIQ=1 ! only liquids.
	if(LOUDER)print*,'KijDevSLE: Kij=',(parm(i),i=1,nParms)
	DO iData=1,nPtsBipDat !	data in module BIPs
		!This needs to be inside the data loop because the temperature may change.
		T=tDat(iData)
		iSolute=IDAT(iData)
		iSolvent=1
		if(iSolute==1)iSolvent=2
		X( iSolute ) = XDAT(iData)
		X( iSolvent)=1-X(iSolute)
		exptXsolute=X( iSolute ) 
		XidSoln( iSolute )=yDat(iData)
		XidSoln( iSolvent)=1-XidSoln(iSolute)
		P=PDAT(iData)
		if(P < zeroTol)then	!guesstimate saturation pressure
!			P=MAX(Pc(1),Pc(2))*2 ! promote liquid root
!			CALL FUGI( T,P,X,NC,LIQ,FUGCL,ZL,ier ) !FUGC=ln(phi)
!			actCo(1:NC)=exp(FUGCL(1:NC)-FUGCPure(1:NC))
			actCo=1 !cancel nonideality.
			Psat1=Pc(1)*10**( 7*(1+acen(1))/3*(1-Tc(1)/T) )
			Psat2=Pc(2)*10**( 7*(1+acen(2))/3*(1-Tc(2)/T) )
			P = ( x(1)*actCo(1)*Psat1+x(2)*actCo(2)*Psat2 )  ! Double to suppress vapor root.
			if(P < 0.1d0)P=0.1d0  
		endif
		if(P <= 0)print*, 'KijDevSLE: PDAT(iData)=P <= 0???'
		do iComp=1,NC
			xOld(1:NC)=zeroTol
			xOld(iComp)=1
			if(LOUDER)print*,'KijDevSLE: Calling fugi for Pure=',iComp
			CALL FUGI( T,P,xOld,NC,LIQ,FUGCL,ZL,ier ) !FUGC=ln(phi)
			fugcPure(iComp)=FUGCL(iComp)
		enddo
		!if(ID(1)==1921 .or. ID(2)==1921)P=pMax ! using a normal guess for P causes BubPL to get stuck in local minimum.
		if(iOptimize>0)then ! 
			if(exptXsolute < zeroTol .or. exptXsolute > 1-zeroTol)then
				ierCalc=3
				print*,'exptX<0? = ',exptXsolute
				pause 'KijDevSLE: exptX < 0?????????????????????????????????????????????????????'
			endif
			initialize=0  ! send the experimental value of x or previous solution if initialize=0.
			if(iData==1.or.iErrCode>0)then
				X( iSolute )=xDat(iData) 
			else
				X( iSolute )=xOld(iSolute) 
			endif
			X(iSolvent)=1-X(iSolute)
			deviate(iData)=0
			ierCalc=0
			CALL FUGI( T,P,X,NC,LIQ,FUGCL,ZL,ier ) !FUGC=ln(phi)
			if(ier(1) > 10)ierCalc=4
			iErrCode=0
			itMax=151
			If(iOptimize==2)then
				Call SLEFL(T,P,NC,Initialize,ITMAX,iSolute,XidSoln,X,iErrCode) !comment this line to use Q&D method.
			elseif(iOptimize==3)then !Wilson model
				ierCalc=0
				X( iSolute ) = XDAT(iData)
				X( iSolvent)=1-X(iSolute)
				!Vmix=X(1)*bVolCc_mol(1)+X(2)*bVolCc_mol(2)
				Om12=bVolCc_mol(2)/bVolCc_mol(1)*exp(-parm(1)/T)
				Om21=bVolCc_mol(1)/bVolCc_mol(2)*exp(-parm(2)/T)
				part= Om12/(X(1)+X(2)*Om12)-Om21/(X(2)+X(1)*Om21)
				actCo(1)=exp( -LOG(X(1)+X(2)*Om12) + X(2)*part )
				actCo(2)=exp( -LOG(X(2)+X(2)*Om21) - X(1)*part )
			else			    	 !NRTL model
				ierCalc=0
				X( iSolute ) = XDAT(iData)
				X( iSolvent)=1-X(iSolute)
				alphaNRTL=0.3
				tau12=(parm(1)/T)
				tau21=(parm(2)/T)
				G12=exp(-alphaNRTL*tau12)
				G21=exp(-alphaNRTL*tau21)
				part1a= tau21*(G21/(X(1)+X(2)*G21))*(G21/(X(1)+X(2)*G21))
				part1b=                                tau12*G12/(X(2)+X(1)*G12)/(X(2)+X(1)*G12)
				part2a= tau12*(G12/(X(2)+X(1)*G12))**2
				part2b=                                tau21*G21/(X(1)+X(2)*G21)**2
				part1=part1a+part1b
				part2=part2a+part2b
				actCo(1)=exp( X(2)*X(2)*part1 )
				actCo(2)=exp( X(1)*X(1)*part2 )
				if(LOUDER)write(*,'(F7.2,f10.4,6f10.3,2f10.3)')T,x(1),parm(1),parm(2),part1a,part2a,part1b,part2b,actCo(1),actCo(2)
			endif
			if(iErrCode>0)ierCalc=4
!			calcXsolute=YDAT(iData)/actCo( iDat(iData) ) !YDAT = ideal solubility.
!			if(LOUDER)print*,'calcX,exptX=',calcXsolute,exptXsolute
!			if(calcXsolute < 1.D-33 )then	! .or. calcXSolute > 1-zeroTol
!				ierCalc=2
!				print*,'calcX,exptX,idX,gam  ',calcXsolute,exptXsolute,YDAT(iData),actCo( iDat(iData) )
!				pause 'KijDevSLE: calcX < 0???'
!			endif
			if(ierCalc==0)then
				if(iOptimize < 3)actCo(1:NC)=exp( FUGCL(1:NC)-fugcPure(1:NC) )
				calcXsolute=YDAT(iData)/actCo(iSolute)	!for quick&dirty Kij optimization, incl Wilson&NRTL
				if(iOptimize==2)calcXsolute=X(iSolute)					!for optimization based on SLEFL
				if(iOptimize >2)itMax=1
				dev1=DLOG(calcXsolute/exptXsolute)
				xOld(iSolute)=calcXsolute ! bootstrap previous temperature.
				!dev1= DLOG( calcXsolute/(1-calcXsolute) )-DLOG( exptXsolute/(1-exptXsolute) ) ! Vladimir's metric, lowers emphasis on very dilute deviations.
				dev1=dev1*100 ! convert to %
				deviate(iData) = dev1 ! use log deviate to minimize bias for bad fits.
				nConverged=nConverged+1 ! Two dev's in above deviate. Consistent with one-sided evaluations when coexistence unavailable.
				!print*,'got one! nConverged=',nConverged
				PAALDK = PAALDK+ABS(dev1) !+ABS(dev2)
				PBIASK = PBIASK+(dev1) !+dev2)
				ssqErr= ssqErr+deviate(iData)*deviate(iData)
				if(LOUDER)print*,'KijDevSLE: Calling fugi.i,xiE,xiC=',iDat(iData),exptXsolute,x(1)
				if(ier(1) > 10)ierCalc=1
				IF(ierCalc.ne.0 .and. LOUDER)WRITE(*,*) 'KijDevSLE: Error in POINT NUMBER = ', iData,ier(1)
				if(ierCalc==0)actCo(1:NC)=exp( fugcL(1:NC)-fugcPure(1:NC) ) !FUGC=ln(phi),
			else
				calcXsolute=0.8686
				dev1=150
				deviate(iData)=150 ! gas+water systems go to nConverged=0 if penalty=200
				!deviate(iData+nPtsBipDat)=150 ! gas+water systems go to nConverged=0 if penalty=200
				if(LAST==1)then
					deviate(iData)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
					!deviate(iData+nPtsBipDat)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
				endif
			endif !last=1 .and. compute activities
			if(LAST==1)WRITE(661,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1,actCo(iSolute),itMax ! write deviations on LastCall
			if(LAST==1.and.LOUDER)WRITE( 6 ,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1 ! write deviations on LastCall
			!WRITE(6 ,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)
		else ! compute LLE @ T if not optimizing.
			itMax=151
			Call SLEFL(T,P,NC,Initialize,ITMAX,iSolute,XidSoln,X,iErrCode)
			calcXSolute=X(iSolute)
			if(iErrCode > 0)cycle
			CALL FUGI( T,P,X,NC,LIQ,FUGCL,ZL,ier ) !FUGC=ln(phi)
			actCo(1:NC)=exp( FUGCL(1:NC)-fugcPure(1:NC) )
			dev1=0
			if(xDat(iData) > 0)then
				dev1=DLOG( x(iSolute)/xDat(iData) )*100
				!if( xDat(iData) > ABS(yDat(iData)) )dev1=DLOG( (1-x(1))/(1-xDat(iData)) )*100
				!dev1=DLOG(  ( 0.5d0-ABS(0.5d0-x(1)) )/( 0.5d0-ABS(0.5d0-xDat(iData)) )  )*100
				nValues=nValues+1
			endif
			deviate(iData) = ( ABS(dev1) ) ! use log deviate to minimize bias for bad fits.
			nConverged=nConverged+1
			PAALDK = PAALDK+deviate(iData)
			PBIASK = PBIASK+(dev1)
			ssqErr= ssqErr+deviate(iData)*deviate(iData)
			if(LAST==1)WRITE(661,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1,actCo(iSolute),itMax ! write deviations on LastCall
			if(LOUDER )WRITE( 6 ,608)T,P,x(1),xDat(iData),yDat(iData),dev1,itMax ! write deviations on LastCall
		endif
	enddo
	close(661)
608	FORMAT( 1X,F8.2,4(1X,F9.6),2(1x,f8.2),i4 )      
617	FORMAT( 1X,i2,F8.2,5(1X,F9.4),(1x,E12.4),i4 )      

	if(nPtsBipDat < 1)pause 'KijVal: nPtsBipDat < 1??? '
	DstdErr=SQRT(DstdErr/(nPtsBipDat))
 	!rmsErr=SQRT( ENORM(nConverged,deviate))
	if(nConverged > 0)then
		rmsErr=SQRT(ssqErr/nConverged)
		PAALDK=PAALDK/nConverged
		PBIASK=PBIASK/nConverged
	else
		rmsErr=8686
		PAALDK=8686
		PBIASK=8686
	endif
	write(*,'(a,i4,2f10.3,7f12.4)')' nConverged,%AALDK,rms%Err,parms=',nConverged,PAALDK,rmsErr,(parm(i),i=1,nparms)
	if(LOUDER)pause 'Check K-values'
	RETURN
	END		 !KijDevSLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Programmed by AV for SpeadKii model
	SUBROUTINE KIIVAL(nc,nParms,parm,deviate,PAADP,pBias,pErrMax,lastCall)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DIMENSION X(NMX) !,Y(NMX),yy(nmx),xx(nmx),ier(20)
	DOUBLE PRECISION parm(nParms),deviate(maxPts)
	character*234 ffile*251
	character*77 errMsg(12)
	LOGICAL LOUDER
	COMMON/eta/etaL,etaV,ZL,ZV
	common/ppVpCoeffs/vpCoeffs(NMX,5)
	DATA iOpt/1/

	errMsg(1)='KiiVal Error: Number of components must equal 1.'
	if(nc.ne.1)then
		if(LOUD)write(*,*)errMsg(1)
		ier=1
		return
	endif
	LOUDER=LOUD
	!LOUDER=.TRUE.
            

	KIJ(1,1)=0
	KTIJ(1,1)=0
	if(lastCall < 0)then !only kijOpt sets last to -1, so that is the only one that affects this.
		write(*,*)'Enter optimization choice: 1=kii, 2=hii'
		read(*,*)iOpt
	endif
	if (nParms==1)then
		KIJ(1,1)=parm(1)
	else
		KIJ(1,1)=parm(1)
		KTIJ(1,1)=parm(2)
	endif
	MAXIT=1111
	PAADP=0
	ssqErr=0
	pBias=0
	DO iData=1,nPtsBipDat
		X(1)=XDAT(iData)
		P=PDAT(iData)
		T=tDat(iData)
		IF(LOUDER)WRITE(*,*)T,P
		IFLAG=0
		ITMAX=MAXIT
		INIT=1    
	    call PsatEar(T,P,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
		!CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)

		IF(ierCode.ne.0  .and. LOUD)WRITE(*,*) 'Error in POINT NUMBER = ', iData

		deviate(iData) = (P-PDAT(iData))/PDAT(iData)*100
		if(ABS(deviate(iData)) > ABS(pErrMax))pErrMax=deviate(iData)
		if(ierCode==0)then
			PAADP = PAADP+ABS(deviate(iData))
			ssqErr= ssqErr+deviate(iData)*deviate(iData)
			pBias=pBias+deviate(iData)
		    !penalize the deviate passed to lmdif, but only compile paadp for converged pts
		else
		    deviate(iData)=10000
        endif
		if(lastCall==1)then
			ffile=TRIM(masterDir)//'\output\tpdat.txt'
			if(DEBUG)ffile='c:\spead\calceos\output\tpdat.txt'
			open(111,file=ffile)
			WRITE(* ,607)tDat(iData),P,deviate(iData),ierCode
			WRITE(111 ,607)tDat(iData),P,deviate(iData),ierCode
		endif
	enddo
	rmsErr=SQRT(ssqErr/nPtsBipDat)
	PAADP=PAADP/nPtsBipDat
	pBias=pBias/nPtsBipDat
	write(*,'(7f10.3)')(parm(i),i=1,nparms),rmsErr,PAADP,pBias
607	FORMAT(1X,F8.4,1X,F10.5,1X,F10.4,1X,i5)      
	RETURN
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE IDACs(NC)
	!C
	!C  PURPOSE:  ITERATIVE LLFLASH EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*1 ANSWER
	!Character*123 outFile
	DoublePrecision xFrac(NMX),gamInf(NMX),fugcLo(NMX),fugcPure(NMX),actCo(NMX),Gmix(20),d2G_dx12(20),dG_dx1(20)
	Integer ierFugi(12)
	!DIMENSION X(NMX),xU(NMX)
	!DIMENSION ZFEED(NMX),kPas(NMX),KINIT(NMX)
	COMMON/eta/etaLo,etaUp,zLo,zUp
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	!data/init/
	if(NC /= 2)then
		pause 'IDACs: Sorry, IDACs only works for NC=2.'
		return
	endif
	TLo=200
	THi=500
	nSteps=30
	PMPa=1
	print*,'Enter TLo(K),THi,nSteps,P(MPa)'
	!read*,TLo,THi,nSteps,pMPa
	deltaT=(THi-TLo)/nSteps
	tKelvin=TLo
	print*,' tKelvin  gamInf(1) gamInf(2) d2Gmin  x1min  xMax1  xMax2  gam2max'
	do while(tKelvin < THi+zeroTol)
		xFrac(1:NC)=zeroTol  ! initialize.
		do iComp=1,NC
			xFrac(iComp)=1
			CALL Fugi(tKelvin,pMPa,xFrac,NC,1,fugcLo,zLo,ierFugi)
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
		!if(LOUDER)print*,'    x1        G/RT        GE/RT        gam1       gam2'
		!if(LOUDER)write(*,'(3F10.3,2F10.3)')0.,0.,0.,gamInf(1),1.0
		do ix=1,19
			xFrac(1)=ix*0.05d0
			xFrac(2)=1-xFrac(1)
			CALL Fugi(tKelvin,pMPa,xFrac,NC,1,fugcLo,zLo,ierFugi)
			actCo(1:NC)=exp( fugcLo(1:NC)-fugcPure(1:NC) )
			Gmix(ix)=xFrac(1)*DLOG(xFrac(1)*actCo(1))+xFrac(2)*DLOG(xFrac(2)*actCo(2))
			GE_RT=xFrac(1)*DLOG(actCo(1))+xFrac(2)*DLOG(actCo(2))
			!if(LOUDER)write(*,'(6F10.3)')xFrac(1),Gmix(ix),GE_RT,actCo(1),actCo(2)
			if(actCo(1) > gam1max)x1Max1=xFrac(1)
			if(actCo(1) > gam1max)gam1max=actCo(1)
			if(actCo(2) > gam2max)x1Max2=xFrac(2)
			if(actCo(2) > gam2max)gam2max=actCo(2)
		enddo
		d2Gmin=1234
		do ix=1,19
			x1=ix*0.05d0
			d2G_dx12(ix)=( Gmix(ix-1)-2*Gmix(ix)+Gmix(ix+1) )/0.0025d0
			dG_dx1=( Gmix(ix-1)-Gmix(ix+1) )/0.1d0
			if(d2G_dx12(ix) < d2Gmin)x1Min=x1
			if(d2G_dx12(ix) < d2Gmin)d2Gmin=d2G_dx12(ix)
			!if(LOUDER)write(*,'(6F10.6)')x1,Gmix(ix),dG_dx1,d2G_dx12(ix)
		enddo

		write(*,'(f9.2,2f10.4,f9.2,5f8.3)')tKelvin,gamInf(1),gamInf(2),d2Gmin,x1min,x1Max1,x1Max2,gam2Max
		tKelvin=tKelvin+deltaT
	enddo
	pause 'IDACs: Success! Read/copy results from screen.'
	return
	end
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LLITER(NC)
	!C
	!C  PURPOSE:  ITERATIVE LLFLASH EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	Character*123 outFile
	DIMENSION X(NMX),xU(NMX)
	DIMENSION ZFEED(NMX),kPas(NMX),KINIT(NMX)
	COMMON/eta/etaLo,etaUp,zLo,zUp
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	!COMMON/KVALUES/kPas
	COMMON/FEED/ZFEED
	!WRITE(*,*)'ENTER DESIRED MAXIMUM FOR NUMBER OF ITERATIONS'
	!READ(*,*)MAXIT
	MAXIT=151
	WRITE(6,*)' '
	outFile=TRIM(masterDir)//'\output\LLIterOut.txt'
	open(601,file=outFile)
	INIT=1
	do J=1,nc
		kPas(J)=KINIT(J)
	enddo        
	!open(61,file='TPXY.txt')
	!write(61,*)' '
	!close(61)
	ANSWER='Y'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
		WRITE(6,*)'ENTER TEMPERATURE(K) AND PRESSURE(MPa)  '
		READ(5,*)T,P
		write(52,*)T,P
		if(NC > 2)then ! For NC=2, K, zFeed, and UOF will be computed from results.
			WRITE(6,*)'ENTER GUESSES FOR xUi/XLiS   '
			READ(5,*)(KINIT(I),I=1,NC)
			WRITE(52,*)'ENTER GUESSES FOR xUi/XLiS   '
			write(52,*)(KINIT(I),I=1,NC)
			WRITE(6,*)' '
			WRITE(6,*)'ENTER FEED COMPOSITIONS      '
			READ(5,*)(ZFEED(I),I=1,NC)
			WRITE(52,*)'ENTER TEMPERATURE(K) AND PRESSURE(MPa)  '
			WRITE(52,*)'ENTER FEED COMPOSITIONS      '
			write(52,*)(ZFEED(I),I=1,NC)
		endif

		!do 6 i=1,1000
		IFLAG=0
		ITMAX=MAXIT
		INIT=2 !check d2G/dx1^2 and search spinodal to determine initial guess.

		!    LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,XU,UOF,iErrCode)
		CALL LLEFL(T,P,NC,INIT,kPas,ITMAX,zFeed,X,xU,UOF,iErrCode)
		if(iErrCode.lt.11.and.iErrCode.gt.0)then
			call PrintErrFlash(iErrCode)
		else 
			if(iErrCode.ge.11)call PrintErrFlash(iErrCode)
			etaLo=etaPasLo
			etaUp=etaPasUp
			zUp=zPasUp
			zLo=zPasLo
			call FlashOut(NC,X,xU,T,P,itmax,Uof) !takes z,eta from common eta
			write(601,'(F8.2,f8.4,4f8.5)')T,P,x(1),xU(1),Uof
		endif

		!  NOTE:  FOR REPEAT LL CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!         PROGRAM BEFORE PROCEEDING.

		INIT=1

		WRITE(6,*)'REPEAT LL CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
	close(601)
182	FORMAT(2X,'COMPONENT IS ',A20,2X,'ID NO. IS ',1X,i5)
183	FORMAT(6E11.4)
659	FORMAT(1X,4(1X,F11.4))
	write(*,*)'Your LL calculations are tabulated in ', TRIM(outFile)
	pause
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LMDevFcn(nPts,nParms,parm,deviate,IFLAG)
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	dimension parm(nParms),deviate(maxPts)
	!C
	!C  SUBROUTINE FCN FOR LMDIFez EXAMPLE.
	nc=2
	last=0
	iflag=last
	iDummy=nPts ! eliminate warning.  deviate must be as maxPts to avoid conflicts later.
	!all we need is the deviate function gotten from kijval below.
	call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,LAST)
	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LMDevLLE(nPts,nParms,parm,deviate,IFLAG)
	USE BIPs, only:maxPts  ! includes nPtsBipDat,TDAT,PDAT,...
	implicit doublePrecision(A-H,K,O-Z)
	dimension parm(*),deviate(maxPts)
	!C
	!C  SUBROUTINE FCN FOR LMDIFez EXAMPLE.
	nc=2
	last=0
	iflag=last  ! disable last indicator when called by LMDif
	iDummy=nPts ! eliminate warning.  deviate must be as maxPts to avoid conflicts later.
	!all we need is the deviate function gotten from kijval below.
	call KijDevLLE(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LAST)
	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LMDevSLE(nPts,nParms,parm,deviate,IFLAG)
	USE BIPs, only:maxPts  ! includes nPtsBipDat,TDAT,PDAT,...
	implicit doublePrecision(A-H,K,O-Z)
	dimension parm(*),deviate(maxPts)
	!C
	!C  SUBROUTINE FCN FOR LMDIFez EXAMPLE.
	nc=2
	last=0
	iflag=last  ! disable last indicator when called by LMDif
	iDummy=nPts ! eliminate warning.  deviate must be as maxPts to avoid conflicts later.
	!all we need is the deviate function gotten from kijval below.
	call KijDevSLE(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LAST)
	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE RegSpeadDevFcn(nPts,nParms,parm,deviate,IFLAG)
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	dimension parm(nParms),deviate(maxPts)
	!C
	!C  SUBROUTINE FCN FOR LMDIFez EXAMPLE.
	iDummy=nPts ! eliminate warning.  deviate must be as maxPts to avoid conflicts later.
	nc=1
	last=0
	iflag=last
	!all we need is the deviate function gotten from kiiVal below.
	call KIIVAL(nc,nParms,parm,deviate,PAADP,pBias,pErrMax,LAST)
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OptiBIPs(NC) !OB option
	!C
	!C  PURPOSE:  BP CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE BINARY INTERACTION COEFFICIENT (KIJ) FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
	!C  PROGRAMMED BY:  JRE 2/96
	!C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
	!C  REFERENCE:NUMERICAL RECIPES
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	INTEGER ier(20)
	DOUBLE PRECISION parm(5),stdErr(5),deviate(maxPts),X(NMX),Y(NMX)
	EXTERNAL LMDevFcn
	character inFile*251,outFile*251

    if(LOUD)then
	    IF(NC.NE.2)pause 'THIS OPTION IS FOR 2 COMPONENTS ONLY'
    end if
	IF(NC.NE.2)RETURN

	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (..\CalcEos\VleData\FN.FT)'
	WRITE(*,*)'FIRST LINE =NPTS, 2ND LINE = TK, PMPA, X,Y, ...'
	READ(*,'(A)')FN
	inFile=TRIM(masterDir)//'\VleData\'//TRIM(FN)
	if(DEBUG)inFile='c:\spead\calceos\VleData\'//TRIM(FN)
	open(51,file=inFile)
	!outFile=TRIM(masterDir)//'\output\KijOut.txt'
	!if(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	!open(67,file=outFile)
	READ(51,*,ERR=861)NPTS
	DO I=1,NPTS
		READ(51,*,ERR=862)tDat(I),PDAT(I),XDAT(I)
		WRITE(*,*)tDat(I),PDAT(I),XDAT(I)
	enddo
	close(51)
	MAXIT=1111
	INIT=1
	WRITE(*,*)'THIS ROUTINE ASSUMES A BINARY SYSTEM FROM MAIN PROGRAM'
	if(iEosOpt.eq.4)then
		write(*,*)'Note: SpeadMD is site-based and must take its HB bips from bipAD.txt and bipDA.txt'
		write(*,*)'Options:  1-kij, 2-kij,Hij'
	elseif(iEosOpt.eq.5.or.iEosOpt.eq.8.or.iEosOpt.eq.9)then	 !AFG 2011 , Added EOS Opt. 9
		write(*,*)'Note: SpeadMD is site-based and must take its HB bips from bipDA.txt (bipAD is inferred from bipDA)'
	elseif(iEosOpt.eq.3)then
		write(*,*)'Options:  1-kij, 2-tauIJ,tauJI, 3-kij,tauIJ,tauJI, '
		write(*,*)'          4-kij,kTij,tauIJ,tauJI, 5-kij,kTij,tauIJ,tauJI,alphaIJ'
	elseif(iEosOpt.eq.7)then
		write(*,*)'Options:	 2-tauIJ,tauJI,3-alphaIJ,tauIJ,tauJI'
	else
		write(*,*)'Options:  1-kij, 2-kij,kTij '
	endif

	READ(5,*)iRegOpt
	nParms=iRegOpt
	WRITE(6,*)'ENTER INITIAL GUESSES'
	READ(5,*)(parm(I),I=1,nParms)                                          
	!C	LWA>5*nParms+NPTS ~ 275
	LWA=275
	!C
	!C  SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
	!C  UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED, 
	!C  THIS IS THE RECOMMENDED SETTING.
	!C
	!c      TOL=SQRT(DPMPAR(1))
	TOL=0.0001
	!C
	!c  factor controls the magnitude of the first step from the initial guess.
	factor=0.01	 !small value starts slow
	if(LOUD)write(*,*)'KIJ,paadp,rmserr'
	CALL LMDifEZ(LMDevFcn,nPts,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
	rmsErr=SQRT( ENORM(nPts,deviate) )
	if(LOUD)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
	if(iEosOpt.le.3.or.iEosOpt.eq.7)then
		if(LOUD)write(* ,*)'   KIJ','      ',' KTIJ','       ',' xsBIP12',' xsBIP21'
	elseif(iEosOpt.le.5.or.iEosOpt.eq.7.or.iEosOpt.eq.8.or.iEosOpt.eq.9)then	 !AFG 2011 , Added EOS Opt. 9
		if(LOUD)write(* ,*)'   KIJ','      ',' HIJ','       ','  kTij'
	else
		if(LOUD)write(* ,*)'   KIJ','      ',' KTIJ' !,'      ',' xsBIP12',' xsBIP21'
	endif
601	format(a,6x,a,7x,a,9x,a)
606 FORMAT(1X,F8.4,5(1X,F10.4))
	if(LOUD)write(* ,606)(parm(i),i=1,nParms)
	if(LOUD)write(* ,*)'StdErr'
	if(LOUD)write(* ,606)(stdErr(i),i=1,nParms)
	if(LOUD)pause 'Converged'
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	open(68,file=outFile,ACCESS='APPEND')
	!OPEN(67,FILE='KIJOUT.txt')
	write(68,601)
	WRITE(68,606)(parm(i),i=1,nParms)
	write(68,*)'StdErr'
	WRITE(68,606)(stdErr(i),i=1,nParms)
	if(LOUD)write(* ,600)
	write(68,600)
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )

	LAST=1 !cause KijVal to print results
	call KIJVAL(nc,nParms,parm,deviate,PAAD1,DstdErr,Last)
607 FORMAT(1X,F8.4,1X,F10.4,1X,F10.8,1X,F10.8,1X,F10.4,1X,i5)      
	DO iPt=1,NPTS
		X(1)=XDAT(iPt)
		X(2)=1-X(1)
		P=PDAT(iPt)
		T=tDat(iPt)
		IFLAG=0
		ITMAX=MAXIT
		INIT=1    
		CALL BUBTL(T,X,NC,INIT,P,ITMAX,Y,ier)
		IF(PDAT(iPt).GT.2.0)IER(1)=1
		if(ier(1).ne.0)then
			CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
		endif
		devi=(t-tDat(iPt))/tDat(iPt)*100
		WRITE(68,607)t,P,X(1),Y(1),devi,ier(1)
		WRITE(6 ,607)t,P,X(1),Y(1),devi,ier(1)
	enddo
	if(LOUD)write(* ,601)
	if(LOUD)write(* ,606)(parm(i),i=1,nParms)
	if(LOUD)write(*,*)'%Err=',PAAD1,' rmsErr = ', rmserr
	if(LOUD)write(*,*)'errCode=',iErrCode
	write(68,*)'%Err=',PAAD1,' rmsErr = ', rmserr
	write(68,*)'errCode=',iErrCode
	CLOSE(68)
	if(LOUD)pause 'These results are tabulated in KijOut.txt.'
	RETURN
861	continue
	if(LOUD)write(*,*)'nPts=',nPts,'?'
	if(LOUD)pause 'Error reading nPts in data file. Check the path and file.'
	close(68)
	close(51)
	return
862	continue
	if(LOUD)pause 'Error reading data line in data file. Check the path and file.'
	close(68)
	close(51)
	return
	END	 ! OptiBIPs

!***********************************************************************
      
	SUBROUTINE PhasEnv(NC)
	!C
	!C  PURPOSE:  COMPUTE phase envelope GIVEN x,y.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character outFile*251
	!DIMENSION 
	DIMENSION X(NMX),Y(NMX),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT

	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	open(61,file=outFile)
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x, 'etaL',6x, 'etaV')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
1000  CONTINUE
	WRITE(6,*)'ENTER COMPOSITION   '
	READ(5,*)(X(I),I=1,NC)
	if(LOUD)write(*,600)
	write(61,600)
	MAXIT=1111
	zero=0

	ITMAX=MAXIT
	tMax=1111
	CALL BootBP(tMax,x,NC,INIT,pMax,ITMAX,y,ier,0)
!c	CALL ErrCheck(ier,iflag)
	if(LOUD)write(* ,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	if(ier(1).eq.0)init=1

!C      PAUSE

	if(LOUD)write(*,*)'Starting dew point curve.'
	write(61,*)'Dew points.'
	do i=1,nc
		Y(i)=X(i)
	enddo

      ITMAX=MAXIT
	pMax=111
	CALL BootDT(tMax,X,NC,INIT,pMax,ITMAX,Y,ier,0)
!c	CALL ErrCheck(ier,iflag)
	if(LOUD)write(* ,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	if(ier(1).eq.0)init=1

86	continue
	if(LOUD)pause 'These results are tabulated in TPXY.txt.'
	CLOSE(61)   
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE PropTable(NC,ier)
	!C
	!C  PURPOSE:  Compile an isobaric property table in increments of T. 
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	!			   uRes_RT, sRes_R, aRes_RT, hRes_RT, CpRes_R, cvRes_R, cmprsblty !cmprsblty=(dP/dRho)T*(1/RT)

	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	LOGICAL LOUDER
	!CHARACTER*1 ANSWER
	CHARACTER*77 OutFile
	DIMENSION fugcL(NMX) !fugcV(NMX),
	DIMENSION X(NMX),Y(NMX),wFrac(NMX) ,gmol(NMX) ! !,VLK(NMX)
	dimension ier(12)
	COMMON/eta/etaL,etaV,zLiqDum,zVapDum
	common/eta2/eta
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	IF(NC > 1)THEN
		if(LOUDER)PAUSE 'SORRY, THIS ROUTINE IS ONLY FOR PURE FLUIDS AT THIS TIME.'
		RETURN
	ENDIF
	LOUDER=LOUD
	LOUDER=.TRUE.
	X(1) = 1
	Y(1) = 1
	rMwAvg=rMw(1)
	wFrac(1) = 1
	OutFile=TRIM(MasterDir)//'\output\PropTable.txt'
	open(601,file=OutFile)
	WRITE(6,*)'ENTER PRESSURE(MPa)  '
	READ(5,*)P !,T
	WRITE(601,*)'PRESSURE(MPa)   '
	write(601,*)P
	print*,' Enter Tstart, Tstop, and Tincrement'
	read*,tStart,tMax,tInc
	tKelvin=tStart-tInc
	print*,     '   T(K)  rho(g/cc)   Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R   '  
	write(601,'(A)')'   T(K)   rho(g/cc)  Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R   SRes_R'    
1000  CONTINUE

		tKelvin=tKelvin+tInc
		!pSatEst = PC(1)*10**(  7*(1+ACEN(1))/3*( 1-TC(1)/tKelvin )  )
		tSatEst=Tc(1)/(   1-(  LOG10(P/Pc(1))*3/( 1+acen(1) )/7  )   )
		LIQ=1
		if(tKelvin > tSatEst)LIQ=0 

		IFLAG=0
		CALL Fugi(tKelvin,P,X,NC,LIQ,fugcL,zFactor,ier)
		rhoG_cc=rMwAvg*P/(zFactor*RGAS*tKelvin)
		DO I=1,6
			IF (ier(I).NE.0)IFLAG=1
		enddo
		IF (IFLAG.EQ.1)WRITE(6,'(a,22I4)')'ier',(ier(I),I=1,12)
		IF (IFLAG.EQ.1)WRITE(6,*)'error in fugi calculation'

		!for Fugi:
		!C  ier = 1 - AT LEAST ONE ERROR
		!C        2 - NOT USED
		!C        3 - NOT USED
		!C        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
		!C        5 - rho IS -VE
		!C        6 - TOO MANY Z ITERATIONS
		!C        7 - eta > 0.53
		!C		11 - goldenZ instead of real Z.

		!if(LOUDER)write(*,*)'DEPARTURE FUNCTIONS: (HLo-Hig)/NkT (SLo-Sig)/Nk (HUp-Hig)/NkT (SUp-Sig)/Nk'
		if(LOUDER)write(* ,610)tKelvin,rhoG_cc,zFactor,aRes_RT,uRes_RT,cmprsblty,CvRes_R,CpRes_R !,DSONK
		WRITE(601,610)tKelvin,rhoG_cc,zFactor,aRes_RT,uRes_RT,cmprsblty,CvRes_R,CpRes_R, sRes_R
		!Check derivatives using FuVtot
		vTotCc=rMwAvg/rhoG_cc
		gmol(1)=1
		isZiter=0 !=>compute derivative props
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGCL,zFactor,iErr)
		pNew = zFactor*rGas*tKelvin/vTotCc
		if(LOUDER)write(* ,612)tKelvin,pNew,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R
		!write(61,612)tKelvin,pNew,zFactor,aRes_RT,uRes_RT,cmprsblty,cvRes_R,CpRes_R,hRes_RT
612	FORMAT(2F9.3,1x,1(F7.4),8(F9.4,1X))      
	if(tKelvin < tMax)GOTO 1000
610	FORMAT(F9.2,1x,2(F8.5),8(F9.4,1X))      
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(2X,'LIQUID X',2X,'Upper Y',1X,'fugacityLo(MPa)',1X,'fugacityUp(MPa)    lnPhiLo',3X,'lnPhiUp')
605 FORMAT(3X,'K-RATIO',7X,'T(K)',3X,'P(MPa)',2X,'etaLo',4X,'etaUp',6X,'ZLo  ',6X,'ZUp')
604	FORMAT(2(3x,F7.4),2(3x,e11.4),2(3x,f12.7))
606	FORMAT(e12.5,3X,F7.2,1X,F7.4,1X,2(F6.4,3X),2(F9.4,1X))
	if(LOUDER)WRITE(*,*)'Output is in:',TRIM(OutFile)
	if(LOUDER)pause 'Success! cf. PropTable.txt'
	close(601)
	RETURN		
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!Proutine
	SUBROUTINE PsatChemPo(tK,pMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCode)
	!COMPUTE VAPOR PRESSURE GIVEN tK AND INITIAL GUESS (IE.pMPa)
	!Ref: Eubank et al, IECR 31:942 (1992)
	!AU:  JRE Aug2010
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!character tabChar*1
	DIMENSION fugc(NMX),xFrac(NMX) 
	DIMENSION IER(12)
	CHARACTER*77 errMsg(0:11) !,errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	errMsg(1)='No spinodal max/min'
	errMsg(2)='Psat iteration did not converge'
	errMsg(3)='Liquid FUGI call failed on last iteration'
	errMsg(4)='Vapor  FUGI call failed on last iteration'
	errMsg(5)='zVap=zLiq on last iteration'
    errMsg(6)='PsatChemPo: FUGI returned T < Tmin error'
	errMsg(7)='PsatChemPo: zLiq < 0 on initial call to Fugi(P~0)'
	errMsg(8)='PsatChemPo: zLiq < 0 on iterative call to Fugi'
	!if(LOUD)write(*,*)'EAR method'
	tMin=TC(1)*0.25
	xFrac(1)=1
	NC=1
	ierCode=0 
	itMax=66
	zCritEff=ZC(1)
	if(iEosOpt==1 .or. iEosOpt==11)zCritEff=0.3074
	if(iEosOpt==2 .or. iEosOpt==4)zCritEff=0.34  !ESD
	rhoCrit=PC(1)/(zCritEff*rGas*TC(1))
	rhoVap=rhoCrit*1.000 !-ve value on input means use default value for etaHi
    if(LOUD)then 
        if(DEBUG) then
	        OPEN(67,file='c:\spead\calceos\output\isotherm.txt') !bonus: write the isotherm from iterations on spinodal
        else
	        OPEN(67,file='spinodal.txt')
        endif
	    write(67,*)' T(K)    rho(mol/cc)   P(MPA)        Z     DA_NKT    DU_NKT    eta'
    endif
	!print*,'Psat: calling spinodal for vapor.'
	call Spinodal(NC,xFrac,tK,0,rhoVap,zVap,aDepVap,dU_NkT,iErrVap)
    if(iErrVap==3)then !first definite call to fuVtot suffices to check for T < Tmin.
        ierCode=6
        goto 86
    endif
	rhoLiq=rhoVap !-ve value on input means use default value for etaLo
	!print*,'Psat: calling spinodal for liquid.'
	call Spinodal(NC,xFrac,tK,1,rhoLiq,zLiq,aDepLiq,dU_NkT,iErrLiq)
	if(LOUD)close(67)
	pMax=rGas*tK*(rhoVap*zVap) 
	pMin=rGas*tK*(rhoLiq*zLiq)
    relDiffP=ABS(pMax-pMin)/pMax
	if(iErrVap.ne.0 .or. iErrLiq.ne.0 .or. relDiffP < 1D-3)then
		ierCode=1
		if(LOUD)write(*,'(a,f8.1,a)')' PsatChemPo: Spinodal error at T = ',tK,' T > Tc?'
		!if(LOUD)pause 
		goto 86
	endif
    !if(pMin < 0)pMin=0
	etaLiq=rhoLiq*bVolCC_mol(1)
	etaVap=rhoVap*bVolCC_mol(1)
	if(LOUD)print*,'Psat after spinodal: pMin,pMax= ',pMin,pMax
	if(LOUD)print*,'Psat after spinodal: etaLiq,etaVap= ',etaLiq,etaVap
	pOld=(pMin+pMax)/2 !when approaching critical, this works well.
	if( (rhoVap/rhoLiq) < 0 .and. LOUD)pause 'Psat error after spinodals: rhoVap/rhoLiq < 0'
	pEuba=rGas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995), assuming 1.2x shifts on V&L relative to spinodals.
	!pOld=pEubank/10  !JRE: I find Eubank overestimates P sometimes. It even get higher than pMax
	if(pMin < 0)then ! if rhoLiq exists at p~0, then Razavi works.
		if(pOld < 0)pOld=1D-5 !If p < 0, then the vapor density goes negative and log(rhoLiq/rhoVap) is indeterminate.
		!NOTE: Don't use pOld=0 when pOld > 0. You might get Golden Error because vdw loop doesn't cross zero!
		print*,'calling fugi for liquid to construct Razavi initial guess.'
		call FUGI(tK,pOld,xFrac,NC,1,FUGC,zLiq,ier) !LIQ=3=>liquid root with no fugc calculation
		if(ier(1)==0)then
			AresLiq=Ares_RT
			rhoLiq=pOld/(zLiq*rGas*tK)			 !AresVap  +  Zvap-1 -ln(Zvap) =AresLiq+ZLiq-1-ln(ZLiq) 
			!rhoLiq=etaL/bVolCC_mol(1)		 =>  !B2*rhoVap+B2*rhoVap+ln(rhoVap/rhoLiq) =AresLiq+0-1
			rhoVap=rhoLiq*EXP(AresLiq-1) !Razavi, FPE, 501:112236(19), Eq 19. Ares_RT from GlobConst
			pSatRazavi=rhoVap*rGas*tK
			!zLiq=pTest/(rhoLiq*rGas*tK) !this takes some of the imprecision off zLiq
			!pEuba2=rGas*tK/(1/rhoVap-1/rhoLiq)*( Ares_RT- 0 +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
			!NOTE: When 1/rhoVap >> 1/rhoLiq, pEuba2=rGas*tK*rhoVap*( Ares_RT-ln(rhoVap/rhoLiq) )=pTest*( Ares_RT - (Ares_RT-1) )
			call FUGI(tK,pSatRazavi,xFrac,NC,0,FUGC,zVap,ier)
			B2cc_mol=(zVap-1)*(zVap*rGas*tK)/pSatRazavi
			rhoVap=rhoVap*exp(-2*B2cc_mol*rhoVap)! => rhoVap=rhoLiq*EXP(Ares_RT-1-2*B2*rhoVap), but Ares-RT changes with FUGI call, so reuse like this.
			pSatRazavi=rhoVap*rGas*tK*(1+B2cc_mol*rhoVap)
			zLiq=pSatRazavi/(rhoLiq*rGas*tK)
			FugcLiq=AresLiq+zLiq-1-LOG(zLiq)
			if(LOUD)print*,'PsatEar: FugcVap,FugcLiq=',Fugc(1),FugcLiq
			if(pSatRazavi < pMax)then
				pOld=pSatRazavi
			else
				if(LOUD)print*,'PsatChemPo: pSatRazavi > pMax?'
				pOld=zeroTol*10
				if(pMax > zeroTol)pOld=pMax/2
			endif
		else
			print*,'PsatChemPo: Error from call to fugi for Razavi. ier=',ier 
		endif ! ier(1)==0
	endif  ! pMax < 0

	fOld = 0.1  ! above gives good guess for P, but assumes chempoV=chempoL.  So set fOld= 0.1 and let iterations sort it out.
	etaLiq=rhoLiq*bVolCC_mol(1)
	!if(LOUD .and. ABS(etaL-etaLiq)/etaLiq > 0.0001 .and. LOUD)pause 'Psat: warning lost precision in rhoLiq'
	if(LOUD)print*,'Initial guess for pSat,etaLiq,fErr=',pOld,etaLiq,fOld

	pMPa=pOld/1.1
	fBest=1234
	pBest=1234
		
	do iter=1,itMax !iterate on pOld according to Eubank criterion
		!print*,'calling fugi for liquid iteration.'
		call FUGI(tK,pMPa,xFrac,NC,3,FUGC,zLiq,ier) !LIQ=3=>liquid root with no fugc calculation
		chemPotL=fugc(1)
		gResLiq=aRes_RT+zLiq-1
		if(ier(1).ne.0)ierCode=3 !declare error but don't stop. if future iterations give valid fugi(), then no worries
		if(zLiq > 0)then
			rhoLiq=pMPa/(zLiq*rGas*tK)
		else
			ierCode=8
			goto 86
		endif
		!print*,'calling fugi for vapor iteration.'
		call Fugi(tK,pMPa,xFrac,NC,2,FUGC,zVap,ier) !LIQ=2=>vapor root with no fugc calculation
		chemPotV=fugc(1)
		gResVap=aRes_RT+zVap-1
		!print*,'computing rhoVap. zVap,rGas,tK=',zVap,rGas,tK
		if(zVap > 0)rhoVap=pMPa/(zVap*rGas*tK)
		if(ier(1).ne.0)ierCode=4 !declare error but don't stop. if future iterations give valid fugi(), then no worries.
		if( (zVap/zLiq) < 0 .and. LOUD)pause 'Psat error : rhoVap/rhoLiq < 0'
		if(ABS(zVap-zLiq) < 1e-3)then
			if(LOUD)write(*,*)'rho,zLiq,zVap',rhoLiq,zVap,zLiq
			if(LOUD)pause 'PsatChemPo: zVap=zLiq'
			ierCode=5
			goto 86
		endif
		!rLogRhoRat=DLOG(rhoVap/rhoLiq)	!for debugging
		!fErr=chemPotL-chemPotV
		!fErr=gDepLiq-gDepVap
		fErr=( gResLiq+zLiq-1-DLOG(zLiq/zVap) )-( gResVap+zVap-1 )
		if( ABS(fErr) < fBest )then
			fBest=ABS(fErr)
			pBest=pMPa
		endif
		if(ABS(fErr-fOld) < 1.e-22)then
			!if(LOUD)write(*,*) 'Error in Psat, fErr ~ fOld'
			ierCode=14
			goto 86
		endif
		!denom=(fErr-fOld)
		!print*,'computing change fOld ,fErr,denom= ',fOld,fErr,denom
		change=fErr*(pMPa-pOld)/( fErr-fOld )
		if(LOUD)write(*,'(a,5e12.4)')' P,%chng,fuL,fuV,fErr',pMPa,change/pMPa*100,chempotL,chemPotV,fErr
		!print*,'pMPa,%change',pMPa,change/pMPa*100
		fOld=fErr
		pOld=pMPa
		IF(ABS(change/pMPa).GT.0.2)then
			change=DSIGN(0.2d0,change)*pMPa !step limiting 
		endif
		!if(LOUD)write(*,'(i3,5f11.5)')iter, PMpa, chemPotL,  chemPotV, change/pMPa*100
		pMPa=pMPa-change
		if(pMPa.lt.1.d-11)pMPa=2.d-11 !don't let the pressure go too low.
		tol=1.D-6
		if(pMPa < 1.D-7)tol=1.D-3
		IF( ABS(change/pMPa) < tol ) exit
	enddo  !iteration on pMPa
	if(iter.ge.itMax)then
		pMPa=pBest
		if(LOUD)print*,'Psat: iterations exceeded. P,fErr=',pMPa,fBest
		ierCode=9
	endif
	if(ierCode.ne.0 .and. pSatRazavi < 0.01)then
		pMPa=pSatRazavi
		ierCode=0 ! At P < 0.01 MPa, pSatRazavi must be pretty good. 
	endif
	!check that FugiTpt did not give warning on last calls
	!pMPa=rGas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
	!pMPa=pOld*( (aDepVap-aDepLiq)/(zVap-zLiq) +DLOG(zVap/zLiq)/(zVap-zLiq) )
	uSatV=DUONKT  ! Last call to FUGI was for vapor.
	if(zVap > 0)rhoVap=pMPa*rMw(1)/(zVap*rGas*tK)
86	continue
	!print*, 'Psat=',pMPa,' Calling Fugi for fugc.'
	if(ierCode==0)call Fugi(tK,pMPa,xFrac,NC,1,FUGC,zLiq,ier) !one last call to fugi for chemPo and debugging.
	if(ier(1)==0)chemPot=FUGC(1) !since chemPotL=chemPotV at pSat       
	uSatL=DUONKT
	if(zLiq > 0)rhoLiq=pMPa*rMw(1)/(zLiq*rGas*tK)
    if(LOUD)write(*,'(f10.2,3f14.10,i2)')tK,pMPa,rhoVap,rhoLiq,ierCode
	if(LOUD)pause 'check Psat'
	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!Proutine
	SUBROUTINE PsatChemPoOld(tK,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode)
	!COMPUTE VAPOR PRESSURE GIVEN tK AND INITIAL GUESS (IE.pMPa)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	!parameter(nTypesMax=9,nWellsMax=11,nPolyOrder=4)
	character tabChar*1
	!common/TptWells/eokWells(nTypesMax,nWellsMax),tptLamda(nWellsMax),sigmaTpt(nTypesMax),idSite(nTypesMax),nDegree(nTypesMax),nTypes,nTptWells
	!common/HbParms/bondVolNm3(nTypesMax),eHbKcal_mol(nTypesMax),idType(nTypesMax),ndHb(nTypesMax),nDs(nTypesMax),nAs(nTypesMax),nTypesHb
	!common/TptParms/a1Coeff(nPolyOrder),a2Coeff(nPolyOrder),zRefCoeff(3),bVolEff(0:3),rMw,idComp !for isochore
	!CHARACTER*20 NAME
	DIMENSION gmol(NMX),fugc(NMX),xFrac(NMX) !NAME(NMX),etaList(maxEtaList),PMpa(maxEtaList) !,ierFugi(20)
	DIMENSION IER(12)
	!CHARACTER*77 errMsg(0:11),errMsgPas
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	!COMMON/NAME/NAME
	!COMMON/fugNum/FUG_num,P_num,dFUG_dN,FUG,dP_dN
	!COMMON/EsdParms/eokP,KCSTAR,DH,C,Q,VX,ND,NDS,NAS
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	data isochore/0/ !for debugging set isochore=1, otherwise isochore=0. This lets you iterate on temperature until pLiq~0
	!ierCode=13 no convergence in 333 iterations
	tabChar=char(9) !ascii tab character
	!vEffNm3=bVolEff(0)/602.22
	tMin=111
	xFrac(1)=1
	NC=1
	ierCode=0 !for debugging set isochore=1, otherwise isochore=0
	if(isochore)then !FYI: JRE copied isochore section from RegSteps Aug2,2010. NOT debugged or tested. 
		!write(6 ,601)idComp,tabChar,idComp,tabChar,idComp,tabChar,(zRefCoeff(iCoeff),tabChar,iCoeff=1,3),(a1Coeff(iCoeff),tabChar,iCoeff=1,nPolyOrder),(a2Coeff(iCoeff),tabChar,iCoeff=1,nPolyOrder) ,&
		!	vEffNm3,tabChar,tMin,tabChar,rMw,tabChar,nTypes,tabChar,(nDegree(iType),tabChar,idType(iType),tabChar,iType=1,nTypes)
!601		format(i5,a1,I5,a1,i11,3(a1,f8.4),4(a1,f11.3),4(a1,f12.1),a1,f10.7,a1,F6.0,a1,F9.3,<1+2*nTypes>(a1,i4),a1,\) !the \ stops the line feed so the compname can be written on the same line.
!		if(LOUD)write(*,*)' ' !end the previous line
		isochore=0 !terminate isochore iteration after the first call
		isZiter=1 !derivatives of P,fug wrt N not required
		write(*,*)'Enter liquid density (g/cc, <0 to stop)'
		read(*,*)rhoLiq
		if(rhoLiq.lt.0)goto 20
		vTotCc=1
		gmol(1)=rhoLiq/rMw(1)
10		continue !correct unit for rho is g/cc
		!eta=rhoLiq/rMw*bVolEff(0)
		!zRef=(  1+eta*( zRefCoeff(1)+eta*(zRefCoeff(2)+eta*zRefCoeff(3)) )  )/(1-eta)**3
		tHi=1e11 !get zRef assuming tHi ~ infinity.
		call FuVtot(isZiter,tHi,vTotCc,gmol,NC,FUGC,zRef,iErr)
		zLiq=zRef
		beps=0
		dBeps=0.05 
		if(LOUD)write(*,*)'1000/T   zCalc'
		do while(zLiq.gt.0)
			beps=beps+dBeps
			!tK=eokWells(1,1)/beps !only works for single site molecule (e.g. xenon)
			!zLiq=ZCalcTpt(rhoLiq,tK,aDepLiq,uDepLiq,ierZ)
			call FuVtot(isZiter,tK,vTotCc,gmol,NC,FUGC,zLiq,iErr)
			tInv=1000/tK
			if(LOUD)write(*,*)tInv,zLiq !,eta
		enddo
		tInvSat=tInv+(zLiq-0)*(tInv-0)/(zLiq-zRef)
		tSatK=1000/tInvSat
		write(*,'(a,f8.1,a)')'Interpolated value of tSatK=',tSatK,'Enter your best guess of tSatK.'
		read(*,*)tSatK
		rhoVap=pMPa*rMw(1)/(rGas*tSatK) !assume pMPa guess is accurate and ig law
		!zVap=ZCalcTpt(rhoVap,tK,aDepVap,uDepVap,ierZ)
		gmol(1)=rhoVap/rMw(1)
		call FuVtot(isZiter,tK,vTotCc,gmol,NC,FUGC,zVap,iErr)
		!aDepVap+zVap-1-ln(zVap) = aDepLiq+zliq-1-ln(zLiq)
		!zLiq = exp(aDepLiq-1-aDepVap)
		pGuess=zLiq*rGas*tK*rhoLiq/rMw(1)
		write(*,*)'rho,tSatK,pSatMPa'
		write(*,*)rhoLiq,tSatK,pGuess
		write(*,*)'Enter another liquid density (g/cc, <0 to stop)'
		read(*,*)rhoLiq
		if(rhoLiq.gt.0)goto 10
	endif !if(isochore)
20	continue
	!pOld=pGuess
	!pOld=6.849716594694031E-003
	!CALL Fugi(tK,pOld,'liq',rhoL,zL,chemPotL,uSatL,ierCodeL)
	!IF(ierCodeL.gt.10)ierCode=11
	!if(ierCode.ne.0)then
		!write(*,*)'Error in Psat.  Fugi failed on first call.' ! ierL,ierV = ' ,ierCodeL,ierCodeV
    if(LOUD)then 
        if(DEBUG) then
	        OPEN(67,file='c:\spead\calceos\output\isotherm.txt') !bonus: write the isotherm from iterations on spinodal
        else
	        OPEN(67,file='spinodal.txt')
        endif
	    write(67,*)' T(K)    rho(mol/cc)   P(MPA)        Z     DA_NKT    DU_NKT    eta'
    endif
	rhoVap=-1 !-ve value on input means use default value for etaHi
	call Spinodal(NC,xFrac,tK,0,rhoVap,zVap,aDepVap,dU_NkT,iErr)
	if(iErr.ne.0)then
		if(LOUD)pause 'PsatChemPo: Spinodal indicates T > Tc.'
		iercode=1
		goto 86
	endif 
	rhoLiq=-1 !-ve value on input means use default value for etaLo
	call Spinodal(NC,xFrac,tK,1,rhoLiq,zLiq,aDepLiq,dU_NkT,iErr) 
	if(LOUD)close(67)
	!endif
	pOld=rGas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
	pOld=zVap*rhoVap*rGas*tK*0.5 ! make it a little less than pMax.
	pOld=0.00136 !scvp at 333K for nbutane
	call Fugi(tK,pOld,xFrac,NC,1,FUGC,zLiq,ier)
	chemPotL=fugc(1)
	aDepLiq=DAONKT
	if(ier(1).ne.0)ierCode=11
	call Fugi(tK,pOld,xFrac,NC,0,FUGC,zVap,ier)
	chemPotV=fugc(1)
	aDepVap=DAONKT
	if(ier(1).ne.0)ierCode=12
	pMPa=pOld*( (aDepVap-aDepLiq)/(zVap-zLiq) )!+DLOG(zVap/zLiq)/(zVap-zLiq) )
	pMPa=pOld*0.5
	iter=0                                  
	fOld=chemPotL-chemPotV
	!if(LOUD)write(*,*)' iter PMpa, chemPotL,  chemPotV, %change '
1002	CONTINUE
		iter=iter+1
		if(iter.gt.333)then
			ierCode=13
			!if(LOUD)write(*,*)'Error in Psat-no convrg'
			goto 86
		endif
		call Fugi(tK,pMPa,xFrac,NC,1,FUGC,zLiq,ier)
		chemPotL=fugc(1)
		aDepLiq=DAONKT
		if(ier(1).ne.0)ierCode=11
		call Fugi(tK,pMPa,xFrac,NC,0,FUGC,zVap,ier)
		chemPotV=fugc(1)
		aDepVap=DAONKT
		if(ier(1).ne.0)ierCode=12
		if(ierCode.ne.0)then
			if(LOUD)write(*,*)'Error in Psat.  Fugi failed on first call.' ! ierL,ierV =',ierCodeL,ierCodeV
		endif
		ierCode=0
		fErr=chemPotL-chemPotV
		if(ABS(fErr-fOld).lt.1.e-33)then
			!if(LOUD)write(*,*) 'Error in Psat, fErr ~ fOld'
			ierCode=14
			goto 86
		endif
		change=fErr*(pMPa-pOld)/(fErr-fOld)
		fOld=fErr
		pOld=pMPa
		IF(ABS(change/pMPa).GT.0.3)then
			change=DSIGN(0.3d0,change)*pMPa !step limiting 
		endif
		!if(LOUD)write(*,'(i3,5f11.5)')iter, PMpa, chemPotL,  chemPotV, change/pMPa*100
		pMPa=pMPa-change
		if(pMPa.lt.1.d-11)pMPa=2.d-11 !don't let the pressure go too low.
	IF(ABS(change/pMPa).GT.1.D-6) GO TO 1002
	!check that FugiTpt did not give warning on last calls
	!pMPa=rGas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	!EAR method of Eubank and Hall (1995)
	!pMPa=pOld*( (aDepVap-aDepLiq)/(zVap-zLiq) +DLOG(zVap/zLiq)/(zVap-zLiq) )
	uSatV=DUONKT
	rhoV=pMPa*rMw(1)/(zV*rGas*tK)
	chemPot=chemPotV !since chemPotL=chemPotV at pSat       
	call Fugi(tK,pMPa,xFrac,NC,1,FUGC,zL,ier)
	uSatL=DUONKT
	rhoL=pMPa*rMw(1)/(zL*rGas*tK)
	!print*, 'tK,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode'
	!print*, tK,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode
86	continue
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE PSITER(NC,iErrCode)
	!C
	!C  PURPOSE:  ITERATIVE CALCULATION OF SOLUBILITY OF SOLVENT IN
	!C    A LIQUID POLYMER BASED ON ASSUMPTION OF ZERO POLYMER IN UPPER PHASE
	! Method:
	! wt fractions are used as the basis for computation. Molar K-Values are converted 
	! to wt fraction for flash computation then converted back for FUGI calculations.
	! Based on vlKwt, we estimate the wFrac of solvents in the polymer. Then we call
	! FUGI to get new vlKwt and compare ther new wFrac's. Iterate until self-consistent.
	!C  PROGRAMMED BY:  JRE 1/96
	!	Revised: jre 12/04
	! Example. TriEtCitrate+PLA50K at 298K, 0.1MPa, kij=0.007, alpha=0.5, SymEsd
	!          KLLwInit=1.1=>x1Wt=0.90909,x1=0.9948,chemPoVap1= -15.4544
	!          1st iteration: chemPoVap1= -15.4522
	!          259 iteration: KLLwt=2.18, x1Wt=0.4589, zLo=0.0209, zUp=0.0100
	!C
	! Notation:
	! yi, xi, zi = WT fractions !!!!
	! xMol, yMol = mole fractions
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	dimension ierFug(12)!,idPoly(NMX)
	dimension FUGC(NMX),chemPoVap(NMX),xOld(NMX),yMol(NMX),xMol(NMX)
	character*77 errMsg(12)
	DIMENSION X(NMX),Y(NMX),Z(NMX),VLK(NMX),rLnVlk(NMX)
	COMMON/eta/etaL,etaV,ZL,ZV
	MAXIT=1111
    iErrCode=0

	WRITE(*,*)'THIS ROUTINE ASSUMES ZERO POLYMER IN THE SOLVENT PHASE'
	WRITE(*,*)'USE FLASH IF xUpolymer>.001'
	write(*,*)'Solvents must be comps 1 thru nSolvComps '
	write(*,*)'Polymers: iPolyComp thru nComps '
	write(*,*)'Note: Kw=K WT BASIS '
	write(*,*)'Note: if convergence fails, try new alpha. '
	INIT=0

	write(*,*)' Enter # of Solvent components (eg. 2 for 2comp mix)'
	read(*,*)nSolvComps
	iPolyComp=nSolvComps+1
	if(nSolvComps.eq.0)then
		iErrCode=1
		errMsg(1)=' Error PsIter:nSolvComps=0'
		goto 86
	endif
	if(iPolyComp.gt.nc)then
		iErrCode=2
		errMsg(2)=' Error PsIter:iPolyComp>nComps'
		goto 86
	endif

	if(LOUD)write(*,*)'INPUT 1 FOR LIQUID SOLVENT OR 0 FOR VAPOR SOLVENT'
	READ(*,*)LIQSOL
	z(iPolyComp)=1
	if(NC.gt.iPolyComp)then
		WRITE(6,*)' Enter  wFracs of polymer comps (solvent free) '
		READ(5,*)(z(J),J=iPolyComp,NC)
	endif
1000  CONTINUE
	ier=0
	WRITE(6,*)'TEMPERATURE(KELVIN) AND PRESSURE (MPa) (-ve to quit)'
	READ(5,*)T,P
	if(T.le.0)return
	if(nSolvComps.eq.1)then
		y(1)=1
	else
		WRITE(6,*)' Enter upper weight fractions of solvent comps '
		READ(5,*)(y(J),J=1,nSolvComps)
	endif
	!convert from wFrac to xFrac
	denomPoly=0
	denomSolv=0
	do j=1,NC
		if(j.lt.iPolyComp)then
			denomSolv=denomSolv+y(j)/rMw(j)
		else
			denomPoly=denomPoly+z(j)/rMw(j)
		endif
	enddo
      IF(INIT.EQ.0)THEN
		do i=iPolyComp,NC
			y(i)=0.1e-33
		enddo
      END IF
	sumy=0
	avgMwVap=0
	do j=1,nSolvComps
		yMol(j)=y(j)/rMw(j)/denomSolv
		sumy=sumy+yMol(j)
		avgMwVap=avgMwVap+yMol(j)*rMw(j)
	enddo
	sumz=0
	avgMwPoly=0
	do j=iPolyComp,nc
		z(j)=z(j)/rMw(j)/denomPoly
		sumz=sumz+z(j)
		avgMwPoly=avgMwPoly+z(j)*rMw(j)
	enddo
	if(ABS(sumy-1).gt.1e-8)then
		if(LOUD)pause 'psiter: sumy .ne. 1'
	endif
	if(ABS(sumz-1).gt.1e-8)then
		if(LOUD)pause 'psiter: sumz .ne. 1'
	endif
            
	!C  GET FUGACITY OF SOLVENT IN THE SOLVENT-RICH PHASE
    CALl Fugi(T,P,yMol,NC,LIQSOL,FUGC,ZV,ierFug)
	do i=1,NC
		chemPoVap(i)=fugc(i)
	enddo
      
	IF(INIT.EQ.0)THEN
		write(*,*)'Enter Picard factor (e.g. 0.5)'
		read(*,*)alpha
		IF(LIQSOL)THEN
			WRITE(*,*)'INPUT INITIAL LL-Kwt FOR SOLVENTS(ie. wtUp/wtLo, e.g. 1.1) '
			READ(*,*)(vlk(i),i=1,nSolvComps)
		ELSE
			do i=1,nSolvComps
				y_x=PC(i)*10**(7*(1+ACEN(i))/3*(1-TC(i)/T))/P
				vlk(i)=y_x*avgMwPoly/2/avgMwVap
				vlk(i)=1.3
			enddo
		END IF
	END IF !INIT==0
	noUpper=1
	do i=1,nSolvComps
		X(i)=Y(i)/vlk(i)
		xOld(i)=x(i)
		IF(x(i).LT.1)noUpper=0
	enddo
	IF(noUpper)THEN
		if(LOUD)write(*,*)'ERROR IN PS, Ksolvent<1, TRY LOWER P OR HIGHER T, K1'
		GOTO 1000
	END IF
	!if(LOUD)write(*,*)'K1,X1,Y1',K1OLD,X(1),Y(1)


	!C SET BLENDING FACTOR FOR PICARD ITERATION
	!ALPHA=+0.4 !This value (-0.1) or zero is necessary for convergence criteria of C3. Octane works well
	!FYI: for TriEtCitrate+PLA50K, alpha=+0.5. Dunno other examples. 
    !Iterate until ABS(chngMax) < 1E-5 where chngMax is the max%dev in wFrac for any component.
	ITER=0
20	ITER=ITER+1
	IF(ITER.GT.2222)then
		if(LOUD)PAUSE 'TOO MANY ITERATIONS IN PSITER'
		if(LOUD)write(*,*)'Error for iter=',iter,'KwSolvents:'
		if(LOUD)write(*,'(5e12.4)')(vlk(i),i=1,nSolvComps)
		iErrCode=4
		goto 86
	endif
	sumx=0
	do i=1,nSolvComps
		if(fugc(i).gt.699)fugc(i)=699 !anything larger than 709 gives overflow. This will effectively make the soly =0.
		sumx=sumx+x(i)
	enddo
	if(sumx.ge.1)then
		errMsg(4)= 'psiter: sumxSolv.ge.1'
		if(LOUD)write(*,*)'Error for iter=',iter,'KwSolvents:'
		if(LOUD)write(*,'(5e12.4)')(vlk(i),i=1,nSolvComps)
		iErrCode=4
		goto 86
	endif
	XPOLY=1-sumx
	DO I=iPolyComp,NC
		X(I)=Z(I)*XPOLY
	enddo
	call WtFracToMolFrac(X,xMol,avgMwLiq,NC) 

	LIQ=1
	CALl Fugi(T,P,xMol,NC,LIQ,FUGC,zLo,ierFug)
	chngMax=0
	sumx=0
	do i=1,nSolvComps
		y_x=EXP(FUGC(i)-chemPoVap(i))
		vlKwt=y_x*avgMwLiq/avgMwVap
		VLK(i)=(ALPHA*vlKwt+(1-ALPHA)*vlk(i))
		X(i)=Y(i)/VLK(i)
		CHNGi=(X(i)-XOLD(I))/X(i)
		if(ABS(CHNGi).gt.chngMax)chngMax=ABS(CHNGi)
		XOLD(i)=X(i)
		sumx=sumx+x(i)
	enddo
	IF(chngMax.GT.1.D-5)GOTO 20
	if(ABS( fugc(1)-chemPoVap(1) ) > 1e-4) then
		if(LOUD)write(*,'(a,f11.4,a,f11.4)')' chemPoUp1=',chemPoVap(1),' chemPoLo1=',fugc(1)
		if(LOUD)write(*,*) 'PS Warning: Upper not equal to lower'
	endif
	noUpper=1
	do i=1,nSolvComps
		if(vlk(i).gt.1)noUpper=0
	enddo
	IF(noUpper)THEN
		errMsg(3)='PS ERROR:ALL LOWER PHASE, TRY LOWER P, HIGHER T'
		iErrCode=3
		GOTO 86
	ENDIF
            
	if(LOUD)write(*,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITER

	! NOW PERFORM SINGLE ITERATION TO ESTIMATE POLYMER IN VAPOR      
	SUMY=0
	DO I=iPolyComp,NC
		VLK(I)=EXP(FUGC(I))	  !this is later divided by exp(fugcVap)
		Y(I)=X(I)*vlk(i)	  !we expect exp(fugcVap(i)) ~ 1
		if(y(i).lt.1e-33)y(i)=1e-33	!if too small, it might mess up fugi.
		SUMY=SUMY+Y(I)
	enddo
	if(sumy.lt.1)then
		do i=1,nSolvComps
			Y(i)=y(i)*(1-SUMY)
		enddo        
	else
		errMsg(5)='ERROR IN PS, ALL VAP, TRY LOWER P?'
		iErrCode=5
		DO I=iPolyComp,NC
			Y(I)=0
		enddo
		GOTO 86
	ENDIF

	call WtFracToMolFrac(Y,yMol,avgMwVap,NC) 

	CALl Fugi(T,P,yMol,NC,LIQSOL,chemPoVap,zUp,ierFug)
	SUMY=0
	DO I=iPolyComp,NC
		rLnVlk(i)=fugc(i)-chemPoVap(i)+LOG(avgMwLiq/avgMwVap)
		VLK(I)=EXP(rLnVlk(i))
		if(VLK(I) > 0.1)VLK(i)=0.1 !Something messed up if polymer goes to upper phase.
		Y(I)=X(I)*VLK(I)
		SUMY=SUMY+Y(I)
	enddo
	if(SUMY > 1)SUMY=1 
	do i=1,nSolvComps
		rLnVlk(i)=LOG(vlk(i))
		Y(i)=Y(I)*(1-SUMY)
	enddo
	if(LOUD)write(*,603)
	if(LOUD)write(*,604)T,P,zLo,zUp
	DO I=1,NC
		WRITE(6,602)NAME(I),ID(I)
		WRITE(6,605)
		WRITE(6,606)X(I),Y(I),VLK(I),rLnVlk(i)/2.3025
		if(LOUD)write(*,*)' '
	enddo
602	FORMAT('  COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(3X,'T(K)',8X,'P(MPa)',7X,' ZL ',6X,' ZV ')
604	FORMAT(1X,2(F10.4,1X),2X,2(F8.4,4X))
605	FORMAT(3X,' wFrLo ',3X,' wFrUp ',2X,'WiUp/WiLo',3X,'log10Kw')
606	FORMAT(1X,(F7.4,1X),E11.4,2X,E11.4,1X,2(F10.4,1X),2X,2(F6.4,4X))

	!C  NOTE:  FOR REPEAT PS CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
	!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
	!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
	!C         PROGRAM BEFORE PROCEEDING.

	INIT=1
86	continue
	if(ier.ne.0)then
		if(LOUD)write(*,*)errMsg(ier)
		do i=1,nSolvComps !must re-initialize in case last guess was Ki<1 for all i
			vlk(i)=1.1
		enddo
	endif

	WRITE(6,*)'REPEAT PS CALCULATIONS? Y/N?'
	READ(*,'(A1)')ANSWER
	IF(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GOTO 1000
	init=0 !we must re-initialize if called again, esp if # of comps changes.
	RETURN
	END

!***********************************************************************
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE PXYT(NC)
	!C
	!C  PURPOSE:  COMPUTE PXY GIVEN T.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character outFile*251,tabChar*1
	DIMENSION X(NMX),Y(NMX),x1(15),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data nx,x1/15,0.001,0.01,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.99,0.999/
	tabChar=char(9) !ascii tab character from key code chart 1
    if(LOUD)then
	    if(nc.gt.2)pause 'error - this option is only for 2 compos'
    end if
	if(nc.gt.2)return
	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	open(61,file=outFile)
	open(62,file='c:\spead\speadchart.dat',STATUS='OLD',IOSTAT=ioFileBad)
	iSpead=1
	if(ioFileBad)iSpead=0
	if(iSpead)then
		write(62,'(i5,1x,a)')id(1),name(1)
		write(62,'(i5,1x,a)')id(2),name(2)
		write(62,*)'T(K)',tabChar,' P(Mpa)',tabChar,' x1',tabChar,' y1 ',tabChar,' etaLiq',tabChar,' etaVap'
		write(62,*)' '
	endif

	WRITE(6,*)'ENTER TEMPERATURE(K)   '
	READ(5,*)tK
	if(LOUD)write(*,600)
	write(61,600)
	MAXIT=1111
	zero=0

	X(1)=x1(1)
	X(2)=1-x1(1)
	ITMAX=MAXIT
	init=0
	T=tK
	CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)
	if(ier(1).ne.0)then
		!CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,1)
	endif
	CALL ErrCheck(ier,iflag)
	if(ier(1).eq.0)init=1
	write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
	DO I=2,(nx-1)
		X(1)=x1(i)
		X(2)=1-x1(i)
		ITMAX=MAXIT
		T=tK
		CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)
		if(ier(1).ne.0)then
			CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,1)
		endif
		CALL ErrCheck(ier,iflag)
		if(ier(1).eq.0)init=1
		write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
		write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
		if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
	enddo
	X(1)=x1(nx)
	X(2)=1-x1(nx)
	ITMAX=MAXIT
	T=tK
	CALL BUBPL(T,X,NC,INIT,P,ITMAX,Y,ier)
	if(ier(1).ne.0)then
		CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,1)
	endif
	CALL ErrCheck(ier,iflag)
	write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
86	continue
	print*,'OutFile=',TRIM(outFile)
	pause 'Check results in OutFile.'
	CLOSE(61)   
	close(62)
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x, 'etaL',6x, 'etaV')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')           
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))      
	RETURN
	END

!***********************************************************************
      SUBROUTINE RegXaXdIo(nComps)
	!C
	!C  PURPOSE:  Regress bonding fractions from FTIR data.
	!C  PROGRAMMED BY:  JRE 3/01
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER fileNameXa*55
	parameter(nPtMx=222)
	external DevXa
	dimension parm(2),stdErr(2),deviate(maxPts)
	common/ RegXaXdData /xFracX(nPtMx,nMx),tKelvinX(nPtMx),xAiAcceptor(nPtMx),nCompsPas,iAcceptor,jDonor
	common/ AVERR / AAPERR,RHOERR,objfun
	write(*,*)'Enter index i for Acceptor of interest, j for Donor'
	read(*,*)iAcceptor,jDonor
	write(*,*)'Enter filename for T(k),xAi,(xFrac(I),I=1,nComps)'
	read(*,*)fileNameXa
	write(*,*)'Enter guess for hIJ'
	read(*,*)bipHij 
	open(55,file=fileNameXa)
	!open(55,file=fileNameXa,ACCESS='READ',ERR=86)
	read(55,*,END=86,ERR=86)nDbaseXa !error here means no file so skip cc search.
    if(LOUD)then
	    if(nDbaseXa.gt.nPtMx)pause 'Error in RegXaXdIo. nDbase>nPtMx'
    end if
	do iPt=1,nDbaseXa
		read(55,*)tKelvinX(iPt),xAiAcceptor(iPt),(xFracX(iPt,iComp),iComp=1,nComps)
	enddo
	close(55)
	nCompsPas=nComps
	parm(1)=bipHij
	nParms=1
	FACTOR=1.d-1
	TOL=1.d-4
	call DevXa(nDbaseXa,nParms,PARM,deviate,IFLAG)
	call LmDifEz(DevXa,nDbaseXa,nParms,PARM,factor,deviate,TOL,INFO,stdErr)
	rmsErr=SQRT( ENORM(nDbaseXa,deviate) )
	if(LOUD)write(*,*)'bipHij,rmsErr=',parm(1),rmserr
	return
86	continue
	if(LOUD)write(*,*)'Error reading XA file'
	return
	end

	SUBROUTINE DevXa(nData,nParms,PARM,deviate,IFLAG)
	USE GlobConst
	USE BIPs  ! incl maxPts
	USE Assoc
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(nPtMx=222)
	dimension ier(12)
	dimension PARM(nParms),deviate(maxPts),fugc(NMX),xFrac(NMX)
	common/ RegXaXdData /xFracX(nPtMx,nMx),tKelvinX(nPtMx),xAiAcceptor(nPtMx),nComps,iAcceptor,jDonor

	HIJ(iAcceptor,jDonor)=parm(1)
	do iPt=1,nData
		tK=tKelvinX(iPt)
		pMPa=0.1
		do iComp=1,nComps
			xFrac(iComp)=xFracX(iPt,iComp)
		enddo
		LIQ=1
		call FuEsdTp(tK,pMPa,xFrac,nComps,LIQ,FUGC,Z,IER)
		deviate(iPt)=xA(1,iAcceptor)-xAiAcceptor(iPt)
	enddo
	if(iflag.eq.0)iflag=0 !dummy to avoid warning
	rmsErr=sqrt( enorm(nData,deviate) )
	return
	end


	!***********************************************************************
	SUBROUTINE RegSpeadIo(NC)
	!C
	!C  PURPOSE:  Regress Pure component paramters for SPEADMD equation of state.
	!C  PROGRAMMED BY:  AV 8/5/07
	!C	Regress Kii for SPEAD, correction for transferability deviations.
	!C
	USE GlobConst
	USE BIPs
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5),TR(NMX)
	Integer idStore(NMX)
	EXTERNAL RegSpeadDevFcn
	character inFile*251,outFile*251,SipFile*251
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/SIPs/KII(NMX,NMX),KTII(NMX,NMX)	! kii = kij(1,1) where Kij is in BIPs.
	!common/FloryWert/ vLiq(nmx)

	!store current values of nc,id.
	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	NC=1
	zero=0
	eightySix= -86
	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in c:\spead\calceos\VleData)'
	WRITE(*,*)'FIRST LINE =NPTS,id 2ND LINE = TK, PMPA'
    !FN='helium.dat'
	READ(*,'(A)')FN
	!inFile=TRIM(masterDir)//'\VpData\'//TRIM(FN)
	inFile=TRIM(FN)
	inFile='c:\spead\calceos\VleData\'//TRIM(FN)
	open(51,file=inFile)
	outFile=TRIM(masterDir)//'\Output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\Output\KijDb.txt'
	open(61,file=outFile)
	SipFile=TRIM(masterDir)//'\Input\SipSpead.txt'
	if(DEBUG)SipFile='c:\spead\calceos\Input\SipSpead.txt'
	open(72,file=SipFile)
	parm(2)=0
	print*,'Enter 1 to fit kii or 2 to fit kii and kTii'
	read*,nParms
	!nParms=2   ! fit kii and kTii
	!nParms=1  ! fit kii,  kTii = 0
    iFlag=0
	ier=0
	do while(ier.eq.0) !loop over data sets in db
		READ(51,*)NPTS,id(1)
		if(NPTS < 1)then
			ier=1 !indicates end of database
			cycle !loop to end will terminate
		endif
		nPtsBipDat=NPTS	 ! passing through BIPs module
		DO I=1,NPTS	!we only get to here if there is no error.  otherwise, cycle takes to enddo while()
			READ(51,*)tDat(I),PDAT(I)
			write(*,'(1x,f8.2,f11.5)')tDat(I),PDAT(I)
			XDAT(I)=1
			YDAT(I)=1
			DO iComp=1,NC
				TR(iComp)=TDAT(I)/TC(iComp)
			ENDDO
        enddo
        rmsBest=12345
        kiiBest=0
        ktiiBest=0
        !bVolStore=bVolCC_mol(1)
100     continue
        print*,'enter guess for kii,ktii ( 86 to stop)'
        read*,(parm(i),i=1,nParms) !,bVolCC_mol(1)
        if(parm(1) .ne. 86)then
        	call RegSpeadDevFcn(nPts,nParms,parm,deviate,IFLAG)
            rmsErr=enorm(nPts,deviate)/nPts
            if(rmsErr < rmsBest)then
                rmsBest=rmsErr
                kiiBest=parm(1)
                ktiiBest=parm(2)
                !bVolBest=bVolCC_mol(1)
            endif
            goto 100
        endif
		parm(1)= kiiBest
		parm(2)= ktiiBest
        !bVolCC_mol(1)=bVolBest
 		MAXIT=1111
		INIT=1
		!C	LWA>5*nParms+NPTS ~ 275
		LWA=275
		!c      TOL=SQRT(DPMPAR(1))
		TOL=0.00001
		!C
		!c  factor controls the magnitude of the first step from the initial guess.
		factor=0.01
		write(*,*)'KII,paadp,rmserr'
		CALL LMDifEZ(RegSpeadDevFcn,nPts,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
		rmsErr=SQRT( ENORM(nPts,deviate) )
		write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		write(* ,'(2e12.4)')(parm(i),i=1,nParms)
		write(*,*)'StdErr'
		write(* ,'(2e12.4)')(stdErr(i),i=1,nParms)
		write(*,*) 'RegSpeadIo: Converged'
		L=1
		call KIIVAL(nc,nParms,parm,deviate,PAADP,pBias,pErrMax,L)
		print*,'KiiVal done.'
		write(* ,607)id(1),(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,pBias,pErrMax,rmsErr,name(1)
		write(61,607)id(1),(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,pBias,pErrMax,rmsErr,name(1)
		write(72,608)id(1),(parm(i),i=1,nParms),name(1)
        exit
	enddo !while(ier.eq.0) loop over data sets in db

	close(51)
	close(61)
	close(72)
	print*,'Success! Output in: ', TRIM(outFile)
	pause 'These results are tabulated in the output file.'
	!restore stored nc,id
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
607	FORMAT(1X,i6,1X,<nParms>F10.7,1X,<nParms>F10.7,f10.2,f10.6,f10.2,f10.2,1x,a20)
608 FORMAT(1X,I6,1X,<nParms>F10.4,1X,A20)
  
	RETURN
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SCFITER(NC)
!C
!C  PURPOSE:  ITERATIVE Supercritical Fluid FUGACITY EVALUATIONS
!C  PROGRAMMED BY:  JRE 2/93
!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*64 Input_File_Name
	DIMENSION X(NMX),fugcL(NMX),ier(12)
	!Dimension fugcV(NMX),Y(NMX),VLK(NMX)
	COMMON/eta/etaL,etaV,ZL,ZV
	common/eta2/eta
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	WRITE(*,*)'ENTER vSolid(cc/mol) for Poynting correction.'
	read(*,*)vSolidCc_mol
	WRITE(52,*)'ENTER vSolid(cc/mol) for Poynting correction.'
	write(52,*)vSolidCc_mol
	WRITE(*,*)'ENTER A,B FOR SOLID SUB PRESSURE:ln(Pmpa)=A-B/TK'
	read(*,*)vpaSolid, vpbSolid
	WRITE(52,*)'ENTER A,B FOR SOLID SUB PRESSURE:ln(Pmpa)=A-B/TK'
	write(52,*)vpaSolid, vpbSolid
	y1CalcOld=0
		
	WRITE(*,*)'What is the file name >:(?'
	WRITE(*,*)'The format of the file is thusly:'
	WRITE(*,*)'First Line:	# of data pointes <tab> # of components'
	WRITE(*,*)'Following Lines are data points:'
	WRITE(*,*)'Temp(K)	Pres(MPa)	xSolid xFirstSolvent xSecondSolvent etc...'
	READ(*,*)Input_File_Name
	WRITE(52,*)'What is the file name >:(?'
	write(52,*)Input_File_Name
	open(69, file="SC_Results.txt")
	!WRITE(79,607)etaL,ZL,fugLi,fugcL,(X(J),J=1,NC),&			fugcL(J),J=1,NC)	
	WRITE(69,*)'T	P	etaL	zL	x(1)'
	open(68, file=Input_File_Name)	
	read(68,*)nMaxIter
!1000  CONTINUE	Disabled when I started reading data from the file
	DO Iter_File_Data=1,nMaxIter
	
	ierScf=0
	!      WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND PRESSURE(MPa)   '
	!      READ(5,*)T,P	Disabled when I started reading data from the file
	READ(68,*)T,P,(X(J),J=1,NC)
	pSatSolid=exp(vpaSolid-vpbSolid/T)
	poynting=exp( (P-pSatSolid)*vSolidCc_mol/rGas/T)
 
	!     WRITE(6,*)'ENTER Guess for SCF COMPOSITIONS (comp1 = solid)'
	!     READ(5,*)(X(J),J=1,NC)	Disabled when I started reading data from the file
	dev=1111
	iter=0
	do while(ierScf.eq.0 .and.abs(dev).gt.1.d-5)
		iter=iter+1
		if(iter.gt.111)ierScf=2
		IFLAG=0
		!CALl Fugi(T,P,X,NC,1,fugcL,ZL,ier)
		CALl Fugi(T,P,X,NC,1,fugcL,ZL,ier)
		DO I=1,12
			IF (ier(I).NE.0)IFLAG=1
		enddo
		IF (IFLAG.EQ.1)then
			WRITE(6,'(a,22I4)')'ier',(ier(I),I=1,12)
			WRITE(6,*)'error in lower phase calculation'
			ierScf=1
		endif
		enhancement= poynting/exp( fugcL(1)	)
		x1Calc=pSatSolid/P*enhancement
		dev=x1Calc-x(1)
		if(LOUD)write(*,*)'iter      y(1),       y1Calc         E       poynt'
		if(LOUD)write(*,'(i4,4e12.4)')iter,x(1),x1Calc,enhancement,poynting
		!pause
		if(x1Calc.gt.0 .and. x1Calc .lt. 1)then
			x(1)=x1Calc
			sumup=0
			do iComp=2,NC
				sumup=sumup+x(iComp)
			enddo
			do iComp=2,NC
				x(iComp)=x(iComp)/sumup*(1-x1Calc)
			enddo
		endif
	enddo !while

	if(ierScf.eq.1 .and. LOUD)write(*,*)'Error returned from FUGI.'
	if(ierScf.eq.2 .and. LOUD)write(*,*)'Error: too many iterations.'

	DHL=DHONKT
	DSL=DSONK
	!for Fugi:
	!C  ier = 1 - AT LEAST ONE ERROR
	!C        2 - NOT USED
	!C        3 - NOT USED
	!C        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!C        5 - rho IS -VE
	!C        6 - TOO MANY Z ITERATIONS
	!C        7 - eta > 0.53
	!C		11 - goldenZ instead of real Z.
	IF(ier(2).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL UPPER PHASE'
	IF(ier(3).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL LOWER PHASE'
	if(ierScf.eq.0 .and. LOUD)then
		DO I=1,1
			WRITE(*,602)NAME(I),ID(I)
			WRITE(69,*)T,P,etaL,ZL,X(1)
			WRITE(*,603)
			fugLi=X(i)*p*EXP(fugcL(i))
			WRITE(*,604)X(I),fugLi,fugcL(I)
			WRITE(*,*)'  '
			WRITE(*,605)
			WRITE(*,606)T,P,etaL,ZL
			write(*,*)' '
			!if(LOUD)PAUSE
		enddo    
		WRITE(*,*)'DEPARTURE FUNCTIONS: (Hscf-Hig)/NkT  (Sscf-Sig)/Nk'
		WRITE(*,610)DHL,DSL
	endif

602   FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603   FORMAT(4X,'  yFrac',3X,'  Fugacity(MPa)',5X,'  lnPhi')
604   FORMAT(1(F12.10),1(F12.10),1(F12.10))
!FORMAT(1(F12.10),2(5x,e11.4),2(1x,E11.4))
605   FORMAT(3X,'T(K)   ',4X,'P(MPa)',3X,'eta',6X,'zFactor')
606   FORMAT(f7.2,3X,F9.3,1X,F8.5,1X,2(F8.5,3X),2(F9.4,1X))
607	FORMAT(f4.3,' ',f7.5,' ',2(f12.9,' '),4(f12.9,' '))
610   FORMAT(21X,4(F10.4,4X))      
!	Disabled when I started reading data from the file
!     WRITE(6,*)'REPEAT SCF CALCULATIONS? Y/N?'
!     READ(*,'(A1)')ANSWER
!      IF(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GOTO 1000
!	Disabled when I started reading data from the file
4412  END DO
	CLOSE(68)
	CLOSE(69)
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE TEITER(NC)
	!
	!  PURPOSE:  ITERATIVE THREE PHASE EVALUATIONS
	!  PROGRAMMED BY:  JRE 2/93
	!
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION S(NMX),KI(NMX)
	DIMENSION X(NMX),X2(NMX),Y(NMX),ier(20)
	DIMENSION ZFEED(NMX)
	COMMON/eta/etaL,etaU,ZL,ZU
	COMMON/FEED/ZFEED
	COMMON/KVALUES/S
	MAXIT=1111
	P=10
	WRITE(*,*)'INITIALIZE PRESSURE FROM RAOULTS LAW? Y/N? '
	READ(*,'(A1)')ANSWER
	INIT=0
	IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
		WRITE(*,*)'ENTER INITIAL GUESS FOR PRESSURE(MPA)'
		READ(*,*)P
		INIT=1
	ENDIF
	WRITE(*,*)'ENTER GUESSES FOR LIQ xUi/xLi  '
	READ(*,*)(KI(I),I=1,NC)
	do J=1,nc
		S(J)=KI(J)
	ENDDO        

1000 CONTINUE
	WRITE(*,*)'ENTER TEMPERATURE(KELVIN) '
	READ(*,*)T
	WRITE(*,*)' '
	WRITE(*,*)'ENTER FEED COMPOSITIONS      '
	READ(*,*)(ZFEED(I),I=1,NC)

	IFLAG=0
	ITMAX=MAXIT
	CALL LLVE(T,zFeed,NC,INIT,P,ITMAX,Y,X,X2,ier,UOF,S)      

	DO I=1,6
		IF (ier(I).NE.0)IFLAG=1
	ENDDO
	IF (IFLAG.EQ.1 .and. LOUD)WRITE(6,'(6(1X,I2))')'ier',(ier(I),I=1,6)
	IF(ier(2).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL UPPER PHASE'
	IF(ier(3).EQ.1 .and. LOUD)WRITE(*,*)'ERROR ALL LOWER PHASE'

	if(LOUD)write(*,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX


	WRITE(6,*)'       T(K)       P(MPa)          ZL           ZU'
	WRITE(6,'(1X,4(1X,F11.4))')T,P,ZL,ZU
	DO I=1,NC
		WRITE(6,182)NAME(I),ID(I)
		WRITE(6,*)' THE L1,L2,V COMPOSITIONS & Yi/XiS ARE: '
		WRITE(6,183)X(I),X2(I),Y(I),S(I)
		WRITE(6,*)'  '
	ENDDO
	WRITE(6,188)UOF

	!  NOTE:  FOR REPEAT TP CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
	!         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
	!         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
	!         PROGRAM BEFORE PROCEEDING.

	INIT=1

	WRITE(6,*)'REPEAT TP CALCULATIONS? Y/N?'
	READ(*,'(A1)')ANSWER
	IF(ANSWER.NE.'N'.AND.ANSWER.NE.'n')GO TO 1000
182	FORMAT(2X,'COMPONENT IS ',A20,2X,'ID NO. IS ',1X,i5)
183	FORMAT(2X,'X= ',E11.4,2X,'X2= ',E11.4,2X,'Y= ',E11.4,2X,'Yi/Xi= ',E11.4)
188	FORMAT(2X,'UOF =  ',F7.5)
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE TXYP(NC)
	!C
	!C  PURPOSE:  COMPUTE TXY GIVEN P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character outFile*251,tabChar*1
	DIMENSION X(NMX),Y(NMX),x1(15),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data nx,x1/15,0.001,0.01,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.99,0.999/
	!data nx,x1/13,0.99,0.95,0.9,.8,.7,.6,.5,.4,.3,.2,.1,0.05,0.01/
	tabChar=char(9) !ascii tab character
	if(nc.ne.2 .and. LOUD)pause 'error - this option is only for 2 compos'
	if(nc.ne.2)return
	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	open(61,file=outFile)
	open(62,file='c:\spead\speadchart.dat',STATUS='OLD',IOSTAT=ioFileBad)
	iSpead=1
	if(ioFileBad)iSpead=0
	if(iSpead)then
		write(62,'(i5,1x,a)')id(1),name(1)
		write(62,'(i5,1x,a)')id(2),name(2)
		write(62,*)'T(K)',tabChar,' P(Mpa)',tabChar,' x1',tabChar,' y1 ',tabChar,' etaLiq',tabChar,' etaVap'
		write(62,*)' '
	endif

	WRITE(6,*)'ENTER PRESSURE(MPa)   '
	READ(5,*)pMPa
	write(*,600)
	write(61,600)
	MAXIT=1111
	zero=0

	X(1)=x1(1)
	X(2)=1-x1(1)
	ITMAX=MAXIT
	P=pMPa
	CALL BUBTL(T,X,NC,INIT,P,ITMAX,Y,ier)
	if(ier(1).ne.0)then
		CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
	endif
	CALL ErrCheck(ier,iflag)
	write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
	if(ier(1).eq.0)init=1

	DO I=2,(nx-1)
		X(1)=x1(i)
		X(2)=1-x1(i)
		ITMAX=MAXIT
		bipKIJ=KIJ(1,2)+KTIJ(1,2)/t
		!write(*,*)'KIJ=',bipKIJ
		P=pMPa
		CALL BUBTL(T,X,NC,INIT,P,ITMAX,Y,ier)
		if(ier(1).ne.0)then
			CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
		endif
		CALL ErrCheck(ier,iflag)
		write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
		write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
		if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
		if(ier(1).eq.0)init=1
	enddo 
	X(1)=x1(nx)
	X(2)=1-x1(nx)
	ITMAX=MAXIT
	P=pMPa
	CALL BUBTL(T,X,NC,INIT,P,ITMAX,Y,ier)
	if(ier(1).ne.0)then
		CALL BootBP(T,X,NC,INIT,P,ITMAX,Y,ier,2)
	endif
	CALL ErrCheck(ier,iflag)
	if(LOUD)write(* ,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	write(61,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)

86	continue
	print*,'OutFile=',TRIM(outFile)
	pause 'Check results in OutFile.'
	CLOSE(61)   
	close(62)
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x, 'etaL',6x, 'etaV')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')           
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))      
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE VLITER(NC)
	!C  PURPOSE:  ITERATIVE VLFLASH EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*1 ANSWER
	DIMENSION X(NMX),Y(NMX)!,VLK(NMX)
	DIMENSION ZFEED(NMX),CALCK(NMX)
	COMMON/eta/etaL,etaV,ZL,ZV

	WRITE(*,*)'INITIALIZE K-RATIOS FROM VAPOR PRESSURE EQUATION? Y/N?'
	READ(*,'(A1)')ANSWER
	WRITE(52,*)'INITIALIZE K-RATIOS FROM VAPOR PRESSURE EQUATION? Y/N?'
	write(52,'(A1)')ANSWER
	INIT=0
	IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
		WRITE(*,*)'ENTER INITIAL GUESS FOR Yi/XiS '
		READ(*,*)(CALCK(I),I=1,NC)
		WRITE(52,*)'ENTER INITIAL GUESS FOR Yi/XiS '
		write(52,*)(CALCK(I),I=1,NC)
		INIT=1
	END IF
	!open(61,file='TPXY.txt')
	!write(61,*)' '
	!close(61)
	ANSWER='Y'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
		MAXIT=1111
		WRITE(6,*)'ENTER TEMPERATURE(KELVIN) & PRESSURE(MPa)  '
		READ(5,*)T,P
		WRITE(6,*)' '
		WRITE(6,*)'ENTER FEED COMPOSITION   '
		READ(5,*)(ZFEED(I),I=1,NC)
		WRITE(52,*)'ENTER TEMPERATURE(KELVIN) & PRESSURE(MPa)  '
		write(52,*)T,P
		WRITE(52,*)' '
		WRITE(52,*)'ENTER FEED COMPOSITION   '
		write(52,*)(ZFEED(I),I=1,NC)

		IFLAG=0
		ITMAX=MAXIT
		CALL VLEFL(T,P,NC,INIT,ZFEED,CALCK,ITMAX,X,Y,VOF,iErrCode)

		IF(iErrCode.ne.0)then
			call PrintErrFlash(iErrCode)
		else
			!if(iErrCode.ne.0)call PrintErrFlash(iErrCode)
			call FlashOut(NC,X,Y,T,P,itmax,vof) !takes z,eta from common eta
		endif

		!C  NOTE:  FOR REPEAT VL CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!C         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!C         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!C         PROGRAM BEFORE PROCEEDING.

		INIT=1

		WRITE(6,*)'REPEAT VL CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
		WRITE(52,*)'REPEAT VL CALCULATIONS? Y/N?'
		write(52,'(A1)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
	if(LOUD)pause 'Results are in TPXY.txt'
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Programmed by AGF 7/10
	!The routine calculates vapor pressure from TrMin to 0.95 then seeks user input.
	subroutine VpIter(NC)
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	character outFile*251
	DIMENSION gmol(NMX)
	LOGICAL LOUDER
	COMMON/eta/etaL,etaV,zLiqDum,zVapDum
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	common/GaussFun/G0,G1,eta0
	LOUDER=LOUD
	LOUDER=.TRUE.
	iComp=1
	if(NC > 1 .and. LOUDER)then
		write(*,*)'VpIter: Only works for a single component. Comp1 or Comp2?'
		read*,iComp
	endif
	gmol(1:NC)=0
	gmol(iComp)=1
	write(*,*)'First calculate curve from TrMin-Tr=0.95, then user input.' 
	write(*,'(a,f7.2,a)')' Tc=',TC(1),'K. Enter TrMin (e.g. 0.45)'
	read(*,*) TrMin
	nSteps=10
	delTr=(0.95d0-TrMin)/nSteps
	iTrMin=INT(TrMin*100/5)*5
	outFile=TRIM(masterDir)//'\output\Binodal.txt'
	open(61,file=outFile)
	if(isTPT)write(6 ,'(a,3F9.3)')' G0,G1,e0 ',G0,G1,eta0
	if(isTPT)write(61,'(a,3F9.3)')'G0,G1,e0 ',G0,G1,eta0
	write(*,*)'  T(K)     PsatMPa     rhoVap(g/cc)    rhoLiq(g/cc)'
	write(61,*)'  T(K)    PsatMPa   rhoVap    rhoLiq(g/cc)'
	do iStep=1,nSteps
		tKelvin=(TrMin+delTr*(iStep-1))*Tc(1)
		call PsatEar(tKelvin,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode)
		if(ierCode.ne.0)then
			if(LOUDER)pause 'Psat returned error code'
			exit
		endif
		!call EqualArea(NC,gmol,tKelvin,Psat,vLCc,vVCc)
		write(* ,'(F8.2,1x,3F9.6,i3)')tKelvin,pMPa,rhoV,rhoL,ierCode
		write(61,'(F8.2,1x,3F9.6,i3)')tKelvin,pMPa,rhoV,rhoL,ierCode
	enddo
	tKelvin=0.96*Tc(1)
	do while(tKelvin.gt.0)
		call PsatEar(tKelvin,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode)
		!call EqualArea(NC,gmol,tKelvin,Psat,vLCc,vVCc)
		write(*,'(F8.2,1x,3F9.6,i3)')tKelvin,pMPa,rhoV,rhoL,ierCode
		if(pMPa > 1.D-3)then
			write(61,'(F8.2,1x,3F9.6,i3)')tKelvin,pMPa,rhoV,rhoL,ierCode
		else
			write(61,'(F8.2,1x,3(E11.4,1x),i3)')tKelvin,pMPa,rhoV,rhoL,ierCode
		endif
		write(*,*)'ENTER Temperature (K). (-ve to stop) '
		READ(*,*) tKelvin
	enddo
	rhoc=Pc(1)*rMw(1)/( Zc(1)*8.31434*Tc(1) )
	etac=rhoc/rMw(1)*bVolCc_Mol(1)
	write(* ,'(F8.2,1x,5F9.6)')tc(1),pc(1),rhoc,rhoc,Zc(1),etac
	write(61,'(F8.2,1x,5F9.6)')tc(1),pc(1),rhoc,rhoc,Zc(1),etac
	close(61) 
	if(LOUDER)write(*,*) 'Your results are in ',TRIM(outFile)
	if(LOUDER)pause
86	return 
	end

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine WtFracToMolFrac(wFrac,xMol,avgMw,nComps)
	USE GlobConst !for rMwPlus
	implicit DoublePrecision(a-h,o-z)
	dimension wFrac(nComps),xMol(nComps)
	denom=0
	do i=1,nComps
		denom=denom+wFrac(i)/rMw(i)
	enddo
	avgMw=0
	do i=1,nComps
		xMol(i)=wFrac(i)/rMw(i)/denom
		avgMw=avgMw+xMol(i)*rMw(i)
	enddo
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine SortTable(table,maxRows,nRows,nCols,keyCol,LoToHi)
	! Purpose: Sort according to keyCol using IMSL routine.
	! 
	USE MSIMSL !for DSROWR()
	Implicit DoublePrecision(a-h,o-z)
	DoublePrecision	table(maxRows,nCols)
	Integer iPerm(maxRows),NI(1),INDKEY(1),NRMISS
	iAbsolute=0	! sort by algebraic values not absolute values
	if(LoToHi > 1 .or. LoToHi < 0)pause 'SortTable: LoToHi must be 0 or 1.'
	iOrder=1-LoToHi	! IMSL swaps the order specifier.
	iRet=2 ! sorted table is returned but not NGROUP or NI. (iRet=0 to return NGROUP and NI).
	nKey=1 ! just sort by one column.
	INDKEY(1)=keyCol
	CALL DSROWR(nRows,nCols,table,maxRows,iAbsolute,IORDR,iRet,nKey,INDKEY,IPERM,NGROUP,NI,NRMISS)!IMSL routine
!LDX — Leading dimension of table exactly as specified in the dimension statement of the calling program.   (Input)
	return
	end
! CALL SROWR (NROW, NCOL, X, LDX, ICOMP, IORDR, IRET, NKEY,INDKEY, IPERM, NGROUP, NI, NRMISS)
! NROW — Number of rows of X.   (Input)
! NCOL — Number of columns of X.   (Input)
! X — NROW by NCOL matrix.   (Input, if IRET = 1; input/output if IRET = 0)
! On input, X contains the matrix to be sorted. If IRET = 0, the output X contains the sorted matrix.
! LDX — Leading dimension of X exactly as specified in the dimension statement of the calling program.   (Input)
! ICOMP — Option giving the method of comparison of the row vectors.   (Input) 
! ICOMP	Action
! 0	Elementwise, by algebraic values
! 1	Elementwise, by absolute values
! IORDR — Option giving the sorting order.   (Input)
! IORDR	Action
! 0	Ascending
! 1	Descending
! IRET — Option to indicate information returned.   (Input)
! IRET	Action
! 0	The sorted X is returned along with NGROUP and NI.
! 1	X is unchanged (detached key sort) and NGROUP and NI are returned.
! 2	The sorted X is returned, but NGROUP and NI are not returned.
! 3	X is unchanged (detached key sort) and NGROUP and NI are not returned.
! NKEY — Number of columns of X on which to sort.   (Input)
! INDKEY — Vector of length NKEY giving the column numbers of X which are to be used in the sort.   (Input)
! IPERM — Permutation vector of length NROW specifying the rearrangement of the rows.   (Output)
! IPERM(I) = J means row I of the sorted X is row J of the unsorted X.
! NGROUP — Number of groups.   (Output, if IRET £ 1)
! The rows of the sorted X are partitioned into groups. A group contains rows that are equal with respect to the method of comparison. NGROUP is the number of groups of different rows.
! NI — Vector of length NGROUP containing the number of rows in each group.   (Output, if IRET £ 1)
! The first NI(1) rows of the sorted X are group number 1. The next NI(2) rows of the sorted X are group number 2. ... The last NI(NGROUP) rows of the sorted X are group number NGROUP. If NGROUP is not known prior to the invocation of this routine, NROW(an upper bound for NGROUP) can be used as the dimension of NI. If IRET ³ 2, NI is not referenced and can be a vector of length one.
! NRMISS — Number of rows that contained NaN in the columns of X used in the sort.   (Output)
! These rows are considered as a separate group from the other NGROUP groups and are put as the last NRMISS rows of the sorted X.
! 
