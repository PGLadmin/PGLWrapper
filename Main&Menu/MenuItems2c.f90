	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ErrCheck(ier,iflag)
	DIMENSION ier(20)
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
	open(641,file=outFile,ACCESS='APPEND')
	if(LOUDER)write(* ,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX
	if(LOUDER)write(* ,'(4(5x,a))')' zLo ',' zUp ','etaLo','etaUp'
	if(LOUDER)write(* ,'(4(f9.6,1x))')zL,zV,etaL,etaV
	if(LOUDER)write(* ,'(a,f8.2,a,f9.4)')' T(K)= ',tKelvin,' P(MPa)=',pMpa
	if(pMpa.lt.1E-3 .or. pMpa > 10E4  .and. LOUDER)write(* ,*)'PMPa=',pMpa
	WRITE(641,*)'REQUIRED NUMBER OF ITERATIONS WAS:',ITMAX
	write(641,'(4(5x,a))')' zLo ',' zUp ','etaLo','etaUp'
	write(641,'(4(f9.6,1x))')zL,zV,etaL,etaV
	write(641,'(a,f8.2,a,f9.4)')' T(K)= ',tKelvin,' P(MPa)=',pMpa
	if(pMpa.lt.1E-3 .or. pMpa.gt.10E4)write(641,*)'PMPa=',pMpa
	avgMwLiq=0
	avgMwVap=0
	do i=1,nComps
		avgMwLiq=avgMwLiq+xFrac(i)*rMw(i)
		avgMwVap=avgMwVap+yFrac(i)*rMw(i)
	enddo

	DO I=1,nComps
		if(LOUDER)write(* ,602)NAME(I),ID(I)
		WRITE(641,602)NAME(I),ID(I)
		VLKI=yFrac(I)/(xFrac(I)+1.D-11)
		wFrVap=yFrac(i)*rMw(i)/avgMwVap
		wFrLiq=xFrac(i)*rMw(i)/avgMwLiq
		if(LOUDER)write(* ,603)
		if(xFrac(i).gt.1E-7)then
			if(LOUDER)write(* ,604)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
			WRITE(641,604)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
		else
			if(LOUDER)write(* ,605)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
			WRITE(641,605)xFrac(I),yFrac(I),VLKI,wFrLiq,wFrVap
		endif
		if(LOUDER)write(* ,*)' '					   
		WRITE(641,603)
		WRITE(641,*)' '					   
	enddo
	if(LOUDER)write(* ,182)VOF
	WRITE(641,182)VOF
	close(641)
182	FORMAT(2X,'UpperMoles/Feed = ',F11.5)
602	FORMAT('  COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(3X,'  LIQUID X  ',3X,'Upper Y',3X,' Yi/Xi ',5X,' wFrLo ',6X,'wFrUp')
604	FORMAT(1X,(F11.8,1X),2(E11.4,1X),F11.8,1X,E11.4)
605	FORMAT(1X,(E11.4,1X),2(E11.4,1X),E11.4,1X,E11.4)
	return
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine BipIo(NC) !note: results passed through USEd BIPs
	USE GlobConst
	USE BIPs
	USE Assoc, only:idLocalType,nTypesTot,aBipAD,aBipDA
	!Integer idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	!DoublePrecision aBipAD(maxTypes,maxTypes),aBipDA(maxTypes,maxTypes) !association bips
	implicit DoublePrecision(a-h,k,o-z)
	character answer*1
	LOGICAL isNRTL
	isNRTL=.FALSE.
	if(iEosOpt==3 .or. iEosOpt==7)isNRTL=.TRUE. 
	DO J=2,NC
		DO I=1,J-1
			if(isNRTL .or. iEosOpt==19)then
				WRITE(*,*)i,j,' - ',j,i
				READ(*,*)xsTau(I,J),xsTau(J,I)
				WRITE(52,*)i,j,'  -  ',j,i
				write(52,*)xsTau(I,J),xsTau(J,I)
				if(isNRTL)then
					xsAlpha(I,J)=0.3d0
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
					xsAlpha(J,I)=xsAlpha(I,J)
				endif
			endif
			if(iEosOpt==7 .or. iEosOpt==19)cycle
			WRITE(*,*)'ENTER KIJ,KTIJ FOR THE PAIR',I,J,'    '
			READ(*,*)KIJ(I,J),KTIJ(I,J)
			WRITE(52,*)' KIJ,KTIJ FOR THE PAIR',I,J,'    '
			write(52,*)KIJ(I,J),KTIJ(I,J)
			KIJ(J,I)=KIJ(I,J)
			KTIJ(J,I)=KTIJ(I,J)
		enddo
		if(LOUD)write(*,601)J,ID(J),(KIJ(J,I),I=1,J-1)
	enddo !J=2,NC

	if(iEosOpt.ge.86.and.itsEven(iEosOpt))then !86=4. disable this option in favor of MEM2 (for now). JRE 20230313
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

	if(iEosOpt==86)then	 !86=5. disable this option in favor of MEM2 (for now). JRE 20230313
		write(*,'(5x,11i8)')(idLocalType(j),j=1,nTypesTot)
		do i=1,nTypesTot
			write(*,'(i5,11f8.4)')idLocalType(i),(aBipDA(i,j),j=1,nTypesTot)
		enddo
		print*,'Enter i,j,BipDA(i,j)'
		read*,i,j,aBipDA(i,j)
		aBipAD(j,i)=aBipDA(i,j)
		aBipAD(i,j)=aBipDA(i,j)	!keep symmetric for now. JRE 20220108
		aBipDA(j,i)=aBipDA(i,j)
	elseif(iEosOpt==10)then
		print*,'BIPio: setting pcSaft. Kij=',KIJ(1,2)
		Call SetParMixPcSaft(1,KIJ(1,2),iErrBip)	
    endif
601 FORMAT(I3,I5,9(1X,F7.4))
602	format(' ENTER xsTau(',i2,',',i2,'),xsTau(',i2,',',i2,') OR ZEROES TO NULLIFY XS EFFECT')
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	DIMENSION X(NMX),FUGC(NMX),ier(20)
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
	open(632,file=outFile)
	vOld1=0
	vOld2=0
	LIQ=2
	do it=130,1,-1
		iBifurcate=0
		tKelvin=it*TC(1)/100
		Call FugiTP( tKelvin,PMPa,X,NC,LIQ,rhoMol_cc,zFactor,aRes,fugc,uRes,iErrF )
		!CALL FUGI(tKelvin,PMPa,X,NC,LIQ,FUGC,zFactor,ier)
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
			write(632,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
			write(* ,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
		else
			write(632,601)it,tKelvin,vNew,ier(1)
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
		Call FugiTP( tKelvin,PMPa,X,NC,LIQ,rhoMol_cc,zFactor,aRes,fugc,uRes,iErrF )
		ier(1)=iErrF
		!CALL FUGI(tKelvin,PMPa,X,NC,LIQ,FUGC,zFactor,ier)
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
			write(632,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
			write(* ,601)it,tKelvin,vNew,ier(1), ' Probable bifurcation.'
		else
			write(632,601)it,tKelvin,vNew,ier(1)
			write(* ,601)it,tKelvin,vNew,ier(1)
		endif
		vOld2=vOld1
		vOld1=vNew
	enddo
	close(632)
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
	open(633,file=outFile)
600	format('   T(K)',4x,'  P(MPA)',5x,'x1Alpha',4x,'x1Beta',2x,'etaAlpha',3x, 'etaBeta')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),6(F9.5,1X))
602	FORMAT(' COMPONENT IS',2X,A20,2X,'ID NO. IS  ',i5)
603	FORMAT(4X,'LIQUID X',4X,'VAPOR Y',3X,'fugcL(MPa)',4X,'fugcV(MPa)')           
604	FORMAT(4X,2(F7.4,4X),2(F7.4,8X))
605 format(1X,F7.2,2X,F10.6,1X,4(f13.11,1X))
606	FORMAT(F12.5,3X,F7.2,1X,F7.4,1X,4(F6.4,3X))
610	FORMAT(21X,4(F10.4,4X))      
	WRITE(6,*)'ENTER PRESSURE(MPa) and starting Temperature (K)'
	READ(5,*)P,T
	if(LOUD)write(*,600)
	write(633,600)
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
			write(633,605)T,P,X(1),xU(1),etaPasLo,etaPasUp
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
	CLOSE(633)   
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

	dimension xFrac(NMX),yFrac(NMX),ier(20)
	dimension idFile(NMX)
	!common/tpxyData/tDat(111),PDAT(111),XDAT(111),npts
	COMMON/eta/etaL,etaV,ZL,ZV

	if(LOUD)write(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (..\VleData\FN.FT)'
	if(LOUD)write(*,*)'FIRST LINE =NPTS,nComps,(id(i),i=1,nComps)'
	if(LOUD)write(*,*)'2ND LINE = TK, PMPA, X1,x2,x3,...'
	READ(*,'(A)')FN
	if(DEBUG)then
		FN=TRIM(PGLinputDir)//TRIM(FN)
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

SUBROUTINE DifCoDb
	!C
	!C  PURPOSE:  ITERATIVE DEWT EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/93
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE CritParmsDb
	Implicit DoublePrecision(A-H,O-Z)
	!CHARACTER*1 ANSWER
	Character*255 inFile,outFile,dumString
	DoublePrecision chemPo(NMX),gmol(NMX)
	inFile=TRIM(PGLInputDir)//'\DifCoDbG&E2010IecrSi.txt'
	open(51,file=inFile)
	outFile=TRIM(PGLInputDir)//'\DifCoOut.txt'
	open(61,file=outFile)
	nComps=1
	gmol(1)=1
	isZiter=0
	nPts=1
	read(51,'(a255)',ioStat=ioErr)dumString
	if(ioErr /= 0)write(*,*)'DifCoDb: Error reading first line of ',TRIM(inFile)
	write(61,*)'ID1,tKelvin,pMpa,rhoG_cc,difCoCm2_s,sRes,pCalc'
	idOld=0
	do while(ioErr==0)
		read(51,'(a255)',ioStat=ioErr)dumString
		if(ioErr /= 0)then
			if(ioErr /= -1)write(*,*)'DifCoDb: error reading line=',TRIM(dumString)
			exit	! end of file
		endif
		read(dumString,*,ioStat=ioErr)ID(1),tKelvin,pMpa,rhoG_cc,difCoCm2_s
		if(ioErr /= 0)then
			write(*,*)'DifCoDb: error reading line=',TRIM(dumString)
			exit	! end of file
		endif
		Call PGLStartup(nComps,5,1,iErrStart)	! iEosOpt=5 means SPEADMD, idOpt=1 means ID(dippr) is USEd.
		if(iErrStart /= 0)then
			write(*,*)' DifCoDb: from PGLStartup, iErrStart,ID(1)=',iErrStart,ID(1)
			cycle
		endif
		if(idOld /= ID(1))then
			write(*,*)' DifCoDb: new ID=',ID(1)
			idOld=ID(1)
		endif
		rhoMol_cc=rhoG_cc/rMwD( CrIndex(ID(1)) )
		vTotCc=1/rhoMol_cc
		call FuTptVtot(isZiter,zFactor,aRes,uRes,chemPo,vTotCc,tKelvin,gmol,nComps,iErr)
		sRes=uRes-aRes
		pCalc=zFactor*rGas*rhoMol_cc*tKelvin
		write(61,form601)ID(1),tKelvin,pMpa,rhoG_cc,difCoCm2_s,sRes,pCalc
	enddo
	close(51)
	close(61)
	write(*,*)'DifCoDb: outFile=',TRIM(outFile)
	pause 'DifCoDb: Success! Check output file.'
	return
end	SUBROUTINE DifCoDb

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
	!Integer ier(20)
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
    x(1)=1
	WRITE(6,*)'ENTER Lower Phase COMPOSITIONS '
	if(NC > 1)READ(5,*)(X(J),J=1,NC)
	WRITE(52,*)'LIQD COMPOSITIONS '
	write(52,*)(X(J),J=1,NC)
	bMix  =SumProduct( NC,x,bVolCC_mol(1:NC) )
	rMwAvg=SumProduct( NC,x,rMw(1:NC) )
	wFrac(1:NC)=x(1:NC)*rMw(1:NC)/rMwAvg
	write(*,*)'Wt Fracs'
	write(*,'(11f8.5)')(wFrac(i),i=1,nc)
	if(initCall)print*,'Fiter: calling Fugi for liq. T(K)=',T
	CALL FugiTP( T,P,X,NC,LIQ,rhoLo,ZL,aRes,fugcL,uRes,iFlag )
	Hres_RT=uRes+ZL-1
	Sres_R=aRes-uRes
	if(initCall)write(*,'(a,f8.2,E12.4)')' Fiter: called Fugi for liq. T(K),ZL=',T,ZL
	IF (IFLAG > 10)then
		WRITE(*,'(a,22I4)')'iErrFugi=',iFlag
		WRITE(*,*)'error in lower phase calculation'
		goto 81 
	endif
	etaLo=rhoLo*bMix
	rhoLo=rhoLo*rMwAvg ! Output as g/cc.
	DHL=hRes_RT  !these are passed through GlobConst
	DSL=sRes_R
	CvL=CvRes_R
	CpL=CpRes_R
	cmL=cmprsblty

	write(6,*)'Upper phase? Enter 0 for vapor, 1 for liquid.'
	read(*,*)iPhaseUp
	write(52,*)'Upper phase? Enter 0 for vapor, 1 for liquid.'
	write(52,*)iPhaseUp
	Y(1)=1
	WRITE(6,*)'ENTER Upper PHASE COMPOSITIONS    '
	if(NC > 1)READ(5,*)(Y(I),I=1,NC)
	WRITE(52,*)'Upper PHASE COMPOSITIONS    '
	write(52,*)(Y(I),I=1,NC)
	!CALL Fugi(T,P,Y,NC,iPhaseUp,fugcV,ZV,ier)
	CALL FugiTP( T,P,Y,NC,iPhaseUp,rhoUp,ZV,aRes,fugcV,uRes,iFlag )
	Hres_RT=uRes+ZV-1
	Sres_R=aRes-uRes
	bMix  =SUM( y(1:NC)*bVolCC_mol(1:NC) )
	rMwAvg=SUM( y(1:NC)*rMw(1:NC) )
	!wFrac(1:NC)=y(1:NC)*rMw(1:NC)/rMwAvg
	etaUp=rhoUp*bMix
	rhoUp=rhoUp*rMwAvg
	IF(iFlag==2 .and. LOUDER)write(*,*)'ERROR ALL UPPER PHASE'
	IF(iFlag==3 .and. LOUDER)write(*,*)'ERROR ALL LOWER PHASE'
	IF (IFLAG > 10 )then
		WRITE(*,'(a,22I4)')'iErrFugi=',iFlag
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
606	FORMAT(e12.5,3X,F7.2,1X,F7.4,1X,2(E11.4,1X),2(F11.8,1X))
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
	!WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND packing fraction   '
	!READ(5,*)T,eta
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND rhoMol_cc   '
	READ(5,*)T,rhoMol_cc

	IFLAG=0
    x(1)=1
    if(NC > 1)then
 		WRITE(6,*)'ENTER COMPOSITIONS '
		READ(5,*)(X(J),J=1,NC)
	endif
	WRITE(52,*)' COMPOSITIONS '
	write(52,*)(X(J),J=1,NC)
	bMix=SUM( x(1:Nc)*bVolCc_mol(1:Nc) )
	if(eta < zeroTol)then
		eta=rhoMol_cc*bMix
	else
		rhoMol_cc=eta/bMix
	endif
	vTotCc=bMix/eta
	WRITE(52,*)'TEMPERATURE(KELVIN),rhoMol_cc, packing fraction   '
	write(52,*)T,rhoMol_cc,eta
	rMwAvg=0
	do i=1,nc
		rMwAvg=rMwAvg+x(i)*rMw(i)
	enddo
	do i=1,nc
		wFrac(i)=x(i)*rMw(i)/rMwAvg
	enddo
	write(*,*)'Wt Fracs'
	write(*,'(11f8.5)')(wFrac(i),i=1,nc)
	if(initCall)print*,'FViter: calling Fugi for liq. T(K),eta=',T,eta
	isZiter=0
	Call FuVtot(isZiter,T,vTotCc,x,NC,FUGCL,ZL,aRes,uRes,iErr) !FUGC here has -lnZ already.
    !FUGCL(1:NC)=FUGCL(1:NC)-DLOG(ZL)
	Hres_RT=uRes+ZL-1
	Sres_R=uRes-aRes
	if(initCall)write(*,'(a,f8.2,E12.4)')' FViter: called Fugi for liq. T(K),ZL=',T,ZL
	etaLo=eta
	rhoLo=1/vTotCc
	PMPa = rhoLo*ZL*Rgas*T
	IF (iErr > 10 )then		! Ignore warnings (ie. iErr < 10)
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
    write(*,'(a,f8.2,f10.7,f10.5)')' T(K),rhoMol_cc,pMPa(calc)= ',T,rhoMol_cc,PMPa
	DO I=1,NC
		write(*,602)NAME(I),ID(I)
		write(6123,602)NAME(I),ID(I)
		activity=fugcL(I)-1
		if(activity > rLogInfinity)then			!polymers can cause this error.
			pause 'Fiter warning: exp overflow.  Setting kRatio to rLogInfinity.'
			activity=rLogInfinity
		endif
		VLK(I)=EXP(activity)
		write(*,603)
		write(6123,603)
		if(x(i) < 1/expLogInf)x(i)=1/expLogInf  !fix if x(i)=0 for one component.
		activity=LOG(x(i)*pMPa)+fugcL(i)	!log of the product is the sum of the logs.
		if(activity > rLogInfinity)then			!polymers can cause this error.
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
	if (iErrCode==0) then
	    outFile=TRIM(masterDir)//'\output\CriticalPoint.txt'
	    open(634,file=outFile)
		if (iEosOpt==5.or.iEosOpt==9) then
			write(634,600)
 			if(LOUD)write(*,*)'Tc=',TC_mix,' etaC=',eta
			if(LOUD)write(*,*)'Pc=',PC_mix,' Zc=',ZC_mix
			write(634,602)TC_mix,PC_mix,ZC_mix,eta
		else
 			write(634,601)
 			if(LOUD)write(*,*)'Tc=',TC_mix,' Pc=',PC_mix
			if(LOUD)write(*,*)'Vc=',VC_mix,' Zc=',ZC_mix
			write(634,602)TC_mix,PC_mix,ZC_mix,VC_mix
		endif
		close(634)
	 	if(LOUD)pause 'Results are saved as output/CriticalPoint.txt'
	endif

600	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, '    EtaC')
601	format(5x, 'Tc (K)', 5x, 'Pc (MPa)',4x, 'Zc', 6x, 'VC (cc/mol)')
602 format(2x,f8.3,5x,f6.4,5x,f5.4,5x,f12.5)

	RETURN
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PURPOSE:  CRITICAL POINT CALCULATION FOR PURE COMPOUNDS					   C
!C  PROGRAMMED BY:  AFG 2010												   C
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
	CALL CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,Acen_pure,TbCalc,iErrCode)	   
	if (iErrCode.eq.0) then
	    outFile=TRIM(masterDir)//'\output\CriticalPoint.txt'
	    open(635,file=outFile)
		if (iEosOpt==5.or.iEosOpt==9) then
			write(635,600)
            eta=bVolCC_mol(1)/VC_Pure
 			if(LOUD)write(*,*)'Tc=',TC_Pure,' etaC=',eta
			if(LOUD)write(*,*)'Pc=',PC_pure,' Zc=',ZC_Pure
			rhoc=PC_pure*rMw(1)/(TC_Pure*Rgas*ZC_Pure)
			if(LOUD)write(*,*)'rhoc(g/cc)=',rhoc
			write(635,602)TC_Pure,PC_PURE,ZC_PURE,eta,rhoc
		else ! must be ESD iEosOpt=2 or 4.
 			write(635,601)
 			!write(*,*)'Tc=',TC_Pure,' Pc=',PC_Pure
			!write(*,*)'Vc=',VC_Pure,' Zc=',ZC_Pure
            write(*,'()')
			write(635,602)TC_Pure,PC_PURE,ZC_PURE,VC_Pure
		endif
		close(635)
        print*,'DB values of Tc,Pc,Zc,acen'
        write(*,'(f8.2,4(1x,f8.4))')Tc(1),Pc(1),Zc(1),acen(1)
	 	print*, 'Results are saved as output/CriticalPoint.txt'
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
	SUBROUTINE CritDB(iErr)		!AUG 10
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE SpeadParms
	USE CritParmsDB
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*123 outfile !,inFile
	Character*77  errMsgPas
	!Integer idComp(nmx),idStore(1234)
	iErr=0
	isZiter=0
	tol=1.d-4
	NC=1
	outFile=TRIM(masterDir)//'\output\ParmsCritEos.txt'
	open(635,file=outFile)
	write(635,*)'    ID        TcEos        PcEos        ZcEos       AcenEos'
	write(635,*)'    ID        TcEos        PcEos        ZcEos       AcenEos'
	do iLine=2,nCritSet
		ID(1)=IDnum(iLine) ! Try every compound in ParmsCrit. If unavailable (ie. iErrStart > 0) in Get(EOS) then cycle.
		Call PGLWrapperStartup(NC,iEosOpt,idCasDb(iLine),iErrStart)
		write(*,611)'ID,Tc',ID(1),Tc(1)
		if(iErrStart > 0)then
			write(*,*)'iLine,iErrStart',iLine,iErrStart,' ',TRIM(errMsgPas)
			cycle
		endif
		CALL CritPure(NC,isZiter,tol,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,Acen_pure,TbCalc,iErrCode)
		if(iErrCode > 0)cycle
		write( * ,601)ID(1),TC_Pure,PC_Pure,ZC_Pure,Acen_pure	   
		write(635,601)ID(1),TC_Pure,PC_Pure,ZC_Pure,Acen_pure	   
	enddo !iLine=1,nDeck
601 format(1x,i6,1x,12E13.5)
611 format(1x,a,i6,1x,12E13.5)
	close(635)
	write(*,*)'CritDB: Done. Output in ',TRIM(outFile)	 	
    RETURN
	END


	SUBROUTINE FREEZINGCURVECALC(NC)
	!
	!  PURPOSE:  iterative solubility calculation using activity coefficients for crystalline solids

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*64 Input_File_Name
	character outFile*251
	DIMENSION fugcL(NMX),X(NMX),ier(20)
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
	DIMENSION X(NMX),fugcL(NMX),ier(20)
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
	USE FugiParts
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER answer*1
	character outFile*251
	DIMENSION X(NMX),gmol(NMX),x1(19)	!,ier(20)
	DIMENSION fugcPure(NMX),fugRepPure(nmx),fugAttPure(nmx),fugAssocPure(nmx),hPure(NMX),sPure(NMX),fugc(NMX)
	DoublePrecision rLnGam(NMX)
	LOGICAL LOUDER
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
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
	open(642,file=outFile)
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
	write(642,*)'  T(K)   Eta'
	write(642,'(f8.2,f8.3)')tKelvin,eta
    isZiter=0
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
		call FuVtot(isZiter,tStore,vTotCc,x,NC,FUGC,zFactor,aRes,uRes,iErr)
		Hres_RT=uRes+zFactor-1
		Sres_R=aRes-uRes
		if(iErr /= 0)then
			print*,' GibbsXsVsEta: Error returned from pure FuVtot. iPure,Z=',iPure,zFactor
			print*,' Try a larger eta?'
			goto 1000
		endif

		fugcPure(iPure)=FUGC(iPure)-Zfactor-LOG(vTotCc)        ! For the total part, we should also subtract lnZ, but that can give ln(-ve). We really want ln(Vi/V), so we initiate with ln(Vi).. 
		fugRepPure(iPure)=FUGRep(iPure)-Zrep       ! Vk/V = Vk/Vk for pure.	   
		fugAttPure(iPure)=FUGAtt(iPure)-Zatt
		fugAssocPure(iPure)=FUGAssoc(iPure)-Zassoc
		hPure(iPure)=Hres_RT
		call FuVtot(isZiter,tInf,vTotCc,x,NC,FUGC,zFactor,aRes,uRes,iErr)
		Hres_RT=uRes+zFactor-1
		Sres_R=aRes-uRes
		sPure(iPure)=aRes-uRes
	enddo
	!if(nc.eq.2 .and. LOUDER)write(6 ,*)'   x1      Z         uXs      -s0Xs        aXs      gIs         V(cc)','     LnGam'
	if(nc.eq.3 .and. LOUDER)write(6 ,*)'   x1      x2      Z         hE      gMix         gE       gIs        V(cc)       LnGam(i)'
	if(nc.eq.2)write(642,*)'   x1      x2      uXs      -s0Xs        aXs      gIs        Z         V(cc)       LnGam(i)'
	if(nc.eq.3)write(642,*)'   x1      x2      x3      hE      gMix         gE       gIs        Z         V(cc)       LnGam(i)'
	write(* ,*)' xi,PMPa,vTotCc,F,lnGam(1,2),lnGamRep(1,2),lnGamAtt(1,2)',',lnGamAssoc(1,2)'
	nX3=nX
	if(nc==2)nX3=1
	do iX3=1,nX3
		do iX1=1,nx
			x(3)=x1(iX3)
			if(nc==2)x(3)=0
			x(1)=x1(iX1)*(1-x(3))
			x(2)=1-(x(1)+x(3))
			sumx=0
			do i=1,NC
                if(LOUDER)then
				    if(x(i).le.0.or.x(i).ge.1)pause 'GibbsXs error in xCalc'
                end if
				sumx=sumx+x(i)
            enddo
            if(LOUDER)then
			    if(ABS(sumx-1) > zeroTol)pause 'GibbsXs error in xSum'
            end if

			gmol(1:nc)=x(1:nc)
			bVolMix=sumproduct(nc,x,bVolCC_mol)
			vTotCc=bVolMix/eta
			!call FuTptVtot(zFactor,aDep,uDep,eta,tKelvin,x,bVolMix,nComps,iErr)
			call FuVtot(isZiter,tKelvin,vTotCc,x,NC,FUGC,zFactor,aRes,uRes,iErr)
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
			hE=uRes-hIs
			aE=aRes-aIs	
			call FuVtot(isZiter,tKelvin,vTotCc,x,NC,FUGC,zFactor,aRes,uRes,iErr)
			!call FuTptVtot(isZiter,zFactor,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErr)
			sDep=aRes-DLOG(zFactor)-uRes
			sE=sDep-sIs !really this is -ve sXs/R
			!if(LOUDER)write(* ,611)(x(i),i=1,nc-1),zFactor,hE,sE,aE,gIs,vTotCc,(rLnGam(i),i=1,nc)
			PMPa=Zfactor*Rgas*tKelvin/vtotcc
			write(* ,612)x(1),PMPa, vTotCc, (rLnGam(i),i=1,nc),(rLnGamRep(i),i=1,nc),(rLnGamAtt(i),i=1,nc),(rLnGamAssoc(i),i=1,nc)
			write(642,611)(x(i),i=1,nc),hE,sE,aE,gIs,zFactor,vTotCc,(rLnGam(i),i=1,nc)
		enddo
	enddo

	CLOSE(642)
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
SUBROUTINE IDACsTP(tKelvin,pMPa,NC,rLnGam,iErr,errMsgPas)
	!C
	!C  PURPOSE:  Gibbs Excess vs. composition at const T,P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE FugiParts
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	!character outFile*251
	character*77 errMsg(0:22),errMsgPas
	!DoublePrecision fugcL(NMX) !X(NMX),x1(0:18),
	!DoublePrecision hPure(NMX),gPure(NMX),Vcc_mol(NMX),fugRepPure(NMX),fugAttPure(NMX),fugAssocPure(NMX)
	DoublePrecision rLnGam(NMX),xFrac(NMX),fugcPure(NMX),fugcLo(NMX)
	LOGICAL LOUDER
	LOUDER=LOUD
	!LOUDER = .TRUE.
	iErr=0
	errMsg( 0)='IDACsTP: No problem!'
	errMsg( 1)='IDACsTP, Warning: tKelvin < 100 K? '
	errMsg(11)='IDACsTP: Sorry, IDACsTP is only for NC=2'
	errMsg(12)='IDACsTP: Sorry, ierFugi > 10.'
    if(NC /= 2)then
	    write(dumpUnit,*) 'error - this option is only for 2 or 3 compos'
		iErr=11
		goto 861
    end if

	!c	write(*,*)'enter file NAME for output'
	!c	read(*,'(A1)')outfile
	!outFile=TRIM(masterDir)//'\output\IDACs.txt'
	!open(643,file=outFile)
	if(tKelvin < 100)then
		write(dumpUnit,*)'IDACsTP: warning T<100K. Sure?(y/n)'
		iErr=1
		goto 861
	endif
	xFrac(1:NC)=zeroTol
	iErr12=0
	do iComp=1,NC
		xFrac(iComp)=1
		Call FugiTP( tKelvin,pMPa,xFrac,NC,1,rhoMol_cc,zLo,aRes,fugcLo,uRes,iErrFugi )
		if(iErrFugi > 10)iErr12=1
		fugcPure(iComp)=fugcLo(iComp)
		if(iComp==1)rLnGam(2)=fugcLo(2)
		if(iComp==2)rLnGam(1)=fugcLo(1)
		xFrac(iComp)=zeroTol
	enddo
	if(iErr12 > 0)iErr=12
	rLnGam(1:NC)=( rLnGam(1:NC)-fugcPure(1:NC) )
	!write(*,*)'OutFile=',TRIM(outfile)
	!pause 'Check results in OutFile.'
861	continue
	errMsgPas=TRIM( errMsg(iErr) )   
	!CLOSE(643)
	RETURN
END !Subroutine IDACsTP()
!***********************************************************************
	SUBROUTINE GibbsXs(NC)
	!C
	!C  PURPOSE:  Gibbs Excess vs. composition at const T,P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE FugiParts
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	CHARACTER answer*1
	character outFile*251
	DoublePrecision X(NMX),x1(0:18),fugcL(NMX)
	DoublePrecision hPure(NMX),gPure(NMX),fugRepPure(NMX),fugAttPure(NMX),fugAssocPure(NMX),fugcPure(NMX)
	DoublePrecision rLnGam(NMX),Vcc_mol(NMX)
	LOGICAL LOUDER
	data nx,x1/18,0.0001,.001,.01,.02,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,.98,.99,.999,0.9999/
	LOUDER=LOUD
	!LOUDER = .TRUE.
    if(NC>3)then
	    write(dumpUnit,*) 'error - this option is only for 2 or 3 compos'
		return
    end if

	!c	write(*,*)'enter file NAME for output'
	!c	read(*,'(A1)')outfile
	outFile=TRIM(masterDir)//'\output\xsGibbs.txt'
	open(643,file=outFile)
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND PRESSURE(MPa), xStep (0 for default)   '
	READ(5,*)T,P,xStep
	if(T < 100)then
		write(*,*)'warning T<100K. Sure?(y/n)'
		read(*,*)answer
		if(answer.eq.'n'.or.answer.eq.'N')goto 1000
	endif
	zero=0
	if(LOUDER)write(* ,'(a,f7.2,a,f7.3)')' T(K)=',T,'P(MPa)=',P
	write(643,*)'  T(K)   P(MPa)'
	write(643,'(f8.2,f8.3)')T,P
	x(1:NC)=zeroTol
	do iPure=1,NC
		x(iPure)=1-zeroTol*(NC-1)
		!call Fugi(t,p,x,NC,1,fugcL,zL,ier)
		call FugiTP( T,P,x,NC,1,rhoMol_cc,zL,aRes,fugcL,uRes,iErrF )
		!bMix=sumproduct(NC,x,bVolCc_mol) ! no need to calculate this. bMix=bVol(iPure)
		!if(LOUDER)write(dumpUnit,form610)' GibbsXs: GpureAssoc=',fugAssocPure(1:NC)
		Vcc_mol(iPure)=1/rhoMol_cc
		fugRepPure(iPure)=fugRep(iPure)-zRep
		fugAttPure(iPure)=fugAtt(iPure)-zAtt
		fugAssocPure(iPure)=fugAssoc(iPure)-zAssoc
		fugcPure(iPure)=fugcL(iPure)
		hPure(iPure)=uRes+zL-1
		if(zL < zeroTol)write(dumpUnit,*)'GibbsXs: 0=zL=',zL
		gPure(iPure)=aRes+zL-1-LOG(zL)
		if(LOUDER)write(*,form611)' i,V,H,G,fAss,uAss=',iPure,Vcc_mol(iPure),hPure(iPure),gPure(iPure),fugAssocPure(iPure),uAssoc
		x(iPure)=zeroTol
	enddo
	if(NC==2 .and. LOUDER)write(6 ,'(a)')'   x1      x2       hE     gMixGE        gE       gIs        LnGam(i)'
	if(NC==3 .and. LOUDER)write(6 ,'(a)')'   x1      x2      x3      hE     gMixGE        gE       gIs        LnGam(i)'
	if(NC==2)write(643,'(a)')'   x1      x2       hE      gMixGE        gE       gIs        LnGam(i)           gMix'
	if(NC==3)write(643,'(a)')'   x1      x2      x3      hE     gMixGE        gE       gIs        LnGam(i)'

	if(bTPT.or.bESD.or.iEosOpt==19)write(* ,'(a,a)')'   x1     eta     hEkJ   vE(cc)    lnGam(1,2)   lnGamRep(1,2)    lnGamAtt(1,2)','  lnGamAssoc(1,2)    '
	nXstep=nX
	if(xStep > zeroTol)nXstep=1+1/xStep
	nX3=nXstep
	xMin=1.D-4
	if(NC==2)nX3=0
	do iX3=0,nX3,1
		do iX1=0,nXstep,1
			x(3)=x1(iX3)
			if(NC > 2 .and. x(3) <   xMin)x(3)=  xMin
			if(NC > 2 .and. x(3) > 1-xMin)x(3)=1-xMin
			if(NC.eq.2)x(3)=0
			if(xStep < zeroTol)then
				x(1)=x1(iX1)*(1-x(3))
				x(2)=1-(x(1)+x(3))
			else
				x(1)=iX1*xStep*(1-x(3))
				if(x(1) < xMin  ) x(1)=  xMin
				if(x(1) > 1-xMin) x(1)=1-xMin
				x(2)=1-(x(1)+x(3))
			endif
			sumx=0
			do i=1,NC
                if(LOUDER)then
				    if(x(i).le.0.or.x(i).ge.1)pause 'GibbsXs error in xCalc'
                end if
				sumx=sumx+x(i)
            enddo
            if(LOUDER)then
			    if(ABS(sumx-1) > zeroTol)pause 'GibbsXs error in xSum'
            end if
			!call Fugi(t,p,x,NC,1,fugcL,zL,ier)
			call FugiTP( T,P,x,NC,1,rhoMol_cc,zL,aRes,fugcL,uRes,iErrF )
			gMix=aRes+zL-1-LOG(zL)
			vTotCc=ZL*Rgas*T/P
			gIs=0
			bMix  =SumProduct( NC,x,bVolCC_mol )
			eta=bMix/vTotCc
			do iComp=1,NC
				gIs=gIs+x(iComp)*LOG(x(iComp))
				if(bESD.or.bTPT.or.iEosOpt==19)then
					rLnGamRep(iComp)=fugRep(iComp)-Zrep*bVolCC_mol(iComp)/bMix-fugRepPure(iComp)
					rLnGamAtt(iComp)=fugAtt(iComp)-Zatt*bVolCC_mol(iComp)/bMix-fugAttPure(iComp)
					rLnGamAssoc(iComp)=fugAssoc(iComp)-Zassoc*bVolCC_mol(iComp)/bMix-fugAssocPure(iComp)
				endif
				rLnGam(iComp)=fugcL(iComp)-fugcPure(iComp) 
				!rLnGam(iComp)=fugc(iComp)-Zfactor*bVolCc_mol(iComp)/bVolMix - LOG(vTotCc)-fugcPure(iComp)
			enddo
			hIs  =SUM( x(1:NC)*hPure(1:NC) ) !SumProduct(NC,x,hPure )
			vIs  =SUM( x(1:NC)*Vcc_mol(1:NC) ) !SumProduct(NC,x,hPure )
			gEgam=SumProduct(NC,x,rLnGam )
			gEmix=gMix - gIs - SumProduct(NC,x,gPure)
			gMix=gEmix+gIs
			hE=uRes+zL-1-hIs
			hEkJ_mol=hE*Rgas*T/1000
			vE=1/rhoMol_cc-vIs
			!if(LOUDER)write(* ,611)(x(i),i=1,NC),hE,gMix,gE,gIs,(rLnGam(i),i=1,NC)
			write(643,611)(x(i),i=1,NC),hE,gMix,gEmix,gIs,(rLnGam(i),i=1,NC),gEgam
			if(bESD .or. bTPT .or. iEosOpt==19)then
				write(dumpUnit,612)x(1),eta, hEkJ_mol, vE, (rLnGam(i),i=1,NC),(rLnGamRep(i),i=1,NC),(rLnGamAtt(i),i=1,NC),(rLnGamAssoc(i),i=1,NC)
			else
				write(dumpUnit,611)x(1:NC),hE,gMix,gEmix,gIs,rLnGam(1:NC)
				write(643,611)x(1:NC),hE,gMix,gEmix,gIs,rLnGam(1:NC)
			endif
		enddo  !iX1
	enddo !iX3

	CLOSE(643)
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
611	format(3(f7.4,1x),5(f10.3,1x),f12.9)    
612	format(2(f7.4,1x),2(f8.3),12(f8.3))    
	END

!***********************************************************************
	SUBROUTINE GibbsXsVtot(NC)
	!C
	!C  PURPOSE:  Gibbs Excess vs. composition at const T,Vtot.	Handy for checking fugacity expressions.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,
	USE FugiParts ! fugRep(nmx),fugAtt(nmx),fugAssoc(nmx),Zrep,Zatt,Zassoc,aRep,aAtt,aAssoc,uAtt,uAssoc
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	CHARACTER answer*1
	character outFile*251
	DoublePrecision X(NMX),x1(0:18),fugcL(NMX)
	DoublePrecision hPure(NMX),gPure(NMX),fugRepPure(NMX),fugAttPure(NMX),fugAssocPure(NMX),fugcPure(NMX)
	DoublePrecision rLnGam(NMX),xTemp(nmx)	!,Vcc_mol(NMX)
	LOGICAL LOUDER
	data nx,x1/18,0.0001,.001,.01,.02,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,.98,.99,.999,0.9999/
	LOUDER=LOUD
	LOUDER = .TRUE.
    if(NC>3)then
	    write(dumpUnit,*) 'error - this option is only for 2 or 3 compos'
		return
    end if

	!c	write(*,*)'enter file NAME for output'
	!c	read(*,'(A1)')outfile
	outFile=TRIM(masterDir)//'\output\xsGibbs.txt'
	open(643,file=outFile)
1000  CONTINUE
	WRITE(6,*)'ENTER TEMPERATURE(KELVIN) AND vTotCc, xStep (0 for default)   '
	READ(5,*)T,vTotCc,xStep
	if(T < 100)then
		write(*,*)'warning T<100K. Sure?(y/n)'
		read(*,*)answer
		if(answer.eq.'n'.or.answer.eq.'N')goto 1000
	endif
	zero=0
	if(LOUDER)write(* ,'(a,f7.2,a,f7.3)')' T(K)=',T,'P(MPa)=',P
	write(643,*)'  T(K)   Vtot(cm3)'
	write(643,'(f8.2,f8.3)')T,vTotCc
	isZiter=0
	xInf=zeroTol
	do iPure=1,NC
		sum=0
		do iComp=1,NC
			x(iComp)=xInf
			sum=sum+x(iComp)
		enddo
		x(iPure)=1-sum
		!call Fugi(t,p,x,NC,1,fugcL,zL,ier)
		Call FuVtot(isZiter,T,vTotCc,x,NC,FUGCL,ZL,aRes,uRes,iErrF)
		P=zL*Rgas*T/vTotCc
		!bMix=sumproduct(NC,x,bVolCc_mol) ! no need to calculate this. bMix=bVol(iPure)
		fugRepPure(iPure)=fugRep(iPure)-zRep
		fugAttPure(iPure)=fugAtt(iPure)-zAtt
		fugAssocPure(iPure)=fugAssoc(iPure)-zAssoc
		fugcPure(iPure)=fugcL(iPure)
		hPure(iPure)=uRes+zL-1
		gPure(iPure)=aRes+zL-1
	enddo
	if(NC==2 .and. LOUDER)write(6 ,'(a)')'   x1      x2       hE     gMixGE        gE       gIs        LnGam(i)'
	if(NC==3 .and. LOUDER)write(6 ,'(a)')'   x1      x2      x3      hE     gMixGE        gE       gIs        LnGam(i)'
	if(NC==2)write(643,'(a)')'   x1      x2       hE      gMixGE        gE       gIs        LnGam(i)           gMix'
	if(NC==3)write(643,'(a)')'   x1      x2      x3      hE     gMixGE        gE       gIs        LnGam(i)'

	if(bTPT.or.bESD)write(* ,'(a,a)')'     x1        aRes      zFactor     fugc(1,2)     fugRep(1,2)      fugAtt(1,2)','    fugAssoc(1,2)    '

	nXstep=nX
	if(xStep > zeroTol)nXstep=1+1/xStep
	nX3=nXstep
	xMin=1.1D-4
	xDiff=1.D-4
	if(NC.eq.2)nX3=0
	do iX3=0,nX3,1
		do iX1=0,nXstep,1
			x(3)=x1(iX3)
			if(NC > 2 .and. x(3) <   xMin)x(3)=  xMin
			if(NC > 2 .and. x(3) > 1-xMin)x(3)=1-xMin
			if(NC.eq.2)x(3)=0
			if(xStep < zeroTol)then
				x(1)=x1(iX1)*(1-x(3))
				x(2)=1-(x(1)+x(3))
			else
				x(1)=iX1*xStep*(1-x(3))
				if(x(1) < xMin  ) x(1)=  xMin
				if(x(1) > 1-xMin) x(1)=1-xMin
				x(2)=1-(x(1)+x(3))
			endif
			sumx=0
			do i=1,NC
                if(LOUDER)then
				    if(x(i).le.0.or.x(i).ge.1)pause 'GibbsXs error in xCalc'
                end if
				sumx=sumx+x(i)
            enddo
            if(LOUDER)then
			    if(ABS(sumx-1) > zeroTol)pause 'GibbsXs error in xSum'
            end if
			xTemp(1)=x(1)+xDiff
			xTemp(2)=x(2) 		 !hold n2=const.
			Call FuVtot(isZiter,T,vTotCc,xTemp,NC,FUGCL,ZL,aResPlus,uRes,iErrF)
			aPlus=aRep
			xTemp(1)=x(1)-xDiff
			xTemp(2)=x(2) 		 !hold n2=const.
			Call FuVtot(isZiter,T,vTotCc,xTemp,NC,FUGCL,ZL,aResMinus,uRes,iErrF)
			aMinus=aRep
			fugcNum=(aPlus*(1+xDiff)-aMinus*(1-xDiff))/(2*xDiff)
			Call FuVtot(isZiter,T,vTotCc,x,NC,FUGCL,ZL,aRes,uRes,iErrF)
			P=zL*Rgas*T/vTotCc
			gMix=aRes+zL-1 !-LOG(zL)
			gIs=0
			bMix  =SumProduct( NC,x,bVolCC_mol )
			eta=bMix/vTotCc
			do iComp=1,NC
				gIs=gIs+x(iComp)*LOG(x(iComp))
				rLnGamRep(iComp)=fugRep(iComp)-Zrep*bVolCC_mol(iComp)/bMix-fugRepPure(iComp)
				rLnGamAtt(iComp)=fugAtt(iComp)-Zatt*bVolCC_mol(iComp)/bMix-fugAttPure(iComp)
				rLnGamAssoc(iComp)=fugAssoc(iComp)-Zassoc*bVolCC_mol(iComp)/bMix-fugAssocPure(iComp)
				rLnGam(iComp)=fugcL(iComp)-fugcPure(iComp) 
				!rLnGam(iComp)=fugc(iComp)-Zfactor*bVolCc_mol(iComp)/bVolMix - LOG(vTotCc)-fugcPure(iComp)
			enddo
			hIs  =SumProduct(NC,x,hPure )
			gEgam=SumProduct(NC,x,rLnGam )
			gEmix=gMix - gIs - SumProduct(NC,x,gPure)
			gMix=gEmix+gIs
			hE=uRes+zL-1-hIs
			!if(LOUDER)write(* ,611)(x(i),i=1,NC),hE,gMix,gE,gIs,(rLnGam(i),i=1,NC)
			write(643,611)(x(i),i=1,NC),hE,gMix,gEmix,gIs,(rLnGam(i),i=1,NC),gEgam
			if(bESD .or. bTPT)then
				write(dumpUnit,612)x(1),aRes, ZL, (fugcL(i),i=1,NC),(fugRep(i),i=1,NC),(fugAtt(i),i=1,NC),(fugAssoc(i),i=1,NC),fugcNum
				!write(dumpUnit,612)x(1),eta, P, (rLnGam(i),i=1,NC),(rLnGamRep(i),i=1,NC),(rLnGamAtt(i),i=1,NC),(rLnGamAssoc(i),i=1,NC)
			else
				write(643,611)x(1:NC),hE,gMix,gEmix,gIs,rLnGam(1:NC)
			endif
		enddo  !iX1
	enddo !iX3

	CLOSE(643)
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
611	format(3(f7.3,1x),5(f10.3,1x),f12.9)    
612	format(3(f10.6,1x),12(f9.4))    
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
    LOGICAL LOUDER
	!COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data tList/0.01,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.1/
    LOUDER=LOUD
    LOUDER=.TRUE.
	iErrCode=0
	errMsg(0)='No Problem in Isochore'
	errMsg(1)='Isochore Error: NC.ne.1'
	errMsg(2)='Isochore Error: itMax exceeded.'
	errMsg(3)='Isochore Error: NC.ne.1'
	if(NC.ne.1)then
		iErrCode=1
		if(LOUDER)write(*,*)'Isochore Error: NC=1 only. Restart CalcEos with number of components = 1'
		goto 86
	endif
	gmol(1)=1				     !Arbitrarily setting 1 mole for basis.
	outFile=TRIM(masterDir)//'\output\Isochore.txt'
	open(644,file=outFile)
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
	if(eta > etaMax)print*,'Isochore: etaMax < eta=',eta
	if(eta > etaMax)goto 10
	write(* ,*)'rho(g/cc) =',rhoLiqG_cc,' V(cc/mol) = ',vTotCc
	write(644,*)'rho(g/cc) =',rhoLiqG_cc,' V(cc/mol) = ',vTotCc
	!if(LOUDER)write(* ,*)'  T(K),      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT,    CvRes/R,     (dP/dRho)/RT'
	!write(644,*)'  T(K),      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT,    CvRes/R,     (dP/dRho)/RT'
	if(LOUDER)write(*,*)'  T(K)    PMPa      Z     Ures/RTc   Ares/RT  (dP/dRho)/RT  CvRes/R  CpRes/R' !   Sdep'    
	write(644,'(a)')'  T(K)    PMPa      Z     Ures/RTc   Ares/RT  (dP/dRho)/RT  CvRes/R  CpRes/R    HRes/RT'    
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
		isZiter=0  ! ln(-ve Z) should not happen for FuVtot.
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,aRes,uRes,iErrFu)
		if(iErrFu)exit
		pNew = zFactor*8.314*tKelvin/vTotCc
		if(bTPT)call NUMDERVS(NC,gmol,tKelvin,1/vTotCc,zFactor,iErrNum)	!Derv props USEd in GlobConst
		if(LOUDER)write(* ,602)tKelvin,pNew,zFactor,uRes*tKelvin/TC(1),aRes,cmprsblty,cvRes_R,CpRes_R
		write(644,602)tKelvin,pNew,zFactor,uRes*tKelvin/TC(1),aRes
	enddo !new temperatures
	print*,'Enter Tstart,Tstop,increment'
	read*,tStart,tStop,tInc
	if(LOUDER)write(*,*)'  T(K)    PMPa      Z     Ures/RTc   Ares/RT  (dP/dRho)/RT  CvRes/R  CpRes/R' !   Sdep'    
	write(644,'(a)')'  T(K)    PMPa      Z     Ures/RTc   Ares/RT  (dP/dRho)/RT  CvRes/R  CpRes/R    HRes/RT'    
	tKelvin=tStart-tInc
	do while(tKelvin <= tStop)
		tKelvin=tKelvin+tInc
		isZiter=0 !=> derivative props moved separate from ln(Z)
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,zFactor,aRes,uRes,iErrFu)
		if(iErrFu)exit
		if(bTPT)call NUMDERVS(NC,gmol,tKelvin,1/vTotCc,zFactor,iErrNum)	!
		pNew = zFactor*rGas*tKelvin/vTotCc
		if(LOUDER)write(* ,602)tKelvin,pNew,zFactor,uRes*tKelvin/TC(1),aRes,cmprsblty,cvRes_R,CpRes_R
		write(644,602)tKelvin,pNew,zFactor,uRes*tKelvin/TC(1),aRes,cmprsblty,cvRes_R,CpRes_R,uRes+zFactor-1  !cmprsblty=(dP/dRho)T*(1/RT)
	enddo
601	format(f11.4,8(' ',f11.4),3(1x,E11.4))
602	FORMAT(2F9.3,1x,1(F7.4),8(f11.4,1X))      
86	continue
	if(LOUDER)pause 'Your results are stored in Isochore.txt'
	close(644)
	!close(666)
	errMsgPas=errMsg(iErrCode)
	if(iErrCode.ne.0)then
		if(LOUDER)write(*,*)errMsgPas
		if(LOUDER)pause
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
	DIMENSION gMol(NMX),fugc(NMX),tList(nList)
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
	open(645,file=outFile)
	!open(666,file='debugData.txt')
	!WRITE(6,*)'ENTER packing fraction   '
	!READ(5,*)eta
	!vTotCc = bVolCc_mol(1)/eta	 !Arbitrarily setting 1 mole for basis.
	!write(6 ,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
	!write(61,*)'Packing fraction =',eta,' V(cc/mol) = ',vTotCc
10	WRITE(6,*)'ENTER Pressure (MPa)   '
	READ(5,*)pMPa
	if(LOUDER)write(* ,*)'  T(K),      rhog/cc,     zFactor,   (A-Aig)/RT,  (U-Uig)/RTc,    eta' !   CvRes/R,     (dP/dRho)/RT'
	write(645,'(a155)')'  T(K),     rhog/cc,     zFactor,   (A-Aig)/RT,  (U-Uig)/RTc,      eta'	     
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
		CALL FugiTP( tKelvin,pMPa,gmol,NC,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErrF )
		if( iErrF > 10 )exit
		rhoNew = rhoMol_cc*rMw(1)
		if(LOUDER)write(* ,601)tKelvin,rhoNew,zFactor,aRes,uRes*tKelvin/TC(1),etaPass
		write(645,601)tKelvin,rhoNew,zFactor,aRes,etaPass,uRes*tKelvin/TC(1),etaPass
	enddo !new temperatures
	if(LOUDER)write(*,*)'  T(K)    rhog/cc    Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes_R/R' !   Sdep'    
	write(645,'(a)')'  T(K)    rhog/cc    Z     Ares/RT   Ures/RT  (dP/dRho)/RT  CvRes/R  CpRes/R    HRes/RT     eta'    
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
		!call FUGI(tKelvin,pMPa,gmol,NC,LIQ,FUGC,zFactor,ierFugi)
		CALL FugiTP( tKelvin,pMPa,gmol,NC,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErrF )
		if(iErrF>10)then
			if(LOUDER)print*,'Isobar: T,ierFugi=',tKelvin,iErrF
			exit
		endif
		if( ABS(zFactor) > 1D-11) rhoMol_cc=pMPa/(zFactor*rGas*tKelvin)
		rhoG_cc = rhoMol_cc*rMw(1)
		eta=rhoMol_cc*bVolCC_mol(1)
		call NUMDERVS(NC,gMol,tKelvin,rhoMol_cc,zFactor,iErrNum)	!It's just way easier and more reliable to get derivative properties numerically.
		if(bTPT.and.initCall)print*,'Isotherm: CvRes_R,pas=',CvRes_R
		if(LOUDER)write(* ,602)tKelvin,rhoNew,zFactor,aRes,uRes,cmprsblty,cvRes_R,CpRes_R,eta
		write(645,602)tKelvin,rhoG_cc,zFactor,aRes,uRes,cmprsblty,cvRes_R,CpRes_R,uRes+zFactor-1,eta  !cmprsblty=(dP/dRho)T*(1/RT)
		initCall=0
	enddo
601	format(f11.4,6(',',f11.4))
602	FORMAT(F9.3,f9.6,1x,1(F7.4),9(F9.4,1X))      
86	continue
	close(645)
	!close(666)
	IF(LOUDER)write(dumpUnit,*)'Output file:',TRIM(outFile)
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
	open(651,file=outFile)
10	write(6,*)'ENTER TEMPERATURE (<0 to stop)'
	read(5,*)tKelvin
	if(tKelvin < 0)goto 86
	tr=tKelvin/TC(1)
	!tKelvin=TC(1)*tr
	write(*  ,*)'Temperature =',tKelvin,' Tc = ',TC(1)
	write(651,*)'Temperature =',tKelvin,' Tc = ',TC(1)
	etaList(1)=1D-11
	dEta=etaMax/100
	do I=2,nEta	!nEta is a parameter
		etaList(I)=(I-1)*dEta
	enddo
	gmol(1)=1
	isZiter=1
	write(* ,600)	!!!!!!!!!!!!!!!!!!!!!!!! Write Header !!!!!!!!!!!!!!!!1
	write(651,600)
	if(bVolCC_mol(1) < 1 .and. LOUD)print*,'Isotherm: 0~bVol=',bVolCc_mol(1)
	if(initCall)print*,'Isotherm: bVol=',bVolCc_mol(1)
	do i=1,nEta   
		if (I==1) write(651,600)
		eta=etaList(i)
		rhoMol_cc=eta/bVolCC_mol(1)
			if(initCall.and.Loud .and. i < 3)print*,'Isotherm: Calling FuVtot. eta,rhoMol_cc=',eta,rhoMol_cc
		call FuVtot(isZiter,tKelvin,1/rhoMol_cc,gmol,NC,FUGC,Z,aRes,uRes,iErrFu)
		if(bTPT)call NUMDERVS(NC,gMol,tKelvin,rhoMol_cc,Z,iErrNum)	!
!		if(bTPT)print*,'Isotherm: cmprsblty,pas=',cmprsblty,cmprsbltyPas
!		if(bTPT)cmprsblty=cmprsbltyPas
		pMPa=Z*Rgas*tKelvin*rhoMol_cc
		write(* ,601) eta,PMpa,Z,uRes,aRes,rhoMol_cc,cmprsblty !cmprsblty from GlobConst
		write(651,601) eta,PMpa,Z,uRes,aRes,rhoMol_cc,cmprsblty
	enddo
	goto 10
86	continue
	close(651)
	errMsgPas=errMsg(ier)
	if(ier.ne.0)then
		write(*,*)errMsgPas
		pause
	endif
	write(*,*) 'Your results are stored in ',TRIM(outFile)
	pause 'Isotherm done!'
600 format(4x,'eta',7x,'PMPa',9x,'Z',9X,'uRes_RT',8x,'aRes_RT',8x,'rho(mol/cc)',9x,'cmprsblty')			
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
	open(652,file=outFile)
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
			CALL FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
			B2grid(I)=(Z-1.d0)/eta
		enddo
		if (J.EQ.1) WRITE(6 ,600)		
		if (J.EQ.1) WRITE(652,600)		
		do k=1,(nEta-1)
			deltaEta(k)=etaList(k+1)-etaList(k)
			B3grid(k)=(B2grid(k+1)-B2grid(k))/deltaEta(k)
		enddo
	!Since Virial Coefficient are defined as eta goes to zero:
		B2=B2grid(1) *bVolCC_mol(1)
		B3=B3grid(1) *bVolCC_mol(1)**2
		WRITE(6  ,601)tKelvin,betaRed,B2,B3	
		WRITE(652,601)tKelvin,betaRed,B2,B3	
	enddo
	close(652)
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
	open(653,file=outFile)
	WRITE(6,*)'ENTER temperature   '
	READ(5,*)trInv
	tKelvin=TC(1)/trInv
	write(653,*)'Temperature =',tKelvin,' Tc = ',TC(1)
	!Note: This is wacky since it is easy to compute Z,U given rho,
	!but I don't have a standard interface for ZCalc of all eosOpt's.
	!I only have a standard interface for fugi.
	xFrac(1)=1
	itMax=111
	LIQ=1
	tol=0.001
	pOld=255	!initialize here so we bootstrap
	write(653,*)'  1/Tr,      pMpa,     zFactor,  (U-Uig)/RTc,   (A-Aig)/RT'
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
		write(653,601)trInv,pNew,zFactor,dUoNkT*tKelvin/TC(1),dAoNkT
	enddo !new temperatures
601	format(f7.4,5(',',f11.4))
86	continue
	if(LOUD)pause 'Your results are stored in Isotherm.txt'
	close(653)
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

	SUBROUTINE HEDB(NC)
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
	character*77 errMsgSub
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	character(len=123) :: dumString,refid !errMsgLookup !errMsg(11),
	!character*11 dum1,dum2,dum3,dum4
	!DoublePrecision deviate(maxPts),parm(5),stdErr(5),temParm(5)
	Integer idStore(NMX) !,idcc(2),QueryNParMix !,iErrGet,idCas(NMX)
	DoublePrecision gPure(NMX),hPure(NMX),Vcc_mol(NMX),x(NMX),TR(NMX),fugcL(NMX)
	LOGICAL LOUDER
	!EXTERNAL LMDevFcn
	character*77 inFile,outFile,errMsgPas !,dumString ,name1,name2
	!character*1 tabChar !,dumString(77),refID(77)
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/FloryWert/ vLiq(nmx)
    !data initial/1/
	!tabChar = char(9)
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
	iAns=1
	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in VleData subdir)'
	WRITE(*,*)'FIRST LINE =NPTS,id1,id2 2ND LINE = TK, PMPA, X,Y, ...'
    FN="\HEMEM2b.txt"
	if(iAns==2)FN="\VLEJaubertAq.txt"
	if(iAns==3)FN="\VLEDannerGess.txt"															
	if(iAns==4)FN="\SolvationSystems.txt"
	if(iAns < 1)then
		write(*,*)'Enter FN.FT (path is assumed to be PGLInputDir.'
		read(*,*)FN
		FN='\'//TRIM(FN)
	endif

!	READ(*,'(A)')FN
	inFile=TRIM(PGLinputDir)//TRIM(FN)
	if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
	open(51,file=inFile)
	write(dumpUnit,*)'inFile=',TRIM(inFile)

	outFile=TRIM(masterDir)//'\output\HEDb.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\HEDb.txt'
	open(661,file=outFile)
	write(661,*)' ' ! purge old contents of KijOut. KijVal uses APPEND
	if(LOUDER)print*,'HEDB: outFile=',TRIM(outFile) 
	if(LOUDER)write(6,*)'    ID1	  ID2  nConv      Kij       stdErr   %PAADP    rmsErrP    dStdErr'

	ier=0
    iSystem=0
	iTimeStart=time()
	paadTot=0
	nConvergedTot=0
	write(6  ,*)'iSystem,id(1),id(2),T,P,kij(1,2),HeCalcMin,HeCalcMax,HeDatMin,HeDatMax,(refid)'
	write(661,*)'iSystem,id(1),id(2),T,P,kij(1,2),HeCalcMin,HeCalcMax,HeDatMin,HeDatMax,(refid)'
	do while(ier==0) !loop over data sets in db
		!READ(51,'(4i5,a,a,a77)',ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
		!pause 'HEDb: Reading new mix'
		READ(51,*,ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
		if(LOUDER)write(*,*)'HEDb: nPts,ids=',nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
        if(ioErr /= 0)then	! shut down and stop the clock.
			timeTot=float(time()-iTimeStart)
			print*,'HEDb: elapsed time(sec)=',timeTot
			if(iSystem > 0)paadTot=paadTot/iSystem
			write(*,'(a,i9,2f10.2)')' Overall nConverged,paad=',nConvergedTot,paadTot
			close(661)
			close(51)
			if(ioErr== -1)then
				pause 'HEDb: end of file.'
			else
				print*,'HEDb: read error. system#,ioErr,File=',iSystem,ioErr,TRIM(inFile)
				pause 'HEDb: check file.'
			endif
			exit ! end of file reached
		endif
        iSystem=iSystem+1
		if(nPtsBipDat.le.0)then
			ier=1 !indicates end of database
			cycle !loop to end will terminate
		endif
		idOpt=1	! USE id
		call PGLStartup(NC,iEosOpt,idOpt,iErrLookup) !(NC,idCas,iErrLookup,errMsgSub)	! idDippr passed by GlobConst. This also sets the class()
		if(iErrLookup > 0)write(*,*)'HEDb: Error from IdCasLookup=',iErrLookup
		write(*,*)'HEDb: New mix, NPTS,ID()=',nPtsBipDat,id(1),id(2)
        if(LOUDER)pause 'Starting new binary mixture'
		if( iErrLookup > 0 )then
			write(*,*)'HEDb: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(661,'(2i10,3f11.4)')id(1),id(2),0,zero,eightySix
			DO I=1,nPtsBipDat	!read the points but do nothing so we can cycle to next
				READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)
			enddo
			cycle !loop to next component but leave ier alone so keep looping.
        endif
		if(bESD)write(*,*)'HEDb:kij=',kij(1,2)
        !if(LOUDER)write(*,'(a,11F10.2)')' Vx=',(vx(i),i=1,NC)
		HeDatMin=  1E11
		HeDatMax= 0
		DO I=1,nPtsBipDat	!we only get to here if there is no error.  otherwise, cycle takes to enddo while()
			READ(51,'(a123)')dumString
			READ(dumString,*)tDat(I),XDAT(I),YDAT(I),PDAT(I)	! place data in module BIPs
			Call ReadRef(dumString,refID) ! ID() or idCas() USEd from GlobConst
			!write(refID,*)dumString(length-14:length)
			if(xDat(i) <   zeroTol)xDat(i)=zeroTol*10	! 10x to avoid y<zeroTol issues.
			if(xDat(i) > 1-zeroTol)xDat(i)=1-zeroTol*10	! 10x to avoid y<zeroTol issues.
			if(ydat(i) < HeDatMin)HeDatMin=ydat(i)
			if( ABS(ydat(i)) > ABS(HeDatMax) )HeDatMax=ydat(i)
			DO iComp=1,NC
				if(Tc(iComp) < zeroTol)pause 'Main: Tc < 0 for some comp???'
				TR(iComp)=TDAT(I)/TC(iComp)
			ENDDO
		enddo !read dat file
		!write(*,'(a,i4,a,a)')' LEN    ,RefID=',length,' ',dumString(length-14:length)
		!write(*,'(a,i4,a,a)')' idStart,RefID=',idStart,' ',TRIM(refID)
		!pause 'check refid'
		T=tdat(1) ! assume constant T,P
		P=pdat(1)
		smallx=1.D-6
		x(1:NC)=smallx
		LIQ=1
		do iPure=1,NC
			x(iPure)=1-smallx*(NC-1)
			!call Fugi(t,p,x,NC,1,fugcL,zL,ier)
			call FugiTP( T,P,x,NC,LIQ,rhoMol_cc,zL,aRes,fugcL,uRes,iErrF )
			!bMix=sumproduct(NC,x,bVolCc_mol) ! no need to calculate this. bMix=bVol(iPure)
			!if(LOUDER)write(dumpUnit,form610)' HeDB: GpureAssoc=',fugAssocPure(1:NC)
			Vcc_mol(iPure)=1/rhoMol_cc
			hPure(iPure)=uRes+zL-1
			if(zL < zeroTol)write(dumpUnit,*)'HeDB: 0>zL=',zL
			gPure(iPure)=aRes+zL-1-LOG(zL)
			!if(LOUDER)write(*,form611)' i,V,H,G,fAss,uAss=',iPure,Vcc_mol(iPure),hPure(iPure),gPure(iPure),fugAssocPure(iPure),uAssoc
			x(iPure)=smallx
		enddo ! iPure
		zero=0
		one=1
		if(LOUDER)write(dumpUnit,*) ' x(1),hRes,V(cc/mol),gRes'
		if(LOUDER)write(*,form600)zero,hPure(2),Vcc_mol(2),gPure(2)
		if(LOUDER)write(*,form600)one ,hPure(1),Vcc_mol(1),gPure(1)

		
		nXstep=21
		dX=1.d0/(nXstep-1)
		LIQ=1
		HeCalcMin=  1E11
		HeCalcMax= 0
		write(dumpUnit,*) ' x(1),hEkJ_mol,vE,gEmix,His,gMix,gIs'
		do iX1=1,nXstep-2
			x(1)=iX1*dX
			x(2)=1-x(1)
			!call Fugi(t,p,x,NC,1,fugcL,zL,ier)
			call FugiTP( T,P,x,NC,LIQ,rhoMol_cc,zL,aRes,fugcL,uRes,iErrF )
			vTotCc=1/rhoMol_cc
			gIs=0
			bMix  =SumProduct( NC,x,bVolCC_mol )
			eta=bMix/vTotCc
			Gis  =SUM( x(1:NC)*gPure(1:NC) ) + SUM(  x(1:NC)*LOG( x(1:NC) )  )
			hIs  =SUM( x(1:NC)*hPure(1:NC) ) !SumProduct(NC,x,hPure )
			vIs  =SUM( x(1:NC)*Vcc_mol(1:NC) ) !SumProduct(NC,x,hPure )
			gMix=aRes+zL-1-LOG(zL)
			gEmix=gMix - gIs
			hRes=uRes+zL-1 
			hE=hRes-hIs
			hEkJ_mol=hE*Rgas*T/1000
			vE=1/rhoMol_cc-vIs
			if(hEkJ_mol < HeCalcMin)HeCalcMin=hEkJ_mol
			if( ABS(hEkJ_mol) > ABS(HeCalcMax) )HeCalcMax=hEkJ_mol
			write(*,'(f7.4,11f11.4)')x(1),hEkJ_mol,vE,gEmix,His,gMix,gIs
		enddo  !iX1
		write(6  ,'( 3i5,7(1x,f8.3),1x,a )')iSystem,id(1),id(2),T,P,kij(1,2),HeCalcMin,HeCalcMax,HeDatMin,HeDatMax,TRIM(refid)
		write(661,'( 3i5,7(1x,f8.3),1x,a )')iSystem,id(1),id(2),T,P,kij(1,2),HeCalcMin,HeCalcMax,HeDatMin,HeDatMax,TRIM(refid)
	enddo !while(ier.eq.0) loop over all data sets in db

	print*,'HEDb: outFile=',TRIM(outFile) 
	pause 'These results are tabulated in HEDb.txt.'
	!restore stored nc,id
	NC=ncStore
	do i=1,NC
		id(i)=idStore(i)
	enddo
	CALL GETCRIT(NC,iErrCrit)
	call IdCasLookup(NC,idCas,ier,errMsgSub)	! idDippr comes from GlobConst. This also resets the class().
	if(iErrCrit.ne.0)then
		if(LOUDER)write(*,*)'Error in Main: ErrorCode from GetCrit = ',iErrCrit
		if(LOUDER)pause
		stop
	endif 
	if(iEosOpt.eq.1)CALL GetPR(NC,iErrGet)
	if(iEosOpt==2)CALL GetEsdCas(NC,idCas,iErrGet)
	if(bESD)CALL GetEsdCas(NC,idCas,iErrGet)
	if(iEosOpt==3)CALL GetPRWS(NC,iErrGet)
	if(bTPT)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
	if(iEosOpt==6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
	if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
	if(iErrGet >0)then
		if(LOUDER)write(*,*)'KijDB Error: failed to restore required pure props.'
		if(LOUDER)write(*,*)'id1,id2',id(1),id(2)
		if(LOUDER)pause
	endif
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
608	FORMAT(1X,3i6,1X,E11.4,1X,E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
      
	RETURN
	END	   ! HEDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MemScedOpt(NC)
	!C
	!C  PURPOSE:  IDAC CALCULATIONS TO DETERMINE THE OPTIMAL
	!C    VALUE FOR THE Acid/Base parameters FOR A SERIES
	!C    OF EXPERIMENTAL DATA SPECIFIED BY USER THROUGH PROMPT.
	!C  NOTE:  OPTIMAL KIJ IS PASSED THROUGH COMMON STATEMENT
	!C  PROGRAMMED BY:  JRE 2/96
	!C  METHOD:   GOLDEN SECTION SEARCH ON A SINGLE PARAMETER
	!C  REFERENCE:NUMERICAL RECIPES
	!C
	USE PortLib !for time()
	USE GlobConst
	USE BIPs  !includes nConverged,maxPts,tDat,...
	USE EsdMem2ParmsDb	!includes idTau,idBeta,...
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(nTau+nBeta+nAlpha+2),stdErr(nTau+nBeta+nAlpha+2),temParm(nTau+nBeta+nAlpha+2) !,TR(NMX)
	Integer idStore(NMX) !,idcc(2),QueryNParMix !,iErrGet,idCas(NMX)
	LOGICAL LOUDER
	EXTERNAL MemScedEval
	character*77 inFile,outFile !,errMsgPas,errMsgSub !,dumString ,name1,name2
	!character*1 tabChar
	!tabChar = char(9)
	LOUDER=LOUD	!This provides local control if desired.
	LOUDER=.TRUE. !Uncomment this line for screen writes from this routine.
	!store current values of nc,id.
	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	NC=2
	zero=0
	eightySix= -86
	iAns=1

	ier=0
    iSystem=0
	iTimeStart=time()
    FN="\IdacRecLazzaroniDb.txt"
	inFile=TRIM(PGLinputDir)//TRIM(FN)
	open(501,file=inFile)
	write(dumpUnit,*)'inFile=',TRIM(inFile)
	read(501,*,ioStat=ioErr)nPts
	do i=1,nPts
		read(501,*,ioStat=ioErr)id2Dat(i),id1Dat(i),tDat(i),gamDat(i)
	enddo
	close(501)
	if(.NOT.isReadEsd)call LoadEsdDb(iErrLoad)	! use ESDMEM2 estimates for initial guesses.
	if(iErrLoad > 0)pause 'MemScedOpt: error loading ESD DB.'
	nParms=nBeta+nAlpha+nTau +1 
	!if(nTau > 0)parm(1:nTau)=0	! vector init.
	do i=1,nBeta
		parm(i)=epsD_kBdb( indexEsd(idBeta(i)) ) ! this will take initial guess from ESDMEM2 database
	enddo
	do i=1,nAlpha
		parm(i+nBeta)=epsA_kBdb( indexEsd(idAlpha(i)) ) ! this will take initial guess from ESDMEM2 database
	enddo
	do i=1,nTau
		j=i+nBeta+nAlpha
		parm(j)=tauDb( indexEsd(idTau(i)) ) ! this will take initial guess from ESDMEM2 database
		if(parm(j) < 1)parm(j)=1			! start with a finite value so scale checking will have some effect.
	enddo
	parm(nParms)=0.02d0 ! scale factor for dTau^2
	iFlag=1	! iFlag=1 means don't print all results. Just get the deviations for converged points.
	call MemScedEval(nPts,nParms,parm,deviate,iFlag)	!nConverged passed by BIPs module.
	rootPts=SQRT(FLOAT(nPts))
	rmsErr0=ENORM(nPts,deviate)/rootPts
	write(*,*)' MemScedOpt: nParms,rmsErr to start=',nParms,rmsErr0
	!pause ' MemScedOpt: check initial error.'
	temparm(1:nParms)=parm(1:nParms)
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	do iter=1,1		! Try increasing or decreasing twice each.
		do i=1,nParms
			factor=1.1d0
			if(isOdd(iter)==1)factor=0.9d0
			temparm(i)=parm(i)*factor
			call MemScedEval(nPts,nParms,temparm,deviate,iFlag)	!nConverged USEd in BIPs.
			rmsErr=( ENORM(nPts,deviate)/rootPts )
			if(rmsErr < rmsErr0)then
				parm(i)=temparm(i)
				rmsErr0=rmsErr
			else
				temparm(i)=parm(i)
				if( ABS(rmsErr-rmsErr0)/rmsErr0 < 1D-5 )then
					idi=idBeta(i)
					if(i > nBeta)idi=idAlpha(i-nBeta)
					if(i > nBeta+nAlpha)idi=idTau(i-nBeta-nAlpha)
					open(661,file=outFile,access='APPEND')
					write( * ,form612)' MemScedOpt: No change for i,idi,parm=',i,idi,parm(i)
					write(661,form612)' MemScedOpt: No change for i,idi,parm=',i,idi,parm(i)
					close(661)
				endif
			endif
		enddo ! nParms
	enddo ! iter
	do i=1,nParms
		write(*,form601)i,parm(i)
	enddo
	pause 'MemScedOpt: Check crude result.'
	TOL=0.00001	 !c      TOL=SQRT(DPMPAR(1))
	factor= 0.100 !c  smaller "factor" means smaller first step from the initial guess.
	print*, 'MemScedOpt: Calling LmDifEz'
	CALL LMDifEZ(MemScedEval,nPts,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
	print*, 'MemScedOpt: Back from LmDifEz'
	iFlag=86	! get the deviations for converged points. Print everything. 
	call MemScedEval(nPts,nParms,parm,deviate,iFlag)	!nConverged passed by BIPs module.
	do i=1,nParms
		idi=idBeta(i)
		if(i > nBeta)idi=idAlpha(i-nBeta)
		if(i > nBeta+nAlpha)idi=idTau(i-nBeta-nAlpha)
		open(661,file=outFile,access='APPEND')
		write( * ,form612)' MemScedOpt: i,idi,parm=',i,idi,parm(i)
		write(661,form612)' MemScedOpt: i,idi,parm=',i,idi,parm(i)
		close(661)
	enddo
	rmsErr= ENORM(nPts,deviate)/rootPts
	write(*,form611)TRIM(outFile),iErrCode,rmsErr
	pause 'MemScedOpt: Success! Check KijOut.txt'
	!restore stored nc,id
	NC=ncStore
	ID(1:NC)=idStore(1:NC)						  
	idOpt=1
	call PGLStartup(NC,iEosOpt,idOpt,iErrStart)	! sets idCas, GetCrit(), Get_Eos_
	if(iErrStart >0)then
		if(LOUDER)write(*,*)'KijDB Error: failed to restore required pure props.'
		if(LOUDER)write(*,*)'id1,id2',id(1),id(2)
		if(LOUDER)pause
	endif
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
606	FORMAT(1X,F8.4,1X,F10.2)      
607	FORMAT(1X,3i6,1X,<nParms>E11.4,1X,<nParms>E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
608	FORMAT(1X,3i6,1X,E11.4,1X,E11.4,f10.2,f10.2,f10.4,1x,a20,1x,a20)
	RETURN
END	   ! MemScedOpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MemScedEval(nPts,nParms,parm,deviate,iFlag)
	USE GlobConst
	USE EsdMem2ParmsDb	!idBeta,... nMemSced,nTau,...
	USE BIPs ! includes nConverged,Kij(),KTij(), ...
	USE VpDb ! for VpCoeffs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	parameter(MAXDAT=4444,MAXPARMS=221)	! From LmDifEzCov2. 
	CHARACTER*77 outfile
	DOUBLEPRECISION, SAVE:: rLnGam(NMX),temParm(maxParms) ! X(NMX),Y(NMX) !,yy(nmx),xx(nmx),
	integer idOld(2) !ier(20),
	DOUBLE PRECISION parm(nParms),deviate(nPts) !,temParm(nParms)
	character*77 errMsg(0:22),errMsgPas !,outFile
	LOGICAL LOUDER
	data nCalls,rmsErr0/1,1E11/
	errMsg(0) =' MemScedEval: no problem.' 
	errMsg(11)=' MemScedEval: experimental gamma < 0???'
	iErr=0
	iPrint=0
	if(iFlag==0)then
		iPrint=1
		LOUDER=.TRUE.
	else
		LOUDER=.FALSE.
	endif

	errMsg(1)='MemScedEval: Number of components must equal 2.'
	if(LOUDER)write(*,*)'MemScedEval: Starting.'
	!if(nTau > 0)tau(1:nTau)=parm(1:nTau)	! start by setting tau=parm for all compds.
	do i=1,nBeta+nAlpha
		if(parm(i) < 10)parm(i)=0		! Let's neglect small values. Just not worth it. < 0 NOT allowed.
		if(parm(i) > 1D4)parm(i)=1D4	! Don't let parm's get too big or expo's might break.
	enddo
	do i=1,nTau
		j=i+nBeta+nAlpha
		tauDb(  indexEsd( idTau(i) )  )=parm(j)
	enddo
	do i=1,nBeta
		epsD_kBdb(  indexEsd( idBeta(i) )  )=parm(i)
	enddo           
	do i=1,nAlpha
		epsA_kBdb(  indexEsd( idAlpha(i) )  )=parm(nBeta+i)
	enddo           
	do i=1,nAssoc
		epsA_kBdb(  indexEsd( idAssoc(i) )  )=2*eTotAssoc(i)-epsD_kBdb(  indexEsd( idAssoc(i) )  )
	enddo           
	NC=2
	KTIJ(1,2)=0								! setting kij=0 means we get gamPhys from ESD defaults.
	HIJ(1,2)=0
	HTIJ(1,2)=0
	xsTau(1,2)=0
	xsTau(2,1)=0
	xsAlpha(1,2)=0.3d0
	xsAlpha(2,1)=xsAlpha(1,2)
	idOld(1:NC)=0
	idOpt=2		! Specifies that we use idCas for PGLStartup.
	iEosOpt=4	! Specifies ESDMEM2 as foundation.
	pMPa=111	! Maybe more than necessary but we want to keep it liquid.
	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	is661open=0
	if(iFlag==86.or.nCalls==1)then
		open(661,file=outFile,ACCESS='APPEND')
		write(661,form600)parm(1:nParms)
		is661open=1
	endif
	rmsErr=0
	nConverged=nPts
	do i=1,nPts
		idCas(1)=id1Dat(i)	! id1Dat and id2Dat USEd from BIPs
		idCas(2)=id2Dat(i)
		if(iFlag < 86)then
			deviate(i)=300
		else
			deviate(i)=0
		endif
		if(gamDat(i) > 0)then
			rLnGamExpt=LOG( gamDat(i) )	! USEd from BIPs
		else
			iErr=11
			cycle
		endif
		if( idCas(1).ne.idOld(1) .or. idCas(2).ne.idOld(2) )then
			idOld(1:NC)=idCas(1:NC)
			call PGLStartup(NC,iEosOpt,idOpt,iErrStart)	! sets idCas, ID, GetCrit(), Get_Eos_
			if(iErrStart.ne.0)pause 'MemScedEval: error from PGLStartup?'
		endif
		dTau=tau(2)-tau(1)
		!if(nTau > 0)dTau=tau(  indexEsd( ID(2) )  )-tau(  indexEsd( ID(1) )  )
		KIJ(1,2)=dTau*dTau/100*0.045d0 +dTau/10*parm(nParms)
		KIJ(2,1)=Kij(1,2)
		!if(LOUDER)write(*,form611)' To IDACsTP. iPt,a1,b1,a2,b2=',i,epsA_kBdb(indexEsd(ID(1))),epsD_kBdb(indexEsd(ID(1))),epsA_kBdb(indexEsd(ID(2))),epsD_kBdb(indexEsd(ID(2)))
		call IDACsTP(tDat(i),pMPa,NC,rLnGam,iErrIdacs,errMsgPas)
		!if(LOUDER)write(*,form612)' From IDACsTP. id1,id2,rLnGam2=',ID(1),ID(2),rLnGam(2)
		if(iErrIdacs.ne.0)then
			nConverged=nConverged-1
			cycle
		endif
		deviate(i)=rLnGam(2)-rLnGamExpt
		rmsErr=rmsErr+deviate(i)*deviate(i)
		if(iFlag==86.or.nCalls==1)then
			write(661,form602)idCas(2),idCas(1),tDat(i),rLnGamExpt,rLnGam(2),deviate(i)
		endif
	enddo
	if(is661open==1)close(661)
	sqArg=rmsErr/nConverged
	if(sqArg < zeroTol)pause ' MemScedEval: SSE < 0?'
	rmsErr=SQRT(rmsErr/nConverged)
	write(*,'(a,2i5,2f15.9,7f12.4)')' nCalls,nConverged,rmsErr=',nCalls,nConverged,rmsErr
	if(rmsErr < rmsErr0)then
		temParm(1:nParms)=parm(1:nParms)
		rmsErr0=rmsErr
		open(661,file=outFile,access='append')
		write(661,form600)rmsErr0,parm(1:nParms)
		close(661)
	endif
	if(nCalls > 1E11)then
		if(LOUDER)pause 'MemScedEval: nCalls > 100. Slow down.'
		if(LOUDER)write(*,form600)parm(1:nParms)
	endif
	if(iFlag==86.or.iPrint==1)then
		open(661,file=outFile,access='append')
		write(661,'(a,2i5,2f10.3,7f12.4)')' nCalls,nConverged,rmsErr=',nCalls,nConverged,rmsErr0
		write(661,form600)temParm(1:nParms)
		if(LOUDER)print*,'KijDB: outFile=',TRIM(outFile) 
		if(LOUDER)print*,'These results are tabulated in outFile.'
		close(661)
	endif
	nCalls=nCalls+1
	if(iFlag==86)iFlag=iErr	! this is the only way to return error info based on LmDiff code.
	RETURN
END ! MemScedEval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER*44 FN
	!PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5),temParm(5),TR(NMX)
	Integer idStore(NMX),idcc(2),QueryNParMix !,iErrGet,idCas(NMX)
	LOGICAL LOUDER
	EXTERNAL LMDevFcn
	character*77 inFile,outFile !,errMsgPas !errMsgSub,dumString ,name1,name2
!	character*1 tabChar
	!COMMON/eta/etaL,etaV,ZL,ZV
	!common/FloryWert/ vLiq(nmx)
    data initial/1/
	!tabChar = char(9)
	LOUDER=LOUD	!This provides local control if desired.
	LOUDER=.TRUE. !Uncomment this line for screen writes from this routine.
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
	print*,'Add +10 to skip opt or 0 to type FN.FT'
	read(*,*)iAns
	iOptimize=1
	if(iAns > 10)then
		iOptimize=0
		iAns=iAns - 10
	elseif(iAns < 0)then
		iOptimize= -1	  ! it means just read the database and rewrite in new format.
		iAns=ABS(iAns)
	endif
    FN="\VLEJaubert.txt"
	if(iAns==2)FN="\VLEJaubertAq.txt"
	if(iAns==3)FN="\VLEDannerGess.txt"															
	if(iAns==4)FN="\SolvationSystems.txt"
	if(iAns < 1)then
		write(*,*)'Enter FN.FT (path is assumed to be PGLInputDir.'
		read(*,*)FN
		FN='\'//TRIM(FN)
	endif
!	READ(*,'(A)')FN
	inFile=TRIM(PGLinputDir)//TRIM(FN)
	if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
	open(51,file=inFile)
	write(dumpUnit,*)'inFile=',TRIM(inFile)

	outFile=TRIM(masterDir)//'\output\KijOut.txt'
	IF(DEBUG)outFile='c:\spead\calceos\output\KijOut.txt'
	open(661,file=outFile)
	write(661,*)' ' ! purge old contents of KijOut. KijVal uses APPEND
	if(iOptimize .ge. 0)close(661) ! purge old contents of KijOut. KijVal uses APPEND.

	outFile=TRIM(masterDir)//'\output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\output\KijDb.txt'
	if(initial==1)then
        initial=0
        open(674,file=outFile)
    else
        open(674,file=outFile,ACCESS='APPEND') ! e.g. JaubertAq is appended to Jaubert if during the same run.
    endif
	if(LOUDER)print*,'KijDB: outFile=',TRIM(outFile) 
	if(LOUDER)write(6,*)'    ID1	  ID2  nConv      Kij       stdErr   %PAADP    rmsErrP    dStdErr'
	write(674,*)'    ID1	  ID2  nConv      Kij       stdErr   %PAADP    rmsErrP    dStdErr'

	if(iEosOpt.eq.5.or.iEosOpt.eq.8.or.iEosOpt.eq.9)then	   !AFG 2011 , Added EOS Opt. 9
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
	do while(ier==0) !loop over data sets in db
		!READ(51,'(4i5,a,a,a77)',ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
		!pause 'KijDb: Reading new mix'
		READ(51,*,ioStat=ioErr)nPtsBipDat,id(1),id(2) !,iRefNum,name1,name2,dumString !,idCas(1),idCas(2)
        if(ioErr /= 0)then
			timeTot=float(time()-iTimeStart)
			print*,'KijDb: elapsed time(sec)=',timeTot
			if(iSystem > 0)paadTot=paadTot/iSystem
			write(*,'(a,i9,2f10.2)')' Overall nConverged,paad=',nConvergedTot,paadTot
			close(674)
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
		if(iOptimize < 0)then
			!write(661,'(i4,2i5,2i6,a)')nPtsBipDat,id(1),id(2),idTrc(1),idTrc(2),' ! Line0:nPts,idDip1,idDip2,idTrc1,idTrc2; Line1-nPts:T(K),P(MPa),x1,y1' 	
		endif !iOptimize < 0

		if(idOpt.eq.2)then
			do i=1,NC
				idcc(i)=id(i)
			enddo
		endif
		write(*,*)'KijDb: New mix, NPTS,ID()=',nPtsBipDat,id(1),id(2)
        !if(LOUDER)pause 'Starting new binary mixture'
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
		idOpt=1	!ID(dippr) is USEd.
		call PGLStartup(NC,iEosOpt,idOpt,iErrStart)	! sets idCas, GetCrit(), Get_Eos_
		if(iErrStart > 0 .or. iOptimize < 0)then
			write(*,*)'KijDB Error: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(674,607)id(1),id(2),0,zero,eightySix
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
			READ(51,*)tDat(I),PDAT(I),XDAT(I),YDAT(I)	! place data in module BIPs
			if(xDat(i) <   zeroTol)xDat(i)=zeroTol*10	! 10x to avoid y<zeroTol issues.
			if(xDat(i) > 1-zeroTol)xDat(i)=1-zeroTol*10	! 10x to avoid y<zeroTol issues.
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
		if(iErrParMix)pause 'KijDb: error returning stored mix parms'
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
					parm(1)= -0.15d0+(iGrid-1)*0.1d0	!search from -0.15 to +0.15
					if(iEosOpt==19)then
						parm(1)=(parm(1)+0.1)*1000	! search from -50 to 250.
						parm(2)=parm(1)			!assume equal values for grid search.
					endif
					call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,0)	!nConverged passed by BIPs module. 0 means enforce penalty for nonconverged points.
					if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
					write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
					if(PAADP < bestAAD)then
						bestAAD=PAADP
						bestKij=parm(1)
					endif
				enddo
				parm(1)=bestKij
				if(nParms > 1)parm(2)=parm(1)
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
		!if(LOUDER)pause 'KijDB: Converged. Last call to KijVal.'
		LastCall=1	! Do not apply penalty for unconverged points on LastCall. Just get the deviations for converged points.
		call KIJVAL(nc,nParms,parm,deviate,PAADP,DstdErr,LastCall)	!nConverged passed by BIPs module.
		if(nConverged > 0)rmsErr=SQRT( DFLOAT(nPtsBipDat)/nConverged )*RmsJre(nPtsBipDat,deviate) !adjust for nConverged.nPts needed to sum over all elements.
		!write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(LOUDER)print*,'KijDB: back from KijVal.'
		if(nParms==1)write(674,608)ID(1),ID(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==2)then
			temParm(1)=parm(2)	! parameter storage is order sensitive for xsTau(i,j). Storage applies ID(1) > ID(2) such that IdBinary=10000*ID(1)+ID(2).
			temParm(2)=parm(1)
			if(ID(1) > ID(2) )temParm(1:2)=parm(1:2)  
			write(674,608)ID(1),ID(2),nConverged,(temParm(i),i=1,2),PAADP,rmsErr,DstdErr,name(1),name(2)
			write(6  ,608)ID(1),ID(2),nConverged,(temParm(i),i=1,2),PAADP,rmsErr,DstdErr,name(1),name(2)
		endif
		write(6 ,'( 3i5,2F10.2,11(1x,E11.4) )')iSystem,id(1),id(2),PAADP,rmsErr,(parm(i),i=1,nParms)
		paadTot=paadTot+paadP
		nConvergedTot=nConvergedTot+nConverged
	enddo !while(ier.eq.0) loop over all data sets in db

	close(51)
	!close(61)
	if(LOUDER)print*,'KijDB: outFile=',TRIM(outFile) 
	if(LOUDER)pause 'These results are tabulated in KijDb.txt.'
	!restore stored nc,id
	NC=ncStore
	ID(1:NC)=idStore(1:NC)
	call PGLStartup(NC,iEosOpt,idOpt,iErrStart)	! sets idCas, GetCrit(), Get_Eos_
	if(iErrStart >0)then
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
	Integer idStore(NMX),QueryNParMix !,sortIndex(maxPts) !,iErrGet,idcc(2),idCas(NMX)
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
	write(*,'(a93)')'Enter 2 to optimize BIPs with flash, 1 to optimize with Ks, or 0 to use stored BIPs and Ks.'
	print*,'Enter -ve to optimize with flash using a custom database file.'
	read*,iOptimize
    FN="\LLeDbPGL6ed96b.txt"
	if(iOptimize < 0)then
		WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in PGLInputDir)'
		READ(*,'(A)')FN
		FN='\'//TRIM(FN)
		iOptimize=2
	endif
	inFile=TRIM(masterDir)//'\input'//TRIM(FN)
	if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
	inUnit=51
	print*,'inFile=',TRIM(inFile)
	open(inUnit,file=inFile)
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
        open(654,file=outFile)
    else
        open(654,file=outFile,ACCESS='APPEND') ! e.g. JaubertAq is appended to Jaubert if during the same run.
    endif
	write(654,*)'   ID1   ID2  nConv      Kij       stdErr    %AALDK    %BIASLK    dStdErr'

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
		print*,'System Header= ',TRIM(dumString)
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
			print*,'KijDbLLE: ID fail! nCritSet???',idCas(1),idCas(2),ID(1),ID(2)
			print*,TRIM(errMsgLookup)
			cycle
		endif
		!print*,'ID=',ID(1),ID(2)
		!call IdCasLookup(NC,idCas,ierLookup,errMsgLookup)	
		!if(LOUDER)pause 'Starting new binary mixture'
		CALL GETCRIT(NC,iErrCrit)
		write(*,'(a,i4,4a)')' KijDbLLE: after GetCrit, iEosOpt,Names=',iEosOpt,' ',TRIM(name(1)),' ',TRIM(name(2))
		if(iErrCrit.ne.0)then
			print*,'ID()=',ID(1),ID(2)
			write(*,form611)' KijDbLLE Error: iErrCit,Tc()=',iErrCrit,Tc(1),Tc(2)
			!pause
			cycle
		endif
		iErrGet=SetNewEos(iEosOpt) !wipe out any instance of other eos.
		print*,'KijDbLLE: after SetNewEos, iEosOpt,nParms=',iEosOpt,nParms
		if(iEosOpt==1)CALL GetPR(NC,iErrGet)
		if(iEosOpt==2)CALL GetEsdCas(NC,idCas,iErrGet)
		if(iEosOpt==4)CALL GetEsdCas(NC,idCas,iErrGet)
		if(iEosOpt==3)CALL GetPRWS(NC,iErrGet)
		if(iEosOpt==5)CALL GetTpt(NC,ID,iErrGet,errMsgPas)	!AFG 2011 , Added EOS Opt. 9
		if(iEosOpt==6)CALL GetFloryWert(NC,ID,iErrGet)!Results placed in common: TptParms, HbParms
		if(iEosOpt==7)CALL GetNRTL (NC,ID,iErrGet)
		if(iEosOpt==10)CALL GetPcSaft(NC,idCas,iErrGet)
		if(iEosOpt==11)CALL GetPrtc(NC,iErrGet)
		if(iEosOpt==17)CALL GetPrLorraine(NC,iErrGet)
		if(iEosOpt==19)CALL GetLsgMem2(NC,ID,iErrGet)
		print*,'KijDbLLE: after Get_EOS_, iEosOpt,bESD,bTPT,bPcSAFT=',iEosOpt,bESD,bTPT,bPcSAFT
		if(iErrGet > 0 .or. iErrCrit > 0 .or. ierLookup > 0)then
			write(*,*)'KijDB Error: failed to get required pure props.'
			write(*,*)'id1,id2',id(1),id(2),name(1),name(2)								 
			write(654,607)id(1),id(2),0,zero,eightySix
			cycle
		endif
		iPt=0
		avgx=0		! LLEFL assumes that xDat < yDat so we need to switch the order if that's not true.
		nx=0
		avgy=0
		ny=0
		do i=1,nLines ! parse the lines into data
			if(var(i,3).ne.1)then
				print*,'var()=',(var(i,j),j=1,nVars)
				pause 'KijDbLLE: XDAT =/= X1???'
			endif
			if( iOptimize ==1 .and. (var(i,4)<0.or.var(i,5)<0) )cycle ! for iOptimize==1, we require coexisting x,y and subcritical.
			!print*,'KijDbLLE: x1Lo,x1Up=',var(i,4),var(i,5)
			iPt=iPt+1
			XDAT(iPt)=var(i,4)
			if(XDAT(iPt)>0)avgx=avgx+XDAT(iPt)
			if(XDAT(iPt)>0)nx=nx+1
			YDAT(iPt)=var(i,5)
			if(YDAT(iPt)>0)avgY=avgY+YDAT(iPt)
			if(YDAT(iPt)>0)nY=nY+1
			TDAT(iPt)=var(i,1)
			PDAT(iPt)=var(i,2)/1000 !convert from kPa to MPa.
			if(TDAT(iPt) > tMax)tMax=TDAT(iPt)
			if(TDAT(iPt) < tMin)tMin=TDAT(iPt)
		enddo
		nPtsBipDat=iPt
		avgx=avgx/nx
		avgy=avgy/ny
		if(nPtsBipDat==0)cycle !omit systems for which no valid points available. e.g., no coexistence when iOptimize==1
		if(avgx > avgy)then
			do i=1,nPtsBipDat
				xTemp=XDAT(i)
				XDAT(i)=YDAT(i)
				YDAT(i)=xTemp
			enddo
		endif
		! All DAT read in for this dataset. Time to regress.
		!write(*,*)'KijDbLLE: New mix, NPTS,ID(),TLast=',nPtsBipDat,id(1),id(2),TDAT(nPtsBipDat)
		nSystems=nSystems+1
		write(dumpUnit ,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax !,Tc(1),Tc(2)
		!write(61,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,Tc(1),Tc(2)
		!write(661,'(3i5,2i15,4f8.2)')nPtsBipDat,id(1),id(2),idCas(1),idCas(2),tMin,tMax,Tc(1),Tc(2)
		if(LOUDER)then
			do i=1,nPtsBipDat
				write(dumpUnit,'(f8.2,f8.4,2e12.4)')TDAT(i),PDAT(i),XDAT(i),YDAT(i)
			enddo
		endif
					
		write(dumpUnit,*)'New mixture. iSystem,id(),nPts=',nSystems,id(1),id(2),nPtsBipDat
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
		if(bESD)nParms=1
		parm(2)=0
		temParm=0	!array initialization.
		!xsTau=0		!array initialization.
		do i=1,nParms
			call QueryParMix(i,temParm(i),iErrParMix)  
			parm(i)=temParm(i) ! this will take initial guess from database, if available.
			if(iErrParMix>0)pause 'KijDbLLE: error from QueryParMix.'
		enddo
		!print*,'KijDbLLE: nParms,Stored BIPs=',nParms,parm(1:nParms)
		if(iErrParMix)print*,'KijDbLLE: error returning stored mix parms'
		write(dumpUnit,form612)' KijDbLLE: IDs,IniParms=',ID(1:2),(temParm(i),i=1,nParms)
		if(iEosOpt==86)then	! 86 means skip this for now. Change 86 to 5 for evaluating solvation effects.
			nParms=2 ! add one to adjust bipDA
			write(*,'(5x,11i7)')(idLocalType(i),i=1,nTypesTot)
			do i=1,nTypesTot
				write(*,'(i5,11f7.4)')idLocalType(i),(aBipDA(i,j),j=1,nTypesTot)
			enddo
			print*,'Enter i,j for site type optimization'
			read*,iDA,jDA
			parm(2)=aBipDA(iDA,jDA)
		endif
		if(iOptimize >0)then  !call LMDif
			MAXIT=111
			INIT=1
			LWA=275	 !C	LWA>5*nParms+NPTS ~ 275
			TOL=0.0001	 !c      TOL=SQRT(DPMPAR(1))
			factor= 0.2 !c  factor=100 enables large first step from the initial guess.
			if(nParms>1)factor=0.2
			write(*,form612)' IDs,ini parms=',ID(1:2),parm(1:nParms)
			!pause 'KijDbLLE: new mixture. check ini.'
			if(LOUDER)write(*,*)'KIJ,paaLDK,rmserr'
			if(ABS(temParm(1)) < zeroTol)then	! signal that parameters do not exist for this compound
				bestKij= 0.1
				bestP2 = 0.1
				bestAAD=12345
				mostConverged=0
				do iGrid=1,6
					parm(1)= (iGrid-1)*0.05d0-0.125	 ! grid from -0.125 to 0.025
					if(nParms>1)then
						parm(1)=(parm(1)+0.1)*1000
						parm(2)=parm(1) !-15
					endif
					call KIJDevLLe(nc,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,0)	!nConverged passed by BIPs module. 0 means enforce penalty for nonconverged points.
					if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
					write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAALDK,rmsErr,DstdErr,name(1),name(2)
					if(nConverged .ge. mostConverged)then
						if(nConverged > mostConverged)then
							bestKij=parm(1)
							if(nParms>1)bestP2=parm(2)
							bestAAD=PAALDK
						elseif(PAALDK < bestAAD)then		!This "else" condition occurs only when nConverged==mostConverged. 
							bestKij=parm(1)					!Only replace best if new PAAD is better.
							if(nParms>1)bestP2=parm(2)
							bestAAD=PAALDK
						endif
						mostConverged=nConverged
					endif
				enddo
				parm(1)=bestKij
				if(nParms>1)parm(2)=bestP2
			elseif(nParms>86)then
				parm(1:NC)=parm(1:NC)+25 ! kick initial guess up a little to force more LLE.
			endif !ParmSrch
			print*,'Initial parm()=',parm(1:nParms)
			!if(LOUDER)
			!pause 'Ready to regress'
			nCoexValues=2*nPtsBipDat ! Coexisting phases have 2 values each.
			CALL LMDifEZ(LMDevLLE,nCoexValues,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
			if(nConverged > 0)rmsErr=SQRT( ENORM(nConverged,deviate) )
			if(LOUDER)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		elseif(ABS(parm(1)) < zeroTol)then
			write(dumpUnit,*)'KijDbLLE: warning, parm(1)=0'
			parm(1)=0.0001d0
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
		call KIJDevLLE(NC,nParms,parm,deviate,PAALDK,DstdErr,PBIASK,LastCall)	!nConverged USEd in BIPs.
		!CALL LMDevLLE(Mdata,Nparms,parm,deviate,IFLAG)	! This function is different from KijDevLLE()
		if(nConverged > 0)rmsErr=SQRT( DFLOAT(nPtsBipDat)/nConverged )*RmsJre(nPtsBipDat,deviate) !adjust for nConverged.nPts needed to sum over all elements.
		!write(* ,607)id(1),id(2),nConverged,(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,rmsErr,DstdErr,name(1),name(2)
		if(nParms==1)write(6 ,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==1)write(654,608)id(1),id(2),nConverged,(parm(i),i=1,1),(StdErr(i),i=1,1),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(654,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		if(nParms==2)write(6 ,609)id(1),id(2),nConverged,(parm(i),i=1,2),PAALDK,PBIASK,DstdErr,name(1),name(2)
		write(6 ,'( 3i5,2F10.2,11(1x,E11.4) )')nSystems,id(1),id(2),PAALDK,PBIASK,(parm(i),i=1,nParms)
	enddo !while(ier.eq.0) loop over all data sets in db
	write(654,*)'nSystem total = ',nSystems
	write( * ,*)' #systems=',nSystems
	!write(661,*)' #systems=',nSystems
	close(51)
	close(654)
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
	inFile=TRIM(masterDir)//'\input\'//TRIM(FN)
	if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
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
	iSumUnit=672
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
			id1Dat(iPt)=NINT(var(i,3))!easier for general read as real, but we need integer for indexing.  
			if(TDAT(iPt) > tMax)tMax=TDAT(iPt)
			if(TDAT(iPt) < tMin)tMin=TDAT(iPt)
			CALL GetTmHfus(id(id1Dat(iPt)),TmKelvin,HfusJ_mol,iErrHfus)
			if( id(1)==1319)then
				print*,'KijDbSle: ID=1319, Hfus,iErrHfus=',HfusJ_mol,iErrHfus
				!pause 'check iErr for 1319'
			endif
			if(iErrHfus.ne.0 .or. TmKelvin < zerotol)then
				print*,'ID(iSolute)=',ID(id1Dat(iPt))
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
			write(dumpUnit,*)'KijDbSLE: warning, parm(1)=0'
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
	write(iSumUnit,*)'nSystem total = ',nSystems
	write( * ,*)' #systems=',nSystems
	!write(661,*)' #systems=',nSystems
	close(51)
	close(iSumUnit)
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
	inFile=TRIM(PGLinputDir)
	WRITE(*,*)'FILE NAME WITH EXPERIMENTAL DATA IN DIR:',TRIM(inFile)
	WRITE(*,*)'FIRST LINE =NPTS, 2ND LINE = TK, PMPA, X, ...'
	FN='Jaubert.txt'
	FN='vlemem2.txt'
	!READ(*,'(A)')FN
	!WRITE(52,*)' NAME OF FILE WITH EXPERIMENTAL DATA (PATH\FN.FT)'
	!write(52,'(1x,A)')FN
	inFile=TRIM(inFile)//'\'//TRIM(FN)
	!if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
	write(*,*)'KijOpt: inFile(full path)=',TRIM(inFile)
	open(51,file=inFile)
	READ(51,*)nPtsBipDat
	pDatMin=1E11
	pDatMax=0
	write(*,*)'     T           P           X1          Y1     '
	DO I=1,nPtsBipDat
		READ(51,*,ERR=861)tDat(I),PDAT(I),XDAT(I),YDAT(I)
		WRITE(*,form600)tDat(I),PDAT(I),XDAT(I),YDAT(I)
		if(PDAT(i) > pDatMax)pDatMax=PDAT(i)
		if(PDAT(i) < pDatMin)pDatMin=PDAT(i)
	enddo
	MAXIT=111
	INIT=1
	WRITE(*,*)'THIS ROUTINE ASSUMES A BINARY SYSTEM FROM MAIN PROGRAM'
	WRITE(6,*)'ENTER RANGE TO OPTIMIZE, (LO, HI)  '
	kij0= -0.15
	kij3=  0.15
	!READ(5,*)KIJ0,KIJ3                                          
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
	write(*,'(a,2f10.2,7f10.5)')' KijOpt: PAAD,rmsLnDev,parm()=',PAADP,rmse,(parm(i),i=1,nParms)

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
	if(iEosOpt==86)then	! change to iEosOpt==4 to optimize solvation energy
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
	if(LOUDER)print*,'KijDevLLE:parms=',(parm(i),i=1,nParms)

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
		if(iOptimize==1)print*,'tDat(i),LOG(CalcK1),LOG(exptK1),LOG(CalcK2),LOG(exptK2),dev1,dev2' 
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
		X(1)=XDAT(iData)
		X(2)=1-X(1)
		Y(1)=YDAT(iData) ! Y IS THE UPPER PHASE COMPOSITION
		Y(2)=1-Y(1)
		deviate(iData)=0
		if( ABS(XDAT(iData)-YDAT(iData)) < zeroTol)cycle ! That can't be for LLE.
		LIQ=1
		MAXIT=333
		if(iOptimize==2)then !  
			initialize=1	! Prefer to use input calcK as initial guess. e.g., previous converged point.
			if(iErrCode > 0 .or. iData==1)then !if there's a problem with previous point then...
				initialize=2 ! use spinodals as initial guess 
				if(xDat(iData) > 0 .and. yDat(iData) > 0 )then	 ! but if the current point represents coexistence, use expt to initialize
					calcK(1)=yDat(iData)/xDat(iData) ! use experimental value if coexistence available. 
					calcK(2)=(1-yDat(iData))/(1-xDat(iData)) ! for soluble compounds, use experimental value if possible.
					initialize=1
				endif
			endif 
			itMax=MAXIT
			Call LLEFL(T,P,NC,initialize,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode)
			if( ABS(x(1)-y(1)) < 0.06d0)iErrCode=86 ! sometimes convergence is indicated when it got too close to consolute.
			if(ITMAX > MAXIT-1)iErrCode=21
			if(iErrCode>0)Call LLEFL(T,P,NC,3,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode) ! If first try of LLEFL fails, try tangent line.
			if( ABS(x(1)-y(1)) < 0.06d0)iErrCode=86 ! sometimes convergence is indicated when it got too close to consolute.
			if(ITMAX > MAXIT-1)iErrCode=21
			dev1=0
			dev2=0 
			if(iErrCode==0)then
				if(xDat(iData) > 0)then
					x1Expt=0.5d0-ABS(0.5d0-xDat(iData))	!automatically takes the smaller
					x1Calc=0.5d0-ABS(0.5d0-x(1))
					dev1=DLOG( x1Calc/x1Expt )*100  
					deviate(iData) = dev1 ! use log deviate to minimize bias for bad fits.
					nValues=nValues+1
					nConverged=nConverged+1 ! Two dev's in above deviate. Consistent with one-sided evaluations when coexistence unavailable.
				endif
				if(yDat(iData) > 0)then
					y1Expt=0.5d0-ABS(0.5d0-yDat(iData))	!automatically takes the smaller
					y1Calc=0.5d0-ABS(0.5d0-y(1))
					dev2=DLOG( y1Calc/y1Expt )*100  
					deviate(iData+nPtsBipDat) = dev2 ! use log deviate to minimize bias for bad fits.
					nValues=nValues+1
					nConverged=nConverged+1 ! Two dev's in above deviate. Consistent with one-sided evaluations when coexistence unavailable.
				endif
				!print*,'got one! nConverged=',nConverged
				PAALDK = PAALDK+ABS(dev1)+ABS(dev2)
				PBIASK = PBIASK+(dev1+dev2)
				ssqErr= ssqErr+dev1*dev1+dev2*dev2
			else
				calcK(1)=1
				calcK(2)=1
				deviate(iData)=500 ! gas+water systems go to nConverged=0 if penalty=200
				deviate(iData+nPtsBipDat)=500 ! gas+water systems go to nConverged=0 if penalty=200
				if(LAST==1)then
					deviate(iData)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
					deviate(iData+nPtsBipDat)=0 ! no penalty for last evaluation. just compute the deviations for those that converged.
				endif
			endif !last=1 .and. compute activities
			if(LAST==1)WRITE(661,607)tDat(iData),x1Calc,x1Expt,y1Calc,y1Expt,dev1,dev2,itMax ! write deviations on LastCall
		elseif(iOptimize==1)then
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
		else ! iOptimize = 0 => compute LLE @ T if not optimizing.
			initialize=1  ! use the value of CALCK from expt to initialize flash iteration..
			if(iData==1)then !If the first point in dataset, we need to start fresh. Dataset should be sorted from LoT to HiT. 
				!initialize=3 !use tangent line
				calcK(1)=0.001	 ! x1=(1-K2)/(K1-K2) => x1=small if K1=large.

				if(xDat(iData) > 0 .and. yDat(iData) > 0)then
					calcK(1)=yDat(iData)/xDat(iData) ! use experimental value if coexistence available. 
					calcK(2)=(1-yDat(iData))/(1-xDat(iData)) ! for soluble compounds, use experimental value if possible. 
				else ! then yDat must be valid
					initialize=2
				endif
			endif
			itMax=151
			Call LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode)
			if(iErrCode>9)Call LLEFL(T,P,NC,3,CALCK,ITMAX,zFeed,X,Y,UOF,iErrCode)
			if( ABS(x(1)-y(1)) < 0.05d0)iErrCode=86 ! sometimes convergence is indicated when iteration crashed too close to consolute.
			if(iErrCode > 0)cycle
			nValues=0
			dev1=0
			dev2=0
			x1Expt=xDat(iData)
			y1Expt=yDat(iData)
			if(xDat(iData) > 0)then
				!dev1=DLOG( x(1)/xDat(iData) )*100
				!if( xDat(iData) > ABS(yDat(iData)) )dev1=DLOG( (1-x(1))/(1-xDat(iData)) )*100
				x1Expt=0.5d0-ABS(0.5d0-xDat(iData))
				x1Calc=0.5d0-ABS(0.5d0-x(1))
				dev1=DLOG( x1Calc/x1Expt )*100  !automatically takes the smaller
				nValues=nValues+1
			endif
			if(yDat(iData) > 0)then
				!dev2=DLOG( y(1)/yDat(iData) )*100
				!if( xDat(iData) < ABS(yDat(iData)) )dev2=DLOG( (1-y(1))/(1-yDat(iData)) )*100
				y1Expt=0.5d0-ABS(0.5d0-yDat(iData))
				y1Calc=0.5d0-ABS(0.5d0-y(1))
				dev2=DLOG( y1Calc/y1Expt )*100  !automatically takes the smaller
				nValues=nValues+1
			endif
			deviate(iData) = ( ABS(dev1)+ABS(dev2) ) ! don't divide by nValues here, Divide by nConverged after enddo..
			nConverged=nConverged+nValues
			PAALDK = PAALDK+deviate(iData) !/nValues
			PBIASK = PBIASK+(dev1+dev2) !/nValues
			ssqErr= ssqErr+deviate(iData)*deviate(iData)
			if(LAST==1)WRITE(661,607)tDat(iData),x1Calc,x1Expt,y1Calc,y1Expt,dev1,dev2,itMax ! write deviations on LastCall
			if(LOUDER )WRITE( 6 ,607)tDat(iData),x1Calc,x1Expt,y1Calc,y1Expt,dev1,dev2,itMax ! write deviations on LastCall
		endif  !iOptimze
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
	if(LOUDER)write(*,form613)' KijDevLLE:ID1,ID2,nCon,%AALDK,parms=',ID(1),ID(2),nConverged,PAALDK,rmsErr,parm(1:nParms)
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
	if(pDatMin > pMax)pMax=pDatMin+0.1 ! e.g. nC2+water the lowest pressure is large. if nConverged=0, LmDif crashes.           
	rmsErr=0
	LIQ=1 ! only liquids.
	if(LOUDER)print*,'KijDevSLE: Kij=',(parm(i),i=1,nParms)
	DO iData=1,nPtsBipDat !	data in module BIPs
		!This needs to be inside the data loop because the temperature may change.
		T=tDat(iData)
		iSolute=id1DAT(iData)
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
!			calcXsolute=YDAT(iData)/actCo( id1Dat(iData) ) !YDAT = ideal solubility.
!			if(LOUDER)print*,'calcX,exptX=',calcXsolute,exptXsolute
!			if(calcXsolute < 1.D-33 )then	! .or. calcXSolute > 1-zeroTol
!				ierCalc=2
!				print*,'calcX,exptX,idX,gam  ',calcXsolute,exptXsolute,YDAT(iData),actCo( id1Dat(iData) )
!				pause 'KijDevSLE: calcX < 0???'
!			endif
			if(ierCalc==0)then
				if(iOptimize < 3)actCo(1:NC)=exp( FUGCL(1:NC)-fugcPure(1:NC) )
				calcXsolute=YDAT(iData)/actCo(iSolute)	!for quick&dirty Kij optimization, incl Wilson&NRTL
				if(iOptimize==2)calcXsolute=X(iSolute)					!for optimization based on SLEFL
				if(iOptimize >2)itMax=1
				dev1=DLOG(calcXsolute/exptXsolute)
				xOld(iSolute)=calcXsolute ! bootstrap previous temperature.
				!dev1= DLOG( calcXsolute/(1-calcXsolute) )-DLOG( exptXsolute/(1-exptXsolute) ) ! Vladimir's metric, 
				dev1=dev1*100 ! convert to %
				deviate(iData) = dev1 ! use log deviate to minimize bias for bad fits.
				nConverged=nConverged+1 ! Two dev's in above deviate. Consistent with one-sided when coexistence unavailable.
				!print*,'got one! nConverged=',nConverged
				PAALDK = PAALDK+ABS(dev1) !+ABS(dev2)
				PBIASK = PBIASK+(dev1) !+dev2)
				ssqErr= ssqErr+deviate(iData)*deviate(iData)
				if(LOUDER)print*,'KijDevSLE: Calling fugi.i,xiE,xiC=',id1Dat(iData),exptXsolute,x(1)
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
					!deviate(iData+nPtsBipDat)=0 ! no penalty for last evaluation. just compute the deviations for converged.
				endif
			endif !last=1 .and. compute activities
			if(LAST==1)WRITE(661,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1,actCo(iSolute),itMax ! 
			if(LAST==1.and.LOUDER)WRITE( 6 ,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1 !LastCall=>writ
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
			if(LAST==1)WRITE(661,617)iSolute,T,P,log(YDAT(iData)),LOG(calcXsolute),LOG(exptXsolute),dev1,actCo(iSolute),itMax 
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
	!C  PURPOSE:  ITERATIVE IDAC EVALUATIONS
	!C  PROGRAMMED BY:  JRE 2/2023
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
	DOUBLEPRECISION X(NMX),xU(NMX)
	DOUBLEPRECISION ZFEED(NMX),kPas(NMX),KINIT(NMX)
	COMMON/eta/etaLo,etaUp,zLo,zUp
	common/EtaLiqPas/etaPasLo,etaPasUp,zPasLo,zPasUp
	!COMMON/KVALUES/kPas
	COMMON/FEED/ZFEED
	!WRITE(*,*)'ENTER DESIRED MAXIMUM FOR NUMBER OF ITERATIONS'
	!READ(*,*)MAXIT
	MAXIT=333
	WRITE(6,*)' '
	outFile=TRIM(masterDir)//'\output\LLIterOut.txt'
	open(601,file=outFile)
	INIT=1
	!open(61,file='TPXY.txt')
	!write(61,*)' '
	!close(61)
	initCall=1
	ANSWER='Y'
    initK=1
	INIT=3 !check d2G/dx1^2 on first call and search spinodal to determine initial guess.
    write(601,*)'T,P,Uof,x(1),xU(1),x(2),xU(2)'
	do while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
1		continue
		WRITE(6,*)'ENTER TEMPERATURE(K) AND PRESSURE(MPa)  '
		READ(5,*)T,P
        if(T < 1)goto 1		! must be a typo
		write(52,*)T,P
		if(NC==2.and.initCall==1)then
			initCall=0
			tFactor=1.05d0
			itMax=22
			Call LLConsolute(NC,T,P,tFactor,itMax,x1Con,Tcon,Gcon,dG_dx1,uncert,iErrCon) !Consolute T is max T for binodal.
			if(iErrCon > 10)then
				iErr=11
				print*,'The guessed temperature is above Tc (consolute T). No LLE here.'
				exit
			endif
			Try2=Tcon/1.01	  !Ignore warnings (iErrCon < 10). Try again starting from the last best estimate of Tcon. 
			tFactor=1.001d0
			itMax=22
			Call LLConsolute(NC,Try2,P,tFactor,itMax,x1Con,Tcon,Gcon,dG_dx1,uncert,iErrCon)
			write(dumpUnit,610)'FYI: upper Tcon,x1con,uncert=',Tcon,x1con,uncert
        elseif(NC > 2)then ! For NC=2, K, zFeed, and UOF will be computed from results.
            if(initK==1)then
				WRITE(6,*)'ENTER GUESSES FOR xUi/XLiS   '
				READ(5,*)(KINIT(I),I=1,NC)
				WRITE(52,*)'ENTER GUESSES FOR xUi/XLiS   '
				write(52,*)(KINIT(I),I=1,NC)
                initK=0
            endif
			WRITE(6,*)' '
			WRITE(6,*)'ENTER FEED COMPOSITIONS      '
			READ(5,*)(ZFEED(I),I=1,NC)
			WRITE(52,*)'ENTER TEMPERATURE(K) AND PRESSURE(MPa)  ',T,P
			WRITE(52,*)'ENTER FEED COMPOSITIONS      '
			write(52,*)(ZFEED(I),I=1,NC)
            X(1:NC)=( zFeed(1:NC)/( 1+(KINIT(1:NC)-1)*0.5d0 )  )		!EL2ed ~Eq. 10.23 assuming U/F=0.5
            sumx=sum(X(1:NC))
            X(1:NC)=X(1:NC)/sumx
            XU(1:NC)=X(1:NC)*KINIT(1:NC)
            sumx=sum(XU(1:NC))
            XU(1:NC)=XU(1:NC)/sumx
		endif

		!do 6 i=1,1000
		IFLAG=0
		ITMAX=MAXIT	! must reinitialize 
        kPas(1:NC)=kInit(1:NC)

		!    LLEFL(T,P,NC,Initialize,CALCK,ITMAX,zFeed,X,XU,UOF,iErrCode)
		CALL LLEFL(T,P,NC,INIT,kPas,ITMAX,zFeed,X,xU,UOF,iErrCode)
		if(iErrCode > 9)then
			call PrintErrFlash(iErrCode)
			write(dumpUnit,611)'LLIter: iErrFL=',iErrCode
			INIT=3
		else 
			etaLo=etaPasLo
			etaUp=etaPasUp
			zUp=zPasUp
			zLo=zPasLo
			call FlashOut(NC,X,xU,T,P,itmax,Uof) !takes z,eta from common eta
			write(601,'(F8.2,2f8.4,9E12.5)')T,P,Uof,x(1),xU(1),x(2),xU(2)
			INIT=1 
			kInit(1:NC)=kPas(1:NC) !Store these for bootstrapping (init=1). Be careful to stay below Tcon.
		endif

		!  NOTE:  FOR REPEAT LL CALCULATIONS IT IS ASSUMED THAT BOOTSTRAPPING
		!         IS DESIRED.  OTHERWISE, THE USER SHOULD ANSWER 'N' TO
		!         THE REPEAT QUESTION BELOW AND GO BACK THROUGH THE MAIN 
		!         PROGRAM BEFORE PROCEEDING.


		WRITE(6,*)'REPEAT LL CALCULATIONS? Y/N?'
		READ(*,'(A1)')ANSWER
	enddo !while(ANSWER.NE.'N'.AND.ANSWER.NE.'n')
	close(601)
182	FORMAT(2X,'COMPONENT IS ',A20,2X,'ID NO. IS ',1X,i5)
183	FORMAT(6E11.4)
659	FORMAT(1X,4(1X,F11.4))
610 format(1x,a,11E12.4)
611 format(1x,a,i8,11E12.4)
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
	if(DEBUG)inFile=TRIM(PGLinputDir)//TRIM(FN)
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
	open(617,file=outFile)
600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x, 'etaL',6x, 'etaV')
601	FORMAT(1X,1(F7.2,2X,F10.6,1X),4(F9.5,1X),i5)
1000  CONTINUE
	WRITE(6,*)'ENTER COMPOSITION   '
	READ(5,*)(X(I),I=1,NC)
	if(LOUD)write(*,600)
	write(617,600)
	MAXIT=1111
	zero=0

	ITMAX=MAXIT
	tMax=1111
	CALL BootBP(tMax,x,NC,INIT,pMax,ITMAX,y,ier,0)
!c	CALL ErrCheck(ier,iflag)
	if(LOUD)write(* ,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	write(617,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	if(ier(1).eq.0)init=1

!C      PAUSE

	if(LOUD)write(*,*)'Starting dew point curve.'
	write(617,*)'Dew points.'
	do i=1,nc
		Y(i)=X(i)
	enddo

      ITMAX=MAXIT
	pMax=111
	CALL BootDT(tMax,X,NC,INIT,pMax,ITMAX,Y,ier,0)
!c	CALL ErrCheck(ier,iflag)
	if(LOUD)write(* ,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	write(617,601)tMax,pMax,X(1),Y(1),etaL,etaV,ier(1)
	if(ier(1).eq.0)init=1

86	continue
	if(LOUD)pause 'These results are tabulated in TPXY.txt.'
	CLOSE(617)   
	RETURN
	END

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE PropTable(NC,iErrF)
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
		!CALL Fugi(tKelvin,P,X,NC,LIQ,fugcL,zFactor,ier)
		CALL FugiTP( tKelvin,P,X,NC,LIQ,rhoMol_cc,zFactor,aRes,fugcL,uRes,iErrF )
		rhoG_cc=rMwAvg*P/(zFactor*Rgas*tKelvin)
		if(iErrF > 0)then
			WRITE(*,601)'PropTable: error in FugiTP calculation. iErrF,T,P,x1=',iErrF,tKelvin,P,x(1)
		endif
601	format(1x,a,i4,8E12.4)

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
		if(LOUDER)write(* ,610)tKelvin,rhoG_cc,zFactor,aRes,uRes,cmprsblty,CvRes_R,CpRes_R !,DSONK
		WRITE(601,610)tKelvin,rhoG_cc,zFactor,aRes,uRes,cmprsblty,CvRes_R,CpRes_R, uRes-aRes
		!Check derivatives using FuVtot
		vTotCc=rMwAvg/rhoG_cc
		gmol(1)=1
		isZiter=0 !=>compute derivative props
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGCL,zFactor,aRes,uRes,iErr)
		pNew = zFactor*rGas*tKelvin/vTotCc
		if(LOUDER)write(* ,612)tKelvin,pNew,zFactor,aRes,uRes,cmprsblty,cvRes_R,CpRes_R
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
	character outFile*251 !,tabChar*1
	DIMENSION X(NMX),Y(NMX),x1(15),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data nx,x1/15,0.001,0.01,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.99,0.999/
!	tabChar=char(9) !ascii tab character from key code chart 1
    if(LOUD)then
	    if(nc.gt.2)pause 'error - this option is only for 2 compos'
    end if
	if(nc.gt.2)return
	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	open(618,file=outFile)
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
	write(618,600)
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
	write(618,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
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
		write(618,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
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
	write(618,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
86	continue
	print*,'OutFile=',TRIM(outFile)
	pause 'Check results in OutFile.'
	CLOSE(618)   
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
	!dimension ier(12)
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
		call FugiESD(tK,pMPa,xFrac,nComps,LIQ,FUGC,rho,Z,aRes,uRes,iErrF)
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
	WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in PGLinputDir\...)'
	WRITE(*,*)'FIRST LINE =NPTS,id 2ND LINE = TK, PMPA'
    !FN='helium.dat'
	READ(*,'(A)')FN
	!inFile=TRIM(masterDir)//'\VpData\'//TRIM(FN)
	inFile=TRIM(FN)
	inFile=TRIM(PGLinputDir)//'\'//TRIM(FN)
	open(51,file=inFile)
	outFile=TRIM(masterDir)//'\Output\KijDb.txt'
	if(DEBUG)outFile='c:\spead\calceos\Output\KijDb.txt'
	open(619,file=outFile)
	SipFile=TRIM(masterDir)//'\Input\SipSpead.txt'
	if(DEBUG)SipFile=TRIM(PGLinputDir)//'\SipSpead.txt'
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
		write(619,607)id(1),(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),PAADP,pBias,pErrMax,rmsErr,name(1)
		write(72,608)id(1),(parm(i),i=1,nParms),name(1)
        exit
	enddo !while(ier.eq.0) loop over data sets in db

	close(51)
	close(619)
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
	Implicit DoublePrecision(A-H,K,O-Z)
	CHARACTER*64 Input_File_Name
	DoublePrecision X(NMX),fugcL(NMX) !,ier(12)
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
		if(iter > 111)ierScf=2
		Call FugiTP( T,P,X,NC,1,rhoMol_cc,ZL,aRes,FugcL,uRes,iFlag )
		IF (IFLAG > 0)then
			WRITE(6,'(a,22I4)')'ScfIter: iErrFugi=',iFlag
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

	!for Fugi:
	!C  ier = 1 - AT LEAST ONE ERROR
	!C        2 - NOT USED
	!C        3 - NOT USED
	!C        4 - ERROR IN ALPHA CALCULATION, SQRT(ALPHA) OR ITERATIONS
	!C        5 - rho IS -VE
	!C        6 - TOO MANY Z ITERATIONS
	!C        7 - eta > etaMax
	!C		11 - goldenZ instead of real Z.
	if(iFlag==1 )write(*,*)'Error returned from FUGI.'
	if(iFlag==6 )write(*,*)'Error: too many iterations.'
	IF(iFlag==4 )WRITE(*,*)'ERROR ALL UPPER PHASE'
	IF(iFlag==5 )WRITE(*,*)'ERROR ALL LOWER PHASE'

	DHL=DHONKT
	DSL=DSONK
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

	SUBROUTINE ThermoCheck(nComps)
	USE GlobConst
	!USE BIPs ! includes nConverged,Kij(),KTij(), ...
	!USE VpDb ! for VpCoeffs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DIMENSION xFrac(NMX),FUGC(NMX) !,yFrac(NMX),yy(nmx),xx(nmx),
	integer LIQ						   
	!DOUBLE PRECISION parm(nParms),deviate(maxPts) !,temParm(nParms)
	character*77 errMsg(12) !,outFile
	!DIMENSION GammaCalc(nmx),GammaExp(nmx),pSatExp(nmx)
	!common/ppVpCoeffs/vpCoeffs(NMX,5)
	!dimension vpCoeffs(nmx,5)
	! DATA iOpt/1/		!JRE 20210324 iOpt=2 was previously to optmize Hij instead of Kij. I gave up on that.
	!common/MEM2parts/FA,FD,betadFA_dBeta,betadFD_dBeta,aAssocPas,uAssocPas,zAssocPas
	errMsg(1)='KijVal Error: Number of components must equal 2.'
	iOpt=1
	write(*,*)'Enter 1 for T,P option or 2 for T,eta option'
	!read(*,*)iOpt
10	continue !re-entry for repeat calculations.
	if(iOpt==1)then
		xFrac(1)=1
		if(nComps==1)then
			write(*,*)'Enter LIQ, T(K), P(MPa) (-ve LIQ to stop)'
			read(*,*)LIQ,tKelvin,pMPa
			if(LIQ < 0)return
		else
			write(*,*)'Enter LIQ,T(K), P(MPa),x(1:nComps-1)...(-ve LIQ to stop)'
			read(*,*)LIQ,tKelvin,pMPa,xFrac(1:nComps-1)
			xFrac(nComps)=sum(xFrac(1:nComps-1))
			if(LIQ < 0)return
		endif
	else
		pause 'ThermoCheck: OOPS! Under development. Please check again later. Sorry.'
		goto 10
		xFrac(1)=1
		if(nComps==1)then
			write(*,*)'Enter T(K), eta'
			read(*,*)tKelvin,eta
		else
			write(*,*)'Enter T(K), eta, x(1:nComps-1)'
			read(*,*)tKelvin,eta,xFrac(1:nComps-1)
			xFrac(nComps)=sum(xFrac(1:nComps-1))
		endif
	endif
	stepSize=5D-4
	if(iOpt==1)then
		CALL FugiTP( tKelvin,pMPa,xFrac,nComps,LIQ,rhoMol_cc,zFactor,aRes_RT,FUGC,uRes_RT,iErrCode )
		if(iErrCode > 0)then
			write(*,*)'ThermoCheck: first call to FUgiTP. iErrCode=',iErrCode
			goto 10
		endif
		write(*,*)' ThermoCheck: rho,aRes_RT',rhoMol_cc,aRes_RT
		!Check that Gdep = SUM(xi*lnPhi_i) 
		Gdep_RT=aRes_RT+zFactor-1 -LOG(zFactor)
		Gcalc=SUM( xFrac(1:nComps)*FUGC(1:nComps) )
		gDev=(Gcalc-Gdep_RT)/Gdep_RT*100
		write(*,form610)' ThermoCheck: Gdep_RT,Gcalc,gDev%',Gdep_RT,Gcalc,gDev
		if( ABS(gDev) > 1.D-2 )pause 'ThermoCheck: OOPS! Gdep_RT is inconsistent!'
		!check that Z = rho*d(Ares/RT)/dRho
		rhoPlus =rhoMol_cc*(1+stepSize)
		rhoMinus=rhoMol_cc*(1-stepSize)
		isZiter=1
		CALL FuVtot(isZiter,tKelvin,1/rhoPlus ,xFrac,nComps,FUGC,Z,aPlus ,uRes,iErrPlus)
		write(*,*)' ThermoCheck: rhoPlus,aPlus',rhoPlus,aPlus
		CALL FuVtot(isZiter,tKelvin,1/rhoMinus,xFrac,nComps,FUGC,Z,aMinus,uRes,iErrMinus)
		write(*,*)' ThermoCheck: rhoMinus,aMinus',rhoMinus,aMinus
		Zcalc=1+rhoMol_cc*(aPlus-aMinus)/(rhoPlus-rhoMinus)
		zDev=(Zcalc-zFactor)   ! don't divide by zFactor here because 
		!write(*,form610)' ThermoCheck: aRes_RT,aPlus,aMinus',aRes_RT,aPlus,aMinus
		write(*,form610)' ThermoCheck: Z,Zcalc,zDev%',zFactor,zCalc,zDev			! Setting tol=5.D-5 makes Z+/- 0.00005.
		if( ABS(zDev) > 5.D-5 )pause 'ThermoCheck: OOPS! zFactor is inconsistent!'	! Imprecision due to iterative nature of MEM2.
		!check that uRes_RT = -T*d(Ares/RT)/dT
		tPlus =tKelvin*(1+stepSize)
		tMinus=tKelvin*(1-stepSize)
		isZiter=0
		CALL FuVtot(isZiter,tPlus ,1/rhoMol_cc,xFrac,nComps,FUGC,Z,aPlus ,uRes,iErrPlus)
		CALL FuVtot(isZiter,tMinus,1/rhoMol_cc,xFrac,nComps,FUGC,Z,aMinus,uRes,iErrMinus)
		Ucalc= -tKelvin*(aPlus-aMinus)/(tPlus-tMinus)
		uDev=(Ucalc-uRes_RT)/uRes_RT*100
		write(*,'(a,2F11.6,3E12.5)')' ThermoCheck: uRes_RT,Ucalc,uDev%',uRes_RT,Ucalc,uDev
		if( ABS(uDev) > 1.D-2 )pause 'ThermoCheck: OOPS! uRes_RT is inconsistent!'
		goto 10
	else
		pause 'ThermoCheck: OOPS! Under development. Please check again later. Sorry.'
		goto 10
	endif
	return
	end


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE TXYP(NC)
	!C
	!C  PURPOSE:  COMPUTE TXY GIVEN P.
	!C  PROGRAMMED BY:  JRE 6/99
	!C
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character outFile*251 !,tabChar*1
	DIMENSION X(NMX),Y(NMX),x1(15),ier(20)
	COMMON/eta/etaL,etaV,ZL,ZV
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	data nx,x1/15,0.001,0.01,0.05,0.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.99,0.999/
	!data nx,x1/13,0.99,0.95,0.9,.8,.7,.6,.5,.4,.3,.2,.1,0.05,0.01/
	!tabChar=char(9) !ascii tab character
	if(nc.ne.2 .and. LOUD)pause 'error - this option is only for 2 compos'
	if(nc.ne.2)return
	outFile=TRIM(masterDir)//'\output\TPXY.txt'
	open(621,file=outFile)
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
	write(621,600)
	MAXIT=1111
	zero=0
    INIT=0

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
	write(621,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)
	if(ier(1).eq.0)init=1

	DO I=2,nx
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
		write(621,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
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
	write(621,601)T,P,X(1),Y(1),etaL,etaV,ier(1)
	if(iSpead)write(62,'(6(f11.4,a1),i5)')T,tabChar,P,tabChar,X(1),tabChar,Y(1),tabChar,etaL,tabChar,etaV,tabChar,ier(1)

86	continue
	print*,'OutFile=',TRIM(outFile)
	pause 'Check results in OutFile.'
	CLOSE(621)   
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
	write(*,*)'First tabulate general curve, then point by point input.'
	write(*,*)' Enter 0 for default (from TrMin-Tr=0.95) or 1 for hi,lo,inc'
	read(*,*)iAns
	if(iAns==0)then 
		write(*,'(a,f7.2,a)')' Tc=',TC(1),'K. Enter TrMin (e.g. 0.45)'
		read(*,*) TrMin
		nSteps=10
		delTr=(0.95d0-TrMin)/nSteps
		iTrMin=INT(TrMin*100/5)*5
		TLo=TrMin*Tc(1)
		delT=delTr*Tc(1)
	else
		write(*,*)' Enter TLo, THi, increment.'
		read(*,*)TLo,THi,delT
		nSteps=(THi-TLo)/delT+0.01
	endif
	outFile=TRIM(masterDir)//'\output\Binodal.txt'
	open(631,file=outFile)
	if(bTPT)write(6 ,'(a,3F9.3)')' G0,G1,e0 ',G0,G1,eta0
	if(bTPT)write(631,'(a,3F9.3)')'G0,G1,e0 ',G0,G1,eta0
	write(*,*)'  T(K)',tabChar,'   PsatMPa',tabChar,'rhoVap',tabChar,'rhoLiq(g/cc)',tabChar,'etaVap',tabChar,'etaLiq'
	write(631,*)'  T(K)',tabChar,'   PsatMPa',tabChar,'rhoVap',tabChar,'rhoLiq(g/cc)',tabChar,'etaVap',tabChar,'etaLiq'
	do iStep=1,nSteps+1
		tKelvin=TLo+delT*(iStep-1)
		call PsatEar(tKelvin,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode)
		if(ierCode > 10)then
			if(LOUDER)pause 'Psat returned error code'
			exit
		endif
		!call EqualArea(NC,gmol,tKelvin,Psat,vLCc,vVCc)
		write(* ,'(F8.2,a1,5(F9.6,a1),i3)')tKelvin,tabChar,pMPa,tabChar,rhoV*rMw(1),tabChar,rhoL*rMw(1),tabChar,rhoV*bVolCC_mol(1),tabChar,rhoL*bVolCC_mol(1),tabChar,ierCode
		write(631,'(F8.2,a1,5(F9.6,a1),i3)')tKelvin,tabChar,pMPa,tabChar,rhoV*rMw(1),tabChar,rhoL*rMw(1),tabChar,rhoV*bVolCC_mol(1),tabChar,rhoL*bVolCC_mol(1),tabChar,ierCode
	enddo
	tKelvin=0.96*Tc(1)
	do while(tKelvin.gt.0)
		call PsatEar(tKelvin,pMPa,chemPot,rhoL,rhoV,uSatL,uSatV,ierCode)
		!call EqualArea(NC,gmol,tKelvin,Psat,vLCc,vVCc)
		write(*,'(F8.2,1x,5F9.6,i3)')tKelvin,pMPa,rhoV*rMw(1),rhoL*rMw(1),rhoV*bVolCC_mol(1),rhoL*bVolCC_mol(1),ierCode
		if(pMPa > 1.D-3)then
			write(631,'(F8.2,a1,5(F9.6,a1),i3)')tKelvin,tabChar,pMPa,tabChar,rhoV*rMw(1),tabChar,rhoL*rMw(1),tabChar,rhoV*bVolCC_mol(1),tabChar,rhoL*bVolCC_mol(1),tabChar,ierCode
		else
			write(631,'(F8.2,a1,5(E11.4,a1),i3)')tKelvin,tabChar,pMPa,tabChar,rhoV*rMw(1),tabChar,rhoL*rMw(1),tabChar,rhoV*bVolCC_mol(1),tabChar,rhoL*bVolCC_mol(1),tabChar,ierCode
		endif
		write(*,*)'ENTER Temperature (K). (-ve to stop) '
		READ(*,*) tKelvin
	enddo
	rhoc=Pc(1)*rMw(1)/( Zc(1)*8.31434*Tc(1) )
	etac=rhoc/rMw(1)*bVolCc_Mol(1)
	write( * ,*)'Tc(K)',tabChar,'PcMPa',tabChar,'rhoc(g/cc)',tabChar,'Zc',tabChar,'etac'
	write(631,*)'Tc(K)',tabChar,'PcMPa',tabChar,'rhoc(g/cc)',tabChar,'Zc',tabChar,'etac'
	write( * ,'(F8.2,a1,5(F9.6,a1))')tc(1),tabChar,pc(1),tabChar,rhoc,tabChar,Zc(1),tabChar,etac
	write(631,'(F8.2,a1,5(F9.6,a1))')tc(1),tabChar,pc(1),tabChar,rhoc,tabChar,Zc(1),tabChar,etac
	close(631) 
	if(LOUDER)write(*,*) 'Your results are in ',TRIM(outFile)
	if(LOUDER)pause
86	return 
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine SortTable(table,maxRows,nRows,nCols,keyCol,LoToHi)
	! Purpose: Sort according to keyCol using NumRecp program indexx().
    ! LoToHi=1 if you want to sort from smallest(least positive) to largest.
	! 
	Implicit DoublePrecision(a-h,o-z)
	DoublePrecision	table(maxRows,nCols),table2(maxRows,nCols),sortVector(maxRows)
    Integer sortIndex(maxRows)
	!Integer iPerm(maxRows),NI(1),INDKEY(1),NRMISS
	!iAbsolute=0	! sort by algebraic values not absolute values
	if(LoToHi > 1 .or. LoToHi < 0)pause 'SortTable: LoToHi must be 0 or 1.'
	!iOrder=1-LoToHi	! IMSL swaps the order specifier.
	!iRet=2 ! sorted table is returned but not NGROUP or NI. (iRet=0 to return NGROUP and NI).
	!nKey=1 ! just sort by one column.
	!INDKEY(1)=keyCol
	!CALL DSROWR(nRows,nCols,table,maxRows,iAbsolute,IORDR,iRet,nKey,INDKEY,IPERM,NGROUP,NI,NRMISS)!IMSL routine
	!LDX — Leading dimension of table exactly as specified in the dimension statement of the calling program.   (Input)
	do i=1,nRows
		sortVector(i)=table(keyCol,i)
    enddo
    call indexx(nRows,sortVector,sortIndex)	! from NumRep
    do i=1,nRows
        do j=1,nCols
			table2(i,j)=table(sortIndex(i),j)
        enddo
    enddo
    table=table2	! copy back to original input.	 
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
! The rows of the sorted X are partitioned into groups. A group contains rows that are equal 
!	with respect to the method of comparison. NGROUP is the number of groups of different rows.
! NI — Vector of length NGROUP containing the number of rows in each group.   (Output, if IRET £ 1)
! The first NI(1) rows of the sorted X are group number 1. The next NI(2) rows of the sorted X are group number 2. ... 
!	The last NI(NGROUP) rows of the sorted X are group number NGROUP. 
!	If NGROUP is not known prior to the invocation of this routine, NROW(an upper bound for NGROUP) can be used as the 
!	dimension of NI. If IRET ³ 2, NI is not referenced and can be a vector of length one.
! NRMISS — Number of rows that contained NaN in the columns of X used in the sort.   (Output)
! These rows are considered as a separate group from the other NGROUP groups and are put as the last NRMISS rows of the sorted X.
! 
