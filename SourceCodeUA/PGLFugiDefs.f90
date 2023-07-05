!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE FugiTP( tKelvin,pMPa,xFrac,nComps,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErr )
	!	Purpose: Call appropriate fugacity routine for EOS specified by iEosOpt USEd in GlobConst.
	!	Input:
	!		tKelvin		Temperature in Kelvins.
	!		pMPa		Pressure in MPa.
	!		xFrac()		Array of mole fractions.
	!		nComps		Number of components in the mixture, maximum of NMX as specified in GlobConst.
	!		LIQ			1 if call is phase is liquid and uRes, FUGC are computed, 3 for liquid without uRes,FUGC; 0,2 for vapor.
	!	Output:
	!		rhoMol_cc	Molar density in mol/cm3.
	!		zFactor		Compressibility factor.
	!		aRes		Residual Helmholtz energy, aRes/RT.
	!		FUGC()		Residual chemical potential array, ln(phi), where phi is the fugacity coefficient.
	!		uRes		Residual internal energy, uRes/RT.
	!		iErr		Error indicator, (0 if no error)	
	USE GlobConst
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision xFrac(NMX),FUGC(NMX)
	DoublePrecision fg(ncomax) ,ft(ncomax) , fp(ncomax) , fx(ncomax,ncomax)	!for PrLorraine
	!DEC$ ATTRIBUTES DLLEXPORT::FugiTP
	!!MS$ ATTRIBUTES DLLEXPORT::FugiTP
    if(LOUD)write(dumpUnit,*)'FUGITP: NC,LIQ=',nComps,LIQ 
    if(LOUD)write(dumpUnit,610)'FUGITP: T,P,gmol=',tKelvin,pMPa,xFrac(1:nComps)
610 format(1x,a,12E12.4)
	if(iEosOpt.eq.0)iEosOpt=1
	if(iEosOpt==1)then
		call FugiPR( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )
	elseif(isESD)then ! iEosOpt=2,4,12,13 
		call FugiESD(tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr)
	elseif(iEosOpt==3)then
	    call FugiPRWS( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,zFactor,iErr )
	elseif(iEosOpt.eq.8)then
		call FuSpeadGamma( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,zFactor,iErr )
	elseif(isTPT)then ! 
 	    call FuTpt( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )	  !AUG 10
	elseif(iEosOpt.eq.7)then
		call FuNRTL( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )
	elseif(iEosOpt==19)then
		call FuLsgMem2( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )
	elseif(isPcSaft)then !iEosOpt=11,13
		call FugiPcSaft(tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr)
 	    !if(LOUD)print*,'Sorry: PcSaft not available at this time'!call FugiPcSaft(T,P,X,NC,LIQ,FUGC,zFactor,ier)	 ! Need to define this
	elseif(iEosOpt.eq.11)then
		call FugiPrTc( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr )	  ! JRE 2020
	elseif(iEosOpt.eq.17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
		ntyp= -1 !vapor. Note: ntyp=0 returns the min fugacity phase.
		if(LIQ==1)ntyp=1
		iErrPRL=0
		call FuPrLorraine_TP(icon,ntyp,mtyp,iErrPRL,tKelvin,pMPa,xFrac,rho,zFactor, &	 ! JRE 20210510
		&                fg,ft,fp,fx,uRes,aRes,CvRes_R,CpRes_R,dpdRho_RT)
		FUGC(1:nComps)=fg(1:nComps)
		!c     mtyp:     (o):      indicator for phase type; ic returns
		!c                       1:   if called with mt = 1/-1 and phase is found
		!c                       -1:  if called with mt = 1/-1 and phase is not found
		!c                       2:   if called with mt = 0 and min g phase is liquid
		!c                       -2:  if called with mt = 0 and min g phase is vapour
		iErr=iErrPRL*10
		if(mtyp  ==  -1)iErr=iErr+1    ! if called with mt = 1/-1 and phase is not found
		if(mtyp  ==   3)iErr=iErr+mTyp ! if called with mt = 1/-1 and phase is not found
	else
	  if(LOUD)pause 'Error in Fugi: unexpected iEosOpt'
	  iErr=10
	endif
    if(LOUD)write(dumpUnit,*)'FUGI: Done. iErr=',iErr 
    if(LOUD)write(dumpUnit,610)'FUGI: zFactor,A,U,FUGC=',zFactor,aRes,uRes,FUGC(1:nComps)
	return
	end	!SUBROUTINE FugiTP() 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE FuVtot(isZiter,tKelvin,vTotCc,gmol,nComps,FUGC,zFactor,aRes,uRes,iErr)
	!	Purpose: Call appropriate fugacity routine for EOS specified by iEosOpt USEd in GlobConst.
	!	Input:
	!		isZiter		1 if call is for a zFactor iteration, where uRes and FUGC are not computed, 0 otherwise.
	!		tKelvin		Temperature (Kelvin)
	!		vTotCc		Total (extensive) volume in cm3.
	!		gmol()		Array of mole numbers, could be mole fractions, but xFrac(i)=gmol(i)/SUM(gmol(1:nComps)) will be computed anyway.
	!		nComps		Number of components in the mixture, maximum of NMX as specified in GlobConst.
	!	Output:
	!		FUGC()		Residual chemical potential array, ln(phi), where phi is the fugacity coefficient.
	!		zFactor		Compressibility factor.
	!		aRes		Residual Helmholtz energy, aRes/RT.
	!		uRes		Residual internal energy, uRes/RT.
	!		iErr		Error indicator, (0 if no error)	
	USE GlobConst !iEosOpt
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(NMX),FUGC(NMX) !,dHkcalMol(NMX),KCSTARp(NMX)
	DoublePrecision ft(3) , fv(3) , fx(3,3)	!for PrLorraine
    Integer initCall,iErr,isZiter
    data initCall/11/
	!DEC$ ATTRIBUTES DLLEXPORT::FuVtot
    if(LOUD)write(dumpUnit,*)'FuVtot: nComps,isZiter=',nComps,isZiter 
    if(LOUD)write(dumpUnit,610)'FuVtot: T,V,gmol=',tKelvin,vTotCc,gmol(1:nComps)
610 format(1x,a,12E12.4)
	if(iEosOpt==1)then
		call FuPrVtot(isZiter,tKelvin,vTotCc,gmol,nComps,FUGC,zFactor,aRes,uRes,iErr)
	!elseif(iEosOpt==4)then ! unused for now.
	elseif(isESD)then ! (iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12) all assume ESD96. 
		call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,nComps,FUGC,zFactor,aRes,uRes,iErr)
	elseif(isTPT)then
		call FuTptVtot(isZiter,zFactor,aRes,uRes,FUGC,vTotCc,tKelvin,gmol,nComps,iErr)	  !AFG 2011
	elseif(iEosOpt==11)then
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,gmol,nComps,FUGC,zFactor,aRes,uRes,iErr)
	elseif(isPcSaft)then
		call FugiPcSaftVtot(nComps,tKelvin,vTotCc,gmol,pMPa,zFactor,aRes,uRes,FUGC,iErr)
	elseif(iEosOpt==17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
        call FuPrLorraine_TV(icon,iErrPRL,tKelvin,vTotCc,gmol,pMPa,zFactor, &
     &                 FUGC,ft,fv,fx,uRes,aRes,CvRes_R,CpRes_R,dpdRho_RT)
		if(iErrPRL.ne.0)iErr=iErrPRL+10
	else
		if(LOUD)print*,'FuVtot: undefined iEosOpt=',iEosOpt
    endif
    if(initCall>0)then !this is it keep DebugDLL.txt down to a reasonable size
        initCall=initCall-1 ! We can set initCall to finite value in data statement, then count down.
        if(LOUD)close(dumpUnit)
        LOUD=.FALSE.
    endif
    if(LOUD)write(dumpUnit,*)'FuVtot: Done. initCall,iErr=',initCall,iErr 
    if(LOUD)write(dumpUnit,610)'FuVtot: zFactor,A,U,FUGC=',zFactor,aRes,uRes,FUGC(1:nComps)
	return
    end

	!***********************************************************************
	SUBROUTINE FUGI( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,zFactor,iErr )
	! This is a dummy routine to patch legacy code through to newer FugiTP(). 
	USE GlobConst
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	!USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision xFrac(*),FUGC(*)
	CALL       FugiTP( tKelvin,pMPa,xFrac,nComps,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErrF )
	iErr=iErrF
	return
	end

