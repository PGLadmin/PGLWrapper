	!C
	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!C
	!C      PURPOSE --- THIS ROUTINE CALLS THE APPROPRIATE ROUTINE FOR COMPUTING
	!C                  FUGACITY COEFFS.  CHOICE IS DICTATED BY iEos FROM common-eosopt
	!C                  iEos = 1 = PR 1976 
	!C                  iEos = 2 = PR SV Wong-Sandler (aicheJ, 38:677, 92)
	!C                  iEos = 3 = ESD
	!C                  	
	!C
	!C      CODED BY:   JRE, 7/8/99
	!C
	!C
	!C
	!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!C
	!C
	SUBROUTINE FUGI( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,zFactor,ier )
	USE GlobConst
	!USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision Xfrac(*),fugc(*)
	!DoublePrecision fg(ncomax) ,ft(ncomax) , fp(ncomax) , fx(ncomax,ncomax)	!for PrLorraine
    integer ier(12) !,ier2(11)
    ier(1)=0 !required if EOS uses iErr instead of ier()
	Call FugiTP( tKelvin,pMPa,xFrac,nComps,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErrF )
	aRes_RT=aRes
	uRes_RT=uRes
	ier(1)=iErrF
	return
	end

	SUBROUTINE FugiTP( tKelvin,pMPa,xFrac,nComps,LIQ,rhoMol_cc,zFactor,aRes,FUGC,uRes,iErrCode )
	USE GlobConst
	USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision Xfrac(*),fugc(*)
	DoublePrecision fg(ncomax) ,ft(ncomax) , fp(ncomax) , fx(ncomax,ncomax)	!for PrLorraine
    integer ier(12) !,ier2(11)
    ier(1)=0 !required if EOS uses iErr instead of ier()
	iErrCode=0.
	if(iEosOpt==0)iEosOpt=1
!	if(iEosOpt==2)then !all other ESD is ESD96.  
!		call FugiESD96( T,P,X,NC,LIQ,FUGC,Z,ier )
	if(isESD)then !all other ESD is ESD96.  
		call FugiESD( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErrF )
	elseif(iEosOpt==1)then
		call FugiPR( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,ier )
	elseif(iEosOpt==3)then
	    call FugiPRWS( T,P,X,NC,LIQ,FUGC,Z,ier )
		aRes=Ares_RT
		uRes=Ures_RT
	elseif(iEosOpt.eq.8)then
		call FuSpeadGamma( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isTPT)then
 	    call FuTpt( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErrF )	  !AUG 10
	elseif(iEosOpt.eq.6)then
	    call FuFloryWert( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(iEosOpt.eq.7)then
		call FuNRTL( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isPcSaft)then 
		call FugiPcSaft(tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErrF)
		if( iErrF.ne.0)ier(1)=1
 	    !if(LOUD)print*,'Sorry: PcSaft not available at this time'!call FugiPcSaft(T,P,X,NC,LIQ,FUGC,Z,ier)	 ! Need to define this
	elseif(iEosOpt.eq.11)then
		call FugiPrTc( tKelvin,pMPa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,ier )	  ! JRE 2020
	elseif(iEosOpt.eq.17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
		ntyp= -1 !vapor. Note: ntyp=0 returns the min fugacity phase.
		if(LIQ==1)ntyp=1
		iErrPRL=0
		call FuPrLorraine_TP(icon,ntyp,mtyp,iErrPRL,tKelvin,pMPa,xFrac,rhoMol_cc,zFactor, &	 ! JRE 20210510
		&                fg,ft,fp,fx,uRes,aRes,CvRes_R,CpRes_R,dpdRho_RT)
		FUGC(1:NC)=fg(1:NC)
		!c     mtyp:     (o):      indicator for phase type; ic returns
		!c                       1:   if called with mt = 1/-1 and phase is found
		!c                       -1:  if called with mt = 1/-1 and phase is not found
		!c                       2:   if called with mt = 0 and min g phase is liquid
		!c                       -2:  if called with mt = 0 and min g phase is vapour
		ier=0
		if(iErrPRL > 20)ier(1)=1 ! indicates desired phase was not found.
		if(iErrPRL > 21)ier(2)=1 ! indicates desired phase was not found.
		if(iErrPRL > 20)ier(3)=1 ! indicates desired phase was not found.
		if(iErrPRL > 20)ier(5)=1 ! indicates desired phase was not found.
		if(iErrPRL > 20)ier(6)=1 ! indicates desired phase was not found.
		if(mtyp  ==  -1)ier(1)=1 ! if called with mt = 1/-1 and phase is not found
		if(mtyp  ==   3)ier(4)=1 ! real part of the imaginary solution
	else
	  if(LOUD)pause 'Error in Fugi: unexpected iEosOpt'
	  iErrCode=1
	endif
	iErrCode=iErrF+10*iErrCode
	return
	end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																		   C
!C Written by AFG, Oct. 2009											   C
!C It is a wrapper for vTot routines									   C
!C			  															   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         	          
	SUBROUTINE FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
	USE GlobConst !iEosOpt
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
    integer ier(12)
	DIMENSION gmol(NC),FUGC(NC) !,dHkcalMol(NMX),KCSTARp(NMX),IER(12)
	DoublePrecision ft(3) , fv(3) , fx(3,3)	!for PrLorraine
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	COMMON/ETA2/ETA
	nComps=NC
!	if(iEosOpt==2)then ! (iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12) all assume ESD96. 
!		call FuEsd96Vtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,Ares,Ures,iErr)
	if(isESD)then ! (iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12) all assume ESD96. 
		call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,Ares,Ures,iErr)
	elseif(isTPT)then
		call FuTptVtot(isZiter,Z,aRes,uRes,FUGC,vTotCc,tKelvin,gmol,nComps,iErr)	  !AFG 2011
	elseif(iEosOpt==1)then
	    !SUBROUTINE FuPrVtot(tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
		call FuPrVtot(isZiter,tKelvin,1/vTotCc,gmol,NC,LIQ,FUGC,Z,Ares,Ures,IER)
	elseif(iEosOpt==11)then
	    !SUBROUTINE FuPrVtot(tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,Ares,Ures,iErr)
	elseif(isPcSaft)then
		!subroutine FugiPcSaftVtot(nComps,tKelvin,vTotCc,gmol,pMPa,zFactor,chemPoRes,iErr)
		call FugiPcSaftVtot(NC,tKelvin,vTotCc,gmol,P,Z,aRes,uRes,FUGC,iErr)
		if( iErr.ne.0)ier(1)=1
	elseif(iEosOpt==17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
        call FuPrLorraine_TV(icon,iErrPRL,tKelvin,vTotCc,gmol,P,Z, &
     &                 FUGC,ft,fv,fx,Ures,Ares,CvRes_R,CpRes_R,dpdRho_RT)
		if(iErrPRL.ne.0)ier(1)=1
	else
		if(LOUD)print*,'FuVtot: undefined iEosOpt=',iEosOpt
    endif
    if(ier(1).ne.0)iErr=1
	return
    end

	!***********************************************************************

