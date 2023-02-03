	SUBROUTINE FUGI( T,P,X,NC,LIQ,FUGC,Z,iErr )
	! This is a dummy routine to patch legacy code through to newer FugiTP(). 
	USE GlobConst
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	!USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision X(*),FUGC(*)
	CALL       FugiTP( T,P,X,NC,LIQ,rhoMol_cc,Z,aRes,FUGC,uRes,iErrF )
	iErr=iErrF
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PURPOSE:	CALL THE APPROPRIATE ROUTINE FOR COMPUTING PROPERTIES GIVEN T,P.   
!	Module Use:
!		iEosOpt:	(GlobConst) indicates the EOS option (cf.Module DLLConst)
!	Input:
!		LIQ:		Phase/property indicator
!			2=>vapor  including zFactor and aRes only. 0=>vapor  also including uRes and FUGC. -2 vapor  with Z,U,A,FUGC plus CvRes, CpRes, (dP/dRho)/(RT) 	 	
!			3=>liquid including zFactor and aRes only. 1=>liquid also including uRes and FUGC. -1 liquid with Z,U,A,FUGC plus CvRes, CpRes, (dP/dRho)/(RT) 	 	
!                  	
	SUBROUTINE FugiTP( T,P,X,NC,LIQ,rhoMol_cc,Z,aRes,FUGC,uRes,iErr )
	USE GlobConst
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision X(*),FUGC(*)
	DoublePrecision fg(ncomax) ,ft(ncomax) , fp(ncomax) , fx(ncomax,ncomax)	!for PrLorraine
    if(LOUD)write(dumpUnit,*)'FUGI: NC,isZiter=',NC,isZiter 
    if(LOUD)write(dumpUnit,610)'FUGI: T,P,gmol=',tKelvin,P,X(1:NC)
610 format(1x,a,12E12.4)
	if(iEosOpt.eq.0)iEosOpt=1
	if(iEosOpt==1)then
		call FugiPR( T,P,X,NC,LIQ,FUGC,rhoMol_cc,Z,aRes,uRes,iErr )
	elseif(isESD)then ! iEosOpt=2,4,12,13 
		call FugiESD(T,P,X,NC,LIQ,FUGC,rhoMol_cc,Z,aRes,uRes,iErr)
	elseif(iEosOpt==3)then
	    call FugiPRWS( T,P,X,NC,LIQ,FUGC,Z,iErr )
	elseif(iEosOpt.eq.8)then
		call FuSpeadGamma( T,P,X,NC,LIQ,FUGC,Z,iErr )
	elseif(isTPT)then ! 
 	    call FuTpt( T,P,X,NC,LIQ,FUGC,rhoMol_cc,Z,aRes,uRes,iErr )	  !AUG 10
	elseif(iEosOpt.eq.7)then
		call FuNRTL( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isPcSaft)then !iEosOpt=11,13
		call FugiPcSaft(T,P,X,NC,LIQ,FUGC,rhoMol_cc,Z,aRes,uRes,iErr)
 	    !if(LOUD)print*,'Sorry: PcSaft not available at this time'!call FugiPcSaft(T,P,X,NC,LIQ,FUGC,Z,ier)	 ! Need to define this
	elseif(iEosOpt.eq.11)then
		call FugiPrTc( T,P,X,NC,LIQ,FUGC,rhoMol_cc,Z,aRes,uRes,iErr )	  ! JRE 2020
	elseif(iEosOpt.eq.17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
		ntyp= -1 !vapor. Note: ntyp=0 returns the min fugacity phase.
		if(LIQ==1)ntyp=1
		iErrPRL=0
		call FuPrLorraine_TP(icon,ntyp,mtyp,iErrPRL,T,P,X,rho,Z, &	 ! JRE 20210510
		&                fg,ft,fp,fx,uRes,aRes,CvRes_R,CpRes_R,dpdRho_RT)
		FUGC(1:NC)=fg(1:NC)
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
    if(LOUD)write(dumpUnit,610)'FUGI: Z,A,U,FUGC=',Z,aRes,uRes,FUGC(1:NC)
	return
	end	!SUBROUTINE FugiTP() 


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																		   C
!C Written by AFG, Oct. 2009											   C
!C It is a wrapper for vTot routines									   C
!C			  															   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         	          
	SUBROUTINE FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
	USE GlobConst !iEosOpt
	USE ESDParms ! for SetParPureEsd and QueryParPureEsd
	IMPLICIT DoublePrecision(A-H,K,O-Z)
	DoublePrecision gmol(NC),FUGC(NC) !,dHkcalMol(NMX),KCSTARp(NMX)
	DoublePrecision ft(3) , fv(3) , fx(3,3)	!for PrLorraine
    Integer initCall,iErr,isZiter
    data initCall/11/
    if(LOUD)write(dumpUnit,*)'FuVtot: NC,isZiter=',NC,isZiter 
    if(LOUD)write(dumpUnit,610)'FuVtot: T,V,gmol=',tKelvin,vTotCc,gmol(1:NC)
610 format(1x,a,12E12.4)
	nComps=NC
	if(iEosOpt==1)then
		call FuPrVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
	!elseif(iEosOpt==4)then ! unused for now.
	elseif(isESD)then ! (iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12) all assume ESD96. 
		call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
	elseif(isTPT)then
		call FuTptVtot(isZiter,Z,aRes,uRes,FUGC,vTotCc,tKelvin,gmol,nComps,iErr)	  !AFG 2011
	elseif(iEosOpt==11)then
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErr)
	elseif(isPcSaft)then
		call FugiPcSaftVtot(NC,tKelvin,vTotCc,gmol,P,Z,aRes,uRes,FUGC,iErr)
	elseif(iEosOpt==17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
        call FuPrLorraine_TV(icon,iErrPRL,tKelvin,vTotCc,gmol,P,Z, &
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
    if(LOUD)write(dumpUnit,610)'FuVtot: Z,A,U,FUGC=',Z,aRes,uRes,FUGC(1:NC)
	return
    end

	!***********************************************************************
