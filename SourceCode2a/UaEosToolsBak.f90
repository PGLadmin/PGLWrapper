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
	SUBROUTINE FUGI(T,P,X,NC,LIQ,FUGC,Z,ier)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
	dimension X(nmx),fugc(nmx)
    integer ier(11) !,ier2(11)
	if(iEosOpt.eq.0)iEosOpt=1
	if(iEosOpt==2.or.iEosOpt==12)then !these are specific to ESD96.  
		call FugiESD(T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(iEosOpt==1)then
		call FugiPR(T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(iEosOpt==3)then
	    call FugiPRWS(T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(iEosOpt==4)then
        call FuEsdTp(T,P,X,NC,LIQ,FUGC,Z,ier)  !AUG 10
	elseif(iEosOpt.eq.8)then
		call FuSpeadGamma  (T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(isTPT)then
 	    call FuTpt(T,P,X,NC,LIQ,FUGC,Z,ier)	  !AUG 10
	elseif(iEosOpt.eq.6)then
	    call FuFloryWert(T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(iEosOpt.eq.7)then
		call FuNRTL(T,P,X,NC,LIQ,FUGC,Z,ier)
	elseif(iEosOpt.eq.10)then
 	  if(LOUD)print*,'Sorry: PcSaft not available at this time'!call FugiPcSaft(T,P,X,NC,LIQ,FUGC,Z,ier)	 ! Need to define this
	elseif(iEosOpt.eq.11)then
 	  call FugiPrTc( T,P,X,NC,LIQ,FUGC,Z,ier)	  ! JRE 2020
	else
	  if(LOUD)pause 'Error in Fugi: unexpected iEosOpt'
	  ier(1)=1
	endif
	return
	end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																		   C
!C Written by AFG, Oct. 2009											   C
!C It is a wrapper for vTot routines									   C
!C			  															   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         	          
	SUBROUTINE FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErr)
	USE GlobConst !iEosOpt
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
    integer ier(11)
	DIMENSION gmol(NMX),FUGC(NMX) !,dHkcalMol(NMX),KCSTARp(NMX),IER(12)
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	COMMON/ETA2/ETA
	nComps=NC
	!	if(iEosOpt==4.or.iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12)then
	if(isESD)then
		call FuEsdVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErr)
	!elseif(iEosOpt==5.or.iEosOpt==8.or.iEosOpt==9.or.iEosOpt==13)then
	elseif(isTPT)then
		call FuTptVtot(isZiter,Z,aDep,uDep,vTotCc,tKelvin,gmol,nComps,iErr)	  !AFG 2011
	elseif(iEosOpt==1)then
	    !SUBROUTINE FuPrVtot(tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
		call FuPrVtot(isZiter,tKelvin,1/vTotCc,gmol,NC,LIQ,FUGC,Z,aDep,uDep,IER)
	elseif(iEosOpt==11)then
	    !SUBROUTINE FuPrVtot(tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aDep,uDep,iErr)
	elseif(iEosOpt==10)then
		!call FuPcSaftVtot(tKelvin,1/vTotCc,gmol,NC,LIQ,FUGC,Z,aDep,uDep,iErr)
	else
		if(LOUD)print*,'FuVtot: undefined iEosOpt=',iEosOpt
    endif
    if(ier(1).ne.0)iErr=1
	return
    end

	!***********************************************************************

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer Function ItsEven(iArg)
    integer iArg
	ItsEven=1-MOD(iArg,2)
	return
	end

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

