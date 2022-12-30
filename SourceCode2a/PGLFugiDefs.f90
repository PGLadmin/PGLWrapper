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
	SUBROUTINE FUGI( T,P,X,NC,LIQ,FUGC,Z,ier )
	USE GlobConst
	USE comflash, only:ncomax ! PrLorraine's module for many constants. ncomax needed to dimension properly
	IMPLICIT DOUBLEPRECISION( A-H,O-Z )
	DoublePrecision X(*),fugc(*)
	DoublePrecision fg(ncomax) ,ft(ncomax) , fp(ncomax) , fx(ncomax,ncomax)	!for PrLorraine
    integer ier(12) !,ier2(11)
    ier(1)=0 !required if EOS uses iErr instead of ier().
	if(iEosOpt.eq.0)iEosOpt=1
	if(iEosOpt==1)then
		call FugiPR( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isESD)then ! iEosOpt=2, 
		call FugiESD( T,P,X,NC,LIQ,FUGC,Z,ier )
	!elseif(iEosOpt==4)then	!unused
	elseif(iEosOpt==3)then
	    call FugiPRWS( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(iEosOpt.eq.8)then
		call FuSpeadGamma( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isTPT)then ! 
 	    call FuTpt( T,P,X,NC,LIQ,FUGC,Z,ier )	  !AUG 10
	!elseif(iEosOpt.eq.6)then  ! unused
	elseif(iEosOpt.eq.7)then
		call FuNRTL( T,P,X,NC,LIQ,FUGC,Z,ier )
	elseif(isPcSaft)then !iEosOpt=11,13
		call FugiPcSaft(NC,T,P,X,LIQ,Z,FUGC,iErr)
		if( iErr.ne.0)ier(1)=1
 	    !if(LOUD)print*,'Sorry: PcSaft not available at this time'!call FugiPcSaft(T,P,X,NC,LIQ,FUGC,Z,ier)	 ! Need to define this
	elseif(iEosOpt.eq.11)then
		call FugiPrTc( T,P,X,NC,LIQ,FUGC,Z,ier )	  ! JRE 2020
	elseif(iEosOpt.eq.17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
		ntyp= -1 !vapor. Note: ntyp=0 returns the min fugacity phase.
		if(LIQ==1)ntyp=1
		iErrPRL=0
		call FuPrLorraine_TP(icon,ntyp,mtyp,iErrPRL,T,P,X,rho,Z, &	 ! JRE 20210510
		&                fg,ft,fp,fx,uRes_RT,aRes_RT,CvRes_R,CpRes_R,dpdRho_RT)
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
    integer ier(12)
	DIMENSION gmol(NC),FUGC(NC) !,dHkcalMol(NMX),KCSTARp(NMX),IER(12)
	DoublePrecision ft(3) , fv(3) , fx(3,3)	!for PrLorraine
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	COMMON/ETA2/ETA
	nComps=NC
	if(iEosOpt==1)then
		call FuPrVtot(isZiter,tKelvin,1/vTotCc,gmol,NC,LIQ,FUGC,Z,aDep,uDep,IER)
	!elseif(iEosOpt==4)then ! unused for now.
	elseif(isESD)then ! (iEosOpt==2.or.iEosOpt==6.or.iEosOpt==12) all assume ESD96. 
		call FuEsd96Vtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aDep,uDep,iErr)
	elseif(isTPT)then
		call FuTptVtot(isZiter,Z,aDep,uDep,FUGC,vTotCc,tKelvin,gmol,nComps,iErr)	  !AFG 2011
	elseif(iEosOpt==11)then
	    !SUBROUTINE FuPrVtot(tAbs,rhoMol_Cc,xFrac,NC,LIQ,FUGC,zFactor,aDep,uDep,IER)
		call FuPrTcVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aDep,uDep,iErr)
	elseif(isPcSaft)then
		!subroutine FugiPcSaftVtot(nComps,tKelvin,vTotCc,gmol,pMPa,zFactor,chemPoRes,iErr)
		call FugiPcSaftVtot(NC,tKelvin,vTotCc,gmol,P,Z,FUGC,iErr)
		if( iErr.ne.0)ier(1)=1
	elseif(iEosOpt==17)then
		icon=2 ! returns fugc and T,V derivatives. See GetPrLorraine for more options. 
        call FuPrLorraine_TV(icon,iErrPRL,tKelvin,vTotCc,gmol,P,Z, &
     &                 FUGC,ft,fv,fx,uRes_RT,aRes_RT,CvRes_R,CpRes_R,dpdRho_RT)
		if(iErrPRL.ne.0)ier(1)=1
	else
		if(LOUD)print*,'FuVtot: undefined iEosOpt=',iEosOpt
    endif
    if(ier(1).ne.0)iErr=1
	return
    end

	!***********************************************************************

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer Function ItsEven(iArg) ! returns 0 for odd and 1 for even.
    integer iArg
	ItsEven=1-MOD(iArg,2)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!int QueryNparPure();    //number of adjustable parameters for the current model and mixture
	integer Function QueryNparPure()
	USE GlobConst, only: iEosOpt !, class
	parameter(nModels=16)
!   cf. GlobConst for iEosOpt's as listed below
!    1     2      3      4      5        6         7         8          9         10      11     12       13         14          15           16
!	'PR','Esd96','PRWS','Esd','SPEAD','FloryWert','NRTL','SpeadGamma','SPEAD11','PcSaft','tcPR','GCESD','GcEsdTb','TffSPEAD','GcPcSaft','GcPcSaft(Tb)'/
	integer nParInert(nModels),nParAssoc(nModels),nParPolar(nModels) ! the # of pure compound parameters for the ith iEosOpt
	!               1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	data nParInert /0,3,0,3,11,0,0,11,0, 3, 4, 3,11, 0, 3, 3/	! e.g. m, sigma, eps/kB. for PcSaft
	data nParAssoc /0,2,0,2, 0,0,0, 0,0, 2, 0, 2, 0, 0, 2, 2/	! GC&Tff EOS's have zero adj par's even if assoc. JRE 20210421
	data nParPolar /0,0,0,0, 0,0,0, 0,0, 2, 0, 0, 0, 0, 0, 0/	! GC&Tff EOS's have zero adj par's even if assoc. JRE 20210421
	! PR & ESD's use Tc,Pc,w. AC models have no pure pars. 
	! SPEAD: 11=A0coeff(3),A1coeff(3),A2Coeff(4),vEffNm3. Does not allow adjustment of Assoc parameters because these must be transferable. 
	! tcPR has alphaN&M, Gc&Tff models have zero adj parameters 
	QueryNparPure=nParInert(iEosOpt)+nParPolar(iEosOpt)+nParAssoc(iEosOpt) 
	!if(class=='assoc')QueryNparPure=QueryNparPure+nParAssoc(iEosOpt)		! m, sigma, eps/kB, Kad, epsHB
	return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParPure(iComp,iParm,value,iErr)
	USE GlobConst      ! iEosOpt,
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

