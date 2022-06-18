!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C Calculation of nonClassical behavior in the vicinity of critical point using White's global RG method.	 C
!C This routine calculates the Helmholtz free energy of density fluctuations in a recursive manner.			 C
!C Adam Bymaster et al. Ind. Eng. Chem. Res. 2008, 47, 6264-6274											 C
!C Jianguo Mi and Yiping Tang, AIChE J, 52: 342-353, 2006													 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE HelmRG(nComps,gmol,tKelvin,eta,bVolMix,aDepNC,zNC,PMpaNC,dP_detaNC,d2P_deta2NC,d3P_deta3NC)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(maxIter=6,maxRhoStep=450,ng=maxRhoStep-1) 
	DIMENSION aEnergy(0:maxIter,0:maxRhoStep),Kcoeff(maxIter),aEnergyRG(maxRhoStep)
	DIMENSION aLong(maxIter,0:maxRhoStep),aShort(maxIter,0:maxRhoStep)
	DIMENSION rhoMesh2(0:maxRhoStep),rhoMesh(0:maxRhoStep),etaSpline(maxRhoStep)
	DIMENSION bfnl(0:maxRhoStep),bfns(0:maxRhoStep),deltaA(maxRhoStep)
	DOUBLEPRECISION mShape(NMX),mShape2,lambda,lambda2 
	DIMENSION gmol(NMX),xFrac(NMX),aBipAD(NMX,NMX),aBipDA(NMX,NMX)
	COMMON/BIPs_SPEAD/aBipAD,aBipDA
	COMMON/HbParms/dHkcalMol(NMX),bondVolNm3Esd(NMX)
	COMMON/HbParms2/ND(NMX),NDS(NMX),NAS(NMX)
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!COMMON/Assoc/eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes),bondVolNm3T(nmx,maxTypes),nDegree(nmx,maxTypes),nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),idType(nmx,maxTypes),localType(maxTypesGlobal),idLocalType(maxTypes),nTypes(NMX),nTypesTot
	common/rg/iFlagRG,iFlagFit,mShape,cutOff(NMX),Phi(NMX),Zi(NMX)

	etaStep=etaMax/maxRhoStep
    if (eta > (etaMax/2.d0)) then
		aDepNC=0.d0
		zNC=0.d0
		PMpaNC=0.d0
		dP_drhoNC=0.d0
		d2P_drho2NC=0.d0
		d3P_drho3NC=0.d0
		return
	endif
!for now it works just for pure components
	totMoles=0.d0
	do i=1,nComps
		totMoles=totMoles+gmol(i)
    enddo
    if(LOUD)then
	    if(totMoles < 1e-22)pause 'HelmRG: totMoles=0'
    end if
	do i=1,nComps
		xFrac(i)=gmol(i)/totMoles
	enddo
	mShape2=mShape(1)*mShape(1)
	sig=(vMolecNm3(1)*6.d0/Pi/mShape(1))**(1.d0/3.d0)  ![=]nm
	sig2=sig*sig
	sig3=sig*sig*sig
 	cutOffL=cutOff(1)*sig !Since the cutOffL param was actually L/sig
	cutOffL2=cutOffL*cutOffL
	cutOffL3=cutOffL*cutoffL2
	eps=Zi(1)*kB  	![=]MPa.nm3
	alpha=16.d0*Pi*eps*sig3/9.d0 	 ![=]MPa.nm6
	zitta2=9.d0*sig2/7.d0 	 ![=]nm2
	const=6.d0/(Pi*mShape(1)*sig3)	 ![=]1/nm3
	vEff=vMolecNm3(1)
	betaR=kB*tKelvin	![=]MPa.nm3
	tTang=3.d0

	rhoStep=etaStep/vEff  !rhoStep[=]1/nm3
	do i=0,ng
		rhoMesh(i)=real(i)*rhoStep
		rhoMesh2(i)=rhoMesh(i)*rhoMesh(i)  
	enddo
	nIter=0  !Means the classical value of HelmHoltz free energy
	aEnergy(nIter,0)=0.d0
	do i=1,ng
		etah=rhoMesh(i)*vEff  
		call FuTptVtotRG(bVolMix,xFrac,etah,tKelvin,nComps,aDepRG,iErr)
		aEnergy(nIter,i)=(aDepRG+dlog(rhoMesh(i))-1.d0)*rhoMesh(i)*betaR   ![=]MPa for chain
	enddo

!Renormalaization Loop.
	do while (nIter.LT.maxIter)	 
		nIter=nIter+1
		lambda=tTang**(nIter-1)*cutOffL
		!lambda=tTang**(nIter)*cutOffL
		lambda2=lambda*lambda
		!lambda2=lambda*lambda/2.d0
		Vn=(lambda/2.d0)**3
		!Vn=lambda**3
		Kcoeff(nIter)=betaR/Vn  !(tTang**(3*nIter))/(Lcutoff3)	  ![=]MPa
		do i=1,ng
			aLong(nIter,i)=aEnergy(nIter-1,i)+alpha*rhoMesh2(i)*mShape2 	 ![=]MPa
			aShort(nIter,i)=aEnergy(nIter-1,i)+alpha*Phi(1)*rhoMesh2(i)*zitta2/lambda2*mShape2  ![=]MPa
			bfnl(i)=aLong(nIter,i)
			bfns(i)=aShort(nIter,i)
		enddo
		aEnergy(nIter,0)=0.d0
		do iRho=1,( (ng+1)/2-1 )!starting from 1 since integration from 0 to 0 would be zero!!
			xkni=1.d0/Kcoeff(nIter) 
			call cgPure(bfnl,bfns,di,xkni,etaStep,iRho,ng) !,Dn,rhoMesh,aAttMesh)
			deltaA(iRho)=-Kcoeff(nIter)*di
			aEnergy(nIter,iRho)=aEnergy(nIter-1,iRho)+deltaA(iRho)
		enddo
		do iRho=((ng+1)/2),ng
			aEnergy(nIter,iRho)=aEnergy(nIter-1,iRho)
		enddo	
	enddo

!The nonClassical part of Helmholtz free energy stored as aEnergy0
	aEnergyRG(1)=0.d0
	etaSpline(1)=0.d0
	do i=2,(ng+1)/2
!aEnergyRG matrix arrays starts from 1 since we pass this matrix to HelmSpline routine
		aEnergyRG(i)=aEnergy(maxIter,(i-1))-aEnergy(0,(i-1)) 
		aEnergyRG(i)=(aEnergyRG(i))/rhoMesh(i-1)/betaR  		
		etaSpline(i)=rhoMesh(i-1)*vEff 
	enddo
	call helmSpline(eta,etaStep,etaSpline,tKelvin,bVolMix,aEnergyRG,(ng+1)/2,aDepNC,zNC,PMpaNC,dP_detaNC,d2P_deta2NC,d3P_deta3NC)
	
602	format(1x,f5.3,2x,f12.7,2x,f12.7,2x,f10.7)

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by Leo Lue																					 C
!C Purpose: numerical integration using trapezoidal rule													 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgPure(bfnl,bfns,di,xkni,dy,imax,ng) 
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension bfnl(0:ng),bfns(0:ng)
      dimension argl(0:ng),args(0:ng)
! *** search for maximum of gnl and ns ***
      aminl=0.d0
      amins=0.d0
      do i=0,imax
        gnl=0.5d0*(bfnl(imax+i)+bfnl(imax-i))-bfnl(imax) 
        gns=0.5d0*(bfns(imax+i)+bfns(imax-i))-bfns(imax) 
        argl(i)=gnl*xkni
        args(i)=gns*xkni
        if (argl(i).lt.aminl) aminl=argl(i)
        if (args(i).lt.amins) amins=args(i)
      end do

! *** numerical integration ***
      qnl=0.5d0*(exp(-argl(0)+aminl)+exp(-argl(imax)+aminl))
      qns=0.5d0*(exp(-args(0)+amins)+exp(-args(imax)+amins))
      do i=1,imax-1
         al=argl(i)-aminl
         as=args(i)-amins
         if (al.lt.30.d0) qnl=qnl+1.d0*exp(-al)
         if (as.lt.30.d0) qns=qns+1.d0*exp(-as)
      end do
      qnl=log(qnl*dy)-aminl	 
      qns=log(qns*dy)-amins	
      di=qns-qnl
	  
      return
      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C The HelmRG subroutine results in a correction term for A at 450 eta values.								 C
!C This routine interpolate between these values for A and it derivatives needed in the critical solver		 C
!C Spline fit is used for interpolation.																	 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE helmSpline(eta,etaStep,etaSpline,tKelvin,bVolMix,fn,ns,aDepNC,zNC,PMpaNC,dP_detaNC,d2P_deta2NC,d3P_deta3NC)
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(NMX=55,AvoNum=602.22d0,rGAs=8.314339637756d0,Pi=3.14159265359796d0)
	PARAMETER(kB=rGas/AvoNum,maxRhoStep=450) 
	DIMENSION etaSpline(maxRhoStep) 
	DIMENSION fn(maxRhoStep),fn1(maxRhoStep),fn11(maxRhoStep),zn(maxRhoStep),zn11(maxRhoStep)
	DIMENSION Pn(maxRhoStep),Pn1(maxRhoStep),Pn11(maxRhoStep),Pn111(maxRhoStep),Pn1111(maxRhoStep),Pn11111(maxRhoStep)

!NonClassical part of A:
	call spline(etaSpline,fn,ns,1.d30,1.d30,fn11)
	call splint(etaSpline,fn,fn11,ns,eta,aDepNC,iErrSplint)

!NonClassical part of Z:
	do i=1,ns-1
		fn1(i)=(fn(i+1)-fn(i))/etaStep-etaStep*(fn11(i)/3.d0+fn11(i+1)/6.d0)
		zn(i)=etaSpline(i)*fn1(i) !+1
	enddo
	i=ns
	fn1(i)=(fn(i)-fn(i-1))/etaStep-etaStep*(fn11(i-1)/3.d0+fn11(i)/6.d0)
	zn(i)=etaSpline(i)*fn1(i) !+1
	call spline(etaSpline,zn,ns,1.d30,1.d30,zn11)
	call splint(etaSpline,zn,zn11,ns,eta,zNC,iErrSplint)

!NonClassical part of P:
	do j=1,ns
		Pn(j)=zn(j)*etaSpline(j)*rGas*tKelvin/bVolMix
	enddo
	call spline(etaSpline,Pn,ns,1.d30,1.d30,Pn11)
	call splint(etaSpline,Pn,Pn11,ns,eta,PMpaNC,iErrSplint)

!NonClassical part of dP/drho:
	do i=1,ns-1
		Pn1(i)=(Pn(i+1)-Pn(i))/etaStep-etaStep*(Pn11(i)/3.d0+Pn11(i+1)/6.d0) !+rGas*tKelvin/bVolMix
	enddo
	i=ns
	Pn1(i)=(Pn(i)-Pn(i-1))/etaStep-etaStep*(Pn11(i-1)/3.d0+Pn11(i)/6.d0) !+rGas*tKelvin/bVolMix
	call spline(etaSpline,Pn1,ns,1.d30,1.d30,Pn111)
	call splint(etaSpline,Pn1,Pn111,ns,eta,dP_detaNC,iErrSplint)

!NonClassical part of d2P/drho2:
	call spline(etaSpline,Pn11,ns,1.d30,1.d30,Pn1111)
	call splint(etaSpline,Pn11,Pn1111,ns,eta,d2P_deta2NC,iErrSplint)

!NonClassical part of d3P/drho3:
	call spline(etaSpline,Pn111,ns,1.d30,1.d30,Pn11111)
	call splint(etaSpline,Pn111,Pn11111,ns,eta,d3P_deta3NC,iErrSplint)

    RETURN 
    END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C Purpose: Fitting the White's parameter for a pure component												 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine FitWhite(NC)
    USE GlobConst
    implicit doubleprecision (A-H,O-Z)
	parameter(nParms=3,nData=3)
	dimension parm(nParms),fvec(nParms),stdErr(nParms) 
	external DevalCrit
	character outFile*80
    character*100 ErrMsg(4)

	ErrMsg(1)='FitWhite Error: NC.ne.1'
	ErrMsg(2)='FitWhite Error: fu option just works with EOS # 9 for now'
	if (iEosOpt.ne.9) then
		write (*,*)ErrMsg(2)
		if(LOUD)pause
		goto 10
	endif
	if (NC.ne.1) then
		write (*,*)ErrMsg(1)
		if(LOUD)pause
		goto 10
	endif
	outFile='c:\MYPROJEX\CalcEos\output\fitWhite.txt'
	open(61,file=outFile)

!Initial guesses for White's RG method parameters
	parm(3)=5.d0	 !Lcutoff	
	parm(2)=550.d0  !eps
	parm(1)=10.d0	!Phi
	tol=1e-5
	factor=0.01D0
	call LmDifEz(DevalCrit,nData,nParms,parm,factor,FVEC,TOL,iErrCode,stdErr)
	write(61,*) ' L/sig  ','   Phi  ','   Zi' 
	write(61,'(3(2x,f7.3))')parm(3),parm(1),parm(2) 
	close(61)

	if(LOUD)pause 'Your results stored in output/fitWhite.txt'
10	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine DevalCrit(nData,nParms,parm,deviate,iFlag)
	implicit doubleprecision (A-H,O-Z)
	parameter(NMX=55,rGAs=8.314339637756d0,AvoNum=602.22d0)
	dimension deviate(nData),parm(nParms)
	doubleprecision mShape(NMX)
	COMMON/ppData/TC(NMX),PC(NMX),ACEN(NMX),ID(NMX)
	COMMON/ppDataPlus/ZC(NMX) !,rMw(nmx)
	COMMON/ETA2/ETA
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	common/rg/iFlagRG,iFlagFit,mShape,cutOff(NMX),Phi(NMX),Zi(NMX)
	iFlag=iFlagRG !just to kill the warning
	NC=1
	gmol=1.d0
	isZiter=0
	cutOff(1)=parm(3)
	phi(1)=parm(1)
	Zi(1)=parm(2)
	toll=1.d-7
	call CritPure(NC,isZiter,toll,TC_Pure,VC_PURE,PC_Pure,ZC_Pure,Acen_pure,iErrCode)
	etaC_PURE=eta
	bVolPure=vMolecNm3(1)*AvoNum
	etaCExp=PC(1)/ZC(1)/rGas/TC(1)*bVolPure
	deviate(1)=dsqrt( ( (TC_Pure-TC(1))/TC(1) )**2 ) 
	deviate(2)=dsqrt( ( (PC_Pure-PC(1))/PC(1) )**2 )*0.1d0	   !The weights are arbitrary
	deviate(3)=dsqrt( ( (etaC_Pure-etaCExp)/etaCExp )**2 )*0.2d0 

	return
	end 

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=1000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  SUBROUTINE splint(xa,ya,y2a,n,x,y,iErr)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
	  iErr=0
      klo=1
      khi=n
1     if (khi-klo > 1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if(h==0)then
		iErr=1
		return
		!pause 'bad xa input in splint'
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      END


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C A version of FuTptVtot special for helmRG routine...														 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine FuTptVtotRG(bVolMix,xFrac,eta,tKelvin,nComps,aDep,iErr)
	USE Assoc !GlobConst+XA,XD,XC
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	character*77 errMsg(11)
	DIMENSION xFrac(NMX)
    DoublePrecision FUGASSOC(NMX),dFUGASSOC_dT(NMX),dFUGASSOC_dRHO(NMX)
    COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/DEPFUN/DUONKT,DAONKT,DSONK,DHONKT
	COMMON/rdf/d2lng,d2g,dlng,dg_deta,dAlph_deta

	iErr=0
	errMsg(1)=' FuTptVtot Error: a1i or a2i > 0'
	errMsg(2)=' FuTptVtot Error: Wertheim returned error.'
	errMsg(3)=' FuTptVtot Error: Error from MixRule. Terminal, sorry.'
	if( (1-eta).lt.1.e-11)then
		write(*,*) 'Error in FuTptVtot: eta > 1'
		if(LOUD)pause
	endif
	call MixRuleRG(xFrac,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,iErrMix)
	if(iErrMix.ne.0)then
		iErr=3
		call BeepMsg(errMsg(iErr))
		return
	endif 
	aAtt=( a1Mix+ a2Mix/tKelvin)/tKelvin 
	if(iErrMix.eq.0)call Wertheim(vMolecNm3,eta,tKelvin,xFrac,nComps,zAssoc,aAssoc,uAssoc,iErrCode) 
	if (zAssoc.NE.0) then
		call WertheimFugc(xFrac,vMolecNm3,tKelvin,nComps,ETA,FUGASSOC,h_nMichelsen,dFUGASSOC_dT,dFUGASSOC_dRHO,dh_dT,dh_dRHO)
		iFlagAssoc=1
	else
		iFlagAssoc=0
	endif
 
	if(iErr.ne.0)then
		write(*,*)errMsg(iErr)
		return
	endif 	
	aDep=a0Mix+aAtt+aAssoc 

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine MixRuleRG(xFrac,eta,nComps,bVolMix,a0Mix,a1Mix,a2Mix,iErr) 
	USE GlobConst, only: avonum,LOUD
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	character*77 errMsg(11)
	DIMENSION xFrac(NMX)
	DIMENSION solParEntro(NMX),solParEnerg(NMX),vMolec(NMX) !,zRefMix(5), !a1Mix(5),a2Mix(5),
	DIMENSION ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	DIMENSION KII(NMX,NMX),KTII(NMX,NMX)
	COMMON/SIPs/KII,KTII
	COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	data initial/1/
	data ks0,ks1/-0.04, 0/		! best I could do based on avg of the entire system AV 11/29/08
	
	iErr=0
	errMsg(1)=' MixRule Error: a1i or a2i > 0'
	errMsg(2)=' MixRule Error: Wertheim returned error.'
	errMsg(3)=' MixRule Error: sqArg = a1i*a1j or a2i*a2j < 0'
	errMsg(4)=' MixRule Error: rMwMix .le. 0.  Check input.'
	errMsg(5)=' MixRule Error: eta out of range.'
	if(eta.ge.1.0 .or. eta.lt.0)then
		iErr=5
		write(*,*)'In TptTerms: eta=',eta
		call BeepMsg(errMsg(iErr))
		return
	endif
	do i=1,nComps
		vMolec(i)=vMolecNm3(i)*AvoNum
	enddo
	if(initial)then
		initial=0
		do iComp=1,nComps
			do jComp=1,nComps
				ks0ij(iComp,jComp)=0
				ks1ij(iComp,jComp)=0
				if(iComp.ne.jComp)ks0ij(iComp,jComp)=ks0
				if(iComp.ne.jComp)ks1ij(iComp,jComp)=ks1
			enddo
		enddo
		etaSto=eta
		etaStd=0.4
		do iComp=1,nComps
			call TptTermsRG(iComp,etaStd,a0i,a1i,a2i,iErrTpt)
			if (iErrTpt.ne.0) then
				iErr=6
				return
			endif
			solParEntro(iComp)=SQRT( a0i*1.987d0*298.d0*etaStd/vMolec(iComp) )
			solParEnerg(iComp)=SQRT(-a1i*1.987d0*etaStd/vMolec(iComp) )
		enddo
		eta=etaSto
	endif
	a0Mix=0
	a1Mix=0
	a2Mix=0
	do iComp=1,nComps
		bVoli=vMolec(iComp)
		call TptTermsRG(iComp,eta,a0i,a1i,a2i,iErrTpt)
		if (iErrTpt.ne.0) then
			iErr=6
			return
		endif
		if(a0i.le.0 .or. a1i.ge.0 .or. a2i.ge.0)then
			write(*,*)'ZCalcTpt error:a0i<0 or a1i or a2i .ge. 0. a0i,a1i,a2i'
			write(*,*)'   eta         a0i        a1i         a2i'
			write(*,'(1x,f8.4,3e13.5)')eta,a0i,a1i,a2i
			write(*,'(a,5e13.5)')' a1Coeff',(a1Coeff(iComp,i),i=1,nTptCoeffs)
			write(*,'(a,5e13.5)')' a2Coeff',(a2Coeff(iComp,i),i=1,nTptCoeffs)
			if(LOUD)pause 'Check your input.'
			iErr=1
			exit !terminate the do loop, Note: "cycle" would go to enddo and keep looping.
		endif
		if(iErr.ne.0)then
			if(LOUD)pause
			return
		endif

		do jComp=1,nComps
			bVolj=vMolec(jComp)
			call TptTermsRG(jComp,eta,a0j,a1j,a2j,iErr)
			if (iErrTpt.ne.0) then
				iErr=6
				return
			endif
		   	vij=(vMolecNm3(iComp)+vMolecNm3(jComp))/2.d0
			chiS=vij*602.22d0/1.987d0/298.d0*(solParEntro(iComp)-solParEntro(jComp))**2
			ksij=(ks0ij(iComp,jComp)+eta*ks1ij(iComp,jComp))*chiS
			aRefij=SQRT( bVoli*bVolj*a0i*a0j )*(1.d0-ksij)/bVolMix  !NOTE: genlzd value for ksij is -0.04
			a0Mix=a0Mix+xFrac(iComp)*xFrac(jComp)*aRefij
			if(nComps==1)then
				bipKij=KII(iComp,jComp)+KTII(iComp,jComp)*eta
			else
				bipKij=KIJ(iComp,jComp)+KTIJ(iComp,jComp)*eta
			endif
			if(iComp==jComp.and.nComps==2)bipKij=0
			if(a1j.ge.0 .or. a2j.ge.0)then
				write(*,*)'VtotTpt error:a1j or a2j .ge. 0. a1j,a2j',a1j,a2j
				if(LOUD)pause 'Check your input.'
				iErr=1
				exit
			endif
			a1ijSq=( bVoli*bVolj * a1i*a1j )
			if(a1ijSq.lt.0)then
				if(LOUD)pause 'VtotTpt Error:sqArg < 0 for a1ijSq'
				iErr=3
				exit
			else
				a1ij=-SQRT(a1ijSq)*(1-bipKij)/bVolMix
			endif
			a1Mix=a1Mix+xFrac(iComp)*xFrac(jComp)*a1ij
			a2ijSq=( bVoli*bVolj * a2i*a2j )
			if(a2ijSq.lt.0)then
				if(LOUD)pause 'ZCalcTpt Error:sqArg < 0 for a2ijSq'
				iErr=3
				exit
			else
				a2ij=-SQRT(a2ijSq)*(1-bipKij)/bVolMix
			endif
			a2Mix=a2Mix+xFrac(iComp)*xFrac(jComp)*a2ij
		enddo  
	enddo 

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	subroutine TptTermsRG(iComp,eta,a0i,a1i,a2i,iErr)
	USE BIPs
	implicit doublePrecision(A-H,K,O-Z)
	character*120 errMsg(4)
	common/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	!common/Assoc/eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes),bondVolNm3T(nmx,maxTypes),nDegree(nmx,maxTypes),nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),idType(nmx,maxTypes),localType(maxTypesGlobal),idLocalType(maxTypes),nTypes(NMX),nTypesTot

	errMsg(1)='TptTermsRG Error: eta out of range.'
	errMsg(2)='TptTermsRG Error: a1 or a2 > 0  '
	iErr=0
	if(eta.ge.1.0 .or. eta.lt.0)then
		write(*,*)'In TptTerms: eta=',eta
		write(*,*)errMsg(1)
		iErr=1
		return
	endif
	if(a1i.gt.0 .or. a2i.gt.0)then
		write(*,*)errMsg(2)
		iErr=2
	endif
	c1=zRefCoeff(iComp,1)+3
	c2=zRefCoeff(iComp,2)-3
	c3=zRefCoeff(iComp,3)+1
	void=1.d0-eta
	void2=void*void
	void3=void*void2
	eta2=eta*eta
	eta3=eta2*eta
	eta4=eta3*eta
	a0i=0.5d0*(c1+c2+c3)/void2-(c2+2*c3)/void-c3*LOG(void)-0.5d0*(c1-c2-3*c3)
	const=2.d0
	ePart1=exp(-const*eta3)
	ePart2=1.d0/(0.2d0+eta) 
	a1i=a1Coeff(iComp,1)*eta+a1Coeff(iComp,2)*eta*ePart1+a1Coeff(iComp,3)*eta3*ePart2+a1Coeff(iComp,4)
	bipKij=KIJ(iComp,iComp)+KTIJ(iComp,iComp)*eta
	denom=(1.d0+50.d0*eta3)
	a2i=(a2Coeff(iComp,1)*eta+a2Coeff(iComp,2)*eta2+a2Coeff(iComp,3)*eta3+a2Coeff(iComp,4)*eta4)/denom
	a1i=a1i*(1.d0-bipKij)	
	a2i=a2i*(1.d0-bipKij)
	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C This routine calculates the Helmholtz free energy using softSAFT EOS. It was used to validate the		 C
!C HelmRG subroutine and it is not used anywhere now in the project											 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    subroutine AsoftSaftRG(tKelvin,rho,aSS,m,sig,eps)

	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(AvoNum=602.22d0,rGAs=8.314472d0,Pi=3.14159265359796d0,kB=rGas/AvoNum)
	DIMENSION a(5,5),x(32),at(8),bt(6),gt(6)
	double precision m

	eok=eps/kB
	a(1,1)=0.49304346593882d0
	a(2,1)=-0.47031983115362d0
	a(3,1)=5.0325486243620d0
	a(4,1)=-7.3633150434385d0
	a(5,1)=2.9043607296043d0
	a(1,2)=2.1528349894745d0
	a(1,3)=-15.955682329017d0
	a(1,4)=24.035999666294d0
	a(1,5)=-8.6437958513990d0
	a(2,2)=1.1471647487376d0
	a(2,3)=37.889828024211d0
	a(2,4)=-84.667121491179d0
	a(2,5)=39.643914108411d0
	a(3,2)=-25.915399226419d0
	a(3,3)=-18.862251310090d0
	a(3,4)=107.63707381726d0
	a(3,5)=-66.602649735720d0
	a(4,2)=51.553565337453d0
	a(4,3)=-40.519369256098d0
	a(4,4)=-38.796692647218d0
	a(4,5)=44.605139198318d0
	a(5,2)=-24.418812869291d0
	a(5,3)= 31.500186765040d0
	a(5,4)=-5.3368920371407d0
	a(5,5)=-9.5183440180133d0
	T=tKelvin
	Ts=T/eok
	Ts2=Ts*Ts
	x(1)=0.8623085097507421d0
	x(2)=2.976218765822098d0
	x(3)=-8.402230115796038d0
	x(4)=0.1054136629203555d0
	x(5)=-0.8564583828174598d0
	x(6)=1.582759470107601d0
	x(7)=0.7639421948305453d0
	x(8)=1.753173414312048d0
	x(9)=2798.291772190376d0
	x(10)=-0.048394220260857657d0
	x(11)=0.9963265197721935d0
	x(12)=-36.98000291272493d0
	x(13)=20.84012299434647d0
	x(14)=83.054021244717285d0
	x(15)=-957.4799715203068d0
	x(16)=-147.7746229234994d0
	x(17)=63.98607852471505d0
	x(18)=16.03993673294834d0
	x(19)=68.05916615864377d0
	x(20)=-2791.29358794945d0
	x(21)=-6.245138304568454d0
	x(22)=-8116.836104958410d0
	x(23)=14.88735559561229d0
	x(24)=-10593.46754655084d0
	x(25)=-113.1607632802822d0
	x(26)=-8867.771540418822d0
	x(27)=-39.86982844450543d0
	x(28)=-4689.27299917261d0
	x(29)=259.3535277438717d0
	x(30)=-2694.523589434903d0
	x(31)=-721.8487631550215d0
	x(32)=172.1802063863269d0
	at(1)=x(1)*Ts+x(2)*sqrt(Ts)+x(3)+x(4)/Ts+x(5)/Ts2
	at(2)=x(6)*Ts+x(7)+x(8)/Ts+x(9)/Ts2
	at(3)=x(10)*Ts+x(11)+x(12)/Ts
	at(4)=x(13)
	at(5)=x(14)/Ts+x(15)/Ts2
	at(6)=x(16)/Ts
	at(7)=x(17)/Ts+x(18)/Ts2
	at(8)=x(19)/Ts2
	bt(1)=x(20)/Ts2+x(21)/(Ts2*Ts)
	bt(2)=x(22)/Ts2+x(23)/(Ts2*Ts2)
	bt(3)=x(24)/Ts2+x(25)/(Ts2*Ts)
	bt(4)=x(26)/Ts2+x(27)/(Ts2*Ts2)
	bt(5)=x(28)/Ts2+x(29)/(Ts2*Ts)
	bt(6)=x(30)/Ts2+x(31)/(Ts2*Ts)+x(32)/(Ts2*Ts2)
	gama=3.d0

    rhos=rho*(sig*sig*sig)   !monomer density
    g=0.d0
	dg=0.d0
    do i=1,5
        do j=1,5
            dg=dg+real(i)*a(i,j)*(rhos**(i+1))*Ts**(1-j)
			g=g+a(i,j)*(rhos**i)*Ts**(1-j)
        enddo
    enddo
    g=1.d0+g
	if (g.LE.0.d0) then
		g=0.001d0
	endif 
    F=exp(-gama*rhos*rhos)
    Gt(1)=(1.d0-F)/(2.d0*gama)
    Gt(2)=-(F*rhos*rhos-2.d0*Gt(1))/(2.d0*gama)
    Gt(3)=-(F*rhos*rhos**3-4.d0*Gt(2))/(2.d0*gama)
    Gt(4)=-(F*rhos*rhos**5-6.d0*Gt(3))/(2.d0*gama)
    Gt(5)=-(F*rhos*rhos**7-8.d0*Gt(4))/(2.d0*gama)
    Gt(6)=-(F*rhos*rhos**9-10.d0*Gt(5))/(2.d0*gama)
    a0=0.d0
    do i=1,8
        a0=a0+(at(i)*rhos**i)/i
    enddo
    do i=1,6
        a0=a0+bt(i)*Gt(i)
    enddo 
	eta=rhos*Pi/6.d0
	!aSS=(4.d0*eta-3.d0*eta*eta)/( (1.d0-eta)**2 )+log(eta)	  !this is the A-repulsive, SS stands for SoftSaft
	aRef=(a0*eok/T) !*m
	aChain=(1.d0-m)*log(g)	   !A/NKT where N is number of chain molecules
	aIdeal=log(rho)-1.d0
    !aSS=(a0*eok/T)*m 		 ! it is the residual part of A, A/NKT for chain
	aSS=aRef+aChain !+aIdeal

    return
    end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																											 C
!C Programmed by AFG 2011																					 C
!C This routine calculates the Helmholtz free energy using softSAFT EOS. It was used to validate the		 C
!C HelmRG subroutine and it is not used anywhere now in the project											 C
!C																											 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    subroutine AsoftSaft(tKelvin,rho,aSS,PSS,m,sig,eps)
 	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(etaStep=2.d0/1000.d0,maxRhoStep=0.8d0/etaStep) 
	PARAMETER(AvoNum=602.22d0,rGAs=8.314472d0,Pi=3.14159265359796d0,kB=0.0138065d0)
	DIMENSION a(5,5),x(32),at(8),bt(6),gt(6),gmol(55)
	double precision m
 
	eok=eps/kB
	a(1,1)=0.49304346593882d0
	a(2,1)=-0.47031983115362d0
	a(3,1)=5.0325486243620d0
	a(4,1)=-7.3633150434385d0
	a(5,1)=2.9043607296043d0
	a(1,2)=2.1528349894745d0
	a(1,3)=-15.955682329017d0
	a(1,4)=24.035999666294d0
	a(1,5)=-8.6437958513990d0
	a(2,2)=1.1471647487376d0
	a(2,3)=37.889828024211d0
	a(2,4)=-84.667121491179d0
	a(2,5)=39.643914108411d0
	a(3,2)=-25.915399226419d0
	a(3,3)=-18.862251310090d0
	a(3,4)=107.63707381726d0
	a(3,5)=-66.602649735720d0
	a(4,2)=51.553565337453d0
	a(4,3)=-40.519369256098d0
	a(4,4)=-38.796692647218d0
	a(4,5)=44.605139198318d0
	a(5,2)=-24.418812869291d0
	a(5,3)= 31.500186765040d0
	a(5,4)=-5.3368920371407d0
	a(5,5)=-9.5183440180133d0
	T=tKelvin
	Ts=T/eok
	Ts2=Ts*Ts
	x(1)=0.8623085097507421d0
	x(2)=2.976218765822098d0
	x(3)=-8.402230115796038d0
	x(4)=0.1054136629203555d0
	x(5)=-0.8564583828174598d0
	x(6)=1.582759470107601d0
	x(7)=0.7639421948305453d0
	x(8)=1.753173414312048d0
	x(9)=2798.291772190376d0
	x(10)=-0.048394220260857657d0
	x(11)=0.9963265197721935d0
	x(12)=-36.98000291272493d0
	x(13)=20.84012299434647d0
	x(14)=83.054021244717285d0
	x(15)=-957.4799715203068d0
	x(16)=-147.7746229234994d0
	x(17)=63.98607852471505d0
	x(18)=16.03993673294834d0
	x(19)=68.05916615864377d0
	x(20)=-2791.29358794945d0
	x(21)=-6.245138304568454d0
	x(22)=-8116.836104958410d0
	x(23)=14.88735559561229d0
	x(24)=-10593.46754655084d0
	x(25)=-113.1607632802822d0
	x(26)=-8867.771540418822d0
	x(27)=-39.86982844450543d0
	x(28)=-4689.27299917261d0
	x(29)=259.3535277438717d0
	x(30)=-2694.523589434903d0
	x(31)=-721.8487631550215d0
	x(32)=172.1802063863269d0
	at(1)=x(1)*Ts+x(2)*sqrt(Ts)+x(3)+x(4)/Ts+x(5)/Ts2
	at(2)=x(6)*Ts+x(7)+x(8)/Ts+x(9)/Ts2
	at(3)=x(10)*Ts+x(11)+x(12)/Ts
	at(4)=x(13)
	at(5)=x(14)/Ts+x(15)/Ts2
	at(6)=x(16)/Ts
	at(7)=x(17)/Ts+x(18)/Ts2
	at(8)=x(19)/Ts2
	bt(1)=x(20)/Ts2+x(21)/(Ts2*Ts)
	bt(2)=x(22)/Ts2+x(23)/(Ts2*Ts2)
	bt(3)=x(24)/Ts2+x(25)/(Ts2*Ts)
	bt(4)=x(26)/Ts2+x(27)/(Ts2*Ts2)
	bt(5)=x(28)/Ts2+x(29)/(Ts2*Ts)
	bt(6)=x(30)/Ts2+x(31)/(Ts2*Ts)+x(32)/(Ts2*Ts2)
	gama=3.d0

    rhos=rho*(sig*sig*sig) !*m		! this is reduced density for monomer
    g=0.d0
	dg=0.d0
    do i=1,5
        do j=1,5
            dg=dg+real(i)*a(i,j)*(rhos**(i+1))*Ts**(1-j)
			g=g+a(i,j)*(rhos**i)*Ts**(1-j)
        enddo
    enddo
    g=1.d0+g
	if (g.LE.0.d0) then
		g=0.001d0
	endif
    aChain=(1.d0-m)*log(g)	   !A/NKT where N is number of chain molecules
    
    F=exp(-gama*rhos*rhos)
    Gt(1)=(1.d0-F)/(2.d0*gama)
    Gt(2)=-(F*rhos*rhos-2.d0*Gt(1))/(2.d0*gama)
    Gt(3)=-(F*rhos*rhos**3-4.d0*Gt(2))/(2.d0*gama)
    Gt(4)=-(F*rhos*rhos**5-6.d0*Gt(3))/(2.d0*gama)
    Gt(5)=-(F*rhos*rhos**7-8.d0*Gt(4))/(2.d0*gama)
    Gt(6)=-(F*rhos*rhos**9-10.d0*Gt(5))/(2.d0*gama)
    a0=0.d0
	P0=0.d0
    do i=1,8
        a0=a0+(at(i)*rhos**i)/real(i)
		P0=P0+at(i)*rhos**(i+1)
    enddo
    do i=1,6
        a0=a0+bt(i)*Gt(i)
		P0=P0+(bt(i)*rhos**(2*i+1))*F
    enddo
	nComps=1
	gmol(1)=1
	eta=rhos*Pi/6.d0  ! based on monomer density
	call HelmRG(nComps,gmol,tKelvin,eta,bVolMix,aDepNC,zNC,PMpaNC,dP_drhoNC,d2P_drho2NC,d3P_drho3NC)
    aSS=(a0*eok/T)*m+aChain+aDepNC   !SS stands for SoftSaft
	Pchain=(1.d0-m)/m*rhos*Ts*(1.d0+dg/(g*rhos))
	PSS=( P0+rhos*Ts+Pchain )*eps/(sig*sig*sig)+PMpaNC
    return
    end

