	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine RegPureEsd2(NC)!note:tc,pc,acen,id & esd parms passed in common. 
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (4/2020) to use CritPure for PGL package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)	 ! max # of EOS parameters to estimate: c,eok,b,Kad,epsAD
	DoublePrecision	parm(3),error(3),stderr(3)
	integer iErr2(NC) ! iErrExactEsd list for each component.
	LOGICAL LOUDER
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDev2
	LOUDER=LOUD
	LOUDER=.TRUE.
	if(NC > 1)then
		pause 'RegpureEsd2 only works for Nc=1. Sorry.'
		return
	endif
	print*,'Initial guess from ParmsEsd? Enter 1 for yes or 0 to use MW guides as guess'
	read(*,*)iAns
	if(iAns==0)then
		c(1)=1+rmw(1)/150
		eokP(1)=359+.66/rmw(1)
		Vx(1)=17*rmw(1)/50
		parm(1)=c(1)
		if(ND(1)==0)call ExactEsd(nC,VX,c,q,eokP,iErr,iErr2)
		if(ND(1)==1)call RegPureDev3b(3,1,parm,error,iflag) !replaces eokP and Vx
	endif
	if(LOUDER)write(*,*)'RegPure1:  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!  general ***********************************
	q(1)=1+1.90476*(c(1)-1)
	nParms=3
	PARM(1)=c(1)
	PARM(2)=eokP(1)
	PARM(3)=Vx(1)
	call RegPureDev2(nData,nParms,parm,error,iflag)
	if(LOUDER)write(*,'(a,3f11.3,i3)')' C,Eok,vX',C(1),eokP(1),VX(1),Nd(1)
	if(LOUDER)write(*,*)'initial devA,devT,devP=',(error(i),i=1,3)
	aceCalc=acen(1)+error(1)
	TcCalc=Tc(1)+error(2)
	PcCalc=Pc(1)+error(3)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
	iflag=0
	iterOpt=0
	rguess=1234
	if(iteropt==0)then
		do while(rguess > 0)
			write(*,*)'enter C (-1 to skip to regression using last guess)'
			read(*,*)rguess,eokP(1),Vx(1)
			if(rguess.gt.0)then
				parm(1)=rguess
				if(ND(1)==1)call RegPureDev3b(3,1,parm,error,iflag) ! use old method as guess for eokP and Vx, given c
				parm(2)=eokP(1)
				parm(3)=vx(1)
				call RegPureDev2(nData,nParms,parm,error,iflag)
				if(LOUDER)write(*,'(a,1x,f9.4,f7.2,f9.4)')' initial devA,devT,devP=',(error(i),i=1,3)
				aceCalc=acen(1)+error(1)
				TcCalc=Tc(1)+error(2)
				PcCalc=Pc(1)+error(3)
				if(LOUDER)write(*,'(a,1x,f8.4,2f9.4,i3)')' c,eok,Vx=',(parm(i),i=1,3),Nd(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
			endif
		enddo
	endif
	q(1)=1+1.90476*(c(1)-1)
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!	pause
	factor=1.d-4
	tol=1.d-4
	nParms=3
	nData=3
	call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info > 4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	iflag=1
	call RegPureDev2(nData,nParms0,parm,error,iflag)
	qShape=1+1.90476*( c(1) -1 )
	rKadNm3=KcStar(1)
	dHkCal_mol=dH(1)*1.987*Tc(1)/1000
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kadNm3    dHkCal/mol nAs   nDs        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	if(LOUDER)write(*,*)'final devA,devT,devP=',(error(i),i=1,3)
	pause 'check errors'
	return
	end
	!******************************************************************************************************************************************
	SUBROUTINE RegPureDev2(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	dimension PARM(nParms),ERROR(nData),x(NMX) !,fugc(NMX) !,y(NMX),ierBp(12)
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	q(1)=1+1.90476*(c(1)-1)
	eokp(1)=Parm(2)
	VX(1)=parm(3)
	bVolCc_mol(1)=VX(1)
	isZiter=1 !use numerical derivatives
	toll=1.d-5
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUD)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	call CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,acen_pure,iErrCode)
	error(1)=acen_pure-acen(1)
	error(2)=TC_Pure-Tc(1)
	error(3)=PC_Pure-Pc(1)
	return
	end
!******************************************************************************************************************************************

	subroutine RegPureEsd2c(dHkJ_molP,rKadStarP,nParms0,iterOpt,iDenOptP,idToReg,nComps)!note:tc,pc,acen,id & esd parms passed in common
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (3/01) to put option in PrEsd package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDev2
	dimension PARM(NPMAX),ERROR(55),stdErr(NMX)
	common/ PUREVP /VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ CONSTK /KIJ(NMX,NMX),INITIAL
	common/ AVERR  /AAPERR,RHOERR,objfun
	common/ CONST  /rKadStar, dHkJ_mol, iDenOpt
	rKadStar=rKadStarP
	dHkJ_mol=dHkJ_molP
	iDenOpt = iDenOptP
	!if(rGas.lt.1.987)rGas=8.314
	do iComp=1,nComps
		if(id(iComp).eq.idToReg)iCompReg=iComp
	enddo

	nC=1 ! within RegPure, we must limit to only 1 compo.
	if(iCompReg.ne.1)then !the comp for regression must be in 1 position so calls to fugi etc sum only over that component
		tcStore=tc(1)
		pcStore=pc(1)
		acenStore=acen(1)
		eokPStore=eokP(1)
		kcStarStore=kcStar(1)
		dHStore=dH(1)
		cShapeStore=c(1)
		QStore=Q(1)
		vxStore=VX(1)
		ndStore=nd(1)
		tc(1)=tc(iCompReg)
		pc(1)=pc(iCompReg)
		acen(1)=acen(iCompReg)
	endif

	c(1)=1+acen(1)*2
	c(1)=1+rmw(1)/150
	u0ok=359+.66/rmw(1)
	v00=17*rmw(1)/50
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/8.314/Tc(1)
	!  general ***********************************
	PARM(1)=c(1)
	PARM(2)=v00
	PARM(3)=359+.66/rmw(1)
	PARM(4)=kcStar(1)
	PARM(5)=dH(1)
	nData=16 !this gives lots of points between tMax and tMin
	nParms=1 !this will force crit match estimate of bVol,eokP
	call RegPureDev2(nData,nParms,parm,error,iflag)
	if(LOUD)write(*,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	if(LOUD)write(*,*)'initial aaperr,sgerr',aaperr,rhoerr
	!outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	!open(66,file=outFile)
	!66 should be open from RegPureIo.
	write(66,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	write(66,*)'initial aaperr,sgerr',aaperr,rhoerr
	objo=objfun
	aaperro=aaperr
	if(iteropt.eq.0)then
		do while(rguess.gt.0)
			write(*,*)'enter C,vx,eok (-1 to skip to regression using last guess)'
			read(*,*)rguess,vguess,u0ok
			if(rguess.gt.0)then
				parm(1)=rguess
				parm(2)=vguess
				parm(3)=u0ok
				r=rguess
				v00=vguess
				call RegPureDev2(nData,nParms0,parm,error,iflag)
				write(*,*)'aaperr,sgerr',aaperr,rhoerr
			endif
		enddo
	else
		!	if(nParms.ge.2)then	!search for starting v00
		!		v00=v00/1.1
		!		parm(2)=v00
		!		call DevPure(nData,nParms,parm,error,iflag)
		!	    write(*,*)'aaperr,sgerr',aaperr,rhoerr
		!	    if(aaperr.lt.aaperro)then
		!	      aaperro=aaperr
		!		  goto 3000
		!	    endif
		!		parm(2)=parm(2)*1.1
		!	end if
		nParms=1 !get starting value by optimizing cShape only.
		NPO=ND(1)
		FACTOR=1.d-4
		TOL=1.D-4
		call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
		parm(2)=VX(1)
		if(nParms0.ge.3)then
			nParms=2
			FACTOR=1.d-4
			TOL=1.D-4
			aaperrB4=aaperr
			if(LOUD)write(*,*)'trying nParms=2'
			call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
			if(aaperr.gt.aaperrB4)then
				if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
				if(LOUD)pause 'Warning at nParms=2: aaperr.gt.aaperrB4'
			endif
		endif
	endif
	NPO=ND(1)
	nParms=nParms0
	qShape=1+1.90476*(C(1)-1)
	dHStar=dHkJ_mol*1000/8.314/tc(1)
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	aaperrB4=aaperr												   
	if(LOUD)write(*,*)'one last call to LmDif, aaperr ',aaperr
	!	pause
	factor=1.d-4
	tol=1.d-4
	parm(3)=eokP(1)
	!	nParms=3
	call LmDifEz(RegPureDev2,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info.gt.4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	if(aaperr.gt.aaperrB4)then
		if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
		if(LOUD)pause 'Warning: aaperr.gt.aaperrB4'
	endif

	!evaluate actual vp error
	iFlag=86
	call RegPureDev2(nData,nParms,parm,error,iflag)
	qShape=1+1.90476*(C(1)-1)
	dHout=dh(1)*Tc(1)*8.314
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	WRITE(66,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	WRITE(66,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
601	FORMAT(I5,1X,2(F9.3,1X),2(F9.3,1X),F9.5,2I4,2(F7.2,1X),f7.4,1x,f7.2,1x,f7.4)
	if(LOUD)PAUSE 'check final results'
86	continue
	if(iCompReg.ne.1)then
		tc(iCompReg)=tc(1)
		pc(iCompReg)=pc(1)
		acen(iCompReg)=acen(1)
		eokP(iCompReg)=eokP(1)
		kcStar(iCompReg)=kcStar(1)
		dH(iCompReg)=dH(1)
		C(iCompReg)=C(1)
		Q(iCompReg)=Q(1)
		VX(iCompReg)=VX(1)
		nd(iCompReg)=nd(1)
		tc(1)=tcStore
		pc(1)=pcStore
		acen(1)=acenStore
		eokP(1)=eokPStore
		kcStar(1)=kcStarStore
		dH(1)=dHStore
		C(1)=cShapeStore
		Q(1)=QStore
		VX(1)=vxStore
		nd(1)=ndStore
	endif
	return
	end

	SUBROUTINE RegPureDev2c(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	INTEGER LIQ
	dimension ier(12)
	dimension PARM(*),ERROR(nData),fugc(NMX),x(NMX),y(NMX),ierBp(12)
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP !we need accurate VP for regression
	common/ CONST / rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	IF(TMAXVP/TC(1).LT.0.8)WRITE(66,*)'Warning TMAX<0.8Tc is TOO LOW'
	IF(TMINVP/TC(1).GT.0.6)WRITE(66,*)'Warning TMIN>0.6Tc is TOO HIGH'
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/Tc(1)/8.314

	!TR7=0.7*TC(1)
	!PR7=EXP( VPA+VPC*LOG(TR7)+VPB/TR7+VPD*TR7**VPE)/PC(1)
	!ACEN(1)=-LOG10(PR7)-1
	!???We must know ACEN from input
	IF(nParms.GE.2)VX(1)    =PARM(2)
	IF(nParms.GE.3)eokP(1)  =PARM(3)
	IF(nParms.GE.4)kcStar(1)=PARM(4)
	IF(nParms.GE.5)dH(1)    =PARM(5)
	Q(1)=1+1.90476*(C(1)-1)
	if(dh(1).gt.111)then
		if(LOUD)write(*,*)'dH is too large.  Spec kJ not J.'
		if(LOUD)pause
		return
	endif
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	if(nParms.ge.3)goto 2000
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev2 during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms.eq.1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms.le.2)eokP(1)=LOG(1.0617+Yc)*Tc(1)
2000 continue		
    !if(LOUD)write(*,'(a,5f10.3)')' c,Vx,eps/k',C(1),VX(1),eokP(1)
    !pause 'check values'
    if(nData==1)then !use the acentric factor to estimate the vp error
        pSat7 = Pc(1)*10**( 7*(1+acen(1))/3*(1-1/0.7) )
        call Fugi(TK,P,X,NC,1,FUGC,ZL,IER)
        fugcL=fugc(1)
        call Fugi(TK,pSat7,X,NC,0,FUGC,ZV,IER)
        ERROR(1)=( EXP( FUGCL-fugc(1) )-1 )*100
        return
    endif 
	tmin=TminVP*1.1
	if(tmin.lt. 0.45*Tc(1))tmin=0.45*Tc(1)
	tMax=tMaxVp*0.9
	if(tMax.gt. 700)tMax=700
	DELT=(TMAX-TMIN)/(nData-1)
	PERR=0
	IPTS=0
	iTempMax=nData-1
	DO IT=1,iTempMax,1
!		TK=TMIN+DELT*(IT-1)
		TK=TMIN+DELT*(IT)
		TR=TK/TC(1)
		P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE )
		if(iFlag.eq.86)then
			pCalc=P
			init=1
			itMax=111
			calL BUBPL(tK,x,nC,init,pCalc,itMax,y,ierBp)
			perri=(pCalc-P)/P*100
			if(LOUD)write(*,'(f7.2,4f10.5)')tK,pCalc,P,perri
			IPTS=IPTS+1
		else
			LIQ=1
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
			fugcL=fugc(1)
			LIQ=0
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZV,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZV,IER)
			fugcV=fugc(1)
			IPTS=IPTS+1
			perri=( EXP( FUGCL-FUGCV )-1 )*100
		endif
		ERROR(IPTS)=perri
		if(kcStar(1).lt.0)error(ipts)=error(ipts)*10000
		if(C(1).lt.1)error(ipts)=error(ipts)*10000
		if(abs(zl-zv).lt.1.d-2)error(ipts)=error(ipts-1)*5
		if(ier(1).ne.0)error(ipts)=error(ipts-1)*10
!		if(abs(u0ok-u0okest)/u0ok.gt.0.03)error(ipts)=error(ipts)*1000
		PERR=PERR+ABS(perri)
		objfun=objfun+error(ipts)*error(ipts)
	enddo

	AAPERR=PERR/(nData-1)

	TK=298
	!P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE ) !
	P=PC(1)*10**( 7/3*(1+ACEN(1))*(1-TC(1)/TK) ) !We don't need super accurate P to estimate rhoLiq.
	if(P.lt.0.1)P=0.1
	LIQ=1
	call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
	!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
	RHOL=RMW(1)*P/8.314/TK/ZL
	IPTS=IPTS+1
	RHOERR=(SG-RHOL)/SG*100
	if(iDenOpt.eq.1)ERROR(nData)=iDenOpt*RHOERR
	if(LOUD)write(*,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	WRITE(66,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	!pause 'check deviations'
	!ERROR(nData)=0
	RETURN
	END

!***********************************************************************
	SUBROUTINE RegPureDev3b(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!INTEGER LIQ
	!dimension ier(12)
	dimension PARM(*),ERROR(nData) ,fugc(NMX),x(NMX) !,y(NMX),ierBp(12)
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	Q(1)=1+1.90476*(C(1)-1)
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev3 during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms==1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms < 3)eokP(1)=LOG(1.0617+Yc)*Tc(1)
	parm(2)=eokP(1)
	parm(3)=VX(1)
	pSat7=Pc(1)*10**( -1-acen(1) )
	T=Tc(1)*0.7D0
	call FugiESD(T,pSat7,X,nC,1,FUGC,rho,zFactor,Ares,Ures,iErrF)
	chemPoLiq=FUGC(1)
	call FugiESD(T,pSat7,X,nC,0,FUGC,rho,zFactor,Ares,Ures,iErrF)
	chemPoVap=FUGC(1)
	error(1)=chemPoLiq-chemPoVap
	if(iFlag==1.and.LOUD)print*,'RegPureDevEsd3: c,chemPoDev',C(1),error(1) !mostly to avoid iFlag warning.
	RETURN
	END

!***********************************************************************
      SUBROUTINE RegPureIo2(NC)
	!C
	!C  PURPOSE:  Regress Pure component parameters of ESD model.
	!C  PROGRAMMED BY:  JRE 3/01
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER fileNameVp*55
	character outFile*251
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ CONST /  rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	print*,'Minimize vp error based on estd rKadStar and dHkJ_mol '  
	print*,'constrained to match Tc. Pc is ~matched implicitly by vp' 
	print*,'or matched explicitly when nParms=1.'
	print*,'best setting is nParms=2 and iteropt=0 => flex to vx.'   
	print*,'kcStar comes from kcStar~kadStar/cShape so kadStar=.025'
	print*,'Data from c:\thermo\esd\VpCcEsdCorr.txt if available.'
	print*,'Enter dHkcal/mol,rKadStar(eg 0.025),iDenOpt(0=ignore SG)'
	read(*,*)dHkJ_mol,rKadStar,iDenOpt
	dHkJ_mol=dHkJ_mol*4.184
	!rKadStar=0.025
	!iDenOpt=0		!0=ignore density in fit, 1=weight sg as one datum
	!initial=0

	WRITE(*,*)'ENTER nParms,iteropt(0=AutoReg, 1=guess then reg)'
	READ(*,*)nParms0,iteropt
	idToReg=id(1)
	ndTemp=1
	iCompToReg=1
	if(nc.ne.1)then
		write(*,*)'Enter dippr id # and Nd of component for regression '
		read(*,*)idToReg,ndTemp
		do iComp=1,nc
			if(id(iComp).eq.idToReg)iCompToReg=iComp
		enddo
		nd(iCompToReg)=ndTemp
	endif
	fileNameVp='c:\thermo\esd\VpCcEsdCorr.txt'
	!write(*,*)'Enter filename for vpdata'
	!read(*,'(a)')fileNameVp
	iFound=0
	open(55,file=fileNameVp,ERR=86)
	read(55,*,END=86,ERR=86)nDbaseVp !error here means no file so skip cc search.
	write(*,*)'Enter chemcad id # of component for regression '
	read(*,*)idToRegCC
	do iVp=1,nDbaseVp
		READ(55,*,END=86)idDb,rmwDb,sgDb,vpADb,vpBDb,vpCDb,vpDDb,vpEDb,tMinVpDb,tMaxVpDb
		if(idDb.eq.idToRegCC)then
			iFound=1
			rmwDum=rmwDb
			sg=sgDb
			vpA=vpADb
			vpB=vpBDb
			vpC=vpCDb
			vpD=vpDDb
			vpE=vpEDb
			tMinVp=tMinVpDb
			tMaxVp=tMaxVpDb
		endif
	enddo
	CLOSE(55)
86	continue !if the dbase doesn't exist, then we must type the data in.
	if(iFound.eq.0)then
		write(*,*)'Compd not found in std dbase.'
		write(*,*)'Enter Mw,SpeGrav,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp'
		read(*,*)rMwDum,sg,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp
	endif
	
	vpA=vpA-LOG(1.E6)!convert dippr vp from Pa to MPa.

	TRLO=tMinVp/TC(iCompToReg)
	TRHI=tMaxVp/TC(iCompToReg)
	outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	open(66,file=outFile)
	IF(TRLO.GT.0.6.OR.TRHI.LT.0.8)THEN
	  if(LOUD)write(* ,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	  WRITE(66,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	END IF

	!  we need a big function here because of switching comp pos., and 
	!initial guesses to converge.
	call RegPureEsd2c(dHkJ_mol,rKadStar,nParms0,iterOpt,iDenOpt,id(iCompToReg),nC)!note:tc,pc,acen,id & esd parms passed in common
	close(66)
	if(LOUD)pause 'Results in RegPureOut.txt'
	return
	end
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE KsVAL2(nc,nParms,parm,tKelvin,deviate,kAADk,kBias,kErrMax,LAST)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DOUBLE PRECISION parm(nParms),deviate(maxPts)
	parameter (nEta=5)
	common/SimData/ xFracc(23),vEffNm3(23),a0xss(20,20),A0Mix(20,20),A1Mix(20,20),A2Mix(20,20),Nps
	!common/XsProps/a0XsSs
	!dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	!common/ksvall/ks0ij,ks1ij
	DIMENSION XXFRAC(NMX)
            
	607	FORMAT(1X,F8.4,1X,F10.4,1X,F10.8,1X,F10.8,1X,F10.4,1X,i5)      

	!KS0IJ(1,2)=0
	!KS1IJ(1,2)=0
	!KS0IJ(2,1)=KS0IJ(1,2)
	!KS1IJ(2,1)=KS1IJ(1,2)
	!KS0IJ(1,1)=0
	!KS1IJ(1,1)=0

		if(last.lt.0)then !only kijOpt sets last to -1, so that is the only one that effects this.
			write(*,*)'Enter optimization choice: 1=ks0, 2=ks0,ks1'
			read(*,*)iOpt
		endif
		if (nParms==1)then
			KS0IJ(1,2)=parm(1)
		else
			KS0IJ(1,2)=parm(1)
			KS1IJ(1,2)=parm(2)
		endif
	MAXIT=1111
	KAADK=0
	ssqErr=0
	kBias=0
	eta=0.d0
	!tKelvin=273
	ikk=1
	nTot=nEta*(NPS-2)
	!Do ikk=1,nTot
	nFracs=NPS-2
		Do iEta=1,nEta
			eta=eta+0.1d0

			DO iData=1,nFracs
			!if(iData==(NPS-1))then
			!	xxFrac(1)=1
			!	xxFrac(2)=0
			!elseif(iData==NPS)then
			!	xxFrac(1)=0
			!	xxFrac(2)=1
			!else
				xxFrac(1)=xFracc(iData)
				xxFrac(2)=1-xxFrac(1)
			!endif
				IFLAG=0
				ITMAX=MAXIT
				INIT=1
				nComps=nc  
				CALL AxsCalc(xxFrac,eta,nComps,tKelvin,a0MixCalc,iErr)
				!a0XsCalc=a0MixCalc-(xxFrac(1)*A0Mix((NPS-1),iEta)+xxFrac(2)*A0Mix(NPS,iEta)) 
			!a0Xs=a0Mix-(xxFrac(1)*

				deviate(ikk) =(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0XsCalc-a0xss(iData,iEta))/a0xss(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))**2! !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100	 !
				if(ABS(deviate(ikk)).gt.ABS(kErrMax))kErrMax=deviate(ikk)
				
				!penalize the deviate passed to lmdif, but only compile paadp for converged pts
				KAADK = KAADK+ABS(deviate(ikk))
				ssqErr= ssqErr+deviate(ikk)*deviate(ikk)
				kBias=kBias+deviate(ikk)
		!	if(ABS(deviate(iData)).gt.ABS(pErrMax))pErrMax=deviate(iData)
			
				!WRITE(67,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)

				ikk=ikk+1
			enddo
		enddo
	!enddo

	rmsErr=SQRT(ssqErr/(nEta*(NPS-2)))
	KAADK=KAADK/(nEta*(NPS-2))
	kBias=kBias/(nEta*(NPS-2))
	if(LOUD)write(*,'(7f10.3)')(parm(i),i=1,nparms),KAADK,rmsErr

	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine RegPureEsdb(NC)!note:tc,pc,acen,id & esd parms passed in common. 
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (4/2020) to use CritPure for PGL package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)	 ! max # of EOS parameters to estimate: c,eok,b,Kad,epsAD
	DoublePrecision	parm(3),error(3),stderr(3)
	integer iErr2(NC) ! iErrExactEsd list for each component.
	LOGICAL LOUDER
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDevb
	LOUDER=LOUD
	LOUDER=.TRUE.
	if(NC > 1)then
		pause 'RegpureEsdb only works for Nc=1. Sorry.'
		return
	endif
	print*,'Initial guess from ParmsEsd? Enter 1 for yes or 0 to use MW guides as guess'
	read(*,*)iAns
	if(iAns==0)then
		c(1)=1+rMw(1)/150
		eokP(1)=359+.66/rMw(1)
		Vx(1)=17*rMw(1)/50
		parm(1)=c(1)
		if(ND(1)==0)call ExactEsd(nC,VX,c,q,eokP,iErr,iErr2)
		if(ND(1)==1)call RegPureDev3(3,1,parm,error,iflag) !replaces eokP and Vx
	endif
	if(LOUDER)write(*,*)'RegPure1:  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!  general ***********************************
	q(1)=1+1.90476*(c(1)-1)
	nParms=3
	PARM(1)=c(1)
	PARM(2)=eokP(1)
	PARM(3)=Vx(1)
	call RegPureDevb(nData,nParms,parm,error,iflag)
	if(LOUDER)write(*,'(a,3f11.3,i3)')' C,Eok,vX',C(1),eokP(1),VX(1),Nd(1)
	if(LOUDER)write(*,*)'initial devA,devT,devP=',(error(i),i=1,3)
	aceCalc=acen(1)+error(1)
	TcCalc=Tc(1)+error(2)
	PcCalc=Pc(1)+error(3)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
	if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
	iflag=0
	iterOpt=0
	rguess=1234
	if(iteropt==0)then
		do while(rguess > 0)
			write(*,*)'enter C (-1 to skip to regression using last guess)'
			read(*,*)rguess,eokP(1),Vx(1)
			if(rguess.gt.0)then
				parm(1)=rguess
				if(ND(1)==1)call RegPureDev3(3,1,parm,error,iflag) ! use old method as guess for eokP and Vx, given c
				parm(2)=eokP(1)
				parm(3)=vx(1)
				call RegPureDevb(nData,nParms,parm,error,iflag)
				if(LOUDER)write(*,'(a,1x,f9.4,f7.2,f9.4)')' initial devA,devT,devP=',(error(i),i=1,3)
				aceCalc=acen(1)+error(1)
				TcCalc=Tc(1)+error(2)
				PcCalc=Pc(1)+error(3)
				if(LOUDER)write(*,'(a,1x,f8.4,2f9.4,i3)')' c,eok,Vx=',(parm(i),i=1,3),Nd(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Exptl Tc,Pc,w=',Tc(1),Pc(1),acen(1)
				if(LOUDER)write(*,'(a,1x,f7.2,2f9.4)')' Calcd Tc,Pc,w=',TcCalc,PcCalc,aceCalc
			endif
		enddo
	endif
	q(1)=1+1.90476*(c(1)-1)
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	!	pause
	factor=1.d-4
	tol=1.d-4
	nParms=3
	nData=3
	call LmDifEz(RegPureDevb,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info > 4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	iflag=1
	call RegPureDevb(nData,nParms0,parm,error,iflag)
	qShape=1+1.90476*( c(1) -1 )
	rKadNm3=KcStar(1)
	dHkCal_mol=dH(1)*1.987*Tc(1)/1000
	if(LOUDER)write(*,*)'  ID       c     q     E/k       Vx    Nd  kadNm3    dHkCal/mol nAs   nDs        '
	if(LOUDER)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	if(LOUDER)write(*,*)'final devA,devT,devP=',(error(i),i=1,3)
	pause 'check errors'
	return
	end
	!******************************************************************************************************************************************
	SUBROUTINE RegPureDevb(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	dimension PARM(nParms),ERROR(nData),x(NMX) !,fugc(NMX) !,y(NMX),ierBp(12)
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	q(1)=1+1.90476*(c(1)-1)
	eokp(1)=Parm(2)
	VX(1)=parm(3)
	bVolCc_mol(1)=VX(1)
	isZiter=1 !use numerical derivatives
	toll=1.d-5
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH        '
	if(LOUD)write(*,'(i9,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,2i3,2f8.2)')ID(1),C(1),q(1),eokP(1),VX(1),ND(1),kcStar(1),dH(1)*1.987*Tc(1)/1000,nAs(1),nDs(1)
	call CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,acen_pure,iErrCode)
	error(1)=acen_pure-acen(1)
	error(2)=TC_Pure-Tc(1)
	error(3)=PC_Pure-Pc(1)
	return
	end
!******************************************************************************************************************************************

	subroutine RegPureEsd(dHkJ_molP,rKadStarP,nParms0,iterOpt,iDenOptP,idToReg,nComps)!note:tc,pc,acen,id & esd parms passed in common
	! THIS ROUTINE minimizes vp error based on est'd rKadStar and dHkJ_mol (spec'd in main),  
	! and constrained to match Tc.  Pc is ~matched implicitly by fitting vp near Tc.
	! or matched explicitly when nParms=1.
	! best setting is nParms=2 and iteropt=0, gives some flex to vx.   
	! kcStar comes from kcStar~.025/cShape because Vseg~bmol/cShape.
	! the calc of Tc is a bit crude: just
	! do 100 weighted successive subst iterations whether you need it or not
	! Routines Req'd: LmDifEz(MinPack),Fugi 
	! Programmed by: JRE (9/96)
	! Revised:  JRE (3/01) to put option in PrEsd package

	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DoublePrecision( A-H, K,O-Z )
	PARAMETER(NPMAX=5)
	!CHARACTER*44 FNVP,FNLD
	EXTERNAL RegPureDev
	dimension PARM(NPMAX),ERROR(55),stdErr(NMX)
	common/ PUREVP /VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ CONSTK /KIJ(NMX,NMX),INITIAL
	common/ AVERR  /AAPERR,RHOERR,objfun
	common/ CONST  /rKadStar, dHkJ_mol, iDenOpt
	rKadStar=rKadStarP
	dHkJ_mol=dHkJ_molP
	iDenOpt = iDenOptP
	!if(rGas.lt.1.987)rGas=8.314
	do iComp=1,nComps
		if(id(iComp).eq.idToReg)iCompReg=iComp
	enddo

	nC=1 ! within RegPure, we must limit to only 1 compo.
	if(iCompReg.ne.1)then !the comp for regression must be in 1 position so calls to fugi etc sum only over that component
		tcStore=tc(1)
		pcStore=pc(1)
		acenStore=acen(1)
		eokPStore=eokP(1)
		kcStarStore=kcStar(1)
		dHStore=dH(1)
		cShapeStore=c(1)
		QStore=Q(1)
		vxStore=VX(1)
		ndStore=nd(1)
		tc(1)=tc(iCompReg)
		pc(1)=pc(iCompReg)
		acen(1)=acen(iCompReg)
	endif

	c(1)=1+acen(1)*2
	c(1)=1+rmw(1)/150
	u0ok=359+.66/rmw(1)
	v00=17*rmw(1)/50
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/8.314/Tc(1)
	!  general ***********************************
	PARM(1)=c(1)
	PARM(2)=v00
	PARM(3)=359+.66/rmw(1)
	PARM(4)=kcStar(1)
	PARM(5)=dH(1)
	nData=16 !this gives lots of points between tMax and tMin
	nParms=1 !this will force crit match estimate of bVol,eokP
	call RegPureDev(nData,nParms,parm,error,iflag)
	if(LOUD)write(*,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	if(LOUD)write(*,*)'initial aaperr,sgerr',aaperr,rhoerr
	!outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	!open(66,file=outFile)
	!66 should be open from RegPureIo.
	write(66,*)'C,vX,Eok',C(1),VX(1),eokP(1)
	write(66,*)'initial aaperr,sgerr',aaperr,rhoerr
	objo=objfun
	aaperro=aaperr
	if(iteropt.eq.0)then
		do while(rguess.gt.0)
			write(*,*)'enter C,vx,eok (-1 to skip to regression using last guess)'
			read(*,*)rguess,vguess,u0ok
			if(rguess.gt.0)then
				parm(1)=rguess
				parm(2)=vguess
				parm(3)=u0ok
				r=rguess
				v00=vguess
				call RegPureDev(nData,nParms0,parm,error,iflag)
				write(*,*)'aaperr,sgerr',aaperr,rhoerr
			endif
		enddo
	else
		!	if(nParms.ge.2)then	!search for starting v00
		!		v00=v00/1.1
		!		parm(2)=v00
		!		call DevPure(nData,nParms,parm,error,iflag)
		!	    write(*,*)'aaperr,sgerr',aaperr,rhoerr
		!	    if(aaperr.lt.aaperro)then
		!	      aaperro=aaperr
		!		  goto 3000
		!	    endif
		!		parm(2)=parm(2)*1.1
		!	end if
		nParms=1 !get starting value by optimizing cShape only.
		NPO=ND(1)
		FACTOR=1.d-4
		TOL=1.D-4
		call LmDifEz(RegPureDev,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
		parm(2)=VX(1)
		if(nParms0.ge.3)then
			nParms=2
			FACTOR=1.d-4
			TOL=1.D-4
			aaperrB4=aaperr
			if(LOUD)write(*,*)'trying nParms=2'
			call LmDifEz(RegPureDev,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
			if(aaperr.gt.aaperrB4)then
				if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
				if(LOUD)pause 'Warning at nParms=2: aaperr.gt.aaperrB4'
			endif
		endif
	endif
	NPO=ND(1)
	nParms=nParms0
	qShape=1+1.90476*(C(1)-1)
	dHStar=dHkJ_mol*1000/8.314/tc(1)
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dH     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.3,f7.3,i3,e11.4,f9.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	aaperrB4=aaperr												   
	if(LOUD)write(*,*)'one last call to LmDif, aaperr ',aaperr
	!	pause
	factor=1.d-4
	tol=1.d-4
	parm(3)=eokP(1)
	!	nParms=3
	call LmDifEz(RegPureDev,nData,nParms,PARM,factor,ERROR,TOL,INFO,stdErr)
	if(info.gt.4 .and. LOUD)write(*,*)'error in lmdif 2nd call'
	if(aaperr.gt.aaperrB4)then
		if(LOUD)write(*,*)'id,aaperr,aaperrB4',id(1),aaperr,aaperrB4
		if(LOUD)pause 'Warning: aaperr.gt.aaperrB4'
	endif

	!evaluate actual vp error
	iFlag=86
	call RegPureDev(nData,nParms,parm,error,iflag)
	qShape=1+1.90476*(C(1)-1)
	dHout=dh(1)*Tc(1)*8.314
	if(LOUD)write(*,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	if(LOUD)write(*,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
	WRITE(66,*)'  ID       c     q     E/k       Vx    Nd  kcStar    dHr     err   '
	WRITE(66,'(i5,2f8.4,f9.2,f7.2,i3,e11.4,f7.2,i4,2f8.2)')ID(1),C(1),qShape,eokP(1),VX(1),NPO,kcStar(1),dHStar,NPO,AAPERR,RHOERR
601	FORMAT(I5,1X,2(F9.3,1X),2(F9.3,1X),F9.5,2I4,2(F7.2,1X),f7.4,1x,f7.2,1x,f7.4)
	if(LOUD)PAUSE 'check final results'
86	continue
	if(iCompReg.ne.1)then
		tc(iCompReg)=tc(1)
		pc(iCompReg)=pc(1)
		acen(iCompReg)=acen(1)
		eokP(iCompReg)=eokP(1)
		kcStar(iCompReg)=kcStar(1)
		dH(iCompReg)=dH(1)
		C(iCompReg)=C(1)
		Q(iCompReg)=Q(1)
		VX(iCompReg)=VX(1)
		nd(iCompReg)=nd(1)
		tc(1)=tcStore
		pc(1)=pcStore
		acen(1)=acenStore
		eokP(1)=eokPStore
		kcStar(1)=kcStarStore
		dH(1)=dHStore
		C(1)=cShapeStore
		Q(1)=QStore
		VX(1)=vxStore
		nd(1)=ndStore
	endif
	return
	end

	SUBROUTINE RegPureDev(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	INTEGER LIQ
	dimension ier(12)
	dimension PARM(*),ERROR(nData),fugc(NMX),x(NMX),y(NMX),ierBp(12)
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP !we need accurate VP for regression
	common/ CONST / rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	IF(TMAXVP/TC(1) < 0.8)WRITE(66,*)'Warning TMAX<0.8Tc is TOO LOW'
	IF(TMINVP/TC(1) > 0.6)WRITE(66,*)'Warning TMIN>0.6Tc is TOO HIGH'
	IF(IFLAG.NE.0)ABC=123	 !this is here so no warning for not using iflag
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	kcStar(1)=rKadStar/C(1)
	dH(1)=dHkJ_mol*1000/Tc(1)/8.314

	!TR7=0.7*TC(1)
	!PR7=EXP( VPA+VPC*LOG(TR7)+VPB/TR7+VPD*TR7**VPE)/PC(1)
	!ACEN(1)=-LOG10(PR7)-1
	!???We must know ACEN from input
	IF(nParms.GE.2)VX(1)    =PARM(2)
	IF(nParms.GE.3)eokP(1)  =PARM(3)
	IF(nParms.GE.4)kcStar(1)=PARM(4)
	IF(nParms.GE.5)dH(1)    =PARM(5)
	Q(1)=1+1.90476*(C(1)-1)
	if(dh(1).gt.111)then
		if(LOUD)write(*,*)'dH is too large.  Spec kJ not J.'
		if(LOUD)pause
		return
	endif
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	if(nParms.ge.3)goto 2000
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms.eq.1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms.le.2)eokP(1)=LOG(1.0617+Yc)*Tc(1)
2000 continue		
    !if(LOUD)write(*,'(a,5f10.3)')' c,Vx,eps/k',C(1),VX(1),eokP(1)
    !pause 'check values'
    if(nData==1)then !use the acentric factor to estimate the vp error
        pSat7 = Pc(1)*10**( 7*(1+acen(1))/3*(1-1/0.7) )
        call Fugi(TK,P,X,NC,1,FUGC,ZL,IER)
        fugcL=fugc(1)
        call Fugi(TK,pSat7,X,NC,0,FUGC,ZV,IER)
        ERROR(1)=( EXP( FUGCL-fugc(1) )-1 )*100
        return
    endif 
	tmin=TminVP*1.1
	if(tmin.lt. 0.45*Tc(1))tmin=0.45*Tc(1)
	tMax=tMaxVp*0.9
	if(tMax.gt. 700)tMax=700
	DELT=(TMAX-TMIN)/(nData-1)
	PERR=0
	IPTS=0
	iTempMax=nData-1
	DO IT=1,iTempMax,1
!		TK=TMIN+DELT*(IT-1)
		TK=TMIN+DELT*(IT)
		TR=TK/TC(1)
		P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE )
		if(iFlag.eq.86)then
			pCalc=P
			init=1
			itMax=111
			calL BUBPL(tK,x,nC,init,pCalc,itMax,y,ierBp)
			perri=(pCalc-P)/P*100
			if(LOUD)write(*,'(f7.2,4f10.5)')tK,pCalc,P,perri
			IPTS=IPTS+1
		else
			LIQ=1
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
			fugcL=fugc(1)
			LIQ=0
			call Fugi(TK,P,X,NC,LIQ,FUGC,ZV,IER)
			!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZV,IER)
			fugcV=fugc(1)
			IPTS=IPTS+1
			perri=( EXP( FUGCL-FUGCV )-1 )*100
		endif
		ERROR(IPTS)=perri
		if(kcStar(1).lt.0)error(ipts)=error(ipts)*10000
		if(C(1).lt.1)error(ipts)=error(ipts)*10000
		if(abs(zl-zv).lt.1.d-2)error(ipts)=error(ipts-1)*5
		if(ier(1).ne.0)error(ipts)=error(ipts-1)*10
!		if(abs(u0ok-u0okest)/u0ok.gt.0.03)error(ipts)=error(ipts)*1000
		PERR=PERR+ABS(perri)
		objfun=objfun+error(ipts)*error(ipts)
	enddo

	AAPERR=PERR/(nData-1)

	TK=298
	!P=EXP( VPA+VPC*LOG(TK)+VPB/TK+VPD*TK**VPE ) !
	P=PC(1)*10**( 7/3*(1+ACEN(1))*(1-TC(1)/TK) ) !We don't need super accurate P to estimate rhoLiq.
	if(P.lt.0.1)P=0.1
	LIQ=1
	call Fugi(TK,P,X,NC,LIQ,FUGC,ZL,IER)
	!call FUGI(TCp,PCp,ACEN,ID,rGas,TK,P,X,NC,LIQ,FUGC,ZL,IER)
	RHOL=RMW(1)*P/8.314/TK/ZL
	IPTS=IPTS+1
	RHOERR=(SG-RHOL)/SG*100
	if(iDenOpt.eq.1)ERROR(nData)=iDenOpt*RHOERR
	if(LOUD)write(*,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	WRITE(66,'(a,9f10.3)')' AAPERR,%RHOERR,parm',AAPERR,RHOERR,(parm(i),i=1,nParms)
	!pause 'check deviations'
	!ERROR(nData)=0
	RETURN
	END

!***********************************************************************
	SUBROUTINE RegPureDev3(nData,nParms,PARM,ERROR,IFLAG)
	!nParms=1 => use critpt for eok and bvol
	!nParms=2 => use critpt for eok
	!nParms>2 => ignore critpt
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!INTEGER LIQ
	!dimension ier(12)
	dimension PARM(*),ERROR(nData) ,fugc(NMX),x(NMX) !,y(NMX),ierBp(12)
	NC=1
	x(1)=1
	C(1)=  PARM(1)
	Q(1)=1+1.90476*(C(1)-1)
	Fhb=EXP(dH(1))-1
	rootc=SQRT(C(1))
 	Zci=(1+1/rootc*(.115-1/rootc*(.186-1/rootc*(.217-.173/rootc))))/3
	bigXc=1-0.06*dH(1)	!Initial guess. This is estimated from Xc~0.7 for alcohols.
	alittle=9.5*q(1)*1.9+4*C(1)*1.7745-1.7745*1.9
	bq=1.7745*1.9*Zci+3*alittle
	BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)* &
		(-bq+SQRT(bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci) )
	Yc=Zci*Zci*Zci/alittle/(BcTest*BcTest)
	ZcOld=Zci
	! compute bVol,eokP that satisfies the critical point.  
	! if nParms.ge.3, then THIS IS SKIPPED
	dev=1234
	iter=0
	! first compute initial guess assuming Xc=const wrt eta
	! Note: this loop is superfluous if Xc=1.
	do while(abs(dev).gt.1.e-6.and.Nd(1).le.1)! this estimate fails for Nd>1
		iter=iter+1
		if(nParms.eq.1)VX(1)=BcTest*rGas*Tc(1)/Pc(1)
		bigBc=VX(1)*pc(1)/rGas/tc(1)
		etac=bigBc/Zci
		if(etac.gt.0.52)then
			if(LOUD)pause 'error in RegPureDev3 during ZcIter.  etac>0.52'
			etac=0.1
			EXIT !break the loop
		endif
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		sqarg=1+4*Nd(1)*alphac
        if(LOUD)then
	        if(sqarg.lt.0)pause 'error in RegPureDev.  sqarg<0'
        end if
	    if(sqarg.lt.0)EXIT
		XAc=2/(1+SQRT(sqarg))
		bigXc=1-Nd(1)*(1-XAc)
		Zci=(bigXc+(1.9-1.7745*Yc)*bigBc)/3
		bq=1.7745*1.9*Zci+3*alittle
		sqarg2=bq*bq+4*alittle*(4*C(1)-1.9)*(9.5*q(1)-1.7745*BigXc)/Zci
        if(LOUD)then
	        if(sqarg2.lt.0)pause 'sqarg2 < 0 in RegPureDev'
        end if
	    if(sqarg2.lt.0)EXIT
		BcTest=Zci*Zci/2/alittle/(4*C(1)-1.9)*(-bq+SQRT(sqarg2) )
		Yc=Zci*Zci*Zci/(alittle*BcTest*BcTest)
		dev=(Zci-ZciOld)/Zci
		ZciOld=Zci
		!bigXcN=3*Zci-(1.9-1.7745*Yc)*BcTest
		!dev=(bigXcN-bigXc)/bigXcN
		!bigXc=bigXcN
		!Yc=Zci*Zci*Zci/alittle/(bigBc*bigBc)
	enddo
	! iterate till dPdEta=d2PdEta2=0
	do iter=1,100  !do 100 weighted successive subst iterations whether you need it or not.
		!if(nParms.eq.1)VX(1)=bigBc*rGas*Tc/Pc
		!bigBc=VX(1)*pc/rGas/tc
		!etac=bigBc/Zci
		alphac=etac*kcStar(1)*Fhb/(1-1.9*etac)
		alphap=alphac/etac/(1-1.9*etac)
		alphapp=alphap*( 1/(etac*(1-1.9*etac))-1/etac+1.9/(1-1.9*etac) )
		sqarg=1+4*Nd(1)*alphac
	    !if(sqarg.lt.0)pause 'error in DevPure.  sqarg<0'
	    if(sqarg.lt.0)CYCLE
		XAc=2/(1+SQRT(sqarg))
		dXdEta= -XAc*XAc*Nd(1)*alphap/SQRT(sqarg)
		d2XdEta2=(-2*XAc*dXdEta*Nd(1)*alphap-XAc*XAc*Nd(1)*alphapp)/	&
			SQRT(sqarg)+2*(XAc*Nd(1)*alphap)**2/sqarg**1.5
		voidfrac=1-1.9*etac
		datt=1+1.7745*Yc*etac
		Zassoc= -Nd(1)*(1-XAc)/voidfrac
		Zci=1+4*C(1)*etac/voidfrac-9.5*q(1)*Yc*etac/datt+Zassoc
		dZdEta=4*C(1)/voidfrac**2-9.5*q(1)*Yc/datt**2+Nd(1)*dXdEta/	&
			voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2
		d2ZdEta2=8*C(1)*1.9/voidfrac**3+2*9.5*q(1)*Yc*1.7745*Yc/datt**3	&
			+Nd(1)*(((-2*1.9*1.9*(1-XAc)/voidfrac+2*1.9*dXdEta)/voidfrac	&
			+d2XdEta2)/voidfrac)
		dBdEta=Zci+etac*dZdEta
		d2BdEta2=2*dZdEta+etac*d2ZdEta2
		Yc=Yc*0.9+0.1*datt*datt/9.5/q(1)*( 4*C(1)/voidfrac**2+Zci/etac+	&
			Nd(1)*dXdEta/voidfrac-Nd(1)*(1-XAc)*1.9/voidfrac**2 )
		etac=etac*0.95-0.05*2*dZdEta/d2ZdEta2
        if(LOUD)then
		    if(etac.gt.0.52)pause 'error in RegPureDev, YcIter.  etac>0.52'
        end if
	ENDDO
	bigBc=2*Zci*Zci/d2ZdEta2/etac
	if(nParms==1)VX(1)=bigBc*rGas*Tc(1)/Pc(1)
	if(nParms < 3)eokP(1)=LOG(1.0617+Yc)*Tc(1)
	parm(2)=eokP(1)
	parm(3)=VX(1)
	pSat7=Pc(1)*10**( -1-acen(1) )
	T=Tc(1)*0.7D0
	call FugiESD(T,pSat7,X,nC,1,FUGC,rho,zFactor,Ares,uRes,iErrF)
	chemPoLiq=FUGC(1)
	call FugiESD(T,pSat7,X,nC,0,FUGC,rho,zFactor,Ares,uRes,iErrF)
	chemPoVap=FUGC(1)
	error(1)=chemPoLiq-chemPoVap
	if(iFlag==1.and.LOUD)print*,'RegPureDevEsd3: c,chemPoDev',C(1),error(1) !mostly to avoid iFlag warning.
	RETURN
	END

!***********************************************************************
      SUBROUTINE RegPureIo(NC)
	!C
	!C  PURPOSE:  Regress Pure component parameters of ESD model.
	!C  PROGRAMMED BY:  JRE 3/01
	!C
	USE GlobConst !NMX, avoNum,RGAS,.. TC,PC,...
	USE EsdParms
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	CHARACTER fileNameVp*55
	character outFile*251
	common/ PUREVP / VPA,VPB,VPC,VPD,VPE,TMINVP,TMAXVP
	common/ CONST /  rKadStar, dHkJ_mol, iDenOpt
	common/ PURELD /SG,rLDA,rLDB,rLDC,rLDD,TMINLD,TMAXLD
	common/ AVERR / AAPERR,RHOERR,objfun
	print*,'Minimize vp error based on estd rKadStar and dHkJ_mol '  
	print*,'constrained to match Tc. Pc is ~matched implicitly by vp' 
	print*,'or matched explicitly when nParms=1.'
	print*,'best setting is nParms=2 and iteropt=0 => flex to vx.'   
	print*,'kcStar comes from kcStar~kadStar/cShape so kadStar=.025'
	print*,'Data from c:\thermo\esd\VpCcEsdCorr.txt if available.'
	print*,'Enter dHkcal/mol,rKadStar(eg 0.025),iDenOpt(0=ignore SG)'
	read(*,*)dHkJ_mol,rKadStar,iDenOpt
	dHkJ_mol=dHkJ_mol*4.184
	!rKadStar=0.025
	!iDenOpt=0		!0=ignore density in fit, 1=weight sg as one datum
	!initial=0

	WRITE(*,*)'ENTER nParms,iteropt(0=AutoReg, 1=guess then reg)'
	READ(*,*)nParms0,iteropt
	idToReg=id(1)
	ndTemp=1
	iCompToReg=1
	if(nc.ne.1)then
		write(*,*)'Enter dippr id # and Nd of component for regression '
		read(*,*)idToReg,ndTemp
		do iComp=1,nc
			if(id(iComp).eq.idToReg)iCompToReg=iComp
		enddo
		nd(iCompToReg)=ndTemp
	endif
	fileNameVp='c:\thermo\esd\VpCcEsdCorr.txt'
	!write(*,*)'Enter filename for vpdata'
	!read(*,'(a)')fileNameVp
	iFound=0
	open(55,file=fileNameVp,ERR=86)
	read(55,*,END=86,ERR=86)nDbaseVp !error here means no file so skip cc search.
	write(*,*)'Enter chemcad id # of component for regression '
	read(*,*)idToRegCC
	do iVp=1,nDbaseVp
		READ(55,*,END=86)idDb,rmwDb,sgDb,vpADb,vpBDb,vpCDb,vpDDb,vpEDb,tMinVpDb,tMaxVpDb
		if(idDb.eq.idToRegCC)then
			iFound=1
			rmw=rmwDb
			sg=sgDb
			vpA=vpADb
			vpB=vpBDb
			vpC=vpCDb
			vpD=vpDDb
			vpE=vpEDb
			tMinVp=tMinVpDb
			tMaxVp=tMaxVpDb
		endif
	enddo
	CLOSE(55)
86	continue !if the dbase doesn't exist, then we must type the data in.
	if(iFound.eq.0)then
		write(*,*)'Compd not found in std dbase.'
		write(*,*)'Enter Mw,SpeGrav,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp'
		read(*,*)rMw,sg,vpA,vpB,vpC,vpD,vpE,tMinVp,tMaxVp
	endif
	
	vpA=vpA-LOG(1.E6)!convert dippr vp from Pa to MPa.

	TRLO=tMinVp/TC(iCompToReg)
	TRHI=tMaxVp/TC(iCompToReg)
	outFile=TRIM(masterDir)//'\output\RegPureOut.txt'
	open(66,file=outFile)
	IF(TRLO.GT.0.6.OR.TRHI.LT.0.8)THEN
	  if(LOUD)write(* ,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	  WRITE(66,'(i4,a,2f7.4)')ID(iCompToReg),' Warning VP RANGE SMALL. TrLo,Hi=',TRLO,TRHI
	END IF

	!  we need a big function here because of switching comp pos., and 
	!initial guesses to converge.
	call RegPureEsd(dHkJ_mol,rKadStar,nParms0,iterOpt,iDenOpt,id(iCompToReg),nC)!note:tc,pc,acen,id & esd parms passed in common
	close(66)
	if(LOUD)pause 'Results in RegPureOut.txt'
	return
	end
!***********************************************************************
	SUBROUTINE RegA0A1A2TptIo(NC)
	!
	! PURPOSE: Regress Tpt terms for Polymer Applications of the SPEADMD EOS.
	!Programmed by: AV 2/17/09
	!
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*44 FN
	PARAMETER (RGOLD=.61803399,CGOLD=.38196602)
	!character*77 errMsgLookup !errMsg(11),
	DoublePrecision deviate(maxPts),parm(5),stdErr(5)
	Integer idStore(NMX)
	EXTERNAL RegTptDevFcn
	character outFile*251
	parameter (nEta=5)
	common/SimData/ xFracc(23),vEffNm3(23),a0xss(20,20),A0Mix(20,20),A1Mix(20,20),A2Mix(20,20),Nps
	!dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	!common/ksvall/ks0ij,ks1ij
	!character junk*5
	!store current values of nc,id.

	ncStore=NC
	do i=1,NC
		idStore(i)=id(i)
	enddo
	zero=0
	eightySix= -86
	!WRITE(*,*)'ENTER NAME OF FILE WITH EXPERIMENTAL DATA (FN.FT in VpData subdir)'
	!WRITE(*,*)'FIRST LINE =NPTS,id 2ND LINE = TK, PMPA'
	!READ(*,'(A)')FN
	!inFile=TRIM(masterDir)//'\Mixes\'//TRIM(FN)
	!if(DEBUG)inFile='c:\spead\calceos\Mixes\'//TRIM(FN)
	!open(51,file=inFile)
	outFile=TRIM(masterDir)//'\Output\TptKs.txt'
	if(DEBUG)outFile='c:\spead\calceos\Output\TptKs.txt'
	open(61,file=outFile)

	!inFile=TRIM(masterDir)//'\MixSummary\'//MixSummary.txt
	!if(DEBUG)inFile='c:\spead\calceos\Mixes\MixSummary.txt'
	!open(51,file=inFile)
	!outFile=TRIM(masterDir)//'\Output\TptKs.txt'
	!if(DEBUG)outFile='c:\spead\calceos\Output\TptKs.txt'
	!open(61,file=outFile)
	nParms=1
	ier=0
	!do while(ier.eq.0) !loop over data sets in db
		do i=1,nParms
			parm(1)=0
			parm(2)=0
		enddo
		!ifound=0	
		!do i=1,1000
		!	READ(51,'(A5)')junk
		!	if(junk=='xFrac')then
		!		ifound=ifound+1
		!		if(ifound==2)exit
		!	endif
		!enddo
		!do i=1,1000
		!	Read(51,'(A500)') READTEXT
		!	READ(READTEXT,'(A4)')junk
		!	if(TRIM(junk)=='Ax,h')exit
		!	Read(readtext,'(4(F6.3,2x))')xFracc(i),wFrac(i),vFrac(i),vEffNm3(i)!,A1Xs2(i),A1Xs3(i),A1Xs4(i),A1Xs5(i),A2Xs1(i),A2Xs2(i),A2Xs3(i),A2Xs4(i),A2Xs5(i)
		!	read(readtext,'(34x,5(f7.3,5x))')A0Xs(i,1),A0Xs(i,2),A0Xs(i,3),A0Xs(i,4),A0Xs(i,5)
			!read(readtext,'(92x,e11.5,x,4(E11.5,x))')A1Xs(i,1),A1Xs(i,2),A1Xs(i,3),A1Xs(i,4),A1Xs(i,5)!,A2Xs1(i),A2Xs2(i),A2Xs3(i),A2Xs4(i),A2Xs5(i)
			!read(readtext,'(152x,e11.5,x,4(E11.5,x))')A2Xs(i,1),A2Xs(i,2),A2Xs(i,3),A2Xs(i,4),A2Xs(i,5)
		!enddo
	!	Call MixSummary(xFracc,vEffNm3,a0xss,A0Mix,A1Mix,A2Mix,Nps)
		MAXIT=1111
		INIT=1
		!C	LWA>5*nParms+NPTS ~ 275
		LWA=275
		!c      TOL=SQRT(DPMPAR(1))
		TOL=0.00001
		!C
		!c  factor controls the magnitude of the first step from the initial guess.
		factor=0.01
		!write(*,*)'KII,paadp,rmserr'
		nTot=(Nps-2)*nEta
		CALL LMDifEZ(RegTptDevFcn,nTot,nParms,parm,factor,deviate,TOL,iErrCode,stdErr)
		rmsErr=SQRT( ENORM(nTot,deviate) )
		if(LOUD)write(*,*)'iErrCode,rmsErr',iErrCode,rmsErr
		if(LOUD)write(* ,'(2f9.4)')(parm(i),i=1,nParms)
		if(LOUD)write(*,*)'StdErr'
		if(LOUD)write(* ,'(2f9.4)')(stdErr(i),i=1,nParms)
		if(LOUD)write(*,*) 'RegA0A1A2TptIo: Converged'
		L=1
        tKelvin=400 !assume this is most common condition.
		call KSVAL(nc,nParms,parm,tKelvin,deviate,KAADK,kBias,kErrMax,L)

		if(LOUD)write(* ,607)id(1),(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),KAADK,kBias,kErrMax,rmsErr
		write(61,607)id(1),(parm(i),i=1,nParms),(StdErr(i),i=1,nParms),KAADK,kBias,kErrMax,rmsErr
	!enddo !while(ier.eq.0) loop over data sets in db

	close(51)
	close(61)
	if(LOUD)pause 'These results are tabulated in TptKs.txt.'
	!restore stored nc,id
	NC=ncStore
	do i=1,NC
		id(i)=idStore(i)
	enddo
	!600	format('   T(K)', 4x, '  P(MPA)', 7x, 'x', 9x, 'y', 7x )
	606	FORMAT(1X,F8.4,1X,F10.2)      
	607	FORMAT(1X,i6,1X,<nParms>F10.4,1X,<nParms>F10.4,f10.2,f10.6,f10.2,f10.2,1x,a20)
	!608 FORMAT(1X,I6,1X,<nParms>F10.4,1X,A20)
  
	RETURN
	END
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine RegTptDevFcn(nTot,nParms,parm,deviate,IFLAG)
	implicit doublePrecision(A-H,K,O-Z)
	dimension parm(nParms),deviate(nTot)
	!C
	!C  SUBROUTINE FCN FOR LMDIFez EXAMPLE.
	nc=2
	last=0
	iflag=last
	!all we need is the deviate function gotten from kijval below.
    tKelvin=400 !assume this is most common condition JRE 20200406
	call KSVAL(nc,nParms,parm,tKelvin,deviate,KAADK,kBias,kErrMax,LAST)
	return
	end


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE KsVAL(nc,nParms,parm,tKelvin,deviate,kAADk,kBias,kErrMax,LAST)
	USE GlobConst
	USE BIPs
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	!CHARACTER*12 FN
	DOUBLE PRECISION parm(nParms),deviate(maxPts)
	parameter (nEta=5)
	common/SimData/ xFracc(23),vEffNm3(23),a0xss(20,20),A0Mix(20,20),A1Mix(20,20),A2Mix(20,20),Nps
	!common/XsProps/a0XsSs
	!dimension ks0ij(NMX,NMX),ks1ij(NMX,NMX)
	!common/ksvall/ks0ij,ks1ij
	DIMENSION XXFRAC(NMX)
            
	607	FORMAT(1X,F8.4,1X,F10.4,1X,F10.8,1X,F10.8,1X,F10.4,1X,i5)      

	!KS0IJ(1,2)=0
	!KS1IJ(1,2)=0
	!KS0IJ(2,1)=KS0IJ(1,2)
	!KS1IJ(2,1)=KS1IJ(1,2)
	!KS0IJ(1,1)=0
	!KS1IJ(1,1)=0

		if(last.lt.0)then !only kijOpt sets last to -1, so that is the only one that effects this.
			write(*,*)'Enter optimization choice: 1=ks0, 2=ks0,ks1'
			read(*,*)iOpt
		endif
		if (nParms==1)then
			KS0IJ(1,2)=parm(1)
		else
			KS0IJ(1,2)=parm(1)
			KS1IJ(1,2)=parm(2)
		endif
	MAXIT=1111
	KAADK=0
	ssqErr=0
	kBias=0
	eta=0.d0
	!tKelvin=273
	ikk=1
	nTot=nEta*(NPS-2)
	!Do ikk=1,nTot
	nFracs=NPS-2
		Do iEta=1,nEta
			eta=eta+0.1d0

			DO iData=1,nFracs
			!if(iData==(NPS-1))then
			!	xxFrac(1)=1
			!	xxFrac(2)=0
			!elseif(iData==NPS)then
			!	xxFrac(1)=0
			!	xxFrac(2)=1
			!else
				xxFrac(1)=xFracc(iData)
				xxFrac(2)=1-xxFrac(1)
			!endif
				IFLAG=0
				ITMAX=MAXIT
				INIT=1
				nComps=nc  
				CALL AxsCalc(xxFrac,eta,nComps,tKelvin,a0MixCalc,iErr)
				!a0XsCalc=a0MixCalc-(xxFrac(1)*A0Mix((NPS-1),iEta)+xxFrac(2)*A0Mix(NPS,iEta)) 
			!a0Xs=a0Mix-(xxFrac(1)*

				deviate(ikk) =(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0XsCalc-a0xss(iData,iEta))/a0xss(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100 !(a0MixCalc-A0Mix(iData,iEta))**2! !(a0MixCalc-A0Mix(iData,iEta))/A0Mix(iData,iEta)*100	 !
				if(ABS(deviate(ikk)).gt.ABS(kErrMax))kErrMax=deviate(ikk)
				
				!penalize the deviate passed to lmdif, but only compile paadp for converged pts
				KAADK = KAADK+ABS(deviate(ikk))
				ssqErr= ssqErr+deviate(ikk)*deviate(ikk)
				kBias=kBias+deviate(ikk)
		!	if(ABS(deviate(iData)).gt.ABS(pErrMax))pErrMax=deviate(iData)
			
				!WRITE(67,607)tDat(iData),P,X(1),Y(1),deviate(iData),ier(1)

				ikk=ikk+1
			enddo
		enddo
	!enddo

	rmsErr=SQRT(ssqErr/(nEta*(NPS-2)))
	KAADK=KAADK/(nEta*(NPS-2))
	kBias=kBias/(nEta*(NPS-2))
	if(LOUD)write(*,'(7f10.3)')(parm(i),i=1,nparms),KAADK,rmsErr

	RETURN
	END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	


	SUBROUTINE CP_Mixture(isZiter,toll,gmol,NC,TC_mix,VC_mix,PC_mix,ZC_mix,iErrCode)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
!C																														   C
!C Written by AGF, Dec. 2009																							   C
!C It calculates the critical point for a given mixture according to Tangent Plane Distance Theory described by Michelsen- C
!C in his book, Thermodynamic models: Fundamentals and computational aspects.											   C
!C Input variables would be number of components and their mole numbers and output variables are TC, PC, VC and ZC 		   C
!C isZiter is a flag fuVtot: 0 means compute all derivatives and fugacity. 1 means just compute ZAU at given T,rho.																					   C
!C																														   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	USE GlobConst
	USE SpeadParms
	IMPLICIT DOUBLEPRECISION(A-H,K,L,O-Z)
    PARAMETER(nV=6)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),Rec_Jac(NMX,NMX)
	DIMENSION VC(NMX),C_V(NMX,NMX),S_V(NMX,NMX),eta_Vc(NMX,NMX),phi_Vc(NMX),phi_Tc(NMX)
	DIMENSION sum_guess(nV),vTotCc_guess(nV),GJ_Jac(2,2),B_GJ(2,1)
	character*200 errMsg(4)
	common/iterC/iter
	COMMON/ETA2/ETA
	COMMON/NC_num_rec/NC_tqli
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)

	ErrMsg(1)='CP_Mixture Error - The determinant of the Jacobian matrix is equal to 0!'
	ErrMsg(2)='CP_Mixture error - Please use "CR" option with EOS # 4, 5 & 9'
	ErrMsg(3)='CP_Mixture warning - Number of iterations exceeded 500 while typically it should be around 10'
	iErrCode=0

!Calculate mole fractions:
	totMoles=0.d0 
	do i=1,nc
		totMoles=totMoles+gmol(i)  
	enddo
	totMoles_old=totMoles
	do i=1,nc
		xFrac(i)=gmol(i)/totMoles
	enddo
	two=2.d0

!Step size for change in T and V:
	if (iEosOpt.eq.4) then
		stepT=1.d6
		stepV=1.d6
	elseif (iEosOpt.eq.5.or.iEosOpt.eq.9) then
		stepT=1.d6
		stepV=1.d6
	else
	 	write(*,*) ErrMsg(2)
		iErrCode=2
		if(LOUD)pause 
		return
	endif

!Initial guess for VC. Six different values are guessed for VC and the best one is picked
	C_Vc=0.1559d0
	sum_Vc=0.d0
	DO I=1,NC
		IF (iEosOpt.eq.4) THEN
			VC(I)=Zc(I)*rGas*TC(I)/PC(I) 
		ELSEIF (iEosOpt.eq.5.or.iEosOpt.eq.9) THEN
  			VC(I)=vMolecNm3(I)*avoNum/0.23d0
		ENDIF
		sum_Vc=sum_Vc+xFrac(I)*(VC(I)**(two/3.d0))
	ENDDO
	DO I=1,NC
		phi_Vc(I)=xFrac(I)*VC(I)**(two/3.d0)/sum_Vc
		DO J=1,NC
			IF (I.ne.J) THEN
				eta_Vc(I,J)=(VC(I)-VC(J))/(VC(I)+VC(J))
				C_V(I,J)=-1.4684d0*eta_Vc(I,J)+C_Vc
				S_V(I,J)=C_V(I,J)*(VC(I)+VC(J))/two
			ELSE
				eta_Vc(I,J)=0.d0
				C_V(I,J)=0.d0
				S_V(I,J)=0.d0
			ENDIF
		ENDDO
	ENDDO
	sum_vc2=0.d0
	sum_vc3=0.d0
	DO I=1,NC
		sum_vc2=sum_vc2+phi_Vc(I)*VC(I)
		DO J=I+1,NC
			IF (I.ne.J) THEN
				sum_Vc3=sum_Vc3+phi_Vc(I)*phi_Vc(J)*S_V(I,J)
			ENDIF
		ENDDO
	ENDDO
	vTotCc_guess(3)=(sum_vc2+sum_vc3)*totMoles					  !Initial guess for V, using volume fraction
	vTotCc_guess(1)=vTotCc_guess(3)-0.1d0*vTotCc_guess(3)
	vTotCc_guess(2)=vTotCc_guess(3)-0.2d0*vTotCc_guess(3)
	vTotCc_guess(4)=vTotCc_guess(3)+0.1d0*vTotCc_guess(3)
	vTotCc_guess(5)=vTotCc_guess(3)+0.2d0*vTotCc_guess(3)
	vTotCc_guess(6)=vTotCc_guess(3)+0.3d0*vTotCc_guess(3)

!Initial guess for TC
	sum_Tc=0.d0
	DO I=1,NC
		sum_Tc=sum_Tc+VC(I)*xFrac(I)
	ENDDO
	TCt=0
	DO I=1,NC
		phi_Tc(I)=VC(I)*xFrac(I)/sum_Tc
		TCt=TCt+phi_Tc(I)*TC(I)
	ENDDO
	IF (iEosOpt.eq.4) THEN
		tKelvin=TCt
	ELSEIF (iEosOpt.eq.5.or.iEosOpt.eq.9) THEN
		tKelvin=TCt+0.05d0*TCt
	ENDIF
	T_old=tKelvin

!This part tests which initial guess for VC is the best:
	DO m=1,nV
		vTotCc=vTotCc_guess(m)
		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_min,dLambda_ds,iErrCode)
		sum_guess(m)=abs(dLambda_ds)+abs(Lambda_min)
	ENDDO
	guess_final=MIN(abs(sum_guess(1)),abs(sum_guess(2)),abs(sum_guess(3)),abs(sum_guess(4)),abs(sum_guess(5)),abs(sum_guess(6)))
	DO I=1,nV
		IF (guess_final.eq.abs(sum_guess(I))) THEN
			vTotCc=vTotCc_guess(I)*1.1d0		! Initial guess for V
			vTotCc_old=vTotCc
		ENDIF
	ENDDO

!----------------------------main iteration--------------------------
	OF=1.d0
	iter=0
	DO WHILE (OF.gt.toll)
		iter=iter+1
 		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_min,dLambda_ds,iErrCode)

		tKelvin=tKelvin+tKelvin/stepT
		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_minUpperT,dLambda_dsUpperT,iErrCode)
		tKelvin=T_old-T_old/stepT
		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_minLowerT,dLambda_dsLowerT,iErrCode)
		dLambda_dT=(Lambda_minUpperT-Lambda_minLowerT)/(two*T_old/stepT)
		dC_dT=(dLambda_dsUpperT-dLambda_dsLowerT)/(two*T_old/stepT)	 !The quantity dlambda/ds is called C in Michelsen's book...

		VtotCc=VtotCc+VtotCc/stepV
		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_minUpperV,dLambda_dsUpperV,iErrCode)
		VtotCc=VtotCc_old-VtotCc_old/stepV
		call numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_minLowerV,dLambda_dsLowerV,iErrCode)
		dLambda_dV=(Lambda_minUpperV-Lambda_minLowerV)/(two*vTotCc_old/stepV)
		dC_dV=(dLambda_dsUpperV-dLambda_dsLowerV)/(two*vTotCc_old/stepV)

		IF (NC.eq.2) THEN
		   	Det_Jac=dLambda_dT*dC_dV-dLambda_dV*dC_dT
			IF (Det_Jac.eq.0) then
				iErrCode=3
				write(*,*) ErrMsg(1)
				write(*,*)'iter= ',iter
				if(LOUD)pause 
				return
			endif
			Rec_Jac(1,1)=dC_dV/Det_Jac
			Rec_Jac(1,2)=-dLambda_dV/Det_Jac
			Rec_Jac(2,2)=dLambda_dT/Det_Jac
			Rec_Jac(2,1)=-dC_dT/Det_Jac
			DeltaT=-Rec_Jac(1,1)*Lambda_min-Rec_Jac(1,2)*dLambda_ds
			DeltaV=-Rec_Jac(2,1)*Lambda_min-Rec_Jac(2,2)*dLambda_ds
		ELSE
			GJ_Jac(1,1)=dLambda_dT
			GJ_Jac(1,2)=dLambda_dV	
			GJ_Jac(2,1)=dC_dT
			GJ_Jac(2,2)=dC_dV
			B_GJ(1,1)=-Lambda_min
			B_GJ(2,1)=-dLambda_ds
			call gaussj(GJ_Jac,2,2,B_GJ,1,1)
			DeltaT=B_GJ(1,1)
			DeltaV=B_GJ(2,1)
		ENDIF
	!Step limiting for V and T:
		if (abs(DeltaV).GT.(vTotCc/10.d0)) then
			DeltaV=DeltaV/(abs(DeltaV))*vTotCc/10.d0
		endif
		if (abs(DeltaT).GT.(tKelvin/10.d0)) then
			DeltaT=DeltaT/(abs(DeltaT))*tKelvin/10.d0
		endif
		OF=Lambda_min**2+dLambda_ds**2		
		tKelvin=tKelvin+DeltaT
		vTotCc=vTotCc+DeltaV*totMoles
		T_old=tKelvin
		vTotCc_old=vTotCc	
		if (iter.GT.100) then
			write(*,*) ' ' 
			write(*,*) ErrMsg(3)
			write(*,*) ' '
			write(*,*) 'The toll variable is too small: toll=',toll,' or'
			write(*,*) 'the system is not well-behaved in the critical region'
			write(*,*) 'If dT=',DeltaT,' and dV=',DeltaV
			write(*,*) 'are small enough for you particular application, you can accept the results'
			write(*,*) 'Otherwise, ignore the CriticalPoint.txt file.' 
			write(*,*) ' ' 
			TC_mix=tKelvin
		    VC_mix=vTotCc/totMoles
			call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
			ZC_mix=Z
			PC_mix=PMpa
			etaC_mix=eta
			return
		endif
	enddo
	TC_mix=tKelvin
	VC_mix=vTotCc/totMoles
	call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
	ZC_mix=Z
	PC_mix=PMpa
	etaC_mix=eta

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							  C
!C	Written by AGF, Dec. 2009																  C
!C	calculation of the critical point for pure components	using ESD, SPEAD and SPEAD11 EOSs C 
!C	There is no input variables, Output variables are TC, PC, ZC and VC 					  C
!C																							  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          !call CritPure(NC,isZiter,toll,TC_Pure,VC_PURE,PC_Pure,ZC_Pure,Acen_pure,iErrCode)
	SUBROUTINE CritPure(NC,isZiter,toll,TC_Pure,VC_Pure,PC_Pure,ZC_Pure,Acen_pure,iErrCode)
	USE GlobConst
	USE BIPs
    IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(N=2)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),VC(NMX) !,ier(11)
	!DIMENSION GJ_Jac(2,2),B_GJ(2,1)
	!DIMENSION KIJ(NMX,NMX),KTIJ(NMX,NMX),HIJ(NMX,NMX),HTIJ(NMX,NMX),xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)
	doubleprecision mShape(NMX),Acen_pure
	character*200 ErrMsg(4) !,outfile*100
	LOGICAL LOUDER
	!COMMON/ETA2/ETA
	!COMMON/BIPs_SPEAD/aBipAD,aBipDA
	!COMMON/BIPs/KIJ,KTIJ,HIJ,HTIJ,xsTau,xsTauT,xsAlpha	 !warning: dim of aBip should be maxTypes, not nmx.
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,iFlagAssoc
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)
	!COMMON/Assoc/eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes),bondVolNm3T(nmx,maxTypes),nDegree(nmx,maxTypes),nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),idType(nmx,maxTypes),localType(maxTypesGlobal),idLocalType(maxTypes),nTypes(NMX),nTypesTot
	common/rg/iFlagRG,iFlagFit,mShape,cutOff(NMX),Phi(NMX),Zi(NMX)
	data initCall/1/
	LOUDER=LOUD
	LOUDER=.TRUE. ! PROVIDES PROSPECT OF LOCAL CONTROL OF CONSOLE FEEDBACK.
	iErrCode=0
	ErrMsg(1)='CrPure error - NC>1 , NC must be equal to 1'
	ErrMsg(2)='CrPure error - Please use "CP" option with EOS # 4, 5 & 9'
	ErrMsg(3)='CrPure warning - Number of iterations exceeded 100 while typically it should be around 7'

	if (NC.ne.1) then
		if(LOUDER)write (*,*) ErrMsg(1)
		iErrCode=1
		if(LOUDER)pause
		return
	endif
	two=2.d0

! Specifying the step size for the iteration on T & V 
		stepT=1.d5
		stepV=1.d5
	IF (iEosOpt.eq.5.or.iEosOpt.eq.9) THEN
        stepT=10/toll
		stepV=10/toll
	ENDIF
	NC=1

! Initial guesses for TC and VC
	VcExp=Zc(1)*rGas*TC(1)/PC(1)
    xFrac(1)=1
    gmol(1)=1
	VC(1)=Zc(1)*rGas*TC(1)/PC(1)
	tKelvin=TC(1)
	IF (iEosOpt==5 .or. iEosOpt==9) THEN
		VC(1)=bVolCC_mol(1)/0.11d0
		tKelvin=TC(1)*1.5 !add a littl to initial guess because TPT2 overestimates.
	ENDIF
	vTotCc=VC(1)

	iErrFlag=1
!	if (iEosOpt.EQ.9.or.iEosOpt.EQ.5.or.iEosOpt.Eq.4.or.iFlagAssoc.Eq.1) iErrFlag=0
!	if(iErrFlag==1)then
!		write(*,*) ErrMsg(2)
!		iErrCode=2
!		if(LOUDER)pause 
!		return
!    endif
!----------------------------main iteration--------------------------
!Basically in the iteration we are minimizing dP/dV and d2P/dV2
!usually the iteration converges after 10 cycles.
	!-------------------------------------------------------------------------------
	!bRef=bVolRef(1,tKelvin)
    !eta=bVolCC_mol(1)/vTotCc
    !tKelvin=5.2             !for debugging
    !eta=0.40469             !for debugging
    if(LOUDER.and.initCall)write(*,'(a,f7.2,a,f7.4,a)')' Default guess is: Tc=',tKelvin,' etac=',eta,' . Enter guess for Tc,etac'
    !read*,tKelvin,eta      !for debugging
    !vTotCc=bVolCC_mol(1)/eta !for debugging
	rho=1/vTotCc
    rhoHi=(1+1/stepV)*rho
    rhoLo=(1-1/stepV)*rho
    isZiter=1               !Simplified to pure numerical derivatives. JRE 20200406
	call FuVtot(isZiter,tKelvin,1/rho,gmol,NC,FUGC,Zbase,aRes,uRes,iErrVtot)
    if(iErrVtot.ne.0)iErrCode=4
    if(iErrCode.ne.0)goto 86
	call FuVtot(isZiter,tKelvin,1/rhoHi,gmol,NC,FUGC,ZHi,aRes,uRes,iErrVtot)
	call FuVtot(isZiter,tKelvin,1/rhoLo,gmol,NC,FUGC,ZLo,aRes,uRes,iErrVtot)
    pBaseJre=zBase*rho*rGas*tKelvin
    pHiJre=zHi*rhoHi*rGas*tKelvin
    pLoJre=zLo*rhoLo*rGas*tKelvin
    dP_dRho=(pHiJre-pLoJre)/(rhoHi-rhoLo)
    d2P_dRho2=(pHiJre+pLoJre-2*pBaseJre)/(rhoHi-rho)**2
	Told=tKelvin
	dPold=dP_dRho
        
	!  P   =  Z*R*T/V
	!dP_dV = R*T/V*dZ_dV - Z*R*T/V^2 = (R*T/V)*(dZ_dV - Z/V)
	!dZ_dV=(ZHi-ZLo)/(2*vTotcc/stepV) 
	!dP_dVcheck = rGas*tKelvin/vTotCc*(dZ_dV-Zbase/vTotcc)
	!d2Z_dV2= (zHi+zLo-2*zBase)/(2*vTotCc/stepV)
	!!d2P_dV2 = R*T*( d2Z/dV2/V - 2*dZ_dV/V^2 + 2*Z/V^3)
	!d2P_dV2=rGas*tKelvin*( d2Z_dV2/vTotCc-2*dZ_dV/vTotCc**2+2*zBase/vTotCc**3 ) 
	!!d2P/dRho2 = -d2P_dV2*vTotCc**2 - dP_dV*(-2/rho**3) = 
	!dP_dRho= -dP_dV*vTotCc**2
	!d2P_dRho2= -d2P_dV2*vTotCc*vTotCc -2*dP_dRho*vTotCc 
	!print*,'dp,dpcheck,tol',dP_dV,dP_dVcheck,toll
	if(LOUDER.and.initCall)write(*,'(f7.2,3F7.4,2e12.4)')tKelvin,eta,PMPa,vTotCc/1000,dP_dRho,d2P_dRho2
	!if(LOUDER)pause 'check initial guess'
	!We calculate d2P/dEta2 numerically, including 5&9, because GaussFun makes 5&9 too complicated. 	
	!tKelvin=tKelvin*1.01
	 if(LOUDER.and.initCall)print*,'bVol,Mw',bVolCC_mol(1),rMw(1)
	Delta_T=33
    rmsBest=1234
    tBest=tKelvin
    rhoBest=rho
	rhoOld=rho
	rho=rho*1.01
	DO iter=1,66 
	!plotting shows that P vs. rho is simple but P vs. V goes concave up, then down, then up.
		d2Pold=d2P_dRho2
		do iterD2P=1,44 !Iterate on rho to find d2P_dRho2=0
            vTotcc=1/rho
            rhoHi=(1+1/stepV)*rho
            rhoLo=(1-1/stepV)*rho
	        call FuVtot(isZiter,tKelvin,1/rho,gmol,NC,FUGC,Zbase,aRes,uRes,iErrVtot)
		    call FuVtot(isZiter,tKelvin,1/rhoHi,gmol,NC,FUGC,ZHi,aRes,uRes,iErrVtot)
		    call FuVtot(isZiter,tKelvin,1/rhoLo,gmol,NC,FUGC,ZLo,aRes,uRes,iErrVtot)
            if(iErrVtot.ne.0)iErrCode=4
            pBaseJre=zBase*rho*rGas*tKelvin
            pHiJre=zHi*rhoHi*rGas*tKelvin
            pLoJre=zLo*rhoLo*rGas*tKelvin
            dP_dRho=(pHiJre-pLoJre)/(rhoHi-rhoLo)
            d2P_dRho2=(pHiJre+pLoJre-2*pBaseJre)/(rhoHi-rho)**2
            if(ABS(d2P_dRho2-d2Pold) < 1E-33)then
                print*,'rhoHi,pHi',rhoHi,pHiJre
                print*,'rho  ,pMi',rho,pBaseJre
                print*,'rhoLo,pLo',rhoLo,pLoJre
                if(LOUDER)pause 'CritPure: (d2P_dRho2-d2Pold) ~0'
                d2P_dRho2=d2Pold+0.1
                stepV=stepV/2 !this error seems to be a problem with roundoff in estimating d2P
            end if
			Delta_Rho=d2P_dRho2/(d2P_dRho2-d2Pold)*(rho-rhoOld)
			if (abs(Delta_Rho) > rho/10)Delta_Rho=Delta_Rho/(abs(Delta_Rho))*rho/10 !Step limiting for V:
            rms=SQRT(dP_dRho*dP_dRho+d2P_dRho2*d2P_dRho2)/1.414
            if(rms < rmsBest)then
                tBest=tKelvin
                rhoBest=rho
                rmsBest=rms
            endif
		    eta=bVolCC_mol(1)/vTotCc
			if(LOUDER.and.initCall)write(*,'(f7.2,3F7.4,3e12.4)')tKelvin,eta,PMPa,rho,dP_dRho,d2P_dRho2 !,rms
			d2Pold=d2P_dRho2
			rhoOld=rho
			rho=rho-Delta_Rho
			if(ABS(Delta_Rho/rho) < 1/stepV .and. ABS(d2P_dRho2) < 0.01)exit !exits do loop.
            if(iterD2P==1)then !this helps if the next guess for T is too similar to the last.
                d2Pold=d2P_dRho2
                rhoOld=rho
                rho=rho*1.001
                cycle
            endif
		enddo
		!if(LOUDER)pause 'check d2P'
        if(iter==1)then  !Establish initial point where d2P=0 but dP is not. Future iterations will be along this inflection.
            dPold=dP_dRho
            tOld=tKelvin
            tKelvin=tOld*1.001
            cycle
        endif
		!vTotCc=1/rho
		!call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Zbase,iErrVtot)
        !if(iErrVtot.ne.0)iErrCode=6
		!dP_dRho= -dP_dV*vTotCc**2 !dP_dV is passed by common/derv1/.
		Delta_T=dP_dRho/(dP_dRho-dPold)*(tKelvin-Told)
		if (abs(Delta_T) > tKelvin/10) Delta_T=Delta_T/(abs(Delta_T))*tKelvin/10 !Step limiting for T:
		dPold=dP_dRho
		Told=tKelvin
		OF=dsqrt(dP_dV*dP_dV)
		tKelvin=tKelvin-Delta_T
		rhoG_cc=rMw(1)/vTotCc
        rms=SQRT(dP_dRho*dP_dRho+d2P_dRho2*d2P_dRho2)/1.414
        if(rms < rmsBest)then
            tBest=tKelvin
            rhoBest=rho
            rmsBest=rms
        endif
		if(LOUDER.and.initCall)write(*,'(f9.4,3F7.4,3e12.4)')tKelvin,eta,PMPa,rhoG_cc,dP_dRho,d2P_dRho2,rms
		IF (iter.GT.22 .or. iErrCode.ne.0 .and. LOUD) THEN
			write(*,*) ' ' 
			write(*,*) 'CritPure: iErrCode= ',iErrCode,' iter = ',iter
			write(*,*) ' ' 
			write(*,*) 'FuVtot error or the tol variable is too small: toll=',toll
			write(*,*) 'or the system is not well-behaved in the critical region'
			write(*,*) 'If dT=',Delta_T,' and dRho=',Delta_Rho
			write(*,*) 'are small enough for your particular application, you can accept the results'
			write(*,*) 'Otherwise, ignore the CriticalPoint.txt file.' 
			write(*,*) 'Best values so far: Tc= ',tBest,'        etac = ',rhoBest*bVolCC_mol(1),'rmsErr=',rmsBest 
			TC_Pure=tBest
			print*,'Enter manual guess for Tc,etac (<0 to terminate)'
			read*,tKelvin,etac
			rho=etac/bVolCC_mol(1)
			if(tKelvin > 0)cycle
			tKelvin=TC_Pure
            rho=rhoBest
			exit
		ENDIF
		if(ABS(Delta_T/tKelvin) < toll)exit !terminate iteration
	ENDDO
    TC_Pure=tBest
    vTotCc=1/rho
	VC_PURE=vTotCc
	call FuVtot(isZiter,tBest,1/rhoBest,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
	ZC_PURE=Z
	PC_PURE=Z*(rhoBest*rGas*tBest)
	etaC_PURE=rhoBest*bVolCC_Mol(1)
	rhoc=rMw(1)*rhoBest
	if(LOUDER)write(*,'(3(a,f7.3))')' Tc=',TC_PURE,'      Pc=',PC_PURE,'     etaC=',etaC_PURE
    if(LOUDER.and.initCall)print*,'CP done. Calling Psat for acen.'
	T7=0.7D0*tBest
	call PsatEar(T7,pMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV,ierCodePsat)
	aceFactor=86  !NOTE: acen() is a vector member of GlobConst
    if(ierCodePsat==0)aceFactor= -log10(pMPa/Pc_Pure) - 1
	acen_pure=aceFactor

	if(LOUDER)write(*,'(5(a,f7.4))')' Zc=',ZC_PURE,'    acen= ',aceFactor,'  rhoc(g/cc)=',rhoc,' Psat7=',pMPa,' T7=',T7
	
	initCall=0
86	RETURN

	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							 C
!C Programmed by AFG 2011																	 C
!C This is a utility subroutine that returns the eigen vectors and the minimum eigen value	 C
!C of the Jacobian matrix of fugacity. For more details refer to Michelsen's book. 			 C
!C																							 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE eigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,eigenVec,LambdaMin,iErrCode)
	USE SpeadParms
	IMPLICIT DOUBLEPRECISION(A-H,K,L,O-Z)
	!PARAMETER(NMX=55)
	CHARACTER*77 errMsg(2)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),Bij(NMX,NMX)
	DIMENSION Lambda(NMX),eigen_vec_mat(NMX),eigenVal(NMX),eigenVec(NMX)
	DIMENSION a_tred2(NC,NC)
	common/iterC/iter
	COMMON/NC_num_rec/NC_tqli
	!COMMON/TptParms/zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX),nTptCoeffs
	COMMON/fugCR/PMpa,dFUG_dN(NMX,NMX),dP_dN(NMX)

	ErrMsg(1)='error eigen - v_norm is equal to zero. Check the eigen vectors...'
	iErrCode=0
	half=1.d0/2.d0
	call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
	DO I=1,NC
		DO J=1,NC
			Bij(I,J)=dFUG_dN(I,J)*sqrt(xFrac(I)*xFrac(J))
		ENDDO
	ENDDO
	IF (NC.eq.2) THEN
		Lambda(1)=half*( Bij(1,1)+Bij(2,2)+sqrt(4.d0*Bij(1,2)*Bij(2,1)+(Bij(1,1)-Bij(2,2))**2)) 
		Lambda(2)=half*( Bij(1,1)+Bij(2,2)-sqrt(4.d0*Bij(1,2)*Bij(2,1)+(Bij(1,1)-Bij(2,2))**2)) 
		LambdaMin=MIN(Lambda(1),Lambda(2))
		eigenVec(2)=1.d0
		eigenVec(1)=Bij(1,2)/(LambdaMin-Bij(1,1))
		v_norm=sqrt(eigenVec(2)**2+eigenVec(1)**2)
		if (v_norm.ne.0.d0) then
			eigenVec(2)=eigenVec(2)/v_norm
			eigenVec(1)=eigenVec(1)/v_norm
		else
			iErrCode=1 
			write(*,*) ErrMsg(1)
			!pause
			return
		endif
	ELSE
		DO I=1,NC
			DO J=1,NC
				a_tred2(I,J)=Bij(I,J)
			ENDDO
		ENDDO
		call tred2(a_tred2,NC,NC,eigenVal,eigen_vec_mat)
		call tqli(eigenVal,eigen_vec_mat,NC,NC,a_tred2)
		call eigsrt(eigenVal,a_tred2,NC,NC)
		LambdaMin=eigenVal(NC)
		DO I=1,NC
			eigenVec(I)=abs(a_tred2(NC,I))
		ENDDO
	ENDIF

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							  C
!C Programmed by AFG 2011																	  C
!C This is a utility subroutine that returns the numerical derivative of lambdaMin calculated C
!C by eigen subroutine in respect to mole fraction.											  C
!C For more details refer to Michelsen's book. 												  C
!C																							  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE numDiffEigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,Lambda_min,dLambda_ds,iErrCode)
	USE GlobConst
    IMPLICIT DOUBLEPRECISION(A-H,K,L,O-Z)
	character*77 errMsg(2)
	DIMENSION gmol(NMX),gmol_old(NMX),xFrac(NMX)
	DIMENSION eigenVec(NMX),eigenVecUp(NMX)

	ErrMsg(1)='error numDiffEigen - totMoles is equal to zero!!'
	iErrCode=0
!Step size for changing gmol	
	two=2.d0
	IF (iEosOpt.eq.4) THEN
		stepEps=1.d5 
	ELSEIF (iEosOpt.eq.5.or.iEosOpt.eq.9) THEN
		stepEps=1.d4 
	ENDIF
	epsilon=1.d0/stepEps

!Calculate lambda (the minimum eigen value of fugacity Jacobian matrix at gmol) 
	if (iErrCode.eq.0) then
		call eigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,eigenVec,LambdaMin,iErrCode)
	else
		return
	endif
	Lambda_min=LambdaMin
	totMolesOld=0.d0
	DO J=1,NC
		totMolesOld=totMolesOld+gmol(J)
	ENDDO

!Calculate lambda at gmol+dg 
	totMoles=0.d0
	DO J=1,NC
		gmol_old(J)=gmol(J)
		gmol(J)=gmol(J)+epsilon*eigenVec(J)*sqrt(gmol(J))
		totMoles=totMoles+gmol(J)
	ENDDO
	do i=1,nc
		if (totMoles.ne.0.d0) then
			xFrac(i)=gmol(i)/totMoles
		else
			iErrCode=2
			write(*,*)ErrMsg(1)
		endif
	enddo
	if (iErrCode.eq.0) then
		call eigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,eigenVecUp,LambdaMin,iErrCode)
	else 
		return
	endif
	Lambda_min_upper=LambdaMin

!Calculate lambda at gmol-dg 
	totMoles=0.d0
	DO J=1,NC
		gmol(J)=gmol_old(J)-epsilon*eigenVec(J)*sqrt(gmol_old(J))
		totMoles=totMoles+gmol(J)
	ENDDO
	do i=1,nc
		if (totMoles.ne.0.d0) then
			xFrac(i)=gmol(i)/totMoles
		else
			iErrCode=2
			write(*,*)ErrMsg(1)
		endif
	enddo
	if (iErrCode.eq.0) then
		call eigen(isZiter,tKelvin,vTotCc,gmol,xFrac,NC,eigenVec,LambdaMin,iErrCode)
	else
		return
	endif
	Lambda_min_lower=LambdaMin

!First derivative of lambda in respect to gmol
	dLambda_ds=(Lambda_min_upper-Lambda_min_lower)/(two*epsilon)		 !is called C in Michelsen's book.
	DO I=1,NC
		gmol(I)=gmol_old(I)
	ENDDO
	totMoles=totMolesOld
	do i=1,nc
		if (totMoles.ne.0.d0) then
			xFrac(i)=gmol(i)/totMoles
		else
			iErrCode=2
			write(*,*)ErrMsg(1)
		endif
	enddo

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C																							 C
!C Programmed by AFG 2011																	 C
!C This routine takes the first derivative of x in respect to T and/or V numerically.		 C
!C Setting stepT=0 means dxdT is not needed. 												 C
!C Setting stepV=0 means dxdV is not needed.												 C
!C Please refer to the code to find the definitions of stepT and stepV. 					 C
!C A typical value of 1.d5 for stepT and stepV is satisfying for most cases.				 C
!C User has to make sure that the quantity "x" is returned properly by FuVtot program. 		 C
!C																							 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   By setting dP_dV in the "x" position, the code knows to operate on dP_dV wherever "x" appears. (JRE, 5/27/2018)
! 		call   numDiff(isZiter,tKelvin,vTotCc,NC,gmol,stepT,0.d0,dP_dV,d2P_dVdT,dxdV)
    subroutine numDiff(isZiter,tKelvin,vTotCc,NC,gmol,stepT,stepV,x,dxdT,dxdV,iErr)
    IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
    PARAMETER(NMX=55)
	character*123 ErrMsg(2)
	DIMENSION gmol(NMX),FUGC(NMX)
	COMMON/Derv1/DFUG_DN_NUM(NMX,NMX),dP_dV,d2P_dV2,d3P_dV3,iFlagAssoc

	ErrMsg(1)='numDiff warning - the quantity "x" is not changing by tKelvin. Make sure that x is returned properly by FuVtot'
	ErrMsg(2)='numDiff warning - the quantity "x" is not changing by VtotCc. Make sure that x is returned properly by FuVtot'
    iErr=0
	two=2.d0
	T_old=tKelvin
	vTotCc_old=vTotCc	
    if (stepT.gt.1.d0) then 
		tKelvin=tKelvin+tKelvin/stepT
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
        if(iErrVtot.ne.0)iErr=1
		x_higherT=x

		tKelvin=T_old-T_old/stepT
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
        if(iErrVtot.ne.0)iErr=iErr+10
        if(iErr.ne.0)return
		x_lowerT=x
		tKelvin=T_old
	
		if (x_higherT.eq.x_lowerT) write(*,*) ErrMsg(1)
		dxdT=(x_higherT-x_lowerT)/(two*tKelvin/stepT)
	else
		dxdT=0.d0 !this is OK when all you want is dxdV
		!pause 'numdiff: stepT < 1???'
	endif
	if (stepV.gt.1) then
		vTotCc=vTotCc+vTotCc/stepV
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
        if(iErrVtot.ne.0)iErr=iErr+100
		x_higherV=x   
		
		vTotCc=vTotCc_old-vTotCc_old/stepV
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,aRes,uRes,iErrVtot)
        if(iErrVtot.ne.0)iErr=iErr+1000
        if(iErr.ne.0)return
		x_lowerV=x
		vTotCc=vTotCc_old
		
		if (x_higherV.eq.x_lowerV) write(*,*) ErrMsg(2)
		dxdV=(x_higherV-x_lowerV)/(two*vTotCc/stepV)
	else
		dxdV=0.d0  !this is OK when all you want is dxdT.
		!pause 'numdiff: stepV < 1???'
	endif
	tKelvin=T_old
	vTotCc=vTotCc_old

    return 
    end

