MODULE Esd2PolarFit
	! Regression utilities for the corrected ESD2-Marshall workflow.
	! Pure saturation, VE, and VLE are deliberately separated by parameter.
	USE GlobConst, only:nmx,zeroTol
	IMPLICIT NONE
	PRIVATE
	PUBLIC SetEsd2PolarMode,FitEsd2PureAtAlpha,CalcEsd2ExcessVolume
	PUBLIC GetEsd2CriticalState,PredictEsd2CriticalFixedParameters
	PUBLIC FitEsd2BetaVE,ScoreEsd2VLE,FitEsd2KijVLE
	PUBLIC FitEsd2BetaKijNested
	PUBLIC FitEsd2AlphaBinaryNested
	PUBLIC FitEsd2AlphaFamilyNested,FitEsd2NonpolarFamily
	PUBLIC FitEsd2HybridFamily
	PUBLIC EvaluateFamilyAlpha
	PUBLIC EvaluateEsd2SystemVE,ConfigureEsd2System,EvaluateEsd2PurePsat,EvaluateEsd2ExcessPoint
	PUBLIC RegisterEsd2PartnerPure
	Integer, SAVE :: nRegisteredPartners=0,registeredPartnerCas(nmx)=0
	DoublePrecision, SAVE :: registeredPartnerQ(nmx)=0.d0,registeredPartnerEps(nmx)=0.d0,registeredPartnerB(nmx)=0.d0
CONTAINS
	SUBROUTINE RegisterEsd2PartnerPure(nPartners,partnerCas,qValues,epsValues,bValues)
		Integer, INTENT(IN) :: nPartners,partnerCas(nPartners)
		DoublePrecision, INTENT(IN) :: qValues(nPartners),epsValues(nPartners),bValues(nPartners)
		nRegisteredPartners=nPartners;registeredPartnerCas=0;registeredPartnerQ=0.d0
		registeredPartnerEps=0.d0;registeredPartnerB=0.d0
		registeredPartnerCas(1:nPartners)=partnerCas
		registeredPartnerQ(1:nPartners)=qValues;registeredPartnerEps(1:nPartners)=epsValues
		registeredPartnerB(1:nPartners)=bValues
	END SUBROUTINE RegisterEsd2PartnerPure

	SUBROUTINE ApplyRegisteredPartner(partnerCas)
		USE EsdParms
		USE GlobConst, only:bVolCc_mol
		Integer, INTENT(IN) :: partnerCas
		Integer i
		do i=1,nRegisteredPartners
			if(registeredPartnerCas(i)/=partnerCas)cycle
			q(2)=registeredPartnerQ(i);c(2)=1.d0+(q(2)-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
			eokP(2)=registeredPartnerEps(i);vx(2)=registeredPartnerB(i);bVolCc_mol(2)=registeredPartnerB(i)
			return
		enddo
	END SUBROUTINE ApplyRegisteredPartner
	SUBROUTINE SetEsd2PolarMode(iPolar,alphaD2,enableMEM2,iErr)
		USE GlobConst, only:iEosOpt
		USE EsdParms
		Integer, INTENT(IN) :: iPolar
		DoublePrecision, INTENT(IN) :: alphaD2
		Logical, INTENT(IN) :: enableMEM2
		Integer, INTENT(OUT) :: iErr
		iErr=0
		if(iEosOpt/=23 .or. iPolar<1 .or. iPolar>nmx .or. alphaD2<0.d0)then
			iErr=11;return
		endif
		alphaPolarD2=0.d0
		alphaPolarD2(iPolar)=alphaD2
		useEsd2Polar=alphaD2>0.d0
		usePolarMEM2=enableMEM2
		disableEsd2MEM2=.not.enableMEM2
		! The polar formulation has one transferable global pair, never fitted.
		ESD2_K10=-0.7d0
		ESD2_K11=2.d0
	END SUBROUTINE SetEsd2PolarMode

	SUBROUTINE FitEsd2PureAtAlpha(nPts,tData,pData,weights,alphaD2,qLower,qUpper, &
	& qFit,epsFit,bFit,logRmse,iErr)
		! Fit q to multipoint Psat.  At every q, epsilon and b are solved so the
		! active (polar or nonpolar) ESD2 reproduces the experimental Tc and Pc.
		USE EsdParms
		USE GlobConst, only:iEosOpt,bVolCc_mol
		Integer, INTENT(IN) :: nPts
		DoublePrecision, INTENT(IN) :: tData(nPts),pData(nPts),weights(nPts)
		DoublePrecision, INTENT(IN) :: alphaD2,qLower,qUpper
		DoublePrecision, INTENT(OUT) :: qFit,epsFit,bFit,logRmse
		Integer, INTENT(OUT) :: iErr
		Integer i,best,iter,iErrLocal
		DoublePrecision grid(7),scores(7),epsGrid(7),bGrid(7),lo,hi,x1,x2,f1,f2
		DoublePrecision eps1,b1,eps2,b2,phi
		iErr=0
		if(iEosOpt/=23 .or. nPts<3 .or. qLower<=0.d0 .or. qUpper<=qLower .or. &
		& MINVAL(tData)<=0.d0 .or. MINVAL(pData)<=0.d0)then
			iErr=11;return
		endif
		CALL SetEsd2PolarMode(1,alphaD2,.FALSE.,iErrLocal)
		best=1
		do i=1,7
			grid(i)=qLower+(qUpper-qLower)*DBLE(i-1)/6.d0
			CALL PureScoreAtQ(grid(i),nPts,tData,pData,weights,scores(i),epsGrid(i),bGrid(i),iErrLocal)
			if(iErrLocal>0)scores(i)=1.d6
			if(scores(i) < scores(best))best=i
		enddo
		!best=MINLOC(scores,1)
		lo=grid(MAX(1,best-1));hi=grid(MIN(7,best+1))
		phi=(DSQRT(5.d0)-1.d0)/2.d0
		x1=hi-phi*(hi-lo);x2=lo+phi*(hi-lo)
		CALL PureScoreAtQ(x1,nPts,tData,pData,weights,f1,eps1,b1,iErrLocal)
		CALL PureScoreAtQ(x2,nPts,tData,pData,weights,f2,eps2,b2,iErrLocal)
		do iter=1,28
			if(ABS(hi-lo)<2.d-3)exit
			if(f1<=f2)then
				hi=x2;x2=x1;f2=f1;eps2=eps1;b2=b1
				x1=hi-phi*(hi-lo)
				CALL PureScoreAtQ(x1,nPts,tData,pData,weights,f1,eps1,b1,iErrLocal)
			else
				lo=x1;x1=x2;f1=f2;eps1=eps2;b1=b2
				x2=lo+phi*(hi-lo)
				CALL PureScoreAtQ(x2,nPts,tData,pData,weights,f2,eps2,b2,iErrLocal)
			endif
		enddo
		if(f1<=f2)then
			qFit=x1;logRmse=f1;epsFit=eps1;bFit=b1
		else
			qFit=x2;logRmse=f2;epsFit=eps2;bFit=b2
		endif
		if(scores(best)<logRmse)then
			qFit=grid(best);logRmse=scores(best);epsFit=epsGrid(best);bFit=bGrid(best)
		endif
		if(logRmse>=1.d5)then;iErr=12;return;endif
		q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
	END SUBROUTINE FitEsd2PureAtAlpha

	SUBROUTINE PureScoreAtQ(qTrial,nPts,tData,pData,weights,score,epsOut,bOut,iErr)
		USE EsdParms
		USE GlobConst, only:bVolCc_mol
		Integer, INTENT(IN) :: nPts
		DoublePrecision, INTENT(IN) :: qTrial,tData(nPts),pData(nPts),weights(nPts)
		DoublePrecision, INTENT(OUT) :: score,epsOut,bOut
		Integer, INTENT(OUT) :: iErr
		Integer i,iPsat
		DoublePrecision sumW,pCalc,chemPot(nmx),rhoL,rhoV,uL,uV,res
		epsOut=MAX(50.d0,eokP(1));bOut=MAX(5.d0,vx(1))
		CALL SolveCriticalAtQ(qTrial,epsOut,bOut,iErr)
		if(iErr>0)then;score=1.d6;return;endif
		q(1)=qTrial;c(1)=1.d0+(qTrial-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsOut;vx(1)=bOut;bVolCc_mol(1)=bOut
		score=0.d0;sumW=0.d0
		do i=1,nPts
			pCalc=pData(i)
			CALL PsatEar(tData(i),pCalc,chemPot,rhoL,rhoV,uL,uV,iPsat)
			if(iPsat>10 .or. pCalc<=0.d0)then
				res=10.d0
			else
				res=DLOG(pCalc/pData(i))
			endif
			score=score+MAX(weights(i),0.d0)*res*res
			sumW=sumW+MAX(weights(i),0.d0)
		enddo
		if(sumW<=zeroTol)sumW=DBLE(nPts)
		score=DSQRT(score/sumW)
	END SUBROUTINE PureScoreAtQ

	SUBROUTINE SolveCriticalAtQ(qTrial,epsValue,bValue,iErr,etaOut)
		USE GlobConst, only:zeroTol
		Integer, INTENT(OUT) :: iErr
		DoublePrecision, INTENT(IN) :: qTrial
		DoublePrecision, INTENT(INOUT) :: epsValue,bValue
		DoublePrecision, INTENT(OUT), OPTIONAL :: etaOut
		Integer iter,iEta,iScale,j,iLocal
		DoublePrecision etaStarts(5),epsScales(3),v(3),trial(3),res(3),resStep(3)
		DoublePrecision jac(3,3),rhs(3),delta(3),steps(3),norm,bestNorm,best(3),scale
		Logical solved
		etaStarts=(/0.08d0,0.13d0,0.20d0,0.30d0,0.40d0/)
		epsScales=(/0.70d0,1.d0,1.35d0/)
		bestNorm=1.d99;best=(/epsValue,bValue,0.15d0/)
		do iScale=1,3
			do iEta=1,5
				v=(/MAX(50.d0,epsValue*epsScales(iScale)),MAX(5.d0,bValue),etaStarts(iEta)/)
				do iter=1,35
					CALL CriticalResidual(qTrial,v(1),v(2),v(3),res,iLocal)
					if(iLocal/=0)exit
					norm=MAXVAL(ABS(res))
					if(norm<bestNorm)then;bestNorm=norm;best=v;endif
					if(norm<2.d-7)exit
					steps=(/MAX(1.d-3,2.d-4*v(1)),MAX(1.d-4,2.d-4*v(2)),MAX(1.d-6,2.d-4*v(3))/)
					do j=1,3
						trial=v;trial(j)=trial(j)+steps(j)
						CALL CriticalResidual(qTrial,trial(1),trial(2),trial(3),resStep,iLocal)
						if(iLocal/=0)exit
						jac(:,j)=(resStep-res)/steps(j)
					enddo
					if(iLocal/=0)exit
					rhs=-res
					CALL SolveLinear3(jac,rhs,delta,solved)
					if(.not.solved)exit
					scale=MIN(1.d0,0.25d0*v(1)/MAX(ABS(delta(1)),zeroTol), &
					&                 0.25d0*v(2)/MAX(ABS(delta(2)),zeroTol), &
					&                 0.25d0*v(3)/MAX(ABS(delta(3)),zeroTol))
					v=v+scale*delta
					v(1)=MAX(50.d0,MIN(2000.d0,v(1)))
					v(2)=MAX(5.d0,MIN(250.d0,v(2)))
					v(3)=MAX(0.025d0,MIN(0.48d0,v(3)))
				enddo
			enddo
		enddo
		epsValue=best(1);bValue=best(2)
		if(PRESENT(etaOut))etaOut=best(3)
		if(bestNorm<2.d-5)then;iErr=0;else;iErr=17;endif
	END SUBROUTINE SolveCriticalAtQ

	SUBROUTINE GetEsd2CriticalState(qValue,epsValue,bValue,etaC,vC,zC,iErr)
		USE GlobConst, only:Tc,Pc,Rgas
		DoublePrecision, INTENT(IN) :: qValue,epsValue,bValue
		DoublePrecision, INTENT(OUT) :: etaC,vC,zC
		Integer, INTENT(OUT) :: iErr
		DoublePrecision epsWork,bWork
		epsWork=epsValue;bWork=bValue
		CALL SolveCriticalAtQ(qValue,epsWork,bWork,iErr,etaC)
		if(iErr/=0 .or. etaC<=0.d0)then
			vC=-9.999d99;zC=-9.999d99;return
		endif
		vC=bWork/etaC
		zC=Pc(1)*vC/(Rgas*Tc(1))
	END SUBROUTINE GetEsd2CriticalState

	SUBROUTINE PredictEsd2CriticalFixedParameters(qValue,epsValue,bValue,alphaD2, &
	& tC,pC,etaC,vC,zC,iErr,tInitial,etaInitial)
		! Predict the unconstrained critical point after changing alpha_p while
		! holding the nonpolar q, epsilon, and b parameters exactly fixed.
		USE GlobConst, only:TcData=>Tc,PcData=>Pc,Rgas
		DoublePrecision, INTENT(IN) :: qValue,epsValue,bValue,alphaD2
		DoublePrecision, INTENT(OUT) :: tC,pC,etaC,vC,zC
		Integer, INTENT(OUT) :: iErr
		DoublePrecision, INTENT(IN), OPTIONAL :: tInitial,etaInitial
		Integer iT,iEta,iter,j,iLocal,modeErr
		DoublePrecision tStarts(6),etaStarts(6),v(2),trial(2),res(2),resStep(2)
		DoublePrecision jac(2,2),rhs(2),delta(2),steps(2),norm,trialNorm
		DoublePrecision bestNorm,bestDistance,distance,best(2),scale,pTrial
		Logical solved
		CALL SetEsd2PolarMode(1,alphaD2,.FALSE.,modeErr)
		if(modeErr/=0 .or. qValue<=0.d0 .or. epsValue<=0.d0 .or. bValue<=0.d0)then
			iErr=11;tC=-9.999d99;pC=-9.999d99;etaC=-9.999d99
			vC=-9.999d99;zC=-9.999d99;return
		endif
		CALL SetPureCandidate(qValue,epsValue,bValue)
		tStarts=(/0.65d0*TcData(1),0.82d0*TcData(1),TcData(1),1.18d0*TcData(1), &
		& 1.38d0*TcData(1),1.65d0*TcData(1)/)
		etaStarts=(/0.055d0,0.085d0,0.125d0,0.18d0,0.27d0,0.39d0/)
		if(PRESENT(tInitial))tStarts(1)=tInitial
		if(PRESENT(etaInitial))etaStarts(1)=etaInitial
		bestNorm=1.d99;bestDistance=1.d99;best=(/TcData(1),0.125d0/)
		do iT=1,SIZE(tStarts)
			do iEta=1,SIZE(etaStarts)
				v=(/tStarts(iT),etaStarts(iEta)/)
				do iter=1,55
					CALL FixedCriticalResidual(qValue,epsValue,bValue,v(1),v(2),res,iLocal)
					if(iLocal/=0)exit
					norm=MAXVAL(ABS(res))
					distance=ABS(v(1)-tStarts(1))/MAX(tStarts(1),1.d0)+ABS(v(2)-etaStarts(1))
					if(norm<bestNorm-1.d-10 .or. (ABS(norm-bestNorm)<1.d-10 .and. distance<bestDistance))then
						bestNorm=norm;bestDistance=distance;best=v
					endif
					if(norm<2.d-9)exit
					steps=(/MAX(2.d-2,2.d-4*v(1)),MAX(1.d-6,2.d-4*v(2))/)
					do j=1,2
						trial=v;trial(j)=trial(j)+steps(j)
						CALL FixedCriticalResidual(qValue,epsValue,bValue,trial(1),trial(2),resStep,iLocal)
						if(iLocal/=0)exit
						jac(:,j)=(resStep-res)/steps(j)
					enddo
					if(iLocal/=0)exit
					rhs=-res
					CALL SolveLinear2(jac,rhs,delta,solved)
					if(.not.solved)exit
					scale=MIN(1.d0,0.20d0*v(1)/MAX(ABS(delta(1)),zeroTol), &
					&                 0.20d0*v(2)/MAX(ABS(delta(2)),zeroTol))
					! Backtrack until the derivative residual decreases.
					do j=1,12
						trial=v+scale*delta
						trial(1)=MAX(0.25d0*TcData(1),MIN(2.d0*TcData(1),trial(1)))
						trial(2)=MAX(0.025d0,MIN(0.48d0,trial(2)))
						CALL FixedCriticalResidual(qValue,epsValue,bValue,trial(1),trial(2),resStep,iLocal)
						if(iLocal==0)then
							trialNorm=MAXVAL(ABS(resStep))
							if(trialNorm<norm)exit
						endif
						scale=0.5d0*scale
					enddo
					if(iLocal/=0 .or. trialNorm>=norm)exit
					v=trial
				enddo
			enddo
		enddo
		tC=best(1);etaC=best(2)
		CALL SetPureCandidate(qValue,epsValue,bValue)
		CALL PurePressureEta(tC,etaC,bValue,pTrial,iLocal)
		pC=pTrial;vC=bValue/etaC
		if(iLocal==0 .and. pC>0.d0 .and. bestNorm<2.d-6)then
			zC=pC*vC/(Rgas*tC);iErr=0
		else
			zC=-9.999d99;iErr=18
		endif
	END SUBROUTINE PredictEsd2CriticalFixedParameters

	SUBROUTINE FixedCriticalResidual(qValue,epsValue,bValue,tKelvin,eta,residual,iErr)
		USE GlobConst, only:Pc
		DoublePrecision, INTENT(IN) :: qValue,epsValue,bValue,tKelvin,eta
		DoublePrecision, INTENT(OUT) :: residual(2)
		Integer, INTENT(OUT) :: iErr
		Integer j,iLocal
		DoublePrecision h,p(5),first,second,pScale
		h=MIN(2.d-4,0.10d0*eta,0.10d0*(0.5d0-eta))
		if(tKelvin<=0.d0 .or. h<=1.d-7 .or. eta-2.d0*h<=0.d0 .or. eta+2.d0*h>=0.5d0)then
			iErr=11;return
		endif
		CALL SetPureCandidate(qValue,epsValue,bValue)
		do j=-2,2
			CALL PurePressureEta(tKelvin,eta+DBLE(j)*h,bValue,p(j+3),iLocal)
			if(iLocal/=0)then;iErr=12;return;endif
		enddo
		first=(p(1)-8.d0*p(2)+8.d0*p(4)-p(5))/(12.d0*h)
		second=(-p(5)+16.d0*p(4)-30.d0*p(3)+16.d0*p(2)-p(1))/(12.d0*h*h)
		pScale=MAX(ABS(Pc(1)),1.d0)
		residual(1)=eta*first/pScale
		residual(2)=eta*eta*second/pScale
		iErr=0
	END SUBROUTINE FixedCriticalResidual

	SUBROUTINE CriticalResidual(qValue,epsValue,bValue,eta,residual,iErr)
		USE GlobConst, only:Tc,Pc
		DoublePrecision, INTENT(IN) :: qValue,epsValue,bValue,eta
		DoublePrecision, INTENT(OUT) :: residual(3)
		Integer, INTENT(OUT) :: iErr
		Integer j,iLocal
		DoublePrecision h,p(5),first,second
		h=2.d-4
		if(eta-2.d0*h<=0.d0 .or. eta+2.d0*h>=0.5d0)then;iErr=11;return;endif
		CALL SetPureCandidate(qValue,epsValue,bValue)
		do j=-2,2
			CALL PurePressureEta(Tc(1),eta+DBLE(j)*h,bValue,p(j+3),iLocal)
			if(iLocal/=0)then;iErr=12;return;endif
		enddo
		first=(p(1)-8.d0*p(2)+8.d0*p(4)-p(5))/(12.d0*h)
		second=(-p(5)+16.d0*p(4)-30.d0*p(3)+16.d0*p(2)-p(1))/(12.d0*h*h)
		residual(1)=p(3)/Pc(1)-1.d0
		residual(2)=eta*first/Pc(1)
		residual(3)=eta*eta*second/Pc(1)
		iErr=0
	END SUBROUTINE CriticalResidual

	SUBROUTINE PurePressureEta(tKelvin,eta,bValue,pMPa,iErr)
		USE GlobConst, only:nmx,Rgas
		DoublePrecision, INTENT(IN) :: tKelvin,eta,bValue
		DoublePrecision, INTENT(OUT) :: pMPa
		Integer, INTENT(OUT) :: iErr
		DoublePrecision n(nmx),fug(nmx),z,a,u,rho
		n=0.d0;n(1)=1.d0;rho=eta/bValue
		CALL FuVtot(1,tKelvin,1.d0/rho,n,1,fug,z,a,u,iErr)
		if(iErr>10 .or. .not.(z>-1.d99 .and. z<1.d99))then;iErr=13;return;endif
		pMPa=z*Rgas*tKelvin*rho
		iErr=0
	END SUBROUTINE PurePressureEta

	SUBROUTINE SolveLinear3(matrix,rhs,solution,solved)
		DoublePrecision, INTENT(IN) :: matrix(3,3),rhs(3)
		DoublePrecision, INTENT(OUT) :: solution(3)
		Logical, INTENT(OUT) :: solved
		DoublePrecision aug(3,4),row(4),factor
		Integer i,j,pivot
		aug(:,1:3)=matrix;aug(:,4)=rhs;solved=.FALSE.
		do i=1,3
			pivot=i
			do j=i+1,3
				if(ABS(aug(j,i))>ABS(aug(pivot,i)))pivot=j
			enddo
			if(ABS(aug(pivot,i))<1.d-14)return
			if(pivot/=i)then;row=aug(i,:);aug(i,:)=aug(pivot,:);aug(pivot,:)=row;endif
			aug(i,:)=aug(i,:)/aug(i,i)
			do j=1,3
				if(j==i)cycle
				factor=aug(j,i);aug(j,:)=aug(j,:)-factor*aug(i,:)
			enddo
		enddo
		solution=aug(:,4);solved=.TRUE.
	END SUBROUTINE SolveLinear3

	SUBROUTINE SolveLinear2(matrix,rhs,solution,solved)
		DoublePrecision, INTENT(IN) :: matrix(2,2),rhs(2)
		DoublePrecision, INTENT(OUT) :: solution(2)
		Logical, INTENT(OUT) :: solved
		DoublePrecision determinant
		determinant=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
		if(ABS(determinant)<1.d-16)then
			solution=0.d0;solved=.FALSE.;return
		endif
		solution(1)=(rhs(1)*matrix(2,2)-matrix(1,2)*rhs(2))/determinant
		solution(2)=(matrix(1,1)*rhs(2)-rhs(1)*matrix(2,1))/determinant
		solved=.TRUE.
	END SUBROUTINE SolveLinear2

	SUBROUTINE SetPureCandidate(qValue,epsValue,bValue)
		USE EsdParms
		USE GlobConst, only:bVolCc_mol
		DoublePrecision, INTENT(IN) :: qValue,epsValue,bValue
		q(1)=qValue;c(1)=1.d0+(qValue-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsValue;vx(1)=bValue;bVolCc_mol(1)=bValue
	END SUBROUTINE SetPureCandidate

	SUBROUTINE CalcEsd2ExcessVolume(tKelvin,pMPa,x1,veCalc,iErr)
		Integer, INTENT(OUT) :: iErr
		DoublePrecision, INTENT(IN) :: tKelvin,pMPa,x1
		DoublePrecision, INTENT(OUT) :: veCalc
		DoublePrecision x(nmx),fug(nmx),rhoMix,rho1,rho2,z,a,u
		Integer iLocal
		iErr=0;veCalc=0.d0
		if(x1<0.d0 .or. x1>1.d0 .or. tKelvin<=0.d0 .or. pMPa<=0.d0)then;iErr=11;return;endif
		x=0.d0;x(1)=x1;x(2)=1.d0-x1
		CALL FugiTP(tKelvin,pMPa,x,2,1,rhoMix,z,a,fug,u,iLocal)
		if(iLocal>10 .or. rhoMix<=0.d0)then;iErr=12;return;endif
		x(1)=1.d0;x(2)=0.d0
		CALL FugiTP(tKelvin,pMPa,x,2,1,rho1,z,a,fug,u,iLocal)
		if(iLocal>10 .or. rho1<=0.d0)then;iErr=13;return;endif
		x(1)=0.d0;x(2)=1.d0
		CALL FugiTP(tKelvin,pMPa,x,2,1,rho2,z,a,fug,u,iLocal)
		if(iLocal>10 .or. rho2<=0.d0)then;iErr=14;return;endif
		veCalc=1.d0/rhoMix-x1/rho1-(1.d0-x1)/rho2
	END SUBROUTINE CalcEsd2ExcessVolume

	SUBROUTINE FitEsd2BetaVE(nPts,tData,pData,x1Data,veData,betaLower,betaUpper,betaFit,rmse,iErr)
		USE BIPs, only:Lij
		Integer, INTENT(IN) :: nPts
		DoublePrecision, INTENT(IN) :: tData(nPts),pData(nPts),x1Data(nPts),veData(nPts)
		DoublePrecision, INTENT(IN) :: betaLower,betaUpper
		DoublePrecision, INTENT(OUT) :: betaFit,rmse
		Integer, INTENT(OUT) :: iErr
		Integer i,best,iter
		DoublePrecision grid(7),score(7),lo,hi,xA,xB,fA,fB,phi
		iErr=0
		if(nPts<3 .or. betaUpper<=betaLower)then;iErr=11;return;endif
		if(MAXVAL(tData)-MINVAL(tData)>0.25d0)then;iErr=12;return;endif ! one VE isotherm only
		best=1
		do i=1,7
			grid(i)=betaLower+(betaUpper-betaLower)*DBLE(i-1)/6.d0
			CALL BetaObjective(grid(i),nPts,tData,pData,x1Data,veData,score(i))
			if(score(i) < score(best))best=i
		enddo
		!best=MINLOC(score);
		lo=grid(MAX(1,best-1));hi=grid(MIN(7,best+1));phi=(DSQRT(5.d0)-1.d0)/2.d0
		xA=hi-phi*(hi-lo);xB=lo+phi*(hi-lo)
		CALL BetaObjective(xA,nPts,tData,pData,x1Data,veData,fA)
		CALL BetaObjective(xB,nPts,tData,pData,x1Data,veData,fB)
		do iter=1,28
			if(ABS(hi-lo)<2.d-5)exit
			if(fA<=fB)then;hi=xB;xB=xA;fB=fA;xA=hi-phi*(hi-lo);CALL BetaObjective(xA,nPts,tData,pData,x1Data,veData,fA)
			else;lo=xA;xA=xB;fA=fB;xB=lo+phi*(hi-lo);CALL BetaObjective(xB,nPts,tData,pData,x1Data,veData,fB);endif
		enddo
		if(fA<=fB)then;betaFit=xA;rmse=DSQRT(fA);else;betaFit=xB;rmse=DSQRT(fB);endif
		if(score(best)<rmse*rmse)then;betaFit=grid(best);rmse=DSQRT(score(best));endif
		Lij(1,2)=betaFit;Lij(2,1)=betaFit
		if(rmse>1.d5)iErr=13
	END SUBROUTINE FitEsd2BetaVE

	SUBROUTINE BetaObjective(beta,nPts,tData,pData,x1Data,veData,mse)
		USE BIPs, only:Lij
		Integer, INTENT(IN) :: nPts
		DoublePrecision, INTENT(IN) :: beta,tData(nPts),pData(nPts),x1Data(nPts),veData(nPts)
		DoublePrecision, INTENT(OUT) :: mse
		Integer i,iErr
		DoublePrecision ve
		Lij(1,2)=beta;Lij(2,1)=beta;mse=0.d0
		do i=1,nPts
			CALL CalcEsd2ExcessVolume(tData(i),pData(i),x1Data(i),ve,iErr)
			if(iErr>0)then;mse=mse+1.d10;else;mse=mse+(ve-veData(i))**2;endif
		enddo
		mse=mse/DBLE(nPts)
	END SUBROUTINE BetaObjective

	SUBROUTINE ScoreEsd2VLE(nPts,mode,tData,pData,x1Data,y1Data,score,pCalc,tCalc,yCalc,iErr,pointStatus)
		! mode=1: isothermal bubble P+y; mode=2: isobaric bubble T+y.
		! mode=3: isobaric bubble T only; mode=4: isothermal bubble P only.
		Integer, INTENT(IN) :: nPts,mode
		DoublePrecision, INTENT(IN) :: tData(nPts),pData(nPts),x1Data(nPts),y1Data(nPts)
		DoublePrecision, INTENT(OUT) :: score,pCalc(nPts),tCalc(nPts),yCalc(nPts)
		Integer, INTENT(OUT) :: iErr
		Integer, INTENT(OUT), OPTIONAL :: pointStatus(nPts)
		Integer i,ii,init,itMax,ier(20)
		DoublePrecision x(nmx),y(nmx),rPrimary,rY,pPrev,yPrev(nmx)
		Logical havePrevious
		iErr=0;score=0.d0
		if(PRESENT(pointStatus))pointStatus=0
		havePrevious=.FALSE.;pPrev=0.d0;yPrev=0.d0
		do ii=1,nPts
			if(mode==1 .or. mode==4)then;i=nPts-ii+1;else;i=ii;endif
			x=0.d0;x(1)=x1Data(i);x(2)=1.d0-x(1);y=0.d0
			pCalc(i)=pData(i);tCalc(i)=tData(i);init=1;itMax=150;ier=0
			if(mode==1 .or. mode==4)then
				if(havePrevious)then
					pCalc(i)=pPrev;y=yPrev;init=2
				elseif(y1Data(i)>=0.d0 .and. y1Data(i)<=1.d0)then
					y(1)=MAX(1.d-10,MIN(1.d0-1.d-10,y1Data(i)));y(2)=1.d0-y(1);init=2
				endif
				CALL BUBPL(tCalc(i),x,2,init,pCalc(i),itMax,y,ier)
				if(ier(1)/=0)then
					! Retry from the supplied pressure without a vapor guess, then PSTART.
					init=1;itMax=200;ier=0;y=0.d0;pCalc(i)=pData(i)
					CALL BUBPL(tCalc(i),x,2,init,pCalc(i),itMax,y,ier)
				endif
				if(ier(1)/=0)then
					init=0;itMax=200;ier=0;y=0.d0
					CALL BUBPL(tCalc(i),x,2,init,pCalc(i),itMax,y,ier)
				endif
				if(ier(1)/=0)CALL BootBPx(tCalc(i),x,2,init,pCalc(i),itMax,y,ier,0)
				if(pCalc(i)>0.d0)rPrimary=DLOG(pCalc(i)/pData(i))/0.01d0
			else
				CALL BUBTL(tCalc(i),x,2,init,pCalc(i),itMax,y,ier)
				rPrimary=(tCalc(i)-tData(i))/1.d0
			endif
			if(ier(1)/=0 .or. ((mode==1.or.mode==4).and.pCalc(i)<=0.d0))then
				rPrimary=25.d0
				if(mode==3 .or. mode==4)then;rY=0.d0;else;rY=25.d0;endif
				iErr=iErr+1;yCalc(i)=0.d0
				if(PRESENT(pointStatus))pointStatus(i)=MAX(1,ABS(ier(1)))
			else
				yCalc(i)=y(1)
				if(mode==1 .or. mode==4)then;pPrev=pCalc(i);yPrev=y;havePrevious=.TRUE.;endif
				if(mode==3 .or. mode==4)then;rY=0.d0;else;rY=(yCalc(i)-y1Data(i))/0.02d0;endif
			endif
			score=score+rPrimary*rPrimary+rY*rY
		enddo
		if(mode==3 .or. mode==4)then;score=score/DBLE(nPts);else;score=score/DBLE(2*nPts);endif
	END SUBROUTINE ScoreEsd2VLE

	SUBROUTINE FitEsd2KijVLE(nPts,mode,tData,pData,x1Data,y1Data,kLower,kUpper,kFit,score,iErr)
		USE BIPs, only:Kij
		Integer, INTENT(IN) :: nPts,mode
		DoublePrecision, INTENT(IN) :: tData(nPts),pData(nPts),x1Data(nPts),y1Data(nPts),kLower,kUpper
		DoublePrecision, INTENT(OUT) :: kFit,score
		Integer, INTENT(OUT) :: iErr
		Integer iter,iLocal
		DoublePrecision lo,hi,xA,xB,fA,fB,phi,pC(nPts),tC(nPts),yC(nPts)
		lo=kLower;hi=kUpper;phi=(DSQRT(5.d0)-1.d0)/2.d0
		xA=hi-phi*(hi-lo);xB=lo+phi*(hi-lo)
		Kij(1,2)=xA;Kij(2,1)=xA;CALL ScoreEsd2VLE(nPts,mode,tData,pData,x1Data,y1Data,fA,pC,tC,yC,iLocal)
		Kij(1,2)=xB;Kij(2,1)=xB;CALL ScoreEsd2VLE(nPts,mode,tData,pData,x1Data,y1Data,fB,pC,tC,yC,iLocal)
		do iter=1,28
			if(ABS(hi-lo)<2.d-5)exit
			if(fA<=fB)then
				hi=xB;xB=xA;fB=fA;xA=hi-phi*(hi-lo);Kij(1,2)=xA;Kij(2,1)=xA
				CALL ScoreEsd2VLE(nPts,mode,tData,pData,x1Data,y1Data,fA,pC,tC,yC,iLocal)
			else
				lo=xA;xA=xB;fA=fB;xB=lo+phi*(hi-lo);Kij(1,2)=xB;Kij(2,1)=xB
				CALL ScoreEsd2VLE(nPts,mode,tData,pData,x1Data,y1Data,fB,pC,tC,yC,iLocal)
			endif
		enddo
		if(fA<=fB)then;kFit=xA;score=fA;else;kFit=xB;score=fB;endif
		Kij(1,2)=kFit;Kij(2,1)=kFit;iErr=0
	END SUBROUTINE FitEsd2KijVLE

	SUBROUTINE FitEsd2BetaKijNested(polarCas,partnerCas,alpha,qFit,epsFit,bFit,nVe,tVe,pVe,xVe,veData, &
	& nVle,vleMode,tVle,pVle,xVle,yVle,kLower,kUpper,betaFit,kFit,veRmse,vleScore,iErr)
		! Preserve the intended parameter separation while treating the coupling correctly:
		! beta*(kij) is fitted only to VE, then kij is selected only from the VLE score.
		! Every VLE evaluation is preceded by a complete binary reconfiguration so that
		! failed flash iterations cannot leak state into the next objective evaluation.
		Integer, PARAMETER :: nGrid=31
		Integer, INTENT(IN) :: polarCas,partnerCas,nVe,nVle,vleMode
		DoublePrecision, INTENT(IN) :: alpha,qFit,epsFit,bFit,tVe(nVe),pVe(nVe),xVe(nVe),veData(nVe)
		DoublePrecision, INTENT(IN) :: tVle(nVle),pVle(nVle),xVle(nVle),yVle(nVle),kLower,kUpper
		DoublePrecision, INTENT(OUT) :: betaFit,kFit,veRmse,vleScore
		Integer, INTENT(OUT) :: iErr
		Integer i,best,iter,iA,iB,iCheck
		DoublePrecision kGrid(nGrid),fGrid(nGrid),betaGrid(nGrid),veGrid(nGrid)
		DoublePrecision lo,hi,xA,xB,fA,fB,betaA,betaB,veA,veB,phi,scoreCheck,betaCheck,veCheck
		iErr=0
		if(nVe<3 .or. nVle<3 .or. kUpper<=kLower)then;iErr=11;return;endif
		best=1
		do i=1,nGrid
			kGrid(i)=kLower+(kUpper-kLower)*DBLE(i-1)/DBLE(nGrid-1)
			CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,kGrid(i),nVe,tVe,pVe,xVe,veData, &
			& nVle,vleMode,tVle,pVle,xVle,yVle,betaGrid(i),veGrid(i),fGrid(i),iA)
			if(iA/=0)fGrid(i)=1.d6
			if(fgrid(i) < fgrid(best))best=i
		enddo
		!best=MINLOC(fGrid)
		if(fGrid(best)>=1.d6)then;iErr=12;return;endif
		betaFit=betaGrid(best);kFit=kGrid(best);veRmse=veGrid(best);vleScore=fGrid(best)
		lo=kGrid(MAX(1,best-1));hi=kGrid(MIN(nGrid,best+1))
		if(hi>lo)then
			phi=(DSQRT(5.d0)-1.d0)/2.d0;xA=hi-phi*(hi-lo);xB=lo+phi*(hi-lo)
			CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,xA,nVe,tVe,pVe,xVe,veData, &
			& nVle,vleMode,tVle,pVle,xVle,yVle,betaA,veA,fA,iA)
			CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,xB,nVe,tVe,pVe,xVe,veData, &
			& nVle,vleMode,tVle,pVle,xVle,yVle,betaB,veB,fB,iB)
			do iter=1,24
				if(ABS(hi-lo)<2.d-5)exit
				if(fA<=fB)then
					hi=xB;xB=xA;fB=fA;betaB=betaA;veB=veA;xA=hi-phi*(hi-lo)
					CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,xA,nVe,tVe,pVe,xVe,veData, &
					& nVle,vleMode,tVle,pVle,xVle,yVle,betaA,veA,fA,iA)
				else
					lo=xA;xA=xB;fA=fB;betaA=betaB;veA=veB;xB=lo+phi*(hi-lo)
					CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,xB,nVe,tVe,pVe,xVe,veData, &
					& nVle,vleMode,tVle,pVle,xVle,yVle,betaB,veB,fB,iB)
				endif
			enddo
			if(iA==0 .and. fA<vleScore)then;kFit=xA;betaFit=betaA;veRmse=veA;vleScore=fA;endif
			if(iB==0 .and. fB<vleScore)then;kFit=xB;betaFit=betaB;veRmse=veB;vleScore=fB;endif
		endif
		! A repeated clean evaluation is a mandatory acceptance test.  This catches
		! the stale 625-point failure plateau that affected the original coordinate fit.
		CALL NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,kFit,nVe,tVe,pVe,xVe,veData, &
		& nVle,vleMode,tVle,pVle,xVle,yVle,betaCheck,veCheck,scoreCheck,iCheck)
		if(iCheck/=0 .or. ABS(scoreCheck-vleScore)>1.d-7*(1.d0+ABS(vleScore)))then
			iErr=13;return
		endif
		betaFit=betaCheck;veRmse=veCheck;vleScore=scoreCheck
	END SUBROUTINE FitEsd2BetaKijNested

	SUBROUTINE NestedKijObjective(polarCas,partnerCas,alpha,qFit,epsFit,bFit,kValue,nVe,tVe,pVe,xVe,veData, &
	& nVle,vleMode,tVle,pVle,xVle,yVle,betaFit,veRmse,vleScore,iErr)
		USE BIPs, only:Kij,Lij
		Integer, INTENT(IN) :: polarCas,partnerCas,nVe,nVle,vleMode
		DoublePrecision, INTENT(IN) :: alpha,qFit,epsFit,bFit,kValue,tVe(nVe),pVe(nVe),xVe(nVe),veData(nVe)
		DoublePrecision, INTENT(IN) :: tVle(nVle),pVle(nVle),xVle(nVle),yVle(nVle)
		DoublePrecision, INTENT(OUT) :: betaFit,veRmse,vleScore
		Integer, INTENT(OUT) :: iErr
		Integer iLocal
		DoublePrecision pCalc(nVle),tCalc(nVle),yCalc(nVle)
		CALL ConfigureEsd2System(polarCas,partnerCas,alpha,qFit,epsFit,bFit,0.d0,kValue,iLocal)
		if(iLocal/=0)then;iErr=21;vleScore=1.d6;return;endif
		CALL FitEsd2BetaVE(nVe,tVe,pVe,xVe,veData,-0.25d0,0.25d0,betaFit,veRmse,iLocal)
		if(iLocal/=0)then;iErr=22;vleScore=1.d6;return;endif
		! Validate VE and evaluate VLE from a fresh, explicitly configured state.
		CALL EvaluateEsd2SystemVE(polarCas,partnerCas,alpha,qFit,epsFit,bFit,betaFit,kValue,nVe, &
		& tVe,pVe,xVe,veData,veRmse,iLocal)
		if(iLocal/=0)then;iErr=23;vleScore=1.d6;return;endif
		CALL ConfigureEsd2System(polarCas,partnerCas,alpha,qFit,epsFit,bFit,betaFit,kValue,iLocal)
		if(iLocal/=0)then;iErr=24;vleScore=1.d6;return;endif
		CALL ScoreEsd2VLE(nVle,vleMode,tVle,pVle,xVle,yVle,vleScore,pCalc,tCalc,yCalc,iLocal)
		! Point failures remain encoded in the score; they do not make the profile disappear.
		iErr=0
	END SUBROUTINE NestedKijObjective

	SUBROUTINE FitEsd2AlphaBinaryNested(nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
	& nVle,vleMode,tVle,pVle,xVle,yVle,alphaLower,alphaUpper,qLower,qUpper, &
	& alphaFit,qFit,epsFit,bFit,betaFit,veRmse,vleScore,iErr)
		! Binary version of beta*(alpha) followed by alpha from VLE, with kij=0.
		USE EsdParms
		USE GlobConst, only:bVolCc_mol
		USE BIPs, only:Kij
		Integer, INTENT(IN) :: nVp,nVe,nVle,vleMode
		DoublePrecision, INTENT(IN) :: tVp(nVp),pVp(nVp),wVp(nVp)
		DoublePrecision, INTENT(IN) :: tVe(nVe),pVe(nVe),xVe(nVe),veData(nVe)
		DoublePrecision, INTENT(IN) :: tVle(nVle),pVle(nVle),xVle(nVle),yVle(nVle)
		DoublePrecision, INTENT(IN) :: alphaLower,alphaUpper,qLower,qUpper
		DoublePrecision, INTENT(OUT) :: alphaFit,qFit,epsFit,bFit,betaFit,veRmse,vleScore
		Integer, INTENT(OUT) :: iErr
		Integer iter,iLocal
		DoublePrecision lo,hi,xA,xB,fA,fB,phi,qA,eA,bA,betaA,veA,qB,eB,bB,betaB,veB
		DoublePrecision seedQ,seedE,seedB
		seedQ=q(1);seedE=eokP(1);seedB=vx(1);lo=alphaLower;hi=alphaUpper
		phi=(DSQRT(5.d0)-1.d0)/2.d0;xA=hi-phi*(hi-lo);xB=lo+phi*(hi-lo)
		CALL NestedAlphaObjective(xA,seedQ,seedE,seedB,nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
		& nVle,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper,fA,qA,eA,bA,betaA,veA,iLocal)
		CALL NestedAlphaObjective(xB,seedQ,seedE,seedB,nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
		& nVle,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper,fB,qB,eB,bB,betaB,veB,iLocal)
		do iter=1,24
			if(ABS(hi-lo)<2.d-3)exit
			if(fA<=fB)then
				hi=xB;xB=xA;fB=fA;qB=qA;eB=eA;bB=bA;betaB=betaA;veB=veA;xA=hi-phi*(hi-lo)
				CALL NestedAlphaObjective(xA,seedQ,seedE,seedB,nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
				& nVle,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper,fA,qA,eA,bA,betaA,veA,iLocal)
			else
				lo=xA;xA=xB;fA=fB;qA=qB;eA=eB;bA=bB;betaA=betaB;veA=veB;xB=lo+phi*(hi-lo)
				CALL NestedAlphaObjective(xB,seedQ,seedE,seedB,nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
				& nVle,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper,fB,qB,eB,bB,betaB,veB,iLocal)
			endif
		enddo
		if(fA<=fB)then
			alphaFit=xA;vleScore=fA;qFit=qA;epsFit=eA;bFit=bA;betaFit=betaA;veRmse=veA
		else
			alphaFit=xB;vleScore=fB;qFit=qB;epsFit=eB;bFit=bB;betaFit=betaB;veRmse=veB
		endif
		CALL SetEsd2PolarMode(1,alphaFit,.FALSE.,iLocal)
		q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
		Kij(1,2)=0.d0;Kij(2,1)=0.d0;iErr=0
	END SUBROUTINE FitEsd2AlphaBinaryNested

	SUBROUTINE NestedAlphaObjective(alpha,seedQ,seedE,seedB,nVp,tVp,pVp,wVp,nVe,tVe,pVe,xVe,veData, &
	& nVle,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper,score,qFit,epsFit,bFit,betaFit,veRmse,iErr)
		USE EsdParms
		USE GlobConst, only:bVolCc_mol
		USE BIPs, only:Kij
		Integer, INTENT(IN) :: nVp,nVe,nVle,vleMode
		DoublePrecision, INTENT(IN) :: alpha,seedQ,seedE,seedB,tVp(nVp),pVp(nVp),wVp(nVp)
		DoublePrecision, INTENT(IN) :: tVe(nVe),pVe(nVe),xVe(nVe),veData(nVe)
		DoublePrecision, INTENT(IN) :: tVle(nVle),pVle(nVle),xVle(nVle),yVle(nVle),qLower,qUpper
		DoublePrecision, INTENT(OUT) :: score,qFit,epsFit,bFit,betaFit,veRmse
		Integer, INTENT(OUT) :: iErr
		DoublePrecision pureScore,pC(nVle),tC(nVle),yC(nVle)
		q(1)=seedQ;eokP(1)=seedE;vx(1)=seedB;bVolCc_mol(1)=seedB
		CALL FitEsd2PureAtAlpha(nVp,tVp,pVp,wVp,alpha,qLower,qUpper,qFit,epsFit,bFit,pureScore,iErr)
		if(iErr>0)then;score=1.d6;veRmse=1.d6;return;endif
		Kij(1,2)=0.d0;Kij(2,1)=0.d0
		CALL FitEsd2BetaVE(nVe,tVe,pVe,xVe,veData,-0.25d0,0.25d0,betaFit,veRmse,iErr)
		if(iErr>0)then;score=1.d6;return;endif
		CALL ScoreEsd2VLE(nVle,vleMode,tVle,pVle,xVle,yVle,score,pC,tC,yC,iErr)
	END SUBROUTINE NestedAlphaObjective

	SUBROUTINE FitEsd2AlphaFamilyNested(nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp, &
	& nVeTotal,veStart,veCount,tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode, &
	& tVle,pVle,xVle,yVle,alphaLower,alphaUpper,qLower,qUpper,doLopo,alphaFit,qFit,epsFit,bFit, &
	& betaFit,partnerScore,balancedScore,lopoAlpha,lopoHeldScore,iErr)
		! All compositions must be ordered with the polar compound as component 1.
		Integer, INTENT(IN) :: nPartners,polarCas,partnerCas(nPartners),nVp,nVeTotal,nVleTotal
		Integer, INTENT(IN) :: veStart(nPartners),veCount(nPartners),vleStart(nPartners),vleCount(nPartners),vleMode(nPartners)
		DoublePrecision, INTENT(IN) :: tVp(nVp),pVp(nVp),wVp(nVp)
		DoublePrecision, INTENT(IN) :: tVe(nVeTotal),pVe(nVeTotal),xVe(nVeTotal),veData(nVeTotal)
		DoublePrecision, INTENT(IN) :: tVle(nVleTotal),pVle(nVleTotal),xVle(nVleTotal),yVle(nVleTotal)
		DoublePrecision, INTENT(IN) :: alphaLower,alphaUpper,qLower,qUpper
		Logical, INTENT(IN) :: doLopo
		DoublePrecision, INTENT(OUT) :: alphaFit,qFit,epsFit,bFit,betaFit(nPartners),partnerScore(nPartners),balancedScore
		DoublePrecision, INTENT(OUT) :: lopoAlpha(nPartners),lopoHeldScore(nPartners)
		Integer, INTENT(OUT) :: iErr
		Integer held,iLocal
		CALL OptimizeFamilyAlpha(0,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,alphaLower,alphaUpper, &
		& qLower,qUpper,alphaFit,qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iErr)
		if(iErr/=0)return
		lopoAlpha=alphaFit;lopoHeldScore=partnerScore
		if(.not.doLopo .or. nPartners<2)return
		do held=1,nPartners
			CALL OptimizeFamilyAlpha(held,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
			& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,alphaLower,alphaUpper, &
			& qLower,qUpper,lopoAlpha(held),qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iLocal)
			if(iLocal/=0)then;iErr=20+held;return;endif
			CALL EvaluateFamilyAlpha(lopoAlpha(held),0,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal, &
			& veStart,veCount,tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle, &
			& qLower,qUpper,qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iLocal)
			lopoHeldScore(held)=partnerScore(held)
		enddo
		! Restore the full-family optimum after the validation loops.
		CALL OptimizeFamilyAlpha(0,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,alphaLower,alphaUpper, &
		& qLower,qUpper,alphaFit,qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iErr)
	END SUBROUTINE FitEsd2AlphaFamilyNested

	SUBROUTINE OptimizeFamilyAlpha(excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal, &
	& veStart,veCount,tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle, &
	& alphaLower,alphaUpper,qLower,qUpper,alphaFit,qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iErr)
		Integer, INTENT(IN) :: excludePartner,nPartners,polarCas,partnerCas(nPartners),nVp,nVeTotal,nVleTotal
		Integer, INTENT(IN) :: veStart(nPartners),veCount(nPartners),vleStart(nPartners),vleCount(nPartners),vleMode(nPartners)
		DoublePrecision, INTENT(IN) :: tVp(nVp),pVp(nVp),wVp(nVp),tVe(nVeTotal),pVe(nVeTotal),xVe(nVeTotal),veData(nVeTotal)
		DoublePrecision, INTENT(IN) :: tVle(nVleTotal),pVle(nVleTotal),xVle(nVleTotal),yVle(nVleTotal)
		DoublePrecision, INTENT(IN) :: alphaLower,alphaUpper,qLower,qUpper
		DoublePrecision, INTENT(OUT) :: alphaFit,qFit,epsFit,bFit,betaFit(nPartners),partnerScore(nPartners),balancedScore
		Integer, INTENT(OUT) :: iErr
		Integer iter,iA,iB
		DoublePrecision lo,hi,xA,xB,fA,fB,phi,qA,eA,bA,qB,eB,bB
		DoublePrecision betaA(nPartners),scoreA(nPartners),betaB(nPartners),scoreB(nPartners)
		DoublePrecision fEnd,qEnd,eEnd,bEnd,betaEnd(nPartners),scoreEnd(nPartners)
		lo=alphaLower;hi=alphaUpper;phi=(DSQRT(5.d0)-1.d0)/2.d0
		xA=hi-phi*(hi-lo);xB=lo+phi*(hi-lo)
		CALL EvaluateFamilyAlpha(xA,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
		& qA,eA,bA,betaA,scoreA,fA,iA)
		CALL EvaluateFamilyAlpha(xB,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
		& qB,eB,bB,betaB,scoreB,fB,iB)
		do iter=1,24
			if(ABS(hi-lo)<2.d-3)exit
			if(fA<=fB)then
				hi=xB;xB=xA;fB=fA;qB=qA;eB=eA;bB=bA;betaB=betaA;scoreB=scoreA;xA=hi-phi*(hi-lo)
				CALL EvaluateFamilyAlpha(xA,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
				& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
				& qA,eA,bA,betaA,scoreA,fA,iA)
			else
				lo=xA;xA=xB;fA=fB;qA=qB;eA=eB;bA=bB;betaA=betaB;scoreA=scoreB;xB=lo+phi*(hi-lo)
				CALL EvaluateFamilyAlpha(xB,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
				& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
				& qB,eB,bB,betaB,scoreB,fB,iB)
			endif
		enddo
		if(fA<=fB)then;alphaFit=xA;balancedScore=fA;qFit=qA;epsFit=eA;bFit=bA;betaFit=betaA;partnerScore=scoreA
		else;alphaFit=xB;balancedScore=fB;qFit=qB;epsFit=eB;bFit=bB;betaFit=betaB;partnerScore=scoreB;endif
		! Explicitly test both bounds because a bounded golden search does not.
		CALL EvaluateFamilyAlpha(alphaLower,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
		& qEnd,eEnd,bEnd,betaEnd,scoreEnd,fEnd,iA)
		if(fEnd<balancedScore)then;alphaFit=alphaLower;balancedScore=fEnd;qFit=qEnd;epsFit=eEnd;bFit=bEnd;betaFit=betaEnd;partnerScore=scoreEnd;endif
		CALL EvaluateFamilyAlpha(alphaUpper,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
		& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
		& qEnd,eEnd,bEnd,betaEnd,scoreEnd,fEnd,iB)
		if(fEnd<balancedScore)then;alphaFit=alphaUpper;balancedScore=fEnd;qFit=qEnd;epsFit=eEnd;bFit=bEnd;betaFit=betaEnd;partnerScore=scoreEnd;endif
		iErr=0;if(balancedScore>=1.d5)iErr=31
	END SUBROUTINE OptimizeFamilyAlpha

	SUBROUTINE EvaluateFamilyAlpha(alpha,excludePartner,nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal, &
	& veStart,veCount,tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
	& qFit,epsFit,bFit,betaFit,partnerScore,balancedScore,iErr)
		USE EsdParms
		USE GlobConst, only:nmx,bVolCc_mol
		USE BIPs, only:Kij
		Integer, INTENT(IN) :: excludePartner,nPartners,polarCas,partnerCas(nPartners),nVp,nVeTotal,nVleTotal
		Integer, INTENT(IN) :: veStart(nPartners),veCount(nPartners),vleStart(nPartners),vleCount(nPartners),vleMode(nPartners)
		DoublePrecision, INTENT(IN) :: alpha,tVp(nVp),pVp(nVp),wVp(nVp),tVe(nVeTotal),pVe(nVeTotal),xVe(nVeTotal),veData(nVeTotal)
		DoublePrecision, INTENT(IN) :: tVle(nVleTotal),pVle(nVleTotal),xVle(nVleTotal),yVle(nVleTotal),qLower,qUpper
		DoublePrecision, INTENT(OUT) :: qFit,epsFit,bFit,betaFit(nPartners),partnerScore(nPartners),balancedScore
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),j,iLocal,i0,i1,nUsed
		DoublePrecision pureScore,veRmse,pC(nVleTotal),tC(nVleTotal),yC(nVleTotal)
		cas=0;cas(1)=polarCas
		CALL PGLWrapperStartup(1,23,cas,iLocal)
		if(iLocal/=0)then;iErr=41;balancedScore=1.d6;return;endif
		CALL FitEsd2PureAtAlpha(nVp,tVp,pVp,wVp,alpha,qLower,qUpper,qFit,epsFit,bFit,pureScore,iLocal)
		if(iLocal/=0)then;iErr=42;balancedScore=1.d6;return;endif
		balancedScore=0.d0;nUsed=0;betaFit=0.d0;partnerScore=0.d0
		do j=1,nPartners
			cas=0;cas(1)=polarCas;cas(2)=partnerCas(j)
			CALL PGLWrapperStartup(2,23,cas,iLocal)
			if(iLocal/=0)then;iErr=42+j;balancedScore=1.d6;return;endif
			CALL ApplyRegisteredPartner(partnerCas(j))
			CALL SetEsd2PolarMode(1,alpha,.FALSE.,iLocal)
			q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
			eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
			Kij(1,2)=0.d0;Kij(2,1)=0.d0
			i0=veStart(j);i1=i0+veCount(j)-1
			CALL FitEsd2BetaVE(veCount(j),tVe(i0:i1),pVe(i0:i1),xVe(i0:i1),veData(i0:i1), &
			& -0.25d0,0.25d0,betaFit(j),veRmse,iLocal)
			if(iLocal/=0)then;iErr=60+j;balancedScore=1.d6;return;endif
			i0=vleStart(j);i1=i0+vleCount(j)-1
			CALL ScoreEsd2VLE(vleCount(j),vleMode(j),tVle(i0:i1),pVle(i0:i1),xVle(i0:i1),yVle(i0:i1), &
			& partnerScore(j),pC(1:vleCount(j)),tC(1:vleCount(j)),yC(1:vleCount(j)),iLocal)
			if(j/=excludePartner)then;balancedScore=balancedScore+partnerScore(j);nUsed=nUsed+1;endif
		enddo
		if(nUsed<1)then;iErr=70;balancedScore=1.d6;else;iErr=0;balancedScore=balancedScore/DBLE(nUsed);endif
	END SUBROUTINE EvaluateFamilyAlpha

	SUBROUTINE FitEsd2NonpolarFamily(nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
	& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,qLower,qUpper, &
	& qFit,epsFit,bFit,betaFit,kijFit,veRmse,partnerScore,iErr)
		USE EsdParms
		USE GlobConst, only:nmx,bVolCc_mol
		USE BIPs, only:Kij,Lij
		Integer, INTENT(IN) :: nPartners,polarCas,partnerCas(nPartners),nVp,nVeTotal,nVleTotal
		Integer, INTENT(IN) :: veStart(nPartners),veCount(nPartners),vleStart(nPartners),vleCount(nPartners),vleMode(nPartners)
		DoublePrecision, INTENT(IN) :: tVp(nVp),pVp(nVp),wVp(nVp),tVe(nVeTotal),pVe(nVeTotal),xVe(nVeTotal),veData(nVeTotal)
		DoublePrecision, INTENT(IN) :: tVle(nVleTotal),pVle(nVleTotal),xVle(nVleTotal),yVle(nVleTotal),qLower,qUpper
		DoublePrecision, INTENT(OUT) :: qFit,epsFit,bFit,betaFit(nPartners),kijFit(nPartners),veRmse(nPartners),partnerScore(nPartners)
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),j,i0,i1,iLocal
		DoublePrecision pureScore
		cas=0;cas(1)=polarCas;CALL PGLWrapperStartup(1,23,cas,iLocal)
		if(iLocal/=0)then;iErr=81;return;endif
		CALL FitEsd2PureAtAlpha(nVp,tVp,pVp,wVp,0.d0,qLower,qUpper,qFit,epsFit,bFit,pureScore,iLocal)
		if(iLocal/=0)then;iErr=82;return;endif
		do j=1,nPartners
			cas=0;cas(1)=polarCas;cas(2)=partnerCas(j);CALL PGLWrapperStartup(2,23,cas,iLocal)
			if(iLocal/=0)then;iErr=82+j;return;endif
			CALL ApplyRegisteredPartner(partnerCas(j))
			CALL SetEsd2PolarMode(1,0.d0,.FALSE.,iLocal)
			q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
			eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
			i0=veStart(j);i1=i0+veCount(j)-1
			Kij(1,2)=0.d0;Kij(2,1)=0.d0;Lij(1,2)=0.d0;Lij(2,1)=0.d0
			CALL FitEsd2BetaKijNested(polarCas,partnerCas(j),0.d0,qFit,epsFit,bFit,veCount(j), &
			& tVe(i0:i1),pVe(i0:i1),xVe(i0:i1),veData(i0:i1),vleCount(j),vleMode(j), &
			& tVle(vleStart(j):vleStart(j)+vleCount(j)-1),pVle(vleStart(j):vleStart(j)+vleCount(j)-1), &
			& xVle(vleStart(j):vleStart(j)+vleCount(j)-1),yVle(vleStart(j):vleStart(j)+vleCount(j)-1), &
			& -0.20d0,0.25d0,betaFit(j),kijFit(j),veRmse(j),partnerScore(j),iLocal)
			if(iLocal/=0)then;iErr=120+j;return;endif
		enddo
		iErr=0
	END SUBROUTINE FitEsd2NonpolarFamily

	SUBROUTINE FitEsd2HybridFamily(nPartners,polarCas,partnerCas,nVp,tVp,pVp,wVp,nVeTotal,veStart,veCount, &
	& tVe,pVe,xVe,veData,nVleTotal,vleStart,vleCount,vleMode,tVle,pVle,xVle,yVle,alphaHalf,qLower,qUpper, &
	& qFit,epsFit,bFit,betaFit,kijFit,veRmse,partnerScore,iErr)
		USE EsdParms
		USE GlobConst, only:nmx,bVolCc_mol
		USE BIPs, only:Kij,Lij
		Integer, INTENT(IN) :: nPartners,polarCas,partnerCas(nPartners),nVp,nVeTotal,nVleTotal
		Integer, INTENT(IN) :: veStart(nPartners),veCount(nPartners),vleStart(nPartners),vleCount(nPartners),vleMode(nPartners)
		DoublePrecision, INTENT(IN) :: tVp(nVp),pVp(nVp),wVp(nVp),tVe(nVeTotal),pVe(nVeTotal),xVe(nVeTotal),veData(nVeTotal)
		DoublePrecision, INTENT(IN) :: tVle(nVleTotal),pVle(nVleTotal),xVle(nVleTotal),yVle(nVleTotal),alphaHalf,qLower,qUpper
		DoublePrecision, INTENT(OUT) :: qFit,epsFit,bFit,betaFit(nPartners),kijFit(nPartners),veRmse(nPartners),partnerScore(nPartners)
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),j,i0,i1,iLocal
		DoublePrecision pureScore
		cas=0;cas(1)=polarCas;CALL PGLWrapperStartup(1,23,cas,iLocal)
		if(iLocal/=0)then;iErr=181;return;endif
		CALL FitEsd2PureAtAlpha(nVp,tVp,pVp,wVp,alphaHalf,qLower,qUpper,qFit,epsFit,bFit,pureScore,iLocal)
		if(iLocal/=0)then;iErr=182;return;endif
		do j=1,nPartners
			cas=0;cas(1)=polarCas;cas(2)=partnerCas(j);CALL PGLWrapperStartup(2,23,cas,iLocal)
			if(iLocal/=0)then;iErr=182+j;return;endif
			CALL ApplyRegisteredPartner(partnerCas(j))
			CALL SetEsd2PolarMode(1,alphaHalf,.FALSE.,iLocal)
			q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
			eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
			Kij(1,2)=0.d0;Kij(2,1)=0.d0;Lij(1,2)=0.d0;Lij(2,1)=0.d0;kijFit(j)=0.d0
			i0=veStart(j);i1=i0+veCount(j)-1
			CALL FitEsd2BetaKijNested(polarCas,partnerCas(j),alphaHalf,qFit,epsFit,bFit,veCount(j), &
			& tVe(i0:i1),pVe(i0:i1),xVe(i0:i1),veData(i0:i1),vleCount(j),vleMode(j), &
			& tVle(vleStart(j):vleStart(j)+vleCount(j)-1),pVle(vleStart(j):vleStart(j)+vleCount(j)-1), &
			& xVle(vleStart(j):vleStart(j)+vleCount(j)-1),yVle(vleStart(j):vleStart(j)+vleCount(j)-1), &
			& -0.20d0,0.25d0,betaFit(j),kijFit(j),veRmse(j),partnerScore(j),iLocal)
			if(iLocal/=0)then;iErr=220+j;return;endif
		enddo
		iErr=0
	END SUBROUTINE FitEsd2HybridFamily

	SUBROUTINE EvaluateEsd2SystemVE(polarCas,partnerCas,alpha,qFit,epsFit,bFit,beta,kijValue,nPts, &
	& tData,pData,xData,veData,rmse,iErr)
		USE EsdParms
		USE GlobConst, only:nmx,bVolCc_mol
		USE BIPs, only:KijMatrix=>Kij,Lij
		Integer, INTENT(IN) :: polarCas,partnerCas,nPts
		DoublePrecision, INTENT(IN) :: alpha,qFit,epsFit,bFit,beta,kijValue
		DoublePrecision, INTENT(IN) :: tData(nPts),pData(nPts),xData(nPts),veData(nPts)
		DoublePrecision, INTENT(OUT) :: rmse
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),i,iLocal
		DoublePrecision veCalc,sse
		cas=0;cas(1)=polarCas;cas(2)=partnerCas
		CALL PGLWrapperStartup(2,23,cas,iLocal)
		if(iLocal/=0)then;iErr=1;rmse=1.d6;return;endif
		CALL ApplyRegisteredPartner(partnerCas)
		CALL SetEsd2PolarMode(1,alpha,.FALSE.,iLocal)
		if(iLocal/=0)then;iErr=2;rmse=1.d6;return;endif
		q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
		KijMatrix(1,2)=kijValue;KijMatrix(2,1)=kijValue;Lij(1,2)=beta;Lij(2,1)=beta
		sse=0.d0
		do i=1,nPts
			CALL CalcEsd2ExcessVolume(tData(i),pData(i),xData(i),veCalc,iLocal)
			if(iLocal/=0)then;iErr=10+i;rmse=1.d6;return;endif
			sse=sse+(veCalc-veData(i))**2
		enddo
		rmse=DSQRT(sse/DBLE(nPts));iErr=0
	END SUBROUTINE EvaluateEsd2SystemVE

	SUBROUTINE ConfigureEsd2System(polarCas,partnerCas,alpha,qFit,epsFit,bFit,beta,kijValue,iErr)
		USE EsdParms
		USE GlobConst, only:nmx,bVolCc_mol
		USE BIPs, only:KijMatrix=>Kij,Lij
		Integer, INTENT(IN) :: polarCas,partnerCas
		DoublePrecision, INTENT(IN) :: alpha,qFit,epsFit,bFit,beta,kijValue
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),iLocal
		cas=0;cas(1)=polarCas;cas(2)=partnerCas
		CALL PGLWrapperStartup(2,23,cas,iLocal)
		if(iLocal/=0)then;iErr=1;return;endif
		CALL ApplyRegisteredPartner(partnerCas)
		CALL SetEsd2PolarMode(1,alpha,.FALSE.,iLocal)
		if(iLocal/=0)then;iErr=2;return;endif
		q(1)=qFit;c(1)=1.d0+(qFit-1.d0)/(ESD2_B0/(ESD2_B0-ESD2_k0))
		eokP(1)=epsFit;vx(1)=bFit;bVolCc_mol(1)=bFit
		KijMatrix(1,2)=kijValue;KijMatrix(2,1)=kijValue;Lij(1,2)=beta;Lij(2,1)=beta;iErr=0
	END SUBROUTINE ConfigureEsd2System

	SUBROUTINE EvaluateEsd2PurePsat(polarCas,alpha,qFit,epsFit,bFit,nPts,tData,pCalc,iErr)
		USE GlobConst, only:nmx
		Integer, INTENT(IN) :: polarCas,nPts
		DoublePrecision, INTENT(IN) :: alpha,qFit,epsFit,bFit,tData(nPts)
		DoublePrecision, INTENT(OUT) :: pCalc(nPts)
		Integer, INTENT(OUT) :: iErr
		Integer cas(nmx),i,iLocal,iPsat
		DoublePrecision chemPot(nmx),rhoL,rhoV,uL,uV
		cas=0;cas(1)=polarCas;CALL PGLWrapperStartup(1,23,cas,iLocal)
		if(iLocal/=0)then;iErr=1;pCalc=-1.d0;return;endif
		CALL SetEsd2PolarMode(1,alpha,.FALSE.,iLocal);CALL SetPureCandidate(qFit,epsFit,bFit)
		do i=1,nPts
			CALL PsatEar(tData(i),pCalc(i),chemPot,rhoL,rhoV,uL,uV,iPsat)
			if(iPsat>10 .or. pCalc(i)<=0.d0)pCalc(i)=-1.d0
		enddo
		iErr=0
	END SUBROUTINE EvaluateEsd2PurePsat

	SUBROUTINE EvaluateEsd2ExcessPoint(tKelvin,pMPa,x1,ve,heKJmol,geRT,seJmolK,iErr)
		USE GlobConst, only:nmx,Rgas
		DoublePrecision, INTENT(IN) :: tKelvin,pMPa,x1
		DoublePrecision, INTENT(OUT) :: ve,heKJmol,geRT,seJmolK
		Integer, INTENT(OUT) :: iErr
		Integer i,iLocal
		DoublePrecision x(nmx),fug(nmx),rho,z,a,u,vPure(2),hPure(2),gPure(2)
		DoublePrecision hMix,gMix,hERT
		if(x1<=0.d0 .or. x1>=1.d0)then;iErr=1;ve=0.d0;heKJmol=0.d0;geRT=0.d0;seJmolK=0.d0;return;endif
		do i=1,2
			x=0.d0;x(i)=1.d0
			CALL FugiTP(tKelvin,pMPa,x,2,1,rho,z,a,fug,u,iLocal)
			if(iLocal>10 .or. rho<=0.d0)then;iErr=10+i;return;endif
			vPure(i)=1.d0/rho;hPure(i)=u+z-1.d0;gPure(i)=a+z-1.d0-DLOG(z)
		enddo
		x=0.d0;x(1)=x1;x(2)=1.d0-x1
		CALL FugiTP(tKelvin,pMPa,x,2,1,rho,z,a,fug,u,iLocal)
		if(iLocal>10 .or. rho<=0.d0)then;iErr=20;return;endif
		ve=1.d0/rho-x(1)*vPure(1)-x(2)*vPure(2)
		hMix=u+z-1.d0;hERT=hMix-x(1)*hPure(1)-x(2)*hPure(2)
		gMix=a+z-1.d0-DLOG(z)
		! Excess Gibbs energy excludes the ideal-mixing entropy.  gMix and
		! gPure here are residual Gibbs energies, so GE/RT is their difference.
		! Subtracting sum(x*ln x) would instead report Gibbs energy of mixing.
		geRT=gMix-x(1)*gPure(1)-x(2)*gPure(2)
		heKJmol=hERT*Rgas*tKelvin/1000.d0;seJmolK=Rgas*(hERT-geRT);iErr=0
	END SUBROUTINE EvaluateEsd2ExcessPoint
END MODULE Esd2PolarFit
