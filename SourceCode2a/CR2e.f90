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
			call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
			ZC_mix=Z
			PC_mix=PMpa
			etaC_mix=eta
			return
		endif
	enddo
	TC_mix=tKelvin
	VC_mix=vTotCc/totMoles
	call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
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
    IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(N=2)
	DIMENSION gmol(NMX),xFrac(NMX),FUGC(NMX),VC(NMX) !,ier(11)
	!DIMENSION GJ_Jac(2,2),B_GJ(2,1)
	DIMENSION KIJ(NMX,NMX),KTIJ(NMX,NMX),HIJ(NMX,NMX),HTIJ(NMX,NMX),xsTau(NMX,NMX),xsTauT(NMX,NMX),xsAlpha(NMX,NMX)
	doubleprecision mShape(NMX),Acen_pure
	character*200 ErrMsg(4) !,outfile*100
	LOGICAL LOUDER
	!COMMON/ETA2/ETA
	COMMON/BIPs_SPEAD/aBipAD,aBipDA
	COMMON/BIPs/KIJ,KTIJ,HIJ,HTIJ,xsTau,xsTauT,xsAlpha	 !warning: dim of aBip should be maxTypes, not nmx.
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
	call FuVtot(isZiter,tKelvin,1/rho,gmol,NC,FUGC,Zbase,iErrVtot)
    if(iErrVtot.ne.0)iErrCode=4
    if(iErrCode.ne.0)goto 86
	call FuVtot(isZiter,tKelvin,1/rhoHi,gmol,NC,FUGC,ZHi,iErrVtot)
	call FuVtot(isZiter,tKelvin,1/rhoLo,gmol,NC,FUGC,ZLo,iErrVtot)
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
	        call FuVtot(isZiter,tKelvin,1/rho,gmol,NC,FUGC,Zbase,iErrVtot)
		    call FuVtot(isZiter,tKelvin,1/rhoHi,gmol,NC,FUGC,ZHi,iErrVtot)
		    call FuVtot(isZiter,tKelvin,1/rhoLo,gmol,NC,FUGC,ZLo,iErrVtot)
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
	call FuVtot(isZiter,tBest,1/rhoBest,gmol,NC,FUGC,Z,iErrVtot)
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
	!half=1.d0/2.d0
	call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
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
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
        if(iErrVtot.ne.0)iErr=1
		x_higherT=x

		tKelvin=T_old-T_old/stepT
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
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
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
        if(iErrVtot.ne.0)iErr=iErr+100
		x_higherV=x   
		
		vTotCc=vTotCc_old-vTotCc_old/stepV
		call FuVtot(isZiter,tKelvin,vTotCc,gmol,NC,FUGC,Z,iErrVtot)
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

