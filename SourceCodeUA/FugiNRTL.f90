!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc PROGRAMED BY AV 06/26/06cccccccccccccccccccccccccccccccccccccccccccccc
!ccccc PURPOSE: CALCULATION OF THERMODYNAMIC PROPERTIES WITH NRTLccccccccc

	SUBROUTINE GetNRTL (nComps,idComp,iErrCode)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=1555)
	character bipFile*234 !,inFile*234
	integer GetBIPs
	DIMENSION idComp(nComps)
	iErrCode=0
!c  note:  bips are passed back through USEd BIPs/
	bipFile=TRIM(PGLinputDir)//'\BipNRTL.txt' ! // is the concatenation operator
     iErrCode=GetBIPs(bipFile,idComp,nComps)
	do i=1,nComps
		tKmin(i)=0.4d0*Tc(i)
	enddo 

	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)

	return                      
	END
	!PROGRAMED BY AV 06/22
	!PURPOSE: CALCULATION OF VAPOR PRESSURES AND FUGACITY COEFFICIENTS BY NRTL ACTIVITY COEFFICIENT MODEL

	subroutine FuNRTL(tK,pMpa,moFrac,nComps,LIQ,FUGC,zFactor,iErr)
!c  use NRTL style of EXPression
	USE GlobConst ! includes vLiq when you call CritParms first.
	USE BIPs
	USE VpDb ! for VpCoeffs
	implicit doubleprecision(A-H,K,O-Z)
	double precision moFrac(Nmx)
	integer kComp
	dimension fugc(nComps),xsGamma(nComps)
	DIMENSION omega(Nmx,Nmx),sumDenom(Nmx),sumNumer(Nmx)
	iErr=0
	!DEBUG=.TRUE.
	vLiqMix=0
	do iComp = 1,nComps
		vLiqMix=vLiqMix+moFrac(iComp)*vLiq(iComp)
	enddo
	zFactor=1
	fugc(iComp)=1
	if(LIQ==0)return ! ideal gas for vapor means... That's all folks!
	zFactor=pMpa*vLiqMix/(rGas*tK)

	sumLog=0
	do jComp=1,nComps
		sumDenom(jComp)=0
		sumNumer(jComp)=0
		do kComp=1,nComps
			omega(kComp,jComp)=EXP(-xsAlpha(kComp,jComp)*xsTau(kComp,jComp)/tK)
			sumDenom(jComp)=sumDenom(jComp)+moFrac(kComp)*omega(kComp,jComp)
			sumNumer(jComp)=sumNumer(jComp)+moFrac(kComp)*omega(kComp,jComp)*xsTau(kComp,jComp)/tK
		enddo
		sumLog=sumLog+moFrac(jComp)*sumDenom(jComp)
	enddo

	xsFreeEn =0
	do kComp=1,nComps
		bigSumK=0.d0
		do jComp=1,nComps
			bigSumK=bigSumK+moFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)*( xsTau(kComp,jComp)/tK-sumNumer(jComp)/sumDenom(jComp) )
		enddo
		xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
		xsFreeEn =xsFreeEn+xsLnGamma*moFrac(kComp)
		fugc(kComp)=xsLnGamma
		xsGamma(kComp)=exp(xsLnGamma)
	enddo
	do kComp=1,nComps
		CALL GetVp(nComps,ID,iErrCode)
		pSat=exp(vpCoeffs(kComp,1)+vpCoeffs(kComp,2)/tK+vpCoeffs(kComp,3)*DLOG(tK)+vpCoeffs(kComp,4)*tK**vpCoeffs(kComp,5))/1000000
		fugc(kComp)=fugC(kComp)+LOG(pSat/pMPa)
	enddo
	return
	end
