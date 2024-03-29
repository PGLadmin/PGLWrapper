!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc PROGRAMED BY AV 06/26/06cccccccccccccccccccccccccccccccccccccccccccccc
!ccccc PURPOSE: CALCULATION OF THERMODYNAMIC PROPERTIES WITH NRTLccccccccc

	SUBROUTINE GetNRTL (nComps,idComp,iErrCode)
	USE GlobConst
	IMPLICIT DOUBLEPRECISION(A-H,K,O-Z)
	PARAMETER(ndb=1555)
	character*234 bipFile !,inFile
	integer GetBIPs
	DIMENSION idComp(nComps)
	iErrCode=0
!c  note:  bips are passed back through USEd BIPs/
	bipFile=TRIM(PGLinputDir)//'\BipNRTL.txt' ! // is the concatenation operator
     iErrCode=GetBIPs(bipFile,idComp,nComps)
	do i=1,nComps
		tKmin(i)=0.3d0*Tc(i)
	enddo

	TcEos(1:nComps)=Tc(1:nComps) !This EOS is consistent with experimental values for critical properties.
	PcEos(1:nComps)=Pc(1:nComps)
	ZcEos(1:nComps)=Zc(1:nComps)

	return                      
	END
	!PROGRAMED BY AV 06/22
	!PURPOSE: CALCULATION OF VAPOR PRESSURES AND FUGACITY COEFFICIENTS BY NRTL ACTIVITY COEFFICIENT MODEL

	subroutine FuNRTL(tKelvin,pMpa,xFrac,nComps,LIQ,FUGC,rhoMol_cc,zFactor,aRes,uRes,iErr)
!c  use NRTL style of EXPression
	USE GlobConst ! includes vLiq when you call CritParms first.
	USE BIPs
	USE VpDb ! for VpCoeffs
	implicit doubleprecision(A-H,K,O-Z)
	double precision xFrac(Nmx)
	integer kComp
	dimension fugc(nComps),xsGamma(nComps)
	DIMENSION omega(Nmx,Nmx),sumDenom(Nmx),sumNumer(Nmx)
	iErr=0
	NC=nComps
	!DEBUG=.TRUE.
	zFactor=1
	fugc(1:NC)=1
	aRes=0
	uRes=0
	if(LIQ==0)return ! ideal gas for vapor means... That's all folks!
	expo=2.d0/7.d0
	vLiqMix=0
	do iComp = 1,nComps
		phi=(1-tKelvin/Tc(iComp))**expo - (1-298/Tc(iComp))**expo  !PGL6ed Eq. 5.2-6, aka Yamada-Gunn adaptation of Rackett eq.
		vSat=vLiq(iComp)*( 0.29056+0.08775*acen(iComp) )**phi
		vLiqMix=vLiqMix+xFrac(iComp)*vSat
	enddo
	rhoMol_cc=1/vLiqMix
	bMix=SUM( xFrac(1:NC)*bVolCC_mol(1:NC) )
	eta=rhoMol_cc*bMix
	zFactor=pMpa*vLiqMix/(rGas*tKelvin)

	sumLog=0
	do jComp=1,nComps
		sumDenom(jComp)=0
		sumNumer(jComp)=0
		do kComp=1,nComps
			omega(kComp,jComp)=EXP(-xsAlpha(kComp,jComp)*xsTau(kComp,jComp)/tKelvin)
			sumDenom(jComp)=sumDenom(jComp)+xFrac(kComp)*omega(kComp,jComp)
			sumNumer(jComp)=sumNumer(jComp)+xFrac(kComp)*omega(kComp,jComp)*xsTau(kComp,jComp)/tKelvin
		enddo
		sumLog=sumLog+xFrac(jComp)*sumDenom(jComp)
	enddo

	xsFreeEn =0
	GmixIdSoln=0
	do kComp=1,nComps
		bigSumK=0.d0
		do jComp=1,nComps
			bigSumK=bigSumK+xFrac(jComp)*omega(kComp,jComp)/sumDenom(jComp)*( xsTau(kComp,jComp)/tKelvin &
			                                                                                    -sumNumer(jComp)/sumDenom(jComp) )
		enddo
		xsLnGamma=sumNumer(kComp)/sumDenom(kComp)+bigSumK
		xsFreeEn =xsFreeEn+xsLnGamma*xFrac(kComp)
		GmixIdSoln =GmixIdSoln+LOG(xFrac(kComp))*xFrac(kComp)
		fugc(kComp)=xsLnGamma
		xsGamma(kComp)=exp(xsLnGamma)
	enddo
	aRes=GmixIdSoln+xsFreeEn !assuming Vxs=0
	do kComp=1,nComps
		CALL GetVp(nComps,ID,iErrCode)
		pSat=exp(vpCoeffs(kComp,1)+vpCoeffs(kComp,2)/tKelvin+vpCoeffs(kComp,3) &
		                                                      *DLOG(tKelvin)+vpCoeffs(kComp,4)*tKelvin**vpCoeffs(kComp,5))/1000000
		fugc(kComp)=fugC(kComp)+LOG(pSat/pMPa)
	enddo
	return
	end
