MODULE ModelSettings ! Keep everything in one place so adding a model is easier. See also QUERYMODEL(PGLDLL.f90).
	parameter(nModels=20)
	integer nParInert(nModels),nParAssoc(nModels),nParPolar(nModels),nParMix(nModels) ! 
	Character*15 EosName(nModels)
	!              1     2       3       4          5          6         7           8              9        10       
	data EosName/'PR','ESD96','PRWS','ESD-MEM2','SPEADMD','Flory-MEM2','NRTL','SpeadGamma-MEM2','SPEAD11','PcSaft',&
			'tcPRq','EgcESD','EgcEsdTb','TffSPEAD','EgcPcSaft','EgcPcSaft(Tb)','tcPR-GE(W)','ESD2','LsgMem2','SptPcSaft'/
	!             11    12      13         14          15          16             17        18       19        20       
	!               1 2 3 4 5  6 7 8  9 10 11 12 13 14 15 16 17 18 19 20
	data nParInert /0,3,0,3,11,0,0,11,0, 3, 4, 3,11, 0, 3, 3, 3, 3, 3, 3/	! e.g. m, sigma, eps/kB. for PcSaft & ESD
	data nParAssoc /0,2,0,2, 0,0,0, 0,0, 2, 0, 2, 0, 0, 2, 2, 2, 2, 2, 2/	! GC&Tff EOS's have zero adj par's 
	data nParPolar /0,0,0,0, 0,0,0, 0,0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/	! GC&Tff EOS's have zero adj par's 
	data nParMix   /1,1,3,1, 1,1,1, 1,1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1/	! Mostly kij.
	! nParMix:
    !    PcSaft(10): lij(Tang-Gross,2010) and kij_assoc are also defined but we don't use them in PGLDLL for now. JRE 20230822
    !    tcPRq(11): kij (for now, could include betaij later). JRE 20210531if(iEosOpt==10)then
    !    PrLorraine(17): age0(1,2) and age0(2,1)
    !    LsgMem2(19): age0(1,2) and age0(2,1) for LSG. MEM2 parameters are assumed transferable.
END MODULE ModelSettings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GlobConst
	!SAVE
	!PUBLIC
	USE ModelSettings, only: nModels, EosName ! defined in UaEosTools.f90
	Implicit NONE
	DoublePrecision pi,twoPi,fourPi,avoNum,Rgas,RgasCal,kB,half,zeroTol
	Integer nmx,nCritSet
	PARAMETER (nmx=55,pi=3.14159265359796d0,twoPi=2.d0*pi, fourPi = 4.d0*pi, half=0.5d0)
	PARAMETER (avoNum=602.214076d0,kB=0.01380649D0,Rgas=avoNum*kB,RgasCal=Rgas/4.184d0,zeroTol=1.D-12)
	          !avoNum[=]cm3/(nm3*mol), kB[=]MPa.nm3/K. cf. PGL6ed, Table 6.1
	PARAMETER (nCritSet=1739) ! This is the number of compounds that should be found in LoadCritParmsDb. 
	!          Change this parameter if you add more compounds.	It must be consistent or you will get a LOAD error.
	!          https://www.nist.gov/si-redefinition 
    !nmx is the max allowed number of Compounds
	!integer :: idComp(nmx),nsTypes(nmx),IDs(nsx),IDsBase(nmx,nsx),siteNum(nmx,maxTypes)
	!integer :: nComps, nsTypesTot,iTPT,iFlagFF,nNormGrid

	DoublePrecision :: Tc(nmx),Tb(nmx), Pc(nmx), ACEN(nmx), Zc(nmx),TcEos(nmx), PcEos(nmx), ACENEos(nmx), ZcEos(nmx), rMw(nmx)
	DoublePrecision bVolCc_mol(nmx)	!vdW molar volume of molecule. 
	DoublePrecision solParm(nmx)    !solubility parameter in (J/cm3)^0.5=MPa^0.5
	DoublePrecision vLiq(nmx)		!liquid molar volume in cm3/mol
	DoublePrecision tKmin(nmx)		!model's minimum recommended temperature, especially for P^vp. Only warning. e.g., SLE OK.
	DoublePrecision CpRes_R, CvRes_R, cmprsblty !cmprsblty=(dP/dRho)T*(1/RT)	= Z+rho*dZ/dRho
	character*255 masterDir,PGLinputDir,dumpFile
	character*30 NAME(nmx)
	character*5 class(nmx) ! Allowed: norml,heavy,polar,assoc,Asso+,gases,siloa,salty,ormet,metal,inorg (cf. PGL6edClasses.xls)
    Logical, SAVE :: LOUD,CheckDLL,isTDE   !LOUD=.TRUE. means writing debug info. CheckDLL=.TRUE. means direct debug info to DebugDLL.txt; isTDE means the call is coming from TDE
	LOGICAL DEBUG, isESD, isTPT, isPcSaft 
	integer ID(nmx), idCas(nmx), idTrc(nmx), iEosOpt, initEos, dumpUnit 
	DoublePrecision etaPass 
    DoublePrecision etaMax  !each EOS has a max value for eta, e.g. PR,TPT: etaMax=1-zeroTol. Set in the Get_ function for the EOS
	DoublePrecision etaPure(nmx) !store eta for each compound at P=0.1MPa and tKmin(K). 
	Character*9  form600 
	Character*12 form601 
	Character*13 form602 
	Character*11 form610 
	Character*14 form611 
	Character*15 form612 
	Character*15 form613 
	data form600,form601,form602/'(12E12.4)','(i7,12E12.4)','(2i7,12E12.4)'/
	data form610,form611,form612,form613/'(a,12E12.4)','(a,i8,12E12.4)','(a,2i8,12E12.4)','(a,3i8,12E12.4)'/ 
	!LOUD = .TRUE.		  !!!!!!!!!!!!!!! YOU CAN'T SET VARIABLES IN A MODULE, ONLY PARAMETERS AND DATA !!!!!!!!!!!!!
	!LOUD = .FALSE.
contains
	integer function SetNewEos(newEosOpt)
	!returns 0 if no error.
	!USE GlobConst ! eliminated because GlobConst contains SetNewEos.
	implicit none
	integer newEosOpt
	SetNewEos=0
	isESD=.FALSE.
	isTPT=.FALSE.
	isPcSaft=.FALSE.
	iEosOpt=newEosOpt
	return
	end function SetNewEos
END MODULE GlobConst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Assoc  ! This module is site-based (similar to Group Contribution (GC) basis). 1st Sums sites per molecule then  molecules.
	USE GlobConst, only:nmx,isTPT,isESD,iEosOpt,LOUD,dumpUnit,dumpFile ! nmx is the maximum number of compounds, typically 55. 
	implicit NONE
	integer maxTypes,nsx,maxTypesGlobal,localPool	
	PARAMETER (maxTypes=44,nsx=maxTypes,maxTypesGlobal=999) 
    parameter(localPool=9999) !this must be long enough to cover all SpeadMd site types (e.g. 1401=methanol hydroxy)
	!maxTypes is the max # of site types for all molecules. maxTypes > sum^NC(count(Types))	
	!nsx = maxTypes (dunno why redundant), nbx=Max Bonds, 
	DoublePrecision eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes)
	DoublePrecision eHbKcal_mol(nmx,maxTypes) ! eHbKcal_mol=(eDonorKcal_mol+eAcceptorKcal_mol)/2
	DoublePrecision bondVolNm3(nmx,maxTypes)
	Integer idType(nmx,maxTypes),nTypes(nmx),nTypesTot !nTypesTot = SUM( nTypes(i) ). Same type on diff molecule is distinct
	Integer nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),nDegree(nmx,maxTypes)
	Integer localType(localPool),idLocalType(maxTypes)	!To accumulate site lists in Wertheim to minimize aBipAd,aBipDa size
	DoublePrecision aBipAD(maxTypes,maxTypes),aBipDA(maxTypes,maxTypes) !association bips
	DoublePrecision XA(nmx,maxTypes),XD(nmx,maxTypes),XC(nmx,maxTypes),cvAssoc
	DoublePrecision rLnGamRep(nmx),rLnGamAtt(nmx),rLnGamAssoc(nmx),bondRate(nmx,maxTypes) ! These replace .../FugiParts/
	LOGICAL LOUDERWert
	 
	!localType is an index of just the types occuring in the current mixture.  
	!e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!idLocalType points back to localType for double-linking. E.g. idLocalType(1)=101.
contains
	Subroutine RdfCalc(rdfContact,dAlpha,eta)
	USE GlobConst
	DoublePrecision denom,eta,denom2,void,void2,void4,dLng,rdfContact,d2Lng,d2g,dAlpha,dg_deta,dAlph_deta
	Integer iRdfOpt,iErrCode
	! Input:
	! eta = packing fraction
	! iEosOpt, isESD, isTPT cf. GlobConst
	! Output:
	! rdfContact= radial distribution function at contact
	! dAlpha=dLn(alpha)/dLn(rho)
	! Miscellaneous:
	! alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/deta) = 1+(eta/rdf)*(dRdf/deta)
	! dLng = (eta/rdf)*(dRdf/deta)     
	! iRdfOpt	- option for characterizing the radial distribution funciton at contact.	
	!			= 0, not specified => error
	!			= 1, ESD form
	!			= 2, Carnahan-Starling form
	!			= 3, activity coefficient model form (rdf=universal constant=g(eta=0.4)
	!			= 4, ESD form, geometric association rule.
	iRdfOpt=0					!ESD non-geometric form
	if(isESD)then
		iRdfOpt=4	!ESD geometric form
	elseif(isTPT)then
		iRdfOpt=2	!CS form
	elseif(iEosOpt==6)then
		iRdfOpt=3	!activity model
	elseif(iEosOpt==19)then
		iRdfOpt=2	!LSG-MEM2 model
	endif

	if(iRdfOpt==0)then
		iErrCode=11
		if(LOUD)write(dumpUnit,*) 'Error in RdfCalc: iRdfOpt not specified'
		if(LOUD)write(dumpUnit,*) '1=ESD form, 2=CS form, 3=(rdf=1), 4=ESD+geometric association'
		return
	endif
	!alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/deta) = 1+(eta/rdf)*(dRdf/deta)
	!dLng = (eta/rdf)*(dRdf/deta) = dLng/dLneta
	denom=1.d0-1.9d0*eta 
	denom2=denom*denom
	void=1.d0-eta
	void2=void*void
	void4=void2*void2    

	dLng=1.9d0*eta/denom
	if(iRdfOpt==2)dLng=eta*( 3.d0/void-1.d0/(2.d0-eta) )
	if(iRdfOpt==3)dLng= -1.d0
	rdfContact=1/denom
	if(iRdfOpt==2)rdfContact=(1.d0-eta/2.d0)/void/void2
	if(iRdfOpt==3)rdfContact=1.d0
	d2Lng=1.9d0*eta/denom2
	d2g=d2Lng-dLng
	dAlpha=1+dLng
	dg_deta=1.9d0/denom2
	dAlph_deta=dg_deta
	if(iRdfOpt==2)dg_deta=(2.5d0-eta)/void4
  	if(iRdfOpt==2)dAlph_deta=3.d0/void2-2.d0/((2.d0-eta)*(2.d0-eta))
	return
	END	Subroutine RdfCalc
END MODULE Assoc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CritParmsDb
	Integer ndb
	Parameter (ndb=3000)
	character*30, STATIC:: NAMED(ndb)	!LoadCrit() loads ParmsCrit database.
	character*5, STATIC:: classDb(ndb) 
	character*11, STATIC:: formDb(ndb) 
	Integer, STATIC:: IDnum(ndb),CrIndex(99999),idCasDb(ndb),nDeckDb ! e.g. TCD(CrIndex(2)) returns Tc of ethane. 
	DoublePrecision, STATIC:: TCD(ndb),PCD(ndb),ACEND(ndb),TbD(ndb),ZCD(ndb),solParmD(ndb),rMwD(ndb),vLiqD(ndb) 
	! LoadCrit uses CrIndex to facilitate lookup. TCD(ndb)=8686. CrIndex()=ndb initially.
END MODULE CritParmsDb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE BIPs	 ! For molecular level binary interaction parameters.
	USE GlobConst, only: nmx 
	integer nConverged,maxPts,nPtsBipDat
	parameter(maxPts=1777) ! this defines the max allowed # of experimental data points in a single binary system.
	DoublePrecision KIJ(nmx,nmx),KTIJ(nmx,nmx),kETAij(nmx,nmx) !usual dispersive BIPs & k^eta_ij(for speadmd)
	DoublePrecision KS0IJ(nmx,nmx),KS1IJ(nmx,nmx)              !entropic BIPs & k^eta_ij
	DoublePrecision HIJ(nmx,nmx),HTIJ(nmx,nmx) !molecular hBonding BIPs for ESD. (Spead aBipAd,aBipDa are site based.)
	DoublePrecision Lij(nmx,nmx) !covolume adjustment.  bVolMix=sum(sum(xi*xj*bij)); bij=(1-Lij)*(bi+bj)/2
	DoublePrecision xsTau(nmx,nmx),xsTauT(nmx,nmx),xsAlpha(nmx,nmx)	!this is for the PRWS/xsNRTL mixing rule.
	DoublePrecision tDat(maxPts),PDAT(maxPts),XDAT(maxPts),YDAT(maxPts) !,deviate(maxPts) 
	!can't include Deviate() cuz it's an argument for LmDif.
	DoublePrecision pDatMin,pDatMax
	Integer iDat(maxPts)  ! sometimes need to indicate whether the data are for comp1 or comp2. e.g. SLE.
END MODULE BIPs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE VpDb
	USE GlobConst, only:nmx !,dumpUnit
	IMPLICIT NONE !DoublePrecision(A-H,O-Z)
	Integer nVpDb
	PARAMETER(nVpDb=2475)
	Integer, STATIC:: IDnum(nVpDb),NUMCOEFFD(nVpDb) ,indexVpDb(9999)
	DoublePrecision, STATIC:: rMINTD(nVpDb) ,VALMIND(nVpDb) ,rMAXTD(nVpDb),VALMAXD(nVpDb),AVGDEVD(nVpDb),vpCoeffsd(nVpDb,5)
	DoublePrecision vpCoeffs(nmx,5)
END MODULE VpDb
