MODULE Assoc  ! This module is site-based (similar to Group Contribution (GC) basis). Sums are over sites per molecule then over molecules.
	USE GlobConst
	implicit NONE
	integer maxTypes,nsx,maxTypesGlobal,localPool	
	PARAMETER (maxTypes=44,nsx=maxTypes,maxTypesGlobal=999) 
    parameter(localPool=9999) !this must be long enough to cover all SpeadMd site types (e.g. 1401=methanol hydroxy)
	!maxTypes is the max # of site types for all molecules. maxTypes > sum^NC(count(Types))	!nsx = maxTypes (dunno why redundant), nbx=Max Bonds, 
	DoublePrecision eHbKcal_mol(nmx,maxTypes),eDonorKcal_mol(nmx,maxTypes),eAcceptorKcal_mol(nmx,maxTypes)
	DoublePrecision bondVolNm3(nmx,maxTypes)
	Integer idType(NMX,maxTypes),nTypes(NMX),nTypesTot !nTypesTot is really the sum of nTypes(i) because the same type on a different molecule is treated distinctly (?)
	Integer nDonors(nmx,maxTypes),nAcceptors(nmx,maxTypes),nDegree(nmx,maxTypes)
	Integer localType(localPool),idLocalType(maxTypes)	!these are to accumulate site lists in Wertheim so the aBipAd,aBipDa arrays don't get too large. cf. AlphaSp
	DoublePrecision aBipAD(maxTypes,maxTypes),aBipDA(maxTypes,maxTypes) !association bips
	DoublePrecision XA(NMX,maxTypes),XD(NMX,maxTypes),XC(NMX,maxTypes),cvAssoc 
	!localType is an index of just the types occuring in the current mixture.  e.g. localType(101)=1 means that the 1st type encountered during reading the dbase was type 101.
	!idLocalType points back to localType for double-linking. E.g. idLocalType(1)=101.
END MODULE Assoc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SpeadParms
	USE GlobConst
	USE Assoc
	DoublePrecision zRefCoeff(NMX,5),a1Coeff(NMX,5),a2Coeff(NMX,5),vMolecNm3(NMX),tKmin(NMX),rMw(NMX)
	Integer nTptCoeffs
	!nnx = Total number of united-atom sites in a molecule.
END MODULE SpeadParms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
