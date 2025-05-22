import numpy as np
from ESD import EsdComp,EsdFluid #,CritProps
from GlobConst import Rgas #,pi,avoNum,kB,SQRT2,RgasCal,zeroTol
#from EsdMix import PvpEar #,EsdFluid
# Test Code
#piTest=pi
CO2Esd=EsdComp(909)
print(f"CO2's bVol={CO2Esd.bVolCC_mol} cc/mol and Mw={CO2Esd.Mw} g/mol from LookupEsdParms.")
EthaneEsd=EsdComp(2) # only idDippr argument is required since the defaults are initialized then replaced by lookup.
print(f"Ethane's Tc={EthaneEsd.TcK} from LookupCritParms.")
print(f"ESD Parms for {EthaneEsd.Name} = {EthaneEsd.cShape},{EthaneEsd.eps_kB},{EthaneEsd.bVolCC_mol}.")
nHeptaneEsd=EsdComp(17)
print(f"{nHeptaneEsd.Name}'s Tb(K)={nHeptaneEsd.TbK}  ")
# Below is optional way to specify CritParms if desired.
# If you use calling arguments in the constructor, it's necessary to specify initial values or get error. These are replaced by function.
nOctaneEsd=EsdComp(27,568.7,2.49,0.3996,398.779,0.2569,250.0,111659,15.44,0.6992,114.23,'norml','C8H18','nOCTANE')
print(f"ESD Parms for {nOctaneEsd.Name} = {nOctaneEsd.qShape},{nOctaneEsd.eps_kB},{nOctaneEsd.bVolCC_mol}.")
print(f"Assoc Parms = {nOctaneEsd.assoc.bondVolNm3},{nOctaneEsd.assoc.epsD_kB},{nOctaneEsd.assoc.epsA_kB} ")
WaterEsd=EsdComp(1921)
print(f"ESD Parms for {WaterEsd.Name} = {WaterEsd.qShape},{WaterEsd.eps_kB},{WaterEsd.bVolCC_mol}.")
print(f"Assoc Parms = {WaterEsd.assoc.bondVolNm3},{WaterEsd.assoc.epsD_kB},{WaterEsd.assoc.epsA_kB} ")
LdenEstG_cc=WaterEsd.liqDenG_cc(373)
PvpMpa=WaterEsd.PvpSc(373)
print(f"At Tb(K)=373, {WaterEsd.Name}'s SC density={LdenEstG_cc:.5g} (g/cc) and SCVP = {PvpMpa:.6g} MPa. ") #:.6g gives 6 sig figs.
LdenScG_cc=WaterEsd.liqDenG_cc(298)
LDenEstG_cc=WaterEsd.EsdDenLiqG_cc(298) # quick test for H2O
PvpSc=WaterEsd.PvpSc(298)
print(f"At T(K)=298, {WaterEsd.Name}'s ~ESD density={LDenEstG_cc:.5g} LdenSc={LdenScG_cc:.5g} (g/cc), and PvpSc = {PvpSc:.6g} MPa. ")
WaterFluid=EsdFluid(WaterEsd)
xFrac=np.zeros(2)
xFrac[0]=1
iErrF,rhoLiq,zLiq,aResLiq,uSatL,chemPot,etaLiq=WaterFluid.ChemPoTP(3,xFrac,298,PvpSc)
iErrF,rhoVap,zVap,aResVap,uSatV,chemPot,etaVap=WaterFluid.ChemPoTP(2,xFrac,298,PvpSc)
if(iErrF > 0):iErr=4 #declare error but don't stop. if future iterations give valid fugi(), then no worries
rhoVap=rhoLiq*np.exp( aResLiq+zLiq-1-aResVap+zVap-1 ) #pVpItic, FPE, 501:112236(19), Eq 19.aResVap =/= 0 for carbo acids
PvpItic=rhoVap*Rgas*298*zVap
print(f"At T(K)=298, {WaterEsd.Name}'s ESD itic density={rhoLiq*WaterEsd.Mw:.5g} (g/cc), and PvpItic = {PvpItic:.6g} MPa. ")
iErr,PvpEar,chemPo1,rhoLiq,rhoVap,uSatL,uSatV= WaterFluid.PvpEar(WaterEsd,298)
print(f"At T(K)=298, {WaterEsd.Name}'s ESD ear  density={rhoLiq*WaterEsd.Mw:.5g} (g/cc), and PvpEar  = {PvpEar :.6g} MPa. ")
c7c8=EsdFluid(nHeptaneEsd,nOctaneEsd)
#print(c7c8.components)
c7DenEstG_cc=nHeptaneEsd.EsdDenLiqG_cc(298) 
c8DenEstG_cc=nOctaneEsd.EsdDenLiqG_cc(298) 
EtOHEsd=EsdComp(1102)
xFrac=[0.5, 0.5]
EtOH_W=EsdFluid(EtOHEsd,WaterEsd)
Vmix=(xFrac[0]*nHeptaneEsd.Mw/c7DenEstG_cc+xFrac[1]*nOctaneEsd.Mw/c8DenEstG_cc)
iErr,Pcalc,zFactor,aRes,uRes,chemPo=c7c8.ChemPoTV(True,xFrac,298,Vmix)
print('iErr,Pcalc,zFactor,aRes for c7c8 from ChemPoTV')
print(iErr,Pcalc,zFactor,aRes)
PMPa=0.1
iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq=c7c8.ChemPoTP(1,xFrac,298,PMPa)
print('iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq for c7c8 from ChemPoTP')
print(iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq)
iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq=EtOH_W.ChemPoTP(1,xFrac,298,PMPa)
print('iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq for EtOH+W from ChemPoTP')
print(iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,etaLiq)


