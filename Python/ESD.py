import numpy as np
import math # for sqrt, log, ... 
from sgtpy.vrmie_pure.psat_saft import psat
from sgtpy.vrmie_pure.ideal import  daideal_drho #, aideal,d2aideal_drho 
#import pandas as pd
#import os
#from ..constants import kb, Na
#from scipy.optimize import root, minimize_scalar
from GlobConst import avoNum,Rgas,RgasCal,zeroTol,kB,PGLInputDir #,pi,SQRT2,
from chempy.util.parsing import formula_to_composition #JRE: to convert chemical formulas to atomCounts in LookupCritParms
#from chempy import Substance
import keyboard
def pauseCheck(force=None):
    if(force==1):input("Press enter to continue.")
    if keyboard.is_pressed('Esc'):
        input("Press enter to continue.")

class CritProps: # Analogous to Car class in Python Crash Course Book.
    """A class for holding the most commonly needed properties of chemical compounds. """
    def __init__(self,idDippr,TcK=190.8686,PcMPa=4.8,acen=0,TbK=111,Zc=0.29,Tmin=88,idCas=1234,solp=6.6,rho298=0.5,Mw=16.04,Class='norml',FORM='CH4',Name='Methane'): # Class "Constructor"
        # NOTE: because this constructor lists all the properties, we must specify all these values whenever we invoke this class.
        self.idDippr=idDippr
        self.TcK=TcK
        self.PcMPa=PcMPa
        self.Zc=Zc
        self.acen=acen
        self.TbK=TbK
        self.Tmin=Tmin
        self.idCas=idCas
        self.solp=solp
        self.rho298=rho298
        self.Mw=Mw
        self.JREClass=Class
        self.FORM=FORM
        self.Name=Name
        self.bAssoc=False # Initialize assuming the compound is not associating.
    def LookupCritParms(self,idDippr):
        # EsdParms={} # a dictionary 
        inFile=PGLInputDir+"ParmsPrTcJaubert.txt"
        #inFile="ParmsPrTcJaubert.txt"
        with open(inFile,encoding='utf-8') as f:
            lines=f.readlines() # reads the entire file as lines.
            # q EOK(K)   b(cc/mo)   KadNm3  eDonEpsK    eAccEpsK    nDegree #eDs    #eAs    idCas   Name
            bThere=False
            i=0
            for line in lines:
                i+=1
                if i==1: continue # skip the first line since its just a header
                part=line.split() # separates out the  
                if int(part[0]) == idDippr:
                    bThere=True
                    self.idDippr=idDippr
                    self.TcK=float(part[1]) # All parts are read as strings. Must convert manually after split().
                    self.PcMPa=float(part[2])
                    self.acen=float(part[3])
                    self.TbK=float(part[4])  
                    self.Zc=float(part[9])  
                    self.Tmin=float(part[10])  
                    self.idCas=int(part[11])  
                    self.solp=float(part[12])  
                    self.rho298=float(part[13])  
                    self.Mw=float(part[14])  
                    self.JREClass=(part[15])  
                    self.FORM=(part[16])
                    try:
                        #f=Substance.from_formula(self.FORM).composition
                        if(self.FORM=='D2'): self.FORM='H2'
                        if(self.FORM=='D2O'): self.FORM='H2O'
                        f = formula_to_composition(self.FORM)
                        self.atomCounts = f #dictionary where atomic number gives number count.
                        self.CHBCFINOSS=np.array([f.get(6,0),f.get(1,0),f.get(35,0),f.get(17,0),f.get(9,0),f.get(53,0),f.get(7,0),f.get(8,0),f.get(16,0),f.get(14,0)] )
                        nCarb=self.CHBCFINOSS[0]
                    except:
                        print("ESD.CritProps: FORM=",self.FORM)
                        self.atomCounts={}
                        self.CHBCFINOSS=np.zeros(8)
                        #composition = Substance.from_formula('D2O').composition
                        #pauseCheck(1)
                    # Outputs: {H': 2, 'O': 1} for H2O. .e.g., nCarb=atomCounts.get('C',0) #returns 0 if no carbons, 
                    #Rgas=kB*avoNum
                    self.Name=(part[17]) 
                    self.bAssoc=False # Initialize assuming the compound is not associating.
                    if self.JREClass=="assoc" or self.JREClass=="Asso+"or self.JREClass=="polar": self.bAssoc=True # Need to include polar for mixes like acetone+chloroform or MEM2 returns immediately.
                    print(f"Compound {idDippr}, {self.Name} found in ParmsPrTcJaubert.txt! bAssoc={self.bAssoc}")
                    break
        return bThere
    def PvpSc(self,tKelvin): #JRE: ShortCut Vapor Pressure Equation. See also PvpEar which uses Equal Area Rule (aka. Maxwell Construction) to solve for Pvp of EOS.
        # log10(Pvp/Pc)=7/3*(1+acen)*(1-T/Tc)
        PvpSc=self.PcMPa*10**( 7*(1+self.acen)/3*(1-self.TcK/tKelvin) )
        return PvpSc
    def liqDenG_cc(self,tKelvin):
        # https://trc.nist.gov/TDE/TDE_Help/Eqns-Pure-DensityLG/Rackett.htm
        # rho_satL = rhoc*Zc^[ -(1-Tr)^(2/7) ] = rhoc*(Pc*Mw/rhoc*R*Tc)^[...] ~EL2ed 9.40
        # ln(rhoL) = ln(rhoc)+[...]*ln(Zc)
        # ln(rhoc) = ln(rho298)-[...@298]*ln(Zc)
        # ln(rhoL) = ln(rho298)+{ [...]-[...@298] }*ln(Zc)
        # [...]-[...@298] = -(1-T/Tc)^2/7 + (1-298/Tc)^2/7 = expo
        expo=(1-298/self.TcK)**(2/7)-(1-tKelvin/self.TcK)**(2/7)
        liqDenG_cc=self.rho298*self.Zc**expo 
        return liqDenG_cc
#################################################################################################
class Assoc: # We can use this as a subclass, e.g., a class within the EsdComp class. Analogous to the Battery class in CrashCourse book.
    """For describing association parameters. Works for ESD or SPEADMD."""
    # This is Ex ~9.6-9.8 in Python Crash Course. 
    def __init__(self): # this constructor has only the self argument, so we don't need to specify anything when invoking.
        # Omitting the calling arguments means you can't pass them in when invoking the class. OK, we use Lookup...
        self.nTypes=2 #This is correct for ESD only.
        self.nDsites=np.zeros(self.nTypes,dtype=int) #nDsites
        self.nAsites=np.zeros(self.nTypes,dtype=int) #nAsites
        self.bondVolNm3=np.zeros(self.nTypes) #bondVolNm3
        self.epsD_kB=np.zeros(self.nTypes) #epsD_kB
        self.epsA_kB=np.zeros(self.nTypes) #epsA_kB
        self.nDegree=np.zeros(self.nTypes,dtype=int) #nDsites
        self.idType=np.zeros(self.nTypes,dtype=int)
        self.nTypes=1 # It's necessary to set as 2 when initializing or np converts it to a scalar, causing an error in MEM2 because it is general. 
        #But really, ESD allows only one type per molecule.
#################################################################################################
class EsdComp(CritProps): # Note how citing CritProps as an argument generates "inheritance."
    """ Build the EsdComp Class "inheriting" the CritProps class. # Ex ~9.12  """ 
    def __init__(self,idDippr,TcK=190.8686,PcMPa=4.6,acen=0,TbK=111,Zc=.29,Tmin=85,idCas=74828,solp=12,rho298=.42,Mw=16,Class='norml',FORM='CH4',Name='METHANE'):
        # setting default values means that these arguments are optional at invocation
        super().__init__(idDippr,TcK,PcMPa,acen,TbK,Zc,Tmin,idCas,solp,rho298,Mw,Class,FORM,Name) # "super()" refers to CritProps
        bThere=False
        self.EOS='ESD'
        self.critical=False
        if self.TcK==190.8686: bThere=self.LookupCritParms(idDippr) # Default Tc means you need to lookup.
        self.assoc=Assoc() #bondVolNm3,epsD_kB,epsA_kB,nTypes,nDsites,nAsites # Example of invoking a subclass
        # if bThere==False:
        # print(f"Compound {self.Name} not found in ParmsEsdMEM2.txt! Estimating parms from Tc,Pc,acen.")
        Wci=self.acen
        if(Wci < 0):
            print("EsdComp constructor: acen < 0. Setting acen=0.")
            Wci=0 #e.g. mercury causes qShape < 0. 
        cShape=1+3.535*Wci+ 0.533*Wci*Wci
        #RooTCinv=1/math.sqrt(cShape)
        k1=1.7745
        qShape=1+(cShape-1)*4/(4-1.9)
        RootQinv=1/math.sqrt(qShape)
        #ZcEsd=1/3+RooTCinv*(.0384+RooTCinv*(-.062+RooTCinv*(0.0723-.0577*RooTCinv)))
        ZcEsd=( 1+RootQinv*(.1451+RootQinv*(-.2046+RootQinv*(0.0323-0*RootQinv))) )/3
        atemp=9.5*qShape*1.9+4*cShape*k1-k1*1.9
        quadB=k1*1.9*ZcEsd+3*atemp
        sqArg=quadB*quadB+4*atemp*(4*cShape-1.9)*(9.5*qShape-k1)/ZcEsd
        Bc=ZcEsd*ZcEsd*( -quadB+math.sqrt(sqArg) )/( 2*atemp*(4*cShape-1.9) )
        Yc=ZcEsd*ZcEsd*ZcEsd/(atemp*Bc*Bc)
        rLnY1=math.log(Yc+1.0617)
        self.cShape=cShape
        self.qShape=qShape
        self.eps_kB=self.TcK*rLnY1
        self.bVolCC_mol=Rgas*self.TcK/self.PcMPa*Bc
        self.ZcEos=ZcEsd
        self.TcEos=self.TcK
        self.PcEos=self.PcMPa
        self.etacEos=Bc/ZcEsd # bP/RT = eta*Z
        self.etaMax=1/1.9
        bThere=self.LookupEsdParms(self.idDippr) # Replace "exact" values if there, but always compute etac and ZcEsd for PvpEar.
        if(bThere):print(f'{self.Name} found in EsdParms! qShape={qShape}')
        # For SGT Computations
        cii=self.cii_correlation() # this could be an array of polynomial(T) coefficients(?), but assume const for now.
        self.cii = np.array(cii, ndmin=1)
        self.Tc=self.TcK        # SgtPy requires Tc in Psat function.
        ######                  End of EsdComp() constructor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def cii_correlation(self, overwrite=False):
        """
        cii_corelation()
        Method that computes the influence parameter of coarsed-grained
        molecules
        AIChE Journal, 62, 5, 1781-1794 (2016)
        Eq. (23)
        Parameters
        ----------
        overwrite : bool
            If true it will overwrite the actual influence parameter.
        Returns
        -------
        cii : float
            correlated influence parameter [J m^5 / mol^2]
        """
        compClass=self.JREClass
        nCarb=(self.CHBCFINOSS[0])
        nSilox=int(  np.sqrt( self.CHBCFINOSS[9]*self.CHBCFINOSS[7] )  ) #if nOxy< 3*nSi, the int() should bring it down to nSi
        nSilox_C1=nSilox/(nCarb+1)
        nSilane=self.CHBCFINOSS[9]-nSilox
        if(nSilane<0):nSilane=0
        nIodo=self.CHBCFINOSS[5]
        nFluoro=(self.CHBCFINOSS[4])
        nChloro=(self.CHBCFINOSS[3])
        nNitro=self.CHBCFINOSS[6]
        nOxy=self.CHBCFINOSS[7]
        nOx_C1=nOxy/(nCarb+1)
        nF_C1=nFluoro/(nCarb+.1)
        #nHy=(self.CHBCFINOSS[1])
        #nSulf=(self.CHBCFINOSS[8])
        #nHalo=self.CHBCFINOSS[2]+self.CHBCFINOSS[3]+self.CHBCFINOSS[4]+self.CHBCFINOSS[5] #BrClFI
        #nHetero_C1=(nHalo+nSilox+nNitro+nOxy)/(nCarb+.1) #Sulfur ~carb, but keep separate
        #nHalo_C1=(nHalo)/(nCarb+.1) 
        #H_C1=nHy/(nCarb+nSilane+.1) #Assume silanes behave like their HC hmomorph.
        SNOCfactor=0
        if(nChloro>.95 or nSilox_C1>.95 or nNitro>.95 or nOx_C1 > 0.45 or nIodo>.95):SNOCfactor=1
        if(nChloro>1.95 or nSilox_C1>1.95 or nNitro>1.95 or nOx_C1 > 0.85):SNOCfactor=2
        if(nSilane>0):SNOCfactor=0
        if(nF_C1>0.4):SNOCfactor=0
        if(compClass=='halom'):SNOCfactor=0
        alpha=1
        q = self.qShape
        eps=self.eps_kB*kB/1E21
        sig3=self.bVolCC_mol/q*6/np.pi/(avoNum) # in nm^3
        sigma=(sig3**(1./3))*1E-9

        ciiCorr = avoNum*1E21* q * (0.12008072630855947 + 2.2197907527439655 * alpha)
        ciiCorr *= ciiCorr
        ciiCorr *= (1.4*eps * sigma**5) #JRE: preliminary calcs showed that correcting the vrMie correlation by 1.4 provided a rough match for methane at Tbp.

# =============================================================================
#         cii = self.ms * (0.12008072630855947 + 2.2197907527439655 * self.alpha)
#         cii *= Na*np.sqrt(self.eps * self.sigma**5)
#         cii **= 2
#         cii = np.array([cii], ndmin=1)
# =============================================================================
        
        cii0=np.sqrt(ciiCorr) #JRE: residual correlation summed to 0 at Tr=0, so use this element to store reference value.
        cii1= -21.521 #JRE: start with norml class as default
        cii2= 112.771
        if(compClass=='heavy'):
            cii1=16.904
            cii2= 2.533
        if(compClass=='assoc' or compClass=='Asso+'):
            cii1=113.595
            cii2=-60.637
        elif(compClass=='polar'):
            cii1=19.134
            cii2=25.046
        cii3= -(cii1+cii2+37) #JRE: for some reason all values sum to 37 at Tr=1
        #cii1=cii2=cii3=0.        
        cii = np.array([cii0,cii1,cii2,cii3] )*(1+0.0412*SNOCfactor) #, ndmin=4)
        #cii *= ciiCorr
        if overwrite:
            self.cii = cii
        return cii
    def SetEsdParms(self,qShape,eps_kB,bVol):
        self.qShape=qShape
        self.eps_kB=eps_kB
        self.bVolCC_mol=bVol
        self.cShape=1+ (qShape-1)*(4-1.9)/4
        if self.bAssoc==True : 
            self.assoc=Assoc()
            self.SetAssocParms() #(bondVolNm3,epsD_kB,epsA_kB,nTypes,nDsites,nAsites)
    def SetAssocParms(self): #bondVolNm3,epsD_kB,epsA_kB,nTypes,nDsites,nAsites):
        self.assoc.bondVolNm3=0.001 # These are default values for alcohols.  
        self.assoc.epsD_kB=2000     # Proper loading would require file reading.
        self.assoc.epsA_kB=2000     # Convention to use lower "assoc" for instance name but upper "Assoc" for Class name.
        self.assoc.nDegree=1
        self.assoc.nTypes=1
        self.assoc.nDsites=[0]      # In principle, the number of donor sites is a list over nTypes. For ESD, nTypes=1.
        self.assoc.nAsites=[0]
    def LookupEsdParms(self,idDippr):
        # EsdParms={} # a dictionary 
        inFile=PGLInputDir+"ParmsEsdMEM2.txt"
        #inFile="ParmsEsdMEM2.txt"
        with open(inFile,encoding='utf-8') as f:
            lines=f.readlines() # reads the entire file as lines.
            # q EOK(K)   b(cc/mo)   KadNm3  eDonEpsK    eAccEpsK    nDegree #eDs    #eAs    idCas   Name
            bThere=False
            for line in lines:
                part=line.split() # separates out the  
                if int(part[0]) == idDippr:
                    bThere=True
                    print(f"Compound {self.Name} found in ParmsEsdMEM2.txt!")
                    self.qShape=float(part[1]) # All parts are read as strings. Must convert manually after split().
                    self.eps_kB=float(part[2])
                    self.bVolCC_mol=float(part[3])
                    self.assoc.bondVolNm3[0]=float(part[4])  
                    self.assoc.epsD_kB[0]=float(part[5])     # Proper loading requires file reading.
                    self.assoc.epsA_kB[0]=float(part[6])
                    self.assoc.nDegree[0]=int(part[7])# Convention to use lower "assoc" for instance name but upper "Assoc" for Class name.
                    self.assoc.nDsites[0]=int(part[8])
                    self.assoc.nAsites[0]=int(part[9])
                    break # exit the for-loop search if found.
        return bThere
    def EsdDenLiqG_cc(self,tKelvin): 
        """ Estimate liquid density according to quadratic solution of the ESD EOS when Z~0.
            Ref. MEM1(1996),EsdPoly(2002) """
        # Z = 1 + q[4η/(1-1.9η) - 9.5Yη/(1+k1Yη)] - (q-1)1.9η/(1-1.9η) - Fassoc^2/(1-1.9η)
        # Z = 1 + 4cη/(1-1.9η) - 9.5qYη/(1+k1Yη) - Fassoc^2/(1-1.9η)
        # 0 = (1-1.9η)*(1+k1Yη) + (4cη-F2)*(1+k1Yη) -9.5qYη(1-1.9η)
        # 0 = 1+k1Yη-1.9η-1.9*k1Y*η^2 + 4cη-F2+4c*k1Yη^2-F2k1Yη - 9.5qYη + 9.5*1.9qYη^2 # Treat F2 as ~constant when alpha is "large" or ~zero.
        # 0 = 1+η(k1Y-1.9+4c-F2k1Y-9.5qY)+η^2*(-1.9k1Y+4c*k1Y+9.5*1.9*qY) ≡ c + bη + aη^2
        # η = { -(-9.5qY+k1Y+4c-F2k1Y) + [ (9.5qY-k1Y-4c+F2k1Y)^2 -4Y*(-1.9k1+4c*k1+9.5*1.9*q) ]½ }/( 2*Y(-1.9k1+4c*k1+9.5*1.9*q) )
        if tKelvin/self.TcK > 0.99: # This method only works for fluids that are well into the liquid range. 
            return 86
        bVol=self.bVolCC_mol
        eta=self.liqDenG_cc(tKelvin)*bVol/self.Mw # start with Rackett to get Fassoc estimate.
        Psat=self.PvpSc(tKelvin)
        Zsat=Psat*bVol/(Rgas*eta*tKelvin)
        if Zsat > 0.35: return 8686 # Something must be wrong if ZLiq > ZcEsd
        Yhb=( np.exp(self.assoc.epsA_kB[0]/tKelvin)-1 )*( np.exp(self.assoc.epsD_kB[0]/tKelvin)-1 )
        sqArg=( eta/bVol*self.assoc.bondVolNm3[0]*Yhb )
        ralph=math.sqrt( sqArg )
        Nd=self.assoc.nDegree[0]
        F2=(  2*Nd*ralph/( 1+math.sqrt(1+4*Nd*ralph*ralph) )  )**2 # EsdPoly below Eq.1
        # if self.bAssoc : F2=1
        cShape=( 4*self.qShape-(self.qShape-1)*1.9 )/4
        # 4c-4 = 2.1(q-1)+4-4 => q = (4c-4)/2.1+1 = (4c-1.9)/(4-1.9)
        bigQ=(4*cShape*eta-F2)/self.qShape
        bigY=math.exp(self.eps_kB/tKelvin)-1.0617
        k1=1.7745
        aQuad= bigY*( 4*cShape*k1+self.qShape*1.9*9.5-1.9*k1 )
        #            ( 4*cShape*k1-1.9*k1+1.9*9.5*self.qShape)*bigY
        bQuad= k1*bigY-1.9+4*cShape-9.5*self.qShape*bigY-F2*k1*bigY
        cQuad= 1-Zsat
        sqArg= bQuad*bQuad-4*aQuad*cQuad
        if sqArg < zeroTol:
            print(f"F^2={F2}, qShape={self.qShape}, bigQ={bigQ}, bigY={bigY}. ") 
            print(f"aQuad={aQuad}, bQuad={bQuad}, sqArg={sqArg}")
            """ For debugging.""" 
            for iEta in range(25,25):
                eta=iEta/100 # Z = [ 1+ η*(k1Y+Qq-9.5qY) + (Qk1+1.9*9.5 )qYη2 ]/(1-1.9eta)/(1+k1Yeta) 
                Zrep=  4*cShape*eta/(1-1.9*eta)
                Zatt= -9.5*self.qShape*bigY*eta/(1+k1*bigY*eta)
                Zassoc = -F2/(1-1.9*eta)
                Zcalc= 1+Zrep+Zatt+Zassoc
                #ZrepAssoc=Zrep+Zassoc
                bigQ=(4*cShape*eta-F2)/self.qShape
                #ZrepAssoc=(4*cShape*eta-F2)/(1-1.9*eta)
                Zcalc=1+(4*cShape*eta-F2)/(1-1.9*eta)-9.5*self.qShape*bigY*eta/(1+k1*bigY*eta)
                Zero=Zcalc*(1-1.9*eta)*(1+k1*bigY*eta)
                """
                ZeroRepAssoc=(1+k1*bigY*eta)*(1-1.9*eta)+(4*cShape*eta-F2)*(1+k1*bigY*eta)
                ZrepAssoc=ZeroRepAssoc/( (1+k1*bigY*eta)*(1-1.9*eta) ) -1 
                ZeroRepAssoc=1+k1*bigY*eta-1.9*eta-1.9*k1*bigY*eta**2+(4*cShape*eta-F2)+4*cShape*k1*bigY*eta**2-F2*k1*bigY*eta
                ZrepAssoc=ZeroRepAssoc/( (1+k1*bigY*eta)*(1-1.9*eta) ) -1 
                ZeroRepAssoc=1-F2+(k1*bigY-1.9+4*cShape-F2*k1*bigY)*eta-1.9*k1*bigY*eta**2+4*cShape*k1*bigY*eta**2
                ZrepAssoc=ZeroRepAssoc/( (1+k1*bigY*eta)*(1-1.9*eta) ) -1 
                ZeroRepAssoc=1-F2+(k1*bigY-1.9+4*cShape-F2*k1*bigY)*eta +(4*cShape*k1-1.9*k1)*bigY*eta*eta
                ZrepAssoc=ZeroRepAssoc/( (1+k1*bigY*eta)*(1-1.9*eta) ) -1 
                ZeroAtt= -9.5*self.qShape*bigY*eta+1.9*9.5*self.qShape*bigY*eta*eta
                """
                Zero=1-F2+(k1*bigY-1.9+4*cShape-F2*k1*bigY-9.5*self.qShape*bigY)*eta +(4*cShape*k1-1.9*k1+1.9*9.5*self.qShape)*bigY*eta*eta
                Zcalc=Zero/(1-1.9*eta)/(1+k1*bigY*eta)
                Zero=(cQuad+eta*bQuad+aQuad*eta*eta) #/(1-1.9*eta)/(1+k1*bigY*eta)
                Zcalc=Zero/(1-1.9*eta)/(1+k1*bigY*eta)
                print(f"eta={eta},Zcalc={Zcalc}, Zrep={Zrep},Zatt={Zatt},Zassoc={Zassoc} ")
            # η = (9.5Y-Q)/(QY+9.5*1.9*bigY); Y > 1 # result from q->inf
            etaLiq=(9.5*bigY-bigQ)/(bigQ*bigY+9.5*1.9*bigY)
            print(f"etaLiqInfChain={etaLiq}")
            denLiqG_cc=etaLiq*self.Mw/bVol # result from q->inf
            #return 868686
        etaLiq= ( -bQuad+math.sqrt(sqArg) )/(2*aQuad)
        print(f"etaLiqQuad={etaLiq}")
        denLiqG_cc=etaLiq*self.Mw/bVol
        return denLiqG_cc
    ##########################################################################################################
    def ChemPoTV(self,bZiter,tKelvin,vCC_mol): 
        """ returns iErr, PMPa,zFactor,aRes,uRes,chemPo """
        iErr=0
        PMPa=zFactor=aRes=uRes=chemPo=0.
        cqFactor=4/(4-1.9) # = 1.90476... where c-1 = (q-1)*cqFactor.
        K1=1.7745
        K2=1.0617
        rho=1/vCC_mol # passing V instead of rho for consistency with mixture code
        bVol=self.bVolCC_mol
        eta=rho*bVol
        voidFrac=1-1.90*eta
        rdfContact=1/voidFrac
        if(eta < zeroTol or voidFrac < zeroTol):
            print(f"ChemPoTV: eta={eta} voidFrac={voidFrac}. Returning iErr=11.")
            iErr=11             
            return iErr, PMPa,zFactor,aRes,uRes,chemPo    
        qShape=self.qShape
        cShape=( 1+(qShape-1)/cqFactor ) # q-1= 1.90476*(c-1) => c = 1+(q-1)/1.90476
        bigY=np.exp(self.eps_kB/tKelvin)-K2
        #K1Yb=K1*bigY*bVol
        nDonTot=self.assoc.nDegree[0]*self.assoc.nDsites[0] #this a vector for single comp in ESD cuz assoc class is general.
        nAccTot=self.assoc.nDegree[0]*self.assoc.nAsites[0]
        bepsA=(self.assoc.epsA_kB[0]/tKelvin)
        bepsD=(self.assoc.epsD_kB[0]/tKelvin)
        YhbA=( np.exp(bepsA)-1 )
        YhbD=( np.exp(bepsD)-1 )
        sqArg=rho*rdfContact*avoNum*self.assoc.bondVolNm3[0]*YhbA
        if(sqArg < 0): #zeroTol is too large for this because HCs=0.00000
            iErr=12
            print(f"ChemPoTV: ralphA calc - rho,sqArg={rho,sqArg}")
            return iErr, PMPa,zFactor,aRes,uRes,chemPo    
        ralphA=math.sqrt( sqArg )
        ralphD=math.sqrt( rho*rdfContact*avoNum*self.assoc.bondVolNm3[0]*YhbD )
        alphAD=ralphA*ralphD
        quadB=1+alphAD*(nAccTot-nDonTot) # Maybe rare, but possible nAcc!=nDon
        FA=(  2*nDonTot*ralphD/( quadB+math.sqrt(quadB*quadB+4*nDonTot*alphAD) )  ) # EsdPoly below Eq.1
        FD=nAccTot*ralphA/(1+FA*ralphA)
        zAssoc= -FA*FD/(1-1.9*eta)
        aAssoc = nAccTot*( -np.log(1+FA*ralphA)+0.5*FA*ralphA/(1+FA*ralphA)) # XA=1/(1+FA*ralphA)
        aAssoc+= nDonTot*( -np.log(1+FD*ralphD)+0.5*FD*ralphD/(1+FD*ralphD)) # XD=1/(1+FD*ralphD)
        dLnAlpha_dLnBeta=np.exp(bepsA)	# Eqs. 44,45
        if(bepsA > 1.E-4):dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsA)/(dLnAlpha_dLnBeta-1)
        betadFD_dBeta= nAccTot*ralphA*0.50*dLnAlpha_dLnBeta/(1+FA*ralphA)
        dLnAlpha_dLnBeta=np.exp(bepsD)
        if(bepsD > 1.E-4):dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsD)/(dLnAlpha_dLnBeta-1)
        betadFA_dBeta = nDonTot*ralphD*0.5*dLnAlpha_dLnBeta/(1+FD*ralphD) 
        uAssoc= -( FA*betadFD_dBeta+FD*betadFA_dBeta ) # Eqs. 43,44
        voidFrac=1-1.90*eta
        denom=voidFrac
        zRep= 4*cShape*eta/denom
        zAtt= -9.5*bigY*qShape*eta/(1+K1*bigY*eta)
        zFactor=(1+zRep+zAtt+zAssoc)
        PMPa=zFactor*Rgas*rho*tKelvin
        aRep= -4/1.9*np.log(voidFrac)*cShape
        aAtt= -9.5*qShape*np.log(1+K1*bigY*eta)/K1
        aRes=aRep+aAtt+aAssoc #-DLOG(Z) !don't subtract log(z) for aRes(T,V). Important for EAR.
        uAtt= -9.5*qShape*(bigY+1)/(1+K1*bigY*eta)
        uRes=uAtt+uAssoc
        chemPo=aRes+zFactor-1
        if( zFactor < zeroTol):  # Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
            iErr=3	  # warning level because another call might produce Z > 0.
        #if(bZiter | iErr>0):return iErr, PMPa,zFactor,aRes,uRes  # don't need the rest if bZiter.
        return iErr,PMPa,zFactor,aRes,uRes,chemPo
    ##########################   End of ChemPoTV(self,bZiter,tKelvin,rhoMol_cc):   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def ChemPoTP(self,LIQ,tKelvin,PMPa): 
        """Returns iErr, rhoMol_cc,zFactor,uRes,aRes,chemPo """
        # LIQ = 0,2 for vapor root, 1,3 for liquid root. 2,3 indicate zIter calculation (no uRes etc)
        bLiq=False
        if(LIQ==1 or LIQ==3):bLiq=True
        iErr=0
        rhoMol_cc=0.
        zFactor=1.
        uRes=0.
        aRes=0.
        chemPo=0
        #NOTE: iErrTmin is checked in FuVtot
        eta=0
        if(bLiq):eta=self.etaMax/1.15
        bVol=self.bVolCC_mol
        bMix=bVol
        Pb_RT=PMPa*bMix/(Rgas*tKelvin)
        #GUESS FOR rho
        eta=Pb_RT  #NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
        if(bLiq or eta>self.etaMax):eta=self.etaMax/1.15
        rho=eta/bMix
        #if(eta > 1/1.9 .and. LOUDER)write(dumpUnit,*)'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
        if(eta < 0):
            print(f"ChemPoTP: initial eta={eta}< 0. Returning iErr=11.")
            iErr=11
            return iErr, rhoMol_cc,zFactor,aRes,uRes,chemPo
        bZiter=True # chemPo calculations are skipped for isZiter=1
        iErr, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,tKelvin,1/rho)
        #Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
        if(iErr > 10):
            return iErr, rhoMol_cc,zFactor,aRes,uRes,chemPo
        etaOld=eta
        errOld=Pb_RT-eta*zFactor
        eta=etaOld/1.150
        #IF (eta < 0 .and. LOUD) write(dumpUnit,31)LIQ
        #if(initial==1.and.LOUD)write(dumpUnit,*)'FugiEsd: initial eta,err',etaOld,errOld
        itMax=77
        errBesteta=1234
        for nIter in range(itMax):
            rho=eta/bMix
            iErr, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,tKelvin,1/rho)
            if(iErr > 10):break
            ERR=Pb_RT-eta*zFactor
            CHNG=ERR/(ERR-errOld)*(eta-etaOld)
            #if(initial==1.and.LOUDER)write(dumpUnit,'(a,2e11.4,3f10.5)')'FugiEsd eta,Z', eta,zFactor
            #if(initial==1.and.LOUDER)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd eta,CHNG,niter',eta,CHNG,niter 
            etaOld=eta
            errOld=ERR
            #  LIMIT THE CHANGE IN Density for liquid.
            if(bLiq and np.abs(CHNG/etaOld) > 0.15):CHNG=np.sign(CHNG)*0.15*etaOld
            #  Low eta must move from zero, so you can't limit its % change
            if(not bLiq and np.abs(CHNG) > 0.020):CHNG=np.sign(CHNG)*0.02
            eta=eta-CHNG
            if(np.abs(CHNG) < errBesteta):
                etaBest=eta
                errBesteta=np.abs(CHNG)
            if(eta < 0 ):
                print(f'ChemPoTP: LIQ={LIQ}, etaV={eta}. Shrinking change.')
                eta=etaOld-CHNG*0.001
            if(eta < 0):
                print(f"ChemPoTP: iterated eta={eta}< 0. Returning iErr=17.")
                iErr=17
                return iErr, rhoMol_cc,zFactor,aRes,uRes,chemPo
            if(eta > self.etaMax):
                print(f'ChemPoTP: LIQ={LIQ}, etaL={eta}. Shrinking change.')
                eta=etaOld-np.sign(CHNG)*0.01*(self.etaMax-etaOld)
            if(np.abs(CHNG) < 1.E-9 and eta > 0):break  # Don't use CHNG/eta. Converge quickly to ideal gas if P->0, ~9 sigfigs if liquid
        ##########################   Iteration Concluded    !!!!!!!!!!!!!!!!!!!!!!!!!
        if(eta < 0 or eta > self.etaMax):
            iErr=17
        rho=eta/bMix
        rhoMol_cc=rho 
        #!  ITERATION ON rho HAS CONCLUDED.  
        #zStore=zFactor   
        bZiter=False
        #One last call to get uRes,chemPo.
        iErrF, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,tKelvin,1/rho)
        if(zFactor < 0):iErrF=11 #JRE: Note that zFactor=0 when computing P=0 root.
        if(iErrF > 0 or nIter > itMax-1 or eta < 0): # if iErr still > 0 on last iteration, then declare error.
            iErr=iErrF 
            if(nIter > itMax-1):iErr+=100
            eta=etaBest
        if(zFactor > zeroTol):chemPo-=np.log(zFactor)	 #!Must subtract ln(Z) when given P as independent variable.
        return iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo #,eta
    ######################   End of Comp.ChemPoTP(LIQ,tKelvin,PMPa)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ##########################   Start SgtPy funs   ############################# 
    def psat(self, T, P0=None, v0=[None, None], Xass0=[None, None],
             full_output=False):
        """
        psat(T, P0)

        Method that computes saturation pressure at fixed T

        Parameters
        ----------

        T : float
            absolute temperature [K]
        P0 : float, optional
            initial value to find saturation pressure [Pa]
        v0: list, optional
            initial guess for liquid and vapor phase, respectively [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
        full_output: bool, optional
            whether to outputs or not all the calculation info.

        Returns
        -------
        psat : float
            saturation pressure [Pa]
        vl : float
            liquid saturation volume [m3/mol]
        vv : float
            vapor saturation volume [m3/mol]
        iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV=self.PvpEar(T)
        vLiq=1/(rhoLiq*1E6)
        vVap=1/(rhoVap*1E6)
        out=PMPa*1E6,vLiq,vVap
        """
        out = psat(self, T, P0, v0, Xass0, full_output)
        return out
    def PvpEar(self,tKelvin):
        """ Returns iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
	    !COMPUTE VAPOR PRESSURE GIVEN tKelvin using the equal area rule.  
	    !Ref: Eubank et al, IECR 31:942 (1992)"""
        errMsg=["0"]*22 #!,errMsgPas
        errMsg[ 1]='PvpEar: Warning P < 1E-8 & iter > 2. Probably T < Tmin.'
        errMsg[ 2]='PvpEar: Warning dP/P > tol. Probably T < Tmin.'
        errMsg[ 3]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[ 4]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[ 9]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[11]='PvpEar: No spinodal max/min'
        errMsg[12]='PvpEar: Psat iteration did not converge'
        errMsg[13]='PvpEar: Liquid FUGI call failed on last iteration'
        errMsg[14]='PvpEar: Vapor  FUGI call failed on last iteration'
        errMsg[15]='PvpEar: zVap=zLiq on last iteration'
        errMsg[16]='PvpEar: ChemPoTP returned T < Tmin error'
        errMsg[17]='PvpEar: zLiq < 0'
        errMsg[18]='PvpEar: Tr > 1-zeroTol'
        errMsg[19]='PvpEar: Critical Error from fugacity calculation'
        errMsg[20]='PvpEar: rhoVap/rhoLiq < 0'
        #tMin=component.TcK*0.25
        #NC=1 # Ltd to single compound.
        iErr=0 
        PMPa=chemPot=rhoLiq=rhoVap=uSatL=uSatV=0.
        zCritEff=self.Zc
        zCritEff=self.ZcEos  # For ESD, this is close enough
        rhoCrit=self.PcMPa/(zCritEff*Rgas*self.TcK)
        rhoVap=rhoCrit*1.000 #-ve value on input means use default value for etaHi
        Tr=tKelvin/self.TcEos
        if(Tr > 1-zeroTol):
            iErr=18
            return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        elif(Tr > 0.85):
            # etaL=2A(1-Tr)-etaV; 2A(1-Tr)- 2etaV= B*(1-Tr)^.5; etaV = 0.1672*(1-Tr)-1.22(1-Tr)^.5/2; Fit for EsdCO2.
            etaVap=self.etacEos+0.1672*(1-Tr)-1.22/2*np.sqrt(1-Tr)
            etaLiq=(2*.1672*(1-Tr)-etaVap)
            rhoVap=etaVap/self.bVolCC_mol
            rhoLiq=etaLiq/self.bVolCC_mol
            pOld=self.PcMPa*10**( 7/3*(1+self.acen)*(1-1/Tr) ) # SCVP
            #pEubank=Rgas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) ) #EAR method of Eubank and Hall (1995)
            #pOld=pEubank/10  #JRE: I find Eubank overestimates P sometimes. It even gets higher than pMax
        else: # If Tr < 0.85, compute pSatItic estimate  FPE,501:112236(19), JCP,110:3043(99)
            pOld=0.01 #If p < 0, then the vapor density goes negative and log(rhoLiq/rhoVap) is indeterminate.
            #NOTE: Don't use pOld=0 when pOld > 0. You might get an error because vdw loop doesn't cross zero#
            iErrF,rhoLiq,zLiq,aResLiq,uSatL,chemPo=self.ChemPoTP(1,tKelvin,pOld)
            if(iErrF > 10):
                iErr=19
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            elif(iErrF > 0):
                iErr=9
            elif(zLiq < zeroTol):
                iErr=17
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
                #AresLiq=Ares_RT						#AresVap  +  Zvap-1 -ln(Zvap) =AresLiq+ZLiq-1-ln(ZLiq)							 
                #rhoLiq=pOld/(zLiq*Rgas*tK)		#=>	#B2*rhoVap+B2*rhoVap+ln(rhoVap/rhoLiq) =AresLiq+0-1  
            rhoVap=rhoLiq*np.exp( aResLiq-1 ) #pSatItic, FPE, 501:112236(19), Eq 19. Ares_RT from GlobConst
            pSatItic=rhoVap*Rgas*tKelvin
            etaVap=rhoVap*self.bVolCC_mol
            #call FuVtot(1,tK,pSatItic,xFrac,NC,0,FUGC,zVap,ier)
            #call FuVtot(1,tK,1/rhoVap,xFrac,NC,FUGC,zVap,aRes,uRes,iErrF) #isZiter=1=>vapor Z with no fugc calculation. zVap >> zero.
            bZiter=False
            iErrF,Pcalc,zVap,aResVap,uSatV,chemPot=self.ChemPoTV(bZiter,tKelvin,1/rhoVap)
            if(iErrF > 10):
                iErr=19
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            elif(iErrF > 0):
                iErr=9
            elif(zLiq < zeroTol):
                iErr=17
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            #B2cc_mol*rhoVap=(zVap-1)
            rhoVap=rhoLiq*np.exp( aResLiq-1-2*(zVap-1) )# => rhoVap=rhoLiq*EXP(AresLiq-1-2*B2*rhoVap).
            pSatItic=rhoVap*Rgas*tKelvin*zVap
            if(pSatItic > zeroTol):
                pTest=pSatItic
                zLiq=pSatItic/(rhoLiq*Rgas*tKelvin) # if Psat > 0 then zLiq > 0
                #FugcLiq=aResLiq+zLiq-1-np.log(zLiq)
            else:
                pTest=zeroTol*10
            #pEuba2=Rgas*tK/(1/rhoVap-1/rhoLiq)*( Ares_RT- 0 +DLOG(rhoLiq/rhoVap) )	#EAR method of Eubank and Hall (1995)
            #NOTE: When 1/rhoVap >> 1/rhoLiq, pEuba2=Rgas*tK*rhoVap*( Ares_RT-ln(rhoVap/rhoLiq) )=pTest*( Ares_RT - (Ares_RT-1) )
            pOld=pTest
        # end if:Tr > 0.85...
        #etaLiq=rhoLiq*self.bVolCC_mol
        if(pOld < zeroTol):pOld=zeroTol
        pBest=1234
        fBest=1234
        itMax=33
        for iter in range(itMax): #iterate on pOld according to Eubank criterion
            iErrF,rhoLiq,zLiq,aResLiq,uSatL,chemPot=self.ChemPoTP(3,tKelvin,pOld)
            if(iErrF > 0):iErr=3 #declare error but don't stop. if future iterations give valid fugi(), then no worries
            iErrF,rhoVap,zVap,aResVap,uSatV,chemPot=self.ChemPoTP(2,tKelvin,pOld)
            if(iErrF > 0):iErr=4 #declare error but don't stop. if future iterations give valid fugi(), then no worries
            if( (rhoVap/rhoLiq) < 0 ):
                iErr=20
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            if(np.abs(zVap-zLiq) < 1e-3):
                iErr=15
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            if(pOld < 1E-8 and iter > 2): # Just use pSatItic estimate in this case. Declare warning.
                #etaLiq=rhoLiq*component.bVolCC_mol
                rhoVap=rhoLiq*np.exp( aResLiq+zLiq-1-aResVap+zVap-1 ) #pSatItic, FPE, 501:112236(19), Eq 19.aDepVap=/=0 for carbo acids
                pTest=rhoVap*Rgas*tKelvin*zVap
                iErr=1 # assign warning in this case. 
                break
            pTest=pOld*( aResLiq-aResVap-np.log(rhoVap/rhoLiq) )/(zVap-zLiq)	#EAR method of Eubank and Hall (1995)
            fErr= aResLiq+zLiq-1 -(aResVap+zVap-1) -np.log(rhoVap/rhoLiq) 
            if( np.abs( fErr) < fBest):	# for very low pSat, like propane.
                fBest=np.abs( fErr)
                pBest=pTest
            change=pTest-pOld
            tol=1E-7
            if(pTest > 0):
                pOld=pTest
            else:
                pOld=pOld/2
            if(np.abs(fErr) < tol):break   
            # end for itMax
        PMPa=pOld
        if(np.abs(change/PMPa) > tol or iter >= itMax):
            iErr=2
            PMPa=pBest
            return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        #iter=0                                  
        #check that Fugi did not give warning on last call
        if(iErr > 10):return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        #if(ierCode==0) call FuVtot(isZiter,tK,1/rhoLiq,xFrac,NC,FUGC,ZLiq,iErrFu) #one last call to fugi for chemPo and debugging.
        chemPot=86.86
        if( zLiq > 0):
            chemPot=aResLiq+zLiq-1 -np.log(zLiq)
        else:
            iErr=20
        return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
##########################   END of PvpEar()   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    def temperature_aux(self, T):
        """
        temperature_aux(T)

        Method that computes temperature dependent parameters.
        It returns the following list:

        temp_aux = [beta, dia, tetha, x0, x03, Fab, epsa]

        Journal of Chemical Physics, 139(15), 1–37 (2013)

        beta: Boltzmann's factor [1/J]
        dia: computed diameter [m] (Eq 7)
        tetha: exp(beta*eps)-1 [Adim] (Below Eq. 63)
        x0: sigma/dia [Adim] (Below Eq. 17)
        x03: x0^3 [Adim]
        Fab: association strength [Adim] (Below Eq. 77)
        epsa: eps / kb / T [Adim]

        Parameters
        ----------

        T : float
            Absolute temperature [K]

        Returns
        -------
        temp_aux : list
             list of computed parameters
        """

        beta = 1 / (T*kB/1E21) #interface with sgtpy requires kB in SI.

        temp_aux = [beta]
        return temp_aux
    
    def logfug_aux(self, temp_aux, P, state, v0=None, Xass0=None):
        """
        logfug_aux(T, P, state, v0, Xass0)

        Method that computes the fugacity coefficient at given
        composition, temperature and pressure.

        Parameters
        ----------
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        logfug : float
            fugacity coefficient
        v : float
            computed volume of the phase [m^3/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        #rhomolecular = avoNum * rho
        beta = temp_aux[0]
        #beta = 1 / (T*kB/1E21) #interface with sgtpy requires kB in SI.
        tKelvin=1/(beta*kB/1E21)
        PMPa=P/1E6
        LIQ=0
        if(state=='L'):LIQ=1
        #ar, Xass = ares(self, rhomolecular, temp_aux, Xass)
        iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo=self.ChemPoTP(LIQ,tKelvin,PMPa)
        #RT = Na/beta
        #Z = P * v / RT
        if(zFactor < zeroTol or rhoMol_cc <zeroTol or iErr>10):
            vM3_mol=8686
            lnphi=8686
        else:
            lnphi = aRes + (zFactor - 1.) - np.log(zFactor)
            vM3_mol=1E-6/rhoMol_cc
        return lnphi, vM3_mol, Xass0

    def muad_aux(self, rho, temp_aux, Xass0=None):
        """
        muad_aux(rho, temp_aux, Xass0)

        Method that computes the adimenstional chemical potential at given
        density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        muad: float
            chemical potential [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """

        #beta = 1 / (T*kB/1E21) #interface with sgtpy requires kB in SI.
        #tKelvin=1/(beta*kB/1E21)
        rhomolecular = rho * avoNum*1E21
        da, Xass = self.dafcn_aux(rhomolecular, temp_aux, Xass0)
        afcn, dafcn = da
        mu = afcn + rhomolecular * dafcn

        return mu, Xass
    def dafcn_aux(self, rho, temp_aux, Xass0=None):
        """
        dafcn_aux(rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its first density derivative.

        Parameters
        ----------
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array
           Helmholtz free energy and its derivative  [J/mol, J m^3/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        #beta = 1 / (T*kB/1E21) #interface with sgtpy requires kB in SI.
        tKelvin=1/(beta*kB/1E21)
        #rhomolecular = rho * avoNum*1E21
        #tKelvin=1/(beta*kB)
        a, Xass = self.dares_drho(rho, tKelvin, Xass0)
        #da = (Z-1)/rhoMolecular
        daIg= daideal_drho(rho, beta)
        a += daIg 
        a *= (avoNum*1E21/beta) #daIg = 1/rhoMolecular
        # FYI:(avoNum*1E21/beta)=Rgas*tKelvin
        # on return da= RT*(Z-1)/(rhoMol_m3*beta) => RT in numerator????
        #a = a/(Rgas*tKelvin)
        return a, Xass0
    def d2afcn_aux(self, rho, temp_aux, Xass0=None):
        """
        dafcn_aux(rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its second density derivative."""
        delRho=rho*0.001
        rhoPlus =rho+delRho
        rhoMinus=rho-delRho
        aPlus, XassPlus = self.dafcn_aux(rhoPlus , temp_aux, Xass0)
        aMinus,XassMinus= self.dafcn_aux(rhoMinus, temp_aux, Xass0)
        aResPlus ,daResPlus =aPlus
        aResMinus,daResMinus=aMinus
        aRes=0.5*(aResPlus+aResMinus)
        daRes=0.5*(daResPlus+daResMinus)
        d2aRes=(daResPlus-daResMinus)/(rhoPlus-rhoMinus)
        a=[aRes,daRes,d2aRes]
        # FYI:(avoNum*1E21/beta)=Rgas*tKelvin
        # on return da= RT*(Z-1)/(rhoMol_m3*beta) => RT in numerator????
        #a = a/(Rgas*tKelvin)
        return a, Xass0
    def dares_drho(self, rho, T, Xass0=None):
        """
        returns         return a[aRes_RT,(Z-1)/rhoMolecular], Xass0
        dares_drho(rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the fluid
        and its first density derivative.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimensionless Helmholtz free energy [Adim, m^3]
        """
        #temp_aux = self.temperature_aux(T)
        #a, Xass = dares_drho(self, rho, temp_aux, Xass0)
        #gmol=np.zeros(2) # force 2D to avoid python declaring 1D array as scalar.
        #gmol[0]=1
        bZiter=False
        rhoMol_cc=rho/(avoNum*1E27) # includes the conversion from m3 to cc.
        iErrF,PMPa,zFactor,aRes_RT,uRes_RT,chemPot=self.ChemPoTV(bZiter,T,1/rhoMol_cc)
        if(iErrF < 11):
            dares_drho=(zFactor-1)/rho #*Rgas*T # Note: rho=rhoMolecular here.
            a=[aRes_RT,dares_drho]
        else:
            a=[8686,8686]
        #a*=(Rgas*T)  # dafcn requires a in Joules...
        return a, Xass0
    def ci(self, T):
        '''
        ci(T)
        Method that evaluates the polynomial for the influence parameters used
        in the SGT theory for surface tension calculations.
        Parameters
        ----------
        T : float
            absolute temperature [K]
        Returns
        -------
        ci: float
            influence parameters [J m5 mol-2]
        '''
        Tr=T/self.TcK
        #ciiCorr=np.polyval(self.cii, Tr)
        cii=self.cii
        pctDevOrig=Tr*( cii[1]+Tr*(cii[2]+Tr*cii[3]) )
        ciiCorr=cii[0]/(  1+pctDevOrig/100  ) #**2
        return ciiCorr*ciiCorr #JRE: here we square it after correcting. 
    def sgt_adim(self, T):
        '''
        sgt_adim(T)

        Method that evaluates adimensional factor for temperature, pressure,
        density, tension and distance for interfacial properties computations
        with SGT.

        Parameters
        ----------
        T : float
        absolute temperature [K]

        Returns
        -------
        Tfactor : float
            factor to obtain dimensionless temperature (K -> K)
        Pfactor : float
            factor to obtain dimensionless pressure (Pa -> Pa/RT)
        rofactor : float
            factor to obtain dimensionless density (mol/m3 -> mol/m3)
        tenfactor : float
            factor to obtain dimensionless surface tension (mN/m)
        zfactor : float
            factor to obtain dimensionless distance  (Angstrom -> m)
        '''
        cii = self.ci(T)  # computing temperature dependent cii

        Tfactor = 1.
        Pfactor = 1.
        rofactor = 1.
        tenfactor = np.sqrt(cii) * 1000  # To give tension in mN/m
        zfactor = 10**-10

        return Tfactor, Pfactor, rofactor, tenfactor, zfactor
    def a0ad_aux(self, rho, temp_aux, Xass0=None):
        """
        a0ad_aux(ro, temp_aux, Xass0)

        Method that computes the adimensional Helmholtz density energy at
        given density and temperature.

        Parameters
        ----------

        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a0ad: float
            Helmholtz density energy [J/m^3]
        Xass : array
            computed fraction of nonbonded sites
        """
        rhomolecular = rho * avoNum*1E21
        a, Xass = self.dafcn_aux(rhomolecular, temp_aux, Xass0)
        a0 = rho*a[0] #dafcn returns [a,da]. It's easier to call than writing new afcn.

        return a0, Xass
    def dOm_aux(self, rho, temp_aux, mu, Psat, Xass0=None):
        """
        dOm_aux(rho, temp_aux, mu, Psat, Xass0)

        Method that computes the adimenstional Thermodynamic Grand potential
        at given density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        mu : float
            adimentional chemical potential at equilibrium
        Psat : float
            adimentional pressure [Pa]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        GPT: float
            Thermodynamic Grand potential [Pa]
        Xass : array
            computed fraction of nonbonded sites
        """
        a0, Xass = self.a0ad_aux(rho, temp_aux, Xass0)
        GPT = a0 - rho*mu + Psat
        return GPT, Xass

    def density_aux(self, temp_aux, P, state, rho0=None, Xass0=None):
        """
        density_aux(T, temp_aux, state, rho0, Xass0)
        Method that computes the density of the fluid at T, P

        Parameters
        ----------
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        rho0 : float, optional
            initial guess to compute density root [mol/m^3]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        density: float
            density [mol/m^3]
        Xass : array
            computed fraction of nonbonded sites

        """
        beta=temp_aux[0]
        tKelvin=1/(beta*kB/1E21)
        LIQ=2
        if(state=='L'):LIQ=3 #density routine does not need other props.
        # LIQ = 0,2 for vapor root, 1,3 for liquid root. 2,3 indicate zIter calculation (no uRes etc)
        PMPa=P/1E6
        iErr, rhoMol_cc,zFactor,uRes,aRes,chemPo=self.ChemPoTP(LIQ,tKelvin,PMPa)
        if iErr < 11:
            rho = rhoMol_cc*1E6
        else:
            rho = 8686
        return rho, Xass0
    def ares(self, rho, T, Xass0=None):
        """
        ares(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the fluid.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: float
           residual dimentionless Helmholtz free energy [Adim]
        """
        #temp_aux = self.temperature_aux(T)
        #a, Xass = ares(self, rho, temp_aux, Xass0)
        a, Xass = self.dares_drho( rho, T)
        #returns a[aRes_RT,(Z-1)/rhoMolecular], Xass
        aRes_RT=a[0]
        return aRes_RT, Xass0


##########################   END of EsdComp Class   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

class EsdFluid(object):
    '''
    class EsdFluid
    Creates an object that contains component(s) [and BIPs] for a fluid modeled with the ESD EOS.

    Parameters
    ----------
    component1 : object
        component created with component class

    Attributes
    ----------
    See EsdComp definition

    Methods
    -------
    add_component : adds a component to the mixture
    Kij : add Kij matrix for ESD EOS
    ci : computes cij matrix at T for SGT
    '''

    def __init__(self, component1,component2=None):
        self.components=[component1]
        if component2 is None:self.nc=1
        else:self.nc=2
        nc=self.nc
        self.etaMax=1/1.9-zeroTol
        if nc == 2: 
            self.components.append(component2)
            ## kij matrix
            self.KIJ0 = np.zeros([self.nc, self.nc])
            self.KIJ1 = np.zeros([self.nc, self.nc])
            id1=component1.idDippr
            id2=component2.idDippr
            inFile="c:/PGLWrapper/input/BipEsdMem2.txt"
            idMix=10000*id1+id2
            if(id1 < id2):idMix=10000*id2+id1
            with open(inFile,encoding='utf-8') as f:
                lines=f.readlines() # reads the entire file as lines.
                # q EOK(K)   b(cc/mo)   KadNm3  eDonEpsK    eAccEpsK    nDegree #eDs    #eAs    idCas   Name
                #bThere=False
                for line in lines:
                    part=line.split() # separates out the  
                    if int(part[0]) == idMix:
                        #bThere=True
                        self.KIJ0[0,1]=part[1]
                        self.KIJ0[1,0]=part[1]
                        self.KIJ1[0,1]=part[2]
                        self.KIJ1[1,0]=part[2]
                        break

    def add_component(self, component):
        """
        add_component(component)

        Method that adds a component to the mixture

        Parameters
        ----------
        component : object
            pure fluid created with component function
        """
        self.components.append(component)
        self.nc+=1

        # kij matrix
        self.KIJ0 = np.zeros([self.nc, self.nc])
        self.KIJ1 = np.zeros([self.nc, self.nc])
    """
    def __add__(self, new_component):
        if isinstance(new_component, component):
            self.add_component(new_component)
        else:
            raise Exception('You can only add components objects to an existing mixture')
        return self
    """
    def ci(self, T):
        """
        Method that computes the matrix of cij interaction parameter for SGT at
        given temperature.

        Parameters
        ----------
        T : float
            absolute temperature [K]

        Returns
        -------
        ci : array_like
            influence parameter matrix at given temperature [J m^5 / mol^2]
        """

        n = len(self.cii)
        ci = np.zeros(n)
        for i in range(n):
            ci[i] = np.polyval(self.cii[i], T)
        self.cij = np.sqrt(np.outer(ci, ci))
        return self.cij
    def RdfCalc(self,eta):
        """Compute the rdf and dlng/dlnEta consistent with the ESD EOS."""
        # Input:
        # eta = packing fraction
        # iEosOpt, bESD, bTPT cf. GlobConst
        # Output:
        # rdfContact= radial distribution function at contact
        # dAlpha=dLn(alpha)/dLn(rho)
        # Miscellaneous:
        # alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/deta) = 1+(eta/rdf)*(dRdf/deta)
        # dLng = (eta/rdf)*(dRdf/deta)     
        # iRdfOpt	- option for characterizing the radial distribution funciton at contact.	
        #			= 0, not specified => error
        #			= 1, ESD form
        #			= 2, Carnahan-Starling form
        #			= 3, activity coefficient model form (rdf=universal constant=g(eta=0.4)
        #			= 4, ESD form, geometric association rule.
        iRdfOpt=4		#ESD geometric form
        #alpha=eta*rdf*kAD*yHB => (eta/alpha)*(dAlpha/deta) = 1+(eta/rdf)*(dRdf/deta)
        #dLng = (eta/rdf)*(dRdf/deta) = dLng/dLneta
        denom=1.0-1.90*eta 
        #denom2=denom*denom
        void=1.0-eta
        void2=void*void
        #void4=void2*void2    

        dLng=1.90*eta/denom
        if(iRdfOpt==2):dLng=eta*( 3/void-1/(2-eta) )
        if(iRdfOpt==3):dLng= -1
        rdfContact=1/denom
        if(iRdfOpt==2):rdfContact=(1-eta/2)/void/void2
        if(iRdfOpt==3):rdfContact=1
        #d2Lng=1.90*eta/denom2
        #d2g=d2Lng-dLng
        dAlpha=1+dLng
        #dg_deta=1.90/denom2
        #dAlph_deta=dg_deta
        #if(iRdfOpt==2):dg_deta=(2.50-eta)/void4
        #if(iRdfOpt==2):dAlph_deta=3/void2-2/((2-eta)*(2-eta))
        return rdfContact, dAlpha

    def MEM2(self,bZiter,tKelvin,rhoMol_cc,xFrac): 
        """
        returns iErr,rhoMol_cc,zAssoc,aAssoc,uAssoc,rLnPhiAssoc
        !  PURPOSE:  COMPUTE THE EXTENTS OF ASSOCIATION (FA,FD) AND properties (zAssoc, lnPhiAssoc,...) given T,rho,x
        !  References: 
        !		Elliott, IECR, 61(42):15724 (2022). doi.org/10.1021/acs.iecr.2c02723
        !       Elliott, JCED, 69:458–471 (2024). doi.org/10.1021/acs.jced.3c00394
        !  INPUT:
        !  bZiter = True if only zAssoc,aAssoc are required (for zFactor iterations), False for lnPhiAssoc,uAssoc
        !  tKelvin = T(K), real*8
        !  xFrac   = component mole fractions, real*8
        !  nComps  = # of components, integer
        !  rhoMol_cc = molar density in gmol/cm3, real*8
        !  OUTPUT:
        !  zAssoc  = association effect on compressibility factor, real*8 
        !  aAssoc  = association effect on residual Helmholtz energy, real*8 
        !  uAssoc  = association effect on residual internal energy, real*8 
        !  rLnPhiAssoc  = association effect on log of fugacity coefficient, real*8
        !  iErr    = error code, integer (see errMsg vector below for code definitions).
        !  Incidentals that may be of interest: 
        !  FA,FD   = THE CHARACTERISTIC ASSOCIATION = (1/XAi-1)/RALPHi, Eqs.7-8
        !  RALPHi  = ROOT OF ALPHA WHERE ALPHAi=rho*bVoli*KADi*Ei*rdfContact, Eq. 5
        """
        global callCounter # "global" is how python makes a variable static. It does not initialize to zero.
        global ralphAmean,ralphDmean,etaOld,rdfOld
        global xOld# 
        global idOld#=[0]*self.nc
        if 'callCounter' not in globals(): callCounter=0 # this initializes callCounter iff it has not been initialized.
        nComps=self.nc
        xOld=np.zeros(nComps)
        ID=np.zeros(nComps)
        idOld=np.zeros(nComps)
        bVol=np.zeros(nComps)
        iErr=0
        zAssoc=0
        uAssoc=0
        aAssoc=0
        rLnPhiAssoc=np.zeros(nComps) # ln(fugacityCoefficient)

        maxTypes=2 # set minimum to 2 so arrays will stay as arrays instead of "1D arrays" being converted to scalar. Necessary to sustain loop structures.
        i= -1
        NoAssoc=True
        for componenti in self.components:
            i+=1
            if(componenti.bAssoc):NoAssoc=False
            ID[i]=(componenti.idDippr)
            bVol[i]=(componenti.bVolCC_mol)
            if np.max(componenti.assoc.nTypes)>maxTypes: maxTypes=np.max(componenti.assoc.nTypes)
        if(NoAssoc):return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
        idDiff=sum( np.abs(ID-idOld) )
        if( callCounter==0 or idDiff > 0) : #None of the following changes if the IDs don't change.
            nTypesMax=maxTypes
            nTypes=np.zeros(nComps,dtype=int)
            idType=np.zeros((nComps,nTypesMax),dtype=int)  #matrices need double parentheses.
            nDegree=np.zeros((nComps,nTypesMax),dtype=int)  #matrices need double parentheses.
            nAcceptors=np.zeros((nComps,nTypesMax),dtype=int)
            nDonors=np.zeros((nComps,nTypesMax),dtype=int)
            XA=np.ones((nComps,nTypesMax)) #matrices need double parentheses. FYI, dtype=np.float64 is equivalent to real*8.
            XD=np.ones((nComps,nTypesMax))
            XC=np.ones((nComps,nTypesMax))
            ralphA=np.zeros((nComps,nTypesMax))             #matrices need double parentheses.
            ralphD=np.zeros((nComps,nTypesMax))
            eAcceptorKcal_mol=np.zeros((nComps,nTypesMax))
            eDonorKcal_mol=np.zeros((nComps,nTypesMax))
            bondVolNm3=np.zeros((nComps,nTypesMax))
            i= -1
            for componenti in self.components:
                i+=1
                nTypes[i]=componenti.assoc.nTypes
                for j in range(componenti.assoc.nTypes):
                    idType[i,j]=componenti.assoc.idType[j]
                    nDegree[i,j]=componenti.assoc.nDegree[j]
                    nAcceptors[i,j]=componenti.assoc.nAsites[j]
                    nDonors[i,j]=componenti.assoc.nDsites[j]
                    eAcceptorKcal_mol[i,j]=componenti.assoc.epsA_kB[j]*RgasCal/1000
                    eDonorKcal_mol[i,j]=componenti.assoc.epsD_kB[j]*RgasCal/1000
                    bondVolNm3[i,j]=componenti.assoc.bondVolNm3[j]
        # Done unpacking components object!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        iErr=0
        zAssoc=0
        uAssoc=0
        aAssoc=0
        #fugAssoc=[0,0]
        #DoublePrecision, SAVE:: xOld(nmx),ralphAmean,ralphDmean,etaOld,rdfOld ! "SAVE" is equivalent to the static attribute
        #Integer, SAVE:: IDold(nmx)  								            ! static is equivalent to the "SAVE" attribute
        bAcid=False
        moreDonors=False
        #Character*123 errMsg(22)
        #data etaOld,rdfOld,ralphAmean,ralphDmean/0.3583,3.0,1.0,1.0/ !This makes these values static for Intel Compiler without STATIC.
        iErr=0
        errMsg=["0"]*22
        errMsg[ 1]=' MEM2: ERROR - NO CNVRG. ' # Warning because future calls (e.g., iterating Z) may converge and remove the issue.
        errMsg[ 2]=' MEM2: Failed bond site balance.' # Warning because some approximations might not require balance. 
        errMsg[11]=' MEM2: rdfContact<1'
        errMsg[12]=' MEM2: sqArg for ralphMean update.'
        errMsg[13]=' MEM2: sqArg(A or D)'
        errMsg[14]=' MEM2: XA,XD, or XC < 0' 
        errMsg[15]=' MEM2: etaOld < 0?' 
        errMsg[16]=' MEM2: Sorry. Derivatives beyond zAssoc and uAssoc have not been implemented yet.'
        errMsg[17]=' MEM2: sqArg < 0 for AA.'
        errMsg[18]=' MEM2: sqArg < 0 for DD.'
        errMsg[19]=' MEM2: sqArg < 0 for CC.'
        if(bZiter < 0):
            iErr=16
        FC=0 #initialize
        FA=0
        FD=0
        zAssoc=0
        aAssoc=0
        uAssoc=0
        rLnPhiAssoc=np.zeros(nComps)
        iErr=0  
        picard=0.830
        Ftol=2.E-7
        #Done with preliminary initializations. Let's write some code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        bVolMix=np.dot(xFrac,bVol)
        eta=rhoMol_cc*bVolMix
        rdfContact,dAlpha=self.RdfCalc(eta)	#dAlpha = dLnAlpha/dLnRho = 1+dLn(g)_dLn(eta)
        if(rdfContact < 1):
            iErr=11
            return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
        ralphAmax=0
        ralphDmax=0
        bAcid=False
        iComp= -1
        for componenti in self.components: # compute the ralphs based in input rho.
            iComp+=1
            for iType in range(componenti.assoc.nTypes): 
                if (componenti.assoc.idType[iType]==1603):
                    bAcid=True
                    epsCC=3*( eAcceptorKcal_mol[iComp,iType]+eDonorKcal_mol[iComp,iType] )/2
                    # ralphC is not subscripted because 1603 is the only acid type and always the same.
                    sqArg=rhoMol_cc*rdfContact*bondVolNm3[iComp,iType]*avoNum*( np.exp(epsCC/tKelvin/RgasCal*1000)-1.0 )
                    if(sqArg < 0):
                        iErr=19
                        return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
                    ralphC=np.sqrt(  sqArg  ) 
                    # Note the "3*" in the exponent, which makes the C bond "special."
                #NOTE: ralphA is initialized to zero, so it stays zero unless nAs > 0.
                if(componenti.assoc.nAsites[iType] > 0):
                    Fwertheim=np.exp(componenti.assoc.epsA_kB[iType]/tKelvin)-1
                    sqArg=rhoMol_cc*rdfContact*componenti.assoc.bondVolNm3[iType]*avoNum*Fwertheim
                    if(sqArg < 0):
                        iErr=17
                        return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
                    ralphA[iComp,iType]=np.sqrt(  sqArg  )
                    if(ralphA[iComp,iType] > ralphAmax):ralphAmax=ralphA[iComp,iType] # Store the max value for later.
                if(componenti.assoc.nDsites[iType] > 0):
                    Fwertheim=np.exp(componenti.assoc.epsD_kB[iType]/tKelvin)-1
                    sqArg=rhoMol_cc*rdfContact*componenti.assoc.bondVolNm3[iType]*avoNum*Fwertheim
                    if(sqArg < 0):
                        iErr=18
                        return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
                    ralphD[iComp,iType]=np.sqrt(  sqArg  )
                    if(ralphD[iComp,iType] > ralphDmax):ralphDmax=ralphD[iComp,iType]
        # Done summing over components to compute ralphs and averages
        FA0=0
        FD0=0
        FC0_ralphC=0
        avgNAS=0
        avgNDS=0
        iComp= -1
        for componenti in self.components:
            iComp+=1
            for j in range(componenti.assoc.nTypes):
                FA0=FA0+xFrac[iComp]*componenti.assoc.nDegree[j]*componenti.assoc.nDsites[j]*ralphD[iComp,j]	
                FD0=FD0+xFrac[iComp]*componenti.assoc.nDegree[j]*componenti.assoc.nAsites[j]*ralphA[iComp,j]
                if(componenti.assoc.nDegree[j]==1603):FC0_ralphC+=xFrac[iComp]*componenti.assoc.nDegree[j] #one C allowed on type 1603 => nCarboxy(i,j)=1 always.
                avgNDS=avgNDS+xFrac[iComp]*componenti.assoc.nDegree[j]*componenti.assoc.nDsites[j]	
                avgNAS=avgNAS+xFrac[iComp]*componenti.assoc.nDegree[j]*componenti.assoc.nAsites[j]
        moreDonors=False
        if(avgNAS < avgNDS):moreDonors=True
        if(FA0==0 or FD0==0):return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc #no need to calculate if either no donors or no acceptors
        ##########################
        newMix=False
        if(callCounter==0):
            newMix=True
            etaOld=0.01 # to bypass etaOld < zeroTol condition
        else: # ie. if callCounter > 0...
            xDiff=np.sum( (xFrac-xOld)**2 )
            idDiff=np.sum( np.abs(ID-idOld) )
        if( idDiff !=0 or xDiff > Ftol/10 or etaOld < zeroTol  ):newMix=True
        if( newMix):
            #Only use default estimate if compounds or composition have changed.
            ralphAmean=ralphAmax
            ralphDmean=ralphDmax
            if(moreDonors):
                ralphDmean=ralphDmean*avgNAS/avgNDS
            else:
                ralphAmean=ralphAmean*avgNDS/avgNAS
        else: #if( SUM(xFrac(1:nComps)) < 0)then
            if(etaOld > 0): # adapt old values.to accelerate Z iterations
                sqArg=eta/etaOld*rdfContact/rdfOld
                if(sqArg < 0):
                    iErr=12
                    return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
                ralphAmean=ralphAmean*np.sqrt(sqArg)
                ralphDmean=ralphDmean*np.sqrt(sqArg)
            else:
                iErr=15
                return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
        #FAold =0 # low guess
        #errOld= -FA0 # ie. when guessing FA=FD=0, sumD=FA0. So, FA(guess)-sumD = 0 - FA0 
        FA    =1 # Take FA0 and FD0 as initial guesses in secant iteration.
        FD    =1 # 
        nIter=0
        #avgFo=(ralphDmean*FD0+ralphAmean*FA0)/2
        itMax=88
        error=123
        nIter= -1
        while(np.abs(error)>Ftol and nIter<itMax):
            nIter+=1
            #if(nIter > 33):picard=0.50
            delFo=(ralphDmean*FD0-ralphAmean*FA0)
            delFA=1+delFo
            delFD=1-delFo
            sqArgA=delFA*delFA+4*ralphAmean*FA0
            sqArgD=delFD*delFD+4*ralphDmean*FD0
            if(sqArgA < zeroTol or sqArgD < zeroTol):
                iErr=13
                return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
            FA = 2*FA0/( delFA+np.sqrt(sqArgA) )
            FD = 2*FD0/( delFD+np.sqrt(sqArgD) )
            sumA=0
            sumD=0
            for i in range(nComps): #If ralphMeans are correct, then FA,FD are correct and we should get sumA=FD and sumD=FA.
                for j in range(nComps):
                    sumA=sumA+xFrac[i]*nDegree[i,j]*nAcceptors[i,j]*ralphA[i,j]/( 1+FA*ralphA[i,j] ) #*FA !=FD=FD0/(1+ralphAmean*FA)
                    sumD=sumD+xFrac[i]*nDegree[i,j]*nDonors[i,j]   *ralphD[i,j]/( 1+FD*ralphD[i,j] ) #*FD !=FA=FA0/(1+ralphDmean*FD) 
            errA=FA-sumD #/FD
            errD=FD-sumA #/FA
            error=np.max([np.abs(errA),np.abs(errD)]) # You need to put brackets around the sequence for np.max() to work.
            #FAold=FA
            #FDold=FD
            #Compute new ralphMeans.
            ralphDmean= picard*(-1+FA0/sumD)/(FD+1E-9)+(1-picard)*ralphDmean #Using new sumD to compute new ralphMeans.
            if(ralphDmean < 0):ralphDmean=(-1+FA0/sumD)/(FD+1E-9)
            if(ralphDmean < 0):ralphDmean=zeroTol
            ralphAmean= picard*(-1+FD0/sumA)/(FA+1E-9)+(1-picard)*ralphAmean #add 1E-9 to avoid possible zero divide if FD->0
            if(ralphAmean < 0):ralphAmean=(-1+FD0/sumA)/(FA+1E-9)
            if(ralphAmean < 0):ralphAmean=zeroTol
            #if(nIter==1):FA1=FA
        #while abs(error) > FTol. 
        if(nIter>itMax-1):iErr=1	#warning level
        #  fAssoc ITERATION HAS CONCLUDED
        if(bAcid):
            sqArg=1+np.sqrt(1+4*FC0_ralphC*ralphC*ralphC)
            if(sqArg < 0):
                iErr=20
                return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
            FC=2*FC0_ralphC*ralphC/( sqArg )
            XCtemp=1/(1+FC*ralphC)
            if(XCtemp < zeroTol):
                iErr=14
                return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
        hAD=0
        hDA=0
        hCC=0
        for i in range(nComps):  # XA and XD USE Assoc.
            for j in range(nTypes[i]):  # XA and XD USE Assoc.
                XA[i,j]= 1/(1+ralphA[i,j]*FA)
                XD[i,j]= 1/(1+ralphD[i,j]*FD)
                hAD=hAD+xFrac[i]*nDegree[i,j]*nAcceptors[i,j]*(1-XA[i,j])
                hDA=hDA+xFrac[i]*nDegree[i,j] * nDonors[i,j] *(1-XD[i,j])
                XC[i,j]= 1
                if(idType[i,j]==1603):
                    XC[i,j]=XCtemp
                    hCC=hCC+xFrac[j]*nDegree[i,j] * (1-XC[i,j])
        zAssoc= -dAlpha*(hAD+hDA+hCC)/2
        for i in range(nComps):  # XA and XD USE Assoc.
            sumLnXi=0.0
            for j in range(nTypes[i]):  # XA and XD USE Assoc.
                if( XA[i,j] < zeroTol or XD[i,j] < zeroTol or XC[i,j] < zeroTol ):
                    iErr=14
                    return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
                sumLnXi=sumLnXi+nDegree[i,j]*(  nAcceptors[i,j]*np.log( XA[i,j] )+nDonors[i,j]*np.log( XD[i,j] )+ np.log( XC[i,j] )  )	#
            aAssoc=aAssoc+xFrac[i]*sumLnXi	# ln(XD)=ln(XA) so *2.
        aAssoc=  aAssoc+(hAD+hDA+hCC)/2  # Eqs. 1 & 40. Validation that dSUM(xi*SumLnXi)/dBeta=0 is given at end of paper.
        xOld=xFrac
        idOld=ID
        etaOld=eta
        rdfOld=rdfContact
        if(bZiter):return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc ###########################################################################
        #aAssoc= 0 !already initialized at the top.
        betadFA_dBeta=0
        betadFD_dBeta=0
        betadFC_dBeta=0
        for i in range(nComps):  # XA and XD USE Assoc.
            sumLnXi=0.0
            for j in range(nTypes[i]):  # XA and XD USE Assoc.
                if( XA[i,j] < zeroTol or XD[i,j] < zeroTol or XC[i,j] < zeroTol ):
                    iErr=14
                sumLnXi+= nDegree[i,j]*(  nAcceptors[i,j]*np.log( XA[i,j] )+nDonors[i,j]*np.log( XD[i,j] )+ np.log( XC[i,j] )  )	#
                bepsA=(eAcceptorKcal_mol[i,j]/RgasCal*1000/tKelvin)
                dLnAlpha_dLnBeta=np.exp(bepsA)	# Eqs. 44,45
                if(bepsA > 1.E-4):dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsA)/(dLnAlpha_dLnBeta-1)
                betadFD_dBeta+= xFrac[i]*nDegree[i,j]*nAcceptors[i,j]*XA[i,j]*ralphA[i,j]*0.50*dLnAlpha_dLnBeta
                bepsD=(eDonorKcal_mol[i,j]/RgasCal*1000/tKelvin)
                dLnAlpha_dLnBeta=np.exp(bepsD)
                if(bepsD > 1.E-4):dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsD)/(dLnAlpha_dLnBeta-1)
                betadFA_dBeta+= xFrac[i]*nDegree[i,j] * nDonors[i,j] *XD[i,j]*ralphD[i,j]*0.5*dLnAlpha_dLnBeta  
                # the 0.5 comes from dLnRalph = 0.5*dLnAlpha
                if(idType[i,j]==1603):
                    bepsC=3*(eDonorKcal_mol[i,j]+eDonorKcal_mol[i,j])/(2*RgasCal/1000*tKelvin)
                    dLnAlpha_dLnBeta=np.exp(bepsC)
                    if(bepsC > 1.E-4):dLnAlpha_dLnBeta=dLnAlpha_dLnBeta*(bepsC)/(dLnAlpha_dLnBeta-1)
                    betadFC_dBeta=betadFC_dBeta+xFrac[i] *nDegree[i,j] *XC[i,j]*ralphC*dLnAlpha_dLnBeta
            #endfor j
            rLnPhiAssoc[i]=sumLnXi-(hAD+hDA+hCC)/2*(dAlpha-1)*bVol[i]/bVolMix # Eq. 42
        #endfor i...
        #aAssocPas=aAssoc
        #(hAD+hDA)/2=hAD = FA*FD by Eq 40.
        uAssoc= -( FA*betadFD_dBeta+FD*betadFA_dBeta+FC*betadFC_dBeta/2 ) # Eqs. 43,44
        if(np.abs(hAD-hDA) > Ftol):
            iErr=2 #Warning
            #if(LOUDER)write(dumpUnit,*)'MEM2: Failed bond site balance.'
        return iErr, zAssoc, aAssoc, uAssoc, rLnPhiAssoc
######################################################### END MEM2 ############################################################################
    def ChemPoTV(self,bZiter,gmol,tKelvin,Vcm3): 
        """Returns iErr, PMPa,zFactor,uRes,aRes,chemPo[]"""
        iErr=0
        cqFactor=4/(4-1.9) # = 1.90476... where c-1 = (q-1)*cqFactor.
        K1=1.7745
        K2=1.0617
        totMols=sum(gmol)
        rho=totMols/Vcm3
        xFrac=[]
        bVol=[]
        cShape=[]
        EOK = np.zeros([self.nc, self.nc])
        Y= np.zeros([self.nc, self.nc])
        qb= np.zeros([self.nc, self.nc])
        YQb= np.zeros([self.nc, self.nc])
        cb= np.zeros([self.nc, self.nc])
        chemPo= np.zeros(self.nc)
        bMix=0
        cMix=0
        K1YbMix=0
        i=0
        for componenti in self.components:
            xFrac.append(gmol[i]/totMols)
            bVol.append(componenti.bVolCC_mol)
            bMix+=xFrac[i]*bVol[i]  # we need to define all bVol before referencing bVol[j]
            cShape.append( 1+(componenti.qShape-1)/cqFactor ) # q-1= 1.90476*(c-1) => c = 1+(q-1)/1.90476
            cMix+=xFrac[i]*cShape[i]
            K1YbMix+=xFrac[i]*K1*( np.exp(componenti.eps_kB/tKelvin)-K2 )*bVol[i]      #1991 form, overwritten if applying 1990 form
            i+=1
        YQbMix=0
        cbMix=0
        i=0
        for componenti in self.components:
            j=0
            for componentj in self.components:
                KijBip=0
                if(i != j):KijBip=self.KIJ0.item(i,j)+self.KIJ1.item(i,j)/tKelvin
                #epsij=np.sqrt(componenti.eps_kB*componentj.eps_kB)*(1-KijBip)
                EOK[i,j]=np.sqrt(componenti.eps_kB*componentj.eps_kB)*(1-KijBip)
                Y[i,j]=np.exp(EOK[i,j]/tKelvin)-K2
                qb[i,j]= (componenti.qShape*bVol[j] + componentj.qShape*bVol[i]) / 2
                YQb[i,j]=qb[i,j]*Y[i,j]
                cb[i,j] = (cShape[i]*bVol[j] + cShape[j]*bVol[i]) / 2 
                # e.g. (x1*c1+x2*c2)*(x1*b1+x2*b2) = x1^2*c1*b1+x1*x2*(c1*b2+c2*b1)+x2^2*b2^2
                YQbMix+=YQb[i,j]*xFrac[i]*xFrac[j]	   
                cbMix += cb[i,j]*xFrac[i]*xFrac[j]    #note: above means <c>=sum(xi*ci) and <b>=sum(xj*bj) and <cb>=<c>*<b> 
                j+=1
            i+=1
        iErrMEM, zAssoc, aAssoc, uAssoc, fugAssoc = self.MEM2(bZiter,tKelvin,rho,xFrac)
        eta=rho*bMix
        voidFrac=1-1.90*eta
        denom=voidFrac
        zRep= 4*cMix*eta/denom
        zAtt= -9.5*YQbMix*rho/(1+K1YbMix*rho)
        zFactor=(1+zRep+zAtt+zAssoc)
        PMPa=zFactor*Rgas*rho*tKelvin
        aRep= -4/1.9*np.log(voidFrac)*cMix
        aAtt= -9.5*YQbMix/K1YbMix*np.log(1+K1YbMix*rho)
        aRes=aRep+aAtt+aAssoc #-DLOG(Z) !don't subtract log(z) for aRes(T,V). Important for EAR.
        uRes=0 # set default value to avoid undefined error.
        if( zFactor < zeroTol):  # Z < 0 is no problem given Vtot because ln(Z) is not relevant.  JRE 20210724
            iErr=3	  # warning level because another call might produce Z > 0.
        if(bZiter or iErr>10):return iErr, PMPa,zFactor,aRes,uRes,chemPo  # don't need the rest if bZiter.
        BdYb_dB=0
        BdYbq_dB=0
        YQbAvg=[]
        #cbAvg=[]
        i= -1
        for componenti in self.components:
            i+=1
	        #ralph[i]=SQRT(alphAD(I,I))	 ! ralph is computed in alpsolEz.
            BdYb_dB+=xFrac[i]*bVol[i]*EOK[i,i]/tKelvin*(Y[i,i]+K2)*K1  # = Beta*d<k1Yb>/dBeta 
            YQbAvg.append(0) # appends the next element in YQbAvg and initiates to zero.
            #cbAvg.append(0) 
            j= -1
            for componentj in self.components:
                j+=1
                BdYbq_dB+=xFrac[j]*xFrac[i]*qb[i,j]*EOK[i,j]/tKelvin*(Y[i,j]+K2) # = Beta*d<Ybq>/dBeta
                YQbAvg[i]+=YQb[i,j]*xFrac[j]
                #cbAvg[i]+=cb[[i,j]]*xFrac[j]
        uAtt= -9.5*YQbMix*rho/(1+K1YbMix*rho)*BdYb_dB/K1YbMix + aAtt*(BdYbq_dB/YQbMix-BdYb_dB/K1YbMix) #FYI:I had omitted the 2nd term previously
        uRes=uAtt+uAssoc
        fugRep=[]
        fugAtt=[]
        i=-1
        for componenti in self.components:
            i+=1
            fugRep.append( aRep*cShape[i]/cMix + zRep*bVol[i]/bMix ) # For pure i, FugRepi= -4ci/1.9*ln(1-1.9eta) + 4ci*eta/(1-1.9eta) ) 
            fugAtt.append( aAtt*( 2*YQbAvg[i]/YQbMix-K1*Y[i,i]*bVol[i]/K1YbMix )+zAtt*K1*Y[i,i]*bVol[i]/K1YbMix ) #91-pres form,
            #fugAssoc[i] returned from MEM2
            chemPo[i]=( fugRep[i]+fugAtt[i]+fugAssoc[i] )  # -np.log(Z)  Don't subtract ln(Z) when given Vtot as independent variable.
            #rLnGamRep[i]=fugRep[i]-zRep*bVol[i]/bMix  ! cf. Bala and Lira (2016), Eqs A6-A14. to correct from constant volume to P.
            #rLnGamAtt[i]=fugAtt[i]-zAtt*bVol[i]/bMix
	        #rLnGamAssoc[i]=fugAssoc[i]-zAssoc*bVol[i]/bMix
        return iErr,PMPa,zFactor,aRes,uRes,chemPo
##########################   End of ChemPoTV   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def ChemPoTP(self,LIQ,xFrac,tKelvin,PMPa): 
        """Returns iErr, rhoMol_cc,zFactor,aRes,uRes,chemPo[], eta """
        # LIQ = 0,2 for vapor root, 1,3 for liquid root.
        bLiq=False
        if(LIQ==1 or LIQ==3):bLiq=True
        iErr=0
        rhoMol_cc=0.
        zFactor=1.
        uRes=0.
        aRes=0.
        chemPo=np.zeros(self.nc)
        #NC=self.nc
        #NOTE: iErrTmin is checked in FuVtot
        eta=0
        if(bLiq):eta=self.etaMax/1.15
        sumx=np.sum( xFrac )
        if(np.abs(sumx-1) > 1e-8):
            iErr=11
            return iErr, rhoMol_cc,zFactor,uRes,aRes,chemPo,eta
        #INITIATE SECANT ITERATION ON rho
        bVol=[]
        bMix=0
        i= -1
        for componenti in self.components:
            i+=1
            bVol.append(componenti.bVolCC_mol)
            bMix+=xFrac[i]*bVol[i]  # we need to define all bVol before referencing bVol[j]
        Pb_RT=PMPa*bMix/(Rgas*tKelvin)
        #GUESS FOR rho
        eta=Pb_RT/1.05  #NOTE: Pb_RT > 1 can happen when Z >>1, like at GPa.
        if(bLiq or eta>self.etaMax):eta=self.etaMax/1.15
        rho=eta/bMix
        #if(eta > 1/1.9 .and. LOUDER)write(dumpUnit,*)'FugiEsd:etaInit > etaMax. P,T=',pMPa,tKelvin 
        bZiter=True # chemPo[] calculations are skipped for isZiter=1
        iErr, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,xFrac,tKelvin,1/rho)
        #Call FuEsdVtot(isZiter,tKelvin,1/rho,xFrac,NC,FUGC,zFactor,Ares,Ures,iErr)
        if(iErr > 10):
            return iErr, rhoMol_cc,zFactor,uRes,aRes,chemPo,eta
        etaOld=eta
        errOld=Pb_RT-eta*zFactor
        eta=etaOld/1.150
        #IF (eta < 0 .and. LOUD) write(dumpUnit,31)LIQ
        #if(initial==1.and.LOUD)write(dumpUnit,*)'FugiEsd: initial eta,err',etaOld,errOld
        itMax=77
        errBesteta=1234
        for nIter in range(itMax):
            rho=eta/bMix
            iErr, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,xFrac,tKelvin,1/rho)
            if(iErr > 10):break
            ERR=Pb_RT-eta*zFactor
            CHNG=ERR/(ERR-errOld)*(eta-etaOld)
            #if(initial==1.and.LOUDER)write(dumpUnit,'(a,2e11.4,3f10.5)')'FugiEsd eta,Z', eta,zFactor
            #if(initial==1.and.LOUDER)write(dumpUnit,'(a,f8.5,e11.4,i3,9f8.3)')'FugiEsd eta,CHNG,niter',eta,CHNG,niter 
            etaOld=eta
            errOld=ERR
            #  LIMIT THE CHANGE IN Density for liquid.
            if(bLiq and np.abs(CHNG/etaOld) > 0.15):CHNG=np.sign(CHNG)*0.15*etaOld
            #  Low eta must move from zero, so you can't limit its % change
            if(not bLiq and np.abs(CHNG) > 0.020):CHNG=np.sign(CHNG)*0.02
            eta=eta-CHNG
            if(np.abs(CHNG) < errBesteta):
                etaBest=eta
                errBesteta=np.abs(CHNG)
            if(eta < 0 or eta > 1/1.9):eta=etaOld-np.sign(CHNG)*0.1*etaOld
            if(np.abs(CHNG) < 1.E-9 and eta > 0):break  # Don't use CHNG/eta. Converge quickly to ideal gas if P->0, ~9 sigfigs if liquid
        ##########################   Iteration Concluded    !!!!!!!!!!!!!!!!!!!!!!!!!
        if(eta < 0 or eta > 1/1.9):
            iErr=17
        rho=eta/bMix
        rhoMol_cc=rho 
        #!  ITERATION ON rho HAS CONCLUDED.  
        #zStore=zFactor   
        bZiter=False
        #One last call to get uRes,chemPo.
        iErrF, Pcalc,zFactor,aRes,uRes,chemPo=self.ChemPoTV(bZiter,xFrac,tKelvin,1/rho)
        if(zFactor < zeroTol):
            iErrF=11
        if(iErrF > 0 or nIter > itMax-1 or eta < 0): # if iErr still > 0 on last iteration, then declare error.
            iErr=iErrF 
            if(nIter > itMax-1):iErr+=14
            eta=etaBest
        chemPo-=np.log(zFactor)	 #!Must subtract ln(Z) when given P as independent variable.
        return iErr,rhoMol_cc,zFactor,aRes,uRes,chemPo,eta
    ##########################   End of ChemPoTP   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def PvpEarMix(self,component,tKelvin):
        """ Returns iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
	    !COMPUTE VAPOR PRESSURE GIVEN tKelvin using the equal area rule.  
	    !Ref: Eubank et al, IECR 31:942 (1992)"""
        errMsg=["0"]*22 #!,errMsgPas
        errMsg[ 2]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[ 3]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[ 4]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[ 9]='PvpEar: Warning from fugacity calculation. Probably T < Tmin.'
        errMsg[11]='PvpEar: No spinodal max/min'
        errMsg[12]='PvpEar: Psat iteration did not converge'
        errMsg[13]='PvpEar: Liquid FUGI call failed on last iteration'
        errMsg[14]='PvpEar: Vapor  FUGI call failed on last iteration'
        errMsg[15]='PvpEar: zVap=zLiq on last iteration'
        errMsg[16]='PvpEar: ChemPoTP returned T < Tmin error'
        errMsg[17]='PvpEar: Calculated Psat < 0.0001 MPa error'
        errMsg[18]='PvpEar: Tr > 1-zeroTol'
        errMsg[19]='PvpEar: Critical Error from fugacity calculation'
        errMsg[20]='PvpEar: rhoVap/rhoLiq < 0'
        #tMin=component.TcK*0.25
        #NC=1 # Ltd to single compound.
        xFrac=np.zeros(2) # force 2D to avoid python declaring 1D array as scalar.
        xFrac[0]=1
        iErr=0 
        PMPa=0
        chemPot=0
        rhoLiq=0
        rhoVap=0
        uSatL=0
        uSatV=0
        zCritEff=component.Zc
        zCritEff=component.ZcEos  # For ESD, this is close enough
        rhoCrit=component.PcMPa/(zCritEff*Rgas*component.TcK)
        rhoVap=rhoCrit*1.000 #-ve value on input means use default value for etaHi
        Tr=tKelvin/component.TcEos
        if(Tr > 1-zeroTol):
            iErr=18
            return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        elif(Tr > 0.85):
            # etaL=2A(1-Tr)-etaV; 2A(1-Tr)- 2etaV= B*(1-Tr)^.5; etaV = 0.1672*(1-Tr)-1.22(1-Tr)^.5/2; Fit for EsdCO2.
            etaVap=component.etacEos+0.1672*(1-Tr)-1.22/2*np.sqrt(1-Tr)
            etaLiq=(2*.1672*(1-Tr)-etaVap)
            rhoLiq=rhoLiq/component.bVolCC_mol
            pOld=component.PcMPa*10**( 7/3*(1+component.acen)*(1-1/Tr) ) # SCVP
            #pEubank=Rgas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) ) #EAR method of Eubank and Hall (1995)
            #pOld=pEubank/10  #JRE: I find Eubank overestimates P sometimes. It even gets higher than pMax
        else: # If Tr < 0.85, compute pSatItic estimate  FPE,501:112236(19), JCP,110:3043(99)
            pOld=0.01 #If p < 0, then the vapor density goes negative and log(rhoLiq/rhoVap) is indeterminate.
            #NOTE: Don't use pOld=0 when pOld > 0. You might get an error because vdw loop doesn't cross zero#
            iErrF,rhoLiq,zLiq,aResLiq,uSatL,chemPo,etaLiq=self.ChemPoTP(1,xFrac,tKelvin,pOld)
            if(iErrF > 10):
                iErr=19
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            elif(iErrF > 0):
                iErr=9
            elif(zLiq < zeroTol):
                iErr=17
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
                #AresLiq=Ares_RT						#AresVap  +  Zvap-1 -ln(Zvap) =AresLiq+ZLiq-1-ln(ZLiq)							 
                #rhoLiq=pOld/(zLiq*Rgas*tK)		#=>	#B2*rhoVap+B2*rhoVap+ln(rhoVap/rhoLiq) =AresLiq+0-1  
            rhoVap=rhoLiq*np.exp( aResLiq-1 ) #pSatItic, FPE, 501:112236(19), Eq 19. Ares_RT from GlobConst
            pSatItic=rhoVap*Rgas*tKelvin
            etaVap=rhoVap*component.bVolCC_mol
            #call FuVtot(1,tK,pSatItic,xFrac,NC,0,FUGC,zVap,ier)
            #call FuVtot(1,tK,1/rhoVap,xFrac,NC,FUGC,zVap,aRes,uRes,iErrF) #isZiter=1=>vapor Z with no fugc calculation. zVap >> zero.
            bZiter=False
            iErrF,Pcalc,zVap,aResVap,uSatV,chemPot=self.ChemPoTV(bZiter,xFrac,tKelvin,1/rhoVap)
            if(iErrF > 10):
                iErr=19
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            elif(iErrF > 0):
                iErr=9
            elif(zLiq < zeroTol):
                iErr=17
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            #B2cc_mol*rhoVap=(zVap-1)
            rhoVap=rhoLiq*np.exp( aResLiq-1-2*(zVap-1) )# => rhoVap=rhoLiq*EXP(AresLiq-1-2*B2*rhoVap).
            pSatItic=rhoVap*Rgas*tKelvin*zVap
            if(pSatItic > zeroTol):
                pTest=pSatItic
                zLiq=pSatItic/(rhoLiq*Rgas*tKelvin) # if Psat > 0 then zLiq > 0
                #FugcLiq=aResLiq+zLiq-1-np.log(zLiq)
            else:
                pTest=zeroTol*10
            #pEuba2=Rgas*tK/(1/rhoVap-1/rhoLiq)*( Ares_RT- 0 +DLOG(rhoLiq/rhoVap) )	#EAR method of Eubank and Hall (1995)
            #NOTE: When 1/rhoVap >> 1/rhoLiq, pEuba2=Rgas*tK*rhoVap*( Ares_RT-ln(rhoVap/rhoLiq) )=pTest*( Ares_RT - (Ares_RT-1) )
            pOld=pTest
        # end if:Tr > 0.85...
        etaLiq=rhoLiq*component.bVolCC_mol
        if(pOld < zeroTol):pOld=zeroTol
        pBest=1234
        fBest=1234
        itMax=33
        for iter in range(itMax): #iterate on pOld according to Eubank criterion
	        #write(dumpUnit,*)'calling fugi for liquid iteration.'
	        #call FUGI(tK,pOld,xFrac,NC,1,FUGC,zLiq,ier) #LIQ=3=>liquid root with no fugc calculation
	        #CALL FugiTP( tK,Pold,xFrac,NC,1,rhoLiq,zLiq,aResLiq,FUGC,uSatL,iErrF )
            iErrF,rhoLiq,zLiq,aResLiq,uSatL,chemPot,etaLiq=self.ChemPoTP(3,xFrac,tKelvin,pOld)
            if(iErrF > 0):iErr=3 #declare error but don't stop. if future iterations give valid fugi(), then no worries
	        #chemPotL=FUGC(1)
	        #aDepLiq=aRes_RT
	        #rhoLiq=pOld/(zLiq*Rgas*tK)
	        #rhoLiq=etaPass/bVolCc_mol(1) # crazy precision is required when etaLiq > 0.98 (and Z->0). e.g. pentanoic acid (1258)
	        #write(dumpUnit,*)'calling fugi for liquid iteration.'
	        #call Fugi(tK,pOld,xFrac,NC,0,FUGC,zVap,ier) #LIQ=2=>vapor root with no fugc calculation
	        #CALL FugiTP( tK,Pold,xFrac,NC,0,rhoVap,zVap,aResVap,FUGC,uSatV,iErrF )
	        #if(iErrF > 0)ierCode=4 #declare warning error but don't stop. if future iterations give valid fugi(), then no worries.
            iErrF,rhoVap,zVap,aResVap,uSatV,chemPot,etaVap=self.ChemPoTP(2,xFrac,tKelvin,pOld)
            if(iErrF > 0):iErr=4 #declare error but don't stop. if future iterations give valid fugi(), then no worries
	        #chemPotV=FUGC(1)
	        #aDepVap=aRes_RT
	        #rhoVap=etaPass/bVolCc_mol(1)
            if( (rhoVap/rhoLiq) < 0 ):
                iErr=20
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            if(np.abs(zVap-zLiq) < 1e-3):
                iErr=15
                return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
            if(pOld < 1E-8 and iter > 2):
                #etaLiq=rhoLiq*component.bVolCC_mol
                rhoVap=rhoLiq*np.exp( aResLiq+zLiq-1-aResVap+zVap-1 ) #pSatItic, FPE, 501:112236(19), Eq 19.aDepVap=/=0 for carbo acids
                pTest=rhoVap*Rgas*tKelvin*zVap
                iErr=1 # assign warning in this case. 
                break
	        #rLogRhoRat=DLOG(rhoVap/rhoLiq)	#for debugging
            pTest=pOld*( aResLiq-aResVap-np.log(rhoVap/rhoLiq) )/(zVap-zLiq)	#EAR method of Eubank and Hall (1995)
	        #fErr= aDepLiq+zLiq-1 -(aDepVap+zVap-1) -DLOG(zLiq/zVap)
            fErr= aResLiq+zLiq-1 -(aResVap+zVap-1) -np.log(rhoVap/rhoLiq) 
            if( np.abs( fErr) < fBest):	# for very low pSat, like propane.
                fBest=np.abs( fErr)
                pBest=pTest
            change=pTest-pOld
            tol=1E-6
            if(pTest > 0):
                pOld=pTest
            else:
                pOld=pOld/10
	        #pOld =Rgas*tK/(1/rhoVap-1/rhoLiq)*( aDepLiq-aDepVap +DLOG(rhoLiq/rhoVap) )	#EAR method of Eubank and Hall (1995)
            if(np.abs(change/pOld) < tol):break   
            # end for itMax
        PMPa=pOld
        if(np.abs(change/PMPa) > tol or iter >= itMax):
            iErr=2
            PMPa=pBest
            return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        #iter=0                                  
        #check that Fugi did not give warning on last call
        if(iErr > 10):return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
        #if(ierCode==0) call FuVtot(isZiter,tK,1/rhoLiq,xFrac,NC,FUGC,ZLiq,iErrFu) #one last call to fugi for chemPo and debugging.
        chemPo[0]=86.86
        if( zLiq > 0):
            chemPo[0]=aResLiq+zLiq-1 -np.log(zLiq)
        else:
            iErr=20
        chemPot=chemPo[0]
        return iErr,PMPa,chemPot,rhoLiq,rhoVap,uSatL,uSatV
##########################   END of PvpEar()   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
