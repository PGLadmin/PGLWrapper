# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""
#print("Hello world")
import numpy as np
#import sys
#sys.path.append('C:\\Users\\ellio\\sgtpy\\vrmie_pure')
import os
print(f"JRE: Current Working Directory: {os.getcwd()}")
# Get the directory of the current script
#script_dir = os.path.dirname(os.path.abspath(__file__))
# Change the working directory to the script's directory
#os.chdir(script_dir)
#print(f"Current Working Directory: {os.getcwd()}")
#eosDir=script_dir+'/vrmie_pure'
#os.chdir(eosDir)
#print(f"EOS Directory (current): {os.getcwd()}")

#from sgtpy import saftvrmie, component #, mixture
#from sgtpy import ESDpy 
from sgtpy.sgt import sgt_pure
import ESD  #JRE: absolute import
import pandas as pd
from GlobConst import PGLInputDir

#"""
#import scipy as sp
#import Cython as cp
#import sgtpy as saft
#"""
#testCode
import keyboard
def pauseCheck(force=None):
    if(force==1):input("Press enter to continue.")
    if keyboard.is_pressed('Esc'):
        input("Press enter to continue.")
testComp=ESD.EsdComp(1997)
atomsTot=testComp.CHBCFINOSS.sum
nCarb=testComp.CHBCFINOSS[0]
nFluoro=testComp.CHBCFINOSS[4]
F_C1=nFluoro/(nCarb+1)
print("nFluoro/(nCarb+1)=",F_C1)
#pauseCheck(1)

#mseg	lrep	epsilon/kB [K]	sigma [A]	cii [J m5 mol-2}	Omega

ULdf=pd.read_csv(PGLInputDir+'ParmsPrTcJaubert.txt',sep=r'[\t,]', skipinitialspace=True,header=0,engine='python')
ULdf.columns = ULdf.columns.str.strip() # Remove leading and trailing whitespace from column names
print(ULdf.columns)
ULdf.set_index('idDippr',inplace=True)
idDip=1

#PPData = ULdf.set_index(ULdf.columns[0]).T.to_dict()
#del ULdf # free up the memory of the ULdf dataframe.
print("Loading stDB....")
stDB=pd.read_excel('c:/users/jarrell.elliot/Downloads/SurfTenDBcleaned.xlsx') #,nrows=44)
print(stDB.head())
print(stDB.columns)
outfile=open('EsdOut.txt',"w")
outfile.close()
with open('EsdOut.txt',"a") as outfile: #The "with" will close the file properly even if there is a crash.
    #print("     T(K)       Psat(kPa)       rhoL(g/cc)       rhoV(g/cc)")
    print("idDippr,acen,acenCalc,T(K),Psat(kPa),rhoL(g/cc),rhoV(g/cc),stCalc,")
    outfile.write("idDippr \t T(K) \t Psat(kPa) \t rhoL(g/cc) \t rhoV(g/cc) \t stCalc \t st(mN/m) \t    Tr  \n")
    i=-1
    idOld=-1
    stDB['stCalc']=0.
    stDB['Tr']=0.
    stDB['bVol']=0.
    stDB['eps_kB']=0.
    stDB['qShape']=0.
    stDB['bondVol']=0.
    stDB['epsD_kB']=0.
    stDB['epsA_kB']=0.
    stDB['pctDev']=0.
    stDB['nHetero_C.1']=0.
    stDB['nHalo_C.1']=0.
    stDB['nF_C.1']=0.
    stDB['H/C.1']=0.
    stDB['fracH']=0.
    stDB['Name']=""
    stDB['JREClass']=""
    stDB['SNOC']=0
    for Row in stDB.itertuples(index=False): #The index is still the row #, so loc(i) should pull the rows sequeentially.
        i+=1
        pauseCheck()
        #if i == 111 : break #JRE comment this to run the whole file.
        idDipprDB,tKelvinDB,STmN_mDB=stDB.loc[i][['ChemID','Temperature (K)','Value (N/m)']]
        if(idDipprDB>9999):idDipprDB -= 20000
        STmN_mDB*=1000
        stCalc=-8686
        try:
            TcK,PcMPa,acen,Mw=ULdf.loc[idDipprDB][['Tc','PcMPa','acen','Mw']]
        except:
            outfile.write(f"{idDipprDB} \t ULdb \t omits this compound.\n")
            print(f"{idDipprDB} \t Compound \t not found in UL dbase.")
            stDB.loc[i,'stCalc']= -8686 #JRE: signals omitted compound
            continue #JRE loop around if we can't find that compound.
        Tr=tKelvinDB/TcK
        P0=PcMPa*10**( 7/3*(1+acen)*(1-1/Tr) )*1E6
        if(idDipprDB != idOld): 
            curComp=ESD.EsdComp(idDipprDB)
            print(idDipprDB,curComp.Name,curComp.JREClass,curComp.atomCounts)
            idOld=idDipprDB
            print("CHBCFINOSS=",curComp.CHBCFINOSS)
        PsatEsd, vLEsd, vVEsd = curComp.psat(T=tKelvinDB, P0=P0) # GCM: omit P0 to use internal guess.
        if(PsatEsd==8686):
            outfile.write(f"{idDipprDB} \t Psat \t failed at T(K)={tKelvinDB}. Tr={Tr}.\n")
            print(f"{idDipprDB}  \t Psat \t failed at T(K)={tKelvinDB}. Tr={Tr} ")
            stDB.loc[i,'stCalc']= -86 #JRE: signals failed Psat
            continue #this cycles to next datapoint.
        PsatEsd/=1000
        rhoLEsd=curComp.Mw/vLEsd/1E6 # convert from mol/m3 to g/cc
        rhoVEsd=curComp.Mw/vVEsd/1E6
        acenCalc= -np.log10(PsatEsd/(1000*PcMPa))-1
        tension = sgt_pure(1/vVEsd,1/vLEsd,tKelvinDB,PsatEsd*1000,curComp)
        if(tension>0):stCalc=tension/1000
        stDB.loc[i,'stCalc']=stCalc
        stDB.loc[i,'Tr']=Tr
        stDB.loc[i,'bVol']=curComp.bVolCC_mol
        stDB.loc[i,'eps_kB']=curComp.eps_kB
        stDB.loc[i,'qShape']=curComp.qShape
        stDB.loc[i,'bondVol']=curComp.assoc.bondVolNm3[0]
        stDB.loc[i,'epsD_kB']=curComp.assoc.epsD_kB[0]
        stDB.loc[i,'epsA_kB']=curComp.assoc.epsA_kB[0]
        stDB.loc[i,'solPrm']=curComp.solp
        stDB.loc[i,'Name']=curComp.Name
        stDB.loc[i,'JREClass']=curComp.JREClass
        atomsTot=curComp.CHBCFINOSS.sum()
        nCarb=(curComp.CHBCFINOSS[0])
        nHy=(curComp.CHBCFINOSS[1])
        nSilox=int(  np.sqrt( curComp.CHBCFINOSS[9]*curComp.CHBCFINOSS[7] )  ) #if nOxy< 3*nSi, the int() should bring it down to nSi
        nSilox_C1=nSilox/(nCarb+1)
        nSilane=curComp.CHBCFINOSS[9]-nSilox
        if(nSilane<0):nSilane=0
        nIodo=curComp.CHBCFINOSS[5]
        nFluoro=(curComp.CHBCFINOSS[4])
        nChloro=(curComp.CHBCFINOSS[3])
        nSulf=(curComp.CHBCFINOSS[8])
        nNitro=curComp.CHBCFINOSS[6]
        nOxy=curComp.CHBCFINOSS[7]
        nOx_C1=nOxy/(nCarb+1)
        nHalo=curComp.CHBCFINOSS[2]+curComp.CHBCFINOSS[3]+curComp.CHBCFINOSS[4]+curComp.CHBCFINOSS[5] #BrClFI
        nHetero_C1=(nHalo+nSilox+nNitro+nOxy)/(nCarb+.1) #Sulfur ~carb, but keep separate
        nHalo_C1=(nHalo)/(nCarb+.1) 
        nF_C1=nFluoro/(nCarb+.1)
        H_C1=nHy/(nCarb+nSilane+.1) #Assume silanes behave like their HC hmomorph.
        fracH=nHy/atomsTot
        SNOCfactor=0
        if(nChloro>.95 or nSilox_C1>.95 or nNitro>.95 or nOx_C1 > 0.45 or nIodo>.95):SNOCfactor=1
        if(nChloro>1.95 or nSilox_C1>1.95 or nNitro>1.95 or nOx_C1 > 0.85):SNOCfactor=2
        if(nSilane>0):SNOCfactor=0
        if(nF_C1>0.4):SNOCfactor=0
        #nHetero_C1=(atomsTot-nCarb-nHy-nSulf*0.99-nSilane*0.98)/(nCarb+.1) #Sulfur ~carb, but keep separate
        stDB.loc[i,'SNOC']=SNOCfactor
        stDB.loc[i,'H/C.1']=H_C1
        stDB.loc[i,'fracH']=fracH
        stDB.loc[i,'nF_C.1']=nF_C1
        stDB.loc[i,'nHalo_C.1']=nHalo_C1
        stDB.loc[i,'nHetero_C.1']=nHetero_C1
        if(stDB.loc[i,'Class']=='normal' or stDB.loc[i,'Class']=='heavy' ):
            if(nHetero_C1<0.15):stDB.loc[i,'Class']='zrmHC'
            #elif(curComp.JREClass=='halom'):stDB.loc[i,'Class']='zrmHC'
            else:stDB.loc[i,'Class']='zrmHt'
        if(tension>0):stDB.loc[i,'pctDev']=(tension-STmN_mDB)/(STmN_mDB+0.5)*100 #0.5 offset reduces outliers in critical region
        #if(tension>0):stDB.loc[i,'pctDev']=(tension-STmN_mDB)/(STmN_mDB)*100 #0.5 offset reduces outliers in critical region
        print('Tr,expt,calc=',"{:10.5f},{:10.4f},{:10.4f}".format(Tr,STmN_mDB,tension))
        
    #    print( '%5d,%5.4f,%5.4f,%5.2f,%5.3f,%5.5f,%5.5f' % 
    #          (idDippr,acen,acenCalc,TKelvin,PsatEsd,rhoLEsd,rhoVEsd) )
        outfile.write( '%5d \t %7.2f \t %7.3f \t %5.5f \t %5.5f \t %10.4f \t %10.4f\t %10.4f \n' % 
              (idDipprDB,tKelvinDB,PsatEsd,rhoLEsd,rhoVEsd,STmN_mDB,tension,Tr) )
    #End of loop of DB
#End of with open()
stDB = stDB[stDB['qShape'].notna() & (stDB['qShape'] > 0)]
stDB = stDB[stDB['stCalc'].notna() & (stDB['stCalc'] > 0)]
stDB = stDB.sort_values(by='Class',ascending=False)
aadHC = stDB.loc[stDB['Class'] == 'zrmHC', 'pctDev'].abs().mean().round(2)
aadHt = stDB.loc[stDB['Class'] == 'zrmHt', 'pctDev'].abs().mean().round(2)
aadPo = stDB.loc[stDB['Class'] == 'polar', 'pctDev'].abs().mean().round(2)
aadAs = stDB.loc[stDB['Class'] == 'assoc', 'pctDev'].abs().mean().round(2)
print("     HC      Hetero       polar      assoc")
print(aadHC,aadHt,aadPo,aadAs)
#stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
stDB = stDB[stDB['stCalc'].notna() & (stDB['stCalc'] > 0)]
stDB.to_excel('stDBwCalcs.xlsx')
classes=stDB['Class'].unique() #Use unique() to compile the unique class names (ie. norml,heavy,polar,assoc)
#UresAB={} #dictionary to store slopes and intercepts
#from scipy.stats import linregress
# Loop through each density group
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharey=True)  # 2 row, 2 columns of plots
for ax, clas in zip(axes.flatten(), classes):
    class_data=stDB[stDB['Class']==clas]
    color='nHetero_C.1'
    if(clas=='assoc' or clas=='zrmHy'):color='H/C.1'
    scatter = ax.scatter(
        class_data['Tr'], class_data['pctDev'], c=class_data[color], cmap='viridis',alpha=0.5
    )
    plt.ylim(-20,20)
    ax.set_title(f'Class {clas}')
    ax.set_xlabel('Tr')
    ax.set_ylabel('pctDev')
fig.colorbar(scatter, ax=axes, orientation='vertical', label='NonC/C')
#plt.tight_layout()
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.1)
plt.show()