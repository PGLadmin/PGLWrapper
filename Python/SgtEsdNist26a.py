# -*- coding: utf-8 -*-
"""
Purpose: Apply the SGTPy code using the ESD EOS. 
References: Mejia, Muller, Chaparro, J. Chem. Inf. Model. 2021, 61, 1244−1250.
            Elliott, Suresh, Dohohue IECR, 1990.
"""
import sys
print("DEBUG:", sys.executable)

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
import ESD  #JRE: absolute import. Requires ESD.py in working dir.
import pandas as pd
from GlobConst import PGLInputDir,DownloadDir,pauseCheck

#"""
#import scipy as sp
#import Cython as cp
#import sgtpy as saft
#"""
testComp=ESD.EsdComp(2684) #1997 is deuterium oxide.
atomsTot=testComp.CHBCFINOSS.sum
nCarb=testComp.CHBCFINOSS[0]
nFluoro=testComp.CHBCFINOSS[4]
F_C1=nFluoro/(nCarb+1) #add one in denom to avoid divide by zero.
print("nFluoro/(nCarb+1)=",F_C1)
#pauseCheck(1)

print("Loading ULdf and IDdb....")
ULdf=pd.read_csv(PGLInputDir+'ParmsPrTcJaubert.txt',sep=r'[\t,]', skipinitialspace=True,header=0,engine='python')
ULdf.columns = ULdf.columns.str.strip() # Remove leading and trailing whitespace from column names
print("ULdf",ULdf.columns)
ULdf.set_index('idDippr',inplace=True)
IDdb=pd.read_csv(PGLInputDir+'IdTrcDipCas.txt',sep=r'[\t]', skipinitialspace=True,header=0,engine='python')
IDdb.columns = IDdb.columns.str.strip() # Remove leading and trailing whitespace from column names

#PPData = ULdf.set_index(ULdf.columns[0]).T.to_dict()
#del ULdf # free up the memory of the ULdf dataframe.
print("Loading stDB. If full dataset be patient...")

# DataBase='SurfTenDbNist26a.xlsx'
# stDB=pd.read_excel(DownloadDir+DataBase) #,nrows=44)
DbSource="DIPPR"
DbSource="NIST"
DataBase='SurfTenDbDIPPR26a.xlsx'
if(DbSource=="NIST"):DataBase='ST-public.txt'
stDB=pd.read_csv(DownloadDir+DataBase,sep=r'[\t]', skipinitialspace=True,header=0,engine='python') 
print(stDB.head(10))
print(stDB.columns)
print("Cleaning stDB....")
valid_dippr = set(ULdf.index)           # Remove rows from stDB that don't have a valid idDippr in ULdf
valid_trc = IDdb[IDdb["idDippr"].isin(valid_dippr)]["idTrc"]
valid_trc = set(valid_trc)
stDB = stDB[stDB["idTrc"].isin(valid_trc)]
stDB["idDippr"] = stDB["idTrc"].map(    # Add the idDippr column to stDB by mapping from IDdb
    IDdb.set_index("idTrc")["idDippr"] )
stDB["JREClass"] = stDB["idDippr"].map(ULdf["Class"])  # Add the JREclass column to stDB by mapping from ULdf
stDB["Name"] = stDB["idDippr"].map(ULdf["Name"])
stDB["Tr"] = stDB["T/K"] / ULdf["Tc"]
TrCut=0.960
#stDB = stDB[stDB['Tr'] < TrCut]         # Uncertainty too high when Tr > TrCut.   
stDB.reset_index(drop=True, inplace=True)   #this is necessary for read_csv and OK otherwise. 
print(f"stDB:len={len(stDB)}")
print(f"ULdf:len={len(ULdf)}")
print(f"stDB cleaned. Tr<{TrCut}, valid idDippr and idTrc. Press enter to start calculations.")
pauseCheck(1)

outfile=open('EsdOut.txt',"w")
outfile.close()
with open('EsdOut.txt',"a") as outfile: #The "with" will close the file properly even if there is a crash.
    #print("     T(K)       Psat(kPa)       rhoL(g/cc)       rhoV(g/cc)")
    print("idDippr,acen,acenCalc,T(K),Psat(kPa),rhoL(g/cc),rhoV(g/cc),stCalc,")
    outfile.write("idDippr \t T(K) \t Psat(kPa) \t rhoL(g/cc) \t rhoV(g/cc) \t stCalc \t st(mN/m) \t    Tr  \n")
    i=-1    # initialize negative so +1 will make the first element "0"
    idOld=-1 # initialize to impossible value for later reference.
    numeric_cols = [
    'stCalc', 'Tr', '1-Tr', 'pctDev', 'bVol', 'eps_kB', 'qShape',
    'bondVol', 'epsD_kB', 'epsA_kB', 'nHetero_C.1', 'nHalo_C.1',
    'nChlo_C.1', 'nF_C.1', 'H/C.1', 'fracH', 'SNOC'  ]
    stDB[numeric_cols] = 0.
    for Row in stDB.itertuples(index=False): #The index is still the row #, so loc(i) should pull the rows sequeentially.
        i+=1
        pauseCheck()
        #if i == 111 : break #JRE comment this to run the whole file.
        #idDipprDB,tKelvinDB,STmN_mDB=stDB.loc[i][['ChemID','Temperature (K)','Value (N/m)']]   #dippr style
        idDipprDB,tKelvinDB,STmN_mDB=stDB.iloc[i][['idDippr','T/K','STN_m']]
        if(idDipprDB>9999):idDipprDB -= 20000 #this corrects review compounds in the DIPPR DB. 
        STmN_mDB*=1000 # convert to mN/m. 
        stCalc=-8686
        try:
            TcK,PcMPa,acen,Mw=ULdf.loc[idDipprDB][['Tc','PcMPa','acen','Mw']]
        except:
            outfile.write(f"{idDipprDB} \t ULdb \t omits this compound.\n")
            print(f"{idDipprDB} \t Compound \t not found in ULorraine dbase.")
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
            stDB.loc[i,'stCalc']= -86       #JRE: signals failed Psat
            continue                        #this cycles to next datapoint.
        PsatEsd/=1000                       # convert from Pa to kPa.
        rhoLEsd=curComp.Mw/vLEsd/1E6        # convert from mol/m3 to g/cc
        rhoVEsd=curComp.Mw/vVEsd/1E6
        acenCalc= -np.log10(PsatEsd/(1000*PcMPa))-1
        tension = sgt_pure(1/vVEsd,1/vLEsd,tKelvinDB,PsatEsd*1000,curComp)
        if(tension>0):stCalc=tension/1000   #convert from Pa-s to mN/m
        stDB.loc[i,'stCalc']=stCalc
        #stDB.loc[i,'Tr']=Tr                #already calculated in the dataframe above.
        stDB.loc[i,'1-Tr']=1-Tr
        stDB.loc[i,'bVol']=curComp.bVolCC_mol
        stDB.loc[i,'eps_kB']=curComp.eps_kB
        stDB.loc[i,'qShape']=curComp.qShape
        stDB.loc[i,'bondVol']=curComp.assoc.bondVolNm3[0]
        stDB.loc[i,'epsD_kB']=curComp.assoc.epsD_kB[0]
        stDB.loc[i,'epsA_kB']=curComp.assoc.epsA_kB[0]
        stDB.loc[i,'solPrm']=curComp.solp
        atomsTot=curComp.CHBCFINOSS.sum()
        nCarb=(curComp.CHBCFINOSS[0])
        nHy=(curComp.CHBCFINOSS[1])
        nHetero=atomsTot-nCarb-nHy
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
        if(nCarb > 1):stDB.loc[i,'nF_C.1']=nF_C1
        stDB.loc[i,'nHalo_C.1']=nHalo_C1
        if(nCarb > 1):stDB.loc[i,'nChlo_C.1']=nChloro/(nCarb) # ignore inorganics 
        stDB.loc[i,'nHetero_C.1']=nHetero_C1
        UpperClass=stDB.loc[i,'JREClass'].upper()   
        stDB.loc[i,'Class']=UpperClass  #For NRMHC and NRMHT, Class is the same.
        if(UpperClass == "NRMSU"):stDB.loc[i,'Class']="NRMHT" #JRE: NRMSU is close to NRMHT.
        Upper2=UpperClass[:2]
        if(Upper2 != "NR"):             #NRMHC and NRMHT remain the same.
            class_map = { "PO":"POLAR","HE":"HEAVY","AS":"ASSOC","HA":"NRMHT","MU":"NRMHT","GA":"NRMHT","IN":"NRMHT","SI":"NRMHT","PR":"NRMHT","OR":"NRMHT" }
            #ToDo: HA=HALOM, PR=PRFLO
            stDB.loc[i, 'Class'] = class_map[Upper2]
        if(tension>0):stDB.loc[i,'pctDev']=(tension-STmN_mDB)/(STmN_mDB)*100 #0.5 offset reduces outliers in critical region
        #Using the calculated value in the denominator here facilitates correlation. Change to expt on final.
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
stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
stDB = stDB.sort_values(by='Class',ascending=False)

tbl = stDB.groupby('Class').agg( AAD=('pctDev', lambda x: x.abs().mean()), BIAS=('pctDev', 'mean'), Pts=('pctDev', 'count'), Comp=('Name', 'nunique') ).round(2).T
tbl.index = ['%AAD', 'BIAS', '#Pts', '#Comp'] 
tbl = tbl[['NRMHC', 'NRMHT', 'POLAR', 'ASSOC', 'HEAVY']]
print(tbl.to_string(
    formatters={ col: lambda x: f"{int(x)}" if x == int(x) else f"{x:.2f}"
    for col in tbl.columns }  )) 

print(PGLInputDir)
#stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
#stDB.to_excel('stDBwCalcs.xlsx')
classes=stDB['Class'].unique().tolist() #Use unique() to compile the unique class names (ie. norml,heavy,polar,assoc)
for class_name in classes:
    subset = stDB[stDB['Class'] == class_name]
    file_name = f'stDB{class_name}.xlsx'
    subset.to_excel(file_name, index=False)
#UresAB={} #dictionary to store slopes and intercepts
#from scipy.stats import linregress
# Loop through each density group
import matplotlib.pyplot as plt
fig, axes = plt.subplots(3, 2, figsize=(8, 8), sharey=True)  # 3 row, 2 columns of plots
vmin = stDB['bVol'].min()
vmax = stDB['bVol'].max()
for ax, clas in zip(axes.flatten(), classes):
    class_data=stDB[stDB['Class']==clas]
    color='bVol'
    #if(clas=='assoc' or clas=='zrmHy'):color='H/C.1'
    scatter = ax.scatter(
        class_data['1-Tr'], class_data['pctDev'], c=class_data[color], cmap='viridis',vmin=vmin, vmax=vmax, alpha=0.5
    )
    ax.set_ylim(-90,90)
    ax.set_title(f'Class {clas}')
    ax.set_xlabel('1-Tr')
    ax.set_ylabel('pctDev')
fig.tight_layout()
fig.colorbar(scatter, ax=axes, orientation='vertical', label='bVol (cm³/mol)')
#plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.4, hspace=0.2)
fig.savefig(r'C:\Users\ellio\OneDrive\Downloads\ClassPlots.svg',dpi=300,bbox_inches='tight')
plt.show()