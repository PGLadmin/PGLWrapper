# -*- coding: utf-8 -*-
"""
Purpose: Apply the SGTPy code using the ESD EOS. 
References: Mejia, Muller, Chaparro, J. Chem. Inf. Model. 2021, 61, 1244−1250.
            Elliott, Suresh, Dohohue IECR, 29:1476-1485 (1990).
            Elliott, JCED, 71:1994-2005 (2026). doi.org/10.1021/acs.jced.5c00675
"""
import sys
print("DEBUG:", sys.executable)
import numpy as np
import os
print(f"JRE: Current Working Directory: {os.getcwd()}")
def BrockBird(TcK,PcMPa,Tb,tKelvin):
    """ Returns the surface tension in mN/m. 
        TcK = critical temperature in K
        PcMPa = critical pressure in MPa
        Tb = normal boiling point in K
        tKelvin = temperature of interest in K
        References: Brock and Bird, J. Chem. Eng. Data, 1962, 7, 233-237. PGL6ed 13.3-3 to 13.3-7.
    """ 
    if(tKelvin>=TcK):return 0.0 #JRE: Above the critical point, there is no surface tension.
    Tr=tKelvin/TcK
    alphac=0.9076*( 1-Tb*np.log(0.101325/PcMPa)/(TcK-Tb) )
    PcBar=PcMPa*10 # convert from MPa to bar
    surfTenDyn_cm=(0.132*alphac-0.279)*(1-Tr)**(11/9)*(TcK*PcBar*PcBar)**(1/3) #JRE: PGL6ed 13.3-3 to 13.3-7
    # TcK*PcMPa*Rg/Rg [=] (J-K/cm3); J-K/cm3 * PcMPa *Rg [=] 1E6 J^2/m^3 *1E-6 N/m2 [=] N^2-m^2/m^3*N/m^2 [=] N^3/m^3
    return surfTenDyn_cm
# Get the directory of the current script
#script_dir = os.path.dirname(os.path.abspath(__file__))
import ESD  #JRE: absolute import. Requires ESD.py in working dir.
import pandas as pd
from GlobConst import PGLInputDir,DownloadDir,pauseCheck,LOUD
print("PGLInputDir=",PGLInputDir)
print("DownloadDir=",DownloadDir)
ULdf=pd.read_csv(PGLInputDir+'ParmsPrTcJaubert.txt',sep=r'[\t,]', skipinitialspace=True,header=0,engine='python')
ULdf.columns = ULdf.columns.str.strip() # Remove leading and trailing whitespace from column names
print(ULdf.columns)
ULdf.set_index('idDippr',inplace=True)
IDdb=pd.read_csv(PGLInputDir+'IdTrcDip.txt',sep=r'[\t]', skipinitialspace=True,header=0,engine='python')
IDdb.columns = IDdb.columns.str.strip() # Remove leading and trailing whitespace from column names
print(  "IDdb duplicates:",             # Previous version of IdTrcDip had duplicate idTrc. 
    IDdb[IDdb["idTrc"].duplicated(keep=False)]
    .sort_values("idTrc")  )

print("Loading stDB. If full dataset be patient...")
DbSource="DIPPR"
#DbSource="NIST"
DataBase='SurfTenDatabase.xlsx'
# stDB=pd.read_excel(DownloadDir+DataBase) #,nrows=44)
if(DbSource=="NIST"):DataBase='ST-public.txt'
if(DbSource=="DIPPR"):DataBase='SurfTenDipprClean.txt'
stDB=pd.read_csv(DownloadDir+DataBase,sep=r'[\t]', skipinitialspace=True,header=0,engine='python') 
stDB.columns = stDB.columns.str.strip() # Remove leading and trailing whitespace from column names
print("stDB",stDB.columns)
#stDB.set_index('idTrc',inplace=True)
print("stDB",stDB.head())
print("stDB",stDB.columns)
if(DbSource=="NIST"):
    valid_dippr = set(ULdf.index)                       # Remove rows from stDB that don't have a valid idDippr in ULdf
    valid_trc = IDdb[IDdb["idDippr"].isin(valid_dippr)]["idTrc"]
    valid_trc = set(valid_trc)
    stDB = stDB[stDB["idTrc"].isin(valid_trc)]
    stDB["idDippr"] = stDB["idTrc"].map(                # Add the idDippr column to stDB by mapping from IDdb
        IDdb.set_index("idTrc")["idDippr"] )
else:
    stDB["idTrc"] = stDB["idDippr"].map(                # Add the idDippr column to stDB by mapping from IDdb
        IDdb.set_index("idDippr")["idTrc"] )
print(f"stDB1:len={len(stDB)}", stDB.head(1))
print(f"ULdf:len={len(ULdf)}", ULdf.head(1))
stDB["JREClass"] = stDB["idDippr"].map(ULdf["Class"])  # Add the JREclass column to stDB by mapping from ULdf
print(f"stDB1:len={len(stDB)}")
stDB["Name"] = stDB["idDippr"].map(ULdf["Name"])    # Add Name column
stDB["Tr"] = stDB["T/K"] / ULdf["Tc"]               # Add Tr column
print(f"stDB2:len={len(stDB)}")                      # Number of points with valid idDippr
TrCut=1.0
#stDB = stDB[stDB['Tr'] < TrCut]                    # Uncertainty too high when Tr > TrCut.  JRE20260923: Wiped out half of DB.
stDB.dropna(subset=["idDippr"], inplace=True)       # delete rows that don't have a Dippr ID. These are mostly inorganics, which we don't want to model.
print(f"stDB3:len={len(stDB)}")                      # This should be same as after dippr filter. 

dropped_rows = stDB[stDB["JREClass"].isna()]
if len(dropped_rows) > 0:
    unique_names = dropped_rows["idCas"].dropna().unique()
    print(f"Unique compds to be dropped: {len(unique_names)}")
    print(dropped_rows[["idDippr","idCas"]].value_counts())
stDB.dropna(subset=["JREClass"], inplace=True)      # delete rows that don't have a JREclass. e.g., ionic liquids.
#stDB.set_index("idDippr", inplace=True)            # change the index to the Dippr ID, which is the key for the ULdf dataframe.
stDB.reset_index(drop=True, inplace=True)           # reset the index to a default integer index
print(f"stDB4:len={len(stDB)}")                     # This should be same as after dippr filter. 
pauseCheck(1)
outfile=open('SurfTenOut.txt',"w")                  # Clear previous output file. 
outfile.close()
with open('SurfTenOut.txt',"a") as outfile:         # The "with" will close the file properly even if there is a crash.
    if(LOUD):print("idTrc \t T(K) \t stCalc \t st(mN/m) \t    Tr  \n")
    outfile.write("idTrc \t T(K) \t stCalc \t st(mN/m) \t    Tr  \n")
    i=-1    # initialize negative so +1 will make the first element "0"
    idOld=-1 # initialize to impossible value for later reference.
    numeric_cols = [ 'stCalc', 'Tr', '1-Tr', 'pctDev', 'bVol', 'eps_kB', 'qShape',
    'bondVol', 'epsD_kB', 'epsA_kB', 'nHetero_C.1', 'nHalo_C.1',
    'nChlo_C.1', 'nF_C.1', 'H/C.1', 'fracH', 'SNOC'  ]
    stDB[numeric_cols] = 0.
    stDB[ 'Class'] = ""                         # Class is a simplified version of JREClass, close to PGL6ed.
    for Row in stDB.itertuples(index=False):    # The index is still the row #, so loc(i) should pull the rows sequentially.
        i+=1
        pauseCheck()
        #if i == 111 : break #JRE comment this to run the whole file.
        idDipprDB,tKelvinDB,STmN_mDB=stDB.loc[i][['idDippr','T/K','STN_m']]
        if(idDipprDB>9999):idDipprDB -= 20000   # This corrects review compounds in the DIPPR DB. 
        STmN_mDB*=1000 # convert to mN/m. 
        stCalc= -8686
        try:
            TcK,PcMPa,acen,Tb,Mw=ULdf.loc[idDipprDB][['Tc','PcMPa','acen','Tb','Mw']]
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
            if(LOUD):print("CHBCFINOSS=",curComp.CHBCFINOSS)
        stCalc = BrockBird(TcK,PcMPa,Tb,tKelvinDB)  # BrockBird returns the surface tension in mN/m [=] Dyne/cm.
        stDB.loc[i,'stCalc']=stCalc
        stDB.loc[i,'Tr']=Tr                         # Done above but doesn't hurt.
        stDB.loc[i,'1-Tr']=1-Tr
        stDB.loc[i,'bVol']=curComp.bVolCC_mol
        stDB.loc[i,'eps_kB']=curComp.eps_kB
        stDB.loc[i,'qShape']=curComp.qShape
        stDB.loc[i,'bondVol']=curComp.assoc.bondVolNm3[0]
        stDB.loc[i,'epsD_kB']=curComp.assoc.epsD_kB[0]
        stDB.loc[i,'epsA_kB']=curComp.assoc.epsA_kB[0]
        stDB.loc[i,'solPrm']=curComp.solp
        #stDB.loc[i,'Name']=curComp.Name                # Already added from ULdb.
        #stDB.loc[i,'JREClass']=curComp.JREClass        # Already added from ULdb.
        atomsTot=curComp.CHBCFINOSS.sum()
        if(atomsTot<1):atomsTot=1 #JRE: avoid divide by zero errors. e.g., Can happen for Nobel gases.
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
        stDB.loc[i,'Class']=UpperClass
        if(UpperClass == "NRMSU"):stDB.loc[i,'Class']="NRMHT" #JRE: NRMSU is close to NRMHT.
        Upper2=UpperClass[:2]
        if(Upper2 != "NR"):             #NRMHC and NRMHT remain the same.
            class_map = { "PO":"POLAR","HE":"HEAVY","AS":"ASSOC","HA":"NRMHT","MU":"NRMHT","GA":"NRMHT","IN":"NRMHT","SI":"NRMHT","PR":"NRMHT","OR":"NRMHT" }
            #ToDo: HA=HALOM, PR=PRFLO
            stDB.loc[i, 'Class'] = class_map[Upper2]
        if(stCalc>0):stDB.loc[i,'pctDev']=(stCalc-STmN_mDB)/(STmN_mDB)*100 #0.5 offset reduces outliers in critical region
        #Using the calculated value in the denominator here facilitates correlation. Change to expt on final.
        if(LOUD):print('Tr,expt,calc=',"{:10.5f},{:10.4f},{:10.4f}".format(Tr,STmN_mDB,stCalc))
        
    #    print( '%5d,%5.4f,%5.4f,%5.2f,%5.3f,%5.5f,%5.5f' % 
    #          (idDippr,acen,acenCalc,TKelvin,PsatEsd,rhoLEsd,rhoVEsd) )
        outfile.write( '%5d \t %7.2f \t %7.3f \t %5.5f \t %5.5f \n' % #\t %10.4f \t %10.4f\t %10.4f \n' % 
              (idDipprDB,tKelvinDB,STmN_mDB,stCalc,Tr) )
    #End of loop over DB
#End of with open()
stDB = stDB[stDB['qShape'].notna() & (stDB['qShape'] > 0)]
stDB = stDB[stDB['stCalc'].notna() & (stDB['stCalc'] > 0)]
stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
stDB = stDB.sort_values(by='Class',ascending=False)

tbl = stDB.groupby('Class').agg( AAD=('pctDev', lambda x: x.abs().mean()), BIAS=('pctDev', 'mean'), Pts=('pctDev', 'count'), Comp=('Name', 'nunique') ).T
tbl.index = ['%AAD', 'BIAS', '#Pts', '#Comp'] 
#tbl.loc['#Pts'] = tbl.loc['#Pts'].astype(int)
#tbl.loc['#Comp'] = tbl.loc['#Comp'].astype(int)
#tbl.loc['%AAD'] = tbl.loc['%AAD'].round(2)
#tbl.loc['BIAS'] = tbl.loc['BIAS'].round(2)
tbl = tbl[['NRMHC', 'NRMHT', 'POLAR', 'ASSOC', 'HEAVY']]

print(  tbl.to_string(
    formatters={ col: lambda x: f"{int(x)}" if x == int(x) else f"{x:.2f}"
    for col in tbl.columns } )  ) 

print(PGLInputDir)
#stDB = stDB[stDB['pctDev'].notna() & (stDB['pctDev'] != 0)]
stDB = stDB[stDB['Tr'] < 0.960]
#stDB.to_excel('stDBwCalcs.xlsx')
classes=stDB['Class'].unique().tolist() #Use unique() to compile the unique class names (ie. norml,heavy,polar,assoc)
nClasses=len(classes)
for class_name in classes:
    subset = stDB[stDB['Class'] == class_name]
    file_name = f'stDB{class_name}.xlsx'
    subset.to_excel(file_name, index=False)
#UresAB={} #dictionary to store slopes and intercepts
#from scipy.stats import linregress
# Loop through each density group
import matplotlib.pyplot as plt
#fig, axes = plt.subplots(3, 2, figsize=(8, 8), sharey=True)  # 3 row, 2 columns of plots
fig, axes = plt.subplots(5, 1, figsize=(4, 13.5), sharey=True)
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
    ax.set_xlim(0,0.6)
    ax.set_title(f'Class {clas}')
    ax.set_xlabel('1-Tr')
    ax.set_ylabel('pctDev')
fig.tight_layout()
fig.colorbar(scatter, ax=axes, orientation='vertical', label='bVol (cm³/mol)')
#plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.4, hspace=0.2)
fig.savefig(r'C:\Users\ellio\OneDrive\Downloads\ClassPlots.svg',dpi=300,bbox_inches='tight')
plt.show()