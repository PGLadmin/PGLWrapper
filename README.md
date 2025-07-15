# The Properties of Gases and Liquids, 6ed (PGL6ed)

## J. Richard Elliott, Vladimir Diky, Thomas Knotts IV, W. Vincent Wilding

Welcome to the PGL6ed GitHub Landing Page! Many of the methods in PGL6ed are sufficiently complicated that software support can be very helpful for readers to implement. Furthermore, we would like to ensure that our evaluations are consistent with the models as conceived by their authors and reproducible by any knowledgeable reader of PGL6ed. Therefore, we provide resources here to aid in those applications. Note that other resources are available at [PGL6ed.byu.edu](http://PGL6ed.byu.edu), primarily tables of group contribution constants for ideal gas properties pertaining to Chapter 4. Please note that all resources are provided here in good spirit but no guarantees of accuracy are intended or implied. If readers encounter errors, [please open a new issue](issues/new). The codes provided are written in Fortran, usually conforming to F90 standards but occasionally implementing Fortran77.

Keywords: Thermodynamics, Equations of State, COSMO-RS, COSMO-SAC, Local Coupled Cluster (LCC) vs. Group Contribution Estimation of Ideal Gas Formation, Critical Properties, Vapor Pressure Correlation and Prediction, Fluid Density, Virial Coefficients, Thermal Properties (H, Cp, compressibility), Vapor-Liquid and Liquid-Liquid and Solid-Liquid Phase Equilibrium (VLE, LLE, SLE), Quasichemical Theory, Infinite Dilution Activity Coefficients (IDACs), Viscosity, Thermal Conductivity, Diffusion, Surface Tension. 

## Resources

The resources on this GitHub site relate primarily to Chapters 6-9 of PGL6ed, focusing on equations of state (EOSs) to correlate and predict thermodynamic properties of pure compounds and mixtures. Chapters 6 and 7 are about pure compounds and mixture properties like vapor pressure, density, enthalpy, and entropy. Chapters 8 and 9 are about vapor-liquid equilibria (VLE), liquid-liquid equilibria (LLE), solid-liquid equilibria (SLE), and infinite dilutions activity coefficients (IDACs). EOSs, including reasonably sophisticated descriptions of mixture interactions, can provide consistent and accurate characterizations of all these properties. Activity models are closely related and generally derive from the same physical conceptions as EOS descriptions, occasionally with additional parameters interjected as in local composition models. The activity models are simple enough to implement without substantial help but we include code here for the NRTL model as an example. Most codes relate to implementations of EOS models. All the resources are organized with PGLWrapper as the main directory.

## Utility files

The [SourceCodeUA](SourceCodeUA) directory also includes many utility files that readers may find useful. These begin with a module file that illustrates how general constants and parameters can be conveyed to the various parts of the code that require them. Fortran programmers unfamiliar with modules can consider these as alternatives to Common statements. The advantage over Common statements is that changing a module does not require changing the rest of the code. For example, USE GlobConst appears in many subroutines because it conveys global constants like the ideal gas constant. If you wanted to add a new constant like RgasCal_mol to complement Rgas (default in Joules/mol-K=MPa-cm3/mol-K), you just need to add RgasCal_mol=Rgas*4.184 to the parameter list in Module GlobConst. The other files in [SourceCodeUA](SourceCodeUA) have filenames that reflect their content. Most of these include literature references at their beginning to facilitate reading about where the equations come from. Relevant literature and comprehensive descriptions are given in PGL6ed of course.

## Model parameters and experimental datasets

The [Input](Input) directory is also important because it provides resources for evaluating the various models as well as input for running the models. For example, [VLEJaubert.txt](Input/VLEJaubert.txt) is our implementation of Jaubert’s VLE database for non-aqueous binary mixtures. [VLEDannerGess.txt](Input/VLEDannerGess.txt) is a similar file provided with permission from Prof. R.P. Danner. [LLeDbPGL6ed96b.txt](Input/LLeDbPGL6ed96b.txt), [SLeDbPGL6ed20220129.txt](Input/SLeDbPGL6ed20220129.txt) are LLE and SLE databases provided by the NIST/TRC group in Boulder, Colorado. The [IdacRecLazzaroniDb.txt](Input/IdacRecLazzaroniDb.txt), [IdacRecJaubert+TdeOrginW.txt](Input/IdacRecJaubert+TdeOrginW.txt), and [IdacRecJaubert+TdeWinOrg.txt](Input/IdacRecJaubert+TdeWinOrg.txt) files provide our databases of IDACs, the first of which is an adaptation from Lazzaroni et al. The water in organic solvent (_WinOrg) file and organic solute in water (_OrginW) files are adaptations of databases of aqueous IDACs as described in Chapter 9.

## Our Mission

We do appreciate that there is much room for improvement in the quality of these codes, extensions to other languages (e.g., Python), and the nature of our organization. We encourage readers to make those improvements and we look forward to learning about that progress. We simply provide what we have as examples of one way to achieve the desired outcomes. In this way, we hope to contribute to our stated mission:
> "to survey the field of physical property estimation in order to promote the best approaches and objectively characterize the quality of the available methods, such that current applications are reliable and future progress of the field is advanced."

## Sample Uses:
One simple use case is to run the PGLEos.exe to generate phase diagrams, regress binary interaction parameters (BIPs), or reproduce the evaluations in PGL6ed. We provide a few examples below.
1. Ethanol+water at 0.1 MPa using ESD-MEM2
     This is a classic system going back to the earliest fermentations and distillations known to man. This example applies the ESD-MEM2 model (iEosOpt=4). The number of components is 2. The (DIPPR) ID numbers are 1102 and 1921. These can be found in .\input\ParmsPrTcJaubert.txt along with the names of the compounds and other key metrics. The screen will display the properties for this model, ending with the default BIP. The default BIP is based on optimization of VLE using the VLeJaubert.txt file as discussed in PGL6ed. You can select Y, to change it or you can change it later using the KN option. (Many common BIPs are called "kij" so KN stands for "kij new". Several other "K_" options pertain to changing BIP values.) For this example, simply hit "enter" to accept the default value of 0.0247. The screen displays a long list of possible calculations. Use the scroll bar to see them all. For this example, select TX and enter the pressure as 0.1. Default values of ethanol mole fractions are displayed with computed temperatures etc. A few lines of sample output are quoted below.
   T(K)      P(MPA)       x         y       etaL      etaV
  371.78    0.100000   0.00100   0.00946   0.37901   0.00032     0
  351.23    0.100000   0.95000   0.95417   0.33750   0.00087     0
  351.14    0.100000   0.99900   0.99906   0.33629   0.00090     0
 OutFile=C:\PGLWrapperRun\output\TPXY.txt
These lines show that the boiling temperatures estimated by the ESD-MEM2 model deviate slightly from the experimental values of 373.1 and 351.4. The result at xE=0.95, yE=0.954 is close to the computed azeotrope. Note that the experimental value for the azeotropic composition is xE=0.895, so this model does not match the azeotrope. (The default BIP minimizes bubble pressure deviations for a database of about 200 points, not just the azeotrope.) The eta values correspond to packing fractions of liquid and vapor. 

2. 1-butanol+water using PC-SAFT (iEosOpt=10)
     This system is interesting because it illustrates LLE as well as VLE. The IDs are: 1105 and 1921. You can use the TX option to check the boiling points of this model. PC-SAFT matches the boiling point of water more accurately than ESD (372.4). Accept the default BIP of 0.0, then select the BI option to generate the (LLE) binodal compositions starting from 270 K and a pressure of 0.4 MPa. The computations end at about 404 K with critical compositions near 0.29. The experimental values are roughly 398 K at x1=0.12, but 404 is close considering the BIP is set to zero. Copy the screen results or open the indicated txt file and paste into your favorite graphing software. Next select the TX option at 0.4 MPa to generate the VLE portion of the diagram. This reaches a minimum at about 408.3 K and x1=0.3. Copy and paste the results into your graphing software to generate the full phase diagram. The pressure of 0.4 MPa was chosen so the VLE and LLE do not collide. At 0.1 MPa, you would need to be careful about interpreting the exact points where the VLE and LLE intersect. 

3. 1-butanol+water using SPEADMD (iEosOpt=5)
    In this example, we vary the BIP to match the experimental value of the critical temperature. The default BIP is -0.0553. Running the BI option at 0.4 MPa and starting at 250 K gives a value of about 415 K. Use the KN option to change the BIP to -0.061, 0. The second value (0 in this case) assigns an optional temperature dependence according to kij=kij0+kij1/T. Now the critical temperature is about 403 K with a composition of x1=0.125. That's pretty close to the experimental value. 

## Sample Projects:

There are three projects that can be built easily using the codes from the PGLAdmin GitHub: PGLDLL, PGLDLLTest, and PGLEos. PGLDLL and its …Test function form the most expected use case. We presume that most readers will have their own flash codes and regression codes for optimizing BIPs and similar parameters. For them, we hope to facilitate incorporating the essential codes from PGL. We have had success calling PGLDLL from the C++ interface to TDE. The plan for the TDEFree version is to include PGLDLL along with other thermo models.  A DLL might not be the easiest way to implement these codes, however. Readers may prefer to simply add these .f90 files into their existing projects. That should be easy, but don’t forget to carry along the PGLWrapper\input directory where the parameter files for the various EOSs are located. We call the main directory PGLWrapper because it "wraps" all the PGL thermo models into a single entity that can be called either through the DLL or through a console application. 

Support for the PGLEos project has been added as an afterthought. We are happy to share these codes, but there are probably many superior codes available. The primary additions from PGLDLL to PGLEos are the MainMenu, MenuItems, FlashSubs, and MinPack. PGLDLL.f90 and PGLWrapper.f90 should not be included in the PGLEos project. The MainMenu provides a crude console interface to a long list of computational options. Most of these are common bubble, dew, and flash routines for VLE, LLE, SLE, and SFE. A few options are available for critical points and parameter regression. We use the LmDif code from MinPack to perform most optimizations. The PGLEos project was used for evaluations of VLE, LLE, IDACs, and SLE in PGL6ed. The databases for these evaluations are included in the PGLWrapper\input directory. The PGLDLL was used for evaluations in Chapter 6 and Chapter 7 through the TDE interface. Those evaluation databases are not available. Readers are welcome to use the PGLEos codes, but please do not expect support when, for example, a difficult flash in the critical region fails to converge. Readers are encouraged to develop their own codes for their own applications. We welcome learning about those developments. 

The following suggestions were developed for a 2025 installation of Windows Visual Studio Community 2019 using the 2025 Fortran Essentials (IFX) without the Math Kernel. This configuration minimizes the impact of the Fortran installation on hard drive space (~2GB). 

To build PGLDLL, add the following as "existing” files:
SourceCodeUA:
0ModuleUAData
FuEsd23a
FuFloryWert2a
FugiLSG-MEM2
FugiNRTL
FugiPr2a
FugiPrTc
fugiPRtcWS
fugiPRWS1a
FuSpeadmd23a
MEM2
PGLDLL.f90
PGLWrapper.f90
PGLFugiDefs
UaDataTools
UaEosTools
Wertheim23a

SourceCode2a:
PcSaftEosGross
PrLorraine

To build the PGLDLLTest project, add the following as "existing” files:
0ModuleUAData
DLLTestMain
(You will need to paste PGLDLL.DLL and PGDLL.lib into the dir where \Debug\PGLDLLTest.exe is built. You may also need to specify PGLDLL.lib in the “Additional Dependencies” of PGLDLLTest(IFX) > properties > Linker > Input.)

To build the PGLEos project, add the following as "existing” files:
SourceCodeUA:
0ModuleUAData
FuEsd23a
FuFloryWert2a
FugiLSG-MEM2
FugiNRTL
FugiPr2a
FugiPrTc
fugiPRtcWS
fugiPRWS1a
FuSpeadmd23a
MainMenu2d
MEM2
MenuItems2c
PGLFugiDefs
UaDataTools
UaEosTools
Wertheim23a

SourceCode2a:
PcSaftEosGross
PrLorraine

FlashSubs:
CR2e
Flashsubs1b

MinPack:
DEIGSRT
DGAUSSJ
DPYTHAG
DTQLI
DTRED2
INDEXX
LmdifEzCov2
MinTools




