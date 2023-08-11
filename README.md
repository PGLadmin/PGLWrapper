# The Properties of Gases and Liquids, 6ed (PGL6ed)

## J. Richard Elliott, Vladimir Diky, Thomas Knotts IV, W. Vincent Wilding

Welcome to the PGL6ed GitHub Landing Page! Many of the methods in PGL6ed are sufficiently complicated that software support can be very helpful for readers to implement. Furthermore, we would like to ensure that our evaluations are consistent with the models as conceived by their authors and reproducible by any knowledgeable reader of PGL6ed. Therefore, we provide resources here to aid in those applications. Note that other resources are available at [PGL6ed.byu.edu](http://PGL6ed.byu.edu), primarily tables of group contribution constants for ideal gas properties pertaining to Chapter 4. Please note that all resources are provided here in good spirit but no guarantees of accuracy are intended or implied. If readers encounter errors, we welcome feedback through the "Issues" link in the upper left of the landing page. The codes provided are written in Fortran, usually conforming to F90 standards but occasionally implementing Fortran77.

## Resources

The resources on this GitHub site relate primarily to Chapters 6-9 of PGL6ed, focusing on equations of state (EOSs) to correlate and predict thermodynamic properties of pure compounds and mixtures. Chapters 6 and 7 are about pure compounds and mixture properties like vapor pressure, density, enthalpy, and entropy. Chapters 8 and 9 are about vapor-liquid equilibria (VLE), liquid-liquid equilibria (LLE), solid-liquid equilibria (SLE), and infinite dilutions activity coefficients (IDACs). EOSs, including reasonably sophisticated descriptions of mixture interactions, can provide consistent and accurate characterizations of all these properties. Activity models are closely related and generally derive from the same physical conceptions as EOS descriptions, occasionally with additional parameters interjected as in local composition models. The activity models are simple enough to implement without substantial help but we include code here for the NRTL model as an example. Most codes relate to implementations of EOS models. All the resources are organized with PGLWrapper as the main directory.

## Utility files

The [SourceCodeUA](SourceCodeUA) directory also includes many utility files that readers may find useful. These begin with a module file that illustrates how general constants and parameters can be conveyed to the various parts of the code that require them. Fortran programmers unfamiliar with modules can consider these as alternatives to Common statements. The advantage over Common statements is that changing a module does not require changing the rest of the code. For example, USE GlobConst appears in many subroutines because it conveys global constants like the ideal gas constant. If you wanted to add a new constant like RgasCal_mol to complement Rgas (default in Joules/mol-K=MPa-cm3/mol-K), you just need to add RgasCal_mol=Rgas*4.184 to the parameter list in Module GlobConst. The other files in [SourceCodeUA](SourceCodeUA) have filenames that reflect their content. Most of these include literature references at their beginning to facilitate reading about where the equations come from. Relevant literature and comprehensive descriptions are given in PGL6ed of course.

## Model parameters and experimental datasets

The [Input](Input) directory is also important because it provides resources for evaluating the various models as well as input for running the models. For example, [VLEJaubert.txt](Input/VLEJaubert.txt) is our implementation of Jaubert’s VLE database for non-aqueous binary mixtures. [VLEDannerGess.txt](Input/VLEDannerGess.txt) is a similar file provided with permission from Prof. R.P. Danner. [LLeDbPGL6ed96b.txt](Input/LLeDbPGL6ed96b.txt), [SLeDbPGL6ed20220129.txt](Input/SLeDbPGL6ed20220129.txt) are LLE and SLE databases provided by the NIST/TRC group in Boulder, Colorado. The [IdacRecLazzaroniDb.txt](Input/IdacRecLazzaroniDb.txt), [IdacRecJaubert+TdeOrginW.txt](Input/IdacRecJaubert+TdeOrginW.txt), and [IdacRecJaubert+TdeWinOrg.txt](Input/IdacRecJaubert+TdeWinOrg.txt) files provide our databases of IDACs, the first of which is an adaptation from Lazzaroni et al. The water in organic solvent (_WinOrg) file and organic solute in water (_OrginW) files are adaptations of databases of aqueous IDACs as described in Chapter 9.

## Our Mission

We do appreciate that there is much room for improvement in the quality of these codes, extensions to other languages (e.g., Python), and the nature of our organization. We encourage readers to make those improvements and we look forward to learning about that progress. We simply provide what we have as examples of one way to achieve the desired outcomes. In this way, we hope to contribute to our stated mission:
> "to survey the field of physical property estimation in order to promote the best approaches and objectively characterize the quality of the available methods, such that current applications are reliable and future progress of the field is advanced."
