# Conversion Notes
## Current Effort
- Carl is working on input file formats and how to read input.
## GMU.f90
### Purpose
GMU.f90 is the main gamma routine that calculates the overall activity coefficient. GMU loads parameters from Aspen, calls routines to load sites, calculates molar volumes, the contributions of physical, wertheim, combinatorial, combinatorial correction and adds them up. Some routines include rudimentary 'print' statments used for debugging within the routines, but not provided normally to the Aspen history file. Only the higher level debugging is provided to the Aspen history file.
#### Variable input/output
- T, P, X, N, IDX - temp, press, mole fractions, number of species, idx is the Aspen indices for the species present in the function call because Aspen calls the routine with only the components in a stream. A simulation might have 10 components, but if only 4 are present in the stream, then the routine is called with N=4 and only those mole fractions.
- KCALC - integer to specify calculation: 1 - only gammas; 2 - only T derivative; 3 - gammas and T derivative
- COVOL - covolume used for calculating combintorial correction.
- KOP - A 10 element integer vector used to determine the options for calculating the activity coefficients
- GAMMA, DGAMMA - the NATURAL log of gamma and d(ln(gamma)/dT)
- NDS, KDIAG, KER, NSUB, NSUP, IDXSUB, IDXSUP, GAMSC, DGAMSC - not used.
##### General comments on I/O
- Your “Get__” functions should be included in your NRTLPA.f90 file for reading the parameters. FuEsd23.f90 or FuSpeadMd23a.f90 provide examples of reading association parameters.
- Probably, sourceCodeUA/FuSpeadMd23a will be more interesting for you. Around line 268 of that code, you will see how site-based parameters are read and stored. Line 370 calls GetAssocBips. FYI: BIPs are “binary interaction parameters.” GetAssocBips() is defined in sourceCodeUA/Wertheim.f90 (wertheim23.f90?) because Wertheim.f90 can be called generically for any TPT1, specifically ESD and SPEADMD in PGLWrapper.
- input/BipDA.txt file lists the BIPs between donors and acceptors. (Acceptors and donors have the same BIPs.) input/ParmsHb4.txt defines the volumes and energies of the donors and acceptors. input/SiteParms2580.txt defines the site types (also defined in ParmsHb4, but maybe it’s better to see the big picture).
- It might be nice if your definitions of site types could be the same as SPEADMD, just to have one less degree of confusion, but I don’t feel strongly about it. If you have a good reason to add a siteType that is not already available in SPEADMD, let me know. I will probably want to add it.
- For example, it would be logical (to me) if you referred to your list of bonding volumes and energies as ParmsHbNRTLPA.txt. BIPs for the physical interactions of NRTL could be called BipNRTLPA.txt. For cross-association BIPs, the filename BipDaNRTLPA.txt makes sense to me, but you may prefer BipDA_NRTLPA.txt.

#### Variables accessed from Aspen memory
- UNIFAC groups are used to setup the sites and association parameters. Since not all species are in a function call, then not all sites are in a fuction call. Use of UNIFAC groups is a limitation of the Aspen GUI.
- UNIQUAC R, Q - used for combinatorial correction.
- nsites - number of sites in Aspen case file
- nsitesp - number of sites present in function call - only the sites on the N species of the function call.
- KAD - kappa value for the association for the sites present, as an nsitesp x nsitesp matrix of interactions. The values are 0 for no interaction.
- EPS - epsilon/k value for the association sites present, as an nsitesp x nsitesp matrix of interactions. The values are 0 for no interaction.
- GMNRTL - user parameter for NRTL in the form tau_ij = a_ij + b_ij/T and alpha_ij is the third element.
#### Debugging output
- Debugging in Aspen is provided in a history file indentified as unit "global_nh". Aspen debugging output is provided if kop(2) > 0.
#### GMU Conversion
Much of the code can be deleted.
- The function DMS_IFCMNC looks up the memory location of Aspen variables in the vectors 'B' and 'IB' depending on whether the variable is real or integer and stores them in the vectors that are implemented. The bulk of the work will be
    - Eliminate Aspen look-ups
    - Replace variable assignments with those for PGLWrapper
    - Convert statements that write to Aspen history files with similar debugging.
## util/loadsites.f90
### Purpose
- Reads Aspen memory to find site present in the Aspen case file and in the function call. The routine loads the structures 'sitep' (type siteinfo) and 'compp' (type species).
- 'sitep' for each site holds the index to the site name and id, the species host index and host mole fraction, the number of occurences on the host, the 'packed id' for the smaller structure with only the sites present.
- 'compp' for each species holds the mole fraction, the number of sites on the species, the site ids for the sites, a flag to indicate if the species self associates.
- aspmx - 0 if a mixture has no self- or cross-association, 1 if association is present.
#### Debugging output
- Debugging in Aspen is provided in a history file indentified as unit "global_nh". Aspen debugging output is provided if kop(2) > 0.
#### loadsites Conversion
- The function DMS_IFCMNC looks up the memory location of Aspen variables in the vectors 'B' and 'IB' depending on whether the variable is real or integer and stores them in the vectors that are implemented. The bulk of the work will be
    - Eliminate Aspen look-ups
    - Replace variable assignments with those for PGLWrapper. The loading of the KAD and EPS matrix can be simplified for PGLWrapper since we have full access to how variables are loaded into memory.
    - Convert statements that write to Aspen history files with similar debugging.
## util/vlu.f90
### Purpose
- Calculate the density of pure components and the density of the mixture assuming ideal mixing
#### Debugging output
- Writes to control panel if Aspen equation is not programmed
#### vlu conversion
- Replace with Rackett/Thodos-Campbell. We will need parameters.
## commoncalc/CalcX.f90
### Purpose
- Calculates X given site KAD, EPS and host mole fractions and mixture density.
#### Debugging output
- Output is provided to the history file and user control panel if the iteration does not converge.
#### CalcX conversion
- Convert statements that write to Aspen history files with similar debugging.