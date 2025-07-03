# NRTL-PA Conversion Project

## Organization
Files are organized in folders according to their use:
- GMU.f90* - Main Aspen interface
- gammaModels
>- phys - physical models
>- nagata1.f90 - variation of Wilson equation by Nagata
>- nrtl.f90 - NRTL
>- scathild.f90 - Scatchard-Hildebrand
>- wilson.f90 - Wilson
- comb - combinatorial
>- flory.f90 - Flory
>- stavgugg.f90 - Staverman-Gugenheim modification
>- correqn1262.f90 - vdw covolume correction
- werth
>- calc_gammaw.f90 - main association routine
>- calc_dAdnk.f90 - calculations the derivatives of Helmholtz wrt n_k
- util - utilities
>- loadsites.f90* - routine to load sites from Aspen interface
>- vlu*.f90 - routine to calculate mixture molar volume assuming ideal mixing using Aspen pure component parameters.
>- sitenspecies.f90 - module to hold 'site' and 'species' information.
- commoncalc - routines used for all Wertheim methods
>- calc_gterms_PCSAFT.f90 - calculate derivative of g
>- calcX.f90* - calculate X for all sites
>- gcalc.f90 - calculate g for all species, calls calc_gterms_PCSAFT if needed
## Conversion
Files marked with * need conversion. See Conversion.md for details on the necessary conversion and progress of the conversion.