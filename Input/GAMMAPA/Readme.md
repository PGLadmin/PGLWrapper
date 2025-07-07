# Instructions for Input Parameters
This folder holds:
- One parameter database files: 'Database - Pure and Site param.xlsx'
- Two parameter input templates with a subset of parameters for a given run
- A set of parameter files for different runs

## Configuring Input Parameter Files
Input parameter files must be tab delimited. This format makes them easy to read and avoids confusion with commas in literature references contained within the files.

The code reads two input files: one for the pure and site parameters and a second for the physical model parameters. Currently, only the NRTL physical parameter is coded.

To create input files, open the Excel template file and copy/paste values to replace those in the template. The original Excel templates can always be recovered from the repository if you make an error, or if you choose, you can create a copy before preparing your input files.

The format is similar for both input files. Each file has header rows that should never be removed as the code expects them to be present. The header rows for each section also include an integer in the first column to indicate how many rows are to be read from each section.

The TRC id numbers are used to identify and lookup species. Always update the correct TRC number if you change the names so that the two input file parameters can be related to each other.

## Preparing Main Input
The main input file includes
- The run options for the 'KOP' vector that specify how the parameters are to be applied. The explanation of KOP is in gammapa.f90.
- The component names and TRC id numbers the pure components for the run
- The standard API volumes (60 F) that can be used instead of temperature-dependent volumes
- The selected EOS covolume that is used for the rdf for the CPA, VDW, and ESD models. Enter a zero if you want to use a constant rdf. The covolume is a required parameter.
- The pure species PCSAFT parameters used for the rdf for that model. If the values are left blank and the correct KOP value is set, then the CPA, VDW, and ESD or constant rdf run options can still be used.
- The literature 'source' from the database for cross referencing of the input parameters.
- The liquid density parameters include the equation number to be used for the temperature-dependent volume parameters. These can be omitted if the API values are used with the correct KOP values.
- The site ids for sites on all species. The electron donor convention is used.
- The site bonding parameters that correspond to the selected rdf model. Note that values fitted to one rdf model should not be applied with a different rdf model especially the preexponential value.

Copy information from the database file and paste into the appropriate table. Any number of species can be used and any number of sites, though the limitations for sites on a given species are below.

When the input file is ready, use 'Save As...' and select a tab separated text file as the file type.

## Preparing Physical Model Input
The physical model input file includes
- The desired NRTL parameters from the database. We recommend always entering species in the database using the smaller id number for the first species to make searching easier. The order is not critical for the code. Recall that parameters regressed with one association rdf model should not be applied with another rdf association model.

When the input file is ready, use 'Save As...' and select a tab separated text file as the file type.

## Limitations on Sites
The current code assumes that a given site type occurs on only one species. However, all sites of a given type on a species share the same parameter values. For example, the -OH on a primary alcohol can have diffent values for the accepting and donating values on each primary alcohol present in a mixture. However, multiple sites of the same type on a given species share the same association values. For ethylene glycol, both alcohol groups share the same association parameter values for the oxygens and protons, though they may differ from the association parameter values on ethanol or 1-propanol.