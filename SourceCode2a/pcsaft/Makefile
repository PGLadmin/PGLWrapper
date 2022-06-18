 FC=  ifort
# FC=  gfortran

# ---- compiling options for ifort --------------  
 FFLAGS= -O3
# ---- ifort: debuging --------------------------
# FFLAGS = -O3 -warn all -check bounds -traceback

# ---- compiling options for gfortran -----------
# FFLAGS= -O3
# ---- gfortran: debugging ----------------------
# FFLAGS= -O3 -Wall -Wextra -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace
# FFLAGS= -O3 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace

LFLAGS= $(FFLAGS)

EXE=    pcsaft

SRC= _mod_basic_variables.f90  levenberg_marquardt.f90  minimizer.f90  toms778.f90  load_pure_and_kij_parameters.f90  input_output.f90  eos_constants.f90  ideal_gas_enthalpy.f90  eos_derivatives.f90  eos_properties.f90  system_definition.f90  density_iteration.f90  matrix_inversion2.f90  Newton_optimizer.f90  hbrd.f90  stability_analysis.f90  stability_analysis_G.f90  flash_algorithm.f90  bubble_dew_point.f90  critical_point.f90  phase_or_property_diagrams_pure.f90  phase_or_property_diagrams_mixture.f90  phase_equilibrium_iteration_mixture.f90  fit_pure_comp_parameters.f90  fit_kij_parameters.f90  main_prog.f90

#
## No need to edit below this line
#
#
#.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core *.mod *genmod.f90
