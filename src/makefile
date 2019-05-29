#
#   $Header: 
#
#SHELL = /bin/sh

# -------------------End-macro-Defs---------------------------
#define all fortan source code

FF = mpif90 
FFLAGS =  
CC = 
CCFLAGS = 
LINK = mpif90
LINKFLAGS = 
USRLIB = 
SYSLIB = 
TARGET = srscd

# 
$(TARGET): mod_derivedtypes.o mod_srscd_constants.o  dprand.o Debug_subroutines.o Defect_attributes.o MeshReader.o ReactionRates.o Initialization_subroutines.o Read_inputs.o Cascade_implantation.o SRSCD_par.o kMC_subroutines.o Deallocate_Lists.o FineMesh_subroutines.o CoarseMesh_subroutines.o Misc_functions.o StrainSubroutines.o Postprocessing_srscd.o PointDefect_concentration.o
		$(LINK) -o $(TARGET) *.o

#
mod_derivedtypes.o: mod_derivedtypes.f90
		$(FF) -c mod_derivedtypes.f90
mod_srscd_constants.o: mod_srscd_constants.f90
	$(FF) -c mod_srscd_constants.f90
#
dprand.o: dprand.f90
	$(FF) -c dprand.f90
Debug_subroutines.o: Debug_subroutines.f90
	$(FF) -c Debug_subroutines.f90
Defect_attributes.o: Defect_attributes.f90
	$(FF) -c Defect_attributes.f90
MeshReader.o: MeshReader.f90
	$(FF) -c MeshReader.f90
ReactionRates.o: ReactionRates.f90
	$(FF) -c ReactionRates.f90
Initialization_subroutines.o: Initialization_subroutines.f90
	$(FF) -c Initialization_subroutines.f90
Read_inputs.o: Read_inputs.f90
	$(FF) -c Read_inputs.f90
Cascade_implantation.o: Cascade_implantation.f90
	$(FF) -c Cascade_implantation.f90
SRSCD_par.o: SRSCD_par.f90
	$(FF) -c SRSCD_par.f90
kMC_subroutines.o: kMC_subroutines.f90
	$(FF) -c kMC_subroutines.f90
Deallocate_Lists.o: Deallocate_Lists.f90
	$(FF) -c Deallocate_Lists.f90
FineMesh_subroutines.o: FineMesh_subroutines.f90
	$(FF) -c FineMesh_subroutines.f90
CoarseMesh_subroutines.o: CoarseMesh_subroutines.f90
	$(FF) -c CoarseMesh_subroutines.f90
Misc_functions.o: Misc_functions.f90
	$(FF) -c Misc_functions.f90
StrainSubroutines.o: StrainSubroutines.f90
	$(FF) -c StrainSubroutines.f90
Postprocessing_srscd.o: Postprocessing_srscd.f90
	$(FF) -c Postprocessing_srscd.f90
PointDefect_concentration.o: PointDefect_concentration.f90
	$(FF) -c PointDefect_concentration.f90

# Cleans

clean:
	rm -f *.o *.mod  $(TARGET)

