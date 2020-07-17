# -------------------Begin-------------------------
#$Header:
#
#SHELL = /bin/sh
# -------------------End---------------------------

#Compiler and Compile options
FC = mpif90
FF = -g

#Source code files, Object files and Target(Executable) file defeine
SRCS =  dprand.f90\
        Debug_subroutines.f90\
        Defect_attributes.f90\
        MeshReader.f90\
        ReactionRates.f90\
        Initialization_subroutines.f90\
        Read_inputs.f90\
        Cascade_implantation.f90\
        MISASCD_par.f90\
        kMC_subroutines.f90\
        Deallocate_Lists.f90\
        FineMesh_subroutines.f90\
        CoarseMesh_subroutines.f90\
        Misc_functions.f90\
        Postprocessing_srscd.f90\
        PointDefect_concentration.f90
OBJS = $(SRCS:.f90=.o)
TARGET = misascd

MODSRCS =   mod_structures.f90\
            mod_constants.f90
MODOBJS = $(MODSRCS:.f90=.o)
MODS = $(MODSRCS:.f90=.mod)

#Generate Target file using Object files
$(TARGET): $(OBJS) $(MODOBJS)
	$(FC) -o $@ $^

#Generate Object files using Source files
%.o: %.f90 $(MODS)
	$(FC) -c $<

$(MODS) $(MODOBJS): $(MODSRCS)
	$(FC) -c $^

# Clean Object files
clean:
	rm $(TARGET) $(OBJS) $(MODOBJS) $(MODS) *.mod

