# -------------------Begin-------------------------
#$Header:
#
#SHELL = /bin/sh
# -------------------End---------------------------

#Compiler and Compile options
FC = mpif90
FF = -g

#Source code files, Object files and Target(Executable) file defeine
SRCS = $(wildcard src_*.f90)
OBJS = $(SRCS:.f90=.o)
TARGET = misascd

MODSRCS =   mod_structures.f90\
            mod_constants.f90\
            mod_globalVariables.f90\
            mod_randdp.f90\
            mod_reactionrates.f90
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
	rm $(TARGET) $(OBJS) $(MODOBJS) $(MODS)

