# $Header: /home/CVS//srscd/src/README.txt,v 1.2 2015/02/17 23:15:38 rdingre Exp $
README file for Parallel Spatially Resolved Stochastic Cluster Dynamics (SRSCD_par)

Aaron Dunn
email: aydunn@sandia.gov, adunn32@me.gatech.edu

******************************************************************************************
Compiling instructions
******************************************************************************************
First load openmpi under version 1.8.
Then, to compile, use the following commands:

cd src
make clean
make

To run, such as Bulk_Implantation_Cascade, use the following commands:

cd ../tests/Bulk_Implantation_Cascade
mpirun -n 8 ../../src/srscd parameters.txt

******************************************************************************************
Description of parameters.txt
******************************************************************************************

This is the main file containing the simulation parameters. It also contains the filenames 
of the other input files and several toggles for various simulation options. The options
that can be toggled in parameters.txt are:

1) MeshType ('uniform' or 'non-uniform')

This indicates to the code whether a uniform cubic mesh is being implemented or if a
non-uniform mesh is being implemented. The case of a non-uniform mesh has not been fully
implemented in the SRSCD code and is not currently operational.

To generate a mesh, run the (separate) MeshGenerator.f90 (using a serial fortran compiler
like f90 or gfortran). To determine the size of the mesh, edit MeshGenInput.txt. Within
MeshGenInput.txt, you can toggle periodic boundary conditions or free surfaces in the z-
direction.

2) DebugRestart ('yes' or 'no')

This option toggles a code restart from a debugging file. This option was created because
of code crashes at very high dose, which took several days to reach. Therefore, given a
debug restart file ('reset.in' is the default filename), you can resume a simulation from
a given point.

The reset.in file has very specific formatting and must include all of the defects in the
simulation, which volume elements they are located in, and other simulation data such
as the dose, the number of implantation events, the elapsed time, etc. Contact Aaron for
more detailed assistance with using this function.

3) NumMaterials (integer number)

This option allows simulation of multiple materials. Each volume element has a material 
ID. There should be one material input file for each material type.

More than 1 material type has not yet been implemented in SRSCD.

4) ImplantType ('cascade' or 'FrenkelPair')

This option determines whether we are implanting cascades or Frenkel pairs. It toggles the
use of an adaptive meshing scheme for cascades as well as causes the program to read in
a list of cascades

5) ImplantScheme ('MonteCarlo' or 'explicit')

This option is only used when ImplantType is 'cascade'. It determines if cascade implantation
is carried out through the Monte Carlo algorithm or if cascades are introduced automatically
with a rate determined by the DPA rate and the volume element size.

The explicit cascade introduction scheme is designed for more efficient parallel operation
but has not been fully investigated.

6) MeshingType ('adaptive' or 'nonAdaptive')

This option toggles the use of an adaptive meshing scheme for damage implantation. 
Currently the code is set up so that adaptive meshing MUST be used with cascade damage. 
Adaptive meshing should not be used with Frenkel pair damage.

7) implantDist ('uniform' or 'nonUniform')

This option has only been applied to the 'nonAdaptive' meshing scheme but should work in
the adaptive meshing scheme as well. It allows a non-uniform spatial distribution of dose
rates and He implantation rates to be created based on an input file. This spatial distribution
is non-uniform in the z-direction only (it is one-dimensional). This option is designed for
use in simulating ion-irradiated thin foils.

8) GrainBoundaries ('yes' or 'no')

This option toggles the ability of grain boundaries that are not explicitly modeled to remove
migrating defects from the simulation.

All other parameters are explained in the comments of the parameters.txt file

!*****************************************************************************************
!Output files and formatting
!*****************************************************************************************

The simulation outputs using a logarithmic timescale: output is generated at:

time=totalTime/200*2^n, with n between 0 and 7, as well as when time=totalTime, for a total
of 8 output points per simulation.

A new output file is created for each repetition of the simulation. SRSCD has an option
to run multiple simulations in a row.

At each output point, information is entered into three output files:
postpr_x.out
rawdat_x.out
totdat_x.out

(here, 'x' refers to the simulation number)

postpr_x.out contains user-defined post-processing information on the number and type 
of defects or reactions in the system (hard-coded into SRSCD)

rawdat_x.out contains a list of all volume elements in the system as well as all defects
present in each volume element. The current method for identifying defect types is as 
follows:

NumHe NumV NumSIA NumSIA_immobile NumDefects

For example, the output:

3 2 0 0 5

Indicates that the defect type in this volume element is a He3V2 cluster, and there are 5 
of these defects in the volume element.

totdat_x.out also contains defect numbers but defects are not separated by volume element
and are instead combined over the entire simulation volume.

An additional output file has been created for non-uniform implantation distributions. The 
output of this analysis is given in:
DepthProfile_x.txt

This gives the He, Vacancy, SIA_mobile, and SIA_sessile concentrations in the material as
a function of z-coordinate. It is only created at the end of the simulation.

!*****************************************************************************************
!General comments
!*****************************************************************************************

The various input files are below:

AnnealedCascades.txt !Contains list of defects in cascades
ExampleMaterial.txt !Contains all allowed reactions and migration/binding energies
ImplantProfile.txt !Contains non-uniform DPA and He/DPA ratio information
MeshInput_New.txt !Contains mesh
reset.in !Reset file for restarting simulations from non-zero DPA

Each input file has a specific format that is currently hard-coded into the SRSCD code.
Questions about the formats of any of these files should be directed to Aaron.

Currently the code has an annealing option, in which the simulation continues with zero
dose rate at a specified temperature for a specified amount of time. The correct initialization
of this annealing period and output of information in a convenient user-defined way is not
optimized because annealing has rarely been used. However, this section of the code does
work.

The simulation automatically divides the mesh between the number of specified processors.
If a processor is left with zero mesh elements, the simulation will return an error.
