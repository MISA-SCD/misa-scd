README file for Parallel Spatially Resolved Stochastic Cluster Dynamics Software -- MISA-SCD

Chendandan
email: chendandan@xs.ustb.edu.cn

******************************************************************************************
Compiling instructions
******************************************************************************************
To compile, use the following commands:

module load [compiler]
cd src
make clean
make

To run under examples/, use the following commands:

cd ../tests/[Example file name]
mpirun -n [number of processes] ../../src/misascd configure.in

******************************************************************************************
Description of configure.in
******************************************************************************************
This is the main file containing the simulation parameters. It also contains the filenames 
of the other input files and several toggles for various simulation options. The options
that can be toggled in configure.in are:

1) ImplantType ('cascade' or 'FrenkelPair')
This option determines whether we are implanting cascades or Frenkel pairs. It toggles the
use of an adaptive meshing scheme for cascades as well as causes the program to read in
a list of cascades

2) ImplantScheme ('MonteCarlo' or 'explicit')
This option is only used when ImplantType is 'cascade'. It determines if cascade implantation
is carried out through the Monte Carlo algorithm or if cascades are introduced automatically
with a rate determined by the DPA rate and the volume element size.

The explicit cascade introduction scheme is designed for more efficient parallel operation
but has not been fully investigated.

3) MeshingType ('adaptive' or 'nonAdaptive')
This option toggles the use of an adaptive meshing scheme for damage implantation. 
Currently the code is set up so that adaptive meshing MUST be used with cascade damage. 
Adaptive meshing should not be used with Frenkel pair damage.

4) implantDist ('uniform' or 'nonUniform')
This option has only been applied to the 'nonAdaptive' meshing scheme but should work in
the adaptive meshing scheme as well. It allows a non-uniform spatial distribution of dose
rates and He implantation rates to be created based on an input file. This spatial distribution
is non-uniform in the z-direction only (it is one-dimensional). This option is designed for
use in simulating ion-irradiated thin foils.

5) GrainBoundaries ('yes' or 'no')
This option toggles the ability of grain boundaries that are not explicitly modeled to remove
migrating defects from the simulation.

6) pointDefect ('yes' or 'no')
This option toggles the mobility of defect clusters.

All other parameters are explained in the comments of the configure.in file

!*****************************************************************************************
!Output files and formatting
!*****************************************************************************************
The simulation outputs using a logarithmic timescale.

A new output file is created for each repetition of the simulation. MISA-SCD has an option
to run multiple simulations in a row.

At each output point, information is entered into two output files:
totdat_x.out

(here, 'x' refers to the simulation number)
totdat_x.out contains defect numbers in the entire simulation volume as well as the statistics for defects.

The current method for identifying defect types is as  follows:
NumCu NumV NumSIA NumSIA_immobile NumDefects

For example, the output:
3 2 0 0 5

Indicates that the defect type in this volume element is a Cu3V2 cluster, and there are 5
of these defects in the volume element.

!*****************************************************************************************
!General comments
!*****************************************************************************************
Each input file has a specific format that is currently hard-coded into the MISA-SCD code.
Questions about the formats of any of these files should be directed to CDD.

Currently the code has an annealing option, in which the simulation continues with zero
dose rate at a specified temperature for a specified amount of time.

The simulation automatically divides the mesh between the number of specified processors.
If a processor is left with zero mesh elements, the simulation will return an error.
