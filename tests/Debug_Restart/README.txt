!****************************************************************************
This is a readme file for the debug restart test test

Description:
In this test, the SRSCD simulation will restart from a set point given in the
file debugRestart.in. DebugRestart.in MUST be in the same format as the file
in this test folder. This file contains the elapsed time, number of
implantation events, number of He implant events, and number of processors to start with in the
simulation. Then it contains a list of all defects present in each mesh
element, with defect type and number given. Note that the simulation does not
include fine meshes in the debug restart, so there are no cascades present at
the beginning of the restart algorithm.

Note that the number of processors you choose must match the number of
processors in the debug restart file. Therefore in the example provided here,
SRSCD must be run with 4 processors.

This functionality is used primarily for debugging problems that occur at
large computation times. It could also be used in the case that you want to
run one simulation up to a certain point, then several simulations starting at
that point (for example, annealing of a given population).
!******************************************************************************

Input files:

MeshGenInput.txt: provides the input for MeshGen, which creates a bulk mesh (80 nm elements, 2x2x2)

Mesh_Bulk.txt: the mesh created by MeshGen

Fe_Cascades.txt: provides cascades for implantation (20 keV cascades)

Fe_Defects.txt: provides list of allowed reactions in Fe, defect types, and defect binding and migration energies

parameters.txt: contains simulation parameters as well as the names of all of the other input files

debugRestart.in: contains the simulation information (elapsed time, number of
cascades) as well as the information on defects present in the simulation at
the beginning of the simulation.

The parameters.txt file contains the main simulation parameters:

Adaptive meshing (5nm fine mesh elements)
Temperature
Dose Rate
Total Dose
Grain Size
Dislocation Density (for sinks)

