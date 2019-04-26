!****************************************************************************
This is a readme file for the isochronal annealing test

Description:
In this test, the SRSCD simulation will implant Frenkel pairs in a bulk Fe
sample at a low (4.5 K) temperature. Then, the simulation will anneal this sample
in steps. 

The parameters chosen for this simulation (matching those of Fu et al. and 
Dalla Torre et al. isochronal annealing simulations) are:

begin annealing at 127 K
Take 44 annealing steps of 300 s each
Multiply temperature by 1.03 at each annealing step

The total dose and intitial annealing temperature are chosen to represent the
remaining defect content after stage 1D annealing (correlated recombination).

This simulation type can be used to identify which defect mechanisms are active at which temperatures.
!****************************************************************************

Input files:

MeshGenInput.txt: provides the input for MeshGen, which creates a bulk mesh (80 nm elements, 2x2x2)

Mesh_IsoAnneal.txt: the mesh created by MeshGen

Fe_Defects.txt: provides list of allowed reactions in Fe, defect types, and defect binding and migration energies

parameters.txt: contains simulation parameters as well as the names of all of the other input files

The parameters.txt file contains the main simulation parameters:

Non-adaptive meshing
Temperature
Dose Rate
Total Dose

Anneal tempeature
Number of annealing steps
total anneal time
amount to multiply anneal temperature by at each annealing step (chosen to be 1.03 here)
