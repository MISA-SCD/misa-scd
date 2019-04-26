!****************************************************************************
This is a readme file for the uniform Frenkel pair implantation in nanoporous test simulation

Description:
This is an example simulation of Frenkel pair implantation in a Fe nanoporous
bulk.
Implantation is assumed to be uniform throughout the material and no helium
is present in this simulation.

The dose rate and implantation temperature are specified in parameters.txt.

This test simulation is intended to provide the user with a simple simulation
with which to begin understanding the various inputs and outputs that SRSCD provides.
!****************************************************************************

Input files:

MeshGenInput.txt: provides the input for MeshGen

Mesh_Nanoporous.txt: the mesh created by MeshGen

Defects_Fe.txt: provides list of allowed reactions in Fe nanoporous bulk, defect types, and defect binding and migration energies

Defects_pore.txt: provides list of allowed reactions in pores (i.e. none)

parameters.txt: contains simulation parameters as well as the names of all of the other input files

The parameters.txt file contains the main simulation parameters:

Non-adaptive meshing
Temperature
Dose Rate (average)
Total Dose (average)

Anneal tempeature
Number of annealing steps = 1 here (only anneal at a single temperature)
total anneal time = 0 s for high temperature implant and 600s for low T implant and annealing
