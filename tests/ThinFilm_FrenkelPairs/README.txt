!****************************************************************************
This is a readme file for the uniform Frenkel pair implantation in Fe thin films test simulation

Description:
This is an example simulation of Frenkel pair implantation in a Fe thin film.
Implantation is assumed to be uniform throughout the material and no helium
is present in this simulation.

The dose rate and implantation temperature are specified in parameters.txt.

The layer thickness is currently set at 100 nm. This can be changed by adjusting
MeshGenInput.txt and running MeshGen. The free surfaces of the thin film are 
located at z=0 and z=100 nm (free surfaces parallel to x-y plane)

This test simulation is intended to provide the user with a simple simulation
with which to begin understanding the various inputs and outputs that SRSCD provides.
!****************************************************************************

Input files:

MeshGenInput.txt: provides the input for MeshGen, which creates a bulk mesh (80 nm elements, 2x2x2)

Mesh_ThinFilm.txt: the mesh created by MeshGen

Fe_Defects.txt: provides list of allowed reactions in Fe, defect types, and defect binding and migration energies

parameters.txt: contains simulation parameters as well as the names of all of the other input files

The parameters.txt file contains the main simulation parameters:

Non-adaptive meshing
Temperature
Dose Rate (average)
Total Dose (average)

Anneal tempeature
Number of annealing steps = 1 here (only anneal at a single temperature)
total anneal time = 0 s for high temperature implant and 600s for low T implant and annealing
