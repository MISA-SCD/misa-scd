This is a readme file for the bulk cascade implantation simulation test

Input files:

MeshGenInput.txt: provides the input for MeshGen, which creates a bulk mesh (80 nm elements, 2x2x2)

Mesh_Bulk.txt: the mesh created by MeshGen

Fe_Cascades.txt: provides cascades for implantation (20 keV cascades)

Fe_Defects.txt: provides list of allowed reactions in Fe, defect types, and defect binding and migration energies

parameters.txt: contains simulation parameters as well as the names of all of the other input files

Description:

This is a bulk cascade implantation in 343 K Fe. The options not used in this simulation are: debug restart, nonuniform implantation profile

The parameters.txt file contains the main simulation parameters:

Adaptive meshing (5nm fine mesh elements)
Temperature
Dose Rate
Total Dose
Grain Size
Dislocation Density (for sinks)

NOTE: if the MonteCarlo option is chosen in parameters.txt, it is suggested
that the simulation be run using srscd2n (or bgsrscd2n)

Other options that can be explored using this input file are:

MonteCarlo vs Explicit cascade implantation - take care of cascade implantation using the regular Monte Carlo scheme or using an explicit scheme
This provides better weak scaling in parallel simulations

Cascade size: currently set at 1000 nm for 20 keV cascades. This gives the probability that defects already in the system will interact with cascade defects.
