!****************************************************************************
This is a readme file for the nonuniform He implantation test

Description:
This simulation is used to match defect accumulation results with experimental
helium ion implantation in thin films, both at low tempeatures and at elevated
temperatures, both with and without annealing.

The main development of this test is the use of a non uniform implantation profile.
This profile contains the helium and dpa rates as a function of depth, as taken
from SRIM calculations. An example file is given here. The format of the file 
MUST be the same as the format given in this example. NOTE that the spacing of
the data points in the nonuniform implantation file does not have to match the
spacing of the mesh elements, as long as the total length of the mesh is less 
than or equal to the total length of the nonuniform implantation data.

Here, we use a 100 nm thick thin film (same thickness in SRIM data). The z-axis
gives depth, and the x- and y-axes are width. These are set to larger values
than the depth in order to have a larger total volume.

The DPA rate and HE/DPA ratio in the parameters.txt file should match the AVERAGE
values given in the nonuniform implantation file.
!****************************************************************************

Input files:

MeshGenInput.txt: provides the input for MeshGen, which creates a bulk mesh (80 nm elements, 2x2x2)

Mesh_ThinFilm.txt: the mesh created by MeshGen

Fe_Defects.txt: provides list of allowed reactions in Fe, defect types, and defect binding and migration energies

ImplantProfile.txt: gives the DPA rate and He implant rate (in at. fraction) as a function of depth.

parameters.txt: contains simulation parameters as well as the names of all of the other input files

The parameters.txt file contains the main simulation parameters:

Non-adaptive meshing
Temperature
Dose Rate (average)
Total Dose (average)

Anneal tempeature
Number of annealing steps = 1 here (only anneal at a single temperature)
total anneal time = 0 s for high temperature implant and 600s for low T implant and annealing
