!****************************************************************************
This is a readme file for He desorption simulations in Fe

Description:
This is a typical simulation type that can be compared to experiment, as He
desorption experiments on irradiated thin foils have been carried out for several
decades.

In a He desorption experiment, a thin foil is irradiated with He ions (creating
damage in the form of Frenkel pairs and He ions) at room temperature or below.

The foil is then annealed and the fraction of helium retained in the foil as 
a function of time is recorded. Annealing is typically only carried out at a 
single (constant) temperature.

In this simulation, helium implantation is carried out uniformly throughout a 
thin foil at room temperature and then annealing at higher temperature is carried
out. The post-processing for the fraction of helium remaining is not explicitly
carried out, but this is simple to do using the outputs that are provided by the
simulation.

NOTE: in this simulation, annealing is carried out in 20 steps of equal increment
(without changing the temperature) in order to output the He composition as a function
of time during annealing. It is better to output on a log scale in time, but this
functionality has not been included in the code yet. This is a trivial change
to the annealing loop in SRSCD.
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
Number of annealing steps = 20 here (only anneal at a single temperature)
Annealing temperature is the same at every step, but we want regular output.
Total anneal time set to 1000 seconds (therefore output every 50 seconds)
