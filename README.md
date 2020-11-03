# MISA-SCD v1.0

README file for Parallel Stochastic Cluster Dynamics Simulation Software -- MISA-SCD v1.0
Developers: [Chendandan](mailto:chendandan_ustb@xs.ustb.edu.cn)

## Compile and Run
Dependency:
1. mpi
2. fortran90

To compile:
```bash
$ module load [compiler]
$ cd $MISA_SCDv1.0_PATH/src
$ make
```

To run:
```bash
$ cd $MISA_SCDv1.0_PATH/examples/[Example file name]
$ mpirun -n [number of processes] ../../src/misascd configure.in
```

## Description of configure.in
This is the main file containing the simulation parameters. It also contains the filenames 
of the other input files and several toggles for various simulation options. It contains four
parts of parameters:
1. **Toggles**, which controls the simulation process.
    * implantType ('cascade' or 'FrenkelPair')
    This option determines whether we are implanting cascades or Frenkel pairs. It toggles the
    use of an adaptive meshing scheme for cascades as well as causes the program to read in
    a list of cascades

    * implantScheme ('MonteCarlo' or 'explicit')
    This option is only used when implantType is 'cascade'. It determines if cascade implantation
    is carried out through the Monte Carlo algorithm or if cascades are introduced automatically
    with a rate determined by the DPA rate and the volume element size.

    * meshingType ('adaptive' or 'nonAdaptive')
    This option toggles the use of an adaptive meshing scheme for cascade implantation.  But adaptive
    meshing should not be used for Frenkel pair implantation.

    * grainBoundaries ('yes' or 'no')
    This option toggles the ability of grain boundaries that are not explicitly modeled to remove
    migrating defects from the simulation.

    * pointDefect ('yes' or 'no')
    This option toggles the mobility of defect clusters.
2. **Simulation parameters**, which defaults to the parameters in the case of irradiation.
If it's a thermal aging simulation, one needs to insert 'agingTime' with a non-zero value
and set 'dpaRate' and 'totalDPA' with zero-value.

3. **Anneal parameters**, which used for annealing simulation.
If it's an anneal simulation, one should set the values and set 'dpaRate' and 'totalDPA' with zero-value.

4. **Output parameters**, which used to control post-processing.
The simulation outputs using a logarithmic timescale. A new output file is created for each repetition of
the simulation. MISA-SCD has an option to run multiple simulations in a row. At each output point, information can be entered into four output files:  
totdat_x.out, defect_x.out, stadat_x.out, xyzdat_x.out

* totdat_x.out: contains defect populations (defect type and their number) and their statistics (number density, average radius and so on) in the whole system.  
* defect_x.out: contains only defect populations (defect type and their number) in the whole system.  
* stadat_x.out: contains only statistics of defects (number density, average radius and so on) in the whole system.  
* xyzdat_x.out: contains defect populations (defect type and their number) and their statistics (number density, average radius and so on) in each mesh.

All other parameters are explained in the comments of the configure.in file

NOTE: the current method for identifying defect types is as  follows:  
NumCu NumV NumSIA NumSIA_immobile NumDefects

For example, the output:  
3 2 0 0 5  
Indicates that the defect type in this volume element is a Cu3V2 cluster, and there are 5 of these defects in the volume element.

## Others
Each input file has a specific format that is currently MISA-SCD v1.0 code. If the format is incorrect, the code will return errors.

The simulation automatically divides the mesh between the number of specified processors. If a processor is left with zero mesh elements, the code will return an error.

The development of the MISA-SCD v1.0 refers to the [SRSCD](https://github.com/aaronydunn/srscd/tree/master/Documents/GradSchool2016/srscd) (Spatially Resolved Stochastic Cluster Dynamics) code.



