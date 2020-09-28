! Created by CDD on 2020/9/28.
!****************************************************************************************
!> Module constants (list of all constants).
!****************************************************************************************
module mod_constants
    implicit none

    integer, parameter :: MASTER=0			 !<Define the master node as ID=0
    integer, parameter :: SPECIES = 4        !<Number of chemical species (typically set to 4: Cu, V, SIA_m, SIA_im)

    !constants
    double precision, parameter :: kboltzmann=8.625d-5	    !<Boltzmann's constant (eV/K)
    double precision, parameter :: pi=3.141592653589793d0	!<Pi
    double precision, parameter :: Zint = 1.2d0				!<Constant representing preference for clustering of interstitials by interstitial clusters (increases clustering cross-section)
    double precision, parameter :: Zv = 1.0d0
    double precision, parameter :: atomSize_Cu = 8.79d-3    !<Cu (nm^3)

    !<input files
    integer, parameter :: PARAFILE = 10     !<Used to read parameter.txtx file
    integer, parameter :: ATTRFILE = 11     !<Used to read Defects.txtx file
    integer, parameter :: MESHFILE = 12     !<Used to read Mesh_*.txt file
    integer, parameter :: CASFILE = 13      !<Used to read cascades.txt File
    !<output file
    integer, parameter :: TOTFILE = 81      !<Used to write totdat.out file, contains defects (type and number) and  statistical data
    integer, parameter :: DEFFILE = 82      !<Used to write defect.out file, contains onle defects (type and number)
    integer, parameter :: STAFILE = 83      !<Used to write stadat.out file, contains only statistical data
    integer, parameter :: XYZFILE = 84      !<Used to write xyzdat.out file, contains defects (type and number) in each mesh

end module mod_constants