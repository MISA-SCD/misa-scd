!****************************************************************************************
!> Module mod_constants (list of globally shared variables and pointers)
!!
!! This module contains the list of all globally shared variables. This includes:
!! 1) Processor and mesh information (code backbone)
!! 2) Reaction lists, defect lists, and cascade lists, both in coarse and fine meshes
!! 3) Global reaction rates
!! 4) Material input information, in the form of derived types (diffusion and binding energies, etc)
!! 5) Universal constants
!! 6) Simulation parameters read in from parameters.txt
!! 7) Other miscellaneous variables used for MPI, debugging, or postprocessing
!****************************************************************************************
module mod_constants

    use mod_structures
    implicit none

    !used for MPI commands
    integer comm                            !<New communication domain: created by MPI_CART_CREAT
    integer dims(3)                         !<Number of processors in x, y, z.
    logical periods(3)                      !<Boundary conditions in x, y, z.  The default is periodic
    integer ierr							!<used for initializing and finalizing MPI
    integer, parameter :: MASTER=0			!<Define the master node as ID=0
    integer, parameter :: maxBufferSize=50	!<Used to define the max size of a send/recieve buffer

    !Processor and mesh information, created in MeshReader.f90
    type(processorData) myProc							!<Contains processor information (id, neighbors)
    type(mesh),allocatable :: myMesh(:)					!<Contains (local) mesh information
    type(boundaryMesh),allocatable ::  myBoundary(:,:)	!<Boundary elements (direction, element #)
    integer numCells						            !<Number of cells in local mesh
    integer numTotal                                    !<Total meshes in the sysytem
    double precision meshLength                         !<Length of a coarse mesh
    integer numx,numy,numz                              !<Number of global meshes in x-direction, y-direction, z-direction
    integer numxLocal,numyLocal,numzLocal               !<Number of local meshes in x-direction, y-direction, z-direction

    !reaction and defect lists
    type(reaction),pointer :: reactionList(:)			!<List of reactions in local (coarse) mesh
    type(defect),pointer :: defectList(:)			    !<List of defects in local (coarse) mesh
    double precision,allocatable :: totalRateVol(:)	    !<Total reaction rate in each volume element
    double precision totalRate						    !<Total reaction rate in this processor
    double precision maxRate						    !<Max reaction rate in all processors

    double precision totalTime                          !<Total time
    double precision elapsedTime                        !<Elapsed time
    integer step                                        !<Current number of time steps

    !List of cascades to be chosen from when cascade implantation occurs
    type(cascadeEvent),pointer :: cascadeList			!<List of cascades (read from file) that can be implanted
    !Fine meshes of cascades implanted in coarse mesh
    type(cascade),pointer :: ActiveCascades				!<List of fine meshes that are active due to recent cascade implantation. Contains defect lists and reaction lists)
    integer numCascades							        !<number of cascades in the cascade input file (used to choose cascades to input in simulation)
    integer numCellsCascade							    !<number of volume elements within a cascade (fine) mesh
    integer numxCascade,numyCascade,numzCascade	        !<number of elements in cascade x-direction, y-direction, z-direction
    !double precision, allocatable :: globalMeshCoord(:,:)
    integer,allocatable :: cascadeConnectivity(:,:) 	!<connectivity matrix for cascade meshes (same for all fine meshes)
    double precision fineLength						    !<length of a cascade volume element (nm)
    double precision cascadeElementVol				    !<volume of a cascade element (nm^3)

    !Material input information, created in MaterialInput.f90
    integer numSpecies							        !<Number of chemical species (typically set to 4: Cu, V, SIA_m, SIA_im)

    type(formationSingle),allocatable :: FormSingle(:)                  !<Parameters for formation of single defects--(numSingleForm)
    type(diffusionSingle),allocatable :: DiffSingle(:)			        !<Parameters for diffusion of single defects--(numSingleDiff)
    type(diffusionFunction), allocatable :: DiffFunc(:)			        !<Parameters for functional forms of diffusion rates for defects--(numFuncDiff)
    type(bindingSingle), allocatable :: BindSingle(:)				    !<Parameters for binding of single defects--(numSingleBind)
    type(bindingFunction), allocatable :: BindFunc(:)				    !<Parameters for functional forms of binding energies for defects--(numFuncBind)
    type(reactionParameters), allocatable :: DissocReactions(:)         !<List of allowed dissociation reactions (and ref. to functional form of reaction rate)--(numDissocReac)
    type(reactionParameters), allocatable :: DiffReactions(:)		    !<List of allowed diffusion reactions (and ref. to functional form of reaction rate)--(numDiffReac)
    type(reactionParameters), allocatable :: SinkReactions(:)		    !<List of allowed sink reactions (and ref. to functional form of reaction rate)--(numSinkReac)
    type(reactionParameters), allocatable :: ImpurityReactions(:)	    !<List of allowed impurity reactions (and ref. to functional form of reaction rate)--(numImpurityReac)
    type(reactionParameters), allocatable :: ClusterReactions(:)	    !<List of allowed clustering reactions (and ref. to functional form of reaction rate)--(numClusterReac)
    type(reactionParameters), allocatable :: ImplantReactions(:)	    !<List of allowed implantation reactions (and ref. to functional form of reaction rate)--(numImplantReac)

    integer :: numSingleForm                                !<Number of single defect formation energy in input file
    integer :: numSingleDiff	                            !<Number of single defect diffusion rates in input file
    integer :: numFuncDiff		                            !<Number of functional forms for diffusion rates in input files
    integer :: numSingleBind                                !<Number of single defect binding energies in input file
    integer :: numFuncBind		                            !<Number of functional forms for binding energies in input files

    integer :: numDissocReac	                            !<Number of dissociation reactions in input file
    integer :: numDiffReac		                            !<Number of diffusion reactions in input file
    integer :: numSinkReac		                            !<Number of sink reactions in input file
    integer :: numImpurityReac	                            !<Number of impurity reactions in input file
    integer :: numClusterReac	                            !<Number of clustering reactions in input file
    integer :: numImplantReac	                            !<Number of implantation reactions in input file (cascade, Frenkel pair currently implemented)

    !constants
    double precision, parameter :: kboltzmann=8.625d-5	    !<Boltzmann's constant (eV/K)
    double precision, parameter :: pi=3.141592653589793d0	!<Pi
    double precision, parameter :: Zint = 1.2d0				!<Constant representing preference for clustering of interstitials by interstitial clusters (increases clustering cross-section)
    double precision, parameter :: Zv = 1.0d0
    !double precision, parameter :: reactionRadius=0.65d0	!<Material parameter used for reaction distances (impacts reaction rates) (nm)
    !double precision, parameter :: lattice = 0.316d0       !<lattice constant (nm) (Fe: 0.2876d0; W: 0.316)
    double precision, parameter :: atomSize_Cu = 8.79d-3    !<Cu (nm^3)

    !Cu solubility CeqCu(T) = exp(DelatS/kB)*exp(-Omega/(kB*T))  Reference: (F. Christien and A. Barbu, 2004)
    double precision ceqV                   !Thermal equilibrium concentration of vacancy
    double precision ceqI                   !Thermal equilibrium concentration of SIA
    double precision concV                  !Vacancy concentration
    double precision concI                  !SIA concentration
    integer numCuCell                       !Initial number of Cu atoms in one mesh
    integer initialNumV                     !Initial number of vacancies in the whole system
    integer initialNumI                     !Initial number of self-interstitial atoms in the whole system
    integer, allocatable :: listVI(:,:)     !List the globalID of the mesh where initial vacancies and self-interstitial atoms are located.

    !simulation parameters, to be computed during simulation
    double precision temperature			!<Temperature (K)
    double precision DPA					!<DPA tracker (not a parameter)
    double precision numDisplacedAtoms		!<number of atoms displaced per cascade, read from cascade file
    double precision totalVolume			!<Volume of single processor's mesh (nm^3)
    double precision systemVol				!<Volume of all processors (global), not including grain boundary elements (nm^3)

    !simulation parameters, to be read during readParameters() in main program
    double precision tempStore				!<Temperature read in (K) - used when temp. changes several times during a simulation
    double precision CuContent              !<The initial content of Cu in iron
    integer numVac                          !<The number of vacancies put in
    integer numInt                          !<The number of interstitials put in
    double precision dpaRate				!<DPA rate in dpa/s
    double precision totalDPA				!<total DPA in simulation
    double precision firr                   !firr = Vconcent / initialCeqv. Radiation enhanced factor
    double precision atomSize				!<atomic volume (nm^3)
    double precision lattice                !<recombination radius (nm)
    double precision burgers				!<magnitude of burgers vector, equal to lattice constant
    double precision reactionRadius         !<recombination radius (nm)
    double precision agingTime              !<Thermal aging time (s)
    double precision meanFreePath			!<mean free path before a defect is absorbed by a grain boundary (AKA avg. grain size)
    double precision dislocationDensity		!<density of dislocations (sinks for point defects)
    double precision impurityDensity		!<denstiy of impurity atoms (traps for SIA loops)
    double precision cascadeVolume			!<Volume of cascade (used for cascade mixing probability)
    double precision cascadeReactionLimit	!<Total reaction rate in a cascade cell to consider it annealed and release cascade back to coarse mesh (s^-1)
    integer max3DInt			            !<largest SIA size that can diffuse in 3D as spherical cluster
    integer numGrains			            !<Number of grains inside polycrystal (default 1)
    integer numSims				            !<Number of times to repeat simulation

    !used for test2
    character(len=20) test3
    integer oneCascadeGCell

    !annealing information
    double precision annealTime				!<Amount of time for anneal (s)
    double precision annealTemp				!<Temperature of anneal stage (K)
    double precision annealTempInc			!<Temperature increment at each annealing step (additive or multipliciative)
    integer	annealSteps			            !<Number of annealing steps
    character(len=20) annealType		    !<('mult' or 'add') toggles additive or multiplicative anneal steps
    logical annealIdentify			        !<(.TRUE. if in annealing phase, .FALSE. otherwise) used to determine how to reset reaction rates (should we include implantation or not)
    integer annealIter                      !<Current number of time steps

    !Toggles
    character(len=20) implantScheme			!<(MonteCarlo or explicit), used to determine if cascades are implanted through Monte Carlo algorithm or explicitly
    character(len=20) implantDist			!<(Uniform or NonUniform), used to determine if defects are implanted uniformly or if DPA rate / He implant rate are given for each volume element
    character(len=20) grainBoundaryToggle	!<Used to determine whether or not we are using grain boundaries to remove defects from simulation
    character(len=20) pointDefectToggle		!<Toggles whether or not we allow HeSIA clusters to form ('yes' or 'no')
    character(len=20) meshType              !<Boundary conditions

    character(len=20) implantType			!<(Frenkel pairs or cascades), used to determine the type of damage in the simulation
    character(len=20) meshingType			!<(adaptive or nonAdaptive), used to determine whether we are simulating cascade implantation with adaptive meshing
    character(len=20) polycrystal			!<(yes or no), used to identify whether or not we have multiple grains in our crystal

    !Output  parameters
    character(len=20) totdatToggle			!<(yes or no), used to toggle whether we output the totdat.out data file
    character(len=20) rawdatToggle			!<(yes or no), used to toggle whether we output the rawdat.out data file
    integer minSCluster                     !<Only n>minSCluster SnVm and Sn clusters are counted
    integer minVoid                         !<Only n>minVoid Vn clusters are counted
    integer minLoop                         !<Only n>minLoop SIAn clusters are counted
    integer minSV                           !<Only (n+m)>minCuV SnVm clusters are counted

    !<input files
    integer, parameter :: PARAFILE = 10                     !<Used to read parameter.txtx file
    integer, parameter :: DEFFILE = 11                      !<Used to read Defects.txtx file
    integer, parameter :: MESHFILE = 12                     !<Used to read Mesh_*.txt file
    integer, parameter :: CASFILE = 13                      !<Used to read cascades.txt File
    !<output file
    integer, parameter :: RAWDAT = 82
    integer, parameter :: TOTDAT = 83

    double precision omega					!<Geometric constant for 3D spherical clustering (see Dunn et al. JNM 2013)
    double precision omega2D				!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision omega1D				!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision omegastar				!<Geometric constant for 3D spherical clustering (see Dunn et al. JNM 2013)
    double precision omegastar1D			!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision omegacircle1D			!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision recombinationCoeff		!<Geometric constant for Frenkel pair recombination (see Dunn et al. JNM 2013)

    !counters
    double precision rateTau(2)             !<Used for collective communication
    integer numImpAnn(2)                    !<Postprocessing: numImpAnn(1) is the num of Frenkel pairs / cascades (local), numImpAnn(2) is the number of annihilation reactions carried out (local)
    integer totalImpAnn(2)                  !<Postprocessing: numImpAnn(1) is the number of implant events across all processors, numImpAnn(2) is the number of annihilation reactions across all processors

    !counters for sink efficiency
    integer numTrapV			            !<Postprocessing: number of vacancies trapped on grain boundary
    integer numTrapSIA			            !<Postprocessing: number of SIAs trapped on grain boundary
    integer numEmitV			            !<Postprocessing: number of vacancies emitted from grain boundary
    integer numEmitSIA			            !<Postprocessing: number of SIAs emitted from grain boundary

contains

end module
