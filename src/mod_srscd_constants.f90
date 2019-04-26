! $Header: /home/CVS//srscd/src/mod_srscd_constants.f90,v 1.21 2015/12/14 21:34:49 aydunn Exp $
!****************************************************************************************
!> Module mod_SRSCD_constants (list of globally shared variables and pointers)
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

module mod_srscd_constants

use DerivedType
implicit none

!Processor and mesh information, created in MeshReader.f90
type(processorData) myProc								!<Contains processor information (id, neighbors)
type(mesh), allocatable :: myMesh(:)					!<Contains (local) mesh information
type(boundaryMesh), allocatable ::  myBoundary(:,:)		!<Boundary elements (direction, element #)
integer numCells										!<Number of cells in local mesh

!reaction and defect lists
type(reaction), pointer :: reactionList(:)				!<List of reactions in local (coarse) mesh
type(defect), pointer :: DefectList(:)					!<List of defects in local (coarse) mesh
double precision totalRate								!<Total reaction rate in this processor
double precision, allocatable :: totalRateVol(:)		!<Total reaction rate in each volume element
double precision maxRate								!<Max reaction rate in all processors

!List of cascades to be chosen from when cascade implantation occurs
type(cascadeEvent), pointer :: cascadeList				!<List of cascades (read from file) that can be implanted

!Fine meshes of cascades implanted in coarse mesh
type(cascade), pointer :: ActiveCascades				!<List of fine meshes that are active due to recent cascade implantation
														!!(Contains defect lists and reaction lists)
integer numCascades										!<number of cascades in the cascade input file (used to choose cascades to input in simulation)
integer numCellsCascade									!<number of volume elements within a cascade (fine) mesh
integer numxCascade										!<number of elements in cascade x-direction
integer numyCascade										!<number of elements in cascade y-direction
integer numzCascade										!<number of elements in cascade z-direction
integer, allocatable :: cascadeConnectivity(:,:) 		!<connectivity matrix for cascade meshes (same for all fine meshes)
double precision fineLength								!<length of a cascade volume element (nm)
double precision cascadeElementVol						!<volume of a cascade element (nm^3)

!Material input information, created in MaterialInput.f90
integer numMaterials						!<Number of material types (eg. copper, niobium or bulk, grain boundary)
integer numSpecies							!<Number of chemical species (typically set to 4: He, V, SIA_glissile, SIA_sessile)

type(diffusionFunction), allocatable :: DiffFunc(:,:)			!<Parameters for functional forms of diffusion rates for defects
type(diffusionSingle), allocatable :: DiffSingle(:,:)			!<Parameters for diffusion of single defects
type(bindingSingle), allocatable :: BindSingle(:,:)				!<Parameters for binding of single defects
type(bindingFunction), allocatable :: BindFunc(:,:)				!<Parameters for functional forms of binding energies for defects
type(reactionParameters), allocatable :: DissocReactions(:,:)	!<List of allowed dissociation reactions (and ref. to functional form of reaction rate)
type(reactionParameters), allocatable :: DiffReactions(:,:)		!<List of allowed diffusion reactions (and ref. to functional form of reaction rate)
type(reactionParameters), allocatable :: SinkReactions(:,:)		!<List of allowed sink reactions (and ref. to functional form of reaction rate)
type(reactionParameters), allocatable :: ImpurityReactions(:,:)	!<List of allowed impurity reactions (and ref. to functional form of reaction rate)
type(reactionParameters), allocatable :: ClusterReactions(:,:)	!<List of allowed clustering reactions (and ref. to functional form of reaction rate)
type(reactionParameters), allocatable :: ImplantReactions(:,:)	!<List of allowed implantation reactions (and ref. to functional form of reaction rate)

integer, allocatable :: numSingleDiff(:)	!<Number of single defect diffusion rates in input file
integer, allocatable :: numFuncDiff(:)		!<Number of functional forms for diffusion rates in input files
integer, allocatable :: numSingleBind(:)	!<Number of single defect binding energies in input file
integer, allocatable :: numFuncBind(:)		!<Number of functional forms for binding energies in input files

integer, allocatable :: numDissocReac(:)	!<Number of dissociation reactions in input file
integer, allocatable :: numDiffReac(:)		!<Number of diffusion reactions in input file
integer, allocatable :: numSinkReac(:)		!<Number of sink reactions in input file
integer, allocatable :: numImpurityReac(:)	!<Number of impurity reactions in input file
integer, allocatable :: numClusterReac(:)	!<Number of clustering reactions in input file
integer, allocatable :: numImplantReac(:)	!<Number of implantation reactions in input file (cascade, Frenkel pair, He currently implemented)

!constants
double precision, parameter :: kboltzmann=8.6173324d-5	!<Boltzmann's constant (eV/K)
double precision, parameter :: pi=3.141592653589793		!<Pi
double precision, parameter :: Zint = 1.2				!<Constant representing preference for clustering of interstitials by interstitial clusters (increases clustering cross-section)
double precision, parameter :: reactionRadius=.5065d0	!<Material parameter used for reaction distances (impacts reaction rates)

!simulation parameters, to be read during readParameters() in main program
double precision temperature			!<Temperature (K)
double precision tempStore				!<Temperature read in (K) - used when temp. changes several times during a simulation
double precision HeDPARatio				!<Helium to dpa ratio (atoms per atom)
double precision DPARate				!<DPA rate in dpa/s
double precision atomsize				!<atomic volume (nm^3)
double precision DPA					!<DPA tracker (not a parameter)
double precision defectDensity			!<total density of defects (?), not sure if this is needed
double precision dislocationDensity		!<density of dislocations (sinks for point defects)
double precision impurityDensity		!<denstiy of impurity atoms (traps for SIA loops)
double precision totalDPA				!<total DPA in simulation
double precision burgers				!<magnitude of burgers vector, equal to lattice constant
double precision numDisplacedAtoms		!<number of atoms displaced per cascade, read from cascade file
double precision meanFreePath			!<mean free path before a defect is absorbed by a grain boundary (AKA avg. grain size)
double precision cascadeReactionLimit	!<Total reaction rate in a cascade cell to consider it annealed and release cascade back to coarse mesh (s^-1)
double precision cascadeVolume			!<Volume of cascade (used for cascade mixing probability)
double precision totalVolume			!<Volume of single processor's mesh (nm^3)
double precision systemVol				!<Volume of all processors (global), not including grain boundary elements (nm^3)

!Sink efficiency parameters
double precision alpha_v				!<Grain boundary sink efficiency for vacancies
double precision alpha_i				!<Grain boundary sink efficiency for interstitials
double precision conc_v					!<Concentration of vacancies found by GB model (used to fit alpha_v)
double precision conc_i					!<Concentration of interstitials found by GB model (used to fit alpha_i)

!annealing information
double precision annealTime				!<Temperature of anneal stage (K)
double precision annealTemp				!<Amount of time for anneal (s)
integer			 annealSteps			!<Number of annealing steps
double precision annealTempInc			!<Temperature increment at each annealing step (additive or multipliciative)
character*20 	 annealType				!<('mult' or 'add') toggles additive or multiplicative anneal steps
logical 		 annealIdentify			!<(.TRUE. if in annealing phase, .FALSE. otherwise) used to determine how to reset reaction rates (should we include implantation or not)

integer numSims							!<Number of times to repeat simulation
integer max3DInt						!<largest SIA size that can diffuse in 3D as spherical cluster
integer SIAPinMin						!<Smallest size of SIA that can pin at HeV clusters
integer numGrains						!<Number of grains inside polycrystal (default 1)

character*20 implantType				!<(Frenkel pairs or cascades), used to determine the type of damage in the simulation
character*20 grainBoundaryToggle		!<Used to determine whether or not we are using grain boundaries to remove defects from simulation
character*20 HeSIAToggle				!<Toggles whether or not we allow HeSIA clusters to form ('yes' or 'no')
character*20 SIAPinToggle				!<Toggles whether or not we allow HeV clusters to pin SIA clusters (without annihilating V+SIA)
character*20 meshingType				!<(adaptive or nonAdaptive), used to determine whether we are simulating cascade implantation with adaptive meshing
character*20 implantScheme				!<(MonteCarlo or explicit), used to determine if cascades are implanted through Monte Carlo algorithm or explicitly
character*20 implantDist				!<(Uniform or NonUniform), used to determine if defects are implanted uniformly or if DPA rate / He implant rate are given for each volume element
character*20 polycrystal				!<(yes or no), used to identify whether or not we have multiple grains in our crystal
character*20 vtkToggle					!<(yes or no), used to toggle whether we want vtk output at each time increment (log scale)
character*20 outputDebug				!<(yes or no), used to toggle whether we want to output a debug restart file at each time increment
character*20 singleElemKMC				!<(yes or no), used to toggle whether we are making one kMC choice per volume element or one kMC choice for the whole processors
character*20 sinkEffSearch				!<(yes or no), used to toggle search for effective sink efficiency
character*20 strainField				!<(yes or no), used to toggle whether we are simulating diffusion in a strain field	
character*20 postprToggle				!<(yes or no), used to toggle whether we output the postpr.out data file
character*20 totdatToggle				!<(yes or no), used to toggle whether we output the totdat.out data file
character*20 rawdatToggle				!<(yes or no), used to toggle whether we output the rawdat.out data file
character*20 xyzToggle					!<(yes or no), used to toggle whether we output an .xyz data file (for visualization)	
character*20 profileToggle				!<(yes or no), used to toggle whether we output a DefectProfile.out data file		

!(hard-coded) constants used for clustering rates
double precision omega					!<Geometric constant for 3D spherical clustering (see Dunn et al. JNM 2013)
double precision omega2D				!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
double precision omega1D				!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
double precision omegastar				!<Geometric constant for 3D spherical clustering (see Dunn et al. JNM 2013)
double precision omegastar1D			!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
double precision omegacircle1D			!<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
double precision recombinationCoeff		!<Geometric constant for Frenkel pair recombination (see Dunn et al. JNM 2013)

!used for MPI commands
integer ierr							!<used for initializing and finalizing MPI
integer, parameter :: MASTER=0			!<Define the master node as ID=0
integer, parameter :: maxBufferSize=50	!<Used to define the max size of a send/recieve buffer

!counters
integer numImplantEvents			!<Postprocessing: number of Frenkel pairs / cascades (local)
integer numHeImplantEvents			!<Postprocessing: number of He implantation events (local)
integer totalImplantEvents			!<Postprocessing: number of implant events across all processors 
integer numHeImplantTotal			!<Postprocessing: number of He implant events across all processors
integer numAnnihilate				!<Postprocessing: number of annihilation reactions carried out

!counters for sink efficiency
integer numTrapV					!<Postprocessing: number of vacancies trapped on grain boundary
integer numTrapSIA					!<Postprocessing: number of SIAs trapped on grain boundary
integer numEmitV					!<Postprocessing: number of vacancies emitted from grain boundary
integer numEmitSIA					!<Postprocessing: number of SIAs emitted from grain boundary

!DEBUG reset parameters
integer numImplantEventsReset		!<For creating restart file (debugging tool, see example): number of cascades/Frenkel pairs
integer numHeImplantEventsReset		!<For creating restart file (debugging tool, see example): number of helium implantaion events
double precision elapsedTimeReset	!<For creating restart file (debugging tool, see example): elapsed time
character*20 debugToggle			!<('yes' or 'no') input parameter indicating whether we are restarting from a file
character*50 restartFileName		!<Name of restart file

!array containing distribution of DPA and He implantation rates
integer numImplantDataPoints							!<For non-uniform implantation, number of input data points through-thickness
double precision, allocatable :: implantRateData(:,:)	!<Data containing implantation rates as a function of depth (for non-uniform implantation)

!Strain-related varialbes
integer numDipole										!<Number of dipole tensors that are read in from a file
type(dipoleTensor), allocatable :: dipoleStore(:)		!<Array of dipole tensors for associated defects, size numDipole
character*50 strainFileName								!<File name of strain field intput data
character*50 dipoleFileName								!<File name of dipole tensor input data

contains

end module
