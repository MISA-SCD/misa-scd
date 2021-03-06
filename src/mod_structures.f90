!***************************************************************************************************
!>Module mod_structures: contains all derived structures created for MISA-SCD1.0.
!>Information read in from input file:
!										Formation energies: formation
!										Diffusivities: diffusion, diffusionFunction
!										Binding energies: binding, bindingFunction
!										Reactions: reactionParameters)
!										Cascade: cascade
!>Information about processors and meshes (processor, mesh, boundaryMesh)
!>Information about defects and reactions (defect, defectUpdateTracker, reaction)
!>Information about cascade (cascadeEvent，cascadeDefect, cascade)
!***************************************************************************************************
module mod_structures
	implicit none
	!**********************************************************************
	!>Constant type: list of defect attributes read in from an input file.
	!**********************************************************************
	type :: formationSingle
		integer, allocatable :: defectType(:)			!<Type of defect of the point defect
		double precision :: Ef							!<Binding energy (eV)
	end type formationSingle

	type :: diffusionSingle
		integer, allocatable :: defectType(:)			!<Type of defect that is allowed to diffuse
		double precision :: D0							!<Diffusion prefactor (nm^2/s)
		double precision :: Em							!<Migration energy (eV)
	end type diffusionSingle

	type :: diffusionFunction
		integer, allocatable :: defectType(:)			!<Type of defect that is allowed to diffuse using this functional form (1's and 0's only).
		integer, allocatable :: min(:)					!<Minimum defect size allowed to diffuse using this functional form
		integer, allocatable :: max(:)					!<Maximum defect size allowed to diffuse using this functional form (-1 indicates infinity)
		double precision, allocatable :: parameters(:)	!<Parameters to input into this functional form (the function is hard-coded into Defect_Attributes.f90)
		integer :: fType								!<ID number of the function type to use to calculate diffusivity
		integer :: numParam								!<Number of parameters needed to input into this functional form
	end type diffusionFunction

	type :: bindingSingle
		integer, allocatable :: defectType(:)			!<Type of defect that can dissociate
		integer, allocatable :: product(:)				!<Type of defect dissociating away (point defect)
		double precision :: Eb							!<Binding energy (eV)
	end type bindingSingle

	type :: bindingFunction
		integer, allocatable :: defectType(:)			!<Type of cluster that is dissociating (1's and 0's only)
		integer, allocatable :: product(:)				!<Type of defect that dissociates from cluster
		integer, allocatable :: min(:)					!<Minimum cluster size allowed to use this functional form for binding energy (size SPECIES)
		integer, allocatable :: max(:)					!<Maximum cluster size allowed to use this functional form for binding energy (size SPECIES)
		double precision, allocatable :: parameters(:)	!<Parameters to input into this functional form (read in from input file)
		integer :: fType								!<ID number of functional form, function is hard coded into Defect_Attributes.f90
		integer :: numParam								!<Number of parameters needed to input into this functional form
	end type bindingFunction

	type :: reactionParameters
		integer, allocatable :: reactants(:,:)			!<Defect types for reactants in this reaction--array size:(SPECIES, numReactants)
		integer, allocatable :: products(:,:)			!<Defect types for products in this reaction--array size:(SPECIES, numProducts)
		integer, allocatable :: min(:)					!<Smallest defect size for reactants (size SPECIES)
		integer, allocatable :: max(:)					!<Largest defect size for reactants (size SPECIES, -1 indicates infinity)
		integer :: numReactants							!<Number of reactants (0, 1, or 2 typically)
		integer :: numProducts							!<Number of products
		integer :: fType								!<ID number of function type used to calculate reaction rate, functions are hard coded in ReactionRates.f90
	end type reactionParameters

	!**********************************************************************
	!>Structures related to processes and meshes
	!**********************************************************************
	type :: processorData
		integer :: taskid								!<ID number of this processor (0=master)
		integer :: numtasks								!<Number of processors in this simulation
		integer :: coords(3)
		double precision :: globalCoord(6)				!<Global boundaries of system (xmin, xmax, ymin, ymax, zmin, zmax)
		double precision :: localCoord(6)				!<Local boundaries of this processor (xmin, xmax, ymin, ymax, zmin, zmax)
		integer :: procNeighbor(6)						!<ID numbers of procesors neighboring this one (right, left, front, back, up, down). Accounts for periodic boundary conditions, if applicable
	end type processorData

	type :: mesh
		double precision :: coordinates(3)				!<Coordinates of center of volume element
		double precision :: length						!<Side length of this volume element (assuming cubic elements)
		double precision :: volume						!<Volume of this volume element
		integer :: globalCell							!<Global ID of the mesh
		integer :: proc									!<Processor ID number that this element is located inside
		integer :: material								!<Material ID number that this element is composed of (currently only set up for one material type)
		integer :: numNeighbors(6)						!<Number of neighbors in each direction (left, right, etc). Could be not equal to 1 in the case of free surfaces or non-uniform mesh.
		integer :: neighbors(6)							!<ID of neighboring volume elements, regardless of if they are in this processor or not. Array size (6)
		integer :: neighborProcs(6)						!<Processor ID numbers of neighboring volume element. Array size (numNeighbors,6)
	end type mesh

	type :: ghostMesh
		integer :: proc									!<Processor ID number of this volume element (typically different from the local processor){
		integer :: material								!<Material ID number of this volume element (currently only set up for one material)
		integer :: localNeighbor						!<ID number of the volume element in the local mesh that borders this volume element
		double precision :: length						!<Side length of this volume element (cubic assumption)
		double precision :: volume						!<Volume of this element
		type(defect), pointer :: defectList				!<List of defects present in this volume element (See defect derived type)
	end type ghostMesh

	!**********************************************************************
	!>defect and reaction
	!**********************************************************************
	type :: defect
		integer, allocatable :: defectType(:) 			!<Array containing the number of particles of each defect species in this defect type.
		integer :: num									!<Number of defects of this type inside this volume element
		integer :: cellNumber							!<Volume element that these defects are located in
		type(defect), pointer :: next					!<Pointer to the next defect in the same volume element (sorted list)
	end type defect

	type :: defectUpdateTracker
		integer, allocatable :: defectType(:) 			!<Type of defect that needs to be updated (see description of defect type for format)
		integer :: num									!<Number of defects of this type in the volume element
		integer :: cellNumber							!<Volume element number that this defect is located in
		integer :: proc									!<Processor ID # that this volume element is located inside
		integer :: dir									!<If defect diffused from a different volume element, indicates direction from which defect came
		integer :: neighbor								!<If defect diffused from a different volume element, indicates volume element number from which defect came
		integer :: cascadeNumber						!<If defect is inside a cascade, indicates cascade ID # that defect is located inside
		type(defectUpdateTracker), pointer :: next		!<Pointer to the next defect that needs to have its reaction rate updated
	end type defectUpdateTracker

	type :: reaction
		integer :: numReactants							!<Number of reactants in this reaction
		integer :: numProducts							!<Number of products in this reaction
		integer, allocatable :: reactants(:,:) 			!<Reactants in this reaction - array size:(SPECIES, numReactants)
		integer, allocatable :: products(:,:)			!<Products in this reaction - array size:(SPECIES, numReactants)
		integer, allocatable :: cellNumber(:) 			!<Cell numbers of defects involved in this reaction. Array size: (numReactants+numProducts). Important for diffusion reactions
		integer, allocatable :: taskid(:) 				!<Processor number of defects involved in this reaction. Array size: (numReactants+numProducts). May not all be the same processor number in the case of diffusion across processor boundaries.
		double precision :: reactionRate				!<Rate (s^-1) for this reaction to be carried out, based on the defects present in the volume
		type(reaction), pointer :: next					!<Pointer to the next reaction in the list in the same volume element (unsorted list as of now)
	end type reaction

	type :: cascadeEvent
		integer :: numDefectsTotal						!<Number of total defects in the cascade list
		integer :: numDisplacedAtoms					!<How much this cascade contributes to the DPA (how many lattice atoms are displaced)
		type(cascadeDefect), pointer :: ListOfDefects	!<Pointer pointing to the list of defects in this cascade
		type(cascadeEvent), pointer:: next				!<Pointer pointing to the next cascade in the list
	end type cascadeEvent

	type :: cascadeDefect
		integer, allocatable :: defectType(:) 			!<Type of defect (SPECIES)
		double precision :: coordinates(3) 				!<Location of defect relative to center of cascade (in nm)
		type(cascadeDefect), pointer :: next			!<Pointer pointing to the next defect in the list
	end type cascadeDefect

	type :: cascade
		integer :: cellNumber							!<Coarse mesh cell number that this cascade is located inside
		integer :: cascadeID							!<ID number of this cascade
		type(defect), pointer :: localDefects(:)		!<Array of defect lists (one for each element in the cascade)
		type(cascade), pointer :: next					!<Pointer pointing to the next cascade currently present in the system
		type(cascade), pointer :: prev					!<Pointer pointing to the previous cascade currently present in the system
		type(reaction), pointer :: reactionList(:) 		!<Array of reaction lists (one for each element in the cascade)
		double precision, allocatable :: totalRate(:)	!<Sum of the rates of all reactions in each cascade element
	end type cascade

	!<PKA spectrum
	type :: cpdf_t
		integer :: size									!number of datas
		double precision, allocatable :: energy(:)		!<PKA energy
		double precision, allocatable :: cpdf(:)		!<Probability between 0 and 1
	end type cpdf_t

	!<cascade fiels
	type :: cascadeFileList		!<List of cascades that are read in from a number of cascade files
		double precision :: temperature					!<temperature of the cascade collision
		double precision :: PKAenergy					!<PKA energy of the cascade collision
		double precision :: averDisAtoms				!<average number of atoms displaced per cascade in this cascade file
		integer :: numCascades							!<number of cascade in this cascade file
		type(cascadeEvent), pointer :: listCascades		!<list of cascades in this cascade file
	end type cascadeFileList

end module
