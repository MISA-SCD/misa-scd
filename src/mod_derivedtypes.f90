! $Header: /home/CVS//srscd/src/mod_derivedtypes.f90,v 1.5 2015/10/09 15:36:46 aydunn Exp $
!***************************************************************************************************
!>mod_DerivedTypes contains the variable derived types created for SRSCD.
!!
!!This module contains derived variable types used in SRSCD. They are roughly divided into
!!defect and reaction lists (pointer lists), processor/mesh variables (used for the structure
!!of the code), cascade information such as lists of defects, and material properties/allowed
!!reactions (including defect diffusion and binding energies, allowed clustering reactions, etc).
!!Some comments on the variable types introduced here are below:
!!
!!
!>General comments:
!!
!!Defects and reactions are all treated as pointers, and a list of defects and reactions is created
!!for each volume element. CellNumber referes to the cell WITHIN A PROCESSOR that the defect is located
!!in.
!!
!!processorData gives the processor ID#, the total number of processors, the boundaries of the 
!!global system as well as the local processor's physical boundaries, and the processor numbers
!!of the neighboring processors (in the 'processor mesh').
!!
!!mesh contains the information on the local mesh within a processor. Each element of mesh has coord-
!!inates, a processor number, a material number, a length (cubic elements only), the number of 
!!neighbors in each direction (for nonuniform meshes), and a list of neighboring elements and their
!!processors.
!!
!!mesh and processor are NOT pointers, these are arrays of a variable type created during the mesh
!!initialization step.
!!
!!DiffusionFunction, diffusionSingle, bindingFunction, bindingSingle, and reactionParameters are material
!!classes that are read in from a file and stored in a particular way in order to quickly find
!!reaction rates of all allowed reactions.
!***************************************************************************************************

module DerivedType
	!Defect and reaction variable classes
	!>Type: defect (pointer list).
	!!
	!!Contains the defect type, the number of said defects in a volume element, the cell number, and
	!!a pointer to the next defect in the list.
	type defect
		integer, allocatable :: defectType(:) !<Array containing the number of particles of each defect species in this defect type.
												!!(Example: here, 3	2 0 0 indicates He_3V_2 and 0 0 20 0 indicates SIA_20 (glissile)
		integer num								!<Number of defects of this type inside this volume element
		integer cellNumber						!<Volume element that these defects are located in
		type(defect), pointer :: next			!<Pointer to the next defect in the same volume element (sorted list)
	end type defect
	
	!This is used in updateDefectList() in order to create a list of defects that need to be updated in updateReactionList
	!>Type: defect update tracker (pointer list)
	!!
	!!Contains a list of defects that need to be updated (along with reaction rates) due to the 
	!!reaction that got carried out during the previous step. Has more information than just the 
	!!defect type in order to make the updateReactionList subroutine easier.
	type defectUpdateTracker
		integer, allocatable :: defectType(:) 	!<Type of defect that needs to be updated (see description of defect type for format)
		integer num								!<Number of defects of this type in the volume element
		integer cellNumber						!<Volume element number that this defect is located in
		integer proc							!<Processor ID # that this volume element is located inside
		integer dir								!<If defect diffused from a different volume element, indicates direction from which defect came
		integer neighbor						!<If defect diffused from a different volume element, indicates volume element number from which defect came
		integer cascadeNumber					!<If defect is inside a cascade, indicates cascade ID # that defect is located inside
		type(defectUpdateTracker), pointer :: next	!<Pointer to the next defect that needs to have its reaction rate updated
	end type defectUpdateTracker
	
	!>Type: reaction (pointer list)
	!!
	!!Contains a list of reactions that can be carried out inside a given volume element,
	!!along with the reactants, products, number of reactants, number of products, and 
	!!reaction rate.
	type reaction
		integer numReactants					!<Number of reactants in this reaction
		integer numProducts						!<Number of products in this reaction
		integer, allocatable :: reactants(:,:) 	!<Reactants in this reaction - array size:(numReactants, numSpecies)
		integer, allocatable :: products(:,:)	!<Products in this reaction - array size:(numReactants, numSpecies)
		integer, allocatable :: cellNumber(:) 	!<Cell numbers of defects involved in this reaction. Array size: (numReactants+numProducts). Important for diffusion reactions
		integer, allocatable :: taskid(:) 		!<Processor number of defects involved in this reaction. Array size: (numReactants+numProducts). May not all be the same processor number in the case of diffusion across processor boundaries.
		double precision reactionRate			!<Rate (s^-1) for this reaction to be carried out, based on the defects present in the volume
		type(reaction), pointer :: next			!<Pointer to the next reaction in the list in the same volume element (unsorted list as of now)
	end type reaction
	
	!Variable classes used for meshing and connectivity
	!>Type: processorData (one variable per processor)
	!!
	!! Contains information about this processor (ID number, coordinates of its boundaries)
	!! as well as its neighbors (ID numbers) and the global system (number of processors, coordinates of entire volume's boundaries)
	type processorData
		integer taskid					!<ID number of this processor (0=master)
		integer numtasks				!<Number of processors in this simulation
		double precision globalCoord(6)	!<Global boundaries of system (xmin, xmax, ymin, ymax, zmin, zmax)
		double precision localCoord(6)	!<Local boundaries of this processor (xmin, xmax, ymin, ymax, zmin, zmax)
		integer procNeighbor(6)			!<ID numbers of procesors neighboring this one (left, right, front, back, up, down). Accounts for periodic boundary conditions, if applicable
	end type processorData
	
	!****************
	!2014/06/16: need to add a variable called 'volume' to mesh class. Initially, volume=length**3, but 
	!the value of volume can be increased or decreased by adding or subtracting fine meshes (adaptive meshing scheme)
	!from the coarse mesh
	!
	!2014/07/22: volume added to mesh and boundaryMesh classes. Volume updated when cascades introduced or destroyed
	!****************
	
	!>Type: mesh (local, one array of this type per processor)
	!!
	!!Contains the mesh and connectivity of the volume elements in the local processor
	!!as well as neighboring processors and material ID numbers. Each variable of this type
	!!represents a single volume element, and a single array of this type is constructed
	!!to represent the entire mesh. This has been partially implemented to allow non-uniform
	!!meshes.
	type mesh
		double precision coordinates(3)				!<Coordinates of center of volume element
		double precision length						!<Side length of this volume element (assuming cubic elements)
		double precision volume						!<Volume of this volume element
		double precision strain(6)					!<Strain tensor at the center of this volume element (e11, e22, e33, e12, e23, e13)
		integer proc								!<Processor ID number that this element is located inside
		integer material							!<Material ID number that this element is composed of (currently only set up for one material type)
		integer numNeighbors(6)						!<Number of neighbors in each direction (left, right, etc). Could be not equal to 1 in the case of free surfaces or non-uniform mesh.
		
		!array sizes: neighbors(direction,num) and neighborProcs(direction,num)
		integer, allocatable :: neighbors(:,:)		!<ID number of neighboring volume elements, regardless of if they are in this processor or not. Array size (numNeighbors, 6)
		integer, allocatable :: neighborProcs(:,:)	!<Processor ID numbers of neighboring volume element. Array size (numNeighbors, 6)
	end type mesh
	
	!>Type: boundary mesh (one array of this type per processor)
	!!
	!!Contains information on volume elements located in other processors that share a boundary
	!!with this processor's domain. Includes cell numbers and volumes as well as lists of 
	!!defects present in neighboring processors' boundary elements. Each variable of this type
	!!represents a single volume element and its defects in a neighboring processor, and an array of these is created
	!!to represent the entire boundary.
	type boundaryMesh
		integer proc						!<Processor ID number of this volume element (typically different from the local processor){
		integer material					!<Material ID number of this volume element (currently only set up for one material)
		integer localNeighbor				!<ID number of the volume element in the local mesh that borders this volume element
		double precision length				!<Side length of this volume element (cubic assumption)
		double precision volume				!<Volume of this element
		type(defect), pointer :: defectList	!<List of defects present in this volume element (See defect derived type)
		double precision strain(6)			!<Strain tensor at the center of this volume element
	end type boundaryMesh
	
	!Variable classes used for material inputs
	!List of diffusion rates that are in functional forms
	
	!>Type: diffusion function (list of functional forms for defect diffusivities read in from input file)
	!!
	!!This derived type contains information read in from an input file concerning which defects
	!!have diffusivity given by a functional form, which functional form to use, and parameters
	!!to input into the functional form. One list (array) of this type is created in each processor.
	type diffusionFunction
		integer, allocatable :: defectType(:)			!<Type of defect that is allowed to diffuse using this functional form (1's and 0's only).
		integer, allocatable :: min(:)					!<Minimum defect size allowed to diffuse using this functional form
		integer, allocatable :: max(:)					!<Maximum defect size allowed to diffuse using this functional form (-1 indicates infinity)
		double precision, allocatable :: parameters(:)	!<Parameters to input into this functional form (the function is hard-coded into Defect_Attributes.f90)
		integer functionType							!<ID number of the function type to use to calculate diffusivity
		integer numParam								!<Number of parameters needed to input into this functional form
	end type diffusionFunction
	
	!List of diffusion rates that are in single-defect form
	
	!>Type: diffusion single (list of single migration energies and diffusion prefactors)
	!!
	!!This derived type contains a list of diffusion energies and prefactors for defects 
	!!read in from an input file. It specifies the defect type and the diffusion parameters.
	!!It is used in SRSCD as an array with size given by a separate input, numDiffSingle
	type diffusionSingle
		integer, allocatable :: defectType(:)	!<Type of defect that is allowed to diffuse
		double precision D						!<Diffusion prefactor (nm^2/s)
		double precision Em						!<Migration energy (eV)
	end type diffusionSingle
	
	!List of binding energies that are in single-defect form
	!>Type: binding single (list of binding energies)
	!!
	!!This derived type contains a list of binding energies for defects read in from an 
	!!input file. It contains both the initial defect type as well as the type of 
	!!defect dissociating from it. It is used in SRSCD as an array with size given by a 
	!!separate input, numBindingSingle
	type bindingSingle
		integer, allocatable :: defectType(:)	!<Type of defect that can dissociate
		integer, allocatable :: product(:)		!<Type of defect dissociating away
		double precision :: Eb					!<Binding energy (eV)
	end type
	
	!List of binding energies that are in functional form
	!>Type: binding function (list of functional forms for defect binding energies read in from input file)
	!!
	!!This derived type contains a list of the functional forms used to calculate the binding
	!!energies of various defect types to clusters. It is used in SRSCD as an array with size
	!!given by an input, numBindingFunc.
	type bindingFunction
		integer, allocatable :: defectType(:)	!<Type of cluster that is dissociating (1's and 0's only)
		integer, allocatable :: product(:)		!<Type of defect that dissociates from cluster
		integer, allocatable :: min(:)			!<Minimum cluster size allowed to use this functional form for binding energy (size numSpecies)
		integer, allocatable :: max(:)			!<Maximum cluster size allowed to use this functional form for binding energy (size numSpecies)
		double precision, allocatable :: parameters(:)	!<Parameters to input into this functional form (read in from input file)
		integer functionType					!<ID number of functional form, function is hard coded into Defect_Attributes.f90
		integer numParam						!<Number of parameters needed to input into this functional form
	end type
	
	!Stored information on reactions allowed in a given material; this is read in from file
	!>Type: reaction parameters (list of allowed reactions, input from file)
	!!
	!!This derived type is used in SRSCD as an array with size given by the number of allowed
	!!reactions that can be carried out in the material model chosen. It contains the parameters
	!!as well as an ID number for the function type that determine the reaction rate for 
	!!a given allowed reaction. Functional forms for reaction rates are hard coded in ReactionRates.f90
	type reactionParameters
		integer, allocatable :: reactants(:,:)	!<Defect types for reactants in this reaction
		integer, allocatable :: products(:,:)	!<Defect types for products in this reaction
		integer, allocatable :: min(:)			!<Smallest defect size for reactants (size numSpecies)
		integer, allocatable :: max(:)			!<Largest defect size for reactants (size numSpecies, -1 indicates infinity)
		integer numReactants					!<Number of reactants (0, 1, or 2 typically)
		integer numProducts						!<Number of products
		integer functionType					!<ID number of function type used to calculate reaction rate, functions are hard coded in ReactionRates.f90
	end type reactionParameters
	
	!Dipole tensor derived type (stored for lookup when computing diffusivity)
	!>Type: dipole tensor (stored dipole tensor data for various defect types, input from file)
	!!
	!!This derived type is used to store information on the dipole tensors
	!!that are used in the simulation to calculate modified diffusivity
	!!in the presence of a strain field. It contains the defect type for which 
	!!the dipole tensor refers as well as the dipole tensor itself (in units of eV).
	type dipoleTensor
		integer, allocatable :: min(:)			!<Smallest defect size for this defect type
		integer, allocatable :: max(:)			!<Largest defect size for this defect type (-1 means infinity)
		double precision equilib(6)				!<Equilibrium state dipole tensor(d11, d22, d33, d12, d23, d13)	
		double precision saddle(6)				!<Saddle point (for migration) dipole tensor
	end type
	
	!List of cascades that is chosen from whenever a cascade is added to the simulation
	!>Type: cascade event (pointer list of read-in cascade data from input file)
	!!
	!!This derived type contains a list of the cascades that are read in from a file
	!!including what defects they contain. Cascades are randomly chosen for implantation
	!!when a cascade event is chosen in the Monte Carlo algorithm.
	type cascadeEvent
		integer NumDefectsTotal			!<Number of total defects in the cascade list
		integer numDisplacedAtoms		!<How much this cascade contributes to the DPA (how many lattice atoms are displaced)
		type(cascadeDefect), pointer :: ListOfDefects	!<Pointer pointing to the list of defects in this cascade
		type(cascadeEvent), pointer:: nextCascade		!<Pointer pointing to the next cascade in the list
	end type cascadeEvent
	
	!Separate defect species used to read in defects from the cascade list text file
	!>Type: cascade defect (pointer list of defect types in a cascade)
	!!
	!!This derived type contains a list of defects as well as their locations within a cascade.
	!!This information is read in from a file and stored for use whenever a cascade event is
	!!chosen in the Monte Carlo algorithm.
	type cascadeDefect
		integer, allocatable :: defectType(:) 	!<Type of defect (numSpecies)
		double precision :: coordinates(3) 		!<Location of defect relative to center of cascade (in nm)
		type(cascadeDefect), pointer :: next	!<Pointer pointing to the next defect in the list
	end type cascadeDefect
	
	!cascade structure used to track presence of cascades in adaptive meshing protocol
	!>Type: cascade (pointer list containing defects and reactions for cascades that are currently
	!!present in the system)
	!!
	!!This derived type contains a pointer list of all of the cascades that are currently 
	!!present in the system, including a fine mesh of all of the defect lists inside the
	!!cascade and all of the reaction lists that can occur inside the cascade. This derived
	!!type is essentially creating a new mesh with a new set of defect lists and reaction
	!!lists inside a coarse mesh element.
	type cascade
		integer cellNumber								!<Coarse mesh cell number that this cascade is located inside
		integer cascadeID								!<ID number of this cascade
		type(defect), pointer :: localDefects(:)		!<Array of defect lists (one for each element in the cascade)
		type(cascade), pointer :: next					!<Pointer pointing to the next cascade currently present in the system
		type(cascade), pointer :: prev					!<Pointer pointing to the previous cascade currently present in the system
		type(reaction), pointer :: reactionList(:) 		!<Array of reaction lists (one for each element in the cascade)
		double precision, allocatable :: totalRate(:)	!<Sum of the rates of all reactions in each cascade element
	end type cascade
end module
