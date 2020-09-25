!***************************************************************************************************
!> Subroutine ReadInputs():  - reads in simulation parameters from parameters.txt
!This subroutine reads in all simulation parameters located in parameters.txt as well file names
!for all other input files (defect attributes, mesh, cascades, implantation, etc).
!***************************************************************************************************
subroutine ReadInputs()
	use mod_structures
	use mod_constants
	implicit none

	character(len=20) :: char
	integer :: i
	character(len=100) :: defectFilename                 !<Filename of defect attributes file
	character(len=100) :: meshFilename                   !<Filename of mesh file
	character(len=100) :: cascadeFilename                !<Filename of cascade file
	logical :: alive1, alive2, flag, flag1

	flag= .false.
	!<read in filename of defectFile
	inquire(file='parameters.txt', exist=alive1)
	if(.not. alive1) then
		write(*,*) 'parameters.txt does not exist'
	else
		open(unit=PARAFILE, file='parameters.txt', status='old',action='read')
	end if

	!*******************************************************
	!<read in toogles
	!*******************************************************
	!<set default valuse for toogles
	implantScheme = 'MonteCarlo'
	implantDist = 'uniform'
	grainBoundaryToggle = 'no'
	pointDefectToggle = 'no'
	meshType = "periodic"

	!<read in filename of defect attributes file
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='defectFile') then
			read(PARAFILE,*) defectFilename
			flag=.true.
		end if
	end do
	flag= .false.

	call readDefectAttributes(defectFilename)

	!<read in mesh file
	do  while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='meshFile') then
			read(PARAFILE,*) meshFilename
			flag=.true.
		end if
	end do
	flag=.false.

	!*************************************************
	!<read in mesh file
	inquire(file=meshFilename, exist=alive2)
	if(.not. alive2) then
		write(*,*) 'mesh file does not exist'
	else
		open(MESHFILE, file=meshFilename, status='old', action='read')
	end if

	do while(flag .eqv. .false.)
		read(MESHFILE,*) char
		if(char=='meshType') then
			read(MESHFILE,*) meshType
			flag=.true.
		end if
	end do
	flag=.false.

	do while(flag .eqv. .false.)
		read(MESHFILE,*) char
		if(char=='length') then
			read(MESHFILE,*) meshLength
			flag=.true.
		end if
	end do
	flag=.false.

	do while(flag .eqv. .false.)
		read(MESHFILE,*) char
		if(char=='numx') then
			read(MESHFILE,*) numx
		else if(char=='numy') then
			read(MESHFILE,*) numy
		else if(char=='numz') then
			read(MESHFILE,*) numz
			flag=.true.
		end if
	end do
	flag=.false.

	close(MESHFILE)
	!*************************************************

	!read in irradiation type
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='implantType') then
			read(PARAFILE,*) implantType
			flag=.true.
		end if
	end do
	flag=.false.

	!read in cascade implantation scheme (Monte Carlo or explicit)
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='implantScheme') then
			read(PARAFILE,*) implantScheme
			flag=.true.
		end if
	end do
	flag=.false.

	!Check to make sure that we only are choosing explicit implantation with cascades
	if(implantScheme=='explicit' .AND. implantType=='FrenkelPair') then
		write(*,*) 'error frenkel pairs with explicit implantation'
	endif

	!read in filename of cascade file
	do while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='cascadeFile') then
			read(PARAFILE,*) cascadeFilename
			flag=.true.
		end if
	end do
	flag=.false.

	if(implantType=='Cascade') then
		call readCascadeList(cascadeFilename)
	else if(implantType=='FrenkelPair') then
		numDisplacedAtoms=1				!Frenkel pair implantation, no need to read cascades, one displaced atom per implant event
	else
		write(*,*) 'error implantType'
	end if

	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='meshingType') then
			read(PARAFILE,*) meshingType
			flag=.TRUE.
		end if
	end do
	flag=.false.

	if(implantType /= 'Cascade' .AND. meshingType == 'adaptive') then
		write(*,*) 'error meshingType, it should be nonAdaptive'
	end if

	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='implantDist') then
			read(PARAFILE,*) implantDist
			flag=.TRUE.
		end if
	end do
	flag=.false.

	!read in grain boundary toggle
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='grainBoundaries') then
			read(PARAFILE,*) grainBoundaryToggle
			flag=.true.
		end if
	end do
	flag=.false.

	!read in point defect toggle
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='pointDefect') then
			read(PARAFILE,*) pointDefectToggle
			flag=.true.
		end if
	end do
	flag=.false.

	!*******************************************************
	!<read in simulartion parameters
	!*******************************************************
	!<set default valuse for toogles
	test3			='no'
	tempStore		=273d0
	CuContent		=0.5d-2
	numVac			=0
	numInt			=0
	dpaRate			=1d-4
	totalDPA		=1d-1
	firr			=1d0
	atomSize		=0d0
	lattice			=0.2876d0	!Fe
	burgers			=0.287d0
	reactionRadius	=0.65d0
	agingTime       =0d0	!2019.04.30 Add

	polycrystal			='no'
	meanFreePath		=330000
	dislocationDensity	=0d0
	impurityDensity		=0d0
	max3DInt			=4
	cascadeVolume		=0d0
	numSims				=1
	numGrains			=1
	cascadeReactionLimit=100d0

	annealTemp		=273d0
	annealTime		=0d0
	annealSteps		=0
	annealType		='add'
	annealTempInc	=0d0

	!Read variables in from file
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(PARAFILE,*) char
		if(char=='start') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		flag1=.FALSE.
		do while(flag1 .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='end') then
				flag1=.TRUE.
				flag=.TRUE.
			else if(char=='test3') then
				flag1=.TRUE.
				read(PARAFILE,*) test3
			else if(char=='temperature') then
				flag1=.TRUE.
				read(PARAFILE,*) tempStore
			else if(char=='CuContent') then
				flag1=.TRUE.
				read(PARAFILE,*) CuContent
			else if(char=='numVac') then
				flag1=.TRUE.
				read(PARAFILE,*) numVac
			else if(char=='numInt') then
				flag1=.TRUE.
				read(PARAFILE,*) numInt
			else if(char=='dpaRate') then
				flag1=.TRUE.
				read(PARAFILE,*) dpaRate
			else if(char=='totalDPA') then
				flag1=.TRUE.
				read(PARAFILE,*) totalDPA
			else if(char=='firr') then
				flag1=.TRUE.
				read(PARAFILE,*) firr
			else if(char=='atomSize') then
				flag1=.TRUE.
				read(PARAFILE,*) atomSize
			else if(char=='lattice') then
				flag1=.TRUE.
				read(PARAFILE,*) lattice
			else if(char=='burgers') then
				flag1=.TRUE.
				read(PARAFILE,*) burgers
			else if(char=='reactionRadius') then
				flag1=.TRUE.
				read(PARAFILE,*) reactionRadius
			else if(char=='annealTemp') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTemp
			else if(char=='annealSteps') then
				flag1=.TRUE.
				read(PARAFILE,*) annealSteps
			else if(char=='agingTime') then		!2019.04.30 Add
				flag1=.TRUE.
				read(PARAFILE,*) agingTime
			else if(char=='annealTime') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTime
			else if(char=='annealType') then
				flag1=.TRUE.
				read(PARAFILE,*) annealType
			else if(char=='annealTempInc') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTempInc
			else if(char=='grainSize') then
				flag1=.TRUE.
				read(PARAFILE,*) meanFreePath
			else if(char=='dislocDensity') then
				flag1=.TRUE.
				read(PARAFILE,*) dislocationDensity
			else if(char=='impurityConc') then
				flag1=.TRUE.
				read(PARAFILE,*) impurityDensity
			else if(char=='max3DInt') then
				flag1=.TRUE.
				read(PARAFILE,*) max3DInt
			else if(char=='cascadeVolume') then
				flag1=.TRUE.
				read(PARAFILE,*) cascadeVolume
			else if(char=='cascRxnLimit') then
				flag1=.TRUE.
				read(PARAFILE,*) cascadeReactionLimit
			else if(char=='numSims') then
				flag1=.TRUE.
				read(PARAFILE,*) numSims
			else if(char=='polycrystal') then
				flag1=.TRUE.
				read(PARAFILE,*) polycrystal
			else if(char=='numGrains') then
				flag1=.TRUE.
				read(PARAFILE,*) numGrains
			else
				write(*,*) 'error parameter: ', char
			end if
		end do
		flag1=.FALSE.
	end do
	flag=.FALSE.

	!*******************************************************
	!<read in output parameters
	!*******************************************************
	!<set default valuse for output parameters
	totdatToggle ='yes'
	defectToggle ='no'
	stadatToggle ='no'
	rawdatToggle ='no'
	minSCluster = 10
	minVoid = 10
	minLoop = 10
	minSV = 10

	do while(flag .eqv. .FALSE.)
		read(PARAFILE,*) char
		if(char=='OutputStart') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		flag1=.FALSE.
		do while(flag1 .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='end') then
				flag1=.TRUE.
				flag=.TRUE.
			else if(char=='totdatToggle') then
				flag1=.TRUE.
				read(PARAFILE,*) totdatToggle
			else if(char=='defectToggle') then
				flag1=.TRUE.
				read(PARAFILE,*) defectToggle
			else if(char=='stadatToggle') then
				flag1=.TRUE.
				read(PARAFILE,*) stadatToggle
			else if(char=='rawdatToggle') then
				flag1=.TRUE.
				read(PARAFILE,*) rawdatToggle
			else if(char=='minSCluster') then
				flag1=.TRUE.
				read(PARAFILE,*) minSCluster
			else if(char=='minVoid') then
				flag1=.TRUE.
				read(PARAFILE,*) minVoid
			else if(char=='minLoop') then
				flag1=.TRUE.
				read(PARAFILE,*) minLoop
			else if(char=='minSV') then
				flag1=.TRUE.
				read(PARAFILE,*) minSV
			else
				write(*,*) 'error'
			end if
		end do
		flag1=.FALSE.
	end do
	flag=.FALSE.

	!***********************************************************************
	!if we are using adaptive meshing, read in the adaptive meshing parameters
	!***********************************************************************
	if(meshingType=='adaptive') then
		flag=.FALSE.

		do while(flag .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='fineStart') then
				flag=.TRUE.
			endif
		end do
		flag=.FALSE.

		do while(flag .eqv. .FALSE.)
			flag1=.FALSE.
			do while(flag1 .eqv. .FALSE.)
				read(PARAFILE,*) char
				if(char=='end') then
					flag1=.TRUE.
					flag=.TRUE.
				else if(char=='fineLength') then
					flag1=.TRUE.
					read(PARAFILE,*) fineLength
				else if(char=='numxFine') then
					flag1=.TRUE.
					read(PARAFILE,*) numxCascade
				else if(char=='numyFine') then
					flag1=.TRUE.
					read(PARAFILE,*) numyCascade
				else if(char=='numzFine') then
					flag1=.TRUE.
					read(PARAFILE,*) numzCascade
				else
					write(*,*) 'error'
				end if
			end do
			flag1=.FALSE.
		end do
		flag=.FALSE.

		cascadeElementVol=fineLength**3d0
		numCellsCascade=numxCascade*numyCascade*numzCascade
	end if

	close(PARAFILE)

	!***********************************************************************
	!clustering rate constants
	!***********************************************************************

	omega=(48d0*pi**2/atomSize**2)**(1d0/3d0) 			!clustering rate parameter for spherical clusters
	omegastar=(4*pi*reactionRadius)/atomSize			!clustering rate parameter modifier due to reaction radius
	omega2D=(4d0*pi/(atomSize*burgers))**(1d0/2d0)		!clustering rate parameter for 1D migrating circular clusters
	omega1D=(9d0*pi/(16d0*atomSize))**(1d0/6d0)			!clustering rate parameter for 1D migrating spherical clusters
	omegastar1D=reactionRadius*(pi/atomSize)**(1d0/2d0)
	!omegastar1D=0d0									!clustering rate parameter modifier due to reaction radius
	omegacircle1D=(1d0/burgers)**(1d0/2d0)				!clustering rate parameter for 1D migrating circular clusters

	recombinationCoeff=4d0*pi*(.4466)/atomSize			!from Stoller et al., not used any longer

	if(myProc%taskid==MASTER) write(*,*) 'cascadeReactionLimit', cascadeReactionLimit

end subroutine

!***************************************************************************************************
!>Subroutine readDefectAttributes(filename). Information read in includes:
!   1) Number of allowed defect types
!   2) Formation energies, diffusion prefactors, migration energies and binding energies
!   3) Allowed reactions
!***************************************************************************************************
subroutine readDefectAttributes(filename)
	use mod_constants
	use mod_structures
	implicit none

	character(len=100), intent(in) :: filename
	character(len=20) :: char
	logical :: flag
	integer :: i,j

	flag=.FALSE.
	open(ATTRFILE, file=filename,action='read', status='old')

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='species') then
			flag=.TRUE.
			read(ATTRFILE,*) numSpecies	!< numSpecies = 4
		end if
	end do
	flag=.FALSE.

	!*******************************************************
	!<Read in formation energies
	!*******************************************************
	do while(flag .eqv. .false.)
		read(ATTRFILE,*) char
		if(char=='formationEnergies') then
			flag=.true.
		end if
	end do
	flag=.false.

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(ATTRFILE,*) numSingleForm
		end if
	end do
	flag=.FALSE.

	allocate(FormSingle(numSingleForm))
	do i=1,numSingleForm
		allocate(FormSingle(i)%defectType(numSpecies))
		read(ATTRFILE,*) (FormSingle(i)%defectType(j),j=1,numSpecies)
		read(ATTRFILE,*) char, FormSingle(i)%Ef
	end do

	!*******************************************************
	!<Read in diffusion prefactors and migration energies
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='diffusionPrefactors') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(ATTRFILE,*) numSingleDiff
		end if
	end do
	flag=.FALSE.

	allocate(DiffSingle(numSingleDiff))
	do i=1,numSingleDiff
		allocate(DiffSingle(i)%defectType(numSpecies))
		read(ATTRFILE,*) (DiffSingle(i)%defectType(j),j=1,numSpecies)
		read(ATTRFILE,*) char, DiffSingle(i)%D, char, DiffSingle(i)%Em
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='numFunction') then
			flag=.TRUE.
			read(ATTRFILE,*) numFuncDiff
		end if
	end do
	flag=.FALSE.

	allocate(DiffFunc(numFuncDiff))
	do i=1,numFuncDiff
		allocate(DiffFunc(i)%defectType(numSpecies))
		read(ATTRFILE,*) (DiffFunc(i)%defectType(j),j=1,numSpecies)	!< read in defectTypes
		allocate(DiffFunc(i)%min(numSpecies))
		allocate(DiffFunc(i)%max(numSpecies))
		read(ATTRFILE,*) char, (DiffFunc(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (DiffFunc(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, DiffFunc(i)%functionType
		read(ATTRFILE,*) char, DiffFunc(i)%numParam
		allocate(DiffFunc(i)%parameters(DiffFunc(i)%numParam))
		if(DiffFunc(i)%numParam /= 0) then
			read(ATTRFILE,*) (DiffFunc(i)%parameters(j),j=1,DiffFunc(i)%numParam)
		end if
	end do

	!*******************************************************
	!<Read in binding energies
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='bindingEnergies') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(ATTRFILE,*) numSingleBind
		end if
	end do
	flag=.FALSE.

	allocate(BindSingle(numSingleBind))
	do i=1,numSingleBind
		allocate(BindSingle(i)%defectType(numSpecies))
		allocate(BindSingle(i)%product(numSpecies))
		read(ATTRFILE,*) (BindSingle(i)%defectType(j),j=1,numSpecies),(BindSingle(i)%product(j),j=1,numSpecies)
		read(ATTRFILE,*) char, BindSingle(i)%Eb
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='numFunction') then
			flag=.TRUE.
			read(ATTRFILE,*) numFuncBind
		end if
	end do
	flag=.FALSE.

	allocate(BindFunc(numFuncBind))
	do i=1,numFuncBind
		allocate(BindFunc(i)%defectType(numSpecies))
		allocate(BindFunc(i)%product(numSpecies))
		read(ATTRFILE,*) (BindFunc(i)%defectType(j),j=1,numSpecies),(BindFunc(i)%product(j),j=1,numSpecies)
		allocate(BindFunc(i)%min(numSpecies))
		allocate(BindFunc(i)%max(numSpecies))
		read(ATTRFILE,*) char, (BindFunc(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (BindFunc(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, BindFunc(i)%functionType
		read(ATTRFILE,*) char, BindFunc(i)%numParam
		allocate(BindFunc(i)%parameters(BindFunc(i)%numParam))
		if(BindFunc(i)%numParam /= 0) then
			read(ATTRFILE,*) (BindFunc(i)%parameters(j),j=1,BindFunc(i)%numParam)
		end if
	end do

	!*******************************************************
	!<Read in reactions
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='dissociation') then
			flag=.TRUE.
			read(ATTRFILE,*) numDissocReac
		end if
	end do
	flag=.FALSE.

	allocate(DissocReactions(numDissocReac))
	do i=1,numDissocReac
		DissocReactions(i)%numReactants=1
		DissocReactions(i)%numProducts=1
		allocate(DissocReactions(i)%reactants(numSpecies,DissocReactions(i)%numReactants))
		allocate(DissocReactions(i)%products(numSpecies,DissocReactions(i)%numProducts))
		read(ATTRFILE,*) (DissocReactions(i)%reactants(j,1),j=1,numSpecies),&
				(DissocReactions(i)%products(j,1),j=1,numSpecies)	!< read in defectType
		allocate(DissocReactions(i)%min(numSpecies))
		allocate(DissocReactions(i)%max(numSpecies))
		read(ATTRFILE,*) char, (DissocReactions(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (DissocReactions(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, DissocReactions(i)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='diffusion') then
			flag=.TRUE.
			read(ATTRFILE,*) numDiffReac
		end if
	end do
	flag=.FALSE.

	allocate(DiffReactions(numDiffReac))
	do i=1,numDiffReac
		DiffReactions(i)%numReactants=1
		DiffReactions(i)%numProducts=1
		allocate(DiffReactions(i)%reactants(numSpecies,DiffReactions(i)%numReactants))
		allocate(DiffReactions(i)%products(numSpecies,DiffReactions(i)%numProducts))
		read(ATTRFILE,*) (DiffReactions(i)%reactants(j,1),j=1,numSpecies),&
				(DiffReactions(i)%products(j,1),j=1,numSpecies)
		allocate(DiffReactions(i)%min(numSpecies))
		allocate(DiffReactions(i)%max(numSpecies))
		read(ATTRFILE,*) char, (DiffReactions(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (DiffReactions(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, DiffReactions(i)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='sinkRemoval') then
			flag=.TRUE.
			read(ATTRFILE,*) numSinkReac
		end if
	end do
	flag=.FALSE.

	allocate(SinkReactions(numSinkReac))

	do i=1,numSinkReac
		SinkReactions(i)%numReactants=1
		SinkReactions(i)%numProducts=0
		allocate(SinkReactions(i)%reactants(numSpecies,SinkReactions(i)%numReactants))
		read(ATTRFILE,*) (SinkReactions(i)%reactants(j,1),j=1,numSpecies)
		allocate(SinkReactions(i)%min(numSpecies))
		allocate(SinkReactions(i)%max(numSpecies))
		read(ATTRFILE,*) char, (SinkReactions(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (SinkReactions(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, SinkReactions(i)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='impurityTrapping') then
			flag=.TRUE.
			read(ATTRFILE,*) numImpurityReac
		end if
	end do
	flag=.FALSE.

	allocate(ImpurityReactions(numImpurityReac))
	do i=1,numImpurityReac
		ImpurityReactions(i)%numReactants=1
		ImpurityReactions(i)%numProducts=1
		allocate(ImpurityReactions(i)%reactants(numSpecies,ImpurityReactions(i)%numReactants))
		allocate(ImpurityReactions(i)%products(numSpecies,ImpurityReactions(i)%numProducts))
		read(ATTRFILE,*) (ImpurityReactions(i)%reactants(j,1),j=1,numSpecies), &
				(ImpurityReactions(i)%products(j,1),j=1,numSpecies)
		allocate(ImpurityReactions(i)%min(numSpecies))
		allocate(ImpurityReactions(i)%max(numSpecies))
		read(ATTRFILE,*) char, (ImpurityReactions(i)%min(j),j=1,numSpecies)
		read(ATTRFILE,*) char, (ImpurityReactions(i)%max(j),j=1,numSpecies)
		read(ATTRFILE,*) char, ImpurityReactions(i)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='clustering') then
			flag=.TRUE.
			read(ATTRFILE,*) numClusterReac
		end if
	end  do
	flag=.FALSE.

	allocate(ClusterReactions(numClusterReac))
	do i=1,numClusterReac
		ClusterReactions(i)%numReactants=2
		ClusterReactions(i)%numProducts=1
		allocate(ClusterReactions(i)%reactants(numSpecies,ClusterReactions(i)%numReactants))
		allocate(ClusterReactions(i)%products(numSpecies,ClusterReactions(i)%numProducts))
		read(ATTRFILE,*) (ClusterReactions(i)%reactants(j,1),j=1,numSpecies),&
				(ClusterReactions(i)%reactants(j,2),j=1,numSpecies)
		allocate(ClusterReactions(i)%min(numSpecies*ClusterReactions(i)%numReactants))
		allocate(ClusterReactions(i)%max(numSpecies*ClusterReactions(i)%numReactants))
		read(ATTRFILE,*) char,(ClusterReactions(i)%min(j),j=1,numSpecies*ClusterReactions(i)%numReactants)
		read(ATTRFILE,*) char,(ClusterReactions(i)%max(j),j=1,numSpecies*ClusterReactions(i)%numReactants)
		do j=1,numSpecies
			ClusterReactions(i)%products(j,1)=ClusterReactions(i)%reactants(j,1)+&
					ClusterReactions(i)%reactants(j,2)
		end do
		read(ATTRFILE,*) char, ClusterReactions(i)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(ATTRFILE,*) char
		if(char=='Implantation') then
			flag=.TRUE.
			read(ATTRFILE,*) numImplantReac
		end if
	end do
	flag=.FALSE.

	allocate(ImplantReactions(numImplantReac))
	do i=1,numImplantReac
		if(i==1) then	!Read in Frenkel pair reaction parameters
			do while(flag .eqv. .FALSE.)
				read(ATTRFILE,*) char
				if(char=='FrenkelPair') then
					flag=.TRUE.
				end if
			end do
			flag=.FALSE.

			ImplantReactions(i)%numReactants=0
			ImplantReactions(i)%numProducts=2
			allocate(ImplantReactions(i)%products(numSpecies,ImplantReactions(i)%numProducts))
			read(ATTRFILE,*) (ImplantReactions(i)%products(j,1),j=1,numSpecies),&
					(ImplantReactions(i)%products(j,2),j=1,numSpecies)
			read(ATTRFILE,*) char, ImplantReactions(i)%functionType

		else if(i==2) then !read in cascade reaction parameters
			do while(flag .eqv. .FALSE.)
				read(ATTRFILE,*) char
				if(char=='Cascade') then
					flag=.TRUE.
				end if
			end do
			flag=.FALSE.

			ImplantReactions(i)%numReactants=-10
			ImplantReactions(i)%numProducts=0
			read(ATTRFILE,*) char, ImplantReactions(i)%functionType
		else
			write(*,*) 'error numImplantReac'
		end if
	end do

	close(ATTRFILE)

end subroutine

!***************************************************************************************************
!> Subroutine readCascadeList(filename): reads defects in a cascade file
!This subroutine reads a list of cascades from a cascade file and stored in derived type cascade.
!***************************************************************************************************
subroutine readCascadeList(filename)
	use mod_constants
	use mod_structures
	implicit none

	character(len=100), intent(in) :: filename
	character(len=20) :: char
	type(cascadeEvent), pointer :: cascadeCurrent
	type(cascadeDefect), pointer :: defectCurrent
	integer :: i, numDefects, j, k
	logical :: flag

	open(CASFILE, file=filename, status='old', action='read')

	allocate(cascadeList)
	cascadeCurrent=>cascadeList
	read(CASFILE,*) numDisplacedAtoms
	read(CASFILE,*)
	read(CASFILE,*) numCascades

	do i=1,numCascades
		read(CASFILE,*)
		read(CASFILE,*) numDefects
		read(CASFILE,*) cascadeCurrent%numDisplacedAtoms
		cascadeCurrent%NumDefectsTotal=numDefects
		allocate(cascadeCurrent%ListOfDefects)
		defectCurrent=>cascadeCurrent%ListOfDefects
		nullify(cascadeCurrent%nextCascade)
		nullify(defectCurrent%next)

		do j=1,numDefects
			allocate(defectCurrent%defectType(numSpecies))
			read(CASFILE,*) (defectCurrent%defectType(k),k=1,numSpecies)
			read(CASFILE,*) (defectCurrent%coordinates(k), k=1,3)

			if(j /= numDefects) then
				allocate(defectCurrent%next)
				nullify(defectCurrent%next%next)
				defectCurrent=>defectCurrent%next
			end if
		end do

		if(i /= numCascades) then
			allocate(cascadeCurrent%nextCascade)
			cascadeCurrent=>cascadeCurrent%nextCascade
		end if
	end do
	nullify(defectCurrent)
	nullify(cascadeCurrent)

	close(CASFILE)

end subroutine
