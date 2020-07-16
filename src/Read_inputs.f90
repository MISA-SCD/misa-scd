!***************************************************************************************************
!> Subroutine ReadInputs():  - reads in simulation parameters from parameters.txt
!This subroutine reads in all simulation parameters located in parameters.txt as well file names
!for all other input files (defect attributes, mesh, cascades, implantation, etc).
!***************************************************************************************************
subroutine ReadInputs()
	use DerivedType
	use mod_constants
	implicit none

	character(len=20) :: char
	integer :: i
	character(len=50) :: defectFilename                 !<Filename of defect attributes file
	character(len=50) :: meshFilename                   !<Filename of mesh file
	character(len=50) :: cascadeFilename                !<Filename of cascade file
	logical alive1, alive2, flag, flag1

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
			read(MESHFILE,*) numx
		else if(char=='numz') then
			read(MESHFILE,*) numx
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
			read(PARAFILE,*) cascadeFile
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
	totdatToggle		='yes'
	rawdatToggle		='no'
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
	use DerivedType
	implicit none

	character(len=50), intent(in) :: filename
	character(len=20) :: char
	logical flag
	integer i, matNum

	flag=.FALSE.
	numMaterials=1
	matNum=1

	!Allocate all of the counters for the number of reactions, etc...
	allocate(numSingleForm(numMaterials))
	allocate(numSingleDiff(numMaterials))
	allocate(numFuncDiff(numMaterials))
	allocate(numSingleBind(numMaterials))
	allocate(numFuncBind(numMaterials))

	allocate(numDiffReac(numMaterials))
	allocate(numClusterReac(numMaterials))
	allocate(numSinkReac(numMaterials))
	allocate(numImplantReac(numMaterials))
	allocate(numDissocReac(numMaterials))
	allocate(numImpurityReac(numMaterials))

	open(DEFFILE, file=filename,action='read', status='old')

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='species') then
			flag=.TRUE.
			read(DEFFILE,*) numSpecies	!< numSpecies = 4
		end if
	end do
	flag=.FALSE.

	!*******************************************************
	!<Read in formation energies
	!*******************************************************
	do while(flag .eqv. .false.)
		read(DEFFILE,*) char
		if(char=='formationEnergies') then
			flag=.true.
		end if
	end do
	flag=.false.

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(DEFFILE,*) numSingleForm(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(FormSingle(numSingleForm(matNum),numMaterials))
	do i=1,numSingleForm(matNum)
		allocate(FormSingle(i,matNum)%defectType(numSpecies))
		read(DEFFILE,*) (FormSingle(i,matNum)%defectType(j),j=1,numSpecies)
		read(DEFFILE,*) char, FormSingle(i,matNum)%Ef
	end do

	!*******************************************************
	!<Read in diffusion prefactors and migration energies
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='diffusionPrefactors') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(DEFFILE,*) numSingleDiff(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(DiffSingle(numSingleDiff(matNum),numMaterials))
	do i=1,numSingleDiff(matNum)
		allocate(DiffSingle(i,matNum)%defectType(numSpecies))
		read(DEFFILE,*) (DiffSingle(i,matNum)%defectType(j),j=1,numSpecies)
		read(DEFFILE,*) char, DiffSingle(i,matNum)%D, char, DiffSingle(i,matNum)%Em
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='numFunction') then
			flag=.TRUE.
			read(DEFFILE,*) numFuncDiff(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(DiffFunc(numFuncDiff(matNum),numMaterials))
	do i=1,numFuncDiff(matNum)
		allocate(DiffFunc(i,matNum)%defectType(numSpecies))
		read(DEFFILE,*) (DiffFunc(i,matNum)%defectType(j),j=1,numSpecies)	!< read in defectTypes
		allocate(DiffFunc(i,matNum)%min(numSpecies))
		allocate(DiffFunc(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (DiffFunc(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (DiffFunc(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, DiffFunc(i,matNum)%functionType
		read(DEFFILE,*) char, DiffFunc(i,matNum)%numParam
		allocate(DiffFunc(i,matNum)%parameters(DiffFunc(i,matNum)%numParam))
		if(DiffFunc(i,matNum)%numParam /= 0) then
			read(DEFFILE,*) (DiffFunc(i,matNum)%parameters(j),j=1,DiffFunc(i,matNum)%numParam)
		end if
	end do

	!*******************************************************
	!<Read in binding energies
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='bindingEnergies') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='numSingle') then
			flag=.TRUE.
			read(DEFFILE,*) numSingleBind(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(BindSingle(numSingleBind(matNum),numMaterials))
	do i=1,numSingleBind(matNum)
		allocate(BindSingle(i,matNum)%defectType(numSpecies))
		allocate(BindSingle(i,matNum)%product(numSpecies))
		read(DEFFILE,*) (BindSingle(i,matNum)%defectType(j),j=1,numSpecies),(BindSingle(i,matNum)%product(j),j=1,numSpecies)
		read(DEFFILE,*) char, BindSingle(i,matNum)%Eb
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='numFunction') then
			flag=.TRUE.
			read(DEFFILE,*) numFuncBind(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(BindFunc(numFuncBind(matNum),numMaterials))
	do i=1,numFuncBind(matNum)
		allocate(BindFunc(i,matNum)%defectType(numSpecies))
		allocate(BindFunc(i,matNum)%product(numSpecies))
		read(DEFFILE,*) (BindFunc(i,matNum)%defectType(j),j=1,numSpecies),(BindFunc(i,matNum)%product(j),j=1,numSpecies)
		allocate(BindFunc(i,matNum)%min(numSpecies))
		allocate(BindFunc(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (BindFunc(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (BindFunc(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, BindFunc(i,matNum)%functionType
		read(DEFFILE,*) char, BindFunc(i,matNum)%numParam
		allocate(BindFunc(i,matNum)%parameters(BindFunc(i,matNum)%numParam))
		if(BindFunc(i,matNum)%numParam /= 0) then
			read(DEFFILE,*) (BindFunc(i,matNum)%parameters(j),j=1,BindFunc(i,matNum)%numParam)
		end if
	end do

	!*******************************************************
	!<Read in reactions
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='dissociation') then
			flag=.TRUE.
			read(DEFFILE,*) numDissocReac(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(DissocReactions(numDissocReac(matNum),numMaterials))
	do i=1,numDissocReac(matNum)
		DissocReactions(i,matNum)%numReactants=1
		DissocReactions(i,matNum)%numProducts=1
		allocate(DissocReactions(i,matNum)%reactants(numSpecies,DissocReactions(i,matNum)%numReactants))
		allocate(DissocReactions(i,matNum)%products(numSpecies,DissocReactions(i,matNum)%numProducts))
		read(DEFFILE,*) (DissocReactions(i,matNum)%reactants(j,1),j=1,numSpecies),&
				(DissocReactions(i,matNum)%products(j,1),j=1,numSpecies)	!< read in defectType
		allocate(DissocReactions(i,matNum)%min(numSpecies))
		allocate(DissocReactions(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (DissocReactions(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (DissocReactions(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, DissocReactions(i,matNum)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='diffusion') then
			flag=.TRUE.
			read(DEFFILE,*) numDiffReac(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(DiffReactions(numDiffReac(matNum),numMaterials))
	do i=1,numDiffReac(matNum)
		DiffReactions(i,matNum)%numReactants=1
		DiffReactions(i,matNum)%numProducts=1
		allocate(DiffReactions(i,matNum)%reactants(numSpecies,DiffReactions(i,matNum)%numReactants))
		allocate(DiffReactions(i,matNum)%products(numSpecies,DiffReactions(i,matNum)%numProducts))
		read(DEFFILE,*) (DiffReactions(i,matNum)%reactants(j,1),j=1,numSpecies),&
				(DiffReactions(i,matNum)%products(j,1),j=1,numSpecies)
		allocate(DiffReactions(i,matNum)%min(numSpecies))
		allocate(DiffReactions(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (DiffReactions(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (DiffReactions(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, DiffReactions(i,matNum)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='sinkRemoval') then
			flag=.TRUE.
			read(DEFFILE,*) numSinkReac(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(SinkReactions(numSinkReac(matNum),numMaterials))

	do i=1,numSinkReac(matNum)
		SinkReactions(i,matNum)%numReactants=1
		SinkReactions(i,matNum)%numProducts=0
		allocate(SinkReactions(i,matNum)%reactants(numSpecies,SinkReactions(i,matNum)%numReactants))
		read(DEFFILE,*) (SinkReactions(i,matNum)%reactants(j,1),j=1,numSpecies)
		allocate(SinkReactions(i,matNum)%min(numSpecies))
		allocate(SinkReactions(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (SinkReactions(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (SinkReactions(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, SinkReactions(i,matNum)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='impurityTrapping') then
			flag=.TRUE.
			read(DEFFILE,*) numImpurityReac(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(ImpurityReactions(numImpurityReac(matNum),numMaterials))
	do i=1,numImpurityReac(matNum)
		ImpurityReactions(i,matNum)%numReactants=1
		ImpurityReactions(i,matNum)%numProducts=1
		allocate(ImpurityReactions(i,matNum)%reactants(numSpecies,ImpurityReactions(i,matNum)%numReactants))
		allocate(ImpurityReactions(i,matNum)%products(numSpecies,ImpurityReactions(i,matNum)%numProducts))
		read(DEFFILE,*) (ImpurityReactions(i,matNum)%reactants(j,1),j=1,numSpecies), &
				(ImpurityReactions(i,matNum)%products(j,1),j=1,numSpecies)
		allocate(ImpurityReactions(i,matNum)%min(numSpecies))
		allocate(ImpurityReactions(i,matNum)%max(numSpecies))
		read(DEFFILE,*) char, (ImpurityReactions(i,matNum)%min(j),j=1,numSpecies)
		read(DEFFILE,*) char, (ImpurityReactions(i,matNum)%max(j),j=1,numSpecies)
		read(DEFFILE,*) char, ImpurityReactions(i,matNum)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='clustering') then
			flag=.TRUE.
			read(DEFFILE,*) numClusterReac(matNum)
		end if
	end  do
	flag=.FALSE.

	allocate(ClusterReactions(numClusterReac(matNum),numMaterials))
	do i=1,numClusterReac(matNum)
		ClusterReactions(i,matNum)%numReactants=2
		ClusterReactions(i,matNum)%numProducts=1
		allocate(ClusterReactions(i,matNum)%reactants(numSpecies,ClusterReactions(i,matNum)%numReactants))
		allocate(ClusterReactions(i,matNum)%products(numSpecies,ClusterReactions(i,matNum)%numProducts))
		read(DEFFILE,*) (ClusterReactions(i,matNum)%reactants(j,1),j=1,numSpecies),&
				(ClusterReactions(i,matNum)%reactants(j,2),j=1,numSpecies)
		allocate(ClusterReactions(i,matNum)%min(numSpecies*ClusterReactions(i,matNum)%numReactants))
		allocate(ClusterReactions(i,matNum)%max(numSpecies*ClusterReactions(i,matNum)%numReactants))
		read(DEFFILE,*) char,(ClusterReactions(i,matNum)%min(j),i=1,numSpecies*ClusterReactions(i,matNum)%numReactants)
		read(DEFFILE,*) char,(ClusterReactions(i,matNum)%max(j),i=1,numSpecies*ClusterReactions(i,matNum)%numReactants)
		do j=1,numSpecies
			ClusterReactions(i,matNum)%products(j,1)=ClusterReactions(i,matNum)%reactants(j,1)+&
					ClusterReactions(i,matNum)%reactants(j,2)
		end do
		read(DEFFILE,*) char, ClusterReactions(i,matNum)%functionType
	end do

	do while(flag .eqv. .FALSE.)
		read(DEFFILE,*) char
		if(char=='Implantation') then
			flag=.TRUE.
			read(DEFFILE,*) numImplantReac(matNum)
		end if
	end do
	flag=.FALSE.

	allocate(ImplantReactions(numImplantReac(matNum),numMaterials))
	do i=1,numImplantReac(matNum)
		if(i==1) then	!Read in Frenkel pair reaction parameters
			do while(flag .eqv. .FALSE.)
				read(DEFFILE,*) char
				if(char=='FrenkelPair') then
					flag=.TRUE.
				end if
			end do
			flag=.FALSE.

			ImplantReactions(i,matNum)%numReactants=0
			ImplantReactions(i,matNum)%numProducts=2
			allocate(ImplantReactions(i,matNum)%products(numSpecies,ImplantReactions(i,matNum)%numProducts))
			read(DEFFILE,*) (ImplantReactions(i,matNum)%products(j,1),j=1,numSpecies),&
					(ImplantReactions(i,matNum)%products(j,2),j=1,numSpecies)
			read(DEFFILE,*) char, ImplantReactions(i,matNum)%functionType

		else if(i==2) then !read in cascade reaction parameters
			do while(flag .eqv. .FALSE.)
				read(DEFFILE,*) char
				if(char=='Cascade') then
					flag=.TRUE.
				end if
			end do
			flag=.FALSE.

			ImplantReactions(i,matNum)%numReactants=-10
			ImplantReactions(i,matNum)%numProducts=0
			read(DEFFILE,*) char, ImplantReactions(i,matNum)%functionType
		else
			write(*,*) 'error numImplantReac'
		end if
	end do

	close(DEFFILE)

end subroutine

!***************************************************************************************************
!> Subroutine readCascadeList(filename): reads defects in a cascade file
!This subroutine reads a list of cascades from a cascade file and stored in derived type cascade.
!***************************************************************************************************
subroutine readCascadeList(filename)
	use mod_constants
	use DerivedType
	implicit none

	character(len=50), intent(in) :: filename
	character(len=20) :: char
	type(cascadeEvent), pointer :: cascadeCurrent
	type(cascadeDefect), pointer :: defectCurrent
	integer i, numDefects, j, k
	logical flag

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
