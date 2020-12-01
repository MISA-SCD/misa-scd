!***************************************************************************************************
!> Subroutine ReadInputs():  - reads in simulation parameters from configure.in
!This subroutine reads in all simulation parameters located in configure.in as well file names
!for all other input files (defect attributes, mesh, cascades, implantation, etc).
!***************************************************************************************************
subroutine ReadInputs()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	character(len=20) :: char
	integer :: i
	character(len=100) :: defectFilename                 !<Filename of defect attributes file
	!character(len=100) :: meshFilename                   !<Filename of mesh file
	character(len=100) :: pkaFilename                    !<Filename of PKA spectrum file
	character(len=100) :: cascadeFilename                !<Filename of cascade file
	logical :: alive1, alive2, flag, flag1

	flag= .false.
	!<read in filename of defectFile
	inquire(file='configure.in', exist=alive1)
	if(.not. alive1) then
		write(*,*) 'configure.in does not exist'
	else
		open(unit=PARAFILE, file='configure.in', status='old',action='read')
	end if

	!*******************************************************
	!<read in toogles
	!*******************************************************
	!<set default valuse for toogles
	implantScheme = 'MonteCarlo'
	implantDist = 'uniform'
	grainBoundaryToggle = 'no'
	pointDefectToggle = 'no'

	!<read in filename of defect attributes file
	do while(flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='defectFile') then
			read(PARAFILE,*) defectFilename
			flag=.true.
		end if
	end do
	flag= .false.

	!<read in mesh file
!	do  while (flag .eqv. .false.)
!		read(PARAFILE,*) char
!		if(char=='meshFile') then
!			read(PARAFILE,*) meshFilename
!			flag=.true.
!		end if
!	end do
!	flag=.false.

	!*************************************************
	!<read in mesh file
!	inquire(file=meshFilename, exist=alive2)
!	if(.not. alive2) then
!		write(*,*) 'mesh file does not exist'
!	else
!		open(MESHFILE, file=meshFilename, status='old', action='read')
!	end if

!	do while(flag .eqv. .false.)
!		read(MESHFILE,*) char
!		if(char=='length') then
!			read(MESHFILE,*) meshLength
!			flag=.true.
!		end if
!	end do
!	flag=.false.

!	do while(flag .eqv. .false.)
!		read(MESHFILE,*) char
!		if(char=='numx') then
!			read(MESHFILE,*) numx
!		else if(char=='numy') then
!			read(MESHFILE,*) numy
!		else if(char=='numz') then
!			read(MESHFILE,*) numz
!			flag=.true.
!		end if
!	end do
!	flag=.false.

!	close(MESHFILE)
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

	!<Whether to use the PKA spectrum
	do while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='PKAspectrum') then
			read(PARAFILE,*) PKAspectrum
			flag=.true.
		end if
	end do
	flag=.false.

	!<read in filename of PKA spectrum file
	do while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='pkaFile') then
			read(PARAFILE,*) pkaFilename
			flag=.true.
		end if
	end do
	flag=.false.

	!<read in the number of cascade files
	do while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='numCascadeFiles') then
			read(PARAFILE,*) numCascadeFiles
			flag=.true.
		end if
	end do
	flag=.false.

	!read in filename of cascade file
	do while (flag .eqv. .false.)
		read(PARAFILE,*) char
		if(char=='cascadeFile') then
			read(PARAFILE,*) cascadeFilename
			flag=.true.
		end if
	end do
	flag=.false.


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

	!do while(flag .eqv. .false.)
	!	read(PARAFILE,*) char
	!	if(char=='implantDist') then
	!		read(PARAFILE,*) implantDist
	!		flag=.TRUE.
	!	end if
	!end do
	!flag=.false.

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
	test3 = 'no'
	tempStore = 273d0
	CuContent = 0.5d-2
	numVac = 0
	numInt = 0
	dpaRate = 1d-4
	totalDPA = 1d-1
	firr = 1d0
	!atomSize		=0d0
	lattice = 0.2876d0	!Fe
	burgers = 0.287d0
	reactionRadius = 0.65d0
	agingTime = 0d0

	polycrystal = 'no'
	grainSize = 330000
	dislocationDensity = 0d0
	impurityDensity = 0d0
	max3DInt = 4
	cascadeVolume = 0d0
	numSims = 1
	numGrains = 1
	cascadeReactionLimit = 100d0

	annealTemp = 273d0
	annealTime = 0d0
	annealSteps = 0
	annealType = 'add'
	annealTempInc = 0d0

	!Read variables in from file
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(PARAFILE,*) char
		if(char=='SimulationStart') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		flag1=.FALSE.
		do while(flag1 .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='SimulationEnd') then
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
			!else if(char=='atomSize') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) atomSize
			else if(char=='lattice') then
				flag1=.TRUE.
				read(PARAFILE,*) lattice
			else if(char=='burgers') then
				flag1=.TRUE.
				read(PARAFILE,*) burgers
			else if(char=='reactionRadius') then
				flag1=.TRUE.
				read(PARAFILE,*) reactionRadius
			!else if(char=='annealTemp') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) annealTemp
			!else if(char=='annealSteps') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) annealSteps
			else if(char=='agingTime') then		!2019.04.30 Add
				flag1=.TRUE.
				read(PARAFILE,*) agingTime
			!else if(char=='annealTime') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) annealTime
			!else if(char=='annealType') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) annealType
			!else if(char=='annealTempInc') then
			!	flag1=.TRUE.
			!	read(PARAFILE,*) annealTempInc
			else if(char=='grainSize') then
				flag1=.TRUE.
				read(PARAFILE,*) grainSize
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
	!<read in anneal parameters
	!*******************************************************
	do while(flag .eqv. .FALSE.)
		read(PARAFILE,*) char
		if(char=='AnnealStart') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		flag1=.FALSE.
		do while(flag1 .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='AnnealEnd') then
				flag1=.TRUE.
				flag=.TRUE.
			else if(char=='annealTemp') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTemp
			else if(char=='annealSteps') then
				flag1=.TRUE.
				read(PARAFILE,*) annealSteps
			else if(char=='annealTime') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTime
			else if(char=='annealType') then
				flag1=.TRUE.
				read(PARAFILE,*) annealType
			else if(char=='annealTempInc') then
				flag1=.TRUE.
				read(PARAFILE,*) annealTempInc
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
	totdatToggle = 'yes'
	defectToggle = 'no'
	stadatToggle = 'no'
	xyzdatToggle = 'no'
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
			if(char=='OutputEnd') then
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
			else if(char=='xyzdatToggle') then
				flag1=.TRUE.
				read(PARAFILE,*) xyzdatToggle
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
	!if(meshingType=='adaptive') then
		!flag=.FALSE.

		do while(flag .eqv. .FALSE.)
			read(PARAFILE,*) char
			if(char=='MeshStart') then
				flag=.TRUE.
			endif
		end do
		flag=.FALSE.

		do while(flag .eqv. .FALSE.)
			flag1=.FALSE.
			do while(flag1 .eqv. .FALSE.)
				read(PARAFILE,*) char
				if(char=='MeshEnd') then
					flag1=.TRUE.
					flag=.TRUE.
				else if(char=='length') then
					flag1=.TRUE.
					read(PARAFILE,*) meshLength
				else if(char=='numx') then
					flag1=.TRUE.
					read(PARAFILE,*) numx
				else if(char=='numy') then
					flag1=.TRUE.
					read(PARAFILE,*) numy
				else if(char=='numz') then
					flag1=.TRUE.
					read(PARAFILE,*) numz
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
	!end if

	close(PARAFILE)

	!*******************************************************
	!<Read other files
	!*******************************************************
	!<read _defects_.txt
	call readDefectAttributes(defectFilename)
	!<read cpdf.*
	if(PKAspectrum == 'yes') then
		call readPKAspectrum(pkaFilename)
	end if

	!<read cascade.txt
	if(implantType=='Cascade') then
		if(numCascadeFiles > 1) then
			call readCascadeFiles(cascadeFilename)
		else
			call readCascadeList(cascadeFilename)
		end if
	else if(implantType=='FrenkelPair') then
		numDisplacedAtoms=1				!Frenkel pair implantation, no need to read cascades, one displaced atom per implant event
	else
		write(*,*) 'error implantType'
	end if

	!***********************************************************************
	!clustering rate constants
	!***********************************************************************
	atomSize=(lattice**3d0)/2d0
	omega=(48d0*pi**2/atomSize**2)**(1d0/3d0) 			!clustering rate parameter for spherical clusters
	omegastar=(4*pi*reactionRadius)/atomSize			!clustering rate parameter modifier due to reaction radius
	omega2D=(4d0*pi/(atomSize*burgers))**(1d0/2d0)		!clustering rate parameter for 1D migrating circular clusters
	omega1D=(9d0*pi/(16d0*atomSize))**(1d0/6d0)			!clustering rate parameter for 1D migrating spherical clusters
	omegastar1D=reactionRadius*(pi/atomSize)**(1d0/2d0)
	!omegastar1D=0d0									!clustering rate parameter modifier due to reaction radius
	omegacircle1D=(1d0/burgers)**(1d0/2d0)				!clustering rate parameter for 1D migrating circular clusters


end subroutine

!***************************************************************************************************
!>Subroutine readDefectAttributes(filename). Information read in includes:
!   1) Number of allowed defect types
!   2) Formation energies, diffusion prefactors, migration energies and binding energies
!   3) Allowed reactions
!***************************************************************************************************
subroutine readDefectAttributes(filename)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	character(len=100), intent(in) :: filename
	character(len=20) :: char
	logical :: flag
	integer :: i,j

	flag=.FALSE.
	open(ATTRFILE, file=filename,action='read', status='old')

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
		allocate(FormSingle(i)%defectType(SPECIES))
		read(ATTRFILE,*) (FormSingle(i)%defectType(j),j=1,SPECIES)
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
		allocate(DiffSingle(i)%defectType(SPECIES))
		read(ATTRFILE,*) (DiffSingle(i)%defectType(j),j=1,SPECIES)
		read(ATTRFILE,*) char, DiffSingle(i)%D0, char, DiffSingle(i)%Em
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
		allocate(DiffFunc(i)%defectType(SPECIES))
		read(ATTRFILE,*) (DiffFunc(i)%defectType(j),j=1,SPECIES)	!< read in defectTypes
		allocate(DiffFunc(i)%min(SPECIES))
		allocate(DiffFunc(i)%max(SPECIES))
		read(ATTRFILE,*) char, (DiffFunc(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (DiffFunc(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, DiffFunc(i)%fType
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
		allocate(BindSingle(i)%defectType(SPECIES))
		allocate(BindSingle(i)%product(SPECIES))
		read(ATTRFILE,*) (BindSingle(i)%defectType(j),j=1,SPECIES),(BindSingle(i)%product(j),j=1,SPECIES)
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
		allocate(BindFunc(i)%defectType(SPECIES))
		allocate(BindFunc(i)%product(SPECIES))
		read(ATTRFILE,*) (BindFunc(i)%defectType(j),j=1,SPECIES),(BindFunc(i)%product(j),j=1,SPECIES)
		allocate(BindFunc(i)%min(SPECIES))
		allocate(BindFunc(i)%max(SPECIES))
		read(ATTRFILE,*) char, (BindFunc(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (BindFunc(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, BindFunc(i)%fType
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
		allocate(DissocReactions(i)%reactants(SPECIES,DissocReactions(i)%numReactants))
		allocate(DissocReactions(i)%products(SPECIES,DissocReactions(i)%numProducts))
		read(ATTRFILE,*) (DissocReactions(i)%reactants(j,1),j=1,SPECIES),&
				(DissocReactions(i)%products(j,1),j=1,SPECIES)	!< read in defectType
		allocate(DissocReactions(i)%min(SPECIES))
		allocate(DissocReactions(i)%max(SPECIES))
		read(ATTRFILE,*) char, (DissocReactions(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (DissocReactions(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, DissocReactions(i)%fType
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
		allocate(DiffReactions(i)%reactants(SPECIES,DiffReactions(i)%numReactants))
		allocate(DiffReactions(i)%products(SPECIES,DiffReactions(i)%numProducts))
		read(ATTRFILE,*) (DiffReactions(i)%reactants(j,1),j=1,SPECIES),&
				(DiffReactions(i)%products(j,1),j=1,SPECIES)
		allocate(DiffReactions(i)%min(SPECIES))
		allocate(DiffReactions(i)%max(SPECIES))
		read(ATTRFILE,*) char, (DiffReactions(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (DiffReactions(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, DiffReactions(i)%fType
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
		allocate(SinkReactions(i)%reactants(SPECIES,SinkReactions(i)%numReactants))
		read(ATTRFILE,*) (SinkReactions(i)%reactants(j,1),j=1,SPECIES)
		allocate(SinkReactions(i)%min(SPECIES))
		allocate(SinkReactions(i)%max(SPECIES))
		read(ATTRFILE,*) char, (SinkReactions(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (SinkReactions(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, SinkReactions(i)%fType
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
		allocate(ImpurityReactions(i)%reactants(SPECIES,ImpurityReactions(i)%numReactants))
		allocate(ImpurityReactions(i)%products(SPECIES,ImpurityReactions(i)%numProducts))
		read(ATTRFILE,*) (ImpurityReactions(i)%reactants(j,1),j=1,SPECIES), &
				(ImpurityReactions(i)%products(j,1),j=1,SPECIES)
		allocate(ImpurityReactions(i)%min(SPECIES))
		allocate(ImpurityReactions(i)%max(SPECIES))
		read(ATTRFILE,*) char, (ImpurityReactions(i)%min(j),j=1,SPECIES)
		read(ATTRFILE,*) char, (ImpurityReactions(i)%max(j),j=1,SPECIES)
		read(ATTRFILE,*) char, ImpurityReactions(i)%fType
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
		allocate(ClusterReactions(i)%reactants(SPECIES,ClusterReactions(i)%numReactants))
		allocate(ClusterReactions(i)%products(SPECIES,ClusterReactions(i)%numProducts))
		read(ATTRFILE,*) (ClusterReactions(i)%reactants(j,1),j=1,SPECIES),&
				(ClusterReactions(i)%reactants(j,2),j=1,SPECIES)
		allocate(ClusterReactions(i)%min(SPECIES*ClusterReactions(i)%numReactants))
		allocate(ClusterReactions(i)%max(SPECIES*ClusterReactions(i)%numReactants))
		read(ATTRFILE,*) char,(ClusterReactions(i)%min(j),j=1,SPECIES*ClusterReactions(i)%numReactants)
		read(ATTRFILE,*) char,(ClusterReactions(i)%max(j),j=1,SPECIES*ClusterReactions(i)%numReactants)
		do j=1,SPECIES
			ClusterReactions(i)%products(j,1)=ClusterReactions(i)%reactants(j,1)+&
					ClusterReactions(i)%reactants(j,2)
		end do
		read(ATTRFILE,*) char, ClusterReactions(i)%fType
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
			allocate(ImplantReactions(i)%products(SPECIES,ImplantReactions(i)%numProducts))
			read(ATTRFILE,*) (ImplantReactions(i)%products(j,1),j=1,SPECIES),&
					(ImplantReactions(i)%products(j,2),j=1,SPECIES)
			read(ATTRFILE,*) char, ImplantReactions(i)%fType

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
			read(ATTRFILE,*) char, ImplantReactions(i)%fType
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
	use mod_globalVariables
	use mod_structures
	implicit none

	character(len=100), intent(in) :: filename
	character(len=20) :: char
	type(cascadeEvent), pointer :: cascadeCurrent
	type(cascadeDefect), pointer :: defectCurrent
	integer :: i, j, k
	logical :: flag

	open(CASFILE, file=filename, status='old', action='read')

	allocate(cascadeList)
	cascadeCurrent=>cascadeList
	read(CASFILE,*) PKAtemperature
	read(CASFILE,*) PKAenergy
	read(CASFILE,*) numDisplacedAtoms
	read(CASFILE,*) numCascades

	do i=1,numCascades
		read(CASFILE,*)
		read(CASFILE,*) cascadeCurrent%numDefectsTotal
		read(CASFILE,*) cascadeCurrent%numDisplacedAtoms
		allocate(cascadeCurrent%ListOfDefects)
		defectCurrent=>cascadeCurrent%ListOfDefects
		nullify(cascadeCurrent%next)
		nullify(defectCurrent%next)

		do j=1, cascadeCurrent%numDefectsTotal
			allocate(defectCurrent%defectType(SPECIES))
			read(CASFILE,*) (defectCurrent%defectType(k),k=1,SPECIES)
			read(CASFILE,*) (defectCurrent%coordinates(k), k=1,3)

			if(j /= cascadeCurrent%numDefectsTotal) then
				allocate(defectCurrent%next)
				nullify(defectCurrent%next%next)
				defectCurrent=>defectCurrent%next
			end if
		end do

		if(i /= numCascades) then
			allocate(cascadeCurrent%next)
			cascadeCurrent=>cascadeCurrent%next
		end if
	end do
	nullify(defectCurrent)
	nullify(cascadeCurrent)

	close(CASFILE)

end subroutine

!***************************************************************************************************
!> Subroutine readEPKAs(filename): reads PKA spectrum file
!The first column is the PKA energy, and the second column is the cpdf
!***************************************************************************************************
subroutine readPKAspectrum(filename)
	use mod_constants
	use mod_structures
	use mod_globalvariables
	implicit none

	character(len=100), intent(in) :: filename       !< filename = '../../inputs/pkas/cpdf.w
	integer :: i

	open(CPDFFILE, file=filename, status='old', action='read')
	read(CPDFFILE, *) EPKAlist%size
	read(CPDFFILE, *)
	allocate(EPKAlist%energy(EPKAlist%size))
	allocate(EPKAlist%cpdf(EPKAlist%size))

	do i=1, EPKAlist%size
		read(CPDFFILE, *) EPKAlist%energy(i), EPKAlist%cpdf(i)
	end do
	close(CPDFFILE)

end subroutine readPKAspectrum

!***************************************************************************************************
!> Subroutine readCascadeFiles(filename): reads cascade defects in a number of cascade files
!This subroutine reads a list of cascades from each cascade file and stored in derived type cascade.
!***************************************************************************************************
subroutine readCascadeFiles(filename)
	use mod_constants
	use mod_structures
	use mod_globalvariables
	implicit none

	character(len=100), intent(in) :: filename       !< filename = '../../inputs/cascades/Fe_*.txt'
	character(len=100) :: systemChars
	character(len=100) :: cascadeFilename
	integer :: fileID, stat, totalDisAtoms, totalCascades, i, j, k
	type(cascadeEvent), pointer :: casCurrent
	type(cascadeDefect), pointer :: casDef

	allocate(cascadeLists(numCascadeFiles))
	totalDisAtoms=0
	totalCascades=0

	systemChars = 'ls -R '//trim(filename)//' > cas.dat'	!write all cascade files' relative path to cas.dat file
	call system(trim(systemChars))

	open(20, file='cas.dat')
	do fileID=1, numCascadeFiles
		!<read filename
		read(20, fmt='(A)', iostat=stat) cascadeFilename
		if(stat /= 0) then
			write(*,*) 'Error: read cascade filename'
		end if
		open(CASFILE, file=cascadeFilename, status='old', action='read')
		!<read number of atoms displaced and number of cascades in this cascade file
		read(CASFILE,*) cascadeLists(fileID)%temperature   !(K)
		read(CASFILE,*) cascadeLists(fileID)%PKAenergy     !(eV)
		read(CASFILE,*) cascadeLists(fileID)%averDisAtoms
		read(CASFILE,*) cascadeLists(fileID)%numCascades
		totalDisAtoms=totalDisAtoms+cascadeLists(fileID)%averDisAtoms*cascadeLists(fileID)%numCascades
		totalCascades=totalCascades+cascadeLists(fileID)%numCascades
		allocate(cascadeLists(fileID)%listCascades)

		!read cascades in this file
		nullify(casCurrent)
		casCurrent=>cascadeLists(fileID)%listCascades
		nullify(casCurrent%next)
		do i=1, cascadeLists(fileID)%numCascades
			read(CASFILE,*)
			read(CASFILE,*) casCurrent%numDefectsTotal
			read(CASFILE,*) casCurrent%numDisplacedAtoms
			allocate(casCurrent%ListOfDefects)
			nullify(casDef)
			casDef=>casCurrent%ListOfDefects
			nullify(casDef%next)

			do j=1,casCurrent%numDefectsTotal
				allocate(casDef%defectType(SPECIES))
				read(CASFILE,*) (casDef%defectType(k),k=1,SPECIES)
				read(CASFILE,*) (casDef%coordinates(k), k=1,3)

				if(j /= casCurrent%numDefectsTotal) then
					allocate(casDef%next)
					casDef=>casDef%next
					nullify(casDef%next)
				end if
			end do

			if(i /= cascadeLists(fileID)%numCascades) then
				allocate(casCurrent%next)
				casCurrent=>casCurrent%next
				nullify(casCurrent%next)
			end if
		end do
		close(CASFILE)
	end do
	numDisplacedAtoms=dble(totalDisAtoms)/dble(totalCascades)

	close(20, status='delete')      !<close  and delete cas.dat

end subroutine readCascadeFiles