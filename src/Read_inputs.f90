
!***************************************************************************************************
!>Subroutine read material input - reads material input information from file
!!
!!Information read in includes:
!!
!!1) Allowed defects and their diffusion rates and binding energies
!!
!!2) Allowed reactions, including single, multi, diffusion, and implantation reactions
!!
!!NOTE: we are double-reading in the values of numSingleDiff(matNum), etc, because it has
!!already been done in the subroutine readReactionListSizes().
!***************************************************************************************************

subroutine readMaterialInput(filename)	!read FeCu_Defects.txt
use DerivedType
use mod_constants
implicit none

integer i, j, count
character*20 :: char
character*50 :: filename
integer, allocatable :: DefectType(:), productType(:)
!type(reaction), pointer :: reactions, reactionCurrent
logical flag
Double precision Diff, Eb
integer matNum

interface
	double precision function findDiffusivity(DefectType)
	integer, allocatable :: DefectType(:)
	end function
	
	double precision function findBinding(DefectType, productType)
	integer, allocatable :: DefectType(:), productType(:)
	end function
end interface

open(80, file=filename,action='read', status='old') !read xx_Defects.txt

flag=.FALSE.

!The following is for the entire parameter set

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='material') then
		flag=.TRUE.
		read(80,*) matNum
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='species') then	
		flag=.TRUE.
		read(80,*) numSpecies
	end if
end do
flag=.FALSE.

!************************************************
!The following is for formation energies parameters only
!************************************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='formationEnergies') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleForm(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='single') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do i=1,numSingleForm(matNum)
	allocate(FormSingle(matNum,i)%defectType(numSpecies))
	read(80,*) (FormSingle(matNum,i)%defectType(j),j=1,numSpecies)
	read(80,*) char, FormSingle(matNum,i)%Ef
end do

!************************************************
!The following is for diffusivity parameters only
!************************************************

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='diffusionPrefactors') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleDiff(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numFunction') then
		flag=.TRUE.
		read(80,*) numFuncDiff(matNum)
	end if
end do
flag=.FALSE.

!allocate(DiffSingle(numSingleDiff))
!allocate(DiffFunc(numFuncDiff))

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='single') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!write(*,*) 'reading single defect diffusion, proc', myProc%taskid
do i=1,numSingleDiff(matNum)
	allocate(DiffSingle(matNum,i)%defectType(numSpecies))
	read(80,*) (DiffSingle(matNum,i)%defectType(j),j=1,numSpecies)
	read(80,*) char, DiffSingle(matNum,i)%D, char, DiffSingle(matNum,i)%Em
end do

!do 16 i=1,numSingleDiff(matNum)
!	write(*,*) 'species', (DiffSingle(matNum,i)%defectType(j),j=1,numSpecies)
!	write(*,*) 'D0', DiffSingle(matNum,i)%D, 'Em', DiffSingle(matNum,i)%Em
!16 continue

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='function') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!write(*,*) 'reading function diffusion values', myProc%taskid
do i=1,numFuncDiff(matNum)
	allocate(DiffFunc(matNum,i)%defectType(numSpecies))
	read(80,*) (DiffFunc(matNum,i)%defectType(j),j=1,numSpecies)	!< read in defectTypes
	count=0
	!We need to know min and max sizes for each defect type included
	allocate(DiffFunc(matNum,i)%min(numSpecies))
	allocate(DiffFunc(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, DiffFunc(matNum,i)%min(j), char, DiffFunc(matNum,i)%max(j)	!< read in min and max
	end do
	read(80,*) char, DiffFunc(matNum,i)%functionType, char, DiffFunc(matNum,i)%numParam	!< read in functionType and numParm
	allocate(DiffFunc(matNum,i)%parameters(DiffFunc(matNum,i)%numParam))
	read(80,*) (DiffFunc(matNum,i)%parameters(j),j=1,DiffFunc(matNum,i)%numParam)	!< read in paramsters
end do

!< test readIn
!do 22 i=1,numFuncDiff
!	write(*,*) 'defect type'
!	write(*,*)	(DiffFunc(i)%defectType(j),j=1,numSpecies)
!	do 24 j=1,numSpecies
!		write(*,*) 'min species', j,'=',DiffFunc(i)%min(j), 'max=', DiffFunction(i)%max(j)
!	24 continue
!	write(*,*) 'Function type', DiffFunc(i)%functionType
!	write(*,*) 'Parameters', (DiffFunc(i)%parameters(j), j=1,DiffFunc(i)%numParam)
!22 continue

!allocate(DefectType(numSpecies))
!write(*,*) 'Test: enter defect type (numSpecies values) to find diffusivity'
!read(*,*) (DefectType(i), i=1,numSpecies)
!
!Diff=findDiffusivity(DefectType)
!write(*,*) 'Diffusivity', Diff

!*****************************************************
!The following is for binding energies parameters only
!*****************************************************

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='bindingEnergies') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!write(*,*) 'reading binding energies', myProc%taskid

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleBind(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numFunction') then
		flag=.TRUE.
		read(80,*) numFuncBind(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='single') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!allocate(BindSingle(numSingleBind))
!allocate(BindFunc(numFuncBind))

do i=1,numSingleBind(matNum)
	allocate(BindSingle(matNum,i)%defectType(numSpecies))
	allocate(BindSingle(matNum,i)%product(numSpecies))
	read(80,*) (BindSingle(matNum,i)%defectType(j),j=1,numSpecies), (BindSingle(matNum,i)%product(j),j=1,numSpecies)
	read(80,*) char, BindSingle(matNum,i)%Eb
end do

!do 30 i=1,numSingleBind
!	write(*,*) 'species', (BindSingle(i)%defectType(j),j=1,numSpecies), 'product', (BindSingle(i)%product(j),j=1,numSpecies)
!	write(*,*) 'Eb', BindSingle(i)%Eb
!30 continue

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='function') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!write(*,*) 'reading function binding values', myProc%taskid
do i=1,numFuncBind(matNum)
	allocate(BindFunc(matNum,i)%defectType(numSpecies))
	allocate(BindFunc(matNum,i)%product(numSpecies))
	read(80,*) (BindFunc(matNum,i)%defectType(j),j=1,numSpecies), (BindFunc(matNum,i)%product(j),j=1,numSpecies)
	allocate(BindFunc(matNum,i)%min(numSpecies))
	allocate(BindFunc(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, BindFunc(matNum,i)%min(j), char, BindFunc(matNum,i)%max(j)
	end do
	read(80,*) char, BindFunc(matNum,i)%functionType, char, BindFunc(matNum,i)%numParam
	allocate(BindFunc(matNum,i)%parameters(BindFunc(matNum,i)%numParam))
	read(80,*) (BindFunc(matNum,i)%parameters(j),j=1,BindFunc(matNum,i)%numParam)
end do

!< test readIn
!allocate(DefectType(numSpecies))
!allocate(productType(numSpecies))
!write(*,*) 'Test: enter defect type (numSpecies values) to find binding energy'
!read(*,*) (DefectType(i), i=1,numSpecies)
!write(*,*) 'Enter product type'
!read(*,*) (productType(i),i=1,numSpecies)

!Eb=findBinding(DefectType, productType)
!write(*,*) 'Binding Energy', Eb

!******************************************************************
!Construct reaction list using migration and binding energies above
!******************************************************************

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='singleDefect') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!********************** dissociation *******************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='dissociation') then
		flag=.TRUE.
		read(80,*) numDissocReac(matNum)
	end if
end do
flag=.FALSE.

!write(*,*) 'reading dissociation reactions', myProc%taskid

!allocate(DissocReactions(numDissocReac))
do i=1,numDissocReac(matNum)
	DissocReactions(matNum,i)%numReactants=1
	DissocReactions(matNum,i)%numProducts=1
	allocate(DissocReactions(matNum,i)%reactants(DissocReactions(matNum,i)%numReactants,numSpecies))
	allocate(DissocReactions(matNum,i)%products(DissocReactions(matNum,i)%numProducts,numSpecies))
	read(80,*) (DissocReactions(matNum,i)%reactants(1,j),j=1,numSpecies),&
		(DissocReactions(matNum,i)%products(1,j),j=1,numSpecies)	!< read in defectType
	allocate(DissocReactions(matNum,i)%min(numSpecies))
	allocate(DissocReactions(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, DissocReactions(matNum,i)%min(j), char, DissocReactions(matNum,i)%max(j)
	end do
	read(80,*) DissocReactions(matNum,i)%functionType
end do


!********************** diffusion *******************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='diffusion') then
		flag=.TRUE.
		read(80,*) numDiffReac(matNum)
	end if
end do
flag=.FALSE.

!write(*,*) 'reading diffusion reactions', myProc%taskid

!allocate(DiffReactions(numDiffReac))
do i=1,numDiffReac(matNum)
	DiffReactions(matNum,i)%numReactants=1
	DiffReactions(matNum,i)%numProducts=1
	allocate(DiffReactions(matNum,i)%reactants(DiffReactions(matNum,i)%numReactants,numSpecies))
	allocate(DiffReactions(matNum,i)%products(DiffReactions(matNum,i)%numProducts,numSpecies))
	read(80,*) (DiffReactions(matNum,i)%reactants(1,j),j=1,numSpecies),&
		(DiffReactions(matNum,i)%products(1,j),j=1,numSpecies)
	allocate(DiffReactions(matNum,i)%min(numSpecies))
	allocate(DiffReactions(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, DiffReactions(matNum,i)%min(j), char, DiffReactions(matNum,i)%max(j)
	end do
	read(80,*) DiffReactions(matNum,i)%functionType
end do

!********************** sinkRemoval *******************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='sinkRemoval') then
		flag=.TRUE.
		read(80,*) numSinkReac(matNum)
	end if
end do
flag=.FALSE.

!write(*,*) 'reading sink removal reactions', myProc%taskid

!allocate(SinkReactions(numSinkReac))
do i=1,numSinkReac(matNum)
	SinkReactions(matNum,i)%numReactants=1
	SinkReactions(matNum,i)%numProducts=0
	allocate(SinkReactions(matNum,i)%reactants(SinkReactions(matNum,i)%numReactants,numSpecies))
	read(80,*) (SinkReactions(matNum,i)%reactants(1,j),j=1,numSpecies)
	allocate(SinkReactions(matNum,i)%min(numSpecies))
	allocate(SinkReactions(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, SinkReactions(matNum,i)%min(j), char, SinkReactions(matNum,i)%max(j)
	end do
	read(80,*) SinkReactions(matNum,i)%functionType
end do


!********************** impurityTrapping *******************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='impurityTrapping') then
		flag=.TRUE.
		read(80,*) numImpurityReac(matNum)
	end if
end do
flag=.FALSE.

!write(*,*) 'reading impurity trapping reactions', myProc%taskid

!allocate(ImpurityReactions(numImpurityReac))
do i=1,numImpurityReac(matNum)
	ImpurityReactions(matNum,i)%numReactants=1
	ImpurityReactions(matNum,i)%numProducts=1
	allocate(ImpurityReactions(matNum,i)%reactants(ImpurityReactions(matNum,i)%numReactants,numSpecies))
	allocate(ImpurityReactions(matNum,i)%products(ImpurityReactions(matNum,i)%numProducts,numSpecies))
	read(80,*) (ImpurityReactions(matNum,i)%reactants(1,j),j=1,numSpecies), &
		(ImpurityReactions(matNum,i)%products(1,j),j=1,numSpecies)
	allocate(ImpurityReactions(matNum,i)%min(numSpecies))
	allocate(ImpurityReactions(matNum,i)%max(numSpecies))
	do j=1,numSpecies
		read(80,*) char, ImpurityReactions(matNum,i)%min(j), char, ImpurityReactions(matNum,i)%max(j)
	end do
	read(80,*) ImpurityReactions(matNum,i)%functionType
end do

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='multipleDefect') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.


!********************** clustering *******************************
do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='clustering') then
		flag=.TRUE.
		read(80,*) numClusterReac(matNum)
	end if
end do
flag=.FALSE.

!write(*,*) 'reading clustering reactions', myProc%taskid

!allocate(ClusterReactions(numClusterReac))
do i=1,numClusterReac(matNum)
	ClusterReactions(matNum,i)%numReactants=2
	ClusterReactions(matNum,i)%numProducts=1
	allocate(ClusterReactions(matNum,i)%reactants(ClusterReactions(matNum,i)%numReactants,numSpecies))
	allocate(ClusterReactions(matNum,i)%products(ClusterReactions(matNum,i)%numProducts,numSpecies))
	read(80,*) (ClusterReactions(matNum,i)%reactants(1,j),j=1,numSpecies),&
		(ClusterReactions(matNum,i)%reactants(2,j),j=1,numSpecies)
	allocate(ClusterReactions(matNum,i)%min(numSpecies*ClusterReactions(matNum,i)%numReactants))
	allocate(ClusterReactions(matNum,i)%max(numSpecies*ClusterReactions(matNum,i)%numReactants))
	do j=1,numSpecies*ClusterReactions(matNum,i)%numReactants
		read(80,*) char, ClusterReactions(matNum,i)%min(j), char, ClusterReactions(matNum,i)%max(j)
	end do
	do j=1,numSpecies
		!ClusterReactions products are reaction-specific, and are not correctly found here. This is just a placeholder.
		ClusterReactions(matNum,i)%products(1,j)=ClusterReactions(matNum,i)%reactants(1,j)+&
												 ClusterReactions(matNum,i)%reactants(2,j)
	end do
	read(80,*) ClusterReactions(matNum,i)%functionType
end do


do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='noDefect') then
		flag=.TRUE.
		read(80,*) numImplantReac(matNum)
	end if
end do
flag=.FALSE.

!allocate(ImplantReactions(numImplantReac))
!********************** FrenkelPair Implantation *******************************
do i=1,numImplantReac(matNum)
	if(i==1) then	!Read in Frenkel pair reaction parameters

		do while(flag .eqv. .FALSE.)
			read(80,*) char
			if(char=='FrenkelPair') then
				flag=.TRUE.
			end if
		end do
		flag=.FALSE.
		
		!write(*,*) 'reading frenkel pair reactions', myProc%taskid
		
		ImplantReactions(matNum,i)%numReactants=0
		ImplantReactions(matNum,i)%numProducts=2
		allocate(ImplantReactions(matNum,i)%products(ImplantReactions(matNum,i)%numProducts,numSpecies))
		read(80,*) (ImplantReactions(matNum,i)%products(1,j),j=1,numSpecies),&
			(ImplantReactions(matNum,i)%products(2,j),j=1,numSpecies)
		read(80,*) ImplantReactions(matNum,i)%functionType

	!else if(i==2) then	!read in He implantation parameters
	
	!	do 61 while(flag .eqv. .FALSE.)
	!		read(80,*) char
	!		if(char=='HeImplant') then
	!			flag=.TRUE.
	!		endif
	!	61 continue
	!	flag=.FALSE.
		
		!Reading He implantation reactions
		
	!	ImplantReactions(matNum,i)%numReactants=0
	!	ImplantReactions(matNum,i)%numProducts=1
	!	allocate(ImplantReactions(matNum,i)%products(ImplantReactions(matNum,i)%numProducts,numSpecies))
	!	read(80,*) (ImplantReactions(matNum,i)%products(1,j),j=1,numSpecies)
	!	read(80,*) ImplantReactions(matNum,i)%functionType
	
	else if(i==2) then !read in cascade reaction parameters

		do while(flag .eqv. .FALSE.)
			read(80,*) char
			if(char=='Cascade') then
				flag=.TRUE.
			end if
		end do
		flag=.FALSE.
		
		!write(*,*) 'reading cascade reactions', myProc%taskid
		
		ImplantReactions(matNum,i)%numReactants=-10
		ImplantReactions(matNum,i)%numProducts=0
		read(80,*) ImplantReactions(matNum,i)%functionType

	else
		write(*,*) 'error numImplantReac'
	end if
end do

close(80)
end subroutine

!***************************************************************************************************
!
!> Subroutine readCascadeList() - reads defects in cascades and locations from a file
!!
!! This subroutine reads a list of cascades from a file (name given in readParameters).
!! This information is stored in derived type cascadeEvent and used whenever a cascade is chosen
!! throughout the duration of the program.
!
!***************************************************************************************************

subroutine readCascadeList()
use mod_constants
use DerivedType

implicit none

type(cascadeEvent), pointer :: cascadeCurrent
type(cascadeDefect), pointer :: defectCurrent
integer i, numDefects, j, k
character*20 char
character*50 filename
logical flag

!read in implantation type
flag=.FALSE.
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='implantType') then
		read(81,*) implantType
		flag=.TRUE.
	end if
end do

!read in cascade implantation scheme (Monte Carlo or explicit)
flag=.FALSE.
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='implantScheme') then
		read(81,*) implantScheme
		flag=.TRUE.
	end if
end do

!Check to make sure that we only are choosing explicit implantation with cascades
if(implantScheme=='explicit' .AND. implantType=='FrenkelPair') then
	write(*,*) 'error frenkel pairs with explicit implantation'
endif

!read in filename of cascade file
flag=.FALSE.
do while (flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='cascadeFile') then
		read(81,*) filename			!read in filename of cascade input file from parameters.txt
		flag=.TRUE.
	end if
end do

flag=.FALSE.
do while (flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='meshingType') then
		read(81,*) meshingType		!read whether we are using adaptive meshing or fine meshing
		flag=.TRUE.
	end if
end do

if(implantType=='Cascade') then

	open(80, file=filename,action='read', status='old')

	allocate(cascadeList)
	cascadeCurrent=>cascadeList
	read(80,*) numDisplacedAtoms
	read(80,*)
	read(80,*) numCascades
	
	do i=1,numCascades
		read(80,*)
		read(80,*) numDefects
		read(80,*) cascadeCurrent%numDisplacedAtoms
		cascadeCurrent%NumDefectsTotal=numDefects
		allocate(cascadeCurrent%ListOfDefects)
		defectCurrent=>cascadeCurrent%ListOfDefects
		nullify(cascadeCurrent%nextCascade)
		nullify(defectCurrent%next)
		
		do j=1,numDefects
			allocate(defectCurrent%defectType(numSpecies))
			read(80,*) (defectCurrent%defectType(k),k=1,numSpecies)
			read(80,*) (defectCurrent%coordinates(k), k=1,3)
			
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
	
	!output read-in cascade list to make sure it works
	!if(myProc%taskid==MASTER) then
	!	cascadeCurrent=>CascadeList
	!	do 12 while(associated(CascadeCurrent))
	!		write(*,*) 'Cascade Event'
	!		defectCurrent=>CascadeCurrent%ListOfDefects
	!		write(*,*) 'total number of defects', cascadeCurrent%NumDefectsTotal
	!		do 13 while(associated(defectCurrent))
	!			write(*,*) 'type', defectCurrent%defectType, 'coordinates', defectCurrent%coordinates
	!			defectCurrent=>defectCurrent%next
	!		13 continue
	!		CascadeCurrent=>CascadeCurrent%nextCascade
	!	12 continue
	!	write(*,*)
	!	read(*,*)
	!	nullify(defectCurrent)
	!	nullify(cascadeCurrent)
	!endif
	
	close(80)

else if(implantType=='FrenkelPair') then
	numDisplacedAtoms=1				!Frenkel pair implantation, no need to read cascades, one displaced atom per implant event
else
	write(*,*) 'error implantType'
endif

end subroutine

!***************************************************************************************************
!
!> Subroutine readImplantData() - reads non-uniform implantation profiles (damage and helium)
!!
!! This subroutine recognizes whether we are in a uniform implantation scheme or a non-uniform implantation
!! scheme. If the implantation is non-uniform, this subroutine reads from a file the implantation
!! profile (in dpa/s for each material point).
!!
!! Note: as constructed, this is hard-coded to read in one-dimensional DPA and He implantation rates
!! (for DPA profiles through the thickness of a material, for example), instead of full 3D DPA distributions.
!! This is typically for thin films with implantation profiles that vary through their depth.
!!
!! This is written such that the input file does not have to have the same number of points with the
!! same z-coordinates as the elements in the mesh. The code will interpolate from the input file what
!! the DPA rate and He implant rate should be in each element. However, if the size of the input file
!! is smaller than the size of the mesh in the z-direction, it will return an error.
!!
!! Inputs: file with DPA and He implant rates in z-coordinates
!! Outputs: information is stored in a global array.
!
!***************************************************************************************************

subroutine readImplantData()
use DerivedType
use mod_constants
implicit none

logical flag
character*20 char
character*50 filename
integer i, j

!read in toggle for non-unifrm distribution
flag=.FALSE.
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='implantDist') then
		read(81,*) implantDist
		flag=.TRUE.
	end if
end do

flag=.FALSE.
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='implantFile') then
		read(81,*) filename
		flag=.TRUE.
	end if
end do

if(implantDist=='Uniform') then
	!Do nothing, no need to read information from file
else if(implantDist=='NonUniform') then
	
	!Read in implantation profile from file and store in a double precision array.
	open(80, file=filename,action='read', status='old')
	
	flag=.FALSE.
	do while(flag .eqv. .FALSE.)
		read(80,*) char
		if(char=='numImplantDataPoints') then
			read(80,*) numImplantDataPoints
			flag=.TRUE.
		end if
	end do
	
	!Allocate the array containing information on DPA rates with numImplantDataPoints rows and
	!3 columns (z-coordinates, DPA rate, He implant rate)

	!2019.05.06 modify: no He
	!allocate(implantRateData(numImplantDataPoints,3))
	allocate(implantRateData(numImplantDataPoints,2))
	
	flag=.FALSE.
	do while(flag .eqv. .FALSE.)
		read(80,*) char
		if(char=='start') then
			flag=.TRUE.
		end if
	end do
	
	!Read implant information from file
	do i=1,numImplantDataPoints
		!2019.05.06 modify: no He
		!read(80,*) (implantRateData(i,j),j=1,3)
		read(80,*) (implantRateData(i,j),j=1,2)
	end do
	
	close(80)
	
else
	write(*,*) 'Error: unknown implantation distribution'
endif

end subroutine

!***************************************************************************************************
!
!> Subroutine selectMaterialInputs() - controlling subroutines that calls other subroutines to read inputs.
!!
!! This subroutine reads the name of the material input file(s) and then reads in relevant material
!! constants (binding and migration energies, number of species, allowed reactions, reaction functional
!! forms, etc)
!
!***************************************************************************************************

subroutine selectMaterialInputs()
use mod_constants
use DerivedType
implicit none

character*20 char
character*50, allocatable :: filename(:)
logical flag
integer i, maxNum

flag=.FALSE.

!read in number of materials
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='numMaterials') then
		read(81,*) numMaterials
		flag=.TRUE.
	end if
end do
flag=.FALSE.

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

allocate(filename(numMaterials))

!Figure out how big to make all of the lists of allowed reactions, etc

do i=1,numMaterials
	!read in filename of mesh file
	do while (flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='materialFile') then
			read(81,*) filename(i)
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	call readReactionListSizes(filename(i)) !many xx_Defects.txt files

end do

!Allocate all lists of binding, diffusion, and reactions
maxNum=0
do i=1,numMaterials
	if(numSingleForm(i) > maxNum) then
		maxNum=numSingleForm(i)
	end if
end do
allocate(FormSingle(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numSingleDiff(i) > maxNum) then
		maxNum=numSingleDiff(i)
	end if
end do
allocate(DiffSingle(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numFuncDiff(i) > maxNum) then
		maxNum=numFuncDiff(i)
	end if
end do
allocate(DiffFunc(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numSingleBind(i) > maxNum) then
		maxNum=numSingleBind(i)
	end if
end do
allocate(BindSingle(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numFuncBind(i) > maxNum) then
		maxNum=numFuncBind(i)
	end if
end do
allocate(BindFunc(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numDiffReac(i) > maxNum) then
		maxNum=numDiffReac(i)
	end if
end do
allocate(DiffReactions(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numClusterReac(i) > maxNum) then
		maxNum=numClusterReac(i)
	endif
end do
allocate(ClusterReactions(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numImplantReac(i) > maxNum) then
		maxNum=numImplantReac(i)
	end if
end do
allocate(ImplantReactions(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numImpurityReac(i) > maxNum) then
		maxNum=numImpurityReac(i)
	end if
end do
allocate(ImpurityReactions(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numSinkReac(i) > maxNum) then
		maxNum=numSinkReac(i)
	end if
end do
allocate(SinkReactions(numMaterials,maxNum))

maxNum=0
do i=1,numMaterials
	if(numDissocReac(i) > maxNum) then
		maxNum=numDissocReac(i)
	end if
end do
allocate(DissocReactions(numMaterials,maxNum))

do i=1,numMaterials
	!these subroutines (located in MaterialInput.f90) read in material parameters.
	call readMaterialInput(filename(i))
end do

!If we are including the strain field, read in the dipoles from a data file
!Must do this after readMaterialInput or numSpecies won't be defined
if(strainField=='yes') then
	call readDipoleTensors(dipoleFileName)
end if

end subroutine

!***************************************************************************************************
!
!> Subroutine readParameters() - reads in simulation parameters from parameters.txt
!!
!! This subroutine reads in all simulation and material parameters located in parameters.txt
!! as well a the file names for all other input files (material input, mesh, cascades, implantation
!! profile, etc). It also reads in toggles for various simulation options. Default values for
!! several options are stored here.
!!
!! Geometric constants used in computing clustering rates are computed at the end of this subroutine.
!
!***************************************************************************************************

subroutine readParameters()
use mod_constants
use DerivedType
implicit none

character*20 char
logical flag, flag2
integer procVol, volume

!Set default values for variables
tempStore		=273d0
CuContent		=0.5d-2
dpaRate			=1d-4
firr			=0d0
atomsize		=1.182d-2
burgers			=0.287d0
totalDPA		=1d-1
alpha_v			=1d0
alpha_i			=1d0

annealTemp		=273d0
agingTime       =0d0	!2019.04.30 Add
annealTime		=0d0
annealSteps		=0
annealType		='add'
annealTempInc	=0d0

grainBoundaryToggle	='no'
pointDefectToggle	='no'
!SIAPinToggle		='no'
polycrystal			='no'
singleElemKMC		='no'
sinkEffSearch		='no'	!< This toggle is mostly legacy code and should be set to ’no’
SIAPinMin			=1
meanFreePath		=330000
dislocationDensity	=0d0
impurityDensity		=0d0
max3DInt			=4
cascadeVolume		=0d0
numSims				=1
numGrains			=1
cascadeReactionLimit=100d0

!Toggles for various output types
postprToggle		='yes'
totdatToggle		='yes'
rawdatToggle		='no'
vtkToggle			='no'
xyzToggle			='no'
outputDebug			='no'
profileToggle		='no'

minCuCluster = 5
minVoid = 5
minLoop = 5
minCuV = 5

!Read variables in from file
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='start') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	flag2=.FALSE.
	do while(flag2 .eqv. .FALSE.)
		read(81,*) char
		if(char=='end') then
			flag2=.TRUE.
			flag=.TRUE.
		else if(char=='temperature') then
			flag2=.TRUE.
			read(81,*) tempStore
		else if(char=='CuContent') then
			flag2=.TRUE.
			read(81,*) CuContent
		else if(char=='dpaRate') then
			flag2=.TRUE.
			read(81,*) DPARate
		else if(char=='firr') then
			flag2=.TRUE.
			read(81,*) firr
		else if(char=='atomSize') then
			flag2=.TRUE.
			read(81,*) atomsize
		else if(char=='burgers') then
			flag2=.TRUE.
			read(81,*) burgers
		else if(char=='totalDPA') then
			flag2=.TRUE.
			read(81,*) totalDPA
		else if(char=='annealTemp') then
			flag2=.TRUE.
			read(81,*) annealTemp
		else if(char=='annealSteps') then
			flag2=.TRUE.
			read(81,*) annealSteps
		else if(char=='agingTime') then		!2019.04.30 Add
			flag2=.TRUE.
			read(81,*) agingTime
		else if(char=='annealTime') then
			flag2=.TRUE.
			read(81,*) annealTime
		else if(char=='annealType') then
			flag2=.TRUE.
			read(81,*) annealType
		else if(char=='annealTempInc') then
			flag2=.TRUE.
			read(81,*) annealTempInc
		else if(char=='grainBoundaries') then
			flag2=.TRUE.
			read(81,*) grainBoundaryToggle
		else if(char=='pointDefect') then
			flag2=.TRUE.
			read(81,*) pointDefectToggle
		else if(char=='grainSize') then
			flag2=.TRUE.
			read(81,*) meanFreePath
		else if(char=='dislocDensity') then
			flag2=.TRUE.
			read(81,*) dislocationDensity
		else if(char=='impurityConc') then
			flag2=.TRUE.
			read(81,*) impurityDensity
		else if(char=='max3DInt') then
			flag2=.TRUE.
			read(81,*) max3DInt
		else if(char=='cascadeVolume') then
			flag2=.TRUE.
			read(81,*) cascadeVolume
		else if(char=='cascRxnLimit') then
			flag2=.TRUE.
			read(81,*) cascadeReactionLimit
!		else if(char=='SIAPinToggle') then
!			flag2=.TRUE.
!			read(81,*) SIAPinToggle
		else if(char=='SIAPinMin') then
			flag2=.TRUE.
			read(81,*) SIAPinMin
		else if(char=='numSims') then
			flag2=.TRUE.
			read(81,*) numSims
		else if(char=='polycrystal') then
			flag2=.TRUE.
			read(81,*) polycrystal
		else if(char=='numGrains') then
			flag2=.TRUE.
			read(81,*) numGrains
!		else if(char=='vtkToggle') then
!			flag2=.TRUE.
!			read(81,*) vtkToggle
!		else if(char=='xyzToggle') then
!			flag2=.TRUE.
!			read(81,*) xyzToggle
!		else if(char=='restartToggle') then
!			flag2=.TRUE.
!			read(81,*) outputDebug
!		else if(char=='postprToggle') then
!			flag2=.TRUE.
!			read(81,*) postprToggle
!		else if(char=='totdatToggle') then
!			flag2=.TRUE.
!			read(81,*) totdatToggle
!		else if(char=='rawdatToggle') then
!			flag2=.TRUE.
!			read(81,*) rawdatToggle
!		else if(char=='profileToggle') then
!			flag2=.TRUE.
!			read(81,*) profileToggle
		else if(char=='singleElemKMC') then
			flag2=.TRUE.
			read(81,*) singleElemKMC
		else if(char=='sinkEffSearch') then
			flag2=.TRUE.
			read(81,*) sinkEffSearch
		else if(char=='alpha_v') then
			flag2=.TRUE.
			read(81,*) alpha_v
		else if(char=='alpha_i') then
			flag2=.TRUE.
			read(81,*) alpha_i
		else if(char=='conc_v') then
			flag2=.TRUE.
			read(81,*) conc_v
		else if(char=='conc_i') then
			flag2=.TRUE.
			read(81,*) conc_i
		else
			write(*,*) 'error readParameters() unrecognized parameter: ', char
		end if
	end do
	flag2=.FALSE.
end do
flag=.FALSE.

!***********************************************************************
!Output parameters
!***********************************************************************
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='OutputStart') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	flag2=.FALSE.
	do while(flag2 .eqv. .FALSE.)
		read(81,*) char
		if(char=='end') then
			flag2=.TRUE.
			flag=.TRUE.
		else if(char=='vtkToggle') then
			flag2=.TRUE.
			read(81,*) vtkToggle
		else if(char=='xyzToggle') then
			flag2=.TRUE.
			read(81,*) xyzToggle
		else if(char=='restartToggle') then
			flag2=.TRUE.
			read(81,*) outputDebug
		else if(char=='postprToggle') then
			flag2=.TRUE.
			read(81,*) postprToggle
		else if(char=='totdatToggle') then
			flag2=.TRUE.
			read(81,*) totdatToggle
		else if(char=='rawdatToggle') then
			flag2=.TRUE.
			read(81,*) rawdatToggle
		else if(char=='profileToggle') then
			flag2=.TRUE.
			read(81,*) profileToggle
		else if(char=='minCuCluster') then
			flag2=.TRUE.
			read(81,*) minCuCluster
		else if(char=='minVoid') then
			flag2=.TRUE.
			read(81,*) minVoid
		else if(char=='minLoop') then
			flag2=.TRUE.
			read(81,*) minLoop
		else if(char=='minCuV') then
			flag2=.TRUE.
			read(81,*) minCuV
		else
			write(*,*) 'error readParameters() unrecognized parameter: '
		end if
	end do
	flag2=.FALSE.
end do
flag=.FALSE.


!***********************************************************************
!if we are using adaptive meshing, read in the adaptive meshing parameters
!***********************************************************************
if(meshingType=='adaptive') then
	flag=.FALSE.
	
	do while(flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='fineStart') then
			flag=.TRUE.
		endif
	end do
	flag=.FALSE.
	
	do while(flag .eqv. .FALSE.)
		flag2=.FALSE.
		do while(flag2 .eqv. .FALSE.)
			read(81,*) char
			if(char=='end') then
				flag2=.TRUE.
				flag=.TRUE.
			else if(char=='fineLength') then
				flag2=.TRUE.
				read(81,*) fineLength
			else if(char=='numxFine') then
				flag2=.TRUE.
				read(81,*) numxCascade
			else if(char=='numyFine') then
				flag2=.TRUE.
				read(81,*) numyCascade
			else if(char=='numzFine') then
				flag2=.TRUE.
				read(81,*) numzCascade
			else
				write(*,*) 'error readParameters() unrecognized parameter'
			end if
		end do
		flag2=.FALSE.
	end do
	flag=.FALSE.
	
	cascadeElementVol=fineLength**3d0
	numCellsCascade=numxCascade*numyCascade*numzCascade
endif

!***********************************************************************
!clustering rate constants
!***********************************************************************

omega=(48d0*pi**2/atomsize**2)**(1d0/3d0) 			!clustering rate parameter for spherical clusters
omegastar=(4*pi*reactionRadius)/atomsize			!clustering rate parameter modifier due to reaction radius
omega2D=(4d0*pi/(atomsize*burgers))**(1d0/2d0)		!clustering rate parameter for 1D migrating circular clusters
omega1D=(9d0*pi/(16d0*atomsize))**(1d0/6d0)			!clustering rate parameter for 1D migrating spherical clusters
omegastar1D=reactionRadius*(pi/atomsize)**(1d0/2d0)
!omegastar1D=0d0										!clustering rate parameter modifier due to reaction radius
omegacircle1D=(1d0/burgers)**(1d0/2d0)				!clustering rate parameter for 1D migrating circular clusters

recombinationCoeff=4d0*pi*(.4466)/atomsize			!from Stoller et al., not used any longer

!***********************************************************************
!initialize counters
!***********************************************************************

!numImplantEvents=0									!counter for number of cascade /Frenkel pairs introduced
!numHeImplantEvents=0								!counter for number of He atoms introduced
!numAnnihilate=0										!counter for number of Vac that are annihilated

if(myProc%taskid==MASTER) write(*,*) 'cascadeReactionLimit', cascadeReactionLimit

end subroutine


!***************************************************************************************************
!>Subroutine reaction list sizes - finds max sizes of reaction parameters for allocation
!!
!!Information read in includes:
!!
!!1) Number of allowed defects
!!
!!2) Number of allowed reactions of each size
!***************************************************************************************************

subroutine readReactionListSizes(filename)	!< filename is 'xx_Defects.txt'
use DerivedType
use mod_constants
implicit none

integer i, j, count
character*20 :: char
character*50 :: filename
integer, allocatable :: DefectType(:), productType(:)
!type(reaction), pointer :: reactions, reactionCurrent
logical flag
Double precision Diff, Eb
integer matNum

open(80, file=filename,action='read', status='old')

flag=.FALSE.

!The following is for the entire parameter set

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='material') then
		flag=.TRUE.
		read(80,*) matNum	!< matNum = 1 or 2
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='species') then	
		flag=.TRUE.
		read(80,*) numSpecies	!< numSpecies = 4
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleForm(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='diffusionPrefactors') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleDiff(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numFunction') then
		flag=.TRUE.
		read(80,*) numFuncDiff(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='bindingEnergies') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numSingle') then
		flag=.TRUE.
		read(80,*) numSingleBind(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='numFunction') then
		flag=.TRUE.
		read(80,*) numFuncBind(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='singleDefect') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='dissociation') then
		flag=.TRUE.
		read(80,*) numDissocReac(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='diffusion') then
		flag=.TRUE.
		read(80,*) numDiffReac(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='sinkRemoval') then
		flag=.TRUE.
		read(80,*) numSinkReac(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='impurityTrapping') then
		flag=.TRUE.
		read(80,*) numImpurityReac(matNum)
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='multipleDefect') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='clustering') then
		flag=.TRUE.
		read(80,*) numClusterReac(matNum)
	end if
end  do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(80,*) char
	if(char=='noDefect') then
		flag=.TRUE.
		read(80,*) numImplantReac(matNum)
	end if
end do
flag=.FALSE.

close(80)
end subroutine
