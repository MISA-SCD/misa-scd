! $Header: /home/CVS//srscd/src/MeshReader.f90,v 1.11 2015/12/10 18:19:18 aydunn Exp $
!>Module Mesh Reader: creates a processor mesh and volume element mesh using input from file
!!
!!This module is responsible for reading the mesh from a file (in a specific format) and creating
!!myProc and myMesh(:), the information in each processor and each mesh, including connectivity and
!!the division of volume elements over processors.
!!
!!There are two formats of mesh accepted at this point: a uniform cubic mesh, and a non-uniform cubic
!!mesh. The first set of subroutines below is associated with a uniform cubic mesh, and the second
!!set is associated with a non-uniform cubic mesh. The code can use either set of subroutines to create
!!a valid myMesh(:) array, but the input file must match the rules for either uniform or nonuniform 
!!mesh input (see rules in input files).
!!
!!The processors are divided among the global volume in such a way as to minimize the surface area 
!!between processors, minimizing the amount of communication necessary between them.
!!
!!5.28.2014: NOTE: the SRSCD_par program is currently unable to accept meshes in which the processor
!!of a neighboring element in a given direction is not equal to the (global) neighboring processor in that same
!!direction (can be caused by uneven meshing and irregular division of elements)

module MeshReader
use mod_srscd_constants
use DerivedType
implicit none

contains

!***************************************************************************************************
!Subroutines for creating uniform mesh. Some subroutines are re-used in non-uniform mesh, but 
!most are different.
!***************************************************************************************************

!>Subroutine read Mesh Uniform - creates processor and mesh files from uniform mesh input
!!
!!This is the main subroutine that controls reading the uniform mesh. It reads in the mesh
!!from a file and divides the volume elements between the processors in the parallel
!!simulation. If a division such that each processor has at least one element is not possible,
!!then the subroutine returns an error. It also creates a 'processor mesh', indicating
!!the bounds of each processor. NOTE: this subroutine is for uniform meshes only, a different
!!subroutine has been created to read in non-uniform meshes (and these have not been 
!!implemented in the rest of the program)
!!
!!Input: filename of mesh file
!!
!!Outputs: myProc (processors with meshes), myMesh, and myBoundary

subroutine readMeshUniform(filename)
use mod_srscd_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), procDivision(3), i, j, k, l, numx, numy, numz, numTotal, matNum
character*20 readIn, meshType
character*50 filename, filename2, filename3
logical flag
double precision volumeFaces(3), totalArea, length, tempCoord(3)
integer element, localElements(3), localElem, numxLocal, numyLocal, numzLocal, globalCell, globalNeighbor, maxElement
integer tracker
double precision, allocatable :: globalMeshCoord(:,:)
double precision, allocatable :: globalStrain(:,:)
integer, allocatable :: globalMeshConnect(:,:), globalMaterial(:)

interface

	subroutine createConnectLocalPeriodic(numx, numy, numz, globalMeshCoord, globalMeshConnect)
	use DerivedType
	integer numx, numy, numz
	double precision, allocatable :: globalMeshCoord(:,:)
	integer, allocatable :: globalMeshConnect(:,:)
	end subroutine
	
	subroutine createConnectLocalFreeSurf(numx, numy, numz, globalMeshCoord, globalMeshConnect)
	use DerivedType
	integer numx, numy, numz
	double precision, allocatable :: globalMeshCoord(:,:)
	integer, allocatable :: globalMeshConnect(:,:)
	end subroutine

end interface

open(80, file=filename,action='read', status='old')
if(strainField=='yes') then
	open(50, file=strainFileName, action='read', status='old')
end if

!Used to read in data from files
flag=.FALSE.

!Step -1: read in mesh type (free surfaces or periodic)
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='meshType') then
		read(80,*) meshType
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!Step 0: read in volume element lengths (true for UNIFORM CUBIC MESH ONLY)
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='length') then
		read(80,*) length
		flag=.TRUE.
	end if
end do
flag=.FALSE.

meshLength = length

!************************************************
!Step 1: read in the global coordinates (max, min)
!************************************************

!The coordinates are the center of volume elements
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='xminmax') then
		read(80,*) myProc%globalCoord(1), myProc%globalCoord(2)
	else if(readIn=='yminmax') then
		read(80,*) myProc%globalCoord(3), myProc%globalCoord(4)
	else if(readIn=='zminmax') then
		read(80,*) myProc%globalCoord(5), myProc%globalCoord(6)
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!***********************************************************************
!Read in the strain field into myMesh%strain(1-6)
!***********************************************************************

if(strainField=='yes') then
	
	!Step -1: skip to start of input data
	do while(flag .eqv. .FALSE.)
		read(50,*) readIn
		if(readIn=='start') then
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	read(50,*)	!Blank line in strain input file
	
end if

!Step 1a: modify global coordinates to include entire simulation volume (not just centers of elements)
!The actual boundary coordinates of the system
do i=1,5,2
	myProc%globalCoord(i)=myProc%globalCoord(i)-length/2d0
	myProc%globalCoord(i+1)=myProc%globalCoord(i+1)+length/2d0
end do

!Debug options
!if(myProc%taskid==MASTER) then
!	write(*,*) 'global coordinates'
!	write(*,*) myProc%globalCoord
!endif

!************************************************
!Step 2: divide processors over gobal volume
!************************************************

!Step 2a: find area of each face of entire volume (x-y, y-z, x-z)
volumeFaces(1)=(myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))	!x-y
volumeFaces(2)=(myProc%globalCoord(4)-myProc%globalCoord(3))*(myProc%globalCoord(6)-myProc%globalCoord(5))	!y-z
volumeFaces(3)=(myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(6)-myProc%globalCoord(5))	!x-z
	
!Step 2b: for each factorization of the number of processors, find the total shared area and minimize
totalArea=dble(myProc%numtasks)*(volumeFaces(1)+volumeFaces(2)+volumeFaces(3))	!This is an upper bound estimate

procDivision(1)=0	!number of processors in x
procDivision(2)=0	!number of processors in y
procDivision(3)=0	!number of processors in z
!Keep factorizations of processors that minimize shared area.
do i=1,myProc%numtasks
	do j=1, myProc%numtasks
		do k=1, myProc%numtasks
			if(i*j*k==myProc%numtasks) then
				!find shared area between processors using this factorization of the number of processors
				if(k*volumeFaces(1)+i*volumeFaces(2)+j*volumeFaces(3) < totalArea .OR. myProc%numtasks==1) then
					procDivision(1)=i
					procDivision(2)=j
					procDivision(3)=k
					totalArea=k*volumeFaces(1)+i*volumeFaces(2)+j*volumeFaces(3)
					!write(*,*) 'new division found', procDivision, totalArea
				end if
				!do nothing, not a factorization of numtasks
			end if
		end do
	end do
end do

if(myProc%taskid==MASTER) then
	write(*,*) 'proc division', procDivision
end if

!step 2c: divide volume among processors according to the factorization above
!**********
!This step locates the boundaries of the local processor within the global volume.
!**********
myProc%localCoord(1)=myProc%globalCoord(1)+&
	mod(myProc%taskid,procDivision(1))*(myProc%globalCoord(2)-myProc%globalCoord(1))/&
	dble(procDivision(1))

myProc%localCoord(2)=myProc%localCoord(1)+&
	(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(procDivision(1))

myProc%localCoord(3)=myProc%globalCoord(3)+&
	mod(myProc%taskid/ProcDivision(1),procDivision(2))*&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))

myProc%localCoord(4)=myProc%localCoord(3)+&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))

myProc%localCoord(5)=myProc%GlobalCoord(5)+&
	mod(myProc%taskid/(ProcDivision(1)*ProcDivision(2)),procDivision(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))

myProc%localCoord(6)=myProc%localCoord(5)+&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))
	
totalVolume=(myProc%localCoord(2)-myProc%localCoord(1))*(myProc%localCoord(4)-myProc%localCoord(3))*&
	(myProc%localCoord(6)-myProc%localCoord(5))

!write(*,*) 'proc', myProc%taskid, 'of', myProc%numtasks
!write(*,*) 'local coordinates', myProc%localCoord

!Step 2d: point processor myProc%procNeighbor(x) to neighboring processors
!********
!This is creating a 'processor mesh' (not the actual mesh) and a connectivity between processors
!based on how the total volume was divided in the previous steps. The connectivity
!in this step is crated by iterating x, then y, then z
!
!This mesh is periodic and does not account for the possibility of free surfaces.
!********
if (myProc%localCoord(1)==myProc%globalCoord(1)) then	!coordinate is at xmin
	myProc%procNeighbor(2)=myProc%taskid+procDivision(1)-1	!left
else
	myProc%procNeighbor(2)=myProc%taskid-1	!left
end if

if (myProc%localCoord(2)==myProc%globalCoord(2)) then	!coordinate is at xmax
	myProc%procNeighbor(1)=myProc%taskid-procDivision(1)+1	!right
else
	myProc%procNeighbor(1)=myProc%taskid+1	!right
end if

if (myProc%localCoord(3)==myProc%globalCoord(3)) then	!coordinate is at ymin
	myProc%procNeighbor(4)=myProc%taskid+procDivision(1)*(procDivision(2)-1)	!back
else
	myProc%procNeighbor(4)=myProc%taskid-procDivision(1)	!back
end if

if (myProc%localCoord(4)==myProc%globalCoord(4)) then	!coordinate is at ymax
	myProc%procNeighbor(3)=myProc%taskid-procDivision(1)*(procDivision(2)-1)	!front
else
	myProc%procNeighbor(3)=myProc%taskid+procDivision(1)	!front
end if

if (myProc%localCoord(5)==myProc%globalCoord(5)) then	!coordinate is at zmin
	myProc%procNeighbor(6)=myProc%taskid+procDivision(1)*procDivision(2)*(procDivision(3)-1)	!down
else
	myProc%procNeighbor(6)=myProc%taskid-procDivision(1)*procDivision(2)	!down
end if

if (myProc%localCoord(6)==myProc%globalCoord(6)) then	!coordinate is at zmax
	myProc%procNeighbor(5)=myProc%taskid-procDivision(1)*procDivision(2)*(procDivision(3)-1)	!up
else
	myProc%procNeighbor(5)=myProc%taskid+procDivision(1)*procDivision(2)	!up
end if

!write(*,*) 'proc', myProc%taskid, 'of', myProc%numtasks
!write(*,*) 'neighbors', myProc%procNeighbor

!Step 3: create mesh in each processor of elements only within that processor's local coordinates
!NOTE: this is an actual mesh of volume elements, not the 'processor mesh' created above.

!Step 3a: find how many volume elements are in global system and which of those are in the local processor's range

!Read the number of elements in x,y,and z directions
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='numx') then
		read(80,*) numx
	else if(readIn=='numy') then
		read(80,*) numy
	else if(readIn=='numz') then
		read(80,*) numz
		flag=.TRUE.
	end if
end do
flag=.FALSE.

numTotal=numx*numy*numz	!total cell in the system
totalMesh = numTotal
totalX = numx
totalY = numy
totalZ = numz

!These arrays create a global mesh and global list of material numbers and coordinates, but are discarded
!once the local mesh is finished.
allocate(globalMaterial(numTotal))
allocate(globalMeshCoord(numTotal,3))
allocate(globalMeshConnect(numTotal,6))
allocate(globalStrain(numTotal,6))	!first 3 coordinates are mesh coordinates, last 6 coordinates are strain tensor


!Create global connectivity (uniform cubic mesh - use same rule as above, count by x then y then z)
if(meshType=='periodic') then
	call createConnectGlobalPeriodicUniform(globalMeshConnect, numx, numy, numz)
else if(meshType=='freeSurfaces') then
	call createConnectGlobalFreeSurfUniform(globalMeshConnect, numx, numy, numz)
end if

!Tells program that we are about to start reading in element coordinates and material numbers
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='elements') then
		flag=.TRUE.
	endif
end do
flag=.FALSE.

!NumxLocal, numyLocal, numzLocal are used to determine the size of the mesh inside the local processor.
!They are calculated during the next step.
numxLocal=0
numyLocal=0
numzLocal=0
element=1
systemVol=0d0

localElements(3)=0
do k=1,numz
	localElements(2)=0
	do j=1,numy
		localElements(1)=0
		do i=1,numx
			!read in the coordinates and material number from file of this element (coordinates are the center of the element)
			read(80,*) globalMeshCoord(element,1),globalMeshCoord(element,2),globalMeshCoord(element,3),&
				globalMaterial(element)
				
			if(strainField=='yes') then
			
				read(50,*) tempCoord(1), tempCoord(2), tempCoord(3), globalStrain(element,1), globalStrain(element,2), &
					globalStrain(element,3), globalStrain(element,4), globalStrain(element,5), globalStrain(element,6)
					
				if(tempCoord(1) /= globalMeshCoord(element,1)) then	!not equal
					write(*,*) 'Error strain mesh does not match normal mesh'
				else if(tempCoord(2) /= globalMeshcoord(element,2)) then
					write(*,*) 'Error strain mesh does not match normal mesh'
				else if(tempCoord(3) /= globalMeshCoord(element,3)) then
					write(*,*) 'Error strain mesh does not match normal mesh'
				end if
			
			end if
				
			if(numMaterials==1) then
				systemVol=systemVol+length**3d0
			else if(globalMaterial(element)==1) then	!only add to system volume if we are NOT at a grain boundary
				systemVol=systemVol+length**3d0
			else
				!Do nothing
			end if
			
			!If this element is within the local bounds of this processor, as defined in step 2
			if(globalMeshCoord(element,1) > myProc%localCoord(1) .AND. globalMeshCoord(element,1) <= myProc%localCoord(2) &
				.AND. globalMeshCoord(element,2) > myProc%localCoord(3) .AND. globalMeshCoord(element,2) <= myProc%localCoord(4) &
				.AND. globalMeshCoord(element,3) > myProc%localCoord(5) .AND. globalMeshCoord(element,3) <= myProc%localCoord(6)) then
				
				!used to count the number of elements in the x-, y-, and z-directions within each processor
				localElements(1)=localElements(1)+1			!add elements in the x-direction
				if(localElements(1)==1) then
					localElements(2)=localElements(2)+1		!create another row in y-direction
				end if
				if(localElements(2)==1 .AND. localElements(1)==1) then
					localElements(3)=localElements(3)+1		!create another row in the z-direction
				end if
				
				!Here, we keep increasing numxLocal, numyLocal, numzLocal until the entire local mesh
				!is read in
				
				if(localElements(1) > numxLocal) then
					numxLocal=localElements(1)
				end if
				if(localElements(2) > numyLocal) then
					numyLocal=localElements(2)
				end if
				if(localElements(3) > numzLocal) then
					numzLocal=localElements(3)
				end if
			end if
			
			element=element+1
		end do
	end do
end do

!If any processors don't have volume elements in them (too many procs), we create an error message
if(numxLocal*numyLocal*numzLocal==0) then
	write(*,*) 'error processors with no volume elements'
	call MPI_ABORT(MPI_COMM_WORLD,ierr)
	stop
end if

!Step 3b: find how many volume elements are in local processor and allocate myMesh accordingly

write(*,*) 'proc', myProc%taskid, 'numx numy numz', numxLocal, numyLocal, numzLocal

allocate(myMesh(numxLocal*numyLocal*numzLocal))

!this variable is used at various points in SRSCD; lets us know the length fo myMesh(:)
numCells=numxLocal*numyLocal*numzLocal	!total cells of this processor
localX=numxLocal
localY=numyLocal
localZ=numzLocal

do i=1,numCells
	myMesh(i)%length=length
	myMesh(i)%volume=length**3d0
	
	!uniform mesh: all elements have 1 neighbor in each direction
	allocate(myMesh(i)%neighbors(6,1))
	allocate(myMesh(i)%neighborProcs(6,1))
	do j=1,6
		myMesh(i)%numNeighbors(j)=1
	end do

end do

!Step 3c: for each volume element in myMesh, assign coordinates and material number (and proc number)
!We can't read from the file directly into myMesh, because we don't know how big to make it until
!we read in all of the volume element information. Therefore, we allocate it and then reread from the
!global mesh back into myMesh.
element=1
localElem=1
do k=1,numz
	do j=1, numy
		do i=1, numx
			if(globalMeshCoord(element,1) > myProc%localCoord(1) .AND. globalMeshCoord(element,1) <= myProc%localCoord(2) &
				.AND. globalMeshCoord(element,2) > myProc%localCoord(3) .AND. globalMeshCoord(element,2) <= myProc%localCoord(4) &
				.AND. globalMeshCoord(element,3) > myProc%localCoord(5) .AND. globalMeshCoord(element,3) <= myProc%localCoord(6)) then
				
				myMesh(localElem)%coordinates(1)=globalMeshCoord(element,1)
				myMesh(localElem)%coordinates(2)=globalMeshCoord(element,2)
				myMesh(localElem)%coordinates(3)=globalMeshCoord(element,3)
				myMesh(localElem)%material=globalMaterial(element)
				myMesh(localElem)%proc=myProc%taskid
				
				if(strainField=='yes') then
					do l=1,6
						myMesh(localElem)%strain(l)=globalStrain(element,l)
					end do
				end if
				
				localElem=localElem+1
			end if
			element=element+1
		end do
	end do
end do

!Step 3d: assign neighbors and processor numbers for neighbors (connectivity in myMesh) - periodic or free surfaces in +/- z

if(meshType=='periodic') then
	call createConnectLocalPeriodicUniform(numxLocal, numyLocal, numzLocal, globalMeshCoord, globalMeshConnect)
else if(meshType=='freeSurfaces') then
	call createConnectLocalFreeSurfUniform(numxLocal, numyLocal, numzLocal, globalMeshCoord, globalMeshConnect)
end if

!This output is optional, used to show the user what the mesh is

!if(myProc%taskid==4) then
!	do 25 i=1,numxLocal*numyLocal*numzLocal
!		write(*,*) 'element', i, 'length', myMesh(i)%length, 'material', myMesh(i)%material
!		write(*,*) 'coords', (myMesh(i)%coordinates(j),j=1,3)
!		write(*,*) 'numNeighbors', (myMesh(i)%numNeighbors(j),j=1,6)
!		write(*,*) 'neighbors and procs'
!		do 26 j=1,6
!			write(*,*) (myMesh(i)%neighbors(j,k),k=1,myMesh(i)%numNeighbors(j))
!			write(*,*) (myMesh(i)%neighborProcs(j,k),k=1,myMesh(i)%numNeighbors(j))
!			write(*,*)
!		26 continue
!	25 continue
!	
!	write(*,*) 'connectivity order +x, -x, +y, -y, +z, -z'
!	write(*,*) 'proc division', procDivision, 'Division order: x, y, z'
!	write(*,*)
!	write(*,*) 'Global coordinates (xmin, xmax, ymin, ymax, zmin, zmax)'
!	write(*,*) myProc%globalCoord
!	write(*,*)
!endif

!***************************************************************************************************
!5/28/2014: initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************
!Using the following global mesh info:
!
!double precision, allocatable :: globalMeshCoord(:,:)
!integer, allocatable :: globalMeshConnect(:,:), globalMaterial(:)

!Step 1: Find the max cell# of any boundary mesh element
maxElement=0
do i=1, numxLocal*numyLocal*numzLocal
	do j=1, 6
		if(myMesh(i)%neighborProcs(j,1) /= myProc%taskid) then	!we are pointing to a different proc
			if(myMesh(i)%neighbors(j,1) > maxElement) then		!searching for the max element number in a neighbor
				maxElement=myMesh(i)%neighbors(j,1)
			end if
		end if
	end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(6,maxElement))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do i=1,maxElement
	do j=1,6
		myBoundary(j,i)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
		myBoundary(j,i)%proc=-10		!default, says that this is not a real element of myBoundary.
	end do
end do

do i=1,numxLocal*numyLocal*numzLocal
	do j=1,6
		if(myMesh(i)%neighborProcs(j,1) == -1) then										!this is a free surface
			!do nothing
		else if(myMesh(i)%neighborProcs(j,1) /= myProc%taskid) then
			myBoundary(j,myMesh(i)%neighbors(j,1))%proc=myMesh(i)%neighborProcs(j,1)	!set proc # of elements in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%length=length						!set length of elements in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%volume=length**3d0					!set volume of elements in myBoundary (changes with cascade addition)
			globalCell=findGlobalCell(myMesh(i)%coordinates,globalMeshCoord)			!find global cell # of element in myBoundary
			globalNeighbor=globalMeshConnect(globalCell,j)								!use global cell # to find material # of element in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%material=globalMaterial(globalNeighbor)	!set material # of elements in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%localNeighbor=i
			
			if(strainField=='yes') then
				do l=1,6
					myBoundary(j,myMesh(i)%neighbors(j,1))%strain(l)=globalStrain(globalNeighbor,l)
				end do
			end if
		end if
	end do
end do

!Close input files
if(strainField=='yes') then
	close(50)
end if
close(80)

end subroutine

!These subroutines creates a global connectivity matrix (elements and their neighbors) based
!on the connectivity scheme of increasing x, then y, then z. (uniform cubic)

!>Subroutine create global connectivity (uniform mesh, periodic boundary conditions)
!!
!!This subroutine creates a global mesh (not dividing between processors) as well as the 
!!connectivity between elements for the case of a uniform cubic mesh.
!!
!!Inputs: numx, numy, numz
!!Output: connectivity

subroutine createConnectGlobalPeriodicUniform(connectivity, numx, numy, numz)
implicit none

integer cell, numElements, numx, numy, numz
integer connectivity(numx*numy*numz, 6)

numElements=numx*numy*numz	!total cells of the system
!************************************************
!periodic boundary condition version
!************************************************
do cell=1,numElements
	if(mod(cell,numx)==0) then !identify cell to the right
		connectivity(cell, 1)=cell-numx+1
	else
		connectivity(cell,1)=cell+1
	end if
	
	if(mod(cell+numx-1,numx)==0) then !identify cell to the left
		connectivity(cell,2)=cell+numx-1
	else
		connectivity(cell,2)=cell-1
	end if
	
	if(mod(cell,numx*numy) > numx*(numy-1) .OR. mod(cell,numx*numy)==0) then
		connectivity(cell,3)=cell-(numx*(numy-1))
	else
		connectivity(cell,3)=cell+numx
	end if
	
	if(mod(cell,numx*numy) <= numx .AND. (mod(cell, numx*numy) /= 0 .OR. numy==1)) then
		connectivity(cell,4)=cell+(numx*(numy-1))
	else
		connectivity(cell,4)=cell-numx
	end if
	
	if(mod(cell,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(cell, numx*numy*numz)==0) then
		connectivity(cell,5)=cell-(numx*numy*(numz-1))
	else
		connectivity(cell,5)=cell+numx*numy
	end if
	
	if(mod(cell,numx*numy*numz) <= numx*numy .AND. (mod(cell,numx*numy*numz) /= 0 .OR. numz==1)) then
		connectivity(cell,6)=cell+(numx*numy*(numz-1))
	else
		connectivity(cell,6)=cell-numx*numy
	end if
end do

end subroutine


!>Subroutine create global connectivity (uniform mesh, free surfaces in z-directions and periodic boundary conditions in other directions)
!!
!!This subroutine creates a global mesh (not dividing between processors) as well as the 
!!connectivity between elements for the case of a uniform cubic mesh with free surfaces.
!!
!!Inputs: numx, numy, numz
!!Output: connectivity

subroutine createConnectGlobalFreeSurfUniform(connectivity, numx, numy, numz)
implicit none

integer cell, numElements, numx, numy, numz
integer connectivity(numx*numy*numz, 6)

numElements=numx*numy*numz

!************************************************
!PBCs in x and y, free in z (cell 0 represents free surface)
!************************************************
do cell=1,numElements
	if(mod(cell,numx)==0) then !identify cell to the right
		connectivity(cell, 1)=cell-numx+1
	else
		connectivity(cell,1)=cell+1
	end if
	
	if(mod(cell+numx-1,numx)==0) then !identify cell to the left
		connectivity(cell,2)=cell+numx-1
	else
		connectivity(cell,2)=cell-1
	end if
	
	if(mod(cell,numx*numy) > numx*(numy-1) .OR. mod(cell,numx*numy)==0) then
		connectivity(cell,3)=cell-(numx*(numy-1))
	else
		connectivity(cell,3)=cell+numx
	end if
	
	if(mod(cell,numx*numy) <= numx .AND. (mod(cell, numx*numy) /= 0 .OR. numy==1)) then
		connectivity(cell,4)=cell+(numx*(numy-1))
	else
		connectivity(cell,4)=cell-numx
	end if
	
	if(mod(cell,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(cell, numx*numy*numz)==0) then
		connectivity(cell,5)=0
	else
		connectivity(cell,5)=cell+numx*numy
	end if
	
	if(mod(cell,numx*numy*numz) <= numx*numy .AND. (mod(cell,numx*numy*numz) /=0 .OR. numz==1)) then
		connectivity(cell,6)=0
	else
		connectivity(cell,6)=cell-numx*numy
	end if
end do

end subroutine

!These subroutines create LOCAL connectivity. They identify the volume element # and processor # of
!neighboring volume elements for each element. The connectivity scheme is the same as in the global
!case, but neighboring processor numbers are used here.

!>Subroutine create local connectivity (uniform mesh, periodic boundary conditions)
!!
!!This subroutine creates a local mesh (for the local processor only) as well as the 
!!connectivity between elements for the case of a uniform cubic mesh. This subroutine
!!has to choose which elements of the global mesh are in the local mesh, and puts those in
!!myMesh along with the element coordinates and their neighbor element ID numbers and processors.
!!
!!Inputs: numx, numy, numz
!!Inputs: numx, numy, numz, global mesh coordinates, global mesh connectivity
!!Output: myMesh

subroutine createConnectLocalPeriodicUniform(numx, numy, numz, globalMeshCoord, globalMeshConnect)
use DerivedType
use mod_srscd_constants

implicit none
include 'mpif.h'
integer numx, numy, numz, cell, numFaces(6)
double precision, allocatable :: globalMeshCoord(:,:)
integer, allocatable :: globalMeshConnect(:,:)
integer globalCell, globalNeighbor, status(MPI_STATUS_SIZE)

!buffer lists to send all information at the end
integer numSendRecv, i, dir, sendList(6*numx*numy*numz, 4)

numSendRecv=0

!************************************************
!periodic boundary condition version
!************************************************
do cell=1,numCells
	if(mod(cell,numx)==0) then !identify cell to the right
		
		!If we are on the right edge of the local mesh, identify the neighboring processor
		myMesh(cell)%neighborProcs(1,1)=myProc%procNeighbor(1)
		
		!If the mesh is only one processor thick, then the neighboring processor is the same as this processor
		if(myMesh(cell)%neighborProcs(1,1)==myProc%taskid) then	
			myMesh(cell)%neighbors(1,1)=cell-numx+1	!use periodic rules from uniform cubic mesh
		else
			!find global cell number of this cell and the neighboring cell.
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,1)
			
			!add these items to sendList, a buffer that is used at the end of this subroutine to communicate
			!with other processors via MPI about which elements have neighbors in other cells
			!(this cell sends information about itself, and the nieghboring cell will do the same in its
			!processor)
			
			numSendRecv=numSendRecv+1		!count the number of items in buffer
			sendList(numSendRecv,1)=cell	!localCellID in this processor
			sendList(numSendRecv,2)=1		!direction: right
			sendList(numSendRecv,3)=globalCell	!the globalCellID of the localCellID
			sendList(numSendRecv,4)=globalNeighbor	!neighbor in the right
		end if
	else
		!if we are still inside the local mesh, don't need to communicate with neighboring cells and 
		!just use the uniform cubic connectivity rules (increase x, then y, then z)
		myMesh(cell)%neighbors(1,1)=cell+1
		myMesh(cell)%neighborProcs(1,1)=myProc%taskid
	end if
	
	!This is the same as above, except that we are pointed to the left instead of the right of the local mesh.
	if(mod(cell+numx-1,numx)==0) then !identify cell to the left
		
		myMesh(cell)%neighborProcs(2,1)=myProc%procNeighbor(2)
		if(myMesh(cell)%neighborProcs(2,1)==myProc%taskid) then
			myMesh(cell)%neighbors(2,1)=cell+numx-1
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,2)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=2	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(2,1)=cell-1
		myMesh(cell)%neighborProcs(2,1)=myProc%taskid
	end if
	
	!Front (+y)
	if(mod(cell,numx*numy) > numx*(numy-1) .OR. mod(cell,numx*numy)==0) then
		
		myMesh(cell)%neighborProcs(3,1)=myProc%procNeighbor(3)
		if(myMesh(cell)%neighborProcs(3,1)==myProc%taskid) then
			myMesh(cell)%neighbors(3,1)=cell-(numx*(numy-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates, globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,3)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=3	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(3,1)=cell+numx
		myMesh(cell)%neighborProcs(3,1)=myProc%taskid
	end if
	
	!Back (-y)
	if(mod(cell,numx*numy) <= numx .AND. (mod(cell, numx*numy) /= 0 .OR. numy==1)) then
		myMesh(cell)%neighborProcs(4,1)=myProc%procNeighbor(4)
		if(myMesh(cell)%neighborProcs(4,1)==myProc%taskid) then
			myMesh(cell)%neighbors(4,1)=cell+(numx*(numy-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,4)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=4	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(4,1)=cell-numx
		myMesh(cell)%neighborProcs(4,1)=myProc%taskid
	end if
	
	!Up (+z)
	if(mod(cell,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(cell, numx*numy*numz)==0) then
		myMesh(cell)%neighborProcs(5,1)=myProc%procNeighbor(5)
		if(myMesh(cell)%neighborProcs(5,1)==myProc%taskid) then
			myMesh(cell)%neighbors(5,1)=cell-(numx*numy*(numz-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,5)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=5	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(5,1)=cell+numx*numy
		myMesh(cell)%neighborProcs(5,1)=myProc%taskid
	end if
	
	!Down (-z)
	if(mod(cell,numx*numy*numz) <= numx*numy .AND. (mod(cell,numx*numy*numz) /= 0 .OR. numz==1)) then
		myMesh(cell)%neighborProcs(6,1)=myProc%procNeighbor(6)
		if(myMesh(cell)%neighborProcs(6,1)==myProc%taskid) then
			myMesh(cell)%neighbors(6,1)=cell+(numx*numy*(numz-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,6)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=6	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(6,1)=cell-numx*numy
		myMesh(cell)%neighborProcs(6,1)=myProc%taskid
	end if
end do

!Now that we have created SendList, the buffer that contains information about which elements in myMesh
!have neighbors on other processors, we send out the cell numbers of these cells to the neighboring procs.
!We label each MPI_SEND with globalCell so that the receiving processor pairs cells correctly using globalNeighbor (see below)

do i=1,numSendRecv
	cell=sendList(i,1)
	dir=sendList(i,2)
	globalCell=sendList(i,3)
	call MPI_SEND(cell, 1, MPI_INTEGER, myMesh(cell)%neighborProcs(dir,1), globalCell, MPI_COMM_WORLD, ierr)
end do

!Second half of send/recieve routine. Here we recieve the neighbor (we sent the cell we were in). 
!We label the MPI_RECV with globalNeighbor to match it to globalCell in the prev. step.
	
do i=1,numSendRecv
	cell=sendList(i,1)
	dir=sendList(i,2)
	globalNeighbor=sendList(i,4)
	call MPI_RECV(myMesh(cell)%neighbors(dir,1), 1, MPI_INTEGER, myMesh(cell)%neighborProcs(dir,1), &
				globalNeighbor, MPI_COMM_WORLD, status, ierr)
end do

end subroutine

!This is the same as the previous subroutine except if the cell is at the free surface, its connectivity cell number
!is 0 and the processor number is -1.

!>Subroutine create local connectivity (uniform mesh, free surfaces in z-direction and periodic boundary conditions in other directions)
!!
!!This subroutine creates a local mesh (for the local processor only) as well as the 
!!connectivity between elements for the case of a uniform cubic mesh with free surfaces. This subroutine
!!has to choose which elements of the global mesh are in the local mesh, and puts those in
!!myMesh along with the element coordinates and their neighbor element ID numbers and processors.
!!
!!Inputs: numx, numy, numz, global mesh coordinates, global mesh connectivity
!!Output: myMesh

subroutine createConnectLocalFreeSurfUniform(numx, numy, numz, globalMeshCoord, globalMeshConnect)
use DerivedType

implicit none
include 'mpif.h'
integer numx, numy, numz, cell, numFaces(6)
double precision, allocatable :: globalMeshCoord(:,:)
integer, allocatable :: globalMeshConnect(:,:)
integer globalCell, globalNeighbor, status(MPI_STATUS_SIZE)

!buffer lists to send all information at the end
integer numSendRecv, i, dir, sendList(6*numx*numy*numz, 4)

numSendRecv=0

!******************************************************************
!free surfaces at z=0 and z=zmax(Global) boundary condition version
!******************************************************************
do cell=1,numCells
	if(mod(cell,numx)==0) then !identify cell to the right
		
		myMesh(cell)%neighborProcs(1,1)=myProc%procNeighbor(1)
		if(myMesh(cell)%neighborProcs(1,1)==myProc%taskid) then
			myMesh(cell)%neighbors(1,1)=cell-numx+1
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,1)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=1	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(1,1)=cell+1
		myMesh(cell)%neighborProcs(1,1)=myProc%taskid
	end if
	
	if(mod(cell+numx-1,numx)==0) then !identify cell to the left
		
		myMesh(cell)%neighborProcs(2,1)=myProc%procNeighbor(2)
		if(myMesh(cell)%neighborProcs(2,1)==myProc%taskid) then
			myMesh(cell)%neighbors(2,1)=cell+numx-1
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,2)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=2	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(2,1)=cell-1
		myMesh(cell)%neighborProcs(2,1)=myProc%taskid
	end if
	
	if(mod(cell,numx*numy) > numx*(numy-1) .OR. mod(cell,numx*numy)==0) then
		
		myMesh(cell)%neighborProcs(3,1)=myProc%procNeighbor(3)
		if(myMesh(cell)%neighborProcs(3,1)==myProc%taskid) then
			myMesh(cell)%neighbors(3,1)=cell-(numx*(numy-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,3)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=3	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(3,1)=cell+numx
		myMesh(cell)%neighborProcs(3,1)=myProc%taskid
	end if
	
	
	if(mod(cell,numx*numy) <= numx .AND. (mod(cell, numx*numy) /= 0 .OR. numy==1)) then
		myMesh(cell)%neighborProcs(4,1)=myProc%procNeighbor(4)
		if(myMesh(cell)%neighborProcs(4,1)==myProc%taskid) then
			myMesh(cell)%neighbors(4,1)=cell+(numx*(numy-1))
		else
			globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
			globalNeighbor=globalMeshConnect(globalCell,4)
			numSendRecv=numSendRecv+1
			sendList(numSendRecv,1)=cell
			sendList(numSendRecv,2)=4	!direction
			sendList(numSendRecv,3)=globalCell
			sendList(numSendRecv,4)=globalNeighbor
		end if
	else
		myMesh(cell)%neighbors(4,1)=cell-numx
		myMesh(cell)%neighborProcs(4,1)=myProc%taskid
	end if
	
	if(mod(cell,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(cell, numx*numy*numz)==0) then
		myMesh(cell)%neighborProcs(5,1)=myProc%procNeighbor(5)
		globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
		globalNeighbor=globalMeshConnect(globalCell,5)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(5,1)=0
			myMesh(cell)%neighborProcs(5,1)=-1
		else
			if(myMesh(cell)%neighborProcs(5,1)==myProc%taskid) then
				myMesh(cell)%neighbors(5,1)=cell-(numx*numy*(numz-1))
			else
				numSendRecv=numSendRecv+1
				sendList(numSendRecv,1)=cell
				sendList(numSendRecv,2)=5	!direction
				sendList(numSendRecv,3)=globalCell
				sendList(numSendRecv,4)=globalNeighbor
			end if
		end if
	else
		myMesh(cell)%neighbors(5,1)=cell+numx*numy
		myMesh(cell)%neighborProcs(5,1)=myProc%taskid
	end if
	
	if(mod(cell,numx*numy*numz) <= numx*numy .AND. (mod(cell,numx*numy*numz) /= 0 .OR. numz==1)) then
		myMesh(cell)%neighborProcs(6,1)=myProc%procNeighbor(6)
		globalCell=findGlobalCell(myMesh(cell)%coordinates,globalMeshCoord)
		globalNeighbor=globalMeshConnect(globalCell,6)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(6,1)=0
			myMesh(cell)%neighborProcs(6,1)=-1
		else
			if(myMesh(cell)%neighborProcs(6,1)==myProc%taskid) then
				myMesh(cell)%neighbors(6,1)=cell+(numx*numy*(numz-1))
			else
				numSendRecv=numSendRecv+1
				sendList(numSendRecv,1)=cell
				sendList(numSendRecv,2)=6	!direction
				sendList(numSendRecv,3)=globalCell
				sendList(numSendRecv,4)=globalNeighbor
			end if
		end if
	else
		myMesh(cell)%neighbors(6,1)=cell-numx*numy
		myMesh(cell)%neighborProcs(6,1)=myProc%taskid
	end if
end do

do i=1,numSendRecv
	cell=sendList(i,1)
	dir=sendList(i,2)
	globalCell=sendList(i,3)
	call MPI_SEND(cell, 1, MPI_INTEGER, myMesh(cell)%neighborProcs(dir,1), globalCell, MPI_COMM_WORLD, ierr)
end do
	
do i=1,numSendRecv
	cell=sendList(i,1)
	dir=sendList(i,2)
	globalNeighbor=sendList(i,4)
	call MPI_RECV(myMesh(cell)%neighbors(dir,1), 1, MPI_INTEGER, myMesh(cell)%neighborProcs(dir,1), &
				globalNeighbor, MPI_COMM_WORLD, status, ierr)
end do

end subroutine

!This finds the global cell number of a cell given its coordinates and the global coordinate list

!>Subroutine find Global Cell
!!
!!This subroutine searches through a list of global cell coordinates and identifies which 
!!global cell ID number the coordinates provided corresponds to
!!
!!Inputs: coordinates of cell we are searching for, list of global cell coordinates
!!Output: cell ID

integer function findGlobalCell(coord, gCoord)
implicit none

double precision coord(3)
double precision, allocatable :: gCoord(:,:)
integer i
logical flag

flag=.FALSE.
i=0
do 10 while (flag .eqv. .FALSE.)
	i=i+1
	if(coord(1)==gCoord(i,1) .AND. coord(2)==gCoord(i,2) .AND. coord(3)==gCoord(i,3)) then
		flag=.TRUE.
	endif
10 continue

findGlobalCell=i
end function


!***************************************************************************************************
!subroutines for creating non-uniform (cubic) mesh.
!***************************************************************************************************

!>Subroutine read Mesh Non Uniform - creates processor and mesh files from non-uniform mesh input
!!
!!This is the main subroutine that controls reading the non-uniform mesh. It reads in the mesh
!!from a file and divides the volume elements between the processors in the parallel
!!simulation. If a division such that each processor has at least one element is not possible,
!!then the subroutine returns an error. It also creates a 'processor mesh', indicating
!!the bounds of each processor. NOTE: this subroutine is for non-uniform meshes only (and these have not been 
!!implemented in the rest of the program)
!!
!!Input: filename of mesh file
!!
!!Outputs: myProc (processors with meshes), myMesh, and myBoundary

!main subroutine that creates mesh (non-uniform cubic)
subroutine readMeshNonUniform(filename)
use mod_srscd_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), procDivision(3), i, j, k, numx, numy, numz, numTotal, matNum
character*20 readIn, meshType
character*50 filename, filename2, filename3
logical flag
double precision volumeFaces(3), totalArea, length
integer elem, localElements(3), localElem, numxLocal, numyLocal, numzLocal, maxNumNeighbors, globalMaxNeighbors
double precision, allocatable :: globalMeshCoord(:,:)
double precision, allocatable :: globalLength(:), procCoordList(:,:)
integer, allocatable :: globalMeshConnect(:,:,:), globalMaterial(:), globalNumNeighbors(:,:)
integer globalCell, globalNeighbor, maxElement

interface

	subroutine createConnectLocal(globalMeshConnect, globalMeshCoord, localElem, procCoordList)
	double precision, allocatable :: globalMeshCoord(:,:)
	double precision, allocatable :: procCoordList(:,:)
	integer, allocatable :: globalMeshConnect(:,:,:)
	integer localElem
	end subroutine
	
	subroutine createConnectGlobalPeriodic(Connect, Coord, NumNeighbors, Length, numTotal, bdryCoord)
	double precision, allocatable :: Coord(:,:), Length(:)
	integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
	double precision bdryCoord(6)
	integer numTotal
	end subroutine
	
	subroutine createConnectGlobalFreeSurf(Connect, Coord, NumNeighbors, Length, numTotal, bdryCoord)
	double precision, allocatable :: Coord(:,:), Length(:)
	integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
	double precision bdryCoord(6)
	integer numTotal
	end subroutine
	
end interface

open(80, file='TestInputNonUniform.txt',action='read', status='old')
flag=.FALSE.

!Step -1: read in mesh type
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='meshType') then
		read(80,*) meshType
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!Step 1: read in the global coordinates (max, min)
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='xminmax') then
		read(80,*) myProc%globalCoord(1), myProc%globalCoord(2)
	else if(readIn=='yminmax') then
		read(80,*) myProc%globalCoord(3), myProc%globalCoord(4)
	else if(readIn=='zminmax') then
		read(80,*) myProc%globalCoord(5), myProc%globalCoord(6)
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!if(myProc%taskid==MASTER) then
!	write(*,*) 'global coordinates'
!	write(*,*) myProc%globalCoord
!endif

!Step 2: divide processors over gobal volume
!NOTE: step to is identical to step 2 in the uniform case. Extended comments can be found in the uniform case.

!Step 2a: find area of each face of entire volume (x-y, y-z, x-z)
volumeFaces(1)=(myProc%globalCoord(2)-myProc%globalCoord(1))*&
	(myProc%globalCoord(4)-myProc%globalCoord(3))
volumeFaces(2)=(myProc%globalCoord(4)-myProc%globalCoord(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5))
volumeFaces(3)=(myProc%globalCoord(2)-myProc%globalCoord(1))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5))
	
!Step 2b: for each factorization of the number of processors, find the total shared area and minimize
totalArea=dble(myProc%numtasks)*(volumeFaces(1)+volumeFaces(2)+volumeFaces(3))
procDivision(1)=0
procDivision(2)=0
procDivision(3)=0

do i=1,myProc%numtasks
	do j=1, myProc%numtasks
		do k=1, myProc%numtasks
			if(i*j*k==myProc%numtasks) then
				!find area using this division
				if(k*volumeFaces(1)+i*volumeFaces(2)+j*volumeFaces(3) < totalArea) then
					procDivision(1)=i
					procDivision(2)=j
					procDivision(3)=k
					totalArea=k*volumeFaces(1)+i*volumeFaces(2)+j*volumeFaces(3)
					!write(*,*) 'new division found', procDivision, totalArea
				end if
			endif
		end do
	end do
end do

!if(myProc%taskid==MASTER) then
!	write(*,*) 'proc division', procDivision
!endif

!step 2c: divide volume among processors
myProc%localCoord(1)=myProc%globalCoord(1)+&
	mod(myProc%taskid,procDivision(1))*(myProc%globalCoord(2)-myProc%globalCoord(1))/&
	dble(procDivision(1))

myProc%localCoord(2)=myProc%localCoord(1)+&
	(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(procDivision(1))

myProc%localCoord(3)=myProc%globalCoord(3)+&
	mod(myProc%taskid/ProcDivision(1),procDivision(2))*&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))

myProc%localCoord(4)=myProc%localCoord(3)+&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))

myProc%localCoord(5)=myProc%GlobalCoord(5)+&
	mod(myProc%taskid/(ProcDivision(1)*ProcDivision(2)),procDivision(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))

myProc%localCoord(6)=myProc%localCoord(5)+&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))
	
totalVolume=(myProc%localCoord(2)-myProc%localCoord(1))*(myProc%localCoord(4)-myProc%localCoord(3))*&
	(myProc%localCoord(6)-myProc%localCoord(5))

!This contains a list of ALL processor coordinate max/mins, needed later
allocate(procCoordList(myProc%numtasks+1,6))
call createProcCoordList(procCoordList, procDivision)

!write(*,*) 'proc', myProc%taskid, 'of', myProc%numtasks
!write(*,*) 'local coordinates', myProc%localCoord

!Step 2d: point processor myProc%procNeighbor(x) to neighboring processors
if (myProc%localCoord(1)==myProc%globalCoord(1)) then	!coordinate is at xmin
	myProc%procNeighbor(2)=myProc%taskid+procDivision(1)-1
else
	myProc%procNeighbor(2)=myProc%taskid-1
end if

if (myProc%localCoord(2)==myProc%globalCoord(2)) then	!coordinate is at xmax
	myProc%procNeighbor(1)=myProc%taskid-procDivision(1)+1
else
	myProc%procNeighbor(1)=myProc%taskid+1
end if

if (myProc%localCoord(3)==myProc%globalCoord(3)) then	!coordinate is at ymin
	myProc%procNeighbor(4)=myProc%taskid+procDivision(1)*(procDivision(2)-1)
else
	myProc%procNeighbor(4)=myProc%taskid-procDivision(1)
end if

if (myProc%localCoord(4)==myProc%globalCoord(4)) then	!coordinate is at ymax
	myProc%procNeighbor(3)=myProc%taskid-procDivision(1)*(procDivision(2)-1)
else
	myProc%procNeighbor(3)=myProc%taskid+procDivision(1)
end if

if (myProc%localCoord(5)==myProc%globalCoord(5)) then	!coordinate is at zmin
	myProc%procNeighbor(6)=myProc%taskid+procDivision(1)*procDivision(2)*(procDivision(3)-1)
else
	myProc%procNeighbor(6)=myProc%taskid-procDivision(1)*procDivision(2)
end if

if (myProc%localCoord(6)==myProc%globalCoord(6)) then	!coordinate is at zmax
	myProc%procNeighbor(5)=myProc%taskid-procDivision(1)*procDivision(2)*(procDivision(3)-1)
else
	myProc%procNeighbor(5)=myProc%taskid+procDivision(1)*procDivision(2)
end if

!write(*,*) 'proc', myProc%taskid, 'of', myProc%numtasks
!write(*,*) 'neighbors', myProc%procNeighbor

!Step 3: create mesh in each processor of elements only within that processor's local coordinates

!Step 3a: read file info into global meshes

!Read total number of elements (global)
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='elements') then
		read(80,*) numTotal
		flag=.TRUE.
	end if
end do
flag=.FALSE.

allocate(globalMeshCoord(numTotal,3))		!coordinates of the center of each element
allocate(globalMaterial(numTotal))			!material type of each element
allocate(globalNumNeighbors(numTotal,6))	!number of neighbors of each element in each direction
											!(directions: +x, -x, +y, -y, +z, -z)
allocate(globalLength(numTotal))			!volume element length

!Count the number of elements that are inside the local processor's bounds
localElem=0
globalMaxNeighbors=0
do elem=1,numTotal
	!read in element coordinates and material number as well as how many neighbors it has in each direction
	read(80,*) (globalMeshCoord(elem,i),i=1,3), globalLength(elem), globalMaterial(elem), (globalNumNeighbors(elem,i), i=1,6)
	
	!GlobalMaxNeighbors is used to create the global connectivity matrix. It is the max number of 
	!neighbors in one direction that any one element has (if elements are of different sizes, one element
	!may have more than one neighbor in a given direction).
	
	do i=1,6
		if(globalNumNeighbors(elem,i) > globalMaxNeighbors) then
			globalMaxNeighbors=globalNumNeighbors(elem,i)
		end if
	end do
	
	!Count the number of elements that are inside the local coordinates.
	if(globalMeshCoord(elem,1) > myProc%localCoord(1) .AND. globalMeshCoord(elem,1) <= myProc%localCoord(2) &
		.AND. globalMeshCoord(elem,2) > myProc%localCoord(3) .AND. globalMeshCoord(elem,2) <= myProc%localCoord(4) &
		.AND. globalMeshCoord(elem,3) > myProc%localCoord(5) .AND. globalMeshCoord(elem,3) <= myProc%localCoord(6)) then
		
		localElem=localElem+1
	end if
end do

!If there are processors with no local volume elements (too many procs), error.
if(localElem==0) then
	write(*,*) 'error processors with no volume elements'
	call MPI_ABORT(MPI_COMM_WORLD,ierr)
	stop
end if

!write(*,*) localElem, 'proc', myProc%taskid

!Step 3b: create local mesh of coordinates
allocate(myMesh(localElem))
!set the global variable numCells, which is the sie of the array myMesh (used at other points in SRSCD code)
numCells=localElem

!Here, we read in the local element coordinates from the global element matrix. We also read in the processor
!number and the number of neighbors. We then allocate the neighbors matrix with the max number of neighbors
!that a given element has (it may have one element neighbor in one direction and multiple neighbors
!in a different direction).

localElem=0
do elem=1,numTotal
	if(globalMeshCoord(elem,1) > myProc%localCoord(1) .AND. globalMeshCoord(elem,1) <= myProc%localCoord(2) &
		.AND. globalMeshCoord(elem,2) > myProc%localCoord(3) .AND. globalMeshCoord(elem,2) <= myProc%localCoord(4) &
		.AND. globalMeshCoord(elem,3) > myProc%localCoord(5) .AND. globalMeshCoord(elem,3) <= myProc%localCoord(6)) then
		
		localElem=localElem+1	!count local elements
		
		myMesh(localElem)%coordinates(1)=globalMeshCoord(elem,1)	!read in coordinates, material, length, etc.
		myMesh(localElem)%coordinates(2)=globalMeshCoord(elem,2)
		myMesh(localElem)%coordinates(3)=globalMeshCoord(elem,3)
		myMesh(localElem)%material=globalMaterial(elem)
		myMesh(localElem)%length=globalLength(elem)
		myMesh(localElem)%volume=globalLength(elem)**3d0
		maxNumNeighbors=0
		!Find the max number of neighbors that this element has in any one direction
		do i=1,6
			myMesh(localElem)%numNeighbors(i)=globalNumNeighbors(elem,i)
			if(myMesh(localElem)%numNeighbors(i) > maxNumNeighbors) then
				maxNumNeighbors=myMesh(localElem)%numNeighbors(i)
			end if
		end do
		myMesh(localElem)%proc=myProc%taskid
		
		!The number of neighbors allowed in all directions is equal to maxNumNeighbors even if
		!the element may have fewer neighbors in some directions. These extra array elements are 
		!ignored for the rest of the program.
		
		allocate(myMesh(localElem)%neighbors(6,maxNumNeighbors))
		allocate(myMesh(localElem)%neighborProcs(6,maxNumNeighbors))
		
	end if
end do

!Step 4: create global and local connectivity matrixes including multiple neighbors to the same direction

!Step 4a: create global connectivity matrix
allocate(globalMeshConnect(numTotal,6,globalMaxNeighbors))

if(meshType=='periodic') then
	call createConnectGlobalPeriodicNonUniform(globalMeshConnect, globalMeshCoord, globalNumNeighbors, globalLength, &
		numTotal, myProc%globalCoord)
else if(meshType=='freeSurfaces') then
	call createConnectGlobalFreeSurfNonUniform(globalMeshConnect, globalMeshCoord, globalNumNeighbors, globalLength, &
		numTotal, myProc%globalCoord)
end if

!if(myProc%taskid==MASTER) then
!	write(*,*) 'Global elements and connectivity'
!	do 19 i=1,numTotal
!		write(*,*) 'element', i, 'length', globalLength(i), 'material', globalMaterial(i)
!		write(*,*) 'coords', (globalMeshCoord(i,j),j=1,3)
!		write(*,*) 'numNeighbors', (globalNumNeighbors(i,j),j=1,6)
!		write(*,*) 'neighbors'
!		do 20 j=1,6
!			write(*,*) (globalMeshConnect(i,j,k),k=1,globalNumNeighbors(i,j))
!		20 continue
!	19 continue
!endif

!Step 4b: create local connectivity myMesh(:)%neighbors(direction,num) and myMesh(:)%neighborProcs(direction,num)

call createConnectLocalNonUniform(globalMeshConnect, globalMeshCoord, localElem, procCoordList)

!OPTIONAL: output to check that the mesh is created correctly

!if(myProc%taskid==MASTER) then
!	write(*,*) 'Local elements and connectivity'
!	do 21 i=1,localElem
!		write(*,*) 'element', i, 'length', myMesh(i)%length, 'material', myMesh(i)%material
!		write(*,*) 'coords', (myMesh(i)%coordinates(j),j=1,3)
!		write(*,*) 'numNeighbors', (myMesh(i)%numNeighbors(j),j=1,6)
!		write(*,*) 'neighbors and procs'
!		do 22 j=1,6
!			write(*,*) (myMesh(i)%neighbors(j,k),k=1,myMesh(i)%numNeighbors(j))
!			write(*,*) (myMesh(i)%neighborProcs(j,k),k=1,myMesh(i)%numNeighbors(j))
!			write(*,*)
!		22 continue
!	21 continue
!	
!	write(*,*) 'connectivity order +x, -x, +y, -y, +z, -z'
!	write(*,*) 'proc division', procDivision, 'Division order: x, y, z'
!	write(*,*)
!	write(*,*) 'Global coordinates (xmin, xmax, ymin, ymax, zmin, zmax)'
!	write(*,*) myProc%globalCoord
!	write(*,*)
!endif

!***************************************************************************************************
!5/28/2014: initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************
!Using the following global mesh info:
!
!double precision, allocatable :: globalMeshCoord(:,:), globalLength(:), procCoordList(:,:)
!integer, allocatable :: globalMeshConnect(:,:,:), globalMaterial(:), globalNumNeighbors(:,:)

!Step 1: Find the max cell# of any boundary mesh element
maxElement=0
do i=1,localElem
	do j=1,6
		do k=1,myMesh(i)%numNeighbors(j)
			if(myMesh(i)%neighborProcs(j,k) /= myProc%taskid) then	!we are pointing to a different proc
				if(myMesh(i)%neighbors(j,k) > maxElement) then		!searching for the max element number in a neighbor
					maxElement=myMesh(i)%neighbors(j,k)
				endif
			endif
		end do
	end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(6,maxElement))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do i=1,maxElement
	do j=1,6
		myBoundary(j,i)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
	end do
end do

!Step 2: initialize myBoundary elements (only the relevant ones) with processor #'s, material #'s, length
do i=1,localElem
	do j=1,6
		do k=1,myMesh(i)%numNeighbors(j)
			if(myMesh(i)%neighborProcs(j,k) /= myProc%taskid) then
				myBoundary(j,myMesh(i)%neighbors(j,k))%proc=myMesh(i)%neighborProcs(j,k)	!set proc # of elements in myBoundary
				globalCell=findGlobalCell(myMesh(i)%coordinates,globalMeshCoord)			!find global cell # of element in myBoundary
				globalNeighbor=globalMeshConnect(globalCell,j,k)								!use global cell # to find material # of element in myBoundary
				myBoundary(j,myMesh(i)%neighbors(j,1))%material=globalMaterial(globalNeighbor)	!set material # of elements in myBoundary
				myBoundary(j,myMesh(i)%neighbors(j,k))%length=globalLength(globalNeighbor)		!set length of elements in myBoundary
				myBoundary(j,myMesh(i)%neighbors(j,k))%volume=globalLength(globalNeighbor)**3d0	!set volume of elements in myBoundary (changes with cascade addition)
				myBoundary(j,myMesh(i)%neighbors(j,1))%localNeighbor=i
			endif
		end do
	end do
end do

close(80)

end subroutine

!This subroutine creates the global connectivity matrix (non-uniform)

!>Subroutine create global connectivity (non-uniform mesh, periodic boundary conditions)
!!
!!This subroutine creates a global mesh (not dividing between processors) as well as the 
!!connectivity between elements for the case of a non-uniform periodic cubic mesh.
!!
!!Inputs: coordinates list, length list, total number of elements, coordinates of boundary
!!Output: connectivity, number of neighbors for each element

subroutine createConnectGlobalPeriodicNonUniform(Connect, Coord, NumNeighbors, Length, numTotal, bdryCoord)
implicit none

double precision, allocatable :: Coord(:,:), Length(:)
integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
integer numTotal, cell, dir, i, searchElem, neighbor
double precision neighborElem, bdryCoord(6)

!********************************************************
!periodic boundary condition version
!********************************************************

do cell=1,numTotal
	do dir=1,6
		neighbor=0
		do searchElem=1,numTotal
			if(dir==1) then
				
				!Here we are using the coordinates of the cell to determine if it is at the edge of the boundary
				!instead of the cell number (in the uniform cubic mesh case).
				
				if(Coord(cell,1)+length(cell)/2d0==bdryCoord(2)) then	
					!periodic
					!We must search through all elements to see if they are the neighbors of other elements. This is because we 
					!don't know how many neighbors each element has. This could be made more efficient in the future if it becomes
					!too computationally expensive, but for now it only has to be carried out once so it should be fine.
					
					if(Coord(searchElem,1)==Coord(cell,1)+(length(cell)+length(searchElem))/2d0-(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1						!count the number of neighbors
						Connect(cell,dir,neighbor)=searchElem	!Add this cell to the connectivity
					end if
				else if(Coord(searchElem,1)==Coord(cell,1)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==2) then
				if(Coord(cell,1)-length(cell)/2d0==bdryCoord(1)) then
					!periodic
					if(Coord(searchElem,1)==Coord(cell,1)-(length(cell)+length(searchElem))/2d0+(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,1)==Coord(cell,1)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==3) then
				if(Coord(cell,2)+length(cell)/2d0==bdryCoord(4)) then
					!periodic
					if(Coord(searchElem,2)==Coord(cell,2)+(length(cell)+length(searchElem))/2d0-(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,2)==Coord(cell,2)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==4) then
				if(Coord(cell,2)-length(cell)/2d0==bdryCoord(3)) then
					!periodic
					if(Coord(searchElem,2)==Coord(cell,2)-(length(cell)+length(searchElem))/2d0+(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,2)==Coord(cell,2)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==5) then
				if(Coord(cell,3)+length(cell)/2d0==bdryCoord(6)) then
					!periodic
					if(Coord(searchElem,3)==Coord(cell,3)+(length(cell)+length(searchElem))/2d0-(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,3)==Coord(cell,3)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==6) then
				if(Coord(cell,3)-length(cell)/2d0==bdryCoord(5)) then
					!periodic
					if(Coord(searchElem,3)==Coord(cell,3)-(length(cell)+length(searchElem))/2d0+(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,3)==Coord(cell,3)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			end if
		end do
		if(NumNeighbors(cell,dir) /= neighbor) then
			!we have found the wrong number of neighbors
			write(*,*) 'error incorrect number of neighbors found', neighbor, 'neighbors', numNeighbors(cell,dir), 'expected'
			write(*,*) 'cell', cell
		end if
	end do
end do

end subroutine


!This is the same as above except that the periodic aspect is removed and cells at boundaries are given '0' neighbors (in +/- z-dir)
!>Subroutine create global connectivity (non-uniform mesh, free surfaces in the z-direction and periodic boundary conditions in other directions)
!!
!!This subroutine creates a global mesh (not dividing between processors) as well as the 
!!connectivity between elements for the case of a non-uniform cubic mesh with free surfaces.
!!
!!Inputs: coordinates list, length list, total number of elements, coordinates of boundary
!!Output: connectivity, number of neighbors for each element

subroutine createConnectGlobalFreeSurfNonUniform(Connect, Coord, NumNeighbors, Length, numTotal, bdryCoord)
implicit none

double precision, allocatable :: Coord(:,:), Length(:)
integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
integer numTotal, cell, dir, i, searchElem, neighbor
double precision neighborElem, bdryCoord(6)

!********************************************************
!free surfaces in +/- z-dir boundary condition version
!********************************************************

do cell=1,numTotal
	do dir=1,6
		neighbor=0
		do searchElem=1,numTotal
			if(dir==1) then
				if(Coord(cell,1)+length(cell)/2d0==bdryCoord(2)) then
					!periodic
					if(Coord(searchElem,1)==Coord(cell,1)+(length(cell)+length(searchElem))/2d0-(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,1)==Coord(cell,1)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==2) then
				if(Coord(cell,1)-length(cell)/2d0==bdryCoord(1)) then
					!periodic
					if(Coord(searchElem,1)==Coord(cell,1)-(length(cell)+length(searchElem))/2d0+(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,1)==Coord(cell,1)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==3) then
				if(Coord(cell,2)+length(cell)/2d0==bdryCoord(4)) then
					!periodic
					if(Coord(searchElem,2)==Coord(cell,2)+(length(cell)+length(searchElem))/2d0-(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,2)==Coord(cell,2)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==4) then
				if(Coord(cell,2)-length(cell)/2d0==bdryCoord(3)) then
					!periodic
					if(Coord(searchElem,2)==Coord(cell,2)-(length(cell)+length(searchElem))/2d0+(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=searchElem
					end if
				else if(Coord(searchElem,2)==Coord(cell,2)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,3)-Coord(cell,3)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==5) then
				if(Coord(cell,3)+length(cell)/2d0==bdryCoord(6)) then
					!periodic
					if(Coord(searchElem,3)==Coord(cell,3)+(length(cell)+length(searchElem))/2d0-(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=0
					end if
				else if(Coord(searchElem,3)==Coord(cell,3)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			else if(dir==6) then
				if(Coord(cell,3)-length(cell)/2d0==bdryCoord(5)) then
					!periodic
					if(Coord(searchElem,3)==Coord(cell,3)-(length(cell)+length(searchElem))/2d0+(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
						
						neighbor=neighbor+1
						Connect(cell,dir,neighbor)=0
					end if
				else if(Coord(searchElem,3)==Coord(cell,3)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,2)-Coord(cell,2)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(searchElem,1)-Coord(cell,1)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(cell,dir,neighbor)=searchElem
				end if
			end if
		end do
		if(NumNeighbors(cell,dir) /= neighbor) then
			!we have found the wrong number of neighbors
			write(*,*) 'error incorrect number of neighbors found', neighbor, 'neighbors', numNeighbors(cell,dir), 'expected'
			write(*,*) 'cell', cell
		end if
	end do
end do

end subroutine

!This subroutine creates a local mesh, including neighbors (connectivity) and neighboring procs

!>Subroutine create local connectivity (non-uniform mesh, either free surfaces or periodic bc's)
!!
!!This subroutine creates a local mesh (for the local processor only) as well as the 
!!connectivity between elements for the case of a non-uniform cubic mesh. This subroutine
!!has to choose which elements of the global mesh are in the local mesh, and puts those in
!!myMesh along with the element coordinates and their neighbor element ID numbers and processors.
!!
!!Inputs: global mesh coordinates, global mesh connectivity, list of local elements, list of coordinate boundaries of each processor
!!Output: myMesh

subroutine createConnectLocalNonUniform(globalMeshConnect, globalMeshCoord, localElem, procCoordList)
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'
double precision, allocatable :: globalMeshCoord(:,:)
double precision, allocatable :: procCoordList(:,:)
integer, allocatable :: globalMeshConnect(:,:,:)
integer localElem, elem, neighbor, dir, globalCell, localNeighbor, globalNeighbor, numSendRecv
integer sendList(6*localElem, 5), i, status(MPI_STATUS_SIZE), neighborProc

!For each local element, find its global cell, find its neighbors (global) in the globalNeighbor variable,
!and if globalNeighbor is not in the local mesh, then send/recieve with the neighboring processor
!in the buffer system (same as uniform mesh)

numSendRecv=0
do elem=1,localElem
	globalCell=findGlobalCell(myMesh(elem)%coordinates,globalMeshCoord)
	do dir=1,6
		do neighbor=1,myMesh(elem)%numNeighbors(dir)
			globalNeighbor=globalMeshConnect(globalCell,dir,neighbor)
			if(globalNeighbor==0) then
				!free surface, set cell number to 0 and proc number to -1
				myMesh(elem)%neighbors(dir,neighbor)=0
				myMesh(elem)%neighborProcs(dir,neighbor)=-1
			else
				!find out if globalNeighbor is in the local mesh. If so, set myMesh(cell)%neighbors and procs
				if(globalMeshCoord(globalNeighbor,1) > myProc%localCoord(1) &
				.AND. globalMeshCoord(globalNeighbor,1) <= myProc%localCoord(2) &
				.AND. globalMeshCoord(globalNeighbor,2) > myProc%localCoord(3) &
				.AND. globalMeshCoord(globalNeighbor,2) <= myProc%localCoord(4) &
				.AND. globalMeshCoord(globalNeighbor,3) > myProc%localCoord(5) &
				.AND. globalMeshCoord(globalNeighbor,3) <= myProc%localCoord(6)) then
					
					localNeighbor=findLocalCell(globalMeshCoord(globalNeighbor,:))	!find cell number of neighbor in local mesh
					myMesh(elem)%neighbors(dir,neighbor)=localNeighbor				
					myMesh(elem)%neighborProcs(dir,neighbor)=myProc%taskid
				!if no, send local mesh number and recieve local mesh number from neigboring processor
				else
					!find the number of the neighboring processor based on the coordinates of the neighboring cell
					neighborProc=findNeighborProc(globalMeshCoord,procCoordList,globalNeighbor)
					myMesh(elem)%neighborProcs(dir,neighbor)=neighborProc
						
					numSendRecv=numSendRecv+1
					
					!Create buffer of information to send/recieve between processors
					sendList(numSendRecv,1)=elem
					sendList(numSendRecv,2)=dir	!direction
					sendList(numSendRecv,3)=globalCell
					sendList(numSendRecv,4)=globalNeighbor
					sendList(numSendRecv,5)=neighbor	!which neighbor number is this
				end if
			end if
		end do
	end do
end do

!Send/recieve the element in the local mesh / the neighbor in the neighboring processor's mesh (respectively).
!MPI_SEND and MPI_RECV functions are matched using the cell number of the sender

do i=1,numSendRecv
	
	elem=sendList(i,1)
	dir=sendList(i,2)
	globalCell=sendList(i,3)
	neighbor=sendList(i,5)
	call MPI_SEND(elem, 1, MPI_INTEGER, myMesh(elem)%neighborProcs(dir,neighbor), globalCell, MPI_COMM_WORLD, ierr)
end do
	
do i=1,numSendRecv
	elem=sendList(i,1)
	dir=sendList(i,2)
	globalNeighbor=sendList(i,4)
	neighbor=sendList(i,5)
	call MPI_RECV(myMesh(elem)%neighbors(dir,neighbor), 1, MPI_INTEGER, myMesh(elem)%neighborProcs(dir,neighbor), &
				globalNeighbor, MPI_COMM_WORLD, status, ierr)
end do

end subroutine			

!Find the cell number in the local mesh using coordinates
!>Function find local cell
!!
!!Finds the cell ID number in the local mesh using the coordinates of the centroid of the cell
!!
!!Input: coordinates
!!Output: cell ID number

integer function findLocalCell(coord)
use DerivedType
use mod_srscd_constants
implicit none

double precision coord(3)
integer i
logical flag

flag = .FALSE.
i=0
do 10 while(flag .eqv. .FALSE.)
	i=i+1
	if(coord(1)==myMesh(i)%coordinates(1) .AND. coord(2)==myMesh(i)%coordinates(2) .AND. coord(3)==myMesh(i)%coordinates(3)) then
		flag=.TRUE.
	endif
10 continue

findLocalCell=i
end function

!create a list of the bounds of each processor (not just the local processor), to be used
!to identify which processor a given cell is inside (necessary for nonuniform meshes)

!>subroutine create processor coordinates list
!!
!!This subroutine creates a list of the bounds (coordinates) of each processor to be used
!!to identify which processor contains a given volume element (necessary for non-uniform
!!mesh creation)
!!
!!Inputs: division of processors (see main subroutine)
!!Output: list of boundaries of each processor

subroutine createProcCoordList(procCoordList, procDivision)
use DerivedType
use mod_srscd_constants
implicit none

integer procDivision(3), i
double precision, allocatable :: procCoordList(:,:)

do 10 i=1,myProc%numtasks
	procCoordList(i,1)=myProc%globalCoord(1)+&
		mod(i-1,procDivision(1))*(myProc%globalCoord(2)-myProc%globalCoord(1))/&
		dble(procDivision(1))
	
	procCoordList(i,2)=procCoordList(i,1)+&
		(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(procDivision(1))
	
	procCoordList(i,3)=myProc%globalCoord(3)+&
		mod((i-1)/ProcDivision(1),procDivision(2))*&
		(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))
	
	procCoordList(i,4)=procCoordList(i,3)+&
		(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(procDivision(2))
	
	procCoordList(i,5)=myProc%GlobalCoord(5)+&
		mod((i-1)/(ProcDivision(1)*ProcDivision(2)),procDivision(3))*&
		(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))
	
	procCoordList(i,6)=procCoordList(i,5)+&
		(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(procDivision(3))
10 continue
end subroutine

!Find the number of the neighboring processor given coordinates within that processor

!>Function find neighboring processor
!! 
!!This function finds the processor ID for a volume element that has ID number==elem in the 
!!GLOBAL connectivity.
!!
!!Inputs: list of global mesh coordinates, list of bounds of each processor's domain, element number
!!Output: processor ID for a given volume element

integer function findNeighborProc(globalMeshCoord,procCoordList,elem)
use DerivedType
use mod_srscd_constants
implicit none

double precision, allocatable :: globalMeshCoord(:,:)
double precision, allocatable :: procCoordList(:,:)
integer elem, i

do 10 i=1,myProc%numtasks
	if(globalMeshCoord(elem,1) .GT. procCoordList(i,1) .AND. globalMeshCoord(elem,1) .LE. procCoordList(i,2) &
	.AND. globalMeshCoord(elem,2) .GT. procCoordList(i,3) .AND. globalMeshCoord(elem,2) .LE. procCoordList(i,4) &
	.AND. globalMeshCoord(elem,3) .GT. procCoordList(i,5) .AND. globalMeshCoord(elem,3) .LE. procCoordList(i,6)) then
		findNeighborProc=i-1
		exit
	endif
10 continue

end function


end module

