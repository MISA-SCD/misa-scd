!**************************************************************************************************
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
!**************************************************************************************************

module MeshReader
use mod_constants
use DerivedType
implicit none

contains

!***************************************************************************************************
!>Subroutine read Mesh Uniform - creates processor and mesh files from uniform mesh input
!
!This is the main subroutine that controls reading the uniform mesh. It reads in the mesh
!from a file and divides the volume elements between the processors in the parallel
!simulation. If a division such that each processor has at least one element is not possible,
!then the subroutine returns an error. It also creates a 'processor mesh', indicating
!the bounds of each processor. NOTE: this subroutine is for uniform meshes only, a different
!subroutine has been created to read in non-uniform meshes (and these have not been
!implemented in the rest of the program)
!
!Input: filename of mesh file
!
!Outputs: myProc (processors with meshes), myMesh, and myBoundary
!***************************************************************************************************

subroutine readMeshUniform(filename)
use mod_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), procDivision(3), i, j, k, l, dir, matNum
character*20 readIn, meshType
character*50 filename, filename2, filename3
logical flag
double precision length, tempCoord(3)
integer element,  localElem
integer tempx1,tempx2,tempy1,tempy2,tempz1,tempz2	!used to determine numxLocal, numyLocal, numzLocal

double precision tempMeshCoord(3), tempStrain(6)
integer tempMaterial

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

!************************************************
!Step 1a: modify global coordinates to include entire simulation volume (not just centers of elements)
!************************************************

!The actual boundary coordinates of the system
do i=1,5,2
	myProc%globalCoord(i)=myProc%globalCoord(i)-length/2d0
	myProc%globalCoord(i+1)=myProc%globalCoord(i+1)+length/2d0
end do

!step 3: divide volume among processors according to the factorization above
!This step locates the boundaries of the local processor within the global volume.
!**********
myProc%localCoord(1)=myProc%globalCoord(1)+myProc%coords(1)*&
		(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(dims(1))	!xmin

myProc%localCoord(2)=myProc%localCoord(1)+&
	(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(dims(1))	!xmax

myProc%localCoord(3)=myProc%globalCoord(3)+myProc%coords(2)*&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))	!ymin

myProc%localCoord(4)=myProc%localCoord(3)+&
	(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))	!ymax

myProc%localCoord(5)=myProc%globalCoord(5)+myProc%coords(3)*&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))	!zmin

myProc%localCoord(6)=myProc%localCoord(5)+&
	(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))	!zmax
	
totalVolume=(myProc%localCoord(2)-myProc%localCoord(1))*(myProc%localCoord(4)-myProc%localCoord(3))*&
	(myProc%localCoord(6)-myProc%localCoord(5))

!Step 4: create mesh in each processor of elements only within that processor's local coordinates
!NOTE: this is an actual mesh of volume elements, not the 'processor mesh' created above.

!Step 4a: find how many volume elements are in global system and which of those are in the local processor's range

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

!Tells program that we are about to start reading in element coordinates and material numbers
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='elements') then
		flag=.TRUE.
	endif
end do
flag=.FALSE.

!Step 4b: find how many volume elements are in local processor and allocate myMesh accordingly
!numxLocal, numyLocal, numzLocal are used to determine the size of the mesh inside the local processor.
numxLocal=0
numyLocal=0
numzLocal=0

!get numxLocal
if(myProc%localCoord(1)==myProc%globalCoord(1) .AND. myProc%localCoord(2)==myProc%globalCoord(2)) then
	numxLocal=numx
else if(myProc%localCoord(1)==myProc%globalCoord(1)) then	!at xmin
	tempx1=0
	tempx2 = myProc%localCoord(2)/(length/2d0)
	if((dble(tempx2)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2)< (dble(tempx2+1)*(length/2d0))) then
		if(mod(tempx2,2)==0) then
			tempx2 = tempx2/2
		else
			tempx2 = (tempx2+1)/2
		end if
	end if
	numxLocal=tempx2
else if(myProc%localCoord(2)==myProc%globalCoord(2)) then	!at xmax
	tempx2=0
	tempx1 = myProc%localCoord(1)/(length/2d0)
	if((dble(tempx1)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(tempx1+1)*(length/2d0))) then
		if(mod(tempx1,2)==0) then
			tempx1 = tempx1/2
		else
			tempx1 = (tempx1+1)/2
		end if
	end if
	numxLocal=numx-tempx1

else	!in the middle
	tempx1 = myProc%localCoord(1)/(length/2d0)
	tempx2 = myProc%localCoord(2)/(length/2d0)
	if((dble(tempx1)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(tempx1+1)*(length/2d0))) then
		if(mod(tempx1,2)==0) then
			tempx1 = tempx1/2
		else
			tempx1 = (tempx1+1)/2
		end if
	end if
	if((dble(tempx2)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2) < (dble(tempx2+1)*(length/2d0))) then
		if(mod(tempx2,2)==0) then
			tempx2 = tempx2/2
		else
			tempx2 = (tempx2+1)/2
		end if
	end if
	numxLocal = tempx2-tempx1

end if

!get numyLocal
if(myProc%localCoord(3)==myProc%globalCoord(3) .AND. myProc%localCoord(4)==myProc%globalCoord(4)) then
	numyLocal=numy
else if(myProc%localCoord(3)==myProc%globalCoord(3)) then	!at ymin
	tempy1=0
	tempy2 = myProc%localCoord(4)/(length/2d0)
	if((dble(tempy2)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(tempy2+1)*(length/2d0))) then
		if(mod(tempy2,2)==0) then
			tempy2 = tempy2/2
		else
			tempy2 = (tempy2+1)/2
		end if
	end if
	numyLocal=tempy2
else if(myProc%localCoord(4)==myProc%globalCoord(4)) then	!at ymax
	tempy2=0
	tempy1 = myProc%localCoord(3)/(length/2d0)
	if((dble(tempy1)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(tempy1+1)*(length/2d0))) then
		if(mod(tempy1,2)==0) then
			tempy1 = tempy1/2
		else
			tempy1 = (tempy1+1)/2
		end if
	end if
	numyLocal=numy-tempy1

else	!in the middle
	tempy1 = myProc%localCoord(3)/(length/2d0)
	tempy2 = myProc%localCoord(4)/(length/2d0)
	if((dble(tempy1)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(tempy1+1)*(length/2d0))) then
		if(mod(tempy1,2)==0) then
			tempy1 = tempy1/2
		else
			tempy1 = (tempy1+1)/2
		end if
	end if
	if((dble(tempy2)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(tempy2+1)*(length/2d0))) then
		if(mod(tempy2,2)==0) then
			tempy2 = tempy2/2
		else
			tempy2 = (tempy2+1)/2
		end if
	end if
	numyLocal = tempy2-tempy1

end if

!get numzLocal
if(myProc%localCoord(5)==myProc%globalCoord(5) .AND. myProc%localCoord(6)==myProc%globalCoord(6)) then
	numzLocal=numz
else if(myProc%localCoord(5)==myProc%globalCoord(5)) then	!at zmin
	tempz1=0
	tempz2 = myProc%localCoord(6)/(length/2d0)
	if((dble(tempz2)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(tempz2+1)*(length/2d0))) then
		if(mod(tempz2,2)==0) then
			tempz2 = tempz2/2
		else
			tempz2 = (tempz2+1)/2
		end if
	end if
	numzLocal=tempz2
else if(myProc%localCoord(6)==myProc%globalCoord(6)) then	!at zmax
	tempz2=0
	tempz1 = myProc%localCoord(5)/(length/2d0)
	if((dble(tempz1)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5)< (dble(tempz1+1)*(length/2d0))) then
		if(mod(tempz1,2)==0) then
			tempz1 = tempz1/2
		else
			tempz1 = (tempz1+1)/2
		end if
	end if
	numzLocal=numz-tempz1

else	!in the middle
	tempz1 = myProc%localCoord(5)/(length/2d0)
	tempz2 = myProc%localCoord(6)/(length/2d0)
	if((dble(tempz1)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5) < (dble(tempz1+1)*(length/2d0))) then
		if(mod(tempz1,2)==0) then
			tempz1 = tempz1/2
		else
			tempz1 = (tempz1+1)/2
		end if
	end if
	if((dble(tempz2)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(tempz2+1)*(length/2d0))) then
		if(mod(tempz2,2)==0) then
			tempz2 = tempz2/2
		else
			tempz2 = (tempz2+1)/2
		end if
	end if
	numzLocal = tempz2-tempz1

end if

!this variable is used at various points in SRSCD; lets us know the length fo myMesh(:)
numCells=numxLocal*numyLocal*numzLocal	!total cells of this processor

!If any processors don't have volume elements in them (too many procs), we create an error message
if(numCells==0) then
	write(*,*) 'error processors with no volume elements'
	call MPI_ABORT(comm,ierr)
	stop
end if

write(*,*) 'proc', myProc%taskid, 'numx numy numz', numxLocal, numyLocal, numzLocal

!Step 3c: Read meshes. And for each volume element in myMesh, assign coordinates and material number (and proc number)
!We can't read from the file directly into myMesh, because we don't know how big to make it until
!we read in all of the volume element information. Therefore, we allocate it and then reread from the
!global mesh back into myMesh.
allocate(myMesh(numCells))
element=0
localElem=0
systemVol=0d0

do k=1,numz
	do j=1,numy
		do i=1,numx
			!read in the coordinates and material number from file of this element (coordinates are the center of the element)
			read(80,*) tempMeshCoord(1),tempMeshCoord(2),tempMeshCoord(3),tempMaterial
			element=element+1

			if(tempMeshCoord(1) > myProc%localCoord(1) .AND. tempMeshCoord(1) <= myProc%localCoord(2) &
					.AND. tempMeshCoord(2) > myProc%localCoord(3) .AND. tempMeshCoord(2) <= myProc%localCoord(4) &
					.AND. tempMeshCoord(3) > myProc%localCoord(5) .AND. tempMeshCoord(3) <= myProc%localCoord(6)) then
				localElem=localElem+1
				myMesh(localElem)%coordinates(1)=tempMeshCoord(1)
				myMesh(localElem)%coordinates(2)=tempMeshCoord(2)
				myMesh(localElem)%coordinates(3)=tempMeshCoord(3)
				myMesh(localElem)%material=tempMaterial
				myMesh(localElem)%proc=myProc%taskid
				myMesh(localElem)%globalCell=element

				myMesh(localElem)%length=length
				myMesh(localElem)%volume=length**3d0

				!uniform mesh: all elements have 1 neighbor in each direction
				allocate(myMesh(localElem)%neighbors(6,1))
				allocate(myMesh(localElem)%neighborProcs(6,1))
				do dir=1,6
					myMesh(localElem)%numNeighbors(dir)=1
				end do
			end if

			if(strainField=='yes') then
			
				read(50,*) tempCoord(1), tempCoord(2), tempCoord(3), tempStrain(1), tempStrain(2), &
						tempStrain(3), tempStrain(4), tempStrain(5), tempStrain(6)
					
				if(tempCoord(1) /= tempMeshCoord(1)) then	!not equal
					write(*,*) 'Error strain mesh does not match normal mesh'
				else if(tempCoord(2) /= tempMeshCoord(2)) then
					write(*,*) 'Error strain mesh does not match normal mesh'
				else if(tempCoord(3) /= tempMeshCoord(3)) then
					write(*,*) 'Error strain mesh does not match normal mesh'
				end if
				if(tempCoord(1) > myProc%localCoord(1) .AND. tempCoord(1) <= myProc%localCoord(2) &
						.AND. tempCoord(2) > myProc%localCoord(3) .AND. tempCoord(2) <= myProc%localCoord(4) &
						.AND. tempCoord(3) > myProc%localCoord(5) .AND. tempCoord(3) <= myProc%localCoord(6)) then
					do l=1,6
						myMesh(localElem)%strain(l)=tempStrain(l)
					end do
				end if
			end if
				
			if(numMaterials==1) then
				systemVol=systemVol+length**3d0
			else if(tempMaterial==1) then	!only add to system volume if we are NOT at a grain boundary
				systemVol=systemVol+length**3d0
			else
				!Do nothing
			end if

		end do
	end do
end do

!Step 3d: assign neighbors and processor numbers for neighbors (connectivity in myMesh) - periodic or free surfaces in +/- z
if(meshType=='periodic') then
	call createConnectLocalPeriodicUniform(length)
else if(meshType=='freeSurfaces') then
	call createConnectLocalFreeSurfUniform(length)
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
!endif

!Close input files
if(strainField=='yes') then
	close(50)
end if
close(80)

end subroutine

!***************************************************************************************************
!>Subroutine initial Mesh Uniform - creates processor and uniform meshes
!
!This is the main subroutine that controls creating the uniform mesh. It creates the mesh
!and divides the volume elements between the processors in the parallel simulation.
!If a division such that each processor has at least one element is not possible,
!then the subroutine returns an error. It also creates a 'processor mesh', indicating
!the bounds of each processor. NOTE: this subroutine is for uniform meshes only, a different
!subroutine has been created to read in non-uniform meshes (and these have not been
!implemented in the rest of the program)
!
!Input: MeshGenInput_XX.txt
!
!Outputs: myProc (processors with meshes), myMesh, and myBoundary
!***************************************************************************************************

subroutine initialMeshUniform(filename)
use mod_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), i, j, k, l, dir, matNum
character*20 readIn, meshType
character*50 filename, filename2, filename3
logical flag
double precision length, tempCoord(3)
integer element,  localElem, grain
integer x, y, z, tempx1,tempx2,tempy1,tempy2,tempz1,tempz2	!used to determine numxLocal, numyLocal, numzLocal

double precision tempMeshCoord(3), tempStrain(6), tempCenter(6)
integer tempMaterial

open(80, file=filename,action='read', status='old')	!MeshGenInput_XX.txt

!if(strainField=='yes') then
!	open(50, file=strainFileName, action='read', status='old')
!end if

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

!Step 0-1: Read the number of elements in x,y,and z directions
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
systemVol = dble(numTotal)*length**3d0

!Tells program that we are about to start reading in element coordinates and material numbers
do while(flag .eqv. .FALSE.)
	read(80,*) readIn
	if(readIn=='end') then
		flag=.TRUE.
	endif
end do
flag=.FALSE.

!***********************************************************************
!Read in the strain field into myMesh%strain(1-6)
!***********************************************************************

!if(strainField=='yes') then

	!Step -1: skip to start of input data
!	do while(flag .eqv. .FALSE.)
!		read(50,*) readIn
!		if(readIn=='start') then
!			flag=.TRUE.
!		end if
!	end do
!	flag=.FALSE.

!	read(50,*)	!Blank line in strain input file

!end if

!************************************************
!Step 1: modify global coordinates to include entire simulation volume (not just centers of elements)
!************************************************

!The  boundary coordinates of the system (xmin,xmax,ymin,ymax,zmin,zmax)
do i=1,5,2
	myProc%globalCoord(i)=0d0
end do
myProc%globalCoord(2)=length*dble(numx)
myProc%globalCoord(4)=length*dble(numy)
myProc%globalCoord(6)=length*dble(numz)

!step 2: divide volume among processors according to the factorization above
!This step locates the boundaries of the local processor within the global volume.
!**********
myProc%localCoord(1)=myProc%globalCoord(1)+&
		mod(myProc%taskid,dims(1))*(myProc%globalCoord(2)-myProc%globalCoord(1))/&
				dble(dims(1))

myProc%localCoord(2)=myProc%localCoord(1)+&
		(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(dims(1))

myProc%localCoord(3)=myProc%globalCoord(3)+&
		mod(myProc%taskid/dims(1),dims(2))*&
				(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))

myProc%localCoord(4)=myProc%localCoord(3)+&
		(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))

myProc%localCoord(5)=myProc%GlobalCoord(5)+&
		mod(myProc%taskid/(dims(1)*dims(2)),dims(3))*&
				(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))

myProc%localCoord(6)=myProc%localCoord(5)+&
		(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))

totalVolume=(myProc%localCoord(2)-myProc%localCoord(1))*(myProc%localCoord(4)-myProc%localCoord(3))*&
		(myProc%localCoord(6)-myProc%localCoord(5))

!Step 4: create mesh in each processor of elements only within that processor's local coordinates
!NOTE: this is an actual mesh of volume elements, not the 'processor mesh' created above.

!Step 4b: find how many volume elements are in local processor and allocate myMesh accordingly
!numxLocal, numyLocal, numzLocal are used to determine the size of the mesh inside the local processor.
numxLocal=0
numyLocal=0
numzLocal=0

!get numxLocal
if(myProc%localCoord(1)==myProc%globalCoord(1) .AND. myProc%localCoord(2)==myProc%globalCoord(2)) then
	tempx1=0
	tempx2=numx
	numxLocal=numx
	tempCenter(1)=length/2d0
	tempCenter(2)=myProc%globalCoord(2)-length/2d0

else if(myProc%localCoord(1)==myProc%globalCoord(1)) then	!at xmin

	tempx1=0
	tempx2 = myProc%localCoord(2)/(length/2d0)
	if((dble(tempx2)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2)< (dble(tempx2+1)*(length/2d0))) then
		if(mod(tempx2,2)==0) then
			tempx2 = tempx2/2
		else
			tempx2 = (tempx2+1)/2
		end if
	end if

	numxLocal=tempx2
	tempCenter(1)=length/2d0
	tempCenter(2)=dble(tempx2)*length-length/2d0

else if(myProc%localCoord(2)==myProc%globalCoord(2)) then	!at xmax
	tempx1 = myProc%localCoord(1)/(length/2d0)
	tempx2=numx
	if((dble(tempx1)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(tempx1+1)*(length/2d0))) then
		if(mod(tempx1,2)==0) then
			tempx1 = tempx1/2
		else
			tempx1 = (tempx1+1)/2
		end if
	end if

	numxLocal=numx-tempx1
	tempCenter(1)=dble(tempx1)*length+length/2d0
	tempCenter(2)=myProc%globalCoord(2)-length/2d0

else	!in the middle
	tempx1 = myProc%localCoord(1)/(length/2d0)
	tempx2 = myProc%localCoord(2)/(length/2d0)
	if((dble(tempx1)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(tempx1+1)*(length/2d0))) then
		if(mod(tempx1,2)==0) then
			tempx1 = tempx1/2
		else
			tempx1 = (tempx1+1)/2
		end if
	end if
	if((dble(tempx2)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2) < (dble(tempx2+1)*(length/2d0))) then
		if(mod(tempx2,2)==0) then
			tempx2 = tempx2/2
		else
			tempx2 = (tempx2+1)/2
		end if
	end if

	numxLocal = tempx2-tempx1
	tempCenter(1)=dble(tempx1)*length+length/2d0
	tempCenter(2)=dble(tempx2)*length-length/2d0

end if

!get numyLocal
if(myProc%localCoord(3)==myProc%globalCoord(3) .AND. myProc%localCoord(4)==myProc%globalCoord(4)) then
	tempy1=0
	tempy2=numy
	numyLocal=numy
	tempCenter(3)=length/2d0
	tempCenter(4)=myProc%globalCoord(4)-length/2d0

else if(myProc%localCoord(3)==myProc%globalCoord(3)) then	!at ymin
	tempy1=0
	tempy2 = myProc%localCoord(4)/(length/2d0)
	if((dble(tempy2)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(tempy2+1)*(length/2d0))) then
		if(mod(tempy2,2)==0) then
			tempy2 = tempy2/2
		else
			tempy2 = (tempy2+1)/2
		end if
	end if
	numyLocal=tempy2
	tempCenter(3)=length/2d0
	tempCenter(4)=dble(tempy2)*length-length/2d0

else if(myProc%localCoord(4)==myProc%globalCoord(4)) then	!at ymax
	tempy1 = myProc%localCoord(3)/(length/2d0)
	tempy2=numy
	if((dble(tempy1)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(tempy1+1)*(length/2d0))) then
		if(mod(tempy1,2)==0) then
			tempy1 = tempy1/2
		else
			tempy1 = (tempy1+1)/2
		end if
	end if
	numyLocal=numy-tempy1
	tempCenter(3)=dble(tempy1)*length+length/2d0
	tempCenter(4)=myProc%globalCoord(4)-length/2d0

else	!in the middle
	tempy1 = myProc%localCoord(3)/(length/2d0)
	tempy2 = myProc%localCoord(4)/(length/2d0)
	if((dble(tempy1)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(tempy1+1)*(length/2d0))) then
		if(mod(tempy1,2)==0) then
			tempy1 = tempy1/2
		else
			tempy1 = (tempy1+1)/2
		end if
	end if
	if((dble(tempy2)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(tempy2+1)*(length/2d0))) then
		if(mod(tempy2,2)==0) then
			tempy2 = tempy2/2
		else
			tempy2 = (tempy2+1)/2
		end if
	end if
	numyLocal = tempy2-tempy1
	tempCenter(3)=dble(tempy1)*length+length/2d0
	tempCenter(4)=dble(tempy2)*length-length/2d0

end if

!get numzLocal
if(myProc%localCoord(5)==myProc%globalCoord(5) .AND. myProc%localCoord(6)==myProc%globalCoord(6)) then
	tempz1=0
	tempz2=numz
	numzLocal=numz
	tempCenter(5)=length/2d0
	tempCenter(6)=myProc%globalCoord(6)-length/2d0

else if(myProc%localCoord(5)==myProc%globalCoord(5)) then	!at zmin
	tempz1=0
	tempz2 = myProc%localCoord(6)/(length/2d0)
	if((dble(tempz2)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(tempz2+1)*(length/2d0))) then
		if(mod(tempz2,2)==0) then
			tempz2 = tempz2/2
		else
			tempz2 = (tempz2+1)/2
		end if
	end if
	numzLocal=tempz2
	tempCenter(5)=length/2d0
	tempCenter(6)=dble(tempz2)*length-length/2d0

else if(myProc%localCoord(6)==myProc%globalCoord(6)) then	!at zmax
	tempz1 = myProc%localCoord(5)/(length/2d0)
	tempz2=numz
	if((dble(tempz1)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5)< (dble(tempz1+1)*(length/2d0))) then
		if(mod(tempz1,2)==0) then
			tempz1 = tempz1/2
		else
			tempz1 = (tempz1+1)/2
		end if
	end if
	numzLocal=numz-tempz1
	tempCenter(5)=dble(tempz1)*length+length/2d0
	tempCenter(6)=myProc%globalCoord(6)-length/2d0

else	!in the middle
	tempz1 = myProc%localCoord(5)/(length/2d0)
	tempz2 = myProc%localCoord(6)/(length/2d0)
	if((dble(tempz1)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5) < (dble(tempz1+1)*(length/2d0))) then
		if(mod(tempz1,2)==0) then
			tempz1 = tempz1/2
		else
			tempz1 = (tempz1+1)/2
		end if
	end if
	if((dble(tempz2)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(tempz2+1)*(length/2d0))) then
		if(mod(tempz2,2)==0) then
			tempz2 = tempz2/2
		else
			tempz2 = (tempz2+1)/2
		end if
	end if
	numzLocal = tempz2-tempz1
	tempCenter(5)=dble(tempz1)*length+length/2d0
	tempCenter(6)=dble(tempz2)*length-length/2d0

end if

!this variable is used at various points in SRSCD; lets us know the length fo myMesh(:)
numCells=numxLocal*numyLocal*numzLocal	!total cells of this processor

!If any processors don't have volume elements in them (too many procs), we create an error message
if(numCells==0) then
	write(*,*) 'error processors with no volume elements'
	call MPI_ABORT(comm,ierr)
	stop
end if

if(myProc%taskid==MASTER) then
	write(*,*) 'numx numy numz', numx, numy, numz
end if

write(*,*) 'proc', myProc%taskid, 'numxLocal numyLocal numzLocal', numxLocal, numyLocal, numzLocal

!Step 3c: Read meshes. And for each volume element in myMesh, assign coordinates and material number (and proc number)
!We can't read from the file directly into myMesh, because we don't know how big to make it until
!we read in all of the volume element information. Therefore, we allocate it and then reread from the
!global mesh back into myMesh.
allocate(myMesh(numCells))

localElem=0

do k=1,numzLocal
	do j=1,numyLocal
		do i=1,numxLocal

			localElem=localElem+1

			x = tempx1+i
			y = tempy1+j
			z = tempz1+k

			myMesh(localElem)%coordinates(1)=tempCenter(1)+(i-1)*length
			myMesh(localElem)%coordinates(2)=tempCenter(3)+(j-1)*length
			myMesh(localElem)%coordinates(3)=tempCenter(5)+(k-1)*length
			myMesh(localElem)%material=1
			myMesh(localElem)%proc=myProc%taskid
			myMesh(localElem)%globalCell=(z-1)*numx*numy+(y-1)*numx+x

			myMesh(localElem)%length=length
			myMesh(localElem)%volume=length**3d0

			!uniform mesh: all elements have 1 neighbor in each direction
			allocate(myMesh(localElem)%neighbors(6,1))
			allocate(myMesh(localElem)%neighborProcs(6,1))
			do dir=1,6
				myMesh(localElem)%numNeighbors(dir)=1
			end do


!			if(strainField=='yes') then

!				read(50,*) tempCoord(1), tempCoord(2), tempCoord(3), tempStrain(1), tempStrain(2), &
!								tempStrain(3), tempStrain(4), tempStrain(5), tempStrain(6)

!				if(tempCoord(1) /= tempMeshCoord(1)) then	!not equal
!					write(*,*) 'Error strain mesh does not match normal mesh'
!				else if(tempCoord(2) /= tempMeshCoord(2)) then
!					write(*,*) 'Error strain mesh does not match normal mesh'
!				else if(tempCoord(3) /= tempMeshCoord(3)) then
!					write(*,*) 'Error strain mesh does not match normal mesh'
!				end if
!				if(tempCoord(1) > myProc%localCoord(1) .AND. tempCoord(1) <= myProc%localCoord(2) &
!						.AND. tempCoord(2) > myProc%localCoord(3) .AND. tempCoord(2) <= myProc%localCoord(4) &
!						.AND. tempCoord(3) > myProc%localCoord(5) .AND. tempCoord(3) <= myProc%localCoord(6)) then
!					do l=1,6
!						myMesh(localElem)%strain(l)=tempStrain(l)
!					end do
!				end if
!			end if

!			if(numMaterials==1) then
!				systemVol=systemVol+length**3d0
!			else if(tempMaterial==1) then	!only add to system volume if we are NOT at a grain boundary
!				systemVol=systemVol+length**3d0
!			else
!				!Do nothing
!			end if

		end do
	end do
end do

!Step 3d: assign neighbors and processor numbers for neighbors (connectivity in myMesh) - periodic or free surfaces in +/- z
if(meshType=='periodic') then
	call createConnectLocalPeriodicUniform(length)
else if(meshType=='freeSurfaces') then
	call createConnectLocalFreeSurfUniform(length)
end if

!Close input files
!if(strainField=='yes') then
!	close(50)
!end if
close(80)

end subroutine

!**************************************************************************************************
!These subroutines create LOCAL connectivity. They identify the volume element # and processor # of
!neighboring volume elements for each element. The connectivity scheme is the same as in the global
!case, but neighboring processor numbers are used here.

!>Subroutine create local connectivity (uniform mesh, periodic boundary conditions)
!
!This subroutine creates a local mesh (for the local processor only) as well as the
!connectivity between elements for the case of a uniform cubic mesh.
!And creates myBoundary with elements that are in neighboring processors that bound this one.
!This subroutine has to choose their neighbor element ID numbers and processors.
!
!Inputs: numxLocal, numyLocal, numzLocal, numx, numy, numz
!Output: myMesh, myBoundary
!**************************************************************************************************

subroutine createConnectLocalPeriodicUniform(length)
use DerivedType
use mod_constants

implicit none
include 'mpif.h'
integer cell, localCell, maxElement
double precision length

!buffer lists to send all information at the end
integer i, dir, tag

integer, allocatable :: sendBuffer(:,:,:), recvBuffer(:,:)
integer materialBuff(numxLocal*numyLocal*numzLocal,6)
integer numSend(6), numRecv(6)

integer sendRequest(6), recvRequest(6)
integer sendStatus(6), recvStatus(6)
integer status(MPI_STATUS_SIZE)

integer tempRecv
logical flagProbe

do dir=1,6
	numSend(dir)=0
	numRecv(dir)=0
end do

!allocate(materialStrain(7,numCells,6))
!************************************************
!periodic boundary condition version
!************************************************
allocate(sendBuffer(6,numxLocal*numyLocal*numzLocal,2))

do cell=1,numCells

	!Right (+x)
	if(mod(cell,numxLocal)==0) then !identify cell to the right
		
		!If we are on the right edge of the local mesh, identify the neighboring processor
		myMesh(cell)%neighborProcs(1,1)=myProc%procNeighbor(1)
		
		!If the mesh is only one processor thick, then the neighboring processor is the same as this processor
		if(myMesh(cell)%neighborProcs(1,1)==myProc%taskid) then	
			myMesh(cell)%neighbors(1,1)=cell-numxLocal+1	!use periodic rules from uniform cubic mesh
		else
			!add these items to sendList, a buffer that is used at the end of this subroutine to communicate
			!with other processors via MPI about which elements have neighbors in other cells
			!(this cell sends information about itself, and the nieghboring cell will do the same in its processor)
			numSend(1)=numSend(1)+1				!count the number of items in buffer
			sendBuffer(1,numSend(1),1)=cell		!localCellID in this processor
			sendBuffer(1,numSend(1),2)=myMesh(cell)%material	!Material ID number that this element is composed of
		end if

	else
		!if we are still inside the local mesh, don't need to communicate with neighboring cells and 
		!just use the uniform cubic connectivity rules (increase x, then y, then z)
		myMesh(cell)%neighbors(1,1)=cell+1
		myMesh(cell)%neighborProcs(1,1)=myProc%taskid
	end if
	
	!Left (-x)
	if(mod(cell+numxLocal-1,numxLocal)==0) then !identify cell to the left
		
		myMesh(cell)%neighborProcs(2,1)=myProc%procNeighbor(2)

		if(myMesh(cell)%neighborProcs(2,1)==myProc%taskid) then
			myMesh(cell)%neighbors(2,1)=cell+numxLocal-1
		else
			numSend(2)=numSend(2)+1
			sendBuffer(2,numSend(2),1)=cell
			sendBuffer(2,numSend(2),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(2,1)=cell-1
		myMesh(cell)%neighborProcs(2,1)=myProc%taskid
	end if
	
	!Front (+y)
	if(mod(cell,numxLocal*numyLocal) > numxLocal*(numyLocal-1) .OR. mod(cell,numxLocal*numyLocal)==0) then
		
		myMesh(cell)%neighborProcs(3,1)=myProc%procNeighbor(3)

		if(myMesh(cell)%neighborProcs(3,1)==myProc%taskid) then
			myMesh(cell)%neighbors(3,1)=cell-(numxLocal*(numyLocal-1))
		else
			numSend(3)=numSend(3)+1
			sendBuffer(3,numSend(3),1)=cell
			sendBuffer(3,numSend(3),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(3,1)=cell+numxLocal
		myMesh(cell)%neighborProcs(3,1)=myProc%taskid
	end if
	
	!Back (-y)
	if(mod(cell,numxLocal*numyLocal) <= numxLocal .AND. (mod(cell, numxLocal*numyLocal) /= 0 .OR. numyLocal==1)) then

		myMesh(cell)%neighborProcs(4,1)=myProc%procNeighbor(4)

		if(myMesh(cell)%neighborProcs(4,1)==myProc%taskid) then
			myMesh(cell)%neighbors(4,1)=cell+(numxLocal*(numyLocal-1))
		else
			numSend(4)=numSend(4)+1
			sendBuffer(4,numSend(4),1)=cell
			sendBuffer(4,numSend(4),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(4,1)=cell-numxLocal
		myMesh(cell)%neighborProcs(4,1)=myProc%taskid
	end if
	
	!Up (+z)
	if(mod(cell,numxLocal*numyLocal*numzLocal) > numxLocal*numyLocal*(numzLocal-1) &
			.OR. mod(cell, numxLocal*numyLocal*numzLocal)==0) then

		myMesh(cell)%neighborProcs(5,1)=myProc%procNeighbor(5)

		if(myMesh(cell)%neighborProcs(5,1)==myProc%taskid) then
			myMesh(cell)%neighbors(5,1)=cell-(numxLocal*numyLocal*(numzLocal-1))
		else
			numSend(5)=numSend(5)+1
			sendBuffer(5,numSend(5),1)=cell
			sendBuffer(5,numSend(5),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(5,1)=cell+numxLocal*numyLocal
		myMesh(cell)%neighborProcs(5,1)=myProc%taskid
	end if
	
	!Down (-z)
	if(mod(cell,numxLocal*numyLocal*numzLocal) <= numxLocal*numyLocal &
			.AND. (mod(cell,numxLocal*numyLocal*numzLocal) /= 0 .OR. numzLocal==1)) then

		myMesh(cell)%neighborProcs(6,1)=myProc%procNeighbor(6)

		if(myMesh(cell)%neighborProcs(6,1)==myProc%taskid) then
			myMesh(cell)%neighbors(6,1)=cell+(numxLocal*numyLocal*(numzLocal-1))
		else
			numSend(6)=numSend(6)+1
			sendBuffer(6,numSend(6),1)=cell
			sendBuffer(6,numSend(6),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(6,1)=cell-numxLocal*numyLocal
		myMesh(cell)%neighborProcs(6,1)=myProc%taskid
	end if

end do

!Now that we have created SendList, the buffer that contains information about which elements in myMesh
!have neighbors on other processors, we send out the cell numbers of these cells to the neighboring procs.
!We label each MPI_SEND with globalCell so that the receiving processor pairs cells correctly using globalNeighbor (see below)
maxElement=0

!*******************************************************
!Send

do dir=1,6
	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		call MPI_ISEND(sendBuffer(dir,1:numSend(dir),:), numSend(dir)*2, MPI_INTEGER, myProc%procNeighbor(dir), &
				200+dir, comm, sendRequest(dir),ierr)
	end if
end do

!Recv
do dir=1,6

	if(mod(dir,2)==0) then
		tag = dir-1
	else
		tag = dir+1
	end if

	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		numRecv(dir)=numSend(dir)
		allocate(recvBuffer(numRecv(dir),2))
		call MPI_IRECV(recvBuffer, numRecv(dir)*2, MPI_INTEGER, myProc%procNeighbor(dir), &
				200+tag, comm, recvRequest(dir), ierr)
		call MPI_WAIT(recvRequest(dir),recvStatus(dir),ierr)

		do i=1, numRecv(dir)
			localCell=sendBuffer(dir,i,1)
			myMesh(localCell)%neighbors(dir,1)=recvBuffer(i,1)
			materialBuff(localCell,dir)=recvBuffer(i,2)

		end do
		deallocate(recvBuffer)
	end if
end do

do dir=1,6
	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		call MPI_WAIT(sendRequest(dir),sendStatus(dir),ierr)
	end if
end do

deallocate(sendBuffer)



!***************************************************************************************************
!Initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************

!Step 1: Find the max cell# of any boundary mesh element
maxElement=0
do i=1, numCells
    do dir=1, 6
        if(myMesh(i)%neighborProcs(dir,1) /= myProc%taskid) then	!we are pointing to a different proc
            if(myMesh(i)%neighbors(dir,1) > maxElement) then		!searching for the max element number in a neighbor
                maxElement=myMesh(i)%neighbors(dir,1)
            end if
        end if
    end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(6,maxElement))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do dir=1,6
	do i=1,maxElement
		myBoundary(dir,i)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
		myBoundary(dir,i)%proc=-10			!default, says that this is not a real element of myBoundary.
        myBoundary(dir,i)%material=0
        myBoundary(dir,i)%length=0d0
        myBoundary(dir,i)%volume=0d0
	end do
end do

do i=1,numCells
    do dir=1,6
        if(myMesh(i)%neighborProcs(dir,1) == -1) then										!this is a free surface
            !do nothing
        else if(myMesh(i)%neighborProcs(dir,1) /= myProc%taskid) then
            myBoundary(dir,myMesh(i)%neighbors(dir,1))%proc=myMesh(i)%neighborProcs(dir,1)	!set proc # of elements in myBoundary
            myBoundary(dir,myMesh(i)%neighbors(dir,1))%length=length						!set length of elements in myBoundary
            myBoundary(dir,myMesh(i)%neighbors(dir,1))%volume=length**3d0					!set volume of elements in myBoundary (changes with cascade addition)
			myBoundary(dir,myMesh(i)%neighbors(dir,1))%material=materialBuff(i,dir)
            myBoundary(dir,myMesh(i)%neighbors(dir,1))%localNeighbor=i

!           if(strainField=='yes') then
!                do l=1,6
!                   myBoundary(dir,myMesh(i)%neighbors(dir,1))%strain(l)=globalStrain(globalNeighbor,l)
!                    myBoundary(dir,myMesh(i)%neighbors(dir,1))%strain(l)=materialStrain(1+l,i,dir)
!               end do
!           end if
        end if
    end do
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

subroutine createConnectLocalFreeSurfUniform(length)
use DerivedType

implicit none
include 'mpif.h'
integer cell, maxElement
double precision length

integer globalCell, globalNeighbor, status(MPI_STATUS_SIZE)

!buffer lists to send all information at the end
integer numSendRecv(6), i,j, dir

integer sendRight(2,numyLocal*numzLocal), recvRight(2,numyLocal*numzLocal)
integer sendLeft(2,numyLocal*numzLocal), recvLeft(2,numyLocal*numzLocal)

integer sendFront(2,numxLocal*numzLocal), recvFront(2,numxLocal*numzLocal)
integer sendBack(2,numxLocal*numzLocal), recvBack(2,numxLocal*numzLocal)

integer sendUp(2,numxLocal*numyLocal), recvUp(2,numxLocal*numyLocal)
integer sendDown(2,numxLocal*numyLocal), recvDown(2,numxLocal*numyLocal)

integer materialBuff(numxLocal*numyLocal*numzLocal)	!materialID of boundary meshes

do i=1,6
	numSendRecv(i)=0
end do

!******************************************************************
!free surfaces at z=0 and z=zmax(Global) boundary condition version
!******************************************************************
do cell=1,numCells
	if(mod(cell,numx)==0) then !identify cell to the right
		
		myMesh(cell)%neighborProcs(1,1)=myProc%procNeighbor(1)
		if(myMesh(cell)%neighborProcs(1,1)==myProc%taskid) then
			myMesh(cell)%neighbors(1,1)=cell-numx+1
		else
			numSendRecv(1)=numSendRecv(1)+1		!count the number of items in buffer
			sendRight(1,numSendRecv(1))=cell	!localCellID in this processor
			sendRight(2,numSendRecv(1))=myMesh(cell)%material	!Material ID number that this element is composed of

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
			numSendRecv(2)=numSendRecv(2)+1
			sendLeft(1,numSendRecv(2))=cell
			sendLeft(2,numSendRecv(2))=myMesh(cell)%material	!Material ID number that this element is composed of

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
			numSendRecv(3)=numSendRecv(3)+1
			sendFront(1,numSendRecv(3))=cell
			sendFront(2,numSendRecv(3))=myMesh(cell)%material	!Material ID number that this element is composed of

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
			numSendRecv(4)=numSendRecv(4)+1
			sendBack(1,numSendRecv(4))=cell
			sendBack(2,numSendRecv(4))=myMesh(cell)%material	!Material ID number that this element is composed of

		end if
	else
		myMesh(cell)%neighbors(4,1)=cell-numx
		myMesh(cell)%neighborProcs(4,1)=myProc%taskid
	end if
	
	if(mod(cell,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(cell, numx*numy*numz)==0) then
		myMesh(cell)%neighborProcs(5,1)=myProc%procNeighbor(5)
		globalCell=myMesh(cell)%globalCell
		globalNeighbor=findgNeighborFreeSurf(globalCell, 5)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(5,1)=0
			myMesh(cell)%neighborProcs(5,1)=-1
		else
			if(myMesh(cell)%neighborProcs(5,1)==myProc%taskid) then
				myMesh(cell)%neighbors(5,1)=cell-(numx*numy*(numz-1))
			else
				numSendRecv(5)=numSendRecv(5)+1
				sendUp(1,numSendRecv(5))=cell
				sendUp(2,numSendRecv(5))=myMesh(cell)%material

			end if
		end if
	else
		myMesh(cell)%neighbors(5,1)=cell+numx*numy
		myMesh(cell)%neighborProcs(5,1)=myProc%taskid
	end if
	
	if(mod(cell,numx*numy*numz) <= numx*numy .AND. (mod(cell,numx*numy*numz) /= 0 .OR. numz==1)) then
		myMesh(cell)%neighborProcs(6,1)=myProc%procNeighbor(6)
		globalCell=myMesh(cell)%globalCell
		globalNeighbor=findgNeighborFreeSurf(globalCell, 6)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(6,1)=0
			myMesh(cell)%neighborProcs(6,1)=-1
		else
			if(myMesh(cell)%neighborProcs(6,1)==myProc%taskid) then
				myMesh(cell)%neighbors(6,1)=cell+(numx*numy*(numz-1))
			else
				numSendRecv(6)=numSendRecv(6)+1
				sendDown(1,numSendRecv(6))=cell
				sendDown(2,numSendRecv(6))=myMesh(cell)%material

			end if
		end if
	else
		myMesh(cell)%neighbors(6,1)=cell-numx*numy
		myMesh(cell)%neighborProcs(6,1)=myProc%taskid
	end if
end do

!***************************************************************
!Send
if(myProc%procNeighbor(1)/=myProc%taskid) then	!right
	call MPI_SEND(sendRight, 2*numyLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(1), 1, comm, ierr)
end if

if(myProc%procNeighbor(2)/=myProc%taskid) then	!left
	call MPI_SEND(sendLeft, 2*numyLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(2), 2, comm, ierr)
end if

if(myProc%procNeighbor(3)/=myProc%taskid) then	!front
	call MPI_SEND(sendFront, 2*numxLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(3), 3, comm, ierr)
end if

if(myProc%procNeighbor(4)/=myProc%taskid) then	!back
	call MPI_SEND(sendBack, 2*numxLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(4), 4, comm, ierr)
end if

if(myProc%procNeighbor(5)/=myProc%taskid) then	!up
	call MPI_SEND(sendUp, 2*numxLocal*numyLocal, MPI_INTEGER, myProc%procNeighbor(5), 5, comm, ierr)
end if

if(myProc%procNeighbor(6)/=myProc%taskid) then	!down
	call MPI_SEND(sendDown, 2*numxLocal*numyLocal, MPI_INTEGER, myProc%procNeighbor(6), 6, comm, ierr)
end if

!************************************************************
!Recv
if(myProc%procNeighbor(1)/=myProc%taskid) then	!right

	call MPI_RECV(recvRight, 2*numyLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(1), 2, comm, status, ierr)

	do i=1, numyLocal*numzLocal
		cell=sendRight(1,i)
		myMesh(cell)%neighbors(1,1)=recvRight(1,i)
		materialBuff(cell)=recvRight(2,i)
	end do

end if

if(myProc%procNeighbor(2)/=myProc%taskid) then	!left

	call MPI_RECV(recvLeft, 2*numyLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(2), 1, comm, status, ierr)

	do i=1, numyLocal*numzLocal
		cell=sendLeft(1,i)
		myMesh(cell)%neighbors(2,1)=recvLeft(1,i)
		materialBuff(cell)=recvLeft(2,i)
	end do
end if

if(myProc%procNeighbor(3)/=myProc%taskid) then	!front

	call MPI_RECV(recvFront, 2*numxLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(3), 4, comm, status, ierr)

	do i=1, numxLocal*numzLocal
		cell=sendFront(1,i)
		myMesh(cell)%neighbors(3,1)=recvFront(1,i)
		materialBuff(cell)=recvFront(2,i)
	end do
end if

if(myProc%procNeighbor(4)/=myProc%taskid) then	!back

	call MPI_RECV(recvBack, 2*numxLocal*numzLocal, MPI_INTEGER, myProc%procNeighbor(4), 3, comm, status, ierr)

	do i=1, numxLocal*numzLocal
		cell=sendBack(1,i)
		myMesh(cell)%neighbors(4,1)=recvBack(1,i)	!cellID
		materialBuff(cell)=recvBack(2,i)	!materialID
	end do
end if

if(myProc%procNeighbor(5)/=myProc%taskid) then	!up

	call MPI_RECV(recvUp, 2*numxLocal*numyLocal, MPI_INTEGER, myProc%procNeighbor(5), 6, comm, status, ierr)

	do i=1, numxLocal*numyLocal
		cell=sendUp(1,i)
		myMesh(cell)%neighbors(5,1)=recvUp(1,i)	!cellID
		materialBuff(cell)=recvUp(2,i)	!materialID
	end do
end if

if(myProc%procNeighbor(6)/=myProc%taskid) then	!down

	call MPI_RECV(recvDown, 2*numxLocal*numyLocal, MPI_INTEGER, myProc%procNeighbor(6), 5, comm, status, ierr)

	do i=1, numxLocal*numyLocal
		cell=sendDown(1,i)
		myMesh(cell)%neighbors(6,1)=recvDown(1,i)	!cellID
		materialBuff(cell)=recvDown(2,i)	!materialID
	end do
end if

!***************************************************************************************************
!Initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************

!Step 1: Find the max cell# of any boundary mesh element
maxElement=0
do i=1, numCells
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

do i=1,numCells
	do j=1,6
		if(myMesh(i)%neighborProcs(j,1) == -1) then										!this is a free surface
			!do nothing
		else if(myMesh(i)%neighborProcs(j,1) /= myProc%taskid .AND. myMesh(i)%neighborProcs(j,1) /= -1) then
			myBoundary(j,myMesh(i)%neighbors(j,1))%proc=myMesh(i)%neighborProcs(j,1)	!set proc # of elements in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%length=length						!set length of elements in myBoundary
			myBoundary(j,myMesh(i)%neighbors(j,1))%volume=length**3d0					!set volume of elements in myBoundary (changes with cascade addition)
			myBoundary(j,myMesh(i)%neighbors(j,1))%material=materialBuff(i)
			myBoundary(j,myMesh(i)%neighbors(j,1))%localNeighbor=i

!           if(strainField=='yes') then
!                do l=1,6
!                   myBoundary(j,myMesh(i)%neighbors(j,1))%strain(l)=globalStrain(globalNeighbor,l)
!                    myBoundary(j,myMesh(i)%neighbors(j,1))%strain(l)=materialStrain(1+l,i,j)
!               end do
!           end if
		end if
	end do
end do


end subroutine

!***************************************************************************************
!>Function find global neighbor (uniform mesh, periodic boundary conditions)
!
!Inputs: globalID, direction
!Output: globalNeighborID of the globalID in the direaction
!***************************************************************************************

integer function findgNeighborPeriodic(globalID, dir)
use DerivedType
use mod_constants
implicit none

integer globalID, dir, neighborID

!************************************************
!periodic boundary condition version
!************************************************

if(dir==1) then
    if(mod(globalID,numx)==0) then !identify cell to the right
        neighborID=globalID-numx+1
    else
        neighborID=globalID+1
    end if
else if(dir==2) then
    if(mod(globalID+numx-1,numx)==0) then !identify cell to the left
        neighborID=globalID+numx-1
    else
        neighborID=globalID-1
    end if
else if(dir==3) then
    if(mod(globalID,numx*numy) > numx*(numy-1) .OR. mod(globalID,numx*numy)==0) then
        neighborID=globalID-(numx*(numy-1))
    else
        neighborID=globalID+numx
    end if
else if(dir==4) then
    if(mod(globalID,numx*numy) <= numx .AND. (mod(globalID, numx*numy) /= 0 .OR. numy==1)) then
        neighborID=globalID+(numx*(numy-1))
    else
        neighborID=globalID-numx
    end if
else if(dir==5) then
    if(mod(globalID,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(globalID, numx*numy*numz)==0) then
        neighborID=globalID-(numx*numy*(numz-1))
    else
        neighborID=globalID+numx*numy
    end if
else if(dir==6) then
    if(mod(globalID,numx*numy*numz) <= numx*numy .AND. (mod(globalID,numx*numy*numz) /= 0 .OR. numz==1)) then
        neighborID=globalID+(numx*numy*(numz-1))
    else
        neighborID=globalID-numx*numy
    end if
end if

findgNeighborPeriodic=neighborID

end function

!***************************************************************************************
!>Function find global neighbor (uniform mesh, free surfaces in z-directions and periodic boundary conditions in other directions)
!
!Inputs: globalID, direction
!Output: globalNeighborID of the globalID in the direaction
!***************************************************************************************

integer function findgNeighborFreeSurf(globalID, dir)
use DerivedType
use mod_constants
implicit none

integer globalID, dir, neighborID

!************************************************
!PBCs in x and y, free in z (cell 0 represents free surface)
!************************************************

if(dir==1) then
    if(mod(globalID,numx)==0) then !identify cell to the right
        neighborID=globalID-numx+1
    else
        neighborID=globalID+1
    end if
else if(dir==2) then
    if(mod(globalID+numx-1,numx)==0) then !identify cell to the left
        neighborID=globalID+numx-1
    else
        neighborID=globalID-1
    end if
else if(dir==3) then
    if(mod(globalID,numx*numy) > numx*(numy-1) .OR. mod(globalID,numx*numy)==0) then
        neighborID=globalID-(numx*(numy-1))
    else
        neighborID=globalID+numx
    end if
else if(dir==4) then
    if(mod(globalID,numx*numy) <= numx .AND. (mod(globalID, numx*numy) /= 0 .OR. numy==1)) then
        neighborID=globalID+(numx*(numy-1))
    else
        neighborID=globalID-numx
    end if
else if(dir==5) then
    if(mod(globalID,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(globalID, numx*numy*numz)==0) then
        neighborID=0
    else
        neighborID=globalID+numx*numy
    end if
else if(dir==6) then
    if(mod(globalID,numx*numy*numz) <= numx*numy .AND. (mod(globalID,numx*numy*numz) /=0 .OR. numz==1)) then
        neighborID=0
    else
        neighborID=globalID-numx*numy
    end if
end if

findgNeighborFreeSurf=neighborID

end function


!*********************************************************************************************************************
!The following subroutine and function are used for non-uniform meshes
!*********************************************************************************************************************

!***************************************************************************************************
!>Subroutine read Mesh Non Uniform - creates processor and mesh files from non-uniform mesh input
!
!This is the main subroutine that controls reading the non-uniform mesh. It reads in the mesh
!from a file and divides the volume elements between the processors in the parallel
!simulation. If a division such that each processor has at least one element is not possible,
!then the subroutine returns an error. It also creates a 'processor mesh', indicating
!the bounds of each processor. NOTE: this subroutine is for non-uniform meshes only (and these have not been
!implemented in the rest of the program)
!
!Input: filename of mesh file
!Outputs: myProc (processors with meshes), myMesh, and myBoundary
!***************************************************************************************************

!main subroutine that creates mesh (non-uniform cubic)
subroutine readMeshNonUniform(filename)
use mod_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), procDivision(3), i, j, k, matNum
character*20 readIn, meshType
character*50 filename, filename2, filename3
logical flag
double precision volumeFaces(3), totalArea, length
integer elem, localElements(3), localElem,  maxNumNeighbors, globalMaxNeighbors
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

	subroutine createConnectGlobalPeriodic(Connect, Coord, NumNeighbors, Length, bdryCoord)
	double precision, allocatable :: Coord(:,:), Length(:)
	integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
	double precision bdryCoord(6)
	end subroutine

	subroutine createConnectGlobalFreeSurf(Connect, Coord, NumNeighbors, Length,  bdryCoord)
	double precision, allocatable :: Coord(:,:), Length(:)
	integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
	double precision bdryCoord(6)
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
	call MPI_ABORT(comm,ierr)
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
		 myProc%globalCoord)
else if(meshType=='freeSurfaces') then
	call createConnectGlobalFreeSurfNonUniform(globalMeshConnect, globalMeshCoord, globalNumNeighbors, globalLength, &
		 myProc%globalCoord)
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

subroutine createConnectGlobalPeriodicNonUniform(Connect, Coord, NumNeighbors, Length,  bdryCoord)
implicit none

double precision, allocatable :: Coord(:,:), Length(:)
integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
integer  cell, dir, i, searchElem, neighbor
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

subroutine createConnectGlobalFreeSurfNonUniform(Connect, Coord, NumNeighbors, Length,  bdryCoord)
implicit none

double precision, allocatable :: Coord(:,:), Length(:)
integer, allocatable :: Connect(:,:,:), NumNeighbors(:,:)
integer  cell, dir, i, searchElem, neighbor
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
use mod_constants
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
	call MPI_SEND(elem, 1, MPI_INTEGER, myMesh(elem)%neighborProcs(dir,neighbor), globalCell, comm, ierr)
end do

do i=1,numSendRecv
	elem=sendList(i,1)
	dir=sendList(i,2)
	globalNeighbor=sendList(i,4)
	neighbor=sendList(i,5)
	call MPI_RECV(myMesh(elem)%neighbors(dir,neighbor), 1, MPI_INTEGER, myMesh(elem)%neighborProcs(dir,neighbor), &
				globalNeighbor, comm, status, ierr)
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
	do while (flag .eqv. .FALSE.)
		i=i+1
		if(coord(1)==gCoord(i,1) .AND. coord(2)==gCoord(i,2) .AND. coord(3)==gCoord(i,3)) then
			flag=.TRUE.
		endif
	end do

	findGlobalCell=i
end function


!Find the cell number in the local mesh using coordinates
!>Function find local cell
!!
!!Finds the cell ID number in the local mesh using the coordinates of the centroid of the cell
!!
!!Input: coordinates
!!Output: cell ID number

integer function findLocalCell(coord)
use DerivedType
use mod_constants
implicit none

double precision coord(3)
integer i
logical flag

flag = .FALSE.
i=0
do while(flag .eqv. .FALSE.)
	i=i+1
	if(coord(1)==myMesh(i)%coordinates(1) .AND. coord(2)==myMesh(i)%coordinates(2) .AND. coord(3)==myMesh(i)%coordinates(3)) then
		flag=.TRUE.
	endif
end do

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
use mod_constants
implicit none

integer procDivision(3), i
double precision, allocatable :: procCoordList(:,:)

do i=1,myProc%numtasks
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
end do
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
use mod_constants
implicit none

double precision, allocatable :: globalMeshCoord(:,:)
double precision, allocatable :: procCoordList(:,:)
integer elem, i

do i=1,myProc%numtasks
	if(globalMeshCoord(elem,1) > procCoordList(i,1) .AND. globalMeshCoord(elem,1) <= procCoordList(i,2) &
	.AND. globalMeshCoord(elem,2) > procCoordList(i,3) .AND. globalMeshCoord(elem,2) <= procCoordList(i,4) &
	.AND. globalMeshCoord(elem,3) > procCoordList(i,5) .AND. globalMeshCoord(elem,3) <= procCoordList(i,6)) then
		findNeighborProc=i-1
		exit
	endif
end do

end function


end module

