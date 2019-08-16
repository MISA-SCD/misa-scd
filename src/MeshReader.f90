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
integer numXmin,numXmax,numYmin,numYmax,numZmin,numZmax	!used to determine numxLocal, numyLocal, numzLocal

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
	numXmin=0
	numXmax = myProc%localCoord(2)/(length/2d0)
	if((dble(numXmax)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2)< (dble(numXmax+1)*(length/2d0))) then
		if(mod(numXmax,2)==0) then
			numXmax = numXmax/2
		else
			numXmax = (numXmax+1)/2
		end if
	end if
	numxLocal=numXmax
else if(myProc%localCoord(2)==myProc%globalCoord(2)) then	!at xmax
	numXmax=0
	numXmin = myProc%localCoord(1)/(length/2d0)
	if((dble(numXmin)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(numXmin+1)*(length/2d0))) then
		if(mod(numXmin,2)==0) then
			numXmin = numXmin/2
		else
			numXmin = (numXmin+1)/2
		end if
	end if
	numxLocal=numx-numXmin

else	!in the middle
	numXmin = myProc%localCoord(1)/(length/2d0)
	numXmax = myProc%localCoord(2)/(length/2d0)
	if((dble(numXmin)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(numXmin+1)*(length/2d0))) then
		if(mod(numXmin,2)==0) then
			numXmin = numXmin/2
		else
			numXmin = (numXmin+1)/2
		end if
	end if
	if((dble(numXmax)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2) < (dble(numXmax+1)*(length/2d0))) then
		if(mod(numXmax,2)==0) then
			numXmax = numXmax/2
		else
			numXmax = (numXmax+1)/2
		end if
	end if
	numxLocal = numXmax-numXmin

end if

!get numyLocal
if(myProc%localCoord(3)==myProc%globalCoord(3) .AND. myProc%localCoord(4)==myProc%globalCoord(4)) then
	numyLocal=numy
else if(myProc%localCoord(3)==myProc%globalCoord(3)) then	!at ymin
	numYmin=0
	numYmax = myProc%localCoord(4)/(length/2d0)
	if((dble(numYmax)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(numYmax+1)*(length/2d0))) then
		if(mod(numYmax,2)==0) then
			numYmax = numYmax/2
		else
			numYmax = (numYmax+1)/2
		end if
	end if
	numyLocal=numYmax
else if(myProc%localCoord(4)==myProc%globalCoord(4)) then	!at ymax
	numYmax=0
	numYmin = myProc%localCoord(3)/(length/2d0)
	if((dble(numYmin)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(numYmin+1)*(length/2d0))) then
		if(mod(numYmin,2)==0) then
			numYmin = numYmin/2
		else
			numYmin = (numYmin+1)/2
		end if
	end if
	numyLocal=numy-numYmin

else	!in the middle
	numYmin = myProc%localCoord(3)/(length/2d0)
	numYmax = myProc%localCoord(4)/(length/2d0)
	if((dble(numYmin)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(numYmin+1)*(length/2d0))) then
		if(mod(numYmin,2)==0) then
			numYmin = numYmin/2
		else
			numYmin = (numYmin+1)/2
		end if
	end if
	if((dble(numYmax)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(numYmax+1)*(length/2d0))) then
		if(mod(numYmax,2)==0) then
			numYmax = numYmax/2
		else
			numYmax = (numYmax+1)/2
		end if
	end if
	numyLocal = numYmax-numYmin

end if

!get numzLocal
if(myProc%localCoord(5)==myProc%globalCoord(5) .AND. myProc%localCoord(6)==myProc%globalCoord(6)) then
	numzLocal=numz
else if(myProc%localCoord(5)==myProc%globalCoord(5)) then	!at zmin
	numZmin=0
	numZmax = myProc%localCoord(6)/(length/2d0)
	if((dble(numZmax)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(numZmax+1)*(length/2d0))) then
		if(mod(numZmax,2)==0) then
			numZmax = numZmax/2
		else
			numZmax = (numZmax+1)/2
		end if
	end if
	numzLocal=numZmax
else if(myProc%localCoord(6)==myProc%globalCoord(6)) then	!at zmax
	numZmax=0
	numZmin = myProc%localCoord(5)/(length/2d0)
	if((dble(numZmin)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5)< (dble(numZmin+1)*(length/2d0))) then
		if(mod(numZmin,2)==0) then
			numZmin = numZmin/2
		else
			numZmin = (numZmin+1)/2
		end if
	end if
	numzLocal=numz-numZmin

else	!in the middle
	numZmin = myProc%localCoord(5)/(length/2d0)
	numZmax = myProc%localCoord(6)/(length/2d0)
	if((dble(numZmin)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5) < (dble(numZmin+1)*(length/2d0))) then
		if(mod(numZmin,2)==0) then
			numZmin = numZmin/2
		else
			numZmin = (numZmin+1)/2
		end if
	end if
	if((dble(numZmax)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(numZmax+1)*(length/2d0))) then
		if(mod(numZmax,2)==0) then
			numZmax = numZmax/2
		else
			numZmax = (numZmax+1)/2
		end if
	end if
	numzLocal = numZmax-numZmin

end if

!this variable is used at various points in SRSCD; lets us know the length fo myMesh(:)
numCells=numxLocal*numyLocal*numzLocal	!total cells of this processor

!If any processors don't have volume elements in them (too many procs), we create an error message
if(numCells==0) then
	write(*,*) 'error processors with no volume elements'
	call MPI_ABORT(comm,ierr)
	stop
end if

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
				allocate(myMesh(localElem)%neighbors(1,6))
				allocate(myMesh(localElem)%neighborProcs(1,6))
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
double precision length
integer localElem, grain
integer x, y, z, numXmin,numXmax,numYmin,numYmax,numZmin,numZmax	!used to determine numxLocal, numyLocal, numzLocal

double precision tempCenter(6)

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
	numXmin=0
	numXmax=numx
	numxLocal=numx
	tempCenter(1)=length/2d0
	tempCenter(2)=myProc%globalCoord(2)-length/2d0

else if(myProc%localCoord(1)==myProc%globalCoord(1)) then	!at xmin

	numXmin=0
	numXmax = myProc%localCoord(2)/(length/2d0)
	if((dble(numXmax)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2)< (dble(numXmax+1)*(length/2d0))) then
		if(mod(numXmax,2)==0) then
			numXmax = numXmax/2
		else
			numXmax = (numXmax+1)/2
		end if
	end if

	numxLocal=numXmax
	tempCenter(1)=length/2d0
	tempCenter(2)=dble(numXmax)*length-length/2d0

else if(myProc%localCoord(2)==myProc%globalCoord(2)) then	!at xmax
	numXmin = myProc%localCoord(1)/(length/2d0)
	numXmax=numx
	if((dble(numXmin)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(numXmin+1)*(length/2d0))) then
		if(mod(numXmin,2)==0) then
			numXmin = numXmin/2
		else
			numXmin = (numXmin+1)/2
		end if
	end if

	numxLocal=numx-numXmin
	tempCenter(1)=dble(numXmin)*length+length/2d0
	tempCenter(2)=myProc%globalCoord(2)-length/2d0

else	!in the middle
	numXmin = myProc%localCoord(1)/(length/2d0)
	numXmax = myProc%localCoord(2)/(length/2d0)
	if((dble(numXmin)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < (dble(numXmin+1)*(length/2d0))) then
		if(mod(numXmin,2)==0) then
			numXmin = numXmin/2
		else
			numXmin = (numXmin+1)/2
		end if
	end if
	if((dble(numXmax)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2) < (dble(numXmax+1)*(length/2d0))) then
		if(mod(numXmax,2)==0) then
			numXmax = numXmax/2
		else
			numXmax = (numXmax+1)/2
		end if
	end if

	numxLocal = numXmax-numXmin
	tempCenter(1)=dble(numXmin)*length+length/2d0
	tempCenter(2)=dble(numXmax)*length-length/2d0

end if

!get numyLocal
if(myProc%localCoord(3)==myProc%globalCoord(3) .AND. myProc%localCoord(4)==myProc%globalCoord(4)) then
	numYmin=0
	numYmax=numy
	numyLocal=numy
	tempCenter(3)=length/2d0
	tempCenter(4)=myProc%globalCoord(4)-length/2d0

else if(myProc%localCoord(3)==myProc%globalCoord(3)) then	!at ymin
	numYmin=0
	numYmax = myProc%localCoord(4)/(length/2d0)
	if((dble(numYmax)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(numYmax+1)*(length/2d0))) then
		if(mod(numYmax,2)==0) then
			numYmax = numYmax/2
		else
			numYmax = (numYmax+1)/2
		end if
	end if
	numyLocal=numYmax
	tempCenter(3)=length/2d0
	tempCenter(4)=dble(numYmax)*length-length/2d0

else if(myProc%localCoord(4)==myProc%globalCoord(4)) then	!at ymax
	numYmin = myProc%localCoord(3)/(length/2d0)
	numYmax=numy
	if((dble(numYmin)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(numYmin+1)*(length/2d0))) then
		if(mod(numYmin,2)==0) then
			numYmin = numYmin/2
		else
			numYmin = (numYmin+1)/2
		end if
	end if
	numyLocal=numy-numYmin
	tempCenter(3)=dble(numYmin)*length+length/2d0
	tempCenter(4)=myProc%globalCoord(4)-length/2d0

else	!in the middle
	numYmin = myProc%localCoord(3)/(length/2d0)
	numYmax = myProc%localCoord(4)/(length/2d0)
	if((dble(numYmin)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < (dble(numYmin+1)*(length/2d0))) then
		if(mod(numYmin,2)==0) then
			numYmin = numYmin/2
		else
			numYmin = (numYmin+1)/2
		end if
	end if
	if((dble(numYmax)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < (dble(numYmax+1)*(length/2d0))) then
		if(mod(numYmax,2)==0) then
			numYmax = numYmax/2
		else
			numYmax = (numYmax+1)/2
		end if
	end if
	numyLocal = numYmax-numYmin
	tempCenter(3)=dble(numYmin)*length+length/2d0
	tempCenter(4)=dble(numYmax)*length-length/2d0

end if

!get numzLocal
if(myProc%localCoord(5)==myProc%globalCoord(5) .AND. myProc%localCoord(6)==myProc%globalCoord(6)) then
	numZmin=0
	numZmax=numz
	numzLocal=numz
	tempCenter(5)=length/2d0
	tempCenter(6)=myProc%globalCoord(6)-length/2d0

else if(myProc%localCoord(5)==myProc%globalCoord(5)) then	!at zmin
	numZmin=0
	numZmax = myProc%localCoord(6)/(length/2d0)
	if((dble(numZmax)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(numZmax+1)*(length/2d0))) then
		if(mod(numZmax,2)==0) then
			numZmax = numZmax/2
		else
			numZmax = (numZmax+1)/2
		end if
	end if
	numzLocal=numZmax
	tempCenter(5)=length/2d0
	tempCenter(6)=dble(numZmax)*length-length/2d0

else if(myProc%localCoord(6)==myProc%globalCoord(6)) then	!at zmax
	numZmin = myProc%localCoord(5)/(length/2d0)
	numZmax=numz
	if((dble(numZmin)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5)< (dble(numZmin+1)*(length/2d0))) then
		if(mod(numZmin,2)==0) then
			numZmin = numZmin/2
		else
			numZmin = (numZmin+1)/2
		end if
	end if
	numzLocal=numz-numZmin
	tempCenter(5)=dble(numZmin)*length+length/2d0
	tempCenter(6)=myProc%globalCoord(6)-length/2d0

else	!in the middle
	numZmin = myProc%localCoord(5)/(length/2d0)
	numZmax = myProc%localCoord(6)/(length/2d0)
	if((dble(numZmin)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5) < (dble(numZmin+1)*(length/2d0))) then
		if(mod(numZmin,2)==0) then
			numZmin = numZmin/2
		else
			numZmin = (numZmin+1)/2
		end if
	end if
	if((dble(numZmax)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < (dble(numZmax+1)*(length/2d0))) then
		if(mod(numZmax,2)==0) then
			numZmax = numZmax/2
		else
			numZmax = (numZmax+1)/2
		end if
	end if
	numzLocal = numZmax-numZmin
	tempCenter(5)=dble(numZmin)*length+length/2d0
	tempCenter(6)=dble(numZmax)*length-length/2d0

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

			x = numXmin+i
			y = numYmin+j
			z = numZmin+k

			myMesh(localElem)%coordinates(1)=tempCenter(1)+(i-1)*length
			myMesh(localElem)%coordinates(2)=tempCenter(3)+(j-1)*length
			myMesh(localElem)%coordinates(3)=tempCenter(5)+(k-1)*length
			myMesh(localElem)%material=1
			myMesh(localElem)%proc=myProc%taskid
			myMesh(localElem)%globalCell=(z-1)*numx*numy+(y-1)*numx+x

			myMesh(localElem)%length=length
			myMesh(localElem)%volume=length**3d0

			!uniform mesh: all elements have 1 neighbor in each direction
			allocate(myMesh(localElem)%neighbors(1,6))
			allocate(myMesh(localElem)%neighborProcs(1,6))
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

integer, allocatable :: send(:,:,:), recv(:,:,:)
integer, allocatable :: sendBuffer(:,:), recvBuffer(:,:)
integer materialBuff(numCells,6)
integer :: numSend(6)=0, numRecv(6)=0

integer sendRequest(6), recvRequest(6)
integer sendStatus(MPI_STATUS_SIZE,6), recvStatus(MPI_STATUS_SIZE,6)
integer status(MPI_STATUS_SIZE)

integer tempRecv
logical flagProbe

!allocate(materialStrain(7,numCells,6))
!************************************************
!periodic boundary condition version
!************************************************
allocate(send(2, numxLocal*numyLocal*numzLocal,6))
allocate(recv(2, numxLocal*numyLocal*numzLocal,6))

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
			send(1,numSend(1),1)=cell		!localCellID in this processor
			send(2,numSend(1),1)=myMesh(cell)%material	!Material ID number that this element is composed of
		end if

	else
		!if we are still inside the local mesh, don't need to communicate with neighboring cells and 
		!just use the uniform cubic connectivity rules (increase x, then y, then z)
		myMesh(cell)%neighbors(1,1)=cell+1
		myMesh(cell)%neighborProcs(1,1)=myProc%taskid
	end if
	
	!Left (-x)
	if(mod(cell+numxLocal-1,numxLocal)==0) then !identify cell to the left
		
		myMesh(cell)%neighborProcs(1,2)=myProc%procNeighbor(2)

		if(myMesh(cell)%neighborProcs(1,2)==myProc%taskid) then
			myMesh(cell)%neighbors(1,2)=cell+numxLocal-1
		else
			numSend(2)=numSend(2)+1
			send(1,numSend(2),2)=cell
			send(2,numSend(2),2)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(1,2)=cell-1
		myMesh(cell)%neighborProcs(1,2)=myProc%taskid
	end if
	
	!Front (+y)
	if(mod(cell,numxLocal*numyLocal) > numxLocal*(numyLocal-1) .OR. mod(cell,numxLocal*numyLocal)==0) then
		
		myMesh(cell)%neighborProcs(1,3)=myProc%procNeighbor(3)

		if(myMesh(cell)%neighborProcs(1,3)==myProc%taskid) then
			myMesh(cell)%neighbors(1,3)=cell-(numxLocal*(numyLocal-1))
		else
			numSend(3)=numSend(3)+1
			send(1,numSend(3),3)=cell
			send(2,numSend(3),3)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(1,3)=cell+numxLocal
		myMesh(cell)%neighborProcs(1,3)=myProc%taskid
	end if
	
	!Back (-y)
	if(mod(cell,numxLocal*numyLocal) <= numxLocal .AND. (mod(cell, numxLocal*numyLocal) /= 0 .OR. numyLocal==1)) then

		myMesh(cell)%neighborProcs(1,4)=myProc%procNeighbor(4)

		if(myMesh(cell)%neighborProcs(1,4)==myProc%taskid) then
			myMesh(cell)%neighbors(1,4)=cell+(numxLocal*(numyLocal-1))
		else
			numSend(4)=numSend(4)+1
			send(1,numSend(4),4)=cell
			send(2,numSend(4),4)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(1,4)=cell-numxLocal
		myMesh(cell)%neighborProcs(1,4)=myProc%taskid
	end if
	
	!Up (+z)
	if(mod(cell,numxLocal*numyLocal*numzLocal) > numxLocal*numyLocal*(numzLocal-1) &
			.OR. mod(cell, numxLocal*numyLocal*numzLocal)==0) then

		myMesh(cell)%neighborProcs(1,5)=myProc%procNeighbor(5)

		if(myMesh(cell)%neighborProcs(1,5)==myProc%taskid) then
			myMesh(cell)%neighbors(1,5)=cell-(numxLocal*numyLocal*(numzLocal-1))
		else
			numSend(5)=numSend(5)+1
			send(1,numSend(5),5)=cell
			send(2,numSend(5),5)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(1,5)=cell+numxLocal*numyLocal
		myMesh(cell)%neighborProcs(1,5)=myProc%taskid
	end if
	
	!Down (-z)
	if(mod(cell,numxLocal*numyLocal*numzLocal) <= numxLocal*numyLocal &
			.AND. (mod(cell,numxLocal*numyLocal*numzLocal) /= 0 .OR. numzLocal==1)) then

		myMesh(cell)%neighborProcs(1,6)=myProc%procNeighbor(6)

		if(myMesh(cell)%neighborProcs(1,6)==myProc%taskid) then
			myMesh(cell)%neighbors(1,6)=cell+(numxLocal*numyLocal*(numzLocal-1))
		else
			numSend(6)=numSend(6)+1
			send(1,numSend(6),6)=cell
			send(2,numSend(6),6)=myMesh(cell)%material
		end if

	else
		myMesh(cell)%neighbors(1,6)=cell-numxLocal*numyLocal
		myMesh(cell)%neighborProcs(1,6)=myProc%taskid
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

		call MPI_SEND(send(:,1:numSend(dir),dir),numSend(dir)*2,MPI_INTEGER,myProc%procNeighbor(dir),200+dir,&
				comm,ierr)

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

	!	tempRecv=0
	!	flagProbe=.FALSE.
	!	call MPI_IPROBE(myProc%procNeighbor(dir), 200+tag,comm,flagProbe,status)

	!	if(flagProbe .eqv. .TRUE.) then
	!		call MPI_GET_COUNT(status,MPI_INTEGER,tempRecv)
	!		numRecv(dir)=tempRecv/2
	!	end if
		numRecv(dir)=numSend(dir)
	!	allocate(recvBuffer(2,numRecv(dir)))
		call MPI_RECV(recv(:,1:numRecv(dir),dir),numRecv(dir)*2,MPI_INTEGER,myProc%procNeighbor(dir),200+tag,&
				comm,status,ierr)
		do i=1, numRecv(dir)
			localCell=send(1,i,dir)
		!	myMesh(localCell)%neighbors(1,dir)=recvBuffer(1,i)
		!	materialBuff(localCell,dir)=recvBuffer(2,i)
			myMesh(localCell)%neighbors(1,dir)=recv(1,i,dir)
			materialBuff(localCell,dir)=recv(2,i,dir)
		end do
	!	deallocate(recvBuffer)
	end if
end do

if(myProc%taskid==MASTER) then
	write(*,*) 'neighbors'
	do i=1,numCells
		write(*,*) myMesh(i)%neighbors
	end do
end if

!***************************************************************************************************
!Initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************

!Step 1: Find the max cell# of any boundary mesh element
maxElement=0
do i=1, numCells
    do dir=1, 6
        if(myMesh(i)%neighborProcs(1,dir) /= myProc%taskid) then	!we are pointing to a different proc
            if(myMesh(i)%neighbors(1,dir) > maxElement) then		!searching for the max element number in a neighbor
                maxElement=myMesh(i)%neighbors(1,dir)
            end if
        end if
    end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(maxElement,6))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do dir=1,6
	do i=1,maxElement
		myBoundary(i,dir)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
		myBoundary(i,dir)%proc=-10			!default, says that this is not a real element of myBoundary.
        myBoundary(i,dir)%material=0
        myBoundary(i,dir)%length=0d0
        myBoundary(i,dir)%volume=0d0
	end do
end do

do i=1,numCells
    do dir=1,6
        if(myMesh(i)%neighborProcs(1,dir) == -1) then										!this is a free surface
            !do nothing
        else if(myMesh(i)%neighborProcs(1,dir) /= myProc%taskid) then
            myBoundary(myMesh(i)%neighbors(1,dir),dir)%proc=myMesh(i)%neighborProcs(1,dir)	!set proc # of elements in myBoundary
            myBoundary(myMesh(i)%neighbors(1,dir),dir)%length=length						!set length of elements in myBoundary
            myBoundary(myMesh(i)%neighbors(1,dir),dir)%volume=length**3d0					!set volume of elements in myBoundary (changes with cascade addition)
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%material=materialBuff(i,dir)
            myBoundary(myMesh(i)%neighbors(1,dir),dir)%localNeighbor=i

!           if(strainField=='yes') then
!                do l=1,6
!                   myBoundary(myMesh(i)%neighbors(1,dir),dir)%strain(l)=globalStrain(globalNeighbor,l)
!                    myBoundary(myMesh(i)%neighbors(1,dir),dir)%strain(l)=materialStrain(1+l,i,dir)
!               end do
!           end if
        end if
    end do
end do

if(myProc%taskid==MASTER) then
	write(*,*) 'maxElement',maxElement
	do dir=1,6
		do i=1,maxElement
			write(*,*) 'dir',dir,'i',i,myBoundary(i,dir)%localNeighbor
		end do
	end do
end if

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
integer cell, maxElement,localCell
double precision length

integer globalCell, globalNeighbor, status(MPI_STATUS_SIZE)

!buffer lists to send all information at the end
integer i,j, dir, tag
integer numSend(6), numRecv(6)

integer, allocatable :: sendBuffer(:,:,:)
integer, allocatable :: recvBuffer(:,:)

integer materialBuff(numxLocal*numyLocal*numzLocal,6)	!materialID of boundary meshes

integer sendRequest(6), recvRequest(6)
integer sendStatus(MPI_STATUS_SIZE,6), recvStatus(MPI_STATUS_SIZE,6)

integer tempRecv
logical flagProbe

numSend(1:6)=0
numRecv(1:6)=0
!******************************************************************
!free surfaces at z=0 and z=zmax(Global) boundary condition version
!******************************************************************
allocate(sendBuffer(2, numxLocal*numyLocal*numzLocal,6))

do cell=1,numCells
	if(mod(cell,numxLocal)==0) then !identify cell to the right
		
		myMesh(cell)%neighborProcs(1,1)=myProc%procNeighbor(1)
		if(myMesh(cell)%neighborProcs(1,1)==myProc%taskid) then
			myMesh(cell)%neighbors(1,1)=cell-numxLocal+1
		else
			numSend(1)=numSend(1)+1				!count the number of items in buffer
			sendBuffer(1,numSend(1),1)=cell		!localCellID in this processor
			sendBuffer(2,numSend(1),1)=myMesh(cell)%material

		end if
	else
		myMesh(cell)%neighbors(1,1)=cell+1
		myMesh(cell)%neighborProcs(1,1)=myProc%taskid
	end if
	
	if(mod(cell+numxLocal-1,numxLocal)==0) then !identify cell to the left
		
		myMesh(cell)%neighborProcs(1,2)=myProc%procNeighbor(2)
		if(myMesh(cell)%neighborProcs(1,2)==myProc%taskid) then
			myMesh(cell)%neighbors(1,2)=cell+numxLocal-1
		else
			numSend(2)=numSend(2)+1
			sendBuffer(1,numSend(2),2)=cell
			sendBuffer(2,numSend(2),2)=myMesh(cell)%material

		end if
	else
		myMesh(cell)%neighbors(1,2)=cell-1
		myMesh(cell)%neighborProcs(1,2)=myProc%taskid
	end if
	
	if(mod(cell,numxLocal*numyLocal) > numxLocal*(numyLocal-1) .OR. mod(cell,numxLocal*numyLocal)==0) then
		
		myMesh(cell)%neighborProcs(1,3)=myProc%procNeighbor(3)
		if(myMesh(cell)%neighborProcs(1,3)==myProc%taskid) then
			myMesh(cell)%neighbors(1,3)=cell-(numxLocal*(numyLocal-1))
		else
			numSend(3)=numSend(3)+1
			sendBuffer(1,numSend(3),3)=cell
			sendBuffer(2,numSend(3),3)=myMesh(cell)%material

		end if
	else
		myMesh(cell)%neighbors(1,3)=cell+numxLocal
		myMesh(cell)%neighborProcs(1,3)=myProc%taskid
	end if
	
	
	if(mod(cell,numxLocal*numyLocal) <= numxLocal .AND. (mod(cell, numxLocal*numyLocal) /= 0 .OR. numyLocal==1)) then
		myMesh(cell)%neighborProcs(1,4)=myProc%procNeighbor(4)
		if(myMesh(cell)%neighborProcs(1,4)==myProc%taskid) then
			myMesh(cell)%neighbors(1,4)=cell+(numxLocal*(numyLocal-1))
		else
			numSend(4)=numSend(4)+1
			sendBuffer(1,numSend(4),4)=cell
			sendBuffer(2,numSend(4),4)=myMesh(cell)%material

		end if
	else
		myMesh(cell)%neighbors(1,4)=cell-numxLocal
		myMesh(cell)%neighborProcs(1,4)=myProc%taskid
	end if
	
	if(mod(cell,numxLocal*numyLocal*numzLocal) > numxLocal*numyLocal*(numzLocal-1) .OR. &
			mod(cell, numxLocal*numyLocal*numzLocal)==0) then

		myMesh(cell)%neighborProcs(1,5)=myProc%procNeighbor(5)
		globalCell=myMesh(cell)%globalCell
		globalNeighbor=findgNeighborFreeSurf(globalCell, 5)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(1,5)=0
			myMesh(cell)%neighborProcs(1,5)=-1
		else
			if(myMesh(cell)%neighborProcs(1,5)==myProc%taskid) then
				myMesh(cell)%neighbors(1,5)=cell-(numxLocal*numyLocal*(numzLocal-1))
			else
				numSend(5)=numSend(5)+1
				sendBuffer(1,numSend(5),5)=cell
				sendBuffer(2,numSend(5),5)=myMesh(cell)%material

			end if
		end if
	else
		myMesh(cell)%neighbors(1,5)=cell+numxLocal*numyLocal
		myMesh(cell)%neighborProcs(1,5)=myProc%taskid
	end if
	
	if(mod(cell,numxLocal*numyLocal*numzLocal) <= numxLocal*numyLocal .AND. &
			(mod(cell,numxLocal*numyLocal*numzLocal) /= 0 .OR. numzLocal==1)) then

		myMesh(cell)%neighborProcs(1,6)=myProc%procNeighbor(6)
		globalCell=myMesh(cell)%globalCell
		globalNeighbor=findgNeighborFreeSurf(globalCell, 6)
		if(globalNeighbor==0) then
			!free surface, set proc id to -1 and cell id to 0 to indicate free surface
			myMesh(cell)%neighbors(1,6)=0
			myMesh(cell)%neighborProcs(1,6)=-1
		else
			if(myMesh(cell)%neighborProcs(1,6)==myProc%taskid) then
				myMesh(cell)%neighbors(1,6)=cell+(numxLocal*numyLocal*(numzLocal-1))
			else
				numSend(6)=numSend(6)+1
				sendBuffer(1,numSend(6),6)=cell
				sendBuffer(2,numSend(6),6)=myMesh(cell)%material

			end if
		end if
	else
		myMesh(cell)%neighbors(1,6)=cell-numxLocal*numyLocal
		myMesh(cell)%neighborProcs(1,6)=myProc%taskid
	end if
end do

!***************************************************************
!Send
do dir=1,6
	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		call MPI_ISEND(sendBuffer(:,1:numSend(dir),dir), numSend(dir)*2, MPI_INTEGER, myProc%procNeighbor(dir), &
				200+dir, comm, sendRequest(dir),ierr)
	end if
end do

!************************************************************
!Recv
do dir=1,6

	if(mod(dir,2)==0) then
		tag = dir-1
	else
		tag = dir+1
	end if

	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		numRecv(dir)=numSend(dir)
		allocate(recvBuffer(2,numRecv(dir)))
		call MPI_IRECV(recvBuffer, numRecv(dir)*2, MPI_INTEGER, myProc%procNeighbor(dir), &
				200+tag, comm, recvRequest(dir), ierr)
		call MPI_WAIT(recvRequest(dir),recvStatus(:,dir),ierr)

		do i=1, numRecv(dir)
			localCell=sendBuffer(1,i,dir)
			myMesh(localCell)%neighbors(1,dir)=recvBuffer(1,i)
			materialBuff(localCell,dir)=recvBuffer(2,i)

		end do
		deallocate(recvBuffer)
	end if
end do

do dir=1,6
	if(myProc%procNeighbor(dir)/=myProc%taskid) then
		call MPI_WAIT(sendRequest(dir),sendStatus(:,dir),ierr)
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
		if(myMesh(i)%neighborProcs(1,dir) /= myProc%taskid) then	!we are pointing to a different proc
			if(myMesh(i)%neighbors(1,dir) > maxElement) then		!searching for the max element number in a neighbor
				maxElement=myMesh(i)%neighbors(1,dir)
			end if
		end if
	end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(maxElement,6))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do i=1,maxElement
	do dir=1,6
		myBoundary(i,dir)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
		myBoundary(i,dir)%proc=-10		!default, says that this is not a real element of myBoundary.
		myBoundary(i,dir)%material=0
		myBoundary(i,dir)%length=0d0
		myBoundary(i,dir)%volume=0d0
	end do
end do

do i=1,numCells
	do dir=1,6
		if(myMesh(i)%neighborProcs(1,dir) == -1) then										!this is a free surface
			!do nothing
		else if(myMesh(i)%neighborProcs(1,dir) /= myProc%taskid .AND. myMesh(i)%neighborProcs(1,dir) /= -1) then
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%proc=myMesh(i)%neighborProcs(1,dir)	!set proc # of elements in myBoundary
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%length=length						!set length of elements in myBoundary
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%volume=length**3d0					!set volume of elements in myBoundary (changes with cascade addition)
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%material=materialBuff(i,dir)
			myBoundary(myMesh(i)%neighbors(1,dir),dir)%localNeighbor=i

!           if(strainField=='yes') then
!                do l=1,6
!                   myBoundary(myMesh(i)%neighbors(1,dir),dir)%strain(l)=globalStrain(globalNeighbor,l)
!                    myBoundary(myMesh(i)%neighbors(1,dir),dir)%strain(l)=materialStrain(1+l,i,dir)
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

integer status(MPI_STATUS_SIZE), i, j, k, matNum
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

	subroutine createProcCoordList(procCoordList)
		implicit none
		double precision, allocatable :: procCoordList(:,:)
	end subroutine

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


!step 1: divide volume among processors
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

!This contains a list of ALL processor coordinate max/mins, needed later
allocate(procCoordList(6,myProc%numtasks+1))
call createProcCoordList(procCoordList)

!write(*,*) 'proc', myProc%taskid, 'of', myProc%numtasks
!write(*,*) 'local coordinates', myProc%localCoord

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

allocate(globalMeshCoord(3,numTotal))		!coordinates of the center of each element
allocate(globalMaterial(numTotal))			!material type of each element
allocate(globalNumNeighbors(6,numTotal))	!number of neighbors of each element in each direction
											!(directions: +x, -x, +y, -y, +z, -z)
allocate(globalLength(numTotal))			!volume element length

!Count the number of elements that are inside the local processor's bounds
localElem=0
globalMaxNeighbors=0
do elem=1,numTotal
	!read in element coordinates and material number as well as how many neighbors it has in each direction
	read(80,*) (globalMeshCoord(i,elem),i=1,3), globalLength(elem), globalMaterial(elem), (globalNumNeighbors(i,elem), i=1,6)

	!GlobalMaxNeighbors is used to create the global connectivity matrix. It is the max number of
	!neighbors in one direction that any one element has (if elements are of different sizes, one element
	!may have more than one neighbor in a given direction).

	do i=1,6
		if(globalNumNeighbors(i,elem) > globalMaxNeighbors) then
			globalMaxNeighbors=globalNumNeighbors(i,elem)
		end if
	end do

	!Count the number of elements that are inside the local coordinates.
	if(globalMeshCoord(1,elem) > myProc%localCoord(1) .AND. globalMeshCoord(1,elem) <= myProc%localCoord(2) &
		.AND. globalMeshCoord(2,elem) > myProc%localCoord(3) .AND. globalMeshCoord(2,elem) <= myProc%localCoord(4) &
		.AND. globalMeshCoord(3,elem) > myProc%localCoord(5) .AND. globalMeshCoord(3,elem) <= myProc%localCoord(6)) then

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
	if(globalMeshCoord(1,elem) > myProc%localCoord(1) .AND. globalMeshCoord(1,elem) <= myProc%localCoord(2) &
		.AND. globalMeshCoord(2,elem) > myProc%localCoord(3) .AND. globalMeshCoord(2,elem) <= myProc%localCoord(4) &
		.AND. globalMeshCoord(3,elem) > myProc%localCoord(5) .AND. globalMeshCoord(3,elem) <= myProc%localCoord(6)) then

		localElem=localElem+1	!count local elements

		myMesh(localElem)%coordinates(1)=globalMeshCoord(1,elem)	!read in coordinates, material, length, etc.
		myMesh(localElem)%coordinates(2)=globalMeshCoord(2,elem)
		myMesh(localElem)%coordinates(3)=globalMeshCoord(3,elem)
		myMesh(localElem)%material=globalMaterial(elem)
		myMesh(localElem)%length=globalLength(elem)
		myMesh(localElem)%volume=globalLength(elem)**3d0
		maxNumNeighbors=0
		!Find the max number of neighbors that this element has in any one direction
		do i=1,6
			myMesh(localElem)%numNeighbors(i)=globalNumNeighbors(i,elem)
			if(myMesh(localElem)%numNeighbors(i) > maxNumNeighbors) then
				maxNumNeighbors=myMesh(localElem)%numNeighbors(i)
			end if
		end do
		myMesh(localElem)%proc=myProc%taskid

		!The number of neighbors allowed in all directions is equal to maxNumNeighbors even if
		!the element may have fewer neighbors in some directions. These extra array elements are
		!ignored for the rest of the program.

		allocate(myMesh(localElem)%neighbors(maxNumNeighbors,6))
		allocate(myMesh(localElem)%neighborProcs(maxNumNeighbors,6))

	end if
end do

!Step 4: create global and local connectivity matrixes including multiple neighbors to the same direction

!Step 4a: create global connectivity matrix
allocate(globalMeshConnect(globalMaxNeighbors,6,numTotal))

if(meshType=='periodic') then
	call createConnectGlobalPeriodicNonUniform(globalMeshConnect, globalMeshCoord, globalNumNeighbors, globalLength, &
		 myProc%globalCoord)
else if(meshType=='freeSurfaces') then
	call createConnectGlobalFreeSurfNonUniform(globalMeshConnect, globalMeshCoord, globalNumNeighbors, globalLength, &
		 myProc%globalCoord)
end if


!Step 4b: create local connectivity myMesh(:)%neighbors(direction,num) and myMesh(:)%neighborProcs(direction,num)

call createConnectLocalNonUniform(globalMeshConnect, globalMeshCoord, localElem, procCoordList)

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
			if(myMesh(i)%neighborProcs(k,j) /= myProc%taskid) then	!we are pointing to a different proc
				if(myMesh(i)%neighbors(k,j) > maxElement) then		!searching for the max element number in a neighbor
					maxElement=myMesh(i)%neighbors(k,j)
				endif
			endif
		end do
	end do
end do

!This tells us how large to allocate myBoundary. NOTE: many elements in myBoundary will be unused. Only
!the elements that represent volume elements on the boundary of myMesh will be used. This represents
!wasted memory for the sake of easier computation

allocate(myBoundary(maxElement,6))	!6 directions, maxElement elements in each direction (more than needed)

!initialize myBoundary with 0 in localNeighbor - signal that myBoundary is not attached to anything
do i=1,maxElement
	do j=1,6
		myBoundary(i,j)%localNeighbor=0	!default, says that this is not a real element of myBoundary.
	end do
end do

!Step 2: initialize myBoundary elements (only the relevant ones) with processor #'s, material #'s, length
do i=1,localElem
	do j=1,6
		do k=1,myMesh(i)%numNeighbors(j)
			if(myMesh(i)%neighborProcs(k,j) /= myProc%taskid) then
				myBoundary(myMesh(i)%neighbors(k,j),j)%proc=myMesh(i)%neighborProcs(k,j)	!set proc # of elements in myBoundary
				globalCell=findGlobalCell(myMesh(i)%coordinates,globalMeshCoord)			!find global cell # of element in myBoundary
				globalNeighbor=globalMeshConnect(globalCell,j,k)								!use global cell # to find material # of element in myBoundary
				myBoundary(myMesh(i)%neighbors(1,j),j)%material=globalMaterial(globalNeighbor)	!set material # of elements in myBoundary
				myBoundary(myMesh(i)%neighbors(k,j),j)%length=globalLength(globalNeighbor)		!set length of elements in myBoundary
				myBoundary(myMesh(i)%neighbors(k,j),j)%volume=globalLength(globalNeighbor)**3d0	!set volume of elements in myBoundary (changes with cascade addition)
				myBoundary(myMesh(i)%neighbors(1,j),j)%localNeighbor=i
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

				if(Coord(1,cell)+length(cell)/2d0==bdryCoord(2)) then
					!periodic
					!We must search through all elements to see if they are the neighbors of other elements. This is because we
					!don't know how many neighbors each element has. This could be made more efficient in the future if it becomes
					!too computationally expensive, but for now it only has to be carried out once so it should be fine.

					if(Coord(1,searchElem)==Coord(1,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1						!count the number of neighbors
						Connect(neighbor,dir,cell)=searchElem	!Add this cell to the connectivity
					end if
				else if(Coord(1,searchElem)==Coord(1,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==2) then
				if(Coord(1,cell)-length(cell)/2d0==bdryCoord(1)) then
					!periodic
					if(Coord(1,searchElem)==Coord(1,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(1,searchElem)==Coord(1,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==3) then
				if(Coord(2,cell)+length(cell)/2d0==bdryCoord(4)) then
					!periodic
					if(Coord(2,searchElem)==Coord(2,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(2,searchElem)==Coord(2,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==4) then
				if(Coord(2,cell)-length(cell)/2d0==bdryCoord(3)) then
					!periodic
					if(Coord(2,searchElem)==Coord(2,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(2,searchElem)==Coord(2,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==5) then
				if(Coord(3,cell)+length(cell)/2d0==bdryCoord(6)) then
					!periodic
					if(Coord(3,searchElem)==Coord(3,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(3,searchElem)==Coord(3,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==6) then
				if(Coord(3,cell)-length(cell)/2d0==bdryCoord(5)) then
					!periodic
					if(Coord(3,searchElem)==Coord(3,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(3,searchElem)==Coord(3,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			end if
		end do
		if(NumNeighbors(dir,cell) /= neighbor) then
			!we have found the wrong number of neighbors
			write(*,*) 'error incorrect number of neighbors found', neighbor, 'neighbors', numNeighbors(dir,cell), 'expected'
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
				if(Coord(1,cell)+length(cell)/2d0==bdryCoord(2)) then
					!periodic
					if(Coord(1,searchElem)==Coord(1,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(1,searchElem)==Coord(1,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==2) then
				if(Coord(1,cell)-length(cell)/2d0==bdryCoord(1)) then
					!periodic
					if(Coord(1,searchElem)==Coord(1,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(2)-bdryCoord(1)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(1,searchElem)==Coord(1,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==3) then
				if(Coord(2,cell)+length(cell)/2d0==bdryCoord(4)) then
					!periodic
					if(Coord(2,searchElem)==Coord(2,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(2,searchElem)==Coord(2,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==4) then
				if(Coord(2,cell)-length(cell)/2d0==bdryCoord(3)) then
					!periodic
					if(Coord(2,searchElem)==Coord(2,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(4)-bdryCoord(3)) &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=searchElem
					end if
				else if(Coord(2,searchElem)==Coord(2,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(3,searchElem)-Coord(3,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==5) then
				if(Coord(3,cell)+length(cell)/2d0==bdryCoord(6)) then
					!periodic
					if(Coord(3,searchElem)==Coord(3,cell)+(length(cell)+length(searchElem))/2d0-(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=0
					end if
				else if(Coord(3,searchElem)==Coord(3,cell)+(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			else if(dir==6) then
				if(Coord(3,cell)-length(cell)/2d0==bdryCoord(5)) then
					!periodic
					if(Coord(3,searchElem)==Coord(3,cell)-(length(cell)+length(searchElem))/2d0+(bdryCoord(6)-bdryCoord(5)) &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then

						neighbor=neighbor+1
						Connect(neighbor,dir,cell)=0
					end if
				else if(Coord(3,searchElem)==Coord(3,cell)-(length(cell)+length(searchElem))/2d0 &
					.AND. dabs(Coord(2,searchElem)-Coord(2,cell)) <= dmax1(length(cell),length(searchElem))/2d0 &
					.AND. dabs(Coord(1,searchElem)-Coord(1,cell)) <= dmax1(length(cell),length(searchElem))/2d0) then
					!we have found a partner
					neighbor=neighbor+1
					Connect(neighbor,dir,cell)=searchElem
				end if
			end if
		end do
		if(NumNeighbors(dir,cell) /= neighbor) then
			!we have found the wrong number of neighbors
			write(*,*) 'error incorrect number of neighbors found', neighbor, 'neighbors', numNeighbors(dir,cell), 'expected'
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
integer sendList(5,6*localElem), i, status(MPI_STATUS_SIZE), neighborProc

!For each local element, find its global cell, find its neighbors (global) in the globalNeighbor variable,
!and if globalNeighbor is not in the local mesh, then send/recieve with the neighboring processor
!in the buffer system (same as uniform mesh)

numSendRecv=0
do elem=1,localElem
	globalCell=findGlobalCell(myMesh(elem)%coordinates,globalMeshCoord)
	do dir=1,6
		do neighbor=1,myMesh(elem)%numNeighbors(dir)
			globalNeighbor=globalMeshConnect(neighbor,dir,globalCell)
			if(globalNeighbor==0) then
				!free surface, set cell number to 0 and proc number to -1
				myMesh(elem)%neighbors(neighbor,dir)=0
				myMesh(elem)%neighborProcs(neighbor,dir)=-1
			else
				!find out if globalNeighbor is in the local mesh. If so, set myMesh(cell)%neighbors and procs
				if(globalMeshCoord(1,globalNeighbor) > myProc%localCoord(1) &
				.AND. globalMeshCoord(1,globalNeighbor) <= myProc%localCoord(2) &
				.AND. globalMeshCoord(2,globalNeighbor) > myProc%localCoord(3) &
				.AND. globalMeshCoord(2,globalNeighbor) <= myProc%localCoord(4) &
				.AND. globalMeshCoord(3,globalNeighbor) > myProc%localCoord(5) &
				.AND. globalMeshCoord(3,globalNeighbor) <= myProc%localCoord(6)) then

					localNeighbor=findLocalCell(globalMeshCoord(:,globalNeighbor))	!find cell number of neighbor in local mesh
					myMesh(elem)%neighbors(neighbor,dir)=localNeighbor
					myMesh(elem)%neighborProcs(neighbor,dir)=myProc%taskid
				!if no, send local mesh number and recieve local mesh number from neigboring processor
				else
					!find the number of the neighboring processor based on the coordinates of the neighboring cell
					neighborProc=findNeighborProc(globalMeshCoord,procCoordList,globalNeighbor)
					myMesh(elem)%neighborProcs(neighbor,dir)=neighborProc

					numSendRecv=numSendRecv+1

					!Create buffer of information to send/recieve between processors
					sendList(1,numSendRecv)=elem
					sendList(2,numSendRecv)=dir	!direction
					sendList(3,numSendRecv)=globalCell
					sendList(4,numSendRecv)=globalNeighbor
					sendList(5,numSendRecv)=neighbor	!which neighbor number is this
				end if
			end if
		end do
	end do
end do

!Send/recieve the element in the local mesh / the neighbor in the neighboring processor's mesh (respectively).
!MPI_SEND and MPI_RECV functions are matched using the cell number of the sender

do i=1,numSendRecv

	elem=sendList(1,i)
	dir=sendList(2,i)
	globalCell=sendList(3,i)
	neighbor=sendList(5,i)
	call MPI_SEND(elem, 1, MPI_INTEGER, myMesh(elem)%neighborProcs(neighbor,dir), globalCell, comm, ierr)
end do

do i=1,numSendRecv
	elem=sendList(1,i)
	dir=sendList(2,i)
	globalNeighbor=sendList(4,i)
	neighbor=sendList(5,i)
	call MPI_RECV(myMesh(elem)%neighbors(neighbor,dir), 1, MPI_INTEGER, myMesh(elem)%neighborProcs(neighbor,dir), &
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
		if(coord(1)==gCoord(1,i) .AND. coord(2)==gCoord(2,i) .AND. coord(3)==gCoord(3,i)) then
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

subroutine createProcCoordList(procCoordList)
use DerivedType
use mod_constants
implicit none

integer i
double precision, allocatable :: procCoordList(:,:)

do i=1,myProc%numtasks
	procCoordList(1,i)=myProc%globalCoord(1)+&
		mod(i-1,dims(1))*(myProc%globalCoord(2)-myProc%globalCoord(1))/&
		dble(dims(1))
	
	procCoordList(2,i)=procCoordList(1,i)+&
		(myProc%globalCoord(2)-myProc%globalCoord(1))/dble(dims(1))
	
	procCoordList(3,i)=myProc%globalCoord(3)+&
		mod((i-1)/dims(1),dims(2))*&
		(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))
	
	procCoordList(4,i)=procCoordList(3,i)+&
		(myProc%globalCoord(4)-myProc%globalCoord(3))/dble(dims(2))
	
	procCoordList(5,i)=myProc%GlobalCoord(5)+&
		mod((i-1)/(dims(1)*dims(2)),dims(3))*&
		(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))
	
	procCoordList(6,i)=procCoordList(5,i)+&
		(myProc%globalCoord(6)-myProc%globalCoord(5))/dble(dims(3))
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
	if(globalMeshCoord(1,elem) > procCoordList(i,1) .AND. globalMeshCoord(1,elem) <= procCoordList(i,2) &
	.AND. globalMeshCoord(2,elem) > procCoordList(i,3) .AND. globalMeshCoord(2,elem) <= procCoordList(i,4) &
	.AND. globalMeshCoord(3,elem) > procCoordList(i,5) .AND. globalMeshCoord(3,elem) <= procCoordList(i,6)) then
		findNeighborProc=i-1
		exit
	endif
end do

end function


end module

