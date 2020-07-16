!***************************************************************************************************
!>Subroutine initialMeshUniform():
!The subroutine creates meshes and sectors, and establishes the mapping relationship among meshes, sectors and processes.
!Each process is divided into eight sectors.
!NOTE: Each sector should have at least two meshes per dimension (x, y, z dimension).
!***************************************************************************************************
subroutine initialMesh()
use mod_constants
use DerivedType

implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE), i, j, k, dir
double precision length
integer localElem
integer x, y, z, numXmin,numXmax,numYmin,numYmax,numZmin,numZmax	!used to determine numxLocal, numyLocal, numzLocal

double precision tempCenter(6)

!open(80, file=filename,action='read', status='old')	!MeshGenInput_XX.txt
length=meshLength
numTotal=numx*numy*numz	!total cell in the system
systemVol = dble(numTotal)*length**3d0

!The  boundary coordinates of the system (xmin,xmax,ymin,ymax,zmin,zmax)
do i=1,5,2
	myProc%globalCoord(i)=0d0
end do
myProc%globalCoord(2)=length*dble(numx)
myProc%globalCoord(4)=length*dble(numy)
myProc%globalCoord(6)=length*dble(numz)

!The  boundary coordinates of the processor (xmin,xmax,ymin,ymax,zmin,zmax)
myProc%localCoord(1)=dble(myProc%coords(1))*myProc%globalCoord(2)/dble(dims(1))    !<xmin
myProc%localCoord(2)=myProc%localCoord(1)+myProc%globalCoord(2)/dble(dims(1))
myProc%localCoord(3)=dble(myProc%coords(2))*myProc%globalCoord(4)/dble(dims(2))
myProc%localCoord(4)=myProc%localCoord(3)+myProc%globalCoord(4)/dble(dims(2))
myProc%localCoord(5)=dble(myProc%coords(3))*myProc%globalCoord(6)/dble(dims(3))
myProc%localCoord(6)=myProc%localCoord(5)+myProc%globalCoord(6)/dble(dims(3))

totalVolume=(myProc%localCoord(2)-myProc%localCoord(1))*(myProc%localCoord(4)-myProc%localCoord(3))*&
        (myProc%localCoord(6)-myProc%localCoord(5))

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

!Step 5: Create meshes. And for each volume element in myMesh, assign coordinates and material number (and proc number)
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
		end do
	end do
end do

!Step 3d: assign neighbors and processor numbers for neighbors (connectivity in myMesh) - periodic or free surfaces in +/- z
if(meshType=='periodic') then
	call createConnectLocalPeriodic(length)
else if(meshType=='freeSurfaces') then
	call createConnectLocalFreeSurf(length)
end if

end subroutine

!**************************************************************************************************
!>Subroutine createConnectLocalPeriodic(length) (periodic boundary conditions)
!It identifies the mesh and processor of neighboring meshes for each mesh.
!The connectivity scheme is the same as in the global case, but neighboring processor numbers are used here.
!**************************************************************************************************
subroutine createConnectLocalPeriodic(length)
use DerivedType
use mod_constants
implicit none
include 'mpif.h'

integer cell, localCell, maxElement
double precision length
!buffer lists to send all information at the end
integer i, dir, tag

integer, allocatable :: send(:,:,:), recv(:,:)
integer, allocatable :: sendBuffer(:,:), recvBuffer(:,:)
integer materialBuff(numCells,6)
integer :: numSend(6)=0, numRecv(6)=0

integer sendRequest(6), recvRequest(6)
integer sendStatus(MPI_STATUS_SIZE,6), recvStatus(MPI_STATUS_SIZE,6)
integer status(MPI_STATUS_SIZE)

integer tempRecv
logical flagProbe

!************************************************
!periodic boundary condition version
!************************************************
allocate(send(2, numxLocal*numyLocal*numzLocal,6))
!allocate(recv(2, numxLocal*numyLocal*numzLocal,6))

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

		call MPI_SEND(send(1,1,dir),numSend(dir)*2,MPI_INTEGER,myProc%procNeighbor(dir),200+dir,&
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

		tempRecv=0
		call MPI_PROBE(myProc%procNeighbor(dir), 200+tag,comm,status,ierr)
		call MPI_GET_COUNT(status,MPI_INTEGER,tempRecv,ierr)
		numRecv(dir)=tempRecv/2

	!	if(flagProbe .eqv. .TRUE.) then
	!		call MPI_GET_COUNT(status,MPI_INTEGER,tempRecv)
	!		numRecv(dir)=tempRecv/2
	!	end if
		!numRecv(dir)=numSend(dir)
		allocate(recv(2,numRecv(dir)))
		call MPI_RECV(recv,numRecv(dir)*2,MPI_INTEGER,myProc%procNeighbor(dir),200+tag,&
				comm,status,ierr)
		do i=1, numRecv(dir)
			localCell=send(1,i,dir)
		!	myMesh(localCell)%neighbors(1,dir)=recvBuffer(1,i)
		!	materialBuff(localCell,dir)=recvBuffer(2,i)
			myMesh(localCell)%neighbors(1,dir)=recv(1,i)
			materialBuff(localCell,dir)=recv(2,i)
		end do
		deallocate(recv)
	end if
end do

!***************************************************************************************************
!Initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************
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

        end if
    end do
end do
deallocate(send)

end subroutine

!**************************************************************************************************
!>Subroutine createConnectLocalFreeSurf(length) (free surfaces in z-direction and periodic boundary conditions in other directions)
!This is the same as the previous subroutine except if the cell is at the free surface, its connectivity cell number
!is 0 and the processor number is -1.
!**************************************************************************************************
subroutine createConnectLocalFreeSurf(length)
use DerivedType
use mod_constants
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
		call MPI_ISEND(sendBuffer(1,1,dir), 2*numSend(dir), MPI_INTEGER, myProc%procNeighbor(dir), &
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
		call MPI_WAIT(recvRequest(dir),recvStatus(1,dir),ierr)

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
		call MPI_WAIT(sendRequest(dir),sendStatus(1,dir),ierr)
	end if
end do

deallocate(sendBuffer)

!***************************************************************************************************
!Initializing myBoundary with elements that are in neighboring processors that bound this one
!***************************************************************************************************
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

		end if
	end do
end do

end subroutine

!***************************************************************************************
!>Function findgNeighborPeriodic(globalID, dir)
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
!>Function findgNeighborFreeSurf(globalID, dir)
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

