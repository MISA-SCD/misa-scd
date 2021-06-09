!***************************************************************************************************
!>Subroutine initialMeshUniform():
!The subroutine creates meshes and sectors, and establishes the mapping relationship among meshes, sectors and processes.
!Each process is divided into eight sectors.
!NOTE: Each sector should have at least two meshes per dimension (x, y, z dimension).
!***************************************************************************************************
subroutine initialMesh()
	use mod_constants
	use mod_globalVariables
	use mod_structures

	implicit none
	include 'mpif.h'

	integer :: status(MPI_STATUS_SIZE), i, j, k, dir, localElem
	double precision :: length
	integer :: x, y, z, numXmin,numXmax,numYmin,numYmax,numZmin,numZmax	!used to determine numxLocal, numyLocal, numzLocal
	double precision :: tempCenter(6)

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
		if((dble(numXmin)*(length/2d0)) <=  myProc%localCoord(1) .AND. myProc%localCoord(1) < &
				(dble(numXmin+1)*(length/2d0))) then
			if(mod(numXmin,2)==0) then
				numXmin = numXmin/2
			else
				numXmin = (numXmin+1)/2
			end if
		end if
		if((dble(numXmax)*(length/2d0)) <=  myProc%localCoord(2) .AND. myProc%localCoord(2) < &
				(dble(numXmax+1)*(length/2d0))) then
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
		if((dble(numYmin)*(length/2d0)) <=  myProc%localCoord(3) .AND. myProc%localCoord(3) < &
				(dble(numYmin+1)*(length/2d0))) then
			if(mod(numYmin,2)==0) then
				numYmin = numYmin/2
			else
				numYmin = (numYmin+1)/2
			end if
		end if
		if((dble(numYmax)*(length/2d0)) <=  myProc%localCoord(4) .AND. myProc%localCoord(4) < &
				(dble(numYmax+1)*(length/2d0))) then
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
		if((dble(numZmin)*(length/2d0)) <=  myProc%localCoord(5) .AND. myProc%localCoord(5) < &
				(dble(numZmin+1)*(length/2d0))) then
			if(mod(numZmin,2)==0) then
				numZmin = numZmin/2
			else
				numZmin = (numZmin+1)/2
			end if
		end if
		if((dble(numZmax)*(length/2d0)) <=  myProc%localCoord(6) .AND. myProc%localCoord(6) < &
				(dble(numZmax+1)*(length/2d0))) then
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
	if(myProc%taskid==MASTER) then
		write(*,*) 'numx numy numz', numx, numy, numz
	end if
	!write(*,*) 'proc', myProc%taskid, 'numxLocal numyLocal numzLocal', numxLocal, numyLocal, numzLocal

	!If any processors don't have volume elements in them (too many procs), we create an error message
	if(numCells==0) then
		write(*,*) 'error processors with no volume elements'
		call MPI_ABORT(comm,ierr)
	end if

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

				myMesh(localElem)%numNeighbors=1
				myMesh(localElem)%neighbors=0
				myMesh(localElem)%neighborProcs=-1
				!allocate(myMesh(localElem)%neighbors(6))
				!allocate(myMesh(localElem)%neighborProcs(6))
				!do dir=1,6
				!	myMesh(localElem)%numNeighbors(dir)=1
				!end do
			end do
		end do
	end do

	!Step 3d: assign neighbors and processor numbers for neighbors (connectivity in myMesh) - periodic or free surfaces
	call createMeshConnect(length)

end subroutine

!**************************************************************************************************
!>Subroutine createMeshConnect(length) (periodic boundary conditions or free surface)
!It identifies the mesh and processor of neighboring meshes for each mesh.
!The connectivity scheme is the same as in the global case, but neighboring processor numbers are used here.
!**************************************************************************************************
subroutine createMeshConnect(length)
	use mod_structures
	use mod_globalVariables
	implicit none
	include 'mpif.h'

	double precision, intent(in) :: length
	integer :: cell, localCell, maxElement, maxBndryCells
	integer :: i, dir, tag
	integer, allocatable :: send(:,:,:), recv(:,:)
	integer :: materialBuff(numCells,6)
	integer :: numSend(6)=0, numRecv(6)=0
	integer :: status(MPI_STATUS_SIZE)
	integer :: tempRecv
	double precision :: commTime1, commTime2

	maxBndryCells = numxLocal*numyLocal
	if((numyLocal*numzLocal) > maxBndryCells) then
		maxBndryCells = numyLocal*numzLocal
	end if
	if((numxLocal*numzLocal) > maxBndryCells) then
		maxBndryCells = numxLocal*numzLocal
	end if

	!allocate(send(2, numxLocal*numyLocal*numzLocal,6))
	allocate(send(2, maxBndryCells,6))

	do cell=1,numCells
		!Right (+x)
		if(mod(cell,numxLocal)==0) then !identify cell to the right
			myMesh(cell)%neighborProcs(1)=myProc%procNeighbor(1)
			if(myProc%procNeighbor(1)==-1) then	!+x free surface
				myMesh(cell)%neighbors(1)=0
			else if(myMesh(cell)%neighborProcs(1)==myProc%taskid) then	!+x periodic
				myMesh(cell)%neighbors(1)=cell-numxLocal+1
			else
				numSend(1)=numSend(1)+1			!count the number of items in buffer
				send(1,numSend(1),1)=cell		!localCellID in this processor
				send(2,numSend(1),1)=myMesh(cell)%material	!Material ID number that this element is composed of
			end if
		else
			myMesh(cell)%neighbors(1)=cell+1
			myMesh(cell)%neighborProcs(1)=myProc%taskid
		end if

		!Left (-x)
		if(mod(cell+numxLocal-1,numxLocal)==0) then !identify cell to the left
			myMesh(cell)%neighborProcs(2)=myProc%procNeighbor(2)
			if(myProc%procNeighbor(2)==-1) then	!-x free surface
				myMesh(cell)%neighbors(2)=0
			else if(myMesh(cell)%neighborProcs(2)==myProc%taskid) then
				myMesh(cell)%neighbors(2)=cell+numxLocal-1
			else
				numSend(2)=numSend(2)+1
				send(1,numSend(2),2)=cell
				send(2,numSend(2),2)=myMesh(cell)%material
			end if
		else
			myMesh(cell)%neighbors(2)=cell-1
			myMesh(cell)%neighborProcs(2)=myProc%taskid
		end if

		!Front (+y)
		if(mod(cell,numxLocal*numyLocal) > numxLocal*(numyLocal-1) .OR. mod(cell,numxLocal*numyLocal)==0) then
			myMesh(cell)%neighborProcs(3)=myProc%procNeighbor(3)
			if(myProc%procNeighbor(3)==-1) then
				myMesh(cell)%neighbors(3)=0
			else if(myMesh(cell)%neighborProcs(3)==myProc%taskid) then
				myMesh(cell)%neighbors(3)=cell-(numxLocal*(numyLocal-1))
			else
				numSend(3)=numSend(3)+1
				send(1,numSend(3),3)=cell
				send(2,numSend(3),3)=myMesh(cell)%material
			end if
		else
			myMesh(cell)%neighbors(3)=cell+numxLocal
			myMesh(cell)%neighborProcs(3)=myProc%taskid
		end if

		!Back (-y)
		if(mod(cell,numxLocal*numyLocal) <= numxLocal .AND. (mod(cell, numxLocal*numyLocal) /= 0 &
				.OR. numyLocal==1)) then
			myMesh(cell)%neighborProcs(4)=myProc%procNeighbor(4)
			if(myProc%procNeighbor(4)==-1) then
				myMesh(cell)%neighbors(4)=0
			else if(myMesh(cell)%neighborProcs(4)==myProc%taskid) then
				myMesh(cell)%neighbors(4)=cell+(numxLocal*(numyLocal-1))
			else
				numSend(4)=numSend(4)+1
				send(1,numSend(4),4)=cell
				send(2,numSend(4),4)=myMesh(cell)%material
			end if

		else
			myMesh(cell)%neighbors(4)=cell-numxLocal
			myMesh(cell)%neighborProcs(4)=myProc%taskid
		end if

		!Up (+z)
		if(mod(cell,numxLocal*numyLocal*numzLocal) > numxLocal*numyLocal*(numzLocal-1) &
				.OR. mod(cell, numxLocal*numyLocal*numzLocal)==0) then
			myMesh(cell)%neighborProcs(5)=myProc%procNeighbor(5)
			if(myProc%procNeighbor(5)==-1) then
				myMesh(cell)%neighbors(5)=0
			else if(myMesh(cell)%neighborProcs(5)==myProc%taskid) then
				myMesh(cell)%neighbors(5)=cell-(numxLocal*numyLocal*(numzLocal-1))
			else
				numSend(5)=numSend(5)+1
				send(1,numSend(5),5)=cell
				send(2,numSend(5),5)=myMesh(cell)%material
			end if

		else
			myMesh(cell)%neighbors(5)=cell+numxLocal*numyLocal
			myMesh(cell)%neighborProcs(5)=myProc%taskid
		end if

		!Down (-z)
		if(mod(cell,numxLocal*numyLocal*numzLocal) <= numxLocal*numyLocal &
				.AND. (mod(cell,numxLocal*numyLocal*numzLocal) /= 0 .OR. numzLocal==1)) then
			myMesh(cell)%neighborProcs(6)=myProc%procNeighbor(6)
			if(myProc%procNeighbor(6)==-1) then
				myMesh(cell)%neighbors(6)=0
			else if(myMesh(cell)%neighborProcs(6)==myProc%taskid) then
				myMesh(cell)%neighbors(6)=cell+(numxLocal*numyLocal*(numzLocal-1))
			else
				numSend(6)=numSend(6)+1
				send(1,numSend(6),6)=cell
				send(2,numSend(6),6)=myMesh(cell)%material
			end if

		else
			myMesh(cell)%neighbors(6)=cell-numxLocal*numyLocal
			myMesh(cell)%neighborProcs(6)=myProc%taskid
		end if
	end do

	maxElement=0
	!*******************************************************
	do dir = 1,6
		if(dir == 1 .OR. dir == 2) then			!<x
			numRecv(dir) = numyLocal*numzLocal
		else if(dir == 3 .OR. dir == 4) then	!<y
			numRecv(dir) = numxLocal*numzLocal
		else if(dir == 5 .OR. dir == 6) then	!<z
			numRecv(dir) = numxLocal*numyLocal
		end if
	end do

	do dir = 1,6

		if(mod(dir,2)==0) then
			tag = dir-1
		else
			tag = dir+1
		end if


		if(myProc%procNeighbor(dir) /= myProc%taskid .AND. myProc%procNeighbor(dir) /= -1) then
			commTime1=MPI_WTIME()
			call MPI_SEND(send(1,1,dir),numSend(dir)*2,MPI_INTEGER,myProc%procNeighbor(dir),dir,comm,ierr)
			commTime2=MPI_WTIME()
			commTimeSum=commTime2-commTime1
		end if

		if(myProc%procNeighbor(tag)/=myProc%taskid .AND. myProc%procNeighbor(tag) /= -1) then
			allocate(recv(2,numRecv(tag)))
			commTime1=MPI_WTIME()
			call MPI_RECV(recv,numRecv(tag)*2,MPI_INTEGER,myProc%procNeighbor(tag),dir,comm,status,ierr)
			commTime2=MPI_WTIME()
			commTimeSum=commTime2-commTime1

			do i=1, numRecv(tag)
				localCell=send(1,i,tag)
				myMesh(localCell)%neighbors(tag)=recv(1,i)
				materialBuff(localCell,tag)=recv(2,i)
			end do

			!if(myProc%taskid == 0) then
			!	do i = 1,numRecv(tag)
			!		write(*,*) 'recvDir', tag, 'numRecv', numRecv(tag), 'recvCell', recv(1,i)
			!	end do
			!end if

			deallocate(recv)
		end if
	end do
	deallocate(send)

	!***************************************************************************************************
	!Initializing myGhost with elements that are in neighboring processors that bound this one
	!***************************************************************************************************
	maxElement=0
	do i=1, numCells
		do dir=1, 6
			if(myMesh(i)%neighborProcs(dir) /= myProc%taskid) then	!we are pointing to a different proc
				if(myMesh(i)%neighbors(dir) > maxElement) then		!searching for the max element number in a neighbor
					maxElement=myMesh(i)%neighbors(dir)
				end if
			end if
		end do
	end do

	allocate(myGhost(maxElement,6))	!6 directions, maxElement elements in each direction (more than needed)
	!initialize myGhost with 0 in localNeighbor - signal that myGhost is not attached to anything
	do dir=1,6
		do i=1,maxElement
			myGhost(i,dir)%localNeighbor=0	!default, says that this is not a real element of myGhost.
			myGhost(i,dir)%proc=-10			!default, says that this is not a real element of myGhost.
			myGhost(i,dir)%material=0
			myGhost(i,dir)%length=0d0
			myGhost(i,dir)%volume=0d0
		end do
	end do

	do i=1,numCells
		do dir=1,6
			if(myMesh(i)%neighborProcs(dir) == -1) then										!this is a free surface
				!do nothing
			else if(myMesh(i)%neighborProcs(dir) /= myProc%taskid) then
				myGhost(myMesh(i)%neighbors(dir),dir)%proc=myMesh(i)%neighborProcs(dir)	!set proc # of elements in myGhost
				myGhost(myMesh(i)%neighbors(dir),dir)%length=length						!set length of elements in myGhost
				myGhost(myMesh(i)%neighbors(dir),dir)%volume=length**3d0					!set volume of elements in myGhost (changes with cascade addition)
				myGhost(myMesh(i)%neighbors(dir),dir)%material=materialBuff(i,dir)
				myGhost(myMesh(i)%neighbors(dir),dir)%localNeighbor=i
			end if
		end do
	end do

end subroutine

!***************************************************************************************
!>Function findgNeighborGID(globalID, dir)
!Inputs: globalID, direction
!Output: globalNeighborID of the globalID in the direaction
!***************************************************************************************
integer function findgNeighborGID(globalID, dir)
	use mod_structures
	use mod_globalVariables
	implicit none

	integer, intent(in) :: globalID, dir
	integer :: neighborID

	if(dir==1) then	!+x
		if(mod(globalID,numx)==0) then !identify cell to the right
			if(periods(1) .eqv. .true.) then	!periodic
				neighborID=globalID-numx+1
			else	!free surface
				neighborID=0
			end if
		else
			neighborID=globalID+1
		end if
	else if(dir==2) then
		if(mod(globalID+numx-1,numx)==0) then !identify cell to the left
			if(periods(1) .eqv. .true.) then
				neighborID=globalID+numx-1
			else
				neighborID=0
			end if
		else
			neighborID=globalID-1
		end if
	else if(dir==3) then
		if(mod(globalID,numx*numy) > numx*(numy-1) .OR. mod(globalID,numx*numy)==0) then	!identify cell to the front
			if(periods(2) .eqv. .true.) then
				neighborID=globalID-(numx*(numy-1))
			else
				neighborID=0
			end if
		else
			neighborID=globalID+numx
		end if
	else if(dir==4) then
		if(mod(globalID,numx*numy) <= numx .AND. (mod(globalID, numx*numy) /= 0 .OR. numy==1)) then	!identify cell to the back
			if(periods(2) .eqv. .true.) then
				neighborID=globalID+(numx*(numy-1))
			else
				neighborID=0
			end if
		else
			neighborID=globalID-numx
		end if
	else if(dir==5) then
		if(mod(globalID,numx*numy*numz) > numx*numy*(numz-1) .OR. mod(globalID, numx*numy*numz)==0) then	!identify cell to up
			if(periods(3) .eqv. .true.) then
				neighborID=globalID-(numx*numy*(numz-1))
			else
				neighborID=0
			end if
		else
			neighborID=globalID+numx*numy
		end if
	else if(dir==6) then
		if(mod(globalID,numx*numy*numz) <= numx*numy .AND. (mod(globalID,numx*numy*numz) /= 0 .OR. numz==1)) then	!identify cell to down
			if(periods(3) .eqv. .true.) then
				neighborID=globalID+(numx*numy*(numz-1))
			else
				neighborID=0
			end if
		else
			neighborID=globalID-numx*numy
		end if
	end if

	findgNeighborGID=neighborID

end function


