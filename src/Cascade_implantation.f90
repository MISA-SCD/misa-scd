
!***************************************************************************************************
!>Subroutine: Choose Cascade
!!Takes list of cascades (read from input file) and chooses one randomly
!!
!!Inputs: CascadeList (global variable)
!!Output: CascadeTemp (pointing at the cascade we want)
!***************************************************************************************************

subroutine chooseCascade(CascadeTemp)
use DerivedType
use randdp
use mod_srscd_constants
implicit none

type(cascadeEvent), pointer :: cascadeTemp
double precision r, atemp

r=dprand()
atemp=0d0
cascadeTemp=>CascadeList

do while(associated(cascadeTemp))
	atemp=atemp+1d0/dble(numCascades)
	if(atemp >= r) then
		exit
	endif
	cascadeTemp=>cascadeTemp%nextCascade
end do
end subroutine

!***************************************************************************************************
!
!> Integer function cascadeCount()
!!
!! This subroutine counts how many cascades are active in the LOCAL mesh (not on all processors)
!! and returns that value
!!
!! Inputs: none
!! Outputs: number of cascades present
!
!***************************************************************************************************

integer function CascadeCount()
use mod_srscd_constants
use DerivedType
implicit none

type(cascade), pointer :: CascadeCurrent
integer count

CascadeCurrent=>ActiveCascades
	
count=0

do 10 while(associated(CascadeCurrent))
	count=count+1
	
	CascadeCurrent=>CascadeCurrent%next
10 continue

CascadeCount=count

end function

!***************************************************************************************************
!
!> subroutine addCascadeExplicit
!!
!! This subroutine takes the place of chooseReaction in the case of explicit cascade implantation.
!! It forces the program to 'choose' a cascade reaction, instead of using the Monte Carlo algorithm
!! to do so.
!!
!! Inputs: none
!! Outputs: reactionCurrent and CascadeCurrent, pointing at cascade reaction.
!
!***************************************************************************************************

subroutine addCascadeExplicit(reactionCurrent)
use DerivedType
use mod_srscd_constants
use randdp
implicit none

type(reaction), pointer :: reactionCurrent
double precision r2, atemp, r2timesa
integer i

!***************************************************************************************************
!Choose from cascades within the coarse mesh (choose which element to implant cascades in)
!***************************************************************************************************

atemp=0d0
r2=dprand()
r2timesa=r2*numCells

outer: do i=1,numCells

	reactionCurrent=>reactionList(i)	!Point reactionCurrent at the cascade implantation reaction in each volume element
	
	atemp=atemp+1d0		!Here we don't have a reaction rate, and each cell is weighted evenly (assuming uniform mesh)
	if(atemp >= r2timesa) then
		exit outer			!exit both loops with reactionCurrent pointing to the randomly chosen reaction
	endif

end do outer

!Checking that reactionCurrent is pointed at the correct reaction
if(implantType=='Cascade') then
	if(reactionCurrent%numReactants==0 .OR. reactionCurrent%numReactants==-10) then
		if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
			numImplantEvents=numImplantEvents+1
		else if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts /= 0) then !He implantation
			write(*,*) 'Error helium implantation in explicit procedure'
		else
			write(*,*) 'Error reaction not allowed ', 'reactants', reactionCurrent%numReactants, &
				'products', reactionCurrent%numProducts, 'rate', reactionCurrent%reactionRate
		end if
	end if
else
	write(*,*) 'Error wrong implant type for explicit procedure'
end if

end subroutine

!***************************************************************************************************
!> logical function cascadeMixingCheck()
!!
!! cascadeMixingCheck checks whether a given defect in the fine mesh interacts with the cascade. It 
!! returns a logical value of true or false.
!!
!! Input: none
!! Output: logical value for whether or not a fine mesh defect combines with a cascade.
!***************************************************************************************************

logical function cascadeMixingCheck()
use mod_srscd_constants
use randdp
implicit none

double precision r1, probability
logical boolean

!step 1: calculate the volume fraction of defectTemp in the cell

!Previous version (n*m cascade mixing checks)
!probability=cascadeVolume/(numCellsCascade*CascadeElementVol*numDisplacedAtoms)

!Current version (n cascade mixing checks)
probability=cascadeVolume/(numCellsCascade*CascadeElementVol)

!write(*,*) 'mixing probability', probability

!step 4: use random number to decide on interaction
r1=dprand()
if(r1 .GT. probability) then
	boolean=.FALSE.
else
	boolean=.TRUE.
endif

cascadeMixingCheck=boolean

end function

!***************************************************************************************************
!
!> Subroutine CascadeUpdateStep(cascadeCell)
!!
!! Carries out the communication necessary when a cascade is created or destroyed in an element
!! that bounds another processor.
!!
!! Step 1: Check to see if cascade was created/destroyed in element that is in the boundary of another
!! 	processor.
!!
!! Step 2: Send message to neighboring processors about whether a cascade was created/destroyed in their
!!	boundaries, and if so the number of defects to recieve in the boundary element
!!
!! Step 3: Send/recieve boundary element defects (just completely re-write boundary element defects in this
!!	step)
!!
!! Step 4: Update reaction rates for all diffusion reactions from elements neighboring the boundary
!!	elements which have been updated
!!
!! Inputs: cascadeCell (integer, 0 if no cascade, otherwise gives number of volume element that cascade
!!	event has occurred in)
!!
!! Outputs: none
!!
!! Actions: see above, sends/recieves information on boundary updates and updates reaction lists.
!
!***************************************************************************************************

subroutine cascadeUpdateStep(cascadeCell)
use DerivedType
use mod_srscd_constants
use ReactionRates
implicit none
include 'mpif.h'

integer cascadeCell

integer i, j, k, tag, cellInfoBuffer(6,5), count 
integer status(MPI_STATUS_SIZE), request
integer, allocatable :: defectBuffer(:,:)
type(defect), pointer :: defectCurrent, defectPrev

integer cellNumber, bndryCellNumber, numDefects, cellVol, defectType(numSpecies)
integer numSend, numRecv, bufferCount, sendBufferCount, sendCount, numBuffersRecv
integer recvInfoBuffer(5)
integer, allocatable :: recvDefectBuffer(:,:)
integer localGrainID, neighborGrainID

interface
	subroutine findDefectInList(defectCurrent, defectPrev, defectType)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer defectType(numSpecies)
	end subroutine
end interface

! Step 1: Check to see if cascade was created/destroyed in element that is in the boundary of another
! 	processor.

if(cascadeCell==0) then

	!do nothing, no cascade
	!write(86,*) 'No Cascade'
	do 22 j=1,6
		cellInfoBuffer(j,1)=0
		cellInfoBuffer(j,2)=0
		cellInfoBuffer(j,3)=0
		cellInfoBuffer(j,4)=0
		cellInfoBuffer(j,5)=0
		
		!Only send/recv if the neighboring proc is different from this one
		if(myProc%procNeighbor(j) .NE. -1 .AND. myProc%procNeighbor(j) .NE. myProc%taskid) then
			call MPI_SEND(cellInfoBuffer(j,:), 5, MPI_INTEGER,myProc%procNeighbor(j),251+j,MPI_COMM_WORLD,ierr)
		endif
		
	22 continue

else
	!write(86,*) 'Cascade Update Step 1: creating and sending buffers'
	do 11 j=1,6
		
		do 12 k=1,myMesh(cascadeCell)%numNeighbors(j)
			
			if(myMesh(cascadeCell)%neighborProcs(j,k) .NE. myProc%taskid .AND. &
				myMesh(cascadeCell)%neighborProcs(j,k) .NE. -1) then
				
				!this element has a neighbor that is on a different processor
				
				!cell number in neighboring proc
				cellInfoBuffer(j,1)=myMesh(cascadeCell)%neighbors(j,k)
				
				!local cell number
				cellInfoBuffer(j,4)=cascadeCell
				
				!Find out how many defects are in this volume element
				defectCurrent=>defectList(cascadeCell)%next
				count=0
				sendBufferCount=0
				
				do 10 while(associated(defectCurrent))
					if(mod(count,maxBufferSize)==0) then	!we have filled up one buffer
						sendBufferCount=sendBufferCount+1
						count=0
					endif
					
					count=count+1
					
					defectCurrent=>defectCurrent%next
				10 continue
				
				cellInfoBuffer(j,2)=count
				cellInfoBuffer(j,5)=sendBufferCount
				
				!Update volume (to send to boundary mesh in neighbor): after cascade addition/deletion,
				!cell volume will have changed
				cellInfoBuffer(j,3)=myMesh(cascadeCell)%volume
				
				!Send info to neighbor: first send size of defectBuffer, then send defectBuffer
				!if(myProc%taskid==MASTER) write(*,*) j, 'sending buffer', 101+j
				
				!write(86,*) 'Sending full buffer, proc', myProc%taskid, 'to', &
				!	myMesh(cascadeCell)%neighborProcs(j,k), 'tag', 251+j
				!write(86,*) 'dir', j, 'numNeighbors', myMesh(cascadeCell)%numNeighbors(j), 'k', k
				!write(86,*) 'cellInfoBuffer', cellInfoBuffer(j,:)
				
!				if(cellInfoBuffer(j,2) .GE. 202) then
!					write(*,*) 'CellInfoBuffer greater than 202'
!					write(*,*) myProc%taskid, 'to', myMesh(cascadeCell)%neighborProcs(j,k)
!				endif
					
				call MPI_SEND(cellInfoBuffer(j,:), 5, MPI_INTEGER,myMesh(cascadeCell)%neighborProcs(j,k),&
					251+j,MPI_COMM_WORLD,ierr)
				
				defectCurrent=>defectList(cascadeCell)%next
				
				do 25 bufferCount=1,sendBufferCount
				
					!All defect buffers have maxBufferSize defects except for the last
					if(bufferCount==sendBufferCount) then
						allocate(defectBuffer(cellInfoBuffer(j,2),numSpecies+1))
						numSend=cellInfoBuffer(j,2)
					else
						allocate(defectBuffer(maxBufferSize,numSpecies+1))
						numSend=maxBufferSize
!						write(*,*) 'Chopped down to ',maxBufferSize,'Proc', myProc%taskid
					endif
					
					!Re-loop through defects and add each type to the defect buffer
					sendCount=0
					
					do 13 while(associated(defectCurrent))
						sendCount=sendCount+1
						
						if(sendCount .GT. maxBufferSize) then
							exit
						endif
						
						do 14 i=1,numSpecies
							defectBuffer(sendCount,i)=defectCurrent%defectType(i)
						14 continue
							
						defectBuffer(sendCount, numSpecies+1)=defectCurrent%num
						defectCurrent=>defectCurrent%next	
					13 continue
					
					!write(86,*) 'Info sent, sending buffer. Size', cellInfoBuffer(j,2)*(numSpecies+1)
					
					call MPI_ISEND(defectBuffer, numSend*(numSpecies+1), MPI_INTEGER, &
						myMesh(cascadeCell)%neighborProcs(j,k), 351*j+bufferCount, MPI_COMM_WORLD,request, ierr)
						
					!write(86,*) 'Buffer sent'
	
					deallocate(defectBuffer)
					
				25 continue
			
			else if(myMesh(cascadeCell)%neighborProcs(j,k) .NE. -1 .AND. &
				myProc%procNeighbor(j) .NE. myProc%taskid) then
				
				!Tag shows that neighboring cell is not in another processor
				cellInfoBuffer(j,1)=0
				cellInfoBuffer(j,2)=0
				cellInfoBuffer(j,3)=0
				cellInfoBuffer(j,4)=0
				cellInfoBuffer(j,5)=0
				
				!write(86,*) 'Sending empty buffer, proc', myProc%taskid, 'to', &
				!	myProc%procNeighbor(j)
				
				call MPI_SEND(cellInfoBuffer(j,:), 5, MPI_INTEGER,myProc%procNeighbor(j),&
					251+j,MPI_COMM_WORLD,ierr)
				
			else if(myMesh(cascadeCell)%neighborProcs(j,k) == -1) then
				
				!Here, even if we are at a free surface, we send a blank buffer set to the other proc
				!because the proc mesh is periodic even when the actual mesh has free surfaces.
				
				cellInfoBuffer(j,1)=0
				cellInfoBuffer(j,2)=0
				cellInfoBuffer(j,3)=0
				cellInfoBuffer(j,4)=0
				cellInfoBuffer(j,5)=0
				
				!write(86,*) 'Sending empty buffer, proc', myProc%taskid, 'to', &
				!	myProc%procNeighbor(j)
				
				call MPI_SEND(cellInfoBuffer(j,:), 5, MPI_INTEGER,myProc%procNeighbor(j),&
					251+j,MPI_COMM_WORLD,ierr)
				
				!Do nothing, free surface
				
			else if(myProc%procNeighbor(j)==myProc%taskid) then
			
				!Do nothing PBCs point towards same proc
				
				!write(86,*) 'doing nothing, periodic BCs'
				
			else
				
				write(*,*) 'error sending in cascadeUpdateStep'
				
			endif
			
		12 continue
	11 continue
endif

!Recieve cellInfoBuffer and boundary element defects, if any
!write(86,*) 'Cascade Update Step 2: Recieving Info'
do 15 i=1,6
	
	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	!if the neighboring proc is not a free surface (no proc)
	if(myProc%procNeighbor(i) .NE. -1 .AND. myProc%procNeighbor(i) .NE. myProc%taskid) then
		
		!write(86,*) 'Recieving info proc', myProc%taskid, 'from', myProc%procNeighbor(i), 'tag', 251+tag
		
		!if(myProc%taskid==MASTER) write(*,*) i, 'receiving', 101+tag
		call MPI_RECV(recvInfoBuffer,5,MPI_INTEGER,myProc%procNeighbor(i),&
			251+tag,MPI_COMM_WORLD,status,ierr)
		
		!write(86,*) 'Info Recvd'
	
		!LOCAL cell number
		cellNumber=recvInfoBuffer(1)
		numDefects=recvInfoBuffer(2)
		cellVol=recvInfoBuffer(3)
		bndryCellNumber=recvInfoBuffer(4)
		numBuffersRecv=recvInfoBuffer(5)
		
!		if(numDefects .GE. 202) then
!			write(*,*) 'Proc', myProc%taskid, 'recieving defects GT 202 from', myProc%procNeighbor(i)
!		endif
		
		!if(myProc%taskid==MASTER) write(*,*) 'local cell number', cellNumber
		
		if(cellNumber .NE. 0) then
			
			!Dividing up buffers in order to not surpass the maximum amount
			!of data that can be sent using MPI_ISEND. All buffers except for
			!the last will have maxBufferSize*(numSpecies+1) integers in them

			!write(86,*) 'Recieving info buffer proc', myProc%taskid, &
			!	'from', myProc%procNeighbor(i)
			
			!Add defects in recvDefectBuffer to correct boundary element (remove all other defects from boundary element first)
			!and change the volume of the bondary mesh element
			myBoundary(i,bndryCellNumber)%volume=cellVol
			
			!remove defects from myBoundary (except for first defect, this is all 0's and is just a placeholder)
			defectCurrent=>myBoundary(i,bndryCellNumber)%defectList%next
			
			do 21 while(associated(defectCurrent))
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
				deallocate(defectPrev%defectType)
				deallocate(defectPrev)
			21 continue
			
			!nullify the %next pointer in the first element of the defect list
			defectCurrent=>myBoundary(i,bndryCellNumber)%defectList
			
			nullify(defectCurrent%next)
			
			do 24 bufferCount=1,numBuffersRecv
				
				!All buffers other than the final one are full (maxBufferSize defects)
				if(bufferCount==numBuffersRecv) then
					numRecv=numDefects
				else
					numRecv=maxBufferSize
				endif
				
				!if(myProc%taskid==MASTER) write(*,*) 'recieving buffer', bufferCount, 'numDefectsRecv', numRecv
					
				allocate(recvDefectBuffer(numRecv,numSpecies+1))
				
				call MPI_IRECV(recvDefectBuffer,numRecv*(numSpecies+1),MPI_INTEGER,myProc%procNeighbor(i),&
					351*tag+bufferCount,MPI_COMM_WORLD,request,ierr)
					
				call MPI_WAIT(request, status, ierr)
				
				!write(86,*) 'Buffer recieved'
			
				!add defects in recvDefectBuffer to boundary element
				do 16 j=1,numRecv
					
					defectCurrent=>myBoundary(i,bndryCellNumber)%defectList
					
					do 17 k=1,numSpecies
						defectType(k)=recvDefectBuffer(j,k)
					17 continue
					
					!if(myProc%taskid==MASTER) write(*,*) 'defect type', defectType, 'num', recvDefectBuffer(j,numSpecies+1)
					
					nullify(defectPrev)
					call findDefectInList(defectCurrent, defectPrev, defectType)
					
					!*********************************************************************************
					!NOTE: the next step should be unnecessary since there should only be one defect of
					!each type and they should all be new in myBoundary(dir,cellNumber), but we will include
					!it anyway for the sake of safety
					!*********************************************************************************
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						
						count=0
						do 18 k=1,numSpecies
							if(defectCurrent%defectType(k)==defectType(k)) then
								count=count+1
							endif
						18 continue
						if(count==numSpecies) then
						!	write(*,*) 'error defect already present in list, cascadeUpdateStep proc', myProc%taskid
						!	write(*,*) 'defectType', defectType
						!	if(myProc%taskid==MASTER) read(*,*)
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=cellNumber
							defectPrev%num=recvDefectBuffer(j,numSpecies+1)
							do 19 k=1,numSpecies
								defectPrev%defectType(k)=defectType(k)
							19 continue
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
					else 			!add a defect to the end of the list
						
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=cellNumber
						defectPrev%num=recvDefectBuffer(j,numSpecies+1)
						do 20 k=1,numSpecies
							defectPrev%defectType(k)=defectType(k)
						20 continue
					endif
					
				16 continue
					
				!memory clear
			
				deallocate(recvDefectBuffer)
				
			24 continue
				
			!Update diffusion rates from LOCAL element into recvInfoBuffer(1) in boundary
			
			!*******************
			!Add Diffusion reactions
			!*******************
			
			!point defectCurrent at defect list in local cell
			defectCurrent=>defectList(cellNumber)
			
			do 23 while(associated(defectCurrent))
				if (myMesh(cellNumber)%numNeighbors(i)==0) then
					write(*,*) 'error myMesh does not have neighbors in this direction'
				endif
				
				if(polycrystal=='yes') then
				
					!Find the grain ID number of the volume element we are in
					localGrainID=myMesh(cellNumber)%material
					
					!Find the grain ID number of the neighboring volume element
					!NOTE: here we don't need to worry about free surfaces, since
					!we are only adding diffusion reactions due to defects that
					!have changed on the boundary of this processor (in another
					!processor, not a free surface)
					if(myProc%procNeighbor(i) .NE. myProc%taskid .AND. &
						myProc%procNeighbor(i) .NE. -1) then
						neighborGrainID=myBoundary(i,myMesh(cellNumber)%neighbors(i,1))%material		
					else
						neighborGrainID=myMesh(myMesh(cellNumber)%neighbors(i,1))%material
					endif
					
					if(localGrainID==neighborGrainID) then
				
						!Allow diffusion between elements in the same grain
						call addDiffusionReactions(cellNumber, bndryCellNumber,&
							myProc%taskid, myProc%procNeighbor(i),i,defectCurrent%defectType)
					
					else
					
						!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
						call addDiffusionReactions(cellNumber, 0, myProc%taskid, -1, i, defectCurrent%defectType)
					
					endif
				
				else	
					!Add diffusion reactions from this cell to neighboring cells
					call addDiffusionReactions(cellNumber, bndryCellNumber,&
						myProc%taskid, myProc%procNeighbor(i),i,defectCurrent%defectType)
				endif
				
				defectCurrent=>defectCurrent%next
			23 continue
		
		endif
		
	else
	
		!Do nothing, free surface
	
	endif
	
15 continue

end subroutine

!***************************************************************************************************
!> subroutine createCascadeConnectivity()
!!
!! This subroutine assigns values to the connectivity matrix (global variable) used for all cascades
!!
!! Input: numxcascade, numycascade, nunmzcascade (global variables) : from parameters.txt
!! Output: cascadeConnectivity (global variable)
!***************************************************************************************************

subroutine createCascadeConnectivity()
use mod_srscd_constants

implicit none

integer cell
!************************************************
!PBCs in x and y, free in z (cell 0 represents free surface)
!************************************************
do cell=1,numCellsCascade
	if(mod(cell,numxcascade)==0) then !identify cell to the right
		!cascadeConnectivity(cell, 1)=cell-numxcascade+1
		cascadeConnectivity(cell, 1)=0	!free in x
	else
		cascadeConnectivity(cell,1)=cell+1
	end if
	
	if(mod(cell+numxcascade-1,numxcascade)==0) then !identify cell to the left
		!cascadeConnectivity(cell,2)=cell+numxcascade-1
		cascadeConnectivity(cell,2)=0	!free in x
	else
		cascadeConnectivity(cell,2)=cell-1
	end if
	
	if(mod(cell,numxcascade*numycascade) > numxcascade*(numycascade-1) .OR. &
		mod(cell,numxcascade*numycascade)==0) then
		cascadeConnectivity(cell,3)=0	!free in y
		!cascadeConnectivity(cell,3)=cell-(numxcascade*(numycascade-1))
	else
		cascadeConnectivity(cell,3)=cell+numxcascade
	end if
	
	if(mod(cell,numxcascade*numycascade) <= numxcascade .AND. mod(cell, numxcascade*numycascade) /= 0) then
		cascadeConnectivity(cell,4)=0	!free in y
		!cascadeConnectivity(cell,4)=cell+(numxcascade*(numycascade-1))
	else
		cascadeConnectivity(cell,4)=cell-numxcascade
	end if
	
	if(mod(cell,numxcascade*numycascade*numzcascade) > numxcascade*numycascade*(numzcascade-1) .OR. &
		mod(cell, numxcascade*numycascade*numzcascade)==0) then
		
		cascadeConnectivity(cell,5)=0	!free in z
	else
		cascadeConnectivity(cell,5)=cell+numxcascade*numycascade
	end if
	
	if(mod(cell,numxcascade*numycascade*numzcascade) <= numxcascade*numycascade .AND. &
		mod(cell,numxcascade*numycascade*numzcascade) /= 0) then
		
		cascadeConnectivity(cell,6)=0	!free in z
	else
		cascadeConnectivity(cell,6)=cell-numxcascade*numycascade
	end if
end do

end subroutine
