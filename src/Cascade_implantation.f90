
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
use mod_constants
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
	end if
	cascadeTemp=>cascadeTemp%nextCascade
end do
end subroutine

!***************************************************************************************************
!> Integer function cascadeCount()
!
! This subroutine counts how many cascades are active in the LOCAL mesh (not on all processors)
! and returns that value
!
! Inputs: none
! Outputs: number of cascades present
!***************************************************************************************************

integer function CascadeCount()
use mod_constants
use DerivedType
implicit none

type(cascade), pointer :: CascadeCurrent
integer count

CascadeCurrent=>ActiveCascades
	
count=0

do while(associated(CascadeCurrent))
	count=count+1
	
	CascadeCurrent=>CascadeCurrent%next
end do

CascadeCount=count

end function

!***************************************************************************************************
!> subroutine addCascadeExplicit
!
! This subroutine takes the place of chooseReaction in the case of explicit cascade implantation.
! It forces the program to 'choose' a cascade reaction, instead of using the Monte Carlo algorithm
! to do so.
!
! Inputs: none
! Outputs: reactionCurrent and CascadeCurrent, pointing at cascade reaction.
!***************************************************************************************************

subroutine addCascadeExplicit(reactionCurrent)
use DerivedType
use mod_constants
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

	if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
!		numImplantEvents=numImplantEvents+1
		numImpAnn(1)=numImpAnn(1)+1
	end if

else
	write(*,*) 'Error wrong implant type for explicit procedure'
end if

end subroutine

!***************************************************************************************************
!> logical function cascadeMixingCheck()
!
! cascadeMixingCheck checks whether a given defect in the fine mesh interacts with the cascade. It
! returns a logical value of true or false.
!
! Input: none
! Output: logical value for whether or not a fine mesh defect combines with a cascade.
!***************************************************************************************************

logical function cascadeMixingCheck()
use mod_constants
use randdp
implicit none

double precision r1, probability
logical boolean

!step 1: calculate the volume fraction of defectTemp in the cell

!Current version (n cascade mixing checks)
probability=cascadeVolume/(numCellsCascade*CascadeElementVol)
!write(*,*) 'mixing probability', probability

!step 4: use random number to decide on interaction
r1=dprand()
if(r1 > probability) then
	boolean=.FALSE.
else
	boolean=.TRUE.
endif

cascadeMixingCheck=boolean

end function

!***************************************************************************************************
!> Subroutine CascadeUpdateStep(cascadeCell)
!
! Carries out the communication necessary when a cascade is created or destroyed in an element
! that bounds another processor.
!
! Step 1: Check to see if cascade was created/destroyed in element that is in the boundary of another
! processor.
!
! Step 2: Send message to neighboring processors about whether a cascade was created/destroyed in their
! boundaries, and if so the number of defects to recieve in the boundary element
!
! Step 3: Send/recieve boundary element defects (just completely re-write boundary element defects in this step)
!
! Step 4: Update reaction rates for all diffusion reactions from elements neighboring the boundary
!	elements which have been updated
!
! Inputs: cascadeCell (integer, 0 if no cascade, otherwise gives number of volume element that cascade
!	event has occurred in)
!
! Outputs: none
!
! Actions: see above, sends/recieves information on boundary updates and updates reaction lists.
!***************************************************************************************************

subroutine cascadeUpdateStep(cascadeCell)
use DerivedType
use mod_constants
use ReactionRates
implicit none
include 'mpif.h'

integer cascadeCell

type(defect), pointer :: defectCurrent, defectPrev

integer i, j, k, tag, count, recvDir
integer cellNumber, bndryCellNumber, numSend, numRecv

integer localGrainID, neighborGrainID

double precision, allocatable :: defectSend(:,:), defectRecv(:,:)

integer status(MPI_STATUS_SIZE)
integer sendRequest, recvRequest
integer sendStatus(MPI_STATUS_SIZE),recvStatus(MPI_STATUS_SIZE)

interface
	subroutine findDefectInList(defectCurrent, defectPrev, defectType)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer defectType(numSpecies)
	end subroutine
end interface

! Step 1: Check to see if cascade was created/destroyed in element that is in the boundary of another processor.
if(cascadeCell==0) then

    numSend=0
    allocate(defectSend(numSpecies+1,numSend))

    do i=1,6

        if(mod(i,2)==0) then
            recvDir=i-1
        else
            recvDir=i+1
        end if

        if(myProc%procNeighbor(i) /= myProc%taskid) then
            call MPI_ISEND(defectSend, (numSpecies+1)*numSend, MPI_DOUBLE_PRECISION, myProc%procNeighbor(i), &
                    900+i, comm, sendRequest, ierr)
            !call MPI_WAIT(sendRequest, sendStatus, ierr)
            !call MPI_SEND(defectSend, (numSpecies+1)*numSend, MPI_DOUBLE_PRECISION, myProc%procNeighbor(i), &
            !        900+i, comm, ierr)
        end if

        if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

            count=0
            call MPI_PROBE(myProc%procNeighbor(recvDir), 900+i, comm,status,ierr)
            call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)

            numRecv=count/(numSpecies+1)
            allocate(defectRecv(numSpecies+1,numRecv))

            !call MPI_IRECV(defectRecv,(numSpecies+1)*numRecv,MPI_DOUBLE_PRECISION,myProc%procNeighbor(recvDir),&
            !        900+i,comm,recvRequest,ierr)
            !call MPI_WAIT(recvRequest, recvStatus, ierr)
            call MPI_RECV(defectRecv,(numSpecies+1)*numRecv,MPI_DOUBLE_PRECISION,myProc%procNeighbor(recvDir),&
                    900+i,comm,status,ierr)

            if(numRecv /= 0) then
                !Update my boundary
                bndryCellNumber=defectRecv(4,1)
                cellNumber=defectRecv(1,1)

                myBoundary(bndryCellNumber,recvDir)%volume=defectRecv(3,1)

                !remove defects from myBoundary (except for first defect, this is all 0's and is just a placeholder)
                defectCurrent=>myBoundary(bndryCellNumber,recvDir)%defectList%next

                !delete exiting defects
                nullify(defectPrev)
                do while(associated(defectCurrent))
                    defectPrev=>defectCurrent
                    defectCurrent=>defectCurrent%next
                    deallocate(defectPrev%defectType)
                    deallocate(defectPrev)
                end do

                !add defects to my boundary
                defectCurrent=>myBoundary(bndryCellNumber,recvDir)%defectList
                do j=2,numRecv
                    nullify(defectCurrent%next)
                    allocate(defectCurrent%next)
                    nullify(defectCurrent%next%next)
                    defectCurrent=>defectCurrent%next
                    allocate(defectCurrent%defectType(numSpecies))
                    defectCurrent%cellNumber=bndryCellNumber
                    defectCurrent%num=defectRecv(numSpecies+1,j)

                    defectCurrent%defectType(:)=defectRecv(1:numSpecies,j)
                end do

                !*******************
                !Add Diffusion reactions
                !*******************
                defectCurrent=>defectList(cellNumber)
                do while(associated(defectCurrent))
                    if (myMesh(cellNumber)%numNeighbors(recvDir)==0) then
                        write(*,*) 'error myMesh does not have neighbors in this direction'
                    end if

                    if(polycrystal=='yes') then

                        !Find the grain ID number of the volume element we are in
                        localGrainID=myMesh(cellNumber)%material

                        !Find the grain ID number of the neighboring volume element
                        !NOTE: here we don't need to worry about free surfaces, since we are only adding diffusion reactions due to defects that
                        !have changed on the boundary of this processor (in another processor, not a free surface)
                        if(myProc%procNeighbor(recvDir)/=myProc%taskid .AND. myProc%procNeighbor(recvDir)/=-1) then
                            neighborGrainID=myBoundary(myMesh(cellNumber)%neighbors(1,recvDir),recvDir)%material
                        else
                            neighborGrainID=myMesh(myMesh(cellNumber)%neighbors(1,recvDir))%material
                        end if

                        if(localGrainID==neighborGrainID) then
                            !Allow diffusion between elements in the same grain
                            call addDiffusionReactions(cellNumber, bndryCellNumber,&
                                    myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                        else
                            !Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
                            call addDiffusionReactions(cellNumber, 0, myProc%taskid, -1, recvDir, defectCurrent%defectType)
                        end if
                    else
                        !Add diffusion reactions from this cell to neighboring cells
                        call addDiffusionReactions(cellNumber, bndryCellNumber,&
                                myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                    end if
                    defectCurrent=>defectCurrent%next
                end do
            end if
            deallocate(defectRecv)
        end if

        if(myProc%procNeighbor(i) /= myProc%taskid) then
            call MPI_WAIT(sendRequest, sendStatus, ierr)
        end if
    end do

else    !cascadeCell /= 0

    numSend=0
	defectCurrent=>defectList(cascadeCell)%next

	outer: do i=1,6
		do k=1,myMesh(cascadeCell)%numNeighbors(i)
			if(myMesh(cascadeCell)%neighborProcs(k,i) /= myProc%taskid .AND. &
					myMesh(cascadeCell)%neighborProcs(k,i) /= -1) then

				do while(associated(defectCurrent))
                    numSend=numSend+1
					defectCurrent=>defectCurrent%next
				end do
				exit outer
			end if
		end do
	end do outer

	if(numSend /=0 ) then
        numSend=numSend+1
		allocate(defectSend(numSpecies+1,numSend))

        defectSend(1,1)=0		!myMesh(cascadeCell)%neighbors(k,dir)
        defectSend(2,1)=0		!numDefects
        defectSend(3,1)=0d0	!myMesh(cascadeCell)%volume
        defectSend(4,1)=0		!cascadeCell
        defectSend(5,1)=0		!useless

		defectCurrent=>defectList(cascadeCell)%next

		do j=2,numSend
            defectSend(1:numSpecies,j)=defectCurrent%defectType(:)
            defectSend(numSpecies+1,j)=defectCurrent%num

			defectCurrent=>defectCurrent%next
		end do
	end if

    !Communication
    do i=1,6

        if(mod(i,2)==0) then
            recvDir=i-1
        else
            recvDir=i+1
        end if

        do k=1,myMesh(cascadeCell)%numNeighbors(i)

            if(myMesh(cascadeCell)%neighborProcs(k,i) /= myProc%taskid .AND. &
                    myMesh(cascadeCell)%neighborProcs(k,i) /= -1) then

                defectSend(1,1)=myMesh(cascadeCell)%neighbors(k,i)
                defectSend(2,1)=count
                defectSend(3,1)=myMesh(cascadeCell)%volume
                defectSend(4,1)=cascadeCell
            else
                numSend = 0
            end if
        end do

        !Send
        if(myProc%procNeighbor(i) /= myProc%taskid) then

            call MPI_ISEND(defectSend, (numSpecies+1)*numSend, MPI_DOUBLE_PRECISION, myProc%procNeighbor(i), &
                    900+i, comm, sendRequest, ierr)
            !call MPI_WAIT(sendRequest, sendStatus, ierr)
            !call MPI_SEND(defectSend, (numSpecies+1)*numSend, MPI_DOUBLE_PRECISION, myProc%procNeighbor(i), &
            !        900+i, comm, ierr)
        end if

        !Recv
        if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

            count=0
            call MPI_PROBE(myProc%procNeighbor(recvDir), 900+i, comm,status,ierr)
            call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)

            numRecv=count/(numSpecies+1)
            allocate(defectRecv(numSpecies+1,numRecv))

            !call MPI_IRECV(defectRecv,(numSpecies+1)*numRecv,MPI_DOUBLE_PRECISION,myProc%procNeighbor(recvDir),&
            !        900+i,comm,recvRequest,ierr)
            !call MPI_WAIT(recvRequest, recvStatus, ierr)
            call MPI_RECV(defectRecv,(numSpecies+1)*numRecv,MPI_DOUBLE_PRECISION,myProc%procNeighbor(recvDir),&
                    900+i,comm,status,ierr)

            if(numRecv /= 0) then

                bndryCellNumber=defectRecv(4,1)
                cellNumber=defectRecv(1,1)

                myBoundary(bndryCellNumber,recvDir)%volume=defectRecv(3,1)

                !remove defects from myBoundary (except for first defect, this is all 0's and is just a placeholder)
                defectCurrent=>myBoundary(bndryCellNumber,recvDir)%defectList%next

                !delete exiting defects
                nullify(defectPrev)
                do while(associated(defectCurrent))
                    defectPrev=>defectCurrent
                    defectCurrent=>defectCurrent%next
                    deallocate(defectPrev%defectType)
                    deallocate(defectPrev)
                end do

                !add defects to my boundary
                defectCurrent=>myBoundary(bndryCellNumber,recvDir)%defectList

                do j=2,numRecv
                    nullify(defectCurrent%next)
                    allocate(defectCurrent%next)
                    nullify(defectCurrent%next%next)
                    defectCurrent=>defectCurrent%next
                    allocate(defectCurrent%defectType(numSpecies))
                    defectCurrent%cellNumber=bndryCellNumber
                    defectCurrent%num=defectRecv(numSpecies+1,j)

                    defectCurrent%defectType(:)=defectRecv(1:numSpecies,j)
                end do

                !*******************
                !Add Diffusion reactions
                !*******************
                defectCurrent=>defectList(cellNumber)
                do while(associated(defectCurrent))
                    if (myMesh(cellNumber)%numNeighbors(recvDir)==0) then
                        write(*,*) 'error myMesh does not have neighbors in this direction'
                    end if

                    if(polycrystal=='yes') then

                        !Find the grain ID number of the volume element we are in
                        localGrainID=myMesh(cellNumber)%material

                        !Find the grain ID number of the neighboring volume element
                        !NOTE: here we don't need to worry about free surfaces, since we are only adding diffusion reactions due to defects that
                        !have changed on the boundary of this processor (in another processor, not a free surface)
                        if(myProc%procNeighbor(recvDir)/=myProc%taskid .AND. myProc%procNeighbor(recvDir)/=-1) then
                            neighborGrainID=myBoundary(myMesh(cellNumber)%neighbors(1,recvDir),recvDir)%material
                        else
                            neighborGrainID=myMesh(myMesh(cellNumber)%neighbors(1,recvDir))%material
                        end if

                        if(localGrainID==neighborGrainID) then
                            !Allow diffusion between elements in the same grain
                            call addDiffusionReactions(cellNumber, bndryCellNumber,&
                                    myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                        else
                            !Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
                            call addDiffusionReactions(cellNumber, 0, myProc%taskid, -1, recvDir, defectCurrent%defectType)
                        end if
                    else
                        !Add diffusion reactions from this cell to neighboring cells
                        call addDiffusionReactions(cellNumber, bndryCellNumber,&
                                myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                    end if
                    defectCurrent=>defectCurrent%next
                end do

            end if
            deallocate(defectRecv)
        end if

        if(myProc%procNeighbor(i) /= myProc%taskid) then
            call MPI_WAIT(sendRequest, sendStatus, ierr)
        end if

    end do
end if

end subroutine

!***************************************************************************************************
!> subroutine createCascadeConnectivity()
!
! This subroutine assigns values to the connectivity matrix (global variable) used for all cascades
!
! Input: numxcascade, numycascade, nunmzcascade (global variables) : from parameters.txt
! Output: cascadeConnectivity (global variable)
!***************************************************************************************************

subroutine createCascadeConnectivity()
use mod_constants

implicit none

integer cell
!************************************************
!PBCs in x and y, free in z (cell 0 represents free surface)
!************************************************
do cell=1,numCellsCascade
	if(mod(cell,numxcascade)==0) then !identify cell to the right
		!cascadeConnectivity(cell, 1)=cell-numxcascade+1
		cascadeConnectivity(1, cell)=0	!free in x
	else
		cascadeConnectivity(1, cell)=cell+1
	end if
	
	if(mod(cell+numxcascade-1,numxcascade)==0) then !identify cell to the left
		!cascadeConnectivity(cell,2)=cell+numxcascade-1
		cascadeConnectivity(2, cell)=0	!free in x
	else
		cascadeConnectivity(2, cell)=cell-1
	end if
	
	if(mod(cell,numxcascade*numycascade) > numxcascade*(numycascade-1) .OR. &
		mod(cell,numxcascade*numycascade)==0) then
		cascadeConnectivity(3, cell)=0	!free in y
		!cascadeConnectivity(cell,3)=cell-(numxcascade*(numycascade-1))
	else
		cascadeConnectivity(3, cell)=cell+numxcascade
	end if
	
	if(mod(cell,numxcascade*numycascade) <= numxcascade .AND. mod(cell, numxcascade*numycascade) /= 0) then
		cascadeConnectivity(4, cell)=0	!free in y
		!cascadeConnectivity(cell,4)=cell+(numxcascade*(numycascade-1))
	else
		cascadeConnectivity(4, cell)=cell-numxcascade
	end if
	
	if(mod(cell,numxcascade*numycascade*numzcascade) > numxcascade*numycascade*(numzcascade-1) .OR. &
		mod(cell, numxcascade*numycascade*numzcascade)==0) then
		
		cascadeConnectivity(5, cell)=0	!free in z
	else
		cascadeConnectivity(5, cell)=cell+numxcascade*numycascade
	end if
	
	if(mod(cell,numxcascade*numycascade*numzcascade) <= numxcascade*numycascade .AND. &
		mod(cell,numxcascade*numycascade*numzcascade) /= 0) then
		
		cascadeConnectivity(6, cell)=0	!free in z
	else
		cascadeConnectivity(6, cell)=cell-numxcascade*numycascade
	end if
end do

end subroutine
