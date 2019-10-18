!*****************************************************************************************
!>double precision generate timestep - chooses a timestep using random number and Monte Carlo algorithm
!!(this is a global timestep)
!*****************************************************************************************

double precision function GenerateTimestep()
use DerivedType
use mod_constants
use randdp
implicit none

double precision r1

r1=dprand()
GenerateTimestep=dlog(1d0/r1)/maxRate

end function

!*****************************************************************************************
!>subroutine choose reaction - chooses a reaction in each processor according to the Monte
!!Carlo algorithm (this is a local reaction)
!!
!!NOTE: to increase computation speed, need these lists to be ordered.
!!
!!NOTE: to increase computation speed, choose an element first and then choose a reaction
!!within that element.
!*****************************************************************************************

subroutine chooseReaction(reactionCurrent, CascadeCurrent)
use DerivedType
use mod_constants
use randdp
implicit none

type(reaction), pointer :: reactionCurrent, reactionTemp
type(cascade), pointer :: cascadeCurrent
type(defect), pointer :: defectTemp
double precision r2, atemp, atemp_cell, r2timesa, atemp_test
integer i, j
double precision totalRateCells

!***************************************************************************************************
!Choose from reactions within the local mesh, including null event, according to kMC algorithm
!***************************************************************************************************

atemp=0d0
atemp_cell=0d0
atemp_test=0d0
r2=dprand()
r2timesa=r2*maxRate
nullify(CascadeCurrent)		!These are default pointed at nothing, indicating null event
nullify(reactionCurrent)	!These are default pointed at nothing, indicating null event

!***********************************************************************
!Choose from reactions within the coarse mesh
!***********************************************************************

outer: do i=1,numCells
	
	!used to choose a volume element that the reaction occurs within
	atemp_cell=atemp_cell+totalRateVol(i)
	
	!search for which volume element we are choosing among
	if(r2timesa <= atemp_cell) then

		!a reaction is chosen in this volume element, so we no longer nullify reactionCurrent
		reactionCurrent=>reactionList(i)
		
		!search for which reaction occurs in that volume element
		do while(associated(reactionCurrent))
			atemp=atemp+reactionCurrent%reactionRate
			if(r2timesa <= atemp) then
				exit outer	!exit both loops with reactionCurrent pointing to the randomly chosen reaction
			endif
			reactionCurrent=>reactionCurrent%next
		end do
		
	else
			atemp=atemp+totalRateVol(i)
	end if
end do outer

atemp_test = atemp_cell

!***********************************************************************
!Choose from reactions within the fine meshes for active cascades
!***********************************************************************
if(r2timesa <= atemp) then
	!do nothing, we have already chosen our reaction
else	!r2timesa > atemp
	CascadeCurrent=>ActiveCascades
	
	outer2: do while(associated(CascadeCurrent))
		
		!For now, just search through all cascade reaction lists for the reaction
		inner2: do i=1, numCellsCascade
		
			!Check if the reaction we are looking for is inside this cascade
			atemp_cell=atemp_cell+CascadeCurrent%totalRate(i)
			
			!If the reaction we are looking for is inside this cascade
			if(r2timesa <= atemp_cell) then
			
				reactionCurrent=>CascadeCurrent%reactionList(i)
				do while(associated(reactionCurrent))
					atemp=atemp+reactionCurrent%reactionRate
					if(r2timesa <= atemp) then
						exit outer2
					endif
					reactionCurrent=>reactionCurrent%next
				end do
						
			else
				
				atemp=atemp+CascadeCurrent%totalRate(i)	!skip this volume element
				
			endif
			
		end do inner2
		
		CascadeCurrent=>CascadeCurrent%next	!Loop over cascades

	end do outer2
end if

!***********************************************************************
!if we have not chosen a reaction at this point, we have a null event
!***********************************************************************

if(.NOT. associated(reactionCurrent)) then
	!write(*,*) 'null event chosen'
else
	!add to DPA and to numImplantEvents if we have an implantation event
	if(implantType=='FrenkelPair') then
		if(reactionCurrent%numProducts==2 .AND. reactionCurrent%numReactants==0) then	!Frenkel pair implantation
            numImpAnn(1)=numImpAnn(1)+1
		end if
	else if(implantType=='Cascade') then
		if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
			numImpAnn(1)=numImpAnn(1)+1
		end if
	else
		write(*,*) 'Error implantType'
	end if
	
	!for post processing: count the number of annihilated point defect (hard coded)
	if(reactionCurrent%numReactants==2) then	!clustering reaction
		if(reactionCurrent%reactants(2,1) /= 0 .AND. reactionCurrent%reactants(3,2) /= 0) then	!V+SIA_mobile
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(3,2))
		else if(reactionCurrent%reactants(2,1) /= 0 .AND. reactionCurrent%reactants(4,2) /= 0) then	!V+SIA_sessile
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(4,2))
		else if(reactionCurrent%reactants(3,1) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then	!SIA_mobile+V
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(3,1),reactionCurrent%reactants(2,2))
		else if(reactionCurrent%reactants(4,1) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then !SIA_sessile+V
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(4,1),reactionCurrent%reactants(2,2))
		end if
	endif
	
	!for post processing: counting trapping and emission reactions from the grain boundary
	if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then	!diffusion/emission reaction
  		if(atemp_cell == atemp_test .AND. reactionCurrent%cellNumber(2) > 0) then	!reactionCurrent%cellNumber(2)<0： diffuse from coarse mesh to fine mesh
			if(myMesh(reactionCurrent%cellNumber(1))%material == 1 .AND. &
				myMesh(reactionCurrent%cellNumber(2))%material == 2) then	!trapping on grain boundary
			
				numTrapV=numTrapV+reactionCurrent%reactants(2,1)
				numTrapSIA=numTrapSIA+reactionCurrent%reactants(3,1)
			
			else if(myMesh(reactionCurrent%cellNumber(1))%material==2 .AND. &
					myMesh(reactionCurrent%cellNumber(2))%material==1) then	!emission from grain boundary
				
				numEmitV=numEmitV+reactionCurrent%reactants(2,1)
				numEmitSIA=numEmitSIA+reactionCurrent%reactants(3,1)
			
			end if
  		end if
	end if
end if

end subroutine


!*****************************************************************************************
!>subroutine choose reaction - chooses a reaction in each processor according to the Monte
!!Carlo algorithm (this is a local reaction)
!!
!!NOTE: to increase computation speed, need these lists to be ordered.
!!
!!NOTE: to increase computation speed, choose an element first and then choose a reaction
!!within that element.
!*****************************************************************************************

subroutine chooseReactionSingleCell(reactionCurrent, CascadeCurrent, cell)
use DerivedType
use mod_constants
use randdp
implicit none

type(reaction), pointer :: reactionCurrent, reactionTemp
type(cascade), pointer :: cascadeCurrent
type(defect), pointer :: defectTemp
double precision r2, atemp, atemp_cell, r2timesa
integer i, j, cell
double precision totalRateCells

!***************************************************************************************************
!Choose from reactions within the local mesh, including null event, according to kMC algorithm
!***************************************************************************************************

atemp=0d0
atemp_cell=0d0
r2=dprand()
r2timesa=r2*maxRate
nullify(CascadeCurrent)		!CsacadeCurrent=NULL
nullify(reactionCurrent)	!reactionCurrent=NULL

!***********************************************************************
!Here we choose from reactions within the coarse mesh
!***********************************************************************
	
!used to choose a volume element that the reaction occurs within
atemp_cell=atemp_cell+totalRateVol(cell)

!search for which volume element we are choosing among
if(r2timesa <= atemp_cell) then

	reactionCurrent=>reactionList(cell)	!a reaction is chosen in this volume element, so we no longer nullify reactionCurrent
	
	!search for which reaction occurs in that volume element
	outer: do while(associated(reactionCurrent))
	
		atemp=atemp+reactionCurrent%reactionRate
		if(r2timesa <= atemp) then
			exit outer		!exit both loops with reactionCurrent pointing to the randomly chosen reaction
		endif
		reactionCurrent=>reactionCurrent%next
		
	end do outer
	
else
	!No reaction chosen in the cell
	atemp=atemp+totalRateVol(cell)	!This is not necessary but I will keep in case adaptive meshing is added to single element KMC later
endif

!NOTE: not allowing adaptive meshing for single element kMC algorithm (no need to search fine meshes for reactions)

!***********************************************************************
!if we have not chosen a reaction at this point, we have a null event
!***********************************************************************

if(.NOT. associated(reactionCurrent)) then
	write(*,*) 'null event chosen'
else
	!add to DPA and to numImplantEvents if we have an implantation event
	if(implantType=='FrenkelPair') then
		if(reactionCurrent%numProducts==2 .AND. reactionCurrent%numReactants==0) then
			!Frenkel pair implantation
!			numImplantEvents=numImplantEvents+1		!LOCAL number of implantation events. Total DPA calculated out using MPI_ALLREDUCE
            numImpAnn(1)=numImpAnn(1)+1
		end if
	else if(implantType=='Cascade') then
		if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
!			numImplantEvents=numImplantEvents+1
            numImpAnn(1)=numImpAnn(1)+1
		end if
	else
		write(*,*) 'Error implantType'
	end if
	
	!for post processing: count annihilation reactions (hard coded)
	if(reactionCurrent%numReactants==2) then	!clustering reaction
		if(reactionCurrent%reactants(2,1) /= 0 .AND. reactionCurrent%reactants(3,2) /= 0) then	!V+SIA_mobile
			
!			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(3,2))
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(3,2))
		
		else if(reactionCurrent%reactants(2,1) /= 0 .AND. reactionCurrent%reactants(4,2) /= 0) then	!V+SIA_sessile
		
!			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(4,2))
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(2,1),reactionCurrent%reactants(4,2))
		
		else if(reactionCurrent%reactants(3,1) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then	!SIA_mobile+V
		
!			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(3,1),reactionCurrent%reactants(2,2))
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(3,1),reactionCurrent%reactants(2,2))
		
		else if(reactionCurrent%reactants(4,1) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then !SIA_sessile+V
		
!			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(4,1),reactionCurrent%reactants(2,2))
            numImpAnn(2)=numImpAnn(2)+min(reactionCurrent%reactants(4,1),reactionCurrent%reactants(2,2))
		
		endif
	endif
endif

end subroutine

!*****************************************************************************************
!> Subroutine update defect list - updates defects according to reaction chosen
!
!This is the most involved and complex subroutine of the entire SRSCD algorithm.
!This subroutine updates defects in the local mesh according to the reaction chosen, and
!communicates with neighboring processors about defects that may have passed into
!a different processor as well as about defects that may have changed on the boundary.
!It also creates a list of defects that have been updated, to inform the next subroutine (updateReactionList(defectUpdate))
!which reactions to update.
!
!Input: reactionCurrent
!Output: defectUpdateCurrent, CascadeCurrent
!*****************************************************************************************

!***************************************************************************************************
!06/04/2014: Discovered bug in updateDefectList: if two reactions send the same defect to the same cell,
!the defect numbers will not both add correctly because defectCurrentUpdate%num is used for the
!defect number, which means that it does not recognize whether another defect has just been added
!in the same cell.
!
!To fix, I am changing the way that the numbers of defects added/subtracted is calculated. Here,
!defectCurrentUpdate%num and buffers(:,:,numSpecies+2) will contain +/- 1, indicating adding
!or removing one of these defects. If a reaction adds or removes more than one of this type of defect,
!multiple defectCurrentUpdate entries will be created.
!
!Thus, the way that defect numbers are tracked needs to be changed in this subroutine, and the way
!that updateReactionList updates reaction rates needs to be checked to make sure that it is using the
!actual defect numbers and not defectCurrentUpdate%num.
!***************************************************************************************************

subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
use mod_constants
use DerivedType
use ReactionRates
use randdp
implicit none
include 'mpif.h'

type(reaction), pointer :: reactionCurrent
type(defect), pointer :: defectCurrent, defectPrev, defectTemp, defectStoreList, defectStore, defectStorePrev, defectTest
type(cascadeDefect), pointer :: cascadeDefectTemp
type(defectUpdateTracker), pointer :: defectUpdateCurrent
type(cascade), pointer :: CascadeCurrent, CascadePrev
type(cascadeEvent), pointer :: cascadeTemp

!Used for cascade
double precision coordinatesTemp(3)
integer cellNumber, mixingEvents, mixingTemp
logical combineBoolean, isCombined
!Used for cascade recombination
double precision r1, atemp

integer i, j, k, l, m, same, products(numSpecies), product2(numSpecies), totalLocalRecv
double precision diffusionRandom
logical flag, flagTemp
integer count

!Function
integer findCellWithCoordinatesFineMesh, chooseRandomCell
logical cascadeMixingCheck

!Used for communication between processors
integer numUpdateSend(6), numUpdateRecv(6)	!the number of defects being sent to (recived from) each processor neighbor
integer numLocalRecv, numBndryRecv			!the number of (loceal/bndry) defects being recieved from each proc neighbor

!NOTE: this final step could be eliminated by keeping the global mesh in each local processor
!(thus each element could be identified as being part of the local mesh of one proc and the boundary of any other procs)
integer numUpdateFinal(6), numUpdateFinalRecv(6)	!number of defects being sent/recieved in final update step
integer recvDir

!create buffers of information to send to neighboring elements (Contains local defects and boundary defects)
integer, allocatable :: firstSend(:,:,:)
integer, allocatable :: firstRecv(:,:,:)
integer, allocatable :: finalSend(:,:,:)	!Only contains local defects
integer, allocatable :: finalBufferRecv(:,:)

integer status(MPI_STATUS_SIZE)
integer sendFirstStatus(MPI_STATUS_SIZE), recvFirstStatus(MPI_STATUS_SIZE)
integer sendFinalStatus(MPI_STATUS_SIZE), recvFinalStatus(MPI_STATUS_SIZE)
integer sendFirstRequest(6), recvFirstRequest(6)
integer sendFinalRequest(6), recvFinalRequest(6)

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
	
	subroutine chooseCascade(CascadeTemp)
		use DerivedType
		implicit none
		type(cascadeEvent), pointer :: CascadeTemp
	end subroutine
	
	subroutine initializeFineMesh(CascadeCurrent)
		use DerivedType
		type(cascade), pointer :: CascadeCurrent
	end subroutine
end interface

!initialize arraies
numUpdateSend=0
numUpdateRecv=0
numUpdateFinal=0
numUpdateFinalRecv=0

!initialize variables
totalLocalRecv=0
mixingEvents=0

!create buffers: size greater than max size needed (at most numReactants and numProducts to change)
!"6" is the number of directions
!"4" is the maximum number of defectts to be updated
!numSpecies+1 = local mesh
!numSpecies+2 = +1/-1
!numSpecies+3 = neighbor mesh (if the defect is in a local mesh), or
!			  =	-100 (if the defect is in a boundry mesh)
allocate(firstRecv(numSpecies+3,4,6))

!if we have chosen a null event
if(.NOT. associated(reactionCurrent)) then

	allocate(firstSend(numSpecies+3,0,6))

	do i=1,6
		if(mod(i,2)==0) then
			recvDir=i-1
		else
			recvDir=i+1
		end if

		!Send
		if(myProc%procNeighbor(i) /= myProc%taskid) then
			call MPI_ISEND(firstSend,0,MPI_INTEGER,myProc%procNeighbor(i),200+i ,comm,sendFirstRequest(i),ierr)
		end if

		!Recv
		if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

			count=0
			call MPI_PROBE(myProc%procNeighbor(recvDir), 200+i,comm,status,ierr)
			call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

			numUpdateRecv(recvDir) = count/(numSpecies+3)

			call MPI_IRECV(firstRecv(1,1,recvDir),(numSpecies+3)*numUpdateRecv(recvDir),MPI_INTEGER,&
					myProc%procNeighbor(recvDir),200+i ,comm,recvFirstRequest(recvDir),ierr)
		end if
	end do
!associated(reactionCurrent)	!if we have chosen an event
else

	!***********************************************************************************************
	!Cascade chosen
	!
	!Initialization of fine mesh: randomly select defects from the coarse mesh into the fine mesh.
	!***********************************************************************************************
	if(reactionCurrent%numReactants==-10) then

        !Communication common defect
        allocate(firstSend(numSpecies+3,0,6))

        do i=1,6
            if(mod(i,2)==0) then
                recvDir=i-1
            else
                recvDir=i+1
            end if

            !Send
            if(myProc%procNeighbor(i) /= myProc%taskid) then
                call MPI_ISEND(firstSend,0,MPI_INTEGER,myProc%procNeighbor(i),200+i ,comm,sendFirstRequest(i),ierr)
            end if

            !Recv
            if(myProc%procNeighbor(recvDir) /= myProc%taskid) then
                count=0
                call MPI_PROBE(myProc%procNeighbor(recvDir), 200+i,comm,status,ierr)
                call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

                numUpdateRecv(recvDir) = count/(numSpecies+3)

                call MPI_IRECV(firstRecv(1,1,recvDir),(numSpecies+3)*numUpdateRecv(recvDir),MPI_INTEGER,&
                        myProc%procNeighbor(recvDir),200+i ,comm,recvFirstRequest(recvDir),ierr)
            end if
        end do

        !Update cascades
        if(meshingType=='adaptive') then
			
			!**************************************************
            !> choose a cascade
			!> create fine mesh
			!> populate with defects from coarse mesh element
			!> update defectUpdate
			!> update localSend if needed
			!> populate with defects from cascade
			!**************************************************

			call chooseCascade(cascadeTemp)	!choose one of the cascades from the list randomly
			cascadeDefectTemp=>cascadeTemp%ListOfDefects

			!initialize the cascade to be added
			CascadeCurrent=>ActiveCascades
			if(.NOT. associated(CascadeCurrent)) then		!no cascade fine meshes are currently in the system
				allocate(ActiveCascades)
				nullify(ActiveCascades%next)
				nullify(ActiveCascades%prev)
				nullify(ActiveCascades%localDefects)
				nullify(ActiveCascades%reactionList)
				ActiveCascades%cellNumber=reactionCurrent%cellNumber(1) !cell number in COARSE MESH of this cascade
				allocate(ActiveCascades%totalRate(numCellsCascade))
				do j=1,numCellsCascade
					ActiveCascades%totalRate(j)=0d0							!reaction rate of all reactions within fine mesh
				end do
				CascadeCurrent=>ActiveCascades
                CascadeCurrent%cascadeID=numImpAnn(1)
				!ActiveCascades%cascadeID=numImpAnn(1)
			else	!cascade fine meshe are already in the system
				j=1
				do while(associated(CascadeCurrent))
					CascadePrev=>CascadeCurrent
					CascadeCurrent=>CascadeCurrent%next
					j=j+1
				end do
				allocate(CascadeCurrent)
				nullify(CascadeCurrent%next)
				nullify(CascadeCurrent%localDefects)
				nullify(CascadeCurrent%reactionList)
				CascadeCurrent%prev=>CascadePrev
				CascadePrev%next=>CascadeCurrent
				CascadeCurrent%cellNumber=reactionCurrent%cellNumber(1)	!cell number in COARSE MESH of this cascade
				allocate(CascadeCurrent%totalRate(numCellsCascade))
				do j=1,numCellsCascade
					CascadeCurrent%totalRate(j)=0d0
				end do
                CascadeCurrent%cascadeID=numImpAnn(1)
			end if
			
			!*******************************************************************
			!when we create a cascade in this coarse mesh element, we reduce its
			!volume by the cascade volume
			!Each time a cascade is added, a new cascade-fine-mesh is created.
			!*******************************************************************
			
			!update volume of couarse mesh element
			myMesh(cascadeCurrent%cellNumber)%volume = myMesh(cascadeCurrent%cellNumber)%volume &
					-CascadeElementVol*dble(numCellsCascade)
			
			!If we add too many cascades to a coarse mesh element, the volume inside can become negative. If this 
			!happens, output an error message. (Only a danger in the case of very high DPA rates and
			!small coarse mesh elements)
			if(myMesh(cascadeCurrent%cellNumber)%volume <= 0d0) then
				write(*,*) 'Error negative coarse mesh volume'
			end if
			
			!*******************************************************************
			!initialize defect list and reaction list in CascadeCurrent. This includes defects
			!placed in the fine mesh from the coarse mesh.
			!*******************************************************************
			call initializeFineMesh(CascadeCurrent)

			!***************************************************************************************
			!Recombination step:
			!Here, we are taking the defects in the initialized fine mesh and combining them with
			!the defects in the cascade according to a probability given by the cascade volume
			!divided by the total fine mesh volume.
			!
			!Create a list of cascade defects that will be added to the simulation after the 
			!recombination step has been carried out
			!
			!This list starts on the second pointer (pointer 1 used for reference only)
			!***************************************************************************************
			
			allocate(defectStoreList)
			allocate(defectStoreList%defectType(numSpecies))
			do i=1, numSpecies
				defectStoreList%defectType(i)=0
			end do
			defectStoreList%num=0
			defectStoreList%cellNumber=0	!unknown

			nullify(defectStoreList%next)
			defectStore=>defectStoreList
			
			!***************************************************************************************
			!Step 1: Create defectStoreList with defects in cascade only
			!Step 2: For each defect in the fine mesh, check whether it combines with cascade
			!Step 3: If yes, remove it from the fine mesh and combine it randomly with one defect in
			!	defectStoreList
			!Step 4: Implant defects in defectStoreList in fine mesh. Remember to correctly update
			!	defect update protocol.
			!***************************************************************************************

			!Step 1:
			do i=1, cascadeTemp%numDefectsTotal
				!**************************************************************
				!Make sure coordinates of all defects in cascade are within fine mesh, using periodic BC's
				!**************************************************************
				coordinatesTemp(1)=cascadeDefectTemp%coordinates(1)
				if(coordinatesTemp(1) > numxcascade*finelength/2d0) then
					coordinatesTemp(1)=coordinatesTemp(1)-numxcascade*finelength
				else if(coordinatesTemp(1) < -numxcascade*finelength/2d0) then
					coordinatesTemp(1)=coordinatesTemp(1)+numxcascade*finelength
				endif
				
				coordinatesTemp(2)=cascadeDefectTemp%coordinates(2)
				if(coordinatesTemp(2) > numycascade*finelength/2d0) then
					coordinatesTemp(2)=coordinatesTemp(2)-numycascade*finelength
				else if(coordinatesTemp(2) < -numycascade*finelength/2d0) then
					coordinatesTemp(2)=coordinatesTemp(2)+numycascade*finelength
				endif
				
				coordinatesTemp(3)=cascadeDefectTemp%coordinates(3)
				if(coordinatesTemp(3) > numzcascade*finelength/2d0) then
					coordinatesTemp(3)=coordinatesTemp(3)-numzcascade*finelength
				else if(coordinatesTemp(3) < -numzcascade*finelength/2d0) then
					coordinatesTemp(3)=coordinatesTemp(3)+numzcascade*finelength
				endif	
				
				!find position cascade is implanted in fine mesh
				cellNumber=findCellWithCoordinatesFineMesh(coordinatesTemp)						
				
				!Create list of defects that need to be added to fine mesh (after casacade-fine mesh
				!mixing has been taken into account)
				allocate(defectStore%next)
				defectStore=>defectStore%next
				allocate(defectStore%defectType(numSpecies))
				nullify(defectStore%next)
				
				do j=1,numSpecies
					defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
				end do
				defectStore%cellNumber=cellNumber
				defectStore%num=1

				if(pointDefectToggle=='yes') then
					if(defectStore%defectType(3) > max3DInt) then
						defectStore%defectType(4) = defectStore%defectType(3)
						defectStore%defectType(3) = 0
					end if
				end if
				
				cascadeDefectTemp=>cascadeDefectTemp%next
			end do

			!combination
			do j=1,numCellsCascade
				defectCurrent=>CascadeCurrent%localDefects(j)
				
				defectPrev=>defectCurrent
				defectTemp=>defectCurrent%next
				do while(associated(defectTemp))
					
					!Step 2:
					!boolean variable, if FALSE then no combination with defectCurrent, if TRUE then combination
					mixingTemp=0
					if(defectTemp%num > 0) then
						do k=1,defectTemp%num
							combineBoolean=cascadeMixingCheck()
						
							if(combineBoolean .eqv. .TRUE.) then
								mixingTemp=mixingTemp+1
							end if
						end do
					end if
					
					!Step 3:
					do k=1,mixingTemp

						!Choose which defect will be combined wtih defectTemp
						nullify(defectStorePrev)
						defectStore=>defectStoreList%next
						
						atemp=0d0
						r1=dprand()
						
						!Move defectStore through defectStoreList until defect chosen for recombination
						do i=1,cascadeTemp%numDefectsTotal
							atemp=atemp+1d0/dble(cascadeTemp%numDefectsTotal)
							if(atemp > r1) then
								exit
							end if
							defectStorePrev=>defectStore
							defectStore=>defectStore%next
						end do
						
						if( .NOT. associated(defectStore)) then
							write(*,*) 'Error DefectStore not associated in cascade mixing'
						else
						    !***********************************************************************
						    !Hard coded: use defect combination rules to combine defectStore and defecTemp
						    !These rules have been transported to a separate subroutine
						    !in order to facilitate hard-coding and keep this subroutine clean.
						    !***********************************************************************
						    do i=1, numSpecies
							    product2(i)=0
						    end do

						    call defectCombinationRules(defectStore%defectType,product2, defectTemp, isCombined)

							if(product2(3)/=0 .OR. product2(4)/=0) then
							    !defectStore is at the middle (or end) of defectStoreList
							    if(associated(defectStorePrev)) then
								    nullify(defectStorePrev%next)
								    allocate(defectStorePrev%next)
								    nullify(defectStorePrev%next%next)
								    defectStorePrev=>defectStorePrev%next
								    allocate(defectStorePrev%defectType(numSpecies))
								    defectStorePrev%cellNumber=defectStore%cellNumber
								    defectStorePrev%num=1
								    do i=1, numSpecies
									    defectStorePrev%defectType(i)=product2(i)
								    end do
								    defectStorePrev%next=>defectStore
							    else	!defectStore is at the beginning of defectStoreList
								    defectStorePrev=>defectStoreList
								    nullify(defectStorePrev%next)
								    allocate(defectStorePrev%next)
								    nullify(defectStorePrev%next%next)
								    defectStorePrev=>defectStorePrev%next
								    allocate(defectStorePrev%defectType(numSpecies))
								    defectStorePrev%cellNumber=defectStore%cellNumber
								    defectStorePrev%num=1
								    do i=1, numSpecies
									    defectStorePrev%defectType(i)=product2(i)
								    end do
								    defectStorePrev%next=>defectStore
							    end if
						    end if

                            if(isCombined .eqv. .TRUE.	) then
                                mixingEvents=mixingEvents+1
                            end if

						end if
	
					end do

                    !NEXT: remove DefectTemp from the fine mesh.
                    !The number of defectTemp removed from the fine mesh is mixingEvents
                    if(mixingEvents < defectTemp%num) then
                        defectTemp%num=defectTemp%num-mixingEvents
                    else    !mixingEvents == efectTemp%num
                        if(defectTemp%num==0) then	!error
                            write(*,*) 'error defect num zero combining with cascade defect'
                        !if defectTemp is in the middle of the list
                        else if(associated(defectPrev) .AND. associated(defectTemp%next)) then
                            defectPrev%next=>defectTemp%next
                            deallocate(defectTemp%defectType)
                            deallocate(defectTemp)
                            defectTemp=>defectPrev
                         !if defectTemp is at the end of the list
                        else if(associated(defectPrev)) then
                            deallocate(defectTemp%defectType)
                            deallocate(defectTemp)
                            defectTemp=>defectPrev
                            nullify(defectPrev%next)
                        !if defectTemp is at the beginning of the list
                        else if(associated(defectTemp%next)) then
                            defectTemp%num=0
                         end if
                    end if

					defectPrev=>defectTemp
					defectTemp=>defectTemp%next
					
				end do
			end do

			defectStore=>defectStoreList%next
			!add all cascade mixing products and normal cascade products to the system
			do while(associated(defectStore))
			
				!**************************************************************
				!Temporary debug: checking that all defects added are admissable defects
				!**************************************************************
				count=0
				do j=1,numSpecies
					if(defectStore%defectType(j) /= 0) then
						count=count+1
					endif
				end do
				if(count > 2) then
					write(*,*) 'Error adding unadmissible defect to cascade'
					write(*,*) defectStore%defectType
					call MPI_ABORT(comm,ierr)
				end if
			
				defectCurrent=>CascadeCurrent%localDefects(defectStore%cellNumber)
				!STEP 2: add DefectTemp+products to the system
				count=0
				do j=1,numSpecies
					if(defectStore%defectType(j)==0) then
						count=count+1
					endif
				end do
				if(count==numSpecies) then
					!do nothing, we have total annihilation
				else
					
					nullify(defectPrev)
					call findDefectInList(defectCurrent, defectPrev, defectStore%defectType)
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						
						count=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==defectStore%defectType(j)) then
								count=count+1
							endif
						end do
						if(count==numSpecies) then
							defectCurrent%num=defectCurrent%num+1

						else		!if the defect is to be inserted in the list

							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectStore%cellNumber
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=defectStore%defectType(j)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectStore%cellNumber
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=defectStore%defectType(j)
						end do
					end if
					
				endif
				defectStore=>defectStore%next
			
			end do
			
			!Final step: set up defectUpdate for all defects in fine mesh
			do i=1,numCellsCascade
				defectTemp=>CascadeCurrent%localDefects(i)%next
				
				do while(associated(defectTemp))
				
					!*******************************************************************************
					!create a new element in defectUpdate and assign all variables
					!
					!Note: this will be done for each member of the cascade, even if more than 
					!one defect of this type is added to this volume element in the fine mesh.
					!
					!This means that the reaction list update will happen multiple times (redundant)
					!for some defects in some volume elements in the fine mesh.
					!*******************************************************************************
					
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(1)
					defectUpdateCurrent%dir=0	!not pointed at a different proc
					allocate(defectUpdateCurrent%defectType(numSpecies))
					do j=1,numSpecies
						defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
					end do
					defectUpdateCurrent%cellNumber=i
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID	!identify which cascade this defect is inside
					defectUpdateCurrent%num=1
					nullify(defectUpdateCurrent%next)
					
					defectTemp=>defectTemp%next
				
				end do
				
			end do

			!*******************************************************************
			!memory erase: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
			!*******************************************************************
			defectStore=>defectStoreList
			
			do while(associated(defectStore))
				defectTemp=>defectStore
				defectStore=>defectStore%next
				deallocate(defectTemp%defectType)
				deallocate(defectTemp)
			end do
			
			nullify(defectStore)
			nullify(defectTemp)
			nullify(defectStoreList)
			nullify(defectCurrent)
			nullify(defectPrev)	

		else if(meshingType=='nonAdaptive') then
			call chooseCascade(CascadeTemp)	!choose one of the cascades from the list randomly
			cascadeDefectTemp=>cascadeTemp%ListOfDefects
			!***************************************************************************************
			!Recombination step:
			!
			!Here, we are taking the defects in the initialized fine mesh and combining them with
			!the defects in the cascade according to a probability given by the cascade volume
			!divided by the total fine mesh volume.
			!
			!Create a list of cascade defects that will be added to the simulation after the 
			!recombination step has been carried out
			!
			!This list starts on the second pointer (pointer 1 used for reference only)
			!***************************************************************************************
			
			allocate(defectStoreList)
			allocate(defectStoreList%defectType(numSpecies))
			do i=1,numSpecies
				defectStoreList%defectType(i)=0
			end do
			defectStoreList%num=0
			defectStoreList%cellNumber=reactionCurrent%cellNumber(1)

			nullify(defectStoreList%next)
			defectStore=>defectStoreList
			
			!***************************************************************************************
			!Step 1: Create defectStoreList with defects in cascade only
			!Step 2: For each defect in the fine mesh, check whether it combines with cascade
			!Step 3: If yes, remove it from the fine mesh and combine it randomly with one defect in
			!	defectStoreList
			!Step 4: Implant defects in defectStoreList in fine mesh. Remember to correctly update
			!	defect update protocol.
			!***************************************************************************************

			!Step 1:
			do i=1, cascadeTemp%numDefectsTotal
				!Create list of defects that need to be added to cell (after casacade
				!mixing has been taken into account)
				
				allocate(defectStore%next)
				defectStore=>defectStore%next
				allocate(defectStore%defectType(numSpecies))
				nullify(defectStore%next)
				
				do j=1,numSpecies
					defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
				end do
				defectStore%cellNumber=reactionCurrent%cellNumber(1)
				defectStore%num=1

				if(pointDefectToggle=='yes') then
					if(defectStore%defectType(3) > max3DInt) then
						defectStore%defectType(4) = defectStore%defectType(3)
						defectStore%defectType(3) = 0
					end if
				end if
				
				cascadeDefectTemp=>cascadeDefectTemp%next
			end do
			
			defectCurrent=>defectList(reactionCurrent%cellNumber(1))	!exited defects in the CoarseMesh
			
			defectPrev=>defectCurrent
			defectTemp=>defectCurrent%next
			
			if(cascadeVolume > 0d0) then
			
				do while(associated(defectTemp))
					
					!Step 3: choose a defect in defectStoreList to combine with the defect that exited in the coarseMesh
					do k=1,defectTemp%num
					!	mixingEvents=mixingEvents+1
						
						!Choose which defect in the cascade will be combined with defectTemp already presetn in the cell
						nullify(defectStorePrev)
						defectStore=>defectStoreList%next
						
						atemp=0d0
						r1=dprand()
						
						!Move defectStore through defectStoreList until defect chosen for recombination
						do i=1,cascadeTemp%numDefectsTotal
							atemp=atemp+1d0/dble(cascadeTemp%numDefectsTotal)
							
							if(atemp > r1) then
								exit
							endif
							defectStorePrev=>defectStore
							defectStore=>defectStore%next
						end do
						
						if( .NOT. associated(defectStore)) then
							write(*,*) 'Error DefectStore not associated in cascade mixing'
						else	!associated(defectStore)

						    !***********************************************************************
						    !Hard coded: use defect combination rules to combine defectStore and defecTemp
						    !These rules have been transported to a separate subroutine
						    !in order to facilitate hard-coding and keep this subroutine clean.
						    !***********************************************************************
						    do i=1, numSpecies
							    product2(i)=0
						    end do

						    call defectCombinationRules(defectStore%defectType, product2,defectTemp, isCombined)

						    if(product2(3)/=0 .OR. product2(4)/=0) then
							    !defectStore is at the middle (or end) of defectStoreList
							    if(associated(defectStorePrev)) then
								    nullify(defectStorePrev%next)
								    allocate(defectStorePrev%next)
								    nullify(defectStorePrev%next%next)
								    defectStorePrev=>defectStorePrev%next
								    allocate(defectStorePrev%defectType(numSpecies))
								    defectStorePrev%cellNumber=defectStore%cellNumber
								    defectStorePrev%num=1
								    do i=1, numSpecies
									    defectStorePrev%defectType(i)=product2(i)
								    end do
								    defectStorePrev%next=>defectStore
							    else	!defectStore is at the beginning of defectStoreList
								    defectStorePrev=>defectStoreList
								    nullify(defectStorePrev%next)
								    allocate(defectStorePrev%next)
								    nullify(defectStorePrev%next%next)
								    defectStorePrev=>defectStorePrev%next
								    allocate(defectStorePrev%defectType(numSpecies))
								    defectStorePrev%cellNumber=defectStore%cellNumber
								    defectStorePrev%num=1
								    do i=1, numSpecies
									    defectStorePrev%defectType(i)=product2(i)
								    end do
								    defectStorePrev%next=>defectStore
							    end if
						    end if

						    if(isCombined .eqv. .TRUE.	) then
                                mixingEvents=mixingEvents+1
                            end if

                        end if
					end do

                    !NEXT: remove DefectTemp from coarse mesh
					if(defectTemp%num > 0) then
						if(mixingEvents < defectTemp%num) then
							defectTemp%num=defectTemp%num-mixingEvents
						else
							!if defectTemp is in the middle of the list
							if(associated(defectPrev) .AND. associated(defectTemp%next)) then
								defectPrev%next=>defectTemp%next
								deallocate(defectTemp%defectType)
								deallocate(defectTemp)
								defectTemp=>defectPrev
								!if defectTemp is at the end of the list
							else if(associated(defectPrev)) then
								deallocate(defectTemp%defectType)
								deallocate(defectTemp)
								defectTemp=>defectPrev
								nullify(defectPrev%next)
								!if defectTemp is at the beginning of the list
							else if(associated(defectTemp%next)) then
								defectTemp%num=0
							end if
						end if
					end if

					defectPrev=>defectTemp
					defectTemp=>defectTemp%next
					
				end do
			end if
			
			defectStore=>defectStoreList%next

			!add all cascade mixing products and normal cascade products to the system
			do while(associated(defectStore))
			
				defectCurrent=>defectList(defectStore%cellNumber)
				!STEP 2: add DefectTemp+products to the system
				count=0
				do j=1,numSpecies
					if(defectStore%defectType(j)==0) then
						count=count+1
					endif
				end do
				if(count==numSpecies) then
					!do nothing, we have total annihilation
				else
					! Updatin defect list with defects contained in cascade defect list
					nullify(defectPrev)
					call findDefectInList(defectCurrent, defectPrev, defectStore%defectType)
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						
						count=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==defectStore%defectType(j)) then
								count=count+1
							endif
						end do
						if(count==numSpecies) then
							defectCurrent%num=defectCurrent%num+1
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectStore%cellNumber
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=defectStore%defectType(j)
							end do
							!if inserted defect is in the middle of the list, point it to the next item in the list
							defectPrev%next=>defectCurrent
						endif
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectStore%cellNumber
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=defectStore%defectType(j)
						end do
					endif
					
				endif
				defectStore=>defectStore%next
			end do
			
			!Final step: set up defectUpdate for all defects in the cell
			defectTemp=>defectList(reactionCurrent%cellNumber(1))%next
			
			do while(associated(defectTemp))
			
				!*******************************************************************************
				!create a new element in defectUpdate and assign all variables
				!
				!Note: this will be done for each member of the cascade, even if more than 
				!one defect of this type is added to this volume element in the cell.
				!
				!This means that the reaction list update will happen multiple times (redundant)
				!for some defects in some volume elements in the cell.
				!*******************************************************************************
				
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=reactionCurrent%taskid(1)
				defectUpdateCurrent%dir=0	!not pointed at a different proc
				allocate(defectUpdateCurrent%defectType(numSpecies))
				do j=1,numSpecies
					defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
				end do
				defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(1)
				defectUpdateCurrent%neighbor=0	!not pointed at a different proc
				defectUpdateCurrent%cascadeNumber=0	!identify which cascade this defect is inside: 0 for coarse mesh
				defectUpdateCurrent%num=1
				nullify(defectUpdateCurrent%next)
				
				defectTemp=>defectTemp%next
			
			end do

			!*******************************************************************
			!memory erase: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
			!*******************************************************************
			defectStore=>defectStoreList
			
			do while(associated(defectStore))
				defectTemp=>defectStore
				defectStore=>defectStore%next
				deallocate(defectTemp%defectType)
				deallocate(defectTemp)
			end do
			
			nullify(defectStore)
			nullify(defectTemp)
			nullify(defectStoreList)
			nullify(defectCurrent)
			nullify(defectPrev)	
			
		else
			write(*,*) 'error meshingType'
		end if

	!***********************************************************************************************
	!Defect update for reactions within the fine mesh.
	!
	!Typically, no communication will be carried out if a reaction is chosen within the fine mesh.
	!The exception is when a defect diffuses from the fine mesh to the coarse mesh.
	!***********************************************************************************************
		
	else if(associated(CascadeCurrent)) then    !Reactions in the fine mesh

        allocate(firstSend(numSpecies+3,reactionCurrent%numProducts,6))
        !***********************************************************************
        !Here, I will assume cubic cells. If a defect migrates from one cell to another,
        !the defect will be given a percent chance of removal from the system
        !based on the cell size (%chance of removal = cell size/grain size)
        !This is to replicate the OKMC practice of removing SIA clusters after they
        !migrate 1 um.
        !***********************************************************************
        flag=.FALSE.
        if(grainBoundaryToggle=='yes') then

            !we have toggled the use of grain boundaries to remove defects from system

            if(reactionCurrent%numReactants==1 .and. reactionCurrent%numProducts==1) then
                diffusionRandom=dprand()

                !randomly choose whether to remove this defect from the system according to the mean free path (parameter) and the
                !length of the volume element that the defect is currently in
                if(diffusionRandom <= fineLength/meanFreePath) then
                    flag=.TRUE.
                    !reactionCurrent%numProducts=0 !remove that defect from system
                    !we can no longer remove the defect in this way because changing the reaction list
                    !is no longer allowed (reaction list is not re-created at each step)
                end if
            end if
        end if

        !Update buffer if needed (only relevant for diffusion from fine to coarse mesh and the coarse mesh is the boundry mesh of neighbor proc)
		!local defects
		if(flag .eqv. .FALSE.) then
			do i=1, reactionCurrent%numProducts
				if(reactionCurrent%taskid(i+reactionCurrent%numReactants) /= -1 .AND. &
						reactionCurrent%cellNumber(i+reactionCurrent%numReactants)==0) then
					do j=1,6	!directions
						if(myMesh(cascadeCurrent%cellNumber)%neighborProcs(1,j) /= myProc%taskid .AND. &
								myMesh(cascadeCurrent%cellNumber)%neighborProcs(1,j)/=-1) then

							!neighboring element not in this proc
							numUpdateSend(j)=numUpdateSend(j)+1

							!species
							do l=1,numSpecies
								firstSend(l,numUpdateSend(j),j)=reactionCurrent%products(l,i)
							end do

							!cell number in local mesh
							firstSend(numSpecies+1,numUpdateSend(j),j)=cascadeCurrent%cellNumber

							!number of defects to be increased by 1
							firstSend(numSpecies+2,numUpdateSend(j),j)=1

							!cell number of neighbor (in different proc)
							firstSend(numSpecies+3,numUpdateSend(j),j)=myMesh(cascadeCurrent%cellNumber)%neighbors(1,j)
						end if
					end do
				end if
			end do
		end if

        !Communication: send products to neighbor proc and receive defects (local and boundry) from neighbor proc
		do i=1,6

			if(mod(i,2)==0) then
				recvDir=i-1
			else
				recvDir=i+1
			end if

			!Send
			if(myProc%procNeighbor(i) /= myProc%taskid) then
				call MPI_ISEND(firstSend(1,1,i),(numSpecies+3)*numUpdateSend(i),MPI_INTEGER,myProc%procNeighbor(i), &
							200+i ,comm,sendFirstRequest(i),ierr)
			end if

			!Recv
			if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

				count=0
				call MPI_PROBE(myProc%procNeighbor(recvDir), 200+i,comm,status,ierr)
				call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

				numUpdateRecv(recvDir) = count/(numSpecies+3)

				call MPI_IRECV(firstRecv(1,1,recvDir),(numSpecies+3)*numUpdateRecv(recvDir),MPI_INTEGER,&
						myProc%procNeighbor(recvDir),200+i,comm,recvFirstRequest(recvDir),ierr)
			end if
		end do

		!Update defects in local mesh
        !Remove reactants from the system
		do i=1, reactionCurrent%numReactants
			
			!create a new element in defectUpdate and assign all variables except for num (will do later)
			allocate(defectUpdateCurrent%next)
			defectUpdateCurrent=>defectUpdateCurrent%next
			defectUpdateCurrent%proc=reactionCurrent%taskid(i)
			defectUpdateCurrent%dir=0	!not pointed at a different proc
			allocate(defectUpdateCurrent%defectType(numSpecies))
			
			do j=1,numSpecies
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(j,i)
			end do
			
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)	!the cell is the fineMeshID
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID	
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i))
			
			do while(associated(defectCurrent)) !search for reactant defect
				same=0
				do j=1,numSpecies
					if(defectCurrent%defectType(j)==reactionCurrent%reactants(j,i)) then
						same=same+1
					endif
				end do
				if(same==numSpecies) then
					exit
				end if
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
			end do
			
			!No need to update local buffer - cascade is always totally contained within the coarse mesh (one processor)
			
			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there fine'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'rate', reactionCurrent%reactionRate
				write(*,*) 'cascade number', cascadeCurrent%cascadeID

				call MPI_ABORT(comm,ierr)
			
			!if there is one defect of this type and it is in the middle of the list, remove it from the list
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next) .AND. associated(defectPrev)) then
				defectPrev%next=>defectCurrent%next !remove that defect type from the system
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)
				defectUpdatecurrent%num=-1			!tell updateReactionList that there are none of this defect left
			
			!if there is one defect of this type and it is at the end of the list then remove it from the list and
			!point defectPrev%next at nothing (nullify)
			else if(defectCurrent%num==1 .AND. associated(defectPrev)) then
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)	!remove the last defect from the system
				nullify(defectPrev%next)
				defectUpdateCurrent%num=-1			!tell updateReactionList that there are none of this defect left
				
			!if there is one defect of this type and it is at the beginning of the list then just make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next)) then !removing first defect from cell i
			!	defectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
				CascadeCurrent%localDefects(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then 							!removing only defect from cell i (single helium) - this is redundant but will keep for now
			!	defectList(reactionCurrent%cellNumber(i))%num=0
				CascadeCurrent%localDefects(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) 'trying to remove a defect that isnt there'
				call MPI_ABORT(comm,ierr)
			else
				!decrease the number of defects by 1 if the number of defects is greater than 1
				defectCurrent%num=defectCurrent%num-1 !remove the defect from the system instead of the entire entry in the list
				defectUpdateCurrent%num=-1	!tell updateReactionList the new number of defects in the cell
			end if
		end do

		!adding products to the system
		if(flag .eqv. .FALSE.) then	!if the diffusion reaction defect did not get removed from the system
			do i=1, reactionCurrent%numProducts
				
				if(reactionCurrent%taskid(i+reactionCurrent%numReactants) == -1) then
					!we have identified a defect that is going to be removed from the system via a free surface
					!Do nothing; no need to add this defect to any lists or update any reaction lists because of it
					!This should not happen because defects in fine meshes should not be able to pass through free surfaces
					write(*,*) 'Error free surface in fine mesh reaction'
				else
				!Update defects and buffer if needed (only relevant for diffusion from fine to coarse mesh)
				
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
					allocate(defectUpdateCurrent%defectType(numSpecies))
					
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					
					nullify(defectPrev)
					do j=1,numSpecies
						products(j)=reactionCurrent%products(j,i)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
					end do
					
					if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants)==0) then
						
						!diffusion from fine mesh to coarse mesh: point defectCurrent at coarse mesh defect list
						defectCurrent=>defectList(CascadeCurrent%cellNumber)
						
						!defectUpdateCurrent tells updateReactionList to update reactions within the coarse mesh
						defectUpdateCurrent%cascadeNumber=0
						defectUpdateCurrent%cellNumber=CascadeCurrent%cellNumber
					else
						
						!reaction within the fine mesh: point defectCurrent at fine mesh defect list
						defectCurrent=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
						
						!defectUpdateCurrent tells updateReactionList to update reactions within fine mesh
						defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID
						defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					end if
	
					! this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)

					!*********************************************
					!Update defects
					!*********************************************
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						end do
						
						if(same==numSpecies) then
							!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+1
						else
							!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectUpdateCurrent%cellNumber
							defectPrev%num=1
							
							do j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(j,i)
							end do
							
							!if inserted defect is in the middle of the list, point it to the next item in the list
							defectPrev%next=>defectCurrent 
						end if
					else
						!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectUpdateCurrent%cellNumber
						defectPrev%num=1
						
						do j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(j,i)
						end do

					end if
				end if
			end do
		else
			!If the defect has been removed due to grain boundary absorption, do nothing
		end if
	
	!***********************************************************************************************
	! Defect update for reactions chosen in the coarse mesh.
	!
	! For defect changes in the boundary of other processors or in the boundary of this processor,
	! local and boundary buffers are created for communication with other processors.
	!***********************************************************************************************

	else	!Reactions in the coarse mesh
        allocate(firstSend(numSpecies+3,reactionCurrent%numReactants+reactionCurrent%numProducts,6))
		!First update local buffer if needed
		do i=1, reactionCurrent%numReactants

			nullify(defectPrev)
			defectCurrent=>defectList(reactionCurrent%cellNumber(i))

			!Check that defectCurrent is pointing towards the reactant
			do while(associated(defectCurrent)) !search for reactant defect
				same=0
				do j=1,numSpecies
					if(defectCurrent%defectType(j)==reactionCurrent%reactants(j,i)) then
						same=same+1
					end if
				end do
				if(same==numSpecies) then
					exit
				end if
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
			end do

			do j=1,6
				if(myMesh(reactionCurrent%cellNumber(i))%neighborProcs(1,j) /= myProc%taskid .AND. &
						myMesh(reactionCurrent%cellNumber(i))%neighborProcs(1,j) /= -1) then	!neighboring element not in this proc

					numUpdateSend(j)=numUpdateSend(j)+1

					do l=1,numSpecies
						firstSend(l,numUpdateSend(j),j)=reactionCurrent%reactants(l,i)
					end do

					!Cell Number in local mesh
					firstSend(numSpecies+1,numUpdateSend(j),j)=reactionCurrent%cellNumber(i)

					!indetify the num of defects to be updated
					if(associated(defectCurrent)) then
						firstSend(numSpecies+2,numUpdateSend(j),j)=-1
					else
						write(*,*) 'Error tried to delete defect that wasnt there and send to neighboring proc'
						write(*,*) 'step', step, 'proc',myProc%taskid,'i',i, 'j',j,&
								'cellNumber', reactionCurrent%cellNumber, 'taskid', reactionCurrent%taskid
						write(*,*) 'numReact',reactionCurrent%numReactants,'numProc', reactionCurrent%numProducts,&
								'reacts', reactionCurrent%reactants, 'products', reactionCurrent%products
						defectTest=>defectList(reactionCurrent%cellNumber(i))
						do while(associated(defectTest))
							write(*,*) defectTest%defectType, defectTest%num, defectTest%cellNumber
							defectTest=>defectTest%next
						end do
					end if

					!Cell number in boundary mesh
					firstSend(numSpecies+3,numUpdateSend(j),j)=myMesh(reactionCurrent%cellNumber(i))%neighbors(1,j)

				end if
			end do
		end do

		!***********************************************************************
		!Here, I will assume cubic cells. If a defect migrates from one cell to another,
		!the defect will be given a percent chance of removal from the system
		!based on the cell size (%chance of removal = cell size/grain size)
		!This is to replicate the OKMC practice of removing SIA clusters after they
		!migrate 1 um.
		!***********************************************************************
		flag=.FALSE.
		if(grainBoundaryToggle=='yes') then

			!we have toggled the use of grain boundaries to remove defects from system
			!7/18/2014: NOTE that if we have diffusion from coarse mesh to fine mesh, the length
			!of the coarse mesh element will be used to calculate the probability of grain boundary
			!absorption, instead of the effective coarse to fine length. This could be changed if desired.

			if(reactionCurrent%numReactants==1 .and. reactionCurrent%numProducts==1) then
				diffusionRandom=dprand()

				!randomly choose whether to remove this defect from the system according to the mean free path (parameter) and the
				!length of the volume element that the defect is currently in
				if(diffusionRandom <= myMesh(reactionCurrent%cellNumber(1))%length/meanFreePath) then
					flag=.TRUE.
					!reactionCurrent%numProducts=0 !remove that defect from system
					!we can no longer remove the defect in this way because changing the reaction list
					!is no longer allowed (reaction list is not re-created at each step)
				end if
			end if
		end if
		!write(*,*) 'step', step, 'proc', myProc%taskid, 'numReact',reactionCurrent%numReactants, &
		!		'numProduct',reactionCurrent%numProducts, 'reactantsCell', reactionCurrent%cellNumber
        if(flag .eqv. .FALSE.) then
            do i=1, reactionCurrent%numProducts
				if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants) < 0) then	!diffuse to fine mesh
					!do nothing
				else

				!diffuse to neighboring processor
                if(reactionCurrent%taskid(i+reactionCurrent%numReactants) /= myProc%taskid .AND. &
                        reactionCurrent%taskid(i+reactionCurrent%numReactants) /= -1) then

                    !we have identified a defect that needs to be added to a neighboring processor and the boundary mesh
					flagTemp=.FALSE.
                    do j=1,6
						if(flagTemp .EQV. .FALSE.) then
                             if(myMesh(reactionCurrent%cellNumber(1))%neighbors(1,j)==&
                                     reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
                                     myMesh(reactionCurrent%cellNumber(1))%neighborProcs(1,j)==&
                                             reactionCurrent%taskid(i+reactionCurrent%numReactants)) then

								 numUpdateSend(j)=numUpdateSend(j)+1

                                 do l=1,numSpecies
                                     firstSend(l,numUpdateSend(j),j)=reactionCurrent%products(l,i)
                                 end do

                                 firstSend(numSpecies+1,numUpdateSend(j),j)=&
                                         reactionCurrent%cellNumber(i+reactionCurrent%numReactants)

                                 firstSend(numSpecies+2,numUpdateSend(j),j)=1	!indicates we are adding one defect

								 firstSend(numSpecies+3,numUpdateSend(j),j)=-100	!indicates the defect is a boundry defect
								 flagTemp=.TRUE.
                             end if
						end if
                    end do
                else
                    !First update local buffer if needed
                    do j=1,6
                        if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(1,j) &
                                /= myProc%taskid .AND. &
								myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(1,j)/=-1) then	!neighboring element not in this proc

							numUpdateSend(j)=numUpdateSend(j)+1

                            !species
                            do l=1,numSpecies
                                firstSend(l,numUpdateSend(j),j)=reactionCurrent%products(l,i)
                            end do

                            !cell number in local mesh
                            firstSend(numSpecies+1,numUpdateSend(j),j)=&
                                    reactionCurrent%cellNumber(i+reactionCurrent%numReactants)

							!number of defects to be increased by 1
							firstSend(numSpecies+2,numUpdateSend(j),j)=1

                            !cell number of neighbor (in different proc)
                            firstSend(numSpecies+3,numUpdateSend(j),j)=&
                                    myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighbors(1,j)
                        end if
                    end do
                end if

				end if
            end do
        end if

        !Communication
        do i=1,6

            if(mod(i,2)==0) then
                recvDir=i-1
            else
                recvDir=i+1
            end if

			!Send
			if(myProc%procNeighbor(i) /= myProc%taskid) then
				call MPI_ISEND(firstSend(1,1,i),(numSpecies+3)*numUpdateSend(i),MPI_INTEGER,myProc%procNeighbor(i),&
						200+i ,comm,sendFirstRequest(i),ierr)
			end if

			!Recv
			if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

				count=0
				call MPI_PROBE(myProc%procNeighbor(recvDir), 200+i,comm,status,ierr)
				call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

				numUpdateRecv(recvDir) = count/(numSpecies+3)
				!totalLocalRecv=totalLocalRecv+numUpdateRecv(recvDir)

				call MPI_IRECV(firstRecv(1,1,recvDir),(numSpecies+3)*numUpdateRecv(recvDir),MPI_INTEGER,&
						myProc%procNeighbor(recvDir),200+i ,comm,recvFirstRequest(recvDir),ierr)
			end if
        end do

		!removing reactants from the system
		do i=1, reactionCurrent%numReactants
			
			!create a new element in defectUpdate and assign all variables except for num (will do later)
			allocate(defectUpdateCurrent%next)
			defectUpdateCurrent=>defectUpdateCurrent%next
			defectUpdateCurrent%proc=reactionCurrent%taskid(i)
			defectUpdateCurrent%dir=0	!not pointed at a different proc
			allocate(defectUpdateCurrent%defectType(numSpecies))
			do j=1,numSpecies
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(j,i)
			end do
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>defectList(reactionCurrent%cellNumber(i))

			!Check that defectCurrent is pointing towards the reactant
			do while(associated(defectCurrent)) !search for reactant defect
				same=0
				do j=1,numSpecies
					if(defectCurrent%defectType(j)==reactionCurrent%reactants(j,i)) then
						same=same+1
					endif
				end do
				if(same==numSpecies) then
					exit
				end if
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
			end do

			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there coarse'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'step',step,'proc',myProc%taskid,'i',i
				write(*,*) 'rate', reactionCurrent%reactionRate
				call MPI_ABORT(comm,ierr)
			
			!if there is one defect of this type and it is in the middle of the list, remove it from the list
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next) .AND. associated(defectPrev)) then
				defectPrev%next=>defectCurrent%next !remove that defect type from the system
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)
				defectUpdatecurrent%num=-1			!tell updateReactionList that there are none of this defect left
			
			!if there is one defect of this type and it is at the end of the list then remove it from the list and
			!point defectPrev%next at nothing (nullify)
			else if(defectCurrent%num==1 .AND. associated(defectPrev)) then
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)	!remove the last defect from the system
				nullify(defectPrev%next)
				defectUpdateCurrent%num=-1			!tell updateReactionList that there are none of this defect left
				
			!if there is one defect of this type and it is at the beginning of the list then just make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next)) then
				defectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates
				defectUpdateCurrent%num=-1
				write(*,*) 'Removing the first defect from the list'
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then
				defectList(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
				write(*,*) 'Removing the only defect from the list'
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) '&&&&trying to remove a defect that isnt there'
				write(*,*) 'step', step,'proc',myProc%taskid,'cell',reactionCurrent%cellNumber(i)
				call MPI_ABORT(comm,ierr)
			else
				!decrease the number of defects by 1 if the number of defects is greater than 1
				!write(86,*) 'decreasing defect num by 1', defectCurrent%num, 'type', defectCurrent%defectType, &
				!	'cell', defectCurrent%cellNumber
				!write(86,*)
				!if(associated(defectPrev))	write(86,*) 'defect', defectPrev%defectType, 'num', defectPrev%num
				!write(86,*) 'defect', defectCurrent%defectType, 'num', defectCurrent%num
				!if(associated(defectCurrent%next)) !write(86,*) 'defect', defectCurrent%next%defectType, &
				!	'num', defectCurrent%next%num
				
				defectCurrent%num=defectCurrent%num-1 !remove the defect from the system instead of the entire entry in the list
				
				if(defectCurrent%num==0) then
					!Error: we should have deleted this defect.
					write(*,*) 'Error did not delete defects when number of defects present equal to zero'
				endif
				
				defectUpdateCurrent%num=-1	!tell updateReactionList the new number of defects in the cell
			endif
		end do

		!adding products to the system
		if(flag .eqv. .FALSE.) then		!if the diffusion reaction defect did not get removed from the system
			
			!***********************************************************************************************
			!EDIT: 5/27/2014: Included the possibility of updating defects in neighboring processors
			!(the information should be stored in reaction%taskid(:))
			!
			!Procedure: if taskid .NE. myProc%taskid, then add to a buffer the number of the proc to send to,
			!the information on what defects are being added and in what volume elements
			!
			!Then send a flag with the number of defects being sent/recieve a flag with number of defects being
			!recieved
			!
			!Then send/recieve information on defects being added to system and add them to the correct volume
			!elements.
			!
			!RE-EDIT: 5/28/2014: Included the boundary mesh: this includes the updating of defects in neighboring
			!processors as well.
			!
			!EDIT: 7/18/2014: Added diffusion from coarse to fine mesh: if reactionCurrent%cellNumber .LT. 0
			!***********************************************************************************************
			
			do i=1, reactionCurrent%numProducts
				
				if(reactionCurrent%taskid(i+reactionCurrent%numReactants) == -1) then
					!we have identified a defect that is going to be removed from the system via a free surface
					!Do nothing; no need to add this defect to any lists or update any reaction lists because of it
					
				else if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants) < 0) then	!diffusion: coarse-to-fine
					!we have identified a defect that is diffusing from the coarse mesh to the fine mesh
					!and therefore neet to treat it differently than diffusion within the coarse mesh.
					!
					!In this case, we do not update the buffers because coarse - to - fine reactions
					!are within the same processor and fine meshes are not within the boundary of 
					!neighboring processors.
					
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
					allocate(defectUpdateCurrent%defectType(numSpecies))
					
					!When coarse-to-fine reactions are chosen, a random cell within the fine mesh
					!is chosen.
					defectUpdateCurrent%cellNumber=chooseRandomCell()	!chooseRandomCell()>0
					
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					
					!In coarse-to-fine reactions, the cascade number is stored in reactionCurrent%cellNumber
					!and identified with a negative. Therefore we make it positive again.
					defectUpdateCurrent%cascadeNumber=-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					
					nullify(defectPrev)
					do j=1,numSpecies
						products(j)=reactionCurrent%products(j,i)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
					end do
					
					!Point CascadeCurrent at the correct cascade
					
					CascadeCurrent=>ActiveCascades
					
					do while(associated(CascadeCurrent))
						
						if(CascadeCurrent%cascadeID==-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)) then
							exit
						endif
						
						CascadeCurrent=>CascadeCurrent%next
					
					end do
					
					!If we went through the previous loop without pointing at a cascade, then we have
					!a problem with placing the coarse-to-fine miration in the correct cascade.
					
					if(.NOT. associated(CascadeCurrent)) then
						write(*,*) 'error coarse-to-fine cascade not correctly pointed'
					endif
					
					!Point defectCurrent at the defect list within the correct element in CascadeCurrent
					defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
					
					!*******************************************
					!Update defects within fine mesh
					!*******************************************

					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							end if
						end do
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+1
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectUpdateCurrent%cellNumber
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(j,i)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectUpdateCurrent%cellNumber
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(j,i)
						end do
					end if

				!Diffusion to neighbor processes
				else if(reactionCurrent%taskid(i+reactionCurrent%numReactants) /= myProc%taskid .AND. &
					reactionCurrent%taskid(i+reactionCurrent%numReactants) /= -1) then
					
					!we have identified a defect that needs to be added to a neighboring processor and the boundary mesh
					
					!***********************************************************************************
					!EDIT: 6/4/2014: we only want to add to the boundary buffer once in this step
					!because this will be used to change only one local defect. However, if the number
					!of processors is large, then one mesh can have a processor bounding it on both sides
					!(assuming periodic boundary conditions). Therefore, we have adjusted the code below
					!to only add to bndrySend once even though we change myBoundary more than once (if
					!the same cell bounds this mesh in two places due to PBCs)
					!
					!This is accomplished using the logical flag. Once the defect is added to the bounary
					!buffer in one direction, it cannot occur again and therefore we won't double-count
					!defect migration.
					!***********************************************************************************
					
					!step 1: find the element in the boundary mesh. The element may bound this processor at more than one place.
					!**************************************************************************************************
					!EDIT: 6/4/2014: this may not be the fastest way to search for which cells in myBoundary to update.
					!EDIT: 6/5/2014: Updated the way that we search through the local mesh for boundary elements which
					!correspond to the one in our reaction. Added %localNeighbor to the boundary class so that we know
					!exactly what the local element is whose neighbor is in the boundary. Now, we know the boundary element
					!number so we just search in all 6 directions for a boundary element that is allocated and that 
					!has a local neighbor in the local mesh. (Avoids searching through all local elements for one
					!with a boundary element corresponding to the boundary element we are looking for)
					!NOTE: old code commented out below
					!NOTE: this both increased the speed of computation and reduced the number of steps. This is surprising
					!as we are not changing anything physical within the system, but I have not found a bug and the results
					!match expected results so far. Will continue to keep an eye out for a bug in this routine.
					!**************************************************************************************************

					!update the direction in defectUpdateCurrent (used for updating reaction lists)
					!create a new element in defectUpdate and assign all variables except for num and dir (will do later)
					do j=1,6
	
						!if myBoundary at this cell has defectList associated, then it is a boundary element for a local element in this cell.
						!NOTE: this will not work in the case of a non-uniform mesh (because there could be more
						!than one local neighbor, will have to change to localNeighbor())
						if(myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%localNeighbor /= 0 &
							.AND. myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%proc == &
							reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
							
							!Create next entry in defectUpdate (used for updating reaction lists)
							allocate(defectUpdateCurrent%next)
							defectUpdateCurrent=>defectUpdateCurrent%next
							defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
							allocate(defectUpdateCurrent%defectType(numSpecies))
							defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							nullify(defectUpdateCurrent%next)
							defectUpdateCurrent%dir=j
							defectUpdateCurrent%num=1	!used for updating reaction lists
							defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
							
							!This information is used to tell us which element in the local
							!processor bounds this element in myBoundary
							defectUpdateCurrent%neighbor=myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)&
								%localNeighbor		
							
							do l=1,numSpecies
								products(l)=reactionCurrent%products(l,i)
								defectUpdateCurrent%defectType(l)=reactionCurrent%products(l,i)
							end do
							
							!Find defect in defect list on myBoundary
							nullify(defectPrev)
							defectCurrent=>myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%defectList
	
							call findDefectInList(defectCurrent, defectPrev, products)
	
							!*****************************************************************
							!Update defects within myBoundary
							!UpdateDefectList: Insert defect in myBoundary(dir,cell)%defectList
							!*****************************************************************

							if(associated(defectCurrent)) then !if we aren't at the end of the list
								same=0
								do l=1,numSpecies
									if(defectCurrent%defectType(l)==products(l)) then
										same=same+1
									endif
								end do
								if(same==numSpecies) then	!if the defect is already present in the list
									defectCurrent%num=defectCurrent%num+1
								else		!if the defect is to be inserted in the list
									nullify(defectPrev%next)
									allocate(defectPrev%next)
									nullify(defectPrev%next%next)
									defectPrev=>defectPrev%next
									allocate(defectPrev%defectType(numSpecies))
									defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
									defectPrev%num=1
									do l=1,numSpecies
										defectPrev%defectType(l)=reactionCurrent%products(l,i)
									end do
									defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
								endif
							else 			!add a defect to the end of the list
								nullify(defectPrev%next)
								allocate(defectPrev%next)
								nullify(defectPrev%next%next)
								defectPrev=>defectPrev%next
								allocate(defectPrev%defectType(numSpecies))
								defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
								defectPrev%num=1
								do l=1,numSpecies
									defectPrev%defectType(l)=reactionCurrent%products(l,i)
								end do
							end if
						end if
					end do
				else
				
				!For reactions that have occurred within this processor in the coarse mesh
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
					allocate(defectUpdateCurrent%defectType(numSpecies))
					defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
					
					nullify(defectPrev)
					do j=1,numSpecies
						products(j)=reactionCurrent%products(j,i)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
					end do
					
					defectCurrent=>defectList(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)

					!***************************************
					!Update defects within coarse mesh
					!***************************************
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						end do
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+1
						else		!if the defect is to be inserted in the list
							
							if(.NOT. associated(defectPrev)) then
								write(*,*) 'Error defectPrev not associated'
								write(*,*) defectCurrent%defectType
								write(*,*) products
							endif
							
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(j,i)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(j,i)
						end do
					end if
				end if
			end do
		else
			!If the defect has been removed due to grain boundary absorption
		end if
	end if
end if	!if associated(reactionCurrent)


do i=1,6

    if(mod(i,2)==0) then
        recvDir=i-1
    else
        recvDir=i+1
    end if

	if(myProc%procNeighbor(i) /= myProc%taskid) then
		call MPI_WAIT(sendFirstRequest(i),sendFirstStatus,ierr)
	end if

	if(myProc%procNeighbor(recvDir) /= myProc%taskid) then
		call MPI_WAIT(recvFirstRequest(recvDir),recvFirstStatus,ierr)
	end if
end do

!totalLocalRecv=1
allocate(finalSend(numSpecies+3,6,6))

!Update defectList
do i=1,6

	numBndryRecv=0
	numLocalRecv=0

	do j=1, numUpdateRecv(i)
		if(firstRecv(numSpecies+3,j,i) == -100) then
			numLocalRecv=numLocalRecv+1
		else
			numBndryRecv=numBndryRecv+1
		end if
	end do

    !Add defects in bndryRecv to myBoundary()
    do j=1,numBndryRecv

        !create a new element in defectUpdate and assign all variables except for num (will do later)
        allocate(defectUpdateCurrent%next)
        defectUpdateCurrent=>defectUpdateCurrent%next
        defectUpdateCurrent%proc=myProc%procNeighbor(i)
        defectUpdateCurrent%dir=i
        allocate(defectUpdateCurrent%defectType(numSpecies))
        defectUpdateCurrent%cellNumber=firstRecv(numSpecies+1,j,i)
        defectUpdateCurrent%num=firstRecv(numSpecies+2,j,i)	!This will be +/- 1 only
        defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
        nullify(defectUpdateCurrent%next)

        !This information is used to tell us which element in the local
        !processor bounds this element in myBoundary
        defectUpdateCurrent%neighbor=firstRecv(numSpecies+3,j,i)

        !point defectCurrent at the defect list in the correct cell of myBoundary
        defectCurrent=>myBoundary(firstRecv(numSpecies+1,j,i),i)%defectList
        if(.NOT. associated(defectCurrent)) then
            write(*,*) 'error myBoundary not allocated correctly'
            write(*,*) 'dir', i, 'cell', firstRecv(numSpecies+1,j,i)
            call MPI_ABORT(comm,ierr)
        end if

        do k=1,numSpecies
            products(k)=firstRecv(k,j,i)
            defectUpdateCurrent%defectType(k)=firstRecv(k,j,i)	!used to update reaction lists
        end do

        !point DefectCurrent at the defect we are looking for (if it is there), otherwise
        !point defectCurrent after and defectPrev before where it should go in defectList.
        nullify(defectPrev)
        call findDefectInList(defectCurrent, defectPrev, products)

        !Next update defects in myBoundary
        if(associated(defectCurrent)) then !if we aren't at the end of the list

            same=0
            do l=1,numSpecies
                if(defectCurrent%defectType(l)==products(l)) then
                    same=same+1
                end if
            end do
            if(same==numSpecies) then	!if the defect is already present in the list
                defectCurrent%num=defectCurrent%num+firstRecv(numSpecies+2,j,i) !this will be +/- 1 only

                !2019.05.26 Debug: if the defect of this type is at the end of the defectList
                !if there is one defect of this type and it is in the middle of the defectList
                if(defectCurrent%num==0 .AND. associated(defectCurrent%next) .AND. associated(defectPrev)) then
                    !delete this defect from the list in myBoundary
                    defectPrev%next=>defectCurrent%next !remove that defect type from the system
                    deallocate(defectCurrent%defectType)
                    deallocate(defectCurrent)
                    nullify(defectCurrent)
                    !if there is one defect of this type and it is 	at the end of the defectList
                else if(defectCurrent%num==0 .AND. associated(defectPrev)) then
                    deallocate(defectCurrent%defectType)
                    deallocate(defectCurrent)
                    nullify(defectCurrent)
                    nullify(defectPrev%next)
                end if
            else		!if the defect is to be inserted in the list

                if(firstRecv(numSpecies+2,j,i) == -1) then
                    write(*,*) 'defectCurrent associated', defectCurrent%defectType, 'num', defectCurrent%num
                    if(associated(defectPrev)) write(*,*) 'defectPrev', defectPrev%defectType, 'num', defectPrev%num
                    write(*,*) 'error in defectUpdate negative defect numbers', myProc%taskid
                    write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
                    write(*,*) (firstRecv(k,j,i),k=1,numSpecies+3)
					write(*,*) 'step',step,'numBndryRecv',numBndryRecv, 'j',j
					defectTemp=>myBoundary(firstRecv(numSpecies+1,j,i),i)%defectList
					do while(associated(defectTemp))
						write(*,*) defectTemp%defectType, defectTemp%num
						defectTemp=>defectTemp%next
					end do

                    call MPI_ABORT(comm,ierr)

                else
                    nullify(defectPrev%next)
                    allocate(defectPrev%next)
                    nullify(defectPrev%next%next)
                    defectPrev=>defectPrev%next
                    allocate(defectPrev%defectType(numSpecies))
                    defectPrev%cellNumber=firstRecv(numSpecies+1,j,i)
                    defectPrev%num=firstRecv(numSpecies+2,j,i)	!This will be +/- 1 only

                    do l=1,numSpecies
                        defectPrev%defectType(l)=firstRecv(l,j,i)
                    end do
                    defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
                end if
            end if
        else 			!add a defect to the end of the list
            if(firstRecv(numSpecies+2,j,i) == -1) then
                write(*,*) 'defectCurrent not associated'
                write(*,*) 'error in defectUpdate negative defect numbers'
                write(*,*) 'j',j,'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
                write(*,*) (firstRecv(k,j,i),k=1,numSpecies+3)
                !write(*,*) 'step',step
                !defectTest=>myBoundary(firstRecv(numSpecies+1,j,i),i)%defectList
                !do while(associated(defectTest))
                !    write(*,*) defectTest%defectType, defectTest%num, defectTest%cellNumber
                !    defectTest=>defectTest%next
                !end do
                call MPI_ABORT(comm,ierr)
            else

                nullify(defectPrev%next)
                allocate(defectPrev%next)
                nullify(defectPrev%next%next)
                defectPrev=>defectPrev%next
                allocate(defectPrev%defectType(numSpecies))
                defectPrev%cellNumber=firstRecv(numSpecies+1,j,i)
                defectPrev%num=firstRecv(numSpecies+2,j,i)	!This will be +/- 1 only

                do l=1,numSpecies
                    defectPrev%defectType(l)=firstRecv(l,j,i)
                end do
            end if
        end if
    end do

    !Add defects in localRecv to defectList()
    do j=numBndryRecv+1,numUpdateRecv(i)

        !create a new element in defectUpdate and assign all variables except for num (will do later)
        allocate(defectUpdateCurrent%next)
        defectUpdateCurrent=>defectUpdateCurrent%next
        defectUpdateCurrent%proc=myProc%taskid
        defectUpdateCurrent%dir=0	!not pointed at a different proc
        defectUpdateCurrent%neighbor=0	!not pointed at a different proc
        allocate(defectUpdateCurrent%defectType(numSpecies))
        defectUpdateCurrent%cellNumber=firstRecv(numSpecies+1,j,i)
        defectUpdateCurrent%num=firstRecv(numSpecies+2,j,i) !This will be +/- 1 only
        defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
        nullify(defectUpdateCurrent%next)

        !point defectCurrent at the defect list in the correct cell of defectList
        defectCurrent=>defectList(firstRecv(numSpecies+1,j,i))

        do k=1,numSpecies
            products(k)=firstRecv(k,j,i)
            defectUpdateCurrent%defectType(k)=firstRecv(k,j,i)
        end do

        !point DefectCurrent at the defect we are looking for (if it is there), otherwise
        !point defectCurrent after and defectPrev before where it should go in defectList.
        nullify(defectPrev)
        call findDefectInList(defectCurrent, defectPrev, products)

        !Next update defects in local mesh
        if(associated(defectCurrent)) then !if we aren't at the end of the list
            same=0
            do l=1,numSpecies
                if(defectCurrent%defectType(l)==products(l)) then
                    same=same+1
                endif
            end do
            if(same==numSpecies) then	!if the defect is already present in the list
                defectCurrent%num=defectCurrent%num+firstRecv(numSpecies+2,j,i)	!This will be +1 only

                if(defectCurrent%num==0) then
                    write(*,*) 'Error zero defects in updateDefectList, step 3 of communication'
                    !		defectPrev%next=>defectCurrent%next !remove that defect type from the system
                    !		deallocate(defectCurrent%defectType)
                    !		deallocate(defectCurrent)
                    !		nullify(defectCurrent)
                end if
            else		!if the defect is to be inserted in the list
                nullify(defectPrev%next)
                allocate(defectPrev%next)
                nullify(defectPrev%next%next)
                defectPrev=>defectPrev%next
                allocate(defectPrev%defectType(numSpecies))
                defectPrev%cellNumber=firstRecv(numSpecies+1,j,i)
                defectPrev%num=firstRecv(numSpecies+2,j,i) !This will be +/- 1 only (should be only +1 here)
                do l=1,numSpecies
                    defectPrev%defectType(l)=firstRecv(l,j,i)
                end do
                defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
            end if
        else 			!add a defect to the end of the list
            nullify(defectPrev%next)
            allocate(defectPrev%next)
            nullify(defectPrev%next%next)
            defectPrev=>defectPrev%next
            allocate(defectPrev%defectType(numSpecies))
            defectPrev%cellNumber=firstRecv(numSpecies+1,j,i)
            defectPrev%num=firstRecv(numSpecies+2,j,i) !This will be +/- 1 only (should be only +1 here)
            do l=1,numSpecies
                defectPrev%defectType(l)=firstRecv(l,j,i)
            end do
        end if

        !*****************************************************
        !Prep for step 4: check to see if element updated has more than one neighbor proc
        !*****************************************************

        do k=1,6

			!if a neighbor of this element is not in this proc or in the proc that just communicated with it
			if (myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(1,k) /= myProc%taskid .AND. &
					myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(1,k) /= myProc%procNeighbor(i) .AND. &
					myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(1,k) /= -1) then

				!Add this defect to a final buffer
				numUpdateFinal(k)=numUpdateFinal(k)+1

				do m=1,numSpecies
					finalSend(m,numUpdateFinal(k),k)=firstRecv(m,j,i)
				end do
				finalSend(numSpecies+1,numUpdateFinal(k),k)=firstRecv(numSpecies+1,j,i)	!local cell number
				finalSend(numSpecies+2,numUpdateFinal(k),k)=firstRecv(numSpecies+2,j,i) !This should be +1 only
				finalSend(numSpecies+3,numUpdateFinal(k),k)=myMesh(firstRecv(numSpecies+1,j,i))%neighbors(1,k)	!cell number in neighboring proc

				if(firstRecv(numSpecies+2,j,i)==-1) then
					write(*,*) 'error: cell in boundary of multiple procs removing defect'
					call MPI_ABORT(comm,ierr)
				end if
			end if
        end do
    end do
end do

if(allocated(firstRecv)) then
	deallocate(firstRecv)
end if

if(allocated(firstSend)) then
    deallocate(firstSend)
end if

!*************
!Step 4: if a local defect is updated due to diffusion to the boundary of another processor,
!	and the local element is also in the boundary of a processor other than the one already 
!	communicated with, let all other processors know of local defects that have been updated
!*************

do i=1,6

    if(mod(i,2)==0) then
        recvDir = i-1
    else
        recvDir = i+1
    end if

    !Send
    if(myProc%procNeighbor(i) /= myProc%taskid) then
		call MPI_SEND(finalSend(1,1,i),(numSpecies+3)*numUpdateFinal(i),MPI_INTEGER,&
				myProc%procNeighbor(i), 400+i ,comm,ierr)
        !call MPI_ISEND(finalSend(1,1,i),(numSpecies+3)*numUpdateFinal(i),MPI_INTEGER,&
        !        myProc%procNeighbor(i), 400+i ,comm,sendFinalRequest(i),ierr)
    end if

    !Recv
    if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

        !Recv Bndry
        count=0
        call MPI_PROBE(myProc%procNeighbor(recvDir), 400+i,comm,status,ierr)
        call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

        numUpdateFinalRecv(recvDir) = count/(numSpecies+3)

        allocate(finalBufferRecv(numSpecies+3,numUpdateFinalRecv(recvDir)))

		call MPI_RECV(finalBufferRecv,(numSpecies+3)*numUpdateFinalRecv(recvDir),MPI_INTEGER,&
				myProc%procNeighbor(recvDir),400+i ,comm,status,ierr)

        !Add defects in finalBufferRecv to myBoundary()
        do j=1,numUpdateFinalRecv(recvDir)

            if(finalBufferRecv(numSpecies+2,j)==-1) then
                write(*,*) 'error deleting defects in finalBufferRecv'
                !call MPI_ABORT(comm,ierr)
            else

                !create a new element in defectUpdate and assign all variables except for num (will do later)
                allocate(defectUpdateCurrent%next)
                defectUpdateCurrent=>defectUpdateCurrent%next
                defectUpdateCurrent%proc=myProc%procNeighbor(recvDir)
                defectUpdateCurrent%dir=recvDir
                allocate(defectUpdateCurrent%defectType(numSpecies))
                defectUpdateCurrent%cellNumber=finalBufferRecv(numSpecies+1,j)	!cell number in boundary mesh
                defectUpdateCurrent%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
                defectUpdateCurrent%neighbor=finalBufferRecv(numSpecies+3,j) !cell number in myMesh
                defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
                nullify(defectUpdateCurrent%next)

                !point defectCurrent at the defect list in the correct cell of myBoundary
                defectCurrent=>myBoundary(finalBufferRecv(numSpecies+1,j),recvDir)%defectList

                do k=1,numSpecies
                    products(k)=finalBufferRecv(k,j)
                    defectUpdateCurrent%defectType(k)=finalBufferRecv(k,j)
                end do

                !point DefectCurrent at the defect we are looking for (if it is there), otherwise
                !point defectCurrent after and defectPrev before where it should go in defectList.
				nullify(defectPrev)
                call findDefectInList(defectCurrent, defectPrev, products)

                !Next update defects in myBoundary
                if(associated(defectCurrent)) then !if we aren't at the end of the list
                    same=0
                    do l=1,numSpecies
                        if(defectCurrent%defectType(l)==products(l)) then
                            same=same+1
                        end if
                    end do
                    if(same==numSpecies) then	!if the defect is already present in the list
                        defectCurrent%num=defectCurrent%num+finalBufferRecv(numSpecies+2,j) !This should be +1 only

                        if(defectCurrent%num==0) then
                            write(*,*) 'Error zero defects in updateDefectList, step 5 of communication'
                        end if
                    else		!if the defect is to be inserted in the list
                        nullify(defectPrev%next)
                        allocate(defectPrev%next)
                        nullify(defectPrev%next%next)
                        defectPrev=>defectPrev%next
                        allocate(defectPrev%defectType(numSpecies))
                        defectPrev%cellNumber=finalBufferRecv(numSpecies+1,j)
                        defectPrev%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
                        do l=1,numSpecies
                            defectPrev%defectType(l)=finalBufferRecv(l,j)
                        end do
                        defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
                    end if
                else 			!add a defect to the end of the list
                    nullify(defectPrev%next)
                    allocate(defectPrev%next)
                    nullify(defectPrev%next%next)
                    defectPrev=>defectPrev%next
                    allocate(defectPrev%defectType(numSpecies))
                    defectPrev%cellNumber=finalBufferRecv(numSpecies+1,j)
                    defectPrev%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
                    do l=1,numSpecies
                        defectPrev%defectType(l)=finalBufferRecv(l,j)
                    end do
                end if
            end if
        end do
        deallocate(finalBufferRecv)
    end if
end do

!do i=1,6
!    if(myProc%procNeighbor(i) /= myProc%taskid) then
!        call MPI_WAIT(sendFinalRequest(i),sendFinalStatus,ierr)
!    end if
!end do

if(allocated(finalSend)) then
    deallocate(finalSend)
end if


end subroutine

!*****************************************************************************************
!> Subroutine update defect list multiple - updates defects according to reactions chosen
!!
!!This is the most involved and complex subroutine of the entire SRSCD algorithm.
!!This subroutine updates defects in the local mesh according to the reaction chosen, and
!!communicates with neighboring processors about defects that may have passed into 
!!a different processor as well as about defects that may have changed on the boundary.
!!It also creates a list of defects that have been updated, to inform the next subroutine
!!which reactions to update.
!!
!!This version of updateDefectList has been created for the case of one KMC domain per volume
!!element, in which many reactions will be chosen during each step. Thus, we have to go through
!!the list pointed to by reactionCurrent rather than just that one reaction.
!*****************************************************************************************

subroutine updateDefectListMultiple(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
use mod_constants
use DerivedType
use ReactionRates
use randdp
implicit none
include 'mpif.h'

type(reaction), pointer :: reactionCurrent
type(defect), pointer :: defectCurrent, defectPrev, defectTemp, defectStoreList, defectStore
type(cascadeDefect), pointer :: cascadeDefectTemp
type(defectUpdateTracker), pointer :: defectUpdate, defectUpdateCurrent
type(cascade), pointer :: CascadeCurrent, CascadePrev
type(cascadeEvent), pointer :: cascadeTemp
double precision coordinates(3), coordinatesTemp(3)
integer cellNumber, mixingEvents, mixingTemp, findCellWithCoordinatesFineMesh
logical combineBoolean, cascadeMixingCheck

integer i, j, k, l, m, same, products(numSpecies), tag, tracker, totalLocalRecv, chooseRandomCell
double precision diffusionRandom
logical flag

!Used for cascade recombination
double precision r1, atemp


!defect update counters
integer :: numUpdateLocal(6)=0			!<the number of defects being sent to each processor neighbor
integer :: numUpdateBndry(6)=0			!<the number of defects being recieved from each processor neighbor

integer count, numDefectsRecv
integer :: numUpdateLocalRecv(6)=0, numUpdateBndryRecv(6)=0	!number of defects being recieved from each proc neighbor

!NOTE: this final step could be eliminated by keeping the global mesh in each local processor
!(thus each element could be identified as being part of the local mesh of one proc and the boundary of 
!any other procs)
integer :: numUpdateFinal(6)=0, numUpdateFinalRecv(6)=0		!number of defects being sent/recieved in final update step

!create buffers of information to send to neighboring elements
integer, allocatable :: localBuffer(:,:,:)		!(direction, count, defectType+cellNumber+numDefects) array containing which processor to send to,
integer, allocatable :: bndryBuffer(:,:,:)		!the number of defects to send this way, and the defect information (defectType, cell number, and new number of defects)
integer, allocatable :: finalBuffer(:,:,:)		
integer, allocatable :: localBufferSend(:,:), localBufferRecv(:,:)	!temporary buffers used for sending/recieving data taken from localBuffer and bndryBuffer
integer, allocatable :: bndryBufferSend(:,:), bndryBufferRecv(:,:)
integer, allocatable :: finalBufferSend(:,:), finalBufferRecv(:,:)

integer status(MPI_STATUS_SIZE)

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
	
	subroutine chooseCascade(CascadeTemp)
		use DerivedType
		implicit none
		type(cascadeEvent), pointer :: CascadeTemp
	end subroutine
	
	subroutine initializeFineMesh(CascadeCurrent)
		use DerivedType
		implicit none
		type(cascade), pointer :: CascadeCurrent
	end subroutine
end interface


!create buffers: size greater than max size needed (max 3 reactants/products per reaction, numCells reactions)
allocate(localBuffer(numSpecies+3,numCells*3,6))
allocate(bndryBuffer(numSpecies+2,numCells*3,6))

!Initialize buffers here instead of inside the loop

do while(associated(reactionCurrent))	!loop through all non-null reactions chosen
	
	!***********************************************************************************************
	!
	!Cascade chosen
	!
	!***********************************************************************************************
	
	
	if(reactionCurrent%numReactants==-10) then 
		
		write(*,*) 'Error cascade implantation in one KMC domain per volume element mode'

	!***********************************************************************************************
	!
	!Defect update for reactions within the fine mesh.
	!
	!***********************************************************************************************
		
	else if(associated(CascadeCurrent)) then
		
		write(*,*) 'Error cascade implantation in one KMC domain per volume element mode'
	
	!***********************************************************************************************
	!
	! Defect update for reactions chosen in the coarse mesh.
	!
	! For defect changes in the boundary of other processors or in the boundary of this processor,
	! local and boundary buffers are created for communication with other processors.
	!
	!***********************************************************************************************

	else
		!Reactions in the coarse mesh

		!removing reactants from the system
		do i=1, reactionCurrent%numReactants
			
			!create a new element in defectUpdate and assign all variables except for num (will do later)
			allocate(defectUpdateCurrent%next)
			defectUpdateCurrent=>defectUpdateCurrent%next
			defectUpdateCurrent%proc=reactionCurrent%taskid(i)
			defectUpdateCurrent%dir=0	!not pointed at a different proc
			allocate(defectUpdateCurrent%defectType(numSpecies))
			do j=1,numSpecies
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(j,i)
			end do
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>defectList(reactionCurrent%cellNumber(i))
			
!			call findDefectInList(defectCurrent, defectPrev, reactionCurrent%reactants(:,i))
			
			!Check that defectCurrent is pointing towards the reactant
!			same=0
!			do j=1,numSpecies
!				if(defectCurrent%defectType(j)==reactionCurrent%reactants(j,i)) then
!					same=same+1
!				endif
!			end do
			
!			if(same /= numSpecies) then
!				write(*,*) 'Proc', myProc%taskid, 'error defectCurrent not pointing to reactants'
!				write(*,*) 'reactants', reactionCurrent%reactants(:,:)
!				write(*,*) 'products', reactionCurrent%products(:,:)
!				write(*,*) 'searching for', reactionCurrent%reactants(:,i)
!				if(associated(defectCurrent)) write(*,*) 'DefectCurrent', defectCurrent%defectType
!				if(associated(defectPrev)) write(*,*) 'DefectPrev', defectPrev%defectType
!			end if
			
			do while(associated(defectCurrent)) !search for reactant defect
				same=0
				do j=1,numSpecies
					if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
						same=same+1
					endif
				end do
				if(same==numSpecies) then
					exit
				endif
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
			end do
			
			!First update local buffer if needed
			do j=1, 6
				do k=1,myMesh(reactionCurrent%cellNumber(i))%numNeighbors(j)
					if(myMesh(reactionCurrent%cellNumber(i))%neighborProcs(k,j) /= myProc%taskid .AND. &
						myMesh(reactionCurrent%cellNumber(i))%neighborProcs(k,j) /= -1) then	!neighboring element not in this proc
						
						numUpdateLocal(j)=numUpdateLocal(j)+1
						
						do l=1,numSpecies
							localBuffer(l,numUpdateLocal(j),j)=reactionCurrent%reactants(l,i)	!defect type
						end do
						
						!Cell Number in local mesh
						localBuffer(numSpecies+1,numUpdateLocal(j),j)=reactionCurrent%cellNumber(i)
						
						!Cell number in boundary mesh
						localBuffer(numSpecies+3,numUpdateLocal(j),j)=myMesh(reactionCurrent%cellNumber(i))%neighbors(k,j)
						
						if(associated(defectCurrent)) then
							localBuffer(numSpecies+2,numUpdateLocal(j),j)=-1
						else
							write(*,*) 'Error tried to delete defect that was not there and send to neighboring proc'
						end if
					end if
				end do
			end do
			
			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there coarse'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'rate', reactionCurrent%reactionRate
				call MPI_ABORT(comm,ierr)
			
			!if there is one defect of this type and it is in the middle of the list, remove it from the list
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next) .AND. associated(defectPrev)) then
				!write(86,*) 'removing defectCurrent from cell', defectCurrent%defectType, 'cell', defectCurrent%cellNumber
				defectPrev%next=>defectCurrent%next !remove that defect type from the system
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)
				defectUpdatecurrent%num=-1			!tell updateReactionList that there are none of this defect left
			
			!if there is one defect of this type and it is at the end of the list then remove it from the list and
			!point defectPrev%next at nothing (nullify)
			else if(defectCurrent%num==1 .AND. associated(defectPrev)) then
				!write(86,*) 'removing defectCurrent from cell', defectCurrent%defectType, 'cell', defectCurrent%cellNumber
				deallocate(defectCurrent%defectType)
				deallocate(defectCurrent)
				nullify(defectCurrent)	!remove the last defect from the system
				nullify(defectPrev%next)
				defectUpdateCurrent%num=-1			!tell updateReactionList that there are none of this defect left
				
			!if there is one defect of this type and it is at the beginning of the list then just make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1 .AND. associated(defectCurrent%next)) then !removing first defect from cell i
				defectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
				defectUpdateCurrent%num=-1
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then 	!removing only defect from cell i (single helium) - this is redundant but will keep for now
				defectList(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) 'trying to remove a defect that isnt there'
				call MPI_ABORT(comm,ierr)
			else
				!decrease the number of defects by 1 if the number of defects is greater than 1
				defectCurrent%num=defectCurrent%num-1 !remove the defect from the system instead of the entire entry in the list
				defectUpdateCurrent%num=-1	!tell updateReactionList the new number of defects in the cell
			endif
		end do
	
		!***********************************************************************
		!Here, I will assume cubic cells. If a defect migrates from one cell to another,
		!the defect will be given a percent chance of removal from the system
		!based on the cell size (%chance of removal = cell size/grain size)
		!This is to replicate the OKMC practice of removing SIA clusters after they
		!migrate 1 um.
		!***********************************************************************
		flag=.FALSE.
		if(grainBoundaryToggle=='yes') then	
		
		!we have toggled the use of grain boundaries to remove defects from system
		!
		!7/18/2014: NOTE that if we have diffusion from coarse mesh to fine mesh, the length
		!of the coarse mesh element will be used to calculate the probability of grain boundary
		!absorption, instead of the effective coarse to fine length. This could be changed if desired.
		
			if(reactionCurrent%numReactants==1 .and. reactionCurrent%numProducts==1) then 
				diffusionRandom=dprand()
				
				!randomly choose whether to remove this defect from the system according to the mean free path (parameter) and the
				!length of the volume element that the defect is currently in
				if(diffusionRandom <= myMesh(reactionCurrent%cellNumber(1))%length/meanFreePath) then
					flag=.TRUE.
					!reactionCurrent%numProducts=0 !remove that defect from system
					!we can no longer remove the defect in this way because changing the reaction list
					!is no longer allowed (reaction list is not re-created at each step)
				endif		
			endif
		endif
	
		!adding products to the system
		if(flag .eqv. .FALSE.) then		!if the diffusion reaction defect did not get removed from the system
			
			!***********************************************************************************************
			!EDIT: 5/27/2014: Included the possibility of updating defects in neighboring processors
			!(the information should be stored in reaction%taskid(:))
			!
			!Procedure: if taskid .NE. myProc%taskid, then add to a buffer the number of the proc to send to,
			!the information on what defects are being added and in what volume elements
			!
			!Then send a flag with the number of defects being sent/recieve a flag with number of defects being
			!recieved
			!
			!Then send/recieve information on defects being added to system and add them to the correct volume
			!elements.
			!
			!RE-EDIT: 5/28/2014: Included the boundary mesh: this includes the updating of defects in neighboring
			!processors as well.
			!
			!EDIT: 7/18/2014: Added diffusion from coarse to fine mesh: if reactionCurrent%cellNumber .LT. 0
			!***********************************************************************************************
			
			do i=1, reactionCurrent%numProducts
				
				if(reactionCurrent%taskid(i+reactionCurrent%numReactants) == -1) then
					!we have identified a defect that is going to be removed from the system via a free surface
					!Do nothing; no need to add this defect to any lists or update any reaction lists because of it
					
				else if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants) < 0) then
					!we have identified a defect that is diffusing from the coarse mesh to the fine mesh
					!and therefore neet to treat it differently than diffusion within the coarse mesh.
					!
					!In this case, we do not update the buffers because coarse - to - fine reactions
					!are within the same processor and fine meshes are not within the boundary of 
					!neighboring processors.
					
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
					allocate(defectUpdateCurrent%defectType(numSpecies))
					
					!When coarse-to-fine reactions are chosen, a random cell within the fine mesh is chosen.
					defectUpdateCurrent%cellNumber=chooseRandomCell()	!cellNumber > 0
					
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					
					!In coarse-to-fine reactions, the cascade number is stored in reactionCurrent%cellNumber
					!and identified with a negative. Therefore we make it positive again.
					defectUpdateCurrent%cascadeNumber=-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					
					nullify(defectPrev)
					do j=1,numSpecies
						products(j)=reactionCurrent%products(j,i)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
					end do
					
					!Point CascadeCurrent at the correct cascade
					
					CascadeCurrent=>ActiveCascades
					
					do while(associated(CascadeCurrent))
						
						if(CascadeCurrent%cascadeID==-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)) then
							exit
						end if
						
						CascadeCurrent=>CascadeCurrent%next
					
					end do
					
					!If we went through the previous loop without pointing at a cascade, then we have
					!a problem with placing the coarse-to-fine miration in the correct cascade.
					
					if(.NOT. associated(CascadeCurrent)) then
						write(*,*) 'error coarse-to-fine cascade not correctly pointed'
					end if
					
					!Point defectCurrent at the defect list within the correct element in CascadeCurrent
					defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
					
					!Next update defects
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						end do
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+1
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectUpdateCurrent%cellNumber
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(j,i)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectUpdateCurrent%cellNumber
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(j,i)
						end do
					end if
					
				else if(reactionCurrent%taskid(i+reactionCurrent%numReactants) /= myProc%taskid .AND. &
					reactionCurrent%taskid(i+reactionCurrent%numReactants) /= -1) then	!coarse-to-coarse diffusion, different processors
					
					!we have identified a defect that needs to be added to a neighboring processor and the boundary mesh
					
					!***********************************************************************************
					!EDIT: 6/4/2014: we only want to add to the boundary buffer once in this step
					!because this will be used to change only one local defect. However, if the number
					!of processors is large, then one mesh can have a processor bounding it on both sides
					!(assuming periodic boundary conditions). Therefore, we have adjusted the code below
					!to only add to bndryBuffer once even though we change myBoundary more than once (if
					!the same cell bounds this mesh in two places due to PBCs)
					!
					!This is accomplished using the logical flag. Once the defect is added to the bounary
					!buffer in one direction, it cannot occur again and therefore we won't double-count
					!defect migration.
					!***********************************************************************************
					
					!step 1: find the element in the boundary mesh. The element may bound this processor at more than one place.
					!**************************************************************************************************
					!EDIT: 6/4/2014: this may not be the fastest way to search for which cells in myBoundary to update.
					!EDIT: 6/5/2014: Updated the way that we search through the local mesh for boundary elements which
					!correspond to the one in our reaction. Added %localNeighbor to the boundary class so that we know
					!exactly what the local element is whose neighbor is in the boundary. Now, we know the boundary element
					!number so we just search in all 6 directions for a boundary element that is allocated and that 
					!has a local neighbor in the local mesh. (Avoids searching through all local elements for one
					!with a boundary element corresponding to the boundary element we are looking for)
					!NOTE: old code commented out below
					!NOTE: this both increased the speed of computation and reduced the number of steps. This is surprising
					!as we are not changing anything physical within the system, but I have not found a bug and the results
					!match expected results so far. Will continue to keep an eye out for a bug in this routine.
					!**************************************************************************************************
					flag=.FALSE.
					
					do j=1,6
						do k=1,myMesh(reactionCurrent%cellNumber(1))%numNeighbors(j)
							if(flag .EQV. .FALSE.) then
								if(myMesh(reactionCurrent%cellNumber(1))%neighbors(k,j)==&
									reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
									myMesh(reactionCurrent%cellNumber(1))%neighborProcs(k,j)==&
									reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
			
									numUpdateBndry(j)=numUpdateBndry(j)+1
									do l=1,numSpecies
										bndryBuffer(l,numUpdateBndry(j),j)=reactionCurrent%products(l,i)
									end do
									bndryBuffer(numSpecies+1,numUpdateBndry(j),j)=&
										reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
									bndryBuffer(numSpecies+2,numUpdateBndry(j),j)=1	!indicates we are adding one defect
									flag=.TRUE.
								endif
							endif
						end do
					end do
					
					!update the direction in defectUpdateCurrent (used for updating reaction lists)
					!create a new element in defectUpdate and assign all variables except for num and dir (will do later)
					do j=1,6
	
						!if myBoundary at this cell has defectList associated, then it is a boundary element for a local element in this cell.
						!NOTE: this will not work in the case of a non-uniform mesh (because there could be more
						!than one local neighbor, will have to change to localNeighbor())
						if(myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%localNeighbor /= 0 &
							.AND. myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%proc == &
							reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
							
							!Create next entry in defectUpdate (used for updating reaction lists)
							allocate(defectUpdateCurrent%next)
							defectUpdateCurrent=>defectUpdateCurrent%next
							defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
							allocate(defectUpdateCurrent%defectType(numSpecies))
							defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							nullify(defectUpdateCurrent%next)
							defectUpdateCurrent%dir=j
							defectUpdateCurrent%num=1	!used for updating reaction lists
							defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
							
							!This information is used to tell us which element in the local
							!processor bounds this element in myBoundary
							defectUpdateCurrent%neighbor=myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)&
								%localNeighbor		
							
							do l=1,numSpecies
								products(l)=reactionCurrent%products(l,i)
								defectUpdateCurrent%defectType(l)=reactionCurrent%products(l,i)
							end do
							
							!Find defect in defect list on myBoundary
							defectCurrent=>myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),j)%defectList
	
							call findDefectInList(defectCurrent, defectPrev, products)
	
							!Insert defect in myBoundary(dir,cell)%defectList
							if(associated(defectCurrent)) then !if we aren't at the end of the list
								same=0
								do l=1,numSpecies
									if(defectCurrent%defectType(l)==products(l)) then
										same=same+1
									endif
								end do
								if(same==numSpecies) then	!if the defect is already present in the list
									defectCurrent%num=defectCurrent%num+1
								else		!if the defect is to be inserted in the list
									nullify(defectPrev%next)
									allocate(defectPrev%next)
									nullify(defectPrev%next%next)
									defectPrev=>defectPrev%next
									allocate(defectPrev%defectType(numSpecies))
									defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
									defectPrev%num=1
									do l=1,numSpecies
										defectPrev%defectType(l)=reactionCurrent%products(l,i)
									end do
									defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
								endif
							else 			!add a defect to the end of the list
								nullify(defectPrev%next)
								allocate(defectPrev%next)
								nullify(defectPrev%next%next)
								defectPrev=>defectPrev%next
								allocate(defectPrev%defectType(numSpecies))
								defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
								defectPrev%num=1
								do l=1,numSpecies
									defectPrev%defectType(l)=reactionCurrent%products(l,i)
								end do
							endif
						endif
					end do
					
				else	!coarse-to-coarse, same processor
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
					allocate(defectUpdateCurrent%defectType(numSpecies))
					defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
					
					nullify(defectPrev)
					do j=1,numSpecies
						products(j)=reactionCurrent%products(j,i)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
					end do
					
					defectCurrent=>defectList(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
	
					!First update local buffer if needed
					do j=1,6
						do k=1,myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%numNeighbors(j)
							if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(k,j)==-1) then
								!do nothing, free surface
							else if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(k,j) &
								/= myProc%taskid) then	!neighboring element not in this proc
							
								numUpdateLocal(j)=numUpdateLocal(j)+1
								
								!species
								do l=1,numSpecies
									localBuffer(l,numUpdateLocal(j),j)=reactionCurrent%products(l,i)
								end do
								
								!cell number in local mesh
								localBuffer(numSpecies+1,numUpdateLocal(j),j)=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
								
								!cell number of neighbor (in different proc)
								localBuffer(numSpecies+3,numUpdateLocal(j),j)=myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))&
									%neighbors(k,j)
								
								!number of defects to be increased by 1
								localBuffer(numSpecies+2,numUpdateLocal(j),j)=1

							end if
						end do
					end do
					
					!Next update defects
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							end if
						end do
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+1
						else		!if the defect is to be inserted in the list
							
							if(.NOT. associated(defectPrev)) then
								write(*,*) 'Error defectPrev not associated'
								write(*,*) defectCurrent%defectType
								write(*,*) products
							endif
							
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							defectPrev%num=1
							do j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(j,i)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
						defectPrev%num=1
						do j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(j,i)
						end do
					end if
					
				end if
				
			end do
		
		else
			!If the defect has been removed due to grain boundary absorption
		end if
	end if
	
	reactionCurrent=>reactionCurrent%next

end do	!do 1000 while(associated(reactionCurrent))

!*********
!Print defects after this step
!*********

!if(myProc%taskid==MASTER) then
!	do 212 i=1,numCells
!		defectCurrent=>defectList(i)
!		do 213 while(associated(defectCurrent))
!			write(*,*)'processor', myProc%taskid,defectCurrent%defectType, defectCurrent%cellNumber, defectCurrent%num
!			defectCurrent=>defectCurrent%next
!		213 continue
!	212 continue
!endif
	
!*************
!Step 2: send/recieve data about local and boundary defects that have changed
!*************

do i=1,6

	if(myProc%procNeighbor(i) /= myProc%taskid) then

		call MPI_SEND(numUpdateLocal(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+i,comm, ierr)	!number of local defects that have changed on boundary of processor i
		call MPI_SEND(numUpdateBndry(i),1,MPI_INTEGER,myProc%procNeighbor(i),i+6,comm, ierr)		!number of defects in the mesh of processor i that have changed (diffusion only)
		
		!EDIT: need to only send the first numUpdateLocal(i) elements of localBuffer(i,:,:)
		
		if(numUpdateLocal(i) /= 0) then

			allocate(localBufferSend(numSpecies+3,numUpdateLocal(i)))
			
			do j=1,numUpdateLocal(i)
				do k=1,numSpecies+3
					localBufferSend(k,j)=localBuffer(k,j,i)
				end do
			end do

			call MPI_SEND(localBufferSend,(numSpecies+3)*numUpdateLocal(i), MPI_INTEGER, &
				myProc%procNeighbor(i),i+12,comm,ierr)
!			if(myProc%taskid==4) then
!				write(*,*) 'dir', i, 'proc sent', myProc%procNeighbor(i), 'numUpdateLocal', numUpdateLocal(i)
!			endif

			deallocate(localBufferSend)

		end if
		
		if(numUpdateBndry(i) /= 0) then

			allocate(bndryBufferSend(numSpecies+2,numUpdateBndry(i)))
			
			do j=1,numUpdateBndry(i)
				do k=1,numSpecies+2
					bndryBufferSend(k,j)=bndryBuffer(k,j,i)
				end do
			end do

			call MPI_SEND(bndryBufferSend,(numSpecies+2)*numUpdateBndry(i), MPI_INTEGER, &
				myProc%procNeighbor(i),i+18,comm,ierr)
			deallocate(bndryBufferSend)

		end if

	end if
end do

!**************
!Step 3: Recieve data about local/bdry defects that have changed and update defectList and myBoundary accordingly
!**************

do i=1,6
	numUpdateFinal(i)=0
end do
totalLocalRecv=0

do i=1,6

	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) /= myProc%taskid) then

		call MPI_RECV(numUpdateBndryRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+tag,comm,status,ierr)	!number of bndry defects that have changed (local to processor i, boundary here)
		call MPI_RECV(numUpdateLocalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+6,comm,status,ierr)	!number of local defects that have changed (bndry to processor i, local here)
		totalLocalRecv=totalLocalRecv+numUpdateLocalRecv(i)
	endif
end do

allocate(finalBuffer(numSpecies+3,totalLocalRecv,6))

do i=1,6
	
	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) /= myProc%taskid) then
		
		if(numUpdateBndryRecv(i) /= 0) then
			
			!Read in defects to update in boundary
			allocate(bndryBufferRecv(numSpecies+3,numUpdateBndryRecv(i)))
			call MPI_RECV(bndryBufferRecv,(numSpecies+3)*numUpdateBndryRecv(i),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+12,comm,status,ierr)

			!Add defects in bndryBufferRecv to myBoundary()
			do j=1,numUpdateBndryRecv(i)

				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%procNeighbor(i)
				defectUpdateCurrent%dir=i	
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=bndryBufferRecv(numSpecies+1,j)
				defectUpdateCurrent%num=bndryBufferRecv(numSpecies+2,j)	!This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)
				
				!This information is used to tell us which element in the local
				!processor bounds this element in myBoundary
				defectUpdateCurrent%neighbor=bndryBufferRecv(numSpecies+3,j)
				
				!point defectCurrent at the defect list in the correct cell of myBoundary
				defectCurrent=>myBoundary(bndryBufferRecv(numSpecies+1,j),i)%defectList
				if(.NOT. associated(defectCurrent)) then
					write(*,*) 'error myBoundary not allocated correctly'
					write(*,*) 'dir', i, 'cell', bndryBufferRecv(numSpecies+1,j)
					call MPI_ABORT(comm,ierr)
				endif
				
				nullify(defectPrev)
				
				do k=1,numSpecies
					products(k)=bndryBufferRecv(k,j)
					defectUpdateCurrent%defectType(k)=bndryBufferRecv(k,j)	!used to update reaction lists
				end do
				
				!point DefectCurrent at the defect we are looking for (if it is there), otherwise
				!point defectCurrent after and defectPrev before where it should go in defectList.
!				if(myProc%taskid==MASTER .OR. myProc%taskid==3) then
!					write(*,*) 'proc', myProc%taskid
!					write(*,*) 'numUpdateBndryRecv', numUpdateBndryRecv(i)
!					write(*,*) 'dir', i, (bndryBufferRecv(k,j),k=1,numSpecies+2)
!				endif
				
				call findDefectInList(defectCurrent, defectPrev, products)

				!Next update defects in myBoundary
				if(associated(defectCurrent)) then !if we aren't at the end of the list
					
					same=0
					do l=1,numSpecies
						if(defectCurrent%defectType(l)==products(l)) then
							same=same+1
						end if
					end do
					if(same==numSpecies) then	!if the defect is already present in the list
						defectCurrent%num=defectCurrent%num+bndryBufferRecv(numSpecies+2,j) !this will be +/- 1 only
						
						if(defectCurrent%num==0) then
							!delete this defect from the list in myBoundary
							defectPrev%next=>defectCurrent%next !remove that defect type from the system
							deallocate(defectCurrent%defectType)
							deallocate(defectCurrent)
							nullify(defectCurrent)
						end if
					else		!if the defect is to be inserted in the list
					
						if(bndryBufferRecv(numSpecies+2,j) == -1) then
							write(*,*) 'defectCurrent associated', defectCurrent%defectType, 'num', defectCurrent%num
							if(associated(defectPrev)) write(*,*) 'defectPrev', defectPrev%defectType, 'num', defectPrev%num
							write(*,*) 'error in defectUpdate negative defect numbers', myProc%taskid
							write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
							write(*,*) (bndryBufferRecv(k,j),k=1,numSpecies+2)
							
							call MPI_ABORT(comm,ierr)
							
						else
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=bndryBufferRecv(numSpecies+1,j)
							defectPrev%num=bndryBufferRecv(numSpecies+2,j)	!This will be +/- 1 only
							
							do l=1,numSpecies
								defectPrev%defectType(l)=bndryBufferRecv(l,j)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
												
					end if
				else 			!add a defect to the end of the list
					
					if(bndryBufferRecv(numSpecies+2,j) == -1) then
						write(*,*) 'defectCurrent not associated'
						write(*,*) 'error in defectUpdate negative defect numbers'
						write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
						write(*,*) (bndryBufferRecv(k,j),k=1,numSpecies+2)
						call MPI_ABORT(comm,ierr)
					else
					
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=bndryBufferRecv(numSpecies+1,j)
						defectPrev%num=bndryBufferRecv(numSpecies+2,j)	!This will be +/- 1 only
											
						do l=1,numSpecies
							defectPrev%defectType(l)=bndryBufferRecv(l,j)
						end do
					end if

				end if
			end do
			deallocate(bndryBufferRecv)
		end if
		
		if(numUpdateLocalRecv(i) /= 0) then
			
			!Read in defects to update in local mesh
			allocate(localBufferRecv(numSpecies+2,numUpdateLocalRecv(i)))
			call MPI_RECV(localBufferRecv,(numSpecies+2)*numUpdateLocalRecv(i),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+18,comm,status,ierr)

			!Add defects in localBufferRecv to defectList()
			do j=1,numUpdateLocalRecv(i)
				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%taskid
				defectUpdateCurrent%dir=0	!not pointed at a different proc
				defectUpdateCurrent%neighbor=0	!not pointed at a different proc
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=localBufferRecv(numSpecies+1,j)
				defectUpdateCurrent%num=localBufferRecv(numSpecies+2,j) !This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)

				!point defectCurrent at the defect list in the correct cell of defectList
				defectCurrent=>defectList(localBufferRecv(numSpecies+1,j))
				nullify(defectPrev)
				
				do k=1,numSpecies
					products(k)=localBufferRecv(k,j)
					defectUpdateCurrent%defectType(k)=localBufferRecv(k,j)
				end do
				
				!point DefectCurrent at the defect we are looking for (if it is there), otherwise
				!point defectCurrent after and defectPrev before where it should go in defectList.

				call findDefectInList(defectCurrent, defectPrev, products)

				!Next update defects in local mesh
				if(associated(defectCurrent)) then !if we aren't at the end of the list
					same=0
					do l=1,numSpecies
						if(defectCurrent%defectType(l)==products(l)) then
							same=same+1
						endif
					end do
					if(same==numSpecies) then	!if the defect is already present in the list
						defectCurrent%num=defectCurrent%num+localBufferRecv(numSpecies+2,j)	!This will be +/- 1 only (should be only +1 here)
						
						if(defectCurrent%num==0) then
							write(*,*) 'Error zero defects in updateDefectList, step 3 of communication'
						endif
					else		!if the defect is to be inserted in the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=localBufferRecv(numSpecies+1,j)
						defectPrev%num=localBufferRecv(numSpecies+2,j) !This will be +/- 1 only (should be only +1 here)
						do l=1,numSpecies
							defectPrev%defectType(l)=localBufferRecv(l,j)
						end do
						defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
					endif
				else 			!add a defect to the end of the list
					nullify(defectPrev%next)
					allocate(defectPrev%next)
					nullify(defectPrev%next%next)
					defectPrev=>defectPrev%next
					allocate(defectPrev%defectType(numSpecies))
					defectPrev%cellNumber=localBufferRecv(numSpecies+1,j)
					defectPrev%num=localBufferRecv(numSpecies+2,j) !This will be +/- 1 only (should be only +1 here)
					do l=1,numSpecies
						defectPrev%defectType(l)=localBufferRecv(l,j)
					end do
				endif
				
				!*****************************************************
				!Prep for step 4: check to see if element updated has more than one neighbor proc
				!*****************************************************
				
				do k=1,6
					do l=1,myMesh(localBufferRecv(numSpecies+1,j))%numNeighbors(k)
						
						!if a neighbor of this element is not in this proc or in the proc that just communicated with it
						if (myMesh(localBufferRecv(numSpecies+1,j))%neighborProcs(l,k) /= myProc%taskid .AND. &
							myMesh(localBufferRecv(numSpecies+1,j))%neighborProcs(l,k) /= myProc%procNeighbor(i) .AND. &
							myMesh(localBufferRecv(numSpecies+1,j))%neighborProcs(l,k) /= -1) then
							
							!Add this defect to a final buffer
							numUpdateFinal(k)=numUpdateFinal(k)+1
							!finalBuffer(1,1,1)=100	!why is this here?
							do m=1,numSpecies
								finalBuffer(m,numUpdateFinal(k),k)=localBufferRecv(m,j)
							end do
							finalBuffer(numSpecies+1,numUpdateFinal(k),k)=localBufferRecv(numSpecies+1,j)	!local cell number
							finalBuffer(numSpecies+2,numUpdateFinal(k),k)=localBufferRecv(numSpecies+2,j) !This should be +1 only
							finalBuffer(numSpecies+3,numUpdateFinal(k),k)=myMesh(localBufferRecv(numSpecies+1,j))%neighbors(l,k)	!cell number in neighboring proc
							
							if(localBufferRecv(numSpecies+2,j)==-1) then
								write(*,*) 'error: cell in boundary of multiple procs removing defect'
								call MPI_ABORT(comm,ierr)
							endif
!							if(myProc%taskid==4) then
!								if(k==5) then
!									write(*,*) 'dir', k, 'final buffer',(finalBuffer(k,numUpdateFinal(k),m),m=1,numSpecies+2), 'add'
!								else if(k==6) then
!									write(*,*) 'dir', k, 'final buffer',(finalBuffer(k,numUpdateFinal(k),m),m=1,numSpecies+2), 'add'
!								endif
!							endif
						end if

					end do
				end do
			end do

			deallocate(localBufferRecv)
		end if
		
	end if
end do

!*************
!Step 4: if a local defect is updated due to diffusion to the boundary of another processor,
!	and the local element is also in the boundary of a processor other than the one already 
!	communicated with, let all other processors know of local defects that have been updated
!*************

do i=1,6
	if(myProc%procNeighbor(i) /= myProc%taskid) then
		call MPI_SEND(numUpdateFinal(i), 1, MPI_INTEGER, myProc%procNeighbor(i), i+24, comm, ierr)	!number of local defects that have changed on boundary of processor i
		
		if(numUpdateFinal(i) /= 0) then
			allocate(finalBufferSend(numSpecies+3,numUpdateFinal(i)))
			
			do j=1,numUpdateFinal(i)
				do k=1,numSpecies+3
					finalBufferSend(k,j)=finalBuffer(k,j,i)
				end do
			end do
			
			call MPI_SEND(finalBufferSend,(numSpecies+3)*numUpdateFinal(i), MPI_INTEGER, &
				myProc%procNeighbor(i),i+30,comm,ierr)
			deallocate(finalBufferSend)
		end if
	end if
end do

!*************
!Step 5: update myBoundary again, based on information sent in step 4
!*************

do i=1,6

	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) /= myProc%taskid) then
		
		
		!number of bndry defects that have changed (local to processor i, boundary here)
		call MPI_RECV(numUpdateFinalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+24,comm,status,ierr)
		
		if(numUpdateFinalRecv(i) /= 0) then
			
			!Read in defects to update in boundary
			allocate(finalBufferRecv(numSpecies+3,numUpdateFinalRecv(i)))
			call MPI_RECV(finalBufferRecv,(numSpecies+3)*numUpdateFinalRecv(i),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+30,comm,status,ierr)

			!Add defects in finalBufferRecv to myBoundary()
			do j=1,numUpdateFinalRecv(i)

				if(finalBufferRecv(numSpecies+2,j)==-1) then
					write(*,*) 'error deleting defects in finalBufferRecv'
					!call MPI_ABORT(comm,ierr)
				else
					
					!create a new element in defectUpdate and assign all variables except for num (will do later)
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=myProc%procNeighbor(i)
					defectUpdateCurrent%dir=i	
					allocate(defectUpdateCurrent%defectType(numSpecies))
					defectUpdateCurrent%cellNumber=finalBufferRecv(numSpecies+1,j)	!cell number in boundary mesh
					defectUpdateCurrent%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
					defectUpdateCurrent%neighbor=finalBufferRecv(numSpecies+3,j) !cell number in myMesh
					defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
					nullify(defectUpdateCurrent%next)
					
					!point defectCurrent at the defect list in the correct cell of myBoundary
					defectCurrent=>myBoundary(finalBufferRecv(numSpecies+1,j),i)%defectList
					nullify(defectPrev)
					
					do k=1,numSpecies
						products(k)=finalBufferRecv(k,j)
						defectUpdateCurrent%defectType(k)=finalBufferRecv(k,j)
					end do
					
					!point DefectCurrent at the defect we are looking for (if it is there), otherwise
					!point defectCurrent after and defectPrev before where it should go in defectList.
	
					call findDefectInList(defectCurrent, defectPrev, products)
	
					!Next update defects in myBoundary
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do l=1,numSpecies
							if(defectCurrent%defectType(l)==products(l)) then
								same=same+1
							end if
						end do
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+finalBufferRecv(numSpecies+2,j) !This should be +1 only
							
							if(defectCurrent%num==0) then
								write(*,*) 'Error zero defects in updateDefectList, step 5 of communication'
							end if
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=finalBufferRecv(numSpecies+1,j)
							defectPrev%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
							do l=1,numSpecies
								defectPrev%defectType(l)=finalBufferRecv(l,j)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=finalBufferRecv(numSpecies+1,j)
						defectPrev%num=finalBufferRecv(numSpecies+2,j) !This should be +1 only
						do l=1,numSpecies
							defectPrev%defectType(l)=finalBufferRecv(l,j)
						end do
					end if
				
				end if
			end do
			
			deallocate(finalBufferRecv)
		end if
	end if

end do

!NOTE: steps 4-5 are necessary in the case that one element is a member of the boundary of 
!more than one processor (ex: the corner of a processor's mesh).

end subroutine

!***************************************************************************************************
!>Subroutine update reaction list - updates reaction rates and creates/deletes reactions according
!!to defect numbers that have changed:
!!
!!This subroutine does the following:
!!
!!1) Find list of affected volume elements (and their procs) - separate into list of neighbors and list of elements involved in reaction
!!
!!1a) Make list of all defect types in each volume element that need to be checked
!!
!!2)	Delete all diffusion reactions associated with reactants and products in neighboring volume elements
!!
!!3) Delete all reactions associated with reactants and products in volume elements associated with reaction
!!
!!4) Add all diffusion reactions associated with reactants and products in neighboring volume elements
!!
!!5) Add all reactions associated with reactants and prodcuts in volume elements associated with reaction
!!
!!In steps 2-5, subtract/add to totalRate each time we add or subtract a reaction from the list
!!
!!EDIT: it may be faster computationally to not delete and re-add reactions, but to just update reaction rates
!!for all appropriate reactions, delete if 0, and add if non-existent.
!!
!!In this way, we are ONLY updating the relevant reaction rates, not all possible reaction rates in the system.
!!
!!As reaction rates are updated, the total reaction rate for the processor is updated as well.
!***************************************************************************************************

!***************************************************************************************************
!06/04/2014: The bug with two defects moving to the same cell at the same time (see discussion before
!updateDefectList) does not affect updateReactionList because defectUpdateCurrent%num is not used 
!to update reaction lists (instead, the actual defect numbers are taken from myMesh(i)%defectList)
!
!Thus, no updates were made to this subroutine for the bug fix.
!***************************************************************************************************

subroutine updateReactionList(defectUpdate)
use DerivedType
use mod_constants
use ReactionRates
implicit none
include 'mpif.h'

type(defectUpdateTracker), pointer :: defectUpdate, defectUpdateCurrent, defectUpdatePrev, defectUpdateNext
type(defect), pointer :: defectCurrent
type(cascade), pointer :: cascadeCurrent
integer numElements, numNeighbors, i, j, k, dir, count, defectTemp(numSpecies), tracker
integer findNumDefectBoundary, findNumDefect
logical flag
integer localGrainID, neighborGrainID

!temporary
double precision temp

defectUpdateCurrent=>defectUpdate%next	!the first is zero

!loop through defects that needs reactions updated in defectUpdate
do while(associated(defectUpdateCurrent))
	
	if(defectUpdateCurrent%num /= 1 .AND. defectUpdateCurrent%num /= -1) then
		!we have an error; all defectUpdateCurrent members should have num = +/-1 (indicates adding or removing)
		write(*,*) 'error defectUpdateCurrent%num not equal to +/- 1', myProc%taskid, defectUpdateCurrent%num
	end if
	
	do i=1,numSpecies
		defectTemp(i)=defectUpdateCurrent%defectType(i)
	end do
			
	!if the defect is within the local mesh, update all relevant reactions
	
	if(defectUpdateCurrent%proc==myProc%taskid) then
		
		!*******************************************************************************************
		! defectUpdateCurrent%cascadeNumber==0 means a defect has changed in the coarse mesh
		! if the defect is within the coarse mesh, then update all relevant reactions in the coarse mesh
		!*******************************************************************************************

		if(defectUpdateCurrent%cascadeNumber==0) then
			
			!Single-defect reactions associated with defects of type defectTemp	
			call addSingleDefectReactions(defectUpdateCurrent%cellNumber,defectTemp)
			
			!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
			!all defect types in the defect list)
			defectCurrent=>defectList(defectUpdateCurrent%cellNumber)
			
			do while(associated(defectCurrent))
				if(defectCurrent%num /= 0) then
					call addMultiDefectReactions(defectUpdateCurrent%cellNumber,defectTemp, defectCurrent%defectType)
				end if
				defectCurrent=>defectCurrent%next
			end do
			
			!NOTE: there is a problem with annihilation reactions: if two defects annihilate and the
			!result is that there are no defects left of those two types in the cell (example, SIA+V->0
			!such that there is nothing left in the cell) after the reaction, addMultiDefectReactions will
			!not remove this reaction from the reaction list because defectCurrent will never point to 
			!one of the two defects. Therefore we must call addMultiDefectReactions using the two defect types
			!listed in defectUpdate
			
			defectUpdateNext=>defectUpdateCurrent%next
			if(associated(defectUpdateNext)) then
				if(defectUpdateNext%cellNumber==defectUpdateCurrent%cellNumber) then
					!we have found two defects that need their reactions updated in the same cell. Thus
					!we update any clustering reactions between them. This may be redundant but prevents
					!us from missing total-annihilation reactions as described above.
					call addMultiDefectReactions(defectUpdateCurrent%cellNumber, defectTemp, defectUpdateNext%defectType)
				end if
			end if
			
			!*******************
			!Diffusion reactions
			!
			!NOTE: diffusion reactions need to be updated in both directions. That is, we need to 
			!update the diffusion reactions in this cell as well as in the neighboring cell
			!(thus taking care of both defect addition and defect removal).
			!
			!Exception: if the neighboring cell is in a different processor, we don't need to 
			!update the neighbor-to-cell diffusion rate because it will be taken care of in the other
			!processor
			!*******************
			do j=1,6
				if (myMesh(defectUpdateCurrent%cellNumber)%numNeighbors(j)==0) then
					write(*,*) 'error myMesh does not have neighbors in this direction'
				end if

				!Add diffusion reactions from this cell to neighboring cells
				if(polycrystal=='yes') then

					!Find the grain ID number of the volume element we are in
					localGrainID=myMesh(defectUpdateCurrent%cellNumber)%material

					!Find the grain ID number of the neighboring volume element
					if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j) /= myProc%taskid .AND. &
							myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j) /= -1) then

						neighborGrainID=myBoundary(myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j),j)%material
					else if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j) == -1) then
						neighborGrainID=localGrainID	!free surface release, don't need to do anything special
					else
						neighborGrainID=myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j))%material
					endif

					if(localGrainID==neighborGrainID) then
						!Allow diffusion between elements in the same grain
						call addDiffusionReactions(defectUpdateCurrent%cellNumber, myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j),&
								myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j),j,defectTemp)
					else
						!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
						call addDiffusionReactions(defectUpdateCurrent%cellNumber, 0, myProc%taskid, -1, j, defectTemp)
					end if
				else
					call addDiffusionReactions(defectUpdateCurrent%cellNumber, myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j),&
							myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j),j,defectTemp)
				end if

				!If the neighboring volume element is in the same processor as this element, add
				!diffusion reactions from neighboring cells into this cell (changes reaction list
				!in neighboring cell)
				!
				!NOTE: don't need to switch j (direction) here even though diffusion is in opposite
				!direction because dir is only used when diffusing into/out of neighbor processor
				if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(1,j)==myProc%taskid) then
					!Don't need to do this for diffusion between different grains
					if(polycrystal=='yes' .AND. myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j))%material &
							== myMesh(defectUpdateCurrent%cellNumber)%material) then
						call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j), defectUpdateCurrent%cellNumber, &
								myProc%taskid, myProc%taskid, j, defectTemp)
					else if(polycrystal=='no') then	!2015.05.20 Noticed that this wasn't here, probably was a bug. Fixed?
						call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(1,j), defectUpdateCurrent%cellNumber, &
								myProc%taskid, myProc%taskid, j, defectTemp)
					end if
				end if
			end do
			
			!***********************************************************************************
			!Diffusion between coarse mesh and fine mesh - coarse to fine
			!
			!For each cascade in the coarse mesh element,  find the number of defects of same type
			!in the fine mesh. Then find the diffusion rate from coarse mesh into fine mesh
			!based on concentration of defects in coarse mesh and fine mesh. Not diffusing
			!into individual cells of fine mesh, as there are too many possibilities (comp. 
			!impractical). Cell will be chosen at random.
			!
			!Separate subroutines created to calculate the number of defects in cascade and to
			!add diffusion reaction
			!
			!To indicate that the diffusion reaction is going from the coarse to the fine mesh,
			!use a negative cell number with the cascadeID in reaction%cellNumber.
			!***********************************************************************************
				
			CascadeCurrent=>ActiveCascades
			
			do while(associated(CascadeCurrent))
			
				if(CascadeCurrent%cellNumber==defectUpdateCurrent%cellNumber) then

					call addDiffusionCoarseToFine(defectUpdateCurrent%cellNumber, myProc%taskid, CascadeCurrent, defectTemp)
				
				endif
				
				CascadeCurrent=>CascadeCurrent%next
				
			end do

		!*******************************************************************************************
		! DefectUpdateCurrent%cascadeNumber .NE. 0 means that the defect that has been updated is inside a cascade.
		! 
		! Therefore we update all relevant reactions within the correct cascade's fine mes
		!*******************************************************************************************
		else
			!Point CascadeCurrent at the cascade associated with the chosen reaction
			CascadeCurrent=>ActiveCascades
			
			do while(associated(CascadeCurrent))
				
				if(CascadeCurrent%cascadeID==defectUpdateCurrent%cascadeNumber) then
					exit
				endif
			
				CascadeCurrent=>CascadeCurrent%next
			end do

			!Check to make sure that we have exited at the correct cascade
			if(.NOT. associated(CascadeCurrent)) then
				write(*,*) 'error CascadeCurrent not associated in updateReactionList'
			end if

			!Single-defect reactions associated with defects of type defectTemp in fine mesh	
			call addSingleDefectReactionsFine(defectUpdateCurrent%cascadeNumber,defectUpdateCurrent%cellNumber,defectTemp)

			!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
			!all defect types in the defect list)
			defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
			do while(associated(defectCurrent))
				if(defectCurrent%num /= 0) then
				
					call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber,&
						defectTemp, defectCurrent%defectType)
				
				end if
					
				defectCurrent=>defectCurrent%next
				
			end do

			!See above for explanation of this section
			defectUpdateNext=>defectUpdateCurrent%next
			if(associated(defectUpdateNext)) then
				if(defectUpdateNext%cellNumber==defectUpdateCurrent%cellNumber) then
					
					call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber, &
						defectTemp, defectUpdateNext%defectType)
						
				end if
			end if
			
			!Diffusion reactions
			do j=1,6

				!Add diffusion reactions from this cell into neighboring cells
				call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber, &
					cascadeConnectivity(j, defectUpdateCurrent%cellNumber),myProc%taskid, myProc%taskid,j,defectTemp)

				!Add diffusion reactions from neighboring cells into this cell. (both directions)
			!	if(cascadeConnectivity(j, defectUpdateCurrent%cellNumber) /= 0) then
			!		call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber, cascadeConnectivity(j, defectUpdateCurrent%cellNumber), &
			!			defectUpdateCurrent%cellNumber,myProc%taskid, myProc%taskid,j,defectTemp)
			!	end if

			end do
		end if

	!if the defect is within the boundary mesh, only update the diffusion reactions for cells touching that boundary element
	else

		!If we are in this section, this means that a defect has changed in the boundary to this 
		!volume element but not IN the volume element. Therefore we only need to update one diffusion
		!reaction in one direction, the direction of the element with the changed defect. No defects
		!have changed in this volume element.

		!*******************************************************************************************
		!2015.04.06 Possible bug: is defectUpdateCurrent%dir pointing in the backwards direction?
		!(from neighboring cell TO local cell instead of vice versa) This may be the case because
		!defectUpdateCurrent was created in a different processor and passed here. Need to check
		!that it was initialized correctly.
		!*******************************************************************************************
		
		if(defectUpdateCurrent%neighbor==-1) then
			write(*,*) 'error neighbor not assigned for diffusion reactions into boundary'
			call MPI_ABORT(comm,ierr)
		end if
		
		if(polycrystal=='yes') then
		
			localGrainID=myMesh(defectUpdateCurrent%neighbor)%material
			
			neighborGrainID=myBoundary(defectUpdateCurrent%cellNumber,defectUpdateCurrent%dir)%material
			
			if(localGrainID==neighborGrainID) then
				
				!Allow diffusion between elements in the same grain
				call addDiffusionReactions(defectUpdateCurrent%neighbor, defectUpdateCurrent%cellNumber, myProc%taskid, &
					defectUpdateCurrent%proc, defectUpdateCurrent%dir, defectTemp)
			
			else
			
				!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
				call addDiffusionReactions(defectUpdateCurrent%neighbor, 0, myProc%taskid, -1, defectUpdateCurrent%dir, &
					defectTemp)
			
			endif
		
		else
		
			call addDiffusionReactions(defectUpdateCurrent%neighbor, defectUpdateCurrent%cellNumber, myProc%taskid, &
				defectUpdateCurrent%proc, defectUpdateCurrent%dir, defectTemp)
			
		end if
		
	end if
	defectUpdateCurrent=>defectUpdateCurrent%next

end do

!Deallocate defectUpdate (memory release)
defectUpdatePrev=>defectUpdate
defectUpdateCurrent=>defectUpdate%next
do while(associated(defectUpdatePrev))
	deallocate(defectUpdatePrev%defectType)
	deallocate(defectUpdatePrev)
	defectUpdatePrev=>defectUpdateCurrent
	if(associated(defectUpdateCurrent)) then
		defectUpdateCurrent=>defectUpdateCurrent%next
	endif
end do

end subroutine
