! $Header: /home/CVS//srscd/src/kMC_subroutines.f90,v 1.12 2015/10/19 19:43:22 aydunn Exp $
!*****************************************************************************************
!>double precision generate timestep - chooses a timestep using random number and Monte Carlo algorithm
!!(this is a global timestep)
!*****************************************************************************************

double precision function GenerateTimestep()
use DerivedType
use mod_srscd_constants
use randdp
implicit none

double precision r1

!***************************************************************************************************
!Here, we generate the GLOBAL timestep using the kMC algorithm.
!Due to the fact that we update totalRate and maxRate after each step, we do not have to calculate
!the total reaction rate at this step.
!***************************************************************************************************

r1=dprand()
GenerateTimestep=1/maxRate*dlog(dble(1d0/r1))

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
use mod_srscd_constants
use randdp
implicit none

type(reaction), pointer :: reactionCurrent, reactionTemp
type(cascade), pointer :: cascadeCurrent
type(defect), pointer :: defectTemp
double precision r2, atemp, atemp_cell, r2timesa
integer i, j
double precision totalRateCells

!***************************************************************************************************
!Choose from reactions within the local mesh, including null event, according to kMC algorithm
!***************************************************************************************************

atemp=0d0
atemp_cell=0d0
r2=dprand()
r2timesa=r2*maxRate
nullify(CascadeCurrent)		!These are default pointed at nothing, indicating null event
nullify(reactionCurrent)	!These are default pointed at nothing, indicating null event

!***********************************************************************
!Here we choose from reactions within the coarse mesh
!***********************************************************************

!2015.04.02: TO DO: use totalRateVol(cell) to choose a volume element within which we choose a reaction

outer: do i=1,numCells
	
	!used to choose a volume element that the reaction occurs within
	atemp_cell=atemp_cell+totalRateVol(i)
	
	!search for which volume element we are choosing among
	if(r2timesa <= atemp_cell) then

		reactionCurrent=>reactionList(i)	!a reaction is chosen in this volume element, so we no longer nullify reactionCurrent
		
		!search for which reaction occurs in that volume element
		do while(associated(reactionCurrent))
			atemp=atemp+reactionCurrent%reactionRate
			if(r2timesa <= atemp) then
				exit outer			!exit both loops with reactionCurrent pointing to the randomly chosen reaction
			endif
			reactionCurrent=>reactionCurrent%next
		end do
		
	else
			atemp=atemp+totalRateVol(i)
	end if
end do outer

!***********************************************************************
!Here we choose from reactions within the fine meshes for active cascades
!***********************************************************************
if(r2timesa <= atemp) then
	!do nothing, we have already chosen our reaction
else
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
		if(reactionCurrent%numProducts==2 .AND. reactionCurrent%numReactants==0) then
			!Frenkel pair implantation
			numImplantEvents=numImplantEvents+1		!LOCAL number of implantation events. Total DPA calculated out using MPI_ALLREDUCE
		else if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts==1) then	!He implantation
			numHeImplantEvents=numHeImplantEvents+1
		end if
	else if(implantType=='Cascade') then
		if(reactionCurrent%numReactants==0 .OR. reactionCurrent%numReactants==-10) then
			if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
				numImplantEvents=numImplantEvents+1
			else if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts .NE. 0) then !He implantation
				numHeImplantEvents=numHeImplantEvents+1
			else
				write(*,*) 'Error reaction not allowed ', 'reactants', reactionCurrent%numReactants, &
					'products', reactionCurrent%numProducts, 'rate', reactionCurrent%reactionRate
				if(associated(CascadeCurrent)) then
					write(*,*) 'Error in cascade', CascadeCurrent%cascadeID
				end if
			end if
		end if
	else
		write(*,*) 'Error implantType'
	end if
	
	!for post processing: count annihilation reactions (hard coded)
	if(reactionCurrent%numReactants==2) then	!clustering reaction
		if(reactionCurrent%reactants(1,2) /= 0 .AND. reactionCurrent%reactants(2,3) /= 0) then	!V+SIA_mobile
			
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,2),reactionCurrent%reactants(2,3))
		
		else if(reactionCurrent%reactants(1,2) /= 0 .AND. reactionCurrent%reactants(2,4) /= 0) then	!V+SIA_sessile
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,2),reactionCurrent%reactants(2,4))
		
		else if(reactionCurrent%reactants(1,3) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then	!SIA_mobile+V
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,3),reactionCurrent%reactants(2,2))
		
		else if(reactionCurrent%reactants(1,4) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then !SIA_sessile+V
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,4),reactionCurrent%reactants(2,2))
		
		endif
	endif
	
	!for post processing: counting trapping and emission reactions from the grain boundary
	if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then	!diffusion/emission reaction
 !! if(reactionCurrent%cellNumber(1) > 0 .AND. reactionCurrent%cellNumber(2) > 0) then
                !write(*,*) 'myProc', myProc%taskid, 'cellNumber(2)', reactionCurrent%cellNumber(2)
		if(myMesh(reactionCurrent%cellNumber(1))%material == 1 .AND. &
			myMesh(reactionCurrent%cellNumber(2))%material == 2) then	!trapping on grain boundary
			
			numTrapV=numTrapV+reactionCurrent%reactants(1,2)
			numTrapSIA=numTrapSIA+reactionCurrent%reactants(1,3)
			
		else if(myMesh(reactionCurrent%cellNumber(1))%material==2 .AND. &
				myMesh(reactionCurrent%cellNumber(2))%material==1) then	!emission from grain boundary
				
			numEmitV=numEmitV+reactionCurrent%reactants(1,2)
			numEmitSIA=numEmitSIA+reactionCurrent%reactants(1,3)
			
		end if
!!  endif
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
use mod_srscd_constants
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
			numImplantEvents=numImplantEvents+1		!LOCAL number of implantation events. Total DPA calculated out using MPI_ALLREDUCE
		else if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts==1) then	!He implantation
			numHeImplantEvents=numHeImplantEvents+1
		end if
	else if(implantType=='Cascade') then
		if(reactionCurrent%numReactants==0 .OR. reactionCurrent%numReactants==-10) then
			if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
				numImplantEvents=numImplantEvents+1
			else if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts .NE. 0) then !He implantation
				numHeImplantEvents=numHeImplantEvents+1
			else
				write(*,*) 'Error reaction not allowed ', 'reactants', reactionCurrent%numReactants, &
					'products', reactionCurrent%numProducts, 'rate', reactionCurrent%reactionRate
				if(associated(CascadeCurrent)) then
					write(*,*) 'Error in cascade', CascadeCurrent%cascadeID
				end if
			end if
		end if
	else
		write(*,*) 'Error implantType'
	end if
	
	!for post processing: count annihilation reactions (hard coded)
	if(reactionCurrent%numReactants==2) then	!clustering reaction
		if(reactionCurrent%reactants(1,2) /= 0 .AND. reactionCurrent%reactants(2,3) /= 0) then	!V+SIA_mobile
			
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,2),reactionCurrent%reactants(2,3))
		
		else if(reactionCurrent%reactants(1,2) /= 0 .AND. reactionCurrent%reactants(2,4) /= 0) then	!V+SIA_sessile
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,2),reactionCurrent%reactants(2,4))
		
		else if(reactionCurrent%reactants(1,3) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then	!SIA_mobile+V
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,3),reactionCurrent%reactants(2,2))
		
		else if(reactionCurrent%reactants(1,4) /= 0 .AND. reactionCurrent%reactants(2,2) /= 0) then !SIA_sessile+V
		
			numAnnihilate=numAnnihilate+min(reactionCurrent%reactants(1,4),reactionCurrent%reactants(2,2))
		
		endif
	endif
endif

end subroutine


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

!*****************************************************************************************
!> Subroutine update defect list - updates defects according to reaction chosen
!!
!!This is the most involved and complex subroutine of the entire SRSCD algorithm.
!!This subroutine updates defects in the local mesh according to the reaction chosen, and
!!communicates with neighboring processors about defects that may have passed into 
!!a different processor as well as about defects that may have changed on the boundary.
!!It also creates a list of defects that have been updated, to inform the next subroutine (updateReactionList(defectUpdate))
!!which reactions to update.
!!
!!Input: reactionCurrent
!!Output: defectUpdateCurrent, CascadeCurrent
!*****************************************************************************************

subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
use mod_srscd_constants
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

!Used for communication between processors
!EDIT: 5/28/2014: delete some of these

!defect update counters
integer numUpdateLocal(6)			!<the number of defects being sent to each processor neighbor
integer numUpdateBndry(6)			!<the number of defects being recieved from each processor neighbor

integer count, numDefectsRecv
integer numUpdateLocalRecv(6), numUpdateBndryRecv(6)	!number of defects being recieved from each proc neighbor

!NOTE: this final step could be eliminated by keeping the global mesh in each local processor
!(thus each element could be identified as being part of the local mesh of one proc and the boundary of 
!any other procs)
integer numUpdateFinal(6), numUpdateFinalRecv(6)		!number of defects being sent/recieved in final update step

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
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
	
	subroutine chooseCascade(CascadeTemp)
	use DerivedType
	type(cascadeEvent), pointer :: CascadeTemp
	end subroutine
	
	subroutine initializeFineMesh(CascadeCurrent)
	use DerivedType
	type(cascade), pointer :: CascadeCurrent
	end subroutine
end interface

!initialize the trackers of number of defects to update
do 9 i=1,6
	numUpdateLocal(i)=0
	numUpdateBndry(i)=0
9 continue

if(associated(reactionCurrent)) then	!if we have not chosen a null event

	!***********************************************************************************************
	!
	!Cascade chosen
	!
	!Initialization of fine mesh, population wtih defects from coarse mesh, and addition of 
	!cascade defects into fine mesh.
	!
	!***********************************************************************************************
	
	
	if(reactionCurrent%numReactants==-10) then 
		
		if(meshingType=='adaptive') then
			
			!**************************************************
			!create fine mesh
			!populate with defects from coarse mesh element
			!update defectUpdate
			!update localBuffer if needed
			!populate with defects from cascade
			!**************************************************
			
			!***************************************************************************************
			!2014/06/16 NOTE: This code is taken from previous SRSCD code and is NOT YET edited for
			!inclusion in SRSCD_par.
			!
			!2014/07/18: This code is now fully integrated into SRSCD_par
			!***************************************************************************************
			call chooseCascade(CascadeTemp)	!choose one of the cascades from the list randomly
			cascadeDefectTemp=>cascadeTemp%ListOfDefects
			
			!This line causes an error because reactionCurrent%numProducts needs to remain =0
			!(unlike in prev. code where all reactions were updated at each step)
			!NOTE: is this even being used in the current code? Can we delete this line without
			!problems? It was used in prev SRSCD for reaction list update but that is being
			!taken care of now by defectUpdate, a pointed list.
			
			!reactionCurrent%numProducts=cascadeTemp%numDefectsTotal
			
			!write(*,*) 'cascade center cell', findCellWithCoordinates(coordinates)
			
			!initialize the cascade to be added
			CascadeCurrent=>ActiveCascades
			if(.NOT. associated(CascadeCurrent)) then		!no cascade fine meshes are currently in the system
				!write(*,*) 'initializing first cascade'
				allocate(ActiveCascades)
				!write(*,*) 'allocated activeCascades'
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
				CascadeCurrent%cascadeID=numImplantEvents
				!write(*,*) 'initialized ActiveCascades'
			else
				j=1
				do 540 while(associated(CascadeCurrent))
					CascadePrev=>CascadeCurrent
					CascadeCurrent=>CascadeCurrent%next
					j=j+1
				540 continue
				!write(*,*) 'initializing cascade', j
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
				CascadeCurrent%cascadeID=numImplantEvents
			endif
			
			!*******************************************************************
			!when we create a cascade in this coarse mesh element, we reduce its
			!volume by the cascade volume
			!*******************************************************************
			
			!update volume of couarse mesh element
			myMesh(cascadeCurrent%cellNumber)%volume=myMesh(cascadeCurrent%cellNumber)%volume-CascadeElementVol*numCellsCascade
			
			!If we add too many cascades to a coarse mesh element, the volume inside can become negative. If this 
			!happens, output an error message. (Only a danger in the case of very high DPA rates and
			!small coarse mesh elements)
			
			if(myMesh(cascadeCurrent%cellNumber)%volume .LE. 0d0) then
				write(*,*) 'Error negative coarse mesh volume'
			endif
			
			!*******************************************************************
			!initialize defect list and reaction list in CascadeCurrent. This includes defects
			!placed in the fine mesh from the coarse mesh.
			!*******************************************************************
			
			call initializeFineMesh(CascadeCurrent)
			
			!output initial mesh populations
!			if(myProc%taskid==MASTER) then
!				do 670 i=1,numCellsCascade
!					write(*,*) 'cascade cell', i
!					defectCurrent=>CascadeCurrent%localDefects(i)
!					do 671 while(associated(defectCurrent))
!						write(*,*) (defectCurrent%defectType(j),j=1,numSpecies), 'num', defectCurrent%num
!						defectCurrent=>defectCurrent%next
!					671 continue
!				670 continue
!			endif
			
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
			nullify(defectStoreList%next)
			defectStore=>defectStoreList
			
			!***************************************************************************************
			!2014.10.07: Updating cascade recombination step for increased efficiency
			!
			!Step 1: Create defectStoreList with defects in cascade only
			!Step 2: For each defect in the fine mesh, check whether it combines with cascade
			!Step 3: If yes, remove it from the fine mesh and combine it randomly with one defect in
			!	defectStoreList
			!Step 4: Implant defects in defectStoreList in fine mesh. Remember to correctly update
			!	defect update protocol.
			!***************************************************************************************

			!Step 1:
			do 590 i=1, cascadeTemp%numDefectsTotal
				
				!**************************************************************
				!Make sure coordinates of all defects in cascade are within fine mesh, using periodic BC's
				!**************************************************************
				coordinatesTemp(1)=cascadeDefectTemp%coordinates(1)
				if(coordinatesTemp(1) .GT. numxcascade*finelength/2d0) then
					coordinatesTemp(1)=coordinatesTemp(1)-numxcascade*finelength
				else if(coordinatesTemp(1) .LT. -numxcascade*finelength/2d0) then
					coordinatesTemp(1)=coordinatesTemp(1)+numxcascade*finelength
				endif
				
				coordinatesTemp(2)=cascadeDefectTemp%coordinates(2)
				if(coordinatesTemp(2) .GT. numycascade*finelength/2d0) then
					coordinatesTemp(2)=coordinatesTemp(2)-numycascade*finelength
				else if(coordinatesTemp(2) .LT. -numycascade*finelength/2d0) then
					coordinatesTemp(2)=coordinatesTemp(2)+numycascade*finelength
				endif
				
				coordinatesTemp(3)=cascadeDefectTemp%coordinates(3)
				if(coordinatesTemp(3) .GT. numzcascade*finelength/2d0) then
					coordinatesTemp(3)=coordinatesTemp(3)-numzcascade*finelength
				else if(coordinatesTemp(3) .LT. -numzcascade*finelength/2d0) then
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
				
				do 556 j=1,numSpecies
					defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
				556 continue
				defectStore%cellNumber=cellNumber
				defectStore%num=1
				
				cascadeDefectTemp=>cascadeDefectTemp%next
			590 continue
			
			do 570 j=1,numCellsCascade
				defectCurrent=>CascadeCurrent%localDefects(j)
				
				defectPrev=>defectCurrent
				defectTemp=>defectCurrent%next
				do 550 while(associated(defectTemp))
					
					!Step 2:
					!boolean variable, if FALSE then no combination with defectCurrent, if TRUE then combination
					mixingTemp=0
					do 571 k=1,defectTemp%num
						combineBoolean=cascadeMixingCheck()
						
						if(combineBoolean .eqv. .TRUE.) then
							mixingTemp=mixingTemp+1
						endif
					571 continue
					
					!Step 3:
					do 572 k=1,mixingTemp
						mixingEvents=mixingEvents+1
						
						!Choose which defect will be combined wtih defectTemp
						defectStore=>defectStoreList%next
						
						atemp=0d0
						r1=dprand()
						
						!Move defectStore through defectStoreList until defect chosen for recombination
						do 573 i=1,cascadeTemp%numDefectsTotal
							atemp=atemp+1d0/dble(cascadeTemp%numDefectsTotal)
							
							if(atemp .GT. r1) then
								exit
							endif
							
							defectStore=>defectStore%next
						573 continue
						
						if( .NOT. associated(defectStore)) then
							write(*,*) 'Error DefectStore not associated in cascade mixing'
						endif
						
						!***********************************************************************
						!Hard coded: use defect combination rules to combine defectStore and defecTemp
						!
						!6/19/2014: These rules have been transported to a separate subroutine
						!in order to facilitate hard-coding and keep this subroutine clean.
						!
						!This subroutine is placed in reactionRates.f90 in accordance with the
						!other clustering combination rules that are hard-coded.
						!***********************************************************************
						
						call defectCombinationRules(defectStore%defectType, defectTemp)
						
						!NEXT: remove DefectTemp from fine mesh
						
						if(defectTemp%num==0) then	!error
							write(*,*) 'error defect num zero combining with cascade defect'
						else if(defectTemp%num==1) then	!remove that defect from the defect list
						
							if(associated(defectPrev)) then !we are not at the beginning of the list
								if(associated(defectTemp%next)) then	!we are not at the end of the list
									defectPrev%next=>defectTemp%next
									deallocate(defectTemp%defectType)
									deallocate(defectTemp)
									defectTemp=>defectPrev
								else
									nullify(defectPrev%next)
									deallocate(defectTemp%defectType)
									deallocate(defectTemp)
									defectTemp=>defectPrev
								endif
							else
								defectTemp%num=0 !beginning of the list is never deleted
							endif
						else
							defectTemp%num=defectTemp%num-1
						endif
	
					572 continue
					
					defectPrev=>defectTemp
					defectTemp=>defectTemp%next
					
				550 continue
			570 continue						
				

!			do 590 i=1, cascadeTemp%numDefectsTotal
			
!				!**************************************************************
!				!Make sure coordinates of all defects in cascade are within fine mesh, using periodic BC's
!				!**************************************************************
!				coordinatesTemp(1)=cascadeDefectTemp%coordinates(1)
!				if(coordinatesTemp(1) .GT. numxcascade*finelength/2d0) then
!					coordinatesTemp(1)=coordinatesTemp(1)-numxcascade*finelength
!				else if(coordinatesTemp(1) .LT. -numxcascade*finelength/2d0) then
!					coordinatesTemp(1)=coordinatesTemp(1)+numxcascade*finelength
!				endif
				
!				coordinatesTemp(2)=cascadeDefectTemp%coordinates(2)
!				if(coordinatesTemp(2) .GT. numycascade*finelength/2d0) then
!					coordinatesTemp(2)=coordinatesTemp(2)-numycascade*finelength
!				else if(coordinatesTemp(2) .LT. -numycascade*finelength/2d0) then
!					coordinatesTemp(2)=coordinatesTemp(2)+numycascade*finelength
!				endif
				
!				coordinatesTemp(3)=cascadeDefectTemp%coordinates(3)
!				if(coordinatesTemp(3) .GT. numzcascade*finelength/2d0) then
!					coordinatesTemp(3)=coordinatesTemp(3)-numzcascade*finelength
!				else if(coordinatesTemp(3) .LT. -numzcascade*finelength/2d0) then
!					coordinatesTemp(3)=coordinatesTemp(3)+numzcascade*finelength
!				endif	
				
!				!find position cascade is implanted in fine mesh
!				cellNumber=findCellWithCoordinatesFineMesh(coordinatesTemp)	
				
!				nullify(defectPrev)
!				do 593 j=1,numSpecies
!					products(j)=cascadeDefectTemp%defectType(j)
!				593 continue
				
!				!write(*,*) 'Inserting defect', (products(j), j=1,numSpecies)
!				!write(*,*) 'Fine Mesh Cell Number', cellNumber, 'Coarse Mesh Cell Number', CascadeCurrent%cellNumber
				
!				!***********************************************************************************
!				!Loop through defects in the fine mesh (existing defects that were placed there during
!				!the fine mesh initialization) and decide if they should combine with this defect
!				!in the cascade. Do this for each defect in the fine mesh once for each cascade defect.
!				!
!				!This is a computationally expensive procedure and slows the simulation down at high
!				!doses. A more efficient way of carrying this out may be helpful in the long run.
!				!
!				!Future work: only decide if each defect in the fine mesh is combined with the 
!				!cascade once per defect in the fine mesh (not once per defect in the fine mesh
!				!per defect in the cascade), then assign each defect chosen this way to a defect
!				!in the cascade randomly. This will take our computation from n*m to n+m.
!				!***********************************************************************************
				
!				do 570 j=1,numCellsCascade
!					defectCurrent=>CascadeCurrent%localDefects(j)
					
!					defectPrev=>defectCurrent
!					defectTemp=>defectCurrent%next
!					do 550 while(associated(defectTemp))
						
!						!boolean variable, if FALSE then no combination with defectCurrent, if TRUE then combination
!						mixingTemp=0
!						do 571 k=1,defectTemp%num
!							combineBoolean=cascadeMixingCheck()
							
!							if(combineBoolean .eqv. .TRUE.) then
!								mixingTemp=mixingTemp+1
!							endif
!						571 continue
						
!						do 572 k=1,mixingTemp
!							mixingEvents=mixingEvents+1
!							!combine cascade defect with defectTemp:
							
!							!STEP 0: calculate the defect created by cascade defect+defectTemp
!							!For now, don't worry about He-SIA cluster creation (cascades only)
							
!							!***********************************************************************
!							!Hard coded: defect combination rules
!							!
!							!6/19/2014: These rules have been transported to a separate subroutine
!							!in order to facilitate hard-coding and keep this subroutine clean.
!							!
!							!This subroutine is placed in reactionRates.f90 in accordance with the
!							!other clustering combination rules that are hard-coded.
!							!***********************************************************************
							
!							call defectCombinationRules(products, defectTemp)
							
!							!STEP 1: remove DefectTemp from fine mesh
							
!							if(defectTemp%num==0) then	!error
!								write(*,*) 'error defect num zero combining with cascade defect'
!							else if(defectTemp%num==1) then	!remove that defect from the defect list
							
!								if(associated(defectPrev)) then !we are not at the beginning of the list
!									if(associated(defectTemp%next)) then	!we are not at the end of the list
!										defectPrev%next=>defectTemp%next
!										deallocate(defectTemp%defectType)
!										deallocate(defectTemp)
!										defectTemp=>defectPrev
!									else
!										nullify(defectPrev%next)
!										deallocate(defectTemp%defectType)
!										deallocate(defectTemp)
!										defectTemp=>defectPrev
!									endif
!								else
!									defectTemp%num=0 !beginning of the list is never deleted
!								endif
!							else
!								defectTemp%num=defectTemp%num-1
!							endif
		
!							!***********************************************************
!							!Note: in order to not include inter-mixing between cascade defects, 
!							!we do not carry out step 2 until step 0-1 have been done for all
!							!defects in cascade
!							!***********************************************************	
!						572 continue			
		
!						defectPrev=>defectTemp
!						defectTemp=>defectTemp%next
!					550 continue
!				570 continue
			
!				!Create list of defects that need to be added to fine mesh (after casacade-fine mesh
!				!mixing has been taken into account)
				
!				allocate(defectStore%next)
!				defectStore=>defectStore%next
!				allocate(defectStore%defectType(numSpecies))
!				nullify(defectStore%next)
				
!				do 556 j=1,numSpecies
!					defectStore%defectType(j)=products(j)
!				556 continue
!				defectStore%cellNumber=cellNumber
!				defectStore%num=1
!				nullify(defectTemp)
				
!				cascadeDefectTemp=>cascadeDefectTemp%next
!			590 continue
			
			defectStore=>defectStoreList%next

			!add all cascade mixing products and normal cascade products to the system
			do 557 while(associated(defectStore))
			
				!**************************************************************
				!Temporary debug: checking that all defects added are admissable defects
				!**************************************************************
				count=0
				do 1000 j=1,numSpecies
					if(defectStore%defectType(j) .NE. 0) then
						count=count+1
					endif
				1000 continue
				if(count .GT. 2) then
					write(*,*) 'Error adding unadmissible defect to cascade'
					write(*,*) defectStore%defectType
					call MPI_ABORT(MPI_COMM_WORLD,ierr)
				endif
				!**************************************************************
				!Temporary debug: checking that all defects added are admissable defects
				!**************************************************************
			
				defectCurrent=>CascadeCurrent%localDefects(defectStore%cellNumber)
				!STEP 2: add DefectTemp+products to the system
				count=0
				do 581 j=1,numSpecies
					if(defectStore%defectType(j)==0) then
						count=count+1
					endif
				581 continue
				if(count==numSpecies) then
					!do nothing, we have total annihilation
				else
					
					nullify(defectPrev)
					call findDefectInList(defectCurrent, defectPrev, defectStore%defectType)
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						
						count=0
						do 582 j=1,numSpecies
							if(defectCurrent%defectType(j)==defectStore%defectType(j)) then
								count=count+1
							endif
						582 continue
						if(count==numSpecies) then
							defectCurrent%num=defectCurrent%num+1
							!write(*,*) 'products', products
							!write(*,*) 'DefectCurrent type', defectCurrent%defectType
						else		!if the defect is to be inserted in the list
							!write(*,*) 'products', products
							!write(*,*) 'DefectCurrent type', defectCurrent%defectType
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectStore%cellNumber
							defectPrev%num=1
							do 54 j=1,numSpecies
								defectPrev%defectType(j)=defectStore%defectType(j)
							54 continue
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
						do 555 j=1,numSpecies
							defectPrev%defectType(j)=defectStore%defectType(j)
						555 continue
					endif
					
				endif
				defectStore=>defectStore%next
			
			557 continue
			
			!Final step: set up defectUpdate for all defects in fine mesh
			do 575 i=1,numCellsCascade
				defectTemp=>CascadeCurrent%localDefects(i)%next
				
				do 576 while(associated(defectTemp))
				
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
					do 577 j=1,numSpecies
						defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
					577 continue
					defectUpdateCurrent%cellNumber=i
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID	!identify which cascade this defect is inside
					defectUpdateCurrent%num=1
					nullify(defectUpdateCurrent%next)
					
					defectTemp=>defectTemp%next
				
				576 continue
				
			575 continue
			
			!*******************************************************************
			!memory erase: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
			!*******************************************************************
			defectStore=>defectStoreList
			
			do 558 while(associated(defectStore))
				defectTemp=>defectStore
				defectStore=>defectStore%next
				deallocate(defectTemp%defectType)
				deallocate(defectTemp)
			558 continue
			
			nullify(defectStore)
			nullify(defectTemp)
			nullify(defectStoreList)
			nullify(defectCurrent)
			nullify(defectPrev)	
		
!			output final defect populations after cascade implantation
!			if(myProc%taskid==MASTER) then
!				do 594 i=1,numCellsCascade
!					write(*,*) 'cascade cell', i
!					defectCurrent=>CascadeCurrent%localDefects(i)
!					do 591 while(associated(defectCurrent))
!						write(*,*) (defectCurrent%defectType(j),j=1,numSpecies), 'num', defectCurrent%num
!						defectCurrent=>defectCurrent%next
!					591 continue
!				594 continue
!				read(*,*)
!			endif
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
			nullify(defectStoreList%next)
			defectStore=>defectStoreList
			
			!***************************************************************************************
			!2014.10.07: Updating cascade recombination step for increased efficiency
			!
			!Step 1: Create defectStoreList with defects in cascade only
			!Step 2: For each defect in the fine mesh, check whether it combines with cascade
			!Step 3: If yes, remove it from the fine mesh and combine it randomly with one defect in
			!	defectStoreList
			!Step 4: Implant defects in defectStoreList in fine mesh. Remember to correctly update
			!	defect update protocol.
			!***************************************************************************************

			!Step 1:
			do 1590 i=1, cascadeTemp%numDefectsTotal
				!Create list of defects that need to be added to cell (after casacade
				!mixing has been taken into account)
				
				allocate(defectStore%next)
				defectStore=>defectStore%next
				allocate(defectStore%defectType(numSpecies))
				nullify(defectStore%next)
				
				do 1556 j=1,numSpecies
					defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
				1556 continue
				defectStore%cellNumber=reactionCurrent%cellNumber(1)
				defectStore%num=1
				
				cascadeDefectTemp=>cascadeDefectTemp%next
			1590 continue
			
			defectCurrent=>defectList(reactionCurrent%cellNumber(1))
			
			defectPrev=>defectCurrent
			defectTemp=>defectCurrent%next
			
			if(cascadeVolume .GT. 0d0) then
			
				do 1550 while(associated(defectTemp))
					
					!Step 3:
					do 1572 k=1,defectTemp%num
						mixingEvents=mixingEvents+1
						
						!Choose which defect in the cascade will be combined with defectTemp already presetn in the cell
						defectStore=>defectStoreList%next
						
						atemp=0d0
						r1=dprand()
						
						!Move defectStore through defectStoreList until defect chosen for recombination
						do 1573 i=1,cascadeTemp%numDefectsTotal
							atemp=atemp+1d0/dble(cascadeTemp%numDefectsTotal)
							
							if(atemp .GT. r1) then
								exit
							endif
							
							defectStore=>defectStore%next
						1573 continue
						
						if( .NOT. associated(defectStore)) then
							write(*,*) 'Error DefectStore not associated in cascade mixing'
						endif
						
						!***********************************************************************
						!Hard coded: use defect combination rules to combine defectStore and defecTemp
						!
						!6/19/2014: These rules have been transported to a separate subroutine
						!in order to facilitate hard-coding and keep this subroutine clean.
						!
						!This subroutine is placed in reactionRates.f90 in accordance with the
						!other clustering combination rules that are hard-coded.
						!***********************************************************************
						
						call defectCombinationRules(defectStore%defectType, defectTemp)
						
						!NEXT: remove DefectTemp from fine mesh
						
						if(defectTemp%num==0) then	!error
							write(*,*) 'error defect num zero combining with cascade defect'
						else if(defectTemp%num==1) then	!remove that defect from the defect list
						
							if(associated(defectPrev)) then !we are not at the beginning of the list
								if(associated(defectTemp%next)) then	!we are not at the end of the list
									defectPrev%next=>defectTemp%next
									deallocate(defectTemp%defectType)
									deallocate(defectTemp)
									defectTemp=>defectPrev
								else
									nullify(defectPrev%next)
									deallocate(defectTemp%defectType)
									deallocate(defectTemp)
									defectTemp=>defectPrev
								endif
							else
								defectTemp%num=0 !beginning of the list is never deleted
							endif
						else
							defectTemp%num=defectTemp%num-1
						endif
	
					1572 continue
					
					defectPrev=>defectTemp
					defectTemp=>defectTemp%next
					
				1550 continue
			
			endif
			
			defectStore=>defectStoreList%next

			!add all cascade mixing products and normal cascade products to the system
			do 1557 while(associated(defectStore))
			
				defectCurrent=>defectList(defectStore%cellNumber)
				!STEP 2: add DefectTemp+products to the system
				count=0
				do 1581 j=1,numSpecies
					if(defectStore%defectType(j)==0) then
						count=count+1
					endif
				1581 continue
				if(count==numSpecies) then
					!do nothing, we have total annihilation
				else
					! Updatin defect list with defects contained in cascade defect list
					nullify(defectPrev)
					call findDefectInList(defectCurrent, defectPrev, defectStore%defectType)
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						
						count=0
						do 1582 j=1,numSpecies
							if(defectCurrent%defectType(j)==defectStore%defectType(j)) then
								count=count+1
							endif
						1582 continue
						if(count==numSpecies) then
							defectCurrent%num=defectCurrent%num+1
							!write(*,*) 'products', products
							!write(*,*) 'DefectCurrent type', defectCurrent%defectType
						else		!if the defect is to be inserted in the list
							!write(*,*) 'products', products
							!write(*,*) 'DefectCurrent type', defectCurrent%defectType
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=defectStore%cellNumber
							defectPrev%num=1
							do 1554 j=1,numSpecies
								defectPrev%defectType(j)=defectStore%defectType(j)
							1554 continue
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
						do 1555 j=1,numSpecies
							defectPrev%defectType(j)=defectStore%defectType(j)
						1555 continue
					endif
					
				endif
				defectStore=>defectStore%next
			
			1557 continue
			
			!Final step: set up defectUpdate for all defects in the cell
			defectTemp=>defectList(reactionCurrent%cellNumber(1))%next
			
			do 1576 while(associated(defectTemp))
			
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
				do 1577 j=1,numSpecies
					defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
				1577 continue
				defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(1)
				defectUpdateCurrent%neighbor=0	!not pointed at a different proc
				defectUpdateCurrent%cascadeNumber=0	!identify which cascade this defect is inside: 0 for coarse mesh
				defectUpdateCurrent%num=1
				nullify(defectUpdateCurrent%next)
				
				defectTemp=>defectTemp%next
			
			1576 continue
							
			!*******************************************************************
			!memory erase: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
			!*******************************************************************
			defectStore=>defectStoreList
			
			do 1558 while(associated(defectStore))
				defectTemp=>defectStore
				defectStore=>defectStore%next
				deallocate(defectTemp%defectType)
				deallocate(defectTemp)
			1558 continue
			
			nullify(defectStore)
			nullify(defectTemp)
			nullify(defectStoreList)
			nullify(defectCurrent)
			nullify(defectPrev)	
			
		else
			write(*,*) 'error meshingType'
		endif

	!***********************************************************************************************
	!
	!Defect update for reactions within the fine mesh.
	!
	!Typically, no communication will be carried out if a reaction is chosen within the fine mesh.
	!The exception is when a defect diffuses from the fine mesh to the coarse mesh.
	!
	!***********************************************************************************************
		
	else if(associated(CascadeCurrent)) then
		
	!Reactions in the fine mesh
		!write(86,*) 'Fine Mesh Reaction Chosen'
		!create buffers: size greater than max size needed (at most numReactants and numProducts to change)
		allocate(localBuffer(6,reactionCurrent%numReactants+reactionCurrent%numProducts,numSpecies+3))
		allocate(bndryBuffer(6,reactionCurrent%numReactants+reactionCurrent%numProducts,numSpecies+2))
		
		!Remove reactants from the system
		!if(myProc%taskid==MASTER) write(*,*) 'Removing reactants'
		!write(86,*) 'Removing reactants'
		do 700 i=1, reactionCurrent%numReactants
			
			!create a new element in defectUpdate and assign all variables except for num (will do later)
			allocate(defectUpdateCurrent%next)
			defectUpdateCurrent=>defectUpdateCurrent%next
			defectUpdateCurrent%proc=reactionCurrent%taskid(i)
			defectUpdateCurrent%dir=0	!not pointed at a different proc
			allocate(defectUpdateCurrent%defectType(numSpecies))
			
			do 701 j=1,numSpecies
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(i,j)
			701 continue
			
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID	
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i))
			
			do 702 while(associated(defectCurrent)) !search for reactant defect
				same=0
				do 703 j=1,numSpecies
					if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
						same=same+1
					endif
				703 continue
				if(same==numSpecies) then
					exit
				endif
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
			702 continue
			
			!No need to update local buffer - cascade is always totally contained within the coarse mesh (one processor)
			
			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there fine'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'rate', reactionCurrent%reactionRate
				write(*,*) 'cascade number', cascadeCurrent%cascadeID
				
				call DEBUGPrintDefects(0)
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
			
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
				DefectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
				defectUpdateCurrent%num=-1
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then 							!removing only defect from cell i (single helium) - this is redundant but will keep for now
				DefectList(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) 'trying to remove a defect that isnt there'
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
			else
				!decrease the number of defects by 1 if the number of defects is greater than 1
				defectCurrent%num=defectCurrent%num-1 !remove the defect from the system instead of the entire entry in the list
				defectUpdateCurrent%num=-1	!tell updateReactionList the new number of defects in the cell
			endif
		700 continue
		
		!Add products to system
		!if(myProc%taskid==MASTER) write(*,*) 'Adding products'
		!***********************************************************************
		!Here, I will assume cubic cells. If a defect migrates from one cell to another,
		!the defect will be given a percent chance of removal from the system
		!based on the cell size (%chance of removal = cell size/grain size)
		!This is to replicate the OKMC practice of removing SIA clusters after they
		!migrate 1 um.
		!***********************************************************************
		!write(86,*) 'Adding Products'
		flag=.FALSE.
		if(grainBoundaryToggle=='yes') then	
		
		!we have toggled the use of grain boundaries to remove defects from system
		
			if(reactionCurrent%numReactants==1 .and. reactionCurrent%numProducts==1) then 
				diffusionRandom=dprand()
				
				!randomly choose whether to remove this defect from the system according to the mean free path (parameter) and the
				!length of the volume element that the defect is currently in
				if(diffusionRandom .LE. fineLength/meanFreePath) then
					flag=.TRUE.
					!reactionCurrent%numProducts=0 !remove that defect from system
					!we can no longer remove the defect in this way because changing the reaction list
					!is no longer allowed (reaction list is not re-created at each step)
				endif		
			endif
		endif
	
		!adding products to the system
		
		if(flag .eqv. .FALSE.) then		
		
		!if the diffusion reaction defect did not get removed from the system
			
			do 704 i=1, reactionCurrent%numProducts
				
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
					do 705 j=1,numSpecies
						products(j)=reactionCurrent%products(i,j)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(i,j)
					705 continue
					
					if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants)==0) then
						
						!diffusion from fine mesh to coarse mesh: point defectCurrent at coarse mesh defect list
					
						defectCurrent=>defectList(CascadeCurrent%cellNumber)
						
						!defectUpdateCurrent tells updateReactionList to update reactions within the coarse mesh
						
						defectUpdateCurrent%cascadeNumber=0
						defectUpdateCurrent%cellNumber=CascadeCurrent%cellNumber
						
					else
						
						!diffusion within the fine mesh: point defectCurrent at fine mesh defect list
						
						defectCurrent=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
						
						!defectUpdateCurrent tells updateReactionList to update reactions within fine mesh
						
						defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID
						defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					
					endif
	
					! this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
					
					if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants)==0) then
						!First update local buffer if needed
						!write(86,*) 'Fine to coarse reaction'
						do 706 j=1,6
							do 707 k=1,myMesh(cascadeCurrent%cellNumber)%numNeighbors(j)
								if(myMesh(cascadeCurrent%cellNumber)%neighborProcs(j,k)==-1) then
									!do nothing, free surface
								else if(myMesh(cascadeCurrent%cellNumber)%neighborProcs(j,k) .NE. myProc%taskid) then	
								
									!neighboring element not in this proc
								
									numUpdateLocal(j)=numUpdateLocal(j)+1
									
									!species
									do 708 l=1,numSpecies
										localBuffer(j,numUpdateLocal(j),l)=reactionCurrent%products(i,l)
									708 continue
									
									!cell number in local mesh
									localBuffer(j,numUpdateLocal(j),numSpecies+1)=cascadeCurrent%cellNumber
									
									!cell number of neighbor (in different proc)
									localBuffer(j,numUpdateLocal(j),numSpecies+3)=myMesh(cascadeCurrent%cellNumber)%neighbors(j,k)
									
									!number of defects to be increased by 1
									localBuffer(j,numUpdateLocal(j),numSpecies+2)=1
									
		!							if(myProc%taskid==4) then
		!								if(j==5) then
		!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
		!								else if(j==6) then
		!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
		!								endif
		!							endif
								endif
							707 continue
						706 continue
					
					endif
					
					!Next update defects
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						
						do 709 j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						709 continue
						
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
							
							do 710 j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(i,j)
							710 continue
							
							!if inserted defect is in the middle of the list, point it to the next item in the list
							
							defectPrev%next=>defectCurrent 
						endif
					else 			
						
						!add a defect to the end of the list
						
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectUpdateCurrent%cellNumber
						defectPrev%num=1
						
						do 711 j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(i,j)
						711 continue

					endif
					
				endif
				
			704 continue
		
		else
			!If the defect has been removed due to grain boundary absorption, do nothing
		endif
	
	!***********************************************************************************************
	!
	! Defect update for reactions chosen in the coarse mesh.
	!
	! For defect changes in the boundary of other processors or in the boundary of this processor,
	! local and boundary buffers are created for communication with other processors.
	!
	!***********************************************************************************************

	else
		!if(myProc%taskid==MASTER) write(*,*) 'Coarse Mesh'
		!Reactions in the coarse mesh
		!create buffers: size greater than max size needed (at most numReactants and numProducts to change)
		allocate(localBuffer(6,reactionCurrent%numReactants+reactionCurrent%numProducts,numSpecies+3))
		allocate(bndryBuffer(6,reactionCurrent%numReactants+reactionCurrent%numProducts,numSpecies+2))
		!removing reactants from the system
		do 10 i=1, reactionCurrent%numReactants
			
			!create a new element in defectUpdate and assign all variables except for num (will do later)
			allocate(defectUpdateCurrent%next)
			defectUpdateCurrent=>defectUpdateCurrent%next
			defectUpdateCurrent%proc=reactionCurrent%taskid(i)
			defectUpdateCurrent%dir=0	!not pointed at a different proc
			allocate(defectUpdateCurrent%defectType(numSpecies))
			do 200 j=1,numSpecies
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(i,j)
			200 continue
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>DefectList(reactionCurrent%cellNumber(i))
			
			call findDefectInList(defectCurrent, defectPrev, reactionCurrent%reactants(i,:))
			
			!Check that defectCurrent is pointing towards the reactant
			same=0
			do 11 j=1,numSpecies
				if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
					same=same+1
				endif
			11 continue
			
			if(same .NE. numSpecies) then
				write(*,*) 'Proc', myProc%taskid, 'error defectCurrent not pointing to reactants'
				write(*,*) 'reactants', reactionCurrent%reactants(:,:)
				write(*,*) 'products', reactionCurrent%products(:,:)
				write(*,*) 'searching for', reactionCurrent%reactants(i,:)
				if(associated(defectCurrent)) write(*,*) 'DefectCurrent', defectCurrent%defectType
				if(associated(defectPrev)) write(*,*) 'DefectPrev', defectPrev%defectType
			endif
			
!			do 11 while(associated(defectCurrent)) !search for reactant defect
!				same=0
!				do 77 j=1,numSpecies
!					if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
!						same=same+1
!					endif
!				77 continue
!				if(same==numSpecies) then
!					exit
!				endif
!				defectPrev=>defectCurrent
!				defectCurrent=>defectCurrent%next
!			11 continue
			
			!First update local buffer if needed
			do 8 j=1,6
				do 7 k=1,myMesh(reactionCurrent%cellNumber(i))%numNeighbors(j)
					if(myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k) .NE. myProc%taskid .AND. &
						myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k) .NE. -1) then	!neighboring element not in this proc
						
						numUpdateLocal(j)=numUpdateLocal(j)+1
						
						do 6 l=1,numSpecies
							localBuffer(j,numUpdateLocal(j),l)=reactionCurrent%reactants(i,l)
						6 continue
						
						!Cell Number in local mesh
						localBuffer(j,numUpdateLocal(j),numSpecies+1)=reactionCurrent%cellNumber(i)
						
						!Cell number in boundary mesh
						localBuffer(j,numUpdateLocal(j),numSpecies+3)=myMesh(reactionCurrent%cellNumber(i))%neighbors(j,k)
						
						if(associated(defectCurrent)) then
							localBuffer(j,numUpdateLocal(j),numSpecies+2)=-1
						else
							write(*,*) 'Error tried to delete defect that wasnt there and send to neighboring proc'
						endif
						
						!write(*,*) 'localBuffer update: dir', j, 'proc', myProc%taskid, 'neighbor', &
						!	myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k)
						!write(*,*) 'buffer:', (localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2)
						
						!write(*,*) 'proc', myProc%taskid, 'neighborProc', myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k)
						!write(*,*) (localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+1)
	!					if(myProc%taskid==4) then
	!						if(j==5) then
	!							write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'remove'
	!						else if(j==6) then
	!							write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'remove'
	!						endif
	!					endif
					endif
				7 continue
			8 continue
			
			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there coarse'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'rate', reactionCurrent%reactionRate
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
			
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
				DefectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
				defectUpdateCurrent%num=-1
				write(*,*) 'Removing the first defect from the list'
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then 							!removing only defect from cell i (single helium) - this is redundant but will keep for now
				DefectList(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
				write(*,*) 'Removing the only defect from the list'
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) 'trying to remove a defect that isnt there'
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
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
		10 continue
	
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
				if(diffusionRandom .LE. myMesh(reactionCurrent%cellNumber(1))%length/meanFreePath) then
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
			
			do 12 i=1, reactionCurrent%numProducts
				
				if(reactionCurrent%taskid(i+reactionCurrent%numReactants) == -1) then
					!we have identified a defect that is going to be removed from the system via a free surface
					!Do nothing; no need to add this defect to any lists or update any reaction lists because of it
					
				else if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .LT. 0) then
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
					defectUpdateCurrent%cellNumber=chooseRandomCell()
					
					nullify(defectUpdateCurrent%next)
					defectUpdateCurrent%num=1		!used for updating reaction lists
					defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					
					!In coarse-to-fine reactions, the cascade number is stored in reactionCurrent%cellNumber
					!and identified with a negative. Therefore we make it positive again.
					defectUpdateCurrent%cascadeNumber=-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
					
					nullify(defectPrev)
					do 713 j=1,numSpecies
						products(j)=reactionCurrent%products(i,j)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(i,j)
					713 continue
					
					!Point CascadeCurrent at the correct cascade
					
					CascadeCurrent=>ActiveCascades
					
					do 712 while(associated(CascadeCurrent))
						
						if(CascadeCurrent%cascadeID==-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)) then
							exit
						endif
						
						CascadeCurrent=>CascadeCurrent%next
					
					712 continue
					
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
					
					!Next update defects
					
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do 714 j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						714 continue
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
							do 715 j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(i,j)
							715 continue
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=defectUpdateCurrent%cellNumber
						defectPrev%num=1
						do 716 j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(i,j)
						716 continue
					endif
					
				else if(reactionCurrent%taskid(i+reactionCurrent%numReactants) .NE. myProc%taskid .AND. &
					reactionCurrent%taskid(i+reactionCurrent%numReactants) .NE. -1) then
					
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
					
					do 225 j=1,6
						do 220 k=1,myMesh(reactionCurrent%cellNumber(1))%numNeighbors(j)
							if(flag .EQV. .FALSE.) then
								if(myMesh(reactionCurrent%cellNumber(1))%neighbors(j,k)==&
									reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
									myMesh(reactionCurrent%cellNumber(1))%neighborProcs(j,k)==&
									reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
			
									numUpdateBndry(j)=numUpdateBndry(j)+1
									do 224 l=1,numSpecies
										bndryBuffer(j,numUpdateBndry(j),l)=reactionCurrent%products(i,l)
									224 continue
									bndryBuffer(j,numUpdateBndry(j),numSpecies+1)=&
										reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
									bndryBuffer(j,numUpdateBndry(j),numSpecies+2)=1	!indicates we are adding one defect
									flag=.TRUE.
								endif
							endif
						220 continue
					225 continue
					
					!update the direction in defectUpdateCurrent (used for updating reaction lists)
					!create a new element in defectUpdate and assign all variables except for num and dir (will do later)
					do 20 j=1,6
	
						!if myBoundary at this cell has defectList associated, then it is a boundary element for a local element in this cell.
						!NOTE: this will not work in the case of a non-uniform mesh (because there could be more
						!than one local neighbor, will have to change to localNeighbor())
						if(myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%localNeighbor .NE. 0 &
							.AND. myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%proc == &
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
							defectUpdateCurrent%neighbor=myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))&
								%localNeighbor		
							
							do 136 l=1,numSpecies
								products(l)=reactionCurrent%products(i,l)
								defectUpdateCurrent%defectType(l)=reactionCurrent%products(i,l)
							136 continue
							
							!Find defect in defect list on myBoundary
							defectCurrent=>myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%defectList
	
							call findDefectInList(defectCurrent, defectPrev, products)
	
							!Insert defect in myBoundary(dir,cell)%defectList
							if(associated(defectCurrent)) then !if we aren't at the end of the list
								same=0
								do 37 l=1,numSpecies
									if(defectCurrent%defectType(l)==products(l)) then
										same=same+1
									endif
								37 continue
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
									do 38 l=1,numSpecies
										defectPrev%defectType(l)=reactionCurrent%products(i,l)
									38 continue
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
								do 39 l=1,numSpecies
									defectPrev%defectType(l)=reactionCurrent%products(i,l)
								39 continue
							endif
						endif
					20 continue
					
	!				do 225 m=1,numCells
	!					do 20 j=1,6
	!						do 220 k=1,myMesh(m)%numNeighbors(j)
	!							
	!							if(myMesh(m)%neighbors(j,k)==reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
	!								myMesh(m)%neighborProcs(j,k)==reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
	!								
	!								!only update the boundary buffer once even if the boundary element bounds the local
	!								!processor at only one place. That is because the boundary buffer becomes a local recv buffer
	!								!in the neighboring processor and we only want to update local defects once.
	!								if(flag .EQV. .FALSE.) then
	!									numUpdateBndry(j)=numUpdateBndry(j)+1
	!									do 224 l=1,numSpecies
	!										bndryBuffer(j,numUpdateBndry(j),l)=reactionCurrent%products(i,l)
	!									224 continue
	!									bndryBuffer(j,numUpdateBndry(j),numSpecies+1)=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
	!									bndryBuffer(j,numUpdateBndry(j),numSpecies+2)=1	!indicates we are adding one defect
	!									flag=.TRUE.
	!								endif
	!									
	!								!update the direction in defectUpdateCurrent (used for updating reaction lists)
	!								!create a new element in defectUpdate and assign all variables except for num and dir (will do later)
	!								allocate(defectUpdateCurrent%next)
	!								defectUpdateCurrent=>defectUpdateCurrent%next
	!								defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
	!								allocate(defectUpdateCurrent%defectType(numSpecies))
	!								defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
	!								nullify(defectUpdateCurrent%next)
	!								defectUpdateCurrent%dir=j
	!								defectUpdateCurrent%num=1	!used for updating reaction lists
	!								
	!								!This information is used to tell us which element in the local
	!								!processor bounds this element in myBoundary
	!								defectUpdateCurrent%neighbor=m
	!							
	!								do 136 l=1,numSpecies
	!									products(l)=reactionCurrent%products(i,l)
	!									defectUpdateCurrent%defectType(l)=reactionCurrent%products(i,l)
	!								136 continue
	!								
	!								!Find defect in defect list on myBoundary
	!								defectCurrent=>myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%defectList
	!
	!								call findDefectInList(defectCurrent, defectPrev, products)
	!
	!								!Insert defect in myBoundary(dir,cell)%defectList and update bndryBuffer
	!	
	!!								if(myProc%taskid==0) then
	!!									if(j==5) then
	!!										write(*,*) 'dir', j, 'migration into bdry',(products(l),l=1,numSpecies), &
	!!											reactionCurrent%cellNumber(i+reactionCurrent%numReactants), 1
	!!									else if(j==6) then
	!!										write(*,*) 'dir', j, 'migration into bdry',(products(l),l=1,numSpecies), &
	!!											reactionCurrent%cellNumber(i+reactionCurrent%numReactants), 1
	!!									endif
	!!								endif
	!								if(associated(defectCurrent)) then !if we aren't at the end of the list
	!									same=0
	!									do 37 l=1,numSpecies
	!										if(defectCurrent%defectType(l)==products(l)) then
	!											same=same+1
	!										endif
	!									37 continue
	!									if(same==numSpecies) then	!if the defect is already present in the list
	!										defectCurrent%num=defectCurrent%num+1
	!									else		!if the defect is to be inserted in the list
	!										nullify(defectPrev%next)
	!										allocate(defectPrev%next)
	!										nullify(defectPrev%next%next)
	!										defectPrev=>defectPrev%next
	!										allocate(defectPrev%defectType(numSpecies))
	!										defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
	!										defectPrev%num=1
	!										do 38 l=1,numSpecies
	!											defectPrev%defectType(l)=reactionCurrent%products(i,l)
	!										38 continue
	!										defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
	!									endif
	!								else 			!add a defect to the end of the list
	!									nullify(defectPrev%next)
	!									allocate(defectPrev%next)
	!									nullify(defectPrev%next%next)
	!									defectPrev=>defectPrev%next
	!									allocate(defectPrev%defectType(numSpecies))
	!									defectPrev%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
	!									defectPrev%num=1
	!									do 39 l=1,numSpecies
	!										defectPrev%defectType(l)=reactionCurrent%products(i,l)
	!									39 continue
	!								endif
	!							endif
	!						220 continue
	!					20 continue
	!				225 continue
					
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
					do 13 j=1,numSpecies
						products(j)=reactionCurrent%products(i,j)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(i,j)
					13 continue
					
					defectCurrent=>DefectList(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
	
					!First update local buffer if needed
					do 40 j=1,6
						do 41 k=1,myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%numNeighbors(j)
							if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(j,k)==-1) then
								!do nothing, free surface
							else if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(j,k) &
								.NE. myProc%taskid) then	!neighboring element not in this proc
							
								numUpdateLocal(j)=numUpdateLocal(j)+1
								
								!species
								do 42 l=1,numSpecies
									localBuffer(j,numUpdateLocal(j),l)=reactionCurrent%products(i,l)
								42 continue
								
								!cell number in local mesh
								localBuffer(j,numUpdateLocal(j),numSpecies+1)=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
								
								!cell number of neighbor (in different proc)
								localBuffer(j,numUpdateLocal(j),numSpecies+3)=myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))&
									%neighbors(j,k)
								
								!number of defects to be increased by 1
								localBuffer(j,numUpdateLocal(j),numSpecies+2)=1
								
	!							if(myProc%taskid==4) then
	!								if(j==5) then
	!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
	!								else if(j==6) then
	!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
	!								endif
	!							endif
							endif
						41 continue
					40 continue
					
					!Next update defects
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do 130 j=1,numSpecies
							if(defectCurrent%defectType(j)==products(j)) then
								same=same+1
							endif
						130 continue
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
							do 14 j=1,numSpecies
								defectPrev%defectType(j)=reactionCurrent%products(i,j)
							14 continue
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
						do 15 j=1,numSpecies
							defectPrev%defectType(j)=reactionCurrent%products(i,j)
						15 continue
					endif
					
				endif
				
			12 continue
		
		else
			!If the defect has been removed due to grain boundary absorption
		endif
	endif

endif	!if associated(reactionCurrent)

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

do 21 i=1,6

	if(myProc%procNeighbor(i) .NE. myProc%taskid) then

		call MPI_SEND(numUpdateLocal(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+i,MPI_COMM_WORLD, ierr)	!number of local defects that have changed on boundary of processor i
		call MPI_SEND(numUpdateBndry(i),1,MPI_INTEGER,myProc%procNeighbor(i),i+6,MPI_COMM_WORLD, ierr)		!number of defects in the mesh of processor i that have changed (diffusion only)
		
		!EDIT: need to only send the first numUpdateLocal(i) elements of localBuffer(i,:,:)
		
		if(numUpdateLocal(i) .NE. 0) then

			allocate(localBufferSend(numUpdateLocal(i),numSpecies+3))
			
			do 23 j=1,numUpdateLocal(i)
				do 24 k=1,numSpecies+3
					localBufferSend(j,k)=localBuffer(i,j,k)
				24 continue
			23 continue

			call MPI_SEND(localBufferSend,numUpdateLocal(i)*(numSpecies+3), MPI_INTEGER, &
				myProc%procNeighbor(i),i+12,MPI_COMM_WORLD,ierr)
!			if(myProc%taskid==4) then
!				write(*,*) 'dir', i, 'proc sent', myProc%procNeighbor(i), 'numUpdateLocal', numUpdateLocal(i)
!			endif

			deallocate(localBufferSend)

		endif
		
		if(numUpdateBndry(i) .NE. 0) then

			allocate(bndryBufferSend(numUpdateBndry(i),numSpecies+2))
			
			do 25 j=1,numUpdateBndry(i)
				do 26 k=1,numSpecies+2
					bndryBufferSend(j,k)=bndryBuffer(i,j,k)
				26 continue
			25 continue

			call MPI_SEND(bndryBufferSend,numUpdateBndry(i)*(numSpecies+2), MPI_INTEGER, &
				myProc%procNeighbor(i),i+18,MPI_COMM_WORLD,ierr)
			deallocate(bndryBufferSend)

		endif

	endif
21 continue

!**************
!Step 3: Recieve data about local/bdry defects that have changed and update defectList and myBoundary accordingly
!**************

do 102 i=1,6
	numUpdateFinal(i)=0
102 continue
totalLocalRecv=0

do 22 i=1,6

	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) .NE. myProc%taskid) then

		call MPI_RECV(numUpdateBndryRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+tag,MPI_COMM_WORLD,status,ierr)	!number of bndry defects that have changed (local to processor i, boundary here)
		call MPI_RECV(numUpdateLocalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+6,MPI_COMM_WORLD,status,ierr)	!number of local defects that have changed (bndry to processor i, local here)
		totalLocalRecv=totalLocalRecv+numUpdateLocalRecv(i)
	endif
22 continue

allocate(finalBuffer(6,totalLocalRecv,numSpecies+3))

do 221 i=1,6
	
	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) .NE. myProc%taskid) then
		
		if(numUpdateBndryRecv(i) .NE. 0) then
			
			!Read in defects to update in boundary
			allocate(bndryBufferRecv(numUpdateBndryRecv(i),numSpecies+3))
			call MPI_RECV(bndryBufferRecv,numUpdateBndryRecv(i)*(numSpecies+3),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+12,MPI_COMM_WORLD,status,ierr)
!			if(myProc%taskid==MASTER) then
!				write(*,*) 'dir', i, 'proc recvd', myProc%procNeighbor(i), 'numUpdateBndryRecv', numUpdateBndryRecv(i)
!			endif
			
			!Add defects in bndryBufferRecv to myBoundary()
			do 27 j=1,numUpdateBndryRecv(i)
!				if(myProc%taskid==0) then
!					if(i==5) then
!						write(*,*) 'dir', i, 'bdry recv',(bndryBufferRecv(j,m),m=1,numSpecies+2)
!					else if(i==6) then
!						write(*,*) 'dir', i, 'bdry recv',(bndryBufferRecv(j,m),m=1,numSpecies+2)
!					endif
!				endif
				
				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%procNeighbor(i)
				defectUpdateCurrent%dir=i	
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=bndryBufferRecv(j,numSpecies+1)
				defectUpdateCurrent%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)
				
				!This information is used to tell us which element in the local
				!processor bounds this element in myBoundary
				defectUpdateCurrent%neighbor=bndryBufferRecv(j,numSpecies+3)
				
				!point defectCurrent at the defect list in the correct cell of myBoundary
				defectCurrent=>myBoundary(i,bndryBufferRecv(j,numSpecies+1))%defectList
				if(.NOT. associated(defectCurrent)) then
					write(*,*) 'error myBoundary not allocated correctly'
					write(*,*) 'dir', i, 'cell', bndryBufferRecv(j,numSpecies+1)
					call MPI_ABORT(MPI_COMM_WORLD,ierr)
				endif
				
				nullify(defectPrev)
				
				do 28 k=1,numSpecies
					products(k)=bndryBufferRecv(j,k)
					defectUpdateCurrent%defectType(k)=bndryBufferRecv(j,k)	!used to update reaction lists
				28 continue
				
				!point DefectCurrent at the defect we are looking for (if it is there), otherwise
				!point defectCurrent after and defectPrev before where it should go in defectList.
!				if(myProc%taskid==MASTER .OR. myProc%taskid==3) then
!					write(*,*) 'proc', myProc%taskid
!					write(*,*) 'numUpdateBndryRecv', numUpdateBndryRecv(i)
!					write(*,*) 'dir', i, (bndryBufferRecv(j,k),k=1,numSpecies+2)
!				endif
				
				call findDefectInList(defectCurrent, defectPrev, products)

				!Next update defects in myBoundary
				if(associated(defectCurrent)) then !if we aren't at the end of the list
					
					same=0
					do 29 l=1,numSpecies
						if(defectCurrent%defectType(l)==products(l)) then
							same=same+1
						endif
					29 continue
					if(same==numSpecies) then	!if the defect is already present in the list
						defectCurrent%num=defectCurrent%num+bndryBufferRecv(j,numSpecies+2) !this will be +/- 1 only
						
						if(defectCurrent%num==0) then
							!delete this defect from the list in myBoundary
							defectPrev%next=>defectCurrent%next !remove that defect type from the system
							deallocate(defectCurrent%defectType)
							deallocate(defectCurrent)
							nullify(defectCurrent)
						endif
					else		!if the defect is to be inserted in the list
					
						if(bndryBufferRecv(j,numSpecies+2) == -1) then
							write(*,*) 'defectCurrent associated', defectCurrent%defectType, 'num', defectCurrent%num
							if(associated(defectPrev)) write(*,*) 'defectPrev', defectPrev%defectType, 'num', defectPrev%num
							write(*,*) 'error in defectUpdate negative defect numbers', myProc%taskid
							write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
							write(*,*) (bndryBufferRecv(j,k),k=1,numSpecies+2)
							
							!call MPI_ABORT(MPI_COMM_WORLD,ierr)
							
						else
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=bndryBufferRecv(j,numSpecies+1)
							defectPrev%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
							
							do 30 l=1,numSpecies
								defectPrev%defectType(l)=bndryBufferRecv(j,l)
							30 continue
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
												
					endif
				else 			!add a defect to the end of the list
					
					if(bndryBufferRecv(j,numSpecies+2) == -1) then
						write(*,*) 'defectCurrent not associated'
						write(*,*) 'error in defectUpdate negative defect numbers'
						write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
						write(*,*) (bndryBufferRecv(j,k),k=1,numSpecies+2)
						!call MPI_ABORT(MPI_COMM_WORLD,ierr)
					else
					
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=bndryBufferRecv(j,numSpecies+1)
						defectPrev%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
											
						do 31 l=1,numSpecies
							defectPrev%defectType(l)=bndryBufferRecv(j,l)
						31 continue
					endif

				endif
			27 continue
			deallocate(bndryBufferRecv)
		endif
		
		if(numUpdateLocalRecv(i) .NE. 0) then
			
			!Read in defects to update in local mesh
			allocate(localBufferRecv(numUpdateLocalRecv(i),numSpecies+2))
			call MPI_RECV(localBufferRecv,numUpdateLocalRecv(i)*(numSpecies+2),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+18,MPI_COMM_WORLD,status,ierr)

			!Add defects in localBufferRecv to defectList()
			do 32 j=1,numUpdateLocalRecv(i)
				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%taskid
				defectUpdateCurrent%dir=0	!not pointed at a different proc
				defectUpdateCurrent%neighbor=0	!not pointed at a different proc
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=localBufferRecv(j,numSpecies+1)
				defectUpdateCurrent%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)

				!point defectCurrent at the defect list in the correct cell of defectList
				defectCurrent=>defectList(localBufferRecv(j,numSpecies+1))
				nullify(defectPrev)
				
				do 33 k=1,numSpecies
					products(k)=localBufferRecv(j,k)
					defectUpdateCurrent%defectType(k)=localBufferRecv(j,k)
				33 continue
				
				!point DefectCurrent at the defect we are looking for (if it is there), otherwise
				!point defectCurrent after and defectPrev before where it should go in defectList.

				call findDefectInList(defectCurrent, defectPrev, products)

				!Next update defects in local mesh
				if(associated(defectCurrent)) then !if we aren't at the end of the list
					same=0
					do 34 l=1,numSpecies
						if(defectCurrent%defectType(l)==products(l)) then
							same=same+1
						endif
					34 continue
					if(same==numSpecies) then	!if the defect is already present in the list
						defectCurrent%num=defectCurrent%num+localBufferRecv(j,numSpecies+2)	!This will be +/- 1 only (should be only +1 here)
						
						if(defectCurrent%num==0) then
							write(*,*) 'Error zero defects in updateDefectList, step 3 of communication'
						endif
					else		!if the defect is to be inserted in the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=localBufferRecv(j,numSpecies+1)
						defectPrev%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only (should be only +1 here)
						do 35 l=1,numSpecies
							defectPrev%defectType(l)=localBufferRecv(j,l)
						35 continue
						defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
					endif
				else 			!add a defect to the end of the list
					nullify(defectPrev%next)
					allocate(defectPrev%next)
					nullify(defectPrev%next%next)
					defectPrev=>defectPrev%next
					allocate(defectPrev%defectType(numSpecies))
					defectPrev%cellNumber=localBufferRecv(j,numSpecies+1)
					defectPrev%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only (should be only +1 here)
					do 36 l=1,numSpecies
						defectPrev%defectType(l)=localBufferRecv(j,l)
					36 continue
				endif
				
				!*****************************************************
				!Prep for step 4: check to see if element updated has more than one neighbor proc
				!*****************************************************
				
				do 100 k=1,6
					do 101 l=1,myMesh(localBufferRecv(j,numSpecies+1))%numNeighbors(k)
						
						!if a neighbor of this element is not in this proc or in the proc that just communicated with it
						if (myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) .NE. myProc%taskid .AND. &
							myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) .NE. myProc%procNeighbor(i) .AND. &
							myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) .NE. -1) then
							
							!Add this defect to a final buffer
							numUpdateFinal(k)=numUpdateFinal(k)+1
							!finalBuffer(1,1,1)=100	!why is this here?
							do 103 m=1,numSpecies
								finalBuffer(k,numUpdateFinal(k),m)=localBufferRecv(j,m)
							103 continue
							finalBuffer(k,numUpdateFinal(k),numSpecies+1)=localBufferRecv(j,numSpecies+1)	!local cell number
							finalBuffer(k,numUpdateFinal(k),numSpecies+2)=localBufferRecv(j,numSpecies+2) !This should be +1 only
							finalBuffer(k,numUpdateFinal(k),numSpecies+3)=myMesh(localBufferRecv(j,numSpecies+1))%neighbors(k,l)	!cell number in neighboring proc
							
							if(localBufferRecv(j,numSpecies+2)==-1) then
								write(*,*) 'error: cell in boundary of multiple procs removing defect'
								call MPI_ABORT(MPI_COMM_WORLD,ierr)
							endif
!							if(myProc%taskid==4) then
!								if(k==5) then
!									write(*,*) 'dir', k, 'final buffer',(finalBuffer(k,numUpdateFinal(k),m),m=1,numSpecies+2), 'add'
!								else if(k==6) then
!									write(*,*) 'dir', k, 'final buffer',(finalBuffer(k,numUpdateFinal(k),m),m=1,numSpecies+2), 'add'
!								endif
!							endif
						endif

					101 continue
				100 continue
			32 continue

			deallocate(localBufferRecv)
		endif
		
	endif
221 continue

!*************
!Step 4: if a local defect is updated due to diffusion to the boundary of another processor,
!	and the local element is also in the boundary of a processor other than the one already 
!	communicated with, let all other processors know of local defects that have been updated
!*************

do 104 i=1,6
	if(myProc%procNeighbor(i) .NE. myProc%taskid) then
		call MPI_SEND(numUpdateFinal(i), 1, MPI_INTEGER, myProc%procNeighbor(i), i+24, MPI_COMM_WORLD, ierr)	!number of local defects that have changed on boundary of processor i
		
		if(numUpdateFinal(i) .NE. 0) then
			allocate(finalBufferSend(numUpdateFinal(i),numSpecies+3))
			
			do 105 j=1,numUpdateFinal(i)
				do 106 k=1,numSpecies+3
					finalBufferSend(j,k)=finalBuffer(i,j,k)
				106 continue
			105 continue
			
			call MPI_SEND(finalBufferSend,numUpdateFinal(i)*(numSpecies+3), MPI_INTEGER, &
				myProc%procNeighbor(i),i+30,MPI_COMM_WORLD,ierr)
			deallocate(finalBufferSend)
		endif
	endif
104 continue

!*************
!Step 5: update myBoundary again, based on information sent in step 4
!*************

do 107 i=1,6

	!We have to switch the tags on MPI_RECV in order for the correct send/recieve pair to be exchanged
	if(i==1 .OR. i==3 .OR. i==5) then
		tag=i+1
	else
		tag=i-1
	endif
	
	if(myProc%procNeighbor(i) .NE. myProc%taskid) then
		
		
		!number of bndry defects that have changed (local to processor i, boundary here)
		call MPI_RECV(numUpdateFinalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+24,MPI_COMM_WORLD,status,ierr)	
		
		if(numUpdateFinalRecv(i) .NE. 0) then
			
			!Read in defects to update in boundary
			allocate(finalBufferRecv(numUpdateFinalRecv(i),numSpecies+3))
			call MPI_RECV(finalBufferRecv,numUpdateFinalRecv(i)*(numSpecies+3),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+30,MPI_COMM_WORLD,status,ierr)

			!Add defects in finalBufferRecv to myBoundary()
			do 108 j=1,numUpdateFinalRecv(i)
				
!				if(myProc%taskid==0) then
!					if(i==5) then
!						write(*,*) 'dir', i, 'final recv',(finalBufferRecv(j,m),m=1,numSpecies+2)
!					else if(i==6) then
!						write(*,*) 'dir', i, 'final recv',(finalBufferRecv(j,m),m=1,numSpecies+2)
!					endif
!				endif
				
				if(finalBufferRecv(j,numSpecies+2)==-1) then
					write(*,*) 'error deleting defects in finalBufferRecv'
					!call MPI_ABORT(MPI_COMM_WORLD,ierr)
				else
					
					!create a new element in defectUpdate and assign all variables except for num (will do later)
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=myProc%procNeighbor(i)
					defectUpdateCurrent%dir=i	
					allocate(defectUpdateCurrent%defectType(numSpecies))
					defectUpdateCurrent%cellNumber=finalBufferRecv(j,numSpecies+1)	!cell number in boundary mesh
					defectUpdateCurrent%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
					defectUpdateCurrent%neighbor=finalBufferRecv(j,numSpecies+3) !cell number in myMesh
					defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
					nullify(defectUpdateCurrent%next)
					
					!point defectCurrent at the defect list in the correct cell of myBoundary
					defectCurrent=>myBoundary(i,finalBufferRecv(j,numSpecies+1))%defectList
					nullify(defectPrev)
					
					do 109 k=1,numSpecies
						products(k)=finalBufferRecv(j,k)
						defectUpdateCurrent%defectType(k)=finalBufferRecv(j,k)
					109 continue
					
					!point DefectCurrent at the defect we are looking for (if it is there), otherwise
					!point defectCurrent after and defectPrev before where it should go in defectList.
	
					call findDefectInList(defectCurrent, defectPrev, products)
	
					!Next update defects in myBoundary
					if(associated(defectCurrent)) then !if we aren't at the end of the list
						same=0
						do 110 l=1,numSpecies
							if(defectCurrent%defectType(l)==products(l)) then
								same=same+1
							endif
						110 continue
						if(same==numSpecies) then	!if the defect is already present in the list
							defectCurrent%num=defectCurrent%num+finalBufferRecv(j,numSpecies+2) !This should be +1 only
							
							if(defectCurrent%num==0) then
								write(*,*) 'Error zero defects in updateDefectList, step 5 of communication'
							endif
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=finalBufferRecv(j,numSpecies+1)
							defectPrev%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
							do 111 l=1,numSpecies
								defectPrev%defectType(l)=finalBufferRecv(j,l)
							111 continue
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						endif
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=finalBufferRecv(j,numSpecies+1)
						defectPrev%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
						do 112 l=1,numSpecies
							defectPrev%defectType(l)=finalBufferRecv(j,l)
						112 continue
					endif
				
				endif
			108 continue
			
			deallocate(finalBufferRecv)
		endif
	endif

107 continue

!NOTE: steps 4-5 are necessary in the case that one element is a member of the boundary of 
!more than one processor (ex: the corner of a processor's mesh).

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
use mod_srscd_constants
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
integer numUpdateLocal(6)			!<the number of defects being sent to each processor neighbor
integer numUpdateBndry(6)			!<the number of defects being recieved from each processor neighbor

integer count, numDefectsRecv
integer numUpdateLocalRecv(6), numUpdateBndryRecv(6)	!number of defects being recieved from each proc neighbor

!NOTE: this final step could be eliminated by keeping the global mesh in each local processor
!(thus each element could be identified as being part of the local mesh of one proc and the boundary of 
!any other procs)
integer numUpdateFinal(6), numUpdateFinalRecv(6)		!number of defects being sent/recieved in final update step

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
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
	
	subroutine chooseCascade(CascadeTemp)
	use DerivedType
	type(cascadeEvent), pointer :: CascadeTemp
	end subroutine
	
	subroutine initializeFineMesh(CascadeCurrent)
	use DerivedType
	type(cascade), pointer :: CascadeCurrent
	end subroutine
end interface

!initialize the trackers of number of defects to update
do i=1,6
	numUpdateLocal(i)=0
	numUpdateBndry(i)=0
end do

!create buffers: size greater than max size needed (max 3 reactants/products per reaction, numCells reactions)
allocate(localBuffer(6,numCells*3,numSpecies+3))
allocate(bndryBuffer(6,numCells*3,numSpecies+2))

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
				defectUpdateCurrent%defectType(j)=reactionCurrent%reactants(i,j)
			end do
			defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i)
			defectUpdateCurrent%neighbor=0	!not pointed at a different proc
			defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
			nullify(defectUpdateCurrent%next)
			
			nullify(defectPrev)
			defectCurrent=>DefectList(reactionCurrent%cellNumber(i))
			
			call findDefectInList(defectCurrent, defectPrev, reactionCurrent%reactants(i,:))
			
			!Check that defectCurrent is pointing towards the reactant
			same=0
			do j=1,numSpecies
				if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
					same=same+1
				endif
			end do
			
			if(same /= numSpecies) then
				write(*,*) 'Proc', myProc%taskid, 'error defectCurrent not pointing to reactants'
				write(*,*) 'reactants', reactionCurrent%reactants(:,:)
				write(*,*) 'products', reactionCurrent%products(:,:)
				write(*,*) 'searching for', reactionCurrent%reactants(i,:)
				if(associated(defectCurrent)) write(*,*) 'DefectCurrent', defectCurrent%defectType
				if(associated(defectPrev)) write(*,*) 'DefectPrev', defectPrev%defectType
			endif
			
!			do 11 while(associated(defectCurrent)) !search for reactant defect
!				same=0
!				do 77 j=1,numSpecies
!					if(defectCurrent%defectType(j)==reactionCurrent%reactants(i,j)) then
!						same=same+1
!					endif
!				77 continue
!				if(same==numSpecies) then
!					exit
!				endif
!				defectPrev=>defectCurrent
!				defectCurrent=>defectCurrent%next
!			11 continue
			
			!First update local buffer if needed
			do j=1, 6
				do k=1,myMesh(reactionCurrent%cellNumber(i))%numNeighbors(j)
					if(myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k) /= myProc%taskid .AND. &
						myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k) /= -1) then	!neighboring element not in this proc
						
						numUpdateLocal(j)=numUpdateLocal(j)+1
						
						do l=1,numSpecies
							localBuffer(j,numUpdateLocal(j),l)=reactionCurrent%reactants(i,l)	!defect type
						end do
						
						!Cell Number in local mesh
						localBuffer(j,numUpdateLocal(j),numSpecies+1)=reactionCurrent%cellNumber(i)
						
						!Cell number in boundary mesh
						localBuffer(j,numUpdateLocal(j),numSpecies+3)=myMesh(reactionCurrent%cellNumber(i))%neighbors(j,k)
						
						if(associated(defectCurrent)) then
							localBuffer(j,numUpdateLocal(j),numSpecies+2)=-1
						else
							write(*,*) 'Error tried to delete defect that was not there and send to neighboring proc'
						endif
						
						!write(*,*) 'localBuffer update: dir', j, 'proc', myProc%taskid, 'neighbor', &
						!	myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k)
						!write(*,*) 'buffer:', (localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2)
						
						!write(*,*) 'proc', myProc%taskid, 'neighborProc', myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j,k)
						!write(*,*) (localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+1)
	!					if(myProc%taskid==4) then
	!						if(j==5) then
	!							write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'remove'
	!						else if(j==6) then
	!							write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'remove'
	!						endif
	!					endif
					endif
				end do
			end do
			
			!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
			if(associated(defectCurrent) .EQV. .FALSE.) then
				write(*,*) 'Tried to delete defect that wasnt there coarse'
				write(*,*) 'reactants', reactionCurrent%reactants
				write(*,*) 'products', reactionCurrent%products
				write(*,*) 'cells', reactionCurrent%CellNumber
				write(*,*) 'rate', reactionCurrent%reactionRate
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
			
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
				DefectList(reactionCurrent%cellNumber(i))%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
				defectUpdateCurrent%num=-1
			
			!if there is one defect of this type and it is the only defect in the list, then make its number 0 (don't remove first defect from list)
			else if(defectCurrent%num==1) then 	!removing only defect from cell i (single helium) - this is redundant but will keep for now
				DefectList(reactionCurrent%cellNumber(i))%num=0
				defectUpdateCurrent%num=-1
			
			!if the defect is in the list but none present, we have chosen a reaction that shouldn't exist
			else if(defectCurrent%num==0) then
				write(*,*) 'trying to remove a defect that isnt there'
				call MPI_ABORT(MPI_COMM_WORLD,ierr)
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
						products(j)=reactionCurrent%products(i,j)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(i,j)
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
								defectPrev%defectType(j)=reactionCurrent%products(i,j)
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
							defectPrev%defectType(j)=reactionCurrent%products(i,j)
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
								if(myMesh(reactionCurrent%cellNumber(1))%neighbors(j,k)==&
									reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
									myMesh(reactionCurrent%cellNumber(1))%neighborProcs(j,k)==&
									reactionCurrent%taskid(i+reactionCurrent%numReactants)) then
			
									numUpdateBndry(j)=numUpdateBndry(j)+1
									do l=1,numSpecies
										bndryBuffer(j,numUpdateBndry(j),l)=reactionCurrent%products(i,l)
									end do
									bndryBuffer(j,numUpdateBndry(j),numSpecies+1)=&
										reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
									bndryBuffer(j,numUpdateBndry(j),numSpecies+2)=1	!indicates we are adding one defect
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
						if(myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%localNeighbor /= 0 &
							.AND. myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%proc == &
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
							defectUpdateCurrent%neighbor=myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))&
								%localNeighbor		
							
							do l=1,numSpecies
								products(l)=reactionCurrent%products(i,l)
								defectUpdateCurrent%defectType(l)=reactionCurrent%products(i,l)
							end do
							
							!Find defect in defect list on myBoundary
							defectCurrent=>myBoundary(j,reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%defectList
	
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
										defectPrev%defectType(l)=reactionCurrent%products(i,l)
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
									defectPrev%defectType(l)=reactionCurrent%products(i,l)
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
						products(j)=reactionCurrent%products(i,j)
						defectUpdateCurrent%defectType(j)=reactionCurrent%products(i,j)
					end do
					
					defectCurrent=>DefectList(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))
	
					!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
					! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
	
					call findDefectInList(defectCurrent, defectPrev, products)
	
					!First update local buffer if needed
					do j=1,6
						do k=1,myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%numNeighbors(j)
							if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(j,k)==-1) then
								!do nothing, free surface
							else if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighborProcs(j,k) &
								/= myProc%taskid) then	!neighboring element not in this proc
							
								numUpdateLocal(j)=numUpdateLocal(j)+1
								
								!species
								do l=1,numSpecies
									localBuffer(j,numUpdateLocal(j),l)=reactionCurrent%products(i,l)
								end do
								
								!cell number in local mesh
								localBuffer(j,numUpdateLocal(j),numSpecies+1)=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
								
								!cell number of neighbor (in different proc)
								localBuffer(j,numUpdateLocal(j),numSpecies+3)=myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))&
									%neighbors(j,k)
								
								!number of defects to be increased by 1
								localBuffer(j,numUpdateLocal(j),numSpecies+2)=1
								
	!							if(myProc%taskid==4) then
	!								if(j==5) then
	!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
	!								else if(j==6) then
	!									write(*,*) 'dir', j, 'buffer contents',(localBuffer(j,numUpdateLocal(j),l),l=1,numSpecies+2), 'add'
	!								endif
	!							endif
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
								defectPrev%defectType(j)=reactionCurrent%products(i,j)
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
							defectPrev%defectType(j)=reactionCurrent%products(i,j)
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

		call MPI_SEND(numUpdateLocal(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+i,MPI_COMM_WORLD, ierr)	!number of local defects that have changed on boundary of processor i
		call MPI_SEND(numUpdateBndry(i),1,MPI_INTEGER,myProc%procNeighbor(i),i+6,MPI_COMM_WORLD, ierr)		!number of defects in the mesh of processor i that have changed (diffusion only)
		
		!EDIT: need to only send the first numUpdateLocal(i) elements of localBuffer(i,:,:)
		
		if(numUpdateLocal(i) /= 0) then

			allocate(localBufferSend(numUpdateLocal(i),numSpecies+3))
			
			do j=1,numUpdateLocal(i)
				do k=1,numSpecies+3
					localBufferSend(j,k)=localBuffer(i,j,k)
				end do
			end do

			call MPI_SEND(localBufferSend,numUpdateLocal(i)*(numSpecies+3), MPI_INTEGER, &
				myProc%procNeighbor(i),i+12,MPI_COMM_WORLD,ierr)
!			if(myProc%taskid==4) then
!				write(*,*) 'dir', i, 'proc sent', myProc%procNeighbor(i), 'numUpdateLocal', numUpdateLocal(i)
!			endif

			deallocate(localBufferSend)

		end if
		
		if(numUpdateBndry(i) /= 0) then

			allocate(bndryBufferSend(numUpdateBndry(i),numSpecies+2))
			
			do j=1,numUpdateBndry(i)
				do k=1,numSpecies+2
					bndryBufferSend(j,k)=bndryBuffer(i,j,k)
				end do
			end do

			call MPI_SEND(bndryBufferSend,numUpdateBndry(i)*(numSpecies+2), MPI_INTEGER, &
				myProc%procNeighbor(i),i+18,MPI_COMM_WORLD,ierr)
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

		call MPI_RECV(numUpdateBndryRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),200+tag,MPI_COMM_WORLD,status,ierr)	!number of bndry defects that have changed (local to processor i, boundary here)
		call MPI_RECV(numUpdateLocalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+6,MPI_COMM_WORLD,status,ierr)	!number of local defects that have changed (bndry to processor i, local here)
		totalLocalRecv=totalLocalRecv+numUpdateLocalRecv(i)
	endif
end do

allocate(finalBuffer(6,totalLocalRecv,numSpecies+3))

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
			allocate(bndryBufferRecv(numUpdateBndryRecv(i),numSpecies+3))
			call MPI_RECV(bndryBufferRecv,numUpdateBndryRecv(i)*(numSpecies+3),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+12,MPI_COMM_WORLD,status,ierr)
!			if(myProc%taskid==MASTER) then
!				write(*,*) 'dir', i, 'proc recvd', myProc%procNeighbor(i), 'numUpdateBndryRecv', numUpdateBndryRecv(i)
!			endif
			
			!Add defects in bndryBufferRecv to myBoundary()
			do j=1,numUpdateBndryRecv(i)
!				if(myProc%taskid==0) then
!					if(i==5) then
!						write(*,*) 'dir', i, 'bdry recv',(bndryBufferRecv(j,m),m=1,numSpecies+2)
!					else if(i==6) then
!						write(*,*) 'dir', i, 'bdry recv',(bndryBufferRecv(j,m),m=1,numSpecies+2)
!					endif
!				endif
				
				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%procNeighbor(i)
				defectUpdateCurrent%dir=i	
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=bndryBufferRecv(j,numSpecies+1)
				defectUpdateCurrent%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)
				
				!This information is used to tell us which element in the local
				!processor bounds this element in myBoundary
				defectUpdateCurrent%neighbor=bndryBufferRecv(j,numSpecies+3)
				
				!point defectCurrent at the defect list in the correct cell of myBoundary
				defectCurrent=>myBoundary(i,bndryBufferRecv(j,numSpecies+1))%defectList
				if(.NOT. associated(defectCurrent)) then
					write(*,*) 'error myBoundary not allocated correctly'
					write(*,*) 'dir', i, 'cell', bndryBufferRecv(j,numSpecies+1)
					call MPI_ABORT(MPI_COMM_WORLD,ierr)
				endif
				
				nullify(defectPrev)
				
				do k=1,numSpecies
					products(k)=bndryBufferRecv(j,k)
					defectUpdateCurrent%defectType(k)=bndryBufferRecv(j,k)	!used to update reaction lists
				end do
				
				!point DefectCurrent at the defect we are looking for (if it is there), otherwise
				!point defectCurrent after and defectPrev before where it should go in defectList.
!				if(myProc%taskid==MASTER .OR. myProc%taskid==3) then
!					write(*,*) 'proc', myProc%taskid
!					write(*,*) 'numUpdateBndryRecv', numUpdateBndryRecv(i)
!					write(*,*) 'dir', i, (bndryBufferRecv(j,k),k=1,numSpecies+2)
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
						defectCurrent%num=defectCurrent%num+bndryBufferRecv(j,numSpecies+2) !this will be +/- 1 only
						
						if(defectCurrent%num==0) then
							!delete this defect from the list in myBoundary
							defectPrev%next=>defectCurrent%next !remove that defect type from the system
							deallocate(defectCurrent%defectType)
							deallocate(defectCurrent)
							nullify(defectCurrent)
						end if
					else		!if the defect is to be inserted in the list
					
						if(bndryBufferRecv(j,numSpecies+2) == -1) then
							write(*,*) 'defectCurrent associated', defectCurrent%defectType, 'num', defectCurrent%num
							if(associated(defectPrev)) write(*,*) 'defectPrev', defectPrev%defectType, 'num', defectPrev%num
							write(*,*) 'error in defectUpdate negative defect numbers', myProc%taskid
							write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
							write(*,*) (bndryBufferRecv(j,k),k=1,numSpecies+2)
							
							call MPI_ABORT(MPI_COMM_WORLD,ierr)
							
						else
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=bndryBufferRecv(j,numSpecies+1)
							defectPrev%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
							
							do l=1,numSpecies
								defectPrev%defectType(l)=bndryBufferRecv(j,l)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
												
					end if
				else 			!add a defect to the end of the list
					
					if(bndryBufferRecv(j,numSpecies+2) == -1) then
						write(*,*) 'defectCurrent not associated'
						write(*,*) 'error in defectUpdate negative defect numbers'
						write(*,*) 'proc', myProc%taskid, 'dir', i, 'neighbor proc', myProc%procNeighbor(i)
						write(*,*) (bndryBufferRecv(j,k),k=1,numSpecies+2)
						call MPI_ABORT(MPI_COMM_WORLD,ierr)
					else
					
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=bndryBufferRecv(j,numSpecies+1)
						defectPrev%num=bndryBufferRecv(j,numSpecies+2)	!This will be +/- 1 only
											
						do l=1,numSpecies
							defectPrev%defectType(l)=bndryBufferRecv(j,l)
						end do
					end if

				end if
			end do
			deallocate(bndryBufferRecv)
		end if
		
		if(numUpdateLocalRecv(i) /= 0) then
			
			!Read in defects to update in local mesh
			allocate(localBufferRecv(numUpdateLocalRecv(i),numSpecies+2))
			call MPI_RECV(localBufferRecv,numUpdateLocalRecv(i)*(numSpecies+2),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+18,MPI_COMM_WORLD,status,ierr)

			!Add defects in localBufferRecv to defectList()
			do j=1,numUpdateLocalRecv(i)
				!create a new element in defectUpdate and assign all variables except for num (will do later)
				allocate(defectUpdateCurrent%next)
				defectUpdateCurrent=>defectUpdateCurrent%next
				defectUpdateCurrent%proc=myProc%taskid
				defectUpdateCurrent%dir=0	!not pointed at a different proc
				defectUpdateCurrent%neighbor=0	!not pointed at a different proc
				allocate(defectUpdateCurrent%defectType(numSpecies))
				defectUpdateCurrent%cellNumber=localBufferRecv(j,numSpecies+1)
				defectUpdateCurrent%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only
				defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
				nullify(defectUpdateCurrent%next)

				!point defectCurrent at the defect list in the correct cell of defectList
				defectCurrent=>defectList(localBufferRecv(j,numSpecies+1))
				nullify(defectPrev)
				
				do k=1,numSpecies
					products(k)=localBufferRecv(j,k)
					defectUpdateCurrent%defectType(k)=localBufferRecv(j,k)
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
						defectCurrent%num=defectCurrent%num+localBufferRecv(j,numSpecies+2)	!This will be +/- 1 only (should be only +1 here)
						
						if(defectCurrent%num==0) then
							write(*,*) 'Error zero defects in updateDefectList, step 3 of communication'
						endif
					else		!if the defect is to be inserted in the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=localBufferRecv(j,numSpecies+1)
						defectPrev%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only (should be only +1 here)
						do l=1,numSpecies
							defectPrev%defectType(l)=localBufferRecv(j,l)
						end do
						defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
					endif
				else 			!add a defect to the end of the list
					nullify(defectPrev%next)
					allocate(defectPrev%next)
					nullify(defectPrev%next%next)
					defectPrev=>defectPrev%next
					allocate(defectPrev%defectType(numSpecies))
					defectPrev%cellNumber=localBufferRecv(j,numSpecies+1)
					defectPrev%num=localBufferRecv(j,numSpecies+2) !This will be +/- 1 only (should be only +1 here)
					do l=1,numSpecies
						defectPrev%defectType(l)=localBufferRecv(j,l)
					end do
				endif
				
				!*****************************************************
				!Prep for step 4: check to see if element updated has more than one neighbor proc
				!*****************************************************
				
				do k=1,6
					do l=1,myMesh(localBufferRecv(j,numSpecies+1))%numNeighbors(k)
						
						!if a neighbor of this element is not in this proc or in the proc that just communicated with it
						if (myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) /= myProc%taskid .AND. &
							myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) /= myProc%procNeighbor(i) .AND. &
							myMesh(localBufferRecv(j,numSpecies+1))%neighborProcs(k,l) /= -1) then
							
							!Add this defect to a final buffer
							numUpdateFinal(k)=numUpdateFinal(k)+1
							!finalBuffer(1,1,1)=100	!why is this here?
							do m=1,numSpecies
								finalBuffer(k,numUpdateFinal(k),m)=localBufferRecv(j,m)
							end do
							finalBuffer(k,numUpdateFinal(k),numSpecies+1)=localBufferRecv(j,numSpecies+1)	!local cell number
							finalBuffer(k,numUpdateFinal(k),numSpecies+2)=localBufferRecv(j,numSpecies+2) !This should be +1 only
							finalBuffer(k,numUpdateFinal(k),numSpecies+3)=myMesh(localBufferRecv(j,numSpecies+1))%neighbors(k,l)	!cell number in neighboring proc
							
							if(localBufferRecv(j,numSpecies+2)==-1) then
								write(*,*) 'error: cell in boundary of multiple procs removing defect'
								call MPI_ABORT(MPI_COMM_WORLD,ierr)
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
		call MPI_SEND(numUpdateFinal(i), 1, MPI_INTEGER, myProc%procNeighbor(i), i+24, MPI_COMM_WORLD, ierr)	!number of local defects that have changed on boundary of processor i
		
		if(numUpdateFinal(i) /= 0) then
			allocate(finalBufferSend(numUpdateFinal(i),numSpecies+3))
			
			do j=1,numUpdateFinal(i)
				do k=1,numSpecies+3
					finalBufferSend(j,k)=finalBuffer(i,j,k)
				end do
			end do
			
			call MPI_SEND(finalBufferSend,numUpdateFinal(i)*(numSpecies+3), MPI_INTEGER, &
				myProc%procNeighbor(i),i+30,MPI_COMM_WORLD,ierr)
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
		call MPI_RECV(numUpdateFinalRecv(i),1,MPI_INTEGER,myProc%procNeighbor(i),tag+24,MPI_COMM_WORLD,status,ierr)	
		
		if(numUpdateFinalRecv(i) /= 0) then
			
			!Read in defects to update in boundary
			allocate(finalBufferRecv(numUpdateFinalRecv(i),numSpecies+3))
			call MPI_RECV(finalBufferRecv,numUpdateFinalRecv(i)*(numSpecies+3),MPI_INTEGER,&
				myProc%procNeighbor(i),tag+30,MPI_COMM_WORLD,status,ierr)

			!Add defects in finalBufferRecv to myBoundary()
			do j=1,numUpdateFinalRecv(i)
				
!				if(myProc%taskid==0) then
!					if(i==5) then
!						write(*,*) 'dir', i, 'final recv',(finalBufferRecv(j,m),m=1,numSpecies+2)
!					else if(i==6) then
!						write(*,*) 'dir', i, 'final recv',(finalBufferRecv(j,m),m=1,numSpecies+2)
!					endif
!				endif
				
				if(finalBufferRecv(j,numSpecies+2)==-1) then
					write(*,*) 'error deleting defects in finalBufferRecv'
					!call MPI_ABORT(MPI_COMM_WORLD,ierr)
				else
					
					!create a new element in defectUpdate and assign all variables except for num (will do later)
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					defectUpdateCurrent%proc=myProc%procNeighbor(i)
					defectUpdateCurrent%dir=i	
					allocate(defectUpdateCurrent%defectType(numSpecies))
					defectUpdateCurrent%cellNumber=finalBufferRecv(j,numSpecies+1)	!cell number in boundary mesh
					defectUpdateCurrent%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
					defectUpdateCurrent%neighbor=finalBufferRecv(j,numSpecies+3) !cell number in myMesh
					defectUpdateCurrent%cascadeNumber=0	!not inside a cascade
					nullify(defectUpdateCurrent%next)
					
					!point defectCurrent at the defect list in the correct cell of myBoundary
					defectCurrent=>myBoundary(i,finalBufferRecv(j,numSpecies+1))%defectList
					nullify(defectPrev)
					
					do k=1,numSpecies
						products(k)=finalBufferRecv(j,k)
						defectUpdateCurrent%defectType(k)=finalBufferRecv(j,k)
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
							defectCurrent%num=defectCurrent%num+finalBufferRecv(j,numSpecies+2) !This should be +1 only
							
							if(defectCurrent%num==0) then
								write(*,*) 'Error zero defects in updateDefectList, step 5 of communication'
							end if
						else		!if the defect is to be inserted in the list
							nullify(defectPrev%next)
							allocate(defectPrev%next)
							nullify(defectPrev%next%next)
							defectPrev=>defectPrev%next
							allocate(defectPrev%defectType(numSpecies))
							defectPrev%cellNumber=finalBufferRecv(j,numSpecies+1)
							defectPrev%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
							do l=1,numSpecies
								defectPrev%defectType(l)=finalBufferRecv(j,l)
							end do
							defectPrev%next=>defectCurrent !if inserted defect is in the middle of the list, point it to the next item in the list
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrev%next)
						allocate(defectPrev%next)
						nullify(defectPrev%next%next)
						defectPrev=>defectPrev%next
						allocate(defectPrev%defectType(numSpecies))
						defectPrev%cellNumber=finalBufferRecv(j,numSpecies+1)
						defectPrev%num=finalBufferRecv(j,numSpecies+2) !This should be +1 only
						do l=1,numSpecies
							defectPrev%defectType(l)=finalBufferRecv(j,l)
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
use mod_srscd_constants
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

defectUpdateCurrent=>defectUpdate%next

!loop through defects that needs reactions updated in defectUpdate
do 10 while(associated(defectUpdateCurrent))
	
	if(defectUpdateCurrent%num .NE. 1 .AND. defectUpdateCurrent%num .NE. -1) then
		!we have an error; all defectUpdateCurrent members should have num = +/-1 (indicates adding or removing)
		write(*,*) 'error defectUpdateCurrent%num not equal to +/- 1', myProc%taskid, defectUpdateCurrent%num
	endif
	
	do 11 i=1,numSpecies
		defectTemp(i)=defectUpdateCurrent%defectType(i)
	11 continue
			
	!if the defect is within the local mesh, update all relevant reactions
	
	if(defectUpdateCurrent%proc==myProc%taskid) then
		
		!*******************************************************************************************
		!
		! defectUpdateCurrent%cascadeNumber==0 means a defect has changed in the coarse mesh
		!
		! if the defect is within the coarse mesh, then update all relevant reactions in the coarse mesh
		!
		!*******************************************************************************************

		if(defectUpdateCurrent%cascadeNumber==0) then
			
			!Single-defect reactions associated with defects of type defectTemp	
			call addSingleDefectReactions(defectUpdateCurrent%cellNumber,defectTemp)
			
			!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
			!all defect types in the defect list)
			defectCurrent=>defectList(defectUpdateCurrent%cellNumber)
			
			do 59 while(associated(defectCurrent))
				if(defectCurrent%num .NE. 0) then
					call addMultiDefectReactions(defectUpdateCurrent%cellNumber,defectTemp, defectCurrent%defectType)
				endif
				defectCurrent=>defectCurrent%next
			59 continue
			
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
				endif
			endif
			
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
			do 60 j=1,6
				if (myMesh(defectUpdateCurrent%cellNumber)%numNeighbors(j)==0) then
					write(*,*) 'error myMesh does not have neighbors in this direction'
				endif
				do 61 k=1,myMesh(defectUpdateCurrent%cellNumber)%numNeighbors(j)
	
	!				if(myProc%taskid==MASTER) then
	!					write(*,*) 'dir', j, 'num1', findNumDefect(defectTemp, defectUpdateCurrent%cellNumber)
	!					if(myProc%taskid==myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k)) then
	!						write(*,*) 'num2_local', findNumDefect(defectTemp,myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k))
	!					else
	!						write(*,*) 'num2_bndry', findNumDefectBoundary(defectTemp,defectUpdateCurrent%cellNumber,j)
	!					endif
	!				endif
					
					!Add diffusion reactions from this cell to neighboring cells
					if(polycrystal=='yes') then

						!Find the grain ID number of the volume element we are in
						localGrainID=myMesh(defectUpdateCurrent%cellNumber)%material

						!Find the grain ID number of the neighboring volume element
						if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k) .NE. myProc%taskid .AND. &
							myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k) .NE. -1) then
						
							neighborGrainID=myBoundary(j,myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k))%material
						else if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k) == -1) then
							neighborGrainID=localGrainID	!free surface release, don't need to do anything special
						else
							neighborGrainID=myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k))%material
						endif

						if(localGrainID==neighborGrainID) then

							!Allow diffusion between elements in the same grain
							call addDiffusionReactions(defectUpdateCurrent%cellNumber, myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k),&
								myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k),j,defectTemp)
						
						else

							!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
							call addDiffusionReactions(defectUpdateCurrent%cellNumber, 0, myProc%taskid, -1, j, defectTemp)													
							
						endif

					else
						
						call addDiffusionReactions(defectUpdateCurrent%cellNumber, myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k),&
							myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k),j,defectTemp)
							
					endif
						
					!If the neighboring volume element is in the same processor as this element, add
					!diffusion reactions from neighboring cells into this cell (changes reaction list
					!in neighboring cell)
					!
					!NOTE: don't need to switch j (direction) here even though diffusion is in opposite
					!direction because dir is only used when diffusing into/out of neighbor processor
					
					if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k)==myProc%taskid) then	
					!Don't need to do this for diffusion between different grains
					if(polycrystal=='yes' .AND. myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k))%material &
						== myMesh(defectUpdateCurrent%cellNumber)%material) then
							
						call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k), defectUpdateCurrent%cellNumber, &
							myProc%taskid, myProc%taskid, j, defectTemp)
						
					else if(polycrystal=='no') then	!2015.05.20 Noticed that this wasn't here, probably was a bug. Fixed?
					
						call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k), defectUpdateCurrent%cellNumber, &
							myProc%taskid, myProc%taskid, j, defectTemp)
							
					endif				
					endif

				61 continue
				
			60 continue
			
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
			
			do 62 while(associated(CascadeCurrent))
			
				if(CascadeCurrent%cellNumber==defectUpdateCurrent%cellNumber) then
					
!					if(myProc%taskid==MASTER) then
!						write(*,*) 'adding coarse to fine diffusion'
!						read(*,*)
!					endif
					
					call addDiffusionCoarseToFine(defectUpdateCurrent%cellNumber, myProc%taskid, CascadeCurrent, defectTemp)
				
				endif
				
				CascadeCurrent=>CascadeCurrent%next
				
			62 continue

		!*******************************************************************************************
		!
		! DefectUpdateCurrent%cascadeNumber .NE. 0 means that the defect that has been updated is inside a cascade.
		! 
		! Therefore we update all relevant reactions within the correct cascade's fine mesh
		!
		!*******************************************************************************************
			
		else

			!Point CascadeCurrent at the cascade associated with the chosen reaction

			CascadeCurrent=>ActiveCascades
			
			do 161 while(associated(CascadeCurrent))
				
				if(CascadeCurrent%cascadeID==defectUpdateCurrent%cascadeNumber) then
					exit
				endif
			
				CascadeCurrent=>CascadeCurrent%next
			161 continue

			!Check to make sure that we have exited at the correct cascade
			if(.NOT. associated(CascadeCurrent)) then
				write(*,*) 'error CascadeCurrent not associated in updateReactionList'
			endif

			!Single-defect reactions associated with defects of type defectTemp in fine mesh	
			call addSingleDefectReactionsFine(defectUpdateCurrent%cascadeNumber,defectUpdateCurrent%cellNumber,defectTemp)

			!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
			!all defect types in the defect list)
			defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
			do 159 while(associated(defectCurrent))
				if(defectCurrent%num .NE. 0) then
				
					call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber,&
						defectTemp, defectCurrent%defectType)
				
				endif
					
				defectCurrent=>defectCurrent%next
				
			159 continue

			!See above for explanation of this section
			defectUpdateNext=>defectUpdateCurrent%next
			if(associated(defectUpdateNext)) then
				if(defectUpdateNext%cellNumber==defectUpdateCurrent%cellNumber) then
					
					call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber, &
						defectTemp, defectUpdateNext%defectType)
						
				endif
			endif
			
			!Diffusion reactions
			do 160 j=1,6
	
	!				if(myProc%taskid==MASTER) then
	!					write(*,*) 'dir', j, 'num1', findNumDefect(defectTemp, defectUpdateCurrent%cellNumber)
	!					if(myProc%taskid==myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j,k)) then
	!						write(*,*) 'num2_local', findNumDefect(defectTemp,myMesh(defectUpdateCurrent%cellNumber)%neighbors(j,k))
	!					else
	!						write(*,*) 'num2_bndry', findNumDefectBoundary(defectTemp,defectUpdateCurrent%cellNumber,j)
	!					endif
	!				endif
				
				!Add diffusion reactions from this cell into neighboring cells as well as diffusion reactions
				!from neighboring cells into this cell. (both directions)

				call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber, defectUpdateCurrent%cellNumber, &
					cascadeConnectivity(defectUpdateCurrent%cellNumber,j),myProc%taskid, myProc%taskid,j,defectTemp)
				
				if(cascadeConnectivity(defectUpdateCurrent%cellNumber,j) .NE. 0) then
					call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber, cascadeConnectivity(defectUpdateCurrent%cellNumber,j), &
						defectUpdateCurrent%cellNumber,myProc%taskid, myProc%taskid,j,defectTemp)
				endif

			160 continue
			
		endif

	!if the defect is within the boundary mesh, only update the diffusion reactions for cells 
	!touching that boundary element
	else

		!If we are in this section, this means that a defect has changed in the boundary to this 
		!volume element but not IN the volume element. Therefore we only need to update one diffusion
		!reaction in one direction, the direction of the element with the changed defect. No defects
		!have changed in this volume element.
		
		!*******************************************************************************************
		!EDIT: 06/05/2014 Updated addDiffusionReactions search in the boundary by adding a %neighbor
		!item in defectUpdate so that we don't have to search for the local element that is the
		!neighbor of the boundary element.
		!*******************************************************************************************
		
		!*******************************************************************************************
		!2015.04.06 Possible bug: is defectUpdateCurrent%dir pointing in the backwards direction?
		!(from neighboring cell TO local cell instead of vice versa) This may be the case because
		!defectUpdateCurrent was created in a different processor and passed here. Need to check
		!that it was initialized correctly.
		!*******************************************************************************************
		
		if(defectUpdateCurrent%neighbor==-1) then
			write(*,*) 'error neighbor not assigned for diffusion reactions into boundary'
			call MPI_ABORT(MPI_COMM_WORLD,ierr)
		endif
		
		if(polycrystal=='yes') then
		
			localGrainID=myMesh(defectUpdateCurrent%neighbor)%material
			
			neighborGrainID=myBoundary(defectUpdateCurrent%dir,defectUpdateCurrent%cellNumber)%material
			
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
			
		endif
		
	endif
	defectUpdateCurrent=>defectUpdateCurrent%next

10 continue

!Deallocate defectUpdate (memory release)
defectUpdatePrev=>defectUpdate
defectUpdateCurrent=>defectUpdate%next
do 20 while(associated(defectUpdatePrev))
	deallocate(defectUpdatePrev%defectType)
	deallocate(defectUpdatePrev)
	defectUpdatePrev=>defectUpdateCurrent
	if(associated(defectUpdateCurrent)) then
		defectUpdateCurrent=>defectUpdateCurrent%next
	endif
20 continue

end subroutine
