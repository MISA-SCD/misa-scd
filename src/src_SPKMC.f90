!*****************************************************************************************
!>double precision generate timestep - chooses a timestep using random number and Monte Carlo algorithm
!(this is a global timestep)
!*****************************************************************************************
double precision function GenerateTimestep()
	use mod_structures
	use mod_globalVariables
	use mod_randdp
	implicit none

	double precision :: r1

	r1=dprand()
	GenerateTimestep=dlog(1d0/r1)/maxRate

end function

!*****************************************************************************************
!>subroutine
!*****************************************************************************************
subroutine chooseImplantReaction(reactionCurrent)
	use mod_structures
	use mod_globalVariables
	implicit none

	type(reaction), pointer, intent(inout) :: reactionCurrent
	integer :: i

	write(*,*) 'step', step, 'oneCascadeGCell',oneCascadeGCell
	nullify(reactionCurrent)	!These are default pointed at nothing, indicating null event

	!***********************************************************************
	!Choose from reactions within the coarse mesh
	!***********************************************************************
	outer: do i=1,numCells

		!search for which volume element we are choosing among
		if(oneCascadeGCell == myMesh(i)%globalCell) then
			!a reaction is chosen in this volume element, so we no longer nullify reactionCurrent
			reactionCurrent=>reactionList(i)
			exit
		end if
	end do outer

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
	end if
end subroutine

!*****************************************************************************************
!>subroutine chooseReaction(reactionCurrent, CascadeCurrent)
!chooses a reaction in each processor according to the Monte Carlo algorithm (this is a local reaction)
!*****************************************************************************************
subroutine chooseReaction(reactionCurrent, CascadeCurrent)
	use mod_structures
	use mod_globalVariables
	use mod_randdp
	implicit none

	type(reaction), pointer, intent(inout) :: reactionCurrent
	type(cascade), pointer, intent(inout) :: cascadeCurrent
	type(reaction), pointer :: reactionTemp
	type(defect), pointer :: defectTemp
	double precision :: r2, atemp, atemp_cell, r2timesa, atemp_test
	integer :: i

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

			reactionCurrent=>reactionList(i)
			!search for which reaction occurs in that volume element
			do while(associated(reactionCurrent))
				atemp=atemp+reactionCurrent%reactionRate
				if(r2timesa <= atemp) then
					exit outer	!exit both loops with reactionCurrent pointing to the randomly chosen reaction
				end if
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
						end if
						reactionCurrent=>reactionCurrent%next
					end do
				else
					atemp=atemp+CascadeCurrent%totalRate(i)	!skip this volume element
				end if
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
!> Subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
!This is the most involved and complex subroutine of the entire algorithm.
!This subroutine updates defects in the local mesh according to the reaction chosen, and
!communicates with neighboring processors about defects that may have passed into
!a different processor as well as about defects that may have changed on the boundary.
!It also creates a list of defects that have been updated, to inform the next subroutine (updateReactionList(defectUpdate))
!which reactions to update.
!***************************************************************************************************
subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
	use mod_globalVariables
	use mod_structures
	use mod_reactionrates
	use mod_randdp
	implicit none
	include 'mpif.h'

	type(reaction), pointer, intent(in) :: reactionCurrent
	type(defectUpdateTracker), pointer, intent(inout) :: defectUpdateCurrent
	type(cascade), pointer, intent(inout) :: CascadeCurrent
	type(defect), pointer :: defectCurrent, defectPrev, defectTemp, defectStoreList, defectStore, defectStorePrev
	type(cascadeDefect), pointer :: cascadeDefectTemp
	type(cascade), pointer :: CascadePrev, CascadeTest
	type(cascadeEvent), pointer :: cascadeTemp
	!Used for cascade
	double precision :: coordinatesTemp(3)
	integer :: cellNumber, mixingEvents, mixingTemp
	logical :: combineBoolean, isCombined
	!Used for cascade recombination
	double precision :: r1, atemp
	integer :: i, j1, j, k, l, m, same, products(numSpecies), product2(numSpecies), totalLocalRecv, count
	double precision :: diffusionRandom
	logical :: flag, flagTemp
	!Used for communication between processors
	integer :: numUpdateSend(6), numUpdateRecv(6)	!the number of defects being sent to (recived from) each processor neighbor
	integer :: numLocalRecv, numBndryRecv			!the number of (loceal/bndry) defects being recieved from each proc neighbor
	!NOTE: this final step could be eliminated by keeping the global mesh in each local processor
	!(thus each element could be identified as being part of the local mesh of one proc and the boundary of any other procs)
	integer :: numUpdateFinal(6), numUpdateFinalRecv(6)	!number of defects being sent/recieved in final update step
	integer :: recvDir
	!create buffers of information to send to neighboring elements (Contains local defects and boundary defects)
	integer, allocatable :: firstSend(:,:,:)
	integer, allocatable :: firstRecv(:,:,:)
	integer, allocatable :: finalSend(:,:,:)	!Only contains local defects
	integer, allocatable :: finalBufferRecv(:,:)
	integer :: status(MPI_STATUS_SIZE)
	integer :: sendFirstStatus(MPI_STATUS_SIZE), recvFirstStatus(MPI_STATUS_SIZE)
	integer :: sendFinalStatus(MPI_STATUS_SIZE), recvFinalStatus(MPI_STATUS_SIZE)
	integer :: sendFirstRequest(6), recvFirstRequest(6)
	integer :: sendFinalRequest(6), recvFinalRequest(6)
	!Function
	integer, external :: findCellWithCoordinatesFineMesh, chooseRandomCell
	logical, external :: cascadeMixingCheck

	!just for testing
	type(defect),pointer :: defectTest
	type(reaction),pointer :: reactionTest
	integer :: countTest

	interface
		subroutine findDefectInList(defectCurrent, defectPrev, products)
			use mod_structures
			use mod_globalVariables
			implicit none
			type(defect), pointer, intent(inout) :: defectCurrent, defectPrev
			integer, intent(in) :: products(numSpecies)
		end subroutine

		subroutine chooseCascade(CascadeTemp)
			use mod_structures
			implicit none
			type(cascadeEvent), pointer, intent(inout) :: CascadeTemp
		end subroutine

		subroutine initializeFineMesh(CascadeCurrent)
			use mod_structures
			type(cascade), pointer, intent(inout) :: CascadeCurrent
		end subroutine
	end interface

	!initialize arraies
	numUpdateSend=0
	numUpdateRecv=0
	numUpdateFinal=0
	numUpdateFinalRecv=0
	!initialize variables
	totalLocalRecv=0

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

	else	!associated(reactionCurrent)	!if we have chosen an event

		!***********************************************************************************************
		!Cascade chosen
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
				call chooseCascade(cascadeTemp)	!cascadeTemp must be associated
				cascadeDefectTemp=>cascadeTemp%ListOfDefects

				!initialize the cascade to be added
				CascadeCurrent=>ActiveCascades
				if(.NOT. associated(CascadeCurrent)) then		!no cascade fine meshes are currently in the system
					allocate(ActiveCascades)
					ActiveCascades%cellNumber=reactionCurrent%cellNumber(1) !cell number in COARSE MESH of this cascade
					ActiveCascades%cascadeID=numImpAnn(1)
					nullify(ActiveCascades%localDefects)
					nullify(ActiveCascades%next)
					nullify(ActiveCascades%prev)
					nullify(ActiveCascades%reactionList)
					allocate(ActiveCascades%totalRate(numCellsCascade))
					do j=1,numCellsCascade
						ActiveCascades%totalRate(j)=0d0							!reaction rate of all reactions within fine mesh
					end do
					CascadeCurrent=>ActiveCascades
					!CascadeCurrent%cascadeID=numImpAnn(1)
				else	!cascade fine meshe are already in the system
					j1=0
					do while(associated(CascadeCurrent))
						CascadePrev=>CascadeCurrent
						CascadeCurrent=>CascadeCurrent%next
						j1=j1+1
					end do
					allocate(CascadeCurrent)
					CascadeCurrent%cellNumber=reactionCurrent%cellNumber(1)	!cell number in COARSE MESH of this cascade
					CascadeCurrent%cascadeID=numImpAnn(1)
					nullify(CascadeCurrent%localDefects)
					nullify(CascadeCurrent%next)
					CascadeCurrent%prev=>CascadePrev
					nullify(CascadeCurrent%reactionList)
					allocate(CascadeCurrent%totalRate(numCellsCascade))
					do j=1,numCellsCascade
						CascadeCurrent%totalRate(j)=0d0
					end do
					CascadePrev%next=>CascadeCurrent
				end if

				!*******************************************************************
				!Each time a cascade is added, a new cascade-fine-mesh is created.
				!*******************************************************************
				!update volume and length of couarse mesh element
				myMesh(cascadeCurrent%cellNumber)%volume = myMesh(cascadeCurrent%cellNumber)%volume &
						-CascadeElementVol*dble(numCellsCascade)
				myMesh(cascadeCurrent%cellNumber)%length = (myMesh(cascadeCurrent%cellNumber)%volume)**(1d0/3d0)

				!If we add too many cascades to a coarse mesh element, the volume inside can become negative. If this
				!happens, output an error message. (Only a danger in the case of very high DPA rates and
				!small coarse mesh elements)
				if(myMesh(cascadeCurrent%cellNumber)%volume <= 0d0) then

					!just for testing
					countTest=0
					CascadeTest=>ActiveCascades%next
					do while(associated(CascadeTest))
						countTest=countTest+1
						CascadeTest=>CascadeTest%next
					end do
					write(*,*) 'Error negative coarse mesh volume    step',step, 'proc', myProc%taskid, &
							'totalCascades',numImpAnn(1), 'No.Cascade', countTest, 'j1', j1-1

					defectTest=>defectList(cascadeCurrent%cellNumber)
					do while(associated(defectTest))
						write(*,*) defectTest%defectType, defectTest%num, defectTest%cellNumber
						defectTest=>defectTest%next
					end do
					write(*,*) '********************************************'
					write(*,*)
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
				nullify(defectStore%next)

				!***************************************************************************************
				!Step 1: Create defectStoreList with defects in cascade only
				!Step 2: For each defect in the fine mesh, check whether it combines with cascade
				!Step 3: If yes, remove it from the fine mesh and combine it randomly with one defect in defectStoreList, then insert the products into  defectStoreList
				!Step 4: Implant defects in defectStoreList in fine mesh. Remember to correctly update defect update protocol.
				!Step 5: Set up defectUpdate for all defects in fine mesh
				!***************************************************************************************
				!Step 1: add defects of cascadeTemp into defectStoreList
				do i=1, cascadeTemp%numDefectsTotal
					!Make sure coordinates of all defects in cascade are within fine mesh, using periodic BC's
					coordinatesTemp(1)=cascadeDefectTemp%coordinates(1)
					if(coordinatesTemp(1) > numxcascade*finelength/2d0) then
						coordinatesTemp(1)=coordinatesTemp(1)-numxcascade*finelength
					else if(coordinatesTemp(1) < -numxcascade*finelength/2d0) then
						coordinatesTemp(1)=coordinatesTemp(1)+numxcascade*finelength
					end if

					coordinatesTemp(2)=cascadeDefectTemp%coordinates(2)
					if(coordinatesTemp(2) > numycascade*finelength/2d0) then
						coordinatesTemp(2)=coordinatesTemp(2)-numycascade*finelength
					else if(coordinatesTemp(2) < -numycascade*finelength/2d0) then
						coordinatesTemp(2)=coordinatesTemp(2)+numycascade*finelength
					end if

					coordinatesTemp(3)=cascadeDefectTemp%coordinates(3)
					if(coordinatesTemp(3) > numzcascade*finelength/2d0) then
						coordinatesTemp(3)=coordinatesTemp(3)-numzcascade*finelength
					else if(coordinatesTemp(3) < -numzcascade*finelength/2d0) then
						coordinatesTemp(3)=coordinatesTemp(3)+numzcascade*finelength
					end if

					!find position cascade is implanted in fine mesh
					cellNumber=findCellWithCoordinatesFineMesh(coordinatesTemp)

					!Create list of defects that need to be added to fine mesh (after casacade-fine mesh
					!mixing has been taken into account)
					allocate(defectStore%next)
					defectStore=>defectStore%next
					allocate(defectStore%defectType(numSpecies))
					do j=1,numSpecies
						defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
					end do
					if(pointDefectToggle=='yes') then
						if(defectStore%defectType(3) > max3DInt) then
							defectStore%defectType(4) = defectStore%defectType(3)
							defectStore%defectType(3) = 0
						end if
					end if
					defectStore%num=1
					defectStore%cellNumber=cellNumber
					nullify(defectStore%next)
					cascadeDefectTemp=>cascadeDefectTemp%next
				end do

				!Step 2: For each initial defect in the fine mesh, check whether it combines with cascade in the defectStoreList
				do j=1,numCellsCascade
					!defectCurrent=>CascadeCurrent%localDefects(j)	!the first defect is 0 0 0 0
					!defectPrev=>defectCurrent
					!defectTemp=>defectCurrent%next

					defectTemp=>CascadeCurrent%localDefects(j)	!the first defect is 0 0 0 0
					defectPrev=>defectTemp
					defectTemp=>defectTemp%next
					do while(associated(defectTemp))

						!Step 2.1: boolean variable, if FALSE then no combination with defectCurrent, if TRUE then combination
						mixingTemp=0
						if(defectTemp%num > 0) then
							do k=1,defectTemp%num
								combineBoolean=cascadeMixingCheck()
								if(combineBoolean .eqv. .TRUE.) then
									mixingTemp=mixingTemp+1
								end if
							end do
						end if

						!Step 2.2: combination
						mixingEvents=0
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
							else if(defectStore%cellNumber==j) then
								!***********************************************************************
								!Hard coded: use defect combination rules to combine defectStore and defecTemp
								!These rules have been transported to a separate subroutine
								!in order to facilitate hard-coding and keep this subroutine clean.
								!***********************************************************************
								do i=1, numSpecies
									product2(i)=0
								end do

								call defectCombinationRules(defectStore%defectType,product2, defectTemp, isCombined)

								if(product2(3)/=0 .OR. product2(4)/=0) then	!product2 is SIAs
									!defectStore is at the middle/end of defectStoreList
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
									else !defectStore is at the beginning of defectStoreList, but not the first
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

				!Step 4: Implant defects in defectStoreList in fine mesh. (add all cascade mixing products and normal cascade products to the system)
				defectStore=>defectStoreList%next
				do while(associated(defectStore))

					!Temporary debug: checking that all defects added are admissable defects
					count=0
					do j=1,numSpecies
						if(defectStore%defectType(j) /= 0) then
							count=count+1
						end if
					end do
					if(count > 2) then
						write(*,*) 'Error adding unadmissible defect to cascade'
						write(*,*) defectStore%defectType
						call MPI_ABORT(comm,ierr)
					end if

					!add DefectTemp+products to the system
					defectCurrent=>CascadeCurrent%localDefects(defectStore%cellNumber)
					count=0
					do j=1,numSpecies
						if(defectStore%defectType(j)==0) then
							count=count+1
						end if
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
								end if
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
							end if
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
					end if

					defectStore=>defectStore%next
				end do

				!Step 5: set up defectUpdate for all defects in fine mesh
				do i=1,numCellsCascade
					defectTemp=>CascadeCurrent%localDefects(i)%next

					do while(associated(defectTemp))
						!*******************************************************************************
						!create a new element in defectUpdate and assign all variables
						!Note: this will be done for each member of the cascade, even if more than
						!one defect of this type is added to this volume element in the fine mesh.
						!This means that the reaction list update will happen multiple times (redundant)
						!for some defects in some volume elements in the fine mesh.
						!*******************************************************************************
						allocate(defectUpdateCurrent%next)
						defectUpdateCurrent=>defectUpdateCurrent%next
						allocate(defectUpdateCurrent%defectType(numSpecies))
						do j=1,numSpecies
							defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
						end do
						defectUpdateCurrent%num=1
						defectUpdateCurrent%cellNumber=i
						defectUpdateCurrent%proc=reactionCurrent%taskid(1)
						defectUpdateCurrent%dir=0	!not pointed at a different proc
						defectUpdateCurrent%neighbor=0	!not pointed at a different proc
						defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID	!identify which cascade this defect is inside
						nullify(defectUpdateCurrent%next)

						defectTemp=>defectTemp%next
					end do
				end do

				!*******************************************************************
				!Free memory: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
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
				!Here, we are taking the defects in the initialized fine mesh and combining them with
				!the defects in the cascade according to a probability given by the cascade volume
				!divided by the total fine mesh volume.
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
				nullify(defectStore%next)

				!***************************************************************************************
				!Step 1: Create defectStoreList with defects in cascade only
				!Step 2: For each defect in the coarse mesh, check whether it combines with cascade
				!Step 3: If yes, remove it from the coarse mesh and combine it randomly with one defect in defectStoreList
				!Step 4: Implant defects in defectStoreList in mesh. Remember to correctly update defect update protocol.
				!***************************************************************************************
				!Step 1: add defects of cascadeTemp into defectStoreList
				do i=1, cascadeTemp%numDefectsTotal
					!Create list of defects that need to be added to cell (after casacade mixing has been taken into account)
					allocate(defectStore%next)
					defectStore=>defectStore%next
					allocate(defectStore%defectType(numSpecies))
					do j=1,numSpecies
						defectStore%defectType(j)=cascadeDefectTemp%defectType(j)
					end do
					if(pointDefectToggle=='yes') then
						if(defectStore%defectType(3) > max3DInt) then
							defectStore%defectType(4) = defectStore%defectType(3)
							defectStore%defectType(3) = 0
						end if
					end if
					defectStore%num=1
					defectStore%cellNumber=reactionCurrent%cellNumber(1)
					nullify(defectStore%next)

					cascadeDefectTemp=>cascadeDefectTemp%next
				end do

				!Step 2: For each exited defect in the coarse mesh, check whether it combines with cascade in the defectStoreList
				!defectCurrent=>defectList(reactionCurrent%cellNumber(1))	!exited defects in the CoarseMesh
				!defectPrev=>defectCurrent
				!defectTemp=>defectCurrent%next

				defectTemp=>defectList(reactionCurrent%cellNumber(1))	!exited defects in the CoarseMesh
				defectPrev=>defectTemp
				defectTemp=>defectTemp%next
				if(cascadeVolume > 0d0) then
					do while(associated(defectTemp))

						!Combination: choose a defect in defectStoreList to combine with the defect that exited in the coarseMesh
						mixingEvents=0
						do k=1,defectTemp%num

							!Choose which defect in the cascade will be combined with defectTemp already present in the cell
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
									!defectStore is at the middle/end of defectStoreList
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

				!Step 4: Implant defects in defectStoreList in coarse mesh. (add all cascade mixing products and normal cascade products to the system)
				defectStore=>defectStoreList%next
				do while(associated(defectStore))

					!add DefectTemp+products to the system
					!defectCurrent=>defectList(defectStore%cellNumber)
					defectCurrent=>defectList(reactionCurrent%cellNumber(1))
					count=0
					do j=1,numSpecies
						if(defectStore%defectType(j)==0) then
							count=count+1
						end if
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
								end if
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
						end if
					end if

					defectStore=>defectStore%next
				end do

				!Final step: set up defectUpdate for all defects in the cell
				defectTemp=>defectList(reactionCurrent%cellNumber(1))%next

				do while(associated(defectTemp))
					!*******************************************************************************
					!create a new element in defectUpdate and assign all variables
					!Note: this will be done for each member of the cascade, even if more than
					!one defect of this type is added to this volume element in the cell.
					!This means that the reaction list update will happen multiple times (redundant)
					!for some defects in some volume elements in the cell.
					!*******************************************************************************
					allocate(defectUpdateCurrent%next)
					defectUpdateCurrent=>defectUpdateCurrent%next
					allocate(defectUpdateCurrent%defectType(numSpecies))
					do j=1,numSpecies
						defectUpdateCurrent%defectType(j)=defectTemp%defectType(j)
					end do
					defectUpdateCurrent%num=1
					defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(1)
					defectUpdateCurrent%proc=reactionCurrent%taskid(1)
					defectUpdateCurrent%dir=0	!not pointed at a different proc
					defectUpdateCurrent%neighbor=0	!not pointed at a different proc
					defectUpdateCurrent%cascadeNumber=0	!identify which cascade this defect is inside: 0 for coarse mesh
					nullify(defectUpdateCurrent%next)

					defectTemp=>defectTemp%next
				end do

				!*******************************************************************
				!Free memory: defectStore, defectStoreList, defectCurrent, defectPrev, defectTemp
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
					if(diffusionRandom <= fineLength/grainSize) then
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
							if(myMesh(cascadeCurrent%cellNumber)%neighborProcs(j) /= myProc%taskid .AND. &
									myMesh(cascadeCurrent%cellNumber)%neighborProcs(j)/=-1) then

								!neighboring element not in this proc
								numUpdateSend(j)=numUpdateSend(j)+1
								do l=1,numSpecies
									firstSend(l,numUpdateSend(j),j)=reactionCurrent%products(l,i)
								end do
								!cell number in local mesh
								firstSend(numSpecies+1,numUpdateSend(j),j)=cascadeCurrent%cellNumber
								!number of defects to be increased by 1
								firstSend(numSpecies+2,numUpdateSend(j),j)=1
								!cell number of neighbor (in different proc)
								firstSend(numSpecies+3,numUpdateSend(j),j)=myMesh(cascadeCurrent%cellNumber)%neighbors(j)
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
					call MPI_ISEND(firstSend(1,1,i),(numSpecies+3)*numUpdateSend(i),MPI_INTEGER,myProc%procNeighbor(i),&
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
						end if
					end do
					if(same==numSpecies) then
						exit
					end if
					defectPrev=>defectCurrent
					defectCurrent=>defectCurrent%next
				end do

				!If reactant defect is not in the list, then we have chosen a reaction that shouldn't exist
				if(associated(defectCurrent) .EQV. .FALSE.) then
					write(*,*) '*************************************'
					write(*,*) 'Tried to delete defect that wasnt there fine  step', step, 'proc', myProc%taskid
					write(*,*) 'reactants', reactionCurrent%reactants
					write(*,*) 'products', reactionCurrent%products
					write(*,*) 'cells', reactionCurrent%CellNumber
					write(*,*) 'rate', reactionCurrent%reactionRate
					write(*,*) 'cascade number', cascadeCurrent%cascadeID
					defectTest=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i))
					do while(associated(defectTest))
						write(*,*) defectTest%defectType, defectTest%num, defectTest%cellNumber
						defectTest=>defectTest%next
					end do
					reactionTest=>CascadeCurrent%reactionList(reactionCurrent%cellNumber(i))
					do while(associated(reactionTest))
						write(*,*) 'numReactants',reactionTest%numReactants, 'numProducts', reactionTest%numProducts
						write(*,*) 'reactants', reactionTest%reactants, 'products',reactionTest%products
						write(*,*) 'fineCell', reactionTest%CellNumber
						reactionTest=>reactionTest%next
					end do

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
						nullify(defectUpdateCurrent%next)
						defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
						defectUpdateCurrent%num=1		!used for updating reaction lists
						defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
						defectUpdateCurrent%neighbor=0	!not pointed at a different proc
						allocate(defectUpdateCurrent%defectType(numSpecies))
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
							defectCurrent=>CascadeCurrent%localDefects(reactionCurrent%cellNumber(i+&
									reactionCurrent%numReactants))

							!defectUpdateCurrent tells updateReactionList to update reactions within fine mesh
							defectUpdateCurrent%cascadeNumber=CascadeCurrent%cascadeID
							defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
						end if

						! this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
						! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
						nullify(defectPrev)
						call findDefectInList(defectCurrent, defectPrev, products)

						!*********************************************
						!Update defects
						!*********************************************
						if(associated(defectCurrent)) then !if we aren't at the end of the list
							same=0
							do j=1,numSpecies
								if(defectCurrent%defectType(j)==products(j)) then
									same=same+1
								end if
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
					if(myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j) /= myProc%taskid .AND. &
							myMesh(reactionCurrent%cellNumber(i))%neighborProcs(j) /= -1) then	!neighboring element not in this proc

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
						firstSend(numSpecies+3,numUpdateSend(j),j)=myMesh(reactionCurrent%cellNumber(i))%neighbors(j)

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
					!if(reactionCurrent%reactants(3,1) > max3DInt) then
					!	diffusionRandom=dprand()
					!	if(diffusionRandom <= fineLength/meanFreePath) then
					!		flag=.TRUE.
					!	end if
					!end if
					diffusionRandom=dprand()

					!randomly choose whether to remove this defect from the system according to the mean free path (parameter) and the
					!length of the volume element that the defect is currently in
					if(diffusionRandom <= myMesh(reactionCurrent%cellNumber(1))%length/grainSize) then
						flag=.TRUE.
						!reactionCurrent%numProducts=0 !remove that defect from system
						!we can no longer remove the defect in this way because changing the reaction list
						!is no longer allowed (reaction list is not re-created at each step)
					end if
				end if
			end if

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
									if(myMesh(reactionCurrent%cellNumber(1))%neighbors(j)==&
											reactionCurrent%cellNumber(i+reactionCurrent%numReactants) .AND. &
											myMesh(reactionCurrent%cellNumber(1))%neighborProcs(j)==&
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
								if(myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%&
										neighborProcs(j) /= myProc%taskid .AND. &
										myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%&
												neighborProcs(j)/=-1) then	!neighboring element not in this proc

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
											myMesh(reactionCurrent%cellNumber(i+reactionCurrent%numReactants))%neighbors(j)
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

				do i=1, reactionCurrent%numProducts

					if(reactionCurrent%taskid(i+reactionCurrent%numReactants) == -1) then
						!we have identified a defect that is going to be removed from the system via a free surface
						!Do nothing; no need to add this defect to any lists or update any reaction lists because of it
					else if(reactionCurrent%cellNumber(i+reactionCurrent%numReactants) < 0) then	!diffusion: coarse-to-fine
						!we have identified a defect that is diffusing from the coarse mesh to the fine mesh
						!and therefore neet to treat it differently than diffusion within the coarse mesh.

						allocate(defectUpdateCurrent%next)
						defectUpdateCurrent=>defectUpdateCurrent%next
						nullify(defectUpdateCurrent%next)
						defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
						allocate(defectUpdateCurrent%defectType(numSpecies))

						!When coarse-to-fine reactions are chosen, a random cell within the fine mesh is chosen.
						defectUpdateCurrent%cellNumber=chooseRandomCell()	!chooseRandomCell()>0
						defectUpdateCurrent%num=1		!used for updating reaction lists
						defectUpdateCurrent%dir=0		!tells updateReactionList that we are not in a boundary mesh
						defectUpdateCurrent%neighbor=0	!not pointed at a different proc

						!In coarse-to-fine reactions, the cascade number is stored in reactionCurrent%cellNumber
						!and identified with a negative. Therefore we make it positive again.
						defectUpdateCurrent%cascadeNumber=-reactionCurrent%cellNumber(i+reactionCurrent%numReactants)

						do j=1,numSpecies
							products(j)=reactionCurrent%products(j,i)
							defectUpdateCurrent%defectType(j)=reactionCurrent%products(j,i)
						end do

						!Point CascadeCurrent at the correct cascade
						CascadeCurrent=>ActiveCascades
						do while(associated(CascadeCurrent))
							if(CascadeCurrent%cascadeID==-reactionCurrent%cellNumber(i+&
									reactionCurrent%numReactants)) then
								exit
							end if
							CascadeCurrent=>CascadeCurrent%next
						end do

						!If we went through the previous loop without pointing at a cascade, then we have
						!a problem with placing the coarse-to-fine miration in the correct cascade.
						if(.NOT. associated(CascadeCurrent)) then
							write(*,*) 'error coarse-to-fine cascade not correctly pointed'
						end if

						!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
						! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
						nullify(defectPrev)
						!Point defectCurrent at the defect list within the correct element in CascadeCurrent
						defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
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

						!update the direction in defectUpdateCurrent (used for updating reaction lists)
						!create a new element in defectUpdate and assign all variables except for num and dir (will do later)
						do j=1,6

							!if myBoundary at this cell has defectList associated, then it is a boundary element for a local element in this cell.
							!NOTE: this will not work in the case of a non-uniform mesh (because there could be more
							!than one local neighbor, will have to change to localNeighbor())
							if(myBoundary(reactionCurrent%cellNumber(i+reactionCurrent%numReactants),&
									j)%localNeighbor /= 0 .AND. myBoundary(reactionCurrent%cellNumber(i+&
									reactionCurrent%numReactants),j)%proc == reactionCurrent%taskid(i+&
									reactionCurrent%numReactants)) then

								!Create next entry in defectUpdate (used for updating reaction lists)
								allocate(defectUpdateCurrent%next)
								defectUpdateCurrent=>defectUpdateCurrent%next
								defectUpdateCurrent%proc=reactionCurrent%taskid(i+reactionCurrent%numReactants)
								allocate(defectUpdateCurrent%defectType(numSpecies))
								defectUpdateCurrent%cellNumber=reactionCurrent%cellNumber(i+&
										reactionCurrent%numReactants)
								nullify(defectUpdateCurrent%next)
								defectUpdateCurrent%dir=j
								defectUpdateCurrent%num=1	!used for updating reaction lists
								defectUpdateCurrent%cascadeNumber=0	!not inside a cascade

								!This information is used to tell us which element in the local
								!processor bounds this element in myBoundary
								defectUpdateCurrent%neighbor=myBoundary(reactionCurrent%cellNumber(i+&
										reactionCurrent%numReactants),j)%localNeighbor

								do l=1,numSpecies
									products(l)=reactionCurrent%products(l,i)
									defectUpdateCurrent%defectType(l)=reactionCurrent%products(l,i)
								end do

								!Find defect in defect list on myBoundary
								nullify(defectPrev)
								defectCurrent=>myBoundary(reactionCurrent%cellNumber(i+&
										reactionCurrent%numReactants),j)%defectList

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
							write(*,*) 'proc',myProc%taskid,defectTemp%defectType, defectTemp%num
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
				if (myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(k) /= myProc%taskid .AND. &
						myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(k) /= myProc%procNeighbor(i) .AND. &
						myMesh(firstRecv(numSpecies+1,j,i))%neighborProcs(k) /= -1) then

					!Add this defect to a final buffer
					numUpdateFinal(k)=numUpdateFinal(k)+1

					do m=1,numSpecies
						finalSend(m,numUpdateFinal(k),k)=firstRecv(m,j,i)
					end do
					finalSend(numSpecies+1,numUpdateFinal(k),k)=firstRecv(numSpecies+1,j,i)	!local cell number
					finalSend(numSpecies+2,numUpdateFinal(k),k)=firstRecv(numSpecies+2,j,i) !This should be +1 only
					finalSend(numSpecies+3,numUpdateFinal(k),k)=myMesh(firstRecv(numSpecies+1,j,i))%neighbors(k)	!cell number in neighboring proc

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

!***************************************************************************************************
!>Subroutine updateReactionList(defectUpdate)
!updates reaction rates and creates/deletes reactions according to defect numbers that have changed:
!This subroutine does the following:
!1) Find list of affected volume elements (and their procs) - separate into list of neighbors and list of elements involved in reaction
!1a) Make list of all defect types in each volume element that need to be checked
!2)	Delete all diffusion reactions associated with reactants and products in neighboring volume elements
!3) Delete all reactions associated with reactants and products in volume elements associated with reaction
!4) Add all diffusion reactions associated with reactants and products in neighboring volume elements
!5) Add all reactions associated with reactants and prodcuts in volume elements associated with reaction
!***************************************************************************************************
subroutine updateReactionList(defectUpdate)
	use mod_structures
	use mod_globalVariables
	use mod_reactionrates
	implicit none
	include 'mpif.h'

	type(defectUpdateTracker), pointer, intent(inout) :: defectUpdate
	type(defectUpdateTracker), pointer :: defectUpdateCurrent, defectUpdatePrev, defectUpdateNext
	type(defect), pointer :: defectCurrent
	type(cascade), pointer :: cascadeCurrent
	integer :: i, j, dir, defectTemp(numSpecies)
	integer :: localGrainID, neighborGrainID

	!just for testing
	type(defect), pointer :: defectTest
	integer :: same

	nullify(defectTest)

	nullify(cascadeCurrent)
	nullify(defectUpdateCurrent)
	nullify(defectUpdatePrev)
	nullify(defectUpdateNext)

	defectUpdateCurrent=>defectUpdate%next	!the first is zero
	do while(associated(defectUpdateCurrent))

		if(defectUpdateCurrent%num /= 1 .AND. defectUpdateCurrent%num /= -1) then
			!we have an error; all defectUpdateCurrent members should have num = +/-1 (indicates adding or removing)
			write(*,*) 'error defectUpdateCurrent%num not equal to +/- 1', myProc%taskid, defectUpdateCurrent%num
		end if
		do i=1, numSpecies
			defectTemp(i)=defectUpdateCurrent%defectType(i)
		end do

		!if the defect is within the local mesh, update all relevant reactions
		if(defectUpdateCurrent%proc==myProc%taskid) then
			!*******************************************************************************************
			! defectUpdateCurrent%cascadeNumber==0 means a defect has changed in the coarse mesh
			! then update all relevant reactions in the coarse mesh
			!*******************************************************************************************
			if(defectUpdateCurrent%cascadeNumber==0) then	!update reaction in coarse mesh

				!Single-defect reactions associated with defects of type defectTemp (Dissociation、sinkRemoval、impurityTrapping)
				call addSingleDefectReactions(defectUpdateCurrent%cellNumber,defectTemp)

				!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType
				!(Scan over all defect types in the defect list)
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
						call addMultiDefectReactions(defectUpdateCurrent%cellNumber, defectTemp, &
								defectUpdateNext%defectType)
					end if
				end if

				!*******************
				!Diffusion reactions
				!NOTE: diffusion reactions need to be updated in both directions. That is, we need to
				!update the diffusion reactions in this cell as well as in the neighboring cell
				!(thus taking care of both defect addition and defect removal).
				!*******************
				do j=1,6
					if (myMesh(defectUpdateCurrent%cellNumber)%numNeighbors(j)==0) then
						write(*,*) 'error myMesh does not have neighbors in this direction'
					end if

					!Add diffusion reactions from this cell to neighboring cells
					if(polycrystal=='yes') then

						localGrainID=myMesh(defectUpdateCurrent%cellNumber)%material

						!Find the grain ID number of the neighboring volume element
						if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j) /= myProc%taskid .AND. &
								myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j) /= -1) then

							neighborGrainID=myBoundary(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j),j)%material
						else if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j) == -1) then
							neighborGrainID=localGrainID	!free surface release, don't need to do anything special
						else
							neighborGrainID=myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j))%material
						endif

						if(localGrainID==neighborGrainID) then
							!Allow diffusion between elements in the same grain
							call addDiffusionReactions(defectUpdateCurrent%cellNumber, &
									myMesh(defectUpdateCurrent%cellNumber)%neighbors(j),&
									myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j),&
									j,defectTemp)
						else
							!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
							call addDiffusionReactions(defectUpdateCurrent%cellNumber,0,myProc%taskid,-1,j,defectTemp)
						end if
					else
						call addDiffusionReactions(defectUpdateCurrent%cellNumber, &
								myMesh(defectUpdateCurrent%cellNumber)%neighbors(j),&
								myProc%taskid, myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j),j,defectTemp)
					end if

					!If the neighboring volume element is in the same processor as this element, add
					!diffusion reactions from neighboring cells into this cell (changes reaction list in neighboring cell)
					!NOTE: don't need to switch j (direction) here even though diffusion is in opposite
					!direction because dir is only used when diffusing into/out of neighbor processor
					if(mod(j,2)==0) then
						dir=j-1
					else
						dir=j+1
					end if

					if(myMesh(defectUpdateCurrent%cellNumber)%neighborProcs(j)==myProc%taskid) then
						!Don't need to do this for diffusion between different grains
						if(polycrystal=='yes' .AND. myMesh(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j))&
								%material == myMesh(defectUpdateCurrent%cellNumber)%material) then
							!call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j), defectUpdateCurrent%cellNumber, &
							!		myProc%taskid, myProc%taskid, j, defectTemp)
							call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j), &
									defectUpdateCurrent%cellNumber, myProc%taskid, myProc%taskid, dir, defectTemp)
						else if(polycrystal=='no') then
							!call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j), defectUpdateCurrent%cellNumber, &
							!		myProc%taskid, myProc%taskid, j, defectTemp)
							call addDiffusionReactions(myMesh(defectUpdateCurrent%cellNumber)%neighbors(j), &
									defectUpdateCurrent%cellNumber, myProc%taskid, myProc%taskid, dir, defectTemp)
						end if
					end if
				end do

				!***********************************************************************************
				!Diffusion between coarse mesh and fine mesh - coarse to fine
				!For each cascade in the coarse mesh element,  find the number of defects of same type
				!in the fine mesh. Then find the diffusion rate from coarse mesh into fine mesh
				!based on concentration of defects in coarse mesh and fine mesh. Not diffusing
				!into individual cells of fine mesh, as there are too many possibilities (comp.
				!impractical). Cell will be chosen at random.
				!***********************************************************************************
				CascadeCurrent=>ActiveCascades
				do while(associated(CascadeCurrent))

					if(CascadeCurrent%cellNumber==defectUpdateCurrent%cellNumber) then
						call addDiffusionCoarseToFine(defectUpdateCurrent%cellNumber,myProc%taskid,&
								CascadeCurrent,defectTemp)
					end if
					CascadeCurrent=>CascadeCurrent%next
				end do

				!*******************************************************************************************
				! DefectUpdateCurrent%cascadeNumber .NE. 0 means that the defect that has been updated is inside a cascade.
				! Therefore we update all relevant reactions within the correct cascade's fine mesh
				!*******************************************************************************************
			else

				!Point CascadeCurrent at the cascade associated with the chosen reaction
				CascadeCurrent=>ActiveCascades

				do while(associated(CascadeCurrent))
					if(CascadeCurrent%cascadeID==defectUpdateCurrent%cascadeNumber) then
						exit
					end if
					CascadeCurrent=>CascadeCurrent%next
				end do

				!Check to make sure that we have exited at the correct cascade
				if(.NOT. associated(CascadeCurrent)) then
					write(*,*) 'error CascadeCurrent not associated in updateReactionList'
				end if

				!Single-defect reactions associated with defects of type defectTemp in fine mesh
				call addSingleDefectReactionsFine(defectUpdateCurrent%cascadeNumber,&
						defectUpdateCurrent%cellNumber,defectTemp)

				!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
				!all defect types in the defect list)
				defectCurrent=>CascadeCurrent%localDefects(defectUpdateCurrent%cellNumber)
				do while(associated(defectCurrent))
					if(defectCurrent%num /= 0) then
						call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, &
								defectUpdateCurrent%cellNumber,defectTemp, defectCurrent%defectType)
					end if
					defectCurrent=>defectCurrent%next
				end do
				!See above for explanation of this section
				defectUpdateNext=>defectUpdateCurrent%next
				if(associated(defectUpdateNext)) then
					if(defectUpdateNext%cellNumber==defectUpdateCurrent%cellNumber) then
						call addMultiDefectReactionsFine(defectUpdateCurrent%cascadeNumber, &
								defectUpdateCurrent%cellNumber, defectTemp, defectUpdateNext%defectType)
					end if
				end if

				!Diffusion reactions
				do j=1,6

					!Add diffusion reactions from this cell into neighboring cells
					call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber,&
							defectUpdateCurrent%cellNumber,cascadeConnectivity(j,defectUpdateCurrent%cellNumber),&
							myProc%taskid, myProc%taskid,j,defectTemp)

					if(mod(j,2)==0) then
						dir=j-1
					else
						dir=j+1
					end if

					!Add diffusion reactions from neighboring cells into this cell. (both directions)
					if(cascadeConnectivity(j, defectUpdateCurrent%cellNumber) /= 0) then
						!call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber, cascadeConnectivity(j, defectUpdateCurrent%cellNumber), &
						!	defectUpdateCurrent%cellNumber,myProc%taskid, myProc%taskid,j,defectTemp)
						call addDiffusionReactionsFine(defectUpdateCurrent%cascadeNumber,&
								cascadeConnectivity(j,defectUpdateCurrent%cellNumber),&
								defectUpdateCurrent%cellNumber,myProc%taskid, myProc%taskid,dir,defectTemp)
					end if
				end do
			end if

			!if the defect is within the boundary mesh, only update the diffusion reactions for cells touching that boundary element
		else

			!If we are in this section, this means that a defect has changed in the boundary to this
			!volume element but not IN the volume element. Therefore we only need to update one diffusion
			!reaction in one direction, the direction of the element with the changed defect. No defects
			!have changed in this volume element.
			if(defectUpdateCurrent%neighbor==-1) then
				write(*,*) 'error neighbor not assigned for diffusion reactions into boundary'
				call MPI_ABORT(comm,ierr)
			end if

			if(polycrystal=='yes') then

				localGrainID=myMesh(defectUpdateCurrent%neighbor)%material
				neighborGrainID=myBoundary(defectUpdateCurrent%cellNumber,defectUpdateCurrent%dir)%material
				if(localGrainID==neighborGrainID) then
					!Allow diffusion between elements in the same grain
					call addDiffusionReactions(defectUpdateCurrent%neighbor, defectUpdateCurrent%cellNumber, &
							myProc%taskid, defectUpdateCurrent%proc, defectUpdateCurrent%dir, defectTemp)
				else
					!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
					call addDiffusionReactions(defectUpdateCurrent%neighbor, 0, myProc%taskid, -1, &
							defectUpdateCurrent%dir, defectTemp)
				end if
			else
				call addDiffusionReactions(defectUpdateCurrent%neighbor, defectUpdateCurrent%cellNumber, &
						myProc%taskid, defectUpdateCurrent%proc, defectUpdateCurrent%dir, defectTemp)
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
		end if
	end do
end subroutine
