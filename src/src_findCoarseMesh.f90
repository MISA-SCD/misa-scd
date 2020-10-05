!***************************************************************************************************
!>Subroutine findDefectInList(defectCurrent, defectPrev, products)
!points defectCurrent at the appropriate defect in a linked list
!Inputs: products
!Outputs: defectCurrent and defectPrev pointed at location in list
!***************************************************************************************************
subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	type(defect), pointer, intent(inout) :: defectCurrent, defectPrev
	integer, intent(in) :: products(SPECIES)
	integer :: same, j

	if(.NOT. associated(defectCurrent)) then
		write(*,*) 'error defectCurrent not associated in findDefectInList'
	end if

	outer:do while(associated(defectCurrent))
		same=0
		inner: do j=1,SPECIES
			if(defectCurrent%defectType(j)==products(j)) then
				same=same+1
			else if(defectCurrent%defectType(j) > products(j)) then
				exit inner
			end if
		end do inner

		if(same==SPECIES) then
			exit outer		!defectCurrent points to the defect that we are trying to add
		else if(same==j-1) then
			if(defectCurrent%defectType(j) > products(j)) then
				exit outer		!defectCurrent points one past the defect we are trying to add and defectPrev points one before
			end if
		else
			defectPrev=>defectCurrent
			defectCurrent=>defectCurrent%next
		end if
	end do outer

end subroutine

!***************************************************************************************************
!> Subroutine findNumDefect(defectType, cellNumber)
!finds the number of defects of a given type inside a certain volume element in the local mesh
!Inputs: defectType, cellNumber
!Outputs: returns number of defects of type in cell
!***************************************************************************************************
integer function findNumDefect(defectType, cellNumber)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	integer, intent(in) :: defectType(SPECIES)
	integer, intent(in) :: cellNumber
	type(defect), pointer :: defectCurrent
	integer :: numDefects, i, count

	numDefects=0
	defectCurrent=>defectList(cellNumber)%next

	do while(associated(defectCurrent))
		count=0
		do i=1,SPECIES
			if(defectType(i)==defectCurrent%defectType(i)) then
				count=count+1
			end if
		end do
		if(count==SPECIES) then
			numDefects=defectCurrent%num
			exit
		else
			defectCurrent=>defectCurrent%next
		end if
	end do

	findNumDefect=numDefects
end function

!***************************************************************************************************
!> Subroutine findNumDefectBoundary(defectType, cellNumber, dir)
!finds the number of defects of a given type inside a certain volume element in the boundary mesh
!Inputs: defectType, cellNumber, direction (used to identify correct element of myGhost)
!Outputs: returns number of defects of type in cell
!***************************************************************************************************
integer function findNumDefectBoundary(defectType, cellNumber, dir)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	integer, intent(in) :: defectType(SPECIES), cellNumber, dir
	type(defect), pointer :: defectCurrent
	integer :: numDefects, i, count

	if(cellNumber==-1) then
		write(*,*) 'trying to find number of defects on free surface'
		numDefects=0
	else
		numDefects=0
		defectCurrent=>myGhost(cellNumber,dir)%defectList

		do while(associated(defectCurrent))
			count=0
			do i=1,SPECIES
				if(defectType(i)==defectCurrent%defectType(i)) then
					count=count+1
				endif
			end do
			if(count==SPECIES) then
				numDefects=defectCurrent%num
				exit
			else
				defectCurrent=>defectCurrent%next
			endif
		end do
	endif

	findNumDefectBoundary=numDefects
end function

!***********************************************************************
!> subroutine countReactionsCoarse(reactionsCoarse)
!counts number of reactions in the coarse mesh (all processors)
!***********************************************************************
subroutine countReactionsCoarse(reactionsCoarse)
	use mod_structures
	use mod_globalVariables
	implicit none
	include 'mpif.h'

	integer, intent(inout) :: reactionsCoarse
	type(Reaction), pointer :: reactionCurrent
	integer :: cell, reactionCounter

	reactionCounter=0

	do cell=1,numCells
		reactionCurrent=>reactionList(cell)
		do while(Associated(reactionCurrent))
			if(reactionCurrent%numReactants==0 .AND. reactionCurrent%numProducts==0) then
				!null reaction, do not count
				reactionCurrent=>reactionCurrent%next
			else
				reactionCounter=reactionCounter+1
				reactionCurrent=>reactionCurrent%next
			end if
		end do
	end do

	call MPI_ALLREDUCE(reactionCounter,reactionsCoarse, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
	
end subroutine

!***************************************************************************************************
!> Subroutine updateImplantRateSingleCell(cell)
!updates defect implantation rates due to change in volume
!***************************************************************************************************
subroutine updateImplantRateSingleCell(cell)
	use mod_updatereactions
	use mod_structures
	use mod_globalVariables
	implicit none

	integer, intent(in) :: cell
	integer :: reac
	type(reaction), pointer :: reactionCurrent

	if(implantType=='FrenkelPair') then
		do reac=1,numImplantReac
			if(ImplantReactions(reac)%numReactants==0 .AND. ImplantReactions(reac)%numProducts==2) then
				exit	!we have found FrenkelPair implantation
			end if
		end do

		totalRate=totalRate-reactionList(cell)%reactionRate
		totalRateVol(cell)=totalRateVol(cell)-reactionList(cell)%reactionRate
		reactionList(cell)%reactionRate=findRate_0th(cell, ImplantReactions(reac))

		totalRate=totalRate+reactionList(cell)%reactionRate
		totalRateVol(cell)=totalRateVol(cell)+reactionList(cell)%reactionRate

	else if(implantType=='Cascade') then
		!search ImplantList for cascade reactions
		do reac=1,numImplantReac
			if(ImplantReactions(reac)%numReactants==-10 .AND. ImplantReactions(reac)%numProducts==0) then
				exit	!we have found cascade implantation
			end if
		end do

		totalRate=totalRate-reactionList(cell)%reactionRate
		totalRateVol(cell)=totalRateVol(cell)-reactionList(cell)%reactionRate

		if(implantScheme=='MonteCarlo') then
			reactionList(cell)%reactionRate=findRate_0th(cell, ImplantReactions(reac))

			totalRate=totalRate+reactionList(cell)%reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionList(cell)%reactionRate
		else if(implantScheme=='explicit') then
			reactionList(cell)%reactionRate=0d0
		else
			write(*,*) 'error implant scheme in reaction list update'
		end if
		nullify(reactionList(cell)%next)
	else
		write(*,*) 'error unknown implantation type'
	end if

end subroutine

!***************************************************************************************************
!> subroutine resetReactionListSingleCell(cell)
!resets an entire reaction list in a single volume element, which used when a cascade is created or
!deleted within that cell, all reaction rates change because volume changes
!***************************************************************************************************
subroutine resetReactionListSingleCell(cell)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	use mod_updatereactions
	implicit none

	integer, intent(in) :: cell
	type(defect), pointer :: defectCurrent, defectUpdate
	type(cascade), pointer :: cascadeCurrent
	integer :: i, j, dir, defectTemp(SPECIES)
	integer :: localGrainID, neighborGrainID

	!First reset the reaction list in this coarse volume element.
	call clearReactionListSingleCell(cell)

	!Only do this step if we are not in the anneal phase (otherwise we wnat reaction rates to remain zero)
	if(annealIdentify .eqv. .FALSE.) then
		!Reset the reaction rate for cascade implantation (changed due to volume change in the cell)
		call updateImplantRateSingleCell(cell)
	end if

	defectUpdate=>defectList(cell)
	do while(associated(defectUpdate))
		do i=1,SPECIES
			defectTemp(i)=defectUpdate%defectType(i)
		end do

		call update_1st_reactions(cell,defectTemp)

		defectCurrent=>defectList(cell)
		do while(associated(defectCurrent))
			call update_2nd_reactions(cell, defectTemp, defectCurrent%defectType)
			defectCurrent=>defectCurrent%next
		end do

		!Diffusion reactions
		do j=1,6

			if (myMesh(cell)%numNeighbors(j)==0) then
				write(*,*) 'error myMesh does not have neighbors in this direction'
			end if

			!Add diffusion reactions from this cell to neighbors and from neighbors to this cell
			if(polycrystal=='yes') then

				!Find the grain ID number of the volume element we are in
				localGrainID=myMesh(cell)%material

				!Find the grain ID number of the neighboring volume element
				if(myMesh(cell)%neighborProcs(j) /= myProc%taskid .AND. myMesh(cell)%neighborProcs(j) /= -1) then
					neighborGrainID=myGhost(myMesh(cell)%neighbors(j),j)%material
				else if(myMesh(cell)%neighborProcs(j) == -1) then
					neighborGrainID=localGrainID
				else
					neighborGrainID=myMesh(myMesh(cell)%neighbors(j))%material
				end if

				if(localGrainID==neighborGrainID) then
					!Allow diffusion between elements in the same grain
					call update_diff_reactions(cell, myMesh(cell)%neighbors(j),&
							myProc%taskid, myMesh(cell)%neighborProcs(j),j,defectTemp)
				else
					!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
					call update_diff_reactions(cell, 0, myProc%taskid, -1, j, defectTemp)
				end if
			else
				call update_diff_reactions(cell, myMesh(cell)%neighbors(j), myProc%taskid, &
						myMesh(cell)%neighborProcs(j), j, defectTemp)
			end if

			if(mod(j,2)==0) then
				dir=j-1
			else
				dir=j+1
			end if

			!Add diffusion reactions from the neighboring cell into this one
			if(myMesh(cell)%neighborProcs(j)==myProc%taskid) then
				if(polycrystal=='yes' .AND. myMesh(myMesh(cell)%neighbors(j))%material == myMesh(cell)%material) then
					!call addDiffusionReactions(myMesh(cell)%neighbors(j), cell, myProc%taskid, &
					!		myProc%taskid, j, defectTemp)
					call update_diff_reactions(myMesh(cell)%neighbors(j),cell,myProc%taskid,myProc%taskid,dir,&
							defectTemp)
				else if(polycrystal=='no') then
					!call addDiffusionReactions(myMesh(cell)%neighbors(j), cell, myProc%taskid, &
					!		myProc%taskid, j, defectTemp)
					call update_diff_reactions(myMesh(cell)%neighbors(j),cell,myProc%taskid,myProc%taskid,dir,&
							defectTemp)
				end if
			end if
		end do

		!***********************************************************************************
		!Diffusion between coarse mesh and fine mesh - coarse to fine
		!***********************************************************************************
		CascadeCurrent=>ActiveCascades
		do while(associated(CascadeCurrent))
			if(CascadeCurrent%cellNumber==cell) then
				call update_diff_coarseToFine(cell, myProc%taskid, CascadeCurrent, defectTemp)
			end if
			CascadeCurrent=>CascadeCurrent%next
		end do
		defectUpdate=>defectUpdate%next
	end do

end subroutine

!***************************************************************************************************
!>Subroutine clearReactionListSingleCell(cellNumber)
!clears reactions from a reaction list in a single volume element
!***************************************************************************************************
subroutine clearReactionListSingleCell(cellNumber)
	use mod_structures
	use mod_globalVariables
	implicit none

	integer, intent(in) :: cellNumber
	type(reaction), pointer :: reactionCurrent, reactionPrev

	reactionCurrent=>reactionList(cellNumber)%next

	do while(associated(reactionCurrent))
		reactionPrev=>reactionCurrent
		reactionCurrent=>reactionCurrent%next

		totalRate=totalRate-reactionPrev%reactionRate
		totalRateVol(cellNumber)=totalRateVol(cellNumber)-reactionPrev%reactionRate

		deallocate(reactionPrev%reactants)
		if(allocated(reactionPrev%products)) then
			deallocate(reactionPrev%products)
		end if
		deallocate(reactionPrev%taskid)
		deallocate(reactionPrev%cellNumber)
		deallocate(reactionPrev)
	end do
	nullify(reactionList(cellnumber)%next)

end subroutine
