!***********************************************************************
!> Subroutine deallocateBoundaryDefectList(): deallocates boundary defect lists
!***********************************************************************
subroutine deallocateBoundaryDefectList()
	use mod_structures
	use mod_globalVariables
	implicit none

	integer :: cell, dir
	type(defect), pointer :: defectCurrent, defectPrev

	do cell=1,numCells
		do dir=1,6
			if(myMesh(cell)%neighborProcs(dir) /= myProc%taskid .AND. myMesh(cell)%neighborProcs(dir) /= -1) then !neighbor cell is in different proc, not free surface (taskid=-1)

				if(myMesh(cell)%neighborProcs(dir) /= myProc%procNeighbor(dir)) then
					write(*,*) 'error neighbor not correct during boundary mesh initialization'
				else
					defectCurrent=>myBoundary(myMesh(cell)%neighbors(dir),dir)%defectList%next

					do while(associated(defectCurrent))
						defectPrev=>defectCurrent
						defectCurrent=>defectCurrent%next
						if(allocated(defectPrev%defectType)) then
							deallocate(DefectPrev%defectType)
						end if
						deallocate(defectPrev)
					end do
				end if
			end if
		end do
	end do

end subroutine

!***********************************************************************
!> Subroutine deallocateCascadeList(): deallocates all stored cascade data
!(read in from file) in cascade list
!***********************************************************************
Subroutine deallocateCascadeList()
	use mod_structures
	use mod_globalVariables
	implicit none

	type(cascadeEvent), pointer :: CascadeTemp, CascadePrev
	type(cascadeDefect), pointer :: DefectTemp, DefectPrev

	CascadeTemp=>CascadeList
	do while(associated(CascadeTemp))
		DefectTemp=>CascadeTemp%listOfDefects
		do while(associated(DefectTemp))
			DefectPrev=>DefectTemp
			DefectTemp=>DefectTemp%next
			if(allocated(defectPrev%defectType)) then
				deallocate(defectPrev%defectType)
			endif
			deallocate(defectPrev)
		end do
		CascadePrev=>CascadeTemp
		CascadeTemp=>CascadeTemp%nextCascade
		deallocate(CascadePrev)
	end do

end subroutine

!***********************************************************************
!> Subroutine deallocateMaterialInput(): deallocates material input data
!(binding and migration energies, etc), include:
!FormSingle(:,:), DiffSingle(:,:), DiffFunc(:,:), BindSingle(:,:), BindFunc(:,:)
!DissocReactions(:,:), DiffReactions(:,:), SinkReactions(:,:), ImpurityReactions(:,:), ClusterReactions(:,:), ImplantReactions(:,:)
!***********************************************************************
subroutine deallocateMaterialInput()
	use mod_structures
	use mod_globalVariables
	implicit none

	integer :: i

	do i=1, numSingleForm
		deallocate(FormSingle(i)%defectType)
	end do

	do i=1,numFuncDiff
		deallocate(DiffFunc(i)%defectType)
		deallocate(DiffFunc(i)%min)
		deallocate(DiffFunc(i)%max)
		deallocate(DiffFunc(i)%parameters)
	end do

	do i=1,numSingleDiff
		deallocate(DiffSingle(i)%defectType)
	end do

	do i=1,numFuncBind
		deallocate(BindFunc(i)%defectType)
		deallocate(BindFunc(i)%product)
		deallocate(BindFunc(i)%min)
		deallocate(BindFunc(i)%max)
		deallocate(BindFunc(i)%parameters)
	end do

	do i=1,numSingleBind
		deallocate(BindSingle(i)%defectType)
		deallocate(BindSingle(i)%product)
	end do

	do i=1,numDissocReac
		deallocate(DissocReactions(i)%reactants)
		deallocate(DissocReactions(i)%products)
		deallocate(DissocReactions(i)%min)
		deallocate(DissocReactions(i)%max)
	end do

	do i=1,numDiffReac
		deallocate(DiffReactions(i)%reactants)
		deallocate(DiffReactions(i)%products)
		deallocate(DiffReactions(i)%min)
		deallocate(DiffReactions(i)%max)
	end do

	do i=1,numSinkReac
		deallocate(SinkReactions(i)%reactants)
		!deallocate(SinkReactions(i)%products)	!no products
		deallocate(SinkReactions(i)%min)
		deallocate(SinkReactions(i)%max)
	end do

	do i=1,numImpurityReac
		deallocate(ImpurityReactions(i)%reactants)
		deallocate(ImpurityReactions(i)%products)
		deallocate(ImpurityReactions(i)%min)
		deallocate(ImpurityReactions(i)%max)
	end do

	do i=1,numClusterReac
		deallocate(ClusterReactions(i)%reactants)
		deallocate(ClusterReactions(i)%products)
		deallocate(ClusterReactions(i)%min)
		deallocate(ClusterReactions(i)%max)
	end do

	do i=1,numImplantReac
		if(allocated(ImplantReactions(i)%reactants)) deallocate(ImplantReactions(i)%reactants)
		if(allocated(ImplantReactions(i)%products)) deallocate(ImplantReactions(i)%products)
		if(allocated(ImplantReactions(i)%min)) deallocate(ImplantReactions(i)%min)
		if(allocated(ImplantReactions(i)%max)) deallocate(ImplantReactions(i)%max)
	end do

	if(allocated(FormSingle)) deallocate(FormSingle)
	if(allocated(DiffFunc)) deallocate(DiffFunc)
	if(allocated(DiffSingle)) deallocate(DiffSingle)
	if(allocated(BindFunc)) deallocate(BindFunc)
	if(allocated(BindSingle)) deallocate(BindSingle)
	if(allocated(DissocReactions)) deallocate(DissocReactions)
	if(allocated(DiffReactions)) deallocate(DiffReactions)
	if(allocated(SinkReactions)) deallocate(SinkReactions)
	if(allocated(ImpurityReactions)) deallocate(ImpurityReactions)
	if(allocated(ClusterReactions)) deallocate(ClusterReactions)
	if(allocated(ImplantReactions)) deallocate(ImplantReactions)

end subroutine

!***********************************************************************
!> Subroutine deallocateDefectList(): deallocates defects in the coarse mesh
!***********************************************************************
subroutine deallocateDefectList()
	use mod_structures
	use mod_globalVariables
	implicit none

	type(defect), pointer :: defectCurrent, defectPrev
	integer :: cell, i, j

	do cell=1,numCells
		defectCurrent=>defectList(cell)%next

		do while(associated(defectCurrent))
			defectPrev=>defectCurrent
			defectCurrent=>defectCurrent%next
			if(allocated(defectPrev%defectType)) then
				deallocate(defectPrev%defectType)
			end if
			deallocate(defectPrev)
		end do

	end do

	do cell=1,numCells
		defectCurrent=>defectList(cell)
		if(allocated(defectCurrent%defectType)) then
			deallocate(defectCurrent%defectType)
		end if
	end do
	deallocate(defectList)

end subroutine

!***********************************************************************
!> Subroutine deallocateReactionList(): deallocates reactions in the coarse mesh
!***********************************************************************
subroutine deallocateReactionList()
	use mod_structures
	use mod_globalVariables
	implicit none

	type(reaction), pointer :: reactionCurrent, reactionPrev
	integer :: cell, i, j

	do cell=1,numCells
		reactionCurrent=>reactionList(cell)%next

		do while(associated(reactionCurrent))

			reactionPrev=>reactionCurrent
			reactionCurrent=>reactionCurrent%next
			if(allocated(reactionPrev%reactants)) then
				deallocate(reactionPrev%reactants)
			end if
			if(allocated(reactionPrev%products)) then
				deallocate(reactionPrev%products)
			end if

			if(allocated(reactionPrev%cellNumber)) then
				deallocate(reactionPrev%cellNumber)
			end if

			if(allocated(reactionPrev%taskid)) then
				deallocate(reactionPrev%taskid)
			end if

			deallocate(reactionPrev)
		end do
	end do

	do cell=1,numCells

		reactionCurrent=>reactionList(cell)
		if(allocated(reactionCurrent%reactants)) then
			deallocate(reactionCurrent%reactants)
		end if

		if(allocated(reactionCurrent%products)) then
			deallocate(reactionCurrent%products)
		end if

		if(allocated(reactionCurrent%cellNumber)) then
			deallocate(reactionCurrent%cellNumber)
		end if

		if(allocated(reactionCurrent%taskid)) then
			deallocate(reactionCurrent%taskid)
		end if
	end do
	deallocate(reactionList)

end subroutine
