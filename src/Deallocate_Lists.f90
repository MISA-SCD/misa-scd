!***********************************************************************
!> Subroutine deallocateBoundaryDefectList(): deallocates boundary defect lists
!***********************************************************************
subroutine deallocateBoundaryDefectList()
	use mod_structures
	use mod_constants
	implicit none

	integer cell, dir
	type(defect), pointer :: defectCurrent, defectPrev

	do cell=1,numCells
		do dir=1,6
			if(myMesh(cell)%neighborProcs(1,dir) /= myProc%taskid .AND. myMesh(cell)%neighborProcs(1,dir) /= -1) then !neighbor cell is in different proc, not free surface (taskid=-1)

				if(myMesh(cell)%neighborProcs(1,dir) /= myProc%procNeighbor(dir)) then
					write(*,*) 'error neighbor not correct during boundary mesh initialization'
				else
					defectCurrent=>myBoundary(myMesh(cell)%neighbors(1,dir),dir)%defectList%next

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
	use mod_constants
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
	use mod_constants
	implicit none

	integer matNum
	integer i

	do matNum=1,numMaterials

		do i=1, numSingleForm(matNum)
			deallocate(FormSingle(i,matNum)%defectType)
		end do

		do i=1,numFuncDiff(matNum)
			deallocate(DiffFunc(i,matNum)%defectType)
			deallocate(DiffFunc(i,matNum)%min)
			deallocate(DiffFunc(i,matNum)%max)
			deallocate(DiffFunc(i,matNum)%parameters)
		end do

		do i=1,numSingleDiff(matNum)
			deallocate(DiffSingle(i,matNum)%defectType)
		end do

		do i=1,numFuncBind(matNum)
			deallocate(BindFunc(i,matNum)%defectType)
			deallocate(BindFunc(i,matNum)%product)
			deallocate(BindFunc(i,matNum)%min)
			deallocate(BindFunc(i,matNum)%max)
			deallocate(BindFunc(i,matNum)%parameters)
		end do

		do i=1,numSingleBind(matNum)
			deallocate(BindSingle(i,matNum)%defectType)
			deallocate(BindSingle(i,matNum)%product)
		end do

		do i=1,numDissocReac(matNum)
			deallocate(DissocReactions(i,matNum)%reactants)
			deallocate(DissocReactions(i,matNum)%products)
			deallocate(DissocReactions(i,matNum)%min)
			deallocate(DissocReactions(i,matNum)%max)
		end do

		do i=1,numDiffReac(matNum)
			deallocate(DiffReactions(i,matNum)%reactants)
			deallocate(DiffReactions(i,matNum)%products)
			deallocate(DiffReactions(i,matNum)%min)
			deallocate(DiffReactions(i,matNum)%max)
		end do

		do i=1,numSinkReac(matNum)
			deallocate(SinkReactions(i,matNum)%reactants)
			!deallocate(SinkReactions(matNum,i)%products)	!no products
			deallocate(SinkReactions(i,matNum)%min)
			deallocate(SinkReactions(i,matNum)%max)
		end do

		do i=1,numImpurityReac(matNum)
			deallocate(ImpurityReactions(i,matNum)%reactants)
			deallocate(ImpurityReactions(i,matNum)%products)
			deallocate(ImpurityReactions(i,matNum)%min)
			deallocate(ImpurityReactions(i,matNum)%max)
		end do

		do i=1,numClusterReac(matNum)
			deallocate(ClusterReactions(i,matNum)%reactants)
			deallocate(ClusterReactions(i,matNum)%products)
			deallocate(ClusterReactions(i,matNum)%min)
			deallocate(ClusterReactions(i,matNum)%max)
		end do

		do i=1,numImplantReac(matNum)
			if(allocated(ImplantReactions(i,matNum)%reactants)) deallocate(ImplantReactions(i,matNum)%reactants)
			if(allocated(ImplantReactions(i,matNum)%products)) deallocate(ImplantReactions(i,matNum)%products)
			if(allocated(ImplantReactions(i,matNum)%min)) deallocate(ImplantReactions(i,matNum)%min)
			if(allocated(ImplantReactions(i,matNum)%max)) deallocate(ImplantReactions(i,matNum)%max)
		end do
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
	use mod_constants
	implicit none

	type(defect), pointer :: defectCurrent, defectPrev
	integer cell, i, j

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
	use mod_constants
	implicit none

	type(reaction), pointer :: reactionCurrent, reactionPrev
	integer cell, i, j

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
