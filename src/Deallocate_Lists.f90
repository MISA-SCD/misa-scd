! $Header: /home/CVS//srscd/src/Deallocate_Lists.f90,v 1.4 2015/05/21 15:14:09 aydunn Exp $
!***********************************************************************
!
!> Subroutine deallocateBoundaryDefectList() - deallocates boundary defect lists
!!
!! This subroutine deallocates all the defects in the boundary mesh at the
!! end of the simulation
!!
!! Inputs: none
!! Outputs: none
!! Action: goes through boundary mesh and deallocates all defects
!
!***********************************************************************

subroutine deallocateBoundaryDefectList()
use DerivedType
use mod_srscd_constants
implicit none

integer cell, dir, k, i
type(defect), pointer :: defectCurrent, defectPrev

do 10 cell=1,numCells
	do 11 dir=1,6
		do 12 k=1,myMesh(cell)%numNeighbors(dir)
			if(myMesh(cell)%neighborProcs(dir,k) .NE. myProc%taskid .AND. &
				myMesh(cell)%neighborProcs(dir,k) .NE. -1) then !neighbor cell is in different proc, not free surface (taskid=-1)
				
				if(myMesh(cell)%neighborProcs(dir,k) .NE. myProc%procNeighbor(dir)) then
					write(*,*) 'error neighbor not correct during boundary mesh initialization'
				else
					!deallocate this boundary element:
					!1) Peruse through defect list and deallocate all defects
					!2) Deallocate entire list
					
					defectCurrent=>myBoundary(dir,myMesh(cell)%neighbors(dir,k))%defectList%next
					
					do 13 while(associated(defectCurrent))
						defectPrev=>defectCurrent
						defectCurrent=>defectCurrent%next
						if(allocated(defectPrev%defectType)) then
							deallocate(DefectPrev%defectType)
						endif
						deallocate(defectPrev)
					13 continue

				endif
			endif
		12 continue
	11 continue
10 continue

end subroutine


!***********************************************************************
!
!> Subroutine deallocateCascadeList() - deallocates all stored cascade data (read in from file) in cascade list
!!
!! This subroutine deallocates all the cascades in the cascade list at the
!! end of the simulation.
!!
!! Inputs: none
!! Outputs: none
!! Action: goes through each cascade list and deallocates all cascades
!! (note: these are the lists used as inputs, not the fine mesh defects)
!
!***********************************************************************

Subroutine deallocateCascadeList()
use DerivedType
use mod_srscd_constants
implicit none

type(cascadeEvent), pointer :: CascadeTemp, CascadePrev
type(cascadeDefect), pointer :: DefectTemp, DefectPrev

CascadeTemp=>CascadeList

do 10 while(associated(CascadeTemp))

	DefectTemp=>CascadeTemp%listOfDefects
	
	do 11 while(associated(DefectTemp))
		
		DefectPrev=>DefectTemp
		DefectTemp=>DefectTemp%next
		
		if(allocated(defectPrev%defectType)) then
			deallocate(defectPrev%defectType)
		endif
		
		deallocate(defectPrev)
	
	11 continue
	
	CascadePrev=>CascadeTemp
	CascadeTemp=>CascadeTemp%nextCascade
	
	deallocate(CascadePrev)

10 continue

end subroutine

!***********************************************************************
!
!> Subroutine deallocateMaterialInput() - deallocates material input data (binding and migration energies, etc)
!!
!! This subroutine deallocates global variables used for calculating
!! allowed reactions:
!!
!!type(diffusionFunction), allocatable :: DiffFunc(:)
!!
!!type(diffusionSingle), allocatable :: DiffSingle(:)
!!
!!type(bindingSingle), allocatable :: BindSingle(:)
!!
!!type(bindingFunction), allocatable :: BindFunc(:)
!!
!!type(reactionParameters), allocatable :: DissocReactions(:), DiffReactions(:), SinkReactions(:)
!!
!!type(reactionParameters), allocatable :: ImpurityReactions(:), ClusterReactions(:), ImplantReactions(:)	
!
!***********************************************************************


subroutine deallocateMaterialInput()
use DerivedType
use mod_srscd_constants
implicit none

integer matNum
integer i

do 20 matNum=1,numMaterials

do 10 i=1,numFuncDiff(matNum)
	deallocate(DiffFunc(matNum,i)%defectType)
	deallocate(DiffFunc(matNum,i)%min)
	deallocate(DiffFunc(matNum,i)%max)
	deallocate(DiffFunc(matNum,i)%parameters)
10 continue

do 11 i=1,numSingleDiff(matNum)
	deallocate(DiffSingle(matNum,i)%defectType)
11 continue

do 12 i=1,numFuncBind(matNum)
	deallocate(BindFunc(matNum,i)%defectType)
	deallocate(BindFunc(matNum,i)%product)
	deallocate(BindFunc(matNum,i)%min)
	deallocate(BindFunc(matNum,i)%max)
	deallocate(BindFunc(matNum,i)%parameters)
12 continue

do 13 i=1,numSingleBind(matNum)
	deallocate(BindSingle(matNum,i)%defectType)
	deallocate(BindSingle(matNum,i)%product)
13 continue

do 14 i=1,numDissocReac(matNum)
	deallocate(DissocReactions(matNum,i)%reactants)
	deallocate(DissocReactions(matNum,i)%products)
	deallocate(DissocReactions(matNum,i)%min)
	deallocate(DissocReactions(matNum,i)%max)
14 continue

do 15 i=1,numDiffReac(matNum)
	deallocate(DiffReactions(matNum,i)%reactants)
	deallocate(DiffReactions(matNum,i)%products)
	deallocate(DiffReactions(matNum,i)%min)
	deallocate(DiffReactions(matNum,i)%max)
15 continue

do 16 i=1,numSinkReac(matNum)
	deallocate(SinkReactions(matNum,i)%reactants)
	!deallocate(SinkReactions(matNum,i)%products)
	deallocate(SinkReactions(matNum,i)%min)
	deallocate(SinkReactions(matNum,i)%max)
16 continue

do 17 i=1,numImpurityReac(matNum)
	deallocate(ImpurityReactions(matNum,i)%reactants)
	deallocate(ImpurityReactions(matNum,i)%products)
	deallocate(ImpurityReactions(matNum,i)%min)
	deallocate(ImpurityReactions(matNum,i)%max)
17 continue

do 18 i=1,numClusterReac(matNum)
	deallocate(ClusterReactions(matNum,i)%reactants)
	deallocate(ClusterReactions(matNum,i)%products)
	deallocate(ClusterReactions(matNum,i)%min)
	deallocate(ClusterReactions(matNum,i)%max)
18 continue

do 19 i=1,numImplantReac(matNum)
!	deallocate(ImplantReactions(matNum,i)%reactants)
!	deallocate(ImplantReactions(matNum,i)%products)
!	deallocate(ImplantReactions(matNum,i)%min)
!	deallocate(ImplantReactions(matNum,i)%max)
19 continue

20 continue

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
!
!> Subroutine deallocateDefectList() - deallocates defects in the coarse mesh
!!
!! This subroutine deallocates all the defects in the coarse mesh at the
!! end of the simulation.
!!
!! Inputs: none
!! Outputs: none
!! Action: goes through each defect list for each volume element and 
!! deallocates all defects
!
!***********************************************************************

subroutine deallocateDefectList()
use DerivedType
use mod_srscd_constants
implicit none

type(defect), pointer :: defectCurrent, defectPrev
integer cell, i, j

do 10 cell=1,numCells
	defectCurrent=>DefectList(cell)%next
	
	do 11 while(associated(defectCurrent))
	
		defectPrev=>defectCurrent
		defectCurrent=>defectCurrent%next
		
		if(allocated(defectPrev%defectType)) then
			deallocate(defectPrev%defectType)
		endif
		
		deallocate(defectPrev)
	
	11 continue

10 continue

do 12 cell=1,numCells
	defectCurrent=>defectList(cell)
	
	if(allocated(defectCurrent%defectType)) then
		deallocate(defectCurrent%defectType)
	endif
12 continue

deallocate(DefectList)

end subroutine

!***********************************************************************
!
!> Subroutine deallocateReactionList() - deallocates reactions in the coarse mesh
!!
!! This subroutine deallocates all the reactions in the coarse mesh at the
!! end of the simulation.
!!
!! Inputs: none
!! Outputs: none
!! Action: goes through each reaction list for each volume element and 
!! deallocates all reactions
!
!***********************************************************************

subroutine deallocateReactionList()
use DerivedType
use mod_srscd_constants
implicit none

type(reaction), pointer :: reactionCurrent, reactionPrev
integer cell, i, j

do 10 cell=1,numCells
	reactionCurrent=>reactionList(cell)%next
	
	do 11 while(associated(reactionCurrent))
	
		reactionPrev=>reactionCurrent
		reactionCurrent=>reactionCurrent%next
		
		if(allocated(reactionPrev%reactants)) then
			deallocate(reactionPrev%reactants)
		endif
		
		if(allocated(reactionPrev%products)) then
			deallocate(reactionPrev%products)
		endif
		
		if(allocated(reactionPrev%cellNumber)) then
			deallocate(reactionPrev%cellNumber)
		endif
		
		if(allocated(reactionPrev%taskid)) then
			deallocate(reactionPrev%taskid)
		endif
		
		deallocate(reactionPrev)
	
	11 continue

10 continue

do 12 cell=1,numCells

	reactionCurrent=>reactionList(cell)
	
	if(allocated(reactionCurrent%reactants)) then
		deallocate(reactionCurrent%reactants)
	endif
	
	if(allocated(reactionCurrent%products)) then
		deallocate(reactionCurrent%products)
	endif
	
	if(allocated(reactionCurrent%cellNumber)) then
		deallocate(reactionCurrent%cellNumber)
	endif
	
	if(allocated(reactionCurrent%taskid)) then
		deallocate(reactionCurrent%taskid)
	endif
	
12 continue

deallocate(reactionList)

end subroutine
