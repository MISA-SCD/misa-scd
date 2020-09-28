!***************************************************************************************************
!> subroutine DEBUGcheckForUnadmissibleDefects(reactionCurrent) - checks to see if any defects are not allowed
!!
!! This subroutine peruses the defect list and cascade defect lists and checks to see if any defects
!! are unadmissible (hard coded information on what constitutes an unadmissible defect; in this case
!! we are concerned about defectType= 5 10 150 3 for example.
!!
!! If an unadmissible defect is found, the subroutine outputs the defect list and the reaction that was
!! chosen at that step.
!!
!! INPUT: reactionCurrent, used to display the reaction that was chosen at this step if needed
!!        step, used to display what step we are on
!!
!! OUTPUT: unadmissible defects, if any (on the screen)
!***************************************************************************************************
subroutine DEBUGCheckForUnadmissible(reactionCurrent)
	use mod_constants
use mod_globalVariables
use mod_structures
implicit none

type(reaction), pointer, intent(inout) :: reactionCurrent
type(defect), pointer :: defectCurrent
type(cascade), pointer :: cascadeCurrent
integer :: i, j, count

interface
	subroutine DEBUGPrintReaction(reactionCurrent)
		use mod_structures
		implicit none
		type(reaction), pointer, intent(in) :: reactionCurrent
	end subroutine
end interface

!Check coarse mesh
do i=1,numCells
	defectCurrent=>defectList(i)
	count=0
	do j=1,SPECIES
		if(defectCurrent%defectType(j) /= 0) then
			count=count+1
		endif
	end do
	
	!No admissible defects have 3 or 4 nonzero numbers in defectType()
	if(count == 3 .OR. count == 4) then
		write(*,*) 'Unadmissible defect found coarse mesh cell', i, 'step', step
		call DEBUGPrintReaction(reactionCurrent)
		call DEBUGPrintDefects()
	endif
	defectCurrent=>defectCurrent%next
end do

!Check fine mesh
CascadeCurrent=>ActiveCascades

do while(associated(CascadeCurrent))
	do i=1,numCellsCascade
		defectCurrent=>CascadeCurrent%localDefects(i)
		count=0
		do j=1,SPECIES
			if(defectCurrent%defectType(j) /= 0) then
				count=count+1
			endif
		end do
		
		!No admissible defects have 3 or 4 nonzero numbers in defectType()
		if(count == 3 .OR. count == 4) then
			write(*,*) 'Unadmissible defect found in fine mesh cell', i, 'step', step, 'Cascade', CascadeCurrent%cascadeID
			call DEBUGPrintReaction(reactionCurrent)
			call DEBUGPrintDefects()
		endif
		defectCurrent=>defectCurrent%next
	end do
	CascadeCurrent=>CascadeCurrent%next
end do

end subroutine

!*****************************************************************************************
!> Subroutine DEBUG Print Defect Update - prints defectUpdate, which is used to identify which
!!defects have been updated (used to update reaction list)
!*****************************************************************************************
subroutine DEBUGPrintDefectUpdate(defectUpdate)
	use mod_constants
use mod_globalVariables
use mod_structures
implicit none

type(defectUpdateTracker), pointer, intent(in) :: defectUpdate
type(defectUpdateTracker), pointer :: defectUpdateCurrent
	
!output list of defects to update (for updating reaction list)
if(myProc%taskid==MASTER) then
	defectUpdateCurrent=>defectUpdate%next
	do while(Associated(defectUpdateCurrent))
		write(*,*) 'update defect'
		write(*,*) defectUpdateCurrent%defectType, 'cell', defectUpdateCurrent%cellNumber, 'num', defectUpdateCurrent%num, &
			'proc', defectUpdateCurrent%proc, 'dir', defectUpdateCurrent%dir
		defectUpdateCurrent=>defectUpdateCurrent%next
	end do
endif

end subroutine

!*****************************************************************************************
!>Subroutine DEBUG Print reaction list - prints all reaction lists at a given Monte Carlo step
!!
!!Prints all reactions in each volume element (using linked reaction list). Have the option
!!to ignore all reactions with rate=0, or to ignore all diffusion reactions, or to ignore
!!implantation reactions (depending on what you are trying to debug)
!!
!!Prints reaction lists in coarse mesh as well as in cascades, if any are present
!*****************************************************************************************
subroutine DEBUGPrintReactionList()
	use mod_constants
use mod_globalVariables
use mod_structures
implicit none

integer :: i, j
type(reaction), pointer :: reactionCurrent
type(cascade), pointer :: CascadeCurrent

!output reaction list
if(myProc%taskid==MASTER) then
	write(*,*) 'processor', myProc%taskid, 'reactions after step', step, 'numCells', numCells
	do i=1, numCells
		write(*,*)  'cell', i
		!!reactionCurrent=>reactionList(i)%next
		!write(*,*) 'neighbors', (myMesh(i)%neighbors(1,j),j=1,6)
		!write(*,*) 'neigProcs', (myMesh(i)%neighborProcs(1,j),j=1,6)
		reactionCurrent=>reactionList(i)
		do while(associated(reactionCurrent))

			!if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then
			!	reactionCurrent=>reactionCurrent%next
			!else
				write(*,*) 'numReactants',reactionCurrent%numReactants,'numProducts',reactionCurrent%numProducts,&
						'cellNumber', reactionCurrent%cellNumber, 'taskid', reactionCurrent%taskid

				do j=1,reactionCurrent%numReactants
					write(*,*) reactionCurrent%reactants(:,j)
				end do

				do j=1,reactionCurrent%numProducts
					write(*,*) reactionCurrent%products(:,j)
				end do

				!write(*,*) 'numProducts', reactionCurrent%numProducts, 'products',reactionCurrent%products

				!do j=1,reactionCurrent%numProducts
				!	write(*,*) (reactionCurrent%products(k,j),k=1,SPECIES)
				!end do

				write(*,*) 'cells and procs'

				do j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
					write(*,*) reactionCurrent%cellNumber(j), reactionCurrent%taskid(j)
				end do

				write(*,*) 'rate', reactionCurrent%reactionRate

				reactionCurrent=>reactionCurrent%next
			!end if
		end do
		write(*,*) '********************************'
	end do
	write(*,*) 'totalRate', totalRate

	write(*,*)
	CascadeCurrent=>ActiveCascades
	do while(associated(CascadeCurrent))
		write(*,*) 'Cascade', CascadeCurrent%CascadeID, 'Coarse cell number', CascadeCurrent%cellNumber
		do i=1,numCellsCascade
			reactionCurrent=>CascadeCurrent%reactionList(i)%next
			if(associated(reactionCurrent)) then
				write(*,*)
				write(*,*)  'cascade cell', i
				write(*,*) 'cascade cell neighbors', cascadeConnectivity(:,i)
	
				do while(associated(reactionCurrent))
					!if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then
						!ignore diffusion reactions for now
					!	reactionCurrent=>reactionCurrent%next
					!else
						write(*,*) 'numReactants',reactionCurrent%numReactants,'numProducts',reactionCurrent%numProducts &
						,'cellNumber', reactionCurrent%cellNumber,'taskid', reactionCurrent%taskid
						!write(*,*) 'numReactants', reactionCurrent%numReactants, 'reactants',reactionCurrent%reactants
						
						do j=1,reactionCurrent%numReactants
							write(*,*) reactionCurrent%reactants(:,j)
						end do

						do j=1,reactionCurrent%numProducts
							write(*,*) reactionCurrent%products(:,j)
						end do
						
						!write(*,*) 'numProducts', reactionCurrent%numProducts, 'products',reactionCurrent%products
						
						!do j=1,reactionCurrent%numProducts
						!	write(*,*) (reactionCurrent%products(k,j),k=1,SPECIES)
						!end do
						
						!write(*,*) 'cells and procs'
						
						!do j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
						!	write(*,*) reactionCurrent%cellNumber(j), reactionCurrent%taskid(j)
						!end do
						
						!write(*,*) 'rate', reactionCurrent%reactionRate
						
						reactionCurrent=>reactionCurrent%next
					!end if
				end do
			end if
		end do
		!write(*,*) 'cascade total rate', CascadeCurrent%totalRate, 'total rate', totalRate
		write(*,*) 'next cascade'
		CascadeCurrent=>CascadeCurrent%next
	end do
	
end if

end subroutine

!*****************************************************************************************
!> Subroutine Debug print defects - prints all defects in all volume elements
!!
!! This subroutine prints all defects in the coarse mesh and / or the cascades in the 
!! master processor only for the purposes of debugging the code.
!!
!!Have the option to skip volume elements that are empty.
!*****************************************************************************************
subroutine DEBUGPrintDefects()
	use mod_constants
use mod_structures
use mod_globalVariables
implicit none

integer :: i, j
type(defect), pointer :: defectCurrent
type(Cascade), pointer :: CascadeCurrent

!Output defects in local mesh
!if(myProc%taskid==MASTER) then
	write(*,*) 'processor', myProc%taskid, 'defects after step', step
	do i=1,numCells
		if(oneCascadeGCell == myMesh(i)%globalCell) then
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				write(*,*) defectCurrent%defectType, defectCurrent%num,defectCurrent%cellNumber
				defectCurrent=>defectCurrent%next
			end do
			exit
		end if

	end do
	write(*,*)
!endif

!Output defects in boundary mesh
!if(myProc%taskid==MASTER) then
!	write(*,*) 'processor', myProc%taskid, 'boundary after step', step
!	do i=1,numCells
!		do j=1,6
!			do k=1,myMesh(i)%numNeighbors(j)
!				if(myMesh(i)%neighborProcs(k,j) /= myProc%taskid .AND. myMesh(i)%neighborProcs(k,j) /= -1) then
!					defectCurrent=>myBoundary(myMesh(i)%neighbors(k,j),j)%defectList%next
!					if(associated(defectCurrent)) then
!						write(*,*) 'proc', myMesh(i)%neighborProcs(k,j),'dir',j, 'element', myMesh(i)%neighbors(k,j)
!					end if
!					do while(associated(defectCurrent))
!						write(*,*) defectCurrent%defectType, defectCurrent%cellNumber, defectCurrent%num
!						defectCurrent=>defectCurrent%next
!					end do
!				end if
!			end do
!		end do
!	end do
!	write(*,*)
!end if

!Output defects in fine mesh
!if(myProc%taskid==MASTER) then
	write(*,*) 'processor', myProc%taskid, 'cascade defects after step', step
	CascadeCurrent=>ActiveCascades
	do while(associated(CascadeCurrent))
		write(*,*) 'CascadeID', CascadeCurrent%cascadeID, 'coarse cell', CascadeCurrent%cellNumber
		do i=1,numCellsCascade
			defectCurrent=>CascadeCurrent%localDefects(i)
			do while(associated(defectCurrent))
				write(*,*) defectCurrent%defectType, defectCurrent%num, defectCurrent%cellNumber
				defectCurrent=>defectCurrent%next
			end do
		end do
		
		!do i=1,numCellsCascade
		!	write(*,*) 'cell',i,'cell total rate',CascadeCurrent%totalRate(i),'processor total rate',totalRate
		!end do

		!write(*,*) 'cascade reaction limit', CascadeReactionLimit
		CascadeCurrent=>CascadeCurrent%next
	end do
	write(*,*)
!end if

end subroutine

!***********************************************************************
!> Subroutine debug print reaction - outputs the reaction chosen at a given step in the 
!master processor. Used for debugging.
!***********************************************************************
subroutine DEBUGPrintReaction(reactionCurrent)
	use mod_constants
use mod_globalVariables
use mod_structures
implicit none

type(reaction), pointer, intent(in) :: reactionCurrent
double precision, external :: totalRateCheck

if(myProc%taskid==MASTER) then	
	write(*,*) 'reaction chosen processor', myProc%taskid, 'step', step
	if(associated(reactionCurrent)) then

		write(*,*) '**************************************choosenReaction*****************************'
		write(*,*) 'numReactants', reactionCurrent%numReactants,'numProducts', reactionCurrent%numProducts

		if(reactionCurrent%numReactants > 0) then
			write(*,*) 'reactants', reactionCurrent%reactants
		end if
		if(reactionCurrent%numProducts > 0) then
			write(*,*) 'products', reactionCurrent%products
		end if
		!if(reactionCurrent%numReactants == -10) then
		!	write(*,*) 'Cascade implantation chosen', numImplantEvents, 'cell', reactionCurrent%cellNumber
		!end if
		write(*,*) 'cells', reactionCurrent%cellNumber, 'procs', reactionCurrent%taskid, 'rate', reactionCurrent%reactionRate
		write(*,*) 'totalRate', totalRate, 'totalRateCell', totalRateVol(reactionCurrent%cellNumber)
		write(*,*) 'totalRateCheck', totalRateCheck(), 'totalRateCellCheck', totalRateVol(reactionCurrent%cellNumber)
	else
		write(*,*) '**************************************choosenReaction*****************************'
		write(*,*) 'null event chosen'
	end if

end if
write(*,*)

end subroutine

!***********************************************************************
!> Subroutine debug print defectList - outputs defects in the coarse meshes of this processor.
!***********************************************************************
subroutine DEBUGPrintDefectList()
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer :: i
	type(defect), pointer :: defectCurrent

	!output reaction list
	if(myProc%taskid==1) then
		write(*,*) 'processor', myProc%taskid, 'reactions after step', step, 'numCells', numCells
		do i=1, numCells
			!!reactionCurrent=>reactionList(i)%next
			defectCurrent=>defectList(i)
			write(*,*) '************defectList********************'
			write(*,*) 'cellNumber', i, 'processor', myProc%taskid
			do while(associated(defectCurrent))
				write(*,*) 'defectType', defectCurrent%defectType, 'numDefects', defectCurrent%num

				defectCurrent=>defectCurrent%next

			end do
		end do
	end if

end subroutine
