! $Header: /home/CVS//srscd/src/Debug_subroutines.f90,v 1.5 2015/06/09 20:51:19 aydunn Exp $
!***************************************************************************************************
!
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
!
!***************************************************************************************************

subroutine DEBUGCheckForUnadmissible(reactionCurrent, step)
use mod_srscd_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionCurrent
type(defect), pointer :: defectCurrent
type(cascade), pointer :: cascadeCurrent

integer i,j,k,count, step

interface
	subroutine DEBUGPrintReaction(reactionCurrent, step)
	use DerivedType
	type(Reaction), pointer :: reactionCurrent
	integer step
	end subroutine
end interface

!Check coarse mesh
do 10 i=1,numCells
	defectCurrent=>defectList(i)
	count=0
	do 11 j=1,numSpecies
		if(defectCurrent%defectType(j) .NE. 0) then
			count=count+1
		endif
	11 continue
	
	!No admissible defects have 3 or 4 nonzero numbers in defectType()
	if(count == 3 .OR. count == 4) then
		write(*,*) 'Unadmissible defect found coarse mesh cell', i, 'step', step
		call DEBUGPrintReaction(reactionCurrent, step)
		call DEBUGPrintDefects(step)
	endif
	defectCurrent=>defectCurrent%next
10 continue

!Check fine mesh
CascadeCurrent=>ActiveCascades

do 12 while(Associated(CascadeCurrent))
	do 13 i=1,numCellsCascade
		defectCurrent=>CascadeCurrent%localDefects(i)
		count=0
		do 14 j=1,numSpecies
			if(defectCurrent%defectType(j) .NE. 0) then
				count=count+1
			endif
		14 continue
		
		!No admissible defects have 3 or 4 nonzero numbers in defectType()
		if(count == 3 .OR. count == 4) then
			write(*,*) 'Unadmissible defect found in fine mesh cell', i, 'step', step, 'Cascade', CascadeCurrent%cascadeID
			call DEBUGPrintReaction(reactionCurrent, step)
			call DEBUGPrintDefects(step)
		endif
		defectCurrent=>defectCurrent%next
	13 continue
	CascadeCurrent=>CascadeCurrent%next
12 continue

end subroutine

!*****************************************************************************************
!> Subroutine DEBUG Print Defect Update - prints defectUpdate, which is used to identify which
!!defects have been updated (used to update reaction list)
!*****************************************************************************************

subroutine DEBUGPrintDefectUpdate(defectUpdate)
use mod_srscd_constants
use DerivedType
implicit none

type(defectUpdateTracker), pointer :: defectUpdate, defectUpdateCurrent
	
!output list of defects to update (for updating reaction list)
if(myProc%taskid==MASTER) then
	defectUpdateCurrent=>defectUpdate%next
	do 18 while(Associated(defectUpdateCurrent))
		write(*,*) 'update defect'
		write(*,*) defectUpdateCurrent%defectType, 'cell', defectUpdateCurrent%cellNumber, 'num', defectUpdateCurrent%num, &
			'proc', defectUpdateCurrent%proc, 'dir', defectUpdateCurrent%dir
		defectUpdateCurrent=>defectUpdateCurrent%next
	18 continue
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

subroutine DEBUGPrintReactionList(step)
use mod_srscd_constants
use DerivedType
implicit none

integer i, j, k, step
type(reaction), pointer :: reactionCurrent
type(cascade), pointer :: CascadeCurrent

!output reaction list
if(myProc%taskid==MASTER) then
	write(*,*) 'processor', myProc%taskid, 'reactions after step', step, 'numCells', numCells
	do 19 i=1, numCells
		!!reactionCurrent=>reactionList(i)%next
		reactionCurrent=>reactionList(i)
		if(associated(reactionCurrent)) then
			write(*,*)  'cell', i
			write(*,*) 'neighbors', (myMesh(i)%neighbors(j,1),j=1,6)
			write(*,*) 'neigProcs', (myMesh(i)%neighborProcs(j,1),j=1,6)

			do 23 while(associated(reactionCurrent))
				!if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then
				!	!ignore diffusion reactions for now
				!	reactionCurrent=>reactionCurrent%next
				!else
					
					write(*,*) 'numReactants', reactionCurrent%numReactants, 'reactants'
					
					do 20 j=1,reactionCurrent%numReactants
						write(*,*) (reactionCurrent%reactants(j,k),k=1,numSpecies)
					20 continue
					
					write(*,*) 'numProducts', reactionCurrent%numProducts, 'products'
					
					do 21 j=1,reactionCurrent%numProducts
						write(*,*) (reactionCurrent%products(j,k),k=1,numSpecies)
					21 continue
					
					write(*,*) 'cells and procs'
					
					do 22 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
						write(*,*) reactionCurrent%cellNumber(j), reactionCurrent%taskid(j)
					22 continue
					
					write(*,*) 'rate', reactionCurrent%reactionRate
					
					reactionCurrent=>reactionCurrent%next
				
				!endif
			23 continue
!			if(myProc%taskid==MASTER) then
!				read(*,*)
!			endif
			write(*,*) '********************************'
		endif
	19 continue
	
	write(*,*) 'Cascades'
	CascadeCurrent=>ActiveCascades
	do 24 while(associated(CascadeCurrent))
		write(*,*) 'Cascade', CascadeCurrent%CascadeID, 'Coarse cell number', CascadeCurrent%cellNumber
		do 25 i=1,numCellsCascade
			reactionCurrent=>CascadeCurrent%reactionList(i)%next
			if(associated(reactionCurrent)) then
				write(*,*)
				write(*,*)  'cascade cell', i
				write(*,*) 'cascade cell neighbors', (cascadeConnectivity(i,j),j=1,6)
	
				do 26 while(associated(reactionCurrent))
					if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then
						!ignore diffusion reactions for now
						reactionCurrent=>reactionCurrent%next
					else
						
						write(*,*) 'numReactants', reactionCurrent%numReactants, 'reactants'
						
						do 27 j=1,reactionCurrent%numReactants
							write(*,*) (reactionCurrent%reactants(j,k),k=1,numSpecies)
						27 continue
						
						write(*,*) 'numProducts', reactionCurrent%numProducts, 'products'
						
						do 28 j=1,reactionCurrent%numProducts
							write(*,*) (reactionCurrent%products(j,k),k=1,numSpecies)
						28 continue
						
						write(*,*) 'cells and procs'
						
						do 29 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
							write(*,*) reactionCurrent%cellNumber(j), reactionCurrent%taskid(j)
						29 continue
						
						write(*,*) 'rate', reactionCurrent%reactionRate
						
						reactionCurrent=>reactionCurrent%next
					
					endif
				26 continue
			endif
			
		25 continue
		write(*,*) 'cascade total rate', CascadeCurrent%totalRate, 'total rate', totalRate
		write(*,*) 'next cascade'
		CascadeCurrent=>CascadeCurrent%next

	24 continue
	
endif

end subroutine

!*****************************************************************************************
!> Subroutine Debug print defects - prints all defects in all volume elements
!!
!! This subroutine prints all defects in the coarse mesh and / or the cascades in the 
!! master processor only for the purposes of debugging the code.
!!
!!Have the option to skip volume elements that are empty.
!*****************************************************************************************

subroutine DEBUGPrintDefects(step)
use DerivedType
use mod_srscd_constants
implicit none

integer i, j, k, step, count
type(defect), pointer :: defectCurrent
type(Cascade), pointer :: CascadeCurrent

!Output defects and boundary in Master (used as a check for code)
!if(myProc%taskid==MASTER) then
!	write(*,*) 'processor', myProc%taskid, 'defects after step', step
!	do 12 i=1,numCells
!		defectCurrent=>defectList(i)%next
!		do 13 while(associated(defectCurrent))
!			write(*,*) defectCurrent%defectType, defectCurrent%cellNumber, defectCurrent%num
!			
!			!Check for defects of type  0 0 0 0
!			count=0
!			do 99 j=1,numSpecies
!				if(defectCurrent%defectType(j)==0) then
!					count=count+1
!				endif
!			99 continue
!			if(count==numSpecies) write(*,*) 'error: zero defect in coarse mesh'
!			!if(count==numSpecies) read(*,*)
!			
!			defectCurrent=>defectCurrent%next
!		13 continue
!	12 continue
!	!read(*,*)
!	write(*,*)
!endif
	
!if(myProc%taskid==MASTER) then
!	write(*,*) 'processor', myProc%taskid, 'boundary after step', step
!	do 14 i=1,numCells
!		do 15 j=1,6
!			do 16 k=1,myMesh(i)%numNeighbors(j)
!				if(myMesh(i)%neighborProcs(j,k) .NE. myProc%taskid) then
!					if(myMesh(i)%neighborProcs(j,k) .NE. -1) then
!						defectCurrent=>myBoundary(j,myMesh(i)%neighbors(j,k))%defectList%next
!						if(associated(defectCurrent)) then
!							write(*,*) 'proc', myMesh(i)%neighborProcs(j,k),'dir',j, 'element', myMesh(i)%neighbors(j,k)
!						endif
!						do 17 while(associated(defectCurrent))
!							write(*,*) defectCurrent%defectType, defectCurrent%cellNumber, defectCurrent%num
!							defectCurrent=>defectCurrent%next
!						17 continue
!					endif
!				endif
!			16 continue
!		15 continue
!	14 continue
!	write(*,*)
!endif

!Output defects in cascdaes
!if(myProc%taskid==MASTER) then
	write(*,*) 'processor', myProc%taskid, 'cascade defects after step', step
	CascadeCurrent=>ActiveCascades
	do 18 while(associated(CascadeCurrent))
		write(*,*) 'Cascade', CascadeCurrent%cascadeID, 'coarse cell', CascadeCurrent%cellNumber, 'defects'
		do 19 i=1,numCellsCascade
			defectCurrent=>CascadeCurrent%localDefects(i)%next
			do 20 while(Associated(defectCurrent))
				write(*,*) defectCurrent%defectType, defectCurrent%cellNumber, defectCurrent%num
				defectCurrent=>defectCurrent%next
			20 continue
		19 continue
		
		do i=1,numCellsCascade
			write(*,*) 'cell', i, 'total rate', CascadeCurrent%totalRate(i), 'total rate', totalRate
		end do

		write(*,*) 'cascade reaction limit', CascadeReactionLimit
		write(*,*)
		CascadeCurrent=>CascadeCurrent%next
	18 continue
!endif

end subroutine

!***********************************************************************
!
!> Subroutine debug print reaction - outputs the reaction chosen at a given step in the 
!!master processor. Used for debugging.
!
!***********************************************************************

subroutine DEBUGPrintReaction(reactionCurrent, step)
use mod_srscd_constants
use DerivedType
implicit none

type(Reaction), pointer :: reactionCurrent
integer step
double precision totalRateCheck

if(myProc%taskid==MASTER) then	
	write(*,*) 'reaction chosen processor', myProc%taskid, 'step', step
	if(associated(reactionCurrent)) then
		if(reactionCurrent%numReactants .GT. 0) then
			write(*,*) 'reactants', reactionCurrent%reactants
		endif
		if(reactionCurrent%numProducts .GT. 0) then
			write(*,*) 'products', reactionCurrent%products
		endif
		if(reactionCurrent%numReactants == -10) then
			write(*,*) 'Cascade implantation chosen', numImplantEvents, 'cell', reactionCurrent%cellNumber
			!read(*,*)
		endif
		write(*,*) 'cells', reactionCurrent%cellNumber, 'procs', reactionCurrent%taskid, 'rate', reactionCurrent%reactionRate
		!write(*,*) 'totalRate', totalRate, 'totalRateCell', totalRateVol(reactionCurrent%cellNumber)
		!write(*,*) 'totalRateCheck', totalRateCheck(), 'totalRateCellCheck', totalRateVol(reactionCurrent%cellNumber)
	else
		write(*,*) 'null event chosen'
	endif
!		read(*,*) !pause once per step (for debugging only)
endif
write(*,*)

end subroutine
