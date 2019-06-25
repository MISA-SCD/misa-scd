!***************************************************************************************************
!>Subroutine Find defect in list - points defectCurrent at the appropriate defect in a linked list
!!
!!This subroutine places defectCurrent on the defect being added to the system if it is already in the system
!!If not, defectCurrent points to the place after the inserted defect and defectPrev points to the place before the inserted defect
!!Ordering: first by He content, then by V content, then by interstitial content
!!
!!Inputs: products
!!Outputs: defectCurrent and defectPrev pointed at location in list
!***************************************************************************************************

subroutine findDefectInList(defectCurrent, defectPrev, products)
use DerivedType
use mod_constants
implicit none

type(defect), pointer :: defectCurrent, defectPrev
integer products(numSpecies), same, j

if(.NOT. associated(defectCurrent)) then
	write(*,*) 'error defectCurrent not associated in findDefectInList'
end if

outer:do while(associated(defectCurrent))
	same=0
	!if(myProc%taskid==MASTER) write(*,*) 'DefectCurrent', defectCurrent%defectType
	inner: do j=1,numSpecies
		if(defectCurrent%defectType(j)==products(j)) then
			same=same+1
		else if(defectCurrent%defectType(j) > products(j)) then
			exit inner
		end if
	end do inner
	
	if(same==numSpecies) then
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
!> Subroutine find number of defects - finds the number of defects of a given type inside a certain volume element in the local mesh
!!
!!This subroutine searches through a linked list for the a defect of a given type
!!and returns the number of defects of that type present in volume element (cellNumber)
!!
!!Inputs: defectType, cellNumber
!!Outputs: returns number of defects of type in cell
!***************************************************************************************************

integer function findNumDefect(defectType, cellNumber)
use DerivedType
use mod_constants
implicit none

type(defect), pointer :: defectCurrent
integer defectType(numSpecies), cellNumber, numDefects, i, count

numDefects=0
defectCurrent=>defectList(cellNumber)

do while(associated(defectCurrent))
	count=0
	do i=1,numSpecies
		if(defectType(i)==defectCurrent%defectType(i)) then
			count=count+1
		endif
	end do
	if(count==numSpecies) then
		numDefects=defectCurrent%num
		exit
	else
		defectCurrent=>defectCurrent%next
	endif
end do

findNumDefect=numDefects
end function

!***************************************************************************************************
!> Subroutine find number of defects boundary - finds the number of defects of a given type inside a certain volume element in the boundary mesh
!!
!!This subroutine searches through a linked list for the a defect of a given type
!!and returns the number of defects of that type present in volume element (cellNumber) - in the boundary mesh
!!
!!Inputs: defectType, cellNumber, direction (used to identify correct element of myBoundary)
!!Outputs: returns number of defects of type in cell
!***************************************************************************************************

integer function findNumDefectBoundary(defectType, cellNumber, dir)
use DerivedType
use mod_constants
implicit none

type(defect), pointer :: defectCurrent
integer defectType(numSpecies), cellNumber, numDefects, i, count, dir

if(cellNumber==-1) then
	write(*,*) 'trying to find number of defects on free surface'
	numDefects=0
else
	numDefects=0
	defectCurrent=>myBoundary(dir, cellNumber)%defectList
	
	do while(associated(defectCurrent))
		count=0
		do i=1,numSpecies
			if(defectType(i)==defectCurrent%defectType(i)) then
				count=count+1
			endif
		end do
		if(count==numSpecies) then
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
!
!> subroutine countReactionsCoarse(reactionsCoarse) - counts number of reactions in the coarse mesh
!!
!! Counts the total number of reactions in the coarse mesh (all processors)
!! and returns the sum.
!!
!! Inputs: none
!! Output: reactionsCoarse, total number of reactions in all coarse
!!         mesh elements including all processors
!
!***********************************************************************

subroutine countReactionsCoarse(reactionsCoarse)

use DerivedType
use mod_constants

implicit none

include 'mpif.h'

type(Reaction), pointer :: reactionCurrent
integer cell
integer reactionsCoarse, reactionCounter, reactionsTemp
integer proc

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
		endif
		
	end do

end do

call MPI_ALLREDUCE(reactionCounter,reactionsCoarse, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
	
end subroutine

!***************************************************************************************************
!
!> Subroutine updateImplantRateSingleCell(cell) - updates defect implantation rates due to change in volume
!!
!! This subroutine updates the implantation rate in a coarse mesh element due to the change in
!! volume from cascade implantation / deletion
!!
!! Input: cell (integer, coarse mesh element id number)
!! Output: updates cascade implantation rate
!
!***************************************************************************************************

subroutine updateImplantRateSingleCell(cell)
use ReactionRates
use DerivedType
use mod_constants
implicit none

integer cell, reac
type(reaction), pointer :: reactionCurrent
integer matNum

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
endif

if(implantType=='FrenkelPair') then
	
	!Do nothing, no fine mesh in Frenkel pair implantation.

else if(implantType=='Cascade') then

	!update reaction rates with cascade implantation. The rate should be given by a function
	!findReactionRate which has as parameters passed to it the cell that this rate is occuring in 
	!and the parameters of the reaction rate as read in from the file (type reactionParameters)
	
	!Part 1: Cascade implantation rate
	
	!search ImplantList for cascade reactions
	do reac=1,numImplantReac(matNum)
		if(ImplantReactions(matNum,reac)%numReactants==-10 .AND. ImplantReactions(matNum,reac)%numProducts==0) then	
			
			!we have found cascade implantation
			exit
		
		endif
	end do
	
	!Update the total reaction rate by subtracting the old reaction rate and adding the new one
	totalRate=totalRate-reactionList(cell)%reactionRate
	totalRateVol(cell)=totalRateVol(cell)-reactionList(cell)%reactionRate
	
	if(implantScheme=='MonteCarlo') then

		reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
		
		!Update total reaction rate for entire processor and for this volume element
		totalRate=totalRate+reactionList(cell)%reactionRate
		totalRateVol(cell)=totalRateVol(cell)+reactionList(cell)%reactionRate

	else if(implantScheme=='explicit') then

		!zero reaction rate for cascade implantation in explicit scheme
		reactionList(cell)%reactionRate=0d0	

	else
		write(*,*) 'error implant scheme in reaction list update'
	endif

	
	nullify(reactionList(cell)%next)


else

	write(*,*) 'error unknown implantation type'

endif

end subroutine

!***************************************************************************************************
!> subroutine resetReactionListSingleCell(cell) - resets an entire reaction list in a single volume element
!
! Resets the reaction list for all reactions within a single cell (used when a cascade is created or
! deleted within that cell, all reaction rates change because volume changes)
!
! Input: cell (integer): cell number
! Outputs: reaction rates for all reactions in a cell
!***************************************************************************************************

subroutine resetReactionListSingleCell(cell)
use mod_constants
use DerivedType
use ReactionRates
implicit none

type(defect), pointer :: defectCurrent, defectUpdate
type(cascade), pointer :: cascadeCurrent
integer numElements, numNeighbors, i, j, k, dir, count, defectTemp(numSpecies), tracker
integer findNumDefectBoundary, findNumDefect
logical flag
integer cell
integer localGrainID, neighborGrainID

!First reset the reaction list in this coarse volume element.
call clearReactionListSingleCell(cell)

!Only do this step if we are not in the anneal phase (otherwise we wnat reaction rates to remain zero)
if(annealIdentify .eqv. .FALSE.) then
	!Reset the reaction rate for cascade implantation (changed due to volume change in the cell)
	call updateImplantRateSingleCell(cell)
end if

defectUpdate=>defectList(cell)

do while(associated(defectUpdate))
	
	do i=1,numSpecies
		defectTemp(i)=defectUpdate%defectType(i)
	end do
	
	!Single-defect reactions associated with defects of type defectTemp	
	call addSingleDefectReactions(cell,defectTemp)

	!Multi-defect reactions associated with defects of type defectTemp and defectCurrent%defectType (Scan over
	!all defect types in the defect list)
	defectCurrent=>defectList(cell)
	do while(associated(defectCurrent))
		call addMultiDefectReactions(cell, defectTemp, defectCurrent%defectType)
		defectCurrent=>defectCurrent%next
	end do

	!Diffusion reactions
	do j=1,6
		
		if (myMesh(cell)%numNeighbors(j)==0) then
			write(*,*) 'error myMesh does not have neighbors in this direction'
		endif
		
		do k=1,myMesh(cell)%numNeighbors(j)

			!Add diffusion reactions from this cell to neighbors and from neighbors to this cell		
			if(polycrystal=='yes') then
			
				!Find the grain ID number of the volume element we are in
				localGrainID=myMesh(cell)%material
				
				!Find the grain ID number of the neighboring volume element
				if(myMesh(cell)%neighborProcs(j,k) /= myProc%taskid .AND. &
					myMesh(cell)%neighborProcs(j,k) /= -1) then
				
					neighborGrainID=myBoundary(j,myMesh(cell)%neighbors(j,k))%material
				else if(myMesh(cell)%neighborProcs(j,k) == -1) then
					neighborGrainID=localGrainID
				else
					neighborGrainID=myMesh(myMesh(cell)%neighbors(j,k))%material
				endif
				
				if(localGrainID==neighborGrainID) then
				
					!Allow diffusion between elements in the same grain
					call addDiffusionReactions(cell, myMesh(cell)%neighbors(j,k),&
						myProc%taskid, myMesh(cell)%neighborProcs(j,k),j,defectTemp)
				
				else
				
					!Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
					call addDiffusionReactions(cell, 0, myProc%taskid, -1, j, defectTemp)													
				
				end if

			else
			
				call addDiffusionReactions(cell, myMesh(cell)%neighbors(j,k), myProc%taskid, &
					myMesh(cell)%neighborProcs(j,k), j, defectTemp)
			
			end if
				
			!Add diffusion reactions from the neighboring cell into this one
			if(myMesh(cell)%neighborProcs(j,k)==myProc%taskid) then
				if(polycrystal=='yes' .AND. myMesh(myMesh(cell)%neighbors(j,k))%material == myMesh(cell)%material) then
					
					call addDiffusionReactions(myMesh(cell)%neighbors(j,k), cell, myProc%taskid, &
						myProc%taskid, j, defectTemp)
						
				else if(polycrystal=='no') then
				
					call addDiffusionReactions(myMesh(cell)%neighbors(j,k), cell, myProc%taskid, &
						myProc%taskid, j, defectTemp)
				
				end if
			end if
		end do
		
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

		if(CascadeCurrent%cellNumber==cell) then
			
!			if(myProc%taskid==MASTER) then
!				write(*,*) 'adding coarse to fine diffusion'
!				read(*,*)
!			endif
			
			call addDiffusionCoarseToFine(cell, myProc%taskid, CascadeCurrent, defectTemp)
		
		endif
		
		CascadeCurrent=>CascadeCurrent%next
		
	end do

	defectUpdate=>defectUpdate%next
	
end do

end subroutine

!***************************************************************************************************
!
!> Subroutine clearReactionListSingleCell(cellNumber) - clears reactions from a reaction list in a single volume element
!!
!! Clears all reactions from a reaction list in a cell, so that these reactions can be reset
!! (used for cascade implantation/removal when all reactions in a coarse mesh element need to be reset)
!!
!! Input: cellNumber
!! Output: none
!! Actions: clears reactions in reactionList(cellNumber)
!
!***************************************************************************************************

subroutine clearReactionListSingleCell(cellNumber)
use DerivedType
use mod_constants
implicit none

integer cellNumber
type(reaction), pointer :: reactionCurrent, reactionPrev

!if(HeDPARatio == 0d0) then

	!leave the first reaction alone (cascade implantation reaction)
	reactionCurrent=>reactionList(cellNumber)%next	

!elseif(HeDPARatio > 0d0) then

	!don't delete first or second reactions (cascade and He implantation)
!	reactionCurrent=>reactionList(cellNumber)%next%next

!else
!	write(*,*) 'Error negative HeDPARatio'
!end if

do while(associated(reactionCurrent))
	reactionPrev=>reactionCurrent
	reactionCurrent=>reactionCurrent%next
	
	!Since we are removing this reaction from the reaction list, subtract it from the total rate.
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

!if(HeDPARatio == 0d0) then
	nullify(reactionList(cellnumber)%next)
!else if(HeDPARatio > 0d0) then
!	nullify(reactionList(cellnumber)%next%next)
!else
!	write(*,*) 'Error negative HeDPARatio'
!end if

end subroutine
