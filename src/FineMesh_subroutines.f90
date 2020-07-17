!***************************************************************************************************
!> Function releaseFineMeshDefects(CascadeCurrent)
!releases fine mesh back to coarse mesh when cascade is annealed
!***************************************************************************************************
subroutine releaseFineMeshDefects(CascadeCurrent)
	use mod_constants
	use mod_structures
	implicit none

	!Data structures used:
	type(cascade), pointer :: CascadeCurrent
	type(defect), pointer :: defectCurrentCoarse, defectCurrentFine, defectPrevCoarse
	type(defect), pointer :: defectCurrent, defectPrev
	type(reaction), pointer :: reactionCurrent, reactionPrev
	integer i, j, count

	interface
		subroutine findDefectInList(defectCurrent, defectPrev, products)
			use mod_structures
			use mod_constants
			type(defect), pointer :: defectCurrent, defectPrev
			integer products(numSpecies)
		end subroutine
	end interface

	nullify(defectCurrentFine)
	nullify(defectCurrentCoarse)
	!Step 0: increase coarse mesh volume
	myMesh(CascadeCurrent%cellNumber)%volume=myMesh(CascadeCurrent%cellNumber)%volume+CascadeElementVol*numCellsCascade
	myMesh(CascadeCurrent%cellNumber)%length=(myMesh(CascadeCurrent%cellNumber)%volume)**(1d0/3d0)
	!Step 1: add defects to coarse mesh
	do i=1,numCellsCascade

		defectCurrentFine=>CascadeCurrent%localDefects(i)
		do while(associated(defectCurrentFine))

			nullify(defectPrevCoarse)
			defectCurrentCoarse=>defectList(CascadeCurrent%cellNumber)

			call findDefectInList(defectCurrentCoarse, defectPrevCoarse, defectCurrentFine%defectType)
			if(associated(defectCurrentCoarse)) then

				count=0
				do j=1,numSpecies
					if(defectCurrentCoarse%defectType(j)==defectCurrentFine%defectType(j)) then
						count=count+1
					end if
				end do

				!case 1: this defect already exists in the coarse mesh
				if(count==numSpecies) then
					!Add the fine mesh defects to the coarse mesh defects
					defectCurrentCoarse%num=defectCurrentCoarse%num+defectCurrentFine%num

				!case 2: defect does not exist in coarse mesh and we are in the middle of coarse list
				else if(associated(defectPrevCoarse)) then
					!Create a new defect in the coarse mesh between defectPrevCoarse and defectCurrentCoarse
					nullify(defectPrevCoarse%next)
					allocate(defectPrevCoarse%next)
					nullify(defectPrevCoarse%next%next)
					defectPrevCoarse=>defectPrevCoarse%next
					allocate(defectPrevCoarse%defectType(numSpecies))
					do j=1,numSpecies
						defectPrevCoarse%defectType(j)=defectCurrentFine%defectType(j)
					end do
					defectPrevCoarse%num=defectCurrentFine%num
					defectPrevCoarse%cellNumber=CascadeCurrent%cellNumber
					defectPrevCoarse%next=>defectCurrentCoarse
				else
					write(*,*) 'Error adding fine mesh to coarse mesh, associated'
				end if

			!case 3: we are at the end of the coarse list
			else if(associated(defectPrevCoarse)) then
				!Create a new defect in the coarse mesh at the end of the list
				nullify(defectPrevCoarse%next)
				allocate(defectPrevCoarse%next)
				nullify(defectPrevCoarse%next%next)
				defectPrevCoarse=>defectPrevCoarse%next
				allocate(defectPrevCoarse%defectType(numSpecies))
				do j=1,numSpecies
					defectPrevCoarse%defectType(j)=defectCurrentFine%defectType(j)
				end do
				defectPrevCoarse%num=defectCurrentFine%num
				defectPrevCoarse%cellNumber=CascadeCurrent%cellNumber
			else
				write(*,*) 'Error adding fine mesh to coarse mesh'

			end if
			defectCurrentFine=>defectCurrentFine%next
		end do
	end do

	!Step 2: delete fine mesh
	nullify(defectCurrent)
	nullify(reactionCurrent)

	!Case 1: CascadeCurrent is in the middle of the list of fine meshes
	if(associated(CascadeCurrent%prev) .AND. associated(CascadeCurrent%next)) then
		!Remove CascadeCurrent from the list of cascades by pointing the previous cascade to the next cascade
		CascadeCurrent%prev%next=>CascadeCurrent%next
		CascadeCurrent%next%prev=>CascadeCurrent%prev

		!Deallocate memory: remove defect list, reaction list, and all other data structures
		do i=1,numCellsCascade
			defectCurrent=>CascadeCurrent%localDefects(i)%next
			reactionCurrent=>CascadeCurrent%reactionList(i)%next

			do while(associated(defectCurrent))
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
				deallocate(defectPrev%defectType)
				deallocate(defectPrev)
			end do

			do while(associated(reactionCurrent))
				reactionPrev=>reactionCurrent
				reactionCurrent=>reactionCurrent%next


				totalRate=totalRate-reactionPrev%reactionRate
				deallocate(reactionPrev%reactants)
				if(allocated(reactionPrev%products)) then
					deallocate(reactionPrev%products)
				end if
				deallocate(reactionPrev%cellNumber)
				deallocate(reactionPrev%taskid)
				deallocate(reactionPrev)
			end do
			deallocate(CascadeCurrent%localDefects(i)%defectType)
		end do

		deallocate(CascadeCurrent%totalRate)
		deallocate(CascadeCurrent%reactionList)
		deallocate(CascadeCurrent%localDefects)
		deallocate(CascadeCurrent)

	!case 2: CascadeCurrent is at the end of the list of fine meshes
	else if(associated(CascadeCurrent%prev)) then

		nullify(CascadeCurrent%prev%next)
		do i=1,numCellsCascade

			defectCurrent=>CascadeCurrent%localDefects(i)%next
			reactionCurrent=>CascadeCurrent%reactionList(i)%next

			do while(associated(defectCurrent))
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
				deallocate(defectPrev%defectType)
				deallocate(defectPrev)
			end do

			do while(associated(reactionCurrent))
				reactionPrev=>reactionCurrent
				reactionCurrent=>reactionCurrent%next

				totalRate=totalRate-reactionPrev%reactionRate
				deallocate(reactionPrev%reactants)
				if(allocated(reactionPrev%products)) then
					deallocate(reactionPrev%products)
				end if
				deallocate(reactionPrev%cellNumber)
				deallocate(reactionPrev%taskid)
				deallocate(reactionPrev)
			end do
			deallocate(CascadeCurrent%localDefects(i)%defectType)
		end do

		deallocate(CascadeCurrent%totalRate)
		deallocate(CascadeCurrent%reactionList)
		deallocate(CascadeCurrent%localDefects)
		deallocate(CascadeCurrent)

	!case 3: CascadeCurrent is at the beginning of the list of fine meshes
	else if(associated(CascadeCurrent%next)) then

		ActiveCascades=>CascadeCurrent%next
		nullify(ActiveCascades%prev)
		do i=1,numCellsCascade

			defectCurrent=>CascadeCurrent%localDefects(i)%next
			reactionCurrent=>CascadeCurrent%reactionList(i)%next

			do while(associated(defectCurrent))
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
				deallocate(defectPrev%defectType)
				deallocate(defectPrev)
			end do

			do while(associated(reactionCurrent))
				reactionPrev=>reactionCurrent
				reactionCurrent=>reactionCurrent%next

				totalRate=totalRate-reactionPrev%reactionRate
				deallocate(reactionPrev%reactants)
				if(allocated(reactionPrev%products)) then
					deallocate(reactionPrev%products)
				end if
				deallocate(reactionPrev%cellNumber)
				deallocate(reactionPrev%taskid)
				deallocate(reactionPrev)
			end do
			deallocate(CascadeCurrent%localDefects(i)%defectType)
		end do

		deallocate(CascadeCurrent%totalRate)
		deallocate(CascadeCurrent%reactionList)
		deallocate(CascadeCurrent%localDefects)
		deallocate(CascadeCurrent)

	!case 4: Only one cascade exists
	else

		nullify(ActiveCascades)
		do i=1,numCellsCascade

			defectCurrent=>CascadeCurrent%localDefects(i)%next
			reactionCurrent=>CascadeCurrent%reactionList(i)%next

			do while(associated(defectCurrent))
				defectPrev=>defectCurrent
				defectCurrent=>defectCurrent%next
				deallocate(defectPrev%defectType)
				deallocate(defectPrev)
			end do

			do while(associated(reactionCurrent))
				reactionPrev=>reactionCurrent
				reactionCurrent=>reactionCurrent%next

				totalRate=totalRate-reactionPrev%reactionRate
				deallocate(reactionPrev%reactants)
				if(allocated(reactionPrev%products)) then
					deallocate(reactionPrev%products)
				end if
				deallocate(reactionPrev%cellNumber)
				deallocate(reactionPrev%taskid)
				deallocate(reactionPrev)
			end do
			deallocate(CascadeCurrent%localDefects(i)%defectType)
		end do

		deallocate(CascadeCurrent%totalRate)
		deallocate(CascadeCurrent%reactionList)
		deallocate(CascadeCurrent%localDefects)
		deallocate(CascadeCurrent)
	end if
		
end subroutine

!***************************************************************************************************
!>function findNumDefectFine(CascadeCurrent, defectType, cellNumber)
!!Finds the number of defects of type defectType in CascadeCurrent%localDefects. If none, returns 0
!***************************************************************************************************
integer function findNumDefectFine(CascadeCurrent, defectType, cellNumber)
	use mod_structures
	use mod_constants
	implicit none

	type(defect), pointer :: defectCurrent
	integer defectType(numSpecies), cellNumber, numDefects, i, count
	type(cascade), pointer :: CascadeCurrent

	numDefects=0
	defectCurrent=>CascadeCurrent%localDefects(cellNumber)

	do while(associated(defectCurrent))

		count=0
		do i=1,numSpecies
			if(defectType(i)==defectCurrent%defectType(i)) then
				count=count+1
			end if
		end do

		if(count==numSpecies) then
			numDefects=defectCurrent%num
			exit
		end if
		defectCurrent=>defectCurrent%next
	end do

	findNumDefectFine=numDefects
end function

!***************************************************************************************************
!>function findNumDefectTotalFine(defectType, CascadeCurrent)
!Finds the number of defects of type defectType in CascadeCurrent%localDefects. If none, returns 0
!***************************************************************************************************
integer function findNumDefectTotalFine(defectType, CascadeCurrent)
	use mod_constants
	use mod_structures
	implicit none

	integer defectType(numSpecies)
	type(cascade), pointer :: CascadeCurrent
	type(defect), pointer :: defectCurrent
	integer cell, count

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_constants
			type(cascade), pointer :: CascadeCurrent
			integer defectType(numSpecies), cell
		end function
	end interface

	count=0
	do cell=1,numCellsCascade
		count=count+findNumDefectFine(CascadeCurrent, defectType, cell)
	end do

	findNumDefectTotalFine=count

end function

!***************************************************************************************************
!> integer function findCellWithCoordinatesFineMesh(coordinates)
! FindCellWithCoordinates returns the cell id number that contains the coordinates given in the
! double precision variable coordinates (for defect implantation in fine mesh at beginning of cascade)
!***************************************************************************************************
integer function findCellWithCoordinatesFineMesh(coordinates)
	use mod_constants
	use mod_structures
	implicit none

	double precision coordinates(3)
	integer i,j,k,cellNumber

	do i=1,numxcascade
		if(coordinates(1)+numxcascade*fineLength/2d0 <= i*fineLength) then
			exit
		end if
	end do

	do j=1,numycascade
		if(coordinates(2)+numycascade*fineLength/2d0 <= j*fineLength) then
			exit
		end if
	end do

	do k=1,numzcascade
		if(coordinates(3)+numzcascade*fineLength/2d0 <= k*fineLength) then
			exit
		end if
	end do

	cellNumber=i+(j-1)*numxcascade+(k-1)*(numxcascade*numycascade)

	findCellWithCoordinatesFineMesh=cellNumber

end function

!***************************************************************************************************
!>Function chooseRandomCell()
!Chooses a cell number at random from the fine mesh. (Uniform distribution)
!***************************************************************************************************
integer function chooseRandomCell()
	use randdp
	use mod_constants
	implicit none

	double precision r, a
	integer i

	r=dprand()
	a=0d0
	do i=1,numxCascade*numyCascade*numzCascade
		a=a+1d0/dble(numxCascade*numyCascade*numzCascade)
		if(a > r) then
			exit
		end if
	end do
	chooseRandomCell=i

end function

!***********************************************************************
!>subroutine countReactionsFine(reactionsFine)
!Counts the total number of reactions in the fine mesh (all processors) and returns the sum (used for postprocessing)
!***********************************************************************
subroutine countReactionsFine(reactionsFine)
	use mod_structures
	use mod_constants
	implicit none
	include 'mpif.h'

	type(Reaction), pointer :: reactionCurrent
	type(cascade), pointer :: CascadeCurrent
	integer cell
	integer reactionsFine, reactionCounter, reactionsTemp
	integer proc

	reactionCounter=0
	CascadeCurrent=>ActiveCascades
	do while(associated(cascadeCurrent))

		do cell=1,numCellsCascade
			reactionCurrent=>CascadeCurrent%reactionList(cell)
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
		CascadeCurrent=>CascadeCurrent%next
	end do

	call MPI_ALLREDUCE(reactionCounter,reactionsFine, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

end subroutine
