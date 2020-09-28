!***************************************************************************************
!>Subroutine initialize vacancy or SIA defect.
!***************************************************************************************
subroutine initializeVIdefect()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	use mod_randdp
	implicit none
	include 'mpif.h'

	integer :: cell, i, maxNumTemp
	double precision :: totalAtoms, rtemp, r1

	!Total atoms in the whole system
	totalAtoms = systemVol/atomSize
	!Number of Cu atoms in every mesh
	numCuCell = anint(CuContent*totalAtoms /dble(numTotal))

	do i=1, numSingleForm
		if(FormSingle(i)%defectType(1)==0 .AND. FormSingle(i)%defectType(2)==1 .AND. &
				FormSingle(i)%defectType(3)==0 .AND. FormSingle(i)%defectType(4)==0) then	!V
			ceqV = dexp(-FormSingle(i)%Ef / (kboltzmann*temperature))
		else if(FormSingle(i)%defectType(1)==0 .AND. FormSingle(i)%defectType(2)==0 .AND. &
				FormSingle(i)%defectType(3)==1 .AND. FormSingle(i)%defectType(4)==0) then	!SIA
			ceqI = dexp(-FormSingle(i)%Ef / (kboltzmann*temperature))
		end if
	end do

	if(numVac==0) then
		initialNumV = nint(ceqV * totalAtoms)
	else
		initialNumV = numVac
		firr=dble(numVac)/systemVol*atomSize
	end if
	if(numInt==0) then
		initialNumI = nint(ceqI * totalAtoms)
	else
		initialNumI = numInt
	end if

	if(initialNumV >= initialNumI) then
		maxNumTemp=initialNumV
	else
		maxNumTemp=initialNumI
	end if

	concV = ceqV
	concI = ceqI

	allocate(listVI(maxNumTemp,2))

	!V
	rtemp = 0d0
	if(myProc%taskid==MASTER) then
		outer1: do i=1, initialNumV
			r1 = dprand()

			inter1: do cell=1, numTotal
				rtemp = rtemp + dble(cell)/dble(numTotal)
				if(r1 <= rtemp) then
					!	tempID = cell-1
					!	x = mod(tempID,numx) +1
					!	tempID = tempID /numx
					!	y = mod(tempID,numy)+1
					!	z = tempID / numz +1
					!coordinate
					!	VcoordinateList(i,1) = meshLength*(x-1)+ meshLength/2d0
					!	VcoordinateList(i,2) = meshLength*(y-1)+ meshLength/2d0
					!	VcoordinateList(i,3) = meshLength*(z-1)+ meshLength/2d0
					!globalCell
					!	VgCellList(i)=cell
					listVI(i,1)=cell

					rtemp = 0d0
					exit inter1
				end if
			end do inter1
		end do outer1
	end if

	!SIA
	rtemp = 0d0
	if(myProc%taskid==MASTER) then
		outer2: do i=1, initialNumI
			r1 = dprand()

			inter2: do cell=1, numTotal
				rtemp = rtemp + dble(cell)/dble(numTotal)
				if(r1 <= rtemp) then

					!globalCell
					!	IgCellList(i)=cell
					listVI(i,2)=cell

					rtemp = 0d0
					exit inter2
				end if
			end do inter2
		end do outer2
	end if

	if(initialNumV > 0 .OR. initialNumI > 0) then
		call MPI_BCAST(listVI, maxNumTemp*2, MPI_INTEGER, MASTER, comm,ierr)
	end if

end subroutine

!*****************************************************************************************
!>Subroutine
!*****************************************************************************************
subroutine initializeOneCascade()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	use mod_randdp
	implicit none
	include 'mpif.h'

	double precision :: rtemp, r1
	integer :: cell

	rtemp = 0d0
	r1 = dprand()
	if(myProc%taskid==MASTER) then
		do cell=1, numTotal
			rtemp = rtemp + dble(cell)/dble(numTotal)
			if(r1 <= rtemp) then
				oneCascadeGCell = cell
				exit
			end if
		end do
	end if
	call MPI_BCAST(oneCascadeGCell, 1, MPI_INTEGER, MASTER, comm,ierr)

end  subroutine

!*****************************************************************************************
!>Subroutine initializeDefectList(): initializes defect list for each volume element
!Begins with a defect with type 0 0 0 0 and num 0. Note that numCells is needed in this subroutine
!*****************************************************************************************
subroutine initializeDefectList()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	integer :: cell, i, j
	type(defect), pointer :: defectCurrent, defectPrev

	nullify(defectCurrent)
	do cell=1,numCells
		allocate(defectList(cell)%defectType(SPECIES))
		do i=1,SPECIES
			defectList(cell)%defectType(i)=0
		end do
		defectList(cell)%num=0
		defectList(cell)%cellNumber=cell
		nullify(defectList(cell)%next)

		defectCurrent=>defectList(cell)

		!SIA_1
		if(initialNumI > 0) then
			do i=1, initialNumI
				if(listVI(i,2)==myMesh(cell)%globalCell) then
					if(defectCurrent%defectType(3)==1) then
						defectCurrent%num = defectCurrent%num +1
					else
						allocate(defectCurrent%next)
						defectCurrent=>defectCurrent%next
						allocate(defectCurrent%defectType(SPECIES))
						do j=1, SPECIES
							defectCurrent%defectType(j) = 0
						end do
						defectCurrent%defectType(3) = 1
						defectCurrent%num=1
						defectCurrent%cellNumber=cell
					end if

				end if
			end do
			nullify(defectCurrent%next)
		end if

		!V_1
		if(initialNumV > 0) then
			do i=1, initialNumV
				if(listVI(i,1)==myMesh(cell)%globalCell) then
					if(defectCurrent%defectType(2)==1) then
						defectCurrent%num = defectCurrent%num +1
					else
						allocate(defectCurrent%next)
						defectCurrent=>defectCurrent%next
						allocate(defectCurrent%defectType(SPECIES))
						do j=1, SPECIES
							defectCurrent%defectType(j) = 0
						end do
						defectCurrent%defectType(2) = 1
						defectCurrent%num=1
						defectCurrent%cellNumber=cell
					end if
				end if
			end do
			nullify(defectCurrent%next)
		end if

		!Cu_1
		if(CuContent > 0d0) then
			allocate(defectCurrent%next)
			defectCurrent=>defectCurrent%next
			allocate(defectCurrent%defectType(SPECIES))
			do i=1, SPECIES
				defectCurrent%defectType(i) = 0
			end do
			defectCurrent%defectType(1) = 1
			defectCurrent%num=numCuCell
			defectCurrent%cellNumber=cell
			nullify(defectCurrent%next)
		end if
	end do
end subroutine

!*****************************************************************************************
!>Subroutine initializeRandomSeeds()
!Uses the computer clock to initialize the random seed of the master processor.
!Then generates several integers in the master processor and uses them to initialize the
!random number seed of the other processors.
!*****************************************************************************************
subroutine initializeRandomSeeds()
	use mod_constants
	use mod_globalVariables
	use mod_structures
	use mod_randdp
	implicit none
	include 'mpif.h'

	integer :: randseed, i, irand
	integer :: status(MPI_STATUS_SIZE)

	if(myProc%taskid == MASTER) then
		call system_clock(Count=randseed)	!return randseed (integer, unit:ms)
		call sdprnd(randseed)
		do i=1,myProc%numtasks-1
			randseed=irand(randseed)
			call mpi_send(randseed, 1, MPI_INTEGER, i, 1000,comm, ierr)
		end do
	else
		call mpi_recv(randseed, 1, MPI_INTEGER, MASTER, 1000, comm, status, ierr)
		call sdprnd(randseed)
	end if

end subroutine

!*****************************************************************************************
!>Subroutine initializeTotalRate()
!finds the total reaction rate in the local processor once reaction lists have been initialized (including implantation reactions)
!*****************************************************************************************
subroutine initializeTotalRate()
	use mod_globalVariables
	use mod_structures
	implicit none

	type(Reaction), pointer :: reactionCurrent
	double precision :: rate, rateCell
	integer :: cell

	rate=0d0
	do cell=1,numCells
		rateCell=0d0
		reactionCurrent=>reactionList(cell)
		do while(associated(reactionCurrent))
			rate=rate+reactionCurrent%reactionRate
			rateCell=rateCell+reactionCurrent%reactionRate
			reactionCurrent=>reactionCurrent%next
		end do
		totalRateVol(cell)=rateCell	!Total reaction rate in this volume element
	end do
	totalRate=rate			!Total reaction rate in the entire processor

end subroutine

!*****************************************************************************************
!>Subroutine initializeReactionList()
!creates a new reaction list for each volume element and initializes implantation reactions (with rates)
!*****************************************************************************************
subroutine initializeReactionList()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	use mod_reactionrates
	implicit none

	integer :: cell, i, j, reac, dir
	type(reaction), pointer :: reactionCurrent
	type(defect), pointer :: defectCurrent
	type(defect), pointer :: defectCurrentTempV
	type(defect), pointer :: defectCurrentTempI

	do cell=1,numCells
		if(totalDPA > 0d0 .AND. dpaRate > 0d0) then

			if(implantType=='FrenkelPair') then
				reactionList(cell)%numReactants=0
				reactionList(cell)%numProducts=2
				allocate(reactionList(cell)%products(SPECIES,reactionList(cell)%numProducts))
				allocate(reactionList(cell)%cellNumber(reactionList(cell)%numProducts))
				allocate(reactionList(cell)%taskid(reactionList(cell)%numProducts))

				!search ImplantList for Frenkel Pair reactions
				do reac=1,numImplantReac
					if(ImplantReactions(reac)%numReactants==0 .AND. &
							ImplantReactions(reac)%numProducts==2) then	!we have found Frenkel pair implantation
						exit	!leave ImplantReactions(reac) pointed at this element of ImplantReactions
					end if
				end do

				do i=1,ImplantReactions(reac)%numProducts
					reactionList(cell)%products(:,i)=ImplantReactions(reac)%products(:,i)
					reactionList(cell)%cellNumber(i)=cell
					reactionList(cell)%taskid(i)=myMesh(cell)%proc
				end do
				reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(reac))
				nullify(reactionList(cell)%next)

			else if(implantType=='Cascade') then
				reactionList(cell)%numReactants=-10 	!This is a flag that the reaction is a cascade (not carried out like a standard reaction)
				reactionList(cell)%numProducts=0		!The products are taken care of separately from the standard defect update procedure
				allocate(reactionList(cell)%cellNumber(1))
				allocate(reactionList(cell)%taskid(1))
				reactionList(cell)%cellNumber(1)=cell
				reactionList(cell)%taskid(1)=myMesh(cell)%proc

				do reac=1,numImplantReac
					if(ImplantReactions(reac)%numReactants==-10 .AND. &
							ImplantReactions(reac)%numProducts==0) then	!we have found cascade implantation
						exit
					end if
				end do

				if(implantScheme=='MonteCarlo') then
					reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(reac))
				else if(implantScheme=='explicit') then
					reactionList(cell)%reactionRate=0d0
				end if
				nullify(reactionList(cell)%next)
			end if
		else	!!no implantation

			reactionList(cell)%numReactants=0
			reactionList(cell)%numProducts=0
			reactionList(cell)%reactionRate=0d0
			allocate(reactionList(cell)%cellNumber(1))
			allocate(reactionList(cell)%taskid(1))
			reactionList(cell)%cellNumber(1)=cell
			reactionList(cell)%taskid(1)=myMesh(cell)%proc
			nullify(reactionList(cell)%next)
		end if
		reactionCurrent=>reactionList(cell)
		defectCurrent=>defectList(cell)

		!Initialize possible reactions of free Cu
		if(CuContent > 0d0) then

			!*******************************************************
			!clustering: Cu+Cu->2Cu
			!*******************************************************
			do reac=1,numClusterReac
				if(ClusterReactions(reac)%reactants(1,1)==1 .AND. ClusterReactions(reac)%reactants(2,1)==0 .AND. &
						ClusterReactions(reac)%reactants(3,1)==0 .AND. ClusterReactions(reac)%reactants(4,1)==0 .AND. &
						ClusterReactions(reac)%reactants(1,2)==1 .AND. ClusterReactions(reac)%reactants(2,2)==0 .AND. &
						ClusterReactions(reac)%reactants(3,2)==0 .AND. ClusterReactions(reac)%reactants(4,2)==0) then
					exit
				end if
			end do

			if(reac <= numClusterReac) then
				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next
				reactionCurrent%numReactants=2
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))
				do j=1, SPECIES
					reactionCurrent%reactants(j,1)=ClusterReactions(reac)%reactants(j,1)
					reactionCurrent%reactants(j,2)=ClusterReactions(reac)%reactants(j,2)
					reactionCurrent%products(j,1)=ClusterReactions(reac)%products(j,1)
				end do
				do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do
				reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
						reactionCurrent%reactants(:,2), cell, ClusterReactions(reac))
				nullify(reactionCurrent%next)
			end if
		end if

		!**********************************************
		!initialize reactions about V
		!*********************************************
		nullify(defectCurrentTempV)
		do while(associated(defectCurrent))
			if(defectCurrent%defectType(1)==0 &
					.AND. defectCurrent%defectType(2)==1 &
					.AND. defectCurrent%defectType(3)==0 &
					.AND. defectCurrent%defectType(4)==0) then
				defectCurrentTempV=>defectCurrent
				exit
			else
				defectCurrent=>defectCurrent%next
			end if
		end do

		if(associated(defectCurrent)) then

			if(CuContent > 0d0) then
				!*******************************************************
				!clustering: Cu+V->CuV
				!*******************************************************
				!search ClusterList for Cu+V->CuV
				do reac=1,numClusterReac
					if(ClusterReactions(reac)%reactants(1,1)==1 .AND. ClusterReactions(reac)%reactants(2,1)==0 .AND. &
							ClusterReactions(reac)%reactants(3,1)==0 .AND. ClusterReactions(reac)%reactants(4,1)==0 .AND. &
							ClusterReactions(reac)%reactants(1,2)==0 .AND. ClusterReactions(reac)%reactants(2,2)==1 .AND. &
							ClusterReactions(reac)%reactants(3,2)==0 .AND. ClusterReactions(reac)%reactants(4,2)==0) then
						exit
					end if
				end do

				if(reac <= numClusterReac) then
					allocate(reactionCurrent%next)
					reactionCurrent=>reactionCurrent%next
					reactionCurrent%numReactants=2
					reactionCurrent%numProducts=1
					allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
					allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
					allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					do j=1, SPECIES
						reactionCurrent%reactants(j,1)=ClusterReactions(reac)%reactants(j,1)
						reactionCurrent%reactants(j,2)=ClusterReactions(reac)%reactants(j,2)
						reactionCurrent%products(j,1)=ClusterReactions(reac)%products(j,1)
					end do
					do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
						reactionCurrent%cellNumber(i)=cell
						reactionCurrent%taskid(i)=myMesh(cell)%proc
					end do
					reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
							reactionCurrent%reactants(:,2), cell, ClusterReactions(reac))
					nullify(reactionCurrent%next)
				end if
			end if

			!*******************************************************************
			!Diffusion: V->V
			!*******************************************************************
			!search DiffList for V->V
			do reac=1,numDiffReac
				if(DiffReactions(reac)%reactants(1,1)==0 .AND. DiffReactions(reac)%reactants(2,1)==1 .AND.&
						DiffReactions(reac)%reactants(3,1)==0 .AND.DiffReactions(reac)%reactants(4,1)==0 .AND.&
						DiffReactions(reac)%products(1,1)==0 .AND. DiffReactions(reac)%products(2,1)==1 .AND. &
						DiffReactions(reac)%products(3,1)==0 .AND. DiffReactions(reac)%products(4,1)==0) then
					exit
				end if
			end do

			if(reac <= numDiffReac) then

				do dir=1, 6
					allocate(reactionCurrent%next)
					reactionCurrent=>reactionCurrent%next

					reactionCurrent%numReactants=1
					reactionCurrent%numProducts=1
					allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
					allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
					allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
					allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))
					do j=1, SPECIES
						reactionCurrent%reactants(j,1)=DiffReactions(reac)%reactants(j,1)
						reactionCurrent%products(j,1)=DiffReactions(reac)%products(j,1)
					end do
					reactionCurrent%cellNumber(1)=cell
					reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(dir)
					reactionCurrent%taskid(1)=myMesh(cell)%proc
					reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(dir)
					reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(:,1),cell,myProc%taskid,&
							myMesh(cell)%neighbors(dir),myMesh(cell)%neighborProcs(dir),dir,DiffReactions(reac))
					nullify(reactionCurrent%next)
				end do
			end if

			!*******************************************************************
			!sinkRemoval: V->0
			!*******************************************************************
			!search sinkRemovalList for V->V
			do reac=1,numSinkReac
				if(SinkReactions(reac)%reactants(1,1)==0 .AND. &
						SinkReactions(reac)%reactants(2,1)==1 .AND. &
						SinkReactions(reac)%reactants(3,1)==0 .AND. &
						SinkReactions(reac)%reactants(4,1)==0) then
					exit
				end if
			end do

			if(reac <= numSinkReac) then

				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next
				reactionCurrent%numReactants=1
				reactionCurrent%numProducts=0
				allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

				do i=1,reactionCurrent%numReactants
					reactionCurrent%reactants(:,i)=SinkReactions(reac)%reactants(:,i)
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do
				reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(:,1), cell, &
						SinkReactions(reac))
				nullify(reactionCurrent%next)
			end if

			if(defectCurrent%num > 1) then
				!*******************************************************
				!clustering: V+V->2V
				!*******************************************************
				!search ClusterList for V+V->2V
				do reac=1,numClusterReac
					if(ClusterReactions(reac)%reactants(1,1)==0 .AND. ClusterReactions(reac)%reactants(2,1)==1 .AND. &
							ClusterReactions(reac)%reactants(3,1)==0 .AND. ClusterReactions(reac)%reactants(4,1)==0 .AND. &
							ClusterReactions(reac)%reactants(1,2)==0 .AND. ClusterReactions(reac)%reactants(2,2)==1 .AND. &
							ClusterReactions(reac)%reactants(3,2)==0 .AND. ClusterReactions(reac)%reactants(4,2)==0) then
						exit
					end if
				end do

				if(reac <= numClusterReac) then

					allocate(reactionCurrent%next)
					reactionCurrent=>reactionCurrent%next

					reactionCurrent%numReactants=2
					reactionCurrent%numProducts=1
					allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
					allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
					allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					do j=1, SPECIES
						reactionCurrent%reactants(j,1)=ClusterReactions(reac)%reactants(j,1)
						reactionCurrent%reactants(j,2)=ClusterReactions(reac)%reactants(j,2)
						reactionCurrent%products(j,1)=ClusterReactions(reac)%products(j,1)
					end do
					do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
						reactionCurrent%cellNumber(i)=cell
						reactionCurrent%taskid(i)=myMesh(cell)%proc
					end do
					reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
							reactionCurrent%reactants(:,2), cell, ClusterReactions(reac))
					nullify(reactionCurrent%next)
				end if
			end if
		end if

		!*******************************************************
		!initialize reactions about SIA
		do while(associated(defectCurrent))
			if(defectCurrent%defectType(1)==0 &
					.AND. defectCurrent%defectType(2)==0 &
					.AND. defectCurrent%defectType(3)==1 &
					.AND. defectCurrent%defectType(4)==0) then
				defectCurrentTempI=>defectCurrent
				exit
			else
				defectCurrent=>defectCurrent%next
			end if
		end do

		if(associated(defectCurrent)) then

			!*******************************************************************
			!Diffusion: SIA->SIA
			!*******************************************************************
			!search DiffList for SIA->SIA
			do reac=1,numDiffReac
				if(DiffReactions(reac)%reactants(1,1)==0 .AND. DiffReactions(reac)%reactants(2,1)==0 .AND. &
						DiffReactions(reac)%reactants(3,1)==1 .AND. DiffReactions(reac)%reactants(4,1)==0 .AND. &
						DiffReactions(reac)%products(1,1)==0 .AND. DiffReactions(reac)%products(2,1)==0 .AND. &
						DiffReactions(reac)%products(3,1)==1 .AND. DiffReactions(reac)%products(4,1)==0) then
					exit
				end if
			end do

			if(reac <= numDiffReac) then

				do dir=1, 6
					allocate(reactionCurrent%next)
					reactionCurrent=>reactionCurrent%next

					reactionCurrent%numReactants=1
					reactionCurrent%numProducts=1
					allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
					allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
					allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					do j=1, SPECIES
						reactionCurrent%reactants(j,1)=DiffReactions(reac)%reactants(j,1)
						reactionCurrent%products(j,1)=DiffReactions(reac)%products(j,1)
					end do
					reactionCurrent%cellNumber(1)=cell
					reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(dir)
					reactionCurrent%taskid(1)=myMesh(cell)%proc
					reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(dir)
					reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(:,1), cell, &
							myProc%taskid, myMesh(cell)%neighbors(dir), myMesh(cell)%neighborProcs(dir), dir, &
							DiffReactions(reac))
					nullify(reactionCurrent%next)
				end do
			end if

			!*******************************************************************
			!sinkRemoval: SIA->0
			!*******************************************************************
			!search sinkRemovalList for SIA->0
			do reac=1,numSinkReac
				if(SinkReactions(reac)%reactants(1,1)==0 .AND. &
						SinkReactions(reac)%reactants(2,1)==0 .AND. &
						SinkReactions(reac)%reactants(3,1)==1 .AND. &
						SinkReactions(reac)%reactants(4,1)==0) then
					exit
				end if
			end do

			if(reac <= numSinkReac) then

				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next
				reactionCurrent%numReactants=1
				reactionCurrent%numProducts=0
				allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

				do i=1,reactionCurrent%numReactants
					reactionCurrent%reactants(:,i)=SinkReactions(reac)%reactants(:,i)
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do
				reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(:,1), cell, &
						SinkReactions(reac))
				nullify(reactionCurrent%next)

			end if

			if(defectCurrent%num > 1) then
				!*******************************************************
				!clustering: SIA+SIA->2SIA
				!*******************************************************
				!search ClusterList for SIA+SIA->2SIA
				do reac=1,numClusterReac
					if(ClusterReactions(reac)%reactants(1,1)==0 .AND. ClusterReactions(reac)%reactants(2,1)==0 .AND. &
							ClusterReactions(reac)%reactants(3,1)==1 .AND. ClusterReactions(reac)%reactants(4,1)==0 .AND. &
							ClusterReactions(reac)%reactants(1,2)==0 .AND. ClusterReactions(reac)%reactants(2,2)==0 .AND. &
							ClusterReactions(reac)%reactants(3,2)==1 .AND. ClusterReactions(reac)%reactants(4,2)==0) then
						exit
					end if
				end do

				if(reac <= numClusterReac) then

					allocate(reactionCurrent%next)
					reactionCurrent=>reactionCurrent%next
					reactionCurrent%numReactants=2
					reactionCurrent%numProducts=1
					allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
					allocate(reactionCurrent%products(SPECIES,reactionCurrent%numProducts))
					allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					do j=1, SPECIES
						reactionCurrent%reactants(j,1)=ClusterReactions(reac)%reactants(j,1)
						reactionCurrent%reactants(j,2)=ClusterReactions(reac)%reactants(j,2)
						reactionCurrent%products(j,1)=ClusterReactions(reac)%products(j,1)
					end do
					if(pointDefectToggle=='yes') then
						if(reactionCurrent%products(3,1) > max3DInt) then
							reactionCurrent%products(4,1) = reactionCurrent%products(3,1)
							reactionCurrent%products(3,1) = 0
						end if
					end if

					do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
						reactionCurrent%cellNumber(i)=cell
						reactionCurrent%taskid(i)=myMesh(cell)%proc
					end do
					reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
							reactionCurrent%reactants(:,2), cell, ClusterReactions(reac))
					nullify(reactionCurrent%next)
				end if
			end if
		end if

		!***********************************
		!V+SIA->0
		!***********************************
		if(associated(defectCurrentTempV) .AND. associated(defectCurrentTempI)) then

			!*******************************************************
			!clustering: V+SIA->0
			!*******************************************************
			!search ClusterList for V+SIA->0
			do reac=1,numClusterReac
				if(ClusterReactions(reac)%reactants(1,1)==0 .AND. ClusterReactions(reac)%reactants(2,1)==1 .AND. &
						ClusterReactions(reac)%reactants(3,1)==0 .AND. ClusterReactions(reac)%reactants(4,1)==0 .AND. &
						ClusterReactions(reac)%reactants(1,2)==0 .AND. ClusterReactions(reac)%reactants(2,2)==0 .AND. &
						ClusterReactions(reac)%reactants(3,2)==1 .AND. ClusterReactions(reac)%reactants(4,2)==0) then
					exit
				end if
			end do

			if(reac <= numClusterReac) then

				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=2
				reactionCurrent%numProducts=0
				allocate(reactionCurrent%reactants(SPECIES,reactionCurrent%numReactants))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants))
				do j=1, SPECIES
					reactionCurrent%reactants(j,1)=ClusterReactions(reac)%reactants(j,1)
					reactionCurrent%reactants(j,2)=ClusterReactions(reac)%reactants(j,2)
				end do
				do i=1,reactionCurrent%numReactants
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do
				!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
				reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
						reactionCurrent%reactants(:,2), cell, ClusterReactions(reac))
				nullify(reactionCurrent%next)
			end if
		end if
	end do

end subroutine

!*****************************************************************************************
!> Subroutine initializeBoundaryDefectList()
!This subroutine initializes the defect lists within the boundary mesh.
!*****************************************************************************************
subroutine initializeBoundaryDefectList()
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	integer :: cell, dir, i, j, gCell, gNeighor
	type(defect), pointer :: defectCurrent
	integer, external :: findgNeighborPeriodic

	do cell=1,numCells
		do dir=1,6
			!neighbor cell is in different proc, not free surface (taskid=-1)
			if(myMesh(cell)%neighborProcs(dir) /= myProc%taskid .AND. &
					myMesh(cell)%neighborProcs(dir) == myProc%procNeighbor(dir) .AND. &
					myMesh(cell)%neighborProcs(dir) /= -1) then

				allocate(myBoundary(myMesh(cell)%neighbors(dir),dir)%defectList)
				defectCurrent=>myBoundary(myMesh(cell)%neighbors(dir),dir)%defectList
				allocate(defectCurrent%defectType(SPECIES))
				do i=1, SPECIES
					defectCurrent%defectType(i)=0
				end do
				defectCurrent%num=0
				defectCurrent%cellNumber=myMesh(cell)%neighbors(dir)
				nullify(myBoundary(myMesh(cell)%neighbors(dir),dir)%defectList%next)

				!SIA_1
				if(initialNumI > 0) then

					do i=1, initialNumI
						gCell=myMesh(cell)%globalCell
						gNeighor=findgNeighborPeriodic(gCell, dir)
						if(listVI(i,2)==gNeighor) then
							if(defectCurrent%defectType(3)==1) then
								defectCurrent%num=defectCurrent%num +1
							else
								allocate(defectCurrent%next)
								defectCurrent=>defectCurrent%next
								allocate(defectCurrent%defectType(SPECIES))
								do j=1, SPECIES
									defectCurrent%defectType(j)=0
								end do
								defectCurrent%defectType(3)=1
								defectCurrent%num=1
								defectCurrent%cellNumber=myMesh(cell)%neighbors(dir)
							end if
						end if
					end do
					nullify(defectCurrent%next)
				end if

				!V_1
				if(initialNumV > 0) then
					do i=1, initialNumV
						gCell=myMesh(cell)%globalCell
						gNeighor=findgNeighborPeriodic(gCell, dir)
						if(listVI(i,1)==gNeighor) then
							if(defectCurrent%defectType(2)==1) then
								defectCurrent%num=defectCurrent%num +1
							else
								allocate(defectCurrent%next)
								defectCurrent=>defectCurrent%next
								allocate(defectCurrent%defectType(SPECIES))
								do j=1, SPECIES
									defectCurrent%defectType(j)=0
								end do
								defectCurrent%defectType(2)=1
								defectCurrent%num=1
								defectCurrent%cellNumber=myMesh(cell)%neighbors(dir)
							end if
						end if
					end do
					nullify(defectCurrent%next)

				end if

				!Cu_1
				if(CuContent > 0d0) then
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(SPECIES))
					do i=1,SPECIES
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%defectType(1)=1
					defectCurrent%num=numCuCell
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir)
					nullify(defectCurrent%next)
				end if
			end if  !myMesh(cell)%neighborProcs(dir) /= myProc%taskid
		end do
	end do

end subroutine

!***************************************************************************************************
!>Subroutine: initializeFineMesh(CascadeCurrent)
!initializes defect and reaction lists in a newly created fine mesh
!***************************************************************************************************
subroutine initializeFineMesh(CascadeCurrent)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	use mod_randdp
	use mod_reactionrates
	implicit none

	type(cascade), pointer, intent(inout) :: CascadeCurrent
	type(defect), pointer :: defectCurrentCoarse, defectCurrentFine
	type(defect), pointer :: defectPrevCoarse, defectPrevFine
	integer :: j, n, k, num, cell, count
	double precision :: volumeRatio, r1, r2, lambda, rstore
	integer :: products(SPECIES)
	integer, external :: binomial, factorial

	interface
		subroutine findDefectInList(defectCurrent, defectPrev, products)
			use mod_constants
			use mod_structures
			use mod_globalVariables
			implicit none
			type(defect), pointer :: defectCurrent, defectPrev
			integer products(SPECIES)
		end subroutine
	end interface

	allocate(CascadeCurrent%localDefects(numCellsCascade))
	allocate(CascadeCurrent%reactionList(numCellsCascade))

	!For each cell, initialize the reaction list (no reaction) and defect list
	do cell=1,numCellsCascade
		allocate(CascadeCurrent%localDefects(cell)%defectType(SPECIES))
		nullify(CascadeCurrent%localDefects(cell)%next)
		do j=1,SPECIES
			CascadeCurrent%localDefects(cell)%defectType(j)=0
		end do
		CascadeCurrent%localDefects(cell)%num=0
		CascadeCurrent%localDefects(cell)%cellNumber=cell   !fineMeshID
		CascadeCurrent%reactionList(cell)%numReactants=0
		CascadeCurrent%reactionList(cell)%numProducts=0
		CascadeCurrent%reactionList(cell)%reactionRate=0d0
		nullify(CascadeCurrent%reactionList(cell)%next)
	end do

	defectCurrentCoarse=>defectList(CascadeCurrent%cellNumber)
	nullify(defectPrevCoarse)
	volumeRatio=CascadeElementVol*dble(numCellsCascade)/(myMesh(CascadeCurrent%cellNumber)%volume)

	do while(associated(defectCurrentCoarse))

		r1=dprand()
		r2=0d0
		n=defectCurrentCoarse%num

		if(n==0) then
			!do nothing
		else
			!*******************************************************************************************
			!finding the number of defects that are moved from coarse mesh to fine mesh requires a
			!binomial distribution. If the number of defects in the coarse mesh is sufficiently small,
			!a full binomial distribution can be used, but if it is too large, it must be approximated
			!by a Poisson distribution and then a Gaussian distrubition becuase large factorials cannot
			!be computed here.
			!*******************************************************************************************
			if(n <= 12) then    !binomial distribution
				!chooses number of defects that are added to fine mesh from coarse mesh of this type
				do k=0,n-1
					r2=r2+dble(binomial(n,k))*dble(volumeRatio**k)*dble((1d0-volumeRatio)**(n-k))
					if(r2 > r1) then
						exit
					end if
				end do
			else    !Poisson distribution, then Gaussian distrubition

				lambda=volumeratio*n
				rstore=0d0
				do k=0,n-1

					!***********************************************************************************
					!For k .LE. 12, a poisson distribution can be calculated. For k>12, factorial(k) is
					!too large and a Gaussian distribution must be assumed.
					!
					!This leads to some errors where k=12 is slightly more common than it should be as
					!a number of defects added from the
					!coarse mesh to the fine mesh, but this problem is almost impossible to overcome
					!considering the size of the numbers that the machine can handle.
					!***********************************************************************************
					if(k <= 12) then    !Poisson distribution
						r2=r2+lambda**(dble(k))*dexp(-lambda)/factorial(k)
						if(r2 > r1) then
							exit
						end if
					else    !Gaussian distrubition
						r2=r2+dexp(-(k-lambda)**(2d0)/(2d0*lambda))/(dsqrt(2*pi*lambda))
						if(rstore-r2 == 0d0) then
							exit
						else if(r2 > r1) then
							exit
						end if
						rstore=r2
					end if
				end do
			end if

			!************************************************************
			!k is the number of defects being deposited into the fine mesh
			!************************************************************
			if(k >  0) then
				do num=1,k

					!***************************************************************
					!Generate the cell within the fine mesh that the defect will be deposited into
					!***************************************************************
					r1=dprand()*dble(numCellsCascade)
					r2=0d0
					do cell=1,numCellsCascade
						r2=r2+dble(cell)
						if(r2 > r1) then
							exit
						end if
					end do

					!***************************************************************
					!Deposit the defect into the fine mesh
					!***************************************************************
					do j=1,SPECIES
						products(j)=defectCurrentCoarse%defectType(j)
					end do
					!write(*,*) 'inserting into fine mesh', (products(j), j=1,SPECIES), 'k', k
					nullify(defectPrevFine)
					defectCurrentFine=>CascadeCurrent%localDefects(cell)

					call findDefectInList(defectCurrentFine, defectPrevFine, products)
					if(associated(defectCurrentFine)) then !if we aren't at the end of the list

						count=0
						!Check to see if this defect already exists in the fine mesh list
						do j=1,SPECIES
							if(defectCurrentFine%defectType(j)==products(j)) then
								count=count+1
							end if
						end do

						if(count==SPECIES) then
							defectCurrentFine%num=defectCurrentFine%num+1
						else		!if the defect is to be inserted in the list
							if(.NOT. associated(defectPrevFine)) then
								write(*,*) 'error defectPrevFine not associated'
							end if
							nullify(defectPrevFine%next)
							allocate(defectPrevFine%next)
							nullify(defectPrevFine%next%next)
							defectPrevFine=>defectPrevFine%next
							allocate(defectPrevFine%defectType(SPECIES))
							defectPrevFine%cellNumber=cell
							defectPrevFine%num=1
							do j=1,SPECIES
								defectPrevFine%defectType(j)=products(j)
							end do
							defectPrevFine%next=>defectCurrentFine
						end if
					else 			!add a defect to the end of the list
						nullify(defectPrevFine%next)
						allocate(defectPrevFine%next)
						nullify(defectPrevFine%next%next)
						defectPrevFine=>defectPrevFine%next
						allocate(defectPrevFine%defectType(SPECIES))
						defectPrevFine%cellNumber=cell
						defectPrevFine%num=1
						do j=1,SPECIES
							defectPrevFine%defectType(j)=products(j)
						end do
					end if
				end do

				!***************************************************************
				!Remove the defect from the coarse mesh
				!***************************************************************
				if(k < n) then
					defectCurrentCoarse%num=defectCurrentCoarse%num-k !remove the defect from the system instead of the entire entry in the list
				else	!defects of this type are all deposited into fine meshes
					if(.NOT. associated(defectCurrentCoarse)) then
						write(*,*) 'Tried to delete defect that wasnt there fine mesh initialization'
					!if defectCurrentCoarse is in the middle of the list
					else if(associated(defectCurrentCoarse%next) .AND. associated(defectPrevCoarse)) then

						defectPrevCoarse%next=>defectCurrentCoarse%next !remove that defect type from the system
						deallocate(defectCurrentCoarse%defectType)
						deallocate(defectCurrentCoarse)
						defectCurrentCoarse=>defectPrevCoarse

					!if defectCurrentCoarse is at the end of the list
					else if(associated(defectPrevCoarse)) then
						deallocate(defectCurrentCoarse%defectType)
						deallocate(defectCurrentCoarse)
						defectCurrentCoarse=>defectPrevCoarse	!remove the last defect from the system
						nullify(defectPrevCoarse%next)

					!if defectCurrentCoarse is at the beginning of the list and there is one of them
					else if(associated(defectCurrentCoarse%next)) then !removing first defect from cell i
						defectCurrentCoarse%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
					end if
				end if
			end if
		end if
		defectPrevCoarse=>defectCurrentCoarse
		defectCurrentCoarse=>defectCurrentCoarse%next
	end do

end subroutine

!***********************************************************************
!> Subroutine annealInitialization()
!This subroutine carries out the following tasks in order to switch from damage introduction to annealing:
!***********************************************************************
subroutine annealInitialization()
	use mod_structures
	use mod_globalVariables
	implicit none

	integer :: cell
	type(reaction), pointer :: reactionCurrent
	type(cascade), pointer :: CascadeCurrent

	!Change temperature to annealTemp
	temperature=annealTemp
	do cell=1,numCells
		reactionCurrent=>reactionList(cell)
		totalRate=totalRate-reactionCurrent%reactionRate
		totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
		reactionCurrent%reactionRate=0d0
	end do

end subroutine
