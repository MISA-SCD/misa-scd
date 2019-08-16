! **************************************************************************************
!>Subroutine initialize vacancy or SIA defect.
!***************************************************************************************
subroutine initializeVIdefect()
use DerivedType
use mod_constants
use randdp
implicit none
include 'mpif.h'

integer cell, i,j
integer tempID,x,y,z
double precision totalAtoms
double precision rtemp, r1


totalAtoms = (myProc%globalCoord(2)-myProc%globalCoord(1))/lattice * &
		(myProc%globalCoord(4)-myProc%globalCoord(3))/lattice * &
		(myProc%globalCoord(6)-myProc%globalCoord(5))/lattice * 2

numCuCell = anint(CuContent*totalAtoms /dble(numTotal))

do i=1, numSingleForm(1)
	if(FormSingle(i,1)%defectType(1)==0 .AND. FormSingle(i,1)%defectType(2)==1 .AND. &
			FormSingle(i,1)%defectType(3)==0 .AND. FormSingle(i,1)%defectType(4)==0) then	!V
		initialCeqv = dexp(-FormSingle(i,1)%Ef / (kboltzmann*temperature))
	else if(FormSingle(i,1)%defectType(1)==0 .AND. FormSingle(i,1)%defectType(2)==0 .AND. &
			FormSingle(i,1)%defectType(3)==1 .AND. FormSingle(i,1)%defectType(4)==0) then	!SIA
		initialCeqi = dexp(-FormSingle(i,1)%Ef / (kboltzmann*temperature))
	end if
end do

initialTotalV = nint(initialCeqv * totalAtoms)
initialTotalI = nint(initialCeqi * totalAtoms)

Vconcent = initialCeqv
SIAconcent = initialCeqi

allocate(VgCellList(initialTotalV))
allocate(IgCellList(initialTotalI))

!V
rtemp = 0d0
if(myProc%taskid==MASTER) then
	outer1: do i=1, initialTotalV
		r1 = dprand()

		inter1: do cell=1, numTotal
			rtemp = rtemp + cell/numTotal
			if(r1 <= rtemp) then
!				tempID = cell-1
!				x = mod(tempID,numx) +1
!				tempID = tempID /numx
!				y = mod(tempID,numy)+1
!				z = tempID / numz +1
				!coordinate
!				VcoordinateList(i,1) = meshLength*(x-1)+ meshLength/2d0
!				VcoordinateList(i,2) = meshLength*(y-1)+ meshLength/2d0
!				VcoordinateList(i,3) = meshLength*(z-1)+ meshLength/2d0
				!globalCell
				VgCellList(i)=cell

				rtemp = 0d0
				exit inter1
			end if
		end do inter1
	end do outer1
end if
if(initialTotalV > 0) then
	call MPI_BCAST(VgCellList, initialTotalV, MPI_INTEGER, MASTER, comm,ierr)
end if

!SIA
rtemp = 0d0
if(myProc%taskid==MASTER) then
	outer2: do i=1, initialTotalI
		r1 = dprand()

		inter2: do cell=1, numTotal
			rtemp = rtemp + cell/numTotal
			if(r1 <= rtemp) then

				!globalCell
				IgCellList(i)=cell

				rtemp = 0d0
				exit inter2
			end if
		end do inter2
	end do outer2
end if

if(initialTotalI > 0) then
	call MPI_BCAST(IgCellList, initialTotalI, MPI_INTEGER, MASTER, comm,ierr)
end if

end subroutine

!*****************************************************************************************
!>Subroutine initialize defect list - initializes defect list for each volume element
!
!Begins with a defect with type 0 0 0 0 and num 0. Note that numCells is needed in this subroutine
!*****************************************************************************************

subroutine initializeDefectList()
use DerivedType
use mod_constants
implicit none

integer cell, i
type(defect), pointer :: defectCurrent, defectPrev

nullify(defectCurrent)

do cell=1,numCells
	allocate(defectList(cell)%defectType(numSpecies))
	do i=1,numSpecies
		defectList(cell)%defectType(i)=0
	end do
	defectList(cell)%num=0
	defectList(cell)%cellNumber=cell
	nullify(defectList(cell)%next)

	defectCurrent=>defectList(cell)

	!Cu_1
	if(CuContent > 0d0) then
		allocate(defectCurrent%next)
		defectCurrent=>defectCurrent%next
		allocate(defectCurrent%defectType(numSpecies))
		do i=1, numSpecies
			defectCurrent%defectType(i) = 0
		end do
		defectCurrent%defectType(1) = 1
		defectCurrent%num=numCuCell
		defectCurrent%cellNumber=cell
		nullify(defectCurrent%next)
	end if

	!V_1
	if(initialTotalV > 0) then
		allocate(defectCurrent%next)
		defectCurrent=>defectCurrent%next
		allocate(defectCurrent%defectType(numSpecies))
		do i=1, numSpecies
			defectCurrent%defectType(i) = 0
		end do
		defectCurrent%defectType(2) = 1
		defectCurrent%num=0
		defectCurrent%cellNumber=cell

		do i=1, initialTotalV
			if(VgCellList(i)==myMesh(cell)%globalCell) then
				defectCurrent%num = defectCurrent%num +1
			end if
		end do

		if(defectCurrent%num > 0) then
			nullify(defectCurrent%next)
		else
			deallocate(defectCurrent%defectType)
			deallocate(defectCurrent)
			nullify(defectCurrent)
		end if
	end if

	defectCurrent=>defectList(cell)
	do while(associated(defectCurrent%next))
		defectCurrent=>defectCurrent%next
	end do

	!SIA_1
	if(initialTotalI > 0) then
		allocate(defectCurrent%next)
		defectCurrent=>defectCurrent%next
		allocate(defectCurrent%defectType(numSpecies))
		do i=1, numSpecies
			defectCurrent%defectType(i) = 0
		end do
		defectCurrent%defectType(3) = 1
		defectCurrent%num=0
		defectCurrent%cellNumber=cell

		do i=1, initialTotalI
			if(IgCellList(i)==myMesh(cell)%globalCell) then
				defectCurrent%num = defectCurrent%num +1
			end if
		end do

		if(defectCurrent%num > 0) then
			nullify(defectCurrent%next)
		else
			deallocate(defectCurrent%defectType)
			deallocate(defectCurrent)
			nullify(defectCurrent)
		end if
	end if

end do

end subroutine

!*****************************************************************************************
!>Subroutine initialize random seeds - seeds random number generator differently for each processor
!!
!!Uses the computer clock to initialize the random seed of the master processor. Then 
!!generates several integers in the master processor and uses them to initialize the 
!!random number seed of the other processors.
!*****************************************************************************************

subroutine initializeRandomSeeds()
use mod_constants
use DerivedType
use randdp
implicit none
include 'mpif.h'

integer randseed, i, irand
integer status(MPI_STATUS_SIZE)

if(myProc%taskid == MASTER) then
	call system_clock(Count=randseed)	!return randseed (integer, unit:ms)
	call sdprnd(randseed)
!	write(*,*)
!	write(*,*) 'random number seed', randseed, 'processor', myProc%taskid
	do i=1,myProc%numtasks-1
		randseed=irand(randseed)
		call mpi_send(randseed, 1, MPI_INTEGER, i, 1000,comm, ierr)
	end do
else
	call mpi_recv(randseed, 1, MPI_INTEGER, MASTER, 1000, comm, status, ierr)
	call sdprnd(randseed)
	!write(*,*)
!	write(*,*) 'random number seed', randseed, 'processor', myProc%taskid
endif

end subroutine

!*****************************************************************************************
!>Subroutine initialize total rate - finds the total reaction rate in the local processor once
!!reaction lists have been initialized (including implantation reactions)
!!
!!This will also initialize the total rate in the case of a debug restart, in which 
!!many reactions are possible at the first step (besides just implantation reactions)

!!Return: totalRateVol(:) Total reaction rate in each volume element of the  local processor
!!		  totalRate		  Total rate of the local processor
!*****************************************************************************************

subroutine initializeTotalRate()
use mod_constants
use DerivedType
implicit none

type(Reaction), pointer :: reactionCurrent
double precision rate, rateCell, rateSingle
integer cell

rate=0d0
rateSingle=0d0

do cell=1,numCells
	rateCell=0d0
	reactionCurrent=>reactionList(cell)
	
	do while(associated(reactionCurrent))
		rate=rate+reactionCurrent%reactionRate
		rateCell=rateCell+reactionCurrent%reactionRate
		reactionCurrent=>reactionCurrent%next
	end do
	
	totalRateVol(cell)=rateCell	!Total reaction rate in this volume element
	
	if(singleElemKMC=='yes') then
		if(rateCell > rateSingle) then
			rateSingle=rateCell
		end if
	end if
	
end do

if(singleElemKMC=='yes') then
	totalRate=rateSingle	!max reaction rate among all volume elements in this processor
else
	totalRate=rate			!Total reaction rate in the entire processor
endif

end subroutine

!*****************************************************************************************
!>Subroutine initialize reaction list - creates a new reaction list for each volume element
!! and initializes implantation reactions (with rates)
!!
!!First reaction in list is either Frenkel pair implantation or cascade implantation, second
!!reaction in the list is Cu clustering, the .
!*****************************************************************************************

subroutine initializeReactionList()
use DerivedType
use mod_constants
use ReactionRates
implicit none

integer cell, i, j, reac, matNum, count
integer dir
type(reaction), pointer :: reactionCurrent
type(defect), pointer :: defectCurrent
type(defect), pointer :: defectCurrentTempV
type(defect), pointer :: defectCurrentTempI

do cell=1,numCells
	!Initialize reaction rates with Frenkel pair implantation. The rate should be given by a function
	!findReactionRate which has as parameters passed to it the cell that this rate is occuring in and
	!the parameters of the reaction rate as read in from the file (type reactionParameters)
	
	!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
	!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
	if(numMaterials==1) then
		matNum=1
	else
		matNum=myMesh(cell)%material
	end if

	if(totalDPA > 0d0 .AND. DPARate > 0d0) then

		if(implantType=='FrenkelPair') then

			!Frenkel pair implantaiton reaction: first in the list
		
			reactionList(cell)%numReactants=0
			reactionList(cell)%numProducts=2
			allocate(reactionList(cell)%products(numSpecies,reactionList(cell)%numProducts))
			allocate(reactionList(cell)%cellNumber(reactionList(cell)%numProducts))
			allocate(reactionList(cell)%taskid(reactionList(cell)%numProducts))
		
			!search ImplantList for Frenkel Pair reactions
			do reac=1,numImplantReac(matNum)
				if(ImplantReactions(reac,matNum)%numReactants==0 .AND. &
						ImplantReactions(reac,matNum)%numProducts==2) then	!we have found Frenkel pair implantation
					exit	!leave ImplantReactions(reac) pointed at this element of ImplantReactions
				end if
			end do
		
			do i=1,ImplantReactions(reac,matNum)%numProducts

				reactionList(cell)%products(:,i)=ImplantReactions(reac,matNum)%products(:,i)

				reactionList(cell)%cellNumber(i)=cell
				reactionList(cell)%taskid(i)=myMesh(cell)%proc
			end do
		
			!Find reaction rate for Frenkel pair implantation using ImplantReactions(reac), which is input info from file.
			reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(reac,matNum))
			nullify(reactionList(cell)%next)

		else if(implantType=='Cascade') then
			!initialize reaction rates with cascade implantation. The rate should be given by a function
			!findReactionRate which has as parameters passed to it the cell that this rate is occuring in
			!and the parameters of the reaction rate as read in from the file (type reactionParameters)

			reactionList(cell)%numReactants=-10 	!This is a flag that the reaction is a cascade (not carried out like a standard reaction)
			reactionList(cell)%numProducts=0		!The products are taken care of separately from the standard defect update procedure
			allocate(reactionList(cell)%cellNumber(1))
			allocate(reactionList(cell)%taskid(1))

			reactionList(cell)%cellNumber(1)=cell
			reactionList(cell)%taskid(1)=myMesh(cell)%proc

			!search ImplantList for cascade reactions
			do reac=1,numImplantReac(matNum)
				if(ImplantReactions(reac,matNum)%numReactants==-10 .AND. &
						ImplantReactions(reac,matNum)%numProducts==0) then	!we have found cascade implantation
					exit
				end if
			end do

			!If implant scheme is Monte Carlo, then this reaction is chosen as part of the regular
			!Monte Carlo algorithm. Otherwise, cascades are added explicitly and this rate is set to 0
			!so that cascade implantation is no longer chosen as part of the MC algorithm.

			if(implantScheme=='MonteCarlo') then
				reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(reac,matNum))
			else if(implantScheme=='explicit') then
				reactionList(cell)%reactionRate=0d0
			endif

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
		!search ClusterList for Cu+Cu->2Cu
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(reac,matNum)%reactants(1,1)==1 .AND. &
					ClusterReactions(reac,matNum)%reactants(2,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(3,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(4,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(1,2)==1 .AND. &
					ClusterReactions(reac,matNum)%reactants(2,2)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(3,2)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(4,2)==0) then
				exit
			end if
		end do

		if(reac <= numClusterReac(matNum)) then
			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
			allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

			reactionCurrent%reactants(:,1)=ClusterReactions(reac,matNum)%reactants(:,1)
			reactionCurrent%reactants(:,2)=ClusterReactions(reac,matNum)%reactants(:,2)

			reactionCurrent%products(:,1)=ClusterReactions(reac,matNum)%products(:,1)

			do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			end do

			!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
			reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
					reactionCurrent%reactants(:,2), cell, ClusterReactions(reac,matNum))

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
			do reac=1,numClusterReac(matNum)
				if(ClusterReactions(reac,matNum)%reactants(1,1)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(1,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,2)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,2)==0) then
					exit
				end if
			end do

			if(reac <= numClusterReac(matNum)) then
				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=2
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))

				reactionCurrent%reactants(:,1)=ClusterReactions(reac,matNum)%reactants(:,1)
				reactionCurrent%reactants(:,2)=ClusterReactions(reac,matNum)%reactants(:,2)

				reactionCurrent%products(:,1)=ClusterReactions(reac,matNum)%products(:,1)

				do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do

				!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
				reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
						reactionCurrent%reactants(:,2), cell, ClusterReactions(reac,matNum))
				nullify(reactionCurrent%next)
			end if

		end if

		!*******************************************************************
		!Diffusion: V->V
		!*******************************************************************
		!search DiffList for V->V
		do reac=1,numDiffReac(matNum)
			if(DiffReactions(reac,matNum)%reactants(1,1)==0 .AND. &
					DiffReactions(reac,matNum)%reactants(2,1)==1 .AND.&
					DiffReactions(reac,matNum)%reactants(3,1)==0 .AND.&
					DiffReactions(reac,matNum)%reactants(4,1)==0 .AND.&
					DiffReactions(reac,matNum)%products(1,1)==0 .AND. &
					DiffReactions(reac,matNum)%products(2,1)==1 .AND. &
					DiffReactions(reac,matNum)%products(3,1)==0 .AND. &
					DiffReactions(reac,matNum)%products(4,1)==0) then
				exit
			end if
		end do

		if(reac <= numDiffReac(matNum)) then

			do dir=1, 6
				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=1
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))


				reactionCurrent%reactants(:,1)=DiffReactions(reac,matNum)%reactants(:,1)
				reactionCurrent%products(:,1)=DiffReactions(reac,matNum)%products(:,1)

				reactionCurrent%cellNumber(1)=cell
				reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(1,dir)
				reactionCurrent%taskid(1)=myMesh(cell)%proc
				reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(1,dir)

				reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(:,1), cell, &
						myProc%taskid, myMesh(cell)%neighbors(1,dir), myMesh(cell)%neighborProcs(1,dir), dir, &
						DiffReactions(reac,matNum))
				nullify(reactionCurrent%next)
			end do

		end if

		!*******************************************************************
		!sinkRemoval: V->0
		!*******************************************************************
		!search sinkRemovalList for V->V
		do reac=1,numSinkReac(matNum)
			if(SinkReactions(reac,matNum)%reactants(1,1)==0 .AND. &
					SinkReactions(reac,matNum)%reactants(2,1)==1 .AND. &
					SinkReactions(reac,matNum)%reactants(3,1)==0 .AND. &
					SinkReactions(reac,matNum)%reactants(4,1)==0) then
				exit
			end if
		end do

		if(reac <= numSinkReac(matNum)) then

			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=0
			allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

			do i=1,reactionCurrent%numReactants

				reactionCurrent%reactants(:,i)=SinkReactions(reac,matNum)%reactants(:,i)

				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			end do

			reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(:,1), cell, &
					SinkReactions(reac,matNum))
			nullify(reactionCurrent%next)

		end if

		if(defectCurrent%num > 1) then
			!*******************************************************
			!clustering: V+V->2V
			!*******************************************************
			!search ClusterList for V+V->2V
			do reac=1,numClusterReac(matNum)
				if(ClusterReactions(reac,matNum)%reactants(1,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,1)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(1,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,2)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,2)==0) then
					exit
				end if
			end do

			if(reac <= numClusterReac(matNum)) then

				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=2
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))



				reactionCurrent%reactants(:,1)=ClusterReactions(reac,matNum)%reactants(:,1)
				reactionCurrent%reactants(:,2)=ClusterReactions(reac,matNum)%reactants(:,2)
				reactionCurrent%products(:,1)=ClusterReactions(reac,matNum)%products(:,1)

				do i=1,reactionCurrent%numReactants+reactionCurrent%numProducts
					reactionCurrent%cellNumber(i)=cell
					reactionCurrent%taskid(i)=myMesh(cell)%proc
				end do
				!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
				reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
						reactionCurrent%reactants(:,2), cell, ClusterReactions(reac,matNum))
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
		do reac=1,numDiffReac(matNum)
			if(DiffReactions(reac,matNum)%reactants(1,1)==0 .AND. &
					DiffReactions(reac,matNum)%reactants(2,1)==0 .AND. &
					DiffReactions(reac,matNum)%reactants(3,1)==1 .AND. &
					DiffReactions(reac,matNum)%reactants(4,1)==0 .AND. &
					DiffReactions(reac,matNum)%products(1,1)==0 .AND. &
					DiffReactions(reac,matNum)%products(2,1)==0 .AND. &
					DiffReactions(reac,matNum)%products(3,1)==1 .AND. &
					DiffReactions(reac,matNum)%products(4,1)==0) then
				exit
			end if
		end do

		if(reac <= numDiffReac(matNum)) then

			do dir=1, 6
				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=1
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))

				reactionCurrent%reactants(:,1)=DiffReactions(reac,matNum)%reactants(:,1)
				reactionCurrent%products(:,1)=DiffReactions(reac,matNum)%products(:,1)

				reactionCurrent%cellNumber(1)=cell
				reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(1,dir)
				reactionCurrent%taskid(1)=myMesh(cell)%proc
				reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(1,dir)

				reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(:,1), cell, &
						myProc%taskid, myMesh(cell)%neighbors(1,dir), myMesh(cell)%neighborProcs(1,dir), dir, &
						DiffReactions(reac,matNum))
				nullify(reactionCurrent%next)

			end do

		end if

		!*******************************************************************
		!sinkRemoval: SIA->0
		!*******************************************************************
		!search sinkRemovalList for SIA->0
		do reac=1,numSinkReac(matNum)
			if(SinkReactions(reac,matNum)%reactants(1,1)==0 .AND. &
					SinkReactions(reac,matNum)%reactants(2,1)==0 .AND. &
					SinkReactions(reac,matNum)%reactants(3,1)==1 .AND. &
					SinkReactions(reac,matNum)%reactants(4,1)==0) then
				exit
			end if
		end do

		if(reac <= numSinkReac(matNum)) then

			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=0
			allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

			do i=1,reactionCurrent%numReactants

				reactionCurrent%reactants(:,i)=SinkReactions(reac,matNum)%reactants(:,i)

				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			end do
			reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(:,1), cell, &
					SinkReactions(reac,matNum))
			nullify(reactionCurrent%next)

		end if

		if(defectCurrent%num > 1) then
			!*******************************************************
			!clustering: SIA+SIA->2SIA
			!*******************************************************
			!search ClusterList for SIA+SIA->2SIA
			do reac=1,numClusterReac(matNum)
				if(ClusterReactions(reac,matNum)%reactants(1,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,1)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,1)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(1,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(2,2)==0 .AND. &
						ClusterReactions(reac,matNum)%reactants(3,2)==1 .AND. &
						ClusterReactions(reac,matNum)%reactants(4,2)==0) then
					exit
				end if
			end do

			if(reac <= numClusterReac(matNum)) then

				allocate(reactionCurrent%next)
				reactionCurrent=>reactionCurrent%next

				reactionCurrent%numReactants=2
				reactionCurrent%numProducts=1
				allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
				allocate(reactionCurrent%products(numSpecies,reactionCurrent%numProducts))
				allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
				allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))


				reactionCurrent%reactants(:,1)=ClusterReactions(reac,matNum)%reactants(:,1)
				reactionCurrent%reactants(:,2)=ClusterReactions(reac,matNum)%reactants(:,2)
				reactionCurrent%products(:,1)=ClusterReactions(reac,matNum)%products(:,1)

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
				!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
				reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
						reactionCurrent%reactants(:,2), cell, ClusterReactions(reac,matNum))
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
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(reac,matNum)%reactants(1,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(2,1)==1 .AND. &
					ClusterReactions(reac,matNum)%reactants(3,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(4,1)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(1,2)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(2,2)==0 .AND. &
					ClusterReactions(reac,matNum)%reactants(3,2)==1 .AND. &
					ClusterReactions(reac,matNum)%reactants(4,2)==0) then
				exit
			end if
		end do

		if(reac <= numClusterReac(matNum)) then

			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=0
			allocate(reactionCurrent%reactants(numSpecies,reactionCurrent%numReactants))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

			reactionCurrent%reactants(:,1)=ClusterReactions(reac,matNum)%reactants(:,1)
			reactionCurrent%reactants(:,2)=ClusterReactions(reac,matNum)%reactants(:,2)

			do i=1,reactionCurrent%numReactants
				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			end do
			!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
			reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(:,1), &
					reactionCurrent%reactants(:,2), cell, ClusterReactions(reac,matNum))
			nullify(reactionCurrent%next)

		end if

	end if

end do

end subroutine

!***************************************************************************************************
!
!> Subroutine initializeDebugRestart() - initializes defect and reaction lists if we are restarting from debug file
!!
!! This subroutine is used to populate the mesh with defects and initialize the DPA at a non-zero value
!! for debugging (restart from a saved point).
!!
!! This subroutine reads in from an input file. The mesh must match the mesh in the actual simulation.
!! The number of processors must also match from the reset file and the current simulation.
!
!***************************************************************************************************

subroutine initializeDebugRestart()
use mod_constants
use DerivedType
implicit none

character*20 char

logical flag

integer procID, numProcs, numDefectTypes
integer defectTypeReset(numSpecies), defectNumReset
integer i, j, k

double precision cellCoord(3)

type(defect), pointer :: defectCurrent, defectPrev

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

flag=.FALSE.

!Read in the name of the debug restart file
if(debugToggle == 'yes') then
	
	open(87, file=restartFileName,action='read', status='old')
	
	!Step 1: read in header information from restart file
	
	do while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='numProcs') then
			read(87,*) numProcs
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	if(numProcs /= myProc%numTasks) then
		write(*,*) 'Error restart file incorrect number of processors'
		if(myProc%taskid==MASTER) read(*,*)
	end if
	
	do while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='numImplantEvents') then
			read(87,*) numImplantEventsReset
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

	do while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='elapsedTime') then
			read(87,*) elapsedTimeReset
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	!Step 2: skip to part of input file that matches this processor
	
	do while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='processor') then
			read(87,*) procID
			if(procID==myProc%taskid) then
				flag=.TRUE.
			end if
		end if
	end do
	flag=.FALSE.
	
	!Step 3: Read in defect information from file to defect list (coarse mesh),
	!noting that the coordinates in the input file must match the coordinates
	!of the coarse mesh element.
	
	do i=1,numCells
		
		do while(flag .eqv. .FALSE.)
			read(87,*) char
			if(char=='coordinates') then
				read(87,*) (cellCoord(j), j=1,3)
				flag=.TRUE.
			end if
		end do
		flag=.FALSE.
		
		do j=1,3
			if(cellCoord(j) /= myMesh(i)%coordinates(j)) then
				write(*,*) 'Error cell coordinates do not match in reset file'
				write(*,*) 'CellCoord', cellCoord(j), 'myMeshCoord', myMesh(i)%coordinates(j)
				if(myProc%taskid==MASTER) read(*,*)
			end if
		end do
		
		do while(flag .eqv. .FALSE.)
			read(87,*) char
			if(char=='numDefectTypes') then
				read(87,*) numDefectTypes
				flag=.TRUE.
			end if
		end do
		flag=.FALSE.
		
		do j=1,numDefectTypes
			read(87,*) (defectTypeReset(k), k=1,numSpecies)
			read(87,*) defectNumReset
			
			defectCurrent=>defectList(i)
			nullify(defectPrev)
			
			call findDefectInList(defectCurrent, defectPrev, defectTypeReset)
			
			!insert defect in list. Since we are starting from an empty list and adding
			!defects in order, we should always have defectPrev pointing at the last 
			!element of the list and defectCurrent pointing at nothing.
			
			if(associated(defectCurrent)) then
				write(*,*) 'error defectCurrent associated in initializeDebugRestart'
			else
				!Create new defect in list
				allocate(defectCurrent)
				allocate(defectCurrent%defectType(numSpecies))
				nullify(defectCurrent%next)
				defectPrev%next=>defectCurrent
				
				do k=1,numSpecies
					defectCurrent%defectType(k)=defectTypeReset(k)
				end do
				defectCurrent%num=defectNumReset
				defectCurrent%cellNumber=i
				
			end if
		end do
		
	end do
	
	close(87)
end if
flag=.FALSE.

end subroutine

!*****************************************************************************************
!> Subroutine initialize boundary defect list - initializes defect lists in the boundary mesh
!!
!!This subroutine initializes the defect lists within the boundary mesh. 
!!
!!NOTE: it is currently not set up to handle situaions where myMesh(cell)%neighborProcs(dir,k) .NE. myProc%procNeighbor(dir)
!!(AKA when the element to the left is a different proc than the proc to the left, due to uneven meshing)
!!This may cause errors when the non-uniform mesh is used.
!*****************************************************************************************

subroutine initializeBoundaryDefectList()
use DerivedType
use mod_constants
use MeshReader
implicit none

integer cell, dir, k, i, j, gCell,gNeighor
type(defect), pointer :: defectCurrent

do cell=1,numCells
	do dir=1,6
		do k=1,myMesh(cell)%numNeighbors(dir)
			!neighbor cell is in different proc, not free surface (taskid=-1)
			if(myMesh(cell)%neighborProcs(k,dir) /= myProc%taskid .AND. &
					myMesh(cell)%neighborProcs(k,dir) == myProc%procNeighbor(dir) .AND. &
					myMesh(cell)%neighborProcs(k,dir) /= -1) then

				!initialize this boundary element:
				!1) Create defectList (allocate the pointer)
				!2) Set first defect to all species 0, num=0, cellNumber=correct cell number in neighboring proc

				allocate(myBoundary(myMesh(cell)%neighbors(k,dir),dir)%defectList)
				defectCurrent=>myBoundary(myMesh(cell)%neighbors(k,dir),dir)%defectList
				allocate(defectCurrent%defectType(numSpecies))
				do i=1, numSpecies
					defectCurrent%defectType(i)=0
				end do
				defectCurrent%num=0
				defectCurrent%cellNumber=myMesh(cell)%neighbors(k,dir)
				nullify(myBoundary(myMesh(cell)%neighbors(k,dir),dir)%defectList%next)

				!Cu_1
				if(CuContent > 0d0) then
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					do i=1,numSpecies
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%defectType(1)=1
					defectCurrent%num=numCuCell
					defectCurrent%cellNumber=myMesh(cell)%neighbors(k,dir)
					nullify(defectCurrent%next)
				end if

				!V_1
				if(initialTotalV > 0) then
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					do i=1, numSpecies
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%defectType(2)=1
					defectCurrent%num=0
					defectCurrent%cellNumber=myMesh(cell)%neighbors(k,dir)

					do i=1, initialTotalV
						gCell=myMesh(cell)%globalCell
						gNeighor=findgNeighborPeriodic(gCell, dir)
						if(VgCellList(i)==gNeighor) then
							defectCurrent%num=defectCurrent%num +1
						end if
					end do

					if(defectCurrent%num > 0) then
						nullify(defectCurrent%next)
					else
						deallocate(defectCurrent%defectType)
						deallocate(defectCurrent)
						nullify(defectCurrent)
					end if
				end if

				!SIA_1
				if(initialTotalI > 0) then
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					do i=1, numSpecies
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%defectType(3)=1
					defectCurrent%num=0
					defectCurrent%cellNumber=myMesh(cell)%neighbors(k,dir)

					do i=1, initialTotalI
						gCell=myMesh(cell)%globalCell
						gNeighor=findgNeighborPeriodic(gCell, dir)
						if(IgCellList(i)==gNeighor) then
							defectCurrent%num=defectCurrent%num +1
						end if
					end do

					if(defectCurrent%num > 0) then
						nullify(defectCurrent%next)
					else
						deallocate(defectCurrent%defectType)
						deallocate(defectCurrent)
						nullify(defectCurrent)
					end if
				end if

			end if  !myMesh(cell)%neighborProcs(dir,k) /= myProc%taskid
		end do  !dir
	end do  !cell
end do

end subroutine

!***************************************************************************************************
!>Subroutine: Initialize Fine Mesh - initializes defect and reaction lists in a newly created fine mesh
!
!Allocates the size of the fine mesh (fine defect list and reaction list included as part of the cascade derived type),
!populates with defects from the coarse mesh and removes those defects from the coarse mesh.
!
!This subroutine includes cascade defect interaction with pre-existing defects.
!
!Inputs: CascadeCurrent (cascade type derived type variable)
!Output: CascadeCurrent (initialized with defects from coarse mesh)
!***************************************************************************************************

subroutine initializeFineMesh(CascadeCurrent)
use DerivedType
use mod_constants
use randdp
use ReactionRates
implicit none

type(cascade), pointer :: CascadeCurrent
type(defect), pointer :: defectCurrentCoarse, defectCurrentFine
type(defect), pointer :: defectPrevCoarse, defectPrevFine

integer i, j, n, k, num, cell, binomial, factorial, count
double precision volumeRatio, r1, r2, lambda, rstore
integer products(numSpecies), matNum

integer reac

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material
end if

!Allocate cascade defect list and reaction list

allocate(CascadeCurrent%localDefects(numCellsCascade))
allocate(CascadeCurrent%reactionList(numCellsCascade))

!For each cell, initialize the reaction list (no reaction) and defect list

do cell=1,numCellsCascade
	!initialize the defect list for each fine mesh
	allocate(CascadeCurrent%localDefects(cell)%defectType(numSpecies))
	nullify(CascadeCurrent%localDefects(cell)%next)
	do j=1,numSpecies
		CascadeCurrent%localDefects(cell)%defectType(j)=0
	end do
	
	CascadeCurrent%localDefects(cell)%num=0
	CascadeCurrent%localDefects(cell)%cellNumber=cell   !fineMeshID

    !initialize the reaction list for each fine mesh
	CascadeCurrent%reactionList(cell)%numReactants=0
	CascadeCurrent%reactionList(cell)%numProducts=0
	CascadeCurrent%reactionList(cell)%reactionRate=0d0
	nullify(CascadeCurrent%reactionList(cell)%next)

end do

!defectCurrentCoarse and defectPrevCoarse will be used to peruse defects in coarse mesh and place in fine mesh
defectCurrentCoarse=>defectList(CascadeCurrent%cellNumber)
nullify(defectPrevCoarse)

!Use random numbers to assign the defects in the coarse element to the fine mesh

!The volume ratio is the ratio of the fine mesh volume to the coarse mesh volume, used to define 
!the probability that a coarse mesh defect is placed in the fine mesh.
volumeRatio=CascadeElementVol*dble(numCellsCascade)/(myMesh(CascadeCurrent%cellNumber)%volume)

!This loop goes through defects in the coarse mesh and decides whether they should be inserted into
!the fine mesh. If yes, the defects are removed from the coarse mesh. UpdateDefect needs to be kept current
!in this loop.

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
					!write(*,*) 'poisson distribution'
					r2=r2+lambda**(dble(k))*dexp(-lambda)/factorial(k)
					!write(*,*) 'r2', r2, lambda**(dble(k))*dexp(-lambda)/factorial(k)
					!write(*,*) 'r1', r1, 'r2', r2, 'k', k, 'n', n
					if(r2 > r1) then
						exit
					end if
				else    !Gaussian distrubition
					!write(*,*) 'gaussian distribution'
					r2=r2+dexp(-(k-lambda)**(2d0)/(2d0*lambda))/(dsqrt(2*pi*lambda))
					!write(*,*) 'r1', r1, 'r2', r2, 'k', k, 'n', n
					if(rstore-r2 == 0d0) then
						!write(*,*) 'gaussian exit'
						exit
					else if(r2 > r1) then
						!write(*,*) 'gaussian exit'
						exit
					end if
					rstore=r2
				end if
			end do
			!write(*,*) 'k', k
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
					endif
				end do
			
				!***************************************************************
				!Deposit the defect into the fine mesh
				!***************************************************************
				do j=1,numSpecies
					products(j)=defectCurrentCoarse%defectType(j)
				end do
				!write(*,*) 'inserting into fine mesh', (products(j), j=1,numSpecies), 'k', k
				nullify(defectPrevFine)
				defectCurrentFine=>CascadeCurrent%localDefects(cell)
				!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
				! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
			
				call findDefectInList(defectCurrentFine, defectPrevFine, products)
			
				if(associated(defectCurrentFine)) then !if we aren't at the end of the list
				
					count=0
					!Check to see if this defect already exists in the fine mesh list
					do j=1,numSpecies
						if(defectCurrentFine%defectType(j)==products(j)) then
							count=count+1
						endif
					end do
				
					if(count==numSpecies) then
				
						defectCurrentFine%num=defectCurrentFine%num+1
					
					else		!if the defect is to be inserted in the list
					
						if(.NOT. associated(defectPrevFine)) then
							write(*,*) 'error defectPrevFine not associated'
						endif
					
						nullify(defectPrevFine%next)
						allocate(defectPrevFine%next)
						nullify(defectPrevFine%next%next)
						defectPrevFine=>defectPrevFine%next
						allocate(defectPrevFine%defectType(numSpecies))
						defectPrevFine%cellNumber=cell
						defectPrevFine%num=1
						do j=1,numSpecies
							defectPrevFine%defectType(j)=products(j)
						end do
						!if inserted defect is in the middle of the list, point it to the next item in the list
						defectPrevFine%next=>defectCurrentFine
					
					endif
				else 			!add a defect to the end of the list
					nullify(defectPrevFine%next)
					allocate(defectPrevFine%next)
					nullify(defectPrevFine%next%next)
					defectPrevFine=>defectPrevFine%next
					allocate(defectPrevFine%defectType(numSpecies))
					defectPrevFine%cellNumber=cell
					defectPrevFine%num=1
					do j=1,numSpecies
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

				!if there is only one element in the list and there is one of them
!				else if(defectCurrentCoarse%num==1) then 	!removing only defect from cell i (single helium) - this is redundant but will keep for now
!				defectCurrentCoarse%num=0

				endif

			end if
		end if
	end if
	
	defectPrevCoarse=>defectCurrentCoarse
	defectCurrentCoarse=>defectCurrentCoarse%next
	
end do

end subroutine
!***************************************************************************************************
!
!> Subroutine initializeMesh() - begins the mesh initialization process at the beginning of the program
!!
!! This subroutine reads the name of the mesh file from the central input file (parameters.txt) and
!! sends it to readMeshUniform() or readMeshNonUniform().
!!
!! NOTE: although the mesh file and creation of mesh files and connectivity differ for uniform/nonuniform
!! meshes, the format of the final global variable created (class mesh, myMesh, in mod_constants)
!! is the same for both readMeshUniform and readMeshNonUniform. Thus the rest of the program can use
!! it either way.
!!
!! Debug tool: if debugRestart=='yes', then we populate the mesh with the defects in the debug file.
!! This allows us to start the simulation further along than the beginning of the simulation.
!
!***************************************************************************************************

subroutine initializeMesh()
use MeshReader
use mod_constants
use DerivedType
implicit none

character*20 char, meshType
character*50 filename, filename2, filename3
logical flag

flag=.FALSE.

!read in filename of mesh file
do  while (flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='meshFile') then
		read(81,*) filename
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!read in whether mesh is uniform or non-uniform
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='meshType') then
		read(81,*) meshType
		flag=.TRUE.
	end if
end do
flag=.FALSE.

!read in whether we have a strain field or not
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='strainField') then
		read(81,*) strainField
		flag=.TRUE.
	end if
end do
flag=.FALSE.

if(strainField=='yes') then
	do while(flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='strainFile') then
			read(81,*) strainFileName
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	do while(flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='dipoleFile') then
			read(81,*) dipoleFileName
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
end if

!these subroutines (located in MeshReader.f90) initialize the mesh and connectivity.
if(meshType=='uniform') then
!	call readMeshUniform(filename)
	call initialMeshUniform(filename)
else if(meshType=='nonUniform') then
	call readMeshNonUniform(filename)
else
	write(*,*) 'error mesh type unknown'
end if

!Check to see if we are starting from a saved point (from a reset file)
do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char == 'debugRestart') then
		read(81,*) debugToggle
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char == 'debugRestartFile') then
		read(81,*) restartFileName
		flag=.TRUE.
	end if
end do
flag=.FALSE.

end subroutine

!***********************************************************************
!> Subroutine annealInitialization() - changes implantation rates and temperature for annealing
!!
!! This subroutine carries out the following tasks in order to switch from
!! damage introduction to annealing:
!!
!! Temperature changed to annealTemp;
!! Damage rate changed to 0;
!! He implant rate changed to 0
!!
!! Inputs: none (uses global variable annealTemp)
!! Outputs: none
!! Actions: prepares system for annealing
!***********************************************************************

subroutine annealInitialization()
use DerivedType
use mod_constants
implicit none

integer cell
type(reaction), pointer :: reactionCurrent
type(cascade), pointer :: CascadeCurrent

!Change temperature to annealTemp

temperature=annealTemp

!Next, reset implant rates. This step will depend on the number of implant
!types originally in the simulation

!Coarse mesh
do cell=1,numCells
	reactionCurrent=>reactionList(cell)
	
	!remove this reaction rate from totalRate
	totalRate=totalRate-reactionCurrent%reactionRate
	totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
	
	!Change reaction rate to 0 for cascades and He implantation
	reactionCurrent%reactionRate=0d0

end do

end subroutine
