! **************************************************************************************
!>Subroutine initialize vacancy or SIA defect.
!***************************************************************************************
subroutine initializeVIdefect()
use DerivedType
use mod_srscd_constants
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
	if(FormSingle(1,i)%defectType(1)==0 .AND. FormSingle(1,i)%defectType(2)==1 .AND. &
			FormSingle(1,i)%defectType(3)==0 .AND. FormSingle(1,i)%defectType(4)==0) then	!V
		initialCeqv = dexp(-FormSingle(1,i)%Ef / (kboltzmann*temperature))
	else if(FormSingle(1,i)%defectType(1)==0 .AND. FormSingle(1,i)%defectType(2)==0 .AND. &
			FormSingle(1,i)%defectType(3)==1 .AND. FormSingle(1,i)%defectType(4)==0) then	!SIA
		initialCeqi = dexp(-FormSingle(1,i)%Ef / (kboltzmann*temperature))
	end if
end do

initialTotalV = nint(initialCeqv * totalAtoms)
initialTotalI = nint(initialCeqi * totalAtoms)

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

call MPI_BCAST(VgCellList, initialTotalV, MPI_INTEGER, MASTER, MPI_COMM_WORLD,ierr)

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

call MPI_BCAST(IgCellList, initialTotalI, MPI_INTEGER, MASTER, MPI_COMM_WORLD,ierr)


end subroutine

!*****************************************************************************************
!>Subroutine initialize defect list - initializes defect list for each volume element
!
!Begins with a defect with type 0 0 0 0 and num 0. Note that numCells is needed in this subroutine
!*****************************************************************************************

subroutine initializeDefectList()
use DerivedType
use mod_srscd_constants
implicit none

integer cell, i
type(defect), pointer :: defectCurrent

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
	allocate(defectCurrent%next)
	defectCurrent=>defectCurrent%next
	allocate(defectCurrent%defectType(numSpecies))
	do i=1, numSpecies
		defectCurrent%defectType(i) = 0
	end do
	defectCurrent%defectType(1) = 1
	defectCurrent%num=numCuCell
	defectCurrent%cellNumber=cell


	!V_1
	allocate(defectCurrent%next)
	defectCurrent=>defectCurrent%next
	allocate(defectCurrent%defectType(numSpecies))
	do i=1, numSpecies
		defectCurrent%defectType(i) = 0
	end do
	defectCurrent%defectType(2) = 1
	defectCurrent%num=0
	defectCurrent%cellNumber=cell

	if(initialTotalV > 0) then
		do i=1, initialTotalV
			if(VgCellList(i)==myMesh(cell)%globalCell) then
				defectCurrent%num = defectCurrent%num +1
			end if
		end do
	end if

	!SIA_1
	allocate(defectCurrent%next)
	defectCurrent=>defectCurrent%next
	allocate(defectCurrent%defectType(numSpecies))
	do i=1, numSpecies
		defectCurrent%defectType(i) = 0
	end do
	defectCurrent%defectType(3) = 1
	defectCurrent%num=0
	defectCurrent%cellNumber=cell

	if(initialTotalI > 0) then
		do i=1, initialTotalI
			if(IgCellList(i)==myMesh(cell)%globalCell) then
				defectCurrent%num = defectCurrent%num +1
			end if
		end do
	end if

	nullify(defectCurrent%next)

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
use mod_srscd_constants
use DerivedType
use randdp
implicit none
include 'mpif.h'

integer randseed, i, irand
integer status(MPI_STATUS_SIZE)

if(myProc%taskid == MASTER) then
	call system_clock(Count=randseed)	!return randseed (integer, unit:ms)
	call sdprnd(randseed)
	write(*,*)
	write(*,*) 'random number seed', randseed, 'processor', myProc%taskid
	do i=1,myProc%numtasks-1
		randseed=irand(randseed)
		call mpi_send(randseed, 1, MPI_INTEGER, i, 1000,MPI_COMM_WORLD, ierr)
	end do
else
	call mpi_recv(randseed, 1, MPI_INTEGER, MASTER, 1000, MPI_COMM_WORLD, status, ierr)
	call sdprnd(randseed)
	!write(*,*)
	write(*,*) 'random number seed', randseed, 'processor', myProc%taskid
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
use mod_srscd_constants
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
use mod_srscd_constants
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
			allocate(reactionList(cell)%products(reactionList(cell)%numProducts,numSpecies))
			allocate(reactionList(cell)%cellNumber(reactionList(cell)%numProducts))
			allocate(reactionList(cell)%taskid(reactionList(cell)%numProducts))
		
			!search ImplantList for Frenkel Pair reactions
			do reac=1,numImplantReac(matNum)
				if(ImplantReactions(matNum,reac)%numReactants==0 .AND. &
						ImplantReactions(matNum,reac)%numProducts==2) then	!we have found Frenkel pair implantation
					exit	!leave ImplantReactions(reac) pointed at this element of ImplantReactions
				end if
			end do
		
			do i=1,ImplantReactions(matNum,reac)%numProducts
				do j=1,numSpecies
					reactionList(cell)%products(i,j)=ImplantReactions(matNum,reac)%products(i,j)
				end do
				reactionList(cell)%cellNumber(i)=cell
				reactionList(cell)%taskid(i)=myMesh(cell)%proc
			end do
		
			!Find reaction rate for Frenkel pair implantation using ImplantReactions(reac), which is input info from file.
			reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
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
				if(ImplantReactions(matNum,reac)%numReactants==-10 .AND. &
						ImplantReactions(matNum,reac)%numProducts==0) then	!we have found cascade implantation
					exit
				end if
			end do

			!If implant scheme is Monte Carlo, then this reaction is chosen as part of the regular
			!Monte Carlo algorithm. Otherwise, cascades are added explicitly and this rate is set to 0
			!so that cascade implantation is no longer chosen as part of the MC algorithm.

			if(implantScheme=='MonteCarlo') then
				reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
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

	if(CuContent > 0d0) then
	!2019.05.04 Add: Initialize possible reactions of free Cu
	!*******************************************************
	!clustering: Cu+Cu->2Cu
	!*******************************************************

		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=2
		reactionCurrent%numProducts=1
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

		!search ClusterList for Cu+Cu->2Cu
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%reactants(1,1)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,1)==1 &
					.AND. ClusterReactions(matNum,reac)%reactants(1,2)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,2)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,3)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,3)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,4)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,4)==0) then
				exit
			end if
		end do

		do i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			do j=1,numSpecies
				reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
				reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
				reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
			end do
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))
		nullify(reactionCurrent%next)

	end if

	!*******************************************************************
	!Diffusion: Cu->Cu
	!*******************************************************************
!	do dir=1, 6
!		allocate(reactionCurrent%next)
!		reactionCurrent=>reactionCurrent%next

!		reactionCurrent%numReactants=1
!		reactionCurrent%numProducts=1
!		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
!		allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
!		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
!		allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

	!search DiffList for Cu->Cu
!		do reac=1,numDiffReac(matNum)
!			if(DiffReactions(matNum,reac)%reactants(1,1)==1 .AND. DiffReactions(matNum,reac)%products(1,1)==1 &
!				.AND. DiffReactions(matNum,reac)%reactants(1,2)==0 .AND. &
!				DiffReactions(matNum,reac)%products(1,2)==0 .AND. &
!				DiffReactions(matNum,reac)%reactants(1,3)==0 .AND. &
!				DiffReactions(matNum,reac)%products(1,3)==0 .AND. &
!				DiffReactions(matNum,reac)%reactants(1,4)==0 .AND. &
!				DiffReactions(matNum,reac)%products(1,4)==0) then
!				exit
!			end if
!		end do

	!do 32 i=1,DiffReactions(matNum,reac)%numReactants+DiffReactions(matNum,reac)%numProducts
!		reactionCurrent%reactants(1,1)=1
!		reactionCurrent%products(1,1)=1
!		do j=2,numSpecies
!			reactionCurrent%reactants(1,j)=DiffReactions(matNum,reac)%reactants(1,j)
!			reactionCurrent%products(1,j)=DiffReactions(matNum,reac)%products(1,j)
!		end do
!		reactionCurrent%cellNumber(1)=cell
!		reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(dir,1)
!		reactionCurrent%taskid(1)=myMesh(cell)%proc
!		reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(dir,1)
	!32 continue
!		reactionCurrent%reactionRate = 0d0
	!reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(1,:), cell, &
	!		myProc%taskid, myMesh(cell)%neighbors(dir,1), myMesh(cell)%neighborProcs(dir,1), dir, &
	!		DiffReactions(matNum,reac))
!		nullify(reactionCurrent%next)

!	end do
	!*******************************************************************
	!Sink trapping: Cu->0
	!*******************************************************************

!	allocate(reactionCurrent%next)
!	reactionCurrent=>reactionCurrent%next

!	reactionCurrent%numReactants=1
!	reactionCurrent%numProducts=0
!	allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
!	allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
!	allocatE(reactionCurrent%taskid(reactionCurrent%numReactants))

!	do reac=1,numSinkReac(matNum)
!		if(SinkReactions(matNum,reac)%reactants(1,1)==1 .AND. SinkReactions(matNum,reac)%reactants(1,2)==0 .AND. &
!			SinkReactions(matNum,reac)%reactants(1,3)==0 .AND. SinkReactions(matNum,reac)%reactants(1,4)==0) then
!			exit
!		end if
!	end do

!	do i=1,SinkReactions(matNum,reac)%numReactants
!		reactionCurrent%reactants(1,1)=1
!		do j=2,numSpecies
!			reactionCurrent%reactants(1,j)=SinkReactions(matNum,reac)%reactants(1,j)
!		end do
!		reactionCurrent%cellNumber(i)=cell
!		reactionCurrent%taskid(i)=myMesh(cell)%proc
!	end do
!	reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(1,:), cell, &
!			SinkReactions(matNum,reac))
!	nullify(reactionCurrent%next)

	!**********************************************
	!initialize reactions about V
	defectCurrent=>defectList(cell)
	do while(associated(defectCurrent))
		if(defectCurrent%defectType(1)==0 &
				.AND. defectCurrent%defectType(2)==1 &
				.AND. defectCurrent%defectType(3)==0 &
				.AND. defectCurrent%defectType(4)==0) then
			exit
		else
			defectCurrent=>defectCurrent%next
		end if
	end do
	defectCurrentTempV=>defectCurrent

	if(associated(defectCurrent)) then
	if(defectCurrent%num > 0) then

		if(CuContent > 0d0) then
			!*******************************************************
			!clustering: Cu+V->CuV
			!*******************************************************

			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

			!search ClusterList for Cu+V->CuV
			do reac=1,numClusterReac(matNum)
				if(ClusterReactions(matNum,reac)%reactants(1,1)==1 .AND. &
						ClusterReactions(matNum,reac)%reactants(2,1)==0 &
						.AND. ClusterReactions(matNum,reac)%reactants(1,2)==0 .AND. &
						ClusterReactions(matNum,reac)%reactants(2,2)==1 .AND. &
						ClusterReactions(matNum,reac)%reactants(1,3)==0 .AND. &
						ClusterReactions(matNum,reac)%reactants(2,3)==0 .AND. &
						ClusterReactions(matNum,reac)%reactants(1,4)==0 .AND. &
						ClusterReactions(matNum,reac)%reactants(2,4)==0) then
					exit
				end if
			end do

			do i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
				do j=1,numSpecies
					reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
					reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
					reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
				end do

				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			end do
			!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
			reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
					reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))
			nullify(reactionCurrent%next)

		end if
		!*******************************************************************
		!Diffusion: V->V
		!*******************************************************************
		do dir=1, 6
			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

		!search DiffList for V->V
			do reac=1,numDiffReac(matNum)
				if(DiffReactions(matNum,reac)%reactants(1,1)==0 .AND. DiffReactions(matNum,reac)%products(1,1)==0 &
						.AND. DiffReactions(matNum,reac)%reactants(1,2)==1 .AND. &
						DiffReactions(matNum,reac)%products(1,2)==1 .AND. &
						DiffReactions(matNum,reac)%reactants(1,3)==0 .AND. &
						DiffReactions(matNum,reac)%products(1,3)==0 .AND. &
						DiffReactions(matNum,reac)%reactants(1,4)==0 .AND. &
						DiffReactions(matNum,reac)%products(1,4)==0) then
					exit
				end if
			end do

			do i=1,DiffReactions(matNum,reac)%numReactants+DiffReactions(matNum,reac)%numProducts
				do j=1,numSpecies
					reactionCurrent%reactants(1,j)=DiffReactions(matNum,reac)%reactants(1,j)
					reactionCurrent%products(1,j)=DiffReactions(matNum,reac)%products(1,j)
				end do
				reactionCurrent%cellNumber(1)=cell
				reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(dir,1)
				reactionCurrent%taskid(1)=myMesh(cell)%proc
				reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(dir,1)
			end do
			reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(1,:), cell, &
				myProc%taskid, myMesh(cell)%neighbors(dir,1), myMesh(cell)%neighborProcs(dir,1), dir, &
				DiffReactions(matNum,reac))
			nullify(reactionCurrent%next)

		end do

		!*******************************************************************
		!sinkRemoval: V->0
		!*******************************************************************
		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=1
		reactionCurrent%numProducts=0
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

		do reac=1,numSinkReac(matNum)
			if(SinkReactions(matNum,reac)%reactants(1,1)==0 .AND. SinkReactions(matNum,reac)%reactants(1,2)==1 .AND. &
					SinkReactions(matNum,reac)%reactants(1,3)==0 .AND. SinkReactions(matNum,reac)%reactants(1,4)==0) then
				exit
			end if
		end do

		do i=1,SinkReactions(matNum,reac)%numReactants
			do j=1,numSpecies
				reactionCurrent%reactants(i,j)=SinkReactions(matNum,reac)%reactants(i,j)
			end do
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(1,:), cell, &
					SinkReactions(matNum,reac))
		nullify(reactionCurrent%next)
	end if

	if(defectCurrent%num > 1) then
		!*******************************************************
		!clustering: V+V->2V
		!*******************************************************

		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=2
		reactionCurrent%numProducts=1
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

		!search ClusterList for V+V->2V
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%reactants(1,1)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,1)==0 &
					.AND. ClusterReactions(matNum,reac)%reactants(1,2)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,2)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,3)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,3)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,4)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,4)==0) then
				exit
			end if
		end do

		do i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			do j=1,numSpecies
				reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
				reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
				reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
			end do
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))
		nullify(reactionCurrent%next)

	end if

	end if

	!*******************************************************
	!initialize reactions about SIA
	defectCurrent=>defectList(cell)
	do while(associated(defectCurrent))
		if(defectCurrent%defectType(1)==0 &
				.AND. defectCurrent%defectType(2)==0 &
				.AND. defectCurrent%defectType(3)==1 &
				.AND. defectCurrent%defectType(4)==0) then
			exit
		else
			defectCurrent=>defectCurrent%next
		end if
	end do
	defectCurrentTempI=>defectCurrent

	if(associated(defectCurrent)) then
	if(defectCurrent%num > 0) then

		!*******************************************************************
		!Diffusion: SIA->SIA
		!*******************************************************************
		do dir=1, 6
			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

			!search DiffList for V->V
			do reac=1,numDiffReac(matNum)
				if(DiffReactions(matNum,reac)%reactants(1,1)==0 .AND. DiffReactions(matNum,reac)%products(1,1)==0 &
						.AND. DiffReactions(matNum,reac)%reactants(1,2)==0 .AND. &
						DiffReactions(matNum,reac)%products(1,2)==0 .AND. &
						DiffReactions(matNum,reac)%reactants(1,3)==1 .AND. &
						DiffReactions(matNum,reac)%products(1,3)==1 .AND. &
						DiffReactions(matNum,reac)%reactants(1,4)==0 .AND. &
						DiffReactions(matNum,reac)%products(1,4)==0) then
					exit
				end if
			end do

			do i=1,DiffReactions(matNum,reac)%numReactants+DiffReactions(matNum,reac)%numProducts
				do j=1,numSpecies
					reactionCurrent%reactants(1,j)=DiffReactions(matNum,reac)%reactants(1,j)
					reactionCurrent%products(1,j)=DiffReactions(matNum,reac)%products(1,j)
				end do
				reactionCurrent%cellNumber(1)=cell
				reactionCurrent%cellNumber(2)=myMesh(cell)%neighbors(dir,1)
				reactionCurrent%taskid(1)=myMesh(cell)%proc
				reactionCurrent%taskid(2)=myMesh(cell)%neighborProcs(dir,1)
			end do
			reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(1,:), cell, &
					myProc%taskid, myMesh(cell)%neighbors(dir,1), myMesh(cell)%neighborProcs(dir,1), dir, &
					DiffReactions(matNum,reac))
			nullify(reactionCurrent%next)

		end do

		!*******************************************************************
		!sinkRemoval: SIA->0
		!*******************************************************************
		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=1
		reactionCurrent%numProducts=0
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants))

		do reac=1,numSinkReac(matNum)
			if(SinkReactions(matNum,reac)%reactants(1,1)==0 .AND. SinkReactions(matNum,reac)%reactants(1,2)==0 .AND. &
					SinkReactions(matNum,reac)%reactants(1,3)==1 .AND. SinkReactions(matNum,reac)%reactants(1,4)==0) then
				exit
			end if
		end do

		do i=1,SinkReactions(matNum,reac)%numReactants
			do j=1,numSpecies
				reactionCurrent%reactants(i,j)=SinkReactions(matNum,reac)%reactants(i,j)
			end do
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(1,:), cell, &
				SinkReactions(matNum,reac))
		nullify(reactionCurrent%next)
	end if

	if(defectCurrent%num > 1) then
		!*******************************************************
		!clustering: SIA+SIA->2SIA
		!*******************************************************

		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=2
		reactionCurrent%numProducts=1
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

		!search ClusterList for V+V->2V
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%reactants(1,1)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,1)==0 &
					.AND. ClusterReactions(matNum,reac)%reactants(1,2)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,2)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,3)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,3)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,4)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,4)==0) then
				exit
			end if
		end do

		do j=1,numSpecies
			reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
			reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
			reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
		end do

		if(pointDefectToggle=='yes') then
			if(reactionCurrent%products(1,3) > max3DInt) then
				reactionCurrent%products(1,4) = reactionCurrent%products(1,3)
				reactionCurrent%products(1,3) = 0
			end if
		end if

		do i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))
		nullify(reactionCurrent%next)

	end if

	end if

	!***********************************
	!V+SIA->0
	if(associated(defectCurrentTempV) .AND. associated(defectCurrentTempI)) then
	if(defectCurrentTempV%num > 0 .AND. defectCurrentTempI%num > 0) then
		!*******************************************************
		!clustering: V+SIA->0
		!*******************************************************
		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=2
		reactionCurrent%numProducts=0
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
		allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

		!search ClusterList for V+V->2V
		do reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%reactants(1,1)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,1)==0 &
					.AND. ClusterReactions(matNum,reac)%reactants(1,2)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,2)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,3)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,3)==1 .AND. &
					ClusterReactions(matNum,reac)%reactants(1,4)==0 .AND. &
					ClusterReactions(matNum,reac)%reactants(2,4)==0) then
				exit
			end if
		end do

		do i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			do j=1,numSpecies
				reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
				reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
			end do
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		end do
		!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))
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
use mod_srscd_constants
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
	use mod_srscd_constants
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
		if(char=='numHeImplantEvents') then
			read(87,*) numHeImplantEventsReset
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
use mod_srscd_constants
use MeshReader
implicit none

integer cell, dir, k, i, j, gCell,gNeighor
type(defect), pointer :: defectCurrent

do cell=1,numCells
	do dir=1,6
		do k=1,myMesh(cell)%numNeighbors(dir)
			if(myMesh(cell)%neighborProcs(dir,k) /= myProc%taskid .AND. &
				myMesh(cell)%neighborProcs(dir,k) /= -1) then !neighbor cell is in different proc, not free surface (taskid=-1)
				
				if(myMesh(cell)%neighborProcs(dir,k) /= myProc%procNeighbor(dir)) then
					write(*,*) 'error neighbor not correct during boundary mesh initialization'
				else
					!initialize this boundary element:
					!1) Create defectList (allocate the pointer)
					!2) Set first defect to all species 0, num=0, cellNumber=correct cell number in neighboring proc
					
!					if(myProc%taskid==MASTER) then
!						write(*,*) 'proc', myProc%taskid, 'cell', cell, 'dir', dir,'neighbor cell', myMesh(cell)%neighbors(dir,k)
!						write(*,*) 'neighbor proc', myBoundary(dir,myMesh(cell)%neighbors(dir,k))%proc
!						!write(*,*) 'length', myBoundary(dir,myMesh(cell)%neighbors(dir,k))%length
!						!write(*,*) 'material', myBoundary(dir,myMesh(cell)%neighbors(dir,k))%material
!					endif

					allocate(myBoundary(dir,myMesh(cell)%neighbors(dir,k))%defectList)
					defectCurrent=>myBoundary(dir,myMesh(cell)%neighbors(dir,k))%defectList
					allocate(defectCurrent%defectType(numSpecies))
					do i=1, numSpecies
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%num=0
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir,k)


					!Cu_1
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					defectCurrent%defectType(1)=1
					do i=2,numSpecies
						defectCurrent%defectType(i)=0
					end do
					defectCurrent%num=numCuCell
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir,k)


					!V_1
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					do i=1, numSpecies
						defectCurrent%defectType(i) = 0
					end do
					defectCurrent%defectType(2) = 1
					defectCurrent%num=0
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir,k)

					if(initialTotalV > 0) then
						do i=1, initialTotalV
							gCell=myMesh(cell)%globalCell
							gNeighor=findgNeighborPeriodic(gCell, dir)
							if(VgCellList(i)==gNeighor) then
								defectCurrent%num = defectCurrent%num +1
							end if
						end do
					end if


					!SIA_1
					allocate(defectCurrent%next)
					defectCurrent=>defectCurrent%next
					allocate(defectCurrent%defectType(numSpecies))
					do i=1, numSpecies
						defectCurrent%defectType(i) = 0
					end do
					defectCurrent%defectType(3) = 1
					defectCurrent%num=0
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir,k)
					if(initialTotalI > 0) then
						do i=1, initialTotalI
							gCell=myMesh(cell)%globalCell
							gNeighor=findgNeighborPeriodic(gCell, dir)
							if(IgCellList(i)==gNeighor) then
								defectCurrent%num = defectCurrent%num +1
							end if
						end do
					end if

					nullify(defectCurrent%next)

				end if  !myMesh(cell)%neighborProcs(dir,k) == myProc%procNeighbor(dir)

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
use mod_srscd_constants
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
	use mod_srscd_constants
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
			
				!***************************************************************
				!Remove the defect from the coarse mesh
				!***************************************************************
!				if(.NOT. associated(defectCurrentCoarse)) then
!					write(*,*) 'Tried to delete defect that wasnt there fine mesh initialization'
			
				!if defectCurrentCoarse is in the middle of the list and there is 1 of them
!				else if(defectCurrentCoarse%num==1 .AND. associated(defectCurrentCoarse%next) .AND. associated(defectPrevCoarse)) then
				
!					defectPrevCoarse%next=>defectCurrentCoarse%next !remove that defect type from the system
!					deallocate(defectCurrentCoarse%defectType)
!					deallocate(defectCurrentCoarse)
!					defectCurrentCoarse=>defectPrevCoarse
			
				!if defectCurrentCoarse is at the end of the list and there is one of them
!				else if(defectCurrentCoarse%num==1 .AND. associated(defectPrevCoarse)) then
!					deallocate(defectCurrentCoarse%defectType)
!					deallocate(defectCurrentCoarse)
!					defectCurrentCoarse=>defectPrevCoarse	!remove the last defect from the system
!					nullify(defectPrevCoarse%next)
			
				!if defectCurrentCoarse is at the beginning of the list and there is one of them
!				else if(defectCurrentCoarse%num==1 .AND. associated(defectCurrentCoarse%next)) then !removing first defect from cell i
!					defectCurrentCoarse%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
			
				!if there is only one element in the list and there is one of them
!				else if(defectCurrentCoarse%num==1) then 	!removing only defect from cell i (single helium) - this is redundant but will keep for now
!					defectCurrentCoarse%num=0
			
				!if we are trying to remove more defects than are actually present (the calculator in prev. step chose k incorrectly)
!				else if(defectCurrentCoarse%num==0) then
!					write(*,*) 'trying to remove a defect that isnt there (fine mesh)'
!					write(*,*) 'k', k, 'num', num, 'num defects coarse', n
!				else
			
				!If there is more than 1 defect of this type, we don't have to to do anything with the pointers
				!and just subtract 1 from the number of defects in the system.
!					defectCurrentCoarse%num=defectCurrentCoarse%num-1 !remove the defect from the system instead of the entire entry in the list
!				endif
			
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
!! meshes, the format of the final global variable created (class mesh, myMesh, in mod_srscd_constants)
!! is the same for both readMeshUniform and readMeshNonUniform. Thus the rest of the program can use
!! it either way.
!!
!! Debug tool: if debugRestart=='yes', then we populate the mesh with the defects in the debug file.
!! This allows us to start the simulation further along than the beginning of the simulation.
!
!***************************************************************************************************

subroutine initializeMesh()
use MeshReader
use mod_srscd_constants
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
use mod_srscd_constants
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
	
!	if(HeDPARatio .GT. 0d0) then
!		reactionCurrent=>reactionCurrent%next
		
		!remove this reaction rate from totalRate
!		totalRate=totalRate-reactionCurrent%reactionRate
!		totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
		
!		reactionCurrent%reactionRate=0d0
!	endif
	
end do

!Fine mesh

!Here, we only zero the first reaction rate if we have He implantation
!(normally, no defects are implanted into the fine mesh except He)

!if(HeDPARatio .GT. 0d0) then

!	CascadeCurrent=>ActiveCascades
	
!	do 11 while(associated(CascadeCurrent))
	
!		do 12 cell=1,numCellsCascade
			
!			reactionCurrent=>CascadeCurrent%reactionList(cell)
			
			!remove this reaction rate from totalRate
!			totalRate=totalRate-reactionCurrent%reactionRate
!			CascadeCurrent%totalRate=CascadeCurrent%totalRate-reactionCurrent%reactionRate
			
!			reactionCurrent%reactionRate=0d0
		
!		12 continue
		
!		CascadeCurrent=>CascadeCurrent%next
		
!	11 continue
	
!endif

end subroutine
