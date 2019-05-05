! $Header: /home/CVS//srscd/src/Initialization_subroutines.f90,v 1.10 2015/10/19 19:43:22 aydunn Exp $

!*****************************************************************************************
!>Subroutine initialize defect list - initializes defect list for each volume element
!!
!!Begins with a defect with type 0 0 0 0 and num 0. Note that numCells is needed in this 
!!subroutine
!*****************************************************************************************

subroutine initializeDefectList()
use DerivedType
use mod_srscd_constants
implicit none

integer cell, i

do 10 cell=1,numCells
	allocate(DefectList(cell)%defectType(numSpecies))
	DefectList(cell)%defectType(1)=1
	do 11 i=2,numSpecies
		DefectList(cell)%defectType(i)=0
	11 continue
	DefectList(cell)%num=CuAtomsEverMesh
	DefectList(cell)%cellNumber=cell
	nullify(DefectList(cell)%next)
10 continue

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

if(myProc%taskid .eq. MASTER) then
	call system_clock(Count=randseed)
	call sdprnd(randseed)
	write(*,*)
	write(*,*) 'random number seed', randseed, 'processor', myProc%taskid
	do 22 i=1,myProc%numtasks-1
		randseed=irand(randseed)
		call mpi_send(randseed, 1, MPI_INTEGER, i, 1000,MPI_COMM_WORLD, ierr)
	22 continue
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

do 10 cell=1,numCells
	rateCell=0d0
	reactionCurrent=>reactionList(cell)
	
	do 11 while(associated(reactionCurrent))
		rate=rate+reactionCurrent%reactionRate
		rateCell=rateCell+reactionCurrent%reactionRate
		reactionCurrent=>reactionCurrent%next
	11 continue
	
	totalRateVol(cell)=rateCell	!Total reaction rate in this volume element
	
	if(singleElemKMC=='yes') then
		if(rateCell .GT. rateSingle) then
			rateSingle=rateCell
		endif
	endif
	
10 continue

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

integer cell, i, j, reac, matNum, dir
type(reaction), pointer :: reactionCurrent

do 10 cell=1,numCells
	!Initialize reaction rates with Frenkel pair implantation. The rate should be given by a function
	!findReactionRate which has as parameters passed to it the cell that this rate is occuring in and
	!the parameters of the reaction rate as read in from the file (type reactionParameters)
	
	!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
	!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
	if(numMaterials==1) then
		matNum=1
	else
		matNum=myMesh(cell)%material
	endif
	
	if(implantType=='FrenkelPair') then

		!Frenkel pair implantaiton reaction: first in the list
		
		reactionList(cell)%numReactants=0
		reactionList(cell)%numProducts=2
		allocate(reactionList(cell)%products(reactionList(cell)%numProducts,numSpecies))
		allocate(reactionList(cell)%cellNumber(reactionList(cell)%numProducts))
		allocate(reactionList(cell)%taskid(reactionList(cell)%numProducts))
		
		!search ImplantList for Frenkel Pair reactions
		do 12 reac=1,numImplantReac(matNum)
			if(ImplantReactions(matNum,reac)%numReactants==0 .AND. ImplantReactions(matNum,reac)%numProducts==2) then	!we have found Frenkel pair implantation
				exit	!leave ImplantReactions(reac) pointed at this element of ImplantReactions
			endif
		12 continue
		
		do 13 i=1,ImplantReactions(matNum,reac)%numProducts
			do 11 j=1,numSpecies
				reactionList(cell)%products(i,j)=ImplantReactions(matNum,reac)%products(i,j)
			
			11 continue
			reactionList(cell)%cellNumber(i)=cell
			reactionList(cell)%taskid(i)=myMesh(cell)%proc
		13 continue
		
		!Find reaction rate for Frenkel pair implantation using ImplantReactions(reac), which is input info from file.
		reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
		nullify(reactionList(cell)%next)

		!2019.05.04 Add: Initialize possible reactions of free Cu

		reactionCurrent=>reactionList(cell)

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
		do 14  reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%numReactants==2 .AND. ClusterReactions(matNum,reac)%numProducts==1) then
				exit
			endif
		14 continue

		do 15 i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			reactionCurrent%reactants(1,1)=1
			reactionCurrent%reactants(2,1)=1
			reactionCurrent%products(1,1)=1
			do 16 j=2,numSpecies
				reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
				reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
				reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
			16 continue
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		15 continue
			!Find reaction rate for Cu clustering using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))

		!*******************************************************************
		!Diffusion: Cu->Cu
		!*******************************************************************
		do 30 dir=1, 6
			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

			!search DiffList for Cu->Cu
			do 31  reac=1,numDiffReac(matNum)
				if(DiffReactions(matNum,reac)%numReactants==1 .AND. DiffReactions(matNum,reac)%numProducts==1) then
					exit
				endif
			31 continue

			do 32 i=1,DiffReactions(matNum,reac)%numReactants+DiffReactions(matNum,reac)%numProducts
				reactionCurrent%reactants(1,1)=1
				reactionCurrent%products(1,1)=1
				do 33 j=2,numSpecies
					reactionCurrent%reactants(1,j)=DiffReactions(matNum,reac)%reactants(1,j)
					reactionCurrent%products(1,j)=DiffReactions(matNum,reac)%products(1,j)
				33 continue
				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			32 continue
			reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(1,:), cell, &
					myProc%taskid, myMesh(cell)%neighbors(1,dir), myMesh(cell)%neighborProcs(1,dir), dir, &
					DiffReactions(matNum,reac))
		30 continue
		!*******************************************************************
		!Sink trapping: Cu->0
		!*******************************************************************

		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=1
		reactionCurrent%numProducts=0
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
		allocatE(reactionCurrent%taskid(reactionCurrent%numReactants))

		do 34  reac=1,numSinkReac(matNum)
			if(SinkReactions(matNum,reac)%numReactants==1 .AND. SinkReactions(matNum,reac)%numProducts==0) then
				exit
			endif
		34 continue

		do 35 i=1,SinkReactions(matNum,reac)%numReactants
			reactionCurrent%reactants(1,1)=1
			do 36 j=2,numSpecies
				reactionCurrent%reactants(1,j)=SinkReactions(matNum,reac)%reactants(1,j)
			36 continue
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		35 continue
		reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(1,:), cell, &
				SinkReactions(matNum,reac))

		nullify(reactionCurrent%next)

!********************************************************************
!2019.05.04
!		!Helium implantation reaction: second in the list

!		if(HeDPARatio .GT. 0d0) then
		
!			reactionCurrent=>reactionList(cell)
			
!			allocate(reactionCurrent%next)
!			reactionCurrent=>reactionCurrent%next
			
!			reactionCurrent%numReactants=0
!			reactionCurrent%numProducts=1
			
!			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
!			allocate(reactionCurrent%cellNumber(reactionCurrent%numProducts))
!			allocatE(reactionCurrent%taskid(reactionCurrent%numProducts))
			
			!search ImplantList for He implant reactions
!			do 14  reac=1,numImplantReac(matNum)
!				if(ImplantReactions(matNum,reac)%numReactants==0 .AND. ImplantReactions(matNum,reac)%numProducts==1) then
!					exit
!				endif
!			14 continue
			
!			do 15 i=1,ImplantReactions(matNum,reac)%numProducts
!				do 16 j=1,numSpecies
!					reactionCurrent%products(i,j)=ImplantReactions(matNum,reac)%products(i,j)
!				16 continue
!				reactionCurrent%cellNumber(i)=cell
!				reactionCurrent%taskid(i)=myMesh(cell)%proc
!			15 continue
			
			!Find reaction rate for He ion implantation using ImplantReactions(reac), which is input from file.
!			reactionCurrent%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
!			nullify(reactionCurrent%next)
			
!		endif
!***********************************************************************
		
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
		do 22 reac=1,numImplantReac(matNum)
			if(ImplantReactions(matNum,reac)%numReactants==-10 .AND. ImplantReactions(matNum,reac)%numProducts==0) then	!we have found cascade implantation
				exit
			endif
		22 continue
		
		!If implant scheme is Monte Carlo, then this reaction is chosen as part of the regular 
		!Monte Carlo algorithm. Otherwise, cascades are added explicitly and this rate is set to 0
		!so that cascade implantation is no longer chosen as part of the MC algorithm.
		
		if(implantScheme=='MonteCarlo') then
			reactionList(cell)%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
		else if(implantScheme=='explicit') then
			reactionList(cell)%reactionRate=0d0
		endif
		
		nullify(reactionList(cell)%next)

		!2019.05.05 Add: Initialize possible reactions of free Cu

		reactionCurrent=>reactionList(cell)

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
		do 41  reac=1,numClusterReac(matNum)
			if(ClusterReactions(matNum,reac)%numReactants==2 .AND. ClusterReactions(matNum,reac)%numProducts==1) then
				exit
			endif
		41 continue

		do 42 i=1,ClusterReactions(matNum,reac)%numReactants+ClusterReactions(matNum,reac)%numProducts
			reactionCurrent%reactants(1,1)=1
			reactionCurrent%reactants(2,1)=1
			reactionCurrent%products(1,1)=1
			do 43 j=2,numSpecies
				reactionCurrent%reactants(1,j)=ClusterReactions(matNum,reac)%reactants(1,j)
				reactionCurrent%reactants(2,j)=ClusterReactions(matNum,reac)%reactants(2,j)
				reactionCurrent%products(1,j)=ClusterReactions(matNum,reac)%products(1,j)
			43 continue
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		42 continue
		!Find reaction rate for Cu clustering  using ClusterReactions(reac), which is input from file.
		reactionCurrent%reactionRate=findReactionRateMultiple(reactionCurrent%reactants(1,:), &
				reactionCurrent%reactants(2,:), cell, ClusterReactions(matNum,reac))

		!*******************************************************************
		!Diffusion: Cu->Cu
		!*******************************************************************

		do 40 dir=1,6
			allocate(reactionCurrent%next)
			reactionCurrent=>reactionCurrent%next

			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants + reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants + reactionCurrent%numProducts))

			!search DiffList for Cu->Cu
			do 44  reac=1,numDiffReac(matNum)
				if(DiffReactions(matNum,reac)%numReactants==1 .AND. DiffReactions(matNum,reac)%numProducts==1) then
					exit
				endif
			44 continue

			do 45 i=1,DiffReactions(matNum,reac)%numReactants+DiffReactions(matNum,reac)%numProducts
				reactionCurrent%reactants(1,1)=1
				reactionCurrent%products(1,1)=1
				do 46 j=2,numSpecies
					reactionCurrent%reactants(1,j)=DiffReactions(matNum,reac)%reactants(1,j)
					reactionCurrent%products(1,j)=DiffReactions(matNum,reac)%products(1,j)
				46 continue
				reactionCurrent%cellNumber(i)=cell
				reactionCurrent%taskid(i)=myMesh(cell)%proc
			45 continue
			reactionCurrent%reactionRate=findReactionRateDiff(reactionCurrent%reactants(1,:), cell, &
					myProc%taskid, myMesh(cell)%neighbors(1,dir), myMesh(cell)%neighborProcs(1,dir), dir, &
					DiffReactions(matNum,reac))
		40 continue
		!*******************************************************************
		!Sink trapping: Cu->0
		!*******************************************************************

		allocate(reactionCurrent%next)
		reactionCurrent=>reactionCurrent%next

		reactionCurrent%numReactants=1
		reactionCurrent%numProducts=0
		allocate(reactionCurrent%reactants(reactionCurrent%numReactants, numSpecies))
		allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants))
		allocatE(reactionCurrent%taskid(reactionCurrent%numReactants))

		do 47  reac=1,numSinkReac(matNum)
			if(SinkReactions(matNum,reac)%numReactants==1 .AND. SinkReactions(matNum,reac)%numProducts==0) then
				exit
			endif
		47 continue

		do 48 i=1,SinkReactions(matNum,reac)%numReactants
			reactionCurrent%reactants(1,1)=1
			do 49 j=2,numSpecies
				reactionCurrent%reactants(1,j)=SinkReactions(matNum,reac)%reactants(1,j)
			49 continue
			reactionCurrent%cellNumber(i)=cell
			reactionCurrent%taskid(i)=myMesh(cell)%proc
		48 continue
		reactionCurrent%reactionRate=findReactionRateSink(reactionCurrent%reactants(1,:), cell, &
				SinkReactions(matNum,reac))

		nullify(reactionCurrent%next)

!********************************************************************
!2019.05.05
!		!Helium implantation reaction: second in the list
		
!		if(HeDPARatio .GT. 0d0) then
		
!			reactionCurrent=>reactionList(cell)
			
!			allocate(reactionCurrent%next)
!			reactionCurrent=>reactionCurrent%next
			
!			reactionCurrent%numReactants=0
!			reactionCurrent%numProducts=1
			
!			allocate(reactionCurrent%products(reactionCurrent%numProducts, numSpecies))
!			allocate(reactionCurrent%cellNumber(reactionCurrent%numProducts))
!			allocatE(reactionCurrent%taskid(reactionCurrent%numProducts))
			
			!search ImplantList for He implant reactions
!			do 24  reac=1,numImplantReac(matNum)
!				if(ImplantReactions(matNum,reac)%numReactants==0 .AND. ImplantReactions(matNum,reac)%numProducts==1) then
!					exit
!				endif
!			24 continue
			
!			do 25 i=1,ImplantReactions(matNum,reac)%numProducts
!				do 26 j=1,numSpecies
!					reactionCurrent%products(i,j)=ImplantReactions(matNum,reac)%products(i,j)
!				26 continue
!				reactionCurrent%cellNumber(i)=cell
!				reactionCurrent%taskid(i)=myMesh(cell)%proc
!			25 continue
			
			!Find reaction rate for He ion implantation using ImplantReactions(reac), which is input from file.
!			reactionCurrent%reactionRate=findReactionRate(cell, ImplantReactions(matNum,reac))
!			nullify(reactionCurrent%next)
			
!		endif
!********************************************************************
	else
		write(*,*) 'error unknown implantation type'
	endif
10 continue

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
	
	do 14 while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='numProcs') then
			read(87,*) numProcs
			flag=.TRUE.
		endif
	14 continue
	flag=.FALSE.
	
	if(numProcs .NE. myProc%numTasks) then
		write(*,*) 'Error restart file incorrect number of processors'
		if(myProc%taskid==MASTER) read(*,*)
	endif
	
	do 15 while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='numImplantEvents') then
			read(87,*) numImplantEventsReset
			flag=.TRUE.
		endif
	15 continue
	flag=.FALSE.
	
	do 16 while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='numHeImplantEvents') then
			read(87,*) numHeImplantEventsReset
			flag=.TRUE.
		endif
	16 continue
	flag=.FALSE.
	
	do 17 while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='elapsedTime') then
			read(87,*) elapsedTimeReset
			flag=.TRUE.
		endif
	17 continue
	flag=.FALSE.
	
	!Step 2: skip to part of input file that matches this processor
	
	do 19 while(flag .eqv. .FALSE.)
		read(87,*) char
		if(char=='processor') then
			read(87,*) procID
			if(procID==myProc%taskid) then
				flag=.TRUE.
			endif
		endif
	19 continue
	flag=.FALSE.
	
	!Step 3: Read in defect information from file to defect list (coarse mesh),
	!noting that the coordinates in the input file must match the coordinates
	!of the coarse mesh element.
	
	do 20 i=1,numCells
		
		do 21 while(flag .eqv. .FALSE.)
			read(87,*) char
			if(char=='coordinates') then
				read(87,*) (cellCoord(j), j=1,3)
				flag=.TRUE.
			endif
		21 continue
		flag=.FALSE.
		
		do 22 j=1,3
			if(cellCoord(j) .NE. myMesh(i)%coordinates(j)) then
				write(*,*) 'Error cell coordinates do not match in reset file'
				write(*,*) 'CellCoord', cellCoord(j), 'myMeshCoord', myMesh(i)%coordinates(j)
				if(myProc%taskid==MASTER) read(*,*)
			endif
		22 continue
		
		do 23 while(flag .eqv. .FALSE.)
			read(87,*) char
			if(char=='numDefectTypes') then
				read(87,*) numDefectTypes
				flag=.TRUE.
			endif
		23 continue
		flag=.FALSE.
		
		do 24 j=1,numDefectTypes
			read(87,*) (defectTypeReset(k), k=1,numSpecies)
			read(87,*) defectNumReset
			
			defectCurrent=>DefectList(i)
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
				
				do 25 k=1,numSpecies
					defectCurrent%defectType(k)=defectTypeReset(k)
				25 continue
				defectCurrent%num=defectNumReset
				defectCurrent%cellNumber=i
				
			endif
		24 continue
		
	20 continue
	
	close(87)
endif
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
implicit none

integer cell, dir, k, i
type(defect), pointer :: defectCurrent

do 10 cell=1,numCells
	do 11 dir=1,6
		do 12 k=1,myMesh(cell)%numNeighbors(dir)
			if(myMesh(cell)%neighborProcs(dir,k) .NE. myProc%taskid .AND. &
				myMesh(cell)%neighborProcs(dir,k) .NE. -1) then !neighbor cell is in different proc, not free surface (taskid=-1)
				
				if(myMesh(cell)%neighborProcs(dir,k) .NE. myProc%procNeighbor(dir)) then
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
					defectCurrent%defectType(1)=1
					do 13 i=2,numSpecies
						defectCurrent%defectType(i)=0
					13 continue
					defectCurrent%num=CuAtomsEverMesh
					defectCurrent%cellNumber=myMesh(cell)%neighbors(dir,k)
					nullify(defectCurrent%next)
				endif
			endif
		12 continue
	11 continue
10 continue

nullify(defectCurrent)

end subroutine

!***************************************************************************************************
!>Subroutine: Initialize Fine Mesh - initializes defect and reaction lists in a newly created fine mesh
!!
!!Allocates the size of the fine mesh (fine defect list and reaction list included as part of the cascade derived type),
!!populates with defects from the coarse mesh and removes those defects from the coarse mesh.
!!
!!This subroutine includes cascade defect interaction with pre-existing defects.
!!
!!Inputs: CascadeCurrent (cascade type derived type variable)
!!Output: CascadeCurrent (initialized with defects from coarse mesh)
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
endif

!Allocate cascade defect list and reaction list

allocate(CascadeCurrent%localDefects(numCellsCascade))
allocate(CascadeCurrent%reactionList(numCellsCascade))

!For each cell, initialize the reaction list (no reaction or He ion implantation) and defect list (single He initialized with num=0)

do 10 cell=1,numCellsCascade
	
	allocate(CascadeCurrent%localDefects(cell)%defectType(numSpecies))
	nullify(CascadeCurrent%localDefects(cell)%next)
	do 11 j=1,numSpecies
		CascadeCurrent%localDefects(cell)%defectType(j)=0
	11 continue
	
	CascadeCurrent%localDefects(cell)%num=0
	CascadeCurrent%localDefects(cell)%cellNumber=cell

	if(HeDPARatio .GT. 0d0) then
	
		CascadeCurrent%reactionList(cell)%numReactants=0
		CascadeCurrent%reactionList(cell)%numProducts=1
		allocate(CascadeCurrent%reactionList(cell)%products(CascadeCurrent%reactionList(cell)%numProducts, numSpecies))
		
		!search ImplantList for He implant reactions
		do 54  reac=1,numImplantReac(matNum)
			if(ImplantReactions(matNum,reac)%numReactants==0 .AND. ImplantReactions(matNum,Reac)%numProducts==1) then
				exit
			endif
		54 continue
		
		do 55 i=1,ImplantReactions(matNum,reac)%numProducts
			do 56 j=1,numSpecies
				CascadeCurrent%reactionList(cell)%products(i,j)=ImplantReactions(matNum,reac)%products(i,j)
			56 continue
			CascadeCurrent%reactionList(cell)%cellNumber(i)=cell
			CascadeCurrent%reactionList(cell)%taskid(i)=myMesh(cell)%proc
		55 continue
		
		!Find reaction rate for He ion implantation using ImplantReactions(reac), which is input from file.
		CascadeCurrent%reactionList(cell)%reactionRate=findReactionRateFine(CascadeCurrent%cellNumber, ImplantReactions(matNum,reac))
		
		!Add this reaction rate to total rate both for the cascade and for the entire volume element
		totalRate=totalRate+cascadeCurrent%reactionList(cell)%reactionRate
		CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionList(cell)%reactionRate		
		
		nullify(CascadeCurrent%reactionList(cell)%next)
	
	else
		
		CascadeCurrent%reactionList(cell)%numReactants=0
		CascadeCurrent%reactionList(cell)%numProducts=0
		CascadeCurrent%reactionList(cell)%reactionRate=0d0
		nullify(CascadeCurrent%reactionList(cell)%next)
		
	endif
10 continue

!defectCurrentCoarse and defectPrevCoarse will be used to peruse defects in coarse mesh and place in fine mesh
defectCurrentCoarse=>DefectList(CascadeCurrent%cellNumber)
nullify(defectPrevCoarse)

!Use random numbers to assign the defects in the coarse element to the fine mesh

!The volume ratio is the ratio of the fine mesh volume to the coarse mesh volume, used to define 
!the probability that a coarse mesh defect is placed in the fine mesh.
volumeRatio=CascadeElementVol*dble(numCellsCascade)/(myMesh(CascadeCurrent%cellNumber)%volume)

!This loop goes through defects in the coarse mesh and decides whether they should be inserted into
!the fine mesh. If yes, the defects are removed from the coarse mesh. UpdateDefect needs to be kept current
!in this loop.

do 12 while(associated(defectCurrentCoarse))
	
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
		
		if(n .LE. 12) then
			!chooses number of defects that are added to fine mesh from coarse mesh of this type
			do 13 k=0,n-1
				r2=r2+dble(binomial(n,k))*dble(volumeRatio**k)*dble((1d0-volumeRatio)**(n-k))
				if(r2 .GT. r1) then
					exit
				endif
			13 continue
		else
			
			lambda=volumeratio*n
			rstore=0d0
			do 44 k=0,n-1
				
				!***********************************************************************************
				!For k .LE. 12, a poisson distribution can be calculated. For k>12, factorial(k) is 
				!too large and a Gaussian distribution must be assumed.
				!
				!This leads to some errors where k=12 is slightly more common than it should be as
				!a number of defects added from the
				!coarse mesh to the fine mesh, but this problem is almost impossible to overcome
				!considering the size of the numbers that the machine can handle.
				!***********************************************************************************
				
				if(k .LE. 12) then
					!write(*,*) 'poisson distribution'
					r2=r2+lambda**(dble(k))*dexp(-lambda)/factorial(k)
					!write(*,*) 'r2', r2, lambda**(dble(k))*dexp(-lambda)/factorial(k)
					!write(*,*) 'r1', r1, 'r2', r2, 'k', k, 'n', n
					if(r2 .GT. r1) then
						exit
					endif
				else
					!write(*,*) 'gaussian distribution'
					r2=r2+dexp(-(k-lambda)**(2d0)/(2d0*lambda))/(dsqrt(2*pi*lambda))
					!write(*,*) 'r1', r1, 'r2', r2, 'k', k, 'n', n
					if(rstore-r2 == 0d0) then
						!write(*,*) 'gaussian exit'
						exit
					else if(r2 .GT. r1) then
						!write(*,*) 'gaussian exit'
						exit
					endif
					rstore=r2
				endif
			44 continue
			!read(*,*)
			!write(*,*) 'k', k
		endif
		
!		if(k .NE. 0) then
!			write(*,*) 'k not 0', k
!			write(*,*) 'defectType', defectCurrentCoarse%defectType
!			!read(*,*)
!		endif
		
		!***********
		!k is the number of defects being deposited into the fine mesh
		!***********
		
		do 20 num=1,k
			
			!***************************************************************
			!Generate the cell within the fine mesh that the defect will be deposited into
			!***************************************************************
			r1=dprand()*dble(numCellsCascade)
			r2=0d0
			do 21 cell=1,numCellsCascade
				r2=r2+dble(cell)
				if(r2 .GT. r1) then
					exit
				endif
			21 continue
			
			!***************************************************************
			!Deposit the defect into the fine mesh
			!***************************************************************
			do 24 j=1,numSpecies
				products(j)=defectCurrentCoarse%defectType(j)
			24 continue
			!write(*,*) 'inserting into fine mesh', (products(j), j=1,numSpecies), 'k', k
			nullify(defectPrevFine)
			defectCurrentFine=>CascadeCurrent%localDefects(cell)
			!this subroutine will move defectCurrent to the place in the list where reactionCurrent%products(i,x) exists OR, if it doesn't exist,
			! defectPrev will point to the defect before the insertion place and defectCurrent will point to the defect after insertion
			
			call findDefectInList(defectCurrentFine, defectPrevFine, products)
			
			if(associated(defectCurrentFine)) then !if we aren't at the end of the list
				
				count=0
				!Check to see if this defect already exists in the fine mesh list
				do 25 j=1,numSpecies
					if(defectCurrentFine%defectType(j)==products(j)) then
						count=count+1
					endif
				25 continue
				
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
					do 14 j=1,numSpecies
						defectPrevFine%defectType(j)=products(j)
					14 continue
					defectPrevFine%next=>defectCurrentFine !if inserted defect is in the middle of the list, point it to the next item in the list
					
				endif
			else 			!add a defect to the end of the list
				nullify(defectPrevFine%next)
				allocate(defectPrevFine%next)
				nullify(defectPrevFine%next%next)
				defectPrevFine=>defectPrevFine%next
				allocate(defectPrevFine%defectType(numSpecies))
				defectPrevFine%cellNumber=cell
				defectPrevFine%num=1
				do 15 j=1,numSpecies
					defectPrevFine%defectType(j)=products(j)
				15 continue
			endif
			
			!***************************************************************
			!Remove the defect from the coarse mesh
			!***************************************************************
			if(.NOT. associated(defectCurrentCoarse)) then
				write(*,*) 'Tried to delete defect that wasnt there fine mesh initialization'
			
			!if defectCurrentCoarse is in the middle of the list and there is 1 of them
			else if(defectCurrentCoarse%num==1 .AND. associated(defectCurrentCoarse%next) .AND. associated(defectPrevCoarse)) then
				
				defectPrevCoarse%next=>defectCurrentCoarse%next !remove that defect type from the system
				deallocate(defectCurrentCoarse%defectType)
				deallocate(defectCurrentCoarse)
				defectCurrentCoarse=>defectPrevCoarse
			
			!if defectCurrentCoarse is at the end of the list and there is one of them
			else if(defectCurrentCoarse%num==1 .AND. associated(defectPrevCoarse)) then
				deallocate(defectCurrentCoarse%defectType)
				deallocate(defectCurrentCoarse)
				defectCurrentCoarse=>defectPrevCoarse	!remove the last defect from the system
				nullify(defectPrevCoarse%next)
			
			!if defectCurrentCoarse is at the beginning of the list and there is one of them
			else if(defectCurrentCoarse%num==1 .AND. associated(defectCurrentCoarse%next)) then !removing first defect from cell i
				defectCurrentCoarse%num=0 !first defect in system never deallocates, it is single helium. set number equal to zero.
			
			!if there is only one element in the list and there is one of them
			else if(defectCurrentCoarse%num==1) then 							!removing only defect from cell i (single helium) - this is redundant but will keep for now
				defectCurrentCoarse%num=0
			
			!if we are trying to remove more defects than are actually present (the calculator in prev. step chose k incorrectly)
			else if(defectCurrentCoarse%num==0) then
				write(*,*) 'trying to remove a defect that isnt there (fine mesh)'
				write(*,*) 'k', k, 'num', num, 'num defects coarse', n
			else
			
			!If there is more than 1 defect of this type, we don't have to to do anything with the pointers
			!and just subtract 1 from the number of defects in the system.
				defectCurrentCoarse%num=defectCurrentCoarse%num-1 !remove the defect from the system instead of the entire entry in the list
			endif
			
		20 continue
	endif
	
	defectPrevCoarse=>defectCurrentCoarse
	defectCurrentCoarse=>defectCurrentCoarse%next
	
12 continue

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
do 10 while (flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='meshFile') then
		read(81,*) filename
		flag=.TRUE.
	endif
10 continue
flag=.FALSE.

!read in whether mesh is uniform or non-uniform
do 11 while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='meshType') then
		read(81,*) meshType
		flag=.TRUE.
	endif
11 continue
flag=.FALSE.

!read in whether we have a strain field or not
do 14 while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char=='strainField') then
		read(81,*) strainField
		flag=.TRUE.
	endif
14 continue
flag=.FALSE.

if(strainField=='yes') then
	do 15 while(flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='strainFile') then
			read(81,*) strainFileName
			flag=.TRUE.
		endif
	15 continue
	flag=.FALSE.
	
	do 16 while(flag .eqv. .FALSE.)
		read(81,*) char
		if(char=='dipoleFile') then
			read(81,*) dipoleFileName
			flag=.TRUE.
		endif
	16 continue
	flag=.FALSE.
endif

!these subroutines (located in MeshReader.f90) initialize the mesh and connectivity.
if(meshType=='uniform') then
	call readMeshUniform(filename)
else if(meshType=='nonUniform') then
	call readMeshNonUniform(filename)
else
	write(*,*) 'error mesh type unknown'
endif

!Check to see if we are starting from a saved point (from a reset file)
do 12 while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char == 'debugRestart') then
		read(81,*) debugToggle
		flag=.TRUE.
	endif
12 continue
flag=.FALSE.

do 13 while(flag .eqv. .FALSE.)
	read(81,*) char
	if(char == 'debugRestartFile') then
		read(81,*) restartFileName
		flag=.TRUE.
	endif
13 continue
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

do 10 cell=1,numCells
	reactionCurrent=>reactionList(cell)
	
	!remove this reaction rate from totalRate
	totalRate=totalRate-reactionCurrent%reactionRate
	totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
	
	!Change reaction rate to 0 for cascades and He implantation
	
	reactionCurrent%reactionRate=0d0
	
	if(HeDPARatio .GT. 0d0) then
		reactionCurrent=>reactionCurrent%next
		
		!remove this reaction rate from totalRate
		totalRate=totalRate-reactionCurrent%reactionRate
		totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
		
		reactionCurrent%reactionRate=0d0
	endif
	
10 continue

!Fine mesh

!Here, we only zero the first reaction rate if we have He implantation
!(normally, no defects are implanted into the fine mesh except He)

if(HeDPARatio .GT. 0d0) then

	CascadeCurrent=>ActiveCascades
	
	do 11 while(associated(CascadeCurrent))
	
		do 12 cell=1,numCellsCascade
			
			reactionCurrent=>CascadeCurrent%reactionList(cell)
			
			!remove this reaction rate from totalRate
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate=CascadeCurrent%totalRate-reactionCurrent%reactionRate
			
			reactionCurrent%reactionRate=0d0
		
		12 continue
		
		CascadeCurrent=>CascadeCurrent%next
		
	11 continue
	
endif

end subroutine
