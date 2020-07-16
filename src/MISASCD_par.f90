! 2019.05: This program is used to simulate I-V-S
!***************************************************************************************************
!>Main program
!!
!!For outline, see main page of manual. This program uses the kinetic Monte Carlo algorithm
!!in a synchronous parallel implementation to simulate radiation damage accumulation and
!!subsequent annealing in metals.
!!
!!Several debugging options can be triggered by un-commenting them in this program, in the
!!case of difficulty.
!***************************************************************************************************
program MISASCD
use DerivedType			!<variable classes for MISASCD
use mod_constants		!<module containing all global variables
use randdp				!<module for double precision random number generation
implicit none
include 'mpif.h'

type(defect), pointer :: defectCurrent									!<used to find defects
type(defectUpdateTracker), pointer :: defectUpdate, defectUpdateCurrent !<used to update reactions
type(reaction), pointer :: reactionCurrent								!<used to find reactions
type(reaction), pointer :: reactionChoiceList, reactionChoiceCurrent	!<used to create a list of chosen reactions (for one KMC domain per volume element case)
type(reaction), pointer :: reactionTemp
type(cascade), pointer :: CascadeCurrent								!<used to find defects/reactions in fine mesh

double precision  GenerateTimestep, TotalRateCheck, rateSingle, tau
integer status(MPI_STATUS_SIZE), sim, outputCounter, nullSteps
integer cascadeCell, i, j, cell
real time1, time2
logical releaseToggle, impCascadeToggle

integer CascadeCount, TotalCascades !<Used to count the number of cascades present in the simulation

character(12) filename, filename2, filename3, filename4, filename5, filename6

interface
    subroutine chooseImplantReaction(reactionCurrent, CascadeCurrent)
        use DerivedType
        implicit none
        type(reaction), pointer :: reactionCurrent
        type(cascade), pointer :: CascadeCurrent
    end subroutine

	subroutine chooseReaction(reactionCurrent, CascadeCurrent)
		use DerivedType
		implicit none
		type(reaction), pointer :: reactionCurrent
		type(cascade), pointer :: CascadeCurrent
	end subroutine

	subroutine addCascadeExplicit(reactionCurrent)
		use DerivedType
		implicit none
		type(reaction), pointer :: reactionCurrent
	end subroutine
	
	subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
		use DerivedType
		implicit none
		type(reaction), pointer :: reactionCurrent
		type(DefectUpdateTracker), pointer :: defectUpdateCurrent
		type(cascade), pointer :: CascadeCurrent
	end subroutine

	subroutine updateReactionList(defectUpdate)
		use DerivedType
		implicit none
		type(DefectUpdateTracker), pointer :: defectUpdate
	end subroutine
	
	subroutine DEBUGPrintReaction(reactionCurrent)
		use DerivedType
		implicit none
		type(reaction), pointer :: reactionCurrent
	end subroutine
	
	subroutine DEBUGPrintDefectUpdate(defectUpdate)
		use DerivedType
		implicit none
		type(defectUpdateTracker), pointer :: defectUpdate
	end subroutine

	subroutine DEBUGcheckForUnadmissible(reactionCurrent)
		use DerivedType
		implicit none
		type(Reaction), pointer :: reactionCurrent
	end subroutine

	subroutine releaseFineMeshDefects(CascadeCurrent)
		use DerivedType
		implicit none
		type(cascade), pointer :: CascadeCurrent
	end subroutine

	double precision function totalRateCascade(CascadeCurrent)
		use DerivedType
		implicit none
		type(Cascade), pointer :: CascadeCurrent
	end function
end interface

call cpu_time(time1)

open(81, file='parameters.txt',action='read', status='old')

!***********************************************************************
!<Initialize MPI interface
!***********************************************************************
periods = .true.
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, myProc%numtasks, ierr)		!read number of processors
call MPI_DIMS_CREATE(myProc%numtasks,3,dims, ierr)
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm,ierr)
call MPI_CART_GET(comm,3,dims,periods,myProc%coords,ierr)
call MPI_CART_RANK(comm,myProc%coords,myProc%taskid,ierr)

call MPI_CART_SHIFT(comm,0,1,myProc%procNeighbor(2),myProc%procNeighbor(1),ierr) !x
call MPI_CART_SHIFT(comm,1,1,myProc%procNeighbor(4),myProc%procNeighbor(3),ierr) !y
call MPI_CART_SHIFT(comm,2,1,myProc%procNeighbor(6),myProc%procNeighbor(5),ierr) !z

if(myProc%taskid==MASTER) then
    write(*,*) 'proc division', dims
end if

!***********************************************************************
!Initialize input parameters
!***********************************************************************
call ReadInputs()
call initialMesh()			!open Mesh_xx.txt file and carry out parallel mesh initialization routine

!<Create fine mesh connectivity
allocate(cascadeConnectivity(6, numCellsCascade))	!< numCellsCascade=numxFine*numyFine*numzFine
call createCascadeConnectivity()

!***********************************************************************
!For running multiple simulations, initialize the defects, reactions, boundary, etc. and loop here.
!***********************************************************************
!< begin the outer loop
do sim=1,numSims

!Reset the temperature to tempStore at the beginning of each simulation
temperature=tempStore	!< Temperature read in (K)
if(myProc%taskid==MASTER) then
	write(*,*) 'Initializing trial', sim, 'proc', myProc%taskid, 'temperature', temperature
end if

!<Initialize output files
if(myProc%taskid==MASTER) then
	filename(1:6)='rawdat'
	filename2(1:6)='totdat'
	!filename4(1:6)='rxnrat'
	!filename5(1:6)='debugg'
	write(unit=filename(7:8), fmt='(I2)') sim
	write(unit=filename2(7:8), fmt='(I2)') sim
	!write(unit=filename4(7:8), fmt='(I2)') sim
	!write(unit=filename5(7:8), fmt='(I2)') sim
	filename(9:12)='.out'
	filename2(9:12)='.out'
	!filename4(9:12)='.out'
	!filename5(9:12)='.log'
	if(rawdatToggle=='yes') open(82, file=filename, action='write', status='Unknown')
	if(totdatToggle=='yes') open(83, file=filename2, action='write', status='Unknown')
	!open(85, file=filename4, action='write', status='Unknown')
	!open(86, file=filename5, action='write', status='Unknown')
end if
!filename5(1:6)='debugg'
!write(unit=filename5(7:8), fmt='(I2)') myProc%taskid
!filename5(9:12)='.log'
!open(86, file=filename5, action='write', status='Unknown')

!***********************************************************************
!<Initialize random number: initialize the random number generators on each processor with a different seed. This
!<is done by creating random integers on the master processor and sending them to the slaves as random number seeds.
!<Initialize defectList, reactionList, boundaryDefectList, totalRate.
!***********************************************************************
call initializeRandomSeeds()		!<set unique random number seeds in each processor
allocate(defectList(numCells))		!<Create list of defects - array
allocate(reactionList(numCells))	!<Create list of reactions - array
allocate(totalRateVol(numCells))	!<Create array of total rates in each volume element
call initializeVIdefect()			!<initialize point defects in the whole system
call initializeDefectList()			!<initialize defects within myMesh
call initializeBoundaryDefectList()	!<initialize defects on boundary of myMesh (in other procs)
call initializeReactionList()		!<initialize reactions within myMesh
call initializeTotalRate()			!<initialize totalRate and maxRate using reactionList(:)
call initializeDebugRestart()		!<input defects into coarse mesh from restart file (for debugging)

if(myProc%taskid==MASTER) then
	write(*,*) 'CuEverMesh', numCuCell
	write(*,*) 'Initial number of V', initialNumV
	write(*,*) 'Initial number of SIA', initialNumI
	write(*,*) 'Equilibrium concentration of V', ceqV
	write(*,*) 'Equilibrium concentration of SIA', ceqI
	write(*,*) 'firr',firr
	write(*,*) 'totalRate', totalRate
end if

!***********************************************************************
!<Restart from a reset file.
!<If we are restarting from a reset file, then we have to start wtih a nonzero dpa and
!<at a nonzero time. All of the 'old' implant events are tracked in the master processor.
!***********************************************************************
if(debugToggle=='yes') then

	if(myProc%taskid==MASTER) then
        numImpAnn(1) = numImplantEventsReset	!<numImplantEventsReset is read in from debug restart file
	else
        numImpAnn(1) = 0
	end if
	numImpAnn(2)=0
	elapsedTime	= elapsedTimeReset

	!Reset the reactions within this cell and diffusion to neighboring cells in the same processor
	do i=1,numCells
		call resetReactionListSingleCell(i)
	end do

	!Create the boundary defects in this processor and send boundary defects to neighboring processor.
	!Update reactions for diffusion into boundary (into other processors)
	!NOTE: this will crash if the number of cells in a neighboring processor is not the same as in the local processor.
	!Easy fix, but putting off until later.
	do i=1,numCells
		call cascadeUpdateStep(i)
	end do
	
	write(*,*) 'Processor', myProc%taskid, 'Defect lists and reaction lists initialized'
else
    numImpAnn(1)=0     	!<Postprocessing: number of Frenkel pairs / cascades (local)
    numImpAnn(2)=0     	!<Postprocessing: number of annihilation reactions carried out
	elapsedTime=0d0		!<The simulation time that has passed
	
	numTrapSIA=0		!<Postprocessing: number of SIAs trapped on grain boundary
	numTrapV=0			!<Postprocessing: number of vacancies trapped on grain boundary
	numEmitSIA=0		!<Postprocessing: number of SIAs emitted from grain boundary
	numEmitV=0			!<Postprocessing: number of vacancies emitted from grain boundary
end if

!<Initialize totalTime
if(totalDPA > 0d0 .AND. dpaRate > 0d0) then
	totalTime=totalDPA/dpaRate	!simulation time
else
	totalTime = agingTime
end if

!<Initialize Counters
step=0
nullSteps=0		!Record the number of steps in which an empty event was selected
TotalCascades=0
outputCounter=0
nullify(ActiveCascades)		! pointer ActiveCascades=NULL
annealIdentify=.FALSE.		!(.TRUE. if in annealing phase, .FALSE. otherwise) used to determine how to reset reaction rates (should we include implantation or not)

!<Initizlize maxRate and tau. rateTau(1)=maxRate, rateTau(2)=tau
rateTau=0d0

!*********************************************************************************************************************
!*********************************************************************************************************************
!**																		                                            **
!** 			                               Damage Accumulation Loop								                **
!** During this process, Frenkerl Pairs or cascades are continuously implanted into system. 						**
!**																													**
!*********************************************************************************************************************
!*********************************************************************************************************************
!<KMC LOOP:
!<1. Calculate total reaction rate and time increment.
!<2, Choose a reaction
!<3. Carry out reaction (update defect)->note all relevant defects that have been changed
!<4. Update reaction lists->use defects deleted and defects added (reactants and products) to update relevant reactions
!<	(first, delete all reactions associated with the REACTANTS AND PRODUCTS in the cell and all DIFFUSION reactions associated with
!<	reactants in neighboring cells. Then add all reaction rates associated with REACTANTS AND PRODUCTS in the cell and diffusion
!<	reactions in neighboring cells.)
!<5. Repeat 1-4.

do while(elapsedTime < totalTime)
!do while(step < 1)
	
	step=step+1
        
	!<Logical variable tells us whether cascade communication step needs to be carried out
	!<(0=no cascade, nonzero=number of volume element where cascade event has happened)
	cascadeCell=0
	impCascadeToggle= .FALSE.
	

	!<find the maximum totalRate in all processors
	call MPI_REDUCE(totalRate,maxRate,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm,ierr)

	if(myProc%taskid==MASTER) then
		rateTau(1)=maxRate
	end if

	!<Generate timestep in the master processor and send it to all other processors
	!<If implantScheme=='explicit', cascade implantation needs to be carried out explicitly
	if(implantScheme=='explicit') then
		if(elapsedTime >= numImpAnn(1)*(numDisplacedAtoms*atomSize)/(totalVolume*dpaRate)) then
			!<Do not generate a timestep in this case; this is an explicit (zero-time) reaction
			if(myProc%taskid==MASTER) then
				rateTau(2)=0d0
			end if
		else
			if(myProc%taskid==MASTER) then
				rateTau(2)=GenerateTimestep()
			end if
		end if
	else if(implantScheme=='MonteCarlo') then
		if(myProc%taskid==MASTER) then
			rateTau(2)=GenerateTimestep()
		end if
	end if

	call MPI_BCAST(rateTau, 2, MPI_DOUBLE_PRECISION, MASTER, comm,ierr)

	maxRate=rateTau(1)	!<update maxRate of the processor
	tau=rateTau(2)		!<update time step of the processor

!*************************************************************************
!	if(mod(step,10)==0) then
!		!Debugging subroutine: outputs the reaction rate in each processor
!		call outputRates(elapsedTime, step)
!	endif
!*************************************************************************
	
	!***********************************************************************************************
	!Choose from reactions in reactionList(:) (local to this processor). Null events possible.
	!If explicit implant scheme has been chosen, implant cascades when elapsed time has passed 
	!a threshold, and do not iterate timestep when doing so.
	!***********************************************************************************************

	!<Initialize defectUpdate. defectUpdate is used to store updated defects
	allocate(defectUpdate)
	allocate(defectUpdate%defectType(numSpecies))
	do i=1,numSpecies
		defectUpdate%defectType(i)=0
	end do
	nullify(defectUpdate%next)
	defectUpdateCurrent=>defectUpdate

	!<choose a reaction in one volume element
	!<If implantScheme=='explicit', cascade implantation needs to be carried out explicitly
	if(implantScheme=='explicit') then
		if(elapsedTime >= numImpAnn(1)*(numDisplacedAtoms*atomSize)/(totalVolume*dpaRate)) then
			!<choose a cell in the peocessor to implant cascade
			call addCascadeExplicit(reactionCurrent)
		else
			call chooseReaction(reactionCurrent, CascadeCurrent)
		end if
	else if(implantScheme=='MonteCarlo') then
		if(test3 == 'yes') then
			call initializeOneCascade()
			call chooseImplantReaction(reactionCurrent, CascadeCurrent)
		else
			call chooseReaction(reactionCurrent, CascadeCurrent)
		end if

	else
		write(*,*) 'error choosing reaction main program'
	end if
		
	!***********************************************************************************************
	!Update defects according to reactions chosen. Communicate defects that have changed on
	!boundaries of other processors and create a list of defects whose reaction rates must be update
	!in this processor due to defects updated in this processor and in neighboring processors.
	!***********************************************************************************************
	if(.NOT. associated(reactionCurrent)) then
		nullSteps=nullSteps+1
	end if

	call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)

	!<If a cascade is chosen, update reaction rates for all defects remaining in coarse mesh element
	if(associated(reactionCurrent)) then
		if(reactionCurrent%numReactants==-10) then

			!Resets reaction list in cell and neighbors (on same processor)
			call resetReactionListSingleCell(reactionCurrent%cellNumber(1))

			!Variable tells us whether cascade communication step needs to be carried out
			!and if so in what volume element in the coarse mesh
			cascadeCell=reactionCurrent%cellNumber(1)
			impCascadeToggle=.TRUE.
		end if
	end if

	!<Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
	elapsedTime=elapsedTime+tau

	call updateReactionList(defectUpdate)

	!if(totalRate < 0d0) then
	!	write(*,*) 'error totalRate less than zero', step
		!call DEBUGPrintDefectList()
	!	call DEBUGPrintReactionList()

		!call MPI_ABORT(comm,ierr)
	!end if

	!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
	!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
	!release the defects into the coarse mesh.
	releaseToggle = .FALSE.
	if(associated(CascadeCurrent)) then
		if(totalRateCascade(CascadeCurrent) < cascadeReactionLimit) then
			
			!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
			!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
			cascadeCell=CascadeCurrent%cellNumber
			releaseToggle=.TRUE.

			!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
			call releaseFineMeshDefects(CascadeCurrent)
			call resetReactionListSingleCell(cascadeCell)
		end if
	end if
	
	!Cascade communication step:
	!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
	!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.
	if(implantType=='Cascade') then
		!call cascadeUpdateStep(cascadeCell)
		call cascadeUpdateStep(releaseToggle, cascadeCell)
	end if

	!****************************************************************************************
	!Every time a cascade is added, we reset the total rate and check how wrong it has become
	!!'Cascade' implant type, check total rate when implantation event is choosen
	!!'FrenkelPair' implant type, check total rate every 1000 steps
	!*****************************************************************************************
	if(associated(reactionCurrent)) then
		if(implantType=='Cascade') then
			!if(reactionCurrent%numReactants==-10) then
			if(impCascadeToggle .EQV. .TRUE.) then	!implanted a cascade
				totalRate=TotalRateCheck()
			end if
		else if(implantType=='FrenkelPair') then
			if(mod(step,1000)==0) then
				totalRate=TotalRateCheck()
			end if
		end if
	end if

	!******************************************
	! Optional: count how many cascades are present at step i and compile to find avg. number of cascades present per step
	!******************************************
	TotalCascades=TotalCascades+CascadeCount()

	!********************************************************************************
	! Output according to outputCounter
	!********************************************************************************
	if(elapsedTime >= totalTime/2.0d6*(2.0d0)**(outputCounter)) then
    !if(elapsedTime >= totalTime/50d0*(outputCounter)) then
	! or if(mod(step,100000)==0) then
		call MPI_REDUCE(numImpAnn,totalImpAnn, 2, MPI_INTEGER, MPI_SUM, 0,comm, ierr)

		if(myProc%taskid==MASTER) then
			DPA=dble(totalImpAnn(1))/(systemVol/(numDisplacedAtoms*atomSize))
			call cpu_time(time2)
			write(*,*)
			write(*,*) 'time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)
			write(83,*) '*********************************************************************************************'
			write(83,*) 'time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)

			if(implantType=='FrenkelPair') then
				write(*,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
				write(83,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
			else if(implantType=='Cascade')	then
				write(*,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
				write(83,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
			else	!Thermal aging
				write(*,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
				write(83,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
			end if

			!Optional: output fraction of steps that are null events
			if(singleElemKMC=='yes') then
				write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
			else
				write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step)
			end if
			write(*,*)
		end if
		
		!Several defect output optionas available.
		if(rawdatToggle=='yes') call outputDefectsXYZ()		!write(82,*): rawdat
		if(totdatToggle=='yes') call outputDefectsTotal()	!write(83,*): totdat.out
		!if(xyzToggle=='yes') call outputDefectsXYZ()	!write(87,*): defect.xyz
		if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)	!write(88,*): VTKout.vtk
		if(outputDebug=='yes') call outputDebugRestart(outputCounter)	!write(88,*): Restart.in

		outputCounter=outputCounter+1

		call MPI_BARRIER(comm,ierr)
	end if

	if(test3 == 'yes') then
		exit
	end if

end do

!***********************************************************************
!Output defects at the end of the implantation loop
!***********************************************************************
call MPI_REDUCE(numImpAnn, totalImpAnn, 2, MPI_INTEGER, MPI_SUM,0, comm, ierr)

if(myProc%taskid==MASTER) then
	DPA=dble(totalImpAnn(1))/(systemVol/(numDisplacedAtoms*atomSize))
	call cpu_time(time2)
	write(*,*)
	write(*,*) 'Final  step'
	write(*,*) 'time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)
	write(83,*) '*********************************************************************************************'
	write(83,*) 'elapsedTime', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)

	if(implantType=='FrenkelPair') then
		write(*,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
		write(83,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
	else if(implantType=='Cascade')	then
		write(*,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
		write(83,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
	else
		write(*,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
		write(83,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
	end if
	write(*,*)
	write(83,*) 'Final  step'
end if

!Final output
if(totdatToggle=='yes') call outputDefectsTotal()
if(rawdatToggle=='yes') call outputDefectsXYZ()
if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
if(outputDebug=='yes') call outputDebugRestart(outputCounter)

!***********************************************************************************************************************
!***********************************************************************************************************************
!**																		                                              **
!** 							              Annealing Loop								     					  **
!**																		   											  **
!***********************************************************************************************************************
!***********************************************************************************************************************

!***************************************************************************
!Annealing loop: change the temperature to annealTemp and change all
!defect production rates to 0 (dpa rate, He implant rate)
!Then carry out simulation normally until elapsedTime=totalTime+annealTime
!***************************************************************************
totalTime=0d0
if(annealTime > 0d0) then
	call annealInitialization()	!Sets up simulation for annealing (change temp, etc)
	annealIdentify=.TRUE.

	outputCounter=0		!Used to output once for each annealing step
	elapsedTime=0d0
	totalTime=annealTime
	step=0
	nullSteps=0

	if(myProc%taskid==MASTER ) then
		write(*,*) '************************************'
		write(*,*) 'Entering Annealing Phase'
		write(83,*) '************************************'
		write(83,*) 'Entering Annealing Phase'
	end if
end if

!begin annealing
!carry out annealing in multiple steps if desired. Each step increases temperature.
!!do annealIter=1,annealSteps	!default value: annealSteps = 1

	!!if(annealType=='mult') then
	!!	temperature=annealTemp*annealTempInc**dble(annealIter-1)
	!!else if(annealType=='add') then
	!!	temperature=annealTemp+dble(annealIter-1)*annealTempInc
	!!else
	!!	write(*,*) 'error unknown anneal type'
	!!end if

	!Possible bug: this will make the reaction rates for implantation non-zero. Not sure why this wasn't a problem before.
	!do i=1,numCells
	!	call resetReactionListSingleCell(i)
	!end do

	totalRate=TotalRateCheck()

	!if annealTime = 0, the following "do" does not execute
	!!do while(elapsedTime < totalTime+dble(annealIter)*annealTime/dble(annealSteps))
	do while(elapsedTime < totalTime)

		step=step+1
	
		!Logical variable tells us whether cascade communication step needs to be carried out
		!(0=no cascade, nonzero=number of volume element where cascade event has happened)
		cascadeCell=0

		call MPI_REDUCE(totalRate,maxRate,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm,ierr)

		if(myProc%taskid==MASTER) then
			rateTau(1)=maxRate
		end if

		if(myProc%taskid==MASTER) then
			rateTau(2)=GenerateTimestep()
			!if(elapsedTime-totalTime+tau > annealTime/dble(annealSteps)*outputCounter) then
			!we have taken a timestep that moves us past this annealing step
			!	rateTau(2)=annealTime/dble(annealSteps)*outputCounter-(elapsedTime-totalTime)
			!end if
		end if

		call MPI_BCAST(rateTau, 2, MPI_DOUBLE_PRECISION, MASTER, comm,ierr)

		maxRate=rateTau(1)	!update maxRate
		tau=rateTau(2)		!update time step

!*************************************************************************
!	if(mod(step,10)==0) then
!		!Debugging subroutine: outputs the reaction rate in each processor
!		call outputRates(step)
!	endif
!*************************************************************************
	
		!***********************************************************************************************
		!Choose from reactions in in reactionList(:) (local to this processor). Null events possible.
		!***********************************************************************************************
		allocate(defectUpdate)
		allocate(defectUpdate%defectType(numSpecies))
		do i=1,numSpecies
			defectUpdate%defectType(i)=0
		end do
		nullify(defectUpdate%next)
		defectUpdateCurrent=>defectUpdate
	

		!choose a reaction in one volume element
		call chooseReaction(reactionCurrent, CascadeCurrent)

		!***********************************************************************************************
		!Update defects according to reactions chosen. Communicate defects that have changed on
		!boundaries of other processors and create a list of defects whose reaction rates must be update
		!in this processor due to defects updated in this processor and in neighboring processors.
		!***********************************************************************************************

		!************
		! Optional: count how many steps are null events
		!************
		if(.NOT. associated(reactionCurrent)) then
			nullSteps=nullSteps+1
		end if
	
		call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
		
		if(associated(reactionCurrent)) then
			if(reactionCurrent%numReactants==-10) then
				write(*,*) 'Error chose casacde implantation in annealing phase'
			end if
		end if

		!Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
		elapsedTime=elapsedTime+tau

		!Update reaction rates for defects involved in reaction chosen
		call updateReactionList(defectUpdate)

		!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
		!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
		!release the defects into the coarse mesh.
		releaseToggle=.FALSE.
		if(associated(CascadeCurrent)) then
			if(totalRateCascade(CascadeCurrent) < cascadeReactionLimit) then
			
				!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
				!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
				cascadeCell=CascadeCurrent%cellNumber
				releaseToggle=.TRUE.

				!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
				call releaseFineMeshDefects(CascadeCurrent)
				call resetReactionListSingleCell(cascadeCell)

				write(*,*) 'taskid', myProc%taskid, 'step', step, 'anneal step', annealIter, &
						'totalRateCascade',totalRateCascade(CascadeCurrent),'release fine meshes'
			end if
		end if
		
		!Cascade communication step:
		!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
		!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.
		call cascadeUpdateStep(releaseToggle,cascadeCell)
		
		!Update the totalRate in order to avoid any truncation error every 1000 steps (this value can be modified)
		if(mod(step,1000)==0) then
			totalRate=TotalRateCheck()
		end if
	
		!******************************************
		! Optional: count how many cascades are present at step i and compile to find avg. number
		! of cascades present per step
		!******************************************
		TotalCascades=TotalCascades+CascadeCount()
	
		!******************************************
		! Output
		!******************************************
		!!if((elapsedTime-totalTime) >= annealTime/dble(annealSteps)*outputCounter) then
		if(elapsedTime >= totalTime/2.0d6*(2.0d0)**(outputCounter)) then

			call MPI_REDUCE(numImpAnn,totalImpAnn, 2, MPI_INTEGER, MPI_SUM,0,comm, ierr)

			if(myProc%taskid==MASTER) then
				DPA=dble(totalImpAnn(1))/(systemVol/(numDisplacedAtoms*atomSize))
				call cpu_time(time2)
				write(*,*)
				write(*,*) 'Anneal time', elapsedTime, 'DPA',DPA, 'steps',step, 'AverageTimeStep', elapsedTime/dble(step)
				!!write(*,*) 'Anneal step', annealIter, 'temperature', temperature
				write(83,*) '*********************************************************************************************'
				write(83,*) 'Anneal time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)
				!!write(83,*) 'Anneal step', annealIter, 'temperature', temperature

				if(implantType=='FrenkelPair') then
					write(*,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
					write(83,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
				else if(implantType=='Cascade')	then
					write(*,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
					write(83,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
				else	!Thermal aging
					write(*,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
					write(83,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
				end if

				!Optional: output fraction of steps that are null events
				if(singleElemKMC=='yes') then
					write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
				else
					write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step)
				end if
				write(83,*)
				write(*,*)
			end if

			!Several defect output optionas available.
			if(totdatToggle=='yes') call outputDefectsTotal()	!write(83,*): totdat.out
			if(rawdatToggle=='yes') call outputDefectsXYZ()		!write(82,*): rawdat
			!if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)		!write(88,*): VTKout.vtk
			if(outputDebug=='yes') call outputDebugRestart(outputCounter)	!write(88,*): Restart.in
		
			outputCounter=outputCounter+1
			!call MPI_BARRIER(comm，ierr)
		end if
	end do	!Anneal time loop

!!end do	!Multiple anneal steps loop

!***********************************************************************
!Output final defect state
!***********************************************************************
if(annealTime > 0d0) then

	call MPI_REDUCE(numImpAnn,totalImpAnn, 2, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

	write(*,*) 'Fraction null steps', dble(nullSteps)/dble(step), 'Proc', myProc%taskid

	if(myProc%taskid==MASTER) then
		DPA=dble(totalImpAnn(1))/(systemVol/(numDisplacedAtoms*atomSize))
		call cpu_time(time2)

		write(*,*) 'Final Defect State'
		write(83,*) 'Final Defect State'

		write(*,*) 'Anneal time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)
		!!write(*,*) 'Anneal step', annealIter, 'temperature', temperature
		write(83,*) '*********************************************************************************************'
		write(83,*) 'Anneal time', elapsedTime, 'DPA', DPA, 'steps', step, 'AverageTimeStep', elapsedTime/dble(step)
		!!write(83,*) 'Anneal step', annealIter, 'temperature', temperature

		if(implantType=='FrenkelPair') then
			write(*,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
			write(83,*) 'FrenkelPairs', totalImpAnn(1), 'computationTime', time2-time1
		else if(implantType=='Cascade')	then
			write(*,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
			write(83,*) 'Cascades', totalImpAnn(1), 'computationTime', time2-time1
		else
			write(*,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
			write(83,*) 'noImplantation', totalImpAnn(1), 'computationTime', time2-time1
		end if
		write(*,*)
		write(83,*) 'Final  step'


		!Optional: output fraction of steps that are null events
		if(singleElemKMC=='yes') then
			write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
		else
			write(83,*) 'Fraction null steps', dble(nullSteps)/dble(step)
		end if

	end if

	!Final output
	if(totdatToggle=='yes') call outputDefectsTotal()
	if(rawdatToggle=='yes') call outputDefectsXYZ()
	if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
	if(outputDebug=='yes') call outputDebugRestart(outputCounter)

end if

!Defect profile is only output during last output step (not at intermediate steps)
!if(profileToggle=='yes') call outputDefectsProfile(sim)

!***********************************************************************
!Final step: release all cascades into coarse mesh
!(deallocate all extra memory here)
!***********************************************************************
call cpu_time(time2)

if(myProc%taskid==MASTER) then
	write(*,*) 'computation time', time2-time1
	write(*,*) 'total steps', step
	write(83,*) 'computation time', time2-time1
	write(83,*) 'total steps', step
	write(*,*) 'Deallocating memory: fine mesh defects and reactions'
end if

do while(associated(ActiveCascades))
	CascadeCurrent=>ActiveCascades

	!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
	!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
	cascadeCell=CascadeCurrent%cellNumber

	!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
	call releaseFineMeshDefects(CascadeCurrent)
	call resetReactionListSingleCell(cascadeCell)
end do

if(myProc%taskid==MASTER) then
	!write(82,*) 'Released all fine mesh defects'
	write(83,*) 'Released all fine mesh defects'
	!write(84,*) 'Released all fine mesh defects'
end if

!call outputDefectsXYZ()
!call outputDefectsTotal(elapsedTime, step)

!Final memory cleanup: deallocate defects and reactions in coarse mesh
if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating memory: coarse mesh defects and reactions'
end if

call deallocateDefectList()
call deallocateReactionList()
call deallocateBoundarydefectList()
deallocate(totalRateVol)
if(allocated(listVI)) then
	deallocate(listVI)
end if

if(myProc%taskid==MASTER) then
	close(82)
	close(83)
	!close(84)
	!close(85)
	!close(86)
	!close(87)
end if

!End of loop for multiple trials
end do

!***********************************************************************
!Last step: deallocating memory
!***********************************************************************
deallocate(cascadeConnectivity)

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating memory: cascade list'
end if

call deallocateCascadeList()

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating material input data'
end if

call deallocateMaterialInput()

close(81)

call MPI_FINALIZE(ierr)		!must be performed at end of simulation to have no error

end program
