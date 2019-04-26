! $Header: /home/CVS//srscd/src/SRSCD_par.f90,v 1.34 2015/12/14 21:34:49 aydunn Exp $

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

program SRSCD
use DerivedType			!variable classes for SRSCD
use MeshReader			!module created for reading in mesh
use mod_srscd_constants		!module containing all global variables
use randdp				!module for double precision random number generation
implicit none
include 'mpif.h'

type(defect), pointer :: defectCurrent									!used to find defects
type(defectUpdateTracker), pointer :: defectUpdate, defectUpdateCurrent !used to update reactions
type(reaction), pointer :: reactionCurrent								!used to find reactions
type(reaction), pointer :: reactionChoiceList, reactionChoiceCurrent	!used to create a list of chosen reactions (for one KMC domain per volume element case)
type(cascade), pointer :: CascadeCurrent								!used to find defects/reactions in fine mesh

double precision elapsedTime, totalTime, tau, GenerateTimestep, TotalRateCheck, rateSingle
integer status(MPI_STATUS_SIZE), step, annealIter, sim, numDefectsRecv, tracker, outputCounter, nullSteps
integer cascadeCell, i, j, k, cell
integer, allocatable :: cellRecvTot(:), defectRecvTot(:,:)
real time1, time2

integer CascadeCount, TotalCascades !<Used to count the number of cascades present in the simulation

character(12) filename, filename2, filename3, filename4, filename5, filename6

double precision rateDiff	!temporary
type(Reaction), pointer :: reactionTemp

!***********************************************************************
!7.2.2015 Adding an iterative search for sink efficiency. Variables below:
!***********************************************************************

double precision, allocatable :: conc_v_store(:), conc_i_store(:)
logical alpha_v_search, alpha_i_search, searchToggle
double precision conc_v_avg, conc_v_stdev, conc_i_avg, conc_i_stdev
double precision computeVConc, computeIConc	!subroutines in postprocessing file
double precision alpha_v_prev, alpha_i_prev
double precision conc_v_prev, conc_i_prev
double precision alpha_temp
integer v_iteration, i_iteration
double precision Residual, Residual_stdev, Residual_prev
double precision Residual_square, Residual_sqrdev

interface
	subroutine chooseReaction(reactionCurrent, CascadeCurrent)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	type(cascade), pointer :: CascadeCurrent
	end subroutine
	
	subroutine chooseReactionSingleCell(reactionCurrent, CascadeCurrent, cell)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	type(cascade), pointer :: CascadeCurrent
	integer cell
	end subroutine
	
	subroutine addCascadeExplicit(reactionCurrent)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	end subroutine
	
	subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	type(DefectUpdateTracker), pointer :: defectUpdateCurrent
	type(cascade), pointer :: CascadeCurrent
	end subroutine
	
	subroutine updateDefectListMultiple(reactionChoiceCurrent, defectUpdateCurrent, CascadeCurrent)
	use DerivedType
	type(reaction), pointer :: reactionChoiceCurrent
	type(defectUpdateTracker), pointer :: defectUpdateCurrent
	type(cascade), pointer :: CascadeCurrent	
	end subroutine
	
	subroutine updateReactionList(defectUpdate)
	use DerivedType
	type(DefectUpdateTracker), pointer :: defectUpdate
	end subroutine
	
	subroutine debugPrintReaction(reactionCurrent, step)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	integer step
	end subroutine
	
	subroutine debugPrintDefectUpdate(defectUpdate)
	use DerivedType
	type(defectUpdateTracker), pointer :: defectUpdate
	end subroutine
	
	subroutine releaseFineMeshDefects(CascadeCurrent)
	use DerivedType
	type(cascade), pointer :: CascadeCurrent
	end subroutine
	
	subroutine DEBUGcheckForUnadmissible(reactionCurrent, step)
	use DerivedType
	type(Reaction), pointer :: reactionCurrent
	integer step
	end subroutine
	
	double precision function totalRateCascade(CascadeCurrent)
	use DerivedType
	type(Cascade), pointer :: CascadeCurrent
	end function
end interface

call cpu_time(time1)

open(81, file='parameters.txt',action='read', status='old')

!Initialize MPI interface

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, myProc%numtasks, ierr)		!read number of processors
call MPI_COMM_RANK(MPI_COMM_WORLD, myProc%taskid, ierr)			!read processor ID of this processor

call initializeMesh()			!open mesh file and carry out parallel mesh initialization routine
call selectMaterialInputs()		!open material input files and read all relevant data (migration and binding energies, etc)
call readCascadeList()			!input cascade list from file to cascadeList (global variable)
call readImplantData()			!input (1-dimensional) non-uniform defect implantation info (DPA rates, He implant rates)
call readParameters()			!read simulation parameters (DPA rate, temperature, etc)

!Create fine mesh connectivity
allocate(cascadeConnectivity(numCellsCascade, 6))
call createCascadeConnectivity()

!***********************************************************************
!7.2.2015 Creating iterative search for sink efficiency, looping here
!***********************************************************************
searchToggle=.TRUE.

!Begin by searching for alpha_i, then switch to alpha_v
alpha_v_search=.FALSE.
alpha_i_search=.TRUE.

allocate(conc_i_store(numSims))
allocate(conc_v_store(numSims))

!Initialize values of alpha_v and alpha_i
if(sinkEffSearch=='yes') then
	conc_i_prev=0d0
	conc_v_prev=0d0
	
	alpha_i_prev=0d0
	alpha_v_prev=0d0
	
	v_iteration=0
	i_iteration=0
endif

do 500 while(searchToggle .eqv. .TRUE.)
!do 501 while(alpha_v_search .eqv. .TRUE.)
!do 502 while(alpha_i_search .eqv. .TRUE.)

!***********************************************************************
!For running multiple simulations, initialize the defects, reactions,
!boundary, etc. and loop here.
!***********************************************************************

do 13 sim=1,numSims

!Reset the temperature to tempStore at the beginning of each simulation
temperature=tempStore

!Initialize output files
if(myProc%taskid==MASTER) then
	filename(1:6)='rawdat'
	filename2(1:6)='totdat'
	filename3(1:6)='postpr'
	!filename4(1:6)='rxnrat'
	!filename5(1:6)='debugg'
	filename6(1:6)='defect'
	write(unit=filename(7:8), fmt='(I2)') sim
	write(unit=filename2(7:8), fmt='(I2)') sim
	write(unit=filename3(7:8), fmt='(I2)') sim
	!write(unit=filename4(7:8), fmt='(I2)') sim
	!write(unit=filename5(7:8), fmt='(I2)') sim
	write(unit=filename6(7:8), fmt='(I2)') sim
	filename(9:12)='.out'
	filename2(9:12)='.out'
	filename3(9:12)='.out'
	!filename4(9:12)='.out'
	!filename5(9:12)='.log'
	filename6(9:12)='.xyz'
	if(rawdatToggle=='yes') open(82, file=filename, action='write', status='Unknown')
	if(totdatToggle=='yes') open(83, file=filename2, action='write', status='Unknown')
	if(postprToggle=='yes') open(84, file=filename3, action='write', status='Unknown')
	!open(85, file=filename4, action='write', status='Unknown')
	!open(86, file=filename5, action='write', status='Unknown')
	if(xyzToggle=='yes') open(87, file=filename6, action='write', status='Unknown')
endif

!filename5(1:6)='debugg'
!write(unit=filename5(7:8), fmt='(I2)') myProc%taskid
!filename5(9:12)='.log'
!open(86, file=filename5, action='write', status='Unknown')

if(myProc%taskid==MASTER) then
	write(*,*) 'Initializing trial', sim, 'proc', myProc%taskid, 'temp', temperature
endif

!NEXT STEPS:
!3. Choose Reaction
!4. Carry out reaction (update defect)->note all relevant defects that have been changed
!5. Update reaction lists->use defects deleted and defects added (reactants and products) to update relevant reactions
!	(first, delete all reactions associated with the REACTANTS AND PRODUCTS in the cell and all DIFFUSION reactions associated with
!	reactants in neighboring cells. Then add all reaction rates associated with REACTANTS AND PRODUCTS in the cell and diffusion 
!	reactions in neighboring cells.)
!6. Calculate total reaction rate and communicate using MPI_GATHER the total reaction rate.
!7. Repeat 2-6

!*******
!Initializing reaction lists and defect lists
!*******

!Here we initialize the random number generators on each processor with a different seed. This 
!is done by creating random integers on the master processor and sending them to the slaves as
!random number seeds.

call initializeRandomSeeds()		!set unique random number seeds in each processor
allocate(DefectList(numCells))		!Create list of defects - array
allocate(reactionList(numCells))	!Create list of reactions - array
allocate(totalRateVol(numCells))	!Create array of total rates in each volume element
call initializeDefectList()			!initialize defects within myMesh
call initializeBoundaryDefectList()	!initialize defects on boundary of myMesh (in other procs)
call initializeReactionList()		!initialize reactions within myMesh
call initializeTotalRate()			!initialize totalRate and maxRate using reactionList(:)
call initializeDebugRestart()		!input defects into coarse mesh from restart file (for debugging)

!call DEBUGPrintReactionList(0)

!******************************************************************
!Initialize Counters
!******************************************************************

if(debugToggle=='yes') then
	
	!If we are restarting from a reset file, then we have to start wtih a
	!nonzero dpa and at a nonzero time. All of the 'old' implant events
	!and He implant events are tracked in the master processor.
	
	if(myProc%taskid==MASTER) then
		numImplantEvents	= numImplantEventsReset
		numHeImplantEvents	= numHeImplantEventsReset
	else
		numImplantEvents	= 0
		numHeImplantEvents	= 0
	endif
	elapsedTime			= elapsedTimeReset
	!Reset the reactions within this cell and diffusion to neighboring
	!cells in the same processor
	do 14 i=1,numCells
		call resetReactionListSingleCell(i)
	14 continue

	!Create the boundary defects in this processor and send boundary
	!defects to neighboring processor. Update reactions for diffusion
	!into boundary (into other processors)
	
	!NOTE: this will crash if the number of cells in a neighboring processor
	!is not the same as in the local processor. Easy fix, but putting
	!off until later.
	
	do 15 i=1,numCells
		call cascadeUpdateStep(i)
	15 continue
	
	write(*,*) 'Processor', myProc%taskid, 'Defect lists and reaction lists initialized'
	
else
	numImplantEvents	= 0
	numHeImplantEvents	= 0
	numAnnihilate		= 0
	elapsedTime			= 0d0
	
	numTrapSIA=0
	numTrapV=0
	numEmitSIA=0
	numEmitV=0
endif

step=0
nullSteps=0

totalTime=totalDPA/DPARate
TotalCascades=0
outputCounter=0
nullify(ActiveCascades)
annealIdentify=.FALSE.

do 10 while(elapsedTime .LT. totalTime)
	
	step=step+1
	
	!Logical variable tells us whether cascade communication step needs to be carried out
	!(0=no cascade, nonzero=number of volume element where cascade event has happened)
	cascadeCell=0
	
	!If we are choosing one reaction per volume element, find the max reaction rate in the volume elements
	!(we don't need to do this otherwise, totalRate is automatically updated each time a reaction is carried out)
	if(singleElemKMC=='yes') then
		if(meshingType=='adaptive') then
			write(*,*) 'Error adaptive meshing not allowed for single element kMC'
		endif
		rateSingle=0d0
		do 20 cell=1,numCells
			if(totalRateVol(cell) .GT. rateSingle) then
				rateSingle=totalRateVol(cell)
			endif
		20 continue
		totalRate=rateSingle
	endif
	
	call MPI_ALLREDUCE(totalRate,maxRate,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

!*************************************************************************
!	if(mod(step,10)==0) then
!		!Debugging subroutine: outputs the reaction rate in each processor
!		call outputRates(step)
!	endif
!*************************************************************************
	
	!***********************************************************************************************
	!Choose from reactions in in reactionList(:) (local to this processor). Null events possible.
	!
	!If explicit implant scheme has been chosen, implant cascades when elapsed time has passed 
	!a threshold, and do not iterate timestep when doing so.
	!***********************************************************************************************
	
	allocate(defectUpdate)
	allocate(defectUpdate%defectType(numSpecies))
	do 1 i=1,numSpecies
		defectUpdate%defectType(i)=0
	1 continue
	nullify(defectUpdate%next)
	defectUpdateCurrent=>defectUpdate
	
	if(singleElemKMC=='yes') then	!choose a reaction in each volume element
	
		if(implantScheme=='explicit') then
		
			write(*,*) 'Error explicit implantation not implemented for single element kMC'
			
		else
			!Generate timestep in the master processor and send it to all other processors
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
			endif
			
			allocate(reactionChoiceList)
			nullify(reactionChoiceList%next)
			reactionChoiceCurrent=>reactionChoiceList
					
			!Choose one reaction in each cell
			do 21 cell=1,numCells
				!Choose reactions here
				call chooseReactionSingleCell(reactionCurrent, CascadeCurrent, cell)

				!Update defects according to the reaction chosen
				!write(*,*) 'cell', cell, 'step', step	
				!call DEBUGPrintReaction(reactionCurrent, step)
			
				!************
				! Optional: count how many steps are null events
				!************
			
				if(.NOT. associated(reactionCurrent)) then
					nullSteps=nullSteps+1
				else
					!Generate list of chosen reactions
					
					allocate(reactionChoiceCurrent%next)
					reactionChoiceCurrent=>reactionChoiceCurrent%next
					
					reactionChoiceCurrent%numReactants=reactionCurrent%numReactants
					allocate(reactionChoiceCurrent%reactants(reactionCurrent%numReactants, numSpecies))
					
					reactionChoiceCurrent%numProducts=reactionCurrent%numProducts
					allocate(reactionChoiceCurrent%products(reactionCurrent%numProducts, numSpecies))
					
					allocate(reactionChoiceCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionChoiceCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					
					reactionChoiceCurrent%reactionRate=reactionCurrent%reactionRate
					
					do 31 i=1,reactionCurrent%numReactants
					
						do 41 j=1,numSpecies
							reactionChoiceCurrent%reactants(i,j)=reactionCurrent%reactants(i,j)
						41 continue
					
						reactionChoiceCurrent%cellNumber(i)=reactionCurrent%cellNumber(i)
						reactionChoiceCurrent%taskid(i)=reactionCurrent%taskid(i)
					
					31 continue
					
					do 32 i=1,reactionCurrent%numProducts
						do 42 j=1,numSpecies
							reactionChoiceCurrent%products(i,j)=reactionCurrent%products(i,j)
						42 continue
						
						reactionChoiceCurrent%cellNumber(i+reactionCurrent%numReactants) = &
							reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							
						reactionChoiceCurrent%taskid(i+reactionCurrent%numReactants) = &
							reactionCurrent%taskid(i+reactionCurrent%numReactants)
							
					32 continue
					
					nullify(reactionChoiceCurrent%next)
					
				endif
			
			!	call DEBUGPrintDefectUpdate(defectUpdate)
			
				!No cascade choices allowed for one KMC domain per element at this point
			
			21 continue
			
			reactionChoiceCurrent=>reactionChoiceList%next
			call updateDefectListMultiple(reactionChoiceCurrent, defectUpdateCurrent, CascadeCurrent)
			
!			call DEBUGPrintDefectUpdate(defectUpdate)
			
			!deallocate memory
			reactionChoiceCurrent=>reactionChoiceList%next
			do 51 while(associated(reactionChoiceCurrent))
				reactionTemp=>reactionChoiceCurrent%next
				deallocate(reactionChoiceCurrent%reactants)
				deallocate(reactionChoiceCurrent%products)
				deallocate(reactionChoiceCurrent%cellNumber)
				deallocate(reactionChoiceCurrent%taskid)
				deallocate(reactionChoiceCurrent)
				
				reactionChoiceCurrent=>reactionTemp
			51 continue
			deallocate(reactionChoiceList)
			
		endif
	
	else	!choose a reaction in one volume element
	
		!If implantScheme=='explicit', cascade implantation needs to be carried out explicitly
		if(implantScheme=='explicit') then
			
			if(elapsedTime .GE. numImplantEvents*(numDisplacedAtoms*atomsize)/(totalVolume*DPARate)) then
				
				call addCascadeExplicit(reactionCurrent)
				
				!Do not generate a timestep in this case; this is an explicit (zero-time) reaction
				
				if(myProc%taskid==MASTER) then
					tau=0d0
				endif
				
			else
			
				call chooseReaction(reactionCurrent, CascadeCurrent)
				
				!Generate timestep in the master processor and send it to all other processors
				
				if(myProc%taskid==MASTER) then
					tau=GenerateTimestep()
				endif
				
			endif
			
		else if(implantScheme=='MonteCarlo') then
			
			call chooseReaction(reactionCurrent, CascadeCurrent)
			
			!Generate timestep in the master processor and send it to all other processors
			
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
			endif
		
		else	
		
			write(*,*) 'error choosing reaction main program'
		
		endif
		
		!***********************************************************************************************
		!Update defects according to reactions chosen. Communicate defects that have changed on 
		!boundaries of other processors and create a list of defects whose reaction rates must be update
		!in this processor due to defects updated in this processor and in neighboring processors.
		!***********************************************************************************************
			
		!call DEBUGPrintReaction(reactionCurrent, step)
		!write(*,*) 'tau', tau, 'totalRate', totalRate
	
		!************
		! Optional: count how many steps are null events
		!************
	
		if(.NOT. associated(reactionCurrent)) then
			nullSteps=nullSteps+1
			if(myProc%numtasks==1) then
				write(*,*) 'Error null event in serial SRSCD'
			endif
		endif
	
		call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)

		!call DEBUGPrintDefectUpdate(defectUpdate)
	
		!If a cascade is chosen, update reaction rates for all defects remaining in coarse mesh element
		
		if(associated(reactionCurrent)) then
			if(reactionCurrent%numReactants==-10) then
				
				!Resets reaction list in cell and neighbors (on same processor)
				call resetReactionListSingleCell(reactionCurrent%cellNumber(1))
				
				!Variable tells us whether cascade communication step needs to be carried out
				!and if so in what volume element in the coarse mesh
				cascadeCell=reactionCurrent%cellNumber(1)
				
			endif
		endif
		!Update reaction rates for defects involved in reaction chosen
	
	endif
	
	!Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
	if(myProc%taskid==MASTER) then
		
		elapsedTime=elapsedTime+tau
		
		!NOTE: we should eliminate this send/recieve pair and only track time in the master to save communication
		do 11 i=1,myProc%numtasks-1
			call MPI_SEND(elapsedTime, 1, MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,ierr)
		11 continue

	else
		!slave processors recieve elapsed time from master
		call MPI_RECV(elapsedTime,1,MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,status,ierr)
	endif
	
	call updateReactionList(defectUpdate)

	if(totalRate .LT. 0d0) then
		write(*,*) 'error totalRate less than zero', step
	endif
	
!	call DEBUGPrintDefects(step)
!	call DEBUGPrintReactionList(step)
!	call DEBUGCheckForUnadmissible(reactionCurrent, step)

	!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
	!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
	!release the defects into the coarse mesh.

	if(associated(CascadeCurrent)) then
		if(totalRateCascade(CascadeCurrent) .LT. cascadeReactionLimit) then
			
			!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
			!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
			cascadeCell=CascadeCurrent%cellNumber

			!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
			call releaseFineMeshDefects(CascadeCurrent)

			call resetReactionListSingleCell(cascadeCell)
			
		endif
	endif
	
	!Cascade communication step:
	!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
	!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.

	call cascadeUpdateStep(cascadeCell)

	!Every time a cascade is added, we reset the total rate and check how wrong it has become
	if(associated(reactionCurrent)) then
		if(implantType=='Cascade') then
			if(reactionCurrent%numReactants==-10) then
				!if(myProc%taskid==MASTER) write(*,*) 'Checking Total Rate'
				totalRate=TotalRateCheck()
			endif
		elseif(implantType=='FrenkelPair') then
			if(mod(step,1000)==0) then
				totalRate=TotalRateCheck()
			endif
		endif
	endif
	
	!******************************************
	! Optional: count how many cascades are present at step i and compile to find avg. number
	! of cascades present per step
	!******************************************
	
	!TotalCascades=TotalCascades+CascadeCount()
	
	!******************************************
	! Output
	!******************************************
	
	if(elapsedTime .GT. totalTime/200d0*(2d0)**(outputCounter)) then
		
		call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		
		DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
			(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))
		
		if(myProc%taskid==MASTER) then
			call cpu_time(time2)
			write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
			write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
				'computation time', time2-time1
			write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
			write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
				'computation time', time2-time1
			
			!Optional: output average number of cascades present per step in local processor
			!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
			
			!Optional: output fraction of steps that are null events
			if(singleElemKMC=='yes') then
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
			else
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
			endif
			
			write(84,*)
			write(*,*)
		endif
		
		!Several defect output optionas available, these outputs should be chosen in input file (currently hard coded)
		
		if(rawdatToggle=='yes') call outputDefects()
		!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
		!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
		if(postprToggle=='yes') then
			if(totdatToggle=='yes') then
				call outputDefectsTotal(elapsedTime,step)
			else
				write(*,*) 'Error outputing postpr.out but not totdat.out'
			endif
		endif
		if(xyzToggle=='yes') call outputDefectsXYZ()
		if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
		if(outputDebug=='yes') call outputDebugRestart(outputCounter, elapsedTime)
		
		outputCounter=outputCounter+1

	endif
	
!	if(mod(step,100000)==0) then
		
!		call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
!		call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		
!		DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
!			(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))
		
!		if(myProc%taskid==MASTER) then
!			call cpu_time(time2)
!			write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
!			write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
!				'computation time', time2-time1
!			write(*,*)
!		endif
		
		!***************************************************************************************
		!2014.10.13
		!
		! Adding debugging section, used to find errors in code that occur rarely
		!
		! Every 100000 steps, we will erase the debug file and open a new one in order to 
		! save memory. The debug file contains written outputs at major locations in the code,
		! in order to see where the code is hanging.
		!
		!***************************************************************************************
		
!		if(myProc%taskid==MASTER) then

!			Close and re-open write file (gets rid of old write data)
!			close(86)
!			open(86, file=filename5, action='write', status='Unknown')
			
!!			Write initial information into write file
!			write(86,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
!			write(86,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He Implant Events', numHeImplantTotal, &
!				'computation time', time2-time1
!			write(86,*)
			
!		endif
		
		
!	endif
	
!	if(mod(step,10)==0) then
!		call outputRates(elapsedTime, step)
!	endif
	
10 continue

!***********************************************************************
!Output defects at the end of the implantation loop
!***********************************************************************

call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))

if(myProc%taskid==MASTER) then
	call cpu_time(time2)

!	write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
!	write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
!				'computation time', time2-time1
!	write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
!	write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
!				'computation time', time2-time1
	
!	!Optional: output average number of cascades present per step in local processor
!	!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
	
!	!Optional: output fraction of steps that are null events
!	write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
	
!	write(84,*) 
!	write(*,*)
endif

!call outputDefects()
!call outputDefectsTotal(elapsedTime, step)
!call outputDefectsProfile(sim)
if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
if(outputDebug=='yes') call outputDebugRestart(outputCounter ,elapsedTime)

!***********************************************************************
!Annealing loop: change the temperature to annealTemp and change all
!defect production rates to 0 (dpa rate, He implant rate)
!
!Then carry out simulation normally until elapsedTime=totalTime+annealTime
!***********************************************************************

if(annealTime .GT. 0d0) then
	call annealInitialization()	!Sets up simulation for annealing (change temp, etc)
	annealIdentify=.TRUE.
endif

outputCounter=1		!Used to output once for each annealing step

if(myProc%taskid==MASTER .AND. annealTime .GT. 0d0) then
	write(*,*) 'Entering Annealing Phase'
!	write(82,*) 'Entering Annealing Phase'
	write(83,*) 'Entering Annealing Phase'
	write(84,*) 'Entering Annealing Phase'
endif


!carry out annealing in multiple steps if desired. Each step increases temperature.
do 101 annealIter=1,annealSteps


if(annealType=='mult') then
	
	temperature=annealTemp*annealTempInc**dble(annealIter-1)

else if(annealType=='add') then

	temperature=annealTemp+dble(annealIter-1)*annealTempInc

else
	write(*,*) 'error unknown anneal type'
endif

!Possible bug: this will make the reaction rates for implantation non-zero. Not sure why this wasn't a problem before.
do 111 i=1,numCells
	call resetReactionListSingleCell(i)
111 continue

totalRate=TotalRateCheck()

if(myProc%taskid==MASTER .AND. annealTime .GT. 0d0) then
	write(*,*) 'Anneal step', annealIter, 'temperature', temperature
	write(84,*) 'Anneal step', annealIter, 'temperature', temperature
	write(83,*) 'Anneal step', annealIter, 'temperature', temperature
endif


do 100 while(elapsedTime .LT. totalTime+dble(annealIter)*annealTime/dble(annealSteps))

	step=step+1
	
	!Logical variable tells us whether cascade communication step needs to be carried out
	!(0=no cascade, nonzero=number of volume element where cascade event has happened)
	cascadeCell=0
	
	!If we are choosing one reaction per volume element, find the max reaction rate in the volume elements
	!(we don't need to do this otherwise, totalRate is automatically updated each time a reaction is carried out)
	if(singleElemKMC=='yes') then
		if(meshingType=='adaptive') then
			write(*,*) 'Error adaptive meshing not allowed for single element kMC'
		endif
		rateSingle=0d0
		do 200 cell=1,numCells
			if(totalRateVol(cell) .GT. rateSingle) then
				rateSingle=totalRateVol(cell)
			endif
		200 continue
		totalRate=rateSingle
	endif
	
	call MPI_ALLREDUCE(totalRate,maxRate,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

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
	do 2 i=1,numSpecies
		defectUpdate%defectType(i)=0
	2 continue
	nullify(defectUpdate%next)
	defectUpdateCurrent=>defectUpdate
	
	if(singleElemKMC=='yes') then	!choose a reaction in each volume element
	
		if(implantScheme=='explicit') then
		
			write(*,*) 'Error explicit implantation not implemented for single element kMC'
			
		else
			!Generate timestep in the master processor and send it to all other processors
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
				if(elapsedTime-totalTime+tau .GT. annealTime/dble(annealSteps)*outputCounter) then
					!we have taken a timestep that moves us past this annealing step
					tau=annealTime/dble(annealSteps)*outputCounter-(elapsedTime-totalTime)
				endif
			endif

			!Choose one reaction in each cell
			do 210 cell=1,numCells
			
				!Choose reactions here
				call chooseReactionSingleCell(reactionCurrent, CascadeCurrent, cell)

				!Update defects according to the reaction chosen
					
			!	call DEBUGPrintReaction(reactionCurrent, step)
			
				!************
				! Optional: count how many steps are null events
				!************
			
				if(.NOT. associated(reactionCurrent)) then
					nullSteps=nullSteps+1
				endif
				
				if(associated(reactionCurrent)) then
					if(reactionCurrent%numReactants==-10) then
						write(*,*) 'Error chose casacde implantation in annealing phase'
					endif
				endif
			
				call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
			
			!	call DEBUGPrintDefectUpdate(defectUpdate)
			
			210 continue

		endif
	
	else	!choose a reaction in one volume element
		
		call chooseReaction(reactionCurrent, CascadeCurrent)
		
		!Generate timestep in the master processor and send it to all other processors
		
		if(myProc%taskid==MASTER) then
			tau=GenerateTimestep()
			if(elapsedTime-totalTime+tau .GT. annealTime/dble(annealSteps)*outputCounter) then
				!we have taken a timestep that moves us past this annealing step
				tau=annealTime/dble(annealSteps)*outputCounter-(elapsedTime-totalTime)
			endif
		endif
	
		!***********************************************************************************************
		!Update defects according to reactions chosen. Communicate defects that have changed on 
		!boundaries of other processors and create a list of defects whose reaction rates must be update
		!in this processor due to defects updated in this processor and in neighboring processors.
		!***********************************************************************************************
	
	!	call DEBUGPrintReaction(reactionCurrent, step)
	
		!************
		! Optional: count how many steps are null events
		!************
		
		if(.NOT. associated(reactionCurrent)) then
			nullSteps=nullSteps+1
			!write(*,*) 'null step', myProc%taskid
		endif
	
		call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent)
		
		if(associated(reactionCurrent)) then
			if(reactionCurrent%numReactants==-10) then
				write(*,*) 'Error chose casacde implantation in annealing phase'
			endif
		endif
	
	endif
	
	!Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
	if(myProc%taskid==MASTER) then
		
		elapsedTime=elapsedTime+tau
		
		!NOTE: we should eliminate this send/recieve pair and only track time in the master to save communication
		do 110 i=1,myProc%numtasks-1
			call MPI_SEND(elapsedTime, 1, MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,ierr)
		110 continue

	else
		!slave processors recieve elapsed time from master
		call MPI_RECV(elapsedTime,1,MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,status,ierr)
	endif

!	call DEBUGPrintDefectUpdate(defectUpdate)

	!Update reaction rates for defects involved in reaction chosen
	
	call updateReactionList(defectUpdate)
	
!	call DEBUGPrintDefects(step)
!	call DEBUGPrintReactionList(step)
!	call DEBUGCheckForUnadmissible(reactionCurrent, step)
	
	!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
	!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
	!release the defects into the coarse mesh.
	
	if(associated(CascadeCurrent)) then
		if(totalRateCascade(CascadeCurrent) .LT. cascadeReactionLimit) then
			
			!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
			!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
			cascadeCell=CascadeCurrent%cellNumber

			!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
			call releaseFineMeshDefects(CascadeCurrent)

			call resetReactionListSingleCell(cascadeCell)
			
		endif
	endif
		
	!Cascade communication step:
	!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
	!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.
	call cascadeUpdateStep(cascadeCell)
		
	!Update the totalRate in order to avoid any truncation error every 1000 steps (this value can be modified)
	
	if(mod(step,1000)==0) then
		totalRate=TotalRateCheck()
	endif
	
	!******************************************
	! Optional: count how many cascades are present at step i and compile to find avg. number
	! of cascades present per step
	!******************************************
	
	!TotalCascades=TotalCascades+CascadeCount()
	
	!******************************************
	! Output
	!******************************************
	
	if((elapsedTime-totalTime) .GE. annealTime/dble(annealSteps)*outputCounter) then
		
		call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		
		DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
			(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))
	
		if(myProc%taskid==MASTER) then
			call cpu_time(time2)
			write(*,*) 'time', elapsedTime, 'anneal time', elapsedTime-totalTime, 'dpa', dpa, 'steps', step
			write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
				'computation time', time2-time1
			write(84,*) 'time', elapsedTime, 'anneal time', elapsedTime-totalTime, 'dpa', dpa, 'steps', step
			write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
				'computation time', time2-time1
			
			!Optional: output average number of cascades present per step in local processor
			!write(*,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
			!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
			
			!Optional: output fraction of steps that are null events
			if(singleElemKMC=='yes') then
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
			else
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
			endif
			
			write(84,*)
			write(*,*)
		endif
		
		if(rawdatToggle=='yes') call outputDefects()
		!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
		!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
		if(postprToggle=='yes') then
			if(totdatToggle=='yes') then
				call outputDefectsTotal(elapsedTime,step)
			else
				write(*,*) 'Error outputing postpr.out but not totdat.out'
			endif
		endif
		if(xyzToggle=='yes') call outputDefectsXYZ()
		if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
		if(outputDebug=='yes') call outputDebugRestart(outputCounter, elapsedTime)
		
		outputCounter=outputCounter+1

	endif
	
!Anneal time loop
100 continue

!Multiple anneal steps loop
101 continue

!***********************************************************************
!Output final defect state
!***********************************************************************

call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))

!write(*,*) 'Fraction null steps', dble(nullSteps)/dble(step), 'Proc', myProc%taskid

if(myProc%taskid==MASTER) then
	call cpu_time(time2)
	
	if(sinkEffSearch=='no') then
		write(*,*) 'Final Defect State'
!		write(82,*) 'Final Defect State'
		write(83,*) 'Final Defect State'
		write(84,*) 'Final Defect State'
		
		write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
		write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
					'computation time', time2-time1
		
		write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
		write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'He implant events', numHeImplantTotal, &
					'computation time', time2-time1
		
		!Optional: output average number of cascades present per step in local processor
		!write(*,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
		!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
		
		!Optional: output fraction of steps that are null events
		if(singleElemKMC=='yes') then
			write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
		else
			write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
		endif
		
		!Optional: output sink efficiency of grain boundary
		if(numMaterials == 2) then
			write(*,*) 'Number of emitted defects', numEmitV, numEmitSIA
			write(*,*) 'Number of trapped defects', numTrapV, numTrapSIA
			write(*,*) 'GB V Sink Eff', 1d0-dble(numEmitV)/dble(numTrapV)
			write(*,*) 'GB SIA Sink Eff', 1d0-dble(numEmitSIA)/dble(numTrapSIA)
			
			write(84,*) 'GB_V_Sink_Eff', 1d0-dble(numEmitV)/dble(numTrapV)
			write(84,*) 'GB_SIA_Sink_Eff', 1d0-dble(numEmitSIA)/dble(numTrapSIA)
		endif
		
		write(84,*) 
		write(*,*)
	
	endif
endif

if(rawdatToggle=='yes') call outputDefects()
!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
if(postprToggle=='yes') then
	if(totdatToggle=='yes') then
		call outputDefectsTotal(elapsedTime,step)
	else
		write(*,*) 'Error outputing postpr.out but not totdat.out'
	endif
endif
if(xyzToggle=='yes') call outputDefectsXYZ()
if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
if(outputDebug=='yes') call outputDebugRestart(outputCounter, elapsedTime)

!Defect profile is only output during last output step (not at intermediate steps)
if(profileToggle=='yes') call outputDefectsProfile(sim)

if(sinkEffSearch=='yes') conc_v_store(sim)=computeVConc()
if(sinkEffSearch=='yes') conc_i_store(sim)=computeIConc()

call cpu_time(time2)

!***********************************************************************
!Final step: release all cascades into coarse mesh
!
!(deallocate all extra memory here)
!***********************************************************************

if(myProc%taskid==MASTER) then
	write(*,*) 'computation time', time2-time1
	write(*,*) 'total steps', step
	
	write(84,*) 'computation time', time2-time1
	write(84,*) 'total steps', step

	write(*,*) 'Deallocating memory: fine mesh defects and reactions'
endif

do 12 while(associated(ActiveCascades))
	CascadeCurrent=>ActiveCascades

	!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
	!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
	cascadeCell=CascadeCurrent%cellNumber

	!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
	call releaseFineMeshDefects(CascadeCurrent)

	call resetReactionListSingleCell(cascadeCell)
	
12 continue

if(myProc%taskid==MASTER) then
!	write(82,*) 'Released all fine mesh defects'
	write(83,*) 'Released all fine mesh defects'
	write(84,*) 'Released all fine mesh defects'
endif

!call outputDefects()
!call outputDefectsTotal(elapsedTime, step)

!Final memory cleanup: deallocate defects and reactions in coarse mesh

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating memory: coarse mesh defects and reactions'
endif

call deallocateDefectList()

call deallocateReactionList()

call deallocateBoundarydefectList()

deallocate(totalRateVol)

if(myProc%taskid==MASTER) then
	close(82)
	close(83)
	close(84)
	close(85)
	!close(86)
	!close(87)
endif

!End of loop for multiple trials
13 continue

!***********************************************************************
!7.2.2015 Adding iterative search for effective sink efficiency
!***********************************************************************
if(sinkEffSearch=='yes') then

!	!Calculate average concentrations among all sims
!	conc_i_avg=0d0
!	do 503 sim=1,numSims
!		conc_i_avg=conc_i_avg+conc_i_store(sim)
!	503 continue
!	conc_i_avg=conc_i_avg/dble(numSims)
	
!	!Calculate standard deviation among all sims
!	conc_i_stdev=0d0
!	do 504 sim=1,numSims
!		conc_i_stdev=conc_i_stdev+(conc_i_store(sim)-conc_i_avg)**2d0
!	504 continue
!	conc_i_stdev=dsqrt(conc_i_stdev/dble(numSims))

!	!Calculate average concentrations among all sims
!	conc_v_avg=0d0
!	do 505 sim=1,numSims
!		conc_v_avg=conc_v_avg+conc_v_store(sim)
!	505 continue
!	conc_v_avg=conc_v_avg/dble(numSims)
	
!	!Calculate standard deviation among all sims
!	conc_v_stdev=0d0
!	do 506 sim=1,numSims
!		conc_v_stdev=conc_v_stdev+(conc_v_store(sim)-conc_v_avg)**2d0
!	506 continue
!	conc_v_stdev=dsqrt(conc_v_stdev/dble(numSims))
	
	!Calculate residual - need residual to be +/- in order to find zero crossing
	Residual=0d0
	do 507 sim=1,numSims
		Residual=Residual+(conc_v_store(sim)-conc_v)+(conc_i_store(sim)-conc_i)
	507 continue
	Residual=Residual/dble(numSims)
	
	!Calculate standard deviation of residual
	Residual_stdev=0d0
	do 508 sim=1,numSims
		Residual_stdev=Residual_stdev+((conc_v_store(sim)-conc_v)+&
			(conc_i_store(sim)-conc_i)-Residual)**2d0
	508 continue
	Residual_stdev=dsqrt(Residual_stdev/dble(numSims))

	!Calculate residual of squares: for convergence criterion
	Residual_square=0d0
	do 509 sim=1,numSims
		Residual_square=Residual_square+dsqrt((conc_v_store(sim)-conc_v)**2d0+(conc_i_store(sim)-conc_i)**2d0)
	509 continue
	Residual_square=Residual_square/dble(numSims)
	
	!Calculate standard deviation of residual of squares
	Residual_sqrdev=0d0
	do 510 sim=1,numSims
		Residual_sqrdev=Residual_sqrdev+(dsqrt((conc_v_store(sim)-conc_v)**2d0+&
			(conc_i_store(sim)-conc_i)**2d0)-Residual_square)**2d0
	510 continue
	Residual_sqrdev=dsqrt(Residual_sqrdev/dble(numSims))
	
	!Change the sign of Residual_square to match the sign of Residual
	Residual_square=sign(Residual_square,Residual)

endif

if(alpha_i_search .eqv. .TRUE.) then
	
if(sinkEffSearch=='no') then
	!Exit all sink efficiency search loops
	searchToggle=.FALSE.
	alpha_i_search=.FALSE.
	alpha_v_search=.FALSE.
else if(sinkEffSearch=='yes') then

	i_iteration=i_iteration+1

	if(alpha_i_search .eqv. .TRUE.) then
		!Check if conc_i_average is within stdev of previous step
		if(dabs(Residual_square) .LT. 1.5d0*Residual_sqrdev .AND. dble(Residual) .LT. Residual_stdev) then
			!we have converged for both values of alpha
			alpha_i_search=.FALSE.
			alpha_v_search=.FALSE.
			write(*,*) 'Both converged: alpha_i', alpha_i, 'alpha_v', alpha_v
			write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
			write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
!			read(*,*)

		else if(dabs(Residual_square) .GE. dabs(Residual_prev) &
			.AND. sign(Residual_square,Residual_prev)==Residual_square &
			.AND. i_iteration .GT. 2) then !we have reached a local minimum in alpha_i
			
			alpha_i_search=.FALSE.
			
			write(*,*) 'Alpha_i converged', alpha_i_prev
			alpha_i=alpha_i_prev
			Residual_square=Residual_prev	!Need to do this because we are going to immediately calculate new value of alpha_v, need to store
											!correct value of residual (not the value given by the most recent alpha_i, but the one previous)
			
			!This forces us to keep iterating until we have converged in one step for both v and i
			if(i_iteration .GT. 1 .OR. v_iteration==0) then
				alpha_v_search=.TRUE.
				v_iteration=0
			endif
		else
			!calculate new value for alpha_i
			alpha_temp=alpha_i
			
			if(i_iteration==1) then
				if(alpha_i==0d0) then
					alpha_i=0.1
				else
					alpha_i=alpha_i/2d0
				endif
			else
				alpha_i=(Residual_square*alpha_i_prev-Residual_prev*alpha_i)/(Residual_square-Residual_prev)
			endif
			
			if(alpha_i .LT. 0d0) then
				alpha_i=0d0
			endif
			
			if(alpha_i .GT. 1d0) then
				alpha_i=1d0
			endif
			
			!store previous step
			alpha_i_prev=alpha_temp
			Residual_prev=Residual_square
			
		endif
		
		!Temp: checking if we are doing it right
		write(*,*) 'i_iteration', i_iteration
		write(*,*) 'alpha_i_old', alpha_i_prev, 'alpha_i', alpha_i
		write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
		write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
		write(*,*)
!		write(*,*) 'conc_i', conc_i, 'conc_i_avg', conc_i_avg, 'conc_i_stdev', conc_i_stdev
!		read(*,*)
	endif

else
	write(*,*) 'Error sinkEffSearch unrecognized'
endif

endif !if(alpha_i_search .eqv. .TRUE.) then


if(alpha_v_search .eqv. .TRUE.) then

if(sinkEffSearch=='no') then
	!Exit all sink efficiency search loops
	searchToggle=.FALSE.
	alpha_i_search=.FALSE.
	alpha_v_search=.FALSE.
else if(sinkEffSearch=='yes') then
	
	!i_iteration=0
	v_iteration=v_iteration+1

	if(alpha_v_search .eqv. .TRUE.) then
		!Check if conc_i_average is within stdev of previous step

		if(dabs(Residual_square) .LT. 1.5d0*Residual_sqrdev .AND. dble(Residual) .LT. Residual_stdev) then
			!we have converged for both values of alpha
			alpha_i_search=.FALSE.
			alpha_v_search=.FALSE.
			write(*,*) 'Both converged: alpha_i', alpha_i, 'alpha_v', alpha_v
			write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
			write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
!			read(*,*)

		else if(dabs(Residual_square) .GE. dabs(Residual_prev) &
			.AND. sign(Residual_square,Residual_prev)==Residual_square &
			.AND. v_iteration .GT. 2) then !we have reached a local minimum in alpha_i
		
			alpha_v_search=.FALSE.
			
			write(*,*) 'Alpha_v converged', alpha_v_prev
			alpha_v=alpha_v_prev

			!This forces us to keep iterating until we have converged in one step for both v and i
			if(v_iteration .GT. 1 .OR. i_iteration==0) then
				alpha_i_search=.TRUE.
				i_iteration=0
			endif
			
		else
			!calculate new value for alpha_v
			alpha_temp=alpha_v
			
			if(v_iteration==1) then
				if(alpha_v==0d0) then
					alpha_v=0.1d0
				else
					alpha_v=alpha_v/2d0
				endif
			else
				alpha_v=(Residual_square*alpha_v_prev-Residual_prev*alpha_v)/(Residual_square-Residual_prev)
			endif
			
			if(alpha_v .LT. 0d0) then
				alpha_v=0d0
			endif
			
			if(alpha_v .GT. 1d0) then
				alpha_v=1d0
			endif
			
			!store previous step
			alpha_v_prev=alpha_temp
			Residual_prev=Residual_square

		endif
		
		write(*,*) 'v_iteration', v_iteration
		write(*,*) 'alpha_v_old', alpha_v_prev, 'alpha_v', alpha_v
		write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
		write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
		write(*,*)
		!write(*,*) 'conc_v', conc_v, 'conc_v_avg', conc_v_avg, 'conc_v_stdev', conc_v_stdev
!		read(*,*)
	endif

else
	write(*,*) 'Error sinkEffSearch unrecognized'
endif

endif !if(alpha_v_search .eqv. .TRUE.) then

if(sinkEffSearch=='yes') then
	if(alpha_i_search .eqv. .FALSE.) then
		if(alpha_v_search .eqv. .FALSE.) then
			searchToggle=.FALSE.
			!we have converged for both alpha_v and alpha_i
			write(*,*) 'Converged for both alphas'
		else
			searchToggle=.TRUE.
		endif
	else
		!we have not yet converged for both alpha_v and alpha_i
		searchToggle=.TRUE.
	endif
endif

500 continue

if(sinkEffSearch=='yes') write(*,*) 'Final alpha_v', alpha_v, 'alpha_i', alpha_i

!***********************************************************************
!Last step: deallocating memory
!***********************************************************************
deallocate(cascadeConnectivity)

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating memory: cascade list'
endif

call deallocateCascadeList()

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating material input data'
endif

call deallocateMaterialInput()

write(*,*) 'Finalizing processor', myProc%taskid

close(81)

call MPI_FINALIZE(ierr)		!must be performed at end of simulation to have no error

end program
