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

double precision  elapsedTime, totalTime, tau, GenerateTimestep, TotalRateCheck, rateSingle
integer status(MPI_STATUS_SIZE),  step, annealIter, sim, numDefectsRecv, tracker, outputCounter, nullSteps
integer cascadeCell, i, j, k, cell
integer, allocatable :: cellRecvTot(:), defectRecvTot(:,:)
real time1, time2

integer CascadeCount, TotalCascades !<Used to count the number of cascades present in the simulation

character(12) filename, filename2, filename3, filename4, filename5, filename6

double precision rateDiff	!temporary
type(Reaction), pointer :: reactionTemp

!2019.05.13 Add: used to change double to integer
!double precision CuAtomsEverMeshTemp
!character*20 CuAtomsEverMeshTemp2

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
	
	subroutine updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent, step)
	use DerivedType
	type(reaction), pointer :: reactionCurrent
	type(DefectUpdateTracker), pointer :: defectUpdateCurrent
	type(cascade), pointer :: CascadeCurrent
	integer step
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

!Initialize input parameters
call initializeMesh()			!open Mesh_xx.txt file and carry out parallel mesh initialization routine
call selectMaterialInputs()		!open material input files (Fe_Defects.txt) and read all relevant data (migration and binding energies, etc)
call readCascadeList()			!input cascade list (Fe_Cascades.txt) from file to cascadeList (global variable)
call readImplantData()			!input (1-dimensional) non-uniform defect implantation info (DPA rates, He implant rates)
call readParameters()			!read simulation parameters (DPA rate, temperature, etc)

!Create fine mesh connectivity
allocate(cascadeConnectivity(numCellsCascade, 6))	!< numCellsCascade=numxFine*numyFine*numzFine
call createCascadeConnectivity()

!***********************************************************************
!7.2.2015 Creating iterative search for sink efficiency, looping here
!***********************************************************************
searchToggle=.TRUE.

!Begin by searching for alpha_i, then switch to alpha_v
alpha_i_search=.TRUE.
alpha_v_search=.FALSE.

allocate(conc_i_store(numSims))
allocate(conc_v_store(numSims))

!Initialize values of alpha_v and alpha_i
if(sinkEffSearch=='yes') then	!< sinkEffSearch was set to 'no' in Read_input.f90
	conc_i_prev=0d0
	conc_v_prev=0d0
	
	alpha_i_prev=0d0
	alpha_v_prev=0d0
	
	v_iteration=0
	i_iteration=0
end if

do while(searchToggle .eqv. .TRUE.)
!do 501 while(alpha_v_search .eqv. .TRUE.)
!do 502 while(alpha_i_search .eqv. .TRUE.)

!***********************************************************************
!For running multiple simulations, initialize the defects, reactions,
!boundary, etc. and loop here.
!***********************************************************************
!< begin the outer loop
do sim=1,numSims

!Reset the temperature to tempStore at the beginning of each simulation
temperature=tempStore	!< Temperature read in (K)

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
end if

!filename5(1:6)='debugg'
!write(unit=filename5(7:8), fmt='(I2)') myProc%taskid
!filename5(9:12)='.log'
!open(86, file=filename5, action='write', status='Unknown')

if(myProc%taskid==MASTER) then
	write(*,*) 'Initializing trial', sim, 'proc', myProc%taskid, 'temperature', temperature
end if

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
allocate(defectList(numCells))		!Create list of defects - array
allocate(reactionList(numCells))	!Create list of reactions - array
allocate(totalRateVol(numCells))	!Create array of total rates in each volume element
call initializeVIdefect()			!2019.05.30 Add
call initializeDefectList()			!initialize defects within myMesh
call initializeBoundaryDefectList()	!initialize defects on boundary of myMesh (in other procs)
call computeVconcent()
call initializeReactionList()		!initialize reactions within myMesh
call initializeTotalRate()			!initialize totalRate and maxRate using reactionList(:)
call initializeDebugRestart()		!input defects into coarse mesh from restart file (for debugging)

if(myProc%taskid==MASTER) then
	write(*,*) 'atomsEverMesh', atomsEverMesh, 'CuEverMesh', CuAtomsEverMesh
	write(*,*) 'initialTotalV', initialTotalV, 'initialTotalSIA', initialTotalSIA
	write(*,*) 'initialCeqv', initialCeqv
	write(*,*) 'initialCeqi', initialCeqi
end if

call DEBUGPrintReactionList(0)		!prints all reaction lists at a given Monte Carlo step
!call DEBUGPrintDefectList(0)
!******************************************************************
!Initialize Counters
!******************************************************************

if(debugToggle=='yes') then		!inpput parameter
	
	!If we are restarting from a reset file, then we have to start wtih a
	!nonzero dpa and at a nonzero time. All of the 'old' implant events
	!and He implant events are tracked in the master processor.
	
	if(myProc%taskid==MASTER) then
		numImplantEvents = numImplantEventsReset		!numImplantEventsReset is read in from debug restart file
	else
		numImplantEvents = 0
	end if
	elapsedTime	= elapsedTimeReset
	!Reset the reactions within this cell and diffusion to neighboring
	!cells in the same processor
	do i=1,numCells
		call resetReactionListSingleCell(i)
	end do

	!Create the boundary defects in this processor and send boundary
	!defects to neighboring processor. Update reactions for diffusion
	!into boundary (into other processors)
	
	!NOTE: this will crash if the number of cells in a neighboring processor
	!is not the same as in the local processor. Easy fix, but putting
	!off until later.
	
	do i=1,numCells
		call cascadeUpdateStep(i)
	end do
	
	write(*,*) 'Processor', myProc%taskid, 'Defect lists and reaction lists initialized'
	
else
	numImplantEvents	= 0		!<Postprocessing: number of Frenkel pairs / cascades (local)
	numAnnihilate		= 0		!<Postprocessing: number of annihilation reactions carried out
	elapsedTime			= 0d0	!The simulation time that has passed
	
	numTrapSIA=0	!<Postprocessing: number of SIAs trapped on grain boundary
	numTrapV=0		!<Postprocessing: number of vacancies trapped on grain boundary
	numEmitSIA=0	!<Postprocessing: number of SIAs emitted from grain boundary
	numEmitV=0		!<Postprocessing: number of vacancies emitted from grain boundary
end if

step=0
nullSteps=0		!Record the number of steps in which an empty event was selected

!2019.04.30 Add
if(totalDPA > 0d0 .AND. DPARate > 0d0) then
	totalTime=totalDPA/DPARate	!simulation time
else
	totalTime = agingTime
end if

TotalCascades=0
outputCounter=0
nullify(ActiveCascades)		! pointer ActiveCascades=NULL
annealIdentify=.FALSE.		!(.TRUE. if in annealing phase, .FALSE. otherwise) used to determine how to reset reaction rates (should we include implantation or not)

do while(elapsedTime < totalTime)
	
	step=step+1
	call computeVconcent()

	!Logical variable tells us whether cascade communication step needs to be carried out
	!(0=no cascade, nonzero=number of volume element where cascade event has happened)
	cascadeCell=0
	
	!If we are choosing one reaction per volume element, find the max reaction rate in the volume elements
	!(we don't need to do this otherwise, totalRate is automatically updated each time a reaction is carried out)
	if(singleElemKMC=='yes') then
		if(meshingType=='adaptive') then
			write(*,*) 'Error adaptive meshing not allowed for single element kMC'
		end if
		rateSingle=0d0
		do cell=1,numCells
			if(totalRateVol(cell) > rateSingle) then
				rateSingle=totalRateVol(cell)
			end if
		end do
		totalRate=rateSingle	!find the maximum totalRate in the local peocessor
	end if
	
	call MPI_ALLREDUCE(totalRate,maxRate,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

!*************************************************************************
!	if(mod(step,10)==0) then
!		!Debugging subroutine: outputs the reaction rate in each processor
!		call outputRates(elapsedTime, step)
!	endif
!*************************************************************************
	
	!***********************************************************************************************
	!Choose from reactions in reactionList(:) (local to this processor). Null events possible.
	!
	!If explicit implant scheme has been chosen, implant cascades when elapsed time has passed 
	!a threshold, and do not iterate timestep when doing so.
	!***********************************************************************************************
	
	allocate(defectUpdate)	!Allocate memory for defectUpdate
	allocate(defectUpdate%defectType(numSpecies))
	do i=1,numSpecies
		defectUpdate%defectType(i)=0	!Initialize
	end do
	nullify(defectUpdate%next)	!defectUpdate%next=NULL
	defectUpdateCurrent=>defectUpdate	!set defectUpdateCurrent point to the defect to be updated
	
	if(singleElemKMC=='yes') then	!choose a reaction in each volume element
	
		if(implantScheme=='explicit') then
			write(*,*) 'Error explicit implantation not implemented for single element kMC'
		else
			!Generate timestep in the master processor and send it to all other processors
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
			end if
			
			allocate(reactionChoiceList)	!type(reaction). Allocate memory for reactionChoiceList
			nullify(reactionChoiceList%next)	!set pointer reactionChoiceList%next=NULL
			reactionChoiceCurrent=>reactionChoiceList	!the first memory of reactionChoiceList hasn't data
					
			!Choose one reaction in each cell
			do cell=1,numCells
				!Choose reactions here
				!Inout: cell.
				!Output: reactionCurrent, CascadeCurrent, numImplantEvents, numHeImplantEvents, numAnnihilate
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
					
					allocate(reactionChoiceCurrent%next)	!Allocate memory for the next
					reactionChoiceCurrent=>reactionChoiceCurrent%next	!point to next
					
					reactionChoiceCurrent%numReactants=reactionCurrent%numReactants
					allocate(reactionChoiceCurrent%reactants(reactionCurrent%numReactants, numSpecies))
					
					reactionChoiceCurrent%numProducts=reactionCurrent%numProducts
					allocate(reactionChoiceCurrent%products(reactionCurrent%numProducts, numSpecies))
					
					allocate(reactionChoiceCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
					allocate(reactionChoiceCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
					
					reactionChoiceCurrent%reactionRate=reactionCurrent%reactionRate
					
					do i=1,reactionCurrent%numReactants
					
						do j=1,numSpecies
							reactionChoiceCurrent%reactants(i,j)=reactionCurrent%reactants(i,j)
						end do
					
						reactionChoiceCurrent%cellNumber(i)=reactionCurrent%cellNumber(i)
						reactionChoiceCurrent%taskid(i)=reactionCurrent%taskid(i)
					
					end do
					
					do i=1,reactionCurrent%numProducts
						do j=1,numSpecies
							reactionChoiceCurrent%products(i,j)=reactionCurrent%products(i,j)
						end do
						
						reactionChoiceCurrent%cellNumber(i+reactionCurrent%numReactants) = &
							reactionCurrent%cellNumber(i+reactionCurrent%numReactants)
							
						reactionChoiceCurrent%taskid(i+reactionCurrent%numReactants) = &
							reactionCurrent%taskid(i+reactionCurrent%numReactants)
							
					end do
					
					nullify(reactionChoiceCurrent%next)		!set pointer reactionChoiceCurrent%next=NULL
					
				end if
			
				!call DEBUGPrintDefectUpdate(defectUpdate)
			
				!No cascade choices allowed for one KMC domain per element at this point
			
			end do
			
			reactionChoiceCurrent=>reactionChoiceList%next
			call updateDefectListMultiple(reactionChoiceCurrent, defectUpdateCurrent, CascadeCurrent)
			
			!call DEBUGPrintDefectUpdate(defectUpdate)
			
			!deallocate memory
			reactionChoiceCurrent=>reactionChoiceList%next
			do while(associated(reactionChoiceCurrent))
				reactionTemp=>reactionChoiceCurrent%next
				deallocate(reactionChoiceCurrent%reactants)
				deallocate(reactionChoiceCurrent%products)
				deallocate(reactionChoiceCurrent%cellNumber)
				deallocate(reactionChoiceCurrent%taskid)
				deallocate(reactionChoiceCurrent)
				
				reactionChoiceCurrent=>reactionTemp
			end do
			deallocate(reactionChoiceList)
			
		end if
	
	else	!choose a reaction in one volume element
	
		!If implantScheme=='explicit', cascade implantation needs to be carried out explicitly
		if(implantScheme=='explicit') then
			
			if(elapsedTime >= numImplantEvents*(numDisplacedAtoms*atomsize)/(totalVolume*DPARate)) then
				
				!Input: none
				!Output: reactionCurrent, pointing at cascade reaction
				call addCascadeExplicit(reactionCurrent)
				
				!Do not generate a timestep in this case; this is an explicit (zero-time) reaction
				
				if(myProc%taskid==MASTER) then
					tau=0d0
				end if
				
			else

				!Input:  none
				!Output: reactionCurrent, CascadeCurrent
				call chooseReaction(reactionCurrent, CascadeCurrent)
				
				!Generate timestep in the master processor and send it to all other processors
				
				if(myProc%taskid==MASTER) then
					tau=GenerateTimestep()
				end if
				
			end if
			
		else if(implantScheme=='MonteCarlo') then
			
			call chooseReaction(reactionCurrent, CascadeCurrent)
			!test
			!call DEBUGPrintReaction(reactionCurrent, step)
			!Generate timestep in the master processor and send it to all other processors
			
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
			end if
		
		else	
		
			write(*,*) 'error choosing reaction main program'
		
		end if
		
		!***********************************************************************************************
		!Update defects according to reactions chosen. Communicate defects that have changed on 
		!boundaries of other processors and create a list of defects whose reaction rates must be update
		!in this processor due to defects updated in this processor and in neighboring processors.
		!***********************************************************************************************
			
		!call DEBUGPrintReaction(reactionCurrent, step)
		!write(*,*) 'tau', tau, 'totalRate', totalRate

		if(.NOT. associated(reactionCurrent)) then
			nullSteps=nullSteps+1
			!if(myProc%numtasks==1) then
			!	write(*,*) 'Error null event in serial SRSCD'
			!endif
		end if

		!Input: reactionCurrent
		!Output: defectUpdateCurrent, CascadeCurrent
		call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent, step)

		!call DEBUGPrintDefectList(step)

		!call DEBUGPrintDefectUpdate(defectUpdate)
	
		!If a cascade is chosen, update reaction rates for all defects remaining in coarse mesh element
		
		if(associated(reactionCurrent)) then
			if(reactionCurrent%numReactants==-10) then
				
				!Resets reaction list in cell and neighbors (on same processor)
				call resetReactionListSingleCell(reactionCurrent%cellNumber(1))
				
				!Variable tells us whether cascade communication step needs to be carried out
				!and if so in what volume element in the coarse mesh
				cascadeCell=reactionCurrent%cellNumber(1)
				
			end if
		end if
		!Update reaction rates for defects involved in reaction chosen
	
	end if
	
	!Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
	if(myProc%taskid==MASTER) then
		
		elapsedTime=elapsedTime+tau
	end if

	call MPI_BCAST(elapsedTime, 1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,ierr)

!*********************************************************
	call updateReactionList(defectUpdate)

	!call DEBUGPrintReactionList(step)

	if(totalRate < 0d0) then
		write(*,*) 'error totalRate less than zero', step
	end if
	
!	call DEBUGPrintDefects(step)
!	call DEBUGPrintReactionList(step)
!	call DEBUGCheckForUnadmissible(reactionCurrent, step)

	!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
	!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
	!release the defects into the coarse mesh.

	if(associated(CascadeCurrent)) then
		if(totalRateCascade(CascadeCurrent) < cascadeReactionLimit) then
			
			!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
			!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
			cascadeCell=CascadeCurrent%cellNumber

			!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
			call releaseFineMeshDefects(CascadeCurrent)

			call resetReactionListSingleCell(cascadeCell)
			
		end if
	end if
	
	!Cascade communication step:
	!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
	!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.

	call cascadeUpdateStep(cascadeCell)

	!****************************************************************************************
	!Every time a cascade is added, we reset the total rate and check how wrong it has become
	!!'Cascade' implant type, check total rate when implantation event is choosen
	!!'FrenkelPair' implant type, check total rate every 1000 steps
	!*****************************************************************************************

	if(associated(reactionCurrent)) then
		if(implantType=='Cascade') then
			if(reactionCurrent%numReactants==-10) then
				!if(myProc%taskid==MASTER) write(*,*) 'Checking Total Rate'
				totalRate=TotalRateCheck()
			end if
		else if(implantType=='FrenkelPair') then
			if(mod(step,1000)==0) then
				totalRate=TotalRateCheck()
			end if
		end if
	end if
	
	!******************************************
	! Optional: count how many cascades are present at step i and compile to find avg. number
	! of cascades present per step
	!******************************************
	
	TotalCascades=TotalCascades+CascadeCount()
	
	!********************************************************************************
	! Output according to outputCounter
	!********************************************************************************
	
	if(elapsedTime >= 1.0d-4*(1.0d1)**(outputCounter)) then
	! or if(mod(step,100000)==0) then
		call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		!call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		
		DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
			(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))
		
		if(myProc%taskid==MASTER) then
			call cpu_time(time2)
			write(*,*)
			write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
			write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
			if(implantType=='FrenkelPair') then
				write(*,*) 'Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
				write(84,*) 'Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
			else if(implantType=='Cascade')	then
				write(*,*) 'Cascades', totalImplantEvents, 'computation time', time2-time1
				write(84,*) 'Cascades', totalImplantEvents, 'computation time', time2-time1
			else
				write(*,*) 'No implantation', totalImplantEvents, 'computation time', time2-time1
				write(84,*) 'No implantation', totalImplantEvents, 'computation time', time2-time1
			end if

			!Optional: output average number of cascades present per step in local processor
			!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
			
			!Optional: output fraction of steps that are null events
			if(singleElemKMC=='yes') then
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
			else
				write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
			end if
			
			write(84,*)
		end if
		
		!Several defect output optionas available, these outputs should be chosen in input file (currently hard coded)
		
		if(rawdatToggle=='yes') call outputDefects(elapsedTime,step)
		!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
		!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
		if(postprToggle=='yes') then
			if(totdatToggle=='yes') then
				call outputDefectsTotal(elapsedTime,step)
			else
				write(*,*) 'Error outputing postpr.out but not totdat.out'
			end if
		end if
		if(xyzToggle=='yes') call outputDefectsXYZ(elapsedTime,step)	!write(87,*): defect.xyz
		if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)	!write(88,*): VTKout.vtk
		if(outputDebug=='yes') call outputDebugRestart(outputCounter, elapsedTime)	!write(88,*): Restart.in
		
		outputCounter=outputCounter+1

		!***************************************************************************************
		! Adding debugging section, used to find errors in code that occur rarely
		!
		! Every 100000 steps, we will erase the debug file and open a new one in order to
		! save memory. The debug file contains written outputs at major locations in the code,
		! in order to see where the code is hanging.
		!***************************************************************************************

		!if(myProc%taskid==MASTER) then

		!	!Close and re-open write file (gets rid of old write data)
		!	close(86)
		!	open(86, file=filename5, action='write', status='Unknown')

		!	!Write initial information into write file
		!	write(86,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
		!	write(86,*) 'Cascades/Frenkel pairs', totalImplantEvents,  'computation time', time2-time1
		!	write(86,*)
		!endif

		!if(mod(step,10)==0) then
		!	call outputRates(elapsedTime, step)
		!endif

	end if

end do

!***********************************************************************
!Output defects at the end of the implantation loop
!***********************************************************************

call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
!call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))

if(myProc%taskid==MASTER) then
	call cpu_time(time2)

	write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
	write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
	write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
	write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
	
	!Optional: output average number of cascades present per step in local processor
	write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
	
	!Optional: output fraction of steps that are null events
	write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
	
	write(84,*)
	write(*,*)
end if

!Final output
call outputDefects(elapsedTime,step)	!write(82,*): rawdat.out
call outputDefectsTotal(elapsedTime, step)	!write(83,*): totdat.out, write(84,*): postpr.out
!call outputDefectsProfile(sim)	!write(99,*): DepthProfile.out
if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
if(outputDebug=='yes') call outputDebugRestart(outputCounter ,elapsedTime)



!***********************************************************************
!Annealing loop: change the temperature to annealTemp and change all
!defect production rates to 0 (dpa rate, He implant rate)
!Then carry out simulation normally until elapsedTime=totalTime+annealTime
!***********************************************************************

if(annealTime > 0d0) then
	call annealInitialization()	!Sets up simulation for annealing (change temp, etc)
	annealIdentify=.TRUE.
end if

outputCounter=1		!Used to output once for each annealing step

if(myProc%taskid==MASTER .AND. annealTime > 0d0) then
	write(*,*) 'Entering Annealing Phase'
	write(82,*) 'Entering Annealing Phase'
	write(83,*) 'Entering Annealing Phase'
	write(84,*) 'Entering Annealing Phase'
end if

!begin annealing
!carry out annealing in multiple steps if desired. Each step increases temperature.
do annealIter=1,annealSteps	!default value: annealSteps = 1


	if(annealType=='mult') then
	
		temperature=annealTemp*annealTempInc**dble(annealIter-1)

	else if(annealType=='add') then

		temperature=annealTemp+dble(annealIter-1)*annealTempInc

	else
		write(*,*) 'error unknown anneal type'
	end if

	!Possible bug: this will make the reaction rates for implantation non-zero. Not sure why this wasn't a problem before.
	do i=1,numCells
		call resetReactionListSingleCell(i)
	end do

	totalRate=TotalRateCheck()

	if(myProc%taskid==MASTER .AND. annealTime > 0d0) then
		write(*,*) 'Anneal step', annealIter, 'temperature', temperature
		write(84,*) 'Anneal step', annealIter, 'temperature', temperature
		write(83,*) 'Anneal step', annealIter, 'temperature', temperature
	end if

	!if annealTime = 0, the following "do" does not execute
	do while(elapsedTime < totalTime+dble(annealIter)*annealTime/dble(annealSteps))

		step=step+1
	
		!Logical variable tells us whether cascade communication step needs to be carried out
		!(0=no cascade, nonzero=number of volume element where cascade event has happened)
		cascadeCell=0
	
		!If we are choosing one reaction per volume element, find the max reaction rate in the volume elements
		!(we don't need to do this otherwise, totalRate is automatically updated each time a reaction is carried out)
		if(singleElemKMC=='yes') then
			if(meshingType=='adaptive') then
				write(*,*) 'Error adaptive meshing not allowed for single element kMC'
			end if
			rateSingle=0d0
			do cell=1,numCells
				if(totalRateVol(cell) > rateSingle) then
					rateSingle=totalRateVol(cell)
				end if
			end do
			totalRate=rateSingle
		end if
	
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
		do i=1,numSpecies
			defectUpdate%defectType(i)=0
		end do
		nullify(defectUpdate%next)
		defectUpdateCurrent=>defectUpdate
	
		if(singleElemKMC=='yes') then	!choose a reaction in each volume element
	
			if(implantScheme=='explicit') then
		
				write(*,*) 'Error explicit implantation not implemented for single element kMC'
			
			else
				!Generate timestep in the master processor and send it to all other processors
				if(myProc%taskid==MASTER) then
					tau=GenerateTimestep()
					if(elapsedTime-totalTime+tau > annealTime/dble(annealSteps)*outputCounter) then
						!we have taken a timestep that moves us past this annealing step
						tau=annealTime/dble(annealSteps)*outputCounter-(elapsedTime-totalTime)
					end if
				end if

				!Choose one reaction in each cell
				do cell=1,numCells
			
					!Choose reactions here
					call chooseReactionSingleCell(reactionCurrent, CascadeCurrent, cell)

					!Update defects according to the reaction chosen
					
					!call DEBUGPrintReaction(reactionCurrent, step)
			
					!************
					! Optional: count how many steps are null events
					!************
			
					if(.NOT. associated(reactionCurrent)) then
						nullSteps=nullSteps+1
					end if
				
					if(associated(reactionCurrent)) then
						if(reactionCurrent%numReactants==-10) then
							write(*,*) 'Error chose casacde implantation in annealing phase'
						end if
					end if
			
					call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent,step)
			
					!call DEBUGPrintDefectUpdate(defectUpdate)
			
				end do

			end if
	
		else	!choose a reaction in one volume element
		
			call chooseReaction(reactionCurrent, CascadeCurrent)
		
			!Generate timestep in the master processor and send it to all other processors
		
			if(myProc%taskid==MASTER) then
				tau=GenerateTimestep()
				if(elapsedTime-totalTime+tau > annealTime/dble(annealSteps)*outputCounter) then
					!we have taken a timestep that moves us past this annealing step
					tau=annealTime/dble(annealSteps)*outputCounter-(elapsedTime-totalTime)
				end if
			end if
	
			!***********************************************************************************************
			!Update defects according to reactions chosen. Communicate defects that have changed on
			!boundaries of other processors and create a list of defects whose reaction rates must be update
			!in this processor due to defects updated in this processor and in neighboring processors.
			!***********************************************************************************************
	
			!call DEBUGPrintReaction(reactionCurrent, step)
	
			!************
			! Optional: count how many steps are null events
			!************
		
			if(.NOT. associated(reactionCurrent)) then
				nullSteps=nullSteps+1
				!write(*,*) 'null step', myProc%taskid
			end if
	
			call updateDefectList(reactionCurrent, defectUpdateCurrent, CascadeCurrent, step)
		
			if(associated(reactionCurrent)) then
				if(reactionCurrent%numReactants==-10) then
					write(*,*) 'Error chose casacde implantation in annealing phase'
				end if
			end if
	
		end if
	
		!Update elapsed time based on tau, generated timestep. If cascade implant chosen in explicit scheme, tau=0d0
		if(myProc%taskid==MASTER) then
		
			elapsedTime=elapsedTime+tau
		end if

		call MPI_BCAST(elapsedTime, 1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,ierr)
		!***********************************************************************************
		!call DEBUGPrintDefectUpdate(defectUpdate)

		!Update reaction rates for defects involved in reaction chosen
	
		call updateReactionList(defectUpdate)
	
!		call DEBUGPrintDefects(step)
!		call DEBUGPrintReactionList(step)
!		call DEBUGCheckForUnadmissible(reactionCurrent, step)
	
		!If we have chosen an event inside a fine mesh, we check the total reaction rate within that
		!fine mesh. If the total rate is less than a set value, we assume the cascade is annealed and
		!release the defects into the coarse mesh.
	
		if(associated(CascadeCurrent)) then
			if(totalRateCascade(CascadeCurrent) < cascadeReactionLimit) then
			
				!Record the coarse mesh cell number of cascade (it is destroyed in releaseFineMeshDefects)
				!Used to reset reaction list and to tell cascadeUpdateStep whether a cascade event has occurred
				cascadeCell=CascadeCurrent%cellNumber

				!Release cascade defects into coarse mesh cell and reset the reaction list within that cell
				call releaseFineMeshDefects(CascadeCurrent)

				call resetReactionListSingleCell(cascadeCell)
			
			end if
		end if
		
		!Cascade communication step:
		!Tell neighbors whether a cascade has occurred in a cell that is a boundary of a neighbor.
		!If so, update boundary mesh (send defects to boundary mesh) and update all diffusion reaction rates.
		call cascadeUpdateStep(cascadeCell)
		
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
	
		if((elapsedTime-totalTime) >= annealTime/dble(annealSteps)*outputCounter) then
		
			call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
			!call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		
			DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
				(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))
	
			if(myProc%taskid==MASTER) then
				call cpu_time(time2)
				write(*,*) 'time', elapsedTime, 'anneal time', elapsedTime-totalTime, 'dpa', dpa, 'steps', step
				write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
				write(84,*) 'time', elapsedTime, 'anneal time', elapsedTime-totalTime, 'dpa', dpa, 'steps', step
				write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
			
				!Optional: output average number of cascades present per step in local processor
				!write(*,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
				!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
			
				!Optional: output fraction of steps that are null events
				if(singleElemKMC=='yes') then
					write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
				else
					write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
				end if
			
				write(84,*)
				write(*,*)
			end if
		
			if(rawdatToggle=='yes') call outputDefects(elapsedTime,step)
			!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
			!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
			if(postprToggle=='yes') then
				if(totdatToggle=='yes') then
					call outputDefectsTotal(elapsedTime,step)
				else
					write(*,*) 'Error outputing postpr.out but not totdat.out'
				end if
			end if
			if(xyzToggle=='yes') call outputDefectsXYZ(elapsedTime,step)
			if(vtkToggle=='yes') call outputDefectsVTK(outputCounter)
			if(outputDebug=='yes') call outputDebugRestart(outputCounter, elapsedTime)
		
			outputCounter=outputCounter+1

		end if
	
	!Anneal time loop
	end do

!Multiple anneal steps loop
end do

!***********************************************************************
!Output final defect state
!***********************************************************************

call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
!call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

DPA=dble(totalImplantEvents)/(((myProc%globalCoord(2)-myProc%globalCoord(1))*(myProc%globalCoord(4)-myProc%globalCoord(3))*&
	(myProc%globalCoord(6)-myProc%globalCoord(5)))/(numDisplacedAtoms*atomsize))

write(*,*) 'Fraction null steps', dble(nullSteps)/dble(step), 'Proc', myProc%taskid

if(myProc%taskid==MASTER) then
	call cpu_time(time2)
	
	if(sinkEffSearch=='no') then
		write(*,*) 'Final Defect State'
		write(82,*) 'Final Defect State'
		write(83,*) 'Final Defect State'
		write(84,*) 'Final Defect State'
		
		write(*,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
		write(*,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
		
		write(84,*) 'time', elapsedTime, 'dpa', dpa, 'steps', step
		write(84,*) 'Cascades/Frenkel pairs', totalImplantEvents, 'computation time', time2-time1
		
		!Optional: output average number of cascades present per step in local processor
		!write(*,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
		!write(84,*) 'Processor ', myProc%taskid, 'Avg. cascades present', dble(TotalCascades)/dble(step)
		
		!Optional: output fraction of steps that are null events
		if(singleElemKMC=='yes') then
			write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step*numCells)
		else
			write(84,*) 'Fraction null steps', dble(nullSteps)/dble(step)
		end if
		
		!Optional: output sink efficiency of grain boundary
		if(numMaterials == 2) then
			write(*,*) 'Number of emitted defects', numEmitV, numEmitSIA
			write(*,*) 'Number of trapped defects', numTrapV, numTrapSIA
			write(*,*) 'GB V Sink Eff', 1d0-dble(numEmitV)/dble(numTrapV)
			write(*,*) 'GB SIA Sink Eff', 1d0-dble(numEmitSIA)/dble(numTrapSIA)
			
			write(84,*) 'GB_V_Sink_Eff', 1d0-dble(numEmitV)/dble(numTrapV)
			write(84,*) 'GB_SIA_Sink_Eff', 1d0-dble(numEmitSIA)/dble(numTrapSIA)
		end if
		
		write(84,*) 
		write(*,*)
	
	end if
end if

if(rawdatToggle=='yes') call outputDefects(elapsedTime,step)
!if(sinkEffSearch=='no' .AND. numMaterials .GT. 1) call outputDefectsBoundary(elapsedTime,step)
!if(sinkEffSearch=='no' .AND. numMaterials==1) call outputDefectsTotal(elapsedTime, step)
if(postprToggle=='yes') then
	if(totdatToggle=='yes') then
		call outputDefectsTotal(elapsedTime,step)
	else
		write(*,*) 'Error outputing postpr.out but not totdat.out'
	end if
end if
if(xyzToggle=='yes') call outputDefectsXYZ(elapsedTime,step)
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
	write(82,*) 'Released all fine mesh defects'
	write(83,*) 'Released all fine mesh defects'
	write(84,*) 'Released all fine mesh defects'
end if

!call outputDefects()
!call outputDefectsTotal(elapsedTime, step)

!Final memory cleanup: deallocate defects and reactions in coarse mesh

if(myProc%taskid==MASTER) then
	write(*,*) 'Deallocating memory: coarse mesh defects and reactions'
end if

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
	close(87)
end if

!End of loop for multiple trials
end do

!***********************************************************************
!7.2.2015 Adding iterative search for effective sink efficiency
!***********************************************************************
if(sinkEffSearch=='yes') then	!default sinkEffSearch = 'no'

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
	do sim=1,numSims
		Residual=Residual+(conc_v_store(sim)-conc_v)+(conc_i_store(sim)-conc_i)
	end do
	Residual=Residual/dble(numSims)
	
	!Calculate standard deviation of residual
	Residual_stdev=0d0
	do sim=1,numSims
		Residual_stdev=Residual_stdev+((conc_v_store(sim)-conc_v)+&
			(conc_i_store(sim)-conc_i)-Residual)**2d0
	end do
	Residual_stdev=dsqrt(Residual_stdev/dble(numSims))

	!Calculate residual of squares: for convergence criterion
	Residual_square=0d0
	do sim=1,numSims
		Residual_square=Residual_square+dsqrt((conc_v_store(sim)-conc_v)**2d0+(conc_i_store(sim)-conc_i)**2d0)
	end do
	Residual_square=Residual_square/dble(numSims)
	
	!Calculate standard deviation of residual of squares
	Residual_sqrdev=0d0
	do sim=1,numSims
		Residual_sqrdev=Residual_sqrdev+(dsqrt((conc_v_store(sim)-conc_v)**2d0+&
			(conc_i_store(sim)-conc_i)**2d0)-Residual_square)**2d0
	end do
	Residual_sqrdev=dsqrt(Residual_sqrdev/dble(numSims))
	
	!Change the sign of Residual_square to match the sign of Residual
	Residual_square=sign(Residual_square,Residual)

end if

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
			if(dabs(Residual_square) < 1.5d0*Residual_sqrdev .AND. dble(Residual) < Residual_stdev) then
				!we have converged for both values of alpha
				alpha_i_search=.FALSE.
				alpha_v_search=.FALSE.
				write(*,*) 'Both converged: alpha_i', alpha_i, 'alpha_v', alpha_v
				write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
				write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
!				read(*,*)

			else if(dabs(Residual_square) >= dabs(Residual_prev) &
				.AND. sign(Residual_square,Residual_prev)==Residual_square &
				.AND. i_iteration > 2) then !we have reached a local minimum in alpha_i
			
				alpha_i_search=.FALSE.
			
				write(*,*) 'Alpha_i converged', alpha_i_prev
				alpha_i=alpha_i_prev
				Residual_square=Residual_prev	!Need to do this because we are going to immediately calculate new value of alpha_v, need to store
											!correct value of residual (not the value given by the most recent alpha_i, but the one previous)
			
				!This forces us to keep iterating until we have converged in one step for both v and i
				if(i_iteration > 1 .OR. v_iteration==0) then
					alpha_v_search=.TRUE.
					v_iteration=0
				end if
			else
				!calculate new value for alpha_i
				alpha_temp=alpha_i
			
				if(i_iteration==1) then
					if(alpha_i==0d0) then
						alpha_i=0.1
					else
						alpha_i=alpha_i/2d0
					end if
				else
					alpha_i=(Residual_square*alpha_i_prev-Residual_prev*alpha_i)/(Residual_square-Residual_prev)
				end if
			
				if(alpha_i < 0d0) then
					alpha_i=0d0
				end if
			
				if(alpha_i > 1d0) then
					alpha_i=1d0
				end if
			
				!store previous step
				alpha_i_prev=alpha_temp
				Residual_prev=Residual_square
			
			end if
		
			!Temp: checking if we are doing it right
			write(*,*) 'i_iteration', i_iteration
			write(*,*) 'alpha_i_old', alpha_i_prev, 'alpha_i', alpha_i
			write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
			write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
			write(*,*)
!			write(*,*) 'conc_i', conc_i, 'conc_i_avg', conc_i_avg, 'conc_i_stdev', conc_i_stdev
!			read(*,*)
		end if

	else
		write(*,*) 'Error sinkEffSearch unrecognized'
	end if

end if !if(alpha_i_search .eqv. .TRUE.) then


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

			if(dabs(Residual_square) < 1.5d0*Residual_sqrdev .AND. dble(Residual) < Residual_stdev) then
				!we have converged for both values of alpha
				alpha_i_search=.FALSE.
				alpha_v_search=.FALSE.
				write(*,*) 'Both converged: alpha_i', alpha_i, 'alpha_v', alpha_v
				write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
				write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
!				read(*,*)

			else if(dabs(Residual_square) >= dabs(Residual_prev) &
				.AND. sign(Residual_square,Residual_prev)==Residual_square &
				.AND. v_iteration > 2) then !we have reached a local minimum in alpha_i
		
				alpha_v_search=.FALSE.
			
				write(*,*) 'Alpha_v converged', alpha_v_prev
				alpha_v=alpha_v_prev

				!This forces us to keep iterating until we have converged in one step for both v and i
				if(v_iteration > 1 .OR. i_iteration==0) then
					alpha_i_search=.TRUE.
					i_iteration=0
				end if
			
			else
				!calculate new value for alpha_v
				alpha_temp=alpha_v
			
				if(v_iteration==1) then
					if(alpha_v==0d0) then
						alpha_v=0.1d0
					else
						alpha_v=alpha_v/2d0
					end if
				else
					alpha_v=(Residual_square*alpha_v_prev-Residual_prev*alpha_v)/(Residual_square-Residual_prev)
				end if
			
				if(alpha_v < 0d0) then
					alpha_v=0d0
				end if
			
				if(alpha_v > 1d0) then
					alpha_v=1d0
				end if
			
				!store previous step
				alpha_v_prev=alpha_temp
				Residual_prev=Residual_square

			end if
		
			write(*,*) 'v_iteration', v_iteration
			write(*,*) 'alpha_v_old', alpha_v_prev, 'alpha_v', alpha_v
			write(*,*) 'Residual', Residual, 'Residual_stdev', Residual_stdev
			write(*,*) 'Residual squares', Residual_square, 'Residual squares dev', Residual_sqrdev
			write(*,*)
			!write(*,*) 'conc_v', conc_v, 'conc_v_avg', conc_v_avg, 'conc_v_stdev', conc_v_stdev
!			read(*,*)
		end if

	else
		write(*,*) 'Error sinkEffSearch unrecognized'
	end if

end if !if(alpha_v_search .eqv. .TRUE.) then

if(sinkEffSearch=='yes') then
	if(alpha_i_search .eqv. .FALSE.) then
		if(alpha_v_search .eqv. .FALSE.) then
			searchToggle=.FALSE.
			!we have converged for both alpha_v and alpha_i
			write(*,*) 'Converged for both alphas'
		else
			searchToggle=.TRUE.
		end if
	else
		!we have not yet converged for both alpha_v and alpha_i
		searchToggle=.TRUE.
	end if
end if

end do

if(sinkEffSearch=='yes') write(*,*) 'Final alpha_v', alpha_v, 'alpha_i', alpha_i

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

write(*,*) 'Finalizing processor', myProc%taskid

close(81)

call MPI_FINALIZE(ierr)		!must be performed at end of simulation to have no error

end program
