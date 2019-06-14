!***************************************************************************************************
!>subroutine outputDefects - outputs raw data for defect populations in various volume elements
!
!Outputs into file: rawdat.out. These contain the complete defect populations per volume element.
!
!Compiles data from local as well as global processors.
!***************************************************************************************************

subroutine outputDefects(elapsedTime,step)
use mod_srscd_constants
use DerivedType
implicit none
include 'mpif.h'

integer buffer, tag, i, j, k,l, status(MPI_STATUS_SIZE), numCellsRecv, numDefectsRecv
!integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount
integer defectCount, recvTemp(numSpecies+1)
type(defect), pointer :: defectCurrent
!double precision coordinatesRecv(3)
double precision, allocatable :: cellDefectSend(:,:)
double precision, allocatable :: cellDefectRecv(:,:)

double precision elapsedTime
integer step

tag=0

if(myProc%taskid==MASTER) then
	!first calculate DPA using MPI_ALLREDUCE
    write(82,*) 'elapsedTime', elapsedTime, '  step', step, 'dpa', DPA
	
	!Write defects in processor 0 to file
	write(82,*) 'processor', myProc%taskid
	do i=1,numCells
		defectCurrent=>defectList(i)%next
		write(82,*) 'cell', i,'globalCell', myMesh(i)%globalID
        write(82,*) 'Cu  ', 'V  ','SIA_m  ', 'SIA_im  ', 'num  '
		do while(associated(defectCurrent))
			write(82,*) defectCurrent%defectType, defectCurrent%num
			defectCurrent=>defectCurrent%next
		end do
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	do i=1,myProc%numtasks-1
		write(82,*) 'processor', i
		call MPI_RECV(numCellsRecv,1,MPI_INTEGER,i,99,MPI_COMM_WORLD,status,ierr)
		
		do j=1,numCellsRecv

            call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,100+2*j,MPI_COMM_WORLD,status,ierr)
            allocate(cellDefectRecv(numSpecies+1,numDefectsRecv+1))

            call MPI_RECV(cellDefectRecv,(numSpecies+1)*(numDefectsRecv+1),MPI_DOUBLE_PRECISION,i,&
                    101+2*j,MPI_COMM_WORLD,status,ierr)

            !!call MPI_RECV(coordinatesRecv,3,MPI_DOUBLE_PRECISION,i,100,MPI_COMM_WORLD,status,ierr)
            write(82,*) 'coordinates',(cellDefectRecv(l,1),l=1,3) , 'cell', j
            write(82,*) 'Cu  ', 'V  ', 'SIA_m  ', 'SIA_im  ', 'num'

            !!call MPI_RECV(cellDefectRecv,numDefectsRecv*(numSpecies+2),MPI_DOUBLE_PRECISION,i,102,MPI_COMM_WORLD,status,ierr)

			do k=2,numDefectsRecv+1
				!!call MPI_RECV(defectTypeRecv,numSpecies,MPI_INTEGER,i,102,MPI_COMM_WORLD,status,ierr)
				!!call MPI_RECV(cellNumberRecv,1,MPI_INTEGER,i,103,MPI_COMM_WORLD,status,ierr)
				!!call MPI_RECV(numRecv,1,MPI_INTEGER,i,104,MPI_COMM_WORLD,status,ierr)
				do l=1, numSpecies+1
					recvTemp(l) = cellDefectRecv(l,k)
				end do
				write(82,*) recvTemp
			end do

			deallocate(cellDefectRecv)
			
			!send a signal telling the other processor to go on to the next cell
			!call MPI_SEND(tag,1,MPI_INTEGER,i,105,MPI_COMM_WORLD,ierr)
		end do
	end do
    write(82,*)

else
		
	call MPI_SEND(numCells, 1, MPI_INTEGER, MASTER, 99, MPI_COMM_WORLD, ierr)
	
	do i=1,numCells
		!!call MPI_SEND(myMesh(i)%coordinates, 3, MPI_DOUBLE_PRECISION, MASTER, 100, MPI_COMM_WORLD, ierr)
		
		!calculate the number of defect types in cell i and send to MASTER.
		defectCount=0
		defectCurrent=>defectList(i)%next
		do while(associated(defectCurrent))
			defectCount=defectCount+1
			defectCurrent=>defectCurrent%next
		end do

		!send number of defects
		call MPI_SEND(defectCount,1,MPI_INTEGER,MASTER,100+2*i,MPI_COMM_WORLD, ierr)
		
		!send actual defect data (loop through defects a second time)
        allocate(cellDefectSend(numSpecies+1,defectCount+1))
        do k=1,3
            cellDefectSend(k,1) = myMesh(i)%coordinates(k)
        end do
        cellDefectSend(k+1,1) = 0
        j=1
		defectCurrent=>defectList(i)%next
		do while(associated(defectCurrent))
			!!call MPI_SEND(defectCurrent%defectType,numSpecies,MPI_INTEGER,MASTER,102,MPI_COMM_WORLD,ierr)
			!!call MPI_SEND(defectCurrent%cellNumber,1,MPI_INTEGER,MASTER,103,MPI_COMM_WORLD,ierr)
			!!call MPI_SEND(defectCurrent%num,1,MPI_INTEGER,MASTER,104,MPI_COMM_WORLD,ierr)
            j=j+1

            do k=1, numSpecies
                cellDefectSend(k,j) = defectCurrent%defectType(k)
            end do

            cellDefectSend(numSpecies+1,j) = defectCurrent%num

			defectCurrent=>defectCurrent%next
		end do

        call MPI_SEND(cellDefectSend, (numSpecies+1)*(defectCount+1), MPI_DOUBLE_PRECISION, MASTER, &
                101+2*i, MPI_COMM_WORLD, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,105,MPI_COMM_WORLD,status,ierr)

		deallocate(cellDefectSend)
	end do
end if



end subroutine

!***********************************************************************
!
!> Subroutine outputDefectsTotal() - outputs the total defects and post processing for the entire volume
!!
!! Outputs a list of all defects in system (defect type, num), regardless
!! of volume element. Compiles such a list using message passing between
!! processors. Used to find total defect populations in simulation volume,
!! not spatial distribution of defects.
!!
!! This subroutine will also automatically output the concentration of
!! vacancies, SIAs, and SIAs of size larger than 16, 20, and 24 defects.
!! (Corresponding to diameters 0.9 nm, 1.0 nm, 1.1 nm)
!
!***********************************************************************

subroutine outputDefectsTotal(elapsedTime, step)
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'

integer buffer, tag, i, j, k, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount, same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
!double precision defectsRecv(numSpecies+1), defectsSend(numSpecies+1)
double precision, allocatable :: defectsRecv(:,:)
double precision, allocatable :: defectsSend(:,:)
!Simulation variables passed by main program
double precision elapsedTime
integer step

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!total number of annihilation reactions
integer numAnnihilateTotal

!Variables for defect counting and concentration
!integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, HeNum, totalHe
integer VClusterNum, SIAClusterNum, LoopSize(3), totalVac, totalSIA, CuClusterNum, totalCu
double precision VMeanRadius, VoidMeanRadius, SIAMeanRadius, LoopMeanRadius, CuMeanRadius
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
end interface

tag=0

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

call MPI_ALLREDUCE(numAnnihilate,numAnnihilateTotal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myProc%taskid==MASTER) then
	!first calculate DPA using MPI_ALLREDUCE
    write(83,*) 'elapsedTime', elapsedTime, '  step', step, 'dpa', DPA
    write(84,*) 'elapsedTime', elapsedTime, '  step', step, 'dpa', DPA
	
	!Create list of defects in processor 0
	do i=1,numCells
		
		!ONLY OUTPUT DEFECTS IN THE BULK (NOT GBS)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
		
			!Do nothing, no output of defects in grain boundaries
			
		else
	
			defectCurrent=>defectList(i)%next	!the first defect in the list is 0
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	
	do i=1,myProc%numtasks-1
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,399,MPI_COMM_WORLD,status,ierr)
		
!		write(*,*) 'recieving ', numDefectsRecv, 'defects'
        !2019.05.14 add
		allocate(defectsRecv(numSpecies+1,numDefectsRecv))
        call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,400,MPI_COMM_WORLD,status,ierr)

		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,400,MPI_COMM_WORLD,status,ierr)
!			write(*,*) 'recieved defect ', j
			
			do k=1,numSpecies
				products(k)=defectsRecv(k,j)
			end do
			
			nullify(defectPrevList)
			
			defectCurrentList=>outputDefectList
			
			call findDefectInList(defectCurrentList, defectPrevList, products)
			
			!Next update defects
			if(associated(defectCurrentList)) then !if we aren't at the end of the list
				same=0
				
				do k=1,numSpecies
					if(defectCurrentList%defectType(k)==products(k)) then
						same=same+1
					endif
				end do
				
				if(same==numSpecies) then	
				
					!if the defect is already present in the list
				
					defectCurrentList%num=defectCurrentList%num+defectsRecv(numSpecies+1,j)
				
				else		
					
					!if the defect is to be inserted in the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
					defectPrevList%num=defectsRecv(numSpecies+1,j)
					
					do k=1,numSpecies
						defectPrevList%defectType(k)=products(k)
					end do
					
					!if inserted defect is in the middle of the list, point it to the next item in the list
					
					defectPrevList%next=>defectCurrentList
				endif
			else if(associated(defectPrevList)) then			
				
				!add a defect to the end of the list
				
				nullify(defectPrevList%next)
				allocate(defectPrevList%next)
				nullify(defectPrevList%next%next)
				defectPrevList=>defectPrevList%next
				allocate(defectPrevList%defectType(numSpecies))
				defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
				defectPrevList%num=defectsRecv(numSpecies+1,j)
				
				do k=1,numSpecies
					defectPrevList%defectType(k)=products(k)
				end do
				
			else
				
				write(*,*) 'error tried to insert defect at beginning of output defect list'

			end if

		end do
        !Signal to other processor to send the next defect
!        call MPI_SEND(tag,1,MPI_INTEGER, i, 405,MPI_COMM_WORLD, ierr)
		deallocate(defectsRecv)
		
	end do
	
	!Output defect list
!	write(*,*) 'Defects ', 'num'
    !Output totdat.out
	write(83,*) 'Cu  ','V  ','SIA_m  ','SIA_im  ','num'
	defectCurrentList=>outputDefectList
	
	!Initialize Defect counters
    CuClusterNum=0
    VClusterNum=0
    SIAClusterNum=0
    totalCu=0
	totalVac=0
	totalSIA=0
	totalVoid=0
	totalLoop=0

    !2019.05.14 Add
    VMeanRadius=0d0
    VoidMeanRadius=0d0
    SIAMeanRadius=0d0
    LoopMeanRadius=0d0
    CuMeanRadius=0d0
	
	do i=1,3
		LoopSize(i)=0
	end do
	VoidNum=0
	
	do while(associated(defectCurrentList))
		
		!Compile statistics for vacancy and SIA concentrations
		if(defectCurrentList%defectType(1) /= 0) then
			
			!HeNum=HeNum+defectCurrentList%num
            CuClusterNum = CuClusterNum + defectCurrentList%num
			
			!totalHe=totalHe+defectCurrentList%defectType(1)*defectCurrentList%num
            totalCu = totalCu + defectCurrentList%defectType(1) * defectCurrentList%num
		
		endif
		
		if(defectCurrentList%defectType(2) /= 0) then
		
			VClusterNum = VClusterNum + defectCurrentList%num
			
			totalVac = totalVac + defectCurrentList%defectType(2) * defectCurrentList%num
			
			if(defectCurrentList%defectType(2) >= 45) then
				VoidNum = VoidNum + defectCurrentList%num
				totalVoid = totalVoid + defectCurrentList%defectType(2) * defectCurrentList%num
			endif
		endif
		
		if(defectCurrentList%defectType(4) /= 0 .OR. defectCurrentList%defectType(3) /= 0) then

            SIAClusterNum = SIAClusterNum + defectCurrentList%num
			
			totalSIA = totalSIA + defectCurrentList%defectType(4) * defectCurrentList%num
			totalSIA = totalSIA + defectCurrentList%defectType(3) * defectCurrentList%num
		
			if(defectCurrentList%defectType(4) >= 16) then
				LoopSize(1)=LoopSize(1)+defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 20) then
				LoopSize(2)=LoopSize(2)+defectCurrentList%num
				totalLoop = totalLoop + defectCurrentList%defectType(4)*defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 24) then
				LoopSize(3)=LoopSize(3)+defectCurrentList%num
			endif
		endif
	
!		write(*,*) defectCurrentList%defectType, defectCurrentList%num
		write(83,*) defectCurrentList%defectType, defectCurrentList%num
		defectCurrentList=>defectCurrentList%next
	end do

    CuMeanRadius = (3*(dble(totalCu)/dble(CuClusterNum))*atomSize/(4*pi))**(1d0/3d0)
    VMeanRadius = (3*(dble(totalVac)/dble(VClusterNum))*atomSize/(4*pi))**(1d0/3d0)
    VoidMeanRadius = (3*(dble(totalVoid)/dble(VoidNum))*atomSize/(4*pi))**(1d0/3d0)
    LoopMeanRadius = ((dble(totalLoop)/dble(LoopSize(2)))*atomSize/(pi*burgers))**(1d0/2d0)

	
!This is now calculated in MeshReader.f90
!	systemVol=(myProc%globalCoord(2)-myProc%globalCoord(1))*&
!			  (myProc%globalCoord(4)-myProc%globalCoord(3))*&
!			  (myProc%globalCoord(6)-myProc%globalCoord(5))*1d-27
	
	
	write(83,*) 'CuClusterNum', CuClusterNum, 'number density (m-3)', dble(CuClusterNum)/(systemVol*1.0d-27), &
			'MeanRadius (nm)', CuMeanRadius
	write(83,*) 'VClusterNum', VClusterNum, 'number density (m-3)', dble(VClusterNum)/(systemVol*1d-27), &
			'VMeanRadius (nm)', VMeanRadius
	write(83,*)	'Voids', VoidNum, 'number density (m-3)', dble(VoidNum)/(systemVol*1d-27), &
			'VoidMeanRadius (nm)', VoidMeanRadius
	write(83,*) 'SIAClusterNum', SIAClusterNum, 'number density (m-3)', dble(SIAClusterNum)/(systemVol*1d-27)
	write(83,*) 'SIA_small (SIANum>=16)',LoopSize(1), 'number density (m-3)', dble(LoopSize(1))/(systemVol*1d-27)
	write(83,*) 'Loops (SIANum>=20)', LoopSize(2), 'number density (m-3)', dble(LoopSize(2))/(systemVol*1d-27), &
			'LoopMeanRadius (nm)', LoopMeanRadius
	write(83,*) 'SIA_large (SIANum>=24)', LoopSize(3), 'number density (m-3)', dble(LoopSize(3))/(systemVol*1d-27)
    write(83,*)

	!Output postpr.out
	write(84,*) 'CuClusterNum', CuClusterNum, 'number density (m-3)', dble(CuClusterNum)/(systemVol*1d-27), &
            'At.%', dble(CuClusterNum)*atomSize/systemVol
	write(84,*) 'VClusterNum', VClusterNum, 'number density (m-3)', dble(VClusterNum)/(systemVol*1d-27), &
            'At.%', dble(VClusterNum)*atomSize/systemVol
	write(84,*)	'Voids', VoidNum, 'number density (m-3)', dble(VoidNum)/(systemVol*1d-27), &
            'At.%', dble(VoidNum)*atomSize/systemVol
	write(84,*) 'SIAClusterNum', SIAClusterNum, 'number density (m-3)', dble(SIAClusterNum)/(systemVol*1d-27), &
            'At.%', dble(SIAClusterNum)*atomSize/systemVol
	write(84,*) 'SIA_small (SIANum>=16)',LoopSize(1), 'number density (m-3)', dble(LoopSize(1))/(systemVol*1d-27), &
            'At.%', dble(LoopSize(1))*atomSize/systemVol
	write(84,*) 'Loops (SIANum>=20)', LoopSize(2), 'number density (m-3)', dble(LoopSize(2))/(systemVol*1d-27), &
            'At.%', dble(LoopSize(2))*atomSize/systemVol
	write(84,*) 'SIA_large (SIANum>=24)', LoopSize(3), 'number density (m-3)', dble(LoopSize(3))/(systemVol*1d-27), &
            'At.%', dble(LoopSize(3))*atomSize/systemVol
    write(84,*)
	
	!*******************************************************************
	!Post processing: calculate relevant information (hard coded) to 
	!output into postprocessing.txt
	!*******************************************************************
	
	!Average vacancy size, average SIA size
	write(84,*) 'AverageCuSize', dble(totalCu)/dble(CuClusterNum)
	write(84,*) 'AverageVSize', dble(totalVac)/dble(VClusterNum)
	write(84,*) 'AverageSIASize', dble(totalSIA)/dble(SIAClusterNum)
	write(84,*) 'AverageVoidSize', dble(totalVoid)/dble(VoidNum)
	write(84,*) 'AverageLoopSize', dble(totalLoop)/dble(LoopSize(2))
    write(84,*)
	
	!Percent vacancies retained/annihilated
	!write(84,*) 'PercentHeRetained', dble(totalHe)/dble(numHeImplantTotal)
	write(84,*) 'PercentVRetained', dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImplantEvents))
	write(84,*) 'PercentVAnnihilated', dble(numAnnihilateTotal)/(dble(numDisplacedAtoms)*dble(totalImplantEvents))
    write(84,*)

	!Computation statistics: number of reactions in coarse/fine mesh, average timestep
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
	write(84,*) 'NumReactionsCoarse', reactionsCoarse
	write(84,*) 'NumReactionsFine', reactionsFine
	write(84,*) 'AverageTimestep', elapsedTime/dble(step)
	write(84,*)

else	!other processors
	numDefectsSend=0
	
	!Create list of defects in processor 
	do i=1,numCells
	
		!ONLY OUTPUT DEFECTS IN BULK (NOT GBs)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
			! Do nothing, no output of defects in grain boundaries
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						end if
					end do
					
					if(same == numSpecies) then
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else    !add a defect to the middle of the list
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
						
						numDefectsSend=numDefectsSend+1
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
					numDefectsSend=numDefectsSend+1
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 399, MPI_COMM_WORLD, ierr)

	defectCurrentList=>outputDefectList%next
	i=0
    !**************************************************
    !2019.05.14 modify
    allocate(defectsSend(numSpecies+1,numDefectsSend))

    do while(associated(defectCurrentList))
        i=i+1

        do j=1,numSpecies
            defectsSend(j,i)=defectCurrentList%defectType(j)
        end do

        defectsSend(numSpecies+1,i)=defectCurrentList%num

        defectCurrentList=>defectCurrentList%next

        if(i==numDefectsSend .AND. associated(defectCurrentList)) then
            write(*,*) 'error outputDefectList size does not match numDefectsSend'
        endif
    end do

    call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 400, MPI_COMM_WORLD, ierr)

	deallocate(defectsSend)

    !This just pauses the sending until the master has recieved all of the above information
!    call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,405,MPI_COMM_WORLD,status,ierr)
    !**************************************************

!	do while(associated(defectCurrentList))
!		i=i+1
	
!		do j=1,numSpecies
!			defectsSend(j)=defectCurrentList%defectType(j)
!		end do
		
!		defectsSend(numSpecies+1)=defectCurrentList%num
		
!		call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 400, MPI_COMM_WORLD, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
!		call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,405,MPI_COMM_WORLD,status,ierr)
		
!		defectCurrentList=>defectCurrentList%next
		
!		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
!			write(*,*) 'error outputDefectList size does not match numDefectsSend'
!		endif
		
!		write(*,*) 'sent defect', i
!	end do


	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
end if

!Deallocate memory

defectCurrentList=>outputDefectList
nullify(defectPrevList)

do while(associated(defectCurrentList))
	defectPrevList=>defectCurrentList
	defectCurrentList=>defectCurrentList%next
	
	if(allocated(defectPrevList%defectType)) deallocate(defectPrevList%defectType)
	
	deallocate(defectPrevList)
end do

nullify(defectCurrentList)
nullify(outputDefectList)

end subroutine

!***********************************************************************
!
!> Subroutine outputDefectsBoundary() - outputs the total defects and post processing on the grain boundary
!!
!! Outputs a list of all defects in the grain boundary (defect type, num), regardless
!! of volume element. Compiles such a list using message passing between
!! processors. Used to find total defect populations in simulation volume,
!! not spatial distribution of defects.
!!
!! This subroutine will also automatically output the concentration of
!! vacancies, SIAs, and SIAs of size larger than 16, 20, and 24 defects.
!! (Corresponding to diameters 0.9 nm, 1.0 nm, 1.1 nm)
!
!***********************************************************************

subroutine outputDefectsBoundary(elapsedTime, step)
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'

integer buffer, tag, i, j, k, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount, same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
!double precision defectsRecv(numSpecies+1), defectsSend(numSpecies+1)

double precision, allocatable :: defectsRecv(:,:), defectsSend(:,:)

double precision atomArea, systemArea

!Simulation variables passed by main program
double precision elapsedTime
integer step

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!total number of annihilation reactions
integer numAnnihilateTotal

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
end interface

tag=0

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

call MPI_ALLREDUCE(numAnnihilate,numAnnihilateTotal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myProc%taskid==MASTER) then
	!first calculate DPA using MPI_ALLREDUCE
	
	write(83,*) 'dpa', DPA
	
	!Create list of defects in processor 0
	do i=1,numCells
		
		!ONLY OUTPUT DEFECTS IN THE GB (NOT BULK)
		if(numMaterials > 1 .AND. myMesh(i)%material == 1) then
		
			!Do nothing, no output of defects in grain boundaries
			
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	
	do i=1,myProc%numtasks-1
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,399,MPI_COMM_WORLD,status,ierr)
		
!		write(*,*) 'recieving ', numDefectsRecv, 'defects'
		allocate(defectsRecv(numSpecies+1,numDefectsRecv))
		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,400,MPI_COMM_WORLD,status,ierr)
		
		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,400,MPI_COMM_WORLD,status,ierr)
!			write(*,*) 'recieved defect ', j
			
			do k=1,numSpecies
				products(k)=defectsRecv(k,j)
			end do
			
			nullify(defectPrevList)
			
			defectCurrentList=>outputDefectList
			
			call findDefectInList(defectCurrentList, defectPrevList, products)
			
			!Next update defects
			if(associated(defectCurrentList)) then !if we aren't at the end of the list
				same=0
				
				do k=1,numSpecies
					if(defectCurrentList%defectType(k)==products(k)) then
						same=same+1
					endif
				end do
				
				if(same==numSpecies) then	
				
					!if the defect is already present in the list
				
					defectCurrentList%num=defectCurrentList%num+defectsRecv(numSpecies+1,j)
				
				else		
					
					!if the defect is to be inserted in the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
					defectPrevList%num=defectsRecv(numSpecies+1,j)
					
					do k=1,numSpecies
						defectPrevList%defectType(k)=products(k)
					end do
					
					!if inserted defect is in the middle of the list, point it to the next item in the list
					
					defectPrevList%next=>defectCurrentList
				endif
			else if(associated(defectPrevList)) then			
				
				!add a defect to the end of the list
				
				nullify(defectPrevList%next)
				allocate(defectPrevList%next)
				nullify(defectPrevList%next%next)
				defectPrevList=>defectPrevList%next
				allocate(defectPrevList%defectType(numSpecies))
				defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
				defectPrevList%num=defectsRecv(numSpecies+1,j)
				
				do k=1,numSpecies
					defectPrevList%defectType(k)=products(k)
				end do
				
			else
				
				write(*,*) 'error tried to insert defect at beginning of output defect list'

			endif
			
			!Signal to other processor to send the next defect
			!call MPI_SEND(tag,1,MPI_INTEGER, i, 405,MPI_COMM_WORLD, ierr)
		end do
		deallocate(defectsRecv)
		
	end do
	
	!Output defect list
!	write(*,*) 'Defects ', 'num'
	write(83,*) 'Defects ', 'num'
	defectCurrentList=>outputDefectList
	
	!Initialize Defect counters
	CuNum=0
	VNum=0
	SIANum=0
	totalCu=0
	totalVac=0
	totalSIA=0
	totalVoid=0
	totalLoop=0
	
	do i=1,3
		LoopSize(i)=0
	end do
	VoidNum=0
	
	do while(associated(defectCurrentList))
		
		!Compile statistics for vacancy and SIA concentrations
		if(defectCurrentList%defectType(1) /= 0) then
			
			CuNum=CuNum+defectCurrentList%num
			
			totalCu=totalCu+defectCurrentList%defectType(1)*defectCurrentList%num
		
		endif
		
		if(defectCurrentList%defectType(2) /= 0) then
		
			VNum=VNum+defectCurrentList%num
			
			totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
			
			if(defectCurrentList%defectType(2) >= 45) then
				VoidNum=VoidNum+defectCurrentList%num
				totalVoid = totalVoid + defectCurrentList%defectType(2)*defectCurrentList%num
			endif
		endif
		
		if(defectCurrentList%defectType(4) /= 0 .OR. defectCurrentList%defectType(3) /= 0) then
		
			SIANum=SIANum+defectCurrentList%num
			
			totalSIA=totalSIA+defectCurrentList%defectType(4)*defectCurrentList%num
			totalSIA=totalSIA+defectCurrentList%defectType(3)*defectCurrentList%num
		
			if(defectCurrentList%defectType(4) >= 16) then
				LoopSize(1)=LoopSize(1)+defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 20) then
				LoopSize(2)=LoopSize(2)+defectCurrentList%num
				totalLoop = totalLoop + defectCurrentList%defectType(4)*defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 24) then
				LoopSize(3)=LoopSize(3)+defectCurrentList%num
			endif
		endif
	
!		write(*,*) defectCurrentList%defectType, defectCurrentList%num
		write(83,*) defectCurrentList%defectType, defectCurrentList%num
		defectCurrentList=>defectCurrentList%next
	end do
	
!This is now calculated in MeshReader.f90
!	systemVol=(myProc%globalCoord(2)-myProc%globalCoord(1))*&
!			  (myProc%globalCoord(4)-myProc%globalCoord(3))*&
!			  (myProc%globalCoord(6)-myProc%globalCoord(5))*1d-27

	atomArea=(9d0*pi*atomsize**2d0/16d0)**(1d0/3d0)
	systemArea=(myProc%globalCoord(4)-myProc%globalCoord(3))*&
				(myProc%globalCoord(6)-myProc%globalCoord(5))
	
	
	write(83,*) 'CuNum', CuNum, 'Conc', dble(CuNum)*atomArea/systemArea
	write(83,*) 'VNum', VNum, 'Conc', dble(VNum)*atomArea/systemArea
	write(83,*)	'Voids', VoidNum, 'Conc', dble(VoidNum)*atomArea/systemArea
	write(83,*) 'SIANum', SIANum, 'Conc', dble(SIANum)*atomArea/systemArea
!	write(83,*) 'SIA_small',LoopSize(1), 'Conc', dble(LoopSize(1))/systemVol
	write(83,*) 'Loops', LoopSize(2), 'Conc', dble(LoopSize(2))*atomArea/systemArea
!	write(83,*) 'SIA_large', LoopSize(3), 'Conc', dble(LoopSize(3))/systemVol
	
	write(84,*) 'CuNum', CuNum, 'Conc', dble(CuNum)*atomArea/systemArea
	write(84,*) 'VNum', VNum, 'Conc', dble(VNum)*atomArea/systemArea
	write(84,*)	'Voids', VoidNum, 'Conc', dble(VoidNum)*atomArea/systemArea
	write(84,*) 'SIANum', SIANum, 'Conc', dble(SIANum)*atomArea/systemArea
!	write(84,*) 'SIA_small',LoopSize(1), 'Conc', dble(LoopSize(1))/systemVol
	write(84,*) 'Loops', LoopSize(2), 'Conc', dble(LoopSize(2))*atomArea/systemArea
!	write(84,*) 'SIA_large', LoopSize(3), 'Conc', dble(LoopSize(3))/systemVol
	
	!*******************************************************************
	!Post processing: calculate relevant information (hard coded) to 
	!output into postprocessing.txt
	!*******************************************************************
	
	!Average vacancy size, average SIA size
	write(84,*) 'AverageCuSize', dble(totalCu)/dble(CuNum)
	write(84,*) 'AverageVSize', dble(totalVac)/dble(VNum)
	write(84,*) 'AverageSIASize', dble(totalSIA)/dble(SIANum)
	write(84,*) 'AverageVoidSize', dble(totalVoid)/dble(VoidNum)
	write(84,*) 'AverageLoopSize', dble(totalLoop)/dble(LoopSize(2))
	
	!Percent vacancies retained/annihilated
!	write(84,*) 'PercentHeRetained', dble(totalHe)/dble(numHeImplantTotal)
	write(84,*) 'PercentVRetained', dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImplantEvents))
	write(84,*) 'PercentVAnnihilated', dble(numAnnihilateTotal)/(dble(numDisplacedAtoms)*dble(totalImplantEvents))
	
	!Computation statistics: number of reactions in coarse/fine mesh, average timestep
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
	write(84,*) 'NumReactionsCoarse', reactionsCoarse
	write(84,*) 'NumReactionsFine', reactionsFine
	write(84,*) 'AverageTimestep', elapsedTime/dble(step)
	
	!Leave blank line
	write(83,*)
	write(84,*)

else
	numDefectsSend=0
	
	!Create list of defects in processor 
	do i=1,numCells
	
		!ONLY OUTPUT DEFECTS IN GBs (NOT BULK)
		if(numMaterials > 1 .AND. myMesh(i)%material == 1) then
			! Do nothing, no output of defects in grain boundaries
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
						
						numDefectsSend=numDefectsSend+1
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
					numDefectsSend=numDefectsSend+1
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 399, MPI_COMM_WORLD, ierr)
	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		
		!call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 400, MPI_COMM_WORLD, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,405,MPI_COMM_WORLD,status,ierr)
		
		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		endif
		
!		write(*,*) 'sent defect', i
	end do
	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 400, MPI_COMM_WORLD, ierr)
	deallocate(defectsSend)
	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
end if

!Deallocate memory

defectCurrentList=>outputDefectList
nullify(defectPrevList)

do while(Associated(defectCurrentList))
	defectPrevList=>defectCurrentList
	defectCurrentList=>defectCurrentList%next
	
	if(allocated(defectPrevList%defectType)) deallocate(defectPrevList%defectType)
	
	deallocate(defectPrevList)
end do

nullify(defectCurrentList)
nullify(outputDefectList)

end subroutine

!***************************************************************************************************
!
!> Subroutine OutputDefectsProfile(sim) - outputs depth profile of defect concentrations
!!
!! This subroutine outputs a list of the z-coordinates in the given system as well as the total
!! concentration of vacancies, interstitials, and helium atoms at each depth. This is meant to 
!! be used by the Gnuplot script DepthProfile.p.
!!
!! Inputs: Defect and mesh information, simulation number (for filename)
!! Outputs: Text document with defect profile
!
!***************************************************************************************************

subroutine outputDefectsProfile(sim)
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'

integer, allocatable :: defectProfileArray(:,:)
integer, allocatable :: defectNumArray(:,:)
integer, allocatable :: tempProfileArray(:,:)
integer, allocatable :: tempNumArray(:,:)
integer zEntry, xEntry, yEntry
integer i, j, k, procID, status(MPI_STATUS_SIZE), sim
type(defect), pointer :: defectCurrent
type(defect), pointer :: defectProfile(:)
type(defect), pointer :: defectListCurrent, defectListPrev
double precision zCoord, xCoord, yCoord
character*20 filename
integer same
integer numDefectsSend, numSend
integer numDefectsRecv, numRecv, defectTypeStore(numSpecies)

integer, allocatable :: buffer(:,:)

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
end interface

!Create an array with z-coordinates given by the z-coordinates of the mesh elements

!First calculate the number of z-coordinates in the entire system

!Assuming uniform mesh (if nonuniform, we cannot use myMesh(1)%length as the length of every element
!numZ=int(((myProc%globalCoord(6)-myProc%globalCoord(5))/myMesh(1)%length))
!numX=int(((myProc%globalCoord(2)-myProc%globalCoord(1))/myMesh(1)%length))
!numY=int(((myProc%globalCoord(4)-myProc%globalCoord(3))/myMesh(1)%length))

!Allocate and initialize array of defect numbers in each processor
allocate(defectProfileArray(numZ,numSpecies))
allocate(defectNumArray(numZ,numSpecies))
allocate(defectProfile(numZ))

allocate(tempProfileArray(numZ,numSpecies))
allocate(tempNumArray(numZ,numSpecies))

do i=1,numz
	allocate(defectProfile(i)%defectType(numSpecies))
	do j=1,numSpecies
		defectProfileArray(i,j)=0
		defectNumArray(i,j)=0
		defectProfile(i)%defectType(j)=0
	end do
	defectProfile(i)%num=0
	defectProfile(i)%cellNumber=i	!Here cellNumber will be used to store the x-coordinate of this cell, not the actual cell ID number
	nullify(defectProfile(i)%next)
end do

!Go through defects in every mesh element and add them to defectProfileArray
do i=1,numCells
	
	!This tells us which index in defectProfileArray we should be adding to
	zEntry=int((myMesh(i)%coordinates(3)+myMesh(i)%length/2d0)/myMesh(i)%length)
	!xEntry=int((myMesh(i)%coordinates(1)+myMesh(i)%length/2d0)/myMesh(i)%length)
	
	defectCurrent=>defectList(i)
	do while(associated(defectCurrent))
		
		!Add the number of defects of each species to defectProfileArray
		do j=1,numSpecies
			
			!For now, only count vacancies with n .GE. 5 (used to count stable proto-voids)
			if(j==2 .AND. defectCurrent%defectType(j) >= 5) then
				defectProfileArray(zEntry,j)=defectProfileArray(zEntry,j)+defectCurrent%defectType(j)*defectCurrent%num
				if(defectCurrent%defectType(j) > 0) then
					defectNumArray(zEntry, j)=defectNumArray(zEntry,j)+defectCurrent%num
				end if
			end if
			
		end do
		
		!Add the current defect to defectProfile
		defectListCurrent=>defectProfile(zEntry)
		nullify(defectListPrev)
		
		call findDefectInList(defectListCurrent, defectListPrev, defectCurrent%defectType)
		
		!Next update defects
		if(associated(defectListCurrent)) then !if we aren't at the end of the list
			same=0
			
			do j=1,numSpecies
				if(defectListCurrent%defectType(j)==defectCurrent%defectType(j)) then
					same=same+1
				endif
			end do
			
			if(same==numSpecies) then	
			
				!if the defect is already present in the list
			
				defectListCurrent%num=defectListCurrent%num+defectCurrent%num
			
			else		
				
				!if the defect is to be inserted in the list
				
				nullify(defectListPrev%next)
				allocate(defectListPrev%next)
				nullify(defectListPrev%next%next)
				defectListPrev=>defectListPrev%next
				allocate(defectListPrev%defectType(numSpecies))
				defectListPrev%cellNumber=zEntry
				defectListPrev%num=defectCurrent%num
				
				do j=1,numSpecies
					defectListPrev%defectType(j)=defectCurrent%defectType(j)
				end do
				
				!if inserted defect is in the middle of the list, point it to the next item in the list
				
				defectListPrev%next=>defectListCurrent
			endif
		else if(associated(defectListPrev)) then			
			
			!add a defect to the end of the list
			
			nullify(defectListPrev%next)
			allocate(defectListPrev%next)
			nullify(defectListPrev%next%next)
			defectListPrev=>defectListPrev%next
			allocate(defectListPrev%defectType(numSpecies))
			defectListPrev%cellNumber=zEntry !no need for cell numbers in outputDefectList
			defectListPrev%num=defectCurrent%num
			
			do j=1,numSpecies
				defectListPrev%defectType(j)=defectCurrent%defectType(j)
			end do
			
		else
			
			write(*,*) 'error tried to insert defect at beginning of output defect list'

		endif
		
		defectCurrent=>defectCurrent%next

	end do

end do
	
!Add up all defect numbers for each processor
if(myProc%taskid==MASTER) then

	do procID=1,myProc%numTasks-1

		call MPI_RECV(tempProfileArray,numz*numSpecies,MPI_INTEGER, procID, 700, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(tempNumArray,numz*numSpecies,MPI_INTEGER, procID, 600,MPI_COMM_WORLD, status, ierr)
		
		do i=1,numz

			!call MPI_RECV(tempArray,numSpecies,MPI_INTEGER, procID, i*700, MPI_COMM_WORLD, status, ierr)
			
			!Add defects of neighboring procs to master processor defect array
			do j=1,numSpecies
				defectProfileArray(i,j)=defectProfileArray(i,j)+tempProfileArray(i,j)
			end do
			
			!call MPI_RECV(tempArray,numSpecies,MPI_INTEGER, procID, i*600,MPI_COMM_WORLD, status, ierr)
			
			do j=1,numSpecies
				defectNumArray(i,j)=defectNumArray(i,j)+tempNumArray(i,j)
			end do
			
			call MPI_RECV(numDefectsRecv, 1, MPI_INTEGER, procID, i*800, MPI_COMM_WORLD, status, ierr)

			allocate(buffer(numDefectsRecv,numSpecies+1))
			call MPI_RECV(buffer, numDefectsRecv*(numSpecies+1), MPI_INTEGER, procID, i*900, MPI_COMM_WORLD, status, ierr)
			
			do j=1,numDefectsRecv
			
				!call MPI_RECV(buffer, numSpecies+1, MPI_INTEGER, procID, i*900+j, MPI_COMM_WORLD, status, ierr)
				
				do k=1,numSpecies
					defectTypeStore(k)=buffer(j,k)
				end do
				
				!Add recieved defect to list of defects in zCoord
				defectListCurrent=>defectProfile(i)
				nullify(defectListPrev)
		
				call findDefectInList(defectListCurrent, defectListPrev, defectTypeStore)
				
				!Next update defects
				if(associated(defectListCurrent)) then !if we aren't at the end of the list
					same=0
					
					do k=1,numSpecies
						if(defectListCurrent%defectType(k)==defectTypeStore(k)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectListCurrent%num=defectListCurrent%num+buffer(j,numSpecies+1)
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectListPrev%next)
						allocate(defectListPrev%next)
						nullify(defectListPrev%next%next)
						defectListPrev=>defectListPrev%next
						allocate(defectListPrev%defectType(numSpecies))
						defectListPrev%cellNumber=i
						defectListPrev%num=buffer(j,numSpecies+1)
						
						do k=1,numSpecies
							defectListPrev%defectType(k)=defectTypeStore(k)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectListPrev%next=>defectListCurrent
					endif
				else if(associated(defectListPrev)) then			
					
					!add a defect to the end of the list
					
					nullify(defectListPrev%next)
					allocate(defectListPrev%next)
					nullify(defectListPrev%next%next)
					defectListPrev=>defectListPrev%next
					allocate(defectListPrev%defectType(numSpecies))
					defectListPrev%cellNumber=i !no need for cell numbers in outputDefectList
					defectListPrev%num=buffer(j,numSpecies+1)
					
					do k=1,numSpecies
						defectListPrev%defectType(k)=defectTypeStore(k)
					end do
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
		
				endif
				
			
			end do

			deallocate(buffer)
		
		end do
	
	end do
	
	!Outputs defectProfileArray to a file

	filename(1:12)='DepthProfile'
	write(unit=filename(13:14), fmt='(I2)') sim
	filename(15:18)='.out'
	
	open(99, file=filename, action='write', status='Unknown')
	
	!write(99,*) 'XCoord ', 'HeConc ', 'VConc ', 'SIAMobileConc ', 'SIASessileConc'
	
	do i=1,numz
		!XCoord=(i-1)*myMesh(1)%length+myMesh(1)%length/2d0
		ZCoord=(i-1)*myMesh(1)%length+myMesh(1)%length/2d0
		
		!write(99,*) XCoord, (dble(defectProfileArray(i,j)/(numY*numz*((1d-9*myMesh(1)%length)**3d0))),&
		!	j=1,3)
		
		!Vacancy concen
		write(99,*) ZCoord, dble(defectProfileArray(i,2)/(numy*numx*((1d-9*myMesh(1)%length)**3d0))), &
			dble(defectNumArray(i,2)/(numy*numx*((1d-9*myMesh(1)%length)**3d0)))
			
		defectListCurrent=>defectProfile(i)
		do while(associated(defectListCurrent))
			
			write(99,*) defectListCurrent%defectType, defectListCurrent%num
			defectListCurrent=>defectListCurrent%next
			
		end do
		
		write(99,*)
	end do
	
	close(99)
	
else
	
	!Send local array to be added to master processor array
	call MPI_SEND(defectProfileArray,numz*numSpecies,MPI_INTEGER, MASTER, 700, MPI_COMM_WORLD, ierr)
	call MPI_SEND(defectNumArray, numz*numSpecies, MPI_INTEGER, MASTER, 600, MPI_COMM_WORLD, ierr)

	do i=1,numz
		!call MPI_SEND(defectProfileArray(i,:),numSpecies,MPI_INTEGER, MASTER, i*700, MPI_COMM_WORLD, ierr)
		
		!call MPI_SEND(defectNumArray(i,:), numSpecies, MPI_INTEGER, MASTER, i*600, MPI_COMM_WORLD, ierr)
		
		!Compute how many defect types are in this list
		numDefectsSend=0
		defectCurrent=>defectProfile(i)%next
		do while(associated(defectCurrent))
			numDefectsSend=numDefectsSend+1
			defectCurrent=>defectCurrent%next
		end do
		
		!Send how many defect types are going to be sent
		call MPI_SEND(numDefectsSend,1,MPI_INTEGER, MASTER, i*800, MPI_COMM_WORLD, ierr)
		
		!Send each defect
		defectCurrent=>defectProfile(i)%next
		numSend=0

		allocate(buffer(numDefectsSend,numSpecies+1))
		
		do while(associated(defectCurrent))
			
			!Create buffer of information to send
			numSend=numSend+1
			do j=1,numSpecies
				buffer(numSend,j)=defectCurrent%defectType(j)
			end do
			buffer(numSend,numSpecies+1)=defectCurrent%num
			
			!Send information on this defect type
			!call MPI_SEND(buffer, numSpecies+1, MPI_INTEGER, MASTER, i*900+numSend, MPI_COMM_WORLD, ierr)
			
			defectCurrent=>defectCurrent%next
		end do
		call MPI_SEND(buffer, numDefectsSend*(numSpecies+1), MPI_INTEGER, MASTER, i*900, MPI_COMM_WORLD, ierr)
		deallocate(buffer)
		
	end do

end if

deallocate(defectProfileArray)
deallocate(defectNumArray)

deallocate(tempProfileArray)
deallocate(tempNumArray)

!Deallocate defect list
do i=1,numz
	deallocate(defectProfile(i)%defectType)
	defectListCurrent=>defectProfile(i)%next

	do while(associated(defectListCurrent))
		defectListPrev=>defectListCurrent
		defectListCurrent=>defectListCurrent%next
		deallocate(defectListPrev%defectType)
		deallocate(defectListPrev)
	end do
end do

deallocate(defectProfile)

end subroutine

!***********************************************************************
!
!> Subroutine outputDefectsXYZ - outputs xyz file of defect concentrations
!!
!! Outputs number of vacancies / SIAs per volume element in .XYZ file
!! (see wikipedia page for documentation on file format)
!!
!! Only outputs defects in the main mesh (ignores fine meshes present)
!
!***********************************************************************

subroutine outputDefectsXYZ(elapsedTime,step)
use mod_srscd_constants
use DerivedType
implicit none
include 'mpif.h'

integer CuCount, VCount, SIACount, numPoints
integer numCellsNeighbor
integer i,k
integer status(MPI_STATUS_SIZE)
!double precision coords(3)
double precision, allocatable ::  xyzRecv(:)
double precision, allocatable ::  xyzSend(:)
integer points(myProc%numTasks), numRecv(myProc%numTasks)

type(defect), pointer :: defectCurrent

double precision elapsedTime
integer step

call MPI_ALLREDUCE(numCells, numPoints, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_GATHERV(numCells,1,MPI_INTEGER,points,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD, ierr)

do i=0, myProc%numTasks-1
	numRecv(i+1)=points(i)*9
end do

allocate(xyzSend(9*numCells))
allocate(xyzRecv(9*numPoints))

do i=1,numCells

	k=(i-1)*9
	xyzSend(k+1) = myMesh(i)%coordinates(1)
	xyzSend(k+2) = myMesh(i)%coordinates(2)
	xyzSend(1+3) = myMesh(i)%coordinates(3)

	xyzSend(k+4) = 0d0  !CuCount
	xyzSend(k+5) = 0d0  !VCount
	xyzSend(k+6) = 0d0  !SIACount

	defectCurrent=>defectList(i)
	do while(associated(defectCurrent))
		xyzSend(k+4) = xyzSend(k+4)+defectCurrent%defectType(1)*defectCurrent%num
		xyzSend(k+5) = xyzSend(k+5)+defectCurrent%defectType(2)*defectCurrent%num
		xyzSend(k+6) = xyzSend(k+6)+defectCurrent%defectType(3)*defectCurrent%num + &
				defectCurrent%defectType(4)*defectCurrent%num

		xyzSend(k+7) = xyzSend(k+4)/(meshLength**3d0)	!Cu concentration
		xyzSend(k+8) = xyzSend(k+5)/(meshLength**3d0)	!V concentration
		xyzSend(k+9) = xyzSend(k+6)/(meshLength**3d0)	!SIA concentration
		defectCurrent=>defectCurrent%next
	end do

end do

call MPI_GATHERV(xyzSend,numRecv,MPI_DOUBLE_PRECISION,xyzRecv,numRecv,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD, ierr)

if(myProc%taskid==MASTER) then
	write(87,*) 'elapsedTime', elapsedTime, '  step', step, 'dpa', DPA
	write(87,*) 'numPoints', numPoints
	write(87,*) 'CellNumber ', 'x ', 'y ', 'z ', 'NumCu ', 'NumV ', 'NumSIA', &
			'CuCon(nm-3) ', 'VCon(nm-3) ', 'SIACon(nm-3)'

	do i=1, numPoints
		k=(i-1)*9
		CuCount=xyzRecv(k+4)
		VCount=xyzRecv(k+5)
		SIACount=xyzRecv(k+6)

		write(87,*) i, xyzRecv(k+1), xyzRecv(k+2), xyzRecv(k+3), CuCount, VCount, SIACount, &
				xyzRecv(k+7),xyzRecv(k+8),xyzRecv(k+9)
	end do
end if

deallocate(xyzSend)
deallocate(xyzRecv)

 !59		format(i4,3(1x,e12.4),3(1x,i6))

end subroutine

!***********************************************************************
!> Subroutine outputDefectsVTK - outputs defect populations in vtk file format
!
! Outputs defects in a VTK formatted file (3D unstructured grid of linear cubes)
!***********************************************************************

subroutine outputDefectsVTK(fileNumber)	!fileNumber = outputCounter
use mod_srscd_constants
use DerivedType
implicit none
include 'mpif.h'

integer fileNumber
integer i,j,k
integer xindex, yindex, zindex
integer numV, numCu, numSIA, GrainID
integer numCellsNeighbor, cellIndex(3)
integer status(MPI_STATUS_SIZE)

!integer, allocatable :: DefectListVTK(:,:,:,:)
!integer, allocatable :: GlobalGrainID(:,:,:)
integer, allocatable :: DefectListGIDVTK(:,:,:,:)
!2019.05.15 Add
integer , allocatable :: defectSend(:,:)
integer, allocatable :: defectRecv(:,:)

double precision xlength, ylength, zlength
type(defect), pointer :: defectCurrent
character(13) :: fileName

if(myProc%taskid==MASTER) then

	!Step 1: Find the dimensions of the (global) mesh (assume cubic elements)

!	numx=(myProc%globalCoord(2)-myProc%globalCoord(1))/myMesh(1)%length
!	numy=(myProc%globalCoord(4)-myProc%globalCoord(3))/myMesh(1)%length
!	numz=(myProc%globalCoord(6)-myProc%globalCoord(5))/myMesh(1)%length

	xlength=myMesh(1)%length
	ylength=myMesh(1)%length
	zlength=myMesh(1)%length

	!Step 2: Allocate DefectListVTK and GlobalGrainID to the same dimensions as global mesh

	allocate(DefectListGIDVTK(numx,numy,numz,4))
	!allocate(GlobalGrainID(numx,numy,numz))

	!Step 3: Fill in DefectListVTK

	do i=1,numCells
	
		!Step 3.1: identify the x,y,z index of this cell
		xindex=((myMesh(i)%coordinates(1)-xlength/2)-myProc%globalCoord(1))/myMesh(1)%length+1
		yindex=((myMesh(i)%coordinates(2)-ylength/2)-myProc%globalCoord(3))/myMesh(1)%length+1
		zindex=((myMesh(i)%coordinates(3)-zlength/2)-myProc%globalCoord(5))/myMesh(1)%length+1
	
		!Step 3.2: count number of vacancies, interstitials, and He atoms in this cell
		defectCurrent=>defectList(i)
	
		numV=0
		numSIA=0
		numCu=0
		GrainID=myMesh(i)%material
	
		do while(associated(defectCurrent))
			numV=numV+defectCurrent%defectType(2)*defectCurrent%num
			numCu=numCu+defectCurrent%defectType(1)*defectCurrent%num
			numSIA=numSIA+defectCurrent%defectType(3)*defectCurrent%num + &
				defectCurrent%defectType(4)*defectCurrent%num
		
			defectCurrent=>defectCurrent%next
		end do
	
		!If we are in a grain boundary or nanoporous simulation, only output defects in the bulk
		if(numMaterials > 1 .AND. GrainID /= 1) then
			numV=0
			numSIA=0
			numCu=0
		end if
	
		!Step 3.3: enter in data point in DefectListVTK
		DefectListGIDVTK(xindex,yindex,zindex,1)=numCu
		DefectListGIDVTK(xindex,yindex,zindex,2)=numV
		DefectListGIDVTK(xindex,yindex,zindex,3)=numSIA
		DefectListGIDVTK(xindex,yindex,zindex,4)=GrainID
		!GlobalGrainID(xindex,yindex,zindex)=GrainID
	
	end do

	!Step 3.4: recieve data from other processors and fill into DefectListVTK
	do i=1,myProc%numTasks-1
	
		call MPI_RECV(numCellsNeighbor, 1, MPI_INTEGER, i, 567, MPI_COMM_WORLD, status, ierr)
		allocate(defectRecv(7,numCellsNeighbor))
		call MPI_RECV(defectRecv, 7*numCellsNeighbor, MPI_INTEGER, i, 571, MPI_COMM_WORLD, status, ierr)

		do j=1,numCellsNeighbor
	
			!call MPI_RECV(cellIndex, 3, MPI_INTEGER, i, 571, MPI_COMM_WORLD, status, ierr)
			!call MPI_RECV(numHe, 1, MPI_INTEGER, i, 568, MPI_COMM_WORLD, status, ierr)
			!call MPI_RECV(numV, 1, MPI_INTEGER, i, 569, MPI_COMM_WORLD, status, ierr)
			!call MPI_RECV(numSIA, 1, MPI_INTEGER, i, 570, MPI_COMM_WORLD, status, ierr)
			!call MPI_RECV(GrainID, 1, MPI_INTEGER, i, 571, MPI_COMM_WORLD, status, ierr)
		
			DefectListGIDVTK(defectRecv(1,j),defectRecv(2,j),defectRecv(3,j),1)=defectRecv(4,j)
			DefectListGIDVTK(defectRecv(1,j),defectRecv(2,j),defectRecv(3,j),2)=defectRecv(5,j)
			DefectListGIDVTK(defectRecv(1,j),defectRecv(2,j),defectRecv(3,j),3)=defectRecv(6,j)
			DefectListGIDVTK(defectRecv(1,j),defectRecv(2,j),defectRecv(3,j),4)=defectRecv(7,j)
			!GlobalGrainID(defectRecv(1,j),defectRecv(2,j),defectRecv(3,j))=defectRecv(7,j)
				
		end do
		deallocate(defectRecv)
	
	end do

	!Step 4: Ouptut data in VTK file format, taking data from DefectListVTK

	fileName(1:6)='VTKout'
	write(unit=fileName(7:9), fmt='(I3)') fileNumber
	fileName(10:13)='.vtk'

	open(88, file=fileName, status='Unknown')

	write(88,'(a)') '# vtk DataFile Version 1.0'
	write(88,'(a)') '3D Unstructured Grid of Linear Cubes'
	write(88,'(a)') 'ASCII'
	write(88,*)
	write(88,'(a)') 'DATASET STRUCTURED_POINTS'
	write(88,'(a)', advance='no') 'DIMENSIONS '
	write(88,'(I4,I4,I4)') numx, numy, numz
	write(88,'(a)', advance='no') 'ORIGIN '
	write(88,'(I4,I4,I4)') 0, 0, 0
	write(88,'(a)', advance='no') 'SPACING '
	write(88,*) xlength, ylength, zlength
	write(88,'(a)', advance='no') 'POINT_DATA '
	write(88,*) numx*numy*numz

	write(88,'(a,a,a,I4)') 'SCALARS ', 'Cu', ' float', 1
	write(88,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'

	do k=1,numz
		do j=1,numy
			do i=1,numx

				write(88,*) DefectListGIDVTK(i,j,k,1)

			end do
		end do
	end do

	write(88,'(a,a,a,I4)') 'SCALARS ', 'Vacancies', ' float', 1
	write(88,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'

	do k=1,numz
		do j=1,numy
			do i=1,numx

				write(88,*) DefectListGIDVTK(i,j,k,2)

			end do
		end do
	end do

	write(88,'(a,a,a,I4)') 'SCALARS ', 'SIAs', ' float', 1
	write(88,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'

	do k=1,numz
		do j=1,numy
			do i=1,numx

				write(88,*) DefectListGIDVTK(i,j,k,3)

			end do
		end do
	end do

	write(88,'(a,a,a,I4)') 'SCALARS ', 'GrainID', ' float', 1
	write(88,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'

	do k=1,numz
		do j=1,numy
			do i=1,numx

				write(88,*) DefectListGIDVTK(i,j,k,4)
	
			end do
		end do
	end do

	close(88)

else

	!Send defect information to master
	call MPI_SEND(numCells, 1, MPI_INTEGER, MASTER, 567, MPI_COMM_WORLD, ierr)

	allocate(defectSend(7,numCells))

	do i=1,numCells
	
		!Step 5.1: identify the x,y,z index of this cell
		!cellIndex(1)=(myMesh(i)%coordinates(1)-myProc%globalCoord(1))/myMesh(1)%length+1
		!cellIndex(2)=(myMesh(i)%coordinates(2)-myProc%globalCoord(3))/myMesh(1)%length+1
		!cellIndex(3)=(myMesh(i)%coordinates(3)-myProc%globalCoord(5))/myMesh(1)%length+1

		defectSend(1,i)=((myMesh(i)%coordinates(1)-myMesh(1)%length/2)-myProc%globalCoord(1))/myMesh(1)%length+1
		defectSend(2,i)=((myMesh(i)%coordinates(2)-myMesh(1)%length/2)-myProc%globalCoord(3))/myMesh(1)%length+1
		defectSend(3,i)=((myMesh(i)%coordinates(3)-myMesh(1)%length/2)-myProc%globalCoord(5))/myMesh(1)%length+1
	
		!Step 5.2: count number of vacancies, interstitials, and He atoms in this cell
		defectCurrent=>defectList(i)
	
		!numV=0
		!numSIA=0
		!numHe=0
		!GrainID=myMesh(i)%material
		defectSend(4,i) = 0	!numCu
		defectSend(5,i) = 0	!numV
		defectSend(6,i) = 0	!numSIA
		defectSend(7,i) = myMesh(i)%material
	
		do while(associated(defectCurrent))
			!numV=numV+defectCurrent%defectType(2)*defectCurrent%num
			!numHe=numHe+defectCurrent%defectType(1)*defectCurrent%num
			!numSIA=numSIA+defectCurrent%defectType(3)*defectCurrent%num + &
			!	defectCurrent%defectType(4)*defectCurrent%num
			defectSend(4,i) = defectSend(4,i) + defectCurrent%defectType(1)*defectCurrent%num
			defectSend(5,i) = defectSend(5,i) + defectCurrent%defectType(2)*defectCurrent%num
			defectSend(6,i) = defectSend(6,i) + defectCurrent%defectType(3)*defectCurrent%num + &
						defectCurrent%defectType(4)*defectCurrent%num
		
			defectCurrent=>defectCurrent%next
		end do
	
		!Step 5.3: Send information to master
		!call MPI_SEND(cellIndex, 3, MPI_INTEGER, MASTER, 571, MPI_COMM_WORLD, ierr)
		!call MPI_SEND(numHe, 1, MPI_INTEGER, MASTER, 568, MPI_COMM_WORLD, ierr)
		!call MPI_SEND(numV, 1, MPI_INTEGER, MASTER, 569, MPI_COMM_WORLD, ierr)
		!call MPI_SEND(numSIA, 1, MPI_INTEGER, MASTER, 570, MPI_COMM_WORLD, ierr)
		!call MPI_SEND(GrainID, 1, MPI_INTEGER, MASTER, 571, MPI_COMM_WORLD, ierr)
	
	end do
	call MPI_SEND(defectSend, 7*numCells, MPI_INTEGER, MASTER, 571, MPI_COMM_WORLD, ierr)

	deallocate(defectSend)

end if

nullify(defectCurrent)

end subroutine

!***************************************************************************************************
!> Subroutine outputRates(step) - outputs reaction rates to a file
!
! Outputs the total reaction rate in each processor to a file. Used for debugging and for parallel
! analysis
!
! Inputs: integer step (the reaction step number)
! Outputs: reaction rate of each processor written in file
!****************************************************************************************************

subroutine outputRates(elapsedTime, step)
use mod_srscd_constants
use DerivedType
implicit none

include 'mpif.h'

integer i, step, status(MPI_STATUS_SIZE)
double precision rate(myProc%numtasks), rateTemp, elapsedTime

if(myProc%taskid==MASTER) then

	!Record master reaction rate
	rate(1)=totalRate
	
	do i=1,myProc%numtasks-1
		
		!Recieve data from other procs
		!record data from other procs in rate()
		
		call MPI_RECV(rateTemp,1,MPI_DOUBLE_PRECISION,i,step,MPI_COMM_WORLD,status,ierr)
		rate(i+1)=rateTemp
		
	end do

	!Output data from other procs to file
	
	write(85,*) 'step', step, 'rates', (rate(i), i=1,myProc%numtasks)
	
else

	!send reaction rate to master proc
	call MPI_SEND(totalRate, 1, MPI_DOUBLE_PRECISION, MASTER, step, MPI_COMM_WORLD, ierr)
	
endif

end subroutine

!****************************************************************************************************
!
!>Subroutine outputDebugRestart(outputCounter) - outputs debug restart.in file
!!
!!Outputs a file restart.in that allows the user to restart the simulation from the
!!point that it ended. This can be used to find errors that occur after significant
!!computation time as well as to run simulations with defects already present.
!
!****************************************************************************************************

subroutine outputDebugRestart(fileNumber, elapsedTime)	!fileNumber = outputCounter
use mod_srscd_constants
use DerivedType
implicit none
include 'mpif.h'

integer status(MPI_STATUS_SIZE)
integer fileNumber
integer cell
integer proc
integer i, j,k
integer numDefectTypes

integer numCellsNeighbor,defectData(numSpecies+1)
double precision coordsAndNumTypesSend(4), coordsAndNumTypesRecv(4)
double precision, allocatable :: defectDataSend(:,:), defectDataRecv(:,:)


double precision elapsedTime
type(defect), pointer :: defectCurrent
character(12) :: fileName

call MPI_ALLREDUCE(numImplantEvents,totalImplantEvents, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
!call MPI_ALLREDUCE(numHeImplantEvents,numHeImplantTotal,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

if(myProc%taskid==MASTER) then
	
	fileName(1:7)='Restart'
	write(unit=fileName(8:9), fmt='(I2)') fileNumber
	fileName(10:12)='.in'
	
	open(88, file=fileName, status='Unknown')
	write(88,*) 'numProcs'
	write(88,*) myProc%numTasks
	write(88,*)
	write(88,*) 'numImplantEvents'
	write(88,*) totalImplantEvents
	write(88,*)
	!write(88,*) 'numHeImplantEvents'
	!write(88,*) numHeImplantTotal
	!write(88,*)
	write(88,*) 'elapsedTime'
	write(88,*) elapsedTime
	write(88,*)
	write(88,*) 'start'
	write(88,*)
	
	!write local defect data
	write(88,*) 'processor'
	write(88,*) myProc%taskid
	write(88,*) 
	
	do cell=1,numCells
	
		write(88,*) 'coordinates'
		write(88,*) (myMesh(cell)%coordinates(i),i=1,3)
		write(88,*)
		
		!Count how many defect types are in this cell
		defectCurrent=>defectList(cell)
		numDefectTypes=0
		
		do while(associated(defectCurrent))
			if(defectCurrent%num > 0) then
				numDefectTypes=numDefectTypes+1
			endif
			defectCurrent=>defectCurrent%next
		end do
		
		!write the defect list into the restart.in file
		write(88,*) 'numDefectTypes'
		write(88,*) numDefectTypes
		
		defectCurrent=>defectList(cell)
		do while(associated(defectCurrent))
			if(defectCurrent%num > 0) then
				write(88,*) (defectCurrent%defectType(i),i=1,numSpecies)
				write(88,*) defectCurrent%num
			endif
			defectCurrent=>defectCurrent%next
		end do
		
		write(88,*) 'end'
		write(88,*)
	
	end do
	
	!Recieve defect data from neighbors
	do proc=1,myProc%numTasks-1
		
		write(88,*) 'processor'
		write(88,*) proc
		write(88,*) 
		
		call MPI_RECV(numCellsNeighbor, 1, MPI_INTEGER, proc, 567, MPI_COMM_WORLD, status, ierr)
		
		do cell=1,numCellsNeighbor
		
			!call MPI_RECV(coordsNeighbor,3,MPI_DOUBLE_PRECISION,proc,2000+cell,MPI_COMM_WORLD,status,ierr)
			call MPI_RECV(coordsAndNumTypesRecv,4,MPI_DOUBLE_PRECISION,proc,1000+cell,MPI_COMM_WORLD,status,ierr)
			write(88,*) 'coordinates'
			write(88,*) (coordsAndNumTypesRecv(i),i=1,3)
			write(88,*)

			numDefectTypes = coordsAndNumTypesRecv(4)
			write(88,*) 'numDefectTypes'
			write(88,*) numDefectTypes

			allocate(defectDataRecv(numSpecies+1,numDefectTypes))
			call MPI_RECV(defectDataRecv,(numSpecies+1)*numDefectTypes,MPI_DOUBLE_PRECISION,proc,3000*cell,MPI_COMM_WORLD,status,ierr)
			
			do i=1,numDefectTypes
				!call MPI_RECV(defectData,numSpecies+1,MPI_INTEGER,proc,3000*cell+i,MPI_COMM_WORLD,status,ierr)
				do k=1, numSpecies
					defectData(k) = defectDataRecv(k,i)
					write(88,*) defectData(k)
				end do
				defectData(numSpecies+1) = defectDataRecv(numSpecies+1,i)
				write(88,*) defectData(numSpecies+1)
			end do
			
			write(88,*) 'end'
			write(88,*)

			deallocate(defectDataRecv)
			
		end do
	
	end do
	
	close(88)

else

	!Send data to master
	
	!Send number of cells to master
	call MPI_SEND(numCells,1,MPI_INTEGER,MASTER,567,MPI_COMM_WORLD,ierr)
	
	do cell=1,numCells
	
		!send coordinates of cell to master
		coordsAndNumTypesSend(1) = myMesh(cell)%coordinates(1)
		coordsAndNumTypesSend(2) = myMesh(cell)%coordinates(2)
		coordsAndNumTypesSend(3) = myMesh(cell)%coordinates(3)

		!call MPI_SEND(myMesh(cell)%coordinates,3,MPI_DOUBLE_PRECISION,MASTER,2000+cell,MPI_COMM_WORLD,ierr)
	
		!Count how many defect types are in the cell
		defectCurrent=>defectList(cell)
		!numDefectTypes=0
		coordsAndNumTypesSend(4)=0
		
		do while(associated(defectCurrent))
			if(defectCurrent%num > 0) then
				!numDefectTypes=numDefectTypes+1
				coordsAndNumTypesSend(4)=coordsAndNumTypesSend(4)+1
			endif
			defectCurrent=>defectCurrent%next
		end do

		!send number of defect types to master
		call MPI_SEND(coordsAndNumTypesSend,4,MPI_DOUBLE_PRECISION,MASTER,1000+cell,MPI_COMM_WORLD,ierr)

		numDefectTypes = coordsAndNumTypesSend(4)
		allocate(defectDataSend(numSpecies+1,numDefectTypes))

		i=1
		defectCurrent=>defectList(cell)
		
		do while(associated(defectCurrent))
			if(defectCurrent%num > 0) then
				do j=1,numSpecies
					defectDataSend(j,i)=defectCurrent%defectType(j)
				end do
				defectDataSend(numSpecies+1,i)=defectCurrent%num
				
				!call MPI_SEND(defectData,numSpecies+1,MPI_INTEGER,MASTER,3000*cell+i,MPI_COMM_WORLD,ierr)
				i=i+1
			endif
			defectCurrent=>defectCurrent%next
		end do
		call MPI_SEND(defectDataSend,(numSpecies+1)*numDefectTypes,MPI_DOUBLE_PRECISION,MASTER,3000*cell,MPI_COMM_WORLD,ierr)

		deallocate(defectDataSend)
	
	end do

endif

end subroutine

!*****************************************************
!function: computeVConc
!*****************************************************

double precision function computeVConc()
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'

integer buffer, tag, i, j, k, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount, same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
!double precision defectsRecv(numSpecies+1), defectsSend(numSpecies+1)

double precision, allocatable :: defectsRecv(:,:)
double precision, allocatable :: defectsSend(:,:)

!Simulation variables passed by main program
double precision elapsedTime
integer step

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!total number of annihilation reactions
integer numAnnihilateTotal

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
end interface

tag=0

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

call MPI_ALLREDUCE(numAnnihilate,numAnnihilateTotal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myProc%taskid==MASTER) then
	
	!Create list of defects in processor 0
	do i=1,numCells
		
		!ONLY OUTPUT DEFECTS IN THE BULK (NOT GBS)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
		
			!Do nothing, no output of defects in grain boundaries
			
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	
	do i=1,myProc%numtasks-1
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,1399,MPI_COMM_WORLD,status,ierr)
		
!		write(*,*) 'recieving ', numDefectsRecv, 'defects'
		allocate(defectsRecv(numSpecies+1,numDefectsRecv))
		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,1400,MPI_COMM_WORLD,status,ierr)

		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,1400,MPI_COMM_WORLD,status,ierr)
!			write(*,*) 'recieved defect ', j
			
			do k=1,numSpecies
				products(k)=defectsRecv(k,j)
			end do
			
			nullify(defectPrevList)
			
			defectCurrentList=>outputDefectList
			
			call findDefectInList(defectCurrentList, defectPrevList, products)
			
			!Next update defects
			if(associated(defectCurrentList)) then !if we aren't at the end of the list
				same=0
				
				do k=1,numSpecies
					if(defectCurrentList%defectType(k)==products(k)) then
						same=same+1
					end if
				end do
				
				if(same==numSpecies) then	
				
					!if the defect is already present in the list
				
					defectCurrentList%num=defectCurrentList%num+defectsRecv(numSpecies+1,j)
				
				else		
					
					!if the defect is to be inserted in the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
					defectPrevList%num=defectsRecv(numSpecies+1,j)
					
					do k=1,numSpecies
						defectPrevList%defectType(k)=products(k)
					end do
					
					!if inserted defect is in the middle of the list, point it to the next item in the list
					
					defectPrevList%next=>defectCurrentList
				end if
			else if(associated(defectPrevList)) then			
				
				!add a defect to the end of the list
				
				nullify(defectPrevList%next)
				allocate(defectPrevList%next)
				nullify(defectPrevList%next%next)
				defectPrevList=>defectPrevList%next
				allocate(defectPrevList%defectType(numSpecies))
				defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
				defectPrevList%num=defectsRecv(numSpecies+1,j)
				
				do k=1,numSpecies
					defectPrevList%defectType(k)=products(k)
				end do
				
			else
				
				write(*,*) 'error tried to insert defect at beginning of output defect list'

			end if
			
			!Signal to other processor to send the next defect
			!call MPI_SEND(tag,1,MPI_INTEGER, i, 1405,MPI_COMM_WORLD, ierr)
		end do

		deallocate(defectsRecv)
		
	end do
	
	defectCurrentList=>outputDefectList
	
	!Initialize Defect counters
	CuNum=0
	VNum=0
	SIANum=0
	totalCu=0
	totalVac=0
	totalSIA=0
	totalVoid=0
	totalLoop=0
	
	do i=1,3
		LoopSize(i)=0
	end do
	VoidNum=0
	
	do while(associated(defectCurrentList))
		
		!Compile statistics for vacancy and SIA concentrations
		if(defectCurrentList%defectType(1) /= 0) then
			
			CuNum=CuNum+defectCurrentList%num
			
			totalCu=totalCu+defectCurrentList%defectType(1)*defectCurrentList%num
		
		end if
		
		if(defectCurrentList%defectType(2) /= 0) then
		
			VNum=VNum+defectCurrentList%num
			
			totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
			
			if(defectCurrentList%defectType(2) >= 45) then
				VoidNum=VoidNum+defectCurrentList%num
				totalVoid = totalVoid + defectCurrentList%defectType(2)*defectCurrentList%num
			end if
		end if
		
		if(defectCurrentList%defectType(4) /= 0 .OR. defectCurrentList%defectType(3) /= 0) then
		
			SIANum=SIANum+defectCurrentList%num
			
			totalSIA=totalSIA+defectCurrentList%defectType(4)*defectCurrentList%num
			totalSIA=totalSIA+defectCurrentList%defectType(3)*defectCurrentList%num
		
			if(defectCurrentList%defectType(4) >= 16) then
				LoopSize(1)=LoopSize(1)+defectCurrentList%num
			end if
			
			if(defectCurrentList%defectType(4) >= 20) then
				LoopSize(2)=LoopSize(2)+defectCurrentList%num
				totalLoop = totalLoop + defectCurrentList%defectType(4)*defectCurrentList%num
			end if
			
			if(defectCurrentList%defectType(4) >= 24) then
				LoopSize(3)=LoopSize(3)+defectCurrentList%num
			end if
		end if
	
		defectCurrentList=>defectCurrentList%next
	end do
	
	computeVConc=dble(Vnum)*atomSize/systemVol
	!do i=1,myProc%numtasks-1
	!	call MPI_SEND(dble(Vnum)*atomsize/systemVol,1,MPI_DOUBLE_PRECISION, i, 1406,MPI_COMM_WORLD, ierr)
	!end do

else
	numDefectsSend=0
	
	!Create list of defects in processor 
	do i=1,numCells
	
		!ONLY OUTPUT DEFECTS IN BULK (NOT GBs)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
			! Do nothing, no output of defects in grain boundaries
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
						
						numDefectsSend=numDefectsSend+1
					end if
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
					numDefectsSend=numDefectsSend+1
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				end if
				
				defectCurrent=>defectCurrent%next
			end do
		
		end if
		
	end do
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 1399, MPI_COMM_WORLD, ierr)

	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		
		!call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 1400, MPI_COMM_WORLD, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,1405,MPI_COMM_WORLD,status,ierr)
		
		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		end if
		
!		write(*,*) 'sent defect', i
	end do

	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 1400, MPI_COMM_WORLD, ierr)

	deallocate(defectsSend)
	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)

	!call MPI_RECV(computeVConc,1,MPI_DOUBLE_PRECISION,MASTER,1406,MPI_COMM_WORLD,status,ierr)
end if

call MPI_BCAST(computeVConc, 1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,ierr)

!Deallocate memory

defectCurrentList=>outputDefectList
nullify(defectPrevList)

do while(Associated(defectCurrentList))
	defectPrevList=>defectCurrentList
	defectCurrentList=>defectCurrentList%next
	
	if(allocated(defectPrevList%defectType)) deallocate(defectPrevList%defectType)
	
	deallocate(defectPrevList)
end do

nullify(defectCurrentList)
nullify(outputDefectList)

end function

!******************************************************
!function: computeIConc()
!******************************************************

double precision function computeIConc()
use DerivedType
use mod_srscd_constants
implicit none

include 'mpif.h'

integer buffer, tag, i, j, k, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount, same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
!double precision defectsRecv(numSpecies+1), defectsSend(numSpecies+1)

double precision, allocatable :: defectsRecv(:,:)
double precision, allocatable :: defectsSend(:,:)

!Simulation variables passed by main program
double precision elapsedTime
integer step

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!total number of annihilation reactions
integer numAnnihilateTotal

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
	use mod_srscd_constants
	type(defect), pointer :: defectCurrent, defectPrev
	integer products(numSpecies)
	end subroutine
end interface

tag=0

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

call MPI_ALLREDUCE(numAnnihilate,numAnnihilateTotal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myProc%taskid==MASTER) then
	
	!Create list of defects in processor 0
	do i=1,numCells
		
		!ONLY OUTPUT DEFECTS IN THE BULK (NOT GBS)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
		
			!Do nothing, no output of defects in grain boundaries
			
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	
	do i=1,myProc%numtasks-1
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,2399,MPI_COMM_WORLD,status,ierr)
		
!		write(*,*) 'recieving ', numDefectsRecv, 'defects'
		allocate(defectsRecv(numSpecies+1,numDefectsRecv))

		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,2400,MPI_COMM_WORLD,status,ierr)
		
		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,2400,MPI_COMM_WORLD,status,ierr)
!			write(*,*) 'recieved defect ', j
			
			do k=1,numSpecies
				products(k)=defectsRecv(k,j)
			end do
			
			nullify(defectPrevList)
			
			defectCurrentList=>outputDefectList
			
			call findDefectInList(defectCurrentList, defectPrevList, products)
			
			!Next update defects
			if(associated(defectCurrentList)) then !if we aren't at the end of the list
				same=0
				
				do k=1,numSpecies
					if(defectCurrentList%defectType(k)==products(k)) then
						same=same+1
					endif
				end do
				
				if(same==numSpecies) then	
				
					!if the defect is already present in the list
				
					defectCurrentList%num=defectCurrentList%num+defectsRecv(numSpecies+1,j)
				
				else		
					
					!if the defect is to be inserted in the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
					defectPrevList%num=defectsRecv(numSpecies+1,j)
					
					do k=1,numSpecies
						defectPrevList%defectType(k)=products(k)
					end do
					
					!if inserted defect is in the middle of the list, point it to the next item in the list
					
					defectPrevList%next=>defectCurrentList
				endif
			else if(associated(defectPrevList)) then			
				
				!add a defect to the end of the list
				
				nullify(defectPrevList%next)
				allocate(defectPrevList%next)
				nullify(defectPrevList%next%next)
				defectPrevList=>defectPrevList%next
				allocate(defectPrevList%defectType(numSpecies))
				defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
				defectPrevList%num=defectsRecv(numSpecies+1,j)
				
				do k=1,numSpecies
					defectPrevList%defectType(k)=products(k)
				end do
				
			else
				
				write(*,*) 'error tried to insert defect at beginning of output defect list'

			endif
			
			!Signal to other processor to send the next defect
			!call MPI_SEND(tag,1,MPI_INTEGER, i, 2405,MPI_COMM_WORLD, ierr)
		end do
		
	end do
	
	defectCurrentList=>outputDefectList
	
	!Initialize Defect counters
	CuNum=0
	VNum=0
	SIANum=0
	totalCu=0
	totalVac=0
	totalSIA=0
	totalVoid=0
	totalLoop=0
	
	do i=1,3
		LoopSize(i)=0
	end do
	VoidNum=0
	
	do while(associated(defectCurrentList))
		
		!Compile statistics for vacancy and SIA concentrations
		if(defectCurrentList%defectType(1) /= 0) then
			
			CuNum=CuNum+defectCurrentList%num
			
			totalCu=totalCu+defectCurrentList%defectType(1)*defectCurrentList%num
		
		endif
		
		if(defectCurrentList%defectType(2) /= 0) then
		
			VNum=VNum+defectCurrentList%num
			
			totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
			
			if(defectCurrentList%defectType(2) >= 45) then
				VoidNum=VoidNum+defectCurrentList%num
				totalVoid = totalVoid + defectCurrentList%defectType(2)*defectCurrentList%num
			endif
		endif
		
		if(defectCurrentList%defectType(4) /= 0 .OR. defectCurrentList%defectType(3) /= 0) then
		
			SIANum=SIANum+defectCurrentList%num
			
			totalSIA=totalSIA+defectCurrentList%defectType(4)*defectCurrentList%num
			totalSIA=totalSIA+defectCurrentList%defectType(3)*defectCurrentList%num
		
			if(defectCurrentList%defectType(4) >= 16) then
				LoopSize(1)=LoopSize(1)+defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 20) then
				LoopSize(2)=LoopSize(2)+defectCurrentList%num
				totalLoop = totalLoop + defectCurrentList%defectType(4)*defectCurrentList%num
			endif
			
			if(defectCurrentList%defectType(4) >= 24) then
				LoopSize(3)=LoopSize(3)+defectCurrentList%num
			endif
		endif
	
		defectCurrentList=>defectCurrentList%next
	end do
	
	computeIConc=dble(SIAnum)*atomSize/systemVol
	!do i=1,myProc%numtasks-1
		!call MPI_SEND(dble(SIAnum)*atomsize/systemVol,1,MPI_DOUBLE_PRECISION, i, 2406,MPI_COMM_WORLD, ierr)
	!end do

else
	numDefectsSend=0
	
	!Create list of defects in processor 
	do i=1,numCells
	
		!ONLY OUTPUT DEFECTS IN BULK (NOT GBs)
		if(numMaterials > 1 .AND. myMesh(i)%material /= 1) then
			! Do nothing, no output of defects in grain boundaries
		else
	
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				
				nullify(defectPrevList)
				
				defectCurrentList=>outputDefectList
				
				call findDefectInList(defectCurrentList, defectPrevList, defectCurrent%defectType)
				
				!Next update defects
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					
					do j=1,numSpecies
						if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
							same=same+1
						endif
					end do
					
					if(same==numSpecies) then	
					
						!if the defect is already present in the list
					
						defectCurrentList%num=defectCurrentList%num+defectCurrent%num
					
					else		
						
						!if the defect is to be inserted in the list
						
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(numSpecies))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num
						
						do j=1,numSpecies
							defectPrevList%defectType(j)=defectCurrent%defectType(j)
						end do
						
						!if inserted defect is in the middle of the list, point it to the next item in the list
						
						defectPrevList%next=>defectCurrentList
						
						numDefectsSend=numDefectsSend+1
					endif
				else if(associated(defectPrevList)) then			
					
					!add a defect to the end of the list
					
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(numSpecies))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num
					
					do j=1,numSpecies
						defectPrevList%defectType(j)=defectCurrent%defectType(j)
					end do
					
					numDefectsSend=numDefectsSend+1
					
				else
					
					write(*,*) 'error tried to insert defect at beginning of output defect list'
	
				endif
				
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 2399, MPI_COMM_WORLD, ierr)

	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		
		!call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 2400, MPI_COMM_WORLD, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,2405,MPI_COMM_WORLD,status,ierr)
		
		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		endif
		
!		write(*,*) 'sent defect', i
	end do

	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 2400, MPI_COMM_WORLD, ierr)
	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
	
	!call MPI_RECV(computeIConc,1,MPI_DOUBLE_PRECISION,MASTER,2406,MPI_COMM_WORLD,status,ierr)
endif

call MPI_BCAST(computeIConc, 1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,ierr)
!Deallocate memory

defectCurrentList=>outputDefectList
nullify(defectPrevList)

do while(Associated(defectCurrentList))
	defectPrevList=>defectCurrentList
	defectCurrentList=>defectCurrentList%next
	
	if(allocated(defectPrevList%defectType)) deallocate(defectPrevList%defectType)
	
	deallocate(defectPrevList)
end do

nullify(defectCurrentList)
nullify(outputDefectList)

end function
