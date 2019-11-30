!***************************************************************************************************
!>subroutine outputDefects - outputs raw data for defect populations in various volume elements
!
! Outputs into file: rawdat.out. These contain the complete defect populations per volume element.
! Compiles data from local as well as global processors.
!***************************************************************************************************

subroutine outputDefectsXYZ()
use mod_constants
use DerivedType
implicit none
include 'mpif.h'

integer i, j, k,l, status(MPI_STATUS_SIZE), defectCount, numRecv
integer numScluster, numVoid, numLoop	!total number of solute/V/SIA/mSnV clusters in per volume element
double precision volTemp

type(defect), pointer :: defectCurrent
integer, allocatable :: cellSend(:,:)
integer, allocatable :: cellRecv(:,:)
integer, allocatable :: numCellRecv(:)

volTemp=atomSize/(meshLength**(3d0))

allocate(numCellRecv(myProc%numtasks))
call MPI_GATHER(numCells,1,MPI_INTEGER,numCellRecv,1,MPI_INTEGER,MASTER,comm,ierr)

if(myProc%taskid/=MASTER) then

	do i=1, numCells

		defectCount=0
		defectCurrent=>defectList(i)%next
		do while(associated(defectCurrent))
			defectCount=defectCount+1
			defectCurrent=>defectCurrent%next
		end do

		allocate(cellSend(numSpecies+1,defectCount+1))

		cellSend(1,1)=defectCount
		cellSend(1,2)=myMesh(i)%globalCell
		cellSend(1,3)=0		!number of Cu clusters
		cellSend(1,4)=0		!number of Cu clusters
		cellSend(1,5)=0		!number of loops

		numScluster=0		!number of Cu clusters
		numVoid=0			!number of void
		numLoop=0			!number of loops

		j=1
		defectCurrent=>defectList(i)%next
		do while(associated(defectCurrent))

			j=j+1
			do k=1, numSpecies
				cellSend(k,j) = defectCurrent%defectType(k)
			end do
			cellSend(numSpecies+1,j) = defectCurrent%num

			if(defectCurrent%defectType(1) > minSCluster) then
				numScluster=numScluster+defectCurrent%num
			else if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2) > minVoid) then
				numVoid=numVoid+defectCurrent%num
			else if(defectCurrent%defectType(3)>minLoop .OR. defectCurrent%defectType(4)>minLoop) then
				numLoop=numLoop+defectCurrent%num
			end if
			defectCurrent=>defectCurrent%next
		end do

		cellSend(1,3)=numScluster
		cellSend(1,4)=numVoid
		cellSend(1,5)=numLoop

		call MPI_SEND(cellSend,(numSpecies+1)*(defectCount+1),MPI_INTEGER,MASTER,100+i,comm,ierr)
		deallocate(cellSend)
	end do

else

    write(82,*) 'time', elapsedTime, 'DPA', DPA, 'steps', step
	
	!Write defects in processor 0 to file
	write(82,*) 'processor', myProc%taskid
	do i=1,numCells
		defectCurrent=>defectList(i)%next
		write(82,*) 'cell', i,'globalCell', myMesh(i)%globalCell
        write(82,*) 'defects (Cu V SIA_m SIA_im num)'

		numScluster=0
		numVoid=0
		numLoop=0

		do while(associated(defectCurrent))
			write(82,*) defectCurrent%defectType, defectCurrent%num
			if(defectCurrent%defectType(1) > minSCluster) then
				numScluster=numScluster+defectCurrent%num
			else if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2) > minVoid) then
				numVoid=numVoid+defectCurrent%num
			else if(defectCurrent%defectType(3)>minLoop .OR. defectCurrent%defectType(4)>minLoop) then
				numLoop=numLoop+defectCurrent%num
			end if
			defectCurrent=>defectCurrent%next
		end do

		write(82,*) 'concentration of point defects (Cu/V/SIA) and clusters (Cu cluster/Void/Loop)'
		write(82,*) dble(numScluster)*volTemp, dble(numVoid)*volTemp, dble(numLoop)*volTemp
		write(82,*)
	end do
	write(82,*)
	write(82,*)
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	do i=1,myProc%numtasks-1
		write(82,*) 'processor', i
		
		do j=1,numCellRecv(i+1)

			numRecv=0
			call MPI_PROBE(i, 100+j,comm,status,ierr)
			call MPI_GET_COUNT(status,MPI_INTEGER,numRecv,ierr)

			numRecv=numRecv/(numSpecies+1)
            allocate(cellRecv(numSpecies+1,numRecv))
            call MPI_RECV(cellRecv,(numSpecies+1)*numRecv,MPI_INTEGER,i,100+j,comm,status,ierr)

			write(82,*) 'cell', j,'globalCell', cellRecv(2,1)
			write(82,*) 'defects (Cu V SIA_m SIA_im num)'

			do k=1,cellRecv(1,1)
				write(82,*) cellRecv(:,k+1)
			end do

			write(82,*) 'concentration of point defects (Cu/V/SIA) and clusters (Cu cluster/Void/Loop)'
			write(82,*) dble(cellRecv(3,1))*volTemp,dble(cellRecv(4,1))*volTemp,dble(cellRecv(5,1))*volTemp
			write(82,*)

			deallocate(cellRecv)
		end do
	end do
    write(82,*)
	write(82,*)

end if

end subroutine

!*******************************************************************************************************
!> Subroutine outputDefectsTotal() - outputs the total defects and post processing for the entire volume
!!
!! Outputs a list of all defects in system (defect type, num), regardless of volume element.
!! Compiles such a list using message passing between processors.
!! Used to find total defect populations in simulation volume, not spatial distribution of defects.
!!
!! This subroutine will also automatically output the concentration of solute, vacancies, and SIAs of size large.
!*******************************************************************************************************

subroutine outputDefectsTotal()
use DerivedType
use mod_constants
implicit none
include 'mpif.h'

integer i, j, k, status(MPI_STATUS_SIZE)
integer products(numSpecies), same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
type(cascade), pointer :: cascadeCurrent
integer, allocatable :: defectsRecv(:,:)
integer, allocatable :: defectsSend(:,:)
integer numDefectsSend
integer, allocatable :: numDefectsRecv(:)

!number of point defects (solute atom, vacancy, self-interstitial atom)
integer pointS, pointV, pointSIA
!total retained vacancies/self-interstitial atoms in the whole system
integer totalVac, totalSIA
!total number of solute atoms/vacancies/SIAs in clusters
integer numS, numV, numSIA, numSV
!total number of solute/V/SIA/mSnV clusters in the whole system
integer numScluster, numVoid, numLoop, numSVcluster
!Number density of clusters
double precision denScluster, denVoid, denLoop, denSVcluster
!Average radius of clusters
double precision radiusScluster, radiusVoid, radiusLoop, radiusSVcluster
!Average size of clusters
double precision sizeScluster, sizeVoid, sizeLoop, sizeSVcluster

double precision VRetained, VAnnihilated, conPointV, conPointSIA

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

allocate(numDefectsRecv(myProc%numtasks))

!Count defects in coarse mesh
numDefectsSend=0
do i=1,numCells

	!Only output defects in bulk (not GBs)
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

!Count defects in fine mesh
if(meshingType=='adaptive') then
	cascadeCurrent=>ActiveCascades
	do while(associated(cascadeCurrent))
		do i=1,numCellsCascade
			defectCurrent=>cascadeCurrent%localDefects(i)%next
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
				endif

				defectCurrent=>defectCurrent%next
			end do
		end do
		cascadeCurrent=>cascadeCurrent%next
	end do
end if

call MPI_GATHER(numDefectsSend,1,MPI_INTEGER,numDefectsRecv,1,MPI_INTEGER,MASTER,comm,ierr)

!Update send buffer
if(myProc%taskid /= MASTER) then

	defectCurrentList=>outputDefectList%next
	allocate(defectsSend(numSpecies+1,numDefectsSend))
	i=0
	do while(associated(defectCurrentList))
		i=i+1
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		defectCurrentList=>defectCurrentList%next
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		end if
	end do
	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_INTEGER, MASTER, 400, comm, ierr)
	deallocate(defectsSend)
else	!MASTER

	!Recv defects from other processors
	do i=1, myProc%numtasks-1
		allocate(defectsRecv(numSpecies+1,numDefectsRecv(i+1)))
		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv(i+1),MPI_INTEGER,i,400,comm,status,ierr)

		do j=1,numDefectsRecv(i+1)

			do k=1,numSpecies
				products(k)=defectsRecv(k,j)
			end do

			nullify(defectPrevList)
			defectCurrentList=>outputDefectList
			call findDefectInList(defectCurrentList, defectPrevList, products)

			!Update outputDefectList
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
		deallocate(defectsRecv)
	end do
end if

!Output totdat.out
if(myProc%taskid==MASTER) then

	!number of point defects
	pointS=0
	pointV=0
	pointSIA=0

	totalVac=0
	totalSIA=0

	!total number of solute (Cu) atoms/vacancies/self-interstitial atoms in clusters
	numS=0
    numV=0
	numSIA=0
	numSV=0

	!total number of solute (Cu)/V/SIA clusters in the whole system
	numScluster=0
	numVoid=0
	numLoop=0
	numSVcluster=0

	!number density of clusters
	denScluster=0d0
	denVoid=0d0
	denLoop=0d0
	denSVcluster=0d0

	!average radius of clusters
	radiusScluster=0d0
	radiusVoid=0d0
	radiusLoop=0d0
	radiusSVcluster=0d0

	!average size of clusters
	sizeScluster=0d0
	sizeVoid=0d0
	sizeLoop=0d0
	sizeSVcluster=0d0

	VRetained=0d0
	VAnnihilated=0d0

	write(83,*) 'defects (Cu, V, SIA_m, SIA_im, num)'
	defectCurrentList=>outputDefectList%next
	do while(associated(defectCurrentList))

		if(defectCurrentList%defectType(1) /= 0) then	!Cu/CuV cluster

			if(defectCurrentList%defectType(1)==1 .AND. defectCurrentList%defectType(2)==0) then
				pointS=defectCurrentList%num
			end if

			if(defectCurrentList%defectType(1) > minSCluster) then
				numS=numS+defectCurrentList%defectType(1)*defectCurrentList%num
				radiusScluster=radiusScluster+dble(defectCurrentList%num)*&
						(3*dble(defectCurrentList%defectType(1))*atomSize/(4*pi))**(1d0/3d0)
				numScluster=numScluster+defectCurrentList%num
			end if

			if((defectCurrentList%defectType(1)+defectCurrentList%defectType(2)) > minSV) then
				numSV=numSV+max(defectCurrentList%defectType(1), defectCurrentList%defectType(2))*defectCurrentList%num
				radiusSVcluster = radiusSVcluster+dble(defectCurrentList%num)*&
						(3*dble(max(defectCurrentList%defectType(1), defectCurrentList%defectType(2)))*&
								atomSize/(4*pi))**(1d0/3d0)
				numSVcluster=numSVcluster+defectCurrentList%num
				if(defectCurrentList%defectType(2) /= 0) then
					totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
				end if
			end if

		else if(defectCurrentList%defectType(2) /= 0) then	!V cluster

			totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num

			if(defectCurrentList%defectType(2)==1) then
				pointV=defectCurrentList%num
			end if

			if(defectCurrentList%defectType(2) > minVoid) then
				numV=numV+defectCurrentList%defectType(2)*defectCurrentList%num
				radiusVoid=radiusVoid+dble(defectCurrentList%num)*&
						(3d0*dble(defectCurrentList%defectType(2))*atomSize/(4d0*pi))**(1d0/3d0)
				numVoid=numVoid+defectCurrentList%num
			end if

		else if(defectCurrentList%defectType(3) /= 0) then	!Loop

			if(defectCurrentList%defectType(3) ==1) then
				pointSIA=defectCurrentList%num
			end if

			if(defectCurrentList%defectType(3) > minLoop) then
				numSIA=numSIA+defectCurrentList%defectType(3)*defectCurrentList%num
				radiusLoop=radiusLoop+dble(defectCurrentList%num)*&
						(dble(defectCurrentList%defectType(3))*atomSize/(pi*burgers))**(1d0/2d0)
				numLoop=numLoop+defectCurrentList%num
			end if

		else if(defectCurrentList%defectType(4) /= 0) then

			if(defectCurrentList%defectType(4) > minLoop) then
				numSIA=numSIA+defectCurrentList%defectType(4)*defectCurrentList%num
				radiusLoop=radiusLoop+dble(defectCurrentList%num)*&
						(dble(defectCurrentList%defectType(4))*atomSize/(pi*burgers))**(1d0/2d0)
				numLoop=numLoop+defectCurrentList%num
			end if

		end if

		!Output defects
		write(83,*) defectCurrentList%defectType, defectCurrentList%num

		defectCurrentList=>defectCurrentList%next
	end do

	denScluster = dble(numScluster)/systemVol
	denVoid = dble(numVoid)/systemVol
	denLoop = dble(numLoop)/systemVol
	denSVcluster = dble(numSVcluster)/systemVol

	!radiusScluster=(3*(dble(numS)/dble(numScluster))*atomSize/(4*pi))**(1d0/3d0)
	!radiusVoid=(3d0*(dble(numV)/dble(numVoid))*atomSize/(4d0*pi))**(1d0/3d0)
    !radiusLoop=((dble(numSIA)/dble(numLoop))*atomSize/(pi*burgers))**(1d0/2d0)
	!radiusSVcluster=(3d0*(dble(numSV)/dble(numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)
	radiusScluster=radiusScluster/dble(numScluster)
	radiusVoid=radiusVoid/dble(numVoid)
	radiusLoop=radiusLoop/dble(numLoop)
	radiusSVcluster=radiusSVcluster/dble(numSVcluster)

	sizeScluster=dble(numS)/dble(numScluster)
	sizeVoid=dble(numV)/dble(numVoid)
	sizeLoop=dble(numSIA)/dble(numLoop)
	!sizeSVcluster=dble(numSV)/dble(numSVcluster)

	conPointV=dble(pointV)/systemVol*atomSize
	conPointSIA=dble(pointSIA)/systemVol*atomSize

	VRetained = dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
	VAnnihilated = dble(totalImpAnn(2))/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))

	!Output totdat.out
	write(83,*)	'numCluster (S/Void/Loop/SV):', numScluster, numVoid, numLoop, numSVcluster
	write(83,*)	'NumberDensity (m-3) (S/Void/Loop/SV):', denScluster*1d27, denVoid*1d27,denLoop*1d27,denSVcluster*1d27
	write(83,*)	'Concentration (S/Void/Loop/SV):', denScluster*atomSize, denVoid*atomSize, denLoop*atomSize,&
			denSVcluster*atomSize
	write(83,*)	'AverageRadius (nm) (S/Void/Loop):', radiusScluster, radiusVoid, radiusLoop
	write(83,*)	'AverageSize (S/Void/Loop):', sizeScluster, sizeVoid, sizeLoop
	write(83,*) 'ConcenPointDefects (V/SIA):', conPointV, conPointSIA
	write(83,*) 'PercentVRetained',VRetained,'PercentVAnnihilated',VAnnihilated
	write(83,*)
	write(83,*)

end if

!Deallocate outputDefectList
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

deallocate(numDefectsRecv)

end subroutine

!***********************************************************************
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
!***********************************************************************

subroutine outputDefectsBoundary()
use DerivedType
use mod_constants
implicit none

include 'mpif.h'

integer buffer, tag, i, j, k, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
integer defectTypeRecv(numSpecies), cellNumberRecv, numRecv, defectCount, same
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList

double precision, allocatable :: defectsRecv(:,:), defectsSend(:,:)

double precision atomArea, systemArea

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
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
				end if
				defectCurrent=>defectCurrent%next
			end do
		
		endif
		
	end do
	
	!Other processors will send information about defects contained in their mesh. This processor 
	!must output all information to data file (can't do it from other processors). Information should
	!be output in the same format for the master and slave processors.
	
	do i=1,myProc%numtasks-1
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,399,comm,status,ierr)

		allocate(defectsRecv(numSpecies+1,numDefectsRecv))
		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,400,comm,status,ierr)
		
		do j=1,numDefectsRecv
			
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
		deallocate(defectsRecv)
		
	end do
	
	!Output defect list
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

		write(83,*) defectCurrentList%defectType, defectCurrentList%num
		defectCurrentList=>defectCurrentList%next
	end do

	atomArea=(9d0*pi*atomSize**2d0/16d0)**(1d0/3d0)
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
	write(84,*) 'PercentVRetained', dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
	write(84,*) 'PercentVAnnihilated', dble(totalImpAnn(2))/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
	
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
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 399, comm, ierr)
	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num

		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		end if

	end do
	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 400, comm, ierr)
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
use mod_constants
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
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

!Allocate and initialize array of defect numbers in each processor
allocate(defectProfileArray(numSpecies,numZ))
allocate(defectNumArray(numSpecies,numZ))
allocate(defectProfile(numZ))

allocate(tempProfileArray(numSpecies,numZ))
allocate(tempNumArray(numSpecies,numZ))

do i=1,numz
	allocate(defectProfile(i)%defectType(numSpecies))
	do j=1,numSpecies
		defectProfileArray(j,i)=0
		defectNumArray(j,i)=0
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
				defectProfileArray(j,zEntry)=defectProfileArray(j,zEntry)+defectCurrent%defectType(j)*defectCurrent%num
				if(defectCurrent%defectType(j) > 0) then
					defectNumArray(j,zEntry)=defectNumArray(j,zEntry)+defectCurrent%num
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

		call MPI_RECV(tempProfileArray,numz*numSpecies,MPI_INTEGER, procID, 700, comm, status, ierr)
		call MPI_RECV(tempNumArray,numz*numSpecies,MPI_INTEGER, procID, 600,comm, status, ierr)
		
		do i=1,numz

			!Add defects of neighboring procs to master processor defect array
			do j=1,numSpecies
				defectProfileArray(j,i)=defectProfileArray(j,i)+tempProfileArray(j,i)
			end do

			do j=1,numSpecies
				defectNumArray(j,i)=defectNumArray(j,i)+tempNumArray(j,i)
			end do
			
			call MPI_RECV(numDefectsRecv, 1, MPI_INTEGER, procID, i*800, comm, status, ierr)

			allocate(buffer(numSpecies+1,numDefectsRecv))
			call MPI_RECV(buffer, (numSpecies+1)*numDefectsRecv, MPI_INTEGER, procID, i*900, comm, status, ierr)
			
			do j=1,numDefectsRecv

				do k=1,numSpecies
					defectTypeStore(k)=buffer(k,j)
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
						defectListCurrent%num=defectListCurrent%num+buffer(numSpecies+1,j)
					
					else		
						
						!if the defect is to be inserted in the list
						nullify(defectListPrev%next)
						allocate(defectListPrev%next)
						nullify(defectListPrev%next%next)
						defectListPrev=>defectListPrev%next
						allocate(defectListPrev%defectType(numSpecies))
						defectListPrev%cellNumber=i
						defectListPrev%num=buffer(numSpecies+1,j)
						
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
					defectListPrev%num=buffer(numSpecies+1,j)
					
					do k=1,numSpecies
						defectListPrev%defectType(k)=defectTypeStore(k)
					end do
					
				else
					write(*,*) 'error tried to insert defect at beginning of output defect list'
				end if

			end do

			deallocate(buffer)
		
		end do
	
	end do
	
	!Outputs defectProfileArray to a file
	filename(1:12)='DepthProfile'
	write(unit=filename(13:14), fmt='(I2)') sim
	filename(15:18)='.out'
	
	open(99, file=filename, action='write', status='Unknown')

	do i=1,numz
		!XCoord=(i-1)*myMesh(1)%length+myMesh(1)%length/2d0
		ZCoord=(i-1)*myMesh(1)%length+myMesh(1)%length/2d0

		!Vacancy concen
		write(99,*) ZCoord, dble(defectProfileArray(2,i)/(numy*numx*((1d-9*myMesh(1)%length)**3d0))), &
			dble(defectNumArray(2,i)/(numy*numx*((1d-9*myMesh(1)%length)**3d0)))
			
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
	call MPI_SEND(defectProfileArray,numz*numSpecies,MPI_INTEGER, MASTER, 700, comm, ierr)
	call MPI_SEND(defectNumArray, numz*numSpecies, MPI_INTEGER, MASTER, 600, comm, ierr)

	do i=1,numz

		!Compute how many defect types are in this list
		numDefectsSend=0
		defectCurrent=>defectProfile(i)%next
		do while(associated(defectCurrent))
			numDefectsSend=numDefectsSend+1
			defectCurrent=>defectCurrent%next
		end do
		
		!Send how many defect types are going to be sent
		call MPI_SEND(numDefectsSend,1,MPI_INTEGER, MASTER, i*800, comm, ierr)
		
		!Send each defect
		defectCurrent=>defectProfile(i)%next
		numSend=0

		allocate(buffer(numSpecies+1,numDefectsSend))
		
		do while(associated(defectCurrent))
			
			!Create buffer of information to send
			numSend=numSend+1
			do j=1,numSpecies
				buffer(j,numSend)=defectCurrent%defectType(j)
			end do
			buffer(numSpecies+1,numSend)=defectCurrent%num

			defectCurrent=>defectCurrent%next
		end do
		call MPI_SEND(buffer, (numSpecies+1)*numDefectsSend, MPI_INTEGER, MASTER, i*900, comm, ierr)
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
!> Subroutine outputDefectsVTK - outputs defect populations in vtk file format
!
! Outputs defects in a VTK formatted file (3D unstructured grid of linear cubes)
!***********************************************************************

subroutine outputDefectsVTK(fileNumber)	!fileNumber = outputCounter
use mod_constants
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
	
		call MPI_RECV(numCellsNeighbor, 1, MPI_INTEGER, i, 567, comm, status, ierr)
		allocate(defectRecv(7,numCellsNeighbor))
		call MPI_RECV(defectRecv, 7*numCellsNeighbor, MPI_INTEGER, i, 571, comm, status, ierr)

		do j=1,numCellsNeighbor

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
	call MPI_SEND(numCells, 1, MPI_INTEGER, MASTER, 567, comm, ierr)

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
		!call MPI_SEND(cellIndex, 3, MPI_INTEGER, MASTER, 571, comm, ierr)
		!call MPI_SEND(numHe, 1, MPI_INTEGER, MASTER, 568, comm, ierr)
		!call MPI_SEND(numV, 1, MPI_INTEGER, MASTER, 569, comm, ierr)
		!call MPI_SEND(numSIA, 1, MPI_INTEGER, MASTER, 570, comm, ierr)
		!call MPI_SEND(GrainID, 1, MPI_INTEGER, MASTER, 571, comm, ierr)
	
	end do
	call MPI_SEND(defectSend, 7*numCells, MPI_INTEGER, MASTER, 571, comm, ierr)

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

subroutine outputRates()
use mod_constants
use DerivedType
implicit none
include 'mpif.h'

integer i,status(MPI_STATUS_SIZE)
double precision rate(myProc%numtasks), rateTemp

call MPI_GATHER(totalRate,1,MPI_DOUBLE_PRECISION,rate,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

if(myProc%taskid==MASTER) then

	write(85,*) 'step', step, 'rates', (rate(i), i=1,myProc%numtasks)

end if

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

subroutine outputDebugRestart(fileNumber)	!fileNumber = outputCounter
use mod_constants
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

type(defect), pointer :: defectCurrent
character(12) :: fileName

if(myProc%taskid==MASTER) then
	
	fileName(1:7)='restart'
	write(unit=fileName(8:9), fmt='(I2)') fileNumber
	fileName(10:12)='.in'
	
	open(88, file=fileName, status='Unknown')
	write(88,*) 'numProcs'
	write(88,*) myProc%numTasks
	write(88,*)
	write(88,*) 'numImplantEvents'
	write(88,*) totalImpAnn(1)
	write(88,*)
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
		
		call MPI_RECV(numCellsNeighbor, 1, MPI_INTEGER, proc, 567, comm, status, ierr)
		
		do cell=1,numCellsNeighbor
		
			!call MPI_RECV(coordsNeighbor,3,MPI_DOUBLE_PRECISION,proc,2000+cell,comm,status,ierr)
			call MPI_RECV(coordsAndNumTypesRecv,4,MPI_DOUBLE_PRECISION,proc,1000+cell,comm,status,ierr)
			write(88,*) 'coordinates'
			write(88,*) (coordsAndNumTypesRecv(i),i=1,3)
			write(88,*)

			numDefectTypes = coordsAndNumTypesRecv(4)
			write(88,*) 'numDefectTypes'
			write(88,*) numDefectTypes

			allocate(defectDataRecv(numSpecies+1,numDefectTypes))
			call MPI_RECV(defectDataRecv,(numSpecies+1)*numDefectTypes,MPI_DOUBLE_PRECISION,proc,3000*cell,comm,status,ierr)
			
			do i=1,numDefectTypes
				!call MPI_RECV(defectData,numSpecies+1,MPI_INTEGER,proc,3000*cell+i,comm,status,ierr)
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
	call MPI_SEND(numCells,1,MPI_INTEGER,MASTER,567,comm,ierr)
	
	do cell=1,numCells
	
		!send coordinates of cell to master
		coordsAndNumTypesSend(1) = myMesh(cell)%coordinates(1)
		coordsAndNumTypesSend(2) = myMesh(cell)%coordinates(2)
		coordsAndNumTypesSend(3) = myMesh(cell)%coordinates(3)

		!call MPI_SEND(myMesh(cell)%coordinates,3,MPI_DOUBLE_PRECISION,MASTER,2000+cell,comm,ierr)
	
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
		call MPI_SEND(coordsAndNumTypesSend,4,MPI_DOUBLE_PRECISION,MASTER,1000+cell,comm,ierr)

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
				
				!call MPI_SEND(defectData,numSpecies+1,MPI_INTEGER,MASTER,3000*cell+i,comm,ierr)
				i=i+1
			endif
			defectCurrent=>defectCurrent%next
		end do
		call MPI_SEND(defectDataSend,(numSpecies+1)*numDefectTypes,MPI_DOUBLE_PRECISION,MASTER,3000*cell,comm,ierr)

		deallocate(defectDataSend)
	
	end do

endif

end subroutine

!*****************************************************
!function: computeVConc
!*****************************************************

double precision function computeVConc()
use DerivedType
use mod_constants
implicit none
include 'mpif.h'

integer i, j, k, same, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList

double precision, allocatable :: defectsRecv(:,:)
double precision, allocatable :: defectsSend(:,:)

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

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
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,1399,comm,status,ierr)

		allocate(defectsRecv(numSpecies+1,numDefectsRecv))
		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,1400,comm,status,ierr)

		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,1400,comm,status,ierr)
			
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
	!	call MPI_SEND(dble(Vnum)*atomSize/systemVol,1,MPI_DOUBLE_PRECISION, i, 1406,comm, ierr)
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
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 1399, comm, ierr)

	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		
		!call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 1400, comm, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,1405,comm,status,ierr)
		
		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		end if

	end do

	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 1400, comm, ierr)

	deallocate(defectsSend)
	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)

	!call MPI_RECV(computeVConc,1,MPI_DOUBLE_PRECISION,MASTER,1406,comm,status,ierr)
end if

call MPI_BCAST(computeVConc, 1, MPI_DOUBLE_PRECISION, MASTER, comm,ierr)

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
use mod_constants
implicit none

include 'mpif.h'

integer i, j, k, same, status(MPI_STATUS_SIZE), numDefectsRecv, numDefectsSend
integer products(numSpecies)
type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList

double precision, allocatable :: defectsRecv(:,:)
double precision, allocatable :: defectsSend(:,:)

!variables for computation statistics
integer reactionsCoarse, reactionsFine

!Variables for defect counting and concentration
integer VNum, SIANum, LoopSize(3), totalVac, totalSIA, CuNum, totalCu
integer totalVoid, totalLoop, VoidNum

interface
	subroutine findDefectInList(defectCurrent, defectPrev, products)
		use DerivedType
		use mod_constants
		implicit none
		type(defect), pointer :: defectCurrent, defectPrev
		integer products(numSpecies)
	end subroutine
end interface

!initialize outputDefectList
allocate(outputDefectList)
allocate(outputDefectList%defectType(numSpecies))
nullify(outputDefectList%next)
do i=1,numSpecies
	outputDefectList%defectType(i)=0
end do
outputDefectList%cellNumber=0
outputDefectList%num=0

!call MPI_ALLREDUCE(numImpAnn,totalImpAnn,2,MPI_INTEGER,MPI_SUM,comm,ierr)

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
	
		call MPI_RECV(numDefectsRecv,1,MPI_INTEGER,i,2399,comm,status,ierr)

		allocate(defectsRecv(numSpecies+1,numDefectsRecv))

		call MPI_RECV(defectsRecv,(numSpecies+1)*numDefectsRecv,MPI_DOUBLE_PRECISION,i,2400,comm,status,ierr)
		
		do j=1,numDefectsRecv
			!call MPI_RECV(defectsRecv,numSpecies+1,MPI_DOUBLE_PRECISION,i,2400,comm,status,ierr)
			
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
		!call MPI_SEND(dble(SIAnum)*atomSize/systemVol,1,MPI_DOUBLE_PRECISION, i, 2406,comm, ierr)
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
		
	call MPI_SEND(numDefectsSend, 1, MPI_INTEGER, MASTER, 2399, comm, ierr)

	allocate(defectsSend(numSpecies+1,numDefectsSend))
	
	defectCurrentList=>outputDefectList%next
	i=0
	
	do while(associated(defectCurrentList))
		i=i+1
	
		do j=1,numSpecies
			defectsSend(j,i)=defectCurrentList%defectType(j)
		end do
		
		defectsSend(numSpecies+1,i)=defectCurrentList%num
		
		!call MPI_SEND(defectsSend, numSpecies+1, MPI_DOUBLE_PRECISION, MASTER, 2400, comm, ierr)
		
		!This just pauses the sending until the master has recieved all of the above information
		!call MPI_RECV(buffer,1,MPI_INTEGER,MASTER,2405,comm,status,ierr)
		
		defectCurrentList=>defectCurrentList%next
		
		if(i==numDefectsSend .AND. associated(defectCurrentList)) then
			write(*,*) 'error outputDefectList size does not match numDefectsSend'
		end if

	end do

	call MPI_SEND(defectsSend, (numSpecies+1)*numDefectsSend, MPI_DOUBLE_PRECISION, MASTER, 2400, comm, ierr)
	
	call countReactionsCoarse(reactionsCoarse)
	call countReactionsFine(reactionsFine)
	
	!call MPI_RECV(computeIConc,1,MPI_DOUBLE_PRECISION,MASTER,2406,comm,status,ierr)
endif

call MPI_BCAST(computeIConc, 1, MPI_DOUBLE_PRECISION, MASTER, comm,ierr)

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
