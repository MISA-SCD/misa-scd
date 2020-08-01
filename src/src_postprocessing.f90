!***************************************************************************************************
!>subroutine outputDefectsXYZ()
!Outputs into file: rawdat.out. These contain the complete defect populations per volume element.
!Compiles data from local as well as global processors.
!***************************************************************************************************
subroutine outputDefectsXYZ()
	use mod_constants
	use mod_structures
	implicit none
	include 'mpif.h'

	integer :: i, j, k, status(MPI_STATUS_SIZE), defectCount, numRecv
	integer :: numScluster, numVoid, numLoop	!total number of solute/V/SIA/mSnV clusters in per volume element
	double precision :: volTemp
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

		write(RAWDAT,*) 'time', elapsedTime, 'DPA', DPA, 'steps', step

		!Write defects in processor 0 to file
		write(RAWDAT,*) 'processor', myProc%taskid
		do i=1,numCells
			defectCurrent=>defectList(i)%next
			write(RAWDAT,*) 'cell', i,'globalCell', myMesh(i)%globalCell
			write(RAWDAT,*) 'defects (Cu V SIA_m SIA_im num)'

			numScluster=0
			numVoid=0
			numLoop=0

			do while(associated(defectCurrent))
				write(RAWDAT,*) defectCurrent%defectType, defectCurrent%num
				if(defectCurrent%defectType(1) > minSCluster) then
					numScluster=numScluster+defectCurrent%num
				else if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2) > minVoid) then
					numVoid=numVoid+defectCurrent%num
				else if(defectCurrent%defectType(3)>minLoop .OR. defectCurrent%defectType(4)>minLoop) then
					numLoop=numLoop+defectCurrent%num
				end if
				defectCurrent=>defectCurrent%next
			end do

			write(RAWDAT,*) 'concentration of point defects (Cu/V/SIA) and clusters (Cu cluster/Void/Loop)'
			write(RAWDAT,*) dble(numScluster)*volTemp, dble(numVoid)*volTemp, dble(numLoop)*volTemp
			write(RAWDAT,*)
		end do
		write(RAWDAT,*)
		write(RAWDAT,*)

		!Other processors will send information about defects contained in their mesh. This processor
		!must output all information to data file (can't do it from other processors). Information should
		!be output in the same format for the master and slave processors.
		do i=1,myProc%numtasks-1
			write(RAWDAT,*) 'processor', i

			do j=1,numCellRecv(i+1)

				numRecv=0
				call MPI_PROBE(i, 100+j,comm,status,ierr)
				call MPI_GET_COUNT(status,MPI_INTEGER,numRecv,ierr)

				numRecv=numRecv/(numSpecies+1)
				allocate(cellRecv(numSpecies+1,numRecv))
				call MPI_RECV(cellRecv,(numSpecies+1)*numRecv,MPI_INTEGER,i,100+j,comm,status,ierr)

				write(RAWDAT,*) 'cell', j,'globalCell', cellRecv(2,1)
				write(RAWDAT,*) 'defects (Cu V SIA_m SIA_im num)'

				do k=1,cellRecv(1,1)
					write(RAWDAT,*) cellRecv(:,k+1)
				end do

				write(RAWDAT,*) 'concentration of point defects (Cu/V/SIA) and clusters (Cu cluster/Void/Loop)'
				write(RAWDAT,*) dble(cellRecv(3,1))*volTemp,dble(cellRecv(4,1))*volTemp,dble(cellRecv(5,1))*volTemp
				write(RAWDAT,*)

				deallocate(cellRecv)
			end do
		end do
		write(RAWDAT,*)
		write(RAWDAT,*)
	end if

end subroutine

!*******************************************************************************************************
!> Subroutine outputDefectsTotal()
!Outputs into file: totdat.out. Outputs a list of all defects in system (defect type, num), regardless of volume element.
!Compiles such a list using message passing between processors.
!Used to find total defect populations in simulation volume, not spatial distribution of defects.
!*******************************************************************************************************
subroutine outputDefectsTotal()
	use mod_structures
	use mod_constants
	implicit none
	include 'mpif.h'

	integer :: i, j, k, status(MPI_STATUS_SIZE)
	integer :: products(numSpecies), same
	type(defect), pointer :: defectCurrent, defectPrevList, defectCurrentList, outputDefectList
	type(cascade), pointer :: cascadeCurrent
	integer, allocatable :: defectsRecv(:,:)
	integer, allocatable :: defectsSend(:,:)
	integer :: numDefectsSend
	integer, allocatable :: numDefectsRecv(:)
	!number of point defects (solute atom, vacancy, self-interstitial atom)
	integer :: pointS, pointV, pointSIA
	!total retained vacancies/self-interstitial atoms in the whole system
	integer :: totalVac, totalSIA
	!total number of solute atoms/vacancies/SIAs in clusters
	integer :: numS, numV, numSIA, numSV
	!total number of solute/V/SIA/mSnV clusters in the whole system
	integer :: numScluster, numVoid, numLoop, numSVcluster
	!Number density of clusters
	double precision :: denScluster, denVoid, denLoop, denSVcluster
	!Average radius of clusters
	double precision :: radiusScluster, radiusVoid, radiusLoop, radiusSVcluster
	!Average size of clusters
	double precision :: sizeScluster, sizeVoid, sizeLoop, sizeSVcluster
	double precision :: VRetained, VAnnihilated, conPointV, conPointSIA

	interface
		subroutine findDefectInList(defectCurrent, defectPrev, products)
			use mod_structures
			use mod_constants
			implicit none
			type(defect), pointer, intent(inout) :: defectCurrent, defectPrev
			integer, intent(in) :: products(numSpecies)
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

		write(TOTDAT,*) 'defects (Cu, V, SIA_m, SIA_im, num)'
		defectCurrentList=>outputDefectList%next
		do while(associated(defectCurrentList))

			if(defectCurrentList%defectType(1) /= 0) then	!Cu/CuV cluster

				if(defectCurrentList%defectType(1)==1 .AND. defectCurrentList%defectType(2)==0) then
					pointS=defectCurrentList%num
				end if

				if(defectCurrentList%defectType(1) > minSCluster) then
					numS=numS+defectCurrentList%defectType(1)*defectCurrentList%num
					radiusScluster=radiusScluster+dble(defectCurrentList%num)*&
							(3*dble(defectCurrentList%defectType(1))*atomSize_Cu/(4*pi))**(1d0/3d0)
					numScluster=numScluster+defectCurrentList%num
				end if

				if(defectCurrentList%defectType(2)>0 .AND. &
						(defectCurrentList%defectType(1)+defectCurrentList%defectType(2)) > minSV) then
					numSV=numSV+(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*defectCurrentList%num
					radiusSVcluster = radiusSVcluster+dble(defectCurrentList%num)*&
							(3*dble(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*&
									atomSize_Cu/(4*pi))**(1d0/3d0)
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
			write(TOTDAT,*) defectCurrentList%defectType, defectCurrentList%num

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
		sizeSVcluster=dble(numSV)/dble(numSVcluster)

		conPointV=dble(pointV)/systemVol*atomSize
		conPointSIA=dble(pointSIA)/systemVol*atomSize

		VRetained = dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
		VAnnihilated = dble(totalImpAnn(2))/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))

		!Output totdat.out
		write(TOTDAT,*)	'numCluster (S/Void/Loop/SV):', numScluster, numVoid, numLoop, numSVcluster
		write(TOTDAT,*)	'NumberDensity (m-3) (S/Void/Loop/SV):',denScluster*1d27,denVoid*1d27,denLoop*1d27,denSVcluster*1d27
		write(TOTDAT,*)	'Concentration (S/Void/Loop/SV):', denScluster*atomSize_Cu, denVoid*atomSize, denLoop*atomSize,&
				denSVcluster*atomSize_Cu
		write(TOTDAT,*)	'AverageRadius (nm) (S/Void/Loop/SV):', radiusScluster, radiusVoid, radiusLoop, radiusSVcluster
		write(TOTDAT,*)	'AverageSize (S/Void/Loop/SV):', sizeScluster, sizeVoid, sizeLoop, sizeSVcluster
		write(TOTDAT,*) 'ConcenPointDefects (V/SIA):', conPointV, conPointSIA
		write(TOTDAT,*) 'PercentVRetained',VRetained,'PercentVAnnihilated',VAnnihilated
		write(TOTDAT,*)
		write(TOTDAT,*)

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
!***********************************************************************
!***********************************************************************
subroutine outputTotal()
	use mod_structures
	use mod_constants
	implicit none
	include 'mpif.h'

	integer :: i
	type(defect), pointer :: defectCurrent
	type(cascade), pointer :: cascadeCurrent
	integer :: sendBuff(8), recvBuff(8)
	!number of point defects (solute atom, vacancy, self-interstitial atom)
	integer :: pointS, pointV, pointSIA
	!total retained vacancies/self-interstitial atoms in the whole system
	integer :: totalVac, totalSIA
	!total number of solute atoms/vacancies/SIAs in clusters
	integer :: numS, numV, numSIA, numSV
	!total number of solute/V/SIA/mSnV clusters in the whole system
	integer :: numScluster, numVoid, numLoop, numSVcluster
	!Number density of clusters
	double precision :: denScluster, denVoid, denLoop, denSVcluster
	!Average radius of clusters
	double precision :: radiusScluster, radiusVoid, radiusLoop, radiusSVcluster
	!Average size of clusters
	double precision :: sizeScluster, sizeVoid, sizeLoop, sizeSVcluster
	double precision :: VRetained, VAnnihilated, conPointV, conPointSIA

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

	do i=1,numCells

		defectCurrent=>defectList(i)%next
		do while(associated(defectCurrent))

			if(defectCurrent%defectType(1) /= 0) then	!Cu/CuV cluster
				if(defectCurrent%defectType(1)==1 .AND. defectCurrent%defectType(2)==0) then
					pointS=defectCurrent%num
				end if

				if(defectCurrent%defectType(1) > minSCluster) then
					numS=numS+defectCurrent%defectType(1)*defectCurrent%num
					radiusScluster=radiusScluster+dble(defectCurrent%num)*&
							(3*dble(defectCurrent%defectType(1))*atomSize/(4*pi))**(1d0/3d0)
					numScluster=numScluster+defectCurrent%num
				end if

				if((defectCurrent%defectType(1)+defectCurrent%defectType(2)) > minSV) then
					numSV=numSV+max(defectCurrent%defectType(1), defectCurrent%defectType(2))*defectCurrent%num
					radiusSVcluster = radiusSVcluster+dble(defectCurrent%num)*&
							(3*dble(max(defectCurrent%defectType(1), defectCurrent%defectType(2)))*&
									atomSize/(4*pi))**(1d0/3d0)
					numSVcluster=numSVcluster+defectCurrent%num
					if(defectCurrent%defectType(2) /= 0) then
						totalVac=totalVac+defectCurrent%defectType(2)*defectCurrent%num
					end if
				end if

			else if(defectCurrent%defectType(2) /= 0) then !V cluster

				totalVac=totalVac+defectCurrent%defectType(2)*defectCurrent%num

				if(defectCurrent%defectType(2)==1) then
					pointV=defectCurrent%num
				end if

				if(defectCurrent%defectType(2) > minVoid) then
					numV=numV+defectCurrent%defectType(2)*defectCurrent%num
					radiusVoid=radiusVoid+dble(defectCurrent%num)*&
							(3d0*dble(defectCurrent%defectType(2))*atomSize/(4d0*pi))**(1d0/3d0)
					numVoid=numVoid+defectCurrent%num
				end if

			else if(defectCurrent%defectType(3) /= 0) then	!Loop

				if(defectCurrent%defectType(3) ==1) then
					pointSIA=defectCurrent%num
				end if

				if(defectCurrent%defectType(3) > minLoop) then
					numSIA=numSIA+defectCurrent%defectType(3)*defectCurrent%num
					radiusLoop=radiusLoop+dble(defectCurrent%num)*&
							(dble(defectCurrent%defectType(3))*atomSize/(pi*burgers))**(1d0/2d0)
					numLoop=numLoop+defectCurrent%num
				end if

			else if(defectCurrent%defectType(4) /= 0) then

				if(defectCurrent%defectType(4) > minLoop) then
					numSIA=numSIA+defectCurrent%defectType(4)*defectCurrent%num
					radiusLoop=radiusLoop+dble(defectCurrent%num)*&
							(dble(defectCurrent%defectType(4))*atomSize/(pi*burgers))**(1d0/2d0)
					numLoop=numLoop+defectCurrent%num
				end if
			end if
			defectCurrent=>defectCurrent%next
		end do
	end do

	!Count defects in fine mesh
	if(meshingType=='adaptive') then
		cascadeCurrent=>ActiveCascades
		do while(associated(cascadeCurrent))
			do i=1,numCellsCascade
				defectCurrent=>cascadeCurrent%localDefects(i)%next
				do while(associated(defectCurrent))

					if(defectCurrent%defectType(1) /= 0) then	!Cu/CuV cluster
						if(defectCurrent%defectType(1)==1 .AND. defectCurrent%defectType(2)==0) then
							pointS=defectCurrent%num
						end if

						if(defectCurrent%defectType(1) > minSCluster) then
							numS=numS+defectCurrent%defectType(1)*defectCurrent%num
							radiusScluster=radiusScluster+dble(defectCurrent%num)*&
									(3*dble(defectCurrent%defectType(1))*atomSize/(4*pi))**(1d0/3d0)
							numScluster=numScluster+defectCurrent%num
						end if

						if((defectCurrent%defectType(1)+defectCurrent%defectType(2)) > minSV) then
							numSV=numSV+max(defectCurrent%defectType(1), defectCurrent%defectType(2))*defectCurrent%num
							radiusSVcluster = radiusSVcluster+dble(defectCurrent%num)*&
									(3*dble(max(defectCurrent%defectType(1), defectCurrent%defectType(2)))*&
											atomSize/(4*pi))**(1d0/3d0)
							numSVcluster=numSVcluster+defectCurrent%num
							if(defectCurrent%defectType(2) /= 0) then
								totalVac=totalVac+defectCurrent%defectType(2)*defectCurrent%num
							end if
						end if

					else if(defectCurrent%defectType(2) /= 0) then !V cluster

						totalVac=totalVac+defectCurrent%defectType(2)*defectCurrent%num

						if(defectCurrent%defectType(2)==1) then
							pointV=defectCurrent%num
						end if

						if(defectCurrent%defectType(2) > minVoid) then
							numV=numV+defectCurrent%defectType(2)*defectCurrent%num
							radiusVoid=radiusVoid+dble(defectCurrent%num)*&
									(3d0*dble(defectCurrent%defectType(2))*atomSize/(4d0*pi))**(1d0/3d0)
							numVoid=numVoid+defectCurrent%num
						end if

					else if(defectCurrent%defectType(3) /= 0) then	!Loop

						if(defectCurrent%defectType(3) ==1) then
							pointSIA=defectCurrent%num
						end if

						if(defectCurrent%defectType(3) > minLoop) then
							numSIA=numSIA+defectCurrent%defectType(3)*defectCurrent%num
							radiusLoop=radiusLoop+dble(defectCurrent%num)*&
									(dble(defectCurrent%defectType(3))*atomSize/(pi*burgers))**(1d0/2d0)
							numLoop=numLoop+defectCurrent%num
						end if

					else if(defectCurrent%defectType(4) /= 0) then

						if(defectCurrent%defectType(4) > minLoop) then
							numSIA=numSIA+defectCurrent%defectType(4)*defectCurrent%num
							radiusLoop=radiusLoop+dble(defectCurrent%num)*&
									(dble(defectCurrent%defectType(4))*atomSize/(pi*burgers))**(1d0/2d0)
							numLoop=numLoop+defectCurrent%num
						end if
					end if

					defectCurrent=>defectCurrent%next
				end do
			end do
			cascadeCurrent=>cascadeCurrent%next
		end do
	end if

	sendBuff(1)=pointV	!number of vacancies
	sendBuff(2)=pointSIA	!number of SIAs
	sendBuff(3)=numS
	sendBuff(4)=numScluster
	sendBuff(5)=numV
	sendBuff(6)=numVoid
	sendBuff(7)=numSIA
	sendBuff(8)=numLoop

	call MPI_REDUCE(sendBuff,recvBuff, 8, MPI_INTEGER, MPI_SUM, 0,comm, ierr)

	!Output totdat.out
	if(myProc%taskid==MASTER) then

		!number of point defects
		!pointS=0
		pointV=recvBuff(1)
		pointSIA=recvBuff(2)

		!totalVac=0
		!totalSIA=0

		!total number of solute (Cu) atoms/vacancies/self-interstitial atoms in clusters
		numS=recvBuff(3)
		numV=recvBuff(5)
		numSIA=recvBuff(7)

		!total number of solute (Cu)/V/SIA clusters in the whole system
		numScluster=recvBuff(4)
		numVoid=recvBuff(6)
		numLoop=recvBuff(8)

		!VRetained=0d0
		!VAnnihilated=0d0



		denScluster = dble(numScluster)/systemVol
		denVoid = dble(numVoid)/systemVol
		denLoop = dble(numLoop)/systemVol
		!denSVcluster = dble(numSVcluster)/systemVol

		radiusScluster=(3*(dble(numS)/dble(numScluster))*atomSize/(4*pi))**(1d0/3d0)
		radiusVoid=(3d0*(dble(numV)/dble(numVoid))*atomSize/(4d0*pi))**(1d0/3d0)
		radiusLoop=((dble(numSIA)/dble(numLoop))*atomSize/(pi*burgers))**(1d0/2d0)
		!radiusSVcluster=(3d0*(dble(numSV)/dble(numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)
		!radiusScluster=radiusScluster/dble(numScluster)
		!radiusVoid=radiusVoid/dble(numVoid)
		!radiusLoop=radiusLoop/dble(numLoop)
		!radiusSVcluster=radiusSVcluster/dble(numSVcluster)

		sizeScluster=dble(numS)/dble(numScluster)
		sizeVoid=dble(numV)/dble(numVoid)
		sizeLoop=dble(numSIA)/dble(numLoop)
		!sizeSVcluster=dble(numSV)/dble(numSVcluster)

		conPointV=dble(pointV)/systemVol*atomSize
		conPointSIA=dble(pointSIA)/systemVol*atomSize

		!VRetained = dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
		!VAnnihilated = dble(totalImpAnn(2))/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))

		!Output totdat.out
		write(TOTDAT,*)	'numCluster (S/Void/Loop):', numScluster, numVoid, numLoop
		write(TOTDAT,*)	'NumberDensity (m-3) (S/Void/Loop):', denScluster*1d27, denVoid*1d27,denLoop*1d27
		write(TOTDAT,*)	'Concentration (S/Void/Loop):', denScluster*atomSize, denVoid*atomSize, denLoop*atomSize
		write(TOTDAT,*)	'AverageRadius (nm) (S/Void/Loop):', radiusScluster, radiusVoid, radiusLoop
		write(TOTDAT,*)	'AverageSize (S/Void/Loop):', sizeScluster, sizeVoid, sizeLoop
		write(TOTDAT,*) 'ConcenPointDefects (V/SIA):', conPointV, conPointSIA
		!write(TOTDAT,*) 'PercentVRetained',VRetained,'PercentVAnnihilated',VAnnihilated
		write(TOTDAT,*)
		write(TOTDAT,*)

	end if

end subroutine
