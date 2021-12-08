!***************************************************************************************************
!>subroutine outputDefectsXYZ()
!Outputs into file: xyzdat.out. These contain the defect populations per volume element.
!***************************************************************************************************
subroutine outputDefectsXYZ()
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none
	include 'mpif.h'

	integer :: numDefects, cell, index, i, j, k
	type(defect), pointer :: defectCurrent
	integer, allocatable :: cellDefects(:)
	integer :: status(MPI_STATUS_SIZE), recvCells, count, numRecv
	integer, allocatable :: sendBuff(:,:), recvBuff(:,:)
	integer :: pointS, pointV, pointSIA		!number of point defects (solute atom, vacancy, self-interstitial atom)
	!integer :: totalVac, totalSIA			!total retained vacancies/self-interstitial atoms in the whole system
	integer :: numS, numV, numSIA, numSV	!total number of solute atoms/vacancies/SIAs in clusters
	integer :: numScluster, numVoid, numLoop, numSVcluster	!total number of solute/V/SIA/mSnV clusters in the whole system
	!double precision :: denScluster, denVoid, denLoop, denSVcluster			!Number density of clusters
	double precision :: radiusScluster, radiusVoid, radiusLoop, radiusSVcluster	!Average radius of clusters

	if(myProc%taskid/=MASTER) then
		numDefects=0
		allocate(cellDefects(numCells))
		do cell=1, numCells
			cellDefects(cell)=0
			nullify(defectCurrent)
			defectCurrent=>defectList(cell)%next
			do while(associated(defectCurrent))
				cellDefects(cell)=cellDefects(cell)+1
				numDefects=numDefects+1
				defectCurrent=>defectCurrent%next
			end do
		end do
		allocate(sendBuff(SPECIES+1,numDefects+numCells))

		index=0
		do cell=1, numCells
			if(cell==1) then
				index=1
			else
				index=index+cellDefects(cell-1)+1
			end if
			sendBuff(1,index)=numCells
			sendBuff(2,index)=myMesh(cell)%globalCell	!global cellID of this cell
			sendBuff(3,index)=cell						!local cellID
			sendBuff(4,index)=cellDefects(cell)			!number of defects in this cell
			sendBuff(5:SPECIES+1,index)=0

			i=0
			nullify(defectCurrent)
			defectCurrent=>defectList(i)%next
			do while(associated(defectCurrent))
				i=i+1
				sendBuff(1:SPECIES,index+i)=defectCurrent%defectType(1:SPECIES)
				sendBuff(SPECIES+1,index+i)=defectCurrent%num
				defectCurrent=>defectCurrent%next
			end do
		end do
		call MPI_SEND(sendBuff,(SPECIES+1)*(numDefects+numCells),MPI_INTEGER,MASTER,100+myProc%taskid,comm,ierr)
		deallocate(sendBuff)
	else
		write(XYZFILE,*) '***********************************************************************'
		write(XYZFILE,*) 'time', elapsedTime, 'steps', step

		!Write defects in processor 0 to file
		do cell=1, numCells
			write(XYZFILE,*) 'processor', myProc%taskid, 'globalCell', myMesh(cell)%globalCell, 'cell', cell
			!number of point defects
			pointS=0
			pointV=0
			pointSIA=0
			!total number of solute (Cu) atoms/vacancies/self-interstitial atoms in clusters
			numV=0
			numSIA=0
			numS=0
			numSV=0
			!total number of solute (Cu)/V/SIA clusters in the whole system
			numVoid=0
			numLoop=0
			numScluster=0
			numSVcluster=0
			!average radius of clusters
			radiusVoid=0d0
			radiusLoop=0d0
			radiusScluster=0d0
			radiusSVcluster=0d0

			nullify(defectCurrent)
			defectCurrent=>defectList(cell)%next
			do while(associated(defectCurrent))
				write(XYZFILE,*) defectCurrent%defectType, defectCurrent%num
				if(defectCurrent%defectType(1) /=0) then	!Cu/Cu_Vac cluster
					if(defectCurrent%defectType(1)==1 .AND. defectCurrent%defectType(2)==0) then
						pointS=pointS+defectCurrent%num
					end if
					if(defectCurrent%defectType(1)>minSCluster .AND. defectCurrent%defectType(2)==0) then	!Cu cluster
						numS=numS+defectCurrent%defectType(1)*defectCurrent%num
						numScluster=numScluster+defectCurrent%num
					else if(defectCurrent%defectType(2)/=0 .AND. (defectCurrent%defectType(1)+defectCurrent%defectType(2))>minSV) then
						numSV=numSV+(defectCurrent%defectType(1)+defectCurrent%defectType(2))*defectCurrent%num
						numSVcluster=numSVcluster+defectCurrent%num
					end if
				else if(defectCurrent%defectType(2)/=0) then	!V cluster
					if(defectCurrent%defectType(2)==1) then
						pointV=pointV+defectCurrent%num
					end if
					if(defectCurrent%defectType(2) > minVoid) then
						numV=numV+defectCurrent%defectType(2)*defectCurrent%num
						numVoid=numVoid+defectCurrent%num
					end if
				else if(defectCurrent%defectType(3)/=0) then
					if(defectCurrent%defectType(3)==1) then
						pointSIA=pointSIA+defectCurrent%num
					end if
					if(defectCurrent%defectType(3) >minLoop) then
						numSIA=numSIA+defectCurrent%defectType(3)*defectCurrent%num
						numLoop=numLoop+defectCurrent%num
					end if
				else if(defectCurrent%defectType(4)/=0) then
					if(defectCurrent%defectType(4) >minLoop) then
						numSIA=numSIA+defectCurrent%defectType(4)*defectCurrent%num
						numLoop=numLoop+defectCurrent%num
					end if
				end if
				defectCurrent=>defectCurrent%next
			end do

			radiusLoop=((dble(numSIA)/dble(numLoop))*atomSize/(pi*burgers))**(1d0/2d0)
			radiusVoid=(3d0*(dble(numV)/dble(numVoid))*atomSize/(4d0*pi))**(1d0/3d0)
			radiusScluster=(3*(dble(numS)/dble(numScluster))*atomSize/(4*pi))**(1d0/3d0)
			radiusSVcluster=(3d0*(dble(numSV)/dble(numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)

			write(XYZFILE,*) 'numCluster (Loop/Void/S/SV/precipitate):'
			write(XYZFILE,*) numLoop, numVoid, numScluster, numSVcluster, (numScluster+numSVcluster)
			write(XYZFILE,*) 'NumberDensity (m-3) (Loop/Void/S/SV/precipitate):'
			write(XYZFILE,*) dble(numLoop)/systemVol*1d27,dble(numVoid)/systemVol*1d27,dble(numScluster)/systemVol*1d27, &
					dble(numSVcluster)/systemVol*1d27, dble(numScluster+numSVcluster)/systemVol*1d27
			write(XYZFILE,*) 'AverageRadius (nm) (Loop/Void/S/SV/precipitate):'
			write(XYZFILE,*) radiusLoop, radiusVoid, radiusScluster, radiusSVcluster, &
					(3d0*(dble(numS+numSV)/dble(numScluster+numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)
			write(XYZFILE,*) 'NumberDensity (m-3) (SIA/V/S):'
			write(XYZFILE,*) dble(pointSIA)/systemVol*1d27, dble(pointV)/systemVol*1d27, dble(pointS)/systemVol*1d27
			write(XYZFILE,*)
		end do
		write(XYZFILE,*)

		!recv defect populations and output then to xyzdat.out
		do i=1,myProc%numtasks-1
			count=0
			call MPI_PROBE(i, 100+i,comm,status,ierr)
			call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)
			numRecv=count/(SPECIES+1)
			allocate(recvBuff(SPECIES+1, numRecv))
			call MPI_RECV(recvBuff,(SPECIES+1)*numRecv,MPI_INTEGER,i,100+i,comm,status,ierr)

			recvCells=recvBuff(1,1)
			do cell=1, recvCells
				if(cell==1) then
					index=1
				else
					index=index+recvBuff(4,index)+1		!recvBuff(4,index): number defects in this cell
				end if
				write(XYZFILE,*) 'processor', i, 'globalCell', recvBuff(2,index), 'cell', recvBuff(3,index)
				!number of point defects
				pointS=0
				pointV=0
				pointSIA=0
				!total number of solute (Cu) atoms/vacancies/self-interstitial atoms in clusters
				numV=0
				numSIA=0
				numS=0
				numSV=0
				!total number of solute (Cu)/V/SIA clusters in the whole system
				numVoid=0
				numLoop=0
				numScluster=0
				numSVcluster=0
				!average radius of clusters
				radiusVoid=0d0
				radiusLoop=0d0
				radiusScluster=0d0
				radiusSVcluster=0d0

				do k=1, recvBuff(4,index)
					write(XYZFILE,*) recvBuff(1:SPECIES+1,index+k)
					if(recvBuff(1,index+k) /=0) then	!Cu/Cu_Vac cluster
						if(recvBuff(1,index+k)==1 .AND. recvBuff(2,index+k)==0) then
							pointS=pointS+recvBuff(SPECIES+1,index+k)
						end if
						if(recvBuff(1,index+k)>minSCluster .AND. recvBuff(2,index+k)==0) then	!Cu cluster
							numS=numS+recvBuff(1,index+k)*recvBuff(SPECIES+1,index+k)
							numScluster=numScluster+recvBuff(SPECIES+1,index+k)
						else if(recvBuff(2,index+k)/=0 .AND. (recvBuff(1,index+k)+recvBuff(2,index+k))>minSV) then
							numSV=numSV+(recvBuff(1,index+k)+recvBuff(2,index+k))*recvBuff(SPECIES+1,index+k)
							numSVcluster=numSVcluster+recvBuff(SPECIES+1,index+k)
						end if
					else if(recvBuff(2,index+k)/=0) then	!V cluster
						if(recvBuff(2,index+k)==1) then
							pointV=pointV+recvBuff(SPECIES+1,index+k)
						end if
						if(recvBuff(2,index+k) > minVoid) then
							numV=numV+recvBuff(2,index+k)*recvBuff(SPECIES+1,index+k)
							numVoid=numVoid+recvBuff(SPECIES+1,index+k)
						end if
					else if(recvBuff(3,index+k)/=0) then
						if(recvBuff(3,index+k)==1) then
							pointSIA=pointSIA+recvBuff(SPECIES+1,index+k)
						end if
						if(recvBuff(3,index+k) >minLoop) then
							numSIA=numSIA+recvBuff(3,index+k)*recvBuff(SPECIES+1,index+k)
							numLoop=numLoop+recvBuff(SPECIES+1,index+k)
						end if
					else if(recvBuff(4,index+k)/=0) then
						if(recvBuff(4,index+k) >minLoop) then
							numSIA=numSIA+recvBuff(4,index+k)*defectCurrent%num
							numLoop=numLoop+recvBuff(SPECIES+1,index+k)
						end if
					end if
				end do

				radiusLoop=((dble(numSIA)/dble(numLoop))*atomSize/(pi*burgers))**(1d0/2d0)
				radiusVoid=(3d0*(dble(numV)/dble(numVoid))*atomSize/(4d0*pi))**(1d0/3d0)
				radiusScluster=(3*(dble(numS)/dble(numScluster))*atomSize/(4*pi))**(1d0/3d0)
				radiusSVcluster=(3d0*(dble(numSV)/dble(numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)

				write(XYZFILE,*) 'numCluster (Loop/Void/S/SV/precipitate):'
				write(XYZFILE,*) numLoop, numVoid, numScluster, numSVcluster, (numScluster+numSVcluster)
				write(XYZFILE,*) 'NumberDensity (m-3) (Loop/Void/S/SV/precipitate):'
				write(XYZFILE,*) dble(numLoop)/systemVol*1d27,dble(numVoid)/systemVol*1d27,dble(numScluster)/systemVol*1d27, &
						dble(numSVcluster)/systemVol*1d27, dble(numScluster+numSVcluster)/systemVol*1d27
				write(XYZFILE,*) 'AverageRadius (nm) (Loop/Void/S/SV/precipitate):'
				write(XYZFILE,*) radiusLoop, radiusVoid, radiusScluster, radiusSVcluster, &
						(3d0*(dble(numS+numSV)/dble(numScluster+numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)
				write(XYZFILE,*) 'NumberDensity (m-3) (SIA/V/S):'
				write(XYZFILE,*) dble(pointSIA)/systemVol*1d27, dble(pointV)/systemVol*1d27, dble(pointS)/systemVol*1d27
				write(XYZFILE,*)
			end do
			write(XYZFILE,*)
			deallocate(recvBuff)
		end do
		write(XYZFILE,*)
		write(XYZFILE,*)
	end if

end subroutine

!*******************************************************************************************************
!> Subroutine outputDefectsTotal(simStatus)
!Outputs into file: totdat.out. Outputs a list of all defects in system (defect type, num), regardless of volume element.
!Compiles such a list using message passing between processors.
!Used to find total defect populations in simulation volume, not spatial distribution of defects.
!*******************************************************************************************************
subroutine outputDefectsTotal(simStatus)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none
	include 'mpif.h'

	integer, intent(in) :: simStatus		!<1: 'irradiation', 0: 'anneal'
	integer :: i, j, k, status(MPI_STATUS_SIZE)
	integer :: products(SPECIES), same
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
	integer :: numS, numV, numSIA, numSV, numCu
	!total number of solute/V/SIA/mSnV clusters in the whole system
	integer :: numScluster, numVoid, numLoop, numSVcluster, numCuCluster
	!Number density of clusters
	double precision :: denScluster, denVoid, denLoop, denSVcluster, denCuCluster
	!Average radius of clusters
	double precision :: radiusScluster, radiusVoid, radiusLoop, radiusSVcluster, radiusCuCluster
	!Average size of clusters
	double precision :: sizeScluster, sizeVoid, sizeLoop, sizeSVcluster, sizeCuCluster
	double precision :: VRetained, VAnnihilated, conPointV, conPointSIA

	interface
		subroutine findDefectInList(defectCurrent, defectPrev, products)
			use mod_constants
			use mod_structures
			use mod_globalVariables
			implicit none
			type(defect), pointer, intent(inout) :: defectCurrent, defectPrev
			integer, intent(in) :: products(SPECIES)
		end subroutine
	end interface

	!initialize outputDefectList
	allocate(outputDefectList)
	allocate(outputDefectList%defectType(SPECIES))
	nullify(outputDefectList%next)
	do i=1,SPECIES
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
				do j=1,SPECIES
					if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
						same=same+1
					end if
				end do

				if(same == SPECIES) then
					!if the defect is already present in the list
					defectCurrentList%num=defectCurrentList%num+defectCurrent%num
				else    !add a defect to the middle of the list

					!if the defect is to be inserted in the list
					nullify(defectPrevList%next)
					allocate(defectPrevList%next)
					nullify(defectPrevList%next%next)
					defectPrevList=>defectPrevList%next
					allocate(defectPrevList%defectType(SPECIES))
					defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
					defectPrevList%num=defectCurrent%num

					do j=1,SPECIES
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
				allocate(defectPrevList%defectType(SPECIES))
				defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
				defectPrevList%num=defectCurrent%num

				do j=1,SPECIES
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
						do j=1,SPECIES
							if(defectCurrentList%defectType(j)==defectCurrent%defectType(j)) then
								same=same+1
							end if
						end do

						if(same == SPECIES) then
							!if the defect is already present in the list
							defectCurrentList%num=defectCurrentList%num+defectCurrent%num
						else    !add a defect to the middle of the list
							!if the defect is to be inserted in the list
							nullify(defectPrevList%next)
							allocate(defectPrevList%next)
							nullify(defectPrevList%next%next)
							defectPrevList=>defectPrevList%next
							allocate(defectPrevList%defectType(SPECIES))
							defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
							defectPrevList%num=defectCurrent%num

							do j=1,SPECIES
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
						allocate(defectPrevList%defectType(SPECIES))
						defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
						defectPrevList%num=defectCurrent%num

						do j=1,SPECIES
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
		allocate(defectsSend(SPECIES+1,numDefectsSend))
		i=0
		do while(associated(defectCurrentList))
			i=i+1
			do j=1,SPECIES
				defectsSend(j,i)=defectCurrentList%defectType(j)
			end do
			defectsSend(SPECIES+1,i)=defectCurrentList%num
			defectCurrentList=>defectCurrentList%next
			if(i==numDefectsSend .AND. associated(defectCurrentList)) then
				write(*,*) 'error outputDefectList size does not match numDefectsSend'
			end if
		end do
		call MPI_SEND(defectsSend, (SPECIES+1)*numDefectsSend, MPI_INTEGER, MASTER, 400, comm, ierr)
		deallocate(defectsSend)
	else	!MASTER

		!Recv defects from other processors
		do i=1, myProc%numtasks-1
			allocate(defectsRecv(SPECIES+1,numDefectsRecv(i+1)))
			call MPI_RECV(defectsRecv,(SPECIES+1)*numDefectsRecv(i+1),MPI_INTEGER,i,400,comm,status,ierr)

			do j=1,numDefectsRecv(i+1)

				do k=1,SPECIES
					products(k)=defectsRecv(k,j)
				end do

				nullify(defectPrevList)
				defectCurrentList=>outputDefectList
				call findDefectInList(defectCurrentList, defectPrevList, products)

				!Update outputDefectList
				if(associated(defectCurrentList)) then !if we aren't at the end of the list
					same=0
					do k=1,SPECIES
						if(defectCurrentList%defectType(k)==products(k)) then
							same=same+1
						end if
					end do

					if(same==SPECIES) then
						!if the defect is already present in the list
						defectCurrentList%num=defectCurrentList%num+defectsRecv(SPECIES+1,j)
					else
						!if the defect is to be inserted in the list
						nullify(defectPrevList%next)
						allocate(defectPrevList%next)
						nullify(defectPrevList%next%next)
						defectPrevList=>defectPrevList%next
						allocate(defectPrevList%defectType(SPECIES))
						defectPrevList%cellNumber=0	!no need for cell numbers in outputDefectList
						defectPrevList%num=defectsRecv(SPECIES+1,j)

						do k=1,SPECIES
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
					allocate(defectPrevList%defectType(SPECIES))
					defectPrevList%cellNumber=0 !no need for cell numbers in outputDefectList
					defectPrevList%num=defectsRecv(SPECIES+1,j)

					do k=1,SPECIES
						defectPrevList%defectType(k)=products(k)
					end do
				else
					write(*,*) 'error tried to insert defect at beginning of output defect list'
				end if
			end do
			deallocate(defectsRecv)
		end do
	end if

	!Output totdat.out, defect.out, satdat.out
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
		numCu=0

		!total number of solute (Cu)/V/SIA clusters in the whole system
		numScluster=0
		numVoid=0
		numLoop=0
		numSVcluster=0
		numCuCluster=0

		!number density of clusters
		denScluster=0d0
		denVoid=0d0
		denLoop=0d0
		denSVcluster=0d0
		denCuCluster=0d0

		!average radius of clusters
		radiusScluster=0d0
		radiusVoid=0d0
		radiusLoop=0d0
		radiusSVcluster=0d0
		radiusCuCluster=0d0

		!average size of clusters
		sizeScluster=0d0
		sizeVoid=0d0
		sizeLoop=0d0
		sizeSVcluster=0d0
		sizeCuCluster=0d0

		VRetained=0d0
		VAnnihilated=0d0

		!<Output time information
		!DPA=dble(totalImpAnn(1))/(systemVol/(numDisplacedAtoms*atomSize))
		DPA = totalImpAnn(3)/(systemVol/atomSize)
		if(totdatToggle=='yes') then
			write(TOTFILE,*)
			write(TOTFILE,*)
			write(TOTFILE,*) '***********************************************************************'
			if(simStatus == 1) then		!<1: 'irradiation', 0: 'anneal'
				write(TOTFILE,*) 'Time', elapsedTime, 'Steps', step
				if(implantType=='FrenkelPair') then
					write(TOTFILE,*) 'FrenkelPairs', totalImpAnn(1), 'DPA', DPA
				else if(implantType=='Cascade')	then
					write(TOTFILE,*) 'Cascades', totalImpAnn(1), 'DPA', DPA
				else	!Thermal aging
					write(TOTFILE,*) 'noImplantation', totalImpAnn(1), 'DPA', DPA
				end if
			else if(simStatus == 0) then
				write(TOTFILE,*) 'Annealing process'
				write(TOTFILE,*) 'Anneal time', elapsedTime, 'Steps', step
			end if
			write(TOTFILE,*) 'Defects (Cu V I_m I_im[defectType] _ [number]):'
		end if
		if(defectToggle=='yes') then
			write(DEFFILE,*)
			write(DEFFILE,*)
			write(DEFFILE,*) '***********************************************************************'
			if(simStatus == 1) then		!<1: 'irradiation', 0: 'anneal'
				write(DEFFILE,*) 'Time', elapsedTime, 'Steps', step
				if(implantType=='FrenkelPair') then
					write(DEFFILE,*) 'FrenkelPairs', totalImpAnn(1), 'DPA', DPA
				else if(implantType=='Cascade')	then
					write(DEFFILE,*) 'Cascades', totalImpAnn(1), 'DPA', DPA
				else	!Thermal aging
					write(DEFFILE,*) 'noImplantation', totalImpAnn(1), 'DPA', DPA
				end if
			else if(simStatus == 0) then
				write(DEFFILE,*) 'Annealing process'
				write(DEFFILE,*) 'Anneal time', elapsedTime, 'Steps', step
			end if
			write(DEFFILE,*) 'Defects (Cu V I_m I_im[defectType] _ [number]):'
		end if
		if(stadatToggle=='yes') then
			write(STAFILE,*)
			write(STAFILE,*)
			write(STAFILE,*) '***********************************************************************'
			if(simStatus == 1) then		!<1: 'irradiation', 0: 'anneal'
				write(STAFILE,*) 'Time', elapsedTime, 'Steps', step
				if(implantType=='FrenkelPair') then
					write(STAFILE,*) 'FrenkelPairs', totalImpAnn(1), 'DPA', DPA
				else if(implantType=='Cascade')	then
					write(STAFILE,*) 'Cascades', totalImpAnn(1), 'DPA', DPA
				else	!Thermal aging
					write(STAFILE,*) 'noImplantation', totalImpAnn(1), 'DPA', DPA
				end if
			else if(simStatus == 0) then
				write(STAFILE,*) 'Annealing process'
				write(STAFILE,*) 'Anneal time', elapsedTime, 'Steps', step
			end if
		end if

		defectCurrentList=>outputDefectList%next
		do while(associated(defectCurrentList))

			if(defectCurrentList%defectType(1) /= 0) then	!Cu/CuV cluster

				if(defectCurrentList%defectType(1)==1 .AND. defectCurrentList%defectType(2)==0) then
					pointS=defectCurrentList%num
				end if

				if(defectCurrentList%defectType(2) == 0) then	!Cu cluster
					if(defectCurrentList%defectType(1) > minSCluster) then
						numS=numS+defectCurrentList%defectType(1)*defectCurrentList%num
						radiusScluster=radiusScluster+dble(defectCurrentList%num)*&
								(3*dble(defectCurrentList%defectType(1))*atomSize/(4*pi))**(1d0/3d0)
						numScluster=numScluster+defectCurrentList%num
					end if
				else	!Cu_Vac cluster
					if((defectCurrentList%defectType(1)+defectCurrentList%defectType(2)) > minSV) then
						numSV=numSV+(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*defectCurrentList%num
						radiusSVcluster = radiusSVcluster+dble(defectCurrentList%num)*&
								(3*dble(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*&
										atomSize/(4*pi))**(1d0/3d0)
						numSVcluster=numSVcluster+defectCurrentList%num
						totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
					end if

				end if

				if(defectCurrentList%defectType(1) > minSCluster) then
					numCu=numCu+defectCurrentList%defectType(1)*defectCurrentList%num
					radiusCuCluster=radiusCuCluster+dble(defectCurrentList%num)*&
							(3*dble(defectCurrentList%defectType(1))*atomSize/(4*pi))**(1d0/3d0)
					numCuCluster=numCuCluster+defectCurrentList%num
				end if

				!if(defectCurrentList%defectType(1) > minSCluster) then
				!	numS=numS+defectCurrentList%defectType(1)*defectCurrentList%num
				!	radiusScluster=radiusScluster+dble(defectCurrentList%num)*&
				!			(3*dble(defectCurrentList%defectType(1))*atomSize_Cu/(4*pi))**(1d0/3d0)
				!	numScluster=numScluster+defectCurrentList%num
				!end if
				!if(defectCurrentList%defectType(2)>0 .AND. &
				!		(defectCurrentList%defectType(1)+defectCurrentList%defectType(2)) > minSV) then
				!	numSV=numSV+(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*defectCurrentList%num
				!	radiusSVcluster = radiusSVcluster+dble(defectCurrentList%num)*&
				!			(3*dble(defectCurrentList%defectType(1)+defectCurrentList%defectType(2))*&
				!					atomSize_Cu/(4*pi))**(1d0/3d0)
				!	numSVcluster=numSVcluster+defectCurrentList%num
				!	if(defectCurrentList%defectType(2) /= 0) then
				!		totalVac=totalVac+defectCurrentList%defectType(2)*defectCurrentList%num
				!	end if
				!end if

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

			!<Output defects
			if(totdatToggle=='yes') then
				write(TOTFILE,*) defectCurrentList%defectType, defectCurrentList%num
			end if
			if(defectToggle=='yes') then
				write(DEFFILE,*) defectCurrentList%defectType, defectCurrentList%num
			end if

			defectCurrentList=>defectCurrentList%next
		end do

		denLoop = dble(numLoop)/systemVol
		denVoid = dble(numVoid)/systemVol
		denScluster = dble(numScluster)/systemVol
		denSVcluster = dble(numSVcluster)/systemVol
		denCuCluster = dble(numCuCluster)/systemVol

		!radiusLoop=((dble(numSIA)/dble(numLoop))*atomSize/(pi*burgers))**(1d0/2d0)
		!radiusVoid=(3d0*(dble(numV)/dble(numVoid))*atomSize/(4d0*pi))**(1d0/3d0)
		!radiusScluster=(3*(dble(numS)/dble(numScluster))*atomSize/(4*pi))**(1d0/3d0)
		!radiusSVcluster=(3d0*(dble(numSV)/dble(numSVcluster))*atomSize/(4d0*pi))**(1d0/3d0)
		radiusLoop=radiusLoop/dble(numLoop)
		radiusVoid=radiusVoid/dble(numVoid)
		radiusScluster=radiusScluster/dble(numScluster)
		radiusSVcluster=radiusSVcluster/dble(numSVcluster)
		radiusCuCluster=radiusCuCluster/dble(numCuCluster)

		sizeLoop=dble(numSIA)/dble(numLoop)
		sizeVoid=dble(numV)/dble(numVoid)
		sizeScluster=dble(numS)/dble(numScluster)
		sizeSVcluster=dble(numSV)/dble(numSVcluster)
		sizeCuCluster=dble(numCu)/dble(numCuCluster)

		conPointV=dble(pointV)/systemVol*atomSize
		conPointSIA=dble(pointSIA)/systemVol*atomSize

		VRetained = dble(totalVac)/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))
		VAnnihilated = dble(totalImpAnn(2))/(dble(numDisplacedAtoms)*dble(totalImpAnn(1)))

		!<Output totdat.out, stadat.out
		if(totdatToggle=='yes') then
			!<writa file
			write(TOTFILE,*) 'numCluster (Loop/Void/S/SV/precipitate/CuCluster):'
			write(TOTFILE,*) numLoop, numVoid, numScluster, numSVcluster, (numScluster+numSVcluster), numCuCluster
			write(TOTFILE,*) 'NumberDensity (m-3) (Loop/Void/S/SV/precipitate/CuCluster):'
			write(TOTFILE,*) denLoop*1d27, denVoid*1d27, denScluster*1d27,  denSVcluster*1d27, &
					(denScluster+denSVcluster)*1d27, denCuCluster*1d27
			write(TOTFILE,*) 'Concentration (Loop/Void/S/SV/precipitate/CuCluster):'
			write(TOTFILE,*) denLoop*atomSize, denVoid*atomSize, denScluster*atomSize, denSVcluster*atomSize, &
			(denScluster+denSVcluster)*atomSize, denCuCluster*atomSize
			write(TOTFILE,*) 'AverageRadius (nm) (Loop/Void/S/SV/precipitate/CuCluster):'
			write(TOTFILE,*) radiusLoop, radiusVoid, radiusScluster, radiusSVcluster, &
			(radiusScluster*dble(numScluster)+radiusSVcluster*dble(numSVcluster))/(dble(numScluster)+dble(numSVcluster)), &
					radiusCuCluster
			write(TOTFILE,*) 'AverageSize (Loop/Void/S/SV/precipitate/CuCluster):'
			write(TOTFILE,*) sizeLoop, sizeVoid, sizeScluster, sizeSVcluster, &
					(dble(numS)+dble(numSV))/(dble(numScluster)+dble(numSVcluster)), sizeCuCluster
			write(TOTFILE,*) 'ConcenPointDefects (V/SIA):', conPointV, conPointSIA
			write(TOTFILE,*) 'PercentVRetained',VRetained,'PercentVAnnihilated',VAnnihilated
		end if
		if(stadatToggle=='yes') then
			write(STAFILE,*) 'numCluster (Loop/Void/S/SV/precipitate):', &
					numLoop, numVoid, numScluster, numSVcluster, (numScluster+numSVcluster)
			write(STAFILE,*) 'NumberDensity (m-3) (Loop/Void/S/SV/precipitate):', &
					denLoop*1d27, denVoid*1d27, denScluster*1d27,  denSVcluster*1d27, (denScluster+denSVcluster)*1d27
			write(STAFILE,*) 'Concentration (Loop/Void/S/SV/precipitate):', &
					denLoop*atomSize, denVoid*atomSize, denScluster*atomSize_Cu, denSVcluster*atomSize_Cu, &
					(denScluster+denSVcluster)*atomSize_Cu
			write(STAFILE,*) 'AverageRadius (nm) (Loop/Void/S/SV/precipitate):', &
					radiusLoop, radiusVoid, radiusScluster, radiusSVcluster, &
					(radiusScluster*dble(numScluster)+radiusSVcluster*dble(numSVcluster))/(dble(numScluster)+dble(numSVcluster))
			write(STAFILE,*) 'AverageSize (Loop/Void/S/SV/precipitate):', &
					sizeLoop, sizeVoid, sizeScluster, sizeSVcluster, (dble(numS)+dble(numSV))/(dble(numScluster)+dble(numSVcluster))
			write(STAFILE,*) 'ConcenPointDefects (V/SIA):', conPointV, conPointSIA
			write(STAFILE,*) 'PercentVRetained',VRetained,'PercentVAnnihilated',VAnnihilated
		end if
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

