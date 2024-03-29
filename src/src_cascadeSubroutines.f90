!***************************************************************************************************
!>Subroutine: chooseCascade(CascadeTemp): chooses one cascade randomly in cascadeList
!!Inputs: CascadeList (global variable)
!!Output: CascadeTemp (pointing at the cascade we want)
!***************************************************************************************************
subroutine chooseCascade(CascadeTemp)
    use mod_structures
    use mod_randdp
    use mod_globalVariables
    implicit none

    type(cascadeEvent), pointer, intent(inout) :: cascadeTemp
    double precision :: r, atemp
    double precision, external :: damageFunc

    r=dprand()
    atemp=0d0
    cascadeTemp=>CascadeList

    do while(associated(cascadeTemp))
        atemp=atemp+1d0/dble(numCascades)
        if(atemp >= r) then
            numImpAnn(3) = numImpAnn(3) + damageFunc(PKAenergy)
            exit
        end if
        cascadeTemp=>cascadeTemp%next
    end do
end subroutine

!****************************************************************************************
!>Subroutine: chooseCascade_files(CascadeTemp): chooses one cascade randomly from a number of cascade files
!Inputs: CascadeList (global variable)
!Output: CascadeTemp (pointing at the cascade we want)
!*****************************************************************************************
subroutine chooseCascade_withFiles(cascadeTemp)
    use mod_structures
    use mod_randdp
    use mod_globalVariables
    implicit none

    type(cascadeEvent), pointer, intent(inout) :: cascadeTemp
    double precision :: r, r1, atemp, atemp1, energy, minEnergy
    integer :: i
    double precision, external :: sample_PKA_energy, damageFunc

    if(PKAspectrum == 'yes') then
        energy = sample_PKA_energy()      !<energy > 0d0

        minEnergy=1.0d8
        outer1: do i=1, numCascadeFiles
            !if(abs(energy-cascadeLists(i)%PKAenergy) <= 5d0) then   !find the cascade file
            if(abs(energy-cascadeLists(i)%PKAenergy) < minEnergy) then   !find the cascade file

                r=dprand()
                atemp=0d0
                cascadeTemp=>cascadeLists(i)%listCascades
                do while(associated(cascadeTemp))
                    atemp=atemp+1d0/dble(cascadeLists(i)%numCascades)
                    if(r<=atemp) then
                        numImpAnn(3) = numImpAnn(3) + damageFunc(cascadeLists(i)%PKAenergy)
                        exit outer1
                    else
                        cascadeTemp=>cascadeTemp%next
                    end if
                end do
            end if
        end do outer1
    else    !not use PKA spectrum
        r1 = dprand()   !<used to choose a cascade file
        atemp1 = 0d0
        outer2: do i=1, numCascadeFiles
            atemp1 =atemp1 +1d0/dble(numCascadeFiles)
            if(r1 <= atemp1) then   !<chosed a cascade file
                r=dprand()
                atemp=0d0
                cascadeTemp=>cascadeLists(i)%listCascades
                do while(associated(cascadeTemp))
                    atemp=atemp+1d0/dble(cascadeLists(i)%numCascades)
                    if(r<=atemp) then   !chosed a cascade from these cascade files
                        numImpAnn(3) = numImpAnn(3) + damageFunc(cascadeLists(i)%PKAenergy)
                        exit outer2
                    else
                        cascadeTemp=>cascadeTemp%next
                    end if
                end do
            end if
        end do outer2
    end if

end subroutine chooseCascade_withFiles

!***************************************************************************************************
!> function sample_PKA_energy(): The PKA energy is sampled by linear interpolation
!***************************************************************************************************
double precision function sample_PKA_energy()
    use mod_structures
    use mod_randdp
    use mod_globalVariables
    implicit none

    double precision :: cpdf, r
    integer :: i, index_PKA
    double precision :: slope

    index_PKA=0
    r=dprand()
    do i=2, EPKAlist%size
        if(r > EPKAlist%cpdf(i-1) .AND. r < EPKAlist%cpdf(i)) then
            index_PKA=i
            exit
        end if
    end do

    slope = (EPKAlist%energy(index_PKA) - EPKAlist%energy(index_PKA-1))/(EPKAlist%cpdf(index_PKA) - EPKAlist%cpdf(index_PKA-1))
    sample_PKA_energy = slope*(r - EPKAlist%cpdf(index_PKA-1)) + EPKAlist%energy(index_PKA-1)

end function

!***************************************************************************************************
!> function damageFunc(Epka): Calculate the number of displaced atoms according to the PKA energy
!***************************************************************************************************
double precision function damageFunc(Epka)
    use mod_constants
    implicit none

    double precision, intent(in) :: Epka
    double precision :: vNRT

    vNRT = 0d0
    if(Epka <= ED) then
        vNRT = 0d0
    else if(Epka <= (2d0*ED/0.8d0)) then
        vNRT = 1d0
    else
        vNRT = 0.8d0*Epka/(2d0*ED)
    end if
    damageFunc = vNRT

end function

!***************************************************************************************************
!> Integer function CascadeCount(): This subroutine counts how many cascades are active in the LOCAL mesh (not on all processors)
! and returns that value
! Inputs: none
! Outputs: number of cascades present
!***************************************************************************************************
integer function CascadeCount()
    use mod_globalVariables
    use mod_structures
    implicit none

    type(cascade), pointer :: CascadeCurrent
    integer :: count

    CascadeCurrent=>ActiveCascades
    count=0

    do while(associated(CascadeCurrent))
        count=count+1
        CascadeCurrent=>CascadeCurrent%next
    end do
    CascadeCount=count
end function

!***************************************************************************************************
!>subroutine addCascadeExplicit(reactionCurrent, CascadeCurrent): This subroutine takes the place of chooseReaction
!in the case of explicit cascade implantation. It forces the program to 'choose' a cascade reaction,
!instead of using the Monte Carlo algorithm to do so.
! Inputs: none
! Outputs: reactionCurrent and CascadeCurrent, pointing at o-th reaction.
!***************************************************************************************************
subroutine addCascadeExplicit(reactionCurrent, CascadeCurrent)
    use mod_structures
    use mod_globalVariables
    use mod_randdp
    implicit none

    type(reaction), pointer, intent(inout) :: reactionCurrent
    type(cascade), pointer, intent(inout) :: cascadeCurrent
    double precision :: r2, atemp, r2timesa
    integer :: i

    atemp=0d0
    r2=dprand()
    r2timesa=r2*numCells
    nullify(reactionCurrent)
    nullify(CascadeCurrent)

    outer: do i=1,numCells
        reactionCurrent=>reactionList(i)
        atemp=atemp+1d0		!Here we don't have a reaction rate, and each cell is weighted evenly (assuming uniform mesh)
        if(atemp >= r2timesa) then
            exit outer		!exit both loops with reactionCurrent pointing to the randomly chosen reaction
        end if
    end do outer

    !Checking that reactionCurrent is pointed at the correct reaction
    if(implantType=='Cascade') then
        if(reactionCurrent%numReactants==-10 .AND. reactionCurrent%numProducts==0) then !Cascade implantation
            numImpAnn(1)=numImpAnn(1)+1
        end if
    else
        write(*,*) 'Error wrong implant type for explicit procedure'
    end if

end subroutine

!***************************************************************************************************
!>logical function cascadeMixingCheck(): cascadeMixingCheck checks whether a given defect in the
!fine mesh interacts with the cascade. It returns a logical value of true or false.
! Input: none
! Output: logical value for whether or not a fine mesh defect combines with a cascade.
!***************************************************************************************************
logical function cascadeMixingCheck()
    use mod_globalVariables
    use mod_randdp
    implicit none

    double precision :: r1, probability
    logical :: boolean

    probability=cascadeVolume/(numCellsCascade*CascadeElementVol)

    r1=dprand()
    if(r1 > probability) then
        boolean=.FALSE.
    else
        boolean=.TRUE.
    endif
    cascadeMixingCheck=boolean

end function

!***************************************************************************************************
!> Subroutine cascadeUpdateStep(releaseToggle, cascadeCell): Carries out the communication necessary
!when a cascade is created or destroyed in an element that bounds another processor.
!Step 1: Check to see if cascade was created/destroyed in element that is in the boundary of another processor.
!Step 2: Send message to neighboring processors about whether a cascade was created/destroyed in their
!boundaries, and if so the number of defects to recieve in the boundary element
!Step 3: Send/recieve boundary element defects (just completely re-write boundary element defects in this step)
!Step 4: Update reaction rates for all diffusion reactions from elements neighboring the boundary elements which have been updated
!Inputs: cascadeCell (integer, 0 if no cascade, otherwise gives number of volume element that cascade event has occurred in)
!Outputs: none
!***************************************************************************************************
subroutine cascadeUpdateStep(releaseToggle, cascadeCell)
    use mod_constants
    use mod_structures
    use mod_globalVariables
    use mod_updatereactions
    implicit none
    include 'mpif.h'

    logical, intent(in) :: releaseToggle
    integer, intent(in) :: cascadeCell
    type(defect), pointer :: defectCurrent, defectPrev
    integer :: j, dir, count, recvDir
    integer :: cellNumber, bndryCellNumber, localGrainID, neighborGrainID
    !Used for communication between processors
    integer :: numSend,numSendTemp, numRecv
    integer, allocatable :: defectSend(:,:), defectRecv(:,:)
    integer :: status(MPI_STATUS_SIZE), sendStatus(MPI_STATUS_SIZE),recvStatus(MPI_STATUS_SIZE)
    integer :: sendRequest, recvRequest
    double precision :: commTime1, commTime2, compTime1, compTime2

    !Fill send buffer
    if(cascadeCell==0) then
        numSend=0
        allocate(defectSend(SPECIES+1,numSend))
    else    !cascadeCell /= 0: Implant a new cascade or release an existing  cascade
        numSend=0
        defectCurrent=>defectList(cascadeCell)%next

        do dir=1,6
            if(myMesh(cascadeCell)%neighborProcs(dir) /= myProc%taskid .AND. &
                    myMesh(cascadeCell)%neighborProcs(dir) /= -1) then

                do while(associated(defectCurrent))
                    numSend=numSend+1
                    defectCurrent=>defectCurrent%next
                end do
                exit
            end if
        end do

        if(numSend /=0 ) then
            numSend=numSend+1
            allocate(defectSend(SPECIES+1,numSend))

            defectSend(1,1)=0		        !myMesh(cascadeCell)%neighbors(1,dir)
            defectSend(2,1)=numSend-1		!numDefects
            defectSend(3,1)=cascadeCell	    !cascadeCell
            defectSend(4,1)=0		        !volume +/- CascadeElementVol*numCellsCascade
            defectSend(5,1)=0		        !useless

            if(releaseToggle .EQV. .TRUE.) then
                defectSend(4,1)=1		        !volume + CascadeElementVol*numCellsCascade
            else
                defectSend(4,1)=-1		        !volume - CascadeElementVol*numCellsCascade
            end if

            defectCurrent=>defectList(cascadeCell)%next

            do j=2,numSend
                defectSend(1:SPECIES,j)=defectCurrent%defectType(:)
                defectSend(SPECIES+1,j)=defectCurrent%num

                defectCurrent=>defectCurrent%next
            end do
        end if
    end if

    do dir=1,6

        if(mod(dir,2)==0) then
            recvDir=dir-1
        else
            recvDir=dir+1
        end if
        numSendTemp = numSend
        if(cascadeCell /= 0) then
            if(myMesh(cascadeCell)%neighborProcs(dir) /= myProc%taskid .AND. &
                    myMesh(cascadeCell)%neighborProcs(dir) /= -1 .AND. numSend /=0)  then
                defectSend(1,1)=myMesh(cascadeCell)%neighbors(dir)
            else
                numSendTemp = 0
            end if
        end if

        !Send
        commTime1 = MPI_WTIME()
        !if(myProc%procNeighbor(dir) /= myProc%taskid) then
            call MPI_ISEND(defectSend, (SPECIES+1)*numSendTemp, MPI_INTEGER, myProc%procNeighbor(dir), &
                    900+dir, comm, sendRequest, ierr)
        !end if

        !if(myProc%procNeighbor(recvDir) /= myProc%taskid) then

            count=0
            call MPI_PROBE(myProc%procNeighbor(recvDir), 900+dir, comm,status,ierr)
            call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)

            numRecv=count/(SPECIES+1)
            allocate(defectRecv(SPECIES+1,numRecv))

            call MPI_RECV(defectRecv,(SPECIES+1)*numRecv,MPI_INTEGER,myProc%procNeighbor(recvDir),&
                    900+dir,comm,status,ierr)

        compTime1 = MPI_WTIME() !标记计算时间的起始位置
            if(numRecv /= 0) then
                !Update my boundary
                bndryCellNumber=defectRecv(3,1)
                cellNumber=defectRecv(1,1)

                if(meshingType=='adaptive') then
                    if(defectRecv(4,1)==-1) then
                        myGhost(bndryCellNumber,recvDir)%volume=myGhost(bndryCellNumber,recvDir)%volume-&
                                CascadeElementVol*dble(numCellsCascade)
                        myGhost(bndryCellNumber,recvDir)%length=(myGhost(bndryCellNumber,recvDir)%volume)&
                                **(1d0/3d0)
                    else if(defectRecv(4,1)==1) then
                        myGhost(bndryCellNumber,recvDir)%volume=myGhost(bndryCellNumber,recvDir)%volume+&
                                CascadeElementVol*dble(numCellsCascade)
                        myGhost(bndryCellNumber,recvDir)%length=(myGhost(bndryCellNumber,recvDir)%volume)&
                                **(1d0/3d0)
                    end if
                end if


                !remove defects from myGhost (except for first defect, this is all 0's and is just a placeholder)
                defectCurrent=>myGhost(bndryCellNumber,recvDir)%defectList%next

                !delete exiting defects
                nullify(defectPrev)
                do while(associated(defectCurrent))
                    defectPrev=>defectCurrent
                    defectCurrent=>defectCurrent%next
                    deallocate(defectPrev%defectType)
                    deallocate(defectPrev)
                end do

                !add defects to my boundary
                defectCurrent=>myGhost(bndryCellNumber,recvDir)%defectList
                do j=2,numRecv
                    nullify(defectCurrent%next)
                    allocate(defectCurrent%next)
                    nullify(defectCurrent%next%next)
                    defectCurrent=>defectCurrent%next
                    allocate(defectCurrent%defectType(SPECIES))
                    defectCurrent%cellNumber=bndryCellNumber
                    defectCurrent%num=defectRecv(SPECIES+1,j)
                    defectCurrent%defectType(:)=defectRecv(1:SPECIES,j)
                end do

                !*******************
                !Add Diffusion reactions
                !*******************
                defectCurrent=>defectList(cellNumber)
                do while(associated(defectCurrent))
                    if (myMesh(cellNumber)%numNeighbors(recvDir)==0) then
                        write(*,*) 'error myMesh does not have neighbors in this direction'
                    end if

                    if(polycrystal=='yes') then

                        localGrainID=myMesh(cellNumber)%material
                        if(myProc%procNeighbor(recvDir)/=myProc%taskid .AND. myProc%procNeighbor(recvDir)/=-1) then
                            neighborGrainID=myGhost(myMesh(cellNumber)%neighbors(recvDir),recvDir)%material
                        else
                            neighborGrainID=myMesh(myMesh(cellNumber)%neighbors(recvDir))%material
                        end if

                        if(localGrainID==neighborGrainID) then
                            !Allow diffusion between elements in the same grain
                            call update_diff_reactions(cellNumber, bndryCellNumber,&
                                    myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                        else
                            !Assume perfect sinks at grain boundaries - treat grain boundaries like free surfaces for now
                            call update_diff_reactions(cellNumber,0,myProc%taskid,-1,recvDir,defectCurrent%defectType)
                        end if
                    else
                        !Add diffusion reactions from this cell to neighboring cells
                        call update_diff_reactions(cellNumber, bndryCellNumber,&
                                myProc%taskid, myProc%procNeighbor(recvDir),recvDir,defectCurrent%defectType)
                    end if
                    defectCurrent=>defectCurrent%next
                end do
            end if
            deallocate(defectRecv)
        !end if
        compTime2 = MPI_WTIME() !标记计算时间的结束位置
        if(myProc%procNeighbor(dir) /= myProc%taskid) then
            call MPI_WAIT(sendRequest, sendStatus, ierr)
        end if
        commTime2 = MPI_WTIME()
        commTimeSum = commTimeSum + (commTime2-commTime1)-(compTime2-compTime1)
        casCommTime = casCommTime + (commTime2-commTime1)-(compTime2-compTime1)
    end do

end subroutine

!***************************************************************************************************
!> subroutine createCascadeConnect(): This subroutine assigns values to the connectivity matrix
!(global variable) used for all cascades
!Input: numxcascade, numycascade, nunmzcascade (global variables) : from configure.in
!Output: cascadeConnectivity (global variable)
!***************************************************************************************************
subroutine  createCascadeConnect()
    use mod_globalVariables
    implicit none

    integer :: cell

    do cell=1,numCellsCascade
        if(mod(cell,numxcascade)==0) then !identify cell to the right
            cascadeConnectivity(1,cell)=cell-numxcascade+1
            !cascadeConnectivity(1, cell)=0	!free in +x
        else
            cascadeConnectivity(1, cell)=cell+1
        end if

        if(mod(cell+numxcascade-1,numxcascade)==0) then !identify cell to the left
            cascadeConnectivity(2,cell)=cell+numxcascade-1
            !cascadeConnectivity(2, cell)=0	!free in -x
        else
            cascadeConnectivity(2, cell)=cell-1
        end if

        if(mod(cell,numxcascade*numycascade) > numxcascade*(numycascade-1) .OR. &
                mod(cell,numxcascade*numycascade)==0) then
            !cascadeConnectivity(3, cell)=0	!free in +y
            cascadeConnectivity(3,cell)=cell-(numxcascade*(numycascade-1))
        else
            cascadeConnectivity(3, cell)=cell+numxcascade
        end if

        if(mod(cell,numxcascade*numycascade) <= numxcascade .AND. mod(cell, numxcascade*numycascade) /= 0) then
            !cascadeConnectivity(4, cell)=0	!free in -y
            cascadeConnectivity(4,cell)=cell+(numxcascade*(numycascade-1))
        else
            cascadeConnectivity(4, cell)=cell-numxcascade
        end if

        if(mod(cell,numxcascade*numycascade*numzcascade) > numxcascade*numycascade*(numzcascade-1) .OR. &
                mod(cell, numxcascade*numycascade*numzcascade)==0) then
            !cascadeConnectivity(5, cell)=0	!free in +z
            cascadeConnectivity(5, cell)=cell-(numxcascade*numycascade*(numzcascade-1))
        else
            cascadeConnectivity(5, cell)=cell+numxcascade*numycascade
        end if

        if(mod(cell,numxcascade*numycascade*numzcascade) <= numxcascade*numycascade .AND. &
                mod(cell,numxcascade*numycascade*numzcascade) /= 0) then
            !cascadeConnectivity(6, cell)=0	!free in -z
            cascadeConnectivity(6, cell)=cell+(numxcascade*numycascade*(numzcascade-1))
        else
            cascadeConnectivity(6, cell)=cell-numxcascade*numycascade
        end if
    end do

end subroutine
