!************************************************
!2019.05.29 Add
!This function is used to compute the vacancy concentration at this time
!************************************************

subroutine computeVconcent()
    use DerivedType
    use mod_srscd_constants
    implicit none
    include 'mpif.h'

    integer i,j
    integer numVacancy, totalVacancy

    type(defect), pointer :: defectCurrent
    type(cascade), pointer :: CascadeCurrent

    nullify(defectCurrent)
    nullify(CascadeCurrent)
    totalVacancy=0d0
    numVacancy = 0d0

    !count the vacancies in fine meshes
    if(totalDPA > 0d0 .AND. DPARate > 0d0) then
        if(implantType=='Cascade' .AND. meshingType=='adaptive') then
            CascadeCurrent=>ActiveCascades

            outer1: do while(associated(CascadeCurrent))

                do j=1, numCellsCascade

                    defectCurrent=>CascadeCurrent%localDefects(j)
                    inter1: do while(associated(defectCurrent))
                        if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2)==1 &
                                .AND. defectCurrent%defectType(3)==0 .AND. defectCurrent%defectType(4)==0) then
                            numVacancy=numVacancy+defectCurrent%num
                            exit inter1
                        else
                            defectCurrent=>defectCurrent%next
                        end if

                    end do inter1

                end do

                CascadeCurrent=>CascadeCurrent%next
            end do outer1

        end if

    end if

    !count the vacancies in Coarse meshes
    outer: do i=1, numCells
        defectCurrent=>defectList(i)%next

        inter: do while(associated(defectCurrent))

            if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2)==1 &
                    .AND. defectCurrent%defectType(3)==0 .AND. defectCurrent%defectType(4)==0) then
                numVacancy=numVacancy+defectCurrent%num
                exit inter
            else
                defectCurrent=>defectCurrent%next
            end if
        end do inter

    end do outer

    call MPI_ALLREDUCE(numVacancy,totalVacancy,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(totalVacancy > 0) then
        Vconcent = totalVacancy/((myProc%globalCoord(2)-myProc%globalCoord(1))/lattice * &
                (myProc%globalCoord(4)-myProc%globalCoord(3))/lattice * &
                (myProc%globalCoord(6)-myProc%globalCoord(5))/lattice * 2)
    else
        Vconcent = initialCeqv
    end if



end subroutine

!************************************************
!This function is used to compute the vacancy concentration at this time
!NOTE:  Adaptive meshes are not considered.
!************************************************

subroutine computeSIAconcent()
    use DerivedType
    use mod_srscd_constants
    implicit none
    include 'mpif.h'

    integer i,j
    integer numSIA, totalSIA

    type(defect), pointer :: defectCurrent
    type(cascade), pointer :: CascadeCurrent

    nullify(defectCurrent)
    nullify(CascadeCurrent)
    totalSIA=0d0
    numSIA = 0d0

    !count the SIAs in fine meshes
    if(totalDPA > 0d0 .AND. DPARate > 0d0) then
        if(implantType=='Cascade' .AND. meshingType=='adaptive') then
            CascadeCurrent=>ActiveCascades

            outer1: do while(associated(CascadeCurrent))

                do j=1, numCellsCascade

                    defectCurrent=>CascadeCurrent%localDefects(j)
                    inter1: do while(associated(defectCurrent))
                        if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2)==0 &
                                .AND. defectCurrent%defectType(3)==1 .AND. defectCurrent%defectType(4)==0) then
                            numSIA=numSIA+defectCurrent%num
                            exit inter1
                        else
                            defectCurrent=>defectCurrent%next
                        end if

                    end do inter1

                end do

                CascadeCurrent=>CascadeCurrent%next
            end do outer1

        end if

    end if

    !count the SIAs in coarse meshes
    outer: do i=1, numCells
        defectCurrent=>defectList(i)%next

        inter: do while(associated(defectCurrent))

            if(defectCurrent%defectType(1)==0 .AND. defectCurrent%defectType(2)==0 &
                    .AND. defectCurrent%defectType(3)==1 .AND. defectCurrent%defectType(4)==0) then
                numSIA=numSIA+defectCurrent%num
                exit inter
            else
                defectCurrent=>defectCurrent%next
            end if
        end do inter

    end do outer

    call MPI_ALLREDUCE(numSIA,totalSIA,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(totalSIA > 0) then
        SIAconcent = totalSIA/((myProc%globalCoord(2)-myProc%globalCoord(1))/lattice * &
                (myProc%globalCoord(4)-myProc%globalCoord(3))/lattice * &
                (myProc%globalCoord(6)-myProc%globalCoord(5))/lattice * 2)
    else
        SIAconcent = initialCeqi
    end if



end subroutine