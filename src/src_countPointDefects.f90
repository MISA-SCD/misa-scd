!************************************************
!2019.05.29 Add
!This function is used to compute vacancy concentration。
!************************************************
subroutine computeVconcent()
    use mod_structures
    use mod_globalVariables
    implicit none
    include 'mpif.h'

    integer :: i,j, numVacancy, totalVacancy
    type(defect), pointer :: defectCurrent
    type(cascade), pointer :: CascadeCurrent

    nullify(defectCurrent)
    nullify(CascadeCurrent)
    totalVacancy=0d0
    numVacancy = 0d0

    !count the vacancies in fine meshes
    if(totalDPA > 0d0 .AND. dpaRate > 0d0) then
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

    call MPI_ALLREDUCE(numVacancy,totalVacancy,1,MPI_INTEGER,MPI_SUM,comm,ierr)

    if(totalVacancy > 0) then
        concV = totalVacancy/((myProc%globalCoord(2)-myProc%globalCoord(1))/lattice * &
                (myProc%globalCoord(4)-myProc%globalCoord(3))/lattice * &
                (myProc%globalCoord(6)-myProc%globalCoord(5))/lattice * 2)
    else
        concV = ceqV
    end if

end subroutine

!************************************************
!This function is used to compute the vacancy concentration at this time
!NOTE:  Adaptive meshes are not considered.
!************************************************
subroutine computeSIAconcent()
    use mod_structures
    use mod_globalVariables
    implicit none
    include 'mpif.h'

    integer :: i,j, numSIA, totalSIA
    type(defect), pointer :: defectCurrent
    type(cascade), pointer :: CascadeCurrent

    nullify(defectCurrent)
    nullify(CascadeCurrent)
    totalSIA=0d0
    numSIA = 0d0

    !count the SIAs in fine meshes
    if(totalDPA > 0d0 .AND. dpaRate > 0d0) then
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

    call MPI_ALLREDUCE(numSIA,totalSIA,1,MPI_INTEGER,MPI_SUM,comm,ierr)

    if(totalSIA > 0) then
        concI = totalSIA/((myProc%globalCoord(2)-myProc%globalCoord(1))/lattice * &
                (myProc%globalCoord(4)-myProc%globalCoord(3))/lattice * &
                (myProc%globalCoord(6)-myProc%globalCoord(5))/lattice * 2)
    else
        concI = ceqI
    end if

end subroutine

!**********************************************************************************
!This function is used to compute the vacancy concentration at this time
!This function will not be used
!**********************************************************************************
double precision function permanentCv()
    use mod_constants
    use mod_structures
    use mod_globalVariables
    implicit none

    double precision :: Kiv, diffV, diffI
    integer :: i

    do i=1,numSingleDiff
        if(DiffSingle(i)%defectType(1)==0 .AND. DiffSingle(i)%defectType(2)==1 .AND. &
                DiffSingle(i)%defectType(3)==0 .AND. DiffSingle(i)%defectType(4)==0) then   !0_V_0_0
            diffV = DiffSingle(i)%D*dexp(-DiffSingle(i)%Em/(kboltzmann*temperature))
        else if(DiffSingle(i)%defectType(1)==0 .AND. DiffSingle(i)%defectType(2)==0 .AND. &
                DiffSingle(i)%defectType(3)==1 .AND. DiffSingle(i)%defectType(4)==0) then   !0_0_I_0
            diffI = DiffSingle(i)%D*dexp(-DiffSingle(i)%Em/(kboltzmann*temperature))
        end if
    end do

    Kiv = 4*pi/atomSize*reactionRadius*(diffV + diffI)
    permanentCv = -dislocationDensity*Zint*diffI/(2*Kiv)+&
            ((dislocationDensity*Zint*diffI/(2*Kiv))**(2d0)+dpaRate*Zint*diffI/(Kiv*diffV))**(1d0/2d0)

end function