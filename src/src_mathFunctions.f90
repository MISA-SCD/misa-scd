!***************************************************************************************************
!> Subroutine factorial
!Factorial function returns n!
!This is limited to numbers n<17 for computer feasibility.
!***************************************************************************************************
integer function factorial(n)
	implicit none

	integer, intent(in) :: n
	integer :: i, temp

	if(n >= 17) then
		write(*,*) 'error factorial too large'
	endif

	temp=1
	do i=1,n
		temp=temp*i
	end do
	factorial=temp
end function

!***************************************************************************************************
!> Subroutine binomial
!Binomial function returns n choose k
!***************************************************************************************************
integer function binomial(n,k)
	implicit none

	integer, intent(in) :: n, k
	integer, external :: factorial

	binomial=factorial(n)/(factorial(k)*factorial(n-k))

end function

!***************************************************************************************************
!> function TotalRateCheck() - checks if our total rate still matches the actual total rate
!This function is used as a diagnostic tool, it calcuates the total rate of all reactions in the
!system and compares it to totalRate, in order to test whether totalRate is being properly updated.
!***************************************************************************************************
double precision function TotalRateCheck()
	use mod_globalVariables
	use mod_structures
	implicit none

	integer :: i
	double precision :: rate, rateCell, rateCascade
	type(reaction), pointer :: reactionCurrent
	type(cascade), pointer :: cascadeCurrent

	CascadeCurrent=>ActiveCascades
	rate=0d0

	!Compute total rate of all reactions in the coarse mesh
	do i=1,numCells
		rateCell=0d0
		reactionCurrent=>reactionList(i)

		do while(associated(reactionCurrent))
			rate=rate+reactionCurrent%reactionRate
			rateCell=rateCell+reactionCurrent%reactionRate
			reactionCurrent=>reactionCurrent%next
		end do

		if(dabs(totalRateVol(i)-rateCell) > 1d0) then
			write(*,*) 'Error: total rate differs significantly from expected in cell', i
			write(*,*) 'TotalRateVol', totalRateVol(i), 'expected value', rateCell
			call MPI_ABORT(comm,ierr)
		end if
		totalRateVol(i)=rateCell			!Also updating the total rate within each volume element
	end do

	!Compute total rate of all reactions in the fine meshes by going through active cascades
	do while(associated(CascadeCurrent))

		do i=1,numCellsCascade
			rateCascade=0d0
			reactionCurrent=>CascadeCurrent%reactionList(i)
			do while(associated(reactionCurrent))
				rate=rate+reactionCurrent%reactionRate
				rateCascade=rateCascade+reactionCurrent%reactionRate
				reactionCurrent=>reactionCurrent%next
			end do

			if(dabs(CascadeCurrent%totalRate(i)-rateCascade) > 1d0) then
				write(*,*) 'Error: total rate differs significantly from expected in cascade', &
						CascadeCurrent%cascadeID
				write(*,*) 'TotalRateCascade', CascadeCurrent%totalRate(i), 'Expected Value', rateCascade
			end if
			CascadeCurrent%totalRate(i)=rateCascade
		end do
		CascadeCurrent=>CascadeCurrent%next
	end do

	!Compare rate to totalRate and return error if different
	!if(myProc%taskid==MASTER) then
	!	if (rate==totalRate) then
	!		write(*,*) 'TotalRate correct', totalRate
	!	else
	!		write(*,*) 'TotalRate incorrect ', 'totalRate=', totalRate, 'actual total rate=', rate
	!	end if
	!end if

	TotalRateCheck=rate
end function

!***************************************************************************************************
!> Function total rate cascade - finds total reaction rate within cascade
!***************************************************************************************************
double precision function totalRateCascade(CascadeCurrent)
	use mod_structures
	use mod_globalVariables
	implicit none

	type(cascade), pointer, intent(in) :: CascadeCurrent
	integer :: i
	double precision :: rateSum

	rateSum=0d0
	do i=1,numCellsCascade
		rateSum=rateSum+CascadeCurrent%totalRate(i)
	end do
	totalRateCascade=rateSum

end function


