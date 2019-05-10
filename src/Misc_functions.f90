! $Header: /home/CVS//srscd/src/Misc_functions.f90,v 1.7 2015/04/07 14:36:22 aydunn Exp $
!***************************************************************************************************
!> Subroutine factorial
!! Factorial function returns n!
!! 
!! Inputs: n (integer)
!! Outputs: factorial (integer)
!!
!! This is limited to numbers n<17 for computer feasibility.
!***************************************************************************************************

integer function factorial(n)
implicit none

integer n, i, temp

if(n >= 17) then
	write(*,*) 'error factorial too large'
endif

temp=1
do i=1,n
	temp=temp*i
	!write(*,*) 'temp', temp
end do
factorial=temp

end function

!***************************************************************************************************
!> Subroutine binomial
!! Binomial function returns n choose k
!! 
!! Inputs: n, k (integers)
!! Outputs: Binomial (integer)
!***************************************************************************************************

integer function binomial(n,k)
implicit none

integer n, k, factorial

!write(*,*) 'binomial', n, k
!write(*,*) 'factorials', factorial(n), factorial(k), factorial(n-k)
!read(*,*) 
binomial=factorial(n)/(factorial(k)*factorial(n-k))

end function

!***************************************************************************************************
!
!> function TotalRateCheck() - checks if our total rate still matches the actual total rate
!!
!! This function is used as a diagnostic tool, it calcuates the total rate of all reactions in the 
!! system and compares it to totalRate, in order to test whether totalRate is being properly updated.
!!
!! Inputs: none
!! Output: total reaction rate (sum of all reaction rates)
!
!***************************************************************************************************

double precision function TotalRateCheck()
use mod_srscd_constants
use DerivedType
implicit none

integer cell, i, j, k
double precision rate, rateCell, rateCascade
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
	endif
	
	totalRateVol(i)=rateCell			!Also updating the total rate within each volume element
end do

!Compute total rate of all reactions in the fine meshes by going through active cascades
!2015.04.02: Need to create something analagous to totalRateVol() for the case of cascades.
!Put it in the cascade derived type definition
do while(associated(CascadeCurrent))

	do i=1,numCellsCascade
		rateCascade=0d0
	
		reactionCurrent=>CascadeCurrent%reactionList(i)
		
		do while(associated(reactionCurrent))
			rate=rate+reactionCurrent%reactionRate
			rateCascade=rateCascade+reactionCurrent%reactionRate			
			reactionCurrent=>reactionCurrent%next
		end do
		
		if(dabs(CascadeCurrent%totalRate(i)-rateCascade) .GT. 1d0) then
			write(*,*) 'Error: total rate differs significantly from expected in cascade', &
				CascadeCurrent%cascadeID
			write(*,*) 'TotalRateCascade', CascadeCurrent%totalRate(i), 'Expected Value', rateCascade
		endif
	
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
!	endif
!endif

TotalRateCheck=rate

end function

!***************************************************************************************************
!> Function total rate cascade - finds total reaction rate within cascade
!!
!! Input: CascadeCurrent%totalRate(:), list of reaction rates within each volume element in the cascade
!! 
!! Output: sum of CascadeCurrent%totalRate(:), total reaction rate for entire cascade
!!
!! This subroutine was added on 2015.04.03 when I switched to tracking the reaction rate in each volume element
!! instead of the entire cascade, for the sake of computational efficiency when choosing reactions.
!***************************************************************************************************

double precision function totalRateCascade(CascadeCurrent)
use DerivedType
use mod_srscd_constants
implicit none

type(cascade), pointer :: CascadeCurrent
integer i
double precision rateSum

rateSum=0

do i=1,numCellsCascade

	rateSum=rateSum+CascadeCurrent%totalRate(i)
	
end do

totalRateCascade=rateSum

end function


