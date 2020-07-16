!*****************************************************************************************
!>double precision function find diffusivity - returns the diffusivity of a given defect type
!!
!! This function looks up the diffusivity of a given defect type using input date from 
!! material parameter input file. It either looks up diffusivity from values in a list or
!! computes diffusivity using diffusivityCompute in the case of a functional form.
!!
!!Input: defect type
!!Output: diffusivity (nm^2/s)
!*****************************************************************************************
double precision function findDiffusivity(matNum, defectType)
use DerivedType
use mod_constants
implicit none

integer defectType(numSpecies)
integer i, j, numSame, matNum
double precision Diff
double precision DiffusivityCompute, diffusivityCu

!Temporary: used as a parameter to vary the diffusivity of all defects on GB
double precision, parameter :: Param=0d0

!***************************************************************************************************
!This function returns the diffusivity of defect DefectType.
!It searches to see if DefectType is listed in DiffSingle, and if not, if it is listed in DiffFunction.
!If it is in neither, it outputs an error message (defect type should not exist)
!***************************************************************************************************

outer1: do i=1,numSingleDiff(matNum)
	numSame=0
	do j=1,numSpecies
		if(defectType(j)==DiffSingle(i,matNum)%defectType(j)) then
			numSame=numSame+1
		end if
	end do

	if (numSame==numSpecies) then
		if(matNum==2) then
			if(DefectType(1)==1 .AND. DefectType(2)==0 .AND. DefectType(3)==0 .AND. DefectType(4)==0) then
				Diff=diffusivityCu(matNum)
				exit outer1
			else
				Diff=DiffSingle(i,matNum)%D*dexp(-(DiffSingle(i,matNum)%Em-Param)/(kboltzmann*temperature))
				exit outer1
			end if

		else
			if(DefectType(1)==1 .AND. DefectType(2)==0 .AND. DefectType(3)==0 .AND. DefectType(4)==0) then

				Diff=diffusivityCu(matNum)
				exit  outer1
			else
				Diff=DiffSingle(i,matNum)%D*dexp(-DiffSingle(i,matNum)%Em/(kboltzmann*temperature))
				exit outer1
			end if
		
		end if
	end if
end do outer1

if(i==numSingleDiff(matNum)+1) then	!did not find defect in single defect list
	do i=1,numFuncDiff(matNum)
		numSame=0
		do j=1,numSpecies
			if(defectType(j)==0 .AND. DiffFunc(i,matNum)%defectType(j)==0) then
				numSame=numSame+1
			else if(defectType(j) /= 0 .AND. DiffFunc(i,matNum)%defectType(j)==1) then
				if(defectType(j) >= DiffFunc(i,matNum)%min(j) .AND. &
						(defectType(j) <= DiffFunc(i,matNum)%max(j) .OR. DiffFunc(i,matNum)%max(j)==-1)) then
					numSame=numSame+1
				end if
			end if
		end do
		if(numSame==numSpecies) then
		
			Diff=DiffusivityCompute(defectType, DiffFunc(i,matNum)%functionType, DiffFunc(i,matNum)%numParam,&
				DiffFunc(i,matNum)%parameters, matNum)
				exit
		end if
	end do
	if(i==numFuncDiff(matNum)+1) then

		Diff=0d0
	end if
end if

findDiffusivity=Diff
end function

!*****************************************************************************************
!>double precision function diffusivityCompute - computes diffusivity using a functional form
!! for defects that don't have their diffusivity given by a value in a list.
!!
!! This function has several hard-coded functional forms for diffusivity, including immobile
!! defects, constant functions, and mobile SIA loops. Additional functional forms can be added
!! as needed.
!*****************************************************************************************
double precision function DiffusivityCompute(DefectType, functionType, numParameters, parameters,matNum)
use mod_constants
use DerivedType
implicit none

integer DefectType(numSpecies)
integer functionType, numParameters, matNum,i
double precision parameters(numParameters)
double precision Diff
double precision D0, Em
double precision diffusivityCu

!***************************************************************************************************
!This function computes diffusivity using functional form and parameters given in the input file
!***************************************************************************************************

if(functionType==1) then	!used for immobile defects
	Diff=0d0
else if(functionType==2) then	!used for constant functions
	Diff=parameters(1)
else if(functionType==3) then	!Mobile defect diffusivity
	D0=parameters(1)+parameters(2)/dble(DefectType(3))**(parameters(3))
	Em=parameters(4)+parameters(5)/dble(DefectType(3))**(parameters(6))
	Diff=D0*dexp(-Em/(kboltzmann*temperature))
else if(functionType==4) then	!Cu diffusivity
	!< Dcu(n) = Dcu(1)/n
	Diff=diffusivityCu(matNum)/dble(DefectType(1))
else if(functionType==5) then	!SIA_m (loop, W)
	D0=parameters(1)*dble(DefectType(3))**(-0.5d0)
	Diff=D0*dexp(-parameters(2)/(kboltzmann*temperature))
else
	write(*,*) 'error incorrect diffusivity function chosen'
endif

DiffusivityCompute=Diff

end function

!**********************************************************************************
!This function is used to compute diffusivity of Cu atom
!**********************************************************************************
double precision function diffusivityCu(matNum)
	use DerivedType
	use mod_constants
	implicit none

	integer DefectType(numSpecies)
	integer matNum,i

	outer: do i=1,numSingleDiff(matNum)
		if(DiffSingle(i,matNum)%defectType(1)==1 .AND. DiffSingle(i,matNum)%defectType(2)==0 .AND. &
				DiffSingle(i,matNum)%defectType(3)==0 .AND. DiffSingle(i,matNum)%defectType(4)==0) then

			if(totalDPA > 0d0 .AND. dpaRate > 0d0) then
			!	diffusivityCu=(DiffSingle(i,matNum)%D*dexp(-DiffSingle(i,matNum)%Em/(kboltzmann*temperature)))* &
			!			(Vconcent/initialCeqv)
				diffusivityCu=(DiffSingle(i,matNum)%D*dexp(-DiffSingle(i,matNum)%Em/(kboltzmann*temperature)))* firr
				exit outer
			else
				diffusivityCu=(DiffSingle(i,matNum)%D*dexp(-DiffSingle(i,matNum)%Em/(kboltzmann*temperature)))
				exit outer
			end if

		end if
	end do outer

end function

!***************************************************************************************************
!This function returns the binding energy of defect DefectType() which releases defect product().
!It searches to see if DefectType and product are listed in BindSingle, and if not, in BindFunc.
!If they are in neither, it outputs an error message (cannot dissociate that product from that defect)
!***************************************************************************************************

!*****************************************************************************************
!>double precision function find binding - returns the binding energy of a given defect type
!!
!! This function looks up the binding energy of a given defect type using input data from 
!! material parameter input file. It either looks up binding energy from values in a list or
!! computes binding energy using bindingCompute in the case of a functional form.
!!
!!Input: defect type, product type (what type of defect dissociates from the cluster)
!!Output: binding energy (eV)
!*****************************************************************************************
double precision function findBinding(matNum, DefectType, productType)
use DerivedType
use mod_constants
implicit none

integer DefectType(numSpecies), productType(numSpecies)
integer i, j, numSame, numSameProduct, matNum
double precision Eb
double precision BindingCompute

!Temporary: used as a parameter to vary the binding energy of all defects on GB
double precision, parameter :: Param=0d0

do i=1,numSingleBind(matNum)
	numSame=0
	numSameProduct=0
	do j=1,numSpecies
		if(DefectType(j)==BindSingle(i,matNum)%defectType(j)) then
			numSame=numSame+1
		end if
		if(productType(j)==BindSingle(i,matNum)%product(j)) then
			numSameProduct=numSameProduct+1
		end if
	end do

	if (numSame==numSpecies .AND. numSameProduct==numSpecies) then
		if(matNum==2) then
			
			Eb=BindSingle(i,matNum)%Eb-Param
			exit
		else
			Eb=BindSingle(i,matNum)%Eb
			exit
		end if
	end if
end do

if(i==numSingleBind(matNum)+1) then	!did not find defect in single defect list
	do i=1,numFuncBind(matNum)
		numSame=0
		numSameProduct=0
		do j=1,numSpecies

			if(DefectType(j)==0 .AND. BindFunc(i,matNum)%defectType(j)==0) then
				numSame=numSame+1
			else if(DefectType(j) /= 0 .AND. BindFunc(i,matNum)%defectType(j)==1) then
				if(DefectType(j) >= BindFunc(i,matNum)%min(j) .AND. &
				(DefectType(j) <= BindFunc(i,matNum)%max(j) .OR. BindFunc(i,matNum)%max(j)==-1)) then

					numSame=numSame+1
				end if
			end if

			if(productType(j)==0 .AND. BindFunc(i,matNum)%product(j)==0) then
				numSameProduct=numSameProduct+1
			else if(productType(j) == 1 .AND. BindFunc(i,matNum)%product(j)==1) then	!used to find dissociation binding energy
				numSameProduct=numSameProduct+1
			end if

		end do

		if(numSame==numSpecies .AND. numSameProduct==numSpecies) then
			
			if(matNum==2) then	!Adjust binding energies on GB
			
				Eb=BindingCompute(DefectType, productType, BindFunc(i,matNum)%functionType, BindFunc(i,matNum)%numParam,&
					BindFunc(i,matNum)%parameters)-Param
				exit
				
			else
			
				Eb=BindingCompute(DefectType, productType, BindFunc(i,matNum)%functionType, BindFunc(i,matNum)%numParam,&
					BindFunc(i,matNum)%parameters)
				exit
				
			endif
			
		endif
	end do
	if(i==numFuncBind(matNum)+1) then
		Eb=0d0
	end if
end if

if(Eb < 0d0) then
	Eb=0d0
end if

findBinding=Eb
end function

!*****************************************************************************************
!>double precision function bindingCompute - computes binding energy using a functional form
!! for defects that don't have their binding energy given by a value in a list.
!!
!! This function has several hard-coded functional forms for binding energy, including vacancy
!! clusters, SIA clusters, He/V clusters, and the activation energy for a sessile-glissile
!!SIA loop transformation
!*****************************************************************************************
double precision function BindingCompute(DefectType, product, functionType, numParameters, parameters)
use mod_constants
implicit none

integer DefectType(numSpecies), product(numSpecies)
integer functionType, numParameters, num, CuNum, VNum, SIANum, i
double precision parameters(numParameters)
double precision Eb

!***************************************************************************************************
!This function computes diffusivity using functional form and parameters given in the input file
!***************************************************************************************************

if(functionType==2) then	!used for Cu cluster dislocation

	CuNum=DefectType(1)
	Eb=parameters(1)*kboltzmann-parameters(2)*kboltzmann*tempStore- &
			(36d0*pi)**(1d0/3d0)*atomSize**(2d0/3d0)*parameters(3)*(dble(CuNum)**(2d0/3d0)-dble(CuNum-1)**(2d0/3d0))

else if(functionType==4) then	!V / SIA cluster dislocation
	num=0
	do i=1,numSpecies
		if(DefectType(i) > num) then
			num=DefectType(i)
			exit
		end if
	end do

	Eb=parameters(1)+(parameters(2)-parameters(1))*(dble(num)**(2d0/3d0)-dble(num-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)

else if(functionType==6) then	!nCumV->Cu+(n-1)CumV

	CuNum=DefectType(1)
	VNum=DefectType(2)
	Eb=parameters(1)+parameters(2)*(dble(CuNum)**(0.85d0)-dble(CuNum+1)**(0.85d0))-&
			parameters(3)*(dble(VNum)**(1d0/3d0)-dble(VNum)**(2d0/3d0))

else if(functionType==7) then	!nCumV->V+nCu(m-1)V

	CuNum=DefectType(1)
	VNum=DefectType(2)
	Eb=parameters(1)-parameters(2)*(dble(VNum)**(1d0/3d0)-dble(VNum+1)**(1d0/3d0))+&
			parameters(3)*(dble(VNum)**(2d0/3d0)-dble(VNum+1)**(2d0/3d0))-parameters(4)*dble(CuNum)*&
			(dble(VNum)**(1d0/3d0)-dble(VNum+1)**(1d0/3d0)+dble(VNum)**(2d0/3d0)-dble(VNum+1)**(2d0/3d0))

else
	write(*,*) 'error incorrect Eb function chosen'
endif

BindingCompute=Eb

end function

!*****************************************************************************************
!>integer function findDefectSize - returns the size of a cluster
!!
!!This function will find the effective size of the defect (hard-coded information), used for determining
!!the radius of the defect (for dissociation and clustering reactions).
!!It returns n, the number of lattice spaces taken up by this defect.
!!
!!NOTE: for He_nV_m clusters, this function returns the larger of m or n
!*****************************************************************************************
integer function findDefectSize(defectType)
use mod_constants
implicit none

integer defectType(numSpecies), max, i

!Hard-coded below and may be changed if the rules for defect size change.
max=0
do i=1, numSpecies
	if(defectType(i) > max) then
		max=defectType(i)
	end if
end do

findDefectSize=max
end function
