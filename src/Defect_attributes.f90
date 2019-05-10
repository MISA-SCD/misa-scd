! $Header: /home/CVS//srscd/src/Defect_attributes.f90,v 1.13 2015/10/09 15:36:46 aydunn Exp $

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

double precision function findDiffusivity(matNum, DefectType)
use DerivedType
use mod_srscd_constants
implicit none

integer defectType(numSpecies)
integer i, j, numSame, matNum
double precision Diff
double precision DiffusivityCompute, permanentCv

!Temporary: used as a parameter to vary the diffusivity of all defects on GB
double precision, parameter :: Param=0d0

!***************************************************************************************************
!This function returns the diffusivity of defect DefectType.
!It searches to see if DefectType is listed in DiffSingle, and if not, if it is listed in DiffFunction.
!If it is in neither, it outputs an error message (defect type should not exist)
!***************************************************************************************************

do i=1,numSingleDiff(matNum)
	numSame=0
	do j=1,numSpecies
		if(DefectType(j)==DiffSingle(matNum,i)%defectType(j)) then
			numSame=numSame+1
		endif
	end do
	if (numSame==numSpecies) then
		if(matNum==2) then
			if(DefectType(1)==1 .AND. DefectType(2)==0 .AND. DefectType(3)==0 .AND. DefectType(4)==0 &
					.AND. DPARate > 0) then
				Diff=DiffSingle(matNum,i)%D*dexp(-(DiffSingle(matNum,i)%Em-Param)/(kboltzmann*temperature)) * &
							(permanentCv(matNum) / initialCeqv)
			else
				Diff=DiffSingle(matNum,i)%D*dexp(-(DiffSingle(matNum,i)%Em-Param)/(kboltzmann*temperature))
			end if
			!Diff=DiffSingle(matNum,i)%D*dexp(-(DiffSingle(matNum,i)%Em-Param)/(kboltzmann*temperature))
			exit
		else
			if(DefectType(1)==1 .AND. DefectType(2)==0 .AND. DefectType(3)==0 .AND. DefectType(4)==0 &
					.AND. DPARate > 0) then
				Diff=DiffSingle(matNum,i)%D*dexp(-DiffSingle(matNum,i)%Em/(kboltzmann*temperature)) * &
					(permanentCv(matNum) / initialCeqv)
			else
				Diff=DiffSingle(matNum,i)%D*dexp(-DiffSingle(matNum,i)%Em/(kboltzmann*temperature))
			end if
			!Diff=DiffSingle(matNum,i)%D*dexp(-DiffSingle(matNum,i)%Em/(kboltzmann*temperature))
			exit
		
		end if
	end if
end do

if(i==numSingleDiff(matNum)+1) then	!did not find defect in single defect list
	do i=1,numFuncDiff(matNum)
		numSame=0
		do j=1,numSpecies
			if(DefectType(j)==0 .AND. DiffFunc(matNum,i)%defectType(j)==0) then
				numSame=numSame+1
			else if(DefectType(j) .NE. 0 .AND. DiffFunc(matNum,i)%defectType(j)==1) then
				if(DefectType(j) .GE. DiffFunc(matNum,i)%min(j)) then
				if(DefectType(j) .LE. DiffFunc(matNum,i)%max(j) .OR. DiffFunc(matNum,i)%max(j)==-1) then
					numSame=numSame+1
				endif
				endif
			endif
		end do
		if(numSame==numSpecies) then
		
			Diff=DiffusivityCompute(DefectType, DiffFunc(matNum,i)%functionType, DiffFunc(matNum,i)%numParam,&
				DiffFunc(matNum,i)%parameters, matNum)
				exit
		endif
	end do
	if(i==numFuncDiff(matNum)+1) then
		write(*,*) 'error defect diffusion not allowed'
		write(*,*) DefectType
		Diff=0d0
	endif
endif

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
use mod_srscd_constants
use DerivedType
implicit none

integer DefectType(numSpecies)
integer functionType, numParameters, matNum
double precision parameters(numParameters)
double precision Diff
double precision D0, Em
double precision permanentCv

!***************************************************************************************************
!This function computes diffusivity using functional form and parameters given in the input file
!***************************************************************************************************

if(functionType==1) then
	!used for immobile defects
	Diff=0d0
else if(functionType==2) then
	!used for constant functions
	Diff=parameters(1)
else if(functionType==3) then
	!Mobile defect diffusivity
	D0=parameters(1)+parameters(2)/dble(DefectType(3))**(parameters(3))
	Em=parameters(4)+parameters(5)/dble(DefectType(3))**(parameters(6))
	
	Diff=D0*dexp(-Em/(kboltzmann*temperature))
else if(functionType==5) then
	!< Dcu(n) = Dcu(1)/n
	if(DPARate > 0d0) then
		Diff=(DiffSingle(matNum,1)%D*dexp(-DiffSingle(1,1)%Em/(kboltzmann*temperature)))/dble(DefectType(1)) * &
				(permanentCv(matNum) / initialCeqv)
	else
		Diff=(DiffSingle(matNum,1)%D*dexp(-DiffSingle(1,1)%Em/(kboltzmann*temperature)))/dble(DefectType(1))
	end if
else
	write(*,*) 'error incorrect diffusivity function chosen'
endif

DiffusivityCompute=Diff

end function

!**********************************************************************************
!This function is used to compute the vacancy concentration in the permanent regime
!**********************************************************************************
double precision function permanentCv(matNum)
	use DerivedType
	use mod_srscd_constants
	implicit none

	double precision Kiv, diffV, diffI
	integer matNum

	diffV = DiffSingle(matNum,2)%D*dexp(-DiffSingle(matNum,2)%Em/(kboltzmann*temperature))
	diffI = DiffSingle(matNum,6)%D*dexp(-DiffSingle(matNum,6)%Em/(kboltzmann*temperature))

	Kiv = 4*pi/atomsize*reactionRadius*(diffV + diffI)

	permanentCv = -dislocationDensity*Zint*diffI/(2*Kiv)+&
			((dislocationDensity*Zint*diffI/(2*Kiv))**(2d0)+DPARate*Zint*diffI/(Kiv*diffV))**(1d0/2d0)

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
use mod_srscd_constants
implicit none

integer DefectType(numSpecies), productType(numSpecies)
integer i, j, numSame, numSameProduct, matNum
double precision Eb
double precision BindingCompute

!Temporary: used as a parameter to vary the binding energy of all defects on GB
double precision, parameter :: Param=0d0

do 10 i=1,numSingleBind(matNum)
	numSame=0
	numSameProduct=0
	do 11 j=1,numSpecies
		if(DefectType(j)==BindSingle(matNum,i)%defectType(j)) then
			numSame=numSame+1
		endif
		if(productType(j)==BindSingle(matNum,i)%product(j)) then
			numSameProduct=numSameProduct+1
		endif
	11 continue
	if (numSame==numSpecies .AND. numSameProduct==numSpecies) then
		if(matNum==2) then
			
			Eb=BindSingle(matNum,i)%Eb-Param
			exit
			
		else
		
			Eb=BindSingle(matNum,i)%Eb
			exit
			
		endif
	endif
10 continue

if(i==numSingleBind(matNum)+1) then	!did not find defect in single defect list
	do 12 i=1,numFuncBind(matNum)
		numSame=0
		numSameProduct=0
		do 13 j=1,numSpecies
			if(DefectType(j)==0 .AND. BindFunc(matNum,i)%defectType(j)==0) then
				numSame=numSame+1
			else if(DefectType(j) .NE. 0 .AND. BindFunc(matNum,i)%defectType(j)==1) then
				if(DefectType(j) .GE. BindFunc(matNum,i)%min(j)) then
				if(DefectType(j) .LE. BindFunc(matNum,i)%max(j) .OR. BindFunc(matNum,i)%max(j)==-1) then
					numSame=numSame+1
				endif
				endif
			endif
			if(productType(j)==0 .AND. BindFunc(matNum,i)%product(j)==0) then
				numSameProduct=numSameProduct+1
			else if(productType(j) == 1 .AND. BindFunc(matNum,i)%product(j)==1) then	!used to find dissociation binding energy
				numSameProduct=numSameProduct+1
			else if(productType(j) .NE. 0 .AND. BindFunc(matNum,i)%product(j)==-1) then	!used to identify de-pinning binding energy
				numSameProduct=numSameProduct+1
			endif
		13 continue
		if(numSame==numSpecies .AND. numSameProduct==numSpecies) then
			
			if(matNum==2) then	!Adjust binding energies on GB
			
				Eb=BindingCompute(DefectType, productType, BindFunc(matNum,i)%functionType, BindFunc(matNum,i)%numParam,&
					BindFunc(matNum,i)%parameters)-Param
				exit
				
			else
			
				Eb=BindingCompute(DefectType, productType, BindFunc(matNum,i)%functionType, BindFunc(matNum,i)%numParam,&
					BindFunc(matNum,i)%parameters)
				exit
				
			endif
			
		endif
	12 continue
	if(i==numFuncBind(matNum)+1) then
		write(*,*) 'error dissociation reaction not allowed'
		write(*,*) DefectType
		write(*,*) ProductType
		Eb=0d0
	endif
endif

if(Eb .LT. 0d0) then
	Eb=0d0
endif

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
use mod_srscd_constants
implicit none

integer DefectType(numSpecies), product(numSpecies)
integer functionType, numParameters, num, HeNum, VNum, SIANum, i
double precision parameters(numParameters)
double precision Eb, Eb_VOnly, Eb_HeV

!***************************************************************************************************
!This function computes diffusivity using functional form and parameters given in the input file
!***************************************************************************************************

if(functionType==1) then
	!used for zero functions
	Eb=0d0
else if(functionType==2) then
	!used for Cu cluster dislocation
	Eb=parameters(1)+parameters(2)*(dble(HeNum)**(0.85d0)-dble(HeNum+1)**(0.85d0))
else if(functionType==3) then
	!Mobile SIA loop diffusivity
	write(*,*) 'error no functionType 3 in BindingCompute'
else if(functionType==4) then
	num=0
	do 10 i=1,numSpecies
		if(DefectType(i) .GT. num) then
			num=DefectType(i)
		endif
	10 continue
!
! Reference:
!
	Eb=parameters(1)+(parameters(2)-parameters(1))*(dble(num)**(2d0/3d0)-dble(num-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)
else if(functionType==5) then
	write(*,*) 'error no functionType 5 in BindingCompute'
else if(functionType==6) then
	HeNum=DefectType(1)
	VNum=DefectType(2)
	Eb=parameters(1)+parameters(2)*(dble(HeNum)**(0.85d0)-dble(HeNum+1)**(0.85d0))-&
			parameters(3)*(dble(VNum)**(1d0/3d0)-dble(VNum)**(2d0/3d0))

	!<if(dble(HeNum)/dble(VNum) .LE. .5d0) then
		!helium cannot dissociate from HeV clusters with mostly V
	!<	Eb=10d0
	!<else
		!Use He_mV_n binding energy (does not apply for low m/n ratios)
		!Eb=parameters(1)-parameters(2)*dlog(dble(HeNum)/dble(VNum))/dlog(10d0)-&
		!	parameters(3)*dlog(dble(HeNum)/dble(VNum))**2d0/(dlog(10d0)**2d0) !Marion and Bulatov, from Terentyev
		
		!EDIT: 2014.10.08: Helium cannot dissociate from HeV clusters at all.
	!<	Eb=10d0
		
		
	!<	if(Eb .LT. 0d0) then
	!<		Eb=0d0
	!<	endif
	
	!<endif
else if(functionType==7) then
	HeNum=DefectType(1)
	VNum=DefectType(2)
	Eb=parameters(1)-parameters(2)*(dble(VNum)**(1d0/3d0)-dble(VNum+1)**(1d0/3d0))+&
			parameters(3)*(dble(VNum)**(2d0/3d0)-dble(VNum+1)**(2d0/3d0))-parameters(3)*dble(HeNum)*&
			(dble(VNum)**(1d0/3d0)-dble(VNum+1)**(1d0/3d0)+dble(VNum)**(2d0/3d0)-dble(VNum+1)**(2d0/3d0))

!	if(dble(HeNum)/dble(VNum) .LE. .5d0) then

!		!use vacancy cluster binding energy
!		Eb=parameters(1)+(parameters(2)-parameters(1))*(dble(VNum)**(2d0/3d0)-dble(VNum-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)
		
!		if(Eb .LT. 0d0) then
!			Eb=0d0
!		endif

!	else
	
		!vacancy cluster binding energy
!
! Reference:
!
	!<	Eb_VOnly=parameters(1)+(parameters(2)-parameters(1))*(dble(VNum)**(2d0/3d0)-dble(VNum-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)
		
	!<	if(Eb_VOnly .LT. 0d0) then
	!<		Eb_VOnly=0d0
	!<	endif

		!He_mV_n binding energy (does not apply for low m/n ratios)
		!Eb_HeV=parameters(3)+parameters(4)*dlog(dble(HeNum)/dble(VNum))/dlog(10d0)+&
		!	parameters(5)*dlog(dble(HeNum)/dble(VNum))**2d0/(dlog(10d0)**2d0) !Marion and Bulatov, from Terentyev
	!<	Eb_HeV=parameters(3)*(dble(HeNum)/dble(VNum))**parameters(4)
		
	!<	if(Eb_HeV .LT. 0d0) then
	!<		Eb_HeV=0d0
	!<	endif
		
		!Modification 2015.03.19: choose the larger binding energy between the V_only functional form and the He_V functional form
	
!		if(Eb_VOnly .GT. Eb_HeV) then
!			Eb=Eb_VOnly+1d0
!		else
!			Eb=Eb_HeV+1d0
!		endif

	!<	Eb=Eb_VOnly+Eb_HeV
		
!	endif
	
else if(functionType==8) then
	
	!SIA sessile - glissile binding energy
	!Using (made-up) empirical functional form
	SIANum=DefectType(4)
	
	!2/3 power law
	!Eb=parameters(1)-parameters(2)*(dble(SIANum)**(2d0/3d0)-dble(SIANum-1)**(2d0/3d0))
	
	!linear binding energy dependence
!
! Reference:
!
	Eb=parameters(1)*SIANum+parameters(2)

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
use mod_srscd_constants
implicit none

integer defectType(numSpecies), max, i

!Hard-coded below and may be changed if the rules for defect size change.
max=0
do 10 i=1, numSpecies
	if(defectType(i) .GT. max) then
		max=defectType(i)
	endif
10 continue

findDefectSize=max
end function

!*****************************************************************************************
!>double precision function findStrainEnergy - Returns the interaction energy of the defect with the local strain field
!!
!!This function takes the double dot product of the defect's dipole tensor with the local
!!strain field and returns that amount (energy, in eV). 
!!
!!NOTE: if the dipole tensor is asymmetric, an averaged strain energy is taken which
!!accounts for all possible orientations of the defect. It is not known if this is the 
!!correct averaging procedure that should be used here.
!*****************************************************************************************

double precision function findStrainEnergy(defectType, cell)
use DerivedType
use mod_srscd_constants
implicit none

integer i, j, same
double precision strainE
integer defectType(numSpecies)
integer cell

strainE=0d0

do 10 i=1,numDipole
	
	!search for defect type in dipole tensor
	same=0
	do 11 j=1,numSpecies
		if(defectType(j) .GE. dipoleStore(i)%min(j) .AND. defectType(j) .LE. dipoleStore(i)%max(j)) then
			same=same+1
		endif
	11 continue
	
	if(same==numSpecies) then
		exit
	endif
	
10 continue

if(i .LE. numDipole) then	!we have identified a dipole tensor
	write(*,*) 'finding strain energy'
	write(*,*) defectType
	
	do 12 j=1,6
		if(j .LE. 3) then
			strainE=strainE+myMesh(cell)%strain(j)*dipoleStore(i)%equilib(j)
		else
			strainE=strainE+2d0*myMesh(cell)%strain(j)*dipoleStore(i)%equilib(j)
		endif
		write(*,*) 'j', j, 'strain', myMesh(cell)%strain(j), 'dipole', dipoleStore(i)%equilib(j)
	12 continue

	read(*,*)
endif

findStrainEnergy=strainE

end function


!*****************************************************************************************
!>double precision function findStrainEnergyBoundary - Returns the interaction energy of the defect with the local strain field in a boundary element
!!
!!This function takes the double dot product of the defect's dipole tensor with the local
!!strain field and returns that amount (energy, in eV). 
!!
!!NOTE: if the dipole tensor is asymmetric, an averaged strain energy is taken which
!!accounts for all possible orientations of the defect. It is not known if this is the 
!!correct averaging procedure that should be used here.
!*****************************************************************************************

double precision function findStrainEnergyBoundary(defectType, dir, cell)
use DerivedType
use mod_srscd_constants
implicit none

integer i, j, same
double precision strainE
integer defectType(numSpecies)
integer cell, dir

strainE=0d0

do 10 i=1,numDipole
	
	!search for defect type in dipole tensor
	same=0
	do 11 j=1,numSpecies
		if(defectType(j) .GE. dipoleStore(i)%min(j) .AND. defectType(j) .LE. dipoleStore(i)%max(j)) then
			same=same+1
		endif
	11 continue
	
	if(same==numSpecies) then
		exit
	endif
	
10 continue

if(i .LE. numDipole) then	!we have identified a dipole tensor

	do 12 j=1,6
		if(j .LE. 3) then
			strainE=strainE+myBoundary(dir,cell)%strain(j)*dipoleStore(i)%equilib(j)
		else
			strainE=strainE+2d0*myBoundary(dir,cell)%strain(j)*dipoleStore(i)%equilib(j)
		endif
	12 continue

endif

findStrainEnergyBoundary=strainE

end function
