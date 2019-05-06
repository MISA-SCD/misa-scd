!***************************************************************************************************
!
!> Subroutine readDipoleTensors() - reads defect dipole tensors in from a file and stores them in global variable
!!
!! This subroutine reads in the dipole tensors of various defects from a file. These dipole tensors are
!! then stored in the global variable dipoleTensor for later reference when calculating the interaction
!! energy between a strain field and various defects (vacancies, interstitials, etc.)
!
!***************************************************************************************************

subroutine readDipoleTensors(filename)
use DerivedType
use mod_srscd_constants
implicit none

character*20 char
character*50 filename
integer i,j
logical flag

open(51, file=filename, action='read', status='old')

flag=.FALSE.

!Skip to where data starts
do while (flag .eqv. .FALSE.)
	read(51,*) char
	if(char=='start') then
		flag=.TRUE.
	end if
end do
flag=.FALSE.

do while (flag .eqv. .FALSE.)
	read(51,*) char
	if(char=='numDipoles') then
		read(51,*) numDipole
		flag=.TRUE.
	end if
end do
flag=.FALSE.

allocate(dipoleStore(numDipole))

do i=1,numDipole
	
	allocate(dipoleStore(i)%min(numSpecies))
	allocate(dipoleStore(i)%max(numSpecies))
	
	do while (flag .eqv. .FALSE.)
		read(51,*) char
		if(char=='min') then
			read(51,*) (dipoleStore(i)%min(j),j=1,numSpecies)
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	do while (flag .eqv. .FALSE.)
		read(51,*) char
		if(char=='max') then
			read(51,*) (dipoleStore(i)%max(j),j=1,numSpecies)
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	do while (flag .eqv. .FALSE.)
		read(51,*) char
		if(char=='equilibrium') then
			read(51,*) (dipoleStore(i)%equilib(j), j=1,6)
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.
	
	do while (flag .eqv. .FALSE.)
		read(51,*) char
		if(char=='saddlePoint') then
			read(51,*) (dipoleStore(i)%saddle(j),j=1,6)
			flag=.TRUE.
		end if
	end do
	flag=.FALSE.

end do

close(51)

end subroutine

!***************************************************************************************************
!
!> Subroutine calculateDeltaEm() - uses dipole tensors at equilibrium position and saddle point to calculate
!! delta Em caused by strain field
!!
!! This subroutine calculates the double dot product of the strain field in a given volume element
!! with the dipole tensor of a given defect type at the equilibrium position and the saddle point
!! position to calculate the change in migration energy of that defect.
!!
!! Initially, this subroutine will be set up without averaging the results over several orientations
!! of the dipole tensor. Once we decide how to best average these results, this will be added to the subroutine. 
!
!***************************************************************************************************

double precision function calculateDeltaEm(cellNumber, defectType)
use DerivedType
use mod_srscd_constants
implicit none

double precision deltaE_equilib, deltaE_saddle
integer i,j
integer cellNumber, defectType(numSpecies), same
integer dipoleID

!Find the defect type in the dipole tensor list (if it is there)
dipoleID=0
do 1 i=1,numDipole
	 
	same=0
	do 2 j=1,numSpecies
	
		if(dipoleStore(i)%min(j) .LE. defectType(j)) then
			if(dipoleStore(i)%max(j) .GE. defectType(j)) then
				same=same+1
			endif
		endif
	
	2 continue
	if(same==numSpecies) then
		dipoleID=i
	endif

1 continue

if(dipoleID == 0) then
	!No dipole tensor for defects of this type
	calculateDeltaEm=0d0
else
	!Here we are NOT taking an average over orientations of the dipole tensor (yet)
	!Delta E = -(strain double dot dipole tensor (in eV))
	deltaE_equilib=-(myMesh(cellNumber)%strain(1)*dipoleStore(dipoleID)%equilib(1)+&
					myMesh(cellNumber)%strain(2)*dipoleStore(dipoleID)%equilib(2)+&
					myMesh(cellNumber)%strain(3)*dipoleStore(dipoleID)%equilib(3)+&
					2d0*myMesh(cellNumber)%strain(4)*dipoleStore(dipoleID)%equilib(4)+&
					2d0*myMesh(cellNumber)%strain(5)*dipoleStore(dipoleID)%equilib(5)+&
					2d0*myMesh(cellNumber)%strain(6)*dipoleStore(dipoleID)%equilib(6))
					
	deltaE_saddle =-(myMesh(cellNumber)%strain(1)*dipoleStore(dipoleID)%saddle(1)+&
					myMesh(cellNumber)%strain(2)*dipoleStore(dipoleID)%saddle(2)+&
					myMesh(cellNumber)%strain(3)*dipoleStore(dipoleID)%saddle(3)+&
					2d0*myMesh(cellNumber)%strain(4)*dipoleStore(dipoleID)%saddle(4)+&
					2d0*myMesh(cellNumber)%strain(5)*dipoleStore(dipoleID)%saddle(5)+&
					2d0*myMesh(cellNumber)%strain(6)*dipoleStore(dipoleID)%saddle(6))
					
	calculateDeltaEm=deltaE_saddle-deltaE_equilib
endif

end function


!*****************************************************************************************
!>double precision function find diffusivity strain - returns the diffusivity of a given defect type
!! inside a given volume element including the impact of the strain field
!!
!! This function looks up the diffusivity of a given defect type using input date from 
!! material parameter input file. It either looks up diffusivity from values in a list or
!! computes diffusivity using diffusivityCompute in the case of a functional form.
!!
!! It then computes the interaction of the strain field inside a given volume element with
!! the defect dipole and modifies the migration energy accordingly.
!!
!! If the resulting migration energy is less than zero, it returns a migration energy of zero.
!!
!!Input: defect type, cell ID number
!!Output: diffusivity (nm^2/s)
!*****************************************************************************************

double precision function findDiffusivityStrain(matNum, DefectType, cellNumber)
use DerivedType
use mod_srscd_constants
implicit none

integer defectType(numSpecies)
integer i, j, numSame, matNum, cellNumber
double precision Diff
double precision DiffusivityComputeStrain
double precision calculateDeltaEm
double precision DeltaEm

!Temporary: used as a parameter to vary the diffusivity of all defects on GB
double precision, parameter :: Param=0.0

!***************************************************************************************************
!This function returns the diffusivity of defect DefectType.
!It searches to see if DefectType is listed in DiffSingle, and if not, if it is listed in DiffFunction.
!If it is in neither, it outputs an error message (defect type should not exist)
!***************************************************************************************************

do 10 i=1,numSingleDiff(matNum)
	numSame=0
	do 11 j=1,numSpecies
		if(DefectType(j)==DiffSingle(matNum,i)%defectType(j)) then
			numSame=numSame+1
		endif
	11 continue
	if (numSame==numSpecies) then
		if(matNum==2) then
		
			Diff=DiffSingle(matNum,i)%D*dexp(-(DiffSingle(matNum,i)%Em-Param)/(kboltzmann*temperature))
			exit
		
		else
			
			DeltaEm=calculateDeltaEm(cellNumber, defectType)
			Diff=DiffSingle(matNum,i)%D*dexp(-(DiffSingle(matNum,i)%Em+DeltaEm)/(kboltzmann*temperature))
			exit
		
		endif
	endif
10 continue

if(i==numSingleDiff(matNum)+1) then	!did not find defect in single defect list
	do 12 i=1,numFuncDiff(matNum)
		numSame=0
		do 13 j=1,numSpecies
			if(DefectType(j)==0 .AND. DiffFunc(matNum,i)%defectType(j)==0) then
				numSame=numSame+1
			else if(DefectType(j) .NE. 0 .AND. DiffFunc(matNum,i)%defectType(j)==1) then
				if(DefectType(j) .GE. DiffFunc(matNum,i)%min(j)) then
				if(DefectType(j) .LE. DiffFunc(matNum,i)%max(j) .OR. DiffFunc(matNum,i)%max(j)==-1) then
					numSame=numSame+1
				endif
				endif
			endif
		13 continue
		if(numSame==numSpecies) then
		
			Diff=DiffusivityComputeStrain(DefectType, DiffFunc(matNum,i)%functionType, DiffFunc(matNum,i)%numParam,&
				DiffFunc(matNum,i)%parameters, cellNumber)
				exit
		endif
	12 continue
	if(i==numFuncDiff(matNum)+1) then
		write(*,*) 'error defect diffusion not allowed'
		write(*,*) DefectType
		Diff=0d0
	endif
endif

findDiffusivityStrain=Diff
end function

!*****************************************************************************************
!>double precision function diffusivityComputeStrain - computes diffusivity using a functional form
!! for defects that don't have their diffusivity given by a value in a list. Includes the effect
!! of strain interacting with defect dipole tensor
!!
!! This function has several hard-coded functional forms for diffusivity, including immobile
!! defects, constant functions, and mobile SIA loops. Additional functional forms can be added
!! as needed.
!!
!! Strain is included by modifying the migration energy Em by an amount equal to the defect dipole
!! tensor double dot product with the strain field at the equilibrium and saddle points
!*****************************************************************************************

double precision function DiffusivityComputeStrain(DefectType, functionType, numParameters, parameters, cellNumber)
use mod_srscd_constants
implicit none

integer DefectType(numSpecies)
integer functionType, numParameters
double precision parameters(numParameters)
double precision Diff
double precision D0, Em
double precision deltaEm
double precision calculateDeltaEm
integer cellNumber

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
	deltaEm=calculateDeltaEm(cellNumber,defectType)
	
	Diff=D0*dexp(-(Em+deltaEm)/(kboltzmann*temperature))
else
	write(*,*) 'error incorrect diffusivity function chosen'
endif

DiffusivityComputeStrain=Diff

end function
