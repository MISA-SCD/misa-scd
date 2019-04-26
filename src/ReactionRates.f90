! $Header: /home/CVS//srscd/src/ReactionRates.f90,v 1.19 2015/10/19 19:43:22 aydunn Exp $
!***************************************************************************************************
!>Module: ReactionRates - adds reactions and calculates rates for the various reaction types in the system.
!!
!!This module contains hard-coded information. It contains the allowed reactions
!!(as read in from the input file) and thier rates.
!!
!!This module also contains addSingleDefectReactions, addMultiDefectReactions, and addDiffusionReactions.
!!These subroutines will, given (a) defect type(s) and a cell number, calculate the reaction rate of
!!all single/multi-defect/diffusion reactions possible in that element and add those rates to the 
!!reaction list in that cell. They shold also delete the reactions that have new reaction rates of 0 in
!!the cell (reactions that are no longer possible after a reaction has occurred).
!
!2014/05/23
!***************************************************************************************************

module ReactionRates
implicit none
contains

!>Subroutine add single defect reactions - adds reactions to a reaction list that require only
!! a single reactant to be carried out. Examples: dissociation, trapping, sinks. Diffusion
!! reactions are not included in this subroutine.
!!
!! Inputs: cell number, defect type (we are adding all single defect reactions associated with 
!! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect type corresponds to an allowed reaction using reaction lists
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect type and number of defects
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.

subroutine addSingleDefectReactions(cell, defectType)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), matNum
integer i, j, k, count, numReactants, numProducts, storeTemp
type(reaction), pointer :: reactionCurrent, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate, totalRateCheck
logical isLegal

!Dissociation reactions. NOTE: the number of reactants and number of products, as well as the form 
!of the products given the reactants, is hard-coded into this section. This type of hard-coding is
!carried out in each section in this module, and is seen here as unavoidable.

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
endif

numReactants=1
numProducts=2
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))


do 10 i=1, numDissocReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do 11 j=1,numSpecies
		if(defectType(j) == 0 .AND. DissocReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. DissocReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. DissocReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. DissocReactions(matNum,i)%max(j)) .OR. DissocReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	11 continue
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do 12 j=1,numSpecies
			
			reactants(1,j)=defectType(j)
			
			products(2,j)=DissocReactions(matNum,i)%products(1,j)
			
			products(1,j)=reactants(1,j)-products(2,j)

		12 continue
		
		!************************************************************************************
		!Hard coded: Frenkel pair creation from HeV cluster with He/V .GT. 5
		!HARD CODED: what to do with dissociation of mobile defects from sessile SIA clusters
		!************************************************************************************
		
		if(products(1,3) .LT. 0) then
			if(products(1,1) .NE. 0 .AND. products(1,2) .NE. 0) then	!Kick-out mechanism
				products(1,3)=0
				products(1,2)=products(1,2)+products(2,3)
			else														!dissociation of mobile defects from sessile SIA clusters
				products(1,3)=0
				products(1,4)=products(1,4)-products(2,3)
			endif
		endif
		
		!***********************************
		!Hard coded: mobilize SIA clusters that decrease in size below max3Dint
		!***********************************
		
		!sessile cluster becomes mobile when it shrinks below max3DInt
		if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt) then
			products(1,3)=products(1,4)
			products(1,4)=0
		endif
		
		reactionCurrent=>reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!*************************************************************************************
		!HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		!HARD CODED: if He/V ratio .GE. 5 then SIA kick-out activated
		!HARD CODED: disallow He1V1 dissociation where the product is V1, to avoid double-counting
		!*************************************************************************************
		
		call checkReactionLegality(numProducts, products, isLegal)
		
		!find the reaction rate. 
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateDissoc(defectType, products, cell, DissocReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!Point reactionSort and reactionSortPrev at the two reactions that reactionCurrent should be
			!inserted between
			
			!Step 1: skip the first 1 or 2 reactions (don't want to change order of implantation reactions)
!			if(HeDPARatio .GT. 0d0) then
!				reactionSort=>reactionList(cell)%next%next
!				reactoinSortPrev=>reactionList(cell)%next
!			else
!				reactionSort=>reactionList(cell)%next
!				reactionSortPrev=>reactionList(cell)
!			endif
			
!			!Step 2: point reactionSort and reactionSortPrev where reactionCurrent should be inserted
!			call sortReactionList(reactionSort, reactionSortPrev, reactionRate)
			
!			!Step 3: create new reaction in between reactionSort and reactionSortPrev
!			nullify(reactionSortPrev%next)
!			allocate(reactionSortPrev%next)
!			reactionCurrent=>reactionSortPrev%next
!			if(associated(reactionSort)) then
!				reactionCurrent%next=>reactionSort
!			else
!				nullify(reactionCurrent%next)
!			endif
			
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=2
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%products(1,j)=products(1,j)
				reactionCurrent%products(2,j)=products(2,j)
			13 continue
			do 14 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			14 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
		
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - dissociation'
		endif

	endif
10 continue

!Sink reactions

numReactants=1
numProducts=0
deallocate(reactants)
deallocate(products)
allocate(reactants(numReactants,numSpecies))

do 20 i=1, numSinkReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this sink reaction
	do 21 j=1,numSpecies
		if(defectType(j) == 0 .AND. SinkReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. SinkReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. SinkReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. SinkReactions(matNum,i)%max(j)) .OR. SinkReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	21 continue
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (sinks)
		do 22 j=1,numSpecies
			reactants(1,j)=defectType(j)
		22 continue
		
		reactionCurrent=>reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateSink(defectType, cell, SinkReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			if(allocated(reactionCurrent%products)) then
				deallocate(reactionCurrent%products)
			endif
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=0
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 23 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
			23 continue
			do 24 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			24 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
		
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		endif

	endif
	
20 continue

!Impurity reactions

numReactants=1
numProducts=1
deallocate(reactants)
if(allocated(products)) then
	deallocate(products)
endif
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 30 i=1, numImpurityReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this impurity reaction
	do 31 j=1,numSpecies
		if(defectType(j) == 0 .AND. ImpurityReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. ImpurityReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. ImpurityReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. ImpurityReactions(matNum,i)%max(j)) .OR. ImpurityReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	31 continue
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (impurities)
		!Impurities change defect types from glissile SIA loops to sesile SIA loops. Therefore,
		!we must change the defectType from 0 0 n 0 to 0 0 0 n. This is hard-coded in here.
		do 32 j=1,numSpecies
			reactants(1,j)=defectType(j)
			if(reactants(1,j) .NE. 0d0) then
				storeTemp=reactants(1,j)
			endif
		32 continue
		
		do 35 j=1,numSpecies
			if(ImpurityReactions(matNum,i)%products(1,j)==1) then
				products(1,j)=storeTemp
			else
				products(1,j)=0
			endif
		35 continue
		
		reactionCurrent=>reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateImpurity(defectType, cell, ImpurityReactions(matNum,i))
		else
			reactionRate=0d0
		endif
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 33 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%products(1,j)=products(1,j)
			33 continue
			do 34 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			34 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		endif

	endif
30 continue

end subroutine

!***************************************************************************************************
!
! Subroutine addSingleDefectReactionsFine
!
! Updates reaction rates associated with defectType in cell in fine mesh of cascade. Single defect
! reactions only (ex. dissociation, trapping at sinks, etc.)
!
! Inputs: cascade (cascade ID number), cell (volume element ID number), defectType
! Outputs: updates single-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!
!***************************************************************************************************

!>Subroutine add single defect reactions fine - adds reactions to a reaction list inside a cascade (fine mesh) that require only
!! a single reactant to be carried out. Examples: dissociation, trapping, sinks. Diffusion
!! reactions are not included in this subroutine.
!!
!! Inputs: cascade ID number, cell number, defect type (we are adding all single defect reactions associated with 
!! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect type corresponds to an allowed reaction using reaction lists
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect type and number of defects (using fine mesh volume element size)
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.

subroutine addSingleDefectReactionsFine(cascadeID, cell, defectType)
use mod_srscd_constants
use DerivedType
implicit none

integer cascadeID, cell, defectType(numSpecies)

type(cascade), pointer :: CascadeCurrent

integer i, j, count, numReactants, numProducts, storeTemp, matNum
type(reaction), pointer :: reactionCurrent, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate, totalRateCheck
logical isLegal

!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this 
!subroutine
CascadeCurrent=>ActiveCascades

do 40 while(associated(CascadeCurrent))

	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	endif
	
	CascadeCurrent=>CascadeCurrent%next
40 continue

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material
endif


!Dissociation reactions. NOTE: the number of reactants and number of products, as well as the form 
!of the products given the reactants, is hard-coded into this section. This type of hard-coding is
!carried out in each section in this module, and is seen here as unavoidable.

numReactants=1
numProducts=2
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 10 i=1, numDissocReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do 11 j=1,numSpecies
		if(defectType(j) == 0 .AND. DissocReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. DissocReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. DissocReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. DissocReactions(matNum,i)%max(j)) .OR. DissocReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	11 continue
		
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)

		do 12 j=1,numSpecies
			
			reactants(1,j)=defectType(j)
			
			products(2,j)=DissocReactions(matNum,i)%products(1,j)
			products(1,j)=reactants(1,j)-products(2,j)
			
		12 continue
		
		!************************************************************************************
		!Hard coded: Frenkel pair creation from HeV cluster with He/V .GT. 5
		!HARD CODED: what to do with dissociation of mobile defects from sessile SIA clusters
		!************************************************************************************
		
		if(products(1,3) .LT. 0) then
			if(products(1,1) .NE. 0 .AND. products(1,2) .NE. 0) then	!Kick-out mechanism
				products(1,3)=0
				products(1,2)=products(1,2)+products(2,3)
			else														!dissociation of mobile defects from sessile SIA clusters
				products(1,3)=0
				products(1,4)=products(1,4)-products(2,3)
			endif
		endif

		!***********************************
		!Hard coded: mobilize SIA clusters that decrease in size below max3Dint
		!***********************************
		
		!sessile cluster becomes mobile when it shrinks below max3DInt
		if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt) then
			products(1,3)=products(1,4)
			products(1,4)=0
		endif
				
		reactionCurrent=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, DissocReactions(matNum,i))
		else
			reactionRate=0d0
		endif
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=2
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%products(1,j)=products(1,j)
				reactionCurrent%products(2,j)=products(2,j)
			13 continue
			do 14 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			14 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - dissociation'
		endif

	endif
10 continue

!Sink reactions

numReactants=1
numProducts=0
deallocate(reactants)
deallocate(products)
allocate(reactants(numReactants,numSpecies))

do 20 i=1, numSinkReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this sink reaction
	do 21 j=1,numSpecies
		if(defectType(j) == 0 .AND. SinkReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. SinkReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. SinkReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. SinkReactions(matNum,i)%max(j)) .OR. SinkReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	21 continue
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (sinks)
		do 22 j=1,numSpecies
			reactants(1,j)=defectType(j)
		22 continue
		
		reactionCurrent=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateSinkFine(CascadeCurrent, defectType, cell, SinkReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			if(allocated(reactionCurrent%products)) then
				deallocate(reactionCurrent%products)
			endif
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=0
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 23 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
			23 continue
			do 24 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			24 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		endif

	endif
	
20 continue

!Impurity reactions

numReactants=1
numProducts=1
deallocate(reactants)
if(allocated(products)) then
	deallocate(products)
endif
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 30 i=1, numImpurityReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this impurity reaction
	do 31 j=1,numSpecies
		if(defectType(j) == 0 .AND. ImpurityReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. ImpurityReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. ImpurityReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. ImpurityReactions(matNum,i)%max(j)) .OR. ImpurityReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	31 continue
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (impurities)
		!Impurities change defect types from glissile SIA loops to sesile SIA loops. Therefore,
		!we must change the defectType from 0 0 n 0 to 0 0 0 n. This is hard-coded in here.
		do 32 j=1,numSpecies
			reactants(1,j)=defectType(j)
			if(reactants(1,j) .NE. 0d0) then
				storeTemp=reactants(1,j)
			endif
		32 continue
		
		do 35 j=1,numSpecies
			if(ImpurityReactions(matNum,i)%products(1,j)==1) then
				products(1,j)=storeTemp
			else
				products(1,j)=0
			endif
		35 continue
		
		reactionCurrent=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateImpurityFine(CascadeCurrent, defectType, cell, ImpurityReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 33 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%products(1,j)=products(1,j)
			33 continue
			do 34 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			34 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		endif
		
	endif
30 continue

end subroutine

!>Subroutine add multi defect reactions - adds reactions to a reaction list that require multiple
!!(only using 2 currently) defects to be carried out. This refers mainly to clustering reactions 
!!or pinning reactions.
!!
!! Inputs: cell number, defect types (we are adding all clustering reactions associated with 
!! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect types both correspond to an allowed clustering reaction reaction using reaction lists
!!
!!1a) Calculate the resulting product (if the reaction is a clustering reaction), using combination 
!!rules. Examples include annihilation of vacancy/interstitial pairs and not allowing He-SIA 
!!clusters to form.
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect types and number of defects 
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.

subroutine addMultiDefectReactions(cell, defectType1, defectType2)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), matNum
type(reaction), pointer :: reactionCurrent, reactionPrev
integer i, j, count, numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate
logical isLegal

integer findNumDefect

!Clustering reactions. NOTE: the number of reactants and number of products, as well as the form 
!of the products given the reactants, is hard-coded into this section. This type of hard-coding is
!carried out in each section in this module, and is seen here as unavoidable.

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
endif


do 10 i=1, numClusterReac(matNum)
	
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	
	do 11 j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
				count=count+1
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
				if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies)) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
						ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
						count=count+1
					endif
				endif
			endif
		else if(defectType1(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
				if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j)) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						count=count+1
					endif
				endif
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
				if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
					ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies) .AND. &
							defectType1(j) .GE. ClusterReactions(matNum,i)%min(j)) then
							count=count+1
						endif
					endif
				endif
			endif
		endif
	11 continue
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction
	
		numReactants=2
		allocate(reactants(numReactants,numSpecies))
		
		!Hard coded: if we are allowing SIA pinning by HeV clusters, then we have 2 products
		if(defectType1(1) .NE. 0 .AND. defectType2(3) .GE. SIAPinMin .AND. SIAPinToggle=='yes') then
			if(defectType1(2)==0 .OR. dble(defectType1(1))/dble(defectType1(2)) .GT. 0.5d0) then
				numProducts=2
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do 121 j=1,numSpecies
					reactants(1,j)=defectType1(j)
					reactants(2,j)=defectType2(j)
					products(1,j)=defectType1(j)
					if(j==3) then
						products(2,j)=0
					else if(j==4) then
						products(2,j)=defectType2(3)
					else
						products(2,j)=defectType2(j)
					endif
				121 continue
			else
				numProducts=1
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do 512 j=1,numSpecies
					reactants(1,j)=defectType1(j)
					reactants(2,j)=defectType2(j)
					products(1,j)=reactants(1,j)+reactants(2,j)
				512 continue
			endif
		else
			numProducts=1
			allocate(products(numProducts,numSpecies))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do 12 j=1,numSpecies
				reactants(1,j)=defectType1(j)
				reactants(2,j)=defectType2(j)
				products(1,j)=reactants(1,j)+reactants(2,j)
			12 continue
		endif
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: He, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of HeV+SIA=>HeSIA (should not be allowed here)
		!*******************************************************************************************
		
		!Vacancy+SIA annihilation - only the larger species remains
		if(products(1,2) .GE. products(1,3)) then
			products(1,2)=products(1,2)-products(1,3)
			products(1,3)=0
		endif
		
		if(products(1,2) .GE. products(1,4)) then
			products(1,2)=products(1,2)-products(1,4)
			products(1,4)=0
		endif
		
		if(products(1,2) .LT. products(1,3)) then
			products(1,3)=products(1,3)-products(1,2)
			products(1,2)=0
		endif
		
		if(products(1,2) .LT. products(1,4)) then
			products(1,4)=products(1,4)-products(1,2)
			products(1,2)=0
		endif
		
		!SIA+SIA clustering
		!two 1D clusters coming together to make a sessile cluster
		if(reactants(1,3) .GT. max3DInt .AND. reactants(2,3) .GT. max3DInt) then
			products(1,4)=products(1,3)
			products(1,3)=0
		endif
		
		!sessile SIA + mobile SIA = sessile SIA
		if(products(1,3) .NE. 0 .AND. products(1,4) .NE. 0) then
			products(1,4)=products(1,3)+products(1,4)
			products(1,3)=0
		endif
		
		if(products(1,1) .NE. 0) then

			!All He-SIA clusters are immobile
			if(products(1,3) .NE. 0) then
				products(1,4)=products(1,3)+products(1,4)
				products(1,3)=0
			endif

		else
		
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt .AND. SIAPinToggle=='no') then
				products(1,3)=products(1,4)
				products(1,4)=0
			endif
		
		endif
		
		!Total Annihilation
		count=0
		do 15 j=1,numSpecies
			if(products(1,j)==0) then
				count=count+1
			endif
		15 continue
		if(count==numSpecies) then
			!we have completely annihilated the defects
			numProducts=0
			deallocate(products)
			allocate(products(numProducts, numSpecies))
		endif
		
		!findReactionInList points reactionCurrent at the reaction if it already exists. If not, reactionCurrent
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		
		reactionCurrent=>reactionList(cell)
		call findReactionInListMultiple(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(matNum,i))
		else
			reactionRate=0d0
!			write(*,*) 'illegal reaction ', 'reactants', reactants, 'products', products
!			if(associated(reactionCurrent)) then
!				write(*,*) 'reactionCurrent associated', reactionCurrent%reactants, 'products', reactionCurrent%products, &
!				'rate', reactionCurrent%reactionRate
!			endif
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate

			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=numProducts
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%reactants(2,j)=reactants(2,j)
				if(numProducts==1) then
					reactionCurrent%products(1,j)=products(1,j)
				else if(numProducts==2) then
					reactionCurrent%products(1,j)=products(1,j)
					reactionCurrent%products(2,j)=products(2,j)
				endif
			13 continue
			do 14 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			14 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
	
	deallocate(reactants)
	deallocate(products)
		
	endif
	
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	
	do 21 j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
				count=count+1
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
				if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j)) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						count=count+1
					endif
				endif
			endif
		else if(defectType1(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
				if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies)) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
						ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
						count=count+1
					endif
				endif
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
				if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
					ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies) .AND. &
							defectType2(j) .GE. ClusterReactions(matNum,i)%min(j)) then
							count=count+1
						endif
					endif
				endif
			endif
		endif
	21 continue
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction
		
		numReactants=2
		allocate(reactants(numReactants,numSpecies))
		
		!Hard coded: if we are allowing SIA pinning by HeV clusters, then we have 2 products
		if(defectType2(1) .NE. 0 .AND. defectType1(3) .GE. SIAPinMin .AND. SIAPinToggle=='yes') then
			if(defectType2(2)==0 .OR. dble(defectType2(1))/dble(defectType2(2)) .GT. 0.5d0) then
				numProducts=2
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (pinning)
				!NOTE: reverse the order of reactants and products so that the reaction is the same
				!(so that findreactioninlistmultiple correctly identifies the reaction)
				do 122 j=1,numSpecies
					reactants(2,j)=defectType1(j)
					reactants(1,j)=defectType2(j)
					products(1,j)=defectType2(j)
					if(j==3) then
						products(2,j)=0
					else if(j==4) then
						products(2,j)=defectType1(3)
					else
						products(2,j)=defectType1(j)
					endif
				122 continue
			else
				numProducts=1
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (clustering)
				!NOTE: reverse the order of reactants and products so that the reaction is the same
				!(so that findreactioninlistmultiple correctly identifies the reaction)
				do 522 j=1,numSpecies
					reactants(2,j)=defectType1(j)
					reactants(1,j)=defectType2(j)
					products(1,j)=reactants(1,j)+reactants(2,j)
				522 continue
			endif
		else
			numProducts=1
			allocate(products(numProducts,numSpecies))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			!NOTE: reverse the order of reactants and products so that the reaction is the same
			!(so that findreactioninlistmultiple correctly identifies the reaction)
			do 22 j=1,numSpecies
				reactants(2,j)=defectType1(j)
				reactants(1,j)=defectType2(j)
				products(1,j)=reactants(1,j)+reactants(2,j)
			22 continue
		endif
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: He, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of HeV+SIA=>HeSIA (should not be allowed here)
		!*******************************************************************************************
		!Vacancy+SIA annihilation - only the larger species remains
		if(products(1,2) .GE. products(1,3)) then
			products(1,2)=products(1,2)-products(1,3)
			products(1,3)=0
		endif
		if(products(1,2) .GE. products(1,4)) then
			products(1,2)=products(1,2)-products(1,4)
			products(1,4)=0
		endif
		if(products(1,2) .LT. products(1,3)) then
			products(1,3)=products(1,3)-products(1,2)
			products(1,2)=0
		endif
		if(products(1,2) .LT. products(1,4)) then
			products(1,4)=products(1,4)-products(1,2)
			products(1,2)=0
		endif
		
		!SIA+SIA clustering
		!two 1D clusters coming together to make a sessile cluster
		if(reactants(1,3) .GT. max3DInt .AND. reactants(2,3) .GT. max3DInt) then
			products(1,4)=products(1,3)
			products(1,3)=0
		endif
		
		!sessile SIA + mobile SIA = sessile SIA
		if(products(1,3) .NE. 0 .AND. products(1,4) .NE. 0) then
			products(1,4)=products(1,3)+products(1,4)
			products(1,3)=0
		endif
		
		if(products(1,1) .NE. 0) then

			!All He-SIA clusters are immobile
			if(products(1,3) .NE. 0) then
				products(1,4)=products(1,3)+products(1,4)
				products(1,3)=0
			endif

		else
		
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt .AND. SIAPinToggle=='no') then
				products(1,3)=products(1,4)
				products(1,4)=0
			endif
		
		endif
		
		!Total Annihilation
		count=0
		do 25 j=1,numSpecies
			if(products(1,j)==0) then
				count=count+1
			endif
		25 continue
		if(count==numSpecies) then
			!we have completely annihilated the defects
			numProducts=0
			deallocate(products)
			allocate(products(numProducts, numSpecies))
		endif
		
		!findReactionInList points reactionCurrent at the reaction if it already exists. If not, reactionCurrent
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		
		reactionCurrent=>reactionList(cell)
		call findReactionInListMultiple(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA left over, need to disallow this reaction
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(matNum,i))
		else
			reactionRate=0d0
!			write(*,*) 'illegal reaction ', 'reactants', reactants, 'products', products
!			if(associated(reactionCurrent)) then
!				write(*,*) 'reactionCurrent associated', reactionCurrent%reactants, 'products', reactionCurrent%products, &
!				'rate', reactionCurrent%reactionRate
!			endif
		endif
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate

			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=numProducts
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 23 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%reactants(2,j)=reactants(2,j)
				if(numProducts==1) then
					reactionCurrent%products(1,j)=products(1,j)
				else if(numProducts==2) then
					reactionCurrent%products(1,j)=products(1,j)
					reactionCurrent%products(2,j)=products(2,j)
				endif
			23 continue
			do 24 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			24 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
	
	deallocate(reactants)
	deallocate(products)
		
	endif
	
10 continue

end subroutine

!***************************************************************************************************
!
! Subroutine addMultiDefectReactionsFine
!
! Updates reaction rates associated with defectType1 and defectType2 in cell in fine mesh of cascade. 
! Multi defect reactions only (ex. clustering)
!
! Inputs: cascade (cascade ID number), cell (volume element ID number), defectType1, defectType2
! Outputs: updates multi-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!
!***************************************************************************************************

!>Subroutine add multi defect reactions (fine mesh) - adds reactions to a reaction list inside a cascade that require multiple
!!(only using 2 currently) defects to be carried out. This refers mainly to clustering reactions 
!!or pinning reactions.
!!
!! Inputs: cascade ID number, cell number, defect types (we are adding all clustering reactions associated with 
!! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect types both correspond to an allowed clustering reaction reaction using reaction lists
!!
!!1a) Calculate the resulting product (if the reaction is a clustering reaction), using combination 
!!rules. Examples include annihilation of vacancy/interstitial pairs and not allowing He-SIA 
!!clusters to form.
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect types and number of defects (and using the cascade volume)
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.

subroutine addMultiDefectReactionsFine(cascadeID, cell, defectType1, defectType2)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), cascadeID, matNum
type(cascade), pointer :: CascadeCurrent
type(reaction), pointer :: reactionCurrent, reactionPrev
integer i, j, count, numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionrate
logical isLegal

!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this 
!subroutine
CascadeCurrent=>ActiveCascades

do 40 while(associated(CascadeCurrent))

	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	endif
	
	CascadeCurrent=>CascadeCurrent%next
40 continue

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material
endif

!Clustering reactions. NOTE: the number of reactants and number of products, as well as the form 
!of the products given the reactants, is hard-coded into this section. This type of hard-coding is
!carried out in each section in this module, and is seen here as unavoidable.

do 10 i=1, numClusterReac(matNum)
	
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	
	do 11 j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
				count=count+1
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
				if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies)) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
						ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
						count=count+1
					endif
				endif
			endif
		else if(defectType1(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
				if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j)) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						count=count+1
					endif
				endif
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
				if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
					ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies) .AND. &
							defectType1(j) .GE. ClusterReactions(matNum,i)%min(j)) then
							count=count+1
						endif
					endif
				endif
			endif
		endif
	11 continue
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction
		
		numReactants=2
		allocate(reactants(numReactants,numSpecies))
		
		!Hard coded: if we are allowing SIA pinning by HeV clusters, then we have 2 products
		if(defectType1(1) .NE. 0 .AND. defectType2(3) .GE. SIAPinMin .AND. SIAPinToggle=='yes') then
			if(defectType1(2)==0 .OR. dble(defectType1(1))/dble(defectType1(2)) .GT. 0.5d0) then
				numProducts=2
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do 121 j=1,numSpecies
					reactants(1,j)=defectType1(j)
					reactants(2,j)=defectType2(j)
					products(1,j)=defectType1(j)
					if(j==3) then
						products(2,j)=0
					else if(j==4) then
						products(2,j)=defectType2(3)
					else
						products(2,j)=defectType2(j)
					endif
				121 continue
			else
				numProducts=1
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do 512 j=1,numSpecies
					reactants(1,j)=defectType1(j)
					reactants(2,j)=defectType2(j)
					products(1,j)=reactants(1,j)+reactants(2,j)
				512 continue
			endif
		else
			numProducts=1
			allocate(products(numProducts,numSpecies))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do 12 j=1,numSpecies
				reactants(1,j)=defectType1(j)
				reactants(2,j)=defectType2(j)
				products(1,j)=reactants(1,j)+reactants(2,j)
			12 continue
		endif
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: He, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of HeV+SIA=>HeSIA (should not be allowed here)
		!*******************************************************************************************
		!Vacancy+SIA annihilation - only the larger species remains
		if(products(1,2) .GE. products(1,3)) then
			products(1,2)=products(1,2)-products(1,3)
			products(1,3)=0
		endif
		if(products(1,2) .GE. products(1,4)) then
			products(1,2)=products(1,2)-products(1,4)
			products(1,4)=0
		endif
		if(products(1,2) .LT. products(1,3)) then
			products(1,3)=products(1,3)-products(1,2)
			products(1,2)=0
		endif
		if(products(1,2) .LT. products(1,4)) then
			products(1,4)=products(1,4)-products(1,2)
			products(1,2)=0
		endif
		
		!SIA+SIA clustering
		!two 1D clusters coming together to make a sessile cluster
		if(reactants(1,3) .GT. max3DInt .AND. reactants(2,3) .GT. max3DInt) then
			products(1,4)=products(1,3)
			products(1,3)=0
		endif
		
		!sessile SIA + mobile SIA = sessile SIA
		if(products(1,3) .NE. 0 .AND. products(1,4) .NE. 0) then
			products(1,4)=products(1,3)+products(1,4)
			products(1,3)=0
		endif
		
		if(products(1,1) .NE. 0) then

			!All He-SIA clusters are immobile
			if(products(1,3) .NE. 0) then
				products(1,4)=products(1,3)
				products(1,3)=0
			endif

		else
		
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt .AND. SIAPinToggle=='no') then
				products(1,3)=products(1,4)
				products(1,4)=0
			endif
		
		endif
		
		!Total Annihilation
		count=0
		do 15 j=1,numSpecies
			if(products(1,j)==0) then
				count=count+1
			endif
		15 continue
		if(count==numSpecies) then
			!we have completely annihilated the defects
			numProducts=0
			deallocate(products)
			allocate(products(numProducts, numSpecies))
		endif
		
		!findReactionInList points reactionCurrent at the reaction if it already exists. If not, reactionCurrent
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.

		reactionCurrent=>CascadeCurrent%reactionList(cell)
		call findReactionInListMultiple(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, ClusterReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=numProducts
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%reactants(2,j)=reactants(2,j)
				if(numProducts==1) then
					reactionCurrent%products(1,j)=products(1,j)
				else if(numProducts==2) then
					reactionCurrent%products(1,j)=products(1,j)
					reactionCurrent%products(2,j)=products(2,j)
				endif
			13 continue
			do 14 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			14 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			!update totalRate
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
		
	deallocate(reactants)
	deallocate(products)

	endif
	
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	
	do 21 j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(matNum,i)%reactants(2,j)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
				count=count+1
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
				if(defectType2(j) .GE. ClusterReactions(matNum,i)%min(j)) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. &
						ClusterReactions(matNum,i)%max(j)==-1) then
						count=count+1
					endif
				endif
			endif
		else if(defectType1(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(2,j) .NE. 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(matNum,i)%reactants(1,j)==0) then
				if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies)) then
					if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
						ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
						count=count+1
					endif
				endif
			else if(defectType2(j) .NE. 0 .AND. ClusterReactions(matNum,i)%reactants(1,j) .NE. 0) then
				if((defectType1(j) .LE. ClusterReactions(matNum,i)%max(j+numSpecies)) .OR. &
					ClusterReactions(matNum,i)%max(j+numSpecies)==-1) then
					if((defectType2(j) .LE. ClusterReactions(matNum,i)%max(j)) .OR. ClusterReactions(matNum,i)%max(j)==-1) then
						if(defectType1(j) .GE. ClusterReactions(matNum,i)%min(j+numSpecies) .AND. &
							defectType2(j) .GE. ClusterReactions(matNum,i)%min(j)) then
							count=count+1
						endif
					endif
				endif
			endif
		endif
	21 continue
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction
		
		numReactants=2
		allocate(reactants(numReactants,numSpecies))
		
		!Hard coded: if we are allowing SIA pinning by HeV clusters, then we have 2 products
		if(defectType2(1) .NE. 0 .AND. defectType1(3) .GE. SIAPinMin .AND. SIAPinToggle=='yes') then
			if(defectType2(2)==0 .OR. dble(defectType2(1))/dble(defectType2(2)) .GT. 0.5d0) then
				numProducts=2
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do 122 j=1,numSpecies
					reactants(2,j)=defectType1(j)
					reactants(1,j)=defectType2(j)
					products(1,j)=defectType2(j)
					if(j==3) then
						products(2,j)=0
					else if(j==4) then
						products(2,j)=defectType1(3)
					else
						products(2,j)=defectType1(j)
					endif
				122 continue
			else
				numProducts=1
				allocate(products(numProducts,numSpecies))
				
				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do 522 j=1,numSpecies
					reactants(2,j)=defectType1(j)
					reactants(1,j)=defectType2(j)
					products(1,j)=reactants(1,j)+reactants(2,j)
				522 continue
			endif
		else
			numProducts=1
			allocate(products(numProducts,numSpecies))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do 22 j=1,numSpecies
				reactants(2,j)=defectType1(j)
				reactants(1,j)=defectType2(j)
				products(1,j)=reactants(1,j)+reactants(2,j)
			22 continue
		endif
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: He, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of HeV+SIA=>HeSIA (should not be allowed here)
		!*******************************************************************************************
		!Vacancy+SIA annihilation - only the larger species remains
		if(products(1,2) .GE. products(1,3)) then
			products(1,2)=products(1,2)-products(1,3)
			products(1,3)=0
		endif
		if(products(1,2) .GE. products(1,4)) then
			products(1,2)=products(1,2)-products(1,4)
			products(1,4)=0
		endif
		if(products(1,2) .LT. products(1,3)) then
			products(1,3)=products(1,3)-products(1,2)
			products(1,2)=0
		endif
		if(products(1,2) .LT. products(1,4)) then
			products(1,4)=products(1,4)-products(1,2)
			products(1,2)=0
		endif
		
		!SIA+SIA clustering
		!two 1D clusters coming together to make a sessile cluster
		if(reactants(1,3) .GT. max3DInt .AND. reactants(2,3) .GT. max3DInt) then
			products(1,4)=products(1,3)
			products(1,3)=0
		endif
		
		!sessile SIA + mobile SIA = sessile SIA
		if(products(1,3) .NE. 0 .AND. products(1,4) .NE. 0) then
			products(1,4)=products(1,3)+products(1,4)
			products(1,3)=0
		endif
		
		if(products(1,1) .NE. 0) then

			!All He-SIA clusters are immobile
			if(products(1,3) .NE. 0) then
				products(1,4)=products(1,3)
				products(1,3)=0
			endif

		else
		
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(1,4) .NE. 0 .AND. products(1,4) .LE. max3DInt .AND. SIAPinToggle=='no') then
				products(1,3)=products(1,4)
				products(1,4)=0
			endif
		
		endif
		
		!Total Annihilation
		count=0
		do 25 j=1,numSpecies
			if(products(1,j)==0) then
				count=count+1
			endif
		25 continue
		if(count==numSpecies) then
			!we have completely annihilated the defects
			numProducts=0
			deallocate(products)
			allocate(products(numProducts, numSpecies))
		endif
		
		!findReactionInList points reactionCurrent at the reaction if it already exists. If not, reactionCurrent
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		
		reactionCurrent=>CascadeCurrent%reactionList(cell)
		call findReactionInListMultiple(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate
		call checkReactionLegality(numProducts, products, isLegal)
		
		if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, ClusterReactions(matNum,i))
		else
			reactionRate=0d0
		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=2
			reactionCurrent%numProducts=numProducts
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 23 j=1,numSpecies
				reactionCurrent%reactants(1,j)=reactants(1,j)
				reactionCurrent%reactants(2,j)=reactants(2,j)
				if(numProducts==1) then
					reactionCurrent%products(1,j)=products(1,j)
				else if(numProducts==2) then
					reactionCurrent%products(1,j)=products(1,j)
					reactionCurrent%products(2,j)=products(2,j)
				endif
			23 continue
			do 24 j=1,reactionCurrent%numReactants+reactionCurrent%numProducts
				reactionCurrent%cellNumber(j)=cell
				reactionCurrent%taskid(j)=myProc%taskid
			24 continue
			reactionCurrent%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			!update totalRate
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
	
	deallocate(reactants)
	deallocate(products)

	endif
	
10 continue

end subroutine

!>Subroutine add diffusion reactions - adds reactions to a reaction list representing diffusion 
!!between volume elements. 
!!
!! Inputs: initial and final cell numbers, initial and final processors (if diffusion occurs between
!!domains of multiple processors), direction (1-6), and defect type (we are adding all diffusion reactions associated with 
!!a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!!in this element and in the neighboring element
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.
!!
!!2015.04.06: Potential bug: if a volume element shares multiple faces with another grain, 
!!we are only creating one reaction for diffusion from the grain to the grain boundary.
!!Need to include diffusion direction in findReactionInListDiff

subroutine addDiffusionReactions(cell1, cell2, proc1, proc2, dir, defectType)
use mod_srscd_constants
use DerivedType
implicit none

integer cell1, cell2, proc1, proc2, defectType(numSpecies), dir, matNum
integer numReactants, numProducts, i, j, k, l, count, neighbor, neighborProc
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionCurrent, reactionPrev
double precision reactionRate

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell1)%material		!For now, assuming that cell1 and cell2 are in the same material.
endif

numReactants=1
numProducts=1
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 10 i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this diffusion reaction
	do 11 j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. DiffReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. DiffReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. DiffReactions(matNum,i)%max(j)) .OR. DiffReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	11 continue
	
	if(count==numSpecies) then	!this defect type is accepted for this diffusion reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do 12 j=1,numSpecies
			reactants(1,j)=defectType(j)
			products(1,j)=defectType(j)
		12 continue

		reactionCurrent=>reactionList(cell1)
		nullify(reactionPrev)
		
		!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)		
		call findReactionInListDiff(reactionCurrent, reactionPrev, reactants, cell1, cell2, proc1, proc2)

		reactionRate=findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, DiffReactions(matNum,i))

!		if(myProc%taskid==3) then
!			write(*,*) 'cell1', cell1, 'cell2', cell2, 'proc1', proc1, 'proc2', proc2, 'dir', dir, 'rate', reactionRate
!		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
			
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 l=1,numSpecies
				reactionCurrent%reactants(1,l)=reactants(1,l)
				reactionCurrent%products(1,l)=products(1,l)
			13 continue
			reactionCurrent%cellNumber(1)=cell1
			reactionCurrent%taskid(1)=proc1
			reactionCurrent%cellNumber(2)=cell2
			reactionCurrent%taskid(2)=proc2
			
			reactionCurrent%reactionRate=reactionRate
			
!			if(myProc%taskid==MASTER) then
!				write(*,*) reactionRate, cell1, cell2, proc1, proc2
!				if(associated(reactionCurrent%next)) then
!					write(*,*) 'reactionCurrent%next associated'
!				endif
!			endif
			
		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		endif
		
	endif
10 continue

end subroutine

!***************************************************************************************************
! Subroutine addDiffusionCoarseToFine
!
! Updates diffusion reaction rates associated wtih defectType from cell1 in coarse mesh to CasacadeCurrent
! fine mesh.
!
! Inputs: cell (cell number in coarse mesh), proc (of coarse mesh cell), CascadeCurrent (contains all
! 	cascade data), defectType(numSpecies)
!
! Outputs: creates new reaction in reactionList for diffusion from coarse mesh to fine mesh
!***************************************************************************************************

!>Subroutine add diffusion reactions from the coarse mesh to a fine mesh element (cascade) - adds reactions to a reaction list representing diffusion 
!!between volume elements. 
!!
!! Inputs: coarse mesh cell number, processor, cascade ID number, and defect type (we are adding all diffusion reactions associated with 
!!a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!!in this element and in the entire fine mesh. A modified reaction distance is used for this reaction rate.
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.

subroutine addDiffusionCoarseToFine(cell, proc, CascadeCurrent, defectType)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, proc, defectType(numSpecies), matNum
type(cascade), pointer :: CascadeCurrent
integer numReactants, numProducts, i, j, k, l, count, numDefectsFine
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionCurrent, reactionPrev
double precision reactionRate

interface
	integer function findNumDefectTotalFine(defectType, CascadeCurrent)
	use mod_srscd_constants
	integer defectType(numSpecies)
	type(cascade), pointer :: CascadeCurrent
	end function
end interface

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
endif
	
numReactants=1
numProducts=1
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 10 i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do 11 j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. DiffReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. DiffReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. DiffReactions(matNum,i)%max(j)) .OR. DiffReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	11 continue
	
	if(count==numSpecies) then	
		!this defect type is accepted for this dissociation reaction
		
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do 12 j=1,numSpecies
			reactants(1,j)=defectType(j)
			products(1,j)=defectType(j)
		12 continue

		reactionCurrent=>reactionList(cell)
		nullify(reactionPrev)
		
		!find the reaction in the reaction list. NOTE: -CascadeCurrent%cascadeID is used in place of cell2
		!to identify that this is a reaction from the coarse mesh to the fine mesh into this cascade.
		call findReactionInListDiff(reactionCurrent, reactionPrev, reactants, cell, -CascadeCurrent%cascadeID, proc, proc)
		
		!Find the total number of defects of type defectType in the fine mesh (all cells)
		numDefectsFine=findNumDefectTotalFine(defectType, CascadeCurrent)
		
		!Find the reaction rate for diffusion from coarse to fine mesh
		reactionRate=findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, DiffReactions(matNum,i))
		
		!Here, we update reactionList by either creating a new reaction or updating the current reaction
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
			
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 l=1,numSpecies
				reactionCurrent%reactants(1,l)=reactants(1,l)
				reactionCurrent%products(1,l)=products(1,l)
			13 continue
			reactionCurrent%cellNumber(1)=cell
			reactionCurrent%taskid(1)=proc
			
			!In coarse-to-fine reactions, the cascade ID number is stored as a negative value in the
			!cell number of the diffusion product (negative is signal that coarse-to-fine reaction
			!is occurring)
			reactionCurrent%cellNumber(2)=-CascadeCurrent%cascadeID
			reactionCurrent%taskid(2)=proc
			
			reactionCurrent%reactionRate=reactionRate
			
!			if(myProc%taskid==MASTER) then
!				write(*,*) reactionRate, cell1, cell2, proc1, proc2
!				if(associated(reactionCurrent%next)) then
!					write(*,*) 'reactionCurrent%next associated'
!				endif
!			endif
			
		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		endif
	endif
10 continue

end subroutine

!***************************************************************************************************
! Subroutine addDiffusionReactionsFine
!
! Updates diffusion reaction rates associated with defectType1 from cell1 to cell2 in fine mesh of cascade. 
!
! Inputs: cascade (cascade ID number), cell1, cell2, (volume element ID numbers), proc1, proc2 (processor numbers),
! dir (diffusion direction), defectType
! Outputs: updates multi-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!***************************************************************************************************

!>Subroutine add diffusion reactions (fine mesh) - adds reactions to a reaction list representing diffusion 
!!between volume elements inside a cascade mesh. 
!!
!!Inputs: cascade ID number, initial and final cell numbers, initial and final processors (if diffusion occurs between
!!domains of multiple processors), direction (1-6), and defect type (we are adding all diffusion reactions associated with 
!!a single defect type, which was a defect that changed in the previous Monte Carlo step)
!!
!!Structure of subroutine:
!!
!!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!!
!!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!!
!!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!!in this element and in the neighboring element (using the fine mesh element volume)
!!
!!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!!reaction is already present and if the reaction rate is nonzero.
!!
!!NOTE: includes possibility of diffusion from fine mesh to coarse mesh.

subroutine addDiffusionReactionsFine(cascadeID, cell1, cell2, proc1, proc2, dir, defectType)
use mod_srscd_constants
use DerivedType
implicit none

integer cascadeID, cell1, cell2, proc1, proc2, dir, defectType(numSpecies), matNum
type(cascade), pointer :: CascadeCurrent
integer numReactants, numProducts, i, j, k, l, count, neighbor, neighborProc
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionCurrent, reactionPrev
double precision reactionRate

!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this 
!subroutine
CascadeCurrent=>ActiveCascades

do 30 while(associated(CascadeCurrent))

	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	endif
	
	CascadeCurrent=>CascadeCurrent%next
30 continue

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material
endif

numReactants=1
numProducts=1
allocate(reactants(numReactants,numSpecies))
allocate(products(numProducts,numSpecies))

do 10 i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do 11 j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(matNum,i)%reactants(1,j) == 0) then
			count=count+1
		else if(defectType(j) .NE. 0 .AND. DiffReactions(matNum,i)%reactants(1,j) .NE. 0) then
			if(defectType(j) .GE. DiffReactions(matNum,i)%min(j)) then
				if((defectType(j) .LE. DiffReactions(matNum,i)%max(j)) .OR. DiffReactions(matNum,i)%max(j)==-1) then
					count=count+1
				endif
			endif
		endif
	11 continue

	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionCurrent at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionCurrent is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (diffusion)
		do 12 j=1,numSpecies
			reactants(1,j)=defectType(j)
			products(1,j)=defectType(j)
		12 continue

		!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)
		reactionCurrent=>CascadeCurrent%reactionList(cell1)
		nullify(reactionPrev)
		
		call findReactionInListDiff(reactionCurrent, reactionPrev, reactants, cell1, cell2, proc1, proc2)
		
		reactionRate=findReactionRateDiffFine(CascadeCurrent, defectType, cell1, proc1, cell2, proc2, dir, DiffReactions(matNum,i))
		
!		if(myProc%taskid==3) then
!			write(*,*) 'cell1', cell1, 'cell2', cell2, 'proc1', proc1, 'proc2', proc2, 'dir', dir, 'rate', reactionRate
!		endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionCurrent) .AND. reactionRate==0d0) then
			
			totalRate=totalRate-reactionCurrent%reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)-reactionCurrent%reactionRate
			
			!deleting reactionCurrent
			reactionPrev%next=>reactionCurrent%next
			deallocate(reactionCurrent%reactants)
			deallocate(reactionCurrent%products)
			deallocate(reactionCurrent%cellNumber)
			deallocate(reactionCurrent%taskid)
			nullify(reactionCurrent%next)
			deallocate(reactionCurrent)
			
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)+reactionRate
			
			!creating new reaction
			allocate(reactionCurrent)
			reactionCurrent%numReactants=1
			reactionCurrent%numProducts=1
			allocate(reactionCurrent%reactants(reactionCurrent%numReactants,numSpecies))
			allocate(reactionCurrent%products(reactionCurrent%numProducts,numSpecies))
			allocate(reactionCurrent%cellNumber(reactionCurrent%numReactants+reactionCurrent%numProducts))
			allocate(reactionCurrent%taskid(reactionCurrent%numReactants+reactionCurrent%numProducts))
			nullify(reactionCurrent%next)
			reactionPrev%next=>reactionCurrent
			do 13 l=1,numSpecies
				reactionCurrent%reactants(1,l)=reactants(1,l)
				reactionCurrent%products(1,l)=products(1,l)
			13 continue
			reactionCurrent%cellNumber(1)=cell1
			reactionCurrent%taskid(1)=proc1
			reactionCurrent%cellNumber(2)=cell2
			reactionCurrent%taskid(2)=proc2
			
			reactionCurrent%reactionRate=reactionRate
			
!			if(myProc%taskid==MASTER) then
!				write(*,*) reactionRate, cell1, cell2, proc1, proc2
!				if(associated(reactionCurrent%next)) then
!					write(*,*) 'reactionCurrent%next associated'
!				endif
!			endif
			
		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionCurrent) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionCurrent) .AND. reactionRate .NE. 0d0) then
			
			!update totalRate
			totalRate=totalRate-reactionCurrent%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)-reactionCurrent%reactionRate+reactionRate
			
			!update reaction rate
			reactionCurrent%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		endif

	endif
10 continue

end subroutine


!***************************************************************************************************
!TEMP example function used to inform the structure of other functions in this module.
!This will be deleted in the final version of the program.
!***************************************************************************************************

!> Function find reaction rate - finds reaction rate for implantation reaction (He, Frenkel pairs, cascades).
!!
!! Inputs: cell ID, reaction parameters (input from file)
!!
!! Calculates reaction rates according to hard-coded formulas. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! NOTE: for the case of non-uniform implantation, reaction rates can be dependent on the z-coordinate
!! of a volume element, using a local DPA rate and helium implantation rate read in from an 
!! input file.

double precision function findReactionRate(cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell
type(reactionParameters) :: reactionParameter

integer numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision Diff, Eb, volume, DPARateLocal, HeImplantRateLocal, zCoord
integer n, numClusters

if(reactionParameter%functionType==1) then
	!***********************************************************************************************
	!This is an example of a function reaction rate. It is not operational currently because all clustering,
	!etc. type reactions are taken care of in other find reaction rate functions
	!***********************************************************************************************
	!calculate Diff and Eb using functions in MaterialInput
	!constants omega, kboltzmann, temperature should be global. Information on defect numbers is global.
	
	!findReactionRate=omega*dble(n)**(4d0/3d0)*Diff*exp(-Eb/(kboltzmann*temperature))*numClusters
	
	write(*,*) 'Error should not be choosing function type 1 in findReactionRate()'

else if(reactionParameter%functionType==10) then	!Frenkel pair implantation
	
	!This is the rate of Frenkel pair implantation events inside a cell with given volume (easy to calculate)
	volume=(myMesh(cell)%length)**3d0
	
	if(implantDist=='Uniform') then
		findReactionRate=volume*DPARate/atomsize
	else if(implantDist=='NonUniform') then
		
		zCoord=myMesh(cell)%coordinates(3)
		DPARateLocal=findDPARateLocal(zCoord)
		findReactionRate=volume*DPARateLocal/atomsize
		
	else
		write(*,*) 'Error implant distribution not recognized'
	endif
	
else if(reactionParameter%functionType==11) then	!cascade implantation
	
	!This is the rate of cascade implantation events inside a cell with a given volume
	volume=(myMesh(cell)%length)**3d0
	
	if(implantDist=='Uniform') then
		findReactionRate=volume*DPARate/(numDisplacedAtoms*atomsize)
	else if(implantDist=='NonUniform') then
		
		zCoord=myMesh(cell)%coordinates(3)
		DPARateLocal=findDPARateLocal(zCoord)
		findReactionRate=volume*DPARateLocal/(numDisplacedAtoms*atomsize)
		
	else
		write(*,*) 'Error implant distribution not recognized'
	endif
	
else if(ReactionParameter%functionType==12) then	!He implantation
	
	!This is the rate of Helium implantation events inside a cell with given volume
	volume=(myMesh(cell)%length)**3d0
	
	if(implantDist=='Uniform') then
		findReactionRate=volume*HeDPARatio*DPARate/atomsize
	else if(implantDist=='NonUniform') then
		
		zCoord=myMesh(cell)%coordinates(3)
		HeImplantRateLocal=findHeImplantRateLocal(zCoord)
		findReactionRate=volume*HeImplantRateLocal/atomsize
		
	else
		write(*,*) 'Error implant distribution not recognized'
	endif
	
else if(ReactionParameter%functionType==13) then	!FP implantation disallowed in grain boundaries

	findReactionRate=0d0
	
else
	write(*,*) 'error function type', ReactionParameter%functionType
endif

end function

!***************************************************************************************************
!
! This subroutine finds the reaction rate for implantation reactions (He only for now) within the
! fine mesh. 
!
! Inputs: cell number in fine mesh, reaction parameters
! Outputs: reaction rate for He implantation or error message
!
!***************************************************************************************************

!> Function find reaction rate fine - finds reaction rate for implantation reaction (Helium or Frenkel pairs only) in a fine mesh.
!! Cascades cannot be implanted inside other cascades
!!
!! Inputs: cell ID, reaction parameters (input from file)
!!
!! Calculates reaction rates according to hard-coded formulas. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.

double precision function findReactionRateFine(cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell
type(reactionParameters) :: reactionParameter

integer numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision Diff, Eb, volume, HeImplantRateLocal, zCoord
integer n, numClusters

if(ReactionParameter%functionType==12) then	!He implantation
	
	!This is the rate of Helium implantation events inside a cell with given volume
	volume=cascadeElementVol
	
	if(implantDist=='Uniform') then
		findReactionRateFine=volume*HeDPARatio*DPARate/atomsize
	else if(implantDist=='NonUniform') then
		
		zCoord=myMesh(cell)%coordinates(3)
		HeImplantRateLocal=findHeImplantRateLocal(zCoord)
		findReactionRateFine=volume*HeImplantRateLocal/atomsize
		
	else
		write(*,*) 'Error implant distribution not recognized'
	endif
	
else
	write(*,*) 'error function type'
endif

end function

!***************************************************************************************************
!This function will return the reaction rate of a function of type given by reactionParameter and
!with reactant given by defectType in cell of myMesh.
!
!reactionParameter%functionType tells the function which functional form to look up to find the 
!reaction rate. The information in this function is mostly hard-coded material-specific information.
!This was seen as inevitable for this program, and this is where the majority of changes should
!occur if we choose to change material parameters.
!***************************************************************************************************

!> Function find reaction rate impurity - finds reaction rate for trapping of SIA loops by impurities (Carbon).
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for SIA loop trapping according to hard-coded formulas. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! In this case, glissile loops become sessile when they interact with impurities. Have the option
!! for 3-D or 1-D diffusion of loops.
!!
!! We have also added the reaction rate for sessile clusters to become mobile again using a 
!! binding energy in this subroutine. The binding energy is read in from an input file.

double precision function findReactionRateImpurity(defectType, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), singleDefectType(numSpecies), num, size
integer productType(numSpecies), i
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff
double precision findDiffusivity, findBinding, Eb
integer findDefectSize, findNumDefect, matNum, grainNum

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell)%material
else
	matNum=myMesh(cell)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell)%material
endif

if(reactionParameter%functionType==3) then

	num=findNumDefect(defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)		!diffusivity of clusters of this type
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar+omega*(dble(size)**(1d0/3d0)))*&
			(Diff)*dble(num)*impurityDensity

else if(reactionParameter%functionType==4) then
	Diff=findDiffusivity(matNum,defectType)
	num=findNumDefect(defectType,cell)
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
		Diff*dble(num)*(impurityDensity**2d0)
		
else if(reactionParameter%functiontype==5) then
	
	!To Do: create reaction rate for sessile clusters becoming mobile again due to a binding energy.
	do 10 i=1,numSpecies
		if(i==3) then
			productType(i)=defectType(4)
		else if(i==4) then
			productType(i)=defectType(3)
		else
			productType(i)=defectType(i)
		endif
	10 continue
	
	num=findNumDefect(defectType,cell)
	Eb=findBinding(matNum,defectType,productType)
	Diff=findDiffusivity(matNum,productType)
	
	reactionRate=omega*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)

else
	write(*,*) 'error impurity trapping function type only admits 3 or 4'
	reactionRate=0d0
endif

findReactionRateImpurity=reactionRate

end function

!***************************************************************************************************
!function findReactionRateImpurityFine
!
!Finds reaction rates for defect clustering with impurities (used here for mobile SIA+impurity -> 
!sessile SIA reactions).
!
!Inputs: defectType(numSpecies), cell, reactionParameter (read from file), CascadeCurrent
!Output: reactionRate
!***************************************************************************************************

!> Function find reaction rate impurity fine - finds reaction rate for trapping of SIA loops by impurities (Carbon) inside a fine mesh (Cascade).
!!
!! Inputs: cascade derived type, cell ID, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for SIA loop trapping according to hard-coded formulas inside a fine mesh (cascade). Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! The cascade derived type must be passed into this subroutine in order to find the number of defects of this type.
!!
!! In this case, glissile loops become sessile when they interact with impurities. Have the option
!! for 3-D or 1-D diffusion of loops.
!!
!! We have also added the reaction rate for sessile clusters to become mobile again using a 
!! binding energy in this subroutine. The binding energy is read in from an input file.

double precision function findReactionRateImpurityFine(CascadeCurrent, defectType, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), singleDefectType(numSpecies), num, size
integer productType(numSpecies), i
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff
type(cascade), pointer :: cascadeCurrent
double precision findDiffusivity, findBinding, Eb
integer findDefectSize, findNumDefect, matNum, grainNum

interface 
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_srscd_constants
	type(cascade), pointer :: CascadeCurrent
	integer defectType(numSpecies), cell
	end function
end interface

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
endif

if(reactionParameter%functionType==3) then

	num=findNumDefectFine(CascadeCurrent,defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)							!diffusivity of clusters of this type
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar+omega*(dble(size)**(1d0/3d0)))*&
			(Diff)*dble(num)*impurityDensity

else if(reactionParameter%functionType==4) then

	Diff=findDiffusivity(matNum,defectType)
	num=findNumDefectFine(CascadeCurrent,defectType,cell)
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
		Diff*dble(num)*(impurityDensity**2d0)

else if(reactionParameter%functiontype==5) then
	
	!To Do: create reaction rate for sessile clusters becoming mobile again due to a binding energy.
	do 10 i=1,numSpecies
		if(i==3) then
			productType(i)=defectType(4)
		else if(i==4) then
			productType(i)=defectType(3)
		else
			productType(i)=defectType(i)
		endif
	10 continue
	
	num=findNumDefectFine(CascadeCurrent,defectType,cell)
	Eb=findBinding(matNum,defectType,productType)
	Diff=findDiffusivity(matNum,productType)
	
	reactionRate=omega*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	
else
	write(*,*) 'error impurity trapping function type only admits 3 or 4'
	reactionRate=0d0
endif

findReactionRateImpurityFine=reactionRate

end function

!> Function find reaction rate dissociation - finds reaction rate for point defects to dissociate from clusters
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies), dissociating defect type
!!
!! Calculates reaction rates for dissociation of a point defect from a cluster. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.

double precision function findReactionRateDissoc(defectType, products, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), products(2,numSpecies), size, num
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff, Eb
double precision findDiffusivity, findBinding
integer findNumDefect, findDefectSize, matNum, grainNum

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell)%material
else
	matNum=myMesh(cell)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell)%material
endif

!dissociation
if(reactionParameter%functionType==1) then
	!dissocation reactions

	Diff=findDiffusivity(matNum,products(2,:))			!diffusivity of the defect dissociating from the cluster

	num=findNumDefect(defectType,cell)			!number of clusters

	size=findDefectSize(defectType)				!Hard-coded, rules for determining which species governs the defect size

	Eb=findBinding(matNum,defectType,products(2,:))

	reactionRate=omega*dble(size)**(4d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	
else
	write(*,*) 'error dissociation function type only admits 1'
	reactionRate=0d0
endif

findReactionRateDissoc=reactionRate

end function

!***************************************************************************************************
!function findReactionRateDissocFine
!
!finds reaction rate for dissociation reactions inside fine mesh.
!
!Inputs: CascadeCurrent, defectType(numSpecies), products(:,:), cell, reactionParameter (read from file)
!Output: reactionRate
!***************************************************************************************************

!> Function find reaction rate dissociation - finds reaction rate for point defects to dissociate from clusters in the fine mesh (cascade)
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies), dissociating defect type, cascade derived type
!!
!! Calculates reaction rates for dissociation of a point defect from a cluster. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! Cascade derived type is used to find the number of defects in the volume element (it is required in a later subroutine)

double precision function findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), products(2,numSpecies), size, num
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff, Eb
type(cascade), pointer :: CascadeCurrent
type(defect), pointer :: defectTemp
double precision findDiffusivity, findBinding
integer findDefectSize, matNum, grainNum

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_srscd_constants
	integer cell, defectType(numSpecies)
	type(cascade), pointer :: CascadeCurrent
	end function
end interface

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
endif

!dissociation
if(reactionParameter%functionType==1) then

	!dissocation reactions
	Diff=findDiffusivity(matNum,products(2,:))								!diffusivity of the defect dissociating from the cluster

	defectTemp=>CascadeCurrent%localDefects(cell)

	num=findNumDefectFine(CascadeCurrent, defectType,cell)			!number of clusters

	size=findDefectSize(defectType)									!Hard-coded, rules for determining which species governs the defect size

	Eb=findBinding(matNum,defectType,products(2,:))						!binding energy of single defect to cluster

	
	reactionRate=omega*dble(size)**(4d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
else
	write(*,*) 'error dissociation function type only admits 1'
	reactionRate=0d0
endif

findReactionRateDissocFine=reactionRate

end function

!> Function find reaction rate sink - finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations)
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for capture of a defect by a sink. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.

double precision function findReactionRateSink(defectType, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), num
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff
double precision findDiffusivity
integer findNumDefect, matNum, grainNum

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell)%material
else
	matNum=myMesh(cell)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell)%material
endif

!sink reaction function type=3
if(reactionParameter%functionType==3) then
	
	num=findNumDefect(defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)		!diffusivity of clusters of this type
	
	if(defectType(3) .NE. 0) then !interstitial defect
		reactionRate=Zint*dislocationDensity*diff*dble(num)
	else
		reactionRate=dislocationDensity*diff*dble(num)
	endif
else
	write(*,*) 'error sink function type only admits 3'
	reactionRate=0d0
endif

findReactionRateSink=reactionRate

end function


!***************************************************************************************************
!function findReactionRateSinkFine
!
!finds reaction rate for trapping of defects at sinks inside fine mesh.
!
!Inputs: CascadeCurrent, defectType(numSpecies), cell, reactionParameter (read from file)
!Output: reactionRate
!***************************************************************************************************

!> Function find reaction rate sink fine - finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations) in the fine mesh (cascade)
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies), cascade derived type
!!
!! Calculates reaction rates for capture of a defect by a sink in the fine mesh (in a cascade). Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.

double precision function findReactionRateSinkFine(CascadeCurrent, defectType, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), num, matNum, grainNum
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff
type(cascade), pointer :: CascadeCurrent
double precision findDiffusivity

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_srscd_constants
	integer cell, defectType(numSpecies)
	type(cascade), pointer :: CascadeCurrent
	end function
end interface

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
endif

!sink reaction function type=3
if(reactionParameter%functionType==3) then
	
	num=findNumDefectFine(CascadeCurrent, defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)							!diffusivity of clusters of this type
	
	if(defectType(3) .NE. 0) then !interstitial defect
		reactionRate=Zint*dislocationDensity*diff*dble(num)
	else
		reactionRate=dislocationDensity*diff*dble(num)
	endif
else
	write(*,*) 'error sink function type only admits 3'
	reactionRate=0d0
endif

findReactionRateSinkFine=reactionRate

end function

!> Function find reaction rate multiple - finds reaction rate for defect clustering
!!
!! Inputs: cell ID, reaction parameters (input from file), 2 defect types (array size numSpecies)
!!
!! Calculates reaction rates for two defects to cluster with each other. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! Several different hard-coded reaction rates are used here, according to the various 
!! reaction rates for clustering between defects depending on their geometry and diffusivity

double precision function findReactionRateMultiple(defectType1, defectType2, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), i, count
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff1, Diff2, size1, size2, num1, num2, vol, Ztemp
integer findNumDefect, findDefectSize, matNum, grainNum
double precision findDiffusivity
double precision rad1, rad2, area

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell)%material
else
	matNum=myMesh(cell)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell)%material
endif

size1=findDefectSize(defectType1)
size2=findDefectSize(defectType2)
Diff1=findDiffusivity(matNum,defectType1)
Diff2=findDiffusivity(matNum,defectType2)
num1=findNumDefect(defectType1,cell)
num2=findNumDefect(defectType2,cell)

!Hard-coded: sink strength bias for SIAs
if(defectType1(3) .NE. 0 .OR. defectType1(4) .NE. 0) then
	Ztemp=1.2d0
else if(defectType2(3) .NE. 0 .OR. defectType2(4) .NE. 0) then
	Ztemp=1.2d0
else
	Ztemp=1.0d0
endif

count=0
do 10 i=1,numSpecies
	if(defectType1(i)==defectType2(i)) then
		count=count+1
	endif
10 continue
if(count==numSpecies) then
	!we have two defects of the same type, have to modify the defect numbers for a defect to combine with itself
	num2=num2-1
endif

vol=myMesh(cell)%volume
area=(myMesh(cell)%length)**2d0	!assuming square elements in grain boundary

if(reactionParameter%functionType==5) then
	
	if(size1+size2 .GT. 4) then	!He-He reactions
		reactionRate=0d0		!hard-coded: limit size of interstitial He clusters to max 4 (no parameters available for larger clusters)
	else
		reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*&
			(Diff1+Diff2)*dble(num1)*dble(num2)*atomsize/vol
	endif
	
else if(reactionParameter%functionType==6) then	!spherical clusters other than He-He clusters

	reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				 *atomsize/vol

else if(reactionParameter%functionType==7) then	!For HeV clusters. NOTE that this is the same as function type 6 (identical for now)

	reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				 *atomsize/vol

else if(reactionParameter%functionType==8) then

	if(defectType1(3) .GT. max3DInt .OR. defectType1(4) .GT. max3DInt) then
		
		!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters
		!in order to have the correct reaction rate.
		
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefect(defectType2,cell)
		num2=findNumDefect(defectType1,cell)
		
	else if(defectType2(3) .GT. max3DInt .OR. defectType2(4) .GT. max3DInt) then
		
		!do nothing, this is the same as the default at the beginning of this subroutine
		
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	endif

	reactionRate=Zint*(omegastar+(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0)))*&
		Diff1*num1*num2*atomsize/vol+&
		(Zint*(omegastar1D+omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
		Diff2*dble(num2)*dble(num1)**(2d0)*(atomsize/vol)**(2d0)
		
else if(reactionParameter%functionType==9) then
	
	reactionRate=(Zint*(omegastar1D+omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0))))**4d0*&
		(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomsize/vol)**(2d0)

else if(reactionParameter%functionType==10) then	!Diffusion rate for 2D defect recombination on the grain boundary

	rad1=(3d0*size1*atomsize/(4d0*pi))**(1d0/3d0)
	rad2=(3d0*size2*atomsize/(4d0*pi))**(1d0/3d0)
	
	if(num1==0 .OR. num2==0) then
	
		reactionRate=0d0
		
	else
	
		reactionRate=	(4d0*pi*dble(num1)*dble(num2)/area)*&
						(Diff1/(2d0*dlog(dsqrt(area/(pi*num2))*(1d0/(rad1+rad2)))-1d0)+&
						 Diff2/(2d0*dlog(dsqrt(area/(pi*num1))*(1d0/(rad1+rad2)))-1d0))
	endif
						 
!	write(*,*) 'defect1', defectType1, 'num', num1
!	write(*,*) 'defect2', defectType2, 'num', num2
!	write(*,*) 'rate', reactionRate

	if(reactionRate .LT. 0d0) then
!		write(*,*) 'Error negative reaction rate', reactionRate
!		write(*,*) 'defect1', defectType1, 'num', num1
!		write(*,*) 'defect2', defectType2, 'num', num2
		
		reactionRate=0d0
	endif

else
	write(*,*) 'error clustering function type only admits 5-9'
	reactionRate=0d0
endif

findReactionRateMultiple=reactionRate
end function


!***************************************************************************************************
!function findReactionRateMultipleFine
!
!finds reaction rate for defect clustering (which takes two defects to occur).
!
!Inputs: CascadeCurrent, defectType1(numSpecies), defectType2(numSpecies), cell, reactionParameter (read from file)
!Output: reactionRate
!***************************************************************************************************

!> Function find reaction rate multiple fine - finds reaction rate for defect clustering in the fine mesh (Cascade)
!!
!! Inputs: cell ID, reaction parameters (input from file), 2 defect types (array size numSpecies), cascade derived type
!!
!! Calculates reaction rates for two defects to cluster with each other inside the fine mesh (in a cascade). Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! Several different hard-coded reaction rates are used here, according to the various 
!! reaction rates for clustering between defects depending on their geometry and diffusivity

double precision function findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), i, count
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff1, Diff2, size1, size2, num1, num2, vol, Ztemp
!integer findNumDefectFine
type(cascade), pointer :: CascadeCurrent
double precision findDefectSize, findDiffusivity
integer matNum, grainNum

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_srscd_constants
	integer cell, defectType(numSpecies)
	type(cascade), pointer :: CascadeCurrent
	end function
end interface

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
endif

size1=findDefectSize(defectType1)
size2=findDefectSize(defectType2)
Diff1=findDiffusivity(matNum,defectType1)
Diff2=findDiffusivity(matNum,defectType2)
num1=findNumDefectFine(CascadeCurrent, defectType1,cell)
num2=findNumDefectFine(CascadeCurrent, defectType2,cell)
vol=cascadeElementVol

!Hard-coded: sink strength bias for SIAs
if(defectType1(3) .NE. 0 .OR. defectType1(4) .NE. 0) then
	Ztemp=1.2d0
else if(defectType2(3) .NE. 0 .OR. defectType2(4) .NE. 0) then
	Ztemp=1.2d0
else
	Ztemp=1.0d0
endif

!Adjust the number of defects in the reaction if both reactants are of the same type

count=0

do 10 i=1,numSpecies
	if(defectType1(i)==defectType2(i)) then
		count=count+1
	endif
10 continue

if(count==numSpecies) then
	!we have two defects of the same type, have to modify the defect numbers for a defect to combine with itself
	num2=num2-1
endif

!list of clustering reaction functional forms

if(reactionParameter%functionType==5) then
	
	if(size1+size2 .GT. 4) then	!He-He reactions
		reactionRate=0d0		!hard-coded: limit size of interstitial He clusters to max 4 (no parameters available for larger clusters)
	else
		reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*&
			(Diff1+Diff2)*dble(num1)*dble(num2)*atomsize/vol
	endif
else if(reactionParameter%functionType==6) then	!spherical clusters other than He-He clusters

	reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				 *atomsize/vol

else if(reactionParameter%functionType==7) then	!For HeV clusters. NOTE that this is the same as function type 6 (identical for now)

	reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				 *atomsize/vol

else if(reactionParameter%functionType==8) then

	if(defectType1(3) .GT. max3DInt .OR. defectType1(4) .GT. max3DInt) then
		
		!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters
		!in order to have the correct reaction rate.
		
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefectFine(CascadeCurrent, defectType2,cell)
		num2=findNumDefectFine(CascadeCurrent, defectType1,cell)
		
	else if(defectType2(3) .GT. max3DInt .OR. defectType2(4) .GT. max3DInt) then
		
		!do nothing, this is the same as the default at the beginning of this subroutine
		
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	endif

	reactionRate=Zint*(omegastar+(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0)))*&
		Diff1*num1*num2*atomsize/vol+&
		(Zint*(omegastar1D+omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
		Diff2*dble(num2)*dble(num1)**(2d0)*(atomsize/vol)**(2d0)
		
else if(reactionParameter%functionType==9) then
	
	reactionRate=(Zint*(omegastar1D+omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0))))**4d0*&
		(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomsize/vol)**(2d0)
else
	write(*,*) 'error clustering function type only admits 5-9'
	reactionRate=0d0
endif

findReactionRateMultipleFine=reactionRate
end function

!> Function find reaction rate diffusion - finds reaction rate for defect diffusion between elements
!!
!! Inputs: cell IDs, processor IDs, diffusion direction, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative).

double precision function findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell1, proc1, cell2, proc2, defectType(numSpecies), num1, num2, dir
type(ReactionParameters) :: reactionParameter
double precision Diff, area1, area2, areaShared, lengthShared, Vol1, Vol2, length1, length2, reactionRate
double precision strainE1, strainE2
integer findNumDefect, findNumDefectBoundary
double precision findStrainEnergy, findStrainEnergyBoundary
double precision findDiffusivity
integer matNum, grainNum, matNeighbor, size
double precision Eb, findBinding
double precision alpha

!If we are in a polycrystalline simulation, reactionParameter will already be
!chosen using matNum=1, and inside this subroutine we use matNum to indicate which
!grain we are inside (rather than which material type we are inside).

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell1)%material
else
	matNum=myMesh(cell1)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell1)%material
endif

if(reactionParameter%functionType==2) then

	Diff=findDiffusivity(matNum,defectType)		!function in MaterialInput
	length1=myMesh(cell1)%length
	area1=length1**2d0
	Vol1=myMesh(cell1)%volume
	num1=findNumDefect(defectType,cell1)
	
	strainE1=findStrainEnergy(defectType, cell1)
	if(strainE1 .NE. 0d0) then
		write(*,*) 'defectType', defectType
		write(*,*) 'strainE1', strainE1
		read(*,*)
	endif
	
	if(proc2==-1) then	!free surface (NOTE: not taking into account strain-assisted diffusion here)
		areaShared=area1
		reactionRate=Diff*areaShared*(dble(num1)/Vol1)/length1
		if(reactionRate .GT. 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		endif
	else	!cell-to-cell diffusion
		!Find various parameters needed for reaction rate
		
		if(proc2==proc1) then
			matNeighbor=myMesh(cell2)%material
			length2=myMesh(cell2)%length
			Vol2=myMesh(cell2)%volume
			strainE2=findStrainEnergy(defectType, cell2)
		else
			matNeighbor=myBoundary(dir,cell2)%material
			length2=myBoundary(dir,cell2)%length
			Vol2=myBoundary(dir,cell2)%volume
			strainE2=findStrainEnergyBoundary(defectType, cell2, dir)
		endif
		
		area2=length2**2d0
		
		!The area of the shared interface between the volume elements is the minimum of the two volume 
		!elements' face areas (assuming cubic elements, nonuniform)
		if(area1 .GT. area2) then
			areaShared=area2
		else
			areaShared=area1
		endif
		
		if(grainNum==matNeighbor) then
			if(proc2==proc1) then
				num2=findNumDefect(defectType,cell2)
			else
				num2=findNumDefectBoundary(defectType,cell2,dir)
			endif
			
			alpha=1d0
			
		else
			
			if(numMaterials .GT. 1) then
				!Sink efficiencies (Default set at 1)
				if(defectType(2) .NE. 0d0) then
					alpha=alpha_v
				else if(defectType(3) .NE. 0d0) then
					alpha=alpha_i
				else
					alpha=1d0
				endif
			else
				alpha=1d0
			endif

			num2=0	!If we are diffusing between two material types, assume perfect sink
		endif
		
		if(strainField=='yes') then
			if(strainE1 .NE. 0d0) then
				write(*,*) 'Strain E1', strainE1, 'Strain E2', strainE2
				write(*,*) 'defect type', defectType
				read(*,*)
			endif
			reactionRate=alpha*Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2+&
				dble(num1)*(strainE1-strainE2)/(Vol1*kboltzmann*temperature))/(length1/2d0 + length2/2d0)
		else
			reactionRate=alpha*Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(length1/2d0 + length2/2d0)
		endif
		
		if(reactionRate .GT. 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		endif
	endif
elseif(reactionParameter%functionType==3) then	!2D diffusion on a plane (rather than 3D diffusion in a volume)
	Diff=findDiffusivity(matNum,defectType)
	length1=myMesh(cell1)%length
	area1=length1**2d0
	num1=findNumDefect(defectType,cell1)
	
	if(proc2==-1) then	!free surface
		lengthShared=length1
		reactionRate=Diff*(lengthShared/length1)*(dble(num1)/area1)
		if(reactionRate .GT. 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		endif
	else	!cell-to-cell diffusion
		!Find various parameters needed for reaction rate
		
		if(proc2==proc1) then
			matNeighbor=myMesh(cell2)%material
			length2=myMesh(cell2)%length
		else
			matNeighbor=myBoundary(dir,cell2)%material
			length2=myBoundary(dir,cell2)%length
		endif
		
		if(grainNum==matNeighbor) then
		
			area2=length2**2d0
			
			!The length of the shared interface between the area elements is the minimum of the two 
			!elements' side lengths (assuming cubic elements, nonuniform)
			if(length1 .GT. length2) then
				lengthShared=length2
			else
				lengthShared=length1
			endif
			
			if(proc2==proc1) then
				num2=findNumDefect(defectType,cell2)
			else
				num2=findNumDefectBoundary(defectType,cell2,dir)
			endif
			
			reactionRate=Diff*(lengthShared/(length1/2d0+length2/2d0))*(dble(num1)/area1-dble(num2)/area2)
			
			if(reactionRate .GT. 0d0) then
				findReactionRateDiff=reactionRate
			else
				findReactionRateDiff=0d0
			endif
		
		else
		
			findReactionRateDiff = 0d0	!Don't let 2D diffusion occur between different material types
			
		endif
		
	endif

elseif(reactionParameter%functionType==4) then	!Dissociation from grain boundary into bulk volume element
	!(treated as a diffusion reaction because we are moving between two volume elements, but the rate is given by a dissociation rate)
	
	if(proc2==proc1) then
		matNeighbor=myMesh(cell2)%material
	else
		matNeighbor=myBoundary(dir,cell2)%material
	endif
	
	if(grainNum==matNeighbor) then
	
		findReactionRateDiff=0d0		!no dissociation from grain boundary to itself
		
	else
		
		!5/31/2015: Diff should be the diffusivity of the defect type IN THE BULK, not on the GB. Therefore,
		!use matNeighbor as input for findDiffusivity instead of matNum.
		
		Diff=findDiffusivity(matNeighbor, defectType)		!diffusivity of the defect dissociating from the GB (in the bulk)
		!Diff=findDiffusivity(matNum,defectType)			!diffusivity of the defect dissociating from the cluster
	
		num1=findNumDefect(defectType,cell1)			!number of clusters
		
		size=1										!not breaking up a cluster, but releasing from grain boundary
		
		Eb=findBinding(matNum,defectType,defectType)	!binding energy of defect to grain boundary
	
		reactionRate=omega*dble(size)**(4d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num1)
		
		if(reactionRate .GT. 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		endif
		
	endif

else
	write(*,*) 'error find reaction rate diffusion'
	findReactionRateDiff=0d0
endif

end function

!***************************************************************************************************
! function findReactionRateCoarseToFine
!
! finds reaction rate for defect diffusion from coarse mesh to fine mesh. Assumes diffusion into 
! entire fine mesh, not a single cell in the fine mesh.
! 
! Inputs: defectType(numSpecies): identification of the defect type, cell: coarse mesh volume element
! number, proc: coarse and fine mesh processor number, numDefectsFine: total number of defects of this
! type in the fine mesh (already calculated), reactionParameter: read from file, contains info on 
! this diffusion reaction.
!***************************************************************************************************

!> Function find reaction rate diffusion coarse to fine - finds reaction rate for defect diffusion between a coarse mesh element and a fine (cascade) mesh
!!
!! Inputs: cell ID, processor ID, reaction parameters (input from file), defect type (array size numSpecies), number of defects of this type in the fine mesh already
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative), with a modified diffusion 
!! distance according to the formula in Dunn et al. (Computational Materials Science 2015)

double precision function findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer defectType(numSpecies), cell, proc, numDefectsFine, num1, num2, matNum, grainNum
type(ReactionParameters) :: reactionParameter

double precision Diff, areaShared, Vol1, Vol2, length1, reactionRate, coarseToFineLength
integer findNumDefect, findNumDefectBoundary
double precision findDiffusivity

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(cell)%material
else
	matNum=myMesh(cell)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(cell)%material
endif

if(reactionParameter%functionType==2) then
	
	!coarse-to-fine mesh diffusion
	!Find various parameters needed for reaction rate
	
	!Length of coarse mesh element
	length1=myMesh(cell)%length

	!*************************************************************************************************
	!Effective length used for diffusion from coarse to fine mesh (assuming fine mesh randomly located
	!within the coarse mesh element, see supporting documents for derivation)
	!
	!NOTE: assuming cubic fine mesh (numxCascade=numyCascade=numzCascade). Can re-derive for non-cubic
	!fine meshes.
	!*************************************************************************************************

	CoarseToFineLength=(length1-numxCascade*fineLength)&
		/(dlog(length1**2d0/(numxCascade*fineLength)**2d0))
	
	!Diffusivity taken from MaterialInput
	Diff=findDiffusivity(matNum,defectType)		
	
	!Volume of coarse mesh element
	Vol1=myMesh(cell)%volume
	
	!Number of defects in coarse mesh element and fine mesh (entire mesh)
	num1=findNumDefect(defectType,cell)
	num2=numDefectsFine
	
	!Area of surface of fine mesh (entire surface)
	areaShared=2d0*(numxCascade*numyCascade+numxCascade*numzCascade+numyCascade*numzCascade)*fineLength**2d0	
	
	!Volume of fine mesh (entire volume)
	Vol2=numxCascade*numyCascade*numzCascade*cascadeElementVol
	
	!************************************************
	!Reaction rate formula using the effective length
	!************************************************
	
	reactionRate=Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(CoarseToFineLength)

	if(reactionRate .GT. 0d0) then
		findReactionRateCoarseToFine=reactionRate
	else
		findReactionRateCoarseToFine=0d0
	endif

else
	write(*,*) 'error find reaction rate diffusion coarse to fine'
	findReactionRateCoarseToFine=0d0
endif

end function

!***************************************************************************************************
!function findReactionRateDiffFine
!
!finds reaction rate for defect diffusion from cell to cell (inside fine mesh). Includes reaction
!rate for diffusion from fine mesh to coarse mesh (if applicable).
!
!Inputs: CascadeCurrent, defectType(numSpecies), cell1, cell2, proc1, proc2, dir, reactionParameter (read from file)
!Output: reactionRate
!***************************************************************************************************

!> Function find reaction rate diffusion fine - finds reaction rate for defect diffusion between elements in the fine mesh (inside a cascade)
!!
!! Inputs: cell IDs, processor IDs, diffusion direction, reaction parameters (input from file), defect type (array size numSpecies), cascade derived type
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative). Includes rate for diffusion
!! from fine mesh to coarse mesh.

double precision function findReactionRateDiffFine(CascadeCurrent, defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
use mod_srscd_constants
use DerivedType
implicit none

integer cell1, proc1, cell2, proc2, defectType(numSpecies), num1, num2, dir, coarseCell
type(ReactionParameters) :: reactionParameter
type(cascade), pointer :: CascadeCurrent
double precision Diff, area1, area2, areaShared, Vol1, Vol2, length1, length2, reactionRate, coarseLength
double precision fineToCoarseLength, coarseVolume
integer findNumDefect, matNum, grainNum
double precision findDiffusivity

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_srscd_constants
	type(cascade), pointer :: CascadeCurrent
	integer cell, defectType(numSpecies)
	end function
end interface

if(polycrystal=='yes') then
	matNum=1
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material	!not worrying about diffusion between multiple material types right now
	grainNum=myMesh(CascadeCurrent%cellNumber)%material
endif

if(reactionParameter%functionType==2) then

	Diff=findDiffusivity(matNum,defectType)		!function in MaterialInput
	
	length1=fineLength
	area1=length1**2d0
	Vol1=cascadeElementVol
	
	num1=findNumDefectFine(CascadeCurrent, defectType,cell1)

	if(proc2==-1) then	
	
		write(*,*) 'error free surface diffusion from inside fine mesh'
	
		!free surface
	
		areaShared=area1
		reactionRate=Diff*areaShared*(dble(num1)/Vol1)/length1
		if(reactionRate .GT. 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		endif
	
	else if(cell2==0) then 
	
	!diffusion from fine mesh to coarse mesh

		!Find information on defects and volume element size in coarse mesh element containing this cascade

		coarseCell=CascadeCurrent%cellNumber
		
		coarseLength=myMesh(coarseCell)%length
		
		coarseVolume=myMesh(coarseCell)%volume
			
		num2=findNumDefect(defectType, coarseCell)

		!***********************************************************************
		!average diffusion distance from cascade element to coarse mesh element
		!(assuming cubic cascade mesh and coarse mesh element)
		!
		!See supporting documentation for derivation of this formula
		!
		!Note: assuming cubic cascade meshes here (used numxCascade for all
		!diffusion directions)
		!***********************************************************************
		
		fineToCoarseLength=(coarseLength-numxCascade*fineLength)&
			/(dlog((coarseLength-(numxCascade-1)*fineLength)**2d0/(fineLength**2d0)))
			
		!Reaction Rate using the average diffusion distance
		
		reactionRate=Diff*(fineLength**2d0)*&
			(dble(num1)/Vol1-dble(num2)/coarseVolume)/(fineToCoarseLength)
		
		if(reactionRate .GT. 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		endif

	else	
	
	!cell-to-cell diffusion
	
	!Find various parameters needed for reaction rate

		if(proc2==proc1) then
			length2=fineLength
		else
			write(*,*) 'error proc-to-proc diffusion inside fine mesh'
		endif
		area2=length2**2d0
		Vol2=cascadeElementVol
		
		!The area of the shared interface between the volume elements is the minimum of the two volume 
		!elements' face areas (assuming cubic elements, nonuniform)
		if(area1 .GT. area2) then
			areaShared=area2
		else
			areaShared=area1
		endif

		if(proc2==proc1) then
			num2=findNumDefectFine(CascadeCurrent, defectType,cell2)
		else
			write(*,*) 'error proc-to-proc diffusion inside fine mesh'
		endif
		
		reactionRate=Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(length1/2d0 + length2/2d0)
		if(reactionRate .GT. 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		endif

	endif
else
	write(*,*) 'error find reaction rate diffusion'
	findReactionRateDiffFine=0d0
endif

end function


!***************************************************************************************************
! subroutine findReactionInList
!
! points reactionCurrent at the reaction in fine mesh with matching reactants and products in cell
! If reaction is not present, reactionCurrent is not associated and reactionPrev points to the end
! of the list.
!
! Inputs: cell, reactants(:,:), products(:,:), numReactants, numProducts
! Outputs: reactionCurrent, reactionPrev (pointers)
!***************************************************************************************************

!>Subroutine find Reaction In List
!!
!!Points reactionCurrent at the reaction in coarse or fine mesh with matching reactants and products in cell
!!If reaction is not present, reactionCurrent is not associated and reactionPrev points to the end
!!of the list.
!!
!!Inputs: cell, reactants(:,:), products(:,:), numReactants, numProducts
!!Outputs: reactionCurrent, reactionPrev (pointers)

subroutine findReactionInList(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
use mod_srscd_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionCurrent, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
integer numReactants, numProducts, count(numReactants+numProducts), i, j, cell
logical flag

do 10 while(associated(reactionCurrent))
	
	if(reactionCurrent%numReactants==numReactants .AND. reactionCurrent%numProducts==numProducts) then
		
		do 11 i=1,numReactants
			count(i)=0
			do 12 j=1,numSpecies
				if(reactionCurrent%reactants(i,j)==reactants(i,j)) then
					count(i)=count(i)+1
				endif
			12 continue
		11 continue
		
		do 13 i=1,numProducts
			count(i+numReactants)=0
			do 14 j=1,numSpecies
				if(reactionCurrent%products(i,j)==products(i,j)) then
					count(i+numReactants)=count(i+numReactants)+1
				endif
			14 continue
		13 continue
		
		flag=.FALSE.
		do 15 i=1,numReactants+numProducts
			if(count(i) .NE. numSpecies) then
				flag=.TRUE.
			endif
		15 continue
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit
		endif
	endif
	
	reactionPrev=>reactionCurrent
	reactionCurrent=>reactionCurrent%next
10 continue
	
end subroutine


!***************************************************************************************************
!This subroutine will point reactionCurrent at the correct diffusion reaction in the reaction 
!list. 
!***************************************************************************************************

!>Subroutine find Reaction In List (diffusion)
!!
!!Points reactionCurrent at the correct diffusion reaction in coarse or fine mesh with matching reactants and products in cells
!!If reaction is not present, reactionCurrent is not associated and reactionPrev points to the end
!!of the list.
!!
!!Inputs: cells, processors, reactants(:,:)
!!Outputs: reactionCurrent, reactionPrev (pointers)

subroutine findReactionInListDiff(reactionCurrent, reactionPrev, reactants, cell1, cell2, proc1, proc2)
use mod_srscd_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionCurrent, reactionPrev
integer, allocatable :: reactants(:,:)
integer cell1, cell2, proc1, proc2, i, j, count(2)
logical flag

do 10 while(associated(reactionCurrent))
	if(reactionCurrent%numReactants==1 .AND. reactionCurrent%numProducts==1) then
		if(reactionCurrent%cellNumber(1)==cell1 .AND. reactionCurrent%cellNumber(2)==cell2) then
			if(reactionCurrent%taskid(1)==proc1 .AND. reactionCurrent%taskid(2)==proc2) then
				count(1)=0
				do 12 j=1,numSpecies
					if(reactionCurrent%reactants(1,j)==reactants(1,j)) then
						count(1)=count(1)+1
					endif
				12 continue	
		
				count(2)=0
				do 14 j=1,numSpecies
					if(reactionCurrent%products(1,j)==reactants(1,j)) then
						count(2)=count(2)+1
					endif
				14 continue
		
				flag=.FALSE.
				do 15 i=1,2
					if(count(i) .NE. numSpecies) then
						flag=.TRUE.
					endif
				15 continue
				
				if(flag .EQV. .FALSE.) then	!we have found the reaction
					exit
				endif
			endif
		endif
	endif
	reactionPrev=>reactionCurrent
	reactionCurrent=>reactionCurrent%next
10 continue

end subroutine


!***************************************************************************************************
!This subroutine will point reactionCurrent at the correct clustering reaction in the reaction 
!list. If it is not present, this subroutine will point to the end of the list.
!***************************************************************************************************

!>Subroutine find Reaction In List Multiple (clustering)
!!
!!Points reactionCurrent at the clustering reaction in the coarse or fine mesh with matching reactants and products in cell
!!If reaction is not present, reactionCurrent is not associated and reactionPrev points to the end
!!of the list.
!!
!!Inputs: cell, reactants(:,:), products(:,:), numReactants, numProducts
!!Outputs: reactionCurrent, reactionPrev (pointers)

subroutine findReactionInListMultiple(reactionCurrent, reactionPrev, cell, reactants, products, numReactants, numProducts)
use mod_srscd_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionCurrent, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
integer numReactants, numProducts, count(numReactants+numProducts), i, j, cell
logical flag

do 10 while(associated(reactionCurrent))
	
	!Currently set up only for clustering reactions
	if(reactionCurrent%numReactants==2 .AND. reactionCurrent%numProducts==numProducts) then
		
		!check if the reactants are the same as the reactionCurrent reactants
		count(1)=0
		count(2)=0
		do 12 j=1,numSpecies
			if(reactionCurrent%reactants(1,j)==reactants(1,j)) then
				count(1)=count(1)+1
			endif
			if(ReactionCurrent%reactants(2,j)==reactants(2,j)) then
				count(2)=count(2)+1
			endif
		12 continue
		
		!check if the products are the same as the reactionCurrent products
		do 13 i=1,numProducts
			count(i+numReactants)=0
			do 14 j=1,numSpecies
				if(reactionCurrent%products(i,j)==products(i,j)) then
					count(i+numReactants)=count(i+numReactants)+1
				endif
			14 continue
		13 continue
		
		flag=.FALSE.
		do 15 i=1,numReactants+numProducts
			if(count(i) .NE. numSpecies) then
				flag=.TRUE.
			endif
		15 continue
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit
		endif
		
		!*************
		!Here we must repeat the above operation but take into account the possibility of the two 
		!reactants being in the wrong order. Thus we switch the indexes in the reactants comparison
		!step.
		!*************
		!check if the reactants are the same as the reactionCurrent reactants
		count(1)=0
		count(2)=0
		do 22 j=1,numSpecies
			if(reactionCurrent%reactants(2,j)==reactants(1,j)) then
				count(1)=count(1)+1
			endif
			if(reactionCurrent%reactants(1,j)==reactants(2,j)) then
				count(2)=count(2)+1
			endif
		22 continue
		
		!check if the products are the same as the reactionCurrent products
		do 23 i=1,numProducts
			count(i+numReactants)=0
			do 24 j=1,numSpecies
				if(reactionCurrent%products(i,j)==products(i,j)) then
					count(i+numReactants)=count(i+numReactants)+1
				endif
			24 continue
		23 continue
		
		flag=.FALSE.
		do 25 i=1,numReactants+numProducts
			if(count(i) .NE. numSpecies) then
				flag=.TRUE.
			endif
		25 continue
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit
		endif
	endif
	
	reactionPrev=>reactionCurrent
	reactionCurrent=>reactionCurrent%next
10 continue
	
end subroutine

!***************************************************************************************************
! defectCombinationRules (subroutine)
!
! Hard-coded information on what happens to defects when they cluster (for example, annihilation, 
! formation of sessile SIA clusters, etc
!
! Inputs: products(numSpecies), defectTemp
! Outputs: products(numSpeices), updated using the combination rules from defectTemp
!***************************************************************************************************

!>defectCombinationRules (subroutine)
!!
!!Hard-coded information on what happens to defects when they cluster (for example, annihilation, 
!!formation of sessile SIA clusters, etc). 
!!
!!Inputs: products(numSpecies), defectTemp
!!Outputs: products(numSpeices), updated using the combination rules from defectTemp
!!
!!NOTE: defect combination rules are primarily input directly into subroutines addMultiDefectReactions
!!and addMultiDefectReactionsFine, and this subroutine is only called in order to get correct
!!defect combination in the case of cascade-defect interactions. In the future, may want to 
!!move all defect combination rules to this subroutine so that they only need to be changed once.

subroutine defectCombinationRules(products, defectTemp)
use derivedType
use mod_srscd_constants
implicit none

integer products(numSpecies)
integer l
type(defect), pointer :: defectTemp
!SIA+SIA clustering
							
!two 1D clusters coming together to make a sessile cluster
if(products(3) .GT. max3DInt .AND. defectTemp%defectType(3) .GT. max3DInt) then
	products(4)=products(3)
	products(3)=0
endif

do 551 l=1,numSpecies
	products(l)=products(l)+defectTemp%defectType(l)
551 continue

!V+SIA recombination
if(products(2) .GE. products(3)) then
	products(2)=products(2)-products(3)
	products(3)=0
else if(products(3) .GE. products(2)) then
	products(3)=products(3)-products(2)
	products(2)=0
endif

if(products(2) .GE. products(4)) then
	products(2)=products(2)-products(4)
	products(4)=0
else if(products(4) .GE. products(2)) then
	products(4)=products(4)-products(2)
	products(2)=0
endif

!sessile+mobile SIA cluster makes sessile cluster
if(products(3) .NE. 0. .AND. products(4) .NE. 0) then
	products(4)=products(3)+products(4)
	products(3)=0
endif

!sessile cluster becomes mobile again when it shrinks below max3DInt
if(products(4) .NE. 0 .AND. products(4) .LE. max3DInt) then
	products(3)=products(4)
	products(4)=0
endif

end subroutine

!***************************************************************************************************
!
! Subroutine checkReactionLegality
!
! This subroutine looks at the products of a combination reaction and checks to see if the reaction
! is allowed (using hard-coded information). If not, the subroutine returns a value of .FALSE. to 
! isLegal.
!
! Inputs: numProducts, products(numProducts, numSpecies)
! Output: isLegal (boolean variable)
!
!***************************************************************************************************

!>Subroutine check Reaction Legality
!!
!!This subroutine looks at the products of a combination reaction and checks to see if the reaction
!!is allowed (using hard-coded information). If not, the subroutine returns a value of .FALSE. to 
!!isLegal.
!!
!!Inputs: numProducts, products(numProducts, numSpecies)
!!Output: isLegal (boolean variable)

subroutine checkReactionLegality(numProducts, products, isLegal)
use mod_srscd_constants
use DerivedType
implicit none

integer numProducts, i
integer products(numProducts, numSpecies)
logical isLegal

isLegal=.TRUE.

!Check for He+SIA

do 10 i=1,numProducts

	!2014.10.08: Allowing He-SIA clusters to form
	
	if(HeSIAToggle=='no') then
		if(products(i,1) .NE. 0 .AND. products(i,3) .NE. 0) then
			isLegal=.FALSE.
		endif
		
		if(products(i,1) .NE. 0 .AND. products(i,4) .NE. 0) then
			isLegal=.FALSE.
		endif
	else if(HeSIAToggle=='yes') then
		!do nothing
	else
		write(*,*) 'Error unknown HeSIAToggle'
	endif
	
	!Check for interstitial He greater than 4
	
	if(products(i,1) .GT. 4 .AND. products(i,2)==0 .AND. &
		products(i,3)==0 .AND. products(i,4)==0) then
		
		isLegal=.FALSE.
	endif

10 continue

!Check for HeV kick-out
if(numProducts==2) then
	if(products(1,1) .NE. 0 .AND. products(2,3) .NE. 0) then
	
		!we are attempting to see if a HeV-> HeV + SIA kickout is legal
		!This reaction is only allowed if He/V .GE. 5. Therefore, we check
		!the reactant's He/V ratio.
		
		if(dble(products(1,1)/(products(1,2)-1)) .LT. 5d0) then
			isLegal=.FALSE.
		endif
	endif
	
	!Check for He1V1 dissociation: only allow one version of this to avoid
	!double-counting the reaction.
	if(products(1,1)==1 .AND. products(1,2)==0 .AND. products(1,3)==0 .AND. products(1,4)==0) then
	if(products(2,1)==0 .AND. products(2,2)==1 .AND. products(2,3)==0 .AND. products(2,4)==0) then
	
			isLegal=.FALSE.
	
	endif
	endif
endif

end subroutine

!***************************************************************************************************
!
! Function findDPARateLocal(zCoord)
!
! This function finds the DPA rate in a non-uniform implantation profile given a z-coordinate of the
! center of the mesh element. An error is returned if the non-uniform implantation profile does not 
! completely surround the mesh.
!
! Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
! Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file
!
!***************************************************************************************************

!>Function find DPA Rate Local (zCoord)
!!
!!This function finds the DPA rate in a non-uniform implantation profile given a z-coordinate of the
!!center of the mesh element. An error is returned if the non-uniform implantation profile does not 
!!completely surround the mesh. This function uses information input from a file.
!!
!!Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
!!Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file

double precision function findDPARateLocal(zCoord)
use DerivedType
use mod_srscd_constants
implicit none

double precision zCoord, xi
integer i

do 10 i=1,numImplantDataPoints
	if(implantRateData(i,1)==zCoord) then
		findDPARateLocal=implantRateData(i,2)
		exit
	elseif(implantRateData(i,1) .GT. zCoord .AND. i .NE. 1) then
	
		!xi is used in interpolaton between points if zCoord does not fall on a coordinate given in the data file
		xi=(zCoord-implantRateData(i-1,1))/(implantRateData(i,1)-implantRateData(i-1,1))
		
		findDPARateLocal=implantRateData(i,2)*xi+implantRateData(i-1,2)*(1d0-xi)	!interpolate between points in data file
		exit
	elseif(implantRateData(i,1) .GT. zCoord .AND. i == 1) then
		
		!The first point in the data file is greater than the first point in the mesh, meaning that
		!we cannot use the data file. Return an error.
		
		write(*,*) 'error DPA rate file starts after mesh in z-direction'
		exit
	endif
10 continue

if(i==numImplantDataPoints+1) then

	!The last point in the data file is smaller than the last point in the mesh, meaning that we cannot
	!use the data file. Return an error.
	
	write(*,*) 'error DPA rate file ends before mesh in z-direction'

endif

end function

!***************************************************************************************************
!
! Function findHeImplantRateLocal(zCoord)
!
! This function finds the He implant rate rate in a non-uniform implantation profile given a z-coordinate of the
! center of the mesh element. An error is returned if the non-uniform implantation profile does not 
! completely surround the mesh.
!
! Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
! Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file
!
!***************************************************************************************************

!>Function find Helium Implantation Rate Local (zCoord)
!!
!!This function finds the He implant rate rate in a non-uniform implantation profile given a z-coordinate of the
!!center of the mesh element. An error is returned if the non-uniform implantation profile does not 
!!completely surround the mesh.
!!
!!Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
!!Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file

double precision function findHeImplantRateLocal(zCoord)
use DerivedType
use mod_srscd_constants
implicit none

double precision zCoord, xi
integer i

do 10 i=1,numImplantDataPoints
	if(implantRateData(i,1)==zCoord) then
		findHeImplantRateLocal=implantRateData(i,3)
		exit
	elseif(implantRateData(i,1) .GT. zCoord .AND. i .NE. 1) then
	
		!xi is used in interpolaton between points if zCoord does not fall on a coordinate given in the data file
		xi=(zCoord-implantRateData(i-1,1))/(implantRateData(i,1)-implantRateData(i-1,1))
		
		findHeImplantRateLocal=implantRateData(i,3)*xi+implantRateData(i-1,3)*(1d0-xi)	!interpolate between points in data file
		exit
	elseif(implantRateData(i,1) .GT. zCoord .AND. i == 1) then
		
		!The first point in the data file is greater than the first point in the mesh, meaning that
		!we cannot use the data file. Return an error.
		
		write(*,*) 'error He implant rate rate file starts after mesh in z-direction'
		exit
	endif
10 continue

if(i==numImplantDataPoints+1) then

	!The last point in the data file is smaller than the last point in the mesh, meaning that we cannot
	!use the data file. Return an error.
	
	write(*,*) 'error He implant rate rate file ends before mesh in z-direction'

endif

end function

end module

