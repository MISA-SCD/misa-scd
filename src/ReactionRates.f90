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
!***************************************************************************************************

module ReactionRates
implicit none
contains

!***************************************************************************************************
!>Subroutine add single defect reactions - adds reactions to a reaction list that require only
! a single reactant to be carried out. Examples: dissociation, trapping, sinks. Diffusion
! reactions are not included in this subroutine.
!
! Inputs: cell number, defect type (we are adding all single defect reactions associated with
! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!
!Structure of subroutine:
!
!1) Identify if the defect type corresponds to an allowed reaction using reaction lists
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect type and number of defects
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!***************************************************************************************************
subroutine addSingleDefectReactions(cell, defectType)
use mod_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), matNum
integer i, j, k, count, numReactants, numProducts, storeTemp
type(reaction), pointer :: reactionUpdate, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate, totalRateCheck
logical isLegal

nullify(reactionUpdate)
nullify(reactionPrev)
!Dissociation reactions. NOTE: the number of reactants and number of products, as well as the form 
!of the products given the reactants, is hard-coded into this section. This type of hard-coding is
!carried out in each section in this module, and is seen here as unavoidable.

!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
end if

!Dissociation reactions
numReactants=1
numProducts=2
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))
do i=1, numDissocReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. DissocReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. DissocReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= DissocReactions(i,matNum)%min(j) .AND. &
                    ((defectType(j) <= DissocReactions(i,matNum)%max(j)) .OR. DissocReactions(i,matNum)%max(j)==-1)) then
                count=count+1
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do j=1,numSpecies
			reactants(j,1)=defectType(j)
			products(j,2)=DissocReactions(i,matNum)%products(j,1)   !point defects
			products(j,1)=reactants(j,1)-products(j,2)
		end do

		!Dissociation of mobile defects from sessile SIA clusters
		if(products(3,1) < 0) then
		!	if(products(1,1) /= 0 .AND. products(2,1) /= 0) then	!Kick-out mechanism. Cu_V dissociates a SIA_m
		!		products(3,1)=0
		!		products(2,1)=products(2,1)+products(3,2)
		!	else    !dissociation of SIA_m from sessile SIA clusters
				products(3,1)=0
				products(4,1)=products(4,1)-products(3,2)
		!	end if
		end if
		!sessile cluster becomes mobile when it shrinks below max3DInt
		if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
			products(3,1)=products(4,1)
			products(4,1)=0
		end if

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)

	!	call checkReactionLegality(numProducts, products, isLegal)
		
		!find the reaction rate.
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateDissoc(defectType, products, cell, DissocReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then	!in the middle
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else    !at the end of the list
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if

			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate

			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=2
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
				reactionUpdate%products(j,2)=products(j,2)
			end do
			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
		
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		else
			write(*,*) 'error updating reaction list - dissociation'
		end if

		exit

	end if
end do

!Sink reactions
deallocate(reactants)
deallocate(products)
numReactants=1
numProducts=0
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numSinkReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this sink reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. SinkReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. SinkReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= SinkReactions(i,matNum)%min(j)) then
				if((defectType(j) <= SinkReactions(i,matNum)%max(j)) .OR. SinkReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (sinks)
		do j=1,numSpecies
			reactants(j,1)=defectType(j)
		end do

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if CuSIA or large He clusters left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateSink(defectType, cell, SinkReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then    !in the middle
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	if(allocated(reactionUpdate%products)) then
			!		deallocate(reactionUpdate%products)
			!	end if
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else    !at the end of the list
			!	deallocate(reactionUpdate%reactants)
			!	if(allocated(reactionUpdate%products)) then
			!		deallocate(reactionUpdate%products)
			!	end if
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			if(allocated(reactionUpdate%products)) then
				deallocate(reactionUpdate%products)
			endif
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=0
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate

			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
			end do
			do j=1,reactionUpdate%numReactants
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
		
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		end if

		exit

	end if
	
end do

!Impurity reactions
deallocate(reactants)
deallocate(products)
numReactants=1
numProducts=1
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numImpurityReac(matNum)
	count=0
	!Check if the defect type is accepted by this impurity reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. ImpurityReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. ImpurityReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= ImpurityReactions(i,matNum)%min(j)) then
				if((defectType(j) <= ImpurityReactions(i,matNum)%max(j)) .OR. ImpurityReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (impurities)
		!Impurities change defect types from mobile SIA loops to sesile SIA loops. Therefore,
		!we must change the defectType from 0 0 n 0 to 0 0 0 n.
		do j=1,numSpecies
			reactants(j,1)=defectType(j)
			if(reactants(j,1) /= 0) then
				storeTemp=reactants(j,1)    !reactan=0 0 n 0, storeTemp=n
			end if
		end do
		do j=1,numSpecies
			if(ImpurityReactions(i,matNum)%products(j,1)==1) then   !0 0 0 1
				products(j,1)=storeTemp
			else
				products(j,1)=0
			end if
		end do

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if CuSIA or large He clusters left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateImpurity(defectType, cell, ImpurityReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then    !in the middle
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else    !at the end pf the list
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=1
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
			end do
			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		else
			write(*,*) 'error updating reaction list - sinks'
		end if
		exit
	end if
end do

end subroutine

!***************************************************************************************************
!>Subroutine add single defect reactions fine - adds reactions to a reaction list inside a cascade (fine mesh) that require only
! a single reactant to be carried out. Examples: dissociation, trapping, sinks. Diffusion
! reactions are not included in this subroutine.
!
! Inputs: cascade ID number, cell number, defect type (we are adding all single defect reactions associated with
! a single defect type, which was a defect that changed in the previous Monte Carlo step)
! Outputs: updates single-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!
!Structure of subroutine:
!
!1) Identify if the defect type corresponds to an allowed reaction using reaction lists
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect type and number of defects (using fine mesh volume element size)
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!***************************************************************************************************
subroutine addSingleDefectReactionsFine(cascadeID, cell, defectType)
use mod_constants
use DerivedType
implicit none

integer cascadeID,cell, defectType(numSpecies)

type(cascade), pointer :: CascadeCurrent

integer i, j, count, numReactants, numProducts, storeTemp, matNum
type(reaction), pointer :: reactionUpdate, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate, totalRateCheck
logical isLegal

nullify(reactionUpdate)
nullify(reactionPrev)

!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this subroutine
CascadeCurrent=>ActiveCascades
do while(associated(CascadeCurrent))
	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	end if
	CascadeCurrent=>CascadeCurrent%next
end do

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
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numDissocReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. DissocReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. DissocReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= DissocReactions(i,matNum)%min(j)) then
				if((defectType(j) <= DissocReactions(i,matNum)%max(j)) .OR. DissocReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
		
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do j=1, numSpecies
			reactants(j,1)=defectType(j)
			products(j,2)=DissocReactions(i,matNum)%products(j,1)
			products(j,1)=reactants(j,1)-products(j,2)
		end do

        !Dissociation of mobile defects from sessile SIA clusters
		if(products(3,1) < 0) then
		!	if(products(1,1) /= 0 .AND. products(2,1) /= 0) then	!Kick-out mechanism. Cu_V dissociates a SIA_m
		!		products(3,1)=0
		!		products(2,1)=products(2,1)+products(3,2)
		!	else    !dissociation of SIA_m from sessile SIA clusters
				products(3,1)=0
				products(4,1)=products(4,1)-products(3,2)
		!	end if
		end if
		!sessile cluster becomes mobile when it shrinks below max3DInt
		if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
			products(3,1)=products(4,1)
			products(4,1)=0
		end if

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, DissocReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then    !in the middle
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else    !at the end of the list
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=2
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
				reactionUpdate%products(j,2)=products(j,2)
			end do
			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - dissociation'
		end if

        exit

	end if
end do

!Sink reactions
deallocate(reactants)
deallocate(products)
numReactants=1
numProducts=0
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numSinkReac(matNum)
	count=0
	!Check if the defect type is accepted by this sink reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. SinkReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. SinkReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= SinkReactions(i,matNum)%min(j)) then
				if((defectType(j) <= SinkReactions(i,matNum)%max(j)) .OR. SinkReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (sinks)
		do j=1, numSpecies
			reactants(j,1)=defectType(j)
		end do

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateSinkFine(CascadeCurrent, defectType, cell, SinkReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then    !in the middle
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    if(allocated(reactionUpdate%products)) then
            !        deallocate(reactionUpdate%products)
            !    endif
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else    !at the end of list
            !    deallocate(reactionUpdate%reactants)
            !    if(allocated(reactionUpdate%products)) then
            !        deallocate(reactionUpdate%products)
            !    endif
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			if(allocated(reactionUpdate%products)) then
				deallocate(reactionUpdate%products)
			endif
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=0
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate

			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
			end do
			do j=1,reactionUpdate%numReactants
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		end if

        exit

	end if
	
end do

!Impurity reactions
deallocate(reactants)
deallocate(products)
numReactants=1
numProducts=1
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numImpurityReac(matNum)
	count=0
	!Check if the defect type is accepted by this impurity reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. ImpurityReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. ImpurityReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= ImpurityReactions(i,matNum)%min(j)) then
				if((defectType(j) <= ImpurityReactions(i,matNum)%max(j)) .OR. ImpurityReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction

		!Create temporary arrays with the defect types associated with this reaction (impurities)
		!Impurities change defect types from glissile SIA loops to sesile SIA loops. Therefore,
		!we must change the defectType from 0 0 n 0 to 0 0 0 n. This is hard-coded in here.
		do j=1,numSpecies
			reactants(j,1)=defectType(j)
			if(reactants(j,1) /= 0) then
				storeTemp=reactants(j,1)
			end if
		end do
		do j=1,numSpecies
			if(ImpurityReactions(i,matNum)%products(j,1)==1) then
				products(j,1)=storeTemp
			else
				products(j,1)=0
			endif
		end do

        !point reactionUpdate at the reaction and reactionPrev at the reaction before it
        !(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
		reactionUpdate=>CascadeCurrent%reactionList(cell)
		call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA or large He clusters left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateImpurityFine(CascadeCurrent, defectType, cell, ImpurityReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then    !in the middle
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else    !at the end of list
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=1
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
			end do
			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - sinks'
		end if

        exit

	end if
end do

end subroutine

!***************************************************************************************************
!>Subroutine add multi defect reactions - adds reactions to a reaction list that require multiple
!(only using 2 currently) defects to be carried out. This refers mainly to clustering reactions or pinning reactions.
!
! Inputs: cell number, defect types (we are adding all clustering reactions associated with
! a single defect type, which was a defect that changed in the previous Monte Carlo step)
!
!Structure of subroutine:
!
!1) Identify if the defect types both correspond to an allowed clustering reaction reaction using reaction lists
!
!1a) Calculate the resulting product (if the reaction is a clustering reaction), using combination
!rules. Examples include annihilation of vacancy/interstitial pairs and not allowing He-SIA
!clusters to form.
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect types and number of defects
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!***************************************************************************************************
subroutine addMultiDefectReactions(cell, defectType1, defectType2)
use mod_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), matNum
type(reaction), pointer :: reactionUpdate, reactionPrev
integer i, j, count, count2, numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionRate
logical isLegal, isLegalTemp

nullify(reactionUpdate)
nullify(reactionPrev)
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
isLegalTemp =.TRUE.

numReactants=2
numProducts=1
allocate(reactants(numSpecies,numReactants))
do i=1, numClusterReac(matNum)

    !*******************************************************
    !defectType1 = ClusterReactions(i,matNum)%reactants(:,1)
    !defectType2 = ClusterReactions(i,matNum)%reactants(:,2)
    !*******************************************************
	count=0
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	do j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
				count=count+1
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
				if(defectType2(j) >= ClusterReactions(i,matNum)%min(j+numSpecies)) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
						ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
						count=count+1
					end if
				end if
			end if
		else if(defectType1(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
				if(defectType1(j) >= ClusterReactions(i,matNum)%min(j)) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						count=count+1
					end if
				end if
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
				if((defectType2(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
					ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						if(defectType2(j) >= ClusterReactions(i,matNum)%min(j+numSpecies) .AND. &
							defectType1(j) >= ClusterReactions(i,matNum)%min(j)) then
							count=count+1
						end if
					end if
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction

        do j=1,numSpecies
            reactants(j,1)=defectType1(j)
            reactants(j,2)=defectType2(j)
        end do

        !CuV+SIA:
        if(defectType1(1)/=0 .AND. defectType1(2)/=0 .AND.  defectType2(3)>defectType1(2)) then
            numProducts=2
            allocate(products(numSpecies,numProducts))
            !Create temporary arrays with the defect types associated with this reaction (SIA pinning)
            do j=1,numSpecies
                if(j==2) then
                    products(j,1)=0
                else
                    products(j,1)=defectType1(j)
                end if
                if(j==3) then
                    products(j,2)=defectType2(3)-defectType1(2)
                else
                    products(j,2)=defectType2(j)
                end if
            end do

		else if(defectType1(1)/=0 .AND. defectType1(2)/=0 .AND.  defectType2(4)>defectType1(2)) then

			numProducts=2
			allocate(products(numSpecies,numProducts))

			!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
			do j=1,numSpecies

				if(j==2) then
					products(j,1)=0
				else
					products(j,1)=defectType1(j)
				end if
				if(j==4) then
					products(j,2)=defectType2(4)-defectType1(2)
				else
					products(j,2)=defectType2(j)
				end if
			end do
		else
			numProducts=1
			allocate(products(numSpecies,numProducts))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do j=1,numSpecies
				products(j,1)=reactants(j,1)+reactants(j,2)
			end do
		end if
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: Cu, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of CuV+SIA=>CuSIA (should not be allowed here)
		!*******************************************************************************************
        if(numProducts==1) then
            !Vacancy+SIA annihilation - only the larger species remains
            if(products(2,1) >= products(3,1)) then
                products(2,1)=products(2,1)-products(3,1)
                products(3,1)=0
            end if
            if(products(2,1) >= products(4,1)) then
                products(2,1)=products(2,1)-products(4,1)
                products(4,1)=0
            end if

            if(products(2,1) < products(3,1)) then
                products(3,1)=products(3,1)-products(2,1)
                products(2,1)=0
            end if
            if(products(2,1) < products(4,1)) then
                products(4,1)=products(4,1)-products(2,1)
                products(2,1)=0
            end if

            !SIA+SIA clustering
            !two 1D clusters coming together to make a sessile cluster
            if(reactants(3,1) > max3DInt .AND. reactants(3,2) > max3DInt) then
                products(4,1)=products(3,1)
                products(3,1)=0
            end if

            !sessile SIA + mobile SIA = sessile SIA
            if(products(3,1) /= 0 .AND. products(4,1) /= 0) then
                products(4,1)=products(3,1)+products(4,1)
                products(3,1)=0
            end if

            !sessile cluster becomes mobile when it shrinks below max3DInt
            if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
                products(3,1)=products(4,1)
                products(4,1)=0
            end if
        else if(numProducts==2) then
            !sessile cluster becomes mobile when it shrinks below max3DInt
            if(products(4,2) /= 0 .AND. products(4,2) <= max3DInt) then
                products(3,2)=products(4,2)
                products(4,2)=0
            end if
        end if

		!onle point defect can move
		if(pointDefectToggle=='yes') then
			if(products(3,1) /= 0 .AND. products(3,1) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if
			if(numProducts==2) then
				if(products(3,2) /= 0 .AND. products(3,2) > max3DInt) then
					products(4,2)=products(3,2)
					products(3,2)=0
				end if
			end if
		end if

		!Total Annihilation
		count2=0
		do j=1,numSpecies
			if(products(j,1)==0) then
				count2=count2+1
			end if
		end do
		if(count2==numSpecies) then
			!we have completely annihilated the defects
			deallocate(products)
			numProducts=0
			allocate(products(numSpecies,numProducts))
		end if
		
		!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		reactionUpdate=>reactionList(cell)
		call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if CuSIA left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	end if
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then	!in the middle
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else    !at the end of list
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if

			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=2
			reactionUpdate%numProducts=numProducts
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%reactants(j,2)=reactants(j,2)
			end do

            if(numProducts==1) then
                reactionUpdate%products=products
            else if(numProducts==2) then
				do j=1, numSpecies
					reactionUpdate%products(j,1)=products(j,1)
					reactionUpdate%products(j,2)=products(j,2)
				end do
            end if

			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
	
	    deallocate(reactants)
	    deallocate(products)

        exit
		
	end if

    !*******************************************************
    !defectType1 = ClusterReactions(i,matNum)%reactants(:,2)
    !defectType2 = ClusterReactions(i,matNum)%reactants(:,1)
    !*******************************************************
	count=0
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	do j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
				count=count+1
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
				if(defectType2(j) >= ClusterReactions(i,matNum)%min(j)) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		else if(defectType1(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
				if(defectType1(j) >= ClusterReactions(i,matNum)%min(j+numSpecies)) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
						ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
						count=count+1
					end if
				end if
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
				if((defectType1(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
					ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						if(defectType1(j) >= ClusterReactions(i,matNum)%min(j+numSpecies) .AND. &
							defectType2(j) >= ClusterReactions(i,matNum)%min(j)) then
							count=count+1
						end if
					end if
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction

		do j=1,numSpecies
			reactants(j,1)=defectType2(j)
			reactants(j,2)=defectType1(j)
		end do

        !SIA+CuV:
        if(defectType2(1)/=0 .AND. defectType2(2)/=0 .AND.  defectType1(3)>defectType2(2)) then
            numProducts=2
            allocate(products(numSpecies,numProducts))
            !Create temporary arrays with the defect types associated with this reaction (SIA pinning)
            do j=1,numSpecies
                if(j==2) then
                    products(j,1)=0
                else
                    products(j,1)=defectType2(j)
                end if
                if(j==3) then
                    products(j,2)=defectType1(3)-defectType2(2)
                else
                    products(j,2)=defectType1(j)
                end if
            end do

		else if(defectType2(1)/=0 .AND. defectType2(2)/=0 .AND.  defectType1(4)>defectType2(2)) then

			numProducts=2
			allocate(products(numSpecies,numProducts))

			do j=1,numSpecies

				if(j==2) then
					products(j,1)=0
				else
					products(j,1)=defectType2(j)
				end if
				if(j==4) then
					products(j,2)=defectType1(4)-defectType2(2)
				else
					products(j,2)=defectType1(j)
				end if
			end do
		else
			numProducts=1
			allocate(products(numSpecies,numProducts))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			!NOTE: reverse the order of reactants and products so that the reaction is the same
			!(so that findreactioninlistmultiple correctly identifies the reaction)
			do j=1, numSpecies
				products(j,1)=reactants(j,1)+reactants(j,2)
			end do
		end if
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: Cu, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of CuV+SIA=>CuSIA (should not be allowed here)
		!*******************************************************************************************
		if(numProducts==1) then
			!Vacancy+SIA annihilation - only the larger species remains
			if(products(2,1) >= products(3,1)) then
				products(2,1)=products(2,1)-products(3,1)
				products(3,1)=0
			end if
			if(products(2,1) >= products(4,1)) then
				products(2,1)=products(2,1)-products(4,1)
				products(4,1)=0
			end if
			if(products(2,1) < products(3,1)) then
				products(3,1)=products(3,1)-products(2,1)
				products(2,1)=0
			end if
			if(products(2,1) < products(4,1)) then
				products(4,1)=products(4,1)-products(2,1)
				products(2,1)=0
			end if

			!SIA+SIA clustering
			!two 1D clusters coming together to make a sessile cluster
			if(reactants(3,1) > max3DInt .AND. reactants(3,2) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if

			!sessile SIA + mobile SIA = sessile SIA
			if(products(3,1) /= 0 .AND. products(4,1) /= 0) then
				products(4,1)=products(3,1)+products(4,1)
				products(3,1)=0
			end if

			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
				products(3,1)=products(4,1)
				products(4,1)=0
			end if
		else if(numProducts==2) then
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,2) /= 0 .AND. products(4,2) <= max3DInt) then
				products(3,2)=products(4,2)
				products(4,2)=0
			end if
		end if

		!onle point defect can move
		if(pointDefectToggle=='yes') then
			if(products(3,1) /= 0 .AND. products(3,1) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if
			if(numProducts==2) then
				if(products(3,2) /= 0 .AND. products(3,2) > max3DInt) then
					products(4,2)=products(3,2)
					products(3,2)=0
				end if
			end if

		end if

		!Total Annihilation
		count2=0
		do j=1,numSpecies
			if(products(j,1)==0) then
				count2=count2+1
			endif
		end do
		if(count2==numSpecies) then
			!we have completely annihilated the defects
			deallocate(products)
			numProducts=0
			allocate(products(numSpecies,numProducts))
		endif
		
		!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		reactionUpdate=>reactionList(cell)
		call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate. HARD CODED: if HeSIA left over, need to disallow this reaction
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	endif
				
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then	!in the middle
				!deleting reactionUpdate
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else
				!deleting reactionUpdate
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if

			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=2
			reactionUpdate%numProducts=numProducts
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%reactants(j,2)=reactants(j,2)
			end do

            if(numProducts==1) then
                reactionUpdate%products=products
            else if(numProducts==2) then
				do j=1, numSpecies
					reactionUpdate%products(j,1)=products(j,1)
					reactionUpdate%products(j,2)=products(j,2)
				end do
            end if

			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		end if
	
	    deallocate(reactants)
	    deallocate(products)

        exit
		
	end if
	
end do

end subroutine

!***************************************************************************************************
!>Subroutine add multi defect reactions (fine mesh) - adds reactions to a reaction list inside a cascade that require multiple
!(only using 2 currently) defects to be carried out. This refers mainly to clustering reactions
!or pinning reactions.
!
! Inputs: cascade ID number, cell number, defect types (we are adding all clustering reactions associated with
! a single defect type, which was a defect that changed in the previous Monte Carlo step)
! Outputs: updates multi-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!
!Structure of subroutine:
!
!1) Identify if the defect types both correspond to an allowed clustering reaction reaction using reaction lists
!
!1a) Calculate the resulting product (if the reaction is a clustering reaction), using combination
!rules. Examples include annihilation of vacancy/interstitial pairs and not allowing He-SIA
!clusters to form.
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect types and number of defects (and using the cascade volume)
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!***************************************************************************************************
subroutine addMultiDefectReactionsFine(cascadeID, cell, defectType1, defectType2)
use mod_constants
use DerivedType
implicit none

integer cascadeID, cell, defectType1(numSpecies), defectType2(numSpecies), matNum
type(cascade), pointer :: CascadeCurrent
type(reaction), pointer :: reactionUpdate, reactionPrev
integer i, j, count, count2, numReactants, numProducts
integer, allocatable :: reactants(:,:), products(:,:)
double precision reactionrate
logical isLegal

nullify(reactionUpdate)
nullify(reactionPrev)
!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this 
!subroutine

CascadeCurrent=>ActiveCascades
do while(associated(CascadeCurrent))

	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	end if
	CascadeCurrent=>CascadeCurrent%next
end do

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
numReactants=2
numProducts=1
allocate(reactants(numSpecies,numReactants))
do i=1, numClusterReac(matNum)

    !*******************************************************
    !defectType1 = ClusterReactions(i,matNum)%reactants(:,1)
    !defectType2 = ClusterReactions(i,matNum)%reactants(:,2)
    !*******************************************************
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	do j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
				count=count+1
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
				if(defectType2(j) >= ClusterReactions(i,matNum)%min(j+numSpecies)) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
						ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
						count=count+1
					end if
				end if
			end if
		else if(defectType1(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
				if(defectType1(j) >= ClusterReactions(i,matNum)%min(j)) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						count=count+1
					end if
				end if
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
				if((defectType2(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
					ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						if(defectType2(j) >= ClusterReactions(i,matNum)%min(j+numSpecies) .AND. &
							defectType1(j) >= ClusterReactions(i,matNum)%min(j)) then
							count=count+1
						end if
					end if
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction

        !CuV+SIA:
        if(defectType1(1)/=0 .AND.  defectType2(3)>defectType1(2)) then
            numProducts=2
            allocate(products(numSpecies,numProducts))
            !Create temporary arrays with the defect types associated with this reaction (SIA pinning)
            do j=1,numSpecies
				reactants(j,1)=defectType1(j)
				reactants(j,2)=defectType2(j)

                if(j==2) then
                    products(j,1)=0
                else
                    products(j,1)=defectType1(j)
                end if
                if(j==3) then
                    products(j,2)=defectType2(3)-defectType1(2)
                else
                    products(j,2)=defectType2(j)
                end if
            end do

		else if(defectType1(1)/=0 .AND.  defectType2(4)>defectType1(2)) then

			numProducts=2
			allocate(products(numSpecies,numProducts))

			!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
			do j=1,numSpecies
				reactants(j,1)=defectType1(j)
				reactants(j,2)=defectType2(j)

				if(j==2) then
					products(j,1)=0
				else
					products(j,1)=defectType1(j)
				end if
				if(j==4) then
					products(j,2)=defectType2(4)-defectType1(2)
				else
					products(j,2)=defectType2(j)
				end if
			end do

		else
			numProducts=1
			allocate(products(numSpecies,numProducts))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do j=1, numSpecies
				reactants(j,1)=defectType1(j)
				reactants(j,2)=defectType2(j)
				products(j,1)=reactants(j,1)+reactants(j,2)
			end do
		end if
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: Cu, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of CuV+SIA=>CuSIA (should not be allowed here)
		!*******************************************************************************************
		if(numProducts==1) then
			!Vacancy+SIA annihilation - only the larger species remains
			if(products(2,1) >= products(3,1)) then
				products(2,1)=products(2,1)-products(3,1)
				products(3,1)=0
			end if
			if(products(2,1) >= products(4,1)) then
				products(2,1)=products(2,1)-products(4,1)
				products(4,1)=0
			end if
			if(products(2,1) < products(3,1)) then
				products(3,1)=products(3,1)-products(2,1)
				products(2,1)=0
			end if
			if(products(2,1) < products(4,1)) then
				products(4,1)=products(4,1)-products(2,1)
				products(2,1)=0
			end if

			!two 1D SIA clusters coming together to make a sessile cluster
			if(reactants(3,1) > max3DInt .AND. reactants(3,2) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			endif

			!sessile SIA + mobile SIA = sessile SIA
			if(products(3,1) /= 0 .AND. products(4,1) /= 0) then
				products(4,1)=products(3,1)+products(4,1)
				products(3,1)=0
			endif

			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
				products(3,1)=products(4,1)
				products(4,1)=0
			end if
		else if(numProducts==2) then
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,2) /= 0 .AND. products(4,2) <= max3DInt) then
				products(3,2)=products(4,2)
				products(4,2)=0
			end if
		end if

		!onle point defect can move
		if(pointDefectToggle=='yes') then
			if(products(3,1) /= 0 .AND. products(3,1) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if
			if(numProducts==2) then
				if(products(3,2) /= 0 .AND. products(3,2) > max3DInt) then
					products(4,2)=products(3,2)
					products(3,2)=0
				end if
			end if
		end if

		!Total Annihilation
		count2=0
		do j=1,numSpecies
			if(products(j,1)==0) then
				count2=count2+1
			endif
		end do
		if(count2==numSpecies) then
			!we have completely annihilated the defects
			deallocate(products)
			numProducts=0
			allocate(products(numSpecies,numProducts))
		endif
		
		!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		reactionUpdate=>CascadeCurrent%reactionList(cell)
		call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .EQV. .TRUE.) then
			reactionRate=findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, ClusterReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	endif
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then
            !    !deleting reactionUpdate
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=2
			reactionUpdate%numProducts=numProducts
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%reactants(j,2)=reactants(j,2)
				if(numProducts==1) then
					reactionUpdate%products(j,1)=products(j,1)
				else if(numProducts==2) then
					reactionUpdate%products(j,1)=products(j,1)
					reactionUpdate%products(j,2)=products(j,2)
				end if

			end do

			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			!update totalRate
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		else
			write(*,*) 'error updating reaction list - clustering'
		end if
		
	    deallocate(reactants)
	    deallocate(products)
		exit
	end if

    !*******************************************************
    !defectType1 = ClusterReactions(i,matNum)%reactants(:,2)
    !defectType2 = ClusterReactions(i,matNum)%reactants(:,1)
    !*******************************************************
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
	!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
	do j=1,numSpecies
		if(defectType1(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,2)==0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
				count=count+1
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
				if(defectType2(j) >= ClusterReactions(i,matNum)%min(j)) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j)) .OR. &
						ClusterReactions(i,matNum)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		else if(defectType1(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,2) /= 0) then
			if(defectType2(j)==0 .AND. ClusterReactions(i,matNum)%reactants(j,1)==0) then
				if(defectType1(j) >= ClusterReactions(i,matNum)%min(j+numSpecies)) then
					if((defectType1(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
						ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
						count=count+1
					end if
				end if
			else if(defectType2(j) /= 0 .AND. ClusterReactions(i,matNum)%reactants(j,1) /= 0) then
				if((defectType1(j) <= ClusterReactions(i,matNum)%max(j+numSpecies)) .OR. &
					ClusterReactions(i,matNum)%max(j+numSpecies)==-1) then
					if((defectType2(j) <= ClusterReactions(i,matNum)%max(j)) .OR. ClusterReactions(i,matNum)%max(j)==-1) then
						if(defectType1(j) >= ClusterReactions(i,matNum)%min(j+numSpecies) .AND. &
							defectType2(j) >= ClusterReactions(i,matNum)%min(j)) then
							count=count+1
						end if
					end if
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect pair is accepted for this clustering reaction

        !SIA+CuV
        if(defectType2(1)/=0 .AND.  defectType1(3)>defectType2(2)) then
			numProducts=2
            allocate(products(numSpecies,numProducts))
            !Create temporary arrays with the defect types associated with this reaction (SIA pinning)
            do j=1,numSpecies
				reactants(j,2)=defectType1(j)
				reactants(j,1)=defectType2(j)

                if(j==2) then
                    products(j,1)=0
                else
                    products(j,1)=defectType2(j)
                end if
                if(j==3) then
                    products(j,2)=defectType1(3)-defectType2(2)
                else
                    products(j,2)=defectType1(j)
                end if
            end do

		else if(defectType2(1)/=0 .AND.  defectType1(4)>defectType2(2)) then

			numProducts=2
			allocate(products(numSpecies,numProducts))

			do j=1,numSpecies
				reactants(j,2)=defectType1(j)
				reactants(j,1)=defectType2(j)

				if(j==2) then
					products(j,1)=0
				else
					products(j,1)=defectType2(j)
				end if
				if(j==4) then
					products(j,2)=defectType1(4)-defectType2(2)
				else
					products(j,2)=defectType1(j)
				end if
			end do
		else
			numProducts=1
			allocate(products(numSpecies,numProducts))
			
			!Create temporary arrays with the defect types associated with this reaction (clustering)
			do j=1, numSpecies
				reactants(j,2)=defectType1(j)
				reactants(j,1)=defectType2(j)
				products(j,1)=reactants(j,1)+reactants(j,2)
			end do
		endif
		
		!*******************************************************************************************
		!Hard-coded: defect combination rules
		!Here we are assuming 4 species: He, V, SIA_mobile, SIA_sessile
		!We have not yet dealt with the case of HeV+SIA=>HeSIA (should not be allowed here)
		!*******************************************************************************************
		if(numProducts==1) then
			!Vacancy+SIA annihilation - only the larger species remains
			if(products(2,1) >= products(3,1)) then
				products(2,1)=products(2,1)-products(3,1)
				products(3,1)=0
			end if
			if(products(2,1) >= products(4,1)) then
				products(2,1)=products(2,1)-products(4,1)
				products(4,1)=0
			end if
			if(products(2,1) < products(3,1)) then
				products(3,1)=products(3,1)-products(2,1)
				products(2,1)=0
			end if
			if(products(2,1) < products(4,1)) then
				products(4,1)=products(4,1)-products(2,1)
				products(2,1)=0
			end if

			!two 1D SIA clusters coming together to make a sessile cluster
			if(reactants(3,1) > max3DInt .AND. reactants(3,2) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if

			!sessile SIA + mobile SIA = sessile SIA
			if(products(3,1) /= 0 .AND. products(4,1) /= 0) then
				products(4,1)=products(3,1)+products(4,1)
				products(3,1)=0
			end if

			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,1) /= 0 .AND. products(4,1) <= max3DInt) then
				products(3,1)=products(4,1)
				products(4,1)=0
			end if
		else if(numProducts==2) then
			!sessile cluster becomes mobile when it shrinks below max3DInt
			if(products(4,2) /= 0 .AND. products(4,2) <= max3DInt) then
				products(3,2)=products(4,2)
				products(4,2)=0
			end if
		end if

		!onle point defect can move
		if(pointDefectToggle=='yes') then
			if(products(3,1) /= 0 .AND. products(3,1) > max3DInt) then
				products(4,1)=products(3,1)
				products(3,1)=0
			end if
			if(numProducts==2) then
				if(products(3,2) /= 0 .AND. products(3,2) > max3DInt) then
					products(4,2)=products(3,2)
					products(3,2)=0
				end if
			end if
		end if

		!Total Annihilation
		count2=0
		do j=1,numSpecies
			if(products(j,1)==0) then
				count2=count2+1
			endif
		end do
		if(count2==numSpecies) then
			!we have completely annihilated the defects
			deallocate(products)
			numProducts=0
			allocate(products(numSpecies,numProducts))
		endif
		
		!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
		!points to nothing and reactionPrev points to the end of the list.
		!NOTE: if order of reactants is backwards, we might not recognize that we have already added
		!this reaction. Thus we could double-add reactions. Will fix later.
		reactionUpdate=>CascadeCurrent%reactionList(cell)
		call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
		
		!find the reaction rate
	!	call checkReactionLegality(numProducts, products, isLegal)
		
	!	if(isLegal .eqv. .TRUE.) then
			reactionRate=findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, ClusterReactions(i,matNum))
	!	else
	!		reactionRate=0d0
	!	endif
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then
            !    !deleting reactionUpdate
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=2
			reactionUpdate%numProducts=numProducts
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%reactants(j,2)=reactants(j,2)
				if(numProducts==1) then
					reactionUpdate%products(j,1)=products(j,1)
				else if(numProducts==2) then
					reactionUpdate%products(j,1)=products(j,1)
					reactionUpdate%products(j,2)=products(j,2)
				end if
			end do

			do j=1,reactionUpdate%numReactants+reactionUpdate%numProducts
				reactionUpdate%cellNumber(j)=cell
				reactionUpdate%taskid(j)=myProc%taskid
			end do
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 adn reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			!update totalRate
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - clustering'
		endif
	
	    deallocate(reactants)
	    deallocate(products)
        exit
	end if
	
end do

end subroutine

!*************************************************************************************************
!>Subroutine add diffusion reactions - adds reactions to a reaction list representing diffusion between volume elements.
!
! Inputs: initial and final cell numbers, initial and final processors (if diffusion occurs between
!domains of multiple processors), direction (1-6), and defect type (we are adding all diffusion reactions associated with
!a single defect type, which was a defect that changed in the previous Monte Carlo step)
!
!Structure of subroutine:
!
!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!in this element and in the neighboring element
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!
!2015.04.06: Potential bug: if a volume element shares multiple faces with another grain,
!we are only creating one reaction for diffusion from the grain to the grain boundary.
!Need to include diffusion direction in findReactionInListDiff
!***********************************************************************************************
subroutine addDiffusionReactions(cell1, cell2, proc1, proc2, dir, defectType)
use mod_constants
use DerivedType
implicit none

integer cell1, cell2, proc1, proc2, defectType(numSpecies), dir, matNum
integer numReactants, numProducts, i, j, count
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionUpdate, reactionPrev
double precision reactionRate

nullify(reactionUpdate)
nullify(reactionPrev)
!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell1)%material		!For now, assuming that cell1 and cell2 are in the same material.
endif

numReactants=1
numProducts=1
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this diffusion reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. DiffReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= DiffReactions(i,matNum)%min(j)) then
				if((defectType(j) <= DiffReactions(i,matNum)%max(j)) .OR. DiffReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	!this defect type is accepted for this diffusion reaction
		!point reactionUpdate at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do j=1, numSpecies
			reactants(j,1)=defectType(j)
			products(j,1)=defectType(j)
		end do

		!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)
		reactionUpdate=>reactionList(cell1)
		nullify(reactionPrev)
		call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)

		reactionRate=findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, DiffReactions(i,matNum))

		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)-reactionUpdate%reactionRate

			!if(associated(reactionUpdate%next)) then	!in the middle
			!	!deleting reactionUpdate
			!	reactionPrev%next=>reactionUpdate%next
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!else
			!	deallocate(reactionUpdate%reactants)
			!	deallocate(reactionUpdate%products)
			!	deallocate(reactionUpdate%cellNumber)
			!	deallocate(reactionUpdate%taskid)
			!	deallocate(reactionUpdate)
			!	nullify(reactionUpdate)
			!	nullify(reactionPrev%next)
			!end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
			
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=1
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
			end do
			reactionUpdate%cellNumber(1)=cell1
			reactionUpdate%taskid(1)=proc1
			reactionUpdate%cellNumber(2)=cell2
			reactionUpdate%taskid(2)=proc2
			
			reactionUpdate%reactionRate=reactionRate

		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell1)=totalRateVol(cell1)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		end if
        exit
	end if
end do

end subroutine

!***************************************************************************************************
!>Subroutine add diffusion reactions from the coarse mesh to a fine mesh element (cascade) -
!adds reactions to a reaction list representing diffusion between volume elements.
!
! Inputs: coarse mesh cell number, processor, cascade ID number, and defect type (we are adding all diffusion reactions associated with
!a single defect type, which was a defect that changed in the previous Monte Carlo step)
! Outputs: creates new reaction in reactionList for diffusion from coarse mesh to fine mesh
!
!Structure of subroutine:
!
!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!in this element and in the entire fine mesh. A modified reaction distance is used for this reaction rate.
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!***************************************************************************************************

subroutine addDiffusionCoarseToFine(cell, proc, CascadeCurrent, defectType)
use mod_constants
use DerivedType
implicit none

integer cell, proc, defectType(numSpecies), matNum
type(cascade), pointer :: CascadeCurrent
integer numReactants, numProducts, i, j, count, numDefectsFine
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionUpdate, reactionPrev
double precision reactionRate

interface
	integer function findNumDefectTotalFine(defectType, CascadeCurrent)
		use mod_constants
		integer defectType(numSpecies)
		type(cascade), pointer :: CascadeCurrent
	end function
end interface

nullify(reactionUpdate)
nullify(reactionPrev)
!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(cell)%material
endif
	
numReactants=1
numProducts=1
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. DiffReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= DiffReactions(i,matNum)%min(j)) then
				if((defectType(j) <= DiffReactions(i,matNum)%max(j)) .OR. DiffReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		end if
	end do
	
	if(count==numSpecies) then	
		!this defect type is accepted for this dissociation reaction
		
		!point reactionUpdate at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (dissociation)
		do j=1, numSpecies
			reactants(j,1)=defectType(j)
			products(j,1)=defectType(j)
		end do

		!find the reaction in the reaction list. NOTE: -CascadeCurrent%cascadeID is used in place of cell2
		!to identify that this is a reaction from the coarse mesh to the fine mesh into this cascade.
        reactionUpdate=>reactionList(cell)
		nullify(reactionPrev)
		call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell, -CascadeCurrent%cascadeID, proc, proc)
		
		!Find the total number of defects of type defectType in the fine mesh (all cells)
		numDefectsFine=findNumDefectTotalFine(defectType, CascadeCurrent)
		
		!Find the reaction rate for diffusion from coarse to fine mesh
		reactionRate=findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, DiffReactions(i,matNum))
		
		!Here, we update reactionList by either creating a new reaction or updating the current reaction
		
		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then
                !deleting reactionUpdate
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)

		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=1
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
			end do
			reactionUpdate%cellNumber(1)=cell
			reactionUpdate%taskid(1)=proc
			
			!In coarse-to-fine reactions, the cascade ID number is stored as a negative value in the
			!cell number of the diffusion product (negative is signal that coarse-to-fine reaction is occurring)
			reactionUpdate%cellNumber(2)=-CascadeCurrent%cascadeID
			reactionUpdate%taskid(2)=proc
			
			reactionUpdate%reactionRate=reactionRate
			
		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!Update total rate (entire processor and this volume element)
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		end if

        exit

	end if
end do

end subroutine

!***************************************************************************************************
!>Subroutine add diffusion reactions (fine mesh) - adds reactions to a reaction list representing diffusion 
!between volume elements inside a cascade mesh.
!
!Inputs: cascade ID number, initial and final cell numbers, initial and final processors (if diffusion occurs between
!domains of multiple processors), direction (1-6), and defect type (we are adding all diffusion reactions associated with
!a single defect type, which was a defect that changed in the previous Monte Carlo step)
! Outputs: updates multi-defect reaction rates in CascadeCurrent%reactionList(cell) for the correct cascade
!
!Structure of subroutine:
!
!1) Identify if the defect type corresponds to an allowed diffusion reaction reaction using reaction lists
!
!2) Find the reaction in the reaction list for this volume element (if it exists) and go to the
!end of the list if not (NOTE: list is unsorted as of now, 03/31/2015)
!
!3) Calculate the reaction rate based on the defect type and numbers of defects of this type
!in this element and in the neighboring element (using the fine mesh element volume)
!
!4) Update the reaction rate / add the reaction / remove the reaction, depending on if the
!reaction is already present and if the reaction rate is nonzero.
!
!NOTE: includes possibility of diffusion from fine mesh to coarse mesh.
!***************************************************************************************************

subroutine addDiffusionReactionsFine(cascadeID, cell1, cell2, proc1, proc2, dir, defectType)
use mod_constants
use DerivedType
implicit none

integer cascadeID, cell1, cell2, proc1, proc2, dir, defectType(numSpecies), matNum
type(cascade), pointer :: CascadeCurrent
integer numReactants, numProducts, i, j, count
integer, allocatable :: reactants(:,:), products(:,:)
type(reaction), pointer :: reactionUpdate, reactionPrev
double precision reactionRate

!CascadeCurrent pointer should be pointing at the cascade whose ID matches the ID number passed into this subroutine
CascadeCurrent=>ActiveCascades
do while(associated(CascadeCurrent))
	if(CascadeCurrent%cascadeID==cascadeID) then
		exit
	endif
	CascadeCurrent=>CascadeCurrent%next
end do

nullify(reactionUpdate)
nullify(reactionPrev)
!In the case of polycrystal simulations, myMesh(cell)%material is the grain ID, not the material number. Therefore
!we must set all values of matNum=1 in this case (only one material type in polycrystal simulations).
if(numMaterials==1) then
	matNum=1
else
	matNum=myMesh(CascadeCurrent%cellNumber)%material
endif

numReactants=1
numProducts=1
allocate(reactants(numSpecies,numReactants))
allocate(products(numSpecies,numProducts))

do i=1, numDiffReac(matNum)
	count=0
	
	!Check if the defect type is accepted by this dissociation reaction
	do j=1,numSpecies
		if(defectType(j) == 0 .AND. DiffReactions(i,matNum)%reactants(j,1) == 0) then
			count=count+1
		else if(defectType(j) /= 0 .AND. DiffReactions(i,matNum)%reactants(j,1) /= 0) then
			if(defectType(j) >= DiffReactions(i,matNum)%min(j)) then
				if((defectType(j) <= DiffReactions(i,matNum)%max(j)) .OR. DiffReactions(i,matNum)%max(j)==-1) then
					count=count+1
				end if
			end if
		endif
	end do

	if(count==numSpecies) then	!this defect type is accepted for this dissociation reaction
		!point reactionUpdate at the reaction and reactionPrev at the reaction before it
		!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
		!to the end of the list)
		
		!Create temporary arrays with the defect types associated with this reaction (diffusion)
		do j=1,numSpecies
			reactants(j,1)=defectType(j)
			products(j,1)=defectType(j)
		end do

		!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)
		reactionUpdate=>CascadeCurrent%reactionList(cell1)
		nullify(reactionPrev)
		call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)
		
		reactionRate=findReactionRateDiffFine(CascadeCurrent, defectType, cell1, proc1, cell2, proc2, dir, DiffReactions(i,matNum))

		!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
		if(associated(reactionUpdate) .AND. reactionRate==0d0) then
			
			totalRate=totalRate-reactionUpdate%reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)-reactionUpdate%reactionRate

            !if(associated(reactionUpdate%next)) then
            !    reactionPrev%next=>reactionUpdate%next
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !else
            !    deallocate(reactionUpdate%reactants)
            !    deallocate(reactionUpdate%products)
            !    deallocate(reactionUpdate%cellNumber)
            !    deallocate(reactionUpdate%taskid)
            !    deallocate(reactionUpdate)
            !    nullify(reactionUpdate)
            !    nullify(reactionPrev%next)
            !end if
			
			!deleting reactionUpdate
			reactionPrev%next=>reactionUpdate%next
			deallocate(reactionUpdate%reactants)
			deallocate(reactionUpdate%products)
			deallocate(reactionUpdate%cellNumber)
			deallocate(reactionUpdate%taskid)
			nullify(reactionUpdate%next)
			deallocate(reactionUpdate)
			
		
		!if reactionRate .NE. 0 and reaction doesn't exist, then create it. Add to totalRate
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			totalRate=totalRate+reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)+reactionRate
			
			!creating new reaction
			allocate(reactionUpdate)
			reactionUpdate%numReactants=1
			reactionUpdate%numProducts=1
			allocate(reactionUpdate%reactants(numSpecies,reactionUpdate%numReactants))
			allocate(reactionUpdate%products(numSpecies,reactionUpdate%numProducts))
			allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
			allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
			nullify(reactionUpdate%next)
			reactionPrev%next=>reactionUpdate
			do j=1, numSpecies
				reactionUpdate%reactants(j,1)=reactants(j,1)
				reactionUpdate%products(j,1)=products(j,1)
			end do

			reactionUpdate%cellNumber(1)=cell1
			reactionUpdate%taskid(1)=proc1
			reactionUpdate%cellNumber(2)=cell2
			reactionUpdate%taskid(2)=proc2
			
			reactionUpdate%reactionRate=reactionRate

		!if reactionRate==0 and reaction doesn't exist, do nothing
		else if(.NOT. associated(reactionUpdate) .AND. reactionRate==0d0) then
			!do nothing
		!if reactionRate .NE. 0 and reaction exists, update the reaction rate. Add/subtract to totalRate
		else if(associated(reactionUpdate) .AND. reactionRate /= 0d0) then
			
			!update totalRate
			totalRate=totalRate-reactionUpdate%reactionRate+reactionRate
			CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)-reactionUpdate%reactionRate+reactionRate
			
			!update reaction rate
			reactionUpdate%reactionRate=reactionRate
		
		else
			write(*,*) 'error updating reaction list - diffusion'
		endif

        exit

	end if
end do

end subroutine

!***************************************************************************************************
!> Function find reaction rate - finds reaction rate for implantation reaction (He, Frenkel pairs, cascades).
!!
!! Inputs: cell ID, reaction parameters (input from file)
!! Outputs: ReactionRate (the reaction rate  of implanation event).
!!
!! Calculates reaction rates according to hard-coded formulas. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!!
!! NOTE: for the case of non-uniform implantation, reaction rates can be dependent on the z-coordinate
!! of a volume element, using a local DPA rate and helium implantation rate read in from an 
!! input file.
!***************************************************************************************************
double precision function findReactionRate(cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell
type(reactionParameters) :: reactionParameter

double precision Diff, Eb, volume, DPARateLocal, HeImplantRateLocal, zCoord
integer n, numClusters

if(reactionParameter%functionType==10) then	!Frenkel pair implantation

	volume=myMesh(cell)%volume
	
	if(implantDist=='Uniform') then
		findReactionRate=volume*dpaRate/atomSize
	else if(implantDist=='NonUniform') then
		zCoord=myMesh(cell)%coordinates(3)
		DPARateLocal=findDPARateLocal(zCoord)
		findReactionRate=volume*DPARateLocal/atomSize
	else
		write(*,*) 'Error implant distribution not recognized'
	endif
	
else if(reactionParameter%functionType==11) then	!Cascade implantation

	volume=myMesh(cell)%volume
	
	if(implantDist=='Uniform') then
		findReactionRate=volume*dpaRate/(numDisplacedAtoms*atomSize)
	else if(implantDist=='NonUniform') then
		zCoord=myMesh(cell)%coordinates(3)
		DPARateLocal=findDPARateLocal(zCoord)
		findReactionRate=volume*DPARateLocal/(numDisplacedAtoms*atomSize)
	else
		write(*,*) 'Error implant distribution not recognized'
	end if

else if(ReactionParameter%functionType==13) then	!Frenkel-Pair implantation disallowed in grain boundaries
	findReactionRate=0d0
else
	write(*,*) 'error function type', ReactionParameter%functionType
end if

end function

!***************************************************************************************************
!> Function find reaction rate impurity - finds reaction rate for trapping of SIA loops by impurities (Carbon).
!
! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies)
!
! Calculates reaction rates for SIA loop trapping according to hard-coded formulas. Have the possibility to create
! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
! function type associated with it, that function type is assigned in the input file.
!
! In this case, glissile loops become sessile when they interact with impurities. Have the option
! for 3-D or 1-D diffusion of loops.
!
! We have also added the reaction rate for sessile clusters to become mobile again using a
! binding energy in this subroutine. The binding energy is read in from an input file.
!***************************************************************************************************
double precision function findReactionRateImpurity(defectType, cell, reactionParameter)
use mod_constants
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

!if(reactionParameter%functionType==3) then
!	num=findNumDefect(defectType,cell)		!number of clusters of this type
!	Diff=findDiffusivity(matNum,defectType)		!diffusivity of clusters of this type
!	size=findDefectSize(defectType)
	
!	reactionRate=(omegastar+omega*(dble(size)**(1d0/3d0)))*(Diff)*dble(num)*impurityDensity

!else if(reactionParameter%functionType==4) then	!impurityTrapping
if(reactionParameter%functionType==4) then	!impurityTrapping
	Diff=findDiffusivity(matNum,defectType)
	num=findNumDefect(defectType,cell)
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
		Diff*dble(num)*(impurityDensity**2d0)

!else if(reactionParameter%functiontype==5) then
	
	!To Do: create reaction rate for sessile clusters becoming mobile again due to a binding energy.
!	do i=1,numSpecies
!		if(i==3) then
!			productType(i)=defectType(4)
!		else if(i==4) then
!			productType(i)=defectType(3)
!		else
!			productType(i)=defectType(i)
!		end if
!	end do
	
!	num=findNumDefect(defectType,cell)
!	Eb=findBinding(matNum,defectType,productType)
!	Diff=findDiffusivity(matNum,productType)
	
!	reactionRate=omega*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
else
	write(*,*) 'error impurity trapping function type only admits 4 or 5'
	reactionRate=0d0
endif

findReactionRateImpurity=reactionRate

end function

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
!***************************************************************************************************
double precision function findReactionRateImpurityFine(CascadeCurrent, defectType, cell, reactionParameter)
use mod_constants
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
	use mod_constants
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

!if(reactionParameter%functionType==3) then

!	num=findNumDefectFine(CascadeCurrent,defectType,cell)		!number of clusters of this type
!	Diff=findDiffusivity(matNum,defectType)							!diffusivity of clusters of this type
!	size=findDefectSize(defectType)
	
!	reactionRate=(omegastar+omega*(dble(size)**(1d0/3d0)))*&
!			(Diff)*dble(num)*impurityDensity

!else if(reactionParameter%functionType==4) then		!impurityTrapping
if(reactionParameter%functionType==4) then		!impurityTrapping
	Diff=findDiffusivity(matNum,defectType)
	num=findNumDefectFine(CascadeCurrent,defectType,cell)
	size=findDefectSize(defectType)
	
	reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
		Diff*dble(num)*(impurityDensity**2d0)

!else if(reactionParameter%functiontype==5) then
	
	!To Do: create reaction rate for sessile clusters becoming mobile again due to a binding energy.
!	do i=1,numSpecies
!		if(i==3) then
!			productType(i)=defectType(4)
!		else if(i==4) then
!			productType(i)=defectType(3)
!		else
!			productType(i)=defectType(i)
!		end if
!	end do
	
!	num=findNumDefectFine(CascadeCurrent,defectType,cell)
!	Eb=findBinding(matNum,defectType,productType)
!	Diff=findDiffusivity(matNum,productType)
	
!	reactionRate=omega*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	
else
	write(*,*) 'error impurity trapping function type only admits 4 or 5'
	reactionRate=0d0
endif

findReactionRateImpurityFine=reactionRate

end function

!**************************************************************************************************************
!> Function find reaction rate dissociation - finds reaction rate for point defects to dissociate from clusters
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies), dissociating defect type
!!
!! Calculates reaction rates for dissociation of a point defect from a cluster. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!**************************************************************************************************************

double precision function findReactionRateDissoc(defectType, products, cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), products(numSpecies,2), size, num
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

	Diff=findDiffusivity(matNum,products(:,2))	!diffusivity of the defect dissociating from the cluster
	num=findNumDefect(defectType,cell)			!number of clusters
	size=findDefectSize(defectType)				!Hard-coded, rules for determining which species governs the defect size
	Eb=findBinding(matNum,defectType,products(:,2))
	if(defectType(3)>max3DInt .OR. defectType(4)/=0) then	!1D SIA
		reactionRate=omega2D*dble(size)**(1d0/2d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	else	!3D defects
		reactionRate=omega*dble(size)**(1d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	end if
!	reactionRate=omega*dble(size)**(4d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	
else
	write(*,*) 'error dissociation function type only admits 1'
	reactionRate=0d0
end if

findReactionRateDissoc=reactionRate

end function

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
!***************************************************************************************************
double precision function findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), products(numSpecies,2), size, num
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff, Eb
type(cascade), pointer :: CascadeCurrent
type(defect), pointer :: defectTemp
double precision findDiffusivity, findBinding
integer findDefectSize, matNum, grainNum

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_constants
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
	Diff=findDiffusivity(matNum,products(:,2))						!diffusivity of the defect dissociating from the cluster
	defectTemp=>CascadeCurrent%localDefects(cell)
	num=findNumDefectFine(CascadeCurrent, defectType,cell)			!number of clusters
	size=findDefectSize(defectType)									!Hard-coded, rules for determining which species governs the defect size
	Eb=findBinding(matNum,defectType,products(:,2))					!binding energy of single defect to cluster
	if(defectType(3)>max3DInt .OR. defectType(4)/=0) then
		reactionRate=omega2D*dble(size)**(1d0/2d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	else
		reactionRate=omega*dble(size)**(1d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
	end if
!	reactionRate=omega*dble(size)**(4d0/3d0)*Diff*dexp(-Eb/(kboltzmann*temperature))*dble(num)
else
	write(*,*) 'error dissociation function type only admits 1'
	reactionRate=0d0
endif

findReactionRateDissocFine=reactionRate

end function

!***************************************************************************************************
!> Function find reaction rate sink - finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations)
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for capture of a defect by a sink. Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!***************************************************************************************************
double precision function findReactionRateSink(defectType, cell, reactionParameter)
use mod_constants
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
if(reactionParameter%functionType==3) then	!sinkRemoval
	
	num=findNumDefect(defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)		!diffusivity of clusters of this type
	
	if(defectType(3) /= 0) then !SIA_m
		reactionRate=Zint*dislocationDensity*diff*dble(num)
	else	!other defect
		reactionRate=dislocationDensity*diff*dble(num)
	end if
else
	write(*,*) 'error sink function type only admits 3'
	reactionRate=0d0
endif

findReactionRateSink=reactionRate

end function

!***************************************************************************************************
!> Function find reaction rate sink fine - finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations) in the fine mesh (cascade)
!!
!! Inputs: cell ID, reaction parameters (input from file), defect type (array size numSpecies), cascade derived type
!!
!! Calculates reaction rates for capture of a defect by a sink in the fine mesh (in a cascade). Have the possibility to create
!! an arbitrary number of unique formulas for computing reaction rates. Each formula has a
!! function type associated with it, that function type is assigned in the input file.
!***************************************************************************************************
double precision function findReactionRateSinkFine(CascadeCurrent, defectType, cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell, defectType(numSpecies), num, matNum, grainNum
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff
type(cascade), pointer :: CascadeCurrent
double precision findDiffusivity

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_constants
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
if(reactionParameter%functionType==3) then	!sinkRemoval
	
	num=findNumDefectFine(CascadeCurrent, defectType,cell)		!number of clusters of this type
	Diff=findDiffusivity(matNum,defectType)							!diffusivity of clusters of this type
	
	if(defectType(3) /= 0) then !SIA_m
		reactionRate=Zint*dislocationDensity*diff*dble(num)
	else
		reactionRate=dislocationDensity*diff*dble(num)
	end if
else
	write(*,*) 'error sink function type only admits 3'
	reactionRate=0d0
end if

findReactionRateSinkFine=reactionRate

end function

!***************************************************************************************************
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
!***************************************************************************************************
double precision function findReactionRateMultiple(defectType1, defectType2, cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), i, count
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff1, Diff2, vol
integer size1, size2, num1, num2
integer findNumDefect, findDefectSize, matNum, grainNum
double precision findDiffusivity
double precision area

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

count=0
do i=1,numSpecies
	if(defectType1(i)==defectType2(i)) then
		count=count+1
	end if
end do
if(count==numSpecies) then
	!we have two defects of the same type, have to modify the defect numbers for a defect to combine with itself
	num2=num2-1
endif

vol=myMesh(cell)%volume
area=(myMesh(cell)%length)**2d0	!assuming square elements in grain boundary

!if(reactionParameter%functionType==5) then	!Cu+Cu

!	reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*&
!	        (Diff1+Diff2)*dble(num1)*dble(num2)*atomSize/vol
    !reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*&
    !        (Diff1+Diff2)*dble(num1)*dble(num2)*atomSize/vol
        
if(reactionParameter%functionType==6) then	!3D-3D

	if((defectType1(3)>0 .AND. defectType1(3) <= max3DInt) .AND. &
			(defectType2(3)>0 .AND. defectType2(3) <= max3DInt)) then	!3D+3D: 3D(SIA) + 3D(SIA)
		reactionRate=Zint*omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				*atomSize/vol
	else
		reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				*atomSize/vol
	end if
	!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
	!			 *atomSize/vol

else if(reactionParameter%functionType==7) then	!3D-1D: Cu/V/CuV + 1D(SIA)

	!reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
	!		*atomSize/vol
	if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

		!if the first defect is the 1D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefect(defectType2,cell)
		num2=findNumDefect(defectType1,cell)

	else if(defectType2(3) > max3DInt .OR. defectType2(4) > max3DInt) then
		!do nothing, this is the same as the default at the beginning of this subroutine
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	end if

	reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1)*dble(num2)*atomSize/vol+&
			(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0))**4d0*&
					Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)

else if(reactionParameter%functionType==8) then	!3D-1D: 3D(SIA) + 1D(SIA)

	if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then
		
		!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefect(defectType2,cell)
		num2=findNumDefect(defectType1,cell)
		
	else if(defectType2(3) > max3DInt .OR. defectType2(4) > max3DInt) then
		!do nothing, this is the same as the default at the beginning of this subroutine
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	end if

	reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1)*dble(num2)*atomSize/vol+&
		(Zint*(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
		Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)
	!reactionRate=Zint*(omegastar+(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0)))*&
	!		Diff1*num1*num2*atomSize/vol+&
	!		(Zint*(omegastar1D+omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
	!				Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)
		
else if(reactionParameter%functionType==9) then	!1D-1D: 1D(SIA) + 1D(SIA)
	
	reactionRate=(Zint*omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0)))**4d0*&
		(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomSize/vol)**(2d0)

!else if(reactionParameter%functionType==10) then	!Diffusion rate for 2D defect recombination on the grain boundary

!	rad1=(3d0*size1*atomSize/(4d0*pi))**(1d0/3d0)
!	rad2=(3d0*size2*atomSize/(4d0*pi))**(1d0/3d0)
	
!	if(num1==0 .OR. num2==0) then
!		reactionRate=0d0
!	else
!		reactionRate=(4d0*pi*dble(num1)*dble(num2)/area)*&
!						(Diff1/(2d0*dlog(dsqrt(area/(pi*num2))*(1d0/(rad1+rad2)))-1d0)+&
!						 Diff2/(2d0*dlog(dsqrt(area/(pi*num1))*(1d0/(rad1+rad2)))-1d0))
!	end if

!	if(reactionRate < 0d0) then
!		reactionRate=0d0
!	end if
else
	write(*,*) 'error clustering function type only admits 6-9'
	reactionRate=0d0
endif

findReactionRateMultiple=reactionRate
end function

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
!***************************************************************************************************
double precision function findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell, defectType1(numSpecies), defectType2(numSpecies), i, count
type(reactionParameters) :: reactionParameter
double precision reactionRate, Diff1, Diff2, size1, size2, num1, num2, vol
!integer findNumDefectFine
type(cascade), pointer :: CascadeCurrent
integer findDefectSize
double precision findDiffusivity
integer matNum, grainNum

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_constants
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
end if

size1=findDefectSize(defectType1)
size2=findDefectSize(defectType2)
Diff1=findDiffusivity(matNum,defectType1)
Diff2=findDiffusivity(matNum,defectType2)
num1=findNumDefectFine(CascadeCurrent, defectType1,cell)
num2=findNumDefectFine(CascadeCurrent, defectType2,cell)
vol=cascadeElementVol
!Adjust the number of defects in the reaction if both reactants are of the same type

count=0
do i=1,numSpecies
	if(defectType1(i)==defectType2(i)) then
		count=count+1
	end if
end do
if(count==numSpecies) then
	!we have two defects of the same type, have to modify the defect numbers for a defect to combine with itself
	num2=num2-1
end if

!list of clustering reaction functional forms
if(reactionParameter%functionType==6) then	!3D-3D

	if((defectType1(3)>0 .AND. defectType1(3) <= max3DInt) .AND. &
			(defectType2(3)>0 .AND. defectType2(3) <= max3DInt)) then	!3D-3D: 3D(SIA) + 3D(SIA)
		reactionRate=Zint*omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				*atomSize/vol
	else
		reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
				*atomSize/vol
	end if

	!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
	!			 *atomSize/vol

else if(reactionParameter%functionType==7) then	!3D-1D: Cu/V/CuV + 1D(SIA)

	!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
	!			 *atomSize/vol
	if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

		!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefectFine(CascadeCurrent, defectType2,cell)
		num2=findNumDefectFine(CascadeCurrent, defectType1,cell)

	else if(defectType2(3) > max3DInt .OR. defectType2(4) > max3DInt) then
		!do nothing, this is the same as the default at the beginning of this subroutine
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	endif

	reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1*num2)*atomSize/vol+&
			(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0))**4d0*&
					Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)


else if(reactionParameter%functionType==8) then	!3D-1D: 3D(SIA) + 1D(SIA)

	if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then
		
		!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
		size1=findDefectSize(defectType2)
		size2=findDefectSize(defectType1)
		Diff1=findDiffusivity(matNum,defectType2)
		Diff2=findDiffusivity(matNum,defectType1)
		num1=findNumDefectFine(CascadeCurrent, defectType2,cell)
		num2=findNumDefectFine(CascadeCurrent, defectType1,cell)
		
	else if(defectType2(3) > max3DInt .OR. defectType2(4) > max3DInt) then
		!do nothing, this is the same as the default at the beginning of this subroutine
	else
		write(*,*) 'error finding 3D-1D reaction rates for non-SIA loops'
		write(*,*) 'defect type 1', defectType1
		write(*,*) 'defect type 2', defectType2
		reactionRate=0d0
	endif

	reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1*num2)*atomSize/vol+&
		(Zint*(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
		Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)
		
else if(reactionParameter%functionType==9) then	!1D-1D: 1D(SIA) + 1D(SIA)
	
	reactionRate=(Zint*omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0)))**4d0*&
		(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomSize/vol)**(2d0)
else
	write(*,*) 'error clustering function type only admits 6-9'
	reactionRate=0d0
end if

findReactionRateMultipleFine=reactionRate
end function

!***************************************************************************************************
!> Function find reaction rate diffusion - finds reaction rate for defect diffusion between elements
!!
!! Inputs: cell IDs, processor IDs, diffusion direction, reaction parameters (input from file), defect type (array size numSpecies)
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative).
!***************************************************************************************************
double precision function findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
use mod_constants
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
else
	matNum=myMesh(cell1)%material	!not worrying about diffusion between multiple material types right now
end if
grainNum=myMesh(cell1)%material

if(reactionParameter%functionType==2) then	!3D diffusion

	Diff=findDiffusivity(matNum,defectType)
	length1=myMesh(cell1)%length
	area1=length1**2d0
	Vol1=myMesh(cell1)%volume
	num1=findNumDefect(defectType,cell1)
	
	strainE1=findStrainEnergy(defectType, cell1)
	if(strainE1 /= 0d0) then
		write(*,*) 'defectType', defectType
		write(*,*) 'strainE1', strainE1
		read(*,*)
	end if
	
	if(proc2==-1) then	!free surface (NOTE: not taking into account strain-assisted diffusion here)
		areaShared=area1
		reactionRate=Diff*areaShared*(dble(num1)/Vol1)/length1
		if(reactionRate > 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		end if
	else	!cell-to-cell diffusion
		!Find various parameters needed for reaction rate
		if(proc2==proc1) then
			matNeighbor=myMesh(cell2)%material
			length2=myMesh(cell2)%length
			Vol2=myMesh(cell2)%volume
			strainE2=findStrainEnergy(defectType, cell2)
		else	!proc2 /= proc1
			matNeighbor=myBoundary(cell2,dir)%material
			length2=myBoundary(cell2,dir)%length
			Vol2=myBoundary(cell2,dir)%volume
			strainE2=findStrainEnergyBoundary(defectType, cell2, dir)
		end if
		
		area2=length2**2d0
		!The area of the shared interface between the volume elements is the minimum of the two volume 
		!elements' face areas (assuming cubic elements, nonuniform)
		if(area1 > area2) then
			areaShared=area2
		else
			areaShared=area1
		end if
		
		if(grainNum==matNeighbor) then
			if(proc2==proc1) then
				num2=findNumDefect(defectType,cell2)
			else
				num2=findNumDefectBoundary(defectType,cell2,dir)
			end if
			
			alpha=1d0
			
		else
			if(numMaterials > 1) then
				!Sink efficiencies (Default set at 1)
				if(defectType(2) /= 0d0) then
					alpha=alpha_v
				else if(defectType(3) /= 0d0) then
					alpha=alpha_i
				else
					alpha=1d0
				end if
			else
				alpha=1d0
			end if
			num2=0	!If we are diffusing between two material types, assume perfect sink
		end if
		
		if(strainField=='yes') then
			if(strainE1 /= 0d0) then
				write(*,*) 'Strain E1', strainE1, 'Strain E2', strainE2
				write(*,*) 'defect type', defectType
				read(*,*)
			endif
			reactionRate=alpha*Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2+&
				dble(num1)*(strainE1-strainE2)/(Vol1*kboltzmann*temperature))/(length1/2d0 + length2/2d0)
		else
			reactionRate=alpha*Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(length1/2d0 + length2/2d0)
		end if
		
		if(reactionRate > 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		end if
	end if
elseif(reactionParameter%functionType==3) then	!2D diffusion on a plane (rather than 3D diffusion in a volume)
	Diff=findDiffusivity(matNum,defectType)
	length1=myMesh(cell1)%length
	area1=length1**2d0
	num1=findNumDefect(defectType,cell1)
	
	if(proc2==-1) then	!free surface
		lengthShared=length1
		reactionRate=Diff*(lengthShared/length1)*(dble(num1)/area1)
		if(reactionRate > 0d0) then
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
			matNeighbor=myBoundary(cell2,dir)%material
			length2=myBoundary(cell2,dir)%length
		endif
		
		if(grainNum==matNeighbor) then
		
			area2=length2**2d0
			
			!The length of the shared interface between the area elements is the minimum of the two 
			!elements' side lengths (assuming cubic elements, nonuniform)
			if(length1 > length2) then
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
			
			if(reactionRate > 0d0) then
				findReactionRateDiff=reactionRate
			else
				findReactionRateDiff=0d0
			end if
		else
			findReactionRateDiff = 0d0	!Don't let 2D diffusion occur between different material types
		end if
	end if

else if(reactionParameter%functionType==4) then	!Dissociation from grain boundary into bulk volume element
	!(treated as a diffusion reaction because we are moving between two volume elements, but the rate is given by a dissociation rate)
	if(proc2==proc1) then
		matNeighbor=myMesh(cell2)%material
	else
		matNeighbor=myBoundary(cell2,dir)%material
	end if
	
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
		
		if(reactionRate > 0d0) then
			findReactionRateDiff=reactionRate
		else
			findReactionRateDiff=0d0
		end if
	end if
else
	write(*,*) 'error find reaction rate diffusion'
	findReactionRateDiff=0d0
endif

end function

!***************************************************************************************************
!> Function find reaction rate diffusion coarse to fine - finds reaction rate for defect diffusion between a coarse mesh element and a fine (cascade) mesh
!!
!! Inputs: cell ID, processor ID, reaction parameters (input from file), defect type (array size numSpecies), number of defects of this type in the fine mesh already
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative), with a modified diffusion 
!! distance according to the formula in Dunn et al. (Computational Materials Science 2015)
!***************************************************************************************************
double precision function findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, reactionParameter)
use mod_constants
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

if(reactionParameter%functionType==2) then	!3D diffusion

	length1=myMesh(cell)%length			!Length of coarse mesh element

	!Effective length used for diffusion from coarse to fine mesh (assuming fine mesh randomly located
	!within the coarse mesh element, see supporting documents for derivation)
	!NOTE: assuming cubic fine mesh (numxCascade=numyCascade=numzCascade). Can re-derive for non-cubic fine meshes.
	CoarseToFineLength=(length1-numxCascade*fineLength)/(dlog(length1**2d0/(numxCascade*fineLength)**2d0))

	Diff=findDiffusivity(matNum,defectType)
	Vol1=myMesh(cell)%volume			!Volume of coarse mesh element
	num1=findNumDefect(defectType,cell)	!Number of defects in coarse mesh element
	num2=numDefectsFine					!Number of defects in fine mesh (entire mesh)
	
	!Area of surface of fine mesh (entire surface)
	areaShared=2d0*(numxCascade*numyCascade+numxCascade*numzCascade+numyCascade*numzCascade)*fineLength**2d0
	!Volume of fine mesh (entire volume)
	Vol2=numxCascade*numyCascade*numzCascade*cascadeElementVol

	reactionRate=Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(CoarseToFineLength)

	if(reactionRate > 0d0) then
		findReactionRateCoarseToFine=reactionRate
	else
		findReactionRateCoarseToFine=0d0
	end if
else
	write(*,*) 'error find reaction rate diffusion coarse to fine'
	findReactionRateCoarseToFine=0d0
endif

end function

!***************************************************************************************************
!> Function find reaction rate diffusion fine - finds reaction rate for defect diffusion between elements in the fine mesh (inside a cascade)
!!
!! Inputs: cell IDs, processor IDs, diffusion direction, reaction parameters (input from file), defect type (array size numSpecies), cascade derived type
!!
!! Calculates reaction rates for diffusion between volume elements. Using the finite-volume
!! version of Fick's law to find rates (using first derivative). Includes rate for diffusion
!! from fine mesh to coarse mesh.
!***************************************************************************************************
double precision function findReactionRateDiffFine(CascadeCurrent, defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
use mod_constants
use DerivedType
implicit none

integer cell1, proc1, cell2, proc2, defectType(numSpecies), num1, num2, dir, coarseCell
type(ReactionParameters) :: reactionParameter
type(cascade), pointer :: CascadeCurrent
double precision Diff, area1, area2, areaShared, Vol1, Vol2, length,length1, length2, reactionRate, coarseLength
double precision fineToCoarseLength, coarseVolume
integer findNumDefect, matNum, grainNum
double precision findDiffusivity

interface
	integer function findNumDefectFine(CascadeCurrent, defectType, cell)
	use mod_constants
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

if(reactionParameter%functionType==2) then	!3D diffusion

	Diff=findDiffusivity(matNum,defectType)		!function in MaterialInput
	
	length1=fineLength
	area1=length1**2d0
	Vol1=cascadeElementVol
	num1=findNumDefectFine(CascadeCurrent, defectType,cell1)

	if(proc2==-1) then	!Diffuse from fine mesh to free surface
	
		write(*,*) 'error free surface diffusion from inside fine mesh'
		areaShared=area1
		reactionRate=Diff*areaShared*(dble(num1)/Vol1)/length1
		if(reactionRate > 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		end if

	else if(cell2==0) then	!fine-to-coarse diffusion

		!Find information on defects and volume element size in coarse mesh element containing this cascade
		coarseCell=CascadeCurrent%cellNumber	!cell number of coarse mesh element
        length=myMesh(coarseCell)%length		!Length of coarse mesh element
		coarseVolume=myMesh(coarseCell)%volume	!Volume of coarse mesh element
		num2=findNumDefect(defectType, coarseCell)	!Number of defects in coarse mesh element

		!average diffusion distance from cascade element to coarse mesh element (assuming cubic cascade mesh and coarse mesh element)
		!See supporting documentation for derivation of this formula
		!Note: assuming cubic cascade meshes here (used numxCascade for all diffusion directions)
		fineToCoarseLength=(length-numxCascade*fineLength)&
			/(dlog((length-(numxCascade-1)*fineLength)**2d0/(fineLength**2d0)))
			
		!Reaction Rate using the average diffusion distance
		reactionRate=Diff*(fineLength**2d0)*&
			(dble(num1)/Vol1-dble(num2)/coarseVolume)/(fineToCoarseLength)
		
		if(reactionRate > 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		end if
	else	!fine-to-fine diffusion

		!Find various parameters needed for reaction rate
		if(proc2==proc1) then
			length2=fineLength
		else
			write(*,*) 'error proc-to-proc diffusion inside fine mesh'
		end if
		area2=length2**2d0
		Vol2=cascadeElementVol
		
		!The area of the shared interface between the volume elements is the minimum of the two volume 
		!elements' face areas (assuming cubic elements, nonuniform)
		if(area1 > area2) then
			areaShared=area2
		else
			areaShared=area1
		end if

		if(proc2==proc1) then
			num2=findNumDefectFine(CascadeCurrent, defectType,cell2)
		else
			write(*,*) 'error proc-to-proc diffusion inside fine mesh'
		endif
		
		reactionRate=Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(length1/2d0 + length2/2d0)
		if(reactionRate > 0d0) then
			findReactionRateDiffFine=reactionRate
		else
			findReactionRateDiffFine=0d0
		end if
	end if
else
	write(*,*) 'error find reaction rate diffusion'
	findReactionRateDiffFine=0d0
end if

end function


!***************************************************************************************************
! subroutine findReactionInList
!
! points reactionUpdate at the reaction in coarse or fine mesh with matching reactants and products in cell
! If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end
! of the list.
!
! Inputs: cell, reactants(:,:), products(:,:), numReactants, numProducts
! Outputs: reactionUpdate, reactionPrev (pointers)
!***************************************************************************************************
subroutine findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
use mod_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionUpdate, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
integer numReactants, numProducts, count(numReactants+numProducts), i, j, cell
logical flag

outer: do while(associated(reactionUpdate))

	if(reactionUpdate%numReactants==numReactants .AND. reactionUpdate%numProducts==numProducts) then
		
		do i=1,numReactants
			count(i)=0
			do j=1,numSpecies
				if(reactionUpdate%reactants(j,i)==reactants(j,i)) then
					count(i)=count(i)+1
				end if
			end do
		end do
		
		do i=1,numProducts
			count(i+numReactants)=0
			do j=1,numSpecies
				if(reactionUpdate%products(j,i)==products(j,i)) then
					count(i+numReactants)=count(i+numReactants)+1
				end if
			end do
		end do
        flag=.FALSE.
		inter: do i=1,numReactants+numProducts
			if(count(i) /= numSpecies) then
				flag=.TRUE.
				exit inter
			end if
		end do inter
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit outer
		end if
	end if
	
	reactionPrev=>reactionUpdate
	reactionUpdate=>reactionUpdate%next
end do outer
	
end subroutine

!***************************************************************************************************
!>Subroutine find Reaction In List (diffusion)
!!
!!Points reactionUpdate at the correct diffusion reaction in coarse or fine mesh with matching reactants and products in cells
!!If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end of the list.
!!
!!Inputs: cells, processors, reactants(:,:)
!!Outputs: reactionUpdate, reactionPrev (pointers)
!***************************************************************************************************
subroutine findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)
use mod_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionUpdate, reactionPrev
integer, allocatable :: reactants(:,:)
integer cell1, cell2, proc1, proc2, i, j, count(2)
logical flag

do while(associated(reactionUpdate))
	if(reactionUpdate%numReactants==1 .AND. reactionUpdate%numProducts==1) then
		if(reactionUpdate%cellNumber(1)==cell1 .AND. reactionUpdate%cellNumber(2)==cell2) then
			if(reactionUpdate%taskid(1)==proc1 .AND. reactionUpdate%taskid(2)==proc2) then
				count(1)=0
				do j=1,numSpecies
					if(reactionUpdate%reactants(j,1)==reactants(j,1)) then
						count(1)=count(1)+1
					end if
				end do
		
				count(2)=0
				do j=1,numSpecies
					if(reactionUpdate%products(j,1)==reactants(j,1)) then
						count(2)=count(2)+1
					endif
				end do
		
				flag=.FALSE.
				do i=1,2
					if(count(i) /= numSpecies) then
						flag=.TRUE.
						exit
					end if
				end do
				
				if(flag .EQV. .FALSE.) then	!we have found the reaction
					exit
				end if
			end if
		end if
	end if
	reactionPrev=>reactionUpdate
	reactionUpdate=>reactionUpdate%next
end do

end subroutine


!***************************************************************************************************
!This subroutine will point reactionUpdate at the correct clustering reaction in the reaction
!list. If it is not present, this subroutine will point to the end of the list.
!***************************************************************************************************

!>Subroutine find Reaction In List Multiple (clustering)
!!
!!Points reactionUpdate at the clustering reaction in the coarse or fine mesh with matching reactants and products in cell
!!If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end
!!of the list.
!!
!!Inputs: cell, reactants(:,:), products(:,:), numReactants, numProducts
!!Outputs: reactionUpdate, reactionPrev (pointers)

subroutine findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
use mod_constants
use DerivedType
implicit none

type(reaction), pointer :: reactionUpdate, reactionPrev
integer, allocatable :: reactants(:,:), products(:,:)
integer numReactants, numProducts, count(numReactants+numProducts), i, j, cell
logical flag

outer: do while(associated(reactionUpdate))
	
	!Currently set up only for clustering reactions
	if(reactionUpdate%numReactants==2 .AND. reactionUpdate%numProducts==numProducts) then
		
		!check if the reactants are the same as the reactionUpdate reactants
		count(1)=0
		count(2)=0
		do j=1,numSpecies
			if(reactionUpdate%reactants(j,1)==reactants(j,1)) then
				count(1)=count(1)+1
			end if
			if(reactionUpdate%reactants(j,2)==reactants(j,2)) then
				count(2)=count(2)+1
			end if
		end do

		!check if the products are the same as the reactionUpdate products
		do i=1,numProducts
			count(i+numReactants)=0
			do j=1,numSpecies
				if(reactionUpdate%products(j,i)==products(j,i)) then
					count(i+numReactants)=count(i+numReactants)+1
				end if
			end do
		end do
		
		flag=.FALSE.
		inter1: do i=1,numReactants+numProducts
			if(count(i) /= numSpecies) then
				flag=.TRUE.
				exit inter1
			end if
		end do inter1
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit outer
		end if
		
		!*************
		!Here we must repeat the above operation but take into account the possibility of the two 
		!reactants being in the wrong order. Thus we switch the indexes in the reactants comparison
		!step.
		!*************
		!check if the reactants are the same as the reactionUpdate reactants
		count(1)=0
		count(2)=0
		do j=1,numSpecies
			if(reactionUpdate%reactants(j,2)==reactants(j,1)) then
				count(1)=count(1)+1
			end if
			if(reactionUpdate%reactants(j,1)==reactants(j,2)) then
				count(2)=count(2)+1
			end if
		end do

		!check if the products are the same as the reactionUpdate products
		do i=1,numProducts
			count(i+numReactants)=0
			do j=1,numSpecies
				if(reactionUpdate%products(j,i)==products(j,i)) then
					count(i+numReactants)=count(i+numReactants)+1
				end if
			end do
		end do
		
		flag=.FALSE.
		inter2: do i=1,numReactants+numProducts
			if(count(i) /= numSpecies) then
				flag=.TRUE.
				exit inter2
			end if
		end do inter2
		
		if(flag .EQV. .FALSE.) then	!we have found the reaction
			exit outer
		end if
	end if
	
	reactionPrev=>reactionUpdate
	reactionUpdate=>reactionUpdate%next
end do outer
	
end subroutine

!***************************************************************************************************
!>defectCombinationRules (subroutine)
!
!Hard-coded information on what happens to defects when they cluster (for example, annihilation,
!formation of sessile SIA clusters, etc).
!
!Inputs: products(numSpecies), product2(numSpecies), defectTemp, isCombined
!Outputs: products(numSpeices), updated using the combination rules from defectTemp
!
!NOTE: defect combination rules are primarily input directly into subroutines addMultiDefectReactions
!and addMultiDefectReactionsFine, and this subroutine is only called in order to get correct
!defect combination in the case of cascade-defect interactions. In the future, may want to
!move all defect combination rules to this subroutine so that they only need to be changed once.
!***************************************************************************************************
subroutine defectCombinationRules(products, product2, defectTemp, isCombined)
use derivedType
use mod_constants
implicit none

integer products(numSpecies), product2(numSpecies)
integer l
type(defect), pointer :: defectTemp
logical isCombined

isCombined = .TRUE.

!SIA+Cu or Cu+SIA, not combine
if(products(1)/=0 .AND. products(2)==0 .AND. &
        (defectTemp%defectType(3)/=0 .OR. defectTemp%defectType(4)/=0)) then    !Cu+SIA
    isCombined=.FALSE.

else if(defectTemp%defectType(1)/=0 .AND. defectTemp%defectType(2)==0 .AND. &
        (products(3)/=0 .OR. products(4)/=0)) then  !SIA+Cu
    isCombined=.FALSE.

else	!Combine

	!two 1D clusters coming together to make a sessile cluster
	if(products(3) > max3DInt .AND. defectTemp%defectType(3) > max3DInt) then
		products(4)=products(3)
		products(3)=0
	end if

	do l=1,numSpecies
		products(l)=products(l)+defectTemp%defectType(l)
	end do

	!CuV+SIA
	if(products(1)/=0 .AND. products(2)/=0 .AND. products(3) > products(2)) then
    	product2(3)=products(3)-products(2)
    	products(2)=0
    	products(3)=0
	else if(products(1)/=0 .AND. products(2)/=0 .AND. products(4) > products(2)) then
    	product2(4)=products(4)-products(2)
    	products(2)=0
    	products(4)=0
	end if

	!V+SIA recombination
	if(products(2) >= products(3)) then
		products(2)=products(2)-products(3)
		products(3)=0
	else if(products(3) > products(2)) then
		products(3)=products(3)-products(2)
		products(2)=0
	end if

	if(products(2) >= products(4)) then
		products(2)=products(2)-products(4)
		products(4)=0
	else if(products(4) > products(2)) then
		products(4)=products(4)-products(2)
		products(2)=0
	end if

	!sessile+mobile SIA cluster makes sessile cluster
	if(products(3) /= 0. .AND. products(4) /= 0) then
		products(4)=products(3)+products(4)
		products(3)=0
	end if

	!sessile cluster becomes mobile again when it shrinks below max3DInt
	if(products(4) /= 0 .AND. products(4) <= max3DInt) then
		products(3)=products(4)
		products(4)=0
	end if
	if(product2(4)/=0 .AND. product2(4) <= max3DInt) then
		product2(3)=product2(4)
		product2(4)=0
	end if

	if(pointDefectToggle=='yes') then
		if(products(3) /= 0 .AND. products(3) > max3DInt) then
			products(4)=products(3)
			products(3)=0
		end if

		if(product2(3) /= 0 .AND. product2(3) > max3DInt) then
			product2(4)=product2(3)
			product2(3)=0
		end if
	end if

end if

end subroutine

!***************************************************************************************************
! Subroutine checkReactionLegality
!
! This subroutine looks at the products of a combination reaction and checks to see if the reaction
! is allowed (using hard-coded information). If not, the subroutine returns a value of .FALSE. to isLegal.
!
! Inputs: numProducts, products(numProducts, numSpecies)
! Output: isLegal (boolean variable)
!***************************************************************************************************
subroutine checkReactionLegality(numProducts, products, isLegal)
use mod_constants
use DerivedType
implicit none

integer numProducts, i
integer products(numSpecies,numProducts)
logical isLegal

isLegal=.TRUE.

!Check for Cu+SIA
do i=1,numProducts
    if(products(1,i)/=0 .AND. (products(3,i)/=0 .OR. products(4,i)/=0)) then
        isLegal=.FALSE.
    end if
end do

!Check for CuV kick-out
if(numProducts==2) then

	!Check for C1V1 dissociation: only allow one version of this to avoid
	!double-counting the reaction.
	if(products(1,1)==1 .AND. products(2,1)==0 .AND. products(3,1)==0 .AND. products(4,1)==0) then
        if(products(1,2)==0 .AND. products(2,2)==1 .AND. products(3,2)==0 .AND. products(4,2)==0) then
	
			isLegal=.FALSE.
	
	    end if
	end if
end if

end subroutine

!***************************************************************************************************
! Function findDPARateLocal(zCoord)
!
! This function finds the DPA rate in a non-uniform implantation profile given a z-coordinate of the
! center of the mesh element. An error is returned if the non-uniform implantation profile does not 
! completely surround the mesh.
!
! Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
! Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file
!***************************************************************************************************
!>Function find DPA Rate Local (zCoord)
!!
!!This function finds the DPA rate in a non-uniform implantation profile given a z-coordinate of the
!!center of the mesh element. An error is returned if the non-uniform implantation profile does not 
!!completely surround the mesh. This function uses information input from a file.
!!
!!Input: zCoord, the z-coordinate of the center of the mesh element we are looking for
!!Output: the DPA rate at that point in the non-uniform implantation profile, based on the input file
!***************************************************************************************************
double precision function findDPARateLocal(zCoord)
use DerivedType
use mod_constants
implicit none

double precision zCoord, xi
integer i

do i=1,numImplantDataPoints
	if(implantRateData(1,i)==zCoord) then
		findDPARateLocal=implantRateData(2,i)
		exit
	elseif(implantRateData(1,i) > zCoord .AND. i /= 1) then
		!First-order interpolation polynomial: L1(x)=(x-x1)/(x0-x1)*y0+(x-x0)/(x1-x0)*y1
		findDPARateLocal=(zCoord-implantRateData(1,i))/(implantRateData(1,i-1)-implantRateData(1,i))* &
				implantRateData(2,i-1) + (zCoord-implantRateData(1,i-1))/(implantRateData(1,i)-implantRateData(1,i-1))*&
				implantRateData(2,i)

		!xi is used in interpolaton between points if zCoord does not fall on a coordinate given in the data file
!		xi=(zCoord-implantRateData(1,i-1))/(implantRateData(1,i)-implantRateData(1,i-1))
!		findDPARateLocal=implantRateData(2,i)*xi+implantRateData(2,i-1)*(1d0-xi)	!interpolate between points in data file
		exit
	elseif(implantRateData(1,i) > zCoord .AND. i == 1) then
		
		!The first point in the data file is greater than the first point in the mesh, meaning that
		!we cannot use the data file. Return an error.
		
		write(*,*) 'error DPA rate file starts after mesh in z-direction'
		exit
	end if
end do

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
use mod_constants
implicit none

double precision zCoord, xi
integer i

do i=1,numImplantDataPoints
	if(implantRateData(1,i)==zCoord) then
		findHeImplantRateLocal=implantRateData(3,i)
		exit
	elseif(implantRateData(1,i) > zCoord .AND. i /= 1) then
	
		!xi is used in interpolaton between points if zCoord does not fall on a coordinate given in the data file
		xi=(zCoord-implantRateData(1,i-1))/(implantRateData(1,i)-implantRateData(1,i-1))
		
		findHeImplantRateLocal=implantRateData(3,i)*xi+implantRateData(3,i-1)*(1d0-xi)	!interpolate between points in data file
		exit
	elseif(implantRateData(1,i) > zCoord .AND. i == 1) then
		
		!The first point in the data file is greater than the first point in the mesh, meaning that
		!we cannot use the data file. Return an error.
		
		write(*,*) 'error He implant rate rate file starts after mesh in z-direction'
		exit
	endif
end do

if(i==numImplantDataPoints+1) then

	!The last point in the data file is smaller than the last point in the mesh, meaning that we cannot
	!use the data file. Return an error.
	
	write(*,*) 'error He implant rate rate file ends before mesh in z-direction'

end if

end function

end module

