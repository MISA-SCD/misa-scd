module mod_reactionrates
	implicit none
contains

!***************************************************************************************************
!>Subroutine addSingleDefectReactions(cell, defectType)
!adds reactions to a reaction list that require only a single reactant to be carried out.
!Examples: dissociation, trapping, sinks. Diffusion reactions are not included in this subroutine.
!***************************************************************************************************
subroutine addSingleDefectReactions(cell, defectType)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES)
	integer i, j, k, count, numReactants, numProducts, storeTemp
	type(reaction), pointer :: reactionUpdate, reactionPrev
	integer, allocatable :: reactants(:,:), products(:,:)
	double precision reactionRate, totalRateCheck
	logical isLegal

	nullify(reactionUpdate)
	nullify(reactionPrev)

	!Dissociation reactions
	numReactants=1
	numProducts=2
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numDissocReac
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. DissocReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. DissocReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= DissocReactions(i)%min(j) .AND. &
						((defectType(j) <= DissocReactions(i)%max(j)) .OR. &
								DissocReactions(i)%max(j)==-1)) then
					count=count+1
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (dissociation)
			do j=1,SPECIES
				reactants(j,1)=defectType(j)
				products(j,2)=DissocReactions(i)%products(j,1)   !point defects
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

			reactionRate=findReactionRateDissoc(defectType, products, cell, DissocReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
		end if
	end do

	!Sink reactions
	deallocate(reactants)
	deallocate(products)
	numReactants=1
	numProducts=0
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numSinkReac
		count=0

		!Check if the defect type is accepted by this sink reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. SinkReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. SinkReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= SinkReactions(i)%min(j)) then
					if((defectType(j) <= SinkReactions(i)%max(j)) .OR. SinkReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (sinks)
			do j=1,SPECIES
				reactants(j,1)=defectType(j)
			end do

			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
			reactionUpdate=>reactionList(cell)
			call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)

			reactionRate=findReactionRateSink(defectType, cell, SinkReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate

				do j=1, SPECIES
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
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numImpurityReac
		count=0
		!Check if the defect type is accepted by this impurity reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. ImpurityReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. ImpurityReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= ImpurityReactions(i)%min(j)) then
					if((defectType(j) <= ImpurityReactions(i)%max(j)) .OR. &
							ImpurityReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do
		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (impurities)
			!Impurities change defect types from mobile SIA loops to sesile SIA loops. Therefore,
			!we must change the defectType from 0 0 n 0 to 0 0 0 n.
			do j=1,SPECIES
				reactants(j,1)=defectType(j)
				if(reactants(j,1) /= 0) then
					storeTemp=reactants(j,1)    !reactan=0 0 n 0, storeTemp=n
				end if
			end do
			do j=1,SPECIES
				if(ImpurityReactions(i)%products(j,1)==1) then   !0 0 0 1
					products(j,1)=storeTemp
				else
					products(j,1)=0
				end if
			end do

			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
			reactionUpdate=>reactionList(cell)
			call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)

			reactionRate=findReactionRateImpurity(defectType, cell, ImpurityReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
!>Subroutine addSingleDefectReactionsFine(cascadeID, cell, defectType)
!adds reactions to a reaction list inside a cascade (fine mesh) that require only a single reactant to be carried out.
!Examples: dissociation, trapping, sinks. Diffusion reactions are not included in this subroutine.
!***************************************************************************************************
subroutine addSingleDefectReactionsFine(cascadeID, cell, defectType)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cascadeID,cell, defectType(SPECIES)
	type(cascade), pointer :: CascadeCurrent

	integer i, j, count, numReactants, numProducts, storeTemp
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

	!Dissociation reactions.
	numReactants=1
	numProducts=2
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numDissocReac
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. DissocReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. DissocReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= DissocReactions(i)%min(j)) then
					if((defectType(j) <= DissocReactions(i)%max(j)) .OR. &
							DissocReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (dissociation)
			do j=1, SPECIES
				reactants(j,1)=defectType(j)
				products(j,2)=DissocReactions(i)%products(j,1)
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

			reactionRate=findReactionRateDissocFine(CascadeCurrent,defectType,products,cell,DissocReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor)
				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
		end if
	end do

	!Sink reactions
	deallocate(reactants)
	deallocate(products)
	numReactants=1
	numProducts=0
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numSinkReac
		count=0
		!Check if the defect type is accepted by this sink reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. SinkReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. SinkReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= SinkReactions(i)%min(j)) then
					if((defectType(j) <= SinkReactions(i)%max(j)) .OR. SinkReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (sinks)
			do j=1, SPECIES
				reactants(j,1)=defectType(j)
			end do

			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
			reactionUpdate=>CascadeCurrent%reactionList(cell)
			call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)

			reactionRate=findReactionRateSinkFine(CascadeCurrent, defectType, cell, SinkReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor)
				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate

				do j=1, SPECIES
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
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numImpurityReac
		count=0
		!Check if the defect type is accepted by this impurity reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. ImpurityReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. ImpurityReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= ImpurityReactions(i)%min(j)) then
					if((defectType(j) <= ImpurityReactions(i)%max(j)) .OR. &
							ImpurityReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction

			!Create temporary arrays with the defect types associated with this reaction (impurities)
			!Impurities change defect types from glissile SIA loops to sesile SIA loops. Therefore,
			!we must change the defectType from 0 0 n 0 to 0 0 0 n. This is hard-coded in here.
			do j=1,SPECIES
				reactants(j,1)=defectType(j)
				if(reactants(j,1) /= 0) then
					storeTemp=reactants(j,1)
				end if
			end do
			do j=1,SPECIES
				if(ImpurityReactions(i)%products(j,1)==1) then
					products(j,1)=storeTemp
				else
					products(j,1)=0
				endif
			end do

			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points to the end of the list)
			reactionUpdate=>CascadeCurrent%reactionList(cell)
			call findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
			reactionRate=findReactionRateImpurityFine(CascadeCurrent, defectType, cell, ImpurityReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor)
				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
!>Subroutine addMultiDefectReactions(cell, defectType1, defectType2)
!adds reactions to a reaction list that require multiple defects to be carried out.
!This refers mainly to clustering reactions or pinning reactions.
!***************************************************************************************************
subroutine addMultiDefectReactions(cell, defectType1, defectType2)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType1(SPECIES), defectType2(SPECIES)
	type(reaction), pointer :: reactionUpdate, reactionPrev
	integer i, j, count, count2, numReactants, numProducts
	integer, allocatable :: reactants(:,:), products(:,:)
	double precision reactionRate
	logical isLegal, isLegalTemp

	nullify(reactionUpdate)
	nullify(reactionPrev)
	isLegalTemp =.TRUE.

	!Clustering reactions.
	numReactants=2
	!numProducts=1
	allocate(reactants(SPECIES,numReactants))
	do i=1, numClusterReac

		!*******************************************************
		!defectType1 = ClusterReactions(i)%reactants(:,1)
		!defectType2 = ClusterReactions(i)%reactants(:,2)
		!*******************************************************
		count=0
		!Check if the defect type is accepted by this dissociation reaction
		!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
		!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
		do j=1,SPECIES
			if(defectType1(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
					count=count+1
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
					if(defectType2(j) >= ClusterReactions(i)%min(j+SPECIES)) then
						if((defectType2(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
								ClusterReactions(i)%max(j+SPECIES)==-1) then
							count=count+1
						end if
					end if
				end if
			else if(defectType1(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
					if(defectType1(j) >= ClusterReactions(i)%min(j)) then
						if((defectType1(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							count=count+1
						end if
					end if
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
					if((defectType2(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
							ClusterReactions(i)%max(j+SPECIES)==-1) then
						if((defectType1(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							if(defectType2(j) >= ClusterReactions(i)%min(j+SPECIES) .AND. &
									defectType1(j) >= ClusterReactions(i)%min(j)) then
								count=count+1
							end if
						end if
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect pair is accepted for this clustering reaction

			do j=1,SPECIES
				reactants(j,1)=defectType1(j)
				reactants(j,2)=defectType2(j)
			end do

			!CuV+SIA:
			if(defectType1(1)/=0 .AND. defectType1(2)/=0 .AND.  defectType2(3)>defectType1(2)) then
				numProducts=2
				allocate(products(SPECIES,numProducts))
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES

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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do j=1,SPECIES
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
			do j=1,SPECIES
				if(products(j,1)==0) then
					count2=count2+1
				end if
			end do
			if(count2==SPECIES) then
				!we have completely annihilated the defects
				deallocate(products)
				numProducts=0
				allocate(products(SPECIES,numProducts))
			end if

			!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
			!points to nothing and reactionPrev points to the end of the list.
			!NOTE: if order of reactants is backwards, we might not recognize that we have already added
			!this reaction. Thus we could double-add reactions. Will fix later.
			reactionUpdate=>reactionList(cell)
			call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, &
					numReactants, numProducts)

			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

				!deleting reactionUpdate
				reactionPrev%next=>reactionUpdate%next
				deallocate(reactionUpdate%reactants)
				if(allocated(reactionUpdate%products)) then
					deallocate(reactionUpdate%products)
				end if
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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
					reactionUpdate%reactants(j,1)=reactants(j,1)
					reactionUpdate%reactants(j,2)=reactants(j,2)
				end do

				if(numProducts==1) then
					reactionUpdate%products=products
				else if(numProducts==2) then
					do j=1, SPECIES
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
		!defectType1 = ClusterReactions(i)%reactants(:,2)
		!defectType2 = ClusterReactions(i)%reactants(:,1)
		!*******************************************************
		count=0
		!Check if the defect type is accepted by this dissociation reaction
		!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
		!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
		do j=1,SPECIES
			if(defectType1(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
					count=count+1
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
					if(defectType2(j) >= ClusterReactions(i)%min(j)) then
						if((defectType2(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							count=count+1
						end if
					end if
				end if
			else if(defectType1(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
					if(defectType1(j) >= ClusterReactions(i)%min(j+SPECIES)) then
						if((defectType1(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
								ClusterReactions(i)%max(j+SPECIES)==-1) then
							count=count+1
						end if
					end if
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
					if((defectType1(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
							ClusterReactions(i)%max(j+SPECIES)==-1) then
						if((defectType2(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							if(defectType1(j) >= ClusterReactions(i)%min(j+SPECIES) .AND. &
									defectType2(j) >= ClusterReactions(i)%min(j)) then
								count=count+1
							end if
						end if
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect pair is accepted for this clustering reaction

			do j=1,SPECIES
				reactants(j,1)=defectType2(j)
				reactants(j,2)=defectType1(j)
			end do

			!SIA+CuV:
			if(defectType2(1)/=0 .AND. defectType2(2)/=0 .AND.  defectType1(3)>defectType2(2)) then
				numProducts=2
				allocate(products(SPECIES,numProducts))
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				do j=1,SPECIES

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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (clustering)
				!NOTE: reverse the order of reactants and products so that the reaction is the same
				!(so that findreactioninlistmultiple correctly identifies the reaction)
				do j=1, SPECIES
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
			do j=1,SPECIES
				if(products(j,1)==0) then
					count2=count2+1
				endif
			end do
			if(count2==SPECIES) then
				!we have completely annihilated the defects
				deallocate(products)
				numProducts=0
				allocate(products(SPECIES,numProducts))
			endif

			!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
			!points to nothing and reactionPrev points to the end of the list.
			!NOTE: if order of reactants is backwards, we might not recognize that we have already added
			!this reaction. Thus we could double-add reactions. Will fix later.
			reactionUpdate=>reactionList(cell)
			call findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, &
					numReactants, numProducts)

			reactionRate=findReactionRateMultiple(defectType1, defectType2, cell, ClusterReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

				!deleting reactionUpdate
				reactionPrev%next=>reactionUpdate%next
				deallocate(reactionUpdate%reactants)
				if(allocated(reactionUpdate%products)) then
					deallocate(reactionUpdate%products)
				end if
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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
					reactionUpdate%reactants(j,1)=reactants(j,1)
					reactionUpdate%reactants(j,2)=reactants(j,2)
				end do

				if(numProducts==1) then
					reactionUpdate%products=products
				else if(numProducts==2) then
					do j=1, SPECIES
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
!>Subroutine addMultiDefectReactionsFine(cascadeID, cell, defectType1, defectType2)
!adds reactions to a reaction list inside a cascade that require multiple defects to be carried out.
!This refers mainly to clustering reactions or pinning reactions.
!***************************************************************************************************
subroutine addMultiDefectReactionsFine(cascadeID, cell, defectType1, defectType2)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cascadeID, cell, defectType1(SPECIES), defectType2(SPECIES)
	type(cascade), pointer :: CascadeCurrent
	type(reaction), pointer :: reactionUpdate, reactionPrev
	integer i, j, count, count2, numReactants, numProducts
	integer, allocatable :: reactants(:,:), products(:,:)
	double precision reactionrate
	logical isLegal

	nullify(reactionUpdate)
	nullify(reactionPrev)

	CascadeCurrent=>ActiveCascades
	do while(associated(CascadeCurrent))

		if(CascadeCurrent%cascadeID==cascadeID) then
			exit
		end if
		CascadeCurrent=>CascadeCurrent%next
	end do

	!Clustering reactions.
	numReactants=2
	numProducts=1
	allocate(reactants(SPECIES,numReactants))
	do i=1, numClusterReac

		!*******************************************************
		!defectType1 = ClusterReactions(i)%reactants(:,1)
		!defectType2 = ClusterReactions(i)%reactants(:,2)
		!*******************************************************
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
		!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
		do j=1,SPECIES
			if(defectType1(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
					count=count+1
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
					if(defectType2(j) >= ClusterReactions(i)%min(j+SPECIES)) then
						if((defectType2(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
								ClusterReactions(i)%max(j+SPECIES)==-1) then
							count=count+1
						end if
					end if
				end if
			else if(defectType1(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
					if(defectType1(j) >= ClusterReactions(i)%min(j)) then
						if((defectType1(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							count=count+1
						end if
					end if
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
					if((defectType2(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
							ClusterReactions(i)%max(j+SPECIES)==-1) then
						if((defectType1(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							if(defectType2(j) >= ClusterReactions(i)%min(j+SPECIES) .AND. &
									defectType1(j) >= ClusterReactions(i)%min(j)) then
								count=count+1
							end if
						end if
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect pair is accepted for this clustering reaction

			!CuV+SIA:
			if(defectType1(1)/=0 .AND.  defectType2(3)>defectType1(2)) then
				numProducts=2
				allocate(products(SPECIES,numProducts))
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do j=1, SPECIES
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
			do j=1,SPECIES
				if(products(j,1)==0) then
					count2=count2+1
				endif
			end do
			if(count2==SPECIES) then
				!we have completely annihilated the defects
				deallocate(products)
				numProducts=0
				allocate(products(SPECIES,numProducts))
			endif

			!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
			!points to nothing and reactionPrev points to the end of the list.
			!NOTE: if order of reactants is backwards, we might not recognize that we have already added
			!this reaction. Thus we could double-add reactions. Will fix later.
			reactionUpdate=>CascadeCurrent%reactionList(cell)
			call findReactionInListMultiple(reactionUpdate,reactionPrev,cell,reactants,products,numReactants,numProducts)

			reactionRate=findReactionRateMultipleFine(CascadeCurrent,defectType1,defectType2,cell,ClusterReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then
				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
		!defectType1 = ClusterReactions(i)%reactants(:,2)
		!defectType2 = ClusterReactions(i)%reactants(:,1)
		!*******************************************************
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		!NOTE: we must check if defectType1 matches with ClusterReactions%reactants(1) and reactants(2)
		!and vice versa with defectType2. We only want to make one reaction rate per pair of reactants.
		do j=1,SPECIES
			if(defectType1(j)==0 .AND. ClusterReactions(i)%reactants(j,2)==0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
					count=count+1
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
					if(defectType2(j) >= ClusterReactions(i)%min(j)) then
						if((defectType2(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							count=count+1
						end if
					end if
				end if
			else if(defectType1(j) /= 0 .AND. ClusterReactions(i)%reactants(j,2) /= 0) then
				if(defectType2(j)==0 .AND. ClusterReactions(i)%reactants(j,1)==0) then
					if(defectType1(j) >= ClusterReactions(i)%min(j+SPECIES)) then
						if((defectType1(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
								ClusterReactions(i)%max(j+SPECIES)==-1) then
							count=count+1
						end if
					end if
				else if(defectType2(j) /= 0 .AND. ClusterReactions(i)%reactants(j,1) /= 0) then
					if((defectType1(j) <= ClusterReactions(i)%max(j+SPECIES)) .OR. &
							ClusterReactions(i)%max(j+SPECIES)==-1) then
						if((defectType2(j) <= ClusterReactions(i)%max(j)) .OR. &
								ClusterReactions(i)%max(j)==-1) then
							if(defectType1(j) >= ClusterReactions(i)%min(j+SPECIES) .AND. &
									defectType2(j) >= ClusterReactions(i)%min(j)) then
								count=count+1
							end if
						end if
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect pair is accepted for this clustering reaction

			!SIA+CuV
			if(defectType2(1)/=0 .AND.  defectType1(3)>defectType2(2)) then
				numProducts=2
				allocate(products(SPECIES,numProducts))
				!Create temporary arrays with the defect types associated with this reaction (SIA pinning)
				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				do j=1,SPECIES
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
				allocate(products(SPECIES,numProducts))

				!Create temporary arrays with the defect types associated with this reaction (clustering)
				do j=1, SPECIES
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
			do j=1,SPECIES
				if(products(j,1)==0) then
					count2=count2+1
				endif
			end do
			if(count2==SPECIES) then
				!we have completely annihilated the defects
				deallocate(products)
				numProducts=0
				allocate(products(SPECIES,numProducts))
			endif

			!findReactionInList points reactionUpdate at the reaction if it already exists. If not, reactionUpdate
			!points to nothing and reactionPrev points to the end of the list.
			!NOTE: if order of reactants is backwards, we might not recognize that we have already added
			!this reaction. Thus we could double-add reactions. Will fix later.
			reactionUpdate=>CascadeCurrent%reactionList(cell)
			call findReactionInListMultiple(reactionUpdate,reactionPrev,cell,reactants,products,numReactants,numProducts)

			reactionRate=findReactionRateMultipleFine(CascadeCurrent,defectType1,defectType2,cell,ClusterReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then
				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell)=CascadeCurrent%totalRate(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
!>Subroutine addDiffusionReactions(cell1, cell2, proc1, proc2, dir, defectType)
!adds reactions to a reaction list representing diffusion between volume elements.
!***********************************************************************************************
subroutine addDiffusionReactions(cell1, cell2, proc1, proc2, dir, defectType)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell1, cell2, proc1, proc2, defectType(SPECIES), dir
	integer numReactants, numProducts, i, j, count
	integer, allocatable :: reactants(:,:), products(:,:)
	type(reaction), pointer :: reactionUpdate, reactionPrev
	double precision reactionRate

	nullify(reactionUpdate)
	nullify(reactionPrev)

	numReactants=1
	numProducts=1
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numDiffReac
		count=0

		!Check if the defect type is accepted by this diffusion reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. DiffReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. DiffReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= DiffReactions(i)%min(j)) then
					if((defectType(j) <= DiffReactions(i)%max(j)) .OR. DiffReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then	!this defect type is accepted for this diffusion reaction
			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
			!to the end of the list)

			!Create temporary arrays with the defect types associated with this reaction (dissociation)
			do j=1, SPECIES
				reactants(j,1)=defectType(j)
				products(j,1)=defectType(j)
			end do

			!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)
			reactionUpdate=>reactionList(cell1)
			nullify(reactionPrev)
			call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)

			reactionRate=findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, DiffReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell1)=totalRateVol(cell1)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
!>Subroutine addDiffusionCoarseToFine(cell, proc, CascadeCurrent, defectType)
!adds reactions to a reaction list representing diffusion between volume elements.
!***************************************************************************************************
subroutine addDiffusionCoarseToFine(cell, proc, CascadeCurrent, defectType)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, proc, defectType(SPECIES)
	type(cascade), pointer :: CascadeCurrent
	integer numReactants, numProducts, i, j, count, numDefectsFine
	integer, allocatable :: reactants(:,:), products(:,:)
	type(reaction), pointer :: reactionUpdate, reactionPrev
	double precision reactionRate

	interface
		integer function findNumDefectTotalFine(defectType, CascadeCurrent)
			use mod_globalVariables
			integer defectType(SPECIES)
			type(cascade), pointer :: CascadeCurrent
		end function
	end interface

	nullify(reactionUpdate)
	nullify(reactionPrev)

	numReactants=1
	numProducts=1
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numDiffReac
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. DiffReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. DiffReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= DiffReactions(i)%min(j)) then
					if((defectType(j) <= DiffReactions(i)%max(j)) .OR. DiffReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			end if
		end do

		if(count==SPECIES) then
			!this defect type is accepted for this dissociation reaction

			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
			!to the end of the list)

			!Create temporary arrays with the defect types associated with this reaction (dissociation)
			do j=1, SPECIES
				reactants(j,1)=defectType(j)
				products(j,1)=defectType(j)
			end do

			!find the reaction in the reaction list. NOTE: -CascadeCurrent%cascadeID is used in place of cell2
			!to identify that this is a reaction from the coarse mesh to the fine mesh into this cascade.
			reactionUpdate=>reactionList(cell)
			nullify(reactionPrev)
			call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell, -CascadeCurrent%cascadeID, &
					proc, proc)

			!Find the total number of defects of type defectType in the fine mesh (all cells)
			numDefectsFine=findNumDefectTotalFine(defectType, CascadeCurrent)

			!Find the reaction rate for diffusion from coarse to fine mesh
			reactionRate=findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, DiffReactions(i))

			!Here, we update reactionList by either creating a new reaction or updating the current reaction

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				!Update total rate (entire processor and this volume element)
				totalRate=totalRate-reactionUpdate%reactionRate
				totalRateVol(cell)=totalRateVol(cell)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
				reactionUpdate%reactionRate=reactionRate
			else
				write(*,*) 'error updating reaction list - diffusion'
			end if
			exit
		end if
	end do

end subroutine

!***************************************************************************************************
!>Subroutine addDiffusionReactionsFine(cascadeID, cell1, cell2, proc1, proc2, dir, defectType)
!adds reactions to a reaction list representing diffusion between volume elements inside a cascade mesh.
!***************************************************************************************************
subroutine addDiffusionReactionsFine(cascadeID, cell1, cell2, proc1, proc2, dir, defectType)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cascadeID, cell1, cell2, proc1, proc2, dir, defectType(SPECIES)
	type(cascade), pointer :: CascadeCurrent
	integer numReactants, numProducts, i, j, count
	integer, allocatable :: reactants(:,:), products(:,:)
	type(reaction), pointer :: reactionUpdate, reactionPrev
	double precision reactionRate

	CascadeCurrent=>ActiveCascades
	do while(associated(CascadeCurrent))
		if(CascadeCurrent%cascadeID==cascadeID) then
			exit
		endif
		CascadeCurrent=>CascadeCurrent%next
	end do

	nullify(reactionUpdate)
	nullify(reactionPrev)

	numReactants=1
	numProducts=1
	allocate(reactants(SPECIES,numReactants))
	allocate(products(SPECIES,numProducts))
	do i=1, numDiffReac
		count=0

		!Check if the defect type is accepted by this dissociation reaction
		do j=1,SPECIES
			if(defectType(j) == 0 .AND. DiffReactions(i)%reactants(j,1) == 0) then
				count=count+1
			else if(defectType(j) /= 0 .AND. DiffReactions(i)%reactants(j,1) /= 0) then
				if(defectType(j) >= DiffReactions(i)%min(j)) then
					if((defectType(j) <= DiffReactions(i)%max(j)) .OR. DiffReactions(i)%max(j)==-1) then
						count=count+1
					end if
				end if
			endif
		end do

		if(count==SPECIES) then	!this defect type is accepted for this dissociation reaction
			!point reactionUpdate at the reaction and reactionPrev at the reaction before it
			!(if reaction does not already exist, reactionUpdate is unallocated and reactionPrev points
			!to the end of the list)

			!Create temporary arrays with the defect types associated with this reaction (diffusion)
			do j=1,SPECIES
				reactants(j,1)=defectType(j)
				products(j,1)=defectType(j)
			end do

			!find the reaction in the reaction list (INCLUDING DIFFUSION DIRECTIONS)
			reactionUpdate=>CascadeCurrent%reactionList(cell1)
			nullify(reactionPrev)
			call findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)

			reactionRate=findReactionRateDiffFine(CascadeCurrent, defectType, cell1, proc1, cell2, proc2, dir, &
					DiffReactions(i))

			!if reactionRate==0 and reaction already exists, then delete it. Subtract from totalRate.
			if(associated(reactionUpdate) .AND. reactionRate==0d0) then

				totalRate=totalRate-reactionUpdate%reactionRate
				CascadeCurrent%totalRate(cell1)=CascadeCurrent%totalRate(cell1)-reactionUpdate%reactionRate

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
				allocate(reactionUpdate%reactants(SPECIES,reactionUpdate%numReactants))
				allocate(reactionUpdate%products(SPECIES,reactionUpdate%numProducts))
				allocate(reactionUpdate%cellNumber(reactionUpdate%numReactants+reactionUpdate%numProducts))
				allocate(reactionUpdate%taskid(reactionUpdate%numReactants+reactionUpdate%numProducts))
				nullify(reactionUpdate%next)
				reactionPrev%next=>reactionUpdate
				do j=1, SPECIES
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
			end if
			exit
		end if
	end do

end subroutine

!***************************************************************************************************
!> Function findReactionRate(cell, reactionParameter)
!finds reaction rate for implantation reaction (Frenkel pairs, cascades).
!***************************************************************************************************
double precision function findReactionRate(cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell
	type(reactionParameters) :: reactionParameter
	double precision Diff, Eb, volume, DPARateLocal, HeImplantRateLocal, zCoord
	integer n, numClusters

	if(reactionParameter%functionType==1) then	!Frenkel pair implantation

		volume=myMesh(cell)%volume

		if(implantDist=='uniform') then
			findReactionRate=volume*dpaRate/atomSize
		else
			write(*,*) 'Error implant distribution not recognized'
		endif

	else if(reactionParameter%functionType==2) then	!Cascade implantation

		volume=myMesh(cell)%volume

		if(implantDist=='uniform') then
			findReactionRate=volume*dpaRate/(numDisplacedAtoms*atomSize)
		else
			write(*,*) 'Error implant distribution not recognized'
		end if

	else if(ReactionParameter%functionType==6) then	!Frenkel-Pair implantation disallowed in grain boundaries
		findReactionRate=0d0
	else
		write(*,*) 'error function type', ReactionParameter%functionType
	end if

end function

!***************************************************************************************************
!> Function findReactionRateImpurity(defectType, cell, reactionParameter)
!finds reaction rate for trapping of SIA loops by impurities (Carbon).
!***************************************************************************************************
double precision function findReactionRateImpurity(defectType, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), num, size
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff
	double precision findDiffusivity, findBinding, Eb
	integer findDefectSize, findNumDefect, grainNum

	if(polycrystal=='yes') then
		grainNum=myMesh(cell)%material
	else
		grainNum=myMesh(cell)%material
	endif

	if(reactionParameter%functionType==14) then	!impurityTrapping
		Diff=findDiffusivity(defectType)
		num=findNumDefect(defectType,cell)
		size=findDefectSize(defectType)

		reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
				Diff*dble(num)*(impurityDensity**2d0)
	else
		write(*,*) 'error impurity trapping function type only admits 4 or 5'
		reactionRate=0d0
	end if

	findReactionRateImpurity=reactionRate

end function

!***************************************************************************************************
!> Function findReactionRateImpurityFine(CascadeCurrent, defectType, cell, reactionParameter)
!finds reaction rate for trapping of SIA loops by impurities (Carbon) inside a fine mesh (Cascade).
!***************************************************************************************************
double precision function findReactionRateImpurityFine(CascadeCurrent, defectType, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), num, size
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff
	type(cascade), pointer :: cascadeCurrent
	double precision findDiffusivity, findBinding, Eb
	integer findDefectSize, findNumDefect, grainNum

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_globalVariables
			type(cascade), pointer :: CascadeCurrent
			integer defectType(SPECIES), cell
		end function
	end interface

	if(polycrystal=='yes') then
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	else
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	end if

	if(reactionParameter%functionType==14) then		!impurityTrapping
		Diff=findDiffusivity(defectType)
		num=findNumDefectFine(CascadeCurrent,defectType,cell)
		size=findDefectSize(defectType)

		reactionRate=(omegastar1D+omegacircle1D*dble(size)**(1d0/2d0)+omega1D)**4d0*&
				Diff*dble(num)*(impurityDensity**2d0)
	else
		write(*,*) 'error impurity trapping function type only admits 4 or 5'
		reactionRate=0d0
	endif

	findReactionRateImpurityFine=reactionRate

end function

!**************************************************************************************************************
!> Function findReactionRateDissoc(defectType, products, cell, reactionParameter)
!finds reaction rate for point defects to dissociate from clusters
!**************************************************************************************************************
double precision function findReactionRateDissoc(defectType, products, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), products(SPECIES,2), size, num
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff, Eb
	double precision findDiffusivity, findBinding
	integer findNumDefect, findDefectSize, grainNum

	if(polycrystal=='yes') then
		grainNum=myMesh(cell)%material
	else
		grainNum=myMesh(cell)%material
	end if

	!dissociation
	if(reactionParameter%functionType==11) then

		Diff=findDiffusivity(products(:,2))	!diffusivity of the defect dissociating from the cluster
		num=findNumDefect(defectType,cell)			!number of clusters
		size=findDefectSize(defectType)				!Hard-coded, rules for determining which species governs the defect size
		Eb=findBinding(defectType,products(:,2))
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
!> Function findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, reactionParameter)
!finds reaction rate for point defects to dissociate from clusters in the fine mesh (cascade)
!***************************************************************************************************
double precision function findReactionRateDissocFine(CascadeCurrent, defectType, products, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), products(SPECIES,2), size, num
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff, Eb
	type(cascade), pointer :: CascadeCurrent
	type(defect), pointer :: defectTemp
	double precision findDiffusivity, findBinding
	integer findDefectSize, grainNum

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_constants
			use mod_globalVariables
			integer cell, defectType(SPECIES)
			type(cascade), pointer :: CascadeCurrent
		end function
	end interface

	if(polycrystal=='yes') then
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	else
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	endif

	!dissociation
	if(reactionParameter%functionType==11) then

		!dissocation reactions
		Diff=findDiffusivity(products(:,2))						!diffusivity of the defect dissociating from the cluster
		defectTemp=>CascadeCurrent%localDefects(cell)
		num=findNumDefectFine(CascadeCurrent, defectType,cell)			!number of clusters
		size=findDefectSize(defectType)									!Hard-coded, rules for determining which species governs the defect size
		Eb=findBinding(defectType,products(:,2))					!binding energy of single defect to cluster
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
!> Function findReactionRateSink(defectType, cell, reactionParameter)
!finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations)
!***************************************************************************************************
double precision function findReactionRateSink(defectType, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), num
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff
	double precision findDiffusivity
	integer findNumDefect, grainNum

	if(polycrystal=='yes') then
		grainNum=myMesh(cell)%material
	else
		grainNum=myMesh(cell)%material
	end if

	!sink reaction function type=12
	if(reactionParameter%functionType==12) then	!sinkRemoval

		num=findNumDefect(defectType,cell)		!number of clusters of this type
		Diff=findDiffusivity(defectType)		!diffusivity of clusters of this type

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
!> Function findReactionRateSinkFine(CascadeCurrent, defectType, cell, reactionParameter)
!finds reaction rate for defects to get absorbed at sinks (typically matrix dislocations) in the fine mesh (cascade)
!***************************************************************************************************
double precision function findReactionRateSinkFine(CascadeCurrent, defectType, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType(SPECIES), num, grainNum
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff
	type(cascade), pointer :: CascadeCurrent
	double precision findDiffusivity

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_constants
			use mod_globalVariables
			integer cell, defectType(SPECIES)
			type(cascade), pointer :: CascadeCurrent
		end function
	end interface

	if(polycrystal=='yes') then
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	else
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	end if

	!sink reaction function type=12
	if(reactionParameter%functionType==12) then	!sinkRemoval

		num=findNumDefectFine(CascadeCurrent, defectType,cell)		!number of clusters of this type
		Diff=findDiffusivity(defectType)							!diffusivity of clusters of this type

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
!> Function findReactionRateMultiple(defectType1, defectType2, cell, reactionParameter)
!finds reaction rate for defect clustering
!***************************************************************************************************
double precision function findReactionRateMultiple(defectType1, defectType2, cell, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType1(SPECIES), defectType2(SPECIES), i, count
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff1, Diff2, vol
	integer size1, size2, num1, num2
	integer findNumDefect, findDefectSize, grainNum
	double precision findDiffusivity
	double precision area

	if(polycrystal=='yes') then
		grainNum=myMesh(cell)%material
	else
		grainNum=myMesh(cell)%material
	endif

	size1=findDefectSize(defectType1)
	size2=findDefectSize(defectType2)
	Diff1=findDiffusivity(defectType1)
	Diff2=findDiffusivity(defectType2)
	num1=findNumDefect(defectType1,cell)
	num2=findNumDefect(defectType2,cell)

	count=0
	do i=1,SPECIES
		if(defectType1(i)==defectType2(i)) then
			count=count+1
		end if
	end do
	if(count==SPECIES) then
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

	if(reactionParameter%functionType==21) then	!3D-3D

		if((defectType1(3)>0 .AND. defectType1(3) <= max3DInt) .AND. &
				(defectType2(3)>0 .AND. defectType2(3) <= max3DInt)) then	!3D+3D: 3D(SIA) + 3D(SIA)
			reactionRate=Zint*omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*&
					dble(num1)*dble(num2)*atomSize/vol
		else
			reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
					*atomSize/vol
		end if
		!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
		!			 *atomSize/vol

	else if(reactionParameter%functionType==22) then	!3D-1D: Cu/V/CuV + 1D(SIA)

		!reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
		!		*atomSize/vol
		if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

			!if the first defect is the 1D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
			size1=findDefectSize(defectType2)
			size2=findDefectSize(defectType1)
			Diff1=findDiffusivity(defectType2)
			Diff2=findDiffusivity(defectType1)
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

		reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1)*dble(num2)&
				*atomSize/vol+(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0))**4d0*&
						Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)

	else if(reactionParameter%functionType==23) then	!3D-1D: 3D(SIA) + 1D(SIA)

		if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

			!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
			size1=findDefectSize(defectType2)
			size2=findDefectSize(defectType1)
			Diff1=findDiffusivity(defectType2)
			Diff2=findDiffusivity(defectType1)
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

		reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1)*dble(num2)&
				*atomSize/vol+(Zint*(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
						Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)
		!reactionRate=Zint*(omegastar+(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0)))*&
		!		Diff1*num1*num2*atomSize/vol+&
		!		(Zint*(omegastar1D+omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0)))**4d0*&
		!				Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)

	else if(reactionParameter%functionType==24) then	!1D-1D: 1D(SIA) + 1D(SIA)

		reactionRate=(Zint*omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0)))**4d0*&
				(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomSize/vol)**(2d0)

	else
		write(*,*) 'error clustering function type only admits 21~24'
		reactionRate=0d0
	endif

	findReactionRateMultiple=reactionRate

end function

!***************************************************************************************************
!> Function findReactionRateMultipleFine(CascadeCurrent, defectType1, defectType2, cell, reactionParameter)
!finds reaction rate for defect clustering in the fine mesh (Cascade)
!***************************************************************************************************
double precision function findReactionRateMultipleFine(CascadeCurrent,defectType1,defectType2,cell,reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell, defectType1(SPECIES), defectType2(SPECIES), i, count
	type(reactionParameters) :: reactionParameter
	double precision reactionRate, Diff1, Diff2, size1, size2, num1, num2, vol
	!integer findNumDefectFine
	type(cascade), pointer :: CascadeCurrent
	integer findDefectSize
	double precision findDiffusivity
	integer grainNum

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_constants
			use mod_globalVariables
			integer cell, defectType(SPECIES)
			type(cascade), pointer :: CascadeCurrent
		end function
	end interface

	if(polycrystal=='yes') then
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	else
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	end if

	size1=findDefectSize(defectType1)
	size2=findDefectSize(defectType2)
	Diff1=findDiffusivity(defectType1)
	Diff2=findDiffusivity(defectType2)
	num1=findNumDefectFine(CascadeCurrent, defectType1,cell)
	num2=findNumDefectFine(CascadeCurrent, defectType2,cell)
	vol=cascadeElementVol
	!Adjust the number of defects in the reaction if both reactants are of the same type

	count=0
	do i=1,SPECIES
		if(defectType1(i)==defectType2(i)) then
			count=count+1
		end if
	end do
	if(count==SPECIES) then
		!we have two defects of the same type, have to modify the defect numbers for a defect to combine with itself
		num2=num2-1
	end if

	!list of clustering reaction functional forms
	if(reactionParameter%functionType==21) then	!3D-3D

		if((defectType1(3)>0 .AND. defectType1(3) <= max3DInt) .AND. &
				(defectType2(3)>0 .AND. defectType2(3) <= max3DInt)) then	!3D-3D: 3D(SIA) + 3D(SIA)
			reactionRate=Zint*omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*&
					dble(num1)*dble(num2)*atomSize/vol
		else
			reactionRate=omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0))*(Diff1+Diff2)*dble(num1)*dble(num2)&
					*atomSize/vol
		end if

		!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
		!			 *atomSize/vol

	else if(reactionParameter%functionType==22) then	!3D-1D: Cu/V/CuV + 1D(SIA)

		!reactionRate=Ztemp*(omegastar+omega*(dble(size1)**(1d0/3d0)+dble(size2)**(1d0/3d0)))*(Diff1+Diff2)*dble(num1)*dble(num2)&
		!			 *atomSize/vol
		if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

			!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
			size1=findDefectSize(defectType2)
			size2=findDefectSize(defectType1)
			Diff1=findDiffusivity(defectType2)
			Diff2=findDiffusivity(defectType1)
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

		reactionRate=(omega*dble(size1)**(1d0/3d0)+omega2D*dble(size2)**(1d0/2d0))*Diff1*dble(num1*num2)*&
				atomSize/vol+(omegacircle1D*dble(size2)**(1d0/2d0)+omega1D*dble(size1)**(1d0/3d0))**4d0*&
						Diff2*dble(num2)*dble(num1)**(2d0)*(atomSize/vol)**(2d0)


	else if(reactionParameter%functionType==23) then	!3D-1D: 3D(SIA) + 1D(SIA)

		if(defectType1(3) > max3DInt .OR. defectType1(4) > max3DInt) then

			!if the first defect is the 1-D diffusing loop, we have to switch the order of the parameters in order to have the correct reaction rate.
			size1=findDefectSize(defectType2)
			size2=findDefectSize(defectType1)
			Diff1=findDiffusivity(defectType2)
			Diff2=findDiffusivity(defectType1)
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

	else if(reactionParameter%functionType==24) then	!1D-1D: 1D(SIA) + 1D(SIA)

		reactionRate=(Zint*omegacircle1D*(dble(size1)**(1d0/2d0)+dble(size2)**(1d0/2d0)))**4d0*&
				(Diff1*dble(num2)+Diff2*dble(num1))*dble(num1*num2)*(atomSize/vol)**(2d0)
	else
		write(*,*) 'error clustering function type only admits 6-9'
		reactionRate=0d0
	end if

	findReactionRateMultipleFine=reactionRate

end function

!***************************************************************************************************
!> Function findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
!finds reaction rate for defect diffusion between elements
!***************************************************************************************************
double precision function findReactionRateDiff(defectType, cell1, proc1, cell2, proc2, dir, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell1, proc1, cell2, proc2, defectType(SPECIES), num1, num2, dir
	type(ReactionParameters) :: reactionParameter
	double precision Diff, area1, area2, areaShared, lengthShared, Vol1, Vol2, length1, length2, reactionRate
	integer findNumDefect, findNumDefectBoundary
	double precision findStrainEnergy, findStrainEnergyBoundary
	double precision findDiffusivity
	integer grainNum, matNeighbor, size
	double precision Eb, findBinding
	double precision alpha

	grainNum=myMesh(cell1)%material
	if(reactionParameter%functionType==13) then	!3D diffusion

		Diff=findDiffusivity(defectType)
		length1=myMesh(cell1)%length
		area1=length1**2d0
		Vol1=myMesh(cell1)%volume
		num1=findNumDefect(defectType,cell1)

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
			else	!proc2 /= proc1
				matNeighbor=myBoundary(cell2,dir)%material
				length2=myBoundary(cell2,dir)%length
				Vol2=myBoundary(cell2,dir)%volume
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
				alpha=1d0
				num2=0	!If we are diffusing between two material types, assume perfect sink
			end if

			reactionRate=alpha*Diff*areaShared*(dble(num1)/Vol1-dble(num2)/Vol2)/(length1/2d0 + length2/2d0)

			if(reactionRate > 0d0) then
				findReactionRateDiff=reactionRate
			else
				findReactionRateDiff=0d0
			end if
		end if
	elseif(reactionParameter%functionType==15) then	!2D diffusion on a plane (rather than 3D diffusion in a volume)
		Diff=findDiffusivity(defectType)
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

	else if(reactionParameter%functionType==16) then	!Dissociation from grain boundary into bulk volume element
		!(treated as a diffusion reaction because we are moving between two volume elements, but the rate is given by a dissociation rate)
		if(proc2==proc1) then
			matNeighbor=myMesh(cell2)%material
		else
			matNeighbor=myBoundary(cell2,dir)%material
		end if

		if(grainNum==matNeighbor) then
			findReactionRateDiff=0d0		!no dissociation from grain boundary to itself
		else

			Diff=findDiffusivity(matNeighbor, defectType)		!diffusivity of the defect dissociating from the GB (in the bulk)
			!Diff=findDiffusivity(defectType)			!diffusivity of the defect dissociating from the cluster

			num1=findNumDefect(defectType,cell1)			!number of clusters
			size=1										!not breaking up a cluster, but releasing from grain boundary
			Eb=findBinding(defectType,defectType)	!binding energy of defect to grain boundary
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
	end if

end function

!***************************************************************************************************
!> Function findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, reactionParameter)
!finds reaction rate for defect diffusion between a coarse mesh element and a fine (cascade) mesh
!***************************************************************************************************
double precision function findReactionRateCoarseToFine(defectType, cell, proc, numDefectsFine, reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer defectType(SPECIES), cell, proc, numDefectsFine, num1, num2, grainNum
	type(ReactionParameters) :: reactionParameter

	double precision Diff, areaShared, Vol1, Vol2, length1, reactionRate, coarseToFineLength
	integer findNumDefect, findNumDefectBoundary
	double precision findDiffusivity

	if(polycrystal=='yes') then
		grainNum=myMesh(cell)%material
	else
		grainNum=myMesh(cell)%material
	end if

	if(reactionParameter%functionType==13) then	!3D diffusion

		length1=myMesh(cell)%length			!Length of coarse mesh element

		!Effective length used for diffusion from coarse to fine mesh (assuming fine mesh randomly located
		!within the coarse mesh element, see supporting documents for derivation)
		!NOTE: assuming cubic fine mesh (numxCascade=numyCascade=numzCascade). Can re-derive for non-cubic fine meshes.
		CoarseToFineLength=(length1-numxCascade*fineLength)/(dlog(length1**2d0/(numxCascade*fineLength)**2d0))

		Diff=findDiffusivity(defectType)
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
	end if

end function

!***************************************************************************************************
!> Function findReactionRateDiffFine(CascadeCurrent,defectType,cell1,proc1,cell2,proc2,dir,reactionParameter)
!finds reaction rate for defect diffusion between elements in the fine mesh (inside a cascade)
!***************************************************************************************************
double precision function findReactionRateDiffFine(CascadeCurrent,defectType,cell1,proc1,cell2,proc2,dir,reactionParameter)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer cell1, proc1, cell2, proc2, defectType(SPECIES), num1, num2, dir, coarseCell
	type(ReactionParameters) :: reactionParameter
	type(cascade), pointer :: CascadeCurrent
	double precision Diff, area1, area2, areaShared, Vol1, Vol2, length,length1, length2, reactionRate, coarseLength
	double precision fineToCoarseLength, coarseVolume
	integer findNumDefect, grainNum
	double precision findDiffusivity

	interface
		integer function findNumDefectFine(CascadeCurrent, defectType, cell)
			use mod_constants
			use mod_globalVariables
			type(cascade), pointer :: CascadeCurrent
			integer cell, defectType(SPECIES)
		end function
	end interface

	if(polycrystal=='yes') then
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	else
		grainNum=myMesh(CascadeCurrent%cellNumber)%material
	end if

	if(reactionParameter%functionType==13) then	!3D diffusion

		Diff=findDiffusivity(defectType)		!function in MaterialInput

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
!subroutine findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
!points reactionUpdate at the reaction in coarse or fine mesh with matching reactants and products in cell
!If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end of the list.
!***************************************************************************************************
subroutine findReactionInList(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	type(reaction), pointer :: reactionUpdate, reactionPrev
	integer, allocatable :: reactants(:,:), products(:,:)
	integer numReactants, numProducts, count(numReactants+numProducts), i, j, cell
	logical flag

	outer: do while(associated(reactionUpdate))

		if(reactionUpdate%numReactants==numReactants .AND. reactionUpdate%numProducts==numProducts) then

			do i=1,numReactants
				count(i)=0
				do j=1,SPECIES
					if(reactionUpdate%reactants(j,i)==reactants(j,i)) then
						count(i)=count(i)+1
					end if
				end do
			end do

			do i=1,numProducts
				count(i+numReactants)=0
				do j=1,SPECIES
					if(reactionUpdate%products(j,i)==products(j,i)) then
						count(i+numReactants)=count(i+numReactants)+1
					end if
				end do
			end do
			flag=.FALSE.
			inter: do i=1,numReactants+numProducts
				if(count(i) /= SPECIES) then
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
!>Subroutine findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)
!Points reactionUpdate at the correct diffusion reaction in coarse or fine mesh with matching reactants and products in cells
!If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end of the list.
!***************************************************************************************************
subroutine findReactionInListDiff(reactionUpdate, reactionPrev, reactants, cell1, cell2, proc1, proc2)
	use mod_constants
	use mod_globalVariables
	use mod_structures
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
					do j=1,SPECIES
						if(reactionUpdate%reactants(j,1)==reactants(j,1)) then
							count(1)=count(1)+1
						end if
					end do

					count(2)=0
					do j=1,SPECIES
						if(reactionUpdate%products(j,1)==reactants(j,1)) then
							count(2)=count(2)+1
						endif
					end do

					flag=.FALSE.
					do i=1,2
						if(count(i) /= SPECIES) then
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
!>Subroutine findReactionInListMultiple(reactionUpdate, reactionPrev, cell, reactants, products, numReactants, numProducts)
!Points reactionUpdate at the clustering reaction in the coarse or fine mesh with matching reactants and products in cell
!If reaction is not present, reactionUpdate is not associated and reactionPrev points to the end of the list.
!***************************************************************************************************
subroutine findReactionInListMultiple(reactionUpdate,reactionPrev,cell,reactants,products,numReactants,numProducts)
	use mod_constants
	use mod_globalVariables
	use mod_structures
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
			do j=1,SPECIES
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
				do j=1,SPECIES
					if(reactionUpdate%products(j,i)==products(j,i)) then
						count(i+numReactants)=count(i+numReactants)+1
					end if
				end do
			end do

			flag=.FALSE.
			inter1: do i=1,numReactants+numProducts
				if(count(i) /= SPECIES) then
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
			do j=1,SPECIES
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
				do j=1,SPECIES
					if(reactionUpdate%products(j,i)==products(j,i)) then
						count(i+numReactants)=count(i+numReactants)+1
					end if
				end do
			end do

			flag=.FALSE.
			inter2: do i=1,numReactants+numProducts
				if(count(i) /= SPECIES) then
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
!>subroutine defectCombinationRules(products, product2, defectTemp, isCombined)
!Hard-coded information on what happens to defects when they cluster (for example, annihilation,
!formation of sessile SIA clusters, etc).
!NOTE: defect combination rules are primarily input directly into subroutines addMultiDefectReactions
!and addMultiDefectReactionsFine, and this subroutine is only called in order to get correct
!defect combination in the case of cascade-defect interactions. In the future, may want to
!move all defect combination rules to this subroutine so that they only need to be changed once.
!***************************************************************************************************
subroutine defectCombinationRules(products, product2, defectTemp, isCombined)
	use mod_constants
	use mod_structures
	use mod_globalVariables
	implicit none

	integer products(SPECIES), product2(SPECIES)
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

		do l=1,SPECIES
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
! Subroutine checkReactionLegality(numProducts, products, isLegal)
! This subroutine looks at the products of a combination reaction and checks to see if the reaction
! is allowed (using hard-coded information). If not, the subroutine returns a value of .FALSE. to isLegal.
!***************************************************************************************************
subroutine checkReactionLegality(numProducts, products, isLegal)
	use mod_constants
	use mod_globalVariables
	use mod_structures
	implicit none

	integer numProducts, i
	integer products(SPECIES,numProducts)
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

end module

