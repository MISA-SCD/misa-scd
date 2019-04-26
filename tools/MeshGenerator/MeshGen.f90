! $Header: /home/CVS//srscd/tools/MeshGenerator/MeshGen.f90,v 1.5 2015/11/12 18:04:58 rdingre Exp $
!***************************************************************************************************
!
! This is the mesh generator program that will output mesh text files that are
! usable in SRSCD_par. The goal of this program is to have an easy and quick way to create uniform
! meshes for SRSCD without changing the (more demanding) inputs in SRSCD.
!
! Inputs: Volume element length (uniform volume elements only), meshType (periodic or free surfaces),
! 		number of elements in x, y, and z-directions.
! 
! Output: Text file that can be read into SRSCD_par with mesh. Each element's center coordinates 
!		as well as the global max/min coordinates in the x, y, and z directions will be given.
!		All materials will be of type 1.
!
!***************************************************************************************************

program MeshGen
implicit none

integer numx, numy, numz, i, j, k, numGrains, grain, n
integer cell, numElements
double precision length
character(20) meshType, temp, filename, iMPALEtoggle, iMPALEfilename
character(20) nanoporous
double precision poreradius, ligamentsize, lengthX, lengthY, lengthZ
double precision lengthUnitCell, centerX, centerY, centerZ
double precision Xunitcell, Yunitcell, Zunitcell, radiustemp
integer numcellX, numcellY, numcellZ, numelementcell
logical flag

!Default value
numGrains=1

!Input file where information on uniform cubic mesh is stored
open(81, file='MeshGenInput.txt',action='read', status='old')

!Part 1: Read in parameters from input file

iMPALEtoggle='no'
nanoporous='no'

flag=.FALSE.

do 10 while(flag .eqv. .FALSE.)

	read(81,*) temp
	
	if(temp=='meshType') then
		read(81,*) meshType
	else if(temp=='length') then
		read(81,*) length
	else if(temp=='numGrains') then
		read(81,*) numGrains
	else if(temp=='numx') then
		read(81,*) numx
	else if(temp=='numy') then
		read(81,*) numy
	else if(temp=='numz') then
		read(81,*) numz
	else if(temp=='filename') then
		read(81,*) filename
	else if(temp=='iMPALEtoggle') then
		read(81,*) iMPALEtoggle
	else if(temp=='iMPALEfilename') then
		read(81,*) iMPALEfilename
	else if(temp=='nanoporous' ) then
	    read(81,*) nanoporous
	else if(temp=='end') then
		flag=.TRUE.
	endif
	
10 continue

if(iMPALEtoggle=='yes') then
	open(83, file=iMPALEfilename, action='read', status='old')
endif

if(nanoporous=='yes') then
	flag=.FALSE.

	do 20 while(flag .eqv. .FALSE.)
	
		read(81,*) temp
		if (temp=='nanoporousStart') then
			flag=.TRUE.
		endif
		
	20 continue
	
	flag=.FALSE.
	
	do 30 while(flag .eqv. .FALSE.)
	
		read(81,*) temp
		
		if(temp=='poreradius') then
			read(81,*) poreradius
		else if(temp=='ligamentsize') then
			read(81,*) ligamentsize		
		else if(temp=='lengthX') then
			read(81,*) lengthX
		else if(temp=='lengthY') then
			read(81,*) lengthY
		else if(temp=='lengthZ') then
			read(81,*) lengthZ
		else if(temp=='end') then
			flag=.TRUE.
		endif
		
	30 continue 
	! Number of repeating unit cells in each direction
	numcellX = INT(lengthX/(2*poreradius+ligamentsize))
	numcellY = INT(lengthY/(2*poreradius+ligamentsize))
	numcellZ = INT(lengthZ/(2*poreradius+ligamentsize))
	! Number of elements in unit cell
	numelementcell = INT((2*poreradius+ligamentsize)/length)
	! Redefine number of volume elements in the entire mesh
	numx = numelementcell*numcellX
	numy = numelementcell*numcellY
	numz = numelementcell*numcellZ
endif

!Part 2: Write output file

open(82, file=filename, action='write', status='unknown')

!First put in comment section

write(82,*) '!*************************************************************************************'
write(82,*) '!This file contains a UNIFORM cubic mesh input for the parallel SRSCD'
write(82,*) '!code. The mesh will be read in via coordinates (x,y,z) of each element center'
write(82,*) '!as well as the material number of each element (eg. material 1=Fe, 2=Cu, etc)'
write(82,*) '!'
write(82,*) '!The first values that will be read in are the number of elements in each direction'
write(82,*) '!(x,y,z).'
write(82,*) '!'
write(82,*) '!*************************************************************************************'
write(82,*)

!Next specify meshType

write(82,*) '!Specify mesh type (currently two types=periodic and freeSurfaces)'
write(82,*) 'meshType'
write(82,*) meshType
write(82,*)

!Next specify uniform cubic mesh element length

write(82,*) '!Specify element lengths (nm)'
write(82,*) 'length'
write(82,*) length
write(82,*)

!Next specify the min and max of each coordinate

write(82,*) '!Specify max and min of each coordinate (nm)'
write(82,*) '!NOTE: the actual max, min will be max, min +/- length/2 because these are assumed'
write(82,*) '!to be the center of the elements and not the edges'
write(82,*) 'xminmax' 
write(82,*) length/2d0, length/2d0+(numx-1)*length
write(82,*) 'yminmax'	
write(82,*) length/2d0, length/2d0+(numy-1)*length
write(82,*) 'zminmax'	
write(82,*) length/2d0, length/2d0+(numz-1)*length
write(82,*)

!Next specify numx, numy, numz

write(82,*) '!Specify number of elements in each direction'
write(82,*) 'numx'	
write(82,*) numx
write(82,*) 'numy'
write(82,*) numy
write(82,*) 'numz'	
write(82,*) numz
write(82,*)

!Next create a list of coordinates of centers of volume elements

write(82,*) '!Coordinates of element centers (x, y, z, elementType) (nm)'
write(82,*) '!NOTE: these must be ordered in the same way as the connectivity matrix (loop x, then y, then z)'
write(82,*) 'elements			!tells computer that the following values are elements'

if(iMPALEtoggle=='no'.AND.nanoporous=='no') then

	grain=1
	n=1
	
	do 11 i=1,numz
		if(dble(i) .GT. dble(n*numz)/dble(numGrains)) then
			grain=grain+1
			n=n+1
		endif
		do 12 j=1,numy
			do 13 k=1,numx
				write(82,*) (k-1)*length+length/2d0, (j-1)*length+length/2d0, (i-1)*length+length/2d0, grain
			13 continue
			write(82,*)
		12 continue
	11 continue

else if(iMPALEtoggle=='yes') then

	read(83,*) numElements
	read(83,*)
	
	if(numElements .NE. numx*numy*numz) then
		write(*,*) 'Error numElements in ', iMPALEfilename, 'does not match MeshGenInput.txt'
	endif
	
	do 14 cell=1,numElements
		read(83,*) k,j,i,grain
		write(82,*) (k-1)*length+length/2d0, (j-1)*length+length/2d0, (i-1)*length+length/2d0, grain
	14 continue

else if(nanoporous=='yes') then

	numGrains = 2
	lengthUnitCell = numelementcell*length

	do 31 i=1,numz
		do 32 j=1,numy
			do 33 k=1,numx
				grain = 1
! Center of volume element
				centerX = (k-1)*length+length/2d0
				centerY = (j-1)*length+length/2d0
				centerZ = (i-1)*length+length/2d0
				Xunitcell = dmod(centerX,lengthUnitCell)
				Yunitcell = dmod(centerY,lengthUnitCell)
				Zunitcell = dmod(centerZ,lengthUnitCell)
				radiustemp = dsqrt((lengthUnitCell/2.d0 - Xunitcell)**2.d0 &
						   + (lengthUnitCell/2.d0 - Yunitcell)**2.d0 &
						   + (lengthUnitCell/2.d0 - Zunitcell)**2.d0)
				if(radiustemp.le.poreradius) grain=2
				write(82,*) (k-1)*length+length/2d0, (j-1)*length+length/2d0, (i-1)*length+length/2d0, grain
			33 continue
			write(82,*)
		32 continue
	31 continue

else

	write(*,*) 'Error unrecognized toggle'
	
endif

end program
