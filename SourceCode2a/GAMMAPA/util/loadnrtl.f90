! loadnrtl.f90
!
!**********************************************
!
!  loadnrtl - loads NRTL parameters
!
!**********************************************

SUBROUTINE LOADNRTL()

use CONSTANTS, only: outfile
use sitenspecies, only: comp, nc
USE PHYS_PARMS, ONLY: aparam, bparam, alpha

integer :: start_index, npair, index1, index2, compid1, compid2
 
! local variables
integer :: i, j, id1, id2
character*500 :: data_line
character*100 :: data_field
character*100 :: filenrtl
LOGICAL :: file_exists = .FALSE.

DO WHILE (.NOT. file_exists)
	print *, 'Enter a the name of the tab-delimited NRTL parameter file:'
	read(*,'(A)') filenrtl
	print *, 'You entered:', trim(filenrtl)
	filenrtl = '..\..\Input\GAMMAPA\'//filenrtl
	INQUIRE(FILE=trim(filenrtl), EXIST=file_exists)
	IF(.NOT. file_exists) print *, 'That file is not found. Check folder and spelling.'
	print *, ' '
END DO

WRITE(outfile,*) 'NRTL file: '//trim(filenrtl)
OPEN(UNIT=1001, FILE=trim(filenrtl))
start_index = 1
READ(UNIT=1001, END=106, FMT='(A)') data_line
call parse_line(start_index,data_line,data_field)
read(data_field,*) npair
READ(UNIT=1001, END=106, FMT='(A)') data_line ! Header

DO i = 1, npair
	start_index = 1
	READ(UNIT=1001, END=106, FMT='(A)') data_line
	call parse_line(start_index,data_line,data_field) !pairid
	call parse_line(start_index,data_line,data_field) !id1
	read(data_field,*) compid1
	call parse_line(start_index,data_line,data_field) !name
	call parse_line(start_index,data_line,data_field) !id2
	read(data_field,*) compid2
	call parse_line(start_index,data_line,data_field) !name
	! find index for components
	call findindex(index1,compid1)
	call findindex(index2,compid2)
	call parse_line(start_index,data_line,data_field) !a12
	read(data_field,*) aparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !b12
	read(data_field,*) bparam(index1,index2)
	call parse_line(start_index,data_line,data_field) !a21
	read(data_field,*) aparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !b21
	read(data_field,*) bparam(index2,index1)
	call parse_line(start_index,data_line,data_field) !alpha
	read(data_field,*) alpha(index1,index2)
	alpha(index2,index1) = alpha(index1,index2)
END DO

106 CLOSE(1001)

END SUBROUTINE LOADNRTL