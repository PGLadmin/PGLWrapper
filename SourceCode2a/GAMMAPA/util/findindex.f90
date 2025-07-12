!****************************************************
! Subroutine FINDINDEX
!
! Find the component index for a given component id in comp structure
!
!***************************************************
SUBROUTINE FINDINDEX(index,compid)

use sitenspecies, only: nc, comp
integer index, compid

do i = 1, nc
 if(comp(i)%id .EQ. compid) EXIT
END DO
index = i
END SUBROUTINE FINDINDEX
!*********************************************************
