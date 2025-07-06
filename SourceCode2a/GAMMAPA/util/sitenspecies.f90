MODULE sitenspecies
    TYPE siteinfo
        INTEGER	:: host = 0	    ! gmu component index for site host. Currently each site on only one host
        CHARACTER(20) :: name = ''    !site name
        INTEGER	:: id = 0	! site id value (e.g. 2 or 502)
        INTEGER	:: noccur = 0	! for use when each site on only one host
        INTEGER	:: packed = 0	! packed id
        REAL*8	:: xhost = 0D0	! mole fraction of host
    END TYPE siteinfo

    TYPE species
        CHARACTER(50):: name
        INTEGER :: id
        LOGICAL	:: selfassoc = .FALSE.	! Flase by default, True if self associating
        LOGICAL	:: dnr = .FALSE.	! Flase by default, True if electron donors on species
        LOGICAL	:: acpt = .FALSE.	! Flase by default, True if electron acceptors on species
        INTEGER	:: nsite = 0		! number of sites on species
        INTEGER	:: sid(5) = (/ 0,0,0,0,0 /)	! site (defined above) indexes to point to the sites on the species, vector size limits the # sites on a species
        REAL*8	:: x = 0D0 		! mole fraction
    END TYPE species

    type (siteinfo), dimension(:), allocatable :: site ! dimensioned later to be the number of sites present
    type (species), dimension(:), allocatable :: comp ! components present
    real*8, dimension(:), allocatable :: x
    real*8, dimension(:,:), allocatable :: KAD, eps
    real*8, dimension(:), allocatable :: PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok ! PCSAFT
    integer :: nc, nsites, aspmx
    ! aspmx =1 if association or solvation is present in mixture (or pure components)

END MODULE sitenspecies
