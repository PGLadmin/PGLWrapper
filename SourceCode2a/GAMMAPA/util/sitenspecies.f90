MODULE sitenspecies
        TYPE siteinfo
                INTEGER	:: host = 0	    ! gmu component index for site host. Currently each site on only one host
                CHARACTER(8) :: name = ''    !site name
                INTEGER	:: idaspen = 0	! aspen group id value (e.g. 4500)
                INTEGER	:: noccur = 0	! for use when each site on only one host
                INTEGER	:: packed = 0	! packed id
                REAL*8	:: xhost = 0D0	! mole fraction of host
        END TYPE siteinfo

        TYPE species
                LOGICAL	:: selfassoc = .FALSE.	! Flase by default, True if self associating
                LOGICAL	:: dnr = .FALSE.	! Flase by default, True if donors on species
                LOGICAL	:: acpt = .FALSE.	! Flase by default, True if acceptors on species
                INTEGER	:: nsite = 0		! number of sites on species
                INTEGER	:: sid(5) = (/ 0,0,0,0,0 /)	! site (defined above) indexes to point to the sites on the species, vector size limits the # sites on a species
                REAL*8	:: x = 0D0 		! mole fraction
        END TYPE species
END MODULE sitenspecies
