MODULE CONSTANTS
	REAL*8, PARAMETER :: avoNum = 602.214076d0, kB = 0.01380649D0, R = avoNum*kB
	INTEGER, PARAMETER :: outfile = 2000 ! unit number for regular output
	INTEGER, PARAMETER :: debugfile = 2001 ! unit for debug file output
	! Units 10, 11, 12 are used for gamma.csv, hxs.csv, vol.csv respectively, within code
	CHARACTER(4), PARAMETER :: missing = '1D35'
END MODULE CONSTANTS
