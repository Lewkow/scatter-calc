
REAL FUNCTION randy(seed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!	randy returns a random number uniformly distributed between
! [0,1]. This is the Park & Miller minimal standard for 
! a random number generator as stated in the book 
! "Exploring Monte Carlo Methods" by Dunn and Shultis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

	IMPLICIT NONE

	INTEGER							:: temp
	INTEGER							:: seed
	INTEGER,PARAMETER		:: a = 16807
	INTEGER,PARAMETER		:: m = 2147483647
	INTEGER,PARAMETER		:: q = 127773
	INTEGER,PARAMETER		:: r = 2836
	
	REAL,PARAMETER			:: minv = 1./m
	
	temp = seed/q
	seed = a*(seed-temp*q) - r*temp
	
	if (seed < 0 ) seed = seed + m
	
	randy = minv*seed

	return

END FUNCTION randy

