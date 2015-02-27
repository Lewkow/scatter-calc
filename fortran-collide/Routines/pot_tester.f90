
SUBROUTINE pot_tester
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Tests the currently used potential by
! writing the potential from POT_I to POT_F
! using POT_N steps. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : POT_I, POT_F, POT_N
	USE physics_constants, ONLY : TOEV

	IMPLICIT NONE

	!---------------------	
	! Internal	
	!---------------------	

	REAL(KIND=8)					:: pot 	! Current potential value
	REAL(KIND=8)					:: R	 	! Current separation distance
	REAL(KIND=8)					:: DR		! Differential separation distance	

	INTEGER								:: i   	! Counter


	DR = (POT_F - POT_I)/DBLE(POT_N-1)

	OPEN(UNIT=70, FILE="../Data/Pot_Test.dat", ACCESS="APPEND")

	DO i=1,POT_N
		R = POT_I + (i-1)*DR
		CALL potential(R,pot)
		WRITE(70,*) R, pot*TOEV
	END DO
	
	CLOSE(70)

END SUBROUTINE
