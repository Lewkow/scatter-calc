
SUBROUTINE depth_finder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Finds the depth of the potential, or the
! lowest point of the potnetial between 
! Ri and Rf. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : Depth

	!------------------------
	! Internal	
	!------------------------

	REAL(KIND=8)								:: Ri
	REAL(KIND=8)								:: Rf
	REAL(KIND=8)								:: DR
	REAL(KIND=8)								:: R
	REAL(KIND=8)								:: low 
	REAL(KIND=8)								:: pot 

	INTEGER											:: N
	INTEGER											:: i

	
	low = 0.0
	Ri  = 1.0 
	Rf  = 10.0 
	N   = 10000 
	DR  = (Rf - Ri)/DBLE(N)	

	DO i=1,N
		R = Ri + i*DR
		CALL potential(R,pot)
		IF (pot .LE. low) THEN
			low = pot
		END IF
	END DO

	Depth = ABS(low)

END SUBROUTINE	
