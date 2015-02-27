
SUBROUTINE eik_phase( b, p )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the phaseshift for the Eikonal method
! by integrating over the range of action for
! a given impact parameter b. Computed serially 
! since each rank will compute for different
! impact parameters.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : K
	USE collision_inputs, ONLY : Mu
	USE grid, ONLY : REnd
	USE eik, ONLY : PhaseInt

	IMPLICIT NONE

	!--------------------------
	! Inputs	
	!--------------------------

	REAL(KIND=8)					:: b     ! Impact parameter

	!--------------------------
	! Outputs	
	!--------------------------

	REAL(KIND=8) 					:: p     ! Phaseshift
	
	!--------------------------
	! Internal	
	!--------------------------

	REAL(KIND=8)					:: delta ! Integration step
	REAL(KIND=8)					:: grand ! Ingegration sum value
	REAL(KIND=8)					:: pot   ! Value of potential  
	REAL(KIND=8)					:: z     ! Current cylindrical dist to core 
	REAL(KIND=8)					:: d     ! Current dist to core
	REAL(KIND=8)					:: ZMax  ! Maximum Z value 
	REAL(KIND=8)					:: coeff ! coeff which changes with integration

	INTEGER								:: i     ! Counter

	IF (b > REnd) THEN
		p = 0.0
	ELSE
		ZMax  = SQRT(REnd**2.0 - b**2.0)
		coeff = -Mu/(2.0*K)
		delta = 2.0*ZMax/DBLE(PhaseInt)
		grand = 0.0
	
		DO i=0,PhaseInt
			z = -ZMax + i*delta
			d = SQRT(b**2.0 + z**2.0)
			CALL potential(d,pot)
			grand = grand + pot
		END DO

		p = grand*coeff*delta

	END IF			

END SUBROUTINE eik_phase


