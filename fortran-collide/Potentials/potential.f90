
SUBROUTINE potential(R,pot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potential wrapper function which call appropriate 
! potential function based on reduced mass of system
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : Proj, Targ
	USE physics_constants, ONLY : TOAMU
	USE current_energy, ONLY : Channel

	IMPLICIT NONE

	!----------------------
	! Inputs	
	!----------------------

	REAL(KIND=8)					:: R

	!----------------------
	! Outputs	
	!----------------------

	REAL(KIND=8)					:: pot
	REAL(KIND=8)					:: dvx
	REAL(KIND=8)					:: d2vx

	IF ( ( (Proj .EQ. 'He4') .OR. (Proj .EQ. 'He3') ) .AND. ( (Targ .EQ. 'He4') .OR. (Targ .EQ. 'He3') ) ) THEN
		CALL pot_HeHe(R,pot)

	ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'H') .OR. (Proj .EQ. 'H' .AND. Targ .EQ. 'He4') ) THEN
		CALL pot_HeH(R,pot)

	ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN
		IF ( Channel .EQ. 'S' ) THEN
			CALL pot_HeO_Pi(R,pot)
		ELSE IF ( Channel .EQ. 'T' ) THEN
			CALL pot_HeO_Sg(R,pot)
		ELSE
			WRITE(*,*) "Channel ", Channel, " is not available!!!"
			WRITE(*,*) "Currently only Singlet PI (S) and Triplet SIGMA (T) channels available"
			STOP
		END IF
			
	ELSE IF ( (Proj .EQ. 'H') .AND. (Targ .EQ. 'H') ) THEN
		IF ( Channel .EQ. 'S' ) THEN
			CALL pot_HH_Singlet(R,pot,dvx,d2vx,1)
		ELSE IF ( Channel .EQ. 'T' ) THEN
			CALL pot_HH_Triplet(R,pot,dvx,d2vx,1)
		ELSE
			WRITE(*,*) "Channel ", Channel, " is not available!!!"
			WRITE(*,*) "Currently only Singlet (S) and Triplet (T) channels available"
			STOP
		END IF
	ELSE
		WRITE(*,*) "Inputs do not match any known potential interaction"
		WRITE(*,*) "                 !!! Terminating NOW !!!"
		STOP 'Inputs do not match any known potential'

	END IF

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE pot_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Writes the current potential being used to the screen
! based on reduced mass of system
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : Proj, Targ
	USE physics_constants, ONLY : TOAMU

	IMPLICIT NONE

	
	WRITE(*,*) "--------------------------------------------------"

  IF ( ( (Proj .EQ. 'He4') .OR. (Proj .EQ. 'He3') ) .AND. ( (Targ .EQ. 'He4') .OR. (Targ .EQ. 'He3') ) ) THEN
		WRITE(*,*) "Potential being used: He+He"

  ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'H') .OR. (Proj .EQ. 'H' .AND. Targ .EQ. 'He4') ) THEN
		WRITE(*,*) "Potential being used: He+H"

  ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN
		WRITE(*,*) "Potential being used: He+O"

  ELSE IF ( (Proj .EQ. 'H') .AND. (Targ .EQ. 'H') ) THEN
		WRITE(*,*) "Potential being used: H+H"

  ELSE
    WRITE(*,*) "Inputs do not match any known potential interaction"
    WRITE(*,*) "                 !!! Terminating NOW !!!"
    STOP 'Inputs do not match any known potential'

  END IF

END SUBROUTINE
