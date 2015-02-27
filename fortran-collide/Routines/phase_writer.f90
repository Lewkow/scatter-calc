
SUBROUTINE phase_end_writer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Writes the number of phases required to converge 
! as a function of collision energy and reduced mass
! to a file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : E, LMax
	USE collision_inputs, ONLY : MU
	USE physics_constants, ONLY : TOEV

	IMPLICIT NONE

	OPEN(UNIT=666, FILE='../Data/PhaseEnds.dat', ACCESS='APPEND')

	WRITE(666,*) MU, E*TOEV, LMax

	CLOSE(666)

END SUBROUTINE phase_end_writer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE phase_writer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Writes the calculated phaseshifts to an output 
!	file. This routine is only called if PHASE_WRITE
! is set to 1 in the SETTINGS.in input file. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : phase, E, LMax	
	USE physics_constants, ONLY : TOEV

	IMPLICIT NONE

	INTEGER													:: L
	CHARACTER(LEN=28)								:: Filename		! Filename to print data to

	WRITE(unit=Filename,fmt="(A, ES7.1, A)") "../Data/Phases_", E*TOEV, "ev.dat"	

	OPEN(UNIT=666, FILE=Filename, ACCESS="APPEND")
	
	DO L=1,LMax
		WRITE(666,*) (L-1), phase(L)
	END DO

	CLOSE(666)

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE eik_phase_writer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates phaseshifts for the current energy
! using the Eikonal Approximation and writes the 
! phaseshifts to a file. This routine is only 
! called if PHASE_WRITE is set to 1 in the 
! SETTINGS.in input file. The output file is of the 
! form ([b , phaseshift]) where 'b' is the impact
! parameter.   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : E
	USE physics_constants, ONLY : TOEV
	
	IMPLICIT NONE

	REAL(KIND=8)										:: b_start		! Initial impact parameter
	REAL(KIND=8)										:: b_end			! Final impact parameter
	REAL(KIND=8)										:: db					! Change in impact parameter
	REAL(KIND=8)										:: b					! Current impact parameter
	REAL(KIND=8)										:: P 					! Current phaseshift

	INTEGER													:: Nb					! Number of impact parameters	
	INTEGER													:: i					! Counter
	
	CHARACTER(LEN=32)								:: Filename		! Filename to print data to


	WRITE(unit=Filename,fmt="(A, ES7.1, A)") "../Data/Eik_Phases_", E*TOEV, "ev.dat"

	OPEN(UNIT=777, FILE=Filename, ACCESS="APPEND")

	Nb      = 1000
	b_start = 0.5
	b_end   = 7.0
	db      = (b_end-b_start)/DBLE(Nb) 		
	
	DO i=0,Nb
		b = b_start + i*db	
		CALL eik_phase(b,P)
		WRITE(777,*) b, P
	END DO

	CLOSE(777)

END SUBROUTINE


