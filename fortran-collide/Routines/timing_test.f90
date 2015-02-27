
SUBROUTINE phase_timing_test(t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Routine to write phase times to file
! Phase_Timings.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	USE mpi_info, ONLY : numrank
	USE collision_inputs, ONLY : Mu
	USE current_energy, ONLY : E
	USE physics_constants, ONLY : TOEV

	!-------------------------------
	! Inputs	
	!-------------------------------

	REAL(KIND=8)							:: t

	OPEN(UNIT=10, FILE="../Data/phase_timings.dat", ACCESS="APPEND")
	
	WRITE(10,*) numrank, Mu, E*TOEV, t

	CLOSE(10)
END SUBROUTINE phase_timing_test

