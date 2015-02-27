
SUBROUTINE call_tests
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine calls the different tests that
! are specified in the input files. This eliminates
! several lines of code in the main program scatter. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE collision_inputs
	USE current_energy, ONLY : E
	USE physics_constants, ONLY : TOEV

	IMPLICIT NONE

	REAL(KIND=8)		:: d_dummy, TC

	INCLUDE 'mpif.h'

!	CALL reduced_coordinates
!	CALL real_imag_write

 	IF (TCS == 1) THEN
  	IF (CALC_TYPE == 0) THEN 	
			CALL write_total_cross_section
!			CALL uni_TCS( E, Mu, d_dummy )
			CALL total_cross_section(TC)
			WRITE(876,*) E*TOEV, TC	
			CALL TCS_HH( E*TOEV, TC )	
			WRITE(877,*) E*TOEV, TC	
		ELSE IF (CALC_TYPE == 1) THEN
	   	CALL eik_total_cross_section
		ELSE IF (CALC_TYPE == 99) THEN
			CALL write_total_cross_section
	  	CALL eik_total_cross_section
		END IF
	END IF

	IF (NUM_TCS .EQ. 1) THEN
		CALL num_total_cross_section
	END IF

	IF (DFCS == 1) THEN
		CALL diffusion_cross_section
	END IF

  IF (DCS == 1) THEN
		IF (CALC_TYPE == 0) THEN
 	  	CALL diff_cross_section
    ELSE IF (CALC_TYPE == 1) THEN
	  	CALL eik_diff_cross_section
    ELSE IF (CALC_TYPE == 99) THEN
 	  	CALL diff_cross_section
    	CALL eik_diff_cross_section
    END IF
  END IF

  IF (HS_TEST == 1) THEN
    CALL hard_sphere_test
  END IF

  IF ( (PROB_DEN == 1) .AND. (rank == 0) ) THEN
    CALL prob_density_write
!    CALL critical_angle_write
  END IF

	IF ( (PROB_DEN_3D == 1) .AND. (rank == 0) ) THEN
    CALL prob_density_3D_write
	END IF

	IF ( (PROB_DEN_FIT == 1) .AND. (rank == 0) ) THEN
!		CALL prob_den_angle_fitter
		CALL write_angle_probability
	END IF

  IF (DC3D == 1) THEN
    CALL diff_cross3d_MPI
  END IF

	IF (AVE_EN_LOSS == 1) THEN
!		CALL energy_loss
		CALL ave_energy_loss
!		CALL ave_percent_energy_loss
	END IF

END SUBROUTINE

