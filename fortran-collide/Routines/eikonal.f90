
SUBROUTINE eikonal
	USE collision_inputs, ONLY : TCS, DCS

	IF (DCS == 1) THEN
		CALL eik_diff_cross_section
	END IF

	IF (TCS == 1) THEN
!		CALL eik_tot_cross_section
	END IF

END SUBROUTINE eikonal
