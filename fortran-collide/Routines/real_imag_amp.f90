
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write real and imaginary amplitudes at 0 deg scattering
! angle to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE real_imag_write
	USE mpi_info, ONLY : rank

	IMPLICIT NONE

	IF (rank .EQ. 0) CALL pw_amp(0.0D0)

END SUBROUTINE real_imag_write

