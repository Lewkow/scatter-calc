
SUBROUTINE eik_diff_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Loops through designated range of scattering
! angles and calculates differential cross section
! useing the Eikonal approximation, prints collision
! energy, angle and differential cross section
! to the file Eik_DCS.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE dcs, ONLY : Amp
	USE collision_inputs, ONLY : ThI, ThF, NTh, FRAME
	USE physics_constants, ONLY : PI, TODEG, TOEV
	USE current_energy, ONLY : E

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!-------------------------------
	! Internal	
	!-------------------------------

	REAL(KIND=8)				:: DT 	! dTheta
	REAL(KIND=8)				:: T  	! CM Theta
	REAL(KIND=8)				:: T_L  ! Lab Theta

	INTEGER							:: i	! Counter

	DT = (ThF - ThI)/DBLE(NTh)

	IF (rank == 0) THEN
		OPEN(UNIT=30, FILE="../Data/Eik_DCS.dat", ACCESS="APPEND")

		DO i=1,NTh
			T = ThI + (i-1)*DT
			CALL eik_amp(T)
			IF (FRAME == 1) THEN
				CALL angle_to_lab(T, T_L)
				T = T_L
			END IF
			WRITE(30,*) E*TOEV, T*TODEG, Amp
		END DO

		CLOSE(30)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE eik_diff_cross_section
