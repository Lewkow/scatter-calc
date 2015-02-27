
SUBROUTINE diff_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Loops through designated range of scattering
! angles and calculates differential cross section
! useing partial wave analysis, prints collision
! energy, angle and differential cross section
! to the file Pw_DCS.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE dcs, ONLY : Amp
	USE collision_inputs, ONLY : ThI, ThF, NTh, FRAME
	USE physics_constants, ONLY : PI, TODEG, TOEV
	USE current_energy, ONLY : E

	IMPLICIT none
	
	INCLUDE 'mpif.h'

	!--------------------------------------
	! Internal	
	!--------------------------------------

	REAL(KIND=8)											:: DT				! dTheta
	REAL(KIND=8)											:: T				! CM Theta
	REAL(KIND=8)											:: T_L			! Lab Theta

	INTEGER														:: i   			! Counter

	CHARACTER(LEN=30)									:: Fname		! Filename to print data to
	CHARACTER(LEN=8)									:: mid

	DT = (ThF - ThI)/DBLE(NTh-1)

	WRITE(unit=Fname,fmt="(A, ES9.3, A)") "../Data/PW_DCS_", E*TOEV, "eV.dat"

	IF (rank == 0) THEN
		OPEN(UNIT=20, FILE=Fname, ACCESS="APPEND")
		WRITE(*,*) "        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        "
		WRITE(*,*) "            Starting DCS Computation             "	
		WRITE(*,*) "        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        "
	END IF

	DO i=1,NTh
		T = ThI + (i-1)*DT
		CALL pw_amp(T)
		IF (rank == 0) THEN
			IF (FRAME == 1) THEN
				CALL angle_to_lab(T, T_L)
				T = T_L
			END IF
			WRITE(20,*) T*TODEG, Amp
		END IF
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	END DO

	CLOSE(20)				

END SUBROUTINE
