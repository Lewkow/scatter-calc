
SUBROUTINE diff_cross3d_MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write to a file, DC_3D.dat, with the format
!				( | E | Theta | Amp | )
!
!	Requires Keys have DC3D set to on to run. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info, ONLY : rank, ierr
	USE collision_inputs, ONLY : NTh, ThI, ThF, FRAME
	USE physics_constants, ONLY : TOEV, TODEG
	USE dcs, ONLY : Amp
	USE current_energy, ONLY : E
	IMPLICIT NONE
	INCLUDE 'mpif.h'

	!--------------------
	! Internal	
	!--------------------

	INTEGER					:: i

	REAL(KIND=8)		:: T
	REAL(KIND=8)		:: DumT
	REAL(KIND=8)		:: DT

	IF (rank == 0) THEN	
		OPEN(UNIT=52, FILE="../Data/DC_3D.dat", ACCESS="APPEND")
	END IF

	DT = (ThF - ThI)/DBLE(NTh-1)

	DO i=1,NTh
		T = ThI + (i-1)*DT	
		CALL pw_amp_MPI(T)
		IF (rank == 0) THEN	
			IF (FRAME == 1) THEN
				CALL amp_to_lab(T,DumT)
				T = DumT
			END IF
			WRITE(52,*) E*TOEV, T*TODEG, Amp
!			WRITE(*,*) E*TOEV, T*TODEG, Amp
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	END DO

	IF (rank == 0) THEN	
		CLOSE(52)
	END IF

END SUBROUTINE
