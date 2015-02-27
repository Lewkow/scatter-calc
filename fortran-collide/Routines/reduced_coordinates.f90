!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Writes reduced coordinates to data file for
! a given energy. Reduced coorinates are given
! in the lab frame and between 0 and max 
! scattering angle for given projectile-target. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE reduced_coordinates
	USE mpi_info 
	USE collision_inputs, ONLY : M1, M2, FRAME
	USE current_energy, ONLY : E
	USE physics_constants, ONLY : TOEV, TOAMU, PI
	USE dcs, ONLY : Amp

	IMPLICIT NONE

	INCLUDE 'mpif.h'
	
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: X, Y
	REAL(KIND=8)	:: Mu, T_in, T_fn, dT, T, T_lab
	INTEGER				:: i, N_t

	IF (rank .EQ. 0) THEN

		OPEN(UNIT=50,FILE='../Data/ReducedCoordinates.dat',ACCESS='APPEND')

		IF (FRAME .NE. 1) FRAME = 1

		Mu   = M1*M2/(M1+M2)
		Mu   = Mu/TOAMU
		N_T  = 100
		T_in = 0.0D0
  	T_fn = 45.0D0
		dT   = (T_fn-T_in)/REAL(N_T-1)

		ALLOCATE(X(N_T),Y(N_T))

		DO i=1,N_t
			T     = T_in + (i-1)*dT
			CALL angle_to_lab(T*PI/180.0D0, T_lab)
			T_lab = T_lab*180.0D0/PI
			X(i)  = E*TOEV*T_lab/Mu
			CALL pw_amp(T*PI/180.0D0)
			Y(i)  = T_lab*SIN(T_lab*PI/180.0D0)*Amp	
			WRITE(*,*) 'Tcm: Tlab: Amp: ', T, T_lab, Amp
		END DO

		DO i=1,N_t
			WRITE(50,*) X(i), Y(i)
		END DO
		DEALLOCATE(X,Y)
		CLOSE(50)

	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE reduced_coordinates

