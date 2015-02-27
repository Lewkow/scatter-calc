
SUBROUTINE hard_sphere_test
	USE current_energy, ONLY : K, E
	USE mpi_info, ONLY : rank
	USE physics_constants, ONLY : TOEV
	INCLUDE 'mpif.h'

	!-----------------------------
	!	Internal	
	!-----------------------------

	REAL(KIND=8)									:: PartWave
	REAL(KIND=8)									:: Theory
	REAL(KIND=8)									:: Error 
	
	IF (rank == 0) THEN	
		CALL total_cross_section(PartWave)
		WRITE(*,*) "PW TC", PartWave
		CALL hard_sphere_TCS(Theory)
		WRITE(*,*) "HS TC", Theory 
		Error = ABS(PartWave - Theory)*100.0/ABS(Theory)
		WRITE(*,*) "E ", E*TOEV, "PW ", PartWave, "TH ", Theory, "Err ", Error
	END IF
	
END SUBROUTINE hard_sphere_test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE hard_sphere_TCS(TC)
	USE current_energy, ONLY : K, E
	USE physics_constants, ONLY : PI	
	
	!-----------------------------
	!	Inputs	
	!-----------------------------

	REAL(KIND=8)								:: TC	
	
	!-----------------------------
	!	Internal	
	!-----------------------------

	INTEGER,PARAMETER						:: MX = 10000
	
	REAL(KIND=8)								:: J	
	REAL(KIND=8)								:: N	
	REAL(KIND=8)								:: small	
	REAL(KIND=8)								:: C	

	INTEGER											:: i
	INTEGER											:: L

	small = 1.0e-10
	i     = 0
	TC    = 0

	DO L=0,(MX-1)
		CALL AJ(K,L,J)
		CALL AN(K,L,N)
		J = J/K
		N = N/K
		C	= (2.0*L+1.0)*J*J/(J*J+N*N)
		TC = TC + C
		IF (C .LT. small) THEN
			WRITE(*,*) "Hard Sphere has converged with L ", L
			EXIT
		END IF
	END DO

	TC = TC*4.0*PI/(K*K)

END SUBROUTINE hard_sphere_TCS


