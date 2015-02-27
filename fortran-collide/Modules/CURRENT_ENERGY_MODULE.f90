
MODULE current_energy

	REAL(KIND=8)															:: K
	REAL(KIND=8)															:: E

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)			:: leg
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)			:: phase
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)			:: S_phase		! Singlet phase
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)			:: T_phase		! Triplet phase
	
	CHARACTER(LEN=1)													:: Channel	

	INTEGER																		:: S_LMax			! Singlet L_Max
	INTEGER																		:: T_LMax			! Triplet L_Max
	INTEGER																		:: LMax
	INTEGER																		:: my_LMax

END MODULE current_energy

