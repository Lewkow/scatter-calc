
MODULE grid

	INTEGER,PARAMETER												:: MaxL = 100000

	INTEGER																	:: NGrid
	
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: R
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: F 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: EF 

	REAL(KIND=8)														:: RStart
	REAL(KIND=8)														:: REnd
	REAL(KIND=8)														:: Step 

END MODULE grid
