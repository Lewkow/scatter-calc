
PROGRAM TCS_data_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program is made to numerically
! integrate a data file containing 2 
! columns of data (| angle | DCS |). 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	INTEGER																	:: N
	INTEGER																	:: IO
	INTEGER																	:: i
	INTEGER																	:: InputLength

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Ang
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: DCS 

	CHARACTER(LEN=50)												:: UserInput
	CHARACTER,ALLOCATABLE,DIMENSION(:) 			:: Filename

	WRITE(*,*) "	TCS INTEGRATOR"
	WRITE(*,*) "------------------"
	WRITE(*,*)
	WRITE(*,*) "	Please Enter the Filename containing data to integrate"
	READ(*,*) UserInput
!	WRITE(*,*) "	You have entered ", UserInput, "!!"	

	i = LEN(UserInput)
	DO WHILE (UserInput(i:i) .EQ. ' ')
		i = i - 1
	END DO
	
	InputLength = i

!	WRITE(*,*) "Length of filename is", InputLength	

	ALLOCATE(Filename(InputLength))
	
	DO i=1,InputLength
		Filename(i:i) = UserInput(i:i)
	END DO

	WRITE(*,*) "	You have entered ", Filename, "!!"

!	OPEN(UNIT=666, FILE=UserInput, ACCESS="APPEND")
	OPEN(UNIT=666, FILE=Filename, ACCESS="APPEND")

!	N = 0
!	DO
!		READ(666, *, IOSTAT = IO)	
!		IF (IO < 0) EXIT
!		N = N+1
!	END DO

!	WRITE(*,*) "Number of lines is ", N	

	CLOSE(666)	
	DEALLOCATE(Filename)

END PROGRAM

