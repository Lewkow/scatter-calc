
SUBROUTINE RStart_finder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find the start of the grid depending on the energy
! of collision. The starting point is the distance
! which is 30% higher in energy than the classical
! turning point. This 30% should account for 
! any quantum tunneling. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE grid, ONLY : RStart
	USE current_energy, ONLY : E

	IMPLICIT NONE

	REAL(KIND=8)									:: r
	REAL(KIND=8)									:: start
	REAL(KIND=8)									:: pot
	REAL(KIND=8)									:: lim
	REAL(KIND=8)									:: dr

	INTEGER												:: j
	INTEGER												:: i

	i 		= 0
	j 		= 0
	start = 10.0
	lim		= 0.3
	dr		= 0.0001
	lim 	= E*(1.0+lim)
	
	DO i=0,10000000 
		r   = start - j*dr
		CALL potential(r,pot)
		IF ( (pot .GE. lim) .OR. (r .LE. 1.0e-5) ) THEN
			EXIT	
		END IF
		j = j+1
	END DO ! i=0

	RStart = r

END SUBROUTINE  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE REnd_finder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find the ending point on the grid. The end is the
! point which has a change in energy smaller than 
! the value of lim. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE grid, ONLY : REnd
	USE physics_constants, ONLY : TOEV

	IMPLICIT NONE

	REAL(KIND=8)								:: ipot
	REAL(KIND=8)								:: iipot
	REAL(KIND=8)								:: dpot 
	REAL(KIND=8)								:: lim 
	REAL(KIND=8)								:: Ri 
	REAL(KIND=8)								:: DR 
	REAL(KIND=8)								:: R 

	INTEGER											:: i

	lim = 0.00000001
	Ri 	= 1.0
	DR  = 0.1
	CALL potential(Ri,ipot)

	DO i=1,100000
		R = Ri + i*DR
		CALL potential(R,iipot)
		dpot = ABS(iipot - ipot)
		IF (dpot .LT. lim) THEN
			REnd = R
			EXIT
		END IF
		ipot = iipot
	END DO 

END SUBROUTINE REnd_finder

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE step_finder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Finds the steplength to use for the grid and also
! finds the number of points to use in the grid
! based on the energy, well depth and starting/ending
! points of the grid. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE grid, ONLY : NGrid, RStart, REnd, Step
  USE current_energy, ONLY : E
  USE collision_inputs, ONLY : Mu, Depth

	IMPLICIT NONE

  REAL(KIND=8)          :: KMax
  REAL(KIND=8)          :: HMax

  KMax = SQRT(2.0*Mu*(E+Depth))
	HMax = (1.0/KMax)/20.0
  
  IF (7.0e-4 > HMax) THEN
    Step = HMax
  ELSE
    Step = 7.0e-4
  END IF

  NGrid = INT( (REnd - RStart)/Step ) + 1

END SUBROUTINE step_finder

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE grid_build
  USE grid, ONLY : R, F, Step, RStart, REnd, NGrid
  USE collision_inputs, ONLY : Mu
	USE mpi_info, ONLY : rank

	IMPLICIT NONE

	INCLUDE 'mpif.h'

  INTEGER               :: i
  REAL(KIND=8)          :: pot

	55 FORMAT(A,I8)
	56 FORMAT(A,ES8.2)

	IF (rank == 0) THEN
		WRITE(*,*) "--------------------------------------------------"	
		WRITE(*,*) "%               Grid Information                 %"	
		WRITE(*,*) "--------------------------------------------------"	
		WRITE(*,55) "Size          ", NGrid
		WRITE(*,56) "Start    (a0) ", RStart
		WRITE(*,56) "End      (a0) ", REnd
		WRITE(*,56) "Stepsize (a0) ", Step
		WRITE(*,*) "--------------------------------------------------"	
		WRITE(*,*)
	END IF  

	CALL potential(RStart,pot)
  R(1) = RStart
  F(1) = -2.0*Mu*pot

  DO i=2,NGrid
    R(i) = R(i-1) + Step
    CALL potential(R(i),pot)
    F(i) = -2.0*Mu*pot
  END DO 

END SUBROUTINE grid_build

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

