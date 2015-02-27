
SUBROUTINE diffusion_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the diffusion cross section using
!	numerical integration of 
! Qd = 2*pi*int_{0}^{pi} Amp(theta)*(1-cos(theta))*sin(theta)*dtheta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE mpi_info
  USE dcs, ONLY : Amp
  USE collision_inputs, ONLY : ThI, ThF, NTh
  USE physics_constants, ONLY : PI, TOEV
	USE current_energy, ONLY : E

	IMPLICIT NONE

  INCLUDE 'mpif.h'

  !-----------------------------------
  ! Internal  
  !-----------------------------------

  INTEGER                           :: NT
  INTEGER                           :: i

  REAL(KIND=8)                      :: DT
  REAL(KIND=8)                      :: T
  REAL(KIND=8)                      :: Dum_T
  REAL(KIND=8)                      :: TI
  REAL(KIND=8)                      :: TF
  REAL(KIND=8)                      :: TC
  REAL(KIND=8)                      :: tstart
  REAL(KIND=8)                      :: tend

  NT = 50000
  TI = 0.0
  TF = PI
  DT = (TF-TI)/DBLE(NT)
  TC = 0.0

  IF (rank == 0) THEN
	
		OPEN(UNIT=40, FILE="../Data/DFCS.dat", ACCESS="APPEND")
    tstart = MPI_WTIME()

    DO i=1,(NT+1)
      T = (i-1)*DT
      CALL pw_amp(T)
      TC = TC + 2.0D0*PI*DT*Amp*(1.0D0-COS(T))*SIN(T)
    END DO

    tend = MPI_WTIME()
    WRITE(*,'(A,ES9.2, F9.4,A)') "Integrated DFCS (a0^2): ", TC, (tend-tstart), "(sec)"
		WRITE(40,*) E*TOEV, TC
		CLOSE(40)
  END IF

  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )



END SUBROUTINE diffusion_cross_section

