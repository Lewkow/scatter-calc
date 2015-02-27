
SUBROUTINE eik_total_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate the total cross section for a given
! collision energy using the Eikonal approximation. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE current_energy, ONLY : K, E
	USE physics_constants, ONLY : TOEV, PI
	USE eik, ONLY : AmpInt, Eik_RStart
	USE grid, ONLY : REnd
	USE collision_inputs, ONLY : IDNT 

	IMPLICIT NONE

	INCLUDE 'mpif.h'
	
	!------------------------------
	! Internal
	!------------------------------

	REAL(KIND=8)										:: TotCross ! Total Cross Section 
	REAL(KIND=8)										:: J        ! Value of bessj0 
	REAL(KIND=8)										:: P        ! Phase
	REAL(KIND=8)										:: b        ! Impact parameter 
	REAL(KIND=8)										:: db       ! Differential impact parameter 
	REAL(KIND=8)										:: chunk    ! Chunk of impact param for rank 
	REAL(KIND=8)										:: tstart   ! Start timer 
	REAL(KIND=8)										:: tend     ! End timer 
	REAL(KIND=8),DIMENSION(numrank)	:: buff			! Buffer used to gather TCS
	
	INTEGER													:: i				! Counter
	INTEGER													:: my_start	! Ranks starting b value 
	INTEGER													:: my_end		! Ranks ending b value 

	chunk    = AmpInt/numrank
	my_start = rank*chunk
	my_end   = (rank+1)*chunk - 1
	db       = (REnd - Eik_RStart)/DBLE(AmpInt) 
	TotCross = 0.0D0	
	tstart   = MPI_WTIME()
	
	IF (rank == (numrank - 1)) THEN
		my_end = AmpInt
	END IF

	DO i=my_start,my_end
		b = Eik_RStart + i*db
		CALL eik_amp_coeff(0.0D0,b,J)
		CALL eik_phase(b,P)
		TotCross = TotCross + b*J*SIN(P)*SIN(P)
	END DO

	IF ( (IDNT == 1) .OR. (IDNT == 2) ) THEN
		TotCross = TotCross*16.0*PI*db
	ELSE IF (IDNT == 0) THEN
		TotCross = TotCross*8.0*PI*db
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )	
	CALL MPI_GATHER( TotCross, 1, MPI_REAL8, &
	& buff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

	IF (rank == 0) THEN
		tend     = MPI_WTIME()
		TotCross = 0.0D0
		DO i=1,numrank
			TotCross = TotCross + buff(i)
		END DO
		OPEN(UNIT=50, FILE="../Data/Eik_TCS.dat", ACCESS="APPEND")
		WRITE(50,*) E*TOEV, TotCross
		WRITE(*,*) "Eikonal OP TCS (a0^2): ", TotCross, (tend-tstart), "(sec)"
		CLOSE(50)
	END IF
	
END SUBROUTINE eik_total_cross_section
