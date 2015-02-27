
PROGRAM collide
  USE mpi_info
  USE collision_inputs
	USE current_energy
	USE physics_constants, ONLY : TODEG, TOEV, PI
	
	IMPLICIT NONE
	
	INCLUDE 'mpif.h'

	REAL(KIND=8)	:: EStart
	REAL(KIND=8)	:: EEnd
	REAL(KIND=8)	:: DE 
	REAL(KIND=8)	:: Start_t 
	REAL(KIND=8)	:: End_t 
	REAL(KIND=8)	:: Full_Start 
	REAL(KIND=8)	:: Full_End 
	REAL(KIND=8)	:: Full_min
	REAL(KIND=8)	:: Full_sec
	REAL(KIND=8)	:: Elapsed

  INTEGER				:: rc
  INTEGER				:: i, L, j 

	PI    = 4.0D0*ATAN(1.0D0)
	TODEG = 180.0D0/PI	

	CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numrank, ierr )
  CALL collision_input_bcast

	!! Write potential to be used
	IF (rank == 0) THEN
		CALL pot_write
	END IF

	IF ( POT_TEST .EQ. 1 .AND. rank .EQ. 0 ) THEN
		Channel = 'T'
		CALL pot_tester	
	END IF

	!! Start Energy Loop	
	IF (rank == 0) THEN	
		Full_Start = MPI_WTIME()
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	
	DO i=1,NumEngy
		E 			= Engy(i)
		K 			= DSQRT(2.0D0*MU*E)
		CurEngy = i

		IF (PROJ .EQ. 'H' .AND. TARG .EQ. 'H') THEN
			!! First channel to compute, singlet
			Channel = 'S'
		END IF

		IF ( (PROJ .EQ. 'He4' .AND. TARG .EQ. 'O') .OR. (PROJ .EQ. 'O' .AND. TARG .EQ. 'He4') ) THEN
			!! First channel to compute, singlet (PI) 
			Channel = 'S'
		END IF

		IF (rank == 0) THEN
			WRITE(*,*)	
			WRITE(*,*) "        !~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!"
			WRITE(*,*) "                 Energy", i
			WRITE(*,*) "        !~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!"
			WRITE(*,*)	
			Start_t = MPI_WTIME()
		END IF

		!! PW Phaseshifts
		IF ( (CALC_TYPE == 0) .OR. (CALC_TYPE == 99) ) THEN

			CALL get_phaseshifts	
			CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

			!! If H+H, get singlet and triplet phaseshifts
			!! and save to s_phase and t_phase
			IF (PROJ .EQ. 'H' .AND. TARG .EQ. 'H') THEN

				S_LMax = LMax
				ALLOCATE( S_phase(S_LMax) )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DO j=1,S_LMax
					S_phase(j) = phase(j)
				END DO
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DEALLOCATE( phase )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				Channel = 'T'
				CALL depth_finder
				CALL get_phaseshifts
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				T_LMax = LMax
				ALLOCATE( T_phase(T_LMax) )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DO j=1,T_LMax
					T_phase(j) = phase(j)
				END DO

			END IF ! H+H	

			!! If He4+O get singlet and triplet phaseshifts
			!! and save to s_phase and t_phase			
			IF (PROJ .EQ. 'He4' .AND. TARG .EQ. 'O' .OR. PROJ .EQ. 'O' .AND. TARG .EQ. 'He4') THEN
				
				S_LMax = LMax
				ALLOCATE( S_phase(S_LMax) )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DO j=1,S_LMax
					S_phase(j) = phase(j)
				END DO
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DEALLOCATE( phase )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				Channel = 'T'
				CALL depth_finder
				CALL get_phaseshifts
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				T_LMax = LMax
				ALLOCATE( T_phase(T_LMax) )
				CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
				DO j=1,T_LMax
					T_phase(j) = phase(j)
				END DO

			END IF ! He4+O
	
		END IF ! PW phaseshift

		IF (CALC_TYPE == 1) THEN
			CALL RStart_finder
			CALL REnd_finder
		END IF

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

		!! Call tests specified by input file
		CALL call_tests

		!! Barrier before deallocating phase array
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

		IF ( CALC_TYPE .EQ. 0 .OR. CALC_TYPE .EQ. 99 ) THEN
			DEALLOCATE(phase)
			IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN
				DEALLOCATE( S_phase )
				DEALLOCATE( T_phase )
			ELSE IF (Proj .EQ. 'He4' .AND. Targ .EQ. 'O' .OR. Proj .EQ. 'O' .AND. Targ .EQ. 'He4') THEN
				DEALLOCATE( S_phase )
				DEALLOCATE( T_phase )
			END IF
		END IF

		IF (rank == 0) THEN
			End_t = MPI_WTIME()
			WRITE(*,*)
			IF ( (End_t-Start_t) .GE. 60.0 ) THEN
				Full_Min = (End_t-Start_t)/60
				Full_Sec = MOD((End_t-Start_t),60.0)
				WRITE(*,"(A,I5,A,I2,A)") "Energy complete in ", INT(Full_Min), " min", &
				& INT(Full_Sec), " sec"
			ELSE	
				WRITE(*,"(A,I2,A)") "Energy complete in ", INT(End_t-Start_t), " sec"
			END IF	
			WRITE(*,*)
			!! write Lend to file
			CALL phase_end_writer
		END IF

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		
	END DO ! i 
	
	IF (rank == 0) THEN
		Full_End = MPI_WTIME()
		Elapsed  = Full_End-Full_Start
		Full_Min = Elapsed/60
		Full_Sec = MOD(Elapsed,60.0) 
		WRITE(*,*)
		WRITE(*,FMT="(A,I5,A,I2,A)") "ASTROSCATT complete in ",INT(Full_Min), " min ", INT(Full_Sec)," sec"
		WRITE(*,FMT="(A,I3,A,I5,A,I2,A)") "Time saved by using", numrank, " processors: ", &
		& INT((Elapsed*(numrank-1))/60)," min ", INT(MOD((Elapsed*(numrank-1)),60.0)), " sec"
		WRITE(*,*)
	END IF

	DEALLOCATE(Engy)

	CALL MPI_FINALIZE(rc)

END PROGRAM collide 
