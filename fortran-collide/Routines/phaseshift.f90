
SUBROUTINE get_phaseshifts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get phaseshifts for a given set of projectile, 
! target and collision energy. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE collision_inputs, ONLY : CurEngy, NumEngy, PHASE_WRITE, PHASE_TIMING, Proj, Targ
	USE current_energy, ONLY : E, K, LMax, phase, Channel
	USE physics_constants, ONLY : TOEV
	USE grid, ONLY : NGrid, MaxL, F, R, EF

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!-------------------------------
	!	Internal	
	!-------------------------------

	REAL(KIND=8)													:: current_phase
	REAL(KIND=8)													:: small_phase
	REAL(KIND=8)													:: start_full
	REAL(KIND=8)													:: end_full 

	REAL(KIND=8),DIMENSION(MaxL/numrank)	:: phase_array

	INTEGER																:: myL
	INTEGER																:: my_count
	INTEGER																:: max_count 
	INTEGER																:: my_max_count 
	INTEGER																:: my_LMax
	INTEGER																:: xx 
	INTEGER																:: L 
	INTEGER																:: i 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	xx					   = 0
	small_phase    = 5.0e-7
	max_count      = 30
	phase_array(:) = 0.0
	
	IF (numrank>max_count) THEN
		my_max_count = 1
	ELSE
		my_max_count = max_count/numrank  
	END IF
	
	IF (rank == 0) THEN
		CALL RStart_finder
		CALL REnd_finder
		CALL step_finder
	END IF ! rank=0

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL phase_bcast

	!-----------------------------
	! Allocate R, F, RF	
	!-----------------------------

	ALLOCATE(R(NGrid))	
	ALLOCATE(F(NGrid))	
	ALLOCATE(EF(NGrid))	

	CALL grid_build
	
	start_full = MPI_WTIME()
	my_count   = 0

	DO L=0,MaxL/numrank
		myL	= rank+L*numrank
		CALL numerov(myL,current_phase)
		phase_array(L+1) = current_phase
		
		IF (ABS(current_phase) <= small_phase) THEN
			my_count = my_count + 1

			IF (my_count .GE. my_max_count) THEN
				my_LMax = myL
				xx = 1	
			END IF	

		END IF

		IF (xx == 1) EXIT
		
	END DO ! L
	
	IF (numrank > 1) THEN
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		CALL my_LMax_collect(L)
		CALL LMax_collect(my_LMax)
		CALL phase_reduce(phase_array)	
	ELSE
		ALLOCATE(phase(my_LMax))
		DO i=1,my_LMax	
			phase(i)= phase_array(i)
		END DO ! i
		LMax  		= L - 1
	END IF ! numrank > 1

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	end_full = MPI_WTIME()
	
	IF (rank == 0) THEN

		WRITE(*,*) "--------------------------------------------------"
		IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN
			IF (Channel .EQ. 'S') THEN
				WRITE(*,*) "%        Singlet Phases Complete                 %"	
			ELSE IF (Channel .EQ. 'T') THEN		
				WRITE(*,*) "%        Triplet Phases Complete                 %"	
			END IF
		ELSE
			WRITE(*,*) "%                Phases Complete                 %"	
		END IF
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,"(A, ES8.2, A)") "Energy             ", E*TOEV, " (eV)"
		IF ( (end_full-start_full) .LT. 60.0D0 ) THEN
			WRITE(*,"(A, F6.2, A)") "Phases Complete in ", end_full-start_full, " sec"
		ELSE
			WRITE(*,"(A,I4,A,F6.2,A)") "Phases Complete in ", &
			& INT((end_full-start_full)/60.0D0), " min ", MOD((end_full-start_full),60.0D0), " sec"
		END IF
		WRITE(*,*) LMax, " Phases Computed"	
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,*)

		IF (PHASE_TIMING == 1) THEN
			CALL phase_timing_test(end_full-start_full)
		END IF

		IF (PHASE_WRITE == 1) THEN
			CALL phase_writer
			CALL eik_phase_writer
		END IF

	END IF 

	!---------------------------------
	! Deallocate memory for EF, F, R	
	!---------------------------------

	DEALLOCATE(EF)
	DEALLOCATE(F)
	DEALLOCATE(R)

END SUBROUTINE get_phaseshifts

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE numerov(L,ph)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate phaseshift for a given L value using
! Numerov's method. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE grid, ONLY : NGrid, EF, F, R, Step
	USE current_energy, ONLY : K

	IMPLICIT NONE

	!---------------------------
	! Inputs
	!---------------------------

	INTEGER						:: L

	!---------------------------
	! Outputs
	!---------------------------

	REAL(KIND=8)			:: ph

	!---------------------------
	! Internal  
	!---------------------------

	REAL(KIND=8)											:: Ksq
	REAL(KIND=8)											:: Lsq
	REAL(KIND=8)											:: HH
	REAL(KIND=8)											:: D1 
	REAL(KIND=8)											:: D2 
	REAL(KIND=8)											:: D4 
	REAL(KIND=8)											:: D12 
	REAL(KIND=8)											:: D24 
	REAL(KIND=8)											:: Grad
	REAL(KIND=8)											:: KR 
	REAL(KIND=8)											:: sine 
	REAL(KIND=8)											:: cose 
	REAL(KIND=8)											:: dsine 
	REAL(KIND=8)											:: dcose 
	REAL(KIND=8)											:: S 
	REAL(KIND=8)											:: C

	REAL(KIND=8),DIMENSION(NGrid+2)		:: Wave 
	
	INTEGER														:: i
	INTEGER														:: j 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Ksq = K*K
	Lsq = DBLE(L*(L+1))
	HH  = (Step*Step)/12.0

	DO i=1,NGrid
		EF(i) = F(i) + Ksq - Lsq/(R(i)*R(i))
	END DO ! i

	Wave(1) = 0.0D0
	Wave(2) = 0.1D0

	DO i=2,(NGrid-1)

		IF ( ABS(Wave(i)) .GT. 1.0e20 ) THEN
			DO j=1,i
				Wave(j) = Wave(j)*1.0e-20
			END DO
		END IF

		IF ( ABS(Wave(i)) .LT. 1.0e-20 ) THEN
			DO j=1,i
				Wave(j) = Wave(j)*1.0e20
			END DO
		END IF	

		Wave(i+1) = ((2.0-10.0*HH*EF(i))*Wave(i)-(1.0+HH*EF(i-1)) &
		& *Wave(i-1))/(1.0+HH*EF(i+1))
	
	END DO !i

	D1    = ( Wave(NGrid-3) - Wave(NGrid-5) )/(2.0*Step)
	D2    = ( Wave(NGrid-2) - Wave(NGrid-6) )/(4.0*Step)
	D4    = ( Wave(NGrid)   - Wave(NGrid-8) )/(8.0*Step)
	D12   = D1 + (D1-D2)/3.0
	D24   = D2 + (D2-D4)/3.0	
	Grad  = D12 + (D12-D24)/15.0
	KR    = K*R(NGrid-4)
	
	CALL AJ(KR,L,sine)
	CALL AN(KR,L,cose)
	CALL AJP(KR,L,dsine)
	CALL ANP(KR,L,dcose)

	dsine = dsine*K
	dcose = dcose*K

	S     = (dsine*Wave(NGrid-4) - sine*Grad)/K
	C     = (dcose*Wave(NGrid-4) - cose*Grad)/K
	ph    = ATAN(S/C)

END SUBROUTINE numerov

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE phase_bcast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Collect grid variables, Step, RStart and REnd 
! from root and broadcast to all other ranks using
! a buffer. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE grid, ONLY : Step, RStart, REnd, NGrid

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	REAL(KIND=8), DIMENSION(3)	:: buff_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF (rank == 0) THEN
		buff_real(1) = Step
		buff_real(2) = RStart
		buff_real(3) = REnd
	END IF

	CALL MPI_BCAST( buff_real, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NGrid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )		

	IF (rank .NE. 0) THEN
		Step   = buff_real(1)
		RStart = buff_real(2)
		REnd   = buff_real(3)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr ) 

END SUBROUTINE phase_bcast

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE phase_reduce(my_phases)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gather individual phase arrays, my_phases,  
! from all ranks and build contiguous phase array, 
! phase, which is broadcast to all ranks. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE current_energy, ONLY : phase, LMax, my_LMax	

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!---------------------------
	! Inputs
	!---------------------------

	REAL(KIND=8),DIMENSION(my_LMax)		:: my_phases

	!---------------------------
	! Internal
	!---------------------------

	REAL(KIND=8),DIMENSION(LMax)			:: buff 

	INTEGER														:: i
	INTEGER														:: j 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ALLOCATE(phase(LMax+1))

	DO i=1,(LMax+1)
		phase(i) = 0.0
	END DO

	CALL MPI_GATHER( my_phases, my_LMax, MPI_REAL8, &
	& buff, my_LMax, MPI_REAL8, 0, MPI_COMM_WORLD, ierr ) 	

	IF (rank == 0) THEN
		DO i=0,(LMax-numrank)
			IF ( MOD(i,numrank) == 0 ) THEN
				j = i/numrank
			END IF
			phase(i+1) = buff( MOD(i,numrank)*my_LMax + j + 1)
		END DO
	END IF	

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( phase, (LMax+1), MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

END SUBROUTINE phase_reduce

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE LMax_collect(L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Collect all ranks maximum L value, LMax, and 
! broadcast that value to all ranks.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE current_energy, ONLY : LMax

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!---------------------------
	! Inputs
	!---------------------------

	INTEGER											:: L

	!---------------------------
	! Internal
	!---------------------------

	INTEGER											:: i

	INTEGER,DIMENSION(numrank)	:: buff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LMax = 0

	CALL MPI_GATHER( L, 1, MPI_INTEGER, buff, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	
	IF (rank == 0) THEN
		DO i=1,numrank
			IF (buff(i) .GT. LMax) THEN
				LMax = buff(i)
			END IF
		END DO
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )	
	CALL MPI_BCAST( LMax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )	

END SUBROUTINE LMax_collect

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE my_LMax_collect(L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Collect all ranks convergent L value
! my_LMax, the number of phases each rank has
! computed, and broadcast the highest convergent 
! L value to the rest of the ranks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE current_energy, ONLY : my_LMax

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!---------------------------
	! Inputs
	!---------------------------

	INTEGER											:: L

	!---------------------------
	! Internal
	!---------------------------

	INTEGER											:: i

	INTEGER,DIMENSION(numrank)	:: buff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	my_LMax = 0

	CALL MPI_GATHER( L, 1, MPI_INTEGER, buff, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	
	IF (rank == 0) THEN
		DO i=1,numrank
			IF (buff(i) .GT. my_LMax) THEN
				my_LMax = buff(i)
			END IF
		END DO
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )	
	CALL MPI_BCAST( my_LMax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )	

END SUBROUTINE my_LMax_collect

