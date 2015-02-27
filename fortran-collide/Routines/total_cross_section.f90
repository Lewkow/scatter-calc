
SUBROUTINE write_total_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section using partial 
! wave phases with the optical theorem.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE collision_inputs, ONLY : IDNT, Proj, Targ
	USE current_energy, ONLY : LMax, K, E, phase
	USE physics_constants, ONLY : PI, TOEV	

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!-------------------------
	!	Internal	
	!-------------------------

	REAL(KIND=8)											:: Cross
	REAL(KIND=8)											:: sphas
	REAL(KIND=8)											:: tstart
	REAL(KIND=8)											:: tend 
	REAL(KIND=8)											:: TC 

	REAL(KIND=8),DIMENSION(numrank)		:: TC_buff 
	
	INTEGER														:: L
	INTEGER														:: chunk 
	INTEGER														:: my_st 
	INTEGER														:: my_fn 

	!!!!!!!!!!!!!!
	!! H+H 
	!!!!!!!!!!!!!!
	IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN	
		IF (rank .EQ. 0) THEN
			tstart = MPI_WTIME()
			CALL HH_TCS(TC)
		END IF
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!!!!!!!!!!!!!!
	!! He + O
	!!!!!!!!!!!!!!
  ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN
		IF (rank .EQ. 0) THEN
			tstart = MPI_WTIME()
	    CALL HeO_TCS(TC)
		END IF
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	ELSE
		chunk = LMax/numrank	
		my_st = rank*chunk
		my_fn = (rank+1)*chunk - 1
		Cross = 0

		IF (rank == 0) THEN
			tstart = MPI_WTIME()
		END IF

		IF (rank == (numrank-1)) THEN
			my_fn = LMax
		END IF

		IF (IDNT == 0) THEN
			!------------------------
			! Non-Ident Particles	
			!------------------------
			DO L=my_st,my_fn
				sphas = DSIN(phase(L+1))
				Cross = Cross+(2.0D0*L+1.0D0)*sphas*sphas
			END DO
		
		ELSE IF (IDNT == 1) THEN
			!------------------------
			! Fermion Particles
			!------------------------
			DO L=my_st,my_fn
				sphas = SIN(phase(L+1))
				Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0-(-1.0)**L)
			END DO
	
		ELSE
			!------------------------
			! Boson Particles
			!------------------------
			DO L=my_st,my_fn
				sphas = SIN(phase(L+1))
				Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0+(-1.0)**L)
			END DO

		END IF

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		CALL MPI_GATHER( Cross, 1, MPI_REAL8, TC_buff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )	

	 	IF ( rank .EQ. 0) THEN
			TC = 0

			DO L=1,numrank
				TC = TC + TC_buff(L)
			END DO

			IF (IDNT == 0) THEN
				TC = TC*4.0*PI/(K*K)
			ELSE
				TC = TC*8.0*PI/(K*K)
			END IF

		END IF		

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

	IF (rank == 0) THEN

		tend = MPI_WTIME()

		OPEN(UNIT=40, FILE="../Data/PW_TCS.dat", ACCESS="APPEND")
		WRITE(40,*) E*TOEV, TC
		WRITE(*,"(A,ES9.2, F9.4,A)") "Optical    TCS (a0^2): ", TC, (tend-tstart), " (sec)"
		CLOSE(40)	

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

END SUBROUTINE write_total_cross_section

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE total_cross_section(TC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section using partial 
! wave phases with the optical theorem.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : IDNT, Proj, Targ
	USE current_energy, ONLY : LMax, K, E, phase
	USE physics_constants, ONLY : PI, TOEV	

	IMPLICIT NONE

	!-------------------------
	!	Outputs	
	!-------------------------

	REAL(KIND=8)											:: TC 

	!-------------------------
	!	Internal	
	!-------------------------

	REAL(KIND=8)											:: Cross
	REAL(KIND=8)											:: sphas
	REAL(KIND=8)											:: tstart
	REAL(KIND=8)											:: tend 

	INTEGER														:: L

	IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN
		CALL HH_TCS(TC)
	ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN
		CALL HeO_TCS(TC)
	ELSE

		Cross = 0.0

		CALL cpu_time(tstart)
		IF (IDNT == 0) THEN
			!------------------------
			! Non-Ident Particles	
			!------------------------
			DO L=0,LMax
				sphas = SIN(phase(L+1))
				Cross = Cross+(2.0*L+1.0)*sphas*sphas
			END DO
		ELSE IF (IDNT == 1) THEN
			!------------------------
			! Fermion Particles
			!------------------------
			DO L=0,LMax
				sphas = SIN(phase(L+1))
				Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0-(-1.0)**L)
			END DO
		ELSE
			!------------------------
			! Boson Particles
			!------------------------
			DO L=0,LMax
				sphas = SIN(phase(L+1))
				Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0+(-1.0)**L)
			END DO
		END IF
	
		CALL cpu_time(tend) 
		IF (IDNT == 0) THEN
			TC = Cross*4.0*PI/(K*K)
		ELSE
			TC = Cross*8.0*PI/(K*K)
		END IF

	END IF
	
END SUBROUTINE total_cross_section

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE total_cross_section_MPI(TC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section using partial 
! wave phases with the optical theorem.  
!						
!									MPI VERSION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE collision_inputs, ONLY : IDNT	
	USE current_energy, ONLY : LMax, K, E, phase
	USE physics_constants, ONLY : PI, TOEV	
	INCLUDE 'mpif.h'

	!-------------------------
	!	Outputs	
	!-------------------------

	REAL(KIND=8)											:: TC 

	!-------------------------
	!	Internal	
	!-------------------------

	REAL(KIND=8)											:: Cross
	REAL(KIND=8)											:: sphas
	REAL(KIND=8)											:: tstart
	REAL(KIND=8)											:: tend 

	REAL(KIND=8),DIMENSION(numrank)		:: TC_buff 
	
	INTEGER														:: L
	INTEGER														:: chunk 
	INTEGER														:: my_st 
	INTEGER														:: my_fn 

	chunk = LMax/numrank	
	my_st = rank*chunk
	my_fn = (rank+1)*chunk - 1
	Cross = 0

	IF (rank == 0) THEN
		tstart = MPI_WTIME()
	END IF

	IF (rank == (numrank-1)) THEN
		my_fn = LMax
	END IF

	IF (IDNT == 0) THEN
		!------------------------
		! Non-Ident Particles	
		!------------------------
		DO L=my_st,my_fn
			sphas = SIN(phase(L+1))
			Cross = Cross+(2.0*L+1.0)*sphas*sphas
		END DO
	
	ELSE IF (IDNT == 1) THEN
		!------------------------
		! Fermion Particles
		!------------------------
		DO L=my_st,my_fn
			sphas = SIN(phase(L+1))
			Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0-(-1.0)**L)
		END DO

	ELSE
		!------------------------
		! Boson Particles
		!------------------------
		DO L=my_st,my_fn
			sphas = SIN(phase(L+1))
			Cross = Cross+(2.0*L+1.0)*sphas*sphas*(1.0+(-1.0)**L)
		END DO

	END IF

	CALL MPI_GATHER( Cross, 1, MPI_REAL8, TC_buff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )	

	IF (rank == 0) THEN
		tend = MPI_WTIME()
		TC   = 0
		DO L=1,numrank
			TC = TC + TC_buff(L)
		END DO
		IF (IDNT == 0) THEN
			TC = TC*4.0*PI/(K*K)
		ELSE
			TC = TC*8.0*PI/(K*K)
		END IF
	END IF

	CALL MPI_BCAST( TC, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE total_cross_section_MPI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE HeO_TCS(TCS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section for He+O collision
! which requires multiple channel calculations, 
! and thus multiple phases, one for each potential. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : S_phase, T_phase, S_LMax, T_LMax, K
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!!! Output
	REAL(KIND=8) 	:: TCS

	!!! Internal
	REAL(KIND=8)	:: T_Amp, S_Amp, S
	INTEGER				:: L

	S_Amp = 0.0D0
	T_Amp = 0.0D0

	!!! Singlet (PI) channel
	DO L=0,S_LMax
		S     = SIN(S_phase(L+1))
		S_Amp = S_Amp + (2.0D0*L+1.0D0)*S*S
	END DO

	!!! Triplet (SIGMA) channel
	DO L=0,T_LMax
		S     = SIN(T_phase(L+1))
		T_Amp = T_Amp + (2.0D0*L+1.0D0)*S*S
	END DO

	TCS = ((4.0D0*PI)/(K*K))*( (2.0D0/3.0D0)*S_Amp + (1.0D0/3.0D0)*T_Amp )	

END SUBROUTINE HeO_TCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE HH_TCS(TCS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section for H+H collision
! which requires multiple channel calculations, 
! and thus multiple phases, one for each potential. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : S_phase, T_phase, S_LMax, T_LMax, K
	USE physics_constants, ONLY : PI
	USE collision_inputs, ONLY : IDNT

	IMPLICIT NONE
	
	!-----------------------------------
	! Output
	!-----------------------------------
	
	REAL(KIND=8)		:: TCS	

	!-----------------------------------
	! Internal 
	!-----------------------------------

	INTEGER					:: L
	REAL(KIND=8)		:: TAmp
	REAL(KIND=8)		:: SAmp
	REAL(KIND=8)		:: S
	REAL(KIND=8)		:: C

	SAmp = 0.0D0
	TAmp = 0.0D0

	IF (IDNT .NE. 0) THEN

		DO L=0,(S_LMax-1)
			S    = SIN(S_phase(L+1))
			IF (MOD(L,2) .EQ. 0) THEN
				C = 1.0D0/8.0D0
			ELSE 
				C = 3.0D0/8.0D0
			END IF
			SAmp = SAmp + C*S*S*(2.0D0*L+1.0D0)
		END DO	

		DO L=0,(T_LMax-1)
			S    = SIN(T_phase(L+1))
			IF (MOD(L,2) .EQ. 0) THEN
				C = 9.0D0/8.0D0
			ELSE
				C = 3.0D0/8.0D0
			END IF
			TAmp = TAmp + C*S*S*(2.0D0*L+1.0D0)
		END DO	

		TCS = (4.0D0*PI/(K*K))*(TAmp + SAmp)

	ELSE
		
		DO L=0,(S_LMax-1)
			S    = SIN(S_phase(L+1))
			SAmp = SAmp + S*S*(2.0D0*L+1.0D0)
		END DO

		DO L=0,(T_LMax-1)
			S    = SIN(T_phase(L+1))
			TAmp = TAmp + S*S*(2.0D0*L+1.0D0)
		END DO

		TCS = (4.0D0*PI/(K*K))*(0.25D0*SAmp + 0.75D0*TAmp)

	END IF		

END SUBROUTINE HH_TCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE num_total_cross_section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the total cross section using phase
! shifts from partial wave analysis. TCS is 
! computed using a numerical integration from 
! theta=0,PI
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE dcs, ONLY : Amp
	USE collision_inputs, ONLY : ThI, ThF, NTh
	USE physics_constants, ONLY : PI, TOEV
	USE current_energy, ONLY : E
	INCLUDE 'mpif.h'

	!-----------------------------------
	! Internal	
	!-----------------------------------

	INTEGER														:: NT
	INTEGER														:: i 

	REAL(KIND=8)											:: DT	
	REAL(KIND=8)											:: T	
	REAL(KIND=8)											:: Dum_T	
	REAL(KIND=8)											:: TI	
	REAL(KIND=8)											:: TF	
	REAL(KIND=8)											:: TC	
	REAL(KIND=8)											:: tstart	
	REAL(KIND=8)											:: tend	
	REAL(KIND=8)											:: y
	REAL(KIND=8)											:: c

	!!! NT must be even to apply Simpson's Rule
	NT = 100000
	TI = 0.0D0
	TF = PI
	DT = TF/DBLE(NT)
	TC = 0.0D0	

	IF (rank == 0) THEN
		tstart = MPI_WTIME()

		DO i=0,NT
			T = i*DT
			CALL pw_amp(T)
			y  = 2.0D0*PI*Amp*DSIN(T)

			IF ( i .EQ. 0 .OR. i .EQ. NT ) THEN
				c = 1.0D0
			ELSE IF ( MOD(i,2) .EQ. 0 ) THEN
				c = 2.0D0
			ELSE 
				c = 4.0D0
			END IF

			TC = TC + y*c 
		END DO

		TC = TC*DT/3.0D0

		tend = MPI_WTIME()
		WRITE(*,"(A,ES10.2,A,F6.2,A)") "Integrated TCS (a0^2): ", TC, " ", (tend-tstart), "(sec)"	
		OPEN(UNIT=600, FILE="../Data/NUM_PW_TCS.dat", ACCESS="APPEND")
		WRITE(600,*) E*TOEV, TC
		CLOSE(600)
	END IF ! rank=0

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE num_total_cross_section

