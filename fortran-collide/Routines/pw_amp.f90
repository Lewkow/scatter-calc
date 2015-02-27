
SUBROUTINE pw_amp_MPI(theta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes scattering amplitude using partial
! wave analysis. Phases must be computed prior
! to computing amplitude. 
!
!								MPI VERSION 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : leg, phase, K, LMax
	USE dcs, ONLY : Amp
	USE collision_inputs, ONLY : FRAME, Proj, Targ	
	USE mpi_info
	USE physics_constants, ONLY : TODEG

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!---------------------------
	!	Inputs	
	!---------------------------

	REAL(KIND=8)										:: theta	
	
	!---------------------------
	!	Internal	
	!---------------------------

	REAL(KIND=8)										:: ReAmp
	REAL(KIND=8)										:: ImAmp
	REAL(KIND=8)										:: Dummy_Amp
	REAL(KIND=8)										:: Coeff 
	REAL(KIND=8)										:: my_ReBuff 
	REAL(KIND=8)										:: my_ImBuff

	REAL(KIND=8),DIMENSION(numrank) :: ReBuff 
	REAL(KIND=8),DIMENSION(numrank) :: ImBuff 

	INTEGER													:: chunk
	INTEGER													:: i 
	INTEGER													:: L 

	CALL legendre(COS(theta),LMax)	

	chunk = LMax/numrank
	ReAmp = 0.0
	ImAmp = 0.0

	DO i=0,(chunk-1)
		L     = rank*chunk + i
		CALL idnt_coeff(L,leg(L+1),phase(L),Coeff)
		ReAmp = ReAmp + Coeff*COS(phase(L))	
		ImAmp = ImAmp + Coeff*SIN(phase(L))	
	END DO
	
	my_ReBuff = ReAmp; 
	my_ImBuff = ImAmp; 

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( my_ReBuff, 1, MPI_REAL8, ReBuff, &
	& 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( my_ImBuff, 1, MPI_REAL8, ImBuff, &
	& 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )	
	
	IF (rank == 0) THEN
		ReAmp = 0.0
		ImAmp = 0.0

		DO i=1,numrank
			ReAmp = ReAmp + ReBuff(i)
			ImAmp = ImAmp + ImBuff(i)
		END DO

		Amp = (1.0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)
		
		IF (FRAME == 1) THEN
			CALL amp_to_lab(theta, Amp, Dummy_Amp)		
			Amp = Dummy_Amp
		END IF
	
	END IF	

	DEALLOCATE(leg)

END SUBROUTINE pw_amp_MPI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE pw_amp(theta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes scattering amplitude using partial
! wave analysis. Phases must be computed prior
! to computing amplitude. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : leg, phase, K, E, LMax, S_phase, T_phase, S_LMax, T_LMax
	USE dcs, ONLY : Amp
	USE collision_inputs, ONLY : FRAME, Proj, Targ, IDNT, RE_IM_WRITE	
	USE physics_constants, ONLY : TODEG, TOEV

	IMPLICIT NONE

	!---------------------------
	!	Inputs	
	!---------------------------

	REAL(KIND=8)										:: theta	
	
	!---------------------------
	!	Internal	
	!---------------------------

	REAL(KIND=8)										:: E_ReAmp
	REAL(KIND=8)										:: O_ReAmp
	REAL(KIND=8)										:: E_ImAmp
	REAL(KIND=8)										:: O_ImAmp
	REAL(KIND=8)										:: ReAmp
	REAL(KIND=8)										:: ImAmp
	REAL(KIND=8)										:: Dummy_Amp
	REAL(KIND=8)										:: Coeff 
	REAL(KIND=8)										:: T_Amp
	REAL(KIND=8)										:: S_Amp
	REAL(KIND=8)										:: S_ReAmp, S_ImAmp, T_ReAmp, T_ImAmp

	INTEGER													:: i 
	INTEGER													:: L 

	IF (RE_IM_WRITE .EQ. 1) OPEN(UNIT=45,FILE="../Data/REAL_IMAG_AMP.dat",ACCESS="APPEND")

	!! H+H collision, multichannel
	IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN
		IF (S_LMax .GT. T_LMax) THEN
			CALL legendre(COS(theta),S_LMax)
		ELSE
			CALL legendre(COS(theta),T_LMax)
		END IF	

		!! Do symmetrized amplitudes if IDNT NE 0
		IF (IDNT .NE. 0) THEN

			E_ReAmp = 0.0D0
			O_ReAmp = 0.0D0
			E_ImAmp = 0.0D0
			O_ImAmp = 0.0D0

			!!! Singlet Channel
			DO i=0,(S_LMax-1)
				L = i
				!!! Even L
				IF (MOD(L,2) .EQ. 0) THEN			
					E_ReAmp = E_ReAmp + (2.0D0*L+1.0D0)*leg(L+1)*COS(S_phase(L+1))*SIN(S_phase(L+1))	
					E_ImAmp = E_ImAmp + (2.0D0*L+1.0D0)*leg(L+1)*SIN(S_phase(L+1))*SIN(S_phase(L+1))	
				!!! Odd L
				ELSE
					O_ReAmp = O_ReAmp + (2.0D0*L+1.0D0)*leg(L+1)*COS(S_phase(L+1))*SIN(S_phase(L+1))	
					O_ImAmp = O_ImAmp + (2.0D0*L+1.0D0)*leg(L+1)*SIN(S_phase(L+1))*SIN(S_phase(L+1))	
				END IF	
			END DO ! i

			S_Amp = (0.5D0/(K*K))*(E_ReAmp*E_ReAmp + E_ImAmp*E_ImAmp) + &
			&     (3.0D0/(2.0D0*K*K))*(O_ReAmp*O_ReAmp + O_ImAmp*O_ImAmp)

			E_ReAmp = 0.0D0
			O_ReAmp = 0.0D0
			E_ImAmp = 0.0D0
			O_ImAmp = 0.0D0

			!!! Triplet Channel
			DO i=0,(T_LMax-1)
				L = i
				!!! Even L
				IF (MOD(L,2) .EQ. 0) THEN			
					E_ReAmp = E_ReAmp + (2.0D0*L+1.0D0)*leg(L+1)*COS(T_phase(L+1))*SIN(T_phase(L+1))	
					E_ImAmp = E_ImAmp + (2.0D0*L+1.0D0)*leg(L+1)*SIN(T_phase(L+1))*SIN(T_phase(L+1))	
				!!! Odd L
				ELSE
					O_ReAmp = O_ReAmp + (2.0D0*L+1.0D0)*leg(L+1)*COS(T_phase(L+1))*SIN(T_phase(L+1))	
					O_ImAmp = O_ImAmp + (2.0D0*L+1.0D0)*leg(L+1)*SIN(T_phase(L+1))*SIN(T_phase(L+1))	
				END IF	
			END DO

			T_Amp = (3.0D0/(2.0D0*K*K))*(E_ReAmp*E_ReAmp + E_ImAmp*E_ImAmp) + &
			&     (0.5D0/(K*K))*(O_ReAmp*O_ReAmp + O_ImAmp*O_ImAmp)

			Amp = 0.25D0*S_Amp + 0.75D0*T_Amp

		!!! Do non-symmetrized amplitudes if IDNT EQ 0	
		ELSE	
			
			!!! S channel
			ReAmp = 0.0D0
			ImAmp = 0.0D0		
	
			Do i=0,(S_LMax-1)
				L = i
	      CALL idnt_coeff(L,leg(L+1),S_phase(L+1),Coeff)
 	     	ReAmp = ReAmp + Coeff*COS(S_phase(L+1))
 	     	ImAmp = ImAmp + Coeff*SIN(S_phase(L+1))
			END DO

			S_Amp = (1.0D0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)

			!!! T channel
			ReAmp = 0.0D0
			ImAmp = 0.0D0		
	
			Do i=0,(T_LMax-1)
				L = i
	      CALL idnt_coeff(L,leg(L+1),T_phase(L+1),Coeff)
 	     	ReAmp = ReAmp + Coeff*COS(T_phase(L+1))
 	     	ImAmp = ImAmp + Coeff*SIN(T_phase(L+1))
			END DO

			T_Amp = (1.0D0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)

			Amp   = 0.25D0*S_Amp + 0.75D0*T_Amp

		END IF ! H+H

	!! He+O collision, multichannel
	ELSE IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN

		IF (S_LMax .GT. T_LMax) THEN
			CALL legendre(COS(theta),S_LMax)
		ELSE
			CALL legendre(COS(theta),T_LMax)
		END IF	

		ReAmp = 0.0D0
		ImAmp = 0.0D0

		!!! Singlet (PI) Channel
		DO i=0,(S_LMax-1)
			L = i
			CALL idnt_coeff(L,leg(L+1),S_phase(L+1),Coeff)
			ReAmp = ReAmp + Coeff*COS(S_phase(L+1))
			ImAmp = ImAmp + Coeff*SIN(S_phase(L+1))
		END DO
	
		S_ReAmp = ReAmp 
		S_ImAmp = ImAmp
		S_Amp   = (1.0D0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)

		ReAmp = 0.0D0
		ImAmp = 0.0D0

		!!! Triplet (SIGMA) Channel
		DO i=0,(T_LMax-1)
			L = i
			CALL idnt_coeff(L,leg(L+1),T_phase(L+1),Coeff)
			ReAmp = ReAmp + Coeff*COS(T_phase(L+1))
			ImAmp = ImAmp + Coeff*SIN(T_phase(L+1))
		END DO

		T_ReAmp = ReAmp
		T_ImAmp = ImAmp
		T_Amp   = (1.0D0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)

		Amp   = (2.0D0/3.0D0)*S_Amp + (1.0D0/3.0D0)*T_Amp	

		IF (RE_IM_WRITE .EQ. 1) WRITE(45,*) E*TOEV, ((2.0D0/3.0D0)*S_ReAmp+(1.0D0/3.0D0)*T_ReAmp), &
		& ((2.0D0/3.0D0)*S_ImAmp+(1.0D0/3.0D0)*T_ImAmp)

	!!! If not multichannel
	ELSE
		ReAmp = 0.0D0
		ImAmp = 0.0D0

		CALL legendre(COS(theta),LMax)	

		DO i=0,(LMax - 1)
			L     =  i
			CALL idnt_coeff(L,leg(L+1),phase(L+1),Coeff)
			ReAmp = ReAmp + Coeff*COS(phase(L+1))	
			ImAmp = ImAmp + Coeff*SIN(phase(L+1))	
		END DO

		IF (RE_IM_WRITE .EQ. 1) WRITE(45,*) E*TOEV, ReAmp, ImAmp
	
		Amp = (1.0D0/(K*K))*(ReAmp*ReAmp + ImAmp*ImAmp)

	END IF ! multi-channel or not
	
	IF (FRAME == 1) THEN
		CALL amp_to_lab(theta, Amp, Dummy_Amp)		
		Amp = Dummy_Amp
	END IF

	DEALLOCATE(leg)

	IF (RE_IM_WRITE .EQ. 1) CLOSE(45)

END SUBROUTINE pw_amp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE idnt_coeff(L, L_leg, L_phase, Coeff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the proper type of coefficient for use
! by pw_amp based on the identical partical 
! symmetry of the collision being computed. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : IDNT

	IMPLICIT NONE

	!---------------------------
	!	Inputs	
	!---------------------------

	INTEGER										:: L

	REAL(KIND=8)							:: L_leg
	REAL(KIND=8)							:: L_phase

	!---------------------------
	!	Outputs	
	!---------------------------

	REAL(KIND=8)							:: Coeff
	
	!---------------------------
	!	Internal	
	!---------------------------

	IF (IDNT == 0) THEN
		Coeff = (2.0*L+1.0)*L_leg*SIN(L_phase)
	ELSE IF (IDNT == 1) THEN
		Coeff = (2.0*L+1.0)*L_leg*SIN(L_phase)*(1.0-(-1.0)**L)
	ELSE
		Coeff = (2.0*L+1.0)*L_leg*SIN(L_phase)*(1.0+(-1.0)**L)
	END IF

END SUBROUTINE idnt_coeff
