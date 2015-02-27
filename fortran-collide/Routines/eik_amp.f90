
SUBROUTINE eik_amp( theta )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the scattering amplitude using the 
! Eikonal approximation. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : M1, M2, Mu, FRAME
	USE current_energy, ONLY : K
	USE eik, ONLY : AmpInt, Eik_RStart
	USE grid, ONLY : REnd	
	USE dcs, ONLY : Amp

	IMPLICIT NONE

!---------------------------
! Inputs
!---------------------------

	REAL(KIND=8)										:: theta    	! Angle to compute ampltude for

!---------------------------
! Internal
!---------------------------

	REAL(KIND=8)  									:: delta    	! Change in b for each integration step
	REAL(KIND=8)										:: RE					! Real amplitude
	REAL(KIND=8)										:: IM					! Imaginary amplitude
	REAL(KIND=8)										:: Dummy_Amp 	! Dummy amplitude
	REAL(KIND=8)										:: b					! Impact parameter
	REAL(KIND=8)										:: J    			! Coefficient dependent on IDNT 
	REAL(KIND=8)										:: P					! Eikonal phase

	INTEGER													:: i        	! Counter

	delta    = (REnd - Eik_RStart)/DBLE(AmpInt)

	RE = 0.0D0
	IM = 0.0D0
	
	DO i=0,AmpInt
		b = Eik_RStart + i*delta	
		CALL eik_amp_coeff(theta,b,J)
		CALL eik_phase(b,P)
		RE = RE + b*J*COS(P)*SIN(P)
		IM = IM + b*J*SIN(P)*SIN(P)
	END DO

	RE  = 2.0*K*delta*RE
	IM  = 2.0*K*delta*IM	
	Amp = RE*RE+IM*IM	

	IF (FRAME == 1) THEN
		CALL amp_to_lab(theta, Amp, Dummy_Amp)
		Amp = Dummy_Amp
	END IF

END SUBROUTINE eik_amp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE eik_amp_MPI( theta )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes the scattering amplitude using the 
! Eikonal approximation. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE collision_inputs, ONLY : M1, M2, Mu, FRAME
	USE current_energy, ONLY : K
	USE eik, ONLY : AmpInt, Eik_RStart
	USE grid, ONLY : REnd	
	USE dcs, ONLY : Amp
	INCLUDE 'mpif.h'

!---------------------------
! Inputs
!---------------------------

	REAL(KIND=8)										:: theta    	! Angle to compute ampltude for

!---------------------------
! Internal
!---------------------------

	REAL(KIND=8)										:: chunk    	! Chunk of impact param for rank
	REAL(KIND=8)										:: my_Amp   	! Ranks amplitude 
	REAL(KIND=8)  									:: delta    	! Change in b for each integration step
	REAL(KIND=8)										:: RE					! Real amplitude
	REAL(KIND=8)										:: IM					! Imaginary amplitude
	REAL(KIND=8)										:: Dummy_Amp 	! Dummy amplitude
	REAL(KIND=8)										:: b					! Impact parameter
	REAL(KIND=8)										:: J    			! Coefficient dependent on IDNT 
	REAL(KIND=8)										:: P					! Eikonal phase
	REAL(KIND=8),DIMENSION(numrank)	:: buff				! Buffer used to gather amplitude by root

	INTEGER													:: i        	! Counter
	INTEGER													:: my_start 	! Ranks starting b value
	INTEGER													:: my_end   	! Ranks ending b value

	
	delta    = (REnd - Eik_RStart)/DBLE(AmpInt)
	chunk    = AmpInt/numrank
	my_start = rank*chunk
	my_end   = (rank+1)*chunk - 1

	IF ( rank == (numrank - 1) ) THEN
		my_end = AmpInt
	END IF

	RE = 0.0D0
	IM = 0.0D0
	
	DO i=my_start,my_end
		b = Eik_RStart + i*delta	
		CALL eik_amp_coeff(theta,b,J)
		CALL eik_phase(b,P)
		RE = RE + b*J*COS(P)*SIN(P)
		IM = IM + b*J*SIN(P)*SIN(P)
	END DO

	RE     = 2.0*K*delta*RE
	IM     = 2.0*K*delta*IM	
	my_Amp = RE*RE+IM*IM	

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( my_Amp, 1, MPI_REAL8, &
	& buff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr ) 

	IF (rank == 0) THEN
		Amp = 0.0
		DO i=1,numrank
			Amp = Amp + buff(i)
		END DO
		IF (FRAME == 1) THEN
			CALL amp_to_lab(theta, Amp, Dummy_Amp)
			Amp = Dummy_Amp
		END IF
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE eik_amp_MPI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE eik_amp_coeff(theta,b,J)
	USE collision_inputs, ONLY : IDNT
	USE current_energy, ONLY : K 
	USE physics_constants, ONLY : PI
!---------------------------
! Inputs
!---------------------------

	REAL(KIND=8)				:: theta	! Angle for amplitude 
	REAL(KIND=8)				:: b    	! Impact parameter 

!---------------------------
! Outputs
!---------------------------

	REAL(KIND=8)				:: J	    ! Coefficient dependent on IDNT

!---------------------------
! Internal 
!---------------------------

	REAL(KIND=8)				:: X1			! Argument for bessj0
	REAL(KIND=8)				:: X2			! Argument for bessj0
	REAL(KIND=8)				:: Y1			! Bessj0  
	REAL(KIND=8)				:: Y2			! Bessj0 

	X1 = 2.0*K*b*SIN(theta/2.0)
	X2 = 2.0*K*b*SIN((PI-theta)/2.0)

	IF (IDNT == 0) THEN
		CALL bessj0(X1,J)
	ELSE IF (IDNT == 1) THEN
		CALL bessj0(X1,Y1)
		CALL bessj0(X2,Y2)
		J = Y1 - Y2
	ELSE IF (IDNT == 2) THEN
		CALL bessj0(X1,Y1)
		CALL bessj0(X2,Y2)
		J = Y1 + Y2
	END IF
	
END SUBROUTINE eik_amp_coeff

