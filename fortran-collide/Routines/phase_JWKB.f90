
SUBROUTINE phase_JWKB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine calculates the phaseshift using
! the semi-classical approximation for high 
! energy scattering. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : K, E, phase
	USE physics_constants, ONLY : TOEV, PI
	USE grid, ONLY : RStart

	IMPLICIT NONE

	REAL(KIND=8)							:: bi
	REAL(KIND=8)							:: bf
	REAL(KIND=8)							:: db
	REAL(KIND=8)							:: b
	REAL(KIND=8)							:: r0
	REAL(KIND=8)							:: I1

	INTEGER										:: Nb	
	INTEGER										:: i


	bi = 1.0/(2.0*K)
	bf = 10.0
	Nb = 100
	db = (bf-bi)/DBLE(Nb-1)

	ALLOCATE(phase(Nb))

	CALL RStart_finder
	r0 = RStart

	OPEN(UNIT=666, FILE="../Data/JWKB_Phases.dat", ACCESS="APPEND")

	DO i=1,Nb
		b = bi + (i-1)*db	
		CALL I1_JWKB(b,I1)
		phase(i) = K*b*(PI/2.0 - r0/b + I1)
		WRITE(666,*) b, phase(i)
	END DO	

	CLOSE(666)

	WRITE(*,*) "K is ", K
	
	DEALLOCATE(phase)

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE I1_JWKB(b,I1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the I1 integral needed for the JWKB
! semi-classical phaseshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE grid, ONLY : RStart
	USE current_energy, ONLY : E

	IMPLICIT NONE

	!!!!!!!!!!!!!!!!!!!
	! Inputs	
	!!!!!!!!!!!!!!!!!!!
	
	REAL(KIND=8)							:: b

	!!!!!!!!!!!!!!!!!!!
	! Outputs	
	!!!!!!!!!!!!!!!!!!!

	REAL(KIND=8)							:: I1

	!!!!!!!!!!!!!!!!!!!
	! Internal	
	!!!!!!!!!!!!!!!!!!!

	REAL(KIND=8)							:: r0			
	REAL(KIND=8)							:: rf			
	REAL(KIND=8)							:: dr			
	REAL(KIND=8)							:: r			
	REAL(KIND=8)							:: grand			
	REAL(KIND=8)							:: pot		
	REAL(KIND=8)							:: C		
	REAL(KIND=8)							:: C1		
	REAL(KIND=8)							:: C2		
	REAL(KIND=8)							:: X		

	INTEGER										:: Nr
	INTEGER										:: i

	Nr 		= 1000
	r0 		= RStart
	rf 		= 25.0
	dr 		= (rf-r0)/DBLE(Nr-1)	
	grand = 0.0	

	DO i=1,Nr	
		r  		= r0 + (i-1)*dr
		CALL potential(r,pot)
		C1 		= pot/E
		C2 		= b*b/(r*r)
		X     = 1.0 - C1 - C2
		IF (X .GT. 0.0) THEN
			C  		= SQRT(X) - 1.0
!		WRITE(*,*) "X ", 1.0-C1-C2, "sqrt(X) ", SQRT(1.0-C1-C2)
			grand = grand + C*dr
		END IF
!		WRITE(*,*) "grand ", grand
	END DO

	I1 = grand/b
!	WRITE(*,*) "b ", b, "I1 ", I1

END SUBROUTINE

