
SUBROUTINE critical_angle_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find the critical angle which has the 
! largest contribution to the total cross
! section.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, TODEG, TOEV
	USE collision_inputs, ONLY : FRAME
	USE current_energy, ONLY : E
	
	IMPLICIT NONE

  !----------------------------
  ! Internal  
  !----------------------------

	REAL(KIND=8)								:: Theta	! Current Angle
	REAL(KIND=8)								:: DT			! Differential Angle
	REAL(KIND=8)								:: TI			! Initial Angle 
	REAL(KIND=8)								:: TF			! Final Angle
	REAL(KIND=8)								:: MaxPD	! Maximum Prob Den 
	REAL(KIND=8)								:: Crit		! Critical Angle 
	REAL(KIND=8)								:: Dum		! Dummy Angle 
	REAL(KIND=8)								:: PD			! Current Prob Den 

	INTEGER											:: NT			! Number of Angles
	INTEGER											:: i			! Counter 

	MaxPD = 0.0
	NT    = 10000	
	TI 		= 0.0
	TF 		= 10.0*PI/180.0	
	DT  	= (TF-TI)/DBLE(NT)

	DO i=0,NT
		Theta = TI + i*DT	
		CALL prob_density(Theta, PD)
	!	WRITE(*,*) Theta*180.0/PI, PD	
		IF (PD .GT. MaxPD) THEN
			MaxPD = PD
			Crit  = Theta
		END IF
	END DO	

	IF (FRAME == 1) THEN
		CALL angle_to_lab(Crit,Dum)
		Crit = Dum
	END IF

	OPEN(UNIT=33, FILE="../Data/CritAngle.dat", ACCESS="APPEND")
	WRITE(33,*) E*TOEV, Crit*TODEG 
	CLOSE(33)

END SUBROUTINE


