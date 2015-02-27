
SUBROUTINE prob_density_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the probability density for a given
! collision. The probability density is given as
!
!				pd(theta) = 2*pi*sin(theta)*Amp/TCS
!
!	Integrating the probability density over all
! angles results in 1.0
!
! Writes the angle and prob density to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, TOEV, TODEG
	USE dcs, ONLY : Amp
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : FRAME, ThI, ThF, NTh

	IMPLICIT NONE

 	!----------------------------
  ! Internal  
  !----------------------------

	REAL(KIND=8)		 :: PD				! Probability Density
	REAL(KIND=8)		 :: TC				! Total Cross Section
	REAL(KIND=8)		 :: Theta		! Current Angle
	REAL(KIND=8)		 :: D_Theta	! Dummy Current Angle
	REAL(KIND=8)		 :: DT				! Differential Angle 
	REAL(KIND=8)		 :: Tot			! Integration sum 
	INTEGER					 :: i				! Counter 
	CHARACTER(LEN=29):: Filename ! Filename to write
	
	DT = (ThF-ThI)/DBLE(NTh)

	CALL total_cross_section(TC)

	WRITE(UNIT=Filename, FMT="(A, ES7.1, A)") "../Data/ProbDen_", E*TOEV, "eV.dat"
	OPEN(UNIT=44, FILE=Filename, ACCESS="APPEND")
	OPEN(UNIT=55, FILE="../Data/Integrated_ProbDen.dat", ACCESS="APPEND")

	Tot = 0.0D0
	
	DO i=0,NTh
		Theta = ThI + i*DT
		CALL pw_amp(Theta)		
		IF (FRAME == 1) THEN
			CALL angle_to_lab(Theta,D_Theta)
			Theta = D_Theta
		END IF	

		PD 	= 2.0*PI*SIN(Theta)*Amp/TC
		Tot = Tot + PD*DT
		WRITE(44,*) Theta*TODEG, PD	
	END DO

	WRITE(55,*) E*TOEV, Tot

	CLOSE(44)
	CLOSE(55)

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE prob_density_3D_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the probability density for a given
! collision. The probability density is given as
!
!				pd(theta) = 2*pi*sin(theta)*Amp/TCS
!
!	Integrating the probability density over all
! angles results in 1.0
!
! Writes the energy, angle and prob density to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, TOEV, TODEG
	USE dcs, ONLY : Amp
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : FRAME, ThI, ThF, NTh

	IMPLICIT NONE

 	!----------------------------
  ! Internal  
  !----------------------------

	REAL(KIND=8)								:: PD				! Probability Density
	REAL(KIND=8)								:: TC				! Total Cross Section
	REAL(KIND=8)								:: Theta		! Current Angle
	REAL(KIND=8)								:: D_Theta	! Dummy Current Angle
	REAL(KIND=8)								:: DT				! Differential Angle 

	INTEGER											:: i				! Counter 

	DT = (ThF-ThI)/DBLE(NTh)

	CALL total_cross_section(TC)
	OPEN(UNIT=44, FILE="../Data/ProbDen3D.dat", ACCESS="APPEND")
	
	DO i=0,NTh
		Theta = ThI + i*DT
		CALL pw_amp(Theta)		
		IF (FRAME == 1) THEN
			CALL angle_to_lab(Theta,D_Theta)
			Theta = D_Theta
		END IF	

		PD = 2.0*PI*SIN(Theta)*Amp/TC
		WRITE(44,*) E*TOEV, Theta*TODEG, PD	
	END DO

	CLOSE(44)


END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE prob_density(Theta, PD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the probability density for a given
! collision. The probability density is given as
!
!				pd(theta) = 2*pi*sin(theta)*Amp/TCS
!
!	Integrating the probability density over all
! angles results in 1.0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, TOEV, TODEG
	USE dcs, ONLY : Amp
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : FRAME

	IMPLICIT NONE

 	!----------------------------
  ! Inputs 
  !----------------------------

	REAL(KIND=8)		:: Theta		! Current Angle
	
	!----------------------------
  ! Outputs  
  !----------------------------
 
	REAL(KIND=8)		:: PD				! Probability Density

	!----------------------------
  ! Internal  
  !----------------------------

	REAL(KIND=8)		:: TC				! Total Cross Section
	REAL(KIND=8)		:: D_Theta	! Dummy Current Angle
	
!	CALL total_cross_section(TC)
	CALL TCS_HH( E*TOEV, TC )
	CALL pw_amp(Theta)		
	IF (FRAME == 1) THEN
		CALL angle_to_lab(Theta,D_Theta)
		Theta = D_Theta
	END IF	

	PD = 2.0*PI*SIN(Theta)*Amp/TC

END SUBROUTINE

