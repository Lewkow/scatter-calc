
MODULE collision_inputs

	CHARACTER(LEN=3)												:: Proj					! Projectile H, He3, He4, O
	CHARACTER(LEN=3)												:: Targ         ! Target H, He3, He4, O
	REAL(KIND=8)														:: M1		   			! Projectile Mass
	REAL(KIND=8)														:: M2		   			! Target Mass
	REAL(KIND=8)														:: Mu		   			! Reduced Mass
	REAL(KIND=8)														:: Depth   			! Depth of Potential Well
	REAL(KIND=8)														:: ThI     			! Initial DCS angle
	REAL(KIND=8)														:: ThF     			! Final DCS angle

	INTEGER																	:: TCS     			! Turn TCS 		on/off (1/0)
	INTEGER																	:: NUM_TCS 			! Turn NUM_TCS	on/off (1/0)
	INTEGER																	:: DFCS    			! Turn DFCS		on/off (1/0)
	INTEGER																	:: DCS     			! Turn DCS 		on/off (1/0)
	INTEGER																	:: DC3D					! TURN DC3D   on/off (1/0)	
	INTEGER																	:: HS_TEST 			! Turn HS_TEST on/off (1/0)
	INTEGER																	:: POT_TEST			! Turn POT_TEST on/off (1/0)
	INTEGER																	:: PROB_DEN			! Turn PROB_DEN on/off (1/0)
	INTEGER																	:: PROB_DEN_FIT	! Turn PROB_DEN_FIT on/off (1/0)
	INTEGER																	:: PROB_DEN_3D	! Turn PROB_DEN on/off (1/0)
	INTEGER																	:: AVE_EN_LOSS	! Turn AVE_EN_LOSS on/off (1/0)
	INTEGER																	:: PHASE_WRITE	! Write Phases on/off (1/0)
	INTEGER																	:: PHASE_TIMING ! Write Phase Times on/off (1/0)
	INTEGER																	:: FRAME   			! Frame to compute in
																								     			! CM -> 0    Lab -> 1
	INTEGER																	:: IDNT    			! Identical particle symmetry
																						         			! NonIdnt -> 0
																									   			! Fermion -> 1   Boson -> 2
	INTEGER																	:: NTh     			! Num DCS angles to compute
	INTEGER																	:: CALC_TYPE		! Type of calculation
																													! PW -> 0  EA -> 1  Both -> 99
	INTEGER																	:: CurEngy 			! Current energy counter
	INTEGER																	:: NumEngy 			! Num energies to compute
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Engy    			! Array to hold energies
	REAL(KIND=8)														:: EngyI				! Initial Energy for Loop
	REAL(KIND=8)														:: EngyF				! Final Energy for Loop
	INTEGER																	:: EngyN				! Number of Energies for Loop

	REAL(KIND=8)														:: POT_I				! Initial R for POT_TEST
	REAL(KIND=8)														:: POT_F				! Final R for POT_TEST
	INTEGER																	:: POT_N				! Number of points for POT_TEST
	INTEGER,PARAMETER												:: RE_IM_WRITE=1


	CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	SUBROUTINE mass_finder
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Finds the mass of the target and 
	! projectile based on the char inputs
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!! Projectile
	IF ( Proj == 'H' ) THEN
		M1 = 1.00782503207
	ELSE IF ( Proj == 'He3' ) THEN
		M1 = 3.0160293191
	ELSE IF ( Proj == 'He4' ) THEN
		M1 = 4.00260325415 
	ELSE IF ( Proj == 'O' ) THEN
		M1 = 15.99491461956
	ELSE
		WRITE(*,*) "Projectile ", Proj, " not known!!!!"
	END IF

	!!! Target
	IF ( Targ == 'H' ) THEN
		M2 = 1.00782503207
	ELSE IF ( Targ == 'He3' ) THEN
		M2 = 3.0160293191
	ELSE IF ( Targ == 'He4' ) THEN
		M2 = 4.00260325415 
	ELSE IF ( Targ == 'O' ) THEN
		M2 = 15.99491461956
	ELSE
		WRITE(*,*) "Target ", Targ, " not known!!!!"
	END IF

	END SUBROUTINE mass_finder

END MODULE collision_inputs
