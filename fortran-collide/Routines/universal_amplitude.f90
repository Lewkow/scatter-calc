!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Use the universal amplitude function to
! obtain the differential cross section
!
! Theta [deg], E [eV], Mu [au]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE uni_amp( E, Mu, theta, DCS )
USE physics_constants, ONLY : PI, TOEV

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E, Mu, Theta

	!! Outputs
	REAL(KIND=8)		:: DCS

	!! Implicit
	REAL(KIND=8)		:: a1, a2, a3, a4, x, tot, dt
	INTEGER					:: i, N

	a1 = -0.136D0
	a2 =  0.993D0
	a3 = -3.042D0
	a4 =  5.402D0

	IF (Theta .EQ. 0.0D0) THEN
		DCS = 0
	ELSE
		x 	= LOG( (E*Theta)/Mu )
		tot = a1*x**3 + a2*x**2 + a3*x + a4 - LOG( Theta*SIN(Theta*PI/180.0D0) )
		DCS = EXP(tot)	
	END IF

END SUBROUTINE uni_amp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Use the universal amplitude function to
! obtain the differential cross section, 
! integrate and obtain the total cross section 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE uni_TCS( E, Mu, TCS )
USE physics_constants, ONLY : PI, TOEV, TOAMU
USE mpi

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E, Mu

	!! Outputs
	REAL(KIND=8)		:: TCS

	!! Implicit
	REAL(KIND=8)		:: x, x_min, x_max, tot, dt, C, t, a, t_min, t_max
	INTEGER					:: i, N

	x_min = 400.0D0
	x_max = 1.0D4
	t_min = ( (Mu/TOAMU)/(E*TOEV) )*(x_min)
	t_max = ( (Mu/TOAMU)/(E*TOEV) )*(x_max) 
	N   	= 100000													! Number of integration steps
	dt  	= PI/REAL(N) 											! Numeric Theta Differential
!	dt  	= (t_max-t_min)/REAL(N) 					! Numeric Theta Differential
	tot 	= 0.0D0	

!	write(*,*) "t_min: t_max: ", t_min, t_max

!	IF ( rank .EQ. 0) THEN
		DO i=0,N
			t = i*dt
			CALL uni_amp(	E*TOEV, Mu/TOAMU, t*180.0D0/PI, a )
			tot = tot + 2.0D0*PI*dt*SIN(t)*a
			write(99,'(6ES10.2)') E*TOEV, Mu/TOAMU, t*180.0D0/PI, LOG((E*TOEV*t*180.0D0/PI)/(Mu/TOAMU)), a, tot
		END DO

		TCS = tot
	
!		WRITE(*,'(A,ES10.2)') "Integrated Universal Amplitude TCS: " , TCS

END SUBROUTINE uni_TCS

