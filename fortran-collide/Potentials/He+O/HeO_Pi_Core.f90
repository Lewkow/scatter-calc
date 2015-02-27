
SUBROUTINE pot_HeO_Pi(R,pot) 
 	USE physics_constants, ONLY : TOEV 
	!--------------------
  ! Inputs  
  !--------------------

  REAL(KIND=8)    :: R

  !--------------------
  ! Outputs 
  !--------------------

  REAL(KIND=8)    :: pot

  !--------------------
  ! Internal  
  !--------------------

	REAL(KIND=8)								:: RR
	REAL(KIND=8)								:: RROld 
	REAL(KIND=8)								:: R4
	REAL(KIND=8)								:: R6
	REAL(KIND=8)								:: R8
	REAL(KIND=8)								:: R10
	REAL(KIND=8)								:: R12
	REAL(KIND=8)								:: R14
	REAL(KIND=8)								:: Alpha
	REAL(KIND=8)								:: Beta
	REAL(KIND=8)								:: Gama
	REAL(KIND=8)								:: Lamda 
	REAL(KIND=8)								:: T1 
	REAL(KIND=8)								:: T2 
	REAL(KIND=8)								:: T3 
	REAL(KIND=8)								:: T31 
	REAL(KIND=8)								:: T4 
	REAL(KIND=8)								:: TOCM 

	REAL(KIND=8)								:: A 
	REAL(KIND=8)								:: aa 
	REAL(KIND=8)								:: z1 
	REAL(KIND=8)								:: z2 

	REAL(KIND=8),DIMENSION(14)	:: Coeff

	IF (R .LT. 2.376) THEN
		A  = 20.29143090894122	
		z1 = 2.0
		z2 = 8.0
		aa = .88534/(z1**0.23+z2**0.23) 
		pot = (A*z1*z2/R)*(0.1818*exp(-3.2*R/aa)+0.5099*exp(-0.9423*R/aa)&
		&+0.2802*exp(-0.4029*R/aa)+0.02817*exp(-0.2016*R/aa))
    pot = pot/TOEV
	ELSE
		TOCM 			= 219474.63137
  	Alpha 		= 1.71319289534039
  	Beta 			= 0.976802973556596
  	Gama 			= -3.81495290352935
  	Lamda 		= -0.393324961747769
  	RR 				= 1.0D0

  	Coeff(1) 	= -510844.473121867
  	Coeff(2) 	= 1116784.56576381
  	Coeff(3) 	= -611547.840972992
  	Coeff(4) 	= 164530.896515184
  	Coeff(5) 	= -24974.4379458087
  	Coeff(6) 	= 1968.3178153098
  	Coeff(7) 	= -36.308471139056
  	Coeff(8) 	= -5.3168142290825
  	Coeff(9) 	= 0.260851069082499
  	Coeff(10) = -1029336.021
  	Coeff(11) = -2886.67864863542
  	Coeff(12) = -2519820024.78251
  	Coeff(13) = 30818852577.5822
  	Coeff(14) = -60898527680.2877

  	T1   			= EXP(-Alpha * (R - Beta))
  	T31  			= TANH(Gama - Lamda * R)
  	T3   			= 0.5 * (1.0 + T31)
  	T2   			= 0.0D0

		DO i=0,8
    	IF (i == 4) THEN
				R4 = RR
			END IF
    
			IF (i == 6) THEN
				R6 = RR
			END IF
    
			IF (i == 8) THEN
				R8 = RR
			END IF

    	T2 		= T2 + RR*Coeff(i+1)
    	RROld = RR
    	RR 		= R*RR
		END DO ! i

  	R10 = R6 * R4;
  	R12 = R8 * R4;
  	R14 = R6 * R8;
  	T4  = Coeff(10)/R6 + Coeff(11)/R8 + Coeff(12)/R10 + Coeff(13)/R12 + Coeff(14)/R14
  	pot = (T1*T2 + T3*T4)/TOCM
	END IF

END SUBROUTINE 
