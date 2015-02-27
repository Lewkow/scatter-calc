
SUBROUTINE pot_HeH_Uni(R,pot)

	!--------------------
	! Inputs	
	!--------------------
		
	REAL(KIND=8)		:: R

	!--------------------
	! Outputs	
	!--------------------
 
	REAL(KIND=8)		:: pot	
 
	!--------------------
	! Internal	
	!--------------------
 
	REAL(KIND=8)		:: Z1
	REAL(KIND=8)		:: Z2
	REAL(KIND=8)		:: a 
	REAL(KIND=8)		:: x 
	REAL(KIND=8)		:: Screen 
 
	Z1 = 2.0
  Z2 = 1.0
  a = 0.88534/(Z1**0.23+Z2**0.23)
  x = R/a
  Screen = 0.1818*EXP(-3.2*x)+0.5099*EXP(-0.9423*x) &
	& +0.2802*Exp(-0.4029*x)+0.02817*EXP(-0.2016*x)
  pot = Z1*Z2*Screen/R

END SUBROUTINE 
