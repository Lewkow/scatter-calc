
SUBROUTINE pot_HeHe(R, Pot)
	USE physics_constants, ONLY : TOEV, TOK
  !--------------------
  ! Inputs  
  !--------------------

  REAL(KIND=8)    						:: R

  !--------------------
  ! Outputs 
  !--------------------

  REAL(KIND=8)    						:: Pot

  !--------------------
  ! Internal  
  !--------------------

	REAL(KIND=8)								:: A
	REAL(KIND=8)								:: a1 
	REAL(KIND=8)								:: a_1 
	REAL(KIND=8)								:: a2 
	REAL(KIND=8)								:: a_2 
	REAL(KIND=8)								:: d1 
	REAL(KIND=8)								:: d2 
	REAL(KIND=8)								:: d3 
	REAL(KIND=8)								:: b 
	REAL(KIND=8)								:: ep_kb 
	REAL(KIND=8)								:: R_ep 
	REAL(KIND=8)								:: sig 
	REAL(KIND=8)								:: fact 
	REAL(KIND=8)								:: ret0 
	REAL(KIND=8)								:: ret1 
	REAL(KIND=8)								:: ret2 
	REAL(KIND=8)								:: atom 
	REAL(KIND=8)								:: aa 
	REAL(KIND=8)								:: Z 

	REAL(KIND=8),DIMENSION(17)	:: C 

	INTEGER											:: n
	INTEGER											:: k 
	INTEGER											:: i 

	A 		=  0.307092338615e7   ! Units K
  a1  	= -0.201651289932e1   ! Units 1/a0
  a_1 	= -0.431646276045e0   ! Units a0
  a2  	= -0.459521265125e-1  ! Units 1/a0/a0
  a_2 	=  0.138539045980e0   ! Units a0*a0
  d1  	=  0.167127323768e-2  ! Dimensionless
  d2  	=  0.178284243205e1   ! Units 1/a0
  d3  	=  0.176635702255e1   ! Dimensionless
  b 		=  0.203625105759e1   ! Units 1/a0
  ep_kb =  10.997898          ! Units K
  R_ep  =  5.608068    				! Units a0
  sig 	=  4.990672    				! Units a0

  C(6)  =  0.4616213781e6   	! Units K*a0^6
  C(8)  =  0.4460565781e7  		! Units K*a0^8
  C(10) =  0.5803352873e8  		! Units K*a0^10
  C(12) =  0.1031677697e10		! Units K*a0^12
  C(14) =  0.2415716766e11 		! Units K*a0^14
  C(16) =  0.7191492488e12 		! Units K*a0^16
	
	ret1	= 0.0D0
	Z			= 2.0D0

	aa		= 0.88534/SQRT(2.0*Z**0.23)

	IF (R .GT. 0.9875) THEN
		ret0	= A*EXP(a1*R+a2*R**(2.0) + a_1*R**(-1.0) &
		& + a_2*R**(-2.0) + d1*SIN(d2*R+d3))
		DO n=3,8
			ret2	= 0.0D0
			DO k=0,(2*n)
				fact = 0.0D0
				IF ( (k == 0) .OR. (k == 1) ) THEN
					fact = 1.0D0
				ELSE
					fact = 1.0D0
					DO i=2,k
						fact = fact*i
					END DO ! i
				END IF ! if (k)
				ret2 = ret2 + (b*R)**k/fact 
			END DO ! k	
			ret1 = ret1 + (C(2*n)/(R**(2*n)))*(1.0D0-EXP(-b*R)*ret2)
		END DO ! n
		Pot = (ret0 - ret1)/TOK ! [K]  
	ELSE
		atom = 0.1818*EXP(-3.2*(R/aa)) + 0.50990*EXP(-0.9423*(R/aa)) &
		& + 0.2802*EXP(-0.4029*(R/aa)) + 0.02817*EXP(-0.2016*(R/aa))
    Pot	=	(Z*Z/R)*atom*10**1.37503/TOEV
	END IF	


END SUBROUTINE 


