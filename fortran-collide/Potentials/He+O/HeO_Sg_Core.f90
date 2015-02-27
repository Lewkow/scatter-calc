Subroutine pot_HeO_Sg(R,pot) 
	USE physics_constants, ONLY : TOEV
! He + O Sigma State Potential at ROHF-UCCSD(ALL-E)/CBS aug-cc-pCVQZ, aug-cc-pCV5Z and aug-cc-pCV6Z
!
!     Input :
!           R in bohr
!     Ouput :
!           pot potential energy in atomic unit
!
!     Function form:
!          Y(x) = Sum{a(i)*x^i}*exp(-alfa*(x-Beta))+0.5*(1+tanh(gama-lamda*x))*Sum{C(2n+6}/X^(2n+6)} n=0,1,2,3,4, i=0,1,...,8
!
!           Factorized to Y(x) = T1*T2+T3*T4
!
!			The core potential has the form of the "universal" screened coulomb potential

	Implicit None
	Integer 					:: Ipower
	Double Precision :: X, FX, DFX, D2FX2
	Double Precision :: R, RR, RROld, R4, R6, R8, R10, R12, R14
	Double Precision :: R6D, R8D, R10D, R12D, R14D
	Double Precision :: Alpha, Beta, gama, Lamda, pot
	Double Precision :: T1, T2, T3, T31, T4, T1D, T2D, T3D, T4D
	Double Precision :: T1DD, T2DD, T3DD, T4DD
	Double Precision :: Coefficient(14)
	DOUBLE PRECISION :: R0, A, B, z1, z2, aa

	Double Precision, Parameter :: Tocm = 219474.631D0

	Data Alpha, Beta /3.56427785004271D+00, 1.15813233432006D+00/
 
	Data gama, Lamda /-8.98488543384186D+00, -1.62100046108703D+00/

	Data Coefficient /-200848.300235467D+00, 363708.167958770D+00, &
	&                    13.6184841600187D+00, 170639.603977519D+00, &
	&                   -80518.9683416413D+00, 47506.0289675006D+00, &
	&                   -2323.44460322074D+00,-1034.83460461637D+00, &
	&                    463.315376221698D+00, -1117125.874D+00,     &
	&                    21861184.8879551D+00, -8978455288.89173D+00,&
	&                    582868156139.390D+00, -10582951626631.9D+00/

	R0 = 2.04433
	RR = 1.0d0

	IF ( R .LE. R0 ) THEN
		A  = 49.99D0
		B  = -1.130801D0
		z1 = 2.0
		z2 = 8.0
		aa = .88534/(z1**0.23+z2**0.23)
		pot = (A*z1*z2/R)*(0.1818*exp(-3.2*R/aa)+0.5099*exp(-0.9423*R/aa)&
		&+0.2802*exp(-0.4029*R/aa)+0.02817*exp(-0.2016*R/aa))
		pot = pot/TOEV + B/TOEV
	ELSE
		T1   = Dexp(-Alpha * (R - Beta))
		T31  = Dtanh(Gama - Lamda * R)
		T3   = 0.5d0 * (1.0d0 + T31)
		T2   = 0.d0
		T2D  = 0.d0
		T2DD = 0.d0
		Do Ipower = 0, 8
			If (Ipower.EQ.4) R4 = RR
			If (Ipower.EQ.6) R6 = RR
			If (Ipower.EQ.8) R8 = RR
			T2 = T2 + RR * Coefficient(Ipower+1)
			RROld = RR
			RR = R * RR
		Enddo
		R10 = R6 * R4
		R12 = R8 * R4
		R14 = R6 * R8
		T4 = Coefficient(10)/R6 + Coefficient(11)/R8 + Coefficient(12)/R10 + Coefficient(13)/R12 + Coefficient(14)/R14
		FX = T1 * T2 + T3 * T4
		FX = FX / tocm
		pot = FX      
	END IF
 
	Return
	END SUBROUTINE
