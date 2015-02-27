
SUBROUTINE pot_HeH(R,pot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! HeH potential using ab initio data points given
! by Peng with cubic spline interpolation. For core,
! R < 0.7 a0, the universal potential is fit to 
! the ab initio data points. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE HeH_cubic_spline

	IMPLICIT NONE

	!--------------------------
	! Input	
	!--------------------------

	REAL(KIND=8)										:: R

	!--------------------------
	! Output	
	!--------------------------

	REAL(KIND=8)										:: Pot

	!--------------------------
	! Internal	
	!--------------------------

	REAL(KIND=8)										:: Z1			! H  charge
	REAL(KIND=8)										:: Z2			! He charge
	REAL(KIND=8)										:: Last		! Last Pot2 value
	REAL(KIND=8)										:: xi		
	REAL(KIND=8)										:: fit		
	REAL(KIND=8)										:: aa		
	REAL(KIND=8)										:: R0		
	REAL(KIND=8)										:: x1	    ! Fitting constant for core	
	REAL(KIND=8)										:: x2			! Fitting constant for core		
	
	INTEGER													:: i			! counter
	INTEGER													:: j			! counter

	Z1 		= 1.0
	Z2		= 2.0
	aa		= 0.88534/SQRT(Z1**(0.23)+Z2**(0.23))

	! Spline
	x1    = 0.6038
	x2		= -0.00897
	R0 		= 0.7

	IF (R .LE. R0) THEN
		fit = 0.1818*EXP(-3.2*R/aa)+0.5099*EXP(-0.9423*R/aa) &
		& +0.2802*EXP(-0.4029*R/aa)+0.02817*EXP(-0.2016*R/aa)
		Pot = (Z1*Z2/R)*fit*x1 + x2
	ELSE
		CALL HeH_cspline(R,Pot)
	END IF

END SUBROUTINE


