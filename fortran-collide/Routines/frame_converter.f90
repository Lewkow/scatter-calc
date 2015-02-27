
SUBROUTINE angle_to_lab(theta_cm, theta_lab)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Converts scattering angle from center of mass (CM)
! to lab frame. 
!
! Both the input theta_cm and the output theta_lab
! are in radians
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : M1, M2
	!------------------------
	! Inputs
	!------------------------

	REAL(KIND=8)				:: theta_cm

	!------------------------
	! Outputs
	!------------------------

	REAL(KIND=8)				:: theta_lab

	!------------------------
	! Internal 
	!------------------------

	REAL(KIND=8)				:: ret

	
	ret 			= SIN(theta_cm)/(M1/M2 + COS(theta_cm))
	theta_lab = DATAN(ret) 

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE amp_to_lab(theta, amp_cm, amp_lab)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Converts scattering amplitude from center of mass
! (CM) to lab frame. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : M1, M2

	IMPLICIT NONE

	!------------------------
	! Inputs
	!------------------------

	REAL(KIND=8)					:: amp_cm
	REAL(KIND=8)					:: theta 

	!------------------------
	! Outputs
	!------------------------

	REAL(KIND=8)					:: amp_lab

	!------------------------
	! Internal 
	!------------------------

	REAL(KIND=8)					:: ret1
	REAL(KIND=8)					:: ret2


	ret1 		= 1.0 + 2.0*(M1/M2)*COS(theta) + (M1/M2)*(M1/M2)
	ret1 		= ret1**1.5
	ret2 		= ABS(1.0 + (M1/M2)*COS(theta))
	amp_lab = amp_cm*(ret1/ret2)

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE energy_to_lab(E_cm, E_lab)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Convert CM energy to lab frame energy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : M1, M2

	IMPLICIT NONE

	REAL(KIND=8)	:: E_cm, E_lab

	E_lab = ((M1+M2)/M2)*E_cm

END SUBROUTINE energy_to_lab

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE energy_to_cm(E_lab, E_cm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Convert lab frame energy to CM frame energy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : M1, M2

	IMPLICIT NONE

	REAL(KIND=8)	:: E_cm, E_lab

	E_cm = (M2/(M1+M2))*E_lab

END SUBROUTINE energy_to_cm

