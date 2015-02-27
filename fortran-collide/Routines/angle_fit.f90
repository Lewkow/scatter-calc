
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write [energy angle probability] for given 
! collision to file which may be used for other
! types of calculations (ie Monte Carlo) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_angle_probability
	USE collision_inputs, ONLY : ThI, ThF, NTh
	USE physics_constants, ONLY : PI, TOEV
	USE current_energy, ONLY : E

	IMPLICIT NONE

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: PD, Theta_now

	REAL(KIND=8)		:: Theta, DT, DA, Tot
	INTEGER					:: i, j, N_tot, N, N_Angles

	OPEN(UNIT=77, FILE="../Data/AngProb3D.dat", ACCESS="APPEND")

	N_tot 	 = 20000								!! total number of angles to sum up to PI
	DT    	 = PI/N_tot		  				!! angular itteration step
	N_Angles = 180*2								!! total number of angles to sum for	
	DA       = PI/REAL(N_Angles) 		!! Difference in angles printed to file

	ALLOCATE(PD(N_tot),Theta_now(N_tot))

	DO i=1,N_tot
		Theta_now(i) = i*DT
		CALL prob_density(Theta_now(i),PD(i))
	END DO

	WRITE(77,*) E*TOEV, 0.0D0, 0.0D0

	DO i=1,N_Angles
		Theta = i*DA
		N     = INT(CEILING(Theta*REAL(N_tot)/PI))
		tot   = 0.0D0
		DO j=1,N
			tot = tot + DT*PD(j)
		END DO
		WRITE(77,*) E*TOEV, Theta*180.0D0/PI, tot
	END DO		

	CLOSE(77)

	DEALLOCATE(PD,Theta_now)

END SUBROUTINE write_angle_probability

!###########################################################
!###########################################################
!###########################################################

SUBROUTINE prob_den_angle_fitter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fits probability density to an analytic 
! function such that y(RND) = theta
! where RND is a random number generated [0-1]
! and theta is a "random" scattering angle. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE collision_inputs, ONLY : ThI, ThF, NTh
	USE LM, ONLY : LM_x, LM_y, LM_N, JAC_ON
	USE physics_constants, ONLY : PI, TOEV
	USE current_energy, ONLY : E

	IMPLICIT NONE

	REAL(KIND=8)	:: PD, Theta, DT, tot
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: PD_A, p, x_tot, y_tot
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: p_Tot 
	INTEGER				:: i, j, cutoff, ci
	INTEGER				:: exp_sum_on, gauss_on, gauss_sum_on, pw_gauss_sum_on

	JAC_ON          = 0	

	exp_sum_on      = 0
	gauss_sum_on    = 1 
	gauss_on				= 0
	pw_gauss_sum_on = 0 

	LM_N = NTh
	DT 	 = (ThF - ThI)/DBLE(NTh-1)	 

	cutoff = 0

	ALLOCATE( PD_A(NTh), x_tot(NTh), y_tot(NTh) )
	
	DO i=1,NTh
		y_tot(i) = ThI + (i-1)*DT
		CALL prob_density(y_tot(i),PD_A(i))
	END DO ! i

	DO i=1,NTh
		Tot = 0.0D0
		DO j=1,i
			Tot = Tot + DT*PD_A(j)
		END DO!j
		x_tot(i) = Tot
		IF ((pw_gauss_sum_on .EQ. 1).AND.(cutoff .EQ. 0).AND.(x_tot(i) .GT. 0.8D0)) cutoff = i 
	END DO!i 	

	WRITE(*,*) "Cutoff Index: ", cutoff, " Angle: ", y_tot(cutoff)

	OPEN(UNIT=666, FILE="../Data/AngleFits.dat", ACCESS="APPEND")

	!!##############################################
	!! Change which function to fit to here
	IF ( exp_sum_on .EQ. 1) THEN
		ALLOCATE(p(4),LM_x(NTh),LM_y(NTh))
		LM_x(:) = x_tot(:)
		LM_y(:) = y_tot(:)
		CALL exponential_sum_fit(NTh,p)	
		WRITE(666,*) E*TOEV, p	
		DEALLOCATE(p, LM_x, LM_y)
	END IF

	IF ( gauss_on .EQ. 1 ) THEN
		ALLOCATE(p(3),LM_x(NTh),LM_y(NTh))
		LM_x(:) = x_tot(:)
		LM_y(:) = y_tot(:)
		CALL gaussian_fit(NTh,p)
		WRITE(666,*) E*TOEV, p
		DEALLOCATE(p, LM_x, LM_y)
	END IF

	IF ( gauss_sum_on .EQ. 1 ) THEN
!		ALLOCATE(p(9),LM_x(NTh),LM_y(NTh))
		ALLOCATE(p(10),LM_x(NTh),LM_y(NTh))
		LM_x(:) = x_tot(:)
		LM_y(:) = y_tot(:)
		CALL gaussian_sum_fit(NTh,p)
		WRITE(666,*) E*TOEV, p
		DEALLOCATE(p, LM_x, LM_y)
	END IF

	IF ( pw_gauss_sum_on .EQ. 1) THEN
		ALLOCATE(p(10),p_tot(10,2))
		ALLOCATE(LM_x(cutoff),LM_y(cutoff))
		LM_x(:) = x_tot(1:cutoff)
		LM_y(:) = y_tot(1:cutoff)
		CALL gaussian_sum_fit(cutoff,p)
		p_tot(:,1) = p
		p(:) = 1.0D0
		DEALLOCATE(LM_x,LM_y)
		ci = NTh - cutoff	
		ALLOCATE(LM_x(ci),LM_y(ci))
		LM_x(:) = x_tot(cutoff+1:NTh)
		LM_y(:) = y_tot(cutoff+1:NTh)
		CALL gaussian_sum_fit(ci,p)
		p_tot(:,2) = p
		WRITE(666,*) E*TOEV, p_tot(:,1), p_tot(:,2)
		DEALLOCATE(p, p_tot, LM_x, LM_y)	
	END IF
	!!##############################################

	DEALLOCATE(PD_A, x_tot, y_tot)
	CLOSE(666)

END SUBROUTINE prob_den_angle_fitter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE exponential_sum_fit(m,p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fit the 4-parameter logistic curve:
!
!     Y = A1 * exp{B1*X} + A2 * exp{B2*x} 
!
! by unweighted least squares.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE LM
IMPLICIT NONE

INTERFACE
  SUBROUTINE exponential_sum_fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE exponential_sum_fcn
END INTERFACE

REAL (dp), ALLOCATABLE  :: fvec(:), fjac(:,:)
INTEGER, PARAMETER      :: n = 4                 ! The number of parameters
INTEGER                 :: info, iostatus, ipvt(n), m
REAL (dp)               :: p(n), tol = 1.0D-03

ALLOCATE( fvec(m), fjac(m,n) )

! Set starting values for parameters.
! A1 = p(1), B1 = p(2), A2 = p(3), B2 = p(4)

p(1) = 1.0D0
p(3) = 1.0D0
p(2) = 1.0D0
p(4) = 1.0D0

CALL lmder1(exponential_sum_fcn, m, n, p, fvec, fjac, tol, info, ipvt)


SELECT CASE (info)
  CASE (:-1)
    WRITE(*, *) 'Users FCN returned INFO = ', -info
  CASE (0)
    WRITE(*, *) 'Improper values for input parameters'
  CASE (4)
    WRITE(*, *) 'Residuals orthogonal to the Jacobian'
    WRITE(*, *) 'There may be an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (5)
    WRITE(*, *) 'Too many calls to FCN'
    WRITE(*, *) 'Either slow convergence, or an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (6:7)
    WRITE(*, *) 'TOL was set too small'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
END SELECT

END SUBROUTINE exponential_sum_fit

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE exponential_sum_fcn(m, n, p, fvec, fjac, iflag)
! Calculate either residuals or the Jacobian matrix.
! A = p(1), B = p(2), C = p(3), D = p(4)
! m = no. of cases, n = no. of parameters (4)

USE LM, ONLY : LM_x, LM_y 
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: p(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag

! Local variables

INTEGER               :: i

IF (iflag == 1) THEN
 fvec = LM_y(1:m) - p(1) * EXP(p(2)*LM_x(1:m)) - p(3) * EXP(p(4)*LM_x(1:m))
ELSE IF (iflag == 2) THEN
  DO i = 1, m
  	fjac(i,1) = -EXP(p(2)*LM_x(i))
    fjac(i,2) = -p(1) * LM_x(i) * EXP(p(2)*LM_x(i))
    fjac(i,3) = -EXP(p(4)*LM_x(i))
    fjac(i,4) = -p(3) * LM_x(i) * EXP(p(4)*LM_x(i)) 
  END DO
END IF

RETURN
END SUBROUTINE exponential_sum_fcn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE gaussian_sum_fit(m,p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fit the 9-parameter logistic curve:
!
!     Y = sum_{i=1}^{i=3} ai*exp(-( (x-bi)/ci )^2) 
!
! by unweighted least squares.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE LM
IMPLICIT NONE

INTERFACE
  SUBROUTINE gaussian_sum_jfcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE gaussian_sum_jfcn

  SUBROUTINE gaussian_sum_fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE gaussian_sum_fcn
END INTERFACE

REAL (dp), ALLOCATABLE  :: fvec(:), fjac(:,:)
!INTEGER, PARAMETER      :: n = 9                 ! The number of parameters
INTEGER, PARAMETER      :: n = 10                 ! The number of parameters
INTEGER                 :: info, iostatus, ipvt(n), m
REAL (dp)               :: p(n), tol = 1.0E-04_dp

ALLOCATE( fvec(m), fjac(m,n) )

! Set starting values for parameters.

p(:) = 1.0D0

IF (JAC_ON .EQ. 1) THEN
	CALL lmder1(gaussian_sum_jfcn, m, n, p, fvec, fjac, tol, info, ipvt)
ELSE
	CALL lmdif1(gaussian_sum_fcn, m, n, p, fvec, tol, info, ipvt)
END IF

SELECT CASE (info)
  CASE (:-1)
    WRITE(*, *) 'Users FCN returned INFO = ', -info
  CASE (0)
    WRITE(*, *) 'Improper values for input parameters'
  CASE (4)
    WRITE(*, *) 'Residuals orthogonal to the Jacobian'
    WRITE(*, *) 'There may be an error in FCN'
    WRITE(*, *) ' Final values p: ', p
  CASE (5)
    WRITE(*, *) 'Too many calls to FCN'
    WRITE(*, *) 'Either slow convergence, or an error in FCN'
    WRITE(*, *) ' Final values of p: ', p
  CASE (6:7)
    WRITE(*, *) 'TOL was set too small'
    WRITE(*, *) ' Final values of p: ', p
END SELECT

DEALLOCATE( fvec, fjac )

END SUBROUTINE gaussian_sum_fit

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE gaussian_sum_fcn(m, n, p, fvec, iflag)
! Calculate residuals 
! m = no. of cases, n = no. of parameters (4)

USE LM, ONLY : LM_x, LM_y 
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: p(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
INTEGER, INTENT(IN OUT)    :: iflag

! Local variables
REAL(KIND=8)						:: exp1, exp2, exp3, f1, f2, f3
REAL(KIND=8),PARAMETER	:: two = 2.0D0, one = 1.0D0
INTEGER               	:: i

!	fvec =  p(1) * EXP( - ( (LM_x(1:m) - p(2))/p(3) )**2 ) &
!       &  P(4) * EXP( - ( (LM_x(1:m) - p(5))/p(6) )**2 ) &
!       &  p(7) * EXP( - ( (LM_x(1:m) - P(8))/p(9) )**2 ) + p(10)
fvec = LM_y(1:m) - p(1) * EXP( - ( (LM_x(1:m) - p(2))/p(3) )**2 ) &
    	         & - P(4) * EXP( - ( (LM_x(1:m) - p(5))/p(6) )**2 ) &
               & - p(7) * EXP( - ( (LM_x(1:m) - P(8))/p(9) )**2 ) - p(10)
RETURN
END SUBROUTINE gaussian_sum_fcn

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE gaussian_sum_jfcn(m, n, p, fvec, fjac, iflag)
! Calculate either residuals or the Jacobian matrix.
! m = no. of cases, n = no. of parameters (4)

USE LM, ONLY : LM_x, LM_y 
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: p(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag

! Local variables
REAL(KIND=8)						:: exp1, exp2, exp3, f1, f2, f3
REAL(KIND=8),PARAMETER	:: two = 2.0D0, one = 1.0D0
INTEGER               	:: i

IF (iflag == 1) THEN
 fvec = LM_y(1:m) - p(1) * EXP( - ( (LM_x(1:m) - p(2))/p(3) )**2 ) &
								& - P(4) * EXP( - ( (LM_x(1:m) - p(5))/p(6) )**2 ) &
                & - p(7) * EXP( - ( (LM_x(1:m) - P(8))/p(9) )**2 ) 
ELSE IF (iflag == 2) THEN
  DO i = 1, m
		f1  			= (LM_x(i) - p(2))/p(3) 
		f2  			= (LM_x(i) - p(5))/p(6) 
		f3  			= (LM_x(i) - p(8))/p(9) 
		f1        = f1*f1
		f2        = f2*f2
		f3        = f3*f3
		exp1 			= EXP( -f1 )
		exp2 			= EXP( -f2 )
		exp3 			= EXP( -f3 )
  	fjac(i,1) = - exp1 
    fjac(i,2) = - (two*p(1)/(p(3)**2))*(LM_x(i)-p(2))*exp1 
    fjac(i,3) = - (two*p(1)/(p(3)**3))*((LM_x(i)-p(2))**2)*exp1
    fjac(i,4) = - exp2 
    fjac(i,5) = - (two*p(4)/(p(6)**2))*(LM_x(i)-p(5))*exp2
    fjac(i,6) = - (two*p(4)/(p(6)**3))*((LM_x(i)-p(5))**2)*exp2
    fjac(i,7) = - exp3
    fjac(i,8) = - (two*p(7)/(p(9)**2))*(LM_x(i)-p(8))*exp3
    fjac(i,9) = - (two*p(7)/(p(9)**3))*((LM_x(i)-p(8))**2)*exp3
!		fjac(i,10)= - one
  END DO
END IF

RETURN
END SUBROUTINE gaussian_sum_jfcn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE gaussian_fit(m,p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fit the 3-parameter logistic curve:
!
!     Y = a*exp(-( (x-b)/c )^2) 
!
! by unweighted least squares.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE LM
IMPLICIT NONE

INTERFACE
  SUBROUTINE gaussian_fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE gaussian_fcn
END INTERFACE

REAL (dp), ALLOCATABLE  :: fvec(:), fjac(:,:)
INTEGER, PARAMETER      :: n = 3                 ! The number of parameters
INTEGER                 :: info, iostatus, ipvt(n), m
REAL (dp)               :: p(n), tol = 1.0D-04

ALLOCATE( fvec(m), fjac(m,n) )

! Set starting values for parameters.

p(:) = 1.0D0

CALL lmder1(gaussian_fcn, m, n, p, fvec, fjac, tol, info, ipvt)


SELECT CASE (info)
  CASE (:-1)
    WRITE(*, *) 'Users FCN returned INFO = ', -info
  CASE (0)
    WRITE(*, *) 'Improper values for input parameters'
  CASE (4)
    WRITE(*, *) 'Residuals orthogonal to the Jacobian'
    WRITE(*, *) 'There may be an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (5)
    WRITE(*, *) 'Too many calls to FCN'
    WRITE(*, *) 'Either slow convergence, or an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (6:7)
    WRITE(*, *) 'TOL was set too small'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
END SELECT

END SUBROUTINE gaussian_fit

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE gaussian_fcn(m, n, p, fvec, fjac, iflag)
! Calculate either residuals or the Jacobian matrix.
! m = no. of cases, n = no. of parameters (4)

USE LM, ONLY : LM_x, LM_y 
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: p(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag

! Local variables
REAL(KIND=8)						:: exp1
REAL(KIND=8),PARAMETER	:: two = 2.0D0, one = 1.0D0
INTEGER               	:: i

IF (iflag == 1) THEN
 fvec = LM_y(1:m) - p(1) * EXP( - ( (LM_x(1:m) - p(2))/p(3) )**2 ) 
ELSE IF (iflag == 2) THEN
  DO i = 1, m
		exp1 			= EXP( - ( (LM_x(i) - p(2))/p(3) )**2 )
  	fjac(i,1) = - exp1 
    fjac(i,2) = - (two*p(1)/p(3)**2)*(LM_x(i)-p(2))*exp1 
    fjac(i,3) = - (two*p(1)/p(3)**3)*((LM_x(i)-p(2))**2)*exp1
  END DO
END IF

RETURN
END SUBROUTINE gaussian_fcn

