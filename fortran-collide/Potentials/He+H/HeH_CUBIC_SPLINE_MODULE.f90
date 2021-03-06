
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module containing all routines and data structures
! needed to compute cubic spline values for data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE HeH_cubic_spline

	IMPLICIT NONE

	INTEGER																	:: Call_spline	! if 1, don't call spline
	INTEGER																	:: N_points

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: y2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: HeH_pot			! potential energy points 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: HeH_sep 			! potential seperation distance

	CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE HeH_array_setup
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	! Sets up the arrays for the He+H splines if
	! not already setup	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	IMPLICIT NONE

	INTEGER		:: i

	! File containing data points for spline
	OPEN(UNIT=66, FILE="../Potentials/He+H/HeH_Points.dat", STATUS="OLD", ACTION="READ")

	! get number of points in file as first line number
	READ(66,*) N_points	

	! allocate arrays
	ALLOCATE( y2(N_points), HeH_pot(N_points), HeH_sep(N_points) )

	DO i=1,N_Points
		READ(66,*) HeH_sep(i), HeH_pot(i)
	END DO	

	END SUBROUTINE HeH_array_setup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE HeH_cspline(x, y)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	! cspline acts as a front end to the user to call 
	! the cubic spline routines. If the spline
	! routine has not yet been called it needs to be 
	! called before calling splinet.	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: x								! x value to find spline to 	

	!! Outputs
	REAL(KIND=8)	:: y								! y value from spline	

	!! Internal
	REAL(KIND=8)	:: yp1, ypn					! end point 2nd derivatives for splines

	yp1 = 1.0D30
	ypn = 1.0D30

	IF (Call_spline /= 1) THEN
		CALL HeH_array_setup
		CALL spline( HeH_sep, HeH_pot, N_points, yp1, ypn, y2 )
		Call_spline = 1
	END IF

	CALL splint( HeH_sep, HeH_pot, y2, N_points, x, y )	
	
	END SUBROUTINE HeH_cspline

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE spline(x,y,n,yp1,ypn,y2)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Numerical Methods routine 
	!
	! Given arrays x(1:n) and y(1:n) containing a tabulated function
	! yi = f(xi), with x1 < x2 < ... < xn, and given values yp1 and ypn 
	! for the first derivative of the interpolating function at points
	! 1 and n, this routine returns an array y2(1:n) of length n which
	! contains the second derivatives of the interpolating function
	! at the tabulated points xi. If yp1 and/or ypn are equal to 
	! 1 x 10^30 or larger, the routine is signaled to set the
	! corresponding boundary condition for a natural spline, with zero
	! second derivative on that boundary. 
	! Parameter: NMAX is teh largest anticipated value of n
	!
	! It is important to understand that the program spline is called
	! only once to process an entire tabualted function in arrays
	! xi and yi. Once this has been done, values of the interpolated 
	! function for any value of x are obtained by calls to the 
	! routine splint (spline interpolation)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER                   :: n
  INTEGER                   :: NMAX
  REAL(KIND=8)              :: yp1
  REAL(KIND=8)              :: ypn
  REAL(KIND=8)              :: x(n)
  REAL(KIND=8)              :: y(n)
  REAL(KIND=8)              :: y2(n)

  PARAMETER (NMAX = 10000)

  INTEGER                   :: i
  INTEGER                   :: k
  REAL(KIND=8)              :: p
  REAL(KIND=8)              :: qn
  REAL(KIND=8)              :: sig
  REAL(KIND=8)              :: un
  REAL(KIND=8)              :: u(NMAX)

  IF (yp1 .GT. .99e30) THEN
    y2(1) = 0
    u(1)  = 0
  ELSE
    y2(1) = -0.5
    u(1)  = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  END IF

DO i=2,n-1
    sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p     = sig*y2(i-1)+2.
    y2(i) = (sig-1.)/p
    u(i)  = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  END DO

  IF (ypn .GT. .99e30) THEN
    qn = 0.
    un = 0.
  ELSE
    qn = 0.5
    un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  END IF

  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

  DO k=n-1,1,-1
    y2(k) = y2(k)*y2(k+1)+u(k)
  END DO

  RETURN

	END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE splint(xa,ya,y2a,n,x,y)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Numerical Methods routine 
	!
	! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate
	! a function, with the xai's in order, and given the array y2a(1:n)
	! which is the output from the spline routine, and given a value
	! of x, this routine returns a cubic-spline interpolated value y. 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER                           :: n
  REAL(KIND=8)                      :: x
  REAL(KIND=8)                      :: y
  REAL(KIND=8)                      :: xa(n)
  REAL(KIND=8)                      :: ya(n)
  REAL(KIND=8)                      :: y2a(n)

  INTEGER                           :: k
  INTEGER                           :: khi
  INTEGER                           :: klo
  REAL(KIND=8)                      :: a
  REAL(KIND=8)                      :: b
  REAL(KIND=8)                      :: h

  klo = 1
  khi = n

1 IF (khi-klo .GT. 1) THEN
    k = (khi+klo)/2
    IF (xa(k) .GT. x) THEN
      khi = k
    ELSE
      klo = k
    END IF
  GOTO 1
  END IF

  h = xa(khi) - xa(klo)
  IF (h .EQ. 0.) THEN
    WRITE(*,*) "bad xa input in splint"
  END IF

  a = (xa(khi)-x)/h
  b = (x-xa(klo))/h
  y = a*ya(klo)+b*ya(khi)+((a**3-1)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

  RETURN

	END SUBROUTINE


END MODULE HeH_cubic_spline

