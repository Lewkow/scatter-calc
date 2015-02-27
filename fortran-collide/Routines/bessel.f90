
SUBROUTINE AJ(X, J, ret)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes X*A_{J}(X) where A_{J}(X) is the Jth
! spherical Bessel function of the first kind
! with an argument of X. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE

	!--------------------------
	!	Inputs	
	!--------------------------
	
	REAL(KIND=8)						:: X
	INTEGER									:: J

	!--------------------------
	!	Outputs	
	!--------------------------

	REAL(KIND=8)						:: ret 

	!--------------------------
	!	Local	
	!--------------------------

	INTEGER									:: IAccuracy
	INTEGER									:: JBig
	INTEGER									:: N  

	REAL(KIND=8)						:: Big
	REAL(KIND=8)						:: FM1 
	REAL(KIND=8)						:: FM2 
	REAL(KIND=8)						:: Temp 
	REAL(KIND=8)						:: FP1 
	REAL(KIND=8)						:: F 
	REAL(KIND=8)						:: FP2 

	IAccuracy = 20
	Big				= 1.0D20
  JBig      = J + INT( SQRT( DBLE(J*IAccuracy) ) )

	IF (J == 0) THEN
		ret = SIN(X)
 	ELSE IF (J == 1) THEN
		ret = SIN(X)/X - COS(X)
	ELSE
		IF (X .GT. (DBLE(J)+0.5D0)) THEN
			! Use upwards recurrance	
!			WRITE(*,*) "Using Upwards recurrance for L", J
			FM1 = SIN(X)/X - COS(X)
			FM2 = SIN(X)
			DO N=2,J
				F   = (2.0D0*N-1.0D0)*FM1/X - FM2
				FM2 = FM1 
				FM1 = F
			END DO
      ret = F
		ELSE
			! Use downwards recurrance	
!			WRITE(*,*) "Using Downwards recurrance for L", J
			FP1 = 1.0D-20
			FP2 = 1.0D-20	
			DO N=JBig,0,-1
        F = (2.0D0*N+3.0D0)*FP1/X - FP2

				IF (N == J) THEN
					Temp = F
				END IF

				IF (ABS(F) .GT. Big) THEN
					FP2 = FP2/Big
					FP1 = FP1/Big
					F   = F/Big

					IF (N .LT. J) THEN
						Temp = Temp/Big
					END IF ! N<J
		
				END IF ! |F| > Big
	
				FP2 = FP1
				FP1 = F

			END DO ! N
			ret = Temp*SIN(X)/F
		END IF ! X > J+0.5
	END IF

END SUBROUTINE AJ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AN(X, J, ret)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes X*AN_{J}(X) where AN_{J}(X) is the Jth
! spherical Bessel function of the second kind
! (Neuman) with an argument of X. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!--------------------------
	!	Inputs	
	!--------------------------
	
	REAL(KIND=8)						:: X
	INTEGER									:: J

	!--------------------------
	!	Outputs	
	!--------------------------

	REAL(KIND=8)						:: ret 

	!--------------------------
	!	Local	
	!--------------------------

	INTEGER									:: N
	REAL(KIND=8)						:: F
	REAL(KIND=8)						:: FM1 
	REAL(KIND=8)						:: FM2

	IF (J == 0) THEN
		ret = -COS(X)
	ELSE IF (J == 1) THEN
		ret = -COS(X)/X - SIN(X)
	ELSE
		FM2 = -COS(X)
		FM1 = -COS(X)/X - SIN(X)
		DO N=2,J
			F   = (2.0*N-1.0)*FM1/X - FM2
			FM2 = FM1
			FM1 = F
		END DO ! N
		ret = F
	END IF

END SUBROUTINE AN	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AJP(X, J, ret)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes d(X*A_{J}(X))/d(X) where A_{J}(X) is the
! Jth spherical Bessel function of the first kind
! with an argument of X. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	IMPLICIT NONE

	!--------------------------
	!	Inputs	
	!--------------------------
	
	REAL(KIND=8)						:: X
	INTEGER									:: J

	!--------------------------
	!	Outputs	
	!--------------------------

	REAL(KIND=8)						:: ret 

	!--------------------------
	!	Locals	
	!--------------------------

	REAL(KIND=8)						:: past
	REAL(KIND=8)						:: pres
	REAL(KIND=8)						:: futr 

	IF (J == 0) THEN
		ret = COS(X)
	ELSE
		CALL AJ(X,J    ,pres)
		CALL AJ(X,(J+1),futr)
		CALL AJ(X,(J-1),past)	
		ret = ( pres/X - (futr - past) )/2.0D0
	END IF

END SUBROUTINE AJP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ANP(X, J, ret)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes d(X*AN_{J}(X))/d(X) where AN_{J}(X) is 
! the Jth spherical Bessel function of the second 
! kind (Neuman) with an argument of X. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!--------------------------
	!	Inputs	
	!--------------------------
	
	REAL(KIND=8)						:: X
	INTEGER									:: J

	!--------------------------
	!	Outputs	
	!--------------------------

	REAL(KIND=8)						:: ret 

	!--------------------------
	!	Locals	
	!--------------------------

	REAL(KIND=8)						:: past
	REAL(KIND=8)						:: pres
	REAL(KIND=8)						:: futr 
	
	IF (J == 0) THEN
		ret = SIN(X)
	ELSE 
		CALL AN(X,J  ,pres)
		CALL AN(X,J+1,futr)
		CALL AN(X,J-1,past)
		ret = ( pres/X - (futr - past) )/2.0D0
	END IF

END SUBROUTINE ANP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE bessj0(X,ans)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computes Bessel function of the first kind of
! order zero. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!--------------------------
	!	Inputs	
	!--------------------------

	REAL(KIND=8)				:: X
	
	!--------------------------
	!	Inputs	
	!--------------------------

	REAL(KIND=8)				:: ans 

	!--------------------------
	!	Internal	
	!--------------------------
	
	REAL(KIND=8)				:: ax
	REAL(KIND=8)				:: z 
	REAL(KIND=8)				:: xx 
	REAL(KIND=8)				:: y 
	REAL(KIND=8)				:: ans1 
	REAL(KIND=8)				:: ans2

	ax = X

  IF (ABS(ax) < 8.0) THEN 
    y=X*X
    ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))))
    ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))))
    ans  = ans1/ans2
 	ELSE     
		z    = 8.0/ax
    y    = z*z
    xx   = ax-0.785398164
    ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)))
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)))
    ans  = SQRT(0.636619772/ax)*(COS(xx)*ans1-z*SIN(xx)*ans2)
	END IF

END SUBROUTINE bessj0

 
