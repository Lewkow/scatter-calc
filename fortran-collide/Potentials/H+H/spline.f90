
!=============================================================

      SUBROUTINE DERIVE(X,Y,YDER,N,YP1,YPN)

! FOR CUBIC SPLINE FITTING PROGRAM

! GIVEN TABLE [X,Y] OF N PAIRS AND FIRST DERIVATIVES
! AT FIRST AND LAST POINTS YP1 AND YPN, THIS SUBROUTINE RETURNS
! YDER, THE SECOND DERIVATIVES AT N POINTS. IF YP1 AND/OR YPN ARE
! UNKNOWN SET TO 1.0E30 OR LARGER AND THEY ARE CALCULATED
! BY THIS SUBROUTINE. U IS A SCRATCH ARRAY.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),Y(N),YDER(N)
      DIMENSION U(N)
      BIG=0.99D30
      IF (N.LT.4) THEN
      WRITE(6,10)
   10 FORMAT('ERROR: MORE THAN 3 DATA POINTS ARE NEEDED')
      STOP
      END IF


! CALCULATE SPLINE COEFFCIENTS 

      XF=X(3)-X(2)
      XB=X(2)-X(1)
      YF=Y(3)-Y(2)
      YB=Y(2)-Y(1)

      IF (YP1.GT.BIG) THEN
        P=2.0D0*XF+XB
        YDER(2)=(XB-XF)/P
        U(2)=6.0D0*(YF-YB*XF/XB)/((XF+XB)*P)
      ELSE
        YDER(1)=-0.5D0
        U(1)=(3.0D0/XB)*(YB/XB-YP1)
        P=2.0D0*(XF+XB)+XB*YDER(1)
        YDER(2)=-XF/P
        U(2)=(6.0D0*(YF/XF-YB/XB)-XB*U(1))/P
      ENDIF

      DO 11 I=3,N-1
        XF=X(I+1)-X(I)
        XB=X(I)-X(I-1)
        YF=Y(I+1)-Y(I)
        YB=Y(I)-Y(I-1)
        P=2.0D0*(XF+XB)+XB*YDER(I-1)
        YDER(I)=-XF/P
        U(I)=(6.0D0*(YF/XF-YB/XB)-XB*U(I-1))/P
11    CONTINUE

      IF (YPN.GT.BIG) THEN
        P=2.0D0*XF+XB
        YDER(N)=(6.0D0*(YF-YB*XF/XB)/(XF+XB)-P*U(N-1))/(XF-XB+P*YDER(N-1))
      ELSE
        YDER(N)=(6.0D0*(YPN-YF/XF)/XF-U(N-1))/(2.0D0+YDER(N-1))
      ENDIF

      DO 12 K=N-1,2,-1
        YDER(K)=YDER(K)*YDER(K+1)+U(K)
12    CONTINUE

      IF (YP1.GT.BIG) THEN
        YDER(1)=YDER(2)-(X(2)-X(1))*(YDER(3)-YDER(2))/(X(3)-X(2))
        YP1=(Y(2)-Y(1))/(X(2)-X(1)) &
     &     -(X(2)-X(1))*(YDER(1)/3.0D0+YDER(2)/6.0D0)
      ELSE
        YDER(1)=YDER(1)*YDER(2)+U(1)
      END IF

      IF (YPN.GT.BIG) THEN
        YPN=YF/XF+XF*(YDER(N-1)/6.0D0+YDER(N)/3.0D0)
      END IF

      RETURN
      END

!=============================================================

      SUBROUTINE SPLINT(X,Y,YDER2,N,XVAL,YVAL)

! CUBIC SPLINE FITTING PROGRAM

! GIVEN TABLE [X,Y] OF N PAIRS AND 2ND DERIVATIVES YDER2,
! THIS SUBROUTINE RETURNS THE FUNCTION VALUE YVAL FOR ARGUMENT XVAL

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),YDER2(N)
      KLO=1
      KHI=N
    1 IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF (X(K).GT.XVAL) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=X(KHI)-X(KLO)
      A=(X(KHI)-XVAL)/H
      B=(XVAL-X(KLO))/H
      YVAL=A*Y(KLO)+B*Y(KHI)+((A**3-A)*YDER2(KLO)+(B**3-B)*YDER2(KHI))*(H**2)/6.0D0
      RETURN
      END

!=============================================================

