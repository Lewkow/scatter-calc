!
!       H---H SInglet Ground State Potentail
!
          SUBROUTINE pot_HH_singlet(X,VX,DVX,D2VX,IOP)

           USE physics_constants, ONLY : TOEV

           IMPLICIT NONE
           
           INTEGER :: IOP

           DOUBLE PRECISION :: X, VX, DVX, D2VX
           DOUBLE PRECISION :: XP1, XP2, XM1, XM2
           DOUBLE PRECISION :: VXP1, VXP2, VXM1, VXM2
           DOUBLE PRECISION :: H
           DOUBLE PRECISION :: MYSPOT 

           EXTERNAL MYSPOT
           

           VX = MYSPOT(X)

           IF (IOP .EQ. 0) RETURN

           IF (X .GT. 15.0D+00) H = 0.05D+00
           IF (X .LT. 0.10D+00) H = X * 0.05D+00    ! move 5% 
           IF ( (X .GE. 0.1D+00) .AND. (X .LE. 15.0D+00) ) H = 0.007D+00

           XP1 = X + H
           XP2 = XP1 + H
           XM1 = X - H
           XM2 = XM1 - H

           VXP1 = MYSPOT(XP1)
           VXP2 = MYSPOT(XP2)
           VXM1 = MYSPOT(XM1)
           VXM2 = MYSPOT(XM2)
           
           DVX = ( 8.0D+00 * (VXP1 - VXM1) + VXM2 - VXP2 ) / (12.0D+00 *H)

           D2VX = ( 16.0D+00 * (VXP1 + VXM1) - 30.0D+00 * VX - VXP2 - VXM2 ) / (12.0D+00 * H * H)

           RETURN

          END SUBROUTINE pot_HH_singlet 
!
!=========================================================================================================
!

          DOUBLE PRECISION FUNCTION MYSPOT(R)

           IMPLICIT NONE
  
           DOUBLE PRECISION :: R
           DOUBLE PRECISION, PARAMETER :: RSHORT = 0.14415114D+00
           DOUBLE PRECISION, PARAMETER :: RSMALL = 11.99D+00
           DOUBLE PRECISION, PARAMETER :: RLARGE = 20.60D+00
           DOUBLE PRECISION :: SHORT

           EXTERNAL SHORT

           IF (R .LT. RSHORT) THEN

              ! USING LENZ-JENSEN Short Range Potential TO GET VERY SHORT RANGE RESULT IF NECESSARY
              MYSPOT = SHORT(R,0)

           ELSE IF ((R .GT. RSMALL) .AND. (R .LT. RLARGE)) THEN

              ! FIXING THE BUMPER AT 12 BOHR IF NECESSARY
              CALL FIXLONG(R,MYSPOT)

           ELSE

              ! USING MJJ's POTENTIA
              CALL CALPOT(R,MYSPOT,0)

           END IF

          END FUNCTION MYSPOT
!
!=========================================================================================================
!
           SUBROUTINE FIXLONG(R,VR)

           IMPLICIT NONE
           INTEGER, PARAMETER :: NP = 253
           INTEGER :: I, ICALL
           
           DOUBLE PRECISION :: W1(253), W2(253), W3(253) 
           DOUBLE PRECISION :: YP1, YPN
           DOUBLE PRECISION :: R, VR
           
           DATA ICALL /0/
           SAVE ICALL
           SAVE  W1, W2, W3

					OPEN(UNIT=65,FILE="../Potentials/H+H/fort.9",STATUS="OLD",ACTION="READ")	

           IF (ICALL .EQ. 0) THEN
              YP1 = 1.0D+50
              YPN = 1.0D+50
              DO I = 1, NP
                 READ (65,*) W1(I), W2(I)
              END DO
              CALL DERIVE(W1,W2,W3,NP,YP1,YPN)
              ICALL = 1
           END IF
           
           CALL SPLINT(W1,W2,W3,NP,R,VR)
              
           RETURN
           END
!
!=========================================================================================================
!

