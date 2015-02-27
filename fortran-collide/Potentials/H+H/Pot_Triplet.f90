!
!      H---H Triplet Ground State Potentail
!
          SUBROUTINE pot_HH_triplet(X,VX,DVX,D2VX,IOP)

					 USE physics_constants, ONLY : TOEV

           IMPLICIT NONE
           
           INTEGER :: IOP

           DOUBLE PRECISION :: X, VX, DVX, D2VX
           DOUBLE PRECISION :: XP1, XP2, XM1, XM2
           DOUBLE PRECISION :: VXP1, VXP2, VXM1, VXM2
           DOUBLE PRECISION :: H
           DOUBLE PRECISION :: MYTPOT 

           EXTERNAL MYTPOT
           

           VX = MYTPOT(X)

           IF (IOP .EQ. 0) RETURN

           IF (X .GT. 15.0D+00) H = 0.05D+00
           IF (X .LT. 0.10D+00) H = X * 0.05D+00    ! move 5% 
           IF ( (X .GE. 0.1D+00) .AND. (X .LE. 15.0D+00) ) H = 0.007D+00

           XP1 = X + H
           XP2 = XP1 + H
           XM1 = X - H
           XM2 = XM1 - H

           VXP1 = MYTPOT(XP1)
           VXP2 = MYTPOT(XP2)
           VXM1 = MYTPOT(XM1)
           VXM2 = MYTPOT(XM2)
           
           DVX = ( 8.0D+00 * (VXP1 - VXM1) + VXM2 - VXP2 ) / (12.0D+00 *H)

           D2VX = ( 16.0D+00 * (VXP1 + VXM1) - 30.0D+00 * VX - VXP2 - VXM2 ) / (12.0D+00 * H * H)

           RETURN

          END SUBROUTINE pot_HH_triplet
!
!=========================================================================================================
!

          DOUBLE PRECISION FUNCTION MYTPOT(R)

           IMPLICIT NONE
  
           DOUBLE PRECISION :: R
           DOUBLE PRECISION, PARAMETER :: RSHORT = 0.20000000D+00
           DOUBLE PRECISION :: SHORT
           DOUBLE PRECISION, PARAMETER :: SCALEFACTOR = 1.16D+00   !  DETERMINED BY NUMBERCAL EXPERIMENT

           EXTERNAL SHORT

           IF (R .LT. RSHORT) THEN

              ! USING LENZ-JENSEN Short Range Potential TO GET VERY SHORT RANGE RESULT IF NECESSARY
              MYTPOT = SHORT(R,0) * SCALEFACTOR
           ELSE

              ! USING MJJ's POTENTIA
              CALL CALPOT(R,MYTPOT,1)

           END IF

          END FUNCTION MYTPOT
!
!=========================================================================================================
!

