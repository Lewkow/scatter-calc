
!   Lenz-Jensen Short Range Potential
!   D J O'Connor and J P Bierack, Nucl. Inst. Meth. B, 15, 14-19 (1986)


         DOUBLE PRECISION FUNCTION SHORT(R,IOP)

       
          IMPLICIT NONE
          INTEGER :: IOP
          DOUBLE PRECISION :: R, SL, SLF, F, PHI, X
          DOUBLE PRECISION :: T, T2
          DOUBLE PRECISION, PARAMETER :: ALPHA = 0.041D+00
          DOUBLE PRECISION, PARAMETER :: BETA = 0.81D+00
          DOUBLE PRECISION, PARAMETER :: Z1 = 1.0D+00
          DOUBLE PRECISION, PARAMETER :: Z2 = 1.0D+00
          DOUBLE PRECISION, PARAMETER :: A = 0.88534D+00

          F = 1.0D+00

          IF (IOP .EQ. 1) F = ALPHA * (DSQRT(Z1) + DSQRT(Z2)) + BETA 

          SL = DSQRT(Z1) + DSQRT(Z2)
          SL = SL**(2.0D+00/3.0D+00)
          SL = A / SL

          SLF = SL * F
          X = R / SLF
          T = 3.11126D+00 * DSQRT(X)
          T2 = T * T
          PHI = 1.0D+00 + T + 0.3344D+00 * T2 + 0.0485D+00 * T * T2 + 2.647D-03 * T2 * T2 
          PHI = DEXP(-T) * PHI

          SHORT = Z1 * Z2 * PHI / R 

          RETURN
         END




 
