      SUBROUTINE NOCROP (NLXM, TRWL, GAI)
      IMPLICIT NONE

      INTEGER NLXM, I1
      REAL GAI, TRWL(NLXM)
      SAVE

      DO I1=1,NLXM
        TRWL(I1) = 0.
      END DO

      GAI   = 0.

      RETURN
      END
