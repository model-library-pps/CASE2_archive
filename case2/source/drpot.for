*----------------------------------------------------------------------*
* SUBROUTINE DRPOT                                                     *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Author : ??                                                          *
* Date   : ??                                                          *
*                                                                      *
* Purpose: This is a simple water balance module that keeps the water  *
*          content of the soil layers at field capacity                *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  determines action of routine                  -     C,I  *
* NLXM    I4  no. of layers as declared in calling program  -     IN   *
* NL      I4  number of layers specified in input file      -      O   *
* TKLX    R4  thickness of soil compartments                m    IN,O  *
* ZRTMS   R4  maximum rooting depth as soil characteristic  -      O   *
* WCADX   R4  volumetric water content airdry               -    IN,O  *
* WCWPX   R4  volumetric water content at wilting point     -    IN,O  *
* WCFCX   R4  volumetric water content at field capacity    -    IN,O  *
* WCSTX   R4  volumetric water content at saturation        -    IN,O  *
* WCLQT   R4  volumetric soil water content per layer       -      O   *
*                                                                      *
* Fatal error checks: none                                             *
* Warnings          : none                                             *
* Subprograms called: none                                             *
* File usage        : none                                             *
*----------------------------------------------------------------------*
      SUBROUTINE DRPOT (ITASK, NLXM , NL,
     &                  TKLX , ZRTMS, WCADX , WCWPX , WCFCX , WCSTX,
     &                  WCLQT)
      IMPLICIT NONE

*     Formal parameters
      INTEGER ITASK,NLXM,NL
      REAL ZRTMS
      REAL TKLX(NLXM) , WCADX(NLXM), WCWPX(NLXM)
      REAL WCFCX(NLXM), WCSTX(NLXM), WCLQT(NLXM)

*     Local variables
      INTEGER IL

      SAVE

      IF (ITASK.EQ.1) THEN

         NL = 1
         DO 10 IL=1,NL
            TKLX(IL)  = 10.
            WCADX(IL) = 0.1
            WCWPX(IL) = 0.12
            WCFCX(IL) = 0.30
            WCSTX(IL) = 0.40
            WCLQT(IL) = WCFCX(IL)
10       CONTINUE
         ZRTMS = 10.
      END IF

      RETURN
      END

