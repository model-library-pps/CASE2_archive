*----------------------------------------------------------------------*
* SUBROUTINE SWSE                                                      *
* Authors: Daniel van Kraalingen                                       *
* Date   : 20-Jun-1994, Version: 1.1                                   *
* Purpose: This subroutine calculates a reduction factor on water      *
*          uptake accounting for the effect of drought.                *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* WCL     R4  Volumetric water content in compartment   cm3/cm3     I  *
* P       R4  Soil water depletion factor                   -       I  *
* WCWET   R4  Volumetric water content at which water   cm3/cm3     I  *
*             stress starts to affect water uptake                     *
* WCWP    R4  Volumetric water content at wilting point cm3/cm3     I  *
* WCFC    R4  Volumetric water content at field cap.    cm3/cm3     I  *
* WCST    R4  Volumetric water content at saturation    cm3/cm3     I  *
* WSE     R4  Reduction factor on soil water uptake         -       O  *
*             as function of drought stress                            *
*                                                                      *
* Fatal error checks: WCWET < WCCR                                     *
* Warnings          : none                                             *
* Subprograms called: FATALERR                                         *
* File usage        : none                                             *
*----------------------------------------------------------------------*

      SUBROUTINE SWSE (WCL,P,WCWET,WCWP,WCFC,WCST,WSE)

      IMPLICIT NONE

*     Formal parameters
      REAL WCL,P,WCWET,WCWP,WCFC,WCST,WSE

*     Local variables
      REAL WCCR, LIMIT
      SAVE

*     Calculation of critical water content, transition point
*     from water-limited to potential transpiration rate
      WCCR = WCWP + (1.-P) * (WCFC - WCWP)

      IF (WCWET.LT.WCCR) CALL FATALERR ('SWSE','WCWET < WCCR')

      IF (WCL.GT.WCWET) THEN
*        Water content larger than optimal
*        growth reduction occurs
         WSE = (WCST-WCL)/(WCST-WCWET)
      ELSE
         IF (WCL.GE.WCCR) THEN
*           Water content is at optimal level
*           no growth reduction
            WSE = 1.
         ELSE IF (WCL.LT.WCWP) THEN
*           Water content is below wilting point
*           growth reduction occurs
            WSE = 0.
         ELSE IF (WCCR.NE.WCWP) THEN
*           Water content is at suboptimal level
*           growth reduction occurs
            WSE = (WCL-WCWP)/(WCCR-WCWP)
         END IF
      END IF

*     Limit reduction between valid range
      WSE = LIMIT (0., 1., WSE)

      RETURN
      END
