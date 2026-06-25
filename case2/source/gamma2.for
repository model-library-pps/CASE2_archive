      REAL FUNCTION GAMMA2 (ALFA,BETA,ISET)

*     Chapter 17 in documentation WOFOST Version 4.1 (1988)

*     This function calculates gamma-distributed numbers for the
*     specified gamma distribution parameters BETA and ALFA.
*     The algorithm used is derived from that of Berman (1971).

*     Subroutines and functions called: RANDOM.
*     Called by routine RNGEN and RNDIS.

*     Author: C. Rappoldt, modified by D. van Kraalingen,

**
      SAVE

      AA = 1./ALFA
      AB = 1./(1.-ALFA)
      TR1 = EXP (-18.42/AA)
      TR2 = EXP (-18.42/AB)
10    RN1 = RANDOM (2,ISET)
      RN2 = RANDOM (2,ISET)
      IF ((RN1-TR1) .GT. 0.) GOTO 20
      S1 = 0.
      GOTO 30
20    CONTINUE
      S1 = RN1**AA
30    CONTINUE
      IF ((RN2-TR2) .GT. 0.) GOTO 40
      S2 = 0.
      GOTO 50
40    CONTINUE
      S2 = RN2**AB
50    CONTINUE
      S12 = S1 + S2
      IF (S12-1.) 60,60,10
60    Z = S1/S12
      RN3 = RANDOM (2,ISET)
      GAMMA2 = -Z*ALOG (RN3)*BETA
      RETURN
      END
