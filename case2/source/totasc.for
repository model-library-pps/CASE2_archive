*======================================================================*
* SUBROUTINE TOTASC                                                    *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Update : April 1995                                                  *
* Version: 1.1                                                         *
* Purpose: This subroutine calculates daily total gross assimilation   *
*          (DTGA) by performing a Gaussian integration over time. At   * 
*          three different times of the day, radiation is computed and * 
*          used to determine assimilation whereafter integration       *
*          takes place.                                                *
*                                                                      *
*                                                                      *
* FORMAL PARAMETERS: (I=input, O=output, C=control, IN=init, T=time)   *
* name   type description                                 units  class *
* ----   ---- -----------                                 -----  ----- *
* IDOY    I4  daynumber since 1 January                     d      T,I *
* LAT     R4  latitude of site                          degrees     I  *
* RDD     R4  daily incoming total global radiation        J/m2/d   I  *
* FRPAR   R4  Fraction PAR in RDD                           -       I  *
* INS     I4  number  of species                              -     I  *
* KDF     R[]  extinction coefficient for leaves            -       I  *
* KS      R[]  extinction coefficient for stems             -       I  *
* AMAX    R[]  actual maximum CO2-assimilation rate for   kg/ha/h   I  *
*             individual leaves                                        *
* EFF     R[]  initial light use efficiency for     kg/ha/h/J m2 s  IN *
*             leaves                                                   *
* LAI     R[]  leaf area index                              ha/ha   I  *
* SAI     R[]  stem area index                              m2/m2   I  *
* HGHT    R[]  total height of a species in the canopy        cm    I  *
* FRABS   R[]  fraction absorbed incoming global radiation    -     O  *
* DTGA    R[]  daily total gross CO2-assimilation         kg/ha/h   O  *
*                                                                      *
* FATAL ERROR CHECKS (execution terminated, message): none             *
*                                                                      *
* SUBROUTINES and FUNCTIONS called: SASTRO, SSKYC, ASSIMC              *
*                                                                      *
* FILE usage: none                                                     *
*                                                                      *
*======================================================================*
      SUBROUTINE TOTASC (IDOY, INS, LAT, RDD, FRPAR, KDF, KS, AMAX, EFF,
     &                   LAI, SAI, HGHT, HGHL,
     &                   FRABS, DTGA)

      IMPLICIT NONE

*     Formal parameters
      INTEGER IDOY, INS

      REAL LAT, RDD, FRPAR
      REAL KDF(INS),  KS(INS),  FRABS(INS)
      REAL AMAX(INS), EFF(INS), DTGA(INS)
      REAL LAI(INS),  SAI(INS), HGHT(INS), HGHL(INS)
	REAL SOLCON,ANGOT,DAYL,DAYLP,DSINB,DSINBE,SINLD,COSLD
	REAL HOUR, SINB, RDPDR , RDPDF

*     Standard local declarations
      INTEGER I1, I2, INGP, IMAX
      PARAMETER (INGP=3, IMAX=20)

      REAL XGAUSS(INGP), WGAUSS(INGP)
      REAL FGROS(IMAX),  DRPABS(IMAX), PARABS(IMAX)

      SAVE
      DATA WGAUSS /0.2778, 0.4444, 0.2778/
      DATA XGAUSS /0.1127, 0.5000, 0.8873/

*     Compute daylength and related data
      CALL SASTRO (IDOY,LAT,
     &             SOLCON,ANGOT,DAYL,DAYLP,DSINB,DSINBE,SINLD,COSLD)

*---- Assimilation set to zero and three different times of
*     the day (HOUR) 

      DO 10 I1=1,INS
         DTGA(I1)  = 0.
         DRPABS(I1) = 0.
10    CONTINUE

      DO 30 I1 = 1,INGP
         HOUR = 12.0 + DAYL * 0.5 * XGAUSS(I1)

*------- At the specified HOUR, radiation is computed and used 
*        to compute assimilation

         CALL SSKYC (HOUR, SOLCON, FRPAR, DSINBE, SINLD, COSLD, RDD,
     &               SINB, RDPDR , RDPDF)

         CALL ASSIMC (INS, SINB, RDPDR, RDPDF, AMAX, EFF,
     &                KDF, KS  , LAI,   SAI,   HGHT, HGHL,
     &              FGROS, PARABS)

*------- Integration of assimilation rate to a daily total (DTGA)
*        and absorbed radiation to a total DRPABS

         DO 20 I2=1,INS
            DTGA(I2)   = DTGA(I2)  + FGROS(I2)  * WGAUSS(I1)
            DRPABS(I2) = DRPABS(I2)+ PARABS(I2) * WGAUSS(I1)
20       CONTINUE                             ! Finish loop over species
30    CONTINUE                          ! Finish loop over the daylength
   
      DO 40 I1=1,INS
         DTGA(I1)  = DTGA(I1)*DAYL
         FRABS(I1) = (DRPABS(I1)*DAYL*3600.*2.) / RDD
40    CONTINUE

      RETURN
      END

*----------------------------------------------------------------------*
* SUBROUTINE SSKYC                                                     *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Authors: Daniel van Kraalingen                                       *
* Date   : 2-Jun-1993, Version 1.0                                     *
* Purpose: This subroutine estimates solar inclination and fluxes of   *
*          diffuse and direct irradiation at a particular time of      *
*          the day.                                                    *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* HOUR    R4  Hour for which calculations should be done     h      I  *
* SOLCON  R4  Solar constant                               W/m2     I  *
* FRPAR   R4  Fraction of total shortwave irradiation        -      I  *
*             that is photosynthetically active (PAR)                  *
* DSINBE  R4  Daily integral of sine of solar height         s      I  *
*             corrected for lower transmission at low                  *
*             elevation                                                *
* SINLD   R4  Intermediate variable from SASTRO              -      I  *
* COSLD   R4  Intermediate variable from SASTRO              -      I  *
* RDD     R4  Daily shortwave radiation                   J/m2/d    I  *
* SINB    R4  Sine of solar inclination at HOUR              -      O  *
* RDPDR   R4  Instantaneous flux of direct photo-          W/m2     O  *
*             synthetically active radiation (PAR)                     *
* RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2     O  *
*             synthetically active irradiation (PAR)                   *
*                                                                      *
* Fatal error checks: RDD <= 0                                         *
* Warnings          : ATMTR > 0.9                                      *
* Subprograms called: ERROR                                            *
* File usage        : none                                             *
*----------------------------------------------------------------------*
      SUBROUTINE SSKYC (HOUR, SOLCON, FRPAR, DSINBE, SINLD, COSLD, RDD,
     &                  SINB, RDPDR , RDPDF)
      IMPLICIT NONE

*     Formal parameters
      REAL HOUR, SOLCON, FRPAR, DSINBE, SINLD, COSLD, RDD
      REAL SINB, RDPDR, RDPDF

*     Local parameters
      REAL ATMTR, FRDIF, SOLHM, RT1

*     Hour on day that solar height is at maximum
      PARAMETER (SOLHM=12.)
      SAVE

      IF (RDD.LE.0.) CALL FATALERR
     &   ('SSKYC','total shortwave irradiation <= zero')

      ATMTR  = 0.
      FRDIF  = 0.
      RDPDF = 0.
      RDPDR = 0.

*     Sine of solar inclination, 0.2617993 is 15 degrees in radians
      SINB = SINLD+COSLD*COS ((HOUR-SOLHM)*0.2617993)

      IF (SINB.GT.0.) THEN
*        Sun is above the horizon

         RT1   = RDD*SINB*(1.+0.4*SINB)/DSINBE
         ATMTR = RT1/(SOLCON*SINB)

         IF (ATMTR.GT.0.9) WRITE (*,'(A,G12.5,A)')
     &       ' WARNING from SSKYC: ATMTR =',ATMTR,', value very large'

         IF (ATMTR.LE.0.22) THEN
            FRDIF = 1.
         ELSE IF (ATMTR.GT.0.22 .AND. ATMTR.LE.0.35) THEN
            FRDIF = 1.-6.4*(ATMTR-0.22)**2
         ELSE
            FRDIF = 1.47-1.66*ATMTR
         END IF

*        Apply lower limit to fraction diffuse
         FRDIF = MAX (FRDIF, 0.15+0.85*(1.-EXP (-0.1/SINB)))

*        Diffuse and direct PAR
         RDPDF = RT1*FRPAR*FRDIF
         RDPDR = RT1*FRPAR*(1.-FRDIF)

      END IF

      RETURN
      END

*======================================================================*
* SUBROUTINE ASSIMC                                                    *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Date   : March 1992                                                  *
* Update : April 1992                                                  *
* Version: 1.1                                                         *
* Purpose: This subroutine (for two or more species in competition)    *
*          performs a Gaussian integration over the depth of the       *
*          canopy for each species; selects five different points      *
*          (m above soil) and computes the leaf area index above each  *
*          point (LAIC), and the leaf area density (LD) and local      *
*          assimilation rate at each point. The integrated variables   *
*          are FGROS and PARABS.                                       *
*                                                                      *
*                                                                      *
* FORMAL PARAMETERS: (I=input, O=output, C=control, IN=init, T=time)   *
* name   type description                                 units  class *
* ----   ---- -----------                                 -----  ----- *
* INS     I4  number of species                             -      I   *
* SINB    R4  sine of solar inclination at HOUR                    I   *
* RDPDR   R4  Instantaneous flux of direct photo-          W/m2    O   *
*             synthetically active radiation (PAR)                     *
* RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2    O   *
*             synthetically active irradiation (PAR)                   *
* AMAX    R[] actual maximum CO2-assimilation rate       kg/ha/h   I   *
*             for individual leaves                                    * 
* EFF     R[]  initial light use efficiency for      kg/ha/h/J m2 s IN *
*             leaves                                                   *
* KDF     R[]  extinction coefficient for leaves             -      I  *
* KS      R[]  extinction coefficient for stems              -      I  *
* LAI     R[]  leaf area index                             ha/ha    I  *
* **** REMOVE SAI
* SAI     R[]  stem area index                             m2/m2    I  *
****
* HGHT    R[]  total height of a species in the canopy       m      I  *
* HGHL    R[]  lower height of a species in the canopy       m      I  *
* FGROS   R[]  canopy assimilation                        kg/ha/h   O  *
* PARABS  R[]  absorbed radiation by species in canopy     J/m/s    O  *
*                                                                      *
* FATAL ERROR CHECKS (execution terminated, message): none             *
*                                                                      *
* SUBROUTINES and FUNCTIONS called: LEAFPA or LEAFRE                   *
*                                                                      *
* FILE usage: none                                                     *
*                                                                      *
*======================================================================*
       SUBROUTINE ASSIMC (INS, SINB, RDPDR, RDPDF, AMAX, EFF,
     &                    KPDFL, KPDFS  , LAI,   SAI,   HGHT, HGHL,
     &                  FGROS, PARABS)
      IMPLICIT NONE

*     Formal parameters
      INTEGER INS
      REAL SINB, RDPDF, RDPDR
      REAL  AMAX(INS), EFF(INS), KPDFL(INS),   KPDFS(INS)
      REAL   LAI(INS), SAI(INS), HGHT(INS), HGHL(INS)
      REAL FGROS(INS), PARABS(INS)


*     Local parameters
      INTEGER I,K,IMAX,IG1,IG2,IGP1,IGP2
      REAL CSLV, RT1, RFLH, RFLS, X 
      REAL EXSDFV, EXSTV, EXDV, AFVPP, FSLLAV

      PARAMETER (CSLV=0.2, IMAX=20)

      REAL FGL(IMAX),   AFT(IMAX)
      REAL LAIC(IMAX),  SAIC(IMAX)
      REAL LD(IMAX),    SD(IMAX)

      REAL KPDRBL(IMAX), KPDRBS(IMAX)
      REAL KPDRTL(IMAX), KPDRTS(IMAX)

      REAL FGRSH(IMAX), FGRSUN(IMAX)
      REAL FGRS(IMAX),  VISSUN(IMAX)
      REAL AFVVL(IMAX), AFVTL(IMAX),  AFVDL(IMAX), AFVSHL(IMAX)
      REAL ABSNL(IMAX)

*     Gauss parameters
      PARAMETER (IGP1=5,  IGP2=3)
      REAL XGAUS1(IGP1),WGAUS1(IGP1), XGAUS2(IGP2), WGAUS2(IGP2)

      SAVE

*     Gaussian points and weights for three and five point integration
      DATA XGAUS1 /0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      DATA WGAUS1 /0.1184635,0.2393144,0.2844444,
     &             0.2393144,0.1184635/
      DATA XGAUS2 /0.1127, 0.5000, 0.8873/
      DATA WGAUS2 /0.2778, 0.4444, 0.2778/

*---- Reflection coefficients of canopy for horizontal (REFVH) and
*     spherical (REFVS) leaves
      RT1  = SQRT(1. - CSLV)
      RFLH = (1. - RT1)/(1. + RT1)
      RFLS = RFLH*2. / (1. + 1.6*SINB)

*---- Extinction coefficients for direct component of direct flux
*     (KPDRB.), total direct flux (KPDRT.) for each species;
*     canopy assimilation and PAR absorption is set to zero

      DO 10 K = 1,INS
*--     Leaves
        KPDRBL(K) = (0.5 / SINB)*KPDFL(K) / (0.8*RT1)
        KPDRTL(K) = KPDRBL(K)*RT1

*--     Stems
        KPDRBS(K)= (0.5 / SINB)*KPDFS(K) / (0.8*RT1)
        KPDRTS(K)= KPDRBS(K)*RT1

        FGROS(K) = 0.
        PARABS(K)= 0.

10    CONTINUE                                       ! Loop over species

*---- Height (m) within canopy is selected; leaf area index 
*     (LAIC, m2/m2) above each Gaussian point and leaf area 
*     density (LD, m2/m3) at each point are calculated in 
*     subroutine LEAFPA and Stem area density in LEAFRE
*     Changed for the use in CASE2: K = 1,1 [not INS]
*     Thus: assimilation only for cocoa plants, not for shade trees. 

      DO 100 K  =1,1           

       ! Start loop over the depth of the canopy
       DO 50  IG1=1,IGP1
        X      = HGHL(K) + XGAUS1(IG1) * (HGHT(K)-HGHL(K))
        EXSDFV = 0.
        EXSTV  = 0.
        EXDV   = 0.
      
        DO 20 I=1,INS

          CALL LEAFPA(HGHT(I),HGHL(I),X,LAI(I),LAIC(I),LD(I))
*          CALL LEAFRE(HGHT(I),X,SAI(I),SAIC(I),SD(I))

*-------- Exponents for light distribution functions: sum of leaf
*         area indices above point X weighted by the extinction
*         coefficients for each species

          EXSDFV = EXSDFV + KPDFL(I)  * LAIC(I) + KPDFS(I)  * SAIC(I)
          EXSTV  = EXSTV  + KPDRTL(I) * LAIC(I) + KPDRTS(I) * SAIC(I)
          EXDV   = EXDV   + KPDRBL(I) * LAIC(I) + KPDRBS(I) * SAIC(I)

20      CONTINUE                                     ! Loop over species

*------ Absorbed fluxes (J/m2 leaf/s) per species at specified
*       height in the canopy: diffuse flux, total direct flux,
*       direct component of direct flux

*------ Leaves
        AFVVL(K) = (1.-RFLH)*RDPDF*KPDFL(K)*EXP(-EXSDFV)
        AFVTL(K) = (1.-RFLS)*RDPDR*KPDRTL(K)*EXP(-EXSTV)
        AFVDL(K) = (1.-CSLV)*RDPDR*KPDRBL(K)*EXP(-EXDV)

*------ Total absorbed diffuse flux (J/m2 leaf/s) and rate of
*       photosynthesis (kg/ha/h)

*------ Shaded leaves
        AFVSHL(K)= AFVVL(K)+AFVTL(K)-AFVDL(K)
        FGRSH(K) = AMAX(K)*(1.-EXP(-AFVSHL(K)*EFF(K)/AMAX(K)))

*------ Direct flux absorbed by sunlit leaves perpendicular to the
*       direct beam (VISPP); instantaneous assimilation of sunlit
*       leaf area (FGRSUN) integrated over the sine of incidence of
*       direct light, assuming a spherical leaf angle distribution
        AFVPP = (1.-CSLV)*RDPDR/SINB

        FGRSUN(K) = 0.
        ABSNL(K)  = 0.
        DO 30 IG2 = 1,IGP2
          VISSUN(K) = AFVSHL(K) + AFVPP     * XGAUS2(IG2)
          FGRS(K)   = AMAX(K)*(1.-EXP(-VISSUN(K)*EFF(K)/AMAX(K)))
          FGRSUN(K) = FGRSUN(K) + FGRS(K)   * WGAUS2(IG2)
          ABSNL(K)  = ABSNL(K)  + VISSUN(K) * WGAUS2(IG2)
30      CONTINUE

*------ Fraction sunlit leaf area (FSLLAV) and local
*       assimilation rate (FGL, kg CO2/ha leaf/h) and absorbed fluxes
        FSLLAV = EXP(-EXDV)
        FGL(K) =  (FSLLAV*FGRSUN(K) + (1. - FSLLAV)*FGRSH(K)) *LD(K)
        AFT(K) =  (FSLLAV*ABSNL(K)  + (1. - FSLLAV)*AFVSHL(K))*LD(K)

*------ Integration of local assimilation rate to canopy
*       assimilation (FGROS) and absorbed light (PARABS)
        FGROS(K)  = FGROS(K) + FGL(K)*WGAUS1(IG1)*(HGHT(K) - HGHL(K))
        PARABS(K) = PARABS(K) + AFT(K)*WGAUS1(IG1) * (HGHT(K)-HGHL(K))

50     CONTINUE                                        ! Loop over depth
100   CONTINUE                                       ! Loop over species

      RETURN
      END

*======================================================================*
* SUBROUTINE LEAFPA                                                    *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Date   : April 1995                                                  *
* Version: 1.1                                                         *
* Purpose: This subroutine assumes a parabolic leaf area distribution; *
*          heights (HGHT, HGHL), a point X and total leaf area index   *
*          (LAI) are input; leaf area index of the canopy above        *
*          point X (LAIC) and the leaf area density (LD) at point X    *
*          are calculated.                                             *
*                                                                      *
*                                                                      *
* FORMAL PARAMETERS: (I=input, O=output, C=control, IN=init, T=time)   *
* name   type description                                 units  class *
* ----   ---- -----------                                 -----  ----- *
* HGHT    R4  total height of a species in the canopy       m     I    *
* HGHL    R4  underside of a species in the canopy          m     I    *
* X       R4  selected height at point X                    m     I    *
* LAI     R4  species leaf area index                      ha/ha  I    *
* LAIC    R4  total species leaf area index above point X  ha/ha  O    *
* LD      R4  leaf area density at point X                 m2/m3  O    *
*                                                                      *
* FATAL ERROR CHECKS (execution terminated, message):                  *
* X < 0.                                                               *
* HGHT < HGHL                                                          *
*                                                                      *
* SUBROUTINES and FUNCTIONS called: ERROR                              *
*                                                                      *
* FILE usage: none                                                     *
*                                                                      *
*======================================================================*
      SUBROUTINE LEAFPA (HGHT,HGHL,X,LAI,
     &                   LAIC,LD)
      IMPLICIT NONE

*     Formal parameters
      REAL HGHT, HGHL, X, LAI, LAIC, LD

*     Local parameters
      REAL HDIF, XPOS

      SAVE

      IF (X .LT. 0. .OR. HGHT .LE. HGHL)
     &   CALL FATALERR ('LEAFPA','Negative X or canopy height')

      IF (X .GT. HGHT) THEN
         LAIC = 0.
         LD   = 0.
      ELSEIF (X .LE. HGHT .AND. X .GE. HGHL ) THEN
         HDIF  = HGHT-HGHL
         XPOS  = X - HGHL
         LAIC = LAI - ((LAI/(HDIF)**3) * XPOS**2 * (3*HDIF - 2*XPOS))
         LD   = (LAI*6./HDIF**3) * XPOS * (HDIF - XPOS)
      ELSEIF (X .GT. 0. .AND. X .LT. HGHL) THEN
         LAIC = LAI
         LD   = 0.
      ENDIF

      RETURN
      END

*======================================================================*
***** NOTE 5-6-2001 *** THIS SUBROUTINE IS NOT USED ANYMORE**** PAZ    *
***** STEM AREA INDEX IS SET TO ZERO                                   * 
* SUBROUTINE LEAFRE                                                    *
************************************************************************
* Date   : March 1992                                                  *
* Version: 1.1                                                         *
* Purpose: This subroutine assumes a rectangular leaf area             *
*          distribution; height (HGHT), a point X and total leaf area  *
*          index (LAI) are input; leaf area index of the canopy above  *
*          point X (LAIC) and the leaf area density (LD) at point X    *
*          are calculated.                                             *
*                                                                      *
*                                                                      *
* FORMAL PARAMETERS: (I=input, O=output, C=control, IN=init, T=time)   *
* name   type description                                 units  class *
* ----   ---- -----------                                 -----  ----- *
* HGHT    R4  total height of a species in the canopy       cm    I    *
* X       R4  selected height at point X                    cm    I    *
* LAI     R4  total leaf area index                        ha/ha  I    *
* LAIC    R4  total leaf area index above point X          ha/ha  O    *
* LD      R4  leaf area density at point X                 m2/m3  O    *
*                                                                      *
* FATAL ERROR CHECKS (execution terminated, message):                  *
* X < 0.                                                               *
* HGHT < 0.                                                            *
*                                                                      *
* SUBROUTINES and FUNCTIONS called: ERROR                              *
*                                                                      *
* FILE usage: none                                                     *
*                                                                      *
*======================================================================*

      SUBROUTINE LEAFRE (HGHT,X,LAI,
     &                   LAIC,LD)
      IMPLICIT NONE

*     Formal parameters
      REAL HGHT, X, LAI, LAIC, LD

      SAVE

      IF (X .LT. 0. .OR. HGHT .LT. 0.) 
     &    CALL FATALERR ('LEAFRE', 'Negative X or canopy height')

      IF (X .LE. HGHT .AND. HGHT .GT. 0.) THEN
         LAIC = LAI * (HGHT - X)/HGHT
         LD   = LAI/HGHT
      ELSE 
         LAIC = 0.
         LD   = 0.
      ENDIF

      RETURN
      END

      INTEGER FUNCTION DATCMP (IYEAR1,IDOY1,IYEAR2,IDOY2)
      IMPLICIT NONE

*     DATCMP = -1 date of simulation is earlier than date from table
*     DATCMP =  0 date of simulation is equal to date from table
*     DATCMP =  1 date of simulation is later than date from table

*     formal parameters
      INTEGER IYEAR1, IDOY1, IYEAR2, IDOY2

*     IYEAR1 - year from table
*     IDOY1  - day from table
*     IYEAR2 - year from simulation
*     IDOY2  - day from simulation

      SAVE

      IF (IYEAR1.LT.1500.OR.IYEAR2.LT.1500.OR.IYEAR1.EQ.IYEAR2) THEN
         IF (IDOY2.LT.IDOY1) THEN
            DATCMP = -1
         ELSE IF (IDOY2.GT.IDOY1) THEN
            DATCMP = 1
         ELSE
            DATCMP = 0
         END IF
      ELSE IF (IYEAR2.LT.IYEAR1) THEN
         IF (IYEAR2.LT.IYEAR1) THEN
            DATCMP = -1
         ELSE IF (IYEAR2.GT.IYEAR1) THEN
            DATCMP = 1
         END IF
      END IF

      RETURN
      END
