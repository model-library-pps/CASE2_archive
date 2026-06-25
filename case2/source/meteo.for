*     To do:
*     Nog een check uitvoeren op weather variabelen in het geval
*     van het WOFOST weersysteem. 
*     check fot tabs in weather files 
*----------------------------------------------------------------------*
* SUBROUTINE METEO                                                     *
*                                                                      *
* Use:   : For CASE2 (Cacao Simulation Engine) version 2.2             *
* Authors: Daniel van Kraalingen / Wouter Gerritsma                    *
* Date   : 12/02/1997, Version: 1.1                                    *
* Purpose: Provides the calling program with weather data (including   *    
*     rainfall) and evapotranspiration values. Different file formats  *
*     can be used depending on the value of IWEATH. Random rainfall    *
*     generation or distrution takes place depending on IWEATH.        *
*     Depending on the value of IWEATH, different variables are used   *
*     to define the data file be used:                                 *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning (unit)                                     class *
* ----   ---- ---------------                                    ----- *
* IYEAR   I4  Year of simulation (y)                                I  *
* IDOY    I4  Day number within year of simulation (d)              I  *
* IWEATH  I4  Options for weather file format (-)                   I  *
* WSTAT   C*  Status code from weather system (-)                   I  *
* IFLAG   I4  Flag handeling error message for weathr.lib (-)       I  *
* WTRDIR  C*  Directory where weather data are stored (-)           I  *
* CNTR    C*  Country name of weather data (-)                      I  *
* ISTN    I4  Station number of weather data (-)                    I  *
* CLFILE  C*  String containing directory and filename (-)          I  *
* IRNDAT  I4  Determines how rainfall data should be obtained          *
*             (0=generated; 1=distributed) (-)                      I  *
* RSETRG  L4  Reset rainfall generator   (only valid if IRNDAT=0)      *
*             (-)                                                   I  *
* RSETRD  L4  Reset rainfall distributor (only valid if IRNDAT=1)      *
*             (-)                                                   I  *
* IUWE    I4  Unit number for weather file (-)                      I  *
* IURA    I4  Unit number by which external rainfall file can be       *
*             opened (-)                                            I  *
* IYEARR  I4  Year for which rainfall data are requested (y)        I  *
* RAFILE  C*  File from which rainfall data is to be obtained (-)   I  *
* ANGA    R4  A value of Angstrom formula (-)                       O  *
* ANGB    R4  B value of Angstrom formula (-)                       O  *
* LAT     R4  Latitude of site (dec.degr.)                          O  *
* LONG    R4  Longitude of site (dec.degrees)                       O  *
* ELEV    R4  Elevation of site (m)                                 O  *
* RDD     R4  Daily shortwave radiation (J.m-2.d)                   O  *
* TMMN    R4  Daily minimum temperature (degrees C)                 O  *
* TMMX    R4  Daily maximum temperature (degrees C)                 O  *
* RAIN    R4  Daily amount of rainfall (mm.d-1)                     O  *
* WN      R4  Average wind speed (m.s-1)                            O  *
* VP      R4  Early morning vapour pressure (kPa)                   O  *
*                                                                      *
* Fatal error checks:                                                  *
* Warnings          :                                                  *
* Subprograms called: libraries TTUTIL & WEATHR                        *
*                     STINFO, CLIMRD, RNREAL                           *
*                     potentially: DBMRD, REPRD, DBMRD, WCURRD         *                                  *
* File usage        :                                                  *
*----------------------------------------------------------------------*
      SUBROUTINE METEO (IYEAR , IDOY  , IWEATH, WSTAT,
     &                  IFLAG , WTRDIR, CNTR  , ISTN  ,
     &                  CLFILE, IRNDAT, RSETRG, RSETRD, IUWE  ,
     &                  IURA  , IYEARR, RAFILE,
     &                  ANGA  , ANGB  , LAT   , LONG  , ELEV  ,
     &                  RDD   , TMMN  , TMMX  , RAIN  , WN    , VP    )

      IMPLICIT NONE

*     formal parameters
      INTEGER IYEAR, IDOY, IWEATH, ISTN, IUWE, IRNDAT , IURA, IYEARR
      INTEGER IFLAG
      REAL LAT , LONG, ELEV, ANGA, ANGB
      REAL TMMX, TMMN, RDD, VP, WN, RAIN  
      CHARACTER*(*) WTRDIR, CNTR, CLFILE, RAFILE, WSTAT

*    
*     CHARACTER*(*) STNAM, GEODIR, DRVDIR, DBMDIR, DBRDIR

      LOGICAL RSETRG, RSETRD

*     local parameters
      INTEGER   ISTAT1, ISTAT2
      REAL      AVRAD , VAPOUR           
      REAL      OANGA , OANGB
      LOGICAL   WTRMES
                             
      SAVE

      IF (IWEATH.EQ.0.OR.IWEATH.EQ.1) THEN

*        monthly climate averages or monthly year averages from WOFOST
*        format, routine has one unit permanently in use (IUWE)

         CALL CLIMRD (WTRDIR, CLFILE, IRNDAT, RSETRG, RSETRD,
     &                 IYEAR, IDOY  , IUWE  , LAT   , ELEV,
     &                 TMMN , TMMX  , AVRAD , VAPOUR, WN  , RAIN)

*        Eigenlijk nog een check uitvoeren op weather variabelen
         WSTAT = '111111'

*        adjust units
         RDD  = AVRAD*1.E6
         VP   = VAPOUR/10.
         
      ELSE IF (IWEATH.EQ.2) THEN 
*       Patch for Angstrom parameters      
        IF ((ANGA.GT.0.).OR.(ANGB.GT.0.)) THEN
           OANGA = ANGA
           OANGB = ANGB
        ENDIF    
      
*-------Open weather file and read station information and return
*       weather data for start day of simulation.
*       Check status of weather system, WTRMES flags if warnings or errors
*       have occurred during the whole simulation. WTRTER flags if the run
*       should be terminated because of missing weather

        CALL STINFO (IFLAG , WTRDIR, ' ', CNTR, ISTN, IYEAR,
     &             ISTAT1, LONG  , LAT, ELEV, ANGA, ANGB)
        CALL WEATHR (IDOY  , ISTAT2, RDD, TMMN, TMMX, VP, WN, RAIN)
        IF (ISTAT1.NE.0.OR.ISTAT2.NE.0) WTRMES = .TRUE.
        WSTAT  = '444444'
        IF (ABS (ISTAT2).GE.111111) THEN
          WRITE (WSTAT,'(I6)') ABS (ISTAT2)
        ELSE IF (ISTAT2.EQ.0) THEN
           WSTAT = '111111'
        END IF

*-------Conversion of total daily radiation from kJ/m2/d to J/m2/d
        RDD = RDD*1.E3

*       Patch for Angstrom parameters      
        IF ((ANGA.EQ.0.).AND.(ANGB.EQ.0.)) THEN
           ANGA = OANGA
           ANGB = OANGB
        ENDIF    
                                    
      ELSE IF (IWEATH.EQ.3) THEN

*        dbmeteo format, routine has two units permanently in use
*        (IUWE+1 and IUWE+2)

*        not within standard WOFOST !
*        CALL DBMRD (DBMDIR, DRVDIR,
*    &               STNAM , IUWE+1, IYEAR, IDOY,
*    &               LAT   , ELEV  ,
*    &               TMMN  , TMMX  , AVRAD, VAPOUR, WN, RAIN)
         CALL FATALERR ('METEO','DBMETEO option not supported')

*        adjust units
         RDD  = AVRAD*1.E6
         VP   = VAPOUR/10.

      ELSE IF (IWEATH.EQ.4) THEN

*       repaired dbmeteo format, routine has one unit permanently in use
*       (IUWE+3)
*       not within standard CASE2
*       CALL REPRD (DBRDIR,
*     &              STNAM, IUWE+3, IYEAR, IDOY,
*     &                LAT, ELEV  ,
*     &               TMMN, TMMX  , AVRAD, VAPOUR, WN, RAIN)
        CALL FATALERR ('METEO','current weather option not supported')

*       adjust units
        RDD  = AVRAD*1.E6
        VP   = VAPOUR/10.

      ELSE IF (IWEATH.EQ.5) THEN

*       read from wgrid.cur weather current year over grids
*       routine has one unit permanently in use (IUWE+4)

*       CALL WCURRD (GEODIR, DRVDIR, STNAM, IUWE+4, IYEAR, IDOY,
*    &                LAT   , ELEV  ,
*    &                TMMN  , TMMX  , AVRAD, VAPOUR, WN , RAIN)
        CALL FATALERR ('METEO','current weather option not supported')

*       adjust units
        RDD  = AVRAD*1.E6
        VP   = VAPOUR/10.

      END IF

*     get rain from external file if necessary
      
      IF ((IWEATH.EQ.0 .OR. IWEATH .EQ. 1).AND. IRNDAT.EQ.2) 
     &   CALL RNREAL (IURA, RAFILE, IYEARR, IDOY, RAIN)

      RETURN
      END

