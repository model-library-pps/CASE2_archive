*----------------------------------------------------------------------*
* SUBROUTINE CLIMRD                                                    *
*                                                                      *
* Use      : For CASE2 (Cacao Simulation Engine) version 2.2           *
* Author(s): Daniel van Kraalingen/Wouter Gerritsma                    *
* Date     : 07-NOV-1996, Version: 1.2                                 *
* Purpose  : Reads weather data from WOFOST format files including the *
*     generation of rainfall. The routine searches the file for the    *
*     requested year and reads weather data, latitude and altitude of  *
*     the site. Years with weather data must be in consecutive order   *
*     on file.                                                         *
*                                                                      *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning (unit)                                     class *
* ----   ---- ---------------                                    ----- *
* CLFILE  C*  Name of climate file (-)                              I  *
* IRNDAT  I4  Determines how rainfall data should be obtained       I  *
*             (0=generated; 1=distributed)(-)                          *
* RSETRG  L4  Reset rainfall generator   (only valid if IRNDAT=0)      *
*             (-)                                                   I  *
* RSETRD  L4  Reset rainfall distributor (only valid if IRNDAT=1)      *
*             (-)                                                   I  *
* IYEAR   I4  Year of simulation (y)                                I  *
* IDOY    I4  Day number within year of simulation (d)              I  *
* IUWE    I4  Unit number for weather file (-)                      I  *
* LAT     R4  Latitude of site (dec.degr.)                          I  *
* ELEV    R4  Elevation of site (m)                                 I  *
* TMMN    R4  Daily minimum temperature (degrees C)                 O  *
* TMMX    R4  Daily maximum temperature (degrees C)                 O  *
* AVRAD   R4  Daily shortwave radiation (kJ.m-2.d)                  O  *
* VAPOUR  R4  Early morning vapour pressure (mbar)                  O  *
* WN      R4  Average wind speed (m.s-1)                            O  *
* RAIN    R4  Daily amount of rainfall (mm.d-1)                     O  *
*                                                                      *
* Fatal error checks:                                                  *
* Warnings          :                                                  *
* Subprograms called:RNDIS, RNGEN, AFGEN (WOFOST subroutine)           * 
* Libraries         :TTUTIL 4.08 (FATALERR)
* File usage        :                                                  *
*----------------------------------------------------------------------*
      SUBROUTINE CLIMRD (WTRDIR, CLFILE, IRNDAT, RSETRG, RSETRD,
     &                   IYEAR , IDOY  , IUWE  , LAT   , ELEV,
     &                   TMMN  , TMMX  , AVRAD , VAPOUR, WN, RAIN)

*     formal parameters
      IMPLICIT NONE
      INTEGER IRNDAT, IYEAR, IDOY, IUWE
      REAL LAT, ELEV, TMMX, TMMN, AVRAD, VAPOUR, WN, RAIN
      CHARACTER*(*) CLFILE, WTRDIR
      LOGICAL RSETRG, RSETRD

*     local parameters 
      REAL DAYNUM, RIDOY
      REAL AFGEN
      REAL DAYNM1(28), DAYNM2(28)
      REAL TMMNTB(28), TMMXTB(28), IRRATB(28)
      REAL VAPPTB(28), WNTB(28), RAINTB(12),RAIND(12),RSERIE(366)
      CHARACTER*80 LINE, CLFILO
      INTEGER ILYEAR, I1, I2, IOS, ILDAY, LEN_TRIM
      LOGICAL OPNFIL,LEAP
      CHARACTER*80 TMPSTR1, TMPSTR2
      SAVE

*     table of day numbers for normal and leap years
      DATA DAYNM1 /1.,0.,15.,0.,45.,0.,74.,0.,105.,0.,135.,0.,166.,0.,
     &        196.,0.,227.,0.,258.,0.,288.,0.,319.,0.,349.,0.,365.,0./
      DATA DAYNM2 /1.,0.,15.,0.,45.,0.,75.,0.,106.,0.,136.,0.,167.,0.,
     &        197.,0.,228.,0.,259.,0.,289.,0.,320.,0.,350.,0.,366.,0./

      DATA CLFILO /' '/, OPNFIL /.FALSE./

*     open new file if new name is different from name of previous call,
*     set start values of local year and day variables

      IF (CLFILE.NE.CLFILO) THEN
         IF (OPNFIL) CLOSE (IUWE) 
         TMPSTR1 = CLFILE
         TMPSTR2 = WTRDIR
         I1 = LEN_TRIM (TMPSTR1)
         I2 = LEN_TRIM (TMPSTR2)
         CLFILE = WTRDIR(1:I2)//TMPSTR1(1:I1)
         CALL LOWERC (CLFILE)
         CALL FOPENG (IUWE,CLFILE,'OLD','FS',0,' ')
         OPNFIL = .TRUE.
         CLFILO = CLFILE
         CALL MOFILP (IUWE)
         ILYEAR = -99
         ILDAY  = -99
      END IF

      IF (IYEAR.NE.ILYEAR) THEN

*        requested year does not match year during previous call

*        assign new day numbers to arrays, these will be used for
*        interpolation later
         LEAP = IYEAR.GT.1500.AND.MOD (IYEAR,4).EQ.0
         DO I1=1,27,2
           IF (LEAP) THEN
             DAYNUM = DAYNM2(I1)
           ELSE
             DAYNUM = DAYNM1(I1)
           END IF
           TMMNTB(I1) = DAYNUM
           TMMXTB(I1) = DAYNUM
           IRRATB(I1) = DAYNUM
           VAPPTB(I1) = DAYNUM
           WNTB(I1)   = DAYNUM
         ENDDO

         IF (IYEAR.LT.ILYEAR) THEN
*          start at top of file if requested year is less than year
*          of previous call
           REWIND (IUWE)
           CALL MOFILP (IUWE)
         END IF

*        read next lines on file, during normal situations when the 
*        requested year (if changed) is one higher than the old year and
*        the year number on file is one higher than in the previous
*        block, these read statements are sufficient to find the
*        requested year

         READ (IUWE,'(A)') LINE
         READ (IUWE,*) ILYEAR,LAT,ELEV

*        if year on file does not match requested year, start search
*        for requested year
20       IF (IYEAR.NE.ILYEAR) THEN
           DO I1=1,13
             READ (IUWE,'(A)',IOSTAT=IOS) LINE
             IF (IOS.NE.0) THEN
               WRITE (*,*) IYEAR
               CALL FATALERR 
     &         ('CLIMRD','cannot find requested year')
               END IF
           ENDDO
           READ (IUWE,*) ILYEAR,LAT,ELEV
         GOTO 20
         END IF

*        requested year is found at this point, read weather data

         DO I1=1,12
           I2 = I1*2+2
           READ (IUWE,*,IOSTAT=IOS) TMMNTB(I2),TMMXTB(I2),IRRATB(I2),
     &           VAPPTB(I2),WNTB(I2),RAINTB(I1),RAIND(I1)
           IF (IOS.NE.0) CALL FATALERR
     &         ('CLIMRD','error while reading climate data')
         ENDDO

* 6.5    inserting monthly values in AFGEN tables
         TMMNTB(2)  = (TMMNTB(4)+TMMNTB(26))/2.
         TMMNTB(28) = TMMNTB(2)
         TMMXTB(2)  = (TMMXTB(4)+TMMXTB(26))/2.
         TMMXTB(28) = TMMXTB(2)
         IRRATB(2)  = (IRRATB(4)+IRRATB(26))/2.
         IRRATB(28) = IRRATB(2)
         VAPPTB(2)  = (VAPPTB(4)+VAPPTB(26))/2.
         VAPPTB(28) = VAPPTB(2)
         WNTB(2)  = (WNTB(4)+WNTB(26))/2.
         WNTB(28) = WNTB(2)
      END IF

*     interpolate weather variables, (this was previously done by
*     subroutine interp)

      RIDOY  = REAL (IDOY)
      TMMN   = AFGEN (TMMNTB, 28, RIDOY)
      TMMX   = AFGEN (TMMXTB, 28, RIDOY)
      AVRAD  = AFGEN (IRRATB, 28, RIDOY)
      VAPOUR = AFGEN (VAPPTB, 28, RIDOY)
      WN   = AFGEN (WNTB, 28, RIDOY)

*     generated rainfall dependent on type of weather data (either monthly
*     climate averages or monthly year averages)

      IF (IRNDAT.EQ.0) THEN
         IF (IDOY.LT.ILDAY.OR.RSETRG)
     &      CALL RNGEN (RSETRG,RAINTB,RAIND,RSERIE)
         RAIN = RSERIE(IDOY)
      ELSE IF (IRNDAT.EQ.1) THEN
         IF (IDOY.LT.ILDAY.OR.RSETRD)
     &      CALL RNDIS (RSETRD,IYEAR,RAINTB,RAIND,RSERIE)
         RAIN = RSERIE(IDOY)
      ELSE IF (IRNDAT.EQ.2) THEN
         RAIN = -99.
      ELSE
         CALL FATALERR ('CLIMRD','illegal rainfall option')
      END IF

*     update local day
      ILDAY  = IDOY

      RETURN
      END
      
