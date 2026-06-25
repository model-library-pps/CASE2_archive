      PROGRAM MAIN
      CALL FSE
      END

*----------------------------------------------------------------------*
* SUBROUTINE MODELS (MODEL2.FOR)                                       *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Authors: Daniel van Kraalingen                                       *
*          Wouter Gerritsma, Liesje Mommer                             *
*          Pieter Zuidema (modifications)                              *
* Date   : 5-Jul-1993, Version: 1.1                                    *
* Update   11 May 1998                                                 * 
* Modifications: February 2002                                          *
* Purpose: This subroutine is the interface routine between the FSE-   *
*          driver and the simulation models. This routine is called    *
*          by the FSE-driver at each new task at each time step. It    *
*          can be used by the user to specify calls to the different   *
*          models that have to be simulated                            *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITO  I4  Unit used for output file                      -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEI1  C*  Name of input file no. 1                       -      I  *
* FILEI2  C*  Name of input file no. 2                       -      I  *
* FILEI3  C*  Name of input file no. 3                       -      I  *
* FILEI4  C*  Name of input file no. 4                       -      I  *
* FILEI5  C*  Name of input file no. 5                       -      I  *
* OUTPUT  L4  Flag to indicate if output should be done      -      I  *
* TERMNL  L4  Flag to indicate if simulation is to stop      -     I/O *
* OUTPUTFQ I4 Code indicating the output frequency           -      I  *
* SOILTYPE I4 Code indicating the selecrted soil type        -      I  *
* DOY     R4  Day number within year of simulation (REAL)    d      I  *
* IDOY    I4  Day number within year of simulation (INTEGER) d      I  *
* YEAR    R4  Year of simulation (REAL)                      y      I  *
* IYEAR   I4  Year of simulation (INTEGER)                   y      I  *
* TIME    R4  Time of simulation                             d      I  *
* STTIME  R4  Start time of simulation                       d      I  *
* FINTIM  R4  Finish time of simulation                      d      I  *
* DELT    R4  Time step of integration                       d      I  *
* ANGA    R4  Regression coefficient in Angstrom formula     -      I  *
* ANGB    R4  Regression coefficient in Angstrom formula     -      I  *
* FRPAR   R4  Fraction PAR in shortwave radiation            -      I  *
* LAT     R4  Latitude of site                            dec.degr. I  *
* LONG    R4  Longitude of site                           dec.degr. I  *
* ELEV    R4  Elevation of site                              m      I  *
* WSTAT   C6  Status code from weather system                -      I  *
* WTRTER  L4  Flag whether weather can be used by model      -      O  *
* RDD     R4  Daily shortwave radiation                   J/m2/d    I  *
* TMMN    R4  Daily minimum temperature                  degrees C  I  *
* TMMX    R4  Daily maximum temperature                  degrees C  I  *
* VP      R4  Early morning vapour pressure                 kPa     I  *
* WN      R4  Average wind speed                            m/s     I  *
* RAIN    R4  Daily amount of rainfall                     mm/d     I  *
*                                                                      *
* Fatal error checks: ETMOD not Pen/Mak/PT; PRODLEVL not 1 or 2;       *
*                     PLTMOD not CACAO/NOCROP                          *
* Warnings          : none                                             *
* Subprograms called: TTUTIL subroutines                               *
*                     CASE2, NOCROP, SETPMD, SETPTD, DRSAHE, DRPOT,    *
*                     INTERCEPT                                        *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &                   FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &                   OUTPUT, TERMNL, OUTPUTFQ,SOILTYPE,
     &                   DOY   , IDOY  , YEAR  , IYEAR,
     &                   TIME  , STTIME, FINTIM, DELT ,
     &                   ANGA  , ANGB  , FRPAR ,
     &                   LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &                   RDD   , TMMN  , TMMX  , VP   , WN, RAIN)
      IMPLICIT NONE

*     Formal parameters
      INTEGER   I1, ITASK, IUNITD, IUNITO, IUNITL, IDOY, IYEAR
      INTEGER   OUTPUTFQ,SOILTYPE
      REAL ANGA, ANGB, TMDI, ETD, ETRD, ETAE, EVSC, DT, FRPAR
      REAL DOY, YEAR, TIME, STTIME, FINTIM, DELT, LAT
      REAL RDD, RF, RFS, GAI,TMMN, TMMX, TMDA, VP, WN 
	REAL RAIN, RAINS, PINT, RAINCU 
      REAL EVSW, EVSWCU, DRAICU, TRWCU
	REAL LONG, ELEV
      REAL MAX
      CHARACTER FILEIT*(*), FILEI1*(*), FILEI2*(*)
      CHARACTER FILEI3*(*), FILEI4*(*), FILEI5*(*)
      LOGICAL   OUTPUT, TERMNL, WTRTER
      CHARACTER WSTAT*6
      
*     Local variables
      INTEGER PRODLEVL, ISURF
      CHARACTER WUSED*6
      CHARACTER*80 ETMOD, PLTMOD
      
*     Water balance declarations
      INTEGER NLXM, NL
      PARAMETER (NLXM=10)
      REAL TRWL(NLXM) , TKL(NLXM) , WCAD(NLXM)   , WCWP(NLXM)
      REAL WCFC(NLXM) , WCST(NLXM), FLXQT(NLXM+1), WCLQT(NLXM)
      REAL FLXCU(NLXM+1), ZRTMS
      REAL WCWET
      LOGICAL GIVEN

      SAVE

*     Code for the use of RDD, TMMN, TMMX, VP, WN, RAIN (in that order)
*     a letter 'U' indicates that the variable is Used in calculations
      DATA WUSED /'------'/
      DATA GIVEN /.FALSE./

*     Check weather data availability
      IF (ITASK.EQ.1.OR.ITASK.EQ.2.OR.ITASK.EQ.4) THEN
         IF (WSTAT(6:6).EQ.'4') THEN
            RAIN       = 0.
            WSTAT(6:6) = '1'
            IF (.NOT.GIVEN) THEN
               WRITE (IUNITL,'(2A)') ' Rain not available,',
     &           ' value set to zero, (patch DvK, Jan 1995)'
               GIVEN = .TRUE.
            END IF
         END IF
         DO I1=1,6
*           Is there an error in the I1-th weather variable ?
            IF (WUSED(I1:I1).EQ.'U' .AND.
     &          WSTAT(I1:I1).EQ.'4') THEN
               WTRTER = .TRUE.
               TERMNL = .TRUE.
               RETURN
            END IF
         END DO
      END IF

      IF (ITASK.EQ.1) THEN
*        ----------------------
*        Initialization section
*        ----------------------

*        Read modules to be used from timer.dat
         CALL RDINIT (IUNITD, IUNITL, FILEIT)
         CALL RDSCHA ('PLTMOD', PLTMOD)
         CALL UPPERC (PLTMOD)
         CALL RDSCHA ('ETMOD' , ETMOD)
         CALL UPPERC (ETMOD)
         CLOSE (IUNITD)

*        Read production level from basic.dat
         CALL RDINIT (IUNITD, IUNITL, FILEI3)
         CALL RDSINT ('PRODLEVL', PRODLEVL)
         CLOSE (IUNITD)

*        Write line to mark start of new run
         WRITE (IUNITO,'(A,76A1)') '*',('=',I1=1,76)

*        Log messages to output file
         WRITE (IUNITO,'(A)') '*'
         WRITE (IUNITO,'(A)') '* FSE driver info:'
         WRITE (IUNITO,'(A,T7,A,I5,A,I4,A)')
     &     '*','Year:',IYEAR,', day:',IDOY,', System start'
         WRITE (IUNITO,'(A)') '*'
         WRITE (IUNITO,'(A)') '* Modules used:'

*        Choose and check evapotranspiration modules
         IF (ETMOD.EQ.'PENMAN') THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','SETPMD: Penman evapotranspiration'
            WUSED(1:5) = 'UUUUU'
         ELSE IF (ETMOD.EQ.'MAKKINK') THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','SETMKD: Makkink evapotranspiration'
            WUSED(1:3) = 'UUU'
         ELSE IF (ETMOD.EQ.'PRIESTLEY TAYLOR') THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','SETPTD: Priestley Taylor evapotranspiration'
            WUSED(1:3) = 'UUU'
         ELSE
            CALL FATALERR
     &      ('MODELS','unknown module name for evapotranspiration')
         END IF

*        Choose and check water balance modules
         IF (PRODLEVL.EQ.1) THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','DRPOT : Water balance for potential situations'
         ELSE IF (PRODLEVL.EQ.2) THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','DRSAHE: Tipping bucket water balance version 1.4'
            WUSED(6:6) = 'U'
         ELSE
            CALL FATALERR
     &      ('MODELS','Wrong value for PRODLEVL')
         END IF

*        Choose and check crop modules
         IF (PLTMOD.EQ.'CACAO'.AND.
     &       PRODLEVL.EQ.1) THEN
            WRITE (IUNITO,'(A,T7,A,/,A,T7,A)')
     &        '*','CASE2:  Cacao at potential production 2.2'
            WUSED(1:3) = 'UUU'
         ELSE IF (PLTMOD.EQ.'CACAO'.AND.
     &            PRODLEVL.EQ.2) THEN
            WRITE (IUNITO,'(A,T7,A,/,A,T7,A)')
     &        '*','CASE2: Cacao at water limited production 2.2'
            WUSED(1:3) = 'UUU'
         ELSE IF (PLTMOD.EQ.'NO CROP') THEN
            WRITE (IUNITO,'(A,T7,A)')
     &        '*','NO CROP: no crop'
         ELSE
            CALL FATALERR
     &      ('MODELS','unknown module name for plant')
         END IF
        
*        Avoid FORCHECK errors
         WCLQT(1) = -99.
         WCST(1)  = -99.

*        Rain interception and throughfall
         CALL INTERCEPT (ITASK, IUNITD, IUNITL, FILEI1,
     &                                  RAIN, RAINS, PINT)
   
*        Water balance module
         IF (PRODLEVL.EQ.2) THEN
           CALL DRSAHE (ITASK , IUNITD, IUNITO, FILEI2, SOILTYPE,
     &                IDOY  , IYEAR , DELT  , OUTPUT,
     &                NLXM  , NL    , EVSC  , RAINS , TRWL,
     &                TKL   , ZRTMS ,
     &                WCAD  , WCWP  , WCFC  , WCST  ,
     &                EVSW  , FLXQT , WCLQT ,
     &                DRAICU, EVSWCU, RAINCU, TRWCU , FLXCU)
         ELSE IF (PRODLEVL.EQ.1) THEN
           CALL DRPOT (ITASK, NLXM , NL,
     &               TKL  , ZRTMS, WCAD , WCWP , WCFC , WCST,
     &               WCLQT)
         END IF
      
*        Crop module
         IF (PLTMOD.EQ.'CACAO') THEN
           CALL CASE2 (PLTMOD, ITASK , IUNITD, IUNITO, IUNITL, FILEI1,
     &            FILEI3, OUTPUT, TERMNL, OUTPUTFQ,
     &            DOY   , IDOY  , IYEAR , DELT  , TIME  , STTIME,
     &            LAT   , FRPAR , RDD   , TMMN  , TMMX  ,
     &            NLXM  , NL    , TRWL  , TKL   ,
     &            WCLQT , WCWP  , WCFC  , WCST  ,
     &            EVSC  , ETRD  , ETAE  , PINT  ,
     &            GAI   , RAIN )
         ELSE IF (PLTMOD.EQ.'NO CROP') THEN						
           CALL NOCROP (NLXM, TRWL, GAI)               
         END IF

      ELSE IF (ITASK.EQ.2) THEN
*        ------------------------
*        Rate calculation section
*        ------------------------

*        Reflection of soil (see: van Laar et al 1992)
         RFS = 0.25*(1.-0.5*WCLQT(1)/WCST(1))

*        Total reflection (see: van Laar et al 1992)
         RF  = RFS*EXP(-0.5*GAI)+0.25*(1.-EXP(-0.5*GAI))

*        Calculate average temperature
         TMDA = (TMMX+TMMN)/2.

         IF (ETMOD.EQ.'PENMAN') THEN
*          Penman evapotranspiration
           CALL SETPMD (IDOY,LAT,ISURF,RF,ANGA,ANGB,TMDI,
     &                  RDD, TMDA, WN,VP, 
     &                  ETD, ETRD,ETAE,  DT)
*           Calculate potential soil evaporation taking into account
*           the standing crop
            EVSC = EXP (-0.5*GAI)*(ETRD+ETAE)
         ELSE IF (ETMOD.EQ.'MAKKINK') THEN
*           Makkink evapotranspiration
*           NOTE: This module is not used in CASE2
*           CALL SETMKD (RDD, TMDA, ETD)
*           Estimate radiation driven and wind and humidity driven part
*           ETRD = 0.75*ETD
*           ETAE = ETD-ETRD
*           Calculate potential soil evaporation taking into account
*           the standing crop
*           EVSC = EXP (-0.5*GAI)*ETD
         ELSE IF (ETMOD.EQ.'PRIESTLEY TAYLOR') THEN
*           Priestley Taylor evapotranspiration
            CALL SETPTD (IDOY,LAT,RF,RDD,TMDA,ETD)
*           Estimate radiation driven and wind and humidity driven part
            ETRD = 0.75*ETD
            ETAE = ETD-ETRD
*           Calculate potential soil evaporation taking into account
*           the standing crop
            EVSC = EXP (-0.5*GAI)*ETD
         END IF
*        Make sure potential soil evaporation is always positive
*        the amount of dew is unreliable anyhow
         EVSC = MAX (EVSC, 0.)
      
*        Rain interception and throughfall
         CALL INTERCEPT (ITASK, IUNITD, IUNITL, FILEI1,
     &                                  RAIN, RAINS, PINT)
   
*        Water balance module
         IF (PRODLEVL.EQ.2) THEN
           CALL DRSAHE (ITASK , IUNITD, IUNITO, FILEI2, SOILTYPE,
     &                IDOY  , IYEAR , DELT  , OUTPUT,
     &                NLXM  , NL    , EVSC  , RAINS , TRWL,
     &                TKL   , ZRTMS ,
     &                WCAD  , WCWP  , WCFC  , WCST  ,
     &                EVSW  , FLXQT , WCLQT ,
     &                DRAICU, EVSWCU, RAINCU, TRWCU , FLXCU)
         ELSE IF (PRODLEVL.EQ.1) THEN
           CALL DRPOT (ITASK, NLXM , NL,
     &               TKL  , ZRTMS, WCAD , WCWP , WCFC , WCST,
     &               WCLQT)
         END IF
      
*        Crop module
         IF (PLTMOD.EQ.'CACAO') THEN
           CALL CASE2 (PLTMOD, ITASK , IUNITD, IUNITO, IUNITL, FILEI1,
     &            FILEI3, OUTPUT, TERMNL, OUTPUTFQ,
     &            DOY   , IDOY  , IYEAR , DELT  , TIME  , STTIME,
     &            LAT   , FRPAR , RDD   , TMMN  , TMMX  ,
     &            NLXM  , NL    , TRWL  , TKL   ,
     &            WCLQT , WCWP  , WCFC  , WCST  ,
     &            EVSC  , ETRD  , ETAE  , PINT  ,
     &            GAI   , RAIN )
         ELSE IF (PLTMOD.EQ.'NO CROP') THEN						
           CALL NOCROP (NLXM, TRWL, GAI)               
         END IF

      ELSE IF (ITASK.EQ.3) THEN
*       --------------------
*       Integration section
*       --------------------
  
*        Water balance module
         IF (PRODLEVL.EQ.2) THEN
           CALL DRSAHE (ITASK , IUNITD, IUNITO, FILEI2, SOILTYPE,
     &                IDOY  , IYEAR , DELT  , OUTPUT,
     &                NLXM  , NL    , EVSC  , RAINS , TRWL,
     &                TKL   , ZRTMS ,
     &                WCAD  , WCWP  , WCFC  , WCST  ,
     &                EVSW  , FLXQT , WCLQT ,
     &                DRAICU, EVSWCU, RAINCU, TRWCU , FLXCU)
         ELSE IF (PRODLEVL.EQ.1) THEN
           CALL DRPOT (ITASK, NLXM , NL,
     &               TKL  , ZRTMS, WCAD , WCWP , WCFC , WCST,
     &               WCLQT)
         END IF
      
*        Crop module
         IF (PLTMOD.EQ.'CACAO') THEN
           CALL CASE2 (PLTMOD, ITASK , IUNITD, IUNITO, IUNITL, FILEI1,
     &            FILEI3, OUTPUT, TERMNL, OUTPUTFQ,
     &            DOY   , IDOY  , IYEAR , DELT  , TIME  , STTIME,
     &            LAT   , FRPAR , RDD   , TMMN  , TMMX  ,
     &            NLXM  , NL    , TRWL  , TKL   ,
     &            WCLQT , WCWP  , WCFC  , WCST  ,
     &            EVSC  , ETRD  , ETAE  , PINT  ,
     &            GAI   , RAIN )
         ELSE IF (PLTMOD.EQ.'NO CROP') THEN						
           CALL NOCROP (NLXM, TRWL, GAI)               
         END IF
       
      ELSE IF (ITASK.EQ.4) THEN
*       ----------------
*       Terminal section
*       ----------------
         WRITE (IUNITO,'(A)') '*'
         WRITE (IUNITO,'(A)') '* FSE driver info:'
         WRITE (IUNITO,'(A,T7,A,I5,A,I4,A)')
     &     '*','Year:',IYEAR,', day:',IDOY,', System end'

*        Water balance module
         IF (PRODLEVL.EQ.2) THEN
           CALL DRSAHE (ITASK , IUNITD, IUNITO, FILEI2, SOILTYPE,
     &                IDOY  , IYEAR , DELT  , OUTPUT,
     &                NLXM  , NL    , EVSC  , RAINS , TRWL,
     &                TKL   , ZRTMS ,
     &                WCAD  , WCWP  , WCFC  , WCST  ,
     &                EVSW  , FLXQT , WCLQT ,
     &                DRAICU, EVSWCU, RAINCU, TRWCU , FLXCU)
         ELSE IF (PRODLEVL.EQ.1) THEN
           CALL DRPOT (ITASK, NLXM , NL,
     &               TKL  , ZRTMS, WCAD , WCWP , WCFC , WCST,
     &               WCLQT)
         END IF

*        Crop module
         IF (PLTMOD.EQ.'CACAO') THEN
           CALL CASE2 (PLTMOD, ITASK , IUNITD, IUNITO, IUNITL, FILEI1,
     &            FILEI3, OUTPUT, TERMNL, OUTPUTFQ,
     &            DOY   , IDOY  , IYEAR , DELT  , TIME  , STTIME,
     &            LAT   , FRPAR , RDD   , TMMN  , TMMX  ,
     &            NLXM  , NL    , TRWL  , TKL   ,
     &            WCLQT , WCWP  , WCFC  , WCST  ,
     &            EVSC  , ETRD  , ETAE  , PINT  ,
     &            GAI   , RAIN )
         ELSE IF (PLTMOD.EQ.'NO CROP') THEN						
           CALL NOCROP (NLXM, TRWL, GAI)               
         END IF

      END IF

      RETURN

      END

*----------------------------------------------------------------------*
* SUBROUTINE INTERCEPT                                                 *
*                                                                      *
* Author : Wouter Gerritsma / Liesje Mommer                            *
* Date   : 1998                                                        *
* Version: 1.0                                                         *
*                                                                      *
* Purpose: Calculate rain intercepted by cocoa crop                    *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEP   C*  Name of file with plant data                   -      I  *
* RAIN    R4  Daily amount of rainfall                     mm/d     I  *
* RAINS   R4  Daily amount of rainfall reaching the soil   mm/d     O  *
* PINT    R4  Daily amount of interecepted rain            mm/d     O  *
*                                                                      *
* Fatal error checks: none                                             *
* Warnings          : none                                             *
* Subprograms called: TTUTIL routines                                  *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE INTERCEPT (ITASK, IUNITD, IUNITL, FILEP,
     &                                  RAIN, RAINS, PINT)
      
      IMPLICIT NONE
      
*     Formal parameters
      INTEGER ITASK, IUNITD, IUNITL
      REAL RAIN, RAINS, PINT
      CHARACTER FILEP*(*)
            
*     Local declarations
      REAL TFALL, STEMFL
      REAL TFALA, TFALB, STFLA, STFLB
      REAL LIMIT
      
      IF (ITASK.EQ.1) THEN
*       ----------------------
*       Initialization section
*       ----------------------
      
*       Initialize input file
        CALL RDINIT (IUNITD, IUNITL, FILEP)

*       Stemflow and troughfall parameters
        CALL RDSREA ('TFALA'  , TFALA  )
        CALL RDSREA ('TFALB'  , TFALB  )
        CALL RDSREA ('STFLA'  , STFLA  )
        CALL RDSREA ('STFLB'  , STFLB  )
 
        CLOSE (IUNITD)
      
      ELSE IF (ITASK.EQ.2) THEN
*       ------------------------
*       Rate calculation section
*       ------------------------
      
*       Troughfall (TFALL), Stemflow (STEMFL) and
*       Rainfall interception (PINT)
        TFALL    = LIMIT (0., RAIN, TFALA*RAIN + TFALB)      ! [mm d-1]
        STEMFL   = LIMIT (0., RAIN, STFLA*RAIN + STFLB)      ! [mm d-1]
        PINT     = MAX(0.,RAIN - TFALL - STEMFL)             ! [mm d-1]

*       Rainfall reaching the soil (RAINS)(= TFALL + STEMFL)
        RAINS    = RAIN-PINT                                 ! [mm d-1]
      
      END IF

      RETURN
      END
