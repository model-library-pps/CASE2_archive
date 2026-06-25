*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*              FORTRAN Simulation Environment (FSE 2.1m)               *
*                            November, 1996                            *
*----------------------------------------------------------------------*
*     Use:     For CASE2 (Cacao Simulation Engine) version 2.2         *
*----------------------------------------------------------------------*
*                                                                      *
*     FSE 2.1m is a simulation environment suited for simulation of    *
*     biological processes in time, such as crop and vegetation growth,*
*     insect population development etc.                               *
*                                                                      *
*     The MAIN program, subroutine FSE and subroutine MODELS are       *
*     programmed by D.W.G. van Kraalingen, DLO Institute for           *
*     Agrobiological and Soil Fertility Research (AB-DLO),             *
*     PO Box 14, 6700 AA, Wageningen, The Netherlands (e-mail:         *
*     d.w.g.van.kraalingen@ab.agro.nl).                                *
*                                                                      *
*     Modifications by Wouter Gerritsma (1998)                         *
*     Modifications by Pieter Zuidema (2002)                           *
*                                                                      *
*     FSE version 2.1 is described in:                                 *
*        Kraalingen, D.W.G. van, 1995. The FSE system for crop         *
*        simulation, version 2.1. Quantitative Approaches in Systems   *
*        Analysis; no.1. DLO Research Institute for Agrobiology and    *
*        Soil Fertility, Wageningen. The C.T. de Wit Graduate School   *
*        for Production Ecology. 70 pp.                                *
*                                                                      *
*     This version has been adapted to handle calls monthly weather    *
*        data files in a WOFOST format                                 * 
*                                                                      *
*     Data files needed for FSE 2.1m:                                  *
*          (excluding data files used by models called from MODELS):   *
*        - CONTROL.DAT (contains file names to be used),               *
*        - timer file whose name is specified in CONTROL.DAT,          *
*        - optionally, a rerun file whose name is specified in         *
*          CONTROL.DAT,                                                *
*        - weather data files as specified in timer file               *
*     Object libraries needed for FSE 2.1m:                            *
*        - TTUTIL (at least version 3.2)                               *
*        - WEATHER (at least version from 17-Jan-1990)                 *
*----------------------------------------------------------------------*
      SUBROUTINE FSE

      USE CHART

      IMPLICIT NONE

*-----Standard declarations for simulation and output control
      INTEGER   ITASK   , INSETS, ISET  , IPFORM, IL, LEN_TRIM
      LOGICAL   OUTPUT  , TERMNL, RDINQR, STRUNF, ENDRNF
      CHARACTER COPINF*1, DELTMP*1
      INTEGER   INPRS   , STRUN , ENDRUN, I1

      INTEGER   IMNPRS
      PARAMETER (IMNPRS=100)
      CHARACTER PRSEL(IMNPRS)*11

*-----Declarations for time control
      INTEGER   IDOY, IYEAR
      REAL      DELT, DOY, FINTIM, PRDEL, STTIME, TIME, YEAR

*-----Declarations for meteo system
      INTEGER   IFLAG    , IUWE  , IWEATH , ISTN, IRNDAT
      INTEGER   IYEARR   , IURA  , NRYEARS
      REAL      ANGA     , ANGB  , ELEV   , LAT , LONG  , FRPAR
      REAL      RDD      , TMMN  , TMMX   , VP  , WN    , RAIN  
      LOGICAL   WTRMES   , WTRTER, RSETRG , RSETRD
      CHARACTER WTRDIR*80, CLFILE*80, RAFILE*80
      CHARACTER CNTR*7,    WSTAT*6,   DUMMY*1

*-----Declarations for soil type and production level
      INTEGER   SOILTYPE, PRODLEVL

*-----Declarations for output
      INTEGER   OUTPUTFQ

*-----Declarations for location seletction
      INTEGER   LOCATION, STRTYR, ENDYR

*     Declarations for higher IWEATH selections      
*     CHARACTER GEODIR*80, DRVDIR*80, DBRDIR*80, DBMDIR*80
*     CHARACTER STNAM*80 , RAFILE*80   
      
*-----Declarations for file names and units
      INTEGER   IUNITR   , IUNITD   , IUNITO   , IUNITL   , IUNITC
      CHARACTER FILEON*80, FILEOL*80
      CHARACTER FILEIC*80, FILEIR*80, FILEIT*80
      CHARACTER FILEI1*80, FILEI2*80, FILEI3*80, FILEI4*80, FILEI5*80
      CHARACTER STRING*80
	INTEGER   SIGLEN

*-----Declarations for observation data facility
      INTEGER   INOD , IOD

      INTEGER   IMNOD
      PARAMETER (IMNOD=100)
      INTEGER   IOBSD(IMNOD)

*-----For communication with OBSSYS routine
      COMMON /FSECM1/ YEAR,DOY,IUNITD,IUNITL,TERMNL

      SAVE

*-----File name for control file and empty strings for input
*     files 1-5. WTRMES flags any messages from the weather system

      DATA FILEIC /'CONTROL.DAT'/
      DATA FILEI1 /' '/, FILEI2 /' '/, FILEI3 /' '/
      DATA FILEI4 /' '/, FILEI5 /' '/
      DATA WTRMES /.FALSE./

      DATA STRUNF /.FALSE./, ENDRNF /.FALSE./

*-----Unit numbers for control file (C), data files (D),
*     output file (O), log file (L) and rerun file (R).
      IUNITC = 10
      IUNITD = 20
      IUNITO = 30
      IUNITL = 40
      IUNITR = 50 

*-----Unit number for WOFOST weather file
      IUWE   = 60

*-----Open control file and read names of normal output file, log file
*     and rerun file (these files cannot be used in reruns)
      CALL RDINIT (IUNITC,0, FILEIC)
      CALL RDSCHA ('FILEON', FILEON)
      CALL RDSCHA ('FILEOL', FILEOL)
      CALL RDSCHA ('FILEIR', FILEIR)

*     check if start run number was found, if there, read it
      IF (RDINQR('STRUN'))  THEN
         CALL RDSINT ('STRUN',STRUN)
         STRUNF = .TRUE.
      END IF
*     check if end run number was found, if there, read it
      IF (RDINQR('ENDRUN')) THEN
         CALL RDSINT ('ENDRUN',ENDRUN)
         ENDRNF = .TRUE.
      END IF
      CLOSE (IUNITC)

*-----Open output file and possibly a log file
      CALL FOPENS  (IUNITO, FILEON, 'NEW', 'DEL')
      IF (FILEOL.NE.FILEON) THEN
         CALL FOPENS  (IUNITL, FILEOL, 'NEW', 'DEL')
      ELSE
         IUNITL = IUNITO
      END IF

c*     initialization of logfile for processing of end_of_run values
c      CALL OPINIT

*-----See if rerun file is present, and if so read the number of rerun
*     sets from rerun file

      CALL RDSETS (IUNITR, IUNITL, FILEIR, INSETS)

*======================================================================*
*======================================================================*
*                                                                      *
*                   Main loop and reruns begin here                    *
*                                                                      *
*======================================================================*
*======================================================================*

      IF (.NOT.ENDRNF) THEN
*        no end run was found in control.dat file
         ENDRUN = INSETS
      ELSE
         ENDRUN = MAX (ENDRUN, 0)
         ENDRUN = MIN (ENDRUN, INSETS)
      END IF

      IF (.NOT.STRUNF) THEN
*        no start run was found in control.dat file
         STRUN = 0
      ELSE
         STRUN = MAX (STRUN, 0)
         STRUN = MIN (STRUN, ENDRUN)
      END IF

      CALL ChartInit (IUNITO+2)

      DO 10 ISET=STRUN,ENDRUN

      WRITE (*,'(A)') '   FSE 2.1m: Initialize model'

*-----Select data set
      CALL RDFROM (ISET, .TRUE.)

      CALL ChartSetRunID (ISET)

*======================================================================*
*                                                                      *
*                        Initialization section                        *
*                                                                      *
*======================================================================*

      ITASK  = 1
      TERMNL = .FALSE.
      WTRTER = .FALSE.
      RSETRG = .TRUE.
      RSETRD = .TRUE.

*-----Read names of timer file and input files 1-5 from control
*     file (these files can be used in reruns)
      CALL RDINIT (IUNITC,IUNITL,FILEIC)
      CALL RDSCHA ('FILEIT', FILEIT)
      IF (RDINQR ('FILEI1')) CALL RDSCHA ('FILEI1', FILEI1)
      IF (RDINQR ('FILEI2')) CALL RDSCHA ('FILEI2', FILEI2)
      IF (RDINQR ('FILEI3')) CALL RDSCHA ('FILEI3', FILEI3)
      IF (RDINQR ('FILEI4')) CALL RDSCHA ('FILEI4', FILEI4)
      IF (RDINQR ('FILEI5')) CALL RDSCHA ('FILEI5', FILEI5)
      CLOSE (IUNITC)

*-----Read information on location,simulation period and soil type
      CALL RDINIT (IUNITD,IUNITL,FILEI3)
      CALL RDSINT ('SOILTYPE', SOILTYPE)
      CALL RDSINT ('LOCATION', LOCATION)
      CALL RDSINT ('IYEAR', IYEAR)
      CALL RDSINT ('NRYEARS', NRYEARS)
      CALL RDSINT ('PRODLEVL', PRODLEVL)
      CLOSE (IUNITD)

*-----Read time, control and miscellaneous weather variables from timer file
      CALL RDINIT (IUNITD  , IUNITL, FILEIT)
      CALL RDSREA ('STTIME', STTIME)
      CALL RDSREA ('DELT'  , DELT  )
      CALL RDSINT ('IPFORM', IPFORM)
      CALL RDSCHA ('COPINF', COPINF)
      CALL RDSCHA ('DELTMP', DELTMP) 
      CALL RDSREA ('FRPAR' , FRPAR)
      IF (RDINQR('IFLAG'))  CALL RDSINT ('IFLAG' , IFLAG)

*-----Read weather parameters depending on the location specified in BASIC.dat 
*     (Added PAZ 1-2002)
*     First for daily or monthly weather (LOCATION <51)
      IF (LOCATION.LT.51) THEN 
*-------read IWEATH
        STRING = 'IWEATH'
  	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) CALL RDSINT (STRING, IWEATH)
        IF (IWEATH .EQ. 0)  CALL FATALERR ('FSE', 
     &   'No information for this LOCATION in timer.dat')
*-------read IRNDAT
        STRING = 'IRNDAT'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        CALL RDSINT (STRING, IRNDAT)
*-------read WTRDIR
        STRING = 'WTRDIR'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) CALL RDSCHA (STRING, WTRDIR)
*-------read CLFILE
        STRING = 'CLFILE'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) CALL RDSCHA (STRING, CLFILE)
*-------read ISTN
        STRING = 'ISTN'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) CALL RDSINT (STRING, ISTN)
*-------read CNTR
        STRING = 'CNTR'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) CALL RDSCHA (STRING, CNTR)
*-------read IRNDAT
        STRING = 'STRTYR'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        CALL RDSINT (STRING, STRTYR)
*-------read ENDYR
        STRING = 'ENDYR'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        CALL RDSINT (STRING, ENDYR)
*     Then for long-term weather (LOCATION > 50)
	ELSE
        STRING = 'CLFILE'
	  SIGLEN = LEN_TRIM(STRING)
        CALL ADDINT(STRING,SIGLEN,LOCATION)
        IF (RDINQR(STRING)) THEN
	    CALL RDSCHA (STRING, CLFILE)
  	  ELSE
	    CALL FATALERR ('FSE', 
     &   'No information for this LOCATION in timer.dat')
        ENDIF
        WTRDIR = 'c:\case2\weather\longterm\'
        IWEATH = 0
        IRNDAT = 0
        STRTYR = 1000
        ENDYR  = 1000
      ENDIF
      CLOSE (IUNITD)

*-----Check whether simulation period as specified in BASIC.dat does not 
*     exceed period of available weather data
      IF (IYEAR.LT.STRTYR.AND.IYEAR.NE.1000.) CALL  FATALERR ('FSE',
     &  'No weather data available for start year specified in BASIC.dat
     &   choose later start year (IYEAR)')
      IF (IYEAR.GT.ENDYR.AND.IYEAR.NE.1000.) CALL  FATALERR ('FSE',
     &  'No weather data available for start year specified in BASIC.dat
     &   choose earlier start year (IYEAR)')
      IF ((IYEAR + NRYEARS - 1).GT.ENDYR.AND.IYEAR.NE.1000.)
     &   CALL  FATALERR ('FSE',
     &  'Simulation period exceeds period for which weather data are   
     &   available: select an earlier start year (IYEAR)                
     &   or a shorter period (NRYEARS) in BASIC.dat')

*-----Calculate FINTIM in days, based on NRYEARS and information on leap years
      FINTIM = 0.
      DO 11 I1 = IYEAR, (IYEAR + NRYEARS - 1)
        IF (I1.GT.1500.AND.MOD (I1,4).EQ.0) THEN
	     FINTIM = FINTIM + 366.
	  ELSE 
	     FINTIM = FINTIM + 365.
        END IF
11    CONTINUE


*-----Read parameter on output frequency from basic.dat and set output
*     parameters accordingly (added PAZ)
      CALL RDINIT (IUNITD  , IUNITL, FILEI3)
      CALL RDSINT ('OUTPUTFQ' , OUTPUTFQ)
      CLOSE (IUNITD)
      CALL RDINIT (IUNITD  , IUNITL, FILEIT)
      IF (OUTPUTFQ.EQ.1) THEN
	  PRDEL = 0.
        CALL RDAINT ('IOBSD' , IOBSD, IMNOD, INOD)
        CALL RDACHA ('PRSEL1',PRSEL,IMNPRS,INPRS)        
      ELSE IF (OUTPUTFQ.EQ.2) THEN
	  PRDEL = 10.
        CALL RDACHA ('PRSEL2',PRSEL,IMNPRS,INPRS)        
      ELSE IF (OUTPUTFQ.EQ.3) THEN 
	  PRDEL = 1.
        CALL RDACHA ('PRSEL3',PRSEL,IMNPRS,INPRS)
      ELSE 
	  CALL FATALERR
     &  ('FSE','Output frequency value (OUTPUTFQ) is invalid')
	END IF
      CLOSE (IUNITD)

*-----See if observation data variable exists, if so read it
*      INOD = 0
*      IF (RDINQR('IOBSD')) THEN
*        CALL RDAINT ('IOBSD' , IOBSD, IMNOD, INOD)
*         IF (IOBSD(1).EQ.0) INOD = 0
*      END IF

*-----See if variable with print selection exists, if so read it
*      INPRS = 0
*      IF (RDINQR('PRSEL')) CALL RDACHA ('PRSEL',PRSEL,IMNPRS,INPRS)


*-----Initialize TIMER and OUTDAT routines
      CALL TIMER2 (ITASK, STTIME, DELT, PRDEL, FINTIM,
     &             IYEAR, TIME  , DOY , IDOY , TERMNL, OUTPUT)
      YEAR = REAL (IYEAR)
      CALL OUTDAT (ITASK, IUNITO, 'TIME', TIME)

      CALL ChartInitialGroup

*-----Open weather file and read station information and return
*     weather data for start day of simulation.
      CALL METEO (IYEAR , IDOY  , IWEATH, WSTAT,
     &            IFLAG , WTRDIR, CNTR  , ISTN  ,
     &            CLFILE, IRNDAT, RSETRG, RSETRD, IUWE  ,
     &            IURA  , IYEARR, RAFILE,
     &            ANGA  , ANGB  , LAT   , LONG  , ELEV  ,
     &            RDD   , TMMN  , TMMX  , RAIN  , WN    , VP    )
 
c*-----initialize OBSSYS routine
c      IF (ITASK.EQ.1) CALL OBSINI

*-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &             FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &             OUTPUT, TERMNL, OUTPUTFQ, SOILTYPE,
     &             DOY   , IDOY  , YEAR  , IYEAR ,
     &             TIME  , STTIME, FINTIM, DELT  , 
     &             ANGA  , ANGB  , FRPAR ,
     &             LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &             RDD   , TMMN  , TMMX  , VP    , WN, RAIN)


*======================================================================*
*                                                                      *
*                      Dynamic simulation section                      *
*                                                                      *
*======================================================================*

      WRITE (*,'(A)') '   FSE 2.1m: DYNAMIC loop'

20    IF (.NOT.TERMNL) THEN

*----------------------------------------------------------------------*
*                     Integration of rates section                     *
*----------------------------------------------------------------------*

      IF (ITASK.EQ.2) THEN

*--------Carry out integration only when previous task was rate
*        calculation

         ITASK = 3

*--------Call routine that handles the different models
         CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &                FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &                OUTPUT, TERMNL, OUTPUTFQ, SOILTYPE,
     &                DOY   , IDOY  , YEAR  , IYEAR ,
     &                TIME  , STTIME, FINTIM, DELT  ,
     &                ANGA  , ANGB  , FRPAR , 
     &                LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &                RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

*--------Turn on output when TERMNL logical is set to .TRUE.
         IF (TERMNL.AND.PRDEL.GT.0.) OUTPUT = .TRUE.

      END IF

*----------------------------------------------------------------------*
*               Calculation of driving variables section               *
*----------------------------------------------------------------------*

      ITASK = 2

*-----Write time of output to screen and file
      CALL OUTDAT (2, 0, 'TIME', TIME)

      CALL ChartNewGroup
      CALL ChartOutputRealScalar('TIME', TIME)

      IF (OUTPUT) THEN
         IF (ISET.EQ.0) THEN
            WRITE (*,'(13X,A,I5,A,F7.2)')
     &        'Default set, Year:', IYEAR, ', Day:', DOY
         ELSE
            WRITE (*,'(13X,A,I3,A,I5,A,F7.2)')
     &        'Rerun set:', ISET, ', Year:', IYEAR, ', Day:', DOY
         END IF
      END IF

*-----Get weather data for new day and flag messages 
      CALL METEO (IYEAR , IDOY  , IWEATH, WSTAT,
     &            IFLAG , WTRDIR, CNTR  , ISTN  ,
     &            CLFILE, IRNDAT, RSETRG, RSETRD, IUWE  ,
     &            IURA  , IYEARR, RAFILE,
     &            ANGA  , ANGB  , LAT   , LONG  , ELEV  ,
     &            RDD   , TMMN  , TMMX  , RAIN  , WN    , VP    )
  
*----------------------------------------------------------------------*
*               Calculation of rates and output section                *
*----------------------------------------------------------------------*

*-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &             FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &             OUTPUT, TERMNL, OUTPUTFQ, SOILTYPE,
     &             DOY   , IDOY  , YEAR  , IYEAR ,
     &             TIME  , STTIME, FINTIM, DELT  , 
     &             ANGA  , ANGB  , FRPAR ,
     &             LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &             RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

      IF (TERMNL.AND..NOT.OUTPUT.AND.PRDEL.GT.0.) THEN
*--------Call model routine again if TERMNL is switched on while
*        OUTPUT was off (this call is necessary to get output to file
*        when a finish condition was reached and output generation
*        was off)
         IF (ISET.EQ.0) THEN
            WRITE (*,'(13X,A,I5,A,F7.2)')
     &        'Default set, Year:', IYEAR, ', Day:', DOY
         ELSE
            WRITE (*,'(13X,A,I3,A,I5,A,F7.2)')
     &        'Rerun set:', ISET, ', Year:', IYEAR, ', Day:', DOY
         END IF
         OUTPUT = .TRUE.
         CALL OUTDAT (2, 0, 'TIME', TIME)
         CALL ChartNewGroup
         CALL ChartOutputRealScalar('TIME', TIME)
         CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &                FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &                OUTPUT, TERMNL, OUTPUTFQ, SOILTYPE,
     &                DOY   , IDOY  , YEAR  , IYEAR ,
     &                TIME  , STTIME, FINTIM, DELT  ,
     &                ANGA  , ANGB  , FRPAR ,
     &                LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &                RDD   , TMMN  , TMMX  , VP    , WN, RAIN)
      END IF

*----------------------------------------------------------------------*
*                             Time update                              *
*----------------------------------------------------------------------*

*-----Check for FINTIM, OUTPUT and observation days
      CALL TIMER2 (ITASK, STTIME, DELT, PRDEL, FINTIM,
     &             IYEAR, TIME  , DOY , IDOY , TERMNL, OUTPUT)
      YEAR = REAL (IYEAR)
      DO 30 IOD=1,INOD,2
         IF (IYEAR.EQ.IOBSD(IOD).AND.IDOY.EQ.IOBSD(IOD+1))
     &       OUTPUT = .TRUE.
30    CONTINUE

      GOTO 20
      END IF

*======================================================================*
*                                                                      *
*                           Terminal section                           *
*                                                                      *
*======================================================================*

      ITASK = 4

      WRITE (*,'(A)') '   FSE 2.1m: Terminate model'

      CALL ChartTerminalGroup

*-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL,
     &             FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5,
     &             OUTPUT, TERMNL, OUTPUTFQ, SOILTYPE,
     &             DOY   , IDOY  , YEAR  , IYEAR ,
     &             TIME  , STTIME, FINTIM, DELT  , 
     &             ANGA  , ANGB  , FRPAR , 
     &             LAT   , LONG  , ELEV  , WSTAT , WTRTER,
     &             RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

*-----Generate output file dependent on option from timer file
      IF (IPFORM.GE.4) THEN
         IF (INPRS.EQ.0) THEN
            CALL OUTDAT (IPFORM, 0, 'Simulation results',0.)
         ELSE
*           Selection of output variables was in timer file
*           write tables according to output selection array PRSEL
            CALL OUTSEL (PRSEL,IMNPRS,INPRS,IPFORM,'Simulation results')
         END IF
      END IF

      IF (WTRTER) THEN
         WRITE (*,'(/,A,/,/,/)')
     &     ' The run was terminated due to missing weather'
         WRITE (IUNITO,'(/,A,/,/,/)')
     &     ' The run was terminated due to missing weather'
         IF (IUNITO.NE.IUNITL) WRITE (IUNITL,'(/,A,/,/,/)')
     &     ' The run was terminated due to missing weather'
      END IF

*-----Delete temporary output file dependent on switch from timer file
      IF (DELTMP.EQ.'Y'.OR.DELTMP.EQ.'y') CALL OUTDAT (99, 0, ' ', 0.)
      
10    CONTINUE

      IF (INSETS.GT.0) CLOSE (IUNITR)

*-----If input files should be copied to the output file,
*     copy rerun file (if present) and timer file and if there, input
*     files 1-5

      IF (COPINF.EQ.'Y'.OR.COPINF.EQ.'y') THEN
         IF (INSETS.GT.0) CALL COPFL2 (IUNITR, FILEIR, IUNITO, .TRUE.)
         CALL COPFL2 (IUNITD, FILEIT, IUNITO, .TRUE.)
         IF (FILEI1.NE.' ') CALL COPFL2 (IUNITD, FILEI1, IUNITO, .TRUE.)
         IF (FILEI2.NE.' ') CALL COPFL2 (IUNITD, FILEI2, IUNITO, .TRUE.)
         IF (FILEI3.NE.' ') CALL COPFL2 (IUNITD, FILEI3, IUNITO, .TRUE.)
         IF (FILEI4.NE.' ') CALL COPFL2 (IUNITD, FILEI4, IUNITO, .TRUE.)
         IF (FILEI5.NE.' ') CALL COPFL2 (IUNITD, FILEI5, IUNITO, .TRUE.)
      END IF

*-----Delete all .TMP files that were created by the RD* routines
*     during simulation
      CALL RDDTMP (IUNITD)

*-----Write to screen which files contain what
      IL = LEN_TRIM (FILEON)
      WRITE (*,'(/,3A)') ' File: ',FILEON(1:IL),
     &  ' contains simulation results'
      WRITE (*,'(2A)') ' File: WEATHER.LOG',
     &  ' contains messages from the weather system'
      IL = LEN_TRIM (FILEOL)
      WRITE (*,'(3A,/)') ' File: ',FILEOL(1:IL),
     &  ' contains messages from the rest of the model'

*-----Write message to screen and output file if warnings and/or errors
*     have occurred from the weather system, pause and wait for return
*     from user to make sure he has seen this message

      IF (WTRMES) THEN

         WRITE (*,'(/,A,/,A,/,A)') ' WARNING from FSE 2.1m:',
     &     ' There have been errors and/or warnings from',
     &     ' the weather system, check file WEATHER.LOG'
         WRITE (IUNITO,'(A,/,A,/,A)') ' WARNING from FSE 2.1m:',
     &     ' There have been errors and/or warnings from',
     &     ' the weather system, check file WEATHER.LOG'

         WRITE (*,'(A)') ' Press <Enter>'
         READ  (*,'(A)') DUMMY

      END IF

*-----Close output file and temporary file of OUTDAT
      CLOSE (IUNITO)
      CLOSE (IUNITO+1)

*-----Close log file (if used)
      IF (FILEOL.NE.FILEON) CLOSE (IUNITL)

*-----Close log file of weather system
      CLOSE (91)

c*-----Write end_of_run values to file
c      CALL OPWRIT (IUNITO)

      RETURN
      END
