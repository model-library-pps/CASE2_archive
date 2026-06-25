* ---------------------------------------------------------- *
* weather -- get weathr data
* ---------------------------------------------------------- *
* PURPOSE
*       WEATHER makes 6 measurements (irradiation, minus temperature,
*       maximum temperature, early morning vapour pressure, mean
*       windspeed and precipitation) available for planth growth
*       modeling programs and other programs requiring a realistic
*       "weather environment" on a daily basis.
* AUTHOR
*       M Verbeek
* DESIGN
*       M Verbeek, D v Kraalingen, C ten Cate.
* (c)
*       Center Agro Biological Research (CABO)
*       P.O. Box 14
*       Wageningen
*       Netherlands
* VERSION
*       1
*       2 : MAXSTN (maximum station code) is 999 in WSDAOP
*       3 : removed string concatenation "//" write variable list
*           (this is was not standard f77)
* DATE
*       17 january 1990
* ------------------------------------------------------------- *

* ------------------------------------------------------------- *
* stinfo -- (re)set weathr system parameters
* ------------------------------------------------------------- *
* PURPOSE
*
*       set (or reset) default values for location of data
*       files and the name of the log file to IPATH and ILOG.
*       Return information about data: the coordinates and altitude
*       of the weather station and in what whay the irradiation
*       data where obtained.
* PARAMETERS
*       name   type        description
*       ----   ----        ---------------------------------------
*       --- in  ---
*       IFLAG   int        output flags
*       IPATH   C*(*)      PATH to data files
*       ILOG    C*(*)      name of logfile
*       ICNT    C*(*)      name of country
*       ISTN    int        code for station
*       IYEAR   int        year of measurments
*        --- out ---
*       LON     real       longitude
*       LAT     real       latitude
*       ALT     real       altitude
*       A       real       first parameter of radiation conversion
*       B       real       second parameter of radiation conversion
* DEFAULTS
*       PATH:   ' '        Default directory
*       LOG:    ' '        Default name = "weather.log"
*       others: -          no defaults
* RETURNS
* ------------------------------------------------------------- *

      SUBROUTINE STINFO (IFLAG,IPATH,ILOG,ICNT,ISTN,IYEAR,
     &                   STAT,LON,LAT,ALT,A,B)

*     IMPLICIT NONE
*     --------- parameters --------
*     --- in ---
      CHARACTER*(*)     IPATH, ILOG, ICNT
      INTEGER           IFLAG, ISTN, IYEAR
*     --- out ---
      INTEGER           STAT
      REAL              LON, LAT, ALT, A, B
*     --------- constants ---------
*                unit var is not used (file closed)
      INTEGER    FSCLSD
      INTEGER    MAXCNT
      CHARACTER  DEFLOG*20
      PARAMETER  (
     &           DEFLOG = 'weather.log' ,
     &           MAXCNT =  6       ,
     &           FSCLSD = -1       )
*     --------- external ----------
      INTEGER  WSILEN
      LOGICAL  WSDAOP, WSATTR, WSRDDA, WSFLGS


*     --------- COMMON: WSCFIL -------------------------------- *
*                       Weathr Subsystem Char FILe names
*                       note: MAXFNM depends on system.
      INTEGER    MAXFNM
      PARAMETER (MAXFNM = 256)

      CHARACTER*(MAXFNM) BPATH, FNAME, LOG
      COMMON   /WSCFIL/  BPATH, FNAME, LOG
*     --------- common wscfil end ----------------------------- *


*     --------- COMMON: WSNFIL  ------------------------------- *
*                       Weathr Subsystem Num  FILe flags and logf. unit.
*     WARNFL  - .TRUE. if to create warings (screen and/or file)
*     ERRFL   - .TRUE. if to create error messages  ( " )

*                Number Of Flags
      INTEGER    NOF
      PARAMETER  (NOF = 4)

      LOGICAL    OF(NOF)
*                unit of logfile
      INTEGER    LUNIT
*                error flag from STINFO to WEATHR (12345 = ok)
      INTEGER    STERR
      LOGICAL    WARNFL, ERRFL

      COMMON /WSNFIL/ OF, LUNIT, WARNFL, ERRFL, STERR
*     --------- common wsnfil end ----------------------------- *

*     --------- local variables ----
      INTEGER   FLAG
*               default logfile name
      LOGICAL   FIRST, DUMMY
*               data unit
      INTEGER   DUNIT
*               datafile identification
      CHARACTER BCNT*(MAXCNT)
      CHARACTER NEWLOG*(MAXFNM)
      INTEGER   BSTN, BYEAR
      REAL      BLON, BLAT, BALT, BA, BB

      SAVE
*     --------- data -------------
      DATA FIRST /.TRUE./
*                 nonsense value
      DATA FLAG   /-54321/
      DATA BCNT   /' '/
      DATA BSTN   /-1/
      DATA BYEAR  /-1/
      DATA BLON   /-199./
      DATA BLAT   /-99./
      DATA BALT   /-99./
      DATA BA     /-99./
      DATA BB     /-99./

*     ------- start module --------
      IF (FIRST) THEN
          FIRST = .FALSE.
*         --- init of variables in common done here
          BPATH = '<no path>'
          FNAME = '<no file>'
          LUNIT = FSCLSD
          LOG   = DEFLOG
       ENDIF

*     --- magic value in STERR if exit ok: STERR = 12345,
*     ---                        on error:       =  -1
      STERR = -1

*     --- check IFLAG, if changed, reset FLAGS
      IF (IFLAG .NE. FLAG) THEN
         IF (.NOT. WSFLGS(IFLAG, STAT)) THEN
*            --- set flag to maximum output and continue
             FLAG = 1111
             DUMMY = WSFLGS(FLAG, STAT)
         ELSE
             FLAG = IFLAG
         ENDIF
      ENDIF

*     --- if ILOG = ' ' the logfile has the default name
      IF (WSILEN(ILOG) .EQ. 0) THEN
         NEWLOG = DEFLOG
      ELSE
         NEWLOG = ILOG
      ENDIF

*     --- see if new log file requested; close old logfile
      IF (NEWLOG .NE. LOG) THEN
         IF (LUNIT .NE. FSCLSD) THEN
            CLOSE(LUNIT)
            LUNIT = FSCLSD
         ENDIF
         LOG = NEWLOG
      ENDIF

*     --- return NIL value's on error
      LON = -199.
      LAT = -99.
      ALT = -99.
      A = -99.
      B = -99.

*     --- if datafile id or PATH changed, initialise buffer
      IF  ((ICNT .NE. BCNT)   .OR. (ISTN .NE. BSTN) .OR.
     &     (IYEAR .NE. BYEAR) .OR. (IPATH .NE. BPATH)) THEN
*        --- open data file, checks parameters
         IF (.NOT. WSDAOP(IPATH,ICNT,ISTN,IYEAR,STAT,DUNIT)) RETURN
*        --- get station location and conversion factors
         IF (.NOT. WSATTR(DUNIT,STAT,LON,LAT,ALT,A,B)) RETURN
*        --- get data in buffer
         IF (.NOT. WSRDDA(DUNIT,LAT, A, B, STAT)) RETURN
*        --- store *valid* parameters in buffer variables
         BPATH = IPATH
         BCNT = ICNT
         BYEAR = IYEAR
         BSTN = ISTN
         BLON = LON
         BLAT = LAT
         BALT = ALT
         BA = A
         BB = B
*        --- close data file
         CLOSE (DUNIT)
         DUNIT = FSCLSD
      ELSE
*        --- buffer already initialised, just return buffer variables
         LON = BLON
         LAT = BLAT
         ALT = BALT
         A   = BA
         B   = BB
      ENDIF

      STAT = 0
      STERR = 12345
      RETURN
      END

* --------------------------------------------------------------------- *
* weathr -- Get weathr data for IDAY from data file opend with STINFO
* --------------------------------------------------------------------- *
* SUBROUTINE WEATHR
*       Returns weather data for a day specified in IDAY from a
*       data file opened with STINFO.
* PARAMETERS
*       name     type     description
*       --- in  ---
*       IDAY     int      day for which weather data are returned
*       --- out ---
*       STAT     int      return status: negative = error,
*                                        zero     = ok,
*                                        positive = warning.
*       RIRRAD    real    irradiation
*       RTMIN     real    minus temperature
*       RTMAX     real    maximum temperature
*       REMPR     real    early morning vapour pressure
*       RMWIND    real    mean wind speed
*       RPRECI    real    precipitation
* DESCRIPTION
*       WEATHR reads six "weather parameters" from the buffer filled by
*       STINFO: irradiation, lowest temperature, highest temperature,
*       early morning vapour pressure and precipitation for the
*       specified IDAY.
* STATUS
*       ISTAT is zero when all data are available.
*       The data returned by WEATHR may or may not be original
*       measurements: the status of these values is returned in
*       ISTAT. Each digit in ISTAT represents one status value for
*       one measurement. The order of digits is the same as in
*       the parameter list.  ISTAT is negative when data are missing
*       or an exception has occured.
* DIGIT FLAGS in ISTAT
*       1   value is an original measurement
*       2   value is an interpolation
*       3   value is an estimate
*       4   value is missing
* EXAMPLES
*       0        All values are OK
*       -444441  All values except precipitation are missing
*       3111111  Irradiation is an estimate.
*       2222223  All values are interpolations, precipiation is an
*                estimate.
*       -1       WEATHR called with wrong day number.
* -------------------------------------------------------------- *

      SUBROUTINE WEATHR(IDAY,
     &           STAT, RIRRAD,RTMIN,RTMAX,REMPR,RMWIND,RPRECI)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
      INTEGER IDAY
*     --- out ---
      REAL      RIRRAD, RTMIN, RTMAX, REMPR, RMWIND, RPRECI
      INTEGER   STAT
*     --------- external ----------
*     void      WSITOA
*     void      WSMESS
      INTEGER   WSILEN


*     --------- COMMON: WSCFIL -------------------------------- *
*                       Weathr Subsystem Char FILe names
*                       note: MAXFNM depends on system.
      INTEGER    MAXFNM
      PARAMETER (MAXFNM = 256)

      CHARACTER*(MAXFNM) PATH, FNAME, LOG
      COMMON   /WSCFIL/  PATH, FNAME, LOG
*     --------- common wscfil end ----------------------------- *

*     --------- COMMON: WSNFIL  ------------------------------- *
*                       Weathr Subsystem Num  FILe flags and logf. unit.
*     OF      - Output File flags
*     OFxxxW  - warning flag,   OFxxxF - fatal error flag
*     OFLOGx  - logfile flag,   OFOUTx - output (screen) flag
*     WARNFL  - .TRUE. if to create warings (screen and/or file)
*     ERRFL   - .TRUE. if to create error messages  ( " )

*                Number Of Flags
      INTEGER    NOF
      PARAMETER (
     &           NOF = 4)

      LOGICAL    OF(NOF)
*                unit of logfile
      INTEGER    LUNIT
*                error flag from STINFO to WEATHR (12345 = ok)
      INTEGER    STERR
      LOGICAL    WARNFL, ERRFL

      COMMON /WSNFIL/ OF, LUNIT, WARNFL, ERRFL, STERR
*     --------- common wsnfil end ----------------------------- *

*     --------- COMMON: WSNBUF ------------------------------- *
*                       Weathr Subsystem Num BUFfer
*     --- BFATTR: attributes of data: DFINT, DFEST, DFUNDE
*     --- BFVALS: values
      INTEGER    NRDAYS, NRDP
      PARAMETER (
     &           NRDP = 6      ,
     &           NRDAYS = 366  )

      REAL       BFVALS(NRDAYS, NRDP)
      INTEGER    BFATTR(NRDAYS, NRDP)

      COMMON /WSNBUF/ BFVALS, BFATTR
*     --------- common wsnbuf end ----------------------------- *

*     --------- constants ---------
*     --- DF: Digit Flags returned in ISTAT
*     --- DFINT : interpolated
*     --- DFEST : estimated
*     --- DFUNDE: undefined
      INTEGER    DFINT, DFEST, DFUNDE
      PARAMETER (
     &           DFINT = 2,
     &           DFEST = 3,
     &           DFUNDE = 4)

*                order of measurements in buffer (and file)
      INTEGER    NIRRAD, NTMIN, NTMAX, NEMPR, NMWIND, NPRECI
      PARAMETER (NIRRAD = 1,
     &           NTMIN = 2,
     &           NTMAX  = 3,
     &           NEMPR = 4,
     &           NMWIND = 5,
     &           NPRECI = 6)

*             wrong day number
      INTEGER SPDAY
*             size of message for WSMESS
      INTEGER MAXMSG
*             internal write wrong
      INTEGER SSWEAT
*             weathr called after error from STINFO or STINFO not called
      INTEGER SOINIT
      PARAMETER (
     &        SOINIT = -21  ,
     &        SSWEAT = -902 ,
     &        MAXMSG = 100  ,
     &        SPDAY = -1    )
*     ------ local variables ------
*             attributes
      INTEGER ATTI
      LOGICAL ISMISS, ISESTI, ISINTP
      CHARACTER MSG*(MAXMSG), MTYPE*(MAXMSG), HLP*10

      SAVE
*     ------- start module --------
*     --- return NIL values on error
      RIRRAD = -99.
      RTMIN = -99.
      RTMAX = -99.
      REMPR = -99.
      RMWIND = -99.
      RPRECI = -99.

*     --- signal from STINFO (12345 "magic value": successful exit)
      IF (STERR .NE. 12345) THEN
         STAT = SOINIT
         CALL WSMESS(
     &        'Error in WEATHR: called after error from STINFO or '//
     &        'STINFO not called', STAT)
         RETURN
      ENDIF

*     --- check day
      IF ((IDAY .LT. 1) .OR. (IDAY .GT. NRDAYS)) THEN
         STAT = SPDAY
         IF (ERRFL) THEN
            CALL WSITOA(IDAY,HLP)
            MSG = 'Error in WEATHR: called with wrong day: '//HLP
            CALL WSMESS(MSG, STAT)
         ENDIF
         RETURN
      ENDIF

*     --- read data from buf, keep track of missing data
      STAT = 111111
      ISMISS = .FALSE.
      ISINTP = .FALSE.
      ISESTI = .FALSE.

      ATTI = BFATTR(IDAY, NIRRAD)
      IF (ATTI .NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (100000 * (ATTI-1))
      ENDIF
      RIRRAD = BFVALS(IDAY, NIRRAD)

      ATTI = BFATTR(IDAY, NTMIN)
      IF (ATTI .NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (10000 * (ATTI-1))
      ENDIF
      RTMIN = BFVALS(IDAY, NTMIN)

      ATTI = BFATTR(IDAY, NTMAX)
      IF (ATTI.NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (1000 * (ATTI-1))
      ENDIF
      RTMAX = BFVALS(IDAY, NTMAX)

      ATTI = BFATTR(IDAY, NEMPR)
      IF (ATTI.NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (100 * (ATTI-1))
      ENDIF
      REMPR = BFVALS(IDAY, NEMPR)

      ATTI = BFATTR(IDAY, NMWIND)
      IF (ATTI .NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (10 * (ATTI-1))
      ENDIF
      RMWIND = BFVALS(IDAY, NMWIND)

      ATTI = BFATTR(IDAY, NPRECI)
      IF (ATTI .NE. 1) THEN
          IF (ATTI .EQ. DFUNDE) ISMISS = .TRUE.
          IF (ATTI .EQ. DFEST) ISESTI = .TRUE.
          IF (ATTI .EQ. DFINT) ISINTP = .TRUE.
          STAT = STAT + (ATTI-1)
      ENDIF
      RPRECI = BFVALS(IDAY, NPRECI)

*     --- check result, create message if necessary and
*         requested.
*         messages:   missing        (printed if ERRFL on)
*                     no missing:    (printed if WARNFL on)
*                       1 -   estimated (not interpolated)
*                       2 -   estimated and interpolated
*                       3 -   interpolated (not estimated)
      IF (STAT .EQ. 111111) THEN
         STAT = 0
*        --- nothing to report
         RETURN
      ELSEIF (ISMISS) THEN
*     --- missing data
         STAT = -STAT
         IF (ERRFL) THEN
            WRITE(MSG, FMT='(3A,I3,1A,I6)', ERR=990)
     &         'Error in WEATHR: missing data, in: ',
     &         FNAME(1:WSILEN(FNAME)),
     &         ' day:', IDAY, ' attr.:', -(STAT)
            CALL WSMESS(MSG, STAT)
         ENDIF
      ELSE
*     --- warning
         IF (WARNFL) THEN
            MTYPE = ' '
            IF ((ISESTI) .AND. (.NOT. ISINTP)) THEN
               MTYPE = 'Warning in WEATHR: estimated'
            ELSEIF ((ISESTI) .AND. (ISINTP)) THEN
               MTYPE = 'Warning in WEATHR: estimated and interpolated'
            ELSE
               MTYPE = 'Warning in WEATHR: interpolated'
            ENDIF
            WRITE(MSG, FMT='(4A,I3,1A,I6)', ERR=990)
     &          MTYPE(1:WSILEN(MTYPE)), ' data, in:',
     &          FNAME(1:WSILEN(FNAME)), ' day: ', IDAY,
     &         ' attr.:', STAT
            CALL WSMESS(MSG, STAT)
         ENDIF
      ENDIF
*     -----------------------------
      RETURN

990   CONTINUE
      STAT = SSWEAT
      CALL WSMESS ('Internal error in WEATHR: ', STAT )
      RETURN
      END

* ------------------------------------------------------------- *
* wsflgs -- set system flags
* ------------------------------------------------------------- *
* DESCRIPTION
*       Verify and set output flags to NEWFLG, store seperate
*       flag values in common field OF.
*       Initializes itself when called first time.
* ------------------------------------------------------------- *
      LOGICAL FUNCTION WSFLGS(NEWFLG, STAT)

*     IMPLICIT NONE
*     -------- parameters ---------
*     --- in ---
*               new flag
      INTEGER   NEWFLG
*     --- out ---
      INTEGER   STAT
*     --------- external ----------
*     void      WSMESS


*     --------- COMMON: WSNFIL  ------------------------------- *
*                       Weathr Subsystem Num  FILe flags and logf. unit.
*     OF      - Output File flags
*     OFxxxW  - warning flag,   OFxxxF - fatal error flag
*     OFLOGx  - logfile flag,   OFOUTx - output (screen) flag
*     WARNFL  - .TRUE. if to create warings (screen and/or file)
*     ERRFL   - .TRUE. if to create error messages  ( " )

      INTEGER    OFOUTF,OFOUTW,OFLOGF,OFLOGW, NOF
      PARAMETER (OFOUTF = 1,
     &           OFOUTW = 2,
     &           OFLOGF = 3,
     &           OFLOGW = 4,
*                Number Of Flags
     &           NOF = 4)

      LOGICAL    OF(NOF)
*                unit of logfile
      INTEGER    LUNIT
*                error flag from STINFO to WEATHR (12345 = ok)
      INTEGER    STERR
      LOGICAL    WARNFL, ERRFL

      COMMON /WSNFIL/ OF, LUNIT, WARNFL, ERRFL, STERR
*     --------- common wsnfil end ----------------------------- *

*     --------- constants ---------
*               Status value : Parameter FLAG wrong
      INTEGER   SPFLAG
      PARAMETER (
     &          SPFLAG  = -2      )
*     --------- local variables ---
*               copy of flag
      INTEGER   F
*               "value" of flag  (1000,100,10,1)
      INTEGER   FI
*               local copy of flag list (OF), OF is set after
*               every thing is checked, to avoid inconsistencies.
      LOGICAL   LOF(NOF)
*
      INTEGER   I
      CHARACTER MSG*80
      LOGICAL   FIRST

      SAVE
*     ----------- data ------------
      DATA FIRST /.TRUE./

*     ------- start module --------
      WSFLGS = .FALSE.
      F = NEWFLG

      IF (FIRST) THEN
*        --- initialise output flags (warning/error for stdout/file)
         DO 50, I = 1, NOF
            OF(I) = .FALSE.
50       CONTINUE
*        --- send errors to output to start with
         OF(OFOUTF) = .TRUE.
         ERRFL = .TRUE.
         WARNFL = .FALSE.
      ENDIF

      WRITE(MSG, '(A29,I10)' ) 'Error in STINFO: wrong flag: ', NEWFLG
      IF (F .LT. 0) THEN
         STAT = SPFLAG
         CALL WSMESS( MSG, STAT )
         RETURN
      ENDIF

*     --- loop for flags
*     --- flag "1000" is OFLOGW: warnings to log file
*     --- ...
*     --- flag "   1" is OFOUTF: fatals to output
      FI = 1000
      DO 100, I = OFLOGW, OFOUTF, -1
         IF ((F - FI) .GE. 0) THEN
            F = F - FI
            IF (F .GT. FI) THEN
*               --- not 1
                STAT = SPFLAG
                CALL WSMESS( MSG, STAT )
                RETURN
            ENDIF
*           --- enable flag I
            LOF(I) = .TRUE.
         ELSE
*           --- disable flag I
            LOF(I) = .FALSE.
         ENDIF
         FI = FI/10
100   CONTINUE

*     --- New flag ok, make it final.
      DO 500, I=1, NOF
         OF(I) = LOF(I)
500   CONTINUE
      ERRFL  = (OF(OFOUTF) .OR. OF(OFLOGF))
      WARNFL = (OF(OFOUTW) .OR. OF(OFLOGW))

      WSFLGS = .TRUE.
      RETURN

      END

* ------------------------------------------------------------- *
* wsdaop -- open new data file, close current (if any)
* ------------------------------------------------------------- *
* DESCRIPTION
*       Checks parameters, constructs a file name , and opens the
*       file.
* REMARKS
*       Although this routine is written in standard fortran, file-
*       name construction is tricky: correct filename syntax is
*       is system dependend.
* Modified: MAXSTN from 99 to 999
* ------------------------------------------------------------- *

      LOGICAL FUNCTION WSDAOP( IPATH, CNT, STN, YEAR,
     &                         STAT,  DUNIT)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
      CHARACTER IPATH*(*)
      CHARACTER CNT*(*)
      INTEGER   STN
      INTEGER   YEAR
*     --- out ---
      INTEGER   STAT
*               data unit
      INTEGER   DUNIT

*     --------- external ----------
*     viod      WSMESS, WSITOA
      INTEGER   WSILEN
      LOGICAL   WSOPEN

*     -------- constants ----------
*               year is wrong
      INTEGER SPYEAR,
*             country is wrong
     &        SPCNT,
*             station is wrong
     &        SPSTN,
*             can't open data file
     &  SFDATA,
*             internal error in wsdaop
     &  SSDAOP,
*             internal error in wsopen
     &  SSOPEN,
*             maximum size of country
     &  MAXCNT,
*             maximum year
     &  MAXYR,
*             minimum year
     &  MINYR,
*             maximum for station code
     &  MAXSTN
      PARAMETER (
     &      SPYEAR = -3    ,
     &      SPCNT  = -4    ,
     &      SPSTN  = -5    ,
     &      SFDATA = -12   ,
     &      SSOPEN = -903  ,
     &      SSDAOP = -900  ,
     &      MAXCNT = 6     ,
     &      MAXYR  = 1999  ,
     &      MINYR  = 1000  ,
     &      MAXSTN = 999   )


*     --------- COMMON: WSCFIL -------------------------------- *
*                       Weathr Subsystem Char FILe names
*                       note: MAXFNM depends on system.
      INTEGER    MAXFNM
      PARAMETER (MAXFNM = 256)

      CHARACTER*(MAXFNM) PATH, FNAME, LOG
      COMMON   /WSCFIL/  PATH, FNAME, LOG
*     --------- common wscfil end ----------------------------- *


*     ------ local variables ------
      INTEGER            SIZE
      INTEGER            LCNT
      INTEGER            IOSS
      CHARACTER*10       S
      CHARACTER*200      MSG,HMSG
*                               data file name, including path
      CHARACTER*(MAXFNM)        FNDAT, LFNAME, HFNAME

      SAVE
*     ------- start module --------
      WSDAOP = .FALSE.
      LFNAME = ' '

*     ---------------------------
*     --- check parameters
*     --- year is 3 digits only
      IF ((YEAR .LT. MINYR) .OR. (YEAR .GT. MAXYR)) THEN
          MSG= 'Error in STINFO: year: "'
          CALL WSITOA(YEAR, S)
          HMSG = MSG
          MSG = HMSG(1:WSILEN(HMSG))//S(1:WSILEN(S))//
     &          '" is out of range'
          STAT = SPYEAR
          CALL WSMESS( MSG, STAT )
          RETURN
      ENDIF

      LCNT = WSILEN(CNT)

*     --- size of country
      IF ((LCNT .LT. 1) .OR. (LCNT .GT. MAXCNT)) THEN
         STAT = SPCNT
         CALL WSMESS(
     &      'Error in STINFO: wrong string size for country', STAT )
         RETURN
      ENDIF

*     --- station
      IF ((STN .LT. 0) .OR. (STN .GT. MAXSTN)) THEN
         STAT = SPSTN
         CALL WSMESS(
     &      'Error in STINFO: station code out of range', STAT )
         RETURN
      ENDIF
*     ---------------------------


*     ---------------------------
*     --- construct name of file (without path), put in LFNAME
*     --- country part
      LFNAME = CNT(1:LCNT)
      SIZE  = LCNT

*     --- get station number part of filename
      CALL WSITOA(STN, S)
      HFNAME = LFNAME
      LFNAME = HFNAME(1:SIZE)//S
      SIZE = WSILEN(LFNAME)

*     --- year comes in extension part of filename
      HFNAME = LFNAME
      LFNAME = HFNAME(1:SIZE) // '.'
      SIZE = SIZE + 1

*     --- get year
      S = ' '
      WRITE(S, ERR= 900, FMT='(I3.3)') (YEAR - 1000)
      HFNAME = LFNAME
      LFNAME = HFNAME(1:SIZE) // S
      SIZE = SIZE + 3
*     ---------------------------

*     ---------------------------
*     --- add path to name
      IF (WSILEN(IPATH) .NE. 0) THEN
          FNDAT = IPATH(1:WSILEN(IPATH))  // LFNAME(1:SIZE)
      ELSE
*         --- allow for "default directory"
          FNDAT = LFNAME
      ENDIF
      SIZE = WSILEN(FNDAT)
*     ---------------------------

*     ---------------------------
*     --- open file for reading
      IF (.NOT. WSOPEN(FNDAT(1:SIZE), 'r', DUNIT, STAT, IOSS)) THEN
*        --- failed, report if not internal error
         IF (STAT .EQ. SSOPEN) RETURN
         CALL WSITOA(IOSS, S)
         MSG = 'Error in STINFO: cannot open: "'//FNDAT(1:SIZE)//
     &         '" (system status ='//S(1:WSILEN(S))//')'
         STAT = SFDATA
         CALL WSMESS(MSG, STAT)
         RETURN
      ENDIF
*     ---------------------------
*     --- make path/filename global for error messages
      FNAME = LFNAME
      PATH  = IPATH
      WSDAOP = .TRUE.
      RETURN

*-------------------------------------------------------------- *
*
*     --- Internal write errors
900   CONTINUE
      STAT = SSDAOP
      CALL WSMESS( 'Internal error in STINFO', STAT )
      RETURN

      END

* ------------------------------------------------------------- *
* wsattr -- get attributes of data for weather station
* ------------------------------------------------------------- *
* DESCRIPTION
*       This routine is called after the file is opened.
*       Comment lines are read until a data line is seen. This line
*       must contain 5 (no more, no less) REAL values which
*       represent (in  this order): longitude, latitude and
*       altitude and irradiaton conversion factors A and B.
*       If no conversion is needed (data are already in irradiation
*       per square meter) then A and B must be zero.
* ------------------------------------------------------------- *

      LOGICAL FUNCTION  WSATTR(DUNIT,STAT,LON,LAT,ALT,A,B)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
*               data file unit
      INTEGER   DUNIT
*     --- out ---
      INTEGER   STAT
      REAL      LON, LAT, ALT
*               radiation conversion constants
      REAL      A, B
*     ------- constants -----------
*               maximum size of input string
      INTEGER   MAXSTR
*               comment char in column 1
      CHARACTER CHRCOM*1
*               signal comment seen
      INTEGER   SMCOMM
*               signal unexpected end of data file
      INTEGER   SFDEOF
*               signal read error in data file
      INTEGER   SFDERR

      PARAMETER (
     &          CHRCOM = '*' ,
     &          SMCOMM = 1   ,
     &          SFDEOF = -14 ,
     &          SFDERR = -15 ,
     &          MAXSTR = 132 )

*     --------- external ----------
*     void      WSMESS

*     ------ local variables ------
      CHARACTER*(MAXSTR) S
*               local copies of location parms
      REAL      LLON, LLAT, LALT, LA, LB
      SAVE
*     ------- start module --------
      WSATTR = .FALSE.

*     --- loop until location data line is read
100   CONTINUE
*        --- read a string
         READ( UNIT=DUNIT, FMT='(A)', END=900, ERR=910 ) S
         IF(S(1:1) .EQ. CHRCOM) THEN
*            --- print if warn-flag(s) is (are) on
             STAT = SMCOMM
             CALL WSMESS( S, STAT)
*            --- if not equal to old value: error in wsmess.
             IF (STAT .NE. SMCOMM) RETURN
             STAT = 0
             GOTO 100
         ENDIF
*     --- end of loop

*     --- not comment, get data again
      BACKSPACE(UNIT=DUNIT)
      READ(UNIT=DUNIT, FMT=*,END=900,ERR=920) LLON, LLAT, LALT,
     &                                        LA, LB
*     --- read was succesfull, copy to global parms.
      LON = LLON
      LAT = LLAT
      ALT = LALT
      A   = LA
      B   = LB
      WSATTR = .TRUE.
      RETURN

* ------------------------------------------------------------- *
*     read errors:

900   CONTINUE
      STAT = SFDEOF
      CALL WSMESS( 'Error in STINFO: unexpected end of file.', STAT )
      RETURN

910   CONTINUE
      STAT = SFDERR
      CALL WSMESS( 'Error in STINFO: unexpected read error.', STAT )
      RETURN

920   CONTINUE
      STAT = SFDERR
      CALL WSMESS(
     &    'Error in STINFO: incorrect geografical data line', STAT )
      RETURN

      END

* ------------------------------------------------------------- *
* wsrdda -- read data in buffer
* ------------------------------------------------------------- *
* DESCRIPTION
*       WSRDDA reads the actual weathr data, i.e. the second part
*       of the datafile. First the buffer is filled with the
*       value for "undefined". Then the data are read.
*       The actual reading is done in WSGREC (Get RECord) witch
*       returns data values and attributes if an attribute line
*       preceded the data record.
*       After reading the available data, the buffer is scanned
*       for undefined values. These are replaced with an
*       interpolated value (if possible). If data are in hours sun
*       then a conversion to irradiation per square meter is
*       is performed.
* ------------------------------------------------------------- *

      LOGICAL FUNCTION  WSRDDA( DUNIT, LAT, A, B,
     &                          STAT )

*     IMPLICIT NONE

*     --------- parameters --------
*     --- in ---
      INTEGER   DUNIT
*               latitude and conversion factors for hours sun to
*               irradiation.
      REAL      LAT, A, B
*     --- out ---
      INTEGER   STAT
*     --------- constants ---------
      REAL      VALUND
      PARAMETER (
     &          VALUND = -99.       )

*     --- DF: Digit Flags returned in ISTAT
*     --- DFINT : interpolated
*     --- DFUNDE: undefined
      INTEGER    DFOK, DFINT, DFUNDE
      PARAMETER (DFOK = 1,
     &           DFINT = 2,
     &           DFUNDE = 4)



*     --------- COMMON: WSNBUF ------------------------------- *
*                       Weathr Subsystem Num BUFfer
*     --- BFATTR: attributes of data: DFOK, DFINT, DFEST, DFUNDE
*     --- BFVALS: values
      INTEGER    NRDAYS, NRDP
      PARAMETER (
     &           NRDP = 6      ,
     &           NRDAYS = 366  )

      REAL       BFVALS(NRDAYS, NRDP)
      INTEGER    BFATTR(NRDAYS, NRDP)

      COMMON /WSNBUF/ BFVALS, BFATTR
*     --------- common wsnbuf end ----------------------------- *

*     --------- external ----------
*               get NRDP data points from file
      LOGICAL   WSGREC
*               convert hours sun to irradiation
*     void      WSCONI
*     -------- local variables ----
      INTEGER   TODAY
      INTEGER   I
      INTEGER   COL
*               list of data for day
      REAL      DATLST(NRDP)
*               list of attributes for day
      INTEGER   ATTLST(NRDP)
*               attributes set in data file with attribute line.
      INTEGER   ATTSET(NRDP)
      REAL      BVAL
      REAL      SLOPE
*               in missing part of column
      LOGICAL   INMIS
*               start of missing part of column
      INTEGER   START

      SAVE
*     ------- start module --------
*     --- initialize data buffer
      DO 100, TODAY = 1, NRDAYS
         DO 110, I = 1, NRDP
            BFVALS(TODAY, I) = VALUND
            BFATTR(TODAY, I) = DFUNDE
110      CONTINUE
100   CONTINUE
*     --- initialize attribute-set list to default
      DO 120, I=1, NRDP
            ATTSET(I) = DFOK
120   CONTINUE

*     --- fill buffer, loop until eof or all days read
200   CONTINUE
*        --- get a line with NRDP data
         IF (.NOT. WSGREC( DUNIT, ATTSET,
     &                     STAT,  TODAY, DATLST, ATTLST)) GOTO 300
         DO 250, COL=1, NRDP
            BFVALS(TODAY, COL) = DATLST(COL)
            BFATTR(TODAY, COL) = ATTLST(COL)
250      CONTINUE
         GOTO 200

*     --- end read loop ---
300   CONTINUE
*     --- abnormal exit if status unequal to eof.
      IF (STAT .NE. 0) RETURN

*     --- see if interpolations needed. (attributes undefined)
*     --- don't interpolate for rainfall: last collumn
      DO 400, COL=1, NRDP-1
         INMIS = .FALSE.
         DO 410, TODAY=1, NRDAYS

            IF ((BFATTR(TODAY,COL) .EQ. DFUNDE) .AND.
     &           (.NOT. INMIS)) THEN
*              --- beginning of missing part
               INMIS = .TRUE.
               START = TODAY

            ELSEIF (INMIS .AND. (BFATTR(TODAY,COL) .NE. DFUNDE)) THEN
*              --- at end of missing part
*              --- if not from day 1 interpolate
               INMIS = .FALSE.
               IF (START .NE. 1) THEN
*                 --- compute slope
                  SLOPE = (BFVALS(TODAY, COL) - BFVALS(START-1,COL))/
     &                    (TODAY-START+1)
                  BVAL = BFVALS(START-1,COL)
                  DO 411, I=START, TODAY-1
                     BFVALS(I, COL) = BVAL + SLOPE*FLOAT(I-START+1)
                     BFATTR(I, COL) = DFINT
411               CONTINUE
               ENDIF
            ENDIF
410      CONTINUE
*        --- we don't have to check for missing at end (inmis=.true.)
*        --- because values are already set to "undefined"
400   CONTINUE

*     --- if data in hours sun, convert to irradiation
      CALL WSCONI(LAT, A, B)

      WSRDDA = .TRUE.
      RETURN

      END

* ------------------------------------------------------------- *
* wsgrec -- get a record: NDPR data with atributes
* ------------------------------------------------------------- *
* DESCRIPTION
*       Datafiles are devided in to parts: a comment section
*       and a data section. This function reads the data section.
*       Two types of data lines are possible: normal data lines
*       and Attribute lines. Attribute lines are marked by the
*       special station value -999 (year and day are ignored
*       but a dummy value must be present). The parameter ATTSET
*       holds the attributes from the last attribute line. The values
*       in ATTLST may differ from ATTSET because data can be missing.
* PARAMETERS
*       name    class  type   description
*       --------------------------------------------------------
*       dunit   in     int    unit of data file
*       attset  in/out int()  attribute's set from last att. line.
*       stat    out    int    exit status (see below)
*       dayrd   out    int    day read
*       datlst  out    real() array (NRDP elements) with data points
*       attlst  out    int()  array (NRDP elements) with attributes
* RETURNS in STAT
*       0       signals success
*       other   signals failure (message is printed)
* RETURNS
*       .FALSE. when done : STAT = 0
*               on error  : STAT = exit error code
*       .TRUE.  got data, more lines expected.
* ------------------------------------------------------------- *

      LOGICAL FUNCTION  WSGREC( DUNIT,
     &                          ATTSET,
     &                          STAT, DAYRD, DATLST, ATTLST)

*     IMPLICIT NONE

      INTEGER   NRDP
      PARAMETER (
     &          NRDP = 6 )

*     -------- parameters ---------
*     --- in ---
*               datafile unit no.
      INTEGER   DUNIT
*     --- in/out ---
*               attributes as specified in datafile.
*               These are in effect until reset in file,
*               and are not affected by missing data.
      INTEGER   ATTSET(NRDP)
*     --- out ---
*               day read
      INTEGER   DAYRD
*               list of data read
      REAL      DATLST(NRDP)
*               list of attributes for this line.
      INTEGER   ATTLST(NRDP)

*     --------- external ----------
*     void      WSMESS

*     --------- constants ---------
      REAL      VALUND
      INTEGER   SFDERR
      INTEGER   SFDEOF

      PARAMETER (
     &          SFDERR = -15       ,
     &          SFDEOF = -14       ,
     &          VALUND = -99.      )
*     --- DF: Digit Flags returned in ISTAT
*     --- DFUNDE: undefined
      INTEGER    DFUNDE
      PARAMETER (DFUNDE = 4)


*     ------ local variables ------
      INTEGER   STAT, I
      INTEGER   YEARRD, STNRD

      SAVE
*     ------- start module --------
      WSGREC = .FALSE.

      READ(UNIT=DUNIT,FMT=*,END=900,ERR=910)
     &             STNRD, YEARRD, DAYRD, DATLST

*     --- Is it a line with data or attributes?
      IF (STNRD .EQ. -999) THEN
*        --- got line with attributes
         DO 100,I=1, NRDP
            ATTSET(I) = NINT(DATLST(I))
100      CONTINUE
*        --- also get next line: MUST have data for these attributes
         READ(UNIT=DUNIT,FMT=*,END=920,ERR=910)
     &             STNRD, YEARRD, DAYRD, DATLST
      ENDIF

*     --- set attributes
      DO 200, I=1, NRDP
         IF (DATLST(I) .LE. VALUND) THEN
            ATTLST(I) = DFUNDE
         ELSE
            ATTLST(I) = ATTSET(I)
         ENDIF
200   CONTINUE

      WSGREC = .TRUE.
      STAT = 0
      RETURN


* ------------------------------------------------------------- *
*     --- labels for read
900   CONTINUE
*     --- eof
      STAT = 0
      RETURN

910   CONTINUE
      STAT =SFDERR
      CALL WSMESS('Error in STINFO: incorrect data file', STAT )
      RETURN

920   CONTINUE
      STAT = SFDEOF
      CALL WSMESS('Error in STINFO: unexpected end of file', STAT )
      RETURN

      END

*
* ------------------------------------------------------------- *
* wsconi -- convert hours sun to irradation
* ------------------------------------------------------------- *
* PUPOSE
*       convert hours sun measurements to irradiation
* DESCRIPTION
*       WSCONI tests if A and B are not zero the measurements
*       are in hours sun and have to be converted to irradiation
*       per square meter.
* ------------------------------------------------------------- *
      SUBROUTINE WSCONI(ILAT, A, B)

*     IMPLICIT NONE

*     --------- parameters --------
*     --- in ---
      REAL A, B, ILAT
*     --------- constants ---------
*                position of irradiation value in buffer and file
      INTEGER    NIRRAD
      PARAMETER (NIRRAD = 1)

*               PI and conversion factor from degrees to radians
      REAL      PI, RAD
*               undefined data
      REAL      VALUND
      PARAMETER (
     &          PI = 3.141592654  ,
     &          RAD = 0.017453292  ,
     &          VALUND  = -99.0        )
*     --------- external ----------

*     --------- COMMON: WSNBUF ------------------------------- *
*                       Weathr Subsystem Num BUFfer
*     --- BFATTR: attributes of data: DFOK, DFINT, DFEST, DFUNDE
*     --- BFVALS: values
      INTEGER    NRDAYS, NRDP
      PARAMETER (
     &           NRDP = 6      ,
     &           NRDAYS = 366  )

      REAL       BFVALS(NRDAYS, NRDP)
      INTEGER    BFATTR(NRDAYS, NRDP)

      COMMON /WSNBUF/ BFVALS, BFATTR
*     --------- common wsnbuf end ----------------------------- *

*     --------- local -----------
*               length of sun day for particalar latitude.
      REAL      DAYL(NRDAYS)
*               angot for this latitude
      REAL      ANGOT(NRDAYS)
*               last used latitude
      REAL      LAT
      REAL      DAY, DSINB, SC, DEC, SINLD, COSLD, AOB
      INTEGER   IDAY
      SAVE
*     --------- data --------------
*               doesn't exist
      DATA      LAT /100/
*     --------- start module ------

*     --- no conversion needed when less then .01.
*     --- (test for less then .01 instead of 0: this is defacto 0)
      IF ((A .LT. 0.01) .AND. (B .LT. 0.01)) RETURN

      IF (ILAT .NE. LAT) THEN
*     --- new latitude, recompute DAYL and ANGOT.
         LAT = ILAT

         DO 100, IDAY=1, 366

            DAY = FLOAT (IDAY)

*           --- declination of the sun as function of daynumber (DAY)
            DEC = -ASIN (SIN(23.45*RAD)*COS(2.*PI*(DAY+10.)/365.))

*           --- SINLD, COSLD and AOB are intermediate variables
            SINLD = SIN (RAD*LAT)*SIN (DEC)
            COSLD = COS (RAD*LAT)*COS (DEC)
            AOB   = SINLD/COSLD

*           ---  daylength
            DAYL(IDAY) = 12.0*(1.+2.*ASIN (AOB)/PI)

*           ---- integral of sine of solar elevation
            DSINB = 3600.*(DAYL(IDAY)*SINLD+24.*COSLD*
     &                            SQRT(1.-AOB*AOB)/PI)

*           --- solar constant (SC) and daily extraterrestrial radiation (ANGOT)
            SC = 1370.*(1.+0.033*COS(2.*PI*DAY/365.))
            ANGOT(IDAY) = SC*DSINB/1000.

100      CONTINUE
      ENDIF

*     --- replace hours sun with radiation from daylength and angot
      DO 200, IDAY=1, 366
            IF (BFVALS(IDAY,NIRRAD) .NE. VALUND) THEN
               BFVALS(IDAY,NIRRAD) =
     &           ANGOT(IDAY)*(A+B*BFVALS(IDAY,NIRRAD)/DAYL(IDAY))
            ENDIF
200   CONTINUE

      RETURN
      END

* ------------------------------------------------------------- *
* wsmess -- (don't) print message on output and/or logfile
* ------------------------------------------------------------- *
* PARAMETERS
*       name    class   type    description
*       instr   in      char*   message to print
*       stat    in/out  integ   status value
* PURPOSE
*       WSMESS is used for printing (or not if flags are
*       disabled) warnings and/or error messages to the
*       standard output and/or log file.
* REMARK
*       STAT is written if WSMESS exits with an error.
* STATUS CODES
*        lowest highest description
*        -------------------------------------------------
*         <<    -111114 data are missing, print this code
*       -111111      -1 error detected (e.g. wrong parameters)
*        1       1      (SMCOMM): print a comment line
*        2       111111 not used
*        111112  333333 warning, print this code
* REMARKS
*       WSMESS is designed so that it is most efficient
*       if all output is disabled.
* -------------------------------------------------------------- *
      SUBROUTINE WSMESS (INSTR, STAT)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
*                string with message (if error)
      CHARACTER  INSTR*(*)
*     --- in/out ---
*                error-code or warning/missing-code
      INTEGER    STAT
*     --------- external ----------
      INTEGER    WSILEN
      LOGICAL    WSOPEN
*     void       WSITOA

*     --------- COMMON: WSNFIL  ------------------------------- *
*                       Weathr Subsystem Num  FILe flags and logf. unit.
*     OF      - Output File flags
*     OFxxxW  - warning flag,   OFxxxF - fatal error flag
*     OFLOGx  - logfile flag,   OFOUTx - output (screen) flag
*     WARNFL  - .TRUE. if to create warings (screen and/or file)
*     ERRFL   - .TRUE. if to create error messages  ( " )

      INTEGER    OFOUTF,OFOUTW,OFLOGF,OFLOGW, NOF
      PARAMETER (OFOUTF = 1,
     &           OFOUTW = 2,
     &           OFLOGF = 3,
     &           OFLOGW = 4,
*                Number Of Flags
     &           NOF = 4)

      LOGICAL    OF(NOF)
*                unit of logfile
      INTEGER    LUNIT
*                error flag from STINFO to WEATHR (12345 = ok)
      INTEGER    STERR
      LOGICAL    WARNFL, ERRFL

      COMMON /WSNFIL/ OF, LUNIT, WARNFL, ERRFL, STERR
*     --------- common wsnfil end ----------------------------- *
*     --------- COMMON: WSCFIL -------------------------------- *
*                       Weathr Subsystem Char FILe names
*                       note: MAXFNM depends on system.
      INTEGER    MAXFNM
      PARAMETER (MAXFNM = 256)

      CHARACTER*(MAXFNM) PATH, FNAME, LOG
      COMMON   /WSCFIL/  PATH, FNAME, LOG
*     --------- common wscfil end ----------------------------- *

*     --------- constants ---------
*                signal comment
      INTEGER    SMCOMM
*                write to (open) logfile failed
      INTEGER    SFLOGW
*                open of logfile failed: unknown mode.
      INTEGER    SSOPEN
*                unknown status code passed to WSMESS
      INTEGER    SSUSTA
*                maximum size for messages
      INTEGER    MAXMSG
*                special unit nbr: not a unit number (file is closed)
      INTEGER    FSCLSD
*                can't open log file
      INTEGER    SFLOGF
*                STINFO not called, or wrong initialisation
      INTEGER    SOINIT

      PARAMETER (
     &           SOINIT = -21     ,
     &           SFLOGF = -13     ,
     &           FSCLSD = -1      ,
     &           MAXMSG = 100     ,
     &           SSOPEN = -903    ,
     &           SFLOGW = -16     ,
     &           SSUSTA = -904    ,
     &           SMCOMM = 1       )

*     ------ local variables ------
      INTEGER    IOSS
      CHARACTER  MSG*(MAXMSG)
      CHARACTER  S*10
      LOGICAL    ISERR
*                status of WSMESS (used with WSOPEN)
      INTEGER    MSTAT

      SAVE
*     ------- start module --------

*     --- get level of status: warning or fatal-error
      ISERR = (STAT .LT. 0)

      IF (STAT .EQ. SOINIT) THEN
*        --- error caused by NOT calling STINFO, must print to output
         PRINT *, INSTR(1:WSILEN(INSTR))
         RETURN
      ENDIF

*     --- create a message if output/logfile enabled
      IF (ISERR) THEN
         IF (ERRFL) THEN
*            --- error output enabled, create message
             MSG = ' '
             IF (STAT .LT. -111111)  THEN
*               --- missing data
                MSG = INSTR
             ELSEIF (STAT .GT. -900) THEN
*               --- some kind of error detected: report INSTR
                MSG = INSTR
             ELSE
*               --- internal error, unknown status
                PRINT *, 'Internal error in WEATHR/STINFO: ',
     &                   SSUSTA
                STAT = SSUSTA
                RETURN
             ENDIF
         ELSE
*            --- error output disabled, nothing to do
             RETURN
         ENDIF

*     --- warning
      ELSE
         IF (WARNFL) THEN
*           --- warning output enabled
            MSG = ' '
            IF (STAT .GT. 111111) THEN
*              --- interpolated or artificial data
               MSG = INSTR
            ELSEIF (STAT .EQ. SMCOMM) THEN
*              --- comment line read, copy to MSG
               MSG = INSTR
            ELSE
*               --- internal error, unknown status
                PRINT *, 'Internal error in WEATHR/STINFO: ',
     &                   SSUSTA
                STAT = SSUSTA
                RETURN
            ENDIF


         ELSE
*           --- warning output disabled, nothing to do
            RETURN
         ENDIF
      ENDIF

*     --- write message to output and/or logfile
      IF (ISERR) THEN
         IF (OF(OFOUTF)) PRINT *, MSG(1:WSILEN(MSG))
         IF (OF(OFLOGF)) THEN
*           --- open logfile if closed
            IF (LUNIT .EQ. FSCLSD) THEN
               IF (.NOT. WSOPEN(LOG, 'w', LUNIT, MSTAT, IOSS))
     &                   GOTO 920
            ENDIF
            WRITE(UNIT=LUNIT, FMT='(1X,A)', ERR=990, IOSTAT=IOSS)
     &            MSG(1:WSILEN(MSG))
         ENDIF

      ELSE
         IF (OF(OFOUTW)) PRINT *, MSG(1:WSILEN(MSG))
         IF (OF(OFLOGW)) THEN
            IF (LUNIT .EQ. FSCLSD) THEN
               IF (.NOT. WSOPEN(LOG, 'w', LUNIT, MSTAT, IOSS))
     &                   GOTO 920
            ENDIF
            WRITE(UNIT=LUNIT, FMT='(1X,A)', ERR=990, IOSTAT=IOSS)
     &            MSG(1:WSILEN(MSG))
         ENDIF
      ENDIF

      RETURN

* ------------------------------------------------------------- *
*     ----error labels

*     --- log open failed: write to output instead
920   CONTINUE
*     --- If not internal error in WSOPEN print a message and put
*     --- code number in STAT. Internal errors are handled in WSOPEN.
      IF (MSTAT .EQ. SSOPEN) THEN
         STAT = MSTAT
         RETURN
      ENDIF
      STAT = SFLOGF
      CALL WSITOA(IOSS, S)
      MSG = 'Error in WEATHR/STINFO: cannot open logfile:"' //
     &      LOG(1:WSILEN(LOG)) // '"' //
     &      ', (system status =' // S(1:WSILEN(S)) // ')'
      PRINT *, MSG
      RETURN
*     -----------------------------

990   CONTINUE
      STAT = SFLOGW
      CALL WSITOA(IOSS, S)
      MSG = 'Error in WEATHR/STINFO: cannot write to logfile:"'//
     &      LOG(1:WSILEN(LOG)) // '"' //
     &      ', (system status =' // S(1:WSILEN(S)) // ')'
      PRINT *, MSG
      RETURN

      END


* ------------------------------------------------------------- *
* wsopen -- open file NAME with MODE, return unit and status
* ------------------------------------------------------------- *
* DESCRIPTION
*       WSOPEN opens the file NAME with access mode MODE.
*       A unit number is assigned to UNUM for accessing the file.
* MODES
*       The following modes are defined:
*       'r'     open file for reading. On multi user systems
*               shareable files can be read.
*       'w'     the file is opend for writing. If the file exists
*               the old version is deleted first.
* RETURNS
*       The function result is .TRUE. if successfull.
*
*       The parameter STAT is used to return errors caused by
*       wrong arguments. IOSS is used to return errors from the
*       environment, e.g. trying to open a file for reading that
*       doesn't exist. These codes are highly system/compiler
*       dependend.
* ------------------------------------------------------------- *

      LOGICAL FUNCTION WSOPEN(NAME, MODE,
     &                        FUNIT, STAT, IOSS)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
      CHARACTER    NAME*(*)
*                  valid modes are: 'r' = readonly; 'w' = writeonly
      CHARACTER    MODE*1
*     --- out ---
*                  file handle
      INTEGER      FUNIT
*                  exit status
      INTEGER      STAT
*                  iostatus returned by system
      INTEGER      IOSS
*     --------- external ----------
      INTEGER      WSILEN
*     void         WSFUN
*     ------ local variables ------
      INTEGER      FSCLSD
      INTEGER      SSOPEN
      PARAMETER (
     &          SSOPEN = -903     ,
     &          FSCLSD = -1       )
*               size of name
      INTEGER   LNAME
      LOGICAL   EXST
      SAVE
*     ------- start module --------
      WSOPEN = .FALSE.
      IOSS = 0
      STAT = 0

*     --- get unit number
      CALL WSFUN(FUNIT)
      LNAME = WSILEN(NAME)

*     --- --- read mode -------------------------------
*     --- open readonly, sequential, existing (old) file
      IF (MODE(1:1) .EQ. 'r') THEN
         OPEN(UNIT=FUNIT,
     &        STATUS='OLD',
     &        IOSTAT=IOSS,
     &        FILE=NAME(1:LNAME))
         IF (IOSS .NE. 0) THEN
            FUNIT = FSCLSD
            RETURN
         ELSE
            WSOPEN = .TRUE.
            RETURN
         ENDIF

*     --- ---- write mode ------------------------
*     --- open new file for sequential writing
      ELSEIF (MODE(1:1) .EQ. 'w') THEN
*        --- to be safe: first delete file if it exists
         INQUIRE(FILE=NAME(1:LNAME), EXIST=EXST)
         IF (EXST) THEN
            OPEN(UNIT=FUNIT,
     &           FILE=NAME(1:LNAME),
     &           STATUS='OLD',
     &           IOSTAT=IOSS )
            CLOSE(UNIT=FUNIT, STATUS='DELETE')
         ENDIF
*        --- now we can open with 'NEW'
         OPEN(UNIT=FUNIT,
     &        FILE = NAME(1:LNAME),
     &        STATUS = 'NEW',
     &        IOSTAT = IOSS )
         IF (IOSS .NE. 0) THEN
            FUNIT = FSCLSD
            RETURN
         ELSE
            WSOPEN = .TRUE.
            RETURN
         ENDIF

*     --- ---- cannot happen ---------------------
      ELSE
           STAT = SSOPEN
           FUNIT = FSCLSD
*          --- must print to output (can't call wsmess)
           PRINT *,
     &       'Error in WEATHER/STINFO: internal error, code:', STAT
           RETURN
      ENDIF
      END

* ------------------------------------------------------------- *
* wsilen -- return significant length of a string
* ------------------------------------------------------------- *
* DESCRIPTION
*       The string is searched from end to begin for the first
*       non-space.
*       If the string is empty 0 is returned.
* HISTORY
*       Author: Daniel van Kraalingen, orginal name: ILEN.
*       Date: Aug 87
* ------------------------------------------------------------- *
      INTEGER FUNCTION WSILEN (STRING)

      CHARACTER*(*) STRING
      SAVE

      DO 10 WSILEN=LEN(STRING),1,-1
         IF (STRING(WSILEN:WSILEN).NE.' ') RETURN
10    CONTINUE

      RETURN
      END

* ------------------------------------------------------------- *
* wsstar -- return first significant character of string
* ------------------------------------------------------------- *
* DESCRIPTION
*       The string is searched from begin to end for the first
*       non-space.
*       If the string is empty (spaces only) 0 is returned.
* HISTORY
*       Author: Daniel van Kraalingen, orginal name: ISTART.
*       Date: Aug 87
* ------------------------------------------------------------- *
      INTEGER FUNCTION WSSTAR (STRING)

      CHARACTER*(*) STRING

      SAVE

      DO 10 WSSTAR=1,LEN(STRING)
         IF (STRING(WSSTAR:WSSTAR).NE.' ') RETURN
10    CONTINUE

      WSSTAR = 0

      RETURN
      END

* ------------------------------------------------------------- *
* wsfun -- get free unit number
* ------------------------------------------------------------- *
* DESCRIPTION
*       wsfun scans for unused unitnumbers, starting with 92
*       down to 40.
*       The numbers are checked with INQUIRE, so that numbers
*       are returned to the "free-list" by using standard CLOSE.
* RETURNS
*        -1 in FD if no valid unit number found.
* BUGS
*       Inconsequent use causes subtle bugs.
* ------------------------------------------------------------- *
      SUBROUTINE WSFUN(FD)

*     IMPLICIT NONE

*     -------- parameters ---------
      INTEGER   FD

*     ------- constants -----------
*               lowest,highest unit number opened:
      INTEGER   LOWUNT, MAXUNT
      PARAMETER (LOWUNT = 40,
     &           MAXUNT = 92)

*     ------ local variables ------
      LOGICAL   ISOPEN
      SAVE
*     ------- start module --------

*     --- fine unused unit number
      DO 10, FD = MAXUNT, LOWUNT, -1
         INQUIRE(UNIT=FD, OPENED=ISOPEN)
         IF (.NOT. ISOPEN) RETURN
10    CONTINUE

      FD = -1
      RETURN

      END

* ------------------------------------------------------------- *
* wsitoa -- convert integer to (left adjusted) ascii in string  *
* ------------------------------------------------------------- *
* DESCRIPTION                                                   *
*       The integer I is written in STRING starting at the      *
*       first postion.                                          *
* ------------------------------------------------------------- *
      SUBROUTINE WSITOA (I, STRING)

*     IMPLICIT NONE

*     -------- parameters ---------
*     --- in ---
      INTEGER   I
*     --- in/out ---
      CHARACTER STRING*(*)
*     --------- external ----------
      INTEGER   WSILEN, WSSTAR

*     ------ local variables ------
      CHARACTER S*80
      INTEGER   IS
      INTEGER   IE
      SAVE
*     ------- start module --------

      WRITE(S, '(I20)' )  I

      IS = WSSTAR(S)
      IE = WSILEN(S)

      STRING = S(IS:IE)
      RETURN

      END
