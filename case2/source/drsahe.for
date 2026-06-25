*----------------------------------------------------------------------*
* SUBROUTINE DRSAHE                                                    *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Author : Daniel van Kraalingen                                       *
* Modifications: Pieter Zuidema (data retrieval from soil.dat)         *
* Date   : August 1994, modification 2002                              *
* Version: 1.4                                                         *
*                                                                      *
* Purpose: Tipping bucket water balance routine                        *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* control                                                              *
* ITASK   I4  determines action of routine                  -     C,I  *
* IUNIT   I4  unit number to be used, see file usage below  -   IN,C,I *
* IUNLOG  I4  unit number in use for LOG FILE               -   IN,C,I *
*             = 0, no log file is used or assumed to exist             *
*             > 0, error messages are written to log file              *
* FILIN   C*  name of file with soil data                   -   IN,C,I *
* NLXM    I4  no. of layers as declared in calling program  -     IN   *
* SOILTYPE I4 type of soil as specified in soil.dat         -      I   *
*                                                                      *
* time variables                                                       *
* TIME    R4  simulation time                               d     T,I  *
* DELT    R4  time step                                     d     T,I  *
* OUTPUT  L4  Flag to indicate if output should be done     -      I   *
*                                                                      *
* dynamic input                                                        *
* EVSC    R4  potential evaporation rate                   mm/d    I   *
* RAIN    R4  rainfall / irrigation rate                   mm/d    I   *
* TRWL    R4  actual transpiration rate per layer          mm/d    I   *
*                                                                      *
* soil description (available after initial call)                      *
* NL      I4  number of layers specified in input file      -      O   *
* TKLX    R4  thickness of soil compartments                m    IN,O  *
* ZRTMS   R4  maximum rooting depth as soil characteristic  -      O   *
* WCADX   R4  volumetric water content airdry               -    IN,O  *
* WCWPX   R4  volumetric water content at wilting point     -    IN,O  *
* WCFCX   R4  volumetric water content at field capacity    -    IN,O  *
* WCSTX   R4  volumetric water content at saturation        -    IN,O  *
*                                                                      *
* dynamic output                                                       *
* EVSW    R4  actual (realized) evaporation rate           mm/d    O   *
* FLXQT   R4  layer boundary fluxes (rates)                mm/d    O   *
* WCLQT   R4  volumetric soil water content per layer       -      O   *
*                                                                      *
* cumulated, derived and help variables                                *
* DRAICU  R4  cumulative drainage by drains                 mm     O   *
* EVSWCU  R4  cumulative evaporation                        mm     O   *
* INFCU   R4  cumulative infiltration                       mm     O   *
* TRWCU   R4  cumulative transpiration                      mm     O   *
* FLXCU   R4  cumulative flux for each layer boundary       mm     O   *
*                                                                      *
* Fatal error checks: many, see docs                                   *
* Warnings          : many, see docs                                   *
* Subprograms called: many from TTUTIL                                 *
* File usage        : IUNIT, IUNIT+1, IUNLOG                           *
*----------------------------------------------------------------------*

      SUBROUTINE DRSAHE (ITASK , IUNIT , IUNLOG, FILIN , SOILTYPE,
     &                   IDOY  , IYEAR , DELT  , OUTPUT,
     &                   NLXM  , NL    , EVSC  , RAIN  , TRWL  ,
     &                   TKLX  , ZRTMS ,
     &                   WCADX , WCWPX , WCFCX , WCSTX,
     &                   EVSW  , FLXQT , WCLQT ,
     &                   DRAICU, EVSWCU, INFCU , TRWCU, FLXCU)
      IMPLICIT NONE

*     Formal parameters
      INTEGER ITASK, IUNIT, IUNLOG, IDOY, IYEAR, SOILTYPE
      CHARACTER FILIN*(*)
      LOGICAL OUTPUT

      REAL DELT  , EVSC  , RAIN , EVSW
      REAL DRAICU, EVSWCU, INFCU, TRWCU, ZRTMS

      INTEGER NLXM, NL
      REAL TRWL(NLXM) , TKLX(NLXM) , WCADX(NLXM)  , WCWPX(NLXM)
      REAL WCFCX(NLXM), WCSTX(NLXM), FLXQT(NLXM+1), WCLQT(NLXM)
      REAL FLXCU(NLXM+1)

*     Local variables
      REAL RAINCU,RNOFCU,RNOFF, INF, WEFF 

*     Variables used to read soil.dat file
      CHARACTER STRING*80
	INTEGER   SIGLEN

*     Soil description arrays
      INTEGER NLLM
      PARAMETER (NLLM=10)
      INTEGER ITYL(NLLM)
      REAL TKL(NLLM)   , TYL(NLLM)  , DEPTH(NLLM)
      REAL WCAD(NLLM)  , WCWP(NLLM) , WCFC(NLLM) , WCST(NLLM)
      REAL WCLQTM(NLLM), WCLCH(NLLM), WCL(NLLM)

*     Array with irrigation data
      INTEGER IMNID, INID
      PARAMETER (IMNID=100)
      REAL IRRTB(IMNID), IRR, IRRCU
      INTEGER YRO,DYO

*     Soil characteristics according to Rijtema/Driessen system
*     The number of soil types defined is NRDTYP
      INTEGER NRDTYP
      PARAMETER (NRDTYP=20)
      REAL MSWCAT(NRDTYP), WCSTT(NRDTYP), MSWCA(NLLM)

*     Soil characteristics according to Van Genuchten system
*     The number of soil types defined is NVGTYP
      INTEGER NVGTYP
      PARAMETER (NVGTYP=2)
      REAL VGWRT(NVGTYP), VGWST(NVGTYP), VGAT(NVGTYP), VGNT(NVGTYP)
      REAL VGA(NLLM)    , VGR(NLLM)    , VGN(NLLM)
      REAL VGM, HLP1, WREL

*     Linear interpolation on user-defined log scale
      REAL PFWC00(NLLM), PFWC01(NLLM), PFWC02(NLLM), PFWC03(NLLM)
      REAL PFWC04(NLLM), PFWC05(NLLM), PFWC06(NLLM), PFWC07(NLLM)
      REAL PFWC08(NLLM), PFWC09(NLLM), PFWC10(NLLM)
      REAL PF(22)

*     Parameters for field capacity, wilting point and airdry
      REAL FIELD, WILTP, AIRDR
      PARAMETER (FIELD = 1.0E2, WILTP = 1.6E4, AIRDR = 1.0E7)

*     Functions
      LOGICAL RDINQR
      REAL LINT, INTGRL

*     Control, switch, temporary and miscellaneous variables
      INTEGER IL, ITMP, DATCMP, I1, I2
      INTEGER SWIT6, SWIT8, SWIT9
      REAL WCUM, WCUMO, TRW, FLOW, CAP, DRAIQT, WCUMCH, CHECK
      REAL TRCH, EVSCL, WCTMP
      REAL VAR(NLLM), RESOIL(NLLM)
      REAL SUM, EES, EVSW2, EVSH, EVSD, RDSLR, DSLR

      SAVE

*     Soil type properties according to the Rijtema/Driessen combination.
*     Data from these tables will be used at the simultaneous occurrence
*     of the following switch values:

*     The 6 soil properties are given in data tables for 20 soils.
*     By defining the soil type TYL(I) for each layer I, a consistent
*     combination of soil properties is selected. Data refer to the
*     twenty standard soil types according to Rijtema (as described
*     by Driessen, 1986).
*
*      1. Coarse sand                11. Fine sandy loam
*      2. Medium coarse sand (mcs)   12. Silt loam
*      3. Medium fine sand           13. Loam
*      4. Fine sand                  14. Sandy clay loam
*      5. Humous loamy mcs           15. Silty clay loam
*      6. Light loamy mcs            16. Clay loam
*      7. Loamy mcs                  17. Light clay
*      8. loamy fine sand            18. Silty clay
*      9. Sandy loam                 19. Heavy clay
*     10. Loess loam                 20. Peat

      DATA MSWCAT  /0.0853, 0.0450, 0.0366, 0.0255, 0.0135,
     &              0.0153, 0.0243, 0.0299, 0.0251, 0.0156,
     &              0.0186, 0.0165, 0.0164, 0.0101, 0.0108,
     &              0.0051, 0.0085, 0.0059, 0.0043, 0.0108/

*     saturated soil moisture content, dimensionless (Rijtema/Driessen)
      DATA WCSTT   /0.3950, 0.3650, 0.3500, 0.3640, 0.4700,
     &              0.3940, 0.3010, 0.4390, 0.4650, 0.4650,
     &              0.5040, 0.5090, 0.5030, 0.4320, 0.4750,
     &              0.4450, 0.4530, 0.5070, 0.5400, 0.8630/

*     Soil type properties according to the van Genuchten system.
*     The 6 soil properties are given in data tables for NVGTYP soils.
*     (for declarations, see above).  By defining the soil type TYL(I)
*     for each layer I, a consistent combination of soil properties
*     is selected. In the case presented here, data refer to
*     a few soil types of the Lovinkhoeve (Peter de Willigen)
*     To be extended by user.
*
*      1. Lovinkhoeve 12b
*      2. Lovinkhoeve 16a

*     TETA-r, dimensionless
      DATA VGWRT / 0.0448,  0.0000/
*     TETA-s, dimensionless
      DATA VGWST / 0.4012,  0.4505/
*     ALPHA in cm-1
      DATA VGAT  / 0.0036,  0.0067/
*     N, dimensionless
      DATA VGNT  / 1.5007,  1.2318/

      IF (DELT.GT.1.) CALL FATALERR ('DRSAHE','delt too large')

      IF (ITASK.EQ.1) THEN

*        ----------------------
*        Initialization section
*        ----------------------

*        Read input file
         CALL RDINIT (IUNIT, IUNLOG, FILIN)

*        Read number of layers and thickness for one of the
*        3 standard Driessen soils specified or 5 user-defined
*        soil types in soil.dat. (Added PAZ, 1-2002) 
*        read number of soil layers
         STRING = 'NL'
 	   SIGLEN = LEN_TRIM(STRING)
         CALL ADDINT(STRING,SIGLEN,SOILTYPE)
         CALL RDSINT (STRING, NL)
         IF (NL.EQ.0) CALL FATALERR 
     &   ('DRSAHE','No information for this SOILTYPE in soil.dat')
*        read number of thickness of soil layers
         STRING = 'TKL'
 	   SIGLEN = LEN_TRIM(STRING)
         CALL ADDINT(STRING,SIGLEN,SOILTYPE)
         CALL RDFREA (STRING, TKL, NLLM, NL)

         IF (NL.GT.NLLM) CALL FATALERR
     &      ('DRSAHE','too many layers defined in data file')
         IF (NLXM.LT.NL) CALL FATALERR
     &      ('DRSAHE','too few layers in external arrays')

*        Read evaporation proportionality factor and switches
*        Note: for CASE2 2.2, the standard values for SWIT8 = 1
*        and for SWIT9 = 2 (Driessen type soil). 
         CALL RDSREA ('EES', EES)
         CALL RDSINT ('SWIT9', SWIT9)
         CALL RDSINT ('SWIT8', SWIT8)

         IF (SWIT9.EQ.1) THEN
*           Moisture characteristics by user-defined parameters

            CALL RDFREA ('WCST', WCST, NLLM, NL)

            IF (SWIT8.EQ.1) THEN

*              Driessen moisture characteristic
               CALL RDFREA ('MSWCA', MSWCA, NLLM, NL)
               DO 10 IL=1,NL
                  WCFC(IL) = WCST(IL)*EXP (-MSWCA(IL)*LOG (FIELD)**2)
                  WCWP(IL) = WCST(IL)*EXP (-MSWCA(IL)*LOG (WILTP)**2)
                  WCAD(IL) = WCST(IL)*EXP (-MSWCA(IL)*LOG (AIRDR)**2)
10             CONTINUE

            ELSE IF (SWIT8.EQ.2) THEN

*              Van Genuchten moisture characteristic
               CALL RDFREA ('VGA', VGA, NLLM, NL)
               CALL RDFREA ('VGR', VGR, NLLM, NL)
               CALL RDFREA ('VGN', VGN, NLLM, NL)

               DO 20 IL=1,NL
                  VGM = 1.-1./VGN(IL)

                  HLP1     = (FIELD*VGA(IL))**VGN(IL)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCFC(IL) = WREL*(WCST(IL)-VGR(IL))+VGR(IL)

                  HLP1     = (WILTP*VGA(IL))**VGN(IL)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCWP(IL) = WREL*(WCST(IL)-VGR(IL))+VGR(IL)

                  HLP1     = (AIRDR*VGA(IL))**VGN(IL)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCAD(IL) = WREL*(WCST(IL)-VGR(IL))+VGR(IL)
20             CONTINUE

            ELSE IF (SWIT8.EQ.3) THEN

*              Linear interpolation on user-defined log scale

*              Read pF values
               CALL RDFREA ('PFWC00', PFWC00, NLLM, NL)
               CALL RDFREA ('PFWC01', PFWC01, NLLM, NL)
               CALL RDFREA ('PFWC02', PFWC02, NLLM, NL)
               CALL RDFREA ('PFWC03', PFWC03, NLLM, NL)
               CALL RDFREA ('PFWC04', PFWC04, NLLM, NL)
               CALL RDFREA ('PFWC05', PFWC05, NLLM, NL)
               CALL RDFREA ('PFWC06', PFWC06, NLLM, NL)
               CALL RDFREA ('PFWC07', PFWC07, NLLM, NL)
               CALL RDFREA ('PFWC08', PFWC08, NLLM, NL)
               CALL RDFREA ('PFWC09', PFWC09, NLLM, NL)
               CALL RDFREA ('PFWC10', PFWC10, NLLM, NL)

*              Set up relative moisture content values
               PF(2)  = 0.
               PF(4)  = 0.1
               PF(6)  = 0.2
               PF(8)  = 0.3
               PF(10) = 0.4
               PF(12) = 0.5
               PF(14) = 0.6
               PF(16) = 0.7
               PF(18) = 0.8
               PF(20) = 0.9
               PF(22) = 1.0

*              Fill array with pf values for subsequent soil layers
               DO 30 IL=1,NL
                  PF(1)  = PFWC00(IL)
                  PF(3)  = PFWC01(IL)
                  PF(5)  = PFWC02(IL)
                  PF(7)  = PFWC03(IL)
                  PF(9)  = PFWC04(IL)
                  PF(11) = PFWC05(IL)
                  PF(13) = PFWC06(IL)
                  PF(15) = PFWC07(IL)
                  PF(17) = PFWC08(IL)
                  PF(19) = PFWC09(IL)
                  PF(21) = PFWC10(IL)
                  WCFC(IL) = WCST(IL)*MAX (0.01, LINT (PF,22,FIELD))
                  WCWP(IL) = WCST(IL)*MAX (0.01, LINT (PF,22,WILTP))
                  WCAD(IL) = WCST(IL)*MAX (0.01, LINT (PF,22,AIRDR))
30             CONTINUE

            ELSE IF (SWIT8.EQ.4) THEN

*              User must specify pf-curve parameters to be read
*              and include error check
               CALL RDFREA ('WCFC', WCFC, NLLM, NL)
               CALL RDFREA ('WCWP', WCWP, NLLM, NL)
               CALL RDFREA ('WCAD', WCAD, NLLM, NL)

            ELSE
               CALL FATALERR ('DRSAHE','Illegal SWIT8 value')
            END IF


         ELSE IF (SWIT9.EQ.2) THEN

*           Physical properties from soil type number
*           Read soil type numbers depending on the value of SOILTYPE
*           (Added PAZ 1-2002)
*           read texture type of soil layers
            STRING = 'TYL'
 	      SIGLEN = LEN_TRIM(STRING)
            CALL ADDINT(STRING,SIGLEN,SOILTYPE)
            CALL RDFREA (STRING, TYL, NLLM, NL)

*           CALL RDFREA ('TYL', TYL, NLLM, NL)

            IF (SWIT8.EQ.1) THEN
*              Driessen moisture characteristic
               DO 40 IL=1,NL
                 ITYL(IL) = NINT (TYL(IL))
                 ITMP = ITYL(IL)
                 WCST(IL) =WCSTT(ITMP)
                 WCFC(IL) =WCST(IL)*EXP(-MSWCAT(ITMP)*LOG(FIELD)**2)
                 WCWP(IL) =WCST(IL)*EXP(-MSWCAT(ITMP)*LOG(WILTP)**2)*1.1
                 WCAD(IL) =WCST(IL)*EXP(-MSWCAT(ITMP)*LOG(AIRDR)**2)
40             CONTINUE

            ELSE IF (SWIT8.EQ.2) THEN
*              Van Genuchten moisture characteristic
               DO 50 IL=1,NL
                  ITYL(IL) = NINT (TYL(IL))
                  ITMP = ITYL(IL)

                  WCST(IL) = VGWST(ITMP)
                  VGM      = 1.-1./VGNT(ITMP)

                  HLP1     = (FIELD*VGAT(ITMP))**VGNT(ITMP)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCFC(IL) = WREL*(WCST(IL)-VGWRT(ITMP))+VGWRT(ITMP)

                  HLP1     = (WILTP*VGAT(ITMP))**VGNT(ITMP)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCWP(IL) = WREL*(WCST(IL)-VGWRT(ITMP))+VGWRT(ITMP)

                  HLP1     = (AIRDR*VGAT(ITMP))**VGNT(ITMP)
                  WREL     = (1.+HLP1)**(-VGM)
                  WCAD(IL) = WREL*(WCST(IL)-VGWRT(ITMP))+VGWRT(ITMP)

50             CONTINUE
            ELSE
               CALL FATALERR
     &            ('DRSAHE','SWIT8 wrong value ; should be 1 or 2')
            END IF
         ELSE

           CALL FATALERR 
     &            ('DRSAHE','SWIT9 wrong value ; should be 1 or 2')

         END IF

*        Initial water contents
         CALL RDSINT ('SWIT6', SWIT6)

         IF (SWIT6.EQ.1) THEN

*           In hydrostatic equilibrium

            WRITE (*,'(2A,/,A)')
     &        ' WARNING from DRSAHE: initial soil moisture',
     &        ' in hydrostatic equilibrium',
     &        ' not implemented. Instead, field capacity is used.'
            IF (IUNLOG.GT.0) WRITE (IUNLOG,'(2A,/,A)')
     &        ' WARNING from DRSAHE: initial soil moisture',
     &        ' in hydrostatic equilibrium',
     &        ' not implemented. Instead, field capacity is used.'

            DO 60 IL=1,NL
               WCLQTM(IL) = WCFC(IL)
60          CONTINUE

         ELSE IF (SWIT6.EQ.2) THEN

*           At observed moisture contents

            CALL RDFREA ('WCLQTM', WCLQTM, NLLM, NL)

            DO 70 IL=1,NL
               IF (WCLQTM(IL).GT.WCFC(IL)) THEN
                  WRITE (*,'(2A)')
     &              ' WARNING from DRSAHE: initial soil moisture',
     &              ' content larger than field capacity'
                  IF (IUNLOG.GT.0) WRITE (IUNLOG,'(2A)')
     &              ' WARNING from DRSAHE: initial soil moisture',
     &              ' content larger than field capacity'
                  WCLQTM(IL) = WCFC(IL)
               ELSE IF (WCLQTM(IL).LT.WCAD(IL)) THEN
                  WRITE (*,'(2A)')
     &              ' WARNING from DRSAHE: initial soil moisture',
     &              ' content less than air dry'
                  IF (IUNLOG.GT.0) WRITE (IUNLOG,'(2A)')
     &              ' WARNING from DRSAHE: initial soil moisture',
     &              ' content less than air dry'
                  WCLQTM(IL) = WCAD(IL)
               END IF
70          CONTINUE

         ELSE IF (SWIT6.EQ.3) THEN

*           At wilting point

            DO 80 IL=1,NL
               WCLQTM(IL) = WCWP(IL)
80          CONTINUE

         ELSE

            CALL FATALERR ('DRSAHE',
     &       'SWIT6 wrong value ; should be 1, 2 or 3')

         END IF

*        Read irrigation data
         IF (RDINQR ('IRRTB')) THEN
            CALL RDAREA ('IRRTB', IRRTB, IMNID, INID)
            IF (MOD (INID,3).NE.0) CALL FATALERR ('DRSAHE',
     &         'number of irrigation data is not multiple of three')
         ELSE
*           No irrigation data
            INID = 0
         END IF

*        End of reading from data file
         CLOSE (IUNIT)

*        Calculate array with depths
         DEPTH(1) = 0.5*TKL(1)
         DO 90 IL=2,NL
            DEPTH(IL) = DEPTH(IL-1)+0.5*TKL(IL-1)+0.5*TKL(IL)
90       CONTINUE

*        Maximum rooting depth as soil characteristic
         ZRTMS = DEPTH(NL)+0.5*TKL(NL)

*        Initialize remaining variables
         WCUM  = 0.
         DO 100 IL=1,NL
            WCLQT(IL) = WCLQTM(IL)
            WCUM      = WCUM+WCLQT(IL)*TKL(IL)*1000.
            FLXCU(IL) = 0.
            FLXQT(IL) = 0.
100      CONTINUE

         FLXCU(NL+1) = 0.
         FLXQT(NL+1) = 0.

         RAINCU = 0.
         IRRCU  = 0.
         RNOFCU = 0.
         INFCU  = 0.
         EVSWCU = 0.
         TRWCU  = 0.
         DRAICU = 0.

         EVSW   = 0.
         DSLR   = 1.

*        Copy soil description arrays to external arrays
         DO 110 IL=1,NL
            TKLX(IL)  = TKL(IL)
            WCADX(IL) = WCAD(IL)
            WCWPX(IL) = WCWP(IL)
            WCFCX(IL) = WCFC(IL)
            WCSTX(IL) = WCST(IL)
110      CONTINUE

*        Set not used elements to zero
         DO 120 IL=NL+1,NLXM
            TKLX(IL)  = 0.
            WCADX(IL) = 0.
            WCWPX(IL) = 0.
            WCFCX(IL) = 0.
            WCSTX(IL) = 0.
120      CONTINUE

         YRO    = 0
         DYO    = 0

      ELSE IF (ITASK.EQ.2) THEN

*        ------------------------
*        Rate calculation section
*        ------------------------

*        Determine rates of change of water balance

*        Check: number of layers for external arrays is large
*        enough to hold data
         IF (NLXM.LT.NL) THEN
            CALL FATALERR ('DRSAHE','too few layers in external arrays')
         END IF

*        Check: evaporation should be positive
         IF (EVSC.LT.0.) THEN
            WRITE (*,'(2A,/,2A)')
     &        ' WARNING from DRSAHE: potential soil evaporation',
     &        ' has negative sign !',
     &        '   To extract water from the soil,',
     &        ' the sign should be positive.'
            EVSC = 0.
         END IF

         EVSCL = -EVSC

*        Check: rainfall should be positive
         IF (RAIN.LT.0.) THEN
            WRITE (*,'(2A,/,2A)')
     &        ' WARNING from DRSAHE: rainfall',
     &        ' has negative sign !',
     &        '   To add water to the soil through rainfall,',
     &        ' the sign should be positive.'
            RAIN = 0.
         END IF

*        Check: transpiration should be negative and 'available'
         DO 190 IL=1,NL
            IF (TRWL(IL).LT.0.) THEN
               WRITE (*,'(2A,/,2A)')
     &           ' WARNING from DRSAHE: transpiration',
     &           ' has negative sign !',
     &           '   To extract water from the soil,',
     &           ' the sign should be positive.'
               TRWL(IL) = 0.
            END IF
190      CONTINUE

*        Set rates of change to zero and make local water status
*        array equal to current water status, reset fluxes
         DO 220 IL=1,NL
            WCLCH(IL) = 0.
            WCL(IL)   = WCLQT(IL)
            FLXQT(IL) = 0.
220      CONTINUE
         FLXQT(NL+1) = 0.

*        Cumulate transpiration
         TRW = 0.
         DO 230 IL=1,NL
            TRW = TRW-TRWL(IL)
230      CONTINUE

*        Effectuate transpiration on local status array
*        (transpiration has a positive value)
         DO 240 IL=1,NL
            TRCH = DELT*TRWL(IL)/(TKL(IL)*1000.)
            IF (WCL(IL).GE.WCWP(IL).AND.(WCL(IL)-TRCH).LT.WCWP(IL)) THEN
               WRITE (*,'(A,/,2A,I3)') ' WARNING from DRSAHE:',
     &           '   Transpiration extracted water below wilting point',
     &           ' in layer:',IL
            END IF
            WCL(IL) = WCL(IL)-TRCH
240      CONTINUE

*        Calculate irrigation application
         IRR = 0.
         IF (INID.GT.0) THEN
*           Irrigation table was present, check table if irrigation
*           should be applied by trying to find an exact date match

            I1 = -1
            I2 = 1
245         IF (I1.NE.0.AND.I2.LE.INID-2) THEN
               I1 = DATCMP (NINT(IRRTB(I2)),NINT(IRRTB(I2+1)),
     &                      IYEAR,IDOY)
               IF (I1.NE.0) I2 = I2+3
            GOTO 245
            END IF

            IF (I1.EQ.0) THEN
*              Exact date match found, are year or day different
*              from previous date (if day and year are equal to year
*              and day from the previous step, delt is smaller than 1
*              and the irrigation should not be reapplied)

               IF (IYEAR.NE.YRO.OR.IDOY.NE.DYO) IRR = IRRTB(I2+2)/DELT
            END IF

            YRO = IYEAR
            DYO = IDOY
         END IF

         RNOFF = MAX (0., 0.15*(RAIN-10.))
         INF   = IRR+RAIN-RNOFF

         EVSH  = MAX (EVSCL, -((WCL(1)-WCAD(1))*TKL(1)*1000.)/DELT-INF)
         EVSD  = MAX (EVSCL, 0.6*EVSCL*(SQRT (DSLR+1.)-SQRT (DSLR))-INF)

         IF (INF.GT.0.5) THEN
            EVSW2  = EVSH
            RDSLR = -(DSLR-1.)/DELT
         ELSE
            EVSW2  = EVSD
            RDSLR = 1.
         END IF

*        Calculate array for exponential extinction of evaporation
         SUM = 0.
         DO 250 IL=1,NL
            VAR(IL) = TKL(IL)*1000.*(WCL(IL)-WCAD(IL))*
     &                EXP (-EES*(DEPTH(IL)-0.25*TKL(IL)))
            SUM     = SUM+VAR(IL)
250      CONTINUE

*        Effectuate evaporation on local status array and calculate
*        the actual soil evaporation
         EVSW = 0.
         DO 260 IL=1,NL
            IF (SUM.GT.0.) THEN
*              Water available somewhere in profile
               RESOIL(IL) = EVSW2*VAR(IL)/SUM
               WCTMP      = MAX (WCAD(IL),
     &                      WCL(IL)+DELT*RESOIL(IL)/(TKL(IL)*1000.))
               RESOIL(IL) = (WCTMP-WCL(IL))*TKL(IL)*1000./DELT
               RESOIL(IL) = MIN (0.,RESOIL(IL))
               EVSW       = EVSW+RESOIL(IL)
               WCL(IL)    = WCTMP
            ELSE
*              Water not available in profile
               RESOIL(IL) = 0.
            END IF
260      CONTINUE


*        Effectuate infiltration on local status array
*        nog een keer naar delt kijken
         FLOW = 0.
         IF (INF.GT.0.) THEN
            FLOW     = INF
            FLXQT(1) = FLOW
            DO 270 IL=1,NL
               CAP = (WCFC(IL)-WCL(IL))*TKL(IL)*1000.
               IF (CAP.LE.FLOW*DELT) THEN
*                 water flow does not fit into compartment
                  WCL(IL) = WCFC(IL)
                  FLOW = FLOW-CAP/DELT
               ELSE
*                 water flow does fit into compartment
                  WCL(IL) = WCL(IL)+FLOW*DELT/(TKL(IL)*1000.)
                  FLOW = 0.
               END IF
               FLXQT(IL+1) = FLOW
270         CONTINUE
         END IF
         DRAIQT = -FLOW

         WCUMCH = 0.
         DO 280 IL=1,NL
            WCLCH(IL) = (WCL(IL) - WCLQT(IL)) / DELT
            WCUMCH    = WCUMCH+WCLCH(IL)*TKL(IL)*1000.
280      CONTINUE

*        Calculate efficiency of rainfall and irrigation application
         WEFF = ABS (TRWCU/(WCUM+DRAICU+EVSWCU+RNOFCU+INFCU))

         IF (OUTPUT) THEN
*           Output integrals
            CALL OUTDAT (2, 0, 'RAINCU', RAINCU)
            CALL OUTDAT (2, 0, 'IRRCU ', IRRCU)
            CALL OUTDAT (2, 0, 'RNOFCU', RNOFCU)
            CALL OUTDAT (2, 0, 'INFCU ', INFCU)
            CALL OUTDAT (2, 0, 'EVSWCU', EVSWCU)
            CALL OUTDAT (2, 0, 'TRWCU ', TRWCU)
            CALL OUTDAT (2, 0, 'DRAICU', DRAICU)
            CALL OUTDAT (2, 0, 'WCUM  ', WCUM)
            CALL OUTARR ('WCLQT', WCLQT, 1, NL)

*           Output rates
*           CALL OUTDAT (2, 0, 'RAIN  ', RAIN)
            CALL OUTDAT (2, 0, 'IRR   ', IRR)
            CALL OUTDAT (2, 0, 'RNOFF ', RNOFF)
            CALL OUTDAT (2, 0, 'INF   ', INF)
            CALL OUTDAT (2, 0, 'EVSW  ', EVSW)
            CALL OUTDAT (2, 0, 'TRW   ', TRW)
            CALL OUTDAT (2, 0, 'DRAIQT', DRAIQT)
            CALL OUTDAT (2, 0, 'WEFF  ' ,WEFF)

         END IF

*        Checks on calculated rate variables
         DO IL=1,NL
            IF (RESOIL(IL).GT.0.) CALL FATALERR
     &         ('DRSAHE','evaporation rate greater than zero')
         END DO
	
      ELSE IF (ITASK.EQ.3) THEN

*        ===========
*        Integration
*        ===========

*        Checks
         IF (NLXM.LT.NL) CALL FATALERR
     &      ('DRSAHE','too few layers in external arrays')
         WCUMO = WCUM

         DO 400 IL=1,NL
*           Water content per layer
            WCLQT(IL) = INTGRL (WCLQT(IL), WCLCH(IL), DELT)
400      CONTINUE

         DO 420 IL=1,NL+1
            FLXCU(IL) = INTGRL (FLXCU(IL), FLXQT(IL), DELT)
420      CONTINUE

*        Cumulative amounts
         RAINCU = INTGRL (RAINCU, RAIN  , DELT)
         IRRCU  = INTGRL (IRRCU , IRR   , DELT)
         RNOFCU = INTGRL (RNOFCU, RNOFF , DELT)
         INFCU  = INTGRL (INFCU , INF   , DELT)
         EVSWCU = INTGRL (EVSWCU, EVSW  , DELT)
         TRWCU  = INTGRL (TRWCU , TRW   , DELT)
         DRAICU = INTGRL (DRAICU, DRAIQT, DELT)

*        Miscellaneous
         WCUM   = INTGRL (WCUM  , WCUMCH, DELT)
         DSLR   = INTGRL (DSLR  , RDSLR , DELT)

*        Check on value of water content per layer
         DO 430 IL=1,NL
            IF (WCLQT(IL).LT.WCAD(IL)-0.01.OR.
     &          WCLQT(IL).LT.0.) CALL FATALERR
     &          ('DRSAHE','water content less than air dry')
            IF (WCLQT(IL).GT.WCFC(IL)+0.01) CALL FATALERR
     &          ('DRSAHE','water content greater than field capacity')
430      CONTINUE

*        Check on correctness of balance, use relative error
         CHECK = (WCUMO-WCUM)+DELT*(INF+DRAIQT+EVSW+TRW)

         IF (ABS (CHECK/(0.5*(WCUMO+WCUM))).GT.0.001) THEN
            CALL FATALERR ('DRSAHE','error in water balance')
         END IF

      ELSE IF (ITASK.EQ.4) THEN
         CONTINUE
      ELSE
         CALL FATALERR ('DRSAHE','wrong ITASK')
      END IF

      RETURN
      END
