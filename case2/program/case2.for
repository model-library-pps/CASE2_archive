*----------------------------------------------------------------------*
*                                                                      *
*  SUBROUTINE CASE2                                                    *
*  Version 2.2, February 2002                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*  Authors: Wouter Gerritsma (CASE2 version 2.1),                      *
*           Liesje Mommer (version February 1999)                      *
*           Pieter Zuidema (as of April 2001)                          *
*  Date   : February 1997, February 1999 (LM), February 2002 (PZ)      *
*  Purpose: This subroutine is the cacao version of SUCROS, to         *
*           calculate growth in situations of water limited production.*
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
* PLTMOD      Name of plant module used (cocoa or no crop)   -      I  *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITO  I4  Unit used for output file                      -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEI1  C*  Name of first file with plant data             -      I  *
* FILEI3  C*  Name of second file with plant data            -      I  *
* OUTPUT  L4  Flag to indicate if output should be done      -      I  *
* TERMNL  L4  Flag to indicate if simulation is to stop      -     I/O *
* OUTPUTFQ L4 Type of output requested (daily, 10-day, annual) -    I  *
* DOY     R4  Day number within year of simulation (REAL)    d      I  *
* IDOY    I4  Day number within year of simulation (INTEGER) d      I  *
* IYEAR   I4  Year of simulation (INTEGER)                   y      I  *
* DELT    R4  Time step of integration                       d      I  *
* TIME    R4  Time of simulation                             d      I  *
* STTIME  R4  Start time of the simulation                   d      I  *
* LAT     R4  Latitude of site                            dec.degr. I  *
* FRPAR   R4  Fraction PAR in shortwave radiation            -      I  *
* RDD     R4  Daily shortwave radiation                   J/m2/d    I  *
* TMMN    R4  Daily minimum temperature                  degrees C  I  *
* TMMX    R4  Daily maximum temperature                  degrees C  I  *
* NLXM    I4  no. of layers as declared in calling program   -      I  *
* NL      I4  number of layers specified in input file       -      I  *
* TRWL[]  R4  actual transpiration rate per layer           mm/d    I  *
* TKL[]   R4  thickness of soil compartments                 m      I  *
* WCLQT[] R4  volumetric soil water content per layer        -      I  *
* WCWPX[] R4  volumetric water content at wilting point      -      I  *
* WCFCX[] R4  volumetric water content at field capacity     -      I  *
* WCSTX[] R4  volumetric water content at saturation         -      I  *
* EVSC    R4  actual (realized) evaporation rate            mm/d    O  *
* ETRD    R4  Radiation driven part of ETPMD               mm/d     O  *
* ETAE    R4  Dryness driven part of ETPMD                 mm/d     O  *
* PINT    R4  Daily amount of interecepted rain            mm/d     O  *
* GAI     R4  Green area index of cocoa and shade crop      -       O  *
* RAIN    R4  Daily amount of rainfall                     mm/d     I  *
*                                                                      *
* Fatal errors: DELT<1; NPL<700; NPL>2500; LAI(2)>3;SWINPUT not 1 or 2;*
*               CUMTKL<1.5; KDF(2)>0.8 or <0.4; HGHT(2)>40;            *
*               HGHT(2)<HGHL(2); HGHT(1)>10; HGHT(1)<HGHL(1);WTOTI<18.5*
*               or >70;AGEIYR<3 or >39;TMAV>40;TMMN.LT.10;CHKPART1<>1  *
*               PR**ID<0;CHKPART2<>1;CHKDIF>0.01; CHK**>0.01           *
* Warnings   :  AGEYR>40;WTOTPP>70;WRES<0;YRAIN<1000                   *
* Subroutines:  TTUTIL functions                                       *
*               TOTASC, LEAF, ROOT, WUPT, POD                          *
* File usage :  IUNITD                                                 *
*----------------------------------------------------------------------*
      SUBROUTINE CASE2  (PLTMOD, ITASK , IUNITD, IUNITO, IUNITL, FILEI1,
     &				   FILEI3, OUTPUT, TERMNL, OUTPUTFQ,
     &                   DOY   , IDOY  , IYEAR , DELT  , TIME  , STTIME,
     &                   LAT   , FRPAR , RDD   , TMMN  , TMMX  , 
     &                   NLXM  , NL    , TRWL  , TKL   ,
     &                   WCLQT , WCWPX , WCFCX , WCSTX ,
     &                   EVSC  , ETRD  , ETAE  , PINT, 
     &                   GAI   , RAIN )
      USE CHART

      IMPLICIT NONE

*     Formal parameters
      INTEGER   ITASK     , IUNITO  , IUNITD, IUNITL, IDOY, IYEAR
      REAL      DOY       , DELT    , TIME  , STTIME
      REAL      LAT       , FRPAR   , RDD   , TMMN  , TMMX
      INTEGER   OUTPUTFQ  , NLXM      , NL
      REAL      TRWL(NLXM), TKL(NLXM)
      REAL      WCLQT(NLXM),WCWPX(NLXM), WCFCX(NLXM), WCSTX(NLXM)
      REAL      EVSC     , ETRD    , ETAE  , PINT
      REAL      GAI      , RAIN  
      LOGICAL   OUTPUT   , TERMNL
      CHARACTER PLTMOD*(*), FILEI1*(*)
	CHARACTER FILEI3*(*)

*     Function table declarations
      INTEGER ITABLE
      PARAMETER (ITABLE=100)

*     Real function declaration from TTUTIL
      REAL LINT, INTGRL, NOTNUL

*     Standard local declarations
      INTEGER NLA, I1

*     Time and age parameters
      REAL TAU, AGEYR
      REAL AGE, AGEI, AGEIYR,IDOYO

*     Initial stand characteristics 
      INTEGER SWINPUT
      REAL AGBIORA, AGBIORB
      REAL NPL, WWDI, WLVI, WPDI, WRTI, WLRTI, WTRTI 
      REAL WTOTI                  !Note that this is on a per-tree basis
      REAL SLAI

*     Species parameter
      INTEGER INS
      PARAMETER (INS = 2)
      REAL  LAI(INS), SAI(INS), HGHT(INS),  HGHL(INS), KDF(INS), KS(INS)
      REAL  AMAX(INS), EFF(INS), DTGA(INS), FRABS(INS), LAIF(INS)
	REAL TRMIS

*     Soil parameter
      REAL TKLTOT

*     Transpiration parameters
      REAL PTRANS, ATRANS, PCEW, PCEWMN, PENMAN, CROPF
	INTEGER IPC 
	PARAMETER (IPC = 10)
      REAL PCEW10(IPC)
      
*     Biomass of organs 
      REAL WTRT, WLRT, WPD, WLV, WWD, WRT 
      REAL WTWURT, LTWURT 
      REAL WTOT, WTOTCUM, WTOTPP, WTOTMIN

*     Dead biomass
      REAL WLVD, WLRTD, WWDD

*     Assimilation parameters
      REAL TMAVD, TMAV
      REAL AMTMP, AMX, AMINIT 
      REAL GPHOT, TNASS
      REAL MAXLAI
      REAL      AMTMPT(ITABLE)
      INTEGER   IAMTMN        

*     Maintenance 
      REAL MAINLRT, MAINTRT, MAINWD, MAINLV, MAINPD, MAINTS, MAINT
      REAL FHRTWD, HRTWDAGE
      REAL Q10, TREF, TEFF     

*     Assimilate requirements
      REAL ASRQTRT, ASRQLRT, ASRQWD, ASRQLV, ASRQPD, ASRQ   
      
*     CO2 production factors
      REAL CO2LV, CO2PD, CO2LRT, CO2TRT, CO2WD
      
*     Carbon content 
      REAL CFTRT, CFLRT, CFWD, CFLV, CFPD  

*     Reserves parameters
      REAL MINRES, MINCON
      REAL WRES1, WRES2, WRES 
      REAL GRES, DRES, DRES1, DRES2
                
*     Growth rates of organs
      REAL GTOT, GTOT1, GTOT2
      REAL GWD1, GWD 
      REAL GTRT, GRT 
      REAL GLRT, GLRT1, GWURT1
      REAL GLV1, GLV
      REAL GPD1, GPD

*     Turnover rates of organs 
      REAL RTOWURT, TOWURT 
      REAL LRTWURTDR, TOLRT
      REAL RTOLV, TOLV
      REAL WDLVDR, TOWD
      REAL RTOPD, TOPD
      REAL TOREQ

*     Decrease rates of organs
      REAL DLV, DLV1, DLV2, DLRT, DWD
      REAL DLV10(IPC) 

*     Root declarations
      REAL CUMTKL, LTRT
      REAL PI
      PARAMETER(PI= 3.1415927)
      REAL FWURT, ATWURT 
      INTEGER NLBM                             
      PARAMETER (NLBM = 10)
      REAL WSERT(NLBM), AWURT(NLBM), WWURT(NLBM)
                        
*     Partitioning to organs
      REAL FTRTRA, FLRTRA, FLVRA, FWDRA, FPDRA 
      REAL FTRTRB, FLRTRB, FLVRB, FWDRB, FPDRB
      REAL PRLVAC, PRWDAC, PRPDAC, PRTRTAC, PRLRTAC
      REAL PRLVID, PRWDID, PRPDID, PRTRTID, PRLRTID, PRSCALE
      REAL FPD, FWD, FLV, FTRT, FLRT, FSCALE

*     Partitioning checks
      REAL CHKPART1, CHKPART2

*     Carbon blance checks
      REAL CHKIN, CHKDIF, CHKFL, CHKLV, CHKWD, CHKLRT, CHKTRT, CHKPD

*     Leaf parameters 
      REAL AVGLVAGE

*     Pod growth and harvest parameters
      REAL YLDPD, YLDBN, WBNCUM, WPDCUM, BHYLD
      REAL FBEANS, HARPODS, PODVALUE, FATCONTENT
      REAL NCONTBN, PCONTBN, KCONTBN, YNLOSS, YPLOSS, YKLOSS
      REAL FMTDUR, FMTA, FMTB, FMTLOS, MOISTC
      REAL TOPD10(1:10)
      INTEGER IPOD
      REAL      ASRQPDTB(ITABLE),CFPDTB(ITABLE)
      INTEGER   IASRQPDN        ,ICFPDN

*     10-Day totals
      REAL D10GTOT,D10YLDPD,D10YLDBN, D10RAIN, D10RDD, D10HARPD  
      REAL D10COUNT

*     Annual totals
      REAL YRDD, YRAIN, YTRANS, YGPHOT
      REAL YGTOTAB, YGTOT
      REAL YHARPD, YRDEFF, YRNEFF
      REAL YHI, YHINCR, YMNBH, BHSUM
	REAL YMNIPOD,IPODSUM, YMNLAI, LAISUM
      REAL YYLDPD, YYLDBN
      REAL YLVD, YWDD

      SAVE

      IF (DELT.LT.1.0) CALL FATALERR
     &   ('CASE2','FE1 - Time step (DELT) too small: <1 day')

      IF (ITASK.EQ.1) THEN       
*        ----------------------
*        Initialization section
*        ----------------------

*        Send title(s) to output file
         CALL OUTCOM (PLTMOD)
         CALL OUTCOM ('CASE2, CAcao Simulation Engine     Version 2.2')
         CALL OUTCOM ('       February 2002')

*        Read plant parameters from plant.dat with UNchangeble values
         CALL RDINIT (IUNITD  , IUNITL, FILEI1)

*        Regression coefficients to relate initial biomass and age
         CALL RDSREA ('AGBIORA' , AGBIORA)       ![kg DW d-1
         CALL RDSREA ('AGBIORB' , AGBIORB)       ![kg DW]

*        Photosynthesis  parameters
         CALL RDSREA ('AMX'   , AMX   )          ![kg CO2 ha-1 leaf h-1] 
         CALL RDSREA ('EFF'   , EFF   )    ![(CO2 ha-1 h-1)/(J m-2 s-1)]
         CALL RDAREA ('AMTMPT', AMTMPT,  ITABLE, IAMTMN)  ![-]
         CALL RDSREA ('AMINIT', AMINIT)          ![-] 
         CALL RDSREA ('MAXLAI', MAXLAI)          ![ha leaf ha-1 ground]
         CALL RDSREA ('KDFL'  , KDF(1))          ![-] 
         CALL RDSREA ('KDFT'  , KS(1))           ![-]

*        Tissue maintenance coefficients
         CALL RDSREA ('MAINLRT', MAINLRT)   ![kg CH2O kg-1 DW lat.roots]
         CALL RDSREA ('MAINTRT', MAINTRT)     ![kg CH2O kg-1 DW taproot]
         CALL RDSREA ('MAINWD', MAINWD)          ![kg CH2O kg-1 DW wood]
         CALL RDSREA ('MAINLV', MAINLV)        ![kg CH2O kg-1 DW leaves]
         CALL RDSREA ('MAINPD', MAINPD)          ![kg CH2O kg-1 DW pods]
         CALL RDSREA ('Q10'   , Q10   )          ![-]
         CALL RDSREA ('TREF'  , TREF  )          ![degr C]

*        Reserves
         CALL RDSREA ('MINCON', MINCON)          ![-]

*        Tissue asssimilate requirements
         CALL RDSREA ('ASRQLRT', ASRQLRT)
	                                  ![kg CH2O kg-1 DW lateral roots]
         CALL RDSREA ('ASRQTRT', ASRQTRT)     ![kg CH2O kg-1 DW taproot]
         CALL RDSREA ('ASRQWD', ASRQWD)          ![kg CH2O kg-1 DW wood]
         CALL RDSREA ('ASRQLV', ASRQLV)        ![kg CH2O kg-1 DW leaves]
         CALL RDAREA ('ASRQPDTB', ASRQPDTB , ITABLE, IASRQPDN ) 
	                  ![kg CH2O kg-1 DW pods]

*        Tissue carbon content
         CALL RDSREA ('CFLRT' , CFLRT )    ![kg C kg-1 DW lateral roots]
         CALL RDSREA ('CFTRT' , CFTRT )          ![kg C kg-1 DW taproot]
         CALL RDSREA ('CFWD'  , CFWD  )          ![kg C kg-1 DW wood]
         CALL RDSREA ('CFLV'  , CFLV  )          ![kg C kg-1 DW leaves]
         CALL RDAREA ('CFPDTB', CFPDTB , ITABLE, ICFPDN )
	                  ![kg C kg-1 DW pods]
   
*        Time control
         CALL RDSREA ('TAU',  TAU)               ![d]

*        Pods and beans parameters
         CALL RDSREA ('NCONTBN',NCONTBN)         ![kg N kg-1 bean DW]
	   CALL RDSREA ('PCONTBN',PCONTBN)         ![kg P kg-1 bean DW]
         CALL RDSREA ('KCONTBN',KCONTBN)         ![kg K kg-1 bean DW]

*        Fermentation parameters
         CALL RDSREA ('FMTA'  , FMTA)            ![d-1]
         CALL RDSREA ('FMTB'  , FMTB)            ![-]
         
*        Partitioning parameters
         CALL RDSREA ('RTOWURT' , RTOWURT)       ![d-1]
         CALL RDSREA ('WDLVDR' , WDLVDR)
	                           ![kg DW dead wood kg-1 DW dead leaves]
         CALL RDSREA ('LRTWURTDR', LRTWURTDR)    ![kg DW dead non-water
	   ! uptaking lateral roots kg-1 DW dead water-uptaking roots]
         CALL RDSREA ('FTRTRA' , FTRTRA)
	                              ![kg DW taproot kg-1 DW whole plant]
         CALL RDSREA ('FLRTRA' , FLRTRA)
	                        ![kg DW lateral roots kg-1 DW whole plant]
         CALL RDSREA ('FLVRA' , FLVRA) ![kg DW leaf kg-1 DW whole plant]
         CALL RDSREA ('FWDRA' , FWDRA) ![kg DW wood kg-1 DW whole plant]
         CALL RDSREA ('FPDRA' , FPDRA)  ![kg DW pod kg-1 DW whole plant]
         CALL RDSREA ('FTRTRB' , FTRTRB)         ![kg DW taproot]
         CALL RDSREA ('FLRTRB' , FLRTRB)         ![kg DW lateral root]    
         CALL RDSREA ('FLVRB' , FLVRB)           ![kg DW leaves kg-1]
         CALL RDSREA ('FWDRB' , FWDRB)           ![kg DW wood]
         CALL RDSREA ('FPDRB' , FPDRB)           ![kg DW pod]
         CALL RDSREA ('AVGLVAGE', AVGLVAGE)      ![d]
         CALL RDSREA ('FWURT', FWURT)            ![-]
         CALL RDSREA ('HRTWDAGE', HRTWDAGE)      ![d]
         CALL RDSREA ('WTOTMIN', WTOTMIN)        ![kg DW plant-1]

*        Close plant file
         CLOSE (IUNITD)

*        Check total soil depth: should be > 1.5 m
         TKLTOT = 0.                             ![m]
         DO I1 = 1,NLXM
	      TKLTOT = TKLTOT + TKL(I1)            ![m]
         ENDDO
         IF (TKLTOT.LT.1.5) CALL FATALERR ('CASE2',
     &     'FE2 - Total thickness of soil layers too low: <1.5 m')

*        Read plant parameters from SECOND data file with changeble values
         CALL RDINIT (IUNITD  , IUNITL, FILEI3)

*        Planting density
         CALL RDSREA ('NPL'   , NPL    )         ![ha-1]

*        Tree height
         CALL RDSREA ('HGHL'  , HGHL(1))         ![m]
         CALL RDSREA ('HGHT'  , HGHT(1))         ![m]

*        Shade tree characteristics
         CALL RDSREA ('SLAI'  , LAI(2))          ![ha leaf ha-1 ground]
         CALL RDSREA ('SKDFL' , KDF(2))          ![-]
         CALL RDSREA ('SHGHL' , HGHL(2) )        ![m]
         CALL RDSREA ('SHGHT' , HGHT(2) )        ![m]

*        Pod-related parameters
         CALL RDSREA ('FMTDUR', FMTDUR)          ![d]
         CALL RDSREA ('MOISTC', MOISTC)          ![-]
         CALL RDSREA ('FBEANS', FBEANS)          ![-]
         CALL RDSREA ('PODVALUE', PODVALUE)      ![-]
         CALL RDSREA ('FATCONTENT',FATCONTENT)   ![-]
      
*        Check planting density: should be between 700 and 2000 ha-1
         IF (NPL.LT.700.) CALL FATALERR ('CASE2',
     &     'FE3 - Planting density too low: < 700 ha-1')
         IF (NPL.GT.2500) CALL FATALERR ('CASE2',
     &  'FE4 - Planting density too high: > 2500 ha-1')

*        Check extiction value for shade trees: should be between 0.4-0.8
         IF ((KDF(2).GT.0.8).OR.(KDF(2).LT.0.4)) CALL FATALERR ('CASE2',
     &  'FE5 - Extinction coefficient of shade trees too 
     &   high or too low')

*        Check value of LAI for shade trees: should not exceed 3 ha ha-1
         IF (LAI(2).GT.3.) CALL FATALERR ('CASE2',
     &   'FE6 - LAI of shade trees too high: >3. ha ha-1')

*        Check value of shade tree:height: should not exceed 40 m
*        and upper height should be larger than lower height
         IF ((HGHT(2).GT.40.).OR.(HGHL(2).GT.40.)) CALL FATALERR 
     &   ('CASE2','FE7 - Shade tree height too high: should be <40 m')
         IF (HGHT(2).LE.HGHL(2)) CALL FATALERR ('CASE2',
     &  'FE8 - Upper shade tree height <= lower height')

*        Check value of shade tree:height: should not exceed 20 m
*        and upper height should be larger than lower height
         IF ((HGHT(1).GT.20.).OR.(HGHL(1).GT.20.)) CALL FATALERR 
     &   ('CASE2','FE9 - Cacao tree height too high: should be <20 m')
         IF (HGHT(1).LE.HGHL(1)) CALL FATALERR ('CASE2',
     &  'FE10 - Upper cacao tree height <= lower height')

*        Initial weight or size of cocoa trees 
*        First determine whether input is age or weight of the plants.
         CALL RDSINT ('SWINPUT' , SWINPUT)       ![-]
         IF (SWINPUT.EQ.2) THEN                                        
*          Input is weight
           CALL RDSREA ('WTOTI'  , WTOTI)        ![kg DW plant-1]
*          Check value of initial size: should be >= 18.5 and <70 kg
           IF (WTOTI.LT.8.5) CALL FATALERR ('CASE2',
     &       'FE11 - Intitial tree biomass too low: < 18.5 kg')
           IF (WTOTI.GE.70.) CALL FATALERR ('CASE2',
     &       'FE12 - Intitial tree biomass too high: > 70.0 kg')
           AGEI = EXP((WTOTI-AGBIORB)/AGBIORA)   ![d]
         ELSEIF (SWINPUT.EQ.1) THEN
*          Input is age
           CALL RDSREA ('AGEIYR'  , AGEIYR   )   ![d]
*          Check value of initial age: should be >= 3 and < 40. yr
           IF (AGEIYR.LT.3.) CALL FATALERR ('CASE2',
     &       'FE13 - Initial tree age too low: < 3 y')
           IF (AGEIYR.GT.39.) CALL FATALERR ('CASE2',
     &       'FE14 - Initial tree age too high: > 39 y')
*          calculate age in days and round to get entire days
           AGEI  = REAL(INT(AGEIYR * 365.))      ![d]
           WTOTI = AGBIORA * LOG(AGEI) + AGBIORB ![kg DW plant-1]
         ELSE
           CALL FATALERR ('CASE2','FE15 - Value of input switch SWINPUT 
     &     is wrong, should be 1 or 2')
         ENDIF 

*        Close plant file
         CLOSE (IUNITD)

*        Initial weight of cocoa tree components (per ha basis)
         PRSCALE =0.
         IF (WTOTI.LT.WTOTMIN) THEN 
           PRSCALE= (FLVRA * WTOTI + FLVRB)+ (FWDRA * WTOTI + FWDRB)+
     &              (FTRTRA *WTOTI + FTRTRB)+(FLRTRA* WTOTI +FLRTRB)+
     &              (FPDRA * WTOTI + FPDRB)     
           WLVI   = (FLVRA*WTOTI+FLVRB)*PRSCALE/WTOTI*NPL  ![kg DW ha-1]
           WWDI   = (FWDRA*WTOTI+FWDRB)*PRSCALE/WTOTI*NPL  ![kg DW ha-1]
           WTRTI  = (FTRTRA*WTOTI+FTRTRB)*PRSCALE/WTOTI*NPL ![kg DW ha-1]
           WLRTI  = (FLRTRA*WTOTI+FLRTRB)*PRSCALE/WTOTI*NPL ![kg DW ha-1]
           WPDI   = 0.     ![kg DW ha-1]
         ELSE
           WLVI   = (FLVRA * WTOTI + FLVRB) * NPL   ![kg DW ha-1]
           WWDI   = (FWDRA * WTOTI + FWDRB) * NPL   ![kg DW ha-1]
           WTRTI  = (FTRTRA * WTOTI + FTRTRB) * NPL ![kg DW ha-1]
           WLRTI  = (FLRTRA * WTOTI + FLRTRB) * NPL ![kg DW ha-1]
           WPDI   = (FPDRA * WTOTI + FPDRB) *NPL    ![kg DW ha-1]
         ENDIF

*        Initialise state variables
         WLV   = WLVI   ![kg DW ha-1]
         WWD   = WWDI   ![kg DW ha-1]
         WTRT  = WTRTI  ![kg DW ha-1]
         WLRT  = WLRTI  ![kg DW ha-1]
         WRT   = WTRT + WLRT                     ![kg DW ha-1]
         WPD   = WPDI   ![kg DW ha-1]
         WTOT   = WLV + WWD + WPD + WRT          ![kg DW ha-1]
         WTOTPP = WTOT / NPL                     ![kg DW plant-1]

*        Plant age in days and years
         AGE = AGEI     ![d]
         AGEYR = AGE / 365.                      ![yr]

*        Initialise subroutines
         CALL ROOT (ITASK, IUNITD,IUNITL,FILEI1,DELT,NL,NLA, NLXM,NLBM, 
     &          I1,TKL, CUMTKL, WCLQT, WCWPX, AGE, NPL, FWURT,
     &          WLRT, WTRT, LTRT, WWURT, WTWURT, LTWURT, 
     &          ATWURT, AWURT, WSERT)
      
         CALL WUPT (ITASK,IUNITD,IUNITL,FILEI1,NLXM,NL, NLA, NLBM,  
     &          TKL   , WCLQT , WCWPX , WCFCX , WCSTX , ATWURT ,
     &          AWURT, ETRD  , ETAE  ,  EVSC  , GAI   ,
     &          TRWL  , PINT, PTRANS, ATRANS, PCEW  ,
     &          PENMAN, CROPF  )

*        Initialise array of PCEW values of the last 10 days, and fill
*        with ten times the same value
         DO I1 = 1,10
           PCEW10(I1) = 1.                       ![-]
         ENDDO
         PCEWMN = SUM(PCEW10(1:10)) / 10.        ![-]

*        Calculate fraction light tranmission through shade canopy,  
*        assuming that cocoa and shade crowns do not overlap.
         TRMIS = EXP(-KDF(2) * LAI(2))

         CALL LEAF (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,AGE,NPL,
     &             TRMIS,GLV,PCEWMN,WTOTPP,
     &             WLV, DLV, DLV1, DLV2, LAI(1))

*        Determine average temperature (TMAV)
         TMAV  = 0.5 * (TMMX + TMMN)             ![degr C]

         CALL POD (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,TMAV,
     &            GPD, WPD, YLDPD, IPOD, BHYLD)

*        Determing initial value for leaf turnover (TOLV), using given
*        maximum leaf age (MAXLAG) and leaf weight (WLV).
         DO I1 = 1,10                                                        
           DLV10(I1) = WLV/AVGLVAGE              ![kg DW ha-1 d-1]
         ENDDO
         TOLV = SUM(DLV10(1:10)) / 10.           ![kg DW ha-1 d-1]

*        Determing initial value for pod turnover (TOPD), using pod 
*        weight per pod category.
         DO I1 = 1,10
          TOPD10(I1) = WPD/REAL(IPOD)            ![kg DW ha-1 d-1]
         ENDDO
         TOPD = SUM(TOPD10(1:10)) / 10.          ![kg DW ha-1 d-1]

*        Set intergrated state variable to zero, for reruns
*        Rates
         GLV     = 0.   ![kg DW ha-1 d-1]
         DLV     = 0.   ![kg DW ha-1 d-1]
         GPD     = 0.   ![kg DW ha-1 d-1]
 	   YLDPD   = 0.   ![kg DW ha-1 d-1]
	   YLDBN   = 0.   ![kg DW ha-1 d-1]
         HARPODS = 0.   ![ha-1 d-1]
         GTRT    = 0.   ![kg DW ha-1 d-1]
         GLRT    = 0.   ![kg DW ha-1 d-1]
         DLRT    = 0.   ![kg DW ha-1 d-1]
         GWD     = 0.   ![kg DW ha-1 d-1]
         DWD     = 0.   ![kg DW ha-1 d-1]

*        Dead biomass 
         WLVD  = 0.     ![kg DW ha-1]
         WLRTD = 0.     ![kg DW ha-1]
         WWDD  = 0.

*        Annual totals        
         YRDD   = 0.    ![MJ m-2 on yearly basis]
         YRAIN  = 0.    ![mm water on yearly basis]
         YTRANS = 0.    ![mm water on yearly basis]
         YYLDPD = 0.    ![kg DW ha-1 y-1]
         YYLDBN = 0.    ![kg dry fermented cocoa ha-1 y-1]
         YGPHOT = 0.    ![kg CH2O ha-1 y-1]
         YGTOT  = 0.    ![kg DW ha-1 y-1]
         YLVD   = 0.    ![kg DW ha-1 y-1]
         YWDD   = 0.    ![kg DW ha-1 y-1]
         YHARPD = 0.    ![y-1]
         YRDEFF = 0.    ![kg DW ha -1 (MJ m-2)-1]
         YRNEFF = 0.    ![kg DW ha -1 mm]
         YHI    = 0.    ![kg DW fermented beans kg-1 DW above ground]
         YHINCR = 0.    ![kg DW fermented beans kg-1 DW above ground]
         BHSUM  = 0.    ![-] 
         YMNBH  = 0.    ![-]
         YNLOSS = 0.    ![kg N ha-1 yr-1]
         YPLOSS = 0.    ![kg P ha-1 yr-1]
         YKLOSS = 0.    ![kg K ha-1 yr-1]
         LAISUM = 0.
         YMNLAI = 0.    ![ha leaf ha-1 ground ]
         IPODSUM= 0.
	   YMNIPOD= 0.    ![d]

*        10-day totals         
         D10GTOT  = 0.  ![kg DW ha-1 10d-1]
         D10YLDPD = 0.  ![kg DW ha-1 10d-1]
         D10YLDBN = 0.  ![kg DW ha-1 10d-1]
         D10RAIN  = 0.  ![mm 10d-1]
         D10RDD   = 0.  ![mm 10d-1]
         D10HARPD = 0.  ![10d-1]
         D10COUNT = 0.

*        Cumulatives
         WTOTCUM = 0.   ![kg DW ha-1]
         WPDCUM = 0.    ![kg DW ha-1]
         WBNCUM = 0.    ![kg DW pods ha-1 yr-1]

*        Variables for carbon balance check
         TNASS  = 0.    ![kg CO2 ha-1]
         CHKDIF = 0.    ![-]
         CHKIN  = 0.    ![kg C ha-1]
         CHKFL  = 0.    ![kg C ha-1]
         CHKLV  = 0.    ![kg C ha-1]
         CHKWD  = 0.    ![kg C ha-1]
         CHKLRT = 0.    ![kg C ha-1]
         CHKTRT = 0.    ![kg C ha-1]
         CHKPD  = 0.    ![kg C ha-1]
        
*        Green Area Index is total leaf area of cocoa and shade crop 
         GAI = 0.       ![m2 green leaves m2 soil]
         DO I1=1,INS
            GAI    = GAI + LAI(I1)             ![m2 green leaf m-2 soil]
         ENDDO

*        Reserve weights
         WRES   = MINCON * WTOT                  ![kg CH2O ha-1]
         MINRES = WRES  ![kg CH2O ha-1]

*        Determine assimilate requirement and C content of pods as 
*        a function of nib fat content
         ASRQPD = LINT(ASRQPDTB, IASRQPDN, FATCONTENT)![kg CH2O kg-1 DW]
         CFPD   = LINT(CFPDTB , ICFPDN  , FATCONTENT)  ![kg C kg-1 DW]

*        CO2 production factors
         CO2LRT = 44./12. * (ASRQLRT*12./30. - CFLRT)  ![kg CO2 kg-1 DW]
         CO2TRT = 44./12. * (ASRQTRT*12./30. - CFTRT)  ![kg CO2 kg-1 DW]
         CO2WD  = 44./12. * (ASRQWD*12./30.  - CFWD)   ![kg CO2 kg-1 DW]
         CO2LV  = 44./12. * (ASRQLV*12./30.  - CFLV)   ![kg CO2 kg-1 DW]
         CO2PD  = 44./12. * (ASRQPD*12./30.  - CFPD)   ![kg CO2 kg-1 DW]

         IDOYO  = IDOY  ![d]


      ELSE IF (ITASK.EQ.2) THEN
*     ------------------------
*     Rate calculation section
*     ------------------------
                    
         AGE = AGEI + TIME - STTIME              ![d]
         AGEYR = AGE / 365.                      ![yr]

*        Weather data
*        Average temperature (TMAV) and day time average (TMAVD)
         TMAV  = 0.5 * (TMMX + TMMN)             ![degr C]
         TMAVD = TMMX - 0.25 * (TMMX-TMMN)       ![degr C]
 
         IF (TMAV.GT.40.) CALL FATALERR
     &     ('CASE2','FE16 - Av. Temp greater than 40 C.')
         IF (TMMN.LT.10.) CALL FATALERR
     &     ('CASE2','FE17 - Minimum temperature below 10 C.')
        
*        Transpiration and water uptake
         CALL ROOT (ITASK,IUNITD,IUNITL,FILEI1,DELT,NL,NLA, NLXM,NLBM,
     &             I1,TKL, CUMTKL, WCLQT, WCWPX, AGE, NPL, FWURT,
     &             WLRT, WTRT, LTRT, 
     &             WWURT, WTWURT, LTWURT, ATWURT, AWURT, WSERT)

         CALL WUPT (ITASK,IUNITD,IUNITL,FILEI1,NLXM,NL, NLA, NLBM,  
     &                  TKL   , WCLQT , WCWPX , WCFCX , WCSTX , ATWURT ,
     &                  AWURT, ETRD  , ETAE  ,  EVSC  , GAI  ,
     &                  TRWL  ,  PINT, PTRANS, ATRANS, PCEW  , 
     &                  PENMAN, CROPF )

*        Calculate mean PCEW value over the last 10 days 
         DO I1 = 9,1,-1
           PCEW10(I1+1) = PCEW10(I1)             ![-]
         ENDDO
         PCEW10(1) = PCEW                        ![-]
         PCEWMN    = SUM(PCEW10(1:10)) / 10.     ![-]

*        Carbohydrate production and respiration
*        Leaf CO2 assimilation
*        Interpolate the temperature correction for the maximum rate
*        of photosynthesis (AMAX(1)) from table (AMPTP)
         AMTMP   = LINT(AMTMPT, IAMTMN, TMAVD)    ![-]
         AMAX(1) = AMX * AMTMP * AMINIT           ![kgCO2 ha-1 h-1]

*        LAI of the cocoa trees is set to the maximum of MAXLAI. At high
*        values of LAI the Gaussian integrations in TOTASC does 
*        not work properly yielding low estimates for photosynthesis
*        (as light availability in large part of the canopy is low) 
         LAIF(1) = MIN(LAI(1), MAXLAI)           ![ha leaf ha-1 ground]
         LAIF(2) = LAI(2)                        ![ha leaf ha-1 ground]

         CALL TOTASC (IDOY, INS, LAT, RDD, FRPAR, KDF, KS,
     &                AMAX, EFF, LAIF, SAI, HGHT, HGHL,
     &                FRABS, DTGA)

*        Carbohydrate production
         GPHOT = DTGA(1) * 30./44. * PCEW        ![kgCH2O ha-1 d-1]

*        Fraction of wood that is heartwood
         IF (AGE.GT.HRTWDAGE) THEN
           FHRTWD = (AGE - HRTWDAGE)/AGE         ![-]
         ELSE 
           FHRTWD = 0.                           ![-]
         END IF

*        Maintenance respiration
         MAINTS = MAINLRT*WLRT + MAINTRT*WTRT*(1-FHRTWD) +
     &            MAINWD*WWD*(1-FHRTWD) + MAINLV*WLV + MAINPD*WPD  
                                                 ![kgCH2O ha-1 d-1]
         TEFF   = Q10**((TMAV - TREF)/10.)       ![-]
         MAINT  = MAINTS * TEFF                  ![kgCH2O ha-1 d-1]
         GRES   = GPHOT - MAINT                  ![kgCH2O ha-1 d-1]
       
*        Growth consists of two parts:
*        1. Replacement of organs with a certain turnover to maintain 
*           biomass levels. This is done for leaves, lateral roots, 
*           wood and fruits. Turnover rates are used. 
*        2. Actual growth of all organs. This is based on fraction 
*           partitioning depending on plant biomass.

*        1.  Replacement
*        1a. Of water-uptaking roots (< 2 mm diameter)
*        Death rate depends on relative death rate and weight of roots.
         TOWURT = RTOWURT * WTWURT               ![kg DW ha-1 d-1]

*        1b. Of coarse lateral roots
*        Calculated as a percentage of loss rate of water-uptaking roots
         TOLRT  = LRTWURTDR * TOWURT             ![kg DW ha-1 d-1]

*        1c. Of leaves
*        First calculate mean leaf loss over the last 10 days  
         DO I1 = 9,1,-1
           DLV10(I1+1) = DLV10(I1)               ![kg DW ha-1 d-1]
         ENDDO
         DLV10(1)    = DLV                       ![kg DW ha-1 d-1]
*        then calculate leaf turnover rate as the mean of leaf loss 
         TOLV = SUM(DLV10(1:10)) / 10.           ![kg DW ha-1 d-1]
        
*        1d. Of wood
*        This is calculated as a percentage of the leaf loss rate
         TOWD = WDLVDR * TOLV                    ![kg DW ha-1 d-1]

*        1e. Of pods
*        Turnover rate of pods is in fact equal to the yield 
*        First calculate mean yield over the last 10 days
         DO I1 = 9,1,-1
           TOPD10(I1+1) = TOPD10(I1)             ![kg DW ha-1 d-1]
         ENDDO
         TOPD10(1)    = YLDPD                    ![kg DW ha-1 d-1] 
         TOPD   = SUM(TOPD10(1:10)) / 10.        ![kg DW ha-1 d-1]

*        Use reserve pool to replace dead biomass of these organs
*        First check whether available reserve mass is more than 
*        the minimum reserve size
         IF ((WRES + GRES * DELT) .GT. (MINRES)) THEN 
           WRES1    = WRES + GRES * DELT - MINRES![kgCH2O ha-1]
         ELSE
           WRES1    = 0                          ![kgCH2O ha-1]
         ENDIF

*        Then check whether available reserve mass is sufficiently large
*        to replace the lost part of all organs.
         IF (WRES1/DELT .GT. ((TOWURT+TOLRT)*ASRQLRT  
     &     + TOLV*ASRQLV*PCEW + TOWD*ASRQWD + TOPD*ASRQPD)) THEN
           GWURT1  = TOWURT                      ![kg DW ha-1 d-1]
           GLRT1   = TOLRT                       ![kg DW ha-1 d-1]
           GLV1    = TOLV * PCEW                 ![kg DW ha-1 d-1]
           GWD1    = TOWD                        ![kg DW ha-1 d-1]
           GPD1    = TOPD                        ![kg DW ha-1 d-1]
           GTOT1    = GWURT1 + GLV1 + GWD1 + GPD1![kg DW ha-1 d-1]
*          If not, reserve weight is distributed proportional to turnover
         ElSEIF ((WRES1/DELT .GT. 0.).AND.(WRES1/DELT .LE. 
     &     ((TOWURT+TOLRT)*ASRQLRT + TOLV*ASRQLV*PCEW + 
     &     TOWD*ASRQWD + TOPD*ASRQPD))) THEN
*          total required CH2O to account for turnover of organs = TOREQ
           TOREQ  = (TOWURT+TOLRT)*ASRQLRT + TOLV*ASRQLV*PCEW 
     &               + TOWD*ASRQWD + TOPD*ASRQPD     
           GWURT1  = ((TOWURT*ASRQLRT)/TOREQ * WRES1/DELT)/ ASRQLRT 
           GLRT1   = ((TOLRT*ASRQLRT) /TOREQ * WRES1/DELT)/ ASRQLRT 
           GWD1    = ((TOWD*ASRQWD)   /TOREQ * WRES1/DELT)/ ASRQWD  
           GLV1    = ((TOLV*ASRQLV*PCEW)/TOREQ * WRES1/DELT)/ASRQLV 
           GPD1    = ((TOPD*ASRQPD)   /TOREQ * WRES1/DELT)/ ASRQPD  
           GTOT1    = GWURT1 + GLV1 + GWD1 + GPD1
		                                       !all [kg DW ha-1 d-1]
         ELSE 
           GWURT1  = 0.                          ![kg DW ha-1 d-1]          
           GLRT1   = 0.                          ![kg DW ha-1 d-1]          
           GLV1    = 0.                          ![kg DW ha-1 d-1]          
           GWD1    = 0.                          ![kg DW ha-1 d-1]          
           GPD1    = 0.                          ![kg DW ha-1 d-1]          
           GTOT1   = 0.                          ![kg DW ha-1 d-1]
         ENDIF

         DRES1 = (GWURT1+GLRT1)*ASRQLRT+GLV1*ASRQLV+GWD1*ASRQWD
     &           +GPD1*ASRQPD                    ![kgCH2O ha-1 d-1]

*        2. Actual growth of all organs. 
*        First determine actual proportions of total DW in  organs
         PRLVAC = WLV / WTOT                     ![-]
         PRWDAC = WWD / WTOT                     ![-]
         PRPDAC = WPD / WTOT                     ![-]
         PRTRTAC = WTRT / WTOT                   ![-]
         PRLRTAC = WLRT / WTOT                   ![-]
         CHKPART1   = PRLVAC + PRWDAC + PRPDAC + PRTRTAC + PRLRTAC ![-]
 
*        Terminate program if CHKPART1 is NOT 1
         IF ((CHKPART1 .GT. 1.005) .OR. (CHKPART1 .LT. 0.995)) 
     &	 CALL FATALERR ('CASE2','FE18 - Actual allometry is wrong: 
     &     sum of proportions (CHKPART1) does not equal 1')
      
*        Then determine "ideal" proportions of total DW in different 
*        organs, based on allometric relations. For leaves and fine 
*        roots, the ideal proportion in modified by the water 
*        availability over the last 10 days. Note that the ideal 
*        proportions are calculated using the total dry weight per plant
*        (WTOTPP), as the regression equations are on a per-plant basis.
*        In two cases: 
*        1. For small plants not bearing pods. The ideal proportion of 
*           biomass in pods PRPDID = 0.
         IF ((WTOTPP).LT.WTOTMIN) THEN
           PRSCALE= (FLVRA*WTOTPP+FLVRB)/WTOTPP*PCEW +
     &              (FWDRA*WTOTPP+FWDRB)/WTOTPP    +
     &              (FTRTRA*WTOTPP+FTRTRB)/WTOTPP  +
     &              (FLRTRA*WTOTPP+FLRTRB)/WTOTPP*(2-PCEW)           ![-]
           PRLVID =((FLVRA*WTOTPP+FLVRB)/WTOTPP*PCEW)/PRSCALE        ![-]
           PRWDID =((FWDRA*WTOTPP+FWDRB)/WTOTPP)/PRSCALE             ![-]
           PRPDID =0.                                                ![-]
           PRTRTID=((FTRTRA*WTOTPP+FTRTRB)/WTOTPP)/PRSCALE           ![-]
           PRLRTID=((FLRTRA*WTOTPP+FLRTRB)/WTOTPP*(2-PCEW))/PRSCALE  ![-]
*        2. For larger plants bearing pods. 
         ELSE
           PRSCALE=(FLVRA*WTOTPP+FLVRB)/WTOTPP*PCEW +
     &             (FWDRA*WTOTPP+FWDRB)/WTOTPP    +
     &             (FPDRA*WTOTPP+FPDRB)/WTOTPP    +
     &             (FTRTRA*WTOTPP+FTRTRB)/WTOTPP  +
     &             (FLRTRA*WTOTPP+FLRTRB)/WTOTPP*(2-PCEW)           ![-]
           PRLVID =((FLVRA*WTOTPP+FLVRB)/WTOTPP * PCEW) / PRSCALE   ![-]
           PRWDID =((FWDRA*WTOTPP+FWDRB)/WTOTPP) / PRSCALE          ![-]
           PRPDID =((FPDRA*WTOTPP+FPDRB)/WTOTPP) / PRSCALE          ![-]
           PRTRTID=((FTRTRA*WTOTPP+FTRTRB)/WTOTPP) / PRSCALE        ![-]
           PRLRTID=((FLRTRA*WTOTPP+FLRTRB)/WTOTPP*(2-PCEW))/PRSCALE ![-]
         ENDIF
 
*       Terminate program if one of ideal proportions is below 0 (this 
*       is possible when total biomass is low.)
        IF ((PRLVID.LT.0.).OR.(PRWDID.LT.0.).OR.(PRTRTID.LT.0.).OR.
     &     (PRLRTID.LT.0.).OR.(PRPDID.LT.0.)) CALL FATALERR
     &     ('CASE2','FE19 - The ideal biomassfraction 
     &     of one or more organs is below 0')

*       Correct partitioning fractions for deviations from the "ideal" 
*       proportions. In two cases: 
*       1. For small plants not bearing pods.
        IF ((WTOTPP).LT.WTOTMIN) THEN
          FSCALE= MAX(0.,FLVRA*(PRLVID-PRLVAC)/PRLVAC) + 
     &            MAX(0.,FWDRA*(PRWDID-PRWDAC)/PRWDAC) + 
     &            MAX(0.,FTRTRA*(PRTRTID-PRTRTAC)/PRTRTAC) +
     &            MAX(0.,FLRTRA *(PRLRTID-PRLRTAC)/PRLRTAC)         ![-]
          FLV   = MAX(0.,(FLVRA *(PRLVID-PRLVAC)/PRLVAC)/FSCALE)    ![-]
          FWD   = MAX(0.,(FWDRA *(PRWDID-PRWDAC)/PRWDAC)/FSCALE)    ![-]
          FPD   = 0.                                                ![-]
          FTRT  = MAX(0.,(FTRTRA*(PRTRTID-PRTRTAC)/PRTRTAC)/FSCALE) ![-]
          FLRT  = MAX(0.,(FLRTRA*(PRLRTID-PRLRTAC)/PRLRTAC)/FSCALE) ![-]
          CHKPART2   = FLV + FWD + FPD + FTRT + FLRT                ![-]
*       2. For larger plants bearing pods.
        ELSE
*         The below if statements prevents PRPDAC from being zero and 
*         PRPDID/PRPDAC from being infinitely large.
          IF (PRPDAC.EQ.0.) PRPDAC = 0.001
            FSCALE= MAX(0.,FLVRA*(PRLVID-PRLVAC)/PRLVAC) + 
     &              MAX(0.,FWDRA*(PRWDID-PRWDAC)/PRWDAC) + 
     &              MAX(0.,FTRTRA*(PRTRTID-PRTRTAC)/PRTRTAC) +
     &              MAX(0.,FLRTRA *(PRLRTID-PRLRTAC)/PRLRTAC) +
     &              MAX(0.,FPDRA *(PRPDID-PRPDAC)/PRPDAC)           ![-]
          IF (FSCALE .EQ. 0.) THEN
            FLV = FLVRA / (FLVRA + FWDRA +FPDRA +FTRTRA +FLRTRA)    ![-]
            FWD = FWDRA / (FLVRA + FWDRA +FPDRA +FTRTRA +FLRTRA)    ![-]
            FPD = FPDRA / (FLVRA + FWDRA +FPDRA +FTRTRA +FLRTRA)    ![-]
            FTRT= FTRTRA / (FLVRA + FWDRA +FPDRA +FTRTRA +FLRTRA)   ![-]
            FLRT= FLRTRA / (FLVRA + FWDRA +FPDRA +FTRTRA +FLRTRA)   ![-]
            CHKPART2 = FLV + FWD + FPD + FTRT + FLRT                ![-]
          ELSE
            FLV    = MAX(0.,(FLVRA  * (PRLVID-PRLVAC)/PRLVAC)/ FSCALE) 
            FWD    = MAX(0.,(FWDRA  * (PRWDID-PRWDAC)/PRWDAC)/ FSCALE) 
            FPD    = MAX(0.,(FPDRA  * (PRPDID-PRPDAC)/PRPDAC)/ FSCALE) 
            FTRT   = MAX(0.,(FTRTRA * (PRTRTID-PRTRTAC)/PRTRTAC)/FSCALE)
            FLRT   = MAX(0.,(FLRTRA * (PRLRTID-PRLRTAC)/PRLRTAC)/FSCALE)
            CHKPART2   = FLV + FWD + FPD + FTRT + FLRT              
		                                                      !all [-]
          ENDIF
        ENDIF

*       Terminate program if CHKPART2 is NOT 1
        IF ((CHKPART2 .GT. 1.005) .OR. (CHKPART2 .LT. 0.995)) 
     &	CALL FATALERR ('CASE2','FE20 - Partitioning is wrong,
     &    total partitioning (CHKPART2) is not equal to 1')
     
*       Determine assimilate requirements for growth
        ASRQ   = ASRQWD*FWD   + ASRQLV*FLV + ASRQPD*FPD +                    
     &           ASRQLRT*FLRT + ASRQTRT*FTRT     ![kg CH2O kg-1 DW]

*       Distribute remaining reserve pool (WRES2) to net organ growth
        IF (WRES1 .EQ. 0.) THEN
          WRES2 = 0.                             ![kg CH2O ha-1]
          DRES2 = 0.                             ![kg CH2O ha-1 d-]
          GTOT2 = 0.                             ![kg DW ha-1 d-1]
        ELSE 
          WRES2 = WRES1 - DRES1*DELT             ![kg CH2O ha-1]
          DRES2 = WRES2 / TAU                    ![kg CH2O ha-1 d-1]
          GTOT2  = DRES2 / ASRQ                  ![kg DW ha-1 d-1]
        ENDIF

*       Total decrease rate of reserves
        DRES   = DRES1 + DRES2                   ![kg CH2O ha-1 d-1]

*       Determine total growth of organs and plant as a total
        GWD    = GWD1 + FWD * GTOT2              ![kg DW ha-1 d-1]
        GLV    = GLV1 + FLV * GTOT2              ![kg DW ha-1 d-1]
        GPD    = GPD1 + FPD * GTOT2              ![kg DW ha-1 d-1]
        GTRT   = FTRT * GTOT2                    ![kg DW ha-1 d-1]
        GLRT   = GWURT1 + GLRT1 + FLRT * GTOT2   ![kg DW ha-1 d-1]
        GRT    = GTRT + GLRT                     ![kg DW ha-1 d-1]
        GTOT   = GTOT1 + GTOT2                   ![kg DW ha-1 d-1]

*       Death rates
        DLRT = TOWURT + TOLRT                    ![kg DW ha-1 d-1]
        DWD  = TOWD                              ![kg DW ha-1 d-1]
              
*       Leaf death rate
        CALL LEAF (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,AGE,NPL,
     &             TRMIS,GLV,PCEWMN,WTOTPP,
     &             WLV, DLV, DLV1, DLV2, LAI(1))

*       Pod yield
        CALL POD (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,TMAV,
     &           GPD, WPD, YLDPD, IPOD, BHYLD)

*       Commercial bean yield (YLDBN) is determined by the fractions of
*       beans per pod (FBEANS), a factor of loss due to the 
*       length of the fermentation process (FMTLOS) and a factor which
*       accounts for the moisture content after drying (MOISTC).
        FMTLOS = FMTDUR * FMTA + FMTB            ![-]
        YLDBN  = YLDPD * FBEANS * FMTLOS * (1.+ MOISTC) 
	                                    ![kg fermented beans ha-1 d-1]
        HARPODS= YLDBN * PODVALUE                ![ha-1 d-1]

*       Calculation of 10-day and annual totals

*       10-day totals are set to zero after each 10 day period
        IF ((TIME-STTIME) .EQ. (D10COUNT*10.+1.)) THEN
          D10GTOT  = 0.                          ![kg DW ha-1 10d-1]
          D10YLDPD = 0.                          ![kg DW ha-1 10d-1]
          D10YLDBN = 0.                          ![kg DW ha-1 10d-1]
          D10RAIN  = 0.                          ![mm 10d-1]
          D10RDD   = 0.                          ![mm 10d-1]
          D10HARPD = 0.                          ![10d-1]
          D10COUNT = D10COUNT + 1.
        ENDIF

*       10-day totals for selected output parameters
        IF (TIME.NE.STTIME) THEN
          D10GTOT  = INTGRL (D10GTOT, GTOT, DELT)  ![kg DW ha-1 10d-1]
          D10YLDPD = INTGRL (D10YLDPD, YLDPD, DELT)![kg DW ha-1 10d-1]
          D10YLDBN = INTGRL (D10YLDBN, YLDBN, DELT)![kg DW ha-1 10d-1]
          D10RAIN  = INTGRL (D10RAIN, RAIN, DELT)  ![mm d10-1]
          D10RDD   = INTGRL (D10RDD, RDD, DELT)    ![J m-2 10d -1]
          D10HARPD = D10YLDBN * PODVALUE           ![10d-1]
        ENDIF

*       Check on rainfall: warning in case annual precipitation is low
        IF (YRAIN.LT.1000..AND.IDOY.EQ.IDOYO.AND.DOY.NE.STTIME)
     &  CALL WARNING('CASE2','WA1 - Very low annual rainfall,<1000 mm')

*       Annual totals are set to zero at the start day of simulations
        IF (IDOY.EQ.IDOYO) THEN
          YGPHOT = 0.                            ![kg CH2O ha-1 y-1]
          YGTOT  = 0.                            ![kg DW ha-1 y-1]
          YGTOTAB= 0.                            ![kg DW ha-1 y-1]
          YYLDPD = 0.                            ![kg DW ha-1 y-1]
          YYLDBN = 0.                            ![kg DW ha-1 y-1]
          YLVD   = 0.                            ![kg DW ha-1 y-1]
          YWDD   = 0.                            ![kg DW ha-1 y-1]
          YRAIN  = 0.                            ![mm y-1]
          YTRANS = 0.                            ![mm y-1]
          YRDD   = 0.                            ![MJ m-2 y-1]
          BHSUM  = 0.                           ![kg DW fermented beans]
          LAISUM = 0.                            ![ha leaf ha-1 ground]
          IPODSUM= 0.                           ![kg DW fermented beans]
          YNLOSS = 0.                            ![kg N ha-1 yr-1]
          YPLOSS = 0.                            ![kg P ha-1 yr-1]
          YKLOSS = 0.                            ![kg K ha-1 yr-1]
        ENDIF
            
*       Annual totals for selected output parameters
        YGPHOT = INTGRL (YGPHOT, GPHOT, DELT)    ![kg CH2O ha-1 y-1]
        YGTOT  = INTGRL (YGTOT, GTOT, DELT)      ![kg DW ha-1 y-1]
        YGTOTAB= INTGRL (YGTOTAB, GTOT-GRT, DELT)  ![kg DW ha-1 y-1]
        YYLDPD = INTGRL (YYLDPD, YLDPD, DELT)    ![kg DW ha-1 y-1]
        YYLDBN = INTGRL (YYLDBN, YLDBN, DELT)    ![kg DW ha-1 y-1]
        YLVD   = INTGRL (YLVD, DLV, DELT)        ![kg DW ha-1 y-1]
        YWDD   = INTGRL (YWDD, DWD, DELT)        ![kg DW ha-1 y-1]
        YRAIN  = INTGRL (YRAIN, RAIN, DELT)      ![mm y-1]
        YRDD   = INTGRL (YRDD, (RDD/1000000), DELT)  ![MJ m-2 y-1]
        YTRANS = INTGRL (YTRANS, ATRANS, DELT)   ![mm y-1]

*       Variables derived from annual yield
        YHARPD = YYLDBN * PODVALUE               ![y-1]
        YNLOSS = YYLDBN * NCONTBN                ![kg N ha-1 yr-1]
        YPLOSS = YYLDBN * PCONTBN                ![kg P ha-1 yr-1]
        YKLOSS = YYLDBN * KCONTBN                ![kg K ha-1 yr-1]

*       Rain and radiation use effiency
        YRDEFF  = YYLDBN / YRDD                ![kg DW ha -1 (MJ m-2)-1]
        YRNEFF  = YYLDBN / NOTNUL(YRAIN)         ![kg DW ha -1 mm]

*       Harvest index annual value
        YHI     = YYLDBN / (WTOT - WRT)          ![-]

*       HINCR is Harvest increment following Cannell 1985
        YHINCR  = YYLDBN / NOTNUL(YGTOTAB)       ![-]

*       Annual average of riping time (IPOD)
        IPODSUM = IPODSUM + IPOD * YLDBN         ![d]
        YMNIPOD = IPODSUM / NOTNUL(YYLDBN)       ![d]

*       Annual average of butter hardness
        BHSUM  = BHSUM + BHYLD * YLDBN          ![kg DW fermented beans]
        YMNBH  = BHSUM / NOTNUL(YYLDBN)          ![-]

*       Annual average of cocoa LAI
        LAISUM = LAISUM + LAI(1)                ![ha leaves ha-1 ground]
        YMNLAI = LAISUM / DOY                   ![ha leaves ha-1 ground]

*       Finish condtions
        IF (AGEYR .GE. 40.) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2',
     &	'WA2 - Maximum tree age (40 yr) is reached')
          TERMNL = .TRUE.
        ENDIF
        IF (WTOTPP  .GE. 70.) THEN
          CALL WARNING ('CASE2',
     &	'WA3 - Maximum tree size (70 kg DW) is reached')
          TERMNL = .TRUE.
        ENDIF
        IF (WRES .LT.    0.) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2','WA4 - Reserves depleted')
          TERMNL = .TRUE.
        ENDIF

*       Output of 
        IF (OUTPUT) THEN

*       Chart output for all three output types (OUTPUTFQ)
           CALL ChartOutputRealScalar('YEAR'  , REAL(IYEAR))
           CALL ChartOutputRealScalar('AGEYR',AGEYR)
           CALL ChartOutputRealScalar('WTOT',WTOT)
           CALL ChartOutputRealScalar('WTOTPP',WTOTPP)
           CALL ChartOutputRealScalar('WRT',WRT)
           CALL ChartOutputRealScalar('WTRT',WTRT)
           CALL ChartOutputRealScalar('WLRT',WLRT)
           CALL ChartOutputRealScalar('WLV',WLV)
           CALL ChartOutputRealScalar('WWD',WWD)
           CALL ChartOutputRealScalar('WPD',WPD)
           CALL ChartOutputRealScalar('WPDCUM',WPDCUM)
           CALL ChartOutputRealScalar('WBNCUM',WBNCUM)

          IF (OUTPUTFQ.EQ.1) THEN
*           For annual output
            CALL ChartOutputRealScalar('YRDD',YRDD)
            CALL ChartOutputRealScalar('YRAIN',YRAIN)
            CALL ChartOutputRealScalar('YTRANS',YTRANS)
            CALL ChartOutputRealScalar('YGPHOT',YGPHOT)
            CALL ChartOutputRealScalar('YGTOT',YGTOT)
            CALL ChartOutputRealScalar('YLVD',YLVD)
            CALL ChartOutputRealScalar('YWDD',YWDD)
            CALL ChartOutputRealScalar('YYLDPD',YYLDPD)
            CALL ChartOutputRealScalar('YYLDBN',YYLDBN)
            CALL ChartOutputRealScalar('YHI', YHI)
            CALL ChartOutputRealScalar('YHINCR', YHINCR)
            CALL ChartOutputRealScalar('YMNBH',YMNBH)
            CALL ChartOutputRealScalar('YMNLAI',YMNLAI)
            CALL ChartOutputRealScalar('YMNIPOD', YMNIPOD)
            CALL ChartOutputRealScalar('YHARPD', YHARPD)
            CALL ChartOutputRealScalar('YRDEFF', YRDEFF)
            CALL ChartOutputRealScalar('YRNEFF', YRNEFF)
            CALL ChartOutputRealScalar('YNLOSS', YNLOSS)
            CALL ChartOutputRealScalar('YPLOSS', YPLOSS)
            CALL ChartOutputRealScalar('YKLOSS', YKLOSS)
          ELSE IF (OUTPUTFQ.EQ.2) THEN
*           For ten-day output
            CALL ChartOutputRealScalar('DOY'   , DOY  )
            CALL ChartOutputRealArray('LAI',LAI,1,1,1,1)
            CALL ChartOutputRealScalar('TMAV'  , TMAV)
            CALL ChartOutputRealScalar('D10RDD',D10RDD)
            CALL ChartOutputRealScalar('D10RAIN',D10RAIN)
            CALL ChartOutputRealScalar('D10YLDPD',D10YLDPD)
            CALL ChartOutputRealScalar('D10YLDBN',D10YLDBN)
            CALL ChartOutputRealScalar('D10GTOT',D10GTOT)
            CALL ChartOutputRealScalar('PCEW',PCEW)
          ELSE IF (OUTPUTFQ.EQ.3) THEN
*           For daily output
            CALL ChartOutputRealScalar('DOY'   , DOY  )
            CALL ChartOutputRealScalar('RAIN'  , RAIN ) 
            CALL ChartOutputRealScalar('RDD'   , RDD )
            CALL ChartOutputRealArray('LAI',LAI,1,1,1,1)
            CALL ChartOutputRealScalar('TMAV'  , TMAV)
            CALL ChartOutputRealArray('FRABS',FRABS,1,1,1,1)
            CALL ChartOutputRealScalar('GPHOT',GPHOT)
            CALL ChartOutputRealScalar('MAINT',MAINT)
            CALL ChartOutputRealScalar('GTOT',GTOT)
            CALL ChartOutputRealScalar('GTOT1',GTOT1)
            CALL ChartOutputRealScalar('GTOT2',GTOT2)
            CALL ChartOutputRealScalar('GRT',GRT)
            CALL ChartOutputRealScalar('GTRT',GTRT)
            CALL ChartOutputRealScalar('GLRT',GLRT) 
            CALL ChartOutputRealScalar('GLV',GLV )
            CALL ChartOutputRealScalar('GWD',GWD)
            CALL ChartOutputRealScalar('GPD',GPD)
            CALL ChartOutputRealScalar('FTRT',FTRT)
            CALL ChartOutputRealScalar('FLRT',FLRT)
            CALL ChartOutputRealScalar('FLV',FLV)
            CALL ChartOutputRealScalar('FWD',FWD)
            CALL ChartOutputRealScalar('FPD',FPD)
            CALL ChartOutputRealScalar('DLRT',DLRT)
            CALL ChartOutputRealScalar('DWD',DWD)
            CALL ChartOutputRealScalar('DLV',DLV)
            CALL ChartOutputRealScalar('DLV1',DLV1)
            CALL ChartOutputRealScalar('DLV2',DLV2)
            CALL ChartOutputRealArray('WWURT',WWURT,1,NLA,1,NLA)
            CALL ChartOutputRealArray('AWURT',AWURT,1,NLA,1,NLA)
            CALL ChartOutputRealScalar('ATWURT',ATWURT)
            CALL ChartOutputRealScalar('LTWURT',LTWURT) 
            CALL ChartOutputRealScalar('WTWURT',WTWURT )
            CALL ChartOutputRealArray('TRWL',TRWL,1,NLA,1,NLA)
            CALL ChartOutputRealScalar('PTRANS',PTRANS)
            CALL ChartOutputRealScalar('ATRANS',ATRANS)
            CALL ChartOutputRealScalar('PCEW',PCEW)
            CALL ChartOutputRealScalar('YLDBN',YLDBN)
            CALL ChartOutputRealScalar('YLDPD',YLDPD)
            CALL ChartOutputRealScalar('BHYLD',BHYLD)
            CALL ChartOutputRealScalar('HARPODS',HARPODS)
          END IF

*          States

           CALL OUTDAT (2, 0, 'AGE'   , AGE)
           CALL OUTDAT (2, 0, 'AGEYR'   , AGEYR)
           CALL OUTDAT (2, 0, 'WTOTCUM', WTOTCUM)
           CALL OUTDAT (2, 0, 'WTOT'   , WTOT)
           CALL OUTDAT (2, 0, 'WTOTPP'   , WTOTPP)
           CALL OUTDAT (2, 0, 'WRT'   , WRT)
           CALL OUTDAT (2, 0, 'WTRT', WTRT)
           CALL OUTDAT (2, 0, 'WLRT', WLRT)
           CALL OUTDAT (2, 0, 'WLV'  , WLV)
           CALL OUTDAT (2, 0, 'WWD'   , WWD)
           CALL OUTDAT (2, 0, 'WPD'   , WPD)
           CALL OUTDAT (2, 0, 'WPDCUM'  , WPDCUM)
           CALL OUTDAT (2, 0, 'WRES'  , WRES)
           CALL OUTDAT (2, 0, 'WRES1'  , WRES1)
           CALL OUTDAT (2, 0, 'WRES2'  , WRES2)
           CALL OUTDAT (2, 0, 'WLVD'  , WLVD)
           CALL OUTDAT (2, 0, 'WWDD'  , WWDD)
           CALL OUTDAT (2, 0, 'WLRTD'  , WLRTD)
           CALL OUTDAT (2, 0, 'TNASS' , TNASS)
           CALL OUTARR ('LAI'   , LAI, 1, INS)
           CALL OUTDAT (2, 0, 'GAI'   , GAI)
           CALL OUTARR ('FRABS' , FRABS, 1, 1)

*          Driving variables and rates
           CALL OUTDAT (2, 0, 'TMAV'  , TMAV)
           CALL OUTDAT (2, 0, 'TMAVD' , TMAVD)
           CALL OUTDAT (2, 0, 'GPHOT' , GPHOT)
           CALL OUTDAT (2, 0, 'MAINT' , MAINT)
           CALL OUTDAT (2, 0, 'GRES'  , GRES)
           CALL OUTDAT (2, 0, 'DRES'  , DRES)
           CALL OUTDAT (2, 0, 'DRES1'  , DRES1)
           CALL OUTDAT (2, 0, 'DRES2'  , DRES2)
           CALL OUTDAT (2, 0, 'MINRES', MINRES)
           CALL OUTDAT (2, 0, 'ASRQ'  , ASRQ)
           CALL OUTDAT (2, 0, 'GTOT'   , GTOT)
           CALL OUTDAT (2, 0, 'GTOT1'   , GTOT1)
           CALL OUTDAT (2, 0, 'GTOT2'   , GTOT2)
           CALL OUTDAT (2, 0, 'GRT'   , GRT)
           CALL OUTDAT (2, 0, 'GTRT', GTRT)
           CALL OUTDAT (2, 0, 'GLRT', GLRT)
           CALL OUTDAT (2, 0, 'GWURT1'   , GWURT1)
           CALL OUTDAT (2, 0, 'GWD'   , GWD)
           CALL OUTDAT (2, 0, 'GLV'   , GLV)
           CALL OUTDAT (2, 0, 'GLV1'   , GLV1)
           CALL OUTDAT (2, 0, 'GPD'   , GPD)
           CALL OUTDAT (2, 0, 'GPD1'   , GPD1)
           CALL OUTDAT (2, 0, 'DLRT', DLRT)
           CALL OUTDAT (2, 0, 'DWD'   , DWD)
           CALL OUTDAT (2, 0, 'DLV'   , DLV)
           CALL OUTDAT (2, 0, 'DLV1'   , DLV1)
           CALL OUTDAT (2, 0, 'DLV2'   , DLV2)
           CALL OUTDAT (2, 0, 'TOREQ'   , TOREQ)
           CALL OUTDAT (2, 0, 'TOWURT'   , TOWURT)
           CALL OUTDAT (2, 0, 'TOLV'   , TOLV)
           CALL OUTDAT (2, 0, 'TOPD'   , TOPD)
           CALL OUTDAT (2, 0, 'FTRT'   , FTRT)
           CALL OUTDAT (2, 0, 'FLRT'   , FLRT)
           CALL OUTDAT (2, 0, 'FLV'   , FLV)
           CALL OUTDAT (2, 0, 'FWD'   , FWD)
           CALL OUTDAT (2, 0, 'FPD'   , FPD)
           CALL OUTDAT (2, 0, 'YLDBN' , YLDBN)
           CALL OUTDAT (2, 0, 'YLDPD' , YLDPD)
           CALL OUTDAT (2, 0, 'HARPODS' ,HARPODS)
           CALL OUTDAT (2, 0, 'WBNCUM', WBNCUM)
           CALL OUTDAT (2, 0, 'PENMAN', PENMAN)
           CALL OUTDAT (2, 0, 'CROPF' , CROPF)

*          Total year values
           CALL OUTDAT (2, 0, 'YRDD'  , YRDD)
           CALL OUTDAT (2, 0, 'YRAIN' , YRAIN)
           CALL OUTDAT (2, 0, 'YTRANS', YTRANS)
           CALL OUTDAT (2, 0, 'YGPHOT', YGPHOT)
           CALL OUTDAT (2, 0, 'YGTOT', YGTOT)
           CALL OUTDAT (2, 0, 'YGTOTAB', YGTOTAB)
           CALL OUTDAT (2, 0, 'YLVD', YLVD)
           CALL OUTDAT (2, 0, 'YWDD', YWDD)
           CALL OUTDAT (2, 0, 'YYLDPD', YYLDPD)
           CALL OUTDAT (2, 0, 'YYLDBN', YYLDBN)
           CALL OUTDAT (2, 0, 'YHI', YHI)
           CALL OUTDAT (2, 0, 'YHINCR', YHINCR)
           CALL OUTDAT (2, 0, 'YMNBH', YMNBH)
           CALL OUTDAT (2, 0, 'YMNLAI', YMNLAI)
           CALL OUTDAT (2, 0, 'YMNIPOD', YMNIPOD)
           CALL OUTDAT (2, 0, 'YHARPD', YHARPD)
           CALL OUTDAT (2, 0, 'YRDEFF', YRDEFF)
           CALL OUTDAT (2, 0, 'YRNEFF', YRNEFF)
           CALL OUTDAT (2, 0, 'YNLOSS', YNLOSS)
           CALL OUTDAT (2, 0, 'YPLOSS', YPLOSS)
           CALL OUTDAT (2, 0, 'YKLOSS', YKLOSS)

*          Ten-day period values
           CALL OUTDAT (2, 0, 'D10RDD'  , D10RDD)
           CALL OUTDAT (2, 0, 'D10RAIN' , D10RAIN)
           CALL OUTDAT (2, 0, 'D10YLDPD', D10YLDPD)
           CALL OUTDAT (2, 0, 'D10YLDBN', D10YLDBN)
           CALL OUTDAT (2, 0, 'D10GTOT', D10GTOT)

*          Parameters from subroutine ROOT             
           CALL OUTDAT (2,0, 'NLA', REAL(NLA))
           CALL OUTDAT (2,0, 'CUMTKL', CUMTKL)
           CALL OUTARR ('WWURT', WWURT, 1, NLA)
           CALL OUTARR ('AWURT', AWURT, 1, NLA)
           CALL OUTARR ('WSERT', WSERT, 1, NLA)
           CALL OUTDAT (2,0, 'ATWURT', ATWURT)
           CALL OUTDAT (2,0, 'LTWURT', LTWURT)
           CALL OUTDAT (2,0, 'WTWURT', WTWURT)
           CALL OUTDAT (2,0, 'LTRT', LTRT)

*          Subroutine WUPT
           CALL OUTARR ('TRWL', TRWL, 1, NLA) 
           CALL OUTDAT (2,0, 'PTRANS', PTRANS)
           CALL OUTDAT (2,0, 'ATRANS', ATRANS)
           CALL OUTDAT (2, 0, 'PINT'  , PINT)
           CALL OUTDAT (2, 0, 'PCEW'  , PCEW)
           CALL OUTDAT (2, 0, 'PCEWMN'  , PCEWMN)

*          Subroutine POD
           CALL OUTDAT (2,0, 'BHYLD',BHYLD) 

           CALL OUTDAT (2, 0, 'DOY'   , DOY  )
           CALL OUTDAT (2, 0, 'YEAR'  , REAL(IYEAR))
           CALL OUTDAT (2, 0, 'RAIN'  , RAIN )
           CALL OUTDAT (2, 0, 'RDD'   , RDD )
           CALL OUTDAT (2, 0, 'ETRD'  , ETRD )
           CALL OUTDAT (2, 0, 'ETAE'  , ETAE )
           CALL OUTDAT (2, 0, 'EVSC'  , EVSC )   

         END IF

      ELSE IF (ITASK.EQ.3) THEN
*       -------------------
*       Integration section
*       -------------------
     
        CALL ROOT (ITASK,IUNITD,IUNITL,FILEI1,DELT,NL, NLA, NLXM, NLBM, 
     &             I1,TKL, CUMTKL, WCLQT, WCWPX, AGE, NPL, FWURT,
     &             WLRT, WTRT, LTRT,  
     &             WWURT, WTWURT, LTWURT, ATWURT, AWURT, WSERT)
        
        CALL LEAF (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,AGE,NPL,
     &             TRMIS,GLV,PCEWMN,WTOTPP,
     &             WLV, DLV, DLV1, DLV2, LAI(1))

        CALL POD (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,TMAV,
     &           GPD, WPD, YLDPD, IPOD, BHYLD)

*       Green Area Index of shade and cocoa trees
        GAI = 0.                                 ![ha leaf ha-1 ground]
        DO I1=1,INS
           GAI    = GAI + LAI(I1)                ![ha leaf ha-1 ground]
        ENDDO    

*       Dry matter
        WLV    = MAX(0.,INTGRL (WLV, GLV-DLV, DELT))      ![kg DW ha-1]
        WWD    = INTGRL (WWD , GWD-DWD, DELT)    ![kg DW ha-1]
        WTRT   = INTGRL (WTRT, GTRT, DELT)       ![kg DW ha-1]
        WLRT   = INTGRL (WLRT, GLRT-DLRT, DELT)  ![kg DW ha-1]
        WRT    = WLRT + WTRT                     ![kg DW ha-1]
        WPD    = MAX(0.,INTGRL (WPD , GPD-YLDPD, DELT))   ![kg DW ha-1]
        WTOT   = WWD   + WLV  + WPD + WTRT + WLRT![kg DW ha-1]
        WTOTPP = WTOT/NPL                        ![kg DW per plant]
        WTOTCUM= INTGRL (WTOTCUM, GTOT, DELT)    ![kg DW ha-1]

*       Dead dry matter
        WLVD   = INTGRL (WLVD, DLV , DELT)       ![kg DW ha-1]
        WWDD   = INTGRL (WWDD, DWD , DELT)       ![kg DW ha-1]
        WLRTD  = INTGRL (WLRTD , DLRT, DELT)     ![kg DW ha-1]

*       Reserves
        WRES   = INTGRL (WRES, GRES-DRES, DELT)  ![kg CH2O ha-1]
        MINRES = MINCON * WTOT

*       Cumulative harvested pods and beans
        WPDCUM = INTGRL (WPDCUM, YLDPD, DELT)    ![kg DW ha-1]
        WBNCUM = INTGRL (WBNCUM, YLDBN, DELT)    ![kg DW ha-1]

*       Carbon balance check - for entire tree
        TNASS  = INTGRL (TNASS,((GPHOT - MAINT - GRES + DRES)
     &                  *44./30.) -
     &                  (GLRT*CO2LRT + GTRT*CO2TRT + GLV*CO2LV            
     &                  + GWD*CO2WD + GPD*CO2PD),DELT)    ![kg CO2 ha-1]
        CHKFL  = TNASS * (12./44.)                        ![kg C ha-1]
        CHKIN  = (WLV + WLVD - WLVI)*CFLV + (WWD + WWDD - WWDI)*CFWD
     &         + (WTRT - WTRTI)*CFTRT + (WLRT + WLRTD - WLRTI)*CFLRT 
     &         + (WPD + WPDCUM - WPDI)*CFPD               ![kg C ha-1]

        CHKDIF = (CHKIN-CHKFL)/NOTNUL(CHKIN)     ![-]

*       Carbon balance check - for plant organs separately

        CHKLV  = GLV*CFLV   + GLV*CO2LV*12./44.  - GLV*ASRQLV*12./30. 
        CHKWD  = GWD*CFWD   + GWD*CO2WD*12./44.  - GWD*ASRQWD*12./30.  
        CHKLRT = GLRT*CFLRT + GLRT*CO2LRT*12./44.- GLRT*ASRQLRT*12./30.
        CHKTRT = GTRT*CFTRT + GTRT*CO2TRT*12./44.- GTRT*ASRQTRT*12./30.
        CHKPD  = GPD*CFPD   + GPD*CO2PD*12./44.  - GPD*ASRQPD*12./30.
                                                 !all [kg C ha-1]

*       Terminate program if CHKDIF is too large
        IF (CHKIN.GT.10. .AND. ABS(CHKDIF) .GT. 0.01) CALL FATALERR
     &     ('CASE2','	FE21 - Carbon balance is wrong
     &     CHKDIF larger than 0.01')

*       Terminate program if CHK** for one of the organs is too large
        IF ((ABS(CHKLV).GT.0.001).OR.(ABS(CHKWD).GT.0.001).OR.
     &     (ABS(CHKLRT).GT.0.001).OR.(ABS(CHKTRT).GT.0.001).OR.
     &     (ABS(CHKPD) .GT.0.001)) CALL FATALERR
     &     ('CASE2','	FE22 - Carbon balance is wrong for one organ, 
     &     CHK** larger than 0.001')

      ELSE IF (ITASK.EQ.4) THEN
*       ----------------
*       Terminal section
*       ----------------
*       No tasks defined for terminal section

      END IF

      RETURN
      END


*     Subroutines

*----------------------------------------------------------------------*
* SUBROUTINE LEAF                                                      *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Authors: Wouter Gerritsma                                            *
*          Pieter Zuidema (as of April 2001)                           *
* Date   : May 1995                                                    *
* Purpose: This subroutine simulates leaf growth and senescence of     *
*          leaf age classes based on the boxcar train without          *
*          dispersion delay techinque                                  *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEI1  C*  Name of file with plant data                  -       I  *
* TERMNL  L4  Flag to indicate if simulation is to stop      -     I/O *
* DELT    R4  Time step of integration                       d      I  *
* AGE     R4  Age of trees                                   d      I  *
* NPL     R$  Density of cacao trees                       ha -1    I  *
* TRMIS   R4  Proportion of light transmitted by shade trees -      I  *
* GLV     R4  Growth rate of leaves                        kg ha-1  I  *
* PCEWMN  R4  Reduction factor for photosynthesis           -       I  *
* WTOTPP  R4  Total dry weight per plant                   kg       I  *
* WLV     R4  Weight of leaves                             kg ha-1  O  *
* DLV     R4  Death rate of leaves                         kg ha-1  O  *
* DLV1    R4  Death rate of leaves, due to ageing          kg ha-1  O  *
* DLV2    R4  Death rate of leaves, due to drought         kg ha-1  O  *
* LAI     R4  Leaf Area Index                              m2 m-2   O  *
*                                                                      *
* Fatal error checks: AVGLVAGE>365; MINLVAGE<0; AVGLVAGE<MINLVAGE      *
* Warnings          : WLEAFTOT <> WLV                                  *
* Subroutines called: several from TTUTIL                              *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE LEAF (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,AGE,NPL,
     &                 TRMIS,GLV,PCEWMN,WTOTPP,
     &                 WLV, DLV, DLV1, DLV2, LAI)

      IMPLICIT NONE

*     Formal parameters
      INTEGER ITASK, IYEAR, IDOY
      LOGICAL TERMNL
      REAL DELT, AGE, TRMIS, GLV, WLV, DLV, DLV1, DLV2, RDLV, ADJLVAGE
      REAL PCEWMN, LAI
      REAL AVGLVAGE, MINLVAGE, NPL 
      REAL SLAR1A, SLAR1B, SLAR2A, SLAR2B, SLAMOD, WTOTPP
      REAL CHKTLV, WLEAFTOT

*     Standard local declarations
      REAL WLEAF(0:365), LA(0:365), SLA(0:365)
      INTEGER ITABLE, I1, ILD, IPLD
      PARAMETER (ITABLE = 100)
      INTEGER   IUNITD, IUNITL
      CHARACTER FILEI1*80

*     Real function declaration from TTUTIL
      REAL LINT, NOTNUL

                              
      IF (ITASK.EQ.1) THEN
*       ----------------------
*       Initialization section
*       ----------------------
        CALL RDINIT (IUNITD  , IUNITL, FILEI1)
        CALL RDSREA ('AVGLVAGE', AVGLVAGE)       ![d]
        CALL RDSREA ('MINLVAGE', MINLVAGE)       ![d]
        CALL RDSREA ('SLAR1A', SLAR1A)           
	                                  ![ha kg-1 leaf DW kg-1 plant DW]
        CALL RDSREA ('SLAR1B', SLAR1B)           ![ha kg-1 leaf DW]    
        CALL RDSREA ('SLAR2A', SLAR2A)           
	                                  ![ha kg-1 leaf DW kg-1 plant DW]
        CALL RDSREA ('SLAR2B', SLAR2B)           ![ha kg-1 leaf DW]    
        CLOSE (IUNITD)

*       Leaf ages  
        ILD = INT(AVGLVAGE)                      ![-]
        ADJLVAGE = NINT(MINLVAGE - PCEWMN*MINLVAGE + PCEWMN*AVGLVAGE)
	                                           ![d]

*       Calculate SLA modifier to account for higher SLA in shaded 
*       conditions. Only when light transmission below shade tree 
*       canopy <0.80 (at higher values there is no effect) 
        IF (TRMIS.LE.0.80) THEN
	    SLAMOD = SLAR2A * LOG(TRMIS) + SLAR2B
        ELSE
	    SLAMOD = 1.
        ENDIF

*       Initialise boxcar train with leaf weights and leaf areas
        LAI   = 0.
        WLEAFTOT= 0.
        DO 10 I1 = 1,ILD
          WLEAF(I1)  = WLV/AVGLVAGE              ![kg DW leaf]
          SLA(I1)    = (SLAR1A * WTOTPP + SLAR1B) * SLAMOD      
		                                       ![ha leaf kg-1 DW leaf]
          LA(I1)     = WLEAF(I1) * SLA(I1)       ![ha leaf ha-1 ground]
          LAI        = LAI + LA(I1)              ![ha leaf ha-1 ground]
          WLEAFTOT   = WLEAFTOT +WLEAF(I1)       ![kg DW leaf]
10      CONTINUE

*       Fatal error checks
        IF (AVGLVAGE.GT.365.) CALL FATALERR
     &     ('CASE2','FE23 - Maximum leaf age greater than one year')
        IF (MINLVAGE.LT.0.) CALL FATALERR
     &     ('CASE2','FE24 - Minimum leaf age less than zero (0)')
        IF (AVGLVAGE.LT.MINLVAGE) CALL FATALERR
     &     ('CASE2','FE25 - Maximum leaf age less than min. leaf age')

*       Check whether leaf weight in boxcars equals leaf weight
        CHKTLV = NOTNUL(WLV) / NOTNUL(WLEAFTOT)
        IF ((CHKTLV.LT.0.95) .OR. (CHKTLV.GT.1.05)) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2 - LEAF',
     &    'WA5 - Sum leaf class weights not equal to total leaf weight')
          TERMNL = .TRUE.
        ENDIF

      ELSE IF (ITASK.EQ.2) THEN
*       ------------------------
*       Rate calculation section
*       ------------------------

*       Characteristics of new leaves
        WLEAF(0)  = GLV * DELT                   ![kg DW ha-1]
        SLA(0)    = (SLAR1A * WTOTPP + SLAR1B) * SLAMOD 
	                                           ![ha leaf kg-1 leaf DW]
        LA(0)     = WLEAF(0)*SLA(0)              ![ha leaf ha-1 ground]

*       Adjust leaf age for water stress sensitivity
        ILD = AVGLVAGE                           ![d]
        ADJLVAGE = NINT(MINLVAGE - PCEWMN*MINLVAGE + PCEWMN*AVGLVAGE)
	                                           ![d]

*       Determine leaf death rates for senescing leaves (DLV1) and 
*       leaves dying due to water shortage (DLV2)
        DLV = 0.                                 ![kg DW ha-1 d-1]
        DLV1 = 0.                                ![kg DW ha-1 d-1]
        DLV2 = 0.                                ![kg DW ha-1 d-1]
        DLV1 = WLEAF(ILD)/DELT                   ![kg DW ha-1 d-1]
        LA(ILD) = 0.                             ![ha leaf ha-1 ground]
        WLEAF(ILD) = 0.                          ![kg DW ha-1]
        RDLV = (1./ADJLVAGE - 1./AVGLVAGE)       ![d-1]
        DLV2 = RDLV * WLV * (AVGLVAGE-1)/AVGLVAGE![kg DW ha-1 d-1]
        DLV  = DLV1 + DLV2                       ![kg DW ha-1 d-1]

      ELSE IF (ITASK.EQ.3) THEN
*       ------------------------
*       Integration section
*       ------------------------
*       Shift all the leaves, weights and areas one class
*       and account for extra leaf loss due to water shortage
        LAI = 0.                                 ![ha leaf ha-1 ground]
        WLEAFTOT = 0.
        DO 30 I1=ILD-1,0,-1
          SLA(I1+1)    = SLA(I1)                 ![ha leaf kg-1 leaf DW]
          WLEAF(I1+1)  = WLEAF(I1) * (1-RDLV)    ![kg DW ha-1]
          LA(I1+1)     = WLEAF (I1) * SLA(I1)    ![ha leaf ha-1 ground]
          LAI          = LAI + LA(I1)            ![ha leaf ha-1 ground]
          WLEAFTOT     = WLEAFTOT + WLEAF(I1)    ![kg DW ha-1]
30      CONTINUE                      
        WLEAF(0)  = 0.                           ![kg DW ha-1]
        LA(0)     = 0.                           ![ha leaf ha-1 ground]
        SLA(0)    = 0.                           ![ha leaf kg-1 leaf DW]
    
        CHKTLV = NOTNUL(WLV + GLV*DELT - DLV*DELT) / NOTNUL(WLEAFTOT)
        IF (WLEAFTOT .GT. 10. .AND. ((CHKTLV.LT.0.95) .OR. 
     &	(CHKTLV.GT.1.05))) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2 - LEAF',
     &    'WA6 - Sum leaf class weights not equal to total lf weight')
          TERMNL = .TRUE.
      ENDIF

      END IF

      RETURN
      END

*----------------------------------------------------------------------*
* SUBROUTINE ROOT                                                      *
*                                                                      *
* Use:     : For CASE2 (Cacao Simulation Engine) version 2.2           *
* Author(s): Liesje Mommer                                             *
*            Pieter Zuidema (as of April 2001)                         *
* Date     : oktober 1998, Version:1.0                                 *
* Purposes:                                                            *
* To calculate the number of layers of which water can be taken up     *
* To calculate the root biomass of cacao, separated in three classes:  *
* finest roots, able to take up water; other fine roots; taproot.      *
* To calculate the distribution of the fine roots in the soil.         *
* To calculate the rootlet surface within a layer.                     *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning (unit)                                     class *
* ----   ---- ---------------                                    ----- *
* PLTMOD  C*  Name of plant module (-)                              ?  *
* ITASK   I4  Task that subroutine should perform (-)               I  *
* IUNITD  I4  Unit number that is used for input files (-)          ?  *
* IUNITL  I4  Unit number that is used for log file (-)             ?  *
* FILEI1  C*  File name with which plant parameters are read (-)    ?  *
* DELT    R4  Time interval of integration (d)                      I  *
* NL      I4  Actual number of soil compartments (-)                I  *
* NLA     I4  No of soil comp. from which water can be extracted(-) I  *
* NLBM    I4  Maximum no of soil comp.                       -      I  *
* NLXM    I4  No. of layers as declared in calling program   -      I  *
* I1      I4  DO-loop counter (<not given>)                         I  *
* TKL     R4  Thicknesses of soil compartments (m)                  I  *
* CUMTKL  R4  Cumulative thickness of rooted soil layers            O  *
* WCLQT[] R4  Volumetric soil water content per layer        (-)    I  *
* WCWPX[] R4  Vol. water content at wilting point per layer (-)     I  *
* AGE     R4  Age of tree (d)                                       I  *
* NPL     R4  Cacao tree density (d)                                I  *
* FWURT   R4  Fraction of lateral roots that may take up water (-)  I  *
* WLRT    R4  Weight of lateral roots (kg DW ha-1)                  I  *
* WTRT    R4  Weight of taproot (kg DW ha-1)                        I  *
* LTRT    R4  Length of taproot (m)                                 I  *
* WWURT   R4  Weight of water-uptaking roots per layer (kg DW ha-1) I  *
* WTWURT  R4  Total weight of water-uptaking roots     (kg DW ha-1) I  *
* LTWURT  R4  Total length of water-uptaking roots         (m ha-1) I  *
* ATWURT  R4  Total area of water-uptaking roots          (m2 ha-1) I  *
* AWURT   R4  Area of water-uptaking roots per layer      (m2 ha-1) I  *
* WSERT   R1  Auxiliary variable to calculate root extension -      O  *
*                                                                      *
* Fatal error checks: none                                             *
* Warnings          : none                                             *
* Subroutines called: TTUTIL functions                                 *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE ROOT (ITASK,IUNITD,IUNITL,FILEI1,DELT,NL,NLA,NLXM,NLBM,
     &              I1, TKL, CUMTKL, WCLQT, WCWPX, AGE, NPL, FWURT,
     &              WLRT, WTRT, LTRT,  
     &              WWURT, WTWURT, LTWURT, ATWURT, AWURT, WSERT)
     
      IMPLICIT NONE
                                    
*     Formal parameters
      INTEGER ITASK 
      REAL DELT
      INTEGER NL, NLA, NLXM, NLBM, I1
      REAL AGE, LTRT
      REAL WLRT, WTWURTR, WTRT
      REAL WTWURT
      REAL ATWURT
      REAL CUMTKL
      REAL LTRTLL 
      REAL GLRT, DLRT
      REAL AWURT(NLBM)
      REAL WWURTR(NLBM) 
      REAL WSERT(NLBM) 
      REAL TKL(NLXM), WCLQT(NLXM), WCWPX(NLXM)  

*     Local declarations          
*     Plant parameters
      REAL NPL

*     Fine root declarations 
*     Total length of roots able to take up water
      REAL LTWURT            
      
*     Regression coefficients for vertical distribution of 
*     water-uptaking roots
      REAL VDWURTRA, VDWURTRB
     
*     Fraction finest roots of all lateral roots, able to take up water,         
*     diameter less then 2 mm. 
      REAL FWURT                                                         

*     Mean diameter and specific root length of the finest root classes
*     (0-1mm and 1-2mm). 
      REAL DIAM1, DIAM2, SPRTL1, SPRTL2

*     Specific weight of wood
      REAL SW                                                           
     
*     Parameter pi
      REAL PI
      PARAMETER (PI= 3.1415927)  
      
      REAL ZERO
      PARAMETER (ZERO=0.)  

*     Array declaration 
      INTEGER NLRT
      PARAMETER (NLRT=10)
      REAL WWURT(1:NLRT), LWURT1(1:NLRT), LWURT2(1:NLRT)

      INTEGER   IUNITD, IUNITL
      CHARACTER FILEI1*80

*     Function table declaration 
*     Variable ITABLE gives the maximum length of the array 
      INTEGER ITABLE 
      PARAMETER (ITABLE = 10)              
      
*     TTUTIL functions
      REAL INTGRL, LINT
      
      IF (ITASK.EQ.1) THEN
*     ----------------------              
*     Initialisation section
*     ----------------------
      
*       Read plant parameters from file 
        CALL RDINIT (IUNITD  , IUNITL, FILEI1)
        CALL RDSREA ('VDWURTRA', VDWURTRA)       ![-]
        CALL RDSREA ('VDWURTRB', VDWURTRB)       ![kg DW ha-1 m-2]
        CALL RDSREA ('SW', SW)                   ![kg ha-1]
        CALL RDSREA ('DIAM1', DIAM1)             ![m]
        CALL RDSREA ('DIAM2', DIAM2)             ![m]
        CALL RDSREA ('SPRTL1', SPRTL1)           ![m kg-1 DW]
        CALL RDSREA ('SPRTL2',SPRTL2)            ![m kg-1 DW]
        CLOSE (IUNITD)
	          
        NLA     = 0                              ![-]        
        CUMTKL  = 0                              ![-]
        LTWURT  = 0.                             ![m ha-1]
        WTWURTR = 0.                             ![kg DW ha-1]
        ATWURT  = 0.                             ![m2 ha-1]
      
*       Arrays are set to zero. 
        DO I1=1, NLXM
          WWURT(I1)  = 0.                        ![kg DW ha-1]
          WWURTR(I1) = 0.                        ![kg DW ha-1]
          LWURT1(I1) = 0.                        ![m ha-1]
          LWURT2(I1) = 0.                        ![m ha-1]
        END DO   
       
        DO I1=1, NLBM
          AWURT(I1) = 0.                         ![m2 ha-1]
        END DO                      
      
*       Total weight of water-uptaking roots
        WTWURT = FWURT * WLRT                    ![kg DW ha-1]

*       Calculation of length of taproot, from known weight.
*       The shape of the taproot is assumed to be that of a cone.
        LTRT = ((WTRT * 1200)/(NPL * SW * PI))**0.3333    ![m]
      
*       Number of layers of which water can be taken up depends on 
*       length of taproot
        DO I1 = 1, NL
           IF (LTRT.GT.CUMTKL) THEN 
              NLA = NLA + 1                      ![-]
              CUMTKL = CUMTKL + TKL(I1)          ![m]
            ELSE 
              NLA = NLA                          ![-]
            ENDIF
        END DO


      ELSE IF (ITASK.EQ.2) THEN
*     ------------------------
*     Rate calculation section
*     ------------------------

*       Calculate vertical distribution of wateruptaking roots. 
*       First determine distribution based on regression, for all
*       except the lowest available layer.       
        CUMTKL  = 0.                             ![m]
        WTWURTR = 0.                             ![kg DW ha-1]
        LTRTLL = 0.                              ![m]
*       - in case the taproot penetrates only one soil layer (NLA=1)
        IF (NLA. EQ .1) THEN
          WWURTR(1)    = (VDWURTRB*(0.5*LTRT)**VDWURTRA)*LTRT
          WTWURTR      = WWURTR(1)
          CUMTKL       = TKL(1)                  ![-]
        ELSE
*       - in case that several layers are available
          DO I1=1,NLA-1
            CUMTKL     = CUMTKL + TKL(I1)        ![m]
            WWURTR(I1) = 0.                      ![kg DW ha-1]
            WWURTR(I1) = (VDWURTRB*(CUMTKL-0.5*TKL(I1))**VDWURTRA)            
     &                    *TKL(I1)               ![kg DW ha-1]
            WTWURTR    = WTWURTR + WWURTR(I1)    ![kg DW ha-1]
          END DO
*         - for last layer: only the extend of taproot in that layer 
          LTRTLL       = LTRT - CUMTKL           ![m]    
          WWURTR(NLA)  = (VDWURTRB*(CUMTKL+0.5*LTRTLL)**VDWURTRA)*LTRTLL
		                                       ![kg DW ha-1]
          WTWURTR      = WTWURTR + WWURTR(NLA)   ![kg DW ha-1]
          CUMTKL       = CUMTKL + TKL(NLA)       ![m]
      ENDIF

*       Then adjust total from regression to actual DW of water-uptaking roots
        WTWURT = FWURT * WLRT                    ![kg DW ha-1]
        DO I1=1, NLA        
          WWURT(I1)  = WWURTR(I1) / WTWURTR * WTWURT  ![kg DW ha-1]
        END DO
      
*       -----------
        DO I1=1, NLA 
          WSERT(I1) = 1.                         ![-]
              IF (WCLQT(I1).LT.WCWPX(I1)) THEN
              WSERT(I1) = 0.                     ![-]
              END IF
        END DO      
  
*       Determine length of water uptaking roots
        LTWURT=0.                                ![m ha-1]
        DO I1=1,NLA            
          LWURT1(I1) = 0.5 * WWURT(I1) * SPRTL1  ![m ha-1]
          LWURT2(I1) = 0.5 * WWURT(I1) * SPRTL2  ![m ha-1]
          LTWURT = LTWURT + LWURT1(I1) + LWURT2(I1) ![m ha-1]
        END DO

*       Determine area of water uptaking roots      
        ATWURT=0.
        DO I1=1, NLA
          AWURT(I1) = (PI*DIAM1)*LWURT1(I1)+                                
     &                (PI*DIAM2)*LWURT2(I1)      ![m2 ha-1]
          ATWURT = ATWURT + AWURT(I1)            ![m2 ha-1]
        END DO
      

      ELSE IF (ITASK.EQ.3) THEN
*     -------------------
*     Integration section
*     -------------------

*       Calculation of length of taproot, from known weight.
*       The shape of the taproot is assumed to be that of a cone.
        LTRT = ((WTRT * 1200)/(NPL * SW * PI))**0.3333

*       Detemine the number of soil layers available to the trees.
        DO I1 = NLA, NL
          IF (LTRT.GT.CUMTKL) THEN 
            NLA = NLA + 1                        ![-]
            CUMTKL = CUMTKL + TKL(I1)            ![m]
          ELSE
            NLA = NLA                            ![-]
          ENDIF
        END DO

      END IF
                   
      RETURN
      
      END
                                   
*----------------------------------------------------------------------*
* SUBROUTINE WUPT                                                      *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Author : Wouter Gerritsma                                            *      
* Date   : June 1995                                                   *
* Version: 1.0                                                         *
*                                                                      *
* Purpose: Calculate potential and actual transpiration, and water     *
*          water uptake from the separate soil layers                  *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEI1   C*  Name of file with plant data                  -      I  *
* NLXM    I4  no. of layers as declared in calling program   -      I  *
* NL      I4  number of layers specified in input file       -      I  *
* NLA     I4  Nr of soil comp. from which water can be extracted(-) I  *
* NLBM    I4  Maximum no of soil comp.                       -      I  *
* TKL[]   R4  thickness of soil compartments                 m      I  *
* WCLQT[] R4  volumetric soil water content per layer        -      I  *
* WCWPX[] R4  volumetric water content at wilting point      -      I  *
* WCFCX[] R4  volumetric water content at field capacity     -      I  *
* WCSTX[] R4  volumetric water content at saturation         -      I  *
* ATWURT  R4  Total area of water-uptaking roots          (m2 ha-1) I  *
* AWURT   R4  Area of water-uptaking roots per layer      (m2 ha-1) I  *
* ETRD    R4  Radiation driven part of ETPMD               mm/d     I  *
* ETAE    R4  Dryness driven part of ETPMD                 mm/d     I  *
* EVSC    R4  actual (realized) evaporation rate           mm/d     I  *
* GAI     R4  total leaf area index                          -      I  *
* ZRT     R4  rooted depth                                   -      I  *
* TRWL[]  R4  Actual transpiration rate per layer          mm/d     O  *
* PINT    R4  Rain intercepted by the canopy               mm/d     I  *
* PTRANS  R4  Potential transpiration rate                 mm/d     O  *
* ATRANS  R4  Actual transpiration rate                    mm/d     O  *
* PCEW    R4  Factor that accounts for reduced                         *
*             photosynthesis due to water stress             -      O  *
* PENMAN  R4  Penman reference value for potential                     *
*             evapotranspiration                           mm/d     O  *
* CROPF   R4  Crop factor for crop water requirement         -      O  *
*                                                                      *
* Fatal error checks: NL>NLBM                                          *
* Warnings          : none                                             *
* Subroutines called: SWSE, many from TTUTIL                           *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE WUPT (ITASK,IUNITD,IUNITL,FILEI1,NLXM,NL, NLA, NLBM,  
     &                  TKL   , WCLQT , WCWPX , WCFCX , WCSTX , ATWURT ,
     &                  AWURT, ETRD  , ETAE  ,  EVSC  , GAI   ,
     &                  TRWL  , PINT, PTRANS, ATRANS, PCEW ,
     &                  PENMAN,  CROPF )

      IMPLICIT NONE

*     Formal parameters

      INTEGER ITASK

      INTEGER NLXM, NL, NLA , NLBM
      REAL TRWL(NLXM) , TKL(NLXM)
      REAL WCLQT(NLXM), WCWPX(NLXM), WCFCX(NLXM), WCSTX(NLXM)
      
      REAL ATWURT
      REAL AWURT(NLBM)
      REAL EVSC, ETRD, ETAE, PTRANS, ATRANS, PCEW
      REAL  PINT, GAI, PENMAN, CROPF

      INTEGER   IUNITD, IUNITL
      CHARACTER FILEI1*80

*     Local declarations
      INTEGER NLLM1
      PARAMETER (NLLM1=10)
      REAL WSEL(NLLM1)
      INTEGER I1
      REAL TRANSC, WCWET, AVAIL, PRWU , P
      REAL LIMIT        
      SAVE

      IF (NL .GT. NLBM) CALL FATALERR
     &   ('CASE2','FE26 - Too many layers in external arrays')

      IF (ITASK .EQ. 1) THEN
*       ----------------------
*       Initialization section
*       ----------------------

        CALL RDINIT (IUNITD  , IUNITL, FILEI1)
        CALL RDSREA ('TRANSC', TRANSC)           ![mm d-1]
        CALL RDSREA ('WCWET' , WCWET )           ![cm3 H2O cm-3 soil]
        CLOSE (IUNITD)
    
*       Transpiration rate per layer is 'zeroed'
        DO 20 I1=1,NL
           TRWL(I1) = 0.                         ![mm d-1]
20      CONTINUE


      ELSE IF (ITASK.EQ.2) THEN
*       ------------------------
*       Rate calculation section
*       ------------------------

*       Transpiration and water uptake

*       Potential transpiration
        PTRANS = MAX(0.,(ETRD*(1.-EXP(-0.5*GAI)) + ETAE*MIN(2.0, GAI)
     &                   - 0.5*PINT))            ![mm d-1]

*       Potential water uptake rate
        PRWU = MAX (0., PTRANS / ATWURT)         ![mm m-2 root area]

*       Soil water depletion factor
        P    = TRANSC / (TRANSC + PTRANS)        ![-]            

*       Calculate actual transpiration (ATRANS) from
        ATRANS = 0.                              ![mm d-1]
        DO 50 I1 = 1, NLA
          CALL SWSE (WCLQT(I1), P, WCWET, WCWPX(I1), WCFCX(I1),
     &               WCSTX(I1), WSEL(I1))
          TRWL(I1) = PRWU * WSEL(I1) * AWURT(I1) ![mm d-1]
          AVAIL = MAX(0., (WCLQT(I1) - WCWPX(I1)) * TKL(I1) * 1000.)
		                                       ![mm]
          IF (TRWL(I1) .GT. AVAIL) TRWL(I1) = AVAIL                        
          ATRANS = ATRANS + ABS(TRWL(I1))        ![mm d-1]
50      CONTINUE

*       Calculate water availability factor
        IF (PTRANS.GT.0.) THEN
          PCEW = ATRANS/PTRANS                   ![-]
        ELSE
          ATRANS = 0.                            ![mm d-1]
          PCEW   = 1.
        END IF

*        WSERT = 1.
*        IF (WCLQT(IRGL).LT.WCWPX(IRGL)) WSERT = 0.

*       Miscellaneous water related variables
        CROPF  = (PTRANS+EVSC)/(ETRD+ETAE)       ![-]
        PENMAN = ETRD+ETAE                       ![mm d-1]
      
      END IF
      RETURN
      END

 
*----------------------------------------------------------------------*
* SUBROUTINE POD                                                       *
*                                                                      *
* Use:     For CASE2 (Cacao Simulation Engine) version 2.2             *
* Authors: Wouter Gerritsma                                            *
*          Pieter Zuidema (as of April 2001)                           *
* Date   : May 1995                                                    *
* Purpose: This subroutine simulates pod growth and development        *
*                                                                      *
* FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
* name   type meaning                                     units  class *
* ----   ---- -------                                     -----  ----- *
* ITASK   I4  Task that subroutine should perform            -      I  *
* IUNITD  I4  Unit that can be used for input files          -      I  *
* IUNITL  I4  Unit used for log  file                        -      I  *
* FILEI1   C*  Name of file with plant data                  -      I  *
* TERMNL  L4  Flag to indicate if simulation is to stop      -     I/O *
* DELT    R4  Time step of integration                       d      I  *
* TMAV    R4  Daily Average Temperature                   degree    I  *
* GPD     R4  Growth rate of pods                  kg DW ha-1 d-1   I  *
* WPD     R4  Weight of pods                           kg DW ha-1   O  *
* YLDPD   R4  Weight of harvested pods             kg DW ha-1 d-1   O  *
* IPOD    R4  Number of pod classes                          -      O  *
* BHYLD   R4  Butter hardness of harvested beans             -      O  *
*                                                                      *
* Fatal error checks: none                                             *
* Warnings          : GPODTOT <> GPD, WPODTOT <> WPD                   *
* Subroutines called: TTUTIL functions                                 *
* File usage        : IUNITD                                           *
*----------------------------------------------------------------------*
      SUBROUTINE POD (ITASK,IUNITD,IUNITL,FILEI1,TERMNL,DELT,TMAV,
     &                GPD, WPD, YLDPD, IPOD, BHYLD)

      IMPLICIT NONE

*     Formal parameters
      LOGICAL TERMNL
      INTEGER ITASK, IYEAR, IDOY
      REAL DELT
      REAL WPD, GPD, YLDPD
      REAL TSS, PGRUSS
      REAL GPODTOT, WPODTOT
      REAL CHKGPOD, CHKWPOD

*     Development rates
      REAL DEVR, DEVRH, DEVRL, TMAV
      REAL DEVRR1A, DEVRR1B, DEVRR2A, DEVRR2B

*     Real function declaration from TTUTIL
      REAL LINT, NOTNUL

*     Butter hardness 
      REAL BHRA, BHRB
      REAL BHYLD

*     Standard local declarations
      REAL WPOD(0:200),STAGE(0:200),SSPOD(0:200),GPOD(0:200),BH(0:200)
      INTEGER ITABLE, I1, IPOD

      PARAMETER (ITABLE = 100)
      REAL      SSTB(ITABLE)
      INTEGER   ISSN        
      INTEGER   IUNITD, IUNITL
      CHARACTER FILEI1*80


      IF (ITASK.EQ.1) THEN
*       ----------------------
*       Initialization section
*       ----------------------
        CALL RDINIT (IUNITD  , IUNITL, FILEI1)
        CALL RDSREA ('DEVRR1A',DEVRR1A)          ![d-1 degr C -1]
        CALL RDSREA ('DEVRR1B',DEVRR1B)          ![d-1]
        CALL RDSREA ('DEVRR2A',DEVRR2A)          ![d-1 degr C -1]
        CALL RDSREA ('DEVRR2B',DEVRR2B)          ![d-1]
        CALL RDAREA ('SSTB'  , SSTB  ,  ITABLE, ISSN ) ![-]
        CALL RDSREA ('BHRA'   , BHRA)            ![-]
        CALL RDSREA ('BHRB'   , BHRB)            ![-]
        CLOSE (IUNITD)

*       Determine inital number of boxcars
        DEVRL = MAX(0., DEVRR1A * TMAV + DEVRR1B)![d-1]     
        DEVRH = MAX(0., DEVRR2A * TMAV + DEVRR2B)![d-1]
        IPOD  = INT(1/MIN(DEVRL,DEVRH)) + 1      ![d]
      
*       Initialise stage distributions
        TSS = 0.                                 ![-]
        DO 10 I1 = 1,IPOD
          STAGE(I1) = REAL(I1)/REAL(IPOD)        ![-]
          WPOD(I1)  = 0.                         ![kg DW ha-1]
          GPOD(I1)  = 0.                         ![kg DW ha-1 d-1]
          SSPOD(I1) = 0.                         ![-]
          IF (STAGE(I1).LT.1.) THEN
            SSPOD(I1) = LINT(SSTB,ISSN,STAGE(I1))![-]
            TSS = TSS + SSPOD(I1)                ![-]
          ELSEIF (STAGE(I1).GE.1.) THEN
            STAGE(I1) = 0.                       ![-]
          ENDIF
10      CONTINUE

*       Initialise boxcar train with pod weights
        DO 20 I1 = 1,IPOD
          WPOD(I1) = WPD * SSPOD(I1)/TSS         ![kg DW ha-1]
          BH(I1)   = 1.6                         ![-]
20      CONTINUE
        STAGE(0) = 0.                            ![-]
        WPOD(0)  = 0.                            ![kg DW ha-1]
        GPOD(0)  = 0.                            ![kg DW ha-1 d-1]

      ELSE IF (ITASK.EQ.2) THEN
*       ------------------------
*       Rate calculation section
*       ------------------------
*       Development rate
        DEVRL = MAX(0., DEVRR1A * TMAV + DEVRR1B)![d-1]     
        DEVRH = MAX(0., DEVRR2A * TMAV + DEVRR2B)![d-1]
        DEVR  = MIN(DEVRL,DEVRH)                 ![d-1]     

*       Growth per unit sink stregth
        PGRUSS = GPD/TSS                         ![kg DW ha-1 d-1]

*       Growth rate of pods per class
        GPODTOT = 0.                             ![kg DW ha-1 d-1]
        DO 30 I1 = 1,IPOD
          IF (SSPOD(I1) .GT. 0.) THEN
            GPOD(I1) = SSPOD(I1)*PGRUSS          ![kg DW ha-1 d-1]
          ELSE
            GPOD(I1) = 0.                        ![kg DW ha-1 d-1]
          ENDIF
        GPODTOT = GPODTOT + GPOD(I1)             ![kg DW ha-1 d-1]
30      CONTINUE

*       Check whether growth in pod classes equals overall pod growth
*       Warning not for very low values of GPODTOT
        CHKGPOD = NOTNUL(GPODTOT)/NOTNUL(GPD)
        IF (GPODTOT.GT.10 .AND. 
     &    ((CHKGPOD.LT.0.99) .OR. (CHKGPOD .GT.1.01))) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2 - POD',
     &    'WA7 - Sum pod class growth not equal to total pod growth')
          TERMNL = .TRUE.
        ENDIF

*       Initialise the first boxcar
        IF (GPD.GT.0.) THEN
          BH(0)    = BHRA*TMAV+BHRB              ![-]
        ELSE
          BH(0)    = 0.                          ![-]
        ENDIF

*       Calculate yield
        YLDPD = 0.
        DO 40 I1 = IPOD,0,-1
          IF (STAGE(I1).GE.1.) THEN
            YLDPD = YLDPD + WPOD(I1)/DELT +GPOD(I1) ![kg DW ha-1 d-1]
            BHYLD = BH(I1)                       ![-]
            WPOD(I1)  = 0.                       ![kg DW ha-1]
            STAGE(I1) = 0.                       ![-]
            SSPOD(I1) = 0.                       ![-]
            BH(I1)    = 0.                       ![-]
            IPOD      = I1 - 1                   ![-]
          ENDIF
40      CONTINUE

      ELSE IF (ITASK.EQ.3) THEN
*       ------------------------
*       Integration section
*       ------------------------

*       Add DW growth, increase stage and move pod classes one boxcar
        TSS     = 0.                             ![-]
        WPODTOT = 0.                             ![kg DW ha-1]
        DO 50 I1 = IPOD,0,-1
          WPOD(I1+1)  = WPOD(I1) + GPOD(I1) * DELT  ![kg DW ha-1]
          STAGE(I1+1) = STAGE(I1)+ DEVR * DELT   ![-]    
          SSPOD(I1+1) = LINT(SSTB,ISSN,STAGE(I1+1)) ![-]
          TSS         = TSS + SSPOD(I1+1)        ![-]
          BH(I1+1)    = BH(I1)                   ![-]
          WPODTOT     = WPODTOT + WPOD(I1+1)     ![kg DW ha-1]
50      CONTINUE
        IPOD = IPOD + 1

*       Give warning if WPODTOT differs from WPD
*       Warning not for very low values of WPODTOT
        CHKWPOD = NOTNUL(WPODTOT)/NOTNUL(WPD+GPD*DELT-YLDPD*DELT) ![-]
        IF (WPODTOT.GT.10 .AND. 
     &	((CHKWPOD.LT.0.99) .OR. (CHKWPOD .GT.1.01))) THEN
          WRITE (*,*) IYEAR,IDOY
          CALL WARNING ('CASE2 - POD',
     &    'WA8 - Sum pod class weight not equal to total pod weight')
          TERMNL = .TRUE.
        ENDIF

      ENDIF
      
      RETURN
      END
