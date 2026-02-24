        
         program atm_chem
 
C     This program was created by James Kasting and modified mainly by
c     Alex Pavlov (AP) and Kara Kerelove (KK). The user-friendly version 
c     of this code was created by Antigona Segura (AS) (2005).
c     Some modifications are identified with the author's initials.
c
c     The code is mostly written in f77 but is compiled in f90 and it 
c     contains some f90 features.
c
C     This program was created from PRESO3.F in June, 1995, for the
C     purpose of calculating ozone levels on planets around stars of
C     different spectral types. It differs from PRESO3 by including 
C     NO3, N2O5, and CL2O2, and by not including the anthropogenic
C     CFC's. Sulfur chemistry was added by Alex Pavlov.
c     
c     The UV fluxes for stars that are not the Sun were provided by
c     Martin Cohen from new IUE (+model) data. 
c
c     Check the notes along the main program before use it. Look for 
c     the word 'NOTE:'
C
C     This program has been updated to the JPL '92 rate constant
C     recommendations. I have not checked all the photolysis cross
C     sections to see if they are up to date.
C
C         THIS PROGRAM IS DESIGNED SPECIFICALLY FOR AN O2 LEVEL
C     OF 1 PAL.  IT DIFFERS FROM THE OLD PRIMO3 BY INCLUDING AN
C     UPDATED ALGORITHM FOR CALCULATING PHOTOLYSIS OF O2 IN THE
C     SCHUMANN-RUNGE BANDS.  THE PHOTOLYSIS ALGORITHMS FOR O2 AND
C     NO MAY NOT BE APPLICABLE TO LOWER O2 LEVELS.  (THIS HAS NOT
C     BEEN CHECKED, BUT IT WOULD BE SURPRISING IF THEY WERE.)
C
C     THIS VERSION OF PRESO3 ALSO CONTAINS THE MODIFICATIONS USED IN 
C     PRIMO2 TO CALCULATE O2 SCHUMANN-RUNGE PHOTOLYSIS BY EXPONENTIAL
C     SUMS AND NO PHOTOLYSIS BY THE MODIFIED CIESLIK AND NICOLET METHOD.
C     THIS EDITION ALSO HAS GIORGI AND CHAMEIDES RAINOUT RATES AND A 
C     MANABE-WETHERALD RELATIVE HUMIDITY DISTRIBUTION. 
C
C       THIS PROGRAM IS A ONE-DIMENSIONAL MODEL OF THE PRIMORDIAL
C     ATMOSPHERE.  THE MIXING RATIOS OF THE LONG-LIVED SPECIES
C     ARE CALCULATED FROM THE EQUATION
C
C     DF/DT  =  (1/N)*D/DZ(KN*DF/DZ) + P/N - LF
C
C     WHERE
C     F = MIXING RATIO (USOL)
C     K = EDDY DIFFUSION COEFFICIENT (EDD)
C     N = TOTAL NUMBER DENSITY (DEN)
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTONS
C     METHOD, E.G. IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.
C
C          THE LIST OF SUBROUTINES IS AS FOLLOWS:
C     (1) GRID   -  SETS UP THE ALTITUDE GRID
C     (2) RATES  -  DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE
C     (3) PHOTO  -  COMPUTES PHOTOLYSIS RATES (CALLS TWOSTR)
C     (4) DENSTY -  COMPUTES ATMOSPHERIC DENSITIES FROM HYDROSTATIC
C                   EQUILIBRIUM AND INITIALIZES ABSORBER PROFILES
C     (5) DIFCO  -  COMPUTES DK = K*N BETWEEN GRID POINTS
C     (6) OUTPUTP -  PRINTS OUT RESULTS
C     (7) DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C                  POINTS BY CALLING CHEMPL
C     (8) CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C                  FOR ONE SPECIES AT ALL GRID POINTS
C     (9) LTNING -  COMPUTES LIGHTNING PRODUCTION RATES FOR O2 AND
C                   N2 BASED ON CHAMEIDES' RESULTS
C    (10) SCALE  -  CAN BE USED TO SCALE THE EQUATIONS BY ARBITRARY
C                   FACTORS TO IMPROVE CONDITIONING OF MATRIX
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (1) TBDY   -  COMPUTES 3-BODY REACTION RATES
C     (2) E1     - EXPONENTIAL INTEGRAL OF ORDER ONE
C

C ***** REACTION LIST *****
C     1)  H2O + O(1D) = 2OH
C     2)  H2 + O(1D) = OH + H
C     3)  H2 + O = OH + H
C     4)  H2 + OH = H2O + H
C     5)  H + O3 = OH + O2
C     6)  H + O2 + M = HO2 + M
C     7)  H + HO2 = H2 + O2
C     8)  H + HO2 = H2O + O
C     9)  H + HO2 = OH + OH
C    10)  OH + O = H + O2
C    11)  OH + HO2 = H2O + O2
C    12)  OH + O3 = HO2 + O2
C    13)  HO2 + O = OH + O2
C    14)  HO2 + O3 = OH + 2O2
C    15)  HO2 + HO2 = H2O2 + O2
C    16)  H2O2 + OH = HO2 + H2O
C    17)  O + O + M = O2 + M
C    18)  O + O2 + M = O3 + M
C    19)  O + O3 = 2O2
C    20)  OH + OH = H2O + O
C    21)  O(1D) + N2 = O(3P) + N2
C    22)  O(1D) + O2 = O(3P) + O2
C    23)  O2 + HV = O(3P) + O(1D)
C    24)  O2 + HV = O(3P) + O(3P)
C    25)  H2O + HV = H + OH
C    26)  O3 + HV = O2 + O(1D)
C    27)  O3 + HV = O2 + O(3P)
C    28)  H2O2 + HV = 2OH
C    29)  CO2 + HV = CO + O(3P)
C    30)  CO + OH = CO2 + H
C    31)  CO + O + M = CO2 + M
C    32)  H + CO + M = HCO + M
C    33)  H + HCO = H2 + CO
C    34)  HCO + HCO = H2CO + CO
C    35)  OH + HCO = H2O + CO
C    36)  O + HCO = H + CO2
C    37)  O + HCO = OH + CO
C    38)  H2CO + HV = H2 + CO
C    39)  H2CO + HV = HCO + H
C    40)  HCO + HV = H + CO
C    41)  H2CO + H = H2 + HCO
C    42)  CO2 + HV = CO + O(1D)
C    43)  H + H + M = H2 + M
C    44)  HCO + O2 = HO2 + CO
C    45)  H2CO + OH = H2O + HCO
C    46)  H + OH + M = H2O + M
C    47)  OH + OH + M = H2O2 + M
C    48)  H2CO + O = HCO + OH
C    49)  H2O2 + O = OH + HO2
C    50)  HO2 + HV = OH + O
C    51)  CH4 + HV  =  1CH2 + H2
C    52)  CH3OOH + HV  =  H3CO + OH
C    53)  N2O + HV  =  N2 + O
C    54)  HNO2 + HV  = NO + OH
C    55)  HNO3 + HV  = NO2 + OH
C    56)  NO + HV  =  N + O
C    57)  NO2 + HV  =  NO + O
C    58)  CH4 + OH  =  CH3 + H2O
C    59)  CH4 + O(1D)  =  CH3 + OH
C    60)  CH4 + O(1D)  =  H2CO + H2
C    61)  1CH2 + CH4  =  2 CH3
C    62)  1CH2 + O2  =  H2CO + O
C    63)  1CH2 + N2  =  3CH2 + N2
C    64)  3CH2 + H2  =  CH3 + H
C    65)  3CH2 + CH4  =  2 CH3
C    66)  3CH2 + O2  =  H2CO + O
C    67)  CH3 + O2 + M  =  CH3O2 + M
C    68)  CH3 + OH  =  H2CO + H2
C    69)  CH3 + O  =  H2CO + H
C    70)  CH3 + O3  =  H2CO + HO2
C    71)  CH3O2 + HO2  =  CH3OOH + O2
C    72)  CH3O2 + CH3O2  =  2 H3CO + O2
C    73)  CH3O2 + NO  =  H3CO + NO2
C    74)  H3CO + O2  =  H2CO + HO2
C    75)  H3CO + O  =  H2CO + OH
C    76)  H3CO + OH  =  H2CO + H2O
C    77)  N2O + O(1D)  =  NO + NO
C    78)  N2O + O(1D)  =  N2 + O2
C    79)  N + O2  =  NO + O
C    80)  N + O3  =  NO + O2
C    81)  N + OH  =  NO + H
C    82)  N + NO  =  N2 + O
C    83)  NO + O3  =  NO2 + O2
C    84)  NO + O + M  =  NO2 + M
C    85)  NO + HO2  =  NO2 + OH
C    86)  NO + OH + M  =  HNO2 + M
C    87)  NO2 + O  =  NO + O2
C    88)  NO2 + OH + M  =  HNO3 + M
C    89)  NO2 + H  =  NO + OH
C    90)  HNO3 + OH  =  H2O + NO3
C    91)  HO2 + NO2 + M  =  HO2NO2 + M
C    92)  HO2NO2 + OH  =  NO2 + H2O + O2
C    93)  HO2NO2 + O  =  NO2 + OH + O2
C    94)  HO2NO2 + M  =  HO2 + NO2 + M
C    95)  HO2NO2 + HV  =  HO2 + NO2
C    96)  CH3OOH + OH  =  CH3O2 + H2O
C    97)  CH3O2 + OH  = H3CO + HO2
C    98)  O3 + NO2  =  O2 + NO3
C    99)  NO2 + NO3  =  NO + NO2 + O2
C   100)  O + NO3  =  O2 + NO2
C   101)  NO + NO3  =  NO2 + NO2
C   102)  OH + NO3  =  HO2 + NO2
C   103)  HO2 + NO3  =  HNO3 + O2
C   104)  NO2 + O + M  =  NO3 + M
C   105)  NO3 + hv  =  NO2 + O
C   106)  NO3 + NO2 + M  =  N2O5 + M
C   107)  N2O5 + hv  =  NO2 + NO3
C   108)  N2O5 + M  =  NO2 + NO3 + M
C   109)  N2O5 + H2O  =  2 HNO3
C************** Added sulfur chemistry ********
C   110)  SO   + HV   =     S    +     O
C   111)  SO2  + HV   =     SO   +     O
C   112)  H2S  + HV   =     HS   +     H
C   113)  SO   + O2   =     O    +     SO2
C   114)  SO   + HO2  =     SO2  +     OH
C   115)  SO   + O    =     SO2
C   116)  SO   + OH   =     SO2  +     H
C   117)  SO2  + OH   =     HSO3
C   118)  SO2  + O    =     SO3
C   119)  SO3  + H2O  =     H2SO4
C   120)  HSO3 + O2   =     HO2  +     SO3
C   121)  HSO3 + OH   =     H2O  +     SO3
C   122)  HSO3 + H    =     H2   +     SO3
C   123)  HSO3 + O    =     OH   +     SO3
C   124)  H2S  + OH   =     H2O  +     HS
C   125)  H2S  + H    =     H2   +     HS
C   126)  H2S  + O    =     OH   +     HS
C   127)  HS   + O    =     H    +     SO
C   128)  HS   + O2   =     OH   +     SO
C   129)  HS   + HO2  =     H2S  +     O2
C   130)  HS   + HS   =     H2S  +     S
C   131)  HS   + HCO  =     H2S  +     CO
C   132)  HS   + H    =     H2   +     S
C   133)  HS   + S    =     H    +     S2
C   134)  S    + O2   =     SO   +     O
C   135)  S    + OH   =     SO   +     H
c-as This reaction was deleted and subtituted by 182a).Note: 182a is from the old order
C   136)  S    + HCO  =     HS   +     CO
c   136a) SO2  + HV   =     S    +     O2
C   137)  S    + HO2  =     HS   +     O2
C   138)  S    + HO2  =     SO   +     OH
C   139)  HS   + H2CO =     H2S  +     HCO
C   140)  SO2  +     HV  =      SO21
C   141)  SO2  +     HV  =      SO23
C   142)  H2SO4 +     HV =       SO2   +   OH  +    OH
C   143)  SO3   +    HV  =      SO2  +     O
C   144)  SO21  +    M   =      SO23 +     M
C   145)  SO21  +    M   =      SO2  +     M
C   146)  SO21  +    HV  =      SO23 +     HV
C   147)  SO21  +    HV  =      SO2  +     HV
C   148)  SO21  +    O2  =      SO3  +     O
C   149)  SO21  +    SO2 =      SO3  +     SO
C   150)  SO23  +    M   =      SO2  +     M
C   151)  SO23  +    HV  =      SO2  +     HV
C   152)  SO23  +    SO2 =      SO3  +     SO
C   153)  SO    +    NO2 =      SO2  +     NO
C   154)  SO    +    O3  =      SO2  +     O2
C   155)  SO2   +    HO2 =      SO3  +     OH
C   156)  HS    +    O3  =      HSO  +     O2
C   157)  HS    +    NO2 =      HSO  +     NO
C   158)  S     +    O3  =      SO   +     O2
C   159)  SO    +    SO  =      SO2  +     S
C   160)  SO3   +    SO  =      SO2  +     SO2
C   161)  S     +    CO2 =      SO   +     CO
C   162)  SO    +    HO2 =      HSO  +     O2
C   163)  SO    +    HCO =      HSO  +     CO
C   164)  H     +    SO  =      HSO
C   165)  HSO   +    HV  =      HS   +     O
C   166)  HSO   +    OH  =      H2O  +     SO
C   167)  HSO   +    H   =      HS   +     OH
C   168)  HSO   +    H   =      H2   +     SO
C   169)  HSO   +    HS  =      H2S  +     SO
C   170)  HSO   +    O   =      OH   +     SO
C   171)  HSO   +    S   =      HS   +     SO
c-as Reactions added for N2O 
c   172)  N2 + O1D = N2O
c   173)  N2O + H = NO + NO + OH
c   174)  N2O + NO = NO2 + N2
c   175)  O1D + CO2 = CO2 + O
C   176)  CO + O1D = CO2
C***********************************************************
C
C***********************************************************
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.  THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS IN FIVE 10-DIGIT COLUMNS
C     STARTING IN COLUMN 11, I.E.
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NSP1 = NSP + 1 (INCLUDES HV)
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C               IS SOLVED
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C        PHOTOLYSIS REACTIONS ARE IDENTIFIED BY THE SYMBOL HV (NOT
C     COUNTED IN EVALUATING NSP).  THREE-BODY REACTIONS ARE WRITTEN
C     IN TWO-BODY FORM, SO THE DENSITY FACTOR MUST BE INCLUDED IN
C     THE RATE CONSTANT.
C        THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICE
C
C     ISPEC(NSP2) = VECTOR CONTAINING THE HOLLERITH NAMES OF THE
C                  CHEMICAL SPECIES.  THE LAST ENTRY MUST BE HV.
C     JCHEM(5,NR) = MATRIX OF CHEMICAL REACTIONS.  THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     ILOSS(2,NSP,NMAX) = MATRIX OF LOSS PROCESSES.  ILOSS(1,I,L)
C                         HOLDS REACTION NUMBER J, ILOSS(2,I,L) HOLDS
C                         REACTANT NUMBER.
C     IPROD(NSP,NMAX) = MATRIX OF PRODUCTION PROCESSES.  IPROD(I,L)
C                       HOLDS REACTION NUMBER J.
C     NUML(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF ILOSS
C     NUMP(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF IPROD
C
c      PARAMETER(NZ=64, NQ=34, NQT=NQ+1)
c      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
c      PARAMETER(NR=220, NF=34)
c      PARAMETER(NSP=55, NSP1=NSP+1, NSP2=NSP+2, NMAX=70)

       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNEQ_LDA.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
C-AJ 03/04/2022 ADD PARAMETER FOR GAUSSIAN INTEGRATION
       INCLUDE 'INCLUDECHEM/parNZA.inc'
c-AJ 12/21/2022 ADD CORRELATED K TABLE
       INCLUDE 'INCLUDECHEM/parCORRK.inc'
C      PARAMETER(NQI=31, NSLS=38, NMAXI=70)
C      PARAMETER(NRI=353, NSPI=78, NSPI1=NSPI+1, NSPI2=NSPI+2)

      PARAMETER(NQI=29, NSLS=36, NMAXI=70)
      PARAMETER(NRI=353, NSPI=76, NSPI1=NSPI+1, NSPI2=NSPI+2)
      
C-AJ 03/10/2022 ADD PARAMETER FOR GAUSSIAN QUARTRATURE
      PARAMETER(nrow=11)      
!C-PL AS I ADDED SGFLUX(NQ),SMFLUX(NQ),VDEP(NQ) TO comBBLOK.inc, delete their statement as follow
      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
     2  ,USAVE(NQ,NZ),R(NZ),U(NQ)!SGFLUX(NQ),SMFLUX(NQ),VDEP(NQ),
      DIMENSION UOLD(29,100),TOLD(100),EDDOLD(100),DENOLD(100),
     2  SO4AEROLD(100),AERSOLOLD(100,3),WFALLOLD(100,3),RPAROLD(100,3)
      DIMENSION DPU(NZ,3),DPL(NZ,3)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)                                    
      DIMENSION TSAV(NZ),TDUM(NZ),TP1(NZ)
      DIMENSION water(NZ),FLOW(NQT),fluxsave(108),sfxsave(10)
C-AJ 08/31/2022 ADD UOLD2 TO READ FROM THE ATMCHEM/IO/atm_composition.dat
      DIMENSION UOLD2(34,100)

      CHARACTER :: STARR*3,DIRDATA*4, AA*11,DIRIO*2
C-AJ 01/24/2022 ADD VARIABLES FOR CH4 LIFETIME CALCULATION
      REAL:: lsCH4_OH,lsCH4_D,lossCH4,lifetimeCH4,lifetimeCH4_1
C-AJ 03/04/2022 ADD PARAMETER FOR GAUSSIAN INTEGRATION
      REAL,DIMENSION(nrow,20)::xi,wi
      INTEGER,DIMENSION(nrow)::NumGau
C-AJ 05/06/2022 ADD OH SOURCE PRINTOUT
      REAL,DIMENSION(NZ):: RateOH1,RateOH25,RateOH59,RatePO3   

C      REAL,DIMENSION(6)::ZYVEC=(/-21.177,-48.608,-76.195,21.177,48.608,
C     2  76.195/)
C      REAL,DIMENSION(7)::ZYVEC2=(/90,66.056,42.138,18.358,-66.056,
C     2  -42.138,-18.358/)
C      REAL,DIMENSION(8)::ZYVEC3=(/79.430,58.296,37.187,16.200,-79.430,
C     2  -58.296,-37.187,-16.200/)
C      REAL,DIMENSION(9)::ZYVEC4=(/90,71.080,52.166,33.277,14.497,
C     2  -71.080,-52.166,-33.277,-14.497/)
C      REAL,DIMENSION(10)::ZYVEC5=(/81.438,64.317,47.202,30.110,13.118,
C     2  -81.438,-64.317,-47.202,-30.110,-13.118/)
C      REAL,DIMENSION(11)::ZYVEC6=(/90,74.363,58.728,43.101,27.494,
C     2  11.978,-74.363,-58.728,-43.101,-27.494,-11.978/)
C      REAL,DIMENSION(12)::ZYVEC7=(/82.806,68.418,54.033,39.655,25.295,
C     2  11.020,-82.806,-68.418,-54.033,-39.655,-25.295,-11.020/)

C      REAL,DIMENSION(6)::WEIGHT=(/0.17132,0.36076,0.46791,0.17132,
C     2  0.36076,0.46791/)
C      REAL,DIMENSION(7)::WEIGHT2=(/0.41796,0.38183,0.27970,0.12948,
C     2  0.38183,0.27970,0.12948/)
C      REAL,DIMENSION(8)::WEIGHT3=(/0.36268,0.31371,0.22238,0.10123,
C     2  0.36268,0.31371,0.22238,0.10123/)
C      REAL,DIMENSION(9)::WEIGHT4=(/0.33024,0.31235,0.26061,0.18065,
C     2  0.08127,0.31235,0.26061,0.18065,0.08127/)
C      REAL,DIMENSION(10)::WEIGHT5=(/0.29552,0.26927,0.21909,0.14945,
C     2  0.29552,0.26927,0.21909,0.14945/)
C      REAL,DIMENSION(11)::WEIGHT6=(/0.27293,0.26280,0.23319,0.18629,
C     2  0.12558,0.05567,0.26280,0.23319,0.18629,0.12558,0.05567/)
C      REAL,DIMENSION(12)::WEIGHT7=(/0.24915,0.23349,0.20317,0.16008,
C     2  0.10694,0.04718,0.24915,0.23349,0.20317,0.16008,0.10694,
C     3  0.04718/)

C-PL ADDED NEW VARIABLES FOR ISOTOPE CALCULATION
C AJ & JL ADD HQO AND HQONO2 
C      DIMENSION DELTAH2CO17(NZ),DELTAO17(NZ),DELTAH2O17(NZ),
C     2 DELTAOH17(NZ),DELTAHO2A17(NZ),DELTAH2O217(NZ),DELTAO3A17(NZ),
C     3 DELTAO3B17(NZ),DELTACO17(NZ),DELTACH3OOHA17(NZ),
C     4 DELTACH3OOHB17(NZ),DELTACH3O2A17(NZ),DELTACH3O2B17(NZ),
C     5 DELTAN2O17(NZ),DELTANO17(NZ),DELTANO217(NZ),DELTAHNO217(NZ),
C     6 DELTAHNO317(NZ),DELTAHO2NO2A17(NZ),DELTAHO2NO2B17(NZ),
C     7 DELTANO317(NZ),DELTAN2O5A17(NZ),DELTAN2O5B17(NZ),DELTAO217(NZ),
C     8 DELTASO17(NZ),DELTASO217(NZ),DELTAH2SO417(NZ),DELTAHSO17(NZ),
C     9 DELTACO217(NZ),DELTAO317(NZ),DELTACH3OOH17(NZ),
C     1 DELTACH3O217(NZ),DELTAHO2NO217(NZ),DELTAN2O517(NZ),
C     2 DELTAHO2B17(NZ),DELTAHO2NO2C17(NZ),DELTAHO217(NZ) 

C      DIMENSION C2DELTAH2CO(NZ),C2DELTAO(NZ),C2DELTAH2O(NZ),
C     2 C2DELTAOH(NZ),C2DELTAHO2A(NZ),C2DELTAH2O2(NZ),C2DELTAO3A(NZ),
C     3 C2DELTAO3B(NZ),C2DELTACO(NZ),C2DELTACH3OOHA(NZ),
C     4 C2DELTACH3OOHB(NZ),C2DELTACH3O2A(NZ),C2DELTACH3O2B(NZ),
C     5 C2DELTAN2O(NZ),C2DELTANO(NZ),C2DELTANO2(NZ),C2DELTAHNO2(NZ),
C     6 C2DELTAHNO3(NZ),C2DELTAHO2NO2A(NZ),C2DELTAHO2NO2B(NZ),
C     7 C2DELTANO3(NZ),C2DELTAN2O5A(NZ),C2DELTAN2O5B(NZ),C2DELTAO2(NZ),
C     8 C2DELTASO(NZ),C2DELTASO2(NZ),C2DELTAH2SO4(NZ),
C     9 C2DELTAHSO(NZ),C2DELTACO2(NZ),C2DELTAO3(NZ),C2DELTACH3OOH(NZ),
C     1 C2DELTACH3O2(NZ),C2DELTAHO2NO2(NZ),C2DELTAN2O5(NZ),
C     2 C2DELTAHO2B(NZ),C2DELTAHO2NO2C(NZ),C2DELTAHO2(NZ) 

      DIMENSION DELTAH2CO17(NZ),DELTAO17(NZ),DELTAH2O17(NZ),
     2 DELTAOH17(NZ),DELTAHO217(NZ),DELTAH2O217(NZ),DELTAO3A17(NZ),
     3 DELTAO3B17(NZ),DELTACO17(NZ),DELTACH3OOHA17(NZ),
     4 DELTACH3OOHB17(NZ),DELTACH3O2A17(NZ),DELTACH3O2B17(NZ),
     5 DELTAN2O17(NZ),DELTANO17(NZ),DELTANO217(NZ),DELTAHNO217(NZ),
     6 DELTAHNO317(NZ),DELTAHO2NO2A17(NZ),DELTAHO2NO2B17(NZ),
     7 DELTANO317(NZ),DELTAN2O5A17(NZ),DELTAN2O5B17(NZ),DELTAO217(NZ),
     8 DELTASO17(NZ),DELTASO217(NZ),DELTAH2SO417(NZ),DELTAHSO17(NZ),
     9 DELTACO217(NZ),DELTAO317(NZ),DELTACH3OOH17(NZ),
     1 DELTACH3O217(NZ),DELTAHO2NO217(NZ),DELTAN2O517(NZ) 

      DIMENSION C2DELTAH2CO(NZ),C2DELTAO(NZ),C2DELTAH2O(NZ),
     2 C2DELTAOH(NZ),C2DELTAHO2(NZ),C2DELTAH2O2(NZ),C2DELTAO3A(NZ),
     3 C2DELTAO3B(NZ),C2DELTACO(NZ),C2DELTACH3OOHA(NZ),
     4 C2DELTACH3OOHB(NZ),C2DELTACH3O2A(NZ),C2DELTACH3O2B(NZ),
     5 C2DELTAN2O(NZ),C2DELTANO(NZ),C2DELTANO2(NZ),C2DELTAHNO2(NZ),
     6 C2DELTAHNO3(NZ),C2DELTAHO2NO2A(NZ),C2DELTAHO2NO2B(NZ),
     7 C2DELTANO3(NZ),C2DELTAN2O5A(NZ),C2DELTAN2O5B(NZ),C2DELTAO2(NZ),
     8 C2DELTASO(NZ),C2DELTASO2(NZ),C2DELTAH2SO4(NZ),
     9 C2DELTAHSO(NZ),C2DELTACO2(NZ),C2DELTAO3(NZ),C2DELTACH3OOH(NZ),
     1 C2DELTACH3O2(NZ),C2DELTAHO2NO2(NZ),C2DELTAN2O5(NZ)

      CHARACTER*5 label1
      CHARACTER*40 CHEMI(5,NRI)
      DIMENSION DELTA17O(NQI),CAP17(NQI),USOLI(NQI,NZ),ZKM(NZ)
      REAL DELTAO3,DELTACH3OOH,DELTACH3O2,DELTAHO2NO2,DELTAN2O5,
     2 CAPO317, CAPCH3OOH17,CAPC3O217,CAPHO2NO217,CAPN2O517
     3 ,CAPHO217,LIGHTNO2I
      COMMON/ISOTOP/ARI(NRI,NZ),ILOSSI(2,NSPI,NMAXI),
     2  JCHEMI(5,NRI),NUMLI(NSPI),NUMPI(NSPI),TPI(NQI),TLI(NQI),
     3  YPI(NQI,NZ),YLI(NQI,NZ),SRI(NQI),TLOSSI(NQI),PHIDEPI(NQI),
     4 ISPECI(NSPI2),DIZ(NSPI2,NZ),IPRODI(NSPI,NMAXI),PSO4AER
     5 ,LBOUNDI(NQI),FLOWI(NQI),FUPI(NQI),CONI(NQI)

      COMMON/RATESIS/RAINGCI(NQI,NZ),DDI(NQI,NZ),DKI(NQI,NZ),
     2  DUI(NQI,NZ),DLI(NQI,NZ),DHUI(NQI,NZ),DHLI(NQI,NZ),HIZ(NQI,NZ),
     3  VDEPI(NQI),RATI(NRI)

      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ
     4  LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LLQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2
C JL 5TH DAY OF THE LUNAR NEW YEAR 2020 
C      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
C     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
C     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,
C     4  LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
C     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
C     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
C     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
C     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
C     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2
C-PL END ADD NEW VARIABLES
     
c Name of the star
      INCLUDE 'INCLUDECHEM/comSTR.inc'
      INCLUDE 'INCLUDECHEM/comFLUXPHOTO.inc'
c DIRP contains DIRCOUP, this variable is defined as a character 
c in comDIRP.inc
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comEBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comGBLOK.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comPRESS1.inc'
      INCLUDE 'INCLUDECHEM/comQBLOK.inc'
      INCLUDE 'INCLUDECHEM/comRBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'
C-AJ 03/04/2022 ADD PARAMETER FOR GAUSSIAN INTEGRATION
      INCLUDE 'INCLUDECHEM/comSBLOKnew.inc'
c-AJ 12/21/2022 ADD new correlated k table
      INCLUDE 'INCLUDECHEM/comCORRK.inc'
C-AJ 01/04/2022 change ALPHA(17,4) TO ALPHAOLD(17,4) in comQBLOK.inc
C ALPHA contains the weights for 12-pt Correlated k table
       data ALPHA/0.023587668193256, 0.053469662997659,
     & 0.080039164271673, 0.101583713361533,
     & 0.116746268269177, 0.124573522906701,
     & 0.124573522906701, 0.116746268269177, 
     & 0.101583713361533, 0.080039164271673,
     & 0.053469662997659, 0.023587668193256/
c in May 6, 2019, JL rip off species of Cl and S. 
c-AJ 08/30/2022 ADD CL BACK
        DATA LH2CO,LO,LH2O,LOH,LHO2,LH2O2,LO3,LH,LH2,LCH4,LCO,
     2  LCH3OOH,LCH3O2,LN2O,LNO,LNO2,LHNO2,LHNO3,LHO2NO2,LNO3,LN2O5,
     3  LCL2O2,LCH3CL,LHOCL,LCL,LCLO,LHCL,LCLONO2,LO2,LH2S,LHS,LSO,
     4  LSO2,LH2SO4,LHSO,LCO2,LSO4AER,LCH21,LCH23,LO1D,LCH3,LH3CO,
     5  LHCO,LN,LNOCL,LCLONO,LCLO2,LCL2,LS,LSO21,LSO23,LHSO3,LSO3,
     6  LS2,LN2/
     7  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     8  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     9  44,45,46,47,48,49,50,51,52,53,54,55/

c in may 17 2019, JL move CO2 from inert to long lived 
C        DATA LH2CO,LO,LH2O,LOH,LHO2,LH2O2,LO3,LH,LH2,LCH4,LCO,
C     2  LCH3OOH,LCH3O2,LN2O,LNO,LNO2,LHNO2,LHNO3,LHO2NO2,LNO3,LN2O5,
C     3  LO2,LH2S,LHS,LSO,LSO2,LH2SO4,LHSO,LCO2,LSO4AER,LCH21,LCH23,
C     4  LO1D,LCH3,LH3CO,LHCO,LN,LS,LSO21,LSO23,LHSO3,LSO3,LS2,LN2/
C     5  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
C     6  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
C     7  44/

C-AP
C
C   NO PREDISSOCIATION COEFFICIENTS (ALLEN AND FREDERICK, 1982)
C
      DATA ANO/-1.790868E+1, -1.924701E-1, -7.217717E-2, 5.648282E-2,
     2  4.569175E-2, 8.353572E-3, 3*0.,
     3  -1.654245E+1, 5.836899E-1, 3.449436E-1, 1.700653E-1,
     4  -3.324717E-2, -4.952424E-2, 1.579306E-2, 1.835462E-2,
     5  3.368125E-3/
C
      DATA BNO/7.836832E+3, -1.549880E+3, 1.148342E+2, -3.777754E+0,
     2  4.655696E-2, 1.297581E+4, -2.582981E+3, 1.927709E+2,
     3  -6.393008E+0, 7.949835E-2/
C
      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/
      DATA RNO2/60*0., .79, .83, .66, .15, 4*0./
C
C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/
C
C ***** SOLUBILITY (GIORGI AND CHAMEIDES) *****
C      DATA H/1.3E+04, 1.0E-99, 1.0E-99, 1.0E+05, 3.3E+04,
C     2       2.0E+05, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99,
C     3       1.0E-03, 2.0E+05, 3.3E+04, 2.5E-02, 1.9E-03,
C     4       7.0E-03, 7.0E+11, 7.0E+11, 7.0E+11, 7.0E-03,
C     5       7.0E+11, 3.2E-4,   0.14,    1.E+5,   1.9E-3,   
c                     o2       
C     6       1.E+4, 7E+11, 9E+3,3.72E-2/
c                                co2 
c-AJ 08/30/2022 ADD CL BACK
C ***** SOLUBILITY (GIORGI AND CHAMEIDES) *****
      DATA H/1.3E+04, 1.0E-99, 1.0E-99, 1.0E+05, 3.3E+04,
     2       2.0E+05, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99,
     3       1.0E-03, 2.0E+05, 3.3E+04, 2.5E-02, 1.9E-03,
     4       7.0E-03, 7.0E+11, 7.0E+11, 7.0E+11, 7.0E-03,
     5       7.0E+11, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99,
     6       1.0E-99, 7.0E+11, 1.0E-99, 3.2E-4,   0.14,      
c                                        o2       
     7       1.E+5, 1.9E-3,1.E+4, 7E+11, 9E+3,3.72E-2/
c                                              co2 
C-AP  I have changed the Henry constant similar to the Archean code
C-AP  only for the sulfur species  from Archean. Note that Henry(SO2) 
C-AP  in the table is higher because we use effective Henry constant
C-AP  from Archean.

c  Temperature from the US Standard Atmosphere 1976. Used when the 
c  code is not coupled to the climate model.
C  JK   Data are estimated above 64 km. I'm adjusting the temperature
C       near the tropopause downward in order to get better 
C       statospheric H2O
      DATA T/288.15, 281.65, 275.15, 268.66, 262.17,
     &       255.68, 249.19, 242.70, 236.21, 229.73,
     &       223.25, 216.77, 216.65, 216.65, 216.65,
     &       216.65, 216.65, 216.65, 216.65, 216.65,
     &       216.65, 217.58, 218.57, 219.57, 220.56, 
     &       221.55, 222.54, 223.54, 224.53, 225.52,
     &       226.51, 227.50, 228.49, 230.97, 233.74, 
     &       236.51, 239.28, 242.05, 244.82, 247.58,
     &       250.35, 253.14, 255.88, 258.64, 261.40,
     &       264.16, 266.96, 269.68, 270.65, 270.65,
     &       270.65, 270.65, 269.03, 266.27, 263.52, 
     &       260.77, 258.02, 255.27, 252.52, 249.77, 
     &       247.02, 244.27, 241.53, 230.78, 230.50,
     &       230.50, 227.00, 225.10, 222.00, 219.60,
     &       216.00, 214.30, 212.00, 210.30, 208.00,
     &       206.40, 204.00, 202.50, 200.00, 198.60,
     &       196.00, 194.70, 192.00, 190.80, 189.00,
     &       186.90, 186.90, 186.90, 186.90, 186.90,
     &       188.00, 190.00, 191.00, 192.50, 194.00,
     &       195.00, 196.50, 198.00, 200.00, 202.00/


      DATA TP1/288.15, 281.65, 275.15, 268.66, 262.17,
     &       255.68, 249.19, 242.70, 236.21, 229.73,
     &       223.25, 216.77, 216.65, 216.65, 216.65,
     &       216.65, 216.65, 216.65, 216.65, 216.65,
     &       216.65, 217.58, 218.57, 219.57, 220.56, 
     &       221.55, 222.54, 223.54, 224.53, 225.52,
     &       226.51, 227.50, 228.49, 230.97, 233.74, 
     &       236.51, 239.28, 242.05, 244.82, 247.58,
     &       250.35, 253.14, 255.88, 258.64, 261.40,
     &       264.16, 266.96, 269.68, 270.65, 270.65,
     &       270.65, 270.65, 269.03, 266.27, 263.52, 
     &       260.77, 258.02, 255.27, 252.52, 249.77, 
     &       247.02, 244.27, 241.53, 230.78, 230.50,
     &       230.50, 227.00, 225.10, 222.00, 219.60,
     &       216.00, 214.30, 212.00, 210.30, 208.00,
     &       206.40, 204.00, 202.50, 200.00, 198.60,
     &       196.00, 194.70, 192.00, 190.80, 189.00,
     &       186.90, 186.90, 186.90, 186.90, 186.90,
     &       188.00, 190.00, 191.00, 192.50, 194.00,
     &       195.00, 196.50, 198.00, 200.00, 202.00/
C
C JK Adjusting this downward in the lower troposphere to get Tsurf=288K
c      DATA T/285.15, 278.65, 272.15, 265.66, 261.17,
c     &       255.68, 249.19, 242.70, 236.21, 225.00,
c     &       210.00, 210.00, 210.00, 210.00, 210.00,
c     &       210.00, 214.00, 214.00, 216.65, 216.65,
c     &       216.65, 217.58, 218.57, 219.57, 220.56, 
c     &       221.55, 222.54, 223.54, 224.53, 225.52,
c     &       226.51, 227.50, 228.49, 230.97, 233.74, 
c     &       236.51, 239.28, 242.05, 244.82, 247.58,
c     &       250.35, 253.14, 255.88, 258.64, 261.40,
c     &       264.16, 266.96, 269.68, 270.65, 270.65,
c     &       270.65, 270.65, 269.03, 266.27, 263.52, 
c     &       260.77, 258.02, 255.27, 252.52, 249.77, 
c     &       247.02, 244.27, 241.53, 230.78, 230.50,
c     &       230.50, 227.00, 225.10, 222.00, 219.60,
c     &       216.00, 214.30, 212.00, 210.30, 208.00,
c     &       206.40, 204.00, 202.50, 200.00, 198.60,
c     &       196.00, 194.70, 192.00, 190.80, 189.00,
c     &       186.90, 186.90, 186.90, 186.90, 186.90,
c     &       188.00, 190.00, 191.00, 192.50, 194.00,
c     &       195.00, 196.50, 198.00, 200.00, 202.00/

C  Temperature profile for O2=0.1 PAL (from Segura et al.,2003, Fig.2)
      DATA TP2/285.15, 278.65, 272.15, 265.66, 261.17,
     &       255.68, 249.19, 242.70, 236.21, 225.00,
     &       210.00, 210.00, 210.00, 210.00, 210.00,
     &       210.00, 214.00, 214.00, 216.65, 216.65,
     &       216.65, 217.58, 218.57, 219.57, 220.56, 
     &       221.55, 222.54, 223.54, 224.53, 225.52,
     &       226.51, 227.50, 228.49, 230.97, 233.74, 
     &       236.51, 239.28, 237.00, 235.00, 233.00,
     &       230.00, 227.00, 224.00, 222.00, 220.00,
     &       218.00, 216.00, 214.00, 212.00, 210.00,
     &       209.00, 208.00, 207.00, 206.00, 205.00, 
     &       204.00, 203.00, 202.00, 201.00, 200.00, 
     &       200.00, 200.00, 200.00, 200.00, 200.00,
     &       200.00, 200.00, 200.00, 200.00, 200.00,
     &       200.00, 200.00, 200.00, 200.00, 200.00,
     &       200.00, 200.00, 200.00, 200.00, 198.60,
     &       196.00, 194.70, 192.00, 190.80, 189.00,
     &       186.90, 186.90, 186.90, 186.90, 186.90,
     &       188.00, 190.00, 191.00, 192.50, 194.00,
     &       195.00, 196.50, 198.00, 200.00, 202.00/


      DATA EDD/64*0., 5.174E+05, 5.674E+05, 6.224E+05, 6.826E+05,
     2  7.487E+05, 8.212E+05, 9.007E+05, 9.879E+05, 28*1.E6/
 
c NOTE: Boundary conditions are defined here
C ***** UPPER BOUNDARY FLUXES *****
      DATA SMFLUX/NQ*0./
C
C ***** EFFUSION VELOCITIES *****
      DATA VEFF/NQ*0./
C
C ***** UPPER BOUNDARY CONDITIONS *****
c      DATA MBOUND/0, 1, 8*0, 1, 24*0/
C JK 6/2019 Ignore downward fluxes of CO and O
c in may 6, JL removed CL               14
C      DATA MBOUND/0, 1, 8*0, 1, 3*0, 0,6*0, 0 ,6*0, 0  ,0/
C                    O       co     NO      O2     CO2  aer
C-AJ 08/30/2022 ADD CL BACK
      DATA MBOUND/0, 1, 8*0, 1, 26*0/
C   0 = CONSTANT EFFUSION VELOCITY (VEFF)
C   1 = CONSTANT FLUX (SMFLUX)
C
C ***** LOWER BOUNDARY FLUXES *****
c NOTE: SGFLUX is redefined later in the main code when the O2 
c mixing ratio is less than 0.21 or when the star is not the Sun
      DATA SGFLUX/NQ*0./
C
C ***** DEPOSITION VELOCITIES (NQ)*****
C      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1., 3*0., 0.2, 1.,
C     2  0., 2*0.2, 6*0.5, 0., 0.5, 1., 3*0.5, 0.02, 1, 3.E-4, 
C     3  3*1./

c in May 6 2018, JL rip off all the CL species 
C      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1.,3.e-4,2*0.,0.2,1.,
c               h2co o  h2o  oh  ho2 h2o2 o3  h        ch3ooh ch3o2
C     2  0., 2*0.2, 5*0.5, 1.e-4,0.02, 1., 3.E-4, 3*1.,0./
c      n2o                o2    h2s   hs   so         co2
C AJ 08/30/2022 ADD CL BACK
      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1.,3.45e-4,2*0.,0.2,1.,
C               h2co o  h2o  oh  ho2 h2o2 o3  h    H2        ch3ooh ch3o2
     2  0., 2*0.2, 6*0.5,0.,0.5,1.,3*0.5, 1.e-4,0.02, 1., 3.E-4, 3*1.,
c      n2o                                o2    h2s   hs   so         
     3  0./
C       co2

c Deposition velocities used for planets around quiescent M dwarfs
c      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1.,2.4e-4,0.,1.2e-4, 0.2,
c     2  1.,0., 2*0.2, 6*0.5, 0., 0.5, 1., 3*0.5, 0.02, 1, 3.E-4, 
c     3  3*1./
C
C ***** LOWER BOUNDARY CONDITIONS (NQT)*****
c NOTE: Mixing ratios for present Earth for H2,CO,CH4, N2O and CH3Cl
c are defined after the model parameters.
c fixed mixing ratios         H2 CH4 CO     N2O 
c      DATA LBOUND/2*0, 1, 5*0, 1, 1, 1, 2*0, 1, 5*0, 2*0, 1, 6*0,1,0/

C-PL FOR THE LOW O2 RUNS
c in May 6, 2019, JL rip off Cl.PENG commented temporarily for present run 
C       DATA LBOUND/2*0, 1, 5*0, 0, 2, 2, 2*0, 2, 5*0, 2*0, 1,6*0,1,0/
c                                             n2o         o2    co2

C-AJ ADD CL BACK
c fixed mixing ratios         H2 CH4 CO     N2O      ch3cl   o2    co2
      DATA LBOUND/2*0, 1, 5*0, 1, 1, 1, 2*0, 1, 8*0, 1, 5*0,1, 6*0, 1,0/
C-AJ 08/30/2022 ADD CL BACK
c fixed surface flux
c       DATA LBOUND/2*0, 1, 5*0, 0, 2, 2, 2*0, 2, 8*0,2,5*0, 1,6*0,1,0/
c

c Boundary conditions for planets around quiescent M dwarfs
c      DATA LBOUND/2*0, 1, 5*0, 0, 1, 0, 2*0, 2, 5*0, 3*0, 2, 12*0/
c Boundary conditions for AD Leo
c      DATA LBOUND/2*0, 1, 5*0, 0, 2, 0, 2*0, 2, 5*0, 3*0, 2, 12*0/

C
C   0 = CONSTANT DEPOSITION VELOCITY (VDEP)
C   1 = CONSTANT MIXING RATIO
C   2 = CONSTANT UPWARD FLUX (SGFLUX)
C
C Modify Henry constants for sensitivity run comparison to Archean atmosphere
      H(LHO2) = 9.E3
      H(LH2O2) = 6.2E5
      H(LH2CO) =  4.25E4

C============== FILE SECTION =============================

      DIRDATA = 'DATA'
      DIRIO = 'IO'

c **** INPUT FILES
c
C-AJ 12/21/2022 read correlated k table
C-AJ 12/21/2022 The data will be read in subroutine srbands.f under CHEM
      OPEN(UNIT=10, file =  DIRDATA//'/k12Gauss.dat', status='old')

c This file contains the chemical reactions used in the code
      OPEN(unit=61, file= DIRDATA//'/primo3s.chm')
c-as Next file contains solar flux and atmospheric data
c-as it MUST be read for all the cases
c-AJ 01/24/2022 replace H2O cross-sections with Ranjan et al., 2020 in photos.pdat  
      OPEN(unit=62, file= DIRDATA//'/photos.pdat', status='old') 
      OPEN(unit=63,file= DIRDATA//'/h2so4.pdat',status='old')
      OPEN(unit=64,file= DIRDATA//'/eddy.pdat',status='old')

c Files with the input parameters
c NOTE: Check this two before runnig the program
       OPEN(unit=65,file= DIRIO//'/input_atmchem.dat')
       OPEN(unit=66,file= DIRIO//'/planet.dat')  

c-as Next file was formerly named atm_chem_coefs.dat.
c-as I has the same format as atm_composition.out   
      OPEN(unit=67, file= DIRIO//'/atm_composition.dat') 

C-AJ 03/10/2022 GUASSIAN DATA
      OPEN(unit=68,file= DIRIO//'/Gaussian_factors.txt')
      
C-AJ 08/31/2022 ADD OLD atm_composition.dat from ATMCHEM
c      OPEN(UNIT=69,file= DIRIO//'/atm_composition2.dat')

c Next file is generated by the climate model contains altitude,
c temperature and water. Formerly called photo_input.dat.
c Only used when ICOUPLE=1
      OPEN(unit=71,file= DIRIO//'/fromClima2Photo.dat') 
      
c-as  Unit 72 is an input and output file that is shared with the climate code
c-as  when ICOUPLE= 1. It is OPEN later in the program to WRITE on it.
c Only used when ICOUPLE=1    
      OPEN(unit=72,file= DIRIO//'/mixing_ratios.dat')

      OPEN(unit=73, file= DIRDATA//'/faruvs.pdat')
      
c  NOTE: IMPORTANT files to read the far UV of all the stars (including 
c  the Sun) and the UV fluxes for stars others than the Sun.
c        74     fluxesKGF_photo.pdat
c        75     M star flux (name it as you like)
c        76     far UV flux (name depends on the star)

C **** OUTPUT FILES
c Main output
      OPEN(unit=90,file= DIRIO//'/outchem.dat')
c Files commented would contain information that is not
c needed for now but the write commands for them still in 
c program. To activate them just remove all c here and in 
c the write commands
c      OPEN(unit=10,file='primo3.plt')
c      OPEN(unit=14,file='primotemp.dat')
c      OPEN(unit=15,file='o3graph.dat')  

c-as MAIN output file MOVED to couple.f
C       OPEN(unit=90,file='outchem.dat')  

c Written on subroutine TWOSTR
c      OPEN(unit=82,file= DIRIO//'/wave_means.dat')

c This file contains the altitude in cm vs the ozone number density (cm^-3)
c      OPEN(unit=83,file= DIRIO//'/O3numdens.out')

c  This is an output file containing altitude, H2O, O3.
c  To be used as input of the climate model, formerly called Pass2SurfMP.dat
      OPEN(unit=84,file= DIRIO//'/fromPhoto2Clima.dat')

c This file contains total hydrogen mixing ratio vs. altitude (cm)
C These OUTPUT files are opened along the program 
      OPEN(unit=86,file= DIRIO//'/total_hydrogen.out')
      OPEN(unit=87,file= DIRIO//'/diffusion_coefs.out')

c In may 9 2019, Jl add the reaction rate table to the code to make life easier 
      OPEN(UNIT=15,file=DIRIO//'/int.rates.out.dat')
      OPEN(UNIT=16,file=DIRIO//'/int.iso_rates.out.dat')     
c in may 16 2019, JL add a clean printout for h2so4 to see possible disarray of input 
      open(UNIT=158, file=DIRIO//'/h2so4.out.dat') 

C AJ 05/08/2022 ADD OUTPUT FOR SOURCE OF OH
      open(UNIT=91,file=DIRIO//'/OH_rates.out')
C AJ 08/27/2022 ADD OUTPUT FOR PLOTS
      open(UNIT=92,file=DIRIO//'/OUTPUT_PLOT.dat')  
C These OUTPUT files are opened along the program 
c	UNIT   	NAME  
c        19     mixing_ratios.dat (main program)

c=================================================================

C********* SET MODEL PARAMETERS *****

C     ZY = SOLAR ZENITH ANGLE (IN DEGREES)
C     AGL = DIURNAL AVERAGING FACTOR FOR PHOTORATES
C     ISEASON = TELLS WHETHER P AND T VARY WITH TIME (THEY DON'T FOR
C               ISEASON < 3)
C     IZYO2 = TELLS WHETHER SOLAR ZENITH ANGLE VARIES WITH TIME (0 SAYS
C             IT DOESN'T; 1 SAYS IT DOES)
C     IO2 = 0 FOR ALLEN AND FREDERICK O2 SCHUMANN-RUNGE COEFFICIENTS
C         = 1 FOR EXPONENTIAL SUM FITS (FOR LOW-O2 ATMOSPHERES)
C     INO = 0 FOR ALLEN AND FREDERICK NO PREDISSOCIATION COEFFICIENTS
C         = 1 FOR MODIFIED CIESLIK AND NICOLET FORMULATION
C     EPSJ = AMOUNT BY WHICH TO PERTURB THINGS FOR JACOBIAN CALCULATION
C     ZTROP = TROPOPAUSE HEIGHT (ABOVE WHICH H2O BEHAVES AS A NONCONDENS
C             ABLE GAS
C     STARR - Character variable to choose a star, it can be:
c             Sun, F2V, K2V,dMV 
c             DO NOT FORGET quotation marks
c     ICOUPLE - 1 = Coupled to the climate model              
c               0 = Not coupled
C     FCO2 = CO2 mixing ratio when ICOUPLE = 0
C     FO2 = O2 mixing ratio when ICOUPLE = 0
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED
C     GauFLAG -0 = Use Gaussian quadrature
c              1 = single solar zenith angle
      read(65,555)
      read(65,*)AA,STARR
      read(65,*)AA,FLUXFAC
      read(65,*)AA,INIT
      read(65,*)AA,TSTOP
      read(65,*)AA,DT
      read(65,*)AA,NSTEPS
      read(65,*)AA,ZY
      read(65,*)AA,AGL
      read(65,*)AA,ISEASON
      read(65,*)AA,IZYO2
      read(65,*)AA,IO2
      read(65,*)AA,INO
      read(65,*)AA,EPSJ
      read(65,*)AA,ZTROP
      read(65,*)AA,FCO2
      read(65,*)AA,FO2   

      read(65,*)AA,GPPOXY
      read(65,*)AA,GPPCDE
      read(65,*)AA,COMBIN
      read(65,*)AA,RECOMB
      read(65,*)AA,GauFLAG
      read(65,*)AA,CLFLAG
  555  format(3/)
      close(65)
      ICOUPLE = 0

*** Read the flux from a star
      call readstar(FLUXFAC)

***** READ THE PLANET PARAMETER DATAFILE *****
      READ(66,502) G,FSCALE,ALB,DELZ,ZTROP,JTROP
 502  FORMAT(F5.1/,F4.2/,F5.3/,E5.1/,E5.1/,I2)
      close(66)
C
***** READ THE ATMOSPHERIC COMPOSITION
c      READ(67,400) UOLD,TOLD,EDDOLD,DENOLD,SO4AEROLD,AERSOLOLD,
c     2  WFALLOLD,RPAROLD
c      DO 3500 J=1,64
c      DO 3501 I=1,NQ
c 3501 USOL(I,J) = UOLD(I,J)
c      T(J) = TOLD(J)
c      EDD(J) = EDDOLD(J)
c      DEN(J) = DENOLD(J)
c      SO4AER(J) = SO4AEROLD(J)
c      AERSOL(J,1) = AERSOLOLD(J,1)
c      WFALL(J,1) = WFALLOLD(J,1)
c 3500 RPAR(J,1) = RPAROLD(J,1)
c
c      DO 3503 J=65,NZ
c      DO 3502 I=1,NQ
c 3502 USOL(I,J) = USOL(I,64)
c      DEN(J) = DENOLD(64)
c      SO4AER(J) = SO4AEROLD(64)
c      AERSOL(J,1) = AERSOLOLD(64,1)
c      WFALL(J,1) = WFALLOLD(64,1)E PHOTO
c 3503 RPAR(J,1) = RPAROLD(64,1)
C
c
c at May 6 2019, JL get ride of CL 
c      READ(67,400) UOLD,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
c 400  FORMAT(1P8E12.5)  13   1970.5     6.10485E-03   0.00000E+00   0.00000E+00   6.05263E+15   0.00000E+00

c      close(67)

C      READ(67,400) Usol,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
C 400  FORMAT(1P8E12.5)
C      close(67)

c-AJ 08/31/2022 ADD CL back read the old atm_composition data
      READ(67,400) USOL,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 400  FORMAT(1P8E12.5)
      close(67)
c      READ(69,400) UOLD2,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
c      close(69)
c      DO J = 1,NZ
c      DO I = 1,21
c      USOL(I,J) = UOLD(I,J)
c      END DO
c      DO I = 22,28
c      USOL(I,J) = UOLD2(I,J)
c      END DO
c      USOL(29,J) = UOLD(22,J)
c      DO I = 30,35
c      USOL(I,J) = UOLD(I-7,J)
c      END DO
c      USOL(36,J) = UOLD(29,J)
c      END DO

c      print 3000, ((usol(i,j),i=1,12),j=1,100,5)
c 3000 format(1p12e12.5)
c      print *
c      print 3000, ((usol(i,j),i=13,24),j=1,100,5)
c      print *
c      print 3000, ((usol(i,j),i=25,36),j=1,100,5)

C-AJ 03/10/2022 Read the Gaussian points and weights for global integration
      call gaussian_data(xi,wi,NumGau)

c      do J = 1,NZ
c      do I = 1,28 
c      USOL(I,J) = Uold(I, J)
c      end do 
c      USOL(29,J) = 3.55E-4 
c      end do 

C-PL CHANGE O2 RATIO HERE TO MAKE SURE FO2 GET THE RIGHT NUMBER 05/24/2019
       DO I=1,NZ
       USOL(LO2,I)=FO2 * USOL(LO2,I)/USOL(LO2,1)
       USOL(LCO2,I)=FCO2 * USOL(LCO2,I)/USOL(LCO2,1)
       O2(I) = USOL(LO2,I)
       CO2(I) = USOL(LCO2,I)
       END DO  

c JL: keep the FO2 and FCO2 factor but be mindful of the rewrite here 
!       FO2  = usol(LO2,1) 
!       FCO2 = USOL(LCO2,1) 

c      USOL(22,J) = 0.21
c      USOL(23,J) = Uold(29,J)
c      USOL(24,J) = Uold(30,J)
c      USOL(25,J) = Uold(31,J)
c      USOL(26,J) = Uold(32,J)
c      USOL(27,J) = Uold(33,J)
c      USOL(28,J) = Uold(34,J)
c      end do 
 

C JK 1/24/16 Change the temperature profile for low-O2 atmospheres
C     DO J=1,NZ
C Use this for O2=0.1 PAL
C     T(J) = TP1(J)
C Use this for O2<0.1 PAL
C      IF(J.GE.16) T(J) = amin1(T(J),200.)
C     ENDDO
C      Print 5000
C 5000 Format(/1x,'T =') 
C      Print 5001, T
C 5001 format(1p10e10.3) PL PRINT OUT TEMPERATURE WE NEED HERE

C-AJ 01/15/2022 Change the temperature profile for low-O2 atmospheres
      IF (FO2 .EQ. 0.021) THEN
      DO J=1,NZ
            T(J) = TP2(J)
      END DO
      END IF
      IF (FO2 .LT. 0.021) THEN
      DO J=1,15
            T(J) = TP2(J)
      END DO
      DO J=16,NZ
            T(J) = 210.
      END DO   
      END IF
      PRINT*,T(100)
      if (ICOUPLE.eq.1) then
	DO J=1, NZ
c Read the temperature and water profiles from the climate code
	  READ(71,*) Z(J),T(J),water(J)
	END DO
      close(71)
      endif
 351  FORMAT(I3,1PE10.3, 1PE12.3, 1PE12.3)

C-KK  Surface mixing ratios to share with the climate code
c  FAR is not needed on this code but it must be transfered to 
c  the climate model.
        READ(72,*) FAR                  !Argon
	READ(72,*) FCH4			!Methane
	READ(72,*) JTROP		!Tropopause layer
        READ(72,*) O3OLD                !Former O3 column depth
      close(72)

C   Initial constant mixing ratios used for present Earth
c       if (STARR=="Sun".and.FO2.gt.0.20) then
c           READ(67,400) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR

       if(INIT.eq.1) goto 77
c      Surface mixing ratios
         if(LBOUND(9).eq.1)FH2 = 5.5E-7
         if(LBOUND(11).eq.1)FCO = 9.0E-8
         if(LBOUND(14).eq.1)FN2O = 3.0E-7
C AJ 08/30/2022 ADD CL BACK
         if(LBOUND(23).eq.1)FCH3CL = 5.0E-10

       DO i = 1, NZ
        if(LBOUND(9).eq.1) then
           USOL(LH2, i) = USOL(LH2,i)*(FH2/USOL(LH2,1))
        endif
       if(LBOUND(10).eq.1) then
           USOL(LCH4,i) = USOL(LCH4,i)*(FCH4/USOL(LCH4,1))
        endif
        if(LBOUND(11).eq.1)then
            USOL(LCO,i) = USOL(LCO,i)*(FCO/USOL(LCO,1))
        endif
        if(LBOUND(14).eq.1) then
            USOL(LN2O,i) = USOL(LN2O,i)*(FN2O/USOL(LN2O,1))
        endif
C-AJ 08/30/2022 ADD CL BACK
        if(LBOUND(23).eq.1)then
            USOL(LCH3CL,i) = USOL(LCH3CL,i)*(FCH3CL/USOL(LCH3CL,1))
        endif
        END DO
   77   continue       
c        endif
c-AJ 12/18/2021 Run differenct methane levels
c       USOL(LCH4,1)=0.8e-6!.5e-4
c       DO J = 2, NZ
c            USOL(LCH4,J) = USOL(LCH4,J)/2.!*5.
c       END DO       
C-KK	Constant flux for long-lived gases. For details see the model 
c       description on: Segura et al. (2003) Astrobiology, 3(4), 689-708.
C-KK	Only used at O2 levels lower than 1 PAL and for planets around 
C-KK    other stars
C-KK	9 - H2   10 - CH4   11 - CO  14 - N2O   23 - CH3Cl
         if(LBOUND(9).eq.2)SGFLUX(9) = -4.55E+09!-4.28E+09!-4.17E+09!-4.05E+09!-4.36E+09!-9.43E+09 ! !-8.61E9
         if(LBOUND(10).eq.2)SGFLUX(10) = 1.14E+11!1.20E+11!1.16E11!1.26E+11!1.80E+11!!2.02E11
c        if(LBOUND(10).eq.2)SGFLUX(10) = 1.5E11  !value for AD Leo
         if(LBOUND(11).eq.2)SGFLUX(11) = 2.21E+11!2.29E+11!2.37E+11!3.02E+11!!3.46E11
         if(LBOUND(14).eq.2)SGFLUX(14) = 1.01E+09!9.93E+08!9.99E+08!9.86E8!9.78E8!9.79E+08!1.14E+09!!1.12E9

c in may 6, 2019, JL removed Cl 
C-AJ 08/30/2022 ADD CL BACK
         if(LBOUND(23).eq.2)SGFLUX(23) = 2.8e8!2.92E8!1.e7

         IF(GauFLAG.EQ.1) THEN
          write(90,*)"SOLAR ZENITH ANGLE (deg) = ", ZY
         ELSE
          write(90,*)"Gaussian points = 8"
         END IF
	 write(90,*)"FIXED SPECIES MIXING RATIOS:"
	 write(90,*)"   O2 = ",FO2," CO2 = ",FCO2
         write(90,*)
         write(90,*)"   CH4 = ",FCH4
c	 print*,"   H2 = ",FH2," CO = ",FCO," N2O = ",FN2O
C-AJ 10/27 ADD
c  310 FORMAT(4(2X,A,2X,1PE12.5))
c       write(90,*)"FIXED SPECIES CONSTANT FLUX:"
c      write(90,310)'H2=',SGFLUX(9),'CO=',SGFLUX(11),
c     2 'N2O= ',SGFLUX(14),'CH4= ',SGFLUX(10)
C-AJ 09/03/2022 Test model sensitivity to CL
        IF(CLFLAG.EQ.0) THEN
         SGFLUX(23) = SGFLUX(23)
         WRITE(90,777)'Standard CH3CL flux = ',SGFLUX(23)
 777     FORMAT(A23,1PE10.3)
        ELSE
         SGFLUX(23) = SGFLUX(23)/30
         WRITE(90,777)'Low CH3CL flux      = ', SGFLUX(23)
        END IF
c
c in May 6, 2019, JL remake species list without cl 

       ISPEC(1) = 4HH2CO
       ISPEC(2) = 1HO
       ISPEC(3) = 3HH2O
       ISPEC(4) = 2HOH
       ISPEC(5) = 3HHO2
       ISPEC(6) = 4HH2O2
       ISPEC(7) = 2HO3
       ISPEC(8) = 1HH
       ISPEC(9) = 2HH2
       ISPEC(10) = 3HCH4
       ISPEC(11) = 2HCO
       ISPEC(12) = 6HCH3OOH
       ISPEC(13) = 5HCH3O2
       ISPEC(14) = 3HN2O
       ISPEC(15) = 2HNO
       ISPEC(16) = 3HNO2
       ISPEC(17) = 4HHNO2
       ISPEC(18) = 4HHNO3
       ISPEC(19) = 6HHO2NO2
       ISPEC(20) = 3HNO3
       ISPEC(21) = 4HN2O5
C-AJ 08/30/2022 ADD CL BACK
       ISPEC(22) = 5HCL2O2
       ISPEC(23) = 5HCH3CL
       ISPEC(24) = 4HHOCL
       ISPEC(25) = 2HCL
       ISPEC(26) = 3HCLO
       ISPEC(27) = 3HHCL
       ISPEC(28) = 6HCLONO2

       ISPEC(29) = 2HO2      
       ISPEC(30) = 3HH2S
       ISPEC(31) = 2HHS
       ISPEC(32) = 2HSO
       ISPEC(33) = 3HSO2
       ISPEC(34) = 5HH2SO4
       ISPEC(35) = 3HHSO
       ISPEC(36) = 3HCO2

C       ISPEC(22) = 2HO2      
C       ISPEC(23) = 3HH2S
C       ISPEC(24) = 2HHS
C       ISPEC(25) = 2HSO
C       ISPEC(26) = 3HSO2
C       ISPEC(27) = 5HH2SO4
C       ISPEC(28) = 3HHSO
C       ISPEC(29) = 3HCO2
C
C   TRIDIAGONAL SOLVER
       ISPEC(37) = 6HSO4AER
c       ISPEC(30) = 6HSO4AER
C
C   Short-lived species
       ISPEC(38) = 4HCH21
       ISPEC(39) = 4HCH23
       ISPEC(40) = 3HO1D
       ISPEC(41) = 3HCH3
       ISPEC(42) = 4HH3CO
       ISPEC(43) = 3HHCO
       ISPEC(44) = 1HN
       ISPEC(45) = 4HNOCL
       ISPEC(46) = 5HCLONO
       ISPEC(47) = 4HCLO2
       ISPEC(48) = 3HCL2
       ISPEC(49) = 1HS
       ISPEC(50) = 4HSO21
       ISPEC(51) = 4HSO23
       ISPEC(52) = 4HHSO3
       ISPEC(53) = 3HSO3
C
C   Inert species
       ISPEC(54) = 2HS2
       ISPEC(55) = 2HN2
       ISPEC(56) = 2HHV
       ISPEC(57) = 1HM
C   Short-lived species      
c       ISPEC(31) = 4HCH21
c       ISPEC(32) = 4HCH23
c       ISPEC(33) = 3HO1D
c       ISPEC(34) = 3HCH3
c       ISPEC(35) = 4HH3CO
c       ISPEC(36) = 3HHCO
c       ISPEC(37) = 1HN
c       ISPEC(38) = 1HS
c       ISPEC(39) = 4HSO21
c       ISPEC(40) = 4HSO23
c       ISPEC(41) = 4HHSO3
c       ISPEC(42) = 3HSO3
C
C   Inert species
c       ISPEC(43) = 2HS2
c       ISPEC(44) = 2HN2
c       ISPEC(45) = 2HHV
c       ISPEC(46) = 1HM
      do i=1,NSP
       NUML(i) = 0
       NUMP(i) = 0
      enddo

C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ (61,200) JCHEM
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      close(61)
 
C NOTE: This command will print the list of chemical reactions in the 
c       output file 
c      PRINT 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
      KJAC = LDA*NEQ
c      PRINT 202,NQ,NZ,KJAC
 202  FORMAT(//1X,'NQ=',I2,5X,'NZ=',I3,5X,'KJAC=',I7)
C     
    
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
c      Print *,'Reading in the reactions'
      DO 5 J=1,NR
c      print *,'J =',J
      DO 5 M=1,5
      IF(JCHEM(M,J).EQ.' ' ) GO TO 5
C-AP Change NSP1 to NSP2 because we added M to the species list
      DO 6 I=1,NSP2
      IF(JCHEM(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      GO TO 25
   5  CONTINUE


C 
C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
c      print *,'Exchanging names for numbers'
      DO 7 M=1,2
      N = 3-M
      DO 7 J=1,NR
c      print *,'J =',J

      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
      NUML(I) = NUML(I) + 1
      IF(NUML(I).GT.NMAX) GO TO 20
      K = NUML(I)
      ILOSS(1,I,K) = J
      ILOSS(2,I,K) = JCHEM(N,J)
   7  CONTINUE
   
C
      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE

c-08/31/2022 debug
c      PRINT *
c      print *, 'Print out the chemical production matrix'
c      DO I=1,NSP
c      PRINT 5001,I
c 5001 FORMAT(/'I =',I3)
c      print 5000,(IPROD(I,K),K=1,30)
c 5000 FORMAT(30(I3,1X))
c      ENDDO

c      print *, 'I =', I
c      print *, 'NUMP =', NUMP
c      print *, 'J=',J
c       print *, 'K=',K
C-AP
C
c** Read the eddy diffusion profile
	do J=1,NZ
	 READ(64,245) EDD(J),Z(J)
	end do 
 245  FORMAT(2(1PE11.3, 1X))
      close(64)
      LTIMES = 0       !Counter for the photorate subroutine
      ISULF = 0
      VOLFLX = 3.E9    !Volcanic flux
c      print *, 'Z before call',Z

      CALL GRID
c      print *, 'Z after call',Z
      DZ = Z(2) - Z(1)

      CALL DENSTY(O2,CO2,P0)
      

      CALL RATES

      CALL DIFCO(O2,CO2)

c  If the star is other than the Sun or there is a flare
c  the fluxes are saved here
      do i=1,108
       fluxsave(i)=FLUX(i)
      enddo     

c  Read the data for the photochemistry and the UV flux of the Sun
      CALL READPHOTO
      if(STARR.ne.'Sun') then
        do i=1,108
         FLUX(i)= fluxsave(i)
        enddo       
      endif
   
C AJ 12/21/2022 read correlated k table from srbands.f
c      print *, 'hello world'
      CALL SRBANDS(T)

c** Water calculation

      JTROP = ZTROP/DZ + 0.01

      CALL PSATRAT(H2O)
c      print *, 'hello world' 
       DO 23 J=1,JTROP 
  23    USOL(LH2O,J) = H2O(J)
 
       
       if(ICOUPLE.eq.1) then
C-KK	This is added to make sure that tropospheric water is being
C-KK	handled consistently. The #s are imported from the climate model.
        DO J = 1, NZ
         USOL(LH2O,J) = water(J)
        END DO
       endif       !end of water calculation


      CALL LTNING(FO2,FCO2,P0)
c      print *, 'After call to LTNING'
      CALL AERTAB
      NZ1 = NZ - 1
      HA = 1.38E-16*T(NZ)/(1.67E-24*28.*980.)


      DTINV = 1./DT
      TIME = 0.
C
C ***** PRINT OUT INITIAL DATA *****
c in may 10 2019, JL common out call for test 

      CALL OUTPUTP(0,NSTEPS,0.,FLOW)

C    
C ***** SET JACOBIAN DIMENSIONING PARAMETERS *****
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   PRINT OUT RESULTS EVERY NPR TIME STEPS
      NPR = NSTEPS
      PRN = NPR
C
C   DO PHOTORATES EVERY MP TIME STEPS
      NPHOT = 0
      MP = 3
      PM = MP
      NN = 0 
 
C
C ***** START THE TIME-STEPPING LOOP *****
C	STARSHIP
      DO 1 N=1,NSTEPS
      TIME = TIME + DT
      NN = NN + 1
      MS = (N-1)/MP
      SM = (N-1)/PM
      IF(NN.EQ.NSTEPS) SM = MS
      IF(SM-MS.GT.0.01) GO TO 18
      IF(N.GT.1 .AND. TIME.LT.1.E4) GO TO 18
c      print *, 'JTROP =', JTROP 
c      print *, 'DEN =', DEN
C
C   STORE ABSORBERS USED TO BLOCK OUT SOLAR UV RADIATION
      DO 35 I=1,NZ
      H2O(I) = ABS(USOL(LH2O,I))
      O3(I) = ABS(USOL(LO3,I))
      O2(I) = USOL(LO2,I) !FO2
      CO2(I) = USOL(LCO2,I) !FCO2
C-AP
      FSO2(I) = ABS(USOL(LSO2,I))
      H2S(I) = ABS(USOL(LH2S,I))
C-AP
      CH4(I) = ABS(USOL(LCH4,I))


c      print *, 'H2O(I) =', H2O 
c      print *, 'O3(I) = ', O3 
c      print *, 'FSO2(I) = ', FSO2
c      print *, 'H2S(I) = ', H2S 
c      print *, 'CH4(I) = ', CH4
  35  CONTINUE
            
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1    
!      CALL DENSTY(O2,CO2,P0)   		
!      CALL DIFCO(O2,CO2)
C-AJ 03/04/2022 ADD FOR GAUSSIAN INTEGRATION
c-AJ 03/10/2022 Find the right right row in the matrix

c-AJ 09/02/2022 ADD FLAG FOR USING GAUSS
      IF (GauFLAG.EQ.1) THEN
       CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)
      ELSE 
      DO i = 1,11
      isave=i
      if(NumGau(i).EQ.NZA) exit
      END DO 
c      print*,NZA
c      print*,'gaussian points = ', NumGau(isave)
c      print*,isave
c isave holds the correct row number for the matrix    
      DO 136 IPHOT=1,NZA
      radians=xi(isave,IPHOT)
      ZY=acos(radians)*180./3.14159
      AGL=.5
c      print*,'radians = ',radians
c      print*,'solar zenith angle = ', ZY
      WEIGHT=wi(isave,IPHOT)

      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)

C-AJ 03/04/2022 ADD FOR GAUSSIAN INTEGRATION      
      DO 137 J=1,NZ
      PO2M(J,IPHOT) = PO2(J)*WEIGHT
      PO2DM(J,IPHOT) = PO2D(J)*WEIGHT
      PO3M(J,IPHOT) = PO3(J)*WEIGHT
      PO3DM(J,IPHOT) = PO3D(J)*WEIGHT
      PH2OM(J,IPHOT) = PH2O(J)*WEIGHT
      PH2O2M(J,IPHOT) = PH2O2(J)*WEIGHT
      PCO2M(J,IPHOT) = PCO2(J)*WEIGHT
      PCO2DM(J,IPHOT) = PCO2D(J)*WEIGHT
      PHCOM(J,IPHOT) = PHCO(J)*WEIGHT
      PH2M(J,IPHOT) = PH2(J)*WEIGHT
      PHO2M(J,IPHOT) = PHO2(J)*WEIGHT
      PCH4M(J,IPHOT) = PCH4(J)*WEIGHT
      PMOOHM(J,IPHOT) = PMOOH(J)*WEIGHT
      PN2OM(J,IPHOT) = PN2O(J)*WEIGHT
      PHNO3M(J,IPHOT) = PHNO3(J)*WEIGHT
      PNOM(J,IPHOT) = PNO(J)*WEIGHT
      PNO2M(J,IPHOT) = PNO2(J)*WEIGHT
      PHNO4M(J,IPHOT) = PHNO4(J)*WEIGHT
      PNO3M(J,IPHOT) = PNO3(J)*WEIGHT
      PN2O5M(J,IPHOT) = PN2O5(J)*WEIGHT
      PSO2M(J,IPHOT) = PSO2(J)*WEIGHT
      PSO21M(J,IPHOT) = PSO21(J)*WEIGHT
      PSO23M(J,IPHOT) = PSO23(J)*WEIGHT
      PSOM(J,IPHOT) = PSO(J)*WEIGHT
      PH2SM(J,IPHOT) = PH2S(J)*WEIGHT
      PH2SO4M(J,IPHOT) = PH2SO4(J)*WEIGHT

C-AJ 08/30/2022 ADD CL BACK
      PCH3CLM(J,IPHOT) = PCH3CL(J)*WEIGHT
      PCL2M(J,IPHOT) = PCL2(J)*WEIGHT
      PCLO2M(J,IPHOT) = PCLO2(J)*WEIGHT
      PHCLM(J,IPHOT) = PHCL(J)*WEIGHT
      PHOCLM(J,IPHOT) = PHOCL(J)*WEIGHT
      PNOCLM(J,IPHOT) = PNOCL(J)*WEIGHT
      PCLONOM(J,IPHOT) = PCLONO(J)*WEIGHT
      PCLONO2M(J,IPHOT) = PCLONO2(J)*WEIGHT
      PCL2O2M(J,IPHOT) = PCL2O2(J)*WEIGHT

 137  CONTINUE
 136  CONTINUE

      DO 138 JJ=1,NZ
C      AR(23,JJ) = PO2DM(JJ,1)+PO2DM(JJ,2)+PO2DM(JJ,3)+PO2DM(JJ,4)
C     2 +PO2DM(JJ,5)+PO2DM(JJ,6)
      AR(23,JJ) = SUM(PO2DM(JJ,1:NZA))
      AR(24,JJ) = SUM(PO2M(JJ,1:NZA))
      AR(25,JJ) = SUM(PH2OM(JJ,1:NZA))
      AR(26,JJ) = SUM(PO3DM(JJ,1:NZA))
      AR(27,JJ) = SUM(PO3M(JJ,1:NZA))
      AR(28,JJ) = SUM(PH2O2M(JJ,1:NZA))
      AR(29,JJ) = SUM(PCO2M(JJ,1:NZA))
      AR(38,JJ) = SUM(PH2M(JJ,1:NZA))
      AR(39,JJ) = SUM(PHCOM(JJ,1:NZA))
      AR(42,JJ) = SUM(PCO2DM(JJ,1:NZA))
      AR(50,JJ) = SUM(PHO2M(JJ,1:NZA))
      AR(51,JJ) = SUM(PCH4M(JJ,1:NZA))
      AR(52,JJ) = SUM(PMOOHM(JJ,1:NZA))
      AR(53,JJ) = SUM(PN2OM(JJ,1:NZA))
      AR(54,JJ) = 1.7E-3
      AR(55,JJ) = SUM(PHNO3M(JJ,1:NZA))
      AR(56,JJ) = SUM(PNOM(JJ,1:NZA))
      AR(57,JJ) = SUM(PNO2M(JJ,1:NZA))
      AR(95,JJ) = SUM(PHNO4M(JJ,1:NZA))
c JL: CL rip cont.
C-AJ 08/30/2022 ADD CL BACK
      AR(101,I) = SUM(PCH3CLM(JJ,1:NZA))
      AR(132,I) = SUM(PCL2M(JJ,1:NZA))
      AR(133,I) = SUM(PCLO2M(JJ,1:NZA))
      AR(134,I) = SUM(PHCLM(JJ,1:NZA))
      AR(135,I) = SUM(PHOCLM(JJ,1:NZA))
      AR(136,I) = SUM(PNOCLM(JJ,1:NZA))
      AR(137,I) = SUM(PCLONOM(JJ,1:NZA))
      AR(138,I) = SUM(PCLONO2M(JJ,1:NZA))
      AR(142,I) = SUM(PCL2O2M(JJ,1:NZA))
      AR(151,I) = SUM(PNO3M(JJ,1:NZA))
      AR(153,I) = SUM(PN2O5M(JJ,1:NZA))

C      AR(105,JJ) = SUM(PNO3M(JJ,1:NZA))
C      AR(107,JJ) = SUM(PN2O5M(JJ,1:NZA))
      AR(156,JJ) = 0.
      AR(157,JJ) = 0.7*SUM(PSO2M(JJ,1:NZA))
      AR(158,JJ) = SUM(PH2SM(JJ,1:NZA))
      AR(186,JJ) = SUM(PSO21M(JJ,1:NZA))
      AR(187,JJ) = SUM(PSO23M(JJ,1:NZA))
      AR(188,JJ) = SUM(PH2SO4M(JJ,1:NZA))
      AR(189,JJ) = 0.
      AR(211,JJ) = SUM(PHO2M(JJ,1:NZA))
      AR(182,JJ) = 0.3*SUM(PSO2M(JJ,1:NZA))

C      AR(110,JJ) = 0.
C      AR(111,JJ) = 0.7*SUM(PSO2M(JJ,1:NZA))
C      AR(112,JJ) = SUM(PH2SM(JJ,1:NZA))
C      AR(140,JJ) = SUM(PSO21M(JJ,1:NZA))
C      AR(141,JJ) = SUM(PSO23M(JJ,1:NZA))
C      AR(142,JJ) = SUM(PH2SO4M(JJ,1:NZA))
C      AR(143,JJ) = 0.
C      AR(165,JJ) = SUM(PHO2M(JJ,1:NZA))
C      AR(136,JJ) = 0.3*SUM(PSO2M(JJ,1:NZA))

 138  CONTINUE

C-AJ 08/31/2022 CHECK THE PHOTOLYSIS RATES IN AR MATRIX
C      PRINT *,'Check Photolysis rates in AR matrix in main code'
C      PRINT 999, (AR(i,100),i=1,NR)
C      PRINT 999, (AR(NR,I),I=1,NZ,5)

      END IF
      PRINT *, 'AR(24,61) IN THE MAIN',AR(24,61)
      CALL AERCON(H2O)
c      print *, ' po3d(1) = ', PO3D(1) 

C
C ****** TIME-DEPENDENT BOUNDARY CONDITIONS ********
C
C  UPPER BOUNDARY
C  Escape of hydrogen: VEFF(H) = (Bi/Ha)/Nt, Nt=total number density
      BOVERH = DI(LH,NZ)*DEN(NZ)/HA
      VEFF(LH) = BOVERH/DEN(NZ)
      BOVERH2 = DI(LH2,NZ)*DEN(NZ)/HA
      VEFF(LH2) = BOVERH2/DEN(NZ)
      BOVERCH4 = DI(LCH4,NZ)*DEN(NZ)/HA
      VEFF(LCH4) = BOVERCH4/DEN(NZ)
      BOVERH2O = DI(LH2O,NZ)*DEN(NZ)/HA
      VEFF(LH2O) = BOVERH2O/DEN(NZ)
C


      IF(NN.EQ.NSTEPS) write(90, 63) BOVERH,VEFF(LH),BOVERH2,VEFF(LH2)
     & ,BOVERCH4,VEFF(LCH4),BOVERH2O,VEFF(LH2O),VEFF(LO2),VEFF(LCO2)
  63  FORMAT(/'Information on upper boundary conditions'/'BOVERH=',
     2  1PE10.3,' VEFF(LH)= ',E10.3,' BOVERH2=',
     3   E10.3,' VEFF(LH2)= ',E10.3,/'BOVERCH4= ',E10.3,' VEFF(LCH4)= '
     4  ,E10.3,' BOVERH2O- ',E10.3,' VEFF(LH2O)= ',E10.3,/'VEFF(LO2)= '
     5  ,E10.3,3x,'VEFF(LCO2)= ',E10.3)
C
C PL 07/2019 Calculate SMFLUX(LCO) 

c-AJ 08/31/2022 PO2 should be AR(24,NZ)now
c      VO2 = (PO2(NZ) + PO2D(NZ)) * HA
c      VCO2 = (PCO2(NZ) + PCO2D(NZ)) * HA
      VO2 = (AR(24,NZ) + AR(23,NZ)) * HA
      VCO2 = (AR(29,NZ) + AR(42,NZ)) * HA
      VEFF(LO2) = VO2
      VEFF(LCO2) = VCO2
      SMFLUX(LCO) = - VCO2*CO2(NZ)*DEN(NZ)   
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ) - 2.*VO2*O2(NZ)*DEN(NZ)
     2  + 2.*VEFF(LCH4)*USOL(LCH4,NZ)*DEN(NZ)

!      SMFLUX(LCO) = - VCO2*CO2(NZ)*DEN(NZ)
!      SMFLUX(LO) = SMFLUX(LCO) - 2.*VO2*O2(NZ)*DEN(NZ)

C
      NMP = NSTEPS - MP
      IF (NN.gt.1 .and. nn.LT.NSTEPS) GO TO 18

      write(90,97)
  97  FORMAT(//1X,'PHOTOLYSIS RATES')
      write(90,98) 
  98  FORMAT(/5X,'Z',7X,'PO2',6X,'PO2D',5X,'PCO2',5X,'PCO2D',4X,
     2  'PH2O',5X,'PO3',6X,'PO3D',5X,'PH2O2',4X,'PHCO',5X,'PH2',
     3  6X,'PHO2')
c AJ 0826/22 REWRITE THE PHOTOLYSIS RATES
c      write(90,99)(Z(I),PO2(I),PO2D(I),PCO2(I),PCO2D(I),PH2O(I),
c     2  PO3(I),PO3D(I),PH2O2(I),PHCO(I),PH2(I),PHO2(I),I=1,NZ,
c     3  3)
      write(90,99)(Z(I),AR(24,I),AR(23,I),AR(29,I),AR(42,I),AR(25,I),
     2  AR(27,I),AR(26,I),AR(28,I),AR(39,I),AR(38,I),AR(50,I),I=1,NZ,3)
  99  FORMAT(2X,1P12E9.2)
      write(90,198)

 198  FORMAT(/5X,'Z',6X,'PCH4',5X,'PCH3OOH',2X,'PN2O',5X,'PHNO3',4X,
     2  'PNO',6X,'PNO2',5X,'PHNO4',4X,'PCCL3F',3X,'PCCL2F2',2X,
     3  'PCCL4',4X,'PCH3CL')
C      write(90,99) (Z(I),PCH4(I),PMOOH(I),PN2O(I),PHNO3(I),PNO(I),
C     2  PNO2(I),PHNO4(I),PCCL3F(I),PCCL2F2(I),PCCL4(I),PCH3CL(I),
C     3  I=1,NZ,3)
      write(90,99) (Z(I),AR(51,I),AR(52,I),AR(53,I),AR(55,I),AR(56,I),
     2  AR(57,I),AR(95,I),PCCL3F(I),PCCL2F2(I),PCCL4(I),AR(101,I),
     3  I=1,NZ,3)
      write(90,197)

 197  FORMAT(/5X,'Z',6X,'PMCCL3',3X,'PCL2',5X,'PHOCL',4X,'PNOCL',4X,
     2  'PCLONO',3X,'PCLONO2',2X,'PCLO2',4X,'PHCL')
C      write(90,199) (Z(I),PMCCL3(I),AR(132,I),AR(135,I),AR(136,I),
C     2  AR(137,I),AR(138,I),AR(133,I),AR(134,I),I=1,NZ,3)
C 199  FORMAT(2X,1P9E9.2)
      write(90,199) (Z(I),PMCCL3(I),AR(132,I),AR(135,I),AR(136,I),
     2  AR(137,I),AR(138,I),AR(133,I),AR(134,I),I=1,NZ,3)
 199  FORMAT(2X,1P9E9.2)
      write(90,298)
 298  format(/5x,'Z',6x,'PNO3',5x,'PN2O5',4x,'PCL2O2',3x,'PSO2',5x,
     2  'PH2S',5x,'PSO21',4x,'PSO23',4x,'PH2SO4')
C      write(90,299) (z(i),pno3(i),pn2o5(i),pcl2o2(i),pso2(i),ph2s(i),
C     2  pso21(i),pso23(i),ph2so4(i),i=1,nz,3)
      write(90,299) (z(i),AR(151,i),AR(153,i),AR(142,i),AR(157,i),
     2  AR(158,i),AR(186,i),AR(187,i),AR(188,i),i=1,nz,3)
 299  format(2x,1p9e9.2)

  18  CONTINUE

      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL SEDMNT(FSULF,IDO)
c      print *,'After SEDMNT'
      DO J=1,NZ
        AERSOL(J,1) = SO4AER(J)*DEN(J)/CONVER(J,1)
      ENDDO                                                               
C

C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,LDA
      DO 17 K=1,NEQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NEQ
  19  RHS(K) = 0.

C
C     (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
c      print *
c      print *,'Before call to DOCHEM'

c      print *, 'USOl =', USOL
      CALL DOCHEM(FVAL,IDO)
c      print *, 'DEN = ', DEN
c      print *, 'DEN = ', DEN
      DO 9 I=1,NQ
      DO 9 J=1,NZ
      K = I + (J-1)*NQ
      RHS(K) = FVAL(I,J)
   9  USAVE(I,J) = USOL(I,J)
C

      DO 3 I=1,NQ
      DO 11 J=1,NZ
      R(J) = EPSJ * ABS(USOL(I,J))
  11  USOL(I,J) = USAVE(I,J) + R(J)
      CALL DOCHEM(FV,0)
C
      DO 12 M=1,NQ
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQ
  12  DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)
C
      DO 10 J=1,NZ
  10  USOL(I,J) = USAVE(I,J)
   3  CONTINUE
C
C   COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO 13 I = 1,NQ
      DO 14 J=2,NZ1
      K = I + (J-1)*NQ
      RHS(K) = RHS(K) - DD(I,J)*USOL(I,J) 
     2  + (DU(I,J) + DHU(I,J))*USOL(I,J+1) 
     3  + (DL(I,J) - DHL(I,J))*USOL(I,J-1)      
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DD(I,J)
      DJAC(KU,K+NQ) = - DU(I,J) - DHU(I,J)
  14  DJAC(KL,K-NQ) = - DL(I,J) + DHL(I,J)
  13  CONTINUE
C
C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQ
      U(K) = USOL(K,1)
      LB = LBOUND(K)
C
C   CONSTANT DEPOSITION VELOCITY
      IF(LB.NE.0) GO TO 16
      RHS(K) = RHS(K) + (DU(K,1) + DHU(K,1))*USOL(K,2) - DU(K,1)*U(K)
     2  - (VDEP(K)/DZ - HI(K,1)/(2.*DZ))*U(K)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(K,1) + VDEP(K)/DZ
     2  - HI(K,1)/(2.*DZ)
      DJAC(KU,K+NQ) = - DU(K,1) - DHU(K,1)
      GO TO 15
C
C   CONSTANT MIXING RATIO
  16  IF(LB.NE.1) GO TO 31
      RHS(K) = 0.
      DO 36 M=1,NQ
      MM = KD + K - M
  36  DJAC(MM,M) = 0.
      DJAC(KU,K+NQ) = 0.
      DJAC(KD,K) = DTINV + DU(K,1)
      GO TO 15
C
C   CONSTANT UPWARD FLUX
C-PL fixed sign error in last term of DJAC from SH
  31  CONTINUE
      RHS(K) = RHS(K) + (DU(K,1) + DHU(K,1))*USOL(K,2) - DU(K,1)*U(K)
     2   + HI(K,1)*U(K)/(2.*DZ) + SGFLUX(K)/DEN(1)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(K,1) - HI(K,1)/(2.*DZ)
      DJAC(KU,K+NQ) = - DU(K,1) - DHU(K,1)
  15  CONTINUE
C
C ***** UPPER BOUNDARY CONDITIONS *****
      DO 30 I=1,NQ
      U(I) = USOL(I,NZ)
      K = I + NZ1*NQ
      MB = MBOUND(I)
C
C   CONSTANT EFFUSION VELOCITY
      IF(MB.NE.0) GO TO 29
      RHS(K) = RHS(K) + (DL(I,NZ) - DHL(I,NZ))*USOL(I,NZ1) 
     2  - DL(I,NZ)*U(I) - (VEFF(I)/DZ + HI(I,NZ)/(2.*DZ))*U(I)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(I,NZ) + VEFF(I)/DZ
     2  + HI(I,NZ)/(2.*DZ)
      DJAC(KL,K-NQ) = - DL(I,NZ) + DHL(I,NZ)
      GO TO 30
C
C   CONSTANT DOWNWARD FLUX
  29  CONTINUE
      RHS(K) = RHS(K) + (DL(I,NZ) - DHL(I,NZ))*USOL(I,NZ1)
     2  - DL(I,NZ)*U(I) - HI(I,NZ)*U(I)/(2.*DZ) 
     3  - SMFLUX(I)/DEN(NZ)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(I,NZ) + HI(I,NZ)/(2.*DZ) 
      DJAC(KL,K-NQ) = - DL(I,NZ) + DHL(I,NZ)
  30  CONTINUE
C
      DO 33 J=1,NZ
      IF(Z(J).GT.ZTROP) GO TO 34
      K = 3 + (J-1)*NQ
      RHS(K) = 0.
      DO 32 M=1,NQ
      MM = M - 3 + KD
  32  DJAC(MM,K) = 0.
      DJAC(KD,K) = DTINV
      DJAC(KU,K+NQ) = 0.
      IF(J.EQ.1) GO TO 33
      DJAC(KL,K-NQ) = 0.
  33  CONTINUE
  34  CONTINUE
C
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INDEX)
      IF(INDEX.NE.0.) write(90,103)N,INDEX
 103  FORMAT(/1X,'N =',I3,5X,'INDEX =',I3)
      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)
C

C   COMPUTE NEW CONCENTRATIONS
      EMAX = 0.
      DO 26 I=1,NQ
      DO 26 J=1,NZ
      K = I + (J-1)*NQ
C-KK	For all runs, inc. standard O2
      IF((I.EQ.LH2S).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHS).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHSO).AND.(Z(J).GT.2.E6)) GOTO 26
      IF((I.EQ.LCH3CL).AND.(Z(J).GT.2.8E6)) GOTO 26
      IF((I.EQ.LH).AND.(Z(J).LT.3.E6)) GOTO 26
      IF (I.EQ.LSO) GOTO 26
      IF (I.EQ.LCL2O2) GOTO 26
      IF((I.EQ.LN2O5).AND.(Z(J).GT.1.E6)) GOTO 26
      IF((I.EQ.LHNO3).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LHO2NO2).AND.(Z(J).GT.5.E6)) GOTO 26
      IF((I.EQ.LHNO2).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LNO3).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LCLONO2).AND.(Z(J).GT.5.E6)) GOTO 26
      IF((I.EQ.LCH3OOH).AND.(Z(J).GT.7.E6)) GOTO 26

C
      REL(I,J) = RHS(K)/USOL(I,J)
      EREL = ABS(REL(I,J))
      EMAX = AMAX1(EMAX,EREL)
      IF(EREL.LT.EMAX) GO TO 26
      IS = I
      JS = J
      UMAX = USOL(I,J)
      RMAX = RHS(K)
  26  USOL(I,J) = USOL(I,J) + RHS(K)

C-PL  DO NOT LET [H2O] CHANGE AGTER CHEMICAL REACTION
C
      DO 4 J=1,NZ
      IF(Z(J).LT.ZTROP) USOL(3,J) = H2O(J)
   4  CONTINUE

C-AP Adding tridiagonal solver
C       TRIDIAGONAL INVERSION *****
      L=1
      I = NQ + L
      IF(I.EQ.LSO4AER) MZ = 50
C-AP      IF(I.EQ.LS8TESTAER) MZ = 40
C-AP      IF(I.EQ.LHCAER) MZ = NZ
      MZ1 = MZ - 1
      MZP1 = MZ + 1
C
C   COMPUTE ADVECTION TERMS FOR PARTICLES
      DPU(1,L) = WFALL(2,L)*DEN(2)/DEN(1)/(2.*DZ)
      DPL(NZ,L) = WFALL(NZ1,L)*DEN(NZ1)/DEN(NZ)/(2.*DZ)
      DO 38 J=2,NZ1
      DPU(J,L) = WFALL(J+1,L)*DEN(J+1)/DEN(J)/(2.*DZ)
  38  DPL(J,L) = WFALL(J-1,L)*DEN(J-1)/DEN(J)/(2.*DZ)
C                                                                              
C
C   TA = LOWER DIAGONAL, TB = DIAGONAL, TC = UPPER DIAGONAL, TY =
C   RIGHT-HAND SIDE
      DO 70 J=1,NZ
      TA(J) = 0.
      TB(J) = 0.
      TC(J) = 0.
  70  TY(J) = 0.
C
      DO 44 J=1,MZ
      TB(J) = YL(I,J)
  44  TY(J) = YP(I,J)/DEN(J)
C
      DO 45 J=2,MZ1
      TA(J) = - DL(I,J) + DPL(J,L)
      TB(J) = TB(J) + DD(I,J)
  45  TC(J) = - DU(I,J) - DPU(J,L)
C                                                                              
C   BOUNDARY CONDITIONS
      TA(MZ) = - DL(I,MZ) + DPL(MZ,L)
      TB(MZ) = TB(MZ) + DL(I,MZ) + 0.5*WFALL(MZ,L)/DZ
      TB(1) = TB(1) + DU(I,1) + (.01 - 0.5*WFALL(1,L))/DZ
      TC(1) = - DU(I,1) - DPU(1,L)
C
      CALL SGTSL(MZ,TA,TB,TC,TY,NFLAG)
C-AP      STOP
      IF (NFLAG.NE.0) write(90,401) N,NFLAG,I
 401  FORMAT(//1X,'TRIDIAGONAL SOLVER FAILED AT N =',I3,2X,
     2  'NFLAG =',I2,2X,'SPECIES #',I2)
C
      IF(I.EQ.LSO4AER) THEN
        DO 59 J=1,MZ
   59     SO4AER(J) = TY(J)
C-AP      ELSEIF(I.EQ.LS8AER) THEN
C-AP        DO 46 J=1,MZ
C-AP   46     S8(J) = TY(J)
C-AP      ELSEIF(I.EQ.LHCAER) THEN
C-AP        DO 60 J=1,MZ
C-AP          HCAER(J) = TY(J)
C-AP   60   CONTINUE
      ENDIF
C
C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION
      IF(I.EQ.LSO4AER .AND. MZ.NE.NZ) THEN
        DO 61 J=MZP1,NZ
          SO4AER(J) = SO4AER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
   61     SO4AER(J) = AMAX1(SO4AER(J),1E-100)
C-AP      ELSEIF(I.EQ.LS8AER .AND. MZ.NE.NZ) THEN
C-AP        DO 47 J=MZP1,NZ
C-AP          S8(J) = S8(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
C-AP  47      S8(J) = AMAX1(S8(J),1E-100)
C-AP      ELSEIF(I.EQ.LHCAER .AND. MZ.NE.NZ) THEN
C-AP        DO 62 J=MZP1,NZ
C-AP          HCAER(J) = HCAER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
C-AP          HCAER(J) = AMAX1(HCAER(J),1E-100)
C-AP   62   CONTINUE
      ENDIF
C-AP   58 CONTINUE

      if(NN.eq.NSTEPS) GOTO 2010
C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
      IF(EMAX.GT.0.20)  DT = 0.7*DTSAVE
      IF(EMAX.GT.0.15)  DT = 0.9*DTSAVE
      IF(EMAX.LT.0.10)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.05)  DT = 1.3*DTSAVE
      IF(EMAX.LT.0.03)  DT = 1.5*DTSAVE
      IF(EMAX.LT.0.01)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.003) DT = 5.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 10.*DTSAVE

c Adjusting time step to assure that TIME<=TSTOP
      TIME1 = TIME + DT
      if(TIME1.ge.TSTOP) then
        DT = ABS(TIME-TSTOP)
        NN = NSTEPS -1
      endif
      DTINV = 1./DT
C
2010  ISP = ISPEC(IS)
      ZMAX = Z(JS)
      IF(SM-MS.GT.0.01) GO TO 317
      write(90,100)N,EMAX,ISP,ZMAX,
     & UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
C-AP
C   COMPUTE ATMOSPHERIC OXIDATION STATE

      DO 42 I=1,NQ
      SR(I) = 0.
      DO 43 J=1,JTROP
  43  SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ
      PHIDEP(I) = VDEP(I)*USOL(I,1)*DEN(1)
  42  TLOSS(I) = SR(I) + PHIDEP(I)


      SR(LSO4AER) = 0.
      DO 48 J=1,JTROP
      SR(LSO4AER) = SR(LSO4AER) + RAINGC(LH2SO4,J)*SO4AER(J)*DEN(J)*DZ
  48  CONTINUE
      PHIDEP(LSO4AER) = (WFALL(1,1) + .01) * SO4AER(1) * DEN(1)
      TLOSS(LSO4AER) = SR(LSO4AER) + PHIDEP(LSO4AER)
C
C   COMPUTE SULFUR BUDGET AND READJUST SO2 (H2S) OUTGASSING RATE IF SO
C   DESIRED (PROGRAM IS SET UP FOR PURE SO2 OUTGASSING)
      SLOSS = TLOSS(LH2S) + TLOSS(LHS) + TLOSS(LSO) +
     2  TLOSS(LSO2) + TLOSS(LH2SO4) + TLOSS(LHSO) + 
     3  TLOSS(LSO4AER)
      SLOSSP = SLOSS - TLOSS(LSO2)
      IF (ISULF.EQ.0 .OR. TIME.LT.1.E6) GO TO 316
      IF (LBOUND(LSO2).EQ.2) SGFLUX(LSO2) = SGFLUX(LSO2) * VOLFLX/SLOSS
      IF (LBOUND(LH2S).EQ.2) SGFLUX(LH2S) = SGFLUX(LH2S) * VOLFLX/SLOSS
 316  CONTINUE
C
 317  CONTINUE
C
C   RETRY TIME STEP IF EMAX EXCEEDS 30 PERCENT
      IF(EMAX.LT.0.3) GO TO 28
      write(90,*)'RETRY TIME STEP BECAUSE EMAX=',EMAX
      DT = 0.5*DTSAVE
      TIME = TIME - DTSAVE
      if(TIME.lt.0) TIME = 0.
      DO 27 I=1,NQ
      DO 27 J=1,NZ
  27  USOL(I,J) = USAVE(I,J)
  28  CONTINUE
C
      NS = N/NPR
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF(SN-NS.GE.1.E-4) GO TO 37 !C-PL prevent print redundant info. for first step 
C
      write(90,*)
      write(90,*)'BEFOREOUT'
      write(90,*) 'N=',N,'NN=',NN,'SN=',SN,'NS=',NS



      CALL OUTPUTP(NN,NSTEPS,TIME,FLOW)
  37  CONTINUE
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
      IF(TIME.GE.TSTOP) NN = NSTEPS - 1
      if(TIME.lt.TSTOP.and.N.eq.NSTEPS) then
       write(*,'(A,I4,A)')'Time not reached after',NSTEPS,
     &  ' of the photochemical model. The run is stopping now.'
       STOP
      endif
          
   1  CONTINUE      
     
C ***** END THE TIME-STEPPING LOOP *****

  22  CONTINUE
      if(ICOUPLe.eq.1) NSTEPS =N

      OPEN(unit=81,file= DIRIO//'/atm_composition.out')
      WRITE(81,399) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 399  FORMAT(1P8E12.5)
      close(81) 

c Transfer results to the climate model (ICOUPLE=1)
      if (ICOUPLE.eq.0) then
       DO 255 I=1,NZ
c Transfer O3 and H2O
        WRITE(84,254) Z(I),PRESS(I),O3(I),H2O(I) 
 254    FORMAT(1PE9.3,3(E10.2))
 255   CONTINUE
       close(84)
      endif
C-KK  Need to find coldtrap by locating water mixing ratio minimum. 

	Jcold = 0

	DO J = 1, NZ
 	   IF (JCOLD .EQ. 0) THEN
           IF (T(J) .LT. T(J+1)) JCOLD = J
	   END IF
	END DO

	OPEN(unit=19,file= DIRIO//'/mixing_ratios.out')
c     Transfer surface mixing ratios
        WRITE(19,*) FAR
	WRITE(19,*) USOL(LCH4,1) 
	WRITE(19,*) FCO2
	WRITE(19,*) O2(1)
	WRITE(19,*) Jcold
        WRITE(19,*) O3COL
       close(19)

C   The following WRITEs out the o3 mixing ratio values so that the 
c   profiles can be compared
C
c      DO 333 I=1,NZ
c        WRITE(15,332) O3(I),Z(I)
c 332  FORMAT(E10.2,1X,1PE9.3)
c 333  CONTINUE 
C   The following WRITEs out the temperature and altitude so that 
C   the profiles can be compared.
C
c      DO 367 I=1,NZ
c	WRITE(14,366) T(I),Z(I)
c 366  FORMAT(1PE9.3,1X,1PE9.3)
c 367  CONTINUE     
c
c   MaKe file for plots
c-as Not needed for now
c      WRITE(10,699)
c 699  format(2x,'Alt',6x,'O',10x,'O3',9x,'CH4',8x,'CO',9x,'N2O',8x,
c     1  'CH3Cl',6x,'n(OH)',6x,'n(O1D)',5x,'PO2',8x,'PO3D')  
c      do 700 i=1,nz
c      zkm = z(i)/1.e5
c      po2tot = po2(i) + po2d(i)
c 700  WRITE(10,701) zkm,usol(2,i),usol(7,i),usol(10,i),usol(11,i),
c     1  usol(14,i),usol(23,i),sl(4,i),sl(31,i),po2tot,po3d(i)
c 701  format(f5.1,1x,1p10e11.4)
      GO TO 21
  20  write(90,300)I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
      GO TO 21
  25  write(90,301)IERR
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)
C
  21  CONTINUE
c      close(82)


C-PL Added the oxygen isotope calculation 
C-SH Added code from Alex Pavlov's version 
c-AJ 05/06/2022 calculate reaction rate of OH
      DO 377 JJ = 1,NZ
      RateOH1(JJ) = AR(1,JJ)*SL(LO1D,JJ)*USOL(LH2O,JJ)*DEN(JJ)
      RateOH25(JJ) = AR(25,JJ)*USOL(LH2O,JJ)*DEN(JJ)
      RateOH59(JJ) = AR(59,JJ)*USOL(LCH4,JJ)*DEN(JJ)*SL(LO1D,JJ)
      RatePO3(JJ) = AR(26,JJ)*USOL(LO3,JJ)*DEN(JJ)
  377  CONTINUE  

c-AJ 05/06/2022 WRITE REACTION RATES OF OH
      write(91,397)
 397  format(/5x,'Z',6X,'Reaction1',3x,'Reaction25',3X,'Reaction59')
      write(91,398) (z(i),RateOH1(I),RateOH25(I),RateOH59(I),I=1,NZ)
 398  format(2x,1p4e9.2)

c-AJ 08/28/2022 WRITE INTO OUPUT FOR Gaussian plot
      write(92,387) 
 387  format(/5x,'Z',9X,'PRESS',6x,'DEN',5x,'ProOH1',4x,'ProOH25',4x,
     2 'ProO3D',4x,'PO2D',7x,'PO2',6x,'PH2O',6x,'PO3D',6x,'PH2O2',5x,
     3 'FO3',6x,'FH2O',6x,'NumO2',6x,'NumOH',6x,'NumO1D')
      write(92,388) (z(I),PRESS(I),DEN(I),RateOH1(I),RateOH25(I),
     2 RatePO3(I),AR(23,I),AR(24,I),AR(25,I),AR(26,I),AR(28,I),
     3 USOL(LO3,I),USOL(LH2O,I),SL(LO2,I),SL(LOH,I),SL(LO1D,I),I=1,NZ)
 388  format(2x,1p16e10.3)     

C-AJ 01/24/2022 CALCULATE LIFETIME OF CH4 IN THE SURFACE
      lsCH4_OH = 4.16E-13*((T(1)/298.)**2.18)*EXP(-1230./T(1))*SL(4,1)
      lsCH4_D = 1.3E-10*SL(33,1) + 7.51E-12*SL(33,1)
      lossCH4 = lsCH4_OH + lsCH4_D
      lifetimeCH4 = 1./lossCH4*1./(3600.*24.*365.)
      lifetimeCH4_1 = SL(10,1)*7.5E5/FLOW(LCH4)*1./(3600.*24.*365.)

      print *,'Surface temperature = ', T(1)
      print'(A,1PE12.5)','CH4 flux = ', FLOW(LCH4)
      print'(A,1PE12.5)','number density of OH = ', SL(4,1)
      print'(A,1PE12.5)','number density of O(1D) = ',SL(33,1)
      print'(A,1PE12.5)','loss CH4 due to OH = ', lsCH4_OH
      print'(A,1PE12.5)','loss CH4 due to O(1D) = ',lsCH4_D
      print'(A,1PE12.5)','loss frequency for CH4 = ', lossCH4
      print'(A,1PE12.5)','CH4 surface lifetime (yr) = ', lifetimeCH4
      print'(A,1PE12.5)','CH4 lifetime (yr) = ', lifetimeCH4_1
c-AJ 01/25/2022 ADD CH4 lifetime in the surface into outchem.dat
      write(90,305) lifetimeCH4
 305  format('CH4 lifetime in the surface (yr)',1PE12.5)
      write(90,306) lifetimeCH4_1
 306  format('CH4 lifetime (yr)',1PE12.5)

C-AJ SKIPE THE ISOTOPE CALCULATION 10/26
      GO TO 6000
      PRINT*, ''
      PRINT *, '----------START O ISOTOPE CALCULATION----------'
      PRINT*, ''
      PRINT 631, GPPOXY,GPPCDE
 631  FORMAT(1X,' Gross Primary Productivity FOR O2 CO2= '
     2 ,1PE8.2,2X,1PE8.2,1X,'cm-2s-1')
      IF (COMBIN.EQ.2) THEN
      PRINT*, ' TURN OFF Combinatorial Effect'
      ELSE 
      PRINT*, ' TURN ON Combinatorial Effect'
      END IF

      IF (RECOMB.EQ.1) THEN
      PRINT*, ' TURN OFF Symmetry Effect'
      ELSE 
      PRINT*, ' TURN ON Symmetry Effect'
      END IF
C-PL COMMENT OUT TILL 46 AS DELTA 33=34=36, SO WE DON'T NEED CAL TWICE HERE 
!      DO J=1,NZ
!       DELTAS3(J) = 0.
!       DELTAS4(J) = 0.
!       DELTAS8(J) = 0.
!      ENDDO
!      PHOTOCOEFF91 = 1. 
!      PHOTOCOEFF1 = 1.!/0.9
!      PHOTOCOEFF2 = 1.!/0.95
!      PHOTOCOEFF38 = 1. 
!      PHOTOCOEFF39 = 1.
!      USOLI = 0. 
!      CALL OXYGEN(PHOTOCOEFF2,PHOTOCOEFF38,PHOTOCOEFF39,
!     2 PHOTOCOEFF91,PHOTOCOEFF1,USOLI,DELTAS4,DELTAS3,DELTAS8)

!      DO J=1,NZ

!      DELTAH2S34(J) = (USOLI(LH2SI,J)-USOL(LH2S,J))/USOL(LH2S,J)*1000 !C-PL THIS IS A TEST FOR MDF S34
!      DELTAHS34(J) = (USOLI(LHSI,J)-USOL(LHS,J))/USOL(LHS,J)*1000 
!      DELTASO34(J) = (USOLI(LSIO,J)-USOL(LSO,J))/USOL(LSO,J)*1000 
!      DELTAHSO34(J) = (USOLI(LHSIO,J)-USOL(LHSO,J))/USOL(LHSO,J)*1000 
!      DELTAS34(J) = (USOLI(LSI,J)-USOL(LS,J))/USOL(LS,J)*1000 
!      DELTAS234(J) = (USOLI(LSIS,J)/2.-USOL(LS2,J))/USOL(LS2,J)*1000 
!      DELTASO234(J) = (USOLI(LSIO2,J)-USOL(LSO2,J))/USOL(LSO2,J)*1000 
!      DELTAH2SO434(J) = (USOLI(LH2SIO4,J)-USOL(LH2SO4,J))/
!     2 USOL(LH2SO4,J)*1000

!      ENDDO
!      DO J=1,10 
!       DELTAS334(J) = DELTAS3(J)
!       DELTAS434(J) = DELTAS4(J)
!       DELTAS834(J) = DELTAS8(J)
!       DELTAS3(J) = 0.
!       DELTAS4(J) = 0.
!       DELTAS8(J) = 0.
!      ENDDO
!      PHOTOCOEFF91 = 1.
!      PHOTOCOEFF1 = 1.!/0.9351 
!      PHOTOCOEFF2 = 1.!/0.3  !SO2 photolysis
!      PHOTOCOEFF38 = 1. 
!      PHOTOCOEFF39 = 1. 

      USOLI = 0.
!      CALL SULFUR(PHOTOCOEFF2,PHOTOCOEFF38,PHOTOCOEFF39,
!     2 PHOTOCOEFF91,PHOTOCOEFF1,USOLI,DELTAS4,DELTAS3,DELTAS8)
C      PRINT *, 'HELLO WORLD'
C      PRINT *, 'USOL18=',USOL(18,1) 
C      PRINT *, 'USOL18II=',USOLI(18,1) 
      CALL OXYGEN(USOLI,FLOW) !USOLI IS THE OUT&IN_PUT VAR


      DO J=1,NZ !C-PL DELTA HERE IS FOR DIFFERENT SPECIES IN DIFFERENT LAYER
      DELTAH2CO17(J) = (USOLI(LH2CQ,J)-USOL(LH2CO,J))
     2 /USOL(LH2CO,J)*1000
      DELTAO17(J) = (USOLI(LQ,J)-USOL(LO,J))/USOL(LO,J)*1000
      DELTAH2O17(J) = (USOLI(LH2Q,J)-USOL(LH2O,J))/USOL(LH2O,J)*1000
      DELTAOH17(J) = (USOLI(LQH,J)-USOL(LOH,J))/USOL(LOH,J)*1000
C AJ 01/27 DELETE /2.
C      DELTAHO2A17(J) = (USOLI(LHOQ,J)-USOL(LHO2,J))/USOL(LHO2,J)*1000
C AJ 02/06 DEBUG
C      DELTAHO2A17(J) = (USOLI(LHOQ,J)/2.-USOL(LHO2,J))/USOL(LHO2,J)*1000
      DELTAHO217(J) = (USOLI(LHOQ,J)/2.-USOL(LHO2,J))/USOL(LHO2,J)*1000
C      DELTAHO217(J) = DELTAHO2A17(J)
C  JL 01/23/2020 ADD HQO AJ 01/27 DELETE /2.
C      DELTAHO2B17(J) = (USOLI(LHQO,J)-USOL(LHO2,J))/USOL(LHO2,J)*1000
C      DELTAHO217(J) = (DELTAHO2A17(J)+DELTAHO2B17(J))/2.

      DELTAH2O217(J) = (USOLI(LH2OQ,J)/2.-USOL(LH2O2,J))
     2 /USOL(LH2O2,J)*1000
      DELTAO3A17(J) = (USOLI(LO2Q,J)/2.-USOL(LO3,J))/USOL(LO3,J)*1000
      DELTAO3B17(J) = (USOLI(LOQO,J)-USOL(LO3,J))/USOL(LO3,J)*1000
      DELTAO317(J)  =  2./3.*DELTAO3A17(J)+ 1./3.*DELTAO3B17(J)    !TOTAL DELTA O3
      DELTACO17(J) = (USOLI(LCQ,J)-USOL(LCO,J))/USOL(LCO,J)*1000
      DELTACH3OOHA17(J) = (USOLI(LCH3OQH,J)-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOHB17(J) = (USOLI(LCH3QOH,J)-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOH17(J)=(DELTACH3OOHA17(J)+DELTACH3OOHB17(J))/2. ! TOTAL DELTA CH3OOH
      DELTACH3O2A17(J) = (USOLI(LCH3OQ,J)-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O2B17(J) = (USOLI(LCH3QO,J)-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O217(J)  = (DELTACH3O2A17(J)+DELTACH3O2B17(J))/2. !TOTAL CH3O2
      DELTAN2O17(J) = (USOLI(LN2Q,J)-USOL(LN2O,J))/USOL(LN2O,J)*1000
      DELTANO17(J) = (USOLI(LNQ,J)-USOL(LNO,J))/USOL(LNO,J)*1000
      DELTANO217(J) = (USOLI(LNOQ,J)/2.-USOL(LNO2,J))/USOL(LNO2,J)*1000
      DELTAHNO217(J) = (USOLI(LHNOQ,J)/2.-USOL(LHNO2,J))
     2 /USOL(LHNO2,J)*1000
      DELTAHNO317(J) = (USOLI(LHNO2Q,J)/3.-USOL(LHNO3,J))
     2 /USOL(LHNO3,J)*1000
C AJ 01/27 DELETE /2.
C      DELTAHO2NO2A17(J) = (USOLI(LHOQNO2,J)-USOL(LHO2NO2,J))
C     2 /USOL(LHO2NO2,J)*1000
C AJ 02/06 DEBUG
      DELTAHO2NO2A17(J) = (USOLI(LHOQNO2,J)/2.-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO2B17(J) = (USOLI(LHO2NOQ,J)/2.-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO217(J) = (DELTAHO2NO2A17(J)+DELTAHO2NO2B17(J))/2. !TOTAL HO2NO2
C  JL 01/23/2020 ADD HQONO2  AJ 01/27 DELETE /2.
C      DELTAHO2NO2C17(J) = (USOLI(LHQONO2,J)-USOL(LHO2NO2,J))
C     2 /USOL(LHO2NO2,J)*1000

c      print*,'HQO_A=',DELTAHO2A17(1) 
c      print*,'HQO_B=',DELTAHO2B17(1) 
C  JL 01/23/2020 REDEFINE HO2NO2 from 1/2 to 1/3
C     DELTAHO2NO217(J) = 1./4.*DELTAHO2NO2A17(J)+1./2.*DELTAHO2NO2B17(J)
C     2 + 1./4.*DELTAHO2NO2C17(J) !TOTAL HO2NO2
      DELTANO317(J) = (USOLI(LNO2Q,J)/3.-USOL(LNO3,J))/USOL(LNO3,J)*1000
      DELTAN2O5A17(J) = (USOLI(LN2O4Q,J)/4.-USOL(LN2O5,J))
     2 /USOL(LN2O5,J)*1000
      DELTAN2O5B17(J) = (USOLI(LN2QO4,J)-USOL(LN2O5,J)) 
     2 /USOL(LN2O5,J)*1000
      DELTAN2O517(J)  = 4./5.*DELTAN2O5A17(J) + 1./5.*DELTAN2O5B17(J) !TOTAL N2O5
      DELTAO217(J) = (USOLI(LOQ,J)/2.-USOL(LO2,J))/USOL(LO2,J)*1000
      DELTASO17(J) = (USOLI(LSQ,J)-USOL(LSO,J))/USOL(LSO,J)*1000
      DELTASO217(J) = (USOLI(LSOQ,J)/2.-USOL(LSO2,J))/USOL(LSO2,J)*1000
      DELTAH2SO417(J) = (USOLI(LH2SO3Q,J)/4.-USOL(LH2SO4,J))
     2 /USOL(LH2SO4,J)*1000
      DELTAHSO17(J) = (USOLI(LHSQ,J)-USOL(LHSO,J))/USOL(LHSO,J)*1000
      DELTACO217(J) = (USOLI(LCOQ,J)/2.-USOL(LCO2,J))/USOL(LCO2,J)*1000

      ENDDO 

C-PL COMMENTED OUT AS THEM CAL. ABOVE      
!      DO J=1,10 
!       DELTAS333(J) = DELTAS3(J)
!       DELTAS433(J) = DELTAS4(J)
!       DELTAS833(J) = DELTAS8(J)
!       DELTAS3(J) = 0.
!       DELTAS4(J) = 0.
!       DELTAS8(J) = 0.
!      ENDDO

      print*, ''
      print*, '--------------out of isotope subroutine---------------'
      print*, ''
    
C-PL COMMENTED OUT      
!      print*, 'sulfur mass balance: '
!      sulf_in = H2SVOLC + S2VOLC !H2S and SO2 are volcanically injected
!      sulf_out = 0.
!      SO4LOSS = 0.

      DO 420 I=1,NQI
      SRI(I) = 0.
      DO 430 J=1,20
  430  SRI(I) = SRI(I) + RAINGCI(I,J)*USOLI(I,J)*DEN(J)*1.e5
      PHIDEPI(I) = VDEPI(I)*USOLI(I,1)*DEN(1)
!      sulf_out = sulf_out+SRI(I)+PHIDEPI(I) C-PL COMMENTED OUT
!!      if(I.eq.LSIS) sulf_out = sulf_out+SRI(I)+PHIDEPI(I)  
C      IF(LBOUNDI(I).EQ.0) THEN 
  420  TLOSSI(I) = SRI(I) + PHIDEPI(I)

C-PL RECALCULATE TLOSSI 07/19
      DO I=1,NQI
      IF (LBOUNDI(I).EQ.0) THEN
      TLOSSI(I) = SRI(I) + PHIDEPI(I)
      ELSE IF (LBOUNDI(I).NE.0.AND.FLOWI(I).LT.0) THEN
      TLOSSI(I) = SRI(I) - FLOWI(I) 
      ELSE     
      TLOSSI(I) = SRI(I)
      END IF
      END DO



      DO J=1,NZ
      SO4LOSS = SO4LOSS + CONSO4(J)*DEN(J)*1.e5
      ENDDO
!      TLOSSI(LSIS) = TLOSSI(LSIS) + RATI(95) !C-PL COMMENTOUT 
!      TLOSSI(LSIS2) = TLOSSI(LSIS2) + RATI(96) 
!      TLOSSI(LSIS3) = TLOSSI(LSIS3) + RATI(97)
!      sulf_out = sulf_out + RATI(95) + RATI(96) + RATI(97)
!      sulf_out = sulf_out + SO4LOSS
  
!!     sulf_out = sulf_out + RATI(88) 
!      PRINT*, ' S8PROD      SO4LOSS    | S     in          out'
!      PRINT 440, RATI(88), SO4LOSS,' | ', sulf_in, sulf_out
!  440  format(1P2E12.4,A3,1P2E12.4)
      print*, ''
      print*, 'RAINOUT RATE, PHIDEP, AND LOWER B.C.'
      print*, '   FOLLOWED BY TP, TL, FUP, FLOW, CON'

      print 1800, (ISPECI(I),I=1,10)
      print 1810, 'SR ',(SRI(I),I=1,10)
      print 1810, 'Phi',(PHIDEPI(I),I=1,10)
      print 1911, 'LBC',(LBOUNDI(I),I=1,10)
      print 1810, 'TPI ',(TPI(I),I=1,10)
      print 1810, 'TLI ',(TLI(I),I=1,10)
      print 1810, 'FUP ',(FUPI(I),I=1,10)
      print 1810, 'FLOW',(FLOWI(I),I=1,10)
      print 1810, 'CON ',(CONI(I),I=1,10)
      print*, ''

      print 1800, (ISPECI(I),I=11,20)
      print 1810, 'SR ',(SRI(I),I=11,20)
      print 1810, 'Phi',(PHIDEPI(I),I=11,20)
      print 1911, 'LBC',(LBOUNDI(I),I=11,20)
      print 1810, 'TPI ',(TPI(I),I=11,20)
      print 1810, 'TLI ',(TLI(I),I=11,20)
      print 1810, 'FUP ',(FUPI(I),I=11,20)
      print 1810, 'FLOW',(FLOWI(I),I=11,20)
      print 1810, 'CON ',(CONI(I),I=11,20)
      print*, ''

C AJ 01/24/2020
      print 1800, (ISPECI(I),I=21,30)
      print 1810, 'SR ',(SRI(I),I=21,30)
      print 1810, 'Phi',(PHIDEPI(I),I=21,30)
      print 1911, 'LBC',(LBOUNDI(I),I=21,30)
      print 1810, 'TPI ',(TPI(I),I=21,30)
      print 1810, 'TLI ',(TLI(I),I=21,30)
      print 1810, 'FUP ',(FUPI(I),I=21,30)
      print 1810, 'FLOW',(FLOWI(I),I=21,30)
      print 1810, 'CON ',(CONI(I),I=21,30)
      print*, ''
C AJ 01/24/2020

      print 1800, (ISPECI(I),I=31,NQI)
      print 1810, 'SR ',(SRI(I),I=31,NQI)
      print 1810, 'Phi',(PHIDEPI(I),I=31,NQI)
      print 1911, 'LBC',(LBOUNDI(I),I=31,NQI)
      print 1810, 'TPI ',(TPI(I),I=31,NQI)
      print 1810, 'TLI ',(TLI(I),I=31,NQI)
      print 1810, 'FUP ',(FUPI(I),I=31,NQI)
      print 1810, 'FLOW',(FLOWI(I),I=31,NQI)
      print 1810, 'CON ',(CONI(I),I=31,NQI)
      print*, ''
 1800 format(8x,10(2x,A8,2x))
 1810 format(A5,1P10E12.4)
 1911 FORMAT(A4,7x,10(I1,11X))
!      print*, 'OXYGEN MASS BALANCE: '
!      OXYGEN_in = SGFLUX(LCO) + SGFLUX(LN2O) ! CO AND N2O ARE injected IN
!      OXYGEN_out  = 0.0
!      DO J=1,NQI
!      OXYGEN_out = OXYGEN_out + SRI(J)
!      END DO
!      OXYGEN_out = OXYGEN_out + PSO4AER
!      OXYGEN_CON = OXYGEN_in - OXYGEN_out
!      PRINT *,'  OXYGEN_IN','  OXYGEN_OUT',  '    OXYGEN_CON'
!      PRINT 440, OXYGEN_in, OXYGEN_out,OXYGEN_CON
!  440  format(1P2E12.4,2x,1P2E12.4,2x,1P2E12.4)
!      PRINT *,''
C-SH the *T variables represent approximate MDF trends for identifying MIF signals in 33S(0.515) and 36S(1.89) species
!      DO J=1,NZ
!      DELTAS333T(J) = 0.515*DELTAS334(J)
!      DELTAH2S33T(J) = 0.515*DELTAH2S34(J)
!      DELTAHS33T(J) = 0.515*DELTAHS34(J)
!      DELTASO33T(J) = 0.515*DELTASO34(J)
!      DELTAHSO33T(J) = 0.515*DELTAHSO34(J)
!      DELTAS33T(J) = 0.515*DELTAS34(J)
!      DELTAS233T(J) = 0.515*DELTAS234(J)
!      DELTAS433T(J) = 0.515*DELTAS434(J)
!      DELTASO233T(J) = 0.515*DELTASO234(J)
!      DELTAH2SO433T(J) = 0.515*DELTAH2SO434(J)
!      enddo
      
C-SH calculate approximate CapDel values for printing
!      CDELTAHSO33 = DELTAHSO33(1) - DELTAHSO33T(1)
!      CDELTAS33 = DELTAS33(1) - DELTAS33T(1)
!      DO J=1,10
!       CDELTAS433(J) = DELTAS433(J) - DELTAS433T(J)
!       CDELTAS333(J) = DELTAS333(J) - DELTAS333T(J)
!      ENDDO
!      DO J=1,NZ
!       CDELTAS233(J) = DELTAS233(J) - DELTAS233T(J)
!       CDELTASO33(J) = DELTASO33(J) - DELTASO33T(J)
!       CDELTASO233(J) = DELTASO233(J) - DELTASO233T(J)
!       CDELTAH2S33(J) = DELTAH2S33(J) - DELTAH2S33T(J)
!       CDELTAHS33(J) = DELTAHS33(J) - DELTAHS33T(J)
!       CDELTAH2SO433(J) = DELTAH2SO433(J) - DELTAH2SO433T(J)
!      ENDDO

C-PL WE ONLY CAL. 17O FRACOXYGEN=0.5305(Crockford et al.,2018)      
!      do i=1,2
      fracoxy = 0.5305
!      if(i.eq.2) fracoxy = 1.89
!      DO J=1,10
!      C2DELTAS(i,J)=DELTAS33(J)-1000.*((1.
!     2 + DELTAS34(J)/1000.)**fracsulf - 1.)
!       C2DELTAS3(i,J) = DELTAS333(J) - 1000.*((1.
!     2 + DELTAS334(J)/1000.)**fracsulf - 1.)
!       C2DELTAS4(i,J) = DELTAS433(J) - 1000.*((1.
!     2 + DELTAS434(J)/1000.)**fracsulf - 1.)
!       C2DELTAS8(i,J) = DELTAS833(J) - 1000.*((1.
!     2 + DELTAS834(J)/1000.)**fracsulf - 1.)
!      ENDDO
      DO J=1,NZ
  ! DELTAHSO33(J)=DELTASO34(J)=DELTASO36(J)
!      C2DELTAHSO(i,J)=DELTAHSO33(J)-1000.*((1.  
!     2 + DELTASO34(J)/1000.)**fracoxy - 1.)
!       C2DELTAS2(i,J) = DELTAS233(J) - 1000.*((1.
!     2 + DELTAS234(J)/1000.)**fracoxy - 1.)
!       C2DELTASO(i,J) = DELTASO33(J)  - 1000.*((1.
!     2 + DELTASO34(J)/1000.)**fracoxy - 1.)
!       C2DELTASO2(i,J) = DELTASO233(J) - 1000.*((1.
!     2 + DELTASO234(J)/1000.)**fracoxy - 1.)
!       C2DELTAH2S(i,J) = DELTAH2S33(J) - 1000.*((1.
!     2 + DELTAH2S34(J)/1000.)**fracoxy - 1.)
!       C2DELTAHS(i,J) = DELTAHS33(J) - 1000.*((1.
!     2 + DELTAHS34(J)/1000.)**fracoxy - 1.)
!       C2DELTAH2SO4(i,J) = DELTAH2SO433(J) - 1000.*((1.
!     2 + DELTAH2SO434(J)/1000.)**fracoxy - 1.)

      C2DELTAH2CO(J)=DELTAH2CO17(J)-1000.*((1.  ! delta17O=delta18O
     2 + DELTAH2CO17(J)/1000.)**fracoxy - 1.)
      C2DELTAO(J)=DELTAO17(J)-1000.*((1.
     2 + DELTAO17(J)/1000.)**fracoxy - 1.)
      C2DELTAH2O(J)=DELTAH2O17(J)-1000.*((1.
     2 + DELTAH2O17(J)/1000.)**fracoxy - 1.)
      C2DELTAOH(J)=DELTAOH17(J)-1000.*((1.
     2 + DELTAOH17(J)/1000.)**fracoxy - 1.)
C AJ LUNAR NEW YEAR'S EVE
C      C2DELTAHO2A(J)=DELTAHO2A17(J)-1000.*((1.
C     2 + DELTAHO2A17(J)/1000.)**fracoxy - 1.)
C      C2DELTAHO2B(J)=DELTAHO2B17(J)-1000.*((1.
C     2 + DELTAHO2B17(J)/1000.)**fracoxy - 1.)

C AJ 02/06 DEBUG
      C2DELTAHO2(J)=DELTAHO217(J)-1000.*((1.
     2 + DELTAHO217(J)/1000.)**fracoxy - 1.)

      C2DELTAH2O2(J)=DELTAH2O217(J)-1000.*((1.
     2 + DELTAH2O217(J)/1000.)**fracoxy - 1.)
      C2DELTAO3A(J)=DELTAO3A17(J)-1000.*((1.
     2 + DELTAO3A17(J)/1000.)**fracoxy - 1.)
      C2DELTAO3B(J)=DELTAO3B17(J)-1000.*((1.
     2 + DELTAO3B17(J)/1000.)**fracoxy - 1.)
      C2DELTAO3(J)=DELTAO317(J)-1000.*((1.
     2 + DELTAO317(J)/1000.)**fracoxy - 1.) ! TOTAL O3
      C2DELTACO(J)=DELTACO17(J)-1000.*((1.
     2 + DELTACO17(J)/1000.)**fracoxy - 1.)
      C2DELTACH3OOHA(J)=DELTACH3OOHA17(J)-1000.*((1.
     2 + DELTACH3OOHA17(J)/1000.)**fracoxy - 1.)
      C2DELTACH3OOHB(J)=DELTACH3OOHB17(J)-1000.*((1.
     2 + DELTACH3OOHB17(J)/1000.)**fracoxy - 1.)
      C2DELTACH3OOH(J)=DELTACH3OOH17(J)-1000.*((1.
     2 + DELTACH3OOH17(J)/1000.)**fracoxy - 1.) !TOTAL CH3OOH
      C2DELTACH3O2A(J)=DELTACH3O2A17(J)-1000.*((1.
     2 + DELTACH3O2A17(J)/1000.)**fracoxy - 1.)
      C2DELTACH3O2B(J)=DELTACH3O2B17(J)-1000.*((1.
     2 + DELTACH3O2B17(J)/1000.)**fracoxy - 1.)
      C2DELTACH3O2(J)=DELTACH3O217(J)-1000.*((1.
     2 + DELTACH3O217(J)/1000.)**fracoxy - 1.)  !TOTAL CH3O2
      C2DELTAN2O(J)=DELTAN2O17(J)-1000.*((1.
     2 + DELTAN2O17(J)/1000.)**fracoxy - 1.)
      C2DELTANO(J)=DELTANO17(J)-1000.*((1.
     2 + DELTANO17(J)/1000.)**fracoxy - 1.)
      C2DELTANO2(J)=DELTANO217(J)-1000.*((1.
     2 + DELTANO217(J)/1000.)**fracoxy - 1.)
      C2DELTAHNO2(J)=DELTAHNO217(J)-1000.*((1.
     2 + DELTAHNO217(J)/1000.)**fracoxy - 1.)
      C2DELTAHNO3(J)=DELTAHNO317(J)-1000.*((1.
     2 + DELTAHNO317(J)/1000.)**fracoxy - 1.)
      C2DELTAHO2NO2A(J)=DELTAHO2NO2A17(J)-1000.*((1.
     2 + DELTAHO2NO2A17(J)/1000.)**fracoxy - 1.)
      C2DELTAHO2NO2B(J)=DELTAHO2NO2B17(J)-1000.*((1.
     2 + DELTAHO2NO2B17(J)/1000.)**fracoxy - 1.)
C AJ LUNAR NEW YEAR'S EVE
C      C2DELTAHO2NO2C(J)=DELTAHO2NO2C17(J)-1000.*((1.
C     2 + DELTAHO2NO2C17(J)/1000.)**fracoxy - 1.)

      C2DELTAHO2NO2(J)=DELTAHO2NO217(J)-1000.*((1.
     2 + DELTAHO2NO217(J)/1000.)**fracoxy - 1.)  !TOTAL HO2NO2
      C2DELTANO3(J)=DELTANO317(J)-1000.*((1.
     2 + DELTANO317(J)/1000.)**fracoxy - 1.)
      C2DELTAN2O5A(J)=DELTAN2O5A17(J)-1000.*((1.
     2 + DELTAN2O5A17(J)/1000.)**fracoxy - 1.)
      C2DELTAN2O5B(J)=DELTAN2O5B17(J)-1000.*((1.
     2 + DELTAN2O5B17(J)/1000.)**fracoxy - 1.)
      C2DELTAN2O5(J)=DELTAN2O517(J)-1000.*((1.
     2 + DELTAN2O517(J)/1000.)**fracoxy - 1.)  !TOATAL N2O5
      C2DELTAO2(J)=DELTAO217(J)-1000.*((1.
     2 + DELTAO217(J)/1000.)**fracoxy - 1.)
      C2DELTASO(J)=DELTASO17(J)-1000.*((1.
     2 + DELTASO17(J)/1000.)**fracoxy - 1.)
      C2DELTASO2(J)=DELTASO217(J)-1000.*((1.
     2 + DELTASO217(J)/1000.)**fracoxy - 1.)
      C2DELTAH2SO4(J)=DELTAH2SO417(J)-1000.*((1.
     2 + DELTAH2SO417(J)/1000.)**fracoxy - 1.)
      C2DELTAHSO(J)=DELTAHSO17(J)-1000.*((1.
     2 + DELTAHSO17(J)/1000.)**fracoxy - 1.)
      C2DELTACO2(J)=DELTACO217(J)-1000.*((1.
     2 + DELTACO217(J)/1000.)**fracoxy - 1.)
      ENDDO
      
      print*, 'delta 17O values with altitude'
      PRINT 102, 'Z(km)',(ISPECI(I),I=1,10)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
      PRINT 104,J-0.5,
C     2 DELTAH2CO17(J),DELTAO17(J),DELTAH2O17(J),
C     3 DELTAOH17(J),DELTAHO2A17(J),DELTAH2O217(J),DELTAO3A17(J),
C     4 DELTAO3B17(J),DELTACO17(J),DELTACH3OOHA17(J)
     2 DELTAH2CO17(J),DELTAO17(J),DELTAH2O17(J),
     3 DELTAOH17(J),DELTAHO217(J),DELTAH2O217(J),DELTAO3A17(J),
     4 DELTAO3B17(J),DELTACO17(J),DELTACH3OOHA17(J)
      ENDDO
      print*, ''
      PRINT 102, 'Z(km)',(ISPECI(I),I=11,20)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 104,J-0.5,
     1 DELTACH3OOHB17(J),DELTACH3O2A17(J),DELTACH3O2B17(J),
     2 DELTAN2O17(J),DELTANO17(J),DELTANO217(J),DELTAHNO217(J),
     3 DELTAHNO317(J),DELTAHO2NO2A17(J),DELTAHO2NO2B17(J)
      ENDDO
      print*, ''
C AJ 01/24/2020
      PRINT 102, 'Z(km)',(ISPECI(I),I=21,30)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 104,J-0.5,
     1 DELTANO317(J),DELTAN2O5A17(J),DELTAN2O5B17(J),DELTAO217(J),
     2 DELTASO17(J),DELTASO217(J),DELTAH2SO417(J),DELTAHSO17(J),
     3 DELTACO217(J) !,DELTAHO2B17(J)
      ENDDO
      print*, ''
C AJ
c      PRINT 102, 'Z(km)',(ISPECI(I),I=31,NQI)
c      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
c       PRINT 104,J-0.5,
c     1 DELTAHO2NO2C17(J) 
c      ENDDO
      print*, ''
      PRINT 102, 'Z(km)','O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I'
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 104,J-0.5,
     1 DELTAO317(J), DELTACH3OOH17(J) ,DELTACH3O217(J),
     2 DELTAHO2NO217(J),  DELTAN2O517(J)
      ENDDO
  102 format(2x,A8,2x,10(A8,2x))
  104 format(2x,F6.1,1x,1P10E10.2)
      
!      print*, '         delta-34S delta-33S'
!      print*, 'species    value     value'
!      PRINT 102,'  H2S  ',DELTAH2S34(1),DELTAH2S33(1)
!      PRINT 102,'  HS   ', DELTAHS34(1),DELTAHS33(1) 
!      PRINT 102,'  SO   ', DELTASO34(1),DELTASO33(1) 
!      PRINT 102,'  SO2  ', DELTASO234(1),DELTASO233(1)
!      PRINT 102,'  H2SO4', DELTAH2SO434(1),DELTAH2SO433(1) 
!      PRINT 102,'  S    ', DELTAS34(1),DELTAS33(1) 
!      PRINT 102,'  S2   ', DELTAS234(1),DELTAS233(1) 
!      PRINT 102,'  S3   ', DELTAS334(1),DELTAS333(1) 
!      PRINT 102,'  S4   ', DELTAS434(1),DELTAS433(1) 
!      PRINT 102,'  S8   ', DELTAS834(1),DELTAS833(1) 
!      PRINT 102,'  HSO  ', DELTAHSO34(1),DELTAHSO33(1)
!  102 format(A8,1P2E10.2)

!      PRINT *,'CDELTAH2S33=', CDELTAH2S33(1), 'LOSS H2S=', TLOSS(LH2S) 
!      PRINT *,'CDELTAHS33=', CDELTAHS33(1), 'LOSS HS=', TLOSS(LHS) 
!      PRINT *,'CDELTAS33=', CDELTAS33, 'LOSS S=', TLOSS(LS) 
!      PRINT *,'CDELTAS233=', CDELTAS233(1), 'LOSS S2=', TLOSS(LS2)*2 
!      PRINT *,'CDELTAS433=', CDELTAS433(1) 
!      PRINT *,'CDELTAS333=', CDELTAS333(1) 
!      PRINT *,'CDELTASO33=', CDELTASO33(1), 'LOSS SO=', TLOSS(LSO) 
!      PRINT *,'CDELTAHSO33=', CDELTAHSO33, 'LOSS HSO=', TLOSS(LHSO) 
!      PRINT *,'CDELTASO233=', CDELTASO233(1), 'LOSS SO2=', TLOSS(LSO2) 
!      PRINT *,'CDELTAH2SO433=', CDELTAH2SO433(1), 'LOSS H2SO4=', SO4LOS
      
!      print*, ''
!      print*, '         precise capdel formula values: '
!      print*, 'CapDel(spec) =   33S       36S       36S/33S'
!      PRINT 1030,'CapDel(H2S) = ', (C2DELTAH2S(i,1),i=1,2),
!     2 C2DELTAH2S(2,1)/C2DELTAH2S(1,1)
!      PRINT 1030,'CapDel(HS) = ', (C2DELTAHS(i,1),i=1,2),
!     2 C2DELTAHS(2,1)/C2DELTAHS(1,1)
!      PRINT 1030,'CapDel(S) = ', (C2DELTAS(i,1),i=1,2),
!     2 C2DELTAS(2,1)/C2DELTAS(1,1)
!      PRINT 1030,'CapDel(S2) = ', (C2DELTAS2(i,1),i=1,2),
!     2 C2DELTAS2(2,1)/C2DELTAS2(1,1)
!      PRINT 1030,'CapDel(S3) = ', (C2DELTAS3(i,1),i=1,2),
!     2 C2DELTAS3(2,1)/C2DELTAS3(1,1)
!      PRINT 1030,'CapDel(S4) = ', (C2DELTAS4(i,1),i=1,2),
!     2 C2DELTAS4(2,1)/C2DELTAS4(1,1)
!      PRINT 1030,'CapDel(S8) = ', (C2DELTAS8(i,1),i=1,2),
!     2 C2DELTAS8(2,1)/C2DELTAS8(1,1)
!      PRINT 1030,'CapDel(SO) = ', (C2DELTASO(i,1),i=1,2),
!     2 C2DELTASO(2,1)/C2DELTASO(1,1)
!      PRINT 1030,'CapDel(HSO) = ', (C2DELTAHSO(i,1),i=1,2),
!     2 C2DELTAHSO(2,1)/C2DELTAHSO(1,1)
!      PRINT 1030,'CapDel(SO2) = ', (C2DELTASO2(i,1),i=1,2),
!     2 C2DELTASO2(2,1)/C2DELTASO2(1,1)
!      PRINT 1030,'CapDel(H2SO4) = ', (C2DELTAH2SO4(i,1),i=1,2),
!     2 C2DELTAH2SO4(2,1)/C2DELTAH2SO4(1,1)
! 1030 format(A16,1P3E10.2)
      
C-SH the following implementation is only APPROXIMATE - does not take into 
! account varying delta profiles with altitude.
!      print*,''
!      PRINT *, 'TOTAL ISOTOPIC BALANCE'
!      do i=1,2
!      label1 = ' 33S '
!      if(i.eq.2) label1 = ' 36S '
!      if(i.eq.2) print*, ''
!      TOTALREDS = C2DELTAH2S(i,1)*TLOSSI(1)+C2DELTAHS(i,1)*
!     2 TLOSSI(2)+C2DELTAS(i,1)*TLOSSI(3) + C2DELTAS2(i,1)*
!     3 TLOSS(8) +C2DELTASO(i,1)*TLOSS(4) + C2DELTAHSO(i,1)*
!     4 TLOSS(7) +C2DELTAS8(i,1)*RATI(88)
!      TOTALOXIS = C2DELTASO2(i,1)*TLOSS(5)
!     2 + C2DELTAH2SO4(i,1)*SO4LOSS + C2DELTAH2SO4(i,1)*TLOSS(6)
!      PRINT 1031, 'TOTALRED',label1,' = ', TOTALREDS
!      PRINT 1031, 'TOTALOXI',label1,' = ', TOTALOXIS 
!      enddo
!      print 1035, 'H2S','HS','S','SIS','SO','HSO','S8','SO2',
!     2 'SO4(aer)','H2SO4'
!      print 1034, 33,C2DELTAH2S(1,1)*TLOSSI(1),C2DELTAHS(1,1)*
!     2 TLOSSI(2),C2DELTAS(1,1)*TLOSSI(3),C2DELTAS2(1,1)*
!     3 TLOSS(8),C2DELTASO(1,1)*TLOSS(4), C2DELTAHSO(1,1)*
!     4 TLOSS(7),C2DELTAS8(1,1)*RATI(88),C2DELTASO2(1,1)*TLOSS(5),
!     5 C2DELTAH2SO4(1,1)*SO4LOSS,C2DELTAH2SO4(1,1)*TLOSS(6)
!      print 1034, 36,C2DELTAH2S(2,1)*TLOSSI(1),C2DELTAHS(2,1)*
!     2 TLOSSI(2),C2DELTAS(2,1)*TLOSSI(3),C2DELTAS2(2,1)*
!     3 TLOSS(8),C2DELTASO(2,1)*TLOSS(4),C2DELTAHSO(2,1)*
!     4 TLOSS(7),C2DELTAS8(2,1)*RATI(88),C2DELTASO2(2,1)*TLOSS(5),
!     5 C2DELTAH2SO4(2,1)*SO4LOSS,C2DELTAH2SO4(2,1)*TLOSS(6)

      goto 2201
C-PL  delete SEVERAL LINES HERE SINCE WE DID'T USE THEM
 2201 continue
    
 1031 format(A8,A5,A3,1PE10.2)
 1034 format(I2,1x,1P12E10.2)
 1035 format(3x,12(A8,2x))
 
      print*, ''
!      DO K=1,2
      label1=' 17O'
!      if(K.eq.2) label1='36S''
      print*, '-------------',label1,' budget ----------'
!      DO 4200 I=1,NQI
! 4200 SRI(I) = 0.
      SO4LOSS1 = 0.
      S8LOSS = 0.
      del_elS = 0.; TLOSS_elS = 0.; TLOSSI_elS = 0.
      elS_cap33 = 0.; elS_cap36 = 0.; elS_frac = 0.
      oxy_in = 0.
      oxy_out = 0.
C-SH I've commandeered the Phidep terms to calculate fractionation balance terms.
!      PHIDEPI(1) = (TLOSSI(1)/TLOSS(LH2S) - 1.)*1000. !DELTA H2S
!      PHIDEPI(2) = (TLOSSI(2)/TLOSS(LHS) - 1.)*1000.  !DELTA HS
!      PHIDEPI(3) = (TLOSSI(3)/TLOSS(LS) - 1.)*1000.   !DELTA S
!      PHIDEPI(4) = (TLOSSI(4)/TLOSS(LSO) - 1.)*1000.  !DELTA SO
!      PHIDEPI(5) = (TLOSSI(5)/TLOSS(LSO2) - 1.)*1000.  !DELTA SO2
!      PHIDEPI(6) = (TLOSSI(6)/TLOSS(LH2SO4) - 1.)*1000.  !DELTA H2SO4
!      PHIDEPI(7) = (TLOSSI(7)/TLOSS(LHSO) - 1.)*1000.  !DELTA HSO
!      PHIDEPI(8) = (TLOSSI(8)/(2.*TLOSS(LS2)) - 1.)*1000.  !DELTA S2
!      PHIDEPI(9) = (TLOSSI(9)/(3.*TLOSS(LS3)) - 1.)*1000.  !DELTA S3
!      PHIDEPI(10) = (TLOSSI(10)/(4.*TLOSS(LS4)) - 1.)*1000.  !DELTA S4
!C-PL DELTA HERE IS THE DELTA VALUE TO THE SEDIMENT RAINOUT+VDEP
      DELTA17O(1)   = (TLOSSI(1)/TLOSS(LH2CO) - 1.)*1000.
      DELTA17O(2)   = (TLOSSI(2)/TLOSS(LO) - 1.)*1000.
      DELTA17O(3)   = 0.0!(TLOSSI(3)/TLOSS(LH2O) - 1.)*1000.
      DELTA17O(4)   = (TLOSSI(4)/TLOSS(LOH) - 1.)*1000.
C JL ON LUNAR NEW YEAR EVE 2020 AJ 01/28 DELETE /2.
      DELTA17O(5)   = (TLOSSI(5)/TLOSS(LHO2) - 1.)*1000.
      DELTA17O(30)   = (TLOSSI(30)/TLOSS(LHO2) - 1.)*1000.
      DELTAHO2 = (DELTA17O(5)+DELTA17O(30))/2. ! TOTAL DELTA HO2

      DELTA17O(6)   = (TLOSSI(6)/(2.*TLOSS(LH2O2))- 1.)*1000.
      DELTA17O(7)   = (TLOSSI(7)/(2.*TLOSS(LO3)) - 1.)*1000. !O2Q
      DELTA17O(8)   = (TLOSSI(8)/TLOSS(LO3) - 1.)*1000. !OQO
      DELTAO3    = 2./3.*DELTA17O(7)+ 1./3.*DELTA17O(8)     !TOTAL DELTA O3
      DELTA17O(9)   = (TLOSSI(9)/TLOSS(LCO) - 1.)*1000.
      DELTA17O(10)   = (TLOSSI(10)/TLOSS(LCH3OOH) - 1.)*1000. !CH3OQH
      DELTA17O(11)   = (TLOSSI(11)/TLOSS(LCH3OOH) - 1.)*1000. !CH3QOH
      DELTACH3OOH = (DELTA17O(10)+DELTA17O(11))/2. ! TOTAL DELTA CH3OOH
      DELTA17O(12)   = (TLOSSI(12)/TLOSS(LCH3O2) - 1.)*1000. !CH3OQ
      DELTA17O(13)   = (TLOSSI(13)/TLOSS(LCH3O2) - 1.)*1000. !CH3QO
      DELTACH3O2  = (DELTA17O(12)+DELTA17O(13))/2. !TOTAL CH3O2
      DELTA17O(14)   = (TLOSSI(14)/TLOSS(LN2O) - 1.)*1000.
      DELTA17O(15)   = (TLOSSI(15)/TLOSS(LNO) - 1.)*1000.
      DELTA17O(16)   = (TLOSSI(16)/(2.*TLOSS(LNO2)) - 1.)*1000.
      DELTA17O(17)   = (TLOSSI(17)/(2.*TLOSS(LHNO2)) - 1.)*1000.
      DELTA17O(18)   = (TLOSSI(18)/(3.*TLOSS(LHNO3)) - 1.)*1000.
      DELTA17O(19)   = (TLOSSI(19)/TLOSS(LHO2NO2) - 1.)*1000. !HOQNO2
      DELTA17O(20)   = (TLOSSI(20)/(2.*TLOSS(LHO2NO2)) - 1.)*1000. !HO2NOQ
      DELTA17O(31)   = (TLOSSI(31)/TLOSS(LHO2NO2) - 1.)*1000. !HQONO2
C AJ 01/28 CHANGE 1./3 TO 1/4. 1/2. 1/4.
      DELTAHO2NO2 = 1./4.*DELTA17O(19)+ 1./2.*DELTA17O(20)
     2 + 1./4.*DELTA17O(31)                                  !TOTAL HO2NO2
      DELTA17O(21)   = (TLOSSI(21)/(3.*TLOSS(LNO3)) - 1.)*1000.
      DELTA17O(22)   = (TLOSSI(22)/(4.*TLOSS(LN2O5)) - 1.)*1000. !N2O4Q
      DELTA17O(23)   = (TLOSSI(23)/TLOSS(LN2O5) - 1.)*1000. !N2QO4
      DELTAN2O5   = 4./5.*DELTA17O(22) + 1./5.*DELTA17O(23) !TOTAL N2O5
      DELTA17O(24)   = (TLOSSI(24)/(2.*TLOSS(LO2)) - 1.)*1000.
!      DELTA17O(24)   = (SRI(24)/(2.*SR(LO2)) - 1.)*1000.
      DELTA17O(25)   = (TLOSSI(25)/TLOSS(LSO) - 1.)*1000.
      DELTA17O(26)   = (TLOSSI(26)/(2.*TLOSS(LSO2)) - 1.)*1000.
      DELTA17O(27)   = (TLOSSI(27)/(4.*TLOSS(LH2SO4)) - 1.)*1000.
      DELTA17O(28)   = (TLOSSI(28)/TLOSS(LHSO) - 1.)*1000.
!      DELTA17O(29)   = (TLOSSI(29)/(2.*TLOSS(LCO2)) - 1.)*1000.
      DELTA17O(29)   = (SRI(29)/(2.*SR(LCO2)) - 1.)*1000. 

C-PL COMPUTE CONSERVATION OF OXYGEN  07/2019 
      OXYDEPI = 0.
      OXYRANI = 0. 
      OXYUPI  = 0.
      OXYLOSI = 0.
C AJ 01/28 ADD (LHQO) +(LHQONO2)
      OXYDEPI = - (FLOWI(LH2CQ) + FLOWI(LQ) + FLOWI(LQH) + FLOWI(LHOQ)
     2  + FLOWI(LH2OQ) + FLOWI(LO2Q) + FLOWI(LOQO) + FLOWI(LCH3OQH)
     3  + FLOWI(LCH3QOH) + FLOWI(LCH3OQ) + FLOWI(LCH3QO) + FLOWI(LNQ) 
     4  + FLOWI(LNOQ) + FLOWI(LHNOQ)+ FLOWI(LHNO2Q) + FLOWI(LHOQNO2) 
     5  + FLOWI(LHO2NOQ) + FLOWI(LNO2Q) + FLOWI(LN2O4Q)
     6  + FLOWI(LN2QO4) + FLOWI(LSQ) +FLOWI(LSOQ) 
     7  + FLOWI(LH2SO3Q) + FLOWI(LHSQ) +FLOWI(LHQO) +FLOWI(LHQONO2))!+ FLOWI(LCOQ) ) !
!      IF (LBOUND(LCQ).EQ.0) OXYDEPI = OXYDEPI - FLOWI(LCO)
!      IF (LBOUND(LN2Q).EQ.0) OXYDEPI = OXYDEPI - FLOWI(LN2O)
      OXYRANI = SRI(LH2CQ) + SRI(LQ) + SRI(LQH) + SRI(LHOQ)
     2  + SRI(LH2OQ) + SRI(LO2Q) + SRI(LOQO) + SRI(LCH3OQH)
     3  + SRI(LCH3QOH) + SRI(LCH3OQ) + SRI(LCH3QO) + SRI(LNQ) 
     4  + SRI(LNOQ) + SRI(LHNOQ)+ SRI(LHNO2Q) + SRI(LHOQNO2) 
     5  + SRI(LHO2NOQ) + SRI(LNO2Q) + SRI(LN2O4Q)
     6  + SRI(LN2QO4) + SRI(LSQ) +SRI(LSOQ) 
     7  + SRI(LH2SO3Q) + SRI(LHSQ) + PSO4AER
     8  + SRI(LCQ) + SRI(LN2Q) + SRI(LOQ) + SRI(LCOQ)
     9  + SRI(LHQO) + SRI(LHQONO2) 
       OXYUPI = FUPI(LH2CQ) + FUPI(LQ) + FUPI(LQH) + FUPI(LHOQ)
     2  + FUPI(LH2OQ) + FUPI(LO2Q) + FUPI(LOQO) + FUPI(LCH3OQH)
     3  + FUPI(LCH3QOH) + FUPI(LCH3OQ) + FUPI(LCH3QO) + FUPI(LNQ) 
     4  + FUPI(LNOQ) + FUPI(LHNOQ)+ FUPI(LHNO2Q) + FUPI(LHOQNO2) 
     5  + FUPI(LHO2NOQ) + FUPI(LNO2Q) + FUPI(LN2O4Q)
     6  + FUPI(LN2QO4) + FUPI(LSQ) +FUPI(LSOQ) 
     7  + FUPI(LH2SO3Q) + FUPI(LHSQ) 
     8  + FUPI(LCQ) + FUPI(LN2Q) + FUPI(LOQ) + FUPI(LCOQ)
     9  + FUPI(LHQO) + FUPI(LHQONO2)
      OXYLOSI = OXYDEPI + OXYRANI + OXYUPI + TPI(LH2Q)-TLI(LH2Q)
 
      OXYPROI = 0.
      OXYPROI = FLOWI(LCQ) + FLOWI(LN2Q)!SGFLUX(LCO) + SGFLUX(LN2O)
      OXYPROI = OXYPROI + 2.*GPPOXY +  FLOWI(LOQ)
     2                  + 2.*GPPCDE +  FLOWI(LCOQ)

      CONOXYI = OXYLOSI - OXYPROI
      print 177, OXYLOSI,OXYPROI,CONOXYI,
     2 CONOXYI/MIN(OXYLOSI,OXYPROI) *100. 
 177  FORMAT(/1X,'CONSERVATION OF OXYGENI:',/5X,'OXYLOSI =',1PE10.3,
     2  2X,'OXYPROI =',E10.3,2X,'CONOXYI =',E10.3,
     3  3X,'ERRI = ',E10.3,' %' )
      print 178, OXYDEPI,OXYRANI,OXYUPI,TPI(LH2Q)-TLI(LH2Q)
 178  FORMAT(/5X,'OXYDEPI =',1PE10.3,2X,'OXYRANI =',E10.3,
     2 2X,'OXYUPI = ',E10.3,2X,'TPI-TLI(H2Q) =',2E10.3)

      PRINT 189, FLOWI(LCQ),FLOWI(LN2Q),2.*GPPOXY+FLOWI(LOQ)
     2 ,2*GPPOXY,2.*GPPCDE+FLOWI(LCOQ),2*GPPCDE
 189  FORMAT(5X,'FLOW(LCQ) =',1PE10.3,2X,'FLOW(LN2Q) =',E10.3,
     2 /5X,'GPP-FLOW(LOQ)=',E10.3,2X,'GPPO2=',E10.3,
     2 2X,'GPP-FLOW(LCOQ)=',E10.3,2X,'GPPCO2=',E10.3)

      print 187, FLOW(LO2)*2.,FLOWI(LOQ)
 187  format(/5X,'2*FLOW(O2) =',1PE12.5,2X,'FLOW(OQ) =',1PE12.5)
      print 188, FLOW(LO2)/SL(LO2,1), FLOW(LCO2)/SL(LCO2,1)
 188  format(5X,'Vdep O2/CO2 from Main code =',1P2E18.8)

      print 289, TPH2OI,TLH2OI,TPH2OI-TLH2OI,TPI(LH2Q)-TPH2OI,
     2 TLI(LH2Q)-TLH2OI,TPI(LH2Q)-TPH2OI-(TLI(LH2Q)-TLH2OI),
     3 TPI(LH2Q),TLI(LH2Q),TPI(LH2Q)-TLI(LH2Q)
 289  FORMAT(/5X,' H2Q     TPI       TLI       TPI-TLI ',
     2       /5X, 'TROPOS',1P3E10.3,
     2       /5X, 'STRATO',1P3E10.3,
     3       /5X, 'TOTAL ',1P3E10.3)
c      DO 4200 I=1,NQI
c 4200 oxy_out = oxy_out + TLOSS(I)*DELTA17O(I)
c      print*, '      weighted delta removal terms'
c      PRINT*, ' oxygen   delta*TLOSS'
c      PRINT 4440,  oxy_out
c 4440 format(10x,1PE10.2)
C-PL COMMENT OUT SINCE THIS IS ABOUT SULFUR
!      DO 4400 I=1,NQI   !C-PL CHANGE LATER MASS BALANCE
! 4400 sulf_out = sulf_out+PHIDEPI(I)*TLOSS(I)
!      DO J=1,NZ
!      S8LOSS = (RATI(88)/(8.*TP(LS8AER)) - 1.)*1000. !DELTA S8
!      SO4LOSS1 = (SO4LOSS/TP(LSO4AER) - 1.)*1000. !DELTA H2SO4 condensation
!      ENDDO
!      sulf_out = sulf_out + SO4LOSS1*SO4LOSS + S8LOSS*RATI(88)
!      SOUT = S2VOLC + H2SVOLC
!      TLOSS_elS = TLOSS(LS)+TLOSS(LS2)*2.+TLOSS(LS3)*3.+TLOSS(LS4)*4.
!     2 +TP(LS8AER)*8.
!      TLOSSI_elS = TLOSSI(3)+TLOSSI(8)+TLOSSI(9)+TLOSSI(10)+RATI(88)
!      del_elS = (TLOSS(LS)*PHIDEPI(3)+TLOSS(LS2)*2.*PHIDEPI(8)+
!     2 TLOSS(LS3)*3.*PHIDEPI(9)+TLOSS(LS4)*4.*PHIDEPI(10)+
!     3 TP(LS8AER)*8.*S8LOSS)/TLOSS_elS  !TOTAL DELTA
!      elS_cap33 = del_elS - 1000.*((1.+del_elS/1000.)**0.515-1.)
!      elS_cap36 = del_elS - 1000.*((1.+del_elS/1000.)**1.89-1.)

C-PL CALCULATE CAP DEL 
      DO I=1,NQI
!      CAP33(I) = PHIDEPI(I) - 1000.*((1.+PHIDEPI(I)/1000.)**0.515-1.)!0.515*PHIDEPI(I)!
!      CAP36(I) = PHIDEPI(I) - 1000.*((1.+PHIDEPI(I)/1000.)**1.89-1.)!1.89*PHIDEPI(I)!
      CAP17(I) = DELTA17O(I)-1000.*((1.+DELTA17O(I)/1000.)**0.5305-1.)
      ENDDO
      CAP17(LH2Q) = 0.0
C-PL  TOTAL CAP DEL NUMBERS FOR SPECIES THAT HAVE 2 ISOTOPIC SPECIES
      CAPO317   = DELTAO3 -1000.*((1.+DELTAO3/1000.)**0.5305-1.) 
      CAPCH3OOH17=DELTACH3OOH-1000.*((1.+DELTACH3OOH/1000.)**0.5305-1.)
      CAPC3O217 = DELTACH3O2 - 1000.*((1.+DELTACH3O2/1000.)**0.5305-1.)
      CAPHO2NO217=DELTAHO2NO2-1000.*((1.+DELTAHO2NO2/1000.)**0.5305-1.)
      CAPN2O517  =DELTAN2O5-1000.*((1.+DELTAN2O5/1000.)**0.5305-1.)   
C JL ON THE LUNAR NEW YEAR EVE 2020 
      CAPHO217  = DELTAHO2-1000.*((1.+DELTAHO2/1000.)**0.5305-1.)  
   
!      CAP33SO4 = SO4LOSS1 - 1000.*((1.+SO4LOSS1/1000.)**0.515-1.)!0.515*SO4LOSS1!
!      CAP36SO4 = SO4LOSS1 - 1000.*((1.+SO4LOSS1/1000.)**1.89-1.)!1.89*SO4LOSS1!
!      CAP33S8 = S8LOSS - 1000.*((1.+S8LOSS/1000.)**0.515-1.)!0.515*S8LOSS!
!      CAP36S8 = S8LOSS - 1000.*((1.+S8LOSS/1000.)**1.89-1.)!1.89*S8LOSS!
      
C-PL COMMENT OUT SINCE THIS IS ABOUT SULFUR    
!      print*, '      weighted delta removal terms'
!      PRINT*, ' S8        SO4      | S   in        out'
!      PRINT 4440, S8LOSS, SO4LOSS1,' | ', sulf_in, sulf_out
! 4440 format(1P2E10.2,A3,1P2E10.2)

!      print*, ''
!      print 1801, (ISPECI(I),I=1,NQI),'SO4con','S8aer','S_n'
!      print 1811, 'delta      ',(PHIDEPI(I),I=1,NQI),SO4LOSS1,S8LOSS
!     2 ,del_elS
!      print 1811, 'CAP-33S    ',(CAP33(I),I=1,NQI),CAP33SO4,CAP33S8
!     2 ,elS_cap33
!!      print 1811, 'CAP33/del',((PHIDEPI(I) - 1000.*((1.+PHIDEPI(I)/
!!     2 1000.)**0.515-1.))/PHIDEPI(I),I=1,NQI),(SO4LOSS1 - 1000.*
!!     3 ((1.+SO4LOSS1/1000.)**0.515-1.))/SO4LOSS1,(S8LOSS - 1000.*
!!     4 ((1.+S8LOSS/1000.)**0.515-1.))/S8LOSS
!      print 1811, 'CAP-36S    ',(CAP36(I),I=1,NQI),CAP36SO4,CAP36S8
!     2 ,elS_cap36
!      print 1811, '36S/33S    ',(CAP36(I)/CAP33(I),I=1,NQI),
!     2 CAP36SO4/CAP33SO4,CAP36S8/CAP33S8,elS_cap36/elS_cap33

      print*, ''
      print 1801, (ISPECI(I),I=1,10)
      print 1811, 'delta      ',(DELTA17O(I),I=1,10)
      print 1811, 'CAP-17O     ',(CAP17(I),I=1,10)

      print*, ''
      print 1811, 'TLOSSI     ',(TLOSSI(I),I=1,10)
      print 1811, 'TLOSS      ',TLOSS(LH2CO),TLOSS(LO),TLOSS(LH2O),
     2 TLOSS(LOH),TLOSS(LHO2),TLOSS(LH2O2),TLOSS(LO3),TLOSS(LO3),
     3 TLOSS(LCO),TLOSS(LCH3OOH)

      print*, ''
      print 1801, (ISPECI(I),I=11,20)
      print 1811, 'delta      ',(DELTA17O(I),I=11,20)
      print 1811, 'CAP-17O    ',(CAP17(I),I=11,20)

      print*, ''
      print 1811, 'TLOSSI     ',(TLOSSI(I),I=11,20)
      print 1811, 'TLOSS      ',TLOSS(LCH3OOH),TLOSS(LCH3O2),
     2 TLOSS(LCH3O2),TLOSS(LN2O),TLOSS(LNO),TLOSS(LNO2),
     3 TLOSS(LHNO2),TLOSS(LHNO3),TLOSS(LHO2NO2),TLOSS(LHO2NO2)
C AJ 01/24/2020
      print*, ''
      print 1801, (ISPECI(I),I=21,30)
      print 1811, 'delta      ',(DELTA17O(I),I=21,30)
      print 1811, 'CAP-17O    ',(CAP17(I),I=21,30)

      print*, ''
      print 1811, 'TLOSSI     ',(TLOSSI(I),I=21,30)
      print 1811, 'TLOSS      ',TLOSS(LNO3),TLOSS(LN2O5),
     2 TLOSS(LN2O5),TLOSS(LO2),TLOSS(LSO),TLOSS(LSO2),
     3 TLOSS(LH2SO4),TLOSS(LHSO),TLOSS(LCO2),TLOSS(LHO2)
C AJ
      print*, ''
      print 1801, (ISPECI(I),I=31,NQI)
      print 1811, 'delta      ',(DELTA17O(I),I=31,NQI)
      print 1811, 'CAP-17O    ',(CAP17(I),I=31,NQI)

      print*, ''
      print 1811, 'TLOSSI     ',(TLOSSI(I),I=31,NQI)
      print 1811, 'TLOSS      ',TLOSS(LHO2NO2)
C-PL SPECIES HAVE TWO ISOTOPIC SPECIES
      print*, ''
      print 1802, 'O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I','HO2I'
 1802 format(18x,6(A10))

      print 1803, 'delta      ',DELTAO3,DELTACH3OOH,
     2 DELTACH3O2,DELTAHO2NO2,DELTAN2O5,DELTAHO2
      print 1803, 'CAP-17O    ',CAPO317,CAPCH3OOH17,
     2 CAPC3O217,CAPHO2NO217,CAPN2O517,CAPHO217

      print*, ''
      print 1803, 'TLOSSI     ',TLOSSI(7)+TLOSSI(8),
     2 TLOSSI(10)+TLOSSI(11),TLOSSI(12)+TLOSSI(13),
     3 TLOSSI(19)+TLOSSI(20),TLOSSI(22)+TLOSSI(23)
      print 1803, 'TLOSS      ',TLOSS(LO3),TLOSS(LCH3OOH),
     2 TLOSS(LCH3O2),TLOSS(LHO2NO2),TLOSS(LN2O5)
 1803 format(A18,1P6E10.3)
CCCC
C-PL COMMENT OUT
!      print*, ''
!      print 1811, 'TLOSSI     ',(TLOSSI(I),I=1,NQI),SO4LOSS,RATI(88)
!     2 ,TLOSSI_elS
!      print 1811, 'TLOSS      ',TLOSS(LH2S),TLOSS(LHS),TLOSS(LS),
!     2 TLOSS(LSO),TLOSS(LSO2),TLOSS(LH2SO4),TLOSS(LHSO),TLOSS(LS2),
!     3 TLOSS(LS3),TLOSS(LS4),TP(LSO4AER),TP(LS8AER),TLOSS_elS
!      print 1811, 'frac(removed)  ',(TLOSSI(I)/SOUT,I=1,NQI),
!     2 SO4LOSS/(S2VOLC+H2SVOLC),RATI(88)/SOUT,TLOSS_elS/SOUT


C-PL COMMENT OUT SINCE WE DON'T HAVE SOUT FOR OXYGEN QUESTION
!      print 1811, 'DEL*TLOSS/Sout ',TLOSS(LH2S)*PHIDEPI(1)/SOUT,
!     2 TLOSS(LHS)*PHIDEPI(2)/SOUT,TLOSS(LS)*PHIDEPI(3)/SOUT,
!     3 TLOSS(LSO)*PHIDEPI(4)/SOUT,TLOSS(LSO2)*PHIDEPI(5)/SOUT,
!     4 TLOSS(LH2SO4)*PHIDEPI(6)/SOUT,TLOSS(LHSO)*PHIDEPI(7)/SOUT,
!     5 2.*TLOSS(LS2)*PHIDEPI(8)/SOUT,3.*TLOSS(LS3)*PHIDEPI(9)/SOUT,
!     6 4.*TLOSS(LS4)*PHIDEPI(10)/SOUT,
!     7 SO4LOSS1*TP(LSO4AER)/SOUT,8.*S8LOSS*TP(LS8AER)/SOUT
!     8 ,del_elS*TLOSS_elS/SOUT
!      print *, "SUM(del*TLOSS) = ",TLOSS(LH2S)*PHIDEPI(1)/SOUT+
!     2 TLOSS(LHS)*PHIDEPI(2)/SOUT+TLOSS(LS)*PHIDEPI(3)/SOUT+
!     3 TLOSS(LSO)*PHIDEPI(4)/SOUT+TLOSS(LSO2)*PHIDEPI(5)/SOUT+
!     4 TLOSS(LH2SO4)*PHIDEPI(6)/SOUT+TLOSS(LHSO)*PHIDEPI(7)/SOUT+
!     5 2.*TLOSS(LS2)*PHIDEPI(8)/SOUT+3.*TLOSS(LS3)*PHIDEPI(9)/SOUT+
!     6 4.*TLOSS(LS4)*PHIDEPI(10)/SOUT+
!     2 SO4LOSS1*TP(LSO4AER)/SOUT+8.*S8LOSS*TP(LS8AER)/SOUT
!      print*, ""
!      print 1811, 'CAPDEL33*TLOS ',TLOSS(LH2S)*CAP33(1)/SOUT,
!     2 TLOSS(LHS)*CAP33(2)/SOUT,TLOSS(LS)*CAP33(3)/SOUT,
!     3 TLOSS(LSO)*CAP33(4)/SOUT,TLOSS(LSO2)*CAP33(5)/SOUT,
!     4 TLOSS(LH2SO4)*CAP33(6)/SOUT,TLOSS(LHSO)*CAP33(7)/SOUT,
!     5 2.*TLOSS(LS2)*CAP33(8)/SOUT,3.*TLOSS(LS2)*CAP33(9)/SOUT,
!     6 4.*TLOSS(LS4)*CAP33(10)/SOUT,
!     2 CAP33SO4*TP(LSO4AER)/SOUT,8.*CAP33S8*TP(LS8AER)/SOUT
!     3 ,elS_cap33*TLOSS_elS/SOUT
!      print*, "SUM(CapDel33*TL) = ",TLOSS(LH2S)*CAP33(1)/SOUT+
!     2 TLOSS(LHS)*CAP33(2)/SOUT+TLOSS(LS)*CAP33(3)/SOUT+
!     3 TLOSS(LSO)*CAP33(4)/SOUT+TLOSS(LSO2)*CAP33(5)/SOUT+
!     4 TLOSS(LH2SO4)*CAP33(6)/SOUT+TLOSS(LHSO)*CAP33(7)/SOUT+
!     5 2.*TLOSS(LS2)*CAP33(8)/SOUT+3.*TLOSS(LS3)*CAP33(9)/SOUT+
!     6 4.*TLOSS(LS4)*CAP33(10)/SOUT+
!     2 CAP33SO4*TP(LSO4AER)/SOUT+8.*CAP33S8*TP(LS8AER)/SOUT
!      print*, ''
!       print 1811, 'CAPDEL36*TLOS ',TLOSS(LH2S)*CAP36(1)/SOUT,
!     2 TLOSS(LHS)*CAP36(2)/SOUT,TLOSS(LS)*CAP36(3)/SOUT,
!     3 TLOSS(LSO)*CAP36(4)/SOUT,TLOSS(LSO2)*CAP36(5)/SOUT,
!     4 TLOSS(LH2SO4)*CAP36(6)/SOUT,TLOSS(LHSO)*CAP36(7)/SOUT,
!     5 2.*TLOSS(LS2)*CAP36(8)/SOUT,3.*TLOSS(LS3)*CAP36(9)/SOUT,
!     6 4.*TLOSS(LS4)*CAP36(10)/SOUT,
!     2 CAP36SO4*TP(LSO4AER)/SOUT,8.*CAP36S8*TP(LS8AER)/SOUT
!     3 ,elS_cap36*TLOSS_elS/SOUT
!      print*, "SUM(CapDel36*TL) = ",TLOSS(LH2S)*CAP36(1)/SOUT+
!     2 TLOSS(LHS)*CAP36(2)/SOUT+TLOSS(LS)*CAP36(3)/SOUT+
!     3 TLOSS(LSO)*CAP36(4)/SOUT+TLOSS(LSO2)*CAP36(5)/SOUT+
!     4 TLOSS(LH2SO4)*CAP36(6)/SOUT+TLOSS(LHSO)*CAP36(7)/SOUT+
!     5 2.*TLOSS(LS2)*CAP36(8)/SOUT+3.*TLOSS(LS3)*CAP36(9)/SOUT+
!     6 4.*TLOSS(LS4)*CAP36(10)/SOUT+
!     2 CAP36SO4*TP(LSO4AER)/SOUT+8.*CAP36S8*TP(LS8AER)/SOUT
      print*, ''
!      enddo
 1801 format(18x,10(A10))
 1811 format(A18,1P10E10.3)!10.2
!     
! 1031 format(A8,A5,A3,1PE10.2)
! 1034 format(I2,1x,1P10E10.2)
! 1035 format(3x,10(A8,2x))
      PRINT*, ''
      print 1032, 'cap-DELTA 17O profiles'
 1032 format(30x,A25)
!      PRINT 1033,'Z(km)','H2S','HS','HSO','SO','SO2','H2SO4',' | ',
!     2 'H2S','HS','HSO','SO','SO2','H2SO4'
! 1033 format(2x,A5,6(A7,3x),A3,6(A7,3x))

!      DO J=1,NZ
!      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
!     2 C2DELTAH2S(1,J),C2DELTAHS(1,J),C2DELTAHSO(1,J),
!     3 C2DELTASO(1,J),C2DELTASO2(1,J),C2DELTAH2SO4(1,J),' | ',
!     4 C2DELTAH2S(2,J),C2DELTAHS(2,J),C2DELTAHSO(2,J),
!     5 C2DELTASO(2,J),C2DELTASO2(2,J),C2DELTAH2SO4(2,J)
! 889  FORMAT(1X,F6.1,1P6E10.2,A3,1P6E10.2)
!      ENDDO    
CCCCCC     
      PRINT 1804,'Z(km)',(ISPECI(I),I=1,10)
 1804 format(16x,11(A10))
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
c     2 C2DELTAH2CO(J),C2DELTAO(J),C2DELTAH2O(J),C2DELTAOH(J),
c     3 C2DELTAHO2A(J), C2DELTAH2O2(J),C2DELTAO3A(J),C2DELTAO3B(J),
c     4 C2DELTACO(J),C2DELTACH3OOHA(J)
     2 C2DELTAH2CO(J),C2DELTAO(J),C2DELTAH2O(J),C2DELTAOH(J),
     3 C2DELTAHO2(J), C2DELTAH2O2(J),C2DELTAO3A(J),C2DELTAO3B(J),
     4 C2DELTACO(J),C2DELTACH3OOHA(J)
 889  FORMAT(18X,F6.1,1P10E10.2)
      ENDDO

      PRINT*, ''
      PRINT 1804,'Z(km)',(ISPECI(I),I=11,20)
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      if(J.le.10.or.modulo(J,4).eq.1)  PRINT 889, ZKM(J),
     2 C2DELTACH3OOHB(J),C2DELTACH3O2A(J),C2DELTACH3O2B(J),
     3 C2DELTAN2O(J), C2DELTANO(J),C2DELTANO2(J),C2DELTAHNO2(J),
     4 C2DELTAHNO3(J),C2DELTAHO2NO2A(J),C2DELTAHO2NO2B(J)
      ENDDO
C AJ 01/24/2020
      PRINT*, ''
      PRINT 1804,'Z(km)',(ISPECI(I),I=21,30)
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
     2 C2DELTANO3(J),C2DELTAN2O5A(J),C2DELTAN2O5B(J),C2DELTAO2(J),
     3 C2DELTASO(J), C2DELTASO2(J),C2DELTAH2SO4(J),C2DELTAHSO(J),
     4 C2DELTACO2(J) !,C2DELTAHO2B(J)
      ENDDO
C AJ

c      PRINT*, ''
c      PRINT 1804,'Z(km)',(ISPECI(I),I=31,NQI)
c      DO J=1,NZ
c      ZKM = Z(J)/1.E5
c      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
c     2 C2DELTAHO2NO2C(J)
c      ENDDO
! PRINT THE SPECIES HAVE TWO ISOTOPIC SPECIES
      PRINT*, ''
      print 1806, 'Z(km)','O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I'
     2 ,'HO2I'
 1806 format(17x,A7,1x,A7,5X,A7,2X,A7,4X,A7,2X,A7,2X,A7)
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
     2 C2DELTAO3(J),C2DELTACH3OOH(J),C2DELTACH3O2(J),
     3 C2DELTAHO2NO2(J),C2DELTAN2O5(J),C2DELTAHO2(J)
      ENDDO

C-AJ SKIP THE ISOTOPE CALCULATION 10/26
 6000 PRINT *,''
      PRINT *,'----------SKIP THE ISOTOPE CALCULATION-----------'
      PRINT *,''


C-PL COMMENT OUT
!      PRINT*, ''
!      print 1032, 'cap-DELTA 33S profiles','cap-DELTA 36S profiles'
! 1032 format(24x,A25,38x,A25)
!      PRINT 1033,'Z(km)','H2S','HS','HSO','SO','SO2','H2SO4',' | ',
!     2 'H2S','HS','HSO','SO','SO2','H2SO4'
! 1033 format(2x,A5,6(A7,3x),A3,6(A7,3x))
!      DO J=1,NZ
!      if(J.le.10.or.modulo(J,4).eq.1) PRINT 889, ZKM(J),
!     2 C2DELTAH2S(1,J),C2DELTAHS(1,J),C2DELTAHSO(1,J),
!     3 C2DELTASO(1,J),C2DELTASO2(1,J),C2DELTAH2SO4(1,J),' | ',
!     4 C2DELTAH2S(2,J),C2DELTAHS(2,J),C2DELTAHSO(2,J),
!     5 C2DELTASO(2,J),C2DELTASO2(2,J),C2DELTAH2SO4(2,J)
! 889  FORMAT(1X,F6.1,1P6E10.2,A3,1P6E10.2)
!      ENDDO
      
!      print*,''
!      PRINT 891, 'CapDel','33S','36S'
!      PRINT 892, 'Z(km)',('S','S2','S3','S4','S8',i=1,2)
!      DO J=1,10
!      PRINT 890, ZKM(J),C2DELTAS(1,J),C2DELTAS2(1,J),C2DELTAS3(1,J),
!     2 C2DELTAS4(1,J),C2DELTAS8(1,J),' | ',C2DELTAS(2,J),C2DELTAS2(2,J)
!     3 ,C2DELTAS3(2,J),C2DELTAS4(2,J),C2DELTAS8(2,J)
!      ENDDO
! 890  FORMAT(1X,F6.1,1P5E10.2,A3,1P5E10.2)
! 891  FORMAT(A6,25x,A3,52x,A3)
! 892  FORMAT(3x,A5,1x,5(A5,5x),' | ',5(A5,5x))
      
      print*, ''
!      print*, 'customized printout for total S outgassing sensitivity'
!      print*, 'S2/3/4/8/elS TLOSS, cap33, cap36/cap33'
!      print 1066,H2SVOLC+S2VOLC,2.*TLOSS(LS2),3.*TLOSS(LS3),
!     & 4.*TLOSS(LS4),8.*TLOSS(LS8AER),TLOSS_elS,CAP33(LSIS),
!     & CAP33(LSIS2),CAP33(LSIS3),CAP33S8,elS_cap33,
!     & CAP33(LSIS)/CAP36(LSIS),CAP33(LSIS2)/CAP36(LSIS2),
!     & CAP33(LSIS3)/CAP36(LSIS3),CAP36S8/CAP33S8,
!     & elS_cap36/elS_cap33
     
      print*, ''
C      print*, 'F(O2), then elemental S/sulfate/SO2/H2S TLOSS, cap33, 
C     & cap36/cap33' !S_volc
C       print 1036,SGFLUX(LO2),TLOSS_elS,TP(LSO4AER),TLOSS(LSO2), !SGFLUX(LO2) !H2SVOLC+S2VOLC
C     & TLOSS(LH2S),elS_cap33,CAP33SO4,CAP33(LSIO2),CAP33(LH2SI),
C     & elS_cap36/elS_cap33,CAP36SO4/CAP33SO4,CAP36(LSIO2)/CAP33(LSIO2),
C     & CAP36(LH2SI)/CAP33(LH2SI)
C-PL COMMENT OUT
!      print 1096,TLOSS_elS,TP(LSO4AER),TLOSS(LSO2),
!     & TLOSS(LH2S),elS_cap33,CAP33SO4,CAP33(LSIO2),CAP33(LH2SI),
!     & elS_cap36/elS_cap33,CAP36SO4/CAP33SO4,CAP36(LSIO2)/CAP33(LSIO2),
!     & CAP36(LH2SI)/CAP33(LH2SI),del_elS,SO4LOSS1,PHIDEPI(LSIO2),
!     & PHIDEPI(LH2SI)
      
!      print*, 'altitude, AERSOL(SO4, S8), RPAR(SO4, S8)'
!!      print*, 'altitude, N(S2/S3/S4), N_sat(S2/S3/S4)'
!!!      print*, 'S2/S3/S4 condensation rates vs. S2/S3/S4 consumption'
!!!      PRINT*, 'then photolysis of S2, S3, and S4 (/cm3/s for all!)'
      DO J=1,50
!      print 1076, ZKM(J),AERSOL(J,1),AERSOL(J,2),RPAR(J,1),RPAR(J,2)
!!      print 1086, ZKM(J),SL(LS2,J),SL(LS3,J),SL(LS4,J),
!!     & (N_Ssvap(I,J),I=2,4)
!!!      print 1086, ZKM(J),AR(231,J)*SL(LS2,J),AR(232,J)*SL(LS3,J),
!!!     & AR(233,J)*SL(LS4,J),
!!!     & AR(225,J)*SL(LS2,J)**2.,AR(226,J)*SL(LS3,J)*SL(LS,J),
!!!     & AR(227,J)*SL(LS4,J)**2.,
!!!     & AR(221,J)*SL(LS2,J),AR(229,J)*SL(LS3,J),AR(228,J)*SL(LS4,J)
      ENDDO
 1036 format(1P5E10.2,1x,1p4e10.2,1x,0p,4(f6.2,1x))
 1066 format(1P6E10.2,1x,1p5e10.2,1x,0p,5(f6.2,1x))
 1076 format(F6.1,1P4E10.2)
 1086 format(F6.1,1P9E10.2)
 1096 format(1P4E10.2,1x,0p,4(f6.1,1x),4(f6.2,1x),4(f6.1,1x))
      print*, ''
      
   !   endif !only print on moderately converged runs

      WRITE(118,1049), '#',0,(L,L=1,NSP)
 1049 FORMAT(1x,A2,1x,53(I8,2x))
      WRITE(118,1037), '#Z(km) ',' DENTOT ',(ISPEC(L),L=1,NSP)
 1037 FORMAT(1x,A6,3x,53(A8,2x))
      DO I=1,10
      WRITE(118,1038) I-0.5,DEN(I),(SL(J,I),J=1,NSP)
      ENDDO
 1038 FORMAT(1x,F6.1,1x,1P53E10.3)
      DO I=11,100,4
      WRITE(118,1038) I-0.5,DEN(I),(SL(J,I),J=1,NSP)
      ENDDO
      
!      PRINT 1050,AR(216,NZ),AR(105,NZ),AR(186,NZ),AR(117,NZ),
!     & AR(118,NZ),AR(222,NZ),AR(119,NZ),AR(122,NZ),AR(123,NZ),
!     & AR(127,NZ),AR(185,NZ),AR(229,NZ),AR(221,NZ),AR(228,NZ)
! 1050 format(1P14E10.2)
      
      PRINT*, "Success!"
!      DO I=1,100
!      PRINT 1039,I-0.5,SL(LS2,I)/DEN(I),SL(LS3,I)/DEN(I),
!     & SL(LS4,I)/DEN(I)
!      ENDDO
 1039 FORMAT('[',F4.1,' , ',1PE10.2,' , ',1PE10.2,' , ',1PE10.2,']')
      m_S = 2.*32.07/6.022E23 !g per S atom
      c_S = ((8.*1.38E-16*T(J))/(3.1415*2.*m_S))**0.5
!      print*,RPAR(10,1),AERSOL(10,1),c_S
      S2_lifetime = c_S*3.1415*RPAR(10,1)*AERSOL(10,1)
!      print*, 'Lifetime of S2 vs. cond. on SO4 (10km): ',S2_lifetime

C-PL PRINT OUT INTEGRATED REACTION RATE FOR EVERY SPECIES 06/12/2019
C
      OPEN(UNIT=90,FILE='OXYGEN_MIF/Reac_List.Oxy')
      READ(90,205) CHEMI
 205  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      DO 702 I=1,NSPI
         ISP = ISPECI(I)
         WRITE(16,703) ISP,TPI(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TPI = ',1PE9.2)
       DO 704 N=1,NRI 
          IF(JCHEMI(3,N).EQ.I .OR. JCHEMI(4,N).EQ.I .OR. 
     2       JCHEMI(5,N).EQ.I)THEN
         IF(RATI(N).NE.0.) WRITE(16,705) N,(CHEMI(J,N),J=1,5),RATI(N)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
      WRITE(16,706) ISP,TLI(I)
 706  FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TLI = ',1PE9.2)
       DO 707 N=1,NRI 
          IF(JCHEMI(1,N).EQ.I .OR. JCHEMI(2,N).EQ.I)THEN
      IF(RATI(N).NE.0.) WRITE(16,705) N,(CHEMI(J,N),J=1,5),RATI(N)
          ENDIF
 707   CONTINUE
 702  CONTINUE


      STOP
      END program atm_chem

c---------------------------------------------------------------


      SUBROUTINE OUTPUTP(N,NSTEPS,TIME,FLOW)
       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'
      INCLUDE 'INCLUDECHEM/comSATBLK.inc'
      INCLUDE 'INCLUDECHEM/comRRATS.inc'
c in may 10 2019, JL include Rblcok for JCHEM
      INCLUDE 'INCLUDECHEM/comRBLOK.inc'
c in may 10 2019, JL include chem; other two are place holder
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)  
      DIMENSION FUP(NQT),FLOW(NQT),CON(NQT),FLUXCH(NQT,NZ)
     2  ,ZF(NZ)
c in may 16 2019, JL include the following to print h2so4.pdat
      DIMENSION TTAB(51),PH2O(51,34),PH2SO4(51,34)
      REAL LIGHTNO2,LIGHTNO2I
C
c in may 10 2019, to include primo3s to this subroutine
      OPEN(unit=61, file= 'DATA'//'/primo3s.chm', status='old')
c in may 16 2019, JL add a clean printout for h2so4 to see possible disarray of input 
c      open(UNIT=1958, file='DIRIO'//'/h2so4.out.dat') 

      JS=N

      ISKIP = 4
      JSKIP = ISKIP
      IF(N.EQ.NSTEPS) ISKIP = 2
      TIMEY = TIME/3600./24./365.
      write(90,100)TIME,TIMEY
 100  FORMAT(/1X,'TIME =',E11.4,5X,'TIMEY =',1pe13.4,1X,'YEARS')
      write(90,101)NPHOT
 101  FORMAT(/1X,'NPHOT =',I3)
C
      write(90,105)
 105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES'/)
      IROW = 12
      LR = NQ/IROW + 1
      RL = FLOAT(NQ)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
 110  FORMAT(/5X,'Z',8X,13(A8,1X))
      DO 20 I=1,3
  20  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
      DO 21 I=4,NZ,ISKIP
  21  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
 120  FORMAT(1X,1P13E9.2)

      write(90,114) (ISPEC(K),K=K1,K2)
 114  FORMAT(5X,'Z',8X,13(A8,1X))
      IF (N.EQ.0) GO TO 8
      write(90,140)
 140  FORMAT(/1X,'TP, TL')
      write(90,145)(TP(K),K=K1,K2)
      write(90,145)(TL(K),K=K1,K2)
 145  FORMAT(10X,1P12E9.2)
   8  CONTINUE
C
C-AP
      write(90,106)
 106  FORMAT(//1X,'MIXING RATIO OF AEROSOL'/)
      write(90,185)
 185  FORMAT(5X,'Z',6X,'SO4AER')
      DO 18 J=1,3
  18  write(90,182) Z(J),SO4AER(J)
      DO 19 J=4,NZ,ISKIP
  19  write(90,182) Z(J),SO4AER(J)
 182  FORMAT(1X,1P2E9.2)
C-AP
      IF (N.EQ.0) RETURN
C
      write(90,183) TP(LSO4AER)
      write(90,184) TL(LSO4AER)
 183  FORMAT(/2X,'TP',6X,1P1E9.2)
 184  FORMAT(2X,'TL',6X,1P1E9.2)
C
      O3_ATMCM = O3COL/2.687e19
      write(90,150) O3COL
 150  FORMAT(//1X,'OZONE COLUMN DEPTH = ',1PE11.4)
      write(90,1150) O3_ATMCM
 1150 FORMAT(1X,'OZONE COLUMN DEPTH = ',F6.3,' ATM CM')
C-AP
      write(90,152) H2SCOL,SO2COL
 152  FORMAT(/1X,'SULFUR COLUMN DEPTHS:  H2S =',1PE10.3,2X,'SO2 =',
     2  E10.3,2X)
C-AP
	DO i = 1, NZ
	 IF (USOL(3,i) .LT. USOL(3,i+1)) THEN
		JCOLD = i
		GOTO 352
	 END IF
	END DO
 352  CONTINUE
      write(90,151) JCOLD, USOL(3,JCOLD)
 151  FORMAT(/1X,I3,' FH2O AT COLD TRAP =',1PE10.3)
      IF(N.LT.NSTEPS) RETURN
C
C ***** PRINT ON LAST ITERATION ONLY *****
      DO 1 I=1,NZ
   1  ZF(I) = Z(I) + 0.5*DZ
C
      DO 3 K=1,NQ
      DO 2 I=1,NZ
   2  SL(K,I) = USOL(K,I)*DEN(I)
!C-PL Added molecular diffusion terms below because they were 
!     missing (7/25/19)
      DO 4 I=1,NZ1
   4  FLUXCH(K,I) = - DK(K,I)*(USOL(K,I+1) - USOL(K,I))/DZ- 0.5*
     2  (HI(K,I)*DEN(I)*USOL(K,I) + HI(K,I+1)*DEN(I+1)*USOL(K,I+1))
   3  CONTINUE

C
      K = LSO4AER
      J = 1
      DO 5 I=1,NZ1
      FLUXCH(K,I) = -DK(K,I)*(SO4AER(I+1) - SO4AER(I))/DZ
     2  - 0.5*(WFALL(I,J)*DEN(I)*SO4AER(I) + WFALL(I+1,J)*DEN(I+1)
     3         *SO4AER(I+1))
   5  CONTINUE
C

!      print *,'Calculating lower boundary flux, FLOW, for O2'
      DO 15 K=1,NQT
      FLOW(K) = FLUXCH(K,1) - (YP(K,1) - YL(K,1)*SL(K,1))*DZ
C      if(k.eq.LO2) print 1500,FLUXCH(K,1),YP(K,1),YL(K,1),SL(K,1),DZ
C 1500 format(1x,1p5e10.3)      
      FUP(K) = FLUXCH(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*SL(K,NZ))*DZ
      CON(K) = TP(K) - TL(K) + FLOW(K) - FUP(K)
 502  format('K=',I2,' FLOW(K)=',1PE10.3,' FUP(K)=',E10.3)
  15  CONTINUE

C-AP
      FLOW(3) = FLUXCH(3,11)
      CON(3) = TP(3) - TL(3) + FLOW(3) - FUP(3)
      DO 6 I=1,10
   6  FLUXCH(3,I) = 0.
C
      write(90,125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)
      DO 9 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1+IROW - 1
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 22 I=1,NZ,ISKIP
       write(90,120) Z(I),(SL(K,I),K=K1,K2)
c      if ((LO3.GE.K1).AND.(LO3.LE.K2)) WRITE(83,355) Z(I),SL(LO3,I)
  22  CONTINUE
   9  CONTINUE
 355  FORMAT(1PE10.3,1PE12.3)
c     close(83)
      write(90,185)
      DO J=1,3
      write(90,182) Z(J),SL(LSO4AER,J)
      END DO
      DO J=4,NZ,ISKIP
      write(90,182) Z(J),SL(LSO4AER,J)
      END DO

      ISKIP = 1
      write(90,155)
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED SPECIES'/)
      ZFL = 0.
      ZFT = ZF(NZ)
     
      DO 10 L=1,LR+1
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR+1) THEN
       K1 = NQT
       K2 = K1
      END IF
C-AP      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,120) ZFL,(FLOW(K),K=K1,K2)
      DO 23 I=1,NZ,ISKIP
  23  write(90,120) ZF(I),(FLUXCH(K,I),K=K1,K2)
      write(90,120) ZFT,(FUP(K),K=K1,K2)
  10  CONTINUE

C-PL RECALCULATE TLOSS   07/2019
      DO I=1,NQ
      IF (LBOUND(I).EQ.0) THEN
      TLOSS(I) = SR(I) + PHIDEP(I)
      ELSE IF (LBOUND(I).NE.0.AND.FLOW(I).LT.0) THEN
      TLOSS(I) = SR(I) -FLOW(I) 
      ELSE     
      TLOSS(I) = SR(I)
      END IF
      END DO

      write(90,175)
 175  FORMAT(/1X,'RAINOUT RATE, PHIDEP, TLOSS AND LOWER B.C.'/)
      write(90,176)
 176  FORMAT(1X,'FOLLOWED BY TP, TL, FUP, FLOW, CON'/)
      DO 13 L=1,LR+1
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
C      IF (L.EQ.LR) K2 = NQT
      IF (L.EQ.LR+1) THEN
       K1 = NQT
       K2 = K1
      END IF
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,145) (SR(K),K=K1,K2)
      write(90,145)(PHIDEP(K),K=K1,K2)
      write(90,145)(TLOSS(K),K=K1,K2) 
      write(90,146) (LBOUND(K),K=K1,K2)
 146  FORMAT(14X,12(I1,8X))
      write(90,145)
      write(90,145) (TP(K),K=K1,K2)
      write(90,145) (TL(K),K=K1,K2)
      write(90,145) (FUP(K),K=K1,K2)
      write(90,145) (FLOW(K),K=K1,K2)
      write(90,145) (CON(K),K=K1,K2)     
  13  CONTINUE

C-AP 175  FORMAT(/1X,'RAINOUT OF LONG-LIVED SPECIES'/)
C-AP      DO 13 L=1,LR
C-AP      K1 = 1 + (L-1)*IROW
C-AP      K2 = K1 + IROW - 1
C-AP      IF (L.EQ.LR) K2 = NQ
C-AP      PRINT 110, (ISPEC(K),K=K1,K2)
C-AP  13  PRINT 145,(SR(K),K=K1,K2)
C-AP

C   COMPUTE CONSERVATION OF SULFUR
!      SULDEP = - (FLOW(LHS) + FLOW(LSO) + FLOW(LH2SO4)
!     2  + FLOW(LHSO) + FLOW(LSO4AER))
!      IF (LBOUND(LSO2).EQ.0) SULDEP = SULDEP - FLOW(LSO2)
!      IF (LBOUND(LH2S).EQ.0) SULDEP = SULDEP - FLOW(LH2S)
!      SULRAN = SR(LH2S) + SR(LHS) + SR(LSO) + SR(LSO2) +
!     2  SR(LH2SO4) + SR(LHSO) + SR(LSO4AER) 
!      SULLOS = SULDEP + SULRAN
!      SULPRO = 0.
!      IF (LBOUND(LSO2).NE.0) SULPRO = SULPRO + FLOW(LSO2)
!      IF (LBOUND(LH2S).NE.0) SULPRO = SULPRO + FLOW(LH2S)
!      SO4LOS = TLOSS(LH2SO4) + TLOSS(LSO4AER)
!      write(90,177) SULLOS,SULPRO,SO4LOS
! 177  FORMAT(/1X,'CONSERVATION OF SULFUR:',/5X,'SULLOS =',1PE10.3,
!     2  2X,'SULPRO =',E10.3,2X,'SO4LOS =',E10.3)

C-PL COMPUTE CONSERVATION OF OXYGEN  07/2019 
      OXYDEP = 0.
      OXYDEP = - (FLOW(LH2CO) + FLOW(LO) + FLOW(LOH) + FLOW(LHO2)*2.
     2  + FLOW(LH2O2)*2. + FLOW(LO3)*3. + FLOW(LCH3OOH)*2.
     3  + FLOW(LCH3O2)*2. + FLOW(LNO) + FLOW(LNO2)*2. + FLOW(LHNO2)*2.
     4  + FLOW(LHNO3)*3. + FLOW(LHO2NO2)*4. + FLOW(LNO3)*3. 
     5  + FLOW(LN2O5)*5. + FLOW(LSO) + FLOW(LSO2)*2. 
     6  + FLOW(LH2SO4) * 4. + FLOW(LHSO)  )
!     7  + FLOW(LCO2)*2. ) !
!      IF (LBOUND(LCO).EQ.0) OXYDEP = OXYDEP - FLOW(LCO)
!      IF (LBOUND(LN2O).EQ.0) OXYDEP = OXYDEP - FLOW(LN2O)
      OXYRAN = SR(LH2CO) + SR(LO) + SR(LOH) + SR(LHO2) * 2.
     2  + SR(LH2O2)*2. + SR(LO3)*3. + SR(LCH3OOH)*2.
     3  + SR(LCH3O2)*2. + SR(LNO) + SR(LNO2)*2. + SR(LHNO2)*2.
     4  + SR(LHNO3)*3. + SR(LHO2NO2)*4. + SR(LNO3)*3. 
     5  + SR(LN2O5)*5. + SR(LSO) + SR(LSO2)*2. 
     6  + SR(LH2SO4) * 4. + SR(LHSO) + TP(LSO4AER) * 4.
     7  + SR(LCO) + SR(LN2O) + SR(LO2)*2. + SR(LCO2) *2.
       OXYUP = FUP(LH2CO) + FUP(LO) + FUP(LOH) + FUP(LHO2) * 2.
     2  + FUP(LH2O2)*2. + FUP(LO3)*3. + FUP(LCH3OOH)*2.
     3  + FUP(LCH3O2)*2. + FUP(LNO) + FUP(LNO2)*2. + FUP(LHNO2)*2.
     4  + FUP(LHNO3)*3. + FUP(LHO2NO2)*4. + FUP(LNO3)*3. 
     5  + FUP(LN2O5)*5. + FUP(LSO) + FUP(LSO2)*2. 
     6  + FUP(LH2SO4) * 4. + FUP(LHSO) 
     7  + FUP(LCO) + FUP(LN2O) + FUP(LO2)*2. + FUP(LCO2) *2.
      OXYLOS = OXYDEP + OXYRAN + OXYUP + TP(LH2O)-TL(LH2O)
      OXYPRO = 0.
      OXYPRO = FLOW(LCO) + FLOW(LN2O)!SGFLUX(LCO) + SGFLUX(LN2O)     
      OXYPRO = OXYPRO + 2.*GPPOXY + FLOW(LO2)*2.
     2                + 2.*GPPCDE + FLOW(LCO2)*2.

      CONOXY = OXYLOS - OXYPRO
      write(90,177)OXYLOS,OXYPRO,CONOXY,CONOXY/MIN(OXYLOS,OXYPRO)*100.
 177  FORMAT(/1X,'CONSERVATION OF OXYGEN:',/5X,'OXYLOS =',1PE10.3,
     2  2X,'OXYPRO =',E10.3,2X,'CONOXY =',E10.3,
     3  3X,'ERR = ',E10.3,'%' )
      write(90,178) OXYDEP,OXYRAN,OXYUP,TP(LH2O)-TL(LH2O)
 178  FORMAT(/5X,'OXYDEP =',1PE10.3,2X,'OXYRAN =',E10.3,
     2 2X,'OXYUP =',E10.3,2X,'TP-TL(H2O) =',2E10.3)

      write(90,189) FLOW(LCO),FLOW(LN2O),2.*GPPOXY+FLOW(LO2)*2.
     2 ,GPPOXY,2.*GPPCDE+FLOW(LCO2)*2,GPPCDE
 189  FORMAT(5X,'FLOW(LCO) =',1PE10.3,2X,'FLOW(LN2O) =',E10.3,
     2 /5X,'GPP-FLOW(LO2)=',E10.3,2X,'GPPOXY =',E10.3,
     3 2X,'GPP-FLOW(LCO2)=',E10.3,2X,'GPPCDE =',E10.3)

      write(90,288) TPH2O,TLH2O,TPH2O-TLH2O,
     2 TP(LH2O)-TPH2O,TL(LH2O)-TLH2O,TP(LH2O)-TPH2O-(TL(LH2O)-TLH2O),
     3 TP(LH2O),TL(LH2O),TP(LH2O)-TL(LH2O)
 288  FORMAT(/5X,' H2O     TP        TL        TP-TL ',
     2       /5X, 'TROPOS',1P3E10.3,
     2       /5X, 'STRATO',1P3E10.3,
     3       /5X, 'TOTAL ',1P3E10.3)

      write(90,285) ABS(FLOW(LO2)),SL(LO2,1)
 285  FORMAT(5x,'FLOW(O2)=',F20.3,/5X,'SL(O2)=',F25.3)
      write(90,286) ABS(FLOW(LCO2)),SL(LCO2,1)
 286  FORMAT(5x,'FLOW(CO2)=',F20.3,/5X,'SL(CO2)=',F25.3)

      write(90,388)
 388  FORMAT(/1X,'RAINOUT RATE, FLOW, FUP,TLOSS')

      write(90,110)  ISPEC(LH2CO),ISPEC(LO),ISPEC(LOH),ISPEC(LHO2),
     2 ISPEC(LH2O2),ISPEC(LO3),ISPEC(LCO),ISPEC(LCH3OOH),ISPEC(LCH3O2),
     3 ISPEC(LN2O)
      write(90,149) 'RAIN',SR(LH2CO),SR(LO),SR(LOH),SR(LHO2)*2.,
     2 SR(LH2O2)*2,SR(LO3)*3.,SR(LCO),SR(LCH3OOH)*2.,SR(LCH3O2)*2.,
     3 SR(LN2O)
      write(90,149) 'FLOW',FLOW(LH2CO),FLOW(LO),FLOW(LOH),
     2 FLOW(LHO2)*2.,FLOW(LH2O2)*2.,FLOW(LO3)*3.,FLOW(LCO),
     3 FLOW(LCH3OOH)*2.,FLOW(LCH3O2)*2.,FLOW(LN2O)
      write(90,149) ' FUP',FUP(LH2CO),FUP(LO),FUP(LOH),
     2 FUP(LHO2)*2.,FUP(LH2O2)*2.,FUP(LO3)*3.,FUP(LCO),
     3 FUP(LCH3OOH)*2.,FUP(LCH3O2)*2.,FUP(LN2O)
      write(90,149) 'TLOS',TLOSS(LH2CO),TLOSS(LO),TLOSS(LOH),
     2 TLOSS(LHO2)*2.,TLOSS(LH2O2)*2.,TLOSS(LO3)*3.,TLOSS(LCO),
     3 TLOSS(LCH3OOH)*2.,TLOSS(LCH3O2)*2.,TLOSS(LN2O)
      write(90,*)

      write(90,110)  ISPEC(LNO),ISPEC(LNO2),ISPEC(LHNO2),ISPEC(LHNO3),
     2 ISPEC(LHO2NO2),ISPEC(LNO3),ISPEC(LN2O5),ISPEC(LO2),ISPEC(LSO),
     3 ISPEC(LSO2)
      write(90,149) 'RAIN',SR(LNO),SR(LNO2)*2.,SR(LHNO2)*2.,
     2 SR(LHNO3)*3.,SR(LHO2NO2)*4.,SR(LNO3)*3.,SR(LN2O5)*5.,
     3 SR(LO2)*2.,SR(LSO),SR(LSO2)*2.
      write(90,149) 'FLOW',FLOW(LNO),FLOW(LNO2)*2.,FLOW(LHNO2)*2.,
     2 FLOW(LHNO3)*3.,FLOW(LHO2NO2)*4.,FLOW(LNO3)*3.,FLOW(LN2O5)*5.
     3 ,FLOW(LO2)*2.,FLOW(LSO),FLOW(LSO2)*2.
      write(90,149) ' FUP',FUP(LNO),FUP(LNO2)*2.,FUP(LHNO2)*2.,
     2 FUP(LHNO3)*3.,FUP(LHO2NO2)*4.,FUP(LNO3)*3.,FUP(LN2O5)*5.,
     3 FUP(LO2)*2.,FUP(LSO),FUP(LSO2)*2.
      write(90,149) 'TLOS',TLOSS(LNO),TLOSS(LNO2)*2.,TLOSS(LHNO2)*2.,
     2 TLOSS(LHNO3)*3.,TLOSS(LHO2NO2)*4.,TLOSS(LNO3)*3.,
     3 TLOSS(LN2O5)*5.,TLOSS(LO2)*2.,TLOSS(LSO),TLOSS(LSO2)*2.
      write(90,*)

      write(90,110)  ISPEC(LH2SO4),ISPEC(LHSO),
     2 ISPEC(LCO2),ISPEC(LSO4AER)
      write(90,149) 'RAIN',SR(LH2SO4)*4.,SR(LHSO),
     2 SR(LCO2)*2.,SR(LSO4AER)*4
      write(90,149) 'FLOW',FLOW(LH2SO4)*4.,FLOW(LHSO),
     2 FLOW(LCO2)*2.,FLOW(LSO4AER)*4.
      write(90,149) ' FUP',FUP(LH2SO4)*4.,
     2 FUP(LHSO),FUP(LCO2)*2.,FUP(LSO4AER)*4.
      write(90,149) 'TLOS',TLOSS(LH2SO4)*4.,TLOSS(LHSO),
     2 TLOSS(LCO2)*2.,TLOSS(LSO4AER)*4.

 149  FORMAT(A5,5X,1P12E9.2)


C-AP
      write(90,179)
  179 FORMAT(/1X,'INTEGRATED REACTION RATES'/)
c      write(90,180) RAT
c  180 FORMAT(1X,1P10E10.3)
C
      WRITE(90,181)
      IROW = 10
      LR = NR/IROW + 1
      RL = FLOAT(NR)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 17 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) THEN
        K2 = NR
        WRITE(90,186) K1,(RAT(K),K=K1,K2),K2
  186   FORMAT(I3,2X,1P2E10.3,82X,I3)
        GO TO 17
      ENDIF
      WRITE(90,180) K1,(RAT(K),K=K1,K2),K2
  180 FORMAT(I3,2X,1P10E10.3,2X,I3)
   17 CONTINUE
      WRITE(90,181)
  181 FORMAT(9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,
     2    '8',9X,'9',8X,'10')
C
      write(90,160)
 160  FORMAT(/1X,'ATMOSPHERIC PARAMETERS AND PH EQ SPECIES')
      NPE = NSP - NQ
      NQ1 = NQ + 1
      LR = NPE/IROW + 1
      RL = FLOAT(NPE)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 12 L=1,LR
      K1 = NQ1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NSP
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 24 I=1,NZ,ISKIP
  24  write(90,120) Z(I),(SL(K,I),K=K1,K2)
  12  CONTINUE
C
      write(90,190)
 190  FORMAT(/1X,'ATMOSPHERIC PARAMETERS')
      write(90,195)
 195  FORMAT(/4X,'Z',9X,'T',9X,'EDD',7X,'DEN')
      write(90,200)(Z(I),T(I),EDD(I),DEN(I),I=1,NZ,ISKIP)
 200  FORMAT(1X,1P4E10.3)
C
      write(90,230)
 230  FORMAT(/1X,'SULFATE AEROSOL PARAMETERS')
      write(90,235)
 235  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'FSULF',4X,
     2  'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'H2SO4S',4X,'H2SO4',5X,
     3  'CONSO4',4X,'CONVER')
      write(90,240) (Z(I),AERSOL(I,1),RPAR(I,1),
     & WFALL(I,1),FSULF(I),
     2  TAUSED(I,1),TAUEDD(I),TAUC(I,1),H2SO4S(I),USOL(LH2SO4,I),
     3  CONSO4(I),CONVER(I,1),I=1,NZ,ISKIP)
 240  FORMAT(1X,1P12E10.3)


c      print *, 'TAUC(I,1) =', TAUC(I,1) 
C
C   Calculate total hydrogen mixing ratio and write to a file
      write(86,122)
 122  format(5x,'zkm',6x,'totH2O',4x,'totH',6x,'totH2',4x,'totCH4',
     2  4x,'tothyd')
      DO I=1,NZ
      zkm = Z(I)/1.E5
      totH2O = USOL(LH2O,I)
      totH = USOL(LH,I)*0.5
      totH2 = USOL(LH2,I)
      totCH4 = USOL(LCH4,I)*2.
      tothyd = totH2O + totH + totH2 + totCH4
      if(I.eq.JTROP) tothytrop = tothyd
      if(I.eq.90) tothyhom = tothyd
      write(86,121) zkm,totH2O,totH,totH2,totCH4,tothyd
      END DO
      close(86)  
 121  FORMAT(1x,1p6e10.3)
C
C Calculate hydrogen escape (in units of H2) based on tothyd at the top
C of the model. Do this both approximately and accurately.
C Approximate
      H2escape1 = 2.5e13*tothyhom
C Accurate
      H2escape2 = FUP(LH2O) + 0.5*FUP(LH) + FUP(LH2) + 2.*FUP(LCH4)
C
      write(90,300)
 300  format(/'Hydrogen escape rate in units of H2 molecules')
      write(90,301) H2escape1
 301  format(5x,'Approximate escape rate =',1pe10.3)
      write(90,302) H2escape2
 302  format(5x,'Escape rate calculated from FUP =',1pe10.3)
      write(90,303) tothyhom
 303  format(5x,'Total hydrogen mixing ratio at the homopause = ',
     2  1pe10.3)
      write(90,304) tothytrop
 304  format(5x,'Total hydrogen mixing ratio at the tropopause =',
     2  1pe10.3)

c in may 9 2019, JL add reaction rate table to the code to make life easier 

      if (N.EQ.NSTEPS) THEN
      REWIND 61
      READ(61,2000)CHEM
2000  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      write(15,2001)(J,(CHEM(M,J),M=1,5),J=1,NR)
2001  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
      DO 702 I=1,NSP
         ISP = ISPEC(I)
         WRITE(15,703) ISP,TP(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TP = ',1PE9.2)
       DO 704 Nj=1,NR  
          IF(JCHEM(3,Nj).EQ.I .OR. JCHEM(4,Nj).EQ.I .OR. 
     2       JCHEM(5,Nj).EQ.I)THEN
           IF(RAT(Nj).NE.0.) WRITE(15,705) Nj,(CHEM(J,Nj),J=1,5),RAT(Nj)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
      WRITE(15,706) ISP,TL(I)
 706     FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TL = ',1PE9.2)
       DO 707 Nj=1,NR 
          IF(JCHEM(1,Nj).EQ.I .OR. JCHEM(2,Nj).EQ.I)THEN
          IF(RAT(Nj).NE.0.) WRITE(15,705) Nj,(CHEM(J,Nj),J=1,5),RAT(Nj)
          ENDIF
 707   CONTINUE
 702  CONTINUE
      close(61)
      END IF

      RETURN
      END
C========================================================================
      subroutine gaussian_data(xi,wi,NumGau)
      parameter(nrow=11)
      dimension xi(nrow,20),wi(nrow,20),NumGau(nrow)
 800  format(2x, I2)
 801  format(F7.5,1X,F7.5)
      DO i = 1,nrow
      read(68,800) NumGau(i)
         DO j = 1,NumGau(i)
         read(68,801) xi(i,j), wi(i,j)
         END DO 
      END DO
      
      RETURN
      END
