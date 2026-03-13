      program atm_chem
 
C     This program was created by James Kasting (JK) and modified in the early
C     2000's by Alex Pavlov (AP) and Kara Krelove (KK). The user-friendly 
C     version of this code was created by Antigona Segura (AS) in 2005.
C     Other students who have modified it include Peng Liu (PL), Jingjun
C     Liu (JL) and Aoshuang Ji (AJ). Some modifications are identified with
C     the author's initials.
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
c     Martin Cohen from new IUE (+ model) data. 
c
c     Check the notes along the main program before using it. Look for 
c     the word 'NOTE:'
C
C     Rate constants are taken from the NIST (National Institute of
C     Standards) Chemical Kinetics Database (https://kinetics.nist.gov/kinetics/).
C     This contains information that used to be published by JPL in
C     their regular monographs on chemical kinetics.
C
C     This version of the program contains a new correlated-k parameterization
C     for O2 photolysis in the Schumann-Runge bands developed by Aoshuang Ji
C     and JK with help from Rafa Fernandez and his group in Argentina, Guillaume
C     Chaverot in France, and Ken Minschwaner at New Mexico Tech. A paper led
C     by Aoshuang is expected to be published in 2024. This replaces the old
C     O2 photolysis routine, which was based on the Allen & Frederick (1982)
C     band model and JK's exponential sum fit to that model. The exponential
C     sum fit did not include temperature dependence, and hence overestimated
C     O2 absorption at temperatures below 300 K.

C     NO PHOTOLYSIS IS TREATED BY THE MODIFIED CIESLIK AND NICOLET METHOD.
C     THIS EDITION ALSO HAS GIORGI AND CHAMEIDES RAINOUT RATES AND A 
C     MANABE-WETHERALD RELATIVE HUMIDITY DISTRIBUTION. 
C
C       THIS PROGRAM IS A ONE-DIMENSIONAL MODEL OF THE PRIMORDIAL
C     ATMOSPHERE.  THE MIXING RATIOS OF THE LONG-LIVED SPECIES
C     ARE CALCULATED FROM THE EQUATION
C
C     DF/DT  =  (1/N)*D/DZ(KN*DF/DZ) + P/N - LF + molecular diffusion terms
C
C     WHERE
C     F = MIXING RATIO (USOL)
C     K = EDDY DIFFUSION COEFFICIENT (EDD)
C     N = TOTAL NUMBER DENSITY (DEN)
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C
C     Note: The original version of this program only included eddy diffusion
C     and was only extended up to 64 km. (64 was the vector length on
C     the Cray 'supercomputer'.) Molecular diffusion was added in the late
C     1990's in order to simulate CH4-rich atmospheres that formed organic
C     haze at high altitudes. The current code extends to 100 km in 1-km
C     increments. A partial description of the code can be found in
C     Catling & Kasting, Atmospheric Evolution on Inhabited and Lifeless
C     Worlds (2017), Appendix B,
C
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTON'S
C     METHOD, i.e., IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.
C
C     Most subroutines are found in separate directories. Two of them
C     are found at the end of this main program.
C
C     Main program subroutines:
C     OUTPUTP -  PRINTS OUT RESULTS
C     GAUSSIAN_DATA -- SETS UP GAUSS POINTS FOR THE SOLAR ZENITH ANGLE
C                      INTEGRATION (IF CHOSEN)
C
C     Subroutines found in subdirectory SETUP:
C     GRID   -  SETS UP THE ALTITUDE GRID
C     RATES  -  DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE.
C               NOTE THAT THE RAINOUT RATE DOES NOT DEPEND ON CLOUD
C               DROPLET PH IN THIS VERSION OF THE MODEL, SO THERE IS
c               NO SUBROUTINE AQUEOUS, WHICH IS PART OF THE LOW-O2 MODEL.
C     DENSTY -  COMPUTES ATMOSPHERIC DENSITIES FROM HYDROSTATIC
C               EQUILIBRIUM AND INITIALIZES ABSORBER PROFILES
C     DIFCO  -  COMPUTES DK = K*N BETWEEN GRID POINTS
C     PSATRAT - CALCULATES THE SATURATION VAPOR PRESSURE OF WATER
C     READPHOTO - READS THE VERY LONG MAIN PHOTOLYSIS DATAFILE AND
C                 THE MUCH SHORTER FAR-UV DATAFILE
C     READSTAR - CAN BE USED TO READ FAR-UV DATAFILES FOR OTHER STARS
C
C     Subroutines and functions found in subdirectory CHEM:
C     DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C              POINTS BY CALLING CHEMPL
C     CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C              FOR ONE SPECIES AT ALL GRID POINTS
C     LTNING - COMPUTES LIGHTNING PRODUCTION RATES FOR O2 AND
C              N2 BASED ON CHAMEIDES' RESULTS
C     PHOTO  -  COMPUTES PHOTOLYSIS RATES (CALLS TWOSTR)
C     TWOSTR -  CALCULATES RADIATIVE TRANSFER USING THE TOON ET AL. (1989)
C               DELTA 2-STREAM METHOD
C     O3PHOT - COMPUTES O(1D) QUANTUM YIELDS IN OZONE PHOTOLYSIS          
C     TBDY   -  FUNCTION THAT COMPUTES 3-BODY REACTION RATES
C
C     Subroutines found in subdirectory MATRIX
C
C     NOT DESCRIBED HERE, BUT THESE ARE THE LINPACK SUBROUTINES NEEDED
C     TO SOLVE THE BANDED JACOBIAN MATRIX USED IN THE LONG-LIVED SPECIES
C     INTEGRATION, AS WELL AS THE TRIDIAGONAL MATRIX USED TO CALCULATE
C     THE VERTICAL DISTRIBUTION OF AEROSOL PARTICLES.
C
C     Subroutines found in subdirectory PRTCL:
C
C     NOT DESCRIBED HERE, BUT THESE ARE THE SUBROUTINES NEEDED TO FIND
C     THE PROPERTIES OF THE SULFATE AEROSOLS. THESE ROUTINES WERE
C     BORROWED FROM BRIAN TOON WHEN HE AND JK WERE BOTH AT NASA AMES.
C
C     Subroutines found in subdirectory OXYGEN_MIF:
C
C     NOT DESCRIBED HERE, BUT THESE ARE THE SUBROUTINES CREATED BY
C     PENG LIU AND JINGUN LIU TO DO THE MASS-INDEPENDENT O ISOTOPE
C     CALCULATIONS DESCRIBED IN P. LIU ET AL. (PNAS, 2021). THESE
C     SUBROUTINES WERE MODIFIED FROM ONES WRITTEN BY ALEX PAVLOV
C     20 YEARS PREVIOUSLY FOR THE SULFUR MIF CALCULATIONS DESCRIBED
C     IN PAVLOV AND KASTING (ASTROBIOLOGY, 2002).
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (2) E1     - EXPONENTIAL INTEGRAL OF ORDER ONE (no longer used)
C
C
C     Note: The reaction list found below does not do anything! It is
C     included simply so that the program can be understood without the
C     accompanying datafile, primo3s.dat, which is found in subdirectory
C     DATA.
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
C  ************* Starting chlorine chemistry ************
C         (but note also reactions 102, 103, 140)
C   101)  CH3Cl + HV  =  CH3 + Cl
C   102)  NO + NO3  =  NO2 + NO2
C   103)  OH + NO3  =  HO2 + NO2
C   104)  CH3Cl + OH  =  CLO + H2COC   105)  CL + O3  =  ClO + O2
C   106)  Cl + H2  = HCl + H
C   107)  Cl + CH4  =  HCl + CH3
C   108)  Cl + CH3Cl  =  CH2Cl + HCl
C   109)  Cl + H2CO  =  HCl + HCO
C   110)  Cl + H2O2  =  HCl + HO2
C   111)  Cl + HO2  =  HCl + O2
C   112)  Cl + HO2  =  ClO + OH
C   113)  Cl + ClONO2  =  Cl + Cl + NO2 (+O)
C   114)  Cl + NO + M  = NOCl + M
C   115)  Cl + NO2 + M  =  ClONO + M
C   116)  Cl + NOCl  =  NO + Cl2
C   117)  Cl + O2 + M  =  ClO2 + M
C   118)  Cl + ClO2  =  Cl2 + O2
C   119)  Cl + ClO2  =  ClO + ClO
C   120)  ClO + O  =  Cl + O2
C   121)  ClO + NO  =  Cl + NO2
C   122)  ClO + NO2 + M  =  ClONO2 + M
C   123)  ClO + HO2  =  HOCl + O2
C   124)  ClO + OH  =  CL + HO2
C   125)  HCl + OH  =  Cl + H2O
C   126)  HOCl + OH  =  ClO + H2O
C   127)  ClONO2 + OH  =  Cl + HO2 + NO2
C   128)  HCl + O  =  Cl + OH
C   129)  HOCl + O  =  ClO + OH
C   130)  ClONO2 + O  =  Cl + O2 + NO2
C   131)  Cl2 + OH  =  HOCl + Cl
C   132)  Cl2 + hv  =  Cl + Cl
C   133)  ClO2 + hv  =  ClO + O
C   134)  HCl + hv  =  H + Cl
C   135)  HOCl + hv  =  OH + Cl
C   136)  NOCl + hv  =  Cl + NO
C   137)  ClONO + hv  =  Cl + NO2
C   138)  CLONO2 + hv  =  Cl + NO3
C   139)  ClO2 + hv  =  Cl + O2
C   140)  HO2 + NO3  =  HNO3 + O2
C   141)  ClO + ClO + M  = Cl2O2 + M
C   142)  Cl2O2 + hv  =  ClO2 + Cl
C   143)  Cl2O2 + M  =  ClO + ClO
C   144)  ClO2 + M  =  Cl + O2
C   145)  Cl + NO3  =  ClO + NO2
C   146)  Cl + HOCl  =  Cl2 + OH
C   147)  ClO + NO3  =  ClONO + O2
C   148)  ClONO + OH  =  HOCl + NO2
C   149)  ClO2 + O  =  ClO + O2
C ************* Additional NOx reactions *********
C   150)  NO2 + O + M  =  NO3 + M
C   151)  NO3 + hv  =  NO2 + O
C   152)  NO3 + NO2 + M  =  N2O5 + M
C   153)  N2O5 + hv  =  NO2 + NO3
C   154)  N2O5 + M  =  NO2 + NO3 + M
C   155)  N2O5 + H2O  =  2 HNO3
C************** Starting sulfur chemistry ********
C   156)  SO   + HV   =     S    +     O
C   157)  SO2  + HV   =     SO   +     O
C   158)  H2S  + HV   =     HS   +     H
C   159)  SO   + O2   =     O    +     SO2
C   160)  SO   + HO2  =     SO2  +     OH
C   161)  SO   + O    =     SO2
C   162)  SO   + OH   =     SO2  +     H
C   163)  SO2  + OH   =     HSO3
C   164)  SO2  + O    =     SO3
C   165)  SO3  + H2O  =     H2SO4
C   166)  HSO3 + O2   =     HO2  +     SO3
C   167)  HSO3 + OH   =     H2O  +     SO3
C   168)  HSO3 + H    =     H2   +     SO3
C   169)  HSO3 + O    =     OH   +     SO3
C   170)  H2S  + OH   =     H2O  +     HS
C   171)  H2S  + H    =     H2   +     HS
C   172)  H2S  + O    =     OH   +     HS
C   173)  HS   + O    =     H    +     SO
C   174)  HS   + O2   =     OH   +     SO
C   175)  HS   + HO2  =     H2S  +     O2
C   176)  HS   + HS   =     H2S  +     S
C   177)  HS   + HCO  =     H2S  +     CO
C   178)  HS   + H    =     H2   +     S
C   179)  HS   + S    =     H    +     S2
C   180)  S    + O2   =     SO   +     O
C   181)  S    + OH   =     SO   +     H
C   182)  SO2  + HV   =     S    +     O2
C   183)  S    + HO2  =     HS   +     O2
C   184)  S    + HO2  =     SO   +     OH
C   185)  HS   + H2CO =     H2S  +     HCO
C   186)  SO2  +     HV  =      SO21
C   187)  SO2  +     HV  =      SO23
C   188)  H2SO4 +     HV =       SO2   +   OH  +    OH
C   189)  SO3   +    HV  =      SO2  +     O
C   190)  SO21  +    M   =      SO23 +     M
C   191)  SO21  +    M   =      SO2  +     M
C   192)  SO21  +    HV  =      SO23 +     HV
C   193)  SO21  +    HV  =      SO2  +     HV
C   194)  SO21  +    O2  =      SO3  +     O
C   195)  SO21  +    SO2 =      SO3  +     SO
C   196)  SO23  +    M   =      SO2  +     M
C   197)  SO23  +    HV  =      SO2  +     HV
C   198)  SO23  +    SO2 =      SO3  +     SO
C   199)  SO    +    NO2 =      SO2  +     NO
C   200)  SO    +    O3  =      SO2  +     O2
C   201)  SO2   +    HO2 =      SO3  +     OH
C   202)  HS    +    O3  =      HSO  +     O2
C   203)  HS    +    NO2 =      HSO  +     NO
C   204)  S     +    O3  =      SO   +     O2
C   205)  SO    +    SO  =      SO2  +     S
C   206)  SO3   +    SO  =      SO2  +     SO2
C   207)  S     +    CO2 =      SO   +     CO
C   208)  SO    +    HO2 =      HSO  +     O2
C   209)  SO    +    HCO =      HSO  +     CO
C   210)  H     +    SO  =      HSO
C   211)  HSO   +    HV  =      HS   +     O
C   212)  HSO   +    OH  =      H2O  +     SO
C   213)  HSO   +    H   =      HS   +     OH
C   214)  HSO   +    H   =      H2   +     SO
C   215)  HSO   +    HS  =      H2S  +     SO
C   216)  HSO   +    O   =      OH   +     SO
C   217)  HSO   +    S   =      HS   +     SO
C ********** Additional N2O and O(1D) reactions *********** 
C   218)  N2 + O1D = N2O
C   219)  N2O + H = NO + NO + OH
C   220)  N2O + NO = NO2 + N2
C   221)  O1D + CO2 = CO2 + O
C   222)  CO + O1D = CO2
C
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
C     The INCLUDE files below are used to make parameters available
C     to all the different subroutines. That way they only need to be
C     defined once, which is as it should be.
       INCLUDE 'INCLUDECHEM/parNZ.inc'	
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNEQ_LDA.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
       INCLUDE 'INCLUDECHEM/parNZA.inc'
       INCLUDE 'INCLUDECHEM/parCORRK.inc'
c-AJ 12/21/2022 THE LAST TWO OF THESE INCLUDE FILES CONTAIN PARAMETERS
C      FOR DIMENSIONING THE SOLAR ZENITH ANGLE ARRAY (NZA) AND THE
C      CORRELATED-K TABLE. 
      
C-AJ 03/10/2022 ADD PARAMETER FOR GAUSSIAN QUADRATURE
C-JK  Parameter nrow is the number of entries in the table, which is 
C     found in the file 'gaussian_factors.txt' in subdirectory IO. It
C     is read in on unit 68. The subroutine that reads it is located at
C     the bottom of this main program file.
      PARAMETER(nrow=11)  
    
      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
     2  ,USAVE(NQ,NZ),R(NZ),U(NQ)
      DIMENSION DPU(NZ,3),DPL(NZ,3)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)                                    
      DIMENSION TSAV(NZ),TDUM(NZ),TP1(NZ)
      DIMENSION water(NZ),FLOW(NQT),fluxsave(108),sfxsave(10)

      CHARACTER :: STARR*3,DIRDATA*4, AA*11,DIRIO*2
C-AJ 01/24/2022 ADD VARIABLES FOR CH4 LIFETIME CALCULATION
      REAL:: lsCH4_OH,lsCH4_D,lossCH4,lifetimeCH4,lifetimeCH4_1
C-AJ 03/04/2022 ADD PARAMETER FOR GAUSSIAN INTEGRATION
      REAL,DIMENSION(nrow,20)::xi,wi
      INTEGER,DIMENSION(nrow)::NumGau
      INTEGER::GauFLAG
C-AJ 05/06/2022 ADD OH SOURCE PRINTOUT
      REAL,DIMENSION(NZ):: RateOH1,RateOH25,RateOH59,RatePO3   

C     The following INCLUDE files contain COMMON blocks used
C     to transfer information between subroutines and the main code.
      INCLUDE 'INCLUDECHEM/comSTR.inc'   	!Name of the star
      INCLUDE 'INCLUDECHEM/comFLUXPHOTO.inc'
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
      INCLUDE 'INCLUDECHEM/comSBLOKnew.inc'	!Parameter for Gaussian integration
      INCLUDE 'INCLUDECHEM/comCORRK.inc'	!Correlated-k tables

C-AJ 01/04/2022 change ALPHA(17,4) TO ALPHAOLD(17,4) in comQBLOK.inc
C ALPHA contains the weights for 12-pt Correlated k table
       data ALPHA/0.173927422568727, 0.326072577431273,
     & 0.326072577431273, 0.173927422568727/

C     The following is an important data statement! Species in the code are
C     referred to by their name starting with an 'L'. This list needs to
C     match the species list stored in the vector ISPEC below. The list
C     of L names must be in the same order as ISPEC.
        DATA LH2CO,LO,LH2O,LOH,LHO2,LH2O2,LO3,LH,LH2,LCH4,LCO,
     2  LCH3OOH,LCH3O2,LN2O,LNO,LNO2,LHNO2,LHNO3,LHO2NO2,LNO3,LN2O5,
     3  LCL2O2,LCH3CL,LHOCL,LCL,LCLO,LHCL,LCLONO2,LO2,LH2S,LHS,LSO,
     4  LSO2,LH2SO4,LHSO,LCO2,LSO4AER,LCH21,LCH23,LO1D,LCH3,LH3CO,
     5  LHCO,LN,LNOCL,LCLONO,LCLO2,LCL2,LS,LSO21,LSO23,LHSO3,LSO3,
     6  LS2,LN2/
     7  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     8  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     9  44,45,46,47,48,49,50,51,52,53,54,55/

C   NO PREDISSOCIATION COEFFICIENTS (ALLEN AND FREDERICK, 1982)
C   Note: No one has looked at these in over 40 years. They may be out
C   of date. WACCM probably does a better job on this. We normally do
C   not care much about NO predissociation, which happens high in the
C   atmosphere.
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
      DATA NUML,NUMP/NSP*0,NSP*0/
C
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
C-AP  I have changed the Henry constants to be similar to those in the 
C     (low-O2) Archean code. This has been done only for the sulfur
C     species. So, the H values for those species are not equal to the
C     actual Henry's Law coefficients. They are 'effective solubilities'.

C  Temperature from the US Standard Atmosphere 1976. Used when the 
C  code is not coupled to the climate model.
C  JK   Data are estimated above 64 km. I'm adjusting the temperature
C       near the tropopause downward in order to get better 
C       statospheric H2O
C replaced by Akhil with the WACCM 0.1pc grid
      DATA T/284.50, 281.29, 275.47, 270.39, 265.00,
     &       258.90, 252.21, 245.06, 237.71, 230.33,
     &       223.25, 216.83, 211.39, 206.96, 203.14,
     &       199.57, 196.16, 193.07, 190.53, 188.65,
     &       187.28, 186.19, 185.22, 184.30, 183.39,
     &       182.51, 181.67, 180.92, 180.30, 179.86,
     &       179.69, 179.75, 179.97, 180.19, 180.23,
     &       180.01, 179.41, 178.53, 177.49, 176.38,
     &       175.35, 174.49, 173.83, 173.38, 173.13,
     &       172.99, 172.92, 172.87, 172.81, 172.71,
     &       172.57, 172.37, 172.09, 171.76, 171.35,
     &       170.89, 170.38, 169.83, 169.28, 168.75,
     &       168.25, 167.79, 167.44, 167.16, 166.93,
     &       166.91, 167.05, 167.19, 167.66, 168.48,
     &       169.30, 170.49, 172.19, 174.07, 176.27,
     &       179.12, 182.34, 185.89, 190.33, 195.37,
     &       200.66, 206.69, 213.20, 219.81, 226.72,
     &       233.85, 241.96, 253.23, 266.67, 282.17,
     &       305.06, 333.27, 362.66, 396.97, 434.94,
     &       472.85, 510.88, 548.97, 586.55, 620.44/
C also replaced by Akhil
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
      DATA TP2/287.71, 281.29, 275.47, 270.39, 265.00,
     &       258.90, 252.21, 245.06, 237.71, 230.33,
     &       223.25, 216.83, 211.39, 206.96, 203.14,
     &       199.57, 196.16, 193.07, 190.53, 188.65,
     &       187.28, 186.19, 185.22, 184.30, 183.39,
     &       182.51, 181.67, 180.92, 180.30, 179.86,
     &       179.69, 179.75, 179.97, 180.19, 180.23,
     &       180.01, 179.41, 178.53, 177.49, 176.38,
     &       175.35, 174.49, 173.83, 173.38, 173.13,
     &       172.99, 172.92, 172.87, 172.81, 172.71,
     &       172.57, 172.37, 172.09, 171.76, 171.35,
     &       170.89, 170.38, 169.83, 169.28, 168.75,
     &       168.25, 167.79, 167.44, 167.16, 166.93,
     &       166.91, 167.05, 167.19, 167.66, 168.48,
     &       169.30, 170.49, 172.19, 174.07, 176.27,
     &       179.12, 182.34, 185.89, 190.33, 195.37,
     &       200.66, 206.69, 213.20, 219.81, 226.72,
     &       233.85, 241.96, 253.23, 266.67, 282.17,
     &       305.06, 333.27, 362.66, 396.97, 434.94,
     &       472.85, 510.88, 548.97, 586.55, 620.44/

C   Eddy diffusion profile from Massie & Hunten, JGR (1981)
      DATA EDD/64*0., 5.174E+05, 5.674E+05, 6.224E+05, 6.826E+05,
     2  7.487E+05, 8.212E+05, 9.007E+05, 9.879E+05, 28*1.E6/
 
C =================== Boundary conditions ==================

C ***** UPPER BOUNDARY CONDITIONS *****
      DATA MBOUND/0, 1, 8*0, 1, 26*0/
C   0 = CONSTANT EFFUSION VELOCITY (VEFF)
C   1 = CONSTANT FLUX (SMFLUX)
C
C ***** EFFUSION VELOCITIES *****
      DATA VEFF/NQ*0./
C
C ***** UPPER BOUNDARY FLUXES *****
      DATA SMFLUX/NQ*0./
C
C
C ***** LOWER BOUNDARY CONDITIONS (NQT)*****
c NOTE: Mixing ratios for present Earth for H2,CO,CH4, N2O and CH3Cl
c are defined after the model parameters.

c fixed mixing ratios         H2 CH4 CO     N2O    ch3cl    o2     co2
      DATA LBOUND/2*0, 1, 5*0, 1, 1, 1, 2*0, 1, 8*0, 1, 5*0,1, 6*0, 1,0/
c fixed surface flux
c       DATA LBOUND/2*0, 1, 5*0, 0, 2, 2, 2*0, 2, 8*0,2,5*0, 1,6*0,1,0/
c
C   0 = CONSTANT DEPOSITION VELOCITY (VDEP)
C   1 = CONSTANT MIXING RATIO
C   2 = CONSTANT UPWARD FLUX (SGFLUX)
C
C ***** DEPOSITION VELOCITIES (NQ)*****
      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1.,4.22e-4,2*0.,0.2,1.,
C               h2co o  h2o  oh  ho2 h2o2 o3  h    H2        ch3ooh ch3o2
     2  0., 2*0.2, 6*0.5,0.,0.5,1.,3*0.5, 1.e-4,0.02, 1., 3.E-4, 3*1.,
c      n2o                                o2    h2s   hs   so         
     3  0./
C       co2

C ***** LOWER BOUNDARY FLUXES *****
c NOTE: SGFLUX is redefined later in the main code when the O2 
c mixing ratio is less than 0.21 or when the star is not the Sun
      DATA SGFLUX/NQ*0./
C
C Modify Henry constants for sensitivity run comparison to Archean atmosphere
      H(LHO2) = 9.E3
      H(LH2O2) = 6.2E5
      H(LH2CO) =  4.25E4

C============== FILE SECTION =============================

      DIRDATA = 'DATA'
      DIRIO = 'IO'

c **** INPUT FILES *****
c
C-AJ 12/21/2022 read correlated k table
C-AJ 12/21/2022 The data will be read in subroutine srbands.f under CHEM
c      OPEN(UNIT=10, file =  DIRDATA//'/k12Gauss.dat', status='old')
c      OPEN(UNIT=10, file =  DIRDATA//'/k4_final0315.dat', status='old')
      OPEN(UNIT=10, file =  DIRDATA//'/k4_0927.dat', status='old')
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

C ***** OUTPUT FILES *****
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

c Written in subroutine TWOSTR
c      OPEN(unit=82,file= DIRIO//'/wave_means.dat')

c This file contains the altitude in cm vs the ozone number density (cm^-3)
c      OPEN(unit=83,file= DIRIO//'/O3numdens.out')

c  This is an output file containing altitude, H2O, O3.
c  To be used as input of the climate model, formerly called Pass2SurfMP.dat
      OPEN(unit=84,file= DIRIO//'/fromPhoto2Clima.dat')

C These OUTPUT files are opened along the program 
c In may 9 2019, Jl add the reaction rate table to the code to make life easier 
      OPEN(UNIT=15,file=DIRIO//'/int.rates.out.dat')
      OPEN(UNIT=16,file=DIRIO//'/int.iso_rates.out.dat')
c This file contains total hydrogen mixing ratio vs. altitude (cm)
      OPEN(unit=86,file= DIRIO//'/total_hydrogen.out')
      OPEN(unit=87,file= DIRIO//'/diffusion_coefs.out')
C AJ 05/08/2022 ADD OUTPUT FOR SOURCE OF OH
      open(UNIT=91,file=DIRIO//'/OH_rates.out')
C AJ 08/27/2022 ADD OUTPUT FOR PLOTS
      open(UNIT=92,file=DIRIO//'/OUTPUT_PLOT.dat')  
c in may 16 2019, JL add a clean printout for h2so4 to see possible disarray of input 
      open(UNIT=158, file=DIRIO//'/h2so4.out.dat') 

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

C  Unit 65 is the file 'input_atmchem.dat' in subdirectory IO
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
c  Unit 66 is the file 'planet.dat' in subdirectory IO
      READ(66,502) G,FSCALE,ALB,DELZ,ZTROP,JTROP
 502  FORMAT(F5.1/,F4.2/,F5.3/,E5.1/,E5.1/,I2)
      close(66)
C
c***** READ THE ATMOSPHERIC COMPOSITION AND AEROSOL PARAMETERS ***********
c   This is what we call the restart file. Unit 67 is the file 'atm_composition.dat'
c   in subdirectory IO. Note that the temperature and eddy diffusion profiles are
c   redefined directly below.

      READ(67,400) USOL,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 400  FORMAT(1P8E12.5)
      close(67)

C-PL CHANGE O2 RATIO HERE TO MAKE SURE FO2 GET THE RIGHT NUMBER 05/24/2019
       DO I=1,NZ
       USOL(LO2,I)=FO2 * USOL(LO2,I)/USOL(LO2,1)
       USOL(LCO2,I)=FCO2 * USOL(LCO2,I)/USOL(LCO2,1)
       O2(I) = USOL(LO2,I)
       CO2(I) = USOL(LCO2,I)
       END DO  

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

c** Read the eddy diffusion profile. Unit 64 is the file 'eddy.pdat'
c   in subdirectory DATA
	do J=1,NZ
	 READ(64,245) EDD(J),Z(J)
	end do 
 245  FORMAT(2(1PE11.3, 1X))
      close(64)

C ================ Used when coupling with the climate code ===============
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
C =================================================================

C ========= Lower boundary conditions for biogenic trace gases ====
!changed by Akhil
       LBOUND(9)  = 1 ! H2 (mixing ratio)
       LBOUND(10) = 1 ! CH4 (mixing ratio)
       LBOUND(14) = 1 ! N2O (mixing ratio)
       LBOUND(23) = 1 ! CH3Cl (mixing ratio)
       LBOUND(1)  = 2 ! CH2O (flux)
       LBOUND(11) = 2 ! CO (flux)
       LBOUND(15) = 2 ! NO (flux)
       LBOUND(33) = 2 ! SO2 (flux)


       if(INIT.eq.1) goto 77
c      Surface mixing ratios
         if(LBOUND(9).eq.1)FH2 = 5E-7
         if(LBOUND(10).eq.1)FCH4 = 0.808E-6
         if(LBOUND(14).eq.1)FN2O = 2.73E-7
	
C AJ 08/30/2022 ADD CL BACK
         if(LBOUND(23).eq.1)FCH3CL = 4.57E-10


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

C-AJ after using the correlated-k coefficient
         ! if(LBOUND(9).eq.2)SGFLUX(9) = -5.57E+09  H2
         ! if(LBOUND(10).eq.2)SGFLUX(10) = 8.54E+10 CH4
         ! if(LBOUND(11).eq.2)SGFLUX(11) = 1.76E+11 CO
         ! if(LBOUND(14).eq.2)SGFLUX(14) = 1.18E+09 N2O
         ! if(LBOUND(23).eq.2)SGFLUX(23) = 2.37E+08 CH3Cl

C      Akhil - WACCM6 CMIP6 emission fluxes
         if(LBOUND(1).eq.2)SGFLUX(1)   = 5.16E+08  ! CH2O
         if(LBOUND(11).eq.2)SGFLUX(11) = 4.98E+10  ! CO
         if(LBOUND(15).eq.2)SGFLUX(15) = 2.37E+09  ! NO
         if(LBOUND(33).eq.2)SGFLUX(33) = 1.2E+08  ! SO2
         if(LBOUND(23).eq.2)SGFLUX(23) = 2.37E+08  ! Default ch3cl

C ============ End biogenic trace gas section ========================

C-AJ 03/10/2022 Read the Gaussian points and weights for global integration
      call gaussian_data(xi,wi,NumGau)
C
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

C-AJ 09/03/2022 Test model sensitivity to CL
        IF(CLFLAG.EQ.0) THEN
         SGFLUX(23) = SGFLUX(23)
         WRITE(90,777)'Standard CH3CL flux = ',SGFLUX(23)
 777     FORMAT(A23,1PE10.3)
        ELSE
         SGFLUX(23) = SGFLUX(23)/30
         WRITE(90,777)'Low CH3CL flux      = ', SGFLUX(23)
        END IF

C ===================== DEFINE THE SPECIES LIST =================

C  Long-lived species
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
C
C   TRIDIAGONAL SOLVER
       ISPEC(37) = 6HSO4AER
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
C   Inert species
       ISPEC(54) = 2HS2
       ISPEC(55) = 2HN2
       ISPEC(56) = 2HHV
       ISPEC(57) = 1HM

      do i=1,NSP
       NUML(i) = 0
       NUMP(i) = 0
      enddo
C
C ================= Set up the chemical reaction network =============

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
      IF(JCHEM(M,J).EQ.1H ) GO TO 5
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
      DO 7 M=1,2
      N = 3-M

      DO 7 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
      NUML(I) = NUML(I) + 1
      IF(NUML(I).GT.NMAX) GO TO 20
      K = NUML(I)
      ILOSS(1,I,K) = J
      ILOSS(2,I,K) = JCHEM(N,J)
   7  CONTINUE

      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE
C ============= End chemical reaction network ===============

      LTIMES = 0       !Counter for the photorate subroutine
      ISULF = 0
      VOLFLX = 3.E9    !Volcanic flux
c      print *, 'Z before call',Z

      CALL GRID
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

      CALL OUTPUTP(0,NSTEPS,0.,FLOW)

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
 
C ================== START THE TIME-STEPPING LOOP ============
      DO 1 N=1,NSTEPS
      TIME = TIME + DT
      NN = NN + 1
      MS = (N-1)/MP
      SM = (N-1)/PM
      IF(NN.EQ.NSTEPS) SM = MS
      IF(SM-MS.GT.0.01) GO TO 18
      IF(N.GT.1 .AND. TIME.LT.1.E4) GO TO 18
C
C ************ Start the photolysis calculation *************

C   STORE ABSORBERS USED TO BLOCK OUT SOLAR UV RADIATION
      DO 35 I=1,NZ
      H2O(I) = ABS(USOL(LH2O,I))
      O3(I) = ABS(USOL(LO3,I))
      O2(I) = USOL(LO2,I)
      CO2(I) = USOL(LCO2,I)
      FSO2(I) = ABS(USOL(LSO2,I))
      H2S(I) = ABS(USOL(LH2S,I))
      CH4(I) = ABS(USOL(LCH4,I))
  35  CONTINUE
            
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
  
C-AJ  03/04/2022 CODE ADDED FOR GAUSSIAN INTEGRATION
C-JK  Code modified by Jim (11/14/2023)

C   Single solar zenith angle case (ZY and AGL are already specified)
      IF(GauFlag.eq.0) GO TO 135
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)
      GO TO 139

C   Gaussian quadrature case (GauFLAG = 1)
 135  CONTINUE
C-JK  NZA is the number of solar zenith angles for Gaussian
C     integration. One can change this by editing the 'parNZA.inc'
C     file in subdirectory INCLUDE. AGL is always 0.5 in this case.
      AGL = 0.5

c-AJ 03/10/2022 Find the right row in the matrix
      DO i = 1,11
      isave=i
      if(NumGau(i).EQ.NZA) exit
      END DO     

C   Loop over solar zenith angles
      DO 136 IPHOT=1,NZA
      radians=xi(isave,IPHOT)
      IF(GauFLAG.eq.0) ZY=acos(radians)*180./3.14159
      WEIGHT=wi(isave,IPHOT)

      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)

C   Copy the photolysis rate vectors into a photolysis matrix     
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
      PCH3CLM(J,IPHOT) = PCH3CL(J)*WEIGHT
      PCL2M(J,IPHOT) = PCL2(J)*WEIGHT
      PCLO2M(J,IPHOT) = PCLO2(J)*WEIGHT
      PHCLM(J,IPHOT) = PHCL(J)*WEIGHT
      PHOCLM(J,IPHOT) = PHOCL(J)*WEIGHT
      PNOCLM(J,IPHOT) = PNOCL(J)*WEIGHT
      PCLONOM(J,IPHOT) = PCLONO(J)*WEIGHT
      PCLONO2M(J,IPHOT) = PCLONO2(J)*WEIGHT
      PCL2O2M(J,IPHOT) = PCL2O2(J)*WEIGHT
 137  CONTINUE   ! End the loop over altitude
 136  CONTINUE   ! End the loop over zenith angles

C-JK  Sum the photolysis rates from quadrature, and recopy the result
C     into the original species photolysis rate vectors, making it look
C     a lot like the old code.

      DO 138 J=1,NZ
      PO2D(J) = SUM(PO2DM(J,1:NZA))
      PO2(J) = SUM(PO2M(J,1:NZA))
      PH2O(J) = SUM(PH2OM(J,1:NZA))
      PO3D(J) = SUM(PO3DM(J,1:NZA))
      PO3(J) = SUM(PO3M(J,1:NZA))
      PH2O2(J) = SUM(PH2O2M(J,1:NZA))
      PCO2(J) = SUM(PCO2M(J,1:NZA))
      PH2(J) = SUM(PH2M(J,1:NZA))
      PHCO(J) = SUM(PHCOM(J,1:NZA))
      PCO2D(J) = SUM(PCO2DM(J,1:NZA))
      PHO2(J) = SUM(PHO2M(J,1:NZA))
      PCH4(J) = SUM(PCH4M(J,1:NZA))
      PMOOH(J) = SUM(PMOOHM(J,1:NZA))
      PN2O(J) = SUM(PN2OM(J,1:NZA))
      PHNO3(J) = SUM(PHNO3M(J,1:NZA))
      PNO(J) = SUM(PNOM(J,1:NZA))
      PNO2(J) = SUM(PNO2M(J,1:NZA))
      PHNO4(J) = SUM(PHNO4M(J,1:NZA))
      PCH3CL(J) = SUM(PCH3CLM(J,1:NZA))
      PCL2(J) = SUM(PCL2M(J,1:NZA))
      PCLO2(J) = SUM(PCLO2M(J,1:NZA))
      PHCL(J) = SUM(PHCLM(J,1:NZA))
      PHOCL(J) = SUM(PHOCLM(J,1:NZA))
      PNOCL(J) = SUM(PNOCLM(J,1:NZA))
      PCLONO(J) = SUM(PCLONOM(J,1:NZA))
      PCLONO2(J) = SUM(PCLONO2M(J,1:NZA))
      PCL2O2(J) = SUM(PCL2O2M(J,1:NZA))
      PNO3(J) = SUM(PNO3M(J,1:NZA))
      PN2O5(J) = SUM(PN2O5M(J,1:NZA))
      PSO2(J) = SUM(PSO2M(J,1:NZA))
      PH2S(J) = SUM(PH2SM(J,1:NZA))
      PSO21(J) = SUM(PSO21M(J,1:NZA))
      PSO23(J) = SUM(PSO23M(J,1:NZA))
      PH2SO4(J) = SUM(PH2SO4M(J,1:NZA))
      PHO2(J) = SUM(PHO2M(J,1:NZA))
 138  CONTINUE

C   Copy the photolysis rates into the reaction rate matrix
 139  CONTINUE
      DO 141 J=1,NZ
      AR(23,J) = PO2D(J)
      AR(24,J) = PO2(J)
      AR(25,J) = PH2O(J)
      AR(26,J) = PO3D(J)
      AR(27,J) = PO3(J)
      AR(28,J) = PH2O2(J)
      AR(29,J) = PCO2(J)
      AR(38,J) = PH2(J)
      AR(39,J) = PHCO(J)
      AR(42,J) = PCO2D(J)
      AR(50,J) = PHO2(J)
      AR(51,J) = PCH4(J)
      AR(52,J) = PMOOH(J)
      AR(53,J) = PN2O(J)
      AR(54,J) = 1.7E-3
      AR(55,J) = PHNO3(J)
      AR(56,J) = PNO(J)
      AR(57,J) = PNO2(J)
      AR(95,J) = PHNO4(J)
      AR(101,J) = PCH3CL(J)
      AR(132,J) = PCL2(J)
      AR(133,J) = PCLO2(J)
      AR(134,J) = PHCL(J)
      AR(135,J) = PHOCL(J)
      AR(136,J) = PNOCL(J)
      AR(137,J) = PCLONO(J)
      AR(138,J) = PCLONO2(J)
      AR(142,J) = PCL2O2(J)
      AR(151,J) = PNO3(J)
      AR(153,J) = PN2O5(J)
      AR(156,J) = 0.
      AR(157,J) = 0.7*PSO2(J)
      AR(158,J) = PH2S(J)
      AR(182,J) = 0.3*PSO2(J)
      AR(186,J) = PSO21(J)
      AR(187,J) = PSO23(J)
      AR(188,J) = PH2SO4(J)
      AR(189,J) = 0.
      AR(211,J) = PHO2(J)
 141  CONTINUE

      CALL AERCON(H2O)
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
      VO2 = (PO2(NZ) + PO2D(NZ)) * HA
      VCO2 = (PCO2(NZ) + PCO2D(NZ)) * HA
      VEFF(LO2) = VO2
      VEFF(LCO2) = VCO2
      SMFLUX(LCO) = - VCO2*CO2(NZ)*DEN(NZ)   
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ) - 2.*VO2*O2(NZ)*DEN(NZ)
     2  + 2.*VEFF(LCH4)*USOL(LCH4,NZ)*DEN(NZ)
C
      NMP = NSTEPS - MP
      IF (nn.LT.NSTEPS) GO TO 18

C   Print out the photolysis rates
      write(90,97)
  97  FORMAT(//1X,'PHOTOLYSIS RATES')
      write(90,98) 
  98  FORMAT(/5X,'Z',7X,'PO2',6X,'PO2D',5X,'PCO2',5X,'PCO2D',4X,
     2  'PH2O',5X,'PO3',6X,'PO3D',5X,'PH2O2',4X,'PHCO',5X,'PH2',
     3  6X,'PHO2')
      write(90,99)(Z(I),PO2(I),PO2D(I),PCO2(I),PCO2D(I),PH2O(I),
     2  PO3(I),PO3D(I),PH2O2(I),PHCO(I),PH2(I),PHO2(I),I=1,NZ,
     3  3)
  99  FORMAT(2X,1P12E9.2)
      write(90,198)
 198  FORMAT(/5X,'Z',6X,'PCH4',5X,'PCH3OOH',2X,'PN2O',5X,'PHNO3',4X,
     2  'PNO',6X,'PNO2',5X,'PHNO4',4X,'PCCL3F',3X,'PCCL2F2',2X,
     3  'PCCL4',4X,'PCH3CL')
      write(90,99) (Z(I),PCH4(I),PMOOH(I),PN2O(I),PHNO3(I),PNO(I),
     2  PNO2(I),PHNO4(I),PCCL3F(I),PCCL2F2(I),PCCL4(I),PCH3CL(I),
     3  I=1,NZ,3)
      write(90,197)
 197  FORMAT(/5X,'Z',6X,'PMCCL3',3X,'PCL2',5X,'PHOCL',4X,'PNOCL',4X,
     2  'PCLONO',3X,'PCLONO2',2X,'PCLO2',4X,'PHCL')
      write(90,199) (Z(I),PMCCL3(I),PCL2(I),PHOCL(I),PNOCL(I),PCLONO(I),
     2  PCLONO2(I),PCLO2(I),PHCL(I),I=1,NZ,3)
 199  FORMAT(2X,1P9E9.2)
      write(90,298)
 298  format(/5x,'Z',6x,'PNO3',5x,'PN2O5',4x,'PCL2O2',3x,'PSO2',5x,
     2  'PH2S',5x,'PSO21',4x,'PSO23',4x,'PH2SO4')
      write(90,299) (z(i),pno3(i),pn2o5(i),pcl2o2(i),pso2(i),ph2s(i),
     2  pso21(i),pso23(i),ph2so4(i),i=1,nz,3)
 299  format(2x,1p9e9.2)
      write(90,299)
  18  CONTINUE
C   End photolysis rate printout

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
      CALL DOCHEM(FVAL,IDO)

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

C  Ignore certain species at certain altitudes when computing
C  EMAX so that they don't limit the time step. (This section
C  of code represents an art, rather than a science. Try not
C  to overuse this. --JK)

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
      IF((I.EQ.LCLONO2).AND.(Z(J).GT.1.E6)) GOTO 26
      IF((I.EQ.LCH3OOH).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LHOCL).AND.(Z(J).GT.1.0E6)) GOTO 26

      
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

C-PL  DO NOT LET [H2O] CHANGE AFTER CHEMICAL REACTION
C
      DO 4 J=1,NZ
      IF(Z(J).LT.ZTROP) USOL(3,J) = H2O(J)
   4  CONTINUE

C-AP Adding tridiagonal solver
C       TRIDIAGONAL INVERSION *****
      L=1   ! Number of aerosol species considered
      I = NQ + L
      IF(I.EQ.LSO4AER) MZ = 50
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

      IF (NFLAG.NE.0) write(90,401) N,NFLAG,I
 401  FORMAT(//1X,'TRIDIAGONAL SOLVER FAILED AT N =',I3,2X,
     2  'NFLAG =',I2,2X,'SPECIES #',I2)
C
      IF(I.EQ.LSO4AER) THEN
        DO 59 J=1,MZ
   59     SO4AER(J) = TY(J)
      ENDIF
C
C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION
      IF(I.EQ.LSO4AER .AND. MZ.NE.NZ) THEN
        DO 61 J=MZP1,NZ
          SO4AER(J) = SO4AER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
   61     SO4AER(J) = AMAX1(SO4AER(J),1E-100)
      ENDIF

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

c Adjusting time step to assure that TIME <= TSTOP
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
C
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
      write(90,*)'BEFORE CALL TO OUTPUT'
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
     
C ================== END THE TIME-STEPPING LOOP ==================

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

      GO TO 21
  20  write(90,300)I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
      GO TO 21
  25  write(90,301)IERR
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)
C
  21  CONTINUE
c      close(82)

c-AJ 05/06/2022 calculate reaction rate of OH
      DO 377 J = 1,NZ
      RateOH1(J) = AR(1,J)*SL(LO1D,J)*USOL(LH2O,J)*DEN(J)
      RateOH25(J) = AR(25,J)*USOL(LH2O,J)*DEN(J)
      RateOH59(J) = AR(59,J)*USOL(LCH4,J)*DEN(J)*SL(LO1D,J)
      RatePO3(J) = AR(26,J)*USOL(LO3,J)*DEN(J)
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
     3 'FO3',6x,'FH2O',6x,'NumO2',6x,'NumOH',6x,'NumO1D',5x,'PN2O',
     4 5x,'PCH3CL')
      write(92,388) (z(I),PRESS(I),DEN(I),RateOH1(I),RateOH25(I),
     2 RatePO3(I),AR(23,I),AR(24,I),AR(25,I),AR(26,I),AR(28,I),
     3 USOL(LO3,I),USOL(LH2O,I),SL(LO2,I),SL(LOH,I),SL(LO1D,I),
     4 AR(53,I),AR(101,I),I=1,NZ)
 388  format(2x,1p18e10.3)     

C-AJ 01/24/2022 CALCULATE LIFETIME OF CH4 IN THE SURFACE
      lsCH4_OH = 4.16E-13*((T(1)/298.)**2.18)*EXP(-1230./T(1))*SL(4,1)
      lsCH4_D = 1.3E-10*SL(33,1) + 7.51E-12*SL(33,1)
      lossCH4 = lsCH4_OH + lsCH4_D
      lifetimeCH4 = 1./lossCH4*1./(3600.*24.*365.)
      lifetimeCH4_1 = SL(10,1)*7.5E5/FLOW(LCH4)*1./(3600.*24.*365.)

      print *,'Temperature z1 z51 z91= ', T(1),T(51),T(91)
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
      
      PRINT*, "Success!"
 1039 FORMAT('[',F4.1,' , ',1PE10.2,' , ',1PE10.2,' , ',1PE10.2,']')
      m_S = 2.*32.07/6.022E23 !g per S atom
      c_S = ((8.*1.38E-16*T(J))/(3.1415*2.*m_S))**0.5
      S2_lifetime = c_S*3.1415*RPAR(10,1)*AERSOL(10,1)

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
 105  FORMAT(//1X,'MIXING RATIOS OF LONG-LIVED SPECIES')
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
 106  FORMAT(/1X,'MIXING RATIO OF AEROSOL')
      write(90,185)
 185  FORMAT(/5X,'Z',6X,'SO4AER')
      DO 18 J=1,3
  18  write(90,182) Z(J),SO4AER(J)
      DO 19 J=4,NZ,ISKIP
  19  write(90,182) Z(J),SO4AER(J)
 182  FORMAT(1X,1P2E9.2)
      write(90,182)
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
 151  FORMAT(/1X,'JCOLD =',I3,' FH2O AT COLD TRAP =',1PE10.3)
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
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES')
      DO 9 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1+IROW - 1
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 22 I=1,NZ,ISKIP
       write(90,120) Z(I),(SL(K,I),K=K1,K2)
  22  CONTINUE
   9  CONTINUE
 355  FORMAT(1PE10.3,1PE12.3)

      write(90,185)
      DO J=1,3
      write(90,182) Z(J),SL(LSO4AER,J)
      END DO

      DO J=4,NZ,ISKIP
      write(90,182) Z(J),SL(LSO4AER,J)
      END DO

      write(90,155)
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED SPECIES')
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
