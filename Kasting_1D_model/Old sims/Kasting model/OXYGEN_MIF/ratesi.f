      SUBROUTINE RATESI
!      INCLUDE 'parameters_add_CO2.txt'
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
c      INCLUDE '../INCLUDECHEM/parNEQ_LDA.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
      INCLUDE '../INCLUDECHEM/parNF.inc'  ! need to check later

      INCLUDE '../INCLUDECHEM/comABLOK.inc' !NEW:JTROP DEL:PMOD
      INCLUDE '../INCLUDECHEM/comBBLOK.inc' !A CONVERT TO AR
      INCLUDE '../INCLUDECHEM/comGBLOK.inc' !NEW:H(NQ)
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRBLOK.inc'
      INCLUDE '../INCLUDECHEM/comSULBLK.inc'
     
C      PARAMETER(NQI=31, NSLS=38, NMAXI=70)
C      PARAMETER(NRI=353, NSPI=78, NSPI1=NSPI+1, NSPI2=NSPI+2)
      PARAMETER(NQI=29, NSLS=36, NMAXI=70)
      PARAMETER(NRI=353, NSPI=76, NSPI1=NSPI+1, NSPI2=NSPI+2)

!      COMMON/RBLOK/A(NR,NZ),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
!     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
!      COMMON/GBLOK/RAIN(NZ),FSAT(NZ),RAINGC(NQ,NZ)
      COMMON/RATESIS/RAINGCI(NQI,NZ),DDI(NQI,NZ),DKI(NQI,NZ),
     2  DUI(NQI,NZ),DLI(NQI,NZ),DHUI(NQI,NZ),DHLI(NQI,NZ),HIZ(NQI,NZ),
     3 VDEPI(NQI),RATI(NRI)

!      COMMON/ABLOK/EDD(NZ),DEN(NZ),DK(NQT,NZ),Z(NZ),T(NZ),G,FSCALE,
!     2  ALB,DELZ,BOVERH,DD(NQT,NZ),DL(NQT,NZ),DU(NQT,NZ),DI(NQT,NZ),
!     3  HI(NQT,NZ),DHU(NQT,NZ),DHL(NQT,NZ),PMOD 

      COMMON/ISOTOP/AI(NRI,NZ),ILOSSI(2,NSPI,NMAXI),
     2  JCHEMI(5,NRI),NUMLI(NSPI),NUMPI(NSPI),TPI(NQI),TLI(NQI),
     3  YPI(NQI,NZ),YLI(NQI,NZ),SRI(NQI),TLOSSI(NQI),PHIDEPI(NQI),
     4 ISPECI(NSPI2),DIZ(NSPI2,NZ),IPRODI(NSPI,NMAXI),PSO4AER
     5 ,LBOUNDI(NQI),FLOWI(NQI),FUPI(NQI),CONI(NQI)

C      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
C     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
C     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,
C     4  LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
C     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
C     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
C     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
C     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
C     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2

      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ,
     4  LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,LIS,
     9  LISO21,LISO23,LIHSO3,LISO3,LIN2
      
      
C ***** Photolysis reactions *****
      DO 1 I=1,NZ
      AI(45,I) = AR(26,I)/2.  ! half O2Q + HV = OQ + O(1D)
      AI(46,I) = AR(26,I)/2. ! hlf O2Q + HV = O2 + Q(1D)
      AI(47,I) = AR(26,I) ! same OQO + HV = OQ + O(1D)

      AI(48,I) = AR(27,I)! same OQO + HV = OQ + O(3P)
      AI(49,I) = AR(27,I)/2.! half O2Q + HV = OQ + O(3P)
      AI(50,I) = AR(27,I)/2.! half O2Q + HV = O2 + Q(3P)

      AI(51,I) = AR(29,I)/2. ! half COQ + HV = CQ + O(3P)
      AI(52,I) = AR(29,I)/2. ! half COQ + HV = CO + Q(3P)

      AI(67,I) = AR(42,I)/2.! half COQ + HV = CQ + O(1D)
      AI(68,I) = AR(42,I)/2.! half COQ + HV = CO + Q(1D)

      AI(279,I)= AR(23,I)/2. !half OQ  +HV =  Q  + O(1D)
      AI(280,I)= AR(23,I)/2. !half OQ  +HV =  O  + Q(1D)

      AI(281,I)= AR(24,I)    ! OQ  + HV =  O  + Q QUESTION      

      AI(282,I)= AR(25,I)    ! H2Q + HV =  H  + OH

      AI(283,I)= AR(28,I)    ! H2OQ+ HV = QH + OH

      AI(284,I)= AR(38,I)    ! H2CQ + HV = H2 + CQ

      AI(285,I)= AR(39,I)    ! H2CQ + HV = HCQ + H

      AI(286,I)= AR(40,I)    ! HCQ + HV = CQ + H
C AJ 01/26
C      AI(287,I)= 0. !AR(50,I)/2. ! half HOQ + HV = QH + O 
C      AI(288,I)= AR(50,I) !/2. ! half HOQ + HV = QH + O

C AJ 02/06 DEBUG
      AI(287,I)= AR(50,I)/2. ! half HOQ + HV = QH + O 
      AI(288,I)= AR(50,I)/2. ! half HOQ + HV = QH + O

 
      AI(289,I)= AR(52,I)    ! CH3OQH + HV = H3CO + QH
      AI(290,I)= AR(52,I)    ! CH3QOH + HV = H3CQ + OH

      AI(291,I)= AR(53,I)    !N2Q + HV = N2 + Q 
 
      AI(292,I)= AR(54,I)/2. !half HNOQ  + HV = NQ + OH 
      AI(293,I)= AR(54,I)/2. !half HNOQ  + HV = NO + QH 

      AI(294,I)= AR(55,I)/3.*2. !half HNO2Q  + HV = NOQ + OH 
      AI(295,I)= AR(55,I)/3.*1. !half HNO2Q  + HV = NO2 + QH 

      AI(296,I)= AR(56,I)    ! NQ +HV =N + Q

      AI(297,I)= AR(57,I)/2. !half NOQ + HV = NQ + O
      AI(298,I)= AR(57,I)/2. !half NOQ + HV = NO + Q

      AI(299,I)= AR(105,I)/3.*2. !half NO2Q + HV = NOQ + O
      AI(300,I)= AR(105,I)/3.*1. !half NO2Q + HV = NO2 + Q

      AI(301,I)= AR(110,I)    ! SQ + HV = S + Q

      AI(302,I)= AR(111,I)/2. !half SOQ + HV = SQ + O
      AI(303,I)= AR(111,I)/2. !half SOQ + HV = SO + Q

      AI(304,I)= AR(136,I)    !SOQ + HV = S + OQ

      AI(305,I)= AR(140,I)    !SOQ + HV = SO1Q 
      AI(306,I)= AR(141,I)    !SOQ + HV = SO3Q

      AI(307,I)= AR(142,I)/2.   !H2SO3Q + HV = SOQ + OH + OH
      AI(308,I)= AR(142,I)/2. !H2SO3Q + HV = SO2 + QH + OH

      AI(309,I)= AR(143,I)/2. !half SO2Q + HV = SOQ + O 
      AI(310,I)= AR(143,I)/2. !half SO2Q + HV = SO2 + Q


      AI(311,I)= AR(146,I)    !SO1Q + HV = SO3Q +HV 
      AI(312,I)= AR(147,I)    !SO1Q + HV = SOQ +HV 

      AI(313,I)= AR(151,I)    !SO3Q + HV = SOQ +HV 
      AI(314,I)= AR(165,I)    !HSQ + HV = HS +Q 

      AI(315,I)= AR(95,I) !same HOQNO2 + HV  =  HOQ + NO2 
      AI(316,I)= AR(95,I) !same HO2NOQ + HV  =  HO2 + NOQ

      AI(317,I)= AR(107,I)/2. !half N2O4Q + HV  =  NOQ + NO3
      AI(318,I)= AR(107,I)/2. !half N2O4Q + HV  =  NO2 + NO2Q 
      AI(319,I)= AR(107,I) !same N2QO4 + HV  =  NO2 + NO2Q 
C AJ 01/26 MOVE TWO PHOTOLYSIS REACTION
C      AI(339,I)= AR(50,I) !HQO + HV = QH + O
C      AI(348,I)= AR(95,I) !HQONO2 + HV = HQO + NO2
C AJ 02/06 DEBUG
      AI(339,I)= 0. !HQO + HV = QH + O
      AI(348,I)= 0. !HQONO2 + HV = HQO + NO2
   1  CONTINUE 
C ***** FILL UP RATE MATRIX *****
      DO 4 I=1,NZ
      AI(1,I) = AR(1,I)! same H2Q + O(1D) = OH + QH
      AI(2,I) = AR(1,I)! same H2O + Q(1D) = OH + QH

      AI(3,I) = AR(2,I)! same H2 + Q(1D) = QH + H
      AI(4,I) = AR(3,I)! same H2 + Q = QH + H 
      AI(5,I) = AR(4,I)! same H2 + QH = H2Q + H

      AI(6,I) = AR(5,I)/2. ! half H +O2Q = QH + O2
      AI(7,I) = AR(5,I)/2.  ! half H + O2Q = OH + OQ
      AI(8,I) = AR(5,I) ! same H + OQO = OH + OQ
C AJ LUNAR NEW YEAR'S EVE
C      AI(9,I) = AR(6,I)/2. ! same H + OQ + M = HOQ + M
C AJ 02/06 DEBUG
      AI(9,I) = AR(6,I) ! same H + OQ + M = HOQ + M

      AI(10,I) = AR(7,I) ! same H + HOQ = H2 + OQ
      AI(11,I) = AR(8,I)/2.! half H + HOQ = H2Q + O
      AI(12,I) = AR(8,I)/2.! half H + HOQ = H2O + Q 
      AI(13,I) = AR(9,I)! same  H + HOQ = QH + OH

      AI(14,I) = AR(10,I)! same QH + O = H + OQ
      AI(15,I) = AR(10,I)! same OH + Q = H + OQ

      AI(16,I) = AR(11,I)! same QH + HO2 = H2Q + O2 
      AI(17,I) = AR(11,I)! same QH + HO2 = H2Q + O2 
C AJ LUNAR NEW YEAR'S EVE
C      AI(18,I) = 0. !AR(12,I) ! same QH + O3 = HOQ + O2
C AJ 02/06 DEBUG
      AI(18,I) = AR(12,I) ! same QH + O3 = HOQ + O2

      AI(19,I) = AR(12,I)/2.! half OH + O2Q = HOQ + O2
      AI(20,I) = AR(12,I)/2.! half OH + O2Q = HO2 + OQ
      AI(21,I) = AR(12,I)!same OH + OQO = HO2 + OQ

      AI(22,I) = AR(13,I)! same HOQ + O = OH + OQ 
      AI(23,I) = AR(13,I)! same HO2 + Q = QH + O2

      AI(24,I) = AR(14,I)/2.*0.0! half HOQ + O3 = QH + 2O2 
      AI(25,I) = AR(14,I)!/2.! half HOQ + O3 = OH + O2 + OQ
      AI(26,I) = AR(14,I)/2.! half HO2 + O2Q = QH + O2 + O2
      AI(27,I) = AR(14,I)/2.! half HO2 + O2Q = OH + O2 + OQ
      AI(28,I) = AR(14,I) ! same HO2 + OQO = OH + O2 + OQ
C AJ LUNAR NEW YEAR'S EVE
      AI(29,I) = AR(15,I)/2.*COMBIN ! ACCELERATED HOQ + HO2 = H2OQ + O2
      AI(30,I) = AR(15,I)/2.*COMBIN !undoubled ! ACCELERATED HOQ + HO2 = H2O2 + OQ 
      AI(31,I) = AR(16,I)/2.! same H2OQ + OH = HOQ + H2O
  
      AI(32,I) = AR(16,I)! same H2O2 + QH = HO2 + H2Q

      AI(33,I) = AR(17,I)*2. ! ACCELERATED Q + O  = OQ 
C-PL-JL add symmetry effect 07/2019 Averaged method
      AI(34,I) = AR(18,I)  !*(1.03+0.92)/2.   !* 0.92 *1.2   ! same Q + O2 + M = O2Q + M
      AI(35,I) = AR(18,I)/2.!*(1.292+1.426)/2. !* 1.45 *1.2   ! half O + OQ + M = O2Q + M
      AI(36,I) = AR(18,I)/2.!*(1.005+1.007)/2.!*1.08 ! half O + OQ + M = OQO + M

      AI(37,I) = AR(19,I) ! same Q + O3 = O2 + OQ
      AI(38,I) = AR(19,I) ! same O + OQO = O2 + OQ
      AI(39,I) = AR(19,I) ! same O + O2Q = O2 + OQ

      AI(40,I) = AR(20,I)/2.*COMBIN ! ACCELERATED QH + OH = H2Q + O
      AI(41,I) = AR(20,I)/2.*COMBIN !undoubled ! ACCELERATED QH + OH = H2O + Q 

      AI(42,I) = AR(21,I) !same Q(1D) + N2 = Q + N2

      AI(43,I) = AR(22,I) !same Q(1D) + O2 = Q + O2 
      AI(44,I) = AR(22,I) !same O(1D) + OQ = O + OQ 
c
      AI(53,I) = AR(30,I)! same CQ + OH = COQ + H
      AI(54,I) = AR(30,I)! same CO + QH = COQ + H 

      AI(55,I) = AR(31,I)! same CQ + O + M = COQ + M
      AI(56,I) = AR(31,I)! same CO + Q + M = COQ + M 

      AI(57,I) = AR(32,I)! same H + CQ + M = HCQ + M

      AI(58,I) = AR(33,I)! same H + HCQ = H2 + CQ 

      AI(59,I) = AR(34,I)/2.*COMBIN! ACCELERATED HCQ + HCO = H2CQ + CO
      AI(60,I) = AR(34,I)/2.*COMBIN!undoubled ! ACCELERATED HCQ + HCO = H2CO + CQ

      AI(61,I) = AR(35,I)! same QH + HCO = H2Q + CO
      AI(62,I) = AR(35,I)! same OH + HCQ = H2O + CQ

      AI(63,I) = AR(36,I)! same Q + HCO = H + COQ 
      AI(278,I) = AR(36,I)! same O + HCQ = H + COQ
 
      AI(64,I) = AR(37,I)! same Q + HCO = QH + CO
      AI(65,I) = AR(37,I)! same Q + HCO = OH + CQ

      AI(66,I) = AR(41,I) ! same H2CQ + H = H2 + HCQ 
c
      AI(69,I) = AR(44,I)! same  HCQ + O2 = HO2 + CQ
C AJ LUNAR NEW YEAR'S EVE
C      AI(70,I) = AR(44,I)/2. ! same  HCO + OQ = HOQ + CO
C AJ 02/06 DEBUG
      AI(70,I) = AR(44,I) ! same  HCO + OQ = HOQ + CO

      AI(71,I) = AR(45,I)! same H2CQ + OH = H2O + HCQ 
      AI(72,I) = AR(45,I)! same H2CO + QH = H2Q + HCO 

      AI(73,I) = AR(46,I)! same H + QH + M = H2Q + M

      AI(74,I) = AR(47,I)*COMBIN!undoubled ! ACCELERATED QH + OH + M = H2OQ + M

      AI(75,I) = AR(48,I) ! same H2CQ + O = HCQ + OH
      AI(76,I) = AR(48,I) ! same H2CO + Q = HCO + QH
C AJ LUNAR NEW YEAR'S EVE
C      AI(77,I) = AR(49,I)/2. ! same  H2OQ + O = OH + HOQ
C AJ 02/06 DEBUG
      AI(77,I) = AR(49,I) ! same  H2OQ + O = OH + HOQ

      AI(78,I) = AR(49,I) ! same  H2O2 + Q = QH + HO2

      AI(79,I) = AR(58,I) ! same CH4 + QH  =  CH3 + H2Q
      AI(80,I) = AR(59,I) ! same CH4 + Q(1D)  =  CH3 + QH
      AI(81,I) = AR(60,I) ! same CH4 + Q(1D)  =  H2CQ + H2

      AI(82,I) = AR(62,I)/2.! half 1CH2 + OQ  =  H2CQ + O
      AI(83,I) = AR(62,I)/2.! half 1CH2 + OQ  =  H2CO + Q 

      AI(84,I) = AR(66,I)/2.!half 3CH2 + OQ  =  H2CQ + O
      AI(85,I) = AR(66,I)/2.!half 3CH2 + OQ  =  H2CO + Q

      AI(86,I) = AR(67,I)/2. ! same CH3 + OQ + M  =  CH3OQ + M
      AI(87,I) = AR(67,I)/2. ! same CH3 + OQ + M  =  CH3QO + M

      AI(88,I) = AR(68,I) ! same CH3 + QH  =  H2CQ + H2
      AI(89,I) = AR(69,I) ! same CH3 + Q  =  H2CQ + H

      AI(90,I) = AR(70,I)/2. !half CH3 + O2Q  =  H2CQ + HO2
      AI(91,I) = AR(70,I)/2. !half CH3 + O2Q  =  H2CO + HOQ
C AJ LUNAR NEW YEAR'S EVE
C      AI(92,I) = AR(70,I)/2.    !same CH3 + OQO  =  H2CO + HOQ
C AJ 02/06 DEBUG
      AI(92,I) = AR(70,I)    !same CH3 + OQO  =  H2CO + HOQ

      AI(93,I) = AR(71,I) !same CH3OQ + HO2  =  CH3OQH + O2   
      AI(94,I) = AR(71,I) !same CH3QO + HO2  =  CH3QOH + O2
      AI(95,I) = AR(71,I) !same CH3O2 + HOQ  =  CH3OOH + OQ

      AI(96,I) = AR(72,I)*COMBIN !ACCELERATED CH3QO + CH3O2  =  H3CO + H3CQ + O2
      AI(97,I) = AR(72,I)*COMBIN !undoubled !ACCELERATED CH3OQ + CH3O2  =  2H3CO + OQ 

      AI(98,I) = AR(73,I) !same CH3OQ + NO  =  H3CO + NOQ
      AI(99,I) = AR(73,I) !same CH3QO + NO  =  H3CQ + NO2 
      AI(100,I)= AR(73,I) !same CH3O2 + NQ  =  H3CO + NOQ

      AI(101,I)= AR(74,I) !same H3CQ + O2  =  H2CQ + HO2
C AJ LUNAR NEW YEAR'S EVE
C      AI(102,I)= AR(74,I)/2. !same H3CO + OQ  =  H2CO + HOQ 
C AJ 02/06 DEBUG
      AI(102,I)= AR(74,I) !same H3CO + OQ  =  H2CO + HOQ 

      AI(103,I)= AR(75,I) !same H3CQ + O  =  H2CQ + OH
      AI(104,I)= AR(75,I) !same H3CO + Q  =  H2CO + QH

      AI(105,I)= AR(76,I) !same H3CQ + OH  =  H2CQ + H2O   
      AI(106,I)= AR(76,I) !same H3CO + QH  =  H2CO + H2Q

      AI(107,I)= AR(77,I) !same N2Q + O(1D)  =  NQ + NO
      AI(108,I)= AR(77,I) !same N2O + Q(1D)  =  NQ + NO

      AI(109,I)= AR(78,I) !same N2Q + O(1D)  =  N2 + OQ
      AI(110,I)= AR(78,I) !same N2O + Q(1D)  =  N2 + OQ 

      AI(111,I)= AR(79,I)/2. !half  N + OQ  =  NQ + O
      AI(112,I)= AR(79,I)/2. !half  N + OQ  =  NO + Q 

      AI(113,I)= AR(80,I)/2. !half  N + O2Q  =  NQ + O2
      AI(114,I)= AR(80,I)/2. !half  N + O2Q  =  NO + OQ
      AI(115,I)= AR(80,I) !same  N + OQO  =  NO + OQ

      AI(116,I)= AR(81,I) !same N + QH  =  NQ + H
      AI(117,I)= AR(82,I) !same N + NQ  =  N2 + Q    

      AI(118,I)= AR(83,I) !same NQ + O3  =  NOQ + O2
      AI(119,I)= AR(83,I)/2. !half NO + O2Q  =  NOQ + O2
      AI(120,I)= AR(83,I)/2. !half NO + O2Q  =  NO2 + OQ
      AI(121,I)= AR(83,I) !same NO + OQO  =  NO2 + OQ

      AI(122,I)= AR(84,I) !same NQ + O + M  =  NOQ + M 
      AI(123,I)= AR(84,I) !same NO + Q + M  =  NOQ + M

      AI(124,I)= AR(85,I) !same NQ + HO2  =  NOQ + OH
C AJ LUNAR NEW YEAR'S EVE
C      AI(125,I)= AR(85,I) !/2. !half NO + HOQ  =  NOQ + OH
C      AI(126,I)= 0. !AR(85,I)/2. !half NO + HOQ  =  NO2 + QH
C AJ 02/06 DEBUG
      AI(125,I)= AR(85,I)  !half NO + HOQ  =  NOQ + OH
      AI(126,I)= AR(85,I)/2. !half NO + HOQ  =  NO2 + QH

      AI(127,I)= AR(86,I) !same NQ + OH + M  =  HNOQ + M
      AI(128,I)= AR(86,I) !same NO + QH + M  =  HNOQ + M 

      AI(129,I)= AR(87,I)/2. !half NOQ + O  =  NQ + O2    
      AI(130,I)= AR(87,I)/2. !half NOQ + O  =  NO + OQ
      AI(131,I)= AR(87,I) !same NO2 + Q  =  NO + OQ 

      AI(132,I)= AR(88,I) !same NOQ + OH + M  =  HNO2Q + M 
      AI(133,I)= AR(88,I) !same NO2 + QH + M  =  HNO2Q + M

      AI(134,I)= AR(89,I)/2. !half NOQ + H  =  NQ + OH
      AI(135,I)= AR(89,I)/2. !half NOQ + H  =  NO + QH

      AI(136,I)= AR(90,I) !same HNO2Q + OH  =  H2O + NO2Q 
      AI(137,I)= AR(90,I) !same HNO3 + QH  =  H2Q + NO3

      AI(138,I)= AR(91,I) !same HOQ + NO2 + M  =  HOQNO2 + M
      AI(139,I)= AR(91,I) !same HO2 + NOQ + M  =  HO2NOQ + M

      AI(140,I)= AR(92,I) !same HOQNO2 + OH  =  NO2 + H2O + OQ
      AI(141,I)= AR(92,I) !same HO2NO2 + QH  =  NO2 + H2Q + O2    
      AI(142,I)= AR(92,I) !same HO2NOQ + OH  =  NOQ + H2O + O2

      AI(143,I)= AR(93,I) !same HOQNO2 + O  =  NO2 + OH + OQ 
      AI(144,I)= AR(93,I) !same HO2NO2 + Q  =  NO2 + QH + O2
      AI(145,I)= AR(93,I) !same HO2NOQ + O  =  NOQ + OH + O2

      AI(146,I)= AR(94,I) !same HOQNO2 + M  =  HOQ + NO2 + M
      AI(147,I)= AR(94,I) !same HO2NOQ + M  =  HO2 + NOQ + M

      AI(148,I)= AR(96,I) !same CH3OQH + OH  =  CH3OQ + H2O
      AI(149,I)= AR(96,I) !same CH3QOH + OH  =  CH3QO + H2O
      AI(150,I)= AR(96,I) !same CH3OOH + QH  =  CH3O2 + H2Q

      AI(151,I)= AR(97,I) !same CH3OQ + OH  = H3CO + HOQ 
      AI(152,I)= AR(97,I) !same CH3O2 + QH  = H3CO + HOQ 
      AI(153,I)= AR(97,I) !same CH3QO + OH  = H3CQ + HO2

      AI(154,I)= AR(98,I)/2. !half O2Q + NO2  =  OQ + NO3 
      AI(155,I)= AR(98,I)/2. !half O2Q + NO2  =  O2 + NO2Q
      AI(156,I)= AR(98,I) !same OQO + NO2  =  OQ + NO3
      AI(157,I)= AR(98,I) !same O3 + NOQ  =  O2 + NO2Q

      AI(158,I)= AR(99,I)/3.  ! /3 NOQ + NO3  =  NQ + NO2 + O2 
      AI(159,I)= AR(99,I)/3.  ! /3 NOQ + NO3  =  NO + NOQ + O2
      AI(160,I)= AR(99,I)/3.  ! /3 NOQ + NO3  =  NO + NO2 + OQ
      AI(161,I)= AR(99,I)/3.  ! /3 NO2 + NO2Q =  NQ + NO2 + O2
      AI(162,I)= AR(99,I)/3.  ! /3 NO2 + NO2Q =  NO + NOQ + O2
      AI(163,I)= AR(99,I)/3.  ! /3 NO2 + NO2Q =  NO + NO2 + OQ

      AI(164,I)= AR(100,I)    !same Q + NO3  =  OQ + NO2 
      AI(165,I)= AR(100,I)/2. !half O + NO2Q  =  OQ + NO2
      AI(166,I)= AR(100,I)/2. !half O + NO2Q  =  O2 + NOQ 

      AI(167,I)= AR(101,I) !same NQ + NO3  =  NO2 + NOQ  
      AI(168,I)= AR(101,I)/2. !half NO + NO2Q  =  NOQ + NO2  
      AI(169,I)= AR(101,I)/2. !half NO + NO2Q  =  NO2 + NOQ 
C AJ LUNAR NEW YEAR'S EVE
      AI(170,I)= AR(102,I) ! QH + NO3  =  HQO + NO2 
C      AI(171,I)= AR(102,I)/3. !half OH + NO2Q  =  HOQ + NO2
C      AI(172,I)= AR(102,I)*2./3. !half OH + NO2Q  =  HO2 + NOQ
C AJ 02/06 DEBUG
      AI(171,I)= AR(102,I)/2. !half OH + NO2Q  =  HOQ + NO2
      AI(172,I)= AR(102,I)/2. !half OH + NO2Q  =  HO2 + NOQ  

      AI(173,I)= AR(103,I) !same HOQ + NO3  =  HNO3 + OQ
      AI(174,I)= AR(103,I) !same HO2 + NO2Q  =  HNO2Q + O2

      AI(175,I)= AR(104,I) !same NOQ + O + M  =  NO2Q + M
      AI(176,I)= AR(104,I) !same NO2 + Q + M  =  NO2Q + M 

      AI(177,I)= AR(106,I)*0.8!/2. !half NO2Q + NO2 + M  =  N2O4Q + M
      AI(178,I)= AR(106,I)*0.2!/2. !half NO2Q + NO2 + M  =  N2QO4 + M
      AI(179,I)= AR(106,I)*0.8!/2. !half NO3 + NOQ + M  =  N2O4Q + M
      AI(180,I)= AR(106,I)*0.2!/2. !half NO3 + NOQ + M  =  N2QO4 + M

      AI(181,I)= AR(108,I)/2. !half N2O4Q + M  =  NOQ + NO3 + M 
      AI(182,I)= AR(108,I)/2. !half N2O4Q + M  =  NO2 + NO2Q + M
      AI(183,I)= AR(108,I) !same N2QO4 + M  =  NO2 + NO2Q + M 

      AI(184,I)= AR(109,I) !same N2O4Q + H2O  = HNO2Q + HNO3
      AI(185,I)= AR(109,I) !same N2QO4 + H2O  = HNO2Q + HNO3
      AI(186,I)= AR(109,I) !same N2O5 + H2Q  =  HNO2Q + HNO3

      AI(187,I)= AR(113,I) !same  SQ   + O2   =     O    +     SOQ
      AI(188,I)= AR(113,I)/2. !half  SO   + OQ   =     Q    +     SO2
      AI(189,I)= AR(113,I)/2. !half  SO   + OQ   =     O    +     SOQ

      AI(190,I)= AR(114,I) !same SQ   + HO2  =     SOQ  +     OH
C AJ LUNAR NEW YEAR'S EVE
C      AI(191,I)= AR(114,I) !/2. !half  SO   + HOQ  =     SOQ  +     OH
C      AI(192,I)= AR(114,I) !/2. !half  SO   + HOQ  =     SO2  +     QH
C AJ 02/06 DEBUG
      AI(191,I)= AR(114,I)/2.  !half  SO   + HOQ  =     SOQ  +     OH
      AI(192,I)= AR(114,I)/2. !half  SO   + HOQ  =     SO2  +     QH

      AI(193,I)= AR(115,I) !same SQ   + O    =     SOQ 
      AI(194,I)= AR(115,I) !same SO   + Q    =     SOQ

      AI(195,I)= AR(116,I) !same SQ   + OH   =     SOQ  +     H 
      AI(196,I)= AR(116,I) !same SO   + QH   =     SOQ  +     H
      AI(197,I)= AR(117,I) !same SOQ  + OH   =     HSO2Q
      AI(198,I)= AR(117,I) !same SO2  + QH   =     HSO2Q
      AI(199,I)= AR(118,I) !same SOQ  + O    =     SO2Q 
      AI(200,I)= AR(118,I) !same SO2  + Q    =     SO2Q 

      AI(201,I)= AR(119,I) !same SO2Q  + H2O =    H2SO3Q 
      AI(202,I)= AR(119,I) !same SO3  + H2Q  =    H2SO3Q 

      AI(203,I)= AR(120,I) !same HSO2Q + O2   =     HO2  +     SO2Q
C AJ LUNAR NEW YEAR'S EVE
C      AI(204,I)= AR(120,I)/2. !same HSO3 + OQ   =     HOQ  +     SO3
C AJ 02/06 DEBUG
      AI(204,I)= AR(120,I) !same HSO3 + OQ   =     HOQ  +     SO3 
 
      AI(205,I)= AR(121,I) !same HSO2Q + OH   =     H2O  +     SO2Q
      AI(206,I)= AR(121,I) !same HSO3 + QH   =     H2Q  +     SO3
 
      AI(207,I)= AR(122,I) !same HSO2Q + H    =     H2   +     SO2Q

      AI(208,I)= AR(123,I) !same HSO2Q + O    =     OH   +     SO2Q 
      AI(209,I)= AR(123,I) !same HSO3 + Q    =     QH   +     SO3

      AI(210,I)= AR(124,I) !same H2S  + QH   =     H2Q  +     HS

      AI(211,I)= AR(126,I) !same H2S  + Q    =     QH   +     HS
      AI(212,I)= AR(127,I) !same HS   + Q    =     H    +     SQ

      AI(213,I)= AR(128,I)/2. !half HS   + OQ   =     QH   +     SO 
      AI(214,I)= AR(128,I)/2. !half HS   + OQ   =     OH   +     SQ

      AI(215,I)= AR(129,I) !same HS   + HOQ  =     H2S  +     OQ 

      AI(216,I)= AR(131,I) !same HS   + HCQ  =     H2S  +     CQ 
      AI(217,I)= AR(134,I)/2. !half S    + OQ   =     SQ   +     O
      AI(218,I)= AR(134,I)/2. !half S    + OQ   =     SO   +     Q

      AI(219,I)= AR(135,I) !same S    + QH   =     SQ   +     H  
      AI(220,I)= AR(137,I) !same  S    + HOQ  =     HS   +     OQ   
C AJ LUNAR NEW YEAR'S EVE
C      AI(221,I)= AR(138,I) !/2. !half S    + HOQ  =     SQ   +     OH  
C      AI(222,I)= AR(138,I) !/2. !half S    + HQO  =     SO   +     QH
C AJ 02/06 DEBUG
      AI(221,I)= AR(138,I)/2. !half S    + HOQ  =     SQ   +     OH  
      AI(222,I)= AR(138,I)/2. !half S    + HQO  =     SO   +     QH  

      AI(223,I)= AR(139,I) !same HS   + H2CQ =     H2S  +     HCQ 
      AI(224,I)= AR(144,I) !same SO1Q  +    M   =      SO3Q +     M
      AI(225,I)= AR(145,I) !same SO1Q  +    M   =      SOQ  +     M 

      AI(226,I)= AR(148,I) !same SO1Q  +    O2  =      SO2Q  +    O
      AI(227,I)= AR(148,I)/2. !half SO21  +    OQ  =      SO2Q  +    O
      AI(228,I)= AR(148,I)/2. !half SO21  +    OQ  =      SO3  +     Q

      AI(229,I)= AR(149,I)/2. !half SO1Q  +    SO2 =      SO2Q  +    SO
      AI(230,I)= AR(149,I)/2. !half SO1Q  +    SO2 =      SO3  +     SQ
      AI(231,I)= AR(149,I)/2. !half SO21  +    SOQ =      SO2Q  +    SO
      AI(232,I)= AR(149,I)/2. !half SO21  +    SOQ =      SO3  +     SQ

      AI(233,I)= AR(150,I) !same SO3Q  +    M   =      SOQ  +     M

      AI(234,I)= AR(152,I)/2. !half SO3Q  +    SO2 =      SO2Q  +    SO 
      AI(235,I)= AR(152,I)/2. !half SO3Q  +    SO2 =      SO3  +     SQ
      AI(236,I)= AR(152,I)/2. !half SO23  +    SOQ =      SO2Q  +    SO 
      AI(237,I)= AR(152,I)/2. !half SO23  +    SOQ =      SO3  +     SQ 

      AI(238,I)= AR(153,I) !same SQ    +    NO2 =      SOQ  +     NO
      AI(239,I)= AR(153,I)/2. !half SO    +    NOQ =      SOQ  +     NO
      AI(240,I)= AR(153,I)/2. !half SO    +    NOQ =      SO2  +     NQ 

      AI(241,I)= AR(154,I) !same SQ    +    O3  =      SOQ  +     O2
      AI(242,I)= AR(154,I)/2. !half SO    +    O2Q  =     SOQ  +     O2
      AI(243,I)= AR(154,I)/2. !half SO    +    O2Q  =     SO2  +     OQ
      AI(244,I)= AR(154,I) !same SO    +    OQO  =     SO2  +     OQ

      AI(245,I)= AR(155,I) !same SOQ   +    HO2 =      SO2Q  +    OH
C AJ LUNAR NEW YEAR'S EVE
C      AI(246,I)= AR(155,I) !/2. !half SO2   +    HOQ =      SO2Q  +    OH
C      AI(247,I)= AR(155,I) !/2. !half SO2   +    HOQ =      SO3  +     QH
C AJ 02/06 DEBUG
      AI(246,I)= AR(155,I)/2. !half SO2   +    HOQ =      SO2Q  +    OH
      AI(247,I)= AR(155,I)/2. !half SO2   +    HOQ =      SO3  +     QH

      AI(248,I)= AR(156,I)/2. !half HS    +    O2Q  =      HSQ  +     O2
      AI(249,I)= AR(156,I)/2. !half HS    +    O2Q  =      HSO  +     OQ 
      AI(250,I)= AR(156,I) !same HS    +    OQO  =      HSO  +     OQ

      AI(251,I)= AR(157,I)/2. !half HS    +    NOQ =      HSQ  +     NO
      AI(252,I)= AR(157,I)/2. !half HS    +    NOQ =      HSO  +     NQ

      AI(253,I)= AR(158,I)/2. !half S     +    O2Q  =      SQ   +     O2
      AI(254,I)= AR(158,I)/2. !half S     +    O2Q  =      SO   +     OQ
      AI(255,I)= AR(158,I) !same S     +    OQO  =      SO   +     OQ

      AI(256,I)= AR(159,I)*COMBIN !undoubled !ACCELERATED SQ   +    SO  =      SOQ  +     S 

      AI(257,I)= AR(160,I) !same SO2Q  +    SO  =      SOQ  +     SO2
      AI(258,I)= AR(160,I) !same SO3   +    SQ  =      SOQ  +     SO2

      AI(259,I)= AR(161,I)/2. !half S     +    COQ =      SQ   +     CO
      AI(260,I)= AR(161,I)/2. !half S     +    COQ =      SO   +     CQ 

      AI(261,I)= AR(162,I) !same SQ    +    HO2 =      HSQ  +     O2
      AI(262,I)= AR(162,I) !same SO    +    HOQ =      HSO  +     OQ

      AI(263,I)= AR(163,I) !same SQ    +    HCO =      HSQ  +     CO
      AI(264,I)= AR(163,I) !same SO    +    HCQ =      HSO  +     CQ 

      AI(265,I)= AR(164,I) !same H     +    SQ  =      HSQ

      AI(266,I)= AR(166,I) !same HSQ   +    OH  =      H2O  +     SQ
      AI(267,I)= AR(166,I) !same HSO   +    QH  =      H2Q  +     SO

      AI(268,I)= AR(167,I) !same HSQ   +    H   =      HS   +     QH 
      AI(269,I)= AR(168,I) !same HSQ   +    H   =      H2   +     SQ 
      AI(270,I)= AR(169,I) !same HSQ   +    HS  =      H2S  +     SQ
 
      AI(271,I)= AR(170,I) !same HSQ   +    O   =      OH   +     SQ
      AI(272,I)= AR(170,I) !same HSO   +    Q   =      QH   +     SO

      AI(273,I)= AR(171,I) !same HSQ   +    S   =      HS   +     SQ
      AI(274,I)= AR(172,I) !same N2 + Q1D = N2Q
      AI(275,I)= AR(173,I) !same N2Q + H = N2 +OH
      AI(276,I)= AR(174,I) !same N2Q + NO = NO2 + N2
      AI(277,I)= AR(174,I) !same N2O + NQ = NOQ + N2 
C-PL ISOTOPE EXCHANGE
      AI(320,I)= 3.4E-12*((300./T(I))**1.1)!Lessard et al 2003   ~1700 times O3 formation rate at 100hPa
      AI(321,I)= 0.5 * AI(320,I)
      AI(322,I)= 2./3.*6.9E-11*EXP(117./T(I)) !Q1D+CO2=COQ+O
      AI(323,I)= 1./3.*6.9E-11*EXP(117./T(I)) !Q1D+CO2=CO2+Q
      AI(324,I)= 1./3.*6.9E-11*EXP(117./T(I)) !O1D+COQ=CO2+Q
      AI(325,I)= 2./3.*6.9E-11*EXP(117./T(I)) !O1D+COQ=COQ+O

C-JL ADD 28 NEW RXN FROM HQO AND HQONO2 
C AJ 02/06 DEBUG
      AI(326,I)= 0. !AR(6,I)*0.5
      AI(327,I)= 0. !AR(7,I)
      AI(328,I)= 0. !AR(8,I)
      AI(329,I)= 0. !AR(9,I)
      AI(330,I)= 0. !AR(11,I)
      AI(331,I)= 0. !AR(12,I)
      AI(332,I)= 0. !AR(13,I)
      AI(333,I)= 0. !AR(14,I)
      AI(334,I)= 0. !AR(15,I) !*0.5
      AI(335,I)= 0. !AR(15,I) !*0.5
      AI(336,I)= 0. !AR(16,I)*0.5
      AI(337,I)= 0. !AR(44,I)*0.5
      AI(338,I)= 0. !AR(49,I)*0.5
c      AI(339,I)= AR(50,I)
      AI(340,I)= 0. !AR(70,I)*0.5
      AI(341,I)= 0. !AR(71,I)
      AI(342,I)= 0. !AR(74,I)*0.5
      AI(343,I)= 0. !AR(85,I)
      AI(344,I)= 0. !AR(91,I)
      AI(345,I)= 0. !AR(92,I)
      AI(346,I)= 0. !AR(93,I)
      AI(347,I)= 0. !AR(94,I)
c      AI(348,I)= AR(95,I)
      AI(349,I)= 0. !AR(103,I)
      AI(350,I)= 0. !AR(120,I)*0.5
      AI(351,I)= 0. !AR(129,I)
      AI(352,I)= 0. !AR(137,I)
      AI(353,I)= 0. !AR(162,I)
   4  CONTINUE
!      print *, 'O3 peng', Z(19),AI(34,19),T(19)
!      print *, 'EXC peng',Z(19),AI(320,19)
      DO 10 I=1,NZ 

        RAINGCI(LH2CQ,I) = RAINGC(LH2CO,I)
        RAINGCI(LQ,I) = RAINGC(LO,I)
        RAINGCI(LH2Q,I) = RAINGC(LH2O,I)
        RAINGCI(LQH,I) = RAINGC(LOH,I)
        RAINGCI(LHOQ,I) = RAINGC(LHO2,I)
        RAINGCI(LH2OQ,I) = RAINGC(LH2O2,I)
        RAINGCI(LO2Q,I) = RAINGC(LO3,I)
        RAINGCI(LOQO,I) = RAINGC(LO3,I)
        RAINGCI(LCQ,I) = RAINGC(LCO,I)
        RAINGCI(LCH3OQH,I) = RAINGC(LCH3OOH,I)
        RAINGCI(LCH3QOH,I) = RAINGC(LCH3OOH,I)
        RAINGCI(LCH3OQ,I) = RAINGC(LCH3O2,I)
        RAINGCI(LCH3QO,I) = RAINGC(LCH3O2,I)
        RAINGCI(LN2Q,I) = RAINGC(LN2O,I)
        RAINGCI(LNQ,I) = RAINGC(LNO,I)
        RAINGCI(LNOQ,I) = RAINGC(LNO2,I)
        RAINGCI(LHNOQ,I) = RAINGC(LHNO2,I)
        RAINGCI(LHNO2Q,I) = RAINGC(LHNO3,I)
        RAINGCI(LHOQNO2,I) = RAINGC(LHO2NO2,I)
        RAINGCI(LHO2NOQ,I) = RAINGC(LHO2NO2,I)
        RAINGCI(LNO2Q,I) = RAINGC(LNO3,I)
        RAINGCI(LN2O4Q,I) = RAINGC(LN2O5,I)
        RAINGCI(LN2QO4,I) = RAINGC(LN2O5,I)
        RAINGCI(LOQ,I) = RAINGC(LO2,I)
        RAINGCI(LSQ,I) = RAINGC(LSO,I)
        RAINGCI(LSOQ,I) = RAINGC(LSO2,I)
        RAINGCI(LH2SO3Q,I) = RAINGC(LH2SO4,I)
        RAINGCI(LHSQ,I) = RAINGC(LHSO,I)
        RAINGCI(LCOQ,I) = RAINGC(LCO2,I)
! AJ 01/23/2020 ADD TWO SPECIES
C        RAINGCI(LHQO,I) = RAINGC(LHO2,I)
C        RAINGCI(LHQONO2,I) = RAINGC(LHO2NO2,I)

  10  CONTINUE
      DO 11 I=1,NZ
C 1
        DUI(LH2CQ,I) = DU(LH2CO,I)
        DLI(LH2CQ,I) = DL(LH2CO,I)
        DDI(LH2CQ,I) = DD(LH2CO,I)
        DHUI(LH2CQ,I) = DHU(LH2CO,I)
        DHLI(LH2CQ,I) = DHL(LH2CO,I)
        HIZ(LH2CQ,I) = HI(LH2CO,I)
        DKI(LH2CQ,I) = DK(LH2CO,I)
C 2
        DUI(LQ,I) = DU(LO,I)
        DLI(LQ,I) = DL(LO,I)
        DDI(LQ,I) = DD(LO,I)
        DHUI(LQ,I) = DHU(LO,I)
        DHLI(LQ,I) = DHL(LO,I)
        HIZ(LQ,I) = HI(LO,I)
        DKI(LQ,I) = DK(LO,I)
C 3
        DUI(LH2Q,I) = DU(LH2O,I)
        DLI(LH2Q,I) = DL(LH2O,I)
        DDI(LH2Q,I) = DD(LH2O,I)
        DHUI(LH2Q,I) = DHU(LH2O,I)
        DHLI(LH2Q,I) = DHL(LH2O,I)
        HIZ(LH2Q,I) = HI(LH2O,I)
        DKI(LH2Q,I) = DK(LH2O,I)
C 4
        DUI(LQH,I) = DU(LOH,I)
        DLI(LQH,I) = DL(LOH,I)
        DDI(LQH,I) = DD(LOH,I)
        DHUI(LQH,I) = DHU(LOH,I)
        DHLI(LQH,I) = DHL(LOH,I)
        HIZ(LQH,I) = HI(LOH,I)
        DKI(LQH,I) = DK(LOH,I)
C 5
        DUI(LHOQ,I) = DU(LHO2,I)
        DLI(LHOQ,I) = DL(LHO2,I)
        DDI(LHOQ,I) = DD(LHO2,I)
        DHUI(LHOQ,I) = DHU(LHO2,I)
        DHLI(LHOQ,I) = DHL(LHO2,I)
        HIZ(LHOQ,I) = HI(LHO2,I)
        DKI(LHOQ,I) = DK(LHO2,I)
C 6
        DUI(LH2OQ,I) = DU(LH2O2,I)
        DLI(LH2OQ,I) = DL(LH2O2,I)
        DDI(LH2OQ,I) = DD(LH2O2,I)
        DHUI(LH2OQ,I) = DHU(LH2O2,I)
        DHLI(LH2OQ,I) = DHL(LH2O2,I)
        HIZ(LH2OQ,I) = HI(LH2O2,I)
        DKI(LH2OQ,I) = DK(LH2O2,I)
C 7
        DUI(LO2Q,I) = DU(LO3,I)
        DLI(LO2Q,I) = DL(LO3,I)
        DDI(LO2Q,I) = DD(LO3,I)
        DHUI(LO2Q,I) = DHU(LO3,I)
        DHLI(LO2Q,I) = DHL(LO3,I)
        HIZ(LO2Q,I) = HI(LO3,I)
        DKI(LO2Q,I) = DK(LO3,I)
C 8
        DUI(LOQO,I) = DU(LO3,I)
        DLI(LOQO,I) = DL(LO3,I)
        DDI(LOQO,I) = DD(LO3,I)
        DHUI(LOQO,I) = DHU(LO3,I)
        DHLI(LOQO,I) = DHL(LO3,I)
        HIZ(LOQO,I) = HI(LO3,I)
        DKI(LOQO,I) = DK(LO3,I)
C 9
        DUI(LCQ,I) = DU(LCO,I)
        DLI(LCQ,I) = DL(LCO,I)
        DDI(LCQ,I) = DD(LCO,I)
        DHUI(LCQ,I) = DHU(LCO,I)
        DHLI(LCQ,I) = DHL(LCO,I)
        HIZ(LCQ,I) = HI(LCO,I)
        DKI(LCQ,I) = DK(LCO,I)
C 10
        DUI(LCH3OQH,I) = DU(LCH3OOH,I)
        DLI(LCH3OQH,I) = DL(LCH3OOH,I)
        DDI(LCH3OQH,I) = DD(LCH3OOH,I)
        DHUI(LCH3OQH,I) = DHU(LCH3OOH,I)
        DHLI(LCH3OQH,I) = DHL(LCH3OOH,I)
        HIZ(LCH3OQH,I) = HI(LCH3OOH,I)
        DKI(LCH3OQH,I) = DK(LCH3OOH,I)
C 11
        DUI(LCH3QOH,I) = DU(LCH3OOH,I)
        DLI(LCH3QOH,I) = DL(LCH3OOH,I)
        DDI(LCH3QOH,I) = DD(LCH3OOH,I)
        DHUI(LCH3QOH,I) = DHU(LCH3OOH,I)
        DHLI(LCH3QOH,I) = DHL(LCH3OOH,I)
        HIZ(LCH3QOH,I) = HI(LCH3OOH,I)
        DKI(LCH3QOH,I) = DK(LCH3OOH,I)
C 12
        DUI(LCH3OQ,I) = DU(LCH3O2,I)
        DLI(LCH3OQ,I) = DL(LCH3O2,I)
        DDI(LCH3OQ,I) = DD(LCH3O2,I)
        DHUI(LCH3OQ,I) = DHU(LCH3O2,I)
        DHLI(LCH3OQ,I) = DHL(LCH3O2,I)
        HIZ(LCH3OQ,I) = HI(LCH3O2,I)
        DKI(LCH3OQ,I) = DK(LCH3O2,I)
C 13
        DUI(LCH3QO,I) = DU(LCH3O2,I)
        DLI(LCH3QO,I) = DL(LCH3O2,I)
        DDI(LCH3QO,I) = DD(LCH3O2,I)
        DHUI(LCH3QO,I) = DHU(LCH3O2,I)
        DHLI(LCH3QO,I) = DHL(LCH3O2,I)
        HIZ(LCH3QO,I) = HI(LCH3O2,I)
        DKI(LCH3QO,I) = DK(LCH3O2,I)
C 14
        DUI(LN2Q,I) = DU(LN2O,I)
        DLI(LN2Q,I) = DL(LN2O,I)
        DDI(LN2Q,I) = DD(LN2O,I)
        DHUI(LN2Q,I) = DHU(LN2O,I)
        DHLI(LN2Q,I) = DHL(LN2O,I)
        HIZ(LN2Q,I) = HI(LN2O,I)
        DKI(LN2Q,I) = DK(LN2O,I)
C 15
        DUI(LNQ,I) = DU(LNO,I)
        DLI(LNQ,I) = DL(LNO,I)
        DDI(LNQ,I) = DD(LNO,I)
        DHUI(LNQ,I) = DHU(LNO,I)
        DHLI(LNQ,I) = DHL(LNO,I)
        HIZ(LNQ,I) = HI(LNO,I)
        DKI(LNQ,I) = DK(LNO,I)
C 16
        DUI(LNOQ,I) = DU(LNO2,I)
        DLI(LNOQ,I) = DL(LNO2,I)
        DDI(LNOQ,I) = DD(LNO2,I)
        DHUI(LNOQ,I) = DHU(LNO2,I)
        DHLI(LNOQ,I) = DHL(LNO2,I)
        HIZ(LNOQ,I) = HI(LNO2,I)
        DKI(LNOQ,I) = DK(LNO2,I)
C 17
        DUI(LHNOQ,I) = DU(LHNO2,I)
        DLI(LHNOQ,I) = DL(LHNO2,I)
        DDI(LHNOQ,I) = DD(LHNO2,I)
        DHUI(LHNOQ,I) = DHU(LHNO2,I)
        DHLI(LHNOQ,I) = DHL(LHNO2,I)
        HIZ(LHNOQ,I) = HI(LHNO2,I)
        DKI(LHNOQ,I) = DK(LHNO2,I)
C 18
        DUI(LHNO2Q,I) = DU(LHNO3,I)
        DLI(LHNO2Q,I) = DL(LHNO3,I)
        DDI(LHNO2Q,I) = DD(LHNO3,I)
        DHUI(LHNO2Q,I) = DHU(LHNO3,I)
        DHLI(LHNO2Q,I) = DHL(LHNO3,I)
        HIZ(LHNO2Q,I) = HI(LHNO3,I)
        DKI(LHNO2Q,I) = DK(LHNO3,I)
C 19
        DUI(LHOQNO2,I) = DU(LHO2NO2,I)
        DLI(LHOQNO2,I) = DL(LHO2NO2,I)
        DDI(LHOQNO2,I) = DD(LHO2NO2,I)
        DHUI(LHOQNO2,I) = DHU(LHO2NO2,I)
        DHLI(LHOQNO2,I) = DHL(LHO2NO2,I)
        HIZ(LHOQNO2,I) = HI(LHO2NO2,I)
        DKI(LHOQNO2,I) = DK(LHO2NO2,I)
C 20
        DUI(LHO2NOQ,I) = DU(LHO2NO2,I)
        DLI(LHO2NOQ,I) = DL(LHO2NO2,I)
        DDI(LHO2NOQ,I) = DD(LHO2NO2,I)
        DHUI(LHO2NOQ,I) = DHU(LHO2NO2,I)
        DHLI(LHO2NOQ,I) = DHL(LHO2NO2,I)
        HIZ(LHO2NOQ,I) = HI(LHO2NO2,I)
        DKI(LHO2NOQ,I) = DK(LHO2NO2,I)
C 21
        DUI(LNO2Q,I) = DU(LNO3,I)
        DLI(LNO2Q,I) = DL(LNO3,I)
        DDI(LNO2Q,I) = DD(LNO3,I)
        DHUI(LNO2Q,I) = DHU(LNO3,I)
        DHLI(LNO2Q,I) = DHL(LNO3,I)
        HIZ(LNO2Q,I) = HI(LNO3,I)
        DKI(LNO2Q,I) = DK(LNO3,I)
C 22
        DUI(LN2O4Q,I) = DU(LN2O5,I)
        DLI(LN2O4Q,I) = DL(LN2O5,I)
        DDI(LN2O4Q,I) = DD(LN2O5,I)
        DHUI(LN2O4Q,I) = DHU(LN2O5,I)
        DHLI(LN2O4Q,I) = DHL(LN2O5,I)
        HIZ(LN2O4Q,I) = HI(LN2O5,I)
        DKI(LN2O4Q,I) = DK(LN2O5,I)
C 23
        DUI(LN2QO4,I) = DU(LN2O5,I)
        DLI(LN2QO4,I) = DL(LN2O5,I)
        DDI(LN2QO4,I) = DD(LN2O5,I)
        DHUI(LN2QO4,I) = DHU(LN2O5,I)
        DHLI(LN2QO4,I) = DHL(LN2O5,I)
        HIZ(LN2QO4,I) = HI(LN2O5,I)
        DKI(LN2QO4,I) = DK(LN2O5,I)
C 24
        DUI(LOQ,I) = DU(LO2,I)
        DLI(LOQ,I) = DL(LO2,I)
        DDI(LOQ,I) = DD(LO2,I)
        DHUI(LOQ,I) = DHU(LO2,I)
        DHLI(LOQ,I) = DHL(LO2,I)
        HIZ(LOQ,I) = HI(LO2,I)
        DKI(LOQ,I) = DK(LO2,I)
C 25
        DUI(LSQ,I) = DU(LSO,I)
        DLI(LSQ,I) = DL(LSO,I)
        DDI(LSQ,I) = DD(LSO,I)
        DHUI(LSQ,I) = DHU(LSO,I)
        DHLI(LSQ,I) = DHL(LSO,I)
        HIZ(LSQ,I) = HI(LSO,I)
        DKI(LSQ,I) = DK(LSO,I)
C 26
        DUI(LSOQ,I) = DU(LSO2,I)
        DLI(LSOQ,I) = DL(LSO2,I)
        DDI(LSOQ,I) = DD(LSO2,I)
        DHUI(LSOQ,I) = DHU(LSO2,I)
        DHLI(LSOQ,I) = DHL(LSO2,I)
        HIZ(LSOQ,I) = HI(LSO2,I)
        DKI(LSOQ,I) = DK(LSO2,I)
C 27
        DUI(LH2SO3Q,I) = DU(LH2SO4,I)
        DLI(LH2SO3Q,I) = DL(LH2SO4,I)
        DDI(LH2SO3Q,I) = DD(LH2SO4,I)
        DHUI(LH2SO3Q,I) = DHU(LH2SO4,I)
        DHLI(LH2SO3Q,I) = DHL(LH2SO4,I)
        HIZ(LH2SO3Q,I) = HI(LH2SO4,I)
        DKI(LH2SO3Q,I) = DK(LH2SO4,I)
C 28
        DUI(LHSQ,I) = DU(LHSO,I)
        DLI(LHSQ,I) = DL(LHSO,I)
        DDI(LHSQ,I) = DD(LHSO,I)
        DHUI(LHSQ,I) = DHU(LHSO,I)
        DHLI(LHSQ,I) = DHL(LHSO,I)
        HIZ(LHSQ,I) = HI(LHSO,I)
        DKI(LHSQ,I) = DK(LHSO,I)
C 29
        DUI(LCOQ,I) = DU(LCO2,I)
        DLI(LCOQ,I) = DL(LCO2,I)
        DDI(LCOQ,I) = DD(LCO2,I)
        DHUI(LCOQ,I) = DHU(LCO2,I)
        DHLI(LCOQ,I) = DHL(LCO2,I)
        HIZ(LCOQ,I) = HI(LCO2,I)
        DKI(LCOQ,I) = DK(LCO2,I)
! AJ 01/23/2020 ADD TWP LONG SPECIES
C 30
C        DUI(LHQO,I) = DU(LHO2,I)
C        DLI(LHQO,I) = DL(LHO2,I)
C        DDI(LHQO,I) = DD(LHO2,I)
C        DHUI(LHQO,I) = DHU(LHO2,I)
C        DHLI(LHQO,I) = DHL(LHO2,I)
C        HIZ(LHQO,I) = HI(LHO2,I)
C        DKI(LHQO,I) = DK(LHO2,I)
C 31
C        DUI(LHQONO2,I) = DU(LHO2NO2,I)
C        DLI(LHQONO2,I) = DL(LHO2NO2,I)
C        DDI(LHQONO2,I) = DD(LHO2NO2,I)
C        DHUI(LHQONO2,I) = DHU(LHO2NO2,I)
C        DHLI(LHQONO2,I) = DHL(LHO2NO2,I)
C        HIZ(LHQONO2,I) = HI(LHO2NO2,I)
C        DKI(LHQONO2,I) = DK(LHO2NO2,I)

  11  CONTINUE
C
      RETURN
      END
