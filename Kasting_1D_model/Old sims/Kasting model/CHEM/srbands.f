
      SUBROUTINE SRBANDS(T)
      
      INCLUDE '../INCLUDECHEM/parCORRK.inc'
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'

      INCLUDE '../INCLUDECHEM/comPRESS1.inc'
      INCLUDE '../INCLUDECHEM/comCORRK.inc'

C     THIS SUBROUTINE USES CORRELATED-K TABLES FOR PARAMETERZING O2
C     PHOTOLYSIS RATES AT THE SR BANDS (175-205 nm)
      
      DIMENSION T(NZ)
C     PL is the log10 of pressure
      DIMENSION TTAB(NumT),RO2_INITIAL(NL,NumT)
      DIMENSION SRO2(KNEW,NL,NumT), SRO2L(KNEW,NL,NumT)
      DIMENSION RO2_1(NL),RO2_2(NL),RO2_3(NL),RO2_4(NL)
C     Temperature and pressure grids
      DATA TTAB/150., 200., 250., 300./
c      DATA PTAB/1., 0.1, 0.01, 0.001, 0.0001/
C     09/08/2023 AJ ADD BRANCHING RATIO FOR SRC AT 150K, 200K, 250K, AND 300K (From Orlando's data)
      DATA RO2_1/7.41E-2,0.49E-2,0.3E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2  0.,0.,0.,0./
      DATA RO2_2/13.01E-2,0.77E-2,0.29E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2  0.,0.,0.,0./
      DATA RO2_3/18.9E-2,1.54E-2,0.39E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2  0.,0.,0.,0./
      DATA RO2_4/24.69E-2,2.93E-2,0.72E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2  0.,0.,0.,0./
C     save RO2_1, RO2_2, RO2_3, RO2_4 into a matrix RO2_INITIAL(NL,NumT)
      DO I = 1,NL
            RO2_INITIAL(I,1)=RO2_1(I)
            RO2_INITIAL(I,2)=RO2_2(I)
            RO2_INITIAL(I,3)=RO2_3(I)
            RO2_INITIAL(I,4)=RO2_4(I)
      END DO    
C     READ THE DATA FILE 'k4_0927.dat' 
 200  FORMAT(10E10.3)
      DO I=1,NumT
c      PRINT *
c      PRINT *,'I = ',I, TTAB(I)
c      DO J=1,NP
C      PRINT *
c      PRINT *,'J =',J, PTAB(J)
      READ(10,100)
 100  FORMAT(///////) !SKIP 8 LINES
      DO L=1,NL       !LOOP OVER WAVELENGTH
      READ(10,101)(SRO2(K,L,I),K=1,4)
 101  FORMAT(6X,6(1X,E14.7))
c      if (I.GT.1. OR. J.GT.1.) GO TO 1
C      print 101,(SRO2(K,L,J,I),K=1,6)
C   1  continue
   
      READ(10,*)
 102  FORMAT(4X,I2,6(1X,1PE14.7))
c      print 102,L,(SRO2(K,L,J,I),K=7,12)
      READ(10,*)
   
      END DO !END L LOOP
c      END DO !END J LOOP  
      END DO !END I LOOP
C
C     CONVERT THE K-COEFFICIENT AND THE PRESSURE GRID TO LOG10 UNITS
c      DO J=1,NP      
c      PTABL(J)=LOG10(PTAB(J))      
c      END DO 
c      print 200,(T(N),N=1,NZ) for debug
C
      DO I=1,NumT
c      DO J=1,NP
      DO L=1,NL 
      DO K=1,KNEW 
      SRO2L(K,L,I)=LOG10(SRO2(K,L,I))
      END DO 
c      END DO
      END DO
      END DO
C      PRINT *,'SRO2L=',SRO2L(1,1,1,NZ)
C
C     PRESS is in CGS units, and calculated in subroutine DENSITY
C     Should be converted to bar in log10
c      DO N=1,NZ
c      P(N)=PRESS(N)/1.E6
c      PL(N)=LOG10(P(N))
c      END DO 
C    
C     START THE INTEPOLATION FROM THE BOTTOM OF THE ATMOSPHERE
      DO N=1,NZ   ! Altitude loop
c      print *
c      print *,'N = ', N
      DO I=1,NumT   ! Temperature loop
            IS=I 
            IF(TTAB(I).GT.T(N))EXIT 
      END DO 
C     T(N) is between TTAB(IS) and TTAB(IS-1)
      FT=(T(N) - TTAB(IS-1))/50. !Assume a 50 K temperature grid
C     FT is the fractional distance between these grid pts
c      DO J=1,NP   ! Pressure loop
c            JS=J
c            IF(PTABL(J).LT.PL(N))EXIT 
c      END DO 
C     PL(N) is between PTABL(JS) and PTABL(JS-1)
c      FP=PL(N) - PTABL(JS)
c      PRINT *,'N =',N,'FT=',FT
c      PRINT *,'N =',N,'FP=',FP
C     using JS here instead of JS-1 to keep FP positive
C     pressure grid is in decreasing order
c-AJ 0129/2023 DEBUG
      IF (N .EQ. 91) THEN 
      PRINT *,'N = ', N    
      PRINT *,'Temp = ', T(N), 'TTAB = ',TTAB(IS)
c      PRINT *,'IS = ', IS
c      PRINT *,'FT = ', FT
c      PRINT *,'P = ',P(N)
c      PRINT *,'PL = ', PL(N),'PTABL = ',PTABL(JS)
c      PRINT *,'JS = ', JS
c      PRINT *,'FP = ', FP
      END IF
C    2 continue

      DO L=1,NL 
c      IF (N .GE. 5) GO TO 3
c      PRINT *, 'L = ', L
c    3 continue

C 09/08/2023 AJ ADD BRANCHING RATIO FOR SRC AT 150K, 200K, 250K, AND 300K (From Orlando's data)
      IF(FT.LE.1) THEN
      RO2(L,N) = RO2_INITIAL(L,IS-1)*(1.-FT)
     2           + RO2_INITIAL(L,IS)*FT
      END IF

      DO K=1,KNEW 
            IF(FT.LE.1) THEN
            BETAL = SRO2L(K,L,IS-1)*(1.-FT)
     2              + SRO2L(K,L,IS)*FT
            BETA(K,L,N) = 10.**BETAL
            END IF

c      IF (N .GT. 1 .OR. L.GT.1 .OR. K .GT. 1) GO TO 5
c            IF (L.EQ.17 .AND. N .EQ. 1) THEN
c            PRINT *
c            PRINT *,'K = ',K
c            PRINT *, 'IS-1 ',SRO2(K,L,IS-1)
c            PRINT *, 'IS ',SRO2(K,L,IS)
c            PRINT *,'BETA = ', BETA(K,L,N)
c            END IF
C      IF (N .GE. 5) GO TO 4
C      PRINT *,'BETA = ', BETA(K,L,N)
C    4 continue

      END DO  !END Gauss loop (K)
      END DO  !END Wavelength loop(L)
      END DO  !END Altitude loop(N)
  104 FORMAT(1X,1P12E9.2)
      DO I = 1,NZ
      DO J = 18,35
      RO2(J,I)=0
      END DO
      END DO
      DO I = 1,NZ
      WRITE(100,*) T(I)
      WRITE(100,104) (RO2(J,I),J=1,35)
      END DO
c      PRINT *
c      PRINT *,'BETA(1,1,N)'
c      PRINT *, 'K = 1'
c      PRINT *,'L = 1'
c      PRINT 300, (N,BETA(1,1,N),N=1,NZ,30)
c      PRINT *,'L = 17'
c      PRINT 300, (N,BETA(1,17,N),N=1,NZ,30)
c      PRINT *, 'K = 12'
c      PRINT *,'L = 1'
c      PRINT 300, (N,BETA(12,1,N),N=1,NZ,30)
c      PRINT *,'L = 17'
c      PRINT 300, (N,BETA(12,17,N),N=1,NZ,30)   
  300 FORMAT(I2,1PE12.5)

      RETURN
      END
