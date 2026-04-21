
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
      DIMENSION TTAB(NT)
      DIMENSION SRO2(KNEW,NL,NT), SRO2L(KNEW,NL,NT)

C     Temperature and pressure grids
      DATA TTAB/150., 200., 250., 300./
c      DATA PTAB/1., 0.1, 0.01, 0.001, 0.0001/

C     READ THE DATA FILE 'k12Gauss.dat' 
 200  FORMAT(10E10.3)
      DO I=1,NT
c      PRINT *
c      PRINT *,'I = ',I, TTAB(I)
c      DO J=1,NP
C      PRINT *
c      PRINT *,'J =',J, PTAB(J)
      READ(10,100)
 100  FORMAT(///////) !SKIP 8 LINES
      DO L=1,NL       !LOOP OVER WAVELENGTH
      READ(10,101)(SRO2(K,L,I),K=1,6)
 101  FORMAT(6X,6(1XE14.7))
c      if (I.GT.1. OR. J.GT.1.) GO TO 1
C      print 101,(SRO2(K,L,J,I),K=1,6)
C   1  continue
   
      READ(10,101)(SRO2(K,L,I),K=7,12)
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
      DO I=1,NT
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
      DO I=1,NT   ! Temperature loop
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
c      IF (N .EQ. 1) THEN 
c      PRINT *,'N = ', N    
c      PRINT *,'Temp = ', T(N), 'TTAB = ',TTAB(IS)
c      PRINT *,'IS = ', IS
c      PRINT *,'FT = ', FT
c      PRINT *,'P = ',P(N)
c      PRINT *,'PL = ', PL(N),'PTABL = ',PTABL(JS)
c      PRINT *,'JS = ', JS
c      PRINT *,'FP = ', FP
c      END IF
C    2 continue

      DO L=1,NL 
c      IF (N .GE. 5) GO TO 3
c      PRINT *, 'L = ', L
c    3 continue

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
