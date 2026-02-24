
      SUBROUTINE SRBANDS(T)
      
      INCLUDE '../INCLUDECHEM/parCORRK.inc'
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'

      INCLUDE '../INCLUDECHEM/comPRESS1.inc'
      INCLUDE '../INCLUDECHEM/comCORRK.inc'

C     THIS SUBROUTINE USES CORRELATED-K TABLES FOR PARAMETERZING O2
C     PHOTOLYSIS RATES AT THE SR BANDS (175-205 nm)
      
      DIMENSION P(NZ), PL(NZ), T(NZ)
C     PL is the log10 of pressure
      DIMENSION TTAB(NT), PTAB(NP), PTABL(NP)
      DIMENSION SRO2(KNEW,NL,NP,NT), SRO2L(KNEW,NL,NP,NT)

C     Temperature and pressure grids
      DATA TTAB/150., 200., 250., 300./
      DATA PTAB/1., 0.1, 0.01, 0.001, 0.0001/

C     READ THE DATA FILE 'k12Gauss.dat' 
 200  FORMAT(10E10.3)
      DO I=1,NT
c      PRINT *
c      PRINT *,'I = ',I, TTAB(I)
      DO J=1,NP
C      PRINT *
c      PRINT *,'J =',J, PTAB(J)
      READ(10,100)
 100  FORMAT(///////) !SKIP 8 LINES
      DO L=1,NL       !LOOP OVER WAVELENGTH
      READ(10,101)(SRO2(K,L,J,I),K=1,6)
 101  FORMAT(6X,6(1XE14.7))
c      if (I.GT.1. OR. J.GT.1.) GO TO 1
C      print 101,(SRO2(K,L,J,I),K=1,6)
C   1  continue
   
      READ(10,101)(SRO2(K,L,J,I),K=7,12)
 102  FORMAT(4X,I2,6(1X,1PE14.7))
c      print 102,L,(SRO2(K,L,J,I),K=7,12)
      READ(10,*)
   
      END DO !END L LOOP
      END DO !END J LOOP  
      END DO !END I LOOP
C
C     CONVERT THE K-COEFFICIENT AND THE PRESSURE GRID TO LOG10 UNITS
      DO J=1,NP      
      PTABL(J)=LOG10(PTAB(J))      
      END DO 
c      print 200,(T(N),N=1,NZ) for debug
C
      DO I=1,NT
      DO J=1,NP
      DO L=1,NL 
      DO K=1,KNEW 
      SRO2L(K,L,J,I)=LOG10(SRO2(K,L,J,I))
      END DO 
      END DO
      END DO
      END DO
C      PRINT *,'SRO2L=',SRO2L(1,1,1,NZ)
C
C     PRESS is in CGS units, and calculated in subroutine DENSITY
C     Should be converted to bar in log10
      DO N=1,NZ
      P(N)=PRESS(N)/1.E6
      PL(N)=LOG10(P(N))
      END DO 
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
      DO J=1,NP   ! Pressure loop
            JS=J
            IF(PTABL(J).LT.PL(N))EXIT 
      END DO 
C     PL(N) is between PTABL(JS) and PTABL(JS-1)
      FP=PL(N) - PTABL(JS)
C      PRINT *,'N =',N,'FT=',FT
c      PRINT *,'N =',N,'FP=',FP
C     using JS here instead of JS-1 to keep FP positive
C     pressure grid is in decreasing order
c-AJ 0129/2023 DEBUG
      IF (N .EQ. 91) THEN 
      PRINT *,'N = ', N    
      PRINT *,'Temp = ', T(N), 'TTAB = ',TTAB(IS)
      PRINT *,'IS = ', IS
      PRINT *,'FT = ', FT
      PRINT *,'P = ',P(N)
      PRINT *,'PL = ', PL(N),'PTABL = ',PTABL(JS)
      PRINT *,'JS = ', JS
      PRINT *,'FP = ', FP
      END IF
C    2 continue

      DO L=1,NL 
c      IF (N .GE. 5) GO TO 3
c      PRINT *, 'L = ', L
c    3 continue

      DO K=1,KNEW 
            IF(FP.GE.0 .AND. FT.LE.1.) THEN
            BETAL = SRO2L(K,L,JS-1,IS-1)*(1.-FT)*FP
     2              + SRO2L(K,L,JS-1,IS)*FT*FP
     3              + SRO2L(K,L,JS,IS-1)*(1.-FT)*(1.-FP)
     4              + SRO2L(K,L,JS,IS)*FT*(1.-FP)
            BETA(K,L,N) = 10.**BETAL
            END IF
            IF(FP.LT.0 .AND. FT.LE.1) THEN
            BETAL = SRO2L(K,L,JS,IS-1)*(1.-FT)
     2              + SRO2L(K,L,JS,IS)*FT
            BETA(K,L,N) = 10.**BETAL
            END IF
            IF(FP.GE.0 .AND. FT.GT.1) THEN
            BETAL = SRO2L(K,L,JS-1,IS)*FP
     2              + SRO2L(K,L,JS,IS)*(1.-FP)
            BETA(K,L,N) = 10.**BETAL
            END IF
c      IF (N .GT. 1 .OR. L.GT.1 .OR. K .GT. 1) GO TO 5
            IF(L.GT.1) GO TO 5
            IF (N .EQ. 91) THEN
            PRINT *
            PRINT *,'K = ',K
            PRINT *, 'JS-1,IS-1 ',SRO2(K,L,JS-1,IS-1)
            PRINT *, 'JS-1,IS ',SRO2(K,L,JS-1,IS)
            PRINT *, 'JS,IS-1 ',SRO2(K,L,JS,IS-1)
            PRINT *, 'JS,IS ',SRO2(K,L,JS,IS)
            PRINT *,'BETA = ', BETA(K,L,N)
            END IF
    5 continue
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
