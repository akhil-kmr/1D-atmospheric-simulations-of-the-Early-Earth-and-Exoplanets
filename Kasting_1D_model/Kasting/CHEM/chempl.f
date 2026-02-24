      SUBROUTINE CHEMPL(D,XP,XL,K)

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc' 
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'    
      
      DIMENSION XP(NZ),XL(NZ),D(NSP2,NZ)

      INCLUDE '../INCLUDECHEM/comRBLOK.inc'
C
      DO 1 I=1,NZ
      XP(I) = 0.
   1  XL(I) = 0.
C  JK 9/1/2022 The commented out lines are a dependency checker
C      print *,'In CHEMPL K =',K
C      print *, 'Start XL loop'

C   LOSS FREQUENCY XL
      NL = NUML(K)
      DO 2 L=1,NL
      J = ILOSS(1,K,L)
      M = ILOSS(2,K,L)
      DO 2 I=1,NZ
c      if(D(M,I).lt.1.e-60.and.i.eq.50) print 500, L,J,D(M,I)
c 500  format('L =',i3,' J =',i3, ' D(M,I) =',1pe10.3)
   2  XL(I) = XL(I) + AR(J,I)*D(M,I)

C      print *, 'Start XP loop'

C   PRODUCTION RATE XP
      NP = NUMP(K)
      DO 3 L=1,NP
      J = IPROD(K,L)
      M = JCHEM(1,J)
      N = JCHEM(2,J)
      DO 3 I=1,NZ
c      if(D(M,I).lt.1.e-60.and.i.eq.50) print 500, L,J,D(M,I)
c      if(D(N,I).lt.1.e-60.and.i.eq.50) print 501, L,J,D(N,I)
c 501  format('L =',i3,' J =',i3, ' D(N,I) =',1pe10.3)
   3  XP(I) = XP(I) + AR(J,I)*D(M,I)*D(N,I)
C 9/2/2022 For dependency checker, put 'stop' in the dochem.f or in the main code after calling the dochem for unperturbed run.
c To pick a layer for dependency checker, it should be neither at the bottom nor the top. Because some species are zero there. 
C
      RETURN
      END
