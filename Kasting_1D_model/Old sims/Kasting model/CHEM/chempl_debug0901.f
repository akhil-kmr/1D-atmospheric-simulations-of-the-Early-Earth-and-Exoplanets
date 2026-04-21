
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
C
C   LOSS FREQUENCY XL
      NL = NUML(K)
      DO 2 L=1,NL
      J = ILOSS(1,K,L)
      M = ILOSS(2,K,L)
      DO 2 I=1,NZ
   2  XL(I) = XL(I) + AR(J,I)*D(M,I)
C
C   PRODUCTION RATE XP
      print 500,K
  500 format(/'K =',i3)
      print *,'In CHEMPL'

      NP = NUMP(K)
      DO 3 L=1,NP
      J = IPROD(K,L)
      M = JCHEM(1,J)
      N = JCHEM(2,J)
      DO 3 I=1,NZ
      if(i.lt.nz) go to 3
      print 501, L,AR(J,I),D(M,I),D(N,I)
  501 format(2x,'L,AR(J,NZ,D(M,NZ),D(N,NZ) =',I3,1p3e10.3)
   3  XP(I) = XP(I) + AR(J,I)*D(M,I)*D(N,I)
C
C      do i=1,nz,9
c      print *,'xp= ',xp(i)
c      end do
      RETURN
      END
