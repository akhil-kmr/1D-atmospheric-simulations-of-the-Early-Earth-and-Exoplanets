      SUBROUTINE DIFCO(O2,CO2)

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comAERBLK.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'
C
      DIMENSION WGT(NQT),HISCALE(NQT),BKMIG(NQT),HD(NQT,NZ),WT(NZ)
     2 ,BKMG(NZ),H(NZ),O2(NZ),CO2(NZ)
C
c      DATA WGT/30., 16., 18., 17., 32., 34., 48.,  1.,  2., 16.,
c     2          28., 48., 47., 44., 30., 46., 47., 63., 79., 62.,
c     3         108.,102., 50., 52., 35., 51., 36., 97., 34., 33.,
c     4          48., 64., 98., 49.,1.E9/
C
c a typo in the old code. 33 instead of 32.
c In May 17 2019, JL moves CO2 into long lived 
C      DATA WGT/30., 16., 18., 17., 33.,34., 48.,  1.,  2., 16.,
C     2          28., 48., 47., 44., 30., 46., 47., 63., 79., 62.,
c               co                  NO                       NO3
C     3         108., 32., 34., 33.,48., 64., 98., 49.,44.,1.E9/
c                    o2                               co2
C-AJ09/01/2022 ADD CL BACK
      DATA WGT/30., 16., 18., 17., 33.,34., 48.,  1.,  2., 16.,
     2          28., 48., 47., 44., 30., 46., 47., 63., 79., 62.,
c               co                  NO                       NO3
     3         108., 102., 50., 52., 35., 51., 36., 97.,32., 34., 
c                                                       o2    
     4          33.,48., 64., 98., 49.,44.,1.E9/
c                                     co2         
      DZ = DELZ
      NZ1 = NZ - 1
      DO J=1,NZ
      WT(J) = O2(J)*32. + CO2(J)*44. + (1.-O2(J)-CO2(J))*28.
      BKMG(J) = 1.38E-16/(1.67E-24*WT(J)*G)
      END DO

C   COMPUTE DIFFUSION LIFETIME AT EACH HEIGHT (H*H/K)
      DO 1 J=1,NZ
       H(J) = BKMG(J) * T(J)
       HSCALE(J) = H(J)
   1   TAUEDD(J) = H(J)*H(J)/EDD(J)
C
C   CALCULATE MOLECULAR DIFFUSION COEFFICIENTS
      DO 2 I=1,NQT
       DO 3 J=1,NZ
   3  DI(I,J) = 1.52E18*(1./WGT(I)+ 1./WT(J))**0.5 *(T(J)**0.5/DEN(J))
   2  CONTINUE
C
C  Write certain diffusion coefficients to a file
      write(87,101)
 101  format(5x,'zkm',6x,'EDD',6x,'DI(LH,J)',2x,'DI(LH2,J)')
      DO J=1,NZ
      zkm = Z(J)/1.e5
      write(87,100) zkm,EDD(J),DI(LH,J),DI(LH2,J)
 100  format(1x,1p4e10.3)
      ENDDO
C
C ***** HD(I,J) = Hi*N at grid step J *****
      DO 4 I=1,NQT
       DO 5 J=1,NZ
        BKMIG(I) = 1.38E-16/(1.67E-24*WGT(I)*G)
        HISCALE(I) = BKMIG(I)*T(J)
        HI(I,J) = DI(I,J)*(1./HISCALE(I) - 1./HSCALE(J))
   5    HD(I,J) = HI(I,J)*DEN(J)
   4  CONTINUE
C
C ***** DK(I,J) = (K+D)*N AT GRID STEP J+1/2 *****
      DO 6 J=1,NZ1
       EDDAV = SQRT(EDD(J)*EDD(J+1))
       DENAV = SQRT(DEN(J)*DEN(J+1))
       DO 7 I=1,NQT
        DK(I,J) = (EDDAV+DI(I,J))*DENAV
   7  CONTINUE
   6  CONTINUE
C
C  STORE CONSTANT JACOBIAN COEFFICIENTS
C
      DO 11 I=1,NQT
      DHU(I,1) = HD(I,2)/DEN(1)/(2.*DZ)
      DHL(I,NZ) = HD(I,NZ1)/DEN(NZ)/(2.*DZ)
      DO 11 J=2,NZ1
      DHU(I,J) = HD(I,J+1)/DEN(J)/(2.*DZ)
  11  DHL(I,J) = HD(I,J-1)/DEN(J)/(2.*DZ)
C
      DZ2 = DZ*DZ
      DO 8 I=1,NQT
       DU(I,1) = DK(I,1)/DEN(1)/DZ2
       DL(I,NZ) = DK(I,NZ1)/DEN(NZ)/DZ2
      DO 8 J=2,NZ1
       DU(I,J) = DK(I,J)/DEN(J)/DZ2
       DL(I,J) = DK(I,J-1)/DEN(J)/DZ2
   8   DD(I,J) = DU(I,J) + DL(I,J)
C

      RETURN
      END
