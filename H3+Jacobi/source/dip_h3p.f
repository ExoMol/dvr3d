      SUBROUTINE DIPD(DIPC,R1,R2,XCOS,NU)
C
C     THIS SUBROUTINE USES THE DIPOLE SURFACE GIVEN BY
C     ROESHE ET AL (J. CHEM. PHYS. 101, 2231 (1994)). 
C     CALCULATE THE COMPONENTS OF THE H3+ DIPOLE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      DIMENSION D(7), DZ(7), DX(7), R(3), B(3)
!jayesh, 14/3/2002
!      COMMON /LOGIC/  ZMORSE, ZNCOR, ZPRINT, ZLPOT, ZTRA, ZSTART, ZEMBED
!      COMMON/MASS/ XMASS(3), G1, G2
       common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
       common /mass/ xmass(3),g1,g2,zembed,zbisc   

      DATA BETA/1.3000/, RE/1.6500/
      DATA D/-0.487805,-0.249076,-0.151777,-0.014064,
     1-0.005687,-0.008458,-0.033219/
C 
      DATA ONE/1.0D0/, TWO/2.0D0/, THREE/3.0D0/, X0/0.0D0/
C
      SR2= SQRT(TWO)
      SR3= SQRT(THREE)
      SR6= SR2*SR3
      R22= R2*R2
      R1C= R1*G1
      RC2= R1 - R1C
      R1C2= R1C*R1C
      RC22= RC2*RC2
      T2= TWO*R1C*R2*XCOS
      T3= TWO*RC2*R2*XCOS
      R(1)= R1
      R(2)= SQRT(R22 + R1C2 - T2)
      R(3)= SQRT(R22 + RC22 + T3)
      BTOT= X0
      XMTOT= X0
      DO 1 M=1,3
      XMTOT= XMTOT + XMASS(M)
      X= R(M)/RE - ONE
      Y= EXP(-BETA*X)
      B(M)= (ONE - Y)/BETA
      BTOT= BTOT + B(M)
1     CONTINUE
      RCM= R2*XMASS(1)/XMTOT
      RCC= R2/THREE
      RDEL= RCM - RCC
      SA= BTOT/SR3
      SZ= (TWO*B(1) - B(2) - B(3))/SR6
      SX= (B(3) - B(2))/SR2
      SZ2= SZ*SZ
      SX2= SX*SX
      SE2= SZ2 + SX2
C
      DZ(1)= SZ
      DX(1)= SX
      DO 2 I=2,4
      II= I-1
      DZ(I)= DZ(II)*SA
      DX(I)= DX(II)*SA
2     CONTINUE
      DZ(5)= SZ2 - SX2
      DZ(6)= DZ(5)*SA
      DZ(7)= SE2*SZ
      DX(5)= -TWO*SZ*SX
      DX(6)= DX(5)*SA
      DX(7)= SE2*SX
C
      DIPZ= X0
      DIPX= X0
      DO 3 I=1,7
      DIPZ= DIPZ + D(I)*DZ(I)
      DIPX= DIPX + D(I)*DX(I)
3     CONTINUE
C
      DIPZ= DIPZ*R2*SR2/SR3  - RDEL
      DIPX= DIPX*R1/SR2
c     this line is in a test for Jim Watson
c     dipx= -dipx
c
      XSIN= SQRT(ONE- XCOS*XCOS)
      IF(.NOT.ZEMBED) GOTO 4
      IF(NU.EQ.0) THEN
         DIPC= DIPZ + XCOS*DIPX
      ELSE
         DIPC= DIPX*XSIN
      ENDIF
      GOTO 5
4     CONTINUE
      IF(NU.EQ.0) THEN
         DIPC= -(DIPX + XCOS*DIPZ)
      ELSE
         DIPC= DIPZ*XSIN
      ENDIF
5     CONTINUE
C
C     DIPC= ONE
      RETURN
      END
C=========================================================== dipcal ===========
C     SUBROUTINE dipd0(R1,R2,XCOS,DIPX,DIPZ)
C written by CRLS on 4th March 1991
C============================================================== DIPD
C This subroutine written by SM7
      SUBROUTINE DIPD0(DIPC,R1,R2,XCOS,NU)
C
C     This subroutine uses the dipole surface given by
C     Botswina et al (J. chem. Phys. 84, 891 (1986)) to
C     calculate the components of the H3+ dipole.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
!      COMMON/MASS/ XMASS(3), G1, G2
      DIMENSION D(7), DZ(7), DX(7), R(3), B(3)
      common /mass/ xmass(3),g1,g2,zembed,zbisc   
      double precision x(3),xe(3),z(3),ze(3)
      DATA BETA/1.3000/, RE/1.6504/
      DATA D/-0.487805,-0.249076,-0.151777,-0.014064,
     1-0.005687,-0.008458,-0.033219/
      DATA ONE/1.0D0/, TWO/2.0D0/, THREE/3.0D0/, X0/0.0D0/
      G1=0.5D0
      XMASS(1)=1D0
      XMASS(2)=1D0
      XMASS(3)=1D0
C
      SR2= SQRT(TWO)
      SR3= SQRT(THREE)
      SR6= SR2*SR3
      R22= R2*R2
      R1C= R1*G1
      RC2= R1 - R1C
      R1C2= R1C*R1C
      RC22= RC2*RC2
      T2= TWO*R1C*R2*XCOS
      T3= TWO*RC2*R2*XCOS
      R(1)= R1
      R(2)= SQRT(R22 + R1C2 - T2)
      R(3)= SQRT(R22 + RC22 + T3)
      BTOT= X0
      XMTOT= X0
      DO 1 M=1,3
      XMTOT= XMTOT + XMASS(M)
      Xold= R(M)/RE - ONE
      Y= EXP(-BETA*Xold)
      B(M)= (ONE - Y)/BETA
      BTOT= BTOT + B(M)
1     CONTINUE
      RCM= R2*XMASS(1)/XMTOT
      RCC= R2/THREE
      RDEL= RCM - RCC
      SA= BTOT/SR3
      SZ= (TWO*B(1) - B(2) - B(3))/SR6
      SX= (B(3) - B(2))/SR2
      SZ2= SZ*SZ
      SX2= SX*SX
      SE2= SZ2 + SX2
C
      DZ(1)= SZ
      DX(1)= SX
      DO 2 I=2,4
      II= I-1
      DZ(I)= DZ(II)*SA
      DX(I)= DX(II)*SA
2     CONTINUE
      DZ(5)= SZ2 - SX2
      DZ(6)= DZ(5)*SA
      DZ(7)= SE2*SZ
      DX(5)= TWO*SZ*SX
      DX(6)= DX(5)*SA
      DX(7)= SE2*SX
C
      DIPZ= X0
      DIPX= X0
      DO 3 I=1,7
      DIPZ= DIPZ + D(I)*DZ(I)
      DIPX= DIPX + D(I)*DX(I)
3     CONTINUE
C
C This gives the dipoles along the axes defined by Botswina et al
C Note that their x and z axes are not mutually orthogonal.
      DIPZ= DIPZ*R2*SR2/SR3  - RDEL
      DIPX= DIPX*R1/SR2
      XSIN= SQRT(ONE- XCOS*XCOS)
C Then convert them to axes where z lies along the scattering coordinate
      DIPZP= DIPZ + XCOS*DIPX
      DIPXP= DIPX*XSIN
C End of bit written by SM7.
C
c finally a conversion to 'lab-fixed' axes...
c first get cartesian components for atom positions in scattering coordinates
      x(1)=0d0
      z(1)=2d0*r2/3d0
      x(2)=-r1*xsin/2d0
      z(2)=-r2/3d0-r1*xcos/2d0
      x(3)=r1*xsin/2d0
      z(3)=-r2/3d0+r1*xcos/2d0
c the equilibrium positions in scattering coordinates are given as follows
      xe(1)=0d0
      ze(1)=re/dsqrt(3d0)
      xe(2)=-re/dsqrt(3d0)
      ze(2)=-re/dsqrt(3d0)/2d0
      xe(3)=re/dsqrt(3d0)
      ze(3)=-re/dsqrt(3d0)/2d0
c the angle is then defined by a cross product \sum_i{\bf r}e_i x M{\bf r}_i=0
c between the equilibrium positions {\bf r}e and the rotated positionsM{\bf r}
      top=0d0
      bottom=0d0
      do 100 i=1,3
        top   =top   +ze(i)*x(i)-xe(i)*z(i)
        bottom=bottom+xe(i)*x(i)+ze(i)*z(i)
100   continue
      tang=top/bottom
      g=atan(tang)
      cosg=cos(g)
      sing=sin(g)
c the dipole can now be referred to these axes
      IF(NU.NE.0) dipC=dipxp*cosg-dipzp*sing
      IF(NU.EQ.0) dipC=dipxp*sing+dipzp*cosg
      RETURN
      END
