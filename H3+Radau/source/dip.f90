      SUBROUTINE  dipd(d0,r1,r2,xcos,nu)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      common /mass/ xmass(3),g1,g2,zembed,zbisc   

      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/,PI/3.1415927D0/


      IF (G1 .EQ. X0) THEN
!        BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
         Q1 = R1
         Q2 = R2
         cost=xcos
      ELSE IF (G2 .EQ. X0) THEN
!        SCATTERING COORDINATES: ATOM 2 = ATOM 3
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         ALPHA= ACOS(XCOS)
         IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
            Q3 = ABS(XX - R2)
            Q2 = (YY + R2)
            COST = -X1
         ELSE IF (XCOS .LE. (TINY - X1)) THEN
            Q3 = (XX + R2)
            Q2 = ABS(YY - R2)
            COST = X1
         ELSE
            Q3 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
            Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
            COST = (Q1**2 + Q2**2 - R1**2) / (X2 * Q1 * Q2)
         ENDIF
         XSIN= SQRT(1.0D0 - XCOS*XCOS)
         BETA= ASIN(XSIN*YY/Q2)
         Q1=R1
      ELSE
!        GENERAL COORDINATES (INCLUDING RADAU): ATOM 1 = ATOM 2
         F1= X1/G1
         F2= X1/G2
         F12= X1 - F1*F2
         P1= R1*(X1-F1)/(G2*F12)
         P2= R2*(X1-F2)/(G1*F12)
         S1= R1-P1
         S2= R2-P2
         Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOS)/(X1-G1)
         Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOS)/(X1-G2)
         Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOS)
         COST = (Q1*Q1 + Q2*Q2 - Q3*Q3)/(X2*Q1*Q2)
!         ALPHA= ASIN(SIN(THETA)*Q2/R2)
!         BETA = ASIN(SIN(THETA)*Q1/R1)
      ENDIF
      THETA = ACOS(COST)

      RR1=Q1
      RR2=SQRT((Q3**2+Q2**2)*0.5d0-Q1**2*0.25d0) 
      xc = -(q3**2-RR2**2-RR1**2*0.25d0)/(RR1*RR2)
      call DIPD0(DIPCZ,DIPCX,RR1,RR2,XC)
      
      if (zbisc) then

         alpha=acos(xcos)
         theta=acos(xc)

         tg=g1*R1*sin(alpha)/(g1*R1*cos(alpha)+R2)
         gamma=atan(tg)
         tgb=(r1-g1*r2)*tan(alpha*0.5d0)/(r1+g1*r2)
         beta=-atan(tgb)

            DIPCA = DIPCX*sin(beta) + DIPCZ * cos(beta)
            DIPCB = DIPCX*cos(beta) + DIPCZ * sin(-beta)

         IF(NU.EQ.0) THEN
            ! perpendicular axis (z)
            D0= DIPCB
         ELSE
            ! bisecting axis (x)
            D0= DIPCA
         ENDIF

         else

            if (zembed) then
               stop
            else
         IF(NU.EQ.0) THEN
            D0 = DIPCZ
         ELSE
            D0 = DIPCX
         ENDIF
            end if
      end if

      return
    end SUBROUTINE dipd

      SUBROUTINE DIPD0(DIPCX,DIPCY,R1,R2,XCOS)
!
!     THIS SUBROUTINE USES THE DIPOLE SURFACE GIVEN BY
!     ROESHE ET AL (J. CHEM. PHYS. 101, 2231 (1994)). 
!     CALCULATE THE COMPONENTS OF THE H3+ DIPOLE.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      DIMENSION D(7), DZ(7), DX(7), R(3), B(3), xmass(3)
!jayesh, 14/3/2002
!       common /mass/ xmass(3),g1,g2,zembed,zbisc   

      DATA BETA/1.3000d0/, RE/1.6500d0/
      DATA D/-0.487805d0,-0.249076d0,-0.151777d0,-0.014064d0,&
      -0.005687d0,-0.008458d0,-0.033219d0/
! 
      DATA ONE/1.0D0/, TWO/2.0D0/, THREE/3.0D0/, X0/0.0D0/
!
      SR2= SQRT(TWO)
      SR3= SQRT(THREE)
      SR6= SR2*SR3
      xmass=1.d0
      
      G1L=xmass(2)/(xmass(2)+xmass(3))
!      G1L=0.5d0

      R22= R2*R2
      R1C= R1*G1L
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
!
      DZ(1)= SZ
      DX(1)= SX
      DO 2 I=2,4
      II= I-1
      DZ(I)= DZ(II)*SA
      DX(I)= DX(II)*SA
2     CONTINUE
!      DX(5)= SZ2 - SX2
!      DZ(5)= TWO*SZ*SX
      DZ(5)= SZ2 - SX2
      DX(5)= -TWO*SZ*SX
      DX(6)= DX(5)*SA
      DZ(6)= DZ(5)*SA
      DX(7)= SE2*SX
      DZ(7)= SE2*SZ
!
      DIPZ= X0
      DIPX= X0
      DO 3 I=1,7
      DIPZ= DIPZ + D(I)*DZ(I)
      DIPX= DIPX + D(I)*DX(I)
3     CONTINUE
!
      DIPZ= DIPZ*R2*SR2/SR3  - RDEL
      DIPX= DIPX*R1/SR2
!     this line is in a test for Jim Watson
!     dipx= -dipx
!
      XSIN= SQRT(ONE- XCOS*XCOS)

            DIPCX= - (DIPX + XCOS*DIPZ) ! z-axis (diatom)
            DIPCY= DIPZ*XSIN            ! x-axis (perp)

!write(98,*)dipcx,dipcy
      RETURN
      END
