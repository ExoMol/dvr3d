      SUBROUTINE DIPD(DIPC,R1,R2,XCOS,NU)
C
C     TRANSFORM GENERALISED COORDINATES TO THOSE FOR PARTICULAR
C     SYSTEM. THIS VERSION TRANSFORMS TO AB2 BONDLENGTH-BONDANGLE
C     COORDINATES. ALLOWANCE MUST BE MADE FOR THE NUMBERING OF THE ATOMS
c     Additionally, the zbisc option is included.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /MASS/ XMASS(3),G1,G2,zembed,zbisc
C
C     (R = r . S = r'. T = theta)
C
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/,PI/3.1415927D0/
C
      IF (G1 .EQ. X0) THEN
C        BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
         Q1 = R1
         Q2 = R2
         THETA = ACOS(XCOS)
      ELSE IF (G2 .EQ. X0) THEN
C        SCATTERING COORDINATES: ATOM 2 = ATOM 3
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         ALPHA= ACOS(XCOS)
         IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
            Q1 = ABS(XX - R2)
            Q2 = (YY + R2)
            COST = -X1
         ELSE IF (XCOS .LE. (TINY - X1)) THEN
            Q1 = (XX + R2)
            Q2 = ABS(YY + R2)
            COST = X1
         ELSE
            Q1 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
            Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
            COST = (Q1**2 + Q2**2 - R1**2) / (X2 * Q1 * Q2)
         ENDIF
         XSIN= SQRT(1.0D0 - XCOS*XCOS)
         THETA = ACOS(COST)
         BETA= ASIN(XSIN*YY/Q2)
      ELSE
C        GENERAL COORDINATES (INCLUDING RADAU): ATOM 1 = ATOM 2
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
         THETA = ACOS(COST)
      ENDIF
C
c      CALL DIPS(DIPX,DIPY,Q1,Q2,THETA)
      CALL DIPS(DIPX,DIPY,Q1,Q2,cost)
C
C     BONDLENGTH-BONDANGLE CO-ORDINATES
C
      IF (G1.EQ.X0) THEN
         GAMMA= THETA/X2
         ycos= COS(GAMMA)
         ysin= SIN(GAMMA)
         IF (ZEMBED) THEN
            IF (NU.EQ.0) THEN
               DIPC= +DIPY*ycos - DIPX*ysin
            ELSE
               DIPC= +DIPX*ycos + DIPY*ysin
            ENDIF
         ELSE
            if (NU.EQ.0) THEN
               DIPC= +DIPY*ycos + DIPX*ysin
            ELSE
               DIPC= -DIPX*ycos + DIPY*ysin
            ENDIF
         ENDIF
C
C     SCATTERING CO-ORDINATES
C
      ELSE if (G2.EQ.X0) THEN
         GAMMA= BETA - THETA/x2
         if (ZEMBED) THEN
            ycos= COS(GAMMA)
            ysin= SIN(GAMMA)
            if (NU.EQ.0) THEN
               DIPC= -DIPX*ysin - DIPY*ycos
            ELSE
               DIPC= +DIPX*ycos - DIPY*ysin
            ENDIF
         ELSE
            DELTA= ALPHA - GAMMA
            ycos= COS(DELTA)
            ysin= SIN(DELTA)
            if (NU.EQ.0) THEN
               DIPC= +DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= +DIPX*ycos - DIPY*ysin
            ENDIF
         ENDIF
C
C     ALL OTHER CO-ORDINATES
C
      ELSE
         IF (ZBISC) THEN
            H1 = G1*Q1
            COSA = (R2*R2 + Q2*Q2 - H1*H1)/(X2*R2*Q2)
            alpha = (acos(xcos)-theta)/x2 - acos(cosa)
            ycos= - COS(ALPHA)
            ysin= + SIN(ALPHA)
            if (NU.EQ.1) THEN
               DIPC= +DIPX*ysin - DIPY*ycos
            ELSE
               DIPC= -DIPX*ycos - DIPY*ysin
            ENDIF
         ELSEIF (ZEMBED) THEN
            H1 = G1*Q1
            COSA = (R2*R2 + H1*H1 - Q2*Q2)/(X2*R2*H1)
            ALPHA= ACOS(COSA)
            ycos= - COS(ALPHA + THETA/X2)
            ysin= + SIN(ALPHA + THETA/X2)
            if (NU.EQ.0) THEN
               DIPC= -DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= -DIPX*ycos - DIPY*ysin
            ENDIF
         ELSE
            H2 = G2*Q2
            COSA = (R1*R1 + H2*H2 - Q1*Q1)/(X2*R1*H2)
            alpha = ACOS(COSA)
            ycos= - COS(alpha + THETA/X2)
            ysin= + SIN(alpha + THETA/X2)
            if (NU.EQ.0) THEN
               DIPC= +DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= -DIPX*ycos + DIPY*ysin
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END
