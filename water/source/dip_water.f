C***********************************************************************
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
               DIPC= +DIPX*ycos + DIPY*ysin
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
C****     SUBROUTINE DIPS    *******************************************
C***********************************************************************
C
C
         SUBROUTINE DIPS(DIPY,DIPZ,R1,R2,THETA)
C
C
C     PURPOSE:
C     -------
C     Components of water molecule dipole moment, for specified water
C     molecule geometries in the Y-Z plane, are returned using a surface
C     fit based on ab-initio points computed by Polyansky et al. (2002)
C     and augmented by relativistic corrections based on a fit of the
C     same function to ab-initio spin-orbit coupling CCSD(T)
C     calculations by Joost van Stralen.
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     NOTES:
C     -----
C     The oxygen atom is considered as fixed at the origin and the H1
C     atom lies in the positive Y-Z quadrant.  The cartesian coordinate
C     system is such that the Z-axis bisects the included bond angle for
C     any geometry.  This subroutine is intended to be called by the
C     triatom program suite of Tennyson et al. (1993).
C
C                                ^  Y
C                                l        H1
C                                l
C                                l
C                - - - - - - - - O - - - - - - - - >  Z
C                                l
C                                l
C                                l         H2
C
C
C     REFERENCES:
C     ----------
C     Schwenke D.W. & Partridge H., 2000
C     Tennyson, J., Miller, S. & Le Sueur, C.R., 1993. Comput. Phys.
C       Comm. 75, 339.
C
C
C     DICTIONARY:
C     ----------
CD    A0BOHR              -  Bohr radius in Aangstroms.
CD    BAREFR              -  reference bond angle in radians.
CD    BLREFR              -  reference bond length in atomic units.
CD    DEBYES              -  number of Debyes in an atomic unit.
CD    DIPY                -  dipole moment component orthogonal to bond
C                            bisector.
CD    DIPYAL              -  all-electron correction to dipole moment
C                            component orthogonal to bond bisector.
C23456789012345678901234567890123456789012345678901234567890123456789012
CD    DIPYRC              -  relativistic correction to dipole moment
C                            component orthogonal to bond bisector.
CD    DIPZ                -  dipole moment component along bond angle
C                            bisector.
CD    DIPZAL              -  all-electron correction to dipole moment
C                            component along bond angle bisector.
C23456789012345678901234567890123456789012345678901234567890123456789012
CD    DIPZRC              -  relativistic correction to dipole moment
C                            component along bond angle bisector.
C23456789012345678901234567890123456789012345678901234567890123456789012
CD    FCBETA              -  scaling for exponetials used in effective
C                            charge calculation.
CD    FPARAM(MPARAM)      -  fit parameters.
CI    IFIRST              -  flag signally first or later call of this
C                            routine.
CD    PRMnnn              -  parameter values as scalars.
CD    R1                  -  separation of oxygen and H1 atoms in Bohrs.
CD    R2                  -  separation of oxygen and H2 atoms in Bohrs.
CD    THETA               -  cosine of bond angle.
C
C
C     SUBROUTINE AND FUNCTION SUBPROGRAMS REQUIRED:
C     --------------------------------------------
C     DIPSAL              -  returns the all electron corrections to
C                            dipole moment components of the water
C                            molecule for specified geometry.
C     DIPSRC              -  returns relativistic corrections to dipole
C                            moment components of the water molecule for
C                            specified geometry.
C
C
C     DATE AND PROGRAMMER:
C     -------------------
C     A.E. Lynas-Gray                                     2003 March 7th
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C-----------------------------------------------------------------------
C
C     IMPORT:  R1     , R2     , THETA
C     EXPORT:  DIPY   , DIPZ
C
C-----------------------------------------------------------------------
C
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C
         PARAMETER  (A0BOHR = 0.52917706D0)
         PARAMETER  (BAREFR = 104.52D0*3.141592653589793/180.0D0)
         PARAMETER  (BLREFR = 0.957624684D0/A0BOHR)
         PARAMETER  (FCBETA = 3.8D-1*A0BOHR**2)
         PARAMETER  (DEBYES = 2.541747578)
         PARAMETER  (MPARAM = 84)
C
C
         DIMENSION  FPARAM(MPARAM)
C
C
         SAVE  IFIRST
         SAVE  PRM001 , PRM002 , PRM003 , PRM004 , PRM005 , PRM006
         SAVE  PRM007 , PRM008 , PRM009 , PRM010 , PRM011 , PRM012
         SAVE  PRM013 , PRM014 , PRM015 , PRM016 , PRM017 , PRM018
         SAVE  PRM019 , PRM020 , PRM021 , PRM022 , PRM023 , PRM024
         SAVE  PRM025 , PRM026 , PRM027 , PRM028 , PRM029 , PRM030
         SAVE  PRM031 , PRM032 , PRM033 , PRM034 , PRM035 , PRM036
         SAVE  PRM037 , PRM038 , PRM039 , PRM040 , PRM041 , PRM042
         SAVE  PRM043 , PRM044 , PRM045 , PRM046 , PRM047 , PRM048
         SAVE  PRM049 , PRM050 , PRM051 , PRM052 , PRM053 , PRM054
         SAVE  PRM055 , PRM056 , PRM057 , PRM058 , PRM059 , PRM060
         SAVE  PRM061 , PRM062 , PRM063 , PRM064 , PRM065 , PRM066
         SAVE  PRM067 , PRM068 , PRM069 , PRM070 , PRM071 , PRM072
         SAVE  PRM073 , PRM074 , PRM075 , PRM076 , PRM077 , PRM078
         SAVE  PRM079 , PRM080 , PRM081 , PRM082 , PRM083 , PRM084
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
         DATA       IFIRST / 0 /
C
C
         DATA (FPARAM(LPARAM),LPARAM=01,MPARAM)
     &   /     8.42992008D-01 ,    -4.75023568D-01 ,     5.88905752D-01
     &   ,    -9.72543418D-01 ,     4.61684495D-01 ,    -5.49846701D-02
     &   ,    -1.53001817D-03 ,    -1.67122588D-01 ,     9.51180682D-02
     &   ,    -3.91257368D-02 ,     7.51228631D-02 ,    -7.51829669D-02
     &   ,     3.28534376D-03 ,    -1.61577165D-01 ,     4.34199750D-01
     &   ,    -1.65720955D-02 ,    -1.68386132D-01 ,    -2.43027329D-01
     &   ,     4.16233838D-01 ,     8.26766491D-02 ,    -1.15110353D-01
     &   ,     2.23781373D-02 ,     7.66085744D-01 ,    -6.40375137D-01
     &   ,    -6.67070970D-02 ,     1.65917203D-02 ,    -1.05564304D-01
     &   ,     3.04661226D-03 ,     1.99622497D-01 ,    -1.90168560D-01
     &   ,     1.28372312D-02 ,    -1.30569965D-01 ,     1.22900084D-01
     &   ,     2.63448544D-02 ,    -2.10397944D-01 ,    -1.92977637D-01
     &   ,     9.75051373D-02 ,     3.68338287D-01 ,    -5.11568069D-01
     &   ,    -1.33396417D-01 ,     3.71872112D-02 ,    -4.27539706D-01
     &   ,     3.81300807D-01 ,    -2.58304298D-01 ,     1.36489511D-01
     &   ,    -3.11779082D-02 ,    -9.25601870D-02 ,    -1.46371782D-01
     &   ,     1.30923809D-02 ,     1.54430047D-01 ,     6.12576865D-02
     &   ,    -1.61683947D-01 ,     1.27978146D-01 ,     5.48809990D-02
     &   ,     1.01166859D-01 ,     3.59075889D-02 ,    -4.56623703D-01
     &   ,     5.62227607D-01 ,    -1.04735762D-01 ,     1.32874161D-01
     &   ,     1.51885718D-01 ,     5.82131624D-01 ,     1.02462292D-01
     &   ,    -7.99435899D-02 ,     3.32773179D-02 ,    -7.47206435D-02
     &   ,    -1.98611431D-02 ,    -1.70582868D-02 ,     4.67161685D-01
     &   ,    -4.82844055D-01 ,     1.92998528D-01 ,     1.25247002D-01
     &   ,    -1.27038658D-01 ,     1.16065912D-01 ,     5.34933817D-04
     &   ,    -1.85741652D-02 ,     6.17188439D-02 ,     1.71201229D-01
     &   ,    -5.97901307D-02 ,     2.57224031D-02 ,    -2.98017502D-01
     &   ,    -5.65484986D-02 ,    -2.46513247D-01 ,     4.58969548D-02/
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement functions for polynomial evaluation:
C
         POLYN1(C0,C1,X) = C0 + C1*X
         POLYN2(C0,C1,C2,X) = C0 + POLYN1(C1,C2,X)*X
         POLYN3(C0,C1,C2,C3,X) = C0 + POLYN2(C1,C2,C3,X)*X
         POLYN4(C0,C1,C2,C3,C4,X) = C0 + POLYN3(C1,C2,C3,C4,X)*X
         POLYN5(C0,C1,C2,C3,C4,C5,X) = C0 + POLYN4(C1,C2,C3,C4,C5,X)*X
         POLYN6(C0,C1,C2,C3,C4,C5,C6,X) = C0
     &                                 + POLYN5(C1,C2,C3,C4,C5,C6,X)*X
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement Function Definitions:
C
C
C     1) Displacement of hydrogen atom from equilibrium (au).
C
         RE(R) =  R - BLREFR
C
C
C     2) Fractional displacement of hydrogen atom from equilibrium.
C
         FE(R) = RE(R)/BLREFR
C
C
C     3) Difference in cosines between bond angle and equilibrium bond
C     angle.
C
         CE(C) = C - COS(BAREFR)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C     4) Compute effective charge using Schwenke & Partridge's (2000)
C        functional form:
C
         Q(R1,R2,C)
     &   = EXP(-FCBETA*RE(R1)**2)*(POLYN6(PRM001,PRM002
     &   ,PRM003,PRM004,PRM005,PRM006,PRM007,FE(R1))
     &   + EXP(-FCBETA*RE(R2)**2)
     &   *(POLYN5(PRM008,PRM009,PRM010,PRM011,PRM012,PRM013,CE(C))*CE(C)
     &   +POLYN5(PRM014,PRM015,PRM016,PRM017,PRM018,PRM019,CE(C))*FE(R2)
     &   +POLYN4(PRM020,PRM021,PRM022,PRM023,PRM024,CE(C))*FE(R2)**2
     &   +POLYN3(PRM025,PRM026,PRM027,PRM028,CE(C))*FE(R2)**3
     &   +POLYN2(PRM029,PRM030,PRM031,CE(C))*FE(R2)**4
     &   +POLYN1(PRM032,PRM033,CE(C))*FE(R2)**5
     &   +PRM034*FE(R2)**6
     &   +POLYN4(PRM035,PRM036,PRM037,PRM038,PRM039,CE(C))*CE(C)*FE(R1)
     &   +POLYN4(PRM040,PRM041,PRM042,PRM043,PRM044,CE(C))*FE(R2)*FE(R1)
     &   +POLYN3(PRM045,PRM046,PRM047,PRM048,CE(C))*FE(R2)**2*FE(R1)
     &   +POLYN2(PRM049,PRM050,PRM051,CE(C))*FE(R2)**3*FE(R1)
     &   +POLYN1(PRM052,PRM053,CE(C))*FE(R2)**4*FE(R1)
     &   +PRM054*FE(R2)**5*FE(R1)
     &   +POLYN3(PRM055,PRM056,PRM057,PRM058,CE(C))*CE(C)*FE(R1)**2
     &   +POLYN3(PRM059,PRM060,PRM061,PRM062,CE(C))*FE(R2)*FE(R1)**2
     &   +POLYN2(PRM063,PRM064,PRM065,CE(C))*FE(R2)**2*FE(R1)**2
     &   +POLYN1(PRM066,PRM067,CE(C))*FE(R2)**3*FE(R1)**2
     &   +PRM068*FE(R2)**4*FE(R1)**2
     &   +POLYN2(PRM069,PRM070,PRM071,CE(C))*CE(C)*FE(R1)**3
     &   +POLYN2(PRM072,PRM073,PRM074,CE(C))*FE(R2)*FE(R1)**3
     &   +POLYN1(PRM075,PRM076,CE(C))*FE(R2)**2*FE(R1)**3
     &   +PRM077*FE(R2)**3*FE(R1)**3
     &   +POLYN1(PRM078,PRM079,CE(C))*CE(C)*FE(R1)**4
     &   +POLYN1(PRM080,PRM081,CE(C))*FE(R2)*FE(R1)**4
     &   +PRM082*FE(R2)**2*FE(R1)**4
     &   +PRM083*CE(C)*FE(R1)**5
     &   +PRM084*FE(R2)*FE(R1)**5))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C
C     Make copies of parameter values on first call only.
C
         IF (IFIRST.EQ.0) THEN
            PRM001 = FPARAM(01)
            PRM002 = FPARAM(02)
            PRM003 = FPARAM(03)
            PRM004 = FPARAM(04)
            PRM005 = FPARAM(05)
            PRM006 = FPARAM(06)
            PRM007 = FPARAM(07)
            PRM008 = FPARAM(08)
            PRM009 = FPARAM(09)
            PRM010 = FPARAM(10)
            PRM011 = FPARAM(11)
            PRM012 = FPARAM(12)
            PRM013 = FPARAM(13)
            PRM014 = FPARAM(14)
            PRM015 = FPARAM(15)
            PRM016 = FPARAM(16)
            PRM017 = FPARAM(17)
            PRM018 = FPARAM(18)
            PRM019 = FPARAM(19)
            PRM020 = FPARAM(20)
            PRM021 = FPARAM(21)
            PRM022 = FPARAM(22)
            PRM023 = FPARAM(23)
            PRM024 = FPARAM(24)
            PRM025 = FPARAM(25)
            PRM026 = FPARAM(26)
            PRM027 = FPARAM(27)
            PRM028 = FPARAM(28)
            PRM029 = FPARAM(29)
            PRM030 = FPARAM(30)
            PRM031 = FPARAM(31)
            PRM032 = FPARAM(32)
            PRM033 = FPARAM(33)
            PRM034 = FPARAM(34)
            PRM035 = FPARAM(35)
            PRM036 = FPARAM(36)
            PRM037 = FPARAM(37)
            PRM038 = FPARAM(38)
            PRM039 = FPARAM(39)
            PRM040 = FPARAM(40)
            PRM041 = FPARAM(41)
            PRM042 = FPARAM(42)
            PRM043 = FPARAM(43)
            PRM044 = FPARAM(44)
            PRM045 = FPARAM(45)
            PRM046 = FPARAM(46)
            PRM047 = FPARAM(47)
            PRM048 = FPARAM(48)
            PRM049 = FPARAM(49)
            PRM050 = FPARAM(50)
            PRM051 = FPARAM(51)
            PRM052 = FPARAM(52)
            PRM053 = FPARAM(53)
            PRM054 = FPARAM(54)
            PRM055 = FPARAM(55)
            PRM056 = FPARAM(56)
            PRM057 = FPARAM(57)
            PRM058 = FPARAM(58)
            PRM059 = FPARAM(59)
            PRM060 = FPARAM(60)
            PRM061 = FPARAM(61)
            PRM062 = FPARAM(62)
            PRM063 = FPARAM(63)
            PRM064 = FPARAM(64)
            PRM065 = FPARAM(65)
            PRM066 = FPARAM(66)
            PRM067 = FPARAM(67)
            PRM068 = FPARAM(68)
            PRM069 = FPARAM(69)
            PRM070 = FPARAM(70)
            PRM071 = FPARAM(71)
            PRM072 = FPARAM(72)
            PRM073 = FPARAM(73)
            PRM074 = FPARAM(74)
            PRM075 = FPARAM(75)
            PRM076 = FPARAM(76)
            PRM077 = FPARAM(77)
            PRM078 = FPARAM(78)
            PRM079 = FPARAM(79)
            PRM080 = FPARAM(80)
            PRM081 = FPARAM(81)
            PRM082 = FPARAM(82)
            PRM083 = FPARAM(83)
            PRM084 = FPARAM(84)
         ENDIF
C
C
C     Ensure that parameter values are never copied again
C
         IFIRST = 1
C
C
C
C     Evaluate fit functions for current estimate of fit parameters.
C
         DIPZ  = (Q(R1,R2,theta)*R1+Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0+theta)/2.0D0)
         DIPY  = (Q(R1,R2,theta)*R1-Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0-theta)/2.0D0)
C
C
C     Dipole moment components are required in atomic units.
C
         DIPZ = DIPZ/DEBYES
         DIPY = DIPY/DEBYES
C
C
C     Augment dipole moment components by relativistic and all electron
C     corrections.
C
         CALL DIPSRC(DIPYRC,DIPZRC,R1,R2,THETA)
C
         CALL DIPSAL(DIPYAL,DIPZAL,R1,R2,THETA)
C
         DIPZ = DIPZ + DIPZAL + DIPZRC
         DIPY = DIPY + DIPYAL + DIPYRC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     END OF SUBROUTINE DIPS
C
         END
C***********************************************************************
C****     SUBROUTINE DIPSAL  *******************************************
C***********************************************************************
C
C
         SUBROUTINE DIPSAL(DIPYAL,DIPZAL,R1,R2,THETA)
C
C
C     PURPOSE:
C     -------
C     All electron - frozen core corrections are returned for components
C     (for specified water molecule geometries in the Y-Z plane) of the
C     water molecule dipole moment, based on a fit of the Schwenke &
C     Partridge (2000) functional form to ab-intio CCSD(T) calculations
C     by Attila G. Csa'sza'.
C
C
C     NOTES:
C     -----
C     The oxygen atom is considered as fixed at the origin and the H1
C     atom lies in the positive Y-Z quadrant.  The cartesian coordinate
C     system is such that the Z-axis bisects the included bond angle for
C     any geometry.  This subroutine is intended to be called by the
C     DIPS routine in the triatom program suite of Tennyson et al.
C     (1993).
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                                ^  Y
C                                l        H1
C                                l
C                                l
C                - - - - - - - - O - - - - - - - - >  Z
C                                l
C                                l
C                                l         H2
C
C
C     REFERENCES:
C     ----------
C     Schwenke D.W. & Partridge H., 2000
C     Tennyson, J., Miller, S. & Le Sueur, C.R., 1993. Comput. Phys.
C       Comm. 75, 339.
C
C
C     DICTIONARY:
C     ----------
CD    A0BOHR              -  Bohr radius in Aangstroms.
CD    BAREFR              -  reference bond angle in radians.
CD    BLREFR              -  reference bond length in atomic units.
CD    DEBYES              -  number of Debyes in an atomic unit.
CD    DIPYAL              -  all-electron correction to dipole moment
C                            component orthogonal to bond bisector.
CD    DIPZAL              -  all-electron correction to dipole moment
C                            component along bond angle bisector.
CD    FCBETA              -  scaling for exponetials used in effective
C                            charge calculation.
CD    FPARAM(MPARAM)      -  fit parameters.
CI    IFIRST              -  flag signally first or later call of this
C                            routine.
CD    PRMnnn              -  parameter values as scalars.
CD    R1                  -  separation of oxygen and H1 atoms in Bohrs.
CD    R2                  -  separation of oxygen and H2 atoms in Bohrs.
CD    THETA               -  bond angle in radians.
C
C
C     DATE AND PROGRAMMER:
C     -------------------
C     A.E. Lynas-Gray                                     2003 March 7th
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C-----------------------------------------------------------------------
C
C     IMPORT:  R1     , R2     , THETA
C     EXPORT:  DIPYAL , DIPZAL
C
C-----------------------------------------------------------------------
C
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C
         PARAMETER  (A0BOHR = 0.52917706D0)
         PARAMETER  (BAREFR = 104.52D0*3.141592653589793/180.0D0)
         PARAMETER  (BLREFR = 0.957624684D0/A0BOHR)
         PARAMETER  (FCBETA = 3.8D-1*A0BOHR**2)
         PARAMETER  (DEBYES = 2.541747578)
         PARAMETER  (MPARAM = 84)
C
C
         DIMENSION  FPARAM(MPARAM)
C
C
         SAVE  IFIRST
         SAVE  PRM001 , PRM002 , PRM003 , PRM004 , PRM005 , PRM006
         SAVE  PRM007 , PRM008 , PRM009 , PRM010 , PRM011 , PRM012
         SAVE  PRM013 , PRM014 , PRM015 , PRM016 , PRM017 , PRM018
         SAVE  PRM019 , PRM020 , PRM021 , PRM022 , PRM023 , PRM024
         SAVE  PRM025 , PRM026 , PRM027 , PRM028 , PRM029 , PRM030
         SAVE  PRM031 , PRM032 , PRM033 , PRM034 , PRM035 , PRM036
         SAVE  PRM037 , PRM038 , PRM039 , PRM040 , PRM041 , PRM042
         SAVE  PRM043 , PRM044 , PRM045 , PRM046 , PRM047 , PRM048
         SAVE  PRM049 , PRM050 , PRM051 , PRM052 , PRM053 , PRM054
         SAVE  PRM055 , PRM056 , PRM057 , PRM058 , PRM059 , PRM060
         SAVE  PRM061 , PRM062 , PRM063 , PRM064 , PRM065 , PRM066
         SAVE  PRM067 , PRM068 , PRM069 , PRM070 , PRM071 , PRM072
         SAVE  PRM073 , PRM074 , PRM075 , PRM076 , PRM077 , PRM078
         SAVE  PRM079 , PRM080 , PRM081 , PRM082 , PRM083 , PRM084
C
C
         DATA       IFIRST / 0 /
C
C
         DATA (FPARAM(LPARAM),LPARAM=01,MPARAM)
     &   /     8.87522590E-04 ,     3.81148129E-04 ,    -4.27389808E-04
     &   ,     2.82407535E-04 ,    -2.52298225E-04 ,    -3.49383481E-05
     &   ,     5.62112182E-05 ,    -1.19848824E-04 ,     3.56289493E-05
     &   ,    -2.14018764E-05 ,     1.61927292E-05 ,    -4.33738933E-05
     &   ,    -2.72791985E-05 ,     8.50149663E-05 ,    -8.21477210E-04
     &   ,    -6.41601873E-05 ,     1.83985045E-04 ,     2.07454921E-03
     &   ,    -2.53058551E-03 ,    -1.29338136E-04 ,     2.94008758E-04
     &   ,     6.48329442E-04 ,    -1.86784845E-03 ,     1.32184662E-03
     &   ,     4.18090858E-05 ,     4.03413840E-04 ,    -9.25142376E-04
     &   ,     1.07434392E-03 ,     2.44981093E-05 ,     5.04454947E-04
     &   ,    -2.55092134E-04 ,    -3.62677092E-05 ,    -6.59114856E-04
     &   ,     2.00275736E-05 ,    -1.58656767E-04 ,     7.60086928E-04
     &   ,    -1.17964961E-03 ,    -7.24685611E-04 ,     1.50384998E-03
     &   ,    -3.00624873E-04 ,     4.11824440E-04 ,     1.26975705E-04
     &   ,    -3.73087154E-04 ,     5.21674054E-04 ,    -1.54370384E-04
     &   ,     4.51390981E-04 ,    -1.02113653E-03 ,     3.04104411E-04
     &   ,    -4.70472105E-06 ,    -1.32201370E-04 ,    -2.53027800E-04
     &   ,     1.41822777E-04 ,    -1.31110733E-04 ,     3.89029519E-05
     &   ,    -8.84054054E-04 ,     2.23726631E-04 ,     1.09040184E-06
     &   ,     5.45179530E-04 ,    -3.41888168E-04 ,     8.58629995E-04
     &   ,    -7.66854500E-04 ,    -2.98567174E-04 ,     3.16300575E-05
     &   ,     7.28875020E-05 ,    -2.82031018E-04 ,    -1.26606756E-04
     &   ,     2.30710924E-04 ,    -2.06779689E-04 ,     3.86956846E-04
     &   ,     1.03125686E-03 ,    -9.61888640E-04 ,     9.75772273E-05
     &   ,     5.28461416E-04 ,    -2.76048842E-04 ,     3.86164931E-04
     &   ,     4.46979713E-04 ,     2.93890567E-04 ,     1.32845191E-04
     &   ,    -7.56434340E-04 ,     5.10016049E-04 ,     3.54172080E-04
     &   ,    -1.32894056E-04 ,     7.84503180E-04 ,     1.21441757E-04/
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement functions for polynomial evaluation:
C
         POLYN1(C0,C1,X) = C0 + C1*X
         POLYN2(C0,C1,C2,X) = C0 + POLYN1(C1,C2,X)*X
         POLYN3(C0,C1,C2,C3,X) = C0 + POLYN2(C1,C2,C3,X)*X
         POLYN4(C0,C1,C2,C3,C4,X) = C0 + POLYN3(C1,C2,C3,C4,X)*X
         POLYN5(C0,C1,C2,C3,C4,C5,X) = C0 + POLYN4(C1,C2,C3,C4,C5,X)*X
         POLYN6(C0,C1,C2,C3,C4,C5,C6,X) = C0
     &                                 + POLYN5(C1,C2,C3,C4,C5,C6,X)*X
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement Function Definitions:
C
C
C     1) Displacement of hydrogen atom from equilibrium (au).
C
         RE(R) =  R - BLREFR
C
C
C     2) Fractional displacement of hydrogen atom from equilibrium.
C
         FE(R) = RE(R)/BLREFR
C
C
C     3) Difference in cosines between bond angle and equilibrium bond
C     angle.
C
         CE(C) = C - COS(BAREFR)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C     4) Compute effective charge using Schwenke & Partridge's (2000)
C        functional form:
C
         Q(R1,R2,C)
     &   = EXP(-FCBETA*RE(R1)**2)*(POLYN6(PRM001,PRM002
     &   ,PRM003,PRM004,PRM005,PRM006,PRM007,FE(R1))
     &   + EXP(-FCBETA*RE(R2)**2)
     &   *(POLYN5(PRM008,PRM009,PRM010,PRM011,PRM012,PRM013,CE(C))*CE(C)
     &   +POLYN5(PRM014,PRM015,PRM016,PRM017,PRM018,PRM019,CE(C))*FE(R2)
     &   +POLYN4(PRM020,PRM021,PRM022,PRM023,PRM024,CE(C))*FE(R2)**2
     &   +POLYN3(PRM025,PRM026,PRM027,PRM028,CE(C))*FE(R2)**3
     &   +POLYN2(PRM029,PRM030,PRM031,CE(C))*FE(R2)**4
     &   +POLYN1(PRM032,PRM033,CE(C))*FE(R2)**5
     &   +PRM034*FE(R2)**6
     &   +POLYN4(PRM035,PRM036,PRM037,PRM038,PRM039,CE(C))*CE(C)*FE(R1)
     &   +POLYN4(PRM040,PRM041,PRM042,PRM043,PRM044,CE(C))*FE(R2)*FE(R1)
     &   +POLYN3(PRM045,PRM046,PRM047,PRM048,CE(C))*FE(R2)**2*FE(R1)
     &   +POLYN2(PRM049,PRM050,PRM051,CE(C))*FE(R2)**3*FE(R1)
     &   +POLYN1(PRM052,PRM053,CE(C))*FE(R2)**4*FE(R1)
     &   +PRM054*FE(R2)**5*FE(R1)
     &   +POLYN3(PRM055,PRM056,PRM057,PRM058,CE(C))*CE(C)*FE(R1)**2
     &   +POLYN3(PRM059,PRM060,PRM061,PRM062,CE(C))*FE(R2)*FE(R1)**2
     &   +POLYN2(PRM063,PRM064,PRM065,CE(C))*FE(R2)**2*FE(R1)**2
     &   +POLYN1(PRM066,PRM067,CE(C))*FE(R2)**3*FE(R1)**2
     &   +PRM068*FE(R2)**4*FE(R1)**2
     &   +POLYN2(PRM069,PRM070,PRM071,CE(C))*CE(C)*FE(R1)**3
     &   +POLYN2(PRM072,PRM073,PRM074,CE(C))*FE(R2)*FE(R1)**3
     &   +POLYN1(PRM075,PRM076,CE(C))*FE(R2)**2*FE(R1)**3
     &   +PRM077*FE(R2)**3*FE(R1)**3
     &   +POLYN1(PRM078,PRM079,CE(C))*CE(C)*FE(R1)**4
     &   +POLYN1(PRM080,PRM081,CE(C))*FE(R2)*FE(R1)**4
     &   +PRM082*FE(R2)**2*FE(R1)**4
     &   +PRM083*CE(C)*FE(R1)**5
     &   +PRM084*FE(R2)*FE(R1)**5))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C
C     Make copies of parameter values on first call only.
C
         IF (IFIRST.EQ.0) THEN
            PRM001 = FPARAM(01)
            PRM002 = FPARAM(02)
            PRM003 = FPARAM(03)
            PRM004 = FPARAM(04)
            PRM005 = FPARAM(05)
            PRM006 = FPARAM(06)
            PRM007 = FPARAM(07)
            PRM008 = FPARAM(08)
            PRM009 = FPARAM(09)
            PRM010 = FPARAM(10)
            PRM011 = FPARAM(11)
            PRM012 = FPARAM(12)
            PRM013 = FPARAM(13)
            PRM014 = FPARAM(14)
            PRM015 = FPARAM(15)
            PRM016 = FPARAM(16)
            PRM017 = FPARAM(17)
            PRM018 = FPARAM(18)
            PRM019 = FPARAM(19)
            PRM020 = FPARAM(20)
            PRM021 = FPARAM(21)
            PRM022 = FPARAM(22)
            PRM023 = FPARAM(23)
            PRM024 = FPARAM(24)
            PRM025 = FPARAM(25)
            PRM026 = FPARAM(26)
            PRM027 = FPARAM(27)
            PRM028 = FPARAM(28)
            PRM029 = FPARAM(29)
            PRM030 = FPARAM(30)
            PRM031 = FPARAM(31)
            PRM032 = FPARAM(32)
            PRM033 = FPARAM(33)
            PRM034 = FPARAM(34)
            PRM035 = FPARAM(35)
            PRM036 = FPARAM(36)
            PRM037 = FPARAM(37)
            PRM038 = FPARAM(38)
            PRM039 = FPARAM(39)
            PRM040 = FPARAM(40)
            PRM041 = FPARAM(41)
            PRM042 = FPARAM(42)
            PRM043 = FPARAM(43)
            PRM044 = FPARAM(44)
            PRM045 = FPARAM(45)
            PRM046 = FPARAM(46)
            PRM047 = FPARAM(47)
            PRM048 = FPARAM(48)
            PRM049 = FPARAM(49)
            PRM050 = FPARAM(50)
            PRM051 = FPARAM(51)
            PRM052 = FPARAM(52)
            PRM053 = FPARAM(53)
            PRM054 = FPARAM(54)
            PRM055 = FPARAM(55)
            PRM056 = FPARAM(56)
            PRM057 = FPARAM(57)
            PRM058 = FPARAM(58)
            PRM059 = FPARAM(59)
            PRM060 = FPARAM(60)
            PRM061 = FPARAM(61)
            PRM062 = FPARAM(62)
            PRM063 = FPARAM(63)
            PRM064 = FPARAM(64)
            PRM065 = FPARAM(65)
            PRM066 = FPARAM(66)
            PRM067 = FPARAM(67)
            PRM068 = FPARAM(68)
            PRM069 = FPARAM(69)
            PRM070 = FPARAM(70)
            PRM071 = FPARAM(71)
            PRM072 = FPARAM(72)
            PRM073 = FPARAM(73)
            PRM074 = FPARAM(74)
            PRM075 = FPARAM(75)
            PRM076 = FPARAM(76)
            PRM077 = FPARAM(77)
            PRM078 = FPARAM(78)
            PRM079 = FPARAM(79)
            PRM080 = FPARAM(80)
            PRM081 = FPARAM(81)
            PRM082 = FPARAM(82)
            PRM083 = FPARAM(83)
            PRM084 = FPARAM(84)
         ENDIF
C
C
C     Ensure that parameter values are never copied again
C
         IFIRST = 1
C
C
C
C     Evaluate fit functions for current estimate of fit parameters.
C
         DIPZAL  = (Q(R1,R2,theta)*R1+Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0+theta)/2.0D0)
         DIPYAL  = (Q(R1,R2,theta)*R1-Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0-theta)/2.0D0)
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     End of SUBROUTINE DIPSAL
C
         END
C***********************************************************************
C****     SUBROUTINE DIPSRC  *******************************************
C***********************************************************************
C
C
         SUBROUTINE DIPSRC(DIPYRC,DIPZRC,R1,R2,THETA)
C
C
C     PURPOSE:
C     -------
C     Relativistic corrections are returned to components
C     (for specified water molecule geometries in the Y-Z plane) of the
C     water molecule dipole moment, based on a fit of the Schwenke &
C     Partridge (2000) functional form to ab-intio spin-orbit coupling
C     CCSD(T) calculations by Joost van Stralen.
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     NOTES:
C     -----
C     The oxygen atom is considered as fixed at the origin and the H1
C     atom lies in the positive Y-Z quadrant.  The cartesian coordinate
C     system is such that the Z-axis bisects the included bond angle for
C     any geometry.  This subroutine is intended to be called by the
C     DIPS routine in the triatom program suite of Tennyson et al.
C     (1993).
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                                ^  Y
C                                l        H1
C                                l
C                                l
C                - - - - - - - - O - - - - - - - - >  Z
C                                l
C                                l
C                                l         H2
C
C
C     REFERENCES:
C     ----------
C     Schwenke D.W. & Partridge H., 2000
C     Tennyson, J., Miller, S. & Le Sueur, C.R., 1993. Comput. Phys.
C       Comm. 75, 339.
C
C
C     DICTIONARY:
C     ----------
CD    A0BOHR              -  Bohr radius in Aangstroms.
CD    BAREFR              -  reference bond angle in radians.
CD    BLREFR              -  reference bond length in atomic units.
CD    DEBYES              -  number of Debyes in an atomic unit.
CD    DIPYRC              -  relativistic correction to dipole moment
C                            component orthogonal to bond bisector.
C23456789012345678901234567890123456789012345678901234567890123456789012
CD    DIPZRC              -  relativistic correction to dipole moment
C                            component along bond angle bisector.
CD    FCBETA              -  scaling for exponetials used in effective
C                            charge calculation.
CD    FPARAM(MPARAM)      -  fit parameters.
CI    IFIRST              -  flag signally first or later call of this
C                            routine.
CD    PRMnnn              -  parameter values as scalars.
CD    R1                  -  separation of oxygen and H1 atoms in Bohrs.
CD    R2                  -  separation of oxygen and H2 atoms in Bohrs.
CD    THETA               -  bond angle in radians.
C
C
C     DATE AND PROGRAMMER:
C     -------------------
C     A.E. Lynas-Gray                                     2003 March 7th
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C-----------------------------------------------------------------------
C
C     IMPORT:  R1     , R2     , THETA
C     EXPORT:  DIPYRC , DIPZRC
C
C-----------------------------------------------------------------------
C
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C
         PARAMETER  (A0BOHR = 0.52917706D0)
         PARAMETER  (BAREFR = 104.52D0*3.141592653589793/180.0D0)
         PARAMETER  (BLREFR = 0.957624684D0/A0BOHR)
         PARAMETER  (FCBETA = 3.8D-1*A0BOHR**2)
         PARAMETER  (DEBYES = 2.541747578)
         PARAMETER  (MPARAM = 84)
C
C
         DIMENSION  FPARAM(MPARAM)
C
C
         SAVE  IFIRST
         SAVE  PRM001 , PRM002 , PRM003 , PRM004 , PRM005 , PRM006
         SAVE  PRM007 , PRM008 , PRM009 , PRM010 , PRM011 , PRM012
         SAVE  PRM013 , PRM014 , PRM015 , PRM016 , PRM017 , PRM018
         SAVE  PRM019 , PRM020 , PRM021 , PRM022 , PRM023 , PRM024
         SAVE  PRM025 , PRM026 , PRM027 , PRM028 , PRM029 , PRM030
         SAVE  PRM031 , PRM032 , PRM033 , PRM034 , PRM035 , PRM036
         SAVE  PRM037 , PRM038 , PRM039 , PRM040 , PRM041 , PRM042
         SAVE  PRM043 , PRM044 , PRM045 , PRM046 , PRM047 , PRM048
         SAVE  PRM049 , PRM050 , PRM051 , PRM052 , PRM053 , PRM054
         SAVE  PRM055 , PRM056 , PRM057 , PRM058 , PRM059 , PRM060
         SAVE  PRM061 , PRM062 , PRM063 , PRM064 , PRM065 , PRM066
         SAVE  PRM067 , PRM068 , PRM069 , PRM070 , PRM071 , PRM072
         SAVE  PRM073 , PRM074 , PRM075 , PRM076 , PRM077 , PRM078
         SAVE  PRM079 , PRM080 , PRM081 , PRM082 , PRM083 , PRM084
C
C
         DATA       IFIRST / 0 /
C
C
         DATA (FPARAM(LPARAM),LPARAM=01,MPARAM)
     &   /    -7.70391431D-04 ,    -4.69740422D-04 ,     6.11087089D-05
     &   ,     4.78438596D-04 ,    -1.12685202D-05 ,    -1.12065063D-04
     &   ,     2.13152598D-05 ,     2.38666544D-04 ,    -1.42028177D-04
     &   ,     2.15814616D-05 ,    -9.79817414D-05 ,     1.62372395D-04
     &   ,    -5.71789205D-05 ,    -2.40074325D-04 ,     8.57470033D-04
     &   ,    -2.04125972D-05 ,    -5.35663676D-05 ,    -2.47124489D-03
     &   ,     2.81054573D-03 ,     1.58077106D-04 ,    -8.44667229D-05
     &   ,    -1.03055337D-03 ,     2.48984690D-03 ,    -1.84013206D-03
     &   ,     4.42688142D-05 ,    -9.19489772D-04 ,     1.21105427D-03
     &   ,    -9.73901129D-04 ,     3.74132833D-05 ,    -3.08710238D-04
     &   ,     1.01051046D-04 ,    -5.23007147D-05 ,     7.02070771D-04
     &   ,     1.02900212D-05 ,     1.67727776D-05 ,    -7.98697816D-04
     &   ,     1.31938688D-03 ,     8.52813944D-04 ,    -1.56458467D-03
     &   ,     6.23775268D-05 ,     7.14333291D-05 ,    -2.67532712D-04
     &   ,     3.46721557D-04 ,    -5.16566506D-04 ,     3.01089778D-04
     &   ,    -1.02985813D-03 ,     1.36199791D-03 ,    -2.36280481D-04
     &   ,     3.13943136D-04 ,    -8.77717248D-05 ,     6.70725945D-04
     &   ,    -2.06068784D-04 ,    -5.33647544D-05 ,    -3.19712388D-04
     &   ,     1.59858353D-03 ,    -5.50020544D-04 ,     3.63587344D-04
     &   ,    -8.30255798D-04 ,     6.64299936D-04 ,    -1.03042321D-03
     &   ,     3.16815305D-04 ,     6.53067313D-04 ,    -2.09874182D-04
     &   ,    -2.80591310D-04 ,     1.78868664D-04 ,    -1.60483294D-04
     &   ,    -5.27635966D-05 ,     2.56790168D-04 ,    -1.07000204D-04
     &   ,    -1.20860175D-03 ,     8.58326152D-04 ,     2.14972752D-04
     &   ,    -5.99867373D-04 ,     1.72207569D-04 ,    -9.83031641D-05
     &   ,    -7.33574678D-04 ,    -4.10908135D-04 ,    -7.07565516D-04
     &   ,     1.60530116D-03 ,    -1.56383263D-04 ,    -2.91009696D-04
     &   ,     4.73564251D-05 ,    -1.05978909D-03 ,     9.25246641D-05/
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement functions for polynomial evaluation:
C
         POLYN1(C0,C1,X) = C0 + C1*X
         POLYN2(C0,C1,C2,X) = C0 + POLYN1(C1,C2,X)*X
         POLYN3(C0,C1,C2,C3,X) = C0 + POLYN2(C1,C2,C3,X)*X
         POLYN4(C0,C1,C2,C3,C4,X) = C0 + POLYN3(C1,C2,C3,C4,X)*X
         POLYN5(C0,C1,C2,C3,C4,C5,X) = C0 + POLYN4(C1,C2,C3,C4,C5,X)*X
         POLYN6(C0,C1,C2,C3,C4,C5,C6,X) = C0
     &                                 + POLYN5(C1,C2,C3,C4,C5,C6,X)*X
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Statement Function Definitions:
C
C
C     1) Displacement of hydrogen atom from equilibrium (au).
C
         RE(R) =  R - BLREFR
C
C
C     2) Fractional displacement of hydrogen atom from equilibrium.
C
         FE(R) = RE(R)/BLREFR
C
C
C     3) Difference in cosines between bond angle and equilibrium bond
C     angle.
C
         CE(C) = C - COS(BAREFR)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C     4) Compute effective charge using Schwenke & Partridge's (2000)
C        functional form:
C
         Q(R1,R2,C)
     &   = EXP(-FCBETA*RE(R1)**2)*(POLYN6(PRM001,PRM002
     &   ,PRM003,PRM004,PRM005,PRM006,PRM007,FE(R1))
     &   + EXP(-FCBETA*RE(R2)**2)
     &   *(POLYN5(PRM008,PRM009,PRM010,PRM011,PRM012,PRM013,CE(C))*CE(C)
     &   +POLYN5(PRM014,PRM015,PRM016,PRM017,PRM018,PRM019,CE(C))*FE(R2)
     &   +POLYN4(PRM020,PRM021,PRM022,PRM023,PRM024,CE(C))*FE(R2)**2
     &   +POLYN3(PRM025,PRM026,PRM027,PRM028,CE(C))*FE(R2)**3
     &   +POLYN2(PRM029,PRM030,PRM031,CE(C))*FE(R2)**4
     &   +POLYN1(PRM032,PRM033,CE(C))*FE(R2)**5
     &   +PRM034*FE(R2)**6
     &   +POLYN4(PRM035,PRM036,PRM037,PRM038,PRM039,CE(C))*CE(C)*FE(R1)
     &   +POLYN4(PRM040,PRM041,PRM042,PRM043,PRM044,CE(C))*FE(R2)*FE(R1)
     &   +POLYN3(PRM045,PRM046,PRM047,PRM048,CE(C))*FE(R2)**2*FE(R1)
     &   +POLYN2(PRM049,PRM050,PRM051,CE(C))*FE(R2)**3*FE(R1)
     &   +POLYN1(PRM052,PRM053,CE(C))*FE(R2)**4*FE(R1)
     &   +PRM054*FE(R2)**5*FE(R1)
     &   +POLYN3(PRM055,PRM056,PRM057,PRM058,CE(C))*CE(C)*FE(R1)**2
     &   +POLYN3(PRM059,PRM060,PRM061,PRM062,CE(C))*FE(R2)*FE(R1)**2
     &   +POLYN2(PRM063,PRM064,PRM065,CE(C))*FE(R2)**2*FE(R1)**2
     &   +POLYN1(PRM066,PRM067,CE(C))*FE(R2)**3*FE(R1)**2
     &   +PRM068*FE(R2)**4*FE(R1)**2
     &   +POLYN2(PRM069,PRM070,PRM071,CE(C))*CE(C)*FE(R1)**3
     &   +POLYN2(PRM072,PRM073,PRM074,CE(C))*FE(R2)*FE(R1)**3
     &   +POLYN1(PRM075,PRM076,CE(C))*FE(R2)**2*FE(R1)**3
     &   +PRM077*FE(R2)**3*FE(R1)**3
     &   +POLYN1(PRM078,PRM079,CE(C))*CE(C)*FE(R1)**4
     &   +POLYN1(PRM080,PRM081,CE(C))*FE(R2)*FE(R1)**4
     &   +PRM082*FE(R2)**2*FE(R1)**4
     &   +PRM083*CE(C)*FE(R1)**5
     &   +PRM084*FE(R2)*FE(R1)**5))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C
C     Make copies of parameter values on first call only.
C
         IF (IFIRST.EQ.0) THEN
            PRM001 = FPARAM(01)
            PRM002 = FPARAM(02)
            PRM003 = FPARAM(03)
            PRM004 = FPARAM(04)
            PRM005 = FPARAM(05)
            PRM006 = FPARAM(06)
            PRM007 = FPARAM(07)
            PRM008 = FPARAM(08)
            PRM009 = FPARAM(09)
            PRM010 = FPARAM(10)
            PRM011 = FPARAM(11)
            PRM012 = FPARAM(12)
            PRM013 = FPARAM(13)
            PRM014 = FPARAM(14)
            PRM015 = FPARAM(15)
            PRM016 = FPARAM(16)
            PRM017 = FPARAM(17)
            PRM018 = FPARAM(18)
            PRM019 = FPARAM(19)
            PRM020 = FPARAM(20)
            PRM021 = FPARAM(21)
            PRM022 = FPARAM(22)
            PRM023 = FPARAM(23)
            PRM024 = FPARAM(24)
            PRM025 = FPARAM(25)
            PRM026 = FPARAM(26)
            PRM027 = FPARAM(27)
            PRM028 = FPARAM(28)
            PRM029 = FPARAM(29)
            PRM030 = FPARAM(30)
            PRM031 = FPARAM(31)
            PRM032 = FPARAM(32)
            PRM033 = FPARAM(33)
            PRM034 = FPARAM(34)
            PRM035 = FPARAM(35)
            PRM036 = FPARAM(36)
            PRM037 = FPARAM(37)
            PRM038 = FPARAM(38)
            PRM039 = FPARAM(39)
            PRM040 = FPARAM(40)
            PRM041 = FPARAM(41)
            PRM042 = FPARAM(42)
            PRM043 = FPARAM(43)
            PRM044 = FPARAM(44)
            PRM045 = FPARAM(45)
            PRM046 = FPARAM(46)
            PRM047 = FPARAM(47)
            PRM048 = FPARAM(48)
            PRM049 = FPARAM(49)
            PRM050 = FPARAM(50)
            PRM051 = FPARAM(51)
            PRM052 = FPARAM(52)
            PRM053 = FPARAM(53)
            PRM054 = FPARAM(54)
            PRM055 = FPARAM(55)
            PRM056 = FPARAM(56)
            PRM057 = FPARAM(57)
            PRM058 = FPARAM(58)
            PRM059 = FPARAM(59)
            PRM060 = FPARAM(60)
            PRM061 = FPARAM(61)
            PRM062 = FPARAM(62)
            PRM063 = FPARAM(63)
            PRM064 = FPARAM(64)
            PRM065 = FPARAM(65)
            PRM066 = FPARAM(66)
            PRM067 = FPARAM(67)
            PRM068 = FPARAM(68)
            PRM069 = FPARAM(69)
            PRM070 = FPARAM(70)
            PRM071 = FPARAM(71)
            PRM072 = FPARAM(72)
            PRM073 = FPARAM(73)
            PRM074 = FPARAM(74)
            PRM075 = FPARAM(75)
            PRM076 = FPARAM(76)
            PRM077 = FPARAM(77)
            PRM078 = FPARAM(78)
            PRM079 = FPARAM(79)
            PRM080 = FPARAM(80)
            PRM081 = FPARAM(81)
            PRM082 = FPARAM(82)
            PRM083 = FPARAM(83)
            PRM084 = FPARAM(84)
         ENDIF
C
C
C     Ensure that parameter values are never copied again
C
         IFIRST = 1
C
C
C
C     Evaluate fit functions for current estimate of fit parameters.
C
         DIPZRC  = (Q(R1,R2,theta)*R1+Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0+theta)/2.0D0)
         DIPYRC  = (Q(R1,R2,theta)*R1-Q(R2,R1,theta)*R2)
     &   *                                SQRT((1.0D0-theta)/2.0D0)
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     End of SUBROUTINE DIPSRC
C
         END
