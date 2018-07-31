C     Dummy main program
      PROGRAM JRH
      CALL DVR3D
      STOP
      END
      SUBROUTINE DVR3D
C
C     PROGRAM               D V R 3 D
C                           ~~~~~~~~~
C     PROGRAM TO DO RO-VIBRATIONAL CALCULATIONS ON TRIATOMIC SYSTEMS
C     USING GENERAL LENGTH, LENGTH, ANGLE COORDINATES AND A CHOICE OF
C     EMBEDDINGS, IN A MULTIDIMENSIONAL DVR IN CHOICE OF COORDINATE
C     ORDERS.
C     THERE ARE various options:
C     (A) CALCULATION AT A FROZEN ANGLE                       ZTWOD  = T
C     (B) ATOM - RIGID DIATOM CALCULATIONS                    NCOORD = 2
C     (C) TRIATOMIC CALCULATIONS IN ALL 3 DIMENSIONS          NCOORD = 3
C     Additionally, for Radau coordinates there is the option to
C     (D) place the z-axis along the bisector                 ZBISC  = T
C     See:
C
C     J.R.HENDERSON & J.TENNYSON, CHEM.PHYS.LETT., 173, 133 (1990).
C     J.R.HENDERSON, PHD THESIS, UNIVERSITY OF LONDON (1990).
C     J.R.HENDERSON, J.TENNYSON & B.T. SUTCLIFFE, J.CHEM.PHYS. 98, 7191 (1993).
C     J. Tennyson & B.T. Sutcliffe, Int. J. Quantum Chem. 42, 941 (1992).
C
C     USE AS FOLLOWS:
C     COMMENTS ON NAMELIST PARAMETERS (& DEFAULTS) IN BLOCK DATA
C     INPUT DATA READ IN SUBROUTINES INSIZE & SETCON
C     THE PROGRAM NEEDS THE FOLLOWING SUBROUTINES:
C     1. (A) IF ZLPOT = .FALSE. (DEFAULT)
C            SUBROUTINE POTV(V,R1,R2,X) SHOULD RETURN THE POTENTIAL V
C            IN HARTREES FOR THE POINT X = COS(THETA) AND
C            BONDLENGTHS R1 & R2 IN BOHR.
C        (B) IF ZLPOT = .TRUE.
C            SUBROUTINE POT(V0,VL,R1,R2) SHOULD RETURN THE LEGENDRE
C            EXPANSION OF THE POTENTIAL IN HARTREE FOR BONDLENGTHS R1 &
C            R2 IN BOHR. V0 IS THE ISOTROPIC TERM AND VL(DIMENSION LPOT)
C            GIVES THE ANISOTROPIC CONTRIBUTIONS.
C        IF (NCOORD .EQ. 2) R1 IS THE FIXED BONDLENGTH
C     2. F02ABF TO DO IN CORE DIAGONALISATION (NAG SUPPLIED ROUTINE).
C     4. ASSEMBLER ROUTINES GTMAIN & CPUSED
C        (FOR WHICH DUMMY FORTRAN ROUTINES PROVIDED - SEE COMMENTS).
C     THE PROGRAMME WORKS IN **** ATOMIC UNITS ***** :
C     1. THE SUBROUTINE POT OR POTV SHOULD RETURN THE POTENTIAL IN
C        HARTREES FOR DISTANCES IN BOHR.
C     2. ALL INPUT IN SETCON IS IN HARTREE OR BOHR EXCEPT
C     3. THE NUCLEAR MASSES ARE READ IN ATOMIC MASS UNITS & CONVERTED.
C     4. THE EIGENVALUES ARE PRINTED IN BOTH HARTREE & WAVENUMBERS.
C
C     THE PROGRAM REQUIRES ABOUT 180K EXCLUDING THE ARRAY
C     SPACE ASKED FOR IN GTMAIN.
C
C***********************************************************************
C
C     THIS IS A DUMMY MAIN PROGRAMME WHICH                         #001
C     HANDLES DYNAMIC ARRAY ALLOCATION.
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      EXTERNAL DYNAM
C
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      NAMELIST/PRT/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DATA NALLOC/0/
C
      WRITE(6,1000)
 1000 FORMAT(1H1,9X,'PROGRAM DVR3D (VERSION OF 15 Feb 1994 (am) )')
C     READ IN NAMELIST INPUT DATA (DEFAULTS IN BLOCK DATA)
      READ(5,PRT)
C
C     READ IN CONTROL PARAMETERS OF PROBLEM.
C
      CALL INSIZE
C
C     DETERMINE REQUIRED STORAGE AREA.
C
      CALL CORE(NCORE)
C
      CALL GTMAIN(DYNAM,NCORE,NALLOC)
 
      WRITE(6,1010)
 1010 FORMAT(/,10X,'INSUFFICIENT STORAGE ******* STOP')
      STOP
      END
      BLOCK DATA
C     STORES DEFAULTS FOR NAMELIST PARAMETERS                       #002
      IMPLICIT LOGICAL (Z)
C
C
C  ZPHAM[F] = T requests printing of the Hamiltonian matrix.
C  ZPRAD[F] = T requests printing of the radial matrix elements.
C  ZP1D [F] = T requests printing of the results of 1D calculations.
C  ZP2D [F] = T requests printing of the results of 2D calculations.
C  ZPMIN[F] = T requests only minimal printing.
C  ZPVEC[F] = T requests printing of the eigenvectors.
C  ZLMAT[F] = T requests printing of the L-matrix.
C  ZCUT[F]  = T final dimension selected using an energy cut-off given
C             by EMAX2.
C           = F final dimension determined by NHAM3.
C  ZMORS1[T]= T use Morse oscillator-like functions for r_1 coordinate;
C           = F use spherical oscillator functions.
C  ZMORS2[T]= T use Morse oscillator-like functions for r_2 coordinate;
C           = F use spherical oscillator functions.
C  ZROT[F]  = F DO VIBRATIONAL PART OF ROTATIONAL CALCULATION BY
C               LOOPING OVER K
C  ZEMBED[T]= T Z AXIS IS ALONG R2, = F Z AXIS IS ALONG R1.
C  ZBISC[F] = T Z axis is along bisector of R1 & R2 (IDIA=-2 only).
C  ZLIN     = T Forces suppresion of functions at last DVR point
C               (ZBISC=T only).
C  ZLADD[T] = T ICREMENT LMAX BY 1 FOR (J,K) BLOCK
C           = F LMAX CONSTANT (=F has a bug), (USED ONLY IF ZROT = .TRUE.)
C  ZLPOT[F] = T potential supplied in POT;
C           = F potential supplied in POTV.
C  ZVEC[F]  = T store the eigenvectors from all the parts of the calculation
C             (1D,2D and 3D) on stream IOUT2.
C             Further information relating to this (arrays IV1 and IV2) is
C             stored on stream IOUT1.
C  ZALL[F]  = T requests no truncation of the intermediate solution.
C  ZTHETA[T]= T let theta be first in the order of solution;
C           = F let theta be last in the order of solution.
C  ZR2R1[T] = T let r_2 come before r_1 in the order of solution;
C           = F let r_1 come before r_2 in the order of solution.
C  ZTRAN[F]= T perform the transformation of the solution coefficients
C             to the expression for the wavefunction amplitudes at the grid
C             points. Store the data on stream IWAVE, ZTRAN = T
C             automatically sets ZVEC = T for non-bisector case.
C  ZQUAD2[T]= T use the DVR quadrature approximation for the integrals of
C             the r_2^{-2} matrix, and hence make its DVR transformation
C             diagonal.
C           = F evaluate the r_2^{-2} integrals fully and perform the
C             full DVR transformation on them.
C             Note that ZQUAD2 = F is not implemented for ZMORS2 = T
C             or for ZTHETA = F.
C  ZPFUN[T] = T store energy levels on stream ILEV
C
C  ILEV[14]     stream for final eigenvalues (formatted).
C  IEIGS1[7]    stream for eigenvalues of the 1D solutions.
C  IVECS1[3]    stream for eigenvectors of the 1D solutions.
C  IEIGS2[2]    stream for eigenvalues of the 2D solutions.
C  IVECS2[4]    stream for eigenvectors of the 2D solutions.
C  IVINT[17]    a scratch stream used for storing intermediate vectors in
C               building the final Hamiltonian.
C  INTVEC[16]   a scratch stream for intermediate storage of the 2D vectors.
C  IBAND[15]    scratch file used for storing bands of the final Hamiltonian.
C  IOUT1[24]    stream for arrays IV1 and IV2, which record the sizes of
C               the truncated vctors. Used when ZVEC = T.
C  IOUT2[25]    stream for the 1D, 2D and 3D vectors for use when ZVEC = T.
C  IWAVE[26]    stores the wavefunction amplitudes at the grid points when
C               ZTRAN = T.
C
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DATA ZPHAM/.FALSE./,ZPRAD/.FALSE./,ZPVEC/.FALSE./,ZROT/.FALSE./,
     1     ZLADD/.true./,ZEMBED/.TRUE./,ZMORS2/.TRUE./,ZLPOT/.FALSE./,
     2     ZPMIN/.FALSE./,ZVEC/.FALSE./,ZQUAD2/.TRUE./,ZCUT/.FALSE./,
     3     ZDIAG/.TRUE./,ZLMAT/.FALSE./,ZWBLK/.TRUE./,ZALL/.FALSE./,
     3     ZP1D/.FALSE./,ZP2D/.FALSE./,ZR2R1/.TRUE./,ZTHETA/.TRUE./,
     3     ZMORS1/.TRUE./,ZTRAN/.FALSE./,ZTWOD/.FALSE./,zbisc/.false./,
     4     IEIGS1/7/,IVECS1/3/,IEIGS2/2/,IVECS2/4/,IVINT/17/,
     5     IBAND/15/,INTVEC/16/,IDIAG1/20/,IDIAG2/21/,IOUT1/24/,
     5     IOUT2/25/,IWAVE/26/,zlin/.false./,ZPFUN/.FALSE./,ILEV/14/
      END
      SUBROUTINE INSIZE
C
C     SET UP COMMON /SIZE/ & WRITE CONTROL PARAMETERS OF PROBLEM    #003
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
C
C     COMMON /SIZE/ STORES CONTROL PARAMETERS FOR THE PROBLEM
C     NPNT1: NUMBER OF (GAUSS-LAGUERRE) DVR POINTS IN R1
C     NPNT2: NUMBER OF (GAUSS-LAGUERRE) DVR POINTS IN R2
C     NMAX1: MAX ORDER OF R1 RADIAL LAGUERRE POLYNOMIAL ( = NPNT1-1)
C     NMAX2: MAX ORDER OF R2 RADIAL LAGUERRE POLYNOMIAL ( = NPNT2-1)
C     MAXLEG:MAX ORDER OF ANGULAR LEGENDRE POLYNOMIAL   ( = NALF -1)
C     MAX2D :UPPER BOUND ON SIZE OF INTERMEDIATE 2D HAMILTONIAN
C     MAX2D2:max2d for smaller block (zbisc=T only)
C     MAX3D :UPPER BOUND ON SIZE OF FULL 3D HAMILTONIAN
C     MAX3D2:max3d for smaller block (zbisc=T only)
C     NALF : NUMBER OF (GAUSS-LEGENDRE) DVR POINTS IN THETA
C     IDVR : NUMBER OF UNIQUE DVR POINTS
C     JROT : JROT IS TOTAL ANGULAR MOMENTUM OF THE MOLECULE
C     KMIN : ZROT=T, KMIN=1 SYM. ROT. BASIS, =0 ANTI-SYM.
C                    KMIN=2 loop over both sym & anti-sym (ZBISC=T Only)
C            ZROT=F, KMIN=FIXED VALUE OF K
C     NPNT : MAX(NPNT1,NPNT2)
C     IDIA : = 1 SCATTERING COORDINATES HETERONUCLEAR DIATOMIC
C            = 2 SCATTERING COORDINATES HOMONUCLEAR DIATOMIC
C            =-1 RADAU      COORDINATES HETRONUCLEAR DIATOMIC
C            =-2 RADAU      COORDINATES HOMONUCLEAR  DIATOMIC
C     IPAR : PARITY OF BASIS - IF IDIA=+/-2: IPAR=0 FOR EVEN & =1 FOR ODD
C     NLIM1: =NMAX1+1*(NMAX1+1+1)/2
C     NLIM2: =NMAX2+1*(NMAX2+1+1)/2
C     NPNTA: THE NUMBER OF DVR POINTS IN
C            THE COORDINATE TO BE TREATED FIRST IN THE DVR SUCCESSIVE
C            DIAGONALISATION-TRUNCATION PROCEDURE
C     NPNTB: THE NUMBER OF DVR POINTS IN THE COORDINATE TO COME SECOND
C     NPNTC: THE NUMBER OF DVR POINTS IN THE COORDINATE TO COME LAST
C     NDIMA: SET EQUAL TO NPNTA AT THE START - USED FOR DIMENSIONING
C     NDIMB: SET EQUAL TO NPNTB AT THE START - USED FOR DIMENSIONING
C     NDIMC: SET EQUAL TO NPNTC AT THE START - USED FOR DIMENSIONING
C     NEVAL: NUMBER OF EIGENVALUES WHICH HAVE TO ACTUALLY BE SUPPLIED
C            AS OUTPUT
C     NCOORD: NUMBER OF VIBRATIONAL COORDINATES EXPLICITLY CONSIDERED
C     IF (NCOORD .NE. 3) SOME OF THE ABOVE ARE DUMMIES, SEE BELOW.
c     max1dv: controls storage space for 1d vectors
C
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      CHARACTER*8 TITLE(9)
C     READ IN CONTROL PARAMETERS OF PROBLEM:
C     NCOORD = 2: ATOM-DIATOM PROBLEM WITH DIATOM RIGID
C     NCOORD = 3: FULL 3-D TRIATOMIC PROBLEM
      READ(5,5) NCOORD
    5 FORMAT(14I5)
C     OTHER PARAMTERS: SEE COMMENTS ABOVE ON /SIZE/ FOR MEANING
C     IF NCOORD=2: ALSO NEED LMAX,LPOT,IDIA,KMIN
C     IF NCOORD=3: ALL PARAMTERS REQUIRED
C
      READ(5,5) NPNT2,JROT,NEVAL,NALF,max2d,max3d,IDIA,
     1          KMIN,NPNT1,IPAR,max3d2,max1dv
      if (jrot .eq. 0) then
         ZEMBED = .TRUE.
         KMIN = 0
         ZROT = .FALSE.
         zbisc = .false.
      endif
      if (idia .eq. -2) then
         ncoord=3
         npnt1=npnt2
         IDVR=NALF
         zmors1=zmors2
         ztheta=.false.
         if (zrot) zbisc=.true.
      else
          if (ztran) ZVEC=.TRUE.
      endif
C
      IF (ZTWOD) THEN
         NCOORD = 3
         IDIA = 1
         NALF = 1
         NALF2= 1
         IDVR = 1
         ZROT=.FALSE.
         ZTHETA=.FALSE.
         max3d=max2d
         neval=min(max2d,neval)
         GOTO 887
      ENDIF
C
C
C     THE GAUSS-LEGENDRE INTEGRATION POINTS ARE ONLY ACUALLY
C     CALCULATED FOR THE HALF-RANGE:
      NALF2 = (NALF+1)/2
C     ARE WE DOING ATOM, ATOM-RIGID DIATOM OR THE FULL PROBLEM?
      IF (NCOORD .EQ. 2) THEN
C        ATOM - RIGID DIATOM CASE, SET DUMMY /SIZE/
         IF (ZR2R1) THEN
            WRITE(6,1010)
 1010 FORMAT(/10X,'ATOM - RIGID DIATOM VIBRATIONAL ANALYSIS WITH:'/)
            NPNT1 = 1
            NMAX1 = 0
         ELSE
            WRITE(6,1011)
 1011 FORMAT(/10X,'FIXED - R2 VIBRATIONAL ANALYSIS WITH:'/)
            NPNT2 = 1
            NMAX2 = 0
         ENDIF
      ELSE
C        FULL CASE: PRINT DATA ABOUT EXTRA RADIAL BASIS
         NCOORD  = 3
         WRITE(6,1020) NPNT1
 1020    FORMAT(/10X,'FULL TRIATOMIC VIBRATIONAL PROBLEM WITH'/
     1       /10X,I5,3X,'RADIAL R1 DVR POINTS USED,')
      ENDIF
C
  887 CONTINUE
C     MAXIMUM ORDER OF POLYNOMIALS;
      MAXLEG = NALF - 1
      NMAX2  = NPNT2- 1
      NMAX1  = NPNT1- 1
C
      IF (ZTHETA) THEN
          NPNTA = NALF/MAX(1,IDIA)
          IF (ZR2R1) THEN
             NPNTB = NPNT2
             NPNTC = NPNT1
          ELSE
             NPNTB = NPNT1
             NPNTC = NPNT2
          ENDIF
      ELSE
          IF (ZR2R1) THEN
             NPNTA = NPNT2
             NPNTB = NPNT1
          ELSE
             NPNTA = NPNT1
             NPNTB = NPNT2
          ENDIF
          NPNTC = NALF/MAX(1,IDIA)
      ENDIF
C
      if (idia .le. -2) then
         if (zrot) then
            max2d  = npnt1*(npnt1+1)/2
            max2d2 = npnt1*(npnt1-1)/2
         else
            MAX2D = npnt1*(npnt1+1-(ipar * 2))/2
            max2d2 = 0
         endif
         IF (.not. ZALL)  MAX3D=MIN(MAX3D,max2d*nalf)
         IF (ZALL) max3d=max2d*nalf
         if (zrot) then
            if (max3d2 .le. 0) max3d2=max3d
            max3d2=min(MAX3D,max2d2*nalf,max3d2)
         endif
         if (max1dv .le. 0) max1dv=max2d*nalf
      else
         IF (ZALL) THEN
            MAX2D = NPNTA*NPNTB
            MAX3D = MAX2D*NPNTC
         ELSE
            MAX2D=MIN(MAX2D,NPNTA*NPNTB)
            MAX3D=MIN(MAX3D,NPNTA*NPNTB*NPNTC)
         ENDIF
         max1dv=max2d
      ENDIF
C
      IF (NEVAL .LE. 0) NEVAL = 10
      NEVAL=MIN(MAX3D,NEVAL)
      IF (ZTWOD) WRITE(6,1023) NPNT1
 1023 FORMAT(/10X,I5,3X,'RADIAL R1 DVR POINTS USED,')
      IF (NCOORD .EQ. 3) WRITE(6,1030) NPNT2,2,NALF,NEVAL,MAX3D
      IF (NCOORD .EQ. 2 .AND. ZR2R1)
     1  WRITE(6,1030) NPNT2,2,NALF,NEVAL,MAX3D
      IF(NCOORD .EQ. 2 .AND. .NOT. ZR2R1)
     1  WRITE(6,1030) NPNT1,1,NALF,NEVAL,MAX3D
 1030 FORMAT(10X,I5,3X,'RADIAL R',I1,' DVR POINTS USED,',
     2      /10X,I5,3X,'ANGULAR DVR POINTS USED, WITH',
     4      /10X,I5,3X,'LOWEST EIGENVECTORS CHOSEN FROM',
     5      /5X,'UP TO',I5,3X,'DIMENSION SECULAR PROBLEM'/)
      READ(5,500)   TITLE
  500 FORMAT(9A8)
      WRITE(6,1040) TITLE
 1040 FORMAT(10X,'TITLE: ',9A8/)
      IF (NCOORD .EQ. 3) THEN
         IF (ZMORS1)  WRITE(6,1050) 1
         IF (.NOT. ZMORS1) WRITE(6,1060) 1
      ENDIF
      IF (NCOORD .EQ. 2 .AND. .NOT. ZR2R1) THEN
         IF (ZMORS1)  WRITE(6,1050) 1
         IF (.NOT. ZMORS1) WRITE(6,1060) 1
      else
         IF (ZMORS2) WRITE(6,1050) 2
         IF (.NOT. ZMORS2) WRITE(6,1060) 2
 1050 FORMAT(10X,'MORSE OSCILLATORS USED FOR R',I1,' BASIS')
 1060 FORMAT(10X,'SPHERICAL OSCILLATORS USED FOR R',I1,' BASIS')
         IF (ZQUAD2) THEN
            WRITE(6,1051)
 1051 FORMAT(/10X,'QUADRATURE APPROXIMATION USED FOR THE R2**(-2) TERMS'
     *       /)
         ELSE
            WRITE(6,1052)
 1052 FORMAT(/10X,'QUADRATURE APPROXIMATION ABANDONED FOR R2**(-2) ',
     *       'TERMS'/)
         ENDIF
      ENDIF
      IF (ZALL) WRITE(6,1067)
 1067 FORMAT(/10X,'ALL SOLUTIONS FROM LOWER DIMENSIONS HAVE BEEN USED')
      IF (ZTHETA) THEN
         IF (ZR2R1) WRITE(6,1042)
 1042    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: THETA -> R2 -> R1')
         IF (.NOT. ZR2R1) WRITE(6,1043)
 1043    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: THETA -> R1 -> R2')
      ELSE
         IF (.NOT.ZTWOD) THEN
            IF (ZR2R1) WRITE(6,1044)
 1044    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: R2 -> R1 -> THETA')
            IF (.NOT. ZR2R1) WRITE(6,1045)
 1045    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: R1 -> R2 -> THETA')
         ELSE
            IF (ZR2R1) WRITE(6,1046)
 1046    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: R2 -> R1')
            IF (.NOT. ZR2R1) WRITE(6,1047)
 1047    FORMAT(10X,'PROBLEM SOLVED IN THE ORDER: R1 -> R2')
         ENDIF
      ENDIF
      IF (ZCUT) THEN
         WRITE(6,1061)
 1061    FORMAT(10X,'FINAL BASIS SELECTED USING ENERGY CUT-OFF')
      ELSE
         IF (ZROT .AND. ZBISC) THEN
            WRITE(6,1062) max3d,max3d2
 1062    FORMAT(/10X,'Final basis comprises',I5,' lowest functions',
     1             ' for EVEN parity Hamiltonian'/
     2           10X,'Final basis comprises',I5,' lowest functions',
     3             ' for ODD  parity Hamiltonian')
         ELSE
           IF (.NOT. ZCUT) WRITE(6,1064) MAX3D
 1064      FORMAT(10X,'FINAL BASIS COMPRISES',I5,' LOWEST FUNCTIONS')
         ENDIF
      ENDIF
      IF (.NOT.ZTWOD) THEN
         IF (ZLMAT) WRITE(6,1065)
 1065    FORMAT(/10X,'PRINTING OF L-MATRIX REQUESTED')
         IF (.NOT.ZLMAT) WRITE(6,1066)
 1066    FORMAT(/10X,'PRINTING OF L-MATRIX NOT REQUESTED')
      ENDIF
      IF (ZP1D) WRITE(6,1071)
 1071 FORMAT(10X,'PRINTING OF 1D EIGENVALUES REQUESTED')
      IF (ZP2D) WRITE(6,1072)
 1072 FORMAT(10X,'PRINTING OF 2D EIGENVALUES REQUESTED')
      IF (ZPHAM) WRITE(6,1070)
 1070 FORMAT(10X,'PRINTING OF HAMILTONIAN MATRIX REQUESTED')
      IF (.NOT.ZPHAM) WRITE(6,1080)
 1080 FORMAT(10X,'PRINTING OF HAMILTONIAN MATRIX NOT REQUESTED')
      IF (ZPRAD) WRITE(6,1090)
 1090 FORMAT(10X,'PRINTING OF RADIAL MATRIX ELEMENTS REQUESTED')
      IF (.NOT.ZPRAD) WRITE(6,1100)
 1100 FORMAT(10X,'PRINTING OF RADIAL MATRIX ELEMENTS NOT REQUESTED')
      IF (ZPVEC) WRITE(6,1110)
 1110 FORMAT(10X,'PRINTING OF EIGENVECTORS REQUESTED')
      IF (.NOT.ZPVEC) WRITE(6,1120)
 1120 FORMAT(10X,'PRINTING OF EIGENVECTORS NOT REQUESTED')
      IF (ZVEC) THEN
         WRITE(6,1130) IOUT2
         WRITE(6,1131) IOUT1
 1130 FORMAT(/10X,'EIGENVALUES & VECTORS WRITTEN TO STREAM IOUT2 ='
     1       ,I4)
 1131 FORMAT(10X,'IMPORTANT INFORMATION ALSO WRITTEN TO STREAM IOUT1 ='
     1       ,I4)
         OPEN(UNIT=IOUT1, FORM='UNFORMATTED')
         OPEN(UNIT=IOUT2, FORM='UNFORMATTED')
      ENDIF
      IF (ZTRAN) THEN
         WRITE(6,1132) IWAVE
 1132 FORMAT(/10X,'WAVEFUNCTION AMPLITUDES WRITTEN TO STREAM IWAVE ='
     1       ,I4)
         OPEN(UNIT=IWAVE, FORM='UNFORMATTED')
      ENDIF
      IF (ABS(JROT) .GT. 1) ZPFUN=.FALSE.
      IF (ZPFUN) THEN
         OPEN(UNIT=ILEV,FORM='FORMATTED')
         IF (JROT .EQ. 0 .AND. MOD(IPAR,2) .EQ. 0) THEN
C           HEADER ON FILE ILEV
            WRITE(ILEV,500) TITLE
            WRITE(6,1134) ILEV
 1134 FORMAT(/10X,'EIGENVALUES      WRITTEN TO START OF STREAM ILEV ='
     1       ,I4)
         ELSE
C           POSITION FILE ILEV
  200       READ(ILEV,*,END=210,ERR=210)
            GOTO 200
  210       CONTINUE
C ******  INCLUSION OF THE FOLLOWING CARD IS MACHINE DEPENDENT ******
            BACKSPACE ILEV
            WRITE(6,1135) ILEV
 1135 FORMAT(/10X,'EIGENVALUES      WRITTEN AT END   OF STREAM ILEV ='
     1       ,I4)
         ENDIF
      ENDIF
C
      IF (IDIA .GT. 0) WRITE(6,1140)
 1140 FORMAT(/10X,'CALCULATION PERFORMED IN SCATTERING COORDINATES')
      IF (IDIA .le. 0) WRITE(6,1150)
 1150 FORMAT(/10X,'CALCULATION PERFORMED IN RADAU COORDINATES')
C
      IF (ZTWOD) GOTO 886
C
      IF (abs(IDIA) .EQ. 2) THEN
         WRITE(6,1180)
 1180    FORMAT(/10X,'DIATOMIC ASSUMED HOMONUCLEAR')
         IF (IPAR .EQ. 1) THEN
            WRITE(6,1190)
 1190       FORMAT(10X,'ODD PARITY FUNCTIONS IN BASIS SET')
         ELSE IF (IPAR .EQ. 0) THEN
            WRITE(6,1200)
 1200       FORMAT(10X,'EVEN PARITY FUNCTIONS IN BASIS SET')
         ELSE
            WRITE(6,1205)
 1205       FORMAT(10X,'ILLEGAL VALUE OF IPAR FOR IDIA = +/-2: STOP')
            STOP
         ENDIF
         if (idia .eq. 2) then
            IDVR=NALF2
            IF (2*IDVR .NE. NALF) GOTO 960
         endif
      ELSE
         WRITE(6,1210)
 1210    FORMAT(/10X,'DIATOMIC ASSUMED HETRONUCLEAR')
         IDVR=NALF
         ipar=0
      ENDIF
      IF (JROT .NE. 0) THEN
         JROT=ABS(JROT)
         IF (ZROT) THEN
            IF (KMIN .NE. 0 .AND. .NOT. ZBISC) KMIN=1
            WRITE(6,1220)
 1220 FORMAT(/10X,'***  VIBRATIONAL PART OF ROT-VIB CALCULATION  ***')
            WRITE(6,1260) JROT
            IF (KMIN .EQ. 1) THEN
               WRITE(6,1270)
 1270 FORMAT(17X,'WITH SYMMETRIC |Jk> + |J-k> FUNCTIONS IN BASIS')
            ELSE IF (KMIN .EQ. 0) THEN
               WRITE(6,1280)
 1280 FORMAT(17X,'WITH ANTI-SYMMETRIC |Jk> - |J-k> FUNCTIONS IN BASIS')
            ELSE
               KMIN=2
               WRITE(6,1275)
 1275 FORMAT(17X,'LOOP OVER SYMMETRIC & ANTI-SYMMETRIC |Jk> FUNCTIONS')
            ENDIF
            IF (ZLADD) WRITE(6,1240)
 1240 FORMAT(15X,'NALF TO BE KEPT CONSTANT WITH K')
            IF (.NOT. ZLADD) WRITE(6,1250)
 1250 FORMAT(15X,'NALF TO DECREASE WITH K')
         ELSE
            WRITE(6,1230) JROT,KMIN
 1230 FORMAT(25X,'J =',I3,'  K =',I3,
     1       /10X,'***  OPTION TO NEGLECT CORIOLIS INTERACTIONS  ***')
         ENDIF
         if (zbisc) then
            zembed=.false.
            write(6,1330)
 1330 FORMAT(/10X,'Z AXIS EMBEDDED ALONG THE BISCETOR OF R1 AND R2')
            IF (zlin) WRITE(6,1340)
 1340 FORMAT(/10X,'Removal of functions with theta = 0 enforced')
         else
            IF (ZEMBED) WRITE(6,1290) 2
 1290       FORMAT(/10X,'Z AXIS EMBEDED ALONG THE R',I1,' COORDINATE')
            IF (.NOT. ZEMBED) WRITE(6,1290) 1
         endif
      ELSE
C        CASE J = 0
         WRITE(6,1260) JROT
 1260    FORMAT(/10X,'J =',I3,' ROTATIONAL STATE')
      ENDIF
C
  886 CONTINUE
C
C     CHECK INPUT PARAMETERS ARE CONSISTENT
      NPNT = MAX(NPNT1,NPNT2)
      NMAX = MAX(NMAX1,NMAX2)
C     DIMENSION OF SQUARE MATRICES STORED IN TRIANGULAR FORM :
      NLIM1 = NPNT1 * (NPNT1+1) / 2
      NLIM2 = NPNT2 * (NPNT2+1) / 2
C
C     Declare DVR sizes for dimensioning the arrays
      NDIMA=NPNTA
      if (idia .le. -2 .and. .not. zrot) ndima=ndima-ipar
      NDIMB=NPNTB
      NDIMC=NPNTC
C
C     STORE PARAMETERS ON DISK FILES REQUESTED
C
      IF (ZVEC)
     1 WRITE (IOUT2) IDIA,IPAR,NPNTA,NPNTB,NPNTC,max2d,max3d,NEVAL
C
      IF (ZMORS2 .AND. .NOT. ZQUAD2) GOTO 961
      IF (.NOT. ZTHETA .AND. .NOT. ZQUAD2) GOTO 962
      IF (idia .eq. -2 .AND. .NOT. ZQUAD2) GOTO 963
 
      RETURN
  960 WRITE(6,970)
  970 FORMAT(//6X,'** NALF MUST BE EVEN WHEN IDIA=2: STOP **')
      STOP
  961 WRITE(6,972)
  972 FORMAT(//6X,'** CAN''T HAVE ZQUAD2 = F WITH ZMORS2 = T : STOP **',
     1        /6X,'               (NOT YET IMPLEMENTED)               ')
      STOP
  962 WRITE(6,973)
  973 FORMAT(//6X,'** CAN''T HAVE ZQUAD2 = F WITH ZTHETA = F : STOP **',
     1        /6X,'               (NOT YET IMPLEMENTED)               ')
      STOP
  963 WRITE(6,974)
  974 FORMAT(//6X,'** CAN''T HAVE ZQUAD2 = F WITH IDIA = -2 : STOP **',
     1        /6X,'               (NOT YET IMPLEMENTED)               ')
      STOP
      END
      SUBROUTINE CORE(NCORE)
C
C     CALCULATE STORAGE REQUIRED IN WORDS (= 8 BYTES)               #005
C     DYNAM: ENTRY POINT FOR LAYOUT OF DYNAMICALLY ALLOCATED STORAGE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      DOUBLE PRECISION I(NALLOC)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      SAVE
      DATA LNAG/2/
C
C     CALCULATE STORAGE POINTERS
C
C
C     ARRAY                STARTS AT         LENGTH (IN WORDS)
C
C     HBL1                   I0                NLIM1
C     HBL2                   I1                NLIM2
C     R2M2                   I2                NLIM2
C     R1M2                   I3                NLIM1
C     WT1                    I4                NPNT1
C     WT2                    I5                NPNT2
C     WALF                   I6                IDVR
C     XALF                   I7                IDVR
C     DNORM1                 I8                NPNT1
C     DNORM2                 I9                NPNT2
C     Y1                     I10               NPNT1
C     Y2                     I11               NPNT2
C     B                      I12               NPNT   +  1
C     C                      I13               NPNT   +  1
C     B1                     I14               NALF
C     CCC1                   I15               NALF
C     PLEG                   I16               (MAXLEG + 1)  *  IDVR
C     PNORM                  I17               (MAXLEG + 1)
C     HAM1                   I18               ndima  *  ndima
C     EIG1                   I19               ndima
C     WORK1                  I20               ndima  *  LNAG
C     BASS1                  I21               NPNT1  *  NPNT1
C     BASS2                  I22               NPNT2  *  NPNT2
C     HAM2                   I23               MAX2D  *  MAX2D
C if symm radau.
C     VECS1D                 I24               max1dv *  ndima
C     EIGS1D                 I25               max1dv
C else
C     VECS1D                 I24               MAX2D  *  NPNTA
C     EIGS1D                 I25               MAX2D
C
C     EIG2                   I26               MAX2D
C     WORK2                  I27               MAX2D  *  LNAG
C
C     HBAND                  I0                MAX2D  *  MAX3D
C
C     CINT                   I28               NPNTA  *  NPNTB  *  MAX2D
C     CINTP                  I29               NPNTA  *  NPNTB  *  MAX2D
C     XLMATR                 I30               IDVR   *  IDVR
C
C     XK1                    I31               NPNT1  *  NPNT1
C     XK2                    I32               NPNT2  *  NPNT2
C     R1M2T                  I33               NPNT1  *  NPNT1
C     R2M2T                  I34               NPNT2  *  NPNT2
C     R1                     I35               NPNT1
C     R2                     I36               NPNT2
C
C     HAM3                   I0                MAX3D  *  MAX3D
C     IV1                    J1                NPNTB  *  NPNTC
C     IV2                    J2                NPNTC
C if symm radau
C     vecs2d                 j7                max3d  *  max2d
C
C     EIGS2D                 J3                MAX3D
C     EVAL                   J4                MAX3D
C     WORK3                  J5                MAX3D  *  LNAG
C     NDIM2D                 J5                NPNTC
C     EIGS2                  J4                MAX2D  *  NPNTC
C
C     and if .ZTRAN. then we want
C     IV1L                   L1                NPNTB  *  NPNTC
C     IV2L                   L2                NPNTC
C     NDIM2L                 L3                NPNTC
C     VECS1L                 L4                MAX2D  *  NPNTA
C     VECS2L                 L5                MAX2D  *  MAX2D
C     VECS3L                 L6                MAX3D
C     PHI                    L7                NPNTA  *  NPNTB
C     PSI                    L8                IDVR   *  NPNT1  *  NPNT2
C     EVALL                  L9                NEVAL
C
C
      ilen=1
      if (idia .eq. -2) ilen=0
 
      I0   =   1
      I1   =   I0   +   NLIM1
      I2   =   I1   +   NLIM2 * ilen
      I3   =   I2   +   NLIM2
      I4   =   I3   +   NLIM1 * ilen
      I5   =   I4   +   NPNT1
      I6   =   I5   +   NPNT2 * ilen
      I7   =   I6   +   IDVR
      I8   =   I7   +   IDVR  * ilen
      I9   =   I8   +   NPNT1
      I10  =   I9   +   NPNT2 * ilen
      I11  =   I10  +   NPNT1
      I12  =   I11  +   NPNT2 * ilen
      I13  =   I12  +   NPNT + 1
      I14  =   I13  +   NPNT + 1
      I15  =   I14  +   NALF
      I16  =   I15  +   NALF
      I17  =   I16  +   (MAXLEG+1) * IDVR
      I18  =   I17  +   (MAXLEG+1)
      I19  =   I18  +   ndima  *  ndima * ilen
      I20  =   I19  +   ndima * ilen
      I21  =   I20  +   ndima  *  LNAG
      I22  =   I21  +   NPNT1  *  NPNT1
      I23  =   I22  +   NPNT2  *  NPNT2 * ilen
      I24  =   I23  +   MAX2D  *  MAX2D * ilen
      I25  =   I24  +   MAX2D *  NPNTA * ilen
      I26  =   I25  +   MAX2D * ilen
      I27  =   I26  +   MAX2D * ilen
      I28  =   I27  +   MAX2D  *  LNAG
C
      K1   =   I0   +   MAX2D  *   MAX3D *  ilen
      I28  =   MAX(K1,I28)
C
      I29  =   I28  +   NPNTA  *  NPNTB  *  MAX2D *  ilen
      I30  =   I29  +   NPNTA  *  NPNTB  *  MAX2D *  ilen
      I31  =   I30  +   IDVR   *  IDVR * ilen
C
      J1   =   I0   +   MAX3D*MAX3D
      J1   =   MAX(J1,I31)
      if (idia .gt. -2) then
         J2   =   J1   + NPNTB * NPNTC
         J3   =   J2   + NPNTC
      else
         J2   =   J1   + NPNTB * NPNTC * 2
         J3   =   J2   + NPNTC * 2
      endif
      J4   =   J3   +   MAX3D
      IF (ZCUT) THEN
         J5   =   J4   +   MAX3D
      ELSE
         J5   =   J4   +   MAX(MAX3D,MAX2D*NPNTC)
      ENDIF
      J6   =   J5   +   MAX(MAX3D*LNAG,NPNTC)
      J7   =   J6
C
      L1   =   I0
      L2   =   L1   +   NPNTB  *  NPNTC  *  ilen
      L3   =   L2   +   NPNTC  *  ilen
      L4   =   L3   +   NPNTC  *  ilen
      L5   =   L4   +   MAX2D  *  NPNTA  *  ilen
      L6   =   L5   +   MAX2D  *  MAX2D  *  ilen
      L7   =   L6   +   MAX3D  *  ilen
      L8   =   L7   +   MAX3D  *  NPNTA  *  NPNTB  *  ilen
      L9   =   L8   +   IDVR   *  NPNT1  *  NPNT2  *  ilen
      L10  =   L9   +   NEVAL  *  ilen
c
      I31  = J6
      IF (ZTRAN) I31  = MAX(J6,L10)
      I32  =   I31  +   NPNT1  *  NPNT1
      I33  =   I32  +   NPNT2  *  NPNT2 * ilen
      I34  =   I33  +   NPNT1  *  NPNT1
      I35  =   I34  +   NPNT2  *  NPNT2 * ilen
      I36  =   I35  +   NPNT1
      I37  =   I36  +   NPNT2 * ilen
      NCORE  =  I37
      if (idia .eq. -2) then
        i7  =  i37
        i18 =  i7 + idvr
        i19 =  i18 + ndima * ndima
        i23 =  i19 + ndima
        i24 =  i23 + max2d * max2d
        i25 =  i24 + max1dv * ndima
        i26 =  i25 + max1dv
        i30 =  i26 + max2d
        j7  =  i30 + idvr  * idvr
        l8  =  j7  + max2d * max3d
        ncore = l8 + max2d * nalf
      endif
C
      WRITE(6,1020)  NCORE
 1020 FORMAT(/,10X,'RUN REQUIRES',I8,' WORDS')
      RETURN
      ENTRY DYNAM(I,NCORE,NALLOC)
      WRITE(6,1030) NALLOC
 1030 FORMAT(/,10X,'ALLOCATED   ',I8,' WORDS')
      IF (NALLOC .LT. NCORE) RETURN
C
      CALL CCMAIN
     &   (I(I0),I(I1),I(I2),I(I3),I(I4),I(I5),I(I6),I(I7),I(I8),I(I9),
     &    I(I10),I(I11),I(I12),I(I13),I(I14),I(I15),I(I16),I(I17),
     &    I(I18),I(I19),I(I20),
     &    I(I21),I(I22),I(I23),I(I24),I(I25),I(I26),I(I27),
     &    I(I0),I(I28),I(I29),I(I30),I(I31),I(I32),I(I33),I(I34),
     &    I(I35),I(I36),
     &    I(I0),I(J1),I(J2),i(j7),I(J3),I(J4),I(J5),I(J5),I(J4),
     &    I(L1),I(L2),I(L3),I(L4),I(L5),I(L6),I(L7),I(L8),I(L9))
      STOP
      END
      SUBROUTINE CCMAIN
     1   (HBL1,HBL2,R2M2,R1M2,WT1,WT2,WALF,XALF,DNORM1,DNORM2,
     1    Y1,Y2,B,C,B1,CCC1,PLEG,PNORM,HAM1,EIG1,WORK1,
     1    BASS1,BASS2,HAM2,VECS1D,EIGS1D,EIG2,WORK2,
     1    HBAND,CINT,CINTP,XLMATR,XK1,XK2,R1M2T,R2M2T,R1,R2,
     1    HAM3,IV1,IV2,vecs2d,EIGS2D,EVAL,WORK3,NDIM2D,EIGS2,
     1    IV1L,IV2L,NDIM2L,VECS1L,VECS2L,VECS3L,PHI,PSI,EVALL)
C
C     SUBROUTINE CCMAIN IS THE 'REAL' MAIN PROGRAMME & CONTAINS     #006
C     THE CALLS TO THE VARIOUS SUBROUTINES WHICH SET & SOLVE THE
C     INTERMEDIATE AND THE FINAL HAMILTONIANS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
C
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /MASS/ XMASS(3),G1,G2
C
      DIMENSION DNORM1(0:NMAX1),NDIM2D(NDIMC),R1M2(NLIM1),
     9          DNORM2(0:NMAX2),R2M2(NLIM2),EIGS2(MAX2D,NDIMC),
     9          BASS1(0:NMAX1,NPNT1),BASS2(0:NMAX2,NPNT2),
     9          Y1(NPNT1),R1(NPNT1),WT1(NPNT1),
     9          Y2(NPNT2),R2(NPNT2),WT2(NPNT2),
     9          B(1),C(NPNT+1),PNORM(0:MAXLEG),
     9          HBL1(NLIM1),HBL2(NLIM2),B1(NALF),
     9          XALF(IDVR),WALF(IDVR),CCC1(NALF),XLMATR(IDVR,IDVR),
     9          PLEG(0:MAXLEG,IDVR),XK1(NPNT1,NPNT1),XK2(NPNT2,NPNT2),
     9          HAM1(NDIMA,NDIMA),HBAND(MAX2D,MAX3D),EIGS1D(max1dv),
     9          EIG1(NDIMA),WORK1(NDIMA),IV1(NDIMC,NDIMB),
     9          VECS1D(max1dv,NDIMA),HAM2(MAX2D,MAX2D),EIG2(MAX2D),
     9          WORK2(MAX2D),IV2(NDIMC),HAM3(MAX3D,MAX3D),
     9          CINT(NDIMA*NDIMB,MAX2D),CINTP(NDIMA*NDIMB,MAX2D),
     9          EIGS2D(MAX3D),EVAL(MAX3D),WORK3(MAX3D),
     9          R1M2T(NPNT1,NPNT1),R2M2T(NPNT2,NPNT2),EVALL(NEVAL),
     9          IV1L(NDIMC,NDIMB),IV2L(NDIMC),NDIM2L(NDIMC),
     9          VECS1L(MAX2D,NDIMA),VECS2L(MAX2D,MAX2D),VECS3L(MAX3D),
     9          PHI(MAX3D,NDIMA,NDIMB),PSI(IDVR,NPNT1,NPNT2),
     9          vecs2d(max2d,max3d)
      DATA X0/0.0D0/,X1/1.0d0/,X2/2.0D0/,X8/8.0D0/,XP5/0.5D0/,
     &     TOLER/1.0D-8/
C
C     ARRAYS ALWAYS IN CORE:
C     HAM3 VIBRATIONAL HAMILTONIAN MATRIX.
C
C     ARRAYS SET UP IN LAGPT:
C     R1M2,R2M2 STORE <R**(-2)> INTEGRATED OVER THE RADIAL BASIS;
C     B & C STORE COEFFICIENTS FOR RECURRENCE RELATIONSHIP;
C     BASS1,BASS2 STORES LAGUERRE POLYNOMIAL BASIS FUNCTIONS;
C     R1 & WT1 POINTS & WEIGHTS FOR NUMERICAL INTEGRATION ON R1;
C     R2 & WT2 POINTS & WEIGHTS FOR NUMERICAL INTEGRATION ON R2;
C     HBL1,HBL2 KINETIC ENERGY MATRICES OF NMAX * NMAX PROBLEMS.
C
C     ARRAYS SET UP IN SETFAC:
C     DNORM1,DNORM2  NORMALISATION CONSTANTS FOR LAGUERRE POLYNOMIALS.
C
C     ARRAYS USED IN DIAGONALISATION:
C     EVAL STORES EIGENVALUES;
C     WORK WORKSPACE FOR DIAGONALISER;
C     HAM3 STORES EIGENVECTORS OF SECULAR PROBLEM.
C
C     open the streams......(most not needed for symmetrised radau case)
      if (idia .gt. -2) then
         OPEN(UNIT=IEIGS1,FORM='UNFORMATTED')
         OPEN(UNIT=IEIGS2,FORM='UNFORMATTED')
         OPEN(UNIT=IVECS1,FORM='UNFORMATTED')
         OPEN(UNIT=IVECS2,FORM='UNFORMATTED')
         OPEN(UNIT=IVINT, FORM='UNFORMATTED')
         OPEN(UNIT=IBAND, FORM='UNFORMATTED')
      endif
      OPEN(UNIT=INTVEC,FORM='UNFORMATTED')
C
C     READ IN MASSES & BASIS SET PARAMETERS
      CALL SETCON(fixcos)
C.....and save them if necessary
      IF (ZVEC) WRITE(IOUT2) ZEMBED,ZMORS1,ZMORS2,ZTHETA,ZR2R1,XMASS,
     1                       G1,G2
C
C     SET UP BINOMIAL AND NORMALISATION ARRAYS
C     (B & C USED FOR WORK SPACE)
C
      IF (ZR2R1) THEN
         CALL SETFAC(BASS2,DNORM1,DNORM2,B,CC1,CC2)
      ELSE
         CALL SETFAC(BASS1,DNORM1,DNORM2,B,CC1,CC2)
      ENDIF
C
C     SET UP POINTS, WEIGHTs & basis FOR NUMERICAL INTEGRATION
C
      IF (NCOORD .EQ. 2) then
C     IN LESS THAN 3-D CASES FIX THE 3-D PARAMTERS
         IF (ZR2R1) THEN
            BASS1(0,1) = X1
            WT1(1) = X1
            R1(1) = RE1
            CALL LAGPT(2,Y2,R2,WT2,B,C,CC2,BASS2,DNORM2,npnt2,nmax2,
     1           zmors2,RE2,DISS2,WE2,BETA2,UR2,A2,IU2)
            HBL2(1) = X0
C           SETUP KINETIC ENERGY & INERTIA INTEGRALS OVER R1
            IF (ZMORS1) THEN
               FKE = BETA1 * BETA1 / (X8 * UR1)
               FKD = BETA1 / X2
               CALL KEINTS(HBL1,FKE,FKD,NLIM1,NMAX1,IU1)
            ELSE
               FKE = BETA1 / (X2 * UR1)
               CALL KEINT2(HBL1,FKE,R1M2,DNORM1,NLIM1,NMAX1,NPNT1,A1)
            ENDIF
         ELSE
            BASS2(0,1) = X1
            WT2(1) = X1
            R2(1) = RE2
            CALL LAGPT(1,Y1,R1,WT1,B,C,CC1,BASS1,DNORM1,npnt1,nmax1,
     1           zmors1,RE1,DISS1,WE1,BETA1,UR1,A1,IU1)
            HBL1(1) = X0
C           SETUP KINETIC ENERGY & INERTIA INTEGRALS OVER R2
            IF (ZMORS2) THEN
               FKE = BETA2 * BETA2 / (X8 * UR2)
               FKD = BETA2 / X2
               CALL KEINTS(HBL2,FKE,FKD,NLIM2,NMAX2,IU2)
            ELSE
               FKE = BETA2 / (X2 * UR2)
               CALL KEINT2(HBL2,FKE,R2M2,DNORM2,NLIM2,NMAX2,NPNT2,A2)
            ENDIF
         ENDIF
      ELSE
         CALL LAGPT(1,Y1,R1,WT1,B,C,CC1,BASS1,DNORM1,npnt1,nmax1,zmors1,
     1           RE1,DISS1,WE1,BETA1,UR1,A1,IU1)
         if (idia .gt. -2) then
            CALL LAGPT(2,Y2,R2,WT2,B,C,CC2,BASS2,DNORM2,npnt2,nmax2,
     1           zmors2,RE2,DISS2,WE2,BETA2,UR2,A2,IU2)
C           SETUP KINETIC ENERGY & INERTIA INTEGRALS OVER R2
            IF (ZMORS2) THEN
              FKE = BETA2 * BETA2 / (X8 * UR2)
              FKD = BETA2 / X2
              CALL KEINTS(HBL2,FKE,FKD,NLIM2,NMAX2,IU2)
            ELSE
              FKE = BETA2 / (X2 * UR2)
              CALL KEINT2(HBL2,FKE,R2M2,DNORM2,NLIM2,NMAX2,NPNT2,A2)
            ENDIF
         ENDIF
C        SETUP KINETIC ENERGY & INERTIA INTEGRALS OVER R1
         IF (ZMORS1) THEN
            FKE = BETA1 * BETA1 / (X8 * UR1)
            FKD = BETA1 / X2
            CALL KEINTS(HBL1,FKE,FKD,NLIM1,NMAX1,IU1)
         ELSE
            FKE = BETA1 / (X2 * UR1)
            CALL KEINT2(HBL1,FKE,R1M2,DNORM1,NLIM1,NMAX1,NPNT1,A1)
         ENDIF
      ENDIF
C
C     WRITE THE QUADRATURE POINTS TO DISK FOR ZVEC = .TRUE.
      IF (ZVEC) THEN
         WRITE (IOUT2) R1
         if (idia .gt. -2) WRITE (IOUT2) R2
      ENDIF
      IF (NCOORD .EQ. 3) THEN
      ENDIF
C
C     take square roots of the weights
      CALL RTWTS(WT1,NPNT1)
      if (idia .gt. -2) CALL RTWTS(WT2,NPNT2)
C
C     SET UP THE TRANSFORMED KINETIC ENERGY INTEGRALS,  T'(HBL) T
C                                                       ~  ~~~  ~
      CALL K1K2(XK1,HBL1,BASS1,WT1,NPNT1,NMAX1,NLIM1)
      if (idia .gt. -2) CALL K1K2(XK2,HBL2,BASS2,WT2,NPNT2,NMAX2,NLIM2)
C
C     ...... and the inertia integrals for spherical oscillators
      IF (.NOT.ZMORS1 .AND. .NOT.ZTWOD)
     1 CALL K1K2(R1M2T,R1M2,BASS1,WT1,NPNT1,NMAX1,NLIM1)
      IF (.NOT.ZMORS2 .AND. .NOT.ZTWOD .and. idia .gt. -2)
     1 CALL K1K2(R2M2T,R2M2,BASS2,WT2,NPNT2,NMAX2,NLIM2)
c     for J > 0, store r**(-2) term for rotlev3
      if (ztran) then
         if (zrot) then
            if (zembed) then
               if (zquad2) then
                  write(iwave) (XP5/(R2(I)*R2(I)*UR2),i=1,npnt2)
               else
                  call outrow(r2m2t,npnt2*npnt2,iwave)
               endif
            else
               write(iwave) (XP5/(R1(I)*R1(I)*UR1),i=1,npnt1)
            endif
         else
            jdia=max(1,idia)
            jstart=kmin
            IF (MOD(JSTART,JDIA) .NE. IPAR) JSTART=JSTART+1
            NANG=(MAXLEG-JSTART)/JDIA+1
            mbass=idvr*npnt1*npnt2
            if (idia .eq. -2) mbass=idvr*max2d
            write(iwave) mbass,jstart,nang,mbass
         endif
         WRITE(iwave) R1
         if (idia .gt. -2) WRITE(iwave) R2
      endif
C
C     Some of the J>0 stuff to get the loop over K correct
      IF (ZROT) THEN
        KD = 1-MIN(KMIN,1)
        KU = JROT
C       If looping over sym & anti-sym, do one extra k=1 block
        IF (KMIN .EQ. 2) KU=KU+1
      ELSE
        KD = KMIN
        KU = KMIN
      ENDIF
      KKZ12 = 0
      KKZ0  = 0
c
      if (idia .eq. -2) then
         max2d1=max2d
         max3d1=max3d
         mx1dv1=max1dv
         mx1dv2=min(max1dv,max2d2*nalf)
      endif
C
C     -------------  start rotational loop here  -------------
C
      DO 40 kk=KD,KU
      IF (KK .LE. JROT) THEN
         KZ=KK
      ELSE
         KZ=1
         KMIN=0
         IPAR=MOD(IPAR+JROT,2)
      ENDIF
c
c     first rewind the scratch files for a calculation with J>0
c     and, if needed, reposition IOUT2 after set up data.
      IF (kk .GT. KD) THEN
        if (idia .gt. -2) then
           REWIND IEIGS1
           REWIND IEIGS2
           REWIND IVECS1
           REWIND IVECS2
           REWIND IVINT
           REWIND IBAND
        endif
        REWIND INTVEC
        IF (ZVEC) THEN
          REWIND IOUT1
          REWIND IOUT2
          do 45 ii=1,4
   45     read(iout2)
        ENDIF
      ENDIF
C
      REALKZ = DBLE(KZ)
C     TSWALF IS THE EXACT SUM OF WEIGHTS FOR GAUSS-JACOBI INTEGRATION
      TSWALF= X2**(KZ+KZ+1)/DBLE(KZ+1)
      DO 30 IA=1,KZ
   30 TSWALF=TSWALF*DBLE(IA)/DBLE(KZ+IA+1)
C
      IF (ZLADD .OR. KZ .EQ. KD) THEN
         NIDVR = IDVR
         NANG  = NALF
         NANG2 = NALF2
         LINCR = KZ
      ELSE
         LINCR = 0
         IF (IDIA .NE. 2) THEN
            NIDVR = IDVR - KZ
            NANG  = NALF - KZ
            NANG2 = (NANG+1)/2
         ELSE IF (IPAR .EQ. 0) THEN
            IF(MOD(KZ,2).EQ.1) KKZ0 = KKZ0 + 2
            NIDVR = IDVR - KKZ0/2
            NANG  = NALF - KKZ0
            NANG2 = (NANG+1)/2
         ELSE
            IF(MOD(KZ,2).EQ.0 .AND. KZ.GT.0) KKZ12 = KKZ12 + 2
            NIDVR = IDVR - KKZ12/2
            NANG  = NALF - KKZ12
            NANG2 = (NANG+1)/2
         ENDIF
      ENDIF
      if (.not. zladd) then
         if (ztheta) then
            npnta=nidvr
         else
            npntc=nidvr
         endif
      endif
C
      IF (ZTWOD) THEN
         XALF(1) = FIXCOS
         GOTO 333
      ENDIF
C
      CALL JACOBI(NANG,NANG2,XALF,WALF,REALKZ,REALKZ,B1,
     1            CCC1,CSWALF,TSWALF)
C
      WRITE(6,1000) NANG,KZ,(XALF(II),WALF(II),II=1,NANG2)
 1000 FORMAT(//I8,' point Gauss-assocaited Legendre integration with',
     1       ' K =',I3//5X,'INTEGRATION POINTS',11X,'WEIGHTS',
     2        /(F23.15,D25.12))
      WRITE(6,1010) CSWALF,TSWALF
 1010 FORMAT(/5X,'COMPUTED SUM OF WEIGHTS',D26.15,
     1       /5X,'EXACT    SUM OF WEIGHTS',D26.15//)
      IF (ABS((CSWALF-TSWALF)/TSWALF) .GT. TOLER) THEN
         WRITE(6,910)
  910    FORMAT(//5X,'POINTS & WEIGHTS IN ERROR, ADJUST ALGORITHM'//)
         STOP
      ENDIF
      CALL ALLPTS(XALF,WALF,NANG,NANG2)
C
      IF (ZVEC)  WRITE(IOUT2) XALF
C
C     SET UP LEGENDRE POLYNOMIALS FOR THE TRANSFORMATION MATRIX
c     save them for rotlev3 or rotlev3b if required
      CALL ASLEG(PLEG,PNORM,MAXLEG,XALF,NIDVR,KZ,LINCR)
      if (ztran) then
        write(iwave) xalf
        write(iwave) kz,maxleg,nidvr,lincr
        write(iwave) ((pleg(i,j)*walf(j),i=0,maxleg),j=1,nidvr)
      endif
C     BUILD THE TRANSFORMED ANGULAR MOMENTUM MATRIX XLMATR;
      IPAR0=0
      IF (IDIA .EQ. 2 .AND. IPAR .EQ. 1) IPAR0=1
      CALL LMATRX(XLMATR,PLEG,WALF,KZ,IPAR0,NIDVR,LINCR)
C
  333 CONTINUE
c
c     For AB2 molecules in Radau coordinates, use separate main
c     driving routine
c
      if (idia .eq. -2) then
         IF (zrot) then
            if (KZ .GT. KD) IPAR=MOD(IPAR+1,2)
            if (ipar .gt. 0) then
               max2d=max2d2
               max3d=max3d2
               max1dv=mx1dv2
            else
               max2d=max2d1
               max3d=max3d1
               max1dv=mx1dv1
            endif
         endif
         call nfmain(xk1,xlmatr,r1,xalf,ham1,eig1,iv1,vecs1d,
     &               eigs1d,work1,ham2,eig2,iv2,vecs2d,eigs2d,ndim2d,
     &               work2,eigs2,ham3,eval,work3,psi,kz)
      else
         call jhmain(xk1,xk2,xlmatr,r1,r2,xalf,ham1,eig1,iv1,vecs1d,
     &               eigs1d,work1,ham2,eig2,iv2,eigs2d,ndim2d,
     &               work2,eigs2,ham3,eval,work3,cint,cintp,hband,
     &               r1m2t,r2m2t,iv1l,iv2l,ndim2l,vecs1l,vecs2l,
     &               vecs3l,phi,psi,evall,kz)
      endif
C
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETCON(fixcos)
C
C     READ IN MASSES & SET CONSTANTS FOR RADIAL BASIS SETS          #007
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
C     SAVE MASSES & G IN CASE THEY ARE NEEDED IN THE POTENTIAL ROUTINE
      COMMON /MASS/ XMASS(3),G1,G2
C     AMTOAU CONVERTS AMU (PROTON MASSES) TO AU (ELECTRON MASSES).
c     DATA AMTOAU/1.8228883D03/
      DATA AMTOAU/1.82284475D03/
      DATA X0,XP5,X1,X4/0.0D0,0.5D0,1.0D0,4.0D0/
C
C     READ COS(THETA) FOR FIXED ANGLE 2-D CALCULATION
C
      READ(5,5) FIXCOS
      IF (ZTWOD) WRITE(6,1088) FIXCOS
 1088 FORMAT(//10X,'TWO-D FIXED ANGLE VIBRATIONAL PROBLEM WITH'
     1       //10X,'***** FIXED VALUE OF COS(THETA) =',F6.2,' *****'/)
C
C     READ MASSES OF THE ATOMS IN ATOMIC MASS UNITS
C
      READ(5,5)     XMASS
    5 FORMAT(3F20.0)
C     READ cut off energies
C     READ PARAMETERS DEFINING ENERGY CUT OFFS FOR EACH BLOCK
      READ(5,5) EMAX1,EMAX2
      IF (.NOT. ZALL) THEN
         IF (ZCUT) THEN
            if (idia .gt. -2) WRITE(6,990) EMAX1,EMAX2
  990       FORMAT(//5X,'CUT-OFF ENERGIES IN WAVENUMBERS:',2D16.8/)
            if (idia .eq. -2) WRITE(6,991)       EMAX2
  991       FORMAT(//5X,'Final cut-off energy in wavenumbers:',2D16.8/)
         ELSE
            if (idia .gt. -2) WRITE(6,992) EMAX1
  992       FORMAT(//5X,'FIRST CUT-OFF ENERGY IN WAVENUMBERS:',1D16.8/)
         ENDIF
      ENDIF
C     SET DEFAULT VALUE OF G1 AND G2
      IF (IDIA .GE. 1) THEN
C        SCATTERING COORDINATES
         G1 = XMASS(2) / (XMASS(2) + XMASS(3))
         G2 = X0
      ELSE
C        RADAU COORDINATES
         A = SQRT(XMASS(3) / (XMASS(1)+XMASS(2)+XMASS(3)))
         B = XMASS(2) / (XMASS(1)+XMASS(2))
         G1 = X1 - A / (A+B-A*B)
         G2 = X1 - A / (X1-B+A*B)
      ENDIF
C     NCOORD = 3: READ PARAMETERS FOR R1 RADIAL BASIS (SEE BELOW)
C     NCOORD = 2: READ FIXED R1 BONDLENGTH RE1. DISS1 & WE1 DUMMY
      READ(5,5)     RE1,DISS1,WE1
C     READ PARAMETERS FOR R2 RADIAL BASIS FUNCTION,
C     FOR MORSE OSCILLATOR FUNCTIONS USE THE FOLLOWING:
C     RE2: EQUILIBRIUM BONDLENGTH OF R2 COORDINATE (IN BOHR)
C     DISS2: DISSOCIATION ENERGY OF THE R2 COORDINATE (IN HARTREE)
C     WE2: FUNDAMENTAL STRETCHING VIBRATION OF R2 (IN HARTREE)
C     FOR SPHERICAL OSCILLATOR FUNCTIONS USE THE FOLLOWING:
C     RE2 : DUMMY
C     DISS2: ORDER OF LAGUERRE POLYNOMIALS USED (DIMENSIONLESS)
C     WE2: FUNDAMENTAL STRETCHING VIBRATION OF R2 (IN HARTREE)
C     ALL ARE TREATED AS VARIATIONALLY OPTIMISABLE PARAMETERS.
      READ(5,5)     RE2,DISS2,WE2
      WRITE(6,1000) XMASS
 1000 FORMAT(//5X,'NUCLEAR MASSES IN AMU:',3F12.6,/)
C     COMPUTE THE EFFECTIVE MOMENTS OF INERTIA
      UR1 = AMTOAU/(G2*G2/XMASS(1)+X1/XMASS(2)+(X1-G2)**2/XMASS(3))
      UR2 = AMTOAU/(X1/XMASS(1)+G1*G1/XMASS(2)+(X1-G1)**2/XMASS(3))
      IF (NCOORD .EQ. 3) GOTO 20
      IF (ZR2R1) THEN
         WRITE(6,1010) RE1,UR1
 1010 FORMAT(/5X,'R1 FIXED BONDLENGTH =',F8.4,' BOHR',
     1            ' & REDUCED MASS =',D16.7,' A.U.'/)
      ELSE
         WRITE(6,1011) RE2,UR2
 1011 FORMAT(/5X,'R2 FIXED BONDLENGTH =',F8.4,' BOHR',
     1            ' & REDUCED MASS =',D16.7,' A.U.'/)
      ENDIF
      IF (ZR2R1) GOTO 30
   20 CONTINUE
      IF (ZMORS1) THEN
          WRITE(6,1020) 1,RE1,DISS1,WE1
 1020 FORMAT(/5X,'MORSE FUNCTION PARAMETERS FOR R',I1,' BASIS',
     1       /5X,'R EQUILIBRIUM =',F8.4,' BOHR, DISSOCIATION ENERGY',
     2    D15.7,' HARTREE &  VIBRATIONAL FREQUENCY =',D15.7,' HARTREE')
         BETA1 = WE1 * SQRT(XP5*UR1/DISS1)
         A1 = X4 * DISS1 / WE1
         IU1 = INT(A1+XP5)
         WRITE(6,1030) UR1,BETA1,A1,IU1
 1030 FORMAT(/5X,'CONSTANTS USED TO CONSTRUCT MORSE OSCILLATORS:',
     1       /5X,'REDUCED MASS =',D16.7,' A.U., BETA =',F8.4,
     2            ' (1/BOHR), A =',D16.7,' AND U =',I5)
      ELSE
          A1=DISS1
          BETA1 = SQRT(WE1 * UR1)
          WRITE(6,1039) A1,WE1,UR1,BETA1
 1039 FORMAT(/5X,'SPHERICAL OSCILLATOR PARAMETERS FOR R1 BASIS:',
     1       /5X,'ALPHA =',F10.5,
     2           ' &  VIBRATIONAL FREQUENCY =',D15.7,' HARTREE',
     3      //5X,'CONSTANTS USED TO CONSTRUCT SPHERICAL OSCILLATORS:',
     4       /5X,'REDUCED MASS =',D16.7,' A.U., BETA =',F12.6,
     5            ' BOHR**-2')
      ENDIF
      IF (NCOORD .EQ. 2 .AND. .NOT. ZR2R1) GOTO 40
      if (idia .eq. -2) goto 40
   30 CONTINUE
      IF (ZMORS2) THEN
         WRITE(6,1020) 2,RE2,DISS2,WE2
         BETA2 = WE2 * SQRT(XP5*UR2/DISS2)
         A2 = X4 * DISS2 / WE2
         IU2 = INT(A2+XP5)
         WRITE(6,1030) UR2,BETA2,A2,IU2
      ELSE
         A2=DISS2
         BETA2 = SQRT(WE2 * UR2)
         WRITE(6,1040) A2,WE2,UR2,BETA2
 1040 FORMAT(/5X,'SPHERICAL OSCILLATOR PARAMETERS FOR R2 BASIS:',
     1       /5X,'ALPHA =',F10.5,
     2           ' &  VIBRATIONAL FREQUENCY =',D15.7,' HARTREE',
     3      //5X,'CONSTANTS USED TO CONSTRUCT SPHERICAL OSCILLATORS:',
     4       /5X,'REDUCED MASS =',D16.7,' A.U., BETA =',F12.6,
     5            ' BOHR**-2')
      ENDIF
   40 CONTINUE
      IF (ZTRAN) THEN
         if (zrot) then
            if (zembed) then
               if (zquad2) then
                  nlim = npnt2
               else
                  nlim = npnt2*npnt2
               endif
            else
               nlim = npnt1
            endif
         else
            nlim = 0
         endif
         zncor=.not.zrot
         if (idia .gt. -2) then
            WRITE(IWAVE) IDIA,IPAR,IDVR,NPNT1,NPNT2,JROT,KMIN,NEVAL,NLIM
            WRITE(IWAVE) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,zncor,ZQUAD2
            WRITE(IWAVE) RE1,DISS1,WE1,RE2,DISS2,WE2
         else
            WRITE(IWAVE) IDIA,IPAR,IDVR,NPNT1,NPNT1,JROT,KMIN,NEVAL,NLIM
            WRITE(IWAVE) ZEMBED,ZMORS1,ZMORS1,XMASS,G1,G2,zncor,ZQUAD2
            WRITE(IWAVE) RE1,DISS1,WE1,RE1,DISS1,WE1
         endif
      ENDIF
      RETURN
      END
      SUBROUTINE SETFAC(BIN,DNORM1,DNORM2,FACT,CC1,CC2)
C
C     SETFAC INITIALISES BINOMIAL ARRAY:                            #021
C       BINOM(I+1,J+1) = I! / (J! * (I-J)!)
C     AND PSEUDO-NORMALISATION ARRAY:
C       DNORM(M) = SQRT((M-1)! * BINOM(NPNT+IU,NPNT-M))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
      DIMENSION DNORM1(0:NMAX1),DNORM2(0:NMAX2),FACT(*),BIN(*)
      DATA XP5/0.5D0/,X1/1.0D0/
      FACT(1) = X1
      COUNT = X1
      DO 10 I=1,NPNT
      FACT(I+1) = COUNT * FACT(I)
      COUNT = COUNT + X1
   10 CONTINUE
      IF (NCOORD .EQ. 2) THEN
        IF (.NOT. ZR2R1) THEN
          IF (ZMORS1) THEN
            ALF=DBLE(IU1)
          ELSE
            ALF=A1+XP5
          ENDIF
          CALL NORMS(DNORM1,BIN,FACT,CC1,ALF,NPNT1,NMAX1)
        ELSE
          IF (ZMORS2) THEN
            ALF=DBLE(IU2)
          ELSE
            ALF=A2+XP5
          ENDIF
          CALL NORMS(DNORM2,BIN,FACT,CC2,ALF,NPNT2,NMAX2)
        ENDIF
      ELSE
        IF (ZMORS1) THEN
          ALF=DBLE(IU1)
        ELSE
          ALF=A1+XP5
        ENDIF
        CALL NORMS(DNORM1,BIN,FACT,CC1,ALF,NPNT1,NMAX1)
        if (idia .gt. -2) then
           IF (ZMORS2) THEN
             ALF=DBLE(IU2)
           ELSE
             ALF=A2+XP5
           ENDIF
           CALL NORMS(DNORM2,BIN,FACT,CC2,ALF,NPNT2,NMAX2)
        ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE NORMS(DNORM,BIN,FACT,CC,ALF,NPNT,NMAX)
C     SET UP FACTORS FOR NORMALISING THE RADIAL BASIS FUNCTIONS     #022
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
C
      DIMENSION DNORM(0:NMAX),BIN(*),FACT(*)
      DATA X1/1.0D0/
      COUNT = DBLE(NPNT) + ALF
C     CC IS EXACT SUM OF WEIGHTS FOR NPNT GAUSS-LAGUERRE INTEGRATION
      CC = FACT(NPNT) / COUNT
C     NORMALISATION ARRAY FOR L(I,ALF): FIRST SET UP BINOMIALS
      NPT1 = NPNT + 1
      BIN(NPT1) = X1
      DO 20 I=1,NPNT
      N = NPT1 - I
      BIN(N) = BIN(N+1) * COUNT / DBLE(I)
      COUNT = COUNT - X1
   20 CONTINUE
      NPT1 = NPT1 + 1
      DO 30 I=0,NMAX
      DNORM(I) = SQRT(BIN(I+1)*FACT(I+1)*FACT(NPT1-I-1))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE LAGPT(ir,Y,R,WT,B,C,CC,BASS,DNORM,npnt,nmax,zmorse,
     1                 RE,DISS,WE,BETA,UR,A,IU)
C
C     SUBROUTINE LAGPT GETS INTEGRATION POINTS AND WEIGHTS FOR      #015
C     NPNT GAUSS LAGUERRE INTEGRATION AND SETS UP BASIS
C     FUNCTIONS AT THE INTEGRATION POINTS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Y), LOGICAL (Z)
      DIMENSION B(NPNT+1),C(NPNT+1),Y(NPNT),R(NPNT),WT(NPNT),
     1          DNORM(0:NMAX),BASS(0:NMAX,NPNT)
      DATA X0,XP5,X1,X2/0.0D0,0.5D0,1.0D0,2.0D0/,TOLER/1.0D-8/
C
      IF (ZMORSE) ALF=DBLE(IU)
      IF (.NOT. ZMORSE) ALF = A + XP5
      ALFM1=ALF-X1
C     SET UP INTEGRATION POINTS AND WEIGHTS.
      CALL LAGUER(NPNT,Y,WT,ALF,B,C,CSX,CSA,TSX,TSA,CC)
      TSA = X1 / (DNORM(0) * DNORM(0))
      WRITE(6,1000) NPNT,ir
 1000 FORMAT(/,I8,' POINT GAUSS-LAGUERRE INTEGRATION',
     1       /,5X,'INTEGRATION POINTS',11X,'WEIGHTS',9X,
     1            'CORRESPONDING R',I1,/)
      DO 60 I=1,NPNT
      IF (ZMORSE) THEN
C         CALCULATE POTENTIAL AT R = RE+BETA(**-1)*LN(A/Y)
          R(I) = RE + LOG(A/Y(I)) / BETA
      ELSE
C         CALCULATE POTENTIAL AT R = SQRT(Y/BETA)
          R(I) = SQRT(Y(I)/BETA)
      ENDIF
      WRITE(6,1010) Y(I),WT(I),R(I)
 1010 FORMAT (F23.15,D25.12,F13.5)
      IF (R(I) .LT. X0) WRITE(6,1015) I
 1015 FORMAT(5X,'***** WARNING: FOR INTEGRATION POINT',I3,
     1       ', R LESS THAN ZERO *****')
C
C     CALCULATE UNNORMALISED LAGUERRE POLYNOMIALS AT Y
C
C     POLYNOMIAL OF ORDER 0
      BASS(0,I) = X1
      IF (NMAX .LT. 1) GOTO 70
C     POLYNOMIAL OF ORDER 1
      AMX = ALF + X1 - Y(I)
      BASS(1,I) = AMX
C     USE RECURRENCE RELATIONSHIPS FOR POLYNOMIALS OF ORDER > 2
C     N * L(N,ALF) = (2*N+ALF-1-X)*L(N-1,ALF) - (N+ALF-1)*L(N-2,ALF)
      EN = X1
      DO 80 N=2,NMAX
      EN = EN + X1
      AMX = AMX + X2
      BASS(N,I) = (AMX * BASS(N-1,I) - (ALFM1+EN) * BASS(N-2,I)) / EN
   80 CONTINUE
   70 CONTINUE
C
      DO 90 N2=0,NMAX
C     NORMALISE POLYNOMIALS
      BASS(N2,I) = BASS(N2,I) * DNORM(N2)
   90 CONTINUE
   60 CONTINUE
C     CHECK THAT THE CORRECT POINTS & WEIGHTS HAVE BEEN GENERATED
      WRITE(6,1020) CSX,CSA,TSX,TSA
 1020 FORMAT(/5X,'COMPUTED SUM OF POINTS',D26.15,' & WEIGHTS',D26.15,
     1       /5X,'EXACT    SUM OF POINTS',D26.15,' & WEIGHTS',D26.15)
      IF (ABS((CSX-TSX)/TSX) .GT. TOLER) GOTO 900
      IF (ABS((CSA-TSA)/TSA) .GT. TOLER) GOTO 900
      RETURN
  900 WRITE(6,910)
  910 FORMAT(//5X,'POINTS & WEIGHTS IN ERROR, ADJUST ALGORITHM',//)
      STOP
      END
      SUBROUTINE LAGUER(NN,X,A,ALF,B,C,CSX,CSA,TSX,TSA,CC)
C
C     CALCULATES POINTS & WEIGHTS FOR GAUSS-LAGUERRE INTEGRATION    #016
C     SEE:
C     "GAUSSIAN QUADRATURE FORMULAS" BY A.H.STROUD & D.SECREST
C      1966, PRENTICE-HALL, P.32.
C     **** VERSION TO AVOID OVERFLOWS (J.T. 25/11/81) ****
C     CALCULATES WEIGHTS DIVIDED BY GAMMA(NN+ALF+1)
C     THIS IS AN INITIALSATION ENTRY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      DIMENSION X(NN),A(NN),B(NN+1),C(NN+1)
      DATA EPS/1.0D-12/,X1/1.0D0/
      CSX=0.0D0
      CSA=0.0D0
      FA=ALF+1.0D0
C     CC = N!                      DENOMINATOR FOR PSEUDO-WEIGHTS: A
C     B(N) = (ALF + 2N -1)             B & C FOR RECURRENCE RELATION
C     C(N) = (N - 1) * ( ALF + N - 1)
      B(1)=FA
      C(1)=0.0D0
      FN=1.0D0
      DO 1 J=2,NN
      FA=FA+2.0D0
      B(J)=FA
      C(J)=FN*(ALF+FN)
    1 FN=FN+1.0D0
      TSX=FN*(ALF+FN)
      XT1=0.0D0
C     FORMULAS FOR INITIAL POINT & STEP CHOSEN BECAUSE THEY WORK!
      XT=(1.0D0+ALF)*(2.0D0+ALF)/(1.0D0+3.0D0*FN+2.0D0*ALF)
      STEP=3.0D0*(1.0D0+ALF)/(1.0D0+3.0D0*FN+ALF)
      CALL LGRECR(PT,DPN,PN1,XT,NN,ALF,B,C)
 
      DO 7 I=1,NN
      IF (I-2) 2,2,4
C     SMALLEST TWO ZEROS: FOUND BY "BRUTE FORCE" SEARCH
    2 XT2 = XT + STEP
      CALL LGRECR(PT2,DPN,PN1,XT2,NN,ALF,B,C)
      IF (SIGN(X1,PT)*SIGN(X1,PT2) .GT. 0.0D0) GOTO 5
      PT = PT2
      XT = 0.5D0 * (XT + XT2)
      GOTO 6
    5 PT = PT2
      XT = XT2
      GOTO 2
C     ALL OTHER ZEROS: FOUND USING FORMULA OF STROUD & SECREST
    4 FI = DBLE(I-2)
      R1 = (1.0D0+2.55D0*FI)/(1.9D0*FI)
      R2 = 1.26D0*FI*ALF/(1.0D0+3.5D0*FI)
      RATIO = (R1+R2)/(1.0D0+0.3D0*ALF)
      XT = XT + RATIO*(XT-XT2)
C
    6 CALL LGROOT(XT,NN,ALF,DPN,PN1,B,C,EPS)
      XT2=XT1
      XT1=XT
      X(I) = XT
      A(I) = CC/DPN/PN1
      CSX = CSX + XT
      CSA = CSA + A(I)
    7 CONTINUE
      RETURN
      END
      SUBROUTINE LGROOT(X,NN,ALF,DPN,PN1,B,C,EPS)
C
C     IMPROVES THE APPROXIMATE ROOT X; IN ADDITION OBTAINS          #017
C          DPN = DERIVATIVE OF P(N) AT X
C          PN1 = VALUE OF P(N-1) AT X
C     THIS ROUTINE IS DUE TO STROUD & SECREST (SEE SUBROUTINE LAGUER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      DIMENSION B(NN+1),C(NN+1)
      DATA ITMAX/10/
      ITER=0
    1 ITER=ITER+1
      CALL LGRECR(P,DPN,PN1,X,NN,ALF,B,C)
      D = P/DPN
      X = X-D
      IF (ABS(D/X) .LE. EPS) RETURN
      IF (ITER - ITMAX) 1,2,2
    2 WRITE(6,100) ITER,D,X
  100 FORMAT(5X,'WARNING: NOCONVERGENCE AFTER',I4,' ITERATIONS',
     1     /,5X,'CURRENT DIFFERENCE',D26.15,' & ROOT',D26.15)
      RETURN
      END
      SUBROUTINE LGRECR(PN,DPN,PN1,X,NN,ALF,B,C)
C
C     USES RECURRENCE RELATIONS TO SET UP POLYNOMIALS               #018
C     THIS ROUTINE IS DUE TO STROUD & SECREST (SEE SUBROUTINE LAGUER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      DIMENSION B(NN+1),C(NN+1)
      P1 = 1.0D0
      P = X - ALF - 1.0D0
      DP1 = 0.0D0
      DP = 1.0D0
      DO 1 J=2,NN
      Q  = (X-B(J))* P-C(J)* P1
      DQ = (X-B(J))*DP-C(J)*DP1 + P
      P1 = P
      P  = Q
      DP1= DP
    1 DP = DQ
      PN = P
      DPN= DP
      PN1= P1
      RETURN
      END
      SUBROUTINE KEINTS(HBL,FKE,FKD,NLIM,NMAX,IU)
C
C     KEINTS CALCULATES ANALYTIC KINETIC ENERGY INTEGRALS OVER R    #012
C     FOR MORSE OSCILLATOR-LIKE FUNCTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
C
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
C
      DIMENSION HBL(NLIM)
C
      DATA X0/0.0D0/
      INDEX = 0
      DO 10 N2=1,NMAX+1
      DO 10 N1=1,N2
      INDEX=INDEX+1
      IF (N1 .EQ. N2) THEN
C        SPECIAL CASE:  N1 = N2
         HBL(INDEX) = FKE * DBLE(2*(N2-1)*(N2+IU)+IU+1)
      ELSE IF (N1+2 .EQ. N2) THEN
C        SPECIAL CASE:  N2 = N1 + 2
         HBL(INDEX) = - FKE *
     1                SQRT(DBLE((IU+N1+1)*(IU+N1))*DBLE((N1+1)*N1))
      ELSE
C        N1+1 = N2 OR N2 > N1 + 2  ALL MATRIX ELEMENTS ARE ZERO
         HBL(INDEX) = X0
      ENDIF
   10 CONTINUE
      IF (.NOT. ZPRAD) RETURN
C     WRITE KINETIC ENERGY INTEGRALS IF REQUESTED
      WRITE(6,500)
  500 FORMAT(//,5X,'RADIAL KINETIC ENERGY MATRIX CALCULATED',
     1             ' ANALYTICALLY',/)
      CALL SYMOUT(HBL,NMAX+1)
      RETURN
      END
      SUBROUTINE KEINT2(HBL,FKE,RM2,DNORM,NLIM,NMAX,NPNT,ALF)
C
C     KEINT2 CALCULATES ANALYTIC KINETIC ENERGY INTEGRALS OVER R2   #013
C     AND MOMENT OF INTERTIA INTEGRAL FOR SPHERICAL OSCILLATOR FUNCTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
C
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
C
      DIMENSION HBL(NLIM),RM2(NLIM),DNORM(*)
C
      DATA X0/0.0D0/,XP5/0.5D0/,X1/1.0D0/
      GAM = FKE / (ALF + XP5)
      DO 10 N1=1,NPNT
   10 GAM = GAM /(DBLE(N1)+ALF+XP5)
      FACT = X1
      FN = X0
      SUM = GAM
      DO 20 N1 = 1,NMAX+1
      INDEX = (N1 * (N1+1)) / 2
      DO 30 N2 = N1,NMAX+1
      RM2(INDEX) = DNORM(N1) * DNORM(N2) * SUM
      INDEX = INDEX + N2
   30 CONTINUE
      GAM = (FN+ALF+XP5) * GAM
      FN = FN + X1
      FACT = FN * FACT
      SUM = SUM + GAM / FACT
   20 CONTINUE
      FACT = ALF * (ALF + X1)
      DO 40 INDEX=1,NLIM
      HBL(INDEX) = - FACT * RM2(INDEX)
   40 CONTINUE
      INDEX = 0
      FN = - X1
      DO 50 N2=1,NMAX+1
      FN = FN + X1
      INDEX=INDEX+N2
C     SPECIAL CASE:  N1 = N2
      HBL(INDEX) = HBL(INDEX) + FKE * (FN+FN+ALF+1.5D0)
C     SPECIAL CASE:  N2 = N1 + 1
      IF (N2 .GT. 1)
     1    HBL(INDEX-1) = HBL(INDEX-1) + FKE * SQRT(FN*(FN+ALF+XP5))
   50 CONTINUE
      IF (.NOT. ZPRAD) RETURN
C     WRITE KINETIC ENERGY & INERTIA INTEGRALS IF REQUESTED
      WRITE(6,500)
  500 FORMAT(//5X,'RADIAL KINETIC ENERGY MATRIX CALCULATED',
     1             ' ANALYTICALLY'/)
      CALL SYMOUT(HBL,NMAX+1)
      WRITE(6,510)
  510 FORMAT(//5X,'MOMENT OF INERTIA MATRIX CALCULATED ANALYTICALLY'/)
      CALL SYMOUT(RM2,NMAX+1)
      RETURN
      END
      SUBROUTINE RTWTS(WT,NPNT)
C
C     TAKE THE ROOTS OF THE RADIAL WEIGHTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      DIMENSION WT(NPNT)
      DO 10 I=1,NPNT
         WT(I) = SQRT(WT(I))
   10 CONTINUE
      RETURN
      END
      SUBROUTINE K1K2(XK,HBL,BASS,WT,NPNT,NMAX,NLIM)
C
C     SET UP THE TRANSFORMED KINETIC ENERGY INTEGRALS,  T'(HBL) T
C                                                       ~  ~~~  ~
C     (NOTE THAT THE RADIAL BASIS FUNCTIONS ARE ALREADY NORMALISED)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
C
      DIMENSION XK(NPNT,NPNT),HBL(NLIM),BASS(0:NMAX,NPNT),WT(NPNT)
C
C     IND(I,J) = MAX(I,J) * (MAX(I,J)-1)/2 + MIN(I,J)
      IND(I,J) = MAX(I,J) * (MAX(I,J)+1)/2 + MIN(I,J) + 1
C
      DO 5 I=1,NPNT
      DO 5 IP=1,I
          XK(I,IP) = 0.0D0
    5 CONTINUE
C
      DO 10 K=1,NPNT
      DO 10 KP=1,K
      WTKKP = WT(K)*WT(KP)
      DO 20 M=0,NMAX
         T = BASS(M,K) * WTKKP
      DO 20 MP=0,NMAX
         IN = IND(M,MP)
         XK(K,KP) = XK(K,KP) + (HBL(IN) * T * BASS(MP,KP))
   20 CONTINUE
      xk(kp,k)=xk(k,kp)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE JACOBI(NN,NN2,X,A,ALF,BTA,B,C,CSA,TSA)
C
C     CALCULATES ZEROS X(I) OF THE NN'TH ORDER JACOBI POLYNOMIAL
C     PN(ALF,BTA) FOR THE SEGMENT (-1,1) & AND CORRESPONDING WEIGHTS
C     FOR GAUSS-JACOBI INTEGRATION. THIS ROUTINE, AND THOSE WHICH
C     FOLLOW, ARE DUE TO A.H.STROUD AND D. SECREST, "GAUSSIAN
C     INTEGRATION FORMULAS", 1966, PRENTICE HALL, PAGE 29.
C     NOTE THAT FOR OUR PURPOSES, ALF= BTA= NU.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(NN),A(NN),B(NN),C(NN)
      DATA X0/0.0D0/,X1/1.0D0/,X2/2.0D0/,EPS/1.0D-12/
      X3= X1 + X2
      X4= X2 + X2
      X8= X4 + X4
      FN= FLOAT(NN)
      CSA= X0
      C(1) = x0
      B(1) = x0
      DO 21 I=2,NN
      XI= FLOAT(I)
      B(I) = x0
      C(I)= X4*(XI-X1)*(ALF+XI-X1)*(BTA+XI-X1)*(ALF+BTA+XI-X1)/
     1       ((ALF+BTA+X2*XI-X1)*(ALF+BTA+X2*XI-X2)*
     2        (ALF+BTA+X2*XI-X2)*(ALF+BTA+X2*XI-X3))
21    CONTINUE
      CC=TSA
      DO 1 I=2,NN
      CC= CC*C(I)
 1    CONTINUE
      DO 12 I=1,NN2
      IF (I .EQ. 1) THEN
C        LARGEST ZERO
         AN= ALF/FN
         BN= BTA/FN
         R1= (X1 + ALF)*(2.78D0/(X4 + FN*FN) +0.768*AN/FN)
         R2= X1 + 1.48D0*AN + 0.96D0*BN + 0.452*AN*AN + 0.83D0*AN*BN
         XT= X1 - R1/R2
      ELSE IF (I .EQ. 2) THEN
C        SECOND ZERO
         R1= (4.1D0 + ALF)/((X1 + ALF)*(X1 + 0.156*ALF))
         R2= X1 + 0.06D0*(FN - X8)*(X1 + 0.12D0*ALF)/FN
         R3= X1 + 0.012*BTA*(X1 + ABS(ALF)/X4)/FN
         RATIO= R1*R2*R3
         XT= XT - RATIO*(X1 - XT)
      ELSE IF (I .EQ. 3) THEN
C        THIRD ZERO
         R1= (1.67D0 + 0.28D0*ALF)/(X1 + 0.37D0*ALF)
         R2= X1 + 0.22D0*(FN - X8)/FN
         R3= X1 + X8*BTA/((6.28D0 + BTA)*FN*FN)
         RATIO= R1*R2*R3
         XT= XT - RATIO*(X(1) - XT)
      ELSE
C        MIDDLE ZEROS
         XT= X3*X(I-1) - X3*X(I-2) + X(I-3)
      ENDIF
C
      CALL ROOT(XT,NN,ALF,BTA,DPN,PN1,B,C,EPS)
      X(I)= XT
      A(I)= CC/(DPN*PN1)
12    CSA= CSA + A(I) + A(I)
      if (2*nn2 .ne. nn) csa=csa-a(nn2)
      RETURN
      END
      SUBROUTINE ROOT(X,NN,ALF,BTA,DPN,PN1,B,C,EPS)
C
C     IMPROVES THE APPROXIMATE ROOT X; IN ADDITION OBTAINS
C          DPN = DERIVATIVE OF P(N) AT X
C          PN1 = VALUE OF P(N-1) AT X.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION B(NN),C(NN)
      ITER= 0
1     ITER= ITER + 1
      CALL RECUR(P,DP,PN1,X,NN,ALF,BTA,B,C)
      D = P/DP
      X = X - D
      IF(ABS(D) - EPS) 3,3,2
2     IF(ITER - 10) 1,3,3
3     DPN= DP
      RETURN
      END
      SUBROUTINE RECUR(PN,DPN,PN1,X,NN,ALF,BTA,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION B(NN),C(NN)
      DATA X0/0.0D0/,X1/1.0D0/,X2/2.0D0/
      P1= X1
      P= X + (ALF-BTA)/(ALF + BTA + X2)
      DP1= X0
      DP= X1
      DO 1 J=2,NN
      Q= (X - B(J))*P - C(J)*P1
      DQ= (X - B(J))*DP + P - C(J)*DP1
      P1= P
      P= Q
      DP1= DP
1     DP= DQ
      PN= P
      DPN= DP
      PN1= P1
      RETURN
      END
      SUBROUTINE ALLPTS(XALF,WALF,NANG,NANG2)
C
C     TAKES THE POINTS & WEIGHTS GENERATED BY LEGPT FOR THE HALF-RANGE
C     AND CREATES NEW ARRAYS FOR THE FULL-RANGE (-1,+1).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
C
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      DIMENSION XALF(IDVR),WALF(IDVR)
C
      SCALE = DBLE(MAX(1,IDIA))
      DO 10 I=1,NANG2
      WALF(I) = SQRT(SCALE*WALF(I))
      IF (IDIA .EQ. 2) GOTO 10
      XALF(NANG+1-I) = XALF(I)
      XALF(I)        =-XALF(I)
      WALF(NANG+1-I) = WALF(I)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE ASLEG(PLEG,PNORM,LMAX,X,NN2,KZ,LINCR)
C
C     CALCULATE POLYNOMIALS 1 TO LMAX AT X = COS(THETA) FOR M = 0 OR 1,
C     USING THE ROUTINE OF PRESS ET AL, PAGE 182,
C     FOR THE POLYNOMIAL PART OF ASSOCIATED LEGENDRE FUNCTIONS.
C     WE HAVE REMOVED SIN(THETA)**2 FOR NU = 1.
C     THIS ENABLES US TO USE JACOBI INTEGRATION WITH ALF = BTA = NU,
C     USING ROUTINES DERIVED FROM BEIDENHARN AND LOUCK.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PLEG(0:LMAX,NN2),X(NN2),PNORM(0:LMAX)
      DATA X1/1.0D0/,X2/2.0D0/
      M = KZ
      IF(M.LT.0.OR.M.GT.LMAX) GOTO 999
      DO 10 I=1,NN2
      PMM = X1
      FACT = X1
      DO 11 J=1,M
          PMM = -PMM * FACT
          FACT = FACT + X2
   11 CONTINUE
C
      PLEG(0,I) = PMM
      PMMP1= X(I)*(M+M+1)*PMM
      PLEG(1,I)= PMMP1
      LL=1
      DO 2 L= 2+M,LMAX+LINCR
      R2LM1 = DBLE(L+L-1)
      RLPMM1= DBLE(L+M-1)
      RLMM  = DBLE(L-M)
      PLL= (X(I)*R2LM1*PMMP1 - RLPMM1*PMM)/RLMM
      PMM= PMMP1
      PMMP1= PLL
      LL=LL+1
      PLEG(LL,I)= PLL
2     CONTINUE
10    CONTINUE
C
C     SET UP THE NORMALISATION CONSTANTS
C     (PNORM)**2 = (2J + 1)/2   *   (J - K)! / (J + K)!
      JSTART = M
      JJ = -1
      DO 13 J = JSTART,LMAX+LINCR
      FACT = X1
      DO 12 I = J-M+1,J+M
         FACTI = DBLE(I)
         FACT = FACT * FACTI
   12 CONTINUE
      RJ = DBLE(J)
      JJ = JJ + 1
      PNORM(JJ) = (RJ + RJ + X1) / (FACT + FACT)
      PNORM(JJ) = SQRT(PNORM(JJ))
   13 CONTINUE
C     NOW NORMALISE THE POLYNOMIALS
      DO 14 I=1,NN2
         JJ = -1
      DO 14 J=JSTART,LMAX+LINCR
         JJ = JJ + 1
         PLEG(JJ,I) = PLEG(JJ,I) * PNORM(JJ)
   14 CONTINUE
      RETURN
999   WRITE(6,200)
200   FORMAT(/,/,5X,'IMPROPER ARGUMENT IN SUBROUTINE ASLEG',/)
      STOP
      END
      SUBROUTINE LMATRX(XLMATR,PLEG,WALF,KZ,IPAR0,NIDVR,LINCR)
C
C     THIS SUBROUTINE SETS UP THE LOWER TRIANGLE OF THE TRANSFORMED
C     ANGULAR MOMENTUM MATRIX L(ALPHA,ALPHA')
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION XLMATR(IDVR,IDVR),PLEG(0:MAXLEG,IDVR),WALF(IDVR)
      DATA X0/0.0D0/,TWO/2.0D0/
      JSTART = KZ
      JDIA=MAX(1,IDIA)
      JJ0=-JDIA
      IF (MOD(JSTART,JDIA) .NE. IPAR0) THEN
         JJ0=JJ0+1
         JSTART=JSTART+1
      ENDIF
      DO 10 K= 1,NIDVR
      TERM = WALF(K)
      DO 10 KP=K,NIDVR
        SUMJ1=X0
        JJ = JJ0
        DO 20 J=JSTART,MAXLEG+LINCR,JDIA
          JJ = JJ + JDIA
          SUMJ1 = SUMJ1 + PLEG(JJ,K) * PLEG(JJ,KP) * DBLE((J+1)*J)
   20   CONTINUE
        XLMATR(KP,K) = SUMJ1 * TERM * WALF(KP)
        xlmatr(k,kp) = xlmatr(kp,k)
   10 CONTINUE
      IF (.NOT. ZLMAT) RETURN
C     WRITE XLMATR IF REQUESTED
      WRITE(6,1010) KZ,IPAR0
 1010 FORMAT(1H1,5X,'L-MATRIX FOR KZ =',I3,', IPAR =',I2/)
      CALL SQOUT(XLMATR,NIDVR)
      RETURN
      END
      subroutine jhmain(xk1,xk2,xlmatr,r1,r2,xalf,ham1,eig1,iv1,vecs1d,
     &                  eigs1d,work1,ham2,eig2,iv2,eigs2d,ndim2d,
     &                  work2,eigs2,ham3,eval,work3,cint,cintp,hband,
     &                  r1m2t,r2m2t,iv1l,iv2l,ndim2l,vecs1l,vecs2l,
     &                  vecs3l,phi,psi,evall,kz)
c
c     This routine controls the DVR calculation in all cases except
c     symmetrised Radau coordinates.
c     Written by James Henderson
c
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
c
      DIMENSION NDIM2D(NDIMC),EIGS2(MAX2D,NDIMC),R1(NPNT1),R2(NPNT2),
     9          XALF(IDVR),XLMATR(IDVR,IDVR),
     9          XK1(NPNT1,NPNT1),XK2(NPNT2,NPNT2),
     9          HAM1(NDIMA,NDIMA),HBAND(MAX2D,MAX3D),EIGS1D(max1dv),
     9          EIG1(NDIMA),WORK1(NDIMA),IV1(NDIMC,NDIMB),
     9          VECS1D(max1dv,NDIMA),HAM2(MAX2D,MAX2D),EIG2(MAX2D),
     9          WORK2(MAX2D),IV2(NDIMC),HAM3(MAX3D,MAX3D),
     9          CINT(NDIMA*NDIMB,MAX2D),CINTP(NDIMA*NDIMB,MAX2D),
     9          EIGS2D(MAX3D),EVAL(MAX3D),WORK3(MAX3D),
     9          R1M2T(NPNT1,NPNT1),R2M2T(NPNT2,NPNT2),EVALL(NEVAL),
     9          IV1L(NDIMC,NDIMB),IV2L(NDIMC),NDIM2L(NDIMC),
     9          VECS1L(MAX2D,NDIMA),VECS2L(MAX2D,MAX2D),VECS3L(MAX3D),
     9          PHI(MAX3D,NDIMA,NDIMB),PSI(IDVR,NPNT1,NPNT2)
c
      TERM  = DBLE(JROT * JROT + JROT - (2 * KZ * KZ))
C     CONSTRUCT THE ONE-DIMENSIONAL HAMILTONIAN MATRIX H1(NPNTA,NPNTA')
C     FOR EACH NPNTB AND NPNTC, THEN SOLVE BY DIAGONALISATION.
      ICALL = 0
      NSUM = 0
      DO 10 IONE = 1,NPNTC
      NHAM2 = 0
      DO 20 ITWO = 1,NPNTB
      IF (ZR2R1) THEN
         CALL MKHAM1(HAM1,XLMATR,IONE,ITWO,TERM,R1,R2,XALF,XK1,XK2)
      ELSE
         CALL MKHAM1(HAM1,XLMATR,ITWO,IONE,TERM,R1,R2,XALF,XK1,XK2)
      ENDIF
C     DIAGONALISE EACH BLOCK, SAVING THE EIGENVALUES AND VECTORS THAT
C     FALL BELOW THE CUT-OFF ENERGY EMAX1.
      CALL DIAG(HAM1,NDIMA,NPNTA,EIG1,WORK1)
      CALL CUT1D(HAM1,EIG1,IV1(IONE,ITWO),EIGS1D,VECS1D,NHAM2,ICALL)
   20 CONTINUE
      WRITE(6,985) NHAM2
  985 FORMAT(5X,' NHAM2 = ',I4)
      NSUM = NSUM + NHAM2
      IF (IONE .EQ. NPNTC) WRITE(6,986) NSUM
  986 FORMAT(/5X,' SUM = ',I5)
      IF (NHAM2 .GT. 0) THEN
C        DUMP THE 1D EIGENAVLUES & VECTORS TO DISK FOR EACH IONE
         CALL OUTROW(EIGS1D,NHAM2,IEIGS1)
         DO 2 ka=1,NPNTA
         IF (ZVEC) CALL OUTROW(VECS1D(1,ka),NHAM2,IOUT2)
         CALL OUTROW(VECS1D(1,ka),NHAM2,IVECS1)
    2    CONTINUE
      ENDIF
   10 CONTINUE
C
      IF (ZVEC) WRITE (IOUT1) IV1
C
      CALL TIMER
C
C     NOW WANT TO MAKE THE TWO-DIMENSIONAL HAM2(NPNTB,NPNTB',I,I'),
C     WHERE I RUNS OVER THE SELECTED HAM1(NPNTA,NPNTA') SOLUTIONS.
C
      REWIND IVECS1
      REWIND IEIGS1
      ICALL = 0
      LOW3D = 0
C     if ZALL then set ZCUT for convenience...
      IF (ZALL) ZCUT = .TRUE.
      DO 31 IONE = 1,NPNTC
      NHAM2 = 0
C     RECALL THE SIZE OF HAM2
      DO 3 ITWO = 1,NPNTB
      NHAM2 = NHAM2 + IV1(IONE,ITWO)
    3 CONTINUE
      NDIM2D(IONE) = NHAM2
C
      IF ( NHAM2 .GT. 0 ) THEN
        CALL MKHAM2(HAM2,EIGS1D,VECS1D,XK1,XK2,IV1,IONE,NHAM2)
        if (.not. ztwod) then
           CALL DIAG(HAM2,MAX2D,NHAM2,EIG2,WORK2)
        else
           CALL DIAG3D(HAM2,NHAM2,eig2,WORK2,kz)
           return
        endif
C
        IF (ZCUT) THEN
          CALL CUT2D(HAM2,EIG2,IV2(IONE),NHAM2,LOW3D,ICALL)
        ELSE
          DO 4 II = 1,NHAM2
          CALL OUTROW(HAM2(1,II),NHAM2,INTVEC)
    4     EIGS2(II,IONE) = EIG2(II)
          LOW3D = MAX3D
        ENDIF
C
      ENDIF
C
      IF (IONE .EQ. NPNTC .AND. ZCUT) WRITE(6,987)  LOW3D,EMAX2
  987 FORMAT(/I14,' EIGENVALUES SELECTED BELOW ',D20.10)
C
   31 CONTINUE
C
C
      IF (.NOT. ZCUT) CALL CHOOSE(EIGS2,NDIM2D,HAM2,IV2,NHAM2,LOW3D)
C
      IF (ZVEC) WRITE (IOUT1) IV2
C
      CALL TIMER
C
      NHAM3 = LOW3D
C
C     SAVE THE REQUIRED BITS TO DISK IF ZDIAG = .FALSE.
C     first open the required disk file
      IF (.NOT. ZDIAG) THEN
        OPEN(UNIT=IDIAG1,FORM='UNFORMATTED')
        OPEN(UNIT=IDIAG2,FORM='UNFORMATTED')
        WRITE(IDIAG1) NPNTA,NPNTB,NPNTC,MAX2D,NHAM3,ZCUT
        CALL IOUTRO(IV1,NDIMB*NDIMC,IDIAG1)
        CALL IOUTRO(IV2,NDIMC,IDIAG1)
C
C       NEED ALSO THE 2-D EIGENVALUES TO BUILD THE FINAL HAMILTONIAN
        IF (ZCUT) THEN
          REWIND IEIGS2
          DO 55 I1 = 1,NPNTC
          IV = IV2(I1)
          IF (IV .GT. 0) CALL GETROW(EIG2,IV,IEIGS2)
          IF (IV .GT. 0) CALL OUTROW(EIG2,IV,IDIAG1)
   55     CONTINUE
        ELSE
          CALL OUTROW(EIGS2,MAX2D*NDIMC,IDIAG1)
        ENDIF
C
      ENDIF
C
      CALL MKHAM3(HAM3,NHAM3,HAM2,NHAM2,CINT,CINTP,IV1,IV2,WORK3,
     1            VECS1D,XK1,XK2,EIGS2D,EIGS2,HBAND,XLMATR,R1,R2,
     1            R1M2T,R2M2T,TERM)
C
      CALL TIMER
C
      CALL DIAG3D(HAM3,NHAM3,EVAL,WORK3,kz)
C
      CALL TIMER
C
C.....finally, can compute the actual wavefunction amplitude at
C.....the grid points if needed.
      IF (ZTRAN) THEN
         CALL TRANS(IV1L,IV2L,NDIM2L,VECS1L,
     1              VECS2L,VECS3L,PHI,PSI,EVALL,NHAM3)
         CALL TIMER
      ENDIF
      return
      end
      SUBROUTINE MKHAM1(HAM1,XLMATR,I1,I2,TERM,R1,R2,XALF,XK1,XK2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
      DIMENSION HAM1(NDIMA,NDIMA),XLMATR(IDVR,IDVR),R1(NPNT1),R2(NPNT2),
     1          XALF(IDVR),XK1(NPNT1,NPNT1),XK2(NPNT2,NPNT2)
      DATA XP5/0.50D0/
C
C
C     ZERO HAM1
      DO 5 K=1,NPNTA
      DO 5 KP=1,K
    5 HAM1(K,KP) = 0.0D0
c     zero rotational excitation term for J=0 cases
      WTERM = 0.0D0
c.....theta first
      IF (ZTHETA) THEN
c
         W1GAMA = XP5 / (R1(I1)*R1(I1)*UR1)
         W2BETA = XP5 / (R2(I2)*R2(I2)*UR2)
         WSUM = W1GAMA + W2BETA
C
         IF (JROT .GT. 0) THEN
            if (zembed) then
C              have TERM * R2**(-2) term
               WTERM = TERM * W2BETA
            else
C              have TERM * R1**(-2) term
               WTERM = TERM * W1GAMA
            endif
         ENDIF
      ENDIF
C
      DO 10 K = 1,NPNTA
      IF (ZTHETA) THEN
         IF (ZLPOT) THEN
            CALL POTVLG(V,R1(I1),R2(I2),XALF(K))
         ELSE
            CALL POTV(V,R1(I1),R2(I2),XALF(K))
         ENDIF
      ELSE
         IF (ZLPOT) THEN
            IF (ZR2R1)THEN
               CALL POTVLG(V,R1(I2),R2(K),XALF(I1))
            ELSE
               CALL POTVLG(V,R1(K),R2(I1),XALF(I2))
            ENDIF
         ELSE
            IF (ZR2R1)THEN
               CALL POTV(V,R1(I2),R2(K),XALF(I1))
            ELSE
               CALL POTV(V,R1(K),R2(I1),XALF(I2))
            ENDIF
         ENDIF
      ENDIF
      IF (ZTHETA) THEN
         HAM1(K,K) = V + WTERM
         DO 11 KP= 1,K
         HAM1(K,KP) = HAM1(K,KP) + XLMATR(K,KP)*WSUM
   11    CONTINUE
      ELSE
         if (jrot .gt. 0) then
            IF (ZEMBED) THEN
              IF (ZR2R1)THEN
                 WTERM = (TERM * XP5) / (R2(K)*R2(K)*UR2)
              ELSE
                 WTERM = (TERM * XP5) / (R2(I1)*R2(I1)*UR2)
              ENDIF
            ELSE
              IF (ZR2R1)THEN
                 WTERM = (TERM * XP5) / (R1(I2)*R1(I2)*UR1)
              ELSE
                 WTERM = (TERM * XP5) / (R1(K)*R1(K)*UR1)
              ENDIF
            ENDIF
         ENDIF
         HAM1(K,K) = V + WTERM
         IF (ZR2R1) THEN
            DO 12 KP= 1,K
            HAM1(K,KP) = HAM1(K,KP) + XK2(K,KP)
   12       CONTINUE
         ELSE
            DO 13 KP= 1,K
            HAM1(K,KP) = HAM1(K,KP) + XK1(K,KP)
   13       CONTINUE
         ENDIF
      ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE MKHAM2(HAM2,EIGS1D,VECS1D,XK1,XK2,IV1,IONE,NHAM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION HAM2(MAX2D,MAX2D),XK2(NPNT2,NPNT2),IV1(NDIMC,NDIMB),
     1          EIGS1D(MAX2D),VECS1D(MAX2D,NDIMA),XK1(NPNT1,NPNT1)
C
C     ZERO LOWER TRIANGLE OF HAM2
      DO 20 I = 1,NHAM2
      DO 20 J = 1,I
   20 HAM2(I,J) = 0.0D0
      DO 30 I3= 1,NPNTA
         CALL GETROW(VECS1D(1,I3),NHAM2,IVECS1)
      DO 40 J = 1,NHAM2
      DO 40 K = 1,J
      HAM2(J,K) = HAM2(J,K) + VECS1D(J,I3)*VECS1D(K,I3)
   40 CONTINUE
   30 CONTINUE
C     MUST NOW MULTIPLY BY XK1 OR XK2
      IVBSM = 0
      DO 50 ITWO = 1,NPNTB
      IVBPSM = 0
      IVB = IV1(IONE,ITWO)
      DO 60 ITWOP = 1,ITWO
      IF (ZTHETA) THEN
        IF (ZR2R1) THEN
           XKTERM = XK2(ITWO,ITWOP)
        ELSE
           XKTERM = XK1(ITWO,ITWOP)
        ENDIF
      ELSE
        IF (ZR2R1) THEN
           XKTERM = XK1(ITWO,ITWOP)
        ELSE
           XKTERM = XK2(ITWO,ITWOP)
        ENDIF
      ENDIF
      IVBP = IV1(IONE,ITWOP)
      DO 70 J  = 1,IVB
      IND1 = IVBSM + J
      DO 70 JP = 1,IVBP
      IND2 = IVBPSM + JP
         HAM2(IND1,IND2) = HAM2(IND1,IND2) * XKTERM
   70 CONTINUE
      IVBPSM = IVBPSM + IVBP
   60 CONTINUE
      IVBSM = IVBSM + IVB
   50 CONTINUE
C
C     NOW ADD THE 1-D EIGENVALUES ALONG THE DIAGONAL
      CALL GETROW(EIGS1D,NHAM2,IEIGS1)
      DO 80 NN = 1,NHAM2
         HAM2(NN,NN) = HAM2(NN,NN) + EIGS1D(NN)
   80 CONTINUE
      RETURN
      END
      SUBROUTINE MKHAM3(HAM3,NHAM3,HAM2,NHAM2,CINT,CINTP,IV1,IV2,WORK3,
     1                  VECS1D,XK1,XK2,EIGS2D,EIGS2,HBAND,XLMATR,
     2                  R1,R2,R1M2T,R2M2T,TERM)
C
C     BUILD THE FINAL 3-D HAMILTONIAN MATRIX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
      DIMENSION CINT(NDIMA*NDIMB,NHAM2),HAM2(MAX2D,MAX2D),IV2(NDIMC),
     1          IV1(NDIMC,NDIMB),VECS1D(MAX2D,NDIMA),
     2          CINTP(NDIMA*NDIMB,NHAM2),HAM3(MAX3D,NHAM3),
     3          XK1(NPNT1,NPNT1),EIGS2D(NHAM3),XK2(NPNT2,NPNT2),
     4          HBAND(MAX2D,MAX3D),EIGS2(MAX2D,NDIMC),
     5          XLMATR(IDVR,IDVR),R1(NPNT1),R2(NPNT2),WORK3(NHAM3),
     6          R1M2T(NPNT1,NPNT1),R2M2T(NPNT2,NPNT2)
      DATA XP5/0.50D0/
C
C     IF ZDIAG = .FALSE. WANT EIGS2D NOW
      IF (.NOT. ZDIAG) THEN
        REWIND IEIGS2
        NEIG = 0
        DO 184 IONE = 1,NPNTC
        IV = IV2(IONE)
        IF (.NOT. ZCUT) THEN
          DO 19 II = 1,IV
          NEIG = NEIG + 1
   19     EIGS2D(NEIG) = EIGS2(II,IONE)
        ENDIF
        IF (IV .GT. 0 .AND. ZCUT) CALL GETROW(EIGS2D,IV,IEIGS2)
  184   CONTINUE
      ENDIF
C
C     FIRST DO THE INTERMEDIATE TRANSFORMATION
      NCINT = NPNTA*NPNTB
      NDIMT = NDIMA*NDIMB
      REWIND IVECS2
      REWIND IVECS1
      DO 10 IONE = 1,NPNTC
C     RECALL THE SIZE OF THE 2-D VECTORS
      NHAM2=0
      DO 3 ITWO=1,NPNTB
    3 NHAM2 = NHAM2 + IV1(IONE,ITWO)
C
      IF (NHAM2 .EQ. 0) GOTO 10
C
      DO 4 JCINT = 1,NHAM2
      DO 4 ICINT = 1,NDIMT
    4 CINT(ICINT,JCINT) = 0.0D0
C
C     BRING BACK THE 1-D VECTORS FOR EACH I1
      DO 23 KK = 1,NPNTA
         IF (IV2(IONE) .GT. 0) THEN
            CALL GETROW(VECS1D(1,KK),NHAM2,IVECS1)
         ELSE
            READ (IVECS1)
         ENDIF
   23 CONTINUE
C
      DO 20 J = 1,IV2(IONE)
C     BRING BACK THE 2-D VECTORS FOR EACH NPNTC
      CALL GETROW(HAM2(1,J),NHAM2,IVECS2)
      IND2 = 0
      DO 30 K = 1,NPNTA
      IND1 = 0
      DO 40 ITWO = 1,NPNTB
      IND2 = IND2 + 1
      DO 50 I = 1,IV1(IONE,ITWO)
      IND1 = IND1 + 1
         CINT(IND2,J) = CINT(IND2,J) + HAM2(IND1,J)*VECS1D(IND1,K)
   50 CONTINUE
   40 CONTINUE
   30 CONTINUE
      IF (.NOT. ZTHETA) THEN
         IN2 = 0
         DO 37 IA = 1,NPNTA
         DO 47 IB = 1,NPNTB
         IF (ZR2R1) THEN
            II1 = IB
            II2 = IA
         ELSE
            II1 = IA
            II2 = IB
         ENDIF
         W1 = XP5 / (R1(II1)*R1(II1)*UR1)
         W2 = XP5 / (R2(II2)*R2(II2)*UR2)
         WSUM = W1 + W2
         IN2 = IN2 + 1
         CINT(IN2,J) = CINT(IN2,J)*SQRT(WSUM)
   47    CONTINUE
   37    CONTINUE
      ENDIF
   20 CONTINUE
C     STORE CINT ON DISK FOR EACH NPNTC
      IF (IV2(IONE) .GT. 0) CALL OUTROW(CINT,NDIMT*IV2(IONE),IVINT)
C
   10 CONTINUE
C
C     NOW DO THE SECOND PART OF THE TRANSFORMATION
C
      REWIND IVINT
      IF (ZTHETA) THEN
      LENGTH = 0
      DO 61 IONE = 1,NPNTC
C.....for each npntc....................................
C.....
C     SET THE HBAND TO ZERO
      DO 12 L = 1,MAX2D
      DO 12 LP= 1,NHAM3
   12 HBAND(L,LP) = 0.0D0
      IVSM = 0
C
      NHAM2 = 0
      DO 5 ITWO = 1,NPNTB
    5 NHAM2 = NHAM2 + IV1(IONE,ITWO)
      IF (NHAM2 .EQ. 0) GOTO 61
C
      IV = IV2(IONE)
      IF (IV .GT. 0) CALL GETROW(CINT,NDIMT*IV,IVINT)
      REWIND IVINT
      IVPSM = 0
      DO 51 IONEP = 1,IONE
      IF (IONE .EQ. IONEP) THEN
        IF (ZQUAD2) THEN
          RM2T = 0.0D0
        ELSE
          D1R2 = R2M2T(IONE,IONEP)
          D2R2 = XP5 / (R2(IONE)*R2(IONE)*UR2)
          RM2T = D1R2 - D2R2
        ENDIF
      ELSE
        IF (ZQUAD2) THEN
          RM2T = 0.0D0
        ELSE
          RM2T = R2M2T(IONE,IONEP)
        ENDIF
      ENDIF
      IF (ZR2R1) THEN
         XKTERM = XK1(IONE,IONEP)
      ELSE
         XKTERM = XK2(IONE,IONEP)
         IF(JROT .GT. 0 .AND. .NOT. ZQUAD2 .AND. ZEMBED)
     1    XKTERM = XKTERM + TERM*RM2T
      ENDIF
      IVP = IV2(IONEP)
      IF (IVP .GT. 0) CALL GETROW(CINTP,NDIMT*IVP,IVINT)
      DO 41 J=1,IV
      IND1 = IVSM + J
      DO 31 JP =1,IVP
      IND2 = IVPSM + JP
C
      IND3 = 0
      DO 81 K = 1,NCINT
      IND3 = IND3 + 1
      HBAND(IND1,IND2) = HBAND(IND1,IND2)
     1                 + CINT(IND3,J)*CINTP(IND3,JP)*XKTERM
   81 CONTINUE
C
      IND3 = 0
      DO 83 NA = 1,NPNTA
      DO 82 K = 1,NPNTB
      IND3 = IND3 + 1
      IND4 = K
      DO 84 NAP = 1,NPNTA
      HBAND(IND1,IND2) = HBAND(IND1,IND2) + CINT(IND3,J)*CINTP(IND4,JP)
     1                 * XLMATR(NA,NAP) * RM2T
      IND4 = IND4 + NPNTB
   84 CONTINUE
   82 CONTINUE
   83 CONTINUE
C
   31 CONTINUE
   41 CONTINUE
      IVPSM = IVPSM + IVP
   51 CONTINUE
C
C
C
      IF (ZDIAG) THEN
         DO 71 JJ=1,IND2
         IF (IV .GT. 0) CALL OUTROW(HBAND(1,JJ),IV,IBAND)
   71    CONTINUE
      ELSE
         DO 93 ISAVE = 1,IV
         LENGTH = LENGTH + 1
         DO 94 JSAVE = 1,IND2
         WORK3(JSAVE)= HBAND(ISAVE,JSAVE)
   94    CONTINUE
         WORK3(LENGTH) = WORK3(LENGTH) + EIGS2D(LENGTH)
         CALL OUTROW(WORK3,LENGTH,IDIAG2)
   93    CONTINUE
      ENDIF
   61 CONTINUE
C
C
      ELSE
C
C
      LENGTH = 0
      DO 66 IONE = 1,NPNTC
C.....for each npntc....................................
C.....
C     SET THE HBAND TO ZERO
      DO 17 L = 1,MAX2D
      DO 17 LP= 1,NHAM3
   17 HBAND(L,LP) = 0.0D0
      IVSM = 0
C
      NHAM2 = 0
      DO 8 ITWO = 1,NPNTB
    8 NHAM2 = NHAM2 + IV1(IONE,ITWO)
      IF (NHAM2 .EQ. 0) GOTO 66
C
      IV = IV2(IONE)
      IF (IV .GT. 0) CALL GETROW(CINT,NDIMT*IV,IVINT)
      REWIND IVINT
      IVPSM = 0
      DO 56 IONEP = 1,IONE
      XKTERM = XLMATR(IONE,IONEP)
      IVP = IV2(IONEP)
      IF (IVP .GT. 0) CALL GETROW(CINTP,NDIMT*IVP,IVINT)
      DO 46 J=1,IV
      IND1 = IVSM + J
      DO 36 JP =1,IVP
      IND2 = IVPSM + JP
C
      IND3 = 0
      DO 86 K = 1,NCINT
      IND3 = IND3 + 1
           HBAND(IND1,IND2) = HBAND(IND1,IND2)
     1                      + CINT(IND3,J)*CINTP(IND3,JP)*XKTERM
   86 CONTINUE
C
   36 CONTINUE
   46 CONTINUE
      IVPSM = IVPSM + IVP
   56 CONTINUE
      IF (ZDIAG) THEN
         DO 76 JJ=1,IND2
         IF (IV .GT. 0) CALL OUTROW(HBAND(1,JJ),IV,IBAND)
   76    CONTINUE
      ELSE
         DO 95 ISAVE = 1,IV
         LENGTH = LENGTH + 1
         DO 96 JSAVE = 1,IND2
         WORK3(JSAVE)= HBAND(ISAVE,JSAVE)
   96    CONTINUE
         WORK3(LENGTH) = WORK3(LENGTH) + EIGS2D(LENGTH)
         CALL OUTROW(WORK3,LENGTH,IDIAG2)
   95    CONTINUE
      ENDIF
   66 CONTINUE
C
      ENDIF
C
      IF (.NOT. ZDIAG) THEN
         WRITE(6,1080)
 1080    FORMAT(/10X,'HAMILTONIAN WRITTEN TO DISK - NOT DIAGONALISED')
         WRITE(6,1081) IDIAG2
 1081    FORMAT(/10X,'HAMILTONIAN BANDS WRITTEN TO STREAM IDIAG =',I4/)
         STOP
      ENDIF
C
      DO 77 I=1,NHAM3
      DO 77 IP=1,I
   77 HAM3(I,IP) = 0.0D0
      REWIND IBAND
      IVSM = 0
      DO 62 IONE = 1,NPNTC
      IVBAND = 0
C
      NHAM2 = 0
      DO 6 ITWO = 1,NPNTB
    6 NHAM2 = NHAM2 + IV1(IONE,ITWO)
      IF (NHAM2 .EQ. 0) GOTO 62
C
      IV = IV2(IONE)
      IVPSM = 0
      DO 52 IONEP = 1,IONE
      IVP = IV2(IONEP)
      DO 42 J=1,IV
      IND1 = IVSM + J
      IND1BD = IVBAND + J
      DO 32 JP =1,IVP
      IND2 = IVPSM + JP
C
   32 CONTINUE
   42 CONTINUE
      IVPSM = IVPSM + IVP
   52 CONTINUE
      DO 72 JJ=1,IND2
      IF (IV .GT. 0) CALL GETROW(HAM3(IVSM + 1,JJ),IV,IBAND)
   72 CONTINUE
      IVSM = IVSM + IV
   62 CONTINUE
C
C
C     NOW PUT THE 2-D EIGENSOLUTIONS ALONG THE DIAGONAL
      REWIND IEIGS2
      NN = 0
      DO 183 IONE = 1,NPNTC
      IV = IV2(IONE)
      IF (.NOT. ZCUT) THEN
      DO 9 II = 1,IV
    9 EIGS2D(II) = EIGS2(II,IONE)
      ENDIF
      IF (IV .GT. 0 .AND. ZCUT) CALL GETROW(EIGS2D,IV,IEIGS2)
      DO 183 I = 1,IV
      NN = NN + 1
         HAM3(NN,NN) = HAM3(NN,NN) + EIGS2D(I)
  183 CONTINUE
C
      IF (ZPHAM) CALL WRTHAM(HAM3,NHAM3)
C
      RETURN
      END
      SUBROUTINE DIAG(HAM,MAXHAM,NHAM,EIG,WORK)
C
C     DIAGONALISE THE APPROPRIATE HAMILTONIAN MATRICES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      DIMENSION HAM(MAXHAM,NHAM),EIG(NHAM),WORK(NHAM)
C
      IFAIL=0
      CALL F02ABF(HAM,MAXHAM,NHAM,EIG,HAM,MAXHAM,WORK,IFAIL)
C
      if (ifail .ne. 0) write(6,100) ifail
      return
100   format(' Diagonalisation has failed with, IFAIL=',i3)
      end
      SUBROUTINE DIAG3D(HAM3,NHAM3,EVAL,WORK3,kz)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION HAM3(MAX3D,NHAM3),EVAL(NHAM3),WORK3(NHAM3)
C     AUTOCM CONVERTS ATOMIC UNITS (HARTREE) TO CM-1.
      DATA AUTOCM/2.19474624D+05/
      DATA X0/0.0D0/
      IF (ZROT) then
         WRITE(6,1040) JROT,KZ
 1040    FORMAT(/10X,'Solutions with J =',I3,' k =',I3)
         IF (idia .EQ. -2) THEN
            IF (IPAR .EQ. 0)  WRITE(6,1025)
 1025       FORMAT(/10X,'EVEN PARITY SOLUTIONS')
            IF (IPAR .EQ. 1)  WRITE(6,1035)
 1035       FORMAT(/10X,'ODD PARITY SOLUTIONS')
         ENDIF
      ENDIF
C
      MEVAL=MIN(NEVAL,NHAM3)
C.....we need this in special cases (NHAM3 < NEVAL)
      IF (ZVEC) WRITE(IOUT2) MEVAL
C
C     DIAGONALISE  THE HAMILTONIAN. THE VECTORS ARE TO OVERWRITE THE
C     HAMILTONIAN.
      IFAIL=1
      nvecln = max3d
      if(idia .eq. -2) nvecln = nham3
      CALL F02ABF(HAM3,nvecln,NHAM3,EVAL,HAM3,nvecln,WORK3,IFAIL)
C     PRINT EIGENVALUES IN ATOMIC UNTIS & WAVENUMBERS
C
      IF (.NOT.ZPMIN) THEN
         WRITE(6,1010) MEVAL
 1010    FORMAT(1H1/5X,'LOWEST',I4,' EIGENVALUES IN HARTREES:',/)
         WRITE(6,1020) (EVAL(I),I=1,MEVAL)
      ENDIF
      IF (ZPFUN) THEN
         IP=JROT-KMIN
         IF (JROT .NE. 0 .AND. IP .NE. 1) GOTO 10
         JDIA=MAX(0,IDIA)
         JPAR=MIN(JDIA,IPAR)
         ISYM=ABS(MIN(0,IDIA))
         IF (IPAR .EQ. 1) ISYM=-ISYM
         WRITE(ILEV,1125) JROT,IP,JDIA,JPAR,ISYM,MEVAL
 1125    FORMAT(6(1X,I6)
         WRITE(ILEV,1126) (EVAL(I),I=1,MEVAL)
 1126    FORMAT(4D20.12)
      ENDIF
   10 CONTINUE
C
C     SAVE THE EIGENVALUES IF NEEDED
      IF (ZVEC) WRITE (IOUT2) (EVAL(I),I=1,MEVAL)
      DO 20 I=1,MEVAL
      WORK3(I) = EVAL(I) * AUTOCM
   20 CONTINUE
      WRITE(6,1030) MEVAL
 1030 FORMAT(//5X,'LOWEST',I4,' EIGENVALUES IN WAVENUMBERS:'/)
      WRITE(6,1020) (WORK3(I),I=1,MEVAL)
 1020 FORMAT(5D24.12/)
C
C     IF REQUESTED PRINT THE EIGENVECTORS
      IF (ZPVEC) THEN
         WRITE(6,1050)
 1050    FORMAT(//'  EIGENVECTORS'/)
         DO 40 I=1,MEVAL
         WRITE(6,1060) (HAM3(J,I),J=1,NHAM3)
 1060    FORMAT(/(1X,10F13.7))
   40    CONTINUE
      ENDIF
C     WRITE THE FINAL VECTORS TO DISK IF REQUIRED
      IF (ZVEC) THEN
         DO 60 L = 1,MEVAL
         CALL OUTROW(HAM3(1,L),NHAM3,IOUT2)
   60    CONTINUE
      ENDIF
      IF (JROT .NE. 0) RETURN
      IF (abs(IDIA) .EQ. 2 .AND. IPAR .EQ. 1) THEN
         II=1
         ezero=x0
         read(5,5,end=55) ezero
    5    format(f20.0)
   55    continue
      ELSE
         E0=WORK3(1)
         II=2
      ENDIF
      WRITE(6,1070)
 1070 FORMAT(///'  BAND ORIGINS IN WAVENUMBERS:'/)
      WRITE(6,1020) (WORK3(I)-E0,I=II,MEVAL)
      RETURN
      END
      SUBROUTINE CHOOSE(EIGS2,NDIM2D,HAM2,IV2,NHAM2,LOW3D)
C
C     THIS ROUTINE CHOOSES THE MAX3D LOWEST EIGENVALUES FROM EIGS2.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION EIGS2(MAX2D,NDIMC),IV2(NDIMC),NDIM2D(NDIMC),
     1          HAM2(MAX2D,MAX2D)
      DATA AUTOCM/2.19474624D+05/
C
      EIGMIN = EIGS2(1,1)
      NHAMSM = 0
      DO 160 I=1,NPNTC
      EIGMIN = MIN(EIGMIN,EIGS2(1,I))
      NHAMSM = NHAMSM + NDIM2D(I)
  160 IV2(I) = 1
      IF (NHAMSM .LT. MAX3D) LOW3D = NHAMSM
      IPT = 1
      DO 200 N=1,LOW3D
  210 IF(IV2(IPT) .LE. NDIM2D(IPT)) THEN
          EIGVIB = EIGS2(IV2(IPT),IPT)
          JPT = IPT
      ELSE
          IPT = IPT + 1
          GOTO 210
      ENDIF
      DO 220 J=IPT+1,NPNTC
      IF(IV2(J) .GT. NDIM2D(J)) GOTO 220
      IF(EIGS2(IV2(J),J) .GE. EIGVIB) GOTO 220
      EIGVIB = EIGS2(IV2(J),J)
      JPT = J
  220 CONTINUE
C     KEEP THE EIGENVALUE:
      IV2(JPT) = IV2(JPT) + 1
  200 CONTINUE
C     STORE THE NUMBER OF EIGENVALUES SELECTED FOR EACH ALPHA
      IV2(1) = IV2(1) - 1
      IVIB = IV2(1)
      WRITE(6,798)
  798 FORMAT(//5X,'SELECTION OUTCOME FOR THE FINAL, CONTRACTED BASIS:')
C     WRITE(6,799)
C 799 FORMAT(5X,'(STARTING FROM THE THETA = 180 END)')
      WRITE(6,800)  IV2(1)
  800 FORMAT(/5X,'NPNTC =  1', ',', ' NO. OF EIGENVECTORS =',I3 /)
      DO 230 I=2,NPNTC
      IV2(I) = IV2(I) - 1
      WRITE(6,900)  I,IV2(I)
  900 FORMAT(5X,'NPNTC =',I3, ',', ' NO. OF EIGENVECTORS =',I3 /)
  230 IVIB = MAX(IVIB,IV2(I))
C
      WRITE(6,998)  LOW3D,EIGMIN,EIGVIB
  998 FORMAT(/I14,' EIGENVALUES SELECTED FROM ',D20.10,' TO',D20.10,
     &    ' HARTREES')
      WRITE(6,999)  LOW3D,EIGMIN*AUTOCM,EIGVIB*AUTOCM
  999 FORMAT(/I14,' EIGENVALUES SELECTED FROM ',D20.10,' TO',D20.10,
     &    ' CM-1')
C
      IF (ZP2D) THEN
         WRITE(6,1051)
 1051    FORMAT(//5X,'2D EIGENVALUES IN WAVENUMBERS:'/)
         DO 32 JONE = 1,NPNTC
         IVJ = IV2(JONE)
         WRITE(6,1050) (EIGS2(I,JONE)*AUTOCM,I=1,IVJ)
 1050    FORMAT(5D24.12/)
         WRITE(6,1052)
 1052    FORMAT(/5X,'------------------------------'/)
   32    CONTINUE
      ENDIF
C
C     SAVE THE VECTORS CHOSEN
      REWIND INTVEC
      DO 31 IONE = 1,NPNTC
      IVM = IV2(IONE)
      NHAM2 = NDIM2D(IONE)
      DO 7 I = 1,NHAM2
    7 CALL GETROW(HAM2(1,I),NHAM2,INTVEC)
      DO 30 IND = 1,IVM
      IF (ZVEC) THEN
         IF (NHAM2 .GT. 0) CALL OUTROW(HAM2(1,IND),NHAM2,IOUT2)
      ENDIF
         IF (NHAM2 .GT. 0) CALL OUTROW(HAM2(1,IND),NHAM2,IVECS2)
   30 CONTINUE
   31 CONTINUE
      RETURN
      END
      SUBROUTINE CUT1D(HAM1,EIG1,IVN,EIGS1D,VECS1D,NHAM2,ICALL)
C
C     THIS ROUTINE SELECTS ALL THE EIGENVALUES THAT ARE LOWER THAN THE
C     THE CUT-OFF ENERGY EMAX1, WHICH IS USER-SUPPLIED IN WAVENUMBERS.
C     THESE EIGENVALUES & THEIR CORRESPONDING VECTORS ARE THEN SAVED
C     TO DISK.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION HAM1(NDIMA,NDIMA),EIG1(NDIMA),
     1          EIGS1D(MAX2D),VECS1D(MAX2D,NDIMA)
      DATA AUTOCM/2.19474624D+05/
C     CHANGE EMAX1 TO HARTREE FOR THE SELECTION
      EMAXAU=EMAX1/AUTOCM
C
      ICALL = ICALL + 1
      IF (.NOT. ZALL) THEN
         IVN = 0
         DO 10 N=1,NPNTA
         IF (EIG1(N) .GT. EMAXAU) GOTO 20
         IVN = IVN + 1
   10    CONTINUE
   20    IV = IVN
      ELSE
         IVN = NPNTA
         IV  = IVN
      ENDIF
C
      IF (ZP1D) THEN
         IF (ICALL .EQ. 1) WRITE(6,1051)
 1051    FORMAT(//5X,'1D EIGENVALUES IN WAVENUMBERS:'/)
         WRITE(6,1050) (EIG1(I)*AUTOCM,I=1,NPNTA)
 1050    FORMAT(5D24.12/)
         WRITE(6,1052)
 1052    FORMAT(/5X,'------------------------------'/)
      ENDIF
C
C     SAVE THE VECTORS AND EIGENVALUES (OVERWRITE FOR EACH GAMMA).
      DO 30 I=1,IV
         NHAM2 = NHAM2 + 1
C
      IF (NHAM2 .GT. MAX2D .AND. .NOT. ZALL)  THEN
         WRITE(6,999)
  999    FORMAT(//6X,'**** CORE EXCEEDED: REDUCE CUT-OFF EMAX1 ****')
         STOP
      ENDIF
C
         EIGS1D(NHAM2)   = EIG1(I)
      DO 31 J=1,NPNTA
         VECS1D(NHAM2,J) = HAM1(J,I)
   31 CONTINUE
   30 CONTINUE
C
      RETURN
      END
      SUBROUTINE CUT2D(HAM2,EIG2,IVM,NHAM2,LOW3D,ICALL)
C
C     THIS ROUTINE SELECTS ALL THE EIGENVALUES THAT ARE LOWER THAN THE
C     THE CUT-OFF ENERGY EMAX2, WHICH IS USER-SUPPLIED IN WAVENUMBERS.
C     THESE EIGENVALUES & THEIR CORRESPONDING VECTORS ARE THEN SAVED
C     TO DISK.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION HAM2(MAX2D,MAX2D),EIG2(NHAM2)
      DATA AUTOCM/2.19474624D+05/
C     CHANGE EMAX2 TO HARTREE FOR THE SELECTION
      EMAXAU=EMAX2/AUTOCM
C
      ICALL = ICALL + 1
      IF (.NOT. ZALL) THEN
         IVM = 0
         DO 10 N=1,NHAM2
         IF(EIG2(N) .GT. EMAXAU) GOTO 20
         IVM = IVM + 1
   10    CONTINUE
   20    LOW3D = LOW3D + IVM
      ELSE
         IVM = NHAM2
         LOW3D = LOW3D + IVM
      ENDIF
C
      IF (LOW3D .GT. MAX3D .AND. .NOT. ZALL)  THEN
         WRITE(6,999)
  999    FORMAT(//6X,'**** CORE EXCEEDED: REDUCE CUT-OFF EMAX2 ****')
         STOP
      ENDIF
C
      IF (ZP2D) THEN
         IF (ICALL .EQ. 1) WRITE(6,1051)
 1051    FORMAT(//5X,'2D EIGENVALUES IN WAVENUMBERS:'/)
         WRITE(6,1050) (EIG2(I)*AUTOCM,I=1,NHAM2)
 1050    FORMAT(5D24.12/)
         WRITE(6,1052)
 1052    FORMAT(/5X,'------------------------------'/)
      ENDIF
C
C     SAVE THE VECTORS AND EIGENVALUES
      IF (IVM .GT. 0) CALL OUTROW(EIG2,IVM,IEIGS2)
      DO 30 IND = 1,IVM
      IF (ZVEC) CALL OUTROW(HAM2(1,IND),NHAM2,IOUT2)
      CALL OUTROW(HAM2(1,IND),NHAM2,IVECS2)
   30 CONTINUE
      RETURN
      END
      subroutine nfmain(hr,htheta,r,theta,ham1,eig1,iv1,vecs1d,
     &                  eigs1d,work1,ham2,eig2,iv2,vecs2d,eigs2d,nv2,
     &                  work2,eigtmp,ham3,eig3,work3,wvfunc,kz)
c
c     This routine controls the DVR calculation in the case of
c     symmetrised Radau coordinates.
c     Written by Nic Fulton, Feb 1993.
c
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
 
      dimension hr(npnt,npnt)
      dimension htheta(nalf,nalf)
      dimension r(npnt)
      dimension theta(nalf)
      dimension ham1(npnt-ipar,npnt-ipar)
      dimension eig1(npnt-ipar)
      dimension iv1(2,npnt,nalf)
      dimension vecs1d(npnt-ipar,max1dv)
      dimension eigs1d(max1dv)
      dimension work1(npnt)
      dimension ham2(max2d,max2d)
      dimension eig2(max2d)
      dimension iv2(2,nalf)
      dimension vecs2d(max2d,max3d)
      dimension eigs2d(max3d)
      dimension nv2(nalf)
      dimension work2(max2d)
      dimension eigtmp(max2d*nalf)
      dimension ham3(max3d,max3d)
      dimension eig3(max3d)
      dimension work3(max3d)
      dimension wvfunc(*)
      DATA X4/4.0D0/,X8/8.0D0/,X16/1.6D1/
 
      IF (jrot .ne. 0) THEN
         TERM  = DBLE(JROT * JROT + JROT - KZ * KZ) / X8
         TERM2 = DBLE(JROT * JROT + JROT - 3 * KZ * KZ) / X4
         IF (ABS(KZ) .EQ. 1) THEN
            TERM3 = DBLE(JROT * JROT + JROT) / X16
            IF (KMIN .GE. 1) TERM3 = -TERM3
         ENDIF
      endif
 
c No need for the 1D diagonalisations as there is no possiblity of
c truncation as the symmetry would be broken.
      write(6,130)
      call nftim('beginning of 2d loop')
      do 30 igamma = 1,nalf
        nv2(igamma) = (npnt*(npnt+1-ipar*2))/2
        nham2 = nv2(igamma)
        call blc2d1(theta(igamma),r,igamma,hr,ham2,nham2,
     1                    term,term2,term3,kz)
        call diag(ham2,nham2,nham2,eig2,work2)
        if (.not. zcut) then
          call choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,
     &                eigtmp,vecs2d,nv2)
          if(igamma .eq. nalf .and. .not. zpmin)
     &        write(6,110) (itmp, iv2(2,itmp),itmp=1,nalf)
        else
          call cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)
          if (.not. zpmin) write(6,110) igamma, iv2(2,igamma)
        endif
30    continue
      call testiv(iv2,nbass)
      nham3 = iv2(1,nalf) + iv2(2,nalf) - 1
      if (.not. zpmin) write(6,120) nham3
      call nftim('end of 2d loop')
      write(6,160)
      call bloc3d(htheta,ham3,eigs1d,vecs1d,iv1,eigs2d,vecs2d,
     &            iv2,nv2,ham1,ham2,nham3,r)
      call nftim('end of 3d Ham building')
      write(6,170)
C
      CALL DIAG3D(HAM3,NHAM3,eig3,WORK3,kz)
      call nftim('end of diagonalising 3d')
      if (ztran) call transr(nv2,iv2,vecs2d,ham3,eig3,r,theta,
     &                      nham3,wvfunc,nbass)
      return
110   format(10x,'For gamma = ',i2,' selected ',i3,' energies.')
120   format(25x,'Total = ',i5)
130   format(5x,'Starting the 2d calculations.')
160   format(5x,'Building the 3d Hamiltonian.')
170   format(5x,'Diagonalising the 3d Hamiltonian.')
      end
 
c***********************************************************************
 
      subroutine blc2d1(xcos,r,igamma,hr,ham2,nham2,
     1                    term,term2,term3,kz)
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      common /split1/ re1,diss1,we1,beta1,ur1,a1,iu1                     
 
      dimension hr(npnt,npnt)
      dimension r(npnt)
      dimension ham2(nham2,nham2)
      data xp5/0.5d0/,x1/1.0d0/

      factr2 = sqrt(xp5)
 
      do 5 j=1,nham2
        do 5 i=1,nham2
          ham2(i,j) = 0.0d0
5     continue

c     Allow for J > 0 case
      if (jrot .ne. 0) then
        FACT =  TERM + TERM2 / (X1 - XCOS)
        IF (KZ .EQ. 1) FACT = FACT + TERM3 * (X1 + XCOS)/(X1-XCOS)
        ia = 0
        do 15 ibeta=1,npnt
          do 25 ialpha=1,ibeta-ipar
            walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
            wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
            wsum = (walpha + wbeta) * fact

            ham2(ialpha+ia,ialpha+ia) = ham2(ialpha+ia,ialpha+ia) + wsum

 25       continue
          ia=ia+ibeta-ipar
          if (ipar .eq. 0) ham2(ia,ia) = ham2(ia,ia) + wsum
 15     continue 
      endif                                                        
                                                                          
 
      q=1.0d0
      if(ipar .eq. 1) q=-1.0d0
      iap=0
      do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ibeta .eq. ibetap)
     &           ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+hr(ialphp,ialpha)
              if(ibeta .eq. ialphp)
     &           ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+q*hr(ibetap,ialpha)
              if(ialpha .eq. ibetap)
     &           ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+q*hr(ialphp,ibeta)
              if(ialpha .eq. ialphp)
     &           ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+hr(ibetap,ibeta)
              if(ialpha .eq. ialphp .and. ibeta .eq. ibetap) then
                call potv(v,r(ialpha),r(ibeta),xcos)
                ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+v
              endif
              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
                call potv(v,r(ibeta),r(ialpha),xcos)
                ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+q*v
              endif
              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=
     &          ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=
     &          ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
      return
      end
 
c***********************************************************************
 
      subroutine choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,
     &                  eigtmp,vecs2d,nv2)
c
c     this routine chooses the max3d lowest eigenvalues from eigs2.
c
      implicit double precision (a-h,o-y), logical (z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
 
      dimension eig2(nham2)
      dimension ham2(nham2,nham2)
      dimension iv2(2,nalf)
      dimension eigs2d(max3d)
      dimension eigtmp(nalf*max2d)
      dimension vecs2d(max2d,max3d)
      dimension nv2(nalf)
      save itotal
      data autocm/2.19474624d+05/
c
      if (igamma .eq. 1) then
         itotal = 0
         rewind intvec
      endif
      do 10 i=1,nham2
        write(intvec) eig2(i)
        write(intvec) (ham2(j,i),j=1,nham2)
10    continue
      iprev=itotal
      inew=nham2
      do 20 ipos = itotal + nham2,1,-1
        if(iprev .ne. 0) then
          if(inew .ne. 0) then
            if(eig2(inew) .lt. eigtmp(iprev)) then
              eigtmp(ipos) = eigtmp(iprev)
              iprev = iprev - 1
            else
              eigtmp(ipos) = eig2(inew)
              inew = inew - 1
            endif
          else
            goto 30
          endif
        else
          eigtmp(ipos) = eig2(inew)
          inew = inew - 1
        endif
20    continue
30    continue
      itotal = itotal + nham2
      if (igamma .eq. nalf) then
        rewind intvec
        if(itotal .le. max3d) then
          emax = eigtmp(itotal)
          write(6,100) itotal,eigtmp(1),emax
        else
          emax = eigtmp(max3d)
          write(6,100) max3d,eigtmp(1),emax
          write(6,110) max3d,eigtmp(1)*autocm,emax*autocm
        endif
        itotal = 1
        do 40 i=1,nalf
          iv2(1,i) = itotal
          ichose=0
          do 50 j=1,nv2(i)
            read(intvec) eig
            if(eig .le. emax) then
              eigs2d(itotal) = eig
              read(intvec) (vecs2d(k,itotal),k=1,nv2(i))
              ichose = ichose + 1
              itotal = itotal + 1
            else
              read(intvec)
            endif
50        continue
          iv2(2,i) = ichose
40      continue
      endif
100   format(10x,'Selecting ',i5,' energies between ',e12.6,
     &       ' and ',e12.6,' in hartrees')
110   format(10x,'Selecting ',i5,' energies between ',e12.6,
     &       ' and ',e12.6,' in wavenumbers')
      return
      end
 
c***********************************************************************
 
      subroutine cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)
c
c     this routine selects all the eigenvalues that are lower than the
c     the cut-off energy emax1, which is user-supplied in wavenumbers.
c     these eigenvalues & their corresponding vectors are then saved
c     in the array vecs1d.
c
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
 
      dimension eig2(nham2)
      dimension ham2(nham2,nham2)
      dimension iv2(2,nalf)
      dimension eigs2d(max3d)
      dimension vecs2d(max2d,max3d)
      save npos
      data autocm /2.19474624D+05/
      if (igamma .eq. 1) npos = 1
 
c     change emax2 to hartree for the selection
      emaxau=emax2/autocm
 
      nvec = 0
      do 10 ialpha = 1,nham2
        if (eig2(ialpha) .lt. emaxau) then
          ntot = npos + nvec
          if (ntot .gt. max3d) then
             write(6,140) max3d, emax2
             stop
          endif
          eigs2d(ntot) = eig2(ialpha)
          do 20 j = 1,nham2
            vecs2d(j,ntot) = ham2(j,ialpha)
20        continue
          nvec = nvec + 1
        endif
10    continue
      iv2(1,igamma) = npos
      iv2(2,igamma) = nvec
      npos = npos + nvec
      return
140   format('Number of 2d eigenvalues greater than ',i4,
     &        ' increase max3d or reduce emax2 which is ',e12.5,
     &        ' cm-1.')
      end
 
      SUBROUTINE TESTIV(IV,nbass)
C
C     Selection vectors for the bisector embedding to ensure that
C     singular region of theta = 0 is not sampled when J > 0.
C     Also calculate which angular grid points are redundant.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      DIMENSION IV(2,nalf)
c
      if (jrot .gt. 0) then
C        First find the extent of the functions in low theta direction
         DO 20 IEND=nalf/2,nalf
         IF (IV(2,IEND) .EQ. 0) GOTO 30
   20    CONTINUE
         IEND=nalf
C        Have we saved functions beyond the point of zero amplitude ?
   30    IOFF=0
         DO 40 I=IEND+1,nalf
C        If so remove them
         IF (IV(2,I) .GT. 0) THEN
            IOFF=IOFF+IV(2,I)
            IV(2,I)=0
            iv(1,i)=iv(1,i-1)
         ENDIF
   40    CONTINUE
         IF (IOFF .GT. 0) WRITE(6,960) IOFF,iv(1,iend)
  960    FORMAT(/5X,'*** Warning:',I4,' Functions removed from theta =',
     1          ' 0 region'/8X,' Basis reset to NHAM3 =',I5)
         IF (IV(2,nalf) .GT. 0) THEN
            IF (zlin) THEN
C              If there are still theta=0 functions, remove them
               WRITE(6,987) IV(2,nalf),iv(1,nalf)-1
               IV(2,nalf)=0
            ELSE
C              Wavefunction has amplitude all the way to theta = 0: STOP
               WRITE(6,950)
               STOP
            ENDIF
         ENDIF
  950    FORMAT(/5X,'BISECTOR EMBEDDING: ',
     1              'WAVEFUNCTION HAS AMPLITUDE FOR THETA = 0. STOP.')
  987    FORMAT(/5X,'*** Warning: ZLIN = T forced the removal of',I4,
     1         ' functions at theta = 0'
     2          /8X,' Basis reset to  NHAM3 =',I5)
       endif
C
      IANG=0
      DO 60 II=1,nalf
   60 IF (IV(2,II) .GT. 0) IANG=IANG+1
      NBASS=IANG*max2d
      IF (ztran) THEN
         WRITE(iwave) IANG,NBASS
         WRITE(iwave) (IV(2,ii),ii=1,nalf)
      ENDIF
      RETURN
      END
 
c***********************************************************************
 
      subroutine bloc3d(htheta,ham3,eigs1d,vecs1d,iv1,eigs2d,
     &                  vecs2d,iv2,nv2,ham1,ham2,nham3,r)
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
 
      dimension htheta(nalf,nalf)
      dimension ham3(nham3,nham3)
      dimension eigs1d(max1dv)
      dimension vecs1d(npnt-ipar,max1dv)
      dimension iv1(2,npnt,nalf)
      dimension eigs2d(max3d)
      dimension vecs2d(max2d,max3d)
      dimension iv2(2,nalf)
      dimension nv2(nalf)
      dimension ham1(npnt-ipar,npnt-ipar)
      dimension ham2(max2d*max2d)
      dimension r(npnt)
C
C Zero ham3
      do 5 j = 1,nham3
        do 5 i = 1,nham3
          ham3(i,j) = 0.0d0
5     continue
      ndim1g = 1
      do 10 igamma = 1,nalf
        if(iv2(2,igamma) .gt. 0) then
          ndim2g = 1
          do 20 igammp = 1,igamma
            if(iv2(2,igammp) .gt. 0) then
              call blc2d2(r,igamma,igammp,htheta,ham2,nv2(igamma))
              call vecmul(vecs2d(1,iv2(1,igamma)),nv2(igamma),
     &                    iv2(2,igamma),vecs2d(1,iv2(1,igammp)),
     &                    nv2(igammp),iv2(2,igammp),max2d,
     &                    ham2,ham3(ndim1g,ndim2g),nham3)
              ndim2g = ndim2g + iv2(2,igammp)
            endif
20        continue
          ndim1g = ndim1g + iv2(2,igamma)
        endif
10    continue
      do 60 i = 1,nham3
        ham3(i,i) = ham3(i,i) + eigs2d(i)
60    continue
      IF (ZPHAM) CALL WRTHAM(HAM3,NHAM3)
C
      return
      end
 
c***********************************************************************
 
      subroutine blc2d2(r,igamma,igammp,htheta,ham2,nham2)
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      common /split1/ re1,diss1,we1,beta1,ur1,a1,iu1
 
      dimension htheta(nalf,nalf)
      dimension r(npnt)
      dimension ham2(nham2,nham2)
      data xp5 /0.5d0/
 
      factr2 = 1.0d0/dsqrt(2.0d0)
 
      do 5 j=1,nham2
        do 5 i=1,nham2
          ham2(i,j) = 0.0d0
5     continue
 
      q=1.0d0
      if(ipar .eq. 1) q=-1.0d0
      iap=0
      do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ialpha .eq. ialphp .and. ibeta .eq. ibetap) then
                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
                wsum = walpha + wbeta
                ham2(ialphp+iap,ialpha+ia)=
     &            ham2(ialphp+iap,ialpha+ia)+htheta(igammp,igamma)*wsum
              endif
              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
                wsum = walpha + wbeta
                ham2(ialphp+iap,ialpha+ia)=
     &           ham2(ialphp+iap,ialpha+ia)+q*htheta(igammp,igamma)*wsum
              endif
              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=
     &          ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=
     &          ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
      return
      end
 
c***********************************************************************
 
      subroutine vecmul(veca,idima,jdima,vecb,idimb,jdimb,nvecln,hama,
     &                                                    hamb,ndim)
 
C This routine does Hb = At Ha B (ummmmm ?????)
 
      implicit double precision(a-h,o-y),logical(z)
      dimension veca(nvecln,jdima)
      dimension vecb(nvecln,jdimb)
      dimension hama(idima,idimb)
      dimension hamb(ndim,jdimb)
 
      do 10 ib = 1,idimb
        do 10 ia = 1,idima
          temp1 = hama(ia,ib)
          if (temp1 .eq. 0.0d0) goto 10
          do 20 ja = 1,jdima
            temp2 = veca(ia,ja) * temp1
            do 30 jb = 1,jdimb
              hamb(ja,jb) = hamb(ja,jb) +  vecb(ib,jb) * temp2
30          continue
20        continue
10    continue
      return
      end
 
c***********************************************************************
 
      subroutine transr(nv2,iv2,vecs2d,ham3,eig3,r,theta,
     &                  nham3,wvfunc,nbass)
      implicit double precision(a-h,o-y),logical(z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      common /split1/ re1,diss1,we1,beta1,ur1,a1,iu1
      common /mass/ xmass(3),g1,g2
 
      dimension nv2(nalf)
      dimension iv2(2,nalf)
      dimension vecs2d(max2d,max3d)
      dimension ham3(nham3,nham3)
      dimension eig3(nham3)
      dimension r(npnt)
      dimension theta(nalf)
      dimension wvfunc(nbass)
 
      call nftim('start of wavefunction generation')
      write(iwave) neval
      call outrow(eig3,neval,iwave)
      do 40 level3=1,neval
        index=1
        do 5 idelta = 1,nbass
          wvfunc(idelta) = 0.0d0
5       continue
        index = 1
        do 30 irs=1,(npnt*(npnt+1-2*ipar))/2
          do 30 igamma=1,nalf
            if (iv2(2,igamma) .eq. 0) goto 30
            do 20 level2=0,iv2(2,igamma)-1
              ivec2=iv2(1,igamma)+level2
              wvfunc(index)=wvfunc(index)+
     &          vecs2d(irs,ivec2)*ham3(ivec2,level3)
20          continue
            index=index+1
30      continue
        write(iwave) wvfunc
40    continue
      call nftim('end of wavefunction generation')
      return
      end
 
c***********************************************************************
 
      SUBROUTINE TRANS(IV1L,IV2L,NDIM2L,VECS1L,
     1                 VECS2L,VECS3L,PHI,PSI,EVALL,NHAM3)
C
C     If ZTRAN then this routine transforms the sets of 1D, 2D and 3D
C     coefficients to PSI, the wavefunction amplitudes at the DVR points
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NPNT1,NPNT2,NALF,NMAX1,NMAX2,MAXLEG,NALF2,IDVR,
     1              NPNT,NLIM1,NLIM2,NEVAL,NCOORD,
     2              JROT,KMIN,IDIA,IPAR,
     3              max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,max1dv,
     3              NDIMA,NDIMB,NDIMC,
     3              EMAX1,EMAX2
      COMMON /OUTP/ ZPHAM,ZPRAD,ZPVEC,ZROT,ZLADD,ZEMBED,ZMORS2,ZLPOT,
     1              ZPMIN,ZVEC,ZQUAD2,ZDIAG,ZLMAT,ZWBLK,zcut,zall,zlin,
     1              ZP1D,ZP2D,ZR2R1,ZTHETA,ZTRAN,zmors1,ztwod,zbisc,
     1              IDIAG1,IDIAG2,IOUT1,IOUT2,iwave,zpfun,ilev,
     2              IEIGS1,IVECS1,IEIGS2,IVECS2,IVINT,IBAND,INTVEC
      COMMON /SPLIT1/ RE1,DISS1,WE1,BETA1,UR1,A1,IU1
      COMMON /SPLIT2/ RE2,DISS2,WE2,BETA2,UR2,A2,IU2
      COMMON /MASS/ XMASS(3),G1,G2
C
C
      DIMENSION VECS1L(MAX2D,NDIMA),VECS2L(MAX2D,MAX2D),VECS3L(MAX3D),
     1          IV1L(NDIMC,NDIMB),NDIM2L(NDIMC),PSI(IDVR,NPNT1,NPNT2),
     1          IV2L(NDIMC),EVALL(NEVAL),PHI(MAX3D,NDIMA,NDIMB)
C
      REWIND IOUT1
      REWIND IOUT2
      REWIND IVECS1
      REWIND IVECS2
C
C     skip header on IOUT2
C
      do 10 i=1,5
   10 read(iout2)
C
C     READ IN IV1 AND IV2 FROM A SEPERATE STREAM (IOUT1)
      CALL IGETRO(IV1L,NDIMB*NDIMC,IOUT1)
      CALL IGETRO(IV2L,NDIMC,IOUT1)
C
C
      DO 32 IONE = 1,NPNTC
      NHAM2 = 0
C     RECALL THE SIZE OF HAM2 FOR EACH NPNTC
      DO 3 ITWO = 1,NPNTB
        NHAM2 = NHAM2 + IV1L(IONE,ITWO)
    3 CONTINUE
   32 NDIM2L(IONE) = NHAM2
C
C     NOW REWRITE THE VECTORS TO DIFFERENT STREAMS
C
      DO 35 IONE = 1,NPNTC
      NHAM2 = NDIM2L(IONE)
      IF (NHAM2 .GT. 0) THEN
         DO 5 KK=1,NPNTA
         CALL GETROW(VECS1L(1,KK),NHAM2,IOUT2)
         IF (IV2L(IONE) .NE. 0) CALL OUTROW(VECS1L(1,KK),NHAM2,IVECS1)
    5    CONTINUE
      ENDIF
   35 CONTINUE
      REWIND IVECS1
C
      DO 31 IONE = 1,NPNTC
      NHAM2 = NDIM2L(IONE)
      IVM = IV2L(IONE)
      IF (NHAM2 .GT. 0) THEN
         DO 30 IND = 1,IVM
         CALL GETROW(VECS2L(1,IND),NHAM2,IOUT2)
         CALL OUTROW(VECS2L(1,IND),NHAM2,IVECS2)
   30    CONTINUE
      ENDIF
   31 CONTINUE
      REWIND IVECS2
C
      READ (IOUT2) MEVAL
      WRITE(IWAVE) MEVAL
C
C     CALL THE EIGENVALUES INTO EVALL
      call getrow(evall,meval,iout2)
      call outrow(evall,meval,iwave)
C
C     NOW READY TO START READING THE 3-D VECTORS
C
C     1-D VECTORS:  IVECS1 (Stream 26)
C     2-D VECTORS:  IVECS2 (Stream 27)
C     3-D VECTORS:  IOUT2  (Stream 25)
C
      IND2 = 0
      DO 58 IC = 1,NPNTC
      NHAM2 = NDIM2L(IC)
C
      DO 56 J  = 1,IV2L(IC)
      IND2 = IND2 + 1
      IF (NHAM2 .GT. 0) CALL GETROW(VECS2L(1,J),NHAM2,IVECS2)
C
      DO 54 IA = 1,NPNTA
      IF (NHAM2.GT.0 .AND. J.EQ.1)
     1 CALL GETROW(VECS1L(1,IA),NHAM2,IVECS1)
      IND1 = 0
      DO 52 IB = 1,NPNTB
      SUM1 = 0.0D0
      DO 50 K  = 1,IV1L(IC,IB)
      IND1 = IND1 + 1
      SUM1 = SUM1 + VECS1L(IND1,IA)*VECS2L(IND1,J)
C
   50 CONTINUE
      PHI(IND2,IA,IB) = SUM1
   52 CONTINUE
C
   54 CONTINUE
C
   56 CONTINUE
C
   58 CONTINUE
C
C
C
      DO 959 LL=1,MEVAL
C
      CALL GETROW(VECS3L,NHAM3,IOUT2)
C
C     SET PSI TO ZERO
      DO 160 IA = 1,IDVR
      DO 160 IB = 1,NPNT1
      DO 160 IC = 1,NPNT2
  160 PSI(IA,IB,IC) = 0.0D0
C
C
C
      IF (ZTHETA) THEN
C
        if (.not. zr2r1) then
C
          IND2 = 0
          DO 158 IC = 1,NPNT2
C
          DO 156 J  = 1,IV2L(IC)
          IND2 = IND2 + 1
C
          DO 154 IA = 1,IDVR
          DO 152 IB = 1,NPNT1
          PSI(IA,IB,IC) = PSI(IA,IB,IC) + VECS3L(IND2)*PHI(IND2,IA,IB)
  152     CONTINUE
C
  154     CONTINUE
C
  156     CONTINUE
C
  158     CONTINUE
C
        else
C
          IND2 = 0
          DO 258 IC = 1,NPNT1
C
          DO 256 J  = 1,IV2L(IC)
          IND2 = IND2 + 1
C
          DO 254 IA = 1,IDVR
          DO 252 IB = 1,NPNT2
          PSI(IA,IC,IB) = PSI(IA,IC,IB) + VECS3L(IND2)*PHI(IND2,IA,IB)
  252     CONTINUE
C
  254     CONTINUE
C
  256     CONTINUE
C
  258     CONTINUE
C
        endif
C
      ELSE
C
        if (.not. zr2r1) then
C
          IND2 = 0
          DO 358 IC = 1,IDVR
C
          DO 356 J  = 1,IV2L(IC)
          IND2 = IND2 + 1
C
          DO 354 IA = 1,NPNT1
          DO 352 IB = 1,NPNT2
          PSI(IC,IA,IB) = PSI(IC,IA,IB) + VECS3L(IND2)*PHI(IND2,IA,IB)
  352     CONTINUE
C
  354     CONTINUE
C
  356     CONTINUE
C
  358     CONTINUE
C
        else
C
          IND2 = 0
          DO 458 IC = 1,IDVR
C
          DO 456 J  = 1,IV2L(IC)
          IND2 = IND2 + 1
C
          DO 454 IA = 1,NPNT2
          DO 452 IB = 1,NPNT1
          PSI(IC,IB,IA) = PSI(IC,IB,IA) + VECS3L(IND2)*PHI(IND2,IA,IB)
  452     CONTINUE
C
  454     CONTINUE
C
  456     CONTINUE
C
  458     CONTINUE
C
        endif
C
      ENDIF
C
C     now save the wavefunction for each level (to stream 28)
      CALL OUTROW(PSI,IDVR*NPNT1*NPNT2,IWAVE)
C
  959 CONTINUE
      RETURN
      END
      SUBROUTINE POTVLG(V,R1,R2,XCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
C
C     NO. OF TERMS IN LEGENDRE EXPANSION MUST BE SET IN THE NEXT 2 LINES
      DIMENSION PLEGV(9),VL(9)
      CALL LEGV(PLEGV,9,XCOS)
      CALL POT(V0,VL,R1,R2)
      V=V0
      DO 20 L=1,9
      V = V + VL(L) * PLEGV(L)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE LEGV(PLEGG,NLEG,X)
C     CALCULATE LEGENDRE POLYNOMIALS 1 TO NLEG AT X                 #029
C     THE POLYNOMIALS ARE NOT NORMALISED.
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      DIMENSION PLEGG(NLEG)
      PLEGG(1) = X
      IF (NLEG .LE. 1) RETURN
      PLEGG(2) = 1.5D0*X*X - 0.5D0
      DO 10 L=3,NLEG
      L1 = L - 1
      PLEGG(L) = (DBLE(L+L1)*X*PLEGG(L1)-DBLE(L1)*PLEGG(L-2))/DBLE(L)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SQOUT(SQMAT,NDIM)
C     PRINT LOWER TRIANGLE OF SQUARE MATRIX
      DOUBLE PRECISION SQMAT(NDIM,NDIM)
      DO 30 I=1,NDIM
      WRITE(6,1020) (SQMAT(I,J),J=1,I)
 1020 FORMAT(10F13.7)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE SYMOUT(SYMMAT,NDIM)
C
C     PRINT OUT LOWER TRIANGLE OF SYMMETRIC MATRICES                #008
C
      DOUBLE PRECISION SYMMAT(1)
      IP=0
    3 LLOW=10*IP+1
      LUP=MIN(LLOW+9,NDIM)
      WRITE(6,4) (I,I=LLOW,LUP)
    4 FORMAT(/,I11,9I13)
      IND0=LLOW*(LLOW+1)/2
      DO 5 I=LLOW,NDIM
      ITOP=MIN(LUP,I)
      KTOP=ITOP-LLOW+IND0
      WRITE(6,7) I,(SYMMAT(J),J=IND0,KTOP)
    7 FORMAT(I4,F12.7,9F13.7)
    5 IND0=IND0+I
      IF(LUP.GE.NDIM) RETURN
      IP=IP+1
      GO TO 3
      END
      SUBROUTINE WRTHAM(HAMIL,NHAM)
C     PRINT HAMILTONIAN MATRIX                                      #011
      DOUBLE PRECISION HAMIL(NHAM,NHAM)
      WRITE(6,1010)
 1010 FORMAT(1H1,5X,'HAMILTONIAN MATRIX'/)
      DO 30 I=1,NHAM
      WRITE(6,1020) (HAMIL(I,J),J=1,I)
 1020 FORMAT(10F13.7)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE GETROW(ROW,NROW,IUNIT)
C     FAST NON-FORMATTED READ
      DOUBLE PRECISION ROW(NROW)
      READ(IUNIT) ROW
      RETURN
      END
      SUBROUTINE IGETRO(IVEC,NSIZE,ISTREAM)
      DIMENSION IVEC(NSIZE)
      READ(ISTREAM) IVEC
      RETURN
      END
      SUBROUTINE OUTROW(ROW,NROW,IUNIT)
C     FAST NON-FORMATTED WRITE                                      #025
      DOUBLE PRECISION ROW(NROW)
      WRITE(IUNIT) ROW
      RETURN
      END
      SUBROUTINE IOUTRO(IVEC,NSIZE,ISTREAM)
      DIMENSION IVEC(NSIZE)
      WRITE(ISTREAM) IVEC
      RETURN
      END
c***********************************************************************
 
      subroutine nftim(text)
      character text*(*)
      write(6,10)
      write(6,*) 'Time at ',text,' is.........'
      call timer
      return
10    format(/)
      end
      SUBROUTINE TIMER
C     PRINTS CURRENT CPU TIME USAGE                                 #030
c     Needs a subroutine which can access the machine clock
      real array(2)
      DATA TIME/0.0/
c     time=etime(array)
      CALL SECOND(TIME)
      MIN = INT(TIME) / 60
      SEC = TIME - REAL(60*MIN)
      WRITE(6,1) MIN,SEC
    1 FORMAT(/I12,' MINS',F6.2,' SECS CPU TIME USED'/)
      RETURN
      END
      SUBROUTINE POT(V0,VL,R1,R)
C
C     CALCULATES ENERGY
C     ESHORT = SIGMA(L) V(L)*P(L)                          L=0,6
C                       V(L)=EXP(-A(L)-B(L)*R-C(L)*R*R)
C     ELONG1 = SIGMA(L) D(L)/R**(L+1)*P(L)                 L=0,6
C     ELONG2 = SIGMA(L) C(L,I)/R**(L+2)*P(I)      L=L1+L2, L1,L2=0,3
C                                                    I=|L1-L2|,L1-L2
C     V & R BOTH IN ATOMIC UNITS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VL(9),RINV(8),A(10),B(10),C(10),D(7)
C
      DATA DAMP/1.0D0/,INAME/0/,LDAMP/2/
      DATA A / 0.1383211600D 01, 0.2957913200D 01, 0.4742029700D 01,
     A         0.1888529900D 01, 0.4414335400D 01, 0.4025649600D 01,
     A         0.5842589900D 01, 0.2616811400D 01, 0.6344657900D 01,
     A        -0.1520228000D 02/
      DATA B /-0.1400070600D 00,-0.1479771600D 01,-0.1811986200D 01,
     B        -0.1287503000D 01,-0.2322971400D 01,-0.2775383200D 01,
     B        -0.3480852900D 01,-0.2655595200D 01,-0.4344980100D 01,
     B         0.6549253700D 01/
      DATA C /-0.2078921600D 00, 0.1161308200D-01, 0.1718068000D-01,
     C        -0.2772849100D-01, 0.7069274200D-01, 0.1377197800D 00,
     C         0.1863111400D 00, 0.5881550400D-02, 0.1529136800D 00,
     C        -0.1302567800D 01/
      DATA D /-0.1         D 01,-0.215135    D 00,-0.3414573000D 01,
     D        -0.3818815000D 01,-0.1584152000D 02,-0.1429374000D 02,
     D        -0.4381719000D 02/
      DATA
     4 C40,C42/-1.05271D1,-3.17   D0/,
     5 C51,C53/-2.062328D1,+3.7320800D0/,
     6 C60,C62,C64/-5.749396D1,-1.068192D2, 1.714139D1/,
     7 C71,C73,C75/-2.028972D2,-7.523207D1,-2.845514D1/,
     8 C80,C82,C84,C86/-4.582015D2,-3.537347D2,-1.126427D2,-1.082786D2/
      DATA DEMPA /-1.515625D0/, DEMPR0 /1.900781D0/
C
C
      DATA ONE/1.0D0/,RSHORT/9.0D0/
C
      IF (INAME .NE. 0) GOTO 1
    1 INAME=1
      R2=R*R
C     LONG RANGE CONTRIBUTIONS SHOULD BE TEMPERED
      DR = R - DEMPR0
      DAMPOT = ONE - DEXP(DEMPA*DR*DR)
      RINV(1)=DAMPOT/R
      DO 30 L=1,7
      RINV(L+1)=RINV(L)/R
   30 CONTINUE
C
C     CALCULATE ENERGY
C
      V0    = C40*RINV(4) + C60*RINV(6) + C80*RINV(8) + RINV(1) * D(1)
      VL(1) = C51*RINV(5) + C71*RINV(7)               + RINV(2) * D(2)
      VL(2) = C42*RINV(4) + C62*RINV(6) + C82*RINV(8) + RINV(3) * D(3)
      VL(3) = C53*RINV(5) + C73*RINV(7)               + RINV(4) * D(4)
      VL(4) =               C64*RINV(6) + C84*RINV(8) + RINV(5) * D(5)
      VL(5) =               C75*RINV(7)               + RINV(6) * D(6)
      VL(6) =                             C86*RINV(8) + RINV(7) * D(7)
C     IF R IN REGION WHERE NO SHORT RANGE CONTRIBUTION, RETURN
      IF (R .GT. RSHORT) GOTO 60
      V0   = V0   + DEXP( A(1) + B(1)*R + C(1)*R2 )
      DO 40 L=2,7
      VL(L-1) = VL(L-1) + DEXP( A(L) + B(L)*R + C(L)*R2 )
   40 CONTINUE
      DO 50 L=8,10
      VL(L-1) =           DEXP( A(L) + B(L)*R + C(L)*R2 )
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
