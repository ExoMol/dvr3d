!     DUMMY MAIN PROGRAM                                   

      CALL ROTLEV
      STOP
      END 

!######################################################################
      SUBROUTINE ROTLEV

!     PROGRAM               R O T L E V 3
!
!     PROGRAM TO DO ROTATIONAL ANALYSIS USING VIBRATIONAL OUTPUT OF
!     DVR3DRJZ RUN WITH ZROT = .TRUE.
!
!     SEE:
!     USE AS FOLLOWS:
!     COMMENTS ON NAMELIST PARAMETERS (& DEFAULTS) IN BLOCK DATA
!     THE PROGRAM NEEDS THE FOLLOWING SUBROUTINES:
!     1. LIMITED CARD INPUT WHICH IS READ IN SUBROUTINE INSIZE
!        AND AN INPUT FILE ON STREAM iwave FROM DVR3D
!     2. F02FJF TO DO ITERATIVE DIAGONALISATION (NAG ROUTINE)
!     or dsyev to do in core diagonalisation (lapack f90 routine).
!     THE PROGRAM WORKS IN **** ATOMIC UNITS ***** :
!     *** VERSION TO PERFORM BASIS SELECTION ON ENERGY,
!         TO LOOP OVER E AND F BLOCKS,
!         USES STORAGE SAVING HAMILTONIAN FILE ISCR,
!         TRANSFORMS VECTORS TO ORIGINAL BASIS.
!     Fortan90 version with dynamic arrays by Max Kostin & Jonathan Tennyson


      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,&
                    zpseg
      NAMELIST/PRT/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,&
                    zpham,zpvec,zdcore,zvec,zpfun,zdiag,ztran,zptra,&
                    zpseg
   
      WRITE(6,1000)
 1000 FORMAT(5x,'PROGRAM ROTLEV3 (VERSION OF March 2002):',/)

!     READ IN NAMELIST INPUT DATA (DEFAULTS IN BLOCK DATA)
      READ(5,PRT)
      OPEN(UNIT=ISCR,FORM='UNFORMATTED')
      call SYSTEM_CLOCK(itime0,irate2,imax2)

!     READ IN CONTROL PARAMETERS OF PROBLEM.
      CALL INSIZE

      CALL select 

      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,'secs CPU time used'/) 

      STOP
      END

!##############################################################################
      BLOCK DATA
!     STORES DEFAULTS FOR NAMELIST PARAMETERS  
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

!     OUTP HOLDS INFORMATION WHICH CONTROLS THE AMOUNT OF PRINTED OUTPUT
!     TOLER: CONVERGENCE TOLERANCE FOR THE ITERATIVE DIAGONALISER
!            TOLER = 0.0 GIVES MACHINE ACCURACY
!     ZPHAM: PRINT MATRIX HAMIL IF ZPHAM = .TRUE.
!     ZPVEC: PRINT EIGENVECTORS IF ZPVEC = .TRUE.
!     THRESH: THRESHOLD FOR PRINTING A COEFFICIENT IF ZPVEC=.TRUE.
!     ZPTRA: PRINT TRANSFORMED VECTORS IF ZPTRA = .TRUE.
!     zdcore: = .true. for in core diagonalisation
!     zdiag:  = .false. do not diagonalise the Hamiltonian matrix.
!     ztran:  = .true. transform eigenvectors back to original basis.
!     zvec:   = .true. eigenvalues and eigenvectors written to disk file.
!     zpfun:  = .true.  eigenvalues concatenated on stream ILEV.
!     ZEMBED: = .TRUE. Z AXIS IS ALONG R2, = .FALSE. Z AXIS IS ALONG R1
!     z1da: = .true., input from DVR1DA, .false. from DVR3DRJ
!     STREAM         HOLDS                              USED IF
!      ILEV    INPUT/OUTPUT OF EIGENVALUES              ZPFUN=.TRUE.
!      iwave   INPUT  EIGENVALUES & EIGENVECTORS        ALWAYS
!      jscr    transformed input vector scratch file    ZTRAN=.TRUE.
!      JVEC    OUTPUT FIRST  SET EIGENVALUES & VECTORS  ZVEC=.TRUE.
!      JVEC2   OUTPUT SECOND SET                        ZVEC=.TRUE.
!      KVEC    OUTPUT FIRST  SET TRANSFORMED VECTORS    ZTRAN=.TRUE.
!      KVEC2   OUTPUT SECOND SET                        ZTRAN=.TRUE.
!      ISCR    HAMILTONIAN FILE                          ALWAYS
!      IRES = 0 NORMAL RUN
!          = 1 RESTART FROM FIRST  CALL TO DGROT
!          = 2 RESTART FROM SECOND CALL TO DGROT
!          < 0 TRANSFORMATION ONLY 
!      IRES  = -1 with kmin=0 or 1 transform one set of vectors  
!                 with kmin=2 transform both sets of vectors            
!      IRES  = -2 with kmin=2  transform 2nd set of vectors
!      RESTART AFTER ZDIAG=FALSE RUN, IWAVE=IRF1 and IRF2 REQUIRED

      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg
      DATA TOLER/0.0D0/,THRESH/0.1D0/,ZPHAM/.FALSE./,ZPVEC/.FALSE./,&
           iwave/26/,ZVEC/.FALSE./,JVEC/3/,JVEC2/2/,ISCR/10/,IRES/0/,&
           ZTRAN/.FALSE./,KVEC/8/,KVEC2/9/,ZPTRA/.FALSE./,jscr/7/,&
           ZPFUN/.FALSE./,ILEV/14/,ZDIAG/.TRUE./,zdcore/.false./,&
           z1da/.false./,IRF1/21/,IRF2/22/,zpseg/.false./
      END

!####################################################################################
      SUBROUTINE INSIZE

!     SET UP COMMON /SIZE/ & WRITE CONTROL PARAMETERS OF PROBLEM 

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

!     COMMON /SIZE/ STORES CONTROL PARAMETERS FOR THE PROBLEM
!     NBASS: MAXIMUM DIMENSION OF ROTATIONAL SECULAR PROBLEM
!     IBASS: ACTUAL  DIMENSION OF ROTATIONAL SECULAR PROBLEM
!     MBASS: MAXIMUM SIZE OF VIBRATIONAL PROBLEM
!     KBASS: SIZE OF ALL VIBRATIONAL PROBLEM SUMMED
!     JROT : TOTAL ANGULAR MOMENTUM OF THE COMPLEX
!     IDIA : =1 HETRONUCLEAR DIATOMIC, =2 HOMONUCLEAR DIATOMIC
!          : =0 RADAU COORDINATES,
!     IPAR : PARITY OF BASIS - IF IDIA=2: IPAR=1 FOR EVEN & =0 FOR ODD
!          : TAKEN ACCORDING TO PARITY OF LMAX
!     KMIN : KMIN=1 FOR SYM. ROTATIONAL BASIS, =0 FOR ANTI-SYM.
!            KMIN>1 LOOP OVER BOTH.
!     NLIM : =npnt1  (ZEMBED=.FALSE.)
!            =npnt2*npnt2  (ZEMBED=T and ZQUAD2=F)
!            =npnt2  (ZEMBED=T and ZQUAD2=T)
!     MEVAL: NUMBER OF EIGENVALUES COMPUTED IN THE VIBRATIONAL PROBLEM
!     NVIB : NUMBER OF VIBRATIONAL EIGENVALUES USED IN ROTATIONAL PROB.
!     NEVAL: NUMBER OF EIGENVALUES WHICH HAVE TO ACTUALLY BE COMPUTED
!     NEVAL2: NEVAL FOR F BLOCK WHEN KMIN>1.
!     KEVAL: NUMBER OF EIGENVECTORS USED FOR ITERATIONS (=NEVAL+4)
!     NBLK : NUMBER OF BLOCKS IN K
!     LOFF : SPACE REQUIRED FOR ALL THE OFF-DIAGONAL BLOCKS
!     LOFF0: SPACE REQUIRED FOR THE LARGEST OFF-DIAGONAL BLOCK
!     npnt1: number of r1 dvr points
!     npnt2: number of r2 dvr points
!     npntt: (maximum) number of theta dvr points

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg
      CHARACTER(len=8) TITLE(9)
      DOUBLE PRECISION, DIMENSION(3) :: xmass
      DATA X0/0.0D0/

      if (zpseg==.true.) then 
      OPEN(UNIT=iwave,FORM='UNFORMATTED',recordtype='segmented')
      open(unit=IRF2,form='unformatted',recordtype='segmented')
      else 
      OPEN(UNIT=iwave,FORM='UNFORMATTED')
      open(unit=IRF2,form='unformatted')     
      end if 


      READ(iwave) IDIA,IPAR,npntt,npnt1,npnt2,JROT,KMIN0,MEVAL,nlim
      read(IWAVE) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,z1da,ZQUAD2

      IF(.NOT. ZDIAG) THEN
         if (zpseg==.true.) then 
            OPEN(UNIT=IRF1,FORM='UNFORMATTED',recordtype='segmented')
         else 
            OPEN(UNIT=IRF1,FORM='UNFORMATTED')
         end if    
         WRITE(IRF1) IDIA,IPAR,npntt,npnt1,npnt2,JROT,KMIN0,MEVAL,nlim
         WRITE(IRF1) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,z1da,ZQUAD2
      ENDIF

      READ(5,5)  NVIB,NEVAL,KMIN,IBASS,NEVAL2
    5 FORMAT(12I5)

!     COMPUTE SIZE OF ROTATIONAL SECULAR PROBLEM

      mbass=npnt1*npnt2*npntt
      NVIB=MIN(NVIB,MEVAL)
      NBLK=JROT
      IF (KMIN .GT. 0) NBLK=NBLK+1
      NBASS=NBLK*NVIB
      IF (IBASS .GT. 0) NBASS=MIN(NBASS,IBASS)
      IF (NEVAL .LE. 0) NEVAL = 10
      NEVAL=MIN(NEVAL,NBASS)
      IF (NEVAL2 .LE. 0) NEVAL2 = NEVAL

      WRITE(6,1000) MEVAL,MBASS,NVIB,NEVAL,NBASS
 1000 FORMAT(/5X,'ROTATIONAL PART OF ROT-VIB CALCULATION  WITH:',&
             /I9,3X,'LOWEST VIBRATIONAL EIGENVECTORS SUPPLIED FROM',&
             /I9,3X,'DIMENSION VIBRATION SECULAR PROBLEM',&
             /I9,3X,'LOWEST VIBRATIONAL EIGENVECTORS ACTUALLY USED',&
             /I9,3X,'LOWEST ROTATIONAL EIGENVECTORS REQUIRED FOR',&
             /I9,3X,'DIMENSION ROTATION SECULAR PROBLEM')
      IF (IBASS .GT. 0) WRITE(6,1005)
 1005 FORMAT(17X,'WITH BASIS SELECTED BY ENERGY ORDERING')

      READ(5,500)   TITLE
  500 FORMAT(9A8)
      WRITE(6,1010) TITLE
 1010 FORMAT(/5X,'TITLE: ',9A8)
      IF (.not.ZDIAG) THEN
         WRITE(6,1012) ISCR,IRF1,IRF2
 1012    FORMAT(/5X,'HAMILTONIAN CONSTRUCTION ONLY REQUESTED',&
                /5X,'TO BE WRITTEN TO STREAM ISCR =',I4,&
                /5X,'RESTART INFORMATION WRITTEN TO STREAMS IRF1=',I4,&
                ' AND IRF2=',I4)
         ZPHAM=.FALSE.
         IRES=0
      else
         if (zdcore) then
            write(6,1013)
 1013 format(/5x,'Diagonalisation performed in core using DSYEV')
         else
            write(6,1014)
 1014 format(/5x,'Diagonalisation performed iteratively using F02FJF')
            IF (TOLER .NE. X0) WRITE(6,1060) TOLER
 1060 FORMAT(5X,'EIGENVALUE CONVERGENCE TOLERANCE, TOLER =',D20.3)
            IF (TOLER .EQ. X0) WRITE(6,1070)
 1070 FORMAT(5X,'EIGENVALUES CONVERGED TO MACHINE ACCURACY')
         endif
      ENDIF
      IF (IRES .NE. 0) WRITE(6,1015) IRES
 1015 FORMAT(/5X,'***** RESTART RUN, IRES =',I2,' *****')
      IF ( IRES .EQ. -1)WRITE(6,1016)
 1016 FORMAT(/5X,'***** TRANSFORMATION ONLY *****')
      IF (ZPHAM) WRITE(6,1020)
 1020 FORMAT(/5X,'PRINTING OF HAMILTONIAN MATRIX REQUESTED')
      IF (.NOT.ZPHAM) WRITE(6,1030)
 1030 FORMAT(/5X,'PRINTING OF HAMILTONIAN MATRIX NOT REQUESTED')
      IF (ZPVEC) WRITE(6,1040) THRESH
 1040 FORMAT(5X,'PRINTING OF EIGENVECTOR COEFFICIENTS GREATER THAN',&
                 ' THRESH =',F5.2,' REQUESTED')
      IF (.NOT.ZPVEC) WRITE(6,1050)
 1050 FORMAT(5X,'PRINTING OF EIGENVECTORS NOT REQUESTED')
      IF (ZTRAN) THEN
          ZVEC = .TRUE.
          IF (ZPTRA) WRITE(6,1052)
 1052     FORMAT(5X,'PRINTING OF TRANSFORMED VECTORS REQUESTED')
          IF (.NOT.ZPTRA) WRITE(6,1053)
 1053     FORMAT(5X,'PRINTING OF TRANSFORMED VECTORS NOT REQUESTED')
      ENDIF
      IF (.not. z1da)  WRITE(6,1051) iwave
 1051 FORMAT(/5X,'DVR3D   DATA   TO BE READ         FROM STREAM ',&
             'IWAVE =',I4)
      IF (z1da)  WRITE(6,1151) iwave
 1151 FORMAT(/5X,'DVR1DA  DATA   TO BE READ         FROM STREAM ',&
             'IWAVE =',I4)
      IF (ZPFUN) WRITE(6,1055) ILEV
 1055 FORMAT(/5X,'EIGENVALUES    TO BE WRITTEN TO END OF STREAM ',&
             'ILEV  =',I4)
      IF (ZVEC) WRITE(6,1054) JVEC
 1054 FORMAT(/5X,'EIGENVALUES & VECTORS TO BE WRITTEN TO STREAM ',&
             'JVEC  =',I4)
      IF (ZVEC .AND. KMIN .GT. 1) WRITE(6,1056) JVEC2
 1056 FORMAT( 5X,'SECOND SET            TO BE WRITTEN TO STREAM ',&
             'JVEC2 =',I4)
      IF (ZTRAN) WRITE(6,1057) jscr
 1057 FORMAT(/5X,'Transformed input vectors stored on    stream ',&
             'JSCR  =',I4)
      IF (ZTRAN) WRITE(6,1058) KVEC
 1058 FORMAT(/5X,'TRANSFORMED VECTORS   TO BE WRITTEN TO STREAM ',&
             'KVEC  =',I4)
      IF (ZTRAN .AND. KMIN .GT. 1) WRITE(6,1059) KVEC2
 1059 FORMAT( 5X,'SECOND SET            TO BE WRITTEN TO STREAM ',&
             'KVEC2 =',I4)
      IF (TOLER .NE. X0) WRITE(6,1066) TOLER
 1066 FORMAT(/5X,'EIGENVALUE CONVERGENCE TOLERANCE, TOLER =',D20.3)
      IF (TOLER .EQ. X0) WRITE(6,1077)
 1077 FORMAT(/5X,'EIGENVALUE CONVERGENCE TO MACHINE ACCURACY')
      IF (IDIA .NE. 2) THEN
         IF (IDIA .EQ. 1) WRITE(6,1080)
 1080    FORMAT(/5X,'DIATOMIC ASSUMED HETRONUCLEAR')
         IF (IDIA .EQ. 0) WRITE(6,1085)
 1085    FORMAT(/5X,'RADAU COORDINATES USED')
      ELSE
         IF (IPAR .EQ. 0) WRITE(6,1100)
 1100    FORMAT(/5X,'DIATOMIC ASSUMED HOMONUCLEAR',&
                /5X,'EVEN PARITY FUNCTIONS IN BASIS SET')
         IF (IPAR .EQ. 1) WRITE(6,1110)
 1110    FORMAT(/5X,'DIATOMIC ASSUMED HOMONUCLEAR',&
                /5X,'ODD PARITY FUNCTIONS IN BASIS SET')
      ENDIF
      WRITE(6,1120) JROT
 1120 FORMAT(/5X,'J =',I3,' ROTATIONAL STATE')
      IF (KMIN .EQ. 0) WRITE(6,1130)
 1130 FORMAT(12X,'WITH ANTI-SYMMETRIC |JK> - |J-K> FUNCTIONS IN BASIS')
      IF (KMIN .EQ. 1) WRITE(6,1140)
 1140 FORMAT(12X,'WITH SYMMETRIC |JK> + |J-K> FUNCTIONS IN BASIS')
      IF (KMIN .GT. 1) WRITE(6,1150)
 1150 FORMAT(12X,'LOOP OVER SYMMETRIC AND ANTI-SYMMETRIC FUNCTIONS')
      IF (KMIN .GT. 0 .AND. KMIN0 .NE. 1) GOTO 960
      IF (ZEMBED) then
         WRITE(6,1160) 2
      else
         WRITE(6,1160) 1
         zquad2=.true.
      endif
 1160 FORMAT(/5X,'Z AXIS EMBEDED ALONG THE R',I1,' COORDINATE')
      RETURN
  960 WRITE(6,970) KMIN,KMIN0
  970 FORMAT(//6X,'**** KMIN =',I3,' BUT KMIN0 =',I3,' STOP ****',/)
      STOP
      END

!#################################################################################
      SUBROUTINE SELECT

!     SUBROUTINE SELECT DETERMINES WHICH VIBRATIONAL BASIS   
!     FUNCTIONS ARE TO BE USED

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                   KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                   kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: EVIB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IV,MVIB

      ALLOCATE(evib(nvib,nblk),iv(nblk),mvib(nblk))

      if (ires .eq. 0)then

!     READ ENERGIES FROM FILE iwave, FIRST SKIP MATRIX ELEMENTS etc
       IF (.NOT. ZDIAG) THEN

       read(IWAVE) RE1,DISS1,WE1,RE2,DISS2,WE2
       write(IRF1) RE1,DISS1,WE1,RE2,DISS2,WE2
!      transfer the DVR points
       read(iwave) dummy
       write(IRF1) dummy
       call getrow(evib,npnt1,iwave)
       call outrow(evib,npnt1,irf1)
       call getrow(evib,npnt2,iwave)
       call outrow(evib,npnt2,irf1)

       ELSE

         do 10 i=1,4
         READ(iwave,END=900)
   10    continue
      ENDIF

!     READ IN NEXT SET OF VIBRATIONAL VECTORS & SET DIAGONAL ELEMENTS
      IPT=0
      KBASS=0
      DO 100 K=1,NBLK
  110 READ(iwave,end=900)
      READ(iwave,end=900) K2,maxleg,idvr,lincr
!     SKIP BACK IF K2=0 AND WE ARE DOING AN F PARITY CALCULATION
      IF (K2 .EQ. 0 .AND. KMIN .EQ. 0) THEN
          READ(iwave,END=900)
          READ(iwave,END=900) meval
          do 120 i=1,meval+1
          READ(iwave,END=900)
  120     continue
          GOTO 110
      ENDIF

      ibass2=idvr*npnt1*npnt2
      read(iwave,end=900)
      read(iwave,end=900) meval
!      write(6,*)K,nvib,meval,iwave
      MVIB(K)=MIN(NVIB,meval)
      CALL GETROW(EVIB(1,K),MVIB(K),iwave)
      do 130 i=1,meval
      READ(iwave,END=900)
  130 continue
      IPT=IPT+MVIB(K)
      KBASS=KBASS+IBASS2
  100 CONTINUE
!     REPOSITION FILE iwave
      REWIND iwave
      READ(iwave)
      READ(iwave)
      READ(iwave)
      IF (IBASS .LE. 0 .OR. IBASS .GE. IPT) THEN
          IBASS=IPT
          IVIB=NVIB
!          write(6,*)'ci sono 300  ',mvib
          GOTO 300
      ENDIF

!     SELECT THE IBASS LOWEST FUNCTIONS

      EMIN=EVIB(1,1)
      DO 160 I=1,NBLK
      EMIN=MIN(EMIN,EVIB(1,I))
      IV(I)=1
  160 continue
      IPT=1
      DO 200 N=1,IBASS
  210 IF (IV(IPT) .LE. MVIB(IPT)) THEN
          EVIBR=EVIB(IV(IPT),IPT)
          JPT=IPT
      ELSE
          IPT=IPT+1
          GOTO 210
      ENDIF
      DO 220 J=IPT+1,NBLK
      IF (IV(J) .GT. MVIB(J)) GOTO 220
      IF (EVIB(IV(J),J) .GE. EVIBR) GOTO 220
      EVIBR=EVIB(IV(J),J)
      JPT=J
  220 CONTINUE
!     KEEP THE BASIS FUNCTION
      IV(JPT)=IV(JPT)+1
  200 CONTINUE
!     STORE THE NUMBER OF FUNCTIONS SELECTED FOR EACH K
      MVIB(1)=IV(1)-1
      IVIB=MVIB(1)

      DO 230 I=2,NBLK
      MVIB(I)=IV(I)-1
      IVIB=MAX(IVIB,MVIB(I))
  230 continue
      WRITE(6,1000) NBASS,EMIN,EVIBR
 1000 FORMAT(/,I9,'FUNCTIONS SELECTED FROM E =',D20.10,' TO',D20.10)
      else  !if ires.ne.0
         read(IRF2) nbass,neval,neval2,ipar,idia,jrot,kmin,nblk
         read(IRF2) mvib
         if(ires .le. -1 .OR. ZTRAN)then
            read(IRF2) kbass
         else
            kbass=1
         endif
         ivib=mvib(1)
         DO 232 I=2,NBLK
           IVIB=MAX(IVIB,MVIB(I))
  232    continue
      endif

!     DETERMINE AND PRINT BASIS SET LABELS

  300 WRITE(6,1010)
 1010 FORMAT(//5X,' BASIS FUNCTIONS SELECTED')
!      write(6,*)'ci sono !'
      IF (KMIN .EQ. 0) THEN
         KZ=0
      ELSE
         KZ=-1
      ENDIF
      ipu=0
      LOFF=0
      LOFF0=0
      DO 310 j=1,NBLK
      IF (J .GT. 1) THEN
         LENG=MVIB(J-1)*MVIB(J)
         LOFF=LOFF+LENG
         LOFF0=MAX(LOFF0,LENG)
      ENDIF
      KZ=KZ+1
      IPD=IPU
      IPU=IPD+MVIB(j)
      if (mvib(j).eq.nvib) then
      WRITE(6,1021) KZ,IPD+1,IPU,1.d2*dble(mvib(j))/dble(nvib)
      else
      WRITE(6,1020) KZ,IPD+1,IPU,1.d2*dble(mvib(j))/dble(nvib)
      end if
  310 CONTINUE
 1020 FORMAT(5X,'K =',I3,', I RUNS FROM',I5,' TO',I5,' ',f6.2,' %')
 1021 FORMAT(5X,'K =',I3,', I RUNS FROM',I5,' TO',I5,' ',f6.2,' % *saturated*')

      DEALLOCATE(evib,iv)

      NBASS=IBASS
      NVIB=IVIB
      call vrmain(mvib,maxleg)

      RETURN
!     ERROR ON UNIT iwave
  900 REWIND iwave
      I=1
  920 READ(iwave,END=930)
      I=I+1
      GOTO 920
  930 WRITE(6,940) iwave,I
  940 FORMAT(/,'   END OF STREAM iwave =',I4,' AT RECORD',I6)
      STOP
      END


!#########################################################################
      SUBROUTINE VRMAIN(MVIB,maxleg)

!     SUBROUTINE VRMAIN IS THE 'REAL' MAIN PROGRAM & CONTAINS     
!     THE CALLS TO THE VARIOUS SUBROUTINES WHICH SET & SOLVE HAMIL

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DIMENSION MVIB(NBLK)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG,eval
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vec
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cvec,dvec
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bvec
      data x0/0.0d0/

    
      if (ztran .and. zpseg) open(unit=jscr,form='unformatted',recordtype='segmented')
      if (ztran) open(unit=jscr,form='unformatted')
      IF (abs(IRES) .EQ. 2) GOTO 20
      if (jrot.eq.1 .and. kmin.eq.0) then
!     J=1f: treat as a special case
         call DSTORE1(eval,kvec2,ezero,2)
         goto 100
      endif

!     SET UP THE HAMILTONIAN MATRIX (NOT FOR RESTART RUNS)

      IF (IRES .EQ. 0) then
         CALL SOLRT(MVIB,maxleg)
         WRITE(6,1050)
 1050    FORMAT(/5X,'HAMILTONIAN CONSTRUCTION COMPLETE')
      endif
      IF(.NOT.ZDIAG) return

      ezero=x0
      read(5,505,end=555) ezero
  505 format(f20.0)
  555 continue
     
      IF (ires .ge. 0) then

!     LOAD THE OFF-DIAGONAL BLOCKS
      IBASS=NBASS
      If (zdcore) then
          keval=nbass
          lwork=max(loff0,3*nbass)
          allocate(vec(nbass,nbass),diag(lwork),eval(keval))
      else
          KEVAL=MIN(IBASS,NEVAL+4)
          allocate(vec(nbass,keval),diag(nbass+loff),eval(ibass))
      endif
      call loadh(diag,mvib,vec,1)

!     DIAGONALISE THE HAMILTONIAN (TWICE IF REQUESTED)

      NOFFD=NBASS+1
      IF (KMIN .EQ. 0) WRITE(6,1000) JROT,IBASS
 1000 FORMAT(5X,'J =',I3,' ROTATIONAL STATE,',I7,' BASIS FUNCTIONS',&
             /12X,'F PARITY, ANTI-SYMMETRIC |JK> - |J-K> FUNCTIONS IN BASIS')
      IF (KMIN .GE. 1) WRITE(6,1010) JROT,IBASS
 1010 FORMAT(5X,'J =',I3,' ROTATIONAL STATE,',I7,' BASIS FUNCTIONS',&
             /12X,'E PARITY, SYMMETRIC |JK> + |J-K> FUNCTIONS IN BASIS')
      IF (IDIA .EQ. 2) THEN
         IF (IPAR .EQ. 0) WRITE(6,1020)
 1020    FORMAT(12X,'EVEN PARITY FUNCTIONS IN BASIS SET')
         IF (IPAR .EQ. 1) WRITE(6,1030)
 1030    FORMAT(12X,'ODD PARITY FUNCTIONS IN BASIS SET')
      ENDIF
     
      CALL DGROT(DIAG,MVIB,eval,VEC,1,NOFFD,ezero,lwork)
      WRITE(6,1060)
 1060 FORMAT(/5X,'DIAGONALISATION COMPLETE')
      ENDIF
      if (ztran .or. zdcore) deallocate(diag,vec)

!     TRANSFORM THE EIGENVECTORS IF REQUESTED
      IF (ZTRAN) THEN 

      ALLOCATE(bvec(mbass,nvib),cvec(nvib, neval),dvec(mbass, neval))

          CALL DSTORE(iwave,jscr,jvec,KVEC,ZPTRA,&
                      MVIB,DVEC,BVEC,CVEC,eval,nblk,1,zpseg)
          DEALLOCATE(bvec,cvec,dvec)
      ENDIF

      IF (KMIN .LE. 1) goto 100

!     DIAGONALISE A SECOND TIME IF KMIN > 1

 20   continue
      IBASS=NBASS-MVIB(1)
      NEVAL=MIN(NEVAL2,IBASS)

      if (jrot.eq.1) then
!     J=1f: treat as a special case
         call DSTORE1(eval,kvec2,ezero,2)
         goto 100
      endif

      if (.not. zdcore) KEVAL=MIN(IBASS,NEVAL+4)
      if (zdcore)       keval=ibass
      if (abs(ires).eq.2) allocate(eval(keval))
      NOFFD=IBASS+MVIB(1)*MVIB(2)+1
      JVEC=JVEC2
      IF (IRES .ge. 0)THEN
         if (ires .eq. 2)then
           ezero=x0
           read(5,505,end=556) ezero
 556      continue
         endif
       WRITE(6,1000) JROT,IBASS
       IF (IDIA .EQ. 2) THEN
         IF (IPAR .EQ. 0) WRITE(6,1020)
         IF (IPAR .EQ. 1) WRITE(6,1030)
       ENDIF

!     LOAD THE HAMILTONIAN MATRIX IF NECESSARY

       IF (ztran .or. IRES.EQ.2 .or. zdcore) then     
          If (zdcore) then
             lwork=max(loff0,3*nbass)  
             allocate(vec(nbass,ibass),diag(lwork))
          else
             allocate(vec(nbass,keval),diag(nbass+loff))
          endif           
          call loadh(diag,mvib,vec,2)
       endif
 
        CALL DGROT(DIAG(MVIB(1)+1),MVIB,eval,VEC,2,NOFFD,ezero,lwork)
        deallocate(vec,diag)
        WRITE(6,1060)
      ENDIF

!     TRANSFORM THE EIGENVECTORS IF REQUESTED

      IF (ZTRAN) THEN
       ALLOCATE(bvec(mbass,nvib),cvec(nvib, neval),dvec(mbass, neval))
          CALL DSTORE(iwave,jscr,JVEC2,KVEC2,ZPTRA,&
                     MVIB,DVEC,BVEC,CVEC,eval,nblk-1,2,zpseg)
          DEALLOCATE(bvec,cvec,dvec)
      ENDIF
 100  DEALLOCATE(eval)
      return
      END
!###########################################################################
      SUBROUTINE WRTHO(DIAG,OFFDG,MVIB)
!     PRINT HAMILTONIAN MATRIX  (out of core version)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y)

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      
      DOUBLE PRECISION, DIMENSION(NBASS) :: DIAG
      DOUBLE PRECISION, DIMENSION(LOFF) :: OFFDG
      DIMENSION MVIB(NBLK)

      WRITE(6,1010) DIAG
 1010 FORMAT(5X,'HAMILTONIAN MATRIX: DIAGONAL ELEMENTS',&
             /(10F13.8))
      IOFF2=0
      DO 10 K=2,NBLK
      IOFF1=IOFF2+1
      IOFF2=IOFF2+MVIB(K-1)*MVIB(K)
      WRITE(6,1020) K-1
 1020 FORMAT(//5X,'OFF-DIAGONAL BLOCK',I3)
      WRITE(6,1030) (OFFDG(I),I=IOFF1,IOFF2)
 1030 FORMAT(10F13.8)
   10 continue
      RETURN
      END      
      subroutine wrthi(hamil,nham)
!     print hamiltonian matrix in core version
      double precision hamil(nham,nham)
      write(6,1010)
 1010 format(5x,'hamiltonian matrix'/)
      do 30 i=1,nham
      write(6,1020) (hamil(i,j),j=1,i)
 1020 format(10f13.7)
   30 continue
      return
      end

!################################################################################
      SUBROUTINE SOLRT(MVIB,maxleg)

!     SUBROUTINE SOLRT SETS UP NON-ZERO PARTS OF HAMILTONIAN       

!     IN SETTING UP THE HAMILTONIAN THE FOLLOWING COUNTERS ARE USED:
!     K RUNS OVER VALUES OF K (PROJECTION OF JROT ON Z AXIS)
!     L RUNS OVER ANGULAR BASIS FUNCTIONS
!     M RUNS OVER R1 RADIAL BASIS FUNCTIONS
!     N RUNS OVER R2 RADIAL BASIS FUNCTIONS
!     IBASS RUNS OVER # OF BASIS FUNCTION IN THE SECULAR PROBLEM

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)

      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DOUBLE PRECISION, DIMENSION(nlim) :: rm2
      DIMENSION MVIB(NBLK)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pleg,dvrvec
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: TCOF,DIAG,OFFDG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COEF1,COEF2

      data x0/0.0d0/,sqrt2/1.4142135623731d0/

      ALLOCATE(TCOF(NVIB),COEF1(MBASS,NVIB),COEF2(MBASS,NVIB),&
               pleg((maxleg+1)*npntt),dvrvec(mbass),diag(nbass),offdg(loff0))

      REWIND ISCR
      if (.not. zdiag) then
         write(IRF2) nbass,neval,neval2,ipar,idia,jrot,kmin,nblk
         write(IRF2)mvib
         write(IRF2)kbass
      endif
!     READ RADIAL MATRIX ELEMENTS FROM FILE

      READ(iwave) rm2
      read(iwave)
      read(iwave)

      JJP1 = JROT * (JROT+1)
      jdia=max(idia,1)
      if (zquad2) n2max=1
      if (.not.zquad2) n2max=npnt2
      nrad=npnt1*npnt2
      IF (KMIN .EQ. 0) K=1
      IF (KMIN .NE. 0) K=0
      MOFF=0
      IDPT=1
!     START A NEW BAND OF THE HAMILTONIAN MATRIX
  300 MOFF=MOFF+1
!     READ IN NEXT SET OF VIBRATIONAL VECTORS & SET DIAGONAL ELEMENTS
  110 READ(iwave)
      READ(iwave)   KZ,maxleg,idvr,lincr
      ibass2=idvr*nrad

!     SKIP BACK IF KZ=0 AND WE ARE DOING AN F PARITY CALCULATION
      IF (ABS(KZ) .LT. K) THEN
         READ(iwave)
         READ(iwave) meval
         do 105 i=1,meval+1
         read(iwave)
  105    continue
         GOTO 110
      ENDIF

      maxlg2=maxleg+lincr
      IF (MOD(maxlg2,JDIA) .NE. IPAR) maxlg2=maxlg2-1

      IF (MVIB(MOFF) .GT. 0) THEN

         call getrow(pleg,(maxleg+1)*idvr,iwave)
         read(iwave) meval
        
!        read diagonal matrix elements

         CALL GETROW(DIAG(IDPT),MVIB(MOFF),iwave)

         IDPT=IDPT+MVIB(MOFF)    
!        read basis vectors and transform to associated Legendres
         call jtran(coef2,mvib(moff),pleg,maxleg,idvr,k,dvrvec,&
                    nrad,ibass2,jdia,lincr,jst2)
      ELSE
         READ(iwave)
         READ(iwave) meval
         do 106 i=1,meval+1
         read(iwave)
  106    continue
      ENDIF

!     COMPUTE OFF DIAGONAL ELEMENTS (NONE FOR FIRST TIME THROUGH)

      IF (MOFF .EQ. 1) GOTO 250
      IF (MN .EQ. 0) GOTO 250
!     ZERO THE NEXT OFF-DIAGONAL BLOCK
      OFFDG = X0
      KKP1=K*(K-1)
      ii=0
      j=0
!     line up the angular functions in the bra and ket
      ijump=(jst2-jst1)/jdia
      jjump=maxlg2/jdia-maxlg1/jdia
      do 200 n1=1,npnt2
      do 205 m1=1,npnt1
      ii=ii+ijump
      do 210 l1=jst2,maxlg1,jdia

      LLP1=L1*(L1+1)
      CJLP = -SQRT(DBLE((LLP1-KKP1)*(JJP1-KKP1)))

!     SPECIAL CASE: K1=0
      IF (K .EQ. 1) CJLP = SQRT2 * CJLP
      ii=ii+1
      j=j+1
      if (.not.zquad2) i=mod(ii-1,nblk1)+1
      do 220 n2=1,n2max
      if (zquad2) then
         i=ii
         IF (ZEMBED) THEN
            COR = CJLP * rm2(N1) 
         ELSE
            COR = CJLP * rm2(M1)
         ENDIF
      else
         cor = cjlp * rm2((n1-1)*npnt2+n2)
      endif

      call dger(mvib(moff), mvib(moff - 1), cor, coef2(j, 1), &
           mbass, coef1(i, 1), mbass, offdg, mvib(moff))

      i=i+nblk1
  220 CONTINUE
  210 CONTINUE
      j=j+jjump
  205 continue
  200 CONTINUE

!     DUMP COMPLETED BLOCK TO UNIT ISCR

      CALL OUTROW(OFFDG,MN,ISCR)

!     NEXT BLOCK: HAVE WE FINISHED?
  250 K=K+1
      IF (K .GT. JROT) GOTO 400
!     size the NEXT OFF-DIAGONAL BLOCK
      MN=MVIB(MOFF)*MVIB(MOFF+1)

!     PREPARE TO READ IN NEXT SET OF VECTORS

      IBASS1=IBASS2
      maxlg1=maxlg2
      jst1=jst2
      nblk1=ibass2/npnt2
      COEF1 = COEF2
      GOTO 300

!     PLACE DIAGONAL ELEMENTS AT END OF SCRATCH FILE

  400 CALL OUTROW(DIAG,NBASS,ISCR)
      DEALLOCATE(TCOF,COEF1,COEF2,pleg,dvrvec,diag,offdg)

      RETURN
      END

!########################################################################
      subroutine jtran(coef,mvib,pleg,maxleg,idvr,kz,dvrvec,nrad,&
                      ibass2,jdia,lincr,jstart)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DOUBLE PRECISION, DIMENSION(0:MAXLEG,idvr) :: PLEG
      DOUBLE PRECISION, DIMENSION(idvr,nrad) :: dvrvec
      DOUBLE PRECISION, DIMENSION(mbass,mvib) :: coef
      DOUBLE PRECISION, DIMENSION(nrad) :: sumk

      DATA X0/0.0D0/

!     DVR3DRJ: TRANSFORM BACK TO THE ORIGINAL FBR-TYPE BASIS IN THE
!     ASSOCIATED LEGENDRE FUNCTIONS

      JSTART=KZ
      JJ0=-JDIA
      IF (MOD(JSTART,JDIA) .NE. IPAR) THEN
         JJ0=JJ0+1
         JSTART=JSTART+1
      ENDIF
      NANG=(MAXLEG+LINCR-JSTART)/JDIA+1

      DO 10 L=1,mvib
!     first read in a new vector
      IF (z1da) then
!        DVR1DA vectors: already as ASSOCIATED LEGENDRE FUNCTIONS
         call getrow(coef(1,l),ibass2,iwave)
         goto 10
      ENDIF

      call getrow(dvrvec,ibass2,iwave)
      JJ=JJ0
      DO 20 J=1,nang         
         sumk = x0
         JJ=JJ+JDIA        
         do 40 k=1,idvr
            do 50 mn=1,nrad
               sumk(mn)=sumk(mn) + dvrvec(k,mn) * pleg(jj,k)
   50       continue
   40    CONTINUE
      ipt=j
      do 60 mn=1,nrad
      coef(ipt,l) = SUMK(mn)
      ipt=ipt+nang
   60 continue
   20 CONTINUE
   10 CONTINUE
!     skip vectors that are not needed
      do 100 l=mvib+1,meval
      read(iwave)
  100 continue
!     If vectors to be back transformed, store for later use
      If (ztran) then
         write(jscr) jstart,nang,nrad*nang
         call outrow(coef,mbass*mvib,jscr)
      endif

      RETURN
      END

!############################################################################
      SUBROUTINE loadh(diag,mvib,hamil,itime)

!     SUBROUTINE loadh loads the Hamiltonian matrix from disk

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                   KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                   kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                   iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                   zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DIMENSION MVIB(NBLK)
      DOUBLE PRECISION, dimension(*) :: diag
      DOUBLE PRECISION, dimension(ibass,ibass) :: hamil

      data x0/0.0d0/

      REWIND ISCR

      if (.not. zdcore) then
         IPT=1+nbass
         DO 10 J=2,NBLK
         LENG=MVIB(J)*MVIB(J-1)
         IF (LENG .GT. 0) CALL GETROW(diag(IPT),LENG,ISCR)
         IPT=IPT+LENG
   10    continue
         CALL GETROW(DIAG,NBASS,ISCR)

!     PRINT HAMILTONIAN MATRIX- IF REQUESTED

         IF (ZPHAM) CALL WRTHO(DIAG,diag(nbass+1),MVIB)
      else

         hamil = x0

         if (itime.gt.1 .and. mvib(1)*mvib(2).gt.0) read(iscr)
         ist2=0
         DO 20 J=itime+1,NBLK
         ist1=ist2
         ist2=ist2+mvib(j-1)
         LENG=MVIB(J)*MVIB(J-1)

         IF (LENG .GT. 0) CALL GETROW(diag,LENG,ISCR)
         ipt=0
         do 30 i1=ist1+1,ist1+mvib(j-1)
         do 40 i2=ist2+1,ist2+mvib(j)
         ipt=ipt+1
         hamil(i1,i2)=diag(ipt)
         hamil(i2,i1)=diag(ipt)
   40    continue
   30    continue
   20    continue
         CALL GETROW(diag,nbass,ISCR)
         n=nbass-ibass
         do 50 i=1,ibass
         hamil(i,i)=diag(n+i)
   50    continue
!     PRINT HAMILTONIAN MATRIX- IF REQUESTED

         IF (ZPHAM) CALL WRTHI(hamil,ibass)

      endif
      REWIND ISCR
      return
      END

!########################################################################
      SUBROUTINE DGiter(DIAG,MVIB,EVAL,VEC,lwork,K1,NOFFD)

!     SUBROUTINE Dgiter solves THE EIGENVALUE PROBLEM:
!          HAMIL * VEC = EVAL * VEC
!     BY USING ITERATIVE NAG ROUTINE F02FJF TO DO DIAGONALISATIONS.

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg
      double precision, external :: vecvec
      external matvec,f02fjz
      save eshift

      DOUBLE PRECISION, DIMENSION(KEVAL) :: EVAL
      DOUBLE PRECISION, DIMENSION(*) :: DIAG
      DOUBLE PRECISION, DIMENSION(IBASS,keval) :: VEC
      DIMENSION MVIB(NBLK)
      DOUBLE PRECISION, DIMENSION(lwork) :: WORK
              
      DATA X0/0.0D0/,X1/1.0D0/,EMAX/1.0D50/

!     CREATE SOME STARTING VECTORS BY USING THE DIAGONAL ELEMENTS
      vec=x0
      ESMALL=-EMAX
      DO 10 I=1,KEVAL
      EBIG=EMAX
      DO 20 J=1,IBASS
      IF (DIAG(J) .GT. EBIG) GOTO 20
      IF (DIAG(J) .LE. ESMALL) GOTO 20
      IND=J
      EBIG=DIAG(J)
   20 CONTINUE
      VEC(IND,I)=X1
      ESMALL=EBIG
   10 CONTINUE

!     SHIFT DIAGONALS TO ENSURE WE GET THE LOWEST EIGENVALULES

      IF (ztran .or. K1 .LE. 1 .OR. IRES .EQ. 2) THEN
         ESHIFT=EBIG
         DO 30 I=1,IBASS
         ESHIFT=MAX(ESHIFT,DIAG(I))
   30    CONTINUE
         DO 40 I=1,IBASS
         DIAG(I)=DIAG(I)-ESHIFT
   40    CONTINUE
      ENDIF

!     DIAGONALISE  THE HAMILTONIAN
      IFAIL=1
      NOITS=10000
      CALL F02FJF(IBASS,NEVAL,KEVAL,NOITS,TOLER,VECVEC,MATVEC,F02FJZ,&
                  KEVAL,VEC,IBASS,EVAL,WORK,LWORK,diag,NOFFD,&
                  MVIB,K1,IFAIL)
      IF (IFAIL .NE. 0) WRITE(6,900) IFAIL
  900 FORMAT(//5X,'F02FJF RETURNED IFAIL =',I3)

      WRITE(6,1000) NOITS
 1000 FORMAT(/5X,'CONVERGENCE AFTER NOITS =',I6,' ITERATIONS')
!     SHIFT EIGENVALUES BACK
      DO 50 I=1,NEVAL
      EVAL(I)=EVAL(I)+ESHIFT
   50 CONTINUE
      return
      end

!###########################################################################
      SUBROUTINE DGROT(DIAG,MVIB,evalcm,VEC,K1,NOFFD,ezero,lwork)

!     SUBROUTINE DIAG SOLVES THE EIGENVALUE PROBLEM:       
!          HAMIL * VEC = EVAL * VEC
!     BY USING ITERATIVE NAG ROUTINE F02FJF TO DO DIAGONALISATIONS.
!     or in core         Lapack routine dsyev

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DOUBLE PRECISION, DIMENSION(neval) :: eval
      DOUBLE PRECISION, DIMENSION(neval) :: evalcm
      DOUBLE PRECISION, DIMENSION(*) :: DIAG
      DOUBLE PRECISION, DIMENSION(IBASS,*) :: VEC
      REAL*8, allocatable :: work(:),evec(:,:)
      integer, allocatable :: iwork(:),iwork2(:)
      DIMENSION MVIB(NBLK),IBIG(IBASS)

!          AUTOCM CONVERTS ATOMIC UNITS (HARTREE) TO CM-1.
      DATA AUTOCM/2.19474624D+05/,x0/0.0d0/

      ip=0 
      if (zdcore) then
         IFAIL=0
         allocate (work(ibass*8),iwork(5*ibass),iwork2(ibass),evec(ibass,neval))
!         call dsyev ('V','U',ibass,vec,ibass,eval,diag,lwork,ifail)
         call DSYEVX('V','I','U',ibass,vec,ibass,0.d0,0.d0,1,neval,&
                     0.d0,nfound,eval,evec,ibass, WORK, 8*ibass, IWORK,&
                     IWORK2, IFAIL )
         IF (IFAIL .NE. 0) WRITE(6,900) IFAIL
         write(6,*)' Found ',nfound,' evalues (out of ',neval,')'
         IF (NEVAL.ne.nfound) stop
         do i=1,neval
            do j=1,ibass
               vec(j,i)=evec(j,i)
            end do
         end do
         deallocate(work,iwork,iwork2,evec)
  900    FORMAT(//5X,'Diagonaliser DSYEV failed with IFAIL =',I3)
      else
         LWORK=3*KEVAL+MAX(KEVAL*KEVAL,NBASS+NBASS)
         call DGiter(DIAG,MVIB,EVAL,VEC,LWORK,K1,NOFFD)
      endif

!     PRINT EIGENVALUES IN ATOMIC UNITS & WAVENUMBERS

      WRITE(6,1000) NEVAL
 1000 FORMAT(//5X,'LOWEST',I4,' EIGENVALUES IN HARTREES',/)
      WRITE(6,1020) EVAL
 1020 FORMAT(5D24.12)
      IF (ZPFUN) THEN
         IF (K1 .EQ. 1) THEN
            OPEN(UNIT=ILEV,FORM='FORMATTED')
            REWIND ILEV
  200       READ(ILEV,*,END=210,ERR=210)
            GOTO 200
  210       CONTINUE
! ***** INCLUSION OF THE FOLLOWING CARD IS MACHINE DEPENDENT *****
!           backspace ilev
         ENDIF
         IP=1-KMIN
         IF (KMIN .GT. 1) IP=K1-1
         WRITE(ILEV,1025) JROT,IP,IDIA,IPAR,0,NEVAL,IBASS
 1025    FORMAT(7I6)
         WRITE(ILEV,1026) EVAL
 1026    FORMAT(4D20.12)
      ENDIF
      IF (ZVEC) THEN
!        WRITE EIGENVALUES, EIGENVECTORS, ETC TO STREAM JVEC
         KZ=KMIN
         IF (KMIN .GT. 1) KZ=2-K1
         if (zpseg==.true.) then 
            OPEN(UNIT=JVEC,FORM='UNFORMATTED',recordtype='segmented')
         else 
            OPEN(UNIT=JVEC,FORM='UNFORMATTED')
         end if 
         REWIND JVEC
         WRITE(JVEC) JROT,KZ,IPAR,NEVAL,IBASS
         WRITE(JVEC) (MVIB(K),K=K1,NBLK)
         CALL OUTROW(EVAL,NEVAL,JVEC)

! GJH CODE TO OUTPUT IN ROWS OF K INTEAD OF ROWS OF I
         mend=0
! loop over k
         do 543 k=k1,nblk
! set begin of k block
         mbeg=mend+1
! set end of k block
         mend=mend+mvib(k)
! write entire row, all of same k
         if(mvib(k).gt.0) write(jvec) ((vec(j,i),j=mbeg,mend),i=1,neval)
543      continue
! END OF GJH CODE.

      ENDIF
      DO 60 I=1,NEVAL
      EVALCM(i) = EVAL(i)*AUTOCM - ezero
   60 CONTINUE
      WRITE(6,1010) NEVAL,ezero
 1010 FORMAT(//5X,'Lowest',I4,' eigenvalues in wavenumbers relative to',&
                  ' EZERO =',D24.12/)
      WRITE(6,1022) EVALCM
 1022 FORMAT(1x,10f13.5/)

      If (IP.eq.0) then
      do i=1,neval
         write(62,1034)i,eval(i),evalcm(i)
      end do
      ELSE
      do i=1,neval
         write(63,1034)i,eval(i),evalcm(i)
      end do
      end if

1034 format(1I4,2f30.20)

      IF (ZPVEC) THEN
          IF (THRESH .LE. X0) THEN
!             PRINT COMPLETE EIGENVECTORS
              WRITE(6,1030)
 1030         FORMAT(//,'    EIGENVECTORS',/)
              DO 70 I=1,NEVAL
              WRITE(6,1040) (VEC(J,I),J=1,IBASS)
 1040         FORMAT(/(1X,10F13.7))
   70         CONTINUE
          ELSE
!             PRINT LARGEST COMPONANTS OF THE EIGENVECTORS
              WRITE(6,1050) THRESH
 1050         FORMAT(//,'EIGENVECTOR COMPONANTS GREATER THAN THRESH =',&
                        F5.2)
              DO 80 I=1,NEVAL
              VBIG=X0
              IPT=0
              DO 90 J=1,IBASS
              VV=ABS(VEC(J,I))
              IF (VV .GT. THRESH) THEN
                  IPT=IPT+1
                  IBIG(IPT)=J
              ENDIF
              IF (IPT .LE. 0 .AND. VV .GT. VBIG) THEN
                  VBIG=VV
                  IBIG(1)=J
              ENDIF
   90         CONTINUE
              WRITE(6,1060) I,(IBIG(J),VEC(IBIG(J),I),J=1,MAX(1,IPT))
 1060         FORMAT('  VECTOR',I3,5(I7,F14.10)/(11X,5(I7,F14.10)))
   80         continue
          ENDIF
      ENDIF
      RETURN
      END

!##################################################################################
 
      SUBROUTINE DSTORE(iwave,jscr,jvec,KVEC,ZPTRA,&
                       MVIB,D,B,C,energy,lblk,itra,zpseg)

!     DSTORE TRANSFORMS THE EIGENVECTORS OF THE SECOND VARIATIONAL #013
!     STEP INTO ONES FOR THE FIRST STEP BASIS AND STORES THE
!     RESULTS IN A FORM SUITABLE FOR program DIPOLE3.

      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt


      DIMENSION MVIB(lblk),NKBAS(lblk),lmin(lblk),lbasis(lblk)
      DOUBLE PRECISION, DIMENSION(MBASS,NVIB) :: B
      DOUBLE PRECISION, DIMENSION(nvib, neval) :: C
      DOUBLE PRECISION, DIMENSION(MBASS, neval) :: D
      DOUBLE PRECISION, DIMENSION(KEVAL) :: ENERGY
      DOUBLE PRECISION, DIMENSION(3) :: XMASS

      DATA X0/0.0D0/,x1/1.0d0/

      WRITE(6,1000) iwave,jscr,jvec,KVEC
 1000 FORMAT(5X,'EIGENVECTOR TRANSFORMATION:',&
             /5X,'INPUT:  DVR3D data,            IWAVE =',I3,&
             /5X,'        ROTLEV3 data (scratch) JSCR  =',I3,&
             /5X,'        ROTLEV3 data (final)   JVEC  =',I3,&
             /5X,'OUTPUT: transformed vectors,   KVEC  =',I3/)
!     Read DVR3D header
      rewind iwave
      read(IWAVE) IDIA,IPAR,IDVR,NPNT1,NPNT2,JROT0,KMIN0,MVAL
      read(IWAVE) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,z1da,ZQUAD2
      read(IWAVE) RE1,DISS1,WE1,RE2,DISS2,WE2
!     READ ROTLEV HEADER
      REWIND JVEC
      READ(JVEC) JROT1,KMIN1,IPAR1,NVAL,IBASS
      READ(JVEC) MVIB

!     CHECK FOR COMPATABILITY
      IF (JROT1 .NE. ABS(JROT0)) THEN
          WRITE(6,900) JROT1,ABS(JROT0)
  900     FORMAT(/5X,'J LEVELS MISMATCHED',&
                 /5X,'JROT1 =',I3,'  JROT0 =',I3)
          STOP
      ENDIF
      IF (KMIN1 .GT. KMIN0) THEN
         WRITE(6,910)
  910    FORMAT(/5X,'KMIN1 AND KMIN0 INCOMPATIBLE')
         STOP
      ENDIF
      IF (IPAR .NE. IPAR1) THEN
         WRITE(6,920) IPAR,IPAR1
  920    FORMAT(/5X,'PARITIES MISMATCHED'/5X,'IPAR= ',I2,'  IPAR1= ',I2)
         STOP
      ENDIF
      IF (ZPTRA) WRITE(6,1010) JROT1,KMIN1,IPAR,NVAL
 1010 FORMAT(/5X,'J =',I3,' KMIN =',I2,'  PARITY =',I2,' NVAL =',I3/)
!     Read basis set parameters from scratch file JSCR
      rewind jscr
      IF (itra .GT. 1) THEN
          READ(jscr)
          READ(jscr)
      ENDIF
      nend=0
      DO 20 K=1,lblk
      read(jscr) lmin(k),lbasis(k),nkbas(k)
      READ(jscr)
      nend=nend+nkbas(k)
   20 CONTINUE
      IF (NEND .GT. KBASS) THEN
         WRITE(6,930) NEND,KBASS
  930    FORMAT(/5X,'ERROR IN SIZE OF TOTAL BASIS',&
                /5X,'NEND =',I7,'  KBASS =',I7/)
         STOP
      ENDIF
!     WRITE HEADER ON NEW FILE
      if (zpseg==.true.) then
          write(235,*)"banane" 
          OPEN(UNIT=KVEC,FORM='UNFORMATTED',recordtype='segmented')
      else 
          OPEN(UNIT=KVEC,FORM='UNFORMATTED')
      end if 
      write(kvec) IDIA,IPAR,idvr,npnt1,npnt2,JROT1,KMIN1,NVAL
      write(kvec) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,.false.
      write(kvec) RE1,DISS1,WE1,RE2,DISS2,WE2
      write(kvec) nend,lmin,lbasis,nkbas
!     transfer the DVR points
      read(iwave)
      CALL getROW(b,NPNT1,iwave)
      CALL OUTROW(b,NPNT1,kvec)
      CALL getROW(b,NPNT2,iwave)
      CALL OUTROW(b,NPNT2,kvec)
      write(kvec) nval
      CALL GETROW(ENERGY,NVAL,JVEC)
      CALL OUTROW(ENERGY,NVAL,KVEC)
!     START TRANSFORMATION STEP
      REWIND jscr
      KK=KMIN1+KMIN1-1


! begining of GJH's code
      IF (itra .gt. 1) THEN
          READ(jscr)
          READ(jscr)
      ENDIF

! loop over k
      DO 80 K=1,lblk

! set dmatrix elements to zero
         d=0.0d0

         IF (MVIB(K) .GT. 0) THEN
! read in B matrix 
            read(jscr)
            CALL GETROW(B,MBASS*MVIB(K),jscr)
! read in C matrix from jvec stream
            READ(jvec) ((C(i,j),i=1,MVIB(K)),j=1, NEVAL)
! multiply matrices to obtain d coefficents.
            call dgemm('N', 'N', nkbas(k), nval, mvib(k), x1, b, &
                 mbass, c, nvib, x0, d, mbass)
         else
            read(jscr)
            read(jscr)
         endif

! next ouput d coefficents in rows of k
         call outrow(d,nkbas(k)*nval,kvec)
         if (zptra) then
            write(6,1020) k
 1020       format(/'Vectors for K-block',i3/)
            do 50 i=1,nval
            write(6,1030) (d(j,i),j=1,nkbas(k))
 1030       format(1X,10F13.7)
   50       continue
         endif
   80 continue   

! end of GJH's code
      WRITE(6,1050) JROT1,KMIN1,IPAR,NVAL,NEND
 1050 FORMAT(//5X,'TRANSFORMATION COMPLETED SUCCESSFULLY ',&
             'FOR STATE: JROT =',I3,', KMIN =',I3,', IPAR =',I3,&
             /5x,'STORED',I4,' VIB-ROT LEVELS',5X,'NEND =',I7,/)
      REWIND JVEC 
      close(unit=kvec)
      RETURN
      END
!##################################################################################
 
      SUBROUTINE DSTORE1(eval,kvec1,ezero,itra)
                       
!     `Transformation' step for J=1f special case      
!     RESULTS IN A FORM SUITABLE FOR program DIPOLE3.

      IMPLICIT DOUBLE PRECISION(A-H,O-Y), LOGICAL(Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      COMMON /OUTP/ toler,thresh,ilev,iwave,jscr,jvec,jvec2,kvec,kvec2,&
                    iscr,ires,irf1,irf2,zpham,zpvec,zdcore,z1da,&
                    zembed,zvec,zquad2,zpfun,zdiag,ztran,zptra,zpseg

      DOUBLE PRECISION, DIMENSION(max(mbass,neval,(npntt+1)**2)) :: B
      DOUBLE PRECISION, DIMENSION(neval) :: eval
      DOUBLE PRECISION, DIMENSION(3) :: XMASS
      DATA AUTOCM/2.19474624D+05/

!     Read DVR3D header
      rewind iwave
      read(IWAVE) IDIA,IPAR,IDVR,NPNT1,NPNT2,JROT,KMIN0,MVAL
      read(IWAVE) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,z1da,ZQUAD2
      read(IWAVE) RE1,DISS1,WE1,RE2,DISS2,WE2
      kmin1=0
!     Read basis set parameters from scratch file JSCR
      rewind jscr
      IF (itra .GT. 1) THEN
          READ(jscr)
          READ(jscr)
      ENDIF
      read(jscr) lmin,lbasis,nkbas
      READ(jscr)
      nend=nkbas
      IF (NEND .GT. KBASS) THEN
         WRITE(6,930) NEND,KBASS
  930    FORMAT(/5X,'ERROR IN SIZE OF TOTAL BASIS',&
                /5X,'NEND =',I7,'  KBASS =',I7/)
         STOP
      ENDIF
      if (ztran) then
!        WRITE HEADER ON NEW FILE
		 if (zpseg==.true.) then 
            OPEN(UNIT=KVEC1,FORM='UNFORMATTED',recordtype='segmented')
         else 
            OPEN(UNIT=KVEC1,FORM='UNFORMATTED')
         end if 
         write(kvec1) IDIA,IPAR,idvr,npnt1,npnt2,JROT,KMIN1,NEVAL
         write(kvec1) ZEMBED,ZMORS1,ZMORS2,XMASS,G1,G2,.false.
         write(kvec1) RE1,DISS1,WE1,RE2,DISS2,WE2
         write(kvec1) nend,lmin,lbasis,nkbas
      endif
!     transfer the DVR pointsn
      read(iwave)
      CALL getROW(b,NPNT1,iwave)
      if (ztran) CALL OUTROW(b,NPNT1,kvec1)
      CALL getROW(b,NPNT2,iwave)
      if (ztran) CALL OUTROW(b,NPNT2,kvec1)
      IF (itra .GT. 1) THEN

          READ(iwave)
          READ(iwave)
          READ(iwave)
          READ(iwave) meval
          do 10 i=1,meval+1
          READ(iwave)
  10      continue
      ENDIF
      CALL getROW(b,idvr,iwave)
      if (ztran) CALL OUTROW(b,idvr,kvec1)
      read(iwave) kz,maxleg,nidvr,lincr
      if (ztran) write(kvec1) kz,maxleg,nidvr,lincr
      nlegs=(maxleg+1)*nidvr
      call getrow(b,nlegs,iwave)
      if (ztran) call outrow(b,nlegs,kvec1)
      read(iwave) meval
      if (ztran) write(kvec1) neval
!     read DVR3DRJZ energies

      CALL GETROW(eval,neval,iwave)
      if (ztran) call outrow(eval,neval,kvec1)
!     PRINT EIGENVALUES IN ATOMIC UNITS & WAVENUMBERS

      WRITE(6,1000) NEVAL
 1000 FORMAT(//5X,'LOWEST',I4,' EIGENVALUES IN HARTREES',/)
      WRITE(6,1020) EVAL
 1020 FORMAT(5D24.12)
      IF (ZPFUN) THEN
         IF (itra .EQ. 1) THEN
            OPEN(UNIT=ILEV,FORM='FORMATTED')
            REWIND ILEV
  200       READ(ILEV,*,END=210,ERR=210)
            GOTO 200
  210       CONTINUE
! ***** INCLUSION OF THE FOLLOWING CARD IS MACHINE DEPENDENT *****
!           backspace ilev
         ENDIF
         IP=1
         WRITE(ILEV,1025) JROT,IP,IDIA,IPAR,0,NEVAL,IBASS
 1025    FORMAT(7I6)
         WRITE(ILEV,1026) EVAL
 1026    FORMAT(4D20.12)
      ENDIF
      DO 60 I=1,NEVAL
      EVAL(i) = EVAL(i)*AUTOCM - ezero
   60 CONTINUE
      WRITE(6,1010) NEVAL,ezero
 1010 FORMAT(//5X,'Lowest',I4,' eigenvalues in wavenumbers relative to',&
                  ' EZERO =',D24.12/)
      WRITE(6,1022) EVAL
 1022 FORMAT(1x,10f13.5/)
      if (.not.ztran) return

      WRITE(6,2000) iwave,jscr,KVEC1
 2000 FORMAT(5X,'EIGENVECTOR DUMMY TRANSFORMATION J=1f case:',&
             /5X,'INPUT:  DVR3D data,            IWAVE =',I3,&
             /5X,'        ROTLEV3 data (scratch) JSCR  =',I3,&
             /5X,'OUTPUT: transformed vectors,   KVEC  =',I3/)
      IF (ZPTRA) WRITE(6,2010) JROT,KMIN1,IPAR,NEVAL,1
 2010 FORMAT(/5X,'J =',I3,' KMIN =',I2,'  PARITY =',I2,' NVAL =',I3/&
              /'Vectors for K-block',i3/)
!     copy across DVR3DRJZ vectors
      do 220 i=1,neval
      call getrow(b,nkbas,iwave)
      call outrow(b,nkbas,kvec1)
      if (zptra) write(6,1030) (b(j),j=1,nkbas)
 1030       format(1X,10F13.7)
  220 continue  

       WRITE(6,2050) JROT,KMIN1,IPAR,NEVAL,NEND
 2050 format(//5x,'Transfer of vectors completed successfully',&
             //5x,'for rotational state: jrot =',i3,', kmin =',i3,&
             ', ipar =',i3,'.  stored',i6,' vib-rot levels',&
             /5x,'nend =',i7)
      close(unit=kvec1)
      RETURN
      END

!###################################################################################
      SUBROUTINE GETROW(ROW,NROW,IUNIT)
!     FAST NON-FORMATTED READ                           
      DOUBLE PRECISION, DIMENSION(NROW) :: ROW
      READ(IUNIT) ROW
      RETURN
      END

!##################################################################################
      SUBROUTINE OUTROW(ROW,NROW,IUNIT)
!     FAST NON-FORMATTED WRITE 
      DOUBLE PRECISION, DIMENSION(NROW) :: ROW                        
      WRITE(IUNIT) ROW
      RETURN
      END

!################################################################################
      FUNCTION VECVEC(IFLAG,N,Z,W,D1,D2,D3,D4)

!     VECVEC RETURNS THE DOT PRODUCT W.Z FOR F02FJF                 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(N) :: W,Z
      VECVEC=dDOT(N,W,1,Z,1)
      RETURN
      END

!################################################################################
      SUBROUTINE MATVEC(IFLAG,N,Z,W,HAMIL,NOFFD,MVIB,K1)

!     MATVEC PERFORMS W = HAMIL * Z FOR F02FJF            
!     WRITTEN TO TAKE ADVANTAGE OF THE BLOCK STRUCTURE
!     NOTE:
!     HAMIL CONTAINS ARRAYS DIAG & OFFDG RELYING ON THEM BEING
!     ADJACENT IN THE DYNAMIC STORE ALLOCATION

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SIZE/ NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,&
                    KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,&
                    kbass,npnt1,npnt2,npntt
      DOUBLE PRECISION, DIMENSION(NBASS) :: W,Z
      DOUBLE PRECISION, DIMENSION(*) :: HAMIL
      DIMENSION MVIB(NBLK)
      parameter (x1=1.0d0)
!     FIRST INCLUDE THE DIAGONAL CONTRIBUTION
      DO 10 I=1,IBASS
      W(I) = Z(I) * HAMIL(I)
   10 CONTINUE
!     THEN THE OFF-DIAGONAL
      IOFF=NOFFD
      I2=1
      DO 20 K=K1,NBLK
      IF (MVIB(K) .LT. 1) GOTO 20
      IF (K .GT. K1) then 
      if (MVIB(K-1) .GE. 1) THEN

         call dgemv('N',MVIB(K), MVIB(K-1), x1, HAMIL(IOFF),&
              MVIB(K), Z(I1),1,x1, W(I2), 1)

          IOFF=IOFF+MVIB(K-1)*MVIB(K)
       endif
      ENDIF
      I1=I2
      I2=I2+MVIB(K)
      IF (K .LT. NBLK) then
      if (MVIB(K+1) .GE. 1) then

         call dgemv('T', mvib(k+1), mvib(k), x1, hamil(ioff),&
              mvib(k+1), z(i2),1,x1 ,w(i1), 1)
         endif
         endif
   20 CONTINUE
      RETURN
      END
