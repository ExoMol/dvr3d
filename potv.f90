
     subroutine potv(V_h,R1,R2,xcos)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y),logical(z)
      COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)
      DATA RZ/1.16   /,RHO/0.00000000  /
      DATA TOANG/0.5291772/, cmtoau/219474.624d0/
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/
      real*8 :: xp(1:300)
      integer :: i_t,npropin
      character(len=10) label
      character(len=20) coeffs_file
      pi=dacos(-1.d0)

      coeffs_file='SO2B2.fit'
!  SO2 B2 state PES . Emil Zak, 6 Dec 2016

!############### Define coordinates #################"
      IF (G1 .EQ. X0) THEN 
         Q1 = R1
         Q2 = R2
         THETA = ACOS(XCOS) 
      ELSE IF (G2 .EQ. X0) THEN 
         XX = R1 * G1
         YY = R1 * (X1 - G1)
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

         THETA = ACOS(COST)
      ELSE 
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
!############### Define coordinates #################"


      open(unit=32,status='old',form='formatted',file=coeffs_file)
      read(32,*) npropin ! Read the number of PES coefficients

      do i = 1,npropin
         read(32,'(A10,I8,ES22.6)') label,i_t, xp(i)
      end do

      close(unit=32)
      !

      !
      R1_ang=Q1*0.5291772d0 ! Into angstroms
      R2_ang=Q2*0.5291772d0 ! Into angstroms
        vp=0
     call fit(vp,R1_ang,R2_ang,THETA,npropin,xp)
    Thup=130.0
    Thlow=85.0
    Thup=Thup*pi/180.0
    Thlow=Thlow*pi/180.0
    
    if(R1_ang .gt. 1.3 .and. R1_ang .lt. 2.0 .and.R2_ang .gt. 1.3 .and. R2_ang .lt. 2.0 .and. THETA .gt. Thlow .and. THETA .lt. Thup) then
        V_h=vp/219474.624d0
        else 
        V_h=4.0
        end if
      end subroutine





