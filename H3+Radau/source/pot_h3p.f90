

      SUBROUTINE POTV(VES,R1,R2,XCOS)     
!
! H3+ global potential due to 
! O.L. Polyansky, R. Prosmiti, W. Klopper and J. Tennyson, 
! Mol. Phys., 98, 261-273 (2000).  Fit 2.
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common /potential/ cv(100), ifit(100), der(100)


      DATA Ediss/0.169361d0/,Esp/0.131126703d0/,m/2/,ga1/1.1d05/

!*******************************************************************
! units: Esp->HARTREE, Ediss->HARTREE, ga0->{HARTREE}-1, ga1->{HARTREE}-3
!*******************************************************************

      e0=0.5d0*(Ediss+Esp)
      de1=e0-Esp 
      ga0=1.14877995d0/de1

      call  sotv(V2,R1,R2,XCOS)
      call RJ(R1,R2,XCOS,R1j,R2j,XCOSj)
      call potv1(V1,R1j,R2j,XCOSj) 
      v22=v2 /219474.6354d0
      call potv11(V4,R1j,R2j,XCOSj,v22) 
      
      v2=v4
      de=v1-e0
      ga=ga0+ga1*(de**m)
      ff=0.5d0*(1+dtanh(ga*de))
      ff1=0.5d0*(1+dtanh(-ga*de))
      do i=1,100
         der(i)=der(i)*ff1
      enddo
      VES=ff*V1+ff1*V2 
      
      return
      end 

!*******************************************************************
!  Long range term VLIM 
!  [C.F.Giese, W.R.Gentry, Phys. Rev. A, 10/6,2156,1974]
!  [W. Kolos, L.Wolniewicz, JCP, 46/4, 1426,1967]
!  [D.G.Truhlar, Int.J. of Quant. Chem.,VI,975,1972]
!    UNITS: HARTREE & BOHR
!
!*******************************************************************

      FUNCTION VLIM(R1,R2,XCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)

!      RHO=R1-1.40083D0
      RHO=R1-1.39990951D0
      A0=2.6091D0+(2.246D0+(0.3181D0-0.1194D0*RHO)*RHO)*RHO
      A2=0.60735D0+(1.3586D0+(0.5573D0-0.3170D0*RHO)*RHO)*RHO
      Q2=0.45886D0+(0.53223D0+(0.03234D0-0.091474D0*RHO)*RHO)*RHO
      P2=1.5D0*XCOS*XCOS-0.5D0
      VLIM=-(A0+A2*P2)/R2**4+Q2*P2/R2**3
!
      RETURN
      end
!c
!c*******************************************************************
!  Long range term VLIM1 
!  [B.H. Bransden, Atomic collision theory,
!   W.A. Benjamin publishers,1970,New York) 
!    UNITS: HARTREE & BOHR
!c
!c*******************************************************************
!c 
      FUNCTION VLIM1(R1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
!c
       RR=R1
       a1=4.5d0
       a2=15.d0
       b1=5.375d0
       VLIM1=-a1/RR**4-a2/RR**6
!c
       RETURN
       end

!c****************************************************************************
!  H3+  Many Body Expansion POTENTIAL and Long range contributions (POTV1) 
!  [J.N. Murrell, S. Carter, S.C. Farantos, P. Huxley, and A.J.C. Varandas,
!   Molecular Potential Energy Functions, Wiley, New York, 1984, 
!   Chapter 11,  pp. 154-160]
!  [D.M. Bishop and S.-K. Shih,  J. Chem. Phys. 64, (1976) 162]
!  [C. Schwartz and R.J. Le Roy, J. Mol. Spectr. 121, (1987) 420] 
!  [R. Schinke, M. Dupuis and W.A. Lester, Jr, J. Chem. Phys. 72, (1980) 3909] 
!  [C.F. Giese and W.R. Gentry, Phys. Rev. A 10, (1974) 2156.]  
!c
!     UNITS: HARTREE & BOHR
!c****************************************************************************
      SUBROUTINE POTV1(Vh3p,xx1,xx2,xx3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(n6=3)
      DIMENSION CA(30),CB(10),C(10)
      dimension q(n6)
      dimension R(n6),dr(n6)
      dimension Pa(n6),Pb(n6),P(n6),VE1(n6),VE2(n6),V1(n6),V2(n6)

      DATA CA/-.487635d0,.460931d0,-1.093472d0,.113334d0,-.170666d0,&
             1.043245d0,-.137693d0,.586113d0,-.678741d0,-.242872d0,&
             .028827d0,-.249638d0,.212290d0,-.257605d0,.235806d0,&
             .082524d0,-.247365d0,.414443d0,-.518085d0,.622019d0,&
             -.216965d0,-.890872d0,-.044282d0,-.221662d0,.159334d0,&
             -.377969d0,.108116d0,.240714d0,-.085919d0,.032833d0/

      DATA CB/ 0.0393d0,-0.0035d0,0.0711d0,0.0269d0,0.0078d0,0.1452d0,&
              0.0097d0,-0.0318d0,0.0234d0,0.0305d0/

      DATA C/ 0.3192d0,-0.7084d0,-0.2253d0,0.9527d0,0.1778d0,0.0597d0,&
             -0.5182d0,-0.1000d0,0.1582d0,0.0426d0/


      DATA DeA,Re1,a1,a2,a3/ 4.73379d0,0.7408d0,&
                            4.38128657d0,6.12491079d0,5.99581549d0/

      DATA DeB,Re2,a4,a5,a6/ 2.792820072d0,1.05835d0,&
                           3.09142103d0,3.4386363d0,2.70834046d0/

       DATA V0A,R0A,g1A/3.782643d0,0.89175599d0,0.820040d0/

      DATA V0B,R0B,g1B/10.2785d0,1.5358d0,1.5027d0/

      DATA V0,R0,g1/0.3878d0,1.1783d0,2.1492d0/
      DATA dpi/3.1415927d0/

!c*****************************
!    CONVERSION  Bohr to angstrom
!c
        q(1)=xx1*0.52917706d0
        q(2)=xx2*0.52917706d0
        q(3)=xx3
        
              CA( 1) =   -.571408    
              CA( 2) =    .506482   
              CA( 3) =  -1.390376  
              CA( 4) =   -.067676 
              CA( 5) =    .269395
              CA( 6) =    .506843    
              CA( 7) =  -1.112170   
              CA( 8) =   1.367337  
              CA( 9) =   -.505933 
              CA(10) =   -.293975    
              CA(11) =   1.067224   
              CA(12) =   -.703635  
              CA(13) =   -.361988 
              CA(14) =   -.794325    
              CA(15) =   -.077245   
              CA(16) =   -.399357  
              CA(17) =    .097876 
              CA(18) =   1.005771    
              CA(19) =  -1.435035   
              CA(20) =   -.500318  
              CA(21) =   1.142223 
              CA(22) =  -3.552668    
              CA(23) =    .013631   
              CA(24) =   -.477215  
              CA(25) =    .176966 
              CA(26) =   -.352025
              CA(27) =    .289533    
              CA(28) =    .201800   
              CA(29) =   -.294577  
              CA(30) =    .170270 

!   CONVERSION TO BOND-LENGTH COORDINATES (R1,R2,R3)

      R(3)=dsqrt(q(2)**2+(0.5d0*q(1))**2+q(1)*q(2)*q(3))
      R(2)=dsqrt(q(2)**2+(0.5d0*q(1))**2-q(1)*q(2)*q(3))
      R(1)=q(1)

!c--------------------------------------------------------------
!c--------------------------------------------------------------
!C
!       TWO-BODY TERMS
!C
!c--------------------------------------------------------------
         sum1=0.0d0
         sum2=0.0d0
         do i=1,n6 
         PA(i)=R(i)-Re1
         PB(i)=R(i)-Re2
         VE1(i)=-DeA*DEXP(-a1*PA(i))
         VE2(i)=-DeB*DEXP(-a4*PB(i))
         V1(i)=(1.0D0+PA(i)*(a1+PA(i)*(a2+PA(i)*a3)))*VE1(i)
         V2(i)=(1.0D0+PB(i)*(a4+PB(i)*(a5+PB(i)*a6)))*VE2(i)
         sum1=sum1+V1(i)
         sum2=sum2+V2(i)
         enddo
         Vaa2=sum1
         Vbb2=sum2
!c----------------------------------------------------------------
!c-------------------------------------------------------------------
!            THREE-BODY TERMS
!C--------------------------------------------------------------------------
            do i=1,n6
               Pa(i)=R(i)-R0A
               Pb(i)=R(i)-R0B
               P(i)=R(i)-R0
            enddo
!c---------------------------------------------------------------------
            prodA=1.d0
            prodB=1.d0
            prod=1.d0
            do i=1,n6
              prodA=prodA*(1.0d0-dtanh(0.5d0*g1A*Pa(i)))
              prodB=prodB*(1.0d0-dtanh(0.5d0*g1B*Pb(i)))
              prod=prod*(1.0d0-dtanh(0.5d0*g1*P(i)))
            enddo
            t1a=prodA
            t1b=prodB
            t1=prod
!C---------------------------------------------------------------------------
            Vaa3=(1.0D0+CA(1)*(Pa(1)+Pa(2)+Pa(3)) +CA(2)*(Pa(1)**2+Pa(2)**2+Pa(3)**2)&
           +CA(3)*(Pa(1)*Pa(2)+Pa(1)*Pa(3)+Pa(2)*Pa(3)) +CA(4)*(Pa(1)**3+Pa(2)**3&
           +Pa(3)**3)  +CA(5)*(Pa(1)**2*Pa(3)+Pa(1)**2*Pa(2)+Pa(1)*Pa(3)**2&
           +Pa(1)*Pa(2)**2+Pa(3)**2*Pa(2)+Pa(3)*Pa(2)**2)+CA(6)*Pa(1)*Pa(2)*Pa(3)&
           +CA(7)*(Pa(1)**4+Pa(2)**4+Pa(3)**4)+CA(8)*(Pa(1)**3*Pa(3)+Pa(1)*Pa(3)**3+Pa(1)**3*Pa(2)&
           +Pa(3)*Pa(2)**3+Pa(3)**3*Pa(2)+Pa(1)*Pa(2)**3)&
       +CA(9)*(Pa(1)**2*Pa(3)**2+Pa(1)**2*Pa(2)**2+Pa(2)**2*Pa(3)**2) +CA(10)*(Pa(1)**2*Pa(2)*Pa(3)+Pa(3)**2*Pa(2)*Pa(1)&
           +Pa(1)*Pa(3)*Pa(2)**2) +CA(11)*(Pa(1)**5+Pa(2)**5+Pa(3)**5)&
       +CA(12)*(Pa(1)**4*Pa(2)+Pa(1)**4*Pa(3)+Pa(1)*Pa(3)**4 +Pa(3)*Pa(2)**4+Pa(3)**4*Pa(2)+Pa(1)*Pa(2)**4)&
       +CA(13)*(Pa(1)**3*Pa(2)**2+Pa(1)**2*Pa(2)**3+Pa(1)**3*Pa(3)**2 +Pa(1)**2*Pa(3)**3+Pa(2)**3*Pa(3)**2+Pa(2)**2*Pa(3)**3)&
       +CA(14)*(Pa(1)**3*Pa(2)*Pa(3)+Pa(1)*Pa(2)*Pa(3)**3  +Pa(1)*Pa(2)**3*Pa(3))&
       +CA(15)*(Pa(1)**2*Pa(2)**2*Pa(3)+Pa(1)**2*Pa(2)*Pa(3)**2 +Pa(1)*Pa(2)**2*Pa(3)**2)&
       +CA(16)*(Pa(1)**6+Pa(2)**6+Pa(3)**6) +CA(17)*(Pa(1)*Pa(2)**5+Pa(1)*Pa(3)**5+Pa(1)**5*Pa(2)&
           +Pa(1)**5*Pa(3)+Pa(2)*Pa(3)**5+Pa(2)**5*Pa(3)) +CA(18)*(Pa(1)**4*Pa(2)**2+Pa(1)**4*Pa(3)**2&
           +Pa(1)**2*Pa(3)**4 +Pa(1)**2*Pa(2)**4+Pa(2)**4*Pa(3)**2  +Pa(2)**2*Pa(3)**4)&
       +CA(19)*(Pa(1)**3*Pa(2)**3+Pa(1)**3*Pa(3)**3&
           +Pa(2)**3*Pa(3)**3) +CA(20)*(Pa(1)**4*Pa(2)*Pa(3)+Pa(1)*Pa(2)**4*Pa(3)&
           +Pa(1)*Pa(2)*Pa(3)**4) +CA(21)*(Pa(1)**3*Pa(2)**2*Pa(3)+Pa(1)**3*Pa(2)*Pa(3)**2&
           +Pa(1)**2*Pa(2)*Pa(3)**3+Pa(1)**2*Pa(2)**3*Pa(3) +Pa(1)*Pa(2)**3*Pa(3)**2+Pa(1)*Pa(2)**2*Pa(3)**3)&
       +CA(22)*(Pa(1)**2*Pa(2)**2*Pa(3)**2) +CA(23)*(Pa(1)**7+Pa(2)**7+Pa(3)**7)&
       +CA(24)*(Pa(1)**2*Pa(2)**5+Pa(1)**2*Pa(3)**5+Pa(1)**5*Pa(2)**2&
           +Pa(1)**5*Pa(3)**2+Pa(2)**2*Pa(3)**5+Pa(2)**5*Pa(3)**2) +CA(25)*(Pa(1)**6*Pa(2)+Pa(1)**6*Pa(3)+Pa(1)*Pa(2)**6&
           +Pa(1)*Pa(3)**6+Pa(2)**6*Pa(3)+Pa(2)*Pa(3)**6) +CA(26)*(Pa(1)**5*Pa(2)*Pa(3)+Pa(1)*Pa(2)*Pa(3)**5&
           +Pa(1)*Pa(2)**5*Pa(3)) +CA(27)*(Pa(1)**4*Pa(2)**3+Pa(1)**4*Pa(3)**3&
           +Pa(1)**3*Pa(3)**4+Pa(1)**3*Pa(2)**4  +Pa(2)**4*Pa(3)**3+Pa(2)**3*Pa(3)**4)&
       +CA(28)*(Pa(1)**4*Pa(2)**2*Pa(3)+Pa(1)**4*Pa(2)*Pa(3)**2&
           +Pa(1)**2*Pa(2)**4*Pa(3)+Pa(1)**2*Pa(2)*Pa(3)**4 +Pa(1)*Pa(2)**4*Pa(3)**2+Pa(1)*Pa(2)**2*Pa(3)**4)&
       +CA(29)*(Pa(1)**3*Pa(2)*Pa(3)**3+Pa(1)**3*Pa(2)**3*Pa(3) +Pa(1)*Pa(2)**3*Pa(3)**3)&
       +CA(30)*(Pa(1)**2*Pa(2)**3*Pa(3)**2 +Pa(1)**3*Pa(2)**2*Pa(3)**2+Pa(1)**2*Pa(2)**2*Pa(3)**3)&
       )*V0A*t1a
!c-----------------------------------------------------------------------
!C---------------------------------------------------------------------------
            Vbb3=(1.0D0+CB(1)*(Pb(1)+Pb(2)+Pb(3))&
           +CB(2)*(Pb(1)**2+Pb(2)**2+Pb(3)**2)&
           +CB(3)*(Pb(1)*Pb(2)+Pb(1)*Pb(3)+Pb(2)*Pb(3))&
           +CB(4)*(Pb(1)**3+Pb(2)**3+Pb(3)**3)&
           +CB(5)*(Pb(1)**2*Pb(3)+Pb(1)**2*Pb(2)+Pb(1)*Pb(3)**2&
           +Pb(1)*Pb(2)**2+Pb(3)**2*Pb(2)+Pb(3)*Pb(2)**2)&
           +CB(6)*Pb(1)*Pb(2)*Pb(3)&
           +CB(7)*(Pb(1)**4+Pb(2)**4+Pb(3)**4)&
           +CB(8)*(Pb(1)**3*Pb(3)+Pb(1)*Pb(3)**3+Pb(1)**3*Pb(2)&
           +Pb(3)*Pb(2)**3+Pb(3)**3*Pb(2)+Pb(1)*Pb(2)**3)&
      +CB(9)*(Pb(1)**2*Pb(3)**2+Pb(1)**2*Pb(2)**2+Pb(2)**2*Pb(3)**2)&
           +CB(10)*(Pb(1)**2*Pb(2)*Pb(3)+Pb(3)**2*Pb(2)*Pb(1)&
           +Pb(1)*Pb(3)*Pb(2)**2))*V0B*t1b
!c-----------------------------------------------------------------------
!C---------------------------------------------------------------------------
            Vab=(1.0D0+C(1)*(P(1)+P(2)+P(3))&
           +C(2)*(P(1)**2+P(2)**2+P(3)**2)&
           +C(3)*(P(1)*P(2)+P(1)*P(3)+P(2)*P(3))&
           +C(4)*(P(1)**3+P(2)**3+P(3)**3)&
           +C(5)*(P(1)**2*P(3)+P(1)**2*P(2)+P(1)*P(3)**2&
           +P(1)*P(2)**2+P(3)**2*P(2)+P(3)*P(2)**2)&
           +C(6)*P(1)*P(2)*P(3)&
           +C(7)*(P(1)**4+P(2)**4+P(3)**4)&
           +C(8)*(P(1)**3*P(3)+P(1)*P(3)**3+P(1)**3*P(2)&
           +P(3)*P(2)**3+P(3)**3*P(2)+P(1)*P(2)**3)&
          +C(9)*(P(1)**2*P(3)**2+P(1)**2*P(2)**2+P(2)**2*P(3)**2)&
           +C(10)*(P(1)**2*P(2)*P(3)+P(3)**2*P(2)*P(1)&
           +P(1)*P(3)*P(2)**2))*V0*t1
!c---------------------------------------------------------------------
!c---------------------------------------------------------------------
!           H+ - H2 potential term (VPTH) for large 
!           separations R>/10.a0
!c-----------------------------------------------------------------------
!C I start my symmetrisation here
             dminr = min (R(1),R(2),R(3))
                 q11  = dsqrt((R(3)**2+R(2)**2)*0.5d0-R(1)**2/4.0d0)
                 q12  = dsqrt((R(1)**2+R(3)**2)*0.5d0-R(2)**2/4.0d0)
                 q13  = dsqrt((R(2)**2+R(1)**2)*0.5d0-R(3)**2/4.0d0)
             if(R(1).eq.dminr) then
                     q(2) = q11
                    q(1) = R(1)
                    q(3) = (R(3)**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    go to 1957
                    endif
             if(R(2).eq.dminr) then
                      q(2) = q12
                        q111=R(1)
                        q112=R(2)
                        q113=R(3)
                    q(1) =q112 
                    q(3) = (q111**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    go to 1957
                    endif
             if(R(3).eq.dminr) then
                    q(2) = q13
                        q111=R(1)
                        q112=R(2)
                        q113=R(3)
                    q(1) =q113 
                    q(3) = (q112**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    endif
 1957              continue    
             xxx1=q(1)/0.52917706d0
             xxx2=q(2)/0.52917706d0
             xxx3=q(3)                 
             RLIM=4.d0*0.52917706d0
             RLIM=4.5d0*0.52917706d0
             RMAX=10.d0*0.52917706d0
             rmmrlim = RMAX - RLIM
             if (dabs(q(2))  .gt. RLIM)  then 
                dr(2)=q(2)-RLIM

                sf1=(dcos(0.5d0*dpi*dr(2)/rmmrlim))**2
                sf2=1.d0-sf1

               Vpth=VLIM(xxx1,xxx2,xxx3)

!         Convert   hartree to eV
               Vpth=Vpth/0.036749d0 
!c
               Vaa3=sf2*Vpth+sf1*Vaa3
               if (q(2)  .gt. RMAX)            Vaa3=Vpth 
            endif
!c------------------------------------------------------
!c------------------------------------------------------
!           H - H2+ potential term (VPTH) for large 
!           separations r>/11.a0
!c-----------------------------------------------------------------------
             RLIM1=5.d0*0.52917706d0
             RMAX1=11.d0*0.52917706d0
             rmmrlim1 = RMAX1 - RLIM1
             if (dabs(q(1)) .gt. RLIM1)  then                   

               Vpth1=VLIM1(xxx1)
!         Convert   hartree to eV
               Vpth1=Vpth1/0.036749d0

                dr(1)=q(1)-RLIM1

                sf1=(dcos(0.5d0*dpi*dr(1)/rmmrlim1))**2
                sf2=1.d0-sf1

               Vbb3=sf2*Vpth1+sf1*Vbb3

               if (q(1)  .gt. RMAX1)            Vbb3=Vpth1
            endif
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
            Vaa=Vaa2+Vaa3
            Vbb=Vbb2+Vbb3
!C-----------------------------------------------------------------------
            Vh3p=0.5d0*(Vaa+Vbb-dsqrt((Vaa-Vbb)**2+4.d0*Vab**2))
!C----------------------------------------------------------------------
!         SHIFT and CONVERTION  eV to HARTREE
!c
          Vh3p=(Vh3p+9.342357355d0)*0.036749d0
!C
      return
      end

      SUBROUTINE SOTV(V,R1,R2,XCOS)    

!     TRANSFORM GENERALISED COORDINATES TO THOSE FOR PARTICULAR
!     SYSTEM. THIS VERSION TRANSFORMS TO AB2 BONDLENGTH-BONDANGLE
!     COORDINATES. ALLOWANCE MUST BE MADE FOR THE NUMBERING OF THE ATOMS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!jayesh 7/10/2002
!      COMMON /MASS/ XMASS(3),G1,G2
      common /mass/ xmass(3),g1,g2,xmassr(3) !consistant with DVR3DRJZ
!     (R = r . S = r'. T = theta)

      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/
      common /potential/ cv(100), ifit(100), der(100)

      IF (G1 .EQ. X0) THEN
!        BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
         Q1 = R1
         Q2 = R2
         Q3= SQRT(R1*R1 + R2*R2 - X2*R1*R2*XCOS)
      ELSE IF (G2 .EQ. X0) THEN
!        SCATTERING COORDINATES: ATOM 2 = ATOM 3
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
            Q1 = ABS(XX - R2)
            Q2 = (YY + R2)
         ELSE IF (XCOS .LE. (TINY - X1)) THEN
            Q1 = (XX + R2)
            Q2 = ABS(YY + R2)
         ELSE
            Q1 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
            Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
         ENDIF
         Q3 = R1
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
      ENDIF


      CALL POTS(V,Q1,Q2,Q3)     
!C
      RETURN
      END

      SUBROUTINE POTS(V,P1,P2,P3)   

!     H3+ POTENTIAL IN 7TH ORDER Morse fit to experimental data 
!     based on the ab initio potential of:
!     R.Rohse,W.Kutzelnigg,R.Jaquet,W.Klopper,in press (1994)
!     UNITS: HARTREE & BOHR
!     asymmetric part

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FT(200)
      common /potential/ cv(100), ifit(100), der(100)
!
      DATA RE/1.650D0/,BET/1.300D0/
      data nv/ 98/,nfmax/12/
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/

      SQ3=SQRT(THREE)
      SQ2=SQRT(TWO)
      FACTOR=BET/RE
      DR1= (P1-RE)
      DR2= (P2-RE)
      DR3= (P3-RE)
       Y1=(ONE-EXP(-FACTOR*DR1))/BET
       Y2=(ONE-EXP(-FACTOR*DR2))/BET
       Y3=(ONE-EXP(-FACTOR*DR3))/BET
      SA=(Y1+Y2+Y3)/SQ3
      SX=(Y3+Y3-Y1-Y2)/(SQ2*SQ3)
      SY=(Y2-Y1)/SQ2
      SX1=(DR3+DR3-DR1-DR2)/(SQ2*SQ3)
      SY1=(DR2-DR1)/SQ2
      QUAD1=SX1**2+SY1**2
      SE1=SQRT(QUAD1)
      phi=acos(sx1/se1)
      npot=1
      ft(1)=one
!c nfmax = 12    
      do 120 norder=1,nfmax
      do 110 n=norder,0,-1
      do 100 k=0,norder-n,3
      m=norder-k-n
      if (mod(m,2) .ne. 0) goto 100
      npot=npot+1
!     ft(npot)=sa**n * se**(m+k) * cos(dble(k)*phi)
      ft(npot)=sa**n * se1**(m+k) * cos(dble(k)*phi)
  100 continue
  110 continue
  120 continue
      V=ZERO
      DO 40 I=1,NV
      V=V+CV(I)*FT(I)
   40 continue 
      DO 103 i=1,100
      der(i)=ft(i)    
  103 continue
      RETURN
      END
      SUBROUTINE POTV11(Vh3p,xx1,xx2,xx3,v22) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(n6=3)
      dimension q(n6)
      dimension R(n6),dr(n6)
      dimension Pa(n6),VE1(n6),V1(n6)
!C
      DATA DeA,Re1,a1,a2,a3/ 4.73379d0,0.7408d0,&
                            4.38128657d0,6.12491079d0,5.99581549d0/
!c
      DATA dpi/3.1415927d0/
!c*****************************
!    CONVERSION  Bohr to angstrom
!c
        q(1)=xx1*0.52917706d0
        q(2)=xx2*0.52917706d0
        q(3)=xx3
!C
!   CONVERSION TO BOND-LENGTH COORDINATES (R1,R2,R3)
!c
      R(3)=dsqrt(q(2)**2+(0.5d0*q(1))**2+q(1)*q(2)*q(3))
      R(2)=dsqrt(q(2)**2+(0.5d0*q(1))**2-q(1)*q(2)*q(3))
      R(1)=q(1)
             dminr = min (R(1),R(2),R(3))
                 q11  = dsqrt((R(3)**2+R(2)**2)*0.5d0-R(1)**2/4.0d0)
                 q12  = dsqrt((R(1)**2+R(3)**2)*0.5d0-R(2)**2/4.0d0)
                 q13  = dsqrt((R(2)**2+R(1)**2)*0.5d0-R(3)**2/4.0d0)
             if(R(1).eq.dminr) then
                     q(2) = q11
                    q(1) = R(1)
                    q(3) = (R(3)**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    go to 1957
                    endif
             if(R(2).eq.dminr) then
                      q(2) = q12
                        q111=R(1)
                        q112=R(2)
                        q113=R(3)
                   R(1) = q112 
                    q(1) = R(1)
                    q(3) = (q111**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    go to 1957
                    endif
             if(R(3).eq.dminr) then
                    q(2) = q13
                        q111=R(1)
                        q112=R(2)
                        q113=R(3)
                   R(1) = q113 
                    q(1) = R(1)
                    q(3) = (q112**2-q(1)**2/4.0d0-q(2)**2)/(q(1)*q(2))
                    endif
 1957              continue    
             xxx1=q(1)/0.52917706d0
             xxx2=q(2)/0.52917706d0
             xxx3=q(3)                 
!c--------------------------------------------------------------
!c---------------------------------------------------------------------
!c---------------------------------------------------------------------
!!           H+ - H2 potential term (VPTH) for large 
!           separations R>/10.a0
!c-----------------------------------------------------------------------
             RLIM=4.d0*0.52917706d0
             RMAX=7.d0*0.52917706d0
                     rmmrlim = RMAX - RLIM
                     sf1 =1.0d0
                     sf2 =0.0d0
                 Vaa3 =     v22         
             if (dabs(q(2))  .gt. RLIM)  then 
                dr(2)=q(2)-RLIM
!               sf1=(dcos(0.5d0*dr(2)))**2
                sf1=(dcos(0.5d0*dpi*dr(2)/rmmrlim))**2
                sf2=1.d0-sf1
!c
!c--------------------------------------------------------------
!c--------------------------------------------------------------

!       TWO-BODY TERMS
!C
!c--------------------------------------------------------------
         sum1=0.0d0
!c!       do i=1,3  
            i=1    
         PA(i)=R(i)-Re1
         VE1(i)=-DeA*DEXP(-a1*PA(i))
         V1(i)=(1.0D0+PA(i)*(a1+PA(i)*(a2+PA(i)*a3)))*VE1(i)
         sum1=sum1+V1(i)
!c!       enddo
         Vaa2=sum1
               Vpth=VLIM(xxx1,xxx2,xxx3)
!         Convert   hartree to eV
               Vpth=Vpth/0.036749d0 
                 deconst = (37170.0+38180.0)/219474.6354d0 /0.036749d0
                 Vpth = deconst + Vaa2 +Vpth
!         Convert  ev  to hartree 
            Vpth= Vpth*0.036749d0      

               Vaa3=sf2*Vpth+sf1*V22 

               if (q(2)  .gt. RMAX)            Vaa3=Vpth 
             endif

            Vh3p= Vaa3                 
      return
      end
      
!*******************************************************************      

      SUBROUTINE RJ(R1r,R2r,XCOSr,R1j,R2j,XCOSj)
!Jayesh 8/10/2002
!This subroutine converts between Radua and Jacobi coordinates

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

!jayesh 7/10/2002
!      COMMON /MASS/ XMASS(3),G1,G2
      common /mass/ xmass(3),g1,g2,xmassr(3) !consistant with DVR3DJRZ

      DATA X0/0.0D0/,X1/1.0D0/,X2/2.0D0/,TINY/9.0D-15/

!If the cooordninates are already Jacobi then psudo Radua coordintaes are
!simply copied to the correct JAocobi variables. Q1, Q2, and Q3 are then 
!calculated for the particular coordinate system
IF (G1 .EQ. X0) THEN
!  BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
   XCOSj = XCOSr
   R1j = R1r
   R2j = R2r

   Q1 = R1j
   Q2 = R2j
   Q3= DSQRT(R1j*R1j + R2j*R2j - X2*R1j*R2j*XCOSj)

ELSEIF (G2.EQ.0.0) THEN
!  SCATTERING COORDINATES: ATOM 2 = ATOM 3
   XCOSj = XCOSr
   R1j = R1r
   R2j = R2r

   XX = R1j * G1
   YY = R1j * (X1 - G1)
   IF (R2j .EQ. X0 .OR. XCOSj .GE. (X1 - TINY)) THEN
      Q1 = DABS(XX - R2j)
      Q2 = (YY + R2j)
   ELSE IF (XCOSj .LE. (TINY - X1)) THEN
      Q1 = (XX + R2j)
      Q2 = DABS(YY + R2j)
   ELSE
      Q1 = DSQRT(XX*XX + R2j*R2j - X2*XX*R2j*XCOSj)
      Q2 = DSQRT(YY*YY + R2j*R2j + X2*YY*R2j*XCOSj)
   ENDIF
   Q3 = R1j


!Else proceed with the conversion to Jacobi
ELSE
   !     GENERAL COORDINATES (INCLUDING RADAU): ATOM 1 = ATOM 2
   F1= X1/G1
   F2= X1/G2
   F12= X1 - F1*F2
   P1= R1r*(X1-F1)/(G2*F12)
   P2= R2r*(X1-F2)/(G1*F12)
   S1= R1r-P1
   S2= R2r-P2
   Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOSr)/(X1-G1)
   Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOSr)/(X1-G2)
   Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOSr)         

   !     Conversion Q1,Q2,Q3 in Jacobi coordinats                  
   R1j=Q3
   R2j=DSQRT((Q2**2+Q1**2)*0.5d0-Q3**2/4.0d0)      
   XCOSj=((Q1**2-Q2**2)/(R2j*Q3))*0.5d0
endif

      return
    END SUBROUTINE RJ

      BLOCK DATA POT_POINTS
      !Jayesh
      !this is a subroutine to read in the potential points
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)      
      common /potential/ cv(100), ifit(100), der(100)
      
      DATA cv/-0.087789852d0, 0.98045682d0, 44864.41269d0, 21476.84833d0, &
           -10851.63608d0, -51590.2709d0, -7204.4524d0, 5829.603497d0,  &
           52139.88524d0, 21815.41705d0, -169.8098455d0, 637.9898341d0, & 
           -32129.45288d0, -28540.41255d0, -2878.831201d0, 363.2206034d0, &
           -2329.206071d0, 15281.97447d0, 21353.68336d0, 9764.910178d0, & 
           399.6016332d0, 56.90514503d0, 4.871752983d0, -8406.165782d0, & 
           -4778.100604d0, -12375.47595d0, -13632.50888d0, -3939.977682d0, & 
           -501.9498256d0, -257.2286076d0, -10.65931215d0, 16786.79531d0, &
           6419.982177d0, 14318.43234d0, 14420.94969d0, 9337.455963d0, &
           1124.258252d0, 866.4016246d0, -37.10320237d0, -26.3356168d0, &
           1.808227334d0, 20981.73391d0, -7708.894559d0, -12335.83598d0, & 
           -16641.8539d0, -12336.83354d0, -802.6889977d0, -1360.451838d0, &
           114.1768807d0, 303.0208632d0, 123.5641465d0, 13.70748335d0, & 
           0.085017887d0, -37925.14149d0, -8082.294389d0, -11659.72527d0, & 
           -3844.084482d0, -5892.174876d0, -1023.092602d0, 1231.687805d0, &
           950.0774708d0, -522.5815229d0, -388.9429905d0, -9.834074368d0, &
           2.756489421d0, -2.363427916d0, -5.380949488d0, -17044.04289d0, &
           6625.762365d0, 19679.34106d0, 36270.14125d0, 32713.47521d0, &
           5834.930375d0, -1089.532875d0, -783.2395159d0, -104.7708236d0, &
           546.8593236d0, -273.7888178d0, -16.56040114d0, -20.8278439d0, &
           -10.33875408d0, -2.580694632d0, 0.013348447d0, 29519.49476d0,&
           5046.99291d0, -4009.477109d0, -23426.87491d0, -21373.81002d0, &
           -5816.439009d0, 790.8732935d0, -672.0810895d0, 442.1395696d0, &
           -329.8193608d0, 340.1812585d0, 15.70723083d0, 30.39865845d0, &
           22.11444588d0, 3.282330882d0,0.0d0, 0.0d0/

       DATA ifit/ 100 * 1/

   END
