      SUBROUTINE potv(V,R1,R2,xcos)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y),logical(z)
      COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)
      dimension rij(3)
      real*8 :: xp(1:300)
      integer :: i_t
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/

      pi=dacos(-1.d0)

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

c Breit surface for H2S.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers

      !
      att1=Q1*0.5291772d0
      att2=Q2*0.5291772d0

      ! call pbbreit(Vpbb,att1,att2,THETA)

cc  H2S lamb corr at the cc level
cc input in A please
cc conversion 0.529177208d0
cc output in wavenumbers

      ! call lambshift(Vlamb,att1,att2,THETA)

c This sub calculates the correction to the h2s BO pes given
c  by the 2 el darwin term.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers

      ! call d2h2s(vdarw,att1,att2,THETA)

c Gaunt surface for H2S.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers

      ! call pbgaunt(Vgaunt,att1,att2,THETA)
      !
      open(unit=32,status='old',form='formatted',file='pot.fit')
      read(32,*) npropin
      do i = 1,npropin
         read(32,*) i_t, xp(i)
      end do
      close(  unit=32)
      !
      if (npropin>size(xp)) then 
        !
        write(6,"('wifin: Too many parameters in pot.fit: ',
     .                    i8,' max = ',i8)") npropin,size(xp)
        stop 'wifin: Too many parameters in pot.fit'
        !
      endif
	!
      !write(6,"(' npropin= ',i8,' max = ',i8)") npropin,size(xp)
	!   
      call poten(xp(1:npropin),vp,att1,att2,THETA)
	!
	!call morbid(att1,att2,THETA,vp)

! This is where the ab initio potential gets defined from the bits and pieces

       !v= (vp +Vpbb+Vlamb+Vgaunt)/219474.624d0
       !v= (vp +Vpbb+Vlamb+vdarw+Vgaunt)/219474.624d0
       v= (vp )/219474.624d0
	!
	!write(6,"('  att1,att2,theta =  ',4f18.4,'  vp = ',f18.8)") 
        !      att1,att2,thet,att1-att2,vp
	!
      end




        SUBROUTINE lambshift(V,r1,r2,z)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

cc  H2S lamb corr at the cc level
cc input in A please
cc conversion 0.529177208d0
cc output in wavenumbers


         ure=0.373851429340111951d0
         dcoste=0.347531934855387636d-1



         xs1=(r1+r2)*ure-1.d0
         xs2=(r1-r2)*ure
         xst=dcos(z)-dcoste

        v=0.3517749424262E+01
     y +0.7967585292444E-03*xst
     y +0.4319615574716E-02*xs1
     y -0.2691250948929E-03*xst**2
     y +0.6586180071120E-03*xs2**2
     y -0.1111735865397E-02*xs1*xst
      v=v+0.2153045911685E-02*xs1**2
     y -0.2690481594546E-02*xst**3
     y -0.2176104087520E-02*xs2**2*xst
     y -0.6974363197770E-03*xs1*xst**2
     y -0.3689966438700E-01*xs1*xs2**2
     y -0.6738725992695E-02*xs1**2*xst
     y -0.1028825525410E-01*xs1**3
     y -0.1353403356769E-02*xst**4
     y -0.2204096713409E-02*xs2**2*xst**2
     y +0.2311483500271E-01*xs2**4
     y +0.6127987300035E-03*xs1*xst**3
     y -0.4858215273930E-02*xs1*xs2**2*xst
     y +0.7514851776139E-03*xs1**2*xst**2
     y +0.1232767777430E+00*xs1**2*xs2**2
     y +0.1628864390992E-01*xs1**3*xst
      v=v+0.2127842392300E-01*xs1**4
     y +0.1557297060788E-01*xst**5
     y -0.1291865163824E-02*xs2**2*xst**3
     y +0.1060606854136E-01*xs1*xst**4
     y -0.6337866495011E-02*xs1*xs2**2*xst**2
     y -0.1343882537761E+00*xs1*xs2**4
     y +0.2175639476592E-01*xs1**2*xst**3
     y +0.1217569172069E-01*xs1**3*xst**2
     y -0.3598492930047E+00*xs1**3*xs2**2
     y -0.1103993647854E+00*xs1**5
     y +0.1146992797608E-01*xst**6
     y +0.6942618444501E-02*xs2**2*xst**4
     y +0.4698702582223E-01*xs2**4*xst**2
     y -0.1005362975968E-01*xs2**6
     y -0.1717030915047E-02*xs1*xst**5
      v=v+0.1928063645074E+00*xs1*xs2**4*xst
     y +0.8134466371036E-02*xs1**2*xst**4
     y +0.4113657301080E+00*xs1**2*xs2**4
     y -0.5435556165923E-01*xs1**3*xst**3
     y -0.3096551228015E-02*xs1**4*xst**2
     y +0.8977467213924E+00*xs1**4*xs2**2
      v=v+0.1870518393999E+00*xs1**6
     y -0.1815202422622E-01*xst**7
     y -0.1098013405862E-01*xs1*xst**6
     y -0.1402280716068E+00*xs1*xs2**4*xst**2
     y -0.3116400901348E-01*xs1**2*xst**5
     y -0.6714023050353E+00*xs1**2*xs2**4*xst
     y -0.9188732367685E-01*xs1**3*xst**4
     y -0.4786248606952E+00*xs1**3*xs2**4
      v=v+0.4014818468875E-01*xs1**4*xst**3
     y +0.4433256971585E+00*xs1**4*xs2**2*xst
     y -0.9503878816037E+00*xs1**5*xs2**2
     y -0.8062699132412E-01*xs1**6*xst
     y -0.1672284004348E-01*xst**8

      v=v*0.01726d0*219474.624d0

      RETURN
      END

        SUBROUTINE pbbreit(V,r1,r2,z)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c Breit surface for H2S.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers


         ure=0.373851429340111951d0
         dcoste=0.347531934855387636d-1

         xs1=(r1+r2)*ure-1.d0
         xs2=(r1-r2)*ure
         xst=dcos(z)-dcoste


         v=0.8765756926097E-01
     y   -0.7470466984977E-06*xst
     y   -0.4591567724873E-03*xs1
     y   +0.2438277378851E-04*xst**2
     y   +0.1059524622662E-02*xs2**2
     y   -0.3189553736298E-05*xs1*xst
     y   +0.1089138001148E-02*xs1**2
     y   +0.7735193302669E-04*xst**3
     y   +0.1406958041130E-04*xs2**2*xst
     y   +0.1630810111481E-04*xs1*xst**2
     y   -0.4542479540109E-02*xs1*xs2**2
     y   +0.1002856632993E-04*xs1**2*xst
     y   -0.1520879954684E-02*xs1**3
          v=v+0.1391996175182E-04*xst**4
     y   -0.4717723013451E-05*xs2**2*xst**2
     y   +0.1558678228909E-02*xs2**4
     y   +0.1437254575428E-04*xs1*xst**3
     y   +0.4118703616626E-05*xs1*xs2**2*xst
     y   -0.8390791525712E-04*xs1**2*xst**2
     y   +0.9097378783447E-02*xs1**2*xs2**2
     y   +0.6306768748460E-04*xs1**3*xst
         v=v+0.1455020958607E-02*xs1**4
     y   -0.3661914850753E-03*xst**5
     y   -0.9407748202265E-04*xs2**2*xst**3
     y   -0.1601893683893E-04*xs1*xst**4
     y   -0.1107885445634E-03*xs1*xs2**2*xst**2
     y   -0.6333559444331E-02*xs1*xs2**4
     y   -0.2190907906757E-03*xs1**2*xst**3
     y   -0.1851060535471E-03*xs1**3*xst**2
     y   -0.1229283119564E-01*xs1**3*xs2**2
     y   -0.6863742081354E-03*xs1**5
     y   -0.2304221502848E-03*xst**6
     y   -0.1158774739975E-03*xs2**2*xst**4
     y   -0.4196870186460E-03*xs2**4*xst**2
     y   +0.3119089521774E-03*xs2**6
        v=v+0.8782036319681E-04*xs1*xst**5
     y   -0.1283168422912E-03*xs1*xs2**4*xst
     y   -0.2998997605629E-04*xs1**2*xst**4
     y   +0.1394552871461E-01*xs1**2*xs2**4
     y   +0.2933772450811E-03*xs1**3*xst**3
     y   +0.4098561367544E-03*xs1**4*xst**2
     y   +0.1069104080091E-01*xs1**4*xs2**2
     y   -0.3326196046609E-03*xs1**6
     y   +0.4539164927449E-03*xst**7
     y   +0.9544310398796E-04*xs1*xst**6
     y   +0.2404606656808E-02*xs1*xs2**4*xst**2
     y   +0.1954265488403E-03*xs1**2*xst**5
     y   +0.1238735531947E-02*xs1**2*xs2**4*xst
     y   +0.7499245400255E-03*xs1**3*xst**4
     y   -0.1397979747135E-01*xs1**3*xs2**4
         v=v+0.6977006962489E-04*xs1**4*xst**3
     y   -0.2347201681002E-02*xs1**4*xs2**2*xst
     y   -0.2837645063101E-02*xs1**5*xs2**2
     y   -0.9643600951545E-03*xs1**6*xst
     y   +0.3831112710055E-03*xst**8

        V=v*219474.636d0

        RETURN
        END

       SUBROUTINE d2h2s(v,r1,r2,z)

      implicit real*8 (a-h,o-z)

c This sub calculates the correction to the h2s BO pes given
c  by the 2 el darwin term.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers

               ure=0.373851429340111951d0
         dcoste=0.347531934855387636d-1



         xs1=(r1+r2)*ure-1.d0
         xs2=(r1-r2)*ure
         xst=dcos(z)-dcoste




       v1=0.3288315669739d-01-0.7071887537443d-05*xst 
     y +0.4746750997890d-04*xs1+0.1814642649677d-05*xst**2 
     y -0.2069467427819d-03*xs2**2+0.7332242095409d-05*xs1*xst      
     y -0.2186179022213d-03*xs1**2+0.2480861522819d-04*xst**3      
     y +0.1737358816430d-04*xs2**2*xst+0.6343875806066d-05*xs1*xst**2      
     y +0.1214158968508d-02*xs1*xs2**2
     y +0.7174021871713d-04*xs1**2*xst
     y +0.4123102067716d-03*xs1**3+0.9304285225731d-05*xst**4
     y +0.2645987433615d-04*xs2**2*xst**2
     y -0.4694606895261d-03*xs2**4-0.4614784740426d-05*xs1*xst**3
     y +0.6735986785697d-04*xs1*xs2**2*xst
       v2=-0.1031568297693d-04*xs1**2*xst**2
     y -0.3018352737328d-02*xs1**2*xs2**2
     y -0.1687936203472d-03*xs1**3*xst-0.5407081511795d-03*xs1**4
     y -0.1486645300254d-03*xst**5
     y +0.2608702101747d-04*xs2**2*xst**3-0.1054721356589d-03*xs1*xst**4
     y +0.4850358295734d-04*xs1*xs2**2*xst**2
     y +0.2354824772719d-02*xs1*xs2**4-0.1994232364960d-03*xs1**2*xst**3
     y -0.1698790405568d-03*xs1**3*xst**2
     y +0.5107201201947d-02*xs1**3*xs2**2+0.7453795906214d-03*xs1**5
     y -0.1048145279829d-03*xst**6-0.5027878777653d-04*xs2**2*xst**4
     y -0.6135483118135d-03*xs2**4*xst**2-0.2786927682760d-03*xs2**6
     y +0.2042132380410d-05*xs1*xst**5
     y -0.2260425245176d-02*xs1*xs2**4*xst
     y -0.4759906023180d-04*xs1**2*xst**4
     y -0.6101027027456d-02*xs1**2*xs2**4
       v3=0.4650204217126d-03*xs1**3*xst**3
     y +0.2376486217449d-03*xs1**4*xst**2
     y -0.7056553518283d-02*xs1**4*xs2**2-0.6813879534093d-03*xs1**6
     y +0.1741840225940d-03*xst**7+0.9478904115519d-04*xs1*xst**6
     y +0.1668130278001d-02*xs1*xs2**4*xst**2
     y +0.3089302549319d-03*xs1**2*xst**5
     y +0.7991472940019d-02*xs1**2*xs2**4*xst
     y +0.8269964331662d-03*xs1**3*xst**4
     y +0.7854115056298d-02*xs1**3*xs2**4
     y -0.2205039469374d-03*xs1**4*xst**3
     y -0.4779484677740d-02*xs1**4*xs2**2*xst
     y +0.4917493052450d-02*xs1**5*xs2**2
     y +0.9271691674515d-03*xs1**6*xst+0.1574429558949d-03*xst**8


       v=(v1+v2+v3)*219474.636d0 

c      write(6,*)'finisco d2 sub'


       RETURN
       END
        SUBROUTINE pbgaunt(V,r1,r2,z)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)


c Gaunt surface for H2S.
c inputs in A (r1,r2) and rad (z)
c output in wavenumbers

         ure=0.373851429340111951d0
         dcoste=0.347531934855387636d-1

         xs1=(r1+r2)*ure-1.d0
         xs2=(r1-r2)*ure
         xst=dcos(z)-dcoste


           v=0.9453750817097E-01
     y   -0.3527497454934E-05*xst
     y   -0.5623398479453E-03*xs1
     y   +0.3184119943909E-04*xst**2
     y   +0.1264983806594E-02*xs2**2
     y   -0.6433822892762E-06*xs1*xst
     y   +0.1293621900960E-02*xs1**2
     y   +0.1153731226002E-03*xst**3
     y   +0.1944371960638E-04*xs2**2*xst
     y   +0.2062421277904E-04*xs1*xst**2
     y   -0.5392912681481E-02*xs1*xs2**2
     y   +0.4215969663932E-04*xs1**2*xst
     y   -0.1813768220418E-02*xs1**3
     y   +0.2231963765567E-04*xst**4
     y   -0.1179435303219E-04*xs2**2*xst**2
          v=v+0.1829273468087E-02*xs2**4
     y   +0.1372457037480E-04*xs1*xst**3
     y   +0.9512308901375E-05*xs1*xs2**2*xst
     y   -0.9855111329499E-04*xs1**2*xst**2
     y   +0.1075732737890E-01*xs1**2*xs2**2
     y   +0.5452750165331E-06*xs1**3*xst
     y   +0.1720769005842E-02*xs1**4
     y   -0.5601955472518E-03*xst**5
     y   -0.1114787862952E-03*xs2**2*xst**3
     y   -0.5837076366327E-04*xs1*xst**4
     y   -0.6509947703250E-04*xs1*xs2**2*xst**2
     y   -0.7017410490921E-02*xs1*xs2**4
     y   -0.3784212684319E-03*xs1**2*xst**3
         v=v-0.2768081471488E-03*xs1**3*xst**2
     y   -0.1342821012540E-01*xs1**3*xs2**2
     y   -0.4023255253942E-03*xs1**5
     y   -0.3576631673448E-03*xst**6
     y   -0.1510735422018E-03*xs2**2*xst**4
     y   -0.5715755491279E-03*xs2**4*xst**2
     y   +0.3721121329871E-03*xs2**6
     y   +0.1125331380038E-03*xs1*xst**5
     y   -0.5909045823656E-03*xs1*xs2**4*xst
     y   -0.7479511751489E-04*xs1**2*xst**4
     y   +0.1414703531371E-01*xs1**2*xs2**4
     y   +0.5764206029177E-03*xs1**3*xst**3
     y   +0.4784970721672E-03*xs1**4*xst**2
     y   +0.8571116522285E-02*xs1**4*xs2**2
     y   -0.1168581477947E-02*xs1**6
     y   +0.6891550629468E-03*xst**7
        v=v+0.1556515276934E-03*xs1*xst**6
     y   +0.2905839729060E-02*xs1*xs2**4*xst**2
     y   +0.3630809675061E-03*xs1**2*xst**5
     y   +0.3296097239931E-02*xs1**2*xs2**4*xst
     y   +0.1318009995887E-02*xs1**3*xst**4
     y   -0.1360157387589E-01*xs1**3*xs2**4
     y   +0.1277729109407E-03*xs1**4*xst**3
     y   -0.3918380101742E-02*xs1**4*xs2**2*xst
     y   +0.1147719312132E-02*xs1**5*xs2**2
     y   -0.1036645894718E-02*xs1**6*xst
     y   +0.5839982090765E-03*xst**8


         V=v*219474.636d0

         RETURN
         END




      SUBROUTINE morbid(q1,q2,th,v)

!
! In here almost all the points have been included, 342 in total
!
!
      implicit real*8 (a-h,o-z)
      double precision r1,r2,th,v
      COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)



           !xp = potparam      

         r1=Q1 !*0.5291772d0
         r2=Q2 !*0.5291772d0
 


        pi = 4.0d0 * atan2(1.0d0,1.0d0)

	  imass = nint(XMASS(3))

        select case(imass)
         case default 
          write (6,"('for this mass there is no MORBID PES supplied')")
	    stop 'check the mass and pes'

	   case (32)
	  

        RE12      = 1.336024d0
        AA1       = 1.658377d0


        V0        =        .00000000d0
        RHOE      = 87.669282*pi/180.0d0
        FA1       =        .00000000d0
        FA2       =   19293.90700000d0
        FA3       =       0.00000000d0
        FA4       =    4239.02200000d0
        FA5       =    6931.27000000d0
        FA6       =    4140.37300000d0
        FA7       =        .00000000d0
        FA8       =       0.00000000d0
        F1A1      =   -3670.331d0
        F2A1      =   -7504.719d0
        F3A1      =       0.0d0
        F4A1      =       0.0d0
        F11       =   39228.206d0
        F1A11     =       0.0d0
        F2A11     =   -3008.372d0
        F3A11     =        .00000000d0
        F13       =    -375.75000000d0
        F1A13     =       0.00000000d0
        F2A13     =    4426.55700000d0
        F3A13     =       0.00000000d0
        F111      =    -898.48500000d0
        F1A111    =        .00000000d0
        F2A111    =   -4720.34500000d0
        F113      =    -301.10300000d0
        F1A113    =        .00000000d0
        F2A113    =       0.00000000d0
        F1111     =     850.09200000d0
        F1A1111   =        .00000000d0
        F1113     =        .00000000d0
        F1A1113   =        .00000000d0
        F1133     =        .00000000d0
        F1A1133   =        .00000000d0


	case(16) ! water 

 



      RE12      = 0.95784800000d0
        AA1       = 2.22600000d0


        V0        =        .00000000d0
        RHOE      = 75.4576d0*pi/180.0d0
        FA1       =         0.000000d0     
        FA2       =         18965.100000d0 
        FA3       =         1649.00000d0   
        FA4       =         4284.00000d0   
        FA5       =         0.00000d0 
        FA6       =         4412.00000d0   
        FA7       =         954.00000d0    
        FA8       =        -4822.00000d0   
        F1A1      =         -6376.00000d0  
        F2A1      =         -3093.00000d0  
        F3A1      =         -6537.00000d0  
        F4A1      =          630.00000d0   
        F11       =         42933.000001d0 
        F1A11     =         -2988.00000d0  
        F2A11     =         -5947.00000d0  
        F3A11     =          0.00000000d0  
        F13       =         -1040.000008d0 
        F1A13     =         6196.00000d0   
        F2A13     =         0.0000000d0   
        F3A13     =         0.000000d0     
        F111      =        -721.00000d0   
        F1A111    =         8236.00000d0   
        F2A111    =         3235.00000d0  
        F113      =         -1041.00000d0  
        F1A113    =          9884.00000d0  
        F2A113    =          5572.00000d0  
        F1111     =          3571.00000d0  
        F1A1111   =        -3885.00000d0   
        F1113     =         -632.00000d0   
        F1A1113   =        -7488.00000d0   
        F1133     =         0.000000d0     
        F1A1133   =       -15008.00000d0   


 
	  

        RE12      = 0.9576257000d0
        AA1       = 2.2260000d0

        V0        =        .00000000d0
        RHOE      = 75.48992000d0*pi/180.0d0
        FA1       =      0.0000d0    
        FA2       =       18902.4000d0  
        FA3       =       1961.5000d0   
        FA4       =       4134.9000d0   
        FA5       =      -1959.6000d0   
        FA6       =       4484.1000d0   
        FA7       =       3961.8000d0   
        FA8       =      -4751.6000d0   
        F1A1      =       -6132.4000d0  
        F2A1      =       -3022.9000d0  
        F3A1      =       -5951.0000d0  
        F4A1      =       1030.0000d0   
        F11       =       42927.83000d0  
        F1A11     =       -3042.6000d0  
        F2A11     =       -3669.0000d0  
        F3A11     =        0.000000d0   
        F13       =       -1046.63000d0  
        F1A13     =       6109.8000d0   
        F2A13     =       0.00000d0     
        F3A13     =       0.0000d0      
        F111      =      -1309.3000d0   
        F1A111    =       1731.000d0    
        F2A111    =      -1423.000d0    
        F113      =       -1255.9000d0  
        F1A113    =        9863.000d0   
        F2A113    =        3712.000d0   
        F1111     =        4159.0000d0  
        F1A1111   =       560.000d0     
        F1113     =       -221.1000d0   
        F1A1113   =      -7238.1000d0   
        F1133     =       0.0000d0      
        F1A1133   =       0.0000d0      




        end select 


        f1 = 0.000000
        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f1a1333 = f1a1113



! calculate potential energy function values
!
      coro=dcos(rhoe)+dcos(th)

      y1=1.00000d+00-dexp(-aa1*(r1-re12))
      y3=1.00000d+00-dexp(-aa1*(r2-re12))

!
! calculate potential energy function values
!

      vp1= fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 + 
     .     fa6*coro**6+fa7*coro**7+fa8*coro**8
      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro
      fe3333= f3333+f1a3333*coro
      fe1113= f1113+f1a1113*coro
      fe1333= f1333+f1a1333*coro
      fe1133= f1133+f1a1133*coro
      vp    =  vp1+fe1*y1+fe3*y3                             
     .        +fe11*y1**2+fe33*y3**2+fe13*y1*y3            
     .        +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3     
     .        +fe133*y1*y3**2                              
     .        +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3   
     .        +fe1333*y1*y3**3+fe1133*y1**2*y3**2

      v = (v0 + vp) ! /219474.6306700d0



        return
        end  subroutine morbid
