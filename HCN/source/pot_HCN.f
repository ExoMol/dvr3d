!    This is the VQZANO+ potential with 3 point mappings,
!      HNC morphing, relativistic (relcor2) and 
!      adiabatic (DBOCcor) corrections OF:-
!
! T. van Mourik, G. J. Harris, O. L. Polyansky, J. Tennyson,
! A. G. Csaszar and P. J. Knowles, J. Chem. Phys. 115, 3706 (2001).
!
!      Written by GJH  14th September 2000
!      Rewritten by GJH 16th March 2001

! A test main program.
c$$$      implicit double precision (a-h,o-z)
c$$$
c$$$      pot =0.0
c$$$      smalr =2.2
c$$$      bigR=3.2
c$$$      xgamma=1.0
c$$$
c$$$      call potv(pot, smallr, bigR, xgamma)
c$$$
c$$$
c$$$      end

! *************************************************

      subroutine potv(pot, smallr, bigR, xgamma)

! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! INPUT
!
! smallr is the C to N bond length.

! bigr is the H to CN center of mass distance.

! xgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0  is linear HCN
!      and an angle of Pi is CNH
!
! OUTPUT
!
! pot the potential energy in Hartree
!
! HCN minimum is at:-
! R=3.1855
! r=2.1785

      implicit double precision (a-h,o-z)
      double precision  coeff(252), bigReC(4), smallreC(3), 
     * BsmallrC(1), BbigRC(3), R_co(5), m_co(9)

      coeff(  1)= 0.1494042105d+00
      coeff(  2)= 0.2179006487D+00
      coeff(  3)=-0.3273920417D+00
      coeff(  4)= 0.8351497054D+00
      coeff(  5)=-0.6619144082D+00
      coeff(  6)= 0.4599713087D+00
      coeff(  7)=-0.2738662958D+00
      coeff(  8)=-0.4440119743D+01
      coeff(  9)= 0.1109277630D+02
      coeff( 10)=-0.1395244884D+02
      coeff( 11)= 0.1101296234D+02
      coeff( 12)=-0.7263133526D+01
      coeff( 13)= 0.2224483967D+01
      coeff( 14)=-0.5673184395D+00
      coeff( 15)=-0.1196408272D+02
      coeff( 16)= 0.3286537933D+02
      coeff( 17)=-0.3404040527D+02
      coeff( 18)= 0.2647031593D+02
      coeff( 19)=-0.1024488163D+02
      coeff( 20)= 0.4328814983D+01
      coeff( 21)= 0.6302400827D+00
      coeff( 22)=-0.1281658459D+02
      coeff( 23)= 0.3501254654D+02
      coeff( 24)=-0.2787051392D+02
      coeff( 25)= 0.2416838074D+02
      coeff( 26)=-0.3077708721D+01
      coeff( 27)= 0.2183263779D+01
      coeff( 28)= 0.2461140633D+01
      coeff( 29)=-0.5094871998D+01
      coeff( 30)= 0.1759796333D+02
      coeff( 31)=-0.7992774963D+01
      coeff( 32)= 0.1048842525D+02
      coeff( 33)= 0.2663014174D+01
      coeff( 34)= 0.9967646599D+00
      coeff( 35)= 0.2099469662D+01
      coeff( 36)=-0.7037402391D+00
      coeff( 37)= 0.3279056072D+01
      coeff( 38)=-0.3060762286D+00
      coeff( 39)= 0.1875680566D+01
      coeff( 40)= 0.9339466095D+00
      coeff( 41)= 0.2784251571D+00
      coeff( 42)= 0.6579222083D+00
      coeff( 43)= 0.4796701431D+01
      coeff( 44)=-0.1318006039D+02
      coeff( 45)= 0.1348409176D+02
      coeff( 46)=-0.1195808506D+02
      coeff( 47)= 0.6425770760D+01
      coeff( 48)=-0.2491237164D+01
      coeff( 49)= 0.6074722409D+00
      coeff( 50)=-0.7898491383D+01
      coeff( 51)= 0.2515898132D+02
      coeff( 52)=-0.2228882599D+02
      coeff( 53)= 0.2319293022D+02
      coeff( 54)=-0.7300096512D+01
      coeff( 55)= 0.4186582088D+01
      coeff( 56)= 0.1468366385D+01
      coeff( 57)=-0.5053308105D+02
      coeff( 58)= 0.1144061890D+03
      coeff( 59)=-0.1330410309D+03
      coeff( 60)= 0.8062403107D+02
      coeff( 61)=-0.4143428802D+02
      coeff( 62)= 0.5329541683D+01
      coeff( 63)=-0.1057870984D+01
      coeff( 64)=-0.6028519058D+02
      coeff( 65)= 0.1080324249D+03
      coeff( 66)=-0.1445681915D+03
      coeff( 67)= 0.5222629547D+02
      coeff( 68)=-0.3854883575D+02
      coeff( 69)=-0.6430528164D+01
      coeff( 70)=-0.2470804155D+00
      coeff( 71)=-0.2972009659D+02
      coeff( 72)= 0.2486521339D+02
      coeff( 73)=-0.6844873810D+02
      coeff( 74)=-0.2519468404D-01
      coeff( 75)=-0.2420088196D+02
      coeff( 76)=-0.9852017403D+01
      coeff( 77)= 0.1558390141D+01
      coeff( 78)=-0.5004767895D+01
      coeff( 79)=-0.3438580275D+01
      coeff( 80)=-0.1220398617D+02
      coeff( 81)=-0.2258141041D+01
      coeff( 82)=-0.8132832527D+01
      coeff( 83)=-0.3991177797D+01
      coeff( 84)= 0.1252798915D+01
      coeff( 85)= 0.2080360413D+02
      coeff( 86)=-0.4752796173D+02
      coeff( 87)= 0.6023006439D+02
      coeff( 88)=-0.4397934341D+02
      coeff( 89)= 0.2710381699D+02
      coeff( 90)=-0.9402210236D+01
      coeff( 91)= 0.3134345770D+01
      coeff( 92)=-0.1111633492D+02
      coeff( 93)= 0.1177921772D+02
      coeff( 94)=-0.3666366196D+02
      coeff( 95)= 0.1091868973D+02
      coeff( 96)=-0.2032316589D+02
      coeff( 97)=-0.1134070456D+00
      coeff( 98)=-0.3911679745D+01
      coeff( 99)=-0.1083814926D+03
      coeff(100)= 0.2490106964D+03
      coeff(101)=-0.2849703674D+03
      coeff(102)= 0.1854083252D+03
      coeff(103)=-0.1098742905D+03
      coeff(104)= 0.3372972488D+02
      coeff(105)=-0.1911353683D+02
      coeff(106)=-0.9813106537D+02
      coeff(107)= 0.2571860962D+03
      coeff(108)=-0.2505248108D+03
      coeff(109)= 0.1899570007D+03
      coeff(110)=-0.1025988770D+03
      coeff(111)= 0.2444090652D+02
      coeff(112)=-0.2244192505D+02
      coeff(113)=-0.4234691143D+01
      coeff(114)= 0.8552397156D+02
      coeff(115)=-0.4999621201D+02
      coeff(116)= 0.1141592178D+03
      coeff(117)=-0.4533333969D+02
      coeff(118)=-0.1452669501D+01
      coeff(119)=-0.1024097919D+02
      coeff(120)= 0.1400469303D+02
      coeff(121)= 0.2283895761D+00
      coeff(122)= 0.5921393037D+00
      coeff(123)= 0.3907504272D+02
      coeff(124)=-0.9498566628D+01
      coeff(125)=-0.6402152061D+01
      coeff(126)=-0.4512357712D+00
      coeff(127)= 0.1556770229D+02
      coeff(128)=-0.4285261154D+02
      coeff(129)= 0.3758481216D+02
      coeff(130)=-0.2836553383D+02
      coeff(131)= 0.5195811749D+01
      coeff(132)=-0.2380146027D+01
      coeff(133)=-0.2525846958D+01
      coeff(134)=-0.6172574234D+02
      coeff(135)= 0.1954359741D+03
      coeff(136)=-0.1905524139D+03
      coeff(137)= 0.1957860260D+03
      coeff(138)=-0.9298155975D+02
      coeff(139)= 0.5106489563D+02
      coeff(140)=-0.1235393524D+02
      coeff(141)=-0.2260525818D+03
      coeff(142)= 0.6543347168D+03
      coeff(143)=-0.5629700317D+03
      coeff(144)= 0.5343278198D+03
      coeff(145)=-0.2099469452D+03
      coeff(146)= 0.1321243896D+03
      coeff(147)=-0.1978376007D+02
      coeff(148)=-0.1678566437D+03
      coeff(149)= 0.5674368286D+03
      coeff(150)=-0.4264934082D+03
      coeff(151)= 0.5247149048D+03
      coeff(152)=-0.1254009094D+03
      coeff(153)= 0.1368169861D+03
      coeff(154)=-0.1946393394D+02
      coeff(155)= 0.6929652691D+01
      coeff(156)= 0.1203537521D+03
      coeff(157)=-0.1582188873D+03
      coeff(158)= 0.2871345215D+03
      coeff(159)=-0.1572170830D+02
      coeff(160)= 0.5963705826D+02
      coeff(161)=-0.9650691986D+01
      coeff(162)= 0.3750868607D+02
      coeff(163)=-0.1728201675D+02
      coeff(164)=-0.4348765564D+02
      coeff(165)= 0.8167091370D+02
      coeff(166)= 0.7986162663D+01
      coeff(167)= 0.3591517687D+01
      coeff(168)=-0.1310849667D+01
      coeff(169)=-0.8493385315D+01
      coeff(170)= 0.2844993019D+02
      coeff(171)=-0.3800553131D+02
      coeff(172)= 0.4474219894D+02
      coeff(173)=-0.2850215530D+02
      coeff(174)= 0.1646527290D+02
      coeff(175)=-0.3970207691D+01
      coeff(176)=-0.1542687836D+03
      coeff(177)= 0.3372832642D+03
      coeff(178)=-0.4464781494D+03
      coeff(179)= 0.2951505432D+03
      coeff(180)=-0.1816697998D+03
      coeff(181)= 0.5316984558D+02
      coeff(182)=-0.6686333179D+01
      coeff(183)=-0.3925593567D+03
      coeff(184)= 0.6971489258D+03
      coeff(185)=-0.9502923584D+03
      coeff(186)= 0.4538467407D+03
      coeff(187)=-0.3279503174D+03
      coeff(188)= 0.6674200439D+02
      coeff(189)=-0.5604239941D+01
      coeff(190)=-0.3418550720D+03
      coeff(191)= 0.4356850586D+03
      coeff(192)=-0.8969200439D+03
      coeff(193)= 0.2099159241D+03
      coeff(194)=-0.2829770813D+03
      coeff(195)= 0.5536326981D+02
      coeff(196)=-0.5499388218D+01
      coeff(197)=-0.6700291443D+02
      coeff(198)= 0.1543539429D+02
      coeff(199)=-0.4875061340D+03
      coeff(200)= 0.2664015007D+02
      coeff(201)=-0.1126098557D+03
      coeff(202)= 0.2384203148D+02
      coeff(203)=-0.3890762568D+01
      coeff(204)= 0.2440213394D+02
      coeff(205)=-0.5005599213D+02
      coeff(206)=-0.1430028381D+03
      coeff(207)= 0.8150361061D+01
      coeff(208)=-0.6253750801D+01
      coeff(209)= 0.1309172034D+01
      coeff(210)=-0.1360354185D+01
      coeff(211)=-0.1142237091D+02
      coeff(212)= 0.3025738716D+02
      coeff(213)=-0.3869804764D+02
      coeff(214)= 0.3254767227D+02
      coeff(215)=-0.1961404419D+02
      coeff(216)= 0.7705432892D+01
      coeff(217)=-0.2599303424D+00
      coeff(218)=-0.7386426544D+02
      coeff(219)= 0.2132664490D+03
      coeff(220)=-0.1990179596D+03
      coeff(221)= 0.1671524963D+03
      coeff(222)=-0.6041564560D+02
      coeff(223)= 0.1448185158D+02
      coeff(224)= 0.1541442037D+01
      coeff(225)=-0.1405283356D+03
      coeff(226)= 0.4659133911D+03
      coeff(227)=-0.2857222900D+03
      coeff(228)= 0.2609895935D+03
      coeff(229)=-0.6865476227D+02
      coeff(230)= 0.1467753887D+02
      coeff(231)= 0.5323611736D+01
      coeff(232)=-0.6710372162D+02
      coeff(233)= 0.4675559082D+03
      coeff(234)=-0.1301919403D+03
      coeff(235)= 0.1887061920D+03
      coeff(236)=-0.4263847733D+02
      coeff(237)= 0.1578407860D+02
      coeff(238)= 0.6043389320D+01
      coeff(239)= 0.4931479263D+02
      coeff(240)= 0.2340046234D+03
      coeff(241)=-0.2837063980D+02
      coeff(242)= 0.7166375732D+02
      coeff(243)=-0.9170268059D+01
      coeff(244)= 0.9936111450D+01
      coeff(245)= 0.2666974306D+01
      coeff(246)= 0.4050671005D+02
      coeff(247)= 0.5192458344D+02
      coeff(248)=-0.1643330765D+02
      coeff(249)= 0.1117124939D+02
      coeff(250)= 0.3955295801D+01
      coeff(251)= 0.1787165642D+01
      coeff(252)= 0.7089445740D-01

      bigReC(1)=3.30884099158515088d0
      bigReC(2)=-2.45634163220171997d0
      bigReC(3)=-1.0989512458800641d0
      bigReC(4)=0.88049730604354981d0

      BbigRC(1)=0.833894529530727957d0
      BbigRC(2)=-0.0635504047336145489d0
      BbigRC(3)=0.284315972604355038d0

      smallreC(1)=2.253505169974654d0
      smallreC(2)=-0.891025354333305986d0
      smallreC(3)=2.22802201493595842d0

      BsmallrC(1)=1.66063327384055204d0
      pot = 0.0d0
      icoefref = 1

! Least Energy Isomerisation PATH COEFFICENTS

      R_co(1) = 2.06021279535746d0
      R_co(2) = 0.206149293036705d0
      R_co(3) = 1.38850960996447d0
      R_co(4) = -0.0678082821258994d0
      R_co(5) = -0.412935160444408d0

! HNC R MOPRHING COEFFIECNTS

      m_co(1) = -0.000303673553151895d0
      m_co(2) = 0.00435465925587477d0
      m_co(3) = 0.00402462914869603d0
      m_co(4) = 0.0139187573949336d0
      m_co(5) = 0.0416697727891190d0
      m_co(6) = 0.0264801392783952d0
      m_co(7) = -0.0216034082004957d0
      m_co(8) = -0.0493168970474053d0
      m_co(9) = -0.0297172762227854d0

      gamma = acos(xgamma)

! CRITCIAL POINT COORDINATE REMAPPING.

      sr = smallr + 0.001635d0 + 
     *     (0.0018646949949d0 * ((sin(gamma)+
     *     (0.131325005d0*sin(gamma*2.0d0)))**4))

      br = bigR + 0.002205d0 -
     *    (0.014287178436d0 * ((sin(gamma)+
     *    (0.131325005d0*sin(gamma*2.0d0)))**4))

      gamma = gamma+ (0.0105990601023d0* 
     *    ((sin(gamma)+(0.131325005d0*sin(gamma*2.0d0)))**16))

      cosgamma = cos(gamma)

      ! BEGIN HNC R MORPHING

      Rmin = 0.0d0
      do j=0, 4
        Rmin = Rmin + (R_co(j+1)*(cosgamma**j))
      enddo

      tempR = bR-Rmin
      shift1=0.0d0

      do i=0, 2
         do j=0, 2

            shift1 = shift1 + (m_co(j+(i*3)+1) * (tempR**i)
     *               * (cosgamma**j)) 
         enddo
      enddo

      bR=bR+(shift1*
     * exp(-((0.6d0*(acos(cosgamma)-3.1415927d0))**6)))

! NEXT THE MAIN FIT CALCULATION

      rLEbigR =  bigReC(1)+(bigReC(2)*cosgamma)
     * +(bigReC(3)*cosgamma*cosgamma)
     * +(bigReC(4)*cosgamma*cosgamma*sr)

      rLEsmallr = smallreC(1)+(smallreC(2)*cosgamma)
     * +(smallreC(3)*cosgamma*cosgamma)

      betabigR= BbigRC(1) + (BbigRC(2) * cosgamma)
     * + (BbigRC(3) * cosgamma * cosgamma)

      betasmallr = BsmallrC(1)

      do i=0, 5
               
         call rmorsecoord(bR, betabigR, rLEbigR, i, bigmorse)

         do j=0, 5

      call rmorsecoord(sr, betasmallr, rLEsmallr, j, smallmorse)

            do k=0, 6

               call plgndr(k, 0, cosgamma, apolynom)

            temp = coeff(icoefref)*bigmorse*smallmorse*apolynom

!       write(6,1001) "      coeff(", icoefref,")=", coeff(icoefref)
! 1001       format(a12,i3,a2,d17.10)
               pot = pot + temp

               icoefref=icoefref+1
            enddo
         enddo
      enddo

! CRITICAL POINT ENERGY SCALING NEXT

      pot = pot +((-1.2381427d-7 *
     * (atan(15.0d0*(gamma-1.2d0))+1.5708d0)) +
     * (1.2115228d-04 *
     * ((sin(gamma)+(0.131325005d0*sin(gamma*2d0)))**16)))

!     lower limits to potential (hole patches)


! general low R patch
      if (bigR .LE. 1.0) then
         pot=100.0d0
      endif

!     Patch for irritating holes at higher low R
      if (cosgamma .le. 0.969d0 .and. cosgamma .gt. 0.071d0) then
!       calculate approxomate LEI
      aprxRbig = 2.113058664103d0 
     * + (cosgamma * 0.1544661511693d0)
     * + (cosgamma * cosgamma * 0.9706585547594d0)
      if (bigR .LE. (aprxRbig-0.8d0)) then
         pot=100.0d0
      endif
      endif

! small r patch
      if (smallr .LE. 1.7d0) then 
       pot = 100.0d0
      endif

!     patch for large smallr
      if (smallr .GE. 3.5d0) then 
        pot = 100.0d0
      endif

!     patch for large R
! general upper limit
      if(bigR .GE. 7.0d0) then
        pot = 100.0d0
      endif

! angular conditions upper limit patchs
      if (bigR .GE. 6.5d0 .AND. cosgamma .le. 0.35d0) then
        pot = 100.0d0
      endif
      if (bigR .GE. 6.0d0 .AND. cosgamma .le. 0.3d0) then
        pot = 100.0d0
      endif
      if (bigR .GE. 5.5d0 .AND. cosgamma .le. 0.23173d0) then
        pot = 100.0d0
      endif

      if (bigR .GE. 5.0d0 .AND. cosgamma .le. 0.12945d0) then
        pot = 100.0
      endif

      if (bigR .GE. 4.8d0 .AND. cosgamma .le. -0.12944d0) then
        pot = 100.0d0
      endif

      pot = pot-4.250289d-5  ! zero potential at HCN minimum

!	call relativistic  correction

      call rel_corr(rcor, smallr, bigR, cosgamma)
      pot = pot + rcor

! 	call DBOC correction

      call DBOC_corr(DBOCcor, smallr, bigR, cosgamma)
       pot = pot + DBOCcor

      return
      end


! ************************************************************************

      subroutine plgndr(l, m, x, pmm)

! subroutine to calculate the asociated legenre polynomials, translated drectly from a C 
! version to fortran version. The C version comes from Numerical recipies in C,  Press et al. 
      implicit double precision (a-h, o-z) 
      double precision fact, pll, pmm, pmmp1, somx2, x, temp
      integer i, ll

      if (m .LT. 0 .OR. m .GT. l .OR. (x*x) .GT. 1.0) then
         write(6,*)  "!ERROR! bad arguments in routine plgndr"
         write(6,*)  m,l,x*x
      endif

      pmm = 1.0d0

      if (m .GT. 0) then
   
         somx2=sqrt((1.0d0 - x) * (1.0d0 + x))
         fact = 1.0d0
         do i=1, m
            pmm = pmm*(-fact * somx2)
            fact = fact+ 2.0d0
         enddo
      endif

      if (l .NE. m) then
         pmmp1 = x * (2.0d0 * dble(m) + 1.0d0) * pmm

         if (l .EQ. (m + 1)) then
            pmm=pmmp1
         else
            do ll = m+2, l
           fm = dble(m)
           fll = ll
           
          pll=((x*((2.0d0*fll)-1.0d0)*pmmp1)-((fll+fm-1.0d0)*pmm))
          temp = (fll-fm)
          pll = pll / temp
            pmm=pmmp1
            pmmp1=pll
            enddo
            pmm=pll
         endif
      endif

!      write(6,*) l, x, pmm

      return
      end

! ********************************************************************

      subroutine rmorsecoord(R, beta, requalib, n, rmor)

      implicit double precision (a-h,o-z)

      ratio = -beta * ((R - requalib) /requalib)
      temp = (1.0d0 - exp(ratio))

      rmor = temp**n

!      write (6,*) temp, rmor, n

      return
      end

! *************************************************

      subroutine DBOC_corr(pot, smallr, bigR, cosgamma)

! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! adiabatic correction surface.
!
! INPUT
! smallr is the C to N bond length.

! bigr is the H to CN center of mass distance.

! cosgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0 coresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH
!
! OUTPUT
!
! pot the adiabatic correction in Hartree
!

      implicit double precision (a-h,o-z)
      double precision  coeff(60)
      EH = 219474.624d0

      coeff(1)=61300.7135899d0
      coeff(2)=895.886162423d0
      coeff(3)=-70565.8990549d0
      coeff(4)=4773.9817209d0
      coeff(5)=8856.65590136d0
      coeff(6)=-80797.3926551d0
      coeff(7)=-1682.1616785d0
      coeff(8)=92565.9091745d0
      coeff(9)=-5602.42892746d0
      coeff(10)=-9198.09369058d0
      coeff(11)=35939.9687435d0
      coeff(12)=929.150781931d0
      coeff(13)=-40169.8574301d0
      coeff(14)=2130.09338513d0
      coeff(15)=2704.67838662d0
      coeff(16)=-5313.95221068d0
      coeff(17)=-157.566257215d0
      coeff(18)=5780.36035164d0
      coeff(19)=-260.350322598d0
      coeff(20)=-167.44015843d0
      coeff(21)=-55370.316994d0
      coeff(22)=4754.39164204d0
      coeff(23)=83030.3122802d0
      coeff(24)=-8281.26766355d0
      coeff(25)=-26969.0665816d0
      coeff(26)=74094.2408883d0
      coeff(27)=-5925.76456913d0
      coeff(28)=-110370.495912d0
      coeff(29)=10527.0804951d0
      coeff(30)=34834.6549419d0
      coeff(31)=-33001.6947725d0
      coeff(32)=2471.70897939d0
      coeff(33)=48649.3467107d0
      coeff(34)=-4437.80010694d0
      coeff(35)=-14790.4934647d0
      coeff(36)=4887.37761221d0
      coeff(37)=-345.374981145d0
      coeff(38)=-7117.25664686d0
      coeff(39)=620.415816354d0
      coeff(40)=2062.69994319d0
      coeff(41)=12554.0937029d0
      coeff(42)=-1936.70511436d0
      coeff(43)=-21736.5239009d0
      coeff(44)=2472.52194725d0
      coeff(45)=9092.30843229d0
      coeff(46)=-16808.5214793d0
      coeff(47)=2498.58479994d0
      coeff(48)=29077.903879d0
      coeff(49)=-3205.61349595d0
      coeff(50)=-12071.3921973d0
      coeff(51)=7489.67926539d0
      coeff(52)=-1077.32941536d0
      coeff(53)=-12914.6546977d0
      coeff(54)=1383.34822742d0
      coeff(55)=5304.19453882d0
      coeff(56)=-1109.80874249d0
      coeff(57)=155.135287705d0
      coeff(58)=1905.01555022d0
      coeff(59)=-198.680122409d0
      coeff(60)=-771.435653363d0


      pot = -840.6019631d0 ! averge value of correction
      icoef = 1

      do i=0, 2
               
         bigmorse=bigR**i

         do j=0, 3

         smallmorse=smallr**j

            do k=0, 4

               apolynom=cosgamma**k

         temp = coeff(icoef)*bigmorse*smallmorse*apolynom

               pot = pot + temp

               icoef=icoef+1
            enddo
         enddo
      enddo

      pot = pot / EH ! CONVERT FORM CM-1 TO HARTREE

      return
      end


! ************************************************************************
 
      subroutine rel_corr(pot, smallr, bigR, cosgamma)
!
! The relativistic correction surface.
!
! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! INPUT
! smallr is the C to N bond length.
!
! bigr is the H to CN center of mass distance.
!
! cosgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0 corresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH
!
! OUTPUT
!
! pot the relativistic correction in Hartree
!

      implicit double precision (a-h,o-z)
      double precision  coeff(80)

      coeff(1)=-0.0805585205815292d0
      coeff(2)=0.0950673303720766d0
      coeff(3)=-0.253364906126069d0
      coeff(4)=-0.084129451710227d0
      coeff(5)=-0.0745419772100667d0
      coeff(6)=0.0477227356064123d0
      coeff(7)=-0.130026022694324d0
      coeff(8)=0.335629886457965d0
      coeff(9)=0.107905706720457d0
      coeff(10)=0.0993733461559782d0
      coeff(11)=-0.0237612599129882d0
      coeff(12)=0.0595297964464055d0
      coeff(13)=-0.149509867823676d0
      coeff(14)=-0.0460734392042994d0
      coeff(15)=-0.0439643487346239d0
      coeff(16)=0.00388436986158473d0
      coeff(17)=-0.00906898001283882d0
      coeff(18)=0.0221602015192701d0
      coeff(19)=0.0065691932071385d0
      coeff(20)=0.00645621168784551d0
      coeff(21)=0.0266090487964399d0
      coeff(22)=-0.115552391578968d0
      coeff(23)=0.252922948126617d0
      coeff(24)=0.0870075200551087d0
      coeff(25)=0.0506809331899964d0
      coeff(26)=-0.041492313771614d0
      coeff(27)=0.156584563160005d0
      coeff(28)=-0.332335667002259d0
      coeff(29)=-0.111669707104665d0
      coeff(30)=-0.06770235689405d0
      coeff(31)=0.0220056679462455d0
      coeff(32)=-0.0710958741686161d0
      coeff(33)=0.146735654274853d0
      coeff(34)=0.047768960584727d0
      coeff(35)=0.0298863818957252d0
      coeff(36)=-0.00375993578145236d0
      coeff(37)=0.0107485671294735d0
      coeff(38)=-0.021583118481192d0
      coeff(39)=-0.00682587066702194d0
      coeff(40)=-0.00437935957298068d0
      coeff(41)=-0.0131085901581377d0
      coeff(42)=0.0460835147852451d0
      coeff(43)=-0.0826635630451006d0
      coeff(44)=-0.0312833158935659d0
      coeff(45)=-0.00904915849992159d0
      coeff(46)=0.0193658221974871d0
      coeff(47)=-0.0618956123071639d0
      coeff(48)=0.107774673538991d0
      coeff(49)=0.040184308696244d0
      coeff(50)=0.0123015451000135d0
      coeff(51)=-0.00973451656341789d0
      coeff(52)=0.0278651342014908d0
      coeff(53)=-0.0471541574047327d0
      coeff(54)=-0.0172250650594388d0
      coeff(55)=-0.00546984649469844d0
      coeff(56)=0.001602231058765d0
      coeff(57)=-0.00418032677019256d0
      coeff(58)=0.00687897584843289d0
      coeff(59)=0.00246635113824095d0
      coeff(60)=0.000805450828032454d0
      coeff(61)=0.00229338173792963d0
      coeff(62)=-0.00610190383854349d0
      coeff(63)=0.00863846660156985d0
      coeff(64)=0.00393643841031911d0
      coeff(65)=0.000202814342030121d0
      coeff(66)=-0.00323246774517316d0
      coeff(67)=0.00812760212687042d0
      coeff(68)=-0.011178025856714d0
      coeff(69)=-0.00505844398525825d0
      coeff(70)=-0.000330146921869386d0
      coeff(71)=0.00154631457145571d0
      coeff(72)=-0.0036279978431523d0
      coeff(73)=0.00484658611284714d0
      coeff(74)=0.00217108938922934d0
      coeff(75)=0.000160901000351014d0
      coeff(76)=-0.000244921320434022d0
      coeff(77)=0.000539999643014766d0
      coeff(78)=-0.000700947801549237d0
      coeff(79)=-0.000311210335395826d0
      coeff(80)=-2.53086031415751d-5

      pot = 0.0472411971073d0
      icoef = 1

      do i=0, 3
               
         bigmorse=bigR**i

         do j=0, 3

         smallmorse=smallr**j

            do k=0, 4

               apolynom=cosgamma**k

            temp = coeff(icoef)*bigmorse*smallmorse*apolynom

               pot = pot + temp

               icoef=icoef+1
            enddo
         enddo
      enddo

      return
      end

