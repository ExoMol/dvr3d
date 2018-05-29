subroutine DIPS(muX,xpx,npropinx, muY, xpy,npropiny, r12, r32, alpha)
implicit none
!==============================================================================================================
!
! INPUT:  r12/ang, r32/ang, alpha /radians
! OUTPUT: mu_Y(a.u.), mu_X(a.u.) 
! mu_Y is the "small" component, perpendicular to the bond-angle bisector  (=0 a.u. at equilibrium)
! mu_X is the "large" component, parallel      to the bond-angle bisector  (=0.730 a.u. at equilibrium)

! In bond embedding both components are non-zero at equilibrium
!==============================================================================================================

 double precision, intent(out) :: muX, muY
 integer  npropinx, npropiny
 double precision xpx(300),xpy(300),xp(300)

      double precision            :: r12,r32,alpha,aa1,pi,re12,alphae,xst,y1,y2,xs1,xs2,v0,vp1,vp2,vp3,re32,aa2
      double precision            :: g1,g2,b1,b2,rhh,vhh,v
      integer(4) :: i 
      !
      pi = 4.0d0 * atan2(1.0d0,1.0d0)



!calculation of muX

      xp=xpx

      re12      = xp(2)
      re32	= xp(248)
      alphae    = xp(3)*pi/180.0d0
      !
      aa1  = xp(4)
      aa2= xp(249)
      b1   = xp(5)
      b2   = xp(6)
      g1   = xp(7)
      g2   = xp(8)
      !

      !y1=1.0d+00-exp(-aa1*(r12-re12))
      !y2=1.0d+00-exp(-aa2*(r32-re32))
      y1=1.0d+00-exp(-aa1*(r12-re12))
      y2=1.0d+00-exp(-aa2*(r32-re32))
      xst=alpha
      !xst=alphae-alpha
	!y1=r12-re12
	!y2=r32-re32  
	!
      xs1 = (y1+y2)*0.5d0
      xs2 = (y1-y2)*0.5d0

  v0= xp(  1)*xs1**0*xs2**0*xst**0

 vp1= xp(  9)*xs1**0*xs2**0*xst**1&
     +xp( 10)*xs1**1*xs2**0*xst**0&
     +xp( 11)*xs1**0*xs2**0*xst**2&
     +xp( 12)*xs1**0*xs2**2*xst**0&
     +xp( 13)*xs1**1*xs2**0*xst**1&
     +xp( 14)*xs1**2*xs2**0*xst**0&
     +xp( 15)*xs1**0*xs2**0*xst**3&
     +xp( 16)*xs1**0*xs2**2*xst**1&
     +xp( 17)*xs1**1*xs2**0*xst**2&
     +xp( 18)*xs1**1*xs2**2*xst**0&
     +xp( 19)*xs1**2*xs2**0*xst**1&
     +xp( 20)*xs1**3*xs2**0*xst**0&
     +xp( 21)*xs1**0*xs2**0*xst**4&
     +xp( 22)*xs1**0*xs2**2*xst**2&
     +xp( 23)*xs1**0*xs2**4*xst**0&
     +xp( 24)*xs1**1*xs2**0*xst**3&
     +xp( 25)*xs1**1*xs2**2*xst**1&
     +xp( 26)*xs1**2*xs2**0*xst**2&
     +xp( 27)*xs1**2*xs2**2*xst**0&
     +xp( 28)*xs1**3*xs2**0*xst**1&
     +xp( 29)*xs1**4*xs2**0*xst**0&
     +xp( 30)*xs1**0*xs2**0*xst**5&
     +xp( 31)*xs1**0*xs2**2*xst**3&
     +xp( 32)*xs1**0*xs2**4*xst**1&
     +xp( 33)*xs1**1*xs2**0*xst**4&
     +xp( 34)*xs1**1*xs2**2*xst**2&
     +xp( 35)*xs1**1*xs2**4*xst**0&
     +xp( 36)*xs1**2*xs2**0*xst**3&
     +xp( 37)*xs1**2*xs2**2*xst**1&
     +xp( 38)*xs1**3*xs2**0*xst**2&
     +xp( 39)*xs1**3*xs2**2*xst**0&
     +xp( 40)*xs1**4*xs2**0*xst**1&
     +xp( 41)*xs1**5*xs2**0*xst**0&
     +xp( 42)*xs1**0*xs2**0*xst**6&
     +xp( 43)*xs1**0*xs2**2*xst**4&
     +xp( 44)*xs1**0*xs2**4*xst**2&
     +xp( 45)*xs1**0*xs2**6*xst**0&
     +xp( 46)*xs1**1*xs2**0*xst**5&
     +xp( 47)*xs1**1*xs2**2*xst**3&
     +xp( 48)*xs1**1*xs2**4*xst**1&
     +xp( 49)*xs1**2*xs2**0*xst**4&
     +xp( 50)*xs1**2*xs2**2*xst**2&
     +xp( 51)*xs1**2*xs2**4*xst**0&
     +xp( 52)*xs1**3*xs2**0*xst**3&
     +xp( 53)*xs1**3*xs2**2*xst**1&
     +xp( 54)*xs1**4*xs2**0*xst**2&
     +xp( 55)*xs1**4*xs2**2*xst**0&
     +xp( 56)*xs1**5*xs2**0*xst**1&
     +xp( 57)*xs1**6*xs2**0*xst**0&
     +xp( 58)*xs1**0*xs2**0*xst**7&
     +xp( 59)*xs1**0*xs2**2*xst**5&
     +xp( 60)*xs1**0*xs2**4*xst**3&
     +xp( 61)*xs1**0*xs2**6*xst**1&
     +xp( 62)*xs1**1*xs2**0*xst**6&
     +xp( 63)*xs1**1*xs2**2*xst**4&
     +xp( 64)*xs1**1*xs2**4*xst**2&
     +xp( 65)*xs1**1*xs2**6*xst**0&
     +xp( 66)*xs1**2*xs2**0*xst**5&
     +xp( 67)*xs1**2*xs2**2*xst**3&
     +xp( 68)*xs1**2*xs2**4*xst**1&
     +xp( 69)*xs1**3*xs2**0*xst**4&
     +xp( 70)*xs1**3*xs2**2*xst**2&
     +xp( 71)*xs1**3*xs2**4*xst**0&
     +xp( 72)*xs1**4*xs2**0*xst**3&
     +xp( 73)*xs1**4*xs2**2*xst**1&
     +xp( 74)*xs1**5*xs2**0*xst**2&
     +xp( 75)*xs1**5*xs2**2*xst**0&
     +xp( 76)*xs1**6*xs2**0*xst**1&
     +xp( 77)*xs1**7*xs2**0*xst**0&
     +xp( 78)*xs1**0*xs2**0*xst**8&
     +xp( 79)*xs1**0*xs2**2*xst**6&
     +xp( 80)*xs1**0*xs2**4*xst**4&
     +xp( 81)*xs1**0*xs2**6*xst**2&
     +xp( 82)*xs1**0*xs2**8*xst**0&
     +xp( 83)*xs1**1*xs2**0*xst**7&
     +xp( 84)*xs1**1*xs2**2*xst**5&
     +xp( 85)*xs1**1*xs2**4*xst**3&
     +xp( 86)*xs1**1*xs2**6*xst**1&
     +xp( 87)*xs1**2*xs2**0*xst**6&
     +xp( 88)*xs1**2*xs2**2*xst**4&
     +xp( 89)*xs1**2*xs2**4*xst**2&
     +xp( 90)*xs1**2*xs2**6*xst**0&
     +xp( 91)*xs1**3*xs2**0*xst**5&
     +xp( 92)*xs1**3*xs2**2*xst**3&
     +xp( 93)*xs1**3*xs2**4*xst**1&
     +xp( 94)*xs1**4*xs2**0*xst**4&
     +xp( 95)*xs1**4*xs2**2*xst**2&
     +xp( 96)*xs1**4*xs2**4*xst**0&
     +xp( 97)*xs1**5*xs2**0*xst**3&
     +xp( 98)*xs1**5*xs2**2*xst**1&
     +xp( 99)*xs1**6*xs2**0*xst**2

 vp2= xp(100)*xs1**6*xs2**2*xst**0&
     +xp(101)*xs1**7*xs2**0*xst**1&
     +xp(102)*xs1**8*xs2**0*xst**0&
     +xp(103)*xs1**0*xs2**0*xst**9&
     +xp(104)*xs1**0*xs2**2*xst**7&
     +xp(105)*xs1**0*xs2**4*xst**5&
     +xp(106)*xs1**0*xs2**6*xst**3&
     +xp(107)*xs1**0*xs2**8*xst**1&
     +xp(108)*xs1**1*xs2**0*xst**8&
     +xp(109)*xs1**1*xs2**2*xst**6&
     +xp(110)*xs1**1*xs2**4*xst**4&
     +xp(111)*xs1**1*xs2**6*xst**2&
     +xp(112)*xs1**1*xs2**8*xst**0&
     +xp(113)*xs1**2*xs2**0*xst**7&
     +xp(114)*xs1**2*xs2**2*xst**5&
     +xp(115)*xs1**2*xs2**4*xst**3&
     +xp(116)*xs1**2*xs2**6*xst**1&
     +xp(117)*xs1**3*xs2**0*xst**6&
     +xp(118)*xs1**3*xs2**2*xst**4&
     +xp(119)*xs1**3*xs2**4*xst**2&
     +xp(120)*xs1**3*xs2**6*xst**0&
     +xp(121)*xs1**4*xs2**0*xst**5&
     +xp(122)*xs1**4*xs2**2*xst**3&
     +xp(123)*xs1**4*xs2**4*xst**1&
     +xp(124)*xs1**5*xs2**0*xst**4&
     +xp(125)*xs1**5*xs2**2*xst**2&
     +xp(126)*xs1**5*xs2**4*xst**0&
     +xp(127)*xs1**6*xs2**0*xst**3&
     +xp(128)*xs1**6*xs2**2*xst**1&
     +xp(129)*xs1**7*xs2**0*xst**2&
     +xp(130)*xs1**7*xs2**2*xst**0&
     +xp(131)*xs1**8*xs2**0*xst**1&
     +xp(132)*xs1**9*xs2**0*xst**0&
     +xp(133)*xs1**0*xs2**0*xst**10&
     +xp(134)*xs1**0*xs2**2*xst**8&
     +xp(135)*xs1**0*xs2**4*xst**6&
     +xp(136)*xs1**0*xs2**6*xst**4&
     +xp(137)*xs1**0*xs2**8*xst**2&
     +xp(138)*xs1**0*xs2**10*xst**0&
     +xp(139)*xs1**1*xs2**0*xst**9&
     +xp(140)*xs1**1*xs2**2*xst**7&
     +xp(141)*xs1**1*xs2**4*xst**5&
     +xp(142)*xs1**1*xs2**6*xst**3&
     +xp(143)*xs1**1*xs2**8*xst**1&
     +xp(144)*xs1**2*xs2**0*xst**8&
     +xp(145)*xs1**2*xs2**2*xst**6&
     +xp(146)*xs1**2*xs2**4*xst**4&
     +xp(147)*xs1**2*xs2**6*xst**2&
     +xp(148)*xs1**2*xs2**8*xst**0&
     +xp(149)*xs1**3*xs2**0*xst**7&
     +xp(150)*xs1**3*xs2**2*xst**5&
     +xp(151)*xs1**3*xs2**4*xst**3&
     +xp(152)*xs1**3*xs2**6*xst**1&
     +xp(153)*xs1**4*xs2**0*xst**6&
     +xp(154)*xs1**4*xs2**2*xst**4&
     +xp(155)*xs1**4*xs2**4*xst**2&
     +xp(156)*xs1**4*xs2**6*xst**0&
     +xp(157)*xs1**5*xs2**0*xst**5&
     +xp(158)*xs1**5*xs2**2*xst**3&
     +xp(159)*xs1**5*xs2**4*xst**1&
     +xp(160)*xs1**6*xs2**0*xst**4&
     +xp(161)*xs1**6*xs2**2*xst**2&
     +xp(162)*xs1**6*xs2**4*xst**0&
     +xp(163)*xs1**7*xs2**0*xst**3&
     +xp(164)*xs1**7*xs2**2*xst**1&
     +xp(165)*xs1**8*xs2**0*xst**2&
     +xp(166)*xs1**8*xs2**2*xst**0&
     +xp(167)*xs1**9*xs2**0*xst**1&
     +xp(168)*xs1**10*xs2**0*xst**0&
     +xp(169)*xs1**0*xs2**0*xst**11&
     +xp(170)*xs1**0*xs2**2*xst**9&
     +xp(171)*xs1**0*xs2**4*xst**7&
     +xp(172)*xs1**0*xs2**6*xst**5&
     +xp(173)*xs1**0*xs2**8*xst**3&
     +xp(174)*xs1**0*xs2**10*xst**1&
     +xp(175)*xs1**1*xs2**0*xst**10&
     +xp(176)*xs1**1*xs2**2*xst**8&
     +xp(177)*xs1**1*xs2**4*xst**6&
     +xp(178)*xs1**1*xs2**6*xst**4&
     +xp(179)*xs1**1*xs2**8*xst**2&
     +xp(180)*xs1**1*xs2**10*xst**0&
     +xp(181)*xs1**2*xs2**0*xst**9

 vp3= xp(182)*xs1**2*xs2**2*xst**7&
     +xp(183)*xs1**2*xs2**4*xst**5&
     +xp(184)*xs1**2*xs2**6*xst**3&
     +xp(185)*xs1**2*xs2**8*xst**1&
     +xp(186)*xs1**3*xs2**0*xst**8&
     +xp(187)*xs1**3*xs2**2*xst**6&
     +xp(188)*xs1**3*xs2**4*xst**4&
     +xp(189)*xs1**3*xs2**6*xst**2&
     +xp(190)*xs1**3*xs2**8*xst**0&
     +xp(191)*xs1**4*xs2**0*xst**7&
     +xp(192)*xs1**4*xs2**2*xst**5&
     +xp(193)*xs1**4*xs2**4*xst**3&
     +xp(194)*xs1**4*xs2**6*xst**1&
     +xp(195)*xs1**5*xs2**0*xst**6&
     +xp(196)*xs1**5*xs2**2*xst**4&
     +xp(197)*xs1**5*xs2**4*xst**2&
     +xp(198)*xs1**5*xs2**6*xst**0&
     +xp(199)*xs1**6*xs2**0*xst**5&
     +xp(200)*xs1**6*xs2**2*xst**3&
     +xp(201)*xs1**6*xs2**4*xst**1&
     +xp(202)*xs1**7*xs2**0*xst**4&
     +xp(203)*xs1**7*xs2**2*xst**2&
     +xp(204)*xs1**7*xs2**4*xst**0&
     +xp(205)*xs1**8*xs2**0*xst**3&
     +xp(206)*xs1**8*xs2**2*xst**1&
     +xp(207)*xs1**9*xs2**0*xst**2&
     +xp(208)*xs1**9*xs2**2*xst**0&
     +xp(209)*xs1**0*xs2**0*xst**12&
     +xp(210)*xs1**0*xs2**2*xst**10&
     +xp(211)*xs1**0*xs2**4*xst**8&
     +xp(212)*xs1**0*xs2**6*xst**6&
     +xp(213)*xs1**0*xs2**8*xst**4&
     +xp(214)*xs1**0*xs2**10*xst**2&
     +xp(215)*xs1**0*xs2**12*xst**0&
     +xp(216)*xs1**1*xs2**0*xst**11&    
     +xp(217)*xs1**1*xs2**2*xst**9&
     +xp(218)*xs1**1*xs2**4*xst**7&
     +xp(219)*xs1**1*xs2**6*xst**5&
     +xp(220)*xs1**1*xs2**8*xst**3&
     +xp(221)*xs1**1*xs2**10*xst**1&
     +xp(222)*xs1**2*xs2**0*xst**10&
     +xp(223)*xs1**2*xs2**2*xst**8&
     +xp(224)*xs1**2*xs2**4*xst**6&
     +xp(225)*xs1**2*xs2**6*xst**4&
     +xp(226)*xs1**2*xs2**8*xst**2&
     +xp(227)*xs1**2*xs2**10*xst**0&
     +xp(228)*xs1**3*xs2**0*xst**9&
     +xp(229)*xs1**3*xs2**2*xst**7&
     +xp(230)*xs1**3*xs2**4*xst**5&
     +xp(231)*xs1**3*xs2**6*xst**3&
     +xp(232)*xs1**3*xs2**8*xst**1&
     +xp(233)*xs1**4*xs2**0*xst**8&
     +xp(234)*xs1**4*xs2**2*xst**6&
     +xp(235)*xs1**4*xs2**4*xst**4&
     +xp(236)*xs1**4*xs2**6*xst**2&
     +xp(237)*xs1**4*xs2**8*xst**0&
     +xp(238)*xs1**5*xs2**0*xst**7&
     +xp(239)*xs1**5*xs2**2*xst**5&
     +xp(240)*xs1**5*xs2**4*xst**3&
     +xp(241)*xs1**5*xs2**6*xst**1&
     +xp(242)*xs1**6*xs2**0*xst**6&
     +xp(243)*xs1**6*xs2**2*xst**4&
     +xp(244)*xs1**6*xs2**4*xst**2&
     +xp(245)*xs1**6*xs2**6*xst**0&
     +xp(246)*xs1**7*xs2**0*xst**5&
     +xp(247)*xs1**7*xs2**2*xst**3



       muX=v0+vp1+vp2+vp3

! Calculation of muY
xp=xpy

  re12      = xp(2)
      re32	= xp(248)
      alphae    = xp(3)*pi/180.0d0
      !
      aa1  = xp(4)
      aa2= xp(249)
      b1   = xp(5)
      b2   = xp(6)
      g1   = xp(7)
      g2   = xp(8)
      !

      !y1=1.0d+00-exp(-aa1*(r12-re12))
      !y2=1.0d+00-exp(-aa2*(r32-re32))
      xst=alpha-alphae
	y1=r12-re12
	y2=r32-re32  
	!
      xs1 = (y1+y2)*0.5d0
      xs2 = (y1-y2)*0.5d0


  v0= xp(  1)*xs1**0*xs2**0*xst**0

 vp1= xp(  9)*xs1**0*xs2**0*xst**1&
     +xp( 10)*xs1**1*xs2**0*xst**0&
     +xp( 11)*xs1**0*xs2**0*xst**2&
     +xp( 12)*xs1**0*xs2**2*xst**0&
     +xp( 13)*xs1**1*xs2**0*xst**1&
     +xp( 14)*xs1**2*xs2**0*xst**0&
     +xp( 15)*xs1**0*xs2**0*xst**3&
     +xp( 16)*xs1**0*xs2**2*xst**1&
     +xp( 17)*xs1**1*xs2**0*xst**2&
     +xp( 18)*xs1**1*xs2**2*xst**0&
     +xp( 19)*xs1**2*xs2**0*xst**1&
     +xp( 20)*xs1**3*xs2**0*xst**0&
     +xp( 21)*xs1**0*xs2**0*xst**4&
     +xp( 22)*xs1**0*xs2**2*xst**2&
     +xp( 23)*xs1**0*xs2**4*xst**0&
     +xp( 24)*xs1**1*xs2**0*xst**3&
     +xp( 25)*xs1**1*xs2**2*xst**1&
     +xp( 26)*xs1**2*xs2**0*xst**2&
     +xp( 27)*xs1**2*xs2**2*xst**0&
     +xp( 28)*xs1**3*xs2**0*xst**1&
     +xp( 29)*xs1**4*xs2**0*xst**0&
     +xp( 30)*xs1**0*xs2**0*xst**5&
     +xp( 31)*xs1**0*xs2**2*xst**3&
     +xp( 32)*xs1**0*xs2**4*xst**1&
     +xp( 33)*xs1**1*xs2**0*xst**4&
     +xp( 34)*xs1**1*xs2**2*xst**2&
     +xp( 35)*xs1**1*xs2**4*xst**0&
     +xp( 36)*xs1**2*xs2**0*xst**3&
     +xp( 37)*xs1**2*xs2**2*xst**1&
     +xp( 38)*xs1**3*xs2**0*xst**2&
     +xp( 39)*xs1**3*xs2**2*xst**0&
     +xp( 40)*xs1**4*xs2**0*xst**1&
     +xp( 41)*xs1**5*xs2**0*xst**0&
     +xp( 42)*xs1**0*xs2**0*xst**6&
     +xp( 43)*xs1**0*xs2**2*xst**4&
     +xp( 44)*xs1**0*xs2**4*xst**2&
     +xp( 45)*xs1**0*xs2**6*xst**0&
     +xp( 46)*xs1**1*xs2**0*xst**5&
     +xp( 47)*xs1**1*xs2**2*xst**3&
     +xp( 48)*xs1**1*xs2**4*xst**1&
     +xp( 49)*xs1**2*xs2**0*xst**4&
     +xp( 50)*xs1**2*xs2**2*xst**2&
     +xp( 51)*xs1**2*xs2**4*xst**0&
     +xp( 52)*xs1**3*xs2**0*xst**3&
     +xp( 53)*xs1**3*xs2**2*xst**1&
     +xp( 54)*xs1**4*xs2**0*xst**2&
     +xp( 55)*xs1**4*xs2**2*xst**0&
     +xp( 56)*xs1**5*xs2**0*xst**1&
     +xp( 57)*xs1**6*xs2**0*xst**0&
     +xp( 58)*xs1**0*xs2**0*xst**7&
     +xp( 59)*xs1**0*xs2**2*xst**5&
     +xp( 60)*xs1**0*xs2**4*xst**3&
     +xp( 61)*xs1**0*xs2**6*xst**1&
     +xp( 62)*xs1**1*xs2**0*xst**6&
     +xp( 63)*xs1**1*xs2**2*xst**4&
     +xp( 64)*xs1**1*xs2**4*xst**2&
     +xp( 65)*xs1**1*xs2**6*xst**0&
     +xp( 66)*xs1**2*xs2**0*xst**5&
     +xp( 67)*xs1**2*xs2**2*xst**3&
     +xp( 68)*xs1**2*xs2**4*xst**1&
     +xp( 69)*xs1**3*xs2**0*xst**4&
     +xp( 70)*xs1**3*xs2**2*xst**2&
     +xp( 71)*xs1**3*xs2**4*xst**0&
     +xp( 72)*xs1**4*xs2**0*xst**3&
     +xp( 73)*xs1**4*xs2**2*xst**1&
     +xp( 74)*xs1**5*xs2**0*xst**2&
     +xp( 75)*xs1**5*xs2**2*xst**0&
     +xp( 76)*xs1**6*xs2**0*xst**1&
     +xp( 77)*xs1**7*xs2**0*xst**0&
     +xp( 78)*xs1**0*xs2**0*xst**8&
     +xp( 79)*xs1**0*xs2**2*xst**6&
     +xp( 80)*xs1**0*xs2**4*xst**4&
     +xp( 81)*xs1**0*xs2**6*xst**2&
     +xp( 82)*xs1**0*xs2**8*xst**0&
     +xp( 83)*xs1**1*xs2**0*xst**7&
     +xp( 84)*xs1**1*xs2**2*xst**5&
     +xp( 85)*xs1**1*xs2**4*xst**3&
     +xp( 86)*xs1**1*xs2**6*xst**1&
     +xp( 87)*xs1**2*xs2**0*xst**6&
     +xp( 88)*xs1**2*xs2**2*xst**4&
     +xp( 89)*xs1**2*xs2**4*xst**2&
     +xp( 90)*xs1**2*xs2**6*xst**0&
     +xp( 91)*xs1**3*xs2**0*xst**5&
     +xp( 92)*xs1**3*xs2**2*xst**3&
     +xp( 93)*xs1**3*xs2**4*xst**1&
     +xp( 94)*xs1**4*xs2**0*xst**4&
     +xp( 95)*xs1**4*xs2**2*xst**2&
     +xp( 96)*xs1**4*xs2**4*xst**0&
     +xp( 97)*xs1**5*xs2**0*xst**3&
     +xp( 98)*xs1**5*xs2**2*xst**1&
     +xp( 99)*xs1**6*xs2**0*xst**2

 vp2= xp(100)*xs1**6*xs2**2*xst**0&
     +xp(101)*xs1**7*xs2**0*xst**1&
     +xp(102)*xs1**8*xs2**0*xst**0&
     +xp(103)*xs1**0*xs2**0*xst**9&
     +xp(104)*xs1**0*xs2**2*xst**7&
     +xp(105)*xs1**0*xs2**4*xst**5&
     +xp(106)*xs1**0*xs2**6*xst**3&
     +xp(107)*xs1**0*xs2**8*xst**1&
     +xp(108)*xs1**1*xs2**0*xst**8&
     +xp(109)*xs1**1*xs2**2*xst**6&
     +xp(110)*xs1**1*xs2**4*xst**4&
     +xp(111)*xs1**1*xs2**6*xst**2&
     +xp(112)*xs1**1*xs2**8*xst**0&
     +xp(113)*xs1**2*xs2**0*xst**7&
     +xp(114)*xs1**2*xs2**2*xst**5&
     +xp(115)*xs1**2*xs2**4*xst**3&
     +xp(116)*xs1**2*xs2**6*xst**1&
     +xp(117)*xs1**3*xs2**0*xst**6&
     +xp(118)*xs1**3*xs2**2*xst**4&
     +xp(119)*xs1**3*xs2**4*xst**2&
     +xp(120)*xs1**3*xs2**6*xst**0&
     +xp(121)*xs1**4*xs2**0*xst**5&
     +xp(122)*xs1**4*xs2**2*xst**3&
     +xp(123)*xs1**4*xs2**4*xst**1&
     +xp(124)*xs1**5*xs2**0*xst**4&
     +xp(125)*xs1**5*xs2**2*xst**2&
     +xp(126)*xs1**5*xs2**4*xst**0&
     +xp(127)*xs1**6*xs2**0*xst**3&
     +xp(128)*xs1**6*xs2**2*xst**1&
     +xp(129)*xs1**7*xs2**0*xst**2&
     +xp(130)*xs1**7*xs2**2*xst**0&
     +xp(131)*xs1**8*xs2**0*xst**1&
     +xp(132)*xs1**9*xs2**0*xst**0&
     +xp(133)*xs1**0*xs2**0*xst**10&
     +xp(134)*xs1**0*xs2**2*xst**8&
     +xp(135)*xs1**0*xs2**4*xst**6&
     +xp(136)*xs1**0*xs2**6*xst**4&
     +xp(137)*xs1**0*xs2**8*xst**2&
     +xp(138)*xs1**0*xs2**10*xst**0&
     +xp(139)*xs1**1*xs2**0*xst**9&
     +xp(140)*xs1**1*xs2**2*xst**7&
     +xp(141)*xs1**1*xs2**4*xst**5&
     +xp(142)*xs1**1*xs2**6*xst**3&
     +xp(143)*xs1**1*xs2**8*xst**1&
     +xp(144)*xs1**2*xs2**0*xst**8&
     +xp(145)*xs1**2*xs2**2*xst**6&
     +xp(146)*xs1**2*xs2**4*xst**4&
     +xp(147)*xs1**2*xs2**6*xst**2&
     +xp(148)*xs1**2*xs2**8*xst**0&
     +xp(149)*xs1**3*xs2**0*xst**7&
     +xp(150)*xs1**3*xs2**2*xst**5&
     +xp(151)*xs1**3*xs2**4*xst**3&
     +xp(152)*xs1**3*xs2**6*xst**1&
     +xp(153)*xs1**4*xs2**0*xst**6&
     +xp(154)*xs1**4*xs2**2*xst**4&
     +xp(155)*xs1**4*xs2**4*xst**2&
     +xp(156)*xs1**4*xs2**6*xst**0&
     +xp(157)*xs1**5*xs2**0*xst**5&
     +xp(158)*xs1**5*xs2**2*xst**3&
     +xp(159)*xs1**5*xs2**4*xst**1&
     +xp(160)*xs1**6*xs2**0*xst**4&
     +xp(161)*xs1**6*xs2**2*xst**2&
     +xp(162)*xs1**6*xs2**4*xst**0&
     +xp(163)*xs1**7*xs2**0*xst**3&
     +xp(164)*xs1**7*xs2**2*xst**1&
     +xp(165)*xs1**8*xs2**0*xst**2&
     +xp(166)*xs1**8*xs2**2*xst**0&
     +xp(167)*xs1**9*xs2**0*xst**1&
     +xp(168)*xs1**10*xs2**0*xst**0&
     +xp(169)*xs1**0*xs2**0*xst**11&
     +xp(170)*xs1**0*xs2**2*xst**9&
     +xp(171)*xs1**0*xs2**4*xst**7&
     +xp(172)*xs1**0*xs2**6*xst**5&
     +xp(173)*xs1**0*xs2**8*xst**3&
     +xp(174)*xs1**0*xs2**10*xst**1&
     +xp(175)*xs1**1*xs2**0*xst**10&
     +xp(176)*xs1**1*xs2**2*xst**8&
     +xp(177)*xs1**1*xs2**4*xst**6&
     +xp(178)*xs1**1*xs2**6*xst**4&
     +xp(179)*xs1**1*xs2**8*xst**2&
     +xp(180)*xs1**1*xs2**10*xst**0&
     +xp(181)*xs1**2*xs2**0*xst**9

 vp3= xp(182)*xs1**2*xs2**2*xst**7&
     +xp(183)*xs1**2*xs2**4*xst**5&
     +xp(184)*xs1**2*xs2**6*xst**3&
     +xp(185)*xs1**2*xs2**8*xst**1&
     +xp(186)*xs1**3*xs2**0*xst**8&
     +xp(187)*xs1**3*xs2**2*xst**6&
     +xp(188)*xs1**3*xs2**4*xst**4&
     +xp(189)*xs1**3*xs2**6*xst**2&
     +xp(190)*xs1**3*xs2**8*xst**0&
     +xp(191)*xs1**4*xs2**0*xst**7&
     +xp(192)*xs1**4*xs2**2*xst**5&
     +xp(193)*xs1**4*xs2**4*xst**3&
     +xp(194)*xs1**4*xs2**6*xst**1&
     +xp(195)*xs1**5*xs2**0*xst**6&
     +xp(196)*xs1**5*xs2**2*xst**4&
     +xp(197)*xs1**5*xs2**4*xst**2&
     +xp(198)*xs1**5*xs2**6*xst**0&
     +xp(199)*xs1**6*xs2**0*xst**5&
     +xp(200)*xs1**6*xs2**2*xst**3&
     +xp(201)*xs1**6*xs2**4*xst**1&
     +xp(202)*xs1**7*xs2**0*xst**4&
     +xp(203)*xs1**7*xs2**2*xst**2&
     +xp(204)*xs1**7*xs2**4*xst**0&
     +xp(205)*xs1**8*xs2**0*xst**3&
     +xp(206)*xs1**8*xs2**2*xst**1&
     +xp(207)*xs1**9*xs2**0*xst**2&
     +xp(208)*xs1**9*xs2**2*xst**0&
     +xp(209)*xs1**0*xs2**0*xst**12&
     +xp(210)*xs1**0*xs2**2*xst**10&
     +xp(211)*xs1**0*xs2**4*xst**8&
     +xp(212)*xs1**0*xs2**6*xst**6&
     +xp(213)*xs1**0*xs2**8*xst**4&
     +xp(214)*xs1**0*xs2**10*xst**2&
     +xp(215)*xs1**0*xs2**12*xst**0&
     +xp(216)*xs1**1*xs2**0*xst**11&    
     +xp(217)*xs1**1*xs2**2*xst**9&
     +xp(218)*xs1**1*xs2**4*xst**7&
     +xp(219)*xs1**1*xs2**6*xst**5&
     +xp(220)*xs1**1*xs2**8*xst**3&
     +xp(221)*xs1**1*xs2**10*xst**1&
     +xp(222)*xs1**2*xs2**0*xst**10&
     +xp(223)*xs1**2*xs2**2*xst**8&
     +xp(224)*xs1**2*xs2**4*xst**6&
     +xp(225)*xs1**2*xs2**6*xst**4&
     +xp(226)*xs1**2*xs2**8*xst**2&
     +xp(227)*xs1**2*xs2**10*xst**0&
     +xp(228)*xs1**3*xs2**0*xst**9&
     +xp(229)*xs1**3*xs2**2*xst**7&
     +xp(230)*xs1**3*xs2**4*xst**5&
     +xp(231)*xs1**3*xs2**6*xst**3&
     +xp(232)*xs1**3*xs2**8*xst**1&
     +xp(233)*xs1**4*xs2**0*xst**8&
     +xp(234)*xs1**4*xs2**2*xst**6&
     +xp(235)*xs1**4*xs2**4*xst**4&
     +xp(236)*xs1**4*xs2**6*xst**2&
     +xp(237)*xs1**4*xs2**8*xst**0&
     +xp(238)*xs1**5*xs2**0*xst**7&
     +xp(239)*xs1**5*xs2**2*xst**5&
     +xp(240)*xs1**5*xs2**4*xst**3&
     +xp(241)*xs1**5*xs2**6*xst**1&
     +xp(242)*xs1**6*xs2**0*xst**6&
     +xp(243)*xs1**6*xs2**2*xst**4&
     +xp(244)*xs1**6*xs2**4*xst**2&
     +xp(245)*xs1**6*xs2**6*xst**0&
     +xp(246)*xs1**7*xs2**0*xst**5&
     +xp(247)*xs1**7*xs2**2*xst**3



       muY=v0+vp1+vp2+vp3
      end subroutine DIPS

