!
! Defining potential energy function 
!
subroutine poten(force,f,r12,r32,theta)
   !
   integer,parameter           :: N = 99
   !
   double precision,intent(in) ::  force(N),r12,r32,theta
   double precision,intent(out)::  f

   double precision            :: xcos,v,v1,v2,v3,v4,v5,v6,v7,v8,alpha

   double precision            :: aa1,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3
   double precision            :: g1,g2,b1,b2,rhh,vhh,th1,th2
   integer                     :: n1
   !
   N1 = size(force)
   !
   alpha = theta
   !
   !write(6,"('r12,r32,alpha = ',3f20.8)") r12,r32,alpha
   !
   re12    = 1.335946936 ! 1.33756337958d0
   alphae  = 92.23305449*3.1415926535897932385/180.d0
   aa1     = 1.8d0
   !
   !re12    = 1.3359007d0
   !alphae  = 92.265883d0*3.1415926535897932385/180.d0
   !aa1     = 1.704d0
   !
   b1   = force(1)
   b2   = force(2)
   g1   = force(3)
   g2   = force(4)
   !
   rhh=sqrt(r12**2+r32**2-2.d0*r12*r32*cos(alpha))
   vhh=b1*exp(-g1*rhh)+b2*exp(-g2*rhh**2)
   !
   !th1 = 0.5d0*( 1.0d0-tanh( 5.0*( r12-0.0d0 ) ) )
   !th3 = 0.5d0*( 1.0d0-tanh( 5.0*( r32-0.0d0 ) ) )
   !
   th1 = 0.5d0*( 1.0d0-tanh( 2.0*( r12-0.0d0 ) ) )
   th3 = 0.5d0*( 1.0d0-tanh( 2.0*( r32-0.0d0 ) ) )
   !
   !vhh = vhh +(200000.d0*(alpha-alphae)**2+1e5)*th1/r12**12+(200000.d0*(alpha-alphae)**2+1e5)*th3/r32**12
   !
   !write(6,"('rhh,vhh = ',2f20.8)") rhh,vhh
   !
   ! calculate potential energy function values
   !
   y1=1.0d+00-exp(-aa1*(r12-re12))
   y2=1.0d+00-exp(-aa1*(r32-re12))
   !
   y3=(cos(alpha)-cos(alphae)) ! *exp(-1.0*(alpha-alphae)**2)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
   !write(6,"('y1,y2,y3 = ',3f20.8)") y1,y2,y3
   !
 v0 = force(5)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1& 
    + force(7)*y1**1*y2**0*y3**0& 
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2& 
    + force(9)*y1**1*y2**0*y3**1& 
    + force(9)*y1**0*y2**1*y3**1& 
    + force(10)*y1**1*y2**1*y3**0& 
    + force(11)*y1**2*y2**0*y3**0& 
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3& 
    + force(13)*y1**1*y2**0*y3**2& 
    + force(13)*y1**0*y2**1*y3**2& 
    + force(14)*y1**1*y2**1*y3**1& 
    + force(15)*y1**2*y2**0*y3**1& 
    + force(15)*y1**0*y2**2*y3**1& 
    + force(16)*y1**2*y2**1*y3**0& 
    + force(16)*y1**1*y2**2*y3**0& 
    + force(17)*y1**3*y2**0*y3**0& 
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then 
  v4 = force(18)*y1**0*y2**0*y3**4& 
    + force(19)*y1**1*y2**0*y3**3& 
    + force(19)*y1**0*y2**1*y3**3& 
    + force(20)*y1**1*y2**1*y3**2& 
    + force(21)*y1**2*y2**0*y3**2& 
    + force(21)*y1**0*y2**2*y3**2& 
    + force(22)*y1**2*y2**1*y3**1& 
    + force(22)*y1**1*y2**2*y3**1& 
    + force(23)*y1**2*y2**2*y3**0& 
    + force(24)*y1**3*y2**0*y3**1& 
    + force(24)*y1**0*y2**3*y3**1& 
    + force(25)*y1**3*y2**1*y3**0& 
    + force(25)*y1**1*y2**3*y3**0& 
    + force(26)*y1**4*y2**0*y3**0& 
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then 
  v5 = force(27)*y1**0*y2**0*y3**5& 
    + force(28)*y1**1*y2**0*y3**4& 
    + force(28)*y1**0*y2**1*y3**4& 
    + force(29)*y1**1*y2**1*y3**3& 
    + force(30)*y1**2*y2**0*y3**3& 
    + force(30)*y1**0*y2**2*y3**3& 
    + force(31)*y1**2*y2**1*y3**2& 
    + force(31)*y1**1*y2**2*y3**2& 
    + force(32)*y1**2*y2**2*y3**1& 
    + force(33)*y1**3*y2**0*y3**2& 
    + force(33)*y1**0*y2**3*y3**2& 
    + force(34)*y1**3*y2**1*y3**1& 
    + force(34)*y1**1*y2**3*y3**1& 
    + force(35)*y1**3*y2**2*y3**0& 
    + force(35)*y1**2*y2**3*y3**0& 
    + force(36)*y1**4*y2**0*y3**1& 
    + force(36)*y1**0*y2**4*y3**1& 
    + force(37)*y1**4*y2**1*y3**0& 
    + force(37)*y1**1*y2**4*y3**0& 
    + force(38)*y1**5*y2**0*y3**0& 
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then 
  v6 = force(39)*y1**0*y2**0*y3**6& 
    + force(40)*y1**1*y2**0*y3**5& 
    + force(40)*y1**0*y2**1*y3**5& 
    + force(41)*y1**1*y2**1*y3**4& 
    + force(42)*y1**2*y2**0*y3**4& 
    + force(42)*y1**0*y2**2*y3**4& 
    + force(43)*y1**2*y2**1*y3**3& 
    + force(43)*y1**1*y2**2*y3**3& 
    + force(44)*y1**2*y2**2*y3**2& 
    + force(45)*y1**3*y2**0*y3**3& 
    + force(45)*y1**0*y2**3*y3**3& 
    + force(46)*y1**3*y2**1*y3**2& 
    + force(46)*y1**1*y2**3*y3**2& 
    + force(47)*y1**3*y2**2*y3**1& 
    + force(47)*y1**2*y2**3*y3**1& 
    + force(48)*y1**3*y2**3*y3**0& 
    + force(49)*y1**4*y2**0*y3**2& 
    + force(49)*y1**0*y2**4*y3**2& 
    + force(50)*y1**4*y2**1*y3**1& 
    + force(50)*y1**1*y2**4*y3**1& 
    + force(51)*y1**4*y2**2*y3**0& 
    + force(51)*y1**2*y2**4*y3**0& 
    + force(52)*y1**5*y2**0*y3**1& 
    + force(52)*y1**0*y2**5*y3**1& 
    + force(53)*y1**5*y2**1*y3**0& 
    + force(53)*y1**1*y2**5*y3**0& 
    + force(54)*y1**6*y2**0*y3**0& 
    + force(54)*y1**0*y2**6*y3**0
 endif


 if (N>54) then 
 v7 = force(55)*y1**0*y2**0*y3**7& 
    + force(56)*y1**1*y2**0*y3**6& 
    + force(56)*y1**0*y2**1*y3**6& 
    + force(57)*y1**1*y2**1*y3**5& 
    + force(58)*y1**2*y2**0*y3**5& 
    + force(58)*y1**0*y2**2*y3**5& 
    + force(59)*y1**2*y2**1*y3**4& 
    + force(59)*y1**1*y2**2*y3**4& 
    + force(60)*y1**2*y2**2*y3**3& 
    + force(61)*y1**3*y2**0*y3**4& 
    + force(61)*y1**0*y2**3*y3**4& 
    + force(62)*y1**3*y2**1*y3**3& 
    + force(62)*y1**1*y2**3*y3**3& 
    + force(63)*y1**3*y2**2*y3**2& 
    + force(63)*y1**2*y2**3*y3**2& 
    + force(64)*y1**3*y2**3*y3**1& 
    + force(65)*y1**4*y2**0*y3**3& 
    + force(65)*y1**0*y2**4*y3**3& 
    + force(66)*y1**4*y2**1*y3**2& 
    + force(66)*y1**1*y2**4*y3**2& 
    + force(67)*y1**4*y2**2*y3**1& 
    + force(67)*y1**2*y2**4*y3**1& 
    + force(68)*y1**4*y2**3*y3**0& 
    + force(68)*y1**3*y2**4*y3**0& 
    + force(69)*y1**5*y2**0*y3**2& 
    + force(69)*y1**0*y2**5*y3**2& 
    + force(70)*y1**5*y2**1*y3**1& 
    + force(70)*y1**1*y2**5*y3**1& 
    + force(71)*y1**5*y2**2*y3**0& 
    + force(71)*y1**2*y2**5*y3**0& 
    + force(72)*y1**6*y2**0*y3**1& 
    + force(72)*y1**0*y2**6*y3**1& 
    + force(73)*y1**6*y2**1*y3**0& 
    + force(73)*y1**1*y2**6*y3**0& 
    + force(74)*y1**7*y2**0*y3**0& 
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then 
 v8 = force(75)*y1**0*y2**0*y3**8& 
    + force(76)*y1**1*y2**0*y3**7& 
    + force(76)*y1**0*y2**1*y3**7& 
    + force(77)*y1**1*y2**1*y3**6& 
    + force(78)*y1**2*y2**0*y3**6& 
    + force(78)*y1**0*y2**2*y3**6& 
    + force(79)*y1**2*y2**1*y3**5& 
    + force(79)*y1**1*y2**2*y3**5& 
    + force(80)*y1**2*y2**2*y3**4& 
    + force(81)*y1**3*y2**0*y3**5& 
    + force(81)*y1**0*y2**3*y3**5& 
    + force(82)*y1**3*y2**1*y3**4& 
    + force(82)*y1**1*y2**3*y3**4& 
    + force(83)*y1**3*y2**2*y3**3& 
    + force(83)*y1**2*y2**3*y3**3& 
    + force(84)*y1**3*y2**3*y3**2& 
    + force(85)*y1**4*y2**0*y3**4& 
    + force(85)*y1**0*y2**4*y3**4& 
    + force(86)*y1**4*y2**1*y3**3& 
    + force(86)*y1**1*y2**4*y3**3& 
    + force(87)*y1**4*y2**2*y3**2& 
    + force(87)*y1**2*y2**4*y3**2& 
    + force(88)*y1**4*y2**3*y3**1& 
    + force(88)*y1**3*y2**4*y3**1& 
    + force(89)*y1**4*y2**4*y3**0& 
    + force(90)*y1**5*y2**0*y3**3& 
    + force(90)*y1**0*y2**5*y3**3& 
    + force(91)*y1**5*y2**1*y3**2& 
    + force(91)*y1**1*y2**5*y3**2& 
    + force(92)*y1**5*y2**2*y3**1& 
    + force(92)*y1**2*y2**5*y3**1& 
    + force(93)*y1**5*y2**3*y3**0& 
    + force(93)*y1**3*y2**5*y3**0& 
    + force(94)*y1**6*y2**0*y3**2& 
    + force(94)*y1**0*y2**6*y3**2& 
    + force(95)*y1**6*y2**1*y3**1& 
    + force(95)*y1**1*y2**6*y3**1& 
    + force(96)*y1**6*y2**2*y3**0& 
    + force(96)*y1**2*y2**6*y3**0& 
    + force(97)*y1**7*y2**0*y3**1& 
    + force(97)*y1**0*y2**7*y3**1& 
    + force(98)*y1**7*y2**1*y3**0& 
    + force(98)*y1**1*y2**7*y3**0& 
    + force(99)*y1**8*y2**0*y3**0& 
    + force(99)*y1**0*y2**8*y3**0
endif

     th1 = 0.5d0*( 1.0d0-tanh( 0.00040d0*( v0+v1+v2-45000.0d0 ) ) )

    !th1 = 0.5d0*( 1.0d0-tanh( 0.0001d0*( v0+v1+v2-50000.0d0 ) ) )


    !f=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh
    !
    f=(v0+v1+v2)+(v3+v4+v5+v6+v7+v8)*th1+vhh
    !
    !if (v0+v1+v2>50000.0) f = v0+v1+v2
    !
    !if (f<0.0) then 
    ! write(6,"(2g18.6,2x,g18.6,2x,11g20.8)") r12,r32,theta/3.1415926535897932385*180.d0,f,v0,v1,v2,v3,v4,v5,v6,v7,v8,vhh
    !endif 


end subroutine poten

