! This program calculates the dipole of:-
! T. van Mourik, G. J. Harris, O. L. Polyansky, J. Tennyson, 
! A. G. Csaszar and P. J. Knowles, J. Chem. Phys. 115, 3706 (2001).
! Written by by GJH March 16th 2001.
!
! 
! Fits are 3*3*5 and 3*3*7 fits to x and z respectively

!      A test main program
 
!      implicit double precision (a-h,o-z)
!      cosg=cos(1.3962634015952)
!      bigr= 2.16667
!      smallr=2.13333 
!      icomp=0

!      should then give dipd=0.499620862d0 for icomp=1
!      for icomp=0 dipd=-0.00166465088d0


!      smallr=2.2 
!      bigr=3.05172d0
!      cosg=-1.0

!      do 10 i=1, 201

!      icomp=0
!      call dipd(dipole0, smallr, bigr, cosg,  icomp )
!      icomp=1
!      call dipd(dipole1, smallr, bigr, cosg,  icomp )

!      write(6,20) cosg, bigr, smallr, dipole0*2.5417662, dipole1*2.5417662
!      cosg=cosg+0.009999d0

! 10    continue
! 20    format(5d12.4)

!      end

! **************************************************************************

      subroutine dipd(dipole, smallr, bigr, cosg,  icomp )

! Jacobi coordinates are used.
! Units of lengths in a0 and dipole in au.
!
! smallr is the C to N bond length. 

! bigr is the H to CN center of mass distance.

! cosg is the cosine of the angle between the bigr and smallr, 
!      an angle of 0 corresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH

! icomp selects the x or z component of the dipole, 
!     if icomp=0 then returns z component of dipole
!     if icomp=1 then the x component of the dipole is returned. Dipole is the output 

! Z dipole component is along C to H bond.
! Xdipole component is perpendicular to z component, H bends in xz plane.


      implicit double precision (a-h,o-z)
      if (icomp .ne. 0 .and. icomp .ne. 1) then 
         write (6, *) "Invalid value for NU"
      else
         if (icomp .eq. 1) dipole=x_comp_dpl(cosg, bigr, smallr)
         if (icomp .eq. 0) dipole=z_comp_dpl(cosg, bigr, smallr)
      endif
!      write(6,*) dipole
      return
      end

! *************************************************************************
      function z_comp_dpl(cosg, bigr, smallr)

      implicit double precision (a-h,o-z)
      common /Z_comp/ coeffz(63)
      dimension B(3), C(3)
! code to calculate the z component of the dipole moment.

! Standard Deviation = 0.008 D

      call set_Z_comp

      B(1) = -0.33194847405555984d0
      B(2) = 0.26846609653383304d0
      B(3) = 1.4784268008920333d0

      C(1) = 0.41155209135581288d0
      C(2) = 0.41038551093500483d0
      C(3) = 0.58913645613028966d0

      temp=0.0
      icoefref=1
      
      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=0, 6
               call plgndr(k, 0, cosg, apolynom)
           temp=temp+(coeffz(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      z_comp_dpl= temp
      return
      end


! ************************************************************************

      function x_comp_dpl(cosg, bigr, smallr)

      implicit double precision (a-h,o-z)
      common /X_comp/ coeffx(45)
      dimension B(3), C(3)
! code to calculate the x component of the dipole moment.
      call set_X_comp

! Standard Deviation = 0.002 D

      temp=0.0
      icoefref=1

      B(1) = -0.47194808781916352d0
      B(2) = -0.54677773064856535d0
      B(3) = -0.46378831330917923d0

      C(1) = 0.0060793434403400061d0
      C(2) = -0.011738771795116187d0
      C(3) = 3.9847070586924296d0

      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=1, 5
               call plgndr(k, 1, cosg, apolynom)
           temp=temp+(coeffx(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      x_comp_dpl= temp

      return
      end


! ************************************************************************

      subroutine plgndr(l, m, x, pmm)

! subroutine to calculate the asociated legenre polynomials, translated drectly from a C 
! version to fortran version. The C version comes from Numerical recipies in C,  Press et al. 

      double precision fact, pll, pmm, pmmp1, somx2, x
      integer i, ll, l, m

      if (m .LT. 0 .OR. m .GT. l .OR. (x*x) .GT. 1.0) then
         write(6,*)  "!ERROR! bad arguments in routine plgndr"
      endif

      pmm = 1.0

      if (m .GT. 0) then
   
         somx2=sqrt((1.0 - x) * (1.0 + x))
         fact = 1.0
         do i=1, m
            pmm = pmm*(-fact * somx2)
            fact = fact+ 2.0
         enddo
      endif

      if (l .NE. m) then
         pmmp1 = x * (2 * m + 1) * pmm

         if (l .EQ. (m + 1)) then
            pmm=pmmp1
         else
       
            do ll = m+2, l
           
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
            enddo
            pmm=pll
         endif
      endif
      return
      end

! ****************************************************************
      subroutine set_X_comp
      implicit double precision (a-h,o-z)
      common /X_comp/ coeffx(45)

      coeffx(1)=3307.26744745775115d0
      coeffx(2)=31514.9042500608293d0
      coeffx(3)=-20148.594904970726d0
      coeffx(4)=546.081750420754701d0
      coeffx(5)=7723.26526967422282d0
      coeffx(6)=-584.587960975932364d0
      coeffx(7)=-32963.3776253884116d0
      coeffx(8)=21191.7974142975256d0
      coeffx(9)=-590.053047331642632d0
      coeffx(10)=-7990.27903779592029d0
      coeffx(11)=-0.031796628743699731d0
      coeffx(12)=-0.00940828117939172576d0
      coeffx(13)=0.00955414571501609397d0
      coeffx(14)=-0.00072965350671948823d0
      coeffx(15)=-0.00340116044617136415d0
      coeffx(16)=-793.754144943871578d0
      coeffx(17)=-3486.98023111115815d0
      coeffx(18)=2421.86358491171674d0
      coeffx(19)=-146.240552443951966d0
      coeffx(20)=-908.472444995532274d0
      coeffx(21)=521.168738604204326d0
      coeffx(22)=3646.28887952151323d0
      coeffx(23)=-2546.78588888870825d0
      coeffx(24)=154.568740785581287d0
      coeffx(25)=939.567725418272133d0
      coeffx(26)=0.00359417943193839101d0
      coeffx(27)=0.00108249483506098995d0
      coeffx(28)=-0.00108048813729769784d0
      coeffx(29)=0.000123617108271128898d0
      coeffx(30)=0.000402299503589065202d0
      coeffx(31)=-2556.9099117606792d0
      coeffx(32)=-28050.9552366028232d0
      coeffx(33)=17760.1943565608219d0
      coeffx(34)=-413.178836891090924d0
      coeffx(35)=-6829.82090935623486d0
      coeffx(36)=106.994087651491977d0
      coeffx(37)=29341.1779083554478d0
      coeffx(38)=-18680.3377835709173d0
      coeffx(39)=449.374983259753014d0
      coeffx(40)=7066.18946387678389d0
      coeffx(41)=0.0282356308831904958d0
      coeffx(42)=0.0083400647971135194d0
      coeffx(43)=-0.00848293720114831184d0
      coeffx(44)=0.000614212002770376261d0
      coeffx(45)=0.00300582095992624047d0


      return
      end


! ***************************************************
      subroutine set_Z_comp
      implicit double precision (a-h,o-z)
      common /Z_comp/ coeffz(63)

      coeffz(1)=-1039.93729474448763d0
      coeffz(2)=353.933029028493396d0
      coeffz(3)=-469.273752106125326d0
      coeffz(4)=5053.36372859932299d0
      coeffz(5)=-174.791252669703362d0
      coeffz(6)=-1790.40179998532812d0
      coeffz(7)=-1737.30128936150733d0
      coeffz(8)=1036.43863494317905d0
      coeffz(9)=-351.840317729231717d0
      coeffz(10)=466.654729016545031d0
      coeffz(11)=-5035.72026003557988d0
      coeffz(12)=174.888204603534292d0
      coeffz(13)=1783.68882790358092d0
      coeffz(14)=1730.16793139302162d0
      coeffz(15)=4.13996041131474449d0
      coeffz(16)=-1.77817321839038074d0
      coeffz(17)=2.46984600203225505d0
      coeffz(18)=-20.8620301949723831d0
      coeffz(19)=0.248557957243695296d0
      coeffz(20)=7.69363017734038415d0
      coeffz(21)=7.82369848101642729d0
      coeffz(22)=-170.934815933980384d0
      coeffz(23)=1038.00304354550415d0
      coeffz(24)=350.15335220496472d0
      coeffz(25)=-1945.5046107913207d0
      coeffz(26)=-192.883669686284673d0
      coeffz(27)=646.020865093223846d0
      coeffz(28)=464.968280321535981d0
      coeffz(29)=170.495969427259945d0
      coeffz(30)=-1035.59570071766628d0
      coeffz(31)=-348.571134507640448d0
      coeffz(32)=1938.59278642478573d0
      coeffz(33)=191.925410021956726d0
      coeffz(34)=-643.632109747355521d0
      coeffz(35)=-463.047969228713786d0
      coeffz(36)=0.60555784179134213d0
      coeffz(37)=-3.61463889815628707d0
      coeffz(38)=-1.65622524995037699d0
      coeffz(39)=8.06544948486997637d0
      coeffz(40)=0.985792239355321115d0
      coeffz(41)=-2.74472229985483116d0
      coeffz(42)=-2.10370116762671985d0
      coeffz(43)=11.6158010664658819d0
      coeffz(44)=-1.34518644018171973d0
      coeffz(45)=-30.6940205962599487d0
      coeffz(46)=12.6697949751099632d0
      coeffz(47)=20.7241897553344583d0
      coeffz(48)=-4.33083732019109644d0
      coeffz(49)=-9.04791261638706548d0
      coeffz(50)=-11.5726193461914656d0
      coeffz(51)=1.34352911201902514d0
      coeffz(52)=30.5762810875374454d0
      coeffz(53)=-12.6246369560126679d0
      coeffz(54)=-20.6427638916496952d0
      coeffz(55)=4.31545198368715627d0
      coeffz(56)=9.01190335420644576d0
      coeffz(57)=-0.0494259888418515414d0
      coeffz(58)=0.00383241424113694307d0
      coeffz(59)=0.132747818418275613d0
      coeffz(60)=-0.0526642790873265692d0
      coeffz(61)=-0.0909424933648184257d0
      coeffz(62)=0.0179744240449498035d0
      coeffz(63)=0.0400004758413293306d0


       return
       end
! ****************************************************************
