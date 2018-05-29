      program DVR3DRJZ
      call dvr3d
      stop
      end

!######################################################################
      subroutine dvr3d
!
!     program               d v r 3 d r j z
!                           ~~~~~~~~~~~~~~~
!     should be cited as:
!         *********** add here ********************
!     program to do ro-vibrational calculations on triatomic systems
!     using general length, length, angle coordinates and a choice of
!     embeddings, in a multidimensional dvr in choice of coordinate
!     orders.
!     there are various options:
!     (a) calculation at a frozen angle                       ztwod  = t
!     (b) atom - rigid diatom calculations                    ncoord = 2
!     (c) triatomic calculations in all 3 dimensions          ncoord = 3
!     additionally, for Radau coordinates there are options to:
!     (d) place the z-axis along the bisector                 zbisc  = t
!     (e) place the z-axis perpendicular to the molecule plane zperp = t
!     see:
!
!     j.r.henderson & j.tennyson, chem.phys.lett., 173, 133 (1990).
!     j.r.henderson, phd thesis, university of london (1990).
!     j.r.henderson, j.tennyson & b.t. sutcliffe, j.chem.phys. 98, 7191 (1993).
!     j. tennyson & b.t. sutcliffe, int. j. quantum chem. 42, 941 (1992).
!
!     use as follows:
!     comments on namelist parameters (& defaults) in block data
!     input data read in subroutines insize & setcon
!     the program needs the following subroutines:
!            subroutine potv(v,r1,r2,x) should return the potential v
!            in hartrees for the point x = cos(theta) and
!            bondlengths r1 & r2 in bohr.
!     2. dsyev to do in core diagonalisation (lapack f90 routine).
!     the programme works in **** atomic units ***** :
!     1. the subroutine potv should return the potential in
!        hartrees for distances in bohr.
!     2. all input in setcon is in hartree or bohr except
!     3. the nuclear masses are read in atomic mass units & converted.
!     4. the eigenvalues are printed in both hartree & wavenumbers.
!
!     Rewritten into fortran 95 by Max Kostin and Jonathan Tennyson
      implicit logical (z)
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      namelist/prt/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common/timing/itime0

      write(6,1000)
 1000 format(5x,'Program DVR3DRJZ (version of June 2010)')

!     read in namelist input data (defaults in block data)
      read(5,prt)

      call SYSTEM_CLOCK(itime0,irate2,imax2)

!     read in control parameters of problem.
      call insize

!     now do the calculation
      call ccmain

      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/) 
      stop
      end

!##############################################################################
      block data
!     stores defaults for namelist parameters  
      implicit logical (z)

!  zpham[f] = t requests printing of the hamiltonian matrix.
!  zprad[f] = t requests printing of the radial matrix elements.
!  zp1d [f] = t requests printing of the results of 1d calculations.
!  zp2d [f] = t requests printing of the results of 2d calculations.
!  zpmin[f] = t requests only minimal printing.
!  zpvec[f] = t requests printing of the eigenvectors.
!  zlmat[f] = t requests printing of the L-matrix.
!  zcut[f]  = t final dimension selected using an energy cut-off given
!             by emax2.
!            = f final dimension determined by nham3.
!  zmors1[t]= t use morse oscillator-like functions for r_1 coordinate;
!           = f use spherical oscillator functions.
!  zmors2[t]= t use morse oscillator-like functions for r_2 coordinate;
!           = f use spherical oscillator functions.
!  zrot[t]  = t do vibrational part of rotational calculation by
!               looping over k
!  zperp[f] = f  z in plane calculation (uses JHMAIN or NFMAIN)
!           = t  z perpendicular calculation (use MKMAIN)
!  zbisc[f] = t  bisector embedding (use NFMAIN), set by the program
!  zembed[t]= t z axis is along r2, = f z axis is along r1.
!               only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.
!  zlin     = t forces suppresion of functions at last dvr point
!               (zbisc=t only).
!  zladd[t] = t NALF kept constant as k increases
!           = f NALF decreases with k (=f has a bug), (only if zrot = .true.)
!  ztwod[f] = t perform 2D calculation only at specified grid point.
!  zvec[f]  = t store the eigenvectors from all the parts of the calculation
!             (1d,2d and 3d) on stream iout2.
!             further information relating to this (arrays iv1 and iv2) is
!             stored on stream iout1.
!  zall[f]  = t requests no truncation of the intermediate solution.
!  ztheta[t]= t let theta be first in the order of solution;
!           = f let theta be last in the order of solution,
!             (used if idia > -2 only).
!  zr2r1[t] = t let r_2 come before r_1 in the order of solution;
!           = f let r_1 come before r_2 in the order of solution.
!             (used if idia > -2 only).
!  ztran[f]= t perform the transformation of the solution coefficients
!             to the expression for the wavefunction amplitudes at the grid
!             points. store the data on stream iwave, ztran = t
!             automatically sets zvec = t for idia > -2.
!  zquad2[t]= t use the dvr quadrature approximation for the integrals of
!             the r_2^{-2} matrix, and hence make its dvr transformation
!             diagonal.
!           = f evaluate the r_2^{-2} integrals fully and perform the
!             full dvr transformation on them.
!             note that zquad2 = f is not implemented for zmors2 = t
!             or for ztheta = f.
!  zdiag[t] = f do not do final diagonalisation, instead the final Hamiltonian
!             matrix is written on units IDIAG1 and IDIAG2. 
!  zpfun[t] = t store energy levels on stream ilev
!  zplot[f] = t prints all needed wavefunctions for plotting purposes
!  ilev[14]     stream for final eigenvalues (formatted).
!  ieigs1[7]    stream for eigenvalues of the 1d solutions.
!  ivecs1[3]    stream for eigenvectors of the 1d solutions.
!  ieigs2[2]    stream for eigenvalues of the 2d solutions.
!  ivecs2[4]    stream for eigenvectors of the 2d solutions.
!  ivint[17]    a scratch stream used for storing intermediate vectors in
!               building the final hamiltonian.
!  intvec[16]   a scratch stream for intermediate storage of the 2d vectors.
!  iband[15]    scratch file used for storing bands of the final hamiltonian.
!  iout1[24]    stream for arrays iv1 and iv2, which record the sizes of
!               the truncated vctors. used when zvec = t.
!  iout2[25]    stream for the 1d, 2d and 3d vectors for use when zvec = t.
!  iwave[26]    stores the wavefunction amplitudes at the grid points when
!               ztran = t.
!
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      data zpham/.false./,zprad/.false./,zpvec/.false./,zrot/.true./,&
           zladd/.true./,zembed/.true./,zmors2/.true./,&
           zpmin/.false./,zvec/.false./,zquad2/.true./,zcut/.false./,&
           zdiag/.true./,zlmat/.false./,zall/.false./,zplot/.false./,&
           zp1d/.false./,zp2d/.false./,zr2r1/.true./,ztheta/.true./,&
           zmors1/.true./,ztran/.false./,ztwod/.false./,zperp/.false./,&
           ieigs1/7/,ivecs1/3/,ieigs2/2/,ivecs2/4/,ivint/17/,&
           iband/15/,intvec/16/,idiag1/20/,idiag2/21/,iout1/24/,&
           iout2/25/,iwave/26/,zlin/.false./,zpfun/.false./,ilev/14/,idip/69/&
           idipd/70/,iwvpb/71/,iplot/40/
      end

!############################################################################
      subroutine insize

!     set up common /size/ & write control parameters of problem  

      implicit real*8 (a-h,o-y), logical (z)

!     common /size/ stores control parameters for the problem
!     npnt1: number of (gauss-laguerre) dvr points in r1
!     npnt2: number of (gauss-laguerre) dvr points in r2
!     nmax1: max order of r1 radial laguerre polynomial ( = npnt1-1)
!     nmax2: max order of r2 radial laguerre polynomial ( = npnt2-1)
!     maxleg:max order of angular legendre polynomial   ( = nalf -1)
!     max2d :upper bound on size of intermediate 2d hamiltonian
!     max2d2:max2d for smaller block (zbisc=t only)
!     max3d :upper bound on size of full 3d hamiltonian
!     max3d2:max3d for smaller block (zbisc=t only)
!     nalf : number of (gauss-legendre) dvr points in theta
!     idvr : number of unique dvr points
!     jrot : jrot is total angular momentum of the molecule
!     kmin : zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
!                    kmin=2 loop over both sym & anti-sym (zbisc=t only)
!            zrot=f, kmin=fixed value of k
!     npnt : max(npnt1,npnt2)
!     idia : = 1 scattering coordinates heteronuclear diatomic
!            = 2 scattering coordinates homonuclear diatomic
!            =-1 radau      coordinates hetronuclear diatomic
!            =-2 radau      coordinates homonuclear  diatomic
!            = 0 radau      coordinates with the z axis perpendicular to
!                the molecular plane.
!     ipar : parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
!     nlim1: =nmax1+1*(nmax1+1+1)/2
!     nlim2: =nmax2+1*(nmax2+1+1)/2
!     npnta: the number of dvr points in
!            the coordinate to be treated first in the dvr successive
!            diagonalisation-truncation procedure
!     npntb: the number of dvr points in the coordinate to come second
!     npntc: the number of dvr points in the coordinate to come last
!     ndima: set equal to npnta at the start - used for dimensioning
!     ndimb: set equal to npntb at the start - used for dimensioning
!     ndimc: set equal to npntc at the start - used for dimensioning
!     neval: number of eigenvalues which have to actually be supplied
!            as output
!     ncoord: number of vibrational coordinates explicitly considered
!     if (ncoord .ne. 3) some of the above are dummies, see below.

      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      character(len=8) title(9)
!     read in control parameters of problem:

!     ncoord = 2: atom-diatom problem with diatom rigid
!     ncoord = 3: full 3-d triatomic problem
      read(5,5) ncoord
    5 format(15i5)
!     other paramters: see comments above on /size/ for meaning
!     if ncoord=2: also need lmax,lpot,idia,kmin
!     if ncoord=3: all paramters required

      read(5,5) npnt2,jrot,neval,nalf,max2d,max3d,idia,&
                kmin,npnt1,ipar,max3d2
      
      zbisc = .false.
      if (jrot .eq. 0) then
         zembed = .true.
         kmin = 0
         zrot = .false.
      endif
      if (idia .eq. -2) then
         ncoord=3
         npnt1=npnt2
         idvr=nalf
         zmors1=zmors2
         ztheta=.false.
         if (jrot .ne. 0.and.zperp.ne..true.) zbisc=.true.
         if (zperp.and.zrot) then
            kmin=1
            if (jrot .eq. 0) kmin=0
         endif
      else
          if (ztran) zvec=.true.
      endif

      if (ztwod) then
         ncoord = 3
         idia = 1
         nalf = 1
         nalf2= 1
         idvr = 1
         zrot=.false.
         ztheta=.false.
         max3d=max2d
         neval=min(max2d,neval)
         goto 887
      endif

!     the gauss-legendre integration points are only acually
!     calculated for the half-range:
      nalf2 = (nalf+1)/2
!     are we doing atom, atom-rigid diatom or the full problem?
      if (ncoord .eq. 2) then
!        atom - rigid diatom case, set dummy /size/
         if (zr2r1) then
            write(6,1010)
 1010 format(/10x,'atom - rigid diatom vibrational analysis with:'/)
            npnt1 = 1
            nmax1 = 0
         else
            write(6,1011)
 1011 format(/5x,'Fixed - r2 vibrational analysis with:'/)
            npnt2 = 1
            nmax2 = 0
         endif
      else
!        full case: print data about extra radial basis
         ncoord  = 3
         write(6,1020) npnt1
 1020    format(/5x,'Full triatomic vibrational problem with'/&
                /5x,i5,3x,'radial r1 dvr points used,')
      endif
      
  887 continue
!     maximum order of polynomials;
      maxleg = nalf - 1
      nmax2  = npnt2- 1
      nmax1  = npnt1- 1

      if (ztheta) then
          npnta = nalf/max(1,idia)
          if (zr2r1) then
             npntb = npnt2
             npntc = npnt1
          else
             npntb = npnt1
             npntc = npnt2
          endif
      else
          if (zr2r1) then
             npnta = npnt2
             npntb = npnt1
          else
             npnta = npnt1
             npntb = npnt2
          endif
          npntc = nalf/max(1,idia)
      endif

      if (idia .eq. -2) then
         if (zrot .and. (jrot+kmin).gt.1) then
            max2d  = npnt1*(npnt1+1)/2
            max2d2 = npnt1*(npnt1-1)/2
         else
            max2d = npnt1*(npnt1+1-(ipar * 2))/2
            max2d2 = 0
         endif
         if (.not. zall)  max3d=min(max3d,max2d*nalf)
         if (zall) max3d=max2d*nalf
         if (zrot) then
            if (max3d2 .le. 0) max3d2=max3d
            max3d2=min(max3d,max2d2*nalf,max3d2)
         endif
      else
         if (zall) then
            max2d = npnta*npntb
            max3d = max2d*npntc
         else
            max2d=min(max2d,npnta*npntb)
            max3d=min(max3d,npnta*npntb*npntc)
         endif
      endif

      if (neval .le. 0) neval = 10
      neval=min(max3d,neval)
      if (ztwod) write(6,1023) npnt1
 1023 format(/5x,i5,3x,'radial r1 dvr points used,')
      if (ncoord .eq. 3) write(6,1030) npnt2,2,nalf,neval,max3d
      if (ncoord .eq. 2 .and. zr2r1)&
          write(6,1030) npnt2,2,nalf,neval,max3d
      if(ncoord .eq. 2 .and. .not. zr2r1)&
          write(6,1030) npnt1,1,nalf,neval,max3d
 1030 format(5x,i5,3x,'radial r',i1,' dvr points used,',&
            /5x,i5,3x,'angular dvr points used, with',&
            /5x,i5,3x,'lowest eigenvectors chosen from',&
            /5x,i5,3x,'maximum dimension secular problem'/)
      if(idia .eq. 2 .and. zperp) then
        write(6,1035)
        stop
 1035  format(/5x,'STOP!!!  ZPERP should be .false. for IDIA=2')
      endif      
! new bit regarding the zplot option:
      if (zplot) then
         if (.not.zperp) write(6,1038)
1038 format(5x,'As of 29th june 2006 zplot works only combined with',&
           /5x,' zperp=TRUE. Apologies for any inconveniences cause.',&
           /5x,' pb. '/)
         read(5,1037)nploti,nplotf,ithre,npth
1037 format(4I5)
         write(6,1036)nploti,nplotf,ithre,npth
 1036 format(5x,'Requested wavefunctions for plotting:', &
            /5x,'waves required from n. ',2x,i5,&
            /5x,'               to   n. ',2x,i5,&
            /5x,'cutting threshold set to',3x,i3,&
            /5x,'finally number of theta points selected:',3x,i3/)
         if (npth.lt.nalf) then
            npth=nalf
       write(6,*)'!#!#: n theta points found too small: increased to ',nalf
         end if
         if (nplotf.eq.0) nplotf=neval
      nplotf=min(nplotf,neval)
      if (nploti.gt.nplotf) then
         write(6,*)'Inconsinstency in call to wf plotting program'
         write(6,*)'Requested to print from wf: '
         write(6,*)'[ ',nploti,' - ',nplotf,' ] '
         write(6,*)'(number of vibrational levels found ',neval,' )'
         stop
      else
      nplot=nplotf-nploti+1
      end if
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      read(5,500)   title
  500 format(9a8)
      write(6,1040) title
 1040 format(5x,'Title: ',9a8/)
      if (ncoord .eq. 3) then
         if (zmors1)  write(6,1050) 1
         if (.not. zmors1) write(6,1060) 1
      endif
      if (ncoord .eq. 2 .and. .not. zr2r1) then
         if (zmors1)  write(6,1050) 1
         if (.not. zmors1) write(6,1060) 1
      else
         if (zmors2) write(6,1050) 2
         if (.not. zmors2) write(6,1060) 2
 1050 format(5x,'Morse oscillators used for r',i1,' basis')
 1060 format(5x,'Spherical oscillators used for r',i1,' basis')
         if (zquad2) then
            write(6,1051)
 1051 format(/5x,'Quadrature approximation used for the r2**(-2) terms'/)
         else
            write(6,1052)
 1052 format(/5x,'Quadrature approximation abandoned for r2**(-2) terms'/)
         endif
      endif
      if (zall) write(6,1067)
 1067 format(/5x,'All solutions from lower dimensions have been used')
      if (ztheta) then
         if (zr2r1) write(6,1042)
 1042    format(5x,'Problem solved in the order: theta -> r2 -> r1')
         if (.not. zr2r1) write(6,1043)
 1043    format(5x,'Problem solved in the order: theta -> r1 -> r2')
      else
         if (.not.ztwod) then
            if (zr2r1) write(6,1044)
 1044    format(5x,'Problem solved in the order: r2 -> r1 -> theta')
            if (.not. zr2r1) write(6,1045)
 1045    format(5x,'Problem solved in the order: r1 -> r2 -> theta')
         else
            if (zr2r1) write(6,1046)
 1046    format(5x,'Problem solved in the order: r2 -> r1')
            if (.not. zr2r1) write(6,1047)
 1047    format(5x,'Problem solved in the order: r1 -> r2')
         endif
      endif
      if (zcut) then
         write(6,1061)
 1061    format(5x,'Final basis selected using energy cut-off')
      else
         if (zrot .and. zbisc .and. (jrot+kmin).gt.1) then
            write(6,1062) max3d,max3d2
 1062    format(/5x,'Final basis comprises',i5,' lowest functions',&
                   ' for even parity hamiltonian'/&
                 5x,'Final basis comprises',i5,' lowest functions',&
                   ' for odd  parity hamiltonian')
         else
           if (.not. zcut) write(6,1064) max3d
 1064      format(5x,'final basis comprises',i5,' lowest functions')
         endif
      endif
      if (.not.ztwod) then
         if (zlmat) write(6,1065)
 1065    format(/5x,'Printing of L-matrix requested')
         if (.not.zlmat) write(6,1066)
 1066    format(/5x,'Printing of L-matrix not requested')
      endif
      if (zp1d) write(6,1071)
 1071 format(5x,'Printing of 1d eigenvalues requested')
      if (zp2d) write(6,1072)
 1072 format(5x,'Printing of 2d eigenvalues requested')
      if (zpham) write(6,1070)
 1070 format(5x,'Printing of hamiltonian matrix requested')
      if (.not.zpham) write(6,1080)
 1080 format(5x,'Printing of hamiltonian matrix not requested')
      if (zprad) write(6,1090)
 1090 format(5x,'Printing of radial matrix elements requested')
      if (.not.zprad) write(6,1100)
 1100 format(5x,'Printing of radial matrix elements not requested')
      if (zpvec) write(6,1110)
 1110 format(5x,'Printing of eigenvectors requested'/)
      if (.not.zpvec) write(6,1120)
 1120 format(5x,'Printing of eigenvectors not requested'/)
      if (zvec) then
         write(6,1130) iout2
         write(6,1131) iout1
 1130 format(5x,'Eigenvalues & vectors   written to stream IOUT2 =',i4)
 1131 format(5x,'Restart information     written to stream IOUT1 =',i4)
         open(unit=iout1, form='unformatted')
         open(unit=iout2, form='unformatted')
      endif
      if (ztran) then
         write(6,1132) iwave
 1132 format(5x,'Wavefunction amplitudes written to stream IWAVE =',i4)
         write(6,1133) idip
 1133 format(5x,'Wavefunctions for dipole written to stream IDIP =',i4)
         write(6,1136) idipd
 1136 format(5x,'Wavefunctions for expect written to stream IDIPD =',i4)
         open(unit=iwave, form='unformatted')
      endif
      if (abs(jrot) .gt. 1) zpfun=.false.
      if (zpfun) then
         open(unit=ilev,form='formatted')
         if (jrot .eq. 0 .and. mod(ipar,2) .eq. 0) then
!           header on file ilev
            write(ilev,500) title
            write(6,1134) ilev
 1134 format(/5x,'Eigenvalues      written to start of stream ilev =',i4)
         else
!           position file ilev
  200       read(ilev,*,end=210,err=210)
            goto 200
  210       continue
! ******  inclusion of the following card is machine dependent ******
!           backspace ilev
            write(6,1135) ilev
 1135 format(/5x,'Eigenvalues      written at end   of stream ilev =',i4)
         endif
      endif
      if (idia .gt. 0) write(6,1140)
 1140 format(/5x,'Calculation performed in scattering coordinates')
      if (idia .le. 0 .and.  .not.zperp) write(6,1150)
 1150 format(/5x,'Calculation performed in Radau coordinates')
      if (idia .le. 0 .and. zperp) write(6,1151)
 1151 format(/5x,'Calculation performed in Radau coordinates with Z axis', &
              /5x,'perpendicular to the plane')
      if (zperp) write(6,1152)
1152  format('zlin not implemented with z-perp')

      if (ztwod) goto 886

      if (abs(idia) .eq. 2) then
         write(6,1180)
 1180    format(/5x,'Diatomic assumed homonuclear')
         if (ipar .eq. 1) then
            write(6,1190)
 1190       format(5x,'Odd parity functions in basis set')
         else if (ipar .eq. 0) then
            write(6,1200)
 1200       format(5x,'Even parity functions in basis set')
         else
            write(6,1205)
 1205       format(5x,'Illegal value of ipar for idia = +/-2: STOP')
            stop
         endif
         if (idia .eq. 2) then
            idvr=nalf2
            if (2*idvr .ne. nalf) goto 960
         endif
      else
         write(6,1210)
 1210    format(/5x,'Diatomic assumed hetronuclear')
         idvr=nalf
         ipar=0
      endif
      if (jrot .ne. 0) then
         jrot=abs(jrot)
         if (zrot) then
            if (kmin .ne. 0 .and. .not. zbisc) kmin=1
            write(6,1220)
 1220 format(/5x,'***  vibrational part of rot-vib calculation  ***')
            write(6,1260) jrot
            if (kmin .eq. 1) then
               write(6,1270)
 1270 format(12x,'with symmetric |Jk> + |J-k> functions in basis')
            else if (kmin .eq. 0) then
               write(6,1280)
 1280 format(12x,'with anti-symmetric |Jk> - |J-k> functions in basis')
            else
               kmin=2
               write(6,1275)
 1275 format(12x,'loop over symmetric & anti-symmetric |jk> functions')
            endif
            if (zladd) write(6,1240)
 1240 format(5x,'Number of angular grid points to be kept constant with k')
            if (.not. zladd) write(6,1250)
 1250 format(/5x,'Nunber of angular grid points to decrease with k')
         else
            write(6,1230) jrot,kmin
 1230 format(5x,'J =',i3,'  k =',i3,&
             /5x,'***  option to neglect coriolis interactions  ***')
            if (abs(kmin) .gt. abs(jrot)) then
               write(6,1235)
 1235 format(5x,'Error: k greater than J. STOP')
               stop
            endif 
         endif
         if (zbisc) then
            zembed=.false.
            write(6,1330)
 1330 format(/5x,'z axis embedded along the biscetor of r1 and r2')
            if (zlin) write(6,1340)
 1340 format(/5x,'Removal of functions with theta = 0 enforced')
         else if (.not.zperp) then
            if (zembed) write(6,1290) 2
 1290       format(/5x,'z axis embedded along the r',i1,' coordinate')
            if (.not. zembed) write(6,1290) 1
         else if (zperp) then 
             write(6,1291) 
 1291       format(/5x,'z axis embedded perp to the molecular plane')
         endif
      else
!        case j = 0
         write(6,1260) jrot
 1260    format(/5x,'J =',i3,' rotational state')
      endif

       write (6,1440)
 1440 format(/5x,'Routine DSYEV to do in core diagonalisation')
  886 continue

!     check input parameters are consistent
      npnt = max(npnt1,npnt2)
      nmax = max(nmax1,nmax2)
!     dimension of square matrices stored in triangular form :
      nlim1 = npnt1 * (npnt1+1) / 2
      nlim2 = npnt2 * (npnt2+1) / 2

!     declare dvr sizes for dimensioning the arrays
      ndima=npnta
      if (idia .le. -2 .and. .not. zrot) ndima=ndima-ipar
      ndimb=npntb
      ndimc=npntc

!     store parameters on disk files requested

      if (zvec) write (iout2) idia,ipar,npnta,npntb,npntc,max2d,max3d,neval

      if (zmors2 .and. .not. zquad2) goto 961
!      if (.not. ztheta .and. .not. zquad2) goto 962
!      if (idia .le. -2 .and. .not. zquad2) goto 963

      return
  960 write(6,970)
  970 format(//6x,'** nalf must be even when idia=2: stop **')
      stop
  961 write(6,972)
  972 format(//6x,'** can''t have zquad2 = f with zmors2 = t : stop **',&
              /6x,'               (not yet implemented)               ')
      stop
  962 write(6,973)
  973 format(//6x,'** can''t have zquad2 = f with ztheta = f : stop **',&
              /6x,'               (not yet implemented)               ')
      stop
  963 write(6,974)
  974 format(//6x,'** can''t have zquad2 = f with idia = -2: stop **',&
              /6x,'               (not yet implemented)               ')
      stop
      end

!##############################################################################
      subroutine ccmain 

!     subroutine ccmain is the 'real' main programme & contains    
!     the calls to the various subroutines which set & solve the
!     intermediate and the final hamiltonians.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common /mass/ xmass(3),g1,g2,xmassr(3)

      REAL*8, ALLOCATABLE, DIMENSION(:) :: dnorm1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r1m2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: r1m2t
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dnorm2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r2m2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dz1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dz2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: bass1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: bass2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: y1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: wt1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: y2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: wt2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: b
      REAL*8, ALLOCATABLE, DIMENSION(:) :: c
      REAL*8, ALLOCATABLE, DIMENSION(:) :: hbl1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: hbl2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xalf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: walf
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pleg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xlmatr
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xk1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xk2 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: r2m2t
      REAL*8, ALLOCATABLE, DIMENSION(:) :: jxcos
      REAL*8, ALLOCATABLE, DIMENSION(:) :: jwalf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sjwalf
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pjac
 
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x8/8.0d0/,xp5/0.5d0/,&
           toler/1.0d-8/

     ALLOCATE(dnorm1(0:nmax1),r1m2(nlim1),dnorm2(0:nmax2),&
               r2m2(nlim2),bass1(0:nmax1,npnt1),bass2(0:nmax2,npnt2),&
               y1(npnt1),r1(npnt1),wt1(npnt1),y2(npnt2),r2(npnt2),&
               wt2(npnt2),b(npnt+1),c(npnt+1),&
                dz1(npnt1,npnt1),dz2(npnt2,npnt2),&
               hbl1(nlim1),hbl2(nlim2),xalf(idvr),walf(idvr),&
               xlmatr(idvr,idvr),pleg (0:maxleg,idvr),xk1(npnt1,npnt1)&
               ,xk2(npnt2,npnt2),r2m2t(npnt2,npnt2),r1m2t(npnt1,npnt1))

    if (zperp) ALLOCATE(sjwalf(idvr),&
                 jxcos(idvr),jwalf(idvr),pjac(0:maxleg,idvr))

!     open the streams......(most not needed for symmetrised radau case)
      if (idia .gt. -2) then
         open(unit=ieigs1,form='unformatted')
         open(unit=ieigs2,form='unformatted')
         open(unit=ivecs1,form='unformatted')
         open(unit=ivecs2,form='unformatted')
         open(unit=ivint, form='unformatted')
         open(unit=iband, form='unformatted')
      endif
      open(unit=intvec,form='unformatted')

!     read in masses & basis set parameters
      call setcon(fixcos)
!.....and save them if necessary
      if (zvec) write(iout2) zembed,zmors1,zmors2,ztheta,zr2r1,xmass,&
                             g1,g2

!     set up binomial and normalisation arrays

      call setfac(dnorm1,dnorm2,cc1,cc2)

!     set up points, weights & basis for numerical integration

      if (ncoord .eq. 2) then
!     in less than 3-d cases fix the 3-d paramters
         if (zr2r1) then
            bass1(0,1) = x1
            wt1(1) = x1
            r1(1) = re1
!            call lagpt(2,y2,r2,wt2,b,c,cc2,bass2,dnorm2,npnt2,nmax2,&
!                       zmors2,re2,beta2,a2,iu2)

! Main change with respect to ALL previous versions: dz are the DVR<->FBR
! unitary transformation matrices for r, that is 
!                    dz_ij = L_j(x_i) sqrt(w_i)
! where L are normalized Laguerre pols, so there is no need for having the 
! weights and norms and whatever hanging around.
! But, in case of need, the weights are outputted from LAGPTNEW
!  wt = dsqrt(w). Those weights do not contain the exponential in themselves
! Also in case of need the polys are calculated and normalized in LAGPTNEW
! but not outputted. Finally, the norms as calculated in setfac can not be
! used with my definition of polys and weights, but are kept for they are used
! independently in keint_ routines.

     CALL LAGPTNEW(2,Y2,R2,WT2,DZ2,npnt2,nmax2,zmors2,RE2,beta2,A2,IU2)
            hbl2(1) = x0
!           setup kinetic energy & inertia integrals over r1
            if (zmors1) then
               fke = beta1 * beta1 / (x8 * ur1)
               call keints(hbl1,fke,nlim1,nmax1,iu1)
            else
               fke = beta1 / (x2 * ur1)
               call keint2(hbl1,fke,r1m2,dnorm1,nlim1,nmax1,npnt1,a1)
            endif
         else
            bass2(0,1) = x1
            wt2(1) = x1
            r2(1) = re2
!            call lagpt(1,y1,r1,wt1,b,c,cc1,bass1,dnorm1,npnt1,nmax1,&
!                       zmors1,re1,we1,a1,iu1)
      CALL LAGPTNEW(1,Y1,R1,WT1,DZ1,npnt1,nmax1,zmors1,RE1,WE1,A1,IU1)
           hbl1(1) = x0
!           setup kinetic energy & inertia integrals over r2
            if (zmors2) then
               fke = beta2 * beta2 / (x8 * ur2)
               call keints(hbl2,fke,nlim2,nmax2,iu2)
            else
               fke = beta2 / (x2 * ur2)
               call keint2(hbl2,fke,r2m2,dnorm2,nlim2,nmax2,npnt2,a2)
            endif
         endif
      else ! if ncoord eq 3 
!         call lagpt(1,y1,r1,wt1,b,c,cc1,bass1,dnorm1,npnt1,nmax1,zmors1,&
!                    re1,beta1,a1,iu1)
       CALL LAGPTNEW(1,Y1,R1,WT1,DZ1,npnt1,nmax1,zmors1,RE1,beta1,A1,IU1)
        if (idia .gt. -2) then
!           call lagpt(2,y2,r2,wt2,b,c,cc2,bass2,dnorm2,npnt2,nmax2,&
!                       zmors2,re2,beta2,a2,iu2)
      CALL LAGPTNEW(2,Y2,R2,WT2,DZ2,npnt2,nmax2,zmors2,RE2,beta2,A2,IU2)
!           setup kinetic energy & inertia integrals over r2
            if (zmors2) then
              fke = beta2 * beta2 / (x8 * ur2)
              call keints(hbl2,fke,nlim2,nmax2,iu2)
            else
              fke = beta2 / (x2 * ur2)
              call keint2(hbl2,fke,r2m2,dnorm2,nlim2,nmax2,npnt2,a2)
            endif
         endif
!        setup kinetic energy & inertia integrals over r1
         if (zmors1) then
            fke = beta1 * beta1 / (x8 * ur1)
            call keints(hbl1,fke,nlim1,nmax1,iu1)
         else
            fke = beta1 / (x2 * ur1)
            call keint2(hbl1,fke,r1m2,dnorm1,nlim1,nmax1,npnt1,a1)
         endif
      endif

!     write the quadrature points to disk for zvec = .true.
      if (zvec) then
         write (iout2) r1
         if (idia .gt. -2) write (iout2) r2
      endif
      if (ncoord .eq. 3) then
      endif

!     take square roots of the weights
!PB:  there is not need for this anymore.....
!      wt1=sqrt(wt1)
!      if (idia .gt. -2) wt2=sqrt(wt2)

!     set up the transformed kinetic energy integrals,  t'(hbl) t
!                                                       ~  ~~~  ~
!      call k1k2(xk1,hbl1,bass1,wt1,npnt1,nmax1,nlim1)
  CALL K1K2NEW(XK1,HBL1,DZ1,NPNT1,NMAX1,NLIM1)

  CALL K1K2NEW(R1M2T,R1M2,DZ1,NPNT1,NMAX1,NLIM1)

      if (idia .gt. -2)&
  CALL K1K2NEW(XK2,HBL2,DZ2,NPNT2,NMAX2,NLIM2)

!         call k1k2(xk2,hbl2,bass2,wt2,npnt2,nmax2,nlim2)

!     ...... and the inertia integrals for spherical oscillators
      if (.not.zmors2 .and. .not.ztwod .and. idia .gt. -2)&
  CALL K1K2NEW(r2m2t,r2m2,DZ2,NPNT2,NMAX2,NLIM2)
!          call k1k2(r2m2t,r2m2,bass2,wt2,npnt2,nmax2,nlim2)

!     some of the j>0 stuff to get the loop over k correct
      if (zrot) then
        kd = 1-min(kmin,1)
        ku = jrot
!       if looping over sym & anti-sym, do one extra k=1 block
        if (kmin .eq. 2) ku=ku+1
      else
        kd = kmin
        ku = kmin
      endif
      kkz12 = 0
      kkz0  = 0
!     for j > 0, store r**(-2) term for rotlev3
      if (ztran) then
         if (ku .gt. kd) then
            if (zembed .and. .not. zperp) then
              if (zquad2) then
                write(iwave) (xp5/(r2(i)*r2(i)*urr2),i=1,npnt2)
              else
                write(iwave) ((ur2*r2m2t(i,j)/urr2,i=1,npnt2),j=1,npnt2)
              endif
            else
              write(iwave) (xp5/(r1(i)*r1(i)*urr1),i=1,npnt1)
            endif
         else
            jdia=max(1,idia)
            jstart=kmin
            if (mod(jstart,jdia) .ne. ipar) jstart=jstart+1
            nang=(maxleg-jstart)/jdia+1
            mbass=idvr*npnt1*npnt2
            if (idia .eq. -2) mbass=idvr*max2d
            write(iwave) mbass,jstart,nang,mbass
         endif
         write(iwave) r1
         if (idia .gt. -2) write(iwave) r2
      endif

      if (idia .eq. -2) then
         max2d1=max2d
         max3d1=max3d
         DEALLOCATE(xk2,r2,r2m2t)
      endif

      DEALLOCATE(dnorm1,r1m2,dnorm2,r2m2,bass1,bass2,&
                 y1,wt1,y2,wt2,b,c,hbl1,hbl2,dz1,dz2)

      k_flag=0;
      iq=0

!     -------------  start rotational loop here  -------------
      do 40 kk=kd,ku
      if (kk .le. jrot) then
         kz=kk
      else
         kz=1
         kmin=0
         if(.not.zperp) ipar=mod(ipar+jrot,2)
      endif

!     first rewind the scratch files for a calculation with j>0
!     and, if needed, reposition iout2 after set up data.
60      if (kk .gt. kd) then
        if (idia .gt. -2) then
           rewind ieigs1
           rewind ieigs2
           rewind ivecs1
           rewind ivecs2
           rewind ivint
           rewind iband
        endif
        rewind intvec
        if (zvec) then
          rewind iout1
          rewind iout2
          do 45 ii=1,4
             read(iout2)
   45     continue
        endif
      endif

      realkz = dble(kz)
!     tswalf is the exact sum of weights for gauss-jacobi integration
      tswalf= x2**(kz+kz+1)/dble(kz+1)
      do 30 ia=1,kz
         tswalf=tswalf*dble(ia)/dble(kz+ia+1)
   30 continue

      if (zladd .or. kz .eq. kd) then
         nidvr = idvr
         nang  = nalf
         nang2 = nalf2
         lincr = kz
      else
         lincr = 0
         if (idia .ne. 2) then
            nidvr = idvr - kz
            nang  = nalf - kz
            nang2 = (nang+1)/2
         else if (ipar .eq. 0) then
            if(mod(kz,2).eq.1) kkz0 = kkz0 + 2
            nidvr = idvr - kkz0/2
            nang  = nalf - kkz0
            nang2 = (nang+1)/2
         else
            if(mod(kz,2).eq.0 .and. kz.gt.0) kkz12 = kkz12 + 2
            nidvr = idvr - kkz12/2
            nang  = nalf - kkz12
            nang2 = (nang+1)/2
         endif
      endif
      if (.not. zladd) then
         if (ztheta) then
            npnta=nidvr
         else
            npntc=nidvr
         endif
      endif

      if (ztwod) then
         xalf(1) = fixcos
         goto 333
      endif

      if(.not. zperp) then
        call jacobi(nang,nang2,xalf,walf,realkz,realkz,cswalf,tswalf)
        write(6,1000) nang,kz,(xalf(ii),walf(ii),ii=1,nang2)
 1000   format(//i8,' point Gauss-associated Legendre integration with',&
             ' k =',i3//5x,'integration points',11x,'weights',&
              //(f23.15,d25.12))
        write(6,1010) cswalf,tswalf
 1010   format(/4x,'Computed sum of weights',d22.15,&
              /4x,'Exact    sum of weights',d22.15//)
          if (abs((cswalf-tswalf)/tswalf) .gt. toler) then
             write(6,910)
  910       format(//5x,'Points & weights in error, adjust algorithm'//)
             stop
          endif
        call allpts(xalf,walf,nang,nang2)
        alf = x0
      else
        realj = DBLE(jrot)
        alf = dsqrt(0.5d0*(((realj)**2+realj)-realkz**2))
        bet = alf
        call jacbasis(xalf,walf,2*nang2,alf,bet)
        xalf=-xalf
        write(6,1020) nang,alf,(xalf(ii),walf(ii),ii=1,nang2)
 1020   format(//i8,' point Gauss-Jacobi integration with',&
             ' alpha = beta =',f7.3//5x,'integration points',11x,'weights',&
              //(f23.15,d25.12))
      endif

      if (zvec) write(iout2) xalf

      if (.not.zperp) then
!     set up Legendre polynomials for the transformation matrix
!     for z in plane case
        call asleg(pleg,maxleg,xalf,nidvr,kz,lincr)
      else
!     set up Jacobi polynomials for the transformation matrix
!     for z perpendicular case
        call jac_basis(nidvr,maxleg,alf,bet,xalf,pleg)
        do  i=1,nidvr
            walf(i) = sqrt(walf(i))
        enddo
      endif

!     save polinomials for rotlev3 or rotlev3b, or rotlev3z
      if (ztran) then
        write(iwave) xalf
        write(iwave) kz,maxleg,nidvr,lincr
        write(iwave) ((pleg(i,j)*walf(j),i=0,maxleg),j=1,nidvr)
        write(iwvpb) ((pleg(i,j)*walf(j),i=0,maxleg),j=1,nidvr)
      endif
      
!     build the transformed angular momentum matrix xlmatr;
      ipar0=0
      if (idia .eq. 2 .and. ipar .eq. 1) ipar0=1
      call lmatrx(xlmatr,pleg,walf,kz,ipar0,nidvr,lincr,alf)


  333 continue

      if (zperp) write(6,1152) kz, iq
 1152 format(/5x,'Calculate block for K =', i3,/5x,'q symmetry =', i3,/)

!     for ab2 molecules in radau coordinates, use separate main
!     driving routine
      if (idia .le. -2) then
         if (zrot) then
            if (kz .gt. kd .and. .not.zperp) ipar=mod(ipar+1,2)
            if (ipar .gt. 0) then
               max2d=max2d2
               max3d=max3d2
            else
               max2d=max2d1
               max3d=max3d1
            endif
         endif
         if (zperp) then
           call mkmain(xk1,r1m2t,xlmatr,r1,xalf,kz)
         else
           call nfmain(xk1,r1m2t,xlmatr,r1,xalf,kz)
         endif    
      else
         call jhmain(xk1,xk2,xlmatr,r1,r2,xalf,r2m2t,kz)
      endif
      
!     calculate block k=1, iq=1 for z-perpendicular case
      if (zperp .and. kk .eq. 1 .and. k_flag .eq. 0) then
        kz=kk
        iq=1
        k_flag=1
        go to 60
      endif
      iq=0
      
   40 continue

      DEALLOCATE(xalf,walf,xlmatr,pleg,xk1,r1)
      if (zperp) DEALLOCATE(jxcos,jwalf,sjwalf,pjac)
      if (idia .gt. -2) DEALLOCATE(r2,xk2,r2m2t)
      return
      end

!#############################################################################
      subroutine setcon(fixcos)

!     read in masses & set constants for radial basis sets          #007

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2
!     save masses & g in case they are needed in the potential routine
      common /mass/ xmass(3),g1,g2,xmassr(3)
!     amtoau converts amu (proton masses) to au (electron masses).
      data amtoau/1.8228883d03/
      data x0,xp5,x1,x4/0.0d0,0.5d0,1.0d0,4.0d0/

!     read cos(theta) for fixed angle 2-d calculation

      read(5,5) fixcos
      if (ztwod) write(6,1088) fixcos
 1088 format(//5x,'two-d fixed angle vibrational problem with'&
             //5x,'***** fixed value of cos(theta) =',f6.2,' *****'/)

!     read masses of the atoms in atomic mass units
!     first vibrational mass...
      read(5,5)     xmass
    5 format(3f20.0)
!     .... then rotational mass
      read(5,5)     xmassr
!     Default rotational mass to vibration mass if it is not set
      if (xmassr(1).le.x0) xmassr = xmass

!     read cut off energies
!     read parameters defining energy cut offs for each block
      read(5,5) emax1,emax2
      if (.not. zall) then
         if (zcut) then
            if (idia .gt. -2) write(6,990) emax1,emax2
  990       format(//5x,'Cut-off energies in wavenumbers:',2d16.8/)
            if (idia .eq. -2) write(6,991)       emax2
  991       format(//5x,'Final cut-off energy in wavenumbers:',2d16.8/)
         else
            if (idia .gt. -2) write(6,992) emax1
  992       format(//5x,'First cut-off energy in wavenumbers:',1d16.8/)
         endif
      endif
!     set default value of g1 and g2
      if (idia .ge. 1) then
!        scattering coordinates
         g1 = xmass(2) / (xmass(2) + xmass(3))
         g2 = x0
      else
!        radau coordinates
         a = sqrt(xmass(3) / (xmass(1)+xmass(2)+xmass(3)))
         b = xmass(2) / (xmass(1)+xmass(2))
         g1 = x1 - a / (a+b-a*b)
         g2 = x1 - a / (x1-b+a*b)
      endif
!     ncoord = 3: read parameters for r1 radial basis (see below)
!     ncoord = 2: read fixed r1 bondlength re1. diss1 & we1 dummy
      read(5,5)     re1,diss1,we1
!     read parameters for r2 radial basis function,
!     for morse oscillator functions use the following:
!     re2: equilibrium bondlength of r2 coordinate (in bohr)
!     diss2: dissociation energy of the r2 coordinate (in hartree)
!     we2: fundamental stretching vibration of r2 (in hartree)
!     for spherical oscillator functions use the following:
!     re2 : dummy
!     diss2: order of laguerre polynomials used (dimensionless)
!     we2: fundamental stretching vibration of r2 (in hartree)
!     all are treated as variationally optimisable parameters.
      read(5,5)     re2,diss2,we2
      write(6,1000) xmass
 1000 FORMAT(/5X,'Vibrational nuclear mass in AMU:',3F12.6)
      if (jrot.ne.0) write(6,1001) xmassr
 1001 FORMAT( 5X,'Rotational  nuclear mass in AMU:',3F12.6/)
!     compute the effective moments of inertia
      ur1 = amtoau/(g2*g2/xmass(1)+x1/xmass(2)+(x1-g2)**2/xmass(3))
      ur2 = amtoau/(x1/xmass(1)+g1*g1/xmass(2)+(x1-g1)**2/xmass(3))
      urr1 = amtoau/(g2*g2/xmassr(1)+x1/xmassr(2)+(x1-g2)**2/xmassr(3))
      urr2 = amtoau/(x1/xmassr(1)+g1*g1/xmassr(2)+(x1-g1)**2/xmassr(3))
 
      if (ncoord .eq. 3) goto 20
      if (zr2r1) then
         write(6,1010) re1,ur1
 1010 format(/5x,'r1 fixed bondlength =',f8.4,' bohr',&
                  ' & reduced mass =',d16.7,' a.u.'/)
      else
         write(6,1011) re2,ur2
 1011 format(/5x,'r2 fixed bondlength =',f8.4,' bohr',&
                  ' & reduced mass =',d16.7,' a.u.'/)
      endif
      if (zr2r1) goto 30
   20 continue
      if (zmors1) then
          write(6,1020) 1,re1,diss1,we1
 1020 format(/5x,'Morse function parameters for r',i1,' basis',&
             /5x,'r equilibrium =',f8.4,' bohr, dissociation energy',&
          d15.7,' hartree &  vibrational frequency =',d15.7,' hartree')
         beta1 = we1 * sqrt(xp5*ur1/diss1)
         a1 = x4 * diss1 / we1
         iu1 = int(a1+xp5)
         write(6,1030) ur1,beta1,a1,iu1
 1030 format(/5x,'Constants used to construct morse oscillators:',&
             /5x,'reduced mass =',d16.7,' a.u., beta =',f8.4,&
                  ' (1/bohr), a =',d16.7,' and u =',i5)
      else
          a1=diss1
          beta1 = sqrt(we1 * ur1)
          write(6,1039) a1,we1,ur1,beta1
 1039 format(/5x,'Spherical oscillator parameters for r1 basis:',&
             /5x,'alpha =',f10.5,&
                 ' &  vibrational frequency =',d15.7,' hartree',&
            //5x,'Constants used to construct spherical oscillators:',&
             /5x,'reduced mass =',d16.7,' a.u., beta =',f12.6,&
                  ' bohr**-2')
      endif
      if (ncoord .eq. 2 .and. .not. zr2r1) goto 40
      if (idia .eq. -2) goto 40
   30 continue
      if (zmors2) then
         write(6,1020) 2,re2,diss2,we2
         beta2 = we2 * sqrt(xp5*ur2/diss2)
         a2 = x4 * diss2 / we2
         iu2 = int(a2+xp5)
         write(6,1030) ur2,beta2,a2,iu2
      else
         a2=diss2
         beta2 = sqrt(we2 * ur2)
         write(6,1040) a2,we2,ur2,beta2
 1040 format(/5x,'Spherical oscillator parameters for r2 basis:',&
             /5x,'alpha =',f10.5,&
                 ' &  vibrational frequency =',d15.7,' hartree',&
            //5x,'Constants used to construct spherical oscillators:',&
             /5x,'reduced mass =',d16.7,' a.u., beta =',f12.6,&
                  ' bohr**-2')
      endif
   40 continue
      if (ztran) then
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
            write(iwave) idia,ipar,idvr,npnt1,npnt2,jrot,kmin,neval,nlim
            write(iwave) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
            write(iwave) re1,diss1,we1,re2,diss2,we2
         else
            write(iwave) idia,ipar,idvr,npnt1,npnt1,jrot,kmin,neval,nlim
            write(iwave) zembed,zmors1,zmors1,xmass,g1,g2,zncor,zquad2
            write(iwave) re1,diss1,we1,re1,diss1,we1
         endif
      endif
      return
      end

!#########################################################################
      subroutine setfac(dnorm1,dnorm2,cc1,cc2)

!     setfac initialises binomial array:                            #021
!       binom(i+1,j+1) = i! / (j! * (i-j)!)
!     and pseudo-normalisation array:
!       dnorm(m) = sqrt((m-1)! * binom(npnt+iu,npnt-m))

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2

      real*8, dimension(0:nmax1) :: dnorm1
      real*8, dimension(0:nmax2) :: dnorm2
      real*8, dimension(npnt+1) :: fact
    
      data xp5/0.5d0/,x1/1.0d0/
      fact(1) = x1
      count = x1
      do 10 i=1,npnt
      fact(i+1) = count * fact(i)
      count = count + x1
   10 continue
      if (ncoord .eq. 2) then
        if (.not. zr2r1) then
          if (zmors1) then
            alf=dble(iu1)
          else
            alf=a1+xp5
          endif
          call norms(dnorm1,fact,cc1,alf,npnt1,nmax1)
        else
          if (zmors2) then
            alf=dble(iu2)
          else
            alf=a2+xp5
          endif
          call norms(dnorm2,fact,cc2,alf,npnt2,nmax2)
        endif
      else
        if (zmors1) then
          alf=dble(iu1)
        else
          alf=a1+xp5
        endif
        call norms(dnorm1,fact,cc1,alf,npnt1,nmax1)
        if (idia .gt. -2) then
           if (zmors2) then
             alf=dble(iu2)
           else
             alf=a2+xp5
           endif
           call norms(dnorm2,fact,cc2,alf,npnt2,nmax2)
        endif
      endif

      return
      end

!##############################################################################
      subroutine norms(dnorm,fact,cc,alf,npnt,nmax)
!     set up factors for normalising the radial basis functions     #022

      implicit real*8 (a-h,o-y)

      real*8, dimension(0:nmax) :: dnorm
      real*8, dimension(npnt+1) :: fact
      real*8, dimension(npnt+1) :: bin

      data x1/1.0d0/
      count = dble(npnt) + alf
!     cc is exact sum of weights for npnt Gauss-Laguerre integration
      cc = fact(npnt) / count
!     normalisation array for l(i,alf): first set up binomials
! bin(i)=g(N+1+a)/(g(N+1-i)*g(a+i+1))
      npt1 = npnt + 1
      bin(npt1) = x1
      do 20 i=1,npnt
      n = npt1 - i
      bin(n) = bin(n+1) * count / dble(i)
      count = count - x1
   20 continue
      npt1 = npt1 + 1
      do 30 i=0,nmax
      dnorm(i) = sqrt(bin(i+1)*fact(i+1)*fact(npt1-i-1))
    30 continue
      return
      end

!############################################################################
      subroutine lagpt(ir,y,r,wt,b,c,cc,bass,dnorm,npnt,nmax,zmorse,&
                       re,beta,a,iu)

!     subroutine lagpt gets integration points and weights for      #015
!     npnt gauss laguerre integration and sets up basis
!     functions at the integration points.

      implicit real*8(a-h,o-y), logical (z)
     
      real*8, dimension(npnt+1) :: b
      real*8, dimension(npnt+1) :: c
      real*8, dimension(npnt) :: y
      real*8, dimension(npnt) :: r
      real*8, dimension(npnt) :: wt
      real*8, dimension(0:nmax) :: dnorm
      real*8, dimension(0:nmax,npnt) :: bass
         
      data x0,xp5,x1,x2/0.0d0,0.5d0,1.0d0,2.0d0/,toler/1.0d-8/

      if (zmorse) alf=dble(iu)
      if (.not. zmorse) alf = a + xp5
      alfm1=alf-x1
!     set up integration points and weights.
      call laguer(npnt,y,wt,alf,b,c,csx,csa,tsx,cc)
      tsa = x1 / (dnorm(0) * dnorm(0))
      write(6,1000) npnt,ir
 1000 format(/,i8,' point Gauss-Laguerre integration',&
             /,5x,'integration points',11x,'weights',9x,&
                  'corresponding r',i1,/)
      do 60 i=1,npnt
      if (zmorse) then
!         calculate potential at r = re+beta(**-1)*ln(a/y)
          r(i) = re + dlog(a/y(i)) / beta
      else
!         calculate potential at r = sqrt(y/beta)
          r(i) = sqrt(y(i)/beta)
      endif
      write(6,1010) y(i),wt(i),r(i)
 1010 format (f23.15,d25.12,f13.5)
      if (r(i) .lt. x0) write(6,1015) i
 1015 format(5x,'***** warning: for integration point',i3,&
             ', r less than zero *****')

!     calculate unnormalised laguerre polynomials at y

!     polynomial of order 0
      bass(0,i) = x1
      if (nmax .lt. 1) goto 70
!     polynomial of order 1
      amx = alf + x1 - y(i)
      bass(1,i) = amx
!     use recurrence relationships for polynomials of order > 2
!     n * l(n,alf) = (2*n+alf-1-x)*l(n-1,alf) - (n+alf-1)*l(n-2,alf)
      en = x1
      do 80 n=2,nmax
      en = en + x1
      amx = amx + x2
      bass(n,i) = (amx * bass(n-1,i) - (alfm1+en) * bass(n-2,i)) / en
   80 continue
   70 continue

      do 90 n2=0,nmax
!     normalise polynomials
      bass(n2,i) = bass(n2,i) * dnorm(n2)
   90 continue
   60 continue

!     check that the correct points & weights have been generated
      write(6,1020) csx,csa,tsx,tsa
 1020 format(/4x,'Computed sum of points',d22.15,' & weights',d22.15,&
             /4x,'Exact    sum of points',d22.15,' & weights',d22.15)
      if (abs((csx-tsx)/tsx) .gt. toler) goto 900
      if (abs((csa-tsa)/tsa) .gt. toler) goto 900
      return
  900 write(6,910)
  910 format(//5x,'points & weights in error, adjust algorithm',//)
      stop
      end

!##########################################################################
      subroutine laguer(nn,x,a,alf,b,c,csx,csa,tsx,cc)

!     calculates points & weights for gauss-laguerre integration    #016
!     see:
!     "gaussian quadrature formulas" by a.h.stroud & d.secrest
!      1966, prentice-hall, p.32.
!     **** version to avoid overflows (j.t. 25/11/81) ****
!     calculates weights divided by gamma(nn+alf+1)
!     this is an initialsation entry

      implicit real*8 (a-h,o-y)
      real*8, dimension(nn) :: x
      real*8, dimension(nn) :: a
      real*8, dimension(nn+1) :: b
      real*8, dimension(nn+1) :: c      

      data eps/1.0d-12/,x1/1.0d0/
      csx=0.0d0
      csa=0.0d0
      fa=alf+1.0d0
!     cc = n!                      denominator for pseudo-weights: a
!     b(n) = (alf + 2n -1)             b & c for recurrence relation
!     c(n) = (n - 1) * ( alf + n - 1)
      b(1)=fa
      c(1)=0.0d0
      fn=1.0d0
      do 1 j=2,nn
         fa=fa+2.0d0
         b(j)=fa
         c(j)=fn*(alf+fn)
         fn=fn+1.0d0
    1 continue
      tsx=fn*(alf+fn)
      xt1=0.0d0
!     formulas for initial point & step chosen because they work!
      xt=(1.0d0+alf)*(2.0d0+alf)/(1.0d0+3.0d0*fn+2.0d0*alf)
      step=3.0d0*(1.0d0+alf)/(1.0d0+3.0d0*fn+alf)
      call lgrecr(pt,dpn,pn1,xt,nn,alf,b,c)

      do 7 i=1,nn
      if (i .gt. 2) goto 4
!     smallest two zeros: found by "brute force" search
    2 xt2 = xt + step
      call lgrecr(pt2,dpn,pn1,xt2,nn,alf,b,c)
      if (dsign(x1,pt)*dsign(x1,pt2) .gt. 0.0d0) goto 5
      pt = pt2
      xt = 0.5d0 * (xt + xt2)
      go to 6
    5 pt = pt2
      xt = xt2
      go to 2
!     all other zeros: found using formula of stroud & secrest
    4 fi = dble(i-2)
      r1 = (1.0d0+2.55d0*fi)/(1.9d0*fi)
      r2 = 1.26d0*fi*alf/(1.0d0+3.5d0*fi)
      ratio = (r1+r2)/(1.0d0+0.3d0*alf)
      xt = xt + ratio*(xt-xt2)

    6 call lgroot(xt,nn,alf,dpn,pn1,b,c,eps)
      xt2=xt1
      xt1=xt
      x(i) = xt
      a(i) = cc/dpn/pn1
      csx = csx + xt
      csa = csa + a(i)
    7 continue
      return
      end

!##############################################################################
      subroutine lgroot(x,nn,alf,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains          #017
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x
!     this routine is due to stroud & secrest (see subroutine laguer)

      implicit real*8 (a-h,o-y)

      real*8, dimension(nn+1) :: b
      real*8, dimension(nn+1) :: c     

      data itmax/10/
      iter=0
    1 iter=iter+1
      call lgrecr(p,dpn,pn1,x,nn,alf,b,c)
      d = p/dpn
      x = x-d
      if (abs(d/x) .le. eps) return
      if (iter .lt. itmax) goto 1 
      write(6,100) iter,d,x
  100 format(5x,'warning: noconvergence after',i4,' iterations',&
             /,5x,'current difference',d26.15,' & root',d26.15)
      return
      end

!###########################################################################
      subroutine lgrecr(pn,dpn,pn1,x,nn,alf,b,c)

!     uses recurrence relations to set up polynomials               #018
!     this routine is due to stroud & secrest (see subroutine laguer)

      implicit real*8 (a-h,o-y)

      real*8, dimension(nn+1) :: b
      real*8, dimension(nn+1) :: c     

      p1 = 1.0d0
      p = x - alf - 1.0d0
      dp1 = 0.0d0
      dp = 1.0d0
      do 1 j=2,nn
         q  = (x-b(j))* p-c(j)* p1
         dq = (x-b(j))*dp-c(j)*dp1 + p
         p1 = p
         p  = q
         dp1= dp
         dp = dq
    1 continue
      pn = p
      dpn= dp
      pn1= p1
      return
      end

!###########################################################################
      subroutine keints(hbl,fke,nlim,nmax,iu)

!     keints calculates analytic kinetic energy integrals over r    #012
!     for morse oscillator-like functions

      implicit real*8 (a-h,o-y), logical (z)

      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nlim) ::  hbl

      data x0/0.0d0/
      index = 0
      do 10 n2=1,nmax+1
        do 20 n1=1,n2
          index=index+1
          if (n1 .eq. n2) then
!             special case:  n1 = n2
              hbl(index) = fke * dble(2*(n2-1)*(n2+iu)+iu+1)
          else if (n1+2 .eq. n2) then
!                  special case:  n2 = n1 + 2
        hbl(index) = - fke * sqrt(dble((iu+n1+1)*(iu+n1))*dble((n1+1)*n1))
          else
!           n1+1 = n2 or n2 > n1 + 2  all matrix elements are zero
            hbl(index) = x0
          endif
   20   continue
   10 continue
      if (.not. zprad) return
!     write kinetic energy integrals if requested
      write(6,500)
  500 format(//,5x,'radial kinetic energy matrix calculated',&
              ' analytically',/)
      call symout(hbl,nmax+1)
      return
      end

!#######################################################################
      subroutine keint2(hbl,fke,rm2,dnorm,nlim,nmax,npnt,alf)

!     keint2 calculates analytic kinetic energy integrals over r2   #013
!     and moment of intertia integral for spherical oscillator functions

      implicit real*8 (a-h,o-y), logical (z)

      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nlim) :: hbl
      real*8, dimension(nlim) :: rm2
      real*8, dimension(*) :: dnorm

      data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/
      gam = fke / (alf + xp5)
      do 10 n1=1,npnt
         gam = gam /(dble(n1)+alf+xp5)
   10 continue  
! gam = g(alf+0.5) / g(a+0.5+N+1)
      fact = x1
      fn = x0
      sum = gam
      do 20 n1 = 1,nmax+1
      index = (n1 * (n1+1)) / 2
      do 30 n2 = n1,nmax+1
      rm2(index) = dnorm(n1) * dnorm(n2) * sum
      index = index + n2
   30 continue
      gam = (fn+alf+xp5) * gam
      fn = fn + x1
      fact = fn * fact
      sum = sum + gam / fact
   20 continue
      fact = alf * (alf + x1)
      do 40 index=1,nlim
      hbl(index) = - fact * rm2(index)
   40 continue
      index = 0
      fn = - x1
      do 50 n2=1,nmax+1
      fn = fn + x1
      index=index+n2
!     special case:  n1 = n2
      hbl(index) = hbl(index) + fke * (fn+fn+alf+1.5d0)
!     special case:  n2 = n1 + 1
      if (n2 .gt. 1)&
          hbl(index-1) = hbl(index-1) + fke * sqrt(fn*(fn+alf+xp5))
   50 continue
      if (.not. zprad) return
!     write kinetic energy & inertia integrals if requested
      write(6,500)
  500 format(//5x,'radial kinetic energy matrix calculated',&
                   ' analytically'/)
      call symout(hbl,nmax+1)
      write(6,510)
  510 format(//5x,'moment of inertia matrix calculated analytically'/)
      call symout(rm2,nmax+1)
      return
      end

!##########################################################################
      subroutine k1k2(xk,hbl,bass,wt,npnt,nmax,nlim)

!     set up the transformed kinetic energy integrals,  t'(hbl) t
!                                                       ~  ~~~  ~
!     (note that the radial basis functions are already normalised)

      implicit real*8 (a-h,o-y), logical (z)

      real*8, dimension(npnt,npnt) :: xk
      real*8, dimension(nlim) :: hbl
      real*8, dimension(0:nmax,npnt) :: bass
      real*8, dimension(npnt) :: wt

      xk = 0.0d0

      do 10 k=1,npnt
        do 20 kp=1,k
          wtkkp = wt(k)*wt(kp)
          do 30 m=0,nmax
            t = bass(m,k) * wtkkp
            do 40 mp=0,nmax
              in = max(m,mp) * (max(m,mp)+1)/2 + min(m,mp) + 1 
              xk(k,kp) = xk(k,kp) + (hbl(in) * t * bass(mp,kp))
   40       continue
   30     continue
        xk(kp,k)=xk(k,kp)
   20   continue
   10 continue
      return
      end

!############################################################################
      subroutine jacobi(nn,nn2,x,a,alf,bta,csa,tsa)

!     calculates zeros x(i) of the nn'th order jacobi polynomial
!     pn(alf,bta) for the segment (-1,1) & and corresponding weights
!     for gauss-jacobi integration. this routine uses a brute force
!     search for zeros due to GJ Harris (2001).
!     note that for our purposes, alf= bta= nu.


      implicit real*8(a-h,o-z)
      real*8, dimension(nn) :: x,a,b,c,xt
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/,&
           eps/1.0d-12/,xstep/1.0d-6/
      fn= dble(nn)
      csa= x0
      c(1) = x0
      b(1) = x0
      do 10 i=2,nn
      xi= float(i)
      b(i) = x0
      c(i)= x4*(xi-x1)*(alf+xi-x1)*(bta+xi-x1)*(alf+bta+xi-x1)/&
             ((alf+bta+x2*xi-x1)*(alf+bta+x2*xi-x2)*&
              (alf+bta+x2*xi-x2)*(alf+bta+x2*xi-x3))
  10 continue
      cc=tsa
      do 15 i=2,nn
       cc= cc*c(i)
   15 continue

! step through the jacobi polynomial and 
! notes a first guess for the posititions of the zeros in the 
! array xt(ii). Zeros are found by looking for a change in sign.

      ii=0
      pm1=x1
      xxx=x1
 30   continue
      call recur(p,dp,pn1,xxx,nn,alf,bta,b,c)

      if (pm1*p .lt. x0) then
         pm1 = -pm1
         ii = ii +1
         xt(ii)=xxx-0.5*xstep
      endif

      if (ii .eq. nn2) then
         do 40 i=1,nn2
         call recur(ptemp,dp,pn1,xt(i),nn,alf,bta,b,c)
40      continue
      else
         xxx=xxx-xstep
         if (xxx .gt. -1.5*xstep) goto 30
         write(6,*) "Incorrect number",ii-1," of zeros found in JACOBI"
         stop
      endif 

       do 20 i=1,nn2
       call root(xt(i),nn,alf,bta,dpn,pn1,b,c,eps)
       x(i)= xt(i)
       a(i)= cc/(dpn*pn1)
       csa= csa + a(i) + a(i)
  20  continue
      if (2*nn2 .ne. nn) csa=csa-a(nn2)
      return
      end

!#####################################################################
      subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x.

      implicit real*8(a-h,o-z)
      real*8, dimension(nn) :: b
      real*8, dimension(nn) :: c
   
      iter= 0
1     iter= iter + 1
      call recur(p,dp,pn1,x,nn,alf,bta,b,c)
      d = p/dp
      x = x - d
      if(abs(d) .le. eps) goto 3 
      if(iter .lt. 10) goto 1 
3     dpn= dp
      return
      end
      subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
      implicit real*8(a-h,o-z)
      real*8, dimension(nn) :: b
      real*8, dimension(nn) :: c
    
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
      p1= x1
      p= x + (alf-bta)/(alf + bta + x2)
      dp1= x0
      dp= x1
      do 1 j=2,nn
        q= (x - b(j))*p - c(j)*p1
        dq= (x - b(j))*dp + p - c(j)*dp1
        p1= p
        p= q
        dp1= dp
        dp= dq
1     continue
      pn= p
      dpn= dp
      pn1= p1
      return
      end

!################################################################
      subroutine allpts(xalf,walf,nang,nang2)

!     takes the points & weights generated by legpt for the half-range
!     and creates new arrays for the full-range (-1,+1).

      implicit real*8 (a-h,o-y)

      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth

      real*8, dimension(idvr) :: xalf,walf

      scale = dble(max(1,idia))
      do 10 i=1,nang2
      walf(i) = sqrt(scale*walf(i))
      if (idia .eq. 2) goto 10
      xalf(nang+1-i) = xalf(i)
      xalf(i)        =-xalf(i)
      walf(nang+1-i) = walf(i)
   10 continue
  
      return
      end

!########################################################################
      subroutine asleg(pleg,lmax,x,nn2,kz,lincr)

!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, page 182,
!     for the polynomial part of associated legendre functions.
!     we have removed sin(theta)**2 for nu = 1.
!     this enables us to use jacobi integration with alf = bta = nu,
!     using routines derived from beidenharn and louck.

      implicit real*8 (a-h,o-z)
      real*8, dimension(0:lmax,nn2) :: pleg
      real*8, dimension(nn2) :: x
      real*8, dimension(0:lmax) :: pnorm

      data x1/1.0d0/,x2/2.0d0/
      m = kz
      if (m.lt.0) goto 999
      do 10 i=1,nn2
      
      !For high J and high k the value pmm can become too large over this loop, i.e. larger than HUGE(X)
      !Therefore in these instances we need pmm to be very small to begin with
      !We later divide by this factor to achieve the originally required number
      !For so2 at J = 200 and k > 90 we find this necessary, so we add this conditional
	if(m.LT.90)then
	  pmm = x1
	else
	  pmm = 1d-250
	end if

      fact = x1
      do 11 j=1,m
          pmm = -pmm * fact
          fact = fact + x2
   11 continue
      
      !Even when we use the above trick to reduce the size of pmm, it can still be too large for subsequent loops
      !We therefore conditionally divide pmm by a large enough factor - we take this to be for k > 150 (for so2)
      if(m.GT.150)then
	pmm = pmm/1.0d200
      end if
      
      pleg(0,i) = pmm
      pmmp1= x(i)*(m+m+1)*pmm
      pleg(1,i)= pmmp1
      
      
      ll=1
      do 2 l= 2+m,lmax+lincr
      r2lm1 = dble(l+l-1)
      rlpmm1= dble(l+m-1)
      rlmm  = dble(l-m)
      pll= (x(i)*r2lm1*pmmp1 - rlpmm1*pmm)/rlmm

      
      pmm= pmmp1
      pmmp1=pll
      ll=ll+1
      pleg(ll,i)= pll
2     continue
10    continue
      
!     set up the normalisation constants
!     (pnorm)**2 = (2j + 1)/2   *   (j - k)! / (j + k)!
      jstart = m
      jj = -1
      do 13 j = jstart,lmax+lincr
      fact = x1
      do 12 i = j-m+1,j+m

	 !After a certain point, square-rooting this factor doesn't reduce the number sufficiently to a coping level
	 !Therefore we take the fourth root instead. We set the cut-off point to be the same as that for pmm, k > 90
	 if(m.LT.90)then
	   facti = sqrt(dble(i))
	 else
	   facti = dble(i)**(1.0d0/4.0d0)
	 end if
        
	 fact = fact * facti

      12 continue
      rj = dble(j)
      jj = jj + 1
      
      if(m.LT.90)then
	   pnorm(jj) = sqrt((rj + rj + x1) /2)/fact
	 else
	   pnorm(jj) = (((rj + rj + x1) /2)**(1.0d0/4.0d0))/fact
      end if
      
      
   13 continue
!     now normalise the polynomials
      do 14 i=1,nn2
         jj = -1
         do 15 j=jstart,lmax+lincr
            jj = jj + 1

            pleg(jj,i) = pleg(jj,i) * pnorm(jj)
	    
	    !This is where we multiply by the large factor initially divided from pmm at the beginning of its loop
	    if(kz.GE.90)then
	      pleg(jj,i) = pleg(jj,i) * 1.0d250 * pnorm(jj)
	    end if
	    
	    !And for the k's which are greater than 150, we need to divide by the reducing factor again
	     if(kz.GT.150)then
	      pleg(jj,i) = pleg(jj,i) * 1.0d200
	    end if

   15    continue
   14 continue


      return
999   write(6,200)
200   format(/,/,5x,'improper argument in subroutine asleg',/)
      stop
      end

!###################################################################
      subroutine lmatrx(xlmatr,pleg,walf,kz,ipar0,nidvr,lincr,alf)

!     this subroutine sets up the lower triangle of the transformed
!     angular momentum matrix l(alpha,alpha')

      implicit real*8 (a-h,o-y), logical (z) 
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(idvr,idvr) :: xlmatr
      real*8, dimension(0:maxleg,idvr) :: pleg
      real*8, dimension(idvr) :: walf

      data x0/0.0d0/

      if(zperp) then      
        jstart = 0
        iparam = jstart
      else
        jstart = kz
        iparam = lincr
      endif

      realkz = kz

      jdia=max(1,idia)
      jj0=-jdia
      if (.not.zperp .and. mod(jstart,jdia) .ne. ipar0) then
          jj0=jj0+1
          jstart=jstart+1
      endif
      do 10 k= 1,nidvr
        term =  walf(k)
        do 11 kp=k,nidvr
          sumj1=x0
          jj = jj0
          do 20 j=jstart,maxleg+iparam,jdia
            jj = jj + jdia
            if (zperp) then
               sumj1 = sumj1 + pleg(jj,k) * pleg(jj,kp) * &
                       dble((j+1+alf)*(j+alf))
            else
               sumj1 = sumj1 + pleg(jj,k) * pleg(jj,kp) * &
                       dble((j+1)*(j))
            endif
   20     continue
        xlmatr(kp,k) = sumj1 * term * walf(kp)
        xlmatr(k,kp) = xlmatr(kp,k)
   11   continue
   10 continue   
      if (.not. zlmat) return
!     write xlmatr if requested
      write(6,1010) kz,ipar0
 1010 format(5x,'L-matrix for kz =',i3,', ipar =',i2/)
      call sqout(xlmatr,nidvr)
      return
      end

!##########################################################################

      subroutine jhmain(xk1,xk2,xlmatr,r1,r2,xalf,r2m2t,kz)

!     this routine controls the dvr calculation in all cases except
!     symmetrised radau coordinates.
!     written by james henderson

      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(npnt1) :: r1
      real*8, dimension(npnt2) :: r2
      real*8, dimension(idvr) :: xalf
      real*8, dimension(idvr,idvr) :: xlmatr
      real*8, dimension(npnt1,npnt1) :: xk1
      real*8, dimension(npnt2,npnt2) :: xk2
      real*8, dimension(npnt2,npnt2) :: r2m2t
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: eigs2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ham1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: hband
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eigs1d
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::vecs1d
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham3
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::cint
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::cintp
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eigs2d
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eval
      REAL*8, ALLOCATABLE, DIMENSION(:) ::evall
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: vecs1l
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::vecs2l
      REAL*8, ALLOCATABLE, DIMENSION(:) ::vecs3l
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::phi
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::psi 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iv1
      INTEGER, ALLOCATABLE, DIMENSION(:) ::ndim2d
      INTEGER, ALLOCATABLE, DIMENSION(:) ::iv2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv1l
      INTEGER, ALLOCATABLE, DIMENSION(:) ::iv2l
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ndim2l

      ALLOCATE(ham1(ndima,ndima),eigs1d(max2d),&
               eig1(ndima),iv1(ndimc,ndimb),vecs1d(max2d,ndima))

      do l1=1,idvr
         do l2=1,idvr
            write(39,*)l1,l2,xlmatr(l1,l2)
         end do
      end do


!      xlmatr=1.d0
      do l=1,idvr
!         xlmatr(l,l)=1.d0
      end do


      term  = dble(jrot * jrot + jrot - (2 * kz * kz))
!     construct the one-dimensional hamiltonian matrix h1(npnta,npnta')
!     for each npntb and npntc, then solve by diagonalisation.
      icall = 0
      nsum = 0
      do 10 ione = 1,npntc
      nham2 = 0
      do 20 itwo = 1,npntb

      if (zr2r1) then
         call mkham1(ham1,xlmatr,ione,itwo,term,r1,r2,xalf,xk1,xk2,kz)
      else
         call mkham1(ham1,xlmatr,itwo,ione,term,r1,r2,xalf,xk1,xk2,kz)
      endif

!     diagonalise each block, saving the eigenvalues and vectors that
!     fall below the cut-off energy emax1.
      call diag_dac(ham1,ndima,npnta,eig1)
      call cut1d(ham1,eig1,iv1(ione,itwo),eigs1d,vecs1d,nham2,icall)
   20 continue
!      write(6,985) nham2
  985 format(5x,' nham2 = ',i4)
      nsum = nsum + nham2
      if (ione .eq. npntc) write(6,986) nsum
  986 format(/5x,' sum = ',i5)
      if (nham2 .gt. 0) then
!        dump the 1d eigenavlues & vectors to disk for each ione
         call outrow(eigs1d,nham2,ieigs1)
         do 2 ka=1,npnta
         if (zvec) call outrow(vecs1d(1,ka),nham2,iout2)
         call outrow(vecs1d(1,ka),nham2,ivecs1)
    2    continue
      endif
   10 continue

      if (zvec) write (iout1) iv1

      DEALLOCATE(ham1,eig1)
      ALLOCATE(ndim2d(ndimc),eigs2(max2d,ndimc),&
               ham2(max2d,max2d),eig2(max2d),iv2(ndimc),eigs2d(max3d))
      
!     now want to make the two-dimensional ham2(npntb,npntb',i,i'),
!     where i runs over the selected ham1(npnta,npnta') solutions.
      rewind ivecs1
      rewind ieigs1
      icall = 0
      low3d = 0
!     if zall then set zcut for convenience...
      if (zall) zcut = .true.
      do 31 ione = 1,npntc
      nham2 = 0
!     recall the size of ham2
      do 3 itwo = 1,npntb
      nham2 = nham2 + iv1(ione,itwo)
    3 continue
      ndim2d(ione) = nham2

      if ( nham2 .gt. 0 ) then
        call mkham2(ham2,eigs1d,vecs1d,xk1,xk2,iv1,ione,nham2)
        if (.not. ztwod) then
           call diag_dac(ham2,max2d,nham2,eig2)
        else
           call diag3d(ham2,nham2,eig2,kz)
           return
        endif

        if (zcut) then
          call cut2d(ham2,eig2,iv2(ione),nham2,low3d,icall)
        else
          do 4 ii = 1,nham2
            call outrow(ham2(1,ii),nham2,intvec)
            eigs2(ii,ione) = eig2(ii)
    4     continue
          low3d = max3d
        endif

     else
        iv2(ione) = 0
     endif
     
      if (ione .eq. npntc .and. zcut) write(6,987)  low3d,emax2
  987 format(/i14,' eigenvalues selected below ',d20.10)

   31 continue

      if (.not. zcut) call choose(eigs2,ndim2d,ham2,iv2,low3d)

      if (zvec) write (iout1) iv2

!      call timer

      nham3 = low3d

!     save the required bits to disk if zdiag = .false.
!     first open the required disk file
      if (.not. zdiag) then
        open(unit=idiag1,form='unformatted')
        open(unit=idiag2,form='unformatted')
        write(idiag1) npnta,npntb,npntc,max2d,nham3,zcut
        call ioutro(iv1,ndimb*ndimc,idiag1)
        call ioutro(iv2,ndimc,idiag1)

!       need also the 2-d eigenvalues to build the final hamiltonian
        if (zcut) then
          rewind ieigs2
          do 55 i1 = 1,npntc
          iv = iv2(i1)
          if (iv .gt. 0) call getrow(eig2,iv,ieigs2)
          if (iv .gt. 0) call outrow(eig2,iv,idiag1)
   55     continue
        else
          call outrow(eigs2,max2d*ndimc,idiag1)
        endif

      endif
      DEALLOCATE(eigs1d,eig2)
      ALLOCATE(hband(max2d,nham3),&
               cint(ndima*ndimb,max2d),cintp(ndima*ndimb,max2d))

      call mkham3(nham3,ham2,cint,cintp,iv1,iv2,ndim2d,&
                  vecs1d,xk1,xk2,eigs2d,eigs2,hband,xlmatr,r1,r2,&
                  r2m2t,term)

      DEALLOCATE(ham2,cint,cintp,vecs1d,hband,iv1)
      ALLOCATE(ham3(nham3,nham3))

      call loadh(ham3,nham3,iv2,ndim2d,eigs2d,eigs2)

      DEALLOCATE(iv2,ndim2d,eigs2d,eigs2)
      ALLOCATE(eval(nham3))

      call diag3d(ham3,nham3,eval,kz)

      DEALLOCATE(ham3,eval)

!.....finally, can compute the actual wavefunction amplitude at
!.....the grid points if needed.
      if (ztran) then
         ALLOCATE(evall(neval),&
                 iv1l(ndimc,ndimb),iv2l(ndimc),ndim2l(ndimc),&
                 vecs1l(max2d,ndima),vecs2l(max2d,max2d),vecs3l(nham3),&
                 phi(nham3,ndima,ndimb),psi(idvr,npnt1,npnt2))

         call trans(iv1l,iv2l,ndim2l,vecs1l,&
                    vecs2l,vecs3l,phi,psi,evall,nham3)

         DEALLOCATE(iv1l,iv2l,ndim2l,vecs1l,vecs2l,vecs3l,phi,psi,evall)
      endif

      return
      end

!#########################################################################
      subroutine mkham1(ham1,xlmatr,i1,i2,term,r1,r2,xalf,xk1,xk2,kz)

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2

      real*8, dimension(ndima,ndima) :: ham1
      real*8, dimension(idvr,idvr) :: xlmatr
      real*8, dimension(npnt1) :: r1
      real*8, dimension(npnt2) :: r2
      real*8, dimension(idvr) :: xalf
      real*8, dimension(npnt1,npnt1) :: xk1
      real*8, dimension(npnt2,npnt2) :: xk2

      data x0/0.0d0/,xp5/0.50d0/,x1/1.0d0/

!     zero ham1
      ham1 = x0
!     zero rotational excitation term for j=0 cases
      wterm = x0
!.....theta first
      if (ztheta) then

         w1gama = xp5 / (r1(i1)*r1(i1)*ur1)
         w2beta = xp5 / (r2(i2)*r2(i2)*ur2)
         wsum = w1gama + w2beta

         if (jrot .gt. 0) then
            if (zembed) then
!              have term * r2**(-2) term
               wterm = term * w2beta * ur2/urr2
            else
!              have term * r1**(-2) term
               wterm = term * w1gama * ur1/urr1
            endif
!     extra NBO term if vib mass .ne. rot mass
            if (kz .gt. 0) then
               s1 = ur1/urr1-x1
               s2 = ur2/urr2-x1
               w3 = dble(kz*kz) * (s1*w1gama + s2*w2beta)
            endif
         endif
      endif

      do 10 k = 1,npnta
      if (ztheta) then
         call potv(v,r1(i1),r2(i2),xalf(k))
!         write(91,1098)r1(i1),r2(i2),xalf(k),v
!         1098 format(3(2x,f16.10),1x,d25.18)
      else
         if (zr2r1)then
            call potv(v,r1(i2),r2(k),xalf(i1))
          else
            call potv(v,r1(k),r2(i1),xalf(i2))
         endif
      endif
      if (ztheta) then
         if (kz .gt. 0) v = v + w3/(x1-xalf(k)**2)
         ham1(k,k) = v + wterm
         do 20 kp= 1,k
         ham1(k,kp) = ham1(k,kp) + xlmatr(k,kp)*wsum
   20    continue    

      else
         if (jrot .gt. 0) then
            if (zembed) then
              if (zr2r1)then
                 wterm = (term * xp5) / (r2(k)*r2(k)*ur2)
              else
                 wterm = (term * xp5) / (r2(i1)*r2(i1)*ur2)
              endif
            else
              if (zr2r1)then
                 wterm = (term * xp5) / (r1(i2)*r1(i2)*ur1)
              else
                 wterm = (term * xp5) / (r1(k)*r1(k)*ur1)
              endif
            endif
         endif
         ham1(k,k) = v + wterm

         if (zr2r1) then
            do 30 kp= 1,k
            ham1(k,kp) = ham1(k,kp) + xk2(k,kp)
   30       continue
         else
            do 40 kp= 1,k
            ham1(k,kp) = ham1(k,kp) + xk1(k,kp)
   40       continue
         endif
      endif
   10 continue

      return
      end

!####################################################################
      subroutine mkham2(ham2,eigs1d,vecs1d,xk1,xk2,iv1,ione,nham2)
      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(max2d,max2d) :: ham2
      real*8, dimension(npnt2,npnt2) :: xk2
      dimension iv1(ndimc,ndimb)
      real*8, dimension(max2d) :: eigs1d
      real*8, dimension(max2d,ndima) :: vecs1d
      real*8, dimension(npnt1,npnt1) :: xk1

!     zero ham2
      ham2 = 0.0d0

      do 10 i3= 1,npnta
         call getrow(vecs1d(1,i3),nham2,ivecs1)
         do 20 j = 1,nham2
           do 30 k = 1,j
             ham2(j,k) = ham2(j,k) + vecs1d(j,i3)*vecs1d(k,i3)
   30      continue
   20    continue
   10 continue
!     must now multiply by xk1 or xk2
      ivbsm = 0
      do 50 itwo = 1,npntb
        ivbpsm = 0
        ivb = iv1(ione,itwo)
        do 60 itwop = 1,itwo
          if (ztheta) then
            if (zr2r1) then
              xkterm = xk2(itwo,itwop)
            else
             xkterm = xk1(itwo,itwop)
            endif
          else
           if (zr2r1) then
           xkterm = xk1(itwo,itwop)
        else
           xkterm = xk2(itwo,itwop)
        endif
      endif
      ivbp = iv1(ione,itwop)
      do 70 j  = 1,ivb
        ind1 = ivbsm + j
        do 77 jp = 1,ivbp
          ind2 = ivbpsm + jp
          ham2(ind1,ind2) = ham2(ind1,ind2) * xkterm
   77   continue
   70 continue
      ivbpsm = ivbpsm + ivbp
   60 continue
      ivbsm = ivbsm + ivb
   50 continue

!     now add the 1-d eigenvalues along the diagonal
      call getrow(eigs1d,nham2,ieigs1)
      do 80 nn = 1,nham2
         ham2(nn,nn) = ham2(nn,nn) + eigs1d(nn)
   80 continue
      return
      end

!################################################################################
      subroutine mkham3(nham3,ham2,cint,cintp,iv1,iv2,ndim2d,&
                        vecs1d,xk1,xk2,eigs2d,eigs2,hband,xlmatr,&
                        r1,r2,r2m2t,term)

!     build the final 3-d hamiltonian matrix.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2

      real*8, dimension(ndima*ndimb,max2d) :: cint
      real*8, dimension(max2d,max2d) :: ham2
      dimension iv1(ndimc,ndimb),iv2(ndimc),ndim2d(npntc)
      real*8, dimension(max2d,ndima) :: vecs1d
      real*8, dimension(ndima*ndimb,max2d) :: cintp
      real*8, dimension(npnt1,npnt1) :: xk1
      real*8, dimension(nham3) :: eigs2d
      real*8, dimension(npnt2,npnt2) :: xk2
      real*8, dimension(max2d,nham3) :: hband
      real*8, dimension(max2d,ndimc) :: eigs2
      real*8, dimension(idvr,idvr) :: xlmatr
      real*8, dimension(npnt1) :: r1
      real*8, dimension(npnt2) :: r2
      real*8, dimension(nham3) :: work3
      real*8, dimension(npnt2,npnt2) :: r2m2t

      data xp5/0.50d0/

!     if zdiag = .false. want eigs2d now
      if (.not. zdiag) then
        rewind ieigs2
        neig = 0
        do 184 ione = 1,npntc
          iv = iv2(ione)
          if (.not. zcut) then
            do 19 ii = 1,iv
              neig = neig + 1
              eigs2d(neig) = eigs2(ii,ione)
   19       continue
          endif
          if (iv .gt. 0 .and. zcut) call getrow(eigs2d,iv,ieigs2)
  184   continue
      endif

!     first do the intermediate transformation
      ncint = npnta*npntb
      ndimt = ndima*ndimb
      rewind ivecs2
      rewind ivecs1

      do 10 ione = 1,npntc
!        recall the size of the 2-d vectors
         nham2=ndim2d(ione)
         if (nham2 .eq. 0) goto 10
         cint = 0.0d0
!     bring back the 1-d vectors for each i1
         do 23 kk = 1,npnta
           if (iv2(ione) .gt. 0) then
             call getrow(vecs1d(1,kk),nham2,ivecs1)
           else
             read (ivecs1)
           endif
   23    continue
         do 20 j = 1,iv2(ione)
!     bring back the 2-d vectors for each npntc
           call getrow(ham2(1,j),nham2,ivecs2)
           ind2 = 0
           do 30 k = 1,npnta
             ind1 = 0
             do 40 itwo = 1,npntb
               ind2 = ind2 + 1
               do 50 i = 1,iv1(ione,itwo)
                 ind1 = ind1 + 1
                 cint(ind2,j) = cint(ind2,j) + ham2(ind1,j)*vecs1d(ind1,k)
   50          continue
   40        continue
   30      continue
           if (.not. ztheta) then
             in2 = 0
             do 37 ia = 1,npnta
               do 47 ib = 1,npntb
                 if (zr2r1) then
                    ii1 = ib
                    ii2 = ia
                 else
                    ii1 = ia
                    ii2 = ib
                 endif
                 w1 = xp5 / (r1(ii1)*r1(ii1)*ur1)
                 w2 = xp5 / (r2(ii2)*r2(ii2)*ur2)
                 wsum = w1 + w2
                 in2 = in2 + 1
                 cint(in2,j) = cint(in2,j)*sqrt(wsum)
   47          continue
   37        continue
           endif
   20    continue
!     store cint on disk for each npntc
         if (iv2(ione) .gt. 0) call outrow(cint,ndimt*iv2(ione),ivint)
   10 continue

!     now do the second part of the transformation
!      endfile ivint


      rewind ivint
      if (ztheta) then
        rm2t = 0.0d0
        length = 0
        do 61 ione = 1,npntc
!        set the hband to zero
         hband = 0.0d0
         ivsm = 0
         nham2 = ndim2d(ione)
         if (nham2 .eq. 0) goto 61
         iv = iv2(ione)
         if (iv .gt. 0) call getrow(cint,ndimt*iv,ivint)
         rewind ivint
         ivpsm = 0
         do 51 ionep = 1,ione
           if (.not.zquad2) then
             if (ione .eq. ionep) then
               d1r2 = r2m2t(ione,ionep)
               d2r2 = xp5 / (r2(ione)*r2(ione)*ur2)
               rm2t = d1r2 - d2r2
             else
               rm2t = r2m2t(ione,ionep)
             endif
           endif
           if (zr2r1) then
             xkterm = xk1(ione,ionep)
           else
             xkterm = xk2(ione,ionep)
             if(jrot .gt. 0 .and. .not. zquad2 .and. zembed) xkterm = xkterm + term*rm2t
           endif
           ivp = iv2(ionep)
           if (ivp .gt. 0) call getrow(cintp,ndimt*ivp,ivint)
           do 41 j=1,iv
             ind1 = ivsm + j
             do 31 jp =1,ivp
               ind2 = ivpsm + jp
               ind3 = 0
               do 81 k = 1,ncint
                  ind3 = ind3 + 1
                  hband(ind1,ind2) = hband(ind1,ind2)&
                                     + cint(ind3,j)*cintp(ind3,jp)*xkterm
   81          continue
               ind3 = 0
               do 83 na = 1,npnta
                  do 82 k = 1,npntb
                    ind3 = ind3 + 1
                    ind4 = k
                    do 84 nap = 1,npnta
                       hband(ind1,ind2) = hband(ind1,ind2) + cint(ind3,j)*cintp(ind4,jp)&
                                           * xlmatr(na,nap) * rm2t
                       ind4 = ind4 + npntb
   84               continue
   82             continue
   83          continue
   31        continue
   41      continue
           ivpsm = ivpsm + ivp
   51    continue
         if (zdiag) then
           do 71 jj=1,ind2
             if (iv .gt. 0) call outrow(hband(1,jj),iv,iband)
   71      continue
           ndim2d(ione)=ind2
         else
           do 93 isave = 1,iv
             length = length + 1
             do 94 jsave = 1,ind2
               work3(jsave)= hband(isave,jsave)
   94        continue
             work3(length) = work3(length) + eigs2d(length)
             call outrow(work3,length,idiag2)
   93      continue
         endif
   61 continue

      else
        length = 0
        do 66 ione = 1,npntc
!         set the hband to zero
          hband=0.0d0
          ivsm = 0
          if (ndim2d(ione) .eq. 0) goto 66
          iv = iv2(ione)
          if (iv .gt. 0) call getrow(cint,ndimt*iv,ivint)
          rewind ivint
          ivpsm = 0
          do 56 ionep = 1,ione
            xkterm = xlmatr(ione,ionep)
            ivp = iv2(ionep)
            if (ivp .gt. 0) call getrow(cintp,ndimt*ivp,ivint)
            do 46 j=1,iv
              ind1 = ivsm + j
              do 36 jp =1,ivp
                ind2 = ivpsm + jp
                ind3 = 0
                do 86 k = 1,ncint
                  ind3 = ind3 + 1
                  hband(ind1,ind2) = hband(ind1,ind2)&
                               + cint(ind3,j)*cintp(ind3,jp)*xkterm
   86           continue
   36         continue
   46       continue
            ivpsm = ivpsm + ivp
   56     continue
          if (zdiag) then
            do 76 jj=1,ind2
              if (iv .gt. 0) call outrow(hband(1,jj),iv,iband)
   76       continue
            ndim2d(ione)=ind2
          else
            do 95 isave = 1,iv
              length = length + 1
              do 96 jsave = 1,ind2
                work3(jsave)= hband(isave,jsave)
   96         continue
            work3(length) = work3(length) + eigs2d(length)
            call outrow(work3,length,idiag2)
   95       continue
          endif
   66  continue
      endif

      if (.not. zdiag) then
         write(6,1080)
 1080    format(/5x,'hamiltonian written to disk - not diagonalised')
         write(6,1081) idiag2
 1081    format(/5x,'hamiltonian bands written to stream idiag =',i4/)
         stop
       endif
       return
       end

      subroutine loadh(ham3,nham3,iv2,ndim2d,eigs2d,eigs2)

!     load the final 3-d hamiltonian matrix.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      dimension iv2(ndimc),ndim2d(npntc)
      real*8, dimension(nham3,nham3) :: ham3
      real*8, dimension(nham3) :: eigs2d
      real*8, dimension(max2d,ndimc) :: eigs2

      ham3 = 0.0d0
      rewind iband
      ivsm = 0
      do 62 ione = 1,npntc
        if (ndim2d(ione) .eq. 0) goto 62
        iv = iv2(ione)
        do 72 jj=1,ndim2d(ione)
          if (iv .gt. 0) call getrow(ham3(ivsm + 1,jj),iv,iband)
   72   continue
        ivsm = ivsm + iv
   62 continue

!     now put the 2-d eigensolutions along the diagonal

      rewind ieigs2
      nn = 0
      do 183 ione = 1,npntc
       iv = iv2(ione)
       if (.not. zcut) then
       do 9 ii = 1,iv
          eigs2d(ii) = eigs2(ii,ione)
    9  continue
       endif
       if (iv .gt. 0 .and. zcut) call getrow(eigs2d,iv,ieigs2)
       do 185 i = 1,iv
         nn = nn + 1
         ham3(nn,nn) = ham3(nn,nn) + eigs2d(i)
  185  continue
  183 continue

      if (zpham) call wrtham(ham3,nham3)

      return
      end

!########################################################################
      SUBROUTINE diag(ham,maxham,nham,eig)

!     diagonalise the appropriate hamiltonian matrices
      implicit real*8 (a-h,o-y), logical (z)
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(maxham,nham) :: ham
      real*8, dimension(nham) :: eig
      real*8, dimension(maxham*3) :: work
      
      ifail=0
      nnham=maxham*3
      call dsyev ('V','L',nham,ham,maxham,eig,work,nnham,ifail)

      if (ifail .ne. 0) write(6,100) ifail
      return
100   format(' diagonalisation has failed with, ifail=',i3)
    END SUBROUTINE diag

!########################################################################
      subroutine diag_dac(ham,maxham,nham,eig)
        ! Diag Divide And Conquer
        ! new routine uses the faster dysevd lapack routine
        ! but does require an extra N**2 workspace,
        ! which for very large runs may be restrictive
        implicit none
        ! external vars
        integer :: maxham,nham,j,i
        real*8 :: ham(maxham,nham),phase
        real*8 :: eig(nham)
        ! internal vars
        integer :: ifail,lwork,liwork
        ! allocatable arrays
        real*8, allocatable ::  work(:)
        integer, allocatable :: iwork(:)

        lwork = 1 + 6*nham + 2*nham**2
        liwork = 3 + 5*nham
        ifail = 0

        allocate(iwork(liwork),work(lwork))

!        write(6,*)'***  ',nham,maxham

        call dsyevd('V','L',nham,ham,maxham,eig,work,lwork,iwork,liwork,ifail)
!!!!!!!!!!!!!!!! pb debugging !!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!        do i=1,nham
!            phase=((-1.d0)**(int(1.d0+sin(real(i)))))
!            phase=1.d0
!            write(6,*)'Phase 3d ... ',i,phase
!           do j=1,nham
!              ham(j,i)=phase*ham(j,i)
!           end do
!        end do

        deallocate(work,iwork)
        
        if (ifail .ne. 0) write(6,100) ifail
        return
100     format(' diagonalisation has failed with, ifail=',i3)
      end subroutine diag_dac

!########################################################################
      subroutine diag_rrr(ham,maxham,nham,eig,zcut,emax2d)
        ! Diag Relativly Robust Representation
        ! even faster then divide and conquer, and uses least workspace
        ! could be programmed to only do neval solns, but i havent here
        implicit none
        ! external vars
        integer :: maxham,nham,max2d,i,j
        real*8 :: ham(maxham,nham)
!        real*8 :: ham_tmp(maxham,nham)
        real*8 :: eig(nham)
        ! internal vars
        character ajob*1
        logical zcut
        real*8 :: vl,vu,abstol,work_query,emax2d,phase
        integer :: il,iu,m,lwork,liwork,ifail,iwork_query
        real*8  :: z(nham,maxham)
        real*8,allocatable :: work(:)
        integer :: isuppz(2*nham)
        integer,allocatable :: iwork(:),ifailv(:)
        
        abstol = 0.0d0
        ifail = 0
        eig=0.d0
        z=0.d0
                        
        ! now do it
        if (.not.zcut ) then
           ajob='A'
!           il=1
!           iu=nham
        ! first ask for work space requirements
        lwork = -1
        liwork = -1
!        call dsyevr('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
!             m,eig,z,maxham,isuppz,work_query,lwork,iwork_query,liwork,ifail)
        allocate (iwork(5*nham))
        allocate (ifailv(nham))
        call dsyevx('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
             m,eig,z,maxham,work_query,lwork,iwork_query,ifailv,ifail)
        if (ifail .ne. 0) write(6,100) ifail
        lwork = int(work_query)
        allocate(work(lwork))
        liwork = iwork_query
!        allocate (iwork(liwork))

        else 
           ajob='V'
           vl=-emax2d
           vu=emax2d

        ! first ask for work space requirements
        lwork = -1
        liwork = -1
!        call dsyevr('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
!             m,eig,z,maxham,isuppz,work_query,lwork,iwork_query,liwork,ifail)
        allocate (iwork(5*nham))
        allocate (ifailv(nham))
        call dsyevx('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
             m,eig,z,maxham,work_query,lwork,iwork_query,ifailv,ifail)
        if (ifail .ne. 0) write(6,100) ifail
        lwork = int(work_query)
        allocate(work(lwork))
        liwork = iwork_query
!        allocate (iwork(liwork))
        end if

!        ham_tmp=ham
!        call dsyevx('V','A','L',nham,ham_tmp,maxham,vl,vu,il,vu,abstol,&
!             m,eig,z,maxham,work,lwork,iwork,ifailv,ifail)
!        do i=1,nham
!           write(69,101)i,eig(i)
!        end do
101     format(I5,d25.18)

!        call dsyevr('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
!             m,eig,z,maxham,isuppz,work,lwork,iwork,liwork,ifail)
        call dsyevx('V',ajob,'L',nham,ham,maxham,vl,vu,il,vu,abstol,&
             m,eig,z,maxham,work,lwork,iwork,ifailv,ifail)
        if (ifail .ne. 0) write(6,100) ifail

        if (zcut) then
        do i=m+1,nham
           eig(i)=2.d0*vu
        end do
        end if

        deallocate(work)
        deallocate(iwork)
        ! copy eigen vectors onto ham
        CALL dcopy(nham*maxham,z,1,ham,1)

        return
100     format(' diagonalisation has failed with, ifail=',i3)
      end subroutine diag_rrr

!########################################################################
      subroutine diag3d(ham3,nham3,eval,kz)

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nham3,nham3) :: ham3
      real*8, dimension(nham3) :: eval
      real*8, dimension(nham3) :: evalcm

!     autocm converts atomic units (hartree) to cm-1.
      data autocm/2.19474624d+05/
      data x0/0.0d0/
      if (zrot) then
         write(6,1040) jrot,kz
 1040    format(/5x,'Solutions with J =',i3,' k =',i3)
         if (idia .eq. -2) then
            if (ipar .eq. 0)  write(6,1025)
 1025       format(/5x,'even parity solutions')
            if (ipar .eq. 1)  write(6,1035)
 1035       format(/5x,'odd parity solutions')
         endif
      endif
!     we need this in special cases (nham3 < neval)
      meval=min(neval,nham3)
      if (zvec) write(iout2) meval

!     diagonalise  the hamiltonian. the vectors are to overwrite the
!     hamiltonian.
      call diag_dac(ham3,nham3,nham3,eval)

! writing J=0 output for dipole in stream idip

      if (jrot.eq.0.or.(jrot.eq.1.and.kz.eq.0)) then
         write(idip)meval                           !6
         write(idip)eval                            !7
         write(idip)kz,ipar,ipar,idvr,ipar          !8
         write(idipd)meval                          !6
         write(idipd)eval                           !7
         write(idipd)kz,ipar,ipar,idvr,ipar         !8
      end if

!     print eigenvalues in atomic untis & wavenumbers

      if (.not.zpmin) then
         write(6,1010) meval
 1010    format(/5x,'lowest',i4,' eigenvalues in hartrees:',/)
         write(6,1020) (eval(i),i=1,meval)
      endif
      if (zpfun) then
         ip=jrot-kmin
         if (jrot .ne. 0 .and. ip .ne. 1) goto 10
         jdia=max(0,idia)
         jpar=min(jdia,ipar)
         isym=abs(min(0,idia))
         if (ipar .eq. 1) isym=-isym
         write(ilev,1125) jrot,ip,jdia,jpar,isym,meval
 1125    format(6i4)
         write(ilev,1126) (eval(i),i=1,meval)
 1126    format(4d20.12)
      endif
   10 continue

!     save the eigenvalues if needed
      if (zvec) write (iout2) (eval(i),i=1,meval)
      do 20 i=1,meval
      evalcm(i) = eval(i) * autocm
   20 continue
      write(6,1030) meval
 1030 format(//5x,'Lowest',i4,' eigenvalues in wavenumbers:'/)
      write(6,1020) (evalcm(i),i=1,meval)
 1020 format(5d24.12/)

!     if requested print the eigenvectors
      if (zpvec) then
         write(6,1050)
 1050    format(//'  Eigenvectors'/)
         do 40 i=1,meval
         write(6,1060) (ham3(j,i),j=1,nham3)
 1060    format((1x,10f13.7))
   40    continue
      endif
!     write the final vectors to disk if required
      if (zvec) then
         do 60 l = 1,meval
         call outrow(ham3(1,l),nham3,iout2)
 1061    format((2I6,d20.10))
   60    continue
      endif
      if (jrot .ne. 0) return
      if (abs(idia) .eq. 2 .and. ipar .eq. 1) then
         ii=1
         ezero=x0
         read(5,5,end=55) ezero
    5    format(f20.0)
   55    continue
      else
         ezero=evalcm(1)
         ii=2
      endif
      write(6,1070)
 1070 format(//'  Band origins in wavenumbers:'/)
      write(6,1021) (evalcm(i)-ezero,i=ii,meval)
 1021 format(5f15.6/)

      return
      end

!########################################################################
      subroutine choose(eigs2,ndim2d,ham2,iv2,low3d)

!     this routine chooses the max3d lowest eigenvalues from eigs2.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(max2d,ndimc) :: eigs2
      dimension iv2(ndimc)
      dimension ndim2d(ndimc)
      real*8, dimension(max2d,max2d) :: ham2

      data autocm/2.19474624d+05/

      eigmin = eigs2(1,1)
      nhamsm = 0
      do 160 i=1,npntc
        eigmin = min(eigmin,eigs2(1,i))
        nhamsm = nhamsm + ndim2d(i)
        iv2(i) = 1
  160 continue
      if (nhamsm .lt. max3d) low3d = nhamsm
      ipt = 1
      do 200 n=1,low3d
  210 if(iv2(ipt) .le. ndim2d(ipt)) then
          eigvib = eigs2(iv2(ipt),ipt)
          jpt = ipt
      else
          ipt = ipt + 1
          goto 210
      endif
      do 220 j=ipt+1,npntc
      if(iv2(j) .gt. ndim2d(j)) goto 220
      if(eigs2(iv2(j),j) .ge. eigvib) goto 220
      eigvib = eigs2(iv2(j),j)
      jpt = j
  220 continue
!     keep the eigenvalue:
      iv2(jpt) = iv2(jpt) + 1
  200 continue
!     store the number of eigenvalues selected for each alpha
      iv2(1) = iv2(1) - 1
      ivib = iv2(1)
      write(6,800)
  800 format(//5x,'Selection outcome for the final, contracted basis:'/)
      write(6,900)  1,iv2(1)
      do 230 i=2,npntc
       iv2(i) = iv2(i) - 1
       write(6,900)  i,iv2(i)
  900  format(5x,'npntc =',i3, ',', ' no. of eigenvectors =',i3 )
       ivib = max(ivib,iv2(i))
  230 continue

      write(6,998)  low3d,eigmin,eigvib
  998 format(/i14,' eigenvalues selected from ',d20.10,' to',d20.10,&
          ' hartrees')
      write(6,999)  low3d,eigmin*autocm,eigvib*autocm
  999 format(/i14,' eigenvalues selected from ',d20.10,' to',d20.10,&
          ' cm-1')


      if (zp2d) then
         write(6,1051)
 1051    format(//5x,'2d eigenvalues in wavenumbers:'/)
         do 32 jone = 1,npntc
         ivj = iv2(jone)
         write(6,1050) (eigs2(i,jone)*autocm,i=1,ivj)
 1050    format(5d24.12/)
         write(6,1052)
 1052    format(/5x,'------------------------------'/)
   32    continue
      endif      

!     save the vectors chosen
      rewind intvec
      do 31 ione = 1,npntc
      ivm = iv2(ione)
      nham2 = ndim2d(ione)
      do 7 i = 1,nham2
        call getrow(ham2(1,i),nham2,intvec)
    7 continue
      do 30 ind = 1,ivm
      if (zvec) then
         if (nham2 .gt. 0) call outrow(ham2(1,ind),nham2,iout2)
      endif
         if (nham2 .gt. 0) call outrow(ham2(1,ind),nham2,ivecs2)
   30 continue
   31 continue

      return
      end

!#####################################################################
      subroutine cut1d(ham1,eig1,ivn,eigs1d,vecs1d,nham2,icall)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax1, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     to disk.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(ndima,ndima) :: ham1
      real*8, dimension(ndima) :: eig1
      real*8, dimension(max2d) :: eigs1d
      real*8, dimension(max2d,ndima) :: vecs1d

      data autocm/2.19474624d+05/
!     change emax1 to hartree for the selection
      emaxau=emax1/autocm

      icall = icall + 1
      if (.not. zall) then
         ivn = 0
         do 10 n=1,npnta
         if (eig1(n) .gt. emaxau) goto 20
         ivn = ivn + 1
   10    continue
   20    iv = ivn
      else
         ivn = npnta
         iv  = ivn
      endif

      if (zp1d) then
         if (icall .eq. 1) write(6,1051)
 1051    format(//5x,'1d eigenvalues in wavenumbers:'/)
         write(6,1050) (eig1(i)*autocm,i=1,npnta)
 1050    format(5d24.12/)
         write(6,1052)
 1052    format(/5x,'------------------------------'/)
      endif

!     save the vectors and eigenvalues (overwrite for each gamma).
      do 30 i=1,iv
         nham2 = nham2 + 1

      if (nham2 .gt. max2d .and. .not. zall)  then
         write(6,999)
  999    format(//6x,'**** core exceeded: reduce cut-off emax1 ****')
         stop
      endif

         eigs1d(nham2)   = eig1(i)
      do 31 j=1,npnta
         vecs1d(nham2,j) = ham1(j,i)
   31 continue
   30 continue

      return
      end

!####################################################################
      subroutine cut2d(ham2,eig2,ivm,nham2,low3d,icall)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax2, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     to disk.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(max2d,max2d) :: ham2
      real*8, dimension(nham2) :: eig2

      data autocm/2.19474624d+05/
!     change emax2 to hartree for the selection
      emaxau=emax2/autocm

      icall = icall + 1
      if (.not. zall) then
         ivm = 0
         do 10 n=1,nham2
         if(eig2(n) .gt. emaxau) goto 20
         ivm = ivm + 1
   10    continue
   20    low3d = low3d + ivm
      else
         ivm = nham2
         low3d = low3d + ivm
      endif

      if (low3d .gt. max3d .and. .not. zall)  then
         write(6,999)
  999    format(//6x,'**** core exceeded: reduce cut-off emax2 ****')
         stop
      endif

      if (zp2d) then
         if (icall .eq. 1) write(6,1051)
 1051    format(//5x,'2d eigenvalues in wavenumbers:'/)
         write(6,1050) (eig2(i)*autocm,i=1,nham2)
 1050    format(5d24.12/)
         write(6,1052)
 1052    format(/5x,'------------------------------'/)
      endif

!     save the vectors and eigenvalues
      if (ivm .gt. 0) call outrow(eig2,ivm,ieigs2)
      do 30 ind = 1,ivm
      if (zvec) call outrow(ham2(1,ind),nham2,iout2)
      call outrow(ham2(1,ind),nham2,ivecs2)
   30 continue
      return
      end

!########################################################################
      subroutine nfmain(hr,r1m2t,htheta,r,theta,kz)

!     this routine controls the dvr calculation in the case of
!     symmetrised radau coordinates.
!     written by nic fulton, feb 1993.

      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(npnt,npnt) :: hr
      real*8, dimension(npnt,npnt) :: r1m2t
      real*8, dimension(nalf,nalf) :: htheta
      real*8, dimension(npnt) :: r
      real*8, dimension(nalf) :: theta
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig2,eigtmp
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::vecs2d
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eigs2d
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham3
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig3
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv2
      INTEGER, ALLOCATABLE, DIMENSION(:) ::nv2
      data x4/4.0d0/,x8/8.0d0/,x16/1.6d1/

      ALLOCATE(ham2(max2d,max2d),eig2(max2d),iv2(2,nalf),&
               vecs2d(max2d,max3d),eigs2d(max3d),nv2(nalf))
      if (.not.zcut) ALLOCATE(eigtmp(nalf*max2d))
      if (jrot .ne. 0) then
         term  = dble(jrot * jrot + jrot - kz * kz) / x8
         term2 = dble(jrot * jrot + jrot - 3 * kz * kz) / x4
         if (abs(kz) .eq. 1) then
            term3 = dble(jrot * jrot + jrot) / x16
            if (kmin .ge. 1 .and. zrot) term3 = -term3
         endif
      endif

! no need for the 1d diagonalisations as there is no possiblity of
! truncation as the symmetry would be broken.
      write(6,130)
      call nftim('beginning of 2d loop')
      do 30 igamma = 1,nalf
  
        nv2(igamma) = (npnt*(npnt+1-ipar*2))/2
        nham2 = nv2(igamma)
        call blc2d1(theta(igamma),r,hr,ham2,nham2,&
                          term,term2,term3,kz)
        call diag(ham2,nham2,nham2,eig2)

        if (.not. zcut) then
          call choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,&
                      vecs2d,nv2,eigtmp)

          if(igamma .eq. nalf .and. .not. zpmin)&
              write(6,110) (itmp, iv2(2,itmp),itmp=1,nalf)
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

      if (.not.zcut) DEALLOCATE(eigtmp)
      ALLOCATE(ham3(nham3,nham3))

      call bloc3d(htheta,r1m2t,ham3,eigs2d,vecs2d,iv2,nv2,ham2,nham3,r)

      DEALLOCATE(ham2,eig2,eigs2d,nv2)
 
      call nftim('end of 3d ham building')
      write(6,170)
      ALLOCATE(eig3(nham3))

      call diag3d(ham3,nham3,eig3,kz)

      call nftim('end of diagonalising 3d')

      if (ztran) call transr(iv2,vecs2d,ham3,eig3,nham3,nbass)

      DEALLOCATE(iv2,vecs2d,ham3,eig3)
      return
110   format(5x,'for gamma = ',i2,' selected ',i3,' energies.')
120   format(25x,'total = ',i5)
130   format(5x,'starting the 2d calculations.')
160   format(5x,'building the 3d hamiltonian.')
170   format(5x,'diagonalising the 3d hamiltonian.')
      end

!    ***********************************************************************

      subroutine blc2d1(xcos,r,hr,ham2,nham2,term,term2,term3,kz)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1

      real*8, dimension(npnt,npnt) :: hr
      real*8, dimension(npnt) :: r
      real*8, dimension(nham2,nham2) :: ham2
      data xp5/0.5d0/,x1/1.0d0/
      factr2 = sqrt(xp5)

      ham2 = 0.0d0

!     allow for j > 0 case
      if (jrot .ne. 0) then
        fact =  term + term2 / (x1 - xcos)
        if (kz .eq. 1) fact = fact + term3 * (x1 + xcos)/(x1-xcos)
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

!      q=0.0d0
      q=1.0d0
      if(ipar .eq. 1) q=-1.0d0
      iap=0
      do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ibeta .eq. ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+hr(ialphp,ialpha)
              if(ibeta .eq. ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*hr(ibetap,ialpha)
              if(ialpha .eq. ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*hr(ialphp,ibeta)
              if(ialpha .eq. ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+hr(ibetap,ibeta)
              if(ialpha .eq. ialphp .and. ibeta .eq. ibetap) then
                call potv(v,r(ibeta),r(ialpha),xcos)
                ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+v
              endif
              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
                call potv(v,r(ibeta),r(ialpha),xcos)
                ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*v
              endif
              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
      return
      end

!     ***********************************************************************

      subroutine choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d,nv2,eigtmp)

!     this routine chooses the max3d lowest eigenvalues from eigs2.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nham2) :: eig2
      real*8, dimension(nham2,nham2) :: ham2
      dimension iv2(2,nalf)
      real*8, dimension(max3d) :: eigs2d
      real*8, dimension(nalf*max2d) :: eigtmp
      real*8, dimension(max2d,max3d) :: vecs2d
      dimension nv2(nalf)
      save itotal
      data autocm/2.19474624d+05/

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

100   format(5x,'selecting ',i5,' energies between ',e12.6,&
             ' and ',e12.6,' in hartrees')
110   format(5x,'selecting ',i5,' energies between ',e12.6,&
             ' and ',e12.6,' in wavenumbers')
      return
      end

!     ***********************************************************************

      subroutine cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax1, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     in the array vecs1d.

      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nham2) :: eig2
      real*8, dimension(nham2,nham2) :: ham2
      dimension iv2(2,nalf)
      real*8, dimension(max3d) :: eigs2d
      real*8, dimension(max2d,max3d) :: vecs2d
      save npos
      data autocm /2.19474624d+05/
      if (igamma .eq. 1) npos = 1

!     change emax2 to hartree for the selection
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
!          write(6,*)'deb ',nvec,ntot,vecs2d(1,ntot),vecs2d(nham2,ntot)
        endif
10    continue
      iv2(1,igamma) = npos
      iv2(2,igamma) = nvec
      npos = npos + nvec
      return
140   format('number of 2d eigenvalues greater than ',i4,&
              ' increase max3d or reduce emax2 which is ',e12.5,&
              ' cm-1.')
      end

!##################################################################
      subroutine testiv(iv,nbass)

!     selection vectors for the bisector embedding to ensure that
!     singular region of theta = 0 is not sampled when j > 0.
!     also calculate which angular grid points are redundant.

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
     
      dimension iv(2,nalf)

      if (jrot .gt. 0) then
!        first find the extent of the functions in low theta direction
         do 20 iend=nalf/2,nalf
         if (iv(2,iend) .eq. 0) goto 30
   20    continue
         iend=nalf
!        have we saved functions beyond the point of zero amplitude ?
   30    ioff=0
         do 40 i=iend+1,nalf
!        if so remove them
         if (iv(2,i) .gt. 0) then
            ioff=ioff+iv(2,i)
            iv(2,i)=0
            iv(1,i)=iv(1,i-1)
         endif
   40    continue

!         zlin = .true.

         if (ioff .gt. 0) write(6,960) ioff,iv(1,iend)
  960    format(/5x,'*** warning:',i4,' functions removed from theta =',&
                ' 0 region'/8x,' basis reset to nham3 =',i5)
         if (iv(2,nalf) .gt. 0) then
            if (zlin) then
!!              if there are still theta=0 functions, remove them
               write(6,987) iv(2,nalf),iv(1,nalf)-1
               iv(2,nalf)=0
            else
!!              wavefunction has amplitude all the way to theta = 0: stop
               write(6,950)
               stop
            endif
         endif
  950    format(/5x,'bisector embedding: ',&
                    'wavefunction has amplitude for theta = 0. stop.')
  987    format(/5x,'*** warning: zlin = t forced the removal of',i4,&
               ' functions at theta = 0'&
                /8x,' basis reset to  nham3 =',i5)
       endif

      iang=0
      do 60 ii=1,nalf
        if (iv(2,ii) .gt. 0) iang=iang+1
   60 continue
      nbass=iang*max2d
      if (ztran) then
         write(iwave) iang,nbass
         write(iwave) (iv(2,ii),ii=1,nalf)
      endif
      return
      end

!***********************************************************************

      subroutine bloc3d(htheta,r1m2t,ham3,eigs2d,vecs2d,iv2,nv2,ham2,nham3,r)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot

      real*8, dimension(nalf,nalf) ::  htheta
      real*8, dimension(nham3,nham3) ::  ham3
      real*8, dimension(max3d) :: eigs2d
      real*8, dimension(max2d,max3d) :: vecs2d
      dimension iv2(2,nalf),nv2(nalf)
      real*8, dimension(max2d,max2d) :: ham2
      real*8, dimension(npnt) :: r

! zero ham3
      ham3 = 0.0d0

      ndim1g = 1
      do 10 igamma = 1,nalf
        if(iv2(2,igamma) .gt. 0) then
          ndim2g = 1
          do 20 igammp = 1,igamma
             ht=htheta(igammp,igamma)
            if(iv2(2,igammp) .gt. 0) then
               if (zquad2) then
              call blc2d2(r,igamma,igammp,ht,ham2,nv2(igamma))
              else
              call blc2d2_quad(r1m2t,r,igamma,igammp,ht,ham2,nv2(igamma))
              end if
              call vecmul(vecs2d(1,iv2(1,igamma)),nv2(igamma),&
                          iv2(2,igamma),vecs2d(1,iv2(1,igammp)),&
                          nv2(igammp),iv2(2,igammp),max2d,&
                          ham2,ham3(ndim1g,ndim2g),nham3)
              ndim2g = ndim2g + iv2(2,igammp)
            endif
20        continue
          ndim1g = ndim1g + iv2(2,igamma)
        endif
10    continue

      do 60 i = 1,nham3
        ham3(i,i) = ham3(i,i) + eigs2d(i)
60    continue
      if (zpham) call wrtham(ham3,nham3)

      return
      end

!    ***********************************************************************

      subroutine blc2d2(r,igamma,igammp,htheta,ham2,nham2)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1

!      real*8, dimension(nalf,nalf) ::  htheta
      real*8, dimension(npnt) :: r
      real*8, dimension(nham2,nham2) :: ham2
      data xp5 /0.5d0/
      factr2 = sqrt(xp5)
    
      ham2 = 0.0d0

!      q=0.0d0
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
                ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+htheta*wsum
              endif
              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
                wsum = walpha + wbeta
                ham2(ialphp+iap,ialpha+ia)=&
                 ham2(ialphp+iap,ialpha+ia)+q*htheta*wsum
              endif
              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
      return
      end

!    ***********************************************************************

      subroutine blc2d2_quad(r1m2t,r,igamma,igammp,htheta,ham2,nham2)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1

      real*8, dimension(npnt,npnt) ::  r1m2t
      real*8, dimension(npnt) :: r
      real*8, dimension(nham2,nham2) :: ham2
      data xp5 /0.5d0/
      factr2 = sqrt(xp5)
    
      ham2 = 0.0d0

!      q=0.0d0
      q=1.0d0
      if(ipar .eq. 1) q=-1.0d0
      iap=0
      do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar


              if(ialpha .eq. ialphp)&
                   ham2(ialphp+iap,ialpha+ia)&
                   =ham2(ialphp+iap,ialpha+ia)&
                   +htheta*r1m2t(ibetap,ibeta)
              if(ibeta .eq. ibetap)&
                   ham2(ialphp+iap,ialpha+ia)&
                   =ham2(ialphp+iap,ialpha+ia)&
                   +htheta*r1m2t(ialphp,ialpha)
              if(ibeta .eq. ialphp)&
                   ham2(ialphp+iap,ialpha+ia)&
                   =ham2(ialphp+iap,ialpha+ia)&
                   +q*htheta*r1m2t(ialpha,ibetap)
              if(ialpha .eq. ibetap)&
                   ham2(ialphp+iap,ialpha+ia)&
                   =ham2(ialphp+iap,ialpha+ia)&
                   +q*htheta*r1m2t(ibeta,ialphp)


!              if(ialpha .eq. ialphp .and. ibeta .eq. ibetap) then
!                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
!                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
!                wsum = walpha + wbeta
!                ham2(ialphp+iap,ialpha+ia)=&
!                  ham2(ialphp+iap,ialpha+ia)+htheta*wsum
!              endif
!              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
!                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
!                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
!                wsum = walpha + wbeta
!                ham2(ialphp+iap,ialpha+ia)=&
!                 ham2(ialphp+iap,ialpha+ia)+q*htheta*wsum
!              endif


              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
      return
      end

!     ***********************************************************************

      subroutine vecmul(veca,idima,jdima,vecb,idimb,jdimb,nvecln,hama,&
                        hamb,ndim)

! this routine does hb = hb + va^T * ha * vb
! could be replaced by blas?

      implicit real*8(a-h,o-y),logical(z)
      real*8, dimension(nvecln,jdima) ::  veca
      real*8, dimension(nvecln,jdimb) :: vecb
      real*8, dimension(idima,idimb) :: hama
      real*8, dimension(ndim,jdimb) :: hamb

      do 10 ib = 1,idimb
        do 20 ia = 1,idima
          temp1 = hama(ia,ib)
          if (temp1 .eq. 0.0d0) goto 20
          do 30 ja = 1,jdima
            temp2 = veca(ia,ja) * temp1
            do 40 jb = 1,jdimb
              hamb(ja,jb) = hamb(ja,jb) +  vecb(ib,jb) * temp2
40          continue
30        continue
20      continue
10    continue
      return
      end

!     ***********************************************************************

      subroutine transr(iv2,vecs2d,ham3,eig3,nham3,nbass)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
     
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /mass/ xmass(3),g1,g2,xmassr(3)

      dimension iv2(2,nalf)
      real*8, dimension(nalf,nalf) :: plegw
      real*8, dimension(max2d,max3d) :: vecs2d
      real*8, dimension(nham3,nham3) :: ham3
      real*8, dimension(nham3) :: eig3
      real*8, dimension(nbass) :: wvfunc
      real*8, dimension(nbass) :: wvfuncdvr

      call nftim('start of wavefunction generation')
      write(iwave) neval
      call outrow(eig3,neval,iwave)

!        write(idipd)(ham3(i,level3),i=1,nham3)  ! line 9 for J=0 dipole 

        if (ztran) then
           rewind(iwvpb)
           read(iwvpb)((plegw(i,j),i=1,nalf),j=1,nalf)
        end if

        nrad=(npnt*(npnt+1-2*ipar))/2
        nkbas=nrad*nalf
      do 40 level3=1,neval
        index=1
        wvfunc = 0.0d0
        wvfuncdvr = 0.0d0
        index = 1    
       do 30 irs=1,nrad
          do 300 igamma=1,nalf
            if (iv2(2,igamma) .eq. 0) goto 300
            do 20 level2=0,iv2(2,igamma)-1
              ivec2=iv2(1,igamma)+level2
              wvfuncdvr(index)=wvfuncdvr(index)+&
                vecs2d(irs,ivec2)*ham3(ivec2,level3)
20          continue
            index=index+1
300       continue
30      continue
        write(iwave) wvfuncdvr
!        write(6,*)'pb debugging ',jrot,kz
!        if (jrot.eq.0.or.(jrot.eq.1.and.kz.eq.0)) then
!        write(idipd) wvfuncdvr   ! line 9 for J=0 dipole 
!        call jtran(wvfunc,wvfuncdvr,plegw,nalf,nrad,nkbas)
!        write(idip) wvfunc       ! line 9 for J=0 dipole 
!        end if
40    continue

      call nftim('end of wavefunction generation')
      return
      end

!     ***********************************************************************

      subroutine transr2(iv2,vecs2d,ham3,eig3,nham3,nbass,kz,rloc)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
     
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /mass/ xmass(3),g1,g2,xmassr(3)

      dimension iv2(2,nalf)
      real*8, dimension(nalf,nalf) :: plegw
      real*8, dimension(max2d,max3d) :: vecs2d
      real*8, dimension(nham3,nham3) :: ham3
      real*8, dimension(nham3) :: eig3
      real*8, dimension(nbass) :: wvfunc
      real*8, dimension(nbass) :: wvfuncdvr
      real*8, dimension(npnt) :: rloc
      real*8, allocatable, dimension(:) :: wr,jxcos0,jwalf0,wmax,r,wf
      real*8, allocatable, dimension(:,:) :: wfs

!     we need this in special cases (nham3 < neval)
      meval=min(neval,nham3)

      call nftim('start of wavefunction generation')
      write(iwave) meval
      call outrow(eig3,meval,iwave)

!        write(idipd)(ham3(i,level3),i=1,nham3)  ! line 9 for J=0 dipole 

        if (ztran) then
           rewind(iwvpb)
           read(iwvpb)((plegw(i,j),i=1,nalf),j=1,nalf)
        end if

        nrad=(npnt*(npnt+1-2*ipar))/2
        levelp=0
        if (zplot) then
           isize=npnt*(npnt+1)*nalf/2
           allocate(wr(npnt),r(npnt))
           allocate(wfs(nplot,isize),wmax(nplot),wf(isize))
           wfs=0.d0
           wmax=0.d0
     call radgrid(zmors1,re1,diss1,we1,npnt,wr,r,idia,xmass)        
            allocate(jxcos0(npth),jwalf0(npth))
            alf = 0.d0
            call jacbasis(jxcos0,jwalf0,npth,alf,alf)

      write(36,1003) npth,alf,(jxcos0(i),jwalf0(i),i=1,npth/2)

1003  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))


      open(unit=iplot,file='wfs/waves.info')

      write(iplot,*)jrot,idia,kpar,iqpar,thresh
      write(iplot,*)zbisc,zperp,zembed,zmors1,zmors2,xmass,g1,g2
      write(iplot,*)re1,diss1,we1,re1,diss1,we1
      write(iplot,*)nmax1
      write(iplot,*)r2
      write(iplot,*)wr**2
      write(iplot,*)nalf
      write(iplot,*)jxcos0
      write(iplot,*)jwalf0
      write(iplot,*)nploti,nplot
      write(iplot,*)(eig3(i),i=nploti,nploti+nplot-1)

           wr=1.d0/wr
        end if

      do 40 level3=1,meval
        index=1
        wvfunc = 0.0d0
        wvfuncdvr = 0.0d0
        index = 1    
       do 30 irs=1,nrad
          do 300 igamma=1,nalf
            if (iv2(2,igamma) .eq. 0) goto 300
            do 20 level2=0,iv2(2,igamma)-1
              ivec2=iv2(1,igamma)+level2
              wvfuncdvr(index)=wvfuncdvr(index)+&
                vecs2d(irs,ivec2)*ham3(ivec2,level3)
20          continue
            index=index+1
300       continue
30      continue
        write(iwave) wvfuncdvr
        if (jrot.eq.0.or.(jrot.eq.1.and.kz.eq.0)) then
        call jtran(wvfunc,wvfuncdvr,plegw,nalf,nrad,nbass)
        write(idipd) wvfuncdvr   ! line 9 for J=0 dipole 
        write(idip) wvfunc       ! line 9 for J=0 dipole 

        if (zplot.and.(level3.ge.nploti)) then

           levelp=levelp+1

      alf=sqrt(dble(jrot*(jrot+1)-kz**2)/2.d0)
      wf=0.d0

      call sumwf(wr,wvfunc,wf,alf,jxcos0,ipar,npnt,nalf,npth,nbass,isize,jwalf0)

      do i=1,isize
         wfs(levelp,i)=wf(i)
      end do

         end if
        end if
40    continue

      if (zplot) then 
         call plotwfs(wfs,nploti,nplot,ithre,isize,wmax) 
            write(iplot,*)wmax
            close(unit=iplot)
         end if

      call nftim('end of wavefunction generation')
      return
      end

!     ***********************************************************************

      subroutine trans(iv1l,iv2l,ndim2l,vecs1l,&
                       vecs2l,vecs3l,phi,psi,evall,nham3)

!     if ztran then this routine transforms the sets of 1d, 2d and 3d
!     coefficients to psi, the wavefunction amplitudes at the dvr points

      implicit real*8 (a-h,o-y), logical (z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
     
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1
      common /split2/ re2,diss2,we2,beta2,ur2,urr2,a2,iu2
      common /mass/ xmass(3),g1,g2,xmassr(3)

      real*8, dimension(max2d,ndima) :: vecs1l
      real*8, dimension(max2d,max2d) :: vecs2l
      real*8, dimension(nham3) :: vecs3l
      dimension iv1l(ndimc,ndimb)
      dimension ndim2l(ndimc)
      real*8, dimension(idvr,npnt1,npnt2) :: psi
      dimension iv2l(ndimc)
      real*8, dimension(neval) :: evall
      real*8, dimension(nham3,ndima,ndimb) :: phi

      rewind iout1
      rewind iout2
      rewind ivecs1
      rewind ivecs2

!     skip header on iout2

      do 10 i=1,5
        read(iout2)
   10 continue

!     read in iv1 and iv2 from a seperate stream (iout1)
      call igetro(iv1l,ndimb*ndimc,iout1)
      call igetro(iv2l,ndimc,iout1)

      do 32 ione = 1,npntc
      nham2 = 0
!     recall the size of ham2 for each npntc
      do 3 itwo = 1,npntb
        nham2 = nham2 + iv1l(ione,itwo)
    3 continue
      ndim2l(ione) = nham2
   32 continue

!     now rewrite the vectors to different streams

      do 35 ione = 1,npntc
      nham2 = ndim2l(ione)
      if (nham2 .gt. 0) then
         do 5 kk=1,npnta
         call getrow(vecs1l(1,kk),nham2,iout2)
         if (iv2l(ione) .ne. 0) call outrow(vecs1l(1,kk),nham2,ivecs1)
    5    continue
      endif
   35 continue
      rewind ivecs1

      do 31 ione = 1,npntc
      nham2 = ndim2l(ione)
      ivm = iv2l(ione)
      if (nham2 .gt. 0) then
         do 30 ind = 1,ivm
         call getrow(vecs2l(1,ind),nham2,iout2)
         call outrow(vecs2l(1,ind),nham2,ivecs2)
   30    continue
      endif
   31 continue
      rewind ivecs2

      read (iout2) meval
      write(iwave) meval

!     call the eigenvalues into evall
      call getrow(evall,meval,iout2)
      call outrow(evall,meval,iwave)

!     now ready to start reading the 3-d vectors

!     1-d vectors:  ivecs1 (stream 26)
!     2-d vectors:  ivecs2 (stream 27)
!     3-d vectors:  iout2  (stream 25)

      ind2 = 0
      do 58 ic = 1,npntc
      nham2 = ndim2l(ic)
      do 56 j  = 1,iv2l(ic)
      ind2 = ind2 + 1
      if (nham2 .gt. 0) call getrow(vecs2l(1,j),nham2,ivecs2)
      do 54 ia = 1,npnta
      if (nham2.gt.0 .and. j.eq.1)&
       call getrow(vecs1l(1,ia),nham2,ivecs1)
      ind1 = 0
      do 52 ib = 1,npntb
      sum1 = 0.0d0
      do 50 k  = 1,iv1l(ic,ib)
      ind1 = ind1 + 1
      sum1 = sum1 + vecs1l(ind1,ia)*vecs2l(ind1,j)
   50 continue
      phi(ind2,ia,ib) = sum1
   52 continue
   54 continue
   56 continue
   58 continue

      do 959 ll=1,meval

      call getrow(vecs3l,nham3,iout2)

!     set psi to zero
       psi = 0.0d0

      if (ztheta) then

        if (.not. zr2r1) then

          ind2 = 0
          do 158 ic = 1,npnt2
          do 156 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 154 ia = 1,idvr
          do 152 ib = 1,npnt1
          psi(ia,ib,ic) = psi(ia,ib,ic) + vecs3l(ind2)*phi(ind2,ia,ib)
  152     continue
  154     continue
  156     continue
  158     continue

        else

          ind2 = 0
          do 258 ic = 1,npnt1
          do 256 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 254 ia = 1,idvr
          do 252 ib = 1,npnt2
          psi(ia,ic,ib) = psi(ia,ic,ib) + vecs3l(ind2)*phi(ind2,ia,ib)
  252     continue
  254     continue
  256     continue
  258     continue

        endif

      else

        if (.not. zr2r1) then

          ind2 = 0
          do 358 ic = 1,idvr
          do 356 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 354 ia = 1,npnt1
          do 352 ib = 1,npnt2
          psi(ic,ia,ib) = psi(ic,ia,ib) + vecs3l(ind2)*phi(ind2,ia,ib)
  352     continue
  354     continue
  356     continue
  358     continue

        else

          ind2 = 0
          do 458 ic = 1,idvr
          do 456 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 454 ia = 1,npnt2
          do 452 ib = 1,npnt1
          psi(ic,ib,ia) = psi(ic,ib,ia) + vecs3l(ind2)*phi(ind2,ia,ib)
  452     continue
  454     continue
  456     continue
  458     continue

        endif

      endif

!     now save the wavefunction for each level (to stream 28)
      call outrow(psi,idvr*npnt1*npnt2,iwave)

  959 continue
      return
      end

!#############################################################################

      subroutine sqout(sqmat,ndim)
!     print lower triangle of square matrix
      real*8 sqmat(ndim,ndim)
      do 30 i=1,ndim
      write(6,1020) (sqmat(i,j),j=1,i)
 1020 format(10f13.7)
   30 continue
      return
      end

!###########################################################################
      subroutine symout(symmat,ndim)

!     print out lower triangle of symmetric matrices                #008

      real*8 symmat(1)
      ip=0
    3 llow=10*ip+1
      lup=min(llow+9,ndim)
      write(6,4) (i,i=llow,lup)
    4 format(/,i11,9i13)
      ind0=llow*(llow+1)/2
      do 5 i=llow,ndim
      itop=min(lup,i)
      ktop=itop-llow+ind0
      write(6,7) i,(symmat(j),j=ind0,ktop)
    7 format(i4,f12.7,9f13.7)
      ind0=ind0+i
    5 continue
      if(lup.ge.ndim) return
      ip=ip+1
      go to 3
      end

      subroutine wrtham(hamil,nham)
!     print hamiltonian matrix                                      #011
      real*8 hamil(nham,nham)
      write(6,1010)
 1010 format(5x,'hamiltonian matrix'/)
      do 30 i=1,nham
      write(6,1020) (hamil(i,j),j=1,i)
 1020 format(10f13.7)
   30 continue
      return
      end

      subroutine getrow(row,nrow,iunit)
!     fast non-formatted read
      real*8 row(nrow)
      read(iunit) row
      return
      end

      subroutine igetro(ivec,nsize,istream)
      dimension ivec(nsize)
      read(istream) ivec
      return
      end

      subroutine outrow(row,nrow,iunit)
!     fast non-formatted write                                      #025
      real*8 row(nrow)
      write(iunit) row
      return
      end

      subroutine ioutro(ivec,nsize,istream)
      dimension ivec(nsize)
      write(istream) ivec
      return
      end
      subroutine nftim(text)
      common/timing/itime0
      character text*(*)
      write(6,10)
      write(6,*) 'Time at ',text,' is.........'
      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/) 
      return
10    format(/)
      end

!########################################################################
      subroutine mkmain(hr,r1m2t,htheta,r,theta,kz,plegw)

!     this routine controls the dvr calculation in the case of
!     symmetrised radau coordinates with z axes perpendicular to the 
!     molecular plane.
!     written by max kostin, 2001.

      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zplot,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,idip,idipd,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,iwvpb,iplot
                                      
      real*8, dimension(npnt,npnt) :: hr
      real*8, dimension(npnt,npnt) :: r1m2t
      real*8, dimension(nalf,nalf) :: htheta
      real*8, dimension(nalf,nalf) :: plegw
      real*8, dimension(npnt) :: r
      real*8, dimension(nalf) :: theta
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig2,eigtmp
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::vecs2d
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eigs2d
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::ham3
      REAL*8, ALLOCATABLE, DIMENSION(:) ::eig3
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv2
      INTEGER, ALLOCATABLE, DIMENSION(:) ::nv2

      data autocm /2.19474624d+05/

      ALLOCATE(ham2(max2d,max2d),eig2(max2d),iv2(2,nalf),&
               vecs2d(max2d,max3d),eigs2d(max3d),nv2(nalf))
      if (.not.zcut) ALLOCATE(eigtmp(nalf*max2d))

      emaxau=emax2/autocm

! no need for the 1d diagonalisations as there is no possiblity of
! truncation as the symmetry would be broken.
      write(6,130)
      call nftim('beginning of 2d loop mkmain')

!      do 30 igamma = nalf,1,-1     
      do 30 igamma = 1,nalf     
        nv2(igamma) = (npnt*(npnt+1-ipar*2))/2      
        nham2 = nv2(igamma)
        call z_blc2d1(theta(igamma),r,hr,ham2,nham2,&
                          term,term2,kz)

!        call diag_dac(ham2,nham2,nham2,eig2)
        call diag_rrr(ham2,nham2,nham2,eig2,zcut,emaxau)

        if (.not. zcut) then
          call choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,&
                      vecs2d,nv2,eigtmp)

          if(igamma .eq. nalf .and. .not. zpmin)&
              write(6,110) (itmp, iv2(2,itmp),itmp=1,nalf)
        else
           call cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)
          if (.not. zpmin) write(6,110) igamma, iv2(2,igamma)

        endif
30    continue

        nbass=max2d*nalf
      if (ztran) then
         write(iwave) nalf,nbass
         write(iwave) (iv2(2,ii),ii=1,nalf)
      endif
!      call testiv(iv2,nbass)

      nham3 = iv2(1,nalf) + iv2(2,nalf) - 1
      if (.not. zpmin) write(6,120) nham3
      call nftim('end of 2d loop')
      write(6,160)

      if (.not.zcut) DEALLOCATE(eigtmp)
      ALLOCATE(ham3(nham3,nham3))

      call bloc3d(htheta,r1m2t,ham3,eigs2d,vecs2d,&
                    iv2,nv2,ham2,nham3,r)

      DEALLOCATE(ham2,eig2,eigs2d,nv2)

      rr=0.d0
      zz=.true.

! writing J=0 output for dipole in stream idip
      if (jrot.eq.0.or.(jrot.eq.1.and.kz.eq.0)) then
     write(idip)ipar,ipar,ipar,npnt1,ipar,jrot,ipar,min(neval,nham3),0,ipar !1
         write(idip)zz,zz,zz,rr,rr,rr,rr,rr,zz                 !2
         write(idip)rr,rr,rr,rr,rr,rr                          !3
         write(idip)                                           !4
         write(idip)r                                          !5
     write(idipd)ipar,ipar,ipar,npnt1,ipar,jrot,ipar,min(neval,nham3),0,ipar !1
         write(idipd)zz,zz,zz,rr,rr,rr,rr,rr,zz                 !2
         write(idipd)rr,rr,rr,rr,rr,rr                          !3
         write(idipd)                                           !4
         write(idipd)r                                          !5
      end if
 
      call nftim('end of 3d ham building')
      write(6,170)
      ALLOCATE(eig3(nham3))

      call diag3d(ham3,nham3,eig3,kz)

      call nftim('end of diagonalising 3d')

      if (ztran) call transr2(iv2,vecs2d,ham3,eig3,nham3,nbass,kz,r)

! writing last line for J=0 output for dipole in stream idip
      if (jrot.eq.0.or.(jrot.eq.1.and.kz.eq.0)) then
         write(idip) theta
         write(idipd) theta
      end if

      DEALLOCATE(iv2,vecs2d,ham3,eig3)
      return
110   format(5x,'for gamma = ',i2,' selected ',i3,' energies.')
120   format(25x,'total = ',i5)
130   format(5x,'starting the 2d calculations.')
160   format(5x,'building the 3d hamiltonian.')
170   format(5x,'diagonalising the 3d hamiltonian.')  
    
      end 
      
!########################################################################     


      subroutine z_blc2d1(xcos,r,hr,ham2,nham2,term,term2,kz)
      implicit real*8(a-h,o-y),logical(z)
      common /size/ npnt1,npnt2,nalf,nmax1,nmax2,maxleg,nalf2,idvr,&
                    npnt,nlim1,nlim2,neval,ncoord,&
                    jrot,kmin,idia,ipar,&
                    max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,&
                    ndima,ndimb,ndimc,iq,emax1,emax2,&
                    nploti,nplotf,nplot,ithre,npth
      common /split1/ re1,diss1,we1,beta1,ur1,urr1,a1,iu1

      real*8, dimension(npnt,npnt) :: hr
      real*8, dimension(npnt) :: r
      real*8, dimension(nham2,nham2) :: ham2
      data xp5/0.5d0/,x1/1.0d0/,x4/4.0d0/
      factr2 = sqrt(xp5)

      ham2 = 0.0d0
      realj = DBLE(jrot)
      realkz = DBLE(kz)

!     allow for j > 0 case
      if (jrot .ne. 0) then
! for q=0(+) term goes with -
! for q=1(-) term goes with +
         if(iq .eq. 0)then
            term= - realj*(realj+x1)/x4
         else
            term= realj*(realj+x1)/x4
         endif
         fact = (realkz**2)/x4
         fact3=0.d0
!         term=1.d0
        if (kz .eq. 1) fact3 = term * (xcos)/(x1-(xcos)**2)
        fact3=0.d0


!-----  Extra NBO term ------
         s=1.0d0-urr1/ur1
         term2=(realj**2+realj-realkz**2)/4.0d0
         fact2=s*term2*x1/(x1-(xcos)**2)
!----------------------------
         ia = 0
         do 15 ibeta=1,npnt
            do 25 ialpha=1,ibeta-ipar
!            walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
!            wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
        walpha = xp5 / (r(ialpha)*r(ialpha)*urr1)
        wbeta  = xp5 / (r(ibeta)*r(ibeta)*urr1)
        wsum   = (walpha + wbeta) * fact
        wsum3  = (walpha + wbeta) * fact3
        wsum2  = (walpha + wbeta) * fact2
        ham2(ialpha+ia,ialpha+ia) = ham2(ialpha+ia,ialpha+ia)+wsum+wsum3+wsum2
25           continue
            ia=ia+ibeta-ipar
!            write(88,*)ia,ibeta,ipar
            if (ipar .eq. 0)ham2(ia,ia) = ham2(ia,ia) + wsum + wsum2 + wsum3
 15     continue
!            stop
      endif

      q=x1
      if(ipar .eq. 1) q=-x1
      iap=0
      do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ibeta .eq. ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+hr(ialphp,ialpha)
              if(ibeta .eq. ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*hr(ibetap,ialpha)
              if(ialpha .eq. ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*hr(ialphp,ibeta)
              if(ialpha .eq. ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+hr(ibetap,ibeta)
              if(ialpha .eq. ialphp .and. ibeta .eq. ibetap) then
                call potv(v,r(ialpha),r(ibeta),xcos)
!                v=0.d0
                ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+v
              endif
              if(ialpha .eq. ibetap .and. ibeta .eq. ialphp) then
                call potv(v,r(ialpha),r(ibeta),xcos)
!                v=0.d0
                ham2(ialphp+iap,ialpha+ia)=&
                  ham2(ialphp+iap,ialpha+ia)+q*v
              endif
              if(ialphp .eq. ibetap) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha .eq. ibeta) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
!            write(88,*)ia,iap
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
!            write(88,*)
!            write(88,*)
!            write(88,*)
      return
      end

!#################################################################

      MODULE constants
      IMPLICIT NONE
      INTEGER, PARAMETER :: real_kind=SELECTED_REAL_KIND(8,40)
      END MODULE constants

      SUBROUTINE jac_basis(nn,nb,alf,bet,x,basis)
      USE constants
      IMPLICIT NONE
      INTEGER :: n,i,nb,nn
      REAL(KIND=real_kind) :: alf, bet
      REAL(KIND=real_kind) :: x(nn),basis(0:nb,nn)
      REAL(KIND=real_kind) :: bass(0:nb,nn),norm(0:nb)

      CALL norms2(norm,nb,alf,bet)
      CALL jacgt(x,bass,alf,bet,nn,nb)

      DO i=0,nb
         DO n=1,nn
            basis(i,n)=bass(i,n)*norm(i)
         END DO
      END DO

      RETURN
      END SUBROUTINE jac_basis

!##################################################################

      SUBROUTINE jacgt(x,bass,alf,bet,nn,nb)
      USE constants
      IMPLICIT NONE
      INTEGER :: n,I,nn,nb
      REAL(KIND=real_kind),EXTERNAL :: gammln
      REAL(KIND=real_kind) :: alf, bet,lmd, x0,x1,x2
      REAL(KIND=real_kind) :: x(nn),bass(0:nb,nn)
      REAL(KIND=real_kind) :: A1n(nb),A2n(nb),A3n(nb),A4n(nb)
      DATA x0,x1,x2/0.0d0,1.0d0,2.0d0/
      lmd=alf+bet+x1
      DO n=1,nb
         A1n(n)=x2*(n+x1)*(n+lmd)*(x2*n+lmd-x1)
         A2n(n)=(x2*n+lmd)*(alf*alf-bet*bet)
         A3n(n)=(x2*n+lmd-x1)*(x2*n+lmd)*(x2*n+lmd+x1)
         A4n(n)=x2*(n+alf)*(n+bet)*(x2*n+lmd+x1)
      END DO   
      bass=x0
      DO 60 I=1,nn
      bass(0,I)=x1
      IF(nb.LT.1) GO TO 70
      bass(1,I)=(alf-bet+(lmd+x1)*x(I))/x2
      DO 80 n=2,nb
         bass(n,I)=((A2n(n-1)+A3n(n-1)*x(I))*bass(n-1,I)&
         &         -A4n(n-1)*bass(n-2,I))/A1n(n-1)
 80   CONTINUE
 70   CONTINUE 
 60   CONTINUE
      RETURN
      END SUBROUTINE jacgt

!##################################################################

      SUBROUTINE gaujac(x,w,n,alf,bet)
      USE constants
      IMPLICIT NONE
      INTEGER :: n
      REAL(KIND=real_kind):: alf,bet,x1,x2,x3
      REAL(KIND=real_kind):: w(n),x(n)
      REAL(KIND=real_kind),PARAMETER :: EPS=3.0D-14
      INTEGER,PARAMETER :: MAXIT=10
      INTEGER :: i,its,j
      REAL(KIND=real_kind)::alfbet,an,bn,r1,r2,r3
      REAL(KIND=real_kind)::c1,c2,c3,p1,p2,p3,pp,temp,z,z1
      REAL(KIND=real_kind),EXTERNAL :: gammln
      DATA x1,x2,x3/1.0d0,2.0d0,3.0d0/
      DO 13 i=1,n
         IF(i==1)THEN
            an=alf/DBLE(n)
            bn=bet/DBLE(n)
            r1=(x1+alf)*(2.78D0/(4.0D0+DBLE(n*n))+0.768D0*an/DBLE(n))
            r2=x1+1.48D0*an+0.96D0*bn+0.452D0*an*an+0.83D0*an*bn
            z =x1-r1/r2
         ELSE IF(i==2)THEN
            r1=(4.1D0+alf)/((x1+alf)*(x1+0.156D0*alf))
            r2=x1+0.06D0*(DBLE(n)-8.0D0)*(1.0D0+0.12D0*alf)/DBLE(n)
            r3=x1+0.012*bet*(x1+0.25D0*DABS(alf))/DBLE(n)
            z=z-(x1-z)*r1*r2*r3
         ELSE IF(i==3)THEN
            r1=(1.67D0+0.28D0*alf)/(x1+0.37D0*alf)
            r2=x1+0.22D0*(DBLE(n)-8.0D0)/DBLE(n)
            r3=x1+8.0D0*bet/((6.28D0+bet)*DBLE(n*n))
            z=z-(x(1)-z)*r1*r2*r3
         ELSE IF(i==n-1)THEN
            r1=(x1+0.235D0*bet)/(0.766D0+0.119D0*bet)
            r2=x1/(x1+0.639D0*(DBLE(n)-4.0D0)/&
            &       (x1+0.71D0*(DBLE(n)-4.0D0)))
            r3=x1/(x1+20.0D0*alf/((7.5D0+alf)*DBLE(n*n)))
            z=z+(z-x(n-3))*r1*r2*r3
         ELSE IF(i==n)THEN
            r1=(x1+0.37D0*bet)/(1.67D0+0.28D0*bet)
            r2=x1/(x1+0.22D0*DBLE(n-8)/DBLE(n))
            r3=x1/(x1+8.0D0*alf/((6.28D0+alf)*DBLE(n*n)))
            z=z+(z-x(n-2))*r1*r2*r3
         ELSE
            z=x3*x(i-1)-x3*x(i-2)+x(i-3)
         ENDIF
         alfbet=alf+bet

         DO 12 its=1,MAXIT
            temp=x2+alfbet
            p1=(alf-bet+temp*z)/x2
            p2=x1
            DO 11 j=2,n
               p3=p2
               p2=p1
               temp=x2*DBLE(j)+alfbet
               c1=x2*DBLE(j)*(DBLE(j)+alfbet)*(temp-x2)
               c2=(temp-x1)*(alf*alf-bet*bet+temp* &
              &               (temp-x2)*z)
               c3=x2*(DBLE(j-1)+alf)*(DBLE(j-1)+bet)*temp
               p1=(c2*p2-c3*p3)/c1
 11         CONTINUE 
             pp=(DBLE(n)*(alf-bet-temp*z)*p1+x2*(DBLE(n)+alf)* &
           &    (DBLE(n)+bet)*p2)/(temp*(x1-z*z))
            z1=z
            z=z1-p1/pp
            IF(ABS(z-z1).LE.EPS) GOTO 1
 12         CONTINUE            
     1   x(i)=z 
         w(i)=DEXP(gammln(alf+DBLE(n))+gammln(bet+DBLE(n))    &
         &    -gammln(DBLE(n+1))-gammln(DBLE(n)+alfbet+x1))&
       &      *temp*x2**alfbet/(pp*p2)

 13    CONTINUE
       RETURN
       END SUBROUTINE gaujac

       FUNCTION gammln(XX)
       USE constants
       IMPLICIT NONE
       INTEGER :: j
       REAL(KIND=real_kind)::GAMMLN,XX
       REAL(KIND=real_kind)::SER,STP,TMP,X,COF(6)
       REAL(KIND=real_kind)::HALF,ONE,FPF
       DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
      &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
       DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
       X=XX-ONE
       TMP=X+FPF
       TMP=(X+HALF)*LOG(TMP)-TMP
       SER=ONE
       DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
!        SER=SER+COF(J)/(X+0.000001d0)
11     CONTINUE
       GAMMLN=TMP+LOG(STP*SER)
       RETURN
       END FUNCTION gammln

!####################################################################

       SUBROUTINE norms2(norm,nn,alf,bet)
       USE constants
       IMPLICIT NONE
       INTEGER :: n,nn
       REAL(KIND=real_kind) :: alf, bet,lmd,x1,x2,norm1
       REAL(KIND=real_kind) :: norm(0:nn)
       REAL(KIND=real_kind) :: a1,a2,a3,a4
       REAL(KIND=real_kind), EXTERNAL :: gammln
       DATA x1,x2/1.0d0,2.0d0/
       lmd=alf+bet+x1
       do n=0,nn
          a1=gammln(DBLE(n+1))
          a2=gammln(DBLE(n)+lmd)
          a3=gammln(DBLE(n)+alf+x1)
          a4=gammln(DBLE(n)+bet+x1)
          norm1=2**(-lmd)*(x2*DBLE(n)+lmd)*exp(a1+a2-a3-a4)
          norm(n)=SQRT(norm1)
       end do
       END SUBROUTINE norms2
!####################################################################

subroutine lagbasis(aa,nlag,dg,wlag,wln,zd,csx)
  implicit none
  !  zd(i,j)=N_{i-1} L^a_{i-1}(x_j) exp{-x_j/2} sqrt{w_j}       
  !  wlag=dsqrt{w}
  !  wln =dln{w}
  !  nlag= number of points
  !  aa=alfa
  !  csx=sum of points
  
  integer :: nlag
  integer :: info,j,i,di, l
  real*8 :: aa,csx,cu2,zc
  real*8 :: xlag(nlag),wlag(nlag),vnor(0:nlag),pl0(0:nlag,nlag)
  real*8 :: dg(nlag),dg1(nlag),zd(nlag,nlag)
  real*8 :: eigen(nlag),work(5*nlag),wln(nlag)
  real*8,external :: gammaln

  do i=0,nlag-1
     di=dble(i)
     dg(i+1)=(2.d0*di+aa+1.d0)
     if (i.ge.1) dg1(i)=-dsqrt((di+aa)*di)
  end do
  
  
  !### then R is diagonalised ############################
  !### checking both zeroes and weights ##################

  CALL  DSTEV('V',nlag,dg,dg1,zd,nlag,work,info)
  
  if (info.ne.0) then
     write(6,*)'Problems with diagonalisation',info
     stop
  endif
  
  cu2=gammaln(aa+1.d0)*0.5d0

  csx=dble(nlag)*(dble(nlag)+aa)

  do i=1,nlag
     wln(i)=dlog(dabs(zd(1,i)))
     wlag(i)=dexp(wln(i))
  end do
  
  do j=1,nlag
     zc=zd(1,j)
     if (zc.lt.0.d0) then
        do i=1,nlag
           zd(i,j)=-zd(i,j)
        end do
     end if
  end do
  
  return
end subroutine lagbasis

!####################################################################

FUNCTION gammaln(xx)
  implicit none
  real*8 :: xx,ser,stp,tmp,x,y,cof(6)
  integer :: j
  real*8 :: gammaln

  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
       -.5395239384953d-5,2.5066282746310005d0/
  
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammaln=tmp+log(stp*ser/x)
  RETURN
END FUNCTION gammaln
!####################################################################
SUBROUTINE LAGPTNEW(ir,Y,R,WT,DZ,npnt,nmax,zmorse,RE,BETA,A,IU)

  !     SUBROUTINE LAGPT GETS INTEGRATION POINTS AND WEIGHTS FOR      #015
  !     NPNT GAUSS LAGUERRE INTEGRATION AND SETS UP BASIS
  !     FUNCTIONS AT THE INTEGRATION POINTS.
  !
  IMPLICIT real*8(A-H,O-Y), LOGICAL (Z)
  dimension Y(NPNT),R(NPNT),WT(NPNT), dz(npnt,npnt),& !wt2(npnt),&
       DNORMNEW(0:NMAX),BASS(0:NMAX,NPNT),wln(npnt)!,f(0:NMAX,NPNT)
  DATA X0,XP5,X1,X2/0.0D0,0.5D0,1.0D0,2.0D0/,TOLER/1.0D-8/
  !
  myid=0
!
!
  IF (ZMORSE) ALF=DBLE(IU)
  IF (.NOT. ZMORSE) ALF = A + XP5
      ALFM1=ALF-X1
!     SET UP INTEGRATION POINTS AND WEIGHTS.
!      CALL LAGUER(NPNT,Y,WT,ALF,B,C,CSX,CSA,TSX,TSA,CC)
     CALL  lagbasis(ALF,npnt,y,wt,wln,dz,tsx)


! please note that wt are the square root of the weights wout the exp factor!!!!

! Note that new weights are defined so that : 
!
!  sum_i w_i = 1 = int 1/(alpha!) rho^alpha exp(-rho) d rho
!
!  but wln_i = 1/2 * log(w_i)
!
!  so sum_i exp(2*wln_i) = 1  (see line 190 below)
!  This is important for patologic cases with very large alpha
!  Also note that the weights are never used in the code.
!  If you need to use them then please check consinstency .

      TSA = X1
      csx=0.d0
      csa=0.d0
      do l=1,npnt
         csx=csx+y(l)
190         csa=csa+dexp(2.d0*wln(l))
      end do
! set up laguerre's normalis 
      do l=0,NPNT-1
         xl=dble(l)
         dd=gammaln(xl+1.d0)-gammaln(xl+alf+1.d0)
         dnormnew(l)=dexp(dd*0.5d0)
      end do
!      csa=csa*(dnormnew(0)**2)

      if (myid .eq. 0) WRITE(6,1000) NPNT,ir
 1000 FORMAT(/,I8,' POINT GAUSS-LAGUERRE INTEGRATION', &
     &       /,5X,'INTEGRATION POINTS',11X,'WEIGHTS',9X,&
     &            'CORRESPONDING R',I1,/)
      DO 60 I=1,NPNT
      IF (ZMORSE) THEN
!         CALCULATE POTENTIAL AT R = RE+BETA(**-1)*LN(A/Y)
          R(I) = RE + LOG(A/Y(I)) / BETA
      ELSE
!         CALCULATE POTENTIAL AT R = SQRT(Y/BETA)
          R(I) = SQRT(Y(I)/BETA)
      ENDIF
!
!     CALCULATE UNNORMALISED LAGUERRE POLYNOMIALS AT Y
!
!     POLYNOMIAL OF ORDER 0
      xep=dexp(-y(I)*0.5d0)
      BASS(0,I) = X1*xep
      IF (NMAX .LT. 1) GOTO 70
!     POLYNOMIAL OF ORDER 1
      AMX = ALF + X1 - Y(I)
      BASS(1,I) = AMX*xep
!     USE RECURRENCE RELATIONSHIPS FOR POLYNOMIALS OF ORDER > 2
!     N * L(N,ALF) = (2*N+ALF-1-X)*L(N-1,ALF) - (N+ALF-1)*L(N-2,ALF)
      EN = X1
      DO 80 N=2,NMAX
      EN = EN + X1
      AMX = AMX + X2
      BASS(N,I) = (AMX * BASS(N-1,I) - (ALFM1+EN) * BASS(N-2,I)) / EN
   80 CONTINUE
   70 CONTINUE
!
      DO 90 N2=0,NMAX
!     NORMALISE POLYNOMIALS
      BASS(N2,I) = BASS(N2,I) * DNORMNEW(N2)
   90 CONTINUE
   60 CONTINUE

 !    do n=0,nmax
 !    do l=1,npnt
 !       f(n,l)=sqrt(beta*y(l)**(ALF+1.d0))*BASS(n,l)
 !    end do
 !    end do
      alf2=0.5d0*(alf+1.d0)
      if (zmorse) then
         do l=1,npnt
            wln(l)=wln(l)-ALF2*log(y(l))
         end do
      else
         do l=1,npnt
            wln(l)=wln(l)-log(sqrt(2.d0))-0.25d0*log(beta)-ALF2*log(y(l))
!            wt(l)=wt(l)/sqrt(sqrt(2.d0)*(beta**0.25d0)*(y(l)**ALF2))
         end do
      end if

     do l=1,npnt
        
      if (myid .eq. 0) WRITE(6,1010) Y(l),WLN(l),R(l)
 1010 FORMAT (F23.15,D25.12,F13.5)
      IF ((R(l) .LT. X0) .and. (myid .eq. 0)) WRITE(6,1015) l
 1015 FORMAT(5X,'***** WARNING: FOR INTEGRATION POINT',I3, &
     &       ', R LESS THAN ZERO *****')
     end do

!     CHECK THAT THE CORRECT POINTS & WEIGHTS HAVE BEEN GENERATED
      if (myid .eq. 0) WRITE(6,1020) CSX,CSA,TSX,TSA
 1020 FORMAT(/5X,'COMPUTED SUM OF POINTS',D26.15,' & WEIGHTS',D26.15,&
     &       /5X,'EXACT    SUM OF POINTS',D26.15,' & WEIGHTS',D26.15)
      IF (ABS((CSX-TSX)/TSX) .GT. TOLER) GOTO 900
      IF (ABS((CSA-TSA)/TSA) .GT. TOLER) GOTO 900
      RETURN
  900 if (myid .eq. 0) WRITE(6,910)
  910 FORMAT(//5X,'(ERROR) POINTS & WEIGHTS IN ERROR, ADJUST ALGORITHM',//)
      STOP
      END
!####################################################################

SUBROUTINE K1K2NEW(XK,HBL,DZ,NPNT,NMAX,NLIM)
  !
  !     SET UP THE TRANSFORMED KINETIC ENERGY INTEGRALS,  T'(HBL) T
  !                                                       ~  ~~~  ~
  !     (NOTE THAT THE RADIAL BASIS FUNCTIONS ARE ALREADY NORMALISED)
  !
  IMPLICIT real*8 (A-H,O-Y), LOGICAL (Z)
  !
  dimension  XK(NPNT,NPNT),HBL(NLIM),DZ(NPNT,NPNT)
  !
  !     IND(I,J) = MAX(I,J) * (MAX(I,J)-1)/2 + MIN(I,J)
  IND(I,J) = MAX(I,J) * (MAX(I,J)+1)/2 + MIN(I,J) + 1
  !
  DO I=1,NPNT
     DO IP=1,I
        XK(I,IP) = 0.0D0
     enddo
  enddo
  
  DO K=1,NPNT
     DO KP=1,K
!        WTKKP = WT(K)*WT(KP)
        DO M=0,NMAX
!           T = BASS(M,K) * WTKKP
           DO MP=0,NMAX
              IN = IND(M,MP)
!              XK(K,KP) = XK(K,KP) + (HBL(IN) * T * BASS(MP,KP))
              XK(K,KP) = XK(K,KP) + (HBL(IN)  * dz(MP+1,KP) * dz(M+1,K))
           enddo
        enddo
        ! this line copies the results onto
        ! the other half of the symmetric matrix
        xk(kp,k)=xk(k,kp)
     enddo
  enddo
  RETURN
END SUBROUTINE K1K2NEW
!#######################################################################
      subroutine jtran(coeffbr,coeffdvr,pleg,nang,nrad,nkbas)
  
      implicit real*8 (a-h,o-y), logical (z)

      REAL*8, DIMENSION(nang,nang) :: pleg
      REAL*8, DIMENSION(nkbas) :: coeffdvr
      REAL*8, DIMENSION(nkbas) :: coeffbr
!      REAL*8, DIMENSION(nang) :: sumk

      data x0/0.0d0/
 
!     transform back to the original fbr-type basis in the
!     associated legendre functions

      ipt=0
      do 50 mn=1,nrad
      do 20 j=1,nang
      sumk=x0
      do 40 k=1,nang
         kk=(mn-1)*nang+k
      sumk=sumk + coeffdvr(kk) * pleg(j,k)
   40 continue
      coeffbr(ipt+j) = sumk
   20 continue
      ipt=ipt+nang
   50 continue

      return
      end
!####################################################################
subroutine radgrid(zmorse,re,diss,we,nr,wt,r,idia,xmass)
implicit none
integer iu,nr,idia
real*8 re,diss,we,a,beta,ur,x4,xp5,x0,x1,xmass(3),g1,g2,b,amtoau
logical zmorse
real*8 r(nr),wt(nr)
real*8, allocatable, dimension(:) :: y
real*8, allocatable, dimension(:,:) :: dz

data x0,xp5,x1,x4/0.0d0,0.5d0,1.0d0,4.0d0/
data amtoau/1.8228883d03/

      if (idia .ge. 1) then
!        scattering coordinates
         g1 = xmass(2) / (xmass(2) + xmass(3))
         g2 = x0
      else
!        radau coordinates
         a = sqrt(xmass(3) / (xmass(1)+xmass(2)+xmass(3)))
         b = xmass(2) / (xmass(1)+xmass(2))
         g1 = x1 - a / (a+b-a*b)
         g2 = x1 - a / (x1-b+a*b)
      endif
 
      write(6,1000) xmass
 1000 FORMAT(/5X,'Vibrational nuclear mass in AMU:',3F12.6)
!     compute the effective moments of inertia

      ur = amtoau/(g2*g2/xmass(1)+x1/xmass(2)+(x1-g2)**2/xmass(3))

zmorse=.false.
      if (zmorse) then
          write(6,1020) re,diss,we
 1020 format(/5x,'Morse function parameters for r basis',&
             /5x,'r equilibrium =',f8.4,' bohr, dissociation energy',&
          d15.7,' hartree &  vibrational frequency =',d15.7,' hartree')
         beta = we * sqrt(xp5*ur/diss)
         a = x4 * diss / we
         iu = int(a+xp5)
         write(6,1030) ur,beta,a,iu
 1030 format(/5x,'Constants used to construct morse oscillators:',&
             /5x,'reduced mass =',d16.7,' a.u., beta =',f8.4,&
                  ' (1/bohr), a =',d16.7,' and u =',i5)
      else
          a=diss
          beta = sqrt(we * ur)
          write(6,1039) a,we,ur,beta
 1039 format(/5x,'Spherical oscillator parameters for r basis:',&
             /5x,'alpha =',f10.5,&
                 ' &  vibrational frequency =',d15.7,' hartree',&
            //5x,'Constants used to construct spherical oscillators:',&
             /5x,'reduced mass =',d16.7,' a.u., beta =',f12.6,&
                  ' bohr**-2')
      endif

allocate(y(nr),dz(nr,nr))

CALL LAGPTNEW(1,Y,R,WT,DZ,nr,nr-1,zmorse,re,beta,a,iu)
deallocate(y,dz)
  return
end subroutine radgrid
!#######################################################################
subroutine plotwfs(wfs,ni,nplot,th0,n3d,wmax)
implicit none 
integer ni,nplot,nst,n3,n2,n1,n0,n3d,i,l,th0,pp,k!,nr,nt
real*8 x,xmax,xmaxw,th,uth,wfs(nplot,n3d),wf(n3d),wmax(nplot)
real*8 sum,ws
character s1a*30,f2*100,state*4

s1a='state'

do i=1,nplot

nst=ni+i-1
n3=nst/1000
n2=(nst-1000*n3)/100
n1=(nst-100*n2-1000*n3)/10
n0=(nst-100*n2-10*n1-1000*n3)
state=char(48+n3)//char(48+n2)//char(48+n1)//char(48+n0)

do l=1,n3d
wf(l)=wfs(i,l)
end do

write(6,*)'====================================='

sum=0.d0
do l=1,n3d
   sum=sum+wf(l) !+wr(n1)*wr(n1)*wr(n2)*wr(n2)*wt(k)*
end do

write(6,*)' Norm (original) .... ',sum

xmax=0.d0
do l=1,n3d
xmax=max(xmax,wf(l))
end do

uth=1.d0/xmax
do l=1,n3d
   wf(l)=wf(l)*uth
end do

th=(10.d0**(th0))
uth=1.d0/th

l=0
do l=1,n3d
   x=wf(l)
   wf(l)=dble(int(wf(l)*uth))*th
end do

x=0.d0
sum=0.d0
do l=1,n3d
   x=max(x,wf(l))
   sum=sum+wf(l)
end do
wmax(i)=x

write(6,*)' Norm (reduced) .... ',sum*xmax
f2='wfs/'//trim(s1a)//'_'//state//'.pb'
write(6,*)trim(f2)
open(unit=21,file=trim(f2))!,form='unformatted')
write(21,*)xmax
write(21,*)wf
close(unit=21)

end do


return
end subroutine plotwfs
!#######################################################################
subroutine sumwf(wr,wf,wfs,alf,jcos,s,nr,nt,nt2,n3dl,n3d,jw)
implicit none
integer nt,nr,l,n1,n2dl,n3dl,i,k2,i0,k1,nst,s,n3,n2,n,n3d,nt2
real*8 jcos(nt2),pleg(0:nt-1,nt2),jx(nt2),wr(nr),jw(nt2)
real*8 sumw,alf,sum
real*8 wf(n3dl),wfs(n3d)

CALL jac_basis(nt2,nt-1,alf,alf,jcos,pleg)

do l=1,nt2
jx(l)=dsqrt(jw(l)*(1.d0-jcos(l)**2)**alf)
end do

do n1=0,nt-1
do l=1,nt2
pleg(n1,l)=pleg(n1,l)*jx(l)!*dsqrt(jwalf(l))
end do
end do

n2dl=nr*(nr+(-1)**s)/2

wfs=0.d0

i0=0
sum=0.d0
do n1=1,nr
do n2=1,n1-s
   n3=n1*(n1-1)/2+n2-1
   do k2=1,nt2

   sumw=0.0
do k1=1,nt
   sumw=sumw+wf(i0+k1)*pleg(k1-1,k2)
end do !k1

      n=n3*nt2+k2
      wfs(n)=sumw**2!*wr(n1)*wr(n2)
   end do !k2
      i0=i0+nt
end do !n2
end do !n1

return
end subroutine sumwf
!########################################################################

subroutine jacbasis(dg,w,n,alf,bet)

  implicit none

  integer :: n
  integer :: info,j,i,di
  real*8 :: xi,alf,bet,a1,a2,a3,a4,ua3,zc,termd,termu,terml
  real*8 :: dg(n),dg1(n),zd(n,n),w(n),vnor(0:n)
  real*8 :: eigen(n),work(5*n)
  real*8,external :: gammaln

      call norms2(vnor,n,alf,bet)
 
  do i=0,n-1
     xi=dble(i)
     a1=2.d0*(xi+1.d0)*(xi+alf+bet+1.d0)*(2.d0*xi+alf+bet)*vnor(i)/vnor(i+1)
     a2=-(2.d0*xi+alf+bet+1.d0)*(alf**2-bet**2)*vnor(i)*vnor(i)
     a3=(2.d0*xi+alf+bet+2.d0)*(2.d0*xi+alf+bet+1.d0)*(2.d0*xi+alf+bet)
     ua3=1.d0/a3
     termu=a1*ua3
     termd=a2/ua3
     dg(i+1)=termd
     if (i.ge.1) then 
     a4=2.d0*(xi+alf)*(xi+bet)*(2.d0*xi+alf+bet+2.d0)*vnor(i)/vnor(i-1)
     terml=a4*ua3
     dg1(i)=terml
     end if
  end do
  
  !### then R is diagonalised ############################
  !### checking both zeroes and weights ##################

  CALL  DSTEV('V',n,dg,dg1,zd,n,work,info)

  if (info.ne.0) then
     write(6,*)'Problems with diagonalisation',info
     stop
  endif
  
  do i=1,n
     w(i)=zd(1,i)**2/vnor(0)**2
  end do
  
  do j=1,n
     zc=zd(1,j)/w(j)
     if (zc.lt.0.d0) then
        do i=1,n
           zd(i,j)=-zd(i,j)
        end do
     end if
  end do
  
  return
end subroutine jacbasis

!####################################################################

