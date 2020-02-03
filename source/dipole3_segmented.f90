      program shallot
      call dipole3b
      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccc
!                                                **001
      subroutine dipole3b
!
!    dipoleb3, and the various subroutines it calls, calculate the
!    transitions between vib-rot eigenfunctions due to the
!    molecular dipole surface.
!
!    transition moments in the molecule fixed z and x axes are
!    calculated separately and then combined to give the line
!    strength, s(f-i).
!
!    the program is based on theory of
!    miller, tennyson and sutcliffe, mol. phys. 66, 429 (1989)
!    updated for dvr wavefunctions and for symmetrised radau coordinates
!
!    dipole3b is drive by programs
!    dvr3drjz and rotlev3b which are extension of the versions in
!    J. Tennyson, J.R. Henderson and N.G. Fulton,
!    Comput. Physomms., 6, 175-198 (1995)
!
!    it has the option to produce output files for program
!    spectra to calculate simulated spectra at a given temperature.
!
!    there are three sources of input for the program:
!
!    a) the dipole surface has to be supplied by the user in
!         subroutine dipd.
!
!    b) the user must also supply the following lines of input:
!
!     l1) namelist /prt/
!         this supplies a number of logical control parameters,
!         whose default values (in brackets) are set by block data
!        b) zprint(f) = .true. supplies extra print out for
!                       debugging purposes.
!        c) zpmin (f) = .true. supplies less  print out for
!                       large runs.
!        d) ztra  (t) = .true. writes out the data needed for
!                       program spectra to calculate simulated
!                       spectra.
!        e) zstart(f) = .true. if we are writing out for spectra
!                       for the first time.
!
!         the namelist can also carry values for the input and
!         output streams, (defaults in brackets):
!        a) iket (11) = input stream for the ket.
!        b) ibra (12) = input stream for the bra.
!        c) itra (13) = output stream to program spectrm (if ztra).
!                       note that for all times other than the
!                       dipole assumes that we have accessed the
!                       permanent dataset or file which has the
!                       data from previous runs and that we are
!                       writing to the end of that file.
!                   ************************************************
!                   **  for the sake of safety you are therefore  **
!                   **  advised to keep one previous edition as   **
!                   **  backup!                                   **
!                   ************************************************
!        d) iscr (24) = scratch file used for restart runs
!     also:
!        a) ires (0) restart parameter
!           ires = 0, normal run
!           ires = 1, restart run
!        b) nblock (1000) number of k --> k' blocks to be attempted
!
!     l2) title  (9a8 format)
!
!     l3) npot, nv1, nv2 (all in i5 format)
!        a) npot  = number of Gauss-Legendre integration points used
!        b) nv1   = number of ket eigenfunctions considered.
!                   if this is input as zero, all available
!                   ket eigenfunctions will be considered when
!                   computing transitions.
!        c) nv2   = as above for the bra.
!
!    l4) the equilibrium position of the molecule is required for
!        rotationless transitions only (not read otherwise), given
!        as r1, r2 (bohr radii) and cos theta in the chosen
!        co-ordinate system.
!
!    c) the bulk of the data required for this program is supplied
!       from programs dvr3d or rotlev3 as follows:
!       (n.b. this data is supplied for both ket and bra)
!
!     l1) mbass, idia, ipar, nlim1, nlim2, j, kmin, neval, zembed,
!         zmors1, zmors2, jrot, ncoord, nmax1, nmax2, lmax
!        a) all of these parameters have the same meaning as in
!           program triatom, to which the user is referred, except:
!        b) j is the absolute value of jrot.
!        c) mbass is the total number of basis functions.
!        d) neval is the number of vib-rot eigenfunctions supplied.
!         the program compares parameters for consistency and takes
!         the maximum value of size parameters where appropriate.
!
!     l3) g1, g2, xmass, re1, diss1, we1, re2, diss2, we2
!         parameters governing the co-ordinate system, the masses,
!         and the morse functions for the radial co-ordinates.
!
!     l4) eval = a one dimensional array holding neval eigenvalues
!                for the vib-rot eigenfunctions.
!
!     l5 - neval+5) d! = a two-dimensional array holding the vectors
!                for the vib-rot eigenfunctions as coefficients
!                relating to the basis function labels supplied in
!                line 2 above.
!
!    input and output on streams iket, ibra and itra are in atomic
!    units in order to keep consistency with other programs in the
!    suite. data printed out at the end of dipole, however, is
!    given in wavenumbers and debye.
!
!    ******* an important note on selection rules. ****************
!
!    because the dipole operator is anti-symmetric, allowed
!    transitions mix wavefunctions of different overall parity.
!    this means that transitions from e (p = 0) to f (p = 1) and
!    vice versa can only occur if j' = j", and conversely if
!    j' = j" +/- 1 only e to e and f to f transitions occur.
!
!    a further selection rule occurs for triatomics which include
!    a homonuclear diatomic (idia = +/-2)
!
!    for idia = 2 (scattering or jacobi coordinates)
!
!    if the wavefunctions have been generated with the molecular
!    z-axis fixed along r2 (zembed= .true.)
!    then j even/odd to j odd/even transitions cannot occur.
!
!    if, however, zembed= .false. (z-axis parallel to r1) then only
!    j even to j odd (and vice versa) transitions occur.
!
!    for idia = -2 (radau coordinates)
!
!    for j' = j" then ipar1 must equal ipar2
!
!    for |j' - j"| = 1, then ipar1 and ipar2 must be different.
!
      implicit double precision (a-h,o-y), logical (z)
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      namelist/prt/ zprint, zpmin, ztra, zstart,&
                    iket, ibra, itra, iscr, ires, nblock
      common /head/ title
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      common/timing/itime0
      character(len=8) title(9)

      write(6,200)
200   format(//,5x,'Program DIPOLE3 (version of January 2004):',/)
      call SYSTEM_CLOCK(itime0,irate2,imax2)
!
!     read in of the logical control parameters zprint, ztra, etc
!
!  *  zprint= .true. gives extra print out - useful for debugging.
!  *  zpmin = .true. gives less  print out - useful for big runs.
!  *  ztra= .true. outputs the information necessary for program
!     spectrum to compute a simulated spectrum.
!  *  the default values of all three parameters are .false.
!
      read(5,prt)
      read(5,100) title
100   format(9a8)
      write(6,202) title
202   format(5x,9a8)

!     read control parameters

      call insize


      call main

      stop
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **002
      block data
      implicit double precision (a-h,o-y), logical (z)
!
!     block data stores the default values of the logical control
!     parameters and the ground zero energy.
!     it also stores the default values of the ket input stream, iket,
!     the bra input stream, ibra, and the output stream for program
!     spectrum, itra.
!
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      data zmors1/.true./, zprint/.false./, ztra/.true./,&
           zmors2/.true./, zpmin /.false./, ires/0/, nblock/1000/,&
           zstart/.false./, iket/11/, ibra/12/, itra/13/, iscr/24/
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **003
      subroutine insize
!
!     subroutine insize reads in the parameters which control the
!     overall size of the problem. some of these are passed from the
!     output streams of programs triatom and rotlev via the input
!     streams iket and ibra.
!     this means that many parameters are duplicated, giving the program
!     the opportunity to check that the bra and ket are compatible.
!     integration parameters for the two radial co-ordinates and the
!     angular functions must be set by the user, and are inputted on
!     stream 5 (five).
!
      implicit double precision (a-h,o-y), logical (z)
!
!     the following size parameters have these functions:
!
!     npnt1: number of gauss-laguerre dvr points for r1
!     npnt2: number of gauss-laguerre dvr points for r2
!     lmax1: max number of angular basis functions in ket
!     lmax2: max number of angular basis functions in bra
!     npot : number of gauss-legendre integration points for dip
!     npnt : max(npnt1,npnt2)
!     nbin : largest binomial coef. required for angular integration(+1)
!     idia : = 1 scattering coordinates heteronuclear diatomic
!            = 2 scattering coordinates homonuclear diatomic
!            = -1 radau coordinates hetronuclear diatomic
!            = -2 radau coordinates homonuclear diatomic
!     ipar : parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
!          : taken from ipar1 and ipar2
!     kmin : kmin=1 for sym. rotational basis, =0 for anti-sym.
!          : for non-coriolis calculations, kmin= k.
!     neval: number of eigenvalues supplied from rotlev or triatom
!     ncoord: number of vibrational coordinates explicitly considered
!     if (ncoord .ne. 3) some of the above are dummies, see below.
!
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
!     save masses, g's and embedding in case they are needed
!     in the dipole routine
      common /mass/ xmass(3),g1,g2,zembed,zbisc

      double precision, dimension(3) :: xm1,xm2

      data toler/1.0d-3/
!
!     read in control parameters of problem:
!     ncoord = 2: atom-diatom problem with diatom rigid
!     ncoord = 3: full 3-d triatomic problem (default value)
      ncoord=3
      write(6,1057) iket,ibra,itra
 1057 format(//5x,'input ket wavefunctions read from stream ',&
             'iket  =',i4,&
             /5x,'input bra wavefunctions read from stream ',&
             'ibra  =',i4,&
             /5x,'output transition data written to stream ',&
             'itra  =',i4)
      if (ires.le.0) write(6,1005) iscr
 1005 format(  5x,'restart data written to           stream ',&
             'iscr  =',i4)
      if (ires.gt.0) write(6,1006) iscr
 1006 format(  5x,'RESTART RUN: input data read from stream ',&
             'iscr  =',i4)

!     read in of data for the ket

      open(unit=iket,form='unformatted',recordtype='segmented')
      read(iket) idia1,ipar1,lmax1,nmax11,nmax12,jrot1,kmin1,neval1
      read(iket) zemb1,zmor11,zmor12,xm1,g11,g21,znco1

!     read in of data for the bra

      open(unit=ibra,form='unformatted',recordtype='segmented')
      read(ibra) idia2,ipar2,lmax2,nmax21,nmax22,jrot2,kmin2,neval2
      read(ibra) zemb2,zmor21,zmor22,xm2,g12,g22,znco2

      open(unit=iscr,form='unformatted',recordtype='segmented')
!      if (ztra)  open(unit=itra,form='unformatted')
      if (ztra)  open(unit=itra,form='unformatted',recordtype='segmented')

!     check the bra and ket are consistent

      if (idia1.ne.idia2) then
          write(6,998) idia1,idia2
998       format(//,5x,'** fatal ** diatomic mismatch',/&
                      5x,'idia1=',i2,'  idia2=',i2,/)
         stop
      else
         idia= idia1

         write(6,888) idia
888      format(/,5x,'diatomic parameter idia  =',i4)
      endif
!cccccccccccccccccccc
      if (jrot1.eq.0.or.jrot2.eq.0) then
         zembed= zemb1
         if (jrot1.eq.0) zembed= zemb2
         if (jrot1.eq.0.and.jrot2.eq.0) then
           write(6,980)
980        format(//,5x,'j = 0 -> 0 not allowed: stop')
           stop
         endif
      else
         if (idia  .gt. -2) then
         if (zemb1.neqv.zemb2) then
             write(6,996)
996          format(/,/,5x,'** fatal ** embedding mismatch',/)
             stop
         else
            zembed= zemb1
         endif
         endif
      endif
!ccccccccccccccccc
      if (idia .gt. -2) then
         if (zembed) then
            if (ipar1.ne.ipar2 .and. idia.eq.2) then
               write(6,997) ipar1,ipar2
997          format(//,5x,'** fatal ** parity mismatch, spin forbidden',&
                   /,5x,'ipar1=',i2,'  ipar2=',i2,/)
               stop
            endif
         else
            if (ipar1.eq.ipar2 .and. idia.eq.2) then
               write(6,997) ipar1,ipar2
               stop
            endif
         endif
         zbisc=.false.
      else
         zbisc=.true.
         zembed=.true.
         if (jrot1 .eq. jrot2) then
            if (ipar1.ne.ipar2) then
               write(6,997) ipar1,ipar2
               stop
            endif
         else
            if (ipar1.eq.ipar2) then
               write(6,997) ipar1,ipar2
               stop
            endif
         endif
      endif
      ipar= ipar1
!cccccccccccccccccccccccccccccccccc
      if (zmor11.neqv.zmor21) then
         write(6,995) 1
995      format(//,5x,'** fatal ** r',i1,' radial function mismatch',/)
         stop
      else
         zmors1= zmor11
      endif
!ccccccccccccccccccccccccccccc
      if (zmor12.neqv.zmor22) then
          write(6,995) 2
          stop
      else
         zmors2= zmor22
      endif
!cccccccccccccccccccccccccccccccccccccc
      write(6,884) ncoord
884   format(/,5x,'number of co-ordinates   =',i4)
!
      if (zbisc) then
         write(6,1330)
 1330    format(/,5x,'z axis embedded along the biscetor of r1 and r2')
      else
         if (zembed) then
           write(6,886) 2
886        format(/,5x,'z axis embedded along r',i1,' co-ordinate'/)
         else
           write(6,886) 1
         endif
      endif
!cccccccccccccccccccccccccccccccc
      if (ncoord .gt. 2) then
         if (zmors1) then
            write(6,885) 1
         else
            write(6,775) 1
         endif
         if (nmax11 .ne. nmax21) then
            write(6,875) 1,nmax11,nmax21
            stop
         else
            npnt1=nmax11
         endif
      else
         npnt1=1
      endif
!cccccccccccccccccccccccccccccccccc
      if (zmors2) then
         write(6,885) 2
      else
         write(6,775) 2
      endif
!cccccccccccccccccccccccccccccccccc
      if (nmax12 .ne. nmax22) then
          write(6,875) 2,nmax12,nmax22
          stop
      else
          npnt2=nmax12
      endif
885   format(5x,'morse functions in r',i1,' radial basis set')
775   format(5x,'spherical functions in r',i1,' radial basis set')
875   format(//,5x,'** fatal ** r2 radial function mismatch',&
             /,5x,i5,' dvr points in bra,',i5,' in ket',/)
!cccccccccccccccccccccccccccccccccccc
      if (zpmin) write(6,2255)
 2255 format(/5x,' Minimum printing requested')
!cccccccccccccccccccccccccccccccccccc
!     check parameters are consistent within limit of toler
      if (abs(g12-g11).gt.toler) then
         write(6,919) g11,g12
919      format(/,5x,'co-ordinate system incompatible',&
               /,2(e18.8,5x),/)
         stop
      else
         g1= g11
      endif
!cccccccccccccccccccccccccccccccccccc
      if (abs(g22-g21).gt.toler) then
         write(6,919) g21,g22
         stop
      else
         g2= g22
      endif
!cccccccccccccccccccccccccccccccccccc
      do 1 i=1,3
      if (abs(xm2(i)-xm1(i)).gt.toler) then
         write(6,918) (xm1(j),xm2(j),j=1,3)
918      format(/,5x,'masses incompatible',&
               /,(2(e18.8,5x)))
         stop
      else
        xmass(i)= xm1(i)
      endif
!cccccccccccccccccccccccccccccccccccc
1     continue
!
      read(ibra) re11,diss11,we11,re12,diss12,we12
      read(iket) re21,diss21,we21,re22,diss22,we22

      if (abs(re21-re11).gt.toler) then
         write(6,917) 1,re11,re21
917      format(/,5x,'re',i1,' parameters incompatible',&
               /,2(e18.8,5x),/)
         stop
      endif
!ccccccccccccccccccccccccccc
      if (abs(diss21-diss11).gt.toler) then
         write(6,916) 1,diss11,diss21
916      format(/,5x,'r',i1,' dissociation energy incompatible',&
               /,2(e18.8,5x),/)
         stop
      endif
!ccccccccccccccccccccccccccc
      if (abs(we21-we11).gt.toler) then
         write(6,915) 1,we11,we21
915      format(/,5x,'r',i1,' morse frequency incompatible',&
               /,2(e18.8,5x),/)
         stop
      endif
!ccccccccccccccccccccccccccc
      if (abs(re22-re12).gt.toler) then
         write(6,917) 2,re12,re22
         stop
      endif
!ccccccccccccccccccccccccccc
      if (abs(diss22-diss12).gt.toler) then
         write(6,916) 2,diss12,diss22
         stop
      endif
      if (abs(we22-we12).gt.toler) then
         write(6,915) 2,we12,we22
         stop
      endif
!cccccccccccccccccccccccccccc
!     are we doing a non coriolis problem?
!     set up jk for the two cases and total basis set sizes
!
!     first correct for the case where j=1f has been done non-coriolis
!     but a full calculation is required
      if (znco1.and..not.znco2) then
         if (abs(jrot1) .gt. 1) then
            write(6,889)
            stop
         else
            kmin1=0
         endif
      endif
!ccccccccccccccccccccccccccccc
      if (znco2.and..not.znco1) then
         if (abs(jrot2) .gt. 1) then
            write(6,889)
            stop
         else
            kmin2=0
         endif
     endif
!ccccccccccccccccccccccccccccc
889   format(//,5x,'attempt to mix no coriolis and coupled calcs: stop')
      if (znco1 .and. znco2) then
         jk1= 1
         jk2= 1
         jrot1=abs(jrot1)
         jrot2=abs(jrot2)
         if (abs(kmin2-kmin1).gt.1) then
            write(6,878)
878         format(//,5x,'k levels differ by more than 1',/)
            stop
         endif
         mblock=1
      else
         if (jrot1.eq.0) kmin1= 1
         jk1= jrot1 + kmin1
         if (jrot2.eq.0) kmin2= 1
         jk2= jrot2 + kmin2
         mblock=jk1+jk2+min(jk1,jk2)-2
      endif
!ccccccccccccccccccccccccccc
      nblock=min(nblock,mblock)
      if (nblock.lt.mblock) write(6,1015) nblock,mblock
 1015 format(/i7,' blocks to be calculated out a maximum of',i4)
!
      if (idia .gt. -2) then
         nrade=npnt1*npnt2
         nrado=nrade
         mbass1=nrade*jk1*lmax1
         mbass2=nrade*jk2*lmax2
      else
         nrade=npnt1*(npnt1+1)/2
         nrado=npnt1*(npnt1-1)/2
         jt=jk1/2
!
         mbass1=(nrade+nrado)*jt
!
         if (2*jt .ne. jk1) then
           if (ipar1 .eq. 0) mbass1=mbass1+nrade
           if (ipar1 .ne. 0) mbass1=mbass1+nrado
         endif
!
         mbass1=mbass1*lmax1
!
         jt=jk2/2
         mbass2=(nrade+nrado)*jt
!
         if (2*jt .ne. jk2) then
           if (ipar2 .eq. 0) mbass2=mbass2+nrade
           if (ipar2 .ne. 0) mbass2=mbass2+nrado
         endif
         mbass2=mbass2*lmax2
       endif
!
!     set other parameters to maximum values
!
      mbass= max(mbass1,mbass2)
      npnt=max(npnt1,npnt2)
!
!     read in of integration parameters and states to be considered
!     zeros will give the default values
!
      read(5,101) npot,nv1,nv2,ibase1,ibase2
101   format(5i5)
       if (ibase1 .le. 0 .or. ibase1 .gt. neval1) ibase1 = 1
       if (ibase2 .le. 0 .or. ibase2 .gt. neval2) ibase2 = 1
!
!     write out data from triatom/rotlev runs
!
      write(6,299)
299   format(/,9x,'parameters passed to dipole for the ket & the bra',/)
      write(6,200) mbass1,mbass2,nmax11,nmax21,nmax12,nmax22,&
                  lmax1,lmax2,jrot1,jrot2,ipar1,ipar2,&
                  kmin1,kmin2,neval1,neval2,nv1,nv2,&
                  ibase1,ibase2
200   format(5x,'total number of basis functions    = ',i6,3x,i6,/&
            5x,'number of r1 radial dvr points     = ',i6,3x,i6,/&
            5x,'number of r2 radial dvr points     = ',i6,3x,i6,/&
            5x,'number of angular functions        = ',i6,3x,i6,/&
            5x,'j, total angular momentum          = ',i6,3x,i6,/&
            5x,'ipar, vibrational parity           = ',i6,3x,i6,/&
            5x,'kmin= 1-p, where j + p is parity,  = ',i6,3x,i6,/&
            5x,'total vib-rot functions available  = ',i6,3x,i6,/&
            5x,'vib-rot functions taken (0 = all)  = ',i6,3x,i6,/&
            5x,'neva1 and neval2 used (all)        = ',i6,3x,i6,/)
!
!     are we doing an allowed transition?
!     stop otherwise.
!
      itot= jrot1 + jrot2 + 1 + kmin1 + kmin2
      if (mod(itot,2) .ne. 0 .or. abs(jrot1 - jrot2) .gt. 1) then
         write(6,9999)
9999     format(/,/,5x,'selection rules violated',/)
         stop
      endif
!
!     reset number of vib-rot functions to be considered
!
      if (nv1.gt.0) neval1=min(nv1,neval1-ibase1+1)
      if (nv1.le.0) neval1=neval1-ibase1+1
      if (nv2.gt.0) neval2= min(nv2,neval2-ibase2+1)
      if (nv2.le.0) neval2=neval2-ibase2+1
      jdia=max(1,idia)
      nbin=jrot1+jrot2+3
!
!     check dimension of angular integration: should be even
!
      if (mod(npot,2).ne.0) npot=npot+1
      nn2= npot/2
! set the number of anglular integartion points needed for the problem
      if (idia.lt.2) then
! radau coordiantes and scattering coordinates without symmetry (hetronuclear)
         ipot = npot
      else
! scattering coordiantes with symmetry (homunuclear)
         ipot = npot / 2
      endif
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **004
      subroutine main
!
      implicit double precision (a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      integer, allocatable, dimension(:):: nbass1,nbass2
      allocate(nbass1(jk1))
      allocate(nbass2(jk2))

       write(*,*) 'Iam in main. size(nbass1)=', size(nbass1)
       write(*,*) 'Iam in main. size(nbass2)=', size(nbass2)
       write(*,*) 'jk1, jk2', jk1, jk2


!     generate the subindex arrays needed for trans

      call genind(nbass1,mbass1,jk1,nbmax1,iket)
      call genind(nbass2,mbass2,jk2,nbmax2,ibra)

      write(*,*) 'checkpoint 0: about to call dmain'
      write(*,*) 'Iam in main. size(nbass1)=', size(nbass1)
      write(*,*) 'Iam in main. size(nbass2)=', size(nbass2)

      call dmain(nbass1,nbass2)
      deallocate(nbass1)
      deallocate(nbass2)

      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **005
      subroutine genind(nbass,mbass,jk,nbmax,ivec)
!
!     this subroutine sets up the sub-index arrays needed for trans
!     and reads in the basis function labels for the vectors.
!     it calculates nbmax, the largest value of nbass, neeeded to
!     dimension the space needed for the d-coefficients.
!
      implicit logical (z)
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2

!       dimension nbass(jk),lmin(jk),lbass(jk)
      integer, dimension(jk) :: nbass
      integer, allocatable, dimension(:) :: lmin, lbass
      allocate( lmin(jk) )
      allocate( lbass(jk) )

      if (zprint) write(6,205) jk,mbass
205   format(//,'   j + kmin =',i3,'   mbass=',i7,/)
!
!     read in the basis function labels
!
      mbass0=-999
      write(*,*) 'Iam in genind. size(nbass)=', size(nbass)
      nbass = -999
      read(ivec, err=333) mbass0,lmin,lbass,nbass
      if (mbass0.gt.mbass) goto 999
!
!     generate the sub-index arrays and  find nbmax
!
      nbmax=nbass(1)
      do 2 k=2,jk
      nbmax=max(nbmax,nbass(k))
2     continue
      if (zprint) then
         write(6,201)
201      format(//,5x,'indices generated by genind',/)
         write(6,202) nbmax
202      format(5x,'nbmax= ',i5,/,5x,'nbass follows',/)
         write(6,*) (nbass(k), k=1,jk)
      endif
      deallocate(lmin)
      deallocate(lbass)
      return
999   write(6,200) ivec,mbass,mbass0
200   format(//,5x,'basis function dimensions in error',&
             /,5x,'unit =',i3,' mbass =',i7,' expected, =',i7,' found')
333   write(6,*)'mbass0 =',mbass0,' not consistent with nbass =',nbass
      stop
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     subroutine dmain(nbass1,nbass2)
!     Effective main program.
!     All data etc for the ket are labelled 1;
!     all data etc for the bra are labelled 2.

      implicit double precision (a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      common /mass/ xmass(3),g1,g2,zembed,zbisc

      !dimension nbass1(jk1)
      !dimension nbass2(jk2)
      integer, dimension(:) :: nbass1(jk1)
      integer, dimension(:) :: nbass2(jk2)
      double precision, allocatable, dimension(:)     :: e1
      double precision, allocatable, dimension(:)     :: e2,xe2
      double precision, allocatable, dimension(:,:)  :: binom
      double precision, allocatable, dimension(:,:) :: dipol
      double precision, allocatable, dimension(:)      :: r1
      double precision, allocatable, dimension(:)      :: r2
      double precision, allocatable, dimension(:)       :: xd ,wtd
      double precision, allocatable, dimension(:,:)  :: tz,tx,sint
      double precision, allocatable, dimension(:) :: dstemp, dc1, dc2, dlower, dmiddle, dupper

      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
      write(*,*) 'checkpoint 1dmain: Allocating'
      write(*,*) 'Iam in dmain. size(nbass1)=', size(nbass1)
      write(*,*) 'Iam in dmain. size(nbass2)=', size(nbass2)

      allocate( e1(neval1)  )
      allocate( e2(neval2)  )
      allocate( xe2(neval2) )
      allocate( binom(nbin,nbin)  )
      allocate( dipol(nrade,ipot) )
      allocate( r1(npnt1) )
      allocate( r2(npnt2) )
      allocate( xd(npot)  )
      allocate( wtd(npot) )
      allocate( tz(neval1,neval2) )
      allocate( tx(neval1,neval2) )
      allocate( sint(neval1,neval2) )

      allocate(dstemp(max(nbmax1*neval1,nbmax2*neval2)))
      allocate(dc1(neval1*max(nrade*ipot,nbmax1)))
      allocate(dc2(neval2*max(nrade*ipot,nbmax2)))
      allocate(dlower(neval2*max(nrade*ipot,nbmax2)))
      allocate(dmiddle(neval2*max(nrade*ipot,nbmax2)))
      allocate(dupper(neval2*max(nrade*ipot,nbmax2)))

      dlower = 0
      dmiddle = 0
      dupper = 0

      write(*,*) 'checkpoint 2 : in dmain'

200   format(///)
!     call to setfac
      call setfac(binom,nbin)

      if (ires.eq.0) then
!     zero tz and tx.....
         tz = x0
         tx = x0
         iblock=0

      else
!     ... or retrieve them for a restart run
         call rdscr(tz,tx,neval1*neval2,iscr,iblock)
      endif

      jblock=0
      kblock=iblock+nblock

!     read in radial dvr grid points (same for bra and ket)
      call getrow(r1,npnt1,ibra)
      if (idia .gt. -2) call getrow(r2,npnt2,ibra)
     write(*,*) 'checkpoint 3 : in dmain'
      if (jk2 .le. 1) then
        read(ibra)
        read(ibra)
        read(ibra)
        if (idia .eq. -2) then
           read(ibra)
           read(ibra)
        endif
!
      endif
           write(*,*) 'checkpoint 4 : in dmain'
        read(iket)
        if (idia .gt. -2) read(iket)
!
      if (jk1 .le. 1) then
        read(iket)
        read(iket)
        read(iket)
        if (idia .eq. -2) then
           read(iket)
           read(iket)
        endif
!
      endif
            write(*,*) 'checkpoint 5 : in dmain'
            write(*,*) 'Iam in dmain. size(nbass1)=', size(nbass1)
            write(*,*) 'Iam in dmain. size(nbass2)=', size(nbass2)
!     read in of energies for the ket
      read(iket) neval
      read(iket)(dum, i=1,ibase1 - 1),e1
      write(6,200)
      write(6,201)
201   format(5x,'First 10 energies for the ket in a.u.',/)
      if (.not.zpmin) write(6,203) (e1(i),i=1,min(10,neval1))

!     read in of energies for the bra

      read(ibra) neval
      read(ibra)(dum, i=1,ibase2 - 1),e2
      write(6,200)
      write(6,202)
202   format(5x,'First 10 energies for the bra in a.u.',/)
      if (.not.zpmin) write(6,203) (e2(i),i=1,min(10,neval2))
203   format(5d24.12,/)

      kbeg1= 0
      kbeg2= 0
      if (znco1.and.znco2) goto 54
      write(*,*) 'checkpoint 6 : in dmain'
      write(*,*) 'Iam in dmain. size(nbass1)=', size(nbass1)
      write(*,*) 'Iam in dmain. size(nbass2)=', size(nbass2)

!cccccccccccccccccccccccccccccc
!     nu = 0 calculation.
!ccccccccccccccccccccccccccccc
      nu = 0

!     call to lagpt

      write(6,6020)
6020  format(/,5x,'lagpt called after')
      call timer
      call lagpt(dipol,r1,r2,xd,wtd,nu)
!     call to trans
      write(6,603)
603   format(/,5x,'trans called after')
      call timer

!     get bra and ket properly lined up

      ks1= 1
      ip=ipar1
!     e to f calculation
      if (kmin1.gt.kmin2) then
         ks1= 2
         ip=1-ip
      endif
!cccccccccccccccccc!

!-----Djedjiga phase correction for Radau coordinate system
      if(idia.eq.-2) then
      ipar11=ipar1
      ipar22=ipar2
      if(kmin1.eq.0) ipar11=ipar1+1
      if(kmin2.eq.0) ipar22=ipar2+1
      endif
!------------------

      do 10 k1= ks1,jk1
      k2= k1 - kmin1 + kmin2
      if (k2.gt.jk2) goto 10
      jblock=jblock+1
      if (jblock.gt.iblock) then
       iblock=iblock+1
       write(*,*) 'checkpoint 7a : in dmain'
       write(*,*) 'size(nbass1)=', size(nbass1)
       write(*,*) 'size(nbass2)=', size(nbass2)
       write(*,*) 'k1, k2', k1, k2
       write(*,*) 'nbass1(k1) =', nbass1(k1)
       write(*,*) 'nbass1(k2) =', nbass2(k2)
       if (nbass1(k1) .eq. 0 .or. nbass2(k2) .eq. 0) then
         write(6,2020) iblock,k1-kmin1,k2-kmin2
 2020   format(/5x,'Block',i4,' k1 =',i3,' to k2 =',i3,' skipped')
      else
           write(*,*) 'checkpoint 7b : in dmain'
         kk=k1-kmin1
         call dsrd(dc1,dstemp,iket,mbass1,nbass1(k1),neval1,&
             k1,kbeg1,jk1,ip,ibase1,xd,kk,nu,ipar1)
         call dsrd(dc2,dstemp,ibra,mbass2,nbass2(k2),neval2,&
           k2,kbeg2,jk2,1-ip,ibase2,xd,kk,nu,ipar2)
         xfac= x1

         if (idia.eq.-2 .and. mod((kk+ipar11)/2+(kk+ipar22)/2,2).ne.0)&
           xfac=-xfac
         call trans(tz,dipol,binom,dc1,dc2,k1,k2,xfac,nu,1)
         call wrscr(tz,tx,neval1*neval2,iscr,iblock)
         write(6,2000) iblock,k1-kmin1,k2-kmin2
 2000    format(/5x,'Block',i4,' k1 =',i3,' to k2 =',i3,' completed')
         if (iblock .ge. kblock) goto 154
       endif
      endif
      ip=1-ip
10    continue

!cccccccccccccccccccccccccccccccccccccccccccc
!     nu = +/-1 calculation.
!ccccccccccccccccccccccccccccccccccccccccccccc
      nu = 1

!     call to lagpt
      write(6,602)
602   format(/,5x,'lagpt called after')
      call lagpt(dipol,r1,r2,xd,wtd,nu)

!     call to trans
      write(6,603)
      call timer

      j1= jk1 - kmin1
      j2= jk2 - kmin2

!     parities for symmetrised radau bisector embedding

      ip=ipar1
      do 11 k1= 1,jk1
      kk1= k1 - kmin1
      if (nbass1(k1).eq.0) goto 110

      if (jblock-iblock.gt.-2) then
	call dsrd(dc1,dstemp,iket,mbass1,&
          nbass1(k1),neval1,k1,kbeg1,jk1,ip,ibase1,xd,kk1,nu,ipar1)
      endif

!cccccccccccccccccccccccccccccccc
!
! 05/09/19: Modification by Dan Underwood
!
! Minimises calls to dsrd subroutine, and stores k-blocks from the bra wavefunction in the memory for faster access.
! For k in ket:
! dupper = k + 1 in bra
! dmiddle = k in bra
! dlower = k - 1 in bra
! The code cycles through these and avoids doing rewinds of the bra file inside the dsrd subroutine, which proves to be prohibitive for large files.
!
!cccccccccccccccccccccccccccccccc
      if (k1.eq.1) then
!	write(*,*) "Start modification"
	
	if (kmin1.eq.0) then
		
		if(kmin2.eq.kmin1) then
			call dsrd(dmiddle,dstemp,ibra,mbass2,nbass2(1),neval2,&
                                1,kbeg2,jk2,1-ip,ibase2,xd,1,1,ipar2)
            if(jk2 .ne. 1)  call dsrd(dupper,dstemp,ibra,mbass2,nbass2(2),neval2,&
                                2,kbeg2,jk2,ip,ibase2,xd,2,1,ipar2)
		else
			call dsrd(dlower,dstemp,ibra,mbass2,nbass2(1),neval2,& 
				1,kbeg2,jk2,ip,ibase2,xd,0,1,ipar2)
			call dsrd(dmiddle,dstemp,ibra,mbass2,nbass2(2),neval2,& 
				2,kbeg2,jk2,1-ip,ibase2,xd,1,1,ipar2)
if(jk2 .ne. 1) call dsrd(dupper,dstemp,ibra,mbass2,nbass2(3),neval2,& 
                   		3,kbeg2,jk2,ip,ibase2,xd,2,1,ipar2)
		endif
	else
		if(kmin2.eq.kmin1) then
			call dsrd(dmiddle,dstemp,ibra,mbass2,nbass2(1),neval2,&
                           1,kbeg2,jk2,1-ip,ibase2,xd,0,1,ipar2)
	if(jk2 .ne. 1) call dsrd(dupper,dstemp,ibra,mbass2,nbass2(2),neval2,&
                           2,kbeg2,jk2,ip,ibase2,xd,1,1,ipar2)
		else
			call dsrd(dupper,dstemp,ibra,mbass2,nbass2(1),neval2,& 
        	           1,kbeg2,jk2,ip,ibase2,xd,1,1,ipar2)
		endif
	endif
      endif

if(jk2 .le. 1) go to 108
!cccccccccccccccccccccccccccccccccccc
!     nu = +1 calculation
!cccccccccccccccccccccccccccccccccccc
      nu= 1
      kk2= kk1 + nu
      k2= kk2 + kmin2
      if (k2.le.jk2) then
       jblock=jblock+1
       if (jblock.gt.iblock) then
        iblock=iblock+1
        if (nbass2(k2).eq.0) then
         write(6,2020) iblock,k1-kmin1,k2-kmin2
        else
         xfac= -x1/sqrt(x2)
         if (kk1.eq.0) xfac= -x1
         if (idia.eq.-2 .and. mod((kk1+ipar11)/2+(kk2+ipar22)/2,2).ne.0)&
           xfac=-xfac
         if (.not. zembed .and. idia .lt. 0) xfac=-xfac
         !call dsrd(dc2,dstemp,ibra,mbass2,nbass2(k2),neval2,&
         !          k2,kbeg2,jk2,ip,ibase2,xd,kk2,nu,ipar2)
         call trans(tx,dipol,binom,dc1,dupper,k1,k2,xfac,nu,ip)
         call wrscr(tz,tx,neval1*neval2,iscr,iblock)
         write(6,2000) iblock,k1-kmin1,k2-kmin2
         if (iblock .ge. kblock) goto 154
        endif
       endif
      endif

108 continue
!cccccccccccccccccccccccccccccccccccc
!     nu = -1 calculation
!cccccccccccccccccccccccccccccccccccc
      nu= -1
      kk2= kk1 + nu
      k2= kk2 + kmin2
      if (k2.ge.1) then
       jblock=jblock+1
       if (jblock.gt.iblock) then
        iblock=iblock+1
        if (nbass2(k2).eq.0) then
         write(6,2020) iblock,k1-kmin1,k2-kmin2
        else
         xfac= x1/sqrt(x2)
         if (kk2.eq.0) xfac= x1
         if (idia.eq.-2 .and. mod((kk1+ipar11)/2+(kk2+ipar22)/2,2).ne.0)&
           xfac=-xfac
         if (.not. zembed .and. idia .lt. 0) xfac=-xfac
         !call dsrd(dc2,dstemp,ibra,mbass2,nbass2(k2),neval2,&
         !          k2,kbeg2,jk2,ip,ibase2,xd,kk2,nu,ipar2)
         call trans(tx,dipol,binom,dc1,dlower,k1,k2,xfac,nu,ip)
         call wrscr(tz,tx,neval1*neval2,iscr,iblock)
         write(6,2000) iblock,k1-kmin1,k2-kmin2
         if (iblock.ge.kblock .and. k1.lt.jk1) goto 154
        endif
       endif
      endif

110   ip=1-ip
      dlower = dmiddle
      dmiddle = dupper
      if (kbeg2.ne.jk2) then
	if (kmin1.eq.kmin2) then
		kk2 = kbeg2+(1 - INT((kmin1 + kmin2)/2))
	else
		kk2 = kbeg2+kmin1
	endif
        call dsrd(dupper,dstemp,ibra,mbass2,nbass2(kbeg2 + 1),neval2,&
                   kbeg2 + 1,kbeg2,jk2,ip,ibase2,xd,kk2,nu,ipar2)
      endif
11    continue

      goto 55
154   if (iblock.ge.mblock) goto 55
      write(6,1540) iblock
1540  format(//i7,' blocks calculated. dipole3 shutting down')
      call timer
      goto 55
54    continue

!     non coriolis coupled case
     write(*,*) 'checkpoint 8 : in dmain'
      nu= abs(kmin2-kmin1)
      nrad=nrado
      if (nu .eq. 1) nrad=nrade

!     call to lagpt
      write(6,602)
      call lagpt(dipol,r1,r2,xd,wtd,nu)
     write(*,*) 'checkpoint 9 : in dmain'
!     call to trans
      write(6,603)

      nu= (kmin2-kmin1)
      k1= kmin1
      k2= kmin2

      call dsrd(dc1,dstemp,ibra,mbass1,nbass1(1),neval1,&
                1,kbeg1,jk1,ipar1,ibase1,xd,k1,nu,ipar1)
      call dsrd(dc2,dstemp,ibra,mbass2,nbass2(1),neval2,&
                1,kbeg2,jk2,ipar2,ibase2,xd,k2,nu,ipar2)

      if (nu.eq.0) then
         xfac= x1
      else if (nu.eq.1) then
         xfac= -x1/sqrt(x2)
         if (k1.eq.0.or.k2.eq.0) xfac= -x1
      else if (nu.eq.-1) then
         xfac= x1/sqrt(x2)
         if (k1.eq.0.or.k2.eq.0) xfac= x1
      endif
      if (idia.eq.-2 .and. mod((k1+ipar1)/2+(k2+ipar2)/2,2).ne.0)&
           xfac=-xfac
      call trans(tx,dipol,binom,dc1,dc2,k1,k2,xfac,nu,ipar1)

!     end of transition dipole moment calculation

!     call to spect
55    write(6,604)
604   format(/,5x,'spect called after')
      call timer
      call spect(tz,tx,e1,e2,sint,xe2)

!     final time
       write(6,605)
605   format(/,5x,'program ended after')
      call timer
      deallocate(dstemp, dc1, dc2)
      stop
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **009
      subroutine setfac(binom,nbin)
!
!     setfa! initialises binomial array:
!        binom(i+1,j+1) = i! / (j! * (i-j)!)

      implicit double precision (a-h,o-y), logical (z)
      double precision, dimension(nbin,nbin) :: binom
      data x1/1.0d0/
      binom(1,1) = x1
      binom(2,1) = x1
      binom(2,2) = x1
      do 10 i=3,nbin
      binom(i,1) = x1
      binom(i,i) = x1
      i1 = i - 1
      do 20 j=2,i1
      binom(i,j) = binom(i1,j-1) + binom(i1,j)
   20 continue
   10 continue
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **011
      subroutine lagpt(d0,r1,r2,xd,wtd,nu)

!     subroutine lagpt obtains values of the dipole at the radial
!     dvr points and angular integration points

      implicit double precision(a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2

      double precision, dimension(*) :: d0
      double precision, dimension(npnt1) :: r1
      double precision, dimension(npnt2) :: r2
      double precision, dimension(npot) :: xd,wtd
      double precision, allocatable, dimension(:) :: b,c

      data x0/0.0d0/,toler/1.0d-8/,&
           x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/

       allocate( b(npot) )
       allocate( c(npot) )

       write(*,*) 'Checkpoint 1 in lagpt'

!     set up points & weights for npot point angular integration

      xnu= x0
      alf= xnu
      bta= alf
      do 10 i=2,npot
      xi= dble(i)
      b(i)= (alf+bta)*(bta-alf)/&
            ((alf+bta+x2*xi)*(alf+bta+x2*xi-x2))
      c(i)= x4*(xi-x1)*(alf+xi-x1)*(bta+xi-x1)*(alf+bta+xi-x1)/&
            ((alf+bta+x2*xi-x1)*(alf+bta+x2*xi-x2)*&
             (alf+bta+x2*xi-x2)*(alf+bta+x2*xi-x3))
10    continue

      call jacobi(npot,xd,wtd,alf,bta,b,c,csa,tsa)
      write(6,1000) npot,0,(xd(i),wtd(i),i=1,nn2)
 1000 format(//,i8,' point gauss-associated legendre integration',&
                  ' with k =',i2,&
             //,5x,'integration points',11x,'weights',&
             /,(f23.15,d25.12))
      write(6,1010) csa,tsa
 1010 format(/,5x,'computed sum of weights',d26.15,&
             /,5x,'exact    sum of weights',d26.15//)
      if (abs((csa-tsa)/tsa) .gt. toler) then
         write(6,940)
  940   format(/,5x,'gauss-legendre weights in error: adjust algorithm')
         stop
      endif
!     define other integration points
      do 20 i=1,nn2
      if (idia .ge.0) then
         j=i+nn2
         xd(j)=-xd(npot-j+1)
         wtd(j)=wtd(npot-j+1)
      else
         j=npot+1-i
         xd(j)=-xd(i)
         wtd(j)=wtd(i)
      endif
   20 continue
!     calculate dipole at (r1,r2,cos\theta)
      if (idia .eq. -2) then
         iadd=1-abs(nu)
      else
         i0=npnt1
      endif
      jdia = max(1, idia)
      ii=0
      do 30 i2=1,npnt2
      if (idia .eq. -2) then
         rr2=r1(i2)
         i0=i2-iadd
      else
         rr2=r2(i2)
      endif
      do 40 i1=1,i0
      do 50 j=1,ipot
      ii=ii+1
      call dipd(d0(ii),r1(i1),rr2,xd(j),nu)
      if (jdia.eq.2) d0(ii)=x2*d0(ii)
      d0(ii)=wtd(j)*d0(ii)
   50 continue
   40 continue
   30 continue
          write(*,*) 'Checkpoint 2 in lagpt'
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine asleg(pleg,pnorm,lmax,x,ipot,m)

!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, numerical recipes, p. 182,
!     for the polynomial part of associated legendre functions.
!     a factor of sin(theta)**m has NOT been removed from all functions.

      implicit double precision (a-h,o-y), logical (z)

      common /mass/ xmass(3),g1,g2,zembed,zbisc
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2

      double precision, dimension(ipot,0:lmax) :: pleg
      double precision, dimension(ipot) :: x
      double precision, dimension(0:lmax) :: pnorm

      data x1/1.0d0/,x2/2.0d0/

      do 10 i=1,ipot

      if (m .lt. 0 .or. abs(x(i)) .gt. x1) then
          write(6,200)
200      format(//5x,'improper argument in subroutine asleg'/)
          stop
      endif

      pmm = x1
      fact = x1
      somx2=sqrt((x1-x(i))*(x1+x(i)))

      do 11 j=1,m
      pmm = -pmm * fact * somx2
      fact = fact + x2
   11 continue
      pleg(i,0)= pmm
      pmmp1= x(i)*dble(m+m+1)*pmm
      pleg(i,1)= pmmp1
      ll=1

!loop ensures that same number of functions is calculated
!for both symmetry and no symmetry
      do 2 l= 2+m,(lmax+m)
      r2lm1 = dble(l+l-1)
      rlpmm1= dble(l+m-1)
      rlmm  = dble(l-m)
      pll= (x(i)*r2lm1*pmmp1 - rlpmm1*pmm)/rlmm

      pmm= pmmp1
      pmmp1= pll
      ll=ll+1
      pleg(i,ll)= pll
2     continue
10    continue

!     set up the normalisation constants
!     (pnorm)**2 = (2j + 1)/2   *   (j - k)! / (j + k)!
      jj = -1
      do 13 j = m,(lmax+m)
      fact = x1
      do 12 i = j-m+1,j+m
      facti = dble(i)
      fact = fact * facti
   12 continue
      jj = jj + 1
      pnorm(jj) = sqrt(dble(j+j+1) / (fact + fact))
   13 continue
!     now normalise the polynomials
      do 14 jj=0,lmax

       do 15 i=1,ipot
       pleg(i,jj) = pleg(i,jj) * pnorm(jj)
15     continue
   14 continue
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                               **018
      subroutine jacobi(nn,x,a,alf,bta,b,c,csa,tsa)

!     calculates zeros x(i) of the nn'th order jacobi polynomial
!     pn(alf,bta) for the segment (-1,1) & and corresponding weights
!     for gauss-jacobi integration. this routine, and those which
!     follow, are due to a.h.stroud and d. secrest, "gaussian
!     integration formulas", 1966, prentice hall, page 29.
!     note that for our purposes, alf= bta= nu.

      implicit double precision(a-h,o-z)
      double precision, dimension(nn) :: x,a,b,c
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/,x6/6.0d0/,&
          x8/8.0d0/,eps/1.0d-12/
      fn= dble(nn)
      nn2= nn/2
      csa= x0
      beta= x1
      if (alf.eq.1) beta = beta/x6
      cc= beta*x2**(alf+bta+x1)
      tsa= cc/x2
      do 10 i=2,nn
      cc= cc*c(i)
   10 continue
      do 20 i=1,nn2
      if (i .eq. 1) then
!         largest zero
      an= alf/fn
      bn= bta/fn
      r1= (x1 + alf)*(2.78d0/(x4 + fn*fn) +0.768*an/fn)
      r2= x1 + 1.48d0*an + 0.96d0*bn + 0.452*an*an + 0.83d0*an*bn
      xt= x1 - r1/r2
      else if (i .eq. 2) then
!         second zero
      r1= (4.1d0 + alf)/((x1 + alf)*(x1 + 0.156*alf))
      r2= x1 + 0.06d0*(fn - x8)*(x1 + 0.12d0*alf)/fn
      r3= x1 + 0.012*bta*(x1 + abs(alf)/x4)/fn
      ratio= r1*r2*r3
      xt= xt - ratio*(x1 - xt)
      else if (i .eq. 3) then
!         third zero
      r1= (1.67d0 + 0.28d0*alf)/(x1 + 0.37d0*alf)
      r2= x1 + 0.22d0*(fn - x8)/fn
      r3= x1 + x8*bta/((6.28d0 + bta)*fn*fn)
      ratio= r1*r2*r3
      xt= xt - ratio*(x(1) - xt)
      else
!         middle zeros
      xt= x3*x(i-1) - x3*x(i-2) + x(i-3)
      endif

      call root(xt,nn,alf,bta,dpn,pn1,b,c,eps)
      x(i)= xt
      a(i)= cc/(dpn*pn1)
      csa= csa + a(i)
20    continue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **019
      subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x.

      implicit double precision(a-h,o-z)
      double precision, dimension(nn) :: b,c
      iter= 0
1     iter= iter + 1
      call recur(p,dp,pn1,x,nn,alf,bta,b,c)
      d = p/dp
      x = x - d
      if (abs(d) .gt. eps .and. iter .lt. 10) goto 1
      dpn= dp
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **020
      subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
      implicit double precision(a-h,o-z)
      double precision, dimension(nn) :: b,c
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
      p1= x1
      p= x + (alf-bta)/(alf + bta + x2)
      dp1= x0
      dp= x1
      do 10 j=2,nn
      q= (x - b(j))*p - c(j)*p1
      dq= (x - b(j))*dp + p - c(j)*dp1
      p1= p
      p= q
      dp1= dp
      dp= dq
  10 continue
      pn= p
      dpn= dp
      pn1= p1
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **021
      subroutine dsrd(d,temp,ivec,mmbass,nbass,ne,kneed,kbeg,&
                     jk,ipar,ibase,xd,kk,nu,jay_ipar)

!     subroutine to read d coefficients from dstore data

      implicit double precision (a-h,o-y), logical(z)
      parameter (iz=1)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common /mass/ xmass(3),g1,g2,zembed,zbisc
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2

      double precision, dimension(ne,max(nrade*ipot,nbass)) :: d
      double precision, dimension(ne,nbass) :: temp
      double precision, dimension(npot) :: xd
      double precision, allocatable, dimension(:,:) :: plegd
      integer, allocatable, dimension(:) :: iv

      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/

      nang=nbass/nrade
      if (ipar.eq.1 .and. zbisc) nang=nbass/nrado
      nrad=nrado
      jdia=max(1,idia)
      if (abs(nu).eq.1 .and. ipar.eq.0) nrad=nrade

      if (jk .gt. 1) then
         kz=kk
!
         if (kneed .le. kbeg) then
           rewind ivec
	   write(*,*) "Rewinded ivec"

           do 10 i=1,kneed+6
           read(ivec)
10        continue
           if (.not. zbisc) read(ivec)
         else
           do 20 i=kbeg,kneed-2
           read(ivec)
20        continue
         endif

!        for nu=0 drop symmetric grid points
         if (nu.eq.0 .and. zbisc .and. ipar.eq.0) then

            read(ivec)(dum,i=1,nbass*(ibase-1)),&
                  ((d(i,j),j=1,nbass),i=1,ne)
            ipt=0
            jpt=0
            do 12 n1=1,npnt1
            do 13 n2=1,n1-1
            do 14 j=1,nang
            ipt=ipt+1
            jpt=jpt+1

            do 15 ie=1,ne
            temp(ie,ipt)=d(ie,jpt)
   15       continue
   14       continue
   13       continue
            jpt=jpt+nang
   12       continue
          else
            read(ivec)(dum,i=1,nbass*(ibase-1)),&
                  ((temp(i,j),j=1,nbass),i=1,ne)
          endif
      else

         rewind ivec

         do 30 i=1,6
         read(ivec)
30       continue
         if (.not. zbisc) read(ivec)
         read(ivec) kz,maxleg,nidvr
         leng=(maxleg+1)*nidvr
         allocate(plegd(0:maxleg,nidvr),iv(nidvr))
         call getrow(plegd,leng,ivec)
         if (zbisc) then
           read(ivec) iang,ibass
           read(ivec) (iv(ii),ii=1,nidvr)
         else
            iang=nidvr
            ibass=mmbass
            do 40 i=1,nidvr
            iv(i)=1
40          continue
         endif
         read(ivec)
         read(ivec)
         call jtran(temp,nrad,ne,plegd,maxleg,nidvr,kz,d,&
                   ivec,ipar,iv,iang,ibass,ibase,nu,d(ne/2,1),jay_ipar)
         deallocate(plegd,iv)
      endif

!     calculate the legendre function
!
      allocate(plegd(ipot,0:(nang*jdia)-1))
      call asleg(plegd,d,(nang*jdia)-1,xd,ipot,kz)

      if (idia.eq.2) then
        if(zembed) then
         !r2 case
           jfirst = mod(kk+ipar1,2) !ipar1 used as in r2 case parity of basis
                                    !same for bra and ket
        else
           jfirst = mod(kk+jay_ipar,2)
        endif
       DO i=1, ipot
          !select polynomials to use
         index = jfirst

         DO ll=0,nang-1
            plegd(i,ll)=plegd(i,index)
            index = index + 2
         END DO
      END DO
      endif


!     evaluate wavefunction at angular integration points

      beta=x1
      d=x0
      ipt=1
      jpt=1
!
      do 60 mn=1,nrad
      call dgemm('n','t',ne,ipot,nang,beta,temp(1,jpt),ne,plegd,&
                ipot,beta,d(1,ipt),ne)
      ipt=ipt+ipot
      jpt=jpt+nang
   60 continue
      kbeg=kneed
      deallocate(plegd)
!     check normalisation if requested
      if (zprint) then
         xnorm=x0
         i=0
         do 70 i=1,nbass
         xnorm=xnorm+d(1,i)**2
   70    continue
          if (idia.eq.2) xnorm = x2 * xnorm
          write(6,*) ' Vector 1 with k =',kk,ipar,ivec,&
                   ' contribution to  normalisation is',xnorm
      endif
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jtran(coef,nrad,mvib,pleg,maxleg,idvr,kz,dvrvec,&
                        ivec,ipar,iv,iang,ibass,ibase,nu,temp,jay_ipar)

      implicit double precision (a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      common /mass/ xmass(3),g1,g2,zembed,zbisc

      double precision, dimension(0:maxleg,idvr) :: pleg
      double precision, dimension(iang,*) :: dvrvec
      double precision, dimension(nrad) :: sumk
      double precision, dimension(mvib,*) :: coef
      dimension iv(idvr)
      double precision, dimension(iang,*) :: temp

      data x0/0.0d0/

!     transform back to the original fbr-type basis in the
!     associated legendre functions
      jstart=kz
      jdia=max(idia,1)
      jj0=-jdia
      if (zembed) then
         if (idia .eq. 2 .and. mod(jstart,2) .ne. ipar1) then !r2 case
           jj0=jj0+1
           jstart=jstart+1
         endif
      else
         if (idia .eq. 2 .and. mod(jstart,2) .ne. jay_ipar) then !r1 case
            jj0=jj0+1
            jstart=jstart+1
         endif
       endif
      nang=(maxleg+kz-jstart)/jdia+1
      do 5 l=1,ibase-1
      read(ivec)
5     continue
      do 10 l=1,mvib
!     first read in a new vector
      if (nu.eq.0 .and. idia.eq.-2 .and. ipar.eq.0) then
         call getrow(temp,ibass,ivec)
         ipt=0
         jpt=0
         do 13 n1=1,npnt1
         do 12 n2=1,n1-1
         jpt=jpt+1
         ipt=ipt+1
         do 11 j=1,iang
         dvrvec(j,ipt)=temp(j,jpt)
   11    continue
   12    continue
         jpt=jpt+1
   13    continue
      else
         call getrow(dvrvec,ibass,ivec)
      endif
      jj=jj0
      do 20 j=1,nang
      sumk=x0
      jj=jj+jdia
      kk=0
      do 40 k=1,idvr
      if (iv(k) .le. 0) goto 40
      kk=kk+1
      do 50 mn=1,nrad
      sumk(mn)=sumk(mn) + dvrvec(kk,mn) * pleg(jj,k)
   50 continue
   40 continue

      ipt=j
      do 60 mn=1,nrad
      coef(l,ipt) = sumk(mn)
      ipt=ipt+nang
   60 continue
   20 continue

   10 continue

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **022
      subroutine trans(t,dipol,binom,dc1,dc2,k1,k2,xfac,nu,ipar)

!     subroutine trans is the main working routine in program dipole.
!     it carries out the necessary angular integrations using 3-j
!     symbols as defined by brink and satchler "angular momentum".
!     trans is called by dmain for each k-k' overlap integral.

!     this part of the program calculates the transition moments
!     leaving out a factor of root$(2j'+1)(2j"+1) which is a
!     common factor which is put in the subroutine spect

!     TRANS uses BLAS rank-1-update routine DGER.
!     Adapted to run in parallel on SGI Origin machines by Greg Harris
!     In this case NCPUS should be set to the number of processors.

      implicit double precision (a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      parameter (NCPUS=1)

      double precision, dimension(neval1,neval2) :: t
      double precision, dimension(*) :: dipol
      double precision, dimension(nbin,nbin) :: binom
      double precision, dimension(neval1,*) :: dc1
      double precision, dimension(neval2,*) :: dc2
      double precision, allocatable, dimension(:,:,:) :: ttemp

      data x0/0.0d0/

      if (ncpus .gt. 1) then
         allocate(ttemp(neval1,neval2,ncpus))
         ttemp=x0
      endif

      if (znco1 .and. znco2) then
         j1= jk1
         j2= jk2
         kk1= kmin1
         kk2= kmin2
      else
         j1= jk1 - kmin1
         j2= jk2 - kmin2
         kk1= k1 - kmin1
         kk2= k2 - kmin2
      endif

!     start the calculation

      x1= threej(j1,1,j2,kk1,nu,-kk2,binom,nbin)*xfac
      if (mod(kk1,2) .ne. 0) x1=-x1
      if (zprint) write(6,*) 'j1, k1, nu, j2, k2 ',j1,kk1,nu,j2,kk2
      if (zprint) write(6,*) 'xfac, x1 =',xfac, x1

      i0=1
      if (abs(nu).eq.1 .and. ipar.eq.0) i0=0
      i=0
      id=0
      n0=npnt2

      do 10 nn1=1,npnt1
      if (idia .eq. -2) n0=nn1-i0
      niter=n0*ipot

! serial case: whole calculation for NCPUS=1 or else the remainder
      if (ncpus .eq. 1) then
         irem=niter
         nparr=0
      else
         irem=mod(niter, ncpus)
         nparr=(niter-irem)/ncpus
      endif

      i2 = i
      id2 = id

      do 101 n1=1,irem
            ii=i2+n1
            idd=id2+n1
            x3=x1*dipol(idd)
      call dger(neval1,neval2,x3,dc1(1,ii),1,dc2(1,ii),1,t,neval1)
 101  continue

      If (nparr .gt. 0) then

! Paralell part
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,j,idd,x3,n1,i2,id2)
         do 20 j=1,ncpus
         i2 = i + irem + ((j-1)*nparr)
         id2 = id + irem + ((j-1)*nparr)
         do 30 n1=1,nparr
            ii=i2+n1
            idd=id2+n1
            x3=x1*dipol(idd)
            call dger(neval1,neval2,x3,dc1(1,ii),1, &
            dc2(1,ii),1,ttemp(1,1,j),neval1)
30    continue
20    continue
!$OMP END DO
!$OMP END PARALLEL
      endif

      i=i+niter
      id=id+niter
!     take care of diagonal case in symmetric Radau
      if (abs(nu).eq.i0 .and. idia .eq. -2) id=id+npot
10    continue

      if (ncpus .gt. 1) then
!  matrix summation
!$OMP PARALLEL
!$OMP DO PRIVATE(ij, jk, ki)
         do 40 jk=1,neval2
         do 41 ki=1,ncpus
         do 42 ij=1,neval1
         t(ij,jk)=t(ij,jk)+ttemp(ij,jk,ki)
42       continue
41       continue
40       continue
!$OMP END DO
!$OMP END PARALLEL
         deallocate(ttemp)
      endif

      return
      end
!                                                **024
      function threej(j1,j2,j3,m1,m2,m3,binom,nbin)

      implicit double precision(a-h,o-z)

      double precision, dimension(nbin,nbin) :: binom
      data zero,one/0.0d0,1.0d0/

      threej = zero
      if (m1+m2+m3 .ne. 0) return
      i1 = -j1+j2+j3+1
      if (i1 .le. 0) return
      i2 = j1-j2+j3+1
      if (i2 .le. 0) return
      i3 =  j1+j2-j3+1
      if (i3 .le. 0) return
      k1 =  j1+m1+1
      if (k1 .le. 0) return
      k2 = j2+m2+1
      if (k2 .le. 0) return
      k3 =  j3+m3+1
      if (k3 .le. 0) return
      l1 = j1-m1+1
      if (l1 .le. 0) return
      l2 = j2-m2+1
      if (l2 .le. 0) return
      l3 = j3-m3+1
      if (l3 .le. 0) return
      n1 = -j1-m2+j3
      n2 = m1-j2+j3
      n3 = j1-j2+m3
      imin = max(-n1,-n2,0)+1
      imax = min(l1,k2,i3)
      if (imin .gt. imax) return
      sign = one

      do 20 i=imin,imax
      sign = -sign
      threej = threej + sign*binom(i1,n1+i)*binom(i2,n2+i)*binom(i3,i)
   20 continue
      threej = threej * sqrt(binom(j2+j2+1,i3)*binom(j1+j1+1,i2)&
             / (binom(j1+j2+j3+2,i3)*dble(j3+j3+1)&
             * binom(j1+j1+1,l1)*binom(j2+j2+1,l2)*binom(j3+j3+1,l3)))
      if (mod(n3+imin,2) .ne. 0) threej = - threej
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **025
      subroutine spect(tz,tx,e1,e2,sint,xe2)

!     subroutine spect calculates s(f-i), the line strength,
!     and controls the printing out of the transition moments,
!     frequencies and nuclear spin effects.
!
!     if ztra is .true. it calls outpt to output energies
!     and line strengths for program spectrum to calculate
!     simulated spectra

      implicit double precision (a-h,o-y), logical (z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      common /head/ title
      common /mass/ xmass(3),g1,g2,zembed,zbisc

      double precision, dimension(neval1,neval2) :: tz,tx
      double precision, dimension(neval1) :: e1
      double precision, dimension(neval2) :: e2
      double precision, dimension(neval1,neval2) :: sint
      double precision, dimension(neval2) :: xe2

      character(len=8)  title(9)
      data autocm/2.19474624d+05/,x0/0.0d0/,&
           autode/2.5417662d0/,&
           detosec/3.136186d-07/
!     autocm converts atomic units to wavenumbers
!     autode converts atomic units to debye
!     detose! converts from s(f-i) in debye**2 to seconds

      if (znco1.and.znco2) then
         j1= jk1
         j2= jk2
      else
         j1= jk1-kmin1
         j2= jk2-kmin2
      endif
      write(6,200)
200   format(///)
      write(6,201)
201   format(//,5x,'*************************************************'&
             //,5x,'print out of dipole transition moments and s(f-1)'&
             //,9x,'frequencies in wavenumbers',&
              /,9x,'transition moments in debye (2.54174a.u.)',&
              /,9x,'s(f-i) in debye**2',&
              /,9x,'einstein a-coefficient in sec-1',//,&
           5x,'*************************************************')

      write(6,200)
      write(6,202) title
202   format(5x,9a8)
      write(6,200)
      write(6,203) j1, kmin1, j2, kmin2, idia, ipar1, ipar2
203   format(5x,'j1=',i4,'  kmin1=',i4,'  j2=',i4,'  kmin2=',i4,&
               '  idia=',i4,'  ipar1=',i4,'  ipar2=',i4,//)
      ezero=x0
      read(5,505,end=555) ezero
  505 format(f20.0)
      write(6,204) ezero
204   format(5x,'ground zero =',e16.8,' cm-1')
  555 continue

      write(6,200)
      write(6,205)
205   format(/,' ie1 ie2   ket energy   bra energy    frequency  z trans', &
      'ition    x transition       dipole       s(f-i)      a-coefficient' ,/)

!     xf is the factor that allows for the root$(2j' + 1)(2j"+1)^G
!     left over from the calculation of the transition moment

      xf= sqrt(dble((2*j1+1)*(2*j2+1)))

!     calculate transition moments, line strengths and a-coefficients

! dd   : transition frequency in cm-1
! dd3  : (transition frequency )**3
!
! sxd  : s(f-i) line strenght
! sxd  = (2*j1+1)*(2*j2+1)*(tz(ie1,ie2)+tx(ie1,ie2))**2
! sxd  : in Debye**2
!
! sint : (2*j1+1)*(2*j2+1)*(tz(ie1,ie2)+tx(ie1,ie2))**2
! sint : in au
! sint is wha will be printed to fort.13
! So fort.13 will contain :
!
!       write(itra) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed
!       gz is ezero  in cm-1

      do 1 ie1=1,neval1
      xe1= e1(ie1)*autocm - ezero

      do 2 ie2=1,neval2

      if (ie1.eq.1) xe2(ie2)= e2(ie2)*autocm - ezero
      dd= xe2(ie2) - xe1

      dd3= abs(dd*dd*dd)
      sx= (tz(ie1,ie2) + tx(ie1,ie2))**2
      sint(ie1,ie2)= sx*xf*xf
      tzd= tz(ie1,ie2)*autode*xf
      txd= tx(ie1,ie2)*autode*xf
      if (.not.zbisc .and. zembed) txd = -txd
      t= abs(tzd + txd)
      sxd= t*t

      if (dd .gt. x0) a= sxd*dd3*detosec/dble(2*j2 + 1)
      if (dd .lt. x0) a= sxd*dd3*detosec/dble(2*j1 + 1)

      if (zpmin .and. max(ie1,ie2).gt.10) goto 2
      write(6,206) ie1,ie2,xe1,xe2(ie2),dd,tzd,txd,t,sxd,a
!206   format(2(i4),3(3x,f10.3),5(2x,e13.6))
206   format(2(i4),3(3x,f14.6),5(2x,es15.8)) ! changed L Lodi 8-Feb-2010
2     continue

!      if (.not.zpmin .or. ie1.le.10) write(6,207)
!207   format(//)

1     continue

! writes in itra e1 and e1 : energy values in au
! then it will write the values of sint ( in au as weel)

      if (ztra) call outpt(e1,e2,sint,ezero)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **026
      subroutine outpt(e1,e2,sint,gz)

!     subroutine outpt outputs the data necessary for program
!     spectrum to simulate laboratory or interstellar spectra.
!     the output data is in atomic units.

      implicit double precision(a-h,o-y), logical(z)
      common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                 npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
                 kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
      common/sym/ idia,ipar1,ipar2,jrot1,jrot2
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      common /mass/ xmass(3),g1,g2,zembed,zbisc

      double precision, dimension(neval1) :: e1
      double precision, dimension(neval2) :: e2
      double precision, dimension(neval1,neval2) :: sint

      if (znco1.and.znco2) then
         j1= jk1
         j2= jk2
      else
         j1= jk1-kmin1
         j2= jk2-kmin2
      endif

!     is this the first write-out?

!      open(unit=itra,form='unformatted',recordtype='segmented')
      if (.not.zstart) then
10       read(itra,end=90)
         goto 10
90       continue
! *****  inclusion of the following card is machine dependent *****
!        backspace itra
      endif

      write(itra) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,&
                 zembed,ibase1,ibase2
      write(itra) e1
      write(itra) e2
      do 20 ie2=1,neval2
      call outrow(sint(1,ie2),neval1,itra)
20    continue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **027
      subroutine getrow(row,nrow,iunit)

      implicit double precision (a-h,o-y)
      double precision, dimension(nrow) :: row
      read(iunit) row
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **028
      subroutine outrow(row,nrow,iunit)

      implicit double precision (a-h,o-y)
      double precision, dimension(nrow) :: row
      write(iunit) row
      return
      end
!cccccccccccccccccccccccccccccccccccccc
      subroutine rdscr(t1,t2,ndim,iscr,iblock)
!     read restart data stored on unit iscr
      implicit double precision (a-h,o-y), logical (z)
      double precision, dimension(ndim) :: t1,t2
      read(iscr) iblock
      read(iscr) t1
      read(iscr) t2
      return
      end
!cccccccccccccccccccccccccccccccccccccccc
      subroutine wrscr(t1,t2,ndim,iscr,iblock)
!     write restart data to unit iscr
      implicit double precision (a-h,o-y), logical (z)
      double precision, dimension(ndim) :: t1,t2
      rewind iscr
      write(iscr) iblock
      write(iscr) t1
      write(iscr) t2
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine timer
!     prints current cpu time usage                                 #030
!     needs a subroutine which can access the machine clock
      common/timing/itime0

      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/)
      return
      end
      SUBROUTINE DIPD(DIPC,R1,R2,XCOS,NU)
!
!     TRANSFORM GENERALISED COORDINATES TO THOSE FOR PARTICULAR
!     SYSTEM. THIS VERSION TRANSFORMS TO AB2 BONDLENGTH-BONDANGLE
!     COORDINATES. ALLOWANCE MUST BE MADE FOR THE NUMBERING OF THE ATOMS
!     Additionally, the zbisc option is included.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), LOGICAL (Z)
      COMMON /MASS/ XMASS(3),G1,G2,zembed,zbisc
      LOGICAL FIRST/.TRUE./
      SAVE FIRST

      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/,PI/3.1415927D0/


      IF (G1 .EQ. X0) THEN       !        BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
         Q1 = R1
         Q2 = R2
         THETA = ACOS(XCOS)
      ELSE IF (G2 .EQ. X0) THEN  !        SCATTERING COORDINATES: ATOM 2 = ATOM 3
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         ALPHA= ACOS(XCOS)
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
         XSIN= SQRT(1.0D0 - XCOS*XCOS)
         THETA = ACOS(COST)
         BETA= ASIN(XSIN*YY/Q2)
      ELSE                     !   GENERAL COORDINATES (INCLUDING RADAU): ATOM 1 = ATOM 2
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

      CALL DIPSA(DIPX,DIPY,Q1,Q2,cost)

!     BONDLENGTH-BONDANGLE CO-ORDINATES
!
      IF (G1.EQ.X0) THEN
         GAMMA= THETA/X2
         ycos= COS(GAMMA)
         ysin= SIN(GAMMA)
         IF (ZEMBED) THEN
            IF (NU.EQ.0) THEN
               DIPC= +DIPY*ycos - DIPX*ysin
            ELSE
               DIPC= +DIPX*ycos + DIPY*ysin
            ENDIF
         ELSE
            if (NU.EQ.0) THEN
               DIPC= +DIPY*ycos + DIPX*ysin
            ELSE
               DIPC= -DIPX*ycos + DIPY*ysin
            ENDIF
         ENDIF

!     SCATTERING CO-ORDINATES

      ELSE if (G2.EQ.X0) THEN
         GAMMA= BETA - THETA/x2
         if (ZEMBED) THEN
            ycos= COS(GAMMA)
            ysin= SIN(GAMMA)
            if (NU.EQ.0) THEN
               DIPC= -DIPX*ysin - DIPY*ycos
            ELSE
               DIPC= +DIPX*ycos - DIPY*ysin
            ENDIF
         ELSE
            DELTA= ALPHA - GAMMA
            ycos= COS(DELTA)
            ysin= SIN(DELTA)
            if (NU.EQ.0) THEN
               DIPC= +DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= +DIPX*ycos - DIPY*ysin
            ENDIF
         ENDIF

!     ALL OTHER CO-ORDINATES

      ELSE
         IF (ZBISC) THEN
            H1 = G1*Q1
            COSA = (R2*R2 + Q2*Q2 - H1*H1)/(X2*R2*Q2)
            alpha = (acos(xcos)-theta)/x2 - acos(cosa)
            ycos= - COS(ALPHA)
            ysin= + SIN(ALPHA)
            if (NU.EQ.1) THEN
               DIPC= +DIPX*ysin - DIPY*ycos
            ELSE
               DIPC= -DIPX*ycos - DIPY*ysin
            ENDIF
         ELSEIF (ZEMBED) THEN
            H1 = G1*Q1
            COSA = (R2*R2 + H1*H1 - Q2*Q2)/(X2*R2*H1)
            ALPHA= ACOS(COSA)
            ycos= - COS(ALPHA + THETA/X2)
            ysin= + SIN(ALPHA + THETA/X2)
            if (NU.EQ.0) THEN
               DIPC= -DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= +DIPX*ycos + DIPY*ysin
            ENDIF
         ELSE
            H2 = G2*Q2
            COSA = (R1*R1 + H2*H2 - Q1*Q1)/(X2*R1*H2)
            alpha = ACOS(COSA)
            ycos= - COS(alpha + THETA/X2)
            ysin= + SIN(alpha + THETA/X2)
            if (NU.EQ.0) THEN
               DIPC= +DIPX*ysin + DIPY*ycos
            ELSE
               DIPC= -DIPX*ycos + DIPY*ysin
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END


