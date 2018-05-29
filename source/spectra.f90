      program spect4

!  this program is designed to handle the output from programs dipole
!  (of the TRIATOM suite) and dipole3 (of the DVR3D suite)
!  for several runs and turn it into data sorted on the transition
!  frequency wif. From this, and input from triatom/dvr3d and rotlev runs,
!  the partition functions and integrated absorption coefficients at
!  the required temperatures.

!  the program structure consists essentially of four main parts:

!  1) the main program, spect, determines the amount of memory needed
!  for the job in hand, and calls subroutine gtmain, a machine
!  dependent routine, which allocates the memory needed. spect also
!  calls the other three main working routines.

!  2) subroutine sortsp. this routine sorts out data from program
!  dipole on ascending transition frequency and may print this out
!  if required.

!  the main input file is "itra", and contains the output from the
!  necessary dipole runs. by default itra is attached to stream 13.
!  itra is read in the main program to determine the size of the
!  problem. the main program also creates a scratch file, item (16),
!  to cut down the amount of data sorted in sortsp.

!  subroutine sortsp creates a scratch file ispe for calculations of
!  spectral intensities, if required. if this file has been saved from
!  a previous run, sortsp may be skipped. ispe defaults to stream 15.

!  3) subroutine pfcalc. this routine calculates the partition function
!  from energy levels calculated by programs triatom and rotlevd at the
!  temperature required for each run.

!  data necessary to calculate the partition function is given on
!  file "ilev". this is defaulted to stream 14.
!  the first eigenenergy on ilev must be the ground state for j=0,
!  (q=0). this is the zero energy of the problem. ilev is read in the
!  main program to determine the size of the problem and in subroutine
!  pfcal! to access the data.

!  4) subroutine spectm. this routine calculates boltzman weighted
!  intensities at the required temperatures and in the frequency ranges
!  determined by the input parameters given by input data lines 4-n.

!  spectm reads the sorted data from sortsp off scratch file ispe.
!  if required it writes out data for a suitable graphics package to
!  data stream iplot (default 20). one option is to produce a spectrum
!  of gaussian profiles by calling subroutine profil. users may wish to
!  substitute this routine with other profile types.
!  lines can also be written to a linelist file, stream ilist (default
!  36), if zlist is .true.


!  input data is supplied to the program on stream 5, as follows:

!  1) namelist prt: four logical parameters:

!                   zout(false) this should be true if the sorted
!                               line strengths are to be written
!                               to the lineprinter. zout is set to
!                               true automatically if zspe is false.

!                   zsort(true) if false subroutine sortsp is skipped.

!                   zspe (true) if false the program stops after sortsp.
!                               units of ispe are atomic units.

!                   zpfun(true) calculates the partition function
!                               from energy levels supplied from
!                               DVR3DRJZ and ROTLEV3/3B.
!                               if zpfun false, the partition function
!                               is set to q read in below. 
!                   the default values of itra(13), ilev(14), ispe(15),
!                   and item(16) may also be reset on this namelist
!                   if required.

!                   gz(0.0)     absolute energy of the ground state
!                               in cm-1 (not used if zpfun=.true.).

!   If zsort=.true. the range of data to be sorted may be set by setting
!                   any or all of wsmin, wsmax, emin, emax, jmax, smin.

!                   emin(-1.d27) the minimum value of e" for which data
!                                is to be written out.

!                   emax(+1.d27) the maximum value of e" for which data
!                                is to be written out. if emax is not
!                                specified all values of e" above emin
!                                are considered.

!                   jmax(-1)    the maximum value of j" for which
!                               transition frequencies, et! are to be
!                               written out in subroutine sortspe.
!                               if jmax is not specified, all j" values
!                               are considered.

!                   smin(0.0d0) The minimum value of the square of the 
!                               transition dipole (D**2).
!                               
!                   smax        The maximum value of the square of the
!                               transition dipole (D**2)        
!
!                   wsmax       The maximum frequency of a transition (cm-1).
!
!                   wsmin       The minimum frequency of a transition (cm-1)

!  judicious use of the parameters on this namelist can enable certain types
!  of transitions to be separated - e.g. fundamental bands, hot bands etc.

!  2) a title, format 9a8.

!  3) ge, go. format 2f10.0   nuclear degeneracy factors for homonuclear
!             diatomic molecules for the j even and j odd states.
!             defaults of 1.0d0 supplied if left totally blank.
!             if the molecule does not contain a homonuclear diatomic
!             part, ge and go are ignored.

!  4) temp, xmin, wmin, wmax, dwl, q. format 6f10.0

!      temp       temperature in k;
!      xmin       the lowest relative intensity to be printed out;
!      wmin and wmax define the transition frequency range in cm-1.
!      dwl        the width of the gaussian profile for the lines (see below).
!      q          is the partition function (only used if ZPFUN=F).
!                 use of default (q=1.0) means the absolute absortpion 
!                  coefficients are not correct. 

!  5) namelist spe. this namelist is read in subroutine spectm.
!     emin1, emax1, emin2, emax2, tinte and jmax control the print out
!     of spectral intensities.

!     tinte set the value below which the intensities are printed out.

!     zemit (.false.) gives emissivities if .true.

!     in addition, there is a parameter zplot(false)  which creates an
!     ascii file suitable for use with a graphics package. this can
!     consist of stick spectra or a gaussian broadened spectrum if
!     zprof (.false.) is .true. at npoint points, it is in:

!     1) frequency in wavenumbers against intensity if zfreq (.true.)
!        is .true.

!     2) wavelength in microns against intensity if zfreq .false.

!     3) either of the above against einstein a-coefficient times
!        spin weighting in zeinst (.false.) is .true.
!        zeinst (.true.) sets zemit to .true.
!     in addition if zene is .true the energy levels and the rotational
!     levels are printed in the output file iplot.

!  if zprof is .true., stream idat takes data from spectm to profil.

!  The program makes use of the dpsort freeware routines, found at
!  http://netlib2.cs.utk.edu/
! 
!  Updated to f90 to use dynamic memory allocation and to preselect transitions
!  by GJH & JT 2001.

      implicit double precision(a-h,o-y), logical(z)
      character(len=8) title(9)
      namelist/prt/ zout, zsort, zspe, zpfun, itra, ilev, ispe, item, &
                    wsmax,wsmin, emin, emax, jmax, smin, gz
      common/logic/ zout, zsort, zspe, zpfun, zembed
      common/base/ ibase1,ibase2
      common/timing/itime0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ee1, ee2, ss
      data autocm/ 2.19474624d+05/, autode/ 2.5417662d0/

!     set time zero for timer
      call SYSTEM_CLOCK(itime0,irate2,imax2)      

!     preset namelist defaults

      zout= .false.
      zsort= .true.
      zspe= .true.
      zpfun= .true.
      itra= 13
      ilev= 14
      ispe= 15
      item= 16
      wsmin= 0.0d0
      wsmax= 1.0d6
      smin=0.0d0
      jmax=500
      emin= -1.0d27
      emax= +1.0d27
      gz=0.0d0

      write(6,*) "Program SPECTRA (version of Feb 2004)"

!     read in logical variables and title
      read(5,prt)
      read(5,100) title
100   format(9a8)
      write(6,202) title
202   format(/5x,9a8/)

!     determine energy zero from levels file if available
      if (zpfun) then
         open(unit=ilev,form='formatted',status='old')
   
         write(6,201) ilev
201     format(5x,'Partition function: levels read from ilev =',i3/)
         read(ilev,*) j
         read(ilev,*) gz
         gz=gz*autocm
      endif
      write(6,203) gz
203   format(5x,'Ground state energy set to',f14.5,' cm-1')
!     determine number of transitions to be processed and create item

      maxr=0
      nr= 0
      nmax1= 0
      nmax2= 0
      open(unit=ispe,form='unformatted')!,recordtype='segmented')
      if(zsort) then
         write(6,1000) wsmin,wsmax,emin,emax,smin,jmax
 1000    format(/5x,'Sorting requested, list preselected using:', &
                /5x,'Minimum frequency,   WSMIN =',d12.5,' cm-1', &
                /5x,'Maximum frequency,   WSMAX =',d12.5,' cm-1', &
                /5x,'Minimum energy,       EMIN =',d12.5,' cm-1', &
                /5x,'Maximum energy,       EMAX =',d12.5,' cm-1', &
                /5x,'Minimum linestrength, SMIN =',d12.5,' D**2', &
                /5x,'Maximum rotations     JMAX =',i5)

         wsmi=wsmin/autocm
         wsma=wsmax/autocm
         emi=(emin+GZ)/autocm
         ema=(emax+GZ)/autocm
         smi=smin/autode**2
         open(unit=itra,form='unformatted')!,recordtype='segmented')
         open(unit=item,form='unformatted')!,recordtype='segmented')
         rewind itra
10       read(itra,end=90) j1,j2,kmin1,kmin2,neval1,neval2, &
                           idia,ipar,isym,gz,zembed,ibase1,ibase2
         maxr=maxr+neval1*neval2
         if (max(j1,j2).gt.jmax) then
            do 9 ie2=1,neval2+2
            read(itra)
9           continue
            goto 10
         endif

         allocate(ee1(neval1), ee2(neval2), ss(neval1))

         ipar1=ipar
         if(neval1.gt.nmax1) nmax1= neval1
         if(neval2.gt.nmax2) nmax2= neval2
         read(itra) ee1
         read(itra) ee2
         if(idia.eq.-2.and.mod(j1,2).eq.1.and.ipar1.eq.0) then
               ipar1=mod(ipar1+j1,2)
         else if(idia.eq.-2.and.mod(j1,2).eq.1.and.ipar1.eq.1) then
               ipar1=0
         endif
         do 1 ie2=1,neval2
         read(itra) ss
         ee=ee2(ie2)
         if (ee.lt.emi .or. ee.gt.ema) goto 1
         do 2 ie1=1,neval1

         if (ee1(ie1).lt.emi) goto 2
         if (ee1(ie1).gt.ema) goto 1
         if (ss(ie1).lt.smi) goto 2
         w= ee - ee1(ie1)
         if (abs(w).ge.wsmi .and. abs(w).le.wsma) then
           if (w.ge.0.0d0) then
             write(item) ipar1,j2,kmin2,ie2,j1,kmin1,ie1, &
                         ee2(ie2),ee1(ie1),w,ss(ie1),ibase2,ibase1
           else
             write(item) ipar1,j1,kmin1,ie1,j2,kmin2,ie2, &
                         ee1(ie1),ee2(ie2),-w,ss(ie1),ibase1,ibase2
           endif
           nr= nr + 1
         endif
  2      continue
  1      continue

         deallocate(ee1, ee2, ss)

         goto 10
90       continue
         close(unit=itra)

         write(6,1010) maxr,itra,nr,item
 1010 format(/,i10,' transition records read from unit ITRA =',i3, &
             /,i10,' selected and      written to unit ITEM =',I3)
         if (nr.le.0) then
            write(6,900) 
  900       format(/'No lines selected: stop')
            stop
         endif
      else
         rewind ispe
         read(ispe)
13       read(ispe,end=93)
         nr= nr+1
         goto 13
93       continue
         write(6,1020) nr,ispe
 1020    format(/5x,'List already sorted:', &
                /i10,' transition records found on unit ISPE =',i3) 
      endif

      call spmain(ilev,ispe,nr,nmax1,nmax2,item,idia,gz)

      stop
      end

! ------------------------------------------------------------------
     
      subroutine spmain(ilev,ispe,nr,nmax1,nmax2,item,idia,gz)

!     this is the effective main program of spectra 
      implicit double precision(a-h,o-y), logical(z)
      common/logic/ zout, zsort, zspe, zpfun, zembed
      
      data dwl/0.0d0/,x1/1.0d0/,x0/0.0d0/

!     sort out the spectrum on frequencies

      if (zsort) call sortsp(nr, item, ispe)
      if (.not.zspe) stop
      rewind ispe
      read(ispe) zembed

!     read in the nuclear spin factors


      read(5,101) ge, go
101   format(6f10.0)
      if(ge.le.x0 .and. go.le.x0) then
        ge= x1
        go= x1
      endif
      write(6,205) ge, go
205   format(/5x,'Even spin factor = ',f6.3,'  odd spin factor = ',f6.3)

!     read in the conditions for spectra required.
! 
      read(5,101, end=92) temp, xmin, wmin, wmax, dwl, q

!     print out run conditions

      write(6,206) temp
206   format(/5x,'Temperature set to: ',f8.2,' K'//)
      if (xmin.ne.x0) then
         write(6,207) xmin
      else
         write(6,208)
      endif
207   format(5x,'Minimum relative intensity required = ',f8.6)
208   format(5x,'All transitions printed out')
      if(wmax.ne.x0) then
         write(6,209) wmin, wmax, dwl
      else
         write(6,210) wmin, dwl
      endif
209   format(5x,'Frequency range from ',f10.3,' cm-1 to ',f10.3,' cm-1',/ &
             5x,'profile half width',f10.3/)
210   format(5x,'Frequency range from ',f10.3,' cm-1 to maximum',/ &
             5x,'profile half width',f10.3/)

!     calculate the partition function

      if (zpfun) then
         nlev=max(nmax1,nmax2)
         call pfcalc(temp,emax,qerr,ge,go,nlev,q,ilev,idia,gz)
         write(6,211) q, qerr
211      format(/5x,'Partition function Q =',d12.5,/, &
                5x,'estimated error in Q =',d12.5,' %')
      else
         if (q .le. x0) q=x1
         write(6,2061) q
2061    format(/5x,' Partition function set to Q =',d12.5)
      endif

!     calculate integrated absorption coefficients
      jdia= abs(idia)
      call spectm(xmin, wmin, wmax, dwl, &
                  ge, go, temp, q, nr, jdia, ispe, gz)

92    continue
       write(6,212)
212   format(/10x,'Job completed successfully')
      call timer
      stop
      end

! ------------------------------------------------------------------

      subroutine sortsp(nr, item, ispe)

!     this subroutine sorts the dipole data on ascending frequencies.
!     data printed out is in cm-1, debye**2 and sec-1.
!     data written to ispe for spectm is in atomic units.

      implicit double precision(a-h,o-y), logical(z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: da
      integer, allocatable, dimension(:,:) :: ia
      integer, allocatable, dimension(:) :: iperm


      common/logic/ zout, zsort, zspe, zpfun, zembed

!      fmem = (64*nr)/1048576.0d0
      fmem = 64.d0*dble(nr)/1048576.0d0 !Changed by Lorenzo Lodi, 10 June 2007
      write(6,300) fmem
300   format(/5x,'Memory required to sort transitions',f10.3,' MB')
      allocate(ia(nr,7), da(nr,4), iperm(nr))

      if(.not.zspe) zout= .true.
      write(ispe) zembed
      ir= 0

!     read in of data from item.
!     input is in atomic units.

      rewind item
10    read(item,end=92,err=998) ipar1, j2,kmin2,ie2, j1,kmin1,ie1, &
                        e2, e1, w, s,ibase2,ibase1

      ir= ir + 1

      ia(ir,1)= ipar1
      ia(ir,2)= j2
      ia(ir,3)= kmin2
      ia(ir,4)= ie2 + ibase2 - 1
      ia(ir,5)= j1
      ia(ir,6)= kmin1
      ia(ir,7)= ie1 + ibase1 - 1
      da(ir,1)= e2
      da(ir,2)= e1
      da(ir,3)= w
      da(ir,4)= s
      goto 10

92    continue

!     now sort on frequencies

      ier = 0

      call dpsort(da(1,3), nr, iperm, 1, ier)
      if(ier .ne. 0) then
         write(6,*) " ERROR returned by DPSORT", ier
         stop
      endif

      do 698 ir=1,7
 
      call ipperm(ia(1,ir), nr, iperm, ier)
      if(ier .ne. 0) write(6,*) " ERROR returned by ipperm ", ier
      if(ir .le. 4) then
         call dpperm(da(1,ir), nr, iperm, ier)
         if(ier .ne. 0) write(6,*) " ERROR returned by dpperm ", ier
      endif   

  698 continue

!     write out transitions for subroutine spectm to ispe.
!     output to ispe is in atomic units.

!     print out results to lineprinter if requested.
!     output to lineprinter is in wavenumbers and debye.

      if (zout) write(6,202)
202   format(//'ipar   j2  p2  i2   j1  p1  i1       e2        e1',9x, &
         'freq       s(f-i)'/)
      i20= 0
      do 4 ir= 1,nr
      write(ispe) (ia(ir,ic), ic=1,7), (da(ir,ic), ic=1,4)
      if (zout) write(6,203) (ia(ir,ic), ic=1,7), (da(ir,ic), ic=1,4)
4     continue
203   format(i3,2x,3(1x,i3),1x,3(1x,i3),3(2x,f9.3),2(3x,e10.3))

      close(unit=item)
      deallocate(ia, da, iperm)

      write(6,204)
204   format(//10x,'Sort completed successfully')
      call timer
      return
998   write(6,898) ir+1
898   format(/,5x,'Failure reading from unit item, record',i6,/)
      stop
       end

! --------------------------------------------------------------

      subroutine pfcalc(temp,emax,qerr,ge,go,nlev,q,ilev,idia,gz)


!     this is a straightforward routine to calculate partition
!     functions from supplied energy levels. the energies are
!     taken from the file supplied by programs triatom and rotlev.

!     the input stream is ilev, which is defaulted to stream 14.

      implicit double precision(a-h,o-y), logical(z)
      common/logic/ zout, zsort, zspe, zpfun, zembed

      double precision, allocatable, dimension(:) :: e
!     program constants
!     autocm converts atomic units to wavenumbers
      data x0/0.0d0/,x1/1.0d0/,&
           hc/ 1.9864476d-16/, &
           bk/ 1.3806581d-16/, & 
           autocm/ 2.19474624d+05/
      allocate(e(nlev))
      
      rewind ilev

      tk= bk*temp
      c1= hc/tk
      q= x0

!     loop over read of energies from ilev

      jmax= 0
      emax= x0
10    read(ilev,838,end=90) j, k, id, ipar, isym, neval
838   format(6i4)
      if(idia.eq.-1 .or. idia.eq.0) then
         ip= isym
      else
         if(idia.eq.2.and..not.zembed) then
            ipar= ipar + j + k + 1
            ipar= mod(ipar,2)
         endif
         ip= ipar
      endif
      write(6,202) neval,j, ip, k
202   format(i8,' levels with j=',i3,' ip=',i2,' kmin=',i2,'  included')
!     check array e is big enough, if not make it bigger...
      if (neval .gt. nlev) then
         deallocate(e)
         nlev=neval      
         allocate(e(nlev))
      endif

      read(ilev,*) (e(i), i=1,neval)
      if(j.gt.jmax) then
         emax= e(1)
         jmax= j
      endif

!     calculate the partition coefficient
!     first: degeneracy factors
      grot=dble(2*j+1)
      if(abs(idia) .eq. 2) then
         if(ip.eq.0 .or. ip.eq.2) then
            gg= ge*grot
         else
            gg= go*grot
         endif
      else
         gg=grot
      endif
      do 3 i=1,neval
      q= q + gg*exp(-(e(i)*autocm-gz)*c1)
3     continue
      goto 10
90    emax= emax*autocm-gz
      write(6,203) gz, emax
203   format(/5x,'Ground zero energy           = ',f13.5,' cm-1',/, &
               5x,'Lowest energy in jmax levels = ',f13.5,' cm-1')

      if(q.le.x0) q= x1
      qend= dble(2*jmax+1)*exp(-emax*c1)
      qerr= 100.0d0*qend/q

      close(unit=ilev)
      deallocate(e)

      return
      end

! ------------------------------------------------------------

      subroutine spectm(xmin, wmin, wmax, dwl, &
                        ge, go, temp, q, nr, jdia, ispe, gz)

!     this is the main working subroutine which calculates the
!     required integrated absorption coefficients.

!     this routine is designed to calculate integrated absorption
!     coefficients using the data produced by subroutine sortsp.
!     it uses the formula:

!           4.162034*10(-19)*gg*w*(exp(-hcw"/kt) - exp(-hcw'/kt))*s(f-i)
!     i(w)= ------------------------------------------------------------
!                                         q

!     where:
!           gg  is the nuclear spin weighting
!               times the degeneracy weighting
!           w   is the transition frequency
!           w"  is the energy of the lower level in wavenumbers
!           w'  is the energy of the upper level in wavenumbers
!           t   is the temperature at which the spectrum is required
!        s(f-i) is the line strength in debye**2 for all magnetic
!               components of the line.
!           q   is the partition function.

!     this gives the intensity in cm./molecule.

!     in many cases the partition function will not be supplied, and it
!     is thus set at unity. this must be taken into account when using
!     integrated absorption coefficients to calculate expected spectral
!     intensities at given column densities.

!     relative absorption coefficients are given relative to the maximum
!     set to unity.

!     if zemit is true, the emissivity, rather than the integrated
!     absoption coefficient, is calculated from:

!               (2j'+1)*gg*hcw*exp(-hcw'/kt)*aif
!         j(w)= --------------------------------
!                             4pi*q

      implicit double precision(a-h,o-y), logical(z)
      double precision, allocatable, dimension(:,:) :: a
      integer, allocatable, dimension(:,:) :: iqnum
      namelist /spe/ emin1,emax1,jmax,zplot,zemit,iplot,zfreq,zeinst, &
                     emin2,emax2,zprof,idat,zene,tinte,zlist,ilist, &
                     zdop,prthr,npoints,xmolm
      common/logic/ zout, zsort, zspe, zpfun, zembed
      common/base/ ibase1,ibase2
!     program constants
      data idat/19/,iplot/20/,ilist/36/, &
           hc/ 1.9864476d-16/, &
           pi/ 3.1415927d0/, & 
           bk/ 1.3806581d-16/, & 
!     autocm converts atomic units to wavenumbers
           autocm/ 2.19474624d+05/, &
!     autode converts atomic units to debye
           autode/ 2.5417662d0/, & 
!     detosc converts from s(f-i) in debye**2 to seconds**(-1)
           detosc/ 3.136186d-07/, & 
           conv  / 4.162034d-19/ 

!     detocg converts from debye to c.g.s. units.

      fmem =  76.d0*dble(nr)/1048576.0d0 !Changed by Lorenzo Lodi, 10 June 2007
      write(6,300) fmem
300   format(/5x,'Memory required to generate spectrum',f10.3,' MB'/)

      allocate(iqnum(nr,7), a(nr,6))

200   format(///)

!     preset namelist parameters

      emin1= -1.0d27
      emax1= +1.0d27
      emin2= -1.0d27
      emax2= +1.0d27
      jmax= -1
      tinte= +1.0d-15
      zplot= .false.
      zlist= .false.
      zprof= .false.
      npoints=3000
      zemit= .false.
      zfreq= .true.
      zeinst= .false.
      zene= .false.
      zdop= .false.
      prthr = 0.1d0
      xmolm = 18.0d0 !molecular mass default (H20).

!     read in namelist parameters

      read(5,spe)
      if(zeinst) zemit=.true.
      write(6,*) "SPECTM output requests:"
      if (zplot) then
         write(6,1010) iplot
 1010    format(5x,'Data for plotting written to stream IPLOT =',i3)
         open(unit=iplot,form='formatted')
         if (zprof) then
            open(unit=idat,form='formatted')
            write(6,1020) npoints,idat
 1020    format(5x,'Full line profile requested at',i5,' points,',&
                 ' scatch  IDAT =',i3)
         endif
      else
         write(6,1030)
 1030    format(5x,'Plot file not requested')
      endif
      if (zlist) then
         write(6,1040) ilist
 1040    format(5x,'Full line list    written to stream ILIST =',i3)
         open(unit=ilist,form='formatted', status='replace')
      else
         write(6,1050)
 1050    format(5x,'Full line list not requested')
      endif
      if (zlist) write(6,1060) prthr
 1060 format(/5x,'Linelist printed using relative threshold PRTHR =', &
        d9.2)

!     determine the spectral range required

      nskip= 0
      ntrans= 0
      nread= 0
      rewind ispe
      read(ispe)
10    read(ispe,end=90,err=99) ipar, j2, k2, ie2, j1, k1, ie1, e2, e1, w
      nread= nread+1
      w= w*autocm
      if(w.lt.wmin) then
         nskip= nskip + 1
         goto 10
      endif
      if(wmax.gt.0.0) then
         if(w.le.wmax) then
            ntrans= ntrans + 1
         else
            goto 90
         endif
      else
         ntrans= ntrans + 1
      endif
      goto 10
90    continue

!     skip unrequired transitions

      rewind ispe
      read(ispe)
      do 1 ir=1,nskip
       read(ispe)
1     continue

!     read in required data from ispe:

      do 2 ir=1,ntrans
      read(ispe) (iqnum(ir,ic),ic=1,7), (a(ir,ic),ic=1,4)
      a(ir,1)= (a(ir,1)*autocm)-gz
      a(ir,2)= (a(ir,2)*autocm)-gz
      a(ir,3)= a(ir,3)*autocm
      a(ir,4)= a(ir,4)*autode*autode
2     continue

!     start the calculation

      c1= hc/(temp*bk)
      if (zemit) goto 40

!     calculate integrated absorption coefficients

      c2= conv/q

     xmax= 0.0d0
      xint= 0.0d0
      if(jdia.eq.2) then
         do 3 ir=1,ntrans
         if(iqnum(ir,1).eq.0.or.iqnum(ir,1).eq.2) then
            gg= ge
         else
            gg= go
         endif
         acur= c2*gg*a(ir,3)*a(ir,4)* &
                  (exp(-a(ir,2)*c1) - exp(-a(ir,1)*c1))
         a(ir,5)= acur
         xint= xint + acur
         if(acur.gt.xmax) xmax= acur
3        continue
      else
         do 4 ir=1,ntrans
         acur= c2*a(ir,3)*a(ir,4)* &
                  (exp(-a(ir,2)*c1) - exp(-a(ir,1)*c1))
         a(ir,5)= acur
         xint= xint + acur

         if(acur.gt.xmax) xmax= acur
4        continue
      endif
      goto 41

!     calculate emissivity coefficients

40    c2= hc*detosc/(4.0d0*pi*q)
      xmax= 0.0d0
      xint= 0.0d0
      if(jdia.eq.2) then
         do 43 ir=1,ntrans
         w= a(ir,3)
         if(iqnum(ir,1).eq.0.or.iqnum(ir,1).eq.2) then
            gg= ge
         else
            gg= go
         endif
         if(.not. zeinst) then
            acur= w*w*w*w*a(ir,4)*c2*gg* &
                  exp(-a(ir,1)*c1)
         else
            acur= w*w*w*gg*a(ir,4)*detosc
         endif
         a(ir,5)= acur
         xint= xint + acur
         if(acur.gt.xmax) xmax= acur
43       continue
      else
         do 44 ir=1,ntrans
         w= a(ir,3)
         if(.not. zeinst) then
            acur= w*w*w*w*a(ir,4)*c2* &
                  exp(-a(ir,1)*c1)
         else
            acur= w*w*w*a(ir,4)*detosc
         endif
         a(ir,5)= acur
         xint= xint + acur
         if(acur.gt.xmax) xmax= acur
44       continue
      endif
41    continue

!     relative absorption or emissivity coefficients

      inc= 0
      neg= 0
      xmax1= 0.0d0
      do 5 ir= 1,ntrans
      if(jmax.gt.-1.and.iqnum(ir,5).gt.jmax) goto 5
      if(a(ir,2).lt.emin1) goto 5
      if(a(ir,2).gt.emax1) goto 5
      if(a(ir,1).lt.emin2) goto 5
      if(a(ir,1).gt.emax2) goto 5
         if(a(ir,5).gt.xmax1) xmax1= a(ir,5)
5     continue

      do 6 ir= 1,ntrans
      a(ir,6)= a(ir,5)/xmax1
      if(a(ir,6).lt.xmin) goto 7
      if(jmax.gt.-1.and.iqnum(ir,5).gt.jmax) goto 7
      if(a(ir,2).lt.emin1) goto 7
      if(a(ir,2).gt.emax1) goto 7
      if(a(ir,1).lt.emin2) goto 7
      if(a(ir,1).gt.emax2) goto 7
         inc= inc + 1
         goto 6
7     continue
         neg= neg + 1
6     continue

      if(jmax.eq.-1) then
         jm=500
      else
         jm= jmax
      endif
      if(.not.zemit) then
         write(6,2021) xmax,xint,emin1,emax1,emin2,emax2,jm,ntrans, &
                     inc,neg
      else
         write(6,2022) xmax,xint,emin1,emax1,emin2,emax2,jm,ntrans, &
                      inc,neg
      endif
2021  format(/5x,'Maximum absorption coefficient = ',e15.6,/, &
             10x,'units are cm./molecule.'//, &
             5x,'Integrated band absorption =',e15.6,//, & 
             5x,'Transitions with lower state energy = ',e13.6,' cm-1', &
             ' to ',e13.6,' cm-1 considered.',/, & 
             5x,'transitions with upper state energy = ',e13.6,' cm-1', &
             ' to ',e13.6,' cm-1 considered.',/, & 
             5x,'Maximum allowed value of J" =',i4,//, &
             i11,' transitions within spectral range',/, &
             i11,' transitions included',/, &
             i11,' transitions neglected')
2022  format(/5x,'Maximum emission coefficient = ',e15.6,/, &
             10x,'units are erg/molecule/sr.'//, &
             5x,'Integrated band emission =',e15.6,//, &
             5x,'Transitions with lower state energy = ',e13.6,' cm-1', &
             ' to ',e13.6,' cm-1 considered.',/, &
             5x,'Transitions with upper state energy = ',e13.6,' cm-1', &
             ' to ',e13.6,' cm-1 considered.',/, &
             5x,'maximum allowed value of j" =',i4,//, &
             i11,' transitions within spectral range',/, &
             i11,' transitions included',/, &
             i11,' transitions neglected')

!     final output

      write(6,*)
      i20= 0
      write(6,203)
203   format(/'ipar  j2 p2 i2  j1 p1 i1    e2         e1   ', &
       '   freq      s(f-i)    abs i(w)   rel i(w)     a(if) '/)
!204   format(i3,2x,2i3,i4,1x,2i3,i4,3f13.6,1x,e16.8,3e11.3)
204   format(i3,2x,2i5,i5,1x,2i5,i5,3f14.6,1x,es16.8,3es15.6) ! Changed by Lorenzo Lodi, 10 June 2007
      do 8 ir=1,ntrans
      if(jmax.gt.-1.and.iqnum(ir,5).gt.jmax) goto 8
      if(a(ir,2).lt.emin1) goto 8
      if(a(ir,2).gt.emax1) goto 8
      if(a(ir,1).lt.emin2) goto 8
      if(a(ir,1).gt.emax2) goto 8
      if(a(ir,6).lt.xmin) goto 8
         w= a(ir,3)
         aa= w*w*w*a(ir,4)*detosc/dble(2*iqnum(ir,2) + 1)
         iqnum(ir,3) = abs(iqnum(ir,3) - 1)
         iqnum(ir,6) = abs(iqnum(ir,6) - 1)
         if (zlist) &
         write(ilist,204) (iqnum(ir,ic),ic=1,7),(a(ir,ic),ic=1,6),aa
         if(a(ir,6).ge.prthr) then
            write(6,204) (iqnum(ir,ic),ic=1,7),(a(ir,ic),ic=1,6),aa
            i20= i20 + 1
         endif
         if(zplot.and..not.zprof) then
            if( a(ir,6) .ge.tinte) then
              if ( zfreq ) then
                  xval= a(ir,3)
               else
                  xval= 10000.0/a(ir,3)
               endif
               if ( .not.zene ) then
                  write(iplot,207) xval, a(ir,6)
               else
                  write(iplot,206) xval, a(ir,6),a(ir,1),a(ir,2), &
                  iqnum(ir,2),iqnum(ir,5)
               endif
            endif
         endif
         if(i20.eq.20) then
            write(6,200)
            i20= 0
         endif
8     continue

      if(zplot.and.zprof) then
         do 9 ir=1,ntrans
         if( zfreq ) then
            a(ir,6)= a(ir,3)
         else
            a(ir,6)= 10000.0/a(ir,3)
         endif
         write(idat,*) a(ir,6), a(ir,5)
9        continue
         fmin= a(1,6)
         fmax= a(ntrans,6)
         call profil(idat,dwl,temp,iplot,fmin,fmax,npoints,&
                         xmolm,zdop,zfreq)
      endif
206   format(1x,f11.5,2x,f13.10,2x,f12.4,2x,f12.4,2x,i3,2x,i3)
207   format(1x,f11.5,2x,f20.10)
      write(6,200)
      write(6,205)
205   format(10x,'Spectrum completed successfully')
      return
99    write(6,*) nread, j2,k2,ie2,j1,k1,ie1,e2,e1,w

      deallocate(iqnum, a)

      stop
      end

! --------------------------------------------------------

      subroutine profil(idat,dwl,temp,iplot,fmin,fmax,npoints,&
                            xmolm,zdop,zfreq)

! This routine calculates a Gaussian or thermal Doppler line profile
! and applies it to the previously calculated integrated absorption 
! coefficient or the emissivity. The frequency dependent total 
! absorption coefficient and emissivity for the molecule are 
! calculated at each frequency or wavelength point. 
!
! Output is to the unit idat.
! Units of the line profiled absorption coefficient are [cm**2 /molecule]
! Units of the line profiled emissivity are [ergs cm /molecule /str].
!   If units of [ergs /micron /molecule /str] are desired, the user
!   must multiply by (10000.0 / (lambda[microns])**2)
!
! ARGUMENTS
! idat     I/O unit, on which frequency/intensity data is stored.
! dwl      Gaussian half width as specified by the user
! temp     Temperature
! ipolt    I/O unit, on which data is to be output.
! fmin     minimum of frequency range   
! fmax     maximum of frequency range
!
! Input arguments from namelist SPE.
!
! npoints  Number of frequency points, default is 3000.
! Zdop     logical parameter, if true then a thermal Doppler
!          line half width is used. The default is false, in which
!          case the user specified value of the half width at half
!          maximum of a Gaussian is used.
! Zfreq    Default .true. units of wavenumber, and calculation 
!          proceeds as discussed below. If .false. units of
!          wavelength in microns. For consistency of intensity units
!          the line profiling is performed only in wavenumbers (cm**-1).
!          But intensity data is output as a function of wavelength.
! xmolm    Mean molecular mass of the species, default is 18.0d0 (H20).
!
! OTHER INPUT
!
!Input from unit idat
!
! w        Wavenumber (cm**-1) or wavelength (microns) of line.
! x        Integrated absorption coefficient (cm / molecule) or
!          emissivity (ergs /s /str) of line.
!
! OUTPUT to unit iplot
!
! wl(npoints)  Wavenumber (cm**-1) or wavelength (microns) of point.
! f(npoints)   Absorption cross section (cm**2 / molecule) or
!              emissivity (ergs cm /s /str) at the point.
!
! The Doppler profile is:
!
!         1        ( m*c*c   )      ( -m*c*c             )
!  DLP = --- * sqrt(---------) * exp(----------- * Dv*Dv )
!         v        ( 2*pi*k*T)      (  2*k*T*v*v         )
!
! Where: v is the central line frequency or wavenumber
!        dv is the displacement from the line centre
!        m is the mass of particle, T is temperature.
!
! The Doppler half width at half maximum is thus:
!
!               ( 2*log(2)*k*T ) 
!  HW = v * sqrt( ------------ )
!               (     m*c*c    )
!
! Note: log(2) is the natural log of 2, ie FORTRAN notation.
!
! the profile can now be written:
!
!          1        (  log(2)  )      ( -Dv*Dv*log(2) )
!  DLP = ---- * sqrt( -------- ) * exp( ------------- )
!         HW        (   pi     )      (  HW*HW        )
!
! the factor sqrt(log(2)) is often dropped from the definition of
! the Doppler half width, so simplifying the all the expressions.
!
! Here, however, we use the rigorous definition of 
! half width at half maximum of the Gaussian.

      implicit double precision (a-h,o-y), logical(z)
      double precision Temp, xmolm, xln2, rtln2, rtln2pi, xlambda 

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: f, wl

      allocate(f(npoints), wl(npoints))

      write(6,*)
      write(6,*) "LINE PROFILING"

      if(zdop) then
        write(6,*) "Thermal Doppler line profiling requested."
        write(6,*) "Molecular mass set to", xmolm
      else
        write(6,*) "Gaussian profiling requested."
        write(6,*) "Half width at half maximum set to ", dwl, " cm**-1"
      endif
      if (zfreq) then
        write(6,*) "Wavenumber units of cm**-1 selected"
      else
        write(6,*) "Wavelength units of microns selected"
      endif

      df= (fmax - fmin)/dble(npoints-1)
! some constants.
      pi= acos(-1.0d0)
      xln2=log(2.0d0)
      rtln2=sqrt(xln2)
      rtln2pi=rtln2/sqrt(pi)

!     zero the intensity array

      do 1 i= 1,npoints
      f(i)= 0.0d0
      wl(i)= fmin + dble(i-1)*df
1     continue
!     read in data from idat

      rewind idat
10    read(idat,*,end=11) w, x

      if(zfreq .eqv. .false.) then
        xlambda=w
        w=10000.0d0/w
      endif

      ! set the half width at half maximum [cm**-1]
      if(zdop) then
       hw=dop_half(temp,w,xmolm)
      else
       hw=dwl
      endif

      fac= rtln2pi/hw  ! Normalisation factor [cm]

      do 12 i=1,npoints
         ! set displacement from line centre in [cm**-1]
      if(zfreq) then
       delw= wl(i) - w
      else
       delw= (10000.0d0/wl(i))-w
      endif

      if(abs(delw).gt.10.0d0*hw) goto 12
      xarg= delw/hw
      xx= x*exp(-xarg*xarg*xln2)  ! Guassian line profile, ....
      f(i)= f(i) + (xx*fac)       ! normalise and bin. 
12    continue
      goto 10

!     output the spectrum

11    continue
      do 20 i=1,npoints

      write(iplot,*) wl(i), f(i)
20    continue

      deallocate(f, wl)

      return
      end

! ----------------------------------------------------------

      function dop_half(T,nu_cm, molm)

! This function calculates the Gaussian thermal Doppler half 
! width at half maximum.

! INPUT:
!        T:        temperature
!        nu_cm:    Wavenumber (cm-1).
!        molm:     Molecular mass
! OUTPUT:
!        dop_half: doppler half width (cm-1).

      implicit none

      double precision T, Nu_cm,molm, dop_half,  C

      data C /3.5682149d-7/

      dop_half=C*sqrt(T/Molm)*nu_cm

      return
      end

! -------------------------------------------------------------

      subroutine getrow(row,nrow,iunit)

      implicit double precision (a-h,o-y)
      dimension row(nrow)
      read(iunit,err=999) row
      return
999   write(6,899) iunit
899   format(/,5x,'fell over in getrow reading unit = ',i6,/)
      stop
      end

! ------------------------------------------------------------

      subroutine outrow(row,nrow,iunit)

      implicit double precision (a-h,o-y)
      dimension row(nrow)
      write(iunit) row
      return
      end
      subroutine timer
!     prints current cpu time usage                                 #030

      implicit double precision (a-h,o-y)
      common/timing/itime0
      write(6,10)
      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/)
      return
10    format(/)
      end


