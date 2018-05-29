!     dummy main program                                           #001
      call rotlev3z
      stop
      end

!#######################################################################
      subroutine rotlev3z

!     program               r o t l e v 3 z
!
!     program to do rotational analysis using vibrational output of
!     dvr3drjz run with idia = -2, zrot = .true., zperp=.false.
!
!     written (originally) by Max Kostin, 2001
!     re-written (so that it works) by Paolo , 2005/2007
!     see:
!     M.A.Kostin, O.L.Polyansky & J.Tennyson, J. Chem. Phys., 116, 7564 (2002)
!
!     use as follows:
!     comments on namelist parameters (& defaults) in block data
!     the program needs the following subroutines:
!     1. limited card input which is read in subroutine insize
!        and files on stream ivec, ivec2 (& ivec1, ivec3) from DVR3DRJZ
!     2. f02fjf to do iterative diagonalisation (nag routine).
!     the program works in **** atomic units ***** 

      implicit double precision (a-h,o-y), logical (z)

      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                namelist/prt/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec & 
                    ,nploti,nplotf,ithre

      write(6,1000)
 1000 format(5x,'Program ROTLEV3Z (version of May 2002)')
!     read in namelist input data (defaults in block data)
      read(5,prt)

!     read in control parameters of problem.
      call insize

!     first for select, then onto main program
      call select 

      stop
      end

!#######################################################################
      block data
!     stores defaults for namelist parameters                       #003
      implicit double precision (a-h,o-y), logical (z)

!     outp holds information which controls the amount of printed output
!     zpvec: print eigenvectors if zpvec = .true.
!     thresh: threshold for printing a coefficient if zpvec=.true.
!     stream         holds                              used if
!      ilev    		input/output of eigenvalues              zpfun=.true.
!      ivec,ivec1      input  eigenvalues & eigenvectors        always
!      ivec2,ivec3    input  eigenvalues & eigenvectors        nktot .gt. 2
!      iscr,iscr_1,iscr-diag,
!      iscr_2,iscr2_1,iscr_f  hamiltonian files                  always


      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
          
      data thresh/0.1d0/,zpham/.false./,zpvec/.false./,idiag/2/,&
           ivec/26/,zvec/.false./,zplot/.false./,jvec/3/,iscr/1/,&
           ivec1/27/,zpfun/.false./,ilev/14/,kvec/8/,kvecpb/9/,&
           ztran/.false./,zptra/.false./,zcut/.false./,nploti/1/,nplotf/0/,&
           ithre/-8/
      end

!########################################################################
      subroutine insize

!     set up common /size/ & write control parameters of problem    #004

      implicit double precision (a-h,o-y), logical (z)

!     common /size/ stores control parameters for the problem
!     nbass: maximum dimension of rotational secular problem
!     ibass: actual dimension of rotational secular problem
!     mbass: maximum size of vibrational problem (excluding linear geom)
!     mbass0: maximum size of vibrational problem (including linear geom)
!     maxblk_even: size of vibrational radial problem (even basis)
!     maxblk_odd: size of vibrational radial problem (odd  basis)
!     ndvr : maximum dimension of theta dvr grid used in vibrational problem
!     iang : maximum number of discrete angles retained in vib. problem
!     npnt : number of gauss-associated legendre grid points requested
!     jrot : total rotational angular momentum
!     if kpar=0 calculate only even K
!          iqpar=0, K=0 even; iqpar=1  K=0 odd
!     if kpar=1 calculate only odd K
!          if qpar=0 q=0("+"), if qpar=1 q=1("-")
!     if iqpar=0 K=0 odd
!     nr : number of dvr points in each radial coordinate
!     meval: number of eigenvalues computed in the vibrational problem
!     nvib : number of vibrational eigenvalues used in rotational prob.
!     neval: number of eigenvalues of interest
!     nktot : number of k values

!      integer, ALLOCATABLE, DIMENSION(:) :: inda1
!      integer, ALLOCATABLE, DIMENSION(:) :: inda2
!      integer, ALLOCATABLE, DIMENSION(:) :: indb1
!      integer, ALLOCATABLE, DIMENSION(:) :: indb2

      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb      
     character(len=8) title(9)

      read(5,5)  nvib,neval,kpar,ibass,iqpar,npnt
    5 format(6i5)
      open(unit=ivec,form='unformatted')
      open(unit=ivec1,form='unformatted')

      rewind ivec
      rewind ivec1

      read(ivec)  idia_1,ipar_1,ndvr_1,nr_1,nr_1,jrot_1,kmin0_1,& 
           meval_1,nlim_1
      read(ivec1) idia_2,ipar_2,ndvr_2,nr_2,nr_2,jrot_2,kmin0_2,&
           meval_2,nlim_2
      nskip=7
      call skipblock(nskip,ivec)
      call skipblock(nskip,ivec1)
      read(ivec) iang_1,ibass_1
      read(ivec1)iang_2,ibass_2
      read(ivec) 
      read(ivec1)
      read(ivec) meval_1
      read(ivec1)meval_2

!     check input data for consistency
      ierr26=abs(ipar_1)
      ierr27=abs(ipar_2-1)
      if (ierr26.gt.0) then
         write(6,850) ivec,ipar
         stop
      endif
      if (ierr27.gt.0) then
        write(6,850) ivec1,ipar_1
        stop
      endif
 850  format(5x,'ERROR: unit =',i3,2x,'has wrong ipar =',i2)
     
!check consistency on all params

      if (iang_1.ne.iang_2) then
         write(6,1200)' iang ',iang_1,iang_2
         stop
      else 
         iang=iang_1
      end if
      if (idia_1.ne.idia_2) then
         write(6,1200)' idia ',idia_1,idia_2
         stop
      else 
         idia=idia_1
      end if
      if (ndvr_1.ne.ndvr_2) then
         write(6,1200)' ndvr ',ndvr_1,ndvr_2
         stop
      else 
         ndvr=ndvr_1
      end if
      if (nr_1.ne.nr_2) then
         write(6,1200)' nr ',nr_1,nr_2
         stop
      else 
         nr=nr_1
      end if
      if (jrot_1.ne.jrot_2) then
         write(6,1200)' jrot ',jrot_1,jrot_2
         stop
      else 
         jrot=jrot_1
      end if
      if (kmin0_1.ne.kmin0_2) then
         write(6,1200)' kmin0',kmin0_1,kmin0_2
         stop
      else 
         kmin0=kmin0_1
      end if
      if (nlim_1.ne.nlim_2) then
         write(6,1200)' nlim ',nlim_1,nlim_2
         stop
      else 
         nlim=nlim_1
      end if

 1200 format(/5x,'Problems with variable',A6,&
             /3x,'stream 26: ',I6,&
             /3x,'stream 27: ',I6,&
             /3x,'i stop here. ')

      if (npnt.le.ndvr) then
         write(6,*)'Npnt set too small...'
         npnt=max(npnt,ndvr)
         write(6,*)'resizing it to ',npnt
      end if



      maxblk_even=nr*(nr+1)/2
      maxblk_odd=nr*(nr-1)/2

      write(6,1120) jrot
 1120 format(/5x,'J =',i3,' rotational state')

      if (kpar .eq. 0) then
        write(6,1)
    1   format(5x,'Calculate EVEN K block')
        if (iqpar .eq. 0) then
           write(6,3)  
    3      format(5x,'Block K=0 contains EVEN parity functions')
        else
           write(6,4)
    4      format(5x,'Block K=0 contains ODD parity functions')
        endif
      else
        write(6,2)
    2   format(5x,'Calculate ODD K block')
        if (iqpar .eq. 0) then
            write(6,33)
            write(6,34)
    33     format(/10x,'Hamiltonian matrix has structure: positive K->q=0')
    34     format(10x,'                                : negative K->q=1')
        else
            write(6,44)
            write(6,45)
    44     format(/10x,'Hamiltonian matrix has structure: positive K->q=1')
    45     format(10x,'                                : negative K->q=0')
        endif
      endif
!
!
!    define k-blocks to be taken into account
      if (kpar.eq.0) then
         kmina=0
         kminb=2
         if (iqpar.eq.0) then
            ipa=ivec
            ipb=ivec1
            mevala=meval_1
            mevalb=meval_2
            ibassa=ibass_1
            ibassb=ibass_2
            iqa=0
            iqb=1
            isa=0
            isb=1
         else
            ipa=ivec1
            ipb=ivec
            mevala=meval_2
            mevalb=meval_1
            ibassa=ibass_2
            ibassb=ibass_1
            iqa=0
            iqb=1
            isa=1
            isb=0
         end if
      else
         kmina=1
         kminb=1
         ipa=ivec
         ipb=ivec1
         mevala=meval_1
         mevalb=meval_2
         ibassa=ibass_1
         ibassb=ibass_2
         isa=0
         isb=1
         if (iqpar.eq.0) then
            iqa=0
            iqb=1
         else
            iqb=0
            iqa=1
         end if
      end if
      nskipka=7+mevala
      nskipkb=7+mevalb
!
!   indicate which k-blocks
! the general structure for file ivec is: 0 1e 1o 2 3 4 5 .....  J
!  which is                               1 2  3  4 5 6 7 ..... J-2    
! define the k-block index as outlined above:
      indk(1)=0
      indk(2)=1
      indk(3)=1
      do k=4,100
         indk(k)=k-2
      end do
!----------------------------------------------
      if (kpar.eq.0) then !even k-blocks
! we build the hamiltonian as follows
!  0e 2e 2o 4e 4o 6e 6o ......
!  1  2  3  4  5  6  7
         inda1(1)=1
         inda2(1)=1
         j=1
         ia=1
         ib=0
         do k=2,jrot,2
         ib=ib+1
         ia=ia+1
         inda1(ia)=k+2
         inda2(ia)=j+1
         indb1(ib)=k+2
         indb2(ib)=j+2
         j=j+2
      end do
      nka=ia
      nkb=ib
      nktot=j
      else if (kpar.eq.1) then
! we build the hamiltonian as follows
!  1a 1b 3a 3b 5a 5b ......
!  1  2  3  4  5  6 
! which k=1 block do we need?
         if (iqpar.eq.0) then
         inda1(1)=2
         indb1(1)=3
         else
         inda1(1)=3
         indb1(1)=2
         end if
         inda2(1)=1
         indb2(1)=2
! now select all remainings odd k-blocks
         i=1
         j=2
         do k=3,jrot,2
            i=i+1
         inda1(i)=k+2
         inda2(i)=j+1
         indb1(i)=k+2
         indb2(i)=j+2
         j=j+2
         end do
      nka=i
      nkb=i
      nktot=j
      end if

      write(6,*)'Stream A assigned to unit ',ipa
      write(6,*)'Stream B assigned to unit ',ipb
      write(6,*)'Number of k-blocks from stream A ...',nka
      write(6,*)'Number of k-blocks from stream B ...',nkb
      write(6,*)'Total number of k-blocks ........',nktot

!     compute size of rotational secular problem

!      if (nvib.ne.0) nvib=0
      if (nvib.eq.0) nvib=max(mevala,mevalb)
      nviba=min(nvib,mevala)
      nvibb=min(nvib,mevalb)


      write(6,1000) meval_1,meval_2,iang*maxblk_even,iang*maxblk_odd,ndvr,nvib,npnt,neval,nbass
 1000 format(/5x,'Rotational part of rot-vib calculation  with:',&
             /i9,3x,'lowest vibrational eigenvectors supplied from stream 26',&
             /i9,3x,'lowest vibrational eigenvectors supplied from stream 27',&
             /i9,3x,'dimension even vibration secular problem,',&
             /i9,3x,'dimension odd vibration secular problem, with',&
             /i9,3x,'angular dvr points,',&
             /i9,3x,'lowest vibrational eigenvectors actually used',&
             /i9,3x,'point gauss-associated legendre integration',&
             /i9,3x,'lowest rotational eigenvectors required for',&
             /i9,3x,'dimension rotation secular problem')
      if (ibass .gt. 0) write(6,1005)
 1005 format(17x,'with basis selected by energy ordering')

      read(5,500)   title
  500 format(9a8)
      write(6,1010) title
 1010 format(/5x,'title: ',9a8)
      write(6,1013)
 1013 format(/5x,'Diagonalisation performed in core using lapack',&
                  ' routine dsyev')
      if (zpvec) write(6,1040) thresh
 1040 format(5x,'Printing of eigenvector coefficients greater than',&
                 ' thresh =',f5.2,' requested')
      if (.not.zpvec) write(6,1050)
 1050 format(5x,'Printing of eigenvectors not requested')
      IF (ZTRAN) THEN
          ZVEC = .TRUE.
          IF (ZPTRA) WRITE(6,1052)
 1052     FORMAT(5X,'PRINTING OF TRANSFORMED VECTORS REQUESTED')
          IF (.NOT.ZPTRA) WRITE(6,1051)
 1051     FORMAT(5X,'PRINTING OF TRANSFORMED VECTORS NOT REQUESTED')
      ENDIF
      write(6,1053) iscr
 1053 format(/5x,'Hamiltonian scratch file to be held on stream ',&
             'iscr  =',i4)
      if (zpfun) write(6,1055) ilev
 1055 format( 5x,'Eigenvalues    to be written to end of stream ',&
             'ilev  =',i4)
      if (zvec) write(6,1054) jvec
 1054 format( 5x,'Eigenvalues & vectors to be written to stream ',&
             'jvec  =',i4)
      if (ztran) write(6,1058) kvec
 1058 format( 5x,'Transformed vectors   to be written to stream ',&
             'kvec  =',i4)
      return
      end
!#######################################################################

      subroutine select
 
!     subroutine select determines which vibrational basis          #007
!     functions are to be used
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
      common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb      
      real*8, allocatable :: eviba(:,:),evibb(:,:),mviba(:),mvibb(:)
      real*8, allocatable :: ea(:),eb(:)
      integer, allocatable :: iva(:),ivb(:)

      integer, DIMENSION(NKTOT) :: IV

      parameter (conv=219474.631d0)

      character(len=4) symm(2)
      data symm/'even','odd '/
!     read energies from file ivec,
!        first skip matrix elements


      write(6,*)'proceed building the off-diag k-blocks...'

      nikb=max(nkb, 1)  ! take care of the case with no B k-block 
      allocate(eviba(mevala,nka),evibb(mevalb,nikb))

      write(6,*)'energy analysis ....'
      write(6,*)'Notation :   k:qs '

! 1. read energies from a
      ipr=ipa
!    skip header  N=5
      rewind(ipr)
      nskip=5
      call skipblock(nskip,ipr)
      j=0
      do ia=1,nka
! skip first j-1 blocks
         i=inda1(ia)
         jp=inda2(ia)
         k2=indk(i)
         do isk=j+1,i-1
      call skipblock(nskipka,ipr)
         end do
! get j-th block      
      read(ipr)                             !6
      read(ipr)                             !7
      read(ipr)                             !8
      read(ipr)                             !9
      read(ipr)                             !10
      read(ipr)                             !11
      read(ipr) (eviba(ii,ia),ii=1,mevala)  !12
      do ii=1,mevala
      read(ipr)                             !13   
      end do
      write(6,190)k2,iqa,isa,eviba(1,ia)*conv,eviba(mevala,ia)*conv
! move 1 to  2     
      j=i
      end do

! 1. read energies from b
      ipr=ipb
!    skip header  N=5
      rewind(ipr)
      nskip=5
      call skipblock(nskip,ipr)
      j=0
      do ib=1,nkb
! skip first j-1 blocks
         i=indb1(ib)
         jp=indb2(ib)
         k2=indk(i)
         do isk=j+1,i-1
      call skipblock(nskipkb,ipr)
         end do
! get j-th block      
      read(ipr)                             !6
      read(ipr)                             !7
      read(ipr)                             !8
      read(ipr)                             !9
      read(ipr)                             !10
      read(ipr)                             !11
      read(ipr) (evibb(ii,ib),ii=1,mevalb)  !12
      do ii=1,mevalb
      read(ipr)                             !13   
      end do
      write(6,190)k2,iqb,isb,evibb(1,ib)*conv,evibb(mevalb,ib)*conv
! move 1 to  2     
      j=i
      end do
190   format('(',I3,':',2I1,'): min ',d20.12,' max ',d20.12)
!--------------------------------------------
!#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
! energy selection of 3d states 

      allocate(iva(nka),ivb(nikb),ea(nka),eb(nikb))
      iva=0
      ivb=0
      nsum=0
      eviba=eviba*conv
      evibb=evibb*conv

      emin=eviba(1,1)
      do ia=1,nka
         emin=min(emin,eviba(1,ia))
      end do
      do ib=1,nkb
         emin=min(emin,evibb(1,ib))
      end do

      emax=eviba(nviba,1)
      do ia=1,nka
         emax=max(emax,eviba(nviba,1))
      end do
      do ib=1,nkb
         emax=max(emax,evibb(nvibb,1))
      end do

      read(5,*)ecut

      nbass=nka*nviba+nkb*nvibb

      if (zcut) then
         write(6,*)'Vib levels selected according to energy cut:'
         write(6,*)' E_cut = ',ecut,'  cm-1 .'
         if (ecut.ge.emax) then
            iva=mevala
            ivb=mevalb
            write(6,*)emax
196 format('Ecut bigger than maximum evalue... ',f20.4)
         else
! case a: zcut=true, selection thus based on energy cut-off
      do ia=1,nka
         ic=0
         do l=1,nviba
            if (eviba(l,ia).gt.ecut.and.ic.eq.0) then
               ic=1
               iva(ia)=l-1
            end if
         end do
         if (ic.eq.0) then
            iva(ia)=nviba
         end if
               nsum=nsum+iva(ia)
      end do
      do ib=1,nkb
         ic=0
         do l=1,nvibb
            if (evibb(l,ib).gt.ecut.and.ic.eq.0) then
               ic=1
               ivb(ib)=l-1
            end if
         end do
         if (ic.eq.0) then
            ivb(ib)=nvibb
         end if
               nsum=nsum+ivb(ib)
      end do
      
      nsum=0
      do ia=1,nka
         nsum=nsum+iva(ia)
      end do
      do ib=1,nkb
         nsum=nsum+ivb(ib)
      end do
      nbass=nsum

      end if
      else
! case b: choosing the lowest nbass states ...... (zcut=false)

         if (ibass.eq.0) then

            write(6,*)'Selecting ALL 3d functions....'
            iva=nviba
            ivb=nvibb

            else 

               nbass=min(nbass,ibass)
               iva=0.d0
               ivb=0.d0

         do ia=1,nka
         ea(ia)=eviba(1,ia)
         end do
         do ib=1,nkb
         eb(ib)=evibb(1,ib)
         end do

      do jl=1,nbass
         
         emina=ea(1)
         i=1
         do ia=2,nka
         emina=min(ea(ia),emina)
         if (emina.eq.ea(ia)) i=ia
         end do

         eminb=eb(1)
         j=1
         do ib=2,nkb
         eminb=min(eb(ib),eminb)
         if (eminb.eq.eb(ib)) j=ib
         end do


         if (emina.lt.eminb) then
            iva(i)=iva(i)+1
            if (iva(i).eq.nviba) then
               ea(i)=1.d300
            else
               ea(i)=eviba(iva(i)+1,i)
            end if
         else 
            ivb(j)=ivb(j)+1
            if (ivb(j).eq.nvibb) then 
               eb(j)=1.d300
            else 
               eb(j)=evibb(ivb(j)+1,j)
            end if
         end if
      end do

      nsum=0
      do ia=1,nka
         nsum=nsum+iva(ia)
      end do
      do ib=1,nkb
         nsum=nsum+ivb(ib)
      end do
      if (nsum.ne.nbass) then
         write(6,*)'Problem in state selection (zcut=F):'
         write(6,*)'Nbass = ',nbass
         write(6,*)'but checksum returned ... ',nsum
      end if
      end if
      end if

      emax=eviba(iva(1),1)
      do ia=1,nka
         if (iva(ia).ne.0) emax=max(emax,eviba(iva(ia),ia))
      end do
      do ib=1,nkb
         if (ivb(ib).ne.0) emax=max(emax,evibb(ivb(ib),ib))
      end do

      write(6,1000) nbass,emin,emax
 1000 format(/i10,' Functions selected from E =',d20.10,' to',d20.10,&
             ' cm-1 ')

      if (nbass.eq.0) then
         write(6,*)'No functions selected...'
         stop
      end if

      open(unit=idiag,form='unformatted')
      do ia=1,nka
         write(idiag)inda2(ia),iva(ia)
         if (iva(ia).gt.0) then
         do i=1,iva(ia)
            write(idiag)eviba(i,ia)
         end do
         end if
      end do
      do ib=1,nkb
         write(idiag)indb2(ib),ivb(ib)
         if (ivb(ib).gt.0) then
         do i=1,ivb(ib)
            write(idiag)evibb(i,ib)
         end do
         end if
      end do
      close(unit=idiag)

      if (neval.eq.0) neval=10
      neval=min(neval,nbass)
 
!     determine  and print basis set labels
 
      write(6,1010)
 1010 format(//5x,' basis functions selected')
      do ia=1,nka
      kz=indk(inda1(ia))
      if (iva(ia).lt.mevala) then
      write(6,1021)kz,iqa,isa,iva(ia),mevala
      else
      write(6,1022)kz,iqa,isa,iva(ia),mevala
      end if
      end do
      do ib=1,nkb
      kz=indk(indb1(ib))
      if (ivb(ib).lt.mevalb) then
      write(6,1021)kz,iqb,isb,ivb(ib),mevalb
      else 
      write(6,1022)kz,iqb,isb,ivb(ib),mevalb
      end if
      end do  
      
 1020 format(5x,'k =',i3,', i runs from',i5,' to',i5,2x,a4)
 1021 format(5x,'(',I3,':',2I1,'):',i5,'/',i5)
 1022 format(5x,'(',I3,':',2I1,'):',i5,'/',i5,' *SATURATED****')
  310 continue

      call vrmain(iva,ivb)

      return
!     error on unit ivec
  900 rewind ivec
      i=1
  920 read(ivec,end=930)
      i=i+1
      goto 920
  930 write(6,940) ivec,i
  940 format(/'   end of stream ivec =',i4,' at record',i6)
      stop
      end

!#####################################################################
      subroutine vrmain(iva,ivb)

!     subroutine vrmain is the 'real' main program & contains       #006
!     the calls to the various subroutines which set & solve hamil
!
      implicit double precision (a-h,o-y), logical (z)

      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
      integer, DIMENSION(nka) :: iva
      integer, DIMENSION(nkb) :: ivb
      integer, DIMENSION(nktot) :: mvib
      real*8, ALLOCATABLE, DIMENSION(:) :: radmee,radmoo,radmeo
      data x0/0.0d0/

      ksize=0
      do ia=1,nka
         ksize=ksize+iva(ia)
      end do
      do ib=1,nkb
         ksize=ksize+ivb(ib)
      end do
      write(6,11) ksize
   11 format(/5x,'Hamiltonian size =',i6)
      neval=min(neval,ksize)

      maxblk_odd=max(1,maxblk_odd)
      maxblk_even=max(1,maxblk_even)

      ALLOCATE(radmee(maxblk_even),radmoo(maxblk_odd),radmeo(maxblk_odd))
!  2d symmetrised radial matrix elements from 1d componants
      call radint(radmee,radmoo,radmeo)

!  set up the hamiltonian matrix (not for restart runs)
      call solrt2(radmee,radmoo,radmeo,iva,ivb,noffblk)
      write(6,1050)
1050  format(/5x,'hamiltonian construction complete')
!      call timer
      DEALLOCATE(radmee,radmoo,radmeo)
 
      ezero=x0
      read(5,505,end=555) ezero
  505 format(f20.0)
  555 continue

 ! compact iva nad ivb into a single array:

      mvib=0
      do ia=1,nka
         mvib(inda2(ia))=iva(ia)
      end do
      do ib=1,nkb
         mvib(indb2(ib))=ivb(ib)
      end do

      call dgrot(1,ezero,ksize,noffblk,mvib)

      write(6,1060)
 1060 format(/5x,'diagonalisation complete')

      return
      end

!#######################################################################
      subroutine radint(radmee,radmoo,radmeo)
 
!     subroutine radint calculates the two-dimensional radial basis
!     functions between two symmetrised orthogonal coordinates.
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb

      real*8, DIMENSION(nr) :: rmb
      real*8, DIMENSION(nr) :: rma
      real*8, DIMENSION(maxblk_even) :: radmee
      real*8, DIMENSION(maxblk_odd) :: radmoo,radmeo
 
!     first read 1d radial matrix elements from file
      rewind(ipa)
      nskip=3
      call skipblock(nskip,ipa)
      read(ipa) (rma(i),i=1,nr)

      rewind(ipb)
      nskip=3
      call skipblock(nskip,ipb)
      read(ipb) (rmb(i),i=1,nr)

!     then use this data to construct the radial matrices
      call mkrad(radmee,nr,rma,0,0,0)
      if (maxblk_odd .gt. 0) then
         call mkrad(radmoo,nr,rma,1,0,0)
         call mkrad(radmeo,nr,rma,1,1,0)
      endif
 
      return
      end

!##############################################################
      subroutine mkrad(radmat,nr,rm2,iq,ip,is)
 
!     construct symmetrised matrix elements with the appropriate parity
 
      implicit double precision (a-h,o-y), logical (z)
      dimension radmat(*),rm2(*)
 
      symm  = dble(1-2*ip)
      sign  = dble(1-2*is)
 
      ipt=0
      do i1=1,nr
         do i2=1,i1-iq
            ipt=ipt+1
            radmat(ipt) = sign * (rm2(i2) + symm*rm2(i1))
         end do
      end do
      return
    end subroutine mkrad

!#######################################################################
      subroutine angin_pl_pl(angmat,pleg1,pleg2,k1,k2,&
                         iv1,iv2,nang1,nang2,angfac)

!     angin_pl_pl calculates the angular integral between blocks k and k+2 
!     for ++ block
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                   ndvr,iang,npnt,maxblk_odd,ibass,&
                   nktot,kpar,iqpar

      real*8, DIMENSION(nang1,nang1) :: pleg1
      real*8, DIMENSION(nang2,nang2) :: pleg2
      integer, DIMENSION(nang1) :: iv1
      integer, DIMENSION(nang2) :: iv2
      real*8, DIMENSION(nang1,nang2) :: fbrmat
      real*8, DIMENSION(nang1,nang2) :: angmat
      real*8, DIMENSION(npnt) :: jxcos,jwalf
      real*8, DIMENSION(nang1,npnt) :: plega
      real*8, DIMENSION(nang2,npnt) :: plegb

      data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/,x2/2.0d0/
 
!     first: set up an npnt gauss-associated legendre quadrature scheme
      angmat=0.d0
      realk1 = dble(k1)
      realk2 = dble(k2)
      npnt2 = (npnt+1)/2
      realj = DBLE(jrot)

      alf = dsqrt(xp5*(realj**2+realj-realk1**2))-x1
      bet = alf

      alft1 = dsqrt(xp5*(realj**2+realj-realk1**2))
      bett1 = alft1

      alft2 = dsqrt(xp5*(realj**2+realj-realk2**2))
      bett2 = alft2

      alfw = (xp5*(alft1+alft2))-x1
      betw = alfw

      call gaujac(jxcos,jwalf,2*npnt2,alfw,betw)

      write(6,1000) npnt,alfw,(jxcos(i),jwalf(i),i=1,npnt2)
1000  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

!     evaluate the polynomials at the quadrature points
         CALL jac_basis(npnt,nang1-1,alft1,bett1,jxcos,plega)
         CALL jac_basis(npnt,nang2-1,alft2,bett2,jxcos,plegb)

!     compute the fbr matrix elements
         fbrmat = x0
         do 40 j=1,nang1
            do 50 k=1,nang2
               sum=x0
               do 30 i=1,npnt
                  term = jwalf(i)*jxcos(i) 
        sum = sum + term * plega(j,i) * plegb(k,i)
30             continue 
        fbrmat(j,k) = sum
50          continue
40       continue

!     final step: transform the fbr matrix elements to the dvr 
      ij=0
      do 60 jp= 1,nang1
         if (iv1(jp) .eq. 0) goto 60
         ij=ij+1
         ik=0
         do 65 kp=1,nang2
            if (iv2(kp) .eq. 0) goto 65
            ik=ik+1
            sum=x0
            do 70 j =1,nang1
               do 75 k=1,nang2
                  sum = sum + pleg1(j,jp) * pleg2(k,kp) * fbrmat(j,k)
75             continue
70          continue
            angmat(ij,ik) = sum * angfac
65       continue
60    continue

      return
      end

!#######################################################################
      subroutine angin_xdg(angmat,pleg1,pleg2,k1,k2,&
                         iv1,iv2,nang1,nang2,angfac)

!     angin_pl_pl calculates the angular integral between blocks k and k+2 
!     for ++ block
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                   ndvr,iang,npnt,maxblk_odd,ibass,&
                   nktot,kpar,iqpar

      real*8, DIMENSION(nang1,nang1) :: pleg1
      real*8, DIMENSION(nang2,nang2) :: pleg2
      integer, DIMENSION(nang1) :: iv1
      integer, DIMENSION(nang2) :: iv2
      real*8, DIMENSION(nang1,nang2) :: fbrmat
      real*8, DIMENSION(nang1,nang2) :: angmat
      real*8, DIMENSION(npnt) :: jxcos,jwalf
      real*8, DIMENSION(nang1,npnt) :: plega
      real*8, DIMENSION(nang2,npnt) :: plegb

      data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/,x2/2.0d0/
 
!     first: set up an npnt gauss-associated legendre quadrature scheme
      angmat=0.d0
      realk1 = dble(k1)
      realk2 = dble(k2)
      npnt2 = (npnt+1)/2
      realj = DBLE(jrot)

      alf = dsqrt(xp5*(realj**2+realj-realk1**2))-x1
      bet = alf

      alft1 = dsqrt(xp5*(realj**2+realj-realk1**2))
      bett1 = alft1

      alft2 = dsqrt(xp5*(realj**2+realj-realk2**2))
      bett2 = alft2

      alfw = (xp5*(alft1+alft2))-x1
      betw = alfw

      call gaujac(jxcos,jwalf,2*npnt2,alfw,betw)

      write(6,1000) npnt,alfw,(jxcos(i),jwalf(i),i=1,npnt2)
1000  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

!     evaluate the polynomials at the quadrature points
         CALL jac_basis(npnt,nang1-1,alft1,bett1,jxcos,plega)
         CALL jac_basis(npnt,nang2-1,alft2,bett2,jxcos,plegb)

!     compute the fbr matrix elements
         fbrmat = x0
         do 40 j=1,nang1
            do 50 k=1,nang2
               sum=x0
               do 30 i=1,npnt
                  term = jwalf(i)*jxcos(i)
                  sum = sum + term * plega(j,i) * plegb(k,i)
30             continue 
        fbrmat(j,k) = sum
50          continue
40       continue

!        stop

!     final step: transform the fbr matrix elements to the dvr 
      ij=0
      do 60 jp= 1,nang1
         if (iv1(jp) .eq. 0) goto 60
         ij=ij+1
         ik=0
         do 65 kp=1,nang2
            if (iv2(kp) .eq. 0) goto 65
            ik=ik+1
            sum=x0
            do 70 j =1,nang1
               do 75 k=1,nang2
                  sum = sum + pleg1(j,jp) * pleg2(k,kp) * fbrmat(j,k)
75             continue
70          continue
            angmat(ij,ik) = sum * angfac
65       continue
60    continue

!      alfw = (xp5*(alft1+alft2))
!      betw = alfw
!      call gaujac(jxcos,jwalf,nang1,alfw,betw)
!      angmat=0.d0
!     do i=1,nang1
!        angmat(i,i)=jxcos(i)*angfac/(1.d0-jxcos(i)**2)
!     end do

!     do i=1,nang1
!     do j=1,nang1
!      write(76,*)i,j,angmat(i,j)
!     end do
!     end do

      return
      end 

!#######################################################################
      subroutine angin_pl_pl_2(angmat,pleg1,pleg2,k1,k2,&
                         iv1,iv2,nang1,nang2,angfac)

!     angin_pl_pl calculates the angular integral between blocks k and k+2 
!     for ++ block
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                   ndvr,iang,npnt,maxblk_odd,ibass,&
                   nktot,kpar,iqpar

      real*8, DIMENSION(nang1,nang1) :: pleg1
      real*8, DIMENSION(nang2,nang2) :: pleg2
      integer, DIMENSION(nang1) :: iv1
      integer, DIMENSION(nang2) :: iv2
      real*8, DIMENSION(nang1,nang2) :: fbrmat
      real*8, DIMENSION(nang1,nang2) :: angmat
      real*8, DIMENSION(npnt) :: jxcos,jwalf
      real*8, DIMENSION(nang1,npnt) :: plega
      real*8, DIMENSION(nang2,npnt) :: plegb

      data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/,x2/2.0d0/
 
!     first: set up an npnt gauss-associated legendre quadrature scheme
      angmat=0.d0
      realk1 = dble(k1)
      realk2 = dble(k2)
      npnt2 = (npnt+1)/2
      realj = DBLE(jrot)

      alf = dsqrt(xp5*(realj**2+realj-realk1**2))
      bet = alf

      alft1 = dsqrt(xp5*(realj**2+realj-realk1**2))
      bett1 = alft1

      alft2 = dsqrt(xp5*(realj**2+realj-realk2**2))
      bett2 = alft2

      alfw = (xp5*(alft1+alft2))
      betw = alfw

      call gaujac(jxcos,jwalf,2*npnt2,alfw,betw)

      write(6,1000) npnt,alfw,(jxcos(i),jwalf(i),i=1,npnt2)
1000  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

!     evaluate the polynomials at the quadrature points
         CALL jac_basis(npnt,nang1-1,alft1,bett1,jxcos,plega)
         CALL jac_basis(npnt,nang2-1,alft2,bett2,jxcos,plegb)

!     compute the fbr matrix elements
         fbrmat = x0
         do 40 j=1,nang1
            do 50 k=1,nang2
               sum=x0
               do 30 i=1,npnt
                  term = jwalf(i)*jxcos(i)/(1.d0-(jxcos(i)**2)) 
        sum = sum + term * plega(j,i) * plegb(k,i)
30             continue 
        fbrmat(j,k) = sum
50          continue
40       continue

!     final step: transform the fbr matrix elements to the dvr 
      ij=0
      do 60 jp= 1,nang1
         if (iv1(jp) .eq. 0) goto 60
         ij=ij+1
         ik=0
         do 65 kp=1,nang2
            if (iv2(kp) .eq. 0) goto 65
            ik=ik+1
            sum=x0
            do 70 j =1,nang1
               do 75 k=1,nang2
                  sum = sum + pleg1(j,jp) * pleg2(k,kp) * fbrmat(j,k)
75             continue
70          continue
            angmat(ij,ik) = sum * angfac
65       continue
60    continue

!            do i=1,nang1
!               do j=1,nang2
!                  write(72,*)i,j,angmat(i,j)
!               end do
!            end do
        angmat=0.d0
       do i=1,npnt
           angmat(i,i) = angfac * jxcos(i)/(1.d0-(jxcos(i)**2)) 
        end do

      return
      end

!########################################################################
!#######################################################################
      subroutine angin_pl_mn(angmat,pleg1,pleg2,k1,k2,&
                         iv1,iv2,nang1,nang2,angfac)

!     angin_pl_pl calculates the angular integral between blocks k and k+2 
!     for ++ block
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                   ndvr,iang,npnt,maxblk_odd,ibass,&
                   nktot,kpar,iqpar

      real*8, DIMENSION(nang1,nang1) :: pleg1
      real*8, DIMENSION(nang2,nang2) :: pleg2
      integer, DIMENSION(nang1) :: iv1 
      integer, DIMENSION(nang2) :: iv2 
      real*8, DIMENSION(nang1,nang2) :: fbrmat
      real*8, DIMENSION(nang1,nang2) :: angmat
      real*8, DIMENSION(npnt) :: jxcos,jwalf
      real*8, DIMENSION(nang1,npnt) :: plega
      real*8, DIMENSION(nang2,npnt) :: plegb

      data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/,x2/2.0d0/
 
!     first: set up an npnt gauss-associated legendre quadrature scheme
      angmat=0.d0
      realk1 = dble(k1)
      realk2 = dble(k2)
      realj = DBLE(jrot)
      npnt2 = (npnt+1)/2

      alft1 = dsqrt(xp5*(realj**2+realj-realk1**2))
      bett1 = alft1

      alft2 = dsqrt(xp5*(realj**2+realj-realk2**2))
      bett2 = alft2

      alfw = (xp5*(alft1+alft2))-xp5
      betw = alfw

      call gaujac(jxcos,jwalf,npnt,alfw,betw)

      write(6,1000) npnt,alfw,(jxcos(i),jwalf(i),i=1,npnt2)
1000  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

!     evaluate the polynomials at the quadrature points
         CALL jac_basis(npnt,nang1-1,alft1,bett1,jxcos,plega)
         CALL jac_basis(npnt,nang2-1,alft2,bett2,jxcos,plegb)

!     compute the fbr matrix elements
         fbrmat = x0
         do 40 j=1,nang1
            do 50 k=1,nang2
               sum=x0
               do 30 i=1,npnt
                  term = - jwalf(i)
        sum = sum + term * plega(j,i) * plegb(k,i)
30             continue 
        fbrmat(j,k) = sum
50          continue
40       continue

!     final step: transform the fbr matrix elements to the dvr 
      ij=0
      do 60 jp= 1,nang1
         if (iv1(jp) .eq. 0) goto 60
         ij=ij+1
         ik=0
         do 65 kp=1,nang2
            if (iv2(kp) .eq. 0) goto 65
            ik=ik+1
            sum=x0
            do 70 j =1,nang1
               do 75 k=1,nang2
                  sum = sum + pleg1(j,jp) * pleg2(k,kp) * fbrmat(j,k)
75             continue
70          continue
            angmat(ij,ik) = sum * angfac
65       continue
60    continue

      return
      end

!#################################################################
      subroutine angin_cor(angmat,k1,iv1,iv2,nang1,nang2,angxt,iq)
 
!     calculates the angular integral between blocks k and k
!     required for the coriolis coupling term for orthogonal (radau)
!     coordinates in the bisector embedding. 
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar

      DIMENSION iv1(ndvr),iv2(ndvr)
      real*8, DIMENSION(ndvr,ndvr) :: fbrmat
      real*8, DIMENSION(iang,iang) :: angmat
      real*8, DIMENSION(npnt) :: jxcos,jwalf,jxcosw,jwalfw
      real*8, DIMENSION(nang1) :: norm
      real*8, DIMENSION(0:1) :: xsign
      real*8, DIMENSION(nang1,npnt) :: plega,plegw

      data x0/0.0d0/,xp5/0.5d0/,x4/4.0d0/
!     evaluate the polynomials at the quadrature points
      angmat=0.d0
      realk = dble(k1)
      realj = DBLE(jrot)
      npnt2 = (npnt+1)/2
      xsign(0)=1.d0
      xsign(1)=-1.d0

! extra-term for k=1
      xx=0.d0
      if (k1.eq.1) xx = angxt
      xy=-xsign(iq)*xp5*dble(k1)

      angmat=x0

      alf = dsqrt(xp5*(realj**2+realj-realk**2))
      bet = alf
      alft = alf-xp5
      bett = alft

      call gaujac(jxcos,jwalf,npnt,alft,bett)
      call gaujac(jxcosw,jwalfw,nang2,alf,bet)

      write(6,1001) npnt,alft,(jxcos(i),jwalf(i),i=1,npnt2)
1001  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

      CALL jac_basis(npnt,nang2-1,alf,bet,jxcos,plega)
      CALL jac_basis(npnt,nang2-1,alf,bet,jxcosw,plegw)
      CALL norms2(norm,nang2-1,alf,bet)
!     now compute the I1 integral
      fbrmat = x0
      fbrmatc = x0
         do 40 j=1,nang1
            do 50 k=1,nang2
               sum=0.d0
               sumc=0.d0
               rk=dble(k-1)
               if (k.gt.1) then
                  fact3=norm(k)/norm(k-1)
                  do i=1,npnt
                     fact1 = - jxcos(i) * (1.d0 + 2.d0 * alf + 2.d0 * rk)* xy 
                     fact2 = 2.d0 * (alf + rk) * fact3 * xy
                     sum  = sum + jwalf(i) * plega(j,i) * & 
                          (( fact1 + xx) * plega(k,i) + fact2 * plega(k-1,i))
                     sumc  = sumc + jwalf(i) * plega(j,i) * & 
                          (fact1 * plega(k,i) + fact2 * plega(k-1,i))
                  end do
               else
                  do i=1,npnt
                     fact1 = - jxcos(i) * (1.d0 + 2.d0 * alf + 2.d0 * rk)* xy 
                     sum  = sum + jwalf(i) * plega(j,i)*(( fact1 + xx) * plega(k,i) )
                     sumc  = sumc + jwalf(i) * plega(j,i)*(fact1 * plega(k,i) )
                  end do
               end if
              fbrmat(j,k)=  sum
50          continue
40       continue

      do j=1,ndvr
         cw=sqrt(jwalfw(j))
         do i=1,ndvr
            plegw(i,j)=plegw(i,j)*cw
         enddo
      enddo

!     second step: transform the fbr matrix elements to the dvr
      i1=0
      do 60 i= 1,nang1
      if (iv1(i) .eq. 0) goto 60
      i1=i1+1
      i2=0
      do 65 ip=1,nang2
      if (iv2(ip) .eq. 0) goto 65
      i2=i2+1
      sum=x0
      do j =1,nang1
      do jp=1,nang2
      sum = sum + plegw(j,i) * plegw(jp,ip)*fbrmat(j,jp)
      end do
      end do
         angmat(i1,i2)=sum
   65 continue
   60 continue
      return
      end

!#################################################################

      subroutine wrtho(diag,offdg,mvib,nbass,nktot)
!     print hamiltonian matrix                                      #009
      implicit double precision (a-h,o-y)
      dimension diag(nbass),offdg(*),mvib(nktot)
      write(6,1010) diag
 1010 format('1',5x,'hamiltonian matrix: diagonal elements',&
             /(10f13.8))
 
      ioff1=1
      ioff2=mvib(1)*mvib(2)
      num=1
      write(6,1030) num,(offdg(i),i=ioff1,ioff2)
 1030 format(//5x,'off-diagonal block number',i3/(10f13.8))
      do 10 k=3,nktot
      ioff1=ioff2+1
      ioff2=ioff2+mvib(k-2)*mvib(k)
      num=num+1
      write(6,1030) num,(offdg(i),i=ioff1,ioff2)
      ioff1=ioff2+1
      ioff2=ioff2+mvib(k-1)*mvib(k)
      num=num+1
      write(6,1030) num,(offdg(i),i=ioff1,ioff2)
   10 continue
      return
      end
      subroutine wrthi(hamil,nham)
!     print hamiltonian matrix in core version
      double precision hamil(nham,nham)
      write(6,1010)
 1010 format(5x,'hamiltonian matrix'/)
      do 30 i=1,nham
      write(6,1020) (hamil(i,j),j=1,i)
 1020 format(10f13.7)
   30 continue
      return
      end

!#######################################################################
!#######################################################################
      subroutine loadh(mvib,hamil,noffblk,ksize,kss)
 
!     subroutine loadh loads the hamiltonian matrix from disk
 
      implicit double precision (a-h,o-y), logical (z)
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
      integer,  DIMENSION(nktot) :: mvib
      integer,  DIMENSION(nktot) :: n0
      real*8, dimension(kss) :: hamil
      real*8, dimension(ksize) :: diagt
      real*8, ALLOCATABLE, DIMENSION(:) :: diag
      real*8, ALLOCATABLE, DIMENSION(:) :: offdg
      data x0/0.0d0/

      hamil=x0
      cm=219474.631d0

      rewind iscr
!      rewind ipa
!      rewind ipb

         n0(1)=1
      do k=2,nktot
         n0(k)=n0(k-1)+mvib(k-1)
      end do

!###########  pozitive K #############
! write diagonal element in hamiltonian 

      write(6,*)'Diagonal elements.....'
      open(unit=idiag,form='unformatted')

      do k=1,nka
         k1=inda1(k)
         k2=indk(k1)
         ip=inda2(k)
         read(idiag)inu,mvv
         do i=1,mvv
            iki=n0(ip)-1+i
            read(idiag)diagt(iki)
         end do
      end do
      do k=1,nkb
         k1=indb1(k)
         k2=indk(k1)
         ip=indb2(k)
         read(idiag)inu,mvv
         do i=1,mvv
            iki=n0(ip)-1+i
            read(idiag)diagt(iki)
         end do
      end do
! skip headers
!      nskip=5
!      call skipblock(nskip,ipa)
!      call skipblock(nskip,ipb)

! first load elements from stream A
!      write(6,*)'Stream A  ....'
!      k0=1
!      do k=1,nka
!         k1=inda1(k)
!         k2=indk(k1)
!         ip=inda2(k)
!         mvib0=mvib(ip)
!         allocate (diag(mvib0))
!! skip from k0 to k1-1        
!         write(6,*)k,': loading block  ',k1,' or k=',k2
!         do i=k0,k1-1
!      call skipblock(nskipka,ipa)            
!         end do
!! load k1 block         
!         call load_diag(diag,mvib0,ipa)
!! put diag into hamil
!         i0=n0(ip)-1
!         do ii=1,mvib0
!            iki=(i0+ii)
!            ik=iki+(iki-1)*(2*ksize-iki)/2
!            hamil(ik)=diag(ii)
!         end do
!         deallocate(diag)
!! set next skipping from k1+1
!         k0=k1+1
!      end do
!! repeat for stream B
!      write(6,*)'Stream B  ....'
!      k0=1
!      do k=1,nkb
!         k1=indb1(k)
!         k2=indk(k1)
!         ip=indb2(k)
!         mvib0=mvib(ip)
!         allocate (diag(mvib0))
!         write(6,*)k,': loading block  ',k1,' or k=',k2
!         do i=k0,k1-1
!      call skipblock(nskipkb,ipb)            
!         end do
!         call load_diag(diag,mvib0,ipb)
!         i0=n0(ip)-1
!         do ii=1,mvib0
!            iki=(i0+ii)
!            ik=iki+(iki-1)*(2*ksize-iki)/2
!            hamil(ik)=diag(ii)
!!            write(6,*)iki,iki*(iki+1)/2,ik,hamil(ik)
!         end do
!         deallocate(diag)
!         k0=k1+1
!      end do



      do i=1,ksize
            ik=i+(i-1)*(2*ksize-i)/2
            hamil(ik)=diagt(i)/cm
!            write(73,*)i,hamil(ik),diagt(i)/cm
      end do

      write(6,*)'off-diag elements ........'

      do ii=1,noffblk
         read(iscr)i,j
         write(6,*)'loading offdiag block ....',ii,i,j
         mvib1=mvib(i)
         mvib2=mvib(j)
         allocate(offdg(mvib1*mvib2))
         read(iscr)offdg
         if (i.eq.j) then
         ll=0
         i0=n0(i)-1
         j0=n0(j)-1
         do ip=1,mvib1
            do jp=1,ip
               ll=(ip-1)*mvib1+jp
               iki=i0+ip
               ikj=j0+jp
               ik=iki+(ikj-1)*(2*ksize-ikj)/2
               hamil(ik)=offdg(ll)+hamil(ik)
!               write(73,*)ip,jp,ik,hamil(ik),offdg(ll)
            end do
         end do
         else
         ll=0
         i0=n0(i)-1
         j0=n0(j)-1
         do ip=1,mvib1
            do jp=1,mvib2
               ll=ll+1
               iki=i0+ip
               ikj=j0+jp
               ik=iki+(ikj-1)*(2*ksize-ikj)/2
               hamil(ik)=offdg(ll)
 !              write(73,*)ip,jp,ik,hamil(ik),offdg(ll)
            end do
         end do
         end if
         deallocate(offdg)
      end do

     return
   end subroutine loadh

!#########################################################################
!#####################################################################
!#########################################################################
     subroutine load_diag(diag,nk,ipr)

      implicit double precision (a-h,o-y), logical (z)

      real*8, DIMENSION(nk) :: diag

! skip k-block header (size=5)
      call skipblock(5,ipr)
! read k-block functions 
      read(ipr)meval
! read lowest nk eigenvalues
      read(ipr)(diag(i),i=1,nk)
! skip all eigenvectors
      call skipblock(meval,ipr)

      return
    end subroutine load_diag 
!#########################################################################
!#####################################################################
!#########################################################################
!######################################################################
!======================================================================
!######################################################################
      subroutine solrt2(radmee,radmoo,radmeo,iva,ivb,noffblk)
 
!     subroutine solrt2 sets up non-zero parts of hamiltonian        #010
!     including the computation of the k dependent angular matrix
!     elements
 
      implicit double precision (a-h,o-y), logical (z)
 
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
      common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
      integer, DIMENSION(nka) :: iva
      integer, DIMENSION(nkb) :: ivb
      real*8, DIMENSION(maxblk_even) :: radmee
      real*8, DIMENSION(maxblk_odd) :: radmoo,radmeo
      real*8, DIMENSION(0:1) :: xsign
      real*8, DIMENSION(iang,iang) :: angmat
      real*8, ALLOCATABLE, DIMENSION(:,:) :: pleg1,pleg2
      real*8, ALLOCATABLE, DIMENSION(:,:) :: plega1,plega2
      real*8, ALLOCATABLE, DIMENSION(:,:) :: plegb1,plegb2
      real*8, ALLOCATABLE, DIMENSION(:,:) :: coef1,coef2
      real*8, ALLOCATABLE, DIMENSION(:,:) :: coefb1,coefb2
      real*8, ALLOCATABLE, DIMENSION(:,:) :: coefa1,coefa2
      real*8, ALLOCATABLE, DIMENSION(:) :: offdg
      integer, ALLOCATABLE, DIMENSION(:) :: iva1,iva2
      integer, ALLOCATABLE, DIMENSION(:) :: iv1,iv2
      integer, ALLOCATABLE, DIMENSION(:) :: ivb1,ivb2


      data x0/0.0d0/,sqrt2/1.4142135623731d0/,xp5/0.5d0/,x4/4.0d0/


!     nktot - number of K

      open(unit=iscr,form='unformatted',status='scratch')
      rewind iscr

      rewind ipa
      rewind ipb
         
      jjp1 = jrot * (jrot+1)
      idpt = 1
      iangsq=iang*iang
      noffblk=0
      xsign(1)=-1.d0
      xsign(0)=1.d0

      mvib0=max(iva(1),ivb(1))
      do ia=2,nka
         mvib0=max(mvib0,iva(ia))
      end do
      do ib=2,nkb
         mvib0=max(mvib0,ivb(ib))
      end do


      write(6,*)'Off-diag k-blocks.....'

! study of the hamiltonian
! k=0 needs special care
!
      write(6,*)'Verbose simulation building of the ham'
      write(6,*)'Searching for non-zero k-overlaps'
      write(6,*)'Notation :   k:qs '




      write(6,*)nka,nkb,nbassa,nbassb
      write(6,*)(inda1(i),i=1,nka)
      write(6,*)(inda2(i),i=1,nka)
      write(6,*)(indb1(i),i=1,nkb)
      write(6,*)(indb2(i),i=1,nkb)

      do ia1=1,nka
         k1=indk(inda1(ia1))
         i1=inda2(ia1)
         iq1=iqa
         is1=isa
         mvib1=iva(ia1)
         do ia2=1,ia1
         k2=indk(inda1(ia2))
         i2=inda2(ia2)
         iq2=iqa
         is2=isa
         mvib2=iva(ia2)
         mn=mvib1*mvib2
         if (abs(k1-k2).eq.2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',iva(ia1),iva(ia2)
         else if (k1.eq.k2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,' diag',iva(ia1),iva(ia2)
         if (k1.eq.1) & 
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  xdg',iva(ia1),iva(ia2)
         end if
         end do
      end do

      do ib1=1,nkb
         k1=indk(inda1(ib1))
         i1=indb2(ib1)
         iq1=iqb
         is1=isb
         mvib1=ivb(ib1)
         do ib2=1,ib1
         k2=indk(indb1(ib2))
         i2=indb2(ib2)
         iq2=iqb
         is2=isb
         mvib2=ivb(ib2)
         mn=mvib1*mvib2
         if (abs(k1-k2).eq.2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',ivb(ib1),ivb(ib2)
         else if (k1.eq.k2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,' diag',ivb(ib1),ivb(ib2)
         if (k1.eq.1) & 
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  xdg',ivb(ib1),ivb(ib2)
         end if
         end do
      end do

      do ia=1,nka
         k1=indk(inda1(ia))
         i1=inda2(ia)
         iq1=iqa
         is1=isa
         mvib1=iva(ia)
         do ib=1,nkb
         k2=indk(inda1(ib))
         i2=indb2(ib)
         iq2=iqb
         is2=isb
         mvib2=ivb(ib)
         mn=mvib1*mvib2
!   coriolis type
         if (k1.eq.k2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  cor',iva(ia),ivb(ib)
         end if
!   rot type
         if (abs(k1-k2).eq.2.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',iva(ia),ivb(ib)
         end if
! xtra rotational term for k1=k2=1
         if (k1.eq.1.and.k2.eq.1.and.mn.ne.0) then
         write(6,1010)k1,iq1,is1,k2,iq2,is2,' xrot',iva(ia),ivb(ib)
         end if

         end do
      end do

1010  format('(',I3,':',2I1,') <-> (',I3,':',2I1,'),   type ...',A5,'. Size:',I5,' x ',I5)

!      goto 901

! effective building of the hamiltonian

      write(6,*)'proceed building the off-diag k-blocks...'
!
! 1. <a|a>

      ipr=ipa
      iq1=iqa
      iq2=iqa
      is1=isa
      is2=isa
      id1=0
      id2=0
      id=isa
      mvl=mevala
      allocate(pleg1(iang,iang),pleg2(iang,iang),iv1(iang),iv2(iang))
      pleg1=0.d0
      pleg2=0.d0
!  due to the particular structure of file ivec, cases k=0,1
!  need to be treated separately
!
!    skip header  N=5
      rewind(ipr)
      nskip=5
      call skipblock(nskip,ipr)
! skip first j-1 blocks
         j=inda1(1)
         jp=inda2(1)
         k2=indk(j)
         mvib2=iva(1)
      allocate(coef2(ibassa,mvib2))
         do isk=1,j-1
      call skipblock(nskipka,ipr)
         end do
! get j-th block      
      read(ipr)                                 !6
      read(ipr)                                 !7
      call getrow(pleg2,iangsq,ipr)             !8 
      read(ipr)                                 !9
      read(ipr) (iv2(ii),ii=1,iang)             !10
      read(ipr)                                 !11
      read(ipr)                                 !12
      call rdcoef(coef2,ibassa,mvib2,ipr)       !13a
      do ii=mvib2+1,mvl
      read(ipr)                                 !13b   
      end do
!-----------------------------------------------------------
      if (k2.eq.1) then !add extra diagonal term ....
         mn = mvib2*mvib2
         if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k2,iq2,is2,k2,iq2,is2,'  xdg',mvib2,mvib2
      angfac = -xsign(iqa)*dble(jjp1)/x4
!      write(6,*)angfac,iq1,iq2,is1,is2
      ALLOCATE(offdg(mn))
      offdg=0.d0
         call angin_xdg(angmat,pleg2,pleg2,k2,k2,&
                          iv2,iv2,iang,iang,angfac)
         if (isa.eq.0) then
         call solofd(mn,radmee,angmat,nr,id,id2,id2,offdg,iang,&
                     iang,iang,mvib2,mvib2,coef2,coef2,&
                     ibassa,ibassa)
         else
         call solofd(mn,radmoo,angmat,nr,id,id2,id2,offdg,iang,&
                     iang,iang,mvib2,mvib2,coef2,coef2,&
                     ibassa,ibassa)
         end if
! output
      noffblk=noffblk+1
         write(iscr)jp,jp
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
   end if
!-----------------------------------------------------------
      do k=2,nka
         i=inda1(k)
         ip=inda2(k)
         k1=indk(i)
         mvib1=iva(k)
      allocate(coef1(ibassa,mvib1))
         do isk=j+1,i-1
      call skipblock(nskipka,ipr)
         end do
! get i-th block      
      read(ipr)                                 !6
      read(ipr)                                 !7
      call getrow(pleg1,iangsq,ipr)             !8 
      read(ipr)                                 !9
      read(ipr) (iv1(ii),ii=1,iang)             !10
      read(ipr)                                 !11
      read(ipr)                                 !12
      call rdcoef(coef1,ibassa,mvib1,ipr)       !13a
      do ii=mvib1+1,mevala
      read(ipr)                                 !13b   
      end do

         mn = mvib1*mvib2
         if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',mvib1,mvib2
         angfac = -sqrt(dble((jjp1-k2*(k2+1))*(jjp1-k2*(k2+3)-2)))/x4
         if (k1 .eq. 0 .or. k2 .eq. 0) angfac =  angfac*sqrt2
         ALLOCATE(offdg(mn))
         offdg=x0
         call angin_pl_pl(angmat,pleg1,pleg2,k1,k2,&
                          iv1,iv2,iang,iang,angfac)
         if (is1.eq.0) then
         call solofd(mn,radmee,angmat,nr,id,id1,id2,offdg,iang,&
                     iang,iang,mvib1,mvib2,coef1,coef2,&
                     ibassa,ibassa)
         else
         call solofd(mn,radmoo,angmat,nr,id,id2,id2,offdg,iang,&
                     iang,iang,mvib1,mvib2,coef1,coef2,&
                     ibassa,ibassa)
         end if
! output
      noffblk=noffblk+1
      write(iscr)ip,jp
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
! move 1 to  2     
      j=i
      jp=ip
      k2=k1
      iv2=iv1
      pleg2=pleg1
      coef0=0.d0
      mvib2=mvib1
      deallocate(coef2)
      allocate(coef2(ibassa,mvib2))
      coef2=coef1
      deallocate(coef1)
      end do
      deallocate(coef2)
      deallocate(pleg1,pleg2,iv1,iv2)
      write(6,*)'<a|a> terminated!!'

! 1. <b|b>

      ipr=ipb
      iq1=iqb
      iq2=iqb
      is1=isb
      is2=isb
      id1=0
      id2=0
      id=isb
      mvl=mevalb
      allocate(pleg1(iang,iang),pleg2(iang,iang),iv1(iang),iv2(iang))
      pleg1=0.d0
      pleg2=0.d0
!  due to the particular structure of file ivec, cases k=0,1
!  need to be treated separately
!
!    skip header  N=5
      rewind(ipr)
      nskip=5
      call skipblock(nskip,ipr)
! skip first j-1 blocks
         j=indb1(1)
         jp=indb2(1)
         k2=indk(j)
         mvib2=ivb(1)
         allocate(coef2(ibassb,mvib2))
         do isk=1,j-1
      call skipblock(nskipkb,ipr)
         end do
! get j-th block      
      read(ipr)                             !6
      read(ipr)                             !7
      call getrow(pleg2,iangsq,ipr)         !8 
      read(ipr)                             !9
      read(ipr) (iv2(ii),ii=1,iang)         !10
      read(ipr)                             !11
      read(ipr)                             !12
      call rdcoef(coef2,ibassb,mvib2,ipr)   !13a
      do ii=mvib2+1,mvl
      read(ipr)                             !13b   
      end do
!-----------------------------------------------------------
      if (k2.eq.1) then !add extra diagonal term ....
         mn = mvib2*mvib2
         if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k2,iq2,is2,k2,iq2,is2,'  xdg',mvib2,mvib2
         angfac = -xsign(iqb)*dble(jjp1)/x4
      ALLOCATE(offdg(mn))
         call angin_xdg(angmat,pleg2,pleg2,k2,k2,&
                          iv2,iv2,iang,iang,angfac)
         if (isb.eq.0) then
         call solofd(mn,radmee,angmat,nr,id,id2,id2,offdg,iang,&
                     iang,iang,mvib2,mvib2,coef2,coef2,&
                     ibassb,ibassb)
         else
         call solofd(mn,radmoo,angmat,nr,id,id2,id2,offdg,iang,&
                     iang,iang,mvib2,mvib2,coef2,coef2,&
                     ibassb,ibassb)
         end if
! output
      noffblk=noffblk+1
         write(iscr)jp,jp
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
   end if
!-----------------------------------------------------------
      do k=2,nkb
         i=indb1(k)
         ip=indb2(k)
         k1=indk(i)
         mvib1=ivb(k)
         allocate(coef1(ibassb,mvib1))
         do isk=j+1,i-1
      call skipblock(nskipkb,ipr)
         end do
! get i-th block      
      read(ipr)                             !6
      read(ipr)                             !7
      call getrow(pleg1,iangsq,ipr)         !8 
      read(ipr)                             !9
      read(ipr) (iv1(ii),ii=1,iang)         !10
      read(ipr)                             !11
      read(ipr)                             !12
      call rdcoef(coef1,ibassb,mvib1,ipr)   !13a
      do ii=mvib1+1,mevalb
      read(ipr)                             !13b   
      end do

         mn = mvib1*mvib2
         if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',mvib1,mvib2
         angfac = -sqrt(dble((jjp1-k2*(k2+1))*(jjp1-k2*(k2+3)-2)))/x4
         if (k1 .eq. 0 .or. k2 .eq. 0) angfac =  angfac*sqrt2
      ALLOCATE(offdg(mn))
         call angin_pl_pl(angmat,pleg1,pleg2,k1,k2,&
                          iv1,iv2,iang,iang,angfac)
         if (is1.eq.0) then
         call solofd(mn,radmee,angmat,nr,id,id2,id1,offdg,iang,&
                     iang,iang,mvib1,mvib2,coef1,coef2,&
                     ibassb,ibassb)
         else
         call solofd(mn,radmoo,angmat,nr,id,id1,id2,offdg,iang,&
                     iang,iang,mvib1,mvib2,coef1,coef2,&
                     ibassb,ibassb)
         end if
! output
      noffblk=noffblk+1
         write(iscr)ip,jp
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
! move 1 to  2      
      jp=ip
      j=i
      k2=k1
      iv2=iv1
      pleg2=pleg1
      mvib2=mvib1
      deallocate(coef2)
      allocate(coef2(ibassb,mvib2))
      coef2=coef1
      deallocate(coef1)
      end do
      deallocate(coef2)
      deallocate(pleg1,pleg2,iv1,iv2)
      write(6,*)'<b|b> terminated!!'
!#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
! more complicate part!! : <a|b>

      iq1=iqa
      iq2=iqb
      is1=isa
      is2=isb
      ida=1-isa
      idb=1-isb
      id=1
      allocate(plega1(iang,iang),plegb1(iang,iang))
      allocate(plega2(iang,iang),plegb2(iang,iang))
      allocate(iva1(iang),ivb1(iang),iva2(iang),ivb2(iang))
      iva1=0.d0
      iva2=0.d0
      ivb1=0.d0
      ivb2=0.d0
      plega1=0.d0
      plega2=0.d0
      plegb1=0.d0
      plegb2=0.d0

! skip headers again
      nskip=5
      rewind(ipa)
      call skipblock(nskip,ipa)
      rewind(ipb)
      call skipblock(nskip,ipb)
! load first k-block for each stream
      ia=1
      ib=1

      i=inda1(ia)
      ip1=inda2(ia)
      k1=indk(i)
      mvib1=iva(ia)
      allocate(coefa1(ibassa,mvib1))
         do isk=1,i-1
      call skipblock(nskipka,ipa)
         end do
! get i-th block      
      read(ipa)                             !6
      read(ipa)                             !7
      call getrow(plega1,iangsq,ipa)        !8 
      read(ipa)                             !9
      read(ipa) (iva1(ii),ii=1,iang)        !10
      read(ipa)                             !11
      read(ipa)                             !12
      call rdcoef(coefa1,ibassa,mvib1,ipa)  !13a
      do ii=mvib1+1,mevala
      read(ipa)                             !13b   
      end do
      i0=i

      j=indb1(ib)
      jp1=indb2(ib)
      k2=indk(j)
      mvib2=ivb(ib)
      allocate(coefb1(ibassb,mvib2))
         do isk=1,j-1
      call skipblock(nskipkb,ipb)
         end do
! get i-th block      
      read(ipb)                             !6
      read(ipb)                             !7
      call getrow(plegb1,iangsq,ipb)        !8 
      read(ipb)                             !9
      read(ipb) (ivb1(ii),ii=1,iang)        !10
      read(ipb)                             !11
      read(ipb)                             !12
      call rdcoef(coefb1,ibassb,mvib2,ipb)  !13a
      do ii=mvib2+1,mevalb
      read(ipb)                             !13b   
      end do
      j0=j

! do extra off-k block if kmina.ne.kminb
      if (k1.ne.k2) then
         mn = mvib1*mvib2  
         if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',mvib1,mvib2
      write(6,*)jp1,ip1
         angfac = sqrt(dble((jjp1-k1*(k1+1))*(jjp1-k1*(k1+3)-2)))/x4
         angfac = xsign(iqb) * angfac
         if (k1 .eq. 0 .or. k2 .eq. 0) angfac =  angfac*sqrt2
      ALLOCATE(offdg(mn))
         call angin_pl_mn(angmat,plegb1,plega1,k2,k1,&
                          ivb1,iva1,iang,iang,angfac)
         call solofd(mn,radmeo,angmat,nr,id,idb,ida,offdg,iang,&
                     iang,iang,mvib2,mvib1,coefb1,coefa1,&
                     ibassb,ibassa)
         write(iscr)jp1,ip1
      noffblk=noffblk+1
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
!  make k1 even to k2      
      ia=ia+1
      i0=i
      i=inda1(ia)
      ip1=inda2(ia)
      k1=indk(i)
      mvib1=iva(ia)
      deallocate(coefa1)
      allocate(coefa1(ibassa,mvib1))
         do isk=i0+1,i-1
      call skipblock(nskipka,ipa)
         end do
! get i-th block      
      read(ipa)                             !6
      read(ipa)                             !7
      call getrow(plega1,iangsq,ipa)        !8 
      read(ipa)                             !9
      read(ipa) (iva1(ii),ii=1,iang)        !10
      read(ipa)                             !11
      read(ipa)                             !12
      call rdcoef(coefa1,ibassa,mvib1,ipa)  !13a
      do ii=mvib1+1,mevala
      read(ipa)                             !13b   
      end do
      end if

      i0=i
      ka1=k1
      kb1=k2
      mviba1=mvib1
      mvibb1=mvib2
      
! calculate first coriolis term on its own
      k1=ka1
      k2=kb1
      mn=mvib1*mvib2
      if (mn.ne.0) then
      ALLOCATE(offdg(mn))
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  cor',mvib1,mvib2
      write(6,*)jp1,ip1
!         coefb1=1.d0
!         coefa1=1.d0
      angfac= dble(jjp1)/x4
         call angin_cor(angmat,k1,ivb1,iva1,iang,iang,angfac,iqb)
         call solofd(mn,radmeo,angmat,nr,id,idb,ida,offdg,iang,&
                     iang,iang,mvib2,mvib1,coefb1,coefa1,&
                     ibassb,ibassa)
         noffblk=noffblk+1
         write(iscr)jp1,ip1
!         write(37,*)'----------------------------------'
!         write(37,*)k1,k2,jjp1
!         write(37,*)ibassa,ibassb,nviba,nvibb
!         write(37,*)jp1,ip1
!         write(37,*)((angmat(ld,jd),ld=1,2),jd=1,2)
!         write(37,*)'coefa1 coefb1  '
!         do ld=1,mvib0
!         write(37,*)ld,coefb1(100,ld),coefa1(100,ld)
!         end do
!         write(37,*)'++++++++++++++++++++++++++++++++++'
!         write(37,*)(offdg(ld),ld=1,mn)
!         write(37,*)'----------------------------------'
         call outrow(offdg,mn,iscr)
         deallocate(offdg)
         end if
! now start the big loop over k

      do jjj=2,nkb

! load two extra k-blocks
! from stream A
      ia=ia+1
      i=inda1(ia)
      ip2=inda2(ia)
      ka2=indk(i)
      mviba2=iva(ia)
      allocate(coefa2(ibassa,mviba2))
         do isk=i0+1,i-1
      call skipblock(nskipka,ipa)
         end do
! get i-th block      
      read(ipa)                             !6
      read(ipa)                             !7
      call getrow(plega2,iangsq,ipa)        !8 
      read(ipa)                             !9
      read(ipa) (iva2(ii),ii=1,iang)        !10
      read(ipa)                             !11
      read(ipa)                             !12
      call rdcoef(coefa2,ibassa,mviba2,ipa) !13a
      do ii=mviba2+1,mevala
      read(ipa)                             !13b   
      end do
      i0=i
! from stream B
      ib=ib+1
      j=indb1(ib)
      jp2=indb2(ib)
      kb2=indk(j)
      mvibb2=ivb(ib)
      allocate(coefb2(ibassb,mvibb2))
         do isk=j0+1,j-1
      call skipblock(nskipkb,ipb)
         end do
! get j-th block      
      read(ipb)                             !6
      read(ipb)                             !7
      call getrow(plegb2,iangsq,ipb)        !8 
      read(ipb)                             !9
      read(ipb) (ivb2(ii),ii=1,iang)        !10
      read(ipb)                             !11
      read(ipb)                             !12
      call rdcoef(coefb2,ibassb,mvibb2,ipb) !13a
      do ii=mvibb2+1,mevalb
      read(ipb)                             !13b   
      end do
      j0=j

!  k11 vs k22  (rot)

      k1=ka1
      k2=kb2
      mvib1=mviba1
      mvib2=mvibb2
      mn = mvib1*mvib2
      if (mn.ne.0) then
      angfac = sqrt(dble((jjp1-k1*(k1+1))*(jjp1-k1*(k1+3)-2)))/x4
      angfac = xsign(iqb) * angfac
      if (k1 .eq. 0 .or. k2 .eq. 0) angfac =  angfac*sqrt2 
      ALLOCATE(offdg(mn))
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',mvib1,mvib2
      write(6,*)jp2,ip1
         call angin_pl_mn(angmat,plegb2,plega1,k2,k1,&
                         ivb2,iva1,iang,iang,angfac)
         call solofd(mn,radmeo,angmat,nr,id,idb,ida,offdg,iang,&
                     iang,iang,mvib2,mvib1,coefb2,coefa1,&
                     ibassb,ibassa)
      noffblk=noffblk+1
         write(iscr)jp2,ip1
         call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
!  k21 vs k12  (rot)

      k1=kb1
      k2=ka2
      mvib1=mvibb1
      mvib2=mviba2
      mn = mvib1*mvib2
      if (mn.ne.0) then
      ALLOCATE(offdg(mn))
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  rot',mvib1,mvib2
      write(6,*)ip2,jp1
      angfac = sqrt(dble((jjp1-k1*(k1+1))*(jjp1-k1*(k1+3)-2)))/x4
      angfac = xsign(iqa) * angfac
      if (k1 .eq. 0 .or. k2 .eq. 0) angfac =  angfac*sqrt2 
      call angin_pl_mn(angmat,plega2,plegb1,k2,k1,&
           iva2,ivb1,iang,iang,angfac)
      call solofd(mn,radmeo,angmat,nr,id,ida,idb,offdg,iang,&
           iang,iang,mvib2,mvib1,coefa2,coefb1,&
           ibassa,ibassb)
      noffblk=noffblk+1
         write(iscr)ip2,jp1
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if
!  k12 vs k22  (cor)

      k1=ka2
      k2=kb2
      mvib1=mviba2
      mvib2=mvibb2
      mn = mvib1*mvib2
      if (mn.ne.0) then
      write(6,*)'calculating '
      write(6,1010)k1,iq1,is1,k2,iq2,is2,'  cor',mvib1,mvib2
      write(6,*)jp2,ip2
      ALLOCATE(offdg(mn))
      call angin_cor(angmat,k1,ivb2,iva2,iang,iang,angfac,iqb)
      call solofd(mn,radmeo,angmat,nr,id,idb,ida,offdg,iang,&
           iang,iang,mvib2,mvib1,coefb2,coefa2,&
           ibassb,ibassa)
      noffblk=noffblk+1
         write(iscr)jp2,ip2
      call outrow(offdg,mn,iscr)
      deallocate(offdg)
      end if

! move k12 -> k1 and k22 -> k2
      ip1=ip2
         ka1=ka2
         plega1=plega2
         iva1=iva2
         mviba1=mviba2
         deallocate(coefa1)
         allocate(coefa1(ibassa,mviba1))
         coefa1=coefa2
         deallocate(coefa2)
      jp1=jp2
         kb1=kb2
         plegb1=plegb2
         ivb1=ivb2
         mvibb1=mvibb2
         deallocate(coefb1)
         allocate(coefb1(ibassb,mvibb1))
         coefb1=coefb2
         deallocate(coefb2)
         end do

      deallocate(coefa1,coefb1)
      deallocate(plega1,plegb1,plega2,plegb2,iva1,ivb1,iva2,ivb2)

901   write(6,*)'Number off-diag blocks calculated ...',noffblk
      return
    end subroutine solrt2 
!#####################################################################
!#####################################################################
      subroutine rdcoef(coef,idim,mev,iv)
!     read first step vector array coef from unit iv
      implicit double precision (a-h,o-y), logical (z)
      dimension coef(idim,mev),ro(idim)
 
      do i=1,mev
         call getrow(ro,idim,iv)
!         ss=ro(1)
!         do ii=1,idim
!            if (abs(ro(ii)).gt.abs(ss)) ss=ro(ii)
!         end do
!         xs=ss/abs(ss)
         xs=1.d0
         ro=xs*ro
         do j=1,idim
            coef(j,i)=ro(j)
         end do
      end do
      return
    end subroutine rdcoef
!########################################################################
      subroutine solofd(mn,radmat,angmat,nr,iq,iq1,iq2,offdg,iang11,&
                 iang1,iang2,mvib1,mvib2,coef1,coef2,ibass1,ibass2)
 
!     construct off-diagonal matrix elements from radial and angular
!     matrix elements and first step vectors.
 
      implicit double precision (a-h,o-y), logical (z)

      real*8, DIMENSION(mn) :: offdg
      real*8, DIMENSION(*) :: coef1
      real*8, DIMENSION(*) :: coef2
      real*8, DIMENSION(iang1,iang2) :: angmat
      real*8, DIMENSION(*) :: radmat
 
      data x0/0.0d0/
!     zero the off-diagonal block
      offdg = x0

      ir1=0
      ir2=0
      ir=0
      do  i1=1,nr
      do  i2=1,i1-iq
      i0=ir1*iang1
      ir1=ir1+1
      j0=ir2*iang2
      ir2=ir2+1
      ir=ir+1
      rm=radmat(ir)

      do  ia1=1,iang1
      i=i0+ia1
      do  ia2=1,iang2
      j=j0+ia2
!     cor is the integral for current basis function
      cor=rm*angmat(ia1,ia2)

!     transform to vibrational basis and add to offdg
      call dger(mvib2,mvib1,cor,coef2(j),ibass2,coef1(i),ibass1,&
            offdg,mvib2)
   end do !ia2
end do !ia1
end do !i2
      ir1=ir1+iq1
      ir2=ir2+iq2
   end do !i1
!     dump completed block to unit iscr

      return
    end subroutine solofd

!##########################################################################
      subroutine dgrot(k1,ezero,ksize,noffblk,mvib)
 
!     subroutine diag solves the eigenvalue problem:                #012
!          hamil * vec = eval * vec
!     by using some lapack routines if they work.
 
      implicit double precision (a-h,o-y), logical (z)
      character JOB*1
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
      common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
      real*8, DIMENSION(ksize) :: EVAL
      real*8, DIMENSION(ksize) :: EVALCM
      real*8, DIMENSION(8*ksize) :: work
      integer, DIMENSION(5*ksize) :: iwork
      integer, DIMENSION(ksize) :: ifail
      real*8, ALLOCATABLE, DIMENSION(:) :: vec
      real*8, ALLOCATABLE, DIMENSION(:,:) :: evec
      integer, DIMENSION(nktot) :: mvib
      integer, dimension(ksize) :: ibig

!          autocm converts atomic units (hartree) to cm-1.
      data autocm/2.19474316d+05/
      data x0/0.0d0/

      ifail=0          
      info=0
      abstol=x0
!     load the hamiltonian matrix elements
      ibass=nbass
      kss=ksize*(ksize+1)/2
      allocate(vec(kss))  

      call loadh(mvib,vec,noffblk,ksize,kss)

      if (zvec) then
         JOB='V'
         allocate(evec(ksize,neval))
         write(6,*)'Job type: evalues & evectors'
         else
         JOB='N'
         allocate(evec(1,1))
         write(6,*)'Job type: evalues only'
      end if

      CALL DSPEVX(JOB,'I','L',ksize,vec,x0,x0,1,neval,&
       ABSTOL,nevals,eval,evec,ksize, WORK, IWORK, IFAIL,INFO )

!      CALL DSPEV(JOB,'L',ksize,vec,eval,evec,ksize, WORK, INFO )
!      nevals=neval

      deallocate(vec)

      IF (INFO .NE. 0) then
         WRITE(6,900) INFO
900  FORMAT(//5X,'DSPEVX RETURNED INFO =',I3)
         WRITE(96,901) IFAIL
901  FORMAT(//5X,'DSPEVX RETURNED IFAIL =',I3)
         stop
      END IF
      write(6,*)'Evalues converged .... ',nevals,' / ',neval

      nevalp=min(nevals,neval)

!     print eigenvalues in atomic units & wavenumbers

      write(6,1000) nevalp
 1000 format(//5x,'lowest',i4,' eigenvalues in hartrees',/)
      write(6,1020) (eval(i),i=1,nevalp)
 1020 format(5d24.12/)
      if (zpfun) then
         if (k1 .eq. 0) then
            open(unit=ilev,form='formatted')
            rewind ilev
  200       read(ilev,*,end=210,err=210)
            goto 200
  210       continue
! ***** inclusion of the following card is machine dependent *****
            backspace ilev
         endif
         ip=mod(jrot+kpar,2)
         write(ilev,1025) jrot,ip,0,0,(2-4*iqpar),nevalp
 1025    format(6i4)
         write(ilev,1026) eval
 1026    format(4d20.12)
      endif

      if (zvec) then
!        write eigenvalues, eigenvectors, etc to stream jvec
         open(unit=jvec,form='unformatted')
         rewind jvec
         write(jvec) jrot,kpar,iqpar,nevalp,ibass
         write(jvec) mvib
         call outrow(eval,nevalp,jvec)

         mbeg=1
         do  k=1,nktot
            mvv=mvib(k)
            if (mvv.ne.0) then
            mend=mbeg+mvv-1
            write(jvec) ((evec(j,i),j=mbeg,mend),i=1,nevalp)
            mbeg=mbeg+mvv
            end if
         enddo
      endif

      do 60 i=1,nevalp
      evalcm(i) = eval(i) * autocm - ezero
   60 continue
      write(6,1010) nevalp,ezero
 1010 format(//5x,'lowest',i4,' eigenvalues in wavenumbers relative to',&
                  ' ezero =',d20.12/)
      write(6,1020) (evalcm(i),i=1,nevalp)

      write(6,*)'For simplicity of use evalues printed also in stream ',62
      do i=1,nevalp
         write(62,1034)i,eval(i),evalcm(i)
      end do
1034 format(1I4,2f30.20)

      if (zpvec) then
          if (thresh .le. x0) then
!             print complete eigenvectors
              write(6,1030)
 1030         format(//'    eigenvectors',/)
              do 70 i=1,neval
              write(6,1040) (evec(j,i),j=1,ibass)
 1040         format(/(1x,10f13.7))
   70         continue
          else
!             print largest components of the eigenvectors
              write(6,1050) thresh
 1050         format(//'eigenvector components greater than thresh =',&
                         f5.2)
              do 80 i=1,neval
              vbig=x0
              ipt=0
              do 90 j=1,ibass
              vv=abs(evec(j,i))
              if (vv .gt. thresh) then
                  ipt=ipt+1
                  ibig(ipt)=j
              endif
              if (ipt .le. 0 .and. vv .gt. vbig) then
                  vbig=vv
                  ibig(1)=j
              endif
   90         continue
              write(6,1060) i,(ibig(j),evec(ibig(j),i),j=1,max(1,ipt))
 1060         format('  vector',i3,5(i7,f14.10)/(11x,5(i7,f14.10)))
   80         continue
          endif
      endif

      deallocate(evec)

!     transform the eigenvectors if requested
      if (ztran) then         
        call dstorepb2(mvib,eval)
      endif

      write(6,*)'Successful end in routine dgrot'

      return !dgrot
      end

!######################################################################
      subroutine getrow(row,nrow,iunit)
!     fast non-formatted read                                       #015
      real*8, DIMENSION(NROW) :: ROW
      read(iunit) row
      return
    end subroutine getrow

!######################################################################
      subroutine outrow(row,nrow,iunit)
!     fast non-formatted write                                      #016
      real*8, DIMENSION(NROW) :: ROW     
      write(iunit) row
      return
    end subroutine outrow

!#######################################################################

      MODULE constants
      IMPLICIT NONE
      integer, PARAMETER :: real_kind=SELECTED_REAL_KIND(8,40)
      END MODULE constants
!#######################################################################

      SUBROUTINE jac_basis(nn,nb,alf,bet,x,basis)
      USE constants
      IMPLICIT NONE
      integer :: n,i,nb,nn
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
      integer :: n,I,nn,nb
      REAL(KIND=real_kind),EXTERNAL :: gammln
      REAL(KIND=real_kind) :: alf, bet,lmd, x0,x1,x2
      REAL(KIND=real_kind) :: x(nn),bass(0:nb,nn)
      REAL(KIND=real_kind), ALLOCATABLE :: A1n(:),A2n(:),A3n(:),&
      &                                    A4n(:)
      DATA x0,x1,x2/0.0d0,1.0d0,2.0d0/
      lmd=alf+bet+x1
      ALLOCATE(A1n(nb),A2n(nb),A3n(nb),A4n(nb))
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
      DEALLOCATE(A1n,A2n,A3n,A4n)
      END SUBROUTINE jacgt

!##################################################################

      SUBROUTINE gaujac(x,w,n,alf,bet)
      USE constants
      IMPLICIT NONE
      integer :: n
      REAL(KIND=real_kind):: alf,bet,x1,x2,x3
      REAL(KIND=real_kind):: w(n),x(n)
      REAL(KIND=real_kind),PARAMETER :: EPS=3.0D-14
      integer,PARAMETER :: MAXIT=10
      integer :: i,its,j
      REAL(KIND=real_kind)::alfbet,an,bn,r1,r2,r3
      REAL(KIND=real_kind)::c1,c2,c3,p1,p2,p3,pp,temp,z,z1
      REAL(KIND=real_kind),EXTERNAL :: gammln
      DATA x1,x2,x3/1.0d0,2.0d0,3.0d0/

!      write(6,*)'entering gaujac.........'
!      write(6,*)alf,alf**2
!      write(6,*)bet,bet**2

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

!####################################################################

       FUNCTION gammln(XX)
       USE constants
       IMPLICIT NONE
       integer :: j
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
11     CONTINUE
       GAMMLN=TMP+LOG(STP*SER)
       RETURN
       END FUNCTION gammln

!####################################################################

       SUBROUTINE norms2(norm,nn,alf,bet)
       USE constants
       IMPLICIT NONE
       integer :: n,nn
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

!###################################################################
!#######################################################################
!#######################################################################
      subroutine dstorepb2(mvib,energy)
 
!     dstore transforms the eigenvectors of the second variational
!     step into ones for the first step basis and stores the
!     results in a form suitable for program dipole3.
!     This version treats each k-block seperately.
!  PB: this is a rude and simple version, that assumes the same mvib
!   for all the k-like blocks, that is, one needs ibass=0 in the input
!   If for any reason this condition is not satisfied, than it tells
!  it and exits.
!  A more elaborate version will be programmed when and if necessary.


      implicit double precision(a-h,o-y), logical(z)
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
                common /pb/ inda1(100),inda2(100),indb1(100),indb2(100),indk(100),&
          iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,&
           nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
      integer, DIMENSION(nktot) :: MVIB
      real*8, DIMENSION(neval) :: ENERGY
      real*8, DIMENSION(nr) :: r
      real*8, ALLOCATABLE, DIMENSION(:) :: w_leg,x_leg,wr,r2,jxcos0,jwalf0
      real*8, ALLOCATABLE, DIMENSION(:) :: wmax
      real*8, ALLOCATABLE, DIMENSION(:,:) :: wfs,wf
      real*8, DIMENSION(3) :: xmass,xmass_1
      integer free

      parameter (iplot=39)

!      write(6,*)'This is a rude version for dstore: pb, 13/07/05.'
!      write(6,*)'Only works if mvib is the same for all k-blocks'
      write(6,*)'This is the final version for dstore: pb, 29/06/06.'
      write(6,*)'hovewer, needs extensive checking....'


      free=0

! set the number of entries in the input files ipa,ipb
! for each k-block

      write(6,1000) ipa,ipb,jvec,kvecpb,kvec
 1000 format(/5x,'eigenvector transformation:',&
             /5x,'input:   dvr3d data (stream A), ivec =',i3,&
             /5x,'input:   dvr3d data (stream B), ivec1 =',i3,&
             /5x,'         rotlev3b data (internal stream), jvec =',i3,&
             /5x,'output:  transformed vectors (DVR), kvecpb =',i3,&
             /5x,'output:  transformed vectors (FBR), kvec =',i3/)
!     read dvr3d header
      rewind ipa
      rewind ipb
      read(ipb) idia_1,ipar_1,idvr_1,npnt1_1,npnt2_1,jrot0_1,kmin0_1,mval_1
      read(ipa) idia1,ipar1,idvr1,npnt1,npnt2,jrot0,kmin0,mval
      read(ipb) zembed_1,zmors1_1,zmors2_1,xmass_1,g1_1,g2_1,zncor_1,zquad2_1
      read(ipa) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
      read(ipb) re1_1,diss1_1,we1_1,re2_1,diss2_1,we2_1
      read(ipa) re1,diss1,we1,re2,diss2,we2
      read(ipb)
      read(ipa)
      read(ipb)
      read(ipa) r
!     read rotlev header
      rewind jvec
      read(jvec) jrot1,kmin1,ipar1,nval,ibass
      read(jvec) mvib

      nval=min(nval,neval)
      IF (ZPTRA) WRITE(6,1010) JROT0_1,KMIN0_1,IPAR_1,NVAL
 1010 FORMAT(/5X,'J =',I3,' KMIN =',I2,'  PARITY =',I2,' NVAL =',I3/)
!     size of final basis?
      nend=0

      nkbas_even=ndvr*(nr+1)*nr/2
      nkbas_odd=ndvr*(nr-1)*nr/2
      nkbas_tot=nr*nr*ndvr
         if (isa.eq.0) nkbas_a=nkbas_even
         if (isa.eq.1) nkbas_a=nkbas_odd
         if (isb.eq.0) nkbas_b=nkbas_even
         if (isb.eq.1) nkbas_b=nkbas_odd
         if ((nkbas_b+nkbas_a).ne.nkbas_tot) then
            write(6,*)'problem in assigning 3d dimension...'
            stop
         end if
         if (kpar .eq. 1) then
            nend=nkbas_a*nka+nkbas_b*nkb
            ipar_new=0
            write(6,*)'Running over odd k values'
         else
            write(6,*)'Running over even k values'
            ipar_new=iqpar
            nend=nkbas_a*nka+nkbas_b*nkb
         endif


         if (zplot) then

            allocate(jxcos0(npnt),jwalf0(npnt))
            alf = 0.d0

            call gaujac(jxcos0,jwalf0,npnt,alf,alf)

      write(36,1003) npnt,alf,(jxcos0(i),jwalf0(i),i=1,npnt/2)

1003  format(//i8,' point gauss-associated legendre integration with ',&
             'ALF =' ,f10.5//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))

            allocate(r2(nr),wr(nr))
      call radgrid(zmors1,re1,diss1,we1,nr,wr,r2,idia,xmass)

      if (abs(r2(1)-r(1)).gt.1.d-13) then
         write(6,*)'Detected inconsinstency in the radial grid'
         write(6,*)'Retrieved from file ',ipa,' : ' 
         write(6,*)r(1)
         write(6,*)'calculated locally: '
         write(6,*)r2(1)
         stop
      end if

      isize=nr*(nr+1)*npnt/2
      if (nplotf.eq.0) nplotf=nval
      nplotf=min(nplotf,nval)
      if (nploti.gt.nplotf) then
         write(6,*)'Inconsinstency in call to wf plotting program'
         write(6,*)'Requested to print from wf: '
         write(6,*)'[ ',nploti,' - ',nplotf,' ] '
         write(6,*)'(number of rotational levels found ',nval,' )'
         stop
      else
         write(6,*)'Requested to print from wf: '
         write(6,*)'[ ',nploti,' - ',nplotf,' ] '
         write(6,*)'(number of rotational levels found ',nval,' )'
      nplot=nplotf-nploti+1
      end if
      allocate(wfs(nplot,isize),wmax(nplot))
      wfs=0.d0
      wmax=0.d0

      open(unit=iplot,file='wfs/waves.info')

      write(iplot,*)jrot,idia,kpar,iqpar,thresh
      write(iplot,*)zbisc,zperp,zembed,zmors1,zmors2,xmass,g1,g2
      write(iplot,*)re1,diss1,we1,re2,diss2,we2
      write(iplot,*)nr
      write(iplot,*)r2
      write(iplot,*)wr**2
      write(iplot,*)npnt
      write(iplot,*)jxcos0
      write(iplot,*)jwalf0
      write(iplot,*)nploti,nplot
      write(iplot,*)(energy(i),i=nploti,nploti+nplot-1)
      wr=1.d0/wr
   else

      write(6,*)'wfs to be plotted not requested.....'
      isize=1
      nplot=1
      allocate(wfs(1,1))
      wfs=0.d0
   end if


      !1e 1o 3e 3o .... for odd K
      !0e 2o 2e ... or 0o 2e 2o ... for even K

!     write header on new file
      open(unit=kvec,form='unformatted')
      open(unit=kvecpb,form='unformatted')

      write(kvec) idia,ipar_new,ndvr,nr,nr,jrot,&
           free,nval,kpar,iqpar
      write(kvec) zembed,zmors1,zmors2,xmass,g1,g2,zncor
      write(kvec) re1,diss1,we1,re2,diss2,we2
      write(kvec) nend,free,free,free
!     transfer the dvr points
      write(kvec) r
      write(kvec) nval

      write(kvecpb) idia,ipar_new,ndvr,nr,nr,jrot,&
           free,nval,kpar,iqpar
      write(kvecpb) zembed,zmors1,zmors2,xmass,g1,g2,zncor
      write(kvecpb) re1,diss1,we1,re2,diss2,we2
      write(kvecpb) nend,free,free,free
!     transfer the dvr points
      write(kvecpb) r
      write(kvecpb) nval

      call getrow(energy,nval,jvec)
      call outrow(energy,nval,kvec)
      call outrow(energy,nval,kvecpb)
 
!--------- start transformation step -----------------

      !start loop over all k-blocks
      ! remember that in file ipa they are ordered:
      ! 0  1(q=0)  1(q=1)  2  3  4  5 .....

      write(6,*)'allocating legendre points and weights, size= ',npnt
      allocate(x_leg(npnt),w_leg(npnt))

      xjp1=dble(jrot*(jrot+1))
         i0=1
         j0=1
         iadd=0
! do extra ka block if necessary.....
      if (kmina.ne.kminb) then
! there is no point to skip as it should be the first block, however
! lets check this
         i=inda1(1)
         k=indk(i)
         mvib0=mvib(inda2(1))
         if (k.ne.0) then
            write(6,*)'Problem in dstorepb2 '
            write(6,*)' i =  ',i
            write(6,*)' k =  ',k
            write(6,*)'but k should be 0 !!!! I stop here.'
            stop
            end if
         allocate(wf(nplot,nkbas_a))
         wf=0.d0
!            if (zplot.and.mvib0.gt.0) then
!               allocate(wf(nplot,nkbas_a))
!               wf=0.d0
!            else
!               allocate(wf(1,1))
!               wf=0.d0
!            end if
         alf=sqrt((xjp1-dble(k**2))/2.d0)
         call gaujac(x_leg,w_leg,npnt,alf,alf)
         call transblock(ipa,kvec,kvecpb,mvib0,nkbas_a,neval,jvec,isa,zplot,nplot,nploti,wf)
         write(kvec)x_leg
         write(kvecpb)x_leg
         i0=2
         iadd=1
         if (zplot.and.mvib0.gt.0) then
            call sumwf(wr,wf,wfs,alf,jxcos0,isa,nr,iang,npnt,nkbas_a,isize,nplot,jwalf0)
         end if
            deallocate(wf)
         end if

      do 50 jj=1,nkb

         i=inda1(jj+iadd)
! iadd accounts for the extra first block if necessary
         ka=indk(i)
         mvib0=mvib(inda2(jj+iadd))

         alf=sqrt((xjp1-dble(ka**2))/2.d0)
         call gaujac(x_leg,w_leg,npnt,alf,alf)

         allocate(wf(nplot,nkbas_a))
         wf=0.d0
!            if (zplot.and.mvib0.gt.0) then
!               allocate(wf(nplot,nkbas_a))
!               wf=0.d0
!            else
!               allocate(wf(1,1))
!               wf=0.d0
!            end if

         do l=i0,i-1
            call skipblock(nskipka,ipa)
         end do
         call transblock(ipa,kvec,kvecpb,mvib0,nkbas_a,neval,jvec,isa,zplot,nplot,nploti,wf)
         write(kvec)x_leg
         write(kvecpb)x_leg
         if (zplot.and.mvib0.gt.0) then
      call sumwf(wr,wf,wfs,alf,jxcos0,isa,nr,iang,npnt,nkbas_a,isize,nplot,jwalf0)
         end if
            deallocate(wf)

         j=indb1(jj)
         kb=indk(j)
         mvib0=mvib(indb2(jj))

         allocate(wf(nplot,nkbas_b))
         wf=0.d0
!            if (zplot.and.mvib0.gt.0) then
!               allocate(wf(nplot,nkbas_b))
!               wf=0.d0
!            else
!               allocate(wf(1,1))
!               wf=0.d0
!            end if

         do l=j0,j-1
            call skipblock(nskipkb,ipb)
         end do
         call transblock(ipb,kvec,kvecpb,mvib0,nkbas_b,neval,jvec,isb,zplot,nplot,nploti,wf)
         write(kvec)x_leg
         write(kvecpb)x_leg
         if (zplot.and.mvib0.gt.0) then
            call sumwf(wr,wf,wfs,alf,jxcos0,isb,nr,iang,npnt,nkbas_b,isize,nplot,jwalf0)
         end if
            deallocate(wf)
         i0=i+1
         j0=j+1
   50 continue

         deallocate(x_leg,w_leg)

         if (zplot) then 
            call plotwfs(wfs,nploti,nplot,ithre,isize,wmax)
            write(iplot,*)wmax
            close(unit=iplot)
         end if
 
      write(6,1050) jrot1,nval,nend
 1050 format(//5x,'transformation completed successfully',&
             //5x,'for rotational state: jrot =',i3,&
             '.  stored',i6,' vib-rot levels',&
             /5x,'nend =',i7)
      rewind jvec
      rewind kvec
      return
      end

!#######################################################################
      subroutine skipblock(nskip,ivpb)

      implicit double precision(a-h,o-y), logical(z)

      do i=1,nskip
         read(ivpb)
         end do
      return
    end subroutine skipblock
!#######################################################################
      subroutine transblock(ivpb,k1,k2,mvib,nkbas,nval,jv,is,z,nwf,ni,wf)

      implicit double precision(a-h,o-y), logical(z)

      real*8, ALLOCATABLE, DIMENSION(:,:) :: pleg
      real*8, ALLOCATABLE, DIMENSION(:,:) :: c
      real*8, ALLOCATABLE, DIMENSION(:,:) :: fbrvec
      real*8, ALLOCATABLE, DIMENSION(:,:) :: dvrvec
      real*8, ALLOCATABLE, DIMENSION(:,:) :: d
      real*8, dimension(nwf,nkbas) :: wf
      integer, ALLOCATABLE, DIMENSION(:) :: iv1

      x1=1.d0
      x0=0.d0
   allocate(d(nkbas,nval))

      if (mvib.gt.0) then
        if (z) wf=x0

   allocate(c(mvib,nval),fbrvec(nkbas,mvib),dvrvec(nkbas,mvib))
   read(jv) ((c(j,i),j=1,mvib),i=1,nval)

      read(ivpb)
!      write(k1)x
!      write(k2)x

      read(ivpb) kz,maxleg,idvr,lincr
      write(k1) kz,is,maxleg,idvr,lincr
      write(k2) kz,is,maxleg,idvr,lincr

      allocate(pleg(idvr,idvr))
      
      call getrow(pleg,idvr*idvr,ivpb)
      read(ivpb) iang1,ibass1

      allocate(iv1(idvr))

      read(ivpb) (iv1(ii),ii=1,idvr)
      read(ivpb) meval
      read(ivpb)

      call jtran2(fbrvec,dvrvec,mvib,pleg,is,iv1,nkbas,ivpb)

!skipping remaining evectors

      nskip=meval-mvib
      call skipblock(nskip,ivpb)
      
      !the code below here performs a standard matrix multiply of the form
      !D = B.C
      !where the matrices are dimensioned d(nkbas(k),nval), 
      !b(nkbas(k),mvib(k)) and c(mvib(k),nval). 
! this is for the FBR  eigenvectors      
      d=0.d0
      call dgemm('N','N',nkbas,nval,mvib,x1,fbrvec,nkbas,c,mvib,&
           x0,d,nkbas)
      do i=1,nval
      write(k1)(d(j,i),j=1,nkbas)
      end do
      if (z) then
         do i=1,nwf
            do j=1,nkbas
               wf(i,j)=d(j,i+ni-1)
            end do
         end do
      end if
      d=0.d0
! this is for the DVR  eigenvectors      
      call dgemm('N','N',nkbas,nval,mvib,x1,dvrvec,nkbas,c,mvib,&
           x0,d,nkbas)
      do i=1,nval
      write(k2)(d(j,i),j=1,nkbas)
      end do
      deallocate(c,dvrvec,fbrvec)
      
      else

         wf=x0

      read(ivpb)
      read(ivpb) kz,maxleg,idvr,lincr
      write(k1) kz,is,maxleg,idvr,lincr
      write(k2) kz,is,maxleg,idvr,lincr

      allocate(pleg(idvr,idvr))
      
      call getrow(pleg,idvr*idvr,ivpb)
      read(ivpb) iang1,ibass1

      allocate(iv1(idvr))

      read(ivpb) (iv1(ii),ii=1,idvr)
      read(ivpb) meval
      read(ivpb)
      nskip=meval
      call skipblock(nskip,ivpb)

         write(6,*)'mvib = ',mvib,' for block: '
         write(6,*)'   Kz   = ',kz
         write(6,*)'   is   = ',is
         write(6,*)' Stream = ',ivpb

         d=x0
      do i=1,nval
      write(k1)(d(j,i),j=1,nkbas)
      end do
      do i=1,nval
      write(k2)(d(j,i),j=1,nkbas)
      end do

      end if
      deallocate(d,pleg,iv1)

      return
    end subroutine transblock 
!#######################################################################
!#######################################################################
      subroutine jtran2(coeffbr,coeffdvr,mvib,pleg,is,iv,nkbas,ifort)
  
      implicit double precision (a-h,o-y), logical (z)
      common /size/ nbass,neval,idia,nr,maxblk_even,jrot,&
                    ndvr,iang,npnt,maxblk_odd,ibass,&
                    nktot,kpar,iqpar
      common /outp/ thresh,zpham,zpvec,zvec,ztran,zptra,zcut,idiag,&
                    zpfun,zplot,ilev,ivec,ivec1,jvec,kvecpb,kvec,iscr &
                    ,nploti,nplotf,ithre 
          
      real*8, DIMENSION(ndvr,iang) :: pleg
      real*8, DIMENSION(nkbas,mvib) :: coeffdvr
      real*8, DIMENSION(nkbas,mvib) :: coeffbr
      real*8, allocatable :: sumt(:,:)
      DIMENSION iv(iang)
      data x0/0.0d0/

      nrad=nr*(nr+1-2*is)/2

      allocate(sumt(iang,nrad))
 
!     transform back to the original fbr-type basis in the
!     associated legendre functions
      do 10 l=1,mvib
!     first read in a new vector
      call getrow(sumt,nkbas,ifort)

      i=0
      do i2=1,nrad
      do i1=1,iang
         i=i+1
         coeffdvr(i,l)=sumt(i1,i2)
      end do
      end do
      ipt=0
      do 50 mn=1,nrad
      do 20 j=1,iang
      sumk=x0
      do 40 k=1,iang
      if (iv(k) .eq. 0) goto 40
      sumk=sumk + sumt(k,mn) * pleg(j,k)
   40 continue
      coeffbr(ipt+j,l) = sumk
   20 continue
      ipt=ipt+iang
   50 continue
   10 continue

      deallocate(sumt)

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
subroutine sumwf(wr,wf,wfs,alf,jcos,s,nr,nt,nt2,n3dl,n3d,nwf,jw)
implicit none
integer nt,nr,l,n1,n2dl,n3dl,i,k2,i0,k1,nst,s,n3,n2,n,nwf,n3d,nt2
real*8 jcos(nt2),pleg(0:nt-1,nt2),jx(nt2),wr(nr),jw(nt2)
real*8 sumw,alf,alf2,sum1,sum2,sum3
real*8 wf(nwf,n3dl),wfs(nwf,n3d),wft(n3d)

CALL jac_basis(nt2,nt-1,alf,alf,jcos,pleg)

do l=1,nt2
jx(l)=dsqrt(jw(l)*(1.d0-jcos(l)**2)**alf)
end do

do n1=0,nt-1
do l=1,nt2
pleg(n1,l)=pleg(n1,l)*jx(l)
end do
end do

!write(71,*)nt2,nt
!do n1=1,nt2
!do n2=1,nt2
!sum3=0.d0
!do l=0,nt-1
!sum3=pleg(l,n1)*pleg(l,n2)+sum3
!end do
!write(71,*)n1,n2,sum3
!end do
!end do

n2dl=nr*(nr+(-1)**s)/2

do nst=1,nwf

wft=0.d0

i0=0
do n1=1,nr
do n2=1,n1-s
   n3=n1*(n1-1)/2+n2-1
   do k2=1,nt2
      n=n3*nt2+k2

   sumw=0.0
do k1=1,nt
   sumw=sumw+wf(nst,i0+k1)*pleg(k1-1,k2)
end do !k1

      wft(n)=sumw!*wr(n1)*wr(n2)
   end do !k2
      i0=i0+nt
end do !n2
end do !n1

do n=1,n3d
wfs(nst,n)=wft(n)**2+wfs(nst,n)
end do ! n
end do ! nst

return
end subroutine sumwf
!########################################################################

subroutine lagbasis(aa,nlag,dg,wlag,wln,zd,csx)
  implicit none
  !  zd(i,j)=N_{i-1} L^a_{i-1}(x_j) exp{-x_j/2} sqrt{w_j}       
  !  wlag=dsqrt{w}
  !  wln =dln{w}
  !  nlag= number of points
  !  aa=alfa
  !  csx=sum of points
  
  integer :: nlag
  integer :: info,j,i,di
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
     wln(i)=cu2+dlog(dabs(zd(1,i)))+0.5d0*dg(i)
     wlag(i)=dexp(wln(i))
  end do
  
  do j=1,nlag
     zc=zd(1,j)/wlag(j)
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
      TSA = X1
      csx=0.d0
      csa=0.d0
      do l=1,npnt
         csx=csx+y(l)
         csa=csa+dexp(-y(l)+2.d0*wln(l))
      end do
! set up laguerre's normalis 
      do l=0,NPNT-1
         xl=dble(l)
         dd=gammaln(xl+1.d0)-gammaln(xl+alf+1.d0)
         dnormnew(l)=dexp(dd*0.5d0)
      end do
      csa=csa*(dnormnew(0)**2)

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
            wt(l)=wt(l)/sqrt(beta*y(l)**ALF2)
         end do
      else
         do l=1,npnt
            wt(l)=wt(l)/sqrt(sqrt(2.d0)*(beta**0.25d0)*(y(l)**ALF2))
         end do
      end if
!     do n1=0,nmax
!        do n2=0,nmax
!           sum=0.d0
!           do l=1,npnt
!              sum=sum+f(n1,l)*f(n2,l)*wt2(l)*wt2(l)
!           end do
!!           write(71,*)n1,n2,sum
!        end do
!     end do

!        write(71,*)ALF,re
     do l=1,npnt
!        write(71,*)l,y(l),r(l)
!        write(71,*)wt(l),wt2(l)
!        wt(l)=wt(l)/(sqrt(2.d0)*beta*y(l)**(ALF+1.d0))
        
      if (myid .eq. 0) WRITE(6,1010) Y(l),WT(l),R(l)
 1010 FORMAT (F23.15,D25.12,F13.5)
      IF ((R(I) .LT. X0) .and. (myid .eq. 0)) WRITE(6,1015) I
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
!#######################################################################
!#######################################################################
!              THE    END      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#######################################################################
!#######################################################################
