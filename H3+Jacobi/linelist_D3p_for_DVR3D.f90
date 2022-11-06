program linelist_DVR3D

 implicit none

 integer j,jf,ji,ilevel,ilevelf,ileveli,igamma,igammaf,igammai,info,j0,i,k,jgamma,kminf,kmini,iparf,ipari,idx
 integer nfiles, jmax, ilines,nenfiles,ifile
 integer nlevels,v_(9),numj(0:200,1:4),itermi,itermf,symv,sym_t,J_t,v_norm(9),igammar,imis,itermi_,itermf_

 real(8),allocatable :: energies(:,:,:),ee1(:),ee2(:),ss(:)
 integer,allocatable :: sym(:),Jktau(:,:),NN(:,:,:),Nlines(:)!,ef(:)

 logical :: energeyfile_do = .false.,record_segmented = .false.

 real(8) :: acoef,abscoef,energy,energyf,energyi,tranfreq,linestr,ZPE,coeff,ztran,xtran,dipole,gg,enermax=1e9,enercutoff=1e9

 character(4) gamma,gammaf,gammai,gammav,jchar,gchar,gammar
 character(50) :: intfilename(400),enrfilename,enrfilenameout,intfilenameout,dscrfilename,filename,enfilename(400)
 character(50) :: record_type = " "
 character(1)  :: ch_t1
 character(len=8)  :: ch_t8
 character(200) :: ch_t200
 character(400) :: form
! character(1), parameter :: efpar(2)=(/'e','f'/)
 character(1), allocatable :: ef(:,:,:),qpar(:,:,:)
 integer        :: unit_f,ilevelmax
 character(5)   :: symm(4) = (/'A1','A2','B1','B2'/)
 integer(4)     :: gns(4) 
 real(8),parameter :: e_thresh = 1d-4,tocm = 2.19474624d+05,todebye=2.5417662d0
 integer(4),parameter :: itra = 13
 integer(4)     :: jrot,kpar,ipar,irottau,meval(2),kronig(2),nsym(0:1,0:1),nasym(0:1,0:1),iasym,isym,Lmax,is,rvsym(0:1,0:1,-2:2)
 integer(4)     :: j1,j2,kmin1,kmin2,neval1,neval2,idia,ibase1,ibase2,ie1,ie2,ipar1,ipar2,p1,p2,j_
 integer(4)     :: ierr
 real(8)        :: gz,energy_tmp
 logical        :: zembed
  
 type numT
      !
      integer,pointer   :: nn(:)
      integer,pointer   :: sym(:)
      integer,pointer   :: J(:)
      real(8),pointer   :: energy(:)
      !
 end type numT

 type(numT),allocatable :: level(:,:)

 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16,detosc=3.136186d-07


      interface
        subroutine get_energies_from_fort14(iunit,jrot,kpar,ipar,meval,lmax,jmax,ierr,enercutoff,enercalc,ef,qpar)
            integer(4),intent(in) :: iunit
            integer,intent(out)   :: jrot,ipar,kpar
            integer,intent(out)  :: meval(2),ierr
            integer(4),intent(in)  :: lmax,jmax
            real*8,intent(in),optional   :: enercutoff
            character(1),intent(inout), optional :: ef(lmax,0:jmax,4),qpar(lmax,0:jmax,4)
            real*8,intent(inout),optional   :: enercalc(lmax,0:jmax,4) 
         end subroutine get_energies_from_fort14
      end interface



    !read number of intensity files
    read*,nfiles
    !
    if (nfiles>400) then 
      print('("Too many files (>400):",i)'),nfiles
      stop 'Too many files'
    endif
    !
    !read intensities filenames
    do i = 1,nfiles
     read*,intfilename(i)
    enddo
    !
    !read number of intensity files
    read*,nenfiles
    !
    if (Nenfiles>400) then 
      print('("Too many en-files (>400):",i)'),Nenfiles
      stop 'Too many en-files'
    endif
    !
    !read intensities filenames
    do i = 1,Nenfiles
     read*,enfilename(i)
    enddo
    !
    !read descr-filename - in
    ! read*,dscrfilename
    !
    !read intensity filename - out
    read*,intfilenameout
    !
    !read energies filename - out 
    read*,enrfilenameout
    !
    !read enermax
    read*,enermax
    !
    !read enercutoff in Eh, absolute 
    read*,enercutoff
    !
    ! read the gns (4 records, one for symmetry)
    read*, gns(1), gns(2), gns(3), gns(4)
    !read record_type
    read(5,"(a50)",end=41),record_type
    !
    if (trim(record_type)=='segmented'.or.trim(record_type)=='SEGMENTED') then
      !
      record_segmented = .true.
      !
    endif 
    !
    41 continue
    !
    !read the energy file 
    !open(unit=11,file=trim(enrfilename))
    !
    i = 0
    numj = 0
    jmax = 0
    nlevels = 0
    !
    ! Define the molec. symmetry group correlatation with the internal DVR3D symmetry indeces.
    !
    ! symmetrical rot. basis functions 
    !
    nsym (0,0) = 1
    nsym (0,1) = 4
    !
    nasym(0,0) = 3
    nasym(0,1) = 2
    !
    ! asymmetrical rot. basis functions 
    !
    !
    nsym (1,0) = 4
    nasym(1,0) = 2
    !
    nsym (1,1) = 1
    nasym(1,1) = 3
    !
    !     ipar irot
    !
    rvsym(0,0,0) = 1
    rvsym(0,1,0) = 2
    rvsym(0,0,1) = 4
    rvsym(0,1,1) = 3
    !
    rvsym(1,0,0) = 3
    rvsym(1,1,0) = 4
    rvsym(1,0,1) = 2
    rvsym(1,1,1) = 1
    !
    Lmax = 0
    !
    print("('Counting energies...')")
    !
    ! go through all energy output files:
    do ifile = 1,Nenfiles
       !
       open(unit=1,file=trim(enfilename(ifile)),err=15)
       !
       print("( 'process ',a)"), trim(enfilename(ifile))
       !
       ilevelmax = 0
       !
       if (jmax>200) stop 'Jmax>200'
       !
      ! do 
         !
         call get_energies_from_fort14(1,jrot,kpar,ipar,meval,lmax,jmax,ierr)
         !print*,kronig, "KRONIG1"
         !!
         !!isym  = nsym(mod(jrot+2,2),ipar)
         !!iasym = nasym(mod(jrot+2,2),ipar)
         !
         isym = rvsym(mod(jrot+2,2),kpar,ipar)
         iasym = rvsym(mod(jrot+2,2),mod(kpar+3,2),mod(ipar+3,2))
         !
         if (numj(jrot,isym)>0) then 
           write(6,"('This J,ipar (1) was already processed: ',i4,i4,2x,a)") jrot,ipar,trim(enfilename(ifile))
           stop 'Error: duplicate J,ipar '
         endif
         !
         if (numj(jrot,iasym)>0) then 
           write(6,"('This J,ipar (2) was already processed: ',i4,i4,2x,a)") jrot,ipar,trim(enfilename(ifile))
           stop 'Error: duplicate J,ipar '
         endif
         !
         numj(jrot,isym)  = meval(1)
         numj(jrot,iasym) = meval(2)
         !
         nlevels = nlevels + sum(meval(:))
         !
         jmax = max(jmax,jrot)
         !
         Lmax = max(Lmax,maxval(meval(:)))
         !
         if (ierr/=0) exit
         !
         cycle
         !
       !enddo
       !
       close(1)
       !
       cycle
       !
     15 write(6,"('Error opennig energy file ',a)") enfilename(ifile)
       stop 'Error opennig energy file' 
       !
    enddo
    !
    print("(' jmax = ',i4)"),jmax
    !
    ! the totoal number of levels
    !
    print("(' nlevels = ',i)"),nlevels
    !
    allocate(sym(nlevels),energies(Lmax,0:jmax,4),Jktau(nlevels,3),NN(Lmax,0:jmax,4),ef(Lmax,0:jmax,4),qpar(Lmax,0:jmax,4),stat=info); 
    if (info/=0) stop 'error: sym,Jktau are out of memory'
    !
    sym = 0
    energies = 0
    NN = 0
    nlevels = 0 
    numj = 0
    !ef=0
    !
    print("('Reading energies...')")
    !
    ! go through all energy output files:
    !idx=1
    do ifile = 1,Nenfiles
      !
       open(unit=1,file=trim(enfilename(ifile)))
       !
       print("( 'process ',a,3x,i)"), enfilename(ifile)
       !
       !ilevelmax = 0
       !
       if (jmax>200) stop 'Jmax>200'
       !
       rewind(1)
       call get_energies_from_fort14(1,jrot,kpar,ipar,meval,Lmax,jmax,ierr,enercutoff,energies,ef,qpar)
       !print*,kronig,"KROONING 2"
       !
       isym = rvsym(mod(jrot+2,2),kpar,ipar)
       iasym = rvsym(mod(jrot+2,2),mod(kpar+3,2),mod(ipar+3,2))
       !
       if (numj(jrot,isym)>0) then 
         write(6,"('This J,ipar (1) was already processed: ',i4,i4,2x,a)") jrot,ipar,trim(enfilename(ifile))
         stop 'Error: duplicate J,ipar '
       endif
       !
       if (numj(jrot,iasym)>0) then 
         write(6,"('This J,ipar (2) was already processed: ',i4,i4,2x,a)") jrot,ipar,trim(enfilename(ifile))
         stop 'Error: duplicate J,ipar '
       endif
       !
       numj(jrot,isym)  = meval(1)
       numj(jrot,iasym) = meval(2)
       !ef(idx:idx+meval(1)-1)=efpar(kronig(1)+1)
       !print*,"first part of J:",idx,idx+meval(1)-1,ef(idx:idx+meval(1)-1)
       !if (jrot .ne. 0) ef(idx+meval(1):idx+meval(1)+meval(2)-1)=efpar(kronig(2)+1)
       !if (jrot .ne. 0) print*,"second part of J:",idx+meval(1),idx+sum(meval(:))-1,ef(idx+meval(1):idx+meval(1)+meval(2)-1)
       !
       nlevels = nlevels + sum(meval(:))
      ! if (jrot .ne. 0 ) then 
      !         idx=idx+sum(meval(:))
      !         print*,"la merda",idx
      ! else 
      !         idx=idx+meval(1)
      !         print*,"le feci"
      ! end if 
       !
       close(1)
       !
    enddo
    !
    ZPE = energies(1,0,1)*tocm
    !
    write(6,"('ZPE = ',f15.6)") ZPE
    !
    energies = energies*tocm-ZPE
    !
    open(unit=12,file=trim(enrfilenameout),action='write',status='replace')
    !
    ! if the file name is none we can skip creating the .states file, otherwise:
    !
    ilines = 0
    !
    ! go through all energy output files:
    !
    energeyfile_do = .true.
    !
    i = 0
    !
    do jrot = 0,jmax
      !print*,"jrot printed:",jrot
      !
      do isym = 1,4  
       ! print*,"isym:",isym,symm(isym)
        !
        do ilevel = 1,numj(jrot,isym)
          !print*,gns(isym)
          !
          if (gns(isym)<1) cycle 
          !
          i = i + 1
          !
          NN(ilevel,jrot,isym) = i
          !
          !prepare and store the .states file
!write(my_fmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,i7,1x,f13.6,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,2i8)"
!write(enunit,my_fmt)  iroot,energyI-intensity%ZPE,nint(intensity%gns(isymI)*( 2.0_rk*jI + 1.0_rk )),nint(jI),&
!                           lande,pm,ef,statename,ivI,(ilambdaI),nint((sigmaI)),nint((omegaI))

          if (energeyfile_do) then
            !
            ! Here is the main print-out of the .states file in the ExoMol format 
            !
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,a1,1x,a1,1x,a4,4(1xi2))") & 
               i,energies(ilevel,jrot,isym),gns(isym)*(2*jrot+1),jrot,qpar(ilevel,jrot,isym),ef(ilevel,jrot,isym),symm(isym),-2,-2,-2,-2
               !
          endif
        enddo
        !
        print("(i4,1x,i3,2x,2i8)"),jrot,isym,numj(jrot,isym),i
        !
      enddo
      !
    enddo
    !
    close(12)

! MP The file works up to here  - here things are fine as well, but I am not properly ure 
    !
    ! The second part: .trans file
    !
    print("(' N of levels =  ',i)"),nlevels
    !
    print("('Generate the Transition file ...')")
    !
    unit_f = 13
    !
    open(unit=unit_f,file=trim(intfilenameout),action='write',status='replace')
    !
    ! Some objects for the statistics 
    allocate(Nlines(0:jmax))
    ilines = 0
    Nlines = 0
    j_ = 0 
    !
    ! go through all intensity output files:
    do i = 1,nfiles
      !
      if (record_segmented) then
        open(unit=itra,file=trim(intfilename(i)),form='unformatted',recordtype='segmented',err=16) ! H2O
      else
        print*,"Opening an unsegmented file"
        open(unit=itra,file=trim(intfilename(i)),form='unformatted',err=16) ! H2S 
      endif
      !
      print("( 'process ',a,3x,i)"), intfilename(i),ilines
      !
      ! total number of lines in the given file
      !
      ilines = 0
      !
      ! number of not-matched states 
      imis = 0
      !
      ! find the j-output 
      do
         !
         read(itra,end=20,err=21) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2
         !write(42,*)  "j1:",j1,"j2:",j2,"kmin1:",kmin1,"kmin2:",kmin2,"neval1:",neval1,"neval2:",neval2,"idia:",idia,"ipar1:",&
         !              ipar1,"ipar2:",ipar2,"gz:", gz,"zembed:",zembed,"ibase1:",ibase1,"ibase2:",ibase2
         !
         !read(itra,end=20,err=21) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed
         !
         !!!!!!!!!!!!!!!
         !write(6,"(9i7,2f16.6)") j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2
         !!!!!!!!!!!!!!!
         !
         ji = j1
         jf = j2
         !
         p1 = abs(1-kmin1)
         p2 = abs(1-kmin2)
         !
         allocate(ee1(neval1), ee2(neval2), ss(neval1),stat=info)
         if (info/=0) stop 'error: ee1,ee2,ss are out of memory'
         !
         read(itra) ee1
         read(itra) ee2
         !
         !!!!!!!!!!!!!!!
         !write(6,"(9i7,2f16.6)") j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,ee1(neval1)*tocm-ZPE,ee2(neval2)*tocm-ZPE
         !!!!!!!!!!!!!!!
         !
         !
         igammai = igamma(J1,p1,ipar1)
         igammaf = igamma(J2,p2,ipar2)
         !
         if ( j_ /= max(j1,j2,j_)) then
            j_ = max(j1,j2,j_)
            write(6,"('J = ',i8,4x,i10)") j_,ilines
         endif 
         !
         energy_tmp = 0
         !
         do ie2=1,neval2
            !
            read(itra) ss 
            !
            energyf=ee2(ie2)*tocm-ZPE
            
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!! exclude !!!!!!!!!!!!!!!
            !if (j1==32.or.j2==32) cycle
            !if (j1==14.or.j2==14) cycle
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !write(6,"('f:',3i7,f16.6)") ie2,j2,igammaf,energyf
            !
            if (abs(energyf-energies(ie2,jf,igammaf))<=e_thresh) then 
              !
              itermf = NN(ie2,jf,igammaf) ;write(54,*)  
              
              !
            elseif ( ee2(ie2)<enercutoff ) then 
              !
              write(6,"('f:',3i7,2f16.6)") ie2,j2,igammaf,energyf,energies(ie2,j2,igammaf)
              !
              itermf = 0
              !
              do is = 1,4
                !
                !is = igammaf
                !
             
                do j = 1,numj(jf,is)
                  !
                  if (abs(energyf-energies(j,jf,is))<=e_thresh) then 
                    !
                    imis = imis + 1
                    itermf = NN(j,jf,is);write(54,*) 
            
                    !
                    exit
                    !
                  endif
                  !
                enddo
              end do 
                !
            endif
            !
            do ie1=1,neval1
               !
               energyi = ee1(ie1)*tocm-ZPE
           
               !
               if (min(energyi,energyf)>enermax) cycle
               !
               !if (abs(1698.165-energyf)<=e_thresh.and.abs(9329.673-energyi)<=e_thresh) then 
               !  continue 
               !endif
               !
               ! Find the i-state by using the energy threshold e_thresh (typically 1e-4)
               ! 
               itermi = 0
               !
               ! Check all states within the same J and symmetry to find the energy match and store the number to itermi
               !
               ilevel = ie1+numj(ji,igammai)-neval1
               !
               if (abs(energyi-energies(ie1,ji,igammai))<=e_thresh) then 
                 !
                 itermi = NN(ie1,ji,igammai)
                 ! 
               elseif (abs(energyi-energies(ilevel,ji,igammai))<=e_thresh) then 
                 !
                 itermi = NN(ilevel,ji,igammai)
                 !
               elseif ( ee1(ie1)<enercutoff ) then
                 !
                 write(6,"('i:',3i7,2f16.6)") ie1,ji,igammai,energyi,energies(ie1,ji,igammai)
                 !
                 !do is = 1,4
                   !
                   is = igammai
                   !
                   do j = 1,numj(ji,is)
                     !
                     if (abs(energyi-energies(j,ji,is))<=e_thresh) then 
                       !
                       imis = imis + 1
                       itermi = NN(j,ji,is)
                       !
                       !igammai = is
                       !
                       exit
                       !
                     endif
                     !
                   enddo
                   !
                 !enddo
                 !
               endif
               !
               tranfreq = energyf-energyi
               !
               itermi_ = itermi
               itermf_ = itermf
               !
               gg = real((2*Jf+1),8) !  gns(igammaf)
               !
               j = -1
               !
               if (tranfreq<0) then 
                 !
                 tranfreq = -tranfreq
                 !
                 j = itermi_
                 itermi_ = itermf_
                 itermf_ = j
                 !
                 gg = real((2*Ji+1),8)
                 !
               endif 
               !
               abscoef= (tranfreq)**3*ss(ie1)*detosc*todebye**2/gg
               !
               Nlines(max(Ji,Jf))= Nlines(max(Ji,Jf)) + 1
               !
               ! if the match is not found - stop here
               !
               !if (itermf==81939.and.itermi==81938.or.(itermi==81939.and.itermf==81938)) then 
               !  continue 
               !endif
               !
               if (itermi_<=0.or.itermf_<=0) then 
                 !
                 if (max(energyf,energyi)>energy_tmp.and.max(ee1(ie1),ee2(ie2))<enercutoff) then 
                   !
                   print("( 'cannot find a match in the energy-file for the ',i,'th transition between levels ',i8,' and ',i8)"), i,ie2,ie1
                   print("( 'Energies: ',2f16.7,' J = ',2i6,' gamma = ',2i4 )"), energyf,energyi,Jf,Ji,igammaf,igammai
                   !
                   energy_tmp = max(energyf,energyi)
                   !
                   if (itermi_<=0.and.j<0.or.(itermf_<=0.and.j>-1)) then
                        print("( 'Energy-1: ',f16.7,' J = ',i6,' gamma = ',i4,e20.12 )"), energyi,Ji,igammai,ee1(ie1)
                   else
                     print("( 'Energy-2: ',f16.7,' J = ',i6,' gamma = ',i4,e20.12 )"), energyf,Jf,igammaf,ee2(ie2)
                   endif
                   !
                   !stop 'cannot find a match in the energy-file'
                   !
                   itermf_ = -1 ; itermi_ = -1
                   !
                   cycle
                   !
                 endif
                 !
               endif
               !
               ! we use this option for large line lists where the transitions from different ranges +100 cm-1
               ! are stored in different files fort.10xxx. These are intermidiate files where the transitions 
               ! are not sorted according the frequencies, but contain a frequency collumn for this sorting. 
               ! For small line lists comment out the next line:
               !
               unit_f = 5001+int(tranfreq/1000)
               !
               ilines = ilines + 1
               !
               write(unit_f,"(2i12,2x,es16.8,2x,f16.6)"),itermf_,itermi_,abscoef,tranfreq
               !
            enddo
            !
         enddo
         !
         deallocate(ee1,ee2,ss) 
         !
         cycle
      20  exit
      
      21  write(6,*) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2
      stop
      
      enddo
      !
      if (imis /= 0) print("( 'missassigned idescr = ',i)"), imis
      !
    16 cycle 
      !
    enddo 
    !
    do i = 0,jmax
     print("(i4,2x,i12)"), i,Nlines(i)
    enddo
    !
    deallocate(sym,energies,Jktau)
    !
    print("( 'Done!')")
    !
    close(1)
    close(13)

end program linelist_DVR3D



    subroutine get_energies_from_fort14(iunit,jrot,kpar,ipar,meval,lmax,jmax,ierr,enercutoff,enercalc,ef,qpar)
      !
      implicit none
      !
      integer(4),intent(in) :: iunit
      integer(4),intent(out)   :: jrot,ipar,kpar
      integer(4),intent(out)  :: meval(2),ierr
      integer(4)              :: i_t(7),isym,iasym,verbose=4,tmp
      integer(4)              :: nsym(0:1,0:1),nasym(0:1,0:1),i2,i,rvsym(0:1,0:1,-2:2),nmax
      !
      character(len=40)     :: label
      character(1), parameter :: efpar(2)=(/'e','f'/),ppar(2)=(/'+','-'/)
      !
      integer(4),intent(in)  :: lmax,jmax
      !
      real*8,intent(in),optional   :: enercutoff
      real*8,intent(inout),optional   :: enercalc(lmax,0:jmax,4)
      character(1),intent(inout),optional :: ef(lmax,0:jmax,4),qpar(lmax,0:jmax,4)
      !
      real*8 :: ener_t,ener(4)
      !
      ! symmetrical rot. basis functions 
      !
      
      ierr = 0 
      !
      nsym(0,0) = 1
      nsym(0,1) = 4
      nsym(1,0) = 3
      nsym(1,1) = 2
      !
      ! asymmetrical rot. basis functions 
      !
      nasym(0,0) = 3
      nasym(0,1) = 2
      nasym(1,0) = 1
      nasym(1,1) = 4
      !
      nsym (0,0) = 1
      nsym (0,1) = 4
      !
      nasym(0,0) = 3
      nasym(0,1) = 2
      !
      ! asymmetrical rot. basis functions 
      !
      nsym (1,0) = 4
      nasym(1,0) = 2
      !
      nsym (1,1) = 1
      nasym(1,1) = 3
      !
      rvsym(0,0,0) = 1
      rvsym(0,1,0) = 2
      rvsym(0,0,1) = 4
      rvsym(0,1,1) = 3
      !
      rvsym(1,0,0) = 3
      rvsym(1,1,0) = 4
      rvsym(1,0,1) = 2
      rvsym(1,1,1) = 1
      !
      !read(iunit,"(a40)",end=31) label
      !print*,"label",label
      !
      meval = 0
      !kronig=0
      kpar=0
      read(iunit,*,end=31) i_t(1:6)
      !print*,"MP i_t",i_t(1:6)
      !
      jrot = i_t(1)
      !
      if (jrot==0) then 
        !
        ipar = 0 ; if (i_t(4)/=0) ipar = 1
        kpar =  i_t(2)
        !
      else
        !
        ipar = 0 ; if (i_t(4)/=0) ipar = 1
        kpar =  i_t(2)
        !
      endif
      !
      !isym  = nsym(mod(jrot+2,2),ipar)
      !iasym = nasym(mod(jrot+2,2),ipar)
      !
      isym = rvsym(mod(jrot+2,2),kpar,ipar)
!      print*,"banane",isym,kpar,ipar,mod(jrot+2,2),rvsym(mod(jrot+2,2),kpar,ipar)
      iasym = rvsym(mod(jrot+2,2),mod(kpar+3,2),mod(ipar+3,2))
!      print*,"lampone",iasym,kpar,ipar,mod(kpar+3,2),mod(ipar+3,2),rvsym(mod(jrot+2,2),mod(kpar+3,2),mod(ipar+3,2))
      !
      ! Number of found solutions
      !
      !kronig(1)=i_t(2)
      meval(1) = i_t(6)
      !
      !
      ! Read the energies from fort.14 
      !
      if (present(enercalc)) then 
        !
        if (size(enercalc(:,jrot,isym))<meval(1)) then 
          !
          write(6,"('get_energies_from_fort14: size of enercalc(:) is too small, meval = ',i8,' maxener =  ',i8)") meval(1),size(enercalc(:,jrot,isym))
          stop 'enercalc is too small'
          !
        endif 
        !
        read(iunit,*) enercalc(1:meval(1),jrot,isym)
        ef(1:meval(1),jrot,isym)=efpar(kpar+1)
        if (mod(jrot+kpar,2) == 0) then 
                qpar(1:meval(1),jrot,isym)=ppar(1)
        else if(mod(jrot+kpar,2)==1) then
                qpar(1:meval(1),jrot,isym)=ppar(2)
        else 
               stop "Error in parity assigment"
        end if 

        !
        nmax = maxloc(enercalc(:,jrot,isym),dim=1,mask=enercalc(:,jrot,isym).le.enercutoff)
        meval(1) = min(nmax,meval(1))
        !
        if (jrot>0) then
          !
          read(iunit,*) i_t(1:6) !read(iunit,"(7i4)",end=14) i_t(1:7)
          !
          ! Number of found solutions
          !
          tmp=i_t(2)
          meval(2) = i_t(6)
          !
          ! Read the energies from fort.14 
          !
          read(iunit,*) enercalc(1:meval(2),jrot,iasym)
          ef(1:meval(2),jrot,iasym)=efpar(tmp+1)
          if (mod(jrot+tmp,2) == 0) then
                qpar(1:meval(2),jrot,iasym)=ppar(1)
          else if(mod(jrot+tmp,2)==1) then
                qpar(1:meval(2),jrot,iasym)=ppar(2)
          else
               stop "Error in parity assigment, second kpar"
          end if
 
          !
        endif 
        !
        nmax = maxloc(enercalc(:,jrot,iasym),dim=1,mask=enercalc(:,jrot,iasym).le.enercutoff)
        meval(2) = min(nmax,meval(2))
        !
      else
        !        
        i2 = mod(meval(1),4)
        !print*, "i2",i2,meval(1)!
        do i = 1,meval(1)-i2,4
          read(iunit,*) ener(1:4)
        enddo
        !print*,"HYPERION"
        !
        if (i2>0) then  
          !
          read(iunit,*) ener(1:i2)
          !
          !print*, "HYPERION II"
        endif
        !
        if (jrot>0) then
          !
          read(iunit,*) i_t(1:6) !read(iunit,"(7i4)",end=14) i_t(1:6)
          !print*,"MP I_T2", i_t(1:6)
          !
          ! Number of found solutions
          !
          !kronig(2)=i_t(2)
          meval(2) = i_t(6)
          !
        endif 
        !
      endif
      !
      !print*, "FINEVOLE",kronig
      return
      !
!   14 continue 
!      meval(2) = 0
!      return 
      !
    31 continue 
      !
      ierr = 1
      return
      !
    end subroutine get_energies_from_fort14
    !
    !
    integer(4) function igamma(J,p,q)
      !
      implicit none
      !
      integer(4),intent(in) :: J,p,q
      integer(4)            :: kmin,ipar
         !
         !p = 0 ; if (kmin==0) p = 1
         !
         !p = kmin
         !q = ipar
         !
         if ( mod(J,2)==0 ) then 
           !
           if (p==0.and.q==0) then
             igamma = 1
           elseif (p==1.and.q==0) then
             igamma = 2
           elseif (p==0.and.q==1) then
             igamma = 4
           elseif (p==1.and.q==1) then
             igamma = 3
           endif 
           !
         else
           !
           if (p==0.and.q==0) then
             igamma = 3
           elseif (p==1.and.q==0) then
             igamma = 4
           elseif (p==0.and.q==1) then
             igamma = 2
           elseif (p==1.and.q==1) then
             igamma = 1
           endif 
           !
         endif
        
        
    end function igamma     
