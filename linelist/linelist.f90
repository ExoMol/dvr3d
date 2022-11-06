program linelist_DVR3D
 use input 
 implicit none

 integer :: iut
 integer, parameter :: wl=500,max_input_lines=500000,cl=80  
 character(len=wl) :: large_fmt, line_buffer
 character(len=cl) :: w,enfilename(400),intfilename(400),label(400),intfilenameout,enrfilenameout,testname
 character(len=4)  :: testunit
 integer j,jf,ji,ilevel,ilevelf,ileveli,igamma,igammaf,igammai,info,j0,i,k,jgamma,kminf,kmini,iparf,ipari,idx,igamma_asym,ic
 integer igamma_r2, igamma_r1,jdx,file_size
 integer nfiles, jmax, ilines,nenfiles,ifile,nlen
 integer nlevels,v_(9),numj(0:200,1:4),itermi,itermf,symv,sym_t,J_t,v_norm(9),igammar,imis,itermi_,itermf_

 real(8),allocatable :: energies(:,:,:),ee1(:),ee2(:),ss(:)
 integer,allocatable :: sym(:),Jktau(:,:),NN(:,:,:),Nlines(:)

 logical :: energeyfile_do = .false.,record_segmented = .false.,presence 

 real(8) :: acoef,abscoef,energy,energyf,energyi,tranfreq,linestr,ZPE,coeff,ztran,xtran,dipole,gg,enermax=1e5,enercutoff=1e5

 character(4) gamma,gammaf,gammai,gammav,jchar,gchar,gammar
 character(cl) :: enrfilename,dscrfilename,filename
 character(cl) :: record_type = " "
 character(1)  :: ch_t1
 character(len=8)  :: ch_t8
 character(200) :: ch_t200
 character(400) :: form
 character(1), allocatable :: ef(:,:,:),qpar(:,:,:)
 integer        :: unit_f,ilevelmax
 character(5)   :: symm(4) = (/'A1','A2','B1','B2'/)
 character(len=40),allocatable :: labx(:)
 integer(4)     :: gns(4) ,dwl 
 real(8),parameter :: e_thresh = 1.00d-1 ,tocm =2.19474624d+05,todebye=2.5417662d0 ,error=0.1d0
 integer(4),parameter :: itra = 13
 integer(4)     :: jrot,kpar,ipar,irottau,meval(2),nsym(0:1,0:1),nasym(0:1,0:1),iasym,isym,Lmax,is,rvsym(0:1,0:1,-2:2)
 integer(4)     :: j1,j2,kmin1,kmin2,neval1,neval2,idia,ibase1,ibase2,ie1,ie2,ipar1,ipar2,p1,p2,j_,nt_lines
 integer(4)     :: ierr,start
 real(8)        :: gz,energy_tmp
 logical        :: eof
 logical        :: zembed,asym,transitions,uvvis
 character(5)      :: symtype
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
  subroutine get_energies_from_fort14(iunit,asym,zembed,jrot,kpar,ipar,meval,lmax,jmax,ierr,enercutoff,enercalc,ef,qpar)
    integer(4),intent(in) :: iunit
    integer,intent(out)   :: jrot,ipar,kpar
    integer,intent(out)  :: meval(2),ierr
    integer(4),intent(in)  :: lmax,jmax
    real*8,intent(in),optional   :: enercutoff
    character(1),intent(inout), optional :: ef(lmax,0:jmax,4),qpar(lmax,0:jmax,4)
    real*8,intent(inout),optional   :: enercalc(lmax,0:jmax,4) 
   logical,intent(in) :: asym,zembed
  end subroutine get_energies_from_fort14
end interface

zembed=.true.
transitions=.true.
uvvis=.false.
  

! Begin input reading 
iut=5

open(unit=iut, status='scratch', action='readwrite')
write(large_fmt, '(A,i0,A)') '(A', wl, ')'
trans_loop: do i=1, max_input_lines
    read(unit=*,fmt=large_fmt,iostat=ierr) line_buffer
    if(ierr /=0) exit trans_loop
    write(iut, '(a)') trim(line_buffer)
    ! This is a hack; I need to know if to echo the input or not
    ! before processing it
    ! The option 'do_not_echo_input' is dealt with here as a special
    ! case
    line_buffer = adjustl(line_buffer) ! remove leading spaces

    do j=1, len(trim(line_buffer)) ! convert to uppercase
       ic = ichar( line_buffer(j:j))
       if( ic >= 97) line_buffer(j:j) = achar(ic-32)
    enddo
enddo trans_loop
rewind(iut)

do

    call read_line(eof,iut) ; if (eof) exit
    call readu(w)
print*,w
    select case(w)
    case ("UVVIS")
        uvvis=.true.
        write(*,*) "Multiple electronic states are considered"
!        !
    case ("NSTATES")
!          !
        call readi(i)
        Nenfiles=i
          if (Nenfiles>400) then 
             print('("Too many state files (>400):",i)'),Nenfiles
             stop 'Too many state files'
          endif
          if (uvvis==.false.) then 
            do i = 1,Nenfiles
               read(iut,*) enfilename(i)
            enddo
          else 
            do i = 1,Nenfiles
               read(iut,*) enfilename(i),label(i)
            enddo
          end if 
!          ! 
    case ("NTRANS")
        call readi(i)
         nfiles=i
         if (nfiles>400) then 
             print('("Too many trans files (>400):",i)'),nfiles
             stop 'Too many state files'
         endif
         do i = 1,nfiles
             read(iut,*) intfilename(i)
            enddo
    case ("OUT_STATE")
        call readu(enrfilenameout)
        !read(iut,*) enrfilenameout
    case ("OUT_TRANS")
        call readu(intfilenameout)
        !read(iut,*) intfilenameout

    case ("SYMMETRY")
        call readu(symtype)
        if (symtype=="ABC") then 
           asym=.True. 
        else if (symtype=="AB2") then 
           asym=.False.
        else 
           stop "Case not included in the program"  
       end if 
    case ("ENERMAX")
       print*,enermax 
       call readf(enermax)
       print*,enermax,2
    case ("ENCUTOFF")
       call readf(enercutoff)
       enercutoff=enercutoff/tocm
    case ("GNS")
        if (symtype=="ABC") then
          call readi(gns(1))
          call readi(gns(2))
        else if (symtype=="AB2") then
          call readi(gns(1))
          call readi(gns(2))
          call readi(gns(3))
          call readi(gns(4))
        else
           stop "The symmetry tipe shoud be defined before GNS "
       end if
    case("SEGMENTED") 
        record_segmented = .true. 
    case ("NOTRANS")
        transitions=.false.
        print*, "Transitions are not considered."
    end select  
 
end do  


! End input reading


if (record_segmented) then
       open(unit=itra,file=trim(intfilename(1)),form='unformatted',recordtype='segmented') ! H2O
else  
       open(unit=itra,file=trim(intfilename(1)),form='unformatted') ! H2S
end if  
      read(itra,err=210) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2
      210  write(6,*) "zembed=",zembed

rewind(itra)
close(itra)

i = 0
numj = 0
jmax = 0
nlevels = 0
!
! Define the molec. symmetry group correlatation with the internal DVR3D symmetry indeces.
! symmetrical rot. basis functions 
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
  call get_energies_from_fort14(1,asym,zembed,jrot,kpar,ipar,meval,lmax,jmax,ierr)
  !
  !
  if (asym==.true.) then
    isym  = igamma_asym(jrot,mod(kpar,2))
    iasym = igamma_asym(jrot,mod(kpar+1,2))
  else if (asym==.false. .and. zembed==.true.) then         
    isym  = igamma_r2(jrot,kpar,ipar)
    iasym = igamma_r2(jrot,mod(kpar+1,2),ipar)
  else 
    isym  = igamma_r1(jrot,kpar,ipar)
    iasym = igamma_r1(jrot,mod(kpar+1,2),ipar)
  end if 
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
allocate(sym(nlevels),energies(Lmax,0:jmax,4),Jktau(nlevels,3),NN(Lmax,0:jmax,4),ef(Lmax,0:jmax,4),qpar(Lmax,0:jmax,4), & 
         labx(nlevels),stat=info); 
if (info/=0) stop 'error: sym,Jktau are out of memory'
!
sym = 0
energies = 0
NN = 0
nlevels = 0 
numj = 0
start=1
!
print("('Reading energies...')")
!
! go through all energy output files:
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
  call get_energies_from_fort14(1,asym,zembed,jrot,kpar,ipar,meval,Lmax,jmax,ierr,enercutoff,energies,ef,qpar)
  !
  !
  if (asym==.true.) then 
    isym  = igamma_asym(jrot,mod(kpar,2))  
    iasym = igamma_asym(jrot,mod(kpar+1,2)) 
else if (asym==.false. .and. zembed==.true.) then         
    isym  = igamma_r2(jrot,kpar,ipar)
    iasym = igamma_r2(jrot,mod(kpar+1,2),ipar)
else 
    isym  = igamma_r1(jrot,kpar,ipar)
    iasym = igamma_r1(jrot,mod(kpar+1,2),ipar)        
  end if 
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
  print*,"Checking meval", meval
  numj(jrot,isym)  = meval(1)
  numj(jrot,iasym) = meval(2)
  !
  nlevels = nlevels + sum(meval(:))
  !
  if (uvvis==.true.) then
    labx(start:nlevels)=label(ifile)
    start=nlevels
  end if 
  !
  close(1)
  !
enddo
!
ZPE =energies(1,0,1)*tocm
!
write(6,"('ZPE = ',f15.6)") ZPE
!
energies = energies*tocm-ZPE
!
! from uppercase to lowercase:
 nlen = len(enrfilenameout)
do i=1,nlen
   ic = ichar(enrfilenameout(i:i))
   if (ic >= 65 .and. ic < 90) enrfilenameout(i:i) = achar(ic+32)
end do 
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


if (asym==.false.) then 
   print*, "Writing the state files for AB2;"
do jrot = 0,jmax
  !print*,"jrot printed:",jrot
  !
  do isym = 1,4  
    !
    jdx=1
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
        if (energeyfile_do) then
          !
          ! Here is the main print-out of the .states file in the ExoMol format  (Commented: we are reporting the char table symmetries)
          ! 
          if (uvvis==.false.) then
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a4,1x,i6,1x,a1,1x,a1,5(1x,a3))") &
              i,energies(ilevel,jrot,isym),gns(isym)*(2*jrot+1),jrot,error,symm(isym), & 
              jdx,ef(ilevel,jrot,isym),qpar(ilevel,jrot,isym),"NaN","NaN","NaN","NaN","NaN"
           jdx=jdx+1
          else 
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a4,1x,a5,1x,i6,1x,a1,1x,a1,5(1x,a3))") &
              i,energies(ilevel,jrot,isym),gns(isym)*(2*jrot+1),jrot,error,symm(isym),labx(i), & 
              jdx,ef(ilevel,jrot,isym),qpar(ilevel,jrot,isym),"NaN","NaN","NaN","NaN","NaN"
           jdx=jdx+1
          end if 
     !
        endif
    enddo
    !
    print("(i4,1x,i3,2x,2i8)"),jrot,isym,numj(jrot,isym),i
    !
  enddo
!
enddo


else 
  print*, "Writing the state files for ABC;"
do jrot = 0,jmax
  !print*,"jrot printed:",jrot
  !
  do isym = 1,4  
    !
    jdx=1 
    print*, "Symmetry ", isym," has", numj(jrot,isym), "levels"
    do ilevel = 1,numj(jrot,isym)

      !
      if (gns(isym)<1) cycle 
        !
        i = i + 1
        !
        NN(ilevel,jrot,isym) = i
        !
        !prepare and store the .states file
        if (energeyfile_do) then
          !
          ! Here is the main print-out of the .states file in the ExoMol format  (Commented: we are reporting the char table symmetries)
          if (uvvis==.false.) then
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a1,1x,a1,5(1xa3))") & 
             i,energies(ilevel,jrot,isym),gns(isym)*(2*jrot+1),jrot,error,qpar(ilevel,jrot,isym),ef(ilevel,jrot,isym),"NaN","NaN","NaN","NaN","NaN"
            jdx=jdx+1
          else 
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a5,1x,a1,1x,a1,5(1xa3))") & 
              i,energies(ilevel,jrot,isym),gns(isym)*(2*jrot+1),jrot,error,labx(i),&
              qpar(ilevel,jrot,isym),ef(ilevel,jrot,isym),"NaN","NaN","NaN","NaN","NaN"
            jdx=jdx+1
          end if   
          !
        endif
    enddo
    !
    print("(i4,1x,i3,2x,2i8)"),jrot,isym,numj(jrot,isym),i
    !
  enddo
!
enddo
end if 


!
 close(12)

if (transitions==.True.) then  
    dwl=maxval(energies(:,:,:))-minval(energies(:,:,:)) ; print*,dwl
    do i=1,int(dwl/1000)
       write(testname,'(A,A,i4,A)') trim(intfilenameout),".",5000+i,".trans" ;
        nlen = len(testname)
       do idx=1,nlen
          ic = ichar(testname(idx:idx))
          if (ic >= 65 .and. ic < 90) testname(idx:idx) = achar(ic+32)
       end do
       open(5000+i,file=testname)
    end do 
! MP The file works up to here  - here things are fine as well, but I am not properly ure 
    !
    ! The second part: .trans file
    !
    print("(' N of levels =  ',i)"),nlevels
    !
    print("('Generate the Transition file ...')")
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
         read(itra,end=20) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2         
         !write(422,*) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2
         !,"gz:", gz,"zembed:",zembed,"ibase1:",ibase1,"ibase2:",ibase2
         !
         !
         ji = j1 ; print*,"Jinitial",ji
         jf = j2 ; print*,"Jfinal  ",jf
         !if (jf > 2 ) exit
         !
         p1 = abs(1-kmin1)
         p2 = abs(1-kmin2)
         !
         allocate(ee1(neval1), ee2(neval2), ss(neval1),stat=info)
         if (info/=0) stop 'error: ee1,ee2,ss are out of memory'
         !
         read(itra) ee1 !; write(54,*)"ee1:",ee1; write(54,*)
         read(itra) ee2 !;write(54,*)"ee2:",ee2; write(54,*)
         !
         !!!!!!!!!!!!!!!
         !write(6,"(9i7,2f16.6)") j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,ee1(neval1)*tocm-ZPE,ee2(neval2)*tocm-ZPE
         !!!!!!!!!!!!!!!
         !
         !
         !
         if (asym==.true.) then 
           igammai = igamma_asym(j1,p1) ; print*,"igamma1:",igamma_asym(j1,p1),j1,p1 
           igammaf = igamma_asym(j2,p2) ; print*,"igamma2:",igamma_asym(j2,p2),j2,p2      
         else if (asym==.false. .and. zembed==.true.) then         
           igammai = igamma_r2(J1,p1,ipar1) ; print*,"igamma1:",igamma_r2(J1,p1,ipar1) 
           igammaf = igamma_r2(J2,p2,ipar2) ; print*,"igamma2:",igamma_r2(J2,p2,ipar2)
else 
           igammai = igamma_r1(J1,p1,ipar1) ; print*,"igamma1:",igamma_r1(J1,p1,ipar1) 
           igammaf = igamma_r1(J2,p2,ipar2) ; print*,"igamma2:",igamma_r1(J2,p2,ipar2)
         end if 
         !
         if ( j_ /= max(j1,j2,j_)) then
            j_ = max(j1,j2,j_)
            write(6,"('J = ',i8,4x,i10)") j_,ilines
         endif 
         !
         energy_tmp = 0
         !unit_f=5001
         !
         do ie2=1,neval2
            !
            read(itra) ss 
            !
            energyf=ee2(ie2)*tocm-ZPE 
            !
            !
            !
            print*,ie2
            if (abs(energyf-energies(ie2,jf,igammaf))<=e_thresh) then 
              !
              itermf = NN(ie2,jf,igammaf) ; print*,"itermf",NN(ie2,jf,igammaf)     
              !
            elseif ( ee2(ie2)<enercutoff ) then 
              !
              write(6,"('f:',3i7,2f16.6)") ie2,j2,igammaf,energyf,energies(ie2,j2,igammaf)
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
                 itermf = NN(j,jf,is)!;write(54,*) 
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
               !
               ! Find the i-state by using the energy threshold e_thresh (typically 1e-4)
               ! 
               itermi = 0
               !
               ! Check all states within the same J and symmetry to find the energy match and store the number to itermi
               !
               ilevel = ie1+numj(ji,igammai)-neval1
               !
               !write(703,*),energyi,energies(ie1,ji,igammai),ie1,ilevel,numj(ji,igammai),neval1
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
                 do is = 1,4
                   !
                   !is = igammai
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
                 enddo
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
               if (itermi_ == 1 .or. itermf_==1 ) write(1024,*) itermf_,itermi_,tranfreq,ie1,abscoef,gg
               !
               Nlines(max(Ji,Jf))= Nlines(max(Ji,Jf)) + 1
               !
               !
               if (itermi_<=0.or.itermf_<=0) then 
                 !
                 if (max(energyf,energyi)>energy_tmp.and.max(ee1(ie1),ee2(ie2))<enercutoff) then 
                   !

                   print("( 'cannot find a match in the energy-file for the ',i,'th transition between levels ',i8,' and ',i8)"), i,ie2,ie1
                   print("( 'Energies: ',2f16.7,' J = ',2i6,' gamma = ',2i4 )"), energyf,energyi,Jf,Ji,igammaf,igammai
                   print("( 'Iterm_i : ',i7,'  Iterm_f : ',i7)"),itermi_,itermf_
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
                unit_f= 5001+int(tranfreq/1000)
                !write(testname,'(A,i4,A)') trim(intfilenameout),unit_f,".trans"
                !inquire(file=testname,exist=presence) ; print*, presence 
                !if (presence==.false. ) then 
                !   open(unit_f,file=testname,action="write")
                !   print*,unit_f, testname 
                !end if 
               !
               ilines = ilines + 1
               !
               !write(unit_f,"(2i12,2x,es16.8,2x,1f16.6)"),itermf_,itermi_,abscoef,tranfreq!, & ! energyf,energyi
               !write(unit_f,"(2i12,2x,es16.8)"),itermf_,itermi_,abscoef!,
                                                       !energies(ie1,ji,igammai),energies(ie2,j2,igammaf)
               write(unit_f,"(4es16.8)")abscoef,tranfreq,ss(ie1),gg
               !
            enddo
            !
         enddo
         !
         deallocate(ee1,ee2,ss) 
         !
         cycle
      20  exit
      
      !21  write(6,*) j1,j2,kmin1,kmin2,neval1,neval2,idia,ipar1,ipar2,gz,zembed,ibase1,ibase2
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

   do i=1,int(dwl/1000)
   INQUIRE(unit=5000+i, SIZE=file_size)
   print*, 5000+i, file_size
   if (file_size==0) then
       close(5000+i,status="delete")
   else
       close(5000+i)
   end if
   end do
else 
 print*,  "Transitions skipped" 
end if 
print("( 'Done!')")
    !
close(1)
close(13)

!do i=1,int(dwl/1000)
!  INQUIRE(unit=5000+i, SIZE=file_size)
!  print*, 5000+i, file_size  
!  if (file_size==0) then 
!      close(5000+i,status="delete")
!  else 
!      close(5000+i)
!  end if 
!end do

end program linelist_DVR3D



    subroutine get_energies_from_fort14(iunit,asym,zembed,jrot,kpar,ipar,meval,lmax,jmax,ierr,enercutoff,enercalc,ef,qpar)
      !
      implicit none
      !
      integer(4),intent(in) :: iunit
      integer(4),intent(out)   :: jrot,ipar,kpar
      integer(4),intent(out)  :: meval(2),ierr
      integer(4)              :: i_t(7),isym,iasym,verbose=4,tmp
      integer(4)              :: nsym(0:1,0:1),nasym(0:1,0:1),i2,i,rvsym(0:1,0:1,-2:2),nmax,igamma,igamma_asym, &
                                 igamma_r1,igamma_r2,idx 
      !
      character(len=40)     :: label
      character(1), parameter :: efpar(2)=(/'e','f'/),ppar(2)=(/'+','-'/),op(4)=(/"o","o","p","p"/)
      !
      integer(4),intent(in)  :: lmax,jmax
      !
      real*8,intent(in),optional   :: enercutoff
      real*8,intent(inout),optional   :: enercalc(lmax,0:jmax,4)
      character(1),intent(inout),optional :: ef(lmax,0:jmax,4),qpar(lmax,0:jmax,4)
      logical,intent(in) :: asym,zembed
      !
      real*8 :: ener_t,ener(4)
      !
      ! symmetrical rot. basis functions 
      !
      ierr = 0 
      !
      !
      !read(iunit,"(a40)",end=31) label
      !print*,"label",label
      !
      meval = 0
      kpar=0
      read(iunit,*,end=31) i_t(1:6)
      !
      !
      jrot = i_t(1) ; print*, "Jrot is ", jrot
      !
      ipar = 0 ; if (i_t(4)/=0) ipar = 1
      kpar =  i_t(2)
      !
      !
      !
      if (asym==.true.) then 
         isym  = igamma_asym(jrot,mod(kpar,2))  
         iasym = igamma_asym(jrot,mod(kpar+1,2)) 
      else if (asym==.false. .and. zembed==.true.) then        
         isym  = igamma_r2(jrot,kpar,ipar)
         iasym = igamma_r2(jrot,mod(kpar+1,2),ipar)
      else 
         isym  = igamma_r1(jrot,kpar,ipar)
         iasym = igamma_r1(jrot,mod(kpar+1,2),ipar)        
      end if 
      !
      ! Number of found solutions
      !
      meval(1) = i_t(6) ; print*,"Number of energy levels read, for ipar=0: ",meval(1)
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
        
        read(iunit,*) enercalc(1:meval(1),jrot,isym) !; print*, "Controllo 1:", size(enercalc(1:meval(1),jrot,isym))
        
        !ef(1:meval(1),jrot,isym)=efpar(kpar+1)   it prints e/f original version
        ef(1:meval(1),jrot,isym)=op(isym)
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
          read(iunit,*) i_t(1:6) !
          !
          ! Number of found solutions
          !
          tmp=i_t(2)
          meval(2) = i_t(6)
          !
          ! Read the energies from fort.14 
          !
          read(iunit,*) enercalc(1:meval(2),jrot,iasym)
          ef(1:meval(2),jrot,iasym)=op(iasym)
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
        !
        if (i2>0) then  
          !
          read(iunit,*) ener(1:i2)
          !
        endif
        !
        if (jrot>0) then
          !print*,"passo da qui"; do idx=1,6
          !read(iunit,*)  i_t(idx) ; print*, "LAMERDA",i_t(idx)
         ! end do 
          read(iunit,*) i_t(1:6) !read(iunit,"(7i4)",end=14) i_t(1:6)
          !
          ! Number of found solutions
          !
          meval(2) = i_t(6); print*,"Number of energy levels read, for ipar=0:",meval(2)
          !
        endif 
        !
      endif
      !
      return
      !
      !
    31 continue 
      !
      ierr = 1
      return
      !
    end subroutine get_energies_from_fort14
    !
    !
     integer(4) function igamma_r2(J,p,q)
      ! Valid for zembed=.True.  z//R
      ! kpar=0 -> "e" ; kpar=1 -> "f"
      ! ipar=0 -> "ortho" ; ipar=1 "para"
      ! There is a flipping between the J even and J odd
      implicit none
      !
      integer(4),intent(in) :: J,p,q
      integer(4)            :: kmin,ipar
         !
         !
         if ( mod(J,2)==0 ) then 
           !
           if (p==0.and.q==0) then ! A1 
             igamma_r2 = 1
           elseif (p==1.and.q==0) then !  A2
             igamma_r2 = 2
           elseif (p==0.and.q==1) then !  B2
             igamma_r2 = 4
           elseif (p==1.and.q==1) then ! B1
             igamma_r2 = 3
           endif 
           !
         else
           !
           if (p==0.and.q==0) then ! A2
             igamma_r2 = 2
           elseif (p==1.and.q==0) then ! A1
             igamma_r2 = 1 
           elseif (p==0.and.q==1) then ! B1 
             igamma_r2 = 3
           elseif (p==1.and.q==1) then ! B2
             igamma_r2 = 4
          endif 
           !
         endif        
    end function igamma_r2    
 
    integer(4) function igamma_r1(J,p,q)
      ! Valid for zembed=.True.  z//r
      ! kpar=0 -> "e" ; kpar=1 -> "f"
      ! ipar+J=0 -> "ortho" ; ipar+J=1 "para"
      ! There is a flipping between the J even and J odd
      implicit none
      !
      integer(4),intent(in) :: J,p,q
      integer(4)            :: kmin,ipar
         !
         !
         if ( mod(J,2)==0 ) then 
           !
           if (p==0.and. q==0) then ! A1 
             igamma_r1 = 1
           elseif (p==1.and. q==0) then !  A2
             igamma_r1 = 3
           elseif (p==0.and. q==1) then !  B2
             igamma_r1 = 4
           elseif (p==1.and. q==1) then ! B1
             igamma_r1 = 2
           endif 
           !
         else
           !
           if (p==0.and. q==0) then ! B1
             igamma_r1 = 3
           elseif (p==1.and. q==0) then ! B2
             igamma_r1 = 1 
           elseif (p==0.and. q==1) then ! A1 
             igamma_r1 = 2
           elseif (p==1.and. q==1) then ! A2
             igamma_r1 = 4
           endif 
           !
         endif        
    end function igamma_r1     
 


    integer(4) function igamma_asym(J,p)
      !
      implicit none
      !
      integer(4),intent(in) :: J,p
         !
         !
         if ( mod(J+p,2)==0 ) igamma_asym = 1
         if ( mod(J+p,2)==1 ) igamma_asym = 2
    end function igamma_asym

