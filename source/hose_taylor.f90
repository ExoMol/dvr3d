program hosetaylor
! THIS PROGRAM CALCULATES |<PSI|PSI>|^2 FOR EACH VALUE OF K. 
! IT THEN PRINTS OUT THE LARGEST COMPONENT AND THE PREDICTED VALUE OF Kc AND Kc
! CURRENTLY BASED ON 'WFNREAD.F90' AND THIS PROGRAM IS OPERATING CORRECTLY FOR RADAU COODINATES. 
!
! USER MUST SUPPLY EZERO!!!!!!!!
!
! OUTPUT FORMAT IS: J, ENERGY, KA, KC, MAX(|<PSI|PSI>|^2),  IPAR, KMIN, MOD(QUANTA OF NU3,2)
!
  common/size/ idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval,jk,ifile

! file is assumed to be attached to unit 26
  ifile=26
  open(unit=ifile,form='unformatted',recordtype='segmented')

! reader header record to determine file structure  
  read(ifile) idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval
  write(*,*) idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval
  If (jrot.eq.0) kmin=1
  jk=jrot+kmin


  If (idia .eq. -2 .and. jk .gt. 1) then
     print *, "Radau, 8 or 9"
      call read_8or9_radau
  elseif (idia .eq. -2) then
     print *, "Radau, 26"
      call read_26_radau
  else
     print *, "ERROR, MUST USE RADAU"
  endif

end program hosetaylor



!###################################################################

subroutine read_26_radau
  implicit double precision (a-h,o-y), logical (z)
  common/size/ idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval,jk,ifile

  integer :: ifile, jk
  integer, allocatable ::  nbass(:), lmin(:), lbass(:), iv(:)
  double precision, allocatable:: e(:), d(:)
  dimension xm(3)
double precision :: component,ezero
character(len=4) :: state
integer :: parity

  read(ifile) zembed,zmors1,zmors1,xmass,g1,g2,zncor,zquad2
  print *, zembed,zmors1,zmors1,xm,g1,g2,zncor,zquad2

  read(ifile) re1,diss1,we1,re1,diss1,we1
  print *, re1,diss1,we1,re1,diss1,we1

! PLEASE SPECIFY EZERO TO OBTAIN ENERGIES WRT ZERO POINT EQUILIBRIUM
  ezero = 0.0d0

  if(ezero .eq. 0.0d0) print *, "ZPE IS ZERO"


  do i=1,5
     read(ifile)
  enddo
  
  read(ifile) iang2,ibass2
!  print *, "iang2,ibass2 ",  iang2,ibass2
  read(ifile)

  read(ifile) meval2
!  print *, "meval2 ",  meval2

!get the enrgy levels
  allocate( e(meval2) )
  call getrow(e,meval2,ifile)
!  print *, "energy levels in Hartrees"
!  print *, e


! D is consists of the values of the wavefunction as expressed on the
! the DVR grid in the 3 coordinates, writen as a 3d array d(NALF,
! NPNT1, NPNT2) NALF, NPNT1, NPNT2 will give the location on the DVR
! grid for that value of the wavefunction.

!get the eigenvectors

if ( (MODULO(jrot,2) .eq. 0 ) ) THEN

    if( (ipar .eq. 0) .and. (kmin .eq. 1) ) then
    state = 'para'
    label = 0
    parity = 1
    else if( (ipar .eq. 0) .and. (kmin .eq. 0) ) then
    state = 'para'
    label = 0
    parity = -1
    else if( (ipar .eq. 1) .and. (kmin .eq. 1) ) then
    state = 'orth'
    label = 1
    parity = 1
    else if( (ipar .eq. 1) .and. (kmin .eq. 0) ) then
    state = 'orth'
    label = 1
    parity = -1
    else
    continue 
    end if

else if ( (MODULO(jrot,2) .eq. 1 ) ) THEN

    if( (ipar .eq. 0) .and. (kmin .eq. 1) ) then
    state = 'orth'
    label = 1
    parity = -1
    else if( (ipar .eq. 0) .and. (kmin .eq. 0) ) then
    state = 'orth'
    label = 1
    parity = 1
    else if( (ipar .eq. 1) .and. (kmin .eq. 1) ) then
    state = 'para'
    label = 0
    parity = -1
    else if( (ipar .eq. 1) .and. (kmin .eq. 0) ) then
    state = 'para'
    label = 0
    parity = 1
    else 
    continue
    end if
else
continue
end if

if (kmin .eq. 0) ka =  1 
if (kmin .eq. 1) ka = 1 - 1

if(ka .eq. 0) then 
kc=jrot
!KA+KC EVEN
if((mod((ka + kc),2) .eq. 0) .and. (state .eq. 'orth')) nu3=1
if((mod((ka + kc),2) .eq. 0) .and. (state .eq. 'para')) nu3=0

if((mod((ka + kc),2) .eq. 1) .and. (state .eq. 'orth')) nu3=0
if((mod((ka + kc),2) .eq. 1) .and. (state .eq. 'para')) nu3=1

else if(ka .gt. 0 ) then
kc1 = jrot - ka
p1=(-1.0d0)**kc1
!----------
kc2 = jrot + 1 - ka
p2=(-1.0d0)**kc2



if(parity .eq. p1) kc=kc1
if(parity .eq. p2) kc=kc2
if((mod((ka + kc),2) .eq. 0) .and. (state .eq. 'orth')) nu3=1
if((mod((ka + kc),2) .eq. 0) .and. (state .eq. 'para')) nu3=0

if((mod((ka + kc),2) .eq. 1) .and. (state .eq. 'orth')) nu3=0
if((mod((ka + kc),2) .eq. 1) .and. (state .eq. 'para')) nu3=1
else
continue 
end if


  allocate( d(ibass2) )
  do j=1, meval2
     call getrow(d,ibass2,ifile)


d=d**2

!WAVE FUNCTION FOR EACH LEVEL IS READ AND SUMMED. THEN REPEAT OVER EACH BLOCK.
component=0.0d0
do i=1,ibass2
component=component+d(i)
end do

write(9,"(i2,2x,f12.5,2x,i2,1x,i2,2x,f8.5,2x,i1,2x,i1,2x,i1)") jrot,(e(j)*2.19474624d+05 - ezero),ka,kc,component,label,kmin,nu3

  enddo

end subroutine read_26_radau

!###################################################################

subroutine read_8or9_radau
  implicit double precision (a-h,o-y), logical (z)
  common/size/ idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval,jk,ifile

  integer :: ifile, jk
  integer, allocatable ::  nbass(:), lmin(:), lbass(:)
  double precision, allocatable:: e(:), d(:)
  dimension xm(3)
double precision, allocatable :: sum(:), component(:,:)
character(len=4) :: state
integer :: parity
double precision,allocatable :: biggest(:)
integer, allocatable :: ka(:),kc(:),kc1(:),kc2(:),p1(:),p2(:),nu3(:)
double precision :: ezero



! PLEASE SPECIFY EZERO TO OBTAIN ENERGIES WRT ZERO POINT EQUILIBRIUM
  ezero = 0.0d0

  if(ezero .eq. 0.0d0) print *, "ZPE IS ZERO"

  read(ifile) zembed,zmorse1,zmorse2,xm,g1,g2,zncor
  write(*,*) zembed,zmorse1,zmorse2,xm,g1,g2,zncor

  read(ifile) re2,diss2,we2,re2,diss2,we2
  write(*,*) re2,diss2,we2,re2,diss2,we2


  allocate( nbass(jk),lmin(jk),lbass(jk) )
  read(ifile) mbass0,lmin,lbass,nbass
  write(*,*) mbass0,lmin,lbass,nbass

  read(ifile)

  !read in the energies, Hartrees
  read(ifile) neval
  write(*,*) neval 
  allocate( e(neval) )
  read(ifile) e


! the array d consists of neval eigenvectors. Thus d(1:nbass(k))
! represents the 1st eigenvector; d( nbass(k)+1:nbass(k) ) the second,
! etc. each eigenvector consists of the values of the wavefunction as
! expressed on the DVR grid in the 3 coordinates, written as a 3D
! array. thus if the eigenvector is placed an array vec(NALF, NPNT1,
! NPNT2), NALF, NPNT1, NPNT2 will give the location on the DVR grid for
! that value of the wavefunction.

  !read in wavefunctions, each k block is read in separtly
!!!!!!
!EACH BLOCK READ IN, ONE AT A TIME
!!!!!!

allocate(sum(neval), component(neval,jk),biggest(neval),ka(neval),kc(neval),kc1(neval),&
kc2(neval),p1(neval),p2(neval),nu3(neval))


do i=1,neval
sum(j)=0.0d0
biggest(j)=0.0d0
end do

if ( (MODULO(jrot,2) .eq. 0 ) ) THEN

    if( (ipar .eq. 0) .and. (kmin .eq. 1) ) then
    state = 'para'
    label = 0
    parity = 1
    else if( (ipar .eq. 0) .and. (kmin .eq. 0) ) then
    state = 'para'
    label = 0
    parity = -1
    else if( (ipar .eq. 1) .and. (kmin .eq. 1) ) then
    state = 'orth'
    label = 1
    parity = 1
    else if( (ipar .eq. 1) .and. (kmin .eq. 0) ) then
    state = 'orth'
    label = 1
    parity = -1
    else
    continue 
    end if

else if ( (MODULO(jrot,2) .eq. 1 ) ) THEN

    if( (ipar .eq. 0) .and. (kmin .eq. 1) ) then
    state = 'orth'
    label = 1
    parity = -1
    else if( (ipar .eq. 0) .and. (kmin .eq. 0) ) then
    state = 'orth'
    label = 1
    parity = 1
    else if( (ipar .eq. 1) .and. (kmin .eq. 1) ) then
    state = 'para'
    label = 0
    parity = -1
    else if( (ipar .eq. 1) .and. (kmin .eq. 0) ) then
    state = 'para'
    label = 0
    parity = 1
    else 
    continue
    end if
else
continue
end if


do k=1, jk
allocate( d(nbass(k)*neval) )
!print *, "K block ", k-1
read(ifile) d
d=d**2

!WAVE FUNCTION FOR EACH LEVEL IS READ AND SUMMED. THEN REPEAT OVER EACH BLOCK.
    do j=1, neval
    component(j,k)=0.0d0
        do i = (j-1)*nbass(k)+1, nbass(k)*(j)
        sum(j) = sum(j) + d(i)
        component(j,k) = component(j,k) + d(i) 
        end do
            if(k .eq. 1) then 
            biggest(j)=component(j,k)
            if (kmin .eq. 0) ka(j) =  k 
            if (kmin .eq. 1) ka(j) = k - 1
            else if ( (k .gt. 1) .and. ( component(j,k) .gt. biggest(j) )) then
            biggest(j)=component(j,k)
            if (kmin .eq. 0) ka(j) = k 
            if (kmin .eq. 1) ka(j) = k - 1
            else
            continue
            end if

if(ka(j) .eq. 0) then 
kc(j)=jrot
!KA+KC EVEN
if((mod((ka(j) + kc(j)),2) .eq. 0) .and. (state .eq. 'orth')) nu3(j)=1
if((mod((ka(j) + kc(j)),2) .eq. 0) .and. (state .eq. 'para')) nu3(j)=0

if((mod((ka(j) + kc(j)),2) .eq. 1) .and. (state .eq. 'orth')) nu3(j)=0
if((mod((ka(j) + kc(j)),2) .eq. 1) .and. (state .eq. 'para')) nu3(j)=1

else if(ka(j) .gt. 0 ) then
kc1(j) = jrot - ka(j)
p1(j)=(-1.0d0)**kc1(j)
!----------
kc2(j) = jrot + 1 - ka(j)
p2(j)=(-1.0d0)**kc2(j)



if(parity .eq. p1(j)) kc(j)=kc1(j)
if(parity .eq. p2(j)) kc(j)=kc2(j)
if((mod((ka(j) + kc(j)),2) .eq. 0) .and. (state .eq. 'orth')) nu3(j)=1
if((mod((ka(j) + kc(j)),2) .eq. 0) .and. (state .eq. 'para')) nu3(j)=0

if((mod((ka(j) + kc(j)),2) .eq. 1) .and. (state .eq. 'orth')) nu3(j)=0
if((mod((ka(j) + kc(j)),2) .eq. 1) .and. (state .eq. 'para')) nu3(j)=1
else
continue 
end if

    end do
deallocate( d )
  enddo

do j=1, neval
write(9,"(i2,2x,f12.5,2x,i2,1x,i2,2x,f8.5,2x,i1,2x,i1,2x,i1)") jrot,(e(j)*2.19474624d+05 - ezero),ka(j),kc(j),biggest(j),label,kmin,nu3(j)
end do




end subroutine read_8or9_radau

!#####################################################################

subroutine getrow(row,nrow,iunit)
  implicit double precision (a-h,o-y)
  double precision, dimension(nrow) :: row      
  read(iunit) row
  return
end subroutine getrow
