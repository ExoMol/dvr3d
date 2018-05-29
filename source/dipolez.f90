
!MODULE DEFINITIONS--------------------------------------------

module input
!  definition of the control input parameters
  save
  !line 1
  !namelist prt
  logical :: zprint ! =T supplies extra print out for debugging purposes
  logical :: ztra   ! =T writes data for spectra to stream ITRA
  logical :: zstart ! =T initiates the output file for the data for SPECTRA
  integer :: iket   ! input stream from DVR3DRJZ/ROTLEV3/ROTLEV3B for the ket (unformated)
  integer :: ibra   ! input stream for the bra (unformmatted)
  integer :: itra   ! output stream to SPECTRA (if ZTRA = T ) (unformated)
  namelist /prt/ zprint, ztra, zstart, iket, ibra, itra    
  
  !line 2
  character(len=72) :: title
  
  !line 3
  integer :: npot   ! number of gauss legendre integration points
  integer :: nv2    ! number of ket eigenfunctions considered
  integer :: nv1    ! number of bra eigenfunctions considered
  integer :: ibase2 ! number of lowest ket eigenfunctions skipped
  integer :: ibase1 ! number of lowest bra eigenfunctions skipped
  
  !line 4
  real(8) :: ezero  ! the ground state of the system in cm-1
end module input

module lists
!linked list definitions for the complicated header reading procedure
!this is just an exercise, feel free to get rid of this mess (but it works ;-) 
  save
  type rgridnode
     real(8) :: rgrid 
     type (rgridnode),pointer :: next
  end type rgridnode
  type eigsnode
     real(8) :: eigs  
     type (eigsnode),pointer :: next
  end type eigsnode
  
  type headnode
     integer :: n3d    ! total number of grid points 
     integer :: nr     ! number of radial grid points
     integer :: Jrot   ! rotational state
     integer :: nval   ! number of eigenvalues
     integer :: kpar   ! laboratory coordinates parity
     integer :: ifpar  ! nuclear permutation parity
     integer :: nk     ! total number of k blocks
     type(rgridnode),pointer :: rgridfirst,rgridlist 
     type(eigsnode),pointer :: eigsfirst,eigslist
     type (headnode),pointer :: next 
  end type headnode
  type(headnode),pointer :: headfirst,headlist 
end module lists

module workdata
  save
  integer :: nk1, nk2, ispar1, ispar2, kz1, kz2, nval, nval1, nval2, nth1, nth2
  integer :: jrot1,jrot2,iqpar1,iqpar2,ifpar1,ifpar2,kpar1,kpar2,n3d1,n3d2,nr1,nr2
  real(8), allocatable :: eigs1(:),eigs2(:)
  real(8), allocatable :: rgrid1(:), rgrid2(:), thgrid1(:), thgrid2(:)  
  
  ! common grid
  integer :: nr,nth
  real(8) :: x,x1,x2
  real(8), allocatable :: thgrid(:),tanth(:),rgrid(:),wt(:)
  real(8), allocatable :: waves1(:,:),waves2(:,:)
  real(8), allocatable :: dipcxa(:),dipcya(:),dipcxb(:),dipcyb(:)
  real(8), allocatable :: pol1(:,:),pol2(:,:)

  ! dipole storage
  real(8), allocatable :: dpba(:,:)!,dpbb(:,:)
  real(8) :: pi,sumpb
end module workdata
! END MODULE DEFINITIONS--------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program dipole
  use input
  use workdata
!  
  implicit none
  integer :: i,channel(2)
  
  pi=acos(-1.0d0)
  sumpb=0.d0

  write(6,*) "program started"  
 
  ! read the control file (or input file)
  call read_input(channel)

  ! read the BRA and KET wavefunction files' headers 
  ! compicated exercise using linked lists
  ! if it bothers you, feel free to change it
  call read_lists(channel)

  ! input verification
  call input_verify
  
  ! check the norms of the wavefunctions
  !  call normcheck(channel(1),nr1,nk1,nval)
  !  call normcheck(channel(2),nr2,nk2,nval)

  !Due to memory restrictions the calculation of the expected values
  !is nested into the k block reading.
  !All dipole calculations are nested in this subroutine.
  call read_kblocks(channel)

  !adding up the real and complex parts of the dipole observable
  !matrix and multipying by the remaining coefficients.
!  call dipole_add
  
  !output the data and calculate relevant quantities
  !using the dipole transition moments
  call spect


  write(6,*) "program finished "
  
end program dipole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END OF MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_input(channel)
use input
use workdata
implicit none
integer :: channel(2)
  ! setting the default values for the prt control group
  zprint=.false.
    ztra=.true.
  zstart=.false.
  iket=11
  ibra=12
  itra=13
  
  read(5,prt)
  read(5,*) title
  read(5,*) npot,nv2,nv1,ibase2,ibase1
  read(5,*) ezero 

  channel(1)=ibra
  channel(2)=iket
  nth=npot
  
  
  
  ! debug verbose
  if(zprint)then
     write(6,*)'DEBUG INFO:'
     write(6,*)'The input values read through standard input are:'
     
     !line 1
     write(6,*)'-----namelist prt data'
     write(6,*)'zprint is',zprint 
     write(6,*)'ztra   is',ztra   
     write(6,*)'zstart is',zstart 
     write(6,*)'iket   is',iket   
     write(6,*)'ibra   is',ibra   
     write(6,*)'itra   is',itra   
     
     !line 2
     write(6,*)'-----title'
     write(6,*)title
     
     !line 3
     write(6,*)'-----basis data'
     write(6,*)'npot   is',npot   
     write(6,*)'nv1    is',nv2 !swapped to correspond to the definition in the paper    
     write(6,*)'nv2    is',nv1 !swapped to correspond to the definition in the paper       
     write(6,*)'ibase1 is',ibase2 !swapped to correspond to the definition in the paper     
     write(6,*)'ibase2 is',ibase1 !swapped to correspond to the definition in the paper    
     
     !line 4
     write(6,*)'-----ground state'
     write(6,*)'ezero is', ezero
  end if
end subroutine read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_lists(channel)
! complicated stuff starts here
use lists
implicit none
integer :: channel(2),i

  ! initiate the general linked list
  allocate(headlist);nullify(headlist%next);headfirst => headlist
  ! read both file's eaders
  do i=1,2
     write(6,*)'reading header in channel', channel(i)
     call read_head(channel(i))
     allocate(headlist%next) ; headlist => headlist%next ; nullify(headlist%next) 
  enddo
  !tranfer the data from the lists into work variables
  call list_transfer 

end subroutine read_lists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine input_verify
use input
use workdata
implicit none

  if(nval1.lt.nv1+ibase1.or.nval2.lt.nv2+ibase2)then
     write(6,*)'nval1.lt.nv2+ibase2.or.nval2.lt.nv1+ibase1 --- the program will stop'
     stop
  end if
  if(kpar1.eq.kpar2) then
     write(6,*)'kpar1.eq.kpar2 -> there are not allowed transitions (d_z=0)!'
     write(6,*)'Program stops here'
     stop
  end if
  if(ifpar1.ne.ifpar2) then
     write(6,*)'ifpar1.ne.ifpar2 -> there are not allowed transitions!'
!     write(6,*)'Program stops here'
!     stop
  end if
  if(nv1.le.0)then
     nv1=nval1-ibase1
  end if
  if(nv2.le.0)then
     nv2=nval2-ibase2
  end if

end subroutine input_verify


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_head(channel)
!here is the complicated stuff
  use input
  use lists
  implicit none
  integer :: channel
  integer :: i,j,mi    ! loop index
  real(8) :: mdp !meaningless double precision value
  real(8),allocatable :: rgrid(:), eigs(:) 
  
!  open(unit=channel,form='unformatted')

  read(channel)mi,mi,mi,headlist%nr,mi,headlist%jrot,mi,headlist%nval,&
       headlist%kpar,headlist%ifpar
  read(channel)
  read(channel)
  read(channel)
  
!   if(zprint)then 
   write(6,*)'-----header in channel',channel
   write(6,*)"Angular momentum quantum number, J (jrot)=", headlist%jrot
   write(6,*)"Number of radial points, nr=",headlist%nr
   write(6,*)"Number of eigenvalues to be obtained, nval=",headlist%nval
   write(6,*)"parity of the K blocks(0 even, 1 odd), kpar=",headlist%kpar
   write(6,*)"parity of the vibrational wave functions, (f=s+q)=", headlist%ifpar
!   endif

  call readrgrid(channel) ! read and store the r grid points in memory 
 
  read(channel)headlist%nval

  call readeigs(channel) ! read and store the energy eigenvalues in memory
  
!  write(6,*)"determining the number of k blocks"

  if((headlist%kpar.eq.1.and.mod(headlist%jrot,2).eq.0).or.(headlist%kpar.eq.0.and.mod(headlist%jrot,2).eq.1)) then
     headlist%nk=headlist%jrot
  elseif((headlist%kpar.eq.1.and.mod(headlist%jrot,2).eq.1).or.(headlist%kpar.eq.0.and.mod(headlist%jrot,2).eq.0)) then
     headlist%nk=headlist%jrot+1
  endif

!  write(6,*)"the number of k blocks is, nk=", headlist%nk

end subroutine read_head

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_kblocks(channel)
! this is where the actual calculation takes place
  use input
  use lists
  use workdata
  implicit none
  integer,intent(in) ::  channel(2)
  integer :: mi ! meaningless integer
  integer :: i  ! 
  integer :: j,k,l,m,n,p !loop indexes
  integer :: mrot1,mrot2 !angular momentum projection on the lab coordinates
  real(8) :: sa,sb
  real(8),external :: s_dipole

     write(6,*)'---READING THE K BLOCKS IN THE BRA AND KET FILES'
     write(6,*)'---PERFORMING THE CALCULATIONS'
     write(6,*)'K_bra  K_ket  S,niu=-1   S,niu=1'  

  ! loop over the BRA k blocks

  do j=1,nk1
     
     read(channel(1))kz1,ispar1,mi,nth1,mi
     n3d1=(nr1*(nr1+1-2*ispar1)/2)*nth1
     iqpar1=abs(ifpar1-ispar1)
     
     if(zprint)then
     write(6,*)'-----the channel',channel(1),' kblock header is:'
     write(6,*)'kz1',kz1
     write(6,*)'ispar1',ispar1
     write(6,*)'nth1',nth1
     write(6,*)'n3d1 is',n3d1
     end if



     allocate(waves1(n3d1,nv1)) !allocation of bra waves 
     call skipblock(ibase1,channel(1))           !!!        
     do l=1,nv1                                  !!!        
        read(channel(1))(waves1(m,l),m=1,n3d1)   !!!        
        !write(6,*)'waves1 for state',l,'is'     !!! read   
        !write(6,*)waves1(:,l)                   !!! kblock 
        write(15)(waves1(m,l),m=1,n3d1)
     end do                                      !!! waves  
     call skipblock(nval1-ibase1-nv1,channel(1)) !!!        
     ! read the angular grid points                                                
     allocate(thgrid1(nth1))
     read(channel(1))(thgrid1(l),l=1,nth1)


     i=0 !counter for the number of reads on file2
     ! loop over the KET kblocks
     do k=1,nk2 
        
            read(channel(2))kz2,ispar2,mi,nth2,mi ; i=i+1
            n3d2=(nr2*(nr2+1-2*ispar2)/2)*nth2
            iqpar2=abs(ifpar2-ispar2)

            if(zprint)then
            write(6,*)'the file 2 kblock part  header is:'
            write(6,*)'kz2',kz2
            write(6,*)'ispar2',ispar2
            write(6,*)'nth2',nth2
            write(6,*)'n3d2 is',n3d2
            end if

               if(nth1.ne.nth2)then
                  write(6,*)nth1,nth2
                  write(6,*)'---the number of angular grid points is not the same in both files'
                  write(6,*)'---the program will stop'
                  stop
               end if

                     call commongrid ! defines the common angular grid
                     call polin      ! calculates the polinimials at grid


                              
               ! calculate part of the the coefficient
               ! and use it as the selection criterion to
               ! read the k blocks in the KET file
                 
               sa=s_dipole(jrot1,jrot2,iqpar1,iqpar2,kz1,kz2,-1) ! niu = -1
               sb=s_dipole(jrot1,jrot2,iqpar1,iqpar2,kz1,kz2,1)  ! niu = 1

               write(6,*)jrot1,jrot2,iqpar1,iqpar2,kz1,kz2
               write(6,110)kz1,kz2,iqpar1,iqpar2,sa

               !if(.true.) then ! FOR NORM CHECK
               
               if(sa.ne.0.0d0.or.sb.ne.0.d0) then
!               if(sa.ne.100000.0d0.or.sb.ne.10000.0d0) then

                  if(.not.allocated(waves2)) then

                     ! read the k block on file 2
                     !write(6,*)'allocating waves'

                     
                     allocate(waves2(n3d2,nv2)) !allocation of ket waves
                     call skipblock(ibase2,channel(2))          !!!
                     do l=1,nv2                                 !!!
                        read(channel(2))(waves2(m,l),m=1,n3d2)  !!!
                        !write(6,*)'waves1 for state',l,'is'    !!! read 
                        !write(6,*)waves1(:,l)                  !!! kblock
                     end do                                     !!! waves
                     call skipblock(nval2-ibase2-nv2,channel(2))!!!
                     ! read the angular grid points                                                
                     allocate(thgrid2(nth2))                              
                     read(channel(2))(thgrid2(m),m=1,nth2)
                     i=i+nval2+1

                                   
                  endif   ! if(.not.allocated(waves2)) then


                call dipole_int(sa,sb) ! calculated the dipole matrix elements
                                    ! real and complex part are separate

               else
!                  write(6,*)'skipping eigenvalues calculation'
               endif
           
               

                    deallocate(pol2,pol1)

            if(.not.allocated(waves2)) then
!               write(6,*)'k block',k,'in file 2 not elected - skipping'
               call skipblock(nval+1,channel(2)); i=i+nval2+1 !skip all records until the next k block
            elseif(allocated(waves2)) then
!               write(6,*)'deallocated waves2 and thgrid2'

               ! deallocate workdata that is no longer necessary
               ! this is done to re-use workspace
               deallocate(waves2,thgrid2)
              
               deallocate(thgrid)
           

! keeping the matrices allocated during the loops
! does the integration over K blocks                 
!               deallocate(dpba)
!               deallocate(dpbb)

            endif

         end do ! nk2
         
         deallocate(waves1,thgrid1)
         if(allocated(pol1)) deallocate(pol1)

         !rewind up to the first k block on file 2 
         call backblock(i,channel(2)) 
         
      end do ! nk1  


100 format(I2,2X,I2,2X,D25.18,2X,D25.18)
110 format(' sdip ',I2,2X,I2,2X,I2,2X,I2,2X,D25.18)
!200   format('( ',I3,',',I3,') <---> ( ',I3,',',I3,')  th ',I3,'  ',f10.4)


    end subroutine read_kblocks
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

subroutine commongrid
use input
use workdata
implicit none
integer :: i

!print*,'nth is',nth

if(allocated(thgrid))deallocate(thgrid)
if(allocated(wt))deallocate(wt)


allocate(thgrid(nth))
allocate(wt(nth))
!print*,'size of thgrid is',size(thgrid)
!print*,'size of wt is',size(wt)


x1=sqrt((jrot1*(jrot1+1)-kz1**2)/2.0d0 )
!print*,'x1 is',x1
x2=sqrt((jrot2*(jrot2+1)-kz2**2)/2.0d0 )
!print*,'x2 is',x2
x=(x1+x2)/2

call gaujac(thgrid,wt,nth,x,x)

if(zprint)then
write(6,*)'-----common angular grid points (gamma)'
write(6,*)'-----common angular grid integration weights' 
do i=1,nth
write(6,109)i,thgrid(i),wt(i)
end do
109 format(I4,2(2x,d20.13))
end if

!stop
end subroutine commongrid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipcpb(r1,r2,t,dx,dy)

  implicit none
  integer :: i,j,k
  real(8) :: r1,r2,t,pi       ! Radau coordinates
  real(8) :: da,db,dx,dy      ! on plane dipole components

! Internal variables:

  real(8) :: beta,alpha,ang3,R3,sr3,dipcx,dipcy,g,tgb

  g = (3.d0 - sqrt(3.d0)) / (3.d0 + sqrt(3.d0))
  pi=4.d0*atan(1.d0)


  ! consider jacobi r and R  opposite to the radau $\theta$
  ! auxiliary to calculate the Jacobi coordinates
             call geom(r1,r2,t,R3,sr3,ang3)
             call DIPD0(DIPCX,DIPCY,sr3,R3,ang3)

             alpha=acos(t)
             tgb=(r1-g*r2)*tan(alpha*0.5d0)/(r1+g*r2)
             beta = - atan(tgb)
!             beta=beta-pi*0.5d0

             dipcx=dipcx
             dipcy=-dipcy

             dx = DIPCY * cos(beta) + DIPCX * sin(beta) 
             dy =-DIPCY * sin(beta) + DIPCX * cos(beta) 

end subroutine dipcpb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipc(i,j)
  use input
  use workdata
  implicit none
  integer :: i,j,k

  real(8) :: dx,dy


!  real(8) :: g,dipcr1,dipcr2,dipcr1t,dipcr2t,tg,tb,ta2,beta,a,b,g1
!  real(8) :: gamma,pbxp,pbzp,pbx,pby,pbxc,alpha,cf,eps,pbyc,tgc,tgd
!  real(8) :: a1x,a1y,a2x,a2y,a3x,a3y,px,py,rx,ry
!  real(8) :: a1jx,a1jy,a2jx,a2jy,a3jx,a3jy,pjx,pjy,rjx,rjy,theta
!  real(8) :: sd1,sdj1,sd2,sdj2,sd3,sdj3


!  real(8),external :: R3_op 
!  real(8),external :: sr1_op
!  real(8),external :: sr3_op
  
!write(6,*)'nr1 is',nr1
!write(6,*)'nr2 is',nr2

if(.not.allocated(dipcxa)) then
  allocate(dipcxa(nth))
  allocate(dipcya(nth))
  allocate(dipcxb(nth))
  allocate(dipcyb(nth))
endif
  
  dipcxa=0.0d0
  dipcya=0.0d0
  dipcxb=0.0d0
  dipcyb=0.0d0

 ! write(6,*)'executing subroutine dipc'

          do k=1,nth
!
!   Here i introduce the value of the dipole at the grid points.
!
             call dipcpb(rgrid(i),rgrid(j),thgrid(k),dx,dy)
             dipcxa(k) = dx
             dipcya(k) = dy
             call dipcpb(rgrid(j),rgrid(i),thgrid(k),dx,dy)
             dipcxb(k) = dx
             dipcyb(k) = dy

          enddo

end subroutine dipc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine polin
  use input
  use workdata
  implicit none
  integer :: i,j,count


  !write(6,*)'entered polin subroutine'
  if(.not.allocated(pol1))then

     allocate(pol1(nth1,nth))
     call jac_basis(nth,nth1-1,x1,x1,thgrid,pol1)

  endif
  
  allocate(pol2(nth2,nth))
  call jac_basis(nth,nth2-1,x2,x2,thgrid,pol2)

  if(zprint)then
     count=0
     !  do i=1,nth1
     do j=1,nth
        count=count+1
        write(10,*)count,pol1(nth1,j)!,pol2(i,j)
     enddo
     !  enddo
  end if
  
end subroutine polin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipole_int(sa,sb)
  use input
  use workdata
!  
  implicit none
  integer :: i,j,k,l,m1,m2,n,pc1,pc2,aux1,aux2,n2r1,n2r2,j1,j2,m,index
  integer :: ispar, isparloop ! common parity
  real(8) :: fac,diag,sumx,sumy,sa,sb,termx,termy,sx,sy,c,sxd,syd
  real(8) :: waves1d(nth,nv1), waves2d(nth,nv2),sum,wf2
  complex(8) :: ci,cq1,cq2,cq,cox,coy

ci= (0.d0,1.d0)
cq1 = ci**(-iqpar1)
cq2 = ci**(-iqpar2)
cq=conjg(cq1)*cq2


  if(ispar1.eq.1.or.ispar2.eq.1) then
       isparloop=1
  else
       isparloop=0
  endif

  if((ispar1+ispar2).eq.1) then
       ispar=1
  else
       ispar=0
  endif

if(zprint)then
write(6,*)'-----ENTERED DIPOLE_INT SUBROUTINE'
write(6,*)'size(waves1)',size(waves1)
write(6,*)'size(waves2)',size(waves2)
write(6,*)'kz1',kz1
write(6,*)'kz2',kz2
write(6,*)'n3d1',n3d1
write(6,*)'n3d2',n3d2
write(6,*)'n3d1*nval1',n3d1*nval1
write(6,*)'n3d2*nval2',n3d2*nval2
write(6,*)'ispar1',ispar1
write(6,*)'ispar2',ispar2
write(6,*)'isparloop',isparloop
write(6,*)'ispar',ispar
write(6,*)'iqpar1',iqpar1
write(6,*)'iqpar2',iqpar2
write(6,*)'nr1',nr1
write(6,*)'nr2',nr2
write(6,*)'nth1',nth1
write(6,*)'nth2',nth2
write(6,*)'ntheta',nth
write(6,*)'sa',sa
end if

   if(.not.allocated(dpba)) then
      allocate(dpba(nval1,nval2))
      dpba=0.0d0
   endif

sy=((-1.d0)**(iqpar1+iqpar2+1)+1.d0)*sa*((-1.d0)**iqpar1)
sx=-((-1.d0)**(iqpar1+iqpar2)+1.d0)*sa
syd=((-1.d0)**(iqpar1+iqpar2+1)+1.d0)*((-1.d0)**iqpar1)
sxd=-((-1.d0)**(iqpar1+iqpar2)+1.d0)
 
pc1=1
pc2=1

n2r1=nth1*nr*(nr+1-2*ispar1)/2
n2r2=nth2*nr*(nr+1-2*ispar2)/2

do n=1,nr
do l=1,n-isparloop
   call dipc(n,l)       ! calculates the dipole surface at grid

   if(n.eq.l) then
      fac = 1.0d0
      diag = 0.0d0
   else
      fac=0.5d0
      if(ispar.eq.1) then
         diag = - 1.0d0
      else
         diag = 1.0d0
      endif
   endif

      waves1d=0.d0
      do j1=1,nv1
       do m=1,nth
          sum=0.d0
          do k=1,nth1
             index=(pc1-1)*nth1+k
             sum=sum+waves1(index,j1)*pol1(k,m)
          end do
          waves1d(m,j1)=sum
       end do
    end do

      waves2d=0.d0
      do j2=1,nv2
       do m=1,nth
          sum=0.d0
          do k=1,nth2
             index=(pc2-1)*nth2+k
             sum=sum+waves2(index,j2)*pol2(k,m)
          end do
          waves2d(m,j2)=sum
       end do
    end do

do j1=1,nv1
do j2=1,nv2

sumx=0.0d0
sumy=0.0d0

do k=1,nth !integration over the common grid points
termx=(dipcxa(k)+diag*dipcxb(k))*fac
termy=(dipcya(k)+diag*dipcyb(k))*fac
wf2=waves1d(k,j1)*waves2d(k,j2)
sumx = sumx + wf2*wt(k)*termx
sumy = sumy + wf2*wt(k)*termy

enddo ! integration over common grid points

!co = (sa*( sumx + ci * sumy) + sb*( sumx - ci * sumy))*cq

cox = sumx * ( sa + sb ) * cq
coy = sumy * ( sa - sb ) * cq * ci 

!c = (sumx*sx + sumy*sy)
dpba(j1,j2)= dpba(j1,j2) + aimag(cox+coy)

enddo!cycle over state j
enddo!cycle over state i

deallocate(dipcxa,dipcya,dipcxb,dipcyb)        
pc1=pc1+1
pc2=pc2+1
enddo !cycle over columns r_grid
    if(ispar1.eq.0.and.isparloop.ne.0) pc1=pc1+1
    if(ispar2.eq.0.and.isparloop.ne.0) pc2=pc2+1
enddo!cycle over lines   r_grid

152 format(2f20.10)   
153 format(2f20.2)   
154 format(2I2,4f8.3)   
155 format(2I2,1f20.10)   
156 format(2I4,2f20.10)   
157 format(6I2,2f20.10)   

end subroutine dipole_int


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dipole_add
  use input
  use workdata
  implicit none
  integer :: i,j,k,m1,m2,tau
  real(8) :: sum1
  integer,parameter :: nbin=30
  real(8) :: binom(nbin,nbin)
  real(8),external :: threej
  real(8) :: d1,d2      ! debugging
   
  if(.not.allocated(dpba))then
     write(6,*)'no transitions are allowed! the program will stop.'
     stop
  end if

  call setfac(binom,nbin)
  
  sum1 = 0.0d0
  do m1=-jrot1,jrot1
     do m2=-jrot2,jrot2
        do tau=-1,1
           sum1=sum1+&
                (threej(jrot1,1,jrot2,m1,-tau,-m2,binom,nbin))**2
        end do
     end do
  end do

if(zprint) then
write(6,*)'the final dipole result is'
write(6,*)'state1   state2   E1-E2  dipole'

  do i=1,nv1
     do j=1,nv2

        dpba(i,j)=dpba(i,j)**2!+dpbb(i,j)**2  ! sum stored in dx to save space      
                
        write(6,*)i,j,eigs1(i)-eigs2(j),dpba(i,j)        
     end do
  end do
  write(6,*)'sum1 in dipole_add is',sum1
endif

  dpba=dpba*sum1
  
end subroutine dipole_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function R3_op(r1,r2,xcos)
  ! radau input
!  use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: R3_op
  real(8) :: q1,q2,q3

  
  call r_to_q(r1,r2,xcos,q1,q2,q3)
  
  R3_op= SQRT((Q1**2+Q2**2)*0.5d0-Q3**2/4.0d0) 
   
end function R3_op

function R2_op(r1,r2,xcos)
  ! radau input
!  use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: R2_op
  real(8) :: q1,q2,q3
  
  call r_to_q(r1,r2,xcos,q1,q2,q3)
  
  R2_op= SQRT((Q3**2+Q1**2)*0.5d0-Q2**2/4.0d0) 
   
end function R2_op

function R1_op(r1,r2,xcos)
  ! radau input
!  use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: R1_op
  real(8) :: q1,q2,q3
  
  call r_to_q(r1,r2,xcos,q1,q2,q3)
  
  R1_op = SQRT((Q2**2+Q3**2)*0.5d0-Q1**2/4.0d0) 
     
end function R1_op


function sr1_op(r1,r2,xcos)
  ! radau input
!  use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: sr1_op
  real(8) :: q1,q2,q3
  
  call r_to_q(r1,r2,xcos,q1,q2,q3)

  sr1_op=q1
   
end function sr1_op

function sr2_op(r1,r2,xcos)
  ! radau input
!  use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: sr2_op
  real(8) :: q1,q2,q3
  
  call r_to_q(r1,r2,xcos,q1,q2,q3)

  sr2_op=q2
   
end function sr2_op

function sr3_op(r1,r2,xcos)
  ! radau input
  ! use input
  implicit none
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: sr3_op
  real(8) :: q1,q2,q3

  
  call r_to_q(r1,r2,xcos,q1,q2,q3)

  sr3_op=q3
   
end function sr3_op
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine geom(r1,r2,xcos,R3_op,sr3_op,ang3)
implicit none
! radau input
!  use input
  real(8),intent(in) :: r1,r2,xcos
  real(8) :: q1,q2,q3
  real(8) :: R3_op
  real(8) :: sr1_op
  real(8) :: sr3_op
  real(8) :: ang3
  
  call r_to_q(r1,r2,xcos,q1,q2,q3)
  
  R3_op= SQRT((Q3**2+Q2**2)*0.5d0-Q1**2*0.25d0) 
  sr3_op=q1
  ! ! Jacobi angle - homonuclear diatom!! 
  ang3 = - (q3**2-R3_op**2-q1**2*0.25d0)/&
                    (sr3_op*R3_op)
end subroutine geom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine r_to_q(r1r,r2r,xcosr,q1,q2,q3)
  implicit none
  
  ! h3p only

  ! coords
  real(8) :: r1r,r2r,xcosr,r1j,r2j,xcosj
  ! bondlengths
  real(8) :: q1,q2,q3
  ! intermediaries
  real(8) :: f1,f2,f12,p1,p2,s1,s2
  
  real(8) :: g1r,g2r,a,b
  real(8) :: x1=1.0d0,x2=2.0d0
 
  A = SQRT(1.d0 / 3.d0)
  B = 1.d0 / 2.d0
  G1r = X1 - A / (A+B-A*B)
  G2r = X1 - A / (X1-B+A*B)
  ! work out bond lengths
  f1=1./g1r
  f2=1./g2r
  F12= X1 - F1*F2
  P1= R1r*(X1-F1)/(G2r*F12)
  P2= R2r*(X1-F2)/(G1r*F12)
  S1= R1r-P1
  S2= R2r-P2
  Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOSr)/(X1-G1r)
  Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOSr)/(X1-G2r)
  Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOSr)    
  
end subroutine r_to_q



!!!!!!!!!!

!!!!!!!!!!!


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function s_dipole(jrota,jrotb,iqpara,iqparb,kza,kzb,niu)
implicit none
integer,parameter :: nbin=30
real(8) :: binom(nbin,nbin),s_dipole,sfact,sd1,sd2,sd3,sd4,sb
integer :: jrota,jrotb,iqpara,iqparb,kza,kzb,niu,kcorr
real(8),external :: threej
real(8),external :: delta

call setfac(binom,nbin)

!kcorr=(kza+kzb)/2
sfact=0.25d0*sqrt(2.d0)/sqrt((1.0d0+delta(kza,0))*(1.0d0+delta(kzb,0)))
sd1=(-1)**(kzb)*  threej(jrota,1,jrotb,-kza,-niu,kzb,binom,nbin) 
sd2=(-1)**iqparb *threej(jrota,1,jrotb,-kza,-niu,-kzb,binom,nbin)
sd3=(-1)**(kza+kzb+iqpara)*threej(jrota,1,jrotb,kza,-niu,kzb,binom,nbin) 
sd4=(-1)**(iqpara+iqparb+kza)*threej(jrota,1,jrotb,kza,-niu,-kzb,binom,nbin) 
!s_dipole=((-1)**(niu+kcorr))*(sd1+sd2+sd3+sd4)*sfact
s_dipole=((-1)**(niu))*(sd1+sd2+sd3+sd4)*sfact

end function s_dipole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function delta(a,b)
implicit none
integer :: a,b

if(a.eq.b)then
delta=1.d0
else
delta=0.d0
endif
end function delta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#######################################################################
      subroutine skipblock(nskip,ivpb)
                                      
      implicit double precision(a-h,o-y), logical(z)
      do i=1,nskip
         read(ivpb)
         end do
      return

      end subroutine skipblock
!#######################################################################

!#######################################################################
      subroutine backblock(nskip,ivpb)

      implicit double precision(a-h,o-y), logical(z)
      do i=1,nskip
         backspace(ivpb)
         end do
      return

      end subroutine backblock
!#######################################################################

!######################################################################





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine list_transfer

  use lists
  use workdata
  implicit none
  integer :: j
!!!

  headlist => headfirst
  
  nr1=headlist%nr

  jrot1=headlist%jrot

  nval1=headlist%nval

  kpar1=headlist%kpar

  ifpar1=headlist%ifpar

  
  ! transfer data in list rgrid list to rgrid array
  allocate(rgrid1(headlist%nr))
  headlist%rgridlist => headlist%rgridfirst
  do j=1,headlist%nr
     rgrid1(j)=headlist%rgridlist%rgrid
     headlist%rgridlist => headlist%rgridlist%next 
  end do
  write(6,*)'rgrid1 is:'
  write(6,*)rgrid1

  ! transfer data in list eigs list to eigs1 array
  allocate(eigs1(headlist%nval))
  headlist%eigslist => headlist%eigsfirst
  do j=1,headlist%nval
     eigs1(j)=headlist%eigslist%eigs
     headlist%eigslist => headlist%eigslist%next 
  end do





  nk1 = headlist%nk
    write(6,*)"nk1 is",nk1

!!! 

  headlist => headlist%next

  nr2=headlist%nr

  jrot2=headlist%jrot

  nval2=headlist%nval

  kpar2=headlist%kpar

  ifpar2=headlist%ifpar

  ! transfer data in list rgrid list to rgrid array
  allocate(rgrid2(headlist%nr))
  headlist%rgridlist => headlist%rgridfirst
  do j=1,headlist%nr
     rgrid2(j)=headlist%rgridlist%rgrid
     headlist%rgridlist => headlist%rgridlist%next 
  end do
  !write(6,*)'rgrid2 is:'
  !write(6,*)rgrid2

  ! transfer data in list eigs list to eigs2 array
  allocate(eigs2(headlist%nval))
  headlist%eigslist => headlist%eigsfirst
  do j=1,headlist%nval
     eigs2(j)=headlist%eigslist%eigs
     headlist%eigslist => headlist%eigslist%next 
  end do


  
  nk2 = headlist%nk
    write(6,*)"nk2 is",nk2



!!!

  if(nval1.ne.nval2)then
     nval=min(nval1,nval2)
     write(6,*)'---the number of wavefunctions on each file is different.'
     write(6,*)'---calculations will use the lowest number:',nval
  else
     nval=nval1
!     write(6,*)'nval is',nval
  end if


  if(any((rgrid1.ne.rgrid2).or.(nr1.ne.nr2)))then
     write(6,*)'---the rgrid points are different.'
     write(6,*)'---the program will stop'
     stop
  else
     allocate(rgrid(headlist%nr))
     rgrid=rgrid1
     nr=nr1
  endif

end subroutine list_transfer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine readrgrid(channel)
  use lists
  implicit none
  integer :: channel 
  integer :: i,j
  real(8),allocatable :: rgrid(:) 

  allocate(rgrid(headlist%nr))
  read(channel)(rgrid(j),j=1,headlist%nr)

  allocate(headlist%rgridlist)  ;  nullify(headlist%rgridlist%next)
  headlist%rgridfirst => headlist%rgridlist
  do j=1,headlist%nr
     headlist%rgridlist%rgrid=rgrid(j)
     allocate(headlist%rgridlist%next) ;  headlist%rgridlist => headlist%rgridlist%next ;  nullify(headlist%rgridlist%next)
  enddo

  !write(6,*)"rgrid values are"  

  headlist%rgridlist => headlist%rgridfirst
  do j=1,headlist%nr
  !write(6,*)headlist%rgridlist%rgrid
  headlist%rgridlist => headlist%rgridlist%next
  enddo

end subroutine readrgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readeigs(channel)
  use lists
  implicit none

  integer :: j,channel
  real(8),allocatable :: eigs(:) 


  allocate(eigs(headlist%nval))
  read(channel)(eigs(j),j=1,headlist%nval)

  allocate(headlist%eigslist)  ;  nullify(headlist%eigslist%next)
  headlist%eigsfirst => headlist%eigslist
  do j=1,headlist%nval
     headlist%eigslist%eigs=eigs(j)
     allocate(headlist%eigslist%next) ;  headlist%eigslist => headlist%eigslist%next ;  nullify(headlist%eigslist%next)
  enddo


!  write(6,*)"eigs values are"  

  headlist%eigslist => headlist%eigsfirst
  do j=1,headlist%nval
!  write(6,*)headlist%eigslist%eigs
  headlist%eigslist => headlist%eigslist%next
  enddo

end subroutine readeigs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!######################################################################

! calls 
!      call gaujac(xt,wt,n,aa,aa)

!      call jac_basis(n,npol,aa1,aa1,xt,pol1)
!      call jac_basis(n,npol,aa2,aa2,xt,pol2)



!##################################################################

      MODULE constants
      IMPLICIT NONE
      INTEGER, PARAMETER :: real_kind=SELECTED_REAL_KIND(8,40)
      END MODULE constants

!##################################################################

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
!            if (i.eq.0) write(33,*)n,basis(i,n),norm(i)
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
      REAL(KIND=real_kind), ALLOCATABLE :: A1n(:),A2n(:),A3n(:),A4n(:)
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

!###################################################################

subroutine normcheck(fi,nr,nk,nval)

  implicit none
  
  integer :: i,j,k,l,m     ! index energy levels
  integer :: mi    ! meaningless integer
  integer :: fi  !file numerator
  integer :: n     !file line counter
  integer :: nr    !number of r grig points
  integer :: nk    !number of k blocks
  integer :: nval !number of states
  integer :: kz,ispar,nth,n3d 

  real(8), allocatable :: waves(:,:), norm(:,:)


     allocate(norm(nval,nval))
     norm=0.0d0
  !write(6,*)'checking the norm'
  n=0
  do i=1,nk
     read(fi)kz,ispar,mi,nth,mi ; n=n+1

     n3d=(nr*(nr+1-2*ispar)/2)*nth

   ! write(6,*)'the file 1 kblock part  header is:'
   ! write(6,*)'kz',kz
   ! write(6,*)'ispar',ispar
   ! write(6,*)'nth',nth
   ! write(6,*)'n3d is',n3d

     allocate(waves(n3d,nval))
     
     do j=1,nval
        read(fi)(waves(l,j),l=1,n3d); n=n+1
        !write(6,*)'waves for state',l,'is'
        !write(6,*)waves(:,l)
     end do

     ! iteratively calculate the norm of the wavefunctions in file 1
     ! the result is printed out after getting the 

     read(fi); n=n+1
     
!     do k=1,nr
!        do l=1,k-ispar
!           do m=1,nth
              
     do l=1,n3d
        do j=1,nval
           do k=1,nval
              norm(j,k) = norm(j,k)+(waves(l,j)*waves(l,k))
           enddo
        enddo
     enddo
     deallocate(waves)
     
  enddo
  
!write(6,*)'The file number',fi,'has norms' 
  
!  do i=1,nval
!     do j=1,nval
!        write(6,*)'norm',i,j,'is',norm(i,j)
!     enddo
!  enddo

  call backblock(n,fi)
  deallocate(norm)


end subroutine normcheck




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **025
subroutine spect
  
  !     subroutine spect calculates s(f-i), the line strength,
  !     and controls the printing out of the transition moments,
  !     frequencies and nuclear spin effects.
  !
  !     if ztra is .true. it calls outpt to output energies
  !     and line strengths for program spectrum to calculate
  !     simulated spectra
  
  ! based on Jayesh's dipole program
  ! adapted by bruno on 25/07/2005  
  ! pb also worked on it
  

  use input
  use workdata
  implicit none
  integer :: nn2,ie1,ie2
  real(8) :: autocm,autode,detosec,x0 
  real(8) :: xf,dd,dd3,sx,sxd,a,cc
  real(8), dimension(nval1,nval2) :: sint
  real(8) :: xe1, xe2
!  real(8), dimension(nval2) :: xe2
!  character(len=8)  title(9)
  
  data autocm/2.19474624d+05/,x0/0.0d0/,&  
       autode/2.5417662d0/,&
       detosec/3.136186d-07/
  !     autocm converts atomic units to wavenumbers
  !     autode converts atomic units to debye
  !     detose! converts from s(f-i) in debye**2 to seconds
  
if(zprint) then 
write(*,*) '-----SUBROUTINE SPECT'
write(*,*) '-----passed data'
write(*,*) 'nval1 is',nval1
write(*,*) 'nval2 is',nval2
end if


  if(.not.allocated(dpba))then
     write(6,*)'no transitions are allowed! the program will stop.'
     stop
  end if



  write(6,200)
200 format(///)
  write(6,201)
201   format(//,5x,'*************************************************'&
            //,5x,'print out of dipole transition moments and s(f-1)'&
            //,9x,'frequencies in wavenumbers',&
            /,9x,'transition moments in debye (2.54174a.u.)',&
            /,9x,'s(f-i) in debye**2',& ! multiply the "dipole" by autode**2
            /,9x,'einstein a-coefficient in sec-1',//,&
            5x,'*************************************************')
  
  write(6,200)
  write(6,202) title
202 format(5x,9a8)
  write(6,200)
  write(6,203) jrot1, jrot2, ifpar1, ifpar2
203 format(5x,'  jrot1=',i4,'  jrot2=',i4,&
         '  ifpar1=',i4,'  ifpar2=',i4,//)
!  ezero=x0
!  read(5,505,end=555) ezero
  write(6,204) ezero
204 format(5x,'ground zero =',e16.8,' cm-1')
  
  write(6,200)
  write(6,205)
205 format(/,' ie1 ie2   ket energy   bra energy    frequency', &
         '       s(f-i)      a-coefficient' ,/)
  
  !     xf is the factor that allows for the root$(2j' + 1)(2j"+1)
  !     left over from the calculation of the transition moment
  
  xf= sqrt(dble((2*jrot1+1)*(2*jrot2+1)))
  
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
  !       write(itra) j1,j2,kmin1,kmin2,nval1,nval2,ifpar1,ifpar2,gz,zembed
  !       gz is ezero  in cm-1
  
          cc=xf*autode

  do 1 ie1=1,nv1
        xe1= eigs1(ie1)*autocm - ezero
     
     do 2 ie2=1,nv2
        xe2= eigs2(ie2)*autocm - ezero

        dd= xe2 - xe1
        dd3= abs(dd*dd*dd)
        
        ! it starts here
        
        sxd= (dpba(ie1,ie2)*cc)**2
        sint(ie1,ie2)= sxd

        if (dd .gt. x0) a= sxd*dd3*detosec/dble(2*jrot2 + 1)
        if (dd .lt. x0) a= sxd*dd3*detosec/dble(2*jrot1 + 1)

        write(6,206) ie1+ibase1,ie2+ibase2,xe1,xe2,dd,sxd,a 
        write(96,209) ie1+ibase1,ie2+ibase2,xe1,xe2,dd,a 
209     format(2(i4),3(3x,f12.5),f16.8)
206     format(2(i4),3(3x,f10.3),(2x,e13.6),(2x,f13.6))
2       continue
!        if (.not.zpmin .or. ie1.le.10) write(6,207)
!        if (ie1.le.nval1) write(6,207)
        write(6,207)
207     format(//)
1       continue
        
        ! writes in itra e1 and e1 : energy values in au
        ! then it will write the values of sint ( in au as weel)
        
        if (ztra) call outpt(eigs1,eigs2,sint,ezero)
        return
        
     end subroutine spect
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !                                                **026
   subroutine outpt(e1,e2,sint,gz)
     
     !     subroutine outpt outputs the data necessary for program
     !     spectrum to simulate laboratory or interstellar spectra.
     !     the output data is in atomic units.
        use input
        use workdata     
        implicit none
        integer :: ie2
!        common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
!             npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,&
!             kmin1,kmin2,jk1,jk2,nval1,nval2,nn2,ibase1,ibase2,ipot
!        common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2
!        common/sym/ idia,ifpar1,ifpar2,jrot1,jrot2
!      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
!      common /mass/ xmass(3),g1,g2,zembed,zbisc
      
      real(8), dimension(nval1) :: e1
      real(8), dimension(nval2) :: e2
      real(8), dimension(nval1,nval2) :: sint
      real(8) :: gz

      
      !     is this the first write-out?
      
!      open(unit=itra,form='unformatted') 
      if (.not.zstart) then
10       read(itra,end=90)
         goto 10
90       continue
         ! *****  inclusion of the following card is machine dependent *****
         !        backspace itra
      endif
      
      write(itra) jrot1,jrot2,nval1,nval2,ifpar1,ifpar2,gz,&
           ibase1,ibase2
      write(itra) e1
      write(itra) e2
      do 20 ie2=1,nval2
      call outrow(sint(1,ie2),nval1,itra)
20    continue
      return
      end subroutine outpt
   !ccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                **028
      subroutine outrow(row,nrow,iunit)

      implicit double precision (a-h,o-y)
      double precision, dimension(nrow) :: row
      write(iunit) row
      return
      end

