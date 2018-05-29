!***********************************************************************!
!***********************************************************************!
!   THE MOST ELEMENTARY CLASS COMES FIRST. THE MOST DEPENDENT CLASS     !
!   COMES LAST. IN THIS FILE, THE LAST IS THE ONE THAT PERFOMS THE      !
!   CONVERSION BETWEEN WAVEFILES AND KBLOCKS.                           !
!***********************************************************************!
!***********************************************************************!

module class_error
implicit none
public
logical :: debug=.true.
contains

  subroutine error(msg)
    implicit none
    character(len=*),intent(in) :: msg
    write(6,*) "(ERROR) ",msg
    stop
  end subroutine error
  
end module class_error


module constants
  implicit none
  save
  ! parameters get compiled in as just numbers
  ! so this module should cost nothing
  INTEGER, PARAMETER :: real_kind=SELECTED_REAL_KIND(8,40)
  !an accuracy tolorence
  double precision,parameter :: toler=1.0d-8

  ! just some numbers >>>
  double precision,parameter :: half=0.5d0,zero=0.0d0,one=1.0d0
  double precision,parameter :: two=2.0d0,four=4.0d0,eight=8.0d0,sixteen=1.6d1
  double precision,parameter :: xp5=0.5d0,x0=0.0d0,x1=1.0d0,x3=3.0d0
  double precision,parameter :: x2=2.0d0,x4=4.0d0,x5=5.0d0,x8=8.0d0,x16=1.6d1
  double precision,parameter :: sqrt2=1.4142135623731d0

  ! wavenumber = autocm * hartree
  double precision,parameter :: autocm=2.19474624d+05

  !block cyclic descriptor array pointers
  integer,parameter :: DESC_LENGTH=9
  integer,parameter :: DESC_DT=1,DESC_CTXT=2,DESC_M=3
  integer,parameter :: DESC_N=4,DESC_MB=5,DESC_NB=6
  integer,parameter :: DESC_RSRC=7,DESC_CSRC=8,DESC_LLD=9

  !default block size for scalapack calculations
  ! recomended values are 32 and 64
  ! block size can be important for scalability
  ! one recomendation for sp2's suggests block size
  ! of 10 to 50 for minimum run time. Work space is minimum
  ! for minimum block size
  integer,parameter :: BLOCKSIZE=64

  !tweakable parameter for rough density of states
  !function that provides loadleveling
  !see rough_equal_desnsity_of_states in pdvr3drz.f90
  double precision,parameter :: roughpower=(2.0d0/7.0d0)

  !numbers of timers in code (see utils.f90)
  integer,parameter :: NTIMERS=7
  !name the timers for readability
  integer,parameter :: TIMEALL=1,TIMECCMAIN=2,TIMEVIBMAIN=3
  integer,parameter :: TIMEDIAG=4,TIMETRANSR=5,TIMEREDIST=6
  integer,parameter :: TIMEBLOC3D=7

  ! file streams
  integer,parameter :: STDERR=0
  integer,parameter :: STDIN=5
  integer,parameter :: STDOUT=6

end module constants

module class_ang_int_pnts
  use class_error
  implicit none
  private
  type, public :: ang_int_pnts
     real(8),pointer :: position(:)
     real(8),pointer :: wheight(:)
  end type ang_int_pnts
  public :: ang_grid ! this is the manual constructor
  !function ang_grid(npot,kz) result(p)
  !  implicit none
  !  type(ang_int_pnts) :: p 
  !  integer :: npot,kz
  !end function ang_grid
  public :: ang_grid_zperp ! this is the manual constructor
  !function ang_grid_zperp(npot,jrot,k) result(p)
  !  implicit none
  !  type(ang_int_pnts) :: p 
  !  integer :: npot,jrot,k
  !end function ang_grid

  public :: ang_grid_destroy ! this is the manual constructor
  !subroutine ang_grid_destroy(p)
  !  type(ang_int_pnts) :: p
  !end subroutine ang_grid_destroy

contains

  subroutine ang_grid_destroy(p)
    type(ang_int_pnts) :: p 
    if(associated(p%position)) deallocate(p%position)
    if(associated(p%wheight)) deallocate(p%wheight)
  end subroutine ang_grid_destroy

  function ang_grid(npot,kz) result(p)
    !     set up points & weights for npot point angular integration
    implicit none
    type(ang_int_pnts) :: p 
    integer :: npot ! number of integration points
     real(8), dimension(npot) :: xd,wtd !grid points position, grid points wheight
    real(8), dimension(npot) :: b,c ! b (recurrence), c(recurrence) 
    real(8) :: xi,x0,x1,x2,x3,x4,toler,csa,tsa,alf,bta,xkz
    integer :: i,j, nn2 ,kz
    data toler/1.0d-8/,x2/2.0d0/

    !STUDY THESE CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
    !      if (jrot > 0) then
    !         ltop = lpot
    !         if (idia .eq. 2 .and. mod(ltop,2) .ne. ipar) ltop=ltop+1
    !         npot=((ltop+2)/2)*2
    !         nn2= npot/2
    !      endif
    
    nn2=npot/2
    allocate(p%position(npot),p%wheight(npot))
    
    !kz=1 
    xkz=dble(kz)
    
    !     tsa is the exact sum of weights for gauss-jacobi integration                                                                                                                                                                                                          
    tsa= x2**(kz+kz+1)/dble(kz+1)                                                                                                                                                                                                                                         
    loop30: do i=1,kz                                                                                                                                                                                                                                                            
       tsa=tsa*dble(i)/dble(kz+i+1)                                                                                                                                                                                                                                  
    end do loop30
   

! paolo's jacobi polinomials
    call jacbasis(xd,wtd,npot,xkz,xkz)
    csa=sum(wtd)
    
!    sum the weights and check 
if(debug) then
!    call jacobi2(npot,nn2,xd,wtd,xkz,xkz,csa,tsa)
    write(6,1000) npot,kz,(xd(i),wtd(i),i=1,nn2)
1000 format(//i8,' point Gauss-associated Legendre integration with',&
          ' k =',i3//5x,'integration points',11x,'weights',&
          //(f23.15,d25.12))
    write(6,1010) csa,tsa
1010 format(/4x,'Computed sum of weights',d22.15,&
          /4x,'Exact    sum of weights',d22.15//)
 end if
    if (abs((csa-tsa)/tsa) .gt. toler) then
       write(6,"(//5x,'Points & weights in error, adjust algorithm'//)")
       stop
    endif
    !     define other integration points
    do i=1,nn2
       j=i+nn2
       xd(j)=-xd(nn2+1-i)
       wtd(j)=wtd(nn2+1-i)
    end do
    do i = 1, 2*nn2
       xd(i) = -xd(i)
    end do
    p%position = xd
    p%wheight  = wtd
  end function ang_grid
  
  function ang_grid_zperp(nang,jrot,k) result(p)
    implicit none
    type(ang_int_pnts) :: p
    integer :: nang,jrot,k
    integer :: nang2
    real(8) :: jjp1, alf
    !setup the angular grid
    ! setup params
    allocate(p%position(nang),p%wheight(nang))
    nang2 = ( nang + 1 ) / 2
    jjp1 = dble( jrot * ( jrot + 1 ) )
    alf = sqrt( ( jjp1 - dble( k**2 ) ) / 2.d0 )
    ! calculate the legendre integration grid points (angular)
    call gaujac( p%position , p%wheight , 2*nang2 , alf , alf )
  end function ang_grid_zperp


subroutine jacbasis(dg,w,n,alf,bet)

  implicit none

  integer :: n
  integer :: info,j,i,di
  real*8 :: xi,alf,bet,a1,a2,a3,a4,ua3,zc,termd,termu,terml
  real*8 :: dg(n),w(n),vnor(-1:n)
  real*8 :: eigen(n)
  real*8,external :: gammaln

real(8), allocatable :: dg1(:),zd(:,:),work(:)
allocate(dg1(n),zd(n,n),work(5*n))

      call norms2(vnor,n,alf,bet)

  do i=0,n-1
     xi=dble(i)
     a1=2.d0*(xi+1.d0)*(xi+alf+bet+1.d0)*(2.d0*xi+alf+bet)*vnor(i)/vnor(i+1)
     a2=-(2.d0*xi+alf+bet+1.d0)*(alf**2-bet**2)*vnor(i)*vnor(i)
     a4=2.d0*(xi+alf)*(xi+bet)*(2.d0*xi+alf+bet+2.d0)*vnor(i)/vnor(i-1)
     a3=(2.d0*xi+alf+bet+2.d0)*(2.d0*xi+alf+bet+1.d0)*(2.d0*xi+alf+bet)
     ua3=1.d0/a3
     termu=a1*ua3
     termd=a2/ua3
     terml=a4*ua3
     dg(i+1)=termd
     if (i.ge.1) dg1(i)=terml
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

       SUBROUTINE norms2(norm,nn,alf,bet)
       USE constants
       IMPLICIT NONE
       INTEGER :: n,nn
       REAL(KIND=real_kind) :: alf, bet,lmd,norm1
       REAL(KIND=real_kind) :: norm(-1:nn)
       REAL(KIND=real_kind) :: a1,a2,a3,a4
       REAL(KIND=real_kind) :: gammln
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


  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ZPERP ROUTINES

SUBROUTINE gaujac(x,w,n,alf,bet)
  USE constants
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=8):: alf,bet
  REAL(KIND=8), intent(out):: w(n),x(n)
  REAL(KIND=8),PARAMETER :: EPS=3.0D-14
  INTEGER,PARAMETER :: MAXIT=10
  INTEGER :: i,its,j
  REAL(KIND=8)::alfbet,an,bn,r1,r2,r3
  REAL(KIND=8)::c1,c2,c3,p1,p2,p3,pp,temp,z,z1
!  REAL(KIND=8),EXTERNAL :: gammln

  DO i=1,n
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

     DO its=1,MAXIT
        temp=x2+alfbet
        p1=(alf-bet+temp*z)/x2
        p2=x1
        DO j=2,n
           p3=p2
           p2=p1
           temp=x2*DBLE(j)+alfbet
           c1=x2*DBLE(j)*(DBLE(j)+alfbet)*(temp-x2)
           c2=(temp-x1)*(alf*alf-bet*bet+temp* &
                &               (temp-x2)*z)
           c3=x2*(DBLE(j-1)+alf)*(DBLE(j-1)+bet)*temp
           p1=(c2*p2-c3*p3)/c1
        enddo
        pp=(DBLE(n)*(alf-bet-temp*z)*p1+x2*(DBLE(n)+alf)* &
             &    (DBLE(n)+bet)*p2)/(temp*(x1-z*z))
        z1=z
        z=z1-p1/pp
        IF(ABS(z-z1).LE.EPS) GOTO 1
     enddo
1    x(i)=z
     w(i)=DEXP(gammln(alf+DBLE(n))+gammln(bet+DBLE(n))    &
          &    -gammln(DBLE(n+1))-gammln(DBLE(n)+alfbet+x1))&
          &      *temp*x2**alfbet/(pp*p2)

  enddo
  RETURN
END SUBROUTINE gaujac


FUNCTION gammln(XX)
  IMPLICIT NONE
  INTEGER :: j
  REAL(KIND=8)::GAMMLN,XX
  REAL(KIND=8)::SER,STP,TMP,X,COF(6)
  REAL(KIND=8)::HALF,ONE,FPF
  DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
       &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
  DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
  X=XX-ONE
  TMP=X+FPF
  TMP=(X+HALF)*LOG(TMP)-TMP
  SER=ONE
  DO J=1,6
     X=X+ONE
     SER=SER+COF(J)/X
     !        SER=SER+COF(J)/(X+0.000001d0)
  enddo
  GAMMLN=TMP+LOG(STP*SER)
  RETURN
END FUNCTION gammln

end module class_ang_int_pnts

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
!***********************************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module class_mass
  implicit none
  save
  real(8) :: g1,g2 ! geometry parameters
  real(8) :: xmass(3)
end module class_mass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module class_geom
  use class_mass
  implicit none
contains
!  subroutine j_to_q(sr,br,xcos,q1,q2,q3,coord)
  subroutine j_to_q(sr,br,xcos,q1,q2,q3)
    implicit none
    ! coords
    real(8),intent(in) :: sr,br,xcos
    real(8) :: sr2,br2,xcos2,brxcos,gmo,gsr,gmosr
    ! bondlengths
    real(8) :: q1,q2,q3
    real(8) :: g
    real(8) :: xxcos


    ! HACK WARNING: the dvr code seems to have the angle on the left side of
    ! the Jacobi coordinate
    xxcos=-xcos
    !
    !
!USING MY DEFINITION OF BONDLENGTHS AS OUTPUT
!    if(coord.eq.'1-23') then
       g=xmass(3)/(xmass(2)+xmass(3))
       br2=br*br
       brxcos=2.0d0*br*xxcos
       gmo=g-1.0d0
       gsr=g*sr
       gmosr=gmo*sr
       q3=sqrt(br2+(gsr-brxcos)*gsr)
       q2=sqrt(br2+(gmosr-brxcos)*gmosr)
       q1=sr 
 !   elseif(coord.eq.'3-12') then
 !      g=xmass(2)/(xmass(1)+xmass(2))
 !      q2=sqrt(br**2+(g*sr)**2-2*g*sr*br*xxcos)
 !      q1=sqrt(br**2-2*br*sr*xxcos*(g-1.0d0)+((g-1.0d0)*sr)**2)
 !      q3=sr 
 !   elseif(coord.eq.'2-31') then
 !      g=xmass(1)/(xmass(3)+xmass(1))
 !      q1=sqrt(br**2+(g*sr)**2-2*g*sr*br*xxcos)
 !      q3=sqrt(br**2-2*br*sr*xxcos*(g-1.0d0)+((g-1.0d0)*sr)**2)
 !      q2=sr 
 !   end if

    !write(6,*)q1,q2,q3

  end subroutine j_to_q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine r_to_q(r1r,r2r,xcosr,q1,q2,q3)
    implicit none !CHECK IF THIS ONLY WORKS FOR HOMONUCLEAR DIATOMIC 12 RADAU
    ! coords
    real(8) :: r1r,r2r,xcosr,r1j,r2j,xcosj
    ! bondlengths
    real(8) :: q1,q2,q3
    ! intermediaries
    real(8) :: f1,f2,f12,p1,p2,s1,s2
    real(8) :: g1r,g2r,a,b
    real(8),parameter :: x1=1.0d0,x2=2.0d0
    A = SQRT(XMASS(3) / (XMASS(1)+XMASS(2)+XMASS(3)))
    B = XMASS(2) / (XMASS(1)+XMASS(2))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine q_to_j(q1,q2,q3,r1j,r2j,xcosj,coord)
    implicit none
    character(len=4),intent(in):: coord !coordinate embedding
    real(8),intent(out) :: r1j,r2j,xcosj!jacobi coordinates, small r , big r and jacobi angle cosine
    real(8),intent(in) :: q1,q2,q3      ! intermediate, bondlength coordinates
    real(8) :: g                        ! diatomic small r mass factor
!intermediate variables:
    real(8) :: q12,q22,q32,g_2

    !these expressions were derived simply using
    !the cosine rule
    !and the bondlength convention defined in SUTC1991:183

    q12=q1*q1
    q22=q2*q2
    q32=q3*q3

    select case (coord)
    case('1-23')
       g=xmass(3)/(xmass(2)+xmass(3))
       g_2=g*g
       r1j = q1
       r2j = SQRT(q32+g_2*q12-g*(q32+q12-q22)) 
       xcosj = (g_2*q12+r2j*r2j-q32)/(2.0d0*r2j*g*q1)
    case('3-12')
       g=xmass(2)/(xmass(1)+xmass(2))
       g_2=g*g
       r1j = q3
       r2j = SQRT(q22+g_2*q32-g*(q22+q32-q12)) 
       xcosj = (g_2*q32+r2j*r2j-q22)/(2.0d0*r2j*g*q3)
    case('2-31')
       g=xmass(1)/(xmass(3)+xmass(1))
       g_2=g*g
       r1j = q2
       r2j = SQRT(q12+g_2*q22-g*(q12+q22-q32)) 
       xcosj = (g_2*q22+r2j*r2j-q12)/(2.0d0*r2j*g*q2)
    case default
       write(6,*)'unrecognised embedding (',coord,') found in subroutine q_to_j' 
       stop
    end select
  end subroutine q_to_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine blba_to_q(r1,r2,xcos,q1,q2,q3)
    real(8),intent(in) :: r1,r2,xcos      !BLBA coordinates
    real(8),intent(out) :: q1,q2,q3       !intermediate, bondlength coordinates
    q1=r1
    q2=r2
    q3=sqrt(r1*r1+r2*r2-2.0d0*r1*r2*xcos)
  end subroutine blba_to_q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine q_to_blba(q1,q2,q3,r1,r2,xcos)
    real(8),intent(out) :: r1,r2,xcos      !BLBA coordinates
    real(8),intent(in) :: q1,q2,q3       !intermediate, bondlength coordinates
    r1=q1
    r2=q2
    xcos=(q1*q1+q2*q2-q3*q3)/(2*q1*q2)
  end subroutine q_to_blba


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine blba_to_j(br1,br2,bxcos,jr1,jr2,jxcos,coord)
    real(8), intent(in) :: br1,br2,bxcos
    character(len=4), intent(in) :: coord
    real(8), intent(out) :: jr1,jr2,jxcos
    real(8) :: q1,q2,q3
    call blba_to_q(br1,br2,bxcos,q1,q2,q3)
    call q_to_j(q1,q2,q3,jr1,jr2,jxcos,coord)
  end subroutine blba_to_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine j_to_blba(jr1,jr2,jxcos,br1,br2,bxcos,coord)
    real(8), intent(out) :: br1,br2,bxcos
    character(len=4), intent(in) :: coord
    real(8), intent(in) :: jr1,jr2,jxcos
    real(8) :: q1,q2,q3
    call j_to_q(jr1,jr2,jxcos,q1,q2,q3)
    call q_to_blba(q1,q2,q3,br1,br2,bxcos)
  end subroutine j_to_blba

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine r_to_j(rr1,rr2,rxcos,jr1,jr2,jxcos,coord)
    real(8), intent(in) :: rr1,rr2,rxcos
    character(len=4), intent(in) :: coord
    real(8), intent(out) :: jr1,jr2,jxcos
    real(8) :: q1,q2,q3
    call r_to_q(rr1,rr2,rxcos,q1,q2,q3)
    call q_to_j(q1,q2,q3,jr1,jr2,jxcos,coord)
  end subroutine r_to_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE RJ(R1r,R2r,XCOSr,Q1,Q2,Q3)
    !Jayesh 8/10/2002
    !This subroutine converts between Radua and Jacobi coordinates
    implicit none
    real(8),intent(in) :: r1r,r2r,xcosr
    real(8) :: r1j,r2j,xcosj
    real(8),parameter :: x0=0.0d0
    real(8),parameter :: x1=1.0d0
    real(8),parameter :: x2=2.0d0
    real(8),parameter :: tiny=9.0D-15

    real(8) :: q1,q2,q3,xx,yy,f1,f2,f12,p1,p2,s1,s2


    !If the cooordninates are already Jacobi then pseudo Radau coordinates are
    !simply copied to the correct Jacobi variables. Q1, Q2, and Q3 are then 
    !calculated for the particular coordinate system
    IF (G1 .EQ. X0) THEN
       !  BONDLENGTH BONDANGLE COORDINATES: ATOM 1 = ATOM 2
       XCOSj = XCOSr
       R1j = R1r
       R2j = R2r

       Q1 = R1j
       Q2 = R2j
       Q3= DSQRT(R1j*R1j + R2j*R2j - X2*R1j*R2j*XCOSj)

    ELSEIF (G2.EQ.0.0) THEN
       !  SCATTERING COORDINATES: ATOM 2 = ATOM 3
       XCOSj = XCOSr
       R1j = R1r
       R2j = R2r

       XX = R1j * G1
       YY = R1j * (X1 - G1)
       IF (R2j .EQ. X0 .OR. XCOSj .GE. (X1 - TINY)) THEN
          Q1 = DABS(XX - R2j)
          Q2 = (YY + R2j)
       ELSE IF (XCOSj .LE. (TINY - X1)) THEN
          Q1 = (XX + R2j)
          Q2 = DABS(YY + R2j)
       ELSE
          Q1 = DSQRT(XX*XX + R2j*R2j - X2*XX*R2j*XCOSj)
          Q2 = DSQRT(YY*YY + R2j*R2j + X2*YY*R2j*XCOSj)
       ENDIF
       Q3 = R1j

       !Else proceed with the conversion to Jacobi
    ELSE
       !     GENERAL COORDINATES (INCLUDING RADAU): ATOM 1 = ATOM 2
       F1= X1/G1
       F2= X1/G2
       F12= X1 - F1*F2
       P1= R1r*(X1-F1)/(G2*F12)
       P2= R2r*(X1-F2)/(G1*F12)
       S1= R1r-P1
       S2= R2r-P2
       Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOSr)/(X1-G1)
       Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOSr)/(X1-G2)
       Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOSr)         

       !     Conversion Q1,Q2,Q3 in Jacobi coordinats                  
       R1j=Q3
       R2j=DSQRT((Q2**2+Q1**2)*0.5d0-Q3**2/4.0d0)      
       XCOSj=((Q2**2-Q1**2)/(R2j*Q3))*0.5d0
    endif

    return
  END SUBROUTINE RJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1!!!!!!!!!!!!!!!!!
end module class_geom



!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
!***********************************************************************!



module class_dvr_grid
  use class_ang_int_pnts
  use class_geom
  use class_mass
  use class_error
  private
  type,public :: dvr_grid
     real(8), pointer :: xmass1 ! vibrational masses of the atoms in atomic units
     real(8), pointer :: xmass2 ! vibrational masses of the atoms in atomic units
     real(8), pointer :: xmass3 ! vibrational masses of the atoms in atomic units
     integer, pointer :: idia
     integer, pointer :: npnt1
     real(8), pointer :: r1(:)
     integer, pointer :: npnt2
     real(8), pointer :: r2(:)
     integer, pointer :: nang
     real(8), pointer :: ang(:)
  end type dvr_grid

  public :: dvr_grid_ang_build
  !function dvr_grid_ang_build(nang,k)result(ang_)
  !  integer :: nang, k
  !  real(8) :: ang_(nang)
  !end function dvr_grid_ang_build
  public :: dvr_grid_ang_build_zperp
  !function dvr_grid_ang_build_zperp(nang,jrot,k)result(ang_)
  !  integer :: nang, jrot,k
  !  real(8) :: ang_(nang)
  public :: dvr_grid_wtd_build
  !function dvr_grid_wtd_build(nang,k)result(wtd_)
  !  integer :: nang, k
  !  real(8) :: wtd_(nang)
  !end function dvr_grid_wtd_build
  public :: dvr_grid_wtd_build_zperp
  !function dvr_grid_wtd_build_zperp(nang,jrot,k)result(wtd_)
  !  integer :: nang, jrot,k
  !  real(8) :: wtd_(nang)
  public :: dvr_grid_write
  !subroutine dvr_grid_write(unit,grid)
  !type(dvr_grid) :: grid
  !integer :: unit
  !end subroutine dvr_grid_write
  public :: dvr_grid_dump
  !subroutine dvr_grid_dump(unit,grid)
  !type(dvr_grid) :: grid
  !integer :: unit
  !end subroutine dvr_grid_dump
  public :: dvr_grid_read
  !subroutine dvr_grid_read(unit,grid)
  !type(dvr_grid) :: grid
  !integer :: unit
  !end subroutine dvr_grid_read
  public :: dvr_grid_destroy
  !subroutine dvr_grid_destroy(grid)
  !type(dvr_grid), pointer :: grid
  !end subroutine dvr_grid_destroy
  public :: dvr_grid_r2_max
  !function dvr_grid_r2_max(grid,coord) result(r2max)
  !  type(dvr_grid) :: grid
  !  character(len=4),intent(in) :: coord 
  !  real(8) :: r2max
  !end function dvr_grid_r2_max

contains

  subroutine dvr_grid_destroy(grid)
    type(dvr_grid), pointer :: grid
    deallocate(grid%xmass1)
    deallocate(grid%xmass2)
    deallocate(grid%xmass3)
    deallocate(grid%idia  )
    deallocate(grid%npnt1 )
    deallocate(grid%r1    )
    deallocate(grid%npnt2 )
    deallocate(grid%r2    )
    deallocate(grid%nang  )
    deallocate(grid%ang   )
    deallocate(grid)
  end subroutine dvr_grid_destroy


  subroutine dvr_grid_write(unit,grid)
    type(dvr_grid) :: grid
    integer :: unit
    write(unit=unit) grid%xmass1 ! vibrational masses of the atoms in atomic units
    write(unit=unit) grid%xmass2 ! vibrational masses of the atoms in atomic units
    write(unit=unit) grid%xmass3 ! vibrational masses of the atoms in atomic units
    write(unit=unit) grid%idia
    write(unit=unit) grid%npnt1 !(1)  number of gauss-laguerre dvr points for r1
    write(unit=unit) grid%r1(:) ! radial grid in r1(npnt1)
    write(unit=unit) grid%npnt2 !(1)  number of gauss-laguerre dvr points for r2
    write(unit=unit) grid%r2(:) ! radial grid in r2(npnt2) 
    write(unit=unit) grid%nang
    write(unit=unit) grid%ang(:) !angular grid points ang_grid(nang)
  end subroutine dvr_grid_write

  subroutine dvr_grid_dump(unit,grid)
    type(dvr_grid) :: grid
    integer :: unit
    write(unit=unit,fmt=*) 'grid%xmass1',grid%xmass1 ! vibrational masses of the atoms in atomic units
    write(unit=unit,fmt=*) 'grid%xmass2',grid%xmass2 ! vibrational masses of the atoms in atomic units
    write(unit=unit,fmt=*) 'grid%xmass3',grid%xmass3 ! vibrational masses of the atoms in atomic units
    write(unit=unit,fmt=*) 'grid%idia  ',grid%idia  
    write(unit=unit,fmt=*) 'grid%npnt1 ',grid%npnt1  !(1)  number of gauss-laguerre dvr points for r1
    write(unit=unit,fmt=*) 'grid%r1(:)'  ! radial grid in r1(npnt1)
    write(unit=unit,fmt=*) grid%r1(:)  ! radial grid in r1(npnt1)
    write(unit=unit,fmt=*) 'grid%npnt2',grid%npnt2  !(1)  number of gauss-laguerre dvr points for r2
    write(unit=unit,fmt=*) 'grid%r2(:)'  ! radial grid in r2(npnt2) 
    write(unit=unit,fmt=*) grid%r2(:)  ! radial grid in r2(npnt2) 
    write(unit=unit,fmt=*) 'grid%nang',grid%nang
    write(unit=unit,fmt=*) 'grid%ang(:)' !angular grid points ang_grid(nang)
    write(unit=unit,fmt=*) grid%ang(:) !angular grid points ang_grid(nang)
  end subroutine dvr_grid_dump



  subroutine dvr_grid_read(unit,grid)
    type(dvr_grid),pointer :: grid
    integer :: unit
    allocate(grid)
    allocate(grid%xmass1)
    read(unit=unit) grid%xmass1 ! vibrational masses of the atoms in atomic units
    allocate(grid%xmass2)
    read(unit=unit) grid%xmass2 ! vibrational masses of the atoms in atomic units
    allocate(grid%xmass3)
    read(unit=unit) grid%xmass3 ! vibrational masses of the atoms in atomic units
    allocate(grid%idia)
    read(unit=unit) grid%idia
    allocate(grid%npnt1) !(1)  number of gauss-laguerre dvr points for r1
    read(unit=unit) grid%npnt1 !(1)  number of gauss-laguerre dvr points for r1
    allocate( grid%r1(grid%npnt1)) ! radial grid in r1(npnt1)
    read(unit=unit) grid%r1 ! radial grid in r1(npnt1)
    allocate(grid%npnt2) !(1)  number of gauss-laguerre dvr points for r1
    read(unit=unit) grid%npnt2 !(1)  number of gauss-laguerre dvr points for r1
    allocate( grid%r2(grid%npnt2)) ! radial grid in r1(npnt1)
    read(unit=unit) grid%r2 ! radial grid in r1(npnt1)
    allocate(grid%nang)
    read(unit=unit) grid%nang
    allocate(grid%ang(grid%nang))
    read(unit=unit) grid%ang(:) !angular grid points ang_grid(nang)
  end subroutine dvr_grid_read

  function dvr_grid_ang_build(nang,k)result(ang_)
    implicit none
    integer :: nang, k
    real(8) :: ang_(nang)
    type(ang_int_pnts) :: ang
    ang = ang_grid(nang,k)
    ang_  = ang%position
    call ang_grid_destroy(ang)
  end function dvr_grid_ang_build

  function dvr_grid_wtd_build(nang,k)result(wtd_)
    implicit none
    integer :: nang, k
    real(8) :: wtd_(nang)
    type(ang_int_pnts) :: ang
    ang = ang_grid(nang,k)
    wtd_  = ang%wheight
    call ang_grid_destroy(ang)
  end function dvr_grid_wtd_build

  function dvr_grid_ang_build_zperp(nang,jrot,k)result(ang_)
    implicit none
    integer :: nang, jrot,k
    real(8) :: ang_(nang)
    type(ang_int_pnts) :: ang
    ang = ang_grid_zperp(nang,jrot,k)
    ang_  = ang%position
    call ang_grid_destroy(ang)
  end function dvr_grid_ang_build_zperp

  function dvr_grid_wtd_build_zperp(nang,jrot,k)result(wtd_)
    implicit none
    integer :: nang, jrot,k
    real(8) :: wtd_(nang)
    type(ang_int_pnts) :: ang
    ang = ang_grid_zperp(nang,jrot,k)
    wtd_  = ang%wheight
    call ang_grid_destroy(ang)
  end function dvr_grid_wtd_build_zperp

  function dvr_grid_r2_max(grid,coord) result(r2max)
    implicit none
    type(dvr_grid) :: grid
    character(len=4),intent(in) :: coord 
    real(8), pointer :: r1,r2,xcos 
    real(8) :: r2max,r2tmp,q1,q2,q3,dummy
    integer :: i,j,k
    xmass(1) = grid%xmass1
    xmass(2) = grid%xmass2
    xmass(3) = grid%xmass3
    r2max=-1.0d0
    if(grid%idia.lt.0) then
       do k=1,grid%npnt2
          do j=1,grid%npnt1
             do i=1,grid%nang
                r2 => grid%r2(k)
                r1 => grid%r1(j)
                xcos => grid%ang(i)
                   call r_to_q(r1,r2,xcos,q1,q2,q3)
                   call q_to_j(q1,q2,q3,dummy,r2tmp,dummy,coord)
                   r2max=max(r2tmp,r2max)
             end do
          end do
       end do
    elseif(grid%idia.gt.0) then
       do k=1,grid%npnt2
          do j=1,grid%npnt1
             do i=1,grid%nang
                r2 => grid%r2(k)
                r1 => grid%r1(j)
                xcos => grid%ang(i)
                   call j_to_q(r1,r2,xcos,q1,q2,q3) 
                   call q_to_j(q1,q2,q3,dummy,r2tmp,dummy,coord)
                   r2max=max(r2tmp,r2max)
             end do
          end do
       end do

    else
       call error('DVR_GRID_R2_MAX: idia.eq.0. Aborting.')
       stop   
    end if
  end function dvr_grid_r2_max

end module class_dvr_grid

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
!***********************************************************************!


module class_kblock
  !COMMON PROPERTIES IN THE WAVEVUNCTION TYPES
  !THIS WILL RESPECT THE SERIAL VERSION OF THE DVR CODES
  !WITH THEIR VARIABLE NAMES
  use class_dvr_grid
  use class_error
  IMPLICIT NONE
  private
  type, public :: kblock
     integer, pointer :: k     ! (1)  kblock label
     integer, pointer :: jk    ! Total number of kblocks of which this one is part.
     integer, pointer :: kpar  ! (1)  kblock parity (postprocessed)
     integer, pointer :: ipar  ! (1)  
     integer, pointer :: jrot  ! (1)   total angular momentum
     integer, pointer :: neval ! (1)  number of eigenvalues
     integer, pointer :: nbass ! size of the total basis (depends on coordinate system and geometry)
!     integer, pointer :: lmin  ! one of the basis function labels
     integer, pointer :: iv(:) ! selected angular basis functions from the 2D step (J=0)
     type(dvr_grid), pointer   :: grid
     real(8), pointer :: eh(:) 
     real(8), pointer :: d(:)
  end type kblock

!KBLOCK CONSTRUCTORS/DESTRUCTORS  

  !KBLOCK I/O
  public :: kblock_write!(block)          write kblock to disk
  !type(kblock) :: block

  public :: kblock_write_head!(block)     write kblock head to disk
  !type(kblock) :: block

  public :: kblock_read!(kz,block)        constroy from file
  !type(kblock),pointer :: block
  !integer :: kz

  public :: kblock_head_read!(unit,block) constroy from file
  !type(kblock),pointer :: block
  !integer :: unit

  public :: kblock_destroy     !(block)   destructor
  !type(kblock),pointer :: block

  public :: kblock_head_destroy!(block)   partial destructor
  !type(kblock),pointer :: block

  public :: kblock_d_read!(kz,block)
  !  type(kblock), pointer :: block
  !  integer :: kz

  public :: kblock_state_read!(kz,block,state,wave) 
  !  integer,intent(in) :: kz, state
  !  type(kblock), pointer :: block
  !  real(8),intent(out) :: wave(block%nbass)
!  public :: kblock_head_dump!(block)         !dump the kblock head in ascii text 

  public :: kblock_head_check
  !subroutine kblock_head_check(jk,kmin)
  !  implicit none
  !  integer, intent(out) :: kmin, jk
  !end subroutine kblock_head_check

contains

  subroutine kblock_write(block,filename)
    implicit none
    type(kblock) :: block
    character(len=20),intent(in),optional :: filename
    if (present(filename)) then
       call kblock_write_head(block,filename)
       call kblock_write_tail(block,filename)
    else
       call kblock_write_head(block)
       call kblock_write_tail(block)
    end if
  end subroutine kblock_write

  subroutine kblock_write_head(block,infilename)
    implicit none
    type(kblock) :: block
    character(len=20),intent(in),optional :: infilename
    integer :: unit
    integer :: idia
    character(len=3) :: k_char
    character(len=200) :: filename
    logical :: exists
    write( k_char ,fmt="(I3.3)" ) block%k 
    if (present(infilename)) then
       filename=trim(infilename)//k_char//'_head.bcs'
    else
       filename='kblock'//k_char//'_head.bcs'
    end if
    inquire(file=trim(filename),exist=exists)
    if(.not.exists) then
       open(unit=100            &
            ,file=trim(filename)&
            ,status='new'       &
            ,form='unformatted')   
       write(unit=100) block%k     !(1)  kblock label
!       print*,'kblock_write_head:block%k', block%k
       write(unit=100) block%jk    !totalnumber of kblocks of which this one is part.
!       print*,'kblock_write_head:block%jk', block%jk
       write(unit=100) block%kpar  !(1)  kblock parity (postprocessed)
!       print*,'kblock_write_head:block%kpar', block%kpar
       write(unit=100) block%ipar  !(1)
!       print*,'kblock_write_head:block%ipar', block%ipar  
       write(unit=100) block%jrot  !(1)   total angular momentum
!       print*,'kblock_write_head:block%jrot', block%jrot  
       write(unit=100) block%neval !(1)  number of eigenvalues
!       print*,'kblock_write_head:block%neval', block%neval  
       write(unit=100) block%nbass ! npnt1*(npnt2-2*kpar)*
!       print*,'kblock_write_head:block%nbass', block%nbass  
       call dvr_grid_write(100,block%grid)
       write(unit=100) block%iv(:)
!       print*,'kblock_write_head:block%iv', block%iv  
!       print*,'kblock_write_head:block%grid%idia',block%grid%idia
!       print*,'kblock_write_head:block%grid%xmass1',block%grid%xmass1
!       print*,'kblock_write_head:block%grid%xmass2',block%grid%xmass2
!       print*,'kblock_write_head:block%grid%xmass3',block%grid%xmass3
!       print*,'kblock_write_head:block%grid%npnt1',block%grid%npnt1
!       print*,'kblock_write_head:block%grid%r1',block%grid%r1
!       print*,'kblock_write_head:block%grid%npnt2',block%grid%npnt2
!       print*,'kblock_write_head:block%grid%r2',block%grid%r2
!       print*,'kblock_write_head:block%grid%nang',block%grid%nang
!       print*,'kblock_write_head:block%grid%ang',block%grid%ang
       write(unit=100) block%eh(:)
!       print*,'kblock_write_head:block%eh', block%eh  
       close(100)
     else
       write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'
    end if
  end subroutine kblock_write_head

  subroutine kblock_write_tail(block,infilename)
    implicit none
    type(kblock) :: block
    character(len=20),intent(in),optional :: infilename
    integer :: i,nbass,neval
    integer :: idia
    character(len=3) :: k_char
    character(len=200) :: filename
    logical :: exists
    real(8),allocatable :: tmp(:)
    write( k_char ,fmt="(I3.3)" ) block%k 
    if (present(infilename)) then
       filename=trim(infilename)//k_char//'_tail.bcs'
    else
       filename='kblock'//k_char//'_tail.bcs'
    end if
    inquire(file=trim(filename),exist=exists)
    if(.not.exists) then
       nbass=block%nbass
       neval=block%neval
       open(unit=101,file=trim(filename),action='write',status='new',form='binary',access='direct',recl=nbass*8)   
       allocate(tmp(nbass))
       do i=1,neval
          tmp = block%d( (i-1)*nbass+1 : i*nbass )          
!          print"('writing from ',I7,'to ',I7)",(i-1)*nbass+1,i*nbass
          write(unit=101,rec=i) tmp
       end do
       deallocate(tmp)
       close(101)
    else
       write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'
    end if
  end subroutine kblock_write_tail

  subroutine kblock_head_read(kz,block,infilename)
    type(kblock), pointer :: block
    character(len=20),intent(in),optional :: infilename
    integer :: kz
    character(len=3) :: k_char
    character(len=200) :: filename
    logical :: exists
    write( k_char ,fmt="(I3.3)" ) kz 
    if (present(infilename)) then
       filename=trim(infilename)//k_char//'_head.bcs'
    else
       filename='kblock'//k_char//'_head.bcs'
    end if
    inquire(file=trim(filename),exist=exists)
    if(exists) then
       open(unit=100            &
            ,file=trim(filename)&
            ,status='old'       &
            ,form='unformatted')   
    allocate(block)
    allocate(block%k)
    read(unit=100) block%k     !(1)  kblock label
    allocate(block%jk)
    read(unit=100) block%jk !totalnumber of kblocks of which this one is part.
    allocate(block%kpar)
    read(unit=100) block%kpar  !(1)  kblock parity (postprocessed)
    allocate(block%ipar)
    read(unit=100) block%ipar  !(1)  
    allocate(block%jrot)  !(1)   total angular momentum
    read(unit=100) block%jrot  !(1)   total angular momentum
    allocate(block%neval)
    read(unit=100) block%neval !(1)  number of eigenvalues
    allocate(block%nbass)
    read(unit=100) block%nbass ! npnt1*(npnt2-2*kpar)*
!bruno
!print*,'block%nbass',block%nbass
!    allocate(block%lmin)
!    read(unit=100) block%lmin ! npnt1*(npnt2-2*kpar)*
    call dvr_grid_read(100,block%grid)
    allocate(block%eh(block%neval))
    allocate( block%iv( block%grid%nang ) )
    read(unit=100) block%iv ! npnt1*(npnt2-2*kpar)*
    read(unit=100) block%eh 
       close(100)
    else
       call error('KBLOCK_HEAD_READ:file '//trim(filename)//' does not exist')
    end if
  end subroutine kblock_head_read

!  subroutine kblock_head_dump(block)
!    type(kblock), pointer :: block
!    character(len=3) :: k_char
!    write(k_char,'(I3.3)') block%k-1
!    open(unit=100,file='kblock'//k_char//'_head_dump.bcs')
!    write(unit=100,fmt=*) 'block%k     ', block%k    !(1)  kblock label
!    write(unit=100,fmt=*) 'block%jk    ', block%jk   !totalnumber of kblocks of which this one is part.
!    write(unit=100,fmt=*) 'block%kpar  ', block%kpar !(1)  kblock parity (postprocessed)
!    write(unit=100,fmt=*) 'block%ipar  ', block%ipar !(1)  
!    write(unit=100,fmt=*) 'block%jrot  ', block%jrot !(1)   total angular momentum
!    write(unit=100,fmt=*) 'block%neval ', block%neval!(1)  number of eigenvalues
!    write(unit=100,fmt=*) 'block%nbass ', block%nbass! npnt1*(npnt2-2*kpar)*
!    call dvr_grid_dump(100,block%grid)
!    write(unit=100,fmt=*) 'block%iv ', block%iv! npnt1*(npnt2-2*kpar)*
!    write(unit=100,fmt=*) 'block%eh ', block%eh
!    close(100)
!  end subroutine kblock_head_dump

subroutine kblock_head_check(jk,kmin)
  implicit none
  integer, intent(out) :: kmin, jk
  logical :: exists
  character(len=200) :: filename
  integer :: unit
  character(len=3) :: k_char
  unit=100
   ! get the number of k blocks from the file
  kmin=0
  write(k_char,fmt="(I3.3)") kmin
  filename='kblock'//k_char//'_head.bcs'
  inquire(file=filename,exist=exists)
  if(exists)then 
     write(6,*)'file '//trim(filename)//' exists'
     open(unit=unit,file=trim(filename),status='old',form='unformatted')
     read(unit=unit)
     read(unit=unit) jk
     close(unit=unit)
     if (jk==1) then 
        write(6,*)'this means fort.26 wavefunction file'
     else
        write(6,*)'this means fort.8 wavefunction file'
     end if
  else
     kmin=1
     write(k_char,fmt="(I3.3)") kmin
     filename='kblock'//k_char//'_head.bcs'
     inquire(file=filename,exist=exists)
     open(unit=unit,file=trim(filename),status='old',form='unformatted')
     read(unit=unit)
     read(unit=unit) jk
     close(unit=unit)
     if(exists)then 
        write(6,*)'file '//trim(filename)//' exists'
        write(6,*)'this means fort.9 wavefunction file'
     else
        call error('cap_matrix_set_size_from_kblock_files: leading kblock file not found!')
     end if
  end if
end subroutine kblock_head_check

  subroutine kblock_d_read(kz,block,infilename)
    type(kblock), pointer :: block
    integer :: kz
    character(len=20),intent(in),optional :: infilename
    character(len=3) :: k_char
    character(len=200) :: filename
    integer :: nbass,neval,i
    real(8),allocatable :: tmp(:)
    nbass=block%nbass
    neval=block%neval
    allocate(tmp(nbass))
    write(k_char,fmt="(I3.3)") kz
    if (present(infilename)) then
       filename=trim(infilename)//k_char//'_tail.bcs'
    else
       filename='kblock'//k_char//'_tail.bcs'
    end if
    open(unit=100,file=trim(filename),status='old',form='binary',access='direct',recl=nbass*8)   
    allocate(block%d(nbass*neval))
    do i=1,neval
       read(unit=100,rec=i) tmp
       block%d((i-1)*nbass+1:i*nbass) = tmp
    end do
    close(100)
  end subroutine kblock_d_read

  subroutine kblock_state_read(kz,block,state,wave,infilename)
    integer,intent(in) :: kz, state
    type(kblock), pointer :: block
    real(8),intent(out) :: wave(block%nbass)
    character(len=20),intent(in),optional :: infilename
    character(len=3) :: k_char
    character(len=200) :: filename
    integer :: nbass,i,neval
    nbass=block%nbass
    neval=block%neval
    write(k_char,fmt="(I3.3)") kz
    if (present(infilename)) then
       filename=trim(infilename)//k_char//'_tail.bcs'
    else
       filename='kblock'//k_char//'_tail.bcs'
    end if
    open(unit=100,file=trim(filename),status='old',form='binary',access='direct',recl=nbass*8)   
       read(unit=100,rec=state) wave
    close(100)
  end subroutine kblock_state_read

  subroutine kblock_read(kz,block,infilename)
    type(kblock), pointer :: block
    integer :: kz
    character(len=20),intent(in),optional :: infilename
    if (present(infilename)) then
       call kblock_head_read(kz,block,infilename)
       call kblock_d_read(kz,block,infilename)
    else
       call kblock_head_read(kz,block)
       call kblock_d_read(kz,block)
    end if
  end subroutine kblock_read

  subroutine kblock_destroy(block)
    type(kblock),pointer :: block
    deallocate(block%k)
    deallocate(block%kpar)
    deallocate(block%jk)
    deallocate(block%ipar)
    deallocate(block%jrot)  !(1)   total angular momentum
    deallocate(block%neval)
    deallocate(block%nbass)
    deallocate(block%iv)
    call dvr_grid_destroy(block%grid)
    deallocate(block%eh)
    if(associated(block%d)) deallocate(block%d)
    deallocate(block)
  end subroutine kblock_destroy

  subroutine kblock_head_destroy(block)
    type(kblock),pointer :: block
    deallocate(block%k)
    deallocate(block%kpar)
    deallocate(block%jk)
    deallocate(block%ipar)
    deallocate(block%jrot)  !(1)   total angular momentum
    deallocate(block%neval)
    deallocate(block%nbass)
!    deallocate(block%lmin)
    call dvr_grid_destroy(block%grid)
    deallocate(block%eh)
    deallocate(block)
  end subroutine kblock_head_destroy

end module class_kblock



!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
!***********************************************************************!




module class_wavefile

  !common properties in all  wavefunction file tipes
  !in accordance with the serial version of the DVR3D suite
  !This structure gives an idea of the arrangement of the data
  !inside the file.
  !
  !Some properties like the total number of k blocks and the
  !file unit and name have been added to this structure and 
  !are not part of the wavefinction file structure
  !
  !It can also be treated as the global data, thus removing the
  !need for any auxiliary class()

  use class_error
  use class_ang_int_pnts
  private
  type,public :: wavefile
!  type line1
     !the declaration order must be kept!!! it will determine the reading order!!
     integer :: idia!   = 1 scattering coordinates hetronuclear diatomic
     !                  = 2 scattering coordinates homonuclear diatomic
     !                  = -1 radau coordinates hetronuclear diatomic
     !                  = -2 radau coordinates homonuclear diatomic
     integer :: ipar!   parity of basis - if |idia|=2: ipar=0 for even & =1 for odd
     integer :: lmax!   max number of angular basis functions in ket
     integer :: npnt1!  number of gauss-laguerre dvr points for r1
     integer :: npnt2!  number of gauss-laguerre dvr points for r2
     integer :: jrot!   total angular momentum
     integer :: kmin!   kmin=1 for sym. rotational basis, =0 for anti-sym. (kparity)
     !                  for non-coriolis calculations, kmin= k.
     integer :: neval!  number of eigenvalues
!  end type line1
!  type line2
     logical :: zembed! used only if j>0 and, zbisc=f and zperp=f     ! = t z-axis ambedded along r2      ! = f z-axis ambedded along r1
     logical :: zmorse1! = t use morse oscilator functions for the r1 coordinate 

     logical :: zmorse2! = t use morse oscilator functions for the r2 coordinate
     ! = f use spherical oscilator functions for the r2 coordinate
     real(8) :: xmass1 ! vibrational masses of the atoms in atomic units
     real(8) :: xmass2 ! vibrational masses of the atoms in atomic units
     real(8) :: xmass3 ! vibrational masses of the atoms in atomic units
     real(8) :: g1    ! generalised coordinate parameters 
     real(8) :: g2    ! (check against idia and xmass(3) ? )
     logical :: zncor ! coriolis problem
!  end type line2
!  type line3
     real(8) :: re1   ! basis function parameters 
     real(8) :: diss1 ! (read the dvr3d paper for reference)
     real(8) :: we1
     real(8) :: re2
     real(8) :: diss2
     real(8) :: we2
!  end type line3
!  type line4
     integer :: mbass0                !maximum size of vibrational problem
      integer, pointer :: lmin(:)
     integer, pointer :: lbass(:)
     integer, pointer :: nbass(:)     !basis size in each kblock
!  end type line4
     ! line 5 (possibly)
     real(8), pointer :: r1(:) ! radial grid in r1(npnt1)
     ! line 6 (if idia > -2 ) yep... they really had to save that 1/2 kbyte >:(
     real(8), pointer :: r2(:) ! radial grid in r2(npnt2) 
     ! line 7 (if idia > -2 ) otherwise 6 ... yeah... that bad...
     ! line 8 (if idia > -2 ) otherwise 7 ... yeah... that bad...
     !integer :: neval
     ! line 9 (if idia > -2 ) otherwise 8 ... yeah... that bad...
     real(8), pointer :: ang(:) ! angular grid points
     real(8), pointer :: wtd(:) ! angular grid wheights

     integer :: kz,maxleg,nidvr,lincr ! not documented quantities

     integer :: iang  ! angular basis size (J=0 only )
     integer :: ibass ! full basis size (J=0 only )
     integer, pointer :: iv(:)  ! 2d vectors selected       !the number of non-zero elements in 
                                ! per grid point (J=0 only )!this vector must = iang2          
     real(8), pointer :: eh(:)  ! egenvalues
     ! line 10 or 9 and the remaing ones
     real(8), pointer :: plegd(:)
     real(8), pointer :: d(:,:) ! (neval,ibass)
     ! alternatively:
     ! the array d consists of neval eigenvectors. thus d(1:nbass(k))
     ! represents the 1st eigenvector; d( nbass(k)+1:2*nbass(k) ) the second,
     ! etc. each eigenvector consists of the values of the wavefunction as
     ! expressed on the dvr grid in the 3 coordinates, written as a 3d
     ! array. thus if the eigenvector is placed an array vec(nalf, npnt1,
     ! npnt2), nalf, npnt1, npnt2 will give the location on the dvr grid for
     ! that value of the wavefunction.
     !the vectors are packed in a triangular form for r1 and r2 so that 
     ! ---> the total number of points is 
     !      nr1*(nr2-2*kpar)/2*nang,
     !      where kpar is 0 for even and 1 for odd. 
     !     real(8),pointer :: ang_grid(:,:) !angular grid points ang_grid(nang,kz)
     ! where jk=jrot+kmin 

     ! additional useful information and derived quantities
     integer :: ifile ! file unit associated with the wavefunction 
     integer :: jk    ! total number of kblocks inside the wavefunction file jk=jrot+kmin
     integer :: ncoord! number of free coordinates: 3 triatomic, 2 atom dissociating from rigid triatomic
     logical :: zbisc ! bisector embedding
     integer :: npnt  ! maximum size of any radial coordinate
     integer :: lpot  ! maximum order of legendre polynomials in angular integration (J>0)
     integer :: nrade ! number of radial basis set in even case
     integer :: nrado ! number of radial basis set in odd case
     integer :: mbass ! maximum size of full basis (counted through all k blocks)
     integer :: npot  ! number of gauss-legendre integration points  (angular coordinate)
     integer :: nn2   ! npot/2 (integer division)
     !genind stuff
     integer :: nbmax ! maximum size of the total basis for all k blocks
     integer :: ndbass ! basis size after defining the angular integration grid
  end type wavefile
  public :: wavefile_open
  !subroutine wavefile_open(file,w)
  !  character(len=200),optional :: file
  !  type(wavefile), pointer :: w
  !end subroutine wavefile_open
  public :: wavefile_close
  !subroutine wavefile_close(w)
  !  type(wavefile), pointer :: w
  !end subroutine wavefile_close
  public :: wavefile_setup
  public :: wavefile_d_integration_setup
  public :: wavefile_destroy
  public :: wavefile_head_dump
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wavefile_open(filename,w)
    implicit none
    character(len=200) :: filename
    type(wavefile), pointer :: w
    allocate(w)
    w%ifile=26 !set the wavefunction unit to 26 throughout the code
    write(6,*)'...opening existing wavefunction file '//trim(filename)
    open(unit=w%ifile,file= trim(filename) ,form='unformatted',status='old')
    write(6,*)'...done!'
    return
  end subroutine wavefile_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wavefile_close(w)
    implicit none
    type(wavefile), pointer :: w
    write(6,*) '...Closing wavefunction file'
    close(unit=w%ifile)
    write(6,*) '...Closed!'
  end subroutine wavefile_close
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine wavefile_setup(w)
    implicit none
    type(wavefile) :: w
    type(ang_int_pnts) :: ang
    !auxiliary variables
    integer :: jt,ltop,i,ii,k
    
    read(w%ifile) w%idia,w%ipar,w%lmax,w%npnt1,w%npnt2,w%jrot,w%kmin,w%neval
    read(w%ifile) w%zembed,w%zmorse1,w%zmorse2,w%xmass1,w%xmass2,w%xmass3,w%g1,w%g2,w%zncor
    read(w%ifile) w%re1,w%diss1,w%we1,w%re2,w%diss2,w%we2
    
    ! determine the type of triatomic calculation ( rigid / non-rigid rotor)
    w%ncoord=3
    if (w%idia .gt. -2) then
       if (w%npnt1 .eq. 1) w%ncoord=2
       w%zbisc=.false.
    else
       w%zbisc=.true.
       w%zembed=.true.
    endif

    ! get the maximum number of radial points in any radial coordinate
    w%npnt = max(w%npnt1,w%npnt2)

    !setup the number of legendre polynimials used in each J case ( = or > 0)
    if ( w%jrot == 0 ) then 
       w%lpot = w%lmax
       w%kmin = 1
    else
       !set a sufficiently high order of legendre polynomials (j>0 gaussian quadrature)
       w%lpot = ( w%lmax + w%jrot + mod( w%lmax + w%jrot , 2 ) )
    end if

    !setup coriolis problem cases
    if (w%zncor) then
       w%jk= 1
       w%jrot=abs(w%jrot)
    else
       w%jk= w%jrot + w%kmin
    endif

    ! determine radial basis and total basis sizes depending on coordinates and 
    ! symmetries
    if (w%idia .gt. -2) then
       w%nrade = w%npnt1 * w%npnt2
       w%nrado = w%nrade
       w%mbass = w%nrade * w%jk * w%lmax
    else
       !triangular pakcking of radial matrix including diagonal (even)
       w%nrade = w%npnt1 * ( w%npnt1 + 1 ) / 2
       !triangular pakcking of radial matrix excluding diagonal (odd)
       w%nrado = w%npnt1 * ( w%npnt1 - 1 ) / 2 
       jt = w%jk / 2
       w%mbass = ( w%nrade + w%nrado ) * ( jt )
       ! test the parity of jk
       if ( 2 * ( jt ) /= w%jk ) then
          if ( w%ipar == 0 ) w%mbass = w%mbass + w%nrade
          if ( w%ipar /= 0 ) w%mbass = w%mbass + w%nrado
       endif
       w%mbass = w%mbass * w%lmax ! this is calculated for checking purposes
    endif
!bruno
!print*, 'ola'
!print*, w%mbass
!stop

    ! setup space for angular integration
    w%npot = w%lpot
    if ( w%jrot > 0 ) then
       ltop = w%lpot
       if ( w%idia == 2 .and. mod( ltop , 2 ) /= w%ipar ) ltop = ltop + 1
       w%npot = ( ( ltop + 2 ) / 2 ) * 2
       w%nn2 = w%npot / 2
    endif
!bruno
!write(6,*)w%npot,w%lpot,w%nn2,ltop
!stop
    !BASIS FUNCTION SETUP

    !read the basis function labels
    allocate( w%lmin( w%jk ) , w%lbass( w%jk ) , w%nbass( w%jk ) )
    read( w%ifile ) w%mbass0, w%lmin, w%lbass, w%nbass

    !check the basis function count for consistency
    if (w%mbass0.gt.w%mbass)  then
       write(6,200) w%ifile,w%mbass,w%mbass0
200    format(//,5x,'basis function dimensions in error',&
            /,5x,'unit =',i3,' mbass =',i5,' expected, =',i5,' found')
       stop
    end if

    !generate the sub-index arrays and  find nbmax
    w%nbmax = w%nbass( 1 )
    do  k=2 , w%jk
       w%nbmax = max( w%nbmax , w%nbass( k ) )
    end do

    !read the radial grid points
    allocate( w%r1( w%npnt1 ) )
    read( w%ifile ) w%r1
    if (w%idia > -2)  then
       allocate( w%r2( w%npnt2 ) ) 
       read( w%ifile ) w%r2
    else if (w%idia == -2) then
       allocate( w%r2( w%npnt2 ) )
       w%r2 = w%r1
    end if

    allocate( w%ang( w%npot ) )
    allocate( w%wtd( w%npot ) )
    !LET THE FILE BRANCHING NIGHTMARE START!!:

!!! case J = 0 !!!   
write(*,*) 'w%jrot ==', w%jrot 
    if ( w%jrot == 0 ) then
       !read the angular grid points
       if ( w%idia > -2 ) then
          read( w%ifile ) w%ang
          if ( w%idia == 2 ) then
             do i=1,w%lmax
                w%ang( 2 * w%lmax + 1 - i ) = - w%ang( i )
             end do
          end if
          read( w%ifile ) 
          read( w%ifile )
       else
          read( w%ifile ) w%ang
          read( w%ifile )
          read( w%ifile )
          read( w%ifile )
          read( w%ifile )
       endif

       !read energy eigenvalues 
       read( w%ifile ) w%neval
       allocate( w%eh( w%neval ) )
       read( w%ifile ) w%eh

       w%ndbass = w%nbass( 1 )
       !inside dsrd

       rewind w%ifile
       do   i=1,6
          read( w%ifile ) 
       end do
       if ( .not. w%zbisc ) read( w%ifile )
       if( w%idia > -2 ) then
          rewind w%ifile
          do i=1,7
             read( w%ifile )
          end do
          read( w%ifile ) w%kz,w%maxleg,w%nidvr,w%lincr
       else
          read( w%ifile ) w%kz,w%maxleg,w%nidvr
       end if

       !leng=(maxleg+1)*nidvr
       !call getrow(plegd,leng,ivec)
       !the legendre polynomials are at this line but need not be read for J=

       read( w%ifile )

       allocate( w%iv( max( w%npot, w%lmax ) ) )
       if ( w%zbisc ) then
          read( w%ifile ) w%iang,w%ibass
          read( w%ifile ) (w%iv(ii),ii=1,w%nidvr)
       else
          w%iang = w%nidvr
          w%ibass = w%mbass
          do i = 1, w%nidvr
             w%iv( i ) = 1
          end do
       end if
       read( w%ifile )   
       read( w%ifile )
       !FROM HERE ON THE KBLOCK CAN BE READ AS:
       !         do l=1, neval
       !           call getrow(temp,ibass,ivec)
       !           do i=1, ibass
       !             d(l,i)=temp(i)
       !           end do
       !         end do

!!! case J > 0 !!!
    elseif( w%jrot > 0 ) then

       !read energy eigenvalues 
       read( w%ifile ) w%neval

       allocate( w%eh( w%neval ) )

 
       read( w%ifile ) w%eh

       ! now the asleg stuff: defining the angular integration grid:
       allocate( ang%position( w%npot ) , ang%wheight( w%npot ) )
       if ( w%idia > -2 ) then 
          ang = ang_grid(w%npot, 0 )
       else 
          ang = ang_grid(w%npot, 0 )
          !this is only for the z perpendicular case... implement input flag
          !ang = ang_grid_zperp(w%npot,w%jrot, 0 )
       end if
       w%ang = ang%position
       w%wtd = ang%wheight
!bruno
!print*,ang%position
!stop
       !inside dsrd
       if ( w%jk == 1 ) then
          rewind w%ifile
          do   i=1,6
             read( w%ifile )
          end do
          if ( .not. w%zbisc ) read( w%ifile )
          if( w%idia > -2 ) then
             rewind w%ifile
             do i=1,7
                read( w%ifile )
             end do
             read( w%ifile ) w%kz,w%maxleg,w%nidvr,w%lincr
          else
             read( w%ifile ) w%kz,w%maxleg,w%nidvr
          end if
          
          !leng=(w%maxleg+1)*w%nidvr
          !call getrow(plegd,leng,ivec)
          w%ndbass = w%nbass( 1 ) * w%npot / w%lmax
          allocate( w%plegd( w%ndbass )  )
          read( w%ifile ) w%plegd
          
          allocate( w%iv( max( w%npot, w%lmax ) ) )
          if ( w%zbisc ) then
             read( w%ifile ) w%iang,w%ibass
             read( w%ifile ) (w%iv(ii),ii=1,w%nidvr)
          else
             w%iang = w%nidvr
             w%ibass = w%ndbass
             do i = 1, w%nidvr
                w%iv( i ) = 1
             end do
          end if
          read( w%ifile )   
          read( w%ifile )
          !FROM HERE ON THE KBLOCK CAN BE READ AS:
          !         do l=1, neval
          !           call getrow(temp,ibass,ivec)
          !           do i=1, ibass
          !             d(l,i)=temp(i)
          !           end do
          !         end do
       else if ( w%jk > 1 ) then
          !FROM HERE ON THE KBLOCKs CAN BE READ AS:
          !w%d=0.0d0
          !read ( w%ifile ) ((d(i,j),j=i,nread),i=1,neval)
          ! meaning one kblock per line

          !unlike the J=0 case, these blocks have to be integrated using 
          !the gaussian quadrature for the angular coordinate.
       end if
       
    end if
 

    


 end subroutine wavefile_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine wavefile_destroy(w)
   implicit none
   type(wavefile), pointer :: w
   if(associated(w%r1))deallocate(w%r1)
   if(associated(w%r2))deallocate(w%r2)
   if(associated(w%ang))deallocate(w%ang)
   if(associated(w%wtd))deallocate(w%wtd)
   if(associated(w%plegd)) deallocate(w%plegd)    
   if(associated(w%iv)) deallocate(w%iv)    
   if(associated(w%lmin))deallocate(w%lmin)
   if(associated(w%lbass))deallocate(w%lbass)
   if(associated(w%nbass))deallocate(w%nbass)
   if(associated(w%eh))deallocate(w%eh)
   if(associated(w%d))deallocate(w%d)
 end subroutine wavefile_destroy
 
 
 subroutine wavefile_d_integration_setup(nk,w)
   implicit none
   integer, intent(in) :: nk ! number of the kblock beying transformed
   type(wavefile) :: w
   
   real(8), allocatable :: plegd(:,:),tmp(:),temp(:,:),wtds(:)
    real(8), allocatable :: plegpb(:,:),pnorm(:)
    integer :: jdia,kz,nrad,nang2,npot2,lnang,i,j,jj,l,ipt,jpt, inpot, iang, ieval, mn
    real(8) :: beta, summ, x1=1.0d0
    

    ! for J = 0 the wavefunctions are in full DVR so they are ready to integrate
    if( w%jrot == 0 ) then
       allocate( tmp( w%nbass( nk ) ) )
       allocate(w%d(w%neval,w%ibass))
       do l=1, w%neval
          read( w%ifile ) tmp
          do i=1, w%ibass
             w%d(l,i)=tmp(i)
          end do
       end do
       deallocate(tmp)
!bruno
!print*, 'd checksum',sum(w%d)
!stop
       return
    end if
!bruno
!print*, 'not J0' 
!stop   
    ! for J > 0 the wavefunctions must be prepared for integration since they are
    ! in FBR in the angular coordinate.
    ! Define the problem size quantities from the wavefile structure
    w%ndbass = w%nbass( nk ) * w%npot / w%lmax
    w%ibass  = w%nbass( nk ) * w%neval
    kz = nk - w%kmin


      !allocate the wavefunctions and read them 
      
    if ( .not. associated(w%d) ) then
       allocate( w%d( w%neval , w%ndbass ) )
    else
       deallocate( w%d)
       allocate( w%d( w%neval , w%ndbass ) )
    end if
    allocate( temp( w%neval , w%ndbass ) )
!bruno
!print*,'size(w%d,1)',size(w%d,1),'size(w%d,2)',size(w%d,2)


    temp=0.0d0
    read( w%ifile ) &
         ( ( temp( i , j ) , j = 1 , w%nbass( nk ) ) , i = 1 , w%neval )
    
!bruno
!print*, 'nread ', w%nbass(nk), ' neval ', w%neval
!print*, 'd checksum', sum(temp)    
!stop        

    w%d = 0.0d0 !initialise the new wavefunction vector to 0
    
    nrad = w%ndbass / w%npot
    jdia = max(1,w%idia)
      lnang = w%lmax*jdia - 1
      if( w%idia == 2 ) lnang = lnang + w%lmin(nk) - kz
      nang2 = w%lmax
      npot2 = w%npot / jdia
      allocate(plegd(w%npot,0:lnang)) 
      allocate(plegpb(w%npot,nang2))
      allocate(pnorm(0:lnang))
      allocate(wtds(w%npot))
      plegd = 0.0d0
      plegpb = 0.0d0
      wtds=sqrt(w%wtd)
      if (w%idia > -2 .and. w%jk == 0 ) kz = 0
      call asleg(plegd,pnorm,lnang,w%ang,w%npot,kz)
      if ( w%idia == 2 ) then
         do i= 1 , w%npot
            jj = w%lmin( nk ) - kz +1
            do j = 1 , nang2
               plegpb(i,j) = wtds( i )*plegd( i , jj - 1 ) 
               jj=jj+jdia
            end do
         end do
      else
         do i = 0 , lnang
            do j = 1 , w%npot
               plegpb(j,i+1) = wtds(j)*plegd(j,i)
            end do
         end do
      end if
      
      !        write(65,*)nang-1,npot
      !        do i=1,nang2
      !          do j=1,npot
      !            write(65,*)i,j,plegpb(i,j)
      !         end do
      !      end do
      
      !        orthonormality test
      !         do k1=1,nang2
      !         do k2=1,nang2
      !         summ = x0
      !         do j=1,npot
      !         summ = summ + plegd(j,k1)*plegd(j,k2)
      !         end do
      !         write(6,*) 'rest: k1,k2 =',k1,k2,'summ plegd  = ', summ
      !         end do
      !         end do
      
      beta=x1
      ipt=0
      jpt=0
      
      ! write(6,*)neval,npot,nang2
      do  mn = 1 , nrad
         do ieval = 1 , w%neval
            do inpot = 1 , w%npot
               summ = 0.0d0
               do iang = 1 , nang2
                  summ = summ + temp( ieval , jpt + iang ) * plegpb( inpot , iang ) 
               end do
               w%d( ieval , ipt + inpot ) = summ
            end do
         end do
!        call dgemm('n','t',neval,npot,nang2,beta,temp(1,jpt),neval,plegpb,&
!                npot,beta,d(1,ipt),neval)
         ipt = ipt + w%npot
         jpt = jpt + nang2
      end do
      deallocate( temp , plegd , plegpb )
!bruno
!print*, 'd checksum', sum(w%d)    
!stop        
 
     return
    end subroutine wavefile_d_integration_setup

      subroutine asleg(pleg,pnorm,lmax,x,nn2,m)

!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, numerical recipes, p. 182,
!     for the polynomial part of associated legendre functions.
!     a factor of sin(theta)**m has NOT been removed from all functions.

      implicit double precision (a-h,o-z)

      double precision, dimension(nn2,0:lmax) :: pleg
      double precision, dimension(nn2) :: x
      double precision, dimension(0:lmax) :: pnorm

      data x1/1.0d0/,x2/2.0d0/

      if (m .lt. 0) goto 999
      do 10 i=1,nn2
      if (abs(x(i)) .gt. x1) goto 999
      pmm = x1
      fact = x1
      somx2=sqrt((x1-x(i))*(x1+x(i)))
      do j=1,m
        pmm = -pmm * fact * somx2
        fact = fact + x2
      end do
      pleg(i,0)= pmm
      pmmp1= x(i)*dble(m+m+1)*pmm
      pleg(i,1)= pmmp1
      ll=1

      do l= 2+m,lmax+m
        r2lm1 = dble(l+l-1)
        rlpmm1= dble(l+m-1)
        rlmm  = dble(l-m)
        pll= (x(i)*r2lm1*pmmp1 - rlpmm1*pmm)/rlmm
        pmm= pmmp1
        pmmp1= pll
        ll=ll+1
        pleg(i,ll)= pll
      end do
10    continue

!     set up the normalisation constants
!     (pnorm)**2 = (2j + 1)/2   *   (j - k)! / (j + k)!
      jj = -1
      do 13 j = m,lmax+m
      fact = x1
      do i = j-m+1,j+m
        facti = dble(i)
        fact = fact * facti
      end do
      jj = jj + 1
      pnorm(jj) = sqrt(dble(j+j+1) / (fact + fact))
   13 continue
!     now normalise the polynomials
      do jj=0,lmax
        do i=1,nn2
          pleg(i,jj) = pleg(i,jj) * pnorm(jj)
        end do
      end do
      return
999   write(6,200)
200   format(//5x,'improper argument in subroutine asleg'/)
      stop
    end subroutine asleg



  subroutine wavefile_head_dump(wavefilename)
    implicit none
    type(wavefile), pointer :: w
    integer :: kz,k
    character(len=3) :: k_char
    character(len=200),intent(in) :: wavefilename
    character(len=200) :: filename
    logical :: exists
    character(len=200) :: format
    filename='wavefile_head_dump.bcs'
    inquire(file=trim(filename),exist=exists)
    if(.not.exists) then
       write(6,*)
       write(6,"(10x,'Printing wave dump file to help setup input file...')")
       call wavefile_open(wavefilename,w)
       call wavefile_setup(w)
!       read(w%ifile) w%idia,w%ipar,w%lmax,w%npnt1,w%npnt2,w%jrot,w%kmin,w%neval
!       read(w%ifile) w%zembed,w%zmorse1,w%zmorse2,w%xmass1,w%xmass2,w%xmass3,w%g1,w%g2,w%zncor
!       read(w%ifile) w%re1,w%diss1,w%we1,w%re2,w%diss2,w%we2

!       open(unit=100            &
!            ,file=trim(filename)&
!            ,status='new'       &
!            ,form='formatted')   

       write(100,*) '--- line1 ---'
       format="('idia =',2x,i2/"//&
            " 'ipar =',3x,i1/"//&
            " 'lmax =',1x,i3/"//&
            " 'npnt1=',1x,i3/"//&
            " 'npnt2=',1x,i3/"//&       
            " 'jrot =',1x,i3/"//&       
            " 'kmin =',2x,i2/"//&       
            " 'neval=',i4)"       
       write(100,fmt=format) w%idia,w%ipar,w%lmax,w%npnt1,w%npnt2,w%jrot,w%kmin,w%neval
       write(100,*) 'jk   = jrot+kmin:'
       write(100,"('jk=',1x,i3)") w%jk
       format="('zmbed  =',3x,L/"//&
            " 'zmorse1=',3x,L/"//&
            " 'zmorse2=',3x,L/"//&
            " 'xmass1 =',1x,f12.8/"//&
            " 'xmass2 =',1x,f12.8/"//&
            " 'xmass3 =',1x,f12.8/"//&
            " 'g1     =',1x,f12.8/"//&       
            " 'g2     =',1x,f12.8/"//&       
            " 'zncor  =',3x,L/)"
       WRITE(100,*)'r1     ='
       do k=1,w%npnt1
          write(100,*)  w%r1(k)
       end do
       WRITE(100,*)'r2     ='
       do k=1,w%npnt2
          write(100,*)  w%r2(k)
       end do
       write(100,fmt=format) w%zembed,w%zmorse1,w%zmorse2,w%xmass1,w%xmass2,w%xmass3,w%g1,w%g2,w%zncor
       format="('re1  =',3x,f12.8/"//&
            " 'diss1=',3x,f12.8/"//&
            " 'we1  =',3x,f12.8/"//&
            " 're2  =',3x,f12.8/"//&
            " 'diss2=',3x,f12.8/"//&
            " 'we2  =',3x,f12.8/)"
       write(100,fmt=format) w%re1,w%diss1,w%we1,w%re2,w%diss2,w%we2
       format="('mbass0  =',3x,I9)"
       write(100,fmt=format) w%mbass0
       do k=1,w%jk
          format="('for k=',I2,', lmin =',I5', lbass =',I5', nbass =',I5)"
          write(100,format) k,w%lmin(k),w%lbass(k),w%nbass(k)
       end do
       format="('neval =',I6)"
       write(100,*)
       write(100,fmt=format) w%neval
       do k=1,w%neval
          write(100,*) w%eh(k)!*eh_to_cm1
       end do
       call wavefile_close(w)
       call wavefile_destroy(w)
    end if
    return
  end subroutine wavefile_head_dump



end module class_wavefile

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
!***********************************************************************!

module class_wavefile_to_kblock
  use class_wavefile
  use class_kblock

contains

  subroutine wavefile_to_kblock_files(wavefilename,blockfilename)
    !this subroutine transforms the wavefunctions from a dvr3d/rotlev
    !calculation into purely dvr wavefunction blocks that can be used for
    !plotting or for calculation of expectation values without
    !complicated "if" constructs, in an object oriented fashion
    !It is based on the xpect3 element of the DVR3D suite, therefore the
    !conversion code is somewhat convoluted
    !the conversion will be made in memory for each kblock file.
    implicit none
    character(len=200),intent(in) :: wavefilename
    character(len=200), optional, intent(in) :: blockfilename
    type(wavefile), pointer :: w
    type(kblock), pointer :: b
    character(len=200) :: filename
    !auxiliary variables
    integer :: jt,ltop,nk
    character(len=3) :: k_char
    logical :: exists
    !allocates the wavefile structure and opens the file,
    !with unit stored in w%ifile    


    write( k_char ,fmt="(I3.3)" ) 0 
    if (present(blockfilename)) then
       filename=trim(blockfilename)//k_char//'_tail.bcs'
    else
       filename='kblock'//k_char//'_tail.bcs'
    end if

    inquire(file=trim(filename),exist=exists)

    if(.not.exists) then

       write( k_char ,fmt="(I3.3)" ) 1 
       if (present(blockfilename)) then
          filename=trim(blockfilename)//k_char//'_tail.bcs'
       else
          filename='kblock'//k_char//'_tail.bcs'
       end if

       inquire(file=trim(filename),exist=exists)

       if(.not.exists) then

          call wavefile_open(wavefilename,w)
          call wavefile_setup(w)

          do nk = 1 , w%jk
             !print*, 'wavefile_to_kblock_files: nk =', nk
             call wavefile_d_integration_setup(nk,w)
!bruno
!print*,'bla'
!stop

             call wavefile_to_kblock_file(w,nk,b)
             !       call kblock_write(b)
             !print*,'wavefile_destroy: before deallocate(b%k)' 
             deallocate(b%k)
             !print*,'wavefile_destroy: before deallocate(b%kpar)' 
             deallocate(b%kpar)
             !print*,'wavefile_destroy: before deallocate(b%grid)' 
             deallocate(b%grid)
             !print*,'wavefile_destroy: before deallocate(b%iv)' 
             nullify(b%iv)
             !print*,'wavefile_destroy: before deallocate(b%d)' 
             !       deallocate(b%d)
             !print*,'wavefile_destroy: before deallocate(b)' 
             deallocate(b)
             !print*,'wavefile_destroy: after deallocate(b)' 
          end do

          call wavefile_destroy(w)    
       else
          write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'
       end if

    else
       write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'
    end if




  end subroutine wavefile_to_kblock_files

  subroutine wavefile_to_kblock_file(w,nk,b,infilename)

    implicit none
    type(wavefile),target,intent(in) :: w
    type(kblock),pointer :: b
    integer,intent(in) :: nk
    character(len=20),intent(in),optional :: infilename

    character(len=3) :: k_char
    character(len=200) :: filename
    logical :: exists
    real(8),allocatable :: tmp(:)
    integer :: kz,i,j,nbass,neval
    real(8) :: c
    
       kz =  nk - w%kmin
       allocate(b)
       allocate(b%k)
       b%k = kz
       b%Jk => w%jk
       allocate(b%kpar)

       b%ipar => w%ipar
       b%jrot => w%jrot

       b%neval => w%neval
       b%nbass => w%ndbass
       allocate(b%grid)
       b%grid%xmass1 => w%xmass1
!print*,'wavefile_to_kblock:b%grid%xmass1',b%grid%xmass1
       b%grid%xmass2 => w%xmass2
!print*,'wavefile_to_kblock:b%grid%xmass2',b%grid%xmass2
       b%grid%xmass3 => w%xmass3
!print*,'wavefile_to_kblock:b%grid%xmass3',b%grid%xmass3
       b%grid%idia  => w%idia
!print*,'wavefile_to_kblock:b%grid%idia',b%grid%idia
       b%grid%npnt1=> w%npnt1
!print*,'wavefile_to_kblock:b%grid%npnt1',b%grid%npnt1
       b%grid%r1 => w%r1
!print*,'wavefile_to_kblock:b%grid%r1',b%grid%r1
       b%grid%npnt2 => w%npnt2
!print*,'wavefile_to_kblock:b%grid%npnt2',b%grid%npnt2
       b%grid%r2=> w%r2
!print*,'wavefile_to_kblock:b%grid%r2',b%grid%r2
       b%grid%nang => w%npot
!print*,'wavefile_to_kblock:b%grid%nang',b%grid%nang
       b%grid%ang => w%ang 
!print*,'wavefile_to_kblock:b%grid%ang',b%grid%ang
!print*,'wavefile_to_kblock: before iv'

       !define the block parity
       if( b%nbass/b%grid%nang == w%nrade ) b%kpar = 0
       if( b%nbass/b%grid%nang == w%nrado ) b%kpar = 1
       if( b%jrot == 0) b%kpar = b%ipar

       if(associated(w%iv))then
          b%iv => w%iv
       else
          allocate( b%iv( b%grid%nang ) )
          b%iv = 1
       end if
!print*,'wavefile_to_kblock: after iv'
       b%eh => w%eh

!       allocate( b%d( w%ndbass * w%neval ) )
!       c = 0.0d0
!print*,'wavefile_to_kblock: before d'
!       do i = 1 , w%neval
!          do j = 1, w%ndbass
!             c = c + 1
!             b%d( c ) = w%d(i,j)
!          end do
!       end do

   call kblock_write_head(b)

    write( k_char ,fmt="(I3.3)" ) b%k 
       filename='kblock'//k_char//'_tail.bcs'

    inquire(file=trim(filename),exist=exists)
    if(.not.exists) then
       nbass=b%nbass
       neval=b%neval
!bruno
!print*,w%ndbass
!print*,b%neval
!stop
       open(unit=101,file=trim(filename),action='write',status='new',form='binary',access='direct',recl=nbass*8)   
       allocate(tmp(w%ndbass))
    
       do i=1,neval
          tmp = w%d( i , 1 : w%ndbass )           
!          print"('writing from ',I7,'to ',I7)",(i-1)*nbass+1,i*nbass
          write(unit=101,rec=i) tmp
       end do
       deallocate(tmp)
       close(101)
    else
       write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'
    end if



!print*,'wavefile_to_kblock: after d'       
  end subroutine wavefile_to_kblock_file




end module class_wavefile_to_kblock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

