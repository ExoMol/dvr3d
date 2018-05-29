
module class_curvs
  use common_parameters
  use class_error
  use class_input 
  use class_complex_h
  private

type, public :: curvs
integer, pointer :: nlambda
integer, pointer :: nstates
real(8), pointer :: eh1
real(8), pointer :: lambda(:)
real(8), pointer :: crvr(:,:)  ! allocated as (nlambda,nstates)
real(8), pointer :: crvi(:,:) 
real(8), pointer :: dcrvi(:,:) 
real(8), pointer :: crvtr(:,:) ! curvature
end type curvs

public :: deltar
public :: curvs_build
public :: curvs_read
public :: curvs_destroy
contains


  subroutine curvs_build(in)
    implicit none
    type(input) :: in
    character(len=200) :: capfilename
    character(len=200) :: diagfilename
    character(len=200) :: diagvecfilename
    character(len=200) :: curvfilename
    integer :: i
    character(len=2) :: ncap_char
    logical :: exists
    integer,parameter ::  capunit     = 120
    integer,parameter ::  diagunit    = 121
    integer,parameter ::  diagvecunit_head = 122
    integer,parameter ::  diagvecunit_tail = 123
    integer,parameter ::  curvunit    = 124
    type(cap_matrix), pointer :: cap

    ! DIAGONALISE THE ComplexH MATRICES
write(6,*)
write(6,"(20x,'----- DIAGONALISING THE ComplexH MATRICES -----')")
write(6,*)

    do i=1,ndeltar
       write(ncap_char,fmt="(I2.2)") i
       capfilename='cap'//ncap_char//'.bcs'
       diagfilename='diag'//ncap_char//'.bcs'
       diagvecfilename='diagvec'//ncap_char
       inquire(file=trim(diagfilename),exist=exists)
       if(.not.exists) then
             open(unit=capunit,file=trim(capfilename),form='unformatted',status='old')
             call cap_matrix_read(capunit,cap)
             write(6,"(10x,'diagonalising CAP matrix with deltar = ',f5.2,' au')"),deltar(i)
             close(capunit)
             open ( unit=diagunit , file = trim(diagfilename), status = "new", form='unformatted')
             if(in%diag%outputvec) then
                open( unit=diagvecunit_head , file = trim(diagvecfilename)//'_head.bcs' , status = "new", form='formatted' )
                write(diagvecunit_head,*) in%cap%molecule
                write(diagvecunit_head,*) in%diag%nlambda
                write(diagvecunit_head,*)  diag_low_state
                write(diagvecunit_head,*)  diag_high_state
                close(diagvecunit_head) 
                open( unit=diagvecunit_tail , file = trim(diagvecfilename)//'_tail.bcs' ,&
                status = "new", &
                form='binary', &
                access='direct', &
                recl=16*(n_diag_states*n_diag_states))
!                open( unit=diagvecunit_tail , file = trim(diagvecfilename)//'_tail.bcs' ,&
!                status = "new", &
!                access='direct', &
!                recl=16*(n_diag_states*n_diag_states))
             end if
             call complex_h_diagonalise(diagunit,diagvecunit_tail,in%diag,cap)
             close(diagunit)
             if(in%diag%outputvec) close(diagvecunit_tail)
             call cap_matrix_destroy(cap)
       else
          write(6,*)'Warning: file ', trim(diagfilename),' already exists -- skipping'
       end if
    end do

    ! DO THE CURVES AND CURVATURES
    do i=1,ndeltar
       write(ncap_char,fmt="(i2.2)") i
       diagfilename='diag'//ncap_char//'.bcs'
       curvfilename='curvs'//ncap_char//'.bcs'
       inquire(file=curvfilename, exist=exists )
       if(exists) then
          write(6,*)'Warning: file ', trim(curvfilename),' already exists -- skipping'        
       else
          write(6,"(10x,'Getting Curves with deltar = ',f5.2,' au')"),deltar(i)
          open ( unit=diagunit , file = trim(diagfilename), status = "old" , action='read', form='unformatted' )
          open ( unit=curvunit , file = trim(curvfilename), status = "new" ,form='unformatted' )
          call curve(diagunit,curvunit,in)
          close(diagunit)
          close(curvunit)
       end if
    end do
    
    ! WRITE THE CURVES FOR PLOTTING (E.G. IN GLUPLOT)
    do i=1,ndeltar
       write(ncap_char,fmt="(i2.2)") i
       curvfilename='curvs'//ncap_char//'.bcs'
       inquire(file='deltar'//ncap_char//'curv_plot.bcs',exist=exists)
       if(exists) then
          write(6,*)'Warning: file deltar'//ncap_char//'curv_plot.bcs already exists -- skipping'        
       else
          write(6,"(10x,'Writing curves for plotting with deltar = ',f5.2,' au')"),deltar(i)
          open ( unit=curvunit , file = trim(curvfilename), status = "old" ,form='unformatted' )
          call curvs_plot(curvunit,i)
          close(curvunit)
       end if
       end do
end subroutine curvs_build

subroutine curvs_read(unit,c)
type(curvs), pointer :: c
integer :: unit
allocate(c)
allocate(c%nlambda)
read(unit) c%nlambda
allocate(c%nstates)
read(unit) c%nstates
allocate(c%eh1)
read(unit) c%eh1
allocate(c%lambda(c%nlambda))
read(unit) c%lambda
allocate(c%crvr(c%nlambda,c%nstates))
read(unit) c%crvr
allocate(c%crvi(c%nlambda,c%nstates))
read(unit) c%crvi
allocate(c%dcrvi(c%nlambda,c%nstates))
read(unit) c%dcrvi
allocate(c%crvtr(c%nlambda,c%nstates))
read(unit) c%crvtr
end subroutine curvs_read

subroutine curvs_write(unit,c)
type(curvs) :: c
integer :: unit
write(unit) c%nlambda
write(unit) c%nstates
write(unit) c%eh1
write(unit) c%lambda
write(unit) c%crvr
write(unit) c%crvi
write(unit) c%dcrvi
write(unit) c%crvtr
end subroutine curvs_write

subroutine curvs_plot(unit,deltar)
!deltar is just for labeling purposes
integer, intent(in) :: unit, deltar 
character(len=2) :: char_deltar
character(len=4) :: char_state
type(curvs),pointer :: c
integer :: i,j
real(8)::max_crvtr
call curvs_read(unit,c)
write(char_deltar,fmt="(I2.2)") deltar
!write(char_state,fmt="(I4.4)") i
open(unit=222,file='deltar'//char_deltar//'curv_plot.bcs',status='new') 
write(6,*)
write(6,"(10x,'Writing curve plots to file ',A)")'deltar'//char_deltar//'curv_plot.bcs'
do i = 1 , c%nstates

max_crvtr=-10000000000000.0d0
do j = 1 , c%nlambda
if(c%crvtr(j,i) > max_crvtr ) max_crvtr = c%crvtr(j,i)
end do

do j = 1 , c%nlambda
write(222,"(I4.4,1X,I5.5,5(1X,ES22.15))") i,j,eh_to_cm1*c%lambda(j),&
c%crvr(j,i)-gstate*eh_to_cm1,c%crvi(j,i),c%dcrvi(j,i),(c%crvtr(j,i)/max_crvtr)**2
end do
write(222,*)
write(222,*)
end do
close(222)

!GNUPLOT
!Index
!The index keyword allows only some of the data sets in a multi-data-set file to be plotted.
!
!Syntax:
!
!     plot 'file' index <m>{{:<n>}:<p>}
!
!
!Data sets are separated by pairs of blank records.
! index m selects only set m;
! index m:n selects sets in the range m to n;
! and index m:n:p selects indices m, m+p, m+2p, etc.,
! but stopping at n.
!Following C indexing, the index 0 is assigned to the first
! data set in the file. Specifying too large an index results
! in an error message. If index is not specified, all sets
! are plotted as a single data set.
!
!Example:
!
!     plot 'file' index 4:5
!
!
!http://www.gnuplot.info/demo/multimsh.html
!Using
!The most common datafile modifier is using.
!
!Syntax:
!
!     plot 'file' using {<entry> {:<entry> {:<entry> ...}}} {'format'}
!
!
!If a format is specified, each datafile record is read 
!using the C library's 'scanf' function, with the specified
!format string. Otherwise the record is read and broken into
!columns at spaces or tabs. A format cannot be specified this
!way for time-format data (instead use set xdata time).
end subroutine curvs_plot



subroutine curvs_destroy(c)
type(curvs), pointer :: c
deallocate(c%nlambda)
deallocate(c%nstates)
deallocate(c%eh1)
deallocate(c%lambda)
deallocate(c%crvr)
deallocate(c%crvi)
deallocate(c%dcrvi)
deallocate(c%crvtr)
deallocate(c)
end subroutine curvs_destroy


  subroutine curve(diagunit,curvunit,in)
    implicit none
    integer :: diagunit,curvunit
    type(input) :: in
     type(curvs), pointer :: c

    integer :: i,j,inu,i0,i1,i2,n,nl,ipar
    real(8) :: R,coe1,coe2,dr1,dr2,di1,di2,ww,conv,hrel,hh,xdummy
    real(8), allocatable:: cr1(:),ci1(:),cr2(:),ci2(:)
    real(8), pointer :: cr(:),ci(:)
    character :: char_state*4
    character(len=150) :: filename
    
    parameter (conv=219474.63067d0)

    
    n=diag_high_state-diag_low_state+1    ! number of states 
    nl=in%diag%nlambda                    ! number of lambdas
    
    ! 
    allocate(c)
    allocate(c%nlambda)
    c%nlambda=nl
    allocate(c%nstates)
    c%nstates=n
    allocate(c%eh1)
    c%eh1=gstate
    allocate(c%lambda(nl))
    allocate(c%crvr(nl,n),c%crvi(nl,n),c%dcrvi(nl,n),c%crvtr(nl,n))

    allocate(cr1(nl),ci1(nl),cr2(nl),ci2(nl))

    do j=1,nl
       do i=1,n
          read(diagunit)xdummy,inu,inu,c%crvr(j,i),c%crvi(j,i)
       end do
    end do

    hrel=in%diag%h*conv
    hh=in%diag%h ! paolo had originally h=1
    coe1=(in%diag%chi-1.d0)/log(in%diag%chi)
    coe2=(in%diag%chi-1.d0)/hh


    do j=1,nl
       c%lambda(j) = in%diag%h*(in%diag%chi**(j-1)-1.d0)/(in%diag%chi-1.d0)
    end do

    !curvature calculation
    do i=1,n
!          i2=i/100
!          i1=(i-i2*100)/10
!          i0=(i-i2*100-i1*10)

          cr => c%crvr(:,i)
          ci => c%crvi(:,i)

       call DRV(1,cr,cr1,1,nl-1,hh)
       call DRV(1,ci,ci1,1,nl-1,hh)
       call DRV(2,cr,cr2,1,nl-1,hh)
       call DRV(2,ci,ci2,1,nl-1,hh)

       do j=1,nl-1
          ww=1.d0/(in%diag%chi**(j-1))
          dr1=cr1(j)*coe1*ww
          di1=ci1(j)*coe1*ww
          c%dcrvi(j,i)=ci1(j)*coe1*ww!/hh
          dr2=ww*ww*coe1*coe1*cr2(j)-cr1(j)*coe1*coe2*ww*ww
          di2=ww*ww*coe1*coe1*ci2(j)-ci1(j)*coe1*coe2*ww*ww
          c%crvtr(j,i)=abs(dr1*di2-dr2*di1)/(dr1**2+di1**2)**(1.5d0)
       end do
    end do

    deallocate(cr1,ci1,cr2,ci2)

    call curvs_write(curvunit,c)
    call curvs_destroy(c)
    
  end subroutine curve


  SUBROUTINE DRV(IK,G,DG,N0,N,H)
    IMPLICIT REAL*8(A-H,O-Z)
    !       El programa calcula DG derivada primera (IK=1) o segunda (IK=2)
    !       de la G. La G es tabulada eventualmente en mas bloques pero
    !       los puntos en que se deriva deben pertenecer todos al mismo
    !       bloque.
    !       H=passo 
    !       N0=indice iniziale
    !       N=numero di punti
    DIMENSION G(1),DG(1)
    K1=N
    K2=N-1
    K3=N-2
    K4=N-3
    K5=N-4
    K6=N-5
    K7=N-6
    HR=1.D0/(60.D0*H)
    H2R=HR/(3.D0*H)
    N1=N0+1
    N2=N1+1
    N3=N2+1
    G1=G(N0)
    G2=G(N1)
    G3=G(N2)
    G4=G(N3)
    G5=G(N0+4)
    G6=G(N0+5)
    G7=G(N0+6)
    GOTO (10,15),IK
10  DG(N0)=HR*(-147.D0*G1+36.D1*G2-45.D1*G3+4.D2*G4-225.D0*G5+&
         72.D0*G6-1.D1*G7)
    DG(N1)=HR*(-1.D1*G1-77.D0*G2+15.D1*G3-5.D1*(G4+G4-G5)-15.D0*G6+&
         G7+G7)
    DG(N2)=HR*(G1+G1-24.D0*G2-35.D0*G3+8.D1*G4-3.D1*G5+8.D0*G6-G7)
    GOTO 20
15  DG(N0)=H2R*(812.D0*G1-3132.D0*G2+5265.D0*G3-508.D1*G4+297.D1*G5&
         -972.D0*G6+137.D0*G7)
    DG(N1)=H2R*(137.D0*G1-147.D0*G2-255.D0*G3+47.D1*G4-285.D0*G5+&
         93.D0*G6-13.D0*G7)
    DG(N2)=H2R*(-13.D0*G1+228.D0*G2-42.D1*G3+2.D2*G4+15.D0*G5-&
         12.D0*G6+G7+G7)
20  DO 60 K=N3,K4
       G1=G(K+1)
       G2=G(K+2)
       G3=G(K+3)
       G4=G(K-1)
       G5=G(K-2)
       G6=G(K-3)
       G7=G(K)
       GOTO(40,50),IK
40     DG(K)=HR*(45.D0*(G1-G4)+9.D0*(G5-G2)+G3-G6)
       GOTO 60
50     DG(K)=H2R*(-49.D1*G7+27.D1*(G1+G4)-27.D0*(G2+G5)+&
            2.D0*(G3+G6))
60     CONTINUE
       G1=G(K7)
       G2=G(K6)
       G3=G(K5)
       G4=G(K4)
       G5=G(K3)
       G6=G(K2)
       G7=G(K1)
       GOTO(80,85),IK
80     DG(K3)=HR*(G1-8.D0*G2+3.D1*G3-8.D1*G4+35.D0*G5+24.D0*G6-G7-G7)
       DG(K2)=HR*(-G1-G1+15.D0*G2-5.D1*(G3-G4-G4)-15.D1*G5+77.D0*G6&
            +1.D1*G7)
       DG(K1)=HR*(1.D1*G1-72.D0*G2+225.D0*G3-4.D2*G4+45.D1*G5-36.D1*G6 &
            +147.D0*G7)
       GOTO 100
85     DG(K3)=H2R*(G1+G1-12.D0*G2+15.D0*G3+2.D2*G4-42.D1*G5+228.D0*G6 &
            -13.D0*G7)
       DG(K2)=H2R*(-13.D0*G1+93.D0*G2-285.D0*G3+47.D1*G4-255.D0*G5 &
            -147.D0*G6+137.D0*G7)
       DG(K1)=H2R*(137.D0*G1-972.D0*G2+297.D1*G3-508.D1*G4+5265.D0*G5 &
            -3132.D0*G6+812.D0*G7)
100    RETURN
    END subroutine drv
    !
  end module class_curvs
