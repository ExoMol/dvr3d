module class_kblock_to_elements
  !COMMON PROPERTIES IN THE WAVEVUNCTION TYPES
  !THIS WILL RESPECT THE SERIAL VERSION OF THE DVR CODES
  !WITH THEIR VARIABLE NAMES
  use class_kblock
  use class_wavefile
  use class_error
  IMPLICIT NONE
  private

  public :: elements_from_kblock_files_mem_split
  !FUNCTION elements_from_kblock_files_mem_split( mem &
  !     , low_state , high_state , n_op,operator ) result( elements )
  !  !SPLIT THE CAP MATRIX CALCULATION INTO
  !  !CHUNKS TO MAXIMIZE AVAILABLE MEMORY USE
  !  !AND MINIMIZE I/O 
  !  real(8): : mem                      ! in MB
  !  integer :: high_state,low_state     ! highest state and lowest state defining the matrix range 
  !  real(8),allocatable :: elements(:,:) ! (N_STATES,N_OP)
  !end function elements_from_kblock_files_mem_split
  public :: diagonal_elements_from_kblock_files
  !SUBROUTINE diagonal_elements_from_kblock_files(low_state,high_state, n_op,operator,elements)
  !  implicit none
  !  integer,intent(in) :: low_state,high_state
  !  real(8),external :: operator
  !  real(8),allocatable :: op_matrix(:,:,:,:) !(nr1,nr2,nxcos,n_op)
  !  real(8),allocatable :: elements(:,:) ! (N_STATES,N_OP)
  !end function diagonal_elements_from_kblock_files 

!PRIVATE ROUTINES  
!  public :: kblock_element
  !  function block_element(block,bra,op_matrix,ket) result(element)
  !    integer :: bra, ket
  !    type(kblock) :: block
  !    real(8) :: op_matrix(:,:,:)
  !    real(8) :: element
  !  end function block_element
!  public :: kblock_element_cross
  !function kblock_element_cross(bra_block,bra,op_matrix,ket_block,ket) result(element)
  !  implicit none
  !  integer :: bra, ket
  !  type(kblock),pointer :: bra_block,ket_block
  !  real(8) :: op_matrix(:,:,:)
  !  real(8) :: element
  !end function kblock_element_cross


!public :: operator_matrix_build!(block,n_op,operator,op_matrix)
!!subroutine operator_matrix_build(block,n_op,operator,op_matrix)
!!!type(kblock),intent(in) :: block
!!integer, intent(in) :: n_op
!!real(8), external :: operator
!!real(8),intent(out) :: op_matrix(:,:,:,:)
!!end subroutine operator_matrix_build


!TESTING ROUTINES
  public :: matrix_element_from_kblock_set
  !function matrix_element(block,bra,op_matrix,ket) result(element)
  !  integer :: bra, ket
  !  type(kblock), pointer :: block(:)
  !  real(8) :: op_matrix(:,:,:)
  !  real(8) :: element
  !end function matrix_element
  public :: elements_from_kblock_files
  !function elements_from_kblock_files(low_state,high_state,op_matrix) result(elements) 
  !  integer :: low_state,high_state
  !  real(8) :: op_matrix(:,:,:)
  !  real(8),allocatable :: elements(:) 
  !end function elements_from_kblock_files

contains



  function diagonal_elements_from_kblock_files(low_state,high_state,n_op,operator) result(elements)
    implicit none
    integer,intent(in) :: low_state,high_state
    integer,intent(in) :: n_op
!    real(8),
    external::operator
    real(8) :: elements(high_state-low_state+1,n_op) 
    real(8),allocatable :: op_matrix(:,:,:,:)
    integer :: l1,l2
    type(kblock), pointer :: block
    integer :: k,kz,kmin,unit,jk,count,nstates
    real(8), allocatable :: kbe(:,:,:)
!    interface 
!       function operator(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos ; integer :: n ; real(8) :: out(n)
!       end function operator
!    end interface

!    real(8),allocatable :: allelements(:) 
!    real(8) :: mem
!    integer :: i,j,count,nstates
    nstates = high_state-low_state+1
!    allocate(allelements(nstates*(nstates+1)/2))
    elements=0.0d0

    call kblock_head_check(jk,kmin)
    allocate(kbe(nstates,jk,n_op))
    do k = 1,jk
       kz=k-(1-kmin)!note that kmin here is defined as the complementary of kmin in rotlev
       print*,'kz',kz
       call kblock_head_read( kz , block )
       allocate(block%d(block%nbass*nstates))
       call readinwaves(low_state,nstates,block)
       allocate( op_matrix( block%grid%nang,block%grid%npnt1 , block%grid%npnt2 ,n_op  ) )
       call operator_matrix_build( block,n_op , operator, op_matrix )
       count=0
       do l1=1,nstates
             count=count+1
             kbe(l1,k,:) = kblock_element(block,l1,n_op,op_matrix,l1)
             elements(count,:) = elements(count,:)  + kbe(l1,k,:)
       end do

!       mem=8000.0d0
!       allelements=elements_from_kblock_files_mem_split( mem &
!       , low_state , high_state , operator )
!
!
!
!       count=0
!       do i=1,nstates
!          do j=1,i
!             count=count+1
!             if (i == j) elements(i) = allelements(count)
!          end do
!       end do

       deallocate(op_matrix)
       call kblock_destroy(block)
    end do
do l2=1,n_op

   write(6,*) "OPERATOR", l2
    do l1=1,nstates
       write(6,*)l1,(kbe(l1,k,l2),k=1,jk)
    end do
 end do
end function diagonal_elements_from_kblock_files



  function elements_from_kblock_files_mem_split( mem &
       , low_state , high_state , n_op,operator ) result( elements )
    !SPLIT THE CAP MATRIX CALCULATION INTO
    !CHUNKS TO MAXIMIZE AVAILABLE MEMORY USE
    !AND MINIMIZE I/O 
    implicit none
    real(8),intent(in) :: mem
    integer,intent(in) :: high_state,low_state
    integer :: n_op
!    real(8), 
external :: operator
    real(8),allocatable :: op_matrix(:,:,:,:)
    real(8),allocatable :: elements(:,:)
    type(kblock),pointer :: bra_block, ket_block
    real(8) :: ms
    character(len=5) :: char_ms
    integer :: nstates,max_n_states,nb,n_bra_states,n_ket_states
    integer :: nsqb,first_bra,first_ket,jk,kmin,k,kz,i,j,prev_count 
    real(8),parameter :: u= 0.00000762939453125000
    real(8),allocatable :: cap_matrix(:,:,:),tmp_cap_matrix(:,:,:) 
!    interface 
!       function operator(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos ; integer :: n ; real(8) :: out(n)
!       end function operator
!    end interface


    nstates = high_state-low_state+1
    allocate( elements( nstates*(nstates+1)/2 , n_op ) )
    elements=0.0d0
    allocate(  cap_matrix( nstates , nstates ,n_op)    )
    allocate( tmp_cap_matrix( nstates , nstates ,n_op) ) 
    cap_matrix=0.0d0
    tmp_cap_matrix=0.0d0
    call kblock_head_check(jk,kmin)
    do k = 1,jk
       kz=k-(1-kmin)!note that kmin here is defined as the complementary of kmin in rotlev
       call kblock_head_read( kz , bra_block )!there are actually the same 
       call kblock_head_read( kz , ket_block )!kblock
       allocate( op_matrix( bra_block%grid%nang , bra_block%grid%npnt1 , bra_block%grid%npnt2, n_op ) )
       call operator_matrix_build( bra_block,n_op , operator, op_matrix )

       ms = bra_block%nbass * u   !size of a state in MB !u= 0.00000762939453125000 MB per real(8) element  
       max_n_states = 0.5d0 * mem / ms
       if(mem .lt. ms ) then
          write(char_ms,"(f5.1)") ms
          call error('ELEMENTS_FROM_KBLOCK_FILES_MEM_SPLIT: provided memory smaller than size of single state ('//char_ms//' MB )')
       end if
       max_n_states = 2 * floor( 0.5d0 * mem / ms ) ! maximum state blocks per BRA / KET
       if (nstates.le.max_n_states) then
          max_n_states = nstates 
       end if
       nb = nstates / max_n_states     ! number of integer blocks (note integer division)


       !Build Large Triangles:
       do i = 1 , nb + 1
          tmp_cap_matrix = 0.0d0
          if ( i .eq. nb + 1 ) then
             n_bra_states = mod ( nstates , max_n_states )
          else
             n_bra_states = max_n_states
          end if
          if (n_bra_states.gt.0) then
             first_bra = low_state + ( i - 1 ) * max_n_states    
             allocate ( bra_block%d( bra_block%nbass*n_bra_states ) )
             call readinwaves ( first_bra , n_bra_states, bra_block )
             call elements_triangle(bra_block,low_state,first_bra &
                  , n_bra_states,n_op, op_matrix , tmp_cap_matrix )
             deallocate ( bra_block%d )                                           !
             cap_matrix = cap_matrix + tmp_cap_matrix   
          end if
       end do !i = 1 , nb

       if (max_n_states.lt.nstates) then
          !Build Small Squares:
          nsqb =  ( nstates - max_n_states ) / (max_n_states/2) ! interger division (remainder truncated)  
          n_bra_states = max_n_states/2
          n_ket_states = max_n_states/2
          allocate ( bra_block%d( bra_block%nbass * n_bra_states ) )
          allocate ( ket_block%d( ket_block%nbass * n_ket_states ) )
          prev_count=1
          do i = 1 , nsqb + 1
             
             if ( i .eq.  nsqb + 1 ) then !this condition deals with the remainder of the matrix
                deallocate( bra_block%d )
                n_bra_states = mod(( nstates - max_n_states ),( max_n_states/2 )) !remainder computed here
                if (n_bra_states.gt.0) allocate ( bra_block%d( bra_block%nbass * n_bra_states )) 
             end if
             if (n_bra_states.gt.0) then ! if this condition ocurrs, the following code is skipped
                first_bra = low_state+(i+1)*(max_n_states/2)
                call readinwaves (first_bra , n_bra_states , bra_block ) !reads a set of eigenvectors sequencially          
                
                if ( prev_count .eq. 1 ) then
                   
                   do j = 1 , i + 1
                      tmp_cap_matrix = 0.0d0
                      first_ket = low_state + ( j - 1 )*(max_n_states/2)
                      if(debug) print*,'elements_from_kblock_files_mem_split:',first_ket , n_ket_states
                      call readinwaves ( first_ket , n_ket_states ,  ket_block )   
                      call elements_square( bra_block , ket_block ,low_state &
                           ,n_bra_states , first_bra &
                           ,n_ket_states , first_ket, n_op, op_matrix &
                           ,tmp_cap_matrix ) 
                      prev_count = j
                      cap_matrix = cap_matrix + tmp_cap_matrix
                   end do
                else if ( prev_count .eq. i ) then
                   do j =  prev_count, 1 , -1
                      tmp_cap_matrix = 0.0d0
                      call elements_square( bra_block , ket_block ,low_state&
                           , n_bra_states , first_bra &
                           , n_ket_states , first_ket, n_op, op_matrix &
                           , tmp_cap_matrix ) 
                      cap_matrix = cap_matrix + tmp_cap_matrix
                      if( j .ne. 1) then
                         first_ket = low_state + ( j - 2 )* (max_n_states/2) 
                         call readinwaves ( first_ket , n_ket_states , ket_block )   
                      end if
                      prev_count = j
                      
                   end do
                endif
             end if
          end do!i = 1 , nsqb + 1
          call kblock_destroy(bra_block)
          call kblock_destroy(ket_block)
       end if
          deallocate(op_matrix)
    end do
       
    !rearrange the cap matrix into a vector of the triangular matrix.


do k=1,n_op
    prev_count=0
    do i=1,nstates
       do j=1,i
          prev_count=prev_count+1
          elements(prev_count,k)=cap_matrix(j,i,k)
       end do
    end do
 end do
deallocate (tmp_cap_matrix,cap_matrix)
!write(unit=6,fmt=*)"printing elements"    
!do i=1,nstates*(nstates+1)/2
!write(unit=66,fmt=*) elements(i)
!end do
!stop
  end function elements_from_kblock_files_mem_split



  subroutine readinwaves ( first_state, n_states, block ) 
    implicit none
    type(kblock),pointer :: block
    integer, intent ( in  ) :: n_states, first_state
    integer                 :: i,nbass
    character(len=3) :: k_char
    real(8), allocatable :: tmp(:)

    write(unit=k_char,fmt="(I3.3)") block%k
    nbass = block%nbass
    allocate(tmp(nbass))
    open(unit=66, &
         file='kblock'//k_char//'_tail.bcs', &
         status="OLD",access="DIRECT",form='binary',recl = nbass*8)
!    open(unit=66, &
!         file='kblock'//k_char//'_tail.bcs', &
!         status="OLD",access="DIRECT",recl = nbass*8)

    do i=1,n_states
       read(66,rec=i+first_state-1) tmp
       block%d((i-1)*nbass+1:i*nbass)=tmp
    enddo
    close(66)
  end subroutine readinwaves

  subroutine elements_triangle(block,low_state,f_state,states,n_op,op_matrix,elements) 
    implicit none
    type(kblock), pointer :: block 
    integer :: low_state,f_state,states
    integer,intent(in) :: n_op
    real(8),intent(in) :: op_matrix(:,:,:,:)
    real(8),intent(inout) :: elements(:,:,:) 
    integer :: l2, l1
    integer :: count,s_state,chunk=1
    s_state=f_state - low_state
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE( l2 ) &
    !$OMP SCHEDULE(DYNAMIC,CHUNK) 

    do l2=1,states
       do l1=1,l2
          elements(l1 + s_state,l2 + s_state,:) = kblock_element(block,l1,n_op,op_matrix,l2)
       end do
    end do

    !$OMP END PARALLEL DO
    return
  end subroutine elements_triangle

  subroutine elements_square( bra_block, ket_block, low_state,&
       bra_states,first_bra,&
       ket_states,first_ket,n_op,op_matrix,elements) 
    implicit none
    type(kblock), pointer :: bra_block ,ket_block 
    integer, intent(in) :: low_state,bra_states,first_bra,ket_states,first_ket
    integer, intent(in) :: n_op
    real(8),intent(in) :: op_matrix(:,:,:,:)
    real(8),intent(inout) :: elements(:,:,:) 
    integer :: l2, l1,chunk=1
    integer :: count,nstates,bra_s_state,ket_s_state
    bra_s_state=first_bra - low_state
    ket_s_state=first_ket - low_state
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE( l1 ) &
    !$OMP SCHEDULE(DYNAMIC,CHUNK) 
    do l1=1,bra_states
       do l2=1,ket_states
          elements(l2 + ket_s_state,l1 + bra_s_state,:) = kblock_element_cross(bra_block,l1,n_op,op_matrix,ket_block,l2)
       end do
    end do
    !$OMP END PARALLEL DO
    return
  end subroutine  elements_square

  function kblock_element(block,bra,n_op,op_matrix,ket) result(element)
    !CALCULATES THE PARTIAL MATRIX ELEMENT FOR A PARTICULAR KBLOCK
    !THE INPUT IS DEFINED AS:
    !
    !BLOCK = VECTOR BLOCK( JK ) OF ALL THE KBLOCKS
    !BRA   = STATE LABEL FOR THE BRA STATE
    !OP_MATRIX = MATRIX WITH DIMENSIONS(NANG,NR1,NR2) WITH OPERATOR CALCULATED AT THE GRID POINTS
    !KET   = STATE LABEL FOR THE KET STATE

    implicit none
    integer,intent(in) :: bra, ket
    type(kblock) :: block
    integer, intent(in) :: n_op
    real(8),intent(in) :: op_matrix(:,:,:,:)
    real(8) :: element(n_op)
    element = kblock_element_cross(block,bra,n_op,op_matrix,block,ket)
  end function kblock_element

  function kblock_element_cross(bra_block,bra,n_op,op_matrix,ket_block,ket) result(element)
    !CALCULATES THE PARTIAL MATRIX ELEMENT FOR A PARTICULAR KBLOCK
    !OF KBLOCKS
    !THE INPUT IS DEFINED AS:
    !
    !BRA_BLOCK = SELF EXPLANATORY
    !KET_BLOCK = SELF EXPLANATORY
    !BRA   = STATE LABEL FOR THE BRA STATE
    ! OP_MATRIX   = MUST BE A REAL FUNCTION OP_MATRIX(R1,R2,COS(THETA)) IN JACOBI COORDINATES
    !KET   = STATE LABEL FOR THE KET STATE

    implicit none
    integer :: bra, ket
    type(kblock),target :: bra_block,ket_block
    integer, intent(in) :: n_op
    real(8),intent(in) :: op_matrix(:,:,:,:)
    real(8) :: element(n_op)
    real(8) :: sum
    real(8) :: op_tmp, op_tmp_t

    integer :: i,j,k,l,count
    real(8), pointer ::r1(:),r2(:),xcos(:)
    real(8), pointer ::ket_d(:),bra_d(:)
    integer, pointer ::nbass,idia,npnt1,npnt2,nang,kpar
    integer :: brapos,ketpos
    integer :: threads,OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,chunk=1,rest
    !the properties following properties are shared by all the k blocks
    !neval jk ipar jrot xmass1 xmass2 xmass3 npnt1 npnt2 r1 r2 eh
!    if (ket_block%nbass.ne.nbass) call error('KBLOCK_ELEMENT_CROSS: ket_block%nbass.ne.nbass')
!    allocate(op(nbass))
    !generate the op_matrix points at the grid:
    idia => bra_block%grid%idia 
    if (idia.ne.ket_block%grid%idia) call error('KBLOCK_ELEMENT_CROSS:bra_block%grid%idia.ne.ket_block%grid%idia.')
    npnt1 => bra_block%grid%npnt1
    npnt2 => bra_block%grid%npnt2
    nang =>  bra_block%grid%nang
    kpar =>  bra_block%kpar
    nbass => bra_block%nbass   
    r1 =>    bra_block%grid%r1
    r2 =>    bra_block%grid%r2
    xcos =>  bra_block%grid%ang !this will have to change in order to be true cross-block
    bra_d =>     bra_block%d
    if (.not.associated(bra_d)) call error('KBLOCK_ELEMENT_CROSS:.not.associated(bra_d)')
    ket_d =>     ket_block%d
    if (.not.associated(ket_d)) call error('KBLOCK_ELEMENT_CROSS:.not.associated(ket_d)')

select case (idia)
case(-2)
do l=1,n_op
    ketpos = (ket-1)*nbass
    brapos = (bra-1)*nbass
    sum=0.0d0
    count=0

    do k=1,npnt1
       do j=1,k-kpar
          do i=1,nang
             op_tmp  =   op_matrix(i,j,k,l)
             op_tmp_t  = op_matrix(i,k,j,l)
             sum = sum + (op_tmp + op_tmp_t)*ket_d( ketpos + i )*bra_d( brapos + i )
          end do
            count=count+nang
!            print*,count
          ketpos = ketpos + nang
          brapos = brapos + nang
       end do
    end do
            
!    if(count.ne.nbass) call error('KBLOCK_ELEMENT_CROSS:count.ne.nbass')
    element(l) =0.5d0*sum    
 end do
case(-1)
do l=1,n_op
    ketpos = (ket-1)*nbass
    brapos = (bra-1)*nbass
    sum=0.0d0
    count=0

    do k=1,npnt2
       do j=1,npnt1
          do i=1,nang
             op_tmp  =   op_matrix(i,j,k,l)
             sum = sum + (op_tmp )*ket_d( ketpos + i )*bra_d( brapos + i )
          end do
            count=count+nang
!            print*,count
          ketpos = ketpos + nang
          brapos = brapos + nang
       end do
    end do
!    if(count.ne.nbass) call error('KBLOCK_ELEMENT_CROSS:count.ne.nbass')
    element(l) = sum    
 end do
case(1)
do l=1,n_op
    ketpos = (ket-1)*nbass
    brapos = (bra-1)*nbass
    sum=0.0d0
    count=0

    do k=1,npnt2
       do j=1,npnt1
          do i=1,nang
             op_tmp  =   op_matrix(i,j,k,l)
             sum = sum + (op_tmp )*ket_d( ketpos + i )*bra_d( brapos + i )
          end do
            count=count+nang
!            print*,count
          ketpos = ketpos + nang
          brapos = brapos + nang
       end do
    end do
!    if(count.ne.nbass) call error('KBLOCK_ELEMENT_CROSS:count.ne.nbass')
    element(l) = sum    
 end do
case(2)
do l=1,n_op
    ketpos = (ket-1)*nbass
    brapos = (bra-1)*nbass
    sum=0.0d0
    count=0

    do k=1,npnt2
       do j=1,npnt1
          do i=1,nang
             !             count=count+1
             op_tmp  =    op_matrix( i , j , k ,l)
             op_tmp_t  =  op_matrix( nang-i+1 , j , k,l )
             sum = sum + ( op_tmp + op_tmp_t ) * ket_d( ketpos + i ) * bra_d( brapos + i )
          end do
          count=count+nang
          ketpos = ketpos + nang
          brapos = brapos + nang
       end do
    end do
    element(l) = 0.5d0*sum    
end do
case default
   call error('KBLOCK_ELEMENT:IDIA value not recognised')
end select

end function kblock_element_cross
  

function matrix_element_from_kblock_set(block,bra,n_op,operator,ket) result(element)
    !WARNING: CURRENT IMPLEMENTATION FOR RADAU BISECTOR COORDINATES ONLY 
    !CALCULATES THE COMPLETE MATRIX ELEMENT FROM THE FULL SET
    !OF KBLOCKS
    !THE INPUT IS DEFINED AS:
    !
    !BLOCK      = VECTOR BLOCK( JK ) OF ALL THE KBLOCKS
    !BRA        = STATE LABEL FOR THE BRA STATE
    !OPERATOR   = MUST BE A REAL FUNCTION OPERATOR(R1,R2,COS(THETA)) IN JACOBI COORDINATES
    !KET        = STATE LABEL FOR THE KET STATE
        implicit none
    integer :: bra, ket
    type(kblock), pointer :: block(:)
    integer, intent(in) :: n_op 
!    real(8), 
external :: operator
    real(8),allocatable :: op_matrix(:,:,:,:)
    real(8) :: element(n_op)
    real(8) :: sum(n_op)
    real(8) :: op_tmp(n_op), op_tmp_t(n_op)
    real(8), pointer :: op(:)
    integer :: i,j,k,l,count
    real(8), pointer :: A(:),B(:),r1,r2,xcos
    integer, pointer :: nbass,neval,idia
!    interface 
!       function operator(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos ; integer :: n ; real(8) :: out(n)
!       end function operator
!    end interface

    !the properties following properties are shared by all the k blocks
    !neval jk ipar jrot xmass1 xmass2 xmass3 npnt1 npnt2 r1 r2 eh
    neval => block(1)%neval
    idia => block(1)%grid%idia !currently only homonuclear diatomics are supported.
                               !abs(idia).eq.-2 homonuclear diatomic Radau coordinates
    sum=0.0d0
    kblockloop: do l=1,block(1)%jk
       allocate( op_matrix( block(l)%grid%nang , block(l)%grid%npnt1 ,&
            block(l)%grid%npnt2 , n_op ) )
       call operator_matrix_build( block(l) , n_op , operator, op_matrix )
       
       sum = sum + kblock_element( block(l) , bra , n_op , op_matrix , ket )
       deallocate(op_matrix)
    end do kblockloop
    element =0.5d0*sum    
  end function matrix_element_from_kblock_set

  function elements_from_kblock_files(low_state,high_state,n_op,operator) result(elements) 
    implicit none
    integer,intent(in) :: low_state,high_state
    integer,intent(in) :: n_op
!    real(8), 
external :: operator
!   interface 
!       function operator(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos ; integer :: n ; real(8) :: out(n)
!       end function operator
!    end interface

    real(8),allocatable :: op_matrix(:,:,:,:)
    real(8),allocatable :: elements(:,:) 
    integer :: l2, l1
    type(kblock), pointer :: block
    integer :: k,kz,kmin,unit,jk,count,nstates
    nstates = high_state-low_state+1
    allocate( elements( nstates*(nstates+1)/2 , n_op ) )
    elements = 0.0d0
    call kblock_head_check(jk,kmin)
    do k = 1,jk
       kz=k-(1-kmin)!note that kmin here is defined as the complementary of kmin in rotlev
       call kblock_read( kz , block )
       allocate( op_matrix( block%grid%nang,block%grid%npnt1 , block%grid%npnt2, n_op ) )
       call operator_matrix_build( block, n_op, operator, op_matrix )
       count=0
       do l1=1,nstates
          print*,k,l1
          do l2=1,l1
             count=count+1
             elements(count,:) = elements(count,:)  + kblock_element(block,l1,n_op,op_matrix,l2)
             !          write(6,"('for bra ',I4,' and ket ',I4,', element is ', d12.3)") l1,l2,element
          end do
       end do
       deallocate(op_matrix)
       call kblock_destroy(block)
    end do
  end function elements_from_kblock_files

!PRIVATE ROUTINES
  subroutine operator_matrix_build(block,n_op,operator,op_matrix)
    implicit none
    type(kblock),intent(in) :: block
    integer, intent(in) :: n_op
!    real(8),
external :: operator
    real(8),intent(out) :: op_matrix(:,:,:,:)
    real(8), pointer::tmp_op(:)
    real(8), pointer :: xcos(:),r1(:),r2(:)
    integer, pointer :: nbass
    integer :: chunk=1,i,j,k

!    interface 
!       function operator(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos ; integer :: n ; real(8) :: out(n)
!       end function operator
!    end interface

    nbass => block%nbass   
    r1 => block%grid%r1
    r2 => block%grid%r2
    xcos => block%grid%ang
    !generate the op_matrix points at the grid:
    
    if(block%jrot == 0 ) then
          do k=1,block%grid%npnt2
             do j=1,block%grid%npnt1
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(SHARED) PRIVATE( i, tmp_op ) &
                !$OMP SCHEDULE(DYNAMIC,CHUNK)
                do i=1,block%grid%nang
                   !selects the angular grid points included in the 3D Hamiltonian
                   if(block%iv(i) /= 0) then
                      
    allocate(tmp_op(n_op))
!                      op_matrix(i,j,k,:)  = operator(r1(j),r2(k),xcos(i),n_op)
                      call operator(r1(j),r2(k),xcos(i),n_op,tmp_op)
                      op_matrix(i,j,k,:) = tmp_op
    deallocate(tmp_op)
                   else
                      op_matrix(i,j,k,:) = 0.0d0
                   end if
                end do
                !$OMP END PARALLEL DO
             end do
          end do
    else
          do k=1,block%grid%npnt2
             do j=1,block%grid%npnt1
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(SHARED) PRIVATE( i, tmp_op ) &
                !$OMP SCHEDULE(DYNAMIC,CHUNK)
                do i=1,block%grid%nang
                   
 !                  op_matrix(i,j,k,:)  = operator(r1(j),r2(k),xcos(i),n_op)
                   allocate(tmp_op(n_op))
                  call operator(r1(j),r2(k),xcos(i),n_op,tmp_op)
                  op_matrix(i,j,k,:) = tmp_op
                  deallocate(tmp_op)
                end do
                !$OMP END PARALLEL DO
             end do
          end do
    end if
    
    
  end subroutine operator_matrix_build


  function dot_product(length,A,B) result(C)
    integer :: length
    real(8) ::A(length),B(length),C
    C=0.0d0
    ! BLAS ROUTINE dgemm
    ! this subroutine does matrix multiplication and addition C = a*( A * B )+b*( C )
    ! where a and b are scalars and A, B and C are matrices 
    !            _ transposes the first vector (A) to become a row
    !           |    _ leaves the second vector (B) as a column 
    !           |   |   _M: number of rows of A   
    !           |   |  |  _N: number of columns B    
    !           |   |  | |  -K: number of columns of A (rows of B) 
    !           |   |  | | |        _a: scalar coefficient
    !           |   |  | | |      |      _A: matrix
    !           |   |  | | |      |     |  _LDA: leading dimension of A
    !           |   |  | | |      |     | |       _B: matrix
    !           |   |  | | |      |     | |      |  _LDB: leading dimension of A
    !           |   |  | | |      |     | |      | |       _b: scalar coefficient
    !           |   |  | | |      |     | |      | |      |      _C: matrix
    !           |   |  | | |      |     | |      | |      |     |  _leading dimension of C
    !           |   |  | | |      |     | |      | |      |     | |
    call dgemm('t','n',1,1,length,1.0d0,A,length,B,length,0.0d0,C,1)
    !this way it is setup to give the dot product between vectors A and B
  end function dot_product

end module class_kblock_to_elements

!program normcheck
!use class_kblock
!use class_kblock_to_elements
!implicit none
!type(kblock),pointer :: kblck
!real(8),allocatable :: elements(:)
!real(8),external :: unity  
!integer :: i,j, count,neval,start_state,end_state,nstates
!
!
!call kblock_head_read(0,kblck)
!neval=kblck%neval
!call kblock_destroy(kblck)
!
!start_state=neval - 10
!end_state=neval
!
!nstates=-(start_state-end_state)+1
!allocate( elements( nstates*(nstates+1)/2 ) )
!elements = elements_from_kblock_files_mem_split(1000.0d0,start_state,end_state,unity) 
!count=0
!do i=start_state,end_state
!do j=start_state,i
!count = count +1
!if(i.eq.j)write(6,*) elements(count) 
!end do
!end do
!
!end program normcheck
!
!function unity(r1,r2,xcos) result(out)
!  implicit none
!  real(8) :: r1,r2,xcos,out
!  out=1.0d0
!end function unity


