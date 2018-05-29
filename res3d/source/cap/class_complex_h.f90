module class_complex_h
  use class_kblock_to_elements
  use class_kblock
  use class_geom
  use class_input
  use common_parameters ! global variables
  use class_error
  implicit none
  private
  type, public:: cap_matrix
     integer, pointer :: low_state
     integer, pointer :: high_state
     real(8), pointer :: eh1 ! evergy of the ground state
     real(8), pointer :: eh(:)
     real(8), pointer :: w(:) ! triangular matrix
  end type cap_matrix
  public :: cap_matrix_build
  !subroutine cap_matrix_build(in)
  !  type(input) ::  in
  !end subroutine cap_matrix_build
  public :: cap_matrix_read
  !subroutine cap_matrix_read(unit,cap)
  !  type(cap_matrix), pointer :: cap
  !  integer :: unit
  !end subroutine cap_matrix_read
  public :: cap_matrix_write
  !subroutine cap_matrix_write(unit,cap)
  !  integer :: unit
  !  type(cap_matrix) :: cap
  !end subroutine cap_matrix_write
  public :: cap_matrix_destroy
  !subroutine cap_matrix_destroy(cap)
  !  type(cap_matrix), pointer :: cap
  !end subroutine cap_matrix_destroy
  public :: complex_h_diagonalise
  !subroutine complex_h_diagonalise(diagunit,diagvecunit,in,cap)
  !  integer, intent(in) :: diagunit
  !  integer, intent(in) :: diagvecunit
  !  type(diag_input), intent(in) :: in
  !  type(cap_matrix), intent(in) :: cap
  !end subroutine cap_matrix_diag
contains

  subroutine cap_matrix_build(in)
    implicit none
    type(input) ::  in
    type(cap_matrix),pointer :: cap
    type(kblock),pointer :: block
    character(len=200) :: filename
    real(8) :: tmp(3)
    integer :: jk, kz,kzth,k

   filename = in%cap%wavefile
    write(6,*)
    write(6,"(20x,'----- BUILDING CAP MATRIX -----')")
    write(6,*)

       !read wavefunction and extract k blocks
       !from class_wavefile
       call wavefile_to_kblock_files(filename)
       !read kblocks and dimension the cap matrix
       !also define parameters from input
       !       call cap_matrix_set_size_from_kblock_files(in%cap,cap)
       tmp=-100000000.0d0
       call kblock_head_check(jk,kz)
       call kblock_head_dump(in,kz)

       do k = 1,jk
          kzth = kz+k-1
          call kblock_head_read(kzth,block)
          call common_parameters_from_kblock_file(in,block)
          call cap_matrix_setup_from_kblock_file(block,cap)
          call kblock_head_destroy(block)
          write(6,*)
          write(6,"(10x,'Setting up rmax for all dissociation channels...')")

          if(rmax(1).gt.tmp(1)) tmp(1)=rmax(1)
          if(rmax(2).gt.tmp(2)) tmp(2)=rmax(2)
          if(rmax(3).gt.tmp(3)) tmp(3)=rmax(3)

       end do

       select case(in%cap%rmaxchoice)
       case('long')
          write(6,"(10x,'Longest CAP rmax chosen for all channels:')")
          rmax=maxval(rmax,1)
       case('short')
          write(6,"(10x,'Shortest CAP rmax chosen for all channels:')")
          rmax=minval(rmax,1)
       case('adapt')
          write(6,"(10x,'Adapted CAP rmax chosen for all channels:')")
       case default
          call error('COMMON_PARAMETERS_FROM_WAVEFILE:unknown rmaxchoice -- aborting')
       end select
       write(6,"(10x,'rmax(1) ',f5.2,' au')")rmax(1)
       write(6,"(10x,'rmax(2) ',f5.2,' au')")rmax(2)
       write(6,"(10x,'rmax(3) ',f5.2,' au')")rmax(3)

    coords: select case(in%cap%coord)
    case("123")
       cap3: select case(in%cap%capchoice)
       case("simple")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple_cap_123)
       case("simple8")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple8_cap_123)
       case("simple3")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple3_cap_123)
       case("manop")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,manop_cap_123)
       case default
          print*,'unrecognised cap... stoping at subroutine cap'
          stop
       end select cap3
    case("12-3")
       cap2: select case(in%cap%capchoice)
       case("simple")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple_cap_12_3)
       case("simple8")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple8_cap_12_3)
       case("manop")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,manop_cap_12_3)
       end select cap2
    case("1-23")
       cap11: select case(in%cap%capchoice)
       case("simple")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple_cap_1_23)
       case("simple8")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple8_cap_1_23)
       case("simple3")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple3_cap_1_23)
       case("manop")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,manop_cap_1_23)
       end select cap11
    case("2-31")
       cap12: select case(in%cap%capchoice)
       case("simple")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple_cap_2_31)
       case("simple8")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple8_cap_2_31)
       case("simple3")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple3_cap_2_31)
       case("manop")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,manop_cap_2_31)
       end select cap12
    case("3-12")
       cap13: select case(in%cap%capchoice)
       case("simple")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple_cap_3_12)
       case("simple8")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple8_cap_3_12)
       case("simple3")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,simple3_cap_3_12)
       case("manop")
             call cap_matrix_deltar_iterate(in%cap%mem,filename,cap,manop_cap_3_12)
       end select cap13
    case default
       call error("CAP_MATRIX_BUILD: unrecognised geometry")
    end select coords
  end subroutine cap_matrix_build

  subroutine cap_matrix_setup_from_kblock_file(block,cap)
    implicit none
    type(input) :: in
    type(kblock), pointer :: block
    type(cap_matrix), pointer :: cap
    allocate(cap)
    allocate(cap%low_state)
    allocate(cap%high_state)
    cap%low_state=cap_low_state
    cap%high_state=cap_high_state
    allocate(cap%eh1)
    cap%eh1=gstate
    allocate(cap%w(n_cap_states*(n_cap_states+1)/2))
    allocate(cap%eh(n_cap_states))
    cap%eh=block%eh(cap_low_state:cap_high_state)
    RETURN
  end subroutine cap_matrix_setup_from_kblock_file

  !private method
  subroutine cap_matrix_deltar_iterate(mem,wavefilename,cap,cap_function)
    implicit none
    real(8),intent(in) :: mem
    character(len=200) :: wavefilename
    type(cap_matrix),pointer :: cap
!    real(8),
    external :: cap_function
    character(len=200) :: filename
    integer :: i
    character(len=2) :: capunit
    logical :: exists
    integer :: deltar_count
    integer :: count
    real(8),allocatable::w(:,:)
!    interface
!       function cap_function(r1,r2,xcos,n) result(out)
!         real(8) :: r1,r2,xcos
!         integer :: n
!         real(8) :: out(n)
!       end function cap_function
!    end interface

    deltar_count=0
    do i=1,ndeltar
       dr => deltar(i)
       write(capunit,fmt="(I2.2)") i
       filename='cap'//capunit//'.bcs'
       inquire(file=trim(filename),exist=exists)
       if(.not.exists) then
          write(6,"(10x,'building cap matrix with deltar= ',F5.2,' au ...')") dr

          deltar_choice(i)=.true.
          deltar_count = deltar_count + 1
       else
          write(6,*)'Warning: file ', trim(filename),' already exists -- skipping'

          deltar_choice(i)=.false.

       end if
    end do

    allocate ( w( n_cap_states * (n_cap_states+1) / 2 , deltar_count) )
    allocate ( deltar_small(deltar_count) )

    count = 0
    do i=1,ndeltar
       if(deltar_choice(i)) then
          count=count+1 ; deltar_small(count)=deltar(i)
       end if
    end do

    if (deltar_count > 0 ) &
         w = elements_from_kblock_files_mem_split( mem , &
         CAP%LOW_STATE , CAP%HIGH_STATE, deltar_count , cap_function )

    count=0
    do i=1,ndeltar
       write(capunit,fmt="(I2.2)") i
       filename='cap'//capunit//'.bcs'

       if(deltar_choice(i)) then
          count=count+1 ; CAP%W = w( : , count )

          open(unit=101,file=trim( filename ),form='unformatted',status='new')
          call cap_matrix_write( 101 , cap )
          close(unit=101)
       end if
    end do
  end subroutine cap_matrix_deltar_iterate

  subroutine cap_matrix_read(unit,cap)
    type(cap_matrix), pointer :: cap
    integer :: unit, nstates
    allocate(cap)
    allocate(cap%low_state)
    read(unit) cap%low_state
    allocate(cap%high_state)
    read(unit) cap%high_state
    nstates=cap%high_state-cap%low_state+1
    allocate(cap%eh1)
    read(unit) cap%eh1
    allocate(cap%eh(nstates))
    read(unit) cap%eh
    allocate(cap%w(nstates*(nstates+1)/2)) ! triangular matrix
    read(unit) cap%w
  end subroutine cap_matrix_read

  subroutine cap_matrix_write(unit,cap)
    type(cap_matrix) :: cap
    integer :: unit
    write(unit) cap%low_state
    write(unit) cap%high_state
    write(unit) cap%eh1
    write(unit) cap%eh
    write(unit) cap%w
  end subroutine cap_matrix_write

  subroutine cap_matrix_destroy(cap)
    type(cap_matrix), pointer :: cap
    deallocate(cap%low_state)
    deallocate(cap%high_state)
    deallocate(cap%eh1)
    deallocate(cap%eh)
    deallocate(cap%w) ! triangular matrix
    deallocate(cap)
  end subroutine cap_matrix_destroy


!
  subroutine complex_h_diagonalise(diagunit,diagvecunit,in,cap)
    implicit none
    integer, intent(in) :: diagunit
    integer, intent(in) :: diagvecunit
    type(diag_input), intent(in) :: in
    type(cap_matrix), intent(in) :: cap
    ! Adapted by Bruno Silva to read from new CAP generation codes
    integer :: m,mi,mf,j,j1,j2,i,i1,i2,info,m0,ldvl,ic
    real(8) ::  ene,thre,emin,sum,lambda
    real(8),allocatable ::  ene_vec(:)
    integer :: mr
    integer ::lwork,md
    real(8), parameter :: conv=219474.63067d0
    !conv -- coversion factor from Eh to cm-1
    integer   ,allocatable  :: ind(:),indp(:),ind2(:),ind3(:)
    real(8)   ,allocatable  ::  eig0(:),w(:,:),dist(:)
    complex(8),allocatable  :: oldh(:,:),origh(:,:),orightmp(:,:)
    complex(8),allocatable  :: eigo(:)

    real(8),pointer  :: rwork(:)
    complex(8),pointer :: ham(:,:),eigl(:),vr(:,:),work(:)
    !    !$OMP threadprivate (rwork,ham,eigl,vr,work)
    integer :: n_vec_mem, b, b_final,mb,b_start, b_end, mem_steps
    complex(8),allocatable :: eig_buffer(:,:), vec_buffer(:,:,:)

    complex(8) :: ci,vl,vnor
    complex(8),allocatable :: vl_vec(:),vnor_vec(:)
    character(len=200) :: filename
    character(len=4) :: char_i
    integer :: count,count2,CHUNK=1,thread,NTHREADS,LRESTART
    logical :: bool_test, exists
    INTERFACE
       INTEGER FUNCTION omp_get_thread_num ()
       END FUNCTION omp_get_thread_num
    END INTERFACE
    INTERFACE
       INTEGER FUNCTION omp_get_num_threads ()
       END FUNCTION omp_get_num_threads
    END INTERFACE

    mr = n_diag_states!diag_high_state - diag_low_state + 1 !number of states used in the calculation
    lwork = 10 * mr    ! parameter used in the diagonaliser
    !md=100          ! maximum of possible contact points between different curves (different lambdas)
    !thre=0.9d0     ! Threshold for comparing different eigenvalues
    md =   in%ccnumber       ! maximum of possible contact points between different curves (different lambdas)
    thre = in%ecthres        ! Threshold for comparing different eigenvalues
    lambda=0.0d0   ! initialise lambda
    !    m  = d%cap%last_state
    !    m0 = d%cap%first_state
    !    mi = d%diag%first_state
    !    mf = d%diag%last_state
    m  = cap_high_state
    m0 = cap_low_state
    mi = diag_low_state
    mf = diag_high_state
    !n_cap_states=m-m0+1
    ci=(0.d0,1.d0) ! define complex number "i".
    ldvl=1         !parameter for the diagonaliser
    allocate( ind(mr),indp(md),ind2(mr),ind3(mr)            )
    allocate( eig0(mr),w(mr,mr),dist(md)        )
    allocate( oldh(mr,mr),origh(mr,mr),orightmp(mr,mr)   )
    allocate( eigo(mr)                            )

    ind  = 0
    indp = 0
    ind2 = 0
    ind3 = 0
    eig0 = 0.0d0
    w    = 0.0d0
    dist = 0.0d0
    oldh = (0.0d0,0.0d0)
    origh = (0.0d0,0.0d0)
    eigo = (0.0d0,0.0d0)
    write(6,*)'         Lambda range (converted to cm-1 from Eh) ....'
    write(6,*)'         l(2). ..........',in%h*conv
    write(6,*)'         N ..............',in%nlambda
    write(6,*)'         l(NLAMBDA) .....',in%h*conv*(in%chi**(in%nlambda-1)-1.d0)/(in%chi-1.d0)

    !READ THE RELEVANT EIGENVALUES:
    j=0
    do i=mi-m0+1,mf-m0+1
       j=j+1
       eig0(j)=cap%eh(i)
    end do

    ! READ THE RELEVANT PART OF THE CAP MATRIX:
    count=0
    count2=0
    j1=0
    do i1=m0,m
       j2=0
       if ((i1).ge.mi.and.(i1).le.mf) then
          j1=j1+1
          bool_test=.true.
       else
          bool_test=.false.
       end if
       do i2=m0,i1
          count=count+1
          if ( bool_test &
               .and. (i2).ge.mi.and.(i2).le.mf ) then
             count2=count2+1
             j2=j2+1
             w(j1,j2)=cap%w(count)
             w(j2,j1)=cap%w(count)
          end if
       end do
    end do
    write(6,*)'         W-matrix read ',count2 ,'triangle elements.... '


    NTHREADS=OMP_GET_NUM_THREADS()
    !CRUDE RESTART FACILITY
    !    call system('mkdir tmp_data')
    !    restart: do i = in%nlambda,1,-1
    !    write(char_i,"(I4.4)") i
    !       inquire(file='tmp_data/tmpdiagvecs'//char_i//'.bcs',exist=exists)
    !       if (exists.and.i.gt.nthreads) then
    !          lrestart=i-nthreads
    !          exit restart
    !       else
    !          lrestart=1
    !       end if
    !    end do restart


    N_VEC_MEM = ( INT(MEM) / 16 )* 1024 / N_DIAG_STATES * 1024  / N_DIAG_STATES
    if ( N_VEC_MEM == 0 ) call error('subroutine cap_matrix_diag: N_VEC_MEM == 0 ')
    MEM_STEPS = 1
    IF (IN%NLAMBDA > N_VEC_MEM) THEN
       MEM_STEPS = IN%NLAMBDA / N_VEC_MEM + 1
       B = N_VEC_MEM
       B_FINAL = MOD ( IN%NLAMBDA ,N_VEC_MEM )
    ELSE IF  (IN%NLAMBDA <= N_VEC_MEM) THEN
       B_FINAL = IN%NLAMBDA
       B = B_FINAL
       MEM_STEPS = 1
    END IF
    ALLOCATE (EIG_BUFFER(N_DIAG_STATES,B))
    ALLOCATE (VEC_BUFFER(N_DIAG_STATES,N_DIAG_STATES,B))
    !MEMORY BLOCK LOOP
    MEM_LOOP: DO MB=1,MEM_STEPS
       eig_buffer = 0.0d0
       vec_buffer = 0.0d0
       B_START = ( MB - 1 ) * B  + 1
       B_END   =  MB * B
       if ( MB == MEM_STEPS ) B_END = ( MB - 1 ) * B + B_FINAL

       !DIAGONALISATION LOOP:

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(SHARED) PRIVATE( i,thread,char_i,lambda,i1,i2,info,vl,ham,eigl,vr,work,rwork ) &
       !$OMP SCHEDULE(DYNAMIC,CHUNK)
       diag:do i= B_START  , B_END

          !print*,'b'

          ! arrays that have to be allocated separately (OMP):
          allocate( ham(mr,mr),eigl(mr),vr(mr,mr),work(lwork),rwork(2*mr) )
          rwork= 0.0d0
          ham  = (0.0d0,0.0d0)
          eigl = (0.0d0,0.0d0)
          work = (0.0d0,0.0d0)
          vr   = (0.0d0,0.0d0)
          !build the complex non-hermitian hamiltonian part
          lambda = in%h * ( in%chi ** (i-1) - 1.d0 ) / ( in%chi - 1.d0 )
          do i1=1,mr                                ! Build the hamiltonian matrix
             do i2=1,mr                             ! to be diagonalised using the
                ham(i1,i2)=-ci*lambda*w(i1,i2)      ! cap matrix and adding to it
             end do! i2=1,mr                        ! the eigenvalues from the DVR3D
             ham(i1,i1)=eig0(i1)+ham(i1,i1)         ! calculation.
          end do! i1=1,mr

          CALL ZGEEV('N','V',mr,ham,mr,eigl, VL, LDVL,vr,mr,&  !Call diagonaliser
               WORK, LWORK, RWORK, INFO )                      !
          if(debug)   write(6,*)'INFO....',info,'   i = ',i
          !thread = omp_get_thread_num()
          !write(char_i,"(I4.4)") i
          !open(unit = 10+thread , file = 'tmp_data/tmpdiagvecs'//char_i//'.bcs', access='direct',recl=8*n_diag_states*n_diag_states)
          !write(unit=10+thread,rec=1) vr
          !close(10+thread)
          !open(unit = 10+thread , file = 'tmp_data/tmpdiageigs'//char_i//'.bcs', access='direct',recl=8*n_diag_states)
          !write(unit=10+thread,rec=1) eigl
          !close(10+thread)
          !print*,'c'
          vec_buffer(:,:,i - B_start + 1 ) = vr(:,:)
          eig_buffer(:,i - B_start + 1 ) = eigl(:)
          deallocate( ham,eigl,vr,work,rwork )
          !print*,'d'
       end do diag
       !$OMP END PARALLEL DO

       allocate( eigl(mr),vr(mr,mr) )

       !sort out the eigenvectors
       vecmatch:do i=B_start,B_end   !(outer loop)

          !SETUP THE EIGEVECTOR RESORTING
          if ( i == 1 ) then
             eigl = (0.0d0,0.0d0)
             vr   = (0.0d0,0.0d0)

             !write(char_i,"(I4.4)") 1
             !open(unit = 300 , file = 'tmp_data/tmpdiagvecs'//char_i//'.bcs', access='direct',recl=8*n_diag_states*n_diag_states)
             !read(unit=300,rec=1) vr
             !close(300)
             !open(unit = 300 , file = 'tmp_data/tmpdiageigs'//char_i//'.bcs', access='direct',recl=8*n_diag_states)
             !read(unit=300,rec=1) eigl
             !close(300)
             !call system('rm -f tmp_data/tmpdiagvecs'//char_i//'.bcs')
             !call system('rm -f tmp_data/tmpdiageigs'//char_i//'.bcs')
             vr(:,:) = vec_buffer(:,:,1)
             eigl(:) = eig_buffer(:,1)

             do i1=1,mr                         ! this loop normalises the eigenvectors
                vl=(0.d0,0.d0)                  ! from each diagonalisation.
                do j=1,mr                       !
                   vl=vl+vr(j,i1)*vr(j,i1)      !
                end do              ! j=1,mr    !
                vnor=1.d0/sqrt(vl)              !
                do j=1,mr                       !
                   origh(j,i1)=vr(j,i1)
                   vr(j,i1)=vr(j,i1)*vnor       !
                end do              ! j=1,mr    !
             end do                 ! i1=1,mr   !
             ! UPDATE THE OLD VECTOR AND WRITE TO FILE THE EIGENVALUES
             do i1=1,mr
                do i2=1,mr
                   oldh(i1,i2)=vr(i1,i2)
                end do! i2=1,mr
                write(unit=diagunit)lambda,i1,1,eigl(i1)*conv    !eigenvalues in cm-1
                eigo(i1)=eigl(i1)
             end do! i1=1,mr
             ! WRITE EIGENVECTORS TO FILE IF REQUESTED
             !Commented out by L Lodi and BC Silva 7- May 2009
!              if (in%outputvec) then
!                 write(unit=diagvecunit,rec=1) vr!(:,:)
!              end if

          else

             !write(char_i,"(I4.4)") i
             !open(unit = 300 , file = 'tmp_data/tmpdiagvecs'//char_i//'.bcs', access='direct',recl=8*n_diag_states*n_diag_states)
             !read(unit=300,rec=1) vr
             !close(300)
             !open(unit = 300 , file = 'tmp_data/tmpdiageigs'//char_i//'.bcs', access='direct',recl=8*n_diag_states)
             !read(unit=300,rec=1) eigl
             !close(300)
             !call system('rm -f tmp_data/tmpdiagvecs'//char_i//'.bcs')
             !call system('rm -f tmp_data/tmpdiageigs'//char_i//'.bcs')
             vr(:,:) = vec_buffer(:,:,i - B_start + 1 )
             eigl(:) = eig_buffer(:,i - B_start + 1 )

             allocate(vl_vec(mr))
             allocate(vnor_vec(mr))
             !$OMP PARALLEL DO &
             !$OMP DEFAULT(SHARED) PRIVATE( i1,j ) &
             !$OMP SCHEDULE(DYNAMIC,CHUNK)
             do i1=1,mr                         ! this loop normalises the eigenvectors
                vl_vec(i1)=(0.d0,0.d0)                  ! from each diagonalisation.
                do j=1,mr                       !
                   vl_vec(i1)=vl_vec(i1)+vr(j,i1)*vr(j,i1)      !
                end do
                vnor_vec(i1)=1.d0/sqrt(vl_vec(i1))              !
                do j=1,mr
                   orightmp(j,i1)=vr(j,i1)                       !
                   vr(j,i1)=vr(j,i1)*vnor_vec(i1)       !
                end do
             end do
             !$OMP END PARALLEL DO
             deallocate(vl_vec)
             deallocate(vnor_vec)
             !resolve the eigenvector correlation between lambda steps
             resolve:do i1=1,mr

                do ic=1,md
                   indp(ic)=0
                   dist(ic)=1.d100
                end do


                ic=0
                ind(i1)=0
                !parelellise loop
                allocate(vl_vec(mr))
                allocate(ene_vec(mr))
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(SHARED) PRIVATE( i2,j ) &
                !$OMP SCHEDULE(DYNAMIC,CHUNK)
                do i2=1,mr
                   vl_vec(i2)=(0.d0,0.d0)
                   do j=1,mr
                      vl_vec(i2)=vl_vec(i2)+(oldh(j,i1))*vr(j,i2)
                   end do
                   ene_vec(i2)=sqrt(dble(conjg(vl_vec(i2))*vl_vec(i2)))
                end do
                !$OMP END PARALLEL DO
                deallocate(vl_vec)


                do i2=1,mr
                   if (ene_vec(i2).gt.thre) then
                      ind(i1)=i2
                      ic=ic+1
                      vl=eigo(i2)-eigl(i1)
                      indp(ic)=i2
                      dist(ic)=dsqrt(dble(conjg(vl)*vl))
                   end if
                end do
                deallocate(ene_vec)

                if (ic.eq.0) then
                   write(6,*)'Warning: Missing contact in '
                   write(6,*)'         eig ...',i1
                   write(6,*)'         lmb ...',i
                else if (ic.gt.1) then
                   write(6,*)'Warning: Multiple contact in '
                   write(6,*)'         eig ........',i1
                   write(6,*)'         lmb ........',i
                   write(6,*)'         Contacts ...',ic
                   write(6,*)'         Resolving multiple contact '
                   write(6,*)(dist(j),j=1,ic)
                   j1=indp(1)
                   do j=2,ic
                      if (dabs(dist(j)).gt.dabs(dist(j-1))) j1=indp(j)
                   end do
                   ind(i1)=j1
                end if

                emin=1.d100
                do i2=1,mr
                   vl=eigo(i1)-eigl(i2)
                   ene=dsqrt(dble(conjg(vl)*vl))
                   if (ene.lt.emin) then
                      ind2(i1)=i2
                      emin=ene
                   end if
                end do


             end do resolve ! do i1=1,mr

             j1=0
             do j=1,mr
                j1=j1+ind(j)
             end do

             j2=mr*(mr+1)/2

             do i1=1,mr
                do i2=1,mr
                   if(ind(i1).eq.0)then
                      oldh(i2,i1)=vr(i2,i1)
                   else
                      oldh(i2,i1)=vr(i2,ind(i1))
                      origh(i2,i1)=orightmp(i2,ind(i1))
                   end if
                end do

                if(ind(i1).eq.0)then
                   write(unit=diagunit)lambda,i1,i,eigl(i1)*conv
                   !             write(unit=diagunit,fmt=*)lambda,i1,i,eigl(i1)*conv
                   eigo(i1)=eigl(i1)
                else
                   !write(87,100)i1,i,eigl(ind(i1))*conv
                   write(unit=diagunit)lambda,i1,i,eigl(ind(i1))*conv
                   !             write(unit=diagunit,fmt=*)lambda,i1,i,eigl(ind(i1))*conv
                   eigo(i1)=eigl(ind(i1))
                end if
             end do

             ! WRITE EIGENVECTORS TO FILE IF REQUESTED
             if (in%outputvec) then
                write(unit=diagvecunit,rec=i) oldh!(:,:)
             end if
          end if
       end do vecmatch ! do i=2,in%nlambda
       deallocate( eigl,vr )
    END DO MEM_LOOP
100 format(2I5,d20.13,2x,d20.13)
102 format(2I5,f12.6)
103 format(3I5,d25.18)
108 format(3I5,d25.18,2x,d25.18)
104 format(' ----------  ',d25.18)
    print*
    print*, '         diagonalisation complete'
    deallocate( ind,indp,ind2,ind3)
    deallocate( eig0,w,dist)
    deallocate( oldh )
    deallocate( eigo )
    call system('rm -rf tmp_data')
  end subroutine complex_h_diagonalise

!  subroutine cap_matrix_diag_ser(diagunit,diagvecunit,in,cap)
!    implicit none
!    integer, intent(in) :: diagunit
!    integer, intent(in) :: diagvecunit
!    type(diag_input), intent(in) :: in
!    type(cap_matrix), intent(in) :: cap
!    ! Adapted by Bruno Silva to read from new CAP generation codes
!    integer :: m,mi,mf,j,j1,j2,i,i1,i2,info,m0,ldvl,ic
!    real(8) ::  ene,thre,emin,sum,lambda
!    integer :: mr
!    integer ::lwork,md
!    real(8), parameter :: conv=219474.63067d0
!    !conv -- coversion factor from Eh to cm-1
!    integer   ,allocatable  :: ind(:),indp(:),ind2(:),ind3(:)
!    real(8)   ,allocatable  ::  eig0(:),w(:,:),rwork(:),dist(:)
!    complex(8),allocatable  :: ham(:,:),oldh(:,:),origh(:,:),orightmp(:,:),eigl(:),work(:)
!    complex(8),allocatable  :: vr(:,:),eigo(:)
!    complex(8) :: ci,vl,vnor
!    character(len=200) :: filename
!    integer :: count,count2
!    logical :: bool_test
!    mr = n_diag_states!diag_high_state - diag_low_state + 1 !number of states used in the calculation
!    lwork = 10 * mr    ! parameter used in the diagonaliser
!    !md=100          ! maximum of possible contact points between different curves (different lambdas)
!    !thre=0.9d0     ! Threshold for comparing different eigenvalues
!    md =   in%ccnumber       ! maximum of possible contact points between different curves (different lambdas)
!    thre = in%ecthres        ! Threshold for comparing different eigenvalues
!    lambda=0.0d0   ! initialise lambda
!    !    m  = d%cap%last_state
!    !    m0 = d%cap%first_state
!    !    mi = d%diag%first_state
!    !    mf = d%diag%last_state
!    m  = cap_high_state
!    m0 = cap_low_state
!    mi = diag_low_state
!    mf = diag_high_state
!    !n_cap_states=m-m0+1
!    ci=(0.d0,1.d0) ! define complex number "i".
!    ldvl=1         !parameter for the diagonaliser
!    allocate( ind(mr),indp(md),ind2(mr),ind3(mr)            )
!    allocate( eig0(mr),w(mr,mr),rwork(2*mr),dist(md)        )
!    allocate( ham(mr,mr),oldh(mr,mr),origh(mr,mr),orightmp(mr,mr),eigl(mr),work(lwork)   )
!    allocate( vr(mr,mr),eigo(mr)                            )
!    ind  = 0
!    indp = 0
!    ind2 = 0
!    ind3 = 0
!    eig0 = 0.0d0
!    w    = 0.0d0
!    rwork= 0.0d0
!    dist = 0.0d0
!    ham  = (0.0d0,0.0d0)
!    oldh = (0.0d0,0.0d0)
!    eigl = (0.0d0,0.0d0)
!    work = (0.0d0,0.0d0)
!    vr   = (0.0d0,0.0d0)
!    eigo = (0.0d0,0.0d0)
!    write(6,*)'         Lambda range (converted to cm-1 from Eh) ....'
!    write(6,*)'         l(2). ..........',in%h*conv
!    write(6,*)'         N ..............',in%nlambda
!    write(6,*)'         l(NLAMBDA) .....',in%h*conv*(in%chi**(in%nlambda-1)-1.d0)/(in%chi-1.d0)
!
!    !READ THE RELEVANT EIGENVALUES:
!    j=0
!    do i=mi-m0+1,mf-m0+1
!       j=j+1
!       eig0(j)=cap%eh(i)
!    end do
!
!    ! READ THE RELEVANT PART OF THE CAP MATRIX:
!    count=0
!    count2=0
!    j1=0
!    do i1=m0,m
!       j2=0
!       if ((i1).ge.mi.and.(i1).le.mf) then
!          j1=j1+1
!          bool_test=.true.
!       else
!          bool_test=.false.
!       end if
!       do i2=m0,i1
!          count=count+1
!          if ( bool_test &
!               .and. (i2).ge.mi.and.(i2).le.mf ) then
!             count2=count2+1
!             j2=j2+1
!             w(j1,j2)=cap%w(count)
!             w(j2,j1)=cap%w(count)
!          end if
!       end do
!    end do
!    write(6,*)'         W-matrix read ',count2 ,'triangle elements.... '
!    lambda=0.0d0
!    do i1=1,mr
!       do i2=1,mr
!          ham(i1,i2)=-ci*lambda*w(i1,i2)
!       end do
!       ham(i1,i1)=eig0(i1)+ham(i1,i1)
!       eigo(i1)=(0.d0,0.d0)
!       eigl(i1)=(0.d0,0.d0)
!    end do
!    ! FIRST DIAGONALIZER CALL TO GET THE REFERENCE EIGENVECTORS
!    !      CALL ZHEEV('NLAMBDA','U',mr,ham,mr,eigl,work,lwork,rwork,INFO )
!    CALL ZGEEV('N','V',mr,ham,mr,eigl, VL, LDVL,vr,mr,&
!         WORK, LWORK, RWORK, INFO )
!    if(debug)  write(6,*)'INFO....',info
!
!    do i1=1,mr                         ! this loop normalises the eigenvectors
!       vl=(0.d0,0.d0)                  ! from each diagonalisation.
!       do j=1,mr                       !
!          vl=vl+vr(j,i1)*vr(j,i1)      !
!       end do              ! j=1,mr    !
!       vnor=1.d0/sqrt(vl)              !
!       do j=1,mr                       !
!          origh(j,i1)=vr(j,i1)
!          vr(j,i1)=vr(j,i1)*vnor       !
!
!       end do              ! j=1,mr    !
!    end do                 ! i1=1,mr   !
!    ! UPDATE THE OLD VECTOR AND WRITE TO FILE THE EIGENVALUES
!    do i1=1,mr
!       do i2=1,mr
!          oldh(i1,i2)=vr(i1,i2)
!       end do! i2=1,mr
!       write(unit=diagunit)lambda,i1,1,eigl(i1)*conv    !eigenvalues in cm-1
!       eigo(i1)=eigl(i1)
!    end do! i1=1,mr
!    ! WRITE EIGENVECTORS TO FILE IF REQUESTED
!    if (in%outputvec) then
!       write(unit=diagvecunit,rec=1) origh!(:,:)
!    end if
!    ! REMAINING DIAGONALISATIONS
!    do i=2,in%nlambda   !(outer loop)
!       lambda = in%h * ( in%chi ** (i-1) - 1.d0 ) / ( in%chi - 1.d0 )
!       do i1=1,mr                                ! Build the hamiltonian matrix
!          do i2=1,mr                             ! to be diagonalised using the
!             ham(i1,i2)=-ci*lambda*w(i1,i2)      ! cap matrix and adding to it
!          end do! i2=1,mr                        ! the eigenvalues from the DVR3D
!          ham(i1,i1)=eig0(i1)+ham(i1,i1)         ! calculation.
!       end do! i1=1,mr
!       CALL ZGEEV('N','V',mr,ham,mr,eigl, VL, LDVL,vr,mr,&  !Call diagonaliser
!            WORK, LWORK, RWORK, INFO )                      !
!      if(debug) write(6,*)'INFO....',info,'   i = ',i
!       do i1=1,mr                         ! this loop normalises the eigenvectors
!          vl=(0.d0,0.d0)                  ! from each diagonalisation.
!          do j=1,mr                       !
!             vl=vl+vr(j,i1)*vr(j,i1)      !
!          end do              ! j=1,mr    !
!          vnor=1.d0/sqrt(vl)              !
!          do j=1,mr
!             orightmp(j,i1)=vr(j,i1)                       !
!             vr(j,i1)=vr(j,i1)*vnor       !
!          end do              ! j=1,mr    !
!       end do                 ! i1=1,mr   !
!       do i1=1,mr
!          do ic=1,md
!             indp(ic)=0
!             dist(ic)=1.d100
!          end do
!          ic=0
!          ind(i1)=0
!          do i2=1,mr
!             vl=(0.d0,0.d0)
!             do j=1,mr
!                vl=vl+(oldh(j,i1))*vr(j,i2)
!             end do
!             ene=sqrt(dble(conjg(vl)*vl))
!             !               write(81,103)i,i1,i2,ene
!             if (ene.gt.thre) then
!                ind(i1)=i2
!                ic=1+ic
!                vl=eigo(i2)-eigl(i1)
!                indp(ic)=i2
!                dist(ic)=dsqrt(dble(conjg(vl)*vl))
!                !                  write(83,*)i,i1,i2,ic,ene,dsqrt(dble(conjg(vl)*vl))
!                !                  write(84,*)i,i1,i2,ic,eigo(i1),eigl(i2)
!             end if
!          end do
!          if (ic.eq.0) then
!             write(6,*)'Warning: Missing contact in '
!             write(6,*)'         eig ...',i1
!             write(6,*)'         lmb ...',i
!          else if (ic.gt.1) then
!             write(6,*)'Warning: Multiple contact in '
!             write(6,*)'         eig ........',i1
!             write(6,*)'         lmb ........',i
!             write(6,*)'         Contacts ...',ic
!             write(6,*)'         Resolving multiple contact '
!             write(6,*)(dist(j),j=1,ic)
!             j1=indp(1)
!             do j=2,ic
!                if (dabs(dist(j)).gt.dabs(dist(j-1))) j1=indp(j)
!             end do
!             ind(i1)=j1
!          end if
!          emin=1.d100
!          do i2=1,mr
!             vl=eigo(i1)-eigl(i2)
!             ene=dsqrt(dble(conjg(vl)*vl))
!             if (ene.lt.emin) then
!                ind2(i1)=i2
!                emin=ene
!             end if
!          end do
!       end do ! do i1=1,mr
!       j1=0
!       do j=1,mr
!          j1=j1+ind(j)
!       end do
!       j2=mr*(mr+1)/2
!       do i1=1,mr
!          do i2=1,mr
!             if(ind(i1).eq.0)then
!                oldh(i2,i1)=vr(i2,i1)
!             else
!                oldh(i2,i1)=vr(i2,ind(i1))
!                origh(i2,i1)=orightmp(i2,ind(i1))
!             end if
!          end do
!          if(ind(i1).eq.0)then
!             write(unit=diagunit)lambda,i1,i,eigl(i1)*conv
!             !             write(unit=diagunit,fmt=*)lambda,i1,i,eigl(i1)*conv
!             eigo(i1)=eigl(i1)
!          else
!             !write(87,100)i1,i,eigl(ind(i1))*conv
!             write(unit=diagunit)lambda,i1,i,eigl(ind(i1))*conv
!             !             write(unit=diagunit,fmt=*)lambda,i1,i,eigl(ind(i1))*conv
!             eigo(i1)=eigl(ind(i1))
!          end if
!       end do
!       ! WRITE EIGENVECTORS TO FILE IF REQUESTED
!       if (in%outputvec) then
!          write(unit=diagvecunit,rec=i) origh!(:,:)
!       end if
!
!    end do ! do i=2,in%nlambda  (outer loop)
!100 format(2I5,d20.13,2x,d20.13)
!102 format(2I5,f12.6)
!103 format(3I5,d25.18)
!108 format(3I5,d25.18,2x,d25.18)
!104 format(' ----------  ',d25.18)
!    print*
!    print*, '         diagonalisation complete'
!    deallocate( ind,indp,ind2,ind3)
!    deallocate( eig0,w,rwork,dist)
!    deallocate( ham,oldh,eigl,work)
!    deallocate( vr,eigo )
!    !    close(86)
!  end subroutine cap_matrix_diag_ser



  subroutine simple_cap_123(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)   :: r1,r2,xcos
    integer,intent(in)   :: n_dr
    real(8),intent(out)  :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    call q_to_j(q1,q2,q3,dummy1,r2j3,dummy2,'3-12')
    count=0
    do i=1,n_dr
          v1= simple_cap(r2j1,rmax(1),deltar_small(i))
          v2= simple_cap(r2j2,rmax(2),deltar_small(i))
          v3= simple_cap(r2j3,rmax(3),deltar_small(i))
         capout(i) = v1 + v2 + v3
    end do
  end subroutine simple_cap_123

  subroutine simple_cap_12_3(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple_cap(r2j1,rmax(1),deltar_small(i))
    v2= simple_cap(r2j2,rmax(2),deltar_small(i))
    capout(i) = v1 + v2

    end do
  end subroutine simple_cap_12_3

  subroutine simple_cap_1_23(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
!    allocate(capout(n_dr))
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    do i=1,n_dr
    v1= simple_cap(r2j1,rmax(1),deltar_small(i))
    capout(i) = v1
    end do
  end subroutine simple_cap_1_23

  subroutine simple_cap_2_31(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple_cap(r2j1,rmax(2),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple_cap_2_31

  subroutine simple_cap_3_12(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= simple_cap(r2j1,rmax(3),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple_cap_3_12

  subroutine simple8_cap_123(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    call q_to_j(q1,q2,q3,dummy1,r2j3,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= simple8_cap(r2j1,rmax(1),deltar_small(i))
    v2= simple8_cap(r2j2,rmax(2),deltar_small(i))
    v3= simple8_cap(r2j3,rmax(3),deltar_small(i))
    capout(i) = v1 + v2 + v3

    end do
  end subroutine simple8_cap_123

  subroutine simple8_cap_12_3(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple8_cap(r2j1,rmax(1),deltar_small(i))
    v2= simple8_cap(r2j2,rmax(2),deltar_small(i))
    capout(i) = v1 + v2

    end do
  end subroutine simple8_cap_12_3

  subroutine simple8_cap_1_23(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count
    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    count=0
    do i=1,n_dr
    v1= simple8_cap(r2j1,rmax(1),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple8_cap_1_23

  subroutine simple8_cap_2_31(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple8_cap(r2j1,rmax(2),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple8_cap_2_31

  subroutine simple8_cap_3_12(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= simple8_cap(r2j1,rmax(3),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple8_cap_3_12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine simple3_cap_123(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    call q_to_j(q1,q2,q3,dummy1,r2j3,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= simple3_cap(r2j1,rmax(1),deltar_small(i))
    v2= simple3_cap(r2j2,rmax(2),deltar_small(i))
    v3= simple3_cap(r2j3,rmax(3),deltar_small(i))
    capout(i) = v1 + v2 + v3

    end do
  end subroutine simple3_cap_123

  subroutine simple3_cap_12_3(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple3_cap(r2j1,rmax(1),deltar_small(i))
    v2= simple3_cap(r2j2,rmax(2),deltar_small(i))
    capout(i) = v1 + v2

    end do
  end subroutine simple3_cap_12_3

  subroutine simple3_cap_1_23(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    do i=1,n_dr
    v1= simple3_cap(r2j1,rmax(1),deltar_small(i))
    capout(i) = v1
    end do
  end subroutine simple3_cap_1_23

  subroutine simple3_cap_2_31(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= simple3_cap(r2j1,rmax(2),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple3_cap_2_31

  subroutine simple3_cap_3_12(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= simple3_cap(r2j1,rmax(3),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine simple3_cap_3_12



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine manop_cap_123(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    call q_to_j(q1,q2,q3,dummy1,r2j3,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1 = manop_cap(r2j1,rmax(1),deltar_small(i))
    v2 = manop_cap(r2j2,rmax(2),deltar_small(i))
    v3 = manop_cap(r2j3,rmax(3),deltar_small(i))
    capout(i) = v1 + v2 + v3
    end do
  end subroutine manop_cap_123

  subroutine manop_cap_12_3(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    call q_to_j(q1,q2,q3,dummy1,r2j2,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= manop_cap(r2j1,rmax(1),deltar_small(i))
    v2= manop_cap(r2j2,rmax(2),deltar_small(i))
    capout(i) = v1 + v2

    end do
  end subroutine manop_cap_12_3

  subroutine manop_cap_1_23(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'1-23')
    count=0
    do i=1,n_dr
    v1= manop_cap(r2j1,rmax(1),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine manop_cap_1_23

  subroutine manop_cap_2_31(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'2-31')
    count=0
    do i=1,n_dr
    v1= manop_cap(r2j1,rmax(2),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine manop_cap_2_31

  subroutine manop_cap_3_12(r1,r2,xcos,n_dr , capout)
    use class_geom ; use common_parameters ; implicit none
    real(8),intent(in)    :: r1,r2,xcos
    integer,intent(in)    :: n_dr
    real(8),intent(out)   :: capout(n_dr)
    real(8) :: q1,q2,q3,r2j,r2j1,r2j2,r2j3,v1,v2,v3
    real(8) :: dummy1,dummy2,ga,gb,gc
    integer :: i,count

    if(idia.lt.0) then
       call r_to_q(r1,r2,xcos,q1,q2,q3) ! radau to bondlength
    else
       call j_to_q(r1,r2,xcos,q1,q2,q3)
    end if
    call q_to_j(q1,q2,q3,dummy1,r2j1,dummy2,'3-12')
    count=0
    do i=1,n_dr
    v1= manop_cap(r2j1,rmax(3),deltar_small(i))
    capout(i) = v1

    end do
  end subroutine manop_cap_3_12

  function simple_cap(r,rmax,deltar) result(cap_out)
    implicit none
    real(8),intent(in) :: r,rmax,deltar
    real(8) :: cap_out, rmin
    rmin=rmax-deltar
    if(r.le.rmin.or.r.gt.rmax) then
       cap_out=0.0d0
    else
       cap_out = ((r-rmin)/deltar)**2
    endif
  end function simple_cap

  function simple8_cap(r,rmax,deltar) result(cap_out)
    implicit none
    real(8),intent(in) :: r,rmax,deltar
    real(8) :: cap_out, rmin
    rmin=rmax-deltar
    if(r.le.rmin.or.r.gt.rmax) then
       cap_out=0.0d0
    else
       cap_out = ((r-rmin)/deltar)**8
    endif
  end function simple8_cap

  function simple3_cap(r,rmax,deltar) result(cap_out)
    implicit none
    real(8),intent(in) :: r,rmax,deltar
    real(8) :: cap_out, rmin, rmaxx
!added by bruno to check the HOCl results
rmaxx=rmax+0.05d0
rmin=rmaxx-deltar
    if(r.le.rmin.or.r.gt.rmax) then
       cap_out=0.0d0
    else
       cap_out = ((r-rmin)/deltar)**3
    endif
  end function simple3_cap

  function manop_cap(r,rmax,deltar) result(cap_out)
    !see J.Chem.Phys,vol 117,page 9552,2002 (David E. Manolopoulos)
    !and J.Chem.Phys,vol 120,page 2247,2004 (David E. Manolopoulos)
    !(specifically the second one for this absorbing potential)
    !
    implicit none
    real(8) :: cap_out
    real(8),intent(in) :: r,deltar,rmax
    real(8) :: x,rmin,rmaxx
    ! value calculated from definite integral (see paper)
    real(8), parameter ::   c=2.6220575542921198105d0

    rmaxx=rmax+0.1 ! to remove the singularity from the grid
    rmin=rmaxx-deltar
    if(r.le.(rmin).or.r.ge.(rmaxx)) then
       cap_out = 0.0d0
    else
       x = (c*(r-(rmin)))/deltar
       cap_out = 4/((c-x)**2) + 4/((c+x)**2) - 8/(c**2)
    endif
    return
  end function manop_cap

end module class_complex_h

