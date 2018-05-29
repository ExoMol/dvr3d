module class_wavefile_to_elements
  use class_error
  use class_wavefile
  use class_kblock
  use class_kblock_to_elements
  implicit none
  private
  public :: elements_from_wavefile
  !function elements_from_wavefile(file,low_state,high_state,operator) result(elements) 
  !  implicit none
  !  character(len=200),intent(in) :: file
  !  integer :: low_state,high_state
  !  real(8), external :: operator
  !  real(8),allocatable :: elements(:) 
  !end function elements_from_wavefile
contains
  function elements_from_wavefile(file,low_state,high_state,operator) result(elements) 
    ! CALCULATES THE SYMMETRIC MATRIX ELEMENTS OF 'OPERATOR' USING THE WAVEFUNCTIONS
    ! OBTAINED FROM THE SERIAL VERSIONS OF EITHER ROTLEV OR DVR3D
    ! THE INPUT IS DEFINED AS:
    !
    ! FILE       = FILE NAME OR PATH TO FILE WHERE THE WAVEFUNCTIONS ARE STORED
    ! LOW_STATE  = INTEGER IDENTIFYING THE LOWEST STATE USED TO BUILD THE MATRIX
    ! HIGH_STATE = INTEGER IDENTIFYING THE HIGHEST STATE USED TO BUILD THE MATRIX
    ! OPERATOR   = MUST BE A REAL FUNCTION OPERATOR(R1,R2,COS(THETA)) IN JACOBI COORDINATES
    !
    !THE OUTPUT OF THIS FUNCTION IS A VECTOR elements( nstates*(nstates+1)/2 )
    !                                 WITH   nstates = high_state-low_state+1
    implicit none
    character(len=200),intent(in) :: file
    integer :: low_state,high_state
    real(8), external :: operator
    real(8), allocatable :: op_matrix(:,:,:)
    real(8),allocatable :: elements(:) 
    real(8),allocatable :: elements_tmp(:,:)
    real(8),allocatable :: norms(:)
    integer :: l2, l1, chunk=1
    type(kblock), pointer :: block
    type(wavefile), pointer :: w
    integer :: k,unit,jk,count,nstates
    ! get the number of k blocks from the file
    nstates = high_state-low_state+1
    allocate(elements(nstates*(nstates+1)/2))
    allocate(elements_tmp(nstates,nstates))
    elements=0.0d0
    elements_tmp=0.0d0
    call wavefile_open(file,w)
    call wavefile_read_header(w)
    jk=w%jk
    do k = 1,jk
       call wavefile_read_kblock_to_mem(k,block,w)
       if (low_state.gt.block%neval.or.&
            high_state.gt.block%neval) then
          print*, 'low_state',low_state
          print*, 'high_state',high_state
          print*, 'block%neval',block%neval
          call error('elements_from_jacobi_wavefile:low_state.gt.block%neval.or.high_state.gt.block%neval')
       end if
       allocate( op_matrix( block%grid%nang,block%grid%npnt1 , block%grid%npnt2   ) )



          call operator_matrix_build( block , operator , op_matrix )


    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE( l2 ) &
    !$OMP SCHEDULE(DYNAMIC,CHUNK) 

          do l2=low_state,high_state
             do l1=low_state,l2
                elements_tmp(l1-low_state+1,l2-low_state+1) =  kblock_element(block,l1,op_matrix,l2)
                if (debug) print*,'elements_tmp(l1,l2)',l1,l2,elements_tmp(l1-low_state+1,l2-low_state+1)
             end do
          end do

    !$OMP END PARALLEL DO
          count=0
          do l2=low_state,high_state
             do l1=low_state,l2
                count=count+1
                elements(count) = elements(count)  + elements_tmp(l1-low_state+1,l2-low_state+1)
                if (debug) print*,'elements(count)',l1,l2,elements(count)
             end do
          end do


       deallocate(op_matrix)
       call kblock_destroy(block)
    end do

    call wavefile_close(w)
    call wavefile_destroy_header(w)

  end function elements_from_wavefile

end module class_wavefile_to_elements
