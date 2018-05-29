module class_stash
  use class_sort
  use class_error
  implicit none
  private
  type,public :: stash_int
     integer :: size=0
     integer, pointer :: value(:)
  end type stash_int

  type,public :: stash_real
     integer :: size=0
     real(8), pointer :: value(:)
  end type stash_real

  interface stash_init
     module procedure stash_init_real, stash_init_int
  end interface

  interface stash_add
     module procedure stash_add_real, stash_add_int
  end interface

  interface stash_destroy
     module procedure stash_destroy_int, stash_destroy_real
  end interface

  interface stash_copy
     module procedure stash_copy_real, stash_copy_int
  end interface

  interface stash_sub
     module procedure stash_sub_real,  stash_sub_int
  end interface

  interface stash_reorder
     module procedure stash_reorder_real, stash_reorder_int
  end interface

  interface stash_sort
     module procedure stash_sort_real, stash_sort_int
  end interface

  interface stash_sum
     module procedure stash_sum_real, stash_sum_int
  end interface

  public :: stash_init!(s)
  public :: stash_add!(s,value)
  public :: stash_sub!(s)
  public :: stash_sort!(s,order)
  public :: stash_sum!(s) result(real(8))
  public :: stash_reorder!(neworder,s)
  public :: stash_copy!(s,t)
  public :: stash_destroy!(s)
  
contains
  
  subroutine stash_init_real(s)
    type(stash_real),intent(inout) :: s
    s%size = 0
  end subroutine stash_init_real

  subroutine stash_init_int(s_i)
    type(stash_int),intent(inout) :: s_i
    s_i%size = 0 
  end subroutine stash_init_int


  subroutine stash_destroy_int(s_i)
    type(stash_int),intent(inout) :: s_i
    if(associated(s_i%value)) deallocate(s_i%value)
    s_i%size = 0
  end subroutine stash_destroy_int
  
  subroutine stash_destroy_real(s)
    type(stash_real),intent(inout) :: s
    if(associated(s%value)) deallocate(s%value)
    s%size=0
  end subroutine stash_destroy_real

  subroutine stash_copy_real(s,t)
    type(stash_real) :: s,t
    if(associated(t%value)) then ;deallocate(t%value) ;t%size=0; end if
       if(s%size.gt.0) then
          allocate(t%value(s%size))
          t%size=s%size
          t%value=s%value 
       end if
  end subroutine stash_copy_real

  subroutine stash_copy_int(s_i,t_i)
    type(stash_int) :: s_i,t_i
    if(associated(t_i%value)) then ;deallocate(t_i%value) ;t_i%size=0; end if
    if(s_i%size.gt.0) then
       allocate(t_i%value(s_i%size))
       t_i%size=s_i%size
       t_i%value=s_i%value 
    end if
  end subroutine stash_copy_int
  
  subroutine stash_add_real(s,value)
    type(stash_real) :: s
    real(8),intent(in) :: value
    real(8), allocatable :: aux_vec(:)

    if (associated(s%value)) then
       allocate(aux_vec(s%size))
       aux_vec = s%value
       deallocate(s%value)
       s%size=s%size+1
       allocate(s%value(s%size))
       s%value(1:(s%size-1))=aux_vec(:)
       s%value(s%size) = value 
       deallocate(aux_vec)

    else
       allocate(s%value(1))
       s%value(1) = value 
       s%size = 1 
    end if
  end subroutine stash_add_real
  
  subroutine stash_add_int(s_i,value_i)
    type(stash_int) :: s_i
    integer,intent(in) :: value_i
    integer, allocatable :: aux_vec(:)
    if (associated(s_i%value))then
       allocate(aux_vec(s_i%size))
       aux_vec = s_i%value
       deallocate(s_i%value)
       s_i%size=s_i%size+1
       allocate(s_i%value(s_i%size))
       s_i%value(1:(s_i%size-1))=aux_vec(:)
       s_i%value(s_i%size) = value_i 
       deallocate(aux_vec)
    else
       allocate(s_i%value(1))
       s_i%value(1) = value_i 
       s_i%size = 1 
    end if
  end subroutine stash_add_int

  subroutine stash_sub_int(s)
    type(stash_int) :: s
    
    real(8), allocatable :: aux_vec(:)
    
    if (associated(s%value))then
       allocate(aux_vec(s%size-1))
       aux_vec = s%value(1:(s%size-1))
       deallocate(s%value)
       s%size=s%size-1
       allocate(s%value(s%size))
       s%value=aux_vec
       deallocate(aux_vec)
    else
       call error('Attempted to remove element from empty stash. Aborting')
    end if
  end subroutine stash_sub_int

  subroutine stash_sub_real(s)
    type(stash_real) :: s
    real(8), allocatable :: aux_vec(:)
    
    if (associated(s%value))then
       allocate(aux_vec(s%size-1))
       aux_vec = s%value(1:(s%size-1))
       deallocate(s%value)
       s%size=s%size-1
       allocate(s%value(s%size))
       s%value=aux_vec
       deallocate(aux_vec)
    else
       call error('Attempted to remove element from empty stash. Aborting')
    end if
  end subroutine stash_sub_real

  
  subroutine stash_sort_int(s,neworder)
    type(stash_int) :: s
    integer,intent(inout) :: neworder(s%size)
    if (associated(s%value)) then
       call sort(s%value,neworder,1,s%size)
    else
       call error('Attempted to sort empty stash. Aborting')
    end if
  end subroutine stash_sort_int

  subroutine stash_sort_real(s,neworder)
    type(stash_real) :: s
    integer,intent(inout) :: neworder(s%size)
    if (associated(s%value)) then
       call sort(s%value,neworder,1,s%size)
    else
       call error('Attempted to sort empty stash. Aborting')
    end if
  end subroutine stash_sort_real
   
  subroutine stash_reorder_int(neworder,s_i)
    type(stash_int) :: s_i
    integer,intent(in) :: neworder(s_i%size)
    integer :: i
    type(stash_int) :: tmp
    if (associated(s_i%value)) then
       call stash_init_int(tmp)
       call stash_copy_int(s_i,tmp)
       do i=1,s_i%size
       s_i%value(i)=tmp%value(neworder(i))   
       end do
    else
       call error('Attempted to reorder empty stash. Aborting')
    end if
  end subroutine stash_reorder_int

  subroutine stash_reorder_real(neworder,s)
    type(stash_real) :: s
    integer,intent(in) :: neworder(s%size)
    integer :: i
    type(stash_real) :: tmp
    if (associated(s%value)) then
       call stash_init_real(tmp)
       call stash_copy_real(s,tmp)
       do i=1,s%size
       s%value(i)=tmp%value(neworder(i))   
       end do
    else
       call error('Attempted to reorder empty stash. Aborting')
    end if
  end subroutine stash_reorder_real

    
  function stash_sum_int(s) result(sum)
    type(stash_int) :: s
    real(8) :: sum
    integer :: i
    sum=0
    if (associated(s%value)) then
       do i=1,s%size
          sum = sum + s%value(i)
       end do
    end if
  end function stash_sum_int

  function stash_sum_real(s) result(sum)
    type(stash_real) :: s
    real(8) :: sum
    integer :: i
    sum=0
    if (associated(s%value)) then
       do i=1,s%size
          sum = sum + s%value(i)
       end do
    end if
  end function stash_sum_real


end module class_stash

! UNCOMMENT TO TEST. TO COMPILE DO:
! ifort -check all -traceback class_stash.f90
!program test_stash
!  use class_stash
!  implicit none
!  type(stash_real),pointer :: s,t
!  type(stash_int),pointer :: s_i,t_i
!
!  integer,allocatable :: order(:)
!
!  integer :: i
!
!!testing the real stash
!
!call stash_add(s,10.1d0)
!call stash_add(s,6.1d0)
!call stash_add(s,3.1d0)
!call stash_add(s,1.1d0)
!call stash_add(s,9.1d0)
!call stash_add(s,11.1d0)  
!call stash_add(s,12.1d0)  
!call stash_add(s,14.1d0)
!!call stash_add(s,0.1d0)    
!
!call stash_copy(s,t)
!
!allocate(order(s%size))
!
!do i=1,s%size
!order(i)=i
!end do
!
!write(6,*) 'unsorted stash'
!do i=1,s%size
!write(6,*) order(i),s%value(i) 
!end do
!
!call stash_sort(s,order)
!
!write(6,*) 'sorted stash'
!do i=1,s%size
!write(6,*) order(i),s%value(i)
!end do
!
!write(6,*) 'stash from sorted indexes'
!call stash_reorder(order,t)
!do i=1,t%size
!write(6,*) order(i),t%value(i)
!end do
!
!write(6,*) 'stash summation'
!write(6,*) stash_sum(t)
!
!
!do while(s%size.gt.0)
!write(6,"('subtracting stash ',I2)") s%size
!call stash_sub(s)
!write(6,*) s%value
!end do
!
!write(6,*) 'destroying stashes'
!call stash_destroy(s)
!call stash_destroy(t)
!
!
!write(6,*) 'testing the integer stash'
!
!
!call stash_add(s_i,10)
!call stash_add(s_i,6)
!call stash_add(s_i,3)
!call stash_add(s_i,1)
!call stash_add(s_i,9)
!call stash_add(s_i,11)  
!call stash_add(s_i,12)  
!call stash_add(s_i,14)
!!call stash_add(s,0.1d0)    
!
!call stash_copy(s_i,t_i)
!
!do i=1,s_i%size
!order(i)=i
!end do
!
!write(6,*) 'unsorted stash'
!do i=1,s_i%size
!write(6,*) order(i),s_i%value(i) 
!end do
!
!call stash_sort(s_i,order)
!
!write(6,*) 'sorted stash'
!do i=1,s_i%size
!write(6,*) order(i),s_i%value(i)
!end do
!
!write(6,*) 'stash from sorted indexes'
!call stash_reorder(order,t_i)
!do i=1,t_i%size
!write(6,*) order(i),t_i%value(i)
!end do
!
!write(6,*) 'stash summation'
!write(6,*) stash_sum(t_i)
!
!
!do while(s_i%size.gt.0)
!write(6,"('subtracting stash ',I2)") s_i%size
!call stash_sub(s_i)
!write(6,*) s_i%value
!end do
!
!write(6,*) 'destroying stashes'
!call stash_destroy(s_i)
!call stash_destroy(t_i)
!
!
!write(6,*) 'test terminated successfuly :)'
!
!end program test_stash
