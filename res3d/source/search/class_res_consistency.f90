Module class_res_consistency
  use common_parameters
  use class_input
  use class_stash
  implicit none
  private
  type res_consist_type
     type(stash_int)  :: res_lbl           ! label of a resonance
     type(stash_real) :: reE,imE,dreE,dimE
     type(stash_real)  :: retol,imtol      !tolerance in the stability point
  end type res_consist_type
  public :: res_consistency
contains

  subroutine res_consistency(in)
    implicit none
    type(input) :: in
    type(res_consist_type),allocatable :: f(:),total(:),inside(:)
    !auxiliary res_consist_type used in averages and final result
    type(res_consist_type) :: tmp_res
    real(8)::rmin, rmin2 !these have to be changed into reals
    integer::nrmin

    real(8)::dummy,rmax2,rmax3

    integer :: i,j,k
    integer :: nstates
    integer :: first_state,last_state
    character(len=2) :: char
    integer :: nfiles, start_unit
    integer :: file_index,count_loop
    !used in reading the data from the files:
    integer :: int_buff
    real(8) :: real_buff1, real_buff2, real_buff3, real_buff4
    !order vector used after swapping the stash
    integer,allocatable :: order(:)
    !repeat spurious resonances removal flag
    logical :: repeat
    !auxiliary stash
    type(stash_real) :: aux_stash_real
    type(stash_int) :: aux_stash_int

    nfiles = ndeltar 
    nstates = diag_high_state-diag_low_state+1
    first_state=1!diag_low_state
    last_state=nstates!diag_high_state
    start_unit=20 !file unit starting point -1 

    !initialise data structures
    allocate(f(nfiles)) 
    do i=1,nfiles
        call stash_init(f(i)%res_lbl )
        call stash_init(f(i)%reE     )
        call stash_init(f(i)%imE     )
        call stash_init(f(i)%retol   )
        call stash_init(f(i)%imtol   )
    end do
    allocate(total(nstates),inside(nstates))

    do i=1,nstates
        call stash_init(total(i)%res_lbl )
        call stash_init(total(i)%reE     )
        call stash_init(total(i)%imE     )
        call stash_init(total(i)%dreE    )
        call stash_init(total(i)%dimE    )
        call stash_init(total(i)%retol     )
        call stash_init(total(i)%imtol     )
        call stash_init(inside(i)%res_lbl )
        call stash_init(inside(i)%reE     )
        call stash_init(inside(i)%imE     )
        call stash_init(inside(i)%dreE    )
        call stash_init(inside(i)%dimE    )
        call stash_init(inside(i)%retol     )
        call stash_init(inside(i)%imtol     )
    end do
        call stash_init(tmp_res%res_lbl )
        call stash_init(tmp_res%reE     )
        call stash_init(tmp_res%imE     )
        call stash_init(tmp_res%dreE    )
        call stash_init(tmp_res%dimE    )
        call stash_init(tmp_res%retol     )
        call stash_init(tmp_res%imtol     )


    ! f is the data structure that keeps the information read from the files
    ! +1 comes from the need of workspace
    file_open:  do i = 1 , nfiles
       file_index = start_unit+i
       write(char,"(I2.2)") i
       open(unit=file_index &
            ,file='res'//char//'.bcs' &
            ,form='formatted'&
            ,status='old')
    end do file_open

    data_read: do i = 1 , nfiles  !do data structure and work space allocation 
       file_index = start_unit+i
j=0
       do                                  ! -This loop runs through the 

!          print*,'i= ',i 
          read(unit=file_index,fmt=*,END=100) &
               int_buff, real_buff1, real_buff2, real_buff3, real_buff4
j=j+1
          call stash_add(f(i)%res_lbl,int_buff ) 
          call stash_add(f(i)%reE,real_buff1   )
          call stash_add(f(i)%imE,real_buff2   )
          call stash_add(f(i)%retol,real_buff3 )
          call stash_add(f(i)%imtol,real_buff4 )
!print"(I5.5,4(1x,F12.5))", int_buff, real_buff1, real_buff2, real_buff3, real_buff4
!print"(I5.5,4(1x,F12.5))",f(i)%res_lbl%value(j)&
!                         ,f(i)%reE%value(j)    &
!                         ,f(i)%imE%value(j)    &
!                         ,f(i)%dreE%value(j)   &
!                         ,f(i)%dimE%value(j)
       end do
100    continue                         
    end do data_read


    
    !group resonances by state label:
    res_match:  do i = 1 , nfiles                    ! this looks at the state label in each file and tries to 
       do j = 1, f(i)%res_lbl%size                   ! mach it to the same label in the following file, 
          scan_loop: do k = 1 , nstates
             if( f(i)%res_lbl%value(j) .eq. k+first_state-1 )then
                call stash_add( total(k)%res_lbl , f(i)%res_lbl%value(j) )
                call stash_add( total(k)%reE     , f(i)%reE%value(j)     )
                call stash_add( total(k)%imE     , f(i)%imE%value(j)     )
                call stash_add( total(k)%dreE    , 0.0d0    )
                call stash_add( total(k)%dimE    , 0.0d0    )
                call stash_add( total(k)%retol    , f(i)%retol%value(j)    )
                call stash_add( total(k)%imtol    , f(i)%imtol%value(j)    )
             end if
          end do scan_loop
       end do
    end do res_match

!do k = 1 , nstates              
!do j = 1, total(k)%res_lbl%size  
!print"(I5.5,4(1x,F12.5))",total(k)%res_lbl%value(j)&
!                         ,total(k)%reE%value(j)    &
!                         ,total(k)%imE%value(j)    &
!                         ,total(k)%dreE%value(j)   &
!                         ,total(k)%dimE%value(j)
!end do
!end do 

    ! deal with spurious resonances
    
    !copy one resonance stash into another...
    do i=1,nstates
       call stash_copy(total(i)%res_lbl,inside(i)%res_lbl) ! copy the contents of total into inside
       call stash_copy(total(i)%reE,inside(i)%reE)
       call stash_copy(total(i)%imE,inside(i)%imE)
       call stash_copy(total(i)%dreE,inside(i)%dreE)
       call stash_copy(total(i)%dimE,inside(i)%dimE)
       call stash_copy(total(i)%retol,inside(i)%retol)
       call stash_copy(total(i)%imtol,inside(i)%imtol)
    end do

!do k = 1 , nstates              
!do j = 1, total(k)%res_lbl%size  
!print"(I5.5,4(1x,F12.5))",inside(k)%res_lbl%value(j)&
!                         ,inside(k)%reE%value(j)    &
!                         ,inside(k)%imE%value(j)    &
!                         ,inside(k)%dreE%value(j)   &
!                         ,inside(k)%dimE%value(j)
!end do
!end do 



repeat = .true.    
    spurious_removal_position: do while(repeat)
       repeat = .false.
    ! Calculate the average value of each quantity for each state    
       if (associated(tmp_res%res_lbl%value)) then
          call stash_destroy(tmp_res%res_lbl)
          call stash_destroy(tmp_res%reE)
          call stash_destroy(tmp_res%imE)
          call stash_destroy(tmp_res%dreE)
          call stash_destroy(tmp_res%dimE)
          call stash_destroy(tmp_res%retol)
          call stash_destroy(tmp_res%imtol)
       end if

    avrg_calc:do i = 1 , nstates  
       if(inside(i)%res_lbl%size.gt.0) then
          call stash_add( tmp_res%res_lbl , inside(i)%res_lbl%value(1)    ) ! this should leave the label unaffected
          call stash_add( tmp_res%reE     , stash_sum( inside(i)%reE   ) / inside(i)%res_lbl%size    )
          call stash_add( tmp_res%imE     , stash_sum( inside(i)%imE   ) / inside(i)%res_lbl%size    )
          call stash_add( tmp_res%dreE    , 0.0d0    ) 
          call stash_add( tmp_res%dimE    , 0.0d0    ) 
          call stash_add( tmp_res%retol    ,stash_sum( inside(i)%retol   ) / inside(i)%res_lbl%size    ) !average the tolerances too
          call stash_add( tmp_res%imtol    ,stash_sum( inside(i)%imtol   ) / inside(i)%res_lbl%size    ) !
       else
          call stash_add( tmp_res%res_lbl , 0    ) ! this should leave the label unaffected
          call stash_add( tmp_res%reE     , 0.0d0    )
          call stash_add( tmp_res%imE     , 0.0d0    )
          call stash_add( tmp_res%dreE    , 0.0d0    )
          call stash_add( tmp_res%dimE    , 0.0d0    )
          call stash_add( tmp_res%retol     , 0.0d0    )
          call stash_add( tmp_res%imtol     , 0.0d0    )
       end if

    end do avrg_calc

!do k = 1 , nstates              
!print"(I5.5,4(1x,F12.5))",tmp_res%res_lbl%value(k)&
!                         ,tmp_res%reE%value(k)    &
!                         ,tmp_res%imE%value(k)    &
!                         ,tmp_res%dreE%value(k)   &
!                         ,tmp_res%dimE%value(k)
!end do 

    
    ! Calculate the distance of each resonance from the average
    dist_calc:do i = 1 , tmp_res%res_lbl%size  
       do j=1,inside(i)%res_lbl%size
          if(tmp_res%res_lbl%value(i).gt.0) then
             inside(i)%dreE%value(j) =  abs( inside(i)%reE%value(j) - tmp_res%reE%value(i) )
          end if
       end do
    end do dist_calc
!do k = 1 , nstates              
!do j = 1, total(k)%res_lbl%size  
!print"(I5.5,4(1x,F12.5))",inside(k)%res_lbl%value(j)&
!                         ,inside(k)%reE%value(j)    &
!                         ,inside(k)%imE%value(j)    &
!                         ,inside(k)%dreE%value(j)   &
!                         ,inside(k)%dimE%value(j)
!end do
!end do 
    
    !sort the resonance stashes by error of position ( real part ) 
    dist_sort:do i = 1 , nstates  
          if(inside(i)%res_lbl%size.gt.0) then
             allocate(order(inside(i)%res_lbl%size))
             do j=1,inside(i)%res_lbl%size
                order(j)=j
             end do
             call stash_sort(inside(i)%dreE,order)
!print*,order
             call stash_reorder(order,inside(i)%reE)
             call stash_reorder(order,inside(i)%imE)
             call stash_reorder(order,inside(i)%dimE)
             call stash_reorder(order,inside(i)%retol)
             call stash_reorder(order,inside(i)%imtol)
             deallocate(order)
          end if
    end do dist_sort

    !Check if not all resonances lie within the real error range, remove the largest error resonance.
    dist_cut:do i = 1 , nstates  
       if(inside(i)%res_lbl%size.gt.0) then 
          if(inside(i)%dreE%value(inside(i)%dreE%size).gt.tmp_res%retol%value(i)) then
             call stash_sub(inside(i)%res_lbl)
             call stash_sub(inside(i)%reE)
             call stash_sub(inside(i)%imE)
             call stash_sub(inside(i)%dreE)
             call stash_sub(inside(i)%dimE)
             call stash_sub(inside(i)%retol)
             call stash_sub(inside(i)%imtol)
             repeat = .true.
          end if
       end if
    end do dist_cut
    
 end do spurious_removal_position

repeat=.true.
    spurious_removal_width: do while(repeat)
       repeat = .false.
    ! Calculate the average value of each quantity for each state    
       if (associated(tmp_res%res_lbl%value)) then
          call stash_destroy(tmp_res%res_lbl)
          call stash_destroy(tmp_res%reE)
          call stash_destroy(tmp_res%imE)
          call stash_destroy(tmp_res%dreE)
          call stash_destroy(tmp_res%dimE)
          call stash_destroy(tmp_res%retol)
          call stash_destroy(tmp_res%imtol)
       end if

    avrg_calc_w:do i = 1 , nstates  

       if(inside(i)%res_lbl%size.gt.0) then
          call stash_add( tmp_res%res_lbl , inside(i)%res_lbl%value(1)    ) ! this should leave the label unaffected
          call stash_add( tmp_res%reE     , stash_sum( inside(i)%reE   ) / inside(i)%res_lbl%size    )
          call stash_add( tmp_res%imE     , stash_sum( inside(i)%imE   ) / inside(i)%res_lbl%size    )
          call stash_add( tmp_res%dreE    , 0.0d0    )
          call stash_add( tmp_res%dimE    , 0.0d0    )
          call stash_add( tmp_res%retol    ,stash_sum( inside(i)%retol   ) / inside(i)%res_lbl%size    ) !average the tolerances too
          call stash_add( tmp_res%imtol    ,stash_sum( inside(i)%imtol   ) / inside(i)%res_lbl%size    ) !
       else
          call stash_add( tmp_res%res_lbl , 0    ) ! this should leave the label unaffected
          call stash_add( tmp_res%reE     , 0.0d0    )
          call stash_add( tmp_res%imE     , 0.0d0    )
          call stash_add( tmp_res%dreE    , 0.0d0    )
          call stash_add( tmp_res%dimE    , 0.0d0    )
          call stash_add( tmp_res%retol    ,0.0d0    ) !average the tolerances too
          call stash_add( tmp_res%imtol    ,0.0d0    ) !
       end if

    end do avrg_calc_w
    
    ! Calculate the distance of each resonance from the average
    dist_calc_w:do i = 1 , tmp_res%res_lbl%size  
       do j=1,inside(i)%res_lbl%size
          if(tmp_res%res_lbl%value(i).gt.0) then
             inside(i)%dimE%value(j) =  abs( inside(i)%imE%value(j) - tmp_res%imE%value(i) )    
             inside(i)%dreE%value(j) =  abs( inside(i)%reE%value(j) - tmp_res%reE%value(i) )            
          end if
       end do
    end do dist_calc_w
    
    !sort the resonance stashes by error of position ( real part ) 
    dist_sort_w:do i = 1 , nstates  
       if(inside(i)%res_lbl%size.gt.0) then
          allocate(order(inside(i)%res_lbl%size))
          do j=1,inside(i)%res_lbl%size
             order(j)=j
          end do
          call stash_sort(inside(i)%dimE,order)
          call stash_reorder(order,inside(i)%reE)
          call stash_reorder(order,inside(i)%imE)
          call stash_reorder(order,inside(i)%dreE)
          call stash_reorder(order,inside(i)%retol)
          call stash_reorder(order,inside(i)%imtol)
          deallocate(order)
       end if
    end do dist_sort_w
    
    !Check if not all resonances lie within the real error range, remove the largest error resonance.
    dist_cut_w:do i = 1 , nstates  
       if(inside(i)%res_lbl%size.gt.0) then
          if(inside(i)%dimE%value(inside(i)%dimE%size).gt.7*tmp_res%imtol%value(i)) then
             call stash_sub(inside(i)%res_lbl)
             call stash_sub(inside(i)%reE)
             call stash_sub(inside(i)%imE)
             call stash_sub(inside(i)%dreE)
             call stash_sub(inside(i)%dimE)
             call stash_sub(inside(i)%retol)
             call stash_sub(inside(i)%imtol)
             repeat = .true.
          end if
       end if
    end do dist_cut_w
    
 end do spurious_removal_width
 
 
 
 !write results to file
    open(unit=10,file='res_consist.bcs',status='new') 
    !     write(10,*)'likelyhood, position, width, pos_error, wid_error' 
    file_write:do i = 1 , nstates

       if(inside(i)%res_lbl%size*100.0d0/nfiles.gt.0 ) then
          !sort the poistion error to get its maximum on output
          allocate(order(inside(i)%dreE%size))
          call stash_sort(inside(i)%dreE,order)
          deallocate(order)
          write(10,'(i4.4,2x,f5.1,2x,f10.3,2x,ES12.3,2x,f7.3,2x,ES12.3)') &
               i+first_state-1 &
               ,inside(i)%res_lbl%size*100.0d0/nfiles&
               ,tmp_res%reE%value(i)    &
               ,tmp_res%imE%value(i)   &
               ,inside(i)%dreE%value(inside(i)%dreE%size)  &
               ,inside(i)%dimE%value(inside(i)%dimE%size)
       end if
    end do file_write

    close(10)
  end subroutine res_consistency


end Module class_res_consistency


