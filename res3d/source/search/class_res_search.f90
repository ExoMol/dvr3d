module class_res_search
  use common_parameters
  use class_curvs
  use class_input
  use class_res_consistency
  use class_stash
  implicit none
  private
  type,public :: res_srch
     type(stash_real) :: l,reE,imE
  end type res_srch
  public :: res_search
  !subroutine res_search(in)
  !type(input) :: in
  !end subroutine res_search
contains

  subroutine res_search(in)

    implicit none
    type(input) :: in
    type(curvs), pointer :: cs
    character(len=200) :: curvsfilename
    character(len=200) :: resfilename
    character(len=2) :: ncap_char
    integer,parameter::curvsunit=130
    integer,parameter::resunit=131
    integer :: i

    do i=1,ndeltar
       write(ncap_char,fmt="(I2.2)") i
       curvsfilename='curvs'//ncap_char//'.bcs'
       resfilename='res'//ncap_char//'.bcs'
       print*,'start resonance search for deltar= ', deltar(i)
       open ( unit=curvsunit , file = trim(curvsfilename) , status='old',form='unformatted' )
       open ( unit=resunit , file = trim(resfilename) ,form='formatted' )
       call curvs_read(curvsunit,cs)
       call res_search_deltar(resunit,cs,in)
       call curvs_destroy(cs)
       close(curvsunit)
       close(resunit)
    end do

    call res_consistency(in)

  end subroutine res_search


  subroutine res_search_deltar(resunit,cs,in)
    implicit none
    integer, intent(in) :: resunit
    type(input), intent(in) :: in
    type(curvs), intent(in) :: cs

    logical :: bool
    integer :: state,step,i,j,k,cnt,cnt_crv,cnt_crvtr,cnt_imEmin,n,odd_n,stt
    type(res_srch) :: crv, crvtr, imEmin
    !    real(8),allocatable :: crv(:,:), crvtr(:,:),imEmin(:,:),aux_vec(:,:)
    real(8),allocatable :: reE(:), imE(:), dimE(:) ! real and imaginary parts of the energy
    real(8),allocatable :: dE(:) ! lambda delta for the energy
    real(8),allocatable :: l(:)  ! lambda value
    real(8),allocatable :: curve(:),lvalue(:) !
    real(8),allocatable :: dEdl(:) !
    real(8) :: immin, immax,remin,remax,retol,imtol
    real(8) :: crv1,crv2,crv3,crv4,crv5,crv6,crv7
    real(8) :: crvtr1,crvtr2,crvtr3,crvtr4,crvtr5,crvtr6,crvtr7
    !character(len=150) :: filename
    character (len=4) :: char_state

    call stash_init(crv%l  ) 
    call stash_init(crv%reE) 
    call stash_init(crv%imE) 
    call stash_init(crvtr%l  ) 
    call stash_init(crvtr%reE) 
    call stash_init(crvtr%imE) 
    call stash_init(imEmin%l  )
    call stash_init(imEmin%reE)
    call stash_init(imEmin%imE)

    n = in%res%resmode
    odd_n = odd(n)
    !    print* , "-----" , n , odd_n
    !   stop
    allocate( reE(odd_n),imE(odd_n),dimE(odd_n)) ! real and imaginary parts of the energy
    allocate( dE(odd_n)                 ) ! lambda delta for the energy
    allocate( l(odd_n)                  ) ! lambda value
    allocate( curve(odd_n),lvalue(odd_n)) !
    allocate( dEdl(odd_n)               ) !

    !filename=trim ( in%char%diagfile )

    !open (unit=21,file="res/"//trim(filename)//"_cm1.bcs",status='new')
    stt=0
    outer_step_loop:    do state = 1,n_diag_states !diag_low_state , diag_high_state
       stt=stt+1

       reE(1) = cs%crvr(1,stt)
       imE(1) = cs%crvi(1,stt)
       dimE(1) = cs%dcrvi(1,stt)
       !PRINT*,reE(1),imE(1)
       !stop
       call lambda(in%diag%h,in%diag%chi,1,l(1)) 
       !print*,'-----', in%diag%h,in%diag%chi,1,l(1)       
       !print*,'-----',ree(1),ime(1),cs%lambda(1)
       !stop
       do step=2,odd_n
          !read(20,*) reE(step), imE(step)
          reE(step) = cs%crvr(step,stt)
          imE(step) = cs%crvi(step,stt)
          dimE(step) = cs%dcrvi(step,stt)

          !          print*,reE(step),imE(step)
          !          stop
          call lambda(in%diag%h,in%diag%chi,step,l(step))
          call deltaE(reE(step-1),reE(step),imE(step-1),imE(step),dE(step))
          dEdl(step) = dE(step)/(l(step)-l(step-1)) 
          !print*,'-----', reE(step-1),reE(step),imE(step-1),imE(step),dE(step)

          !read(22,*) lvalue(step),curve(step)
          lvalue(step) = cs%lambda(step)
          curve(step)  = cs%crvtr(step,stt)

          !print*, '-----',  lvalue(step),curve(step)
       end do
       !stop
       !from here on the buffered data is used to calculate the differences
       cnt_crv=0
       cnt_crvtr=0
       cnt_imEmin=0
       immin = 100000000000.0d0
       immax = -100000000000.0d0
       remin = 100000000000.0d0
       remax = -100000000000.0d0
       do step= odd_n + 1 , in%diag%nlambda
          !print*,'----', l
          do i = 1 , odd_n - 1
             l(i)=l(i+1)
             lvalue(i)=lvalue(i+1)
             dE(i)=dE(i+1)
             reE(i)=reE(i+1)
             imE(i)=imE(i+1)
             dimE(i)=dimE(i+1)
             curve(i)=curve(i+1)
             dEdl(i)=dEdl(i+1)
          end do
          !print*,'----', l
          !stop

          !read(20,*) reE( odd_n ), imE( odd_n )
          reE( odd_n ) = cs%crvr(step,stt)
          imE( odd_n ) = cs%crvi(step,stt)
          dimE( odd_n ) = cs%dcrvi(step,stt)
          !read(22,*) lvalue( odd_n ), curve( odd_n )
          lvalue( odd_n ) = cs%lambda(step)
          curve( odd_n )  = cs%crvtr(step,stt)

          call lambda(in%diag%h,in%diag%chi,step,l( odd_n ))
          call deltaE(reE( odd_n -1 ),reE( odd_n ),imE( odd_n - 1 ),imE( odd_n ),dE( odd_n ))
          dEdl(odd_n) = dE(odd_n)!/(l(odd_n)-l(odd_n-1)) 

          !get the maximum values for the real part and the complex part
          !to setup the resonance detection window: small for spoting
          ! large for checking consistency.

          if(reE(odd_n).gt.remax)remax=reE(odd_n)
          if(reE(odd_n).lt.remin)remin=reE(odd_n)
          if(imE(odd_n).gt.immax)immax=imE(odd_n)
          if(imE(odd_n).lt.immin)immin=imE(odd_n)

          ! test the point concentration
          bool=.true.
          do i = n-n/10, n - 2  

             bool = bool .and. dEdl(odd_n-i) &
                  .gt.dEdl(odd_n-1-i) &
                  .and.dEdl(3+i) &
                  .lt.dEdl(2+i)
          end do
          if(.not.bool) then
             crv1=0.0D0
             crv2=0.0D0
             crv3=0.0D0
             crv4=0.0D0
             crv5=0.0D0
             crv6=0.0D0
             crv7=0.0D0
             do i=1,(2*n)/21
                crv1=crv1+dEdl(i+1+6*(2*n)/21)
                crv2=crv2+dEdl(i+1+7*(2*n)/21)
                crv3=crv3+dEdl(i+1+9*(2*n)/21)
                crv4=crv4+dEdl(i+1+10*(2*n)/21)
                crv5=crv3+dEdl(i+1+11*(2*n)/21)
                crv6=crv6+dEdl(i+1+12*(2*n)/21)
                crv7=crv7+dEdl(i+1+13*(2*n)/21)
             end do
             if( crv1 > crv2 .and. crv2 > crv3 .and.&
                  crv3 > crv4 .and. crv4 < crv5 .and.&
                  crv5 < crv6 .and. crv6 < crv7 ) bool=.true.
          end if
          if (.not.bool) then
             crv1=0.0D0
             crv2=0.0D0
             crv3=0.0D0
             crv4=0.0D0
             crv5=0.0D0
             crv6=0.0D0
             crv7=0.0D0
             do i=1,(2*n)/7
                crv1=crv1+dEdl(i+1)
                crv2=crv2+dEdl(i+1+(2*n)  /7)
                crv3=crv3+dEdl(i+1+2*(2*n)/7)
                crv4=crv4+dEdl(i+1+3*(2*n)/7)
                crv5=crv3+dEdl(i+1+4*(2*n)/7)
                crv6=crv6+dEdl(i+1+5*(2*n)/7)
                crv7=crv7+dEdl(i+1+6*(2*n)/7)
             end do
             if( crv1 > crv2 .and. crv2 > crv3 .and. &
                 crv3 > crv4 .and. crv4 < crv5 .and. &
                 crv5 < crv6 .and. crv6 < crv7 ) bool=.true.
          end if
          
          if(bool) then           
             call stash_add(crv%l  ,l(n)  )
             call stash_add(crv%reE,reE(n))
             call stash_add(crv%imE,imE(n))
          end if

          ! test the curvature peak

!          bool=.true.
!          do i = 0, n - 20  
!             bool = bool .and. curve(odd_n-i).lt.curve(odd_n-1-i) &
!                  .and.curve(2+i).gt.curve(3+i)
!          end do
          bool = .false.
          crvtr1=0.0D0
          crvtr2=0.0D0
          crvtr3=0.0D0
          crvtr4=0.0D0
          crvtr5=0.0D0
          crvtr6=0.0D0
          crvtr7=0.0D0
          do i=1,odd_n/21
             crvtr1=crvtr1+curve(i+6*odd_n/21)
             crvtr2=crvtr2+curve(i+7*odd_n/21)
             crvtr3=crvtr3+curve(i+8*odd_n/21)
             crvtr4=crvtr4+curve(i+9*odd_n/21)
             crvtr5=crvtr5+curve(i+10*odd_n/21)
             crvtr6=crvtr6+curve(i+11*odd_n/21)
             crvtr7=crvtr7+curve(i+12*odd_n/21)
          end do
          if( crvtr1 < crvtr2 .and. crvtr2 < crvtr3 .and.&
               crvtr3 < crvtr4 .and. crvtr4 > crvtr5 .and.&
               crvtr5 > crvtr6 .and. crvtr6 > crvtr7) bool=.true.
          if(.not.bool) then            
             crvtr1=0.0D0
             crvtr2=0.0D0
             crvtr3=0.0D0
             crvtr4=0.0D0
             crvtr5=0.0D0
             crvtr6=0.0D0
             crvtr7=0.0D0
             do i=1,odd_n/7
                crvtr1=crvtr1+curve(i)
                crvtr2=crvtr2+curve(i+odd_n/7)
                crvtr3=crvtr3+curve(i+2*odd_n/7)
                crvtr4=crvtr4+curve(i+3*odd_n/7)
                crvtr5=crvtr5+curve(i+4*odd_n/7)
                crvtr6=crvtr6+curve(i+5*odd_n/7)
                crvtr7=crvtr7+curve(i+6*odd_n/7)
             end do
          if( crvtr1 < crvtr2 .and. crvtr2 < crvtr3 .and.&
              crvtr3 < crvtr4 .and. crvtr4 > crvtr5 .and.&
              crvtr5 > crvtr6 .and. crvtr6 > crvtr7) bool=.true.
          end if
          if (bool) then
             call stash_add(crvtr%l  ,l(n)  )
             call stash_add(crvtr%reE,reE(n))
             call stash_add(crvtr%imE,imE(n))
          end if

          !spot the maxima/minima in the complex part
          bool=.true.
          do i = n-n/10, n - 2  
             bool = bool .and.(    (dimE(odd_n-i).gt.0.0d0.and.dimE(2+i).lt.0.0d0)&
                  .or.&
                  (dimE(odd_n-i).lt.0.0d0.and.dimE(2+i).gt.0.0d0)    )
          end do

          if(bool) then           
             call stash_add(imEmin%l  ,l(n)  )
             call stash_add(imEmin%reE,reE(n))
             call stash_add(imEmin%imE,imE(n))
          end if



       end do   


       !setup the relative tolerance for spotting resonances:
       retol = in%res%re_thres
       imtol = in%res%im_thres*abs(immax-immin)/100.0d0
       if (imtol.gt.abs(in%res%resupcut-in%res%reslocut)) imtol = abs(in%res%resupcut-in%res%reslocut)


       !filter out resonances
       if(crv%l%size.gt.0 .and. crvtr%l%size.gt.0 .and. imEmin%l%size.gt.0) then
          !if(allocated(imEmin)) then

          test_imEmin: do k=1,imEmin%l%size

             if(abs(imEmin%imE%value(k)).lt.in%res%resupcut &
                  .and. abs(imEmin%imE%value(k)).gt.in%res%reslocut) then

                do i = 1, crv%l%size

                   if(abs( abs( crv%imE%value(i) ) - abs( imEmin%imE%value(k) ) ) .lt.imtol &
                        .and. abs( abs(crv%reE%value(i)) - abs(imEmin%reE%value(k))) .lt. retol) then

                      do j = 1, crvtr%l%size

                         if(abs(abs(crvtr%imE%value(j))-abs(imEmin%imE%value(k))) .lt. imtol.and.&
                              abs(abs(crvtr%reE%value(j))-abs(imEmin%reE%value(k))) .lt. retol.and.&
                              abs(abs(crv%imE%value(i))-abs(crvtr%imE%value(j))) .lt. imtol.and.&
                              abs(abs(crv%reE%value(i))-abs(crvtr%reE%value(j))) .lt. retol ) then


                            write(unit=resunit,fmt="(I4.4,2X,F15.7,2x,G17.8,2(2X,D9.2),2(2X,D12.5))") &
                                 state,crv%reE%value(i)-gstate*eh_to_cm1,crv%imE%value(i) &
                                 ,retol &
                                 ,imtol &
                                 ,crv%l%value(i),crvtr%l%value(j)
                            exit test_imEmin

                         end if
                      end do
                   end if
                end do
             end if
          end do test_imEmin
       end if

       call stash_destroy(crv%l     )
       call stash_destroy(crv%reE   )
       call stash_destroy(crv%imE   )
       call stash_destroy(crvtr%l   )
       call stash_destroy(crvtr%reE )
       call stash_destroy(crvtr%imE )
       call stash_destroy(imEmin%l  )
       call stash_destroy(imEmin%reE)
       call stash_destroy(imEmin%imE)

    end do outer_step_loop



  end subroutine res_search_deltar

  subroutine deltaE(reE1,reE2,imE1,imE2,d) 
    implicit none
    real(8) :: reE1,reE2,imE1,imE2,d
    d=( dsqrt( (reE2-reE1)*(reE2-reE1) + (imE2-imE1)*(imE2-imE1) ) ) 
  end subroutine deltaE

  subroutine lambda(h,chi,i,l) 
    implicit none
    real(8) :: conv=219474.63067
    integer :: i
    real(8) :: l,h,chi
    l=h*conv*(chi**dble(i)-1.0d0)/(chi-1.0d0) 
  end subroutine lambda

  subroutine phase(ree1,ree2,ime1,ime2,ph)
    implicit none
    real(8),intent(in) :: ree1,ree2,ime1,ime2
    real(8),intent(out) :: ph
    ph=atan((ime1-ime2)/(ree1-ree2))
  end subroutine phase

  integer function odd(i)
    !returns the ith odd number
    integer :: i
    odd=2*i+1
  end function odd

end module class_res_search
