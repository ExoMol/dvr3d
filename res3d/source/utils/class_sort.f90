module class_sort
implicit none
private
interface sort
module procedure sort_int
module procedure sort_real
end interface

public :: sort
contains
subroutine sort_real(a,ind,l,r)
integer,intent(in) :: l,r !l = leftmost index, r = righmost index
real(8),intent(inout) :: a(r-l+1)
integer,intent(inout) :: ind(r-l+1)
call quicksort(a,ind,l,r) 
call insertionsort(a,ind,l,r)
end subroutine sort_real

subroutine sort_int(a,ind,l,r)
integer,intent(in) :: l,r !l = leftmost index, r = righmost index
integer,intent(inout) :: a(r-l+1)
integer,intent(inout) :: ind(r-l+1)
call quicksort_int(a,ind,l,r) 
call insertionsort_int(a,ind,l,r)
end subroutine sort_int


recursive subroutine quicksort(a,ind,l,r)
  integer,intent(in) :: l,r !l = leftmost index, r = righmost index
  real(8),intent(inout) :: a(:)
  integer,intent(inout) :: ind(:)
  
  integer :: M=4,i,j,v

  if ((r-l)>M) then 
     i = (r+l)/2
     if (a(l)>a(i)) then ;call swap_real(a,l,i);call swap_int(ind,l,i);end if
     if (a(l)>a(r)) then ;call swap_real(a,l,r);call swap_int(ind,l,r);end if  
     if (a(i)>a(r)) then ;call swap_real(a,i,r);call swap_int(ind,i,r);end if
     j = r-1
     call swap_real(a,i,j)
     call swap_int(ind,i,j)
     i = l
     v = a(j)
     outer: do
        inner1:do 
           i=i+1 ; if(i == 0) exit inner1
           if(a(i).ge.v) exit inner1
        end do inner1
        inner2:do 
           j=j-1 ; if(j == 0) exit inner2
           if(a(j).le.v) exit inner2
        end do inner2
        if (j.lt.i) exit outer
     end do outer
     call swap_real(a,i,r-1) ; call swap_int(ind,i,r-1)
     call quicksort(a,ind,l,j)
     call quicksort(a,ind,i+1,r)
  end if
end subroutine quicksort

recursive subroutine quicksort_int(a,ind,l,r)
  integer,intent(in) :: l,r !l = leftmost index, r = righmost index
  integer,intent(inout) :: a(:)
  integer,intent(inout) :: ind(:)
  
  integer :: M=4,i,j,v

  if ((r-l)>M) then 
     i = (r+l)/2
     if (a(l)>a(i)) then ;call swap_int(a,l,i);call swap_int(ind,l,i);end if
     if (a(l)>a(r)) then ;call swap_int(a,l,r);call swap_int(ind,l,r);end if  
     if (a(i)>a(r)) then ;call swap_int(a,i,r);call swap_int(ind,i,r);end if
     j = r-1
     call swap_int(a,i,j)
     call swap_int(ind,i,j)
     i = l
     v = a(j)
     outer: do
        inner1:do 
           i=i+1 ; if(a(i).ge.v) exit inner1
        end do inner1
        inner2:do 
           j=j-1 ; if(a(j).le.v) exit inner2
        end do inner2
        if (j.lt.i) exit outer
     end do outer
     call swap_int(a,i,r-1) ; call swap_int(ind,i,r-1)
     call quicksort_int(a,ind,l,j)
     call quicksort_int(a,ind,i+1,r)
  end if
end subroutine quicksort_int


subroutine swap_real(a,i,j)
real(8),intent(inout) :: a(:)
integer,intent(in) :: i,j 

real(8) :: T
T=a(i)
a(i)=a(j)
a(j)=T
end subroutine swap_real

subroutine swap_int(a,i,j)
integer,intent(inout) :: a(:)
integer,intent(in) :: i,j 

integer :: T
T=a(i)
a(i)=a(j)
a(j)=T
end subroutine swap_int

subroutine insertionsort(a,ind,lo0,hi0)
  integer,intent(in) :: lo0, hi0
  real(8),intent(inout) :: a(:)
  integer,intent(inout) :: ind(:)
  integer :: i,j,v_int
  real(8) :: v
  
  do i=lo0+1,hi0
 
     v=a(i)
     v_int=ind(i)
     j=i

     loop: do while(j.gt.lo0.and.a(j-1).gt.v)
        a(j)=a(j-1)
        ind(j)=ind(j-1)
        j=j-1
        if(j==lo0) exit loop ! this makes sure that a(j-1) is not evaluated when j=lo0
     end do loop
     a(j)=v
     ind(j)=v_int
  end do

end subroutine insertionsort

subroutine insertionsort_int(a,ind,lo0,hi0)
  integer,intent(in) :: lo0, hi0
  integer,intent(inout) :: a(:)
  integer,intent(inout) :: ind(:)
  integer :: i,j,v_int
  real(8) :: v
  
  do i=lo0+1,hi0
 
     v=a(i)
     v_int=ind(i)
     j=i

     loop: do while(j.gt.lo0.and.a(j-1).gt.v)
        a(j)=a(j-1)
        ind(j)=ind(j-1)
        j=j-1
        if(j==lo0) exit loop ! this makes sure that a(j-1) is not evaluated when j=lo0
     end do loop
     a(j)=v
     ind(j)=v_int
  end do

end subroutine insertionsort_int

end module class_sort


!UNCOMMENT TO TEST. COMPILE DOING:
!ifort -check all -traceback class_sort.f90

!program sort_test
!use class_sort
!implicit none
!real(8) :: array(20),tmp(20)
!integer :: l=1,r=20,i
!integer :: indexes(20)
!
!
!array(1)  =234.11d0
!array(2)  =324.11d0
!array(3)  =4.11d0
!array(4)  =6.11d0
!array(5)  =7.11d0
!array(6)  =24.11d0
!array(7)  =783.11d0
!array(8)  =9876.11d0
!array(9)  =567.11d0
!array(10) =57.11d0
!array(11) =36.11d0
!array(12) =45.11d0
!array(13) =27.11d0
!array(14) =83.11d0
!array(15) =552.11d0
!array(16) =786.11d0
!array(17) =628.11d0
!array(18) =5226.11d0
!array(19) =45.11d0
!array(20) =8.11d0
!
!tmp = array
!
!do i=1,20
!indexes(i)=i
!end do
!
!do i=1,20
!write(6,*) indexes(i),array(i)
!end do
!
!
!call sort(array,indexes,l,r)
!
!write(6,*) 'sorted array:'
!do i=1,20
!write(6,*) indexes(i),array(i)
!end do
!write(6,*) 'array from sorted indexes'
!do i=1,20
!write(6,*) indexes(i),tmp(indexes(i))
!end do
!
!
!end program sort_test



