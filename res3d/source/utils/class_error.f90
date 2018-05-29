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
