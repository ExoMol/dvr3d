
module class_cap_input
  use class_error
  implicit none
  private
  type, public :: cap_input
     character(len=7) :: molecule  ! name of the molecule (for labeling purposes only)
     character(len=20) :: capchoice! function that defines the CAP shape
     character(len=4) :: coord     ! string containing indicating CAPed channels (check case)
     character(len=5) :: rmaxchoice  ! string indicating rmax choice,'long' or 'short' or 'adapt'
     real(8) :: drmax,drmin        ! maximum and minimum space covered by the CAP (all channels)
     integer :: ndr                ! number of intermediate distances that need to be checked
     real(8) :: lowest_energy      ! BAND ORIGINS first state to build CAP matrix in cm-1 (inclusive)
     real(8) :: highest_energy     ! BAND ORIGINS last state to  build CAP matrix in cm-1(inclusive)
     real(8) :: ground_energy     ! BAND ORIGINS last state to  build CAP matrix in cm-1(inclusive)
     character(len=200) :: wavefile! files from DVR3D (wavefunction)
     real(8) :: mem                ! available memory in megabytes
  end type cap_input

  public :: read_cap_input
  public :: check_cap_input

contains

  subroutine read_cap_input(c)
    implicit none
    type(cap_input), intent(out) :: c
    write(6,"(20X,'----- READING RESONANCE INPUT ----')")
    read(5,*) !comment line
    read(5,*) c%molecule
    read(5,*) c%capchoice                     ! function that defines the CAP shape
    read(5,*) c%coord                         ! string containing indicating CAPed channels (check case)
    read(5,*) c%rmaxchoice
    read(5,*) c%drmin,c%drmax,c%ndr           ! units of drmax/min a0
    read(5,*) c%lowest_energy,c%highest_energy! lowest and highest energy of interest (cm -1)
    read(5,*) c%ground_energy
    read(5,*) c%wavefile                      ! files from DVR3D (wavefunction)
    read(5,*) c%mem
  end subroutine read_cap_input

  subroutine check_cap_input(c)
    implicit none
    type(cap_input), intent(in) :: c
    write(6,*)
    write(6,"(10X,'Input file was read as:')")!comment line
    write(6,*) !comment line
    write(6,"(10x,'Molecule label = ',A7)") c%molecule
    write(6,"(10x,'CAP shape      = ',A20)") c%capchoice  ! function that defines the CAP shape
    write(6,"(10x,'CAP channel    = ',A4)") c%coord       ! string containing indicating CAPed channels (check case)
    capcase:  select case(c%capchoice)
    case("simple")
       write(6,"(A)",ADVANCE="yes") "           =>  Using Jacobi Simple CAP - dissociation of"
    case("simple3")
       write(6,"(A)",ADVANCE="yes") "           =>  Using Jacobi Simple (3rd power) CAP  - dissociation of"
    case("simple8")
       write(6,"(A)",ADVANCE="yes") "           =>  Using Jacobi Simple (8th power) CAP  - dissociation of"
    case("manop")
       write(6,"(A)",ADVANCE='yes') "           =>  Using Jacobi Manolopoulos CAP - dissociation of"
    case default
       call error("CHECK_CAP_INPUT:unrecognized CAP form... aborting")
    end select capcase
    coords: select case(c%coord)
    case("1-23")
       write(6,*) "              atom 1 from diatom 23"
    case("2-31")
       write(6,*) "              atom 2 from diatom 31"
    case("3-12")
       write(6,*) "              atom 3 from diatom 12"
    case("12-3")
       write(6,*) "              atom 1 from diatom 23 and atom 2 from diatom 31 (c2v molecules)"
    case("123")
       write(6,*) "              atom 1, atom 2 and atom 3 "
    case default
       call error( 'CHECK_CAP_INPUT:unrecognized coord... aborting' )
    end select coords
    write(6,"(10x,'CAPs end - rmax - position, relative to each channels grid  = ',A5)") c%rmaxchoice       ! string containing indicating CAPed channels (check case)
    write(6,"(10x,'CAP length ranging from ',f5.2,' to ',f5.2,' au in ',I2,' steps')") c%drmin,c%drmax,c%ndr
    write(6,"(10x,'CAP Energy window from ',f9.2,' to ',f9.2,' cm^-1 ')") c%lowest_energy,c%highest_energy           ! last state to  build CAP matrix (inclusive)
    write(6,"(10x,'with potential equilibrium at ',f9.2,' cm^-1 ')") c%ground_energy
    write(6,"(10x,'Wave function filename is: ',A20)") trim(c%wavefile)  ! files from DVR3D (wavefunction)
    write(6,"(10x,'Available memory defined as: ',f10.1,' MB')") c%mem
    if(c%mem.le.0) then
       write(6,"(10x,'   => calculation will be fully done in memory')")
    elseif(c%mem.gt.0) then
       write(6,"(10x,'   => cap matrix will be built using kblock scratch')")
    end if
    write(6,*)
    if(c%drmax.le.c%drmin)                  call error('CHECK_CAP_INPUT:c%drmax.le.c%drmin')
    if(c%ndr.le.0)                          call error('CHECK_CAP_INPUT:c%ndr.le.0')
    if(c%lowest_energy.ge.c%highest_energy) call error('CHECK_CAP_INPUT:c%lowest_energy.ge.c%highest_energy')
    if(c%lowest_energy.lt.c%ground_energy)  call error('CHECK_CAP_INPUT:c%lowest_energy.lt.c%ground_energy')

  end subroutine check_cap_input

end module class_cap_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module class_diag_input
  use class_error
  implicit none
  private
  type, public :: diag_input
     real(8) :: lowest_energy,highest_energy ! BAND ORIGIN
     real(8) :: chi,h           ! lambda exponential function progression parameters
     integer :: nlambda         ! number of lambda points to build a trajectory
     logical :: outputvec       ! integer that is set either to T/F indicating if eigevectors from diag are stored
     integer :: ccnumber        ! maximum???? number of crossings (temporary parameter)
     real(8) :: ecthres         ! threshold for correlating eigenvectors while tracking eigenvectors in diag
  end type diag_input

  public :: read_diag_input
  public :: check_diag_input

contains

  subroutine read_diag_input(diag)
    type(diag_input),intent(out) :: diag
    read(5,*)  !this skips the comment line
    read(5,*)  diag%lowest_energy,diag%highest_energy
    read(5,*)  diag%chi, diag%h, diag%nlambda ! should become optional if manop is chosen
    read(5,*)  diag%ccnumber, diag%ecthres
    read(5,*)  diag%outputvec
  end subroutine read_diag_input

  subroutine check_diag_input(d)
    implicit none
    type(diag_input), intent(in) :: d
    real(8), parameter :: eh_to_cm1=219474.63067d0
    write(6,*)  !this skips the comment line
    write(6,"(10x,'Diagonalisation energy window from ',f9.2,' to ',f9.2,' cm^-1 ')")  d%lowest_energy,d%highest_energy
    write(6,"(10x,'Lambda function parameters, Xi = ',f5.2,&
    ' and h = ',E9.2,' calculated at ',I4,' points')")  d%chi, d%h, d%nlambda
    write(6,*)'           => Lambda range (converted to cm-1 from Eh) ....'
    write(6,*)'              l(2). ..........',d%h*eh_to_cm1
    write(6,*)'              N ..............',d%nlambda
    write(6,*)'              l(NLAMBDA) .....',d%h*eh_to_cm1*(d%chi**(d%nlambda-1)-1.d0)/(d%chi-1.d0)

    write(6,"(10x,'Max. num of crossings allowed is ',I3,' with eigenvector correlation min. = ',F3.1)")  d%ccnumber, d%ecthres
    write(6,"(10X,'Vector output is set to ',L)")  d%outputvec

    if(d%lowest_energy.ge.d%highest_energy) call error('CHECK_DIAG_INPUT:d%lowest_energy.ge.d%highest_energy')

  end subroutine check_diag_input


end module class_diag_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module class_res_input
use class_error
  implicit none
private
  type, public :: res_input
     real(8) :: re_thres       !real part of energy window for all criteria
     real(8) :: im_thres       !imaginary part of energy window for all criteria
     integer :: resmode         ! number of points in peak >= 3 !!
     real(8) :: resupcut          !resonance maximum width value cutoff (wavenumbers)
     real(8) :: reslocut          !resonance minimum width value cutoff (wavenumbers)
  end type res_input

public :: read_res_input
public :: check_res_input

contains

  subroutine read_res_input(r)
    type(res_input), intent(out) :: r
    read(5,*)  !this skips the comment line
    read(5,*)  r%resmode
    read(5,*)  r%resupcut, r%reslocut
    read(5,*)  r%re_thres, r%im_thres
  end subroutine read_res_input

  subroutine check_res_input(r)
    type(res_input), intent(in) :: r
    write(6,*)  !this skips the comment line
    write(6,"(10x,'Number of eigenvalues used to find resonant point = ',I2)")  r%resmode
    write(6,"(10x,'Resonance half width window from ',D9.3,' to ',D9.3,' cm^-1 ')")  r%resupcut, r%reslocut
    write(6,"(10x,'Real error threshold = ',D9.3)")  r%re_thres
    write(6,"(10x,'Imaginary error threshold = ',D9.3)")  r%im_thres
  if(r%resupcut.le.r%reslocut)call error('CLASS_RES_INPUT:r%resupcut.le.r%reslocut')
  end subroutine check_res_input

end module class_res_input

!!!!!!!!!!!!!!!!!!!!!!!!
!Inheriting classes
!!!!!!!!!!!!!!!!!!!!!!!!
module class_input
  use class_error
  use class_cap_input
  use class_diag_input
  use class_res_input
  implicit none
  type, public :: input
     type (cap_input)  :: cap
     type (diag_input) :: diag
     type (res_input) :: res
  end type input

!public :: cap_input
! inherited type
!public :: diag_input
! inherited type

public :: read_input
!  subroutine read_input(data)
!    !this subroutine reads the data from the input file
!    !and dumps it into the object of type input_type
!    type (input),intent(out) :: data
!  end subroutine read_input
contains
  subroutine read_input(data)
    !this subroutine reads the data from the input file
    !and dumps it into the object of type input_type
    implicit none
    integer :: idia,ipar,max_ang_bas,nmaxr1,nmaxr2,jrot,kmin1,neval
    type (input),intent(out) :: data
    type (cap_input)  :: cap
    type (diag_input) :: diag
    type (res_input)  :: res
    call read_cap_input(cap)
    call check_cap_input(cap)
    call read_diag_input(diag)
    call check_diag_input(diag)
    call read_res_input(res)
    call check_res_input(res)
    data%cap  = cap
    data%diag = diag
    data%res  = res
!DO A FINAL CAP AND DIAG ENERGY INPUT CHECK:
if (cap%highest_energy .lt. diag%highest_energy ) call error('READ_INPUT:cap%highest_energy .lt. diag%highest_energy')
if (cap%lowest_energy .gt. diag%lowest_energy ) call error('READ_INPUT:cap%lowest_energy .gt. diag%lowest_energy')
  end subroutine read_input



end module class_input

module common_parameters
  use class_error
  use class_input
  use class_kblock
  use class_wavefile
  use class_wavefile_to_kblock
  use class_dvr_grid
  implicit none
  public
  save
  integer :: idia
  integer :: ndeltar
  real(8), pointer :: dr
  real(8), allocatable, target :: deltar(:)
  real(8), allocatable:: deltar_small(:)
  logical,allocatable :: deltar_choice(:)
  real(8) :: rmax(3)
  character(len=4) :: coord
  character(len=20) :: capchoice
  character(len=200) :: filelabel
  integer :: cap_low_state
  integer :: cap_high_state
  integer :: n_cap_states
  integer :: diag_low_state
  integer :: diag_high_state
  integer :: n_diag_states
  integer :: nlambda
  real(8) :: gstate
  real(8), parameter :: eh_to_cm1=219474.63067d0
  real(8) :: mem
contains

  subroutine common_parameters_from_kblock_file(in,block,uni)
    implicit none
    type(input) :: in
    type(kblock), pointer :: block
    integer,optional :: uni
    integer :: unit
    real(8) :: eh_min , eh_max
    integer :: k,i,j
    integer :: lowest_state,highest_state
    character(len=200) :: wavefilename
    real(8) :: tmp1,tmp2,tmp3
    if(.not.present(uni)) then
       unit=6
    else
       unit=uni
    end if
    idia=block%grid%idia
    coord=in%cap%coord
    capchoice=in%cap%capchoice
    WRITE(UNIT,*)
    write(unit,"(10x,'Setting up the cap matrix energy window...')")
    eh_min=in%cap%lowest_energy/eh_to_cm1
    eh_max=in%cap%highest_energy/eh_to_cm1
    gstate=in%cap%ground_energy/eh_to_cm1
    loop_cap_min: do i=1,block%neval
       lowest_state=i
       !       if((block%eh(i)-block%eh(1)).ge.eh_min) exit loop_cap_min
       if((block%eh(i)-gstate).ge.eh_min) exit loop_cap_min
    end do loop_cap_min
    loop_cap_max: do i=block%neval,1,-1
       highest_state=i
       !       if((block%eh(i)-block%eh(1)).le.eh_max) exit loop_cap_max
       if((block%eh(i)-gstate).le.eh_max) exit loop_cap_max
    end do loop_cap_max

    if(highest_state.le.lowest_state)then
       write(unit,*)eh_min+gstate,eh_max+gstate,lowest_state,highest_state
       call error('CAP_MATRIX_SET_SIZE: highest_state.le.lowest_state')
    end if
    cap_low_state  = lowest_state
    cap_high_state = highest_state
    n_cap_states=cap_high_state-cap_low_state+1
    write(unit,*)
    write(unit,"(10x,'ground state energy     = ',F8.1,'cm^-1')")gstate*eh_to_cm1
    write(unit,"(10x,'lowest energy selected  = ',F8.1,'cm^-1')")block%eh(cap_low_state)*eh_to_cm1
    write(unit,"(10x,'highest energy selected = ',F8.1,'cm^-1')")block%eh(cap_high_state)*eh_to_cm1
    write(unit,"(10x,'total number of states  = ',I4)"          )n_cap_states
    write(unit,"(10x,'   starting at            ',I4)"          )cap_low_state
    write(unit,"(10x,'   ending at              ',I4)"          )cap_high_state
    write(unit,*)
    write(unit,"(10x,'Setting up the matrix diagonalisation energy window...')")
    eh_min=in%diag%lowest_energy/eh_to_cm1
    eh_max=in%diag%highest_energy/eh_to_cm1
    loop_diag_min: do i=1,block%neval
       lowest_state=i
       !       if((block%eh(i)-block%eh(1)).ge.eh_min) exit loop_diag_min
       if((block%eh(i)-gstate).ge.eh_min) exit loop_diag_min

    end do loop_diag_min
    loop_diag_max: do i=block%neval,1,-1
       highest_state=i
       !       if((block%eh(i)-block%eh(1)).le.eh_max) exit loop_diag_max
       if((block%eh(i)-gstate).le.eh_max) exit loop_diag_max
    end do loop_diag_max

    if(highest_state.le.lowest_state)then
       write(unit,*)eh_min+gstate,eh_max+gstate,lowest_state,highest_state
       call error('CAP_MATRIX_SET_SIZE: highest_state.le.lowest_state')
    end if
    diag_low_state  = lowest_state
    diag_high_state = highest_state
    n_diag_states=diag_high_state-diag_low_state+1
    write(unit,*)
    write(unit,"(10x,'ground state energy     = ',F8.1,'cm^-1')")gstate*eh_to_cm1
    write(unit,"(10x,'lowest energy selected  = ',F8.1,'cm^-1')")block%eh(diag_low_state)*eh_to_cm1
    write(unit,"(10x,'highest energy selected = ',F8.1,'cm^-1')")block%eh(diag_high_state)*eh_to_cm1
    write(unit,"(10x,'total number of states  = ',I4)"          )n_diag_states
    write(unit,"(10x,'   starting at            ',I4)"          )diag_low_state
    write(unit,"(10x,'   ending at              ',I4)"          )diag_high_state

    if(debug) then
       write(unit,*)'**  DEBUG INFO  **'
       write(unit,*)'common_parameters_from_kblock_file:'
       write(unit,*)'cap_low_state  ', cap_low_state
       write(unit,*)'cap_high_state ', cap_high_state
       write(unit,*)'diag_low_state ', diag_low_state
       write(unit,*)'diag_high_state', diag_high_state
       write(unit,*)'** END DEBUG INFO **'
    end if
    write(unit,*)
    write(unit,"(10x,'Setting up the cap deltar values...')")
    ndeltar=in%cap%ndr
    if (.not.allocated(deltar_choice))allocate( deltar_choice(ndeltar) )
    WRITE(UNIT,*)
    write(unit,"(10x,'The selected ',i2,' CAP ranges are: ')")ndeltar
    if(.not.allocated(deltar))allocate(deltar(ndeltar))
    if (ndeltar.eq.1) then
       deltar(1)=in%cap%drmin
       write(unit,"(10x,f5.2,' au')")deltar(1)
    elseif (ndeltar.gt.1) then
       do i=1,ndeltar
          deltar(i)=dble(i-1)*(in%cap%drmax-in%cap%drmin)/dble(ndeltar-1)+in%cap%drmin
          write(unit,"(10x,f5.2,' au')")deltar(i)
       end do
    end if
    !write(unit,*)
    !write(unit,"(10x,'Setting up rmax for all dissociation channels...')")
    rmax(1)=dvr_grid_r2_max(block%grid,'1-23')
    rmax(2)=dvr_grid_r2_max(block%grid,'2-31')
    rmax(3)=dvr_grid_r2_max(block%grid,'3-12')
    write(unit,"(10x,'for kz = ',I2,' :')")block%k
    write(unit,"(10x,'rmax(1) for 1-23 = ',f5.2,' au')")rmax(1)
    write(unit,"(10x,'rmax(2) for 2-31 = ',f5.2,' au')")rmax(2)
    write(unit,"(10x,'rmax(3) for 3-12 = ',f5.2,' au')")rmax(3)
    do i=1,3
       do j=1,ndeltar
          if(rmax(i).le.deltar(j))call error('COMMON_PARAMETERS_FROM_KBLOCK_FILE:rmax(i).le.deltar(j) -- review dump file')
       end do
    end do
    mem = in%cap%mem
    RETURN
  end subroutine common_parameters_from_kblock_file

  subroutine kblock_head_dump(in,kz)
    implicit none
    type(input) ::  in
    type(kblock), pointer :: block
    integer :: kz,i
    character(len=3) :: k_char
    character(len=200) :: filename
    logical :: exists
    call kblock_head_read(kz,block)
    write( k_char ,fmt="(I3.3)" ) block%k
    filename='kblock'//k_char//'_head_dump.bcs'
    inquire(file=trim(filename),exist=exists)
    if(.not.exists) then
       write(6,*)
       write(6,"(10x,'Printing kblock dump file to help setup input file...')")

       open(unit=100            &
            ,file=trim(filename)&
            ,status='new'       &
            ,form='formatted')
       write(unit=100,fmt=*) 'block%k    ',block%k     !(1)  kblock label
       write(unit=100,fmt=*) 'block%jk   ',block%jk    !totalnumber of kblocks of which this one is part.
       write(unit=100,fmt=*) 'block%kpar ',block%kpar  !(1)  kblock parity (postprocessed)
       write(unit=100,fmt=*) 'block%ipar ',block%ipar  !(1)
       write(unit=100,fmt=*) 'block%jrot ',block%jrot  !(1)   total angular momentum
       write(unit=100,fmt=*) 'block%neval',block%neval !(1)  number of eigenvalues
       write(unit=100,fmt=*) 'block%nbass',block%nbass ! npnt1*(npnt2-2*kpar)*
       call dvr_grid_dump(100,block%grid)
       write(unit=100,fmt=*) 'block%eh(:)*eh_to_cm1'
       do i=1,block%neval
          write(unit=100,fmt=*) i, block%eh(i)*eh_to_cm1
       end do
       call common_parameters_from_kblock_file(in,block,100)
       close(100)
       call kblock_head_destroy(block)
       write(6,*)
       write(6,"(10x,'wrote ',A,' to disk')")trim(filename)
       write(6,"(10x,'review the file to adjust input')")
       write(6,*)
       write(6,"(10x,'Re-run resonance to start the calculation')")
       write(6,"(10x,'(without deleting the file)')")
       stop
    else
       write(6,"(10x,A,' already exists ')")trim(filename)
       write(6,"(10x,' -- proceeding with calculation')")
    end if
    return
  end subroutine kblock_head_dump

end module common_parameters







