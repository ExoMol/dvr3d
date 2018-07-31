
  module fit_module


    implicit none
    !
    integer, parameter    :: verbose  = 3 ! Verbosity level    
    !
    public fitting
    !
    private
    !
    integer, parameter :: enermax=6000    ! maximal number of calculated energies for a given J and a given symmetry 
    !integer, parameter :: obsmax=5000     ! maximal number of experimental energies
    integer, parameter :: pointsmax=7000  ! maximal number of potential data points 
    !
    integer, parameter :: f_inp = 5, f_out = 6   ! read/write units 
    !
    double precision,save :: pi                  ! the pi number 
    !
    integer,parameter          ::  findif = 3    ! Dinite difference differentiation with 2 or 3 points
    !
    double precision,parameter :: stadev_best=1e-04,stab_best=1e-12  ! best standard error and srability 
                                                 ! the fit will finish when these reached. 
    !
    double precision,parameter :: fitfactordeltax=0.001   ! parameter for the finite differncies differention 
    !
    character(len=70) :: deriv_type = 'hellman  '  ! hellman or direct - how we calculate derivatives. 
                                                 ! direct means the finite diferentiation of energies wrt parameters.
    !
    integer,parameter :: nofititer=6             ! the Jacobi matrix can be evaluated only each nofititer+1 step. For all 
                                                 ! iterations in between the the same jacobi matrix is used 
    integer,parameter :: f_en =206               ! All computed energies (not only those that have obs. counterparts) will be written into here.
    integer,parameter :: f_pot=207               ! The current potential parameters will be written here. 
    integer,parameter :: f_potpoints =208        ! Computed at the current iteration step potential energy points are printed here. 
    integer,parameter :: f_res=209               ! it is where we write the fitting results 
    character(len=200)::  char_job(10)           ! The derectives from the input file are stored in this array. Will be used to create all job-input files.
    character(len=200)::  char_vib_job(100)      ! The vib. derectives from the input file are stored in this array. Will be used to create all job-input files.
    character(len=200)::  char_rot_job(100,2)      ! The rot. derectives from the input file are stored in this array. Will be used to create all job-input files.
    character(len=1)  :: mark
    !
    integer,parameter :: sym2parity(4) = (/0,1,0,1/) ! parity versus symmetry
    !
    ! global variables, as they appear in the DVR3D input file.
    !
    integer ::   npnt2_,jrot_,neval_,nalf_,max2d_,max3d_,idia_,&
                 kmin_,npnt1_,ipar_,max3d2_
    integer ::   ibass_ = 0, nvib_, neval2_ = 0
    integer      :: npnt2x(100),nevalx(100),jrot_card(100),icard_j(0:100)
    !
    integer ::   jmax = 10    ! maximal value up to whicj all j = 0,1,2,..jmax will be evaluated. 
    !
    real(8) ::   ezero_ = -1                     ! zero point energy as computed at J=0,A1 
    !
    ! type to deal with the calculated DVR3D energies 
    !
    type  FTdvr3dT
      real(8),pointer   ::  energy(:)   ! to store the calculated energies for  given (J,symmetry)
      real(8),pointer   ::  derj(:,:)   ! to store the calculated derivatives 
    end type FTdvr3dT

    !
   
    contains
     !
     subroutine fitting 


      integer           :: parmax                  ! total number of potential parameters 
      integer           :: jrot                    ! current value of the rotational quantum number 
      character(len=6)  :: pes_type   = 'dvr3d'    ! Can be used to switch between potential functiion. 
      !
      character(len=30) :: dvr3d_exe               ! DVR3d exe-file. 
      character(len=30) :: rotlev_exe              ! DVR3d exe-file. 
      character(len=30) :: xpect_exe               ! DVR3d exe-file. 

      integer           :: alloc,info, ierror      ! error state variables 

      integer           :: npts,en_npts,pot_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used. 
      !
      integer           :: nrow,ncol               ! loop-variables
      integer           :: numpar                  ! number of varying parameters 
      integer           :: ndigits                 ! number of rounded digits in the parameters output. 
      integer           :: itmax                   ! Number of fitting iterations.
      !
      ! objects to store the pot. parameters, geometries, values of poten. function.  
      !
      double precision,allocatable :: parold(:),local(:,:),pot_values(:),potparam(:)
      !
      !integer,allocatable:: ipower(:,:)            ! powers of the potential function expansion. Can be used 
      !                                             ! to form the potential function os just as a parameter name.
      character(len=8),allocatable :: nampar(:)    ! parameter names 
      
      integer,allocatable :: ivar(:)               ! switch for the parameters fit: ivar(i) = 0 for the parameters that do not vary.
      !
      integer,allocatable :: J_obs(:),sym_obs(:),N_obs(:),quanta_obs(:,:) ! Arrays where we 
                                                   ! store the information about obs. energies: 
                                                   ! J, symmetry, number in the block, vib. quantum numbers
      !
      integer             :: jrot_list(1:100)      ! list of j-s to be processed. 
      !
      ! objects to store the energies 
      double precision,allocatable :: enercalc(:),ener_obs(:),eps(:)
      character(len=2),allocatable      :: gamma_obs(:)
      !
      type(FTdvr3dT)                    :: dvr3d(4)   ! caclualted term DVR3D values;
      !
      ! objects to needed for the fitting procedure: jacobi matrix, derivatives, matrices to solve the 
      ! the ststem of linear equations, standard errors and so on. 
      double precision,allocatable :: rjacob(:,:),al(:,:),bl(:),dx(:),ai(:,:),sterr(:),derj0(:)
      double precision,allocatable :: derj(:,:,:),sigma(:)

      double precision,allocatable :: wt_bit(:),wtall(:) ! weight factors - only for the energies and total.
      !
      !
      double precision :: wtsum,fit_factor,r1,r2,theta,ssq,rms,sum_sterr,conf_int
      double precision :: stadev_old,stability,stadev,tempx,deltax,v,potright,potleft
      double precision :: ssq1,ssq2,rms1,rms2
      logical          :: still_run
      integer          :: iener1,iener2,j,meval,i,l
      integer          :: iJ,ipar,isym,jsym,iener,irow,icolumn,ncards,icard,nsym(0:1,0:1),nasym(0:1,0:1),iasym
      !
      !     For the details on the robust fit see J.K.G. Watson, JMS, 219, 326 (2003).
      !     The idea is optimize the weigth factors iteratevely depending on the fitting results, as  follows:
      !     w_i = 1 / sigma_i^2+alpha*eps_i^2 
      !     where sigma_i represent the experimental data precision, alpha is the Watson parameter (can be also adjsuted) and 
      !     eps_i = Obs - Calc deviatons. 
      !
      double precision :: robust_fit = 0 
      !                   robust_fit = 0 corresponds to the standard fit, the robust fittting procedure is switched off. 
      !                   robust_fit <> =0, the robust fittting procedure is on.
      !                   robust_fit also defines the target accuracy of the fit. For example for robust_fit  = 0.001 
      !                   all data with (obs-calc) < 0.001 will be considered as relatively bad (i.e. outliers ) 
      !                   and their weights will we reduced. 
      !
      ! different character objects 
      !
      character(len=68)  :: fmt1,fmt2   ! to have nice parameters output   
      character(len=2)   :: fmt0,ch2    
      character(len=1)   :: ch1
      character(len=6)   :: fit_type    ! to switch between fitting methods. 
      character(len=7)   :: out_name    ! file name for outputs
      character(len=7)   :: ch_tmp      ! parameter names 
      character(len=80)  :: char_tmp    ! to read a current line in the input file.
      character(len=200) :: char_job_   
      !
      ! parameters needed for the dgelss lapack routine: 
      !
      double precision, allocatable :: Tsing(:,:)
      integer                       :: rank0,jlistmax,neval_t(4),lwork
      double precision,allocatable  :: wspace(:)
      !
      integer                     :: itime0,itime2,itime,irate2,imax2 ! time control variables 
      !
      integer                     :: fititer   ! iteration variable
      !
      logical                     :: ifopen,isys
      !
      logical                     ::  SYSTEMQQ  ! system function for  calling extyernal programs, 
                                                ! platform dependent.  

       !
       ! START PROGRAM HERE 
       !
       if (verbose>=2) write(f_out,"(/'The fitting procedure starts...'/)")
       !
       ! read the beginning time. 
       !
       call SYSTEM_CLOCK(itime0,irate2,imax2)
       !
       !dir = 'C:\sergei\programs\fiiting_dvr\myprogs\fit\'
       !
       ! pi constant 
       !
       pi = 4.0d0 * atan2(1.0d0,1.0d0)
       !
       ! Read the input file, count all energies, parameters, data points, so that we can 
       ! allocate all needed arrays and only then actuly read the input again. 
       !
       ! First 10 lines  stand for the DVR3D vibrational input (ipar=0), which will be 
       ! used to generate all other input files needed for DVR3D. 
       !
       !
       read(f_inp,"(a200)")  char_job(1)
       read(f_inp,"(a200)")  char_job(2)
       !
       char_job_(1:3)=""
       !
       char_vib_job(:)="xxx"
       !
       ncards = 0
       icard_j(:) = 99
       !
       ! prepare the dvr3d object
       !
       do isym = 1,4
	 dvr3d(isym)%energy => null()
	 dvr3d(isym)%derj => null()
       enddo
       !
       do 
          !
          read(f_inp,"(a200)") char_job_ 
          if ( char_job_(1:3)=="---") exit
          !
          ! Read the enviroment constants
          !
          read(char_job_,"(11i5)") npnt2_,jrot_,neval_,nalf_,max2d_,max3d_,idia_,&
                            kmin_,npnt1_,ipar_,max3d2_

          if (jrot_>=0.and.jrot_<=100) then
             !
             ncards = ncards + 1
             !
             jrot_card(ncards) = jrot_
             icard_j(jrot_) = ncards
             !
             char_vib_job(ncards) = char_job_
             npnt2x(ncards) = npnt2_
             nevalx(ncards) = neval_
             !
             if (ncards==1) then  
               npnt2x(:) = npnt2_
               nevalx(:) = neval_
             endif
             !
             if (ncards==1.and.jrot_/=0) then
                !
                write (f_out,"('illegal value J ',i,' in the firts vib. card (not zero)')") jrot_
                stop 'illegal value J in the first card'
                !
             endif 
             !
          else
             !
             write (f_out,"('illegal value J in input:',i)")  jrot_
             stop 'illegal value J in input'
             !
          endif 
          !
       enddo 

       do i = 3,9
         read(f_inp,"(a200)") char_job(i)
       enddo 
       !
       ! skip a line 
       read(f_inp,"(a80)") char_tmp
       !

       !read(f_inp,"(a200)")  char_job(11)

       icard = 1 
       do 
          !
          icard = icard + 1
          !
          read(f_inp,"(a200)") char_job_ 
          if ( char_job_(1:3)=="---") exit
          !
          jrot_  = jrot_card(icard)
          !
          char_rot_job(icard,1) = char_job_
          !
          read(f_inp,"(a200)") char_job_ 
          if ( char_job_(1:3)=="---") then 
             write (f_out,"('missing second line in the rot. card ',i,' for j = ',i)") icard,jrot_
             stop 'missing second line in the rot. card'
          endif 
          !
          char_rot_job(icard,2) = char_job_
          !
          ! Read the enviroment constants
          !
          read(char_job_,"(11i5)") nvib_,neval_
          !
          if (nvib_/=nevalx(icard)) then
             !
             write (f_out,"('value nvib = ',i,' in the rot. card ',i,' for j = ',i)") nvib_,icard,jrot_
             write (f_out,"('different from nvib = ',i,' from vib. card above')") nevalx(icard)
             stop 'illegal value nvib in rot. card'
             !
          endif 
          !
       enddo 
       !
       ! skip a line 
       read(f_inp,"(a80)") char_tmp
       !
       ! Read the list of j-values to be evaluated. 
       !
       ! Count the totl number of records (j-s) on the line:
       !
       i = 1
       j = 0 
       !
       ch1 = ' '
       !
       do while (i<80.and.ch1/='<'.and.ch1/='=')
         read(char_tmp(i:i),"(a1)") ch1
         i = i + 1 

         do while (ch1==' ')
           read(char_tmp(i:i),"(a1)") ch1
           i = i + 1 
         enddo 

         do while (ch1/=' '.and.ch1/='<'.and.ch1/='=')
           read(char_tmp(i:i),"(a1)") ch1
           read(char_tmp(i-1:i),"(a2)") ch2
           j = j + 1 
           i = i + 1 
           !
           if (ch2(1:1)/=' '.and.ch2(2:2)/=' ') j = j - 1
           !
         enddo 
         !
       enddo 
       !
       i = i-1 
       jlistmax = j
       jrot_list = 0 
       !
       ! if the list is not presented  - computed all J-s upto J=jmax, 
       ! for which prepare the list of J-s. 
       !
       if (i>=80) then
         do  jrot = 0,jmax
           jrot_list(jrot+1) = jrot 
         enddo 
       else
         read(char_tmp(1:i),*) (jrot_list(j),j=1,jlistmax)
         do  j = 1,jlistmax
           jrot = jrot_list(j)
           !
           icard = icard_j(jrot)
           !
           char_job_ = char_vib_job(icard)
           !
           if (char_job_(1:3)=="xxx") then
              !
              write (f_out,"('an input vib. card  for J = ',i,' was not specified')") jrot_
              stop 'input: illegal J in the  list -J values-'
              !
           endif 
           !
         enddo 
       endif
       ! 
       ! Maximal number of iterations.
       !
       read(f_inp,*) itmax
       !
       ! fit_factor defines the ralative importance of the obs. energies vs. the potential energies 
       ! in the simultaneous fit: 
       ! fit_factor = weight_ener/weight_pot
       !
       read(f_inp,*) fit_factor 
       !
       ! Switch on/off the robust fitting
       !
       read(f_inp,*) robust_fit
       !
       ! The fitting methods, e.g. dgelss or linur, or the method to solve the linear system. 
       !
       read(f_inp,"(a6)") fit_type
       !
       ! The name for the output file
       !
       read(f_inp,"(a6)") out_name
       !
       ! Can be used to switch between potential functiion. 
       !
       read(f_inp,"(a6)") pes_type
       !
       ! Read the names of the exe-files: dvr3d, rotlev, and xpect.
       !
       read(f_inp,"(a30)") dvr3d_exe
       read(f_inp,"(a30)") rotlev_exe
       read(f_inp,"(a30)") xpect_exe
       !
       write(f_out,"('itmax = ',i5)") itmax
       write(f_out,"('fit_factor = ',e12.4)") fit_factor
       write(f_out,"('fit_type = ',a6)") trim(fit_type)
       write(f_out,"('pes_type = ',a6)") trim(pes_type)
       if (robust_fit>0) then 
          write(f_out,"('Robust fit will be used with the Watson parameter = ',f14.8)") robust_fit
          !
       end if 
       write(f_out,"('J-s = ',100i)")  (jrot_list(j),j=1,jlistmax)
       write(f_out,"('dvr3d  exe file is ',a30)") dvr3d_exe
       write(f_out,"('rotlev exe file is ',a30)") rotlev_exe
       write(f_out,"('xpect  exe file is ',a30)") xpect_exe
       !
       ! skip a line 
       read(f_inp,"(a80)") char_tmp
       !
       ! Number of the potential parameters 
       read(f_inp,*) parmax
       !
       ! Allocate all arrays that relate to the potential parameters 
       !
       allocate (parold(parmax),potparam(parmax),nampar(parmax),ivar(parmax),&
                 al(parmax,parmax),ai(parmax,parmax),bl(parmax),dx(parmax),sterr(parmax),&
                 Tsing(parmax,parmax),stat=alloc)
       if (alloc/=0) then
         write (f_out,"(' Error ',i,' initializing potparam objects')") alloc
         stop 'potparam - alloc'
       end if
       !
       do i=1,parmax
          !
          read(f_inp,"(a80)") char_tmp
          read(char_tmp(1:8),"(a8)") nampar(i)
          read(char_tmp(9:80),*) ivar(i),potparam(i)
          !write(nampar(i),"('pot_',3i1)") ipower(1:3,i)
          !
       end do
       !
       ! check the the next line in the input, has to be "---"
       !
       read(f_inp,"(a80)") char_tmp
       if (char_tmp(1:3)/="---") then 
         write (f_out,"('Number of parameters do not agree with number of lines :',i10)")  parmax
         stop 'wrong number of parameters or inputlines'
       endif 
       !
       ! Input Obs. energies
       !
       ! Count the obs. term values. 
       !
       npts = 0
       char_tmp(1:3)=""
       do while( char_tmp(1:3)/="---" )
          !
          npts=npts+1
          read (f_inp,"(a80)") char_tmp
          !
       enddo 
       !
       ! Number of energies  
       en_npts = npts-1
       !
       write(f_out,"('Number of obs. data points: ',i9)") en_npts
       !
       ! Allocate the obs. energy related arrays
       !
       allocate (ener_obs(en_npts),enercalc(en_npts),J_obs(en_npts),&
                 sym_obs(en_npts),gamma_obs(en_npts),N_obs(en_npts),&
                 quanta_obs(3,en_npts),stat=alloc)
       if (alloc/=0) then
         write (f_out,"(' Error ',i,' initializing obs. energy related arrays')") alloc
         stop 'obs. energy arrays - alloc'
       end if
       !
       ! Potential energy points 
       !
       ! Count them:
       !
       pot_npts = 0
       char_tmp(1:3)=""
       do while( char_tmp(1:3)/="---" )
         !
         pot_npts=pot_npts+1
         read (f_inp,"(a80)") char_tmp
         !
       enddo 
       !
       pot_npts=pot_npts-1
       !
       write(f_out,"('Number of potential energy data points: ',i9)") pot_npts
       !
       npts = pot_npts + en_npts
       !
       allocate (wtall(npts),wt_bit(npts),&
                local(3,pot_npts),pot_values(pot_npts),stat=alloc)
       !
       if (robust_fit>0) allocate (sigma(npts),stat=alloc)
       !
       if (alloc/=0) then
         write (f_out,"(' Error ',i,' initializing weights and potpoints')") alloc
         stop 'weights - alloc'
       end if

       !
       ! Now we go back the beginning of the input file, search for the line whether the 
       ! obs. energies start and read all of them (en_npts). 
       ! After that we do the same with the potential energy points. 
       !
       rewind(f_inp)
       !
       char_tmp(1:3)=""
       i = 0
       do while( char_tmp(1:3)/="---".or.i/=5)
         !
         read (f_inp,"(a80)") char_tmp
         if (char_tmp(1:3)=="---") i = i+1
         !
       enddo 

       !
       ! Read the obs. data
       !
       do i=1,en_npts
          !
          read (f_inp,*) J_obs(i),gamma_obs(i),N_obs(i), &
                         ener_obs(i),quanta_obs(1:3,i),wtall(i)
          !
          if (all(jrot_list(:)/=J_obs(i))) wtall(i) = 0
          !
          select case(trim(gamma_obs(i))) 
            !
          case default 
            !
            write(f_out,"('Illegal symmetry ',2a,' in the obs. energy section')") gamma_obs(i)
            stop 'Illegal symmetry '
            !
          case("A1")
            sym_obs(i)  = 1
          case("A2")
            sym_obs(i)  = 2
          case("B1")
            sym_obs(i)  = 3
          case("B2")
            sym_obs(i)  = 4
            !
          end select
          !
       enddo 
       !
       ! Define the molec. symmetry group correlatation with the internal DVR3D symmetry indeces.
       !
       ! symmetrical rot. basis functions 
       !
       nsym(0,0) = 1
       nsym(0,1) = 4
       nsym(1,0) = 3
       nsym(1,1) = 2
       !
       ! asymmetrical rot. basis functions 
       !
       nasym(0,0) = 3
       nasym(0,1) = 2
       nasym(1,0) = 1
       nasym(1,1) = 4
       !
       ! Skip a line
       !
       read (f_inp,"(a80)") char_tmp
       !
       ! Read the potential energy data
       !
       do i=1,pot_npts
         !
         read (f_inp,*) local(1:3,i),pot_values(i),wtall(en_npts+i)
         !
       enddo 
       !
       close(f_inp,status="keep")
       !
       ! write the potential parameters into the file, where  
       ! we store the current values 
       !
       open(f_pot,file='pot.fit',status='replace')
       !
       write(f_pot,*) parmax
       !
       do i=1,parmax
          write(f_pot,"(i8,e24.12,4x,a8)") ivar(i),potparam(i),nampar(i)
       end do
       !
       close(f_pot)
       !
       ! ztran has to be .false. for vibrational (J=0) problem,
       ! and .true. for rovibrational calcualtions. 
       !
       !if (itmax/=0.and.Jmax==0.and.ztran) then 
       !  write (f_inp,"(' for nonzero itmax and jmax ztran must be false')") 
       !  stop 'please disactivate ztran for fitting'
       !endif 
       !
       !if (Jmax/=0.and..not.ztran) then 
       !  write (f_inp,"(' for jmax/=0 ztran must be true')") 
       !  stop 'please activate ztran for J/=0'
       !endif 
       !
       ! Adjust the obs. and pot. weight factor to be
       ! fit_factor =w_ener/w_pot
       ! 
       ! 1. Normalizing the pot. weight factor sum(weights(i)) = number of data 
       !
       !
       wtsum = sum(wtall(1:en_npts))
       wtall(1:en_npts) = wtall(1:en_npts)/wtsum
       !
       wtsum = sum(wtall(en_npts+1:npts))
       wtall(en_npts+1:npts) = wtall(en_npts+1:npts)/wtsum
       !
       ! 2. Factorizing the obs. weights by the factor "fit_factor":
       !
       wtall(1:en_npts) = wtall(1:en_npts)*fit_factor
       !
       ! 3. And normilizing the all weight factors. 
       !
       wtsum = sum(wtall(1:npts))
       wtall(1:npts) = wtall(1:npts)/wtsum
       ! 
       ! Count how many data points actually will be fitted. 
       !
       nused=0
       wtsum=0.0d0
       !
       wt_bit = 0 
       !
       do i=1,npts
         if (wtall(i) .gt. 0.0d0) then 
           nused=nused+1
           !
           wt_bit(i) = 1.0d0
           !
         endif
       enddo
       !
       ! sigma = exp. data precision for robust fit 
       !
       if (robust_fit>0) then
         !
         sigma = 1.0d0
         do i=1,npts
           if (wtall(i)>0) sigma(i) = sigma(i)/sqrt(wtall(i))*robust_fit
         enddo
         !
         !wtsum = 1.0d0 ! sqrt(sum(sigma(1:en_npts)**2))
         !sigma(:) = sigma(:)*robust_fit/wtsum
         !
       endif 
       !
       write(f_out,"('Number of data points used in the fit: ',i9)") nused
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       allocate (derj0(parmax),rjacob(npts,parmax),eps(npts),stat=alloc)
       if (alloc/=0) then
           write (f_out,"(' Error ',i,' initializing fitting objects')") alloc
           stop 'derj0 - alloc'
       end if
       !
       !
       ! The last object to allocate - the lapack related work array
       !
       lwork = 50*parmax
       !
       allocate (wspace(lwork),stat=alloc)
       if (alloc/=0) then 
        write(f_out,"('wspace - out of memory')")  
        stop 'wspace - out of memory'
       endif
       !
       ! we will write the obs.-calc. data into the standard output file (unit=6) 
       !
       open(f_res,file=trim(out_name)//'.res',status='replace')
       !
       ! prepare the file to write all computed energies 
       !
       open(f_en,file=trim(out_name)//'.en',status='replace')
       !
       ! The fitting loop is about to start. 
       !
       ! fititer will count the iterations. 
       !
       fititer = 0
       !
       ! The initial parameters to be remembered. 
       !
       parold=potparam
       !
       ! Create the input file for the expectaion values jobs, 
       ! to be used in connaction with the program 'xpext3.x'
       !
       call create_the_xpect_job_file(parmax,ivar)
       !
       ! The outer fitting loop - allows to restart the fitting in case 
       ! we decide to remove some of the varying parameters. 
       ! This option is working together with fit_type ='linur'
       !
       still_run = .true.
       outer_loop: do while (still_run)  
          !
          ! Initial values for the standard error and  stability.
          !
          stadev_old = 1e10
          stability = 1e10
          stadev    = 1e10
          !
          numpar  = 0
          !
          ! Count the actually varying number of parameters:
          !
          numpar  = 0
          do i=1,parmax
            if (ivar(i) .gt. 0) numpar=numpar+1
          enddo 
          !
          rjacob = 0 
          !
          ! Prepaper the first run:
          !
          close(14)
          write(14,"('start')")
          close(14)
          !
          open(1,file='dvrvib.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
          open(1,file='dvrrot.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
          open(1,file='xpect.out',status='replace')  ; write(1,"('Start fitting...')") ; close(1)
          !
          isys = systemqq('rm dvrvib.out dvrrot.out xpect.out')
          !
          ! The loop starts here. 
          !
          do while( fititer.le.itmax .and. stadev.ge.stadev_best.and. stability.ge.stab_best)
            !
            fititer = fititer + 1
            !
            ! If fit_factor is set zero - there would be no fit to the experimental energies. 
            !
            if (fit_factor<1e-12) rjacob(1:en_npts,:) = 0
            !
            write(f_out,"(/'Iteration = ',i8)") fititer
            !
            do j = 1,jlistmax
               !
               jrot = jrot_list(j)
               !
               ! here we store the number of computed term values for each symmetry block
               neval_t = 0
               !
               ! Two parites have to be evaluated.
               !
               do ipar = 0,1
                 !
                 !isym = 2*ipar+1 ! ; if (ipar==1.and.jrot==0) isym = 4
                 !
                 isym  = nsym(mod(jrot+2,2),ipar)
                 iasym = nasym(mod(jrot+2,2),ipar)
                 !
                 ! Create the job-input file to be used with the dvr3d.x program
                 ! in the J=0 calculations for the givem parity ipar. 
                 !
                 call create_the_vib_job_file(jrot,ipar,meval)
                 !
                 isys = systemqq('rm fort.14')
                 !
                 !isys = systemqq(trim(dir)//'dvr3drjz+ser.exe < '//trim(dir)//'dvr0.job') 
                 !
                 if (verbose>=3) write(f_out,"(/'calling DVR3D program for J = ',i4,', ipar =  ',i2'...')") jrot,ipar
                 !
                 isys = systemqq(dvr3d_exe//'< dvrvib.inp >> dvrvib.out')
                 !
                 if (jrot>0) then
                    !
                    ! The rot. energies will be written from the very first line. 
                    !
                    write(14,"('start')") ;  close(14)
                    isys = systemqq('rm fort.14')
                    !
                    isys = systemqq('cp fort.26 fort.4')
                    !
                    call create_the_rot_job_file(ipar,jrot,meval)
                    !
                    if (verbose>=3) write(f_out,"('callling rotlev program ...')") 
                    ! 
                    isys = systemqq(rotlev_exe//'<dvrrot.inp >>dvrrot.out')
                    !
                 endif 
                 !
                 if (associated(dvr3d(isym)%energy)) deallocate(dvr3d(isym)%energy)
                 if (associated(dvr3d(isym)%derj)) deallocate(dvr3d(isym)%derj)
                 !
                 allocate (dvr3d(isym)%energy(meval),dvr3d(isym)%derj(meval,parmax),stat=alloc)
                 if (alloc/=0) then
                   write (f_out,"(' Error ',i,' initializing dvr3d(isym)%energy')") alloc
                   stop 'dvr3d%energy - alloc'
                 end if
                 !
                 ! Read the energies from fort.14 into dvr3d(isym)%energy.
                 !
                 call get_energies_from_fort14(jrot,ipar,0,meval,dvr3d(isym  )%energy(:))
                 !
                 ! Number of calculated energies.
                 !
                 neval_t(isym) = meval
                 !
                 ! the antysemtric rotational functions:
                 !
                 if (jrot/=0) then 
                   !
                   if (associated(dvr3d(iasym)%energy)) deallocate(dvr3d(iasym)%energy)
                   if (associated(dvr3d(iasym)%derj)) deallocate(dvr3d(iasym)%derj)
                   !
                   allocate (dvr3d(iasym)%energy(meval),dvr3d(iasym)%derj(meval,parmax),stat=alloc)
                   if (alloc/=0) then
                     write (f_out,"(' Error ',i,' initializing dvr3d(iasym)%energy')") alloc
                     stop 'dvr3d%energy - alloc'
                   end if
                   call get_energies_from_fort14(jrot,ipar,1,meval,dvr3d(iasym)%energy(:))
                   !
                   ! Number of calculated energies.
                   !
                   neval_t(iasym) = meval
                   !
                 endif
                 !
                 ! Zero point energy:
                 !
                 if (jrot==0.and.ipar==0) then 
                   ezero_ = dvr3d(1)%energy(1)
                   isys = systemqq('cp fort.14 fort.57')
                 endif 
                 !
                 ! The derivatives of the energy wrt parameters will be evalueted 
                 ! only if 1) it is a fitting job, not just energy calculations (j/=0) and
                 !         2) it is a energy fitting job, not just a fit to the ab initio points 
                 !            (fit_factor is not zero). 
                 !
                 if (itmax.ge.1.and.fit_factor>1e-12.and.trim(deriv_type)=='hellman'.and.mod(fititer+2,3)==0) then 
                   !
                   ! J=0 part is processed differently than J/=0. 
                   ! We need xpect3.exe to be called only once. 
                   !
                   if (jrot==0) then 
                     !
                     isys = systemqq('mv fort.26 fort.11')
                     !
                     if (ipar/=0) isys = systemqq('cp fort.28.em fort.28')
                     !
                     if (verbose>=3) write(f_out,"('calling xpect program for J=0...')")
                     !
                     isys = systemqq(xpect_exe//'< xpect.inp >> xpect.out')
                     !
                     if (ipar==0) isys = systemqq('mv fort.28 fort.28.em')
                     !
                     ! The derivatives are stored in fort.12, we read them into dvr3d(isym)%derj.
                     !
                     call read_jacobi_matrix(meval,parmax,dvr3d(isym)%derj(1:meval,1:parmax) )
                     !
                     if (jrot==0.and.ipar==0) then 
                       derj0(1:parmax) = dvr3d(1)%derj(1,1:parmax)
                     endif 
                     !
                   else
                     !
                     !isys = systemqq('cp fort.26 fort.4')
                     !
                     !if (verbose>=3) write(f_out,"('Call rotlev program ...')") 
                     !
                     !isys = systemqq(rotlev_exe//' <dvrrot.inp >>dvrrot.out')
                     !
                     isys = systemqq('cp fort.28.em fort.28')
                     isys = systemqq('mv fort.8 fort.11')
                     !
                     ! 1. xpect3 #1 
                     !
                     if (verbose>=3) write(f_out,"('calling xpect program for |jk>+|j-k>...')") 
                     !
                     isys = systemqq(xpect_exe//'< xpect.inp >> xpect.out')
                     !
                     call read_jacobi_matrix(meval,parmax,dvr3d(isym)%derj(1:meval,1:parmax) )
                     !
                     isys = systemqq('mv fort.9 fort.11')
                     !
                     ! 2. xpect3 #2 
                     !
                     if (verbose>=3) write(f_out,"('calling xpect program for |jk>-|j-k>...')") 
                     !
                     isys = systemqq(xpect_exe//'< xpect.inp >> xpect.out')
                     !
                     ! The derivatives are stored in fort.12, we read them into dvr3d(iasym)%derj.
                     !
                     call read_jacobi_matrix(meval,parmax,dvr3d(iasym)%derj(1:meval,1:parmax) )
                     !
                   endif 
                   !
                   ! Jacobi matrix stores derivatives only for energies that have obs. counterparts 
                   ! in input with weight /= 0, i.e. participating in the fit. 
                   !
                   do nrow = 1,en_npts
                     !
                     iJ = J_obs(nrow) ; jsym = sym_obs(nrow) !; ipar = mod(isym,2)
                     !
                     if (iJ==Jrot.and.(isym==jsym.or.iasym==jsym)) then 
                       !
                       iener  = n_obs(nrow) ;  if (iener>meval.or.isym>4.or.iasym>4) cycle 
                       !
                       ! Use only the derivatives for parameters with ivar/=0, i.e. varying paramaters. 
                       !
                       ncol=0
                       do  i=1,parmax
                         if (ivar(i) .ne. 0) then
                            !
                            ncol=ncol+1
                            !
                            rjacob(nrow,ncol) = dvr3d(jsym)%derj(iener,i)-derj0(i)
                            !
                          endif
                          ! 
                       enddo 
                       !
                     endif 
                     !
                   enddo
                   !
                 endif 
                 !
               enddo
               !
               ! Print out the calc. and obs.-calc., i.e. result of the fit. 
               ! Only fitted energies are printed. 
               !
               write(f_out,"(/1X,100('-'),/'| ## |  N |  J | sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    v1 v2 v3  |',/1X,100('-'))")
               !
               eps = 0
               !
               do nrow = 1,en_npts
                 !
                 iJ = J_obs(nrow) ; isym = sym_obs(nrow) 
                 !
                 if (iJ==Jrot) then 
                   !
                   iener  = n_obs(nrow) ;  if (iener>neval_t(isym).or.isym>4) cycle ! iener  = 1
                   !
                   enercalc(nrow) = dvr3d(isym)%energy(iener)-ezero_
                   !
                   eps(nrow) = ener_obs(nrow)-enercalc(nrow)
                   !
                   write (f_out,"(4i5,' ',3f13.4,2x,e8.2,5x,3(i3))") &
                      nrow,iener,iJ,isym,enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
                   wtall(nrow),quanta_obs(1:3,nrow)
                   !
                 endif 
                 !
               enddo ! --- nrow
               !
               ! Printing all calculated term values. If the obs. counterpats exist, 
               ! the obs.-calc. are printed as well. 
               ! This list can be used to identify the obs-calc pairs, i.e. the number of 
               ! a term value as it appear in the list of calculated energies. 
               !
               write(f_en,"(/'Iteration = ',i4)") fititer
               !
               write(f_en,"(/1X,100('-'),/'| ## |  N |  J | Sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    v1 v2 v3  |',/1X,100('-'))")
               !
               do isym = 1,4
                  !
                  do irow = 1,neval_t(isym)
                     !
                     if (irow==1.or.dvr3d(isym)%energy(irow)-ezero_>sqrt(epsilon(1.0))) then 
                        !
                        nrow = 1
                        !
                        do while (nrow/=en_npts+1.and..not.(irow==n_obs(nrow).and.Jrot==J_obs(nrow).and.isym==sym_obs(nrow)) )
                          !
                          nrow = nrow+1
                          !
                        enddo
                        !
                        if (nrow<en_npts) then
			  !
			  mark = " "
			  !
			  if (abs(eps(nrow))>1.0) mark = "!"
                          !
                          write(f_en,"(4i5,' ',3f13.4,2x,e8.2,2x,a1,2x,3(i3))") irow,n_obs(nrow),Jrot,sym_obs(nrow),enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
                          wtall(nrow),mark,quanta_obs(1:3,nrow)
                          !
                        else
                          !
                          write(f_en,"(4i5,' ',3f13.4,2x,e8.2,5x,3(i3))") irow,0,Jrot,isym,0.0,dvr3d(isym)%energy(irow)-ezero_,0.0,0.0
                          !
                        endif 
                        !
                     endif 
                     !
                  enddo
                  !
               enddo 
               !
            enddo 
            !
            ! Alternative way of calculating the derivatives  - with the finite 
            ! differencies. It is essentially slower and we use it only for 
            ! the testing of the xpect3 derivativies.
            !
            if (trim(deriv_type)/='hellman'.and.itmax.ge.1.and.fit_factor>1e-12) then
              !
              call finite_diff_of_energy(rjacob,potparam)
              !
            endif 
            !
            ! Here the potential energy section starts. 
            !
            do nrow=1,pot_npts
              !
              !  we assume here that the angles are written in degrees, 
              !  not in radians. I.e. we have to convert the degrees to in radians.
              !     
              r1   =local(1,nrow)
              r2   =local(2,nrow)
              theta=local(3,nrow)*pi/1.8d+02
              !
              ! Call the potential function. It has to be included into the "poten" subroutine
              !
              call poten(potparam,v,r1,r2,theta)
              !
              ! eps - epsilon = ab initio energies - calculated pot. energies,
              ! where we comntinue counting the fitting data points starting with en_npts - 
              ! obs. data. 
              !
              eps(nrow+en_npts) = pot_values(nrow)-v
              !
              ! Calculate derivatives with respect to parameters 
              ! using the finite diff. method. 
              !
              if (itmax.ge.1.and.(mod(fititer-1,nofititer+1).eq.0)) then
                !
                ncol=0
                numpar = 0
                do  i=1,parmax
                  if (ivar(i) .ne. 0) then
                     !
                     ncol=ncol+1
                     tempx=potparam(i)
                     deltax=fitfactordeltax*abs(tempx)
                     if (deltax .le. 1e-15) deltax=1e-5
                     !
                     potparam(i)=tempx+deltax
                     call poten(potparam,potright,r1,r2,theta)
                     !
                     potparam(i)=tempx-deltax
                     !
                     call poten(potparam,potleft,r1,r2,theta)
                     !
                     potparam(i)=tempx
                     rjacob(nrow+en_npts,ncol)=(potright-potleft)/(2.0d0*deltax)
                     !
                  endif
                enddo ! --- ncol
                !
                numpar = ncol
                !
              endif     
              !
            enddo  ! ---  nrow
            !
            ! ssq  - weighted rms**2, rms  - root mean square deviation. 
            !
            ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
            rms=sqrt(sum(eps(1:npts)*eps(1:npts))/npts)
            !
            ! Prepare the linear system a x = b as in the Newton fitting approach.  
            !
            if (itmax.ge.1) then
               !----- form the a and b matrix ------c
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow    
                   al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))
                   al(icolumn,irow)=al(irow,icolumn)
                 enddo
               enddo
               !
               ! form B matrix 
               do irow=1,numpar      
                 bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
               enddo   
               !
               ! Two types of the linear solver are availible: 
               ! 1. linur (integrated into the program, from Ulenikov Oleg)
               ! 2. dgelss - Lapack routine (recommended).
               !
               select case (trim(fit_type)) 
               !
               case default
                 !
                 write (f_out,"('fit_type ',a,' unknown')") trim(fit_type)
                 stop 'fit_type unknown'
                 !
               case('linur') 
                 !
                 call linur(numpar,numpar,al(1:numpar,1:numpar),bl(1:numpar),dx(1:numpar),ierror)
                 !
                 ! In case of dependent parameters  "linur" reports an error = ierror, 
                 ! which is a number of the dependent parameter. We can remove this paramter 
                 ! from the fit and set its value to zero. And start the iteration again. 
                 !
                 if (ierror.ne.0) then 
                   ncol=0
                   do i=1,parmax
                      if (ivar(i) .ne. 0) then
                        ncol=ncol+1
                        if  ( ncol.eq.ierror ) then 
                            ivar(i) = 0
                            potparam(i) = parold(i)
                            write(f_out,"(i,'-th is out - ',a8)") i,nampar(i)
                        endif 
                      endif 
                   enddo 
                   cycle outer_loop    
                  endif 
                  !
               case ('dgelss')
                 !
                 ai = al 
                 call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bl(1:numpar),numpar,Tsing,1.D-12,RANK0,wspace,lwork,info)
                 !
                 if (info/=0) then
                   write(f_out,"('dgelss:error',i)") info
                   stop 'dgelss'
                 endif
                 !
                 dx = bl
                 !
               end select 
               !
               !----- update the parameter values ------!
               !
               ncol=0
               do i=1,parmax
                if (ivar(i) > 0) then
                     ncol=ncol+1
                     potparam(i)=potparam(i)+dx(ncol)
                  endif
               enddo
               !
               if (robust_fit>0) then
                 !
                 call robust_fitting(numpar,sigma(1:npts),eps(1:npts),wtall(1:npts))
                 !
                 ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
                 !
               endif 
               !
               ! write the potential parameters into the file, where  
               ! we store the current values 
               !
               open(f_pot,file='pot.fit',status='replace')
               !
               write(f_pot,*) parmax
               !
               do i=1,parmax
                  write(f_pot,"(i8,e24.12,4x,a8)") ivar(i),potparam(i),nampar(i)
               end do
               !
               close(f_pot)
               !
               ! Estimate standard deviation error. 
               !
               if ( nused.ne.numpar ) then 
                 stadev=dsqrt(ssq/float(nused-numpar))
               else 
                 stadev=dsqrt(ssq/nused)
               endif
               !
               ! Estimate the standard errors for each parameter using 
               ! the inverse matrix of a. 
               !
               call invmat(al,ai,numpar,parmax)
               !
               sum_sterr=0.d0
               ncol = 0 
               do i=1,parmax
                  if (ivar(i) > 0) then
                      ncol=ncol+1
                     if (nused.eq.numpar) then  
                        sterr(ncol)=0
                     else
                        sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                        sum_sterr=sum_sterr+abs(sterr(ncol)/potparam(i))
                     endif
                   endif
               enddo    
               !
               sum_sterr=sum_sterr/numpar 
               !
               ! This is how we define stability of the fit:
               ! as a relative change of stadev comparing with the step before. 
               !
               stability=abs( (stadev-stadev_old)/stadev )
               stadev_old=stadev
               !
            else
               !
               stadev=dsqrt(ssq/nused)
               !
            endif
            !
            ! Print the updated parameters. 
            !
            write(f_out,"(/'Potential paramters:')")
            !
            do i=1,parmax
              write (f_out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),potparam(i)
            enddo
            !
            ! Print the potential energy points into a separate unit. 
            !
            inquire(f_potpoints,opened=ifopen)
            if ( ifopen ) then
               rewind(f_potpoints)
            else
               open (f_potpoints,file=trim(out_name)//'.pot',status='replace' )
            endif
            !
            write(f_potpoints,"(1h1,5x,12('*'),' ab initio points ',  &
                 12('*')// &
                 4x,'r1',5x,'r2',5x,'theta',7x,'ab initio PES',3x, &
                 'cal.PES',3x,'a-c',3x,'weight'/)")
            !
            do nrow=1,pot_npts
              !
              r1   =local(1,nrow)
              r2   =local(2,nrow)
              theta=local(3,nrow)
              !
               v = pot_values(nrow)-eps(nrow+en_npts)
               write (f_potpoints,"(2f8.4,' ',f10.2,' ',3f12.2,e10.2)") r1,r2,theta, &
                      pot_values(nrow),v, &
                      eps(nrow+en_npts),wtall(nrow+en_npts) 
              !
            enddo
            !
            ! Output some staqtistics and results 
            !
            !  only if we are fitting:  
            !
            if (itmax.ne.0) then
              !
              !still_run = .false.
              ! 
              write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
              l = 0 
              do i=1,parmax
                if (ivar(i) .ne. 0) then
                   l=l+1
                   ndigits = 0
                   conf_int = sterr(l)
                   do while (conf_int.le.10.0)
                     ndigits = ndigits +1 
                     conf_int = conf_int*10
                   enddo
                   !
                   !if (ndigits<1) ndigits = 12
                   !
                   if (conf_int>1e8) conf_int = 0 
                   !
                   write(fmt0,"(i2)") ndigits
                   fmt0 = adjustl(fmt0)
                   fmt1 = '(a8,i4,2x,f22.'//fmt0//'    ''    ('',i14,'')'')'
                   write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i),nint(conf_int)
                else 
                   ndigits =2
                   if (potparam(i).ne.0.0) ndigits = 8

                   write(fmt0,"(i2)") ndigits
                   fmt0 = adjustl(fmt0)
                   fmt1 = '(a8,I4,2x,f22.'//fmt0//')'
                   write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i)
                endif
              enddo  ! --- i
              !
            endif 
            !
            still_run = .false.
            !
            rewind(f_res)
            !
            write(f_res,"(/1X,100('-'),/'| ## |  N |  J | tau|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    v1 v2 v3  |',/1X,100('-'))")
            !
            do nrow = 1,en_npts
               !
               iJ = J_obs(nrow) ; isym = sym_obs(nrow) 
               !
               iener  = n_obs(nrow) ;  if (iener>neval_t(isym).or.isym>4) cycle ! iener  = 1
               !
               if (iJ<=maxval(jrot_list(:))) then 
                  !
                  write (f_res,"(4i5,' ',3f13.4,2x,e8.2,5x,3(i3))") &
                      nrow,iener,iJ,isym,enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
                      wtall(nrow),quanta_obs(1:3,nrow)
                  !
               endif 
               !
            enddo
            !
            ! Print out the ssq for the rovib. energies and pot. data points separetely:
            !
            ssq1 = 0 ; ssq2 = 0 
            !
            wtsum = sum(wt_bit(1:en_npts))
            !
            if (wtsum/=0) ssq1 = sqrt( sum(eps(1:en_npts)**2*dble(wt_bit(1:en_npts)))/wtsum )
            !
            wtsum = sum(wt_bit(1+en_npts:npts))
            !
            if (wtsum/=0) ssq2 = sqrt( sum(eps(1+en_npts:npts)**2*dble(wt_bit(1+en_npts:npts)))/wtsum )

            rms1=sqrt(sum(eps(1:en_npts)**2)/en_npts)
            rms2=sqrt(sum(eps(1+en_npts:npts)**2)/pot_npts)

            !
            write (f_out,6552) fititer,nused,numpar,stadev,ssq1,ssq2,stability
            write (f_res,6553) fititer,nused,numpar,stadev,rms1,rms2,stability
            !
            !
          enddo  ! --- fititer
          !
       enddo outer_loop
       !
       ! Print out the final information: paramteres and eneries.
       !
       rewind(f_res)
       !
       write (f_res,"(/36('-')/' Parameter | W |      Value        |',/36('-')/)")
       !
       do i=1,parmax
          write (f_res,"(A8,1X,I4,2X,e18.8)") nampar(i),ivar(i),potparam(i)
       enddo
       !
       write(f_res,"(/1X,100('-'),/'| ## |  N |  J | tau|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    v1 v2 v3  |',/1X,100('-'))")
       !
       do nrow = 1,en_npts
          !
          iJ = J_obs(nrow) ; isym = sym_obs(nrow) 
          !
          !iener  = n_obs(nrow)
          !
          iener  = n_obs(nrow) ;  if (iener>neval_t(isym).or.isym>4) cycle ! iener  = 1
          !
          if (iJ<=maxval(jrot_list(:))) then 
             !
             write (f_res,"(4i5,' ',3f13.4,2x,e8.2,5x,3(i3))") &
                 nrow,iener,iJ,isym,enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
                 wtall(nrow),quanta_obs(1:3,nrow)
             !
          endif
          !
       enddo
       !
       !
       do nrow=1,pot_npts
         !
         r1   =local(1,nrow)
         r2   =local(2,nrow)
         theta=local(3,nrow)
         !
         v = pot_values(nrow)-eps(nrow+en_npts)
         write (f_res,"(2f8.4,' ',f10.2,' ',3f12.2,e10.2)") r1,r2,theta, &
                 pot_values(nrow),v, &
                 eps(nrow+en_npts),wtall(nrow+en_npts) 
         !
       enddo
       !
       ssq1 = 0 ; ssq2 = 0 
       !
       wtsum = sum(wt_bit(1:en_npts))
       !
       if (wtsum/=0) ssq1 = sqrt( sum(eps(1:en_npts)**2*dble(wt_bit(1:en_npts)))/wtsum )
       !
       wtsum = sum(wt_bit(1+en_npts:npts))
       !
       if (wtsum/=0) ssq2 = sqrt( sum(eps(1+en_npts:npts)**2*dble(wt_bit(1+en_npts:npts)))/wtsum )
       !
       write (f_res,6551) fititer,nused,numpar,stadev,ssq1,ssq2
       write (f_res,6550) fititer,nused,numpar,stadev,rms,stability

  

6550   format(/3X,67('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms      | Stability |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E10.5,' | ',E10.5,'  |  ',&
            E8.3,' |',/3X,67('-')/)


6551   format(/3X,67('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq1     |   ssq2    |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E10.5,' | ',E10.5,'  |  ',&
            E8.3,' |',/3X,67('-')/)

6552   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq_ener |   ssq_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E10.5,' | ',E10.5,'  |  ',&
            E8.3,' |',E10.3,' |',/3X,80('-')/)

6553   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms_ener |   rms_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E10.5,' | ',E10.5,'  |  ',&
            E8.3,' |',E10.3,' |',/3X,80('-')/)

       !
       call SYSTEM_CLOCK(itime2,irate2,imax2)
       !
       itime=(itime2-itime0)/irate2
       write(f_out,"(/i10,' secs CPU time used'/)") itime
       !

  contains 


    subroutine finite_diff_of_energy(rjacob,potparam)
      !
      real(8) :: rjacob(:,:),potparam(:)
      double precision,allocatable :: enerleft(:),enerright(:)
      !
      ! calculate derivatives with respect to parameters
       allocate (enerright(en_npts),enerleft(en_npts),stat=alloc)
       if (alloc/=0) then
           write (f_out,"(' Error ',i,' initializing enerright objects')") alloc
           stop 'enerright - alloc'
       end if
      !
      !
      ncol=0
      !
      do  i=1,parmax
        if (ivar(i) > 0) then
          !
          write( f_res,"('jacob:',i4,'-th param')") i
          !
          ncol=ncol+1
          tempx=potparam(i)
          deltax=fitfactordeltax*abs(tempx)
          if (deltax .le. 1e-15) deltax=1e-6
          !
          select case ( findif )
          case default
               write (f_out,"(' Bad finite diff. parameter',i)") findif
               stop 'Bad finite diff. parameter'
          case (2)
             !
             enerright = enercalc
             enerleft  = enercalc
             !
             if (mod(fititer,2).eq.0) then 
               !
               potparam(i)=tempx+deltax
               !
               do j = 1,jlistmax
                  !
                  jrot = jrot_list(j)
                  !
                  do nrow = 1,en_npts
                     !
                     iJ = J_obs(nrow) ; isym = sym_obs(nrow) 
                     !
                     if (iJ==Jrot) then 
                        !
                        iener  = n_obs(nrow) ; if (iener>neval_t(isym).or.isym>4) iener  = 1
                        !
                        enerright(nrow) = dvr3d(isym)%energy(iener)-ezero_
                        !
                     endif
                     !
                  enddo 
                  ! 
               enddo 
               !
             endif   ! --- flag_fit
             !
             if (mod(fititer,2).eq.1) then 
               !
               potparam(i)=tempx-deltax
               !
               do j = 1,jlistmax
                  !
                  jrot = jrot_list(j)
                  !
                  do nrow = 1,en_npts
                     !
                     iJ = J_obs(nrow) ; isym = sym_obs(nrow) 
                     !
                     if (iJ==Jrot) then 
                        !
                        iener  = n_obs(nrow) ; if (iener>neval_t(isym).or.isym>4) iener  = 1
                        !
                        enerleft(nrow) = dvr3d(isym)%energy(iener)-ezero_
                        !
                     endif
                     !
                  enddo ! --- nrow
                  ! 
               enddo 
               !
             endif   ! --- flag_fit
             !
             rjacob(1:en_npts,ncol)=(enerright(1:en_npts)-enerleft(1:en_npts))/(deltax)
             !
          case (3)
             !
             potparam(i)=tempx+deltax
             !
             ! write the potential parameters into the file, where  
             ! we store the current values 
             !
             open(f_pot,file='pot.fit',status='replace')
             !
             write(f_pot,*) parmax
             !
             do l=1,parmax
                write(f_pot,"(i8,e24.12,4x,a8)") ivar(l),potparam(l),nampar(l)
             end do
             !
             close(f_pot)
             !
             do ipar = 0,1
               !
               isym = 2*ipar+1 ; if (ipar==1.and.jrot==0) isym = 4
               !
               ! Create the job-input file to be used with the dvr3d.x program
               ! in the J=0 calculations for the givem parity ipar. 
               !
               call create_the_vib_job_file(jrot,ipar,meval)
               !
               isys = systemqq('rm fort.14')
               !
               !isys = systemqq(trim(dir)//'dvr3drjz+ser.exe < '//trim(dir)//'dvr0.job') 
               !
               isys = systemqq(dvr3d_exe//'< dvrvib.inp >> dvrvib.out')
               !
               if (jrot>0) then
                  !
                  isys = systemqq('cp fort.26 fort.4')
                  !
                  call create_the_rot_job_file(ipar,jrot,meval)
                  ! 
                  isys = systemqq(rotlev_exe//'<dvrrot.inp >>dvrrot.out')
                  !
               endif 
               !
               if (associated(dvr3d(isym)%energy)) deallocate(dvr3d(isym)%energy)
               !
               allocate (dvr3d(isym)%energy(meval),stat=alloc)
               if (alloc/=0) then
                 write (f_out,"(' Error ',i,' initializing dvr3d(isym)%energy')") alloc
                 stop 'dvr3d%energy - alloc'
               end if
               !
               ! Read rhe energies from fort.14 into dvr3d(isym)%energy.
               !
               call get_energies_from_fort14(jrot,ipar,0,meval,dvr3d(isym  )%energy(:))
               !
             enddo
             !
             do j = 1,jlistmax
                !
                jrot = jrot_list(j)
                ! 
                do nrow = 1,en_npts
                   !
                   iJ = J_obs(nrow) ; isym = sym_obs(nrow)
                   !
                   if (iJ==Jrot) then 
                      !
                      iener  = n_obs(nrow) ; if (iener>neval_t(isym).or.isym>4) cycle
                      !
                      enerright(nrow) = dvr3d(isym)%energy(iener)
                      !
                   endif
                   !
                enddo ! --- nrow
                ! 
             enddo 
             !
             potparam(i)=tempx-deltax
             !
             ! write the potential parameters into the file, where  
             ! we store the current values 
             !
             open(f_pot,file='pot.fit',status='replace')
             !
             write(f_pot,*) parmax
             !
             do l=1,parmax
                write(f_pot,"(i8,e24.12,4x,a8)") ivar(l),potparam(l),nampar(l)
             end do
             !
             close(f_pot)
             !
             !
             do ipar = 0,1
               !
               isym = 2*ipar+1 ; if (ipar==1.and.jrot==0) isym = 4
               !
               ! Create the job-input file to be used with the dvr3d.x program
               ! in the J=0 calculations for the givem parity ipar. 
               !
               call create_the_vib_job_file(jrot,ipar,meval)
               !
               isys = systemqq('rm fort.14')
               !
               !isys = systemqq(trim(dir)//'dvr3drjz+ser.exe < '//trim(dir)//'dvr0.job') 
               !
               isys = systemqq(dvr3d_exe//'< dvrvib.inp >> dvrvib.out')
               !
               if (jrot>0) then
                  !
                  isys = systemqq('cp fort.26 fort.4')
                  !
                  call create_the_rot_job_file(ipar,jrot,meval)
                  ! 
                  isys = systemqq(rotlev_exe//'<dvrrot.inp >>dvrrot.out')
                  !
               endif 
               !
               if (associated(dvr3d(isym)%energy)) deallocate(dvr3d(isym)%energy)
               !
               allocate (dvr3d(isym)%energy(meval),stat=alloc)
               if (alloc/=0) then
                 write (f_out,"(' Error ',i,' initializing dvr3d(isym)%energy')") alloc
                 stop 'dvr3d%energy - alloc'
               end if
               !
               ! Read rhe energies from fort.14 into dvr3d(isym)%energy.
               !
               call get_energies_from_fort14(jrot,ipar,0,meval,dvr3d(isym  )%energy(:))
               !
             enddo
             !
             do j = 1,jlistmax
                !
                jrot = jrot_list(j)
                !
                do nrow = 1,en_npts
                   !
                   iJ = J_obs(nrow) ; isym = sym_obs(nrow)
                   !
                   if (iJ==Jrot) then 
                      !
                      iener  = n_obs(nrow) ; if (iener>neval_t(isym).or.isym>4) cycle
                      !
                      enerleft(nrow) = dvr3d(isym)%energy(iener)
                      !
                   endif
                   !
                enddo ! --- nrow
                ! 
             enddo 
             !
             rjacob(1:en_npts,ncol)=(enerright(1:en_npts)-enerleft(1:en_npts))/(2.0d0*deltax)
             !
             rjacob(2:en_npts,ncol) = rjacob(2:en_npts,ncol) - rjacob(1,ncol)
             !
          end select
          !
          potparam(i)=tempx
          !
          ! write the potential parameters into the file, where  
          ! we store the current values 
          !
          open(f_pot,file='pot.fit',status='replace')
          !
          write(f_pot,*) parmax
          !
          do l=1,parmax
             write(f_pot,"(i8,e24.12,4x,a8)") ivar(l),potparam(l),nampar(l)
          end do
          !
          close(f_pot)
          !
        endif
      enddo ! --- ncol
      !
      deallocate (enerright,enerleft)
      !
    end subroutine finite_diff_of_energy


  end subroutine fitting
  !



  !
  subroutine linur(dimen,npar,coeff,constant,solution,error)

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur


!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  integer,intent(in)           :: npar,dimen
  double precision,intent(in)  :: al(npar,npar)
  double precision,intent(out) :: ai(npar,npar)
  double precision             :: h(npar),p,q
  integer                      :: i1,i2,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine invmat
!------------------------------------------!

    subroutine create_the_vib_job_file(jrot,ipar,neval)
      !
      integer,intent(in) :: jrot,ipar
      integer,intent(out):: neval
      integer            :: icard
      integer            :: npnt2_t,jrot_t,neval_t,nalf_t,max2d_t,max3d_t,idia_t,&
                                           kmin_t,npnt1_t,ipar_t,max3d2_t
      character(len=200)  ::  char_job_(10)
      !
      open(1,file='dvrvib.inp')
      !
      char_job_(1:10) = char_vib_job(1:10)
      !
      write(1,"(a200)") char_job(1)
      write(1,"(a20)") char_job( 2)
      !
      icard = icard_j(jrot)
      !
      if (jrot==0.and.ipar==0) then
        !
        !read(char_job(3),"(11i5)") npnt2_,jrot_,neval_,nalf_,max2d_,max3d_,idia_,&
        !                    kmin_,npnt1_,ipar_,max3d2_
        !
        write(1,"(a50)") char_vib_job(icard)
        !
        neval = nevalx(icard)
        !
        if (verbose>=4) write(f_out,"(/'neval (J=0,A1) =  ',i)") neval
         
      else
        !
        ! Change the size of the basis set according to the value of jrot
        !
        !neval = max(jrot,1)*neval_
        !
        neval = nevalx(icard)
        !
        !npnt2x(jrot)
        !

        read(char_vib_job(icard),"(11i5)") npnt2_t,jrot_t,neval_t,nalf_t,max2d_t,max3d_t,idia_t,&
                                           kmin_t,npnt1_t,ipar_t,max3d2_t
        !
        if (nevalx(icard)/=neval_t.or.npnt2x(icard)/=npnt2_t.or.jrot/=jrot_t) then
           !
           write (f_out,"('create_the_vib_job_file: inconsistence with input parameters for icard = :',i)") icard
           write (f_out,"('nevalx(icard)/=neval_tor.npnt2x(icard)/=npnt2_t.or.jrot/=jrot_',6i)") &
           nevalx(icard),neval_t,npnt2x(icard),npnt2_t,jrot,jrot_t
           stop 'create_the_vib_job_file: inconsistence in parameters'
           !
        endif 
        !
        write(1,"(11i5)") npnt2_t,jrot,neval,nalf_t,max2d_t,max3d_t,idia_t,&
                          kmin_t,npnt1_t,ipar,max3d2_t
        !
      endif
      !
      write(1,"(a41)") char_job( 3)
      write(1,"(a74)") char_job( 4)
      write(1,"(a56)") char_job( 5)
      write(1,"(a56)") char_job( 6)
      write(1,"(a53)") char_job( 7)
      write(1,"(a53)") char_job( 8)
      write(1,"(a53)") char_job( 9)
      !
      close(1,status='keep')
      !
    end subroutine create_the_vib_job_file



    subroutine get_energies_from_fort14(jrot,ipar,irottau,meval,enercalc)
      !
      integer,intent(in)   :: jrot,ipar,irottau
      integer,intent(out)  :: meval
      integer              :: i_t(6)
      !
      real*8,intent(out)   :: enercalc(:)
      !
      open(2,file='fort.14',action='read',status='old')
      !
      if (ipar==0.and.jrot==0) read(2,*) 
      !
      read(2,*) i_t(1:6)
      !
      ! Number of found solutions
      !
      !
      meval = i_t(6)
      !
      if (verbose>=4) write(f_out,"(/'meval (found energies) =  ',i7,/'jrot,ipar,irottau = ',3i3)") meval,jrot,ipar,irottau
      !
      if (size(enercalc)<meval) then 
        !
        write(f_out,"('get_energies_from_fort14: size of enercalc(:) is too small, meval = ',i8,' maxener =  ',i8)") meval,size(enercalc)
        stop 'enercalc is too small'
        !
      endif 
      !
      ! Read the energies from fort.14 
      !
      read(2,*) enercalc(1:meval)
      !
      ! read the energies for anti-symmetric rotational functions. 
      !
      if (irottau==1) then
        !
        read(2,*) i_t(1:6)
        !
        ! Number of found solutions
        !
        meval = i_t(6)
        !
        if (size(enercalc)<meval) then 
          !
          write(f_out,"('get_energies_from_fort14: size enercalc too small (asym), meval= ',i8,' maxener= ',i8)") & 
                         meval,size(enercalc)
          stop 'enercalc asym is too small'
          !
        endif 
        !
        ! Read the energies from fort.14 
        !
        read(2,*) enercalc(1:meval)
        !
      endif 
      !
      enercalc(1:meval) = enercalc(1:meval)*219474.624
      !
      close(2,status='keep')
      !
    end subroutine get_energies_from_fort14
    !
    !
    subroutine create_the_rot_job_file(ipar,jrot,neval)
	  !
      integer,intent(in)   :: ipar,jrot
      integer,intent(inout) :: neval
      integer              :: nvib,kmin,i_t(6),icard,nvib_t,neval_t,kmin_t,ibass_t,neval2_t
      !
      !open(1,file='fort.14',action='read',status='old')
      !
      !if (ipar==0) read(2,*) 
      !
      !read(2,*) i_t(1:6)
      !
      !close (1)
      !
      nvib = neval ! i_t(6)
      !
      open(1,file='dvrrot.inp')
      !
      icard = icard_j(jrot)
      !
      write(1,"(a200)") char_rot_job(icard,1)
      !
      read(char_rot_job(icard,2),"(11i5)") nvib_t,neval_t,kmin_t,ibass_t,neval2_t
      !
      write(1,"(11i5)") nvib,neval_t,kmin_t,ibass_t,neval2_t
      !
      ! write(1,"(a200)") char_rot_job(icard,2)
      !
      !nvib = nevalx(jrot)
      !
      neval = neval_t
      !
      !write(1,"(11i5)") nvib,neval,kmin_,ibass_,neval2_
      write(1,*) 
      write(1,"(f10.4)") ezero_
      !
      close(1,status='keep')
      !
    end subroutine create_the_rot_job_file
    !
    subroutine create_the_xpect_job_file(nparams,ivar)
	  !
      integer,intent(in)  :: nparams,ivar(nparams)
      integer             :: lpot,npropin,nprt,nv1,i
      !
      open(1,file='xpect.inp')
      !
      write(1,"('&PRT zform=.true., /')")
      write(1,"('Expectation values, fort.11,properties=3')")
      !nvib = neval_
      !
      lpot = nalf_ ! 38
      npropin = nparams
      nprt = 0
      nv1 = 0
      !
      write(1,"(4i5)") lpot,npropin,nprt,nv1
      do i = 1,npropin 
        !
        if (ivar(i)>0) then
          !
          write(1,"(i5)") 1
          !
        else
          !
          write(1,"(i5)") 0
          !    
        endif
        !
      enddo
      !
      close(1,status='keep')
      !
    end subroutine create_the_xpect_job_file
    !
    !
    subroutine read_jacobi_matrix(enermax,parmax,derj)
	  !
      integer,intent(in)  :: enermax,parmax
      real(8),intent(out) :: derj(enermax,parmax)
      integer             :: ie,ipar
      !
      open(1, file='fort.12',status='old')
      !
      !derj(1,:) = 0 
      !
      do ie = 1,enermax 
        !
        read(1,*)
        !
        do ipar = 1,parmax 
          !
          read(1,*) derj(ie,ipar)
          !    
        enddo
        !
      enddo
      !
      close(1,status='keep')
      !
    end subroutine read_jacobi_matrix

    subroutine robust_fitting(numpar,sigma,eps,wt)

      integer,intent(in) :: numpar
      !real(8),intent(inout) :: a_wats
      real(8),intent(in)    :: sigma(:),eps(:)
      real(8),intent(inout) :: wt(:)
      !
      integer            :: npts,i,nrow,nused
      real               :: da1,da2,wtsum,da
      !
      npts = size(sigma)
      !
      nused = 0 
      do i=1,npts
        if (wt(i)>0) nused=nused+1
      enddo
      !
      !Watson alpha-parameter
      ! 
      !do i = 1,-1
      !  !
      !  da1 = 0
      !  da2 = 0
      !  !
      !  do nrow=1,npts
      !    if (wt(nrow)>small_) then 
      !      da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
      !      da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
      !    endif 
      !  enddo 
      !  !
      !  da =( da1 -  real(nused-numpar,rk) )/da2
      !  a_wats = a_wats + da
      !  !
      !  if (a_wats<sqrt(small_)) a_wats = 1e-3+real(i,rk)*1e-2
      !  !
      !enddo
      ! 
      !  adjusting the weights  
      ! 
      do nrow=1,npts
         if (wt(nrow)>0) wt(nrow) = 1.0d0/( sigma(nrow)**2 + eps(nrow)**2 )
      enddo 
      ! 
      ! "re-normalizing" the weight factors
      !
      wtsum = sum(wt(1:npts))
      wt(1:npts) = wt(1:npts)/wtsum
      !
    end subroutine robust_fitting
 

  end module fit_module
  !
  program driver
    use fit_module

    call fitting()

  end program driver


