!>hold main variable globally, all other modules will need this module included
!>to access the filament/particle arrays and other important variables
!>also contains important routines for initialising random number generator and
!>reading in runfile
module cdata
  use omp_lib
  !**********VORTEX FILAMENT******************************************************
  !>our main structure which holds vortex points
  !!@param x position of the vortex point  
  !!@param u velocity of the vortex point
  !!@param u1 @param u2 stored velocities for Adams-Bashforth
  !!@param ghosti @param ghostb ghost particles for periodic b.c's
  !!@param ghostii @param ghostbb ghost particles for periodic b.c's  
  !!@param infront @param behind flag to make points an orientated filament
  !!@param closest closest particle, used in reconnections
  !!@param closestd separation between closest particle and particle 
  !!@param pinnedi @param pinnedb if the particle is pinned to boundary
  !!@param wpinned which boundary are we pinned to ?
  type qvort 
    real :: x(3)
    real :: u(3), u1(3), u2(3)
    real :: ghosti(3), ghostb(3)
    integer :: infront, behind 
    integer :: closest
    real :: closestd
    integer :: wpinned(3)
    logical :: pinnedi=.false., pinnedb=.false.
  end type
  !>main filament vector
  type(qvort), allocatable :: f(:) 
  !>number of vortex points in the simulation - size of f is pcount
  integer :: pcount
  !**************TIME PARAMS*******************************************************
  !>time held globally
  real :: t=0. 
  !>current timestep
  integer :: itime 
  !>integer loop starts from (altered by reading in stored data - i.e. restarting)
  integer :: nstart=1 
  !***********DIAGNOSTIC INFO******************************************************
  !>total number of reconnections
  integer :: recon_count=0 
  !>total number of reconnections with boundaries
  integer :: wall_recon_count=0 
  !>total number of particle removals due to contraction of filament
  integer :: remove_count=0 
  !>total length of filaments
  real :: total_length
  !>average separation of the vortex points 
  real :: avg_sep
  !>maximum velocity
  real :: maxu
  !>maximum velocity change - acceleration x dt 
  real :: maxdu
  !>mean curvature
  real :: kappa_bar 
  real :: kappa_min, kappa_max !min/max curvature 
  !>number of evalulations when using tree algorithm
  integer :: tree_eval=0
  !***********CONSTANTS************************************************************
  !some constants - precompute for speed
  real, parameter :: pi=3.14159265358979324
  real, parameter :: rootpi=1.77245385090551
  real, parameter :: one_half = (1./2.)
  real, parameter :: three_twos=(3./2.)
  real, parameter :: twenty_three_twelve=(23./12.)
  real, parameter :: four_thirds=(4./3.)
  real, parameter :: five_twelths=(5./12.)
  !***********RUN.IN***************************************************************
  !parameters from run.in, given protected status so treated like parameters
  !by routines in the rest of the code...
  !--------main parameters-please set these in run.in-----------------------------
  integer, protected :: nsteps, shots, recon_shots=1
  integer, protected :: init_pcount
  real, protected ::  delta
  !>timestep
  real :: dt
  real, protected :: box_size=0., cylind_r=0.
  real, protected :: quant_circ=9.97E-4 !He-4 by default
  real , protected :: corea=8.244023E-9 !He-4 by default
  !>are boundaries periodic, open or solid?  
  character(len=40) :: boundary_x='periodic' !set to periodic by default
  character(len=40) :: boundary_y='periodic'
  character(len=40) :: boundary_z='periodic'
  integer :: n_periodic=0 !number of periodic dimensions
  integer, allocatable :: periodic_loop_array(:,:)
  integer :: n_mirror=0 !number of solid dimensions
  integer, allocatable :: mirror_loop_array(:,:)
  logical :: cylindrical_boundaries=.false. !used if we have cylindrical boundaries
  !key arguements that must be set
  character(len=30), protected :: velocity, initf, boundary
  !order of derivatives
  !-----------arguements used by initial.mod/initial_cond.mod-------------
  integer, protected :: line_count=1
  real, protected :: scale_factor=1 ! in terms of box size
  real, protected :: rotation_factor=1 !1=2*pi, 0=0
  !how much we translate random_loops initial condition by 
  real, protected :: loop_translate(3)=1.!for separate xyz components
  !--------the following parameters add special features-------------------------
  !---------------------tree-algorithm-------------------------------------------
  real, protected :: tree_theta=0.
  logical, protected :: tree_print=.false.
  logical, protected :: tree_extra_correction=.true.
  !------------normal fluid component--------------------------------------------
  character(len=30), protected :: normal_velocity='zero'
  real, protected :: alpha(2)=0. !mutual friction coefficients
  real, protected :: norm_vel_xflow=0.5
  !----------------------------code testing---------------------------------------
  logical, protected :: NAN_test=.true.!test for NANs in arrays
  logical, protected :: overide_timestep_check=.false.!do not perform initial dt test
  integer :: mloop !for looping over multiple initial conditions
  !----------------------------warnings---------------------------------------
  integer, private :: max_warn_count=5 !maximum warning count code can survive
  integer, private :: warn_count=0 !number of warnings
  !------------------------------openmp--------------------------------------------
  logical,private :: serial_run=.false.
  integer,private :: qvort_nproc=0
  contains
  !*************************************************************************************************  
  !>read the file run.in obtaining all parameters at runtime, avoiding the need to recompile the code
  subroutine read_run_file()
    implicit none
    ! input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0

    open(fh, file='run.in')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '     ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)
          select case (label)
          case ('nsteps')
             !number of steps, enter a natural number
             read(buffer, *, iostat=ios) nsteps !number of steps to take
          case ('shots')
             !number often to print to file, enter a natural number
             read(buffer, *, iostat=ios) shots !how often to print to file
          case ('recon_shots')
             !how often to try reconnection algorithm, enter a natural number
             read(buffer, *, iostat=ios) recon_shots !how often to perform reconnection algorithm
          case ('pcount')
             !initial number of particles, can be overwritten by intial condition, enter natural number
             read(buffer, *, iostat=ios) init_pcount !initial particle count
          case ('dt')
             !size of timestep, cheked based on delta
             read(buffer, *, iostat=ios) dt !timestep value
          case ('delta')
             !resolution, real number
             read(buffer, *, iostat=ios) delta !spatial resolution
          case ('quant_circ')
             !quatum of circulation, not necessary to set, accepts real #
             read(buffer, *, iostat=ios) quant_circ !quantum of circulation
          case ('corea')
             read(buffer, *, iostat=ios) corea !size of vortex core
          case ('box_size')
             !size of box must be>0, real number
             read(buffer, *, iostat=ios) box_size !size of periodic box
          case ('cylind_r')
             read(buffer, *, iostat=ios) cylind_r !radius of cylinder
          case ('velocity')
             !velocity field options are LIA, BS, Tree
             read(buffer, *, iostat=ios) velocity !BS/LIA/Tree
          case ('boundary_x')
             read(buffer, *, iostat=ios) boundary_x !open/periodic/mirror
          case ('boundary_y')
             read(buffer, *, iostat=ios) boundary_y !open/periodic/mirror
          case ('boundary_z')
             read(buffer, *, iostat=ios) boundary_z !open/periodic/mirror
          case ('tree_theta')
             read(buffer, *, iostat=ios) tree_theta !tree code, opening angle
          case ('tree_print')
             read(buffer, *, iostat=ios) tree_print !print the tree mesh
          case ('tree_extra_correction')
             read(buffer, *, iostat=ios) tree_extra_correction !extra correction suggested by B&H
          case ('normal_velocity')
             read(buffer, *, iostat=ios) normal_velocity !zero/xflow/ABC/KS   
          case ('norm_vel_xflow')
             read(buffer, *, iostat=ios) norm_vel_xflow !counterflow velocity       
          case ('alpha')
             read(buffer, *, iostat=ios) alpha !mutual friction
          case ('initf')
             read(buffer, *, iostat=ios) initf !initial setup of filaments   
          case ('scaling_factor')
             read(buffer, *, iostat=ios) scale_factor !for scaling random loops 
          case ('loop_translate')
             read(buffer, *, iostat=ios) loop_translate !for scaling random loops 
          case ('rotation_factor')
             read(buffer, *, iostat=ios) rotation_factor !for rotating random loops      
          case ('line_count')
             read(buffer, *, iostat=ios) line_count !used in certain intial conditions         
          case ('NAN_test')
             read(buffer, *, iostat=ios) NAN_test !test for NANs
          case ('max_warn_count')
             read(buffer, *, iostat=ios) max_warn_count !maximum number of warning counts we survive
          case ('overide_timestep_check')
             read(buffer, *, iostat=ios) overide_timestep_check !no timestep check 
          case ('serial_run')
             read(buffer, *, iostat=ios) serial_run !do not run in parallel
          case ('qvort_nproc')
             read(buffer, *, iostat=ios) qvort_nproc !specify number of processes
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do
  end subroutine
  !*************************************************************************************************  
  !>print error message to screen and stop the run
  !!provide a location of error and the message to print
  subroutine fatal_error(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    character(len=300) :: message_file
    write (*,*) '-------------------------FATAL ERROR-------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) "FYI: t= ", t, "iteration number= ", itime 
    write (*,*) '-------------------------------------------------------------'
    stop
  end subroutine
  !*************************************************************************************************  
  !>print a warning message to screen - will not stop the code
  !!provide a location of error and the message to print
  subroutine warning_message(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    write (*,*) '-------------------------WARNING----------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) '------------------------------------------------------------'
    !increment warning count
    warn_count=warn_count+1
    if (warn_count>max_warn_count) then
      call fatal_error('warning count','maximum number of warning counts exceeded')
    end if
  end subroutine
  !*************************************************************************************************  
  !>generate a new random seed, or read one in from ./data if restating code
  subroutine init_random_seed()
     !CREATE A NEW RANDOM SEED, UNLESS RESTARTING CODE
     integer :: i, n=16, clock
     integer, allocatable :: seed(:)
     logical :: seed_exists=.false.
     allocate(seed(n))
     
     call system_clock(count=clock)
     
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     !read the seed in from file if possible
     inquire(file='./data/seed.dat',exist=seed_exists)
     if (seed_exists) then
       write(*,*)'reading in random seed from ./data'
       open(unit=37,file='./data/seed.dat',form='unformatted')
         read(37) seed
       close(37)
     else
       write(*,*)'generating new random seed saving to ./data'
       open(unit=37,file='./data/seed.dat',form='unformatted')
         write(37) seed
       close(37)
     end if
     call random_seed(put = seed)
     deallocate(seed)
  end subroutine 
  !*************************************************
  subroutine init_openmp
    implicit none
    write(*,'(a)') ' ------------------------OPENMP-----------------------' 
    if (serial_run) then
      call omp_set_num_threads(1) 
      write(*,'(a)') 'serial run, running on one process'
    else
      if (qvort_nproc>0) then
        call omp_set_num_threads(qvort_nproc) 
      else
        qvort_nproc=omp_get_max_threads()
      end if
      write(*,'(a,i2.2,a)') ' parallel run, running on ', qvort_nproc, ' processes'
    end if
  end subroutine
  !**************************************************************************************************  
end module
