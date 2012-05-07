!> Module which contains all the routines/call to routines which setup and
!!restart the code.
module initial
  use cdata
  use initial_cond
  use boundary
  use output
  use normal_fluid
  use boundary
  contains
  !*************************************************************************
  !>Prints and sets up intial conditions - will give warnings/errors if there
  !>are any conflicting options set in run.in
  subroutine init_setup()
    implicit none
    logical :: restart !can we restart or not?
    write(*,'(a)') ' ---------------------VORTEX PARAMETERS------------------' 
    write(*,'(a,f9.7)') ' quantum of circulation is:', quant_circ 
    write(*,'(a,e9.3)') ' core size is:', corea 
    !check particle separation has been set
    if (delta<epsilon(0.)) call fatal_error('init.mod','delta must be set in run.in')
    !check that the particle count (pcount) has been set
    if (init_pcount<5) then
      call fatal_error('init.mod','you must set enough intial particles')
    end if 
    !we must check the timestep is sufficient to resolve the motion
    !based on the smallest separation possible in the code
    write(*,'(a)') ' ---------------------TIME-STEP--------------------'
    if (overide_timestep_check) then
      write(*,'(a,e10.4,a)') ' forcing dt to be: ', dt, ' overiding timestep check'
    else 
      call timestep_check !initial.mod
    end if
    write(*,'(a,i3.2,a)') ' outputting filament information every ', shots, ' time-steps'
    !periodic bounday conditions?
    write(*,'(a)') ' ---------------------BOUNDARY CONDITIONS--------------------' 
    if (box_size>0.) then
      write(*,'(a,f8.3)') ' box size:', box_size
      !what is the boundary
      select case(boundary_x)
        case('periodic')
          write(*,'(a)') ' boundary x: periodic'
          n_periodic=n_periodic+1
        case('open')
          write(*,'(a)') ' boundary x: open'
        case('solid')
          write(*,'(a)') ' boundary x: solid'
          n_mirror=n_mirror+1
        case('cylind')
          !check the y boundary is also set to cylinder
          if (boundary_y/='cylind') call fatal_error('init_setup:',&
          'boundary_y must also be set to cylind')
          !now check that the cylinder radius has been set
          if (cylind_r>epsilon(0.)) then 
            write(*,'(a)') ' boundary x/y: running with cylindrical boundary'
            write(*,'(a,f9.4)') ' radius of cylinder: ', cylind_r
            !to speed up logical evaluations if we are running with cylindrical 
            !boundaries we activate the following logical variable
            cylindrical_boundaries=.true.
          else
            call fatal_error('init_setup:', 'cylinder radius not set')
          end if
        case default
          call fatal_error('init_setup:', 'incorrect boundary_x parameter')
      end select
      select case(boundary_y)
        case('periodic')
          write(*,'(a)') ' boundary y: periodic'
          n_periodic=n_periodic+1
        case('open')
          write(*,'(a)') ' boundary y: open'
        case('solid')
          write(*,'(a)') ' boundary y: solid'
          n_mirror=n_mirror+1
        case('cylind')
          !check the x boundary is also set to cylinder
          if (boundary_x/='cylind') call fatal_error('init_setup:',&
          'boundary_x must also be set to cylind')
        case default
          call fatal_error('init_setup:', 'incorrect boundary_y parameter')
      end select
      select case(boundary_z)
        case('periodic')
          write(*,'(a)') ' boundary z: periodic'
          n_periodic=n_periodic+1
        case('open')
          write(*,'(a)') ' boundary z: open'
        case('solid')
          write(*,'(a)') ' boundary z: solid'
          n_mirror=n_mirror+1
        case default
          call fatal_error('init_setup:', 'incorrect boundary_z parameter')
      end select
      !now create the periodic_loop_array
      call create_periodic_loop_array !boundary.mod
      call create_mirror_loop_array !boundary.mod
    else
      call fatal_error('init_setup:', 'box size is less than or equal zero')
    end if
    write(*,'(a)') ' --------------------INITIAL CONDITIONS--------------------' 
    !check if we can restart the code
    inquire(file="./data/var.dat", exist=restart)
    if (restart) then
      call data_restore !init.mod
    else
      pcount=init_pcount
      allocate(f(pcount)) !main vector allocated
      !choose the correct setup routine based on the value of initf in run.in
      select case(initf)
        case('single_loop')
          call setup_single_loop !initial_cond.mod
        case('random_loops')
          call setup_random_loops !initial_loop.mod
        case('half_loop')
          call setup_half_loop !initial_loop.mod
        case default
          call fatal_error('cdata.mod:init_setup', &
                         'invalid choice for initf parameter') !cdata.mod
      end select
    end if
    !enforce boundary conditions
    call enforce_boundary
    !print initial filament to file
    call initial_printf !output.mod
    write(*,'(a)') ' ---------------------VELOCITY CALCULATION----------------------' 
    !print information about the velocity field to screen
    select case(velocity)
      case('LIA')
        write(*,*) 'using local induction approximation - scales like O(N)'
      case('BS')
        write(*,*) 'using full Biot-Savart integral - scales like O(N^2)'
      case('Tree')
        write(*,*) 'using tree approximation to Biot-Savart integral - scales like O(NlogN)'
     case default
       print*, 'correct value for velocity in run.in has not been set'
       print*, 'options are: LIA, BS, Tree'
       call fatal_error('init.mod:init_setup', & 
        'correct value for "velocity" in run.in has not been set') !cdata.mod
    end select
    write(*,'(a)') ' ---------------------Tree Algorithm----------------------' 
    if (tree_theta>0) then
      !check that the boundary conditions are not open or the code should not run, unless we adapt it!
      if ((boundary_x=='open').or.(boundary_y=='open').or.(boundary_z=='open')) then
        call fatal_error('init.mod:init_setup', & 
          'tree approximation cannot run with open boundaries in its current form') !cdata.mod
      end if
      write(*,*) 'using tree algorithms for reconnection routine - scales like O(NlogN)'
      if (tree_extra_correction) write(*,*) 'making extra correction to tree algorithm'
    else
      write(*,*) 'using brute force reconnection routine - scales like O(N^2)'
    end if
    write(*,'(a)') ' --------------------NORMAL FLUID--------------------' 
    call setup_normal_fluid !normal_fluid.mod
  end subroutine
  !**********************************************************************
  !>restart the code code periodically writes all the main variables to a file
  !>./data/dump.dat, if this is detected at startup then the code will restart
  !>from that data file
  subroutine data_restore
    !restart the code
    implicit none
    integer :: dummy_itime 
    open(unit=63,file="./data/var.dat",FORM='unformatted')
      read(63) pcount
      read(63) recon_count
      read(63) dummy_itime
      read(63) t
      allocate(f(pcount))
      read(63) f
      write(*,*) 'restored vortex filament'
    close(63)
    write(*,*) 'data read in from dump file at t=', t
    nstart=dummy_itime+1
  end subroutine
  !****************************************************************
  !> check the timestep is OK if we have a vortex filament only
  subroutine timestep_check
    implicit none
    real :: delta_min, dt_max
    delta_min=delta/2.
    dt_max=((delta_min)**2)/(quant_circ*log(delta_min*1E8/pi))
    if (dt<dt_max) then
      write(*,'(a,e10.4)') ' dt is below maximum possible dt:', dt_max
    else
      write(*,'(a,e10.4)') ' warning set dt below ', dt_max
      call fatal_error('initial.mod:timestep_check','dt is too large')
    end if
  end subroutine
end module

