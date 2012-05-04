!>The main program, runs the qvort code
program run
  use cdata 
  use initial
  use output
  use boundary
  use timestep
  use line
  use diagnostic
  use reconnection
  use tree
  implicit none
  integer :: i
  call banner_print !output.mod
  call init_random_seed !cdata.mod
  !read in parameters
  call read_run_file !cdata.mod
  call init_setup !initial.mod
  call init_openmp !cdata.mod
  !print dimensions info for matlab plotting
  call print_dims !output.mod
  !begin time loop
  write(*,*) 'setup complete: beginning time loop'
  do itime=nstart, nsteps
    !-------------------------ghost points-------------------------
    call ghostp !boundary.mod
    !---------------------build tree routines----------------------
    if (tree_theta>0) then
      call construct_tree !tree.mod
    end if
    !---------------------velocity operations----------------------
    call pmotion !timestep.mod
    !---------------------line operations--------------------------
    call pinsert !line.mod
    if (mod(itime, recon_shots)==0) then
      if (tree_theta>0) then
        !we may need to empty the tree and then redraw it at this point
        call pclose_tree !tree.mod
      else
        call pclose !line.mod
      end if
      call precon !reconnection.mod
      call premove !line.mod
      call wall_recon
    end if
    !-------------------boundary conditions------------------------
    call enforce_boundary !boundary.mod
    !---------------once all algorithms have been run--------------
    t=t+dt  !increment the time
    !---------------------diagnostic info--------------------------
    call calculate_diagnostics !diagnostics.mod
    !--------------now do all data output--------------------------
    if (mod(itime,shots)==0) then
      call data_dump !output.mod
      call print_info !output.mod
    end if
    !-------------------remove the tree----------------------------
    if (tree_theta>0) then
      call empty_tree(vtree) !empty the tree to avoid a memory leak
      deallocate(vtree%parray) ; deallocate(vtree)
      nullify(vtree) !just in case!
    end if
    !--------------------final sanity checks----------------------
    if (NAN_test) call NAN_finder !general.mod
  end do !close main loop
end program




