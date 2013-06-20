!> Module which contains initial conditions (loops) for filament (set in run.in) for more 
!!information on the options see \ref INIT
module initial_cond
  use cdata
  use general
  contains
  !*************************************************************************
  !>set up a single loop in the x-y plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$
  subroutine setup_single_loop
    implicit none
    real :: velocity, ring_energy
    real :: radius
    integer :: i 
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    velocity=(quant_circ/(4*pi*radius))*(log(8*radius/corea)-.5)
    ring_energy=0.5*(quant_circ**2)*radius*(log(8*radius/corea)-2.)
    write(*,*) 'initf: single loop, radius of loop:', radius
    write(*,*) 'velocity should be:', velocity
    write(*,*) 'energy should be:', ring_energy
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount)+box_size(2)/2
      f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do   
  end subroutine
  !*************************************************************************
   !>set up a half loop at the x-y plane touching the solid boundary at x=BoxSize/2,
   !>it's size is dictated by the initial number of particles set and the size of \f$\delta\f$
   !>RISTO's first attemp to make a new initial condition.
   subroutine setup_half_loop
     implicit none
     real :: velocity, ring_energy
     real :: radius
     integer :: i
     radius=(0.75*pcount*delta)/pi !75% of potential size
     velocity=(quant_circ/(4*pi*radius))*(log(8*radius/corea)-.5)
     ring_energy=0.5*(quant_circ**2)*radius*(log(8*radius/corea)-2.)
     write(*,*) 'initf: half loop, radius of loop:', radius
     write(*,*) 'velocity should be:', velocity
     write(*,*) 'energy should be:', ring_energy
     !check the boundary conditions are OK
     if (boundary_x=='solid') then
       !loop over particles setting spatial and 'loop' position
       do i=1, pcount
         f(i)%x(1)=0.5*box_size(1)-radius*sin(pi*real(i-1)/(pcount-1))
         f(i)%x(2)=radius*cos(pi*real(i-1)/(pcount-1))
         f(i)%x(3)=0.
         if (i==1) then
           f(i)%pinnedb=.true. ; f(i)%behind=i ; f(i)%wpinned=(/1,0,0/)
           f(i)%infront=i+1
         else if (i==pcount) then
          f(i)%pinnedi=.true. ; f(i)%infront=i ; f(i)%wpinned=(/1,0,0/)
          f(i)%behind=i-1
         else
           f(i)%behind=i-1 ; f(i)%infront=i+1
         end if
         !zero the stored velocities
         f(i)%u1=0. ; f(i)%u2=0.
       end do
     else
      call fatal_error('init.mod:setup_half_loop', &
      'incorrect boundary conditions')     
     end if 
   end subroutine
!*************************************************************************
   !>set up a half loop at the z-y plane touching the solid boundary at z=BoxSize/2,
   !>it's size is dictated by the initial number of particles set and the size of \f$\delta\f$
   subroutine setup_half_loopz
     implicit none
     real :: velocity, ring_energy
     real :: radius
     integer :: i
     radius=(0.75*pcount*delta)/pi !75% of potential size
     velocity=(quant_circ/(4*pi*radius))*(log(8*radius/corea)-.5)
     ring_energy=0.5*(quant_circ**2)*radius*(log(8*radius/corea)-2.)
     write(*,*) 'initf: half loop, radius of loop:', radius
     write(*,*) 'velocity should be:', velocity
     write(*,*) 'energy should be:', ring_energy
     !check the boundary conditions are OK
     if (boundary_z=='solid') then
       !loop over particles setting spatial and 'loop' position
       do i=1, pcount
         f(i)%x(1)=0.
         f(i)%x(2)=radius*cos(pi*real(i-1)/(pcount-1))
         f(i)%x(3)=0.5*box_size(3)-radius*sin(pi*real(i-1)/(pcount-1))
         if (i==1) then
           f(i)%pinnedb=.true. ; f(i)%behind=i ; f(i)%wpinned=(/0,0,1/)
           f(i)%infront=i+1
         else if (i==pcount) then
          f(i)%pinnedi=.true. ; f(i)%infront=i ; f(i)%wpinned=(/0,0,1/)
          f(i)%behind=i-1
         else
           f(i)%behind=i-1 ; f(i)%infront=i+1
         end if
         !zero the stored velocities
         f(i)%u1=0. ; f(i)%u2=0.
       end do
     else
      call fatal_error('init.mod:setup_half_loop', &
      'incorrect boundary conditions')     
     end if 
   end subroutine
  !*************************************************************************
  subroutine setup_random_loops
    implicit none
    real :: loop_radius
    real :: anglex,angley,anglez
    real,dimension(3)::translate
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: loop_size
    integer:: loop_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_random_loops', &
      'you have not set a value for line_count in run.in')
    end if
    if (mod(pcount,line_count)/=0) then
      call fatal_error('init.mod:setup_random_loops', &
      'pcount/line_count is not an integer')
    end if
    loop_size=int(pcount/line_count)
    loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
    write(*,'(a,i5.1,a)') ' drawing ', line_count, ' random loops in the box'
    write(*,'(a,i5.1,a)') ' each loop contains ', loop_size, ' particles'
    write(*,'(a,f7.4)') ' radius of each loop: ', loop_radius
    if (rotation_factor<1.) then
      write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
    end if
    if (scale_factor<1.) then
      write(*,'(a,3f13.4,a)') ' loops occupy ', 100*loop_translate, '% of box volume'
    end if
    do i=1, line_count
      call random_number(anglex)
      call random_number(angley)
      call random_number(anglez)
      call random_number(translate)
      anglex=anglex*2*pi*rotation_factor
      angley=angley*2*pi*rotation_factor
      anglez=anglez*2*pi*rotation_factor
      translate=loop_translate*((box_size*translate-box_size/2.)-loop_radius)
        
      do j=1, loop_size

        loop_position=j+(i-1)*loop_size

        dummy_xp_1(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
        dummy_xp_1(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
        dummy_xp_1(3)=0.

        dummy_xp_2(1)=dummy_xp_1(1)
        dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
        dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

        dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
        dummy_xp_3(2)=dummy_xp_2(2)
        dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

        dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
        dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
        dummy_xp_4(3)=dummy_xp_3(3)
    
        f(loop_position)%x(1)=dummy_xp_4(1)+translate(1)
        f(loop_position)%x(2)=dummy_xp_4(2)+translate(2)
        f(loop_position)%x(3)=dummy_xp_4(3)+translate(3)

        if(j==1) then
          f(loop_position)%behind=i*loop_size
          f(loop_position)%infront=loop_position+1
        else if (j==loop_size) then
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=(i-1)*loop_size+1
        else
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=loop_position+1
        end if
        f(loop_position)%u1=0. ; f(loop_position)%u2=0.
      end do
    end do
  end subroutine
  !*************************************************************************
  !>lines from the lop of the box to the bottom, helical waves added 
  subroutine setup_hayder_wave
    implicit none
    real :: wave_number(5), amp(5), phase(5)
    integer :: pcount_required
    integer :: i, k
    !test run.in parameters, if wrong program will exit
     if (boundary_z=='periodic') then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size(3)/(0.75*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_wave_line', &
      'periodic boundary conditions required')
    end if
    write(*,'(a)') ' drawing a line from -z to +z'
    write(*,'(a)') '5 helical wave pertubations '
    do i=1, pcount
      f(i)%x(1)=0. 
      f(i)%x(2)=0. 
      f(i)%x(3)=box_size(3)/2.-box_size(3)*real(2*i-1)/(2.*pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. 
    end do
    amp=(/5.5,4.5,3.8,4.0,4.5/)
    wave_number=(/5.,4.,3.,2.,1./)
    phase=(/.3422,.6343,.2213,.989342,.5345/)
    !amp=(/1.,0.,0.,0.,0./)
    !wave_number=(/5.,4.,3.,2.,1./)
    !phase=0.
    amp=amp*wave_amp
    phase=2*pi*phase
    do k=1, 5
      print*, k,' wavenumber ',wave_number(k),' amp ', amp(k), ' phase ', phase(k)
      do i=1, pcount
        f(i)%x(1)=f(i)%x(1)+&
                amp(k)*cos(wave_number(k)*2.*pi*real(2.*i-1)/(2.*pcount)+phase(k))
        f(i)%x(2)=f(i)%x(2)+&
                amp(k)*sin(wave_number(k)*2.*pi*real(2.*i-1)/(2.*pcount)+phase(k))
      end do 
    end do
  end subroutine
end module
