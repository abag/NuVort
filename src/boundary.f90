!>all the routines used to enforce boundary conditions should be contained within this 
!!module. Three boundary conditions are supported, open, periodic and solid
module boundary
  use cdata
  use general
  contains 
  !******************************************************************
  !>dummy routine, calls get_ghost_p below, set's the ghost particles for all
  !>points (except empty ones marked with f(i)%infront=0.
  subroutine ghostp
    implicit none
    integer :: i
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particles
      call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb)
    end do
    !$omp end parallel do
  end subroutine
  !******************************************************************
  !>set the ghost particles, essentially these are the positions of the
  !!particles infront/behind, if they are at the other side of the box 
  !!due to periodic b.c. this position must be adjusted we simply subract
  !!(or add) the box size onto their position if the distance between a particle
  !!and it's neighbour is too large to be explained (even by terrbile resolution!)
  subroutine get_ghost_p(i,ginfront,gbehind)
    implicit none
    integer, intent(IN) :: i
    real :: ginfront(3), gbehind(3)
    real :: r_tmp
    if (f(i)%pinnedi) then
      if (sum(abs(f(i)%wpinned))==2) then
        !spherical reflection
        r_tmp=get_radius(f(f(i)%behind)%x)
        ginfront(:)=f(f(i)%behind)%x*(2*cylind_r-r_tmp)/r_tmp
      else if (sum(abs(f(i)%wpinned))==1) then
        !cartesian reflection
        ginfront(:)=f(f(i)%behind)%x(:)+2*abs(f(i)%wpinned)*(0.5*f(i)%wpinned*box_size-f(f(i)%behind)%x(:))
      end if
      gbehind(:)=f(f(i)%behind)%x(:)
    else if (f(i)%pinnedb) then
      ginfront(:)=f(f(i)%infront)%x(:)
      if (sum(abs(f(i)%wpinned))==2) then 
        !spherical reflection
        r_tmp=get_radius(f(f(i)%infront)%x)
        ginfront(:)=f(f(i)%infront)%x*(2*cylind_r-r_tmp)/r_tmp
      else if (sum(abs(f(i)%wpinned))==1) then
        !cartesian reflection
        gbehind(:)=f(f(i)%infront)%x(:)+2*abs(f(i)%wpinned)*(0.5*f(i)%wpinned*box_size-f(f(i)%infront)%x(:))
      end if
    else
      ginfront(:)=f(f(i)%infront)%x(:)
      gbehind(:)=f(f(i)%behind)%x(:)
    end if
    !if periodic then must do more in here
    if (boundary_x=='periodic') then
      !we must ensure that ginfront/gbehind is not on the other side of the box
      !---------------------x------------------------------
      if ((f(i)%x(1)-ginfront(1))>(box_size/2.)) then
        ginfront(1)=ginfront(1)+box_size
      elseif ((f(i)%x(1)-ginfront(1))<(-box_size/2.)) then
        ginfront(1)=ginfront(1)-box_size
      end if
      if ((f(i)%x(1)-gbehind(1))>(box_size/2.)) then
        gbehind(1)=gbehind(1)+box_size
      elseif ((f(i)%x(1)-gbehind(1))<(-box_size/2.)) then
        gbehind(1)=gbehind(1)-box_size
      end if
    end if
    if (boundary_y=='periodic') then
      !---------------------y------------------------------
      if ((f(i)%x(2)-ginfront(2))>(box_size/2.)) then
        ginfront(2)=ginfront(2)+box_size
      elseif ((f(i)%x(2)-ginfront(2))<(-box_size/2.)) then
        ginfront(2)=ginfront(2)-box_size
      end if
      if ((f(i)%x(2)-gbehind(2))>(box_size/2.)) then
        gbehind(2)=gbehind(2)+box_size
      elseif ((f(i)%x(2)-gbehind(2))<(-box_size/2.)) then
        gbehind(2)=gbehind(2)-box_size
      end if
    end if
    if (boundary_z=='periodic') then
      !---------------------z------------------------------
      if ((f(i)%x(3)-ginfront(3))>(box_size/2.)) then
        ginfront(3)=ginfront(3)+box_size
      elseif ((f(i)%x(3)-ginfront(3))<(-box_size/2.)) then
        ginfront(3)=ginfront(3)-box_size
      end if
      if ((f(i)%x(3)-gbehind(3))>(box_size/2.)) then
        gbehind(3)=gbehind(3)+box_size
      elseif ((f(i)%x(3)-gbehind(3))<(-box_size/2.)) then
        gbehind(3)=gbehind(3)-box_size
      end if
    end if
  end subroutine
  !******************************************************************
  !>if we have periodic boundary conditions if a point leaves one side 
  !>of the box, reinsert it on the opposite side. If solid boundaries,
  !>we denote a particle on the boundary with pinnedi or pinnedb being true
  !>due to timestepping the particle will move away from the boundaries, this 
  !>routine simply resets the position.
  subroutine enforce_boundary()
    implicit none
    real :: theta_store
    integer :: i
    integer :: pinned_component
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !-------------x------------------  
      if (boundary_x=='periodic') then   
        if (f(i)%x(1)>(box_size/2.)) then
          f(i)%x(1)=f(i)%x(1)-box_size
        else if (f(i)%x(1)<(-box_size/2.)) then
          f(i)%x(1)=f(i)%x(1)+box_size
        end if
      end if
      !-------------y------------------
      if (boundary_y=='periodic') then
        if (f(i)%x(2)>(box_size/2.)) then
          f(i)%x(2)=f(i)%x(2)-box_size
        else if (f(i)%x(2)<(-box_size/2.)) then
          f(i)%x(2)=f(i)%x(2)+box_size
        end if
      end if
      !-------------z------------------
      if (boundary_z=='periodic') then
        if (f(i)%x(3)>(box_size/2.)) then
          f(i)%x(3)=f(i)%x(3)-box_size
        else if (f(i)%x(3)<(-box_size/2.)) then
          f(i)%x(3)=f(i)%x(3)+box_size
        end if
      end if
      !----------solid boundaries----------
      if (f(i)%pinnedi.or.f(i)%pinnedb) then
        !set the particle position back to the boundary 
        if (sum(abs(f(i)%wpinned))==2) then
          !is this particle 
          theta_store=atan2(f(i)%x(2),f(i)%x(1))
          !set the radius to be cylind_r
          f(i)%x(1)=cylind_r*cos(theta_store)
          f(i)%x(2)=cylind_r*sin(theta_store)
        else if (sum(abs(f(i)%wpinned))==1) then
          !we use maxloc to find out if the particle is either
          !on the x,y or z boundaries
          pinned_component=maxloc(abs(f(i)%wpinned),1)
          !enforcing zero flux at boundaries
          if (f(i)%wpinned(pinned_component)>0) then
            f(i)%x(pinned_component)=0.5*box_size
          else if (f(i)%wpinned(pinned_component)<0) then
            f(i)%x(pinned_component)=-0.5*box_size
          end if
        end if
      end if
    end do
    !$omp end parallel do
  end subroutine
  !*****************************************************************
  !> this routine is used called by initial.mod to create an array of 'shifts'
  !> of the vortex configuration needed for periodic boundary conditions. This
  !> array is looped over in timestep.mod
  subroutine create_periodic_loop_array 
    implicit none
    integer :: i, j, k
    integer :: ii, jj, kk
    integer :: counter=0
    n_periodic=3**n_periodic-1
    allocate(periodic_loop_array(n_periodic,3)) 
    do i=-1,1 ; do j=-1,1 ; do k=-1,1
      if (boundary_x=='periodic') then
        ii=i 
      else
        if (i/=0) cycle
        ii=0
      end if  
      if (boundary_y=='periodic') then
        jj=j
      else 
        if (j/=0) cycle
        jj=0
      end if 
      if (boundary_z=='periodic') then
        kk=k
      else 
        if (k/=0) cycle
        kk=0
      end if
      if (ii==0.and.jj==0.and.kk==0) cycle
      counter=counter+1
      periodic_loop_array(counter,1)=ii
      periodic_loop_array(counter,2)=jj
      periodic_loop_array(counter,3)=kk
    end do ; end do ; end do
  end subroutine
  !*****************************************************************
  !> this routine is used called by initial.mod to create an array of 'reflections'
  !> of the vortex configuration needed for solid boundary conditions. This
  !> array is looped over in timestep.mod
  subroutine create_mirror_loop_array 
    implicit none
    integer :: i, j, k
    integer :: ii, jj, kk
    integer :: counter=0
    n_mirror=3**n_mirror-1
    allocate(mirror_loop_array(n_mirror,3)) 
    do i=-1,1 ; do j=-1,1 ; do k=-1,1
      if (boundary_x=='solid') then
        ii=i 
      else
        if (i/=0) cycle
        ii=0
      end if  
      if (boundary_y=='solid') then
        jj=j
      else 
        if (j/=0) cycle
        jj=0
      end if 
      if (boundary_z=='solid') then
        kk=k
      else 
        if (k/=0) cycle
        kk=0
      end if
      if (ii==0.and.jj==0.and.kk==0) cycle
      counter=counter+1
      mirror_loop_array(counter,1)=ii
      mirror_loop_array(counter,2)=jj
      mirror_loop_array(counter,3)=kk
    end do ; end do ; end do
  end subroutine
end module
