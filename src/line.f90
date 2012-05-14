!>all routines which alter the geometry of the vortex filament/flux tube should
!>be contained in this module.
!>The main routines here insert and remove particles to maintain a roughly
!>constant resolution along the filaments.
module line
  use cdata
  use general
  use boundary
  use reconnection
  contains
  !>insert new points to maintain the resolution of the line,
  !>this will mean the point separation lies between \f$\delta\f$
  !>and \f$\delta/2\f$ to ensure the curvature is not affected by this the new point
  !>\f$ \mathbf{s}_{i'} \f$ is positioned between i, i+1 at position
  !>\f[ \mathbf{s}_{i'}=\frac{1}{2}(\mathbf{s}_i+\mathbf{s}_{i+1})+\left( \sqrt{R^2_{i'}
  !!-\frac{1}{4}\ell_{i+1}^2}-R_{i'} \right)\frac{\mathbf{s}_{i'}''}{|\mathbf{s}_{i'}''|},\f]
  !! where \f$R_{i'}=|\mathbf{s}_{i'}''|^{-1}\f$.
  subroutine pinsert
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: disti, curv, f_ddot(3)
    integer :: par_new
    integer :: old_pcount, i
    !check that we have particles
    if (count(mask=f(:)%infront>0)==0) then
      call fatal_error('line.mod:pinsert','vortex line length is 0, run over')
    end if
    old_pcount=pcount
    total_length=0. !zero this
    do i=1, old_pcount
      if (i>size(f)) then
        !bug check
        print*, 'I think there is a problem' ; exit
      end if
      if (f(i)%infront==0) cycle !empty particle
      if (f(i)%pinnedi) cycle !ignore pinned points infront
      !get the distance between the particle and the one infront
      disti=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
      total_length=total_length+disti !measure total length of filaments
      if (disti>delta) then 
        !we need a new particle 
        !the first step is assess where to put the particle?
        !1. is there an empty slot in out array?
        if (minval(f(:)%infront)==0) then
          !there is an empty slot in the array
          !now find its location minloc is fortran intrinsic function
          par_new=minloc(f(:)%infront,1)
        else
          !2. we must resize the array - 
          !increment by 1 (speed) - use a dummy array tmp
          allocate(tmp(size(f)+1)) ; tmp(:)%infront=0 !0 the infront array
          !copy accross information
          tmp(1:size(f)) = f
          !now deallocate tmp and transfer particles back to f
          !move_alloc is new intrinsic function (fortran 2003)
          call move_alloc(from=tmp,to=f)
          !pcount+1 must be free - so put the particle there
          par_new=pcount+1 ; pcount=pcount+1 !increase the particle count
        end if
        !insert a new particle between i and infront using curvature
        !get second derivative at i
        call get_deriv_2(i,f_ddot) !general.mod
        curv=sqrt(dot_product(f_ddot,f_ddot)) !get the curvature
        if (curv>1E-5) then !if curvature very small linear interpolation OK
          curv=curv**(-1) !actually we want the inverse
          if (curv**2-0.25*disti**2>0.) then
          !could be negative, avoid this
          f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)+&
                       (sqrt(curv**2-0.25*disti**2)-curv)*curv*f_ddot
          else
            !linear interpolation
            f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)
          end if
        else
          !linear interpolation
          f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)
        end if
        !average the current velocity
        f(par_new)%u=0.5*(f(i)%u+f(f(i)%infront)%u) 
        !zero the older velocities
        f(par_new)%u1=0. ; f(par_new)%u2=0.
        !zero the pinned conditions
        f(par_new)%pinnedi=.false. ; f(par_new)%pinnedb=.false.
        f(par_new)%wpinned=0
        !set correct infront and behinds & ghostzones
        f(par_new)%behind=i ; f(par_new)%infront=f(i)%infront
        call get_ghost_p(par_new,f(par_new)%ghosti, f(par_new)%ghostb) !periodic.mod
        f(f(i)%infront)%behind=par_new
        call get_ghost_p(f(i)%infront,f(f(i)%infront)%ghosti, f(f(i)%infront)%ghostb) !periodic.mod         
        f(i)%infront=par_new
        call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb) !periodic.mod         
      end if
    end do
    !calculate average separation of particles
    avg_sep=total_length/(old_pcount-count(mask=f(:)%infront==0))
  end subroutine
  !*************************************************************************
  !>remove points along the filament if they are compressed to the point where
  !>the separation between the point i and i+2 is less than \f$\delta\f$
  !>if phonon emission is set to true in run.in then particles with high 
  !>curvature are removed, smoothing the loop to mimic dissipation at large k
  subroutine premove
    implicit none
    real :: distii
    integer :: infront, tinfront
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      if (f(i)%pinnedi) cycle !ignore pinned points
      infront=f(i)%infront
      if (f(infront)%pinnedi) cycle !ignore one before of pinned points
      !get the distance between the particle and the one twice infront
      distii=distf(i,f(f(i)%infront)%infront)
      !see if we can remove points
      if (distii<.99*delta) then
        !print to file the curvature of this particle
        tinfront=f(f(i)%infront)%infront
        !remove the particle at infront
        f(tinfront)%behind=i ; f(i)%infront=tinfront
        call clear_particle(infront) !general.mod
        !check the size of the new loop
        call loop_killer(i) !reconnection.mod
        remove_count=remove_count+1
      end if
      !also check for two points on the boundary and remove
      if (f(i)%pinnedb.and.f(f(i)%infront)%pinnedi) then
        call clear_particle(f(i)%infront) ; call clear_particle(i) !general.mod
      end if
    end do
  end subroutine
  !******************************************************************
  !>find the closest particle to i using N^2 operation this is
  !>done by looping over all particles and caling distf from general.mod
  !\todo we need to do particles on the boundary as well (periodic)
  subroutine pclose
    implicit none
    integer :: i, j
    real :: dist
    !$omp parallel do private(i,j,dist)
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      f(i)%closestd=100. !arbitrarily high
      !do not test pinned points
      if ((f(j)%pinnedi).or.(f(j)%pinnedb)) cycle
      do j=1, pcount
        if (f(j)%infront==0) cycle !empty particle
        if ((i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j).and. &
        !the above line ensures we do not reconnect with particles infront/behind
            (f(j)%pinnedi.eqv..false.).and.(f(j)%pinnedb.eqv..false.)) then
        !the above line ensures we do not reconnect with pinned points
          dist=distf(i,j)
          if (dist<f(i)%closestd) then
           f(i)%closest=j
           f(i)%closestd=dist
          end if
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine
end module
