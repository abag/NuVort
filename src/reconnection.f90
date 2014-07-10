!>all reconnection routines are contained in this module this then feeds into line.mod
module reconnection
  use cdata
  use general
  implicit none
  contains
  !******************************************************************
  !>reconnect filaments if they become too close - default scheme
  !>omp is not used here at present as this is a sensitive routine which 
  !>changes other particles
  subroutine precon
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: l_before, l_after !to check line length before and after
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !loop over i, j a helper
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      !now we determine if we can reconnect
      if ((f(i)%closestd<delta/2.).and.(f(i)%closestd>epsilon(1.))) then
        j=f(i)%closest
        !another saftery check
        if (j==pari) cycle ; if (j==parb) cycle ; if (j==0) cycle
        !these two could have reconnected earlier in this case j will be empty
        if (f(j)%infront==0) cycle
        parji=f(j)%infront ; parjb=f(j)%behind
        !we can reconnect based on distance
        !now check whether parallel
        tangent1=norm_tanf(i) ; tangent2=norm_tanf(j) !general.mod
        dot_val=dot_product(tangent1,tangent2) !intrinsic function
        if ((dot_val>0.9)) then
          !we cannot reconnect as filaments are parallel          
        else
          !now we must check the line length before and after
          l_before=dist_gen(f(i)%x,f(i)%ghosti)+dist_gen(f(i)%x,f(i)%ghostb)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(j)%x,f(j)%ghostb)
          l_after=dist_gen(f(i)%x,f(i)%ghostb)+dist_gen(f(i)%x,f(j)%x)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(pari)%x,f(parjb)%x)
          if (l_after<=l_before) then
            !reconnect the filaments
            recon_count=recon_count+1 !keep track of the total # of recons
            !print the time of the recon and its location
            open(unit=61,file='./data/recon_location.log',position='append')
              write(61,*) t, 0.5*(f(i)%x+f(j)%x), acos(dot_val), (l_before-l_after)
            close(61)
            !set correct behind_infront
            f(parjb)%infront=pari
            f(pari)%behind=parjb
            f(i)%infront=j
            f(j)%behind=i
            !check the size of these new loops
            call loop_killer(pari) ; call loop_killer(i)
          end if
        end if 
      end if
    end do
  end subroutine
  !******************************************************************
  !>reconnect filaments if they become too close - default scheme
  !>omp is not used here at present as this is a sensitive routine which 
  !>changes other particles
  subroutine wall_recon
    implicit none
    real :: dist_wall
    integer :: pari, parb
    integer :: i, j !we must loop over all particles
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !perform a quick cleaning operation, we could create pairs of particles
      !which are stuck to the boundaries
      if (f(i)%pinnedb.and.f(f(i)%infront)%pinnedi) then
        !remove these points and cycle
        call clear_particle(f(i)%infront) !obv clear particle infront first
        call clear_particle(i) 
      end if
      if ((f(i)%pinnedi).or.(f(i)%pinnedb)) cycle !obviously don't test pinned points
      !see if we are close to a boundary 
      ! x - boundary
      !------------------------first check cartesian box---------------------
      if (boundary_x=='solid') then
        if (f(i)%x(1)>0.5*(box_size(1)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ; f(pari)%wpinned=(/1,0,0/)
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/1,0,0/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        else if (f(i)%x(1)<-0.5*(box_size(1)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ; f(pari)%wpinned=(/-1,0,0/) 
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/-1,0,0/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        end if 
      end if
      ! y - boundary
      if (boundary_y=='solid') then
        if (f(i)%x(2)>0.5*(box_size(2)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ; f(pari)%wpinned=(/0,1,0/)
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/0,1,0/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        else if (f(i)%x(2)<-0.5*(box_size(2)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ; f(pari)%wpinned=(/0,-1,0/) 
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/0,-1,0/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        end if 
      end if
      ! z - boundary
      if (boundary_z=='solid') then
        if (f(i)%x(3)>0.5*(box_size(3)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ;f(pari)%wpinned=(/0,0,1/) 
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/0,0,1/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        else if (f(i)%x(3)<-0.5*(box_size(3)-delta)) then
          pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
          f(pari)%pinnedb=.true. ; f(pari)%behind=pari ;f(pari)%wpinned=(/0,0,-1/)
          f(parb)%pinnedi=.true. ; f(parb)%infront=parb ; f(parb)%wpinned=(/0,0,-1/)
          call clear_particle(i) !general.mod
          !we must test if we have 'double pinned the particle infront or behind'
          if (f(pari)%pinnedi.and.f(pari)%pinnedb) call clear_particle(pari)
          if (f(parb)%pinnedi.and.f(parb)%pinnedb) call clear_particle(parb)
          wall_recon_count=wall_recon_count+1 !keep track of the total # of wall recons
        end if 
      end if
    end do
  end subroutine
  !**************************************************
  !>removes loops with less than 4 particles this is needed
  !>to ensure derivatives can be calculated correctly pinned points 
  !>are ignored, we remove small loops on the boundary in premove
  subroutine loop_killer(particle)
    implicit none
    integer :: particle, next
    integer :: store_next
    integer :: i, counter
    counter=1
    next=particle
    do i=1, pcount
      !do not checked pinned points we check these in a
      !separate routine
      if (f(next)%pinnedi.or.f(next)%pinnedb) return   
      next=f(next)%infront
      if (next==particle) exit  
      counter=counter+1
    end do
    ! If loop is too small destroy
    if (counter<4) then
      next=particle 
      do i=1, pcount
        store_next=f(next)%infront
        call clear_particle(next) !general.mod
        next=store_next
        if (next==particle) then
          exit  
        end if
      end do
    end if
  end subroutine
end module
