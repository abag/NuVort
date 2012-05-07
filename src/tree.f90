!>all the tree routines used by the code are contained within this module,
!>pointers are used to create a linked list, be very careful to avoid memory
!>leaks!
module tree
   ! CONTAINS ALL THE ROUTINES USED TO EVALUATE THE BIOT-SAVART INTEGRAL
   ! USING TREE CODE METHODS - A WARNING, POINTERS ARE USED HEAVILY HERE
   ! SO BE VERY CAREFUL, MEMORY LEAKS ARE A DISTINCT POSSIBILITY!
   use cdata
   use general
   !>the tree structure, uses pointers
   !>@param posx @param posy @param posz bottem left corner of box
   !>@param width the width of the mesh 
   !>@param fbl child cell front bottom left, @param fbr child cell front bottom right
   !>@param ftl child cell front top left @param ftr child cell front top right
   !>@param bbl child cell back bottom left, @param bbr child cell back bottom right
   !>@param btl child cell back top left @param btr child cell back top right
   !>@param all particles in the cell
   !>@param centx @param centy @param centz the center of 'mass' of the cell
   !>@param circ the total circulation (vector) of the cell
   type node
      real :: posx,posy,posz 
      real :: width 
      integer :: pcount 
      type (node), pointer :: fbl, fbr !front (y) bottom left, right (child)
      type (node), pointer :: ftl, ftr !front (y) top left, right (child)
      type (node), pointer :: bbl, bbr !back (y) bottom left, right (child)
      type (node), pointer :: btl, btr !back (y) top left, right (child)
      type (qvort), allocatable :: parray(:) 
      real :: centx, centy, centz 
      real :: vec_cent(3) !center in vector form
      real :: circ(3) !total circulation vector
   end type node
   !>the vortex point tree
   type (node), pointer :: vtree 
   !>how many evaluations are required to calculate the tree vel field
   integer, protected :: eval_counter
  contains
  !>An array to set the main node of the tree and then call the recurisive 
  !>subroutine build_tree this is called in the main program
  subroutine construct_tree()
    allocate(vtree)
    allocate(vtree%parray(pcount))
    !set up the main node, not a 100% sure if this needs doing as will never 
    !be used for evaluations
    vtree%posx=-box_size/2. ; vtree%posy=-box_size/2. ; vtree%posz=-box_size/2.
    vtree%width=box_size
    vtree%pcount=pcount
    vtree%parray=f
    vtree%centx=sum(f(:)%x(1), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    vtree%centy=sum(f(:)%x(2), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    vtree%centz=sum(f(:)%x(3), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    vtree%vec_cent=(/vtree%centx,vtree%centy,vtree%centz/)
    vtree%circ=0. !0 the circulation
    !now nullify its children - for safety
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !the recursive routine which builds the tree
  end subroutine
  !**************************************************************************
  !>called by above, build the tree, continue to subdivide until a cell is empty
  !>or contains only one particle - octree structure used each new cell is created
  !>by the routine below new_tree
  recursive subroutine build_tree(vtree)
    implicit none
    type(node), pointer :: vtree
    character (len=40) :: print_file
    if (vtree%pcount>1) then
      !divide up into 8 children
      !front, bottom, left (x,y,z)
      call new_tree(vtree%fbl, vtree%posx, vtree%posy, vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, bottom, right (x+w,y,z)
      call new_tree(vtree%fbr, vtree%posx+(vtree%width/2), vtree%posy, vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, top, left (x,y,z+w)
      call new_tree(vtree%ftl, vtree%posx, vtree%posy,vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, top, right (x+w,y,z+w)
      call new_tree(vtree%ftr, vtree%posx+(vtree%width/2), vtree%posy,vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, bottom, left (x,y+w,z)
      call new_tree(vtree%bbl, vtree%posx, vtree%posy+(vtree%width/2), vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, bottom, right (x+w,y+w,z)
      call new_tree(vtree%bbr, vtree%posx+(vtree%width/2), vtree%posy+(vtree%width/2), vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, top, left (x,y+w,z+w)
      call new_tree(vtree%btl, vtree%posx, vtree%posy+(vtree%width/2), vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, top, right (x+w,y+w,z+w)
      call new_tree(vtree%btr, vtree%posx+(vtree%width/2), vtree%posy+(vtree%width/2), vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
    else
      !print to file to test the algorithm works
      if ((mod(itime,shots)==0).and.tree_print) then
        write(unit=print_file,fmt="(a,i3.3,a)")"./data/tree_mesh",itime/shots,".log"
        open(unit=97,file=print_file,action='write', position='append')
          write(97,*) vtree%posx, vtree%posx+vtree%width, vtree%posy, &
                     vtree%posy+vtree%width, vtree%posz, vtree%posz+vtree%width, vtree%pcount
        close(97)
      end if
    end if
  end subroutine
  !********************************************************************************
  !>allocate a new cell, populate with particles and calculate centre of mass, 
  !>total circulation vector
  subroutine new_tree(vtree,posx,posy,posz,width, loc_pcount,parray)
    !create a new tree
    implicit none
    integer, intent(IN) :: loc_pcount
    real, intent(IN) :: posx, posy, posz, width
    type(qvort) :: parray(loc_pcount)
    type(node), pointer :: vtree
    integer :: i, counter
    allocate(vtree)
    allocate(vtree%parray(loc_pcount))
    !set up this node of the tree with the input info
    vtree%posx=posx ; vtree%posy=posy ; vtree%posz=posz
    vtree%width=width ; vtree%centx=0. ; vtree%centy=0.
    vtree%centz=0. ; vtree%circ=0. 
    vtree%parray(:)%infront=0 !0 the infront marker
    counter=0 !0 the counter
    do i=1, loc_pcount
      !empty points are not part of the tree
      if (parray(i)%infront==0) cycle
      !neither are points which go into the wall
      if (parray(i)%pinnedi) cycle !this marks vorex segment which goes into the wall
      if((parray(i)%x(1)>=posx).and.(parray(i)%x(1)<(posx+width)).and. &
         (parray(i)%x(2)>=posy).and.(parray(i)%x(2)<(posy+width)).and. &
         (parray(i)%x(3)>=posz).and.(parray(i)%x(3)<(posz+width))) then
        counter=counter+1
        vtree%parray(counter)=parray(i)
        !add the position of the particle to the position of the mesh
        vtree%centx=vtree%centx+parray(i)%x(1)
        vtree%centy=vtree%centy+parray(i)%x(2)
        vtree%centz=vtree%centz+parray(i)%x(3)
        !also account for the circulation
        vtree%circ(1)=vtree%circ(1)+parray(i)%ghosti(1)-parray(i)%x(1)
        vtree%circ(2)=vtree%circ(2)+parray(i)%ghosti(2)-parray(i)%x(2)
        vtree%circ(3)=vtree%circ(3)+parray(i)%ghosti(3)-parray(i)%x(3)
      end if
    end do
    !take the average of the positions
    if (counter/=0) then
     vtree%centx=vtree%centx/counter
     vtree%centy=vtree%centy/counter
     vtree%centz=vtree%centz/counter
    end if
    !store in vector form
    vtree%vec_cent=(/vtree%centx,vtree%centy,vtree%centz/)
    vtree%pcount=counter
    !nullify the children
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !recall the routine that builds the tree
  end subroutine
  !********************************************************************************
  !>empty out the tree structure to avoid a memory leak - repeatedly calls dead_tree
  recursive subroutine empty_tree(vtree)
    implicit none
    type (node), pointer :: vtree
    if(vtree%pcount>1) then
      !kill all it's 8 children
      call dead_tree(vtree%fbl)
      call dead_tree(vtree%fbr)
      call dead_tree(vtree%ftl)
      call dead_tree(vtree%ftr)
      call dead_tree(vtree%bbl)
      call dead_tree(vtree%bbr)
      call dead_tree(vtree%btl)
      call dead_tree(vtree%btr)
    end if
  end subroutine
  !********************************************************************************
  !>deallocate a tree
  subroutine dead_tree(vtree)
    !chop down that tree (metaphorically) 
    implicit none
    type (node), pointer :: vtree
    call empty_tree(vtree) !kill its children (if it has any)
    !print*, vtree%pcount
    deallocate(vtree%parray)
    deallocate(vtree)
    nullify(vtree)
  end subroutine
  !********************************************************************************
  !>find the closest particle using the tree structure, theoretically should be NlogN
  !>repeatedly calls closest_tree
  subroutine pclose_tree
    !loop over all the particles and find the closest one to i using the tree mesh
    implicit none
    integer :: i
    !$omp parallel do private(i)
    do i=1, pcount
      f(i)%closestd=100. !arbitrarily high
      if (f(i)%infront==0) cycle !empty particle
      !do not test if you are pinned 
      if (f(i)%pinnedi.or.f(i)%pinnedb) cycle 
      call closest_tree(i,vtree)
    end do
    !$omp end parallel do
  end subroutine
  !********************************************************************************
  !>recurisve routine used by pclose_tree to find the closest particle to i
  recursive subroutine closest_tree(i,vtree)
    !find the nearest particle to i
    implicit none
    integer, intent(IN) :: i
    type (node), pointer :: vtree
    real :: theta, dist !tree distance
    integer :: j
     if (f(i)%infront==0) return !zero particle exit
     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((f(i)%x(1)-vtree%centx)**2+(f(i)%x(2)-vtree%centy)**2+(f(i)%x(3)-vtree%centz)**2)
     theta=vtree%width/dist
     if (vtree%pcount==1.or.theta<.5) then
       !see if the distance if less than the minimum distance
       if (dist<f(i)%closestd) then
         j=f(vtree%parray(1)%infront)%behind
         if (vtree%pcount>1) then
           !skip this the particle must be isolated
           !probably will not happen in practice
           f(i)%closest=0
         else if ((i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j)) then
           !test to see if j is pinned 
           !note we have excluded pinnnedi when building the tree 
           !but I leave this here anyway
           if (f(j)%pinnedi.or.f(j)%pinnedb) then
             !do nothing its pinned
           else
             f(i)%closest=j ; f(i)%closestd=dist !OK not pinned
           end if 
         end if
       end if
     else
       !open the box up and use the child cells
       call closest_tree(i,vtree%fbl) ; call closest_tree(i,vtree%fbr)
       call closest_tree(i,vtree%ftl) ; call closest_tree(i,vtree%ftr)
       call closest_tree(i,vtree%bbl) ; call closest_tree(i,vtree%bbr)
       call closest_tree(i,vtree%btl) ; call closest_tree(i,vtree%btr)
     end if
   end subroutine
   !************************************************************************
   !>get the tree code approximation to the biot savart integral, giving
   !>the induced velocity at particle i due to all other vortices, accepts 
   !>an arguement shift, which allows you to shift the position of the
   !>vortices for periodic boundary conditions
   recursive subroutine tree_walk(i,vtree,shift,u)
     implicit none
     integer, intent(IN) :: i
     real, intent(IN) :: shift(3) !shift tree (periodicity)
     type (node), pointer :: vtree !the tree
     real :: u(3) !the velocity at particle i
     real :: vect(3) !helper vector
     real :: dist, theta !helper variables
     real :: geo_dist !distance between centre of vorticity and centre of cell
     real :: a_bs, b_bs, c_bs, u_bs(3) !helper variables 
     integer :: j=0 !the particle in the tree mesh
     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((f(i)%x(1)-(vtree%centx+shift(1)))**2+&
               (f(i)%x(2)-(vtree%centy+shift(2)))**2+&
               (f(i)%x(3)-(vtree%centz+shift(3)))**2)
     if (tree_extra_correction) then
       geo_dist=sqrt((vtree%centx-(vtree%posx+vtree%width/2.))**2+&
                     (vtree%centy-(vtree%posy+vtree%width/2.))**2+&
                     (vtree%centz-(vtree%posz+vtree%width/2.))**2)
       if (dist-geo_dist>0.) then
         theta=vtree%width/(dist-geo_dist)
       else
         theta=vtree%width/dist 
       end if
     else
       theta=vtree%width/dist 
     end if
     if (vtree%pcount==1.or.theta<tree_theta) then
       !use the contribution of this cell
       if (vtree%pcount==1) then!this will be particles very close to i
         !check that the particle is not itself or the particle behind
         j=f(vtree%parray(1)%infront)%behind
         if (maxval(abs(shift))<epsilon(0.)) then
           !if the box has not been shifted (i.e. not periodic wrap) then
           !we should not take the contribution from the desingularised part
           if (j==i) return
           if (j==f(i)%behind) return
         end if
         if (dist<epsilon(0.)) then
           call fatal_error('tree.mod:tree_walk', & 
           'singularity in BS (tree) velocity field - &
           this is normally caused by having recon_shots too large') !cdata.mod
         end if
       end if
       eval_counter=eval_counter+1 !increment this counter
       vect(1)=((vtree%centx+shift(1))-f(i)%x(1)) 
       vect(2)=((vtree%centy+shift(2))-f(i)%x(2)) 
       vect(3)=((vtree%centz+shift(3))-f(i)%x(3))
       a_bs=dist**2
       b_bs=2.*dot_product(vect,vtree%circ)
       c_bs=dot_product(vtree%circ,vtree%circ) 
       if (4*a_bs*c_bs-b_bs**2<epsilon(0.)) return !avoid 1/0
       u_bs=cross_product(vect,vtree%circ)
       u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
       u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
       u=u+u_bs
     else
       !open the box up and use the child cells
       call tree_walk(i,vtree%fbl,shift,u) ; call tree_walk(i,vtree%fbr,shift,u)
       call tree_walk(i,vtree%ftl,shift,u) ; call tree_walk(i,vtree%ftr,shift,u)
       call tree_walk(i,vtree%bbl,shift,u) ; call tree_walk(i,vtree%bbr,shift,u)
       call tree_walk(i,vtree%btl,shift,u) ; call tree_walk(i,vtree%btr,shift,u)
     end if 
   end subroutine
   !************************************************************************
   !>get the tree code approximation to the mirrored biot savart integral, giving
   !>the induced velocity at particle i due to all the image vortices reflected in a
   !> specified plane, reflect
   recursive subroutine tree_walk_mirror(i,vtree,reflect,u)
     implicit none
     integer, intent(IN) :: i !the particle we are interested in
     integer, intent(IN) :: reflect(3) !reflect the tree (solid bounds.)
     type (node), pointer :: vtree !the tree
     real :: u(3) !the velocity at particle i
     real :: vect(3) !helper vector
     real :: dist, theta !helper variables
     real :: geo_dist !distance between centre of vorticity and centre of cell
     real :: a_bs, b_bs, c_bs, u_bs(3) !helper variables 
     integer :: j=0 !the particle in the tree mesh
     real :: reflect_cent(3) !helper for reflected center of mass
     real :: s_img(3), s_gi_img(3) !dummy variables
     if (vtree%pcount==0) return !empty box no use
     !before we reflect get the distance between the geometric centre of the cell
     !and the centre of mass
     geo_dist=sqrt((vtree%centx-(vtree%posx+vtree%width/2.))**2+&
                   (vtree%centy-(vtree%posy+vtree%width/2.))**2+&
                   (vtree%centz-(vtree%posz+vtree%width/2.))**2)
     !reflect the centre of mass in the appropriate plane
     s_img=vtree%vec_cent+vtree%circ
     s_gi_img=vtree%vec_cent
     s_gi_img=s_gi_img+2*abs(reflect)*(0.5*reflect*box_size-s_gi_img)
     s_img=s_img+2*abs(reflect)*(0.5*reflect*box_size-s_img)
     dist=dist_gen(f(i)%x,s_img)
     !make the correction suggested by barnes and hut?
     if (tree_extra_correction) then
       if (dist-geo_dist>0.) then
         theta=vtree%width/(dist-geo_dist)
       else
         theta=vtree%width/dist 
       end if
     else
       theta=vtree%width/dist 
     end if
     theta=vtree%width/dist 
     if (vtree%pcount==1.or.theta<tree_theta) then
       !use the contribution of this cell
       if (vtree%pcount==1) then !this will be particles very close to i
         !check that the particle can be used
         j=f(vtree%parray(1)%infront)%behind
         !the line below is probably not needed as particles which go into 
         !the wall are not included when building the tree
         if (f(j)%pinnedi) return!ignore vorex segment which goes into the wall
         !if pinned behind do not take the contribution from your own image as it 
         if ((sum(f(i)%wpinned-reflect)==0).and.(sum(f(i)%wpinned)<2)) then
           if ((f(i)%pinnedb).and.(i==j)) return!has already been used in the local part
           if ((f(i)%pinnedi).and.(j==f(i)%behind)) return!has already been used in the local part
         end if
         if (dist<epsilon(0.)) then
           call fatal_error('tree.mod:tree_walk', & 
           'singularity in BS (tree) velocity field - &
           this is normally caused by having recon_shots too large') !cdata.mod
         end if
       end if
       a_bs=dist**2
       b_bs=2.*dot_product((s_img-f(i)%x),(s_gi_img-s_img))
       c_bs=dist_gen_sq(s_img,s_gi_img) !distance sqd between s_img, s_gi_img
       if (4*a_bs*c_bs-b_bs**2<epsilon(0.)) return !avoid 1/0
       u_bs=cross_product((s_img-f(i)%x),(s_gi_img-s_img))
       u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
       u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
       u=u+u_bs
     else
       !open the box up and use the child cells
       call tree_walk_mirror(i,vtree%fbl,reflect,u) ; call tree_walk_mirror(i,vtree%fbr,reflect,u)
       call tree_walk_mirror(i,vtree%ftl,reflect,u) ; call tree_walk_mirror(i,vtree%ftr,reflect,u)
       call tree_walk_mirror(i,vtree%bbl,reflect,u) ; call tree_walk_mirror(i,vtree%bbr,reflect,u)
       call tree_walk_mirror(i,vtree%btl,reflect,u) ; call tree_walk_mirror(i,vtree%btr,reflect,u)
     end if 
   end subroutine
end module
