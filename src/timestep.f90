!>timestepping and velocity routines are contained in this module
module timestep
  use cdata
  use general
  use normal_fluid
  use tree
  contains
  !*************************************************
  !>implement adams bashforth time-stepping scheme to move particles
  !>\f[
  !! \mathbf{s}_{i}^{n+1}=\mathbf{s}_{i}^{n}+\frac{\Delta t}{12}(23\mathbf{u}_{i}^{n}
  !! -16\mathbf{u}_{i}^{n-1}+5\mathbf{u}_{i}^{n-2})+\mathcal{O}(\Delta t^4)
  !>\f]
  !>adaptive timestep can be set in run.in in which case a comparison between 2nd and 3rd
  !>order scheme is used to estimate error and hence adjust the timestep
  subroutine pmotion()
    implicit none
    real :: u(3) !dummy variable used to store velocities
    integer :: i
    !now loop over all points and get the velocity field
    !$omp parallel do private(i,u)
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call calc_velocity(u,i)
      f(i)%u(:)=u(:) !store the velocity for time-step
    end do
    !$omp end parallel do
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      if (maxval(abs(f(i)%u1))==0) then
        f(i)%x(:)=f(i)%x(:)+dt*f(i)%u(:) !euler
      else if (maxval(abs(f(i)%u2))==0) then
        !first order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+three_twos*dt*f(i)%u(:)-one_half*dt*f(i)%u1(:)
      else 
        !second order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+twenty_three_twelve*dt*f(i)%u(:)-four_thirds*dt*f(i)%u1(:)+five_twelths*dt*f(i)%u2(:)
      end if
      f(i)%u2(:)=f(i)%u1(:)  
      f(i)%u1(:)=f(i)%u(:)
    end do
    !$omp end parallel do
    call randomise_normal_fluid !normal_fluid
  end subroutine
  !*************************************************
  !>get the velocity of each particle subject to the superfluid velocity
  !>plus any normal fluid/forcing
  subroutine calc_velocity(u,i)
    implicit none
    integer, intent(IN) :: i
    real :: u(3), u_norm(3), u_mf(3), u_bs(3), u_mir(3)!velocities
    real :: cov(3), cov_vel(3) !centre of vorticity
    real :: curv, beta !LIA
    real :: hyperviscous_alpha !for hyperviscosity
    real :: f_dot(3), f_ddot(3) !first and second derivs
    integer :: j !used to loop in periodic cases
    !what scheme are we using? (LIA/BS/Tree)
    select case(velocity)
      case('LIA')
        !use the local induction approximation
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        !calculate the curvature
        curv=sqrt(dot_product(f_ddot,f_ddot))
        if (curv<epsilon(0.)) then
          curv=epsilon(0.) !we must check for zero curvature
        else
          curv=1./curv
        end if
        !caluculate beta based on the curvature
        beta=(quant_circ/(4.*pi))*log(4.6*curv/corea)
        u=beta*cross_product(f_dot,f_ddot) !general.mod
      case('BS')
        !full (nonlocal) biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log((1.1/corea)*sqrt(dist_gen(f(i)%x,f(i)%ghosti)*dist_gen(f(i)%x,f(i)%ghostb)))
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        !now we do the non-local part
        !-------------------cartesian box----------------
          u_bs=0. !always 0 before calling the routine
          call biot_savart(i,u_bs)
          u=u+u_bs
          !account for periodic boundary conditions
          !we must shift the mesh in required directions
          u_bs=0. !zero u_bs
          do j=1, n_periodic
            call biot_savart_shift(i,u_bs,periodic_loop_array(j,:))
          end do 
          u=u+u_bs
          !and now solid boundaries via image vortices
          u_mir=0. !zero u_mir
          do j=1, n_mirror
            call biot_savart_mirror(i,u_mir,mirror_loop_array(j,:))
          end do 
          u=u+u_mir
      case('Tree')
        !tree approximation to biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log((1.1/corea)*sqrt(dist_gen(f(i)%x,f(i)%ghosti)*dist_gen(f(i)%x,f(i)%ghostb)))
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        !now walk the tree to get the non-local contribution
          u_bs=0. !zero u_bs
          call tree_walk(i,vtree,(/0.,0.,0./),u_bs) !tree.mod
          u=u+u_bs
          !now account for periodic  boundary conditons
          u_bs=0. !zero u_bs
          do j=1, n_periodic
            call tree_walk(i,vtree,periodic_loop_array(j,:),u_bs)
          end do 
          u=u+u_bs
          !and now solid boundaries via mirrored tree
          u_mir=0. !zero u_mir
          do j=1, n_mirror
            call tree_walk_mirror(i,vtree,mirror_loop_array(j,:),u_mir)
          end do 
          u=u+u_mir
    end select
    f(i)%u_sup=u ! store superfluid velocity
    !-------------------------normal fluid--------------------------
    !check that either of the mutual friction coefficients are >0
    if ((abs(alpha(1))>epsilon(0.)).or.(abs(alpha(2))>epsilon(0.))) then
      if (t<normal_fluid_cutoff) then !cutoff time set in run.in
        if (t<t_zero_normal_fluid) then
          call get_normal_velocity(f(i)%x,u_norm) !normal_fluid.mod
        else
          u_norm=0.
        end if
        ! \todo this could be improved calculating same thing twice
        u_mf=alpha(1)*cross_product(f_dot,(u_norm-u))- &
        alpha(2)*cross_product(f_dot,cross_product(f_dot,(u_norm-u)))
        u=u+u_mf !this way we store the mutual friction velocity
        f(i)%u_mf=u_mf !store it for printing
      end if
    end if
    u=u+super_velocity
  end subroutine
  !**************************************************************************
  !>the desingularised biot savart integral
  !>\f[
  !>\frac{d\mathbf{s}_i}{dt}=\frac{\Gamma}{4\pi} \ln \left(\frac{\sqrt{\ell_i
  !>\ell_{i+1}}}{a}\right)\mathbf{s}_i' \times \mathbf{s}_i'' 
  !>+\frac{\Gamma}{4 \pi} \oint_{\cal L'} \frac{(\mathbf{s}_i-\mathbf{r}) }
  !>{\vert \mathbf{s}_i - \mathbf{r} \vert^3}
  !>\times {\bf d}\mathbf{r}
  !>\f] note the LIA part is calculated in calc_velocity
  subroutine biot_savart(i,u)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j !needed to loop over all particles
    do j=1, pcount
      !check that the particle is not empty/i/f(i)%behind
      if ((f(j)%infront==0).or.(i==j).or.(f(i)%behind==j)) cycle
      if (f(j)%pinnedi) cycle !ignore vorex segment which goes into the wall
      a_bs=distfsq(j,i) !distance squared between i and j
      b_bs=2.*dot_product((f(j)%x-f(i)%x),(f(j)%ghosti-f(j)%x))
      c_bs=dist_gen_sq(f(j)%ghosti,f(j)%x) !distance sqd between j, j+1
      !add non local contribution to velocity vector
      if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
      u_bs=cross_product((f(j)%x-f(i)%x),(f(j)%ghosti-f(j)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of j
    end do
  end subroutine
  !**************************************************************************
  !>as above but shifts the particles (by a vector shift) for periodicity
  subroutine biot_savart_shift(i,u,shift)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables
    real :: shift(3) !moves all the particles for periodic_bc
    integer :: j !needed to loop over all particles
    do j=1, pcount
      !check that the particle is not empty or the particle behind
      !if ((f(j)%infront==0).or.(i==j).or.(f(i)%behind==j)) cycle
      if ((f(j)%infront==0).or.(f(i)%behind==j)) cycle
      if (f(j)%pinnedi) cycle!ignore vorex segment which goes into the wall
      a_bs=dist_gen_sq(f(i)%x,f(j)%x+shift) !distance squared between i and j+shift
      b_bs=2.*dot_product((f(j)%x+shift-f(i)%x),(f(j)%ghosti-f(j)%x))
      c_bs=dist_gen_sq(f(j)%ghosti,f(j)%x) !distance sqd between j, j+1
      !add non local contribution to velocity vector
      if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
      u_bs=cross_product((f(j)%x+shift-f(i)%x),(f(j)%ghosti-f(j)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of j
    end do
  end subroutine
  !**************************************************************************
  !>as above but reflects the particles for image vortices
  subroutine biot_savart_mirror(i,u,reflect)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables
    real :: s_img(3), s_gi_img(3) !dummy variables
    integer  :: reflect(3) !for reflecting 
    integer :: j !needed to loop over all particles
    do j=1, pcount
      !check that the particle is not empty or goes into the wall
      if ((f(j)%infront==0).or.(f(j)%pinnedi)) cycle
      !if pinned behind do not take the contribution from your own image as it 
      if ((sum(f(i)%wpinned-reflect)==0).and.(sum(f(i)%wpinned)<2)) then
        if ((f(i)%pinnedb).and.(i==j)) cycle!has already been used in the local part
        if ((f(i)%pinnedi).and.(j==f(i)%behind)) cycle!has already been used in the local part
      end if
      !we now create dummy variables reflecting f(j)%x and f(j)%ghosti in the required
      !axis and storing the results in s_img, s_gi_img respectively
      s_gi_img=f(j)%x+2*abs(reflect)*(0.5*reflect*box_size-f(j)%x)
      s_img=f(j)%ghosti+2*abs(reflect)*(0.5*reflect*box_size-f(j)%ghosti)
      !note above s_gi_img comes from f(j)%x and s_img from the ghost infront
      !as we have a reflection, hence orientation if flipped
      a_bs=dist_gen_sq(f(i)%x,s_img) !distance squared between i and s_img
      b_bs=2.*dot_product((s_img-f(i)%x),(s_gi_img-s_img))
      c_bs=dist_gen_sq(s_img,s_gi_img) !distance sqd between s_img, s_gi_img
      !add non local contribution to velocity vector
      if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
      u_bs=cross_product((s_img-f(i)%x),(s_gi_img-s_img))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of j
    end do
  end subroutine
end module
