!>normal fluid component in the equation of motion, if the filament acts as a vortex
!>velocity enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$
!>if the filament is a material line then the normaly velocity is used directly.
!>for further information on the normal fluid options see \ref NF
module normal_fluid
  use cdata
  use general
    real, private :: urms_norm
    !>normal fluid mesh - all grid based normal velocities use this (e.g. Navier Stokes future)
    !>@param x the position on the grid
    !>@param u the velocity
    !>@param v a scalar used in the compressible flow case
    !>@param div the divergence of the velocity field at x
    type norm_fluid_grid
      private
      real :: x(3)
      real :: u(3)
      real :: v
      real :: div
    end type
    !use the abbreviation nfm (normal fluid mesh)
    !>the size of the normal fluid mesh, put in run.in soon
    integer, private, parameter :: nfm_size=128
    real, private :: nfm_res(3), nfm_inv_res(3) !resolution/inv resolution of normal fluid mesh
    !>the normal fluid mesh
    type(norm_fluid_grid), allocatable, private :: nfm(:,:,:)
    real, allocatable, private :: laizet_u(:,:,:,:)
    real, allocatable, private :: xflow_noise(:)
    integer, private :: xflow_noise_size=100000
    integer, private :: xflow_noise_counter=1
    real, allocatable, private :: xp(:),yp(:),zp(:)
    contains 
    !************************************************************
    !>setup everything needed to use a normal fluid this may seem redundant
    !>at the moment but useful for more complicated normal velocities
    subroutine setup_normal_fluid
      !in here we print to file the timescale of the flow
      implicit none
      write(*,*) 'normal fluid velocity field is: ', trim(normal_velocity)
      select case(normal_velocity)
        case('xflow')
          write(*,'(a,f6.3)') ' u(x)=', norm_vel_xflow
        case('shear_xflow')
          write(*,'(a,f6.3,a,f6.3)') ' A= ', norm_vel_xflow,' w= ', norm_shear_omega
          call setup_gen_normalf !normal_fluid.mod
        case('channel')
          write(*,*) 'channel flow with scaling amplitude A'
          write(*,'(a,f6.3,a,f6.3)') ' A= ', norm_vel_xflow
          call setup_gen_normalf !normal_fluid.mod
        case('noisey_xflow')
          write(*,*) 'noisey counter flow with scaling amplitude A'
          write(*,'(a,f6.3,a,f6.3)') ' A= ', norm_vel_xflow
          write(*,*) 'loading in noise from:  ', noise_file
          call load_noise  !normal_fluid.mod
          call setup_gen_normalf !normal_fluid.mod
        case('ychannel')
          write(*,*) 'y channel flow in x direction with scaling amplitude A'
          write(*,'(a,f6.3,a,f6.3)') ' A= ', norm_vel_xflow
          call setup_gen_normalf !normal_fluid.mod
        case('laizetDNS')
          write(*,*) "DNS from Sylvain's code"
          call  setup_laizet
          call setup_gen_normalf !normal_fluid.mod
        case('zero')
        case default
          call fatal_error('normal_fluid.mod:get_normal_fluid', &
          'correct parameter for normal_veloctity not set')
      end select
      write(*,'(a,f6.4,a,f6.4)') ' mutual friction coefficients alpha=',alpha(1),' alpha`=',alpha(2)
      if (normal_fluid_cutoff<1000.) then
        write(*,'(a,f6.3)') ' normal fluid is turned off (and T=0) when t=',normal_fluid_cutoff
      end if
      if (t_zero_normal_fluid<1000.) then
        write(*,'(a,f6.3)') ' normal fluid is zeroed (finite temp. remains) when t=',t_zero_normal_fluid
      end if
    end subroutine
    !************************************************************
    subroutine get_normal_velocity(x,u)
      implicit none
      real, intent(IN) :: x(3) !position of the particle
      real, intent(OUT) :: u(3) !velocity at x
      real :: r, phi, theta !used to convert to polar coords
      real :: u_r, u_theta !used to convert to polar coords
      real :: cossum, ak !helper variable for channel flow
      integer :: peri, perj, perk !used to loop in periodic cases
      integer :: ii ! for looping
      u=0. ! a safety check , 0 the field before we begin!
      select case(normal_velocity)
        case('zero')
          u=0. ! no flow
        case('xflow')
          u=0. ; u(1)=norm_vel_xflow !flow in the x direction
        case('noisey_xflow')
          u=0. ; u(1)=xflow_noise(xflow_noise_counter)*norm_vel_xflow
          !print*, u,xflow_noise_counter
        case('shear_xflow')
          u=0. ; u(1)=norm_vel_xflow*sin(2*pi*x(3)/box_size(3)+t*norm_shear_omega) !flow in the x direction
        case('ychannel')
          u=0. ; u(1)=norm_vel_xflow*(1.-(2.*x(2)/box_size(2))**2) !flow in the x direction
        case('channel')
          u=0. ; cossum=0.
          do ii=1, 10
            ak=(2*ii-1)*(pi/2);
            cossum=cossum+4*((-1)**ii/(ak**3))*(cosh(ak*2*x(2)/box_size(2))*cos(ak*2*x(3)/box_size(3)))/cosh(ak)
          end do
          u(1)=norm_vel_xflow*(1.-(2*x(3)/box_size(3))**2+cossum)*1.2;
        case('laizetDNS')
          call laizet_interp(x,u)
      end select
      u=u+super_velocity
    end subroutine    
    !**********************************************************
    !>set up simple normal fluid models (ones with an analytic form)
    !>and print this field (discretized on a cubic mesh) to a binary file
    subroutine setup_gen_normalf
      implicit none
      integer :: i, j, k
      !print dimensions to file for matlab
      open(unit=77,file='./data/nfm_dims.log',status='replace')
        write(77,*) nfm_size
      close(77)
      allocate(nfm(nfm_size,nfm_size,nfm_size))
      nfm_res=(real(box_size)/nfm_size)
      nfm_inv_res=1./nfm_res
      !$omp parallel do private(i,j,k)
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        nfm(k,j,i)%x(1)=nfm_res(1)*real(2*i-1)/2.-(box_size(1)/2.)
        nfm(k,j,i)%x(2)=nfm_res(2)*real(2*j-1)/2.-(box_size(2)/2.)
        nfm(k,j,i)%x(3)=nfm_res(3)*real(2*k-1)/2.-(box_size(3)/2.)
      end do ; end do ; end do
      !$omp end parallel do
      urms_norm=0. !0 the root mean squared velocity
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        !get the velocity field - shearing wave
        call get_normal_velocity(nfm(k,j,i)%x,nfm(k,j,i)%u)
        urms_norm=urms_norm+(nfm(k,j,i)%u(1)**2+nfm(k,j,i)%u(2)**2+nfm(k,j,i)%u(3)**2)
      end do ; end do ; end do
      urms_norm=sqrt(urms_norm/(nfm_size**3))
      print*, urms_norm
      write(*,'(a)') ' velocity field calculated, printing to ./data/norm_init_mesh.dat'
      open(unit=92,file='./data/norm_init_mesh.dat',form='unformatted',status='replace',access='stream')
        write(92) nfm(nfm_size/2,nfm_size/2,1:nfm_size)%x(1)
        write(92) nfm(:,:,:)%u(1)
        write(92) nfm(:,:,:)%u(2)
        write(92) nfm(:,:,:)%u(3)
      close(92)
    end subroutine
    !**********************************************************
    subroutine load_noise
      implicit none
      allocate(xflow_noise(xflow_noise_size))
      open(unit=37,file=noise_file)
        read(37,*) xflow_noise
      close(37)
    end subroutine
    !**********************************************************
    subroutine setup_laizet
      implicit none
      integer :: i
      allocate(laizet_u(128,129,84,3))
      allocate(xp(128))
      allocate(yp(129))
      allocate(zp(84))
      open(unit=37,file='./channel/channel_out.dat',form='unformatted',access='stream')
        read(37) laizet_u(:,:,:,1)
        read(37) laizet_u(:,:,:,2)
        read(37) laizet_u(:,:,:,3)
      close(37)
      laizet_u=laizet_u/2.
      print*, 'max u velocity: ', maxval(laizet_u(:,:,:,1))
      !print*, maxval(laizet_u(:,:,:,2))
      !print*, maxval(laizet_u(:,:,:,3))
      open(unit=37,file='./channel/yp.dat')
        read(37,*) yp
      close(37)
      do i=1,128
        xp(i)=(i-1)*(1./128)*4*pi-2*pi
      end do
      do i=1,84
        zp(i)=(i-1)*(1./84)*4*pi/3-4*pi/6
      end do
      yp=yp/20 
      xp=xp/20
      zp=zp/20
      !check box size! - should be [4\pi,2,4\pi/3]
      write(*,*) 'please check box sizes - it should be:  [4\pi,2,4\pi/3]/20'
    end subroutine
    !**********************************************************
    subroutine laizet_interp(x,u)
      implicit none
      real, intent(IN) :: x(3)
      real, intent(OUT) :: u(3)
      real :: xd,yd,zd
      real,dimension(3) :: c00,c10,c01,c11,c0,c1
      integer :: nx, nnx, ny, nny, nz, nnz
      !find nearest and next nearest points in x,y,z
      !-------------X--------------
      nx=minloc(abs(x(1)-xp),1)
      if (x(1)>xp(nx)) then
        nnx=nx+1
      else
        nx=nx-1
        nnx=nx+1
      end if
      if (nnx==129) nnx=1 !periodicity
      
      !print*, 'x', xp(nx), x(1), xp(nnx)
      !-------------Y--------------
      ny=minloc(abs(x(2)-yp),1)
      if (x(2)>yp(ny)) then
        nny=ny+1
      else
        ny=ny-1
        nny=ny+1
      end if
      !print*, 'y', yp(ny), x(2), yp(nny)
      !-------------Z--------------
      nz=minloc(abs(x(3)-zp),1)
      if (x(3)>zp(nz)) then
        nnz=nz+1
      else
        nz=nz-1
        nnz=nz+1
      end if
      if (nnz==85) nnz=1 !periodicity
      !print*, 'z', zp(nz), x(3), zp(nnz)
      !now we can do the interpolation
      xd=(x(1)-xp(nx))/((1./128)*4*pi)
      yd=(x(2)-yp(ny))/(yp(nny)-yp(ny))
      zd=(x(3)-zp(nz))/((1./84)*4*pi/3)
      !interpolate in x
      c00=laizet_u(nx,ny,nz,:)*(1-xd)+laizet_u(nnx,ny,nz,:)*xd
      c10=laizet_u(nx,nny,nz,:)*(1-xd)+laizet_u(nnx,nny,nz,:)*xd
      c01=laizet_u(nx,ny,nnz,:)*(1-xd)+laizet_u(nnx,ny,nnz,:)*xd
      c11=laizet_u(nx,nny,nnz,:)*(1-xd)+laizet_u(nnx,nny,nnz,:)*xd
      !interpolate in y
      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd
      !interpolate in z
      u=c0*(1-zd)+c1*zd
      !print*, 'interp', u
      !print*, 'nn', laizet_u(nx,ny,nz,:)
    end subroutine
    !**********************************************************
    subroutine randomise_normal_fluid
      implicit none
      select case(normal_velocity)
        case('noisey_xflow')
        xflow_noise_counter=xflow_noise_counter+1 !increment counter
        if (xflow_noise_counter>xflow_noise_size) xflow_noise_counter=1 !reset
      end select
    end subroutine
end module

