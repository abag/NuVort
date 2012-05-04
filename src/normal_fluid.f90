!>normal fluid component in the equation of motion, if the filament acts as a vortex
!>velocity enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$
!>if the filament is a material line then the normaly velocity is used directly.
!>for further information on the normal fluid options see \ref NF
module normal_fluid
  use cdata
  use general
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
        case('zero')
        case default
          call fatal_error('normal_fluid.mod:get_normal_fluid', &
          'correct parameter for normal_veloctity not set')
      end select
      write(*,'(a,f6.4,a,f6.4)') ' mutual friction coefficients alpha=',alpha(1),' alpha`=',alpha(2) 
    end subroutine
    !************************************************************
    subroutine get_normal_velocity(x,u)
      implicit none
      real, intent(IN) :: x(3) !position of the particle
      real, intent(OUT) :: u(3) !velocity at x
      real :: r, phi, theta !used to convert to polar coords
      real :: u_r, u_theta !used to convert to polar coords
      integer :: peri, perj, perk !used to loop in periodic cases
      u=0. ! a safety check , 0 the field before we begin!
      select case(normal_velocity)
        case('zero')
          u=0. ! no flow
        case('xflow')
          u=0. ; u(1)=norm_vel_xflow !flow in the x direction
      end select
    end subroutine    
end module

