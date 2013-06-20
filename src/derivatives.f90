!>routines to calculate spatial derivatives
module derivatives
  use cdata
  contains
  !> get the first derivative at position i using a second order adaptive mesh
  !>finite difference scheme:
  !!\f[
  !!\frac{d \mathbf{s}_i}{d \xi}=\frac{\ell_{i-1}\mathbf{s}_{i+1}+(\ell_{i+1}-\ell_{i-1})\mathbf{s}_i+\ell_{i+1}\mathbf{s}_{i-1}}
  !!  {2\ell_{i+1}\ell_{i-1}}+{\cal O}(\ell^2)
  !! \f]
  subroutine get_deriv_1(i,s_dot)
    implicit none
    integer, intent(IN) :: i
    real :: s_dot(3), disti, distb
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      s_dot(1:3)=distb*f(i)%ghosti(1:3)+ &
                (disti-distb)*f(i)%x(1:3)- &
                 disti*f(i)%ghostb(1:3)
      s_dot=s_dot/(2.*distb*disti)      
  end subroutine
  !*********************************************************************
  !>get the second derivative at position i using a second order adaptive 
  !>mesh finite difference scheme:
  !!\f[
  !!\frac{d^2 \mathbf{s}_i}{d \xi^2}=\frac{2\mathbf{s}_{i+1}}{\ell_{i+1}(\ell_{i+1}+\ell_{i-  1})}-
  !!\frac{2\mathbf{s}_i}{\ell_{i+1}\ell_{i-1}}+\frac{2\mathbf{s}_{i-1}}{\ell_{i-1}(\ell_{i+  1}+\ell_{i-1})}+{\cal O}(\ell^2)
  !! \f]
  subroutine get_deriv_2(i,s_ddot)
    implicit none
    integer, intent(IN) :: i
    real :: s_ddot(3), disti, distb
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      s_ddot(1:3)=2.*(f(i)%ghosti(1:3)/(disti*(disti+distb))- &
                      f(i)%x(1:3)/(disti*distb)+ &
                      f(i)%ghostb(1:3)/(distb*(disti+distb)))
  end subroutine
  !*********************************************************************
  !>calculate the distance between points \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  !>just pythagoras
  real function dist_gen(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen=sqrt((a(1)-b(1))**2+&
                  (a(2)-b(2))**2+&
                  (a(3)-b(3))**2)
  end function
end module

