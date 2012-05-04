!>general routines e.g. vector products, distance evaluation
!>this is used by all other modules (excluding cdata/derivatives)
module general
  use cdata
  use derivatives
  contains
  !*********************************************************************
  !>empty all a points information - used if a point is removed
  !>(e.g. contraction of filament)
  subroutine clear_particle(i)
    implicit none
    integer, intent(IN) :: i
    f(i)%x=0.
    f(i)%u=0. ; f(i)%u1=0. ; f(i)%u2=0.
    f(i)%ghosti=0. ; f(i)%ghostb=0.
    f(i)%infront=0 ; f(i)%behind=0
    f(i)%closest=0 ; f(i)%closestd=0.
    f(i)%pinnedi=.false. ; f(i)%pinnedb=.false.
    f(i)%wpinned=0
  end subroutine
  !*********************************************************************
  !>calculate the distance between points in the f vector
  !!The distance between \f$(x_1,y_1,z_1)\f$ and \f$(x_2,y_2,z_2)\f$ is 
  !!\f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$.
  real function distf(i,j)
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distf=sqrt((f(i)%x(1)-f(j)%x(1))**2+&
               (f(i)%x(2)-f(j)%x(2))**2+&
               (f(i)%x(3)-f(j)%x(3))**2)
  end function
  !*********************************************************************
  !>calculate the squared distance between particles in the f vector
  !!\f${(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$.
  real function distfsq(i,j)
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distfsq=(f(i)%x(1)-f(j)%x(1))**2+&
            (f(i)%x(2)-f(j)%x(2))**2+&
            (f(i)%x(3)-f(j)%x(3))**2
  end function
  !*********************************************************************
  !>calculate the curvature at the particle i: \f$|\mathbf{s}''|\f$
  real function curvature(i)
    use Cdata
    implicit none
    integer, intent(IN) :: i
    real :: fddot(3)
    call get_deriv_2(i,fddot)
    curvature=sqrt(dot_product(fddot,fddot))
  end function
  !*********************************************************************
  !>calculate the squared distance between points \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  real function dist_gen_sq(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen_sq=(a(1)-b(1))**2+&
                (a(2)-b(2))**2+&
                (a(3)-b(3))**2
  end function
  !*********************************************************************
  !>calculate the angle between vectors \f$\mathbf{a}\f$ and  \f$\mathbf{b}\f$
  !!\f[ \theta=\cos^{-1} \frac{ \mathbf{a} \cdot \mathbf{b}}{ab} \f]
  real function vector_angle(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    vector_angle=acos(dot_product(a,b)/ &
                 (sqrt(dot_product(a,a))*sqrt(dot_product(b,b))))
  end function
  !*********************************************************************
  !>call finite difference to calculate the tangent vector at point i
  !!and then normalise, used to check angle between vortices when reconnecting.
  function tangentf(i)
    use Cdata
    implicit none
    real, dimension(3) :: tangentf
    integer, intent(IN) :: i
    real :: dist
    !get distance between particle and ghost particle (infront)
    dist=dist_gen(f(i)%x,f(i)%ghosti)
    !now determine the vector (low order finite diff.)
    tangentf(:)=(f(i)%ghosti(:)-f(i)%x(:))/dist
  end function
  !*********************************************************************
  !>calculate the normalised tangent vector at point i low order fininte diff.
  function norm_tanf(i)    
    use Cdata
    implicit none
    real, dimension(3) :: norm_tanf
    integer, intent(IN) :: i
    real :: length !length of vector
    !now determine the vector (finite diff.)
    call get_deriv_1(i,norm_tanf)
    !calculate length of vector
    length=sqrt(norm_tanf(1)**2+norm_tanf(2)**2+norm_tanf(3)**2)
    norm_tanf(:)=norm_tanf(:)/length !normalise
  end function
  !*********************************************************************
  !>calculate the normalised normal vector at point i:
  !>\f$\hat{\mathbf{s}''}\f$
  function normalf(i)
    use Cdata
    implicit none
    real, dimension(3) :: normalf
    integer, intent(IN) :: i
    real :: length !length of vector
    !get the second derivative
    call get_deriv_2(i,normalf)
    !calculate length of vector
    length=sqrt(normalf(1)**2+normalf(2)**2+normalf(3)**2)
    normalf(:)=normalf(:)/length !normalise
  end function
  !*********************************************************************
  !>calculate the normalised binormal vector at point i:
  !>\f$\hat{\mathbf{s}'} \times \hat{\mathbf{s}''}\f$
  function binormalf(i)
    use Cdata
    implicit none
    real, dimension(3) :: s_dot, s_ddot, binormalf
    integer, intent(IN) :: i
    real :: length !length of vector
    !get the first/second derivative
    call get_deriv_1(i,s_dot)
    call get_deriv_2(i,s_ddot)
    binormalf=cross_product(s_dot,s_ddot)
    !calculate length of vector
    length=sqrt(binormalf(1)**2+binormalf(2)**2+binormalf(3)**2)
    binormalf(:)=binormalf(:)/length !normalise
  end function
  !*********************************************************************
  !>calculate the cross product of \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  function cross_product(a,b)
    implicit none
    real, dimension(3) :: cross_product
    real, dimension(3), intent(IN) :: a, b
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
  end function
  !*********************************************************************
  !>calculate the L2-norm of input a
  !!\f[ |mathbf{x}|=\sqrt{x_1^2+x_2^2+x_3^2} \f]
  real function vector_norm(a)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a
    vector_norm=(sqrt(dot_product(a,a)))
  end function
  !*********************************************************************
  !> a routine to find any NANs in the positions of all the allocated main
  !> arrays
  subroutine NAN_finder
    implicit none
    integer :: i
    if (allocated(f)) then
      if (any(isnan(f%x(1))).or.any(isnan(f%x(2))).or.any(isnan(f%x(3)))) then
        do i=1, pcount
          if ((isnan(f(i)%x(1))).or. &
              (isnan(f(i)%x(2))).or. &
              (isnan(f(i)%x(3)))) then
            write(*,*) 'I have found a NAN in the f%x array'
            write(*,*) 'location, i=',i
            write(*,*) 'f(i)%x=',f(i)%x
            write(*,*) 'f(i)', f(i)
          end if
        end do
        call fatal_error('run.x','there is a NAN in the (filament) f%x array')
      end if
    end if
  end subroutine
  !*********************************************************************
  !> a routine purely used for code testing finds empty particles 
  subroutine zero_finder(location)
    implicit none
    integer :: i, zcount=0
    character(len=*) :: location
    do i=1, pcount
      if (f(i)%infront==0) zcount=zcount+1
    end do
    write(*,*) 'zero finder called at', trim(location)
    write(*,*) 'zero count= ', zcount
    write(*,*) 'ending run...' ; stop 
  end subroutine
end module

