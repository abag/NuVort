!> diagnostic routines purely for vortex filament
module diagnostic
  use cdata
  use general
  use normal_fluid
  contains
  !>dummy routine to call all diagnostic routines
  subroutine calculate_diagnostics()
    implicit none
    !call friction_force
    if (mod(itime, shots)==0) then
      call velocity_info !diagnostics.mod
      call curv_info !diagnostics.mod
    end if
  end subroutine
  !*************************************************
  !>get the maximum velocity and change of velocity
  !>on the filament
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    allocate(uinfo(pcount,2))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(:)%u2(1))**2+&
                    (f(:)%u1(2)-f(:)%u2(2))**2+&
                    (f(:)%u1(3)-f(:)%u2(3))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
    deallocate(uinfo)
    mean_u_sup=sum(sqrt(f(:)%u_sup(1)**2+f(:)%u_sup(2)**2+f(:)%u_sup(3)**2))/count(mask=f(:)%infront>0)
    mean_u_mf=sum(sqrt(f(:)%u_mf(1)**2+f(:)%u_mf(2)**2+f(:)%u_mf(3)**2))/count(mask=f(:)%infront>0)
  end subroutine 
  !*************************************************
  !>caculate the mean, min, max curvature of the filament
  subroutine curv_info()
    implicit none
    real, allocatable :: curvi(:)
    integer :: i, j
    allocate(curvi(pcount)) !allocate this array pcount size
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) then
        curvi(i)=0. !check for 'empty' particles
      else
        curvi(i)=curvature(i) !general.mod
      end if
      f(i)%curv=curvi(i) !store this for printing
    end do
    !$omp end parallel do
    !compute the average/min/max of this array
    kappa_bar=sum(curvi)/count(mask=f(:)%infront>0)
    kappa_max=maxval(curvi)
    kappa_min=minval(curvi,mask=curvi>0)
    deallocate(curvi) !deallocate helper array
  end subroutine
  !*************************************************
  subroutine friction_force
    integer :: i
    real :: f_dot(3)  
    real :: u_norm(3)
    real :: disti
    real :: force(3), force2(3), force3(3), force4(3)
    real :: max_force(pcount)
    real, parameter :: D=2.79E-5
    real, parameter :: Dt=-4.14E-6
    force=0. ; force2=0. ; force3=0. ; force4=0.
    open(unit=76,file='./test_friction.log')
     do i=1, pcount
      if (f(i)%infront==0) cycle

        call get_deriv_1(i,f_dot)
        call get_normal_velocity(f(i)%x,u_norm)
       disti=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
       force=force+cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf)))*disti
       force2=force2+cross_product(f_dot,cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf))))*disti
       force3=force3+abs(D*cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf)))*disti+&
              Dt*cross_product(f_dot,cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf))))*disti)
       force4=force4+abs(D*cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf)))*disti)+&
              abs(Dt*cross_product(f_dot,cross_product(f_dot,(u_norm-(f(i)%u-f(i)%u_mf))))*disti)
       max_force(i)=sqrt(dot_product(force3,force3))
    end do
   print*, force
   print*, force2
   force=force*D
   force2=force2*Dt
   print*, sqrt(dot_product(force+force2,force+force2))/(0.1*0.1*0.2)
   !print*, sqrt(dot_product(abs(force)+abs(force2),abs(force)+abs(force2)))/(0.1*0.1*0.2)
   print*, sqrt(dot_product(force3,force3))/(0.1*0.1*0.2)
   print*, maxval(force3)
   !print*, sqrt(dot_product(force4,force4))/(0.1*0.1*0.2)
   close(76)
   stop
  end subroutine
end module
