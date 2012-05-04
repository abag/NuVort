!> diagnostic routines purely for vortex filament
module diagnostic
  use cdata
  use general
  contains
  !>dummy routine to call all diagnostic routines
  subroutine calculate_diagnostics()
    implicit none
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
    end do
    !$omp end parallel do
    !compute the average/min/max of this array
    kappa_bar=sum(curvi)/count(mask=f(:)%infront>0)
    kappa_max=maxval(curvi)
    kappa_min=minval(curvi,mask=curvi>0)
    deallocate(curvi) !deallocate helper array
  end subroutine
end module
