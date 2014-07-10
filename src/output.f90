!>The main output routines from the code, whilst specifc modules may contain there own
!>output, the main routines for printing structures should be in here.
module output
  use cdata
  use diagnostic
  contains
  !**********************************************************************
  !> print various dimensions/flags to file for matlab to read in for plotting
  subroutine print_dims()
    implicit none
    open(unit=77,file='./data/dims.log',status='replace')
      write(77,*) delta, '%resolution - \delta' 
      write(77,*) box_size(1), '%x box size'
      write(77,*) box_size(2), '%y box size'
      write(77,*) box_size(3), '%z box size'
    close(77)
  end subroutine
  !**********************************************************************
  !>print information to screen/file
  subroutine print_info()
    implicit none
    call printf(itime/shots) !output.mod
    open(unit=78,file='data/ts.log',position='append')
    if (itime==shots) then
      write(*,'(a)') '--var--------t--------pcount--------recon-----wall_recon---avg_d-----&
                length--------maxu---------maxdu-------curv------removed'
      write(78,*) '%--var--------t-------pcount--------recon-----wall_recon----avg_d-----&
                   length--------maxu---------maxdu------curv------removed-----tree_eval'
    end if
    write(*,'(i6.4,f13.7,i10.1,i13.1,i13.1,f7.4,f13.6,f13.5,f13.5,f10.2,i13.1)') &
    itime/shots,t,count(mask=f(:)%infront>0),recon_count,wall_recon_count,avg_sep/delta,&
    total_length,maxu,maxdu,kappa_bar,&
    remove_count
    write(78,'(i6.4,f13.7,i10.1,i13.1,i13.1,f7.4,f13.6,f13.5,f13.5,f10.2,i13.1,i13.1)') &
itime/shots,t,count(mask=f(:)%infront>0),recon_count,wall_recon_count,avg_sep/delta,&
total_length,maxu,maxdu,kappa_bar,&
remove_count, tree_eval
    close(78)
    open(unit=79,file='data/curvature.log',position='append')
      write(79,*) kappa_bar, kappa_min, kappa_max
    close(79)
    open(unit=79,file='data/super_vs_friction.log',position='append')
      write(79,*) t, mean_u_sup, mean_u_mf
  end subroutine
  !**********************************************************************
  !>print the f (filament) array as (un)formatted data for use with gnuplot/matlab
  subroutine printf(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (filenumber==10000) call warning_message('output.mod','run out of filenumbers to print var to')
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/var",filenumber,".log"
    open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
      write(98) t
      write(98) pcount
      write(98) f(:)%x(1)
      write(98) f(:)%x(2)
      write(98) f(:)%x(3)
      write(98) f(:)%infront
      write(98) f(:)%u(1)
      write(98) f(:)%u(2)
      write(98) f(:)%u(3)
      write(98) f(:)%u_mf(1)
      write(98) f(:)%u_mf(2)
      write(98) f(:)%u_mf(3)
      write(98) f(:)%curv
      write(98) f(:)%stretch
    close(98)
  end subroutine
  !**********************************************************************
  !>print the initial f (filament) array as formatted data for use with gnuplot/matlab
  subroutine initial_printf
    implicit none
    integer :: i !for looping
    open(unit=98,file='./data/var_initial.log',status='replace')
      do i=1, pcount
        write(98,*) f(i)%x(1:3), f(i)%infront
      end do
    close(98)
  end subroutine
  !**********************************************************************
  !>store everything needed to restart the code
  !>\todo a few diagnostics are not dumped and resume from 0 please fix
  subroutine data_dump
    implicit none
    open(unit=53,file="./data/var.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) recon_count
      write(53) itime
      write(53) t
      write(53) f
    close(53)
  end subroutine
  !**************************************************
  subroutine banner_print()
    implicit none
    integer :: today(3), now(3)
    character(len=30) :: user_name, host_name
    call getenv('USER', user_name)
    call hostnm(host_name)
    call idate(today) ! today(1)=day, (2)=month, (3)=year
    call itime(now)   ! now(1)=hour, (2)=minute, (3)=second
    write(*,*) "                                        " 
    write(*,*) "        _   _     __     __         _   "  
    write(*,*) "       | \ | |_   \ \   / /__  _ __| |_ "
    write(*,*) "       |  \| | | | \ \ / / _ \| '__| __|"
    write(*,*) "       | |\  | |_| |\ V / (_) | |  | |_ "
    write(*,*) "       |_| \_|\__,_| \_/ \___/|_|   \__|"
    write(*,*) "                                        "
    write(*,*) "                                        " 
    write(*,*) 'user info: ', trim(user_name), '@', trim(host_name)
    write ( *, 10 )  today(1), today(2), today(3), now
    10    format ( ' date ', i2.2, '/', i2.2, '/', i4.4, '; time ', &
                  i2.2, ':', i2.2, ':', i2.2 )
    write(*,*) "-------------------------------------" 
  end subroutine
end module
