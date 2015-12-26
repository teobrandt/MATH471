program homework7
  use mpi
  use xycoord
  implicit none
  integer :: ierr, nprocs, myid
  integer :: status(MPI_STATUS_SIZE)
  
  integer, parameter :: nx = 91
  integer, parameter :: ny = 101

  integer :: nr,ns,nrl,nsl,i,j

  integer :: px, py, ix_off, iy_off ! what is my x-proc index and what is my offset
  integer :: p_left, p_right, px_max, py_max
  integer :: p_down, p_up
  integer :: nxl,nyl !local size
  integer :: remx,remy
  real(kind = 8) :: hr,hs
  real(kind = 8), dimension(:),  allocatable :: r,s
  real(kind = 8), dimension(:,:),allocatable :: t,x,y,u,up,um,ur,us
  real(kind = 8), dimension(:,:),allocatable :: xr,xs,yr,ys,rx,ry,sx,sy,jac,a
  integer :: xint_sum, yint_sum
  CHARACTER(7) :: charit

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

  nr = nx
  ns = ny
  hr = 2.d0/dble(nr-1)
  hs = 2.d0/dble(ns-1)

  ! Label the processes from 1 to px_max
  ! You will have to improve on this and split the
  ! grid in px_max*py_max = nprocs
  ! in such a way that communication is minimized
  call min_domain(nprocs, px_max, py_max)
  call get_id(myid, px, px_max, py, py_max)
  
  ! Split up the grid in the "x-direction"
  nxl = nx/px_max
  remx = nx-nxl*px_max
  nyl = ny/py_max
  remy = ny-nyl*py_max
  if (px .le. remx) then
     nxl = nxl + 1
     ix_off = (px-1)*nxl
  else
     ix_off = (remx)*(nxl+1) + (px-(remx+1))*nxl
  end if
  if (py .le. remy) then
     nyl = nyl + 1
     iy_off = (py-1)*nyl
  else
     iy_off = (remy)*(nyl+1) + (py-(remy+1))*nyl
  end if
  ! write(*,*) 'Crazy arithmetic! : ',px,nxl,ix_off
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_Reduce(nxl,xint_sum,1,&
       MPI_INTEGER,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(nyl,yint_sum,1,&
       MPI_INTEGER,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr)
  if(myid == 0) then
     if (py_max*nx .ne. xint_sum) then
        write(*,*) 'Something is wrong with the number of points in x-direction: ',&
             nx,xint_sum
     end if
     if (px_max*ny .ne. yint_sum) then
        write(*,*) 'Something is wrong with the number of points in y-direction: ',&
             ny,yint_sum
     end if
  end if
  nrl = nxl
  nsl = nyl

  ! who are my neighbours? 
  ! NOTE THAT WE ARE MAPPING TO THE 1D INDEXING HERE!!!!
  p_left  = myid - py_max
  p_right = myid + py_max
  p_up = myid - 1
  p_down = myid + 1
  if (px .eq. px_max) p_right = MPI_PROC_NULL
  if (px .eq. 1) p_left = MPI_PROC_NULL  
  if (py .eq. py_max) p_down = MPI_PROC_NULL
  if (py .eq. 1) p_up = MPI_PROC_NULL

  ! Allocate memory for the various arrays
  allocate(r(nrl),s(nsl),u(0:nrl+1,0:nsl+1),t(0:nrl+1,0:nsl+1),up(0:nrl+1,0:nsl+1),um(0:nrl+1,0:nsl+1),&
       ur(0:nrl+1,0:nsl+1),us(0:nrl+1,0:nsl+1))
  allocate(x(0:nrl+1,0:nsl+1),y(0:nrl+1,0:nsl+1))
  allocate(xr(0:nrl+1,0:nsl+1),yr(0:nrl+1,0:nsl+1),xs(0:nrl+1,0:nsl+1),ys(0:nrl+1,0:nsl+1)&
       ,rx(0:nrl+1,0:nsl+1),ry(0:nrl+1,0:nsl+1),sx(0:nrl+1,0:nsl+1),sy(0:nrl+1,0:nsl+1)&
       ,jac(0:nrl+1,0:nsl+1),a(0:nrl+1,0:nsl+1))
  
  ! Testing communication here
  do i = 1, nrl
     do j = 1, nsl
        t(i, j) = myid + 1
     end do
  end do

  if(myid == 0) then 
     print *, "Before comm tests.."
  end if

  ! send to the top receive from the bottom
  call MPI_Sendrecv(t(0:nrl+1, 1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       t(0:nrl+1, nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(t(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
       t(0:nrl+1, 0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(t(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       t(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(t(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       t(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! End comm tests

  if(myid == 0) then
     print *, "..finished comm tests."
  end if

  do i = 1,nrl
     r(i) = -1.d0 + dble(i-1+ix_off)*hr
  end do
  do i = 1,nsl
     s(i) = -1.d0 + dble(i-1)*hs
  end do

  do j = 1,nsl
     do i = 1,nrl
        x(i,j) = x_coord(r(i),s(j))
        y(i,j) = y_coord(r(i),s(j))
     end do
  end do

  ! At a boundary we extrapolate
  x(1:nrl,0) = 2.d0*x(1:nrl,1)-x(1:nrl,2)
  y(1:nrl,0) = 2.d0*y(1:nrl,1)-y(1:nrl,2)
  x(1:nrl,nsl+1) = 2.d0*x(1:nrl,nsl)-x(1:nrl,nsl-1)
  y(1:nrl,nsl+1) = 2.d0*y(1:nrl,nsl)-y(1:nrl,nsl-1)

  ! Only proc 1 and px_max have boundaries
  if(px .eq. 1) then
     x(0,:) = 2.d0*x(1,:)-x(2,:)
     y(0,:) = 2.d0*y(1,:)-y(2,:)
  end if
  if(px .eq. px_max) then
     x(nrl+1,:) = 2.d0*x(nrl,:)-x(nrl-1,:)
     y(nrl+1,:) = 2.d0*y(nrl,:)-y(nrl-1,:)
  end if

  ! send to the top receive from the bottom
  call MPI_Sendrecv(x(0:nrl+1, 1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       x(0:nrl+1, nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(0:nrl+1, 1),nrl+2,MPI_DOUBLE_PRECISION,p_up,124,&
       y(0:nrl+1, nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,124,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(x(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,125,&
       x(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,125,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,&
       y(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(x(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,127,&
       x(0:nrl+1, 0),nrl+2,MPI_DOUBLE_PRECISION,p_up,127,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,128,&
       y(0:nrl+1, 0),nrl+2,MPI_DOUBLE_PRECISION,p_up,128,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(x(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,129,&
       x(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,129,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,130,&
       y(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,130,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(myid == 0) then
     print *, "Through #1"
  end if

  ! Now compute the metric
  xr(1:nrl,:) = (x(2:nrl+1,:)-x(0:nrl-1,:))/(2.d0*hr)
  yr(1:nrl,:) = (y(2:nrl+1,:)-y(0:nrl-1,:))/(2.d0*hr)
  xs(:,1:nsl) = (x(:,2:nsl+1)-x(:,0:nsl-1))/(2.d0*hs)
  ys(:,1:nsl) = (y(:,2:nsl+1)-y(:,0:nsl-1))/(2.d0*hs)

  ! At a boundary we extrapolate
  xr(1:nrl,0) = 2.d0*xr(1:nrl,1)-xr(1:nrl,2)
  xs(1:nrl,0) = 2.d0*xs(1:nrl,1)-xs(1:nrl,2)
  yr(1:nrl,0) = 2.d0*yr(1:nrl,1)-yr(1:nrl,2)
  ys(1:nrl,0) = 2.d0*ys(1:nrl,1)-ys(1:nrl,2)
  xr(1:nrl,nsl+1) = 2.d0*xr(1:nrl,nsl)-xr(1:nrl,nsl-1)
  xs(1:nrl,nsl+1) = 2.d0*xs(1:nrl,nsl)-xs(1:nrl,nsl-1)
  yr(1:nrl,nsl+1) = 2.d0*yr(1:nrl,nsl)-yr(1:nrl,nsl-1)
  ys(1:nrl,nsl+1) = 2.d0*ys(1:nrl,nsl)-ys(1:nrl,nsl-1)

  ! Only proc 1 and px_max have boundaries
  if(px .eq. 1) then
     xr(0,:) = 2.d0*xr(1,:)-xr(2,:)
     xs(0,:) = 2.d0*xs(1,:)-xs(2,:)
     yr(0,:) = 2.d0*yr(1,:)-yr(2,:)
     ys(0,:) = 2.d0*ys(1,:)-ys(2,:)
  end if
  if(px .eq. px_max) then
     xr(nrl+1,:) = 2.d0*xr(nrl,:)-xr(nrl-1,:)
     xs(nrl+1,:) = 2.d0*xs(nrl,:)-xs(nrl-1,:)
     yr(nrl+1,:) = 2.d0*yr(nrl,:)-yr(nrl-1,:)
     ys(nrl+1,:) = 2.d0*ys(nrl,:)-ys(nrl-1,:)
  end if

  ! send to the top recieve from the bottom
  call MPI_Sendrecv(xr(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,223,&
       xr(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,223,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,224,&
       yr(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,224,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(xr(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,225,&
       xr(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,225,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,226,&
       yr(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,226,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(xr(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,227,&
       xr(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,227,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,228,&
       yr(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,228,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(xr(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,229,&
       xr(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,229,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,230,&
       yr(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,230,MPI_COMM_WORLD,status,ierr)

  ! send to the top recieve from the bottom
  call MPI_Sendrecv(xs(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,231,&
       xs(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,231,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,232,&
       ys(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,232,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(xs(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,233,&
       xs(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,233,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,234,&
       ys(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,234,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(xs(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,235,&
       xs(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,235,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,236,&
       ys(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,236,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(xs(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,237,&
       xs(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,237,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,238,&
       ys(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,238,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(myid == 0) then
     print *, "Through #2"
  end if

  jac = xr*ys-xs*yr
  
  rx =  ys/jac
  ry = -xs/jac
  sx = -yr/jac
  sy =  xr/jac

  ! send to the top recieve from the bottom
  call MPI_Sendrecv(u(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       u(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(u(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       u(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(u(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
       u(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(u(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       u(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(myid == 0) then
     print *, "Through #3"
  end if
  

  u =  x + y 
  ! Now compute the derivative in the interior
  ur(1:nrl,1:nsl) = (u(2:nrl+1,1:nsl)-u(0:nrl-1,1:nsl))/(2.d0*hr)
  us(1:nrl,1:nsl) = (u(1:nrl,2:nsl+1)-u(1:nrl,0:nsl-1))/(2.d0*hs)
  up = rx*ur+sx*us
  WRITE(charit,"(I7.7)") myid
  call  printdble2d(t(0:nrl+1,0:nsl+1),nrl+2,nsl+2,'t'//charit//'.txt')
  call  printdble2d(x(1:nrl,1:nsl),nrl,nsl,'x'//charit//'.txt')
  call  printdble2d(y(1:nrl,1:nsl),nrl,nsl,'y'//charit//'.txt')
  call  printdble2d(up(1:nrl,1:nsl),nrl,nsl,'up'//charit//'.txt')
  
  
  ! Test Boundary Conditions
  ! Test Forcing
  ! Test Computation of right hand side 
  !do i = ..
  !   do j = 
  !      call mms(u_MMS,x(i,j),y(i,j),t,m)
  !        utt(i,j) = u_mms(0,0,2)
  !        uxx(i,j) = u_mms(2,0,0)
  !        uexact(i,j)= u_mms(0,0,0)


! compute error
!write(*,*) maxval(abs(u-uexact)) 
 ! Finally do some..
  ! Time-stepping
  ! Output
  ! 

  ! Deallocate
  call mpi_finalize(ierr)
end program homework7

subroutine min_domain(n_procs, px_max, py_max)
  implicit none
  real (kind=8) :: div
  integer, intent(in) :: n_procs
  integer, intent(out) :: px_max, py_max
  integer :: i

  do i = 1, ceiling(sqrt(real(n_procs)))
     div = real(n_procs) / real(i)
     if (div .eq. floor(div)) then
        px_max = i
        py_max = n_procs / i
     end if
  end do
end subroutine min_domain

subroutine get_id(my_id, px, px_max, py, py_max)
  implicit none
  integer, intent(in) :: my_id, px_max, py_max
  integer, intent(out) :: px, py

  px = (my_id / px_max) + 1
  py = mod(my_id, py_max) + 1

  write (*,*) my_id, px, py
end subroutine get_id

subroutine printdble2d(u,nx,ny,str)
  implicit none
  integer, intent(in) :: nx,ny
  real(kind = 8), intent(in) :: u(nx,ny)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=1,ny,1
     do i=1,nx,1
        write(2,fmt='(E24.16)',advance='no') u(i,j)
     end do
     write(2,'()')
  end do
  close(2)
end subroutine printdble2d
