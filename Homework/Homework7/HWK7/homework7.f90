program homework7
  use mpi
  use xycoord
  implicit none

  integer :: ierr, nprocs, myid
  integer :: status(MPI_STATUS_SIZE)
  integer, parameter :: nx = 160
  integer, parameter :: ny = 160
  real(kind = 8), parameter :: t_end = 10.d0
  real(kind = 8), parameter :: c = 1.d0 !This is the speed of the fluid (ie.speed of sound)
  integer :: nr,ns,nrl,nsl,i,j,it,nt
  integer :: px,py,ix_off,iy_off ! what is my x-proc index and what is my offset
  integer :: p_left,p_right,px_max
  integer :: p_down,p_up,py_max
  integer :: nxl,nyl !local size
  integer :: remx,remy
  real(kind = 8) :: hr,hs,dt,a,t
  real(kind = 8), dimension(:),  allocatable :: r,s
  real(kind = 8), dimension(:,:),allocatable :: x,y,u,f,up,um,ur,us
  real(kind = 8), dimension(:,:),allocatable :: xr,xs,yr,ys,rx,ry,sx,sy,jac,u_tt,Ju_tt
  real(kind = 8), dimension(:,:),allocatable :: u_ex,u_ttex,error_tave
  integer :: xint_sum, yint_sum
  CHARACTER(7) :: charit!,charstep
  real(kind = 8) :: mpit, totit
  real(kind = 8), dimension(:), allocatable :: locerr,toterr
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  
  !===================================================================================
  !Setting up the grid parameters and communication paths
  !===================================================================================
  mpit = MPI_Wtime()

  nr = nx!The number of points in the r direction is the same as in x
  ns = ny!"..."s"..."y"..."
  hr = 2.d0/dble(nr-1)!The step size: hr=2/(nr-1) ex. nr=50, hr~0.041
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
  
  !===================================================================================
  !Generate grid in r,s coordinate system from x,y
  !===================================================================================
  ! Allocate memory for the various arrays
  allocate(r(nrl),s(nsl),u(0:nrl+1,0:nsl+1),up(0:nrl+1,0:nsl+1),um(0:nrl+1,0:nsl+1),&
       ur(0:nrl+1,0:nsl+1),us(0:nrl+1,0:nsl+1),error_tave(0:nrl+1,0:nsl+1))
  allocate(x(0:nrl+1,0:nsl+1),y(0:nrl+1,0:nsl+1))
  allocate(xr(0:nrl+1,0:nsl+1),yr(0:nrl+1,0:nsl+1),xs(0:nrl+1,0:nsl+1),ys(0:nrl+1,0:nsl+1)&
       ,rx(0:nrl+1,0:nsl+1),ry(0:nrl+1,0:nsl+1),sx(0:nrl+1,0:nsl+1),sy(0:nrl+1,0:nsl+1)&
       ,jac(0:nrl+1,0:nsl+1),f(0:nrl+1,0:nsl+1))

  !What are these offsets for?? These offsets account for the ghostpoints since the calc u_tt occurs in
  !the interior of the grid and there is a boarder of "ghost points"
  !r goes from -1 to 1, s goes fro -1 to 1
  do i = 1,nrl
     r(i) = -1.d0 + dble(i-1+ix_off)*hr
  end do
  do i = 1,nsl
     s(i) = -1.d0 + dble(i-1+iy_off)*hs
  end do
  !write(*,*) 'r',r!,'s',s

  !Here we call the x and y curvalinear coordinates
  do j = 1,nsl
     do i = 1,nrl
        x(i,j) = x_coord(r(i),s(j))
        y(i,j) = y_coord(r(i),s(j))
     end do
  end do
  !write(*,*) 'x',x(1:nrl,1)!,'y',y

  ! At a boundary we extrapolate
  x(1:nrl,0) = 2.d0*x(1:nrl,1)-x(1:nrl,2)
  y(1:nrl,0) = 2.d0*y(1:nrl,1)-y(1:nrl,2)
  x(1:nrl,nsl+1) = 2.d0*x(1:nrl,nsl)-x(1:nrl,nsl-1)
  y(1:nrl,nsl+1) = 2.d0*y(1:nrl,nsl)-y(1:nrl,nsl-1)
  !write(*,*) 'x',x(1,1:nrl)

  ! Only proc 1 and px_max have boundaries
  if(px .eq. 1) then
     x(0,:) = 2.d0*x(1,:)-x(2,:)
     y(0,:) = 2.d0*y(1,:)-y(2,:)
  end if
  if(px .eq. px_max) then
     x(nrl+1,:) = 2.d0*x(nrl,:)-x(nrl-1,:)
     y(nrl+1,:) = 2.d0*y(nrl,:)-y(nrl-1,:)
  end if

  !===================================================================================
  !Communication using MPI
  !===================================================================================
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

  !===================================================================================
  !Compute the derivatives of x,y w.r.t r,s using finite difference formulas
  !===================================================================================
  ! Now compute the metric
  !write(*,*) x(2:nrl+1,:)

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
  call MPI_Sendrecv(ys(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,236,&
       ys(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,236,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(xs(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,237,&
       xs(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,237,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,238,&
       ys(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,238,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !===================================================================================
  !Pre-Space Discretization
  !===================================================================================
  do i = 0,nrl+1
     do j = 0,nsl+1
        jac(i,j) = xr(i,j)*ys(i,j)-xs(i,j)*yr(i,j)
     end do
  end do
  ! send to the top recieve from the bottom
  call MPI_Sendrecv(jac(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       jac(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(jac(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       jac(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(jac(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
       jac(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(jac(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       jac(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !write(*,*) 'jac',jac
  !stop
  do i = 0,nrl+1
     do j = 0,nsl+1
        rx(i,j) = ys(i,j)/jac(i,j)
        ry(i,j) = -xs(i,j)/jac(i,j)
        sx(i,j) = -yr(i,j)/jac(i,j)
        sy(i,j) = xr(i,j)/jac(i,j)
     end do
  end do

    ! send to the top recieve from the bottom
  call MPI_Sendrecv(rx(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,223,&
       rx(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,223,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ry(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,224,&
       ry(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,224,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(rx(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,225,&
       rx(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,225,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ry(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,226,&
       ry(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,226,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(rx(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,227,&
       rx(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,227,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ry(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,228,&
       ry(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,228,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(rx(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,229,&
       rx(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,229,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ry(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,230,&
       ry(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,230,MPI_COMM_WORLD,status,ierr)

  ! send to the top recieve from the bottom
  call MPI_Sendrecv(sx(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,231,&
       sx(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,231,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(sy(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,232,&
       sy(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,232,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(sx(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,233,&
       sx(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,233,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(sy(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,234,&
       sy(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,234,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(sx(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,235,&
       sx(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,235,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(sy(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,236,&
       sy(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,236,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(sx(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,237,&
       sx(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,237,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(sy(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,238,&
       sy(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,238,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !===================================================================================
  !u(x,y,t) = ...
  !===================================================================================
  do i = 0,nrl+1
     do j = 0,nsl+1
        u(i,j) = sin(-x(i,j))*sin(y(i,j)) !For MMStrig
        !u(i,j) = (x(i,j)+x(i,j)*x(i,j))*&
             !(y(i,j)+y(i,j)*y(i,j))*&
             !(0.d0) !For MMStrig (t = 0)
     end do
  end do
  !===================================================================================
  !Boundary Values and Initial Conditions - none for MMS check
  !===================================================================================
  !u(0:nrl+1,0) = 0
  !u(0,0:nsl+1) = 0
  !u(nrl+1,0:nsl+1) = 0
  !u(0:nrl+1,nsl+1) = 0
  
  !write(*,*) u(1,1:nsl)
  !stop
  
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

  ! Now compute the derivative in the interior
  ur(1:nrl,1:nsl) = (u(2:nrl+1,1:nsl)-u(0:nrl-1,1:nsl))/(2.d0*hr)
  us(1:nrl,1:nsl) = (u(1:nrl,2:nsl+1)-u(1:nrl,0:nsl-1))/(2.d0*hs)
  !do i = 1,nrl
  !   ur(i,nsl) = 0.d0
  !   ur(i,0) = 0.d0
  !end do
  ! Now comute the derivative at the boundaries
  ur(1:nrl,0) = 2.d0*ur(1:nrl,1)-ur(1:nrl,2)
  us(1:nrl,0) = 2.d0*us(1:nrl,1)-us(1:nrl,2)
  ur(1:nrl,nsl+1) = 2.d0*ur(1:nrl,nsl)-ur(1:nrl,nsl-1)
  us(1:nrl,nsl+1) = 2.d0*us(1:nrl,nsl)-us(1:nrl,nsl-1)
  ! send to the top recieve from the bottom
  call MPI_Sendrecv(ur(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,231,&
       ur(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,231,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(us(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,232,&
       us(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,232,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(ur(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,233,&
       ur(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,233,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(us(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,234,&
       us(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,234,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(ur(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,235,&
       ur(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,235,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(us(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,236,&
       us(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,236,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(ur(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,237,&
       ur(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,237,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(us(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,238,&
       us(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,238,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i = 0,nrl+1
     do j = 0,nsl+1
        up(i,j) = rx(i,j)*ur(i,j)+sx(i,j)*us(i,j)
     end do
  end do
  ! send to the top recieve from the bottom
  call MPI_Sendrecv(up(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       up(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(up(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       up(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(up(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
       up(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(up(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       up(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  WRITE(charit,"(I7.7)") myid
  call  printdble2d(x(1:nrl,1:nsl),nrl,nsl,'x'//charit//'.txt')
  call  printdble2d(y(1:nrl,1:nsl),nrl,nsl,'y'//charit//'.txt')
  call  printdble2d(up(1:nrl,1:nsl),nrl,nsl,'up'//charit//'.txt')

  !===================================================================================
  !Test Boundary Conditions
  !===================================================================================

  !===================================================================================
  !Test Forcing
  !===================================================================================
  do i = 0,nrl+1
     do j = 0,nsl+1
        f(i,j) = 2.d0*c*c*u(i,j) !For MMStrig
        !f(i,j) = 0 !For MMSpoly
     end do
  end do
  !write(*,*) 'u',u(1,1:nsl)
  !write(*,*) 'f',f(1,1:nsl)
  !stop
  call MPI_Sendrecv(f(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
       f(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(f(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       f(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the bottom recieve from the top
  call MPI_Sendrecv(f(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
       f(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(f(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       f(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  !===================================================================================
  !Test Computation of right hand side
  !===================================================================================
  allocate(u_tt(0:nrl-1,0:nsl-1),Ju_tt(0:nrl+1,0:nsl+1))
  a = 1.d0
  call compute_rhs(Ju_tt,u,f,rx,ry,sx,sy,c,hr,hs,nrl,nsl,jac)
  do i = 1,nrl
     do j = 1,nsl
        u_tt(i-1,j-1) = Ju_tt(i,j)/jac(i,j)
     end do
  end do
  
  write(*,*) 'Error:', maxval(abs(u_tt))!,'utt is:', u_tt(2:nrl-1,1)
  call MPI_Reduce(mpit,totit,1,&
       MPI_REAL8,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr)
  if(myid == 0) then
     write(*,*) 'time', totit / real(px_max * py_max)
  end if
  !stop
  !stop
  !write(*,*) 'What is the error', max(hr,hs),maxval(abs(u_tt)2.d0*sin(x(1:nrl,1:nsl))*sin(y(1:nrl,1:nsl)))
  !stop
  !===================================================================================
  !Time Discretization
  !===================================================================================
  do i = 0,nrl+1
     do j = 0,nsl+1
        um(i,j) = u(i,j)
     end do
  end do

  !dt = 0.1d0*minval(sqrt(abs(jac)))*min(hr,hs)
  dt = 0.0075d0
  nt = floor(t_end/dt)+1
  dt = t_end/dble(nt)
  allocate(u_ex(1:nrl,1:nsl),u_ttex(1:nrl,1:nsl))
  allocate(locerr(1:nt),toterr(1:nt))
  do it = 1,nt
     t = dt*dble(it-1)
     !================================================================================
     !Set Boundaries
     !================================================================================
     u(1,:) = 2.d0*sin(5.d0*dt*dble(it-1))
     u(nrl,:) = 0.d0
     u(:,1) = 0.d0
     u(:,nsl) = 0.d0
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
     !================================================================================
     !Compute Forcing
     !================================================================================
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

     do i = 0,nrl+1
        do j = 0,nsl+1
           !f(i,j) = 2.d0*c*c*u(i,j) !For MMStrig
           f(i,j) = 0 !For MMSpoly
        end do
     end do
     !write(*,*) 'u',u(1,1:nsl)
     !write(*,*) 'f',f(1,1:nsl)
     !stop
     call MPI_Sendrecv(f(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
          f(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
     ! send to the left recieve from the right
     call MPI_Sendrecv(f(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
          f(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
     ! send to the bottom recieve from the top
     call MPI_Sendrecv(f(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
          f(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
     ! send to the right recieve from the left
     call MPI_Sendrecv(f(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
          f(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     !================================================================================
     !Compute RHS
     !================================================================================
     call compute_rhs(Ju_tt,u,f,rx,ry,sx,sy,c,hr,hs,nrl,nsl,jac)
     do i = 1,nrl
        do j = 1,nsl
           u_tt(i-1,j-1) = Ju_tt(i,j)/jac(i,j)
        end do
     end do
     call MPI_Sendrecv(u_tt(0:nrl-1,1),nrl,MPI_DOUBLE_PRECISION,p_up,123,&
          u_tt(0:nrl-1,nsl-1),nrl,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
     ! send to the left recieve from the right
     call MPI_Sendrecv(u_tt(1,0:nsl-1),nsl,MPI_DOUBLE_PRECISION,p_left,124,&
          u_tt(nrl-1,0:nsl-1),nsl,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
     ! send to the bottom recieve from the top
     call MPI_Sendrecv(u_tt(0:nrl-1,nsl-1),nrl,MPI_DOUBLE_PRECISION,p_down,125,&
          u_tt(0:nrl-1,0),nrl,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
     ! send to the right recieve from the left
     call MPI_Sendrecv(u_tt(nrl-1,0:nsl-1),nsl,MPI_DOUBLE_PRECISION,p_right,126,&
          u_tt(0,0:nsl-1),nsl,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !================================================================================
     !Update the solution at the next time level
     !================================================================================
     up(1:nrl,1:nsl) = 2.d0*u(1:nrl,1:nsl) - um(1:nrl,1:nsl) + dt*dt*u_tt
     !================================================================================
     !Shift the time levels and save output
     !================================================================================
     !call MMSpoly(u_ex,u_ttex,x,y,t,nrl,nsl)
     !call MMStrig(u_ex,u_ttex,x,y,t,nrl,nsl)
     
     do i = 1,nrl
        do j = 1,nsl
           um(i,j) = u(i,j)
           u(i,j) = up(i,j)
     !      error_tave(i,j) = abs(u_tt(i-1,j-1)-u_ttex(i,j))
        end do
     end do
     !locerr(it) = sum(error_tave(1:nrl,1:nsl))
     !call MPI_Reduce(locerr(it),toterr(it),1,&
     !     MPI_REAL8,MPI_SUM,&
     !     0,MPI_COMM_WORLD,ierr)
     !if(myid == 0) then
     !   toterr(it) = toterr(it) / real(nx * ny)
     !end if
     call MPI_Sendrecv(um(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_up,123,&
          um(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_down,123,MPI_COMM_WORLD,status,ierr)
     ! send to the left recieve from the right
     call MPI_Sendrecv(um(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
          um(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
     ! send to the bottom recieve from the top
     call MPI_Sendrecv(um(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_down,125,&
          um(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_up,125,MPI_COMM_WORLD,status,ierr)
     ! send to the right recieve from the left
     call MPI_Sendrecv(um(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
          um(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
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
     !================================================================================
     !Compute error and display
     !================================================================================
     WRITE(charit,"(I7.7)") it
     call  printdble2d(u(1:nrl,1:nsl),nrl,nsl,'u'//charit//'.txt')
     !call compute_rhs(u_tt,u,f,rx,ry,sx,sy,a,hr,hs,nrl,nsl,jac)
     !call MMStrig(u_ex,u_ttex,x,y,t,nrl,nsl)
     
  end do
  
  !write(*,*) 'time avg error: ',(sum(toterr)/dble(nt))!, u_ex
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
  !
  !call uttsolver(ucheck,rx,ry,sx,sy,a,hr,hs,jac,nrl,nsl,u_tt)
  !write(*,*) u_tt
  mpit = MPI_Wtime() - mpit
  call MPI_Reduce(mpit,totit,1,&
       MPI_REAL8,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr)
  if(myid == 0) then
     write(*,*) 'time', totit / real(px_max * py_max)
  end if
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
        if(abs(u(i,j)).lt. 1d-16) then
           write(2,fmt='(E24.16)',advance='no') 0.d0
        else
           write(2,fmt='(E24.16)',advance='no') u(i,j)
        end if
     end do
     write(2,'()')
  end do
  close(2)
end subroutine printdble2d
!compute_rhs(Ju_tt,u,f,rx,ry,sx,sy,c,hr,hs,nrl,nsl,jac)
subroutine compute_rhs(Ju_tt,u,f,rx,ry,sx,sy,c,hr,hs,nrl,nsl,jac)
  implicit none
  integer :: nrl,nsl
  real(kind = 8), dimension(0:nrl+1,0:nsl+1), intent(in) :: rx,ry,sx,sy,u,f,jac
  real(kind = 8) :: Ju_tt(0:nrl+1,0:nsl+1)
  real(kind = 8) :: c
  !real(kind = 8), dimension(:,:), allocatable, intent(inout) :: u_tt
  !real(kind = 8) :: u_tt(0:nrl-1,0:nsl-1)
  integer :: i,j
  real(kind = 8), intent(in) :: hr,hs
  real(kind = 8), dimension(:,:), allocatable :: Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8

  allocate(Q1(0:nrl+1,0:nsl+1))
  allocate(Q2(0:nrl+1,0:nsl+1))
  allocate(Q3(0:nrl+1,0:nsl+1))
  allocate(Q4(0:nrl+1,0:nsl+1))
  allocate(Q5(0:nrl+1,0:nsl+1))
  allocate(Q6(0:nrl+1,0:nsl+1))
  allocate(Q7(0:nrl+1,0:nsl+1))
  allocate(Q8(0:nrl+1,0:nsl+1))
  !allocate(Ju_tt(0:nrl+1,0:nsl+1))
  !allocate(u_tt(0:nrl-1,0:nsl-1))
  !Dimension of Ju_tt chosen to compensate for offset of ghost points

  do i = 0,nrl+1
     do j = 0,nsl+1
        Q1(i,j) = jac(i,j)*rx(i,j)*c*c*rx(i,j)
        Q2(i,j) = jac(i,j)*rx(i,j)*c*c*sx(i,j)
        Q3(i,j) = jac(i,j)*ry(i,j)*c*c*ry(i,j)
        Q4(i,j) = jac(i,j)*ry(i,j)*c*c*sy(i,j)
        Q5(i,j) = jac(i,j)*sx(i,j)*c*c*rx(i,j)
        Q6(i,j) = jac(i,j)*sx(i,j)*c*c*sx(i,j)
        Q7(i,j) = jac(i,j)*sy(i,j)*c*c*ry(i,j)
        Q8(i,j) = jac(i,j)*sy(i,j)*c*c*sy(i,j)
     end do
  end do

  !===================================================================================
  !Compute interior points using finite difference equations
  !===================================================================================
  do i = 1,nrl
     do j = 1,nsl
        Ju_tt(i,j) = (1.d0/(2.d0*hr*hr))*(u(i+1,j)*(Q1(i+1,j)+Q1(i,j))-&
             u(i,j)*(Q1(i+1,j)+2.d0*Q1(i,j)+Q1(i-1,j))+u(i-1,j)*(Q1(i,j)+Q1(i-1,j)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(4.d0*hr*hs))*(Q2(i,j+1)*(u(i+1,j+1)-&
             u(i-1,j+1))+Q2(i,j-1)*(u(i-1,j-1)-u(i+1,j-1)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(2.d0*hr*hr))*(u(i+1,j)*(Q3(i+1,j)+&
             Q3(i,j))-u(i,j)*(Q3(i+1,j)+2.d0*Q3(i,j)+Q3(i-1,j))+u(i-1,j)*(Q3(i,j)+Q3(i-1,j)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(4.d0*hr*hs))*(Q4(i,j+1)*(u(i+1,j+1)-&
             u(i-1,j+1))+Q4(i,j-1)*(u(i-1,j-1)-u(i+1,j-1)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(4.d0*hr*hs))*(Q5(i+1,j)*(u(i+1,j+1)-&
             u(i+1,j-1))+Q5(i-1,j)*(u(i-1,j-1)-u(i-1,j+1)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(2.d0*hs*hs))*(u(i,j+1)*(Q6(i,j+1)+&
             Q6(i,j))-u(i,j)*(Q6(i,j+1)+2.d0*Q6(i,j)+Q6(i,j-1))+u(i,j-1)*(Q6(i,j)+Q6(i,j-1)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(4.d0*hr*hs))*(Q7(i+1,j)*(u(i+1,j+1)-&
             u(i+1,j-1))+Q7(i-1,j)*(u(i-1,j-1)-u(i-1,j+1)))
        Ju_tt(i,j) = Ju_tt(i,j) + (1.d0/(2.d0*hs*hs))*(u(i,j+1)*(Q8(i,j+1)+&
             Q8(i,j))-u(i,j)*(Q8(i,j+1)+2.d0*Q8(i,j)+Q8(i,j-1))+u(i,j-1)*(Q8(i,j)+Q8(i,j-1)))
        Ju_tt(i,j) = Ju_tt(i,j) + jac(i,j)*f(i,j)
     end do
  end do

  !do i = 0,nrl-1
  !   do j = 0,nsl-1
  !      u_tt(i,j) = Ju_tt(i+1,j+1)/jac(i+1,j+1)
  !   end do
  !end do

  !deallocate(Q1)
  !deallocate(Q2)
  !deallocate(Q3)
  !deallocate(Q4)
  !deallocate(Q5)
  !deallocate(Q6)
  !deallocate(Q7)
  !deallocate(Q8)
  !deallocate(Ju_tt)

end subroutine compute_rhs

subroutine MMStrig(u_ex,u_ttex,x,y,t,nrl,nsl)
  implicit none
  integer :: nrl,nsl,i,j
  real(kind=8), dimension(0:nrl+1,0:nsl+1), intent(in) :: x,y
  real(kind=8) :: u_ex(1:nrl,1:nsl),u_ttex(1:nrl,1:nsl)
  real(kind=8) :: t

  do i = 1,nrl
     do j = 1,nsl
        u_ex(i,j) = sin(t-x(i,j))*sin(y(i,j))
        u_ttex(i,j) = 0
     end do
  end do
end subroutine MMStrig

subroutine MMSpoly(u_ex,u_ttex,x,y,t,nrl,nsl)
  implicit none
  integer :: nrl,nsl,i,j
  real(kind=8), dimension(0:nrl+1,0:nsl+1), intent(in) :: x,y
  real(kind=8) :: u_ex(1:nrl,1:nsl),u_ttex(1:nrl,1:nsl)
  real(kind=8) :: t

  do i = 1,nrl
     do j = 1,nsl
        u_ex(i,j) = (x(i,j)+x(i,j)*x(i,j))*&
             (y(i,j)+y(i,j)*y(i,j))*&
             (t+t*t)
        u_ttex(i,j) = 0
     end do
  end do

end subroutine MMSpoly
