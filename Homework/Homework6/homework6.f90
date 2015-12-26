program hwk6
  ! 
  ! 
  ! 
  implicit none
  integer :: nx,i,nsteps,it
  real(kind = 8) :: hx,xl,xr,dt,tend,cfl,t,amp
  real(kind = 8), dimension(:), allocatable :: x,u,v,a
  ! For Runge-Kutta, use the second dimension as 1 = u, 2 = v.
  real(kind = 8), dimension(:,:), allocatable :: k1,k2,k3,k4
  real(kind = 8), dimension(:), allocatable :: up,vp
  real(kind = 8), parameter :: pi = acos(-1.d0) 
  CHARACTER(7) :: charit

  nx = 100
  xl = -1.0d0
  xr =  2.d0
  tend = 5.d0
  cfl = 0.1d0

  amp = 0.1d0

  ! Allocate memory for the various arrays
  allocate(x(0:nx),a(0:nx),u(0:nx),v(0:nx))
  allocate(k1(nx-1,2),k2(nx-1,2),k3(nx-1,2),k4(nx-1,2))
  allocate(up(0:nx),vp(0:nx))

  ! Set up grid, initial data and a(x)
  hx = (xr-xl)/dble(nx)
  do i = 0,nx
     x(i) = xl + dble(i)*hx
     ! initial data
     u(i) = exp(-36.d0*x(i)*x(i))!cos(x(i))
     v(i) = 0.d0
     a(i) = 1.d0+amp*sin((x(i)+1.d0)**4)
  end do
  call printdble1d(x,nx,'x.txt')
  call printdble1d(a,nx,'a.txt')
  dt = cfl*hx
  nsteps = floor(tend/dt) + 1
  dt = tend/dble(nsteps)

  ! Save the initial data to file
  it = 0
  WRITE(charit,"(I7.7)") it
  call printdble1d(u,nx,'u'//charit//'.txt')
  call printdble1d(v,nx,'v'//charit//'.txt')

  do it = 1,nsteps
     t = dble(it-1)*dt

     ! This is Heun's method which is NOT STABLE!!! for this problem.
     ! Replace it with RK4.
     ! Stage 1, compute k1
     up = u
     vp = v
     call rhside(k1(1:nx-1,1),k1(1:nx-1,2),up,vp,t,a,nx,hx)
     ! Stage 2, compute k2
     t = t + 0.5*dt
     up(1:nx-1) = u(1:nx-1) + 0.5d0*dt*k1(:,1)
     vp(1:nx-1) = v(1:nx-1) + 0.5d0*dt*k1(:,2)
     call rhside(k2(1:nx-1,1),k2(1:nx-1,2),up,vp,t,a,nx,hx)

     ! Combine the stages to find the approximate solution at 
     ! "t = t_n + dt"

     t = dble(it)*dt
     u(1:nx-1) = u(1:nx-1) &
          + (dt/2.d0)*(k1(:,1)+k2(:,1))
     v(1:nx-1) = v(1:nx-1) &
          + (dt/2.d0)*(k1(:,2)+k2(:,2))

     WRITE(charit,"(I7.7)") it
     call printdble1d(u,nx,'u'//charit//'.txt')
     call printdble1d(v,nx,'v'//charit//'.txt')
     
  end do
  
  deallocate(x,a,u,v,k1,k2,k3,k4,up,vp)

end program hwk6
