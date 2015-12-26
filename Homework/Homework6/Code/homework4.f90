program hwk4
  use xycoord ! use the module xycoord to set the mapping 
  implicit none
  integer :: nr,ns,i,j
  double precision :: integral
  real(kind = 8) :: hr,hs
  real(kind = 8), external :: f
  real(kind = 8), dimension(:), allocatable :: r,s
  real(kind = 8), dimension(:,:), allocatable :: u,ur,us,xcr,xcs,ycr,ycs,sx,rx,sy,ry
  
  real(kind = 8), dimension(:,:), allocatable :: xc,yc
  ! Number of grid points for the cartesian grid
  nr = 50
  ns = 50
  
  ! Allocate memory for the various arrays
  allocate(r(0:nr),s(0:ns),u(0:nr,0:ns),ur(0:nr,0:ns),us(0:nr,0:ns),xcr(0:nr,0:ns),ycr(0:nr,0:ns),xcs(0:nr,0:ns),ycs(0:nr,0:ns))
  allocate(xc(0:nr,0:ns),yc(0:nr,0:ns),sx(0:nr,0:ns),sy(0:nr,0:ns),rx(0:nr,0:ns),ry(0:nr,0:ns))
  
  hr = 2.d0/dble(nr)
  hs = 2.d0/dble(ns)
  do i = 0,nr
     r(i) = -1.d0 + dble(i)*hr
  end do
  do i = 0,ns
     s(i) = -1.d0 + dble(i)*hs
  end do

  do j = 0,ns
     do i = 0,nr
        xc(i,j) = x_coord(r(i),s(j))
        yc(i,j) = y_coord(r(i),s(j))
     end do
  end do
  call  printdble2d(xc,nr,ns,'x.txt')
  call  printdble2d(yc,nr,ns,'y.txt')
  ! 
  do j = 0,ns
     do i = 0,nr
        !u(i,j) = sin(r(i))*cos(s(j))
        !u(i,j) = sin(xc(i,j)) * cos(yc(i,j))
        !u(i,j) = exp(xc(i,j) + yc(i,j))
        !u(i,j) = (r(i))**2 + (s(j))**2
        u(i,j) = xc(i,j)**2 + yc(i,j)**2
     end do
  end do
  ! Differentiate in the r-direction
  do i = 0,ns
     call differentiate(u(0:nr,i),ur(0:nr,i),hr,nr)
  end do
  do i = 0,ns
     call differentiate(xc(0:nr,i),xcr(0:nr,i),hr,nr)
  end do
  do i = 0,ns
     call differentiate(yc(0:nr,i),ycr(0:nr,i),hr,nr)
  end do

  ! Differentiate in the s-direction
  do i = 0,nr
     call differentiate(u(i,0:ns),us(i,0:ns),hs,ns)
  end do
  do i = 0,nr
     call differentiate(xc(i,0:ns),xcs(i,0:ns),hs,ns)
  end do
  do i = 0,nr
     call differentiate(yc(i,0:ns),ycs(i,0:ns),hs,ns)
  end do

  do i = 0,nr
     do j = 0,ns
        sx(i,j) = -ycr(i,j)/(xcr(i,j)*ycs(i,j)-xcs(i,j)*ycr(i,j))
        rx(i,j) = ycs(i,j)/(xcr(i,j)*ycs(i,j)-xcs(i,j)*ycr(i,j))
        sy(i,j) = xcr(i,j)/(xcr(i,j)*ycs(i,j)-xcs(i,j)*ycr(i,j))
        ry(i,j) = -xcs(i,j)/(xcr(i,j)*ycs(i,j)-xcs(i,j)*ycr(i,j))
     end do
  end do

  !sx = -yr/(xr*ys-xs*yr)
  !rx = ys/(xr*ys-xs*yr)
  !sy = xr/(xr*ys-xs*yr)
  !ry = -xs/(xr*ys-xs*yr)

  call trap_2D(f,-1.d0,1.d0,-1.d0,1.d0,integral,nr,ns)
  print *, integral
end program hwk4

!real(kind = 8) function f(r, s)
!  implicit none
!  real(kind = 8), intent(in) :: r,s
!  real(kind = 8) :: u,ur,us,x,xr,xs,y,yr,ys,hr,hs
!  integer :: nr,ns

!  nr = 30
!  ns = 60

!  hr = 2.d0 / dble(nr)
!  hs = 2.d0 / dble(ns)
  
!  x = r + 0.1d0*s
!  call differentiate(x,xr,hr,nr)
!  call differentiate(x,xs,hs,ns)
!  y = s
!  call differentiate(y,yr,hr,nr)
!  call differentiate(y,ys,hs,ns)
!  u = sin(x)*cos(y)
!  call differentiate(u,ur,hr,nr)
!  call differentiate(u,us,hs,ns)
!  f = r + s
!    f = ((ys/(xr*ys-xs*yr))*ur+(-yr/(xr*ys-xs*yr))*us+(-xs/(xr*ys-xs*yr))*ur+(xr/(xr*ys-xs*yr))*us-(cos(x)-sin(y)))**2
!end function f

real(kind = 8) function f(r, s)
  implicit none
  real(kind = 8), intent(in) :: r,s
  
  f = 1
end function f

subroutine trap_2D(f,a,b,c,d,integral,nx,ny)
implicit none
double precision f, a, b, c, d, integral, sum
double precision hx, hy, x, y
integer nx, ny
integer i, j

hx = (b-a)/dfloat(nx-1)
hy = (d-c)/dfloat(ny-1)
   
! calculate the corner's terms
sum = f(a,c)+f(a,d)+f(b,c)+f(b,d)

! calculate single sums
do i=2,nx-1
   x = a + hx*(i-1)
   sum = sum + 2.0*(f(x,c)+f(x,d))
end do

do j=2, ny-1
   y = c + hy*(j-1)
   sum = sum + 2.0*(f(a,y)+f(b,y))
end do

! calculate the double sum
do i = 2,nx-1
   x = a + hx*(i-1)
   do j = 2,ny-1
      y = c + hy*(j-1)
      sum = sum + 4.0*f(x,y)
   end do
end do
integral = 0.25*hx*hy*sum
return
end subroutine trap_2D
