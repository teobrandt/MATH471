!
! A super-basic example of solving an equation
! using Newton's method
!
! This is the template file used for the scripted version  
!
program newton
  
  implicit none
  double precision :: f,fp,x,dx,nx,snx
  integer :: iter
  
  ! Here we try to find the solution to the
  ! f(x) = 0
  
  x = -0.5d0
  dx = 1.d0
  iter = 0
  do while(abs(dx) > 10E-15)
     iter = iter + 1
     f = ffun(x)
     fp = fpfun(x)
     dx = -f/fp
     x = x + dx
     nx = abs(x)/abs(x-dx)
     snx = nx/abs(x-dx)
     write(*, 950) ' sin(x)+cos(x*x) ', ' & ', iter, ' & ', x, ' & ', dx, ' & ', nx, ' & ', snx, ' \\'
     950 format(A20, A3, I2.2, 4(A3, E24.16), A3)
  end do

contains

  double precision function ffun(x)
    implicit none
    double precision :: x

    ffun = sin(x)+cos(x*x)

  end function ffun

  double precision function fpfun(x)
    implicit none
    double precision :: x

    fpfun = cos(x)-2.d0*x*sin(x*x)

  end function fpfun

end program newton
