!What shoud be a super-basic example of Trapazoidal rule approximation
!(except this is my first fortran attempt...)
!
!
program trapappx

  implicit none
  double precision :: error_n1,h,y,x_l,x_r,x_i,y_i,N,y_appx,y_out,error_n0,r,rr
  real(kind = 8), external :: f
  integer :: i,iter
  integer :: end_criteria

  !Here we try to approximate the integration
  ! I = int(-1,1){exp(cos(kx))}dx
  !Where k = pi or pi^2

  x_l = -1.d0
  x_r = 1.d0
  N = 2.d0!Start with 2 because the terms x_0 and x_1 are defined above as x_l and x_r, respectivly
  y = (f(x_l)+f(x_r))/2.d0
  h = (x_r-x_l)/N
  error_n1 = 1.d0!ambiguous assignment of error value so that the loop runs
  y_appx = 100!again, ambiguous assignment to run loop
  do iter = 1, 1325
     end_criteria = iter
     y_i = 0.d0
     y_out = 0.d0
     do i = 1,end_criteria
        x_i = x_l + (dble(i)*h)
        y_i = f(x_i)
        y_out = y_out + y_i
     end do
     error_n0 = error_n1
     error_n1 = abs(y_appx - h*(y + y_out))
     r = error_n1/error_n0
     r = error_n1/(error_n0*error_n0)
     y_appx = 0.d0
     y_appx = h*(y + y_out)
     N = N + 1.d0
     h = (x_r-x_l)/N
     write(*,*) iter, y_appx,  error_n1

  end do
  !print *, y_appx
  !print *, iter

 ! do i = 1,n
 !    h = (x_r-x_l)/N
 !    N = N + 1
 !    x_i = x_l+dble(i)*h
 !    y_i = y_i+f(x_i)

end program trapappx

real(kind = 8) function f(x)
  implicit none
  real(kind = 8), intent(in) :: x
  real(kind = 8), parameter :: pi = 3.141592653589793238462643383279d0
  real(kind = 8):: k
  k = pi * pi
  f = exp(cos(k*x))
end function f
