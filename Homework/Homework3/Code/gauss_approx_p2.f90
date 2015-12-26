program gauss_approx_p2
  implicit none
  real(kind = 8), parameter :: pi = acos(-1.d0)
  real(kind = 8), dimension(:), allocatable :: w, x, y
  real(kind = 8) :: Integral_value, last_val, error_n
  real(kind = 8) :: dE = 0.d0
  integer :: n

  Integral_value = 1.d0
  
  do n = 1, 1325
     last_val = Integral_value
     
     allocate(x(0:n))
     allocate(y(0:n))
     allocate(w(0:n))

     call lglnodes(x,w,n)     
     y = exp(cos(pi*x))
     Integral_value = sum(y*w)

     deallocate(w)
     deallocate(x)
     deallocate(y)

     dE = abs(Integral_value - last_val) ! relative absolute error
          
     write(*,*) n, Integral_value, dE
  end do

  
end program gauss_approx_p2
