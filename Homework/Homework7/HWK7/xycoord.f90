module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
contains
  
  real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
    !x_coord = -2.d0+ (2.d0+r+0.2d0*sin(5.d0*s))*cos(s) 
    !x_coord = r-3.d0*s
    x_coord = 4.d0*r+2.d0*sin(r)*0.05d0*sin(r)
  end function x_coord

  real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s
    !y_coord = (2.d0+r+0.3d0*sin(5.d0*s))*sin(s) 
    !y_coord = s+r
    y_coord = 1.d0*sin(3.d0*r)+2.d0*sin(s)
  end function y_coord
    
end module xycoord
