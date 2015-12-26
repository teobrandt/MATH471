module MMS
  real(kind = 8), parameter :: pi = acos(-1.d0)
contains
  
  subroutine MMStrig(u_ex,u_tt,x,y,t,nrl,nsl)
    implicit none
    integer :: nrl,nsl
    real(kind=8), dimension(0:nrl+1,0:nsl+1), intent(in) :: x,y
    real(kind=8) :: u_ex(0:nrl-1),u_tt(0:nrl-1)
    real(kind=8), intent(in) :: t
    
    do i = 0,nrl-1
       do j = 0,nsl-1
          u_ex(i,j) = sin(t-x(i+1,j+1))*sin(y(i+1,j+1))
          u_tt(i,j) = 0
       end do
    end do
  end subroutine MMStrig

  subroutine MMSpoly(u_ex,u_tt,x,y,t,nrl,nsl)
    implicit none
    integer :: nrl,nsl
    real(kind=8), dimension(0:nrl+1,0:nsl+1), intent(in) :: x,y
    real(kind=8) :: u_ex(0:nrl-1),u_tt(0:nrl-1)
    real(kind=8), intent(in) :: t
    
    do i = 0,nrl-1
       do j = 0,nsl-1
          u_ex(i,j) = (x(i+1,j+1)+x(i+1,j+1)*x(i+1,j+1))&
               *(y(i+1,j+1)+y(i+1,j+1)*y(i+1,j+1))&
               *(t(i+1,j+1)+t(i+1,j+1)*t(i+1,j+1))
          u_tt(i,j) = 0
       end do
    end do
    
  end subroutine MMSpoly
    
end module MMS
