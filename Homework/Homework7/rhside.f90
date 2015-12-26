subroutine rhside(fu,fv,u,v,t,a,nx,hx)
  ! This subroutine computes an approximation to 
  ! the time derivative (i.e. the right hand side) 
  ! of the wave equation on the form
  ! 
  ! u_t = v_t, 
  ! v_t = (a(x)*u_x)_x.
  !
  ! with dirichlet boundary conditions at both ends. 
  !
  implicit none
  integer, intent(in) :: nx
  real(kind = 8), intent(in) :: t,hx
  real(kind = 8), dimension(0:nx), intent(in) :: u,v,a
  real(kind = 8), dimension(1:nx-1), intent(inout) :: fu,fv
  real(kind = 8), dimension(0:nx-1) :: ap_ave
  real(kind = 8), dimension(0:nx-1) :: ux
  real(kind = 8), dimension(1:nx-1) :: uxx
  real(kind = 8) :: hxi
  integer :: i
  ! Avoid division.
  hxi = 1.d0/hx
  ! Compute D_+(u) and E_{1/2}(a)
  do i = 0,nx-1
     ap_ave(i) = 0.5d0*(a(i+1)+a(i))
     ux(i) = hxi*(u(i+1)-u(i))
  end do
  ! Compute D_-(E_{1/2}(a) D_+(u))
  do i = 1,nx-1
     uxx(i) = hxi*(ap_ave(i)*ux(i)-ap_ave(i-1)*ux(i-1))
  end do
  
  fu = v(1:nx-1)  
  fv = uxx

end subroutine rhside
