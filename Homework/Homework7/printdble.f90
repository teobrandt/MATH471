subroutine printdble1d(u,nx,str)
  implicit none
  integer, intent(in) :: nx
  real(kind = 8), intent(in) :: u(0:nx)
  character(len=*), intent(in) :: str
  integer :: i
  open(2,file=trim(str),status='unknown')
  do i=0,nx,1
     write(2,fmt='(E24.16)',advance='no') u(i)
  end do
  close(2)
end subroutine printdble1d

subroutine printdble2d(u,nx,ny,str)
  implicit none
  integer, intent(in) :: nx,ny
  real(kind = 8), intent(in) :: u(0:nx,0:ny)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=0,ny,1
     do i=0,nx,1
        write(2,fmt='(E24.16)',advance='no') u(i,j)
     end do
     write(2,'()')
  end do
  close(2)
end subroutine printdble2d
