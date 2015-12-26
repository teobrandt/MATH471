program getmatrix
  ! 
  ! 
  implicit none
  integer :: nx,i,N
  real(kind = 8) :: hx,xl,xr,t,amp
  real(kind = 8), dimension(:), allocatable :: x,u,v,a,col
  real(kind = 8), dimension(:,:), allocatable :: k1,AMAT,acopy
  ! Stuff for computing the eigenvalues
  integer :: LDA, LDVL, LDVR, LWORK, INFO
  real(kind = 8), dimension(:), allocatable :: WR,WI,WORK
  real(kind = 8), dimension(:,:), allocatable :: VL,VR

  nx = 100
  xl = -1.0d0
  xr =  2.d0
  amp = 0.0d0
  
  ! Allocate memory for the various arrays
  allocate(x(0:nx),a(0:nx),u(0:nx),v(0:nx))
  allocate(k1(nx-1,2),col(2*(nx-1)),AMAT(2*(nx-1),2*(nx-1)))
  ! Set up grid, initial data and a(x)
  hx = (xr-xl)/dble(nx)
  do i = 0,nx
     x(i) = xl + dble(i)*hx
     a(i) = 1.d0+amp*sin((x(i)+1.d0)**4)
  end do
  
  u = 0.d0
  v = 0.d0
  t = 0.d0
  do i = 1,2*(nx-1)
     col = 0.d0;
     col(i) = 1.d0
     u(1:nx-1) = col(1:nx-1)
     v(1:nx-1) = col(nx:2*nx-2)
     call rhside(k1(1:nx-1,1),k1(1:nx-1,2),u,v,t,a,nx,hx)
     AMAT(1:nx-1,i)    = k1(1:nx-1,1)
     AMAT(nx:2*nx-2,i) = k1(1:nx-1,2)
  end do
  call printdble2d(Amat,2*nx-3,2*nx-3,'amat.txt')
  
  ! This is an example of a 

  N = 2*(nx-1)
  LDA = N
  LDVL = 1
  LDVR = 1
  LWORK = 10*N
  allocate(Acopy(N,N),WR(N),WI(n),VL(LDVL,N),VR(LDVL,N),WORK(LWORK))
  Acopy = AMAT
  call dgeev ('N','N', N, Acopy, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
  
  call printdble1d(WR,N-1,'lam_re.txt')
  call printdble1d(WI,N-1,'lam_im.txt')

  deallocate(ACOPY,WR,WI,VL,VR,WORK)

  ! Make up a right hand side B by doing matmul(AMAT,YEX) 
  ! For some YEX, then: 
  ! 1. Factor A = PLU using the Lapack routine DGETRF, NOTE that A is overwritten!  
  ! 2. Solve  A*Y = PLU*Y = B  DGETRS
  ! 3. Make sure Y-YEX is really small. 
  
  deallocate(x,a,k1,amat,u,v)

end program getmatrix
