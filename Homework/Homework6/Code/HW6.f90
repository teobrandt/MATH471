program HW6
  !$ use omp_lib
  implicit none
  integer :: nr, ns, i, j
  real(kind = 8), parameter :: pi = 3.141592653589793238462643383279d0
  real(kind = 8), dimension(:), allocatable :: rc,sc
  real(kind = 8), dimension(:,:), allocatable :: X,Y,R,S
  real(kind = 8), dimension(:,:), allocatable :: u,ux,uy,uex_x,uex_y
  real(kind = 8), dimension(:,:), allocatable :: ur,Xr,Yr,us,Xs,Ys,Rx,Ry,Sx,Sy
  real(kind = 8), dimension(:,:), allocatable :: integrand,RAD,JAC
  real(kind = 8), dimension(:), allocatable :: I1
  real(kind = 8) :: hr,hs,begin,end,error,Ie
  !Set the number of cores to be used for parallel processing
  !$ call OMP_set_num_threads(1)
  do nr = 200,200
  
     ns = nr
  !Allocate statements will be done ahead of time so that computational time includes assignments
     allocate(rc(0:nr),sc(0:ns),I1(0:ns))
     allocate(integrand(0:nr,0:ns),JAC(0:nr,0:ns))
     allocate(X(0:nr,0:ns),Y(0:nr,0:ns),R(0:nr,0:ns),S(0:nr,0:ns),RAD(0:nr,0:ns))
     allocate(u(0:nr,0:ns),uex_x(0:nr,0:ns),uex_y(0:nr,0:ns),ux(0:nr,0:ns),uy(0:nr,0:ns))
     allocate(ur(0:nr,0:ns),Xr(0:nr,0:ns),Yr(0:nr,0:ns),us(0:nr,0:ns),Xs(0:nr,0:ns),Ys(0:nr,0:ns))
     allocate(Rx(0:nr,0:ns),Ry(0:nr,0:ns),Sx(0:nr,0:ns),Sy(0:nr,0:ns))
  
!Begin timing the computational part of the program
     !$ begin = omp_get_wtime()
!Set up the cartesian grid spacing
     hr = 2.d0/dble(nr)
     hs = 2.d0/dble(ns)
!Define grids in r and s
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,nr
        rc(i) = -1.d0 + dble(i)*hr
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,ns
        sc(i) = -1.d0 + dble(i)*hs
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,nr
        !$OMP PARALLEL DO PRIVATE(j)
        do j = 0,ns
           R(i,j) = rc(i)
           S(i,j) = sc(j)
           X(i,j) = (rc(i)+2.d0)*cos(sc(j)*pi/2.d0)
           Y(i,j) = (rc(i)+2.d0)*sin(sc(j)*pi/2.d0)
        end do
        !$OMP END PARALLEL DO
     end do
     !$OMP END PARALLEL DO
     call printdble2d(X,nr,ns,'x.txt')
     call printdble2d(Y,nr,ns,'y.txt')

     !$OMP PARALLEL DO PRIVATE(j)
     do j = 0,ns
        !$OMP PARALLEL DO PRIVATE(i)
        do i = 0,nr
           RAD = sqrt(X(i,j)**2.d0+Y(i,j)**2)
           u(i,j) = -(RAD(i,j)*(RAD(i,j)-3))
           uex_x = 2.d0*X(i,j)
           uex_y = 2.d0*Y(i,j)
        end do
        !$OMP END PARALLEL DO
     end do
     !$OMP END PARALLEL DO

  ! Differentiate in the r-direction
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,ns
        call differentiate(u(0:nr,i),ur(0:nr,i),hr,nr)
        call differentiate(X(0:nr,i),Xr(0:nr,i),hr,nr)
        call differentiate(Y(0:nr,i),Yr(0:nr,i),hr,nr)
     end do
     !$OMP END PARALLEL DO

  ! Differentiate in the s-direction
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,nr
        call differentiate(u(i,0:ns),us(i,0:ns),hs,ns)
        call differentiate(X(i,0:ns),Xs(i,0:ns),hs,ns)
        call differentiate(Y(i,0:ns),Ys(i,0:ns),hs,ns)
     end do
     !$OMP END PARALLEL DO 
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,nr
        !$OMP PARALLEL DO PRIVATE(j)
        do j = 0,ns
           Sx(i,j) = -Yr(i,j)/(Xr(i,j)*Ys(i,j)-Xs(i,j)*Yr(i,j))
           Rx(i,j) = Ys(i,j)/(Xr(i,j)*Ys(i,j)-Xs(i,j)*Yr(i,j))
           Sy(i,j) = Xr(i,j)/(Xr(i,j)*Ys(i,j)-Xs(i,j)*Yr(i,j))
           Ry(i,j) = -Xs(i,j)/(Xr(i,j)*Ys(i,j)-Xs(i,j)*Yr(i,j))
           ux(i,j) = Rx(i,j)*ur(i,j)+Sx(i,j)*us(i,j)
           uy(i,j) = Ry(i,j)*ur(i,j)+Sy(i,j)*us(i,j)
           JAC(i,j) = Xr(i,j)*Ys(i,j)-Xs(i,j)*Yr(i,j)
        end do
        !$OMP END PARALLEL DO
     end do
     !$OMP END PARALLEL DO

  ! Constructing the integrand
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 0,nr
        !$OMP PARALLEL DO PRIVATE(j)
        do j = 0,ns
           Integrand(i,j) = u(i,j)*JAC(i,j)
        end do
        !$OMP END PARALLEL DO
     end do
     !$OMP END PARALLEL DO
  
!TIME TO INTEGRATE!!!!!
     !$OMP PARALLEL DO PRIVATE(j)
     do j = 1,ns
        I1(j) = hr*sum(Integrand(0:nr,j))
     end do
     !$OMP END PARALLEL DO
!SECOND INTEGRATION IN THE r DIRECTION
     error = hr*sum(I1(0:ns))
     Ie = (8.d0/3.d0)*pi
     error = abs(error - Ie)
     !$ end = omp_get_wtime()
     !$ write(*,*) end - begin, nr
     !write(*,*) '10 number of strings (20x20)'
     deallocate(rc,sc,I1)
     deallocate(integrand)
     deallocate(X,Y,R,S)
     deallocate(u,uex_x,uex_y,ux,uy)
     deallocate(ur,Xr,Yr,us,Xs,Ys)
     deallocate(Rx,Ry,Sx,Sy)

  end do
  write(*,*) '1 string(s)', error/100
end program HW6
