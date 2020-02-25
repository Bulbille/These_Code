! Transformee de Fourier complexe-complexe 2d

program tjmccfft2d

  implicit none

  integer, parameter              :: m = 1024, n = 1024
  complex, dimension(0:n-1,0:m-1) :: x
  real, dimension(100+2*(n+m))    :: table
  real, dimension(4096*max(n,m))  :: work

  integer :: isys
  integer :: i, j

  ! On prepare le tableau d'entree
  do j = 0,m-1
    do i = 0,n-1
      x(i,j) = cmplx(j*j*j+i*i,i*i*i+j*j)
    end do
  end do

  call ccfft2d(0,n,m,1.d0,x,n,x,n,table,work,isys)
  call ccfft2d(1,n,m,1.d0,x,n,x,n,table,work,isys)

end program tjmccfft2d
