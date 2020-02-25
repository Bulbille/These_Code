! Transformee de Fourier complexe-complexe 2d

program tjmccfft2d

  implicit none

  integer, parameter :: m = 18

  integer, parameter :: n = 24
  complex, dimension(0:n-1,0:m-1) :: x, xx
  real, dimension(0:n-1,0:m-1,2) :: rx
  equivalence ( x, rx )

  integer, parameter :: ntable = 100+2*(n+m)
  real, dimension(0:ntable-1) :: table

  ! Les routines jm ont besoin de 4*n*m. D'ou tronconnage
  integer, parameter :: nwork = 512*max(n,m)
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k, l
  real :: twopi
  complex :: s

  twopi = 2 * acos(real(-1))

  ! On prepare le tableau d'entree
  call random_number(rx)
  xx = x

  isys = 0
  scale = 1./real(n)

  isign = 0
  call ccfft2d(isign,n,m,scale,x,n,x,n,table,work,isys)
  isign = 1
  print *,'jmccfft2d ',n,m,isign,scale
  call ccfft2d(isign,n,m,scale,x,n,x,n,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') x

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  ! On reprepare le tableau d'entree
  x = xx
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n-1
      s = 0
      do l = 0,m-1
        do k = 0,n-1
           s = s+   cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n))) &
           &      * cmplx(cos(twopi*j*l/real(m)),isign*sin(twopi*j*l/real(m))) &
           &      * x(k,l)
        end do
      end do
      write(11,'(2e25.12)') s*scale
    end do
  end do

end program tjmccfft2d
