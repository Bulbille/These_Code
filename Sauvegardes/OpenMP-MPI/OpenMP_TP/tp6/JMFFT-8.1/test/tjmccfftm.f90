! Transformee de Fourier complexe-complexe multiple

program tjmccfftm

  implicit none

  integer, parameter :: m = 7

  integer, parameter :: n = 36
  complex, dimension(0:n-1,0:m-1) :: x, xx
  real, dimension(0:n-1,0:m-1,2) :: rx
  equivalence ( x, rx )

  integer, parameter :: ntable = 100+2*n
  real, dimension(0:ntable-1) :: table

  integer, parameter :: nwork = 4*m*n
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k
  real :: twopi
  complex :: s

  twopi = 2 * acos(real(-1))

  ! On prepare le tableau d'entree
  call random_number( rx ) 
  xx = x

  isys = 0
  scale = 1./real(n)

  isign = 0
  call ccfftm(isign,n,m,scale,x,n,x,n,table,work,isys)
  isign = 1
  print *,'jmccfftm ',n,m,isign,scale
  call ccfftm(isign,n,m,scale,x,n,x,n,table,work,isys)

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
      do k = 0,n-1
        s = s+ cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))*x(k,j)
      end do
      write(11,'(2e25.12)') s*scale
    end do
  end do

end program tjmccfftm
