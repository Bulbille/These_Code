! Transformee de Fourier reel-complexe multiple

program tjmrfftmlt2

  implicit none

  integer, parameter :: m = 7

  integer, parameter :: n = 36
  real, dimension(0:n+1,0:m-1) :: x, xx

  integer, parameter :: ntrigs = 2*n
  real, dimension(0:ntrigs-1) :: trigs

  integer, parameter :: nifax = 19
  real, dimension(0:nifax-1) :: ifax

  integer, parameter :: nwork = 2*m*n
  real, dimension(0:nwork-1) :: work

  integer :: isign
  integer :: i, j, k
  real :: twopi
  complex :: s
  integer :: inc, jump

  twopi = 2 * acos(real(-1))

  ! On prepare le tableau d'entree
  call random_number( x )
  xx = x

  call fftfax(n,ifax,trigs)
  isign = -1
  print *,'jmrfftmlt ',n,m,isign
  inc = 1
  jump = n+2
  call rfftmlt(x,work,trigs,ifax,inc,jump,n,m,isign)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') x(0:2*(n/2+1)-1,0:m-1)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')

  ! On reprepare le tableau d'entree
  x = xx
  do j = 0,m-1
    do i = 0,n/2
      s = 0
      do k = 0,n-1
        s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))*x(k,j)
      end do
      write(11,'(2e25.12)') s/real(n)
    end do
  end do

end program tjmrfftmlt2
