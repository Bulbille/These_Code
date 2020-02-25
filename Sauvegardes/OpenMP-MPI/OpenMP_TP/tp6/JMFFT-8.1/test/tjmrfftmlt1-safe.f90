! Transformee de Fourier complexe-reel multiple

program tjmrfftmlt1

  implicit none

  integer, parameter :: m = 7

  integer, parameter :: n = 36
  complex, dimension(0:n/2,0:m-1) :: x, xx
  real, dimension(0:n/2,0:m-1,2) :: rx
  equivalence ( x, rx )
  real, dimension(0:2*(n/2+1)-1,0:m-1) :: y
  equivalence(x(0,0),y(0,0))

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

  ! On prepare le tableau d'entree sans forcer a 0 les termes necessaires
  call random_number( rx )
  xx = x

  call fftfax(n,ifax,trigs)
  isign = +1
  print *,'jmrfftmlt ',n,m,isign
  inc = 1
  jump = 2*(n/2+1)
  call rfftmlt(x,work,trigs,ifax,inc,jump,n,m,isign)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') y(0:n-1,0:m-1)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')

  ! On reprepare le tableau d'entree en forcant a 0 les termes necessaires
  x = xx
  x(0,:) = cmplx( real( x(0,:) ), 0 )
  if ( mod( n, 2 ) == 0 ) x(n/2,:) = cmplx( real( x(n/2,:) ), 0 )
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n-1
      s = 0
      do k = 0,n/2
        s = s+cmplx(cos(twopi*i*k/real(n)),sin(twopi*i*k/real(n)))*x(k,j)
      end do
      do k = n/2+1,n-1
        s = s+cmplx(cos(twopi*i*k/real(n)),sin(twopi*i*k/real(n)))* &
        & conjg(x(n-k,j))
      end do
      write(11,'(e25.12)') real(s)
    end do
  end do
  close(11)

end program tjmrfftmlt1
