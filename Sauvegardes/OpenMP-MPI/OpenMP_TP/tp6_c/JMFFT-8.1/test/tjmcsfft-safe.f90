! Transformee de Fourier complexe-reelle simple

program tjmcsfft

  implicit none

  integer, parameter :: n = 8
  complex, dimension(0:n/2) :: x, xx
  real, dimension(0:n/2,2) :: rx
  equivalence ( x, rx )
  real, dimension(0:n-1) :: y

  ! Pour stocker les cosinus et les sinus
  ! En fait 100+2*n suffisent pour les routines jm
  integer, parameter :: ntable = 100+4*n
  real, dimension(0:ntable-1) :: table

  ! En fait, 2*n suffisent pour les routines jm
  integer, parameter :: nwork = 4+4*n
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, k
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree sans forcer a 0 les termes necessaires
  call random_number( rx )
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call csfft(isign,n,scale,x,y,table,work,isys)
  isign = 1
  print *,'jmcsfft ',n,isign,scale
  call csfft(isign,n,scale,x,y,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') y
  close(10)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree en forcant a 0 les termes necessaires
  x = xx
  x(0) = cmplx( real( x(0) ), 0 )
  if ( mod( n, 2 ) == 0 ) x(n/2) = cmplx( real( x(n/2) ), 0 )
  xx = x
  ! Et on calcule
  do i = 0,n-1
    s = 0
    do k = 0,n/2
      s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))*x(k)
    end do
    do k = n/2+1,n-1
      s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))* &
      & conjg(x(n-k))
    end do
    write(11,'(2e25.12)') real(s*scale)
  end do
  close(11)

end program tjmcsfft
