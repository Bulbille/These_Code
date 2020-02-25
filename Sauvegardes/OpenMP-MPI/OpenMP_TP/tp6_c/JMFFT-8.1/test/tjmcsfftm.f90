! Transformee de Fourier complexe-reelle multiple

program tjmcsfftm

  implicit none

  integer, parameter :: m = 1
  integer, parameter :: n = 8
  complex, dimension(0:n/2,0:m-1) :: x, xx
  real, dimension(0:n/2,0:m-1,2) :: rx
  equivalence ( x, rx )
  real, dimension(0:n-1,0:m-1) :: y

  ! Pour stocker les cosinus et les sinus
  integer, parameter :: ntable = 100+2*n
  real, dimension(0:ntable-1) :: table

  ! En fait, les routines jm n'ont besoin que de 2*n*m
  integer, parameter :: nwork = (2*n+4)*m
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree en forcant a 0 les termes necessaires
  call random_number( rx )
  x(0,:) = cmplx( real( x(0,:) ), 0 )
  if ( mod( n, 2 ) == 0 ) x(n/2,:) = cmplx( real( x(n/2,:) ), 0 )
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call csfftm(isign,n,m,scale,x,n/2+1,y,n,table,work,isys)
  isign = 1
  print *,'jmcsfftm ',n,m,isign,scale
  call csfftm(isign,n,m,scale,x,n/2+1,y,n,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') y
  close(10)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree
  x = xx
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n-1
      s = 0
      do k = 0,n/2
        s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))*x(k,j)
      end do
      do k = n/2+1,n-1
        s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))* &
        & conjg(x(n-k,j))
      end do
      write(11,'(2e25.12)') real(s*scale)
    end do
  end do
  close(11)

end program tjmcsfftm
