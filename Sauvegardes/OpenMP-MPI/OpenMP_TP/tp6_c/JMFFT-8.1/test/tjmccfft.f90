! Transformee de Fourier complexe-complexe 1d

program tjmccfft

  implicit none

  integer, parameter :: n = 36
  complex, dimension(0:n-1) :: x, xx
  real, dimension(0:n-1,2) :: rx
  equivalence( x, rx )

  ! Pour stocker les cosinus et les sinus
  ! En fait, pour les routines jm, 100+2*n suffisent
  integer, parameter :: ntable = 100+8*n
  real, dimension(0:ntable-1) :: table

  ! En fait, pour les routines jm, 4*n suffisent
  integer, parameter :: nwork = 8*n
  real, dimension(nwork) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree
  call random_number( rx )
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call ccfft(isign,n,scale,x,x,table,work,isys)
  isign = 1
  print *,'jmccfft ',n,0,isign,scale
  call ccfft(isign,n,scale,x,x,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') x

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  x = xx
  do i = 0,n-1
    s = 0
    do j = 0,n-1
      s = s + cmplx(cos(twopi*i*j/real(n)),isign*sin(twopi*i*j/real(n)))*x(j)
    end do
    write(11,'(2e25.12)') s*scale
  end do

end program tjmccfft
