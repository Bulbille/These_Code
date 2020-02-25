! Transformee de Fourier reelle-complexe 1d

program tjmscfft

  implicit none

  integer, parameter :: n = 12
  real, dimension(0:n-1) :: x, xx
  complex, dimension(0:n/2) :: y

  ! Pour stocker les cosinus et les sinus
  integer, parameter :: ntable = 100+2*n
  real, dimension(0:ntable-1) :: table

  ! En fait, les routines jm n'ont besoin que de 2*n
  integer, parameter :: nwork = 4+4*n
  real, dimension(nwork) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree
  call random_number( x )
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call scfft(isign,n,scale,x,y,table,work,isys)
  isign = 1
  print *,'jmscfft ',n,isign,scale
  call scfft(isign,n,scale,x,y,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') y

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree
  x = xx
  ! Et on calcule
  do i = 0,n/2
    s = 0
    do j = 0,n-1
      s = s + cmplx(cos(twopi*i*j/real(n)),isign*sin(twopi*i*j/real(n)))*x(j)
    end do
    write(11,'(2e25.12)') s*scale
  end do

end program tjmscfft
