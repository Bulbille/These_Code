! Transformee de Fourier reelle-complexe 1d multiple

program tjmscfftm

  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 16
  real, dimension(0:n-1,0:m-1) :: x, xx
  complex, dimension(0:n/2,0:m-1) :: y

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

  ! On prepare le tableau d'entree
  call random_number( x )
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call scfftm(isign,n,m,scale,x,n,y,n/2+1,table,work,isys)
  isign = 1
  print *,'jmscfftm ',n,m,isign,scale
  call scfftm(isign,n,m,scale,x,n,y,n/2+1,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') y

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree
  x = xx
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n/2
      s = 0
      do k = 0,n-1
        s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n)))*x(k,j)
      end do
      write(11,'(2e25.12)') s*scale
    end do
  end do

end program tjmscfftm
