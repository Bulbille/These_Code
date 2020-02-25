! Transformee de Fourier en sinus multiple

program tjmsinfftmlt

  implicit none

  integer, parameter :: m = 4

  integer, parameter :: n = 8
  real(8), dimension(0:n*m-1) :: x, xx

  integer, parameter :: ntrigs = 3*n
  real(8), dimension(0:ntrigs-1) :: trigs
  integer, parameter :: nfax = 19
  integer, dimension(0:nfax-1) :: ifax

  integer, parameter :: nwork = 2*m*n + m*(n+2)
  real(8), dimension(0:nwork-1) :: work

  integer :: inc, jump

  real(8) :: scale

  integer :: i, j, k
  real(8) :: rij
  real(8) :: pi

  ! On genere x
  call random_number(x)

  ! On force l'imparité
  x(0:m-1) = 0

  ! On sauvegarde en vue du calcul exact
  xx = x

  ! On calcule une fois pour toutes la factorisation
  ! et les tableaux trigonométriques
  call sinfftfax(n,ifax,trigs)

  ! Appel de jmsinfftmlt
  inc = m
  jump = 1
  call sinfftmlt(x,work,trigs,ifax,inc,jump,n,m)

  ! On affiche les resultats
  scale = sqrt( real(2,kind=8) / n )

  ! On imprime la sortie dans le fichier temp1
  open( 11, file = 'temp1', form = 'formatted', status = 'unknown' )
  do j = 0, m-1
    do i = 0, n-1
      write( 11, '(F15.8, I8, I8 )' ) scale*x(j+i*m), i, j
    end do
  end do

  ! On reinitialise x
  x = xx

  ! Gestion de pi
  pi = acos( real(-1,kind=8) )

  ! La valeur exacte
  open( 12, file = 'temp2', form = 'formatted', status = 'unknown' )
  do j = 0, m-1
    do i = 0, n-1
      rij = 0
      do k = 0, n-1
        rij = rij + x(j+k*m) * sin( ( (k+1)*(i+1)*pi ) / n )
      end do
      write( 12, '(F15.8, I8, I8 )' ) scale*rij, i, j
    end do
  end do

end program tjmsinfftmlt
