! Transformee de Fourier complexe-complexe 3d

program tjmccfft3d

  implicit none

  integer, parameter :: n = 9
  integer, parameter :: m = 4
  integer, parameter :: l = 18
  complex, dimension(0:n-1,0:m-1,0:l-1) :: x, xx
  real, dimension(0:n-1,0:m-1,0:l-1,2) :: rx
  equivalence( x, rx )

  integer, parameter :: ntable = 100+2*(n+m+l)
  real, dimension(0:ntable-1) :: table

  ! Les routines jm ont besoin de 4*n*m*l. D'ou tronconnage
  integer, parameter :: nwork = 512*max(n,m,l)
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k
  integer :: i1, j1, k1
  real :: twopi
  complex :: s

  twopi = 2 * acos(real(-1))

  ! On prepare le tableau d'entree
  call random_number(rx)
  xx = x

  isys = 0
  scale = 1./real(n)

  isign = 0
  call ccfft3d(isign,n,m,l,scale,x,n,m,x,n,m,table,work,isys)
  isign = 1
  print *,'jmccfft3d ',n,m,l,isign,scale
  call ccfft3d(isign,n,m,l,scale,x,n,m,x,n,m,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(2e25.12)') x

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  ! On reprepare le tableau d'entree
  x = xx
  ! Et on recalcule
  do k = 0,l-1
    do j = 0,m-1
      do i = 0,n-1
        s = 0
        do k1 = 0,l-1
          do j1 = 0,m-1
            do i1 = 0,n-1
               s = s + cmplx( &
               &               cos(twopi*i*i1/real(n)), &
               &         isign*sin(twopi*i*i1/real(n))) &
               &   *   cmplx( &
               &               cos(twopi*j*j1/real(m)), &
               &         isign*sin(twopi*j*j1/real(m))) &
               &   *   cmplx( &
               &               cos(twopi*k*k1/real(l)), &
               &         isign*sin(twopi*k*k1/real(l))) &
               &   * x(i1,j1,k1)
            end do
          end do
        end do
        write(11,'(2e25.12)') s*scale
      end do
    end do
  end do

end program tjmccfft3d
