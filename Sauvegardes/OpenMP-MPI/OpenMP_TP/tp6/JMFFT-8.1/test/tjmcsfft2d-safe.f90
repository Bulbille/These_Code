! Transformee de Fourier complexe-reelle multiple

program tjmcsfft2d

  implicit none

  integer, parameter :: m = 1
  integer, parameter :: n = 8
  complex, dimension(0:n/2,0:m-1) :: x, xx
  real, dimension(0:n/2,0:m-1,2) :: rx
  equivalence ( x, rx )
  real, dimension(0:n-1+2,0:m-1) :: y

  ! Pour stocker les cosinus et les sinus
  integer, parameter :: ntable = 100+2*(n+m)
  real, dimension(0:ntable-1) :: table

  ! Les routines jm ont besoin de 2*2*(n/2+1)*m. D'ou tronconnage.
  integer, parameter :: nwork = 512*max(n,m)
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k, l
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree sans forcer a 0 les termes necessaires
  call random_number( rx )
  ! On force la symetrie hermitienne pour i=0 et i=n/2
  do j = m/2+1,m-1
    x(0,j) = conjg(x(0,m-j))
  end do
  !JMT- x(0,0) = cmplx(real(x(0,0)),0)
  !JMT- if (mod(m,2) == 0) x(0,m/2) = cmplx(real(x(0,m/2)),0)
  if (mod(n,2) == 0) then
    do j = m/2+1,m-1
      x(n/2,j) = conjg(x(n/2,m-j))
    end do
    !JMT- x(n/2,0) = cmplx(real(x(n/2,0)),0)
    !JMT- if (mod(m,2) == 0) x(n/2,m/2) = cmplx(real(x(n/2,m/2)),0)
  end if
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call csfft2d(isign,n,m,scale,x,n/2+1,y,n+2,table,work,isys)
  isign = 1
  print *,'jmcsfft2d ',n,m,isign,scale
  call csfft2d(isign,n,m,scale,x,n/2+1,y,n+2,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') y(0:n-1,0:m-1)
  close(10)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree, mais on force a 0 les termes necessaires
  x = xx
  x(0,0) = cmplx(real(x(0,0)),0)
  if (mod(m,2) == 0) x(0,m/2) = cmplx(real(x(0,m/2)),0)
  if (mod(n,2) == 0) then
    x(n/2,0) = cmplx(real(x(n/2,0)),0)
    if (mod(m,2) == 0) x(n/2,m/2) = cmplx(real(x(n/2,m/2)),0)
  end if
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n-1
      s = 0
      do l = 0,m-1
        do k = 0,n/2
          s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n))) &
          &    *cmplx(cos(twopi*j*l/real(m)),isign*sin(twopi*j*l/real(m))) &
          &    *x(k,l)
        end do
        do k = n/2+1,n-1
          s = s+cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n))) &
          &    *cmplx(cos(twopi*j*l/real(m)),isign*sin(twopi*j*l/real(m))) &
          &    *conjg(x(n-k,mod(m-l,m)))
        end do
      end do
      write(11,'(2e25.12)') real(s*scale)
    end do
  end do
  close(11)

end program tjmcsfft2d
