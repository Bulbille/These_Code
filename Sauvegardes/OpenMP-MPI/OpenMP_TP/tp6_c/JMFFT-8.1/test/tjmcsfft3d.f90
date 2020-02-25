! Transformee de Fourier complexe-reelle multiple

program tjmcsfft3d

  implicit none

  integer, parameter :: m = 1
  integer, parameter :: n = 8
  integer, parameter :: l = 8
  complex, dimension(0:n/2,0:m-1,0:l-1) :: x, xx
  real, dimension(0:n/2,0:m-1,0:l-1,2) :: rx
  equivalence ( x, rx )
  real, dimension(0:n-1+2,0:m-1,0:l-1) :: y

  ! Tableau temporaire pour forcer la symetrie hermitienne
  complex, dimension(0:m-1,0:l-1) :: z

  ! Pour stocker les cosinus et les sinus
  integer, parameter :: ntable = 100+2*(n+m+l)
  real, dimension(0:ntable-1) :: table

  ! Les routines jm ont besoin de 2*2*(n/2+1)*m*l. D'ou tronconnage.
  integer, parameter :: nwork = 512*max(n,m,l)
  real, dimension(0:nwork-1) :: work

  integer :: isign
  real :: scale
  integer :: isys
  integer :: i, j, k, i1, j1, k1
  real :: twopi
  complex :: s

  ! On prepare le tableau d'entree (avec symetrie hermitienne)
  call random_number( rx )
  ! On force la symetrie hermitienne pour i=0 et i=n/2
  do i = 0,n/2,n/2
    do k = 0,l-1
      do j = 0,m-1
        z(j,k) = ( x(i,j,k) + conjg(x(i,mod(m-j,m),mod(l-k,l))) ) / 2
      end do
    end do
    do k = 0,l-1
      do j = 0,m-1
        x(i,j,k) = z(j,k)
      end do
    end do
  end do
  xx = x

  scale = 1.
  isys = 0

  isign = 0
  call csfft3d(isign,n,m,l,scale,x,n/2+1,m,y,n+2,m,table,work,isys)
  isign = 1
  print *,'jmcsfft3d ',n,m,l,isign,scale
  call csfft3d(isign,n,m,l,scale,x,n/2+1,m,y,n+2,m,table,work,isys)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') y(0:n-1,0:m-1,0:l-1)
  close(10)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1))
  ! On reprepare le tableau d'entree (avec symetrie hermitienne)
  x = xx
  ! Et on calcule
  do k = 0,l-1
    do j = 0,m-1
      do i = 0,n-1
        s = 0
        do k1 = 0,l-1
          do j1 = 0,m-1
            do i1 = 0,n/2
              s = s+cmplx(cos(twopi*i*i1/real(n)),isign*sin(twopi*i*i1/real(n))) &
              &    *cmplx(cos(twopi*j*j1/real(m)),isign*sin(twopi*j*j1/real(m))) &
              &    *cmplx(cos(twopi*k*k1/real(l)),isign*sin(twopi*k*k1/real(l))) &
              &    *x(i1,j1,k1)
            end do
            do i1 = n/2+1,n-1
              s = s+cmplx(cos(twopi*i*i1/real(n)),isign*sin(twopi*i*i1/real(n))) &
              &    *cmplx(cos(twopi*j*j1/real(m)),isign*sin(twopi*j*j1/real(m))) &
              &    *cmplx(cos(twopi*k*k1/real(l)),isign*sin(twopi*k*k1/real(l))) &
              &    *conjg(x(n-i1,mod(m-j1,m),mod(l-k1,l)))
            end do
          end do
        end do
        write(11,'(2e25.12)') real(s*scale)
      end do
    end do
  end do
  close(11)

end program tjmcsfft3d
