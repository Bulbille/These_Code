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
  integer :: j2, k2

  ! On prepare deux tableaux d'entree
  ! x en forcant la symetrie hermitienne, sauf termes a parties reelles nulles
  ! xx en forcant la symetrie hermitienne, avec termes a parties reelles nulles
  call random_number( rx )
  xx = x
  ! On force la symetrie hermitienne pour i=0 et i=n/2
  do k = 0,l-1
    k2 = l-k
    if (k == 0) k2 = 0
    do j = 0,m-1
      j2 = m-j
      if ( j== 0) j2 = 0
      if ( ( j > j2 ) .or. ( j == j2 .and. k > k2 ) ) then
        x(0,j,k) = conjg( x( 0,j2,k2 ) )
        xx(0,j,k) = conjg( xx( 0,j2,k2 ) )
      end if
      ! Note : le test suivant engloble les cas 0 et moitie si parite
      !JMT- if ( (j==j2) .and. (k==k2) ) x(0,j,k) = real( x(0,j,k) )
      if ( (j==j2) .and. (k==k2) ) xx(0,j,k) = real( xx(0,j,k), 8 )
    end do
  end do
  if ( mod(n,2) == 0 ) then
    do k = 0,l-1
      k2 = l-k
      if (k == 0) k2 = 0
      do j = 0,m-1
        j2 = m-j
        if ( j== 0) j2 = 0
        if ( ( j > j2 ) .or. ( j == j2 .and. k > k2 ) ) then
          x(n/2,j,k) = conjg( x( n/2,j2,k2 ) )
          xx(n/2,j,k) = conjg( xx( n/2,j2,k2 ) )
        end if
        ! Note : le test suivant engloble les cas 0 et moitie si parite
        !JMT- if ( (j==j2) .and. (k==k2) ) x(n/2,j,k) = real( x(n/2,j,k) )
        if ( (j==j2) .and. (k==k2) ) xx(n/2,j,k) = real( xx(n/2,j,k), 8 )
      end do
    end do
  end if

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
  x = xx
  open(11,file='temp2',status='unknown',form='formatted')
  twopi = 2 * acos(real(-1,8))
  ! Et on calcule
  do k = 0,l-1
    do j = 0,m-1
      do i = 0,n-1
        s = 0
        do k1 = 0,l-1
          do j1 = 0,m-1
            do i1 = 0,n/2
              s = s+cmplx(cos(twopi*i*i1/real(n,8)),isign*sin(twopi*i*i1/real(n,8))) &
              &    *cmplx(cos(twopi*j*j1/real(m,8)),isign*sin(twopi*j*j1/real(m,8))) &
              &    *cmplx(cos(twopi*k*k1/real(l,8)),isign*sin(twopi*k*k1/real(l,8))) &
              &    *x(i1,j1,k1)
            end do
            do i1 = n/2+1,n-1
              s = s+cmplx(cos(twopi*i*i1/real(n,8)),isign*sin(twopi*i*i1/real(n,8))) &
              &    *cmplx(cos(twopi*j*j1/real(m,8)),isign*sin(twopi*j*j1/real(m,8))) &
              &    *cmplx(cos(twopi*k*k1/real(l,8)),isign*sin(twopi*k*k1/real(l,8))) &
              &    *conjg(x(n-i1,mod(m-j1,m),mod(l-k1,l)))
            end do
          end do
        end do
        write(11,'(2e25.12)') real(s*scale,8)
      end do
    end do
  end do
  close(11)

end program tjmcsfft3d
