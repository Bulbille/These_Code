subroutine sinfftmlt(x,work,trigs,ifax,inc,jump,n,m)

  implicit none

  ! Arguments
  integer, intent(in) ::  m, n
  integer :: inc, jump
  real(kind=8), dimension(0:(m-1)*jump+(n-1)*inc) :: x
  real(kind=8), dimension(0:2*m*n+m*(n+2)-1) :: work
  real(kind=8), dimension(0:3*n-1) :: trigs
  integer, dimension(0:18) :: ifax

  ! Variables locales
  integer :: ntrigs, nwork
  real(kind=8) :: pi
  real(kind=8) :: half
  character(len=*), parameter :: nomsp = 'SINFFTMLT'
  integer :: incx, jumpx, isign
  integer :: i, j
  real(kind=8) :: s, t, u

  ! Gestion de pi
  pi = acos(real(-1,kind=8))
  half = 1/real(2,kind=8)

  ! Gestion de table
  ntrigs = 3*n

  ! Gestion de work (dimension pour jmrfftmlt)
  nwork = 2*m*n

  ! On prepare le tableau d'entree
  work(:) = 0
  do j = 0,m-1
    work(nwork+j*(n+2)) = 0
  end do
  if (n > 16 .or. m < 8) then
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = 1,n-1
        s = trigs(2*n+i)
        t = x(j*jump+(i-1)*inc)
        u = x(j*jump+(n-i-1)*inc)
        work(nwork+i+j*(n+2)) = s*(t+u)+half*(t-u)
      end do
    end do
  else
    do i = 1,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        s = trigs(2*n+i)
        t = x(j*jump+(i-1)*inc)
        u = x(j*jump+(n-i-1)*inc)
        work(nwork+i+j*(n+2)) = s*(t+u)+half*(t-u)
      end do
    end do
  end if

  ! On appelle le sous-programme de transformee de Fourier
  isign = -1
  incx = 1
  jumpx = n+2
  call rfftmlt(work(nwork),work,trigs,ifax,incx,jumpx,n,m,isign)

  ! On reconstitue x
  ! Note : Il faut tenir compte des particularites de rfftmlt, qui met un
  !        facteur 1/n par defaut et qui prend une exponentielle negative
  !        Ceci ne s'applique pas bien sur au terme sauvegarde

  ! D'abord les indices impairs
  if (n/2 > 8 .or. m < 8) then
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = 0,n/2-1
        x(j*jump+(2*i+1)*inc) = -n*work(nwork+2*i+3+j*(n+2))
      end do
    end do
  else
    do i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        x(j*jump+(2*i+1)*inc) = -n*work(nwork+2*i+3+j*(n+2))
      end do
    end do
  end if

  ! Ensuite les indices pairs
  ! Cas particulier indice 0
  do j = 0,m-1
    x(j*jump) = n*half*work(nwork+j*(n+2))
  end do
  ! Cas general : recurrence
  if (m > 16 .or. n/2 < 8) then
    do i = 1,(n-1)/2
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        x(j*jump+2*i*inc) = x(j*jump+(2*i-2)*inc) + n*work(nwork+j*(n+2)+2*i)
      end do
    end do
  else
    do j = 0,m-1
      ! Attention : Il faut vectoriser une recurrence
      ! En attendant, voici une version scalaire
      do i = 1,(n-1)/2
        x(j*jump+2*i*inc) = x(j*jump+(2*i-2)*inc) + n*work(nwork+j*(n+2)+2*i)
      end do
    end do
  end if

end subroutine sinfftmlt
