subroutine cosfftmlt(x,work,trigs,ifax,inc,jump,n,m)

  implicit none

  ! Arguments
  integer, intent(in) ::  m, n
  integer :: inc, jump
  real(kind=8), dimension(0:(m-1)*jump+n*inc) :: x
  real(kind=8), dimension(0:2*m*n+m+m*(n+2)-1) :: work
  real(kind=8), dimension(0:4*n-1) :: trigs
  integer, dimension(0:18) :: ifax

  ! Variables locales
  integer :: ntrigs, nwork
  real(kind=8) :: pi
  real(kind=8) :: half
  character(len=*), parameter :: nomsp = 'COSFFTMLT'
  integer :: incx, jumpx, isign
  integer :: i, j

  ! Gestion de pi et de half
  pi = acos(real(-1,kind=8))
  half = 1/real(2,kind=8)

  ! Gestion de table
  ntrigs = 4*n

  ! Gestion de work (dimension pour jmrfftmlt)
  nwork = 2*m*n

  ! On calcule par anticipation le premier terme de rang impair
  do j = 0, m-1
    work(nwork+j) = half * ( x(j*jump) - x(j*jump+n*inc) )
  end do
  if (m > 16 .or. n < 8) then
    do i = 1, n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0, m-1
        work(nwork+j) = work(nwork+j) + trigs(3*n+i) * x(j*jump+i*inc)
      end do
    end do
  else
    do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = 1, n-1
        work(nwork+j) = work(nwork+j) + trigs(3*n+i) * x(j*jump+i*inc)
      end do
    end do
  end if

  ! On prepare le tableau d'entree
  if (n > 16 .or. m < 8) then
    do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = 0, n-1
        work(nwork+m+i+j*(n+2)) = &
           half * (x(j*jump+i*inc) + x(j*jump+(n-i)*inc) ) &
          - trigs(2*n+i) * ( x(j*jump+i*inc) - x(j*jump+(n-i)*inc) )
      end do
    end do
  else
    do i = 0, n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0, m-1
        work(nwork+m+i+j*(n+2)) = &
           half * (x(j*jump+i*inc) + x(j*jump+(n-i)*inc) ) &
          - trigs(2*n+i) * ( x(j*jump+i*inc) - x(j*jump+(n-i)*inc) )
      end do
    end do
  end if

  ! On appelle le sous-programme de transformee de Fourier
  isign = -1
  incx = 1
  jumpx = n+2
  call rfftmlt(work(nwork+m),work,trigs,ifax,incx,jumpx,n,m,isign)

  ! On reconstitue x
  ! Note : Il faut tenir compte des particularites de rfftmlt, qui met un
  !        facteur 1/n par defaut et qui prend une exponentielle negative
  !        Ceci ne s'applique pas bien sur au terme sauvegarde

  ! Traitement des indices pairs
  if (n/2 > 16 .or. m < 8) then
    do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = 0, n, 2
        x(j*jump+i*inc) = n*work(nwork+m+i+j*(n+2))
      end do
    end do
  else
    do i = 0, n, 2
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0, m-1
        x(j*jump+i*inc) = n*work(nwork+m+i+j*(n+2))
      end do
    end do
  end if

  ! Traitement des indices impairs
  do j = 0, m-1
    x(j*jump+inc) = work(nwork+j)
  end do
  if (m > 16 .or. n/2 < 8) then
    do i = 3, n, 2
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0, m-1
        x(j*jump+i*inc) = x(j*jump+(i-2)*inc) - n*work(nwork+m+i+j*(n+2))
      end do
    end do
  else
    do j = 0, m-1
      ! Attention : Il faut vectoriser une recurrence
      ! En attendant, voici une version scalaire
      do i = 3, n, 2
        x(j*jump+i*inc) = x(j*jump+(i-2)*inc) - n*work(nwork+m+i+j*(n+2))
      end do
    end do
  end if

end subroutine cosfftmlt
