! $Header: /opt/cvsroot/jmfft/lib/jmccfft2d.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine ccfft2d(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: n, m, ldx, ldy
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*ldx*m-1) :: x
  real(kind=8), intent(out), dimension(0:2*ldy*m-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*(n+m)-1) :: table
  real(kind=8), intent(inout), dimension(0:512*max(n,m)-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i, j
  integer :: ioff
  integer :: ntable, nwork
  integer :: nfact, mfact
  integer, dimension(0:99) :: fact
  integer :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  logical :: debut, fin
  character(len=*), parameter :: nomsp = 'CCFFT2D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n) call jmerreur2(nomsp,9,ldx,n)
  if (ldy < n) call jmerreur2(nomsp,14,ldy,n)

  ! Gestion de table
  ntable = 100+2*(n+m)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    table(0:mfact-1) = fact(0:mfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0  ,n)
    call jmtable(table,ntable,100+2*n,m)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    fact(0:mfact-1) = nint(table(0:mfact-1))
  end if

  ! Gestion de work
  !nwork = 4*n*m
  !nwork = 512*max(n,m)
  call jmgetnwork(nwork,512*max(n,m),4*max(n,m))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .true.
  do

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup(m,4*n,nwork,debut,.true.,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! Note : On copie en transposant
    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = jdeb,jfin
        work(j-jdeb+m_temp*i)     =       scale*x(2*i  +j*2*ldx)
        work(j-jdeb+m_temp*(n+i)) = isign*scale*x(2*i+1+j*2*ldx)
      end do
    end do
    ioff = 0

    ! Attention : ioff1 est peut-etre modifie en sortie
    call jmccm1d(m_temp,n,fact,100,0    ,table,ntable,100+0  ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = jdeb,jfin
        y(2*i  +j*2*ldy) = work(ioff+j-jdeb+m_temp*i)
        y(2*i+1+j*2*ldy) = work(ioff+j-jdeb+m_temp*(n+i))
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

  ! On fait les T.F. sur l'autre dimension
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(n,4*m,nwork,debut,.true.,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        work(i-ideb+n_temp*j)     = y(2*i  +j*2*ldy)
        work(i-ideb+n_temp*(m+j)) = y(2*i+1+j*2*ldy)
      end do
    end do
    ioff = 0

    call jmccm1d(n_temp,m,fact,100,nfact,table,ntable,100+2*n,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        y(2*i  +j*2*ldy) =       work(ioff+i-ideb+n_temp*j)
        y(2*i+1+j*2*ldy) = isign*work(ioff+i-ideb+n_temp*(m+j))
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

end subroutine ccfft2d
