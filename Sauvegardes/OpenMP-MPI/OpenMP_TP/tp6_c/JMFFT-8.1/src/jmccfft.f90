! $Header: /opt/cvsroot/jmfft/lib/jmccfft.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine ccfft(isign,n,scale,x,y,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: n
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*n-1) :: x
  real(kind=8), intent(out), dimension(0:2*n-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*n-1) :: table
  real(kind=8), intent(inout), dimension(0:4*n-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i
  integer :: ioff
  integer :: ntable, nwork
  integer :: nfact
  integer, dimension(0:99) :: fact
  character(len=*), parameter :: nomsp = 'CCFFT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 4*n

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  do i = 0,n-1
    work(i)   =       scale* x(2*i)
    work(n+i) = isign*scale* x(2*i+1)
  end do
  ioff = 0

  ! On appelle le sous-programme principal
  call jmccm1d(1,n,fact,100,0,table,ntable,100+0,work,nwork,ioff)

  ! On recopie dans le tableau d'arrivee
  do i = 0,n-1
    y(2*i)   =         work(ioff  +i)
    y(2*i+1) = isign * work(ioff+n+i)
  end do

end subroutine ccfft
