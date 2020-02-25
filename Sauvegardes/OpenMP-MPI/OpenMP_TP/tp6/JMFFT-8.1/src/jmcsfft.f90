! $Header: /opt/cvsroot/jmfft/lib/jmcsfft.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine csfft(isign,n,scale,x,y,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: n
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*(n/2)+1) :: x
  real(kind=8), intent(out),  dimension(0:n-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*n-1) :: table
  real(kind=8), intent(inout), dimension(0:2*n-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i
  integer :: ntable, nwork
  integer :: nfact
  integer, dimension(0:99) :: fact
  integer :: dimx, debx, incx, jumpx
  integer :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'CSFFT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (mod(n,2) /= 0) call jmerreur1(nomsp,24,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n

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

  ! On appelle le sous-programme principal
  dimx = 2*(n/2)+2 ; debx = 0 ; incx = 1 ; jumpx = 0
  dimy = n         ; deby = 0 ; incy = 1 ; jumpy = 0
  call jmcsm1dxy(1,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpy,y,dimy,deby,incy,jumpy,isign,scale)

end subroutine csfft
