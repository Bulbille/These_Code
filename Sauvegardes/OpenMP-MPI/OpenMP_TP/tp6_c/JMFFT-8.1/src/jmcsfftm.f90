! $Header: /opt/cvsroot/jmfft/lib/jmcsfftm.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine csfftm(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: m, n, ldx, ldy
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*ldx-1,0:m-1) :: x
  real(kind=8), intent(out), dimension(0:ldy-1,0:m-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*n-1) :: table
  real(kind=8), intent(inout), dimension(0:2*n*m-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i
  integer :: ntable, nwork
  integer :: nfact
  integer, dimension(0:99) :: fact
  integer :: dimx, debx, incx, jumpx
  integer :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'CSFFTM'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n/2+1) call jmerreur2(nomsp,10,ldx,n/2+1)
  if (ldy < n    ) call jmerreur2(nomsp,14,ldy,n    )
  if (mod(n,2) /= 0 .and. mod(m,2) /= 0) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n*m

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
  dimx = 2*ldx*m ; debx = 0 ; incx = 1 ; jumpx = 2*ldx
  dimy = ldy*m   ; deby = 0 ; incy = 1 ; jumpy = ldy
  call jmcsm1dxy(m,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

end subroutine csfftm
