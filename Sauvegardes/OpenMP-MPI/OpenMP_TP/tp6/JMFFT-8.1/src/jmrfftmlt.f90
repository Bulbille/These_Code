! $Header: /opt/cvsroot/jmfft/lib/jmrfftmlt.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine rfftmlt(x,work,trigs,ifax,incx,jumpx,n,m,isign)

  implicit none

  ! Constantes pour les arguments
  integer, parameter :: nfax = 19

  ! Arguments
  integer, intent(in) :: incx, jumpx, n, m, isign
  real(kind=8), intent(inout), dimension(0:m*(n+2)-1) :: x
  real(kind=8), intent(out), dimension(0:2*n*m-1) :: work
  real(kind=8), intent(in), dimension(0:2*n-1) :: trigs
  integer, intent(in), dimension(0:nfax-1) :: ifax

  ! Variables locales
  integer :: ntrigs, nwork
  real(kind=8) :: scale
  integer :: dimx, debx
  integer :: signe
  real(kind=8) :: scale_temp
  character(len=*), parameter :: nomsp = 'RFFTMLT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,1,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (mod(n,2) /= 0 .and. mod(m,2) /= 0) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de trigs
  ntrigs = 2*n

  ! Gestion de work
  nwork = 2*n*m

  if (isign == -1) then

    ! On appelle le sous-programme principal
    scale = real(1,kind=8)/real(n,kind=8)
    dimx = m*(n+2)
    debx = 0
    signe = -1
    call jmscm1dxy(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork, &
    & x,dimx,debx,incx,jumpx,x,dimx,debx,incx,jumpx,signe,scale)

  else

    ! On appelle le sous-programme principal
    dimx = m*(n+2)
    debx = 0
    signe = 1
    scale_temp = real(1,kind=8)
    call jmcsm1dxy(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork, &
    & x,dimx,debx,incx,jumpx,x,dimx,debx,incx,jumpx,signe,scale_temp)

  end if

end subroutine rfftmlt
