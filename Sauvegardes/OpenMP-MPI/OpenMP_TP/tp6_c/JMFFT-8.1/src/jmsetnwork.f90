! $Header: /opt/cvsroot/jmfft/lib/jmsetnwork.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmsetnwork(nwork)

  ! Subroutine appelee par l'utilisateur pour augmenter le nwork
  ! des routines 2d et 3d

  implicit none

  ! Arguments
  integer, intent(in) :: nwork

  ! Variables locales
  character(len=*), parameter :: nomsp = 'JMSETNWORK'
  integer :: nwork2

  if (nwork <= 0) then
    call jmerreur1(nomsp,4,nwork)
  end if

  nwork2 = nwork
  call jmgetsetnwork(nwork2,'s')

end subroutine jmsetnwork
