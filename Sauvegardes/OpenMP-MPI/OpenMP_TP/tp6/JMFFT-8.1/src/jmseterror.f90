! $Header: /opt/cvsroot/jmfft/lib/jmseterror.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmseterror(arret)

  implicit none

  ! Arguments
  integer, intent(in) :: arret

  ! Variables locales
  integer :: arret2

  arret2 = arret
  call jmgetseterror(arret2,'s')

end subroutine jmseterror
