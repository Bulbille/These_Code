! $Header: /opt/cvsroot/jmfft/lib/jmseterreur.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmseterreur(arret)

  implicit none

  ! Arguments
  logical, intent(in) :: arret

  ! Variables locales
  logical :: arret2

  arret2 = arret
  call jmgetseterreur(arret2,'s')

end subroutine jmseterreur
