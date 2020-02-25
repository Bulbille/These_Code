! $Header: /opt/cvsroot/jmfft/lib/jmgeterreur.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgeterreur(arret)

  implicit none

  ! Arguments
  logical, intent(out) :: arret

  ! Variables locales

  call jmgetseterreur(arret,'g')

end subroutine jmgeterreur
