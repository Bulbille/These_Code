! $Header: /opt/cvsroot/jmfft/lib/jmgetseterreur.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetseterreur(arret,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  implicit none

  ! Arguments
  logical, intent(inout) :: arret
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  logical, save :: arret_last = .true.
  !$OMP THREADPRIVATE(arret_last)

  if (type == 's') then
    arret_last = arret
  else if (type == 'g') then 
    arret = arret_last
  end if

end subroutine jmgetseterreur
