! $Header: /opt/cvsroot/jmfft/lib/jmgetsetnwork.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetsetnwork(nwork,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  implicit none

  ! Arguments
  integer, intent(inout) :: nwork
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer, save :: nwork_last = -1
  !$OMP THREADPRIVATE(nwork_last)

  if (type == 's') then
    nwork_last = nwork
  else if (type == 'g') then 
    nwork = nwork_last
  end if

end subroutine jmgetsetnwork
