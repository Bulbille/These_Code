! $Header: /opt/cvsroot/jmfft/lib/jmgetsetcode.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetsetcode(code,type)

  ! Subroutine qui permet de stocker le dernier code de retour obtenu
  ! Ceci evite de recourir a un common ...

  implicit none

  ! Arguments
  integer, intent(inout) :: code
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer, save :: code_last = 0
  !$OMP THREADPRIVATE(code_last)

  if (type == 's') then
    code_last = code
  else if (type == 'g') then 
    code = code_last
  end if

end subroutine jmgetsetcode
