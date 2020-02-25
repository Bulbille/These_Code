! $Header: /opt/cvsroot/jmfft/lib/jmsetcode.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmsetcode(code)

  implicit none

  ! Arguments
  integer, intent(in) :: code

  ! Variables locales
  integer :: errcode

  errcode = code
  call jmgetsetcode(errcode,'s')

end subroutine jmsetcode
