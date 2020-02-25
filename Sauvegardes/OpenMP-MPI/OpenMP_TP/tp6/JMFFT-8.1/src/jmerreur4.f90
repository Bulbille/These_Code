! $Header: /opt/cvsroot/jmfft/lib/jmerreur4.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmerreur4(nomsp,code,var1,var2,var3,var4)

  implicit none

  ! Arguments
  character(len=*), intent(in) :: nomsp
  integer, intent(in) :: code
  integer, intent(in) :: var1, var2, var3,var4

  ! Variables locales
  integer :: arret
  character(len=80) :: message

  call jmgeterror(arret)
  if (arret == 1) then
    call jmgetmessage(code,message)
    print *,'JMFFT Erreur dans ',trim(nomsp),' : ',trim(message), &
    & ' (',var1,var2,var3,var4,')'
    stop 1
  else
    call jmsetcode(code)
  end if

end subroutine jmerreur4
