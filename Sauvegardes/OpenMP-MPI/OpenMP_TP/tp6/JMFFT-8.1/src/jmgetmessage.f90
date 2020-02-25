! $Header: /opt/cvsroot/jmfft/lib/jmgetmessage.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetmessage(code,message)

  implicit none

  ! Arguments
  integer, intent(in) :: code
  character(len=*), intent(out) :: message

  ! Variables locales
  integer, parameter :: mm = 26
  character(len=34), dimension(0:mm-1) :: messages = (/ &
  & "Pas d'erreur                     ",                &
  & "Isign doit etre egal a -1 ou 1   ",                &
  & "Isign doit etre egal a 0, -1 ou 1",                &
  & "Nombres premiers trop grands     ",                &
  & "Nwork negatif ou nul             ",                &
  & "Nwork trop petit                 ",                &
  & "Tronconnage impossible           ",                &
  & "Trop de facteurs premiers        ",                &
  & "l doit etre >= 1                 ",                &
  & "ldx doit etre >= n               ",                &
  & "ldx doit etre >= n/2+1           ",                &
  & "ldx1 doit etre >= n              ",                &
  & "ldx1 doit etre >= n/2+1          ",                &
  & "ldx2 doit etre >= m              ",                &
  & "ldy doit etre >= n               ",                &
  & "ldy doit etre >= n+2             ",                &
  & "ldy doit etre >= n/2+1           ",                &
  & "ldy1 doit etre >= n              ",                &
  & "ldy1 doit etre >= n+2            ",                &
  & "ldy1 doit etre >= n/2+1          ",                &
  & "ldy2 doit etre >= m              ",                &
  & "m doit etre >= 1                 ",                &
  & "m ou n doit etre pair            ",                &
  & "n doit etre >= 1                 ",                &
  & "n doit etre pair                 ",                &
  & "n ou m ou l doit etre pair       "                 &
  & /)

  if (code < 0 .or. code >= mm) then
    print *,'JMFFT GETMESSAGE Code invalide : ',code
  end if

  message = messages(code)

end subroutine jmgetmessage
