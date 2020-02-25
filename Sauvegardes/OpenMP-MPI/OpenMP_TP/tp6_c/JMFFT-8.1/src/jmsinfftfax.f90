! $Header: /opt/cvsroot/jmfft/lib/jmsinfftfax.f90,v 1.3 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine sinfftfax(n,ifax,trigs)

  implicit none

  ! Constantes pour les arguments
  integer, parameter :: nfax = 19

  ! Arguments
  integer, intent(in) :: n
  real(kind=8), dimension(0:3*n-1), intent(out) :: trigs
  integer, dimension(0:nfax-1), intent(out) :: ifax

  ! Variables locales
  integer :: ntrigs 
  integer :: ifin
  character(len=*), parameter :: nomsp = 'SINFFTFAX'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ntrigs = 3*n

  ! Factorisation de n dans ifax
  call jmfact(n,ifax,nfax,0,ifin)

  ! Preparation des tables
  call jmtableSin(trigs,ntrigs,0,n)

end subroutine sinfftfax
