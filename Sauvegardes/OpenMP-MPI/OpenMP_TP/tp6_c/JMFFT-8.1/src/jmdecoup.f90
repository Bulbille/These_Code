! $Header: /opt/cvsroot/jmfft/lib/jmdecoup.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

! Pour tronconner de facon a tenir dans le nwork disponible

subroutine jmdecoup(n,nr,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

  implicit none

  ! Arguments
  integer, intent(in) :: n, nr, nwork
  logical, intent(in) :: debut, mpair
  integer, intent(out) :: n_temp, ideb, nwork_temp
  integer, intent(inout) :: ifin
  logical, intent(out) :: fin

  ! Variables locales
  character(len=*), parameter :: nomsp = 'JMDECOUP'

  ! n*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n en morceaux pour tenir

  ! Gestion de debut
  if (debut) then
    ideb = 0
  else
    ideb = ifin+1
  end if

  ! Gestion de n_temp et ifin
  n_temp = nwork/nr
  ! Si m impair, on doit eviter que n_temp soit impair (routine cs et sc)
  if (.not.mpair .and. mod(n_temp,2) /= 0) n_temp = n_temp-1
  ifin = min(ideb+n_temp-1,n-1)
  n_temp = ifin-ideb+1
  ! On verifie que n_temp n'est pas nul
  if (n_temp <= 0) then
    call jmerreur3(nomsp,6,n,nr,nwork)
  end if
  nwork_temp = n_temp*nr

  ! Gestion de fin
  if (ifin == n-1) then
    fin = .true.
  else
    fin = .false.
  end if

end subroutine jmdecoup
