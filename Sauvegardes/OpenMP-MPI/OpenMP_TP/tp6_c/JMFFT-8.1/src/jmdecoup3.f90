! $Header: /opt/cvsroot/jmfft/lib/jmdecoup3.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

! Pour tronconner en dimension 3 de facon a tenir dans le nwork disponible

subroutine jmdecoup3(n,m,nmr,nwork,debut,lpair,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)

  implicit none

  ! Arguments
  integer, intent(in) :: n, m, nmr, nwork
  logical, intent(in) :: debut, lpair
  integer, intent(out) :: nmtemp, nwork_temp
  integer, intent(out)   :: ideb, jdeb
  integer, intent(inout) :: ifin, jfin
  logical, intent(out) :: fini

  ! Variables locales
  integer :: ijdeb, ijfin
  character(len=*), parameter :: nomsp = 'JMDECOUP3'

  ! n*m*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n et m en morceaux pour tenir

  ! Gestion de debut
  if (debut) then
    ideb = 0
    jdeb = 0
  else
    if (ifin < n-1) then
      ideb = ifin+1
      jdeb = jfin
    else
      ideb = 0
      jdeb = jfin+1
    end if
  end if

  ! Gestion de nmtemp
  nmtemp = nwork/nmr
  ! Si l impair, on doit eviter que nmtemp soit impair (routine cs et sc)
  if (.not.lpair .and. mod(nmtemp,2) /= 0) nmtemp = nmtemp-1
  ! Pour simplifier, on passe par des indices 2d
  ijdeb = ideb+jdeb*n
  ijfin = min(ijdeb+nmtemp-1,n*m-1)
  nmtemp = ijfin-ijdeb+1
  ! On verifie que nmtemp n'est pas nul
  if (nmtemp <= 0) then
    call jmerreur4(nomsp,6,n,m,nmr,nwork)
  end if
  nwork_temp = nmtemp*nmr

  ! On deduit ifin et jfin de ijfin
  jfin = ijfin/n
  ifin = ijfin-n*jfin

  ! Gestion de fin
  if (ifin == n-1 .and. jfin == m-1) then
    fini = .true.
  else
    fini = .false.
  end if

end subroutine jmdecoup3
