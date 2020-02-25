! $Header: /opt/cvsroot/jmfft/lib/jmgetnwork.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetnwork(nwork,nwork_def,nwork_min)

  ! On recupere la valeur de nwork si elle a ete augmentee par l'utilisateur
  ! Sinon on prend la valeur par defaut
  ! Il s'agit du nwork des routines 2d et 3d

  implicit none

  ! Arguments
  integer, intent(out) :: nwork
  integer, intent(in)  :: nwork_def, nwork_min

  ! Variables locales
  integer :: nwork_loc
  character(len=*), parameter :: nomsp = 'JMGETNWORK'

  call jmgetsetnwork(nwork_loc,'g')

  ! Valeur par defaut
  if (nwork_loc == -1) then
    nwork = nwork_def
  ! Valeur invalide (trop petite)
  else if (nwork_loc < nwork_min) then
    call jmerreur2(nomsp,5,nwork_loc,nwork_min)
  ! Valeur correcte
  else
    nwork = nwork_loc
  end if

end subroutine jmgetnwork
