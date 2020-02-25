! $Header: /opt/cvsroot/jmfft/lib/jmccm1d.f90,v 1.3 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
  implicit none

  ! Arguments
  integer, intent(in) :: m, n
  integer, intent(in) :: nfact, ifact
  integer, intent(in), dimension(0:nfact-1) :: fact
  integer, intent(in) :: ntable,itable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: nterms
  integer :: np, pp, lastnp, premier
  integer :: nprod, nprod1, nprod2
  integer :: n2, p2, n3, p3, n5, p5
  integer :: i

  logical, save :: copyright = .false.
  !$OMP THREADPRIVATE(copyright)

!$OMP SINGLE
  if (.not.copyright) then
    copyright = .true.
    print *,' '
    print *,'************************************************************'
    print *,'*       Portable Fourier transforms by JMFFTLIB            *'
    print *,'* Author : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr   *'
    print *,'************************************************************'
    print *,' '
  end if
!$OMP END SINGLE

  ! On recupere les facteurs
  nterms = fact(ifact)
  n2 = fact(ifact+1)
  p2 = fact(ifact+2)
  n3 = fact(ifact+3)
  p3 = fact(ifact+4)
  n5 = fact(ifact+5)
  p5 = fact(ifact+6)
  nprod = n2*n3*n5
  do i = 7,nterms-1,2
    nprod = nprod*fact(ifact+i+1)**fact(ifact+i)
  end do

  ! On fait n3*n5 T.F. de n2 (qui est en puissances de 2)
  if (n2 /= 1) then
    call jmccm1d2(p2,n2,m*(nprod/n2),table,ntable,itable,n,n/n2,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n2 /= 1 .and. nprod /= n2) then
    call jmcctranspcs(m,n,n2,nprod/n2,table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n5*n2 T.F. de n3 (en puissances de 3)
  if (n3 /= 1) then
    call jmccm1d3(p3,n3,m*(nprod/n3),table,ntable,itable,n,n/n3,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n3 /= 1 .and. nprod /= n3) then
    call jmcctranspcs(m*n2,n,n3,nprod/(n2*n3), &
    & table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n2*n3 T.F. de n5 (en puissances de 5)
  if (n5 /= 1) then
    call jmccm1d5(p5,n5,m*(nprod/n5),table,ntable,itable,n,n/n5,work,nwork,ioff)
  end if

  ! On transpose s'il y a lieu (si on a fait quelque chose et s'il reste des
  ! termes a traiter
  if (n5 /= 1 .and. nprod /= n5 .and. nterms > 7) then
    call jmcctranspcs(m*n2*n3,n,n5,nprod/(n2*n3*n5), &
    & table,ntable,itable,work,nwork,ioff)
  end if
  nprod1 = m*n2*n3
  nprod2 = n2*n3*n5
  lastnp = n5

  ! On passe aux nombres premiers autres que 2, 3 et 5
  do i = 7,nterms-1,2

    pp = fact(ifact+i)
    premier = fact(ifact+i+1)
    np = premier**pp

    call jmccm1dp(premier,pp,m*(nprod/np), &
    & table,ntable,itable,n,n/np,work,nwork,ioff)

    nprod1 = nprod1 * lastnp
    nprod2 = nprod2 * np
    if (np /= 1 .and. nprod /= np .and. nterms > i+1) then
      call jmcctranspcs(nprod1,n,np,nprod/nprod2, &
      & table,ntable,itable,work,nwork,ioff)
    end if
    lastnp = np

  end do

end subroutine jmccm1d
