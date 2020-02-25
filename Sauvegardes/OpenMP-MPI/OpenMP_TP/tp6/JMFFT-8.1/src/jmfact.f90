! $Header: /opt/cvsroot/jmfft/lib/jmfact.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmfact(n,fact,nfact,ideb,ifin)

  implicit none

  ! Arguments
  integer, intent(in) :: n, nfact, ideb
  integer, intent(out) :: ifin
  integer, intent(inout), dimension(0:nfact-1) :: fact

  ! Variables locales
  integer :: m
  integer :: n2, p2, n3, p3, n5, p5
  character(len=*), parameter :: nomsp = 'JMFACT'
  ! Nombres premiers
  integer, parameter :: npremiers = 7
  integer, dimension(0:npremiers-1) :: premiers = (/7,11,13,17,19,23,29/)
  integer :: ip, premier, pp, np

  m = n

  ! Etude des puissances de deux
  p2 = 0
  n2 = 1
  do
    if (mod(m,2) == 0) then
      p2 = p2+1
      n2 = n2*2
      m  = m/2
    else
      exit
    end if
  end do
  ifin = ideb+3
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+1) = n2
  fact(ideb+2) = p2

  ! Etude des puissances de trois
  p3 = 0
  n3 = 1
  do
    if (mod(m,3) == 0) then
      p3 = p3+1
      n3 = n3*3
      m  = m/3
    else
      exit
    end if
  end do
  ifin = ifin+2
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+3) = n3
  fact(ideb+4) = p3

  ! Etude des puissances de cinq
  p5 = 0
  n5 = 1
  do
    if (mod(m,5) == 0) then
      p5 = p5+1
      n5 = n5*5
      m  = m/5
    else
      exit
    end if
  end do
  ifin = ifin+2
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+5) = n5
  fact(ideb+6) = p5

  ! On met a jour le nombre de termes
  fact(ideb) = 7

  ! Si on a fini
  if (n2*n3*n5 == n) return

  ! Il reste maintenant des facteurs premiers bizarres
  ! On va boucler tant qu'on n'a pas fini ou tant qu'on n'a pas epuise la liste

  do ip = 0,npremiers-1

    premier = premiers(ip)

    pp = 0
    np = 1
    do
      if (mod(m,premier) == 0) then
        pp = pp+1
        np = np*premier
        m  = m/premier
      else
        exit
      end if
    end do
    ifin = ifin+2
    if (ifin > nfact) &
    & call jmerreur2(nomsp,7,nfact,ifin)
    fact(ifin-2) = pp
    fact(ifin-1) = premier
    fact(ideb) = fact(ideb) + 2

    ! Si le nombre est completement factorise, inutile de continuer
    if (m == 1) exit

  end do

  ! On regarde si la factorisation est terminee
  if (m == 1) then
    return
  else
    call jmerreur1(nomsp,3,n)
  end if
  
end subroutine jmfact
