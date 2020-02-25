! $Header: /opt/cvsroot/jmfft/lib/jmcfftmlt.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine cfftmlt(xr,xi,work,trigs,ifax,inc,jump,n,m,isign)

  implicit none

  ! Constantes pour les arguments
  integer, parameter :: nfax = 19

  ! Arguments
  integer, intent(in) :: inc, jump, n, m, isign
  real(kind=8), intent(inout), dimension(0:m*n-1) :: xr, xi
  real(kind=8), intent(out), dimension(0:4*n*m-1) :: work
  real(kind=8), intent(in), dimension(0:2*n-1) :: trigs
  integer, intent(in), dimension(0:nfax-1) :: ifax

  ! Variables locales
  integer :: ntrigs, nwork
  integer :: ioff
  integer :: i, j
  character(len=*), parameter :: nomsp = 'CFFTMLT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,1,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)

  ! Gestion de table
  ntrigs = 2*n

  ! Gestion de work
  nwork = 4*n*m

  if (isign == 1) then

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        work(j+m*i)     = xr(jump*j+inc*i)
        work(j+m*(n+i)) = xi(jump*j+inc*i)
      end do
    end do

  else ! isign = -1

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        work(j+m*i)     =  xr(jump*j+inc*i)
        work(j+m*(n+i)) = -xi(jump*j+inc*i)
      end do
    end do

  end if

  ! On appelle le sous-programme principal
  ioff = 0
  call jmccm1d(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork,ioff)

  ! On recopie le tableau de travail dans le tableau de sortie
  if (isign == 1) then

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        xr(jump*j+inc*i) = work(ioff+j+m*i)
        xi(jump*j+inc*i) = work(ioff+j+m*(n+i))
      end do
    end do

  else

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        xr(jump*j+inc*i) =  work(ioff+j+m*i)
        xi(jump*j+inc*i) = -work(ioff+j+m*(n+i))
      end do
    end do

  end if

end subroutine cfftmlt
