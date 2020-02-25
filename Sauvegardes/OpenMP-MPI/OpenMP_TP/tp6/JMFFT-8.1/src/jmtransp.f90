! $Header: /opt/cvsroot/jmfft/lib/jmtransp.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmtransp(n,m,l,work,nwork,ioff)

  implicit none

  ! Arguments
  integer, intent(in) :: n, m, l
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: ioff1, ioff2
  integer :: ij, k

  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! On transpose (nm)(l) en (l)(nm) en distinguant les parties reelles et im.
  if (m*n >= 16 .or. l < 8) then

    do k = 0,l-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do ij = 0,m*n-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      end do
    end do

  else

    do ij = 0,m*n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do k = 0,l-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      end do
    end do

  end if

  ioff = ioff2

end subroutine jmtransp
