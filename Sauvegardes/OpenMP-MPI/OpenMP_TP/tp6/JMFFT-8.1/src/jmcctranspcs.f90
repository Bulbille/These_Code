! $Header: /opt/cvsroot/jmfft/lib/jmcctranspcs.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmcctranspcs(m,n,n2,n3,table,ntable,itable,work,nwork,ioff)

  ! Cette subroutine permute le contenu du tableau work de la facon suivante
  ! On considere qu'a l'origine ce tableau est en (m,n3,n2)
  ! On doit transposer le terme d'ordre (k,j,i) en (k,i,j)
  ! On en profite pour faire les multiplications par wij
  ! Le role de n est seulement d'attaquer les bonnes valeurs du tableau table
  ! (il y a un ecart de n entre les cos et les sin, et le stride entre
  !  les cos est de n/(n2*n3)
  ! Note : le sens de n2 et n3 ici n'a rien a voir avec celui de jmccm1d

  implicit none

  ! Arguments
  integer, intent(in) :: m, n
  integer, intent(in) :: n2, n3
  integer, intent(in) :: ntable,itable
  real(kind=8), intent(in),  dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout),  dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: i, j, k
  real(kind=8) :: t, u, c, s
  integer :: ioff1, ioff2
  integer :: is

  ! Gestion des offsets
  ioff1 = ioff
  ioff2 = nwork/2-ioff

  ! Gestion du stride
  is = n/(n2*n3)

  if ( m >= 16 .or. (n2 < 8 .and. n3 < 8) ) then

    do i = 0,n2-1
      do j = 0,n3-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,m-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else if ( n2 >= 16 .or. n3 < 8 ) then

    do j = 0,n3-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n2-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else

    do i = 0,n2-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,n3-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  end if

  ioff = ioff2

end subroutine jmcctranspcs
