! $Header: /opt/cvsroot/jmfft/lib/jmscm1d.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmscm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  implicit none

  ! Arguments
  integer, intent(in) :: m, n
  integer, intent(in) :: nfact, ifact
  integer, intent(inout), dimension(0:nfact-1) :: fact
  integer, intent(in) :: ntable,itable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: ioff1, ioff2
  integer :: i, j
  real(kind=8) :: t, u, v, w
  real(kind=8) :: c, s
  integer :: is, it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue
    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    end if
        
    ! On fait m/2 t.f. complexes de longueur n
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On regenere le resultat
    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue les indices pairs et impairs selon n

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1
    ioff1 = nwork/2 - ioff2

    ! Maintenant, il faut reconstituer la t.f. reelle

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    end if

  end if

  ioff = ioff1

end subroutine jmscm1d
