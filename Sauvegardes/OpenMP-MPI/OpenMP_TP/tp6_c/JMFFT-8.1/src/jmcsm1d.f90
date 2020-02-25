! $Header: /opt/cvsroot/jmfft/lib/jmcsm1d.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmcsm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  implicit none

  ! Arguments
  integer, intent(in) :: m, n
  integer, intent(in) :: nfact, ifact
  integer, intent(inout), dimension(0:nfact-1) :: fact
  integer, intent(in) :: ntable,itable
  real(kind=8), intent(inout), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: ioff1, ioff2
  integer :: i, j
  real(kind=8) :: t, u, v, w, tt, uu, vv, ww
  real(kind=8) :: c, s
  integer :: it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue

    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
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
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    end if

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On reconstitue

    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2*(uu+vv)
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

    ! On reconstitue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        end do
      end do

    end if

  end if

  ioff = ioff1

end subroutine jmcsm1d
