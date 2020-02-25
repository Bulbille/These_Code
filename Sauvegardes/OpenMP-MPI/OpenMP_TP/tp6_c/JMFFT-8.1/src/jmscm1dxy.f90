! $Header: /opt/cvsroot/jmfft/lib/jmscm1dxy.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

! Variante de jmscm1d ou on fournit x en entree et en sortie
subroutine jmscm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  implicit none

  ! Arguments
  integer, intent(in) :: m, n
  integer, intent(in) :: nfact, ifact
  integer, intent(inout), dimension(0:nfact-1) :: fact
  integer, intent(in) :: ntable,itable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(in) :: dimx, debx, incx, jumpx
  integer, intent(in) :: dimy, deby, incy, jumpy
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in),  dimension(0:dimx-1) :: x
  real(kind=8), intent(out), dimension(0:dimy-1) :: y
  integer, intent(in) :: isign

  ! Variables locales
  integer :: i, j
  real(kind=8) :: t, u, v, w
  real(kind=8) :: c, s
  integer :: is, it
  integer :: ioff

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
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        end do
      end do

    end if

    ! On fait m/2 t.f. complexes de longueur n
    ioff = 0
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On regenere le resultat
    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2
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
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2
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
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

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
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2 + c*(u+w)/2 - s*(v-t)/2)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2 + c*(v-t)/2 + s*(u+w)/2)
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
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2 + c*(u+w)/2 - s*(v-t)/2)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2 + c*(v-t)/2 + s*(u+w)/2)
        end do
      end do

    end if

  end if

end subroutine jmscm1dxy
