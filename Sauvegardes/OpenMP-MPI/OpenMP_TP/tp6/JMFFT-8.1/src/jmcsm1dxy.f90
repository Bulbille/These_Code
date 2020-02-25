! $Header: /opt/cvsroot/jmfft/lib/jmcsm1dxy.f90,v 1.3 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

! Variante de jmcsm1d ou on fournit x en entree et en sortie
subroutine jmcsm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

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
  real(kind=8), intent(in),  dimension(0:dimx-1) :: x
  real(kind=8), intent(out), dimension(0:dimy-1) :: y
  integer, intent(in) :: isign
  real(kind=8), intent(in) :: scale

  ! Variables locales
  integer :: i, j
  real(kind=8) :: t, u, v, w, tt, uu, vv, ww
  real(kind=8) :: c, s
  integer :: it
  integer :: ioff

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue
    ! Note : Les cas i = 0 et i = n/2 (s'il y a lieu), sont traites a part

    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 1,(n-1)/2
        it = n-i
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 1,(n-1)/2
          it = n-i
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    end if

    ! Cas particulier i = 0
    do j = 0,m/2-1
      work(       +j) = scale*x(debx+jumpx*(j    ))
      work(nwork/4+j) = scale*x(debx+jumpx*(j+m/2))
    end do

    ! Cas particulier i = n/2 (si n est pair seulement)
    if ( mod(n,2) == 0 ) then
      do j = 0,m/2-1
        work(        n/2*m/2+j) = scale*x(debx+incx*n+jumpx*(j    ))
        work(nwork/4+n/2*m/2+j) = scale*x(debx+incx*n+jumpx*(j+m/2))
      end do
    end if

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    ioff = 0
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On reconstitue

    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue
    ! Note : Le cas i = 0 sera traite a part

    if (m >= 16 .or. n/2 < 8) then

      do i = 1,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) =  2*(tt-ww)
          work(nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 1,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = t+v
          uu = u-w
          vv = c*(t-v)+s*(u+w)
          ww = c*(u+w)-s*(t-v)
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) = tt-ww
          work(nwork/4+m*i+j) = -uu-vv
        end do
      end do

    end if

    ! Cas particulier i = 0
    do j = 0,m-1
      t = scale*x(debx       +j*jumpx)
      v = scale*x(debx+n*incx+j*jumpx)
      work(        j) = v+t
      work(nwork/4+j) = v-t
    end do

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

    ! On reconstitue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        end do
      end do

    else

!dir$ ivdep
      do j = 0, m-1
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        end do
      end do

    end if

  end if

end subroutine jmcsm1dxy
