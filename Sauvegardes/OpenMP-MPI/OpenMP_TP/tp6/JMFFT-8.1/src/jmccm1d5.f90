! $Header: /opt/cvsroot/jmfft/lib/jmccm1d5.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1d5(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  implicit none

  ! Arguments
  integer, intent(in) :: p, n, m
  integer, intent(in) :: ntable,itable,ntable2, mtable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: i, k, jl
  real(kind=8) :: x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,t0,t1,t2,t3,t4,u0,u1,u2,u3,u4
  real(kind=8) :: c1, s1, c2, s2, c3, s3, c4, s4
  integer :: ioff1, ioff2

  ! Gestion des constantes cosinus
  real(kind=8), save :: twopi5
  !$OMP THREADPRIVATE(twopi5)
  real(kind=8), save :: ctwopi51, ctwopi52, ctwopi53, ctwopi54
  !$OMP THREADPRIVATE(ctwopi51,ctwopi52,ctwopi53,ctwopi54)
  real(kind=8), save :: stwopi51, stwopi52, stwopi53, stwopi54
  !$OMP THREADPRIVATE(stwopi51,stwopi52,stwopi53,stwopi54)
  logical, save :: first = .true.
  !$OMP THREADPRIVATE(first)

  ! On recupere cos et sin de 2*pi/5
  if (first) then
    first = .false.
    twopi5   = 2*acos(real(-1,kind=8))/real(5,kind=8)
    ctwopi51 = cos(  twopi5)
    stwopi51 = sin(  twopi5)
    ctwopi52 = cos(2*twopi5)
    stwopi52 = sin(2*twopi5)
    ctwopi53 = cos(3*twopi5)
    stwopi53 = sin(3*twopi5)
    ctwopi54 = cos(4*twopi5)
    stwopi54 = sin(4*twopi5)
  end if

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*5**(p-i-1) >= 16 .or. 5**i < 8) then

      do k = 0,5**i-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*5**(p-i-1)-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        end do

      end do

    else

      do jl = 0,m*5**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,5**i-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

end subroutine jmccm1d5
