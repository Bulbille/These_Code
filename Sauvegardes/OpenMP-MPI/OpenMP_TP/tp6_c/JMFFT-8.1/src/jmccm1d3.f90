! $Header: /opt/cvsroot/jmfft/lib/jmccm1d3.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1d3(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

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
  real(kind=8) :: x1, x2, x3, y1, y2, y3, t1, t2, t3, u1, u2, u3
  real(kind=8) :: c2, s2, c3, s3
  integer :: it1,iu1,it2,iu2,it3,iu3
  integer :: jt1,ju1,jt2,ju2,jt3,ju3
  real(kind=8) :: r,s,t,u

  ! Gestion des constantes cosinus
  real(kind=8), save :: ctwopi3, stwopi3
  !$OMP THREADPRIVATE(ctwopi3)
  !$OMP THREADPRIVATE(stwopi3)
  logical, save :: first = .true.
  !$OMP THREADPRIVATE(first)
  integer :: ioff1, ioff2

  ! On recupere cos et sin de 2*pi/3
  if (first) then
    first = .false.
    ctwopi3 = -real(1,kind=8)/real(2,kind=8)
    stwopi3 = sqrt(real(3,kind=8))/real(2,kind=8)
  end if

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*3**(p-i-1) >= 16 .or. 3**i < 8) then

      do k = 0,3**i-1

        ! Les sinus et cosinus
        c2 = table(itable+        mtable*  3**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
        c3 = table(itable+        mtable*2*3**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

        ! Les indices
        it1 = ioff1        +m*(k*3**(p-i)             )
        iu1 = ioff1+nwork/4+m*(k*3**(p-i)             )
        it2 = ioff1        +m*(k*3**(p-i)+  3**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*3**(p-i)+  3**(p-i-1))
        it3 = ioff1        +m*(k*3**(p-i)+2*3**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1))
        jt1 = ioff2        +m*( k        *3**(p-i-1))
        ju1 = ioff2+nwork/4+m*( k        *3**(p-i-1))
        jt2 = ioff2        +m*((k+  3**i)*3**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+  3**i)*3**(p-i-1))
        jt3 = ioff2        +m*((k+2*3**i)*3**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+2*3**i)*3**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*3**(p-i-1)-1

          r = (c2*work(it2+jl))-(s2*work(iu2+jl))
          s = (c2*work(iu2+jl))+(s2*work(it2+jl))
          t = (c3*work(it3+jl))-(s3*work(iu3+jl))
          u = (c3*work(iu3+jl))+(s3*work(it3+jl))
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          work(jt1+jl) = x1 + r + t
          work(ju1+jl) = y1 + s + u
          work(jt2+jl) = x1 + ctwopi3*(r+t) - stwopi3*(s-u)
          work(ju2+jl) = y1 + ctwopi3*(s+u) + stwopi3*(r-t)
          work(jt3+jl) = x1 + ctwopi3*(r+t) + stwopi3*(s-u)
          work(ju3+jl) = y1 + ctwopi3*(s+u) - stwopi3*(r-t)

        end do

      end do

    else

      do jl = 0,m*3**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,3**i-1

          t1 = work(ioff1+jl        +m*(k*3**(p-i)             ))
          u1 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)             ))
          t2 = work(ioff1+jl        +m*(k*3**(p-i)+  3**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+  3**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*3**(p-i)+2*3**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1)))

          ! Les sinus et cosinus
          c2 = table(itable+        mtable*  3**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
          c3 = table(itable+        mtable*2*3**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

          ! On premultiplie
          x1 = t1
          y1 = u1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi3
          work(ioff2+jl        +m*( k        *3**(p-i-1))) = &
          & x1 + x2                    + x3
          work(ioff2+jl+nwork/4+m*( k        *3**(p-i-1))) = &
          & y1 + y2                    + y3
          work(ioff2+jl        +m*((k+  3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2-stwopi3*y2 + ctwopi3*x3+stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+  3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2+stwopi3*x2 + ctwopi3*y3-stwopi3*x3
          work(ioff2+jl        +m*((k+2*3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2+stwopi3*y2 + ctwopi3*x3-stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+2*3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2-stwopi3*x2 + ctwopi3*y3+stwopi3*x3

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

end subroutine jmccm1d3
