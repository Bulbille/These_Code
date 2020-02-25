! $Header: /opt/cvsroot/jmfft/lib/jmccm1d4.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1d4(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

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
  real(kind=8) :: x0,x1,x2,x3,y0,y1,y2,y3,t0,t1,t2,t3,u0,u1,u2,u3
  real(kind=8) :: x0px2,x0mx2,x1px3,x1mx3
  real(kind=8) :: y0py2,y0my2,y1py3,y1my3
  real(kind=8) :: c1, s1, c2, s2, c3, s3
  integer :: ioff1, ioff2
  integer :: it0,iu0,it1,iu1,it2,iu2,it3,iu3
  integer :: jt0,ju0,jt1,ju1,jt2,ju2,jt3,ju3

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*4**(p-i-1) >= 16 .or. 4**i < 8) then

      do k = 0,4**i-1

        ! Les sinus et cosinus
        c1 = table(itable+        mtable*  4**(p-i-1)*k)
        s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
        c2 = table(itable+        mtable*2*4**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
        c3 = table(itable+        mtable*3*4**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

        ! Les indices
        it0 = ioff1        +m*(k*4**(p-i)             )
        iu0 = ioff1+nwork/4+m*(k*4**(p-i)             )
        it1 = ioff1        +m*(k*4**(p-i)+  4**(p-i-1))
        iu1 = ioff1+nwork/4+m*(k*4**(p-i)+  4**(p-i-1))
        it2 = ioff1        +m*(k*4**(p-i)+2*4**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1))
        it3 = ioff1        +m*(k*4**(p-i)+3*4**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1))
        jt0 = ioff2        +m*( k        *4**(p-i-1))
        ju0 = ioff2+nwork/4+m*( k        *4**(p-i-1))
        jt1 = ioff2        +m*((k+  4**i)*4**(p-i-1))
        ju1 = ioff2+nwork/4+m*((k+  4**i)*4**(p-i-1))
        jt2 = ioff2        +m*((k+2*4**i)*4**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+2*4**i)*4**(p-i-1))
        jt3 = ioff2        +m*((k+3*4**i)*4**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+3*4**i)*4**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*4**(p-i-1)-1

          x0px2 = work(it0+jl) + (c2*work(it2+jl)-s2*work(iu2+jl))
          x0mx2 = work(it0+jl) - (c2*work(it2+jl)-s2*work(iu2+jl))
          y0py2 = work(iu0+jl) + (c2*work(iu2+jl)+s2*work(it2+jl))
          y0my2 = work(iu0+jl) - (c2*work(iu2+jl)+s2*work(it2+jl))
          x1px3 = (c1*work(it1+jl)-s1*work(iu1+jl))+(c3*work(it3+jl)-s3*work(iu3+jl))
          x1mx3 = (c1*work(it1+jl)-s1*work(iu1+jl))-(c3*work(it3+jl)-s3*work(iu3+jl))
          y1py3 = (c1*work(iu1+jl)+s1*work(it1+jl))+(c3*work(iu3+jl)+s3*work(it3+jl))
          y1my3 = (c1*work(iu1+jl)+s1*work(it1+jl))-(c3*work(iu3+jl)+s3*work(it3+jl))

          ! Il reste a multiplier par les twopi4
          work(jt0+jl) = (x0px2)+(x1px3)
          work(jt2+jl) = (x0px2)-(x1px3)
          work(ju0+jl) = (y0py2)+(y1py3)
          work(ju2+jl) = (y0py2)-(y1py3)
          work(jt1+jl) = (x0mx2)-(y1my3)
          work(jt3+jl) = (x0mx2)+(y1my3)
          work(ju1+jl) = (y0my2)+(x1mx3)
          work(ju3+jl) = (y0my2)-(x1mx3)

        end do

      end do

    else

      do jl = 0,m*4**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,4**i-1

          t0 = work(ioff1+jl        +m*(k*4**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*4**(p-i)+  4**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+  4**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*4**(p-i)+2*4**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*4**(p-i)+3*4**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  4**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
          c2 = table(itable+        mtable*2*4**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
          c3 = table(itable+        mtable*3*4**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi4
          work(ioff2+jl        +m*( k        *4**(p-i-1))) = x0+x1+x2+x3
          work(ioff2+jl+nwork/4+m*( k        *4**(p-i-1))) = y0+y1+y2+y3
          work(ioff2+jl        +m*((k+  4**i)*4**(p-i-1))) = x0-y1-x2+y3
          work(ioff2+jl+nwork/4+m*((k+  4**i)*4**(p-i-1))) = y0+x1-y2-x3
          work(ioff2+jl        +m*((k+2*4**i)*4**(p-i-1))) = x0-x1+x2-x3
          work(ioff2+jl+nwork/4+m*((k+2*4**i)*4**(p-i-1))) = y0-y1+y2-y3
          work(ioff2+jl        +m*((k+3*4**i)*4**(p-i-1))) = x0+y1-x2-y3
          work(ioff2+jl+nwork/4+m*((k+3*4**i)*4**(p-i-1))) = y0-x1-y2+x3

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

end subroutine jmccm1d4
