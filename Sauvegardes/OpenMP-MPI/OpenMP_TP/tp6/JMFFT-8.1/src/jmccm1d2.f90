! $Header: /opt/cvsroot/jmfft/lib/jmccm1d2.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1d2(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  implicit none

  ! Arguments
  integer, intent(in) :: p, n, m
  integer, intent(in) :: ntable,itable,ntable2,mtable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: k, jl
  integer :: it1,iu1,it2,iu2
  integer :: jt1,ju1,jt2,ju2
  real(kind=8) :: x1, x2, y1, y2
  real(kind=8) :: c, s
  integer :: ioff1, ioff2

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  if (mod(p,2)==0) then

    ! Si p est pair, on peut travailler entierement en base 4
    call jmccm1d4(p/2,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff1)
    ioff = ioff1

  else

    ! On fait les premieres etapes en base 4
    call jmccm1d4(p/2,n,2*m,table,ntable,itable,ntable2,mtable*2,work,nwork,ioff1)
    ioff2 = nwork/2-ioff1
    ! On fait la derniere etape en base 2
    if (m >= 16 .or. 2**(p-1) < 8) then
      do k = 0,2**(p-1)-1

        ! Les sinus et cosinus
        c = table(itable+        mtable*k)
        s = table(itable+ntable2+mtable*k)

        ! Les indices
        it1 = ioff1        +m*(k*2  )
        iu1 = ioff1+nwork/4+m*(k*2  )
        it2 = ioff1        +m*(k*2+1)
        iu2 = ioff1+nwork/4+m*(k*2+1)
        jt1 = ioff2        +m*( k          )
        ju1 = ioff2+nwork/4+m*( k          )
        jt2 = ioff2        +m*((k+2**(p-1)))
        ju2 = ioff2+nwork/4+m*((k+2**(p-1)))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m-1
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          x2 = work(it2+jl)
          y2 = work(iu2+jl)
          work(jt1+jl) = x1 + ( x2*c - y2*s )
          work(ju1+jl) = y1 + ( x2*s + y2*c )
          work(jt2+jl) = x1 - ( x2*c - y2*s )
          work(ju2+jl) = y1 - ( x2*s + y2*c )
        end do
      end do
    else
      do jl = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,2**(p-1)-1
          x1 = work(ioff1+jl        +m*(k*2  ))
          y1 = work(ioff1+jl+nwork/4+m*(k*2  ))
          x2 = work(ioff1+jl        +m*(k*2+1))
          y2 = work(ioff1+jl+nwork/4+m*(k*2+1))
          ! Les sinus et cosinus
          c = table(itable+        mtable*k)
          s = table(itable+ntable2+mtable*k)
          work(ioff2+jl        +m*( k          )) = x1 + ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*( k          )) = y1 + ( x2*s + y2*c )
          work(ioff2+jl        +m*((k+2**(p-1)))) = x1 - ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*((k+2**(p-1)))) = y1 - ( x2*s + y2*c )
        end do
      end do
    end if

    ioff = ioff2

  end if

end subroutine jmccm1d2
