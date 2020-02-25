! $Header: /opt/cvsroot/jmfft/lib/jmccm1dp.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmccm1dp(p,q,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! On fait m t.f. d'ordre q en base p (p**q)
  ! Note : n n'est pas utilise ici

  implicit none

  ! Arguments
  integer, intent(in) :: p, q, m
  integer, intent(in) :: ntable,itable,ntable2, mtable
  real(kind=8), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=8), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: i, k, jl, jp, kp
  real(kind=8) :: ck, sk, tk, uk, cpjk, spjk
  integer :: pqq, pi, pqi, pqii
  integer :: ikpr, ikpi, ijpr, ijpi
  integer :: itr, iti, jtr, jti
  integer :: ioff1, ioff2
  real(kind=8) :: c11, c12, c21, c22

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Pour le calcul des cos(2*pi/p)
  pqq = p**(q-1)

  ! Boucle sur les etapes
  do i = 0, q-1

    pi   = p**i
    pqi  = p**(q-i)
    pqii = p**(q-i-1)

    do k = 0,pi-1

      do jp = 0,p-1

        do jl = 0,m*pqii-1

          ijpr = ioff2 + jl + m*((k+jp*pi)*pqii)
          ijpi = ijpr + nwork/4

          work(ijpr) = 0
          work(ijpi) = 0

        end do

      end do

      do kp = 0,p-1

        itr = itable+mtable*kp*pqii*k
        iti = itr + ntable2
        ck = table(itr)
        sk = table(iti)

        do jp = 0,p-1

          ! Gymanstique infernale pour recuperer cos(2*pi/p) etc
          jtr = itable+mtable*pqq*mod(jp*kp,p)
          jti = jtr + ntable2
          cpjk = table(jtr)
          spjk = table(jti)
          c11 = (cpjk*ck-spjk*sk)
          c12 = (cpjk*sk+spjk*ck)
          c21 = (cpjk*sk+spjk*ck)
          c22 = (cpjk*ck-spjk*sk)

!dir$ ivdep
!ocl novrec
!cdir nodep
          do jl = 0,m*pqii-1

            ikpr = ioff1+jl+m*(k*pqi+kp*pqii)
            ikpi = ikpr + nwork/4
            tk = work(ikpr)
            uk = work(ikpi)

            ijpr = ioff2+jl+m*((k+jp*pi)*pqii)
            ijpi = ijpr + nwork/4

            work(ijpr) = work(ijpr) + tk*c11-uk*c12
            work(ijpi) = work(ijpi) + tk*c21+uk*c22

          end do

        end do

      end do

    end do

    ! On alterne les offsets
    ioff1 = nwork/2 - ioff1
    ioff2 = nwork/2 - ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

end subroutine jmccm1dp
