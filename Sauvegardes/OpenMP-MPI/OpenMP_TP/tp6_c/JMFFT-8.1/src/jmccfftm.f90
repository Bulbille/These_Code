! $Header: /opt/cvsroot/jmfft/lib/jmccfftm.f90,v 1.2 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine ccfftm(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: n, m, ldx, ldy
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*ldx*m-1) :: x
  real(kind=8), intent(out), dimension(0:2*ldy*m-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*n-1) :: table
  real(kind=8), intent(inout), dimension(0:4*n*m-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i, j
  integer :: ioff
  integer :: ntable, nwork
  integer :: nfact
  integer, dimension(0:99) :: fact
  character(len=*), parameter :: nomsp = 'CCFFTM'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n) call jmerreur2(nomsp,9,ldx,n)
  if (ldy < n) call jmerreur2(nomsp,14,ldy,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 4*n*m

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    do j = 0,m-1
      work(j+m*i)     =       scale*x(2*i  +2*ldx*j)
      work(j+m*(n+i)) = isign*scale*x(2*i+1+2*ldx*j)
    end do
  end do

  ! On appelle le sous-programme principal
  ioff = 0
  call jmccm1d(m,n,fact,100,0,table,ntable,100+0,work,nwork,ioff)

  ! On recopie le tableau de travail dans le tableau de sortie
  do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    do j = 0,m-1
      y(2*i  +2*ldy*j) =       work(ioff+j+m*i)
      y(2*i+1+2*ldy*j) = isign*work(ioff+j+m*(n+i))
    end do
  end do

end subroutine ccfftm
