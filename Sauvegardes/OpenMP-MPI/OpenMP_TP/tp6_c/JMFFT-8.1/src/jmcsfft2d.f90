! $Header: /opt/cvsroot/jmfft/lib/jmcsfft2d.f90,v 1.3 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine csfft2d(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: m, n, ldx, ldy
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*ldx*m-1) :: x
  real(kind=8), intent(out), dimension(0:ldy*m-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*(n+m)-1) :: table
  real(kind=8), intent(inout), dimension(0:512*max(n,m)-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i, j
  integer :: ntable, nwork, ioff
  integer :: nfact, mfact
  integer, dimension(0:99) :: fact
  integer :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  logical :: debut, fin
  integer :: dimy, deby, incy, jumpy
  integer :: signe
  real(kind=8) :: scale_temp
  logical :: npair, mpair
  character(len=*), parameter :: nomsp = 'CSFFT2D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair et mpair
  npair = ( mod(n,2) == 0 )
  mpair = ( mod(m,2) == 0 )

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n/2+1) call jmerreur2(nomsp,10,ldx,n/2+1)
  if (ldy < n+2  ) call jmerreur2(nomsp,15,ldy,n+2  )
  if ( .not.npair .and. .not.mpair) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*(n+m)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    table(0:mfact-1) = fact(0:mfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0  ,n)
    call jmtable(table,ntable,100+2*n,m)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    fact(0:mfact-1) = nint(table(0:mfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m
  !nwork = 512*max(n,m)
  call jmgetnwork(nwork,512*max(n,m),4*max(n,m))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(n/2+1,4*m,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail sans permuter
    ! les dimensions (n en premier) pour faire d'abord la tf sur m
    ! On en profite pour premultiplier et pour tenir compte du signe
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        work(             n_temp*j+i-ideb) =       scale*x(2*i  +2*ldx*j)
        work(nwork_temp/4+n_temp*j+i-ideb) = isign*scale*x(2*i+1+2*ldx*j)
      end do
    end do

    ! Precautions a posteriori de symetrie hermitienne
    ! En fait on ne s'occupe que des termes qui doivent etre reels
    ! (on ne peut rien faire pour forcer la symetrie hermitienne des autres)

    if ( ideb == 0 ) then
      ! Pour j = 0
      work(nwork_temp/4) = 0
      ! Pour j = m/2
      if ( mod(m,2) == 0 ) then
        work(nwork_temp/4+n_temp*m/2) = 0
      end if
    end if

    if ( ifin == n/2 .and. mod(n,2) == 0 ) then
      ! Pour j = 0
      work(nwork_temp/4+ifin-ideb) = 0
      ! Pour j = m/2
      if ( mod(m,2) == 0 ) then
        work(nwork_temp/4+n_temp*m/2+ifin-ideb) = 0
      end if
    end if

    ! On fait la FFT complexe -> complexe sur la deuxieme dimension (m)
    ioff = 0
    call jmccm1d(n_temp,m,fact,100,nfact,table,ntable,100+2*n,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
         y(2*i  +ldy*j) = work(ioff+             n_temp*j+i-ideb)
         y(2*i+1+ldy*j) = work(ioff+nwork_temp/4+n_temp*j+i-ideb)
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

  ! On fait les T.F. sur l'autre dimension
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(m,2*n,nwork,debut,npair,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On fait la FFT complexe -> reel sur le premiere dimension (n)
    dimy = ldy*m   ; deby = jdeb*ldy   ; incy = 1 ; jumpy = ldy
    signe = 1
    scale_temp = real(1,kind=8)
    call jmcsm1dxy(m_temp,n,fact,100,0,table,ntable,100+0, &
    & work,nwork_temp, &
    & y,dimy,deby,incy,jumpy,y,dimy,deby,incy,jumpy,signe,scale_temp)

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

end subroutine csfft2d
