! $Header: /opt/cvsroot/jmfft/lib/jmcsfft3d.f90,v 1.3 2004/04/01 15:48:31 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine csfft3d(isign,n,m,l,scale,x,ldx1,ldx2,y,ldy1,ldy2,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: m, n, l, ldx1, ldx2, ldy1, ldy2
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:2*ldx1*ldx2*l-1) :: x
  real(kind=8), intent(out), dimension(0:ldy1*ldy2*l-1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*(n+m+l)-1) :: table
  real(kind=8), intent(inout), dimension(0:512*max(n,m,l)-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  integer :: i, j, k
  integer :: ntable, nwork, ioff
  integer :: nfact, mfact, lfact
  integer, dimension(0:99) :: fact
  integer :: ideb, ifin, jdeb, jfin, kdeb, kfin, i1, i2, j1, j2
  integer :: nltemp, nmtemp, mltemp, nwork_temp, iwork, iwork0
  logical :: debut, fini
  logical :: npair, mpair, lpair
  character(len=*), parameter :: nomsp = 'CSFFT3D'
  logical :: jHermit, kHermit

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair
  npair = (mod(n,2) == 0)
  mpair = (mod(m,2) == 0)
  lpair = (mod(l,2) == 0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (l < 1) call jmerreur1(nomsp,8,l)
  if (ldx1 < n/2+1) call jmerreur2(nomsp,12,ldx1,n/2+1)
  if (ldy1 < n+2  ) call jmerreur2(nomsp,18,ldy1,n+2  )
  if (ldx2 < m) call jmerreur2(nomsp,13,ldx2,m)
  if (ldy2 < m) call jmerreur2(nomsp,20,ldy2,m)
  if (.not.mpair .and. .not.npair .and. .not.lpair) &
  & call jmerreur3(nomsp,25,n,m,l)

  ! Gestion de table
  ntable = 100+2*(n+m+l)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    call jmfact(l,fact,100,mfact,lfact)
    table(0:lfact-1) = fact(0:lfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0      ,n)
    call jmtable(table,ntable,100+2*n    ,m)
    call jmtable(table,ntable,100+2*(n+m),l)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    lfact = nint(table(mfact)) + mfact
    fact(0:lfact-1) = nint(table(0:lfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m*l
  !nwork = 512*max(n,m,l)
  call jmgetnwork(nwork,512*max(n,m,l),4*max(n,m,l))

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup3(n/2+1,m,4*l,nwork,debut,.true.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites

    do k = 0,l-1

      iwork = 0
      kHermit = ( k == 0 ) .or. ( mod( l, 2 ) == 0 .and. k == l/2 )

      do j = jdeb,jfin

        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin

        iwork0 = iwork
        jHermit = ( j == 0 ) .or. ( mod( m, 2 ) == 0 .and. j == m/2 )

!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2

          work(             iwork+k*nmtemp) = &
          &       scale*x(2*i  +2*ldx1*j+2*ldx1*ldx2*k)
          work(nwork_temp/4+iwork+k*nmtemp) = &
          & isign*scale*x(2*i+1+2*ldx1*j+2*ldx1*ldx2*k)
          iwork = iwork+1

        end do

        ! Precautions a posteriori de symetrie hermitienne
        ! En fait on ne s'occupe que des termes qui doivent etre reels
        ! (on ne peut rien faire pour forcer la symetrie hermitienne des autres)
        ! On doit examiner 8 termes
        if ( i1 == 0 .and. jHermit .and. kHermit ) &
          work(nwork_temp/4+iwork0+k*nmtemp) = 0
        if ( mod( n,2) == 0 .and. i2 == n/2 .and. jHermit .and. kHermit ) &
          work(nwork_temp/4+iwork0+i2-i1+k*nmtemp) = 0

      end do

    end do

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    call jmccm1d(nmtemp,l,fact,100,mfact,table,ntable,100+2*(n+m),work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +ldy1*j+ldy1*ldy2*k) = work(ioff+             iwork+k*nmtemp)
          y(2*i+1+ldy1*j+ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(n/2+1,l,4*m,nwork,debut,.true.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+j*nltemp) = &
          & y(2*i  +ldy1*j+ldy1*ldy2*k)
          work(nwork_temp/4+iwork+j*nltemp) = &
          & y(2*i+1+ldy1*j+ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    call jmccm1d(nltemp,m,fact,100,nfact,table,ntable,100+2*n    ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +ldy1*j+ldy1*ldy2*k) = work(ioff             +iwork+j*nltemp)
          y(2*i+1+ldy1*j+ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(m,l,4*(n/2+1),nwork,debut,npair,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do i = 0,n/2
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          work(             iwork+i*mltemp) = y(2*i  +ldy1*j+ldy1*ldy2*k)
          work(nwork_temp/4+iwork+i*mltemp) = y(2*i+1+ldy1*j+ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    call jmcsm1d(mltemp,n,fact,100,0    ,table,ntable,100+0      ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do i = 0,n-1
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          y(i+ldy1*j+ldy1*ldy2*k) = work(ioff+iwork+i*mltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

end subroutine csfft3d
