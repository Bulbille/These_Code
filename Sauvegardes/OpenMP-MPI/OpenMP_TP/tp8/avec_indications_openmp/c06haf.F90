! Goal     : Sinus FFT and inverse
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

SUBROUTINE C06HAF(m,n,x,init,table,work,ifail)
  IMPLICIT NONE

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  INTEGER, INTENT(in)                                        :: m, n
  REAL(kind=prec), INTENT(inout), DIMENSION(0:m*n-1)         :: x
  CHARACTER(len=1), INTENT(in)                               :: init
  REAL(kind=prec), DIMENSION(0:100+3*n-1), INTENT(inout)     :: table
  REAL(kind=prec), DIMENSION(0:4*(n/2+1)*m-1), INTENT(inout) :: work
  INTEGER, INTENT(inout)                                     :: ifail

  ! Variables locales
  INTEGER, SAVE          :: lastn = -1
  LOGICAL, SAVE          :: first = .TRUE.
  INTEGER                :: ntable, nwork
  INTEGER                :: n2, p2, n3, p3, n5, p5
  INTEGER                :: isign
  REAL(kind=prec)        :: scale
  INTEGER                :: i, j, ii
  REAL(kind=prec), SAVE  :: pi
  INTEGER                :: ioff
  REAL(kind=prec)        :: s, t, u

  ! Gestion de pi
  IF (first) pi = ACOS(-1.0_prec)

  ! Gestion de table
  ntable = 100+3*n

  ! Gestion de work
  nwork = 4*(n/2+1)*m

  ! Gestion de ifail
  ifail = 0

  ! Test sur M et N
  IF (m < 1) THEN
    ifail = 1
    RETURN
  ELSE IF (n < 1) THEN
    ifail = 2
    RETURN
  END IF

  ! Test sur isign
  IF (init == 'i' .OR. init == 'I') THEN
    !$OMP CRITICAL
    first = .FALSE.
    lastn = n
    !$OMP END CRITICAL
    ! On factorise
    CALL jmfact(n,n2,p2,n3,p3,n5,p5)
    table(0) = n2
    table(1) = p2
    table(2) = n3
    table(3) = p3
    table(4) = n5
    table(5) = p5
    table(99) = SQRT(2.0_prec/REAL(n,KIND=PREC))
    ! On genere la table normale
    CALL jmtable(table,ntable,100+0,n)
    ! On genere la table des sinus des i*pi/n
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    DO i = 0,n-1
      table(100+2*n+i) = SIN(pi*REAL(i,KIND=PREC)/REAL(n,KIND=PREC))
    END DO

  ELSE IF (init == 's' .OR. init == 'S') THEN

    ! On verifie que le precedent appel etait bon
    IF (first) THEN
      ifail = 4
      RETURN
    ELSE IF (lastn /= n) THEN
      ifail = 5
      RETURN
    END IF

  ELSE IF (init == 'r' .OR. init == 'R') THEN

    n2 = NINT(table(0))
    n3 = NINT(table(2))
    n5 = NINT(table(4))
    IF (n2*n3*n5 /= n) THEN
      ifail = 5
       RETURN
    END IF

  ELSE

    ifail = 3
    RETURN

  END IF

  ! On recupere la factorisation et scale
  n2 = NINT(table(0))
  p2 = NINT(table(1))
  n3 = NINT(table(2))
  p3 = NINT(table(3))
  n5 = NINT(table(4))
  p5 = NINT(table(5))
  scale = table(99)

  ! On forme le tableau de travail
  ! Note : on ne prend que les n-1 premiers elements de x
  ! Cas particulier i = 0
  DO j = 0,m-1
    work(j) = 0
  END DO
  IF (m > 16 .OR. n < 8) THEN
    DO i = 1,n-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
      DO j = 0,m-1
        s = table(100+2*n+i)
        t = x(m*(i-1)+j)
        u = x(m*(n-i-1)+j)
        work(m*i+j) = scale * (s*(t+u)+0.5_prec*(t-u))
      END DO
    END DO
  ELSE
    DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
      DO i = 1,n-1
        s = table(100+2*n+i)
        t = x(m*(i-1)+j)
        u = x(m*(n-i-1)+j)
        work(m*i+j) = scale * (s*(t+u)+0.5_prec*(t-u))
      END DO
    END DO
  END IF

  ! On appelle le sous-programme de transformee de Fourier
  ioff = 0
  ! On appelle le sous-programme principal
  CALL jmscm1d(m,n,n2,p2,n3,p3,n5,p5,table,ntable,100+0,work,nwork,ioff)

  ! On reconstitue x
  ! D'abord les indices impairs
  IF (m > 16 .OR. n/2 < 8) THEN
    DO i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
      DO j = 0,m-1
        x(m*(2*i+1)+j) = work(ioff+nwork/4+m*(i+1)+j)
      END DO
    END DO
  ELSE
    DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
      DO i = 0,n/2-1
        x(m*(2*i+1)+j) = work(ioff+nwork/4+m*(i+1)+j)
      END DO
    END DO
  END IF
  ! Ensuite les indices pairs
  ! Cas particulier indice 0
  DO j = 0,m-1
    x(j) = work(ioff+j)/2
  END DO
  ! Cas general : recurrence
  IF (m > 16 .OR. n/2 < 8) THEN
    DO i = 1,n/2-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
      DO j = 0,m-1
        x(m*(2*i)+j) = x(m*(2*i-2)+j) + work(ioff+m*i+j)
      END DO
    END DO
  ELSE
    DO j = 0,m-1

      ! Attention : Il faut vectoriser une recurrence
      ! En attendant, voici une version scalaire
      DO i = 1,n/2-1
        x(m*(2*i)+j) = x(m*(2*i-2)+j) + work(ioff+m*i+j)
      END DO
    END DO
  END IF

END SUBROUTINE C06HAF


subroutine jmfact(n,n2,p2,n3,p3,n5,p5)
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(out) :: n2, p2, n3, p3, n5, p5

  ! Variables locales
  integer :: m

  m = n

  ! Etude des puissances de deux
  p2 = 0
  n2 = 1
  do
    if (mod(m,2) == 0) then
      p2 = p2+1
      n2 = n2*2
      m  = m/2
    else
      exit
    end if
  end do

  ! Etude des puissances de trois
  p3 = 0
  n3 = 1
  do
    if (mod(m,3) == 0) then
      p3 = p3+1
      n3 = n3*3
      m  = m/3
    else
      exit
    end if
  end do

  ! Etude des puissances de cinq
  p5 = 0
  n5 = 1
  do
    if (mod(m,5) == 0) then
      p5 = p5+1
      n5 = n5*5
      m  = m/5
    else
      exit
    end if
  end do

  ! On verifie a posteriori la factorisation
  if (n2*n3*n5 /= n) then
    print *,'JMMFACT: ',n,' non factorizable in power of 2, 3 or 5'
    stop 1
  end if

end subroutine jmfact


subroutine jmtable(table,ntable,itable,n)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in)                              :: ntable, itable, n
  real(kind=prec), intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
  real(kind=prec), save :: twopi
  real(kind=prec)       :: temp, temp1, temp2

  integer :: i

  twopi = 2.0_prec * acos(-1.0_prec)

  ! Calcul des sinus et cosinus

  ! Si n est impair, pas d'astuce
  if (mod(n,2) /= 0) then
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do i = 0,n-1
      table(itable+  i) = cos(twopi*real(i,kind=prec)/real(n,kind=prec))
      table(itable+n+i) = sin(twopi*real(i,kind=prec)/real(n,kind=prec))
    end do

  ! Si n est multiple de 4, astuces en serie
  else if (mod(n,4) == 0) then
    ! On se debarrasse des cas limite
    table(itable+      0) =  1
    table(itable+n+    0) =  0
    table(itable+    n/4) =  0
    table(itable+n+  n/4) =  1
    table(itable+    n/2) = -1
    table(itable+n+  n/2) =  0
    table(itable+  3*n/4) =  0
    table(itable+n+3*n/4) = -1
    ! Cas general
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do i = 1,n/4-1
      temp = cos(twopi*real(i,kind=prec)/real(n,kind=prec))
      table(itable+    i)     =  temp
      table(itable+    n/2-i) = -temp
      table(itable+    n/2+i) = -temp
      table(itable+    n-i)   =  temp
      table(itable+n+  n/4+i) =  temp
      table(itable+n+  n/4-i) =  temp
      table(itable+n+3*n/4+i) = -temp
      table(itable+n+3*n/4-i) = -temp
    end do

  ! Si n est simplement multiple de 2 (sans etre multiple de 4)
  else
    ! On se debarrasse des cas limite
    table(itable+      0) =  1
    table(itable+n+    0) =  0
    table(itable+    n/2) = -1
    table(itable+n+  n/2) =  0
    ! Cas general
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do i = 1,n/2-1
      temp1 = cos(twopi*real(i,kind=prec)/real(n,kind=prec))
      table(itable+      i) =  temp1
      table(itable+  n/2+i) = -temp1
      temp2 = sin(twopi*real(i,kind=prec)/real(n,kind=prec))
      table(itable+n+    i) =  temp2
      table(itable+n+n/2+i) = -temp2
    end do

  end if

end subroutine jmtable


subroutine jmscm1d(m,n,n2,p2,n3,p3,n5,p5,table,ntable,itable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in)                               :: m, n, n2, p2, n3, p3, n5, p5
  integer, intent(in)                               :: ntable,itable
  real(kind=prec), intent(in), dimension(0:ntable-1)   :: table
  integer, intent(in)                               :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout)                            :: ioff

  ! Variables locales
  integer      :: ioff1, ioff2
  integer      :: i, j
  real(kind=prec) :: t, u, v, w
  real(kind=prec) :: c, s
  integer      :: is, it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue
    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do j = 0,m/2-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do i = 0,n-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    end if
        
    ! On fait m/2 t.f. complexes de longueur n
    call jmccm1d(m/2,n,n2,p2,n3,p3,n5,p5,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On regenere le resultat
    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
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
!CDIR NODEP
        do j = 0, m-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do i = 0, n/2-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    call jmccm1d(m,n,n2/2,p2-1,n3,p3,n5,p5,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! Maintenant, il faut reconstituer la t.f. reelle

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do j = 0,m-1
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do i = 0,n/2
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    end if

  end if

  ioff = ioff1

end subroutine jmscm1d



subroutine jmccm1d(m,n,n2,p2,n3,p3,n5,p5,table,ntable,itable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in) :: m, n, n2, p2, n3, p3, n5, p5
  integer, intent(in) :: ntable,itable
  real(kind=prec), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! On fait n3*n5 T.F. de n2 (qui est en puissances de 2)
  if (n2 /= 1) then
    call jmccm1d2(p2,n2,m*n3*n5,table,ntable,itable,n,n/n2,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n2 /= 1 .and. n3*n5 /= 1) then
    call jmcctranspcs(m,n,n2,n3*n5,table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n5*n2 T.F. de n3 (en puissances de 3)
  if (n3 /= 1) then
    call jmccm1d3(p3,n3,m*n5*n2,table,ntable,itable,n,n/n3,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n3 /= 1 .and. n5*n2 /= 1) then
    call jmcctranspcs(m*n2,n,n3,n5,table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n2*n3 T.F. de n5 (en puissances de 5)
  if (n5 /= 1) then
    call jmccm1d5(p5,n5,m*n2*n3,table,ntable,itable,n,n/n5,work,nwork,ioff)
  end if

end subroutine jmccm1d


subroutine jmccm1d2(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in) :: p, n, m
  integer, intent(in) :: ntable,itable,ntable2,mtable
  real(kind=prec), intent(in), dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: i, k, jl
  integer :: it1,iu1,it2,iu2
  integer :: jt1,ju1,jt2,ju2
  real(kind=prec) :: x1, x2, y1, y2
  real(kind=prec) :: c, s
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
!CDIR NODEP
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
!CDIR NODEP
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


subroutine jmccm1d3(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in)                               :: p, n, m
  integer, intent(in)                               :: ntable,itable,ntable2, mtable
  real(kind=prec), intent(in), dimension(0:ntable-1)   :: table
  integer, intent(in)                               :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout)                            :: ioff

  ! Gestion des constantes cosinus
  real(kind=prec), save :: ctwopi3, stwopi3
  logical, save      :: first = .true.

  ! Variables locales
  integer :: i, k, jl
  real(kind=prec) :: x1, x2, x3, y1, y2, y3, t1, t2, t3, u1, u2, u3
  real(kind=prec) :: c2, s2, c3, s3
  integer :: it1,iu1,it2,iu2,it3,iu3
  integer :: jt1,ju1,jt2,ju2,jt3,ju3
  real(kind=prec) :: cx2px3,cy2py3,sx2mx3,sy2my3 

  integer :: ioff1, ioff2

  ! On recupere cos et sin de 2*pi/3
  if (first) then
    !$OMP CRITICAL
    first   = .false.
    ctwopi3 = -1.0_prec/2.0_prec
    stwopi3 = sqrt(3.0_prec)/2.0_prec
    !$OMP END CRITICAL
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
!CDIR NODEP
        do jl = 0,m*3**(p-i-1)-1

          ! On premultiplie

          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          x2 = (c2*work(it2+jl))-(s2*work(iu2+jl))
          y2 = (c2*work(iu2+jl))+(s2*work(it2+jl))
          x3 = (c3*work(it3+jl))-(s3*work(iu3+jl))
          y3 = (c3*work(iu3+jl))+(s3*work(it3+jl))

          ! On calcule les quantites redondantes
          cx2px3 = x1 + ctwopi3*(x2+x3)
          cy2py3 = y1 + ctwopi3*(y2+y3)
          sx2mx3 = stwopi3*(x2-x3)
          sy2my3 = stwopi3*(y2-y3)

          work(jt1+jl) = x1 +      x2 +      x3
          work(ju1+jl) = y1 +      y2 +      y3
          work(jt2+jl) = cx2px3 - sy2my3
          work(ju2+jl) = cy2py3 + sx2mx3
          work(jt3+jl) = cx2px3 + sy2my3
          work(ju3+jl) = cy2py3 - sx2mx3

        end do

      end do

    else

      do jl = 0,m*3**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!CDIR NODEP
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


subroutine jmccm1d4(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in)                               :: p, n, m
  integer, intent(in)                               :: ntable,itable,ntable2, mtable
  real(kind=prec), intent(in), dimension(0:ntable-1)   :: table
  integer, intent(in)                               :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout)                            :: ioff

  ! Variables locales
  integer :: i, k, jl
  real(kind=prec) :: x0,x1,x2,x3,y0,y1,y2,y3,t0,t1,t2,t3,u0,u1,u2,u3
  real(kind=prec) :: x0px2,x0mx2,x1px3,x1mx3
  real(kind=prec) :: y0py2,y0my2,y1py3,y1my3
  real(kind=prec) :: c1, s1, c2, s2, c3, s3
  integer :: ioff1, ioff2
  integer :: tt
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
!CDIR NODEP
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
!CDIR NODEP
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


subroutine jmccm1d5(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in)                               :: p, n, m
  integer, intent(in)                               :: ntable,itable,ntable2, mtable
  real(kind=prec), intent(in), dimension(0:ntable-1)   :: table
  integer, intent(in)                               :: nwork
  real(kind=prec), intent(inout), dimension(0:nwork-1) :: work
  integer, intent(inout)                            :: ioff

  ! Gestion des constantes cosinus
  real(kind=prec), save :: twopi5
  real(kind=prec), save :: ctwopi51, ctwopi52, ctwopi53, ctwopi54
  real(kind=prec), save :: stwopi51, stwopi52, stwopi53, stwopi54
  logical, save      :: first = .true.

  ! Variables locales
  integer :: i, k, jl
  real(kind=prec) :: x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,t0,t1,t2,t3,t4,u0,u1,u2,u3,u4
  real(kind=prec) :: c1, s1, c2, s2, c3, s3, c4, s4
  integer :: ioff1, ioff2

  ! On recupere cos et sin de 2*pi/5
  if (first) then
    !$OMP CRITICAL
    first = .false.
    twopi5   = 2.0_prec*acos(-1.0_prec)/5.0_prec
    ctwopi51 = cos(  twopi5)
    stwopi51 = sin(  twopi5)
    ctwopi52 = cos(2.0_prec*twopi5)
    stwopi52 = sin(2.0_prec*twopi5)
    ctwopi53 = cos(3.0_prec*twopi5)
    stwopi53 = sin(3.0_prec*twopi5)
    ctwopi54 = cos(4.0_prec*twopi5)
    stwopi54 = sin(4.0_prec*twopi5)
    !$OMP END CRITICAL
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
!CDIR NODEP
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
!CDIR NODEP
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


subroutine jmcctranspcs(m,n,n2,n3,table,ntable,itable,work,nwork,ioff)
  implicit none

  ! Arguments
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(12)
  integer, intent(in) :: m, n
  integer, intent(in) :: n2, n3
  integer, intent(in) :: ntable,itable
  real(kind=prec), intent(in),  dimension(0:ntable-1) :: table
  integer, intent(in) :: nwork
  real(kind=prec), intent(inout),  dimension(0:nwork-1) :: work
  integer, intent(inout) :: ioff

  ! Variables locales
  integer :: i, j, k
  real(kind=prec) :: t, u, c, s
  integer :: ioff1, ioff2
  integer :: is

  ! Gestion des offsets
  ioff1 = ioff
  ioff2 = nwork/2-ioff

  ! Gestion du stride
  is = n/(n2*n3)

  if ( m >= 16 .or. (n2 < 8 .and. n3 < 8) ) then

    do i = 0,n2-1
      do j = 0,n3-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do k = 0,m-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else if ( n2 >= 16 .or. n3 < 8 ) then

    do j = 0,n3-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do i = 0,n2-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else

    do i = 0,n2-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
        do j = 0,n3-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  end if

  ioff = ioff2

end subroutine jmcctranspcs
