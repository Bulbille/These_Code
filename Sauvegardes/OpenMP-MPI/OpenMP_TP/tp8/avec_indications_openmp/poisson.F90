!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: Fortran -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!
!! poisson.f90 --- TP8 : résolution d'un problème de poisson 2D par une méthode
!!             --- mixte (différences finies en X et FFT en sinus en Y)
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Fri Feb  5 14:11:00 1999
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 14:13:58 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program poisson
!$ use OMP_LIB
  implicit none

  integer, parameter                  :: ni=2501, nj=4001
  integer                             :: i, j, ir, t_debut, t_fin, code
  integer                             :: niter_min, niter_max
!$ integer                            :: nb_taches
  real(kind=8), dimension(ni,nj)      :: a_inf, a_diag, a_sup, b, u, u_a
  real(kind=8), dimension(ni)         :: x
  real(kind=8), dimension(nj)         :: y, vp, temp
  real(kind=8)                        :: Lx, Ly, hx, hy, ecart, pi, cx, cy
  real(kind=8)                        :: temps, t_cpu_0, t_cpu_1, t_cpu
  real(kind=8), dimension(100+3*nj)   :: table
  real(kind=8), dimension(4*(nj/2+1)) :: work
  external                            :: C06HAF

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution poisson en parallele avec ",i2," threads")',nb_taches

  ! Calcul des parametres physiques.
  Lx = 1.0_8        ! longueur du domaine en X
  Ly = 1.0_8        ! longueur du domaine en Y
  pi = acos(-1.0_8) ! Calcul de Pi

  ! Calcul des parametres de la grille.
  hx   = Lx / real(ni-1, kind=8) ! pas d'espace en X
  hy   = Ly / real(nj-1, kind=8) ! pas d'espace en Y
  cx   = 1.0_8 / (hx*hx)
  cy   = 1.0_8 / (hy*hy)
  x(:) = (/ (real(i-1, kind=8)*hx, i=1,ni) /) ! coordonnees en X
  y(:) = (/ (real(j-1, kind=8)*hy, j=1,nj) /) ! coordonnees en Y

  ! Initialisation du temps CPU.
  call cpu_time(t_cpu_0)

  ! Initialisation du temps de restitution.
  call system_clock(count=t_debut, count_rate=ir)

  ! Initialisation de la solution a calculer "u", de la solution
  ! analytique "u_a" et des conditions limites.
  u(:,:) = 0.0_8
  !$OMP ...............................................................
  !$ nb_taches = ......................................................
  !$OMP ...............................................................
  do j = 1, nj
    do i = 1, ni
      u_a(i,j) = sin(pi * y(j)) * cos(pi * x(i))
    end do
    u(1,j)  = u_a(1,j)
    u(ni,j) = u_a(ni,j)
  end do
  !$OMP ...............................................................
  !$OMP ...............................................................
  do i = 1, ni
    u(i,1)  = u_a(i,1)
    u(i,nj) = u_a(i,nj)
  end do
  !$OMP ...............................................................
  !$OMP ...............................................................

  ! 1) Calcul des valeurs propres "vp" associees à l'operateur - d2 /dy2.
  ! 2) Initialisation de la matrice tridiagonale associée a l'operateur symetrique
  !    -(d2 /dx2 + d2 /dy2) dans la base des vecteurs propres.
  ! 3) Initialisation du second membre "b" en tenant compte des conditions limites.
  a_inf(:,:) = 0.0_8 ; a_diag(:,:) = 0.0_8 ; a_sup(:,:) = 0.0_8
  vp(:) = 0.0_8 ; b(:,:) = 0.0_8
  !$OMP ...............................................................
  !$OMP ...............................................................
  do j = 2, nj-1
    vp(j) = 4.0_8 * cy * sin(0.5*real(j-1,kind=8)*pi/real(nj-1,kind=8))**2
    do i = 2, ni-1
      a_inf(i,j) = - cx
      a_diag(i,j)=   2.0_8 * cx + vp(j)
      a_sup(i,j) = - cx
      b(i,j)     = - cy * (u_a(i,j+1) - 2.0_8 * u_a(i,j) + u_a(i,j-1)) &
                   - cx * (u_a(i-1,j) - 2.0_8 * u_a(i,j) + u_a(i+1,j))
    end do
    b(ni-1,j) = b(ni-1,j) + cx * u(ni,j)
    b(2,j)    = b(2,j)    + cx * u(1,j)
  end do
  !$OMP ...............................................................
  !$OMP ...............................................................
  do i = 2, ni-1
    b(i,2)    = b(i,2)    + cy * u(i,1)
    b(i,nj-1) = b(i,nj-1) + cy * u(i,nj)
  end do
  !$OMP ...............................................................
  !$OMP ...............................................................

  ! Initialisation des coefficients trigonométriques en vue d'une FFT en sinus.
  work(:)=0.0_8 ; code=0
  call C06HAF( 1, nj-1, vp(2), 'i', table, work, code )

  ! Calcul du second membre "b" dans la base des vecteurs propres (FFT en sinus).
  !$OMP ...............................................................
  code=0 ; temp(:)=0.0_8 ; work(:)=0.0_8
  !$OMP ...............................................................
  do i = 2, ni-1
    temp(2:nj-1) = b(i,2:nj-1)
    call C06HAF(1, nj-1, temp(2), 'r', table, work, code)
    b(i,2:nj-1) = temp(2:nj-1)
  end do
  !$OMP ..............................................................
  !$OMP ..............................................................

  ! Résolution de nj-2 systèmes lineaires indépendants
  ! de taille (ni-2,ni-2) dans la base des vecteurs propres.
  niter_min=ni ; niter_max=0
  !$OMP ..............................................................
  call gradient_conjugue(ni, nj, a_inf, a_diag, a_sup, b, u, niter_min, niter_max)
  !$OMP .............................................................

  ! Calcul de la solution physique dans la base canonique (FFT inverse).
  ! On effectue ni-2 FFT independantes de nj-2 séquences.
  !$OMP .............................................................
  code=0 ;  temp(:)=0.0_8 ; work(:)=0.0_8
  !$OMP .............................................................
  do i = 2, ni-1
    temp(2:nj-1) = b(i,2:nj-1)
    call C06HAF(1, nj-1, temp(2), 'r', table, work, code)
    u(i,2:nj-1) = temp(2:nj-1)
  end do
  !$OMP ............................................................
  !$OMP ............................................................

  ! Temps elapsed final
  call system_clock(count=t_fin, count_rate=ir)
  temps=real(t_fin - t_debut,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Norme de l'ecart entre la solution calculee et la solution analytique
  ecart = maxval(abs(u(:,:) - u_a(:,:))) / real(ni*nj, kind=8)

  ! Impression du resultat
  print '(/)'  
  print '(3X,"Min. & Max. iteration numbers : ",1X,2(I4,3X),/,   &
         & 3X,"Norm(|u - u_a|)               : ",2X,1PE10.3,//,   &
         & 1X,"Analytic solution : ",3(1PE9.2,1X),"... ",3(1PE9.2,1X),/, &
         & 1X,"Computed solution : ",3(1PE9.2,1X),"... ",3(1PE9.2,1X),//)',&
         & niter_min,niter_max,ecart, &
         & u_a(2,2),u_a(3,3),u_a(4,4),u_a(ni-3,nj-3),u_a(ni-2,nj-2),u_a(ni-1,nj-1), &
         & u(2,2),u(3,3),u(4,4),u(ni-3,nj-3),u(ni-2,nj-2),u(ni-1,nj-1)

  print '(//,3X,"Temps elapsed              : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU                  : ",1PE10.3," sec.",//)', &
           & temps,t_cpu

end program poisson
