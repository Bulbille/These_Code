!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! cg.f90 --- TP4 : resolution d'un systeme lineaire symetrique
!!        --- par une methode de gradient conjugue preconditionne
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Jan  9 16:05:48 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 10:39:12 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gradient_conjugue
!$ use OMP_LIB
  implicit none

! Dimension par defaut de la taille des matrices
#ifndef VAL_N
#define VAL_N 4201
#endif

  integer, parameter           :: n=VAL_N
  integer                      :: i, j, ir, t0, t1, iteration=0
  real(kind=8), dimension(n,n) :: a
  real(kind=8), dimension(n)   :: x, b, pr, z, r, p, q
  real(kind=8)                 :: norme, rho, rho_ancien, alpha, beta, gamma
  real                         :: temps, t_cpu_0, t_cpu_1, t_cpu
!$ integer                     :: nb_taches

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution cg en parallele avec ",i2," threads")',nb_taches

  ! Initialisation de la matrice et du second membre
  call random_number(a)
  call random_number(b)

  ! Temps CPU de calcul initial
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t0, count_rate=ir)

  !$OMP PARALLEL

  ! On symmetrise la matrice
  !$OMP WORKSHARE
  forall (i=2:n, j=1:n-1, i>j) a(j,i) = a(i,j)

  ! On muscle la diagonale principale
  forall (i=1:n) a(i,i) = a(i,i) + 800.0_8

  ! On choisit le preconditionnement Jacobi=Diagonale(a)
  forall (i=1:n) pr(i) = a(i,i)

  ! Solution initiale
  x(:) = 1.0_8

  ! Ecart initial
  r(:) = b(:) - matmul(a(:,:),x(:))
  !$OMP END WORKSHARE

  ! Resolution par la methode du gradient conjugue
  do
     !$OMP SINGLE
     rho = 0.0_8 ; gamma = 0.0_8
     iteration = iteration + 1
     !$OMP END SINGLE

     !$OMP WORKSHARE
     z(:) = r(:) / pr(:)
     !$OMP END WORKSHARE

     !$OMP DO REDUCTION(+:rho) SCHEDULE(RUNTIME)
     do i = 1, n
        rho = rho + r(i)*z(i)
     end do
     !$OMP END DO

     if (iteration > 1) then
        !$OMP SINGLE
        beta = rho / rho_ancien
        !$OMP END SINGLE

        !$OMP WORKSHARE
        p(:) = z(:) + beta*p(:)
        !$OMP END WORKSHARE
     else
        !$OMP WORKSHARE
        p(:) = z(:)
        !$OMP END WORKSHARE
     end if

     !$OMP DO REDUCTION(+:gamma) SCHEDULE(RUNTIME)
     do i = 1, n
        q(i)  = sum(a(i,:)*p(:))
        gamma = gamma  + p(i)*q(i)
     end do
     !$OMP END DO

     !$OMP SINGLE
     alpha = rho / gamma
     !$OMP END SINGLE

     !$OMP WORKSHARE
     x(:) = x(:) + alpha * p(:)
     r(:) = r(:) - alpha * q(:)

     ! Test de convergence
     norme = maxval( abs( r(:) ) )/real(n,kind=8)
     !$OMP END WORKSHARE
     if( (norme <= 10*epsilon(1.0_8)) .or. (iteration >= n) ) exit
     !$OMP SINGLE
     rho_ancien = rho
     !$OMP END SINGLE
  end do
!$OMP END PARALLEL

  ! Temps elapsed final.
  call system_clock(count=t1, count_rate=ir)
  temps=real(t1 - t0,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final.
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Impression du resultat
  print '(//,3X,"Taille systeme      : ",I6,/  &
           &,3X,"Iterations          : ",I4,/       &
           &,3X,"Norme du residu     : ",1PE10.3,/  &
           &,3X,"Temps elapsed       : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU           : ",1PE10.3," sec.",//)',n,iteration,norme,temps,t_cpu

end program gradient_conjugue
