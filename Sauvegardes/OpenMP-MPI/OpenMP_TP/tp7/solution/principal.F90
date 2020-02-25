!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! principal.f90 --- TP7 : résolution d'un système linéaire général multiple
!!               --- avec la méthode du Bi_CGSTAB.
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Jan  9 16:05:48 2001
!! Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
!! Dern. mod. le   : Thu Aug 01 14:45:30 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program principal
!$ use OMP_LIB
  implicit none
  integer, parameter             :: lda=VAL_LDA, ldb=lda, ldx=lda, n=VAL_N, m=VAL_M
  integer                        :: i, ir, t0, t1, nb_iter_min, nb_iter_max
  real(kind=8), dimension(lda,n) :: a
  real(kind=8), dimension(ldb,m) :: b
  real(kind=8), dimension(ldx,m) :: x
  real(kind=8), dimension(n)     :: pr
  real(kind=8)                   :: norme
  real(kind=8)                   :: temps, t_cpu_0, t_cpu_1, t_cpu
!$ integer                       :: nb_taches

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution bi_cgstab en parallele avec ",i2," threads")',nb_taches

  ! On initialise la matrice et le second membre
  call random_number(a)
  call random_number(b)

  ! On muscle la diagonale principale
  forall (i=1:n) a(i,i) = a(i,i) + 100.0_8

  ! On choisi le preconditionnement Jacobi=Diagonale(a)
  forall (i=1:n) pr(i) = a(i,i)

  ! Solution initiale.
  x(:,:) = 1.0_8

  ! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t0, count_rate=ir)

  ! Resolution du systeme lineaire multiple Ax=b en parallèle.
  !$OMP PARALLEL
  call bi_cgstab(n, m, a, lda, pr, b, ldb, x, ldx, norme, nb_iter_min, nb_iter_max)
  !$OMP END PARALLEL

  ! Temps elapsed final.
  call system_clock(count=t1, count_rate=ir)
  temps=real(t1 - t0,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final.
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Impression du resultat.
  print '(//,3X,"Nb iterations min. et max. : ",I4,1X,I4,/ &
           &,3X,"Norme                      : ",1PE10.3,/  &
           &,3X,"Temps elapsed              : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU                  : ",1PE10.3," sec.",//)', &
           & nb_iter_min,nb_iter_max,norme,temps,t_cpu

end program principal
