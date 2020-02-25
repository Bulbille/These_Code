!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! jacobi.f90 --- TP2 : resolution d'un systeme lineaire par la methode de jacobi
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Jan  9 16:07:17 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Mon Jan 21 16:04:34 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program jacobi
!$ use OMP_LIB
  implicit none

! Dimension par defaut de la taille des matrices
#ifndef VAL_N
#define VAL_N 1201
#endif
#ifndef VAL_D
#define VAL_D 800
#endif


  integer, parameter           :: n=VAL_N, diag=VAL_D
  integer                      :: i, j, ir, t0, t1, iteration=0
  real(kind=8), dimension(n,n) :: a
  real(kind=8), dimension(n)   :: x, x_courant, b
  real(kind=8)                 :: norme, temps, t_cpu_0, t_cpu_1, t_cpu
!$ integer                     :: nb_taches

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution jacobi en parallele avec ",i2," threads")',nb_taches

  ! Initialisation de la matrice et du second membre
  call random_number(a)
  call random_number(b)

  ! On muscle la diagonale principale de la matrice
  forall (i=1:n) a(i,i)=a(i,i)+diag

  ! Solution initiale
  x(:) = 1.0_8

  ! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t0, count_rate=ir)

  ! Resolution par la methode de Jacobi
  Jaco : do
     iteration = iteration + 1

     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do i = 1, n
        x_courant(i) = 0.
        do j = 1, i-1
           x_courant(i) = x_courant(i) + a(i,j)*x(j)
        end do
        do j = i+1, n
           x_courant(i) = x_courant(i) + a(i,j)*x(j)
        end do
        x_courant(i) = (b(i) - x_courant(i))/a(i,i)
     end do
     !$OMP END PARALLEL DO

     ! Test de convergence
     norme = maxval( abs(x(:) - x_courant(:)) )/real(n,kind=8)
     if( (norme <= epsilon(1.0_8)) .or. (iteration >= n) ) exit Jaco

     x(:) = x_courant(:)
  end do Jaco

  ! Temps elapsed final
  call system_clock(count=t1, count_rate=ir)
  temps=real(t1 - t0,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Impression du resultat
  print '(//,3X,"Taille du systeme   : ",I5,/       &
           &,3X,"Iterations          : ",I4,/       &
           &,3X,"Norme               : ",1PE10.3,/  &
           &,3X,"Temps elapsed       : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU           : ",1PE10.3," sec.",//)',n,iteration,norme,temps,t_cpu

end program jacobi
