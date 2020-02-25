!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! pi2opt.F90 --- TP3 : calcul de Pi par intégration numérique.
!!
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Wed Feb 14 15:22:09 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 10:40:56 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pi
  !                  __
  !  But : calcul de || par la methode des rectangles (point milieu).
  !
  !                   / 1
  !                  |       4            __          
  !                  |   ----------  dx = ||
  !                  |    1 + x**2
  !                 / 0

!$ use OMP_LIB
  implicit none

! Dimension par defaut de la taille des matrices
#ifndef VAL_N
#define VAL_N 30000000
#endif

  integer, parameter :: n=VAL_N
  real(kind=8)       :: f, x, a, h, Pi_estime, Pi_calcule, ecart, Pi_calcule_loc
  real(kind=8)       :: temps, t_cpu_0, t_cpu_1, t_cpu
  integer            :: i, ir, t1, t2, k
!$ integer           :: nb_taches
 
  ! Fonction instruction a integrer
  f(a) = 4.0_8 / ( 1.0_8 + a*a )

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution pi en parallele avec ",i2," threads")',nb_taches

  ! Valeur estimée de Pi
  Pi_estime = acos(-1.0_8)

  ! Longueur de l'intervalle d'integration.
  h = 1.0_8 / real(n,kind=8)

  ! Temps CPU de calcul initial
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t1, count_rate=ir)

  ! Boucle artificielle a ne pas toucher
  do k=1,100
     
     ! Calcul de Pi
     Pi_calcule = 0.0_8
     Pi_calcule_loc = 0.0_8
     !$OMP PARALLEL PRIVATE(x) FIRSTPRIVATE(Pi_calcule_loc)
     !$OMP DO SCHEDULE(RUNTIME)
     do i = 1, n
        x = h * ( real(i,kind=8) - 0.5_8 )
        Pi_calcule_loc = Pi_calcule_loc + f(x)
     end do
     !$OMP END DO
     !$OMP ATOMIC
     Pi_calcule = Pi_calcule + Pi_calcule_loc
     !$OMP END PARALLEL
     Pi_calcule = h * Pi_calcule
     
  enddo

  ! Temps elapsed final
  call system_clock(count=t2, count_rate=ir)
  temps=real(t2 - t1,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Ecart entre la valeur estimee et la valeur calculee de Pi.
  ecart = abs(Pi_estime - Pi_calcule)

  ! Impression du resultat.
  print '(//,3X,"Nombre d''intervalles       : ",I10,/  &
           &,3X,"| Pi_estime - Pi_calcule | : ",1PE10.3,/  &
           &,3X,"Temps elapsed              : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU                  : ",1PE10.3," sec.",//)',n,ecart,temps,t_cpu

end program pi
