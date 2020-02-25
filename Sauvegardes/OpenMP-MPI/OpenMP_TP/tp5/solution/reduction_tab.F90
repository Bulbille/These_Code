!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! reduction_tab.f90 --- TP5 : réduction d’un tableau
!! 
!! Auteur          : Pierre-François LAVALLEE (CNRS/IDRIS) <lavallee@idris.fr>
!! Créé le         : Tue Jun 28 18:34:00 2011
!! Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
!! Dern. mod. le   : Thu Jul 25 10:53:34 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program reduction_tableaux
  !$ use OMP_LIB
  implicit none

! Dimension par defaut de la taille des matrices
#ifndef VAL_NMOLEC
#define VAL_NMOLEC 10000
#endif
#ifndef VAL_NMOL
#define VAL_NMOL 10000
#endif
#ifndef VAL_N
#define VAL_N 10
#endif

  integer, parameter :: nmolec=VAL_NMOLEC, nmol=VAL_NMOL, n=VAL_N
  integer            :: i, j, k, t0, t1, ir
  real(kind=kind(1.d0)), dimension(nmol,n,nmolec) :: tab
  real(kind=kind(1.d0)), dimension(nmol)          :: tab1, tab2, tab1c, tab2c
  real(kind=kind(1.d0))                           :: err, temps, t_cpu_0, t_cpu_1, t_cpu
!$ integer                                        :: nb_taches

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution reduction_tab en parallele avec ",i2," threads")',nb_taches

  tab2(1:nmol) = 0

  !$OMP PARALLEL PRIVATE(tab1)

  ! Initialisation du tableau
  ! First-touch pour garantir un fonctionnement optimal sur les systèmes NUMA
  !$OMP DO SCHEDULE(STATIC)
  do k=1,nmolec
    do j=1,n
      do i=1,nmol
        tab(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP SINGLE
  ! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t0, count_rate=ir)
  !$OMP END SINGLE

  ! Nid de boucle a paralleliser
  !$OMP DO REDUCTION(+:tab2) SCHEDULE(STATIC)
  do k=1,nmolec
     tab1(1:nmol) = 0

     do j=1,n
        do i=1,nmol
           tab1(i) = tab1(i) + tab(i,j,k)
        enddo
     enddo

     tab2(1:nmol) = tab2(1:nmol) + 2*tab1(1:nmol)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Temps elapsed final
  call system_clock(count=t1, count_rate=ir)
  temps=real(t1 - t0,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Verification des resultats
  tab2c(1:nmol) = 0
  do k=1,nmolec
     tab1c(1:nmol) = 0

     do j=1,n
        do i=1,nmol
           tab1c(i) = tab1c(i) + tab(i,j,k)
        enddo
     enddo

     tab2c(1:nmol) = tab2c(1:nmol) + 2*tab1c(1:nmol)
  enddo
  
  err=maxval(abs(tab2c-tab2)/abs(tab2c))

  ! Impression du resultat
  print '(//,3X,"Temps elapsed       : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU           : ",1PE10.3," sec.",/ &
           &,3X,"Erreur relative     : ",1PE10.3," ",//)',temps,t_cpu,err

end program reduction_tableaux
