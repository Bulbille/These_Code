!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prod_mat.f90 --- TP1 : produit de matrices : C = A * B
!!
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Wed Feb 14 15:22:09 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 10:39:56 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program prod_mat
!$ use OMP_LIB
  implicit none

! Dimension par defaut de la taille des matrices
#ifndef VAL_M
#define VAL_M 701
#endif
#ifndef VAL_N
#define VAL_N 801
#endif

  integer, parameter           :: m=VAL_M, n=VAL_N
  real(kind=8), dimension(m,n) :: a
  real(kind=8), dimension(n,m) :: b
  real(kind=8), dimension(m,m) :: c
  integer                      :: i, j, k, ir, t1, t2
  real(kind=8)                 :: temps, t_cpu_0, t_cpu_1, t_cpu
!$ integer                     :: nb_taches

  !$OMP PARALLEL
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ print '(//,3X,"Execution prod_mat en parallele avec ",i2," threads")',nb_taches

  ! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t1, count_rate=ir)

  !$OMP......................................................
  ! Intialisation des matrices A, B et C.
  !$OMP......................................................
  do j = 1, n
    do i = 1, m
      a(i,j) = real(i+j,kind=8)
    end do
  end do
  !$OMP......................................................

  !$OMP......................................................
  do j = 1, m
    do i = 1, n
      b(i,j) = real(i-j,kind=8)
    end do
  end do
  !$OMP......................................................

  !$OMP......................................................
  do j = 1, m
    do i = 1, m
      c(i,j) = 0.0_8
    end do
  end do
  !$OMP......................................................

  ! Produit de matrices
  !$OMP......................................................
  do j = 1, m
    do k = 1, n
      do i = 1, m
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do
  !$OMP......................................................
  !$OMP......................................................

  ! Temps elapsed final
  call system_clock(count=t2, count_rate=ir)
  temps=real(t2 - t1,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Impression du resultat.
  print '(//,3X,"Valeurs de m et n   : ",I5,I5/,              &
           & 3X,"Temps elapsed       : ",1PE10.3," sec.",/, &
           & 3X,"Temps CPU           : ",1PE10.3," sec.",/, &
           & 3X,"Resultat partiel    : ",2(1PE10.3,1X)," ... ",2(1PE10.3,1X),//)', &
           m,n,temps,t_cpu,c(2,2),c(3,3),c(m-2,m-2),c(m-1,m-1)

end program prod_mat
