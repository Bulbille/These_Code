!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dependance.F90 --- TP9 : boucles avec dépendances
!!
!! Auteur          : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Créé le         : Wed Nov  5 10:27:35 2008
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 13:54:50 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program dependance
  implicit none

  integer, parameter             :: nx=VAL_NX, ny=VAL_NY
  real(kind=8), dimension(nx,ny) :: V, W
  integer                        :: i, j, k, ir, t1, t2
  real(kind=8)                   :: temps, t_cpu_0, t_cpu_1, t_cpu, norme

  ! Intialisation des matrices V et W
  do i = 1, nx
    do j = 1, ny
      V(i,j) = real(i+j,kind=8)/(nx+ny)
      W(i,j) = real(i+j,kind=8)/(nx+ny)
    end do
  end do

  ! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

  ! Temps elapsed de reference.
  call system_clock(count=t1, count_rate=ir)

  ! Boucles avec dependance
  do j = 2, ny
     do i = 2, nx
        V(i,j) =( V(i,j) + V(i-1,j) + V(i,j-1))/3
     end do
  end do

  ! Temps elapsed final
  call system_clock(count=t2, count_rate=ir)
  temps=real(t2 - t1,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  ! Verification de la justesse de la parallelisation
  ! Ne pas modifier cette partie SVP
  do j = 2, ny
     do i = 2, nx
        W(i,j) =( W(i,j) + W(i-1,j) + W(i,j-1))/3
     end do
  end do
  norme = 0.0
  do j = 2, ny
     do i = 2, nx
        norme = norme + (V(i,j)-W(i,j))*(V(i,j)-W(i,j))
     end do
  end do

  ! Impression du resultat.
  print '(//,3X,"Valeurs de nx et ny : ",I5,I5/,              &
           & 3X,"Temps elapsed       : ",1PE10.3," sec.",/, &
           & 3X,"Temps CPU           : ",1PE10.3," sec.",/, &
           & 3X,"Norme (PB si /= 0)  : ",1PE10.3,//)', &
           nx,ny,temps,t_cpu,norme

end program dependance
