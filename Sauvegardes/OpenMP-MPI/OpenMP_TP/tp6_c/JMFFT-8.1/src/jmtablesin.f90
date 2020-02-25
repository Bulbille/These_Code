subroutine jmTableSin(table,ntable,itable,n)

  ! On doit seulement calculer en plus les sinus de pi/n

  implicit none

  ! Arguments
  integer, intent(in) :: ntable, itable, n
  real(kind=8), intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
  real(kind=8) :: pi
  integer :: i
  integer :: lasti

  pi = acos(real(-1,kind=8))

  ! On commence par calculer les premiers termes
  call jmtable(table,ntable,itable,n)

  ! Pour les indices pairs, c'est facile (prendre dans la table precedente)
!dir$ ivdep
!ocl novrec
!cdir nodep
  do i = 0, (n-1)/2
    table( itable + 2*n + 2*i ) = table( itable + n + i )
  end do

  ! Pour les indices impairs, on fait ca en deux temps
!dir$ ivdep
!ocl novrec
!cdir nodep
  do i = 1, n/2, 2
    table( itable + 2*n + i ) = sin( pi * i / real( n, 8 ) )
  end do
  lasti = i
!dir$ ivdep
!ocl novrec
!cdir nodep
  do i = lasti, n-1, 2
    table( itable + 2*n + i ) = table( itable + 2*n + n-i )
  end do

end subroutine jmTableSin
