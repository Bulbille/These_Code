subroutine jmTableSinCos(table,ntable,itable,n)

  ! On doit calculer en plus les sinus et les cosinus de pi/n

  implicit none

  ! Arguments
  integer, intent(in) :: ntable, itable, n
  real(kind=8), intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
  real(kind=8) :: pi
  integer :: i
  integer :: lasti

  pi = acos(real(-1,kind=8))

  ! On commence par calculer les premiers termes et les sinus
  call jmtable(table,ntable,itable,n)
  call jmTableSin(table,ntable,itable,n)

  ! On va calculer les cosinus avec des astuces trigonometriques

  ! Si n est pair, on peut transformer les cosinus en sinus
  if ( mod(n,2) == 0) then

    ! On commence par les bornes
    table( itable + 3*n       ) = 1
    table( itable + 3*n + n/2 ) = 0

    ! On calcule sur le premier quadrant
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1, n/2 - 1
      table( itable + 3*n + i ) =   table( itable + 2*n + i + n/2 )
    end do

    ! Puis sur le second
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = n/2 + 1, n-1
      table( itable + 3*n + i ) = - table( itable + 2*n + i - n/2 )
    end do

  else

    ! On ne peut pas calculer les cosinus a l'aide des sinus
    ! On fait tout simplement comme dans jmTableSin

  ! Pour les indices pairs, c'est facile (prendre dans la table precedente)
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 0, (n-1)/2
      table( itable + 3*n + 2*i ) = table( itable + i )
    end do

  ! Pour les indices impairs, on fait ca en deux temps
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1, n/2, 2
      table( itable + 3*n + i ) = cos( pi * i / real( n, 8 ) )
    end do
    lasti = i
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = lasti, n-1, 2
      table( itable + 3*n + i ) = - table( itable + 3*n + n-i )
    end do

  end if

end subroutine jmTableSinCos
