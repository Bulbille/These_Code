! Test des sous-programmes de generation des sinus-cosinus
! pour les T.F. en sinus et cosinus

program tTable

  implicit none

  integer, parameter :: nmax = 32
  integer, parameter :: nmin = nmax/2
  integer, parameter :: ntable = nmax*4
  real(8), dimension(0:ntable-1) :: table
  real(8), dimension(0:ntable-1) :: tableRef
  real(8) :: pi
  integer :: n
  integer :: i

  pi = acos( -1._8 )

  ! Verification de jmtable
  do n = nmin, nmax
    table(:) = 0
    call jmtable( table, 2*n, 0, n )
    do i = 0, n-1
      tableRef(i)   = cos( 2 * pi * i / real( n, 8 ) )
      tableRef(n+i) = sin( 2 * pi * i / real( n, 8 ) )
    end do
    print *, maxval( abs( table(0:2*n-1) - tableRef(0:2*n-1) ) )
  end do
  print *

  ! Verification de jmTableSin
  do n = nmin, nmax
    table(:) = 0
    call jmTableSin( table, 3*n, 0, n )
    do i = 0, n-1
      tableRef(i)     = cos( 2 * pi * i / real( n, 8 ) )
      tableRef(n+i)   = sin( 2 * pi * i / real( n, 8 ) )
      tableRef(2*n+i) = sin(     pi * i / real( n, 8 ) )
    end do
    print *, maxval( abs( table(0:3*n-1) - tableRef(0:3*n-1) ) )
  end do
  print *

  ! Verification de jmTableSinCos
  do n = nmin, nmax
    table(:) = 0
    call jmTableSinCos( table, 4*n, 0, n )
    do i = 0, n-1
      tableRef(i)     = cos( 2 * pi * i / real( n, 8 ) )
      tableRef(n+i)   = sin( 2 * pi * i / real( n, 8 ) )
      tableRef(2*n+i) = sin(     pi * i / real( n, 8 ) )
      tableRef(3*n+i) = cos(     pi * i / real( n, 8 ) )
    end do
    print *, maxval( abs( table(0:4*n-1) - tableRef(0:4*n-1) ) )
  end do

end program tTable
