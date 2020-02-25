program strassen
  implicit none
  integer(kind=8), parameter     :: n=2048, threshold=128
  real(kind=8), dimension(:,:)   :: A(n,n), B(n,n), C(n,n)
  real(kind=8)                   :: temps, t_cpu_0, t_cpu_1, t_cpu, erreur
  integer                        :: ir, t0, t1

  ! Initialisation des matrices A et B
  call random_number(A)
  call random_number(B)

  ! Temps CPU de reference
  call cpu_time(t_cpu_0)
  ! Temps elapsed de reference.
  call system_clock(count=t0, count_rate=ir)

  ! Calcul C=A*B par la methode recursive de Strassen
  call strassen_multiply(A, B, C, n)

  ! Temps elapsed final.
  call system_clock(count=t1, count_rate=ir)
  temps=real(t1 - t0,kind=8)/real(ir,kind=8)

  ! Temps CPU de calcul final.
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

  call calcul_erreur(A, B, C, n, erreur)

  print '(//,3X,"Erreur                     : ",1PE10.3,/ &
           &,3X,"Temps elapsed              : ",1PE10.3," sec.",/ &
           &,3X,"Temps CPU                  : ",1PE10.3," sec.",//)', &
       & erreur,temps,t_cpu

contains

recursive subroutine strassen_multiply(AA, BB, CC, nn)
  integer, intent(in)   :: nn
  real(8), intent(in)  :: AA(nn,nn), BB(nn,nn)
  real(8), intent(out) :: CC(nn,nn)
  real(8), dimension(nn/2,nn/2) :: A11, A21, A12, A22, B11, B21, B12, B22
  real(8), dimension(nn/2,nn/2) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7

  if(iand(nn,1) /= 0 .OR. nn < threshold) then
     CC = matmul(AA,BB)
  else
     A11 = AA(1:nn/2,1:nn/2)
     A21 = AA(nn/2+1:nn,1:nn/2)
     A12 = AA(1:nn/2,nn/2+1:nn)
     A22 = AA(nn/2+1:nn,nn/2+1:nn)
     B11 = BB(1:nn/2,1:nn/2)
     B21 = BB(nn/2+1:nn,1:nn/2)
     B12 = BB(1:nn/2,nn/2+1:nn)
     B22 = BB(nn/2+1:nn,nn/2+1:nn)
     call strassen_multiply(A11+A22, B11+B22, Q1, nn/2)
     call strassen_multiply(A21+A22, B11, Q2, nn/2)
     call strassen_multiply(A11, B12-B22, Q3, nn/2)
     call strassen_multiply(A22, -B11+B21, Q4, nn/2)
     call strassen_multiply(A11+A12, B22, Q5, nn/2)
     call strassen_multiply(-A11+A21, B11+B12, Q6, nn/2)
     call strassen_multiply(A12-A22, B21+B22, Q7, nn/2)
     CC(1:nn/2,1:nn/2) = Q1+Q4-Q5+Q7
     CC(nn/2+1:nn,1:nn/2) = Q2+Q4
     CC(1:nn/2,nn/2+1:nn) = Q3+Q5
     CC(nn/2+1:nn,nn/2+1:nn) = Q1+Q3-Q2+Q6
  end if
  return
end subroutine strassen_multiply

end program strassen
