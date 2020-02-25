subroutine calcul_erreur(AA, BB, CC, nn, error)
  implicit none
  integer, intent(in)  :: nn
  real(8), intent(in)  :: AA(nn,nn), BB(nn,nn), CC(nn,nn)
  real(8), intent(out) :: error

  error = sqrt(sum((CC-matmul(AA,BB))**2)) / nn
end subroutine calcul_erreur
