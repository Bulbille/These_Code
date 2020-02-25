!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! bi_cgstab.f90 --- TP7 : méthode du bi-gradient conjugué stabilisée
!!               --- (Réf. Templates for the Solution of Linear Systems. p.27)
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Jan  9 16:05:48 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 11:12:23 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bi_cgstab(n, m, a, lda, pr, b, ldb, x, ldx, norme, nb_iter_min, nb_iter_max)
  implicit none

  ! Variables passées en argument.
  integer, intent(in)                           :: lda, ldb, ldx, n, m
  real(kind=8), intent(in), dimension(lda,n)    :: a
  real(kind=8), intent(in), dimension(n)        :: pr
  real(kind=8), intent(in), dimension(ldb,m)    :: b
  real(kind=8), intent(inout), dimension(ldx,m) :: x
  real(kind=8), intent(out)                     :: norme
  integer, intent(out)                          :: nb_iter_min, nb_iter_max

  ! Variables locales a la procedure.
  real(kind=8), dimension(n) :: r, r_tilde, p, p_hat, s, s_hat, v, t
  real(kind=8)               :: norme_locale, rho, rho_old, alpha, beta, omega
  integer                    :: j, iteration

  nb_iter_min=n ; nb_iter_max=0 ; norme=0.0_8

  ! Resolution de m systemes lineaires independants.
  !$OMP.............................................................
  do j = 1, m

    ! Initialisation du compteur sur le nombre d'iterations.
    iteration = 0

    ! Ecart initiale.
    r(:)       = b(1:n,j) - matmul(a(1:n,1:n),x(1:n,j))
    r_tilde(:) = r(:)

    ! Resolution par la methode du bi-gradient conjugue stabilisee.
    BiCGSTAB : do
       iteration = iteration + 1

       rho = dot_product(r_tilde(:), r(:))
       if(rho == 0.0_8) then
         print *," La procedure Bi-CGSTAB a echoue, rho=0 !"
         exit BiCGSTAB
       end if
       if (iteration == 1) then
         p(:) = r(:)
       else
         beta = (rho / rho_old) * (alpha / omega)
         p(:) = r(:) + beta * (p(:) - omega * v(:))
       end if

       p_hat(:) = p(:) / pr(:)
       v(:) = matmul(a(1:n,1:n), p_hat(:))
       alpha = rho / dot_product(r_tilde(:),v(:))
       s(:) = r(:) - alpha * v(:)

       ! Test de convergence 1
       norme_locale = maxval( abs( s(:) ) )/real(n,kind=8)
       if( (norme_locale <= 10*epsilon(1.0_8)) .or. (iteration >= n) ) then
         x(1:n,j) = x(1:n,j) + alpha * p_hat(:)
         exit BiCGSTAB
       end if

       s_hat(:) = s(:) / pr(:)
       t(:) = matmul(a(1:n,1:n), s_hat(:))
       omega = dot_product(t(:), s(:)) / dot_product(t(:),t(:))
       if ( omega == 0.0_8 ) then
          print *," La procedure Bi-CGSTAB a echoue, omega=0 !"
          exit BiCGSTAB
       end if
       x(1:n,j) = x(1:n,j) + alpha * p_hat(:) + omega * s_hat(:)
       r(:) = s - omega * t(:)

       ! Test de convergence 2
       norme_locale = maxval( abs( r(:) ) )/real(n,kind=8)
       if( (norme_locale <= 10*epsilon(1.0_8)) .or. (iteration >= n) ) exit BiCGSTAB

       rho_old = rho
    end do BiCGSTAB

    !$OMP.............................................................
    nb_iter_min = min(nb_iter_min,iteration)
    !$OMP.............................................................
    nb_iter_max = max(nb_iter_max,iteration)
    !$OMP.............................................................
    norme       = max(norme,norme_locale)
  end do
  !$OMP.............................................................

end subroutine bi_cgstab
