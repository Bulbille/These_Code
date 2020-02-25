!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! cg.f90 --- TP8 : résolution d'un système linéaire symétrique
!!        --- par une méthode de gradient conjugué préconditionnée
!!
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Jan  9 16:05:48 2001
!! Dern. mod. par  : Pierre-Francois Lavallee <lavallee@idris.fr>
!! Dern. mod. le   : Tue Dec 22 14:09:16 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gradient_conjugue(ni, nj, a_inf, a_diag, a_sup, b, u, niter_min, niter_max)
  implicit none

  integer                           :: ni, nj, i, j, k
  integer                           :: niter, niter_min, niter_max
  real(kind=8), dimension(ni,nj)    :: a_inf, a_diag, a_sup, b, u
  real(kind=8), dimension(ni)       :: pr, r, z, p, q
  real(kind=8)                      :: norme, rho, rho_ancien, alpha, beta, gamma

  p(:) = 0.0_8

  !$OMP .....................................................................
  do j = 2, nj-1

     ! Preconditionnement du type Jacobi=Diagonale(a)
     pr(2:ni-1) = a_diag(2:ni-1,j)

     ! Ecart initiale
     do i = 2, ni-1
       r(i) = b(i,j) - a_inf(i-1,j)*u(i-1,j) &
                     - a_diag(i,j)*u(i,j)    &
                     - a_sup(i+1,j)*u(i+1,j)
     end do

     ! Solution initiale
     b(2:ni-1,j) = u(2:ni-1,j)

     niter=0
     norme = 1.0_8
     do while ( (norme > 10*epsilon(1.0_8)) .and. (niter < ni) )
       niter = niter + 1
       rho=0.0_8 ; gamma=0.0_8

        do i = 2, ni-1
           z(i) = r(i) / pr(i)
           rho = rho + r(i)*z(i)
        end do

        if (niter > 1) then
           beta = rho / rho_ancien

           do i = 2, ni-1
              p(i) = z(i) + beta*p(i)
           end do
        else
           do i = 2, ni-1
              p(i) = z(i)
           end do
        end if

        do i = 2, ni-1
           q(i)  = a_inf(i-1,j)*p(i-1) + a_diag(i,j)*p(i) + a_sup(i+1,j)*p(i+1)
           gamma = gamma + p(i)*q(i)
        end do

        alpha = rho / gamma

        do i = 2, ni-1
           b(i,j) = b(i,j) + alpha * p(i)
           r(i)   = r(i) - alpha * q(i)
        end do

        ! Test de convergence
        norme = maxval( abs( r(2:ni-1) ) )/real(ni,kind=8)
        rho_ancien = rho
     end do ! while

     !$OMP .....................................................................
     niter_min=min(niter_min,niter)
     !$OMP .....................................................................
     niter_max=max(niter_max,niter)
  end do
  !$OMP .........................................................................

end subroutine gradient_conjugue
