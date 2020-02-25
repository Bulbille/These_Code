! Transformee de Fourier complexe-complexe multiple

program tjmcfftmlt

  implicit none

  integer, parameter :: m = 7

  integer, parameter :: n = 36
  real, dimension(0:n-1,0:m-1) :: ar, ai, ar2, ai2

  integer, parameter :: ntrigs = 2*n
  real, dimension(0:ntrigs-1) :: trigs

  integer, parameter :: nifax = 19
  real, dimension(0:nifax-1) :: ifax

  integer, parameter :: nwork = 4*m*n
  real, dimension(0:nwork-1) :: work

  integer :: isign
  integer :: i, j, k
  real :: twopi
  complex :: s
  integer :: inc, jump

  twopi = 2 * acos(real(-1))

  ! On prepare le tableau d'entree
  call random_number( ar )
  call random_number( ai )
  ar2 = ar
  ai2 = ai

  call cftfax(n,ifax,trigs)
  isign = 1
  print *,'jmcfftmlt ',n,m,isign
  inc = 1
  jump = n
  call cfftmlt(ar,ai,work,trigs,ifax,inc,jump,n,m,isign)

  ! On imprime le tableau de sortie
  open(10,file='temp1',status='unknown',form='formatted')
  write(10,'(e25.12)') ((ar(i,j),ai(i,j),i=0,n-1),j=0,m-1)

  ! Ce qu'il faut trouver
  open(11,file='temp2',status='unknown',form='formatted')

  ! On reprepare le tableau d'entree
  ar = ar2
  ai = ai2
  ! Et on calcule
  do j = 0,m-1
    do i = 0,n-1
      s = 0
      do k = 0,n-1
        s = s + cmplx(cos(twopi*i*k/real(n)),isign*sin(twopi*i*k/real(n))) &
            & * cmplx(ar(k,j),ai(k,j))
      end do
      write(11,'(e25.12)') s
    end do
  end do
  close(11)

end program tjmcfftmlt
