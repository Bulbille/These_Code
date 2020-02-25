! $Header: /opt/cvsroot/jmfft/lib/jmtable.f90,v 1.2 2004/04/01 15:48:32 teuler v8 $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS, teuler@lcp.u-psud.fr
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmtable(table,ntable,itable,n)

  implicit none

  ! Arguments
  integer, intent(in) :: ntable, itable, n
  real(kind=8), intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
 ! real(kind=8), save :: twopi
  real(kind=8) :: twopi
  real(kind=8) :: temp, temp1, temp2

  integer :: i

  twopi = 2 * acos(real(-1,kind=8))

  ! Calcul des sinus et cosinus

  ! Si n est multiple de 4, astuces en serie
  if (mod(n,4) == 0) then
    ! On se debarrasse des cas limite
    table(itable+      0) =  1
    table(itable+n+    0) =  0
    table(itable+    n/4) =  0
    table(itable+n+  n/4) =  1
    table(itable+    n/2) = -1
    table(itable+n+  n/2) =  0
    table(itable+  3*n/4) =  0
    table(itable+n+3*n/4) = -1
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/4-1
      temp = cos(twopi*real(i,kind=8)/real(n,kind=8))
      table(itable+    i)     =  temp
      table(itable+    n/2-i) = -temp
      table(itable+    n/2+i) = -temp
      table(itable+    n-i)   =  temp
      table(itable+n+  n/4+i) =  temp
      table(itable+n+  n/4-i) =  temp
      table(itable+n+3*n/4+i) = -temp
      table(itable+n+3*n/4-i) = -temp
    end do

  ! Si n est simplement multiple de 2 (sans etre multiple de 4)
  else if (mod(n,2) == 0) then
    ! On se debarrasse des cas limite
    table(itable+    0) =  1
    table(itable+n+  0) =  0
    table(itable+  n/2) = -1
    table(itable+n+n/2) =  0
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/2-1
      temp1 = cos(twopi*real(i,kind=8)/real(n,kind=8))
      table(itable+      i) =  temp1
      table(itable+  n/2+i) = -temp1
      temp2 = sin(twopi*real(i,kind=8)/real(n,kind=8))
      table(itable+n+    i) =  temp2
      table(itable+n+n/2+i) = -temp2
    end do

  ! Si n est impair
  else
    ! On se debarrasse des cas limite
    table(itable+  0) =  1
    table(itable+n+0) =  0
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/2
      temp1 = cos(twopi*real(i,kind=8)/real(n,kind=8))
      table(itable+    i) =  temp1
      table(itable+  n-i) =  temp1
      temp2 = sin(twopi*real(i,kind=8)/real(n,kind=8))
      table(itable+n+  i) =  temp2
      table(itable+n+n-i) = -temp2
    end do

  end if

end subroutine jmtable
