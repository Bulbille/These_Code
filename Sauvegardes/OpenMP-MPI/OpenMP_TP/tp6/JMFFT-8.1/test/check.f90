      program check

      implicit none

      real(kind=8), parameter ::  tol = 1.e-9
      real(kind=8) :: x, y
      integer :: ios1, ios2
      integer :: irc
      integer :: nl = 0

      open(10,file='temp1',form='formatted',status='old')
      open(11,file='temp2',form='formatted',status='old')

      irc = 0

      do

        nl = nl + 1

        read(10,*,iostat=ios1) x
        read(11,*,iostat=ios2) y

        if (ios1 /= 0 .and. ios2 /= 0) then
          exit
        else if (ios1 /= 0 .xor. ios2 /= 0) then
          irc = 1
          exit
        ! Mixture de test relatif et de test absolu
        else if ( abs(x-y) <= max(tol*(abs(x)+abs(y)),tol) ) then
          continue
        else
          irc = 1
          exit
        end if

      end do

      if (irc == 0) then
        print *,'Test OK'
      else
        print *,'Problemes ligne', nl
        stop 1
      end if

      end program check
