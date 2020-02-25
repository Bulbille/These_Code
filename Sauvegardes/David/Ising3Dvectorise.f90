!----------------------------------------------------------------------------------
! PROGRAM de résolution Ising 3D par algo Wolff
!----------------------------------------------------------------------------------

program Ising3D
  implicit none

  include 'global3D.var' 				!include global variables

  ! theoretical estimations of parameters for the 3D Ising model
  double precision, parameter 	:: Betac = 0.2216544    ! cf article, Beta critique
  double precision, parameter   :: Ksi0  = 0.501        ! correlation length for Ising 3D on cubic lattice above Tc
  double precision, parameter   :: nu    = 0.6301       ! critical exponent linking casimir force and reduced system size

  ! the system matrix
  integer(1),allocatable	:: Cfg(:)               ! configuration of the spins 
  integer,allocatable           :: ppv(:)               ! nearest neighbors of the spins
  integer                       :: A                    ! size of confining boundaries
  integer 			:: x, y, z, nth, neigh  ! variables of position

  ! parameters of MC steps
  integer	 		:: nbsf         	! number of spin flip in a MC step
  integer                       :: nbcf                 ! nbr of cluster flip in MC step, following the article
  double precision              :: C

  ! for data files writing
  character(len=7) 		:: Lstr='0000000'

  ! outputs parameters
  double precision 		:: Energy, Emoy,E2moy, Mmoy, M2moy, M4moy, MCsize, Khi, Cv,Binder
  integer 			:: Magnetization,Csize

  ! computation parameters
  integer                       :: Ti                    ! current temperature step in MC
  double precision              :: T, Beta               ! current temperature in MC
  integer 			:: bin, step, flip, stepSize,L,therm
  double precision 		:: Ebin, Mbin, Khibin, Cvbin

  ! compute error with bootstrap method
  double precision, allocatable :: ResE(:),ResM(:),ResCsize(:)
  double precision              :: errbootE,errbootM,errbootE2,errbootM2,errbootM4,errbootCsize

  ! program time measurement
  real				:: t1,t2

  ! for random numbers generation
  integer 			:: idum=-1
  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface

  ! open the input file
  open(10,file='input3D.dat')
  read(10,*) LX 				!size of the spin matrix
  read(10,*) LY 				!size of the spin matrix
  read(10,*) LZ 				!size of the spin matrix
  read(10,*) Tmin 				!max scanned temperature
  read(10,*) Tmax 				!main scanned temperature
  read(10,*) nT 				!number of temp
  read(10,*) J
  read(10,*) binnb
  read(10,*) Nboot
  close(10)

  ! initialize Results tables
  allocate(ResE(1:binnb))
  allocate(ResM(1:binnb))
  allocate(ResCsize(1:binnb))

  do stepSize = 1,6

     L = stepSize*15
     LX=L
     LY=L
     LZ=L

     ! System dimensions
     N = L*L*L
     A = L*L

     ! initialize SpinMatrix and ppv
     allocate(Cfg(0:N-1))
     allocate(ppv(1:6*N))

     ! Initial Configuration
     Cfg = 1

     ! set the ppv table
     do z=0,LZ-1
        do y=0,LY-1
           do x=0,LX-1
              nth = x + LX*y + LX*LY*z
              ppv((6*nth)+1) = nth - x + modulo(x-1,LX)
              ppv((6*nth)+2) = nth - x + modulo(x+1,LX)
              ppv((6*nth)+3) = nth - LX*y + LX*modulo(y-1,LY)
              ppv((6*nth)+4) = nth - LX*y + LX*modulo(y+1,LY)
              ppv((6*nth)+5) = nth - LX*LY*z + LX*LY*modulo(z-1,LZ)
              ppv((6*nth)+6) = nth - LX*LY*z + LX*LY*modulo(z+1,LZ)
           enddo
        enddo
     enddo

     write(Lstr, '(i7)') N
     open(unit=1,form='formatted',status='replace',action='write',file='Ising3D_'//trim(adjustl(Lstr))//'.dat')
     write(1,'(a6,es14.5,a7,i10,a7,i10,a7,i10)') '# J = ', J, '  LX = ', LX, '  LY = ', LY, '  LZ = ', LZ
     write(1,'(16a14)') '#      T      ','       E      ','      dE      ','       M      ','      dM      ',& 
          '      Cv      ','       dCv    ','      Khi     ','      dKhi    ','   Binder     ','   dBinder    ',&
          '     MCsize   ','   dMCsize    ','     L        ','     N        ','    calctime  '
     write(1,'(16a14)') '#      1      ','       2      ','       3      ','       4      ','       5      ',&
          '       6      ','       7      ','        8     ','       9      ','      10      ','      11      ',&
          '       12     ','      13      ','      14      ','      15      ','      16      '
     close(1)

     ! set the nbr of Metropolis spin flips in a MC step
     nbsf = 3*LX*LY*LZ                           

     ! start Monte Carlo
     ! chose current temperature
     do Ti=0,nT
        T = Tmin + real(Ti)*(Tmax - Tmin)/real(nT)
        Beta = 1.d0/T

        ! start counting computation time
        call CPU_TIME( t1 )

        ! thermalization
        therm = floor(real(N)/1d2)
        do step=1,therm
           call ClusterFlipWithCsizeVect(Cfg, ppv, Beta, idum, Csize )
        enddo

        ! chosen # of calls to Wolff in a MC step depending on T following an empirical law
        if (T < 3.5) then
           write(*,*) 'ERROR : nbcf set to arbitrary 1 value'
           C=1
           nbcf=1
        else if ((T<=4.5).AND.(T>=3.5)) then
           C = -(56d-2)*T**2 + (38d-1)*T - (556d-2)
           nbcf = floor(1d0/C)
        else if ((T>4.5).AND.(T<4.6)) then 
           C = 5d-2
           nbcf = 20
        else if ((T>=4.6).AND.(T<=5.5)) then
           C = (12d-5)*(1d0/(T-(45d-1)))
           nbcf = floor(1d0/C)
        else if (T>5.5) then
           write(*,*) 'ERROR : nbcf set to arbitrary 1000 value'
           C=1d-3
           nbcf=1000
        end if

        Emoy=0d0
        E2moy=0d0
        Mmoy=0d0
        M2moy=0d0
        M4moy=0d0
        MCsize=0d0

        ! MC
        do bin=1,binnb
           do flip=1,nbcf
              call ClusterFlipWithCsizeVect(Cfg, ppv, Beta, idum, Csize )
           enddo
           Energy = 0.d0; Magnetization = 0.d0
           do nth=0,N-1
                    Energy = Energy + Cfg(nth)*(Cfg(ppv(6*nth+2))+Cfg(ppv(6*nth+4))+Cfg(ppv(6*nth+6)))
                    Magnetization = Magnetization + Cfg(nth)
           enddo
           Energy = J*Energy

           Emoy=Emoy+Energy
           E2moy=E2moy+Energy**2
           Mmoy=Mmoy+real(abs(Magnetization))
           M2moy=M2moy+real(Magnetization)**2
           M4moy=M4moy+real(Magnetization)**4
           MCsize=MCsize+real(Csize)

           ResM(bin) = real(abs(Magnetization))
           ResE(bin)= Energy
           ResCsize(bin)= real(Csize)

        enddo

        Emoy=Emoy/real(binnb)
        E2moy=E2moy/real(binnb)
        Mmoy=Mmoy/real(binnb)
        M2moy=M2moy/real(binnb)
        M4moy=M4moy/real(binnb)
        MCsize=MCsize/real(binnb)

        Cv  = (Beta**2)*(E2moy-(Emoy**2))
        Khi = (1d0/real(N))*(M2moy-(Mmoy**2))
        Binder = 1d0-(M4moy/(3*(M2moy**2)))

        call bootstrapMean(ResM, idum, errbootM)
        call bootstrapFluctuation(ResM, idum, errbootM2)
        call bootstrapBinder(ResM, idum, errbootM4)
        call bootstrapMean(ResE, idum, errbootE)
        call bootstrapFluctuation(ResE, idum, errbootE2)
        call bootstrapMean(ResCsize, idum, errbootCsize)

        errbootM2=(1d0/real(N))*errbootM2
        errbootE2=(Beta**2)*errbootE2

        call CPU_TIME(t2)
        write(*,*) 'temps d''execution', t2-t1

        open(unit=22,form='formatted',status='old',action='write',position='append',&
             file='Ising3D_'//trim(adjustl(Lstr))//'.dat')    
        write(22,'(16es14.5)') T,Emoy/real(N),errbootE/real(N),Mmoy/real(N),errbootM/real(N),Cv/real(N),errbootE2/real(N),&
             Khi/real(N),errbootM2/real(N),Binder,errbootM4,MCsize/real(N),errbootCsize/real(N),real(L),real(N), t2-t1
        close(22)

     enddo
     deallocate(Cfg)
     deallocate(ppv)
  enddo
  deallocate(ResE)
  deallocate(ResM)
  deallocate(ResCsize)

end program Ising3D


!---------------------------------------------------
!
!---------------------------------------------------
function ran2(idum)
  implicit none
  integer, parameter :: IM1=2147483563, IM2=2147483399, IA1=40014, IA2=40692, IQ1=53668,&
       IQ2=52774, IR1=12211, IR2=3791, NTAB=32
  integer, parameter :: IMM1 = IM1-1, NDIV=1+IMM1/NTAB
  double precision, parameter :: AM=1./IM1, EPS=1.2e-7, RNMX=1.-EPS
  integer :: idum, j, k
  double precision :: ran2
  integer :: idum2=123456789,iy=0
  integer, dimension(NTAB) :: iv=0

  if ( idum <= 0 ) then
     idum = max( -idum, 1 )
     idum2 = idum
     do j = NTAB+8,1,-1
        k = idum/IQ1
        idum = IA1*( idum - k*IQ1 ) - k*IR1
        if ( idum < 0 ) idum = idum + IM1
        if (j <= NTAB) iv(j) = idum
     end do
     iy = iv(1)
  endif
  k = idum/IQ1
  idum = IA1*( idum - k*IQ1 ) - k*IR1
  if ( idum < 0 ) idum = idum + IM1
  k = idum2/IQ2
  idum2 = IA2*( idum2 - k*IQ2 ) - k*IR2
  if ( idum2 < 0 ) idum2 = idum2 + IM2
  j = 1 + iy/NDIV
  iy = iv(j) - idum2
  iv(j) = idum
  if( iy < 1 ) iy = iy + IMM1
  ran2 = min( AM*iy, RNMX )
  return
end function ran2


!---------------------------------------------------
!
!---------------------------------------------------
subroutine SpinFlip(Cfg, Beta, idum)
  implicit none

  include 'global3D.var' 				!include global variables

  integer(1), intent(inout) :: Cfg(0:LX-1,0:LY-1,0:LZ-1)
  double precision, intent(in) :: Beta
  integer, intent(inout) :: idum

  double precision :: DeltaE
  integer :: x, y, Z

  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface

  x = floor( LX*ran2(idum) );  y = floor( LY*ran2(idum) ); z = floor( LZ*ran2(idum) )

  DeltaE = Cfg(modulo(x+1,LX),y,z) + Cfg(modulo(x-1,LX),y,z) + Cfg(x,modulo(y+1,LY),z) + Cfg(x,modulo(y-1,LY),z) + &
           Cfg(x,y,modulo(z+1,LZ)) + Cfg(x,y,modulo(z-1,LZ))
  DeltaE = -2.*J*Cfg(x,y,z)*DeltaE

  if ( ran2(idum) < exp( -Beta*DeltaE ) ) then
     Cfg(x,y,z) = -Cfg(x,y,z)
  endif
end subroutine SpinFlip


!---------------------------------------------------
!
!---------------------------------------------------
subroutine ClusterFlipWithCsizeVect(Cfg, neighbor, Beta, idum, clustersize )
  implicit none

  include 'global3D.var' 

  integer(1), intent(inout) :: Cfg(0:(LX*LY*LZ)-1)
  integer, intent(inout) :: neighbor(1:6*(LX*LY*LZ))
  double precision, intent(in) :: Beta
  integer, intent(inout) :: idum
  integer,intent(inout) :: clustersize

  double precision, dimension(-1:1) :: p
  double precision :: DeltaE
  integer :: nb, s, spinnb, nbn, sn, nnnb
  integer, dimension(LX*LY*LZ) :: Cluster
  integer(1), dimension(LX*LY*LZ) :: ClusterSign

  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface


  p(1) = 1. - exp( 2*Beta*J ); p(-1) = 1. - exp( -2*Beta*J )

  nb = floor( LX*LY*LZ*ran2(idum) );  s = Cfg(nb)

  clustersize = 1; spinnb = 1
  Cluster(clustersize) = nb; ClusterSign(clustersize) = s
  !flip
  Cfg(nb) = -s

  do nnnb=1,6
     nbn = neighbor((6*nb)+nnnb); sn = Cfg(nbn)
     if ( ran2(idum) < p(s*sn) ) then
        clustersize = clustersize + 1;
        Cluster(clustersize) = nbn; ClusterSign(clustersize) = sn
        !flip
        Cfg(nbn) = -sn
     endif
  enddo

  do while ( spinnb < clustersize )
     spinnb = spinnb + 1
     nb = Cluster(spinnb); s = ClusterSign(spinnb)
     do nnnb=1,6
        nbn = neighbor((6*nb)+nnnb); sn = Cfg(nbn)
        if ( ran2(idum) < p(s*sn) ) then
           clustersize = clustersize + 1;
           Cluster(clustersize) = nbn; ClusterSign(clustersize) = sn
           !flip
           Cfg(nbn) = -sn
        endif
     enddo
  enddo

end subroutine ClusterFlipWithCsizeVect


!-------------------------------------------------------------------------------------------
! Programms to compute error using bootstrap method
!-------------------------------------------------------------------------------------------

subroutine bootstrapMean(Res, idum, errboot)

! Programm to compute error using bootstrap method : Res contains the binnb values of 
! a calculated variable (exemple Magnetization or Energy) and the programm calculates 
! errboot the error on the mean calculated using the binnb values

  include 'global3D.var' 				!include global variables

  double precision, intent(inout) :: Res(1:binnb)
  integer, intent(inout) :: idum
  double precision, intent(inout) :: errboot

  ! compute error with bootstrap method
  double precision              :: meanboot
  double precision              :: meanmeanboot
  integer                       :: stepb, stepb2,i

  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface

  ! bootstrap method
  errboot = 0d0
  meanmeanboot = 0d0
  do stepb = 1,Nboot
     meanboot = 0d0
     do stepb2 = 1,binnb
        i = 1 + floor((binnb-1)*(ran2(idum)))
        meanboot = meanboot + Res(i)
     end do
     errboot = errboot + (meanboot/real(binnb))**2
     meanmeanboot = meanmeanboot + meanboot/real(binnb)
  end do
  errboot = errboot/real(Nboot)
  meanmeanboot = meanmeanboot/real(Nboot)
  errboot = errboot - meanmeanboot**2
  errboot = sqrt(errboot)

end subroutine bootstrapMean


subroutine bootstrapFluctuation(Res, idum, errboot)

! Programm to compute error using bootstrap method : Res contains the binnb values of 
! a calculated variable (exemple Magnetization or Energy). The programm calculates 
! errboot the error on mean(Res^2)-(mean(Res))^2 calculated using the binnb values

  include 'global3D.var' 				!include global variables

  double precision, intent(inout) :: Res(1:binnb)
  integer, intent(inout) :: idum
  double precision, intent(inout) :: errboot

  ! compute error with bootstrap method
  double precision              :: meanA,meanA2
  double precision              :: meanmeanboot
  integer                       :: stepb, stepb2,i

  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface

  ! bootstrap method
  errboot = 0d0
  meanmeanboot = 0d0
  do stepb = 1,Nboot
     meanA = 0d0
     meanA2= 0d0
     do stepb2 = 1,binnb
        i = 1 + floor((binnb-1)*(ran2(idum)))
        meanA = meanA + Res(i)
        meanA2= meanA2 + (Res(i)**2)
     end do
     errboot = errboot + ((meanA2/real(binnb))-(meanA/real(binnb))**2)**2
     meanmeanboot = meanmeanboot + ((meanA2/real(binnb))-(meanA/real(binnb))**2)
  end do
  errboot = errboot/real(Nboot)
  meanmeanboot = meanmeanboot/real(Nboot)
  errboot = errboot - meanmeanboot**2
  errboot = sqrt(errboot)

end subroutine bootstrapFluctuation



subroutine bootstrapBinder(Res, idum, errboot)

! Programm to compute error using bootstrap method : Res contains the binnb values of 
! the magnetization. The programm calculates 
! errboot the error on the Binder cumulant 1-(mean(M^4)/3*(mean(M^2)^2))

  include 'global3D.var' 				!include global variables

  double precision, intent(inout) :: Res(1:binnb)
  integer, intent(inout) :: idum
  double precision, intent(inout) :: errboot

  ! compute error with bootstrap method
  double precision              :: meanA4,meanA2
  double precision              :: meanmeanboot
  integer                       :: stepb, stepb2,i

  interface
     function ran2(idum)
       double precision :: ran2
       integer, intent(inout) :: idum
     end function ran2
  end interface

  ! bootstrap method
  errboot = 0d0
  meanmeanboot = 0d0
  do stepb = 1,Nboot
     meanA4 = 0d0
     meanA2= 0d0
     do stepb2 = 1,binnb
        i = 1 + floor((binnb-1)*(ran2(idum)))
        meanA4 = meanA4 + (Res(i)**4)
        meanA2= meanA2 + (Res(i)**2)
     end do
     errboot = errboot + ( 1d0 - (meanA4/(3*(meanA2**2))) )**2
     meanmeanboot = meanmeanboot + ( 1d0 - (meanA4/(3*(meanA2**2))) )
  end do
  errboot = errboot/real(Nboot)
  meanmeanboot = meanmeanboot/real(Nboot)
  errboot = errboot - meanmeanboot**2
  errboot = sqrt(errboot)

end subroutine bootstrapBinder
