! NOTE: This program is similar to treebranch, but only computes the mean
! pressure and flow throughout the tree. This is stored for every alpha beta branch
program treebranch_mean

use f90_tools
use tree_pres_mean
implicit none

!integer, parameter :: tmstps = 1024
!real(lng), parameter:: trm_rst = 0.1


integer            :: j,k,output
real(lng)          :: t,p,df,rho,mu,Lr,Lr2,qc,g,a,PeriodND,p_ref
real(lng), allocatable :: p0(:),pnew(:),qnew(:),Freq(:),Omega(:)
complex(lng), allocatable :: P0_omega(:),p0C(:),Pbranch(:),Qbranch(:)
complex(lng), allocatable :: Ptree(:,:), Qtree(:,:)
character (len=30) :: fn, buffer

include 'readdata.f90'  ! nah 28.8.08 reads parameter values from a file.

allocate (p0(1:tmstps),pnew(1:tmstps),qnew(1:tmstps),Freq(tmstps+1),Omega(tmstps+1))
allocate (P0_omega(1:tmstps),p0C(1:tmstps),Pbranch(1:tmstps),Qbranch(1:tmstps))
allocate (Ptree(0:Maxgen,1:tmstps),Qtree(0:Maxgen,1:tmstps))

allocate (Computed(0:Maxgen,0:Maxgen)) !  NAH
allocate (p0_network(0:MaxBranches,1:5),q0_network(0:MaxBranches,1:5)) !  MJC
allocate (pL_network(0:MaxBranches,1:5),qL_network(0:MaxBranches,1:5)) !  MJC
!allocate (jbranches(0:Maxgen))         !  NAH

!allocate (radStat(0:maxStat))
!allocate (pStat(0:maxStat))

! Physical parameters
rho    = 1.055_lng                ! Density of blood [g/cm^3].
mu     = 0.032_lng                ! Viscosity of blood [g/cm/s].

Lr     = 1.0            ! Characteristic radius of the vessels in the tree [cm].
Lr2    = Lr**2          ! = sq(Lr) The squared radius [cm2].
g      = 981.0_lng          ! The gravitational force [cm/s^2].
qc      = 10.0_lng*Lr2       ! The characteristic flow [cm^3/s].
!tmstps = 1024
!df     = 1.0/Period
PeriodND =Period*qc/Lr**3
df     = 1.0/PeriodND                         ! Frequency interval.
Freq   = (/ (j*df, j=-tmstps/2, tmstps/2) /) ! Frequency-vector (abscissae). 
Omega  = 2.0*pi*Freq                           ! Freq.-vector scaled by a factor 2pi.
trm_rst = trm_rst!*qc/rho/g/Lr


write(*,*) '***************************************'
write(*,*) 'Omega(N/2+1) = ', Omega(tmstps/2 + 1)



localmax = 0
branches = 0
!index_in = 0


! THIS IS WHERE THE CODE READS IN THE PRESSURE AT A TERMINAL VESSEL
p_ref = 2.0*1333.220_lng
open (unit=20,file='pterm.dat')
do j=1, tmstps
  read(20,*) p
  p0(j) = p
end do
close (unit=20)

p0(:) = (p0(:)*1333.220_lng-p_ref)!/rho/g/Lr

!! WRITE TO FILE
!open (unit=2,file='p0.2d', status = 'replace')
!do j = 1,tmstps
!  write(2,*) p0(j)
!end do
!close (unit=2)
!
p0C = p0  ! creating a complex number with 0 imaginary part


P0_omega = FFTshift(bitreverse(FFT(p0C))) !Transform P to frequency domain

write(*,*) 'P0_omega(N/2+1)', P0_omega(tmstps/2 + 1)

! CALL comp_pres IN TRES_PRES.F90
output = comp_pres(P0_omega(tmstps/2 + 1),trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,qc,g)


end program treebranch_mean
