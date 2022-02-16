program treebranch
!**************  Modified to use parameter data file and allocatable arrays. nah 28.8.08

use f90_tools
use tree_pres
implicit none

integer            :: j,k
real(lng)          :: t,p,df,rho,mu,Lr,Lr2,Fr2,q,g,a
real(lng), allocatable :: p0(:),pnew(:),qnew(:),Freq(:),Omega(:)
complex(lng), allocatable :: P0_omega(:),p0C(:),Pbranch(:),Qbranch(:),Ptree(:,:)
character (len=30) :: fn, buffer

include 'readdata.f90'  ! nah 28.8.08 reads parameter values from a file.

allocate (p0(1:tmstps),pnew(1:tmstps),qnew(1:tmstps),Freq(tmstps+1),Omega(tmstps+1))
allocate (P0_omega(1:tmstps),p0C(1:tmstps),Pbranch(1:tmstps),Qbranch(1:tmstps),Ptree(0:Maxgen,1:tmstps))
!allocate (Qtree(0:Maxgen,1:tmstps))

allocate (Computed(0:Maxgen,0:Maxgen)) !  NAH
allocate (jbranches(0:Maxgen))         !  NAH

allocate (radStat(0:maxStat))
allocate (pStat(0:maxStat))

! Physical parameters
df     = 1/Period                            ! Frequency interval.
Freq   = (/ (j*df, j=-tmstps/2, tmstps/2) /) ! Frequency-vector (abscissae). 
Omega  = 2*pi*Freq                           ! Freq.-vector scaled by a factor 2pi.

rho    = 1.055				! Density of blood [g/cm^3].
mu     = 0.049				! Viscosity of blood [g/cm/s].

Lr     = 0.1            ! Characteristic radius of the vessels in the tree [cm].
Lr2    = Lr**2          ! = sq(Lr) The squared radius [cm2].
g      = 981.0          ! The gravitational force [cm/s^2].
q      = 10.0*Lr2       ! The characteristic flow [cm^3/s].
!Fr2    = (q**(1/2))/g/(Lr**5)! The squared Froudes number.
 trm_rst = trm_rst*q/rho/g/Lr

localmax = 0
branches = 0

!For 'p4e.2d' r_root = 0.55, r_min = 0.005 cm

! THIS IS WHERE THE CODE READS IN THE PRESSURE AT A TERMINAL VESSEL
!open (unit=2,file='p4e.dat')
open (unit=2,file='p_mouse.dat')
do j=1, tmstps
  read(2,*) t, p
  p0(j) = p
end do
close (unit=2)

p0(:) = p0(:)*1333.22/rho/g/Lr

! WRITE TO FILE
open (unit=2,file='p0.2d', status = 'replace')
do j = 1,tmstps
  write(2,*) p0(j)
end do
close (unit=2)
  
p0C = p0  ! creating a complex number with 0 imaginary part


P0_omega = FFTshift(bitreverse(FFT(p0C))) !Transform P to frequency domain

write(*,*) 'P0_omega(N/2+1)', P0_omega(tmstps/2 + 1)

! CALL comp_pres IN TRES_PRES.F90
! Get alpha side predictions
Ptree = comp_pres(tmstps,Omega,P0_omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,1,1)

! Added by MJC: Need to compute the flow as well. Using the average admittance does NOT work
!Qtree = comp_pres(tmstps,Omega,P0_omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,1,2)

write(*,*)'alpha branches = ',branches
do j = 0, branches + 1 ! nah 24/6/08 - to print pressures at end of final vessels.
      
  do k = 1, tmstps
    Pbranch(k) = Ptree(j,k)
    Qbranch(k) = Pbranch(k)/Computed(j,0)
!    write(*,*) 'ERROR NOT HERE'
  end do


  pnew = real(IFFT(bitreverse(FFTshift(Pbranch))),lng)
  qnew = real(IFFT(bitreverse(FFTshift(Qbranch))),lng)

  a = j
  if (a < 10) then
    write (buffer,'(I1)') floor(a)
  else 
    write (buffer,'(I2)') floor(a)
  end if
  fn = 'p' // trim(buffer) // '_alpha.2d'
  open (unit=21, file=fn, status='replace') 
  do k=1,tmstps
    write (21,*) pnew(k)*rho*g*Lr/1332.20, qnew(k)*q
  end do
  close(unit=21)
end do


! beta branch calculations
Ptree = comp_pres(tmstps,Omega,P0_omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,2,1)
! Added by MJC: Need to compute the flow as well. Using the average admittance does NOT work
!Qtree = comp_pres(tmstps,Omega,P0_omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,2,2)

write(*,*)'beta branches = ',branches
do j = 0, branches + 1
do k = 1, tmstps
Pbranch(k) = Ptree(j,k)
Qbranch(k) = Pbranch(k)/Computed(0,j)
end do
pnew = real(IFFT(bitreverse(FFTshift(Pbranch))),lng)
qnew = real(IFFT(bitreverse(FFTshift(Qbranch))),lng)
a = j
if (a < 10) then
write (buffer,'(I1)') floor(a)
else
write (buffer,'(I2)') floor(a)
end if
fn = 'p' // trim(buffer)// '_beta.2d'
open (unit=21, file=fn, status='replace')
do k=1,tmstps
write (21,*) pnew(k)*rho*g*Lr/1333.22, qnew(k)*q
end do
close(unit=21)
end do

end program treebranch
