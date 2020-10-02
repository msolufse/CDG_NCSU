!  MODULE ROOT_IMP

!***************************************************************************
!*                                                                         *
!* Module: root_imp.f90                                                    *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!* Contains subroutines and functions needed to compute the impedance at   *
!* the root of a structured tree (with constant asymmetry-ratio and a      *
!* geometry where the size of all branches depends of the radius of the    *
!* branch at the root of the tree.                                         *
!*                                                                         *
!***************************************************************************
module root_imp
use f90_tools     ! Contains several tools used for when computing FFT.
implicit none

private
public impedance  ! The subroutine called from the c-program arteries.cxx
                  ! (the driver routine)
public impedance_init
public impedance_close


! A temporary matrix for storing root impedances in parts of the
! structured tree, those which are repeated because of the constant
! asymmetry ratio.

complex(lng), allocatable   :: Z_om(:)


contains

!***************************************************************************
!*                                                                         *
!* Subroutine: impedance                                                   *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!*                                                                         *
!* Compute the Impedance as a function of space and time at the root of a  *
!* structured tree.                                                        *
!*                                                                         *
!* This subroutine is called from the c++ subroutine arteries.cxx whenever *
!* a new end-tube (of the tree of larger arteries) is constructed. We then *
!* use the resulting impedances to get a time-dependent boundary-condition.*
!* It is called by:                                                        *
!*                                                                         *
!* Impedance (tmstps,Tper,ff*,rho,mu,r_root,r_min,z_xt,Lr,Fr2,q,g)         *
!*                                                                         *
!* where:                                                                  *
!*                                                                         *
!* tmstps    Number of sample points in time (Must be a power of 2)!       *
!* Tper      The period of the simulation, s.                              *
!* ff*       Constants to model Eh/r = ff1*exp(ff2*r)+ff3, g/cm/s^2.       *
!* ff*       Constants to model Eh/r=-ff1*tanh(ff2*(r-ff3))+ff4, g/cm/s^2. *
!* rho       Blood density, g/cm^3.                                        *
!* mu        Bloodplasma viscosity, g/cm/s.                                *
!* r_root    The radius of the root of the structured tree, cm.            *
!* r_min     The terminal radius of the structured tree, cm.               *
!* z_xt      Will return the impedances at the root of the structured tree *
!*           in the time domain.                                           *
!* Lr        Characteristic length scale (radius), cm.                     *
!* Fr2       The squared Froude number.                                    *
!* q         Characteristic flow-rate, cm^3/s.                             *
!* g         The gravitational force, cm/s^2.                              *
!*                                                                         *
!***************************************************************************
subroutine impedance (tmstps,Period,r_root,r_min,y_xt,R1,R2,CT,LT)
implicit none

  integer,   intent(in)      :: tmstps
  real(lng), intent(in)      :: Period,r_root,r_min,R1,R2,CT,LT
  real(lng)                  ::  y_xt(tmstps)
  !real(lng)                  ::  z_xt(tmstps)

  integer                    :: j
  real(lng)                  :: df, Freq(tmstps+1), Omega(tmstps+1), Z_om0,Z_om1,Zc,Zc1
  complex(lng)               :: ii, Z_hat(tmstps), Y_hat(tmstps)

  integer, parameter                  :: nbuf = 2, f1 = 10
  character (len=90)                  :: fn
  character (len=40), dimension(nbuf) :: buffer  ! Temporary strings
  integer k

  ! Physical parameters
  df     = 1/Period                            ! Frequency interval. 
  Freq   = (/ (j*df, j=-tmstps/2, tmstps/2) /) ! Frequency-vector (abscissae). 
  Omega  = 2*pi*Freq                           ! Freq.-vector scaled by a factor 2pi.
  
  !ii = cmplx(0,1);
  ii = cmplx(0.0_lng,1.0_lng,lng)
  
  !Z_om = (R1+R2+ii*Omega*CT*R1*R2)/(1+ii*Omega*CT*R2)
   Z_om = R1+((R2+ii*Omega*LT)/(1+ii*R2*Omega*CT-Omega**2*LT*CT))
   !Z_om = (R1+R2+ii*Omega*(CT*R1*R2+LT)-R1*Omega**2*LT*CT)/(1+ii*Omega*CT*R2-Omega*Omega*LT*CT)


 
  Z_om0 = real(Z_om(tmstps/2+1),lng)   ! Dirty hack, that makes Z_om real
  Z_om1 = real(Z_om(1),lng)                              ! first at the lowest frequency.

  Zc1 = 0
  do j = 5,11
     Zc = Zc1+ real(Z_om(tmstps/2+j))
     Zc1 = Zc
  end do
  Zc = Zc/7
  

!write(*,*) '  Z(0)  ','   R1  ', '   R2  ', '    RT   ','    CT  ','  root_R'
!write(*,'(6F8.3)') Z_om0 , R1,  R2,  R1+R2, CT, r_root
  ! Transform P back to the time domain.
  ! Divide by tmstps to approximate continuous inv. Fourier transform.
  ! In particular, amplitude must be independent of resolution.
 ! z_xt   = real(IFFT(bitreverse(FFTshift(Z_om)/Period)),lng)

  
  Z_hat = Z_om
  Y_hat = 1/Z_om
  y_xt   = real(IFFT(bitreverse(FFTshift(Y_hat)/Period)),lng)

  write (buffer(1),'(I4)') floor(1000*r_root)
  write (buffer(2),'(I4)') floor(100*r_min)
  do k = 1, nbuf
    buffer(k) = adjustl(buffer(k))
  end do
  fn = 'Zhat' // trim(buffer(1)) // '_' // trim(buffer(2))

  open (f1, file=fn, action='write') 
 
  do k=1,tmstps
    ! write (f1,'(3F26.16)') Omega(k)/Lr**3*q, Z_hat(k)*rho*g*Lr/q
      write (f1,'(3F26.16)') Omega(k), Z_hat(k)

  end do 
  close(f1)

  return
  
end subroutine impedance


subroutine impedance_init (tmstps)

  integer,   intent(in)  :: tmstps
 
  allocate(Z_om(tmstps))
  
end subroutine impedance_init


subroutine impedance_close 

  deallocate(Z_om)
    
end subroutine impedance_close

end module root_imp



