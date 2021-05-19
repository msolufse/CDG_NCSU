module tree_pres
use f90_tools
implicit none

! Define these values apriori
integer, parameter:: Maxgen = 100
integer, parameter:: tmstps = 1024!8192
!real(lng), parameter:: Period = 0.85

  real(lng), save :: Period, trm_rst, ff1, ff2, ff3, r_root, r_min, alpha, beta, lrr!, tmstps
!  integer, save   :: Maxgen, tmstps, maxStat, Ctmstps
  complex(lng), allocatable :: Computed(:,:)
  integer, allocatable :: jbranches(:)


real(lng), parameter :: Hct = 0.45

complex(lng), allocatable   :: Z_om(:)

! integer localmax, branchesj(0:Maxgen), branches
integer localmax, branches

contains


!***************************************************************************
!* Function: Z0func                                                        *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            *
!***************************************************************************
recursive function Z0func (omega_k, alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g) result (Z0)
implicit none

  integer, intent(in)   :: alpha_pow, beta_pow
  real(lng), intent(in) :: omega_k,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g
  ! ADDED BY MJC
  real(lng)    :: alpha_update, beta_update
  integer      :: generations
real(lng)    :: nu, r,diameter,eta_D,mu_D,mu1, r_d, l, A, A_d, D, wom
real(lng)    ::  Hct_term, C_exp !,Hct
!  real(lng)    :: beta, alpha
  complex(lng) :: i, g_omega, c_omega, kappa, Z0, ZL, Zl_0, Zr_0, t1, t2

  ! Physical constants.
  i     = cmplx(0.0_lng,1.0_lng,lng) ! The complex unit.

!  alpha  = ((asym**(expo/2)+1.0)**(-1/expo))**alpha_pow  ! Scaling parameter.
!  beta = (sqrt(asym)*(asym**(expo/2)+1.0)**(-1/expo))**beta_pow   ! do.
  ! Added by MJC
    alpha_update = alpha**alpha_pow
    beta_update  = beta**beta_pow
  generations = alpha_pow + beta_pow  ! The current generation.
  r_d  = alpha_update*beta_update*r_root     ! Radius at root, cm.
  A_d  = pi*r_d**2             ! Cross-sectional area, cm^2.
  r    = r_d                   !b Radius at root, dimension-less.
  A    = A_d                   ! Cross-sectional area, dimension-less.
! ADDED BY MJC
    if (r<0.005) then
            l=15.75*(r**1.10) ! From Qureshi 2014
        else
            l = 1.79*(r**0.47) ! From Qureshi 2014
    !        l = lrr*r
        endif

    !l = 12.4*(r**1.10) ! from Olufsen 2012

    diameter = 2.0*r*(10.0**4.0)
    eta_D = 3.2+6*exp(-0.085*diameter)-2.44*exp(-0.06*(diameter**0.645))

    ! Terms that depend on hematocrit
    C_exp = (0.8+exp(-0.075*diameter))*(1.0/(1.0+10**(-11.0) * diameter**12.0) - 1.0) + (1.0/(1.0+10**(-11.0) * diameter**12.0))
    !    Hct = 0.45
    Hct_term = ((1.0 - Hct)**C_exp - 1.0)/((1.0-0.45)**C_exp - 1.0)
    mu_D  = (1.0 + (eta_D - 1.0)*Hct_term*(diameter/(diameter - 1.1))**2.0)*(diameter/(diameter - 1.1))**2.0
    mu1 = mu*mu_D/3.2

  nu   = mu1/rho                ! Kinematic blood viscosity, m^2/s.
  D    = 1.0/(ff1*exp(ff2*r_d)+ff3)*3*A_d/2 ! Distensibility.
  wom  = r_d*sqrt(omega_k/nu)  ! Womersleys parameter.
  ! Temporary functions of r.
  if (wom > 3.0) then 
    g_omega = sqrt(D*A/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)) 
    c_omega = sqrt(A/D/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom))
  else
    if (wom > 2.0) then
      g_omega = sqrt(D*A/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
      c_omega = sqrt(A/D/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
    else
      if (wom == 0) then
        g_omega = 0.0
        c_omega = 0.0
      else
        g_omega = sqrt(D*A/rho)*sqrt(i*wom**2/8+wom**4/48)
        c_omega = sqrt(A/D/rho)*sqrt(i*wom**2/8+wom**4/48)
      end if
    end if
  end if
 
  ! Temporary function of omega_k. 
  if (omega_k /= 0) then
    kappa  = omega_k*l/c_omega
  else
    kappa   = 0.0
  end if

  ! Determine the resistance of the root of the terminal branches.
  ! if (generations >= Maxgen) 
  if (r <= r_min) then
    if (generations >= localmax) then
      localmax = generations
    end if
    if (generations >= Maxgen) then
      write (*,*) 'Generations level exceeded'
    end if
    if (omega_k==0) then
        ZL = trm_rst!*q/rho/g/Lr
    else
        ZL=0.0
    endif
  else
    ! Get Z0 recursively at reduced frequencies.
    if (abs(Computed(alpha_pow+1, beta_pow)) /= 0.0) then
      Zl_0 = Computed(alpha_pow+1, beta_pow)
    else  
      Zl_0 = Z0func (omega_k,alpha_pow+1,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)
    end if

    if (abs(Computed(alpha_pow, beta_pow+1)) /= 0.0) then
      Zr_0 = Computed(alpha_pow, beta_pow+1)
    else
      Zr_0 = Z0func (omega_k,alpha_pow,beta_pow+1,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)
    end if

!    ZL = 1.0/(1.0/Zl_0 + 1.0/Zr_0)
    ZL = (Zl_0*Zr_0)/(Zl_0+Zr_0)

  end if
 
  ! Finally get Z0(omega) as theoretically derived.
  if (g_omega /= 0.0) then 
    t1   = i*sin(kappa)/g_omega + cos(kappa)*ZL
    t2   = cos(kappa) + i*g_omega*sin(kappa)*ZL
    Z0   = (t1/t2)
  else
    Z0  = (8.0*mu*l/(A*r**2.0) + ZL)
  end if

  ! This should be non-dimensional
  Computed(alpha_pow,beta_pow) = Z0!*qc/rho/g/Lr
!    write(*,*) 'Z0',Z0,alpha_pow,beta_pow
end function Z0func

!==============================================================================================
!====================function computes the pressure at the end of each vessel==================
function PLfunc (omega_k,Z0,P0,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,pq) result (result)

  integer, intent(in)   :: alpha_pow, beta_pow, pq
  real(lng), intent(in) :: omega_k,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g
  complex(lng), intent(in) :: P0, Z0
  
! ADDED BY MJC
  real(lng)    :: alpha_update, beta_update
  real(lng)    :: nu, r,diameter,eta_D,mu_D,mu1, r_d, l, A, A_d, D, wom
 real(lng)    :: Hct_term, C_exp!,Hct

!  real(lng)    :: beta, alpha
  complex(lng) :: i, g_omega, c_omega, kappa, Q0, PL,QL,result
  complex(lng) :: ZL,t1,t2

  ! Physical constants.
  i  = cmplx(0.0_lng,1.0_lng,lng) ! The complex unit.

 ! Added by MJC
    alpha_update = alpha**alpha_pow
    beta_update  = beta**beta_pow
  r_d  = alpha_update*beta_update*r_root     ! Radius at root, cm.
  A_d  = pi*r_d**2             ! Cross-sectional area, cm^2.
  r    = r_d                   ! Radius at root, dimension-less.
  A    = A_d                   ! Cross-sectional area, dimension-less.


if (r<0.005) then
        l=15.75*(r**1.10) ! From Qureshi 2014
    else
        l = 1.79*(r**0.47) ! From Qureshi 2014
!        l = lrr*r
    endif

!l = 12.4*(r**1.10) ! from Olufsen 2012

diameter = 2.0*r*(10**4)
    eta_D = 3.2+6*exp(-0.085*diameter)-2.44*exp(-0.06*(diameter**0.645))
    ! Terms that depend on hematocrit
    C_exp = (0.8+exp(-0.075*diameter))*(1.0/(1.0+10**(-11.0) * diameter**12.0) - 1.0) + (1.0/(1.0+10**(-11.0) * diameter**12.0))
!    Hct = 0.45
Hct_term = ((1.0 - Hct)**C_exp - 1.0)/((1.0-0.45)**C_exp - 1.0)
mu_D  = (1.0 + (eta_D - 1.0)*Hct_term*(diameter/(diameter - 1.1))**2.0)*(diameter/(diameter - 1.1))**2.0
!mu_D  = (1.0 + (eta_D - 1.0)*(diameter/(diameter - 1.1))**2.0)*(diameter/(diameter - 1.1))**2.0

    mu1 = mu*mu_D/3.2


  nu   = mu1/rho                ! Kinematic blood viscosity, m^2/s.
  D    = 1.0/(ff1*exp(ff2*r_d)+ff3)*3.0*A_d/2.0 ! Distensibility.
  wom  = r_d*sqrt(omega_k/nu)  ! Womersleys parameter.
!if (omega_k == 0) write(*,*) 'PLfunc',alpha_update,beta_update,r_d, A, l, diameter, eta_D, mu1, wom
  ! Temporary functions of r.
  if (wom > 3.0) then 
    g_omega = sqrt(D*A/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)) 
    c_omega = sqrt(A/D/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom))
  else
    if (wom > 2.0) then
      g_omega = sqrt(D*A/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
      c_omega = sqrt(A/D/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
    else
      if (wom == 0) then
        g_omega = 0.0
        c_omega = 0.0
      else
        g_omega = sqrt(D*A/rho)*sqrt(i*(wom**2.0)/8.0+(wom**4.0)/48.0)
        c_omega = sqrt(A/D/rho)*sqrt(i*(wom**2.0)/8.0+(wom**4.0)/48.0)
      end if
    end if
  end if
 
  ! Temporary function of omega_k. 
   if (omega_k /= 0) then
    kappa  = omega_k*l/c_omega
  else
    kappa   = 0.0
  end if

  Q0 = P0 / Z0

! ADD AN IF STATEMENT HERE FOR PRESSURE VS FLOW
  if (kappa == 0.0) then
    PL = P0 - (8.0*mu1*l*Q0)/(A*r**2.0)
  else
    PL = P0*cos(kappa) - i*Q0/(g_omega)*sin(kappa)!*rho*g*Lr/qc
  end if
    if (pq==1) then
        result = PL
    else ! Return flow at L
    !! MJC: Compute ZL
        if (g_omega /= 0.0) then
        ! Original
!          t1   = -i*Q0*sin(kappa)/g_omega + P0*cos(kappa)
!          t2   = Q0*cos(kappa) - i*g_omega*P0*sin(kappa)
!          ZL   = (t1/t2)!*qc/rho/g/Lr
        ! MJC 1
!            t1 = i/(g_omega*Z0) - cos(kappa)
!            t2 = g_omega*sin(kappa) - cos(kappa)/Z0
!            ZL = t1/t2
            ! MJC 2
            t1 = i/(g_omega)*sin(kappa) - cos(kappa)*Z0
            t2 = g_omega*cos(kappa)*Z0 - cos(kappa)
            ZL = t1/t2
        else
!            ZL  = (8.0*mu1*l/(A*r**2.0))!*qc/rho/g/Lr
            ZL  = Z0 - (8.0*mu1*l/(A*r**2.0))!*qc/rho/g/Lr
        end if
!    write(*,*)"Z0",Z0,"ZL",ZL
     QL = PL/ZL
    result = QL
    end if
end function PLfunc


!***************************************************************************
!* Function: comp_pres                                                     *
!* Version: 2.0                                                            *
!* By: NA Hill, Uni Glasgow                                                *
!* Edited: MJ Colebank, NC State University                                *
!***************************************************************************
function comp_pres (N,Omega,P_om,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,qc,g,ab,pq) result (output)
implicit none

  integer, intent(in)       :: N, ab, pq
  real(lng), intent(in)     :: Omega(:),trm_rst,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,r_min
  complex(lng), intent(in)  :: P_om(N)
  complex(lng)              :: temp, temp2, tempP(N), Z_om(N), P0, Z0, Pspecial(N)
  complex(lng)              :: P(0:Maxgen,1:N),Q(0:Maxgen,1:N),output(0:Maxgen,1:N)
  complex(lng)              :: Pin(0:Maxgen,1:N),Qin(0:Maxgen,1:N)
  integer                   :: k, j, kk
  real(lng)                 :: PspReal

! Compute the pressure at the root, i.e. the 0th generation
  do k = 1, N
    P(0,k) = P_om(k) !P0 at root of tree
  end do

! Initialize Pspecial to be the 0th generation
  do k = 2, N
    Pspecial(k-1) = P_om(k)
  end do
  Pspecial(N) = conjg(P_om(1))
! MJC
!do k = 1, N
!  Pspecial(k) = P_om(k)
!end do
!Pspecial(N) = conjg(P_om(1))

  ! Omega contains tmstps+1 frequencies because when computing the
  ! impedance it is easier to invert it when it is computed for all positive
  ! frequencies and since the interval goes from [-N/2/Tper:N/2/Tper-df],
  ! we include the frequency N/2/Tper in Omega and from this we get
  ! Z(-N/2/Tper) we then end up throwing pushing Z(N/2/Tper) out.
 
  Z_om = cmplx (0.0_lng,0.0_lng,lng) ! Initialize Z_om with zeros

  ! For all the positive frequencies compute the impedance.
  ! Since we know that the system is self-adjoint we don't
  ! have to compute the negative frequencies as well. 
  do k = N/2+1, N+1
    Computed  = 0.0_lng ! For each frequency make sure we don't carry
                        ! any values with us from the computations for the
                        ! previous frequency, so make sure Computed is zero.
    ! Since Z_om only has N places leave result for the frequency k at
    ! Z_om(k-1) we will later make up for this.
    Z_om(k-1) = Z0func (Omega(k),0,0,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)
!    write(*,*) "Impedance: ", Z_om(k-1), " at k=",k

    if (ab == 1) then
        do j = 0, Maxgen-1
            if (Computed(j,0) /= 0) then
                branches = j
                Z0 = Computed(j,0)!*qc/rho/g/Lr
                if (j == 0) then
                    P0 = Pspecial(k-1)
                else
                    P0 = P(j-1,k-1)
                end if
                Pin(j,k-1) = P0
                Qin(j,k-1) = P0/Z0
                P(j,k-1) = PLfunc(Omega(k),Z0,P0,j,0,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,1)
                Q(j,k-1) = PLfunc(Omega(k),Z0,P0,j,0,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,2)
                    if (k==N/2+1) then
                        write(*,*) 'P0',Pin(j,k-1),'Q0',Qin(j,k-1)
                        write(*,*) 'PL',P(j,k-1),'QL',Q(j,k-1)
!                        write(*,*) 'Z0',Z_om(k-1)
                    end if
            end if
        end do
    else if (ab == 2) then
        do j = 0, Maxgen-1
            if (Computed(0,j) /= 0) then
                branches = j
                Z0 = Computed(0,j)!*qc/rho/g/Lr
                if (j == 0) then
                    P0 = Pspecial(k-1)
                else
                    P0 = P(j-1,k-1)
                end if
                Pin(j,k-1) = P0
                Qin(j,k-1) = P0/Z0
                P(j,k-1) = PLfunc(Omega(k),Z0,P0,0,j,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,1)
                Q(j,k-1) = PLfunc(Omega(k),Z0,P0,0,j,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,2)
            end if
        end do
    end if
  end do

  temp = Z_om (N/2)
  
  ! Apply self-adjointness 
  Z_om(1:N/2)   = conjg(flipud(Z_om(N/2+1:N)))

  ! Shift the results for the positive frequencies one to the right
  ! to make up for the above.
  Z_om(N/2+1:N) = eoshift(Z_om(N/2+1:N),-1)    

  ! Insert Z(0,0) as described in the theoretical derivation
  Z_om(N/2+1) = temp

  open (unit=1, file='Z0om.2d', status = 'replace')
  do k = 1, N
    write(1,*) real(Z_om(k)), imag(Z_om(k))
  end do
  close (unit=1)

if (pq==1) then
  do j = 0, Maxgen-1
    do k = 1, N
      tempP(k) = P(j,k)
!      tempP(k) = Pin(j,k)
    end do
    
    temp2 = tempP(N/2)
    tempP(1:N/2)   = conjg(flipud(tempP(N/2+1:N)))
    tempP(N/2+1:N) = eoshift(tempP(N/2+1:N),-1)
    tempP(N/2+1)   = temp2

    do k= 1, N
      output(j,k) = tempP(k)
    end do
  end do

  open (unit=1, file='P0om.2d', status = 'replace')
  do k = 1, N
    write(1,*) real(output(0,k)), imag(output(0,k))
  end do
  close (unit=1)

  open (unit=1, file='P1om.2d', status = 'replace')
  do k = 1, N
    write(1,*) real(output(1,k)), imag(output(1,k))
  end do
  close (unit=1)
else
    do j = 0, Maxgen-1
      do k = 1, N
        tempP(k) = Q(j,k)
!        tempP(k) = Qin(j,k)
      end do
      
      temp2 = tempP(N/2)
      tempP(1:N/2)   = conjg(flipud(tempP(N/2+1:N)))
      tempP(N/2+1:N) = eoshift(tempP(N/2+1:N),-1)
      tempP(N/2+1)   = temp2

      do k= 1, N
        output(j,k) = tempP(k)
      end do
    end do

    open (unit=1, file='Q0om.2d', status = 'replace')
    do k = 1, N
      write(1,*) real(output(0,k)), imag(output(0,k))
    end do
    close (unit=1)

    open (unit=1, file='Q1om.2d', status = 'replace')
    do k = 1, N
      write(1,*) real(output(1,k)), imag(output(1,k))
    end do
    close (unit=1)
end if

end function comp_pres

end module tree_pres


