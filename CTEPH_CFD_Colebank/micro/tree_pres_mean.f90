module tree_pres_mean
use f90_tools
implicit none

! Define these values apriori
integer, parameter:: Maxgen = 100
integer, parameter:: MaxBranches = 100000
integer, parameter:: tmstps = 1024!8192
!real(lng), parameter:: Period = 0.85

  real(lng), save :: Period, trm_rst, ff1, ff2, ff3, r_root, r_min, alpha, beta, lrr!, tmstps
!  integer, save   :: Maxgen, tmstps, maxStat, Ctmstps
  complex(lng), allocatable :: Computed(:,:)
  real(lng), allocatable :: p0_network(:,:),q0_network(:,:)
  real(lng), allocatable :: pL_network(:,:),qL_network(:,:)
!  integer, allocatable :: jbranches(:)


real(lng), parameter :: Hct = 0.45

complex(lng), allocatable   :: Z_om(:)

! integer localmax, branchesj(0:Maxgen), branches
integer localmax, branches!, index_in

contains


!***************************************************************************
!* Function: Z0func                                                        *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            *
!***************************************************************************
recursive function Z0func (alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g) result (Z0)
implicit none

  integer, intent(in)   :: alpha_pow, beta_pow
  real(lng), intent(in) :: ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g
  ! ADDED BY MJC
  real(lng)    :: alpha_update, beta_update
  integer      :: generations
  real(lng)    :: nu, r,diameter,eta_D,mu_D,mu1, r_d, l, A, A_d, D, wom
  real(lng)    ::  Hct_term, C_exp !,Hct
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
  ! Temporary functions of r.
  
  ! Determine the resistance of the root of the terminal branches.
  ! if (generations >= Maxgen) 
  if (r <= r_min) then
    if (generations >= localmax) then
      localmax = generations
    end if
    if (generations >= Maxgen) then
      write (*,*) 'Generations level exceeded'
    end if
        ZL = trm_rst!*q/rho/g/Lr
  else
    ! Get Z0 recursively at reduced frequencies.
    if (abs(Computed(alpha_pow+1, beta_pow)) /= 0.0) then
      Zl_0 = Computed(alpha_pow+1, beta_pow)
    else  
      Zl_0 = Z0func (alpha_pow+1,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)
    end if

    if (abs(Computed(alpha_pow, beta_pow+1)) /= 0.0) then
      Zr_0 = Computed(alpha_pow, beta_pow+1)
    else
      Zr_0 = Z0func (alpha_pow,beta_pow+1,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)
    end if

!    ZL = 1.0/(1.0/Zl_0 + 1.0/Zr_0)
    ZL = (Zl_0*Zr_0)/(Zl_0+Zr_0)

  end if
 
  ! Compute pressure and flow in the network
  Z0  = (8.0*mu*l/(A*r**2.0) + ZL)
!  write(*,*) 'Z0',Z0, 'alpha',alpha_pow, 'beta',beta_pow

  ! This should be non-dimensional

  Computed(alpha_pow,beta_pow) = Z0!*qc/rho/g/Lr
!    write(*,*) 'Z0',Z0,alpha_pow,beta_pow
end function Z0func

!==============================================================================================
!====================function computes the pressure at the end of each vessel==================
function PLfunc (Z0,P0,Q0,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,pq) result (result)

  integer, intent(in)   :: alpha_pow, beta_pow, pq
  real(lng), intent(in) :: ff1,ff2,ff3,rho,mu,r_root
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
    endif

!l = 12.4*(r**1.10) ! from Olufsen 2012

diameter = 2.0*r*(10**4)
    eta_D = 3.2+6*exp(-0.085*diameter)-2.44*exp(-0.06*(diameter**0.645))
    ! Terms that depend on hematocrit
    C_exp = (0.8+exp(-0.075*diameter))*(1.0/(1.0+10**(-11.0) * diameter**12.0) - 1.0) + (1.0/(1.0+10**(-11.0) * diameter**12.0))
!    Hct = 0.45
    Hct_term = ((1.0 - Hct)**C_exp - 1.0)/((1.0-0.45)**C_exp - 1.0)
    mu_D  = (1.0 + (eta_D - 1.0)*Hct_term*(diameter/(diameter - 1.1))**2.0)*(diameter/(diameter - 1.1))**2.0

    mu1 = mu*mu_D/3.2


  nu   = mu1/rho                ! Kinematic blood viscosity, m^2/s.
  D    = 1.0/(ff1*exp(ff2*r_d)+ff3)*3.0*A_d/2.0 ! Distensibility.


  Q0 = P0 / Z0

! ADD AN IF STATEMENT HERE FOR PRESSURE VS FLOW
    PL = P0 - (8.0*mu1*l*Q0)/(A*r**2.0)
    if (pq==1) then
        result = PL
    else ! Return flow at L
    !! MJC: Compute ZL
     ZL  = Z0 - (8.0*mu1*l/(A*r**2.0))!*qc/rho/g/Lr
     QL = PL/ZL
    result = QL
    end if
end function PLfunc

!==============================================================================================
!====================function computes the pressure and flow at every branch==================
recursive function PQ_network (Z_in,P0,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,index_in) result (result)

  integer, intent(in)   :: alpha_pow, beta_pow, index_in
  real(lng), intent(in) :: ff1,ff2,ff3,rho,mu,r_root
complex(lng), intent(in) :: P0, Z_in(0:Maxgen,0:Maxgen)
  
! ADDED BY MJC
  real(lng)    :: alpha_update, beta_update
  real(lng)    :: nu, r,diameter,eta_D,mu_D,mu1, r_d, l, A, A_d, D, wom
 real(lng)    :: Hct_term, C_exp!,Hct
 integer   :: index

!  real(lng)    :: beta, alpha
  complex(lng) :: i, g_omega, c_omega, kappa, Q0,Z0, PL,QL,result
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

    mu1 = mu*mu_D/3.2

    ! Start recurive function call
    ! If Z0 is not defined, don't call loop
    if (Computed(alpha_pow,beta_pow)/=0) then
    !First, calculate the flow at the inlet
    Z0 = Computed(alpha_pow,beta_pow)
    Q0 = P0 / Z0
!    write(*,*) 'Z0',Z0,'Q0',Q0
    ! Now that we know pressure and flow at the inlet, compute P/Q at L
    PL = PLfunc (Z0,P0,Q0,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,1)
    QL = PLfunc (Z0,P0,Q0,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,2)

    ! Now recursively go down the network and update all appropriate vectors/matrices
    index = index_in
    p0_network(index,1) = index
    p0_network(index,2) = alpha_pow
    p0_network(index,3) = beta_pow
    p0_network(index,4) = r_d
    p0_network(index,5) = P0

    pL_network(index,1) = index
    pL_network(index,2) = alpha_pow
    pL_network(index,3) = beta_pow
    pL_network(index,4) = r_d
    pL_network(index,5) = PL

    q0_network(index,1) = index
    q0_network(index,2) = alpha_pow
    q0_network(index,3) = beta_pow
    q0_network(index,4) = r_d
    q0_network(index,5) = Q0

    qL_network(index,1) = index
    qL_network(index,2) = alpha_pow
    qL_network(index,3) = beta_pow
    qL_network(index,4) = r_d
    qL_network(index,5) = QL

!    write(*,*) 'index',index,'Z0',Z_in(alpha_pow,beta_pow)
    if (r_d>=r_min) then
        
        index = index+1
        index = PQ_network (Z_in,PL,alpha_pow+1,beta_pow,ff1,ff2,ff3,rho,mu,r_root,index)

        index = index+1
        index = PQ_network (Z_in,PL,alpha_pow,beta_pow+1,ff1,ff2,ff3,rho,mu,r_root,index)
    endif

endif
result = index
end function PQ_network


!***************************************************************************
!* Function: comp_pres                                                     *
!* Version: 2.0                                                            *
!* By: NA Hill, Uni Glasgow                                                *
!* Edited: MJ Colebank, NC State University                                *
!***************************************************************************
function comp_pres (P_om,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,qc,g) result (output)
implicit none

  real(lng), intent(in)     :: trm_rst,ff1,ff2,ff3,rho,mu,r_root,Lr,qc,g,r_min
  complex(lng)              :: P_om, Z_om(1:MaxBranches,5)
  complex(lng)              :: P0, Z0
  integer                   :: output
  integer                   :: k, j, kk, index_in
  real(lng)                 :: PspReal
  character (len=30)        :: fn


  ! Omega contains tmstps+1 frequencies because when computing the
  ! impedance it is easier to invert it when it is computed for all positive
  ! frequencies and since the interval goes from [-N/2/Tper:N/2/Tper-df],
  ! we include the frequency N/2/Tper in Omega and from this we get
  ! Z(-N/2/Tper) we then end up throwing pushing Z(N/2/Tper) out.
 
  Z_om = cmplx (0.0_lng,0.0_lng,lng) ! Initialize Z_om with zeros

  ! For all the positive frequencies compute the impedance.
  ! Since we know that the system is self-adjoint we don't
  ! have to compute the negative frequencies as well.

    ! This will generate all impedance values
    index_in = 0
    P0 = P_om
    Z_om = Z0func(0,0,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,qc,g)

    index_in = PQ_network(Z_om,P0,0,0,ff1,ff2,ff3,rho,mu,r_root,index_in)


! Write everything to file
write(*,*) 'Program finished: number of branches ', index_in

fn = 'p0_tree'
open (unit=21, file=fn, status='replace')
do j = 0, index_in
    write (21,*) p0_network(j,1), p0_network(j,2), p0_network(j,3), &
                p0_network(j,4), p0_network(j,5)/1333.22 + 2.0
end do
close(unit=21)

fn = 'q0_tree'
open (unit=22, file=fn, status='replace')
do j = 0, index_in
    write (22,*) q0_network(j,1), q0_network(j,2), q0_network(j,3), &
                 q0_network(j,4), q0_network(j,5)*10.0
end do
close(unit=22)


fn = 'pL_tree'
open (unit=23, file=fn, status='replace')
do j = 0, index_in
    write (23,*) pL_network(j,1), pL_network(j,2), pL_network(j,3), &
                pL_network(j,4), pL_network(j,5)/1333.22 + 2.0
end do
close(unit=23)

fn = 'qL_tree'
open (unit=24, file=fn, status='replace')
do j = 0, index_in
    write (24,*) qL_network(j,1), qL_network(j,2), qL_network(j,3), &
                 qL_network(j,4), qL_network(j,5)*10.0
end do
close(unit=24)

end function comp_pres

end module tree_pres_mean


