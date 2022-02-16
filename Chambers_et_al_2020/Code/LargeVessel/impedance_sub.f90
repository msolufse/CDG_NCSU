subroutine impedance_driver(tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g,alpha_in,beta_in,lrr1,lrr2,term_ID)
  use f90_tools
  use root_imp 
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g
  real(lng), intent(in)  :: alpha_in,beta_in,lrr1,lrr2
  real(lng), intent(out) :: y_xt(tmstps)
  integer, intent(in)    :: term_ID


  call impedance (tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g,alpha_in,beta_in,lrr1,lrr2, term_ID)
end subroutine impedance_driver
