subroutine impedance_driver(tmstps,Period,r_root,r_min,y_xt,R1,R2,CT,LT)
  use f90_tools
  use root_imp 
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,r_root,r_min,R1,R2,CT,LT
  real(lng), intent(out) :: y_xt(tmstps)
  !write(6,*) 'R1=', R1,'R2=',R2,'CT=',CT,'LT=',LT
  call impedance (tmstps,Period,r_root,r_min,y_xt,R1,R2,CT,LT)
end subroutine impedance_driver
