subroutine external_driver(tmstps,Period,rho,mu,r_root,r_min,yE11,yE12,yE21,yE22,&
                            Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,&
                            asym,expo,lrrA,lrrV)
  use f90_tools
  use ext_new_match
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,rho,mu,r_root,r_min,Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV
  real(lng), intent(out) :: yE11(tmstps),yE12(tmstps),yE21(tmstps),yE22(tmstps)
  integer :: j
!write(*,*)'Impedancedriver started'    ! JAM
!write(*,*)'timesteps',tmstps    ! JAM
  do j = 1, tmstps
    yE11(j) = 0.0
    yE12(j) = 0.0
    yE21(j) = 0.0
    yE22(j) = 0.0
!write(*,*)j,' Impedance initiated'    ! JAM
  end do
!write(*,*)'New match called'    ! JAM
  call external (tmstps,Period,rho,mu,r_root,r_min,yE11,yE12,yE21,yE22,&
                  Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV)
end subroutine external_driver
