
subroutine boundary_condition_c
  use variables
  implicit none
  real(8) :: tend
  
  ! Uc(1,0:Ny,1)=0.25d0
  !x方向境界条件
  !上下流端水位一定
  ! Uc(1:1,2:Ny-1,1)=0.1d0
  ! Uc(1:1,2:Ny-1,4)=0.3d0
  ! Uc(1:3,2:Ny-1,2)=0.430887d0
  select case(boundtype)
  case(1)
    !x方向上下流端の壁でuh=0
    Uc(0,0:Ny+1,2)=0.d0
    Ec(0,0:Ny+1,1)=0.d0
    Ec(0,0:Ny+1,3)=0.d0
    Ec(0,0:Ny+1,4)=0.d0
    Fc(0,0:Ny+1,2)=0.d0
    Uc(Nx+1,0:Ny+1,2)=0.d0
    Ec(Nx+1,0:Ny+1,1)=0.d0
    Ec(Nx+1,0:Ny+1,3)=0.d0
    Ec(Nx+1,0:Ny+1,4)=0.d0
    Fc(Nx+1,0:Ny+1,2)=0.d0
    !x方向側壁でuvhの差が0
    Fc(0:Nx+1,0,2)=Fc(0:Nx+1,1,2)
    Fc(0:Nx+1,Ny+1,2)=Fc(0:Nx+1,Ny,2)
    Ec(0:Nx+1,0,3)=Ec(0:Nx+1,1,3)
    Ec(0:Nx+1,Ny+1,3)=Ec(0:Nx+1,Ny,3)
    !y方向上下流端の壁でvh=0
    Uc(0:Nx+1,0,3)=0.0d0
    Ec(0:Nx+1,0,3)=0.0d0
    Fc(0:Nx+1,0,1)=0.0d0
    Fc(0:Nx+1,0,2)=0.0d0
    Fc(0:Nx+1,0,4)=0.0d0
    Uc(0:Nx+1,Ny+1,3)=0.0d0
    Ec(0:Nx+1,Ny+1,3)=0.0d0
    Fc(0:Nx+1,Ny+1,1)=0.0d0
    Fc(0:Nx+1,Ny+1,2)=0.0d0
    Fc(0:Nx+1,Ny+1,4)=0.0d0
    !y方向側壁でuvhの差が0
    Fc(0,0:Ny,2)=Fc(1,0:Ny,2)
    Fc(Nx+1,0:Ny,2)=Fc(Nx,0:Ny,2)
    Ec(0,0:Ny,3)=Ec(1,0:Ny,3)
    Ec(Nx+1,0:Ny,3)=Ec(Nx,0:Ny,3)
  case(2)
    !x方向上流端の壁でh,c=const
    Uc(0,1:Ny,1)=0.1d0
    Uc(0,1:Ny,4)=0.0d0
    !x方向上流端の壁でuh=0
    Uc(0,0:Ny+1,2)=0.d0
    Ec(0,0:Ny+1,1)=0.d0
    Ec(0,0:Ny+1,3)=0.d0
    Ec(0,0:Ny+1,4)=0.d0
    Fc(0,0:Ny+1,2)=0.d0
    !x方向側壁でuvhの差が0
    Fc(0:Nx+1,0,2)=Fc(0:Nx+1,1,2)
    Fc(0:Nx+1,Ny+1,2)=Fc(0:Nx+1,Ny,2)
    Ec(0:Nx+1,0,3)=Ec(0:Nx+1,1,3)
    Ec(0:Nx+1,Ny+1,3)=Ec(0:Nx+1,Ny,3)
    !y方向側壁でuvhの差が0
    Fc(0,0:Ny+1,2)=Fc(1,0:Ny+1,2)
    Fc(Nx+1,0:Ny+1,2)=Fc(Nx,0:Ny+1,2)
    Ec(0,0:Ny+1,3)=Ec(1,0:Ny+1,3)
    Ec(Nx+1,0:Ny+1,3)=Ec(Nx,0:Ny+1,3)
    !y方向上下流端の壁でvh=0
    Uc(0:Nx+1,0,3)=0.0d0
    Ec(0:Nx+1,0,3)=0.0d0
    Fc(0:Nx+1,0,1)=0.0d0
    Fc(0:Nx+1,0,2)=0.0d0
    Fc(0:Nx+1,0,4)=0.0d0
    Uc(0:Nx+1,Ny+1,3)=0.0d0
    Ec(0:Nx+1,Ny+1,3)=0.0d0
    Fc(0:Nx+1,Ny+1,1)=0.0d0
    Fc(0:Nx+1,Ny+1,2)=0.0d0
    Fc(0:Nx+1,Ny+1,4)=0.0d0
  case(3)
    !x方向上流端の壁でh,c=const
    Uc(0:1,1:Ny,1)=hup
    Uc(0:1,1:Ny,2)=qup/Bini
    ! Uc(1,1:Ny,2)=Uc(2,1:Ny,2)
    Uc(0:1,1:Ny,3)=0.d0
    Uc(0,1:Ny,4)=hup*cup !土砂ありの場合
    Uc(0,1:Ny,5)=Zbini(0,1:Ny) !土砂ありの場合
    !下流端
    Uc(Nx:Nx+1,1:Ny,1)=hdown
    Uc(Nx:Nx+1,1:Ny,2)=qup/Bini
    Uc(Nx:Nx+1,1:Ny,3)=0.d0
    ! Uc(Nx+1,1:Ny,3)=hup*cup
    ! Uc(Nx+1,1:Ny,3)=Zbini(Nx+1,1:Ny)
    ! Uc(Nx+1,1:Ny,4)=hup*cup

    ! Uc(Nx,1:Ny,1)=Uc(Nx-1,1:Ny,1)
    ! Uc(Nx,1:Ny,2)=Uc(Nx-1,1:Ny,2)
    ! Uc(Nx+1,1:Ny,1)=Uc(Nx,1:Ny,1)
    ! Uc(Nx+1,1:Ny,2)=Uc(Nx,1:Ny,2)
    ! Uc(Nx+1,1:Ny,3)=0.d0
    ! Uc(Nx+1,1:Ny,4)=hup*cup

    ! Uc(Nx+1,1:Ny,1)=hup
    ! Uc(Nx+1,1:Ny,2)=qup/Bini
    ! Uc(Nx+1,1:Ny,3)=0.d0
    
    !両岸の境界条件は，差分時に使う？
    ! Uc(0,0:Ny+1,4)=hup*cup !土砂ありの場合
    ! Ec(0,0:Ny+1,2)=(Qup/Bini/hup)**2*hup+0.5d0*g*hup
    ! Fc(0,0:Ny+1,3)=0.5d0*g*hup
    !x方向上流端の壁でuh=0
    ! Ec(0,0:Ny+1,1)=qup/Bini
    ! Ec(0,0:Ny+1,3)=0.d0
    ! Ec(0,0:Ny+1,4)=0.d0
    ! Fc(0,0:Ny+1,2)=0.d0
    !x方向側壁でuvhの差が0
    ! Fc(0:Nx+1,0,2)=Fc(0:Nx+1,1,2)
    ! Fc(0:Nx+1,Ny+1,2)=Fc(0:Nx+1,Ny,2)
    ! Ec(0:Nx+1,0,3)=Ec(0:Nx+1,1,3)
    ! Ec(0:Nx+1,Ny+1,3)=Ec(0:Nx+1,Ny,3)
    !y方向側壁でuvhの差が0
    ! Fc(0,0:Ny+1,2)=Fc(1,0:Ny+1,2)
    ! Fc(Nx+1,0:Ny+1,2)=Fc(Nx,0:Ny+1,2)
    ! Ec(0,0:Ny+1,3)=Ec(1,0:Ny+1,3)
    ! Ec(Nx+1,0:Ny+1,3)=Ec(Nx,0:Ny+1,3)
    !y方向上下流端の壁でvh=0
    ! Uc(0:Nx+1,0,3)=0.0d0
    ! Ec(0:Nx+1,0,3)=0.0d0
    ! Fc(0:Nx+1,0,1)=0.0d0
    ! Fc(0:Nx+1,0,2)=0.0d0
    ! Fc(0:Nx+1,0,4)=0.0d0
    ! Uc(0:Nx+1,Ny+1,3)=0.0d0
    ! Ec(0:Nx+1,Ny+1,3)=0.0d0
    ! Fc(0:Nx+1,Ny+1,1)=0.0d0
    ! Fc(0:Nx+1,Ny+1,2)=0.0d0
    ! Fc(0:Nx+1,Ny+1,4)=0.0d0
    !x方向上流端の壁でh,c=const,uh=0
    ! Uc(0,0:Ny+1,1)=hup
    ! Uc(0,0:Ny+1,2)=qup/Bini
    ! Uc(0,0:Ny+1,3)=0.0d0
    ! Uc(0,0:Ny+1,4)=hup*cup !土砂ありの場合
    ! Uc(0,0:Ny+1,5)=Uc(0,0:Ny+1,5)
    ! Ec(0,0:Ny+1,1)=Qup/Bini
    ! Ec(0,0:Ny+1,2)=(Qup/Bini/hup)**2*hup+0.5d0*g*hup
    ! Ec(0,0:Ny+1,3)=0.d0
    ! Ec(0,0:Ny+1,4)=Qup/Bini*Cup
    ! Ec(0,0:Ny+1,5)=0.d0
    ! Fc(0,0:Ny+1,1)=0.d0
    ! Fc(0,0:Ny+1,2)=0.d0
    ! Fc(0,0:Ny+1,3)=0.5d0*g*hup
    ! Fc(0,0:Ny+1,4)=0.d0
    ! Fc(0,0:Ny+1,5)=0.d0
    ! !x方向側壁(左岸)でuvhの差が0
    ! Fc(0:Nx+1,0,2)=Fc(0:Nx+1,1,2)
    ! Fc(0:Nx+1,Ny+1,2)=Fc(0:Nx+1,Ny,2)
    ! Ec(0:Nx+1,0,3)=Ec(0:Nx+1,1,3)
    ! Ec(0:Nx+1,Ny+1,3)=Ec(0:Nx+1,Ny,3)
    ! !y方向側壁でuvhの差が0
    ! Fc(0,0:Ny+1,2)=Fc(1,0:Ny+1,2)
    ! Fc(Nx+1,0:Ny+1,2)=Fc(Nx,0:Ny+1,2)
    ! Ec(0,0:Ny+1,3)=Ec(1,0:Ny+1,3)
    ! Ec(Nx+1,0:Ny+1,3)=Ec(Nx,0:Ny+1,3)
    ! !y方向上下流端の壁でvh=0
    ! Uc(0:Nx+1,0,3)=0.0d0
    ! Ec(0:Nx+1,0,3)=0.0d0
    ! Fc(0:Nx+1,0,1)=0.0d0
    ! Fc(0:Nx+1,0,2)=0.0d0
    ! Fc(0:Nx+1,0,4)=0.0d0
    ! Uc(0:Nx+1,Ny+1,3)=0.0d0
    ! Ec(0:Nx+1,Ny+1,3)=0.0d0
    ! Fc(0:Nx+1,Ny+1,1)=0.0d0
    ! Fc(0:Nx+1,Ny+1,2)=0.0d0
    ! Fc(0:Nx+1,Ny+1,4)=0.0d0
  case(4)
    !x方向上流端の壁でh,c=const
    Uc(0,1:Ny,1)=Uc(1,1:Ny,1) !hup
    Uc(0,1:Ny,2)=Uc(1,1:Ny,2) !qup/Bini
    Uc(0:1,1:Ny,3)=0.d0 !0.d0
    Uc(0,1:Ny,4)=Uc(1,1:Ny,4) !hup*cup !土砂ありの場合
    Uc(0:1,1:Ny,5)=Zbini(0:1,1:Ny) !土砂ありの場合
    !下流端
    Uc(Nx+1,1:Ny,1:5)=Uc(Nx,1:Ny,1:5)
  case(5)
    !x方向上流端の壁でh,c=const
    Tend=20.d0
    Uc(0,1:Ny,1:5)=Uc(1,1:Ny,1:5)
    ! Uc(0,1:Ny,1)=0.d0 !Uc(1,1:Ny,1) !hup
    ! Uc(0,1:Ny,2)=Uc(0,1:Ny,2) !qup/Bini
    ! Uc(0,1:Ny,3)=0.d0 !0.d0
    !下流端
    Uc(Nx+1,1:Ny,1:5)=Uc(Nx,1:Ny,1:5)
    if(T<Tend)then
      Uc(2:3,46:50,1)=hup !0.04d0
      Uc(2:3,46:50,2)=qup !1.d0/1000.d0
      Uc(2:3,46:50,3)=0.d0
      Uc(2:3,46:50,4)=hup*cup !0.04d0*0.3d0
      Uc(2:3,46:50,5)=Zbini(2:3,46:50)
    else
      Uc(2:3,46:50,1)=0.d0 !0.04d0
      Uc(2:3,46:50,2)=0.d0 !1.d0/1000.d0
      Uc(2:3,46:50,3)=0.d0
      Uc(2:3,46:50,4)=0.d0 !0.04d0*0.3d0
      Uc(2:3,46:50,5)=Zbini(2:3,46:50)
    endif
  case(6)
    Uc(0,1:Ny,1:5)=Uc(1,1:Ny,1:5)
    Uc(Nx+1,1:Ny,1:5)=Uc(Nx,1:Ny,1:5)
    Uc(1:Nx,0,1:5)=Uc(1:Nx,1,1:5)
    Uc(1:Nx,Ny+1,1:5)=Uc(1:Nx,Ny,1:5)
  case(7)
    if(rank_l(my_rank)<0)Uc(0,1:Ny,1:5)=Uc(1,1:Ny,1:5)
    if(rank_r(my_rank)<0)Uc(Nx+1,1:Ny,1:5)=Uc(Nx,1:Ny,1:5)
    if(rank_u(my_rank)<0)Uc(1:Nx,0,1:5)=Uc(1:Nx,1,1:5)
    if(rank_d(my_rank)<0)Uc(1:Nx,Ny+1,1:5)=Uc(1:Nx,Ny,1:5)
  end select


end subroutine boundary_condition_c
