
subroutine boundary_condition_p
  use variables
  implicit none
  real(8) Tend
  ! Uc(1,0:Ny,1)=0.25d0
  !x方向境界条件
  !上下流端水位一定
  ! Uc(1:1,2:Ny-1,1)=0.1d0
  ! Uc(1:1,2:Ny-1,4)=0.3d0
  ! Uc(1:3,2:Ny-1,2)=0.430887d0
  select case(boundtype)
  case(1)
    !x方向上下流端の壁でuh=0
    Up(0,0:Ny+1,2)=0.d0
    Ep(0,0:Ny+1,1)=0.d0
    Ep(0,0:Ny+1,3)=0.d0
    Ep(0,0:Ny+1,4)=0.d0
    Fp(0,0:Ny+1,2)=0.d0
    Up(Nx+1,0:Ny+1,2)=0.d0
    Ep(Nx+1,0:Ny+1,1)=0.d0
    Ep(Nx+1,0:Ny+1,3)=0.d0
    Ep(Nx+1,0:Ny+1,4)=0.d0
    Fp(Nx+1,0:Ny+1,2)=0.d0
    !x方向側壁でuvhの差が0
    Fp(0:Nx+1,0,2)=Fp(0:Nx+1,1,2)
    Fp(0:Nx+1,Ny+1,2)=Fp(0:Nx+1,Ny,2)
    Ep(0:Nx+1,0,3)=Ep(0:Nx+1,1,3)
    Ep(0:Nx+1,Ny+1,3)=Ep(0:Nx+1,Ny,3)
    !y方向上下流端の壁でvh=0
    Up(0:Nx+1,0,3)=0.0d0
    Ep(0:Nx+1,0,3)=0.0d0
    Fp(0:Nx+1,0,1)=0.0d0
    Fp(0:Nx+1,0,2)=0.0d0
    Fp(0:Nx+1,0,4)=0.0d0
    Up(0:Nx+1,Ny+1,3)=0.0d0
    Ep(0:Nx+1,Ny+1,3)=0.0d0
    Fp(0:Nx+1,Ny+1,1)=0.0d0
    Fp(0:Nx+1,Ny+1,2)=0.0d0
    Fp(0:Nx+1,Ny+1,4)=0.0d0
    !y方向側壁でuvhの差が0
    Fp(0,0:Ny,2)=Fp(1,0:Ny,2)
    Fp(Nx+1,0:Ny,2)=Fp(Nx,0:Ny,2)
    Ep(0,0:Ny,3)=Ep(1,0:Ny,3)
    Ep(Nx+1,0:Ny,3)=Ep(Nx,0:Ny,3)
  case(2)
    !x方向上流端の壁でh,c=const
    Up(0,1:Ny,1)=0.1d0
    Up(0,1:Ny,4)=0.0d0
    !x方向上流端の壁でuh=0
    Up(0,0:Ny+1,2)=0.d0
    Ep(0,0:Ny+1,1)=0.d0
    Ep(0,0:Ny+1,3)=0.d0
    Ep(0,0:Ny+1,4)=0.d0
    Fp(0,0:Ny+1,2)=0.d0
    !x方向側壁でuvhの差が0
    Fp(0:Nx+1,0,2)=Fp(0:Nx+1,1,2)
    Fp(0:Nx+1,Ny+1,2)=Fp(0:Nx+1,Ny,2)
    Ep(0:Nx+1,0,3)=Ep(0:Nx+1,1,3)
    Ep(0:Nx+1,Ny+1,3)=Ep(0:Nx+1,Ny,3)
    !y方向側壁でuvhの差が0
    Fp(0,0:Ny+1,2)=Fp(1,0:Ny+1,2)
    Fp(Nx+1,0:Ny+1,2)=Fp(Nx,0:Ny+1,2)
    Ep(0,0:Ny+1,3)=Ep(1,0:Ny+1,3)
    Ep(Nx+1,0:Ny+1,3)=Ep(Nx,0:Ny+1,3)
    !y方向上下流端の壁でvh=0
    Up(0:Nx+1,0,3)=0.0d0
    Ep(0:Nx+1,0,3)=0.0d0
    Fp(0:Nx+1,0,1)=0.0d0
    Fp(0:Nx+1,0,2)=0.0d0
    Fp(0:Nx+1,0,4)=0.0d0
    Up(0:Nx+1,Ny+1,3)=0.0d0
    Ep(0:Nx+1,Ny+1,3)=0.0d0
    Fp(0:Nx+1,Ny+1,1)=0.0d0
    Fp(0:Nx+1,Ny+1,2)=0.0d0
    Fp(0:Nx+1,Ny+1,4)=0.0d0
  case(3)
    !x方向上流端の壁でh,c=const
    Up(0:1,1:Ny,1)=hup
    Up(0:1,1:Ny,2)=qup/Bini
    Up(0:1,1:Ny,3)=0.d0
    ! Up(1,1:Ny,2)=Up(2,1:Ny,2)
    Up(0,1:Ny,4)=hup*cup !土砂ありの場合
    Uc(0,1:Ny,5)=Zbini(0,1:Ny) !土砂ありの場合
    !下流端
    Up(Nx:Nx+1,1:Ny,1)=hdown
    Up(Nx:Nx+1,1:Ny,2)=qup/Bini
    Up(Nx:Nx+1,1:Ny,3)=0.d0

    ! Up(Nx+1,1:Ny,4)=hup*cup

    ! Up(Nx,1:Ny,1)=Up(Nx-1,1:Ny,1)
    ! Up(Nx,1:Ny,2)=Up(Nx-1,1:Ny,2)
    ! Up(Nx+1,1:Ny,1)=Up(Nx,1:Ny,1)
    ! Up(Nx+1,1:Ny,2)=Up(Nx,1:Ny,2)
    ! Up(Nx+1,1:Ny,3)=0.d0
    ! Up(Nx+1,1:Ny,4)=hup*cup
    
    !x方向上流端の壁でh,c=const
    ! Up(0,1:Ny,1)=hup
    ! Up(0,1:Ny,4)=hup*cup !土砂ありの場合
    ! Up(0,0:Ny+1,2)=qup/Bini
    ! Up(0,0:Ny+1,4)=hup*cup !土砂ありの場合
    ! Ep(0,0:Ny+1,2)=(Qup/Bini/hup)**2*hup+0.5d0*g*hup
    ! Fp(0,0:Ny+1,3)=0.5d0*g*hup
    !x方向上流端の壁でuh=0
    ! Ep(0,0:Ny+1,1)=qup/Bini
    ! Ep(0,0:Ny+1,3)=0.d0
    ! Ep(0,0:Ny+1,4)=0.d0
    ! Fp(0,0:Ny+1,2)=0.d0
    ! !x方向側壁でuvhの差が0
    ! Fp(0:Nx+1,0,2)=Fp(0:Nx+1,1,2)
    ! Fp(0:Nx+1,Ny+1,2)=Fp(0:Nx+1,Ny,2)
    ! Ep(0:Nx+1,0,3)=Ep(0:Nx+1,1,3)
    ! Ep(0:Nx+1,Ny+1,3)=Ep(0:Nx+1,Ny,3)
    ! !y方向側壁でuvhの差が0
    ! Fp(0,0:Ny+1,2)=Fp(1,0:Ny+1,2)
    ! Fp(Nx+1,0:Ny+1,2)=Fp(Nx,0:Ny+1,2)
    ! Ep(0,0:Ny+1,3)=Ep(1,0:Ny+1,3)
    ! Ep(Nx+1,0:Ny+1,3)=Ep(Nx,0:Ny+1,3)
    ! !y方向上下流端の壁でvh=0
    ! Up(0:Nx+1,0,3)=0.0d0
    ! Ep(0:Nx+1,0,3)=0.0d0
    ! Fp(0:Nx+1,0,1)=0.0d0
    ! Fp(0:Nx+1,0,2)=0.0d0
    ! Fp(0:Nx+1,0,4)=0.0d0
    ! Up(0:Nx+1,Ny+1,3)=0.0d0
    ! Ep(0:Nx+1,Ny+1,3)=0.0d0
    ! Fp(0:Nx+1,Ny+1,1)=0.0d0
    ! Fp(0:Nx+1,Ny+1,2)=0.0d0
    ! Fp(0:Nx+1,Ny+1,4)=0.0d0
  case(4)
    !x方向上流端の壁でh,c=const
    Up(0,1:Ny,1)=Up(1,1:Ny,1) !hup
    Up(0:1,1:Ny,2)=Up(0:1,1:Ny,2)!0.d0 !qup/Bini
    Up(0:1,1:Ny,3)=0.d0
    Up(0,1:Ny,4)=Up(1,1:Ny,4) !hup*cup !土砂ありの場合
    Up(0:1,1:Ny,5)=Zbini(0:1,1:Ny) !土砂ありの場合
    !下流端
    ! Up(Nx:Nx+1,1:Ny,1)=hdown
    Up(Nx+1,1:Ny,1:5)=Up(Nx,1:Ny,1:5)
    ! Up(Nx:Nx+1,1:Ny,3)=0.d0
  case(5)
        !x方向上流端の壁でh,c=const
    Tend=20.d0
    Up(0,1:Ny,1:5)=Up(1,1:Ny,1:5)
    ! Up(0,1:Ny,1)=0.d0 !Uc(1,1:Ny,1) !hup
    ! Up(0,1:Ny,2)=Up(0,1:Ny,2) !qup/Bini
    ! Up(0,1:Ny,3)=0.d0 !0.d0
    !下流端
    Up(Nx+1,1:Ny,1:5)=Up(Nx,1:Ny,1:5)
    if(T<Tend)then
      Up(2:3,46:50,1)=hup !0.04d0
      Up(2:3,46:50,2)=qup !1.d0/1000.d0
      Up(2:3,46:50,3)=0.d0
      Up(2:3,46:50,4)=hup*cup !0.04d0*0.3d0
      Up(2:3,46:50,5)=Zbini(2:3,46:50)
    else
      Up(2:3,46:50,1)=0.d0 !0.04d0
      Up(2:3,46:50,2)=0.d0 !1.d0/1000.d0
      Up(2:3,46:50,3)=0.d0
      Up(2:3,46:50,4)=0.d0 !0.04d0*0.3d0
      Up(2:3,46:50,5)=Zbini(2:3,46:50)
    endif
    case(6)
      Up(0,1:Ny,1:5)=Up(1,1:Ny,1:5)
      Up(Nx+1,1:Ny,1:5)=Up(Nx,1:Ny,1:5)
      Up(1:Nx,0,1:5)=Up(1:Nx,1,1:5)
      Up(1:Nx,Ny+1,1:5)=Up(1:Nx,Ny,1:5)
    case(7)
      if(rank_l(my_rank)<0)Up(0,1:Ny,1:5)=Up(1,1:Ny,1:5)
      if(rank_r(my_rank)<0)Up(Nx+1,1:Ny,1:5)=Up(Nx,1:Ny,1:5)
      if(rank_u(my_rank)<0)Up(1:Nx,0,1:5)=Up(1:Nx,1,1:5)
      if(rank_d(my_rank)<0)Up(1:Nx,Ny+1,1:5)=Up(1:Nx,Ny,1:5)
  end select
end subroutine boundary_condition_p