
subroutine Initial_condition
  use variables
  implicit none
  real(8) :: kx,ky,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf!,ustar0
  character :: dummy
  h=0.d0
  ca=0.d0
  u=0.d0
  v=0.d0
  iero=0.d0
  Uc=0.d0
  Ec=0.d0
  Fc=0.d0
  Cc=0.d0
  Qxc=0.d0
  Qyc=0.d0
  Up=0.d0
  Ep=0.d0
  Fp=0.d0
  Cp=0.d0
  Qxp=0.d0
  Qyp=0.d0
  Zb=100.d0

  select case (inittype)
  case(1)
    DO I=1,Nx+1
      if(I<=Nx/2)then
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)  !勾配0
      else
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)  !勾配0
      endif
    ENDDO
    h(0:int(Nx/2),0:Ny+1)=0.5d0
    h(int(Nx/2)+1:Nx+1,0:Ny+1)=0.05d0
    Ca(0:Nx+1,0:Ny+1)=0.d0
  case(2)
    DO I=1,Nx+1
      if(I<=Nx/2)then
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx  !勾配sliに設定
      else
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx*0.01d0!*0.01d0  !勾配sliに設定
      endif
    ENDDO
    h(1:Nx,1:Ny)=hup
    Ca(1:Nx,1:Ny)=Cup
    u(1:Nx,1:Ny)=Qup/Bini/hup
    v(1:Nx,1:Ny)=0.d0
  case(3)
    DO I=1,Nx+1
      ! Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx  !勾配sliに設定
      if(I<=Nx/2)then
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx  !勾配sliに設定
      else
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx*0.01d0!*0.01d0  !勾配sliに設定
      endif
    ENDDO
    zbini=zb
    ! zb(101:105,11:20)=zb(101:105,11:20)+0.5d0
    h(1:Nx,1:Ny)=hup
    ! h(1:int(Nx/2),1:Ny)=hup
    ! h(int(Nx/2)+1:Nx,1:Ny)=hdown
    ! h(101:105,11:20)=0.d0
    Ca(0:Nx+1,0:Ny+1)=Cup
    u(1:Nx,1:Ny)=Qup/Bini/hup
    ! u(100:106,11:20)=0.d0
    ! u(1:int(Nx/2),1:Ny)=Qup/Bini/hup
    ! u(int(Nx/2)+1:Nx,1:Ny)=Qup/Bini/hdown
    v(1:Nx,1:Ny)=0.d0
    write(*,*) "u(200,10)=",u(200,10),qup/Bini/hup
    zbmin(1:Nx,1:Ny)=0.d0
  case(4)
    DO I=1,Nx+1
      if(I<=Nx/2)then
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-sli*dx  !勾配sliに設定
      else
        Zb(I,0:Ny+1)=Zb(I-1,0:Ny+1)-slidown*dx!*0.01d0  !勾配sliに設定
      endif
    ENDDO
    zbini=zb
    h(1:Nx,1:Ny)=0.d0
    Ca(0:Nx+1,0:Ny+1)=0.d0
    h(10:20,1:Ny)=hini
    Ca(10:20,1:Ny)=cini
    v(1:Nx,1:Ny)=0.d0
    u(1:Nx,1:Ny)=0.d0
    ! h(1:int(Nx/2),1:Ny)=hup
    ! h(int(Nx/2)+1:Nx,1:Ny)=hdown
    ! u(1:int(Nx/2),1:Ny)=Qup/Bini/hup
    ! u(int(Nx/2)+1:Nx,1:Ny)=Qup/Bini/hdown
    zbmin(1:Nx,1:Ny)=0.d0
    write(*,*) "u(200,10)=",u(200,10)
  case(5)
    DO J=1,Ny
      read(1001,*) (Zb(I,J),I=1,Nx)
    ENDDO
    Zb(0,1:Ny)=Zb(1,1:Ny)
    Zb(Nx+1,1:Ny)=Zb(Nx,1:Ny)
    Zb(0:Nx,0)=Zb(0:Nx,1)
    Zb(0:Nx,Ny+1)=Zb(0:Nx,Ny)
    zbini=zb
    ! zbmin=zb
    zbmin(1:275,1:Ny)=Zb
    zbmin(276:Nx,1:Ny)=0.d0
    ! zbmin(5:155,1:Ny)=zb(5:155,1:Ny)-0.1d0
    h=0.d0
    Ca=0.d0
    v=0.d0
    u=0.d0
  case(6)
    DO J=1,Ny
    ! DO J=Ny,1,-1
      read(1001,*) (Zb(I,J),I=1,Nx)
    ENDDO
    DO I=1,5
      read(1002,*) dummy
    ENDDO
    h=0.d0
    Ca=0.d0
    DO J=1,Ny
    ! DO J=Ny,1,-1
      read(1002,*) (h(I,J),I=1,Nx)
      DO I=1,Nx
        if(h(I,J)>hcr)then
          Ca(I,J)=Cini
        ! else
        !   h(I,J)=hcr
        ENDIF
      ENDDO
    ENDDO
    v=0.d0
    u=0.d0
    Zb(0,1:Ny)=Zb(1,1:Ny)
    Zb(Nx+1,1:Ny)=Zb(Nx,1:Ny)
    Zb(0:Nx,0)=Zb(0:Nx,1)
    Zb(0:Nx,Ny+1)=Zb(0:Nx,Ny)
    zbini=zb
    ! zbmin=zb
    zbmin=Zb-1.d0
  end select

  ! call dry_or_wet(h,zb,u,v,iswet,hcr,Nx,Ny,1)
  call dry_or_wet_f(h,zb,u,v,iswet,hcr,Nx,Ny)
  ! where(h>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  !calc U
  Uc(0:Nx+1,0:Ny+1,1)=h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,2)=u(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,3)=v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,4)=Ca(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,5)=Zb(0:Nx+1,0:Ny+1)
  !calc E
  Ec(0:Nx+1,0:Ny+1,1)=u(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Ec(0:Nx+1,0:Ny+1,2)=u(0:Nx+1,0:Ny+1)**2*h(0:Nx+1,0:Ny+1)+0.5d0*g*h(0:Nx+1,0:Ny+1)**2
  Ec(0:Nx+1,0:Ny+1,3)=u(0:Nx+1,0:Ny+1)*v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Ec(0:Nx+1,0:Ny+1,4)=ca(0:Nx+1,0:Ny+1)*u(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Ec(0:Nx+1,0:Ny+1,5)=0.d0
  !calc F
  Fc(0:Nx+1,0:Ny+1,1)=v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Fc(0:Nx+1,0:Ny+1,2)=u(0:Nx+1,0:Ny+1)*v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Fc(0:Nx+1,0:Ny+1,3)=v(0:Nx+1,0:Ny+1)**2*h(0:Nx+1,0:Ny+1)+0.5d0*g*h(0:Nx+1,0:Ny+1)**2
  Fc(0:Nx+1,0:Ny+1,4)=ca(0:Nx+1,0:Ny+1)*v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Fc(0:Nx+1,0:Ny+1,5)=0.d0


  DO J=1,Ny
    !$omp parallel do default(shared),private(I,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf)
    DO I=1,Nx
      Ec(I,J,1)=u(I,J)*h(I,J)
      Ec(I,J,2)=u(I,J)**2*h(I,J)+0.5d0*g*h(I,J)**2
      Ec(I,J,3)=u(I,J)*v(I,J)*h(I,J)
      Ec(I,J,4)=ca(I,J)*u(I,J)*h(I,J)
      Ec(I,J,5)=0.d0
      Fc(I,J,1)=v(I,J)*h(I,J)
      Fc(I,J,2)=u(I,J)*v(I,J)*h(I,J)
      Fc(I,J,3)=v(I,J)**2*h(I,J)+0.5d0*g*h(I,J)**2
      Fc(I,J,4)=ca(I,J)*v(I,J)*h(I,J)
      Fc(I,J,5)=0.d0
      if(iswet(I,J)==1)then
        if(ca(I,J)>=0.4d0*cstar.and.h(I,J)/dm<30.d0)then
          Sfx(I,J)=u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
          /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
          Sfy(I,J)=v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
          /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
        elseif(ca(I,J)>0.01d0.and.h(I,J)/dm<30.d0)then
          Sfx(I,J)=1.d0/0.49d0*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
          Sfy(I,J)=1.d0/0.49d0*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
        else
          Sfx(I,J)=nm**2*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
          Sfy(I,J)=nm**2*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
        endif
      else
        Sfx(I,J)=0.d0
        Sfy(I,J)=0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I+1,J)==1)then
        S0x(I,J)=-(Zb(I+1,J)-Zb(I,J))/Dx
      else
        S0x(I,J)=0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I,J+1)==1)then
        S0y(I,J)=-(Zb(I,J+1)-Zb(I,J))/Dy !collectorなので前進差分
      else
        S0y(I,J)=0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)
      if(iswet(I+1,J)==1.and.iswet(I-1,J)==1)then
        uhxx=(eps*(u(I+1,J)*h(I+1,J)-u(I,J)*h(I,J))/dx-eps*(u(I,J)*h(I,J)-u(I-1,J)*h(I-1,J))/dx)/dx
      else
        uhxx=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J-1)==1)then
        uhyy=(eps*(u(I,J+1)*h(I,J+1)-u(I,J)*h(I,J))/dy-eps*(u(I,J)*h(I,J)-u(I,J-1)*h(I,J-1))/dy)/dy
      else
        uhyy=0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)
      if(iswet(I+1,J)==1.and.iswet(I-1,J)==1)then
        vhyy=(eps*(v(I,J+1)*h(I,J+1)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I,J-1)*h(I,J-1))/dy)/dy
      else
        vhyy=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J-1)==1)then
        vhxx=(eps*(v(I+1,J)*h(I+1,J)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I-1,J)*h(I-1,J))/dx)/dx
      else
        vhxx=0.d0
      endif
      !calc water surface gradient!要検討
      if(u(I,J)>0.d0)then
        thewx=(-Zb(I+1,J)-H(I+1,J)+Zb(I,J)+H(I,J))/dx
      elseif(u(I,J)<0.d0)then
        thewx=(-Zb(I,J)-H(I,J)+Zb(I-1,J)+H(I-1,J))/dx
      else
        thewx=0.d0
      endif
      if(v(I,J)>0.d0)then
        thewy=(-Zb(I,J+1)-H(I,J+1)+Zb(I,J)+H(I,J))/dy
      elseif(v(I,J)<0.d0)then
        thewy=(-Zb(I,J)-H(I,J)+Zb(I,J-1)+H(I,J-1))/dy
      else
        thewy=0.d0
      endif
      if(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)>0.d0)then
        thewxy=atan((sin(atan(thewx))*u(I,J)+sin(atan(thewy))*v(I,J))/&
        dsqrt(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)))  !nakagawa
      else
        thewxy=0.d0
      endif
      !calc Cainf
      if(abs(tan(thewxy))>=tan(phi))then
        cainf=0.9d0*Cstar
      elseif(abs(tan(thewxy))>0.138d0)then
        cainf=min(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))),0.9d0*Cstar)
      elseif(abs(tan(thewxy))>0.03d0)then
        cainf=6.7d0*(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))))**2
      else
        taustar=rho/(sigma-rho)*h(I,J)*tan(thewxy)/dm
        taustarc=0.04d0*10.d0**(1.72*tan(thewxy))
        alphc=sqrt(2.d0*(0.425d0-sigma*tan(thewxy)/(sigma-rho))/(1.d0-sigma*tan(thewxy)/(sigma-rho)))
        if(taustar>0.d0)then
          cainf=(1.d0+5.d0*tan(thewxy))*tan(thewxy)/(sigma/rho-1.d0)*&
          (1.d0-alphc**2*taustarc/taustar)*(1.d0-alphc*sqrt(taustarc/taustar))
        else
          cainf=0.d0
        endif
      ENDIF
      if(cainf>=Ca(I,J))then
        iero(I,J)=delero*(cainf-Ca(I,J))/(Cstar-cainf)*&
        sqrt(u(I,J)**2+v(I,J)**2)*h(I,J)/dm
      else
        iero(I,J)=deldep*(cainf-Ca(I,J))/Cstar*sqrt(u(I,J)**2+v(I,J)**2)
      endif
      if(caltype==1)then
        Cc(I,J,1)=iero(I,J)
      elseif(caltype==0)then
        Cc(I,J,1)=0.d0
      else
        write(*,*) "error"
      endif
      Cc(I,J,2)=g*h(I,J)*(S0x(I,J)-Sfx(I,J))+uhxx+uhyy
      Cc(I,J,3)=g*h(I,J)*(S0y(I,J)-Sfy(I,J))+vhxx+vhyy
      Cc(I,J,4)=iero(I,J)*Cstar
      Cc(I,J,5)=-iero(I,J)
    ENDDO
    !$omp end parallel do
  enddo
  DO k=1,kmax
    DO J=1,Ny
      !$omp parallel do default(shared),private(I,kx,ky)
      DO I=1,Nx
        ! if(h(I-1,J)>0.d0.and.h(I+1,J)>0.d0)then
        if(I>=1.and.I<=Nx)then
          ! kx=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)/Dx
          kx=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(Sfx(I,J)))*h(I,J)/Dx
          Qxc(I,J,K)=0.125d0*kx*(Uc(I+1,J,K)-2.d0*Uc(I,J,K)+Uc(I-1,J,K))
        else
          Qxc(I,J,K)=0.d0
        endif
        ! if(h(I,J-1)>0.d0.and.h(I,J+1)>0.d0)then
        if(J>1.and.J<Ny)then
          ! ky=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)/Dy
          ky=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(Sfy(I,J)))*h(I,J)/Dy
          Qyc(I,J,K)=0.125d0*ky*(Uc(I,J+1,K)-2.d0*Uc(I,J,K)+Uc(I,J+1,K))
        else
          Qxc(I,J,K)=0.d0
        endif
      ENDDO
      !$omp end parallel do
    ENDDO
  ENDDO

  Up=Uc
  Ep=Ec
  Fp=Fc
  Cp=Cc
  Qxp=Qxc
  Qyp=Qyc


  ! !calc Sf
  ! DO J=0,Ny+1
  !   DO I=0,nx+1
  !     if(h(I,J)>Hcr)then
  !       if(ca(I,J)>=0.4d0*cstar.and.h(I,J)/dm<30.d0)then
  !         Sfx(I,J)=u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
  !         /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
  !         Sfy(I,J)=v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
  !         /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
  !       elseif(ca(I,J)>0.01d0.and.h(I,J)/dm<30.d0)then
  !         Sfx(I,J)=1.d0/0.49d0*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
  !         Sfy(I,J)=1.d0/0.49d0*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
  !       else
  !         Sfx(I,J)=nm**2*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
  !         Sfy(I,J)=nm**2*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
  !       endif
  !     else
  !       Sfx(I,J)=0.d0
  !       Sfy(I,J)=0.d0
  !     endif
  !   ENDDO
  ! ENDDO
  ! !calc Sf0
  ! DO J=1,Ny
  !   DO I=1,Nx
  !     S0x(I,J)=-(Zb(I+1,J)-Zb(I-1,J))/2.d0/Dx
  !     S0y(I,J)=-(Zb(I,J+1)-Zb(I,J-1))/2.d0/Dy
  !   enddo
  ! enddo
  ! S0x(0,1:Ny)=S0x(1,1:Ny)
  ! S0x(Nx+1,1:Ny)=S0x(Nx,1:Ny)
  ! S0x(1:Nx,0)=S0x(1:Nx,1)
  ! S0x(1:Nx,Ny+1)=S0x(1:Nx,Ny)
  ! S0x(0,0)=0.d0
  ! S0x(Nx+1,0)=0.d0
  ! S0x(0,Ny+1)=0.d0
  ! S0x(Nx+1,Ny+1)=0.d0
  ! S0y(0,1:Ny)=S0y(1,1:Ny)
  ! S0y(Nx+1,1:Ny)=S0y(Nx,1:Ny)
  ! S0y(1:Nx,0)=S0y(1:Nx,1)
  ! S0y(1:Nx,Ny+1)=S0y(1:Nx,Ny)
  ! S0y(0,0)=0.d0
  ! S0y(Nx+1,0)=0.d0
  ! S0y(0,Ny+1)=0.d0
  ! S0y(Nx+1,Ny+1)=0.d0
  ! write(*,*) "here ok"
  ! !calc reynolds stress
  ! DO J=1,Ny
  !   DO I=1,Nx
  !     eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)
  !     uhxx=(eps*(u(I+1,J)*h(I+1,J)-u(I,J)*h(I,J))/dx-eps*(u(I,J)*h(I,J)-u(I-1,J)*h(I-1,J))/dx)/dx
  !     uhyy=(eps*(u(I,J+1)*h(I,J+1)-u(I,J)*h(I,J))/dy-eps*(u(I,J)*h(I,J)-u(I,J-1)*h(I,J-1))/dy)/dy
  !     eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)
  !     vhyy=(eps*(v(I,J+1)*h(I,J+1)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I,J-1)*h(I,J-1))/dy)/dy
  !     vhxx=(eps*(v(I+1,J)*h(I+1,J)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I-1,J)*h(I-1,J))/dx)/dx
  !     if(u(I,J)>0.d0)then
  !       thewx=(-Zb(I+1,J)-H(I+1,J)+Zb(I,J)+H(I,J))/dx
  !     elseif(u(I,J)<0.d0)then
  !       thewx=(-Zb(I,J)-H(I,J)+Zb(I-1,J)+H(I-1,J))/dx
  !     else
  !       thewx=0.d0
  !     endif
  !     if(v(I,J)>0.d0)then
  !       thewy=(-Zb(I,J+1)-H(I,J+1)+Zb(I,J)+H(I,J))/dy
  !     elseif(v(I,J)<0.d0)then
  !       thewy=(-Zb(I,J)-H(I,J)+Zb(I,J-1)+H(I,J-1))/dy
  !     else
  !       thewy=0.d0
  !     endif
  !     if(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)>0.d0)then
  !       thewxy=atan((sin(atan(thewx))*u(I,J)+sin(atan(thewy))*v(I,J))/&
  !       dsqrt(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)))  !nakagawa
  !     else
  !       thewxy=0.d0
  !     endif
  !     !calc water surface gradient
  !     if(abs(tan(thewxy))>=tan(phi))then
  !       cainf=0.9d0*Cstar
  !     elseif(abs(tan(thewxy))>0.138d0)then
  !       cainf=min(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))),0.9d0*Cstar)
  !     elseif(abs(tan(thewxy))>0.03d0)then
  !       cainf=6.7d0*(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))))**2
  !     else
  !       taustar=rho/(sigma-rho)*h(I,J)*tan(thewxy)/dm
  !       taustarc=0.04d0*10.d0**(1.72*tan(thewxy))
  !       alphc=sqrt(2.d0*(0.425d0-sigma*tan(thewxy)/(sigma-rho))/(1.d0-sigma*tan(thewxy)/(sigma-rho)))
  !       cainf=(1.d0+5.d0*tan(thewxy))*tan(thewxy)/(sigma/rho-1.d0)*&
  !       (1.d0-alphc**2*taustarc/taustar)*(1.d0-alphc*sqrt(taustarc/taustar))
  !     ENDIF
  !     if(cainf>=Ca(I,J))then
  !       iero(I,J)=delero*(cainf-Ca(I,J))/(Cstar-cainf)*&
  !       sqrt(u(I,J)**2+v(I,J)**2)*h(I,J)/dm
  !     else
  !       iero(I,J)=deldep*(cainf-Ca(I,J))/Cstar*sqrt(u(I,J)**2+v(I,J)**2)
  !     endif
  !     !calc C
  !     if(caltype==1)then
  !       Cc(I,J,1)=iero(I,J)
  !     elseif(caltype==0)then
  !       Cc(I,J,1)=0
  !     else
  !       write(*,*) "error"
  !     endif
  !     Cc(I,J,2)=g*h(I,J)*(S0x(I,J)-Sfx(I,J))+uhxx+uhyy
  !     Cc(I,J,3)=g*h(I,J)*(S0y(I,J)-Sfy(I,J))+vhxx+vhyy
  !     Cc(I,J,4)=iero(I,J)*Cstar
  !     Cc(I,J,5)=-iero(I,J)
  !   ENDDO
  ! ENDDO
  ! Cc(0,0:Ny+1,1:5)=0.d0
  ! Cc(Nx+1,0:Ny+1,1:5)=0.d0
  ! Cc(0:Nx+1,0,1:5)=0.d0
  ! Cc(0:Nx+1,Ny+1,1:5)=0.d0

  ! DO J=0,Ny+1
  !   DO I=0,Nx+1
  !     !calc kxky for Q
  !     DO k=1,kmax
  !     !calc Q
  !       if(I>1.and.I<Nx)then
  !         kx=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)/Dx
  !         Qxc(I,J,K)=0.125d0*kx*(Uc(I+1,J,K)-2.d0*Uc(I,J,K)+Uc(I-1,J,K))
  !       else
  !         Qxc(I,J,K)=0.d0
  !       endif
  !       if(J>1.and.J<Ny)then
  !         ky=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)/Dy
  !         Qyc(I,J,K)=0.125d0*ky*(Uc(I,J+1,K)-2.d0*Uc(I,J,K)+Uc(I,J+1,K))
  !       else
  !         Qyc(I,J,K)=0.d0
  !       endif
  !     ENDDO
  !   ENDDO
  ! ENDDO

end subroutine Initial_condition