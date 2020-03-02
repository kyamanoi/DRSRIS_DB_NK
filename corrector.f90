subroutine corrector(alphain,betain)
  use variables
  implicit none
  real(8),intent(in) ::alphain,betain
  real(8) alpha,beta,kx,ky,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf
  integer :: ktemp
  real(8) taux,tauy,kd,kf,kkf,kkd,ee,fb,phis,phii,phiw,cs,cbar,eta,hh,rr,kds,kfs,ct

  DO k=1,kmax
    DO J=1,Ny
      !$omp parallel do default(shared),private(I,alpha,beta)
      DO I=1,Nx
        if(K==2)then !x方向の計算
          if(Iswet(I,J)==0)then
            Uc(I,J,K)=0.d0
          ! elseif(Iswet(I+1,J)==0.or.Iswet(I-1,J)==0)then!上端
          elseif(Iswet(I+1,J)==0.and.Up(I+1,J,5)>Up(I,J,1)+Up(I,J,5))then!上端
          ! elseif(Iswet(I+1,J)==0)then!上端
            Uc(I,J,K)=0.d0
          elseif(Iswet(I,J+1)==0)then!即岸
            Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
            &+0.5d0*Dt*Cp(I,J,K)
          else
            Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
            &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
          endif
          ! if(Uc(I,J,K)*Up(I,J,K)<0.d0 .and. abs(Uc(I,J,K)/Uc(I,J,1))>100.d0)then
          !   write(*,*) "u is reversed at ", my_rank,I,J, Uc(I,J,K)/Uc(I,J,1),Up(I,J,K)/Up(I,J,1)
          !   Uc(I,J,K)=0.d0
          ! endif
        elseif(k==3)then
          if(Iswet(I,J)==0)then
            Uc(I,J,K)=0.d0
          ! elseif(Iswet(I,J+1)==0.or.Iswet(I,J-1)==0)then!上端
          elseif(Iswet(I,J+1)==0.and.Up(I,J+1,5)>Up(I,J,1)+Up(I,J,5))then!上端
          ! elseif(Iswet(I,J+1)==0)then!上端
            Uc(I,J,K)=0.d0
          elseif(Iswet(I+1,J)==0)then!即岸
            Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))&
            &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
          else
            Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
            &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
          endif
          ! if(Uc(I,J,K)*Up(I,J,K)<0.d0 .and. abs(Uc(I,J,K)/Uc(I,J,1))>100.d0)then
          !   write(*,*) "v is reversed at ", my_rank,I,J, Uc(I,J,K)/Uc(I,J,1),Up(I,J,K)/Up(I,J,1)
          !   Uc(I,J,K)=0.d0
          ! endif
          ! elseif(k==1)then
        !   if(Iswet(I,J+1)==0)then!上端
        !     Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
        !     &+0.5d0*Dt*Cp(I,J,K)
        !   elseif(Iswet(I+1,J)==0)then!即岸
        !     Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))&
        !     &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
        !   else
        !     Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
        !     &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
        !   endif
        else
          Uc(I,J,K)=0.5d0*(Uc(I,J,K)+Up(I,J,K))-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))&
          &-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))+0.5d0*Dt*Cp(I,J,K)
        endif
        ! if(K==5.and.Uc(I,J,K)<Zbini(I,J))then
        !   write(*,*) "Zb is smaller than zbini (corrector)"
        !   write(*,*) T,Dt
        !   write(*,*) I,J,Zbini(I,J),Uc(I,J,K)
        !   write(*,*) "zb=",Up(0,J,K),Up(1,J,K),Up(2,J,K),Up(3,J,K),Up(4,J,K),Up(5,J,K)
        !   write(*,*) Qxp(I+1,J,K),Qxp(I,J,K)

        ! if(K==2.and.Uc(I,J,2)/Uc(I,J,1)<-1000.d0)then
        !   write(*,*) "too small uc",Uc(I,J,K)
        !   write(*,*) I,J,K,T,Dt
        !   write(*,*) "Up=",Up(I,J,K)
        !   write(*,*) "uc1",Uc(I,J,2)/Uc(I,J,1)
        !   write(*,*) "uc2",Up(I,J,2)/Up(I,J,1)
        !   write(*,*) "term1",-0.5d0*Dt/Dx*((Ep(I+1,J,K)-Ep(I,J,K))+(Qxp(I+1,J,K)-Qxp(I,J,K)))
        !   write(*,*) "term2",-0.5d0*Dt/Dy*((Fp(I,J+1,K)-Fp(I,J,K))+(Qyp(I,J+1,K)-Qyp(I,J,K)))
        !   write(*,*) Cp(I-1,J-1,K),Cp(I,J-1,K),Cp(I+1,J-1,K)
        !   write(*,*) Cp(I-1,J,K),Cp(I,J,K),Cp(I+1,J,K)
        !   write(*,*) Cp(I-1,J+1,K),Cp(I,J+1,K),Cp(I+1,J+1,K)
        !   read(*,*)
        ! endif
        !   read(*,*)
        ! endif
        ! if(K==4.and.Uc(I,J,4)/Uc(I,J,1)<0.d0)then
        !   write(*,*) "ca is smaller than zero",Uc(I,J,K)
        !   write(*,*) I,J,K,T,Dt
        !   read(*,*)
        ! endif

        if(K==1.and.Uc(I,J,K)<0.d0)then
          ! write(*,*) Uc(I,J,K)
          ! write(*,*) "uc1 is smaller than zero in corrector"
          ! write(*,*) Iswet(I-1,J+1),Iswet(I,J+1),Iswet(I+1,J+1)
          ! write(*,*) Iswet(I-1,J),Iswet(I,J),Iswet(I+1,J)
          ! write(*,*) Iswet(I-1,J-1),Iswet(I,J-1),Iswet(I+1,J-1)
          ! write(*,*) "up1(h)"
          ! write(*,*) up(I-1,J+1,1),up(I,J+1,1),up(I+1,J+1,1)
          ! write(*,*) up(I-1,J,1),up(I,J,1),up(I+1,J,1)
          ! write(*,*) up(I-1,J-1,1),up(I,J-1,1),up(I+1,J-1,1)
          ! write(*,*) "up2(qx)"
          ! write(*,*) up(I-1,J+1,2),up(I,J+1,2),up(I+1,J+1,2)
          ! write(*,*) up(I-1,J,2),up(I,J,2),up(I+1,J,2)
          ! write(*,*) up(I-1,J-1,2),up(I,J-1,2),up(I+1,J-1,2)
          ! write(*,*) "up3(qy)"
          ! write(*,*) up(I-1,J+1,3),up(I,J+1,3),up(I+1,J+1,3)
          ! write(*,*) up(I-1,J,3),up(I,J,3),up(I+1,J,3)
          ! write(*,*) up(I-1,J-1,3),up(I,J-1,3),up(I+1,J-1,3)
          ! write(*,*) "Zb"
          ! write(*,*) zb(I-1,J+1),zb(I,J+1),zb(I+1,J+1)
          ! write(*,*) zb(I-1,J),zb(I,J),zb(I+1,J)
          ! write(*,*) zb(I-1,J-1),zb(I,J-1),zb(I+1,J-1)

          ! write(*,*) Uc(I,J,1), Uc(I,J,2), Uc(I,J,3)
          ! write(*,*) Uc(I+1,J,1), Uc(I+1,J,2), Uc(I+1,J,3)
          ! write(*,*) Uc(I,J+1,1), Uc(I,J+1,2), Uc(I,J+1,3)
          ! write(*,*) I,J,K,T,Dt
          ! call fwrite_xyz("h",h)
          ! call fwrite_xyz("u",u)
          ! call fwrite_xyz("v",v)
          ! call fwrite_xyz("zb",zb)
          ! call fwrite_xyz("ca",ca)
          ! read(*,*)
          Uc(I,J,K)=0.d0
        endif
        if(k==4)then
          if(Uc(I,J,K)<0.d0)then
            Uc(I,J,K)=0.d0
          elseif(Uc(I,J,1)>0.d0.and.Uc(I,J,K)/Uc(I,J,1)>cstar*0.9d0)then
            Uc(I,J,k)=cstar*0.9d0*Uc(I,J,1)
          endif
        endif
      ENDDO
      !$omp end parallel do
    ENDDO
  ENDDO

  DO J=0,Ny+1
    !$omp parallel do default(shared),private(I)
    DO I=0,Nx+1
      if(Uc(I,J,1)>Hcr)then
        h(I,J)=Uc(I,J,1)
        u(I,J)=Uc(I,J,2)/Uc(I,J,1)
        v(I,J)=Uc(I,J,3)/Uc(I,J,1)
        Ca(I,J)=Uc(I,J,4)/Uc(I,J,1)
        if(Ca(I,J)>cstar)then
          write(*,*)"ca is larger than cster (corrector) at ",I,J,Ca(I,J),cstar
        endif
      else
        h(I,J)=Uc(I,J,1)
        u(I,J)=0.d0
        v(I,J)=0.d0
        ca(I,J)=0.d0
      endif
      zb(I,J)=Uc(I,J,5)
    ENDDO
    !$omp end parallel do
  ENDDO

  if(PETOT>1)then
    call mpi_boundary_2d_y(my_rank,h(0:Nx+1,0:Ny+1))
    call mpi_boundary_2d_x(my_rank,h(0:Nx+1,0:Ny+1))
  endif
  call dry_or_wet(h,zb,u,v,iswet,hcr,Nx,Ny,1)
  ! call dry_or_wet_f(h,zb,u,v,iswet,hcr,Nx,Ny)
  Uc(0:Nx+1,0:Ny+1,1)=h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,2)=u(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,3)=v(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,4)=Ca(0:Nx+1,0:Ny+1)*h(0:Nx+1,0:Ny+1)
  Uc(0:Nx+1,0:Ny+1,5)=Zb(0:Nx+1,0:Ny+1)

  ! where(h>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  DO J=1,Ny
    !$omp parallel do default(shared),&
    !$omp private(I,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf,&
    !$omp taux,tauy,kd,kf,kkf,kkd,ee,fb,phis,phii,phiw,cs,cbar,eta,hh,rr,kds,kfs,ct)
    DO I=1,Nx
      if(iswet(I,J)==1.and.iswet(I+1,J)==1)then
        thewx=(-Zb(I+1,J)-H(I+1,J)+Zb(I,J)+H(I,J))/dx
      else
        thewx=0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I,J+1)==1)then
      thewy=(-Zb(I,J+1)-H(I,J+1)+Zb(I,J)+H(I,J))/dy
      else
        thewy=0.d0
      endif
      if(abs(u(I,J)**2+v(I,J)**2)>0.d0)then
        thewxy=atan((thewx*u(I,J)+thewy*v(I,J))/sqrt(u(I,J)**2+v(I,J)**2))
      else
        thewxy=0.d0
      endif
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
        !takahashi
        if(ca(I,J)>=0.4d0*cstar.and.h(I,J)/dm<30.d0)then
        ! if(abs(tan(thewxy))>0.138d0)then
          Sfx(I,J)=u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
          /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
          Sfy(I,J)=v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2&
          /(8.d0*h(I,J)**3.d0*(ca(I,J)+(1.d0-ca(I,J))*rho/sigma)*((Cstar/ca(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
        elseif(ca(I,J)>0.01d0.and.h(I,J)/dm<30.d0)then
        ! elseif(abs(tan(thewxy))>0.03d0)then
          Sfx(I,J)=1.d0/0.49d0*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
          Sfy(I,J)=1.d0/0.49d0*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/g
        else
          Sfx(I,J)=nm**2*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
          Sfy(I,J)=nm**2*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
        endif
        !egashira
        ! kd=0.0828d0
        ! kf=0.2d0
        ! ee=0.85d0
        ! cs=0.5d0*Cstar
        ! cbar=ca(I,J)
        ! hh=h(I,J)
        ! if(cbar>cs)then
        !   kkd=kd*sigma/rho*(1-ee**2)*cbar**(1.d0/3.d0)
        !   kkf=kf*(1.d0-cbar)**(5.d0/3.d0)*cbar**(2.d0/3.d0)
        !   fb=25.d0/4.d0*(kkd+kkf)*(hh/dm)**(-2.d0)
        ! else
        !   eta=sqrt(kf)*((1.d0-cs)/cs)**(1.d0/3.d0)*dm
        !   phiw=1.d0/karman*sqrt(1.d0-cbar/cs)*(eta/hh+1.d0-cbar/cs)*(hh/dm)**(-1)*log((eta/hh+1.d0-cbar/cs)/(eta/hh)) &
        !   & - 1/karman*(1.d0-cbar/cs)**(1.5d0)*(hh/dm)**(-1)
        !   rr=1.d0+cs/cstar*((1.d0-(cs/cstar)**0.2d0)*((sigma/rho-1)*cbar+1.d0)-1.d0)
        !   kds=kd*sigma/rho*(1.d0-ee**2)*cs**(1.d0/3.d0)
        !   kfs=kf*(1.d0-cs)**(5.d0/3.d0)*cs**(2.d0/3.d0)
        !   if(rr==0.d0)then
        !     phis=1.d0/2.d0/sqrt(kds+kfs)*(cbar/cs)**2*sqrt(1-cbar/cs)
        !     phii=1.d0/sqrt(kds+kfs)*cbar/cs*(1.d0-cbar/cs)**1.5d0
        !   else
        !     phis=4.d0/15.d0/(rr**2*sqrt(kds+kfs)) &
        !     & *((1.d0-cbar/cs)**2.5d0-((1.d0-cbar/cs)+rr*cbar/cs)**1.5d0*((1.d0-cbar/cs)-1.5d0*rr*cbar/cs))
        !     phii=-2.d0/3.d0*(1.d0-cbar/cs)/(rr*sqrt(kds+kfs))*((1-cbar/cs)**2.5d0-((1.d0-cbar/cs)+rr*cbar/cs)**1.5d0)
        !   endif
        !   fb=(1.d0-(cs/cstar)**0.2d0)*((sigma/rho-1.d0)*Cbar+1.d0)*(phis+phii+phiw)**(-2)*(hh/dm)**(-2)
        ! endif
        ! taux=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zb(I+1,J)-Zb(I,J))/Dx))*tan(phi)
        ! tauy=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zb(I,J+1)-Zb(I,J))/Dy))*tan(phi)
        ! sfx(I,J)=(taux+rho*fb*u(I,J)*sqrt(u(I,J)**2+v(I,J)**2))/rho/hh
        ! sfy(I,J)=(tauy+rho*fb*v(I,J)*sqrt(u(I,J)**2+v(I,J)**2))/rho/hh
      !egashiratakahashimix
        ! kd=0.0828d0
        ! kf=0.2d0
        ! ee=0.85d0
        ! cs=0.5d0*Cstar
        ! cbar=ca(I,J)
        ! hh=h(I,J)
        ! if(ca(I,J)>=0.4d0*cstar.and.h(I,J)/dm<30.d0)then
        !   kkd=kd*sigma/rho*(1-ee**2)*cbar**(1.d0/3.d0)
        !   kkf=kf*(1.d0-cbar)**(5.d0/3.d0)*cbar**(2.d0/3.d0)
        !   fb=25.d0/4.d0*(kkd+kkf)*(hh/dm)**(-2.d0)
        !   taux=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zb(I+1,J)-Zb(I,J))/Dx))*tan(phi)
        !   tauy=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zb(I,J+1)-Zb(I,J))/Dy))*tan(phi)
        !   sfx(I,J)=(taux+rho*fb*u(I,J)*sqrt(u(I,J)**2+v(I,J)**2))/rho/g
        !   sfy(I,J)=(tauy+rho*fb*v(I,J)*sqrt(u(I,J)**2+v(I,J)**2))/rho/g
        ! elseif(ca(I,J)>0.01d0.and.h(I,J)/dm<30.d0)then
        !   Sfx(I,J)=1.d0/0.49d0*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/hh
        !   Sfy(I,J)=1.d0/0.49d0*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)*dm**2/h(I,J)**3.d0/hh
        ! else
        !   Sfx(I,J)=nm**2*u(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
        !   Sfy(I,J)=nm**2*v(I,J)*sqrt(v(I,J)**2+u(I,J)**2)/h(I,J)**(4.d0/3.d0)
        ! endif
      else
        Sfx(I,J)=0.d0
        Sfy(I,J)=0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I+1,J)==1)then
        S0x(I,J)=-(Zb(I+1,J)-Zb(I,J))/Dx
      else
        S0x(I,J)=-(Zb(I+1,J)-Zb(I,J))/Dx !0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I,J+1)==1)then
        S0y(I,J)=-(Zb(I,J+1)-Zb(I,J))/Dy !collectorなので前進差分
      else
        S0y(I,J)=-(Zb(I,J+1)-Zb(I,J))/Dy  !0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)
      if(iswet(I+1,J)==1.and.iswet(I,J)==1.and.iswet(I-1,J)==1)then
        uhxx=(eps*(u(I+1,J)*h(I+1,J)-u(I,J)*h(I,J))/dx-eps*(u(I,J)*h(I,J)-u(I-1,J)*h(I-1,J))/dx)/dx
      else
        uhxx=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J)==1.and.iswet(I,J-1)==1)then
        uhyy=(eps*(u(I,J+1)*h(I,J+1)-u(I,J)*h(I,J))/dy-eps*(u(I,J)*h(I,J)-u(I,J-1)*h(I,J-1))/dy)/dy
      else
        uhyy=0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)
      if(iswet(I+1,J)==1.and.iswet(I,J)==1.and.iswet(I-1,J)==1)then
        vhyy=(eps*(v(I,J+1)*h(I,J+1)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I,J-1)*h(I,J-1))/dy)/dy
      else
        vhyy=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J)==1.and.iswet(I,J-1)==1)then
        vhxx=(eps*(v(I+1,J)*h(I+1,J)-v(I,J)*h(I,J))/dy-eps*(v(I,J)*h(I,J)-v(I-1,J)*h(I-1,J))/dx)/dx
      else
        vhxx=0.d0
      endif
      !calc water surface gradient!要検討
      ! if(u(I,J)>0.d0)then
      !   thewx=(-Zb(I+1,J)-H(I+1,J)+Zb(I,J)+H(I,J))/dx
      ! elseif(u(I,J)<0.d0)then
      !   thewx=(-Zb(I,J)-H(I,J)+Zb(I-1,J)+H(I-1,J))/dx
      ! else
      !   thewx=0.d0
      ! endif
      ! if(v(I,J)>0.d0)then
      !   thewy=(-Zb(I,J+1)-H(I,J+1)+Zb(I,J)+H(I,J))/dy
      ! elseif(v(I,J)<0.d0)then
      !   thewy=(-Zb(I,J)-H(I,J)+Zb(I,J-1)+H(I,J-1))/dy
      ! else
      !   thewy=0.d0
      ! endif
      ! if(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)>0.d0)then
      !   thewxy=atan((sin(atan(thewx))*u(I,J)+sin(atan(thewy))*v(I,J))/&
      !   dsqrt(u(I,J)**2*cos(atan(thewx))**2+v(I,J)**2*cos(atan(thewy)**2)))  !nakagawa
      ! else
      !   thewxy=0.d0
      ! endif
      !calc Cainf
      if(abs(tan(thewxy))>=tan(phi))then
        cainf=0.9d0*Cstar
      elseif(abs(tan(thewxy))>0.138d0)then
      ! elseif(ca(I,J)>=0.4d0*cstar.and.h(I,J)/dm<30.d0)then
        cainf=min(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))),0.9d0*Cstar)
      elseif(abs(tan(thewxy))>0.03d0)then
      ! elseif(ca(I,J)>0.01d0.and.h(I,J)/dm<30.d0)then
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
        if(Zb(I,J)>Zbmin(I,J))then
          iero(I,J)=min(delero*(cainf-Ca(I,J))/(Cstar-cainf)*&
          sqrt(u(I,J)**2+v(I,J)**2)*h(I,J)/dm,(zb(I,J)-Zbmin(I,J))/Dt)
          ! iero(I,J)=delero*(cainf-Ca(I,J))/(Cstar-cainf)*&
          ! sqrt(u(I,J)**2+v(I,J)**2)*h(I,J)/dm
        else
          iero(I,J)=0.d0
        endif
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
      if((Uc(I,J,2)+Dt*(g*h(I,J)*(S0x(I,J)-Sfx(I,J))+uhxx+uhyy))*Uc(I,J,2)<0.d0)then !avoid reverse flow
        Sfx(I,J)=S0x(I,J)-(-2.d0*Uc(I,J,2)/Dt-uhxx-uhyy)/(g*h(I,J))
      endif
      if((Uc(I,J,3)+Dt*(g*h(I,J)*(S0y(I,J)-Sfy(I,J))+vhxx+vhyy))*Uc(I,J,3)<0.d0)then !avoid reverse flow
        Sfy(I,J)=S0y(I,J)-(-2.d0*Uc(I,J,3)/Dt-vhxx-vhyy)/(g*h(I,J))
      endif
      Cc(I,J,2)=g*h(I,J)*(S0x(I,J)-Sfx(I,J))+uhxx+uhyy
      Cc(I,J,3)=g*h(I,J)*(S0y(I,J)-Sfy(I,J))+vhxx+vhyy
      Cc(I,J,4)=iero(I,J)*Cstar
      Cc(I,J,5)=-iero(I,J)
      ! if(abs(iero(I,J))>0.d0)then
      !   write(*,*) "iero has an value"
      !   write(*,*) I,J,Ca(I,J),iero(I,J)
      !   read(*,*)
      ! endif
    ENDDO
    !$omp end parallel do
  enddo
  DO k=1,kmax
    DO J=1,Ny
      !$omp parallel do default(shared),private(I,kx,ky)
      DO I=1,Nx
        ! if(h(I-1,J)>0.d0.and.h(I+1,J)>0.d0)then
        ! if(I>=1.and.I<=Nx)then
          ! kx=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0x(I,J)+Sfx(I,J)))*h(I,J)/Dx
          kx=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(Sfx(I,J)))*h(I,J)/Dx
          ! if(K==5.and.Zb(I,J)<Zbmin(I,J))then
          !   kx=0.d0
          ! endif
          Qxc(I,J,K)=0.125d0*kx*(Uc(I+1,J,K)-2.d0*Uc(I,J,K)+Uc(I-1,J,K))
        ! else
          ! Qxc(I,J,K)=0.d0
        ! endif
        ! if(h(I,J-1)>0.d0.and.h(I,J+1)>0.d0)then
        ! if(J>1.and.J<Ny)then
          ! ky=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(S0y(I,J)+Sfy(I,J)))*h(I,J)/Dy
          ky=8.d0*Kv(K)*sqrt(g*h(I,J)*abs(Sfy(I,J)))*h(I,J)/Dy
          ! if(K==5.and.Zb(I,J)<Zbmin(I,J))then
          !   ky=0.d0
          ! endif
          Qyc(I,J,K)=0.125d0*ky*(Uc(I,J+1,K)-2.d0*Uc(I,J,K)+Uc(I,J+1,K))
        ! else
          ! Qxc(I,J,K)=0.d0
        ! endif
        ! if(k==5.and.Qxc(I,J,K)/=0.d0)then
        !   write(*,*) "bbbbb"
        !   write(*,*) I,J,K
        !   read(*,*)
        ! endif
      ENDDO
      !$omp end parallel do
    ENDDO
  ENDDO


end subroutine corrector