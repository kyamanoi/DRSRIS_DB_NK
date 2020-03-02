subroutine predictor(alphain,betain)
  use variables
  implicit none
  real(8),intent(in) ::alphain,betain
  real(8) alpha,beta,kx,ky,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf
  ! real(8),dimension(0:nx+1,0:ny+1) :: htemp,utemp,vtemp,catemp,zbtemp,sfxtemp,sfytemp,s0xtemp,s0ytemp,ierotemp
  real(8) taux,tauy,kd,kf,kkf,kkd,ee,fb,phis,phii,phiw,cs,cbar,eta,hh,rr,kds,kfs,ct

  DO k=1,kmax
    DO J=1,Ny
      !$omp parallel do default(shared),private(I,alpha,beta)
      DO I=1,Nx
        if(K==2)then !x方向の計算
          if(Iswet(I,J)==0)then
            Up(I,J,K)=0.d0
          ! elseif(Iswet(I-1,J)==0.or.Iswet(I+1,J)==0)then !上端
          elseif(Iswet(I-1,J)==0.and.Uc(I-1,J,5)>Uc(I,J,1)+Uc(I,J,5))then !上端
          ! elseif(Iswet(I-1,J)==0)then !上端
            Up(I,J,K)=0.d0
          elseif(Iswet(I,J-1)==0)then !即岸
            Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
            &+Dt*Cc(I,J,K)
          else
            Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
            &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
          endif
          ! if(Up(I,J,K)*Uc(I,J,K)<0.d0 .and. abs(Up(I,J,K)/Up(I,J,1))>100.d0)then
          !   write(*,*) "u is reversed at ", my_rank,I,J, Up(I,J,K)/Up(I,J,1),Uc(I,J,K)/Uc(I,J,1)
          !   Up(I,J,K)=0.d0
          ! endif
        elseif(k==3)then
          if(Iswet(I,J)==0)then
            Up(I,J,K)=0.d0
          ! elseif(Iswet(I,J-1)==0.or.Iswet(I,J+1)==0)then !上端
          elseif(Iswet(I,J-1)==0.and.Uc(I,J-1,5)>Uc(I,J,1)+Uc(I,J,5))then !上端
          ! elseif(Iswet(I,J-1)==0)then !上端
            Up(I,J,K)=0.d0
          elseif(Iswet(I-1,J)==0)then !即岸
            Up(I,J,K)=Uc(I,J,K)&
            &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
          else
            Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
            &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
          endif
          ! if(Up(I,J,K)*Uc(I,J,K)<0.d0 .and. abs(Up(I,J,K)/Up(I,J,1))>100.d0)then
          !   write(*,*) "v is reversed at ", my_rank,I,J, Up(I,J,K)/Up(I,J,1),Uc(I,J,K)/Uc(I,J,1)
          !   Up(I,J,K)=0.d0
          ! endif
          ! elseif(k==1)then
        !   if(Iswet(I,J-1)==0)then !上端
        !     Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
        !     &+Dt*Cc(I,J,K)
        !   elseif(Iswet(I-1,J)==0)then !即岸
        !     Up(I,J,K)=Uc(I,J,K)&
        !     &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
        !   else
        !     Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
        !     &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
        !   endif
        else
          Up(I,J,K)=Uc(I,J,K)-Dt/Dx*((Ec(I,J,K)-Ec(I-1,J,K))-(Qxc(I,J,k)-Qxc(I-1,J,K)))&
          &-Dt/Dy*((Fc(I,J,K)-Fc(I,J-1,K))-(Qyc(I,J,K)-Qyc(I,J-1,K)))+Dt*Cc(I,J,K)
        endif
        ! if(K==5.and.Up(I,J,K)<Zbini(I,J))then
        !   write(*,*) "Zb is smaller than zbini (predictor)"
        !   write(*,*) T,Dt
        !   write(*,*) I,J,Zbini(I,J),Up(I,J,K)
        !   write(*,*) Qxc(I-1,J,K),Qxc(I,J,k),Qxc(I+1,J,K)
          
        !   read(*,*)
        ! endif

        if(K==1.and.Up(I,J,K)<0.d0)then
          ! write(*,*) "up1 is smaller than zero in predictor"
          ! write(*,*) Up(I,J,K)
          ! write(*,*) I,J,K,T,Dt
          ! write(*,*) "iswqet(h)"
          ! write(*,*) Iswet(I-1,J+1),Iswet(I,J+1),Iswet(I+1,J+1)
          ! write(*,*) Iswet(I-1,J),Iswet(I,J),Iswet(I+1,J)
          ! write(*,*) Iswet(I-1,J-1),Iswet(I,J-1),Iswet(I+1,J-1)
          ! write(*,*) "uc1(h)"
          ! write(*,*) uc(I-1,J+1,1),uc(I,J+1,1),uc(I+1,J+1,1)
          ! write(*,*) uc(I-1,J,1),uc(I,J,1),uc(I+1,J,1)
          ! write(*,*) uc(I-1,J-1,1),uc(I,J-1,1),uc(I+1,J-1,1)
          ! write(*,*) "uc2(qx)"
          ! write(*,*) uc(I-1,J+1,2),uc(I,J+1,2),uc(I+1,J+1,2)
          ! write(*,*) uc(I-1,J,2),uc(I,J,2),uc(I+1,J,2)
          ! write(*,*) uc(I-1,J-1,2),uc(I,J-1,2),uc(I+1,J-1,2)
          ! write(*,*) "uc3(qy)"
          ! write(*,*) uc(I-1,J+1,3),uc(I,J+1,3),uc(I+1,J+1,3)
          ! write(*,*) uc(I-1,J,3),uc(I,J,3),uc(I+1,J,3)
          ! write(*,*) uc(I-1,J-1,3),uc(I,J-1,3),uc(I+1,J-1,3)
          ! write(*,*) "Zbtemp"
          ! write(*,*) zbtemp(I-1,J+1),zbtemp(I,J+1),zbtemp(I+1,J+1)
          ! write(*,*) zbtemp(I-1,J),zbtemp(I,J),zbtemp(I+1,J)
          ! write(*,*) zbtemp(I-1,J-1),zbtemp(I,J-1),zbtemp(I+1,J-1)

          ! call fwrite_xyz("h",h)
          ! call fwrite_xyz("u",u)
          ! call fwrite_xyz("v",v)
          ! call fwrite_xyz("zb",zb)
          ! call fwrite_xyz("ca",ca)
          ! read(*,*)
          Up(I,J,K)=0.d0
        endif
        if(k==4)then
          if(Up(I,J,K)<0.d0)then
            Up(I,J,K)=0.d0
          elseif(Up(I,J,1)>0.d0.and.Up(I,J,K)/Up(I,J,1)>cstar*0.9d0)then
            Up(I,J,k)=cstar*0.9d0*Up(I,J,1)
          endif
        endif

        ! if(K==4.and.Up(I,J,4)/Up(I,J,1)<0.d0)then
        !   write(*,*) "ca is smaller than zero (pred)",Up(I,J,K)
        !   write(*,*) I,J,K,T,Dt
        !   read(*,*)
        ! endif

      ENDDO
      !$omp end parallel do
    ENDDO
  ENDDO

  ! write(*,*) iswet(100,10),iswet(101,10),iswet(102,10)
  ! write(*,*) iswet(100,11),iswet(101,11),iswet(102,11)
  ! write(*,*) iswet(100,12),iswet(101,12),iswet(102,12)

  ! write(*,*) "h_c(101,11)=",Uc(101,11,1)
  ! write(*,*) "qx_c(101,11)=",Uc(101,11,2)
  ! write(*,*) "qy_c(101,11)=",Uc(101,11,3)

  ! write(*,*) "cc1(101,11)=",Cc(101,11,1)
  ! write(*,*) "cc2(101,11)=",Cc(101,11,2)
  ! write(*,*) "cc3(101,11)=",Cc(101,11,3)

  ! write(*,*) "h_pred(101,11)=",Up(101,11,1)
  ! write(*,*) "qx_pred(101,11)=",Up(101,11,2)
  ! write(*,*) "qy_pred(101,11)=",Up(101,11,3)

  DO J=0,Ny+1
    !$omp parallel do default(shared),private(I)!,htemp,utemp,vtemp,sfxtemp,sfytemp,kx,ky,eps,uhxx,uhyy,vhxx,vhyy)
    DO I=0,Nx+1
      if(Up(I,J,1)>Hcr)then
        htemp(I,J)=Up(I,J,1)
        utemp(I,J)=Up(I,J,2)/Up(I,J,1)
        vtemp(I,J)=Up(I,J,3)/Up(I,J,1)
        Catemp(I,J)=Up(I,J,4)/Up(I,J,1)
        if(Catemp(I,J)>cstar)then
          write(*,*)"catemp is larger than cster (predictor) at ",I,J,Catemp(I,J),cstar
        endif
      else
        htemp(I,J)=Up(I,J,1)
        utemp(I,J)=0.d0
        vtemp(I,J)=0.d0
        Catemp(I,J)=0.d0
      endif
      if(caltype==1)then
        zbtemp(I,J)=Up(I,J,5)
      else
        zbtemp(I,J)=Uc(I,J,5)
      endif
    enddo
    !$omp end parallel do
  enddo
  if(PETOT>1)then
    call mpi_boundary_2d_y(my_rank,htemp(0:Nx+1,0:Ny+1))
    call mpi_boundary_2d_x(my_rank,htemp(0:Nx+1,0:Ny+1))
  endif
  call dry_or_wet(htemp,zbtemp,utemp,vtemp,iswet,hcr,Nx,Ny,0)

  ! call dry_or_wet_b(htemp,zbtemp,utemp,vtemp,iswet,hcr,Nx,Ny)
  Up(0:Nx+1,0:Ny+1,1)=htemp(0:Nx+1,0:Ny+1)
  Up(0:Nx+1,0:Ny+1,2)=utemp(0:Nx+1,0:Ny+1)*htemp(0:Nx+1,0:Ny+1)
  Up(0:Nx+1,0:Ny+1,3)=vtemp(0:Nx+1,0:Ny+1)*htemp(0:Nx+1,0:Ny+1)
  Up(0:Nx+1,0:Ny+1,4)=Catemp(0:Nx+1,0:Ny+1)*htemp(0:Nx+1,0:Ny+1)
  Up(0:Nx+1,0:Ny+1,5)=Zbtemp(0:Nx+1,0:Ny+1)
  
  ! where(htemp>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  DO J=1,Ny
    !$omp parallel do default(shared),&
    !$omp private(I,kx,ky,eps,uhxx,uhyy,vhxx,vhyy,thewx,thewy,thewxy,taustar,taustarc,alphc,cainf,&
    !$omp taux,tauy,kd,kf,kkf,kkd,ee,fb,phis,phii,phiw,cs,cbar,eta,hh,rr,kds,kfs,ct)
    DO I=1,Nx
      if(iswet(I,J)==1.and.iswet(I-1,J)==1)then
        thewx=(-Zbtemp(I,J)-Htemp(I,J)+Zbtemp(I-1,J)+Htemp(I-1,J))/dx
      else
        thewx=0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I,J-1)==1)then
        thewy=(-Zbtemp(I,J)-Htemp(I,J)+Zbtemp(I,J-1)+Htemp(I,J-1))/dy
      else
        thewy=0.d0
      endif
      if(abs(utemp(I,J)**2+vtemp(I,J)**2)>0.d0)then
        thewxy=atan(abs((thewx*utemp(I,J)+thewy*vtemp(I,J))/sqrt(utemp(I,J)**2+vtemp(I,J)**2)))
      else
        thewxy=0.d0
      endif
      !calc E
      Ep(I,J,1)=utemp(I,J)*htemp(I,J)
      Ep(I,J,2)=utemp(I,J)**2*htemp(I,J)+0.5d0*g*htemp(I,J)**2
      Ep(I,J,3)=utemp(I,J)*vtemp(I,J)*htemp(I,J)
      Ep(I,J,4)=catemp(I,J)*utemp(I,J)*htemp(I,J)
      Ep(I,J,5)=0.d0
      !calc F
      Fp(I,J,1)=vtemp(I,J)*htemp(I,J)
      Fp(I,J,2)=utemp(I,J)*vtemp(I,J)*htemp(I,J)
      Fp(I,J,3)=vtemp(I,J)**2*htemp(I,J)+0.5d0*g*htemp(I,J)**2
      Fp(I,J,4)=catemp(I,J)*vtemp(I,J)*htemp(I,J)
      Fp(I,J,5)=0.d0
      !calc Sfxtemp
      if(iswet(I,J)==1)then
        ! takahashi
        if(catemp(I,J)>=0.4d0*cstar.and.htemp(I,J)/dm<30.d0)then
        ! if(abs(tan(thewxy))>0.138d0)then
          Sfxtemp(I,J)=utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2&
          /(8.d0*htemp(I,J)**3.d0*(catemp(I,J)+(1.d0-catemp(I,J))*rho/sigma)*((Cstar/catemp(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
          Sfytemp(I,J)=vtemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2&
          /(8.d0*htemp(I,J)**3.d0*(catemp(I,J)+(1.d0-catemp(I,J))*rho/sigma)*((Cstar/catemp(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
        elseif(catemp(I,J)>0.01d0.and.htemp(I,J)/dm<30.d0)then
        ! elseif(abs(tan(thewxy))>0.03d0)then
          Sfxtemp(I,J)=1.d0/0.49d0*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2/htemp(I,J)**3.d0/g
          Sfytemp(I,J)=1.d0/0.49d0*vtemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2/htemp(I,J)**3.d0/g
        else
          Sfxtemp(I,J)=nm**2*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)/htemp(I,J)**(4.d0/3.d0)
          Sfytemp(I,J)=nm**2*vtemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)/htemp(I,J)**(4.d0/3.d0)
        endif
        !egashira
        ! kd=0.0828d0
        ! kf=0.2d0
        ! ee=0.85d0
        ! cs=0.5d0*Cstar
        ! cbar=catemp(I,J)
        ! hh=htemp(I,J)
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
        ! taux=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zbtemp(I+1,J)-Zbtemp(I,J))/Dx))*tan(phi)
        ! tauy=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zbtemp(I,J+1)-Zbtemp(I,J))/Dy))*tan(phi)
        ! sfxtemp(I,J)=(taux+rho*fb*utemp(I,J)*sqrt(utemp(I,J)**2+vtemp(I,J)**2))/rho/hh
        ! sfytemp(I,J)=(tauy+rho*fb*vtemp(I,J)*sqrt(utemp(I,J)**2+vtemp(I,J)**2))/rho/hh
        !egashiratakahashimix
        ! kd=0.0828d0
        ! kf=0.2d0
        ! ee=0.85d0
        ! cs=0.5d0*Cstar
        ! cbar=catemp(I,J)
        ! hh=htemp(I,J)
        ! if(catemp(I,J)>=0.4d0*cstar.and.htemp(I,J)/dm<30.d0)then
        !   kkd=kd*sigma/rho*(1-ee**2)*cbar**(1.d0/3.d0)
        !   kkf=kf*(1.d0-cbar)**(5.d0/3.d0)*cbar**(2.d0/3.d0)
        !   fb=25.d0/4.d0*(kkd+kkf)*(hh/dm)**(-2.d0)
        !   taux=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zbtemp(I+1,J)-Zbtemp(I,J))/Dx))*tan(phi)
        !   tauy=(cbar/cstar)**0.2d0*(sigma-rho)*cbar*g*hh*cos(atan(-(Zbtemp(I,J+1)-Zbtemp(I,J))/Dy))*tan(phi)
        !   sfxtemp(I,J)=(taux+rho*fb*utemp(I,J)*sqrt(utemp(I,J)**2+vtemp(I,J)**2))/rho/g
        !   sfytemp(I,J)=(tauy+rho*fb*vtemp(I,J)*sqrt(utemp(I,J)**2+vtemp(I,J)**2))/rho/g
        ! elseif(catemp(I,J)>0.01d0.and.htemp(I,J)/dm<30.d0)then
        !   Sfxtemp(I,J)=1.d0/0.49d0*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2/htemp(I,J)**3.d0/hh
        !   Sfytemp(I,J)=1.d0/0.49d0*vtemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2/htemp(I,J)**3.d0/hh
        ! else
        !   Sfxtemp(I,J)=nm**2*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)/htemp(I,J)**(4.d0/3.d0)
        !   Sfytemp(I,J)=nm**2*vtemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)/htemp(I,J)**(4.d0/3.d0)
        ! endif
      else
        Sfxtemp(I,J)=0.d0
        Sfytemp(I,J)=0.d0
      endif
      ! if(sfxtemp(I,J)>10000)then
      !   write(*,*) "too large sfxtemp",sfxtemp(I,J),T
      !   write(*,*) "debrisflow",utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2&
      !   /(8.d0*htemp(I,J)**3.d0*(catemp(I,J)+(1.d0-catemp(I,J))*rho/sigma)*((Cstar/catemp(I,J))**(1.d0/3.d0)-1.d0)**2.d0)/g
      !   write(*,*) "debrisflow",utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2&
      !   /(8.d0*htemp(I,J)**3.d0*(0.54d0+(1.d0-0.54d0)*rho/sigma)*((Cstar/0.54d0)**(1.d0/3.d0)-1.d0)**2.d0)/g
      !   write(*,*) "soryujo",1.d0/0.49d0*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)*dm**2/htemp(I,J)**3.d0/g
      !   write(*,*) "water",nm**2*utemp(I,J)*sqrt(vtemp(I,J)**2+utemp(I,J)**2)/htemp(I,J)**(4.d0/3.d0)
      !   write(*,*) "htemp",htemp(I,J)
      !   write(*,*) "utemp",utemp(I,J)
      !   write(*,*) "vtemp",vtemp(I,J)
      !   write(*,*) "catemp",catemp(I,J)
      !   write(*,*) "at",I,J
      !   write(*,*) catemp(0,19),catemp(1,19),catemp(2,19)
      !   write(*,*) catemp(0,20),catemp(1,20),catemp(2,20)
      !   write(*,*) catemp(0,21),catemp(1,21),catemp(2,21)
      !   write(*,*) "iswet",iswet(I,J),htemp(I,J)
      !   call fwrite_xyz("h",h+zb-zbini)
      !   call fwrite_xyz("u",u)
      !   call fwrite_xyz("v",v)
      !   call fwrite_xyz("zb",zb)
      !   call fwrite_xyz("dzb",zb-zbini)
      !   call fwrite_xyz("ca",ca)
  
      !   read(*,*)
      ! endif
      if(iswet(I,J)==1.and.iswet(I-1,J)==1)then
        S0xtemp(I,J)=-(Zbtemp(I,J)-Zbtemp(I-1,J))/Dx
      else
        S0xtemp(I,J)=-(Zbtemp(I,J)-Zbtemp(I-1,J))/Dx !0.d0
      endif
      if(iswet(I,J)==1.and.iswet(I,J-1)==1)then
        S0ytemp(I,J)=-(Zbtemp(I,J)-Zbtemp(I,J-1))/Dy !predictorなので後進差分
      else
        S0ytemp(I,J)=-(Zbtemp(I,J)-Zbtemp(I,J-1))/Dy !0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*htemp(I,J)*abs(S0xtemp(I,J)+Sfxtemp(I,J)))*htemp(I,J)
      if(iswet(I+1,J)==1.and.iswet(I,J)==1.and.iswet(I-1,J)==1)then
        uhxx=(eps*(utemp(I+1,J)*htemp(I+1,J)-utemp(I,J)*htemp(I,J))/dx-eps*(utemp(I,J)*htemp(I,J)-utemp(I-1,J)*htemp(I-1,J))/dx)/dx
      else
        uhxx=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J)==1.and.iswet(I,J-1)==1)then
        uhyy=(eps*(utemp(I,J+1)*htemp(I,J+1)-utemp(I,J)*htemp(I,J))/dy-eps*(utemp(I,J)*htemp(I,J)-utemp(I,J-1)*htemp(I,J-1))/dy)/dy
      else
        uhyy=0.d0
      endif
      eps=1.d0/6.d0*karman*sqrt(g*htemp(I,J)*abs(S0ytemp(I,J)+Sfytemp(I,J)))*htemp(I,J)
      if(iswet(I+1,J)==1.and.iswet(I,J)==1.and.iswet(I-1,J)==1)then
        vhyy=(eps*(vtemp(I,J+1)*htemp(I,J+1)-vtemp(I,J)*htemp(I,J))/dy-eps*(vtemp(I,J)*htemp(I,J)-vtemp(I,J-1)*htemp(I,J-1))/dy)/dy
      else
        vhyy=0.d0
      endif
      if(iswet(I,J+1)==1.and.iswet(I,J)==1.and.iswet(I,J-1)==1)then
        vhxx=(eps*(vtemp(I+1,J)*htemp(I+1,J)-vtemp(I,J)*htemp(I,J))/dy-eps*(vtemp(I,J)*htemp(I,J)-vtemp(I-1,J)*htemp(I-1,J))/dx)/dx
      else
        vhxx=0.d0
      endif
      !calc water surface gradient !要検討
      ! if(utemp(I,J)>0.d0)then
      !   thewx=(-Zbtemp(I+1,J)-Htemp(I+1,J)+Zbtemp(I,J)+Htemp(I,J))/dx
      ! elseif(utemp(I,J)<0.d0)then
      !   thewx=(-Zbtemp(I,J)-Htemp(I,J)+Zbtemp(I-1,J)+Htemp(I-1,J))/dx
      ! else
      !   thewx=0.d0
      ! endif
      ! if(vtemp(I,J)>0.d0)then
      !   thewy=(-Zbtemp(I,J+1)-Htemp(I,J+1)+Zbtemp(I,J)+Htemp(I,J))/dy
      ! elseif(vtemp(I,J)<0.d0)then
      !   thewy=(-Zbtemp(I,J)-Htemp(I,J)+Zbtemp(I,J-1)+Htemp(I,J-1))/dy
      ! else
      !   thewy=0.d0
      ! endif
      ! if(utemp(I,J)**2*cos(atan(thewx))**2+vtemp(I,J)**2*cos(atan(thewy)**2)>0.d0)then
      !   thewxy=atan((sin(atan(thewx))*u(I,J)+sin(atan(thewy))*v(I,J))/&
      !   dsqrt(utemp(I,J)**2*cos(atan(thewx))**2+vtemp(I,J)**2*cos(atan(thewy)**2)))  !nakagawa
      ! else
      !   thewxy=0.d0
      ! endif
      !calc Cainf
      if(abs(tan(thewxy))>=tan(phi))then
        cainf=0.9d0*Cstar
      elseif(abs(tan(thewxy))>0.138d0)then
      ! elseif(catemp(I,J)>=0.4d0*cstar.and.htemp(I,J)/dm<30.d0)then
        cainf=min(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))),0.9d0*Cstar)
      elseif(abs(tan(thewxy))>0.03d0)then
      ! elseif(catemp(I,J)>0.01d0.and.htemp(I,J)/dm<30.d0)then
        cainf=6.7d0*(rho*abs(tan(thewxy))/((sigma-rho)*(tan(phi)-abs(tan(thewxy)))))**2
      else
        taustar=rho/(sigma-rho)*htemp(I,J)*tan(thewxy)/dm
        taustarc=0.04d0*10.d0**(1.72*tan(thewxy))
        alphc=sqrt(2.d0*(0.425d0-sigma*tan(thewxy)/(sigma-rho))/(1.d0-sigma*tan(thewxy)/(sigma-rho)))
        if(taustar>0.d0)then
          cainf=(1.d0+5.d0*tan(thewxy))*tan(thewxy)/(sigma/rho-1.d0)*&
          (1.d0-alphc**2*taustarc/taustar)*(1.d0-alphc*sqrt(taustarc/taustar))
        else
          cainf=0.d0
        endif
      ENDIF
      if(cainf>=Catemp(I,J))then
        if(Zbtemp(I,J)>Zbmin(I,J))then
          ierotemp(I,J)=min(delero*(cainf-Catemp(I,J))/(Cstar-cainf)*&
          sqrt(utemp(I,J)**2+vtemp(I,J)**2)*htemp(I,J)/dm,(Zbtemp(I,J)-Zbmin(I,J))/Dt)
          ! ierotemp(I,J)=delero*(cainf-Catemp(I,J))/(Cstar-cainf)*&
          ! sqrt(utemp(I,J)**2+vtemp(I,J)**2)*htemp(I,J)/dm
        else
          ierotemp(I,J)=0.d0
        endif
        ! write(*,*) "erosion",Cainf,Catemp(I,J),I,J,ierotemp(I,J),htemp(I,J)
        ! read(*,*)
      else
        ierotemp(I,J)=deldep*(cainf-Catemp(I,J))/Cstar*sqrt(utemp(I,J)**2+vtemp(I,J)**2)
      endif
      if(caltype==1)then
        Cp(I,J,1)=ierotemp(I,J)
      elseif(caltype==0)then
        Cp(I,J,1)=0.d0
      else
        write(*,*) "error"
      endif
      if(iswet(I,J)==1)then
        ! if((Up(I,J,2)+Dt*(g*htemp(I,J)*(S0xtemp(I,J)-Sfxtemp(I,J))+uhxx+uhyy))*Up(I,J,2)<0.d0)then !avoid reverse flow
        if((0.5d0*Up(I,J,2)+0.5d0*Uc(I,J,2)+0.5d0*Dt*(g*htemp(I,J)*(S0xtemp(I,J)-Sfxtemp(I,J))+uhxx+uhyy))*Up(I,J,2)<0.d0)then !avoid reverse flow
            Sfxtemp(I,J)=S0xtemp(I,J)-(-2.d0*Up(I,J,2)/Dt-uhxx-uhyy)/(g*htemp(I,J))
        endif
        ! if((Up(I,J,3)+Dt*(g*htemp(I,J)*(S0ytemp(I,J)-Sfytemp(I,J))+vhxx+vhyy))*Up(I,J,3)<0.d0)then !avoid reverse flow
        if((0.5d0*Up(I,J,3)+0.5d0*Uc(I,J,3)+0.5d0*Dt*(g*htemp(I,J)*(S0ytemp(I,J)-Sfytemp(I,J))+vhxx+vhyy))*Up(I,J,3)<0.d0)then !avoid reverse flow
            Sfytemp(I,J)=S0ytemp(I,J)-(-2.d0*Up(I,J,3)/Dt-vhxx-vhyy)/(g*htemp(I,J))
        endif
      endif
      Cp(I,J,2)=g*htemp(I,J)*(S0xtemp(I,J)-Sfxtemp(I,J))+uhxx+uhyy
      Cp(I,J,3)=g*htemp(I,J)*(S0ytemp(I,J)-Sfytemp(I,J))+vhxx+vhyy
      Cp(I,J,4)=ierotemp(I,J)*Cstar
      Cp(I,J,5)=-ierotemp(I,J)
      ! if(abs(ierotemp(I,J))>0.d0)then
      !   write(*,*) "ierotemp has an value"
      !   write(*,*) I,J,Catemp(I,J),ierotemp(I,J)
      !   read(*,*)
      ! endif
    ENDDO
    !$omp end parallel do
  ENDDO

  DO k=1,kmax
    DO J=1,Ny
      !$omp parallel do default(shared),private(I,kx,ky)
      DO I=1,Nx
        ! if(htemp(I-1,J)>0.d0.and.htemp(I+1,J)>0.d0)then
        ! if(I>=1.and.I<=Nx)then !上下流端には値がある
          ! kx=8.d0*Kv(K)*sqrt(g*htemp(I,J)*abs(S0xtemp(I,J)+Sfxtemp(I,J)))*htemp(I,J)/Dx
          kx=8.d0*Kv(K)*sqrt(g*htemp(I,J)*abs(Sfxtemp(I,J)))*htemp(I,J)/Dx
          ! if(K==5.and.Zbtemp(I,J)<Zbmin(I,J))then
          !   kx=0.d0
          ! endif
          Qxp(I,J,K)=0.125d0*kx*(Up(I+1,J,K)-2.d0*Up(I,J,K)+Up(I-1,J,K))
        ! else
          ! Qxp(I,J,K)=0.d0
        ! endif
        ! if(htemp(I,J-1)>0.d0.and.htemp(I,J+1)>0.d0)then
        ! if(J>1.and.J<Ny)then
          ! ky=8.d0*Kv(K)*sqrt(g*htemp(I,J)*abs(S0ytemp(I,J)+Sfytemp(I,J)))*htemp(I,J)/Dy
          ky=8.d0*Kv(K)*sqrt(g*htemp(I,J)*abs(Sfytemp(I,J)))*htemp(I,J)/Dy
          ! if(K==5.and.Zbtemp(I,J)<Zbmin(I,J))then
          !   ky=0.d0
          ! endif
          Qyp(I,J,K)=0.125d0*ky*(Up(I,J+1,K)-2.d0*Up(I,J,K)+Up(I,J+1,K))
        ! else
          ! Qyp(I,J,K)=0.d0
        ! endif
        ! if(k==5.and.Qxp(I,J,K)/=0.d0)then
        !   write(*,*) "aaaaa"
        !   write(*,*) I,J,K
        !   write(*,*) Qxp(I,J,K)
        !   write(*,*) 0.125d0*kx*(Uc(I+1,J,K)-2.d0*Uc(I,J,K)+Uc(I-1,J,K))
        !   write(*,*) 0.125d0*kx*(Up(I+1,J,K)-2.d0*Up(I,J,K)+Up(I-1,J,K))
        !   read(*,*)
        ! endif

      ENDDO
    !$omp end parallel do
    ENDDO
  ENDDO

end subroutine predictor
