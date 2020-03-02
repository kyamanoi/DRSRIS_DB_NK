program maccormack
  !$ use omp_lib
  use variables
  use mpi
  implicit none
  Integer ierr,ii,ranktemp,lefttemp,righttemp,uptemp,downtemp
  character dummy
  character(1000) ftemp,dirname,command
  !$ real(8) :: time(10),timetot(10)
  ! allocate(Uc(Nx,Ny,3),Ec(Nx,Ny,3),Fc(Nx,Ny,3),Cc(Nx,Ny,3),Up(Nx,Ny,3),Ep(Nx,Ny,3),Fp(Nx,Ny,3),Cp(Nx,Ny,3))
  ! allocate(u(Nx,Ny),v(Nx,Ny),h(Nx,Ny),S0x(Nx,Ny),Sfx(Nx,Ny),S0y(Nx,Ny),Sfy(Nx,Ny))
  ! open(1001,file='./input_hanransuiro/dem.asc')
  ! open(1002,file='./input_dem/ldepth.asc')
  dirname="."
  call MPI_INIT    (ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr )
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr )
  write(*,*) "PETOT=",PETOT
  allocate(rank_r(0:PETOT-1),rank_u(0:PETOT-1),rank_l(0:PETOT-1),rank_d(0:PETOT-1))

  ! write(command,"('mkdir -p ',A,'/output/rank',i4.4,'/asc')")trim(dirname),my_rank
  write(command,"('mkdir -p ',A,'/output/rank',i5.5,'/asc')")trim(dirname),my_rank
  call system(command)
  ! write(command,"('mkdir -p ',A,'/output/rank',i4.4,'/xyz')")trim(dirname),my_rank
  write(command,"('mkdir -p ',A,'/output/rank',i5.5,'/xyz')")trim(dirname),my_rank
  call system(command)
  
  open(1000,file="around_rank.dat",STATUS='old',action='read')
  read(1000,*)
  DO II=0,PETOT-1
    read(1000,*) ranktemp,lefttemp,righttemp,uptemp,downtemp
    ! if(my_rank==ranktemp)then
      rank_l(II)=lefttemp
      rank_r(II)=righttemp
      rank_u(II)=uptemp
      rank_d(II)=downtemp
    ! endif
  enddo

  !   do ii=0,PETOT-1
  !   if(mod(ii+1,N_mpi_x)==0)then
  !     rank_r(ii)=-1
  !   else
  !     rank_r(ii)=ii+1
  !   endif
  !   if(mod(ii+1,N_mpi_x)==1)then
  !     rank_l(ii)=-1
  !   else
  !     rank_l(ii)=ii-1
  !   endif
  !   if(ii<N_mpi_x)then
  !     rank_u(ii)=-1
  !   else
  !     rank_u(ii)=ii-N_mpi_x
  !   endif
  !   if(ii>N_mpi_x*(N_mpi_y-1)-1)then
  !     rank_d(ii)=-1
  !   else
  !     rank_d(ii)=ii+N_mpi_x
  !   endif
  ! enddo

  ! write(ftemp,"( A,'/input/rank',I4.4,'/dem.asc')") trim(dirname),my_rank
  write(ftemp,"( A,'/input/rank',I5.5,'/dem.asc')") trim(dirname),my_rank
  ! write(ftemp,"( A,'/input_5mres/rank',I4.4,'/dem.asc')") trim(dirname),my_rank
  open(1001,file=ftemp,STATUS='old',action='read')
  ! write(ftemp,"( A,'/input/rank',I4.4,'/ldepth.asc')") trim(dirname),my_rank
  write(ftemp,"( A,'/input/rank',I5.5,'/ldepth.asc')") trim(dirname),my_rank
  ! write(ftemp,"( A,'/input/rank',I4.4,'/ldepth_points.asc')") trim(dirname),my_rank
  open(1002,file=ftemp,STATUS='old',action='read')
  

  read(1001,*) dummy, Nx
  read(1001,*) dummy, Ny
  read(1001,*) dummy, xllcorner
  read(1001,*) dummy, yllcorner
  read(1001,*) dummy, dx
  dy=dx

  allocate(Uc(0:Nx+1,0:Ny+1,5),Ec(0:Nx+1,0:Ny+1,5),Fc(0:Nx+1,0:Ny+1,5),Cc(0:Nx+1,0:Ny+1,5),Up(0:Nx+1,0:Ny+1,5),Ep(0:Nx+1,0:Ny+1,5)&
  ,Fp(0:Nx+1,0:Ny+1,5),Cp(0:Nx+1,0:Ny+1,5),Qxc(0:Nx+1,0:Ny+1,5),Qyc(0:Nx+1,0:Ny+1,5),Qxp(0:Nx+1,0:Ny+1,5),Qyp(0:Nx+1,0:Ny+1,5),&
  Uctemp(0:Nx+1,0:Ny+1,5))
  allocate(u(0:Nx+1,0:Ny+1),v(0:Nx+1,0:Ny+1),h(0:Nx+1,0:Ny+1),S0x(0:Nx+1,0:Ny+1),Sfx(0:Nx+1,0:Ny+1),S0y(0:Nx+1,0:Ny+1),&
  Sfy(0:Nx+1,0:Ny+1),xix(0:Nx+1,0:Ny+1),xiy(0:Nx+1,0:Ny+1),zb(0:Nx+1,0:Ny+1),qb(0:Nx+1,0:Ny+1),ca(0:Nx+1,0:Ny+1),&
  iero(0:Nx+1,0:Ny+1),zbini(0:Nx+1,0:Ny+1),Zbmin(0:Nx+1,0:Ny+1))
  allocate(TVDx(1:Nx),TVDy(1:Ny))
  allocate(iswet(0:Nx+1,0:Ny+1))
  allocate(htemp(0:Nx+1,0:Ny+1),utemp(0:Nx+1,0:Ny+1),vtemp(0:Nx+1,0:Ny+1),catemp(0:Nx+1,0:Ny+1),zbtemp(0:Nx+1,0:Ny+1),&
  sfxtemp(0:Nx+1,0:Ny+1),sfytemp(0:Nx+1,0:Ny+1),s0xtemp(0:Nx+1,0:Ny+1),s0ytemp(0:Nx+1,0:Ny+1),ierotemp(0:Nx+1,0:Ny+1))
  allocate(wlmax(0:Nx+1,0:Ny+1))
  ! kv(1:3)=2.5d0 !2.5d0 !10.d0
  ! kv(4:5)=1.d0
  Dt=Dtini
  Bini=Ny*Dy
  T=0.d0
  step=1

  if(caltype==1)then
    kmax=5
  elseif(caltype==0)then
    kmax=3
  endif
  call Initial_condition
  if(PETOT>1)then
    do K=1,kmax
      call mpi_boundary_2d_x(my_rank,Uc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Uc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_x(my_rank,Ec(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Ec(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_x(my_rank,Fc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Fc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_x(my_rank,Cc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Cc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_x(my_rank,Qxc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Qxc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_x(my_rank,Qyc(0:Nx+1,0:Ny+1,K))
      call mpi_boundary_2d_y(my_rank,Qyc(0:Nx+1,0:Ny+1,K))
    enddo
    call mpi_boundary_int_2d_x(my_rank,iswet(0:Nx+1,0:Ny+1))
    call mpi_boundary_int_2d_y(my_rank,iswet(0:Nx+1,0:Ny+1))
  endif
  wlmax=0.d0
  DO while (T<Tmax)
    !1:backward, 0:forward
    !$ time(1) = omp_get_wtime()
    call boundary_condition_c
    !$ time(2) = omp_get_wtime()
    !$ timetot(2) =timetot(2)+time(2)-time(1)
    ! call TVD
    ! if(step==1)then
    !   call predictor(0.d0,0.d0)
    !   call corrector(1.d0,1.d0)
    ! elseif(step==2)then
    !   call predictor(0.d0,1.d0)
    !   call corrector(1.d0,0.d0)
    ! elseif(step==3)then
    !   call predictor(1.d0,1.d0)
    !   call corrector(0.d0,0.d0)
    ! elseif(step==4)then
    !   call predictor(1.d0,0.d0)
    !   call corrector(0.d0,1.d0)
    ! endif
    ! if(step<4)then
    !   step=step+1
    ! else
    !   step=1
    ! endif
    call predictor(1.d0,1.d0)
    !$ time(3) = omp_get_wtime()
    !$ timetot(3) =timetot(3)+time(3)-time(2)
    ! call check_cfl(0.5d0)
    if(PETOT>1)then
      do K=1,kmax
        call mpi_boundary_2d_x(my_rank,Up(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Up(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Ep(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Ep(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Fp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Fp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Cp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Cp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Qxp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Qxp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Qyp(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Qyp(0:Nx+1,0:Ny+1,K))
      enddo
      call mpi_boundary_int_2d_x(my_rank,iswet(0:Nx+1,0:Ny+1))
      call mpi_boundary_int_2d_y(my_rank,iswet(0:Nx+1,0:Ny+1))
    endif
    call boundary_condition_p
    !$ time(4) = omp_get_wtime()
    !$ timetot(4) =timetot(4)+time(4)-time(3)
    call corrector(0.d0,0.d0)
    !$ time(5) = omp_get_wtime()
    !$ timetot(5) =timetot(5)+time(5)-time(4)
    if(PETOT>1)then
      do K=1,kmax
        call mpi_boundary_2d_x(my_rank,Uc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Uc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Ec(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Ec(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Fc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Fc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Cc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Cc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Qxc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Qxc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_x(my_rank,Qyc(0:Nx+1,0:Ny+1,K))
        call mpi_boundary_2d_y(my_rank,Qyc(0:Nx+1,0:Ny+1,K))
      enddo
      call mpi_boundary_int_2d_x(my_rank,iswet(0:Nx+1,0:Ny+1))
      call mpi_boundary_int_2d_y(my_rank,iswet(0:Nx+1,0:Ny+1))
    endif
    ! call check_cfl(fs_Dt)
    where(wlmax(0:Nx+1,0:Ny+1)<Uc(0:Nx+1,0:Ny+1,1))wlmax(0:Nx+1,0:Ny+1)=Uc(0:Nx+1,0:Ny+1,1)
    !$ time(6) = omp_get_wtime()
    !$ timetot(6) =timetot(6)+time(6)-time(5)
    call IO
    call check_cfl(fs_Dt)
    !$ time(7) = omp_get_wtime()
    !$ timetot(7) =timetot(7)+time(7)-time(6)
    ! write(*,*) "h=",h(2,10)
    ! write(*,*) "qx=",u(2,10)*h(2,10)
    ! write(*,*) "Ca=",Ca(2,10)
    ! write(*,*) "Zb=",Zb(2,10)
    ! read(*,*)
    T=T+Dt
  ENDDO

  contains
  subroutine IO
    if(T==0)then
    ! write(ftemp,"( A,'/output/rank',I4.4,'/Hcenter.dat')") trim(dirname),my_rank
    ! open(2001,file=ftemp)
    ! write(ftemp,"( A,'/output/rank',I4.4,'/Ucenter.dat')") trim(dirname),my_rank
    ! open(2002,file=ftemp)
    ! write(ftemp,"( A,'/output/rank',I4.4,'/Vcenter.dat')") trim(dirname),my_rank
    ! open(2003,file=ftemp)
    ! write(ftemp,"( A,'/output/rank',I4.4,'/Ccenter.dat')") trim(dirname),my_rank
    ! open(2004,file=ftemp)
    ! write(ftemp,"( A,'/output/rank',I4.4,'/Zbcenter.dat')") trim(dirname),my_rank
    ! open(2005,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/h.asc')") trim(dirname),my_rank
      write(ftemp,"( A,'/output/rank',I5.5,'/h.asc')") trim(dirname),my_rank
      open(2101,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/u.asc')") trim(dirname),my_rank
      ! write(ftemp,"( A,'/output/rank',I5.5,'/u.asc')") trim(dirname),my_rank
      ! open(2102,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/v.asc')") trim(dirname),my_rank
      ! write(ftemp,"( A,'/output/rank',I5.5,'/v.asc')") trim(dirname),my_rank
      ! open(2103,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/zb.asc')") trim(dirname),my_rank
      ! write(ftemp,"( A,'/output/rank',I5.5,'/zb.asc')") trim(dirname),my_rank
      ! open(2104,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/dzb.asc')") trim(dirname),my_rank
      write(ftemp,"( A,'/output/rank',I5.5,'/dzb.asc')") trim(dirname),my_rank
      open(2105,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/ca.asc')") trim(dirname),my_rank
      write(ftemp,"( A,'/output/rank',I5.5,'/ca.asc')") trim(dirname),my_rank
      open(2106,file=ftemp)
      ! write(ftemp,"( A,'/output/rank',I4.4,'/wlmax.asc')") trim(dirname),my_rank
      write(ftemp,"( A,'/output/rank',I5.5,'/wlmax.asc')") trim(dirname),my_rank
      open(2107,file=ftemp)
    endif
    if(mod(T,300.d0)<mod(T-Dt,300.d0))then
    ! if(mod(T,1.d0)<mod(T-Dt,1.d0))then
        ! if(mod(T+Dt*0.00001d0,1.d0)<Dt*0.1d0)then
      call fwrite_asc_time(2101,h)
      ! call fwrite_asc_time(2102,u)
      ! call fwrite_asc_time(2103,v)
      ! call fwrite_asc_time(2104,zb)
      call fwrite_asc_time(2105,zb-zbini)
      call fwrite_asc_time(2106,ca)
      call fwrite_asc_time(2107,wlmax)
      ! call fwrite_asc("h",h)!+zb-zbini)
      ! call fwrite_asc("u",u)
      ! call fwrite_asc("v",v)
      ! call fwrite_asc("zb",zb)
      ! call fwrite_asc("dzb",zb-zbini)
      ! call fwrite_asc("ca",ca)
      ! call fwrite_asc("eta",h+zb)
      ! call fwrite_asc("s0x",s0x)
      ! call fwrite_asc("s0y",s0y)
      ! call fwrite_asc("sfx",sfx)
      ! call fwrite_asc("sfy",sfy)
      ! call fwrite_asc_int("iw",iswet)
      ! call fwrite_asc("htemp",htemp)!+zbtemp-zbini)
      ! call fwrite_asc("utemp",utemp)
      ! call fwrite_asc("vtemp",vtemp)
      ! call fwrite_asc("zbtemp",zbtemp)
      ! call fwrite_asc("dzbtemp",zbtemp-zbini)
      ! call fwrite_asc("catemp",catemp)
      ! call fwrite_asc("s0xtemp",s0xtemp)
      ! call fwrite_asc("s0ytemp",s0ytemp)
      ! call fwrite_asc("sfxtemp",sfxtemp)
      ! call fwrite_asc("sfytemp",sfytemp)
      ! call fwrite_asc("etatemp",htemp+zbtemp)
    endif
    ! if(mod(T+Dt*0.00001d0,1.d0)<Dt*0.1d0)then
      ! call fwrite_xyz("h",h) !+zb-zbini)
      ! call fwrite_xyz("u",u)
      ! call fwrite_xyz("v",v)
      ! call fwrite_xyz("zb",zb)
      ! call fwrite_xyz("dzb",zb-zbini)
      ! call fwrite_xyz("ca",ca)
    ! endif
    if(mod(T,1.d0)<mod(T-Dt,1.d0))then
    ! if(mod(T+Dt*0.00001d0,1.d0)<Dt*0.001d0)then
      ! write(2001,'(f16.6,3500f16.6)') T,(h(I,int(Ny/2)+1),I=1,Nx)
      ! write(2002,'(f16.6,3500f16.6)') T,(u(I,int(Ny/2)+1),I=1,Nx)
      ! write(2003,'(f16.6,3500f16.6)') T,(v(I,int(Ny/2)+1),I=1,Nx)
      ! write(2004,'(f16.6,3500f16.6)') T,(ca(I,int(Ny/2)+1),I=1,Nx)
      ! write(2005,'(f16.6,3500f16.6)') T,(zb(I,int(Ny/2)+1),I=1,Nx)
      if(my_rank==0)then
        write(*,*) T,"/",Tmax
        write(*,*) "Dt=1/2^",log(1.d0/dt)/log(2.d0)
        write(*,*) "cfl=",Cfl
        ! write(*,*) "total_water=",sum(Uc(:,:,1))
        !$ write(*,*) "timetimetimetime"
        !$ write(*,*) "boundaryC",timetot(2)
        !$ write(*,*) "predictor",timetot(3)
        !$ write(*,*) "boundaryP",timetot(4)
        !$ write(*,*) "corrector",timetot(5)
        !$ write(*,*) "IO",timetot(6)
      endif
    endif
  end subroutine IO
end program maccormack

! subroutine TVD
!   use variables
!   implicit none
!   real(8) :: nu,rplus,rminus
  
!   nu=maxval(u+sqrt(g*h))*Dt/Dx
!   DO K=1,3
!     DO J=1,Ny
!       DO I=1,Nx
!         rplus=Uc(I,J,K)
!       enddo
!     ENDDO
!   ENDDO
  
!   contains
!   function gpm(rin,nuin)
!     real(8) gpm,rin,nuin
!     if(rin>=0.d0)then
!       gpm=0.5d0*nuin*(1.d0-abs(nuin)*(1.d0-min(2.d0*rin,1.d0))) 
!     else
!       gpm=0.5d0*nuin*(1.d0-abs(nuin)) 
!     endif
!   end function gpm

! end subroutine TVD
