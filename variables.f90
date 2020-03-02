module variables
  implicit none
  real(8),allocatable,dimension(:,:,:) :: Uc,Ec,Fc,Cc,Up,Ep,Fp,Cp,Qxc,Qyc,Qxp,Qyp,Uctemp
  real(8),allocatable,dimension(:,:) :: u,v,h,S0x,Sfx,S0y,Sfy,xix,xiy,zb,zbini,zbmin,qb,ca,iero,&
  htemp,utemp,vtemp,catemp,zbtemp,sfxtemp,sfytemp,s0xtemp,s0ytemp,ierotemp,wlmax
  real(8),allocatable,dimension(:) :: TVDx,TVDy
  real(8),parameter :: &
  nm=0.07d0,&
  karman=0.4d0,&
  g=9.80665d0,&
  Tmax=7201.d0,&
  Hcr=0.05d0,&
  sli=18.d0/180.d0*4.d0*atan(1.d0),& !=0.324919696232906d0,& !0.04d0,& !0.324919696232906d0,& !0.04d0,&
  slidown=2.d0/180.d0*4.d0*atan(1.d0),& !0.324919696232906d0,& !0.04d0,& !0.06992681194351d0,&!  0.04d0 !0.086749973905571d0 !=2.d0
  qup=0.0001d0,&!0.0136d0,& !0.005408d0,& !0.0136d0,&
  hup=0.02d0,&!0.01d0 ,&!0.021790608238721d0,& !0.02179d0,& !0.01162311d0,& !0.021790608238721d0,&
  cup=0.3d0,&!0.d0,& !0.35d0,&!0.35d0,& 
  hdown=0.01d0,& !0.021790608238721d0*10.d0**(3.d0/5.d0),& !0.021790608238721d0,&
  hini=0.1d0,&
  cini=0.5d0,&
  Dtini=1.d0/2.d0**10,&
  Dtmin=1.d0/2.d0**15,&
  dtkeeptime=0.d0,& !2.d0,&
  Dtmax=0.25d0,& ! 1.d0,& !0.25d0
  fs_Dt=0.05d0,& !0.2d0 ,&!0.004d0,& 0.1d0
  rho=1.d0,&
  Con=0.4d0,&
  cstar=0.6d0,&
  sigma=2.65d0,&
  dm=0.05d0,& !0.00342d0,&
  phi=35.d0/180.d0*4.d0*atan(1.d0),&
  deldep=0.05d0,&
  delero=0.0007d0
! nm=0.05d0,&
! karman=0.4d0,&
! g=9.80665d0,&
! Tmax=1200.d0,&
! Hcr=0.002d0,&
! sli=18.d0/180.d0*4.d0*atan(1.d0),& !=0.324919696232906d0,& !0.04d0,& !0.324919696232906d0,& !0.04d0,&
! slidown=2.d0/180.d0*4.d0*atan(1.d0),& !0.324919696232906d0,& !0.04d0,& !0.06992681194351d0,&!  0.04d0 !0.086749973905571d0 !=2.d0
! qup=0.0015d0,&!0.0136d0,& !0.005408d0,& !0.0136d0,&
! hup=0.02d0,&!0.01d0 ,&!0.021790608238721d0,& !0.02179d0,& !0.01162311d0,& !0.021790608238721d0,&
! cup=0.3d0,&!0.d0,& !0.35d0,&!0.35d0,& 
! hdown=0.01d0,& !0.021790608238721d0*10.d0**(3.d0/5.d0),& !0.021790608238721d0,&
! hini=0.1d0,&
! cini=0.3d0,&
! Dtini=1.d0/2.d0**10,&
! rho=1.d0,&
! Con=0.4d0,&
! cstar=0.6d0,&
! sigma=2.65d0,&
! dm=0.00342d0,&
! phi=36.d0/180.d0*4.d0*atan(1.d0),&
! deldep=0.05d0,&
! delero=0.0007d0


  real(8),dimension(5) ::kv=(/2.5d0,2.5d0,2.5d0,1.d0,1.d0/) !=2.d0
  ! real(8),dimension(5) ::kv=(/0.d0,0.d0,0.d0,1.d0,1.d0/) !=2.d0
  ! real(8),dimension(5) ::kv=(/10.d0,10.d0,10.d0,1.d0,1.d0/) !=2.d0

  real(8) :: T,Dt,Bini,Cfl,Dx,Dy,xllcorner,yllcorner
  !Dx=0.02d0,Dy=0.02d0
  integer I,J,K,step,kmax
  integer,allocatable,dimension(:,:) :: iswet
  integer,parameter :: caltype=1,& !1:debrisflow,0:water
  boundtype=7,& !1:水槽，type=2:下流端開放水路，定常供給,3:跳水実験,4:土石流
  inittype=6 !1:dambrake,2:一定水深,3:跳水実験,4:土石流,5:file_read
  integer :: Nx,Ny !Nx=375,Ny=95
  integer my_rank,PETOT
  integer,dimension(:),allocatable,save::rank_r,rank_u,rank_l,rank_d
  integer :: N_mpi_x=90,N_mpi_y=320
  ! integer :: N_mpi_x=75,N_mpi_y=96
  ! integer :: N_mpi_x=2,N_mpi_y=2
  ! deldep=0.d0,delero=0.d0
  integer :: iserror=0
end module variables

!parameter for flume inudnnation
! nm=0.05d0,&
! karman=0.4d0,&
! g=9.80665d0,&
! Tmax=1200.d0,&
! Hcr=0.002d0,&
! sli=18.d0/180.d0*4.d0*atan(1.d0),& !=0.324919696232906d0,& !0.04d0,& !0.324919696232906d0,& !0.04d0,&
! slidown=2.d0/180.d0*4.d0*atan(1.d0),& !0.324919696232906d0,& !0.04d0,& !0.06992681194351d0,&!  0.04d0 !0.086749973905571d0 !=2.d0
! qup=0.0001d0,&!0.0136d0,& !0.005408d0,& !0.0136d0,&
! hup=0.02d0,&!0.01d0 ,&!0.021790608238721d0,& !0.02179d0,& !0.01162311d0,& !0.021790608238721d0,&
! cup=0.3d0,&!0.d0,& !0.35d0,&!0.35d0,& 
! hdown=0.01d0,& !0.021790608238721d0*10.d0**(3.d0/5.d0),& !0.021790608238721d0,&
! hini=0.1d0,&
! cini=0.3d0,&
! Dtini=1.d0/2.d0**8,&
! rho=1.d0,&
! Con=0.4d0,&
! cstar=0.6d0,&
! sigma=2.65d0,&
! dm=0.00342d0,&
! phi=35.d0/180.d0*4.d0*atan(1.d0),&
! deldep=0.05d0,&
! delero=0.0007d0

