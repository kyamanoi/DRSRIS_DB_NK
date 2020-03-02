subroutine check_cfl(fs)
  use variables
  use mpi
  implicit none
  real(8) Dtcr,maxvel,maxvelall
  real(8),intent(in) :: fs
  integer,dimension(2) :: mloc
  integer ierr

  ! Dtcr=Dx/(maxval(sqrt(u(1:Nx,1:Ny)**2+v(1:Nx,1:Ny)**2)+sqrt(2.d0*g*h(1:Nx,1:Ny)))) !for single node

  maxvel=(maxval(sqrt(u(1:Nx,1:Ny)**2+v(1:Nx,1:Ny)**2)+sqrt(2.d0*g*h(1:Nx,1:Ny))))
    if(maxvel>Dx/Dtmin*Fs)then
      write(*,*) "maxvel modified!!,rank=",my_rank,"T=",T
      maxvel=0.d0
      iserror=1
      u=0.d0
      v=0.d0
      h=0.d0
      Ca=0.d0
      Uc(:,:,1:4)=0.d0
      Ec(:,:,1:4)=0.d0
      Fc(:,:,1:4)=0.d0
      Cc=0.d0
      Qxc=0.d0
      Qyc=0.d0
      iswet=0.d0
      Up(:,:,1:4)=0.d0
      Ep(:,:,1:4)=0.d0
      Fp(:,:,1:4)=0.d0
      Cp=0.d0
      Qxp=0.d0
      Qyp=0.d0
    endif
  if(PETOT>1)then
    call mpi_allreduce(maxvel,maxvelall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  else
    maxvelall=maxvel
  endif
  Dtcr=Dx/maxvelall
  Cfl=Dt/Dtcr

  if(T>dtkeeptime)then
    if(Dt>fs*Dtcr)then
      do while (Dt>Dtcr*Fs)
        ! write(*,*) "dt become 1/2"
        Dt=Dt/2.d0
        if(Dt<Dtmin)then
          if(maxvel==maxvelall)then
            write(*,*) "too small dt. rank=",my_rank,Dt
            write(*,*) "maxvel=",maxvelall,maxvel
            mloc=maxloc(sqrt(u(1:Nx,1:Ny)**2+v(1:Nx,1:Ny)**2)+sqrt(2.d0*g*h(1:Nx,1:Ny)))
            write(*,*) my_rank,mloc(1),mloc(2)
            write(*,*) "                        ",u(mloc(1),mloc(2) - 1),"                "
            write(*,*) u(mloc(1)-1,mloc(2)),u(mloc(1),mloc(2)),u(mloc(1)+1,mloc(2))
            write(*,*) "                        ",u(mloc(1),mloc(2) + 1),"                "
            write(*,*) 
            write(*,*) "                        ",v(mloc(1),mloc(2) - 1),"                "
            write(*,*) v(mloc(1)-1,mloc(2)),v(mloc(1),mloc(2)),v(mloc(1)+1,mloc(2))
            write(*,*) "                        ",v(mloc(1),mloc(2) + 1),"                "
            write(*,*) 
            write(*,*) "                        ",h(mloc(1),mloc(2) - 1),"                "
            write(*,*) h(mloc(1)-1,mloc(2)),h(mloc(1),mloc(2)),h(mloc(1)+1,mloc(2))
            write(*,*) "                        ",h(mloc(1),mloc(2) + 1),"                "
            write(*,*) 
            write(*,*) "                        ",sfx(mloc(1),mloc(2) - 1),"                "
            write(*,*) sfx(mloc(1)-1,mloc(2)),sfx(mloc(1),mloc(2)),sfx(mloc(1)+1,mloc(2))
            write(*,*) "                        ",sfx(mloc(1),mloc(2) + 1),"                "
            write(*,*) 
            write(*,*) "                        ",sfy(mloc(1),mloc(2) - 1),"                "
            write(*,*) sfy(mloc(1)-1,mloc(2)),sfy(mloc(1),mloc(2)),sfy(mloc(1)+1,mloc(2))
            write(*,*) "                        ",sfy(mloc(1),mloc(2) + 1),"                "
          endif
          ! stop
        endif
      enddo
    endif

    if(Dt<fs*Dtcr*0.5)then
      do while (Dt<fs*Dtcr*0.5)
        ! write(*,*) "dt become 2"
        Dt=Dt*2.d0
      enddo
    endif
    if(DT>Dtmax)then
      Dt=Dtmax
    endif
  endif

  ! if(Dt>fs*Dtcr)then
  !   write(*,*) "Dt is too large!!!"
  !   write(*,*) "cfl=",Cfl
  !   write(*,*) Dtcr,Dt,fs,T
  !   write(*,*) maxval(sqrt(u**2+v**2)),maxval(sqrt(2.d0*g*h))
  !   write(*,*) maxval(u(1:Nx,1:Ny)),maxval(v(1:Nx,1:Ny)),maxval(h(1:Nx,1:Ny))
  !   call fwrite_xyz("h",h)
  !   call fwrite_xyz("u",u)
  !   call fwrite_xyz("v",v)
  !   call fwrite_xyz("zb",zb)
  !   call fwrite_xyz("dzb",zb-zbini)
  !   call fwrite_xyz("ca",ca)
  !   write(*,*) "file wrote"
  !   ! read(*,*)
  !   do while (Dt>Dtcr*Fs)
  !     write(*,*) "dt become 1/2"
  !     Dt=Dt/2.d0
  !   enddo
  ! endif

end subroutine check_cfl
  