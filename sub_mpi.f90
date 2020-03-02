subroutine mpi_boundary_2d_x(RANK,Tarray)
  use variables
  use mpi
  implicit none
  integer RANK,i1,i2,ierr
  integer,dimension(MPI_STATUS_SIZE) :: istat1
  real(8),dimension(0:Nx+1,0:Ny+1) :: Tarray
  real(8),dimension(Ny) :: send1d,send1d2,recv1d,recv1d2

  if(rank_l(RANK)<0)then !leftest
    call MPI_iRecv(recv1d(1),Ny,MPI_REAL8,rank_r(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Ny
      send1d(I)=Tarray(Nx,I)
    enddo
    call MPI_Send(send1d(1),Ny,MPI_REAL8,rank_r(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    do I=1,Ny
      Tarray(Nx+1,I)=recv1d(I)
    ENDDO
  elseif(rank_r(RANK)<0)then !rightest
    call MPI_iRecv(recv1d2(1),Ny,MPI_REAL8,rank_l(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Ny
      send1d2(I)=Tarray(1,I)
    enddo
    call MPI_Send(send1d2(1),Ny,MPI_REAL8,rank_l(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    DO I=1,Ny
      Tarray(0,I)=recv1d2(I)
    enddo
  else
    call MPI_iRecv(recv1d(1),Ny,MPI_REAL8,rank_r(RANK),0,MPI_COMM_WORLD,i1,ierr)
    call MPI_iRecv(recv1d2(1),Ny,MPI_REAL8,rank_l(RANK),0,MPI_COMM_WORLD,i2,ierr)
    DO I=1,Ny
      send1d(I)=Tarray(Nx,I)
      send1d2(I)=Tarray(1,I)
    enddo
    call MPI_Send(send1d(1),Ny,MPI_REAL8,rank_r(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Send(send1d2(1),Ny,MPI_REAL8,rank_l(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    call MPI_Wait(i2,Istat1,ierr)
    do I=1,Ny
      Tarray(Nx+1,I)=recv1d(I)
      Tarray(0,I)=recv1d2(I)
    enddo
  endif
!~   CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine mpi_boundary_2d_x

subroutine mpi_boundary_2d_y(RANK,Tarray)
  use variables
  use mpi
  implicit none
  integer RANK,i1,i2,ierr
  real(8),dimension(0:Nx+1,0:Ny+1) :: Tarray
  real(8),dimension(Nx) :: send1d,recv1d,send1d2,recv1d2
  integer,dimension(MPI_STATUS_SIZE) :: istat1


  if(rank_u(RANK)<0)then !top
    call MPI_iRecv(recv1d(1),Nx,MPI_REAL8,rank_d(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Nx
      send1d(I)=Tarray(I,Ny)
    enddo
    call MPI_Send(send1d(1),Nx,MPI_REAL8,rank_d(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    do I=1,Nx
      Tarray(I,Ny+1)=recv1d(I)
    ENDDO
  elseif(rank_d(RANK)<0)then !bottom
    call MPI_iRecv(recv1d2(1),Nx,MPI_REAL8,rank_u(RANK),0,MPI_COMM_WORLD,I1,ierr)
    DO I=1,Nx
      send1d2(I)=Tarray(I,1)
    enddo
    call MPI_Send(send1d2(1),Nx,MPI_REAL8,rank_u(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    DO I=1,Nx
      Tarray(I,0)=recv1d2(I)
    ENDDO
  else
    call MPI_iRecv(recv1d(1),Nx,MPI_REAL8,rank_d(RANK),0,MPI_COMM_WORLD,I1,ierr)
    call MPI_iRecv(recv1d2(1),Nx,MPI_REAL8,rank_u(RANK),0,MPI_COMM_WORLD,I2,ierr)
    DO I=1,Nx
      send1d(I)=Tarray(I,Ny)
      send1d2(I)=Tarray(I,1)
    enddo
    call MPI_Send(send1d(1),Nx,MPI_REAL8,rank_d(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Send(send1d2(1),Nx,MPI_REAL8,rank_u(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I2,Istat1,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    do I=1,Nx
      Tarray(I,Ny+1)=recv1d(I)
      Tarray(I,0)=recv1d2(I)
    ENDDO
  endif
!~   CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  
end subroutine mpi_boundary_2d_y



subroutine mpi_boundary_int_2d_x(RANK,Tarray)
  use variables
  use mpi
  implicit none
  integer RANK,i1,i2,ierr
  integer,dimension(MPI_STATUS_SIZE) :: istat1
  integer,dimension(0:Nx+1,0:Ny+1) :: Tarray
  integer,dimension(Ny) :: send1d,send1d2,recv1d,recv1d2

  if(rank_l(RANK)<0)then !leftest
    call MPI_iRecv(recv1d(1),Ny,MPI_INTEGER,rank_r(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Ny
      send1d(I)=Tarray(Nx,I)
    enddo
    call MPI_Send(send1d(1),Ny,MPI_INTEGER,rank_r(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    do I=1,Ny
      Tarray(Nx+1,I)=recv1d(I)
    ENDDO
  elseif(rank_r(RANK)<0)then !rightest
    call MPI_iRecv(recv1d2(1),Ny,MPI_INTEGER,rank_l(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Ny
      send1d2(I)=Tarray(1,I)
    enddo
    call MPI_Send(send1d2(1),Ny,MPI_INTEGER,rank_l(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    DO I=1,Ny
      Tarray(0,I)=recv1d2(I)
    enddo
  else
    call MPI_iRecv(recv1d(1),Ny,MPI_INTEGER,rank_r(RANK),0,MPI_COMM_WORLD,i1,ierr)
    call MPI_iRecv(recv1d2(1),Ny,MPI_INTEGER,rank_l(RANK),0,MPI_COMM_WORLD,i2,ierr)
    DO I=1,Ny
      send1d(I)=Tarray(Nx,I)
      send1d2(I)=Tarray(1,I)
    enddo
    call MPI_Send(send1d(1),Ny,MPI_INTEGER,rank_r(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Send(send1d2(1),Ny,MPI_INTEGER,rank_l(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(i1,Istat1,ierr)
    call MPI_Wait(i2,Istat1,ierr)
    do I=1,Ny
      Tarray(Nx+1,I)=recv1d(I)
      Tarray(0,I)=recv1d2(I)
    enddo
  endif
!~   CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine mpi_boundary_int_2d_x

subroutine mpi_boundary_int_2d_y(RANK,Tarray)
  use variables
  use mpi
  implicit none
  integer RANK,i1,i2,ierr
  integer,dimension(0:Nx+1,0:Ny+1) :: Tarray
  integer,dimension(Nx) :: send1d,recv1d,send1d2,recv1d2
  integer,dimension(MPI_STATUS_SIZE) :: istat1


  if(rank_u(RANK)<0)then !top
    call MPI_iRecv(recv1d(1),Nx,MPI_INTEGER,rank_d(RANK),0,MPI_COMM_WORLD,i1,ierr)
    DO I=1,Nx
      send1d(I)=Tarray(I,Ny)
    enddo
    call MPI_Send(send1d(1),Nx,MPI_INTEGER,rank_d(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    do I=1,Nx
      Tarray(I,Ny+1)=recv1d(I)
    ENDDO
  elseif(rank_d(RANK)<0)then !bottom
    call MPI_iRecv(recv1d2(1),Nx,MPI_INTEGER,rank_u(RANK),0,MPI_COMM_WORLD,I1,ierr)
    DO I=1,Nx
      send1d2(I)=Tarray(I,1)
    enddo
    call MPI_Send(send1d2(1),Nx,MPI_INTEGER,rank_u(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    DO I=1,Nx
      Tarray(I,0)=recv1d2(I)
    ENDDO
  else
    call MPI_iRecv(recv1d(1),Nx,MPI_INTEGER,rank_d(RANK),0,MPI_COMM_WORLD,I1,ierr)
    call MPI_iRecv(recv1d2(1),Nx,MPI_INTEGER,rank_u(RANK),0,MPI_COMM_WORLD,I2,ierr)
    DO I=1,Nx
      send1d(I)=Tarray(I,Ny)
      send1d2(I)=Tarray(I,1)
    enddo
    call MPI_Send(send1d(1),Nx,MPI_INTEGER,rank_d(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Send(send1d2(1),Nx,MPI_INTEGER,rank_u(RANK),0,MPI_COMM_WORLD,ierr)
    call MPI_Wait(I2,Istat1,ierr)
    call MPI_Wait(I1,Istat1,ierr)
    do I=1,Nx
      Tarray(I,Ny+1)=recv1d(I)
      Tarray(I,0)=recv1d2(I)
    ENDDO
  endif
!~   CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  
end subroutine mpi_boundary_int_2d_y
