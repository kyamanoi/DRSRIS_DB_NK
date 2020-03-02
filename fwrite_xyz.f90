
subroutine fwrite_xyz(name,outarray)
  use variables
  implicit none
  character(*) name
  character(len=1000) fname,timechar
  character(len=13) fmt
  character(len=4) rankname
  ! real(8),dimension(Nx,Ny) :: outarray
  real(8),dimension(0:Nx+1,0:Ny+1) :: outarray
  real(8) xmargin,ymargin

  xmargin=mod(my_rank,n_mpi_x)*dx*Nx
  ymargin=int(my_rank/n_mpi_x)*dy*Ny
  write(rankname,'(I4.4)') my_rank
  
  write(timechar,'(F0.4)') T
  if(T==0.d0)timechar="0.0000"
  fname= "./output/rank"//rankname//"/xyz/"//trim(name)//trim(timechar)//".xyz"
  OPEN(5001,FILE=fname,STATUS='UNKNOWN')
  write(5001,'(A30)')"# X          Y         Z"

  ! do J=1,Ny
  !   DO I=1,Nx
  !     write(5001,'(3f16.6)') (real(I)-0.5d0)*Dx,(real(J)-0.5d0)*Dy,outarray(I,J)
  !   ENDDO
  ! enddo
  DO I=1,Nx
    do J=1,Ny
      ! write(5001,'(3f16.6)') (real(I)-0.5d0)*Dx,(real(J)-0.5d0)*Dy,outarray(I,J)
      write(5001,'(3f16.6)') (real(I)-0.5d0)*Dx+xmargin,(real(J)-0.5d0)*Dy+ymargin,outarray(I,J)
    ENDDO
    write(5001,*) 
  enddo

  flush(5001)
  close(5001)
end subroutine fwrite_xyz

subroutine fwrite_txyz(fnum,outarray)
  use variables
  implicit none
  integer,intent(in) :: fnum
  ! real(8),dimension(Nx,Ny) :: outarray
  real(8),dimension(0:Nx+1,0:Ny+1) :: outarray
  real(8) xmargin,ymargin

  xmargin=mod(my_rank,n_mpi_x)*dx*Nx
  ymargin=int(my_rank/n_mpi_x)*dy*Ny
  write(fnum,'(A30)')"# T         X          Y         Z"
  DO I=1,Nx
    do J=1,Ny
      write(fnum,'(4f16.6)') T,(real(I)-0.5d0)*Dx+xmargin,(real(J)-0.5d0)*Dy+ymargin,outarray(I,J)
    ENDDO
    write(fnum,*) 
  enddo

  flush(fnum)
  close(fnum)
end subroutine fwrite_txyz
