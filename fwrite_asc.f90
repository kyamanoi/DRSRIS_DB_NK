
subroutine fwrite_asc(name,outarray)
  use variables
  implicit none
  character(*) name
  character(len=1000) fname,timechar
  character(len=13) fmt
  character(len=4) rankname
  ! real(8),dimension(Nx,Ny) :: outarray
  real(8),dimension(0:Nx+1,0:Ny+1) :: outarray
  write(rankname,'(I4.4)') my_rank

  write(fmt,'(''('',I0,''E18.6e3)'')') Nx
  write(timechar,'(F0.4)') T
  if(T==0.d0)timechar="0.0000"
  ! fname= "./output/"//trim(name)//trim(timechar)//".asc"
  fname= "./output/rank"//rankname//"/asc/"//trim(name)//trim(timechar)//".asc"
  OPEN(5001,FILE=fname,STATUS='UNKNOWN')
  write(5001,'(A5,I16)')"ncols",Nx
  write(5001,'(A5,I16)')"nrows",Ny
  write(5001,'(A9,f16.6)')"xllcorner",xllcorner
  write(5001,'(A9,f16.6)')"yllcorner",yllcorner
  write(5001,'(A8,f16.6)')"cellsize",Dx
  ! do J=Ny,1,-1
  do J=1,Ny
    write(5001,fmt) (outarray(I,J),I=1,Nx)
  enddo
  flush(5001)
  close(5001)
end subroutine fwrite_asc


subroutine fwrite_asc_int(name,outarray)
  use variables
  implicit none
  character(*) name
  character(len=1000) fname,timechar
  character(len=13) fmt
  character(len=4) rankname
  ! real(8),dimension(Nx,Ny) :: outarray
  integer,dimension(0:Nx+1,0:Ny+1) :: outarray
  write(rankname,'(I4.4)') my_rank

  write(fmt,'(''('',I0,''I5)'')') Nx
  write(timechar,'(F0.4)') T
  if(T==0.d0)timechar="0.0000"
  ! fname= "./output/"//trim(name)//trim(timechar)//".asc"
  fname= "./output/rank"//rankname//"/asc/"//trim(name)//trim(timechar)//".asc"

  OPEN(5001,FILE=fname,STATUS='UNKNOWN')
  write(5001,'(A5,I16)')"ncols",Nx
  write(5001,'(A5,I16)')"nrows",Ny
  write(5001,'(A9,f16.6)')"xllcorner",xllcorner
  write(5001,'(A9,f16.6)')"yllcorner",yllcorner
  write(5001,'(A8,f16.6)')"cellsize",Dx
  ! do J=Ny,1,-1
  do J=1,Ny
    write(5001,fmt) (outarray(I,J),I=1,Nx)
  enddo
  flush(5001)
  close(5001)
end subroutine fwrite_asc_int


subroutine fwrite_asc_time(fnum,outarray)
  use variables
  implicit none
  character(len=13) fmt
  character(len=4) rankname
  integer,intent(in) :: fnum
  real(8),dimension(0:Nx+1,0:Ny+1) :: outarray
  write(rankname,'(I4.4)') my_rank

  write(fmt,'(''('',I0,''E18.6e3)'')') Nx
  write(fnum,*)
  write(fnum,'(A5 f16.6)')"time",T
  write(fnum,'(A5,I16)')"ncols",Nx
  write(fnum,'(A5,I16)')"nrows",Ny
  write(fnum,'(A9,f16.6)')"xllcorner",xllcorner
  write(fnum,'(A9,f16.6)')"yllcorner",yllcorner
  write(fnum,'(A8,f16.6)')"cellsize",Dx
  do J=1,Ny
    write(fnum,fmt) (outarray(I,J),I=1,Nx)
  enddo
  flush(fnum)
end subroutine fwrite_asc_time

subroutine fwrite_asc_int_time(fnum,outarray)
  use variables
  implicit none
  character(len=13) fmt
  character(len=4) rankname
  integer,intent(in) :: fnum
  integer,dimension(0:Nx+1,0:Ny+1) :: outarray
  write(rankname,'(I4.4)') my_rank

  write(fmt,'(''('',I0,''I5)'')') Nx
  write(fnum,*)
  write(fnum,'(A5 f16.6)')"time",T
  write(fnum,'(A5,I16)')"ncols",Nx
  write(fnum,'(A5,I16)')"nrows",Ny
  write(fnum,'(A9,f16.6)')"xllcorner",xllcorner
  write(fnum,'(A9,f16.6)')"yllcorner",yllcorner
  write(fnum,'(A8,f16.6)')"cellsize",Dx
  do J=1,Ny
    write(fnum,fmt) (outarray(I,J),I=1,Nx)
  enddo
  flush(fnum)
end subroutine fwrite_asc_int_time
