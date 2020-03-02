subroutine dry_or_wet(hin,zbin,uin,vin,iswet,hcr,Nx,Ny,isforward)
  implicit none
  integer I,J
  integer,intent(in) :: Nx,Ny,isforward
  real(8),intent(in) :: hcr
  real(8),dimension(0:Nx+1,0:Ny+1),intent(in) :: zbin
  real(8),dimension(0:Nx+1,0:Ny+1),intent(inout) :: hin,uin,vin
  real(8),dimension(0:Nx+1,0:Ny+1) :: eta
  integer,dimension(0:Nx+1,0:Ny+1),intent(out) :: iswet

  open(3001,file='./output/log.log')

  ! where(hin>hcr)
  !   eta=hin+zbin
  ! elsewhere
  !   eta=-9999.d0
  ! endwhere

  ! where(hin>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  DO I=0,Nx+1
    DO J=0,Ny+1
      if(hin(I,J)>hcr)then
        eta(I,J)=hin(I,J)+zbin(I,J)
        iswet(I,J)=1
      else
        eta(I,J)=-9999
        iswet(I,J)=0
      endif
    enddo
  enddo

  !Liang et al., employing international journal for numerical methods in fluids, 2007
  DO I=1,Nx
    DO J=1,Ny
      if(iswet(I,J)==0)then
        if(iswet(I-1,J)==0.and.iswet(I,J-1)==0.and.iswet(I+1,J)==0.and.iswet(I,J+1)==0)then
          cycle
        endif
        if(zbin(I,J)>max(eta(I-1,J),eta(I,J-1),eta(I+1,J),eta(I,J+1)))then
          cycle
        endif
        if(eta(I+1,J)==max(eta(I-1,J),eta(I,J-1),eta(I+1,J),eta(I,J+1)).and.eta(I+1,J)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I+1,J)>2.d0*Hcr)then
            hin(I+1,J)=hin(I+1,J)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I+1,J)=eta(I+1,J)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type3",I,J
          else
            cycle
          endif
        elseif(eta(I,J+1)==max(eta(I-1,J),eta(I,J-1),eta(I+1,J),eta(I,J+1)).and.eta(I,J+1)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I,J+1)>2.d0*Hcr)then
            hin(I,J+1)=hin(I,J+1)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I,J+1)=eta(I,J+1)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type4",I,J
          else
            cycle
          endif
        elseif(eta(I-1,J)==max(eta(I-1,J),eta(I,J-1),eta(I+1,J),eta(I,J+1)).and.eta(I-1,J)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I-1,J)>2.d0*Hcr)then
            hin(I-1,J)=hin(I-1,J)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I-1,J)=eta(I-1,J)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type1",I,J
          else
            cycle
          endif
        elseif(eta(I,J-1)==max(eta(I-1,J),eta(I,J-1),eta(I+1,J),eta(I,J+1)).and.eta(I,J-1)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I,J-1)>2.d0*Hcr)then
            hin(I,J-1)=hin(I,J-1)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I,J-1)=eta(I,J-1)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type2",I,J
          else
            cycle
          endif
        endif
      endif
      if(hin(I,J)<Hcr.and.iswet(I,J)==1)then
        write(*,*) "wetdry error"
        write(*,*) hin(I,J),iswet(I,J)
        read(*,*)
      endif
    ENDDO
  enddo

  where(iswet==0)
    uin=0.d0
    vin=0.d0
  endwhere

  ! if(isforward==1)then
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I+1,J).and.hin(I+1,J)>hcr.and.uin(I+1,J)-uin(I,J)<0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J+1).and.hin(I,J+1)>hcr.and.vin(I,J+1)-vin(I,J)<0.d0)))then
  !       !   iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! else
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I-1,J).and.hin(I-1,J)>hcr.and.uin(I,J)-uin(I-1,J)>0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J-1).and.hin(I,J-1)>hcr.and.vin(I,J)-vin(I,J-1)>0.d0)))then
  !         ! iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! endif
end subroutine dry_or_wet

subroutine dry_or_wet_f(hin,zbin,uin,vin,iswet,hcr,Nx,Ny)
  implicit none
  integer I,J
  integer,intent(in) :: Nx,Ny
  real(8),intent(in) :: hcr
  real(8),dimension(0:Nx+1,0:Ny+1),intent(in) :: zbin
  real(8),dimension(0:Nx+1,0:Ny+1),intent(inout) :: hin,uin,vin
  real(8),dimension(0:Nx+1,0:Ny+1) :: eta
  integer,dimension(0:Nx+1,0:Ny+1),intent(out) :: iswet

  open(3001,file='./output/log.log')

  ! where(hin>hcr)
  !   eta=hin+zbin
  ! elsewhere
  !   eta=-9999.d0
  ! endwhere

  ! where(hin>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  DO I=0,Nx+1
    DO J=0,Ny+1
      if(hin(I,J)>hcr)then
        eta(I,J)=hin(I,J)+zbin(I,J)
        iswet(I,J)=1
      else
        eta(I,J)=-9999
        iswet(I,J)=0
      endif
    enddo
  enddo

  !Liang et al., employing international journal for numerical methods in fluids, 2007
  DO I=1,Nx
    DO J=1,Ny
      if(iswet(I,J)==0)then
        if(iswet(I+1,J)==0.and.iswet(I,J+1)==0)then
          cycle
        endif
        if(zbin(I,J)>max(eta(I+1,J),eta(I,J+1)))then
          cycle
        endif
        if(eta(I+1,J)==max(eta(I+1,J),eta(I,J+1)).and.eta(I+1,J)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I+1,J)>2.d0*Hcr)then
            hin(I+1,J)=hin(I+1,J)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I+1,J)=eta(I+1,J)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type3",I,J
          else
            cycle
          endif
        elseif(eta(I,J+1)==max(eta(I+1,J),eta(I,J+1)).and.eta(I,J+1)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I,J+1)>2.d0*Hcr)then
            hin(I,J+1)=hin(I,J+1)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I,J+1)=eta(I,J+1)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type4",I,J
          else
            cycle
          endif
        endif
      endif
      if(hin(I,J)<Hcr.and.iswet(I,J)==1)then
        write(*,*) "wetdry error"
        write(*,*) hin(I,J),iswet(I,J)
        read(*,*)
      endif
    ENDDO
  enddo

  where(iswet==0)
    uin=0.d0
    vin=0.d0
  endwhere

  ! if(isforward==1)then
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I+1,J).and.hin(I+1,J)>hcr.and.uin(I+1,J)-uin(I,J)<0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J+1).and.hin(I,J+1)>hcr.and.vin(I,J+1)-vin(I,J)<0.d0)))then
  !       !   iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! else
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I-1,J).and.hin(I-1,J)>hcr.and.uin(I,J)-uin(I-1,J)>0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J-1).and.hin(I,J-1)>hcr.and.vin(I,J)-vin(I,J-1)>0.d0)))then
  !         ! iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! endif
end subroutine dry_or_wet_f

subroutine dry_or_wet_b(hin,zbin,uin,vin,iswet,hcr,Nx,Ny)
  implicit none
  integer I,J
  integer,intent(in) :: Nx,Ny
  real(8),intent(in) :: hcr
  real(8),dimension(0:Nx+1,0:Ny+1),intent(in) :: zbin
  real(8),dimension(0:Nx+1,0:Ny+1),intent(inout) :: hin,uin,vin
  real(8),dimension(0:Nx+1,0:Ny+1) :: eta
  integer,dimension(0:Nx+1,0:Ny+1),intent(out) :: iswet

  open(3001,file='./output/log.log')

  ! where(hin>hcr)
  !   eta=hin+zbin
  ! elsewhere
  !   eta=-9999.d0
  ! endwhere

  ! where(hin>hcr)
  !   iswet=1
  ! elsewhere
  !   iswet=0
  ! endwhere

  DO I=0,Nx+1
    DO J=0,Ny+1
      if(hin(I,J)>hcr)then
        eta(I,J)=hin(I,J)+zbin(I,J)
        iswet(I,J)=1
      else
        eta(I,J)=-9999
        iswet(I,J)=0
      endif
    enddo
  enddo

  !Liang et al., employing international journal for numerical methods in fluids, 2007
  DO I=1,Nx
    DO J=1,Ny
      if(iswet(I,J)==0)then
        if(iswet(I-1,J)==0.and.iswet(I,J-1)==0)then
          cycle
        endif
        if(zbin(I,J)>max(eta(I-1,J),eta(I,J-1)))then
          cycle
        endif
        if(eta(I-1,J)==max(eta(I-1,J),eta(I,J-1)).and.eta(I-1,J)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I-1,J)>2.d0*Hcr)then
            hin(I-1,J)=hin(I-1,J)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I-1,J)=eta(I-1,J)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type1",I,J
          else
            cycle
          endif
        elseif(eta(I,J-1)==max(eta(I-1,J),eta(I,J-1)).and.eta(I,J-1)>zbin(I,J)+2.d0*Hcr)then
          if(hin(I,J-1)>2.d0*Hcr)then
            hin(I,J-1)=hin(I,J-1)-Hcr
            hin(I,J)=hin(I,J)+Hcr
            eta(I,J-1)=eta(I,J-1)-Hcr
            iswet(I,J)=1
            write(3001,*) "modif_type2",I,J
          else
            cycle
          endif
        endif
      endif
      if(hin(I,J)<Hcr.and.iswet(I,J)==1)then
        write(*,*) "wetdry error"
        write(*,*) hin(I,J),iswet(I,J)
        read(*,*)
      endif
    ENDDO
  enddo

  where(iswet==0)
    uin=0.d0
    vin=0.d0
  endwhere

  ! if(isforward==1)then
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I+1,J).and.hin(I+1,J)>hcr.and.uin(I+1,J)-uin(I,J)<0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J+1).and.hin(I,J+1)>hcr.and.vin(I,J+1)-vin(I,J)<0.d0)))then
  !       !   iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! else
  !   DO I=1,Nx
  !     DO J=1,Ny
  !       if(hin(I,J)>hcr)then
  !         iswet(I,J)=1
  !       ! elseif(hin(I,J)<Hcr&
  !       ! &.and.((eta(I,J)<eta(I-1,J).and.hin(I-1,J)>hcr.and.uin(I,J)-uin(I-1,J)>0.d0)&
  !       ! &.and.(eta(I,J)<eta(I,J-1).and.hin(I,J-1)>hcr.and.vin(I,J)-vin(I,J-1)>0.d0)))then
  !         ! iswet(I,J)=1
  !       else
  !         iswet(I,J)=0
  !       endif
  !     ENDDO
  !   ENDDO
  ! endif
end subroutine dry_or_wet_b
