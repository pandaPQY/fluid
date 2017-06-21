module pencil_fft
!use parameters
  implicit none
  save
  integer,parameter :: n=512
  integer(4),parameter :: NULL=0
  integer(8) planx,plany,planz,iplanx,iplany,iplanz

  ! fft arrays
  real        r3(n,n,n),r0(n,n)
  complex     c3(n/2,n,n)
  real        rxyz(n*1+2  ,n,n)
  complex     cxyz(n*1/2+1,n,n)
  complex     cyyyxz(n,1,1,n/2+1,n)
  complex     cyyxz(n,     1,n/2+1,n)
  complex     czzzxy(n,1,1,n/2+1,n)

  ! communication coarrays
  complex ctransfer1(n/2,n,1)
  complex ctransfer2(n,n/2+1,1)
  complex ctransfer3(n,n,1,1)
  complex ctransfer4(n/2+1,n,1)

  ! equivalence statements
  equivalence(cyyyxz,cyyxz,r3,c3)
  equivalence(czzzxy,rxyz,cxyz)

  contains

  subroutine pencil_fft_forward
    implicit none
    save
    sync all
    call c2x
    call sfftw_execute(planx)
    call x2y
    call sfftw_execute(plany)
    call y2z
    call sfftw_execute(planz)
    call z2y
    call y2x
  endsubroutine

  subroutine pencil_fft_backward
    implicit none
    save
    sync all
    call x2y
    call y2z
    call sfftw_execute(iplanz)
    call z2y
    call sfftw_execute(iplany)
    call y2x
    call sfftw_execute(iplanx)
    call x2c
    r3=r3/(n*1)/(n*1)/(n*1)
  endsubroutine

  subroutine c2x
    implicit none
    save
    integer(8) i0,i1,i2,islab
!    do islab=1,n ! loop over cells in z, extract slabs
!      ctransfer1=c3(:,:,islab) ! 1 slabs of c3 copied to ctransfer1
        ! i1=mod()
        cxyz(1:n/2,:,:)=c3!ctransfer1
!    enddo
  endsubroutine

  subroutine x2y
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,n ! loop over z
!        ctransfer2=transpose(cxyz(:,:,islab))
        cyyxz(:,1,:,islab)=transpose(cxyz(:,:,islab))!ctransfer2
    enddo
  endsubroutine

  subroutine y2z
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,n/2+1 ! loop over slices in x direction
!        ctransfer3=transpose(cyyyxz(:,:,:,islab,:))
        czzzxy(:,1,1,islab,:)=transpose(cyyyxz(:,1,1,islab,:))!ctransfer3
    enddo
  endsubroutine

  subroutine z2y
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,n/2+1 ! loop over slices in x direction
!        ctransfer3=transpose(czzzxy(:,:,:,islab,:))
        cyyyxz(:,1,1,islab,:)=transpose(czzzxy(:,1,1,islab,:))!ctransfer3
    enddo
  endsubroutine

  subroutine y2x
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,n ! loop over z
 !       ctransfer4=transpose(cyyxz(:,:,:,islab))
        cxyz(:,:,islab)=transpose(cyyxz(:,1,:,islab))!ctransfer4
    enddo
  endsubroutine

  subroutine x2c
    implicit none
    save
    integer(8) i0,i1,i2,islab
!    do islab=1,n
!        ctransfer1=cxyz(:,:,islab)
        c3=cxyz(1:n/2,:,:)!ctransfer1
!    enddo
  endsubroutine

  subroutine create_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_plan_many_dft_r2c(planx,1,n*1,n*n,cxyz,NULL,1,n*1+2,cxyz,NULL,1,n*1/2+1,FFTW_MEASURE)
    call sfftw_plan_many_dft_c2r(iplanx,1,n*1,n*n,cxyz,NULL,1,n*1/2+1,cxyz,NULL,1,n*1+2,FFTW_MEASURE)
    call sfftw_plan_many_dft(plany,1,n*1,(n/2+1)*n,cyyxz,NULL,1,n*1,cyyxz,NULL,1,n*1,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplany,1,n*1,(n/2+1)*n,cyyxz,NULL,1,n*1,cyyxz,NULL,1,n*1,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(planz,1,n*1,(n/2+1)*n,czzzxy,NULL,1,n*1,czzzxy,NULL,1,n*1,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplanz,1,n*1,(n/2+1)*n,czzzxy,NULL,1,n*1,czzzxy,NULL,1,n*1,FFTW_BACKWARD,FFTW_MEASURE)
  endsubroutine

  subroutine destroy_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_destroy_plan(planx)
    call sfftw_destroy_plan(iplanx)
    call sfftw_destroy_plan(plany)
    call sfftw_destroy_plan(iplany)
    call sfftw_destroy_plan(planz)
    call sfftw_destroy_plan(iplanz)
  endsubroutine
endmodule
