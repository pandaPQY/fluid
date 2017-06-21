module parameters
!use pencil_fft
  integer, parameter :: ng=512,nx=ng,ny=ng,nz=ng,neq=4
  integer, parameter :: nbin=40
  integer, parameter :: ns=-150  !!!!! shifted cells
  real, dimension(neq,nx,ny,nz) :: u
  character(4) npd
  character(30) fn
  character(4) bt
  integer, parameter :: nt_last=0
  integer, parameter :: nt_this=nt_last+150           !!!!!!!!!!!!!!!!!!!!!!!change
  logical, parameter :: backup =0
  integer, parameter :: id=256  !xyz
  character(1),parameter :: ax='y' !xyz
  character(4) face
  character(5),parameter :: sam='_t2'              !!!!!!!!!!!!!!!!!!!!!!change
  character(5),parameter :: sam2='t2'              !!!!!!!!!!!!!!!!!!!!!!change
  character(10),parameter :: sam3='_t2.dat'              !!!!!!!!!!!!!!!!!!!!!!change
  character(2) tail
  real,dimension(neq,nx,ny,nz) :: curl,curlsm,uc4,ue
  logical,parameter :: output_u=0
  logical,parameter :: output_curl=0
  logical,parameter :: output_pdfcurl=0
  logical,parameter :: output_pdfcurlsm=0
  logical,parameter :: output_findcurl=0
  logical,parameter :: output_findcurlsm=0
  logical,parameter :: output_pdf_ebmode=0
  logical,parameter :: output_ebmode=0
  logical,parameter :: output_backup=0
  logical,parameter :: output_shift=1
  logical,parameter :: output_shiftsm=0
contains 
subroutine openfile
  if (nt_last<50) tail=''
  if (nt_last==50) tail='_2'
  if (nt_last==100) tail='_3'
  write(npd,'(i4)') ng
  write(face,'(i4)') id
  print*,'write in:'
!   fn='ucrl'//ax//'_'//trim(adjustl(npd))//'_'//trim(adjustl(face))//ax//'_t10.dat'
!   print*,fn
!   open(553,file=trim(fn),status='replace',access='stream')

if (output_curl==1) then
   fn='ucrl'//ax//'_'//trim(adjustl(npd))//'_'//trim(adjustl(face))//ax//trim(sam3)
   print*,fn
   open(521,file=trim(fn),status='replace',access='stream')
endif
if (output_u==1) then
   fn='u0_'//trim(adjustl(npd))//'_'//trim(adjustl(face))//ax//trim(sam3)
   print*,fn
   open(522,file=trim(fn),status='replace',access='stream')
endif
if (output_pdfcurl==1) then
   fn='pdfcurl_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(519,file=trim(fn),status='replace',access='stream')
endif
if (output_pdfcurlsm==1) then
   fn='pdfcurl_'//trim(adjustl(npd))//'_sm'//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(518,file=trim(fn),status='replace',access='stream')
endif
if (output_findcurl==1) then
   fn='findcurl'//'_'//trim(adjustl(npd))//'_ex'//trim(sam2)//trim(tail)//'.dat'
   print*,fn
   open(516,file=trim(fn),status='replace',access='stream')
endif
if (output_findcurlsm==1) then
   fn='findcurl'//'_'//trim(adjustl(npd))//'_sm'//trim(sam2)//trim(tail)//'.dat'
   print*,fn
   open(517,file=trim(fn),status='replace',access='stream')
endif
if (output_pdf_ebmode==1) then
   fn='pdfv'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(513,file=trim(fn),status='replace',access='stream')
   fn='pdfve'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(514,file=trim(fn),status='replace',access='stream')
   fn='pdfvb'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(515,file=trim(fn),status='replace',access='stream')
endif
if (output_ebmode==1) then
   fn='v'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(550,file=trim(fn),status='replace',access='stream')
   fn='ve'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(551,file=trim(fn),status='replace',access='stream')
   fn='vb'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(552,file=trim(fn),status='replace',access='stream')
endif
if (output_shift==1) then
   fn='u0shift-150'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(560,file=trim(fn),status='replace',access='stream')
   fn='ucrlshift-150'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(561,file=trim(fn),status='replace',access='stream')
endif
if (output_shiftsm==1) then
   fn='ucrlshift20sm'//'_'//trim(adjustl(npd))//trim(sam)//trim(tail)//'.dat'
   print*,fn
   open(562,file=trim(fn),status='replace',access='stream')
endif
  if (backup == 1) then
     print*, 'backup, from',nt_last,' to',nt_this
     write(bt,'(i4)') nt_last
     fn='backup_u'//trim(adjustl(npd))//'_'//trim(adjustl(bt))//trim(sam3)
     print*,'last time:',fn
     open(527,file=trim(fn),status='old',access='stream')
     read(527) u
     close(527)
  else
     call init(u,nx,ny,nz)
     print*,'from the beginnig to',nt_this
  endif
  write(bt,'(i4)') nt_this
if (output_backup==1) then
   fn='backup_u'//trim(adjustl(npd))//'_'//trim(adjustl(bt))//trim(sam3)
   print*,'this time:',fn
   open(520,file=fn,status='replace',access='stream')
endif
end subroutine openfile

  subroutine initrandom
    implicit none
    complex ctmp(nx,ny,nz)
    real rtmp(2,nx/2,ny,nz)
    integer nn(3)
    integer i,j,k,ii,jj,kk,k2
    real ak,var,pi
    nn=(/nx,ny,nz /)
    write(*,*)'in the initrandom'
    pi=4*atan(1.)
    call random_seed()
    call random_number(rtmp)
    ctmp(::2,:,:)=sqrt(-2.*log(rtmp(1,:,:,:)))*cos(2*pi*rtmp(2,:,:,:))
    ctmp(2::2,:,:)=sqrt(-2.*log(rtmp(1,:,:,:)))*sin(2*pi*rtmp(2,:,:,:))
    call fourn(ctmp,nn,3,1)
!$omp parallel do default(none) shared(ctmp) private(k,kk,j,jj,i,ii,k2,ak)
    do k=1,nz
       kk=k-1
       if (k>nz/2)kk=kk-nz
       do j=1,ny
          jj=j-1
          if (j>ny/2)jj=jj-ny
          do i=1,nx
             ii=i-1
             if (ii>nx/2) ii=ii-nx
             k2=ii**2+jj**2+kk**2
             ak=sqrt(k2*1.)/max(nx,ny,nz)
             if (k2 .eq. 0) then
                ctmp(i,j,k)=0
                cycle
             end if
             ctmp(i,j,k)=ctmp(i,j,k)/ak**1.5
             if (ak>0.2) ctmp(i,j,k)=0
!             if (ak>0.2) print*,'ak = ',ak
          end do
       end do
    end do
    call fourn(ctmp,nn,3,-1)
    u(1,:,:,:)=real(ctmp)
    u(2:,:,:,:)=0
    var=sum(u(1,:,:,:)**2/(nx*ny*nz))

!!!!!!!!!!!perturbation
    u(1,:,:,:)=1+u(1,:,:,:)/sqrt(var)/10   !/20
    !u(1,:,:,:)=1+u(1,:,:,:)/sqrt(var)/40    

    write(*,*) 'initial condition: maxu = ',maxval(u(1,:,:,:)),' minu = ',minval(u(1,:,:,:))
  end subroutine initrandom

  subroutine init(u,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(4,nx,ny,nz)
  real(8) var,tmp1,tmp2,ave
  integer,parameter :: until=100
!    pi=4*atan(1.)
    u=0
    u(1,:,:,:)=1.
    u(1,1:until,:,:)   =u(1,1:until,:,:)    +0.1
    u(1,until+1:nx,:,:)=u(1,until+1:nx,:,:) -0.1
    u(1,:,1:until,:)   =u(1,:,1:until,:)    +0.2
    u(1,:,until+1:ny,:)=u(1,:,until+1:ny,:) -0.2
    u(1,:,:,1:until)   =u(1,:,:,1:until)    +0.3
    u(1,:,:,until+1:nz)=u(1,:,:,until+1:nz) -0.3
    u(2:4,:,:,:)=0
!    u(1,:,:,:)=u(1,:,:,:)/sqrt(3.)
    tmp1=(sum(real(u(1,:,:,:),8))/nx/ny/nz)**2
    tmp2=sum(real(u(1,:,:,:),8)**2)/nx/ny/nz
    ave=sum(real(u(1,:,:,:),8))/nx/ny/nz
    var=sqrt(tmp2-tmp1)
    print*,'test: var=',var,tmp2,tmp1,ave
!    u(1,:,:,:)=1.+(u(1,:,:,:)-1.)/var*0.1
    tmp1=(sum(real(u(1,:,:,:),8))/nx/ny/nz)**2
    tmp2=sum(real(u(1,:,:,:),8)**2)/nx/ny/nz
    ave=sum(real(u(1,:,:,:),8))/nx/ny/nz
    var=sqrt(tmp2-tmp1)
    print*,'test after trans: var=',var,tmp2,tmp1,ave
  end subroutine init

end module parameters
