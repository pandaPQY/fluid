module calculation
contains
subroutine getuc4(u,uc4,nx,ny,nz)
implicit none
integer nx,ny,nz
real u(4,nx,ny,nz)
real(8) u8(4,nx,ny,nz)
real(8), dimension(nx,ny,nz) :: rho,rhog2,v2 !x
real(8), dimension(4,nx,ny,nz) :: u4,uc8 !x
real uc4(4,nx,ny,nz)
u4=0
u8=u
v2=3*sum(u8(2:4,:,:,:)**2,1)/4 !x
rho=2*sqrt(u8(1,:,:,:)**2-v2)-u8(1,:,:,:)    ! u=T^00, v=T^0k !x
rhog2=rho/4+3*u8(1,:,:,:)/4              !    v=u(2,:)/rhog2 !x
u4(2:4,:,:,:) = 1.5*u8(2:4,:,:,:)*spread(rho**0.25/sqrt(rhog2),1,3)!u4=n^i !x
uc8=u4/spread(sqrt(rho),1,4)             !uc4=n^i/average(n)
uc8(1:3,:,:,:)=uc8(2:4,:,:,:)
uc8(4,:,:,:)=sqrt(sum(uc8(1:3,:,:,:)**2,1))
uc4=uc8
end subroutine getuc4

subroutine curlv(uc4,curl,nx,ny,nz)
implicit none
integer nx,ny,nz
real,dimension(4,nx,ny,nz) :: uc4,curl
integer i,j,k,ip,jp,kp
real vyx,vxy,vzy,vyz,vzx,vxz,curlx,curly,curlz

!$omp parallel do private(i,j,k,ip,jp,kp,vzx,vyx,vxy,vzy,vxz,vyz,curlx,curly,curlz) shared(uc4,curl)
  do k=1,nz !z                   
     do j=1,ny !y
        do i=1,nx
           kp=mod(k,nz)+1
           jp=mod(j,ny)+1
           ip=mod(i,nx)+1 !y
       vzx=(uc4(3,ip,j,k)-uc4(3,i,j,k)+uc4(3,ip,jp,k)-uc4(3,i,jp,k) + &
           uc4(3,ip,j,kp)-uc4(3,i,j,kp)+uc4(3,ip,jp,kp)-uc4(3,i,jp,kp) )/4
       vyx=(uc4(2,ip,j,k)-uc4(2,i,j,k)+uc4(2,ip,jp,k)-uc4(2,i,jp,k) + &
           uc4(2,ip,j,kp)-uc4(2,i,j,kp)+uc4(2,ip,jp,kp)-uc4(2,i,jp,kp) )/4
       vxy=(uc4(1,i,jp,k)-uc4(1,i,j,k)+uc4(1,ip,jp,k)-uc4(1,ip,j,k) + &
           uc4(1,i,jp,kp)-uc4(1,i,j,kp)+uc4(1,ip,jp,kp)-uc4(1,ip,j,kp) )/4
       vzy=(uc4(3,i,jp,k)-uc4(3,i,j,k)+uc4(3,ip,jp,k)-uc4(3,ip,j,k) + &
           uc4(3,i,jp,kp)-uc4(3,i,j,kp)+uc4(3,ip,jp,kp)-uc4(3,ip,j,kp) )/4
       vxz=(uc4(1,i,j,kp)-uc4(1,i,j,k)+uc4(1,ip,j,kp)-uc4(1,ip,j,k) + &
           uc4(1,i,jp,kp)-uc4(1,i,jp,k)+uc4(1,ip,jp,kp)-uc4(1,ip,jp,k) )/4
       vyz=(uc4(2,i,j,kp)-uc4(2,i,j,k)+uc4(2,ip,j,kp)-uc4(2,ip,j,k) + &
           uc4(2,i,jp,kp)-uc4(2,i,jp,k)+uc4(2,ip,jp,kp)-uc4(2,ip,jp,k) )/4
       curlx=vzy-vyz
       curly=vxz-vzx
       curlz=vyx-vxy
       curl(1,i,j,k)=curlx
       curl(2,i,j,k)=curly
       curl(3,i,j,k)=curlz
       curl(4,i,j,k)=curlx**2+curly**2+curlz**2
       enddo
    enddo
  enddo
end subroutine curlv

subroutine findcurl(curl,nx,ny,nz)
integer nx,ny,nz
real curl(4,nx,ny,nz)
write(516) maxval(curl(1,:,:,:)),maxloc(curl(1,:,:,:)),minval(curl(1,:,:,:)),minloc(curl(1,:,:,:))
write(516) maxval(curl(2,:,:,:)),maxloc(curl(2,:,:,:)),minval(curl(2,:,:,:)),minloc(curl(2,:,:,:))
write(516) maxval(curl(3,:,:,:)),maxloc(curl(3,:,:,:)),minval(curl(3,:,:,:)),minloc(curl(3,:,:,:))
write(516) maxval(curl(4,:,:,:)),maxloc(curl(4,:,:,:)),minval(curl(4,:,:,:)),minloc(curl(4,:,:,:))
print*,'exmax',maxval(curl(1,:,:,:)),maxval(curl(2,:,:,:)),maxval(curl(3,:,:,:)),maxval(curl(4,:,:,:))
print*,'exmin',minval(curl(1,:,:,:)),minval(curl(2,:,:,:)),minval(curl(3,:,:,:)),minval(curl(4,:,:,:))
end subroutine findcurl

subroutine smooth(u,usm,nx,ny,nz)
integer nx,ny,nz
real,dimension(4,nx,ny,nz) :: u,usm,tmp,tmp2,tmp3
double precision utmp(4,0:nx+2,0:ny+2,0:nz+2),tmp4(3)
double precision,dimension(4,nx,ny,nz) :: usm8
integer i,j,k

utmp(:,1:nx,1:ny,1:nz)=u
utmp(:,0,0,0)=u(:,nx,ny,nz)
utmp(:,0,1:ny,0)=u(:,nx,:,nz)
utmp(:,1:nx,0,0)=u(:,:,ny,nz)
utmp(:,0,0,1:nz)=u(:,nx,ny,:)
utmp(:,0,1:ny,1:nz)=u(:,nx,:,:)
utmp(:,1:nx,0,1:nz)=u(:,:,ny,:)
utmp(:,1:nx,1:ny,0)=u(:,:,:,nz)
utmp(:,nx+1:nx+2,ny+1:ny+2,nz+1:nz+2)=u(:,1:2,1:2,1:2)
utmp(:,nx+1:nx+2,1:ny,1:nz)=u(:,1:2,:,:)
utmp(:,1:nx,ny+1:ny+2,1:nz)=u(:,:,1:2,:)
utmp(:,1:nx,1:ny,nz+1:nz+2)=u(:,:,:,1:2)
utmp(:,nx+1:nx+2,ny+1:ny+2,1:nz)=u(:,1:2,1:2,:)
utmp(:,nx+1:nx+2,1:ny,nz+1:nz+2)=u(:,1:2,:,1:2)
utmp(:,1:nx,ny+1:ny+2,nz+1:nz+2)=u(:,:,1:2,1:2)
utmp(:,0,ny+1:ny+2,nz+1:nz+2)=u(:,nx,1:2,1:2)
utmp(:,nx+1:nx+2,0,nz+1:nz+2)=u(:,1:2,ny,1:2)
utmp(:,nx+1:nx+2,ny+1:ny+2,0)=u(:,1:2,1:2,nz)
utmp(:,0,ny+1:ny+2,1:nz)=u(:,nx,1:2,:)
utmp(:,nx+1:nx+2,0,1:nz)=u(:,1:2,ny,:)
utmp(:,0,1:ny,nz+1:nz+2)=u(:,nx,:,1:2)
utmp(:,nx+1:nx+2,1:ny,0)=u(:,1:2,:,nz)
utmp(:,1:nx,0,nz+1:nz+2)=u(:,:,ny,1:2)
utmp(:,1:nx,ny+1:ny+2,0)=u(:,:,1:2,nz)
utmp(:,0,0,nz+1:nz+2)=u(:,nx,ny,1:2)
utmp(:,0,ny+1:ny+2,0)=u(:,nx,1:2,nz)
utmp(:,nx+1:nx+2,0,0)=u(:,1:2,ny,nz)
usm=0 
usm8=0
!$omp parallel do private(i,j,k,tmp4) shared(utmp,usm)
do k=1,nz
   do j=1,ny
      do i=1,nx
         tmp4(1)=sum(utmp(1,i-1:i+2,j-1:j+2,k-1:k+2))
         tmp4(2)=sum(utmp(2,i-1:i+2,j-1:j+2,k-1:k+2))
         tmp4(3)=sum(utmp(3,i-1:i+2,j-1:j+2,k-1:k+2))
         usm(1:3,i,j,k)=real(tmp4,4)
      enddo
   enddo
enddo
!!!$omp parallel do collapse(3) private(i,j,k,tmp,tmp2,tmp3) reduction(+:usm8)
!do k=1,4
!   do j=1,4
!      do i=1,4
!         tmp=cshift(u,3-i,2)
!         tmp2=cshift(tmp,3-j,3)
!         tmp3=cshift(tmp2,3-k,4)
!         usm8=usm8+real(tmp3,8)
!      enddo
!   enddo
!enddo
usm=usm/4**3
print*,'smooth',maxval(usm(1:3,:,:,:)),minval(usm(1:3,:,:,:)),maxval(sum(usm(1:3,:,:,:),1)),minval(sum(usm(1:3,:,:,:),1))
end subroutine smooth

subroutine findcurlsm(curl3,nx,ny,nz)
integer nx,ny,nz
real curl3(4,nx,ny,nz)
!!$omp end parallel do
curl3(4,:,:,:)=sum(curl3(1:3,:,:,:)**2,1)
write(517) maxval(curl3(1,:,:,:)),maxloc(curl3(1,:,:,:)),minval(curl3(1,:,:,:)),minloc(curl3(1,:,:,:))
write(517) maxval(curl3(2,:,:,:)),maxloc(curl3(2,:,:,:)),minval(curl3(2,:,:,:)),minloc(curl3(2,:,:,:))
write(517) maxval(curl3(3,:,:,:)),maxloc(curl3(3,:,:,:)),minval(curl3(3,:,:,:)),minloc(curl3(3,:,:,:))
write(517) maxval(curl3(4,:,:,:)),maxloc(curl3(4,:,:,:)),minval(curl3(4,:,:,:)),minloc(curl3(4,:,:,:))
print*,'smmax',maxval(curl3(1,:,:,:)),maxval(curl3(2,:,:,:)),maxval(curl3(3,:,:,:)),maxval(curl3(4,:,:,:))
print*,'smmin',minval(curl3(1,:,:,:)),minval(curl3(2,:,:,:)),minval(curl3(3,:,:,:)),minval(curl3(4,:,:,:))
end subroutine findcurlsm

subroutine pdf_curl(curl,nx,ny,nz,nbin)
integer nx,ny,nz,nbin
real curl(4,nx,ny,nz)
real(8) pdf(4,2,nbin)
integer i
do i=1,4
   call calcpdf(curl(i,:,:,:),pdf,nx,ny,nz,nbin)
enddo
write(519) pdf
print*,'pdf',sum(pdf(4,2,:))
end subroutine pdf_curl

subroutine pdf_curl_sm(usm,nx,ny,nz,nbin)
integer nx,ny,nz,nbin
real usm(4,nx,ny,nz)
real(8) pdfsm(4,2,nbin)
integer i
do i=1,4
   call calcpdf(usm(i,:,:,:),pdfsm,nx,ny,nz,nbin)
enddo
write(518) pdfsm
print*,'pdf',sum(pdfsm(4,2,:))
end subroutine pdf_curl_sm

subroutine writepdfv(uc4,ue,nx,ny,nz,check1,check2)
use pencil_fft
implicit none
integer nx,ny,nz
!integer nn(3)
real,dimension(4,nx,ny,nz) :: uc4
complex,dimension(3,nx/2+1,ny,nz) :: ctmp,ctmpe,ctmpee
integer i,j,k,ii,jj,kk,i_dim
real,dimension(4,nx,ny,nz) :: ue,ub,uee,ubb
integer,parameter :: nbin=40
real(8),dimension(4,2,nbin) :: pdfv,pdfve,pdfvb
logical check1,check2
!nn=(/nx,ny,nz/)
!print*,'before fourn'
!call create_penfft_plan
!print*,'maxmin=',maxval(uc4(1:3,:,:,:)),minval(uc4(1:3,:,:,:))
do i=1,3
  r3=uc4(i,:,:,:)
  call pencil_fft_forward
  ctmp(i,:,:,:)=cxyz
enddo
print*,'fourn done'
!!$omp parallel do default(none) shared(ctmp,ctmpe,nx,ny,nz) private(k,kk,j,jj,i,ii,i_dim) 
do k=1,nz
   kk=mod(k+nz/2-1,nz)-nz/2
   do j=1,ny
      jj=mod(j+ny/2-1,ny)-ny/2
      do i=1,nx/2+1
         ii=mod(i+nx/2-1,nx)-nx/2
         do i_dim=1,3
            call del(ctmp(1:3,i,j,k),ctmpe(i_dim,i,j,k),nx,ii,jj,kk,i_dim)
         enddo
      end do
   end do
end do
ctmpe(:,1,1,1)=0
print*,'del done'
do i=1,3
   cxyz=ctmpe(i,:,:,:)
   call pencil_fft_backward
   ue(i,:,:,:)=r3
enddo
!call destroy_penfft_plan
print*,'fourn back'
ue(1:3,:,:,:)=real(ue(1:3,:,:,:))
print*,'real'
ue(4,:,:,:)=sqrt(ue(1,:,:,:)**2+ue(2,:,:,:)**2+ue(3,:,:,:)**2)
print*,'ue4'
!print*,'maxmin e=',maxval(ue(1:3,:,:,:)),minval(ue(1:3,:,:,:))
ub(1:3,:,:,:)=real(uc4(1:3,:,:,:),8)-real(ue(1:3,:,:,:),8)
do i=1,3
print*,'maxmin  =',i,maxval(uc4(i,:,:,:)),minval(uc4(i,:,:,:))
print*,'maxmin e=',i,maxval(ue(i,:,:,:)),minval(ue(i,:,:,:))
print*,'maxmin b=',i,maxval(ub(i,:,:,:)),minval(ub(i,:,:,:))
enddo
ub(4,:,:,:)=sqrt(ub(1,:,:,:)**2+ub(2,:,:,:)**2+ub(3,:,:,:)**2)
if (check2==1) then
print*,'output ebmode'
!write(550) uc4(:,nx/2-10:nx/2+10,ny/2-10:ny/2+10,nz/2-10:nz/2+10)
!write(551) ue(:,nx/2-10:nx/2+10,ny/2-10:ny/2+10,nz/2-10:nz/2+10)
!write(552) ub(:,nx/2-10:nx/2+10,ny/2-10:ny/2+10,nz/2-10:nz/2+10)
endif
end subroutine writepdfv

subroutine del(v,v2,n,k1,k2,k3,i_dim)
implicit none
integer n
complex v(3),vv(3),v2
integer k1,k2,k3,i_dim,dim1,dim2,dim3,k(3)
complex ekx(3),ekx2(3)
real(8) delta12,delta13
real(8) laplas(3),delta(3)
real(8),parameter :: pi=4*atan(1.)
dim1=i_dim
dim2=mod(dim1+1,3)+1
dim3=mod(dim2+1,3)+1
k(dim1)=k1
k(dim2)=k2
k(dim3)=k3
vv(dim1)=v(1)
vv(dim2)=v(2)
vv(dim3)=v(3)
laplas=2.*cos(2*pi*k/n)-2.
delta=-1.*sin(2*pi*k(1)/n)*sin(2*pi*k/n)
delta(1)=laplas(1)

if (sum(laplas)==0.) then
   v2=0
else
   v2=sum(delta*vv)/sum(laplas)
endif
end subroutine del

subroutine calcpdf(u,pdf,nx,ny,nz,nbin)
implicit none
integer nx,ny,nz
real u(nx,ny,nz)
integer nbin
real(8) pdf(2,nbin)
real maxu,minu,interval
integer i,j,k,ibin
pdf=0
minu=minval(u)
maxu=maxval(u)
interval=(maxu-minu)/nbin
!!$omp parallel do private(ibin) shared(nbin,pdf,minu,interval)
do ibin=1,nbin
   pdf(1,ibin)=minu+interval*(ibin-1)
enddo
!print*,'prepared'
!!$omp parallel do private(i,j,k,ibin) shared(u,minu,interval,nx,ny,nz) reduction(+:pdf)
do k=1,nz
   do j=1,ny
      do i=1,nx
         ibin=real(u(i,j,k)-minu,4)/interval+1
         pdf(2,ibin)=pdf(2,ibin)+1
      enddo
   enddo
enddo
print*,'pdf',sum(pdf(2,:))
end subroutine calcpdf

subroutine shift(u,us,ns,nx,ny,nz)
implicit none
integer nx,ny,nz,ns
real,dimension(nx,ny,nz) :: u,us
integer i,j,k
!$omp parallel do private(j,k) shared(nz,ny,us,u,ns)
do k=1,nz
   do j=1,ny
      us(:,j,k)=cshift(u(:,j,k),ny+ny/2-j-k+1+ns,1)
   enddo
enddo
!!$omp end parallel do
!us=cshift(us,ns/2,2)
!us=cshift(us,ns/2,3)
end subroutine shift

subroutine  writeshift(u,uc4,ns,nx,ny,nz,check1,check2)
implicit none
integer nx,ny,nz,ns 
real,dimension(4,nx,ny,nz) :: u,uc4,curl,curlsm
real,dimension(nx,ny,nz) :: us,us1,us2
real,dimension(ny,nz) :: output
integer i 
logical check1,check2
call shift(u(1,:,:,:),us,ns,nx,ny,nz)
call curlv(uc4,curl,nx,ny,nz)
curl(4,:,:,:)=sum(curl(1:3,:,:,:),1)
call shift(curl(4,:,:,:),us1,ns,nx,ny,nz)
if (check1==1) then
   write(560) us(1,:,:)
   output=1./sqrt(3.)*us1(1,:,:)
   write(561) sign(sqrt(abs(output)),output)
endif
if (check2==1) then
   call smooth(curl,curlsm,nx,ny,nz)
   curlsm(4,:,:,:)=sum(curlsm(1:3,:,:,:),1)
   call shift(curlsm(4,:,:,:),us2,ns,nx,ny,nz)
   output=1./sqrt(3.)*us2(1,:,:)
   print*,'output',maxval(output),minval(output)
   write(562) sign(sqrt(abs(output)),output)
endif
end subroutine writeshift

end module calculation
