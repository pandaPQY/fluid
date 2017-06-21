module writefile
contains
subroutine writecurlx(uc4,ny,nz)
integer ny,nz
real, dimension(4,2,ny,nz) :: uc4
integer i,j,k,ip,jp,kp
real vzy,vzx,curlx
i=1
ip=i+1
do k=1,nz
   kp=mod(k,nz)+1
   do j=1,ny
      jp=mod(j,ny)+1
     vzy=(uc4(3,i,jp,k)-uc4(3,i,j,k)+uc4(3,ip,jp,k)-uc4(3,ip,j,k) + &
         uc4(3,i,jp,kp)-uc4(3,i,j,kp)+uc4(3,ip,jp,kp)-uc4(3,ip,j,kp) )/4
     vyz=(uc4(2,i,j,kp)-uc4(2,i,j,k)+uc4(2,ip,j,kp)-uc4(2,ip,j,k) + &
         uc4(2,i,jp,kp)-uc4(2,i,jp,k)+uc4(2,ip,jp,kp)-uc4(2,ip,jp,k) )/4
     curlx=vzy-vyz
     write(521) sign(sqrt(abs(curlx)),curlx)
  enddo
enddo
end subroutine writecurlx

subroutine writecurly(uc4,nx,nz)
integer nx,nz
real, dimension(4,nx,2,nz) :: uc4    !y
integer i,j,k,ip,jp,kp 
real vxz,vzx,curly
j=1
jp=j+1
do k=1,nz
   kp=mod(k,nz)+1
      do i=1,nx
         ip=mod(i,nx)+1
     vzx=(uc4(3,ip,j,k)-uc4(3,i,j,k)+uc4(3,ip,jp,k)-uc4(3,i,jp,k) + &
         uc4(3,ip,j,kp)-uc4(3,i,j,kp)+uc4(3,ip,jp,kp)-uc4(3,i,jp,kp) )/4
     vxz=(uc4(1,i,j,kp)-uc4(1,i,j,k)+uc4(1,ip,j,kp)-uc4(1,ip,j,k) + &
         uc4(1,i,jp,kp)-uc4(1,i,jp,k)+uc4(1,ip,jp,kp)-uc4(1,ip,jp,k) )/4
     curly=vxz-vzx
     write(521) sign(sqrt(abs(curly)),curly)
  enddo
enddo
end subroutine writecurly

subroutine writecurly2(uc4,nx,nz)
integer nx,nz
real, dimension(4,nx,2,nz) :: uc4    !y
integer i,j,k,ip,jp,kp
real vxz,vzx,curly
j=1
jp=j+1
do k=1,nz
   kp=mod(k,nz)+1
      do i=1,nx
         ip=mod(i,nx)+1
     vzx=(uc4(3,ip,j,k)-uc4(3,i,j,k)+uc4(3,ip,jp,k)-uc4(3,i,jp,k) + &
         uc4(3,ip,j,kp)-uc4(3,i,j,kp)+uc4(3,ip,jp,kp)-uc4(3,i,jp,kp) )/4
     vxz=(uc4(1,i,j,kp)-uc4(1,i,j,k)+uc4(1,ip,j,kp)-uc4(1,ip,j,k) + &
         uc4(1,i,jp,kp)-uc4(1,i,jp,k)+uc4(1,ip,jp,kp)-uc4(1,ip,jp,k) )/4
     curly=vxz-vzx
     write(553) sign(sqrt(abs(curly)),curly)
  enddo
enddo
!print*,'curl:',maxval(curly),minval(curly)
end subroutine writecurly2

subroutine writecurlz(uc4,nx,ny)
integer nx,ny
real, dimension(4,nx,ny,2) :: uc4    !z
integer i,j,k,ip,jp,kp
real vyx,vxy
k=1!z
kp=2!z
do j=1,ny!y
   jp=mod(j,ny)+1!y
   do i=1,nx        !x
      ip=mod(i,nx)+1           !x
  vyx=(uc4(2,ip,j,k)-uc4(2,i,j,k)+uc4(2,ip,jp,k)-uc4(2,i,jp,k) + &
      uc4(2,ip,j,kp)-uc4(2,i,j,kp)+uc4(2,ip,jp,kp)-uc4(2,i,jp,kp) )/4
  vxy=(uc4(1,i,jp,k)-uc4(1,i,j,k)+uc4(1,ip,jp,k)-uc4(1,ip,j,k) + &
      uc4(1,i,jp,kp)-uc4(1,i,j,kp)+uc4(1,ip,jp,kp)-uc4(1,ip,j,kp) )/4
  curlz=vyx-vxy
  write(521) sign(sqrt(abs(curlz)),curlz)
enddo
enddo
end subroutine writecurlz

subroutine writegrad(u,n,ax,id)
implicit none
integer n
real u(4,n,n,n),gradu(3,n,n,n)
integer i,j,k,id
character(1) ax
gradu(1,:,:,:)=u(1,:,:,:)-cshift(u(1,:,:,:),-1,1)
gradu(2,:,:,:)=u(1,:,:,:)-cshift(u(1,:,:,:),-1,2)
gradu(3,:,:,:)=u(1,:,:,:)-cshift(u(1,:,:,:),-1,3)
if (ax=='y') then
   write(530) gradu(1,:,id,:)
   write(531) gradu(2,:,id,:)
   write(532) gradu(3,:,id,:)
   write(533) sum(gradu(:,:,id,:)**2,1)
else if (ax=='x') then
   write(530) gradu(1,id,:,:)
   write(531) gradu(2,id,:,:)
   write(532) gradu(3,id,:,:)
   write(533) sum(gradu(:,id,:,:)**2,1)
else if (ax=='z') then
   write(530) gradu(1,:,:,id)
   write(531) gradu(2,:,:,id)
   write(532) gradu(3,:,:,id)
   write(533) sum(gradu(:,:,:,id)**2,1)
endif

end subroutine writegrad

end module writefile
