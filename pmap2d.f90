  integer, parameter :: nx=512,nz=nx,nt=151
  real, dimension(nx,nz) :: u0
  character*4 fn3
  character(20) fnl
  character(4) npd
  logical,parameter :: v=1
  character(1), parameter :: ax='y'
  integer cp,check
  write(npd,'(i4)') nx
  if (v == 1) then
     if (ax == 'x') fnl='u0_'//trim(adjustl(npd))//'_129x.dat'
     if (ax == 'y') fnl='u0_'//trim(adjustl(npd))//'_256y.dat'
     if (ax == 'z') fnl='u0_'//trim(adjustl(npd))//'_103z.dat'
  else
     if (ax == 'x') fnl='ucrl_'//trim(adjustl(npd))//'_129x.dat'
     if (ax == 'y') fnl='ucrl_'//trim(adjustl(npd))//'_256y.dat'
     if (ax == 'z') fnl='ucrl_'//trim(adjustl(npd))//'_103z.dat'
  endif
  
  open(12,file=trim(adjustl(fnl)),access='stream')
  rmax=0
  rmin=1000
!  do it=1,nt
!  read(12) u0
!  rmax=max(rmax,maxval(u0))
!  rmin=min(rmin,minval(u0))
!  enddo
!  rewind(12)
!  write(*,*) rmax,rmin
  do it=1,nt
  write(fn3,'(i4.4)')it-1
  read(12) u0
!  if (it>=220) then
     cp=0
!     write(fn3,'(i3.3)')it-220
     if (ax == 'x') u0=cshift(u0,-153,2)
     if (ax == 'y') u0=cshift(cshift(u0,-127,1),-153,2)
     if (ax == 'z') u0=cshift(u0,-127,1)
!     if (ax == 'x') u0=cshift(u0,-40,2)
!     if (ax == 'y') u0=cshift(cshift(u0,-17,1),-40,2)
!     if (ax == 'z') u0=cshift(u0,-17,1)
     if (it == 118) cp=1
     if (v == 1 ) then
        call pmap('movies'//'/u0'//'.'//fn3//'.pgm',u0,nx,nz,1,cp,check)
     else
        call pmap('movies'//'/ucrl'//ax//'.'//fn3//'.pgm',u0,nx,nz,1,cp,check)
     endif
!  endif
  enddo

contains
  subroutine pmap(fn,rmap1,nx,ny,iscale,cp,check)
  real rmap(nx,ny),rmap1(nx,ny)
  integer*2, dimension(nx,ny) :: imap
  integer*1, dimension(nx,ny) :: imap1
  integer*2, dimension(nx,nz) :: tmp
  character(len=*):: fn
  integer npix,mypos,cp,i,check

  npix=min(ny/2-1,nx/2-1,300)
  
  
  rmap=rmap1
  iscale1=iscale
  do while (iscale1 > 1)      
     rmap=sign((sqrt(abs(rmap))),rmap)
     iscale1=iscale1-1
  end do
  rmax=maxval(rmap)
  rmin=minval(rmap)
  write(*,*) trim(fn),rmax,rmin
  imap=255-255*(rmap-rmin)/(rmax-rmin)
  imap1=127*(rmap-rmin)/(rmax-rmin)
!  tmp=imap
!  imap=0
!imap(239-1:239+1,216-1:216+1)=0
!imap(454-2:454+2,60-2:60+2)=0
!imap(194-3:194+3,33-3:33+3)=0
imap(256,256)=0
  if (cp ==1 ) then
     check=check+1
!     if (check==1) imap(239-1:239+1,216-1:216+1)=0
!     if (check==2) imap(454-2:454+2,60-2:60+2)=0
!     if (check==3) imap(194-3:194+3,33-3:33+3)=0
     imap(1:20,118)=0
     imap(1:20,139)=0
     do i=1,20
        imap(i,118+i)=0
     enddo

  endif
!  imap(239,128)=0  !x
!  imap(1:65,207:257)=tmp(1:65,207:257)      !y
!  imap(1:65,88:168)=tmp(1:65,88:168)        !z
  open(10,file=fn)
  write(10,'(2hP5)')
  write(10,*)nx,ny
  write(10,*) 255
!  write(10,*) 127
  INQUIRE(UNIT=10, POS=mypos)
  close(10)
  open(10,file=fn, access='stream',position='append')
!  write(10,pos=mypos) int(imap,1)
  write(10) int(imap,1)
  close(10)
end subroutine pmap

end program
