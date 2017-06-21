!ifort -axMIC-AVX512 -O3 -qopenmp rad3d.f90 fourn.for  -o rad3d.x
program rad3d
use parameters
use calculation
use writefile
use pencil_fft
call openfile
  cfl=0.7
  write(*,*) 'u: ',u(:,1,1,1)
  write(*,*) 'maxmin u:',maxval(u(1,:,:,:)),minval(u(1,:,:,:))
  t=nt_last
if (output_pdf_ebmode==1 .or. output_ebmode==1) call create_penfft_plan
  do it=nt_last+1,nt_this
     do i=1,n/1000+1
        dt=timestep()
        call sweepx
        call sweepy
        call sweepz
        call sweepz
        call sweepy
        call sweepx
        t=t+dt
     end do
     write(*,*) t,dt,maxval(u(1,:,:,:)),minval(u(1,:,:,:))
!     write(25) maxval(u(1,:,:,:)),maxloc(u(1,:,:,:)),minval(u(1,:,:,:)),minloc(u(1,:,:,:))
!     if (it>=0) then
     call getuc4(u,uc4,nx,ny,nz)
     if (output_pdf_ebmode==1 .or. output_ebmode==1) call writepdfv(uc4,ue,nx,ny,nz,output_pdf_ebmode,output_ebmode)
     if (output_pdfcurl==1 .or. output_findcurl==1 .or. &
          output_pdfcurlsm==1 .or. output_findcurlsm==1 ) then
        call curlv(uc4,curl,nx,ny,nz)
        if (output_pdfcurl==1) call pdf_curl(curl,nx,ny,nz,nbin)
        if (output_findcurl==1) call findcurl(curl,nx,ny,nz)
        if (output_pdfcurlsm==1 .or. output_findcurlsm==1) then
           call smooth(curl,curlsm,nx,ny,nz)
           if (output_pdfcurlsm==1) call pdf_curl_sm(curlsm,nx,ny,nz,nbin)
           if (output_findcurlsm==1) call findcurlsm(curlsm,nx,ny,nz)
        endif
     endif
     if (output_u==1) then
             if (ax=='y') then
           write(522) u(1,:,id,:)
        else if (ax=='x') then 
           write(522) u(1,id,:,:)
        else if (ax=='z') then
           write(522) u(1,:,:,id) 
        endif
     endif
     if (output_curl==1) then
             if (ax=='y') then
           call writecurly(uc4(:,:,id:id+1,:),nx,nz)
        else if (ax=='x') then 
           call writecurlx(uc4(:,id:id+1,:,:),ny,nz)
        else if (ax=='z') then
           call writecurlz(uc4(:,:,:,id:id+1),nx,ny)
        endif
     endif
     if (output_shift==1 .or. output_shiftsm==1) call writeshift(u,uc4,ns,nx,ny,nz,output_shift,output_shiftsm)
  end do
if (output_pdf_ebmode==1 .or. output_ebmode==1) call destroy_penfft_plan
if (output_backup==1)  write(520) u
  
contains
  
!  subroutine init
!!    pi=4*atan(1.)
!    u=0
!    u(1,:,:,:)=1
!    u(1,1:50,:,:)=u(1,1:50,:,:)+1
!    u(1,:,1:50,:)=u(1,:,1:50,:)+1
!    u(1,:,:,1:50)=u(1,:,:,1:50)+1
!!    u(1,:,1,:)=1+0.1*spread((/(sin(2*pi*(i-0.5)/nx),i=1,nx)/),2,nz)&
!!    +0.1*spread((/(sin(2*pi*(i-0.5)/nz),i=1,nz)/),1,nx)
!    u(2,:,:,:)=0
!  end subroutine init

  subroutine pspect
    implicit none
    complex ctmp4(nx,ny,nz,4)
    integer nn(3)
    integer i,j,k,ii,jj,kk,k2,ik
    real ak,var,pi,ak3(3)
    real, dimension(nx+ny+nz):: pk,pkv,pkdiv,w
    nn=(/nx,ny,nz /)
    
    !$omp parallel do default(shared)
    do i=1,4
       ctmp4(:,:,:,i)=u(i,:,:,:)
!       call fourn(ctmp4(:,:,:,i),nn,3,1)
    end do
    ctmp4=ctmp4/(nx*ny*nz)
    pk=0
    pkdiv=0
    pkv=0
    w=0
!$omp parallel do default(none) shared(ctmp4) private(k,kk,j,jj,i,ii,k2,ak,ak3,ik) reduction(+:pk,pkv,pkdiv,w)
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
             ak=sqrt(k2*1.)
             ik=nint(ak)
             if (ik < 1) cycle
             ak3=(/ ii,jj,kk/)/ak
             pk(ik)=pk(ik)+abs(ctmp4(i,j,k,1))**2*ak**3
             pkdiv(ik)=pkdiv(ik)+abs(sum(ctmp4(i,j,k,2:)*ak3))**2*ak**5
             pkv(ik)=pkv(ik)+sum(abs(ctmp4(i,j,k,2:))**2)*ak**5
             w(ik)=w(ik)+1
          end do
       end do
    end do
    where(w>0) 
       pk=pk/w
       pkv=pkv/w
       pkdiv=pkdiv/w
    end where
    do k=nx+ny+nz,1,-1
       if (w(k)>0) exit
    enddo
    write(*,*) it,k
    write(18) pk(:k),pkdiv(:k),pkv(:k)
    flush(18)

  end subroutine pspect

  real function timestep()
!!$omp workshare
    timestep= cfl/maxval(1./sqrt(3.)+sqrt(sum(u(2:4,:,:,:)**2,1)))   !!!! 1/2
!!$omp end workshare
    return
  end function timestep
    
  subroutine sweepx
    implicit none
    integer j,k
    real u1d(neq,nx)

    !$omp parallel do default(shared) private(j,k,u1d)
    do k=1,nz
       do j=1,ny
          u1d=u(:,:,j,k)
          call relaxing(u1d,nx)
          u(:,:,j,k)=u1d
       end do
    end do
    return
  end subroutine sweepx

  subroutine sweepy
    implicit none
    integer i,k
    real u1d(neq,ny)

    !$omp parallel do default(shared) private(i,k,u1d)
    do k=1,nz
       do i=1,nx
          u1d((/1,3,2,4/),:)=u(:,i,:,k)
          call relaxing(u1d,ny)
          u(:,i,:,k) = u1d((/1,3,2,4/),:)
       end do
    end do
    return
  end subroutine sweepy
  subroutine sweepz
    implicit none
    integer i,j
    real u1d(neq,nz)
    !$omp parallel do default(shared) private(i,j,u1d)
    do j=1,ny
       do i=1,nx
          u1d((/1,4,3,2/),:)=u(:,i,j,:)
          call relaxing(u1d,nz)
          u(:,i,j,:)=u1d((/1,4,3,2/),:)
       end do
    end do
    return
  end subroutine sweepz
  
  recursive subroutine relaxing(u,n)
    implicit none
    integer n
    real cmax
    real, dimension(n) :: c
    real, dimension(neq,n) :: u,u1,w,fu,fr,fl,dfl,dfr

    cmax=averageflux(u,w,c,n)
    fr=(u*spread(c,1,neq)+w)/2
    fl=cshift(u*spread(c,1,neq)-w,1,2)/2
    fu=fr-fl
    u1=u-(fu-cshift(fu,-1,2))*dt/2

    cmax=averageflux(u1,w,c,n)
    fr=(u1*spread(c,1,neq)+w)/2
    dfl=(fr-cshift(fr,-1,2))/2
    dfr=cshift(dfl,1,2)
    call minmod(fr,dfl,dfr,n)

    fl=cshift(u1*spread(c,1,neq)-w,1,2)/2
    dfl=(cshift(fl,-1,2)-fl)/2
    dfr=cshift(dfl,1,2)
    call minmod(fl,dfl,dfr,n)

    fu=fr-fl
    u=u-(fu-cshift(fu,-1,2))*dt
    return
  end subroutine relaxing
  recursive real function averageflux(u,w,c,n)
    implicit none
    integer n
    real, dimension(n)::rho,rhog2,v2,v,w(neq,n),c(n)
    real u(neq,n)

    ! relativistic ideal fluid
    ! stress-energy tensor T=(\rho+p)u^\mu u^\nu + p g^\mu\nu
    ! metric convention g=diag(-1,1,1,1)
    ! Lorenz factor 1/\gamma^2=1-\beta^2
    ! 4-velocity u0=\gamma,u1=\gamma\beta, u_\mu u^\mu=-1
    ! equation of state: p=rho/3
    ! conservation law: T^00,0 + T^01,1=0 (continuity)
    ! T^01,0+T^11,1 = 0 (momentum)
    ! define: U0=T^00=(\rho+p)\gamma^2-p, U1=T^01 = (\rho+p)\gamma^2\beta
    ! simplify: U0=\rho(4\gamma^2-1)/3,U1=4\rho\gamma^2\beta/3
    ! U0=\rho(4/(1-\beta^2)-1)/3,U1=4\rho\beta/3/(1-\beta^2)
    ! U1/U0=4\beta/3+O(\beta^3)
    ! T11=\rho(4\gamma^2\beta^2+1)/3=U1^2/U0+p+O(beta^3)
    
    
    
    v2=3*sum(u(2:4,:)**2,1)/4
    rho=2*sqrt(u(1,:)**2-v2)-u(1,:)
    rhog2=rho/4+3*u(1,:)/4
!    v=u(2,:)/rhog2
    v = 1.5*u(2,:)*rho**0.25/sqrt(rhog2)
    c=(1/sqrt(3.)+abs(v))
    averageflux=maxval(c)
    w(1,:)=u(2,:)
    w(2,:)=3*u(2,:)**2/(3*u(1,:)+rho)+rho/3
    w(3,:)=3*u(2,:)*u(3,:)/(3*u(1,:)+rho)
    w(4,:)=3*u(2,:)*u(4,:)/(3*u(1,:)+rho)
    return

  end function averageflux
  
  subroutine vanleer(f,a,b,n)
    implicit none
    integer n
    real, dimension(neq,n) :: f,a,b,c
    
    c=a*b
    where(c>0) f=f+2*c/(a+b)
    return
  end subroutine vanleer
  recursive subroutine minmod(f,a,b,n)
    implicit none
    integer n
    real, dimension(neq,n) :: f,a,b

    f=f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))/2.
    return
  end subroutine minmod

end program
