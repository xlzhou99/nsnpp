
  program main
! second order
      use Input_var
      implicit double precision (a-h,o-z)
    
  real (kind=8),dimension((n2-1)*n2/2)       :: ed,en
  real (kind=8),dimension(0:n2)              :: aa,ww,wd,xn,wtn,wn
  real (kind=8),dimension(0:m2)              :: wm,xm,wtm
  real (kind=8),dimension((n2+1)*(n2+1))     :: an,bn
  real (kind=8),dimension((m2+1)*(m2+1))     :: am,bm
  real (kind=8),dimension(4*(n2+1)*(m2+1))   :: wkd,wkp,wkp0
  real (kind=8),dimension(0:n2,9)            :: ud,vd
  real (kind=8),dimension(0:n2,0:n2)         :: dn
  real (kind=8),dimension(0:m2,0:m2)         :: dm
  real (kind=8),dimension(0:m2,0:m2)         :: em
  real (kind=8),dimension(0:n2,0:m2)  :: u2,v2,p2,f1,f2,fc1,fc2,uex,vex,pex
  real (kind=8),dimension(0:n2,0:m2)  :: ux,vx,ux1,vx1
  real (kind=8),dimension(0:n2,0:m2)  :: u,v,p,u1,v1,p1,up1,un1,phi1,px,py
  real (kind=8),dimension(0:n2,0:m2)  :: un2,up2,uz2,uz1,phi2,fp,fn,unex,upex,phiex
  real (kind=8),dimension(0:n2,0:m2)  :: up2_log,up1_log,un2_log,un1_log

  real (kind=8),dimension(0:n2,9)     :: ud1,vd1
  real (kind=8),dimension(2)          :: rlm,rlp 
  real (kind=8),dimension(6)          :: a
  real (kind=8),dimension(10,10)  :: ord_data,ord_ns,ord_npp
  real*4 tt(2)
  character*16 fname

	call cpu_time(start1)

      pi2=pi*pi
      call leinitz(n,xn,wtn,an,bn,dn)
      call leinitz(m,xm,wtm,am,bm,dm)
      call ledeig(n,wd,ed,wkp)
      call leneig(n,m,wn,en,aa,ww,wkp)
	call cpu_time(finish1)
        print*,"time cost for prepare LGL data is ",finish1-start1
		print*,"iord=",iord
		print*,"dt=",dt
		print*,"T=", tmax
		print*,"Re=", re

  ! stop
      xr = 1.d0
      xl = -1.d0
      yr = 1.d0
      yl = -1.d0
    eps=1.d0; z1=1.d0;z2=-1.d0
!--- the coefficients to Mapping the basis element to objective. -----
      rlm(1) = (xr-xl)/2.d0!  rlm, r:'right'; l:'left'; m : 'minus(-)',
      rlp(1) = (xr+xl)/2.d0!  rlp, r:'right'; l:'left'; p : 'plus(+)'.
      rlm(2) = (yr-yl)/2.d0
      rlp(2) = (yr+yl)/2.d0
!**** set up the parameters
    a=1.d0
    be  = 1.d0/(h*h)
    alpha = re*a(5)*a(5)


	call cpu_time(start) 
    CC1=1.1D0
    t0=0.d0
    t=t0
        do j =0, m
        do i =0, n
            x = xn(i)
            y = xm(j)  
            tmp=cos(pi*x)*cos(pi*y)
            u2(i,j)=pi*sin(2*pi*y)*(sin(pi*x))**2
            v2(i,j)=-pi*sin(2*pi*x)*(sin(pi*y))**2
            p2(i,j)=sin(pi*x)*sin(pi*y)
            up2(i,j)=cc1+tmp
            un2(i,j)=cc1-tmp
            phi2(i,j)=tmp/pi2
       enddo
    enddo

    call nltermc_ns(n,m,n2,m2,h,u2,v2,ux,vx,dn,dm)
      do j=0,m
        do i=0,n
            up2_log(i,j)=log(up2(i,j))
            un2_log(i,j)=log(un2(i,j))
        enddo
      enddo

    eng_ns2=.5d0*(dleproduct_z(n,m,u2,u2,wtn,wtm,2)+dleproduct_z(n,m,v2,v2,wtn,wtm,2))
    eng_npp2=pnp_energy(n,m,n2,up2,un2,phi2,wtn,wtm,dn,dm)    
      tmp=0.0d0
     do j=0,m
        do i=0,n
          cc=up2(i,j)*(up2_log(i,j)-1.d0)+un2(i,j)*(un2_log(i,j)-1.d0)+ &
            & 0.5*( z1*up2(i,j)+z2*un2(i,j))*phi2(i,j)
          tmp=tmp+cc*wtn(i)*wtm(j)
        enddo
     enddo
     eng_npp=tmp
      sv2=sqrt(eng_npp+c0)
      imax=(tmax-t0+1.D-9)/dt
!**** start the main iteration on time ************
!*
!* u2, v2, p2 will be the current value
!* u1, v1, p1 will be the previous value
!* ux, vx     will be the current nonlinear terms
!* ux1, vx1   will be the previous nonlinear terms
!***************************************************
        
      t1=dtime(tt)
      do it=1,imax
         t=t0+it*dt
         if(abs(u2(1,1)+u2(3,4)+v2(2,5)).gt.500.d0) then
            print*,'overflow at t=',t
            stop
         endif
    f1=0.d0;f2=0.d0;fp=0.d0;fn=0.d0
     if((iord==1).or.(it==1))then
      fc1=u2
      fc2=v2
      else
        fc1=2.d0*u2-u1
        fc2=2.d0*v2-v1
      endif
      iwkd=max(it-1,1) 
      call sav_pnp2D_2nd(fc1,fc2,fp,fn,up1,un1,phi1,up1_log,un1_log,& 
                        & up2,un2,phi2,up2_log,un2_log, eps,wtn,an,bn,dn,& 
                        & wtm,am,bm,dm,wd,ed, wn,en,aa,ww,wkp,wkp0,cc_B1,res,eng_npp,iwkd,it,t)
		uz2=0.d0
      call sav_ns2D(cc_B1,res,eng_npp,up2,un2,uz2,phi2,f1,f2, &
                 & u2,v2,p2,ux,vx,u1,v1,p1,ux1,vx1, sv2, sv1, & 
                & wtn,an,bn,dn,wtm,am,bm,dm,wd,ed,wn,en,aa,ww,wkd,wkp0,iwkd,it,t,imax)

!======================================================================     
!     ..............  END  OF the iteration
    enddo
	
    call cpu_time(finish)
    print*, 'main iteration cpu time ', finish - start
    print*,'Average cputime in each iteration',(finish-start)/(imax)
	print*, start, '&', finish,'&',finish - start,'&',(finish-start)/(imax),'&',(finish-start)/(imax)/(N**3)
	

!********* print out the errors ******************
  
     open (10,file='data.m')
     write(10,*) 'u=['
      do j=0,m
         write (10,100) (u2(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'v=['
      do j=0,m
         write (10,100) (v2(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'p=['
      do j=0,m
         write (10,100) (p2(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'up=['
      do j=0,m
         write (10,100) (up2(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'un=['
      do j=0,m
         write (10,100) (un2(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'phi=['
      do j=0,m
         write (10,100) (phi2(i,j),i=0,n)
      enddo
      write(10,*) '];'

      write(10,*) 'x=['
      write (10,100) (xn(i),i=0,n)
      write(10,*) '];'

      write(10,*) 'y=['
      write (10,100) (xm(i),i=0,m)
      write(10,*) '];'
!   do j=0,m
!         do i=0,n
!            print*,i,j, pex(i,j)-p2(i,j)
!         enddo
!      enddo



      call  error( n, m, wtn,wtm, u1,u2, resu, res2)
      call  error( n, m, wtn,wtm, v1,v2, resv, res2)
      print *, 'Vel-Error L2 =',sqrt(resu**2+resv**2)
      call  error( n, m, wtn,wtm, p1,p2, res1, res2)
      print *, 'p Error  L2 =', res1 ! , 'max =',res2
      call  error( n, m, wtn,wtm, up1,up2, res1, res2)
      print *, 'c1 Error L2 =', res1 ! , 'max =',res2
      call  error( n, m, wtn,wtm, un1,un2, res1, res2)
      print *, 'c2 Error L2 =', res1  !, 'max =',res2
      call  error( n, m, wtn,wtm, phi1,phi2, res1, res2)
      print *, 'phi Error L2 =', res1 ! , 'max =',res2

      call ledxc(n,m,u2,n2,f1,n2,dn)
      call ledyc(n,m,v2,n2,f2,n2,dm)
      do j=0,m
         do i=0,n
            u(i,j) =  f1(i,j)+f2(i,j)
         enddo
      enddo
      v=0
      call  error( n, m, wtn,wtm, u,v, res1, res2)
      print*,'L2 norm of  div u=',res1

        print*,"iord=",iord
      print*,"dt=",dt
      print*,"T=", tmax
	  print*,"N=", N



  stop "done!"





!    print*,"iw=", iw
  100  format (4(1pg14.6),'...')
  200  format(50(1pg12.5,2x))


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     ..............  END  OF  PROGRAMME .......................
!  
!
  stop
  end program main
!
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

 !======================================================  

      subroutine sav_ns2D(cc_B1,res,eng_npp,up2,un2,uz2,phi2,f,g, &
                   & u2,v2,p2,ux,vx,u1,v1,p1,ux1,vx1,sv2,sv1,  & 
                  & wtn,an,bn,dn, wtm,am,bm,dm,wd,ed,wn,en,aa,ww,wkd,wkp0,iwkd,it,t,imax)
 !* Navier-Stokes equations in 2 Dimensional domain
 !*   i.e.  u_t+u*u_x+v*u_y+p_x-\nu*(u_{xx}+u_{yy})=f
 !*         v_t+u*v_x+v*v_y+p_y-\nu*{v_{xx}+v_{yy}}=g
!*with nonhomogeneous Dirichlet boundary condtions
 !****************************************************
!***************************************************** !
      use Input_var

      implicit double precision (a-h,o-z)    
    
      real (kind=8),dimension(0:n2,0:m2)         :: u1,v1,p1,u2,v2,p2,f,g,up2,un2,phi2,uz2,ux,vx,ux1,vx1
      real (kind=8),dimension((n2-1)*n2/2)       :: ed,en
      real (kind=8),dimension(0:n2)              :: aa,ww,wd,xn,wtn,wn
      real (kind=8),dimension(0:m2)              :: wm,xm,wtm
      real (kind=8),dimension((n2+1)*(n2+1))     :: an,bn
      real (kind=8),dimension((m2+1)*(m2+1))     :: am,bm
      real (kind=8),dimension(4*(n2+1)*(m2+1))   :: wkd,wkp0
      real (kind=8),dimension(0:n2,9)            :: ud,vd
      real (kind=8),dimension(0:n2,0:n2)         :: dn
      real (kind=8),dimension(0:m2,0:m2)         :: dm
      real (kind=8),dimension(0:n2,9)     :: ud1,vd1
      real (kind=8),dimension(6)          :: a!x \in (a(2),a(1)), y \in (a(4),a(3)),
      real (kind=8),dimension(0:10)       :: d,b
      real*4 tt(2)
      real (kind=8),allocatable :: ut1(:,:),ut2(:,:),vt1(:,:),vt2(:,:),uti(:,:),vti(:,:),uxbar(:,:),vxbar(:,:)
      real (kind=8),allocatable :: u(:,:),v(:,:),p(:,:),wk1(:,:),wk2(:,:),wk3(:,:),wk4(:,:),fc1(:,:),fc2(:,:)


      allocate(wk3(0:n2,0:m2),wk4(0:n2,0:m2),ut1(0:n2,0:m2),ut2(0:n2,0:m2),vt1(0:n2,0:m2),vt2(0:n2,0:m2),uti(0:n2,0:m2),vti(0:n2,0:m2))
      allocate(u(0:n,0:m),v(0:n,0:m),uxbar(0:n,0:m),vxbar(0:n,0:m),p(0:n,0:m),wk1(0:n2,0:m2),wk2(0:n2,0:m2))
      allocate(fc1(0:n,0:m),fc2(0:n,0:m))



!**** set up the parameters
    a=1
!    h   = a(6)/a(5)! a(5)=(\Delta x)/2, a(6)=(\Delta y)/2
    be  = 1.d0/(h*h)
    alpha = re*a(5)*a(5)
    b3  = 1.d0/dt
    b1=3.d0/2.d0/dt
    a3=b3*alpha
    a1=b1*alpha

!********************************************************* 
! set some parameters:
    z1=1.d0; z2=-1.d0; z3=2.d0
    wk3=(z1*up2+z2*un2+z3*uz2)
    call ledxc(n,m,phi2,n2,wk1,n2,dn)
    call ledyc(n,m,phi2,n2,wk2,n2,dm)
!       
    do j =0, m
           do i =0, n
               x = xn(i)
               y = xm(j)
               fc1(i,j)= wk3(i,j)*wk1(i,j) 
               fc2(i,j)= wk3(i,j)*wk2(i,j)
            end do
            end do 

        if((iord==1).or.(it==1))then
             rt=b3; rm= rt*alpha
            u=u2; v=v2
            uxbar=ux; vxbar=vx
        else
            rt=3.d0/2.d0*b3; rm= rt*alpha
            u=2.d0*u2-.5d0*u1; v=2.d0*v2-.5d0*v1
            uxbar=2.d0*ux-ux1; vxbar=2.d0*vx-vx1
        endif
        call ledxc(n,m,p2,n2,wk1,n2,dn)
        call ledyc(n,m,p2,n2,wk2,n2,dm)
!!        do j=0,m
!!        do j=0,m
!!            do i=0,n
!!              u(i,j)=alpha*( f(i,j)-fc1(i,j)+b3*u(i,j) - wk1(i,j))
!!              v(i,j)=alpha*( g(i,j)-fc2(i,j)+b3*v(i,j) - wk2(i,j))
!!              wk1(i,j)=-alpha*uxbar(i,j)
!!              wk2(i,j)=-alpha*vxbar(i,j)
!!
!!              u1(i,j)=u2(i,j)
!!              v1(i,j)=v2(i,j)
!!              p1(i,j)=p2(i,j)
!!              ux1(i,j)=ux(i,j)
!!              vx1(i,j)=vx(i,j)
!!           enddo
!!        enddo
        
        do j=0,m
               do i=0,n
                  f(i,j)=alpha*( f(i,j)+b3*u(i,j) - wk1(i,j)/a(5) )
                  g(i,j)=alpha*( g(i,j)+b3*v(i,j) - wk2(i,j)/h/a(5) )
                  wk3(i,j)=-alpha*(uxbar(i,j)+fc1(i,j))
                  wk4(i,j)=-alpha*(vxbar(i,j)+fc2(i,j))
                  ux1(i,j)=ux(i,j)
                  vx1(i,j)=vx(i,j)
                  u1(i,j)=u2(i,j)
                  v1(i,j)=v2(i,j)
                  p1(i,j)=p2(i,j)
               enddo
            enddo

!******  solving the first equations for u and v*******
!vec(\tilde{u}^{n+1}_1)=(u2,v2); vec(\tilde{u}^{n+1}_2)=(u,v)



    call leptos2(n,m,f,n2,p,bn,bm)
    call leptos2(n,m,g,n2,p,bn,bm)
    call leptos2(n,m,wk3,n2,p,bn,bm)
    call leptos2(n,m,wk4,n2,p,bn,bm)

    call ledhmz2(n,m,rm,be,f,n2,ut1,n2,ed,wd,wkd,iwkd) 
    call ledhmz2(n,m,rm,be,g,n2,vt1,n2,ed,wd,wkd,2)
    call ledhmz2(n,m,rm,be,wk3,n2,ut2,n2,ed,wd,wkd,2)
    call ledhmz2(n,m,rm,be,wk4,n2,vt2,n2,ed,wd,wkd,2)
      uxbar=uxbar+fc1
      vxbar=vxbar+fc2
    call leptos2(n,m,uxbar,n2,p,bn,bm)
    call leptos2(n,m,vxbar,n2,p,bn,bm)
    d3=dlepsw(n,m,1,uxbar,ut1,n2,wtn,wtm)+dlepsw(n,m,1,vxbar,vt1,n2,wtn,wtm)
    d4=dlepsw(n,m,1,uxbar,ut2,n2,wtn,wtm)+dlepsw(n,m,1,vxbar,vt2,n2,wtn,wtm)

    ccc1=dt/sqrt(eng_npp+c0)/2.d0
    if((iord==1).or.(it==1))then
!   s=(exp(t/T_ex)*d3+sav2/dt)/( (T_ex+dt)/(T_ex*dt)-exp(2*t/T_ex)*d4 )*exp(t/T_ex)
    s=( sv2+ccc1*(cc_B1+d3) )/( sqrt(eng_npp+c0)+ccc1*(res-d4) )

    else
      s=(2.d0*sv2-0.5*sv1+ccc1*(cc_B1+d3) )/(3.d0/2.d0*sqrt(eng_npp+c0)+ccc1*(res-d4) )
      
    endif
    sv1=sv2
    sv2=s*dsqrt(eng_npp+c0)
 !  print*,'s =',s, abs(sv2-dsqrt(eng_npp+c0))
 !   print*, 'abs(s-1)='
 !   print*, abs(1-s)
! 
 !   print*, abs(sv2-dsqrt(eng_npp+c0))
         do j=0,m
            do i=0,n
              uti(i,j)=ut1(i,j)+s*ut2(i,j)
              vti(i,j)=vt1(i,j)+s*vt2(i,j)
            enddo
         enddo
    call ledx(n,m,uti,n2,wk1,n2)
    call ledy(n,m,vti,n2,wk2,n2)

          do j=0,m
            do i=0,n
              f(i,j)=(wk1(i,j)+wk2(i,j)/h)
               wk1(i,j)=-rt*f(i,j)

            enddo
         enddo
!****** compute \phi^{n+1}:s=(p^{n+1}-p^n) ********
!**** p2:=\phi^{n+1}_1
!**** p:=\phi^{n+1}_2
         call lenhmz2(n,m,0.d0,be,wk1,n2,p2,n2,en,wn,aa,ww,wkp0,iwkd)

         call ledx(n,m,p2,n2,wk1,n2)
         call ledy(n,m,p2,n2,wk2,n2)

         do j=0,m
            do i=0,n
               uxbar(i,j)=uti(i,j)-wk1(i,j)/rt
               vxbar(i,j)=vti(i,j)-wk2(i,j)/h/rt
             !  p2(i,j)=p2(i,j)+p(i,j)!-v(i,j)/re
            enddo
         enddo
!
!**** compute d(4):=\nu \int_\Omega \nabla \bar{{\mb u}}^{n+1} and d(6):=A_2 ****
        call leptos2(n,m,fc1,n2,wk1,bn,bm)
        call leptos2(n,m,fc2,n2,wk1,bn,bm)
        cc_nps= dlepsw(n,m,1,uti,fc1,n2,wtn,wtm)+dlepsw(n,m,1,vti,fc2,n2,wtn,wtm)
        wk2=p1 
        call leptos2(n,m,wk2,n2,wk1,bn,bm)
    if  ((it.eq.1).or.(iord==1)) then
        do j=0,m
          do i=0,n
              u2(i,j)=uxbar(i,j)
              v2(i,j)=vxbar(i,j)
              p2(i,j)=wk2(i,j)+ p2(i,j)
              u(i,j)=u2(i,j)
              v(i,j)=v2(i,j)
             ! p(i,j)=p2(i,j)
          enddo
        enddo
        else
        do j=0,m
          do i=0,n
              u2(i,j)=uxbar(i,j)
              v2(i,j)=vxbar(i,j)
              p2(i,j)=wk2(i,j)+ p2(i,j)!-f(i,j)/alpha
          enddo
        enddo
        endif


    if((iord==2).and.(it==imax))then
        do j=0,m
          do i=0,n
              p2(i,j)=p2(i,j)-f(i,j)/alpha
          enddo
        enddo
        endif

        do i=0,n
            p2(0,i)  =0.d0
            p2(i,0)  =0.d0
            p2(n-1,i)=0.d0
            p2(i,n-1)=0.d0
            p2(n,i)  =0.d0
            p2(i,n)  =0.d0
          enddo
    call lestop2(n,m,u2,n2,wk1,an,am)
    call lestop2(n,m,v2,n2,wk1,an,am)
    call lestop2(n,m,p2,n2,wk1,an,am)
    call nltermc_ns(n,m,n2,m2,h,u2,v2,ux,vx,dn,dm)
    phi2=s*phi2

    deallocate(ut1,ut2,vt1,vt2,uti,vti,u,v,uxbar,vxbar,wk1,wk2,wk3,wk4,p,fc1,fc2)
    return
  stop
  end subroutine sav_ns2D

!****************************************************
 





    subroutine  sav_pnp2d_2nd(fc1,fc2,f1,f2,u1,v1,p1,u1_log,v1_log,&
                        & u2,v2,p2,u2_log,v2_log,eps,wtn,an,bn,dn, & 
                        &  wtm,am,bm,dm,wd,ed, wn,en,aa,ww,wkp,wkp0,cc_B1,res,c_eng,iwkd,it,t)
! 3 ions
        use Input_var


    implicit double precision (a-h,o-z)    

  real (kind=8),dimension(0:n2,0:m2)         :: u2,v2,p2,f1,f2,fc1,fc2,u1,v1,p1
  real (kind=8),dimension(0:n2,0:m2)         :: u1_log,v1_log,u2_log,v2_log
  real (kind=8),dimension(4*(n2+1)*(m2+1))   :: wkp0,wkp
  real (kind=8),dimension((n2-1)*n2/2)       :: ed,en
  real (kind=8),dimension(0:n2)              :: aa,ww,wd,xn,wtn,wn
  real (kind=8),dimension(0:m2)              :: wm,xm,wtm
  real (kind=8),dimension((n2+1)*(n2+1))     :: an,bn
  real (kind=8),dimension((m2+1)*(m2+1))     :: am,bm
  real (kind=8),dimension(0:n2,0:n2)         :: dn
  real (kind=8),dimension(0:m2,0:m2)         :: dm
  real (kind=8),dimension(6)          :: a!x \in (a(2),a(1)), y \in (a(4),a(3)),
  real*4 tt(2)
  character*16 fname
  real (kind=8),allocatable :: u(:,:),v(:,:),p(:,:),g1(:,:),g2(:,:),wk1(:,:),wk2(:,:)

  allocate(g1(0:n2,0:m2),g2(0:n2,0:m2),wk1(0:n2,0:m2),wk2(0:n2,0:m2))
  allocate(u(0:n2,0:m2),v(0:n2,0:m2),p(0:n2,0:m2))

    z1=1.d0;z2=-1.d0
    D1=1.d0; D2=1.d0
!**** set up the parameters
    a=1
!    h   = a(6)/a(5)! a(5)=(\Delta x)/2, a(6)=(\Delta y)/2
    be  = 1.d0/(h*h)
    alpha = a(5)*a(5)
    b3  = 1.d0/dt
    a3=b3*alpha

    if ((iord==1).or.(it==1)) then
      rm=a3
      do j=0,m
          do i=0,n
               u(i,j)=u2_log(i,j)
               v(i,j)=v2_log(i,j)
               p(i,j)=p2(i,j)
           enddo
        enddo

              call pnp_f2(n,m,n2,h,z1,z2,u,v,p,fc1,fc2,g1,g2,dn,dm)
         do j=0,m
          do i=0,n
            g1(i,j)=g1(i,j)+b3*u2_log(i,j)  +f1(i,j)/u2(i,j)!
            g2(i,j)=g2(i,j)+b3*v2_log(i,j)  +f2(i,j)/v2(i,j)!
          enddo
        enddo
      else
        rm=a3*3.d0/2.d0
          do j=0,m
            do i=0,n
               p(i,j)=2.d0*p2(i,j)-p1(i,j)
               u(i,j)=2.d0*u2_log(i,j)-u1_log(i,j)
               v(i,j)=2.d0*v2_log(i,j)-v1_log(i,j)
           enddo
        enddo
      call pnp_f2(n,m,n2,h,z1,z2,u,v,p,fc1,fc2,g1,g2,dn,dm)
      do j=0,m
          do i=0,n
            g1(i,j)=g1(i,j)+b3*(2.d0*u2_log(i,j)-.5*u1_log(i,j)) +f1(i,j)/(2*u2(i,j)-u1(i,j))!
            g2(i,j)=g2(i,j)+b3*(2.d0*v2_log(i,j)-.5*v1_log(i,j)) +f2(i,j)/(2*v2(i,j)-v1(i,j))!
          enddo
        enddo

    endif
        u1=u2
        v1=v2
        p1=p2
        u1_log=u2_log
        v1_log=v2_log
        
        call leptos2(n,m,g1,n2,v,bn,bm)
        call leptos2(n,m,g2,n2,v,bn,bm)
         call lenhmz2(n,m,rm,be,g1,n2,u2_log,n2,en,wn,aa,ww,wkp,iwkd)
         call lenhmz2(n,m,rm,be,g2,n2,v2_log,n2,en,wn,aa,ww,wkp,2)
         call lestop2(n,m,u2_log,n2,v,an,am)
         call lestop2(n,m,v2_log,n2,v,an,am)
!      

        do j=0,m
          do i=0,n  
            u2(i,j)=exp(u2_log(i,j))
            v2(i,j)=exp(v2_log(i,j))
  !          print*,i,u2(i,j),v2(i,j)
          enddo
        enddo

!        stop "message"
            cc1=dleproduct_z(n,m,u1,f,wtn,wtm,1)/dleproduct_z(n,m,u2,f,wtn,wtm,1)
            cc2=dleproduct_z(n,m,v1,f,wtn,wtm,1)/dleproduct_z(n,m,v2,f,wtn,wtm,1)
             do j=0,m
          do i=0,n           
            u2(i,j)=cc1*u2(i,j)
            v2(i,j)=cc2*v2(i,j)
            u(i,j)=z1*u2(i,j)+z2*v2(i,j)  ! +f3(i,j)
          enddo
        enddo
        
        do j=0,m
        do i=0,n
            u2_log(i,j)=log(u2(i,j))
            v2_log(i,j)=log(v2(i,j))
        enddo
      enddo

       call leptos2(n,m,u,n2,v,bn,bm)

        call lenhmz2(n,m,0.d0,1.d0,u,n2,p2,n2,en,wn,aa,ww,wkp0,iwkd)
        call lestop2(n,m,p2,n2,v,an,am)
        p2=p2/eps

          do j=0,m
          do i=0,n
            g1(i,j)=u2_log(i,j)+z1*p2(i,j)!\mu_1
            g2(i,j)=v2_log(i,j)+z2*p2(i,j)!\mu_2
          enddo
        enddo

!*** compute \sum_i \int_\Omega \mu_i^{n+1}*f_i **********************************
       cc_B1= dleproduct_z(n,m,g1,f1,wtn,wtm,2)+dleproduct_z(n,m,g2,f2,wtn,wtm,2)

    tmp=0.0d0
     do j=0,m
        do i=0,n
          cc=u2(i,j)*(u2_log(i,j)-1.d0)+v2(i,j)*(v2_log(i,j)-1.d0)+ &
            & 0.5*( z1*u2(i,j)+z2*v2(i,j))*p2(i,j)
          tmp=tmp+cc*wtn(i)*wtm(j)
        enddo
     enddo
     c_eng=tmp
!*** compute (c1, |\nabla\mu_1|^2)+(c2,|\nabla\mu_2|^2) ***
        res1= pnp_edr(n,m,u2,n2,g1,wtn,wtm,dn,dm)
        res2= pnp_edr(n,m,v2,n2,g2,wtn,wtm,dn,dm)
        res=res1+res2
!**********************************************************
    
        deallocate(u,v,p,g1,g2,wk1,wk2)
        return
       stop 
    end subroutine   sav_pnp2d_2nd







    subroutine pnp_f2(n,m,n2,h,d1,d2,u,v,p,fc1,fc2,g1,g2,dn,dm)
!*****************************************************
!* computer 
!*    u_x*u_x + 1/(h*h)*u_y*u_y +  d1*((u v_x+u v_y/h)_x+(u v_x+u v_y/h)_y)& 
!*    +d1*(v_xx+1/(h*h)*v_yy)
!*        input:  u,w,v:  physical value, nochange on output.
!*        output: wk: physical value 
!*
!*****************************************************
   
    implicit double precision (a-h, o-z)
    real (kind=8),dimension(0:n,0:m) ::u,v,p,fc1,fc2,g1,g2
      real (kind=8),dimension(0:n,0:n)  :: dn
       real (kind=8),dimension(0:n,0:n)  :: dm
  real (kind=8),allocatable :: ux(:,:),uy(:,:),px(:,:),py(:,:),pxx(:,:),pyy(:,:)
    
    allocate(ux(0:n,0:m),uy(0:n,0:m),px(0:n,0:m),py(0:n,0:m),pxx(0:n,0:m),pyy(0:n,0:m))

!
        call ledxc(n,m,p,n2,px,n2,dn)
        call ledyc(n,m,p,n2,py,n2,dm)
!
        call ledxc(n,m,px,n2,pxx,n2,dn)
      call ledyc(n,m,py,n2,pyy,n2,dm)
!
       call ledxc(n,m,u,n2,ux,n2,dn)
        call ledyc(n,m,u,n2,uy,n2,dm)
!
    do j=0,m
        do i=0,n
        g1(i,j)= ux(i,j)**2+uy(i,j)**2  + & 
            & d1*( ux(i,j)*px(i,j)+uy(i,j)*py(i,j)+pxx(i,j)+pyy(i,j) )
        enddo
        enddo
!
    call ledxc(n,m,v,n2,ux,n2,dn)
        call ledyc(n,m,v,n2,uy,n2,dm)
!
    do j=0,m
        do i=0,n
        g2(i,j)= ux(i,j)**2+uy(i,j)**2  + & 
            & d2*( ux(i,j)*px(i,j)+uy(i,j)*py(i,j)+pxx(i,j)+pyy(i,j) ) 
        enddo
        enddo

        do j=0,m
        do i=0,n
          ux(i,j)=u(i,j)*fc1(i,j)
          uy(i,j)=u(i,j)*fc2(i,j)
        enddo
        enddo
        call ledxc(n,m,ux,n2,px,n2,dn)
        call ledyc(n,m,uy,n2,py,n2,dm)
        do j=0,m
        do i=0,n
          g1(i,j)=g1(i,j)-(px(i,j)+py(i,j))
        enddo
        enddo

        do j=0,m
        do i=0,n
          ux(i,j)=v(i,j)*fc1(i,j)
          uy(i,j)=v(i,j)*fc2(i,j)
        enddo
        enddo
        call ledxc(n,m,ux,n2,px,n2,dn)
        call ledyc(n,m,uy,n2,py,n2,dm)
        do j=0,m
        do i=0,n
          g2(i,j)=g2(i,j)-(px(i,j)+py(i,j))
        enddo
        enddo
    deallocate(ux,uy,px,py,pxx,pyy)

       return
       stop
    
        end subroutine pnp_f2
  !================================================================  


  double precision function pnp_energy(n,m,n2,u,v,w,wtn,wtm,dn,dm)
!********************************************
!*  Compute the discrete Energy of PNP equations
!*
!*   Input:  u,v,w: contains the physical-value of L-G-L points
!*
!*          wtn,wtm contain the weights relative to the discrete products
!*          wtn and wtm can be obtained by calling leinit.f
!*
!*   Output: Value of the discrete Energy of PNP equations!
!*******************************************
  implicit double precision (a-h, o-z)

     real (kind=8),dimension(0:n)       :: wtn
       real (kind=8),dimension(0:m)       :: wtm
       real (kind=8),dimension(0:n,0:n)  :: dn
       real (kind=8),dimension(0:m,0:m)  :: dm
    real (kind=8),dimension(0:n2,0:m) :: u,v,w

  real (kind=8),allocatable :: ux(:,:),uy(:,:)
    
    allocate(ux(0:n2,0:m),uy(0:n2,0:m))

  call ledxc(n,m,w,n2,ux,n2,dn)
  call ledyc(n,m,w,n2,uy,n2,dm)

   tmp=0.d0

     do j=0,m
        do i=0,n
          cc=u(i,j)*(log(u(i,j))-1.d0)+v(i,j)*(log(v(i,j))-1.d0)+ &
            & 0.5*(ux(i,j)**2+uy(i,j)**2)
          tmp=tmp+cc*wtn(i)*wtm(j)
        enddo
     enddo
     pnp_energy=tmp

  deallocate(ux,uy)


      return
      end function pnp_energy
!================================================================




        double precision function pnp_edr(n,m,u,n2,v,wtn,wtm,dn,dm)
!*****************************************************
!* compute
!*    res=(u,|\nabla v|^2)
!*        input:  u,v:  physical value, nochange on output.
!*        output: wk: physical value 
!*
!*****************************************************
   
    implicit double precision (a-h, o-z)
    real (kind=8),dimension(0:n2,0:m) ::u,v
    real (kind=8),dimension(0:n)       :: wtn
       real (kind=8),dimension(0:m)       :: wtm
       real (kind=8),dimension(0:n,0:n)  :: dn
       real (kind=8),dimension(0:m,0:m)  :: dm
  real (kind=8),allocatable :: v_y(:,:),wk(:,:)
!
    allocate(v_y(0:n2,0:m),wk(0:n2,0:m))
!
    call ledxc(n,m,v,n2,wk,n2,dn)
    call ledyc(n,m,v,n2,v_y,n2,dm)
        res=0.d0
         do j=0,m
            do i=0,n
              cc=wtn(i)*wtm(j)
            res=res+u(i,j)*(wk(i,j)**2+v_y(i,j)**2)*cc
          enddo
        enddo
        pnp_edr=res
  deallocate(v_y,wk)
  return
  stop

  end function pnp_edr



        double precision function dleproduct_z(n,m,u,v,wtn,wtm,ikw)
!********************************************
!*  Compute the discrete inner product of two given functions
!*
!*   Input:  u,v: contains the physical-value of L-G-L points
!*
!*          wtn,wtm contain the weights relative to the discrete products
!*          wtn and wtm can be obtained by calling leinit.f
!*
!*   Output: res: Value of the discrete inner product of u and v!
!*******************************************
        implicit double precision (a-h,o-z)
        real (kind=8),dimension(0:n)       :: wtn
        real (kind=8),dimension(0:m)       :: wtm
        real (kind=8),dimension(0:n,0:m)   :: u,v
!
        res = 0.

  if (ikw==1)then
     do j =0, m
        do i =0, n
          res = res + u(i,j)*wtn(i)*wtm(j)
        end do
        end do
     else


        do j =0, m
        do i =0, n
          res = res + u(i,j)*v(i,j)*wtn(i)*wtm(j)
        end do
        end do
    endif
    dleproduct_z=res
!
        return
  end function dleproduct_z

