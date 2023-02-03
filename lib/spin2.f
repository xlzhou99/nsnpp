      program spin2
****************************************************
*   Second order BDF projection method for time-dependent
*   Navier-Stokes equations: Application to spin-down problem !
*  
* Notes:
*  1. Must be compiled with LAPACK and BLAS (available from www.netlib.org)
*  2. Reynolds # = 2 * re
*  3. Real time = 2 * t
*  4. The results can be viewed in matlab by typing "stream".
*
* Author: Jie Shen
*
* References:
*  on spectral method:
*   Jie Shen, SIAM J. Sci. Comput., v. 15, 1489-1505 (1994)
*  on projection method:
*   J.L. Guermond, J. Shen. "On the error estimates of rotational
*    pressure-correction projection methods" Math. Comput., 
*    to appear, available from   my web page
*****************************************************
      parameter (n2=128,m2=128,nn=(n2+1)*(n2+1),mm=(m2+1)*(m2+1))
      parameter (nm=(n2+1)*(m2+1))
      implicit double precision (a-h,o-z)
      dimension u(0:n2,0:m2),v(0:n2,0:m2),ux(0:n2,0:m2)
      dimension vx(0:n2,0:m2),u1(0:n2,0:m2),u2(0:n2,0:m2)
      dimension v1(0:n2,0:m2),v2(0:n2,0:m2)
      dimension ed((n2-1)*n2/2),en((n2-1)*n2/2),aa(0:n2),ww(0:n2)
      dimension wd(0:n2),wn(0:n2),wk1(0:n2,0:m2),wk2(0:n2,0:m2)
      dimension xn(0:n2),xm(0:m2),an(nn),bn(nn),am(mm),bm(mm),wkd(4*nm)
      dimension p1(0:n2,0:m2),p2(0:n2,0:m2),wkp(4*nm)
      dimension ux1(0:n2,0:m2),vx1(0:n2,0:m2)
      dimension wtn(0:n2),wtm(0:m2)
      real*4 tt(2)
      character*16 fname
      pi=acos(-1.0d0)
      open (20,file='data')
      read (20,*)  n, m, re, dt, tmax,h,irestart,fname,ifinal
      open(50,file=fname)

      call leinit(n,xn,wtn,an,bn)
      call leinit(m,xm,wtm,am,bm)
      call ledeig(n,wd,ed,wkp)
      call leneig(n,m,wn,en,aa,ww,wkp)

**** set up the parameters
      al=3*re/(2*dt)
      be=1.d0/(h*h)
      a2=re/(2*dt)
      a3=re/dt


      if(irestart.eq.0) then
****  Determine the initial condition for the velocity ***
         do j=0,m
            do i=0,n
               wk1(i,j)=0.d0
            enddo
         enddo
         wk1(0,0)=-2.d0
         call ledhmz2(n,m,0.d0,be,wk1,n2,wk2,n2,ed,wd,wkd,1)
         call ledy(n,m,wk2,n2,u2,n2)
         call ledx(n,m,wk2,n2,v2,n2)
*******Inital data at time step 0 **************
         do j=0,m
            do i=0,n
               p2(i,j)=0.d0
               u2(i,j)=u2(i,j)/h
               v2(i,j)=-v2(i,j)
               u(i,j)=u2(i,j)
               v(i,j)=v2(i,j)
            enddo
         enddo
***** compute the nonlinear terms ******
         call nlterm(n,m,h,n2,m2,u,v,ux,vx,wk1,wk2,an,am,bn,bm)

      else
         read (50,*) n,m,re0,dt0,t0
         print*,'using previous data with n1=',m,'  t0=',t0
         if ((re.ne.re0).or.(dt.ne.dt0)) then
            print*,'WARNING: different Reynolds number or time step'
            stop
         endif
         read (50,*) ((u2(i,j),i=0,n),j=0,m)
         read (50,*) ((v2(i,j),i=0,n),j=0,m)
         read (50,*) ((p2(i,j),i=0,n),j=0,m)
         read (50,*) ((u1(i,j),i=0,n),j=0,m)
         read (50,*) ((v1(i,j),i=0,n),j=0,m)
         read (50,*) ((ux(i,j),i=0,n),j=0,m)
         read (50,*) ((vx(i,j),i=0,n),j=0,m)
         read (50,*) ((ux1(i,j),i=0,n),j=0,m)
         read (50,*) ((vx1(i,j),i=0,n),j=0,m)
      endif
      imax=(tmax-t0+1.D-9)/dt


**** start the main iteration on time *******
*
* u2, v2, p2 will be the current value
* u1, v1,    will be the previous value
* ux, vx,    will be the current nonlinear terms
* ux1, vx1   will be the previous nonlinear terms
*********************************************
      t1=dtime(tt)
      do it=1,imax
         t=t0+it*dt
         if(abs(u2(1,1)+u2(3,4)+v2(2,5)).gt.500.d0) then
            print*,'overflow at t=',t
            stop
         endif

******* compute u= -re p_x   + a2*(4u^n-u^{n-1}) + nonlinear term, 
******          v= -re p_y/h + a2*(4v^n-v^{n-1}) + nonlinear term, 


         call ledx(n,m,p2,n2,wk1,n2)
         call ledy(n,m,p2,n2,wk2,n2)

         if ((irestart.ne.1).and.(it.eq.1)) then
            do j=0,m
               do i=0,n
                  u(i,j)=-re*(ux(i,j)+wk1(i,j))  +a3*u2(i,j)
                  v(i,j)=-re*(vx(i,j)+wk2(i,j)/h)+a3*v2(i,j)
                  ux1(i,j)=ux(i,j)
                  vx1(i,j)=vx(i,j)
                  u1(i,j)=u2(i,j)
                  v1(i,j)=v2(i,j)
               enddo
            enddo
         else
            do j=0,m
               do i=0,n
                  u(i,j)=-re*(2*ux(i,j)-ux1(i,j)+wk1(i,j))
     1                   +a2*(4*u2(i,j)-u1(i,j))
                  v(i,j)=-re*(2*vx(i,j)-vx1(i,j)+wk2(i,j)/h)
     1                   +a2*(4*v2(i,j)-v1(i,j))
                  ux1(i,j)=ux(i,j)
                  vx1(i,j)=vx(i,j)
                  u1(i,j)=u2(i,j)
                  v1(i,j)=v2(i,j)
               enddo
            enddo
         endif

******  solving the first equations for u and v*******

         if ((irestart.ne.1).and.(it.eq.1)) then
            rm=a3
         else
            rm=al
         endif
         iwkd=max(it-1,1)
         call ledhmz2(n,m,rm,be,u,n2,u2,n2,ed,wd,wkd,iwkd)

         call ledhmz2(n,m,rm,be,v,n2,v2,n2,ed,wd,wkd,2)


***** solving the equation for (p^{n+1}-p^n) ******

******* compute u2_x + v2_y/h *******
         call ledx(n,m,u2,n2,wk1,n2)
         call ledy(n,m,v2,n2,wk2,n2)
         do j=0,m
            do i=0,n
               u(i,j)=wk1(i,j)+wk2(i,j)/h
               v(i,j)=u(i,j)
            enddo
         enddo

****** compute -2*dt/3 (p^{n+1}-p^n) ********

         call lenhmz2(n,m,0.d0,be,u,n2,p1,n2,en,wn,aa,ww,wkp,it)

*********** update u^{n+1}, v^{n+1} and p^{n+1} ******
         call ledx(n,m,p1,n2,wk1,n2)
         call ledy(n,m,p1,n2,wk2,n2)
         if ((irestart.ne.1).and.(it.eq.1)) then
            tmp=1.d0/dt
         else
            tmp=3.d0/(2*dt)
         endif
         do j=0,m
            do i=0,n
               u2(i,j)=u2(i,j)+wk1(i,j)
               v2(i,j)=v2(i,j)+wk2(i,j)/h
               p2(i,j)=p2(i,j)-tmp*p1(i,j)-v(i,j)/re
*               p2(i,j)=p2(i,j)-tmp*p1(i,j)
            enddo
         enddo

      do j=0,m
         do i=0,n
            u(i,j)=u2(i,j)
            v(i,j)=v2(i,j)
         enddo
      enddo
      call nlterm(n,m,h,n2,m2,u,v,ux,vx,wk1,wk2,an,am,bn,bm)

      enddo         
*** end of the main iteration on t ****

      t1=dtime(tt)
      print*,'cpu elapsed=',tt(1)
************************************
      print*,'t=',t,'   dt=',dt,' re=',re
      rewind (50)
      write (50,*) n,m,re,dt,t
         write (50,*) ((u2(i,j),i=0,n),j=0,m)
         write (50,*) ((v2(i,j),i=0,n),j=0,m)
         write (50,*) ((p2(i,j),i=0,n),j=0,m)
         write (50,*) ((u1(i,j),i=0,n),j=0,m)
         write (50,*) ((v1(i,j),i=0,n),j=0,m)
         write (50,*) ((ux(i,j),i=0,n),j=0,m)
         write (50,*) ((vx(i,j),i=0,n),j=0,m)
         write (50,*) ((ux1(i,j),i=0,n),j=0,m)
         write (50,*) ((vx1(i,j),i=0,n),j=0,m)

      if(ifinal.eq.1) then
         write(35,*) ((u2(i,j),i=0,n),j=0,m)
      else
         read (35,*) ((wk2(i,j),i=0,n),j=0,m)
         do j=0,m
            do i=0,n
               wk2(i,j)=wk2(i,j)-u2(i,j)
            enddo
         enddo
         bl2=dlel2no(n,m,1,wk2,n2,wtn,wtm)
         print*,'L2 error of velocity=',bl2
      endif
      do j=0,m
         do i=0,n
            wk1(i,j)=u2(i,j)
            wk2(i,j)=v2(i,j)
         enddo
      enddo
      call lestop2(n,m,wk1,n2,u,an,am)
      call lestop2(n,m,wk2,n2,u,an,am)

      open (10,file='stream.m')

      write(10,*) 'u=['
      do j=0,m
         write (10,100) (wk1(i,j),i=0,n)
      enddo
      write(10,*) '];'
      write(10,*) 'v=['
      do j=0,m
         write (10,100) (wk2(i,j),i=0,n)
      enddo
      write(10,*) '];'

**** Compute the vorticity and the stream function*****
      call ledx(n,m,v2,n2,wk1,n2)
      call ledy(n,m,u2,n2,wk2,n2)
      do j=0,m
         do i=0,n
            wk1(i,j)=(wk1(i,j)-wk2(i,j)/h)*.5d0
            u(i,j)=wk1(i,j)
         enddo
      enddo
      call lestop2(n,m,u,n2,wk2,an,am)
      if(ifinal.eq.1) then
         write(35,*) ((u(i,j),i=0,n),j=0,m)
      else
         read (35,*) ((wk2(i,j),i=0,n),j=0,m)
         bmax=0.d0
         do j=0,m
            do i=0,n
               wk2(i,j)=wk2(i,j)-u(i,j)
               if(bmax.lt.abs(wk2(i,j))) bmax=abs(wk2(i,j))
            enddo
         enddo
         bl2=dlel2no(n,m,0,wk2,n2,wtn,wtm)
         print*,'Max error of vorticity=',bmax,'  L2 error=',bl2
      endif

      bmax=0.d0
      do j=0,m
         do i=0,n
            if(bmax.lt.abs(u(i,j))) bmax=abs(u(i,j))
         enddo
      enddo
      print*,'Max norm of vorticity=',bmax
      write(10,*) 'vor=['
      do j=0,m
         write (10,100) (u(i,j),i=0,n)
      enddo
      write(10,*) '];'
*** Stream function ****
      call ledhmz2(n,m,0.d0,be,wk1,n2,wk2,n2,ed,wd,wkd,1)
      call lestop2(n,m,wk2,n2,wk1,an,am)
      if(ifinal.eq.1) then
         write(35,*) ((wk2(i,j),i=0,n),j=0,m)
      else
         read (35,*) ((wk1(i,j),i=0,n),j=0,m)
         bmax=0.d0
         do j=0,m
            do i=0,n
               wk1(i,j)=wk1(i,j)-wk2(i,j)
               if(bmax.lt.abs(wk1(i,j))) bmax=abs(wk1(i,j))
            enddo
         enddo
         bl2=dlel2no(n,m,0,wk1,n2,wtn,wtm)
         print*,'Max error of Stream=',bmax,'  L2 error=',bl2
      endif

      bmax=0.d0
      do j=0,m
         do i=0,n
            if(bmax.lt.abs(wk2(i,j))) bmax=abs(wk2(i,j))
         enddo
      enddo
      print*,'Max norm of stream fn (by solving a Poisson eqn)=',bmax
      write(10,*) 'str=['
      do j=0,m
         write (10,100) (wk2(i,j),i=0,n)
      enddo
      write(10,*) '];'
******* compute div u= u2_x + v2_y/h *******
      call ledx(n,m,u2,n2,wk1,n2)
      call ledy(n,m,v2,n2,wk2,n2)
      do j=0,m
         do i=0,n
            u(i,j)=wk1(i,j)+wk2(i,j)/h
         enddo
      enddo
      call lestop2(n,m,u,n2,wk1,an,am)
      divl2= dlel2no(n,m,0,u,n2,wtn,wtm)
      print*,'L2 norm of div u=',divl2
*      print*,'div u at the two corners=',u(0,0),u(n,0)
      write(30,200) ((u(i,j),j=0,m),i=0,n)

      write(10,*) 'x=['
      write (10,100) (xn(i),i=0,n)
      write(10,*) '];'

      write(10,*) 'y=['
      write (10,100) (xm(i),i=0,m)
      write(10,*) '];'
      write(10,*) 'contour(x,y,str)'

  100  format (4(1pg14.6),'...')
  200  format(50(1pg12.5,2x))
       end
