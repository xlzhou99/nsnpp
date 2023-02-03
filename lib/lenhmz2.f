      subroutine lenhmz2(n1,m1,rm,beta,f,nf,u,nu,e,w,aa,ww,ws,init)
***************************************
* Solve the 2D Helmholtz equation: rm*u - u_xx - beta u_yy= f with
*     homogeneous Neumann B.C.
*
*  Input:  f: f(0:nf, 0:m1) coefficients in Legendre of right hand side f
*          nf: nf>=n1
*        e,w: Output from leneig.f, contain eigenfunctions and eigenvalues
*             of -u_{xx}= \lambda u.
*         ws: working array. There are initialized with init=1, they can be
*             repeatedly used with init.ne.1.
*         rm: constant of the Helmholtz equation
*       init: initialization flag: init=1, initialize ws.
*
*  Output: u, coefficients in Legendre of the solution. (nu >= n1)
*
* Remark: Most efficient when m1 >= n1 !  
*
* Basis: phi_i = ww(i) (L_i - aa(i) L_{i+2}) i=0,1,..., n1-2
*        (phi'_i,phi'_j)= Delta_{ij} for i,j >= 1, ww(0)=1.
*
***************************************
      implicit double precision (a-h,o-z)
      dimension e(n1/2-1,n1-2),w(n1-2),ws(3*(m1-2),0:n1-2),u(0:nu,0:m1)
      dimension f(0:nf,0:m1),aa(0:0),ww(0:0)
      n=n1-2
      m=m1-2
      n2=n/2

      call lestof2b(n1,m1,f,nf,u,nu,aa,ww)

***  Take care of the components corresponding to first row and first column
      if (rm.eq.0.d0) then
         u(0,0)=0.d0
         do j=1,m
            u(0,j)=u(0,j)/(2*beta)
         enddo
         do i=1,n
            u(i,0)=.5d0*u(i,0)
         enddo
      else
         u(0,0)=u(0,0)*.25d0/rm
         call lenhmz0(m,u(0,1),nu+1,2*rm,2*beta,ws(1,0),aa,ww,init)
         call lenhmz0(n,u(1,0),1,2*rm,2.d0,f,aa,ww,1)
      endif

*** Store the first column and row in the last column and row of F !
      do i=1,n
         f(i-1,m)=u(i,0)
      enddo
      do j=1,m
         f(n,j-1)=u(0,j)
      enddo
      f(n,m)=u(0,0)

***F=E^t*U *** 

      do j=1,m
         do i=2,n,2
            f(i/2-1,j-1)=u(i,j)
         enddo
         do i=1,n,2
            f(n2+(i-1)/2,j-1)=u(i,j)
         enddo
      enddo
      call dgemm('t','n',n2,m,n2,1.d0,e,n2,f,nf+1,0.d0,u,nu+1)
      do j=0,m-1
         do i=0,n2-1
            f(i,j)=f(n2+i,j)
            f(n2+i,j)=u(i,j)
         enddo
      enddo
      call dgemm('t','n',n2,m,n2,1.d0,e(1,n2+1),n2,f,nf+1
     1          ,0.d0,u,nu+1)

      do j=0,m-1
         do i=0,n2-1
            f(i,j)=f(n2+i,j)
            f(n2+i,j)=u(i,j)
         enddo
      enddo

****
      do i=1,n
         rm1=1+rm*w(i)
      call lenhmz0(m,f(i-1,0),nf+1,rm1,w(i)*beta,ws(1,i),aa,ww,init)
      enddo


***** U=E*F *****
      call dgemm('n','n',n2,m,n2,1.d0,e,n2,f,nf+1,0.d0,u,nu+1)
      do j=0,m-1
         do i=0,n2-1
            tmp=f(n2+i,j)
            f(n2+i,j)=u(i,j)
            u(i,j)=tmp
         enddo
      enddo
      call dgemm('n','n',n2,m,n2,1.d0,e(1,n2+1),n2,u,nu+1
     1           ,0.d0,f,nf+1)

      do j=0,m-1
         do i=0,n2-1
            u(2*i+1,j)=f(n2+i,j)
            u(2*i,j)=f(i,j)
         enddo
      enddo

**** rearrange u to its original position
      do j=m-1,0,-1
         do i=n-1,0,-1
            u(i+1,j+1)=u(i,j)
         enddo
      enddo
      do i=1,n
         u(i,0)=f(i-1,m)
      enddo
      do j=1,m
         u(0,j)=f(n,j-1)
      enddo
      u(0,0)=f(n,m)

      call ftos2b(n1,m1,u,nu,f,aa,ww)

      return
      end

