       subroutine ledhmz2(n1,m1,rm,beta,f,nf,u,nu,e,w,ws,init)
***************************************
* Solve the 2D Helmholtz equation: rm*u - u_xx - beta u_yy= f
* with homogeneous Dirichlet B.C. in IxI.
*
*  Input:  f: f(0:nf, 0:m1) coefficients in Legendre of right hand side f
*         nf: nf>=n1
*        e,w: Output from ledeig.f, contain eigenfunctions and eigenvalues
*             of -u_{xx}= \lambda u.
*         ws: working array. There are initialized with init=1, they can be
*             repeatedly used with init.ne.1.
*         rm,beta: constants in the Helmholtz equation
*       init: initialization flag: init=1, initialize ws.
*
*  Output: u, coefficients in Legendre of the solution. (nu >= n1)
*
* Remark: Most efficient when m1 >= n1 !  
***************************************
      implicit double precision (a-h,o-z)
      dimension e(n1/2,n1-1),w(0:n1-2)
      dimension ws(3*(m1-1),0:n1-2),u(0:nu,0:m1),f(0:nf,0:m1)
      n=n1-2
      m=m1-2
      n2=n/2+1

      call lestof2(n1,m1,f,nf,u,nu)

**** u=A_1^{-1/2} u ***
      do i=0,n
         f(i,0)=1.d0/dsqrt(4*i+6.d0)
      enddo
      do j=0,m
         do i=0,n
            u(i,j)=f(i,0)*u(i,j)
         enddo
      enddo

***F=E^t*U ***
      do j=0,m
         do i=0,n,2
            f(i/2,j)=u(i,j)
         enddo
         do i=1,n,2
            f(n2+(i-1)/2,j)=u(i,j)
         enddo
      enddo
      call dgemm('t','n',n2,m+1,n2,1.d0,e,n2,f,nf+1,0.d0,u,nu+1)
      do j=0,m
         do i=0,n/2-1
            f(i,j)=f(n2+i,j)
            f(n/2+i,j)=u(i,j)
         enddo
         f(n,j)=u(n/2,j)
      enddo
      call dgemm('t','n',n/2,m+1,n/2,1.d0,e(1,n2+1),n2,f,nf+1
     1          ,0.d0,u,nu+1)

      do j=0,m
         do i=0,n/2
            f(i,j)=f(n/2+i,j)
         enddo
         do i=1,n/2
            f(n2-1+i,j)=u(i-1,j)
         enddo
      enddo

      do i=0,n
         rm1=1+rm*w(i)
         call ledhmz1(m,f(i,0),nf+1,rm1,w(i)*beta,ws(1,i),init)
      enddo

***** U=E*F *****
      call dgemm('n','n',n2,m+1,n2,1.d0,e,n2,f,nf+1,0.d0,u,nu+1)
      do j=0,m
         f(n/2,j)=u(0,j)
         do i=1,n/2
            tmp=f(n/2+i,j)
            f(n/2+i,j)=u(i,j)
            u(i-1,j)=tmp
         enddo
      enddo
      call dgemm('n','n',n/2,m+1,n/2,1.d0,e(1,n2+1),n2,u,nu+1
     1           ,0.d0,f,nf+1)

      do j=0,m
         do i=1,n2
            u(2*(i-1),j)=f(n/2+i-1,j)
         enddo
         do i=1,n/2
            u(2*i-1,j)=f(i-1,j)
         enddo
      enddo

**** u=A_1^{-1/2} u ***
      do i=0,n
         f(i,0)=1.d0/dsqrt(4*i+6.d0)
      enddo
      do j=0,m
         do i=0,n
            u(i,j)=f(i,0)*u(i,j)
         enddo
      enddo
 
      call ftos2(n1,m1,u,nu,f)
      return
      end

