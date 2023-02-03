      subroutine lenhmz1(n,f,inf,rm,alpha,a,aa,ww,init)
***************************************
* Solve the 1D Helmholtz equation: rm*u-alpha* u_xx= f!
*  with homogeneous Neumann B.C. 
*
*  Input:  f: f_{i*inf} =(RHS, \phi_i),i=0,1,...,n 
*       init: initialization flag: init=1, initialize a
*          a: working spaces
*  Output: f, f_{i*inf}: coefficients in \phi_i of the solution.
***************************************** 
      implicit double precision (a-h, o-z)
      dimension a(0:n,3),f(0:n*inf),aa(0:n),ww(0:n)

      if (rm.eq.0.d0) then

         f(0)=0.d0
         do j=1,n
            f(j*inf)=f(j*inf)/alpha
         enddo

      else

      if(init.eq.1) then

**** set up the tridiagonal matrix a ****
      do i=0,n-2
         a(i,3)=-rm*aa(i)/(i+2.5d0)*ww(i)*ww(i+2)
         a(i+2,1)=a(i,3)
      enddo
      
      do i=1,n
      a(i,2)=alpha+rm*(1.d0/(i+.5d0)+aa(i)*aa(i)/(i+2.5d0))*ww(i)*ww(i)
      enddo
      a(0,2)=rm*2

**** LU decomposition *****

         do i=2,n
            a(i,1)=a(i,1)/a(i-2,2)
            a(i,2)=a(i,2)-a(i,1)*a(i-2,3)
         enddo

      endif


**** Ly=f ****
      do j=2,n
         f(j*inf)=f(j*inf)-a(j,1)*f((j-2)*inf)
      enddo

**** Ux=y *****

      f(n*inf)=f(n*inf)/a(n,2)
      f((n-1)*inf)=f((n-1)*inf)/a(n-1,2)

      do j=n-2,0,-1
         f(j*inf)=(f(j*inf)-a(j,3)*f((j+2)*inf))/a(j,2)
      enddo

      endif

      return
      end
