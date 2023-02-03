      subroutine ledhmz1(n,f,inf,rm,alpha,a,init)
***************************************
* Solve the 1D Helmholtz equation: rm*u-alpha* u_xx= f!
*  with homogeneous Dirichlet B.C.
*  Input:  f: f_{i*inf} =(RHS, \phi_i),i=0,1,...,n
*       init: initialization flag: init=1, initialize a
*          a: working spaces
*  Output: f, f_{i*inf}: coefficients in \phi_i of the solution.
***************************************** 
      implicit double precision (a-h, o-z)
      dimension a(0:n,3),f(0:0)
      if(init.eq.1) then

**** set up the a *****
         do i=0,n
         a(i,1)=-rm/(i+.5d0)
         a(i,3)=-rm/(i+2.5d0)
         a(i,2)=alpha*(4*i+6)-a(i,1)-a(i,3)
      enddo

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

      return
      end



