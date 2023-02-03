      subroutine lenhmz0(n,f,inf,rm,alpha,a,aa,ww,init)
***************************************
* "Solve" the 1D  equation: rm*u-alpha* u_xx= f without the term phi_0 !
*  with homogeneous Neumann B.C.
*  This routine is only to be used in lenhmz2.f ! Use lenhmz1.f to solve
*  the 1-D Helmholtz equation above !
*
*  Input:  f: f_{i*inf} =(RHS, \phi_i),i=1,...,n
*       init: initialization flag: init=1, initialize a
*          a: working spaces
*  Output: f, f_{i*inf}: coefficients in \phi_i of the solution.
***************************************** 
      implicit double precision (a-h, o-z)
      dimension a(n,3),f(0:(n-1)*inf),aa(0:n),ww(0:n)
      if(init.eq.1) then

**** set up the tridiagonal matrix a ****
      do i=1,n-2
         a(i,3)=-rm*aa(i)/(i+2.5d0)*ww(i)*ww(i+2)
         a(i+2,1)=a(i,3)

      enddo
      
      do i=1,n
      a(i,2)=alpha+rm*(1.d0/(i+.5d0)+aa(i)*aa(i)/(i+2.5d0))*ww(i)*ww(i)
      enddo

**** LU decomposition *****

         do i=3,n
            a(i,1)=a(i,1)/a(i-2,2)
            a(i,2)=a(i,2)-a(i,1)*a(i-2,3)
         enddo

      endif
      
**** Ly=f ****
      do j=3,n
         f((j-1)*inf)=f((j-1)*inf)-a(j,1)*f((j-3)*inf)
      enddo

**** Ux=y *****

      f((n-1)*inf)=f((n-1)*inf)/a(n,2)
      f((n-2)*inf)=f((n-2)*inf)/a(n-1,2)

      do j=n-2,1,-1
         f((j-1)*inf)=(f((j-1)*inf)-a(j,3)*f((j+1)*inf))/a(j,2)
      enddo


      return
      end
