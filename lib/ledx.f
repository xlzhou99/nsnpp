      subroutine ledx(n,m,u,nu,v,nv)
*************************
* compute the first derivative w.r.t. the first variable for 1-D and 2-D
* Legendre series.
*     Input: u, spectral values in legendre series, unchanged at output
*    output: v=u_x, spectral values in legendre series
*
* Rmk: for 1-D problem, call the subroutine with m=0!
*************************  
      implicit double precision (a-h, o-z)
      dimension u(0:nu,0:m),v(0:nv,0:m)
      do j=0,m
         v(n,j)=0.d0
         v(n-1,j)=(2*n-1)*u(n,j)
         do i=n-2,0,-1
            v(i,j)=(2*i+1)*(v(i+2,j)/(2*i+5)+u(i+1,j))
         enddo
      enddo
      return
      end
