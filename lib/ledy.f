      subroutine ledy(n,m,u,nu,v,nv)
*************************
* compute the first derivative w.r.t. the second variable for 1-D and 2-D
* Legendre series.
*     Input: u, spectral values in legendre series, unchanged at output
*    output: v=u_y, spectral values in legendre series
*
*************************  
      implicit double precision (a-h, o-z)
      dimension u(0:nu,0:m),v(0:nv,0:m)
      do i=0,n
         v(i,m)=0.d0
         v(i,m-1)=(2*m-1)*u(i,m)
      enddo
      do j=m-2,0,-1
         do i=0,n
            v(i,j)=(2*j+1)*(v(i,j+2)/(2*j+5)+u(i,j+1))
         enddo
      enddo
      return
      end
