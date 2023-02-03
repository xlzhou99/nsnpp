      double precision function dlel2no(n,m,isp,a,na,wtn,wtm)
********************************************
*  Compute the L^2 norm of a given function
*
*   Input: isp=1, then a contains the Legendre-coefficients
*
*          isp.ne.1, then a contains the values at the L-G-L points,
*          wtn,wtm contain the weights relative to the discrete products
*          wtn and wtm can be obtained by calling leinit.f
*
*   Output: Value of the l^2 or discrete quadratic norm
*******************************************
      implicit double precision (a-h,o-z)
      dimension a(0:na,0:m),wtn(0:n),wtm(0:m)
      tmp=0.d0

      if(isp.eq.1) then

         do j=0,m
            t1=1.d0/(j+.5d0)
            do i=0,n
               tmp=tmp+a(i,j)*a(i,j)*t1/(i+.5d0)
            enddo
         enddo
         dlel2no=sqrt(tmp)

      else
         
         do j=0,m
            do i=0,n
               tmp=tmp+a(i,j)*a(i,j)*wtn(i)*wtm(j)
            enddo
         enddo
         dlel2no=sqrt(tmp)

      endif
      return
      end
