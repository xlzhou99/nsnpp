      double precision function dlepsw(n,m,isp,a,b,na,wtn,wtm)
!********************************************
!*  Compute the L^2 inner product of two given function
!*
!*   Input: isp=1, then a contains the Legendre-coefficients
!*
!*          isp.ne.1, then a contains the values at the L-G-L points,
!*          wtn,wtm contain the weights relative to the discrete products
!*          wtn and wtm can be obtained by calling leinit.f
!*
!* Rmk: 1-D case can be called with m=0 !
!*
!*   Output: Value of the l^2 or discrete quadratic product
!*******************************************

      implicit double precision (a-h,o-z)
      dimension a(0:na,0:m),b(0:na,0:m),wtn(0:n),wtm(0:m)
      tmp=0.d0

      if (m.eq.0) then

      if(isp.eq.1) then
         do i=0,n
            tmp=tmp+a(i,0)*b(i,0)/(i+.5d0)
         enddo
         dlepsw=tmp

      else
         
         do i=0,n
            tmp=tmp+a(i,0)*b(i,0)*wtn(i)
         enddo
         dlepsw=tmp

      endif

      return

      endif

      if(isp.eq.1) then

         do j=0,m
            t1=1.d0/(j+.5d0)
            do i=0,n
               tmp=tmp+a(i,j)*b(i,j)*t1/(i+.5d0)
            enddo
         enddo
         dlepsw=tmp

      else
         
         do j=0,m
            do i=0,n
               tmp=tmp+a(i,j)*b(i,j)*wtn(i)*wtm(j)
            enddo
         enddo
         dlepsw=tmp
      endif
      return
      end
