      subroutine lestop2(n,m,s,ns,f,an,am)
*************************************
*** 2D Legendre transfrom  ****
*** Input: an(L_k(x_j)), am(L_m(y_j)): output from leinit.f
*** Input : s (spectral value) 
***        an,am: output from leinit.f
***         f: working array of dimension >=(n+1)*(m+1)
*** output: s (physical value)
***
*** Rmk. Only correct if n and m are even !
*************************************
      implicit double precision (a-h,o-z)
      dimension s(0:ns,0:m),an(0:n/2,0:n),am(0:m/2,0:m),f(0:n,0:m)
      data one/1.d0/,zero/0.d0/

      if (mod(n,2)+mod(m,2).ne.0) then
         print*,'n or m is not even !'
         stop
      endif

      do i=0,m
         do j=0,n/2-1
            f(j,i)=s(2*j,i)
            f(n/2+1+j,i)=s(2*j+1,i)
         enddo
         f(n/2,i)=s(n,i)
      enddo


      call  dgemm('t','n',n/2+1,m+1,n/2+1,one,an,n/2+1,f,n+1
     1            ,zero,s,ns+1)

      call  dgemm('t','n',n/2,m+1,n/2,one,an(0,n/2+1),n/2+1,f(n/2+1,0)
     1            ,n+1,zero,s(n/2+1,0),ns+1)

****** We used the fact that L_j(x(n/2))=0 for j odd !

      do i=0,m
         do j=0,n/2-1
            f(j,i)=s(j,i)+s(n/2+1+j,i)
            f(n-j,i)=s(j,i)-s(n/2+1+j,i)
         enddo
         f(n/2,i)=s(n/2,i)
      enddo



      do i=0,m/2-1
         do j=0,n
            s(j,i)=f(j,2*i)
            s(j,m/2+1+i)=f(j,2*i+1)
         enddo
      enddo
      do j=0,n
         s(j,m/2)=f(j,m)
      enddo

      call  dgemm('n','n',n+1,m/2+1,m/2+1,one,s,ns+1,am
     1            ,m/2+1,zero,f,n+1)

      call  dgemm('n','n',n+1,m/2,m/2,one,s(0,m/2+1),ns+1,am(0,m/2+1)
     1            ,m/2+1,zero,f(0,m/2+1),n+1)
      
      do i=0,m/2-1
         do j=0,n
              s(j,i)=f(j,i)+f(j,m/2+i+1)
            s(j,m-i)=f(j,i)-f(j,m/2+i+1)
         enddo
      enddo
      do j=0,n
         s(j,m/2)=f(j,m/2)
      enddo

      return 
      end
