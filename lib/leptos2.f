      subroutine leptos2(n,m,f,nf,s,bn,bm)
*************************************
*** 2D Legendre trbnsfrom  ****
*** Input: bn,bm (L_k(x_j)*w_j/(\gbmma_k)): output from leinit.f
*** Input:  f (physical value) 
***        bn,bm: output from leinit.f
***         s: working array of dimension >=(n+1)*(m+1)
*** output: f (spectral value)
*
* Rmk. Only correct if n and m are even !
*************************************
      implicit double precision (a-h,o-z)
      dimension f(0:nf,0:m),s(0:n,0:m),bn(0:n/2,0:n),bm(0:m/2,0:m)
      data one/1.d0/,zero/0.d0/

      if (mod(n,2)+mod(m,2).ne.0) then
         print*,'n or m is not even !'
         stop
      endif

      do i=0,m
         do j=0,n/2-1
            s(j,i)=f(j,i)+f(n-j,i)
            s(n/2+j+1,i)=f(j,i)-f(n-j,i)
         enddo
         s(n/2,i)=f(n/2,i)
      enddo

****** We used the fact that L_j(x(n/2))=0 for j odd !


      call  dgemm('n','n',n/2+1,m+1,n/2+1,one,bn,n/2+1,s,n+1
     1            ,zero,f,nf+1)

      call  dgemm('n','n',n/2,m+1,n/2,one,bn(0,n/2+1),n/2+1,s(n/2+1,0)
     1            ,n+1,zero,f(n/2+1,0),nf+1)

      do j=0,m

      do i=0,n/2-1
         s(2*i,j)=f(i,j)
         s(2*i+1,j)=f(n/2+i+1,j)
      enddo
      s(n,j)=f(n/2,j)

      enddo

      do j=0,m/2-1
         do i=0,n
            f(i,j)=s(i,j)+s(i,m-j)
            f(i,m/2+j+1)=s(i,j)-s(i,m-j)
         enddo
      enddo
      do i=0,n
         f(i,m/2)=s(i,m/2)
      enddo

      call  dgemm('n','t',n+1,m/2+1,m/2+1,one,f,nf+1,bm
     1            ,m/2+1,zero,s,n+1)

      call  dgemm('n','t',n+1,m/2,m/2,one,f(0,m/2+1),nf+1,bm(0,m/2+1)
     1            ,m/2+1,zero,s(0,m/2+1),n+1)

      do i=0,m/2-1
         do j=0,n
            f(j,2*i)=s(j,i)
            f(j,2*i+1)=s(j,m/2+i+1)
         enddo
      enddo
      do j=0,n
         f(j,m)=s(j,m/2)
      enddo

      return
      end
