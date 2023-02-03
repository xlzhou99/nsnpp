      subroutine leinit(n,xjac,wt,a,b)
*******************
*** Initialization for Legendre-Galerkin ******
***
***  Output:    a: to be used for lestop.f 
***             b: to be used for leptos.f 
***                (only a will be used in letrfm1.f)
***             xjac:  Gauss-Lobatto points
***             wt:     weights for the discrete inner product
***
***  More precisely, a, b are defined as follows:
****  for k even,  A(k/2   ,j)=L_k(x_j), B(k/2    ,j)=L_k(x_j)*w_j/(gamma_k) 
****  for k odd,  A((k-1)/2,n/2+1+j)=L_k(x_j)
****            , B((k-1)/2,n/2+1+j)=L_k(x_j)*w_j/(gamma_k) 
******************* 
      implicit double precision (a-h, o-z)
      dimension xjac(0:n),wt(0:n),a(0:n/2,0:n),b(0:n/2,0:n)
      alpha=0.d0
      beta =0.d0
      call jacobl(n,alpha,beta,xjac)

***   compute L_0(x_j) and L_1(x_j) ***         

      do i=0,n/2
         b(i,0)=1.d0
         b(i,1)=xjac(i)
      enddo

**** compute the rest of L_k(x_j)  ---> b_jk  *****
      
      do i=0,n/2
         do j=1,n-1
            b(i,j+1)=((2*j+1.d0)*xjac(i)*b(i,j)-j*b(i,j-1))/(j+1.d0)
         enddo
      enddo
**** compute the weights ****

      a1=2.d0/(n*(n+1.d0))
      do i=0,n/2
         wt(i)=a1/(b(i,n)**2)
         wt(n-i)=wt(i)
      enddo

****** rearrange a so that:  for j even, a(j/2,k) =L_j(x_k) *****
******                       for j odd,  a((j-1)/2,n/2+1+k)= L_j(x_k)
      do k=0,n/2-1
         do j=0,n,2
            a(j/2,k)=b(k,j)
         enddo
         do j=1,n,2
            a((j-1)/2,n/2+k+1)=b(k,j)
         enddo
      enddo
      do j=0,n,2
         a(j/2,n/2)=b(n/2,j)
      enddo

**** compute: for k even,  B(k/2,j)=L_k(x_j)*w_j/(\gamma_k) ****
****          for k odd,   B((k-1)/2,n/2+1+j)=L_k(x_j)*w_j/(\gamma_k) ****
      do j=0,n/2
         do i=0,n-1,2
            b(i/2,j)=(i+.5d0)*a(i/2,j)*wt(j)
         enddo
      enddo

      do j=0,n/2
         b(n/2,j)=.5d0*n*a(n/2,j)*wt(j)
      enddo

      do j=0,n/2-1
         do i=1,n-1,2
            b((i-1)/2,n/2+1+j)=(i+.5d0)*a((i-1)/2,n/2+1+j)*wt(j)
         enddo
      enddo


      return
      end
