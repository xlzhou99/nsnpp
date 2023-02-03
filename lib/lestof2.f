      subroutine lestof2(n,m,s,ns,ss,nss)
***********************************
*** Purpose: compute ss(i,j)=(s,\phi_i(x)\phi_j(y))
***          where \phi_i(x) = L_i(x) - L_{i+2}(x).
***
*** Input:    n,m: size of the problem
***           s: coefficients of the Legendre expansion
*** Output:
***           ss: ss(i,j)=(s,\phi_i(x)\phi_j(y))., nss>=n-2
***********************************
      implicit double precision (a-h, o-z)
      dimension s(0:ns,0:m),ss(0:nss,0:(m-2))

      do j=0,m-2
         b1=1.d0/(j+.5d0)
         b2=1.d0/(j+2.5d0)
         do  i=0,n-2
            a1=1.d0/(i+.5d0)
            a2=1.d0/(i+2.5d0)
            ss(i,j)=a1*(s(i,j)*b1-s(i,j+2)*b2)
     1             +a2*(s(i+2,j+2)*b2-s(i+2,j)*b1)
         enddo
      enddo
      
      return
      end
