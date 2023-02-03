      subroutine ftos2(n,m,ss,nss,s)
******Transfrom from \phi to L *******
*  Purpose: compute the coefficients of the Legendre or Chebyshev
*   expansion from the coefficients of expansion in  \Phi_k(x)\phi_j(y),
*           where \phi_k (x)= L_k(x) - L_{k+2} (x) ,
*             or  \phi_k (x)= T_k(x) - T_{k+2} (x) . 
*
* Input:  n,m:  size of the problem
*          ss: spectral value in \Phi_k(x)\phi_j(y), nss>=n
*           s: working array of dimension >=(n+1)*(m-1)
* Output:  ss: spectral value in L_k(x)L_j(y) or T_k(x)T_j(y)
**************************************
      implicit double precision (a-h, o-z)
      dimension s(0:n,0:m-2),ss(0:nss,0:m)
      do i=0,m-2
         s(0,i)=ss(0,i)
         s(1,i)=ss(1,i)
         s(n-1,i)=-ss(n-3,i)
         s(n,i)=-ss(n-2,i)
         do j=2,n-2
            s(j,i)=ss(j,i)-ss(j-2,i)
         enddo
      enddo

      do i=0,n
         ss(i,0)=s(i,0)
         ss(i,1)=s(i,1)
         ss(i,m-1)=-s(i,m-3)
         ss(i,m)=-s(i,m-2)
         do j=2,m-2
            ss(i,j)=s(i,j)-s(i,j-2)
         enddo
      enddo

      return
      end
