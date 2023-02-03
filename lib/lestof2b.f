      subroutine lestof2b(n,m,s,ns,ss,nss,aa,ww)
***********************************
*   For Neumann boundary conditions !
*
*** computing ss=(s,\phi_i(x)\phi_j(y)) with phi_i=ww(i)(l_i-aa(i)L_{i+2))
*** Input:
***           s: spectral value in L (the righthand side f)
*** Output:
***           ss: Spectral value in \phi
***********************************
      implicit double precision (a-h, o-z)
      dimension s(0:ns,0:m),ss(0:nss,0:m),aa(0:0),ww(0:0)

      do i=0,n-2
         a1=1.d0/(i+.5d0)
         a2=1.d0/(i+2.5d0)*aa(i)
         do j=0,m-2
            b1=1.d0/(j+.5d0)
            b2=1.d0/(j+2.5d0)*aa(j)
            ss(i,j)=(a1*(s(i,j)*b1-s(i,j+2)*b2)
     1              +a2*(s(i+2,j+2)*b2-s(i+2,j)*b1))*(ww(i)*ww(j))
         enddo
      enddo
      
      return
      end
