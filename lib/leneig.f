      subroutine leneig(n1,m1,d,z,aa,ww,wk)
**************************************************
**** Compute the eigenvalues and eigenvectors of 
****  -D_xx with homogeneous Neumann boundary condition
****  using Legendre-Galerkin approximation 
*
* Basis: phi_i = ww(i) (L_i - aa(i) L_{i+2}) i=0,1,..., n1-2
*        (phi'_i,phi'_j)= Delta_{ij} for i,j >= 1 !
* Input: n1: size of the program
*        wk: working array of dimension >= 3*n1
* Output: d: eigenvalues
*         z: eigenvectors
*        aa,ww: as described above of size max(n1,m1)+1
*
* The output of this routine is to be used in lenhmz2.f !
**************************************************
      implicit double precision (a-h,o-z)
      dimension d(n1-2),z((n1-2)/2,n1-2),wk(1),aa(0:0),ww(0:0)

      n=n1-2
      n2=n/2
      do 1 i=1,max(n1,m1)
         aa(i)=i*(i+1.d0)/((i+2.d0)*(i+3.d0))
 1    ww(i)=1.d0/dsqrt(aa(i)*(4*i+6.d0))
      aa(0)=0.d0
      ww(0)=1.d0
****** Even modes ********
      do 10 i=1,n2
         c1=1.d0/(2*i+.5d0)
         c2=1.d0/(2*i+2.5d0)
         d(i)=(c1+c2*aa(2*i)**2)*ww(2*i)*ww(2*i)
 10   wk(i)=-c2*aa(2*i)*ww(2*i)*ww(2*i+2)
      call dstev('V',n2,d,wk,z,n2,wk(n2+1),info)

****** Odd modes *******
      do 20 i=1,n2
         c1=1.d0/(2*i-.5d0)
         c2=1.d0/(2*i+1.5d0)
         d(n2+i)=(c1+c2*aa(2*i-1)**2)*ww(2*i-1)**2
 20   wk(i)=-c2*aa(2*i-1)*ww(2*i-1)*ww(2*i+1)

      call dstev('V',n2,d(n2+1),wk,z(1,n2+1),n2,wk(n2+1),info)

      return
      end

