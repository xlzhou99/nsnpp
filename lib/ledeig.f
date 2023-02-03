      subroutine ledeig(n1,d,z,wk)
**************************************************
* Compute the eigenvalues and eigenvectors of the problem
*  "-u_xx = \lambda u,  u(-1)=u(1)=0"
*  using Legendre-Galerkin approximation 
*
*  Input:    n1: size of the problem
*            wk: working array
*  Output:    d: The eigenvalues
*             z: The eigenvectors
*
*  Rmk. The output of this routine is to be used in ledhmz2.f !
**************************************************
      implicit double precision (a-h,o-z)
      dimension d(n1-1),z(n1/2,n1-1),wk(n1/2*(n1/2+2))

      n=n1-2
      n2=n1/2
      
****** Even modes ********
      do i=1,n2
         c1=1.d0/(2*(i-1)+.5d0)
         c2=1.d0/(2*(i-1)+2.5d0)
         c3=1.d0/dsqrt(8*(i-1)+6.d0)
         d(i)=(c1+c2)*(c3*c3)
         wk(i)=-c2*c3/dsqrt(8*i+6.d0)
      enddo


      call dstev('V',n2,d,wk,z,n2,wk(n2+1),info)

******* Odd modes *******
      do i=1,n/2
         c1=1.d0/(2*i-.5d0)
         c2=1.d0/(2*i+1.5d0)
         c3=1.d0/dsqrt(8*i+2.d0)
         d(n2+i)=(c1+c2)*(c3*c3)
         wk(i)=-c2*c3/dsqrt(8*i+10.d0)
      enddo

      call dstev('V',n2-1,d(n2+1),wk,z(1,n2+1),n2,wk(n2),info)

      return
      end
