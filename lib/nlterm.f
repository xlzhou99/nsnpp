      subroutine nlterm(n,m,h,n2,m2,u,v,ux,vx,uy,vy,an,am,bn,bm)
*****************************************************
* computer the non-linear term of the Navier-Stokes equation
*  Nlterm:   u u_x +  v u_y/h,   u v_x +  v v_y/h
*        input:  u,v:  spectral value, destroyed on output.
*        output: ux,vx: spectral value of non-linear term  
*
*****************************************************
      implicit double precision (a-h,o-z)
      dimension u(0:n2,0:m2),v(0:n2,0:m2),ux(0:n2,0:m2),vx(0:n2,0:m2)
      dimension uy(0:n2,0:m2),vy(0:n2,0:m2)
      dimension an(0:n,0:n),am(0:m,0:m),bn(0:n,0:n),bm(0:m,0:m)


      call ledx(n,m,u,n2,ux,n2)
      call ledy(n,m,u,n2,uy,n2)

      call lestop2(n,m,u, n2,vx,an,am)
      call lestop2(n,m,ux,n2,vx,an,am)
      call lestop2(n,m,uy,n2,vx,an,am)

      call ledx(n,m,v,n2,vx,n2)
      call lestop2(n,m,vx,n2,vy,an,am)
      

      do j=0,m
         do i=0,n
            ux(i,j)=u(i,j)*ux(i,j)
            vx(i,j)=u(i,j)*vx(i,j)
         enddo
      enddo

      call ledy(n,m,v,n2,vy,n2)
      call lestop2(n,m,v, n2,u,an,am)
      call lestop2(n,m,vy,n2,u,an,am)

      do j=0,m
         do i=0,n
            ux(i,j)=ux(i,j)+v(i,j)*uy(i,j)/h
            vx(i,j)=vx(i,j)+v(i,j)*vy(i,j)/h
         enddo
      enddo

      call leptos2(n,m,ux,n2,uy,bn,bm)
      call leptos2(n,m,vx,n2,uy,bn,bm)

      return
      end
