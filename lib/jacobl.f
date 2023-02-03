      subroutine jacobl(n,alpha,beta,xjac)
******************
* Computing the Gauss-Lobatto points for Jacoby polynomials
******************
*
* n: degree of approximation
*
* alpha, beta: Parameters in jacoby weight
*
* xjac: output containing Gauss-Lobatto points
*
      implicit double precision (a-h,o-z)
      dimension xjac(1)
      common /jacpar/alp,bet,rv
      data kstop/20/,one/1.0d0/,onen/-1.0d0/
      data eps/1.0d-18/
      alp=alpha
      bet=beta
      rv=1+alp
      np=n+1

      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,one)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,onen)
      det=pnp*pnm1m-pnm*pnm1p
      rp=-pnp1p
      rm=-pnp1m
      a=(rp*pnm1m-rm*pnm1p)/det
      b=(rm*pnp-rp*pnm)/det

      xjac(1)=1.0d0
      nh=(n+1)/2

      pi=dacos(-1.0d0)
      dth=pi/(2*n+1)
      cd=cos(2*dth)
      sd=sin(2*dth)
      cs=cos(dth)
      ss=sin(dth)

      do 39 j=2,nh
         x=cs
         do 29 k=1,kstop
            call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
            poly=pnp1+a*pn+b*pnm1
            pder=pdnp1+a*pdn+b*pdnm1
            recsum=0.0d0
            jm=j-1
            do 27 i=1,jm
               recsum=recsum+one/(x-xjac(i))
 27         continue
 28         continue
            delx=-poly/(pder-recsum*poly)
            x=x+delx
            if(dabs(delx).lt.eps) goto 30
 29      continue
 30      continue
         xjac(j)=x
         cssave=cs*cd-ss*sd
         ss=cs*sd+ss*cd
         cs=cssave
 39   continue
      xjac(np)=-1.0d0
      npp=n+2

      do 49 i=2,nh
         xjac(npp-i)=-xjac(i)
 49   continue
      if (n.ne.2*(n/2)) return
      xjac(nh+1)=0.0d0
      return
      end


      subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)

      implicit double precision (a-h,o-z)
      common /jacpar/alp,bet,rv
      apb=alp+bet
      poly=1.0
      pder=0.0
      if (n.eq.0) return
      polylst=poly
      pderlst=pder
      poly=rv*x
      pder=rv
      if (n.eq.1) return
      do 19 k=2,n
         a1=2*k*(k+apb)*(2*k+apb-2)
         a2=(2*k+apb-1)*(alp**2-bet**2)
         b3=2*k+apb-2
         a3=b3*(b3+1)*(b3+2)
         a4=2*(k+alp-1)*(k+bet-1)*(2.*k+apb)
         polyn=((a2+a3*x)*poly-a4*polylst)/a1
         pdern=((a2+a3*x)*pder-a4*pderlst+a3*poly)/a1
         psave=polylst
         pdsave=pderlst
         polylst=poly
         poly=polyn
         pderlst=pder
         pder=pdern
 19   continue
      polym1=polylst
      pderm1=pderlst
      polym2=psave
      pderm2=pdsave
      return
      end

