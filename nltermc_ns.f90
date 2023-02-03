  subroutine nltermc_ns(n,m,n2,m2,h,u,v,ux,vx,dn,dm)
!*****************************************************
!* computer the non-linear term of the Navier-Stokes equation
!*  Nlterm:   u u_x +  v u_y/h,   u v_x +  v v_y/h
!*        input:  u,v:  physical value, nochange on output.
!*        output: ux,vx: physical value of non-linear term  
!*
!*****************************************************
  implicit double precision (a-h, o-z)

     real (kind=8),dimension(0:n2)       :: wtn
       real (kind=8),dimension(0:m2)       :: wtm
       real (kind=8),dimension(0:n2,0:n2)  :: dn
       real (kind=8),dimension(0:m2,0:m2)  :: dm
       real (kind=8),dimension(0:n2,0:m2)  ::u,v,ux,vx
  real (kind=8),allocatable :: wk1(:,:),wk2(:,:)

  allocate(wk1(0:n2,0:m2),wk2(0:n2,0:m2))
  call ledxc(n,m,u,n2,wk1,n2,dn)
  call ledyc(n,m,u,n2,wk2,n2,dm)

    do j=0,m
      do i=0,n
        ux(i,j)=u(i,j)*wk1(i,j)+v(i,j)*wk2(i,j)/h
      enddo
    enddo
    call ledxc(n,m,v,n2,wk1,n2,dn)
  call ledyc(n,m,v,n2,wk2,n2,dm)
    do j=0,m
      do i=0,n
        vx(i,j)=u(i,j)*wk1(i,j)+v(i,j)*wk2(i,j)/h
      enddo
    enddo
    deallocate(wk1,wk2)
    return
    end subroutine nltermc_ns