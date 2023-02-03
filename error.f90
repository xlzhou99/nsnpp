    subroutine error( n, m, ron, rom, f,g, res1, res2)
! input: f,g: physical values
    implicit double precision(a-h,o-z)
    dimension ron(0:n), rom(0:m), f(0:n,0:m), g(0:n,0:m)
    res1 = 0.
    res2 = 0.	
     !   cmoy=f(2,2)-g(2,2)
    do j = 0,m
    do i = 0,n
        x = f(i,j)-g(i,j)-cmoy
        res1 = res1 + ron(i)*rom(j)*x*x
        res2 = max(res2,abs(x))
        end do
        end do
    res1 = sqrt(res1)

    return
    end

