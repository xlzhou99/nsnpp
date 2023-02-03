	module Input_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		implicit none
        integer, parameter      :: n        = 64
        integer, parameter      :: m        = n, mn=m*n       
        integer, parameter      :: n2       = n, nn=(n2+1)*(n2+1)
        integer, parameter      :: m2       = m, mm=(m2+1)*(m2+1), nm=(n2+1)*(m2+1)
        integer, parameter      :: iord     = 2!index of order

        real (kind=8),parameter :: re       = 1.d0

        real (kind=8),parameter :: dt       = 0.0001d0 
        real (kind=8),parameter :: tmax     = 1.0d0
        real (kind=8),parameter :: T_ex     = tmax
        real (kind=8),parameter :: h        = 1.d0
	    real (kind=8),parameter :: pi        = dacos(-1.0d0)
        integer, parameter      :: irestart = 0

        real (kind=8),parameter :: c0       = 10.d0


        character(len=99)       :: fname1   = "re400"
        integer, parameter      :: ifinal   = 1
    end module Input_var
