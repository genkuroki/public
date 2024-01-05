program sumofsinc
implicit none
integer, parameter :: dp = kind(0.d0)
real(dp) :: t1, t2, r

call cpu_time(t1)
r = f(100000000)
call cpu_time(t2)

print *, "Time", t2-t1
print *, r

contains

    real(dp) function f(N) result(r)
    integer, intent(in) :: N
    integer :: i
    r = 0
    do i = 1, N
        r = r + sin(real(i,dp)) / real(i,dp)
    end do
    r = 1 + 2 * r
    end function

end program
