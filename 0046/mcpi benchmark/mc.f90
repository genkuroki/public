program main
  implicit none
  integer :: seed, i, n
  integer :: num_inside, num_points
  real(8) :: x, y, r2

  seed = 20231226
  num_points = 1000000000
              
  num_inside = 0
  do i = 1, num_points
    call LCGs(seed, x)
    call LCGs(seed, y)
    r2 = x**2 + y**2
    if(r2<1d0)num_inside = num_inside + 1
  end do

  write(*,*)"pi=",4*dble(num_inside)/num_points


  contains
! Linear congruential generators (LCGs)
! Parameters are provided by Park and Miller
! See https://c-faq.com/lib/rand.html
    subroutine LCGs(seed, rand_num)
      implicit none
      integer,parameter :: a = 48271
      integer,parameter :: m = 2147483647
      integer,parameter :: q = m/a
      integer,parameter :: r = mod(m,a)
      integer,intent(inout) :: seed
      real(8),intent(out) :: rand_num
      integer :: hi, lo, test

      hi = seed/q
      lo = mod(seed, q)
      test = a * lo - r * hi

      if(test>0)then
        seed = test
      else
        seed = test + m
      end if

      rand_num = dble(seed)/m

    end subroutine LCGs
end program main


