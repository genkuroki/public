program main
  implicit none
  real(8) :: x, v
  real(8),allocatable :: xt(:), vt(:)
  real(8) :: mass, k, dt
  integer :: it, nt

  mass = 1d0
  k    = 1d0
  dt   = 1d-2
  nt   = 100000000

  allocate(xt(0:nt), vt(0:nt))

  x = 0d0
  v = 1d0


  do it = 0, nt
    xt(it) = x
    vt(it) = v

    call Runge_Kutta_4th(x,v,dt,mass,k)

  end do

  open(20,file="result_fortran.out")
  do it = nt-1000, nt
    write(20,"(3e26.16e3)")it*dt, xt(it), vt(it)
  end do
  close(20)


  contains 

    subroutine Runge_Kutta_4th(x,v,dt,mass,k)
      implicit none
      real(8),intent(inout) :: x, v
      real(8),intent(in) :: dt, mass, k
      real(8) :: x1,x2,x3,x4,v1,v2,v3,v4

! RK1
      x1 = v
      v1 = force(x, mass, k)

! RK2
      x2 = v+0.5d0*dt*v1
      v2 = force(x+0.5d0*x1*dt, mass, k)

! RK3
      x3 = v+0.5d0*dt*v2
      v3 = force(x+0.5d0*x2*dt, mass, k)

! RK4
      x4 = v+dt*v3
      v4 = force(x+x3*dt, mass, k)

      x = x + (x1+2d0*x2+2d0*x3+x4)*dt/6d0
      v = v + (v1+2d0*v2+2d0*v3+v4)*dt/6d0

    end subroutine Runge_Kutta_4th

    real(8) function force(x,mass,k)
      implicit none
      real(8),intent(in) :: x, mass, k

      force = -x*k/mass

    end function force

end program main
