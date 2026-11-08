!=======================================================================
! Chapter 4 -- SPHERICAL BESSEL FUNCTIONS and the SPHERICAL WELL (F2008)
! Part A: j_l by Miller downward recurrence (normalised by j_0 or j_1),
!         y_l by the stable upward recurrence, small-x series, zeros;
!         checks: cross product x^2 (j_l y_{l-1} - j_{l-1} y_l) = 1 and
!         j_l(x) = sqrt(pi/2x) J_{l+1/2}(x) for l = 0 (sin x / x)
! Part B: u = r R,  u'' + [2E - l(l+1)/r^2] u = 0,  u(1) = 0
!         Numerov + bisection; naive start vs Frobenius values
! Output: screen tables + ch4_f_conv.dat
! Compile: gfortran -O2 -o ch4 ch4_spherical_bessel_well.f90 && ./ch4
!=======================================================================
module sph_mod
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
contains
  real(dp) function j_series(l, x)
    integer,  intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: df, t
    integer  :: k
    df = 1.0_dp
    do k = 1, 2*l+1, 2
       df = df*k
    end do
    j_series = 0.0_dp; t = 1.0_dp
    do k = 0, 29
       j_series = j_series + t
       t = -t*x*x/(2*(k + 1)*(2*l + 2*k + 3))
    end do
    j_series = x**l/df*j_series
  end function j_series

  subroutine j_miller(L, x, jv)
    integer,  intent(in)  :: L
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: jv(0:L)
    real(dp), allocatable :: f(:)
    real(dp) :: j0, j1, scale
    integer  :: M, k
    if (x < 1.0d-3) then
       do k = 0, L
          jv(k) = j_series(k, x)
       end do
       return
    end if
    M = L + int(x) + 20 + int(sqrt(40*(L + x)))
    allocate(f(0:M+1)); f = 0.0_dp; f(M) = 1.0d-30
    do k = M, 1, -1
       f(k-1) = (2*k + 1)/x*f(k) - f(k+1)
       if (abs(f(k-1)) > 1.0d250) f(k-1:M+1) = f(k-1:M+1)*1.0d-250
    end do
    j0 = sin(x)/x; j1 = sin(x)/x**2 - cos(x)/x
    if (abs(j0) > abs(j1)) then
       scale = j0/f(0)
    else
       scale = j1/f(1)
    end if
    jv = f(0:L)*scale
  end subroutine j_miller

  subroutine y_up(L, x, yv)
    integer,  intent(in)  :: L
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: yv(0:L)
    integer :: k
    yv(0) = -cos(x)/x
    if (L >= 1) yv(1) = -cos(x)/x**2 - sin(x)/x
    do k = 1, L-1
       yv(k+1) = (2*k + 1)/x*yv(k) - yv(k-1)
    end do
  end subroutine y_up

  subroutine j_up(L, x, jv)
    integer,  intent(in)  :: L
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: jv(0:L)
    integer :: k
    jv(0) = sin(x)/x
    if (L >= 1) jv(1) = sin(x)/x**2 - cos(x)/x
    do k = 1, L-1
       jv(k+1) = (2*k + 1)/x*jv(k) - jv(k-1)
    end do
  end subroutine j_up

  real(dp) function jl(l, x)
    integer,  intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: jv(0:l+1)
    call j_miller(l+1, x, jv); jl = jv(l)
  end function jl

  subroutine jl_zeros(l, S, z)
    integer,  intent(in)  :: l, S
    real(dp), intent(out) :: z(S)
    real(dp) :: x, x2, prev, cur, zz, dz, d, jv(0:l+1)
    integer  :: found, it
    found = 0
    x = merge(l + 0.5_dp, 1.0_dp, l > 0); prev = jl(l, x)
    do while (found < S)
       x2 = x + 0.1_dp; cur = jl(l, x2)
       if (sign(1.0_dp, cur) /= sign(1.0_dp, prev)) then
          zz = 0.5_dp*(x + x2)
          do it = 1, 50
             call j_miller(l+1, zz, jv)
             if (l > 0) then
                d = jv(l-1) - (l + 1)/zz*jv(l)
             else
                d = -jv(1)
             end if
             dz = jv(l)/d; zz = zz - dz
             if (abs(dz) < 1.0d-15*zz) exit
          end do
          found = found + 1; z(found) = zz
       end if
       x = x2; prev = cur
    end do
  end subroutine jl_zeros
end module sph_mod

module swell_mod
  use sph_mod
  implicit none
contains
  real(dp) function frob(r, l, E)
    real(dp), intent(in) :: r, E
    integer,  intent(in) :: l
    real(dp) :: a, s, p
    integer  :: k
    s = 0.0_dp; a = 1.0_dp; p = 1.0_dp
    do k = 0, 39
       s = s + a*p
       a = -2*E*a/((2*k + 2)*(2*k + 2*l + 3)); p = p*r*r
    end do
    frob = r**(l+1)*s
  end function frob

  real(dp) function simpson(y, h, n)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    simpson = h/3*(y(0) + y(n) + 4*sum(y(1:n-1:2)) + 2*sum(y(2:n-2:2)))
  end function simpson

  subroutine shoot(E, l, N, naive, u)
    real(dp), intent(in)  :: E
    integer,  intent(in)  :: l, N
    logical,  intent(in)  :: naive
    real(dp), intent(out) :: u(0:N)
    real(dp) :: h, r(0:N), c(0:N)
    integer  :: j, js
    h = 1.0_dp/N; r = [(j*h, j=0,N)]
    c(0) = 1.0_dp; c(1:N) = 1 + h*h*(2*E - l*(l+1)/r(1:N)**2)/12
    u = 0.0_dp
    if (naive) then
       js = 0; u(1) = h**(l+1)
    else
       js = 1; u(0) = 0.0_dp; u(1) = frob(r(1), l, E); u(2) = frob(r(2), l, E)
    end if
    do j = js+1, N-1
       u(j+1) = ((12 - 10*c(j))*u(j) - c(j-1)*u(j-1))/c(j+1)
    end do
  end subroutine shoot

  subroutine eigen(l, Eex, N, naive, E, u)
    integer,  intent(in)  :: l, N
    real(dp), intent(in)  :: Eex
    logical,  intent(in)  :: naive
    real(dp), intent(out) :: E, u(0:N)
    real(dp) :: a, b, fa, mid
    integer  :: it
    a = 0.85_dp*Eex; b = 1.15_dp*Eex
    call shoot(a, l, N, naive, u); fa = u(N)
    do it = 1, 64
       mid = 0.5_dp*(a + b)
       call shoot(mid, l, N, naive, u)
       if (sign(1.0_dp, u(N)) == sign(1.0_dp, fa)) then
          a = mid; fa = u(N)
       else
          b = mid
       end if
    end do
    E = 0.5_dp*(a + b)
    call shoot(E, l, N, naive, u)
    u = u/sqrt(simpson(u*u, 1.0_dp/N, N))
  end subroutine eigen
end module swell_mod

program ch4
  use sph_mod
  use swell_mod
  implicit none
  integer, parameter :: NN = 400
  real(dp) :: jv(0:30), ju(0:30), yv(0:30), z(0:3,4), x, cp, E, u(0:NN), uex(0:NN), Eex, e1(7), e2(7), r
  real(dp) :: xc(3) = [0.5_dp, 3.0_dp, 20.0_dp]
  real(dp), allocatable :: uu(:)
  integer  :: i, k, l, s, Ns(7) = [20, 40, 80, 160, 320, 640, 1280]
  integer  :: llist(6) = [0,0,1,1,2,3], slist(6) = [1,2,1,2,1,1]

  print '(a)', repeat('=',78)
  print '(a)', 'PART A  Spherical Bessel functions (Fortran)'
  print '(a)', repeat('=',78)
  call j_miller(30, 1.0_dp, jv); call j_up(30, 1.0_dp, ju)
  print '(a)', ' j_l(1): upward vs Miller (relative difference)'
  do l = 5, 30, 5
     print '(a,i3,a,es24.15,a,es10.2)', '   l=', l, '  Miller=', jv(l), '  rel.diff upward=', abs(ju(l)/jv(l) - 1)
  end do
  do k = 1, 3
     x = xc(k)
     call j_miller(25, x, jv(0:25)); call y_up(25, x, yv(0:25))
     cp = 0.0_dp
     do l = 1, 25
        cp = max(cp, abs(x*x*(jv(l)*yv(l-1) - jv(l-1)*yv(l)) - 1))
     end do
     print '(a,f5.1,a,es9.1)', ' cross product x=', x, ':  max_l |x^2 (j_l y_{l-1} - j_{l-1} y_l) - 1| = ', cp
  end do
  print '(a,es24.15,a)', ' j_3(1e-4) (series) = ', jl(3, 1.0d-4), '   (exact 9.523809518518520e-15)'
  do l = 0, 3
     call jl_zeros(l, 4, z(l,:))
     print '(a,i1,a,4f16.12)', ' zeros j_', l, ':', z(l,:)
  end do

  print '(/a)', repeat('=',78)
  print '(a)', 'PART B  Spherical infinite well, Numerov, N = 400'
  print '(a)', repeat('=',78)
  print '(a2,a3,a17,a17,a10,a16,a12)', 'l', 's', 'E Numerov', 'z^2/2', '|dE|', 'norm', 'max|u-ue|'
  do k = 1, 6
     l = llist(k); s = slist(k); Eex = z(l,s)**2/2
     call eigen(l, Eex, NN, .false., E, u)
     uex(0) = 0.0_dp
     do i = 1, NN
        r = i/real(NN,dp)
        uex(i) = r*jl(l, z(l,s)*r)*sqrt(2.0_dp)/abs(jl(l+1, z(l,s)))
     end do
     if (sum(u*uex) < 0) u = -u
     print '(i2,i3,2f17.12,es10.1,f16.12,es12.1)', l, s, E, Eex, abs(E - Eex), simpson(u*u, 1.0_dp/NN, NN), &
          maxval(abs(u - uex))
  end do

  open(10, file='ch4_f_conv.dat', status='replace')
  print '(/a)', 'Error of E_{l,1} vs N: naive start | Frobenius values at first nodes'
  do l = 0, 2
     Eex = z(l,1)**2/2
     do i = 1, 7
        allocate(uu(0:Ns(i)))
        call eigen(l, Eex, Ns(i), .true.,  E, uu); e1(i) = abs(E - Eex)
        call eigen(l, Eex, Ns(i), .false., E, uu); e2(i) = abs(E - Eex)
        write(10,'(i3,i6,3es16.7)') l, Ns(i), 1.0_dp/Ns(i), e1(i), e2(i)
        deallocate(uu)
     end do
     print '(a,i1,a,7es9.1)', ' l=', l, ' naive    :', e1
     print '(a,i1,a,7es9.1)', ' l=', l, ' Frobenius:', e2
  end do
  close(10)
end program ch4
