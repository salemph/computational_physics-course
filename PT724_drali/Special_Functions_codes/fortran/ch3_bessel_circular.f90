!=======================================================================
! Chapter 3 -- BESSEL FUNCTIONS J_n and the CIRCULAR INFINITE WELL (F2008)
! Part A: J_n(x) by Miller's backward recurrence (normalised with
!         J0 + 2 sum J_2k = 1), forward recurrence, power series,
!         zeros by McMahon + Newton; independent check with the Fortran
!         2008 intrinsic bessel_jn and the Jacobi-Anger identity
! Part B: u = sqrt(r) R,  u'' + [2E - (m^2-1/4)/r^2] u = 0,  u(1) = 0
!         Numerov outward shooting + bisection, Frobenius start
! Output: screen tables + ch3_f_conv.dat
! Compile: gfortran -O2 -o ch3 ch3_bessel_circular.f90 && ./ch3
!=======================================================================
module bessel_mod
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
contains
  ! J_0..J_nmax(x), x > 0, Miller backward recurrence from even index M
  subroutine j_miller(nmax, x, Jv, Mstart)
    integer,  intent(in)  :: nmax
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: Jv(0:nmax)
    integer,  intent(in), optional :: Mstart
    real(dp), allocatable :: f(:)
    real(dp) :: nrm
    integer  :: M, k
    if (present(Mstart)) then
       M = Mstart
    else
       M = 2*((max(nmax, int(x)) + int(sqrt(64.0_dp*max(nmax, int(x), 1))) + 10)/2)
    end if
    M = M + mod(M, 2)
    allocate(f(0:M+1))
    f = 0.0_dp; f(M) = 1.0d-30
    do k = M, 1, -1
       f(k-1) = 2*k/x*f(k) - f(k+1)
       if (abs(f(k-1)) > 1.0d250) f(k-1:M+1) = f(k-1:M+1)*1.0d-250
    end do
    nrm = f(0) + 2*sum(f(2:M:2))
    Jv = 0.0_dp
    Jv(0:min(nmax, M)) = f(0:min(nmax, M))/nrm
  end subroutine j_miller

  real(dp) function J(n, x)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: Jv(0:n+1)
    if (x == 0.0_dp) then
       J = merge(1.0_dp, 0.0_dp, n == 0); return
    end if
    call j_miller(n+1, x, Jv); J = Jv(n)
  end function J

  real(dp) function j_series(n, x)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: term
    integer  :: k
    term = (x/2)**n/gamma(n + 1.0_dp); j_series = 0.0_dp
    do k = 0, 300
       j_series = j_series + term
       term = -term*(x/2)**2/((k + 1)*(k + 1 + n))
       if (abs(term) < 1.0d-17*abs(j_series) .and. k > x) exit
    end do
  end function j_series

  subroutine bessel_zeros(n, S, z)
    integer,  intent(in)  :: n, S
    real(dp), intent(out) :: z(S)
    real(dp) :: bb, mu, x, dx, Jv(0:n+2), d
    integer  :: k, it
    mu = 4.0_dp*n*n
    do k = 1, S
       bb = (k + n/2.0_dp - 0.25_dp)*pi
       x = bb - (mu - 1)/(8*bb) - 4*(mu - 1)*(7*mu - 31)/(3*(8*bb)**3)
       do it = 1, 50
          call j_miller(n+2, x, Jv)
          if (n > 0) then
             d = (Jv(n-1) - Jv(n+1))/2
          else
             d = -Jv(1)
          end if
          dx = Jv(n)/d; x = x - dx
          if (abs(dx) < 1.0d-15*x) exit
       end do
       z(k) = x
    end do
  end subroutine bessel_zeros
end module bessel_mod

module well_mod
  use bessel_mod
  implicit none
contains
  real(dp) function frob(r, m, E)
    real(dp), intent(in) :: r, E
    integer,  intent(in) :: m
    real(dp) :: a, s, p
    integer  :: k
    s = 0.0_dp; a = 1.0_dp; p = 1.0_dp
    do k = 0, 39
       s = s + a*p
       a = -2*E*a/((2*k + 2)*(2*k + 2 + 2*m)); p = p*r*r
    end do
    frob = r**(m + 0.5_dp)*s
  end function frob

  real(dp) function simpson(y, h, n)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    simpson = h/3*(y(0) + y(n) + 4*sum(y(1:n-1:2)) + 2*sum(y(2:n-2:2)))
  end function simpson

  subroutine shoot(E, m, N, rs, u)
    real(dp), intent(in)  :: E, rs
    integer,  intent(in)  :: m, N
    real(dp), intent(out) :: u(0:N)
    real(dp) :: h, r(0:N), c(0:N)
    integer  :: j, js
    h = 1.0_dp/N; r = [(j*h, j=0,N)]
    c(0) = 1.0_dp; c(1:N) = 1 + h*h*(2*E - (m*m - 0.25_dp)/r(1:N)**2)/12
    js = max(1, nint(rs/h))
    u = 0.0_dp
    do j = 0, js+1
       u(j) = frob(r(j), m, E)
    end do
    do j = js+1, N-1
       u(j+1) = ((12 - 10*c(j))*u(j) - c(j-1)*u(j-1))/c(j+1)
    end do
  end subroutine shoot

  subroutine eigen(m, Ex, N, rs, E, u)
    integer,  intent(in)  :: m, N
    real(dp), intent(in)  :: Ex, rs
    real(dp), intent(out) :: E, u(0:N)
    real(dp) :: a, b, fa, mid
    integer  :: it
    a = 0.85_dp*Ex; b = 1.15_dp*Ex
    call shoot(a, m, N, rs, u); fa = u(N)
    do it = 1, 64
       mid = 0.5_dp*(a + b)
       call shoot(mid, m, N, rs, u)
       if (sign(1.0_dp, u(N)) == sign(1.0_dp, fa)) then
          a = mid; fa = u(N)
       else
          b = mid
       end if
    end do
    E = 0.5_dp*(a + b)
    call shoot(E, m, N, rs, u)
    u = u/sqrt(simpson(u*u, 1.0_dp/N, N))
  end subroutine eigen
end module well_mod

program ch3
  use bessel_mod
  use well_mod
  implicit none
  integer, parameter :: NN = 400
  real(dp) :: Jv(0:80), Jf(0:30), z(0:3,5), x, s, E, u(0:NN), uex(0:NN), r, Eex, e1(7), e2(7), maxd
  real(dp), allocatable :: uu(:)
  integer  :: i, k, m, n, sidx, Ns(7) = [20, 40, 80, 160, 320, 640, 1280]
  real(dp) :: xja(3) = [0.5_dp, 7.3_dp, 25.0_dp]
  integer  :: mlist(7) = [0,0,0,1,1,2,3], slist(7) = [1,2,3,1,2,1,1]

  print '(a)', repeat('=',78)
  print '(a)', 'PART A  Bessel functions J_n (Fortran); reference = intrinsic bessel_jn'
  print '(a)', repeat('=',78)
  maxd = 0.0_dp
  do i = 1, 400
     x = 0.1_dp*i
     call j_miller(10, x, Jv(0:10))
     do n = 0, 10
        maxd = max(maxd, abs(Jv(n) - bessel_jn(n, x)))
     end do
  end do
  print '(a,es9.2)', ' max |J_miller - bessel_jn| for n<=10, 0.1<=x<=40: ', maxd
  print '(a)', ' power series vs Miller (relative difference):'
  print '(a,2es24.15,es10.2)', '   J_0(30): ', j_series(0, 30.0_dp), J(0, 30.0_dp), abs(j_series(0,30.0_dp)/J(0,30.0_dp)-1)
  print '(a,2es24.15,es10.2)', '   J_3(50): ', j_series(3, 50.0_dp), J(3, 50.0_dp), abs(j_series(3,50.0_dp)/J(3,50.0_dp)-1)
  call j_miller(30, 1.0_dp, Jv(0:30))
  Jf(0) = Jv(0); Jf(1) = Jv(1)
  do n = 1, 29
     Jf(n+1) = 2*n/1.0_dp*Jf(n) - Jf(n-1)
  end do
  print '(a)', ' forward recurrence at x = 1, relative error vs Miller:'
  do n = 5, 30, 5
     print '(a,i3,a,es10.2,a,es10.2)', '   n=', n, '  forward ', abs(Jf(n)/Jv(n) - 1), &
          '   Miller vs bessel_jn ', abs(Jv(n)/bessel_jn(n, 1.0_dp) - 1)
  end do
  do k = 1, 3
     x = xja(k)
     call j_miller(80, x, Jv)
     s = Jv(0)
     do i = 1, 40
        s = s + 2*(-1)**i*Jv(2*i)
     end do
     print '(a,f5.1,a,es9.1)', ' Jacobi-Anger x=', x, ':  J0 - 2J2 + 2J4 - ... - cos x = ', s - cos(x)
  end do
  do m = 0, 3
     call bessel_zeros(m, 5, z(m,:))
     print '(a,i1,a,5f16.12)', ' zeros J_', m, ':', z(m,:)
  end do
  print '(a,4es9.1)', ' check |bessel_jn(m, j_m1)|:', (abs(bessel_jn(m, z(m,1))), m=0,3)

  print '(/a)', repeat('=',78)
  print '(a)', 'PART B  Circular infinite well, Numerov, N = 400, r_s = 0.125'
  print '(a)', repeat('=',78)
  print '(a2,a3,a17,a17,a10,a16,a12,a10,a10)', 'm', 's', 'E Numerov', 'j^2/2', '|dE|', 'norm', 'max|u-ue|', '<r^2>', 'exact'
  do k = 1, 7
     m = mlist(k); sidx = slist(k); Eex = z(m,sidx)**2/2
     call eigen(m, Eex, NN, 0.125_dp, E, u)
     do i = 0, NN
        r = i/real(NN,dp)
        uex(i) = sqrt(r)*bessel_jn(m, z(m,sidx)*r)*sqrt(2.0_dp)/abs(bessel_jn(m+1, z(m,sidx)))
     end do
     if (sum(u*uex) < 0) u = -u
     print '(i2,i3,2f17.12,es10.1,f16.12,es12.1,2f10.6)', m, sidx, E, Eex, abs(E - Eex), &
          simpson(u*u, 1.0_dp/NN, NN), maxval(abs(u - uex)), &
          simpson(u*u*[((i/real(NN,dp))**2, i=0,NN)], 1.0_dp/NN, NN), 1/3.0_dp + 2*(m*m - 1)/(3*z(m,sidx)**2)
  end do

  open(10, file='ch3_f_conv.dat', status='replace')
  print '(/a)', 'Error of E_{m,1} vs N: start at first nodes | at r_s = 1/8'
  do m = 0, 2
     Eex = z(m,1)**2/2
     do i = 1, 7
        allocate(uu(0:Ns(i)))
        call eigen(m, Eex, Ns(i), 0.0_dp, E, uu);   e1(i) = abs(E - Eex)
        call eigen(m, Eex, Ns(i), 0.125_dp, E, uu); e2(i) = abs(E - Eex)
        write(10,'(i3,i6,3es16.7)') m, Ns(i), 1.0_dp/Ns(i), e1(i), e2(i)
        deallocate(uu)
     end do
     print '(a,i1,a,7es9.1)', ' m=', m, ' first nodes:', e1
     print '(a,i1,a,7es9.1)', ' m=', m, ' r_s = 1/8  :', e2
  end do
  close(10)
end program ch3
