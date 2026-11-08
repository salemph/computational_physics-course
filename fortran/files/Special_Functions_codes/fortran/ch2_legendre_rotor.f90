!=======================================================================
! Chapter 2 -- LEGENDRE FUNCTIONS and the ANGULAR EQUATION  (Fortran 2008)
! Part A: P_l by Bonnet recurrence vs explicit sum, orthonormal associated
!         functions Pbar_l^m by the stable l-recurrence, Gauss-Legendre
!         rule by Newton, orthonormality and quadrature-error checks
! Part B: u = sqrt(sin t) Theta,  u'' + [lam + 1/4 + (1/4-m^2)/sin^2 t] u = 0
!         Numerov from both poles, Frobenius start at t_s (or first nodes),
!         Wronskian matching at t = pi/2, lam = l(l+1)
! Output: screen tables + ch2_f_conv.dat
! Compile: gfortran -O2 -o ch2 ch2_legendre_rotor.f90 && ./ch2
!=======================================================================
module legendre_mod
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
contains
  subroutine leg_rec(L, x, P)
    integer,  intent(in)  :: L
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: P(0:L)
    integer :: k
    P(0) = 1.0_dp
    if (L >= 1) P(1) = x
    do k = 1, L-1
       P(k+1) = ((2*k + 1)*x*P(k) - k*P(k-1))/(k + 1)
    end do
  end subroutine leg_rec

  real(dp) function leg_sum(l, x)                 ! explicit sum (cancellation)
    integer,  intent(in) :: l
    real(dp), intent(in) :: x
    integer :: k
    leg_sum = 0.0_dp
    do k = 0, l/2
       leg_sum = leg_sum + (-1)**k*exp(log_gamma(l+1.0_dp) - log_gamma(k+1.0_dp) - log_gamma(l-k+1.0_dp) &
                 + log_gamma(2*l-2*k+1.0_dp) - log_gamma(l+1.0_dp) - log_gamma(l-2*k+1.0_dp))*x**(l-2*k)
    end do
    leg_sum = leg_sum/2.0_dp**l
  end function leg_sum

  ! orthonormal associated Legendre functions Pbar_l^m, l = m..L
  subroutine plm_bar(L, m, x, P)
    integer,  intent(in)  :: L, m
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: P(0:L)
    real(dp) :: s, p0, al, alm1
    integer  :: k
    P = 0.0_dp
    s = sqrt(max(1.0_dp - x*x, 0.0_dp))
    p0 = 1.0_dp/sqrt(2.0_dp)
    do k = 1, m
       p0 = -sqrt((2*k + 1.0_dp)/(2*k))*s*p0
    end do
    P(m) = p0
    if (L > m) P(m+1) = sqrt(2*m + 3.0_dp)*x*p0
    do k = m+2, L
       al   = sqrt((4.0_dp*k*k - 1)/(real(k,dp)**2 - m*m))
       alm1 = sqrt((4.0_dp*(k-1)**2 - 1)/(real(k-1,dp)**2 - m*m))
       P(k) = al*(x*P(k-1) - P(k-2)/alm1)
    end do
  end subroutine plm_bar

  subroutine gauss_legendre(n, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: x(n), w(n)
    real(dp) :: z, dz, P(0:n), dPn
    integer  :: i, it
    do i = 1, n
       z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
       do it = 1, 100
          call leg_rec(n, z, P)
          dPn = n*(z*P(n) - P(n-1))/(z*z - 1)
          dz = P(n)/dPn; z = z - dz
          if (abs(dz) < 1.0d-15) exit
       end do
       call leg_rec(n, z, P)
       dPn = n*(z*P(n) - P(n-1))/(z*z - 1)
       x(n+1-i) = z; w(n+1-i) = 2/((1 - z*z)*dPn*dPn)
    end do
  end subroutine gauss_legendre
end module legendre_mod

module rotor_mod
  use legendre_mod
  implicit none
  real(dp), parameter :: b(0:7) = [1.0_dp, 1/3.0_dp, 1/15.0_dp, 2/189.0_dp, 1/675.0_dp, &
                                   2/10395.0_dp, 1382/58046625.0_dp, 4/1403325.0_dp]
contains
  real(dp) function frob(t, m, lam)                ! u = t^{m+1/2} sum a_k t^{2k}
    real(dp), intent(in) :: t, lam
    integer,  intent(in) :: m
    real(dp) :: a(0:23), s, q, LL
    integer  :: k, j
    q = 0.25_dp - m*m; LL = lam + 0.25_dp
    a(0) = 1.0_dp
    do k = 1, 23
       s = -LL*a(k-1)
       do j = 1, min(k, 7)
          s = s - q*b(j)*a(k-j)
       end do
       a(k) = s/(2*k*(2*k + 2*m))
    end do
    s = 0.0_dp
    do k = 0, 23
       s = s + a(k)*t**(2*k)
    end do
    frob = t**(m + 0.5_dp)*s
  end function frob

  real(dp) function simpson(y, h, n)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    simpson = h/3*(y(0) + y(n) + 4*sum(y(1:n-1:2)) + 2*sum(y(2:n-2:2)))
  end function simpson

  subroutine shoot(lam, m, N, ts, D, u)
    real(dp), intent(in)  :: lam, ts
    integer,  intent(in)  :: m, N
    real(dp), intent(out) :: D, u(0:N)
    real(dp) :: h, t(0:N), c(0:N), pL(0:N), pR(0:N), sL, sR
    integer  :: j, js, Mm, jg
    h = pi/N; t = [(j*h, j=0,N)]; t(N) = pi
    c = 1.0_dp
    c(1:N-1) = 1 + h*h*(lam + 0.25_dp + (0.25_dp - m*m)/sin(t(1:N-1))**2)/12
    js = max(1, nint(ts/h)); Mm = N/2
    pL = 0.0_dp; pR = 0.0_dp
    do j = 0, js+1
       pL(j) = frob(t(j), m, lam); pR(N-j) = frob(t(j), m, lam)
    end do
    do j = js+1, Mm
       pL(j+1) = ((12 - 10*c(j))*pL(j) - c(j-1)*pL(j-1))/c(j+1)
    end do
    do j = N-js-1, Mm, -1
       pR(j-1) = ((12 - 10*c(j))*pR(j) - c(j+1)*pR(j+1))/c(j-1)
    end do
    sL = maxval(abs(pL(0:Mm+1))); sR = maxval(abs(pR(Mm-1:N)))
    pL = pL/sL; pR = pR/sR
    D = pL(Mm+1)*pR(Mm) - pL(Mm)*pR(Mm+1)
    jg = Mm; if (abs(pR(Mm+1)) > abs(pR(Mm))) jg = Mm + 1
    u(0:jg) = pL(0:jg); u(jg+1:N) = pR(jg+1:N)*pL(jg)/pR(jg)
  end subroutine shoot

  subroutine eigen(l, m, N, ts, lam, u)
    integer,  intent(in)  :: l, m, N
    real(dp), intent(in)  :: ts
    real(dp), intent(out) :: lam, u(0:N)
    real(dp) :: a, bb, fa, fm, mid, D
    integer  :: it
    a = l*(l+1) - l - 0.5_dp; bb = l*(l+1) + l + 1.5_dp
    call shoot(a, m, N, ts, fa, u)
    do it = 1, 64
       mid = 0.5_dp*(a + bb)
       call shoot(mid, m, N, ts, fm, u)
       if (sign(1.0_dp, fm) == sign(1.0_dp, fa)) then
          a = mid; fa = fm
       else
          bb = mid
       end if
    end do
    lam = 0.5_dp*(a + bb)
    call shoot(lam, m, N, ts, D, u)
    u = u/sqrt(simpson(u*u, pi/N, N))
  end subroutine eigen
end module rotor_mod

program ch2
  use legendre_mod
  use rotor_mod
  implicit none
  integer, parameter :: NN = 400
  real(dp) :: P(0:200), xq(30), wq(30), G(18,18), Pq(0:20), u(0:NN), ex(0:NN), lam, c2, c2ex, mx, t
  real(dp), allocatable :: uu(:)
  integer  :: i, j, k, l, m, Ns(7) = [20, 40, 80, 160, 320, 640, 1280]
  integer  :: ll(9) = [0,1,2,1,2,3,2,3,5], mm(9) = [0,0,0,1,1,1,2,2,3], lconv(3) = [2, 2, 3]
  real(dp) :: e1(7), e2(7), xA(3) = [0.3_dp, 0.3_dp, 0.9_dp]
  integer  :: lA(3) = [10, 40, 100]

  print '(a)', repeat('=',78)
  print '(a)', 'PART A  Legendre functions (Fortran)'
  print '(a)', repeat('=',78)
  do k = 1, 3
     l = lA(k); t = xA(k)
     call leg_rec(l, t, P(0:l))
     print '(a,i4,a,f4.1,a,es24.15,a,es24.15,a,es9.2)', ' l=', l, ' x=', t, '  Bonnet=', P(l), '  sum=', &
          leg_sum(l, t), '  rel.diff=', abs(leg_sum(l, t)/P(l) - 1)
  end do
  call plm_bar(200, 180, 0.5_dp, P)
  print '(a,es24.15,a)', ' Pbar_200^180(0.5) = ', P(200), '   (mpmath: 7.231608361615227e-02)'
  call gauss_legendre(30, xq, wq)
  G = 0.0_dp
  do k = 1, 30
     call plm_bar(20, 3, xq(k), Pq)
     do i = 3, 20
        do j = 3, 20
           G(i-2,j-2) = G(i-2,j-2) + wq(k)*Pq(i)*Pq(j)
        end do
     end do
  end do
  do i = 1, 18
     G(i,i) = G(i,i) - 1.0_dp
  end do
  print '(a,es9.2)', ' orthonormality Pbar_l^3 (l=3..20), 30-pt Gauss-Legendre: max|G-I| = ', maxval(abs(G))
  print '(a)', ' Gauss-Legendre errors, int_{-1}^{1} e^x dx and int 1/(1+25x^2) dx:'
  do k = 2, 30, 4
     call gauss_legendre(k, xq(1:k), wq(1:k))
     print '(a,i3,2es11.2)', '   n=', k, abs(sum(wq(1:k)*exp(xq(1:k))) - (exp(1.0_dp) - exp(-1.0_dp))), &
          abs(sum(wq(1:k)/(1 + 25*xq(1:k)**2)) - 0.4_dp*atan(5.0_dp))
  end do

  print '(/a)', repeat('=',78)
  print '(a)', 'PART B  Angular equation, Numerov, N = 400, t_s = 0.4'
  print '(a)', repeat('=',78)
  print '(a2,a3,a18,a7,a10,a16,a12,a10,a10)', 'l', 'm', 'lambda', 'l(l+1)', '|dlam|', 'norm', 'max|u-ue|', &
       '<cos^2>', 'exact'
  do k = 1, 9
     l = ll(k); m = mm(k)
     call eigen(l, m, NN, 0.4_dp, lam, u)
     do j = 0, NN
        t = j*pi/NN
        call plm_bar(l, m, cos(t), P(0:l))
        ex(j) = sqrt(abs(sin(t)))*P(l)
     end do
     if (sum(u*ex) < 0) u = -u
     mx = maxval(abs(u - ex))
     c2 = simpson(u*u*cos([(j*pi/NN, j=0,NN)])**2, pi/NN, NN)
     c2ex = (2.0_dp*l*(l+1) - 2*m*m - 1)/((2*l - 1)*(2*l + 3))
     print '(i2,i3,f18.12,i7,es10.1,f16.12,es12.1,2f10.7)', l, m, lam, l*(l+1), abs(lam - l*(l+1)), &
          simpson(u*u, pi/NN, NN), mx, c2, c2ex
  end do

  open(10, file='ch2_f_conv.dat', status='replace')
  print '(/a)', 'Error of lambda vs N: Frobenius start at first nodes | at t_s = 0.4'
  do m = 0, 2
     do i = 1, 7
        allocate(uu(0:Ns(i)))
        call eigen(lconv(m+1), m, Ns(i), 0.0_dp, lam, uu); e1(i) = abs(lam - lconv(m+1)*(lconv(m+1)+1))
        call eigen(lconv(m+1), m, Ns(i), 0.4_dp, lam, uu); e2(i) = abs(lam - lconv(m+1)*(lconv(m+1)+1))
        write(10,'(i3,i6,3es16.7)') m, Ns(i), pi/Ns(i), e1(i), e2(i)
        deallocate(uu)
     end do
     print '(a,i1,a,7es9.1)', ' m=', m, ' first nodes:', e1
     print '(a,i1,a,7es9.1)', ' m=', m, ' t_s = 0.4  :', e2
  end do
  close(10)
end program ch2
