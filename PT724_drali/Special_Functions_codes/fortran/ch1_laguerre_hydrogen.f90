!=======================================================================
! Chapter 1 -- LAGUERRE POLYNOMIALS and the HYDROGEN ATOM  (Fortran 2008)
! Part A: L_n^(a)(x) by recurrence and by explicit sum, Gauss-Laguerre
!         nodes/weights by Newton, orthogonality check, quadrature errors
! Part B: radial Schroedinger equation for hydrogen (atomic units, Z=1)
!         u'' = [l(l+1)/r^2 - 2/r - 2E] u ,   E_n = -1/(2n^2)
!         Numerov, two-sided shooting, Frobenius start, Wronskian matching
! Output: screen tables + ch1_f_states.dat, ch1_f_conv.dat
! Compile: gfortran -O2 -o ch1 ch1_laguerre_hydrogen.f90 && ./ch1
!=======================================================================
module laguerre_mod
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
contains
  ! L_0..L_n^(a)(x) by (k+1)L_{k+1} = (2k+1+a-x)L_k - (k+a)L_{k-1}
  subroutine lag_rec(n, a, x, L)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a, x
    real(dp), intent(out) :: L(0:n)
    integer :: k
    L(0) = 1.0_dp
    if (n >= 1) L(1) = 1.0_dp + a - x
    do k = 1, n-1
       L(k+1) = ((2*k + 1 + a - x)*L(k) - (k + a)*L(k-1))/(k + 1)
    end do
  end subroutine lag_rec

  real(dp) function lag(n, a, x)
    integer,  intent(in) :: n
    real(dp), intent(in) :: a, x
    real(dp) :: L(0:max(n,1))
    call lag_rec(max(n,1), a, x, L); lag = L(n)
  end function lag

  ! explicit sum  sum_k (-1)^k C(n+a,n-k) x^k/k!   (suffers cancellation)
  real(dp) function lag_sum(n, a, x)
    integer,  intent(in) :: n
    real(dp), intent(in) :: a, x
    integer :: k
    lag_sum = 0.0_dp
    do k = 0, n
       lag_sum = lag_sum + (-1)**k*exp(log_gamma(n+a+1) - log_gamma(real(n-k+1,dp)) &
                 - log_gamma(a+k+1) - log_gamma(real(k+1,dp)))*x**k
    end do
  end function lag_sum

  ! Gauss-Laguerre: Newton on L_n^(a) with Stroud-Secrest initial guesses
  subroutine gauss_laguerre(n, a, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a
    real(dp), intent(out) :: x(n), w(n)
    real(dp) :: z, dz, L(0:n+1), dL
    integer  :: i, it, ai
    do i = 1, n
       if (i == 1) then
          z = (1 + a)*(3 + 0.92_dp*a)/(1 + 2.4_dp*n + 1.8_dp*a)
       else if (i == 2) then
          z = x(1) + (15 + 6.25_dp*a)/(1 + 0.9_dp*a + 2.5_dp*n)
       else
          ai = i - 2
          z = x(i-1) + ((1 + 2.55_dp*ai)/(1.9_dp*ai) + 1.26_dp*ai*a/(1 + 3.5_dp*ai)) &
              *(x(i-1) - x(i-2))/(1 + 0.3_dp*a)
       end if
       do it = 1, 100
          call lag_rec(n, a, z, L(0:n))
          dL = (n*L(n) - (n + a)*L(n-1))/z
          dz = L(n)/dL; z = z - dz
          if (abs(dz) < 1.0d-15*max(1.0_dp, z)) exit
       end do
       x(i) = z
       call lag_rec(n+1, a, z, L)
       w(i) = exp(log_gamma(n+a+1) - log_gamma(real(n+1,dp)))*z/((n+1)**2*L(n+1)**2)
    end do
  end subroutine gauss_laguerre
end module laguerre_mod

module hydrogen_mod
  use laguerre_mod
  implicit none
contains
  ! normalised exact radial function u = r R_nl, using the Laguerre recurrence
  real(dp) function u_exact(n, l, r)
    integer,  intent(in) :: n, l
    real(dp), intent(in) :: r
    real(dp) :: rho, nrm
    rho = 2*r/n
    nrm = sqrt((2.0_dp/n)**3*exp(log_gamma(real(n-l,dp)) - log_gamma(real(n+l+1,dp)))/(2*n))
    u_exact = r*nrm*exp(-rho/2)*rho**l*lag(n-l-1, real(2*l+1,dp), rho)
  end function u_exact

  ! Frobenius series  u = r^{l+1} sum a_k r^k
  real(dp) function frob(r, l, E)
    real(dp), intent(in) :: r, E
    integer,  intent(in) :: l
    real(dp) :: a(0:39), s, p
    integer  :: k
    a(0) = 1.0_dp; a(1) = -2.0_dp/(2*l + 2)
    do k = 2, 39
       a(k) = (-2*a(k-1) - 2*E*a(k-2))/(k*(k + 2*l + 1))
    end do
    s = 0.0_dp; p = 1.0_dp
    do k = 0, 39
       s = s + a(k)*p; p = p*r
    end do
    frob = r**(l+1)*s
  end function frob

  real(dp) function simpson(y, h, n)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    simpson = h/3*(y(0) + y(n) + 4*sum(y(1:n-1:2)) + 2*sum(y(2:n-2:2)))
  end function simpson

  ! two-sided shooting; naive = .true. -> textbook start u0=0, u1=h^(l+1)
  subroutine shoot(E, l, h, N, naive, D, u)
    real(dp), intent(in)  :: E, h
    integer,  intent(in)  :: l, N
    logical,  intent(in)  :: naive
    real(dp), intent(out) :: D, u(0:N)
    real(dp) :: r(0:N), c(0:N), pL(0:N), pR(0:N), rt, sL, sR
    integer  :: j, m, js
    r = [(j*h, j=0,N)]
    c(0) = 1.0_dp
    c(1:N) = 1 + h*h*(2*E + 2/r(1:N) - l*(l+1)/r(1:N)**2)/12
    rt = (1 + sqrt(max(1 + 2*E*l*(l+1), 0.0_dp)))/(-2*E)
    m = min(max(nint(rt/h), nint(0.25_dp/h)+4), N-3)
    pL = 0.0_dp
    if (naive) then
       js = 0; pL(1) = h**(l+1)
    else
       js = max(1, nint(0.25_dp/h))                    ! Frobenius start at r_s = 0.25
       do j = 0, js+1
          pL(j) = frob(r(j), l, E)
       end do
    end if
    do j = js+1, m
       pL(j+1) = ((12 - 10*c(j))*pL(j) - c(j-1)*pL(j-1))/c(j+1)
    end do
    pR = 0.0_dp; pR(N-1) = 1.0d-30
    do j = N-1, m, -1
       pR(j-1) = ((12 - 10*c(j))*pR(j) - c(j+1)*pR(j+1))/c(j-1)
    end do
    sL = maxval(abs(pL(0:m+1))); sR = maxval(abs(pR(m-1:N)))
    pL = pL/sL; pR = pR/sR
    D = pL(m+1)*pR(m) - pL(m)*pR(m+1)
    u(0:m) = pL(0:m); u(m+1:N) = pR(m+1:N)*pL(m)/pR(m)
  end subroutine shoot

  subroutine eigen(n, l, h, naive, E, u, Nout)
    integer,  intent(in)  :: n, l
    real(dp), intent(in)  :: h
    logical,  intent(in)  :: naive
    real(dp), intent(out) :: E
    real(dp), allocatable, intent(out) :: u(:)
    integer,  intent(out) :: Nout
    real(dp) :: a, b, fa, fm, mid, Ex, D
    integer  :: it
    Ex = -0.5_dp/n**2
    Nout = nint(40.0_dp*n/h); Nout = Nout + mod(Nout, 2)
    allocate(u(0:Nout))
    a = 1.25_dp*Ex; b = 0.75_dp*Ex
    call shoot(a, l, h, Nout, naive, fa, u)
    do it = 1, 60
       mid = 0.5_dp*(a + b)
       call shoot(mid, l, h, Nout, naive, fm, u)
       if (sign(1.0_dp, fm) == sign(1.0_dp, fa)) then
          a = mid; fa = fm
       else
          b = mid
       end if
    end do
    E = 0.5_dp*(a + b)
    call shoot(E, l, h, Nout, naive, D, u)
    u = u/sqrt(simpson(u*u, h, Nout))
  end subroutine eigen
end module hydrogen_mod

program ch1
  use laguerre_mod
  use hydrogen_mod
  implicit none
  real(dp) :: xq(40), wq(40), Lm(0:12), G(0:12,0:12), hn, maxo, E, x
  real(dp), allocatable :: u(:)
  real(dp) :: hs(5) = [0.08_dp, 0.04_dp, 0.02_dp, 0.01_dp, 0.005_dp], errn(5), errf(5)
  real(dp) :: rexp, vexp, mx, rr
  integer  :: i, j, k, n, l, ng, cases(2,3)
  integer  :: nlist(6) = [1, 2, 2, 3, 3, 3], llist(6) = [0, 0, 1, 0, 1, 2]
  real(dp) :: xs(4) = [0.5_dp, 5.0_dp, 20.0_dp, 60.0_dp]

  print '(a)', repeat('=',78)
  print '(a)', 'PART A  Laguerre polynomials (Fortran)'
  print '(a)', repeat('=',78)
  print '(a)', 'L_30^(1)(x): recurrence vs explicit sum (difference shows the cancellation of the sum)'
  do i = 1, 4
     x = xs(i)
     print '(a,f5.1,a,es24.15,a,es24.15,a,es9.2)', ' x=', x, '  rec=', lag(30, 1.0_dp, x), &
          '  sum=', lag_sum(30, 1.0_dp, x), '  rel.diff=', abs(lag_sum(30,1.0_dp,x)/lag(30,1.0_dp,x) - 1)
  end do
  call gauss_laguerre(40, 2.0_dp, xq, wq)
  G = 0.0_dp
  do k = 1, 40
     call lag_rec(12, 2.0_dp, xq(k), Lm)
     do i = 0, 12
        do j = 0, 12
           G(i,j) = G(i,j) + wq(k)*Lm(i)*Lm(j)
        end do
     end do
  end do
  maxo = 0.0_dp
  do i = 0, 12
     do j = 0, 12
        hn = sqrt(exp(log_gamma(i+3.0_dp) - log_gamma(i+1.0_dp))*exp(log_gamma(j+3.0_dp) - log_gamma(j+1.0_dp)))
        maxo = max(maxo, abs(G(i,j)/hn - merge(1.0_dp, 0.0_dp, i == j)))
     end do
  end do
  print '(a,es9.2)', ' orthogonality a=2, n,m<=12, 40-point Gauss-Laguerre: max|G/h - delta| = ', maxo
  print '(a)', ' Gauss-Laguerre errors for int e^-x cos x = 1/2 and int e^-x sqrt(x) = sqrt(pi)/2:'
  do n = 4, 40, 6
     call gauss_laguerre(n, 0.0_dp, xq(1:n), wq(1:n))
     print '(a,i3,a,es10.2,a,es10.2)', '   n=', n, '  cos: ', abs(sum(wq(1:n)*cos(xq(1:n))) - 0.5_dp), &
          '   sqrt: ', abs(sum(wq(1:n)*sqrt(xq(1:n))) - sqrt(pi)/2)
  end do

  print '(/a)', repeat('=',78)
  print '(a)', 'PART B  Hydrogen atom, Numerov h = 0.01, R = 40 n, Frobenius start'
  print '(a)', repeat('=',78)
  print '(a2,a3,a17,a14,a10,a16,a12,a10,a10,a11)', 'n', 'l', 'E Numerov', 'E exact', '|dE|', 'norm', &
       'max|u-ue|', '<r>', '<r>ex', '<V>/2E'
  open(10, file='ch1_f_states.dat', status='replace')
  do i = 1, 6
     n = nlist(i); l = llist(i)
     call eigen(n, l, 0.01_dp, .false., E, u, ng)
     mx = 0.0_dp
     if (sum(u*[(u_exact(n,l,k*0.01_dp), k=0,ng)]) < 0) u = -u     ! fix arbitrary sign
     do j = 0, ng
        rr = j*0.01_dp
        mx = max(mx, abs(u(j) - u_exact(n, l, rr)))
     end do
     rexp = simpson(u*u*[(j*0.01_dp, j=0,ng)], 0.01_dp, ng)
     vexp = -simpson([0.0_dp, u(1:ng)**2/[(j*0.01_dp, j=1,ng)]], 0.01_dp, ng)
     print '(i2,i3,f17.12,f14.10,es10.1,f16.12,es12.1,f10.5,f10.5,f11.7)', n, l, E, -0.5_dp/n**2, &
          abs(E + 0.5_dp/n**2), simpson(u*u, 0.01_dp, ng), mx, rexp, (3*n*n - l*(l+1))/2.0_dp, vexp/(2*E)
     write(10,'(2i3,2es24.15)') n, l, E, abs(E + 0.5_dp/n**2)
     deallocate(u)
  end do
  close(10)

  cases = reshape([1,0, 2,1, 3,2], [2,3])
  open(11, file='ch1_f_conv.dat', status='replace')
  print '(/a)', 'Eigenvalue error vs h: naive start u1 = h^(l+1)  |  Frobenius start'
  do k = 1, 3
     n = cases(1,k); l = cases(2,k)
     do i = 1, 5
        call eigen(n, l, hs(i), .true.,  E, u, ng); errn(i) = abs(E + 0.5_dp/n**2); deallocate(u)
        call eigen(n, l, hs(i), .false., E, u, ng); errf(i) = abs(E + 0.5_dp/n**2); deallocate(u)
        write(11,'(2i3,3es16.7)') n, l, hs(i), errn(i), errf(i)
     end do
     print '(a,i1,a,i1,a,5es10.2)', ' (', n, ',', l, ') naive    : ', errn
     print '(a,i1,a,i1,a,5es10.2)', ' (', n, ',', l, ') Frobenius: ', errf
  end do
  close(11)
end program ch1
