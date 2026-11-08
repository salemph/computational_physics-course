!=======================================================================
! Problem 3 -- QUANTUM HARMONIC OSCILLATOR with the Numerov method
! Units hbar = m = omega = 1.  V = x^2/2,  exact E_n = n + 1/2
! psi'' = -k^2(x) psi ,  k^2 = 2E - x^2
!
! Two-sided shooting: outward from -L, inward from +L, matched at the
! node m nearest to the right classical turning point sqrt(2E).
! Matching function = discrete Wronskian (no poles):
!     D(E) = psiL(m+1) psiR(m) - psiL(m) psiR(m+1)
! Bracketing by an energy scan, refinement by bisection, labelling by
! node counting.  Verification: Simpson normalisation, comparison with
! Hermite functions, orthogonality, virial theorem <V> = E/2.
! Error analysis: E_n error vs h (order 4), vs box half-width L.
! Output : screen tables + p3_states.dat, p3_conv_h.dat, p3_conv_L.dat
! Compile: gfortran -O2 -o p3 p3_harmonic_oscillator.f90 && ./p3
!=======================================================================
module oscillator
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
  integer,  parameter :: maxst = 10
contains
  pure elemental function V(x)
    real(dp), intent(in) :: x
    real(dp) :: V
    V = 0.5_dp*x*x
  end function V

  function simpson(y, h, n) result(s)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    real(dp) :: s
    s = h/3.0_dp*( y(0) + y(n) + 4.0_dp*sum(y(1:n-1:2)) + 2.0_dp*sum(y(2:n-2:2)) )
  end function simpson

  ! exact normalised eigenfunction via the Hermite recurrence
  function psi_exact(n, x) result(p)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: p, H0, H1, H2, fact
    integer  :: j
    H0 = 1.0_dp; H1 = 2.0_dp*x; H2 = H0
    if (n == 0) then
       H2 = H0
    else if (n == 1) then
       H2 = H1
    else
       do j = 2, n
          H2 = 2.0_dp*x*H1 - 2.0_dp*(j-1)*H0
          H0 = H1; H1 = H2
       end do
    end if
    fact = 1.0_dp
    do j = 2, n
       fact = fact*j
    end do
    p = pi**(-0.25_dp)/sqrt(2.0_dp**n*fact)*H2*exp(-x*x/2.0_dp)
  end function psi_exact

  !---------------------------------------------------------------
  ! shooting for trial energy E on grid x(0:N); returns D and the
  ! assembled (unnormalised) wavefunction
  !---------------------------------------------------------------
  subroutine shoot(E, x, h, N, D, psi)
    integer,  intent(in)  :: N
    real(dp), intent(in)  :: E, x(0:N), h
    real(dp), intent(out) :: D, psi(0:N)
    real(dp) :: c(0:N), pL(0:N), pR(0:N), sL, sR
    integer  :: m, j
    c = 1.0_dp + h*h*2.0_dp*(E - V(x))/12.0_dp
    m = 0
    do j = 0, N                                 ! first node >= turning point
       if (x(j) >= sqrt(max(2.0_dp*E, 0.0_dp))) then
          m = j; exit
       end if
    end do
    m = max(2, min(m, N-2))
    pL = 0.0_dp; pR = 0.0_dp
    pL(0) = 0.0_dp; pL(1) = 1.0d-10             ! outward
    do j = 1, m
       pL(j+1) = ((12.0_dp - 10.0_dp*c(j))*pL(j) - c(j-1)*pL(j-1))/c(j+1)
    end do
    pR(N) = 0.0_dp; pR(N-1) = 1.0d-10           ! inward
    do j = N-1, m, -1
       pR(j-1) = ((12.0_dp - 10.0_dp*c(j))*pR(j) - c(j+1)*pR(j+1))/c(j-1)
    end do
    sL = maxval(abs(pL(0:m+1))); sR = maxval(abs(pR(m-1:N)))
    pL = pL/sL; pR = pR/sR
    D = pL(m+1)*pR(m) - pL(m)*pR(m+1)           ! discrete Wronskian
    psi(0:m)   = pL(0:m)
    psi(m+1:N) = pR(m+1:N)*(pL(m)/pR(m))
  end subroutine shoot

  integer function nodes(psi, N)
    integer,  intent(in) :: N
    real(dp), intent(in) :: psi(0:N)
    real(dp) :: thr, last
    integer  :: j
    thr = 1.0d-6*maxval(abs(psi)); nodes = 0; last = 0.0_dp
    do j = 0, N
       if (abs(psi(j)) > thr) then
          if (last*psi(j) < 0.0_dp) nodes = nodes + 1
          last = psi(j)
       end if
    end do
  end function nodes

  !---------------------------------------------------------------
  ! all eigenstates with E < Emax on [-L,L] with N intervals
  !---------------------------------------------------------------
  subroutine eigenstates(L, N, Emax, nfound, nlab, Eval, states)
    real(dp), intent(in)  :: L, Emax
    integer,  intent(in)  :: N
    integer,  intent(out) :: nfound, nlab(maxst)
    real(dp), intent(out) :: Eval(maxst), states(0:N, maxst)
    real(dp) :: x(0:N), psi(0:N), h, E, Eold, D, Dold, a, b, fa, mid, fm, s, dE
    integer  :: j, nn
    h = 2.0_dp*L/N; dE = 0.05_dp
    x = [(-L + j*h, j=0,N)]
    nfound = 0
    Eold = 0.05_dp; call shoot(Eold, x, h, N, Dold, psi)
    E = Eold
    do while (E + dE < Emax)
       E = E + dE
       call shoot(E, x, h, N, D, psi)
       if (sign(1.0_dp,D) /= sign(1.0_dp,Dold)) then
          a = E - dE; b = E; fa = Dold
          do while (b - a > 1.0d-14*max(1.0_dp, a))   ! bisection
             mid = 0.5_dp*(a + b)
             call shoot(mid, x, h, N, fm, psi)
             if (sign(1.0_dp,fm) == sign(1.0_dp,fa)) then
                a = mid; fa = fm
             else
                b = mid
             end if
          end do
          call shoot(0.5_dp*(a+b), x, h, N, fm, psi)
          psi = psi/sqrt(simpson(psi*psi, h, N))      ! normalisation
          nn = nodes(psi, N)
          s = 0.0_dp
          do j = 0, N
             s = s + psi(j)*psi_exact(nn, x(j))
          end do
          if (s < 0.0_dp) psi = -psi
          nfound = nfound + 1
          nlab(nfound) = nn; Eval(nfound) = 0.5_dp*(a+b); states(:,nfound) = psi
       end if
       Dold = D
    end do
  end subroutine eigenstates

  function order(h, e, m) result(p)
    integer,  intent(in) :: m
    real(dp), intent(in) :: h(m), e(m)
    real(dp) :: p, lx(m), ly(m)
    lx = log(h); ly = log(e)
    p = (m*sum(lx*ly) - sum(lx)*sum(ly)) / (m*sum(lx*lx) - sum(lx)**2)
  end function order

  real(dp) function Efor(nwant, nfound, nlab, Eval)
    integer,  intent(in) :: nwant, nfound, nlab(maxst)
    real(dp), intent(in) :: Eval(maxst)
    integer :: j
    Efor = huge(1.0_dp)
    do j = 1, nfound
       if (nlab(j) == nwant) Efor = Eval(j)
    end do
  end function Efor
end module oscillator

program harmonic_oscillator
  use oscillator
  implicit none
  integer, parameter :: NL = 1400, nh = 8, nLs = 8
  integer  :: Ns(nh) = [70, 140, 280, 560, 1120, 2240, 4480, 8960]
  real(dp) :: Ls(nLs) = [3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp, 5.0_dp, 5.5_dp, 6.0_dp, 7.0_dp]
  integer  :: wanted(3) = [0, 2, 4]
  real(dp), allocatable :: st(:,:), x(:)
  real(dp) :: Eval(maxst), h, Vexp, ex, maxerr, hs(nh), errs(3,nh), errL(3,nLs), omax, ov
  real(dp) :: L, E1, E2, ER
  integer  :: nlab(maxst), nf, i, j, q, N

  L = 7.0_dp; h = 2.0_dp*L/NL
  allocate(st(0:NL, maxst), x(0:NL))
  x = [(-L + j*h, j=0,NL)]
  call eigenstates(L, NL, 6.0_dp, nf, nlab, Eval, st)
  print '(a)', repeat('=',78)
  print '(a)', 'Harmonic oscillator  V = x^2/2, box [-7, 7], N = 1400 (h = 0.01)'
  print '(a)', repeat('=',78)
  print '(a2,a18,a8,a11,a16,a13,a11,a9)', 'n', 'E_numerov', 'E_exact', '|dE|', 'norm', &
       'max|psi-ex|', '<V>', 'E/2'
  open(10, file='p3_states.dat', status='replace')
  do i = 1, nf
     maxerr = 0.0_dp
     do j = 0, NL
        ex = psi_exact(nlab(i), x(j)); maxerr = max(maxerr, abs(st(j,i) - ex))
     end do
     Vexp = simpson(st(:,i)**2*V(x), h, NL)
     print '(i2,f18.12,f8.1,es11.2,f16.12,es13.2,f11.6,f9.4)', nlab(i), Eval(i), nlab(i)+0.5_dp, &
          abs(Eval(i) - nlab(i) - 0.5_dp), simpson(st(:,i)**2, h, NL), maxerr, Vexp, Eval(i)/2
  end do
  do j = 0, NL
     write(10,'(7es16.7)') x(j), st(j,1:nf)
  end do
  close(10)
  omax = 0.0_dp
  do i = 1, nf
     do q = 1, nf
        ov = simpson(st(:,i)*st(:,q), h, NL)
        if (i == q) ov = ov - 1.0_dp
        omax = max(omax, abs(ov))
     end do
  end do
  print '(a,es9.2)', 'max |<psi_i|psi_j> - delta_ij| = ', omax
  deallocate(st, x)

  ! ---------------- convergence in h
  open(11, file='p3_conv_h.dat', status='replace')
  print '(/a)', 'Eigenvalue error vs step (L = 7):'
  print '(a6,a9,3a12)', 'N', 'h', '|dE_0|', '|dE_2|', '|dE_4|'
  do i = 1, nh
     N = Ns(i); hs(i) = 2.0_dp*L/N
     allocate(st(0:N, maxst))
     call eigenstates(L, N, 6.0_dp, nf, nlab, Eval, st)
     do q = 1, 3
        errs(q,i) = abs(Efor(wanted(q), nf, nlab, Eval) - wanted(q) - 0.5_dp)
     end do
     print '(i6,f9.5,3es12.2)', N, hs(i), errs(:,i)
     write(11,'(4es16.7)') hs(i), errs(:,i)
     deallocate(st)
  end do
  close(11)
  print '(a,3f6.2)', 'observed order p (n=0,2,4): ', (order(hs(1:5), errs(q,1:5), 5), q=1,3)

  allocate(st(0:280, maxst))
  call eigenstates(L, 140, 6.0_dp, nf, nlab, Eval, st(0:140,:))
  E1 = Efor(4, nf, nlab, Eval)
  call eigenstates(L, 280, 6.0_dp, nf, nlab, Eval, st)
  E2 = Efor(4, nf, nlab, Eval)
  ER = (16.0_dp*E2 - E1)/15.0_dp
  print '(/a,f13.10,a,f13.10,a,f13.10,a,es8.1)', 'Richardson n=4: E(0.1)=', E1, '  E(0.05)=', E2, &
       '  E_R=', ER, '  err=', abs(ER - 4.5_dp)
  deallocate(st)

  ! ---------------- truncation error vs L
  open(12, file='p3_conv_L.dat', status='replace')
  print '(/a)', 'Box-truncation error (h = 0.005):'
  print '(a5,3a12)', 'L', '|dE_0|', '|dE_2|', '|dE_4|'
  do i = 1, nLs
     N = nint(2.0_dp*Ls(i)/0.005_dp); N = N + mod(N,2)
     allocate(st(0:N, maxst))
     call eigenstates(Ls(i), N, 6.0_dp, nf, nlab, Eval, st)
     do q = 1, 3
        errL(q,i) = abs(Efor(wanted(q), nf, nlab, Eval) - wanted(q) - 0.5_dp)
     end do
     print '(f5.1,3es12.2)', Ls(i), errL(:,i)
     write(12,'(4es16.7)') Ls(i), errL(:,i)
     deallocate(st)
  end do
  close(12)
end program harmonic_oscillator
