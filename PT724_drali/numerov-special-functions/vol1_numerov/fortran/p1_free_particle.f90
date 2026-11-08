!=======================================================================
! Problem 1 -- FREE PARTICLE solved with the Numerov method (Fortran 90)
! Units hbar = m = 1 :  psi'' = -k^2 psi ,  k^2 = 2(E - V),  V = 0
!
! Part A : initial-value problem psi(0)=0, k=2  -> exact psi = sin(kx)
!          global error vs h, exact start vs Taylor start psi(h)=kh
! Part B : free particle in a box [0,L], psi(0)=psi(L)=0
!          shooting + bisection, Simpson normalisation, comparison with
!          sqrt(2/L) sin(n pi x/L), convergence of E_n,  Richardson
! Output : screen tables + p1_partA.dat, p1_partB.dat, p1_states.dat
! Compile: gfortran -O2 -o p1 p1_free_particle.f90 && ./p1
!=======================================================================
module numerov_tools
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238_dp
contains
  !---------------------------------------------------------------
  ! Numerov recursion, forward direction, real wavefunction
  ! psi(n+1) c(n+1) = (12 - 10 c(n)) psi(n) - c(n-1) psi(n-1),
  ! c = 1 + h^2 k2 / 12
  !---------------------------------------------------------------
  subroutine numerov(k2, psi, h, n)
    integer,  intent(in)    :: n
    real(dp), intent(in)    :: k2(0:n), h
    real(dp), intent(inout) :: psi(0:n)       ! psi(0), psi(1) given on entry
    real(dp) :: c(0:n)
    integer  :: i
    c = 1.0_dp + h*h*k2/12.0_dp
    do i = 1, n-1
       psi(i+1) = ((12.0_dp - 10.0_dp*c(i))*psi(i) - c(i-1)*psi(i-1)) / c(i+1)
    end do
  end subroutine numerov

  !---------------------------------------------------------------
  ! composite Simpson rule on 0..n (n even)
  !---------------------------------------------------------------
  function simpson(y, h, n) result(s)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(0:n), h
    real(dp) :: s
    s = h/3.0_dp*( y(0) + y(n) + 4.0_dp*sum(y(1:n-1:2)) + 2.0_dp*sum(y(2:n-2:2)) )
  end function simpson

  !---------------------------------------------------------------
  ! least-squares slope of log(e) vs log(h)  -> observed order
  !---------------------------------------------------------------
  function order(h, e, m) result(p)
    integer,  intent(in) :: m
    real(dp), intent(in) :: h(m), e(m)
    real(dp) :: p, lx(m), ly(m)
    lx = log(h); ly = log(e)
    p = (m*sum(lx*ly) - sum(lx)*sum(ly)) / (m*sum(lx*lx) - sum(lx)**2)
  end function order
end module numerov_tools

module box_problem
  use numerov_tools
  implicit none
contains
  ! value of psi(L) for trial energy E
  function shoot_box(E, L, n) result(f)
    real(dp), intent(in) :: E, L
    integer,  intent(in) :: n
    real(dp) :: f, h, k2(0:n), psi(0:n)
    h = L/n
    k2 = 2.0_dp*E
    psi(0) = 0.0_dp; psi(1) = 1.0d-3*h
    call numerov(k2, psi, h, n)
    f = psi(n)
  end function shoot_box

  ! n-th eigenvalue + normalised eigenfunction
  subroutine box_eigen(nst, L, n, E, psi)
    integer,  intent(in)  :: nst, n
    real(dp), intent(in)  :: L
    real(dp), intent(out) :: E, psi(0:n)
    real(dp) :: Eex, Ea, Eb, a, b, fa, fb, fm, m, best, h, k2(0:n)
    integer  :: i, nscan = 200
    Eex = (nst*pi/L)**2/2.0_dp
    Ea = 0.5_dp*Eex*(1.0_dp - 1.0_dp/nst) + 1.0d-3
    Eb = Eex*(1.0_dp + 1.0_dp/nst)
    best = huge(1.0_dp)
    fa = shoot_box(Ea, L, n)
    do i = 1, nscan                      ! bracket: sign change nearest Eex
       b  = Ea + (Eb - Ea)*i/nscan
       fb = shoot_box(b, L, n)
       a  = Ea + (Eb - Ea)*(i-1)/nscan
       if (sign(1.0_dp,fa) /= sign(1.0_dp,fb)) then
          if (abs(a - Eex) < best) then
             best = abs(a - Eex); E = a
          end if
       end if
       fa = fb
    end do
    a = E; b = E + (Eb - Ea)/nscan; fa = shoot_box(a, L, n)
    do while (b - a > 1.0d-13*max(1.0_dp, abs(a)))      ! bisection
       m = 0.5_dp*(a + b); fm = shoot_box(m, L, n)
       if (sign(1.0_dp,fm) == sign(1.0_dp,fa)) then
          a = m; fa = fm
       else
          b = m
       end if
    end do
    E = 0.5_dp*(a + b)
    h = L/n; k2 = 2.0_dp*E
    psi(0) = 0.0_dp; psi(1) = 1.0d-3*h
    call numerov(k2, psi, h, n)
    psi = psi/sqrt(simpson(psi**2, h, n))               ! normalisation
    if (sum(psi*sin(nst*pi*[(i*h, i=0,n)]/L)) < 0.0_dp) psi = -psi
  end subroutine box_eigen
end module box_problem

program free_particle
  use numerov_tools
  use box_problem
  implicit none
  integer, parameter :: nA = 8, nB = 7
  integer  :: NsA(nA) = [25, 50, 100, 200, 400, 800, 1600, 3200]
  integer  :: NsB(nB) = [20, 40, 80, 160, 320, 640, 1280]
  real(dp) :: hA(nA), eA(nA), eT(nA), hB(nB), eB(3, nB)
  real(dp), allocatable :: x(:), psi(:), k2(:), ex(:), st(:,:)
  real(dp) :: k, X0, h, E, Eex, Eh, Eh2, ER
  integer  :: i, j, ng, s
  character(len=*), parameter :: fmtA = '(i6,f10.5,es18.3,f8.2,es18.3,f8.2)'

  !------------------------------ PART A
  k = 2.0_dp; X0 = 10.0_dp
  print '(a)', repeat('=',70)
  print '(a)', 'PART A : free particle, psi(0)=0, k = 2, exact solution sin(kx)'
  print '(a)', repeat('=',70)
  print '(a6,a10,a18,a8,a18,a8)', 'N', 'h', 'err(exact start)', 'ratio', 'err(Taylor start)', 'ratio'
  open(10, file='p1_partA.dat', status='replace')
  do i = 1, nA
     ng = NsA(i); h = X0/ng
     allocate(x(0:ng), psi(0:ng), k2(0:ng))
     x = [(j*h, j=0,ng)]; k2 = k*k
     psi(0) = 0.0_dp; psi(1) = sin(k*h)        ! exact start
     call numerov(k2, psi, h, ng)
     eA(i) = maxval(abs(psi - sin(k*x)))
     psi(0) = 0.0_dp; psi(1) = k*h             ! Taylor start
     call numerov(k2, psi, h, ng)
     eT(i) = maxval(abs(psi - sin(k*x)))
     hA(i) = h
     if (i == 1) then
        print fmtA, ng, h, eA(i), 0.0, eT(i), 0.0
     else
        print fmtA, ng, h, eA(i), eA(i-1)/eA(i), eT(i), eT(i-1)/eT(i)
     end if
     write(10,'(3es16.7)') h, eA(i), eT(i)
     deallocate(x, psi, k2)
  end do
  close(10)
  print '(a,f6.2,a)', 'observed order (exact start)  p = ', order(hA, eA, nA), '  (theory 4)'
  print '(a,f6.2,a)', 'observed order (Taylor start) p = ', order(hA, eT, nA), '  (theory 2)'

  !------------------------------ PART B
  print '(/a)', repeat('=',70)
  print '(a)', 'PART B : free particle in a box L = 1  (E_n = n^2 pi^2 / 2)'
  print '(a)', repeat('=',70)
  print '(a3,a18,a18,a11,a16,a15)', 'n', 'E_numerov', 'E_exact', 'rel.err', 'norm', 'max|psi-ex|'
  ng = 200; h = 1.0_dp/ng
  allocate(x(0:ng), psi(0:ng), ex(0:ng), st(0:ng,4))
  x = [(j*h, j=0,ng)]
  do s = 1, 4
     call box_eigen(s, 1.0_dp, ng, E, psi)
     Eex = (s*pi)**2/2.0_dp
     ex  = sqrt(2.0_dp)*sin(s*pi*x)
     st(:,s) = psi
     print '(i3,2f18.12,es11.3,f16.12,es15.3)', s, E, Eex, abs(E-Eex)/Eex, &
          simpson(psi**2,h,ng), maxval(abs(psi-ex))
  end do
  print '(a)', 'overlap matrix <psi_i|psi_j>  (should be identity):'
  do i = 1, 4
     write(*,'(4es11.2)') (simpson(st(:,i)*st(:,j),h,ng), j=1,4)
  end do
  open(11, file='p1_states.dat', status='replace')
  do j = 0, ng
     write(11,'(5es16.7)') x(j), st(j,1:4)
  end do
  close(11)
  deallocate(x, psi, ex, st)

  print '(/a)', 'Convergence of the eigenvalues in the box:'
  open(12, file='p1_partB.dat', status='replace')
  do i = 1, nB
     ng = NsB(i); hB(i) = 1.0_dp/ng
     allocate(psi(0:ng))
     do s = 1, 3
        call box_eigen(s, 1.0_dp, ng, E, psi)
        Eex = (s*pi)**2/2.0_dp
        eB(s,i) = abs(E-Eex)/Eex
     end do
     write(12,'(4es16.7)') hB(i), eB(:,i)
     deallocate(psi)
  end do
  close(12)
  do s = 1, 3
     write(*,'(a,i1,a,7es10.2,a,f5.2)') '  n=', s, ': ', eB(s,:), ' -> order', order(hB(1:5), eB(s,1:5), 5)
  end do

  ! Richardson extrapolation (error ~ h^4 -> weight 16/15, -1/15)
  allocate(psi(0:80))
  call box_eigen(1, 1.0_dp, 40, Eh, psi(0:40))
  call box_eigen(1, 1.0_dp, 80, Eh2, psi)
  ER = (16.0_dp*Eh2 - Eh)/15.0_dp
  print '(/a,f15.12,a,f15.12,a,f15.12,a,es9.2)', 'Richardson n=1: E(h)=', Eh, ' E(h/2)=', Eh2, &
       ' E_R=', ER, ' err=', abs(ER - pi**2/2)
end program free_particle
