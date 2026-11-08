!=======================================================================
! Problem 2 -- RECTANGULAR POTENTIAL BARRIER with the Numerov method
! Units hbar = m = 1.  V = V0 for 0 < x < a, 0 elsewhere.
! psi'' = -k^2(x) psi,  k^2 = 2(E - V)
!
! Backward integration of the COMPLEX Numerov recursion from the right,
! where psi = exp(ikx) (transmitted wave, F = 1), to the left, where
! psi = A exp(ikx) + B exp(-ikx).   T = 1/|A|^2,  R = |B|^2/|A|^2.
! Checks: R + T = 1 (flux conservation) and T vs analytic formula.
! Three treatments of the potential discontinuity (mode):
!   1 = misaligned : barrier edges inside grid cells
!   2 = aligned    : edges on grid nodes, V(edge) = V0/2
!   3 = corrected  : aligned + one-sided V in neighbouring stencils
!                    + jump term delta = -h^3/12 [k^2] psi'(edge)
! Output : screen tables + p2_T.dat, p2_conv.dat, p2_wave.dat
! Compile: gfortran -O2 -o p2 p2_square_barrier.f90 && ./p2
!=======================================================================
module barrier
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: V0 = 5.0_dp, aw = 1.0_dp, Lpad = 2.0_dp
  complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
contains
  function T_exact(E) result(T)
    real(dp), intent(in) :: E
    real(dp) :: T, kap, q
    if (abs(E - V0) < 1.0d-12) then
       T = 1.0_dp/(1.0_dp + V0*aw*aw/2.0_dp)
    else if (E < V0) then
       kap = sqrt(2.0_dp*(V0 - E))
       T = 1.0_dp/(1.0_dp + V0**2*sinh(kap*aw)**2/(4.0_dp*E*(V0 - E)))
    else
       q = sqrt(2.0_dp*(E - V0))
       T = 1.0_dp/(1.0_dp + V0**2*sin(q*aw)**2/(4.0_dp*E*(E - V0)))
    end if
  end function T_exact

  subroutine scatter(E, h, mode, T, R, wavefile)
    real(dp), intent(in)  :: E, h
    integer,  intent(in)  :: mode
    real(dp), intent(out) :: T, R
    character(len=*), intent(in), optional :: wavefile
    real(dp), allocatable :: x(:), Vx(:), k2(:), c(:), jump(:)
    logical,  allocatable :: edge(:)
    complex(dp), allocatable :: psi(:)
    complex(dp) :: rhs, dpsi, M(2,2), det, Ac, Bc
    real(dp) :: k, xl, xr, shift, cp, cm
    integer  :: N, i, j

    k = sqrt(2.0_dp*E)
    shift = 0.0_dp
    if (mode == 1) shift = 0.37_dp*h
    xl = -Lpad + shift; xr = aw + Lpad + shift
    N = nint((xr - xl)/h)
    allocate(x(0:N), Vx(0:N), k2(0:N), c(0:N), psi(0:N), edge(0:N), jump(0:N))
    x = [(xl + i*h, i=0,N)]
    Vx = 0.0_dp; edge = .false.; jump = 0.0_dp
    where (x > 0.0_dp .and. x < aw) Vx = V0
    if (mode /= 1) then                       ! edge nodes get the mean value
       do i = 0, N
          if (abs(x(i)) < 1.0d-3*h .or. abs(x(i) - aw) < 1.0d-3*h) then
             Vx(i) = V0/2.0_dp; edge(i) = .true.
          end if
       end do
    end if
    k2 = 2.0_dp*(E - Vx)
    c  = 1.0_dp + h*h*k2/12.0_dp
    if (mode == 3) then
       do i = 1, N-1
          if (edge(i)) jump(i) = 2.0_dp*(E - Vx(i+1)) - 2.0_dp*(E - Vx(i-1))
       end do
    end if

    psi(N)   = exp(iu*k*x(N))                 ! transmitted wave, F = 1
    psi(N-1) = exp(iu*k*x(N-1))
    do j = N-1, 1, -1                         ! backward Numerov
       cp = c(j+1); cm = c(j-1)
       if (mode == 3 .and. .not. edge(j)) then
          if (edge(j+1)) cp = 1.0_dp + h*h*2.0_dp*(E - Vx(j))/12.0_dp
          if (edge(j-1)) cm = 1.0_dp + h*h*2.0_dp*(E - Vx(j))/12.0_dp
       end if
       rhs = (12.0_dp - 10.0_dp*c(j))*psi(j) - cp*psi(j+1)
       if (mode == 3 .and. edge(j)) then
          dpsi = (-3.0_dp*psi(j) + 4.0_dp*psi(j+1) - psi(j+2))/(2.0_dp*h)
          rhs  = rhs - h**3/12.0_dp*jump(j)*dpsi
       end if
       psi(j-1) = rhs/cm
    end do

    ! solve [e^{ikx0} e^{-ikx0}; e^{ikx1} e^{-ikx1}] (A,B) = (psi0, psi1)
    M(1,1) = exp(iu*k*x(0)); M(1,2) = exp(-iu*k*x(0))
    M(2,1) = exp(iu*k*x(1)); M(2,2) = exp(-iu*k*x(1))
    det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
    Ac  = ( M(2,2)*psi(0) - M(1,2)*psi(1))/det
    Bc  = (-M(2,1)*psi(0) + M(1,1)*psi(1))/det
    T = 1.0_dp/abs(Ac)**2
    R = abs(Bc)**2/abs(Ac)**2

    if (present(wavefile)) then
       open(20, file=wavefile, status='replace')
       do i = 0, N
          write(20,'(2es16.7)') x(i), abs(psi(i)/Ac)**2
       end do
       close(20)
    end if
    deallocate(x, Vx, k2, c, psi, edge, jump)
  end subroutine scatter

  function order(h, e, m) result(p)
    integer,  intent(in) :: m
    real(dp), intent(in) :: h(m), e(m)
    real(dp) :: p, lx(m), ly(m)
    lx = log(h); ly = log(e)
    p = (m*sum(lx*ly) - sum(lx)*sum(ly)) / (m*sum(lx*lx) - sum(lx)**2)
  end function order
end module barrier

program square_barrier
  use barrier
  implicit none
  real(dp) :: Elist(11) = [0.5_dp, 1.0_dp, 2.5_dp, 4.0_dp, 4.9_dp, 5.1_dp, 6.0_dp, &
                           8.0_dp, 9.935_dp, 12.0_dp, 15.0_dp]
  real(dp) :: hs(7) = [0.1_dp, 0.05_dp, 0.025_dp, 0.0125_dp, 0.00625_dp, 0.003125_dp, 0.0015625_dp]
  real(dp) :: Es(2) = [2.0_dp, 8.0_dp]
  real(dp) :: err(3,7), fl(7), T, R, E
  integer  :: i, ie, md

  print '(a)', repeat('=',74)
  print '(a)', 'Square barrier  V0 = 5,  a = 1   (hbar = m = 1),  h = 0.01, corrected scheme'
  print '(a)', repeat('=',74)
  print '(a6,2a15,a11,a13,a12)', 'E', 'T_numerov', 'T_exact', '|dT|', 'R', 'R+T-1'
  do i = 1, 11
     call scatter(Elist(i), 0.01_dp, 3, T, R)
     print '(f6.3,2f15.10,es11.2,f13.8,es12.2)', Elist(i), T, T_exact(Elist(i)), &
          abs(T - T_exact(Elist(i))), R, R + T - 1.0_dp
  end do

  open(11, file='p2_T.dat', status='replace')
  do i = 1, 400
     E = 0.05_dp + (20.0_dp - 0.05_dp)*(i-1)/399.0_dp
     call scatter(E, 0.01_dp, 3, T, R)
     write(11,'(3es16.7)') E, T, T_exact(E)
  end do
  close(11)

  open(12, file='p2_conv.dat', status='replace')
  do ie = 1, 2
     do i = 1, 7
        do md = 1, 3
           call scatter(Es(ie), hs(i), md, T, R)
           err(md,i) = abs(T - T_exact(Es(ie)))
           if (md == 3) fl(i) = abs(R + T - 1.0_dp)
        end do
        write(12,'(f6.2,5es16.7)') Es(ie), hs(i), err(:,i), fl(i)
     end do
     print '(/a,f4.1,a,f15.12)', 'Convergence of T at E = ', Es(ie), '   T_exact = ', T_exact(Es(ie))
     print '(a10,a13,a14,a13,a15)', 'h', 'misaligned', 'aligned V0/2', 'corrected', '|R+T-1| corr.'
     do i = 1, 7
        print '(f10.6,es13.2,es14.2,es13.2,es15.1)', hs(i), err(:,i), fl(i)
     end do
     print '(a,3f6.2)', 'observed order (misaligned, aligned, corrected): ', &
          order(hs(1:5), err(1,1:5), 5), order(hs(1:5), err(2,1:5), 5), order(hs(1:5), err(3,1:5), 5)
  end do
  close(12)

  call scatter(2.0_dp, 0.005_dp, 3, T, R, 'p2_wave_E2.dat')
  call scatter(8.0_dp, 0.005_dp, 3, T, R, 'p2_wave_E8.dat')
end program square_barrier
