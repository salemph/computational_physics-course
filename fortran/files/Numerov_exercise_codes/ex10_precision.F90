!=======================================================================
! Exercise 10 -- ROUND-OFF: harmonic-oscillator ground state E_0 in
! single precision (real(4)) versus double precision (real(8)).
! The same source is compiled twice with the C preprocessor:
!   gfortran -O2 -cpp -DWP=4 -o ex10_sp ex10_precision.F90
!   gfortran -O2 -cpp -DWP=8 -o ex10_dp ex10_precision.F90
! Method: two-sided Numerov shooting with the discrete-Wronskian
! matching function D(E) (identical to Problem 3), box [-L, L], L = 6.
! Output: h, |E_0 - 0.5| for N = 12 ... 76800  (file ex10_wpX.dat)
!=======================================================================
module prec
  implicit none
  integer, parameter :: wp = WP
end module prec

module ho_shoot
  use prec
  implicit none
contains
  function Dfun(E, L, N) result(D)
    real(wp), intent(in) :: E, L
    integer,  intent(in) :: N
    real(wp) :: D, h, x, c(0:N), pL(0:N), pR(0:N), sL, sR
    integer  :: j, m
    h = 2.0_wp*L/N
    do j = 0, N
       x = -L + j*h
       c(j) = 1.0_wp + h*h*(2.0_wp*E - x*x)/12.0_wp
    end do
    m = nint((sqrt(2.0_wp*E) + L)/h)          ! node at the right turning point
    m = max(2, min(m, N-2))
    pL = 0.0_wp; pR = 0.0_wp
    pL(1) = 1.0e-10_wp
    do j = 1, m
       pL(j+1) = ((12.0_wp - 10.0_wp*c(j))*pL(j) - c(j-1)*pL(j-1))/c(j+1)
    end do
    pR(N-1) = 1.0e-10_wp
    do j = N-1, m, -1
       pR(j-1) = ((12.0_wp - 10.0_wp*c(j))*pR(j) - c(j+1)*pR(j+1))/c(j-1)
    end do
    sL = maxval(abs(pL(0:m+1))); sR = maxval(abs(pR(m-1:N)))
    D = (pL(m+1)/sL)*(pR(m)/sR) - (pL(m)/sL)*(pR(m+1)/sR)
  end function Dfun

  function ground_state(L, N) result(E)
    real(wp), intent(in) :: L
    integer,  intent(in) :: N
    real(wp) :: E, a, b, fa, fm, mid
    integer  :: it
    a = 0.3_wp; b = 0.7_wp                    ! brackets only E_0 = 0.5
    fa = Dfun(a, L, N)
    do it = 1, 200                            ! bisection down to machine precision
       mid = 0.5_wp*(a + b)
       if (mid <= a .or. mid >= b) exit
       fm = Dfun(mid, L, N)
       if (sign(1.0_wp, fm) == sign(1.0_wp, fa)) then
          a = mid; fa = fm
       else
          b = mid
       end if
    end do
    E = 0.5_wp*(a + b)
  end function ground_state
end module ho_shoot

program ex10
  use prec
  use ho_shoot
  implicit none
  integer, parameter :: nN = 13
  integer  :: Ns(nN) = [12, 24, 48, 96, 192, 384, 768, 1536, 3072, 6144, 12288, 24576, 49152]
  real(wp) :: L, E, h
  integer  :: i
  character(len=20) :: fname
  L = 6.0_wp
  write(fname,'(a,i1,a)') 'ex10_wp', wp, '.dat'
  open(10, file=fname, status='replace')
  print '(a,i2,a,es10.3)', 'working precision: real(', wp, ')   epsilon = ', epsilon(1.0_wp)
  print '(a8,a12,a22,a14)', 'N', 'h', 'E_0', '|E_0 - 0.5|'
  do i = 1, nN
     h = 2.0_wp*L/Ns(i)
     E = ground_state(L, Ns(i))
     print '(i8,es12.4,f22.16,es14.4)', Ns(i), h, E, abs(E - 0.5_wp)
     write(10,'(2es16.7)') real(h,8), real(abs(E - 0.5_wp),8)
  end do
  close(10)
end program ex10
