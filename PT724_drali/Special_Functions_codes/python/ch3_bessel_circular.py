#!/usr/bin/env python3
"""
Chapter 3 -- BESSEL FUNCTIONS J_n (cylindrical) and the CIRCULAR INFINITE WELL
==============================================================================
Part A : computation of J_n(x)
         * power series                                  DLMF §10.2  (cancellation)
         * forward recurrence  J_{n+1} = (2n/x) J_n - J_{n-1}   DLMF §10.6 (unstable n > x)
         * Miller backward recurrence + normalisation
           J_0 + 2 sum J_{2k} = 1                        DLMF §3.6(iii), §10.12
         * trapezoid on the Bessel integral               DLMF §10.9 (geometric convergence)
         * zeros j_{n,s}: McMahon asymptotics + Newton   DLMF §10.21
         * verification: mpmath 50-digit reference, Jacobi-Anger identity
           cos x = J_0 - 2 J_2 + 2 J_4 - ...  (DLMF §10.12), scipy.special.jv
Part B : particle in a circular infinite well of radius 1 (hbar = m = 1)
         psi = R(r) e^{i m phi},  u = sqrt(r) R :
         u'' + [ 2E - (m^2 - 1/4)/r^2 ] u = 0 ,  u(1) = 0 ,  u ~ r^{m+1/2}
         exact  E_{m,s} = j_{m,s}^2 / 2 ,  R = J_m(j r) sqrt(2)/|J_{m+1}(j)|
         Numerov outward shooting with Frobenius start at r_s, bisection on u(1)
         * verification: normalisation int R^2 r dr = 1, comparison with J_m,
           orthogonality, energies vs zeros computed in Part A
         * error analysis: start at first nodes (order 2m, no convergence for m=0)
           vs fixed r_s (order 4); Richardson
Outputs: tables + ch3_functions.png, ch3_partA_errors.png, ch3_well.png, ch3_convergence.png
"""
import numpy as np
from math import factorial, lgamma, exp
import mpmath as mp
from scipy.special import jv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cd6"]
mp.mp.dps = 50


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


# ============================================================ PART A
def J_series(n, x, kmax=200):
    """J_n(x) = sum_k (-1)^k (x/2)^{2k+n} / (k! (k+n)!)"""
    x = float(x); s = 0.0
    term = (x / 2) ** n / factorial(n)
    for k in range(kmax):
        s += term
        term *= -(x / 2) ** 2 / ((k + 1) * (k + 1 + n))
        if abs(term) < 1e-17 * abs(s) and k > x: break
    return s


def J_forward(nmax, x, J0, J1):
    J = np.zeros(nmax + 1); J[0], J[1] = J0, J1
    for n in range(1, nmax):
        J[n + 1] = 2 * n / x * J[n] - J[n - 1]
    return J


def miller_start(nmax, x, digits=16):
    """starting index for Miller's algorithm (Numerical Recipes style + margin)"""
    return 2 * ((max(nmax, int(x)) + int(np.sqrt(digits * 4 * max(nmax, int(x), 1))) + 10) // 2)


def J_miller(nmax, x, M=None):
    """J_0..J_nmax(x) by backward recurrence from M (even), normalised with
       J_0 + 2 (J_2 + J_4 + ...) = 1.  Valid for x > 0."""
    if M is None: M = miller_start(nmax, x)
    M += M % 2
    f = np.zeros(M + 2); f[M + 1] = 0.0; f[M] = 1e-30
    for n in range(M, 0, -1):
        f[n - 1] = 2 * n / x * f[n] - f[n + 1]
        if abs(f[n - 1]) > 1e250:                     # rescale to avoid overflow
            f[n - 1:] *= 1e-250
    norm = f[0] + 2 * f[2:M + 1:2].sum()
    return f[:nmax + 1] / norm


def J_trapezoid(n, x, K):
    """J_n(x) = (1/pi) int_0^pi cos(n t - x sin t) dt with K-panel trapezoid"""
    t = np.linspace(0, np.pi, K + 1)
    f = np.cos(n * t - x * np.sin(t))
    return (f.sum() - 0.5 * (f[0] + f[-1])) / K


def J(n, x):
    """recommended: Miller for x>0, J_n(0)=delta"""
    if x == 0: return 1.0 if n == 0 else 0.0
    return J_miller(n + 1, x)[n]


def dJ(n, x):
    Jm = J_miller(n + 2, x)
    return (Jm[n - 1] - Jm[n + 1]) / 2 if n > 0 else -Jm[1]


def bessel_zeros(n, S):
    """first S positive zeros of J_n: McMahon start, Newton refinement"""
    z = []
    mu = 4 * n * n
    for s in range(1, S + 1):
        b = (s + n / 2 - 0.25) * np.pi
        x = b - (mu - 1) / (8 * b) - 4 * (mu - 1) * (7 * mu - 31) / (3 * (8 * b) ** 3)
        for _ in range(50):
            Jm = J_miller(n + 2, x)
            d = (Jm[n - 1] - Jm[n + 1]) / 2 if n > 0 else -Jm[1]
            dx = Jm[n] / d; x -= dx
            if abs(dx) < 1e-15 * x: break
        z.append(x)
    return np.array(z)


def partA():
    print("=" * 78)
    print("PART A  Bessel functions of the first kind J_n(x)")
    print("=" * 78)
    print("J_n(x): series / Miller / mpmath")
    print(f"{'n':>3} {'x':>6} {'Miller':>24} {'mpmath':>24} {'rel.err Miller':>15} {'rel.err series':>15}")
    for n, x in ((0, 1.0), (5, 10.0), (0, 30.0), (3, 50.0), (20, 5.0), (50, 1.0)):
        ref = float(mp.besselj(n, x))
        jm = J(n, x); js = J_series(n, x)
        print(f"{n:3d} {x:6.1f} {jm:24.15e} {ref:24.15e} {abs(jm-ref)/abs(ref):15.2e} {abs(js-ref)/abs(ref):15.2e}")
    # forward recurrence instability at x = 1
    x = 1.0
    Jm = J_miller(30, x)
    Jf = J_forward(30, x, Jm[0], Jm[1])
    ref = np.array([float(mp.besselj(n, x)) for n in range(31)])
    print("\nForward vs backward recurrence at x = 1 (relative error):")
    for n in (2, 5, 10, 15, 20, 25, 30):
        print(f"  n={n:2d}: forward {abs(Jf[n]-ref[n])/abs(ref[n]):9.2e}   Miller {abs(Jm[n]-ref[n])/abs(ref[n]):9.2e}   J_n = {ref[n]:.3e}")
    # Jacobi-Anger identity check
    for x in (0.5, 7.3, 25.0):
        Jm = J_miller(80, x)
        s = Jm[0] + 2 * sum((-1) ** k * Jm[2 * k] for k in range(1, 40))
        print(f"  Jacobi-Anger x={x:5.1f}: J0 - 2J2 + 2J4 - ... - cos x = {s-np.cos(x):.1e}")
    xs = np.linspace(0.1, 40, 50)
    print(f"  max |J_miller - scipy jv| for n<=10, 0.1<x<40: "
          f"{max(np.max(abs(J_miller(10, x) - jv(np.arange(11), x))) for x in xs):.1e}")
    z = {n: bessel_zeros(n, 5) for n in range(4)}
    for n in range(4):
        refz = [float(mp.besseljzero(n, s)) for s in range(1, 6)]
        print(f"  zeros J_{n}: " + ", ".join(f"{v:.12f}" for v in z[n]) + f"   max err {max(abs(z[n]-refz)):.1e}")

    # error-analysis data
    xx = np.linspace(0.5, 60, 240)
    ref0 = np.array([float(mp.besselj(0, x)) for x in xx])
    e_ser = [abs(J_series(0, x) - r) for x, r in zip(xx, ref0)]
    e_mil = [abs(J(0, x) - r) for x, r in zip(xx, ref0)]
    Ms = np.arange(2, 62, 2)
    x = 10.0; r5 = float(mp.besselj(5, x))
    e_M = [abs(J_miller(5, x, M=M)[5] - r5) / abs(r5) if M > 5 else np.nan for M in Ms]
    Ks = np.arange(2, 41, 1)
    e_K = [abs(J_trapezoid(5, x, K) - r5) / abs(r5) + 1e-18 for K in Ks]
    nn = np.arange(31)
    e_f = abs(Jf - ref) / abs(ref); e_b = abs(J_miller(30, 1.0) - ref) / abs(ref)
    print("\nMiller start index M (x=10, n=5): " + "  ".join(f"M={M}:{e:.1e}" for M, e in zip(Ms, e_M) if M in (10, 16, 20, 26, 30, 40)))
    print("Trapezoid panels K (x=10, n=5):   " + "  ".join(f"K={K}:{e:.1e}" for K, e in zip(Ks, e_K) if K in (5, 10, 15, 20, 25, 30)))

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    xg = np.linspace(0.001, 20, 600)
    for n in range(4):
        ax[0].plot(xg, [J(n, v) for v in xg], color=C[n], lw=2, label=f"$J_{n}$")
        ax[0].plot(z[n][:3], 0 * z[n][:3], "o", ms=4, color=C[n])
    ax[0].axhline(0, color="#888888", lw=.8)
    ax[0].set(xlabel="x", title=r"(a) $J_n(x)$ and its first zeros $j_{n,s}$")
    ax[0].legend(frameon=False, fontsize=8, ncol=4)
    from scipy.special import yv
    for n in range(3):
        ax[1].plot(xg[30:], yv(n, xg[30:]), color=C[n], lw=2, label=f"$Y_{n}$")
    ax[1].set(ylim=(-2, 0.8), xlabel="x", title=r"(b) Second solution $Y_n(x)$ (singular at 0)")
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch3_functions.png", dpi=170)

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.7))
    ax[0].semilogy(xx, np.array(e_ser) + 1e-18, color=C[1], lw=1.5, label="power series")
    ax[0].semilogy(xx, np.array(e_mil) + 1e-18, color=C[0], lw=1.5, label="Miller")
    ax[0].set(xlabel="x", ylabel="absolute error", title=r"(a) $J_0(x)$: series cancellation")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].semilogy(nn[1:], e_f[1:] + 1e-18, "s-", color=C[1], lw=1.5, ms=4, label="forward recurrence")
    ax[1].semilogy(nn[1:], e_b[1:] + 1e-18, "o-", color=C[0], lw=1.5, ms=4, label="backward (Miller)")
    ax[1].set(xlabel="order n", ylabel="relative error", title="(b) Recurrence direction, x = 1")
    ax[1].legend(frameon=False, fontsize=8)
    ax[2].semilogy(Ms, e_M, "o-", color=C[0], lw=1.5, ms=4, label="Miller, start index M")
    ax[2].semilogy(Ks, e_K, "s-", color=C[2], lw=1.5, ms=4, label="trapezoid, K panels")
    ax[2].set(xlabel="M or K", ylabel="relative error", title=r"(c) Convergence for $J_5(10)$")
    ax[2].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch3_partA_errors.png", dpi=170)
    return z


# ============================================================ PART B
def frobenius(r, m, E, nterm=40):
    """u = r^{m+1/2} sum a_k r^{2k},  2k(2k+2m) a_k = -2E a_{k-1}"""
    r = np.asarray(r, float)
    s = np.zeros_like(r); a = 1.0; p = np.ones_like(r)
    for k in range(nterm):
        s += a * p
        a = -2 * E * a / ((2 * k + 2) * (2 * k + 2 + 2 * m)); p = p * r * r
    return r ** (m + 0.5) * s


def simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def shoot(E, m, N, rs):
    h = 1.0 / N; r = h * np.arange(N + 1)
    k2 = np.zeros(N + 1); k2[1:] = 2 * E - (m * m - 0.25) / r[1:] ** 2
    c = 1 + h * h * k2 / 12
    js = max(1, int(rs / h + 0.5))
    u = np.zeros(N + 1); u[:js + 2] = frobenius(r[:js + 2], m, E)
    for n in range(js + 1, N):
        u[n + 1] = ((12 - 10 * c[n]) * u[n] - c[n - 1] * u[n - 1]) / c[n + 1]
    return u, r


def eigen(m, s, N, rs=0.125, zeros=None):
    Ex = zeros[m][s - 1] ** 2 / 2
    a, b = Ex * 0.85, Ex * 1.15
    fa = shoot(a, m, N, rs)[0][-1]
    for _ in range(64):
        mid = 0.5 * (a + b); fm = shoot(mid, m, N, rs)[0][-1]
        if np.sign(fm) == np.sign(fa): a, fa = mid, fm
        else: b = mid
    E = 0.5 * (a + b)
    u, r = shoot(E, m, N, rs)
    u /= np.sqrt(simpson(u * u, r[1]))
    j = zeros[m][s - 1]
    ex = np.sqrt(r) * np.array([J(m, j * v) for v in r]) * np.sqrt(2) / abs(J(m + 1, j))
    if np.dot(u, ex) < 0: u = -u
    return E, u, r, ex


def partB(z):
    print("\n" + "=" * 78)
    print("PART B  Circular infinite well, radius 1:  Numerov in u = sqrt(r) R,  N = 400, r_s = 0.125")
    print("=" * 78)
    # scan, m = 0
    Eg = np.linspace(0.5, 80, 800)
    F = np.array([shoot(E, 0, 400, 0.125)[0][-1] for E in Eg])
    roots = [0.5 * (Eg[i] + Eg[i + 1]) for i in np.where(np.sign(F[:-1]) != np.sign(F[1:]))[0]]
    print("  scan m=0: sign changes of u(1;E) near " + ", ".join(f"{v:.2f}" for v in roots)
          + "   (exact j_{0,s}^2/2 = " + ", ".join(f"{v**2/2:.2f}" for v in z[0][:4]) + ")")
    print(f"\n{'m':>2} {'s':>2} {'E Numerov':>17} {'j^2/2 (Part A)':>17} {'|dE|':>9} {'norm':>15} {'max|u-u_ex|':>12} {'<r^2>':>10} {'exact':>10}")
    store = {}
    for (m, s) in ((0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (2, 1), (3, 1)):
        E, u, r, ex = eigen(m, s, 400, zeros=z)
        Ex = z[m][s - 1] ** 2 / 2
        store[m, s] = (u, r, ex)
        print(f"{m:2d} {s:2d} {E:17.12f} {Ex:17.12f} {abs(E-Ex):9.1e} {simpson(u*u,r[1]):15.12f} {np.max(abs(u-ex)):12.1e} {simpson(u*u*r*r, r[1]):10.6f} {1/3+2*(m*m-1)/(3*z[m][s-1]**2):10.6f}")
    ov = max(abs(simpson(store[0, a][0] * store[0, b][0], 1 / 400)) for a in (1, 2, 3) for b in (1, 2, 3) if a < b)
    print(f"  orthogonality (m=0, s=1,2,3): max |<u_s|u_s'>| = {ov:.1e}")

    Ns = np.array([20, 40, 80, 160, 320, 640, 1280])
    hs = 1.0 / Ns
    err = {}
    for m in (0, 1, 2):
        Ex = z[m][0] ** 2 / 2
        for rs, lab in ((0.0, "first nodes"), (0.125, "r_s")):
            err[m, lab] = [abs(eigen(m, 1, N, rs, z)[0] - Ex) for N in Ns]
    print("\nError of E_{m,1} vs h: Frobenius start at first nodes  vs  at fixed r_s = 0.125")
    print(f"{'N':>6} {'h':>9}" + "".join(f" {'m=%d nodes' % m:>12} {'m=%d r_s' % m:>11}" for m in (0, 1, 2)))
    for i, N in enumerate(Ns):
        print(f"{N:6d} {hs[i]:9.6f}" + "".join(f" {err[m,'first nodes'][i]:12.2e} {err[m,'r_s'][i]:11.2e}" for m in (0, 1, 2)))
    for m in (0, 1, 2):
        for lab in ("first nodes", "r_s"):
            e = np.array(err[m, lab]); ok = e > 1e-12
            print(f"  m={m} {lab:12s}: observed order {np.polyfit(np.log(hs[ok][1:6]), np.log(e[ok][1:6]), 1)[0]:.2f}")
    E1 = eigen(0, 1, 160, 0.125, z)[0]; E2 = eigen(0, 1, 320, 0.125, z)[0]
    ER = (16 * E2 - E1) / 15
    print(f"Richardson m=0: E(1/160)={E1:.12f} E(1/320)={E2:.12f} E_R={ER:.12f} err={abs(ER - z[0][0]**2/2):.1e}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    for i, (m, s) in enumerate(((0, 1), (0, 2), (1, 1), (2, 1), (1, 2))):
        u, r, ex = store[m, s]
        R = np.zeros_like(u); R[1:] = u[1:] / np.sqrt(r[1:]); R[0] = R[1] if m == 0 else 0
        ax[0].plot(r, R, color=C[i], lw=2, label=f"m={m}, s={s}")
        ax[0].plot(r[1::16], ex[1::16] / np.sqrt(r[1::16]), "o", ms=2.5, color="#333333")
    ax[0].set(xlabel="r", ylabel="R(r)", title=r"(a) Radial functions (dots: $\sqrt{2}J_m(jr)/|J_{m+1}(j)|$)")
    ax[0].legend(frameon=False, fontsize=7, ncol=2)
    # 2D density of (m=1,s=2) cos(phi) standing wave
    xg = np.linspace(-1, 1, 241); X, Y = np.meshgrid(xg, xg)
    Rr = np.hypot(X, Y); Ph = np.arctan2(Y, X)
    j = z[1][1]
    psi = np.where(Rr <= 1, jv(1, j * Rr) * np.cos(Ph), np.nan)
    im = ax[1].imshow(psi ** 2, extent=(-1, 1, -1, 1), origin="lower", cmap="Blues")
    ax[1].set(title=r"(b) $|\psi|^2$ for m=1, s=2  ($\cos\phi$ standing wave)", xlabel="x", ylabel="y")
    style(ax[0])
    fig.tight_layout(); fig.savefig("ch3_well.png", dpi=170)

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for i, m in enumerate((0, 1, 2)):
        ax.loglog(hs, err[m, "first nodes"], "s--", color=C[i], lw=1.5, ms=5, label=f"m={m}: start at first nodes")
        ax.loglog(hs, err[m, "r_s"], "o-", color=C[i], lw=2, ms=5, label=fr"m={m}: Frobenius start at $r_s$=1/8")
    ax.loglog(hs, err[0, "r_s"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax.loglog(hs, err[1, "first nodes"][0] * (hs / hs[0]) ** 2, "-.", color="#888888", label=r"$\propto h^2$")
    ax.set(xlabel="h = 1/N", ylabel=r"$|E_{m,1}-j_{m,1}^2/2|$", title="Circular well: the start fault (order $\\min(4,2m)$)\nand its cure")
    ax.legend(frameon=False, fontsize=7); style(ax)
    fig.tight_layout(); fig.savefig("ch3_convergence.png", dpi=170)


if __name__ == "__main__":
    z = partA(); partB(z)
