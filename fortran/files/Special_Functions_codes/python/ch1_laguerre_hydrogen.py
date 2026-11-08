#!/usr/bin/env python3
"""
Chapter 1 -- LAGUERRE POLYNOMIALS and the HYDROGEN ATOM (Numerov)
=================================================================
Part A : computation of the generalised Laguerre polynomials L_n^(a)(x)
         * three-term recurrence  (stable)               DLMF §18.9
         * explicit finite sum    (cancellation!)          DLMF §18.5
         * Gauss-Laguerre nodes/weights: Newton on L_n^(a), and Golub-Welsch
         * verification: orthogonality by quadrature, mpmath reference
         * error analysis: sum vs recurrence, quadrature convergence
Part B : hydrogen atom, radial equation (hbar = m_e = e = 4 pi eps0 = 1, Z = 1)
         u'' = [ l(l+1)/r^2 - 2/r - 2E ] u ,  u = r R(r)
         exact  E_n = -1/(2n^2),  R_nl ~ e^{-r/n} (2r/n)^l L_{n-l-1}^{(2l+1)}(2r/n)
         Numerov, two-sided shooting, Frobenius start at r_s, Wronskian matching
         * verification: normalisation, comparison with the Laguerre-based
           exact R_nl (computed with Part A), <r>, virial <V> = 2E, orthogonality
         * error analysis: start at first grid nodes vs fixed r_s, h^4 law
Outputs: tables + ch1_functions.png, ch1_partA_errors.png, ch1_hydrogen.png,
         ch1_convergence.png
"""
import numpy as np
from math import lgamma, exp, factorial
import mpmath as mp
from scipy.linalg import eigh_tridiagonal
from scipy.special import eval_genlaguerre
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cd6"]
mp.mp.dps = 50


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


# ============================================================ PART A
def laguerre_rec(n, a, x):
    """L_0..L_n^(a)(x) by the three-term recurrence (returns array shape (n+1, len(x)))
       (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}"""
    x = np.asarray(x, float)
    L = np.zeros((n + 1,) + x.shape)
    L[0] = 1.0
    if n >= 1:
        L[1] = 1.0 + a - x
    for k in range(1, n):
        L[k + 1] = ((2 * k + 1 + a - x) * L[k] - (k + a) * L[k - 1]) / (k + 1)
    return L


def laguerre_sum(n, a, x):
    """explicit sum  L_n^(a)(x) = sum_k (-1)^k C(n+a, n-k) x^k / k!"""
    x = np.asarray(x, float)
    s = np.zeros_like(x)
    for k in range(n + 1):
        binom = exp(lgamma(n + a + 1) - lgamma(n - k + 1) - lgamma(a + k + 1))
        s += (-1) ** k * binom * x ** k / factorial(k)
    return s


def gauss_laguerre_newton(n, a=0.0, tol=1e-15):
    """nodes & weights of n-point Gauss-Laguerre rule: Newton on L_n^(a)
       with the Stroud-Secrest / Numerical Recipes initial guesses."""
    x = np.zeros(n); w = np.zeros(n)
    for i in range(n):
        if i == 0:
            z = (1 + a) * (3 + 0.92 * a) / (1 + 2.4 * n + 1.8 * a)
        elif i == 1:
            z = x[0] + (15 + 6.25 * a) / (1 + 0.9 * a + 2.5 * n)
        else:
            ai = i - 1
            z = x[i - 1] + ((1 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * a / (1 + 3.5 * ai)) \
                * (x[i - 1] - x[i - 2]) / (1 + 0.3 * a)
        for _ in range(100):
            L = laguerre_rec(n, a, np.array([z]))[:, 0]
            dL = (n * L[n] - (n + a) * L[n - 1]) / z       # derivative, DLMF §18.9
            dz = L[n] / dL
            z -= dz
            if abs(dz) < tol * max(1.0, z):
                break
        x[i] = z
        Lp1 = laguerre_rec(n + 1, a, np.array([z]))[n + 1, 0]
        w[i] = exp(lgamma(n + a + 1) - lgamma(n + 1)) * z / ((n + 1) ** 2 * Lp1 ** 2)
    return x, w


def gauss_laguerre_gw(n, a=0.0):
    """Golub-Welsch: eigenvalues of the symmetric Jacobi matrix"""
    k = np.arange(n)
    d = 2 * k + 1 + a
    e = -np.sqrt((k[1:]) * (k[1:] + a))
    x, V = eigh_tridiagonal(d, e)
    w = exp(lgamma(a + 1)) * V[0] ** 2
    return x, w


def partA():
    print("=" * 78)
    print("PART A  Generalised Laguerre polynomials L_n^(a)(x)")
    print("=" * 78)
    # --- 1. values vs high-precision reference
    xs = np.array([0.5, 5.0, 20.0, 60.0])
    print("Recurrence vs explicit sum vs mpmath (50 digits), n = 30, a = 1")
    print(f"{'x':>6} {'recurrence':>24} {'mpmath':>24} {'rel.err rec':>12} {'rel.err sum':>12}")
    for x in xs:
        ref = float(mp.laguerre(30, 1, x))
        r = laguerre_rec(30, 1, np.array([x]))[30, 0]
        s = laguerre_sum(30, 1, np.array([x]))[0]
        print(f"{x:6.1f} {r:24.15e} {ref:24.15e} {abs(r-ref)/abs(ref):12.2e} {abs(s-ref)/abs(ref):12.2e}")
    # --- 2. Gauss-Laguerre rules
    xn, wn = gauss_laguerre_newton(20)
    xg, wg = gauss_laguerre_gw(20)
    print(f"\n20-point Gauss-Laguerre: max|x_newton - x_GW| = {np.max(abs(xn-xg)):.2e}, "
          f"max|w_newton - w_GW|/w = {np.max(abs(wn-wg)/wg):.2e}, sum w = {wn.sum():.15f}")
    # --- 3. orthogonality by quadrature  int e^-x x^a L_m L_n = Gamma(n+a+1)/n! delta
    a = 2.0
    xq, wq = gauss_laguerre_newton(40, a)
    L = laguerre_rec(12, a, xq)
    G = (L * wq) @ L.T
    hn = np.array([exp(lgamma(n + a + 1) - lgamma(n + 1)) for n in range(13)])
    print(f"orthogonality (a=2, n,m<=12, 40-pt rule): max |G_mn/h_n - delta| = "
          f"{np.max(abs(G/np.sqrt(np.outer(hn,hn)) - np.eye(13))):.2e}")
    print(f"agreement with scipy eval_genlaguerre (n<=12, 40 nodes): "
          f"{np.max(abs(L - np.array([eval_genlaguerre(n, a, xq) for n in range(13)]))/(1+abs(L))):.2e}")

    # --- error analysis figure data
    xx = np.linspace(0.1, 100, 300)
    ref = np.array([float(mp.laguerre(30, 1, x)) for x in xx])
    e_rec = abs(laguerre_rec(30, 1, xx)[30] - ref) / abs(ref)
    e_sum = abs(laguerre_sum(30, 1, xx) - ref) / abs(ref)
    # quadrature convergence
    ns = np.arange(2, 41, 2)
    f1 = lambda x: np.cos(x)                    # int e^-x cos x = 1/2
    f2 = lambda x: 1 / (1 + x)                  # int e^-x/(1+x) = e E1(1)
    f3 = lambda x: np.sqrt(x)                   # int e^-x sqrt x = sqrt(pi)/2
    ex2 = float(mp.e * mp.e1(1))
    err = {k: [] for k in ("cos", "1/(1+x)", "sqrt")}
    for n in ns:
        xq, wq = gauss_laguerre_newton(n)
        err["cos"].append(abs(wq @ f1(xq) - 0.5))
        err["1/(1+x)"].append(abs(wq @ f2(xq) - ex2))
        err["sqrt"].append(abs(wq @ f3(xq) - np.sqrt(np.pi) / 2))
    print("\nGauss-Laguerre quadrature errors:")
    for i, n in enumerate(ns):
        if n in (2, 6, 10, 20, 30, 40):
            print(f"  n={n:2d}: cos {err['cos'][i]:.2e}   1/(1+x) {err['1/(1+x)'][i]:.2e}   sqrt(x) {err['sqrt'][i]:.2e}")
    sl = np.polyfit(np.log(ns[5:]), np.log(err["sqrt"][5:]), 1)[0]
    print(f"  sqrt(x): algebraic decay ~ n^{sl:.2f}")

    # ---- figures
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    x = np.linspace(0, 15, 600)
    L = laguerre_rec(5, 0, x)
    for n in range(6):
        ax[0].plot(x, L[n], color=C[n], lw=2, label=f"$L_{n}$")
    ax[0].set(ylim=(-10, 20), xlabel="x", title=r"(a) Laguerre polynomials $L_n(x)$ ($\alpha=0$)")
    ax[0].legend(frameon=False, fontsize=8, ncol=3)
    for n in range(4):
        y = np.exp(-x / 2) * laguerre_rec(n, 0, x)[n]
        ax[1].plot(x, y, color=C[n], lw=2, label=f"$e^{{-x/2}}L_{n}$")
    ax[1].set(xlabel="x", title="(b) Laguerre functions (orthonormal on $[0,\\infty)$)")
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch1_functions.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].semilogy(xx, e_sum + 1e-17, color=C[1], lw=1.5, label="explicit sum")
    ax[0].semilogy(xx, e_rec + 1e-17, color=C[0], lw=1.5, label="three-term recurrence")
    ax[0].set(xlabel="x", ylabel="relative error", title=r"(a) $L_{30}^{(1)}(x)$ vs 50-digit reference")
    ax[0].legend(frameon=False, fontsize=8)
    for i, (k, mk) in enumerate(zip(err, ("o-", "s-", "^-"))):
        ax[1].semilogy(ns, np.array(err[k]) + 1e-18, mk, color=C[i], lw=2, ms=4, label=fr"$\int_0^\infty e^{{-x}}\,{k}\,dx$" if k != "sqrt" else r"$\int_0^\infty e^{-x}\sqrt{x}\,dx$")
    ax[1].set(xlabel="number of nodes n", ylabel="absolute error", title="(b) Gauss–Laguerre quadrature convergence")
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch1_partA_errors.png", dpi=170)


# ============================================================ PART B
def R_exact(n, l, r):
    """normalised hydrogen radial function from the Laguerre recurrence (Part A)"""
    rho = 2 * r / n
    norm = np.sqrt((2.0 / n) ** 3 * exp(lgamma(n - l) - lgamma(n + l + 1)) / (2 * n))
    return norm * np.exp(-rho / 2) * rho ** l * laguerre_rec(n - l - 1, 2 * l + 1, rho)[n - l - 1]


def frobenius(r, l, E, nterm=40):
    """u(r) = r^{l+1} sum a_k r^k,  a_k k(k+2l+1) = -2 a_{k-1} - 2E a_{k-2}"""
    a = [1.0, -2.0 / (2 * l + 2)]
    for k in range(2, nterm):
        a.append((-2 * a[k - 1] - 2 * E * a[k - 2]) / (k * (k + 2 * l + 1)))
    s = np.zeros_like(r); p = np.ones_like(r)
    for ak in a:
        s += ak * p; p = p * r
    return r ** (l + 1) * s


def simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def shoot(E, l, h, R, rs):
    """returns D(E), assembled u, matching index.  rs = start radius (0 -> first nodes)"""
    N = int(round(R / h)); r = h * np.arange(N + 1)
    k2 = np.zeros(N + 1)
    k2[1:] = 2 * E + 2 / r[1:] - l * (l + 1) / r[1:] ** 2
    c = 1 + h * h * k2 / 12
    naive = rs < 0                                          # textbook start u0=0, u1=h^{l+1}
    js = 0 if naive else max(1, int(rs / h + 0.5))
    # outer classical turning point  E = -1/r + l(l+1)/(2r^2)
    rt = (1 + np.sqrt(max(1 + 2 * E * l * (l + 1), 0))) / (-2 * E)
    m = int(min(max(round(rt / h), js + 3), N - 3))
    pL = np.zeros(N + 1)
    if naive:
        pL[1] = h ** (l + 1)                                # c_0 u_0 taken as 0 (wrong for l<=1)
    else:
        pL[:js + 2] = frobenius(r[:js + 2], l, E)          # series up to the start nodes
    for n in range(js + 1, m + 1):
        pL[n + 1] = ((12 - 10 * c[n]) * pL[n] - c[n - 1] * pL[n - 1]) / c[n + 1]
    pR = np.zeros(N + 1); pR[N - 1] = 1e-30
    for n in range(N - 1, m - 1, -1):
        pR[n - 1] = ((12 - 10 * c[n]) * pR[n] - c[n + 1] * pR[n + 1]) / c[n - 1]
    sL = np.max(abs(pL[:m + 2])); sR = np.max(abs(pR[m - 1:]))
    pL /= sL; pR /= sR
    D = pL[m + 1] * pR[m] - pL[m] * pR[m + 1]
    u = np.concatenate([pL[:m + 1], pR[m + 1:] * pL[m] / pR[m]])
    return D, u, r


def eigen(n, l, h, R=None, rs=0.25):
    """bisection on D(E) in a bracket of +-15% around a scan-located root"""
    Ex = -0.5 / n ** 2
    R = 40.0 * n if R is None else R                        # box grows with n (u ~ r^n e^{-r/n})
    a, b = Ex * 1.25, Ex * 0.75
    fa = shoot(a, l, h, R, rs)[0]; fb = shoot(b, l, h, R, rs)[0]
    assert np.sign(fa) != np.sign(fb), (n, l, h)
    for _ in range(60):
        mid = 0.5 * (a + b); fm = shoot(mid, l, h, R, rs)[0]
        if np.sign(fm) == np.sign(fa): a, fa = mid, fm
        else: b = mid
    E = 0.5 * (a + b)
    _, u, r = shoot(E, l, h, R, rs)
    u /= np.sqrt(simpson(u * u, h))
    ex = r * R_exact(n, l, r)
    if np.dot(u, ex) < 0: u = -u
    return E, u, r


def scan_levels(l, h=0.01, R=160.0, rs=0.25):
    """locate all levels of given l with E < -0.02 by scanning D(E) on a log grid"""
    Eg = -np.exp(np.linspace(np.log(0.6), np.log(0.02), 500))
    D = np.array([shoot(E, l, h, R, rs)[0] for E in Eg])
    return [0.5 * (Eg[i] + Eg[i + 1]) for i in np.where(np.sign(D[:-1]) != np.sign(D[1:]))[0]]


def partB():
    print("\n" + "=" * 78)
    print("PART B  Hydrogen atom: radial Schroedinger equation, Numerov, h = 0.01, R = 40 n")
    print("=" * 78)
    for l in (0, 1, 2):
        lev = scan_levels(l)
        print(f"  scan l={l}: sign changes near E = " + ", ".join(f"{E:.4f}" for E in lev))
    print(f"\n{'n':>2} {'l':>2} {'E Numerov':>17} {'E exact':>13} {'|dE|':>9} {'norm':>15} {'max|u-u_ex|':>12}"
          f" {'<r>':>9} {'<r>exact':>9} {'<V>/2E':>10}")
    states = {}
    for (n, l) in ((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 3)):
        E, u, r = eigen(n, l, 0.01)
        h = r[1]
        ex = r * R_exact(n, l, r)
        rexp = simpson(u * u * r, h)
        Vexp = -simpson(np.concatenate([[0.0], u[1:] ** 2 / r[1:]]), h)        # <V> = -<1/r>
        states[n, l] = (u, r)
        print(f"{n:2d} {l:2d} {E:17.12f} {-0.5/n**2:13.10f} {abs(E+0.5/n**2):9.1e} {simpson(u*u,h):15.12f} "
              f"{np.max(abs(u-ex)):12.1e} {rexp:9.5f} {(3*n*n-l*(l+1))/2:9.5f} {Vexp/(2*E):10.7f}")
    for l in (0, 1):
        ns = [n for (n, ll) in states if ll == l]
        def ov(a, b):
            ua, ub = states[a, l][0], states[b, l][0]; k = min(len(ua), len(ub))
            return simpson(ua[:k - (k + 1) % 2] * ub[:k - (k + 1) % 2], 0.01)
        O = max(abs(ov(a, b)) for a in ns for b in ns if a < b)
        print(f"  orthogonality l={l}: max |<u_n|u_n'>| (n != n') = {O:.1e}")

    # ---------------- convergence: start strategies
    hs = np.array([0.08, 0.04, 0.02, 0.01, 0.005])
    cases = [(1, 0), (2, 1), (3, 2)]
    err = {}
    for (n, l) in cases:
        for rs, lab in ((-1.0, "naive"), (0.25, "Frobenius")):
            err[n, l, lab] = [abs(eigen(n, l, h, rs=rs)[0] + 0.5 / n ** 2) for h in hs]
    print("\nEigenvalue error vs h  (naive start u0=0,u1=h^(l+1)  |  Frobenius start)")
    print(f"{'h':>8}" + "".join(f"  {'(%d,%d) naive' % c:>14} {'(%d,%d) Frob.' % c:>13}" for c in cases))
    for i, h in enumerate(hs):
        print(f"{h:8.4f}" + "".join(f"  {err[c[0],c[1],'naive'][i]:14.2e} {err[c[0],c[1],'Frobenius'][i]:13.2e}" for c in cases))
    for c in cases:
        for lab in ("naive", "Frobenius"):
            e = np.array(err[c[0], c[1], lab]); ok = e > 2e-12
            p = np.polyfit(np.log(hs[ok][:5]), np.log(e[ok][:5]), 1)[0]
            print(f"  (n,l)={c}: {lab:12s} observed order {p:.2f}")
    # Richardson
    E1 = eigen(2, 1, 0.04)[0]; E2 = eigen(2, 1, 0.02)[0]
    print(f"\nRichardson (2p): E(0.04)={E1:.12f}  E(0.02)={E2:.12f}  E_R={(16*E2-E1)/15:.12f}  err={abs((16*E2-E1)/15+0.125):.1e}")
    # box truncation
    print("Box truncation (h=0.01), n=3,l=0:")
    for R in (25, 30, 40, 50, 60, 80, 100, 120):
        print(f"  R={R:3d}: |dE| = {abs(eigen(3, 0, 0.01, R=float(R))[0] + 0.5/9):.2e}")

    # ---------------- figures
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    for i, (n, l) in enumerate(((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2))):
        u, r = states[n, l]; sel = r < 25
        ax[0].plot(r[sel], u[sel] ** 2, color=C[i], lw=2, label=f"{n}{'spd'[l]}")
        ax[0].plot(r[sel][::60], (r * R_exact(n, l, r))[sel][::60] ** 2, "o", ms=2.5, color="#333333")
    ax[0].set(xlim=(0, 25))
    ax[0].set(xlabel="r (Bohr radii)", ylabel=r"$|u_{nl}(r)|^2=r^2|R_{nl}|^2$",
              title="(a) Radial densities (lines: Numerov, dots: Laguerre exact)")
    ax[0].legend(frameon=False, fontsize=8, ncol=2)
    for i, (n, l) in enumerate(((1, 0), (2, 1), (3, 2), (3, 0))):
        u, r = states[n, l]
        ax[1].semilogy(r, abs(u - r * R_exact(n, l, r)) + 1e-18, color=C[i], lw=1.2, label=f"{n}{'spd'[l]}")
    ax[1].set(xlim=(0, 60), ylim=(1e-15, 1e-8), xlabel="r", ylabel=r"$|u_{num}-u_{exact}|$", title="(b) Pointwise error, h = 0.01")
    ax[1].legend(frameon=False, fontsize=8, ncol=4)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch1_hydrogen.png", dpi=170)

    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    for i, c in enumerate(cases):
        ax.loglog(hs, err[c[0], c[1], "naive"], "s--", color=C[i], lw=1.5, ms=5, label=f"{c[0]}{'spd'[c[1]]}: naive start $u_1=h^{{l+1}}$")
        ax.loglog(hs, err[c[0], c[1], "Frobenius"], "o-", color=C[i], lw=2, ms=5, label=f"{c[0]}{'spd'[c[1]]}: Frobenius start")
    ax.loglog(hs, err[1, 0, "Frobenius"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax.loglog(hs, err[1, 0, "naive"][0] * (hs / hs[0]) ** 2, "-.", color="#888888", label=r"$\propto h^2$")
    ax.set(xlabel="h", ylabel=r"$|E_{nl}+1/2n^2|$", title="Hydrogen: convergence and the effect of the start")
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2); style(ax)
    fig.tight_layout(); fig.savefig("ch1_convergence.png", dpi=170)


if __name__ == "__main__":
    partA(); partB()
