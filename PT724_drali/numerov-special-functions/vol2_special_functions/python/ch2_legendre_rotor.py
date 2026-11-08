#!/usr/bin/env python3
"""
Chapter 2 -- LEGENDRE FUNCTIONS and the ANGULAR SCHROEDINGER EQUATION (Numerov)
===============================================================================
Part A : Legendre polynomials P_l(x) and associated functions P_l^m(x)
         * Bonnet recurrence for P_l                      DLMF §14.10, §18.9
         * normalised associated functions  Pbar_l^m  by the stable
           recurrence in l (no factorials, no overflow)    DLMF §14.10, §14.7
         * explicit sum for P_l                            DLMF §18.5 (cancellation)
         * Gauss-Legendre nodes and weights (Newton)       DLMF §3.5(v)
         * verification: orthonormality, mpmath reference, scipy lpmv
         * error analysis: sum vs recurrence, unnormalised overflow,
           Gauss-Legendre vs trapezoid convergence
Part B : angular equation (rigid rotor, hbar = I = 1, E = lambda/2)
         -(1/sin t) d/dt (sin t dTheta/dt) + m^2/sin^2 t Theta = lambda Theta
         Liouville transform  u = sqrt(sin t) Theta :
         u'' + [ lambda + 1/4 + (1/4 - m^2)/sin^2 t ] u = 0 ,  u ~ t^{m+1/2}
         exact lambda = l(l+1), Theta = Pbar_l^m(cos t)
         Numerov, two-sided shooting from both poles, Frobenius start,
         Wronskian matching at t = pi/2
         * verification: normalisation int u^2 dt = 1, comparison with Pbar_l^m,
           <cos^2 t>, parity, orthogonality
         * error analysis: start at first nodes (order 2m!) vs fixed t_s (order 4)
Outputs: tables + ch2_functions.png, ch2_partA_errors.png, ch2_rotor.png, ch2_convergence.png
"""
import numpy as np
from math import lgamma, exp, comb
import mpmath as mp
from scipy.special import lpmv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cd6"]
mp.mp.dps = 50


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


# ============================================================ PART A
def legendre_rec(L, x):
    """P_0..P_L(x) by Bonnet:  (l+1) P_{l+1} = (2l+1) x P_l - l P_{l-1}"""
    x = np.asarray(x, float)
    P = np.zeros((L + 1,) + x.shape); P[0] = 1.0
    if L >= 1: P[1] = x
    for l in range(1, L):
        P[l + 1] = ((2 * l + 1) * x * P[l] - l * P[l - 1]) / (l + 1)
    return P


def legendre_sum(l, x):
    """explicit sum  P_l = 2^-l sum_k (-1)^k C(l,k) C(2l-2k,l) x^{l-2k}"""
    x = np.asarray(x, float)
    return sum((-1) ** k * comb(l, k) * comb(2 * l - 2 * k, l) * x ** (l - 2 * k)
               for k in range(l // 2 + 1)) / 2.0 ** l


def plm_bar(L, m, x):
    """orthonormal associated Legendre functions Pbar_l^m(x), l = m..L
       (int_{-1}^{1} Pbar^2 dx = 1, Condon-Shortley phase included)
       Pbar_m^m     = -sqrt((2m+1)/(2m)) sqrt(1-x^2) Pbar_{m-1}^{m-1},  Pbar_0^0 = 1/sqrt 2
       Pbar_{m+1}^m = sqrt(2m+3) x Pbar_m^m
       Pbar_l^m     = a_l ( x Pbar_{l-1}^m - Pbar_{l-2}^m / a_{l-1} ),
                      a_l = sqrt((4l^2-1)/(l^2-m^2))"""
    x = np.asarray(x, float)
    s = np.sqrt(np.maximum(1 - x * x, 0.0))
    p = np.full_like(x, 1 / np.sqrt(2.0))
    for k in range(1, m + 1):
        p = -np.sqrt((2 * k + 1) / (2 * k)) * s * p
    out = np.zeros((L + 1,) + x.shape)
    out[m] = p
    if L > m:
        out[m + 1] = np.sqrt(2 * m + 3) * x * p
    a = lambda l: np.sqrt((4 * l * l - 1) / (l * l - m * m))
    for l in range(m + 2, L + 1):
        out[l] = a(l) * (x * out[l - 1] - out[l - 2] / a(l - 1))
    return out


def plm_unnorm(L, m, x):
    """textbook unnormalised recurrence (starts with (2m-1)!!: overflows for large m)"""
    x = np.asarray(x, float)
    with np.errstate(all="ignore"):
        dfact = np.prod(np.arange(1, 2 * m, 2, dtype=float)) if m > 0 else 1.0
        out = np.zeros((L + 1,) + x.shape)
        out[m] = (-1) ** m * dfact * (1 - x * x) ** (m / 2)
        if L > m: out[m + 1] = x * (2 * m + 1) * out[m]
        for l in range(m + 1, L):
            out[l + 1] = ((2 * l + 1) * x * out[l] - (l + m) * out[l - 1]) / (l - m + 1)
    return out


def gauss_legendre(n, tol=1e-15):
    """Newton on P_n with the Tricomi initial guesses cos(pi (i-1/4)/(n+1/2))"""
    x = np.zeros(n); w = np.zeros(n)
    for i in range(n):
        z = np.cos(np.pi * (i + 0.75) / (n + 0.5))
        for _ in range(100):
            P = legendre_rec(n, np.array([z]))[:, 0]
            dP = n * (z * P[n] - P[n - 1]) / (z * z - 1)
            dz = P[n] / dP; z -= dz
            if abs(dz) < tol: break
        P = legendre_rec(n, np.array([z]))[:, 0]
        dP = n * (z * P[n] - P[n - 1]) / (z * z - 1)
        x[i], w[i] = z, 2 / ((1 - z * z) * dP * dP)
    return x[::-1], w[::-1]


def partA():
    print("=" * 78)
    print("PART A  Legendre polynomials and associated Legendre functions")
    print("=" * 78)
    print("P_l(x): Bonnet recurrence vs explicit sum vs mpmath (50 digits)")
    print(f"{'l':>4} {'x':>6} {'recurrence':>24} {'mpmath':>24} {'rel.err rec':>12} {'rel.err sum':>12}")
    for l, x in ((10, 0.3), (40, 0.3), (60, 0.05), (100, 0.9)):
        ref = float(mp.legendre(l, x))
        r = legendre_rec(l, np.array([x]))[l, 0]; s = legendre_sum(l, np.array([x]))[0]
        print(f"{l:4d} {x:6.2f} {r:24.15e} {ref:24.15e} {abs(r-ref)/abs(ref):12.2e} {abs(s-ref)/abs(ref):12.2e}")
    print("\nNormalised Pbar_l^m vs unnormalised recurrence vs mpmath, x = 0.5")
    for l, m in ((5, 2), (30, 10), (120, 100), (200, 180)):
        nrm = mp.sqrt((2 * l + 1) / mp.mpf(2) * mp.factorial(l - m) / mp.factorial(l + m))
        ref = nrm * mp.legenp(l, m, 0.5)                     # mpmath includes (-1)^m
        pb = plm_bar(l, m, np.array([0.5]))[l, 0]
        pu = plm_unnorm(l, m, np.array([0.5]))[l, 0]
        pu_n = float(nrm) * pu if np.isfinite(pu) else np.nan
        print(f"  l={l:3d} m={m:3d}: Pbar = {pb: .15e}  mpmath = {float(ref): .15e}  rel.err {abs(pb-float(ref))/abs(float(ref)):.1e}"
              f"   unnormalised: {'overflow (inf)' if not np.isfinite(pu) else f'{pu_n: .6e}'}")
    xg, wg = gauss_legendre(30)
    P = plm_bar(20, 3, xg)[3:]
    G = (P * wg) @ P.T
    print(f"\northonormality of Pbar_l^3, l=3..20, 30-point Gauss-Legendre: max|G-I| = {np.max(abs(G-np.eye(len(G)))):.1e}")
    ls = np.arange(3, 21)
    scip = np.array([lpmv(3, l, xg) * np.sqrt((2 * l + 1) / 2 * exp(lgamma(l - 2) - lgamma(l + 4))) for l in ls])
    print(f"agreement with scipy lpmv (normalised): {np.max(abs(P - scip)):.1e}")
    print(f"30-point rule: sum w = {wg.sum():.15f}, int x^58 dx exact? err = {abs(wg@xg**58 - 2/59):.1e}")

    # error data
    xx = np.linspace(-0.999, 0.999, 400)
    ref40 = np.array([float(mp.legendre(40, x)) for x in xx])
    e_rec = abs(legendre_rec(40, xx)[40] - ref40) / (abs(ref40) + 1e-3)
    e_sum = abs(legendre_sum(40, xx) - ref40) / (abs(ref40) + 1e-3)
    ns = np.arange(2, 41, 2)
    tests = {r"$e^x$": (np.exp, np.e - 1 / np.e), r"$1/(1+25x^2)$": (lambda x: 1 / (1 + 25 * x * x), 0.4 * np.arctan(5)),
             r"$|x|^3$": (lambda x: abs(x) ** 3, 0.5)}
    eg = {k: [] for k in tests}; et = []
    for n in ns:
        xq, wq = gauss_legendre(n)
        for k, (f, ex) in tests.items():
            eg[k].append(abs(wq @ f(xq) - ex) + 1e-18)
        xt = np.linspace(-1, 1, n); ht = xt[1] - xt[0]
        ft = np.exp(xt); et.append(abs(ht * (ft.sum() - 0.5 * (ft[0] + ft[-1])) - (np.e - 1 / np.e)))
    print("\nGauss-Legendre quadrature errors (n nodes):")
    for i, n in enumerate(ns):
        if n in (2, 4, 8, 12, 20, 40):
            print(f"  n={n:2d}: " + "  ".join(f"{k.strip('$')}: {eg[k][i]:.1e}" for k in tests) + f"   trapezoid e^x: {et[i]:.1e}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    x = np.linspace(-1, 1, 500)
    P = legendre_rec(5, x)
    for l in range(6):
        ax[0].plot(x, P[l], color=C[l], lw=2, label=f"$P_{l}$")
    ax[0].set(xlabel="x", title=r"(a) Legendre polynomials $P_l(x)$")
    ax[0].legend(frameon=False, fontsize=8, ncol=3)
    for i, (l, m) in enumerate(((2, 0), (2, 1), (2, 2), (4, 2), (6, 3))):
        ax[1].plot(x, plm_bar(l, m, x)[l], color=C[i], lw=2, label=fr"$\bar P_{l}^{{{m}}}$")
    ax[1].set(xlabel="x", title=r"(b) Orthonormal associated functions $\bar P_l^m(x)$")
    ax[1].legend(frameon=False, fontsize=8, ncol=3)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch2_functions.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].semilogy(xx, e_sum + 1e-18, color=C[1], lw=1.2, label="explicit sum")
    ax[0].semilogy(xx, e_rec + 1e-18, color=C[0], lw=1.2, label="Bonnet recurrence")
    ax[0].set(xlabel="x", ylabel="error / (|P|+10$^{-3}$)", title=r"(a) $P_{40}(x)$ vs 50-digit reference")
    ax[0].legend(frameon=False, fontsize=8)
    for i, k in enumerate(tests):
        ax[1].semilogy(ns, eg[k], "o-", color=C[i], lw=2, ms=4, label=f"Gauss–Legendre, {k}")
    ax[1].semilogy(ns, et, "s--", color="#888888", lw=1.5, ms=4, label=r"trapezoid, $e^x$")
    ax[1].set(ylim=(1e-17, 1), xlabel="number of nodes n", ylabel="absolute error", title=r"(b) $\int_{-1}^{1} f\,dx$: quadrature convergence")
    ax[1].legend(frameon=False, fontsize=7)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch2_partA_errors.png", dpi=170)


# ============================================================ PART B
INV_SIN2 = [1.0, 1 / 3, 1 / 15, 2 / 189, 1 / 675, 2 / 10395, 1382 / 58046625, 4 / 1403325]


def frobenius(t, m, lam, nterm=24):
    """u = t^{m+1/2} sum a_k t^{2k};  1/sin^2 t = sum b_j t^{2j-2}
       2k(2k+2m) a_k = -(1/4-m^2) sum_{j>=1} b_j a_{k-j} - (lam+1/4) a_{k-1}"""
    q = 0.25 - m * m; L = lam + 0.25
    b = INV_SIN2 + [0.0] * nterm   # higher b_j negligible for t <= 0.5 at double precision
    a = [1.0]
    for k in range(1, nterm):
        s = -L * a[k - 1] - q * sum(b[j] * a[k - j] for j in range(1, k + 1))
        a.append(s / (2 * k * (2 * k + 2 * m)))
    t = np.asarray(t, float)
    return t ** (m + 0.5) * sum(ak * t ** (2 * k) for k, ak in enumerate(a))


def simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def shoot(lam, m, N, ts):
    """N even; grid t_j = j pi/N.  ts = start angle (0 -> first nodes)."""
    h = np.pi / N; t = h * np.arange(N + 1); t[N] = np.pi
    k2 = np.zeros(N + 1)
    k2[1:N] = lam + 0.25 + (0.25 - m * m) / np.sin(t[1:N]) ** 2
    c = 1 + h * h * k2 / 12
    js = max(1, int(ts / h + 0.5))
    M = N // 2
    pL = np.zeros(N + 1); pL[:js + 2] = frobenius(t[:js + 2], m, lam)
    for n in range(js + 1, M + 1):
        pL[n + 1] = ((12 - 10 * c[n]) * pL[n] - c[n - 1] * pL[n - 1]) / c[n + 1]
    pR = np.zeros(N + 1); pR[N - js - 1:] = frobenius(np.maximum(np.pi - t[N - js - 1:], 0.0), m, lam)
    for n in range(N - js - 1, M - 1, -1):
        pR[n - 1] = ((12 - 10 * c[n]) * pR[n] - c[n + 1] * pR[n + 1]) / c[n - 1]
    sL = np.max(abs(pL[:M + 2])); sR = np.max(abs(pR[M - 1:]))
    pL /= sL; pR /= sR
    D = pL[M + 1] * pR[M] - pL[M] * pR[M + 1]
    # glue at M+1 or M, whichever is larger in |pR| (odd states vanish at pi/2)
    j = M if abs(pR[M]) > abs(pR[M + 1]) else M + 1
    u = np.concatenate([pL[:j + 1], pR[j + 1:] * pL[j] / pR[j]])
    return D, u, t


def eigen(l, m, N, ts=0.4):
    lex = l * (l + 1)
    a, b = lex - l - 0.5, lex + l + 1.5            # brackets only lambda = l(l+1)
    fa = shoot(a, m, N, ts)[0]
    assert np.sign(fa) != np.sign(shoot(b, m, N, ts)[0])
    for _ in range(64):
        mid = 0.5 * (a + b); fm = shoot(mid, m, N, ts)[0]
        if np.sign(fm) == np.sign(fa): a, fa = mid, fm
        else: b = mid
    lam = 0.5 * (a + b)
    _, u, t = shoot(lam, m, N, ts)
    h = t[1]
    u /= np.sqrt(simpson(u * u, h))
    ex = np.sqrt(np.abs(np.sin(t))) * plm_bar(l, m, np.cos(t))[l]
    if np.dot(u, ex) < 0: u = -u
    return lam, u, t, ex


def partB():
    print("\n" + "=" * 78)
    print("PART B  Angular equation / rigid rotor: Numerov in theta, N = 400 (h = pi/400), t_s = 0.4")
    print("=" * 78)
    # scan to show the spectrum for m = 1
    N = 400
    lg = np.linspace(0.05, 45, 900)
    D = np.array([shoot(L, 1, N, 0.4)[0] for L in lg])
    roots = [0.5 * (lg[i] + lg[i + 1]) for i in np.where(np.sign(D[:-1]) != np.sign(D[1:]))[0]]
    print("  scan m=1: sign changes of D(lambda) near " + ", ".join(f"{r:.2f}" for r in roots) + "  (exact 2, 6, 12, 20, 30, 42)")
    print(f"\n{'l':>2} {'m':>2} {'lambda Numerov':>18} {'l(l+1)':>7} {'|dlam|':>9} {'norm':>15} {'max|u-u_ex|':>12} {'<cos^2>':>10} {'exact':>10}")
    store = {}
    for (l, m) in ((0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (3, 1), (2, 2), (3, 2), (5, 3)):
        lam, u, t, ex = eigen(l, m, N)
        h = t[1]
        c2 = simpson(u * u * np.cos(t) ** 2, h)
        c2ex = (2 * l * (l + 1) - 2 * m * m - 1) / ((2 * l - 1) * (2 * l + 3))
        store[l, m] = (u, t, ex)
        print(f"{l:2d} {m:2d} {lam:18.12f} {l*(l+1):7d} {abs(lam-l*(l+1)):9.1e} {simpson(u*u,h):15.12f} {np.max(abs(u-ex)):12.1e} {c2:10.7f} {c2ex:10.7f}")
    ov = max(abs(simpson(store[a, 1][0] * store[b, 1][0], np.pi / N)) for a in (1, 2, 3) for b in (1, 2, 3) if a < b)
    print(f"  orthogonality m=1: max |<u_l|u_l'>| = {ov:.1e};   parity u(pi-t) = (-1)^(l+m) u(t): "
          f"max err {max(np.max(abs(store[l,m][0][::-1] - (-1)**(l+m)*store[l,m][0])) for (l,m) in store):.1e}")

    Ns = np.array([20, 40, 80, 160, 320, 640, 1280])
    hs = np.pi / Ns
    err = {}
    for m, l in ((0, 2), (1, 2), (2, 3)):
        for ts, lab in ((0.0, "first nodes"), (0.4, "t_s = 0.4")):
            err[m, lab] = [abs(eigen(l, m, N, ts)[0] - l * (l + 1)) for N in Ns]
    print("\nEigenvalue error vs h: Frobenius start at first two grid nodes  vs  at fixed t_s = 0.4")
    print(f"{'N':>6} {'h':>8}" + "".join(f" {'m=%d nodes' % m:>12} {'m=%d t_s' % m:>11}" for m in (0, 1, 2)))
    for i, N in enumerate(Ns):
        print(f"{N:6d} {hs[i]:8.5f}" + "".join(f" {err[m,'first nodes'][i]:12.2e} {err[m,'t_s = 0.4'][i]:11.2e}" for m in (0, 1, 2)))
    for m in (0, 1, 2):
        for lab in ("first nodes", "t_s = 0.4"):
            e = np.array(err[m, lab]); ok = e > 1e-12
            print(f"  m={m} {lab:12s}: observed order {np.polyfit(np.log(hs[ok][1:6]), np.log(e[ok][1:6]), 1)[0]:.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
    for i, (l, m) in enumerate(((0, 0), (1, 0), (2, 0), (2, 1), (3, 2))):
        u, t, ex = store[l, m]
        ax[0].plot(t, u, color=C[i], lw=2, label=f"l={l}, m={m}")
        ax[0].plot(t[::16], ex[::16], "o", ms=2.5, color="#333333")
    ax[0].set(xlabel=r"$\theta$", ylabel=r"$u(\theta)=\sqrt{\sin\theta}\,\Theta(\theta)$",
              title=r"(a) Numerov (lines) vs $\sqrt{\sin\theta}\,\bar P_l^m(\cos\theta)$ (dots)")
    ax[0].legend(frameon=False, fontsize=7, ncol=2)
    th = np.linspace(0, 2 * np.pi, 400)
    for i, (l, m) in enumerate(((1, 0), (2, 0), (2, 1), (3, 2))):
        rr = abs(plm_bar(l, m, np.cos(th))[l]) ** 2
        ax[1].plot(rr * np.sin(th), rr * np.cos(th), color=C[i], lw=2, label=fr"$|Y_{{{l}{m}}}|^2$")
    ax[1].set_aspect("equal"); ax[1].set(xlabel="x-z plane", title=r"(b) Polar plots of $|\bar P_l^m(\cos\theta)|^2$")
    ax[1].legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch2_rotor.png", dpi=170)

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for i, m in enumerate((0, 1, 2)):
        ax.loglog(hs, err[m, "first nodes"], "s--", color=C[i], lw=1.5, ms=5, label=f"m={m}: start at first nodes")
        ax.loglog(hs, err[m, "t_s = 0.4"], "o-", color=C[i], lw=2, ms=5, label=fr"m={m}: Frobenius start at $\theta_s$=0.4")
    ax.loglog(hs, err[1, "t_s = 0.4"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax.loglog(hs, err[1, "first nodes"][0] * (hs / hs[0]) ** 2, "-.", color="#888888", label=r"$\propto h^2$")
    ax.set(xlabel=r"$h=\pi/N$", ylabel=r"$|\lambda-l(l+1)|$", title="Angular equation: the order of convergence is $\\min(4, 2m)$\nunless the start is moved away from the pole")
    ax.legend(frameon=False, fontsize=7); style(ax)
    fig.tight_layout(); fig.savefig("ch2_convergence.png", dpi=170)


if __name__ == "__main__":
    partA(); partB()
