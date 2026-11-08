#!/usr/bin/env python3
"""
Chapter 4 -- SPHERICAL BESSEL FUNCTIONS and the SPHERICAL INFINITE WELL
=======================================================================
Part A : computation of j_l(x), y_l(x)
         * closed forms j_0 = sin x / x, j_1, y_0, y_1           DLMF §10.49
         * upward recurrence f_{l+1} = (2l+1)/x f_l - f_{l-1}   DLMF §10.51
           (stable for y_l, UNSTABLE for j_l when l > x)
         * Miller downward recurrence for j_l, normalised by j_0 = sin x / x
         * small-x series  j_l ~ x^l/(2l+1)!! (1 - x^2/(2(2l+3)) + ...)  DLMF §10.53
         * zeros z_{l,s}: bisection bracket + Newton, j_l' = j_{l-1} - (l+1) j_l / x
         * verification: cross product j_l y_{l-1} - j_{l-1} y_l = x^{-2} (DLMF §10.50),
           mpmath reference, scipy spherical_jn
Part B : particle in a spherical infinite well of radius 1 (hbar = m = 1)
         psi = R(r) Y_lm,  u = r R :
         u'' + [ 2E - l(l+1)/r^2 ] u = 0 ,  u(1) = 0 ,  u ~ r^{l+1}  (analytic!)
         exact  E_{l,s} = z_{l,s}^2/2,  R = sqrt(2) j_l(z r)/|j_{l+1}(z)|
         Numerov outward shooting, bisection on u(1)
         * verification: normalisation, comparison with j_l, orthogonality, <r^2>
         * error analysis: here the first-node start is harmless (u is analytic),
           in contrast with Chapters 2-3; h^4 law, Richardson
Outputs: tables + ch4_functions.png, ch4_partA_errors.png, ch4_well.png, ch4_convergence.png
"""
import numpy as np
import mpmath as mp
from scipy.special import spherical_jn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cd6"]
mp.mp.dps = 50


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


def ref_j(l, x):
    return float(mp.sqrt(mp.pi / (2 * mp.mpf(x))) * mp.besselj(l + mp.mpf(1) / 2, x))


def ref_y(l, x):
    return float(mp.sqrt(mp.pi / (2 * mp.mpf(x))) * mp.bessely(l + mp.mpf(1) / 2, x))


# ============================================================ PART A
def j_up(L, x):
    j = np.zeros(L + 1); j[0] = np.sin(x) / x
    if L >= 1: j[1] = np.sin(x) / x ** 2 - np.cos(x) / x
    for l in range(1, L):
        j[l + 1] = (2 * l + 1) / x * j[l] - j[l - 1]
    return j


def y_up(L, x):
    y = np.zeros(L + 1); y[0] = -np.cos(x) / x
    if L >= 1: y[1] = -np.cos(x) / x ** 2 - np.sin(x) / x
    for l in range(1, L):
        y[l + 1] = (2 * l + 1) / x * y[l] - y[l - 1]
    return y


def j_series(l, x, nterm=30):
    df = np.prod(np.arange(1, 2 * l + 2, 2, dtype=float))
    s, t = 0.0, 1.0
    for k in range(nterm):
        s += t
        t *= -x * x / (2 * (k + 1) * (2 * l + 2 * k + 3))
    return x ** l / df * s


def j_miller(L, x, M=None):
    """downward recurrence from M, normalised with j_0 = sin x / x"""
    if x < 1e-3:                                    # tiny x: series is exact to round-off
        return np.array([j_series(l, x) for l in range(L + 1)])
    if M is None:
        M = L + int(x) + 20 + int(np.sqrt(40 * (L + x)))
    f = np.zeros(M + 2); f[M] = 1e-30
    for l in range(M, 0, -1):
        f[l - 1] = (2 * l + 1) / x * f[l] - f[l + 1]
        if abs(f[l - 1]) > 1e250: f[l - 1:] *= 1e-250
    j0 = np.sin(x) / x
    # normalise with the larger of j_0, j_1 (j_0 vanishes at x = k pi)
    j1 = np.sin(x) / x ** 2 - np.cos(x) / x
    scale = j0 / f[0] if abs(j0) > abs(j1) else j1 / f[1]
    return f[:L + 1] * scale


def jl(l, x):
    return j_miller(l + 1, x)[l]


def jl_zeros(l, S):
    """first S zeros: bracket by scanning (step 0.1), then Newton"""
    zs = []
    x = l + 0.5 if l > 0 else 1.0
    prev = jl(l, x)
    while len(zs) < S:
        x2 = x + 0.1; cur = jl(l, x2)
        if np.sign(cur) != np.sign(prev):
            z = 0.5 * (x + x2)
            for _ in range(50):
                jj = j_miller(l + 1, z)
                d = (jj[l - 1] - (l + 1) / z * jj[l]) if l > 0 else -jj[1]
                dz = jj[l] / d; z -= dz
                if abs(dz) < 1e-15 * z: break
            zs.append(z)
        x, prev = x2, cur
    return np.array(zs)


def partA():
    print("=" * 78)
    print("PART A  Spherical Bessel functions j_l, y_l")
    print("=" * 78)
    print(f"{'l':>3} {'x':>6} {'j_l Miller':>24} {'mpmath':>24} {'rel.err Miller':>15} {'rel.err upward':>15}")
    for l, x in ((2, 5.0), (5, 1.0), (10, 2.0), (20, 30.0), (30, 1.0), (3, 1e-4)):
        ref = ref_j(l, x); jm = jl(l, x); ju = j_up(l, x)[l]
        print(f"{l:3d} {x:6.0e} {jm:24.15e} {ref:24.15e} {abs(jm-ref)/abs(ref):15.2e} {abs(ju-ref)/abs(ref):15.2e}")
    print("\ny_l by upward recurrence (stable):")
    for l, x in ((5, 1.0), (20, 3.0), (30, 1.0)):
        ref = ref_y(l, x); yu = y_up(l, x)[l]
        print(f"  l={l:2d} x={x}: y_l = {yu:.15e}  mpmath {ref:.15e}  rel.err {abs(yu-ref)/abs(ref):.1e}")
    print("\nCross product  j_l y_{l-1} - j_{l-1} y_l - 1/x^2  (independent check, DLMF §10.50):")
    for x in (0.5, 3.0, 20.0):
        jj = j_miller(25, x); yy = y_up(25, x)
        cp = max(abs((jj[l] * yy[l - 1] - jj[l - 1] * yy[l]) * x * x - 1) for l in range(1, 26))
        print(f"  x={x:5.1f}: max over l=1..25 of |x^2 (j_l y_(l-1) - j_(l-1) y_l) - 1| = {cp:.1e}")
    xs = np.linspace(0.05, 40, 80)
    print(f"  max |j_miller - scipy spherical_jn| (l<=15, 0.05<x<40) = "
          f"{max(np.max(abs(j_miller(15, x) - spherical_jn(np.arange(16), x))) for x in xs):.1e}")
    z = {l: jl_zeros(l, 4) for l in range(4)}
    for l in range(4):
        refz = [float(mp.findroot(lambda t: mp.sqrt(mp.pi/(2*t))*mp.besselj(l+0.5, t), v)) for v in z[l]]
        print(f"  zeros j_{l}: " + ", ".join(f"{v:.12f}" for v in z[l]) + f"   max err vs mpmath {max(abs(z[l]-refz)):.1e}")

    # error data
    ls = np.arange(0, 31)
    ref1 = np.array([ref_j(l, 1.0) for l in ls])
    e_up = abs(j_up(30, 1.0) - ref1) / abs(ref1)
    e_dn = abs(j_miller(30, 1.0) - ref1) / abs(ref1)
    xx = np.logspace(-6, 0, 120)
    ref3 = np.array([ref_j(3, x) for x in xx])
    j3_closed = (15 / xx ** 4 - 6 / xx ** 2) * np.sin(xx) - (15 / xx ** 3 - 1 / xx) * np.cos(xx)  # DLMF §10.49
    e_cl = abs(j3_closed / xx * 1 - ref3) / abs(ref3)
    e_cl = abs(((15 / xx ** 3 - 6 / xx) * np.sin(xx) - (15 / xx ** 2 - 1) * np.cos(xx)) / xx - ref3) / abs(ref3)
    e_mi = np.array([abs(jl(3, x) - r) / abs(r) for x, r in zip(xx, ref3)])
    print("\nj_3 closed trigonometric form vs Miller at small x (relative error):")
    for x in (1e-1, 1e-2, 1e-3, 1e-4):
        i = np.argmin(abs(xx - x))
        print(f"  x={xx[i]:.0e}: closed form {e_cl[i]:.1e}   Miller/series {e_mi[i]:.1e}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    xg = np.linspace(0.01, 15, 500)
    for l in range(4):
        ax[0].plot(xg, [jl(l, v) for v in xg], color=C[l], lw=2, label=f"$j_{l}$")
        ax[0].plot(z[l][:3], 0 * z[l][:3], "o", ms=4, color=C[l])
    ax[0].axhline(0, color="#888888", lw=.8)
    ax[0].set(xlabel="x", title=r"(a) $j_l(x)$ and zeros $z_{l,s}$")
    ax[0].legend(frameon=False, fontsize=8, ncol=4)
    for l in range(3):
        ax[1].plot(xg[20:], [y_up(l, v)[l] for v in xg[20:]], color=C[l], lw=2, label=f"$y_{l}$")
    ax[1].set(ylim=(-1.5, 0.5), xlabel="x", title=r"(b) $y_l(x)$ (singular at the origin)")
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch4_functions.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].semilogy(ls, e_up + 1e-18, "s-", color=C[1], lw=1.5, ms=4, label="upward recurrence")
    ax[0].semilogy(ls, e_dn + 1e-18, "o-", color=C[0], lw=1.5, ms=4, label="downward (Miller)")
    ax[0].set(xlabel="l", ylabel="relative error", title=r"(a) $j_l(1)$: direction of recurrence")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].loglog(xx, e_cl + 1e-18, color=C[1], lw=1.5, label="closed trigonometric form")
    ax[1].loglog(xx, e_mi + 1e-18, color=C[0], lw=1.5, label="Miller / series")
    ax[1].set(xlabel="x", ylabel="relative error", title=r"(b) $j_3(x)$ for small x: cancellation")
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch4_partA_errors.png", dpi=170)
    return z


# ============================================================ PART B
def frobenius(r, l, E, nterm=40):
    """u = r^{l+1} sum a_k r^{2k},  2k(2k+2l+1) a_k = -2E a_{k-1}"""
    r = np.asarray(r, float); s = np.zeros_like(r); a = 1.0; p = np.ones_like(r)
    for k in range(nterm):
        s += a * p
        a = -2 * E * a / ((2 * k + 2) * (2 * k + 2 * l + 3)); p = p * r * r
    return r ** (l + 1) * s


def simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def shoot(E, l, N, rs, naive=False):
    h = 1.0 / N; r = h * np.arange(N + 1)
    k2 = np.zeros(N + 1); k2[1:] = 2 * E - l * (l + 1) / r[1:] ** 2
    c = 1 + h * h * k2 / 12
    if naive:                                   # u0 = 0, u1 = h^{l+1}, c0 u0 = 0
        js = 0; u = np.zeros(N + 1); u[1] = h ** (l + 1)
    else:
        js = max(1, int(rs / h + 0.5))
        u = np.zeros(N + 1); u[:js + 2] = frobenius(r[:js + 2], l, E)
    for n in range(js + 1, N):
        u[n + 1] = ((12 - 10 * c[n]) * u[n] - c[n - 1] * u[n - 1]) / c[n + 1]
    return u, r


def eigen(l, s, N, z, rs=0.0, naive=False):
    Ex = z[l][s - 1] ** 2 / 2
    a, b = Ex * 0.85, Ex * 1.15
    fa = shoot(a, l, N, rs, naive)[0][-1]
    for _ in range(64):
        mid = 0.5 * (a + b); fm = shoot(mid, l, N, rs, naive)[0][-1]
        if np.sign(fm) == np.sign(fa): a, fa = mid, fm
        else: b = mid
    E = 0.5 * (a + b)
    u, r = shoot(E, l, N, rs, naive)
    u /= np.sqrt(simpson(u * u, r[1]))
    zz = z[l][s - 1]
    ex = r * np.array([jl(l, zz * v) if v > 0 else (1.0 if l == 0 else 0.0) * zz for v in r]) * np.sqrt(2) / abs(jl(l + 1, zz))
    ex[0] = 0.0
    if np.dot(u, ex) < 0: u = -u
    return E, u, r, ex


def r2_exact(l, z):
    """<r^2> = int_0^1 j_l(zr)^2 r^4 dr / int_0^1 j_l(zr)^2 r^2 dr, by 50-digit quadrature"""
    f = lambda t: mp.besselj(l + mp.mpf(1) / 2, z * t) ** 2 / t       # = (2 z/pi) j_l^2 * t^... common factor
    num = mp.quad(lambda t: f(t) * t ** 4, [0, 1]); den = mp.quad(lambda t: f(t) * t ** 2, [0, 1])
    return float(num / den)


def order(hs, e):
    """observed order, fitted only on the range where the error still drops (>6x per halving)"""
    idx = [0]
    for i in range(len(e) - 1):
        if e[i] / e[i + 1] > 6: idx.append(i + 1)
        else: break
    idx = idx[1:] if len(idx) > 3 else idx
    return np.polyfit(np.log(np.array(hs)[idx]), np.log(np.array(e)[idx]), 1)[0]


def partB(z):
    print("\n" + "=" * 78)
    print("PART B  Spherical infinite well, radius 1: Numerov in u = r R, N = 400, Frobenius start at first nodes")
    print("=" * 78)
    Eg = np.linspace(0.5, 80, 800)
    F = np.array([shoot(E, 1, 400, 0.0)[0][-1] for E in Eg])
    roots = [0.5 * (Eg[i] + Eg[i + 1]) for i in np.where(np.sign(F[:-1]) != np.sign(F[1:]))[0]]
    print("  scan l=1: sign changes of u(1;E) near " + ", ".join(f"{v:.2f}" for v in roots)
          + "  (exact " + ", ".join(f"{v**2/2:.2f}" for v in z[1][:4]) + ")")
    print(f"\n{'l':>2} {'s':>2} {'E Numerov':>17} {'z^2/2 (Part A)':>17} {'|dE|':>9} {'norm':>15} {'max|u-u_ex|':>12} {'<r^2>':>9} {'exact':>9}")
    store = {}
    for (l, s) in ((0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (3, 1)):
        E, u, r, ex = eigen(l, s, 400, z)
        Ex = z[l][s - 1] ** 2 / 2
        zz = z[l][s - 1]
        store[l, s] = (u, r, ex)
        print(f"{l:2d} {s:2d} {E:17.12f} {Ex:17.12f} {abs(E-Ex):9.1e} {simpson(u*u,r[1]):15.12f} {np.max(abs(u-ex)):12.1e} "
              f"{simpson(u*u*r*r,r[1]):9.6f} {r2_exact(l, zz):9.6f}")
    ov = max(abs(simpson(store[0, 1][0] * store[0, 2][0], 1 / 400)), abs(simpson(store[1, 1][0] * store[1, 2][0], 1 / 400)))
    print(f"  orthogonality (same l, different s): max |<u|u'>| = {ov:.1e}")

    Ns = np.array([20, 40, 80, 160, 320, 640, 1280])
    hs = 1.0 / Ns
    err = {}
    for l in (0, 1, 2):
        Ex = z[l][0] ** 2 / 2
        err[l, "naive"] = [abs(eigen(l, 1, N, z, naive=True)[0] - Ex) for N in Ns]
        err[l, "Frobenius"] = [abs(eigen(l, 1, N, z)[0] - Ex) for N in Ns]
    print("\nError of E_{l,1} vs h: naive start (u0=0, u1=h^(l+1))  vs  Frobenius values at the first nodes")
    print(f"{'N':>6} {'h':>9}" + "".join(f" {'l=%d naive' % l:>12} {'l=%d Frob.' % l:>11}" for l in (0, 1, 2)))
    for i, N in enumerate(Ns):
        print(f"{N:6d} {hs[i]:9.6f}" + "".join(f" {err[l,'naive'][i]:12.2e} {err[l,'Frobenius'][i]:11.2e}" for l in (0, 1, 2)))
    for l in (0, 1, 2):
        for lab in ("naive", "Frobenius"):
            print(f"  l={l} {lab:10s}: observed order {order(hs, err[l, lab]):.2f}")
    E1 = eigen(1, 1, 80, z)[0]; E2 = eigen(1, 1, 160, z)[0]; ER = (16 * E2 - E1) / 15
    print(f"Richardson l=1: E(1/80)={E1:.12f} E(1/160)={E2:.12f} E_R={ER:.12f} err={abs(ER-z[1][0]**2/2):.1e}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    for i, (l, s) in enumerate(((0, 1), (0, 2), (1, 1), (2, 1), (3, 1))):
        u, r, ex = store[l, s]
        ax[0].plot(r, u, color=C[i], lw=2, label=f"l={l}, s={s}")
        ax[0].plot(r[::16], ex[::16], "o", ms=2.5, color="#333333")
    ax[0].set(xlabel="r", ylabel="u(r) = r R(r)", title=r"(a) Numerov vs $\sqrt{2}\,r\,j_l(zr)/|j_{l+1}(z)|$ (dots)")
    ax[0].legend(frameon=False, fontsize=7, ncol=2)
    lev = sorted([(z[l][s] ** 2 / 2, l, s + 1) for l in range(4) for s in range(3) if z[l][s] ** 2 / 2 < 60])
    for E, l, s in lev:
        ax[1].hlines(E, l - 0.3, l + 0.3, color=C[l], lw=2.5)
        ax[1].text(l + 0.33, E, f"{s}{'spdf'[l]}", fontsize=7, va="center")
    ax[1].set(xticks=range(4), xticklabels=["s (l=0)", "p (l=1)", "d (l=2)", "f (l=3)"], ylabel="E",
              title="(b) Level scheme $E=z_{l,s}^2/2$ (cf. nuclear shell model)")
    for a_ in ax: style(a_)
    fig.tight_layout(); fig.savefig("ch4_well.png", dpi=170)

    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    for i, l in enumerate((0, 1, 2)):
        ax.loglog(hs, err[l, "naive"], "s--", color=C[i], lw=1.5, ms=5, label=f"l={l}: naive start $u_1=h^{{l+1}}$")
        ax.loglog(hs, err[l, "Frobenius"], "o-", color=C[i], lw=2, ms=5, label=f"l={l}: Frobenius values at first nodes")
    ax.loglog(hs, err[0, "Frobenius"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax.loglog(hs, err[1, "naive"][0] * (hs / hs[0]) ** 3, "-.", color="#888888", label=r"$\propto h^3$")
    ax.set(xlabel="h = 1/N", ylabel=r"$|E_{l,1}-z_{l,1}^2/2|$", title="Spherical well: analytic u = r R needs no special start\n(only the naive $u_1$ for l = 1 is wrong)")
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2); style(ax)
    fig.tight_layout(); fig.savefig("ch4_convergence.png", dpi=170)


if __name__ == "__main__":
    z = partA(); partB(z)
