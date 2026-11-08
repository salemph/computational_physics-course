#!/usr/bin/env python3
"""
Problem 3 -- QUANTUM HARMONIC OSCILLATOR with the Numerov method
================================================================
Units: hbar = m = omega = 1.   V(x) = x^2 / 2
psi''(x) = -k^2(x) psi(x),  k^2(x) = 2 [E - V(x)] = 2E - x^2
Exact:  E_n = n + 1/2,   psi_n(x) = (pi^-1/4 / sqrt(2^n n!)) H_n(x) exp(-x^2/2)

Eigenvalue method: two-sided shooting ("outward + inward" integration)
  * outward Numerov from x = -L (psi = 0, next point tiny) up to node m+1
  * inward  Numerov from x = +L (psi = 0, next point tiny) down to node m-1
  * m = grid node closest to the right classical turning point x_t = sqrt(2E)
    (integrating INTO a classically forbidden region is unstable: the growing
     exponential swamps the decaying one; we therefore always integrate from
     the forbidden region TOWARDS the allowed one)
  * matching function (discrete Wronskian / Casoratian, has no poles):
        D(E) = psiL[m+1]*psiR[m] - psiL[m]*psiR[m+1]
    D(E) = 0  <=>  the two pieces are proportional  <=>  E is an eigenvalue
  * eigenvalues bracketed by scanning D(E), refined by bisection,
    labelled by counting nodes of the assembled wavefunction.
Verification: normalisation (Simpson), overlap with exact psi_n, orthogonality,
              virial theorem <V> = E/2, Richardson extrapolation.
Error analysis: discretisation error vs h (order 4), truncation error vs box L,
                round-off floor for very small h.
Outputs: tables + p3_wavefunctions.png, p3_convergence.png
"""
import numpy as np
from math import factorial
from numpy.polynomial.hermite import hermval
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def V(x):
    return 0.5 * x * x


def simpson(y, h):
    n = len(y) - 1
    assert n % 2 == 0, "Simpson needs an even number of intervals"
    return h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


def psi_exact(n, x):
    coef = np.zeros(n + 1); coef[n] = 1.0
    return np.pi ** -0.25 / np.sqrt(2.0 ** n * factorial(n)) * hermval(x, coef) * np.exp(-x * x / 2)


def make_grid(L, N):
    x = np.linspace(-L, L, N + 1)
    return x, x[1] - x[0]


def shoot(E, x, h):
    """Return (D, psi_assembled, m) for trial energy E."""
    N = len(x) - 1
    k2 = 2.0 * (E - V(x))
    c = 1.0 + h * h * k2 / 12.0
    xt = np.sqrt(max(2.0 * E, 0.0))                     # right turning point
    m = int(np.clip(np.searchsorted(x, xt), 2, N - 2))
    # outward: 0 .. m+1
    pL = np.zeros(m + 2)
    pL[0], pL[1] = 0.0, 1e-10
    for n in range(1, m + 1):
        pL[n + 1] = ((12.0 - 10.0 * c[n]) * pL[n] - c[n - 1] * pL[n - 1]) / c[n + 1]
    # inward: N .. m-1
    pR = np.zeros(N + 1)
    pR[N], pR[N - 1] = 0.0, 1e-10
    for n in range(N - 1, m - 1, -1):
        pR[n - 1] = ((12.0 - 10.0 * c[n]) * pR[n] - c[n + 1] * pR[n + 1]) / c[n - 1]
    # scale to avoid overflow (positive factors do not change the sign of D)
    sL = np.max(np.abs(pL)); sR = np.max(np.abs(pR[m - 1:]))
    pL /= sL; pR /= sR
    D = pL[m + 1] * pR[m] - pL[m] * pR[m + 1]           # discrete Wronskian
    psi = np.empty(N + 1)
    psi[:m + 1] = pL[:m + 1]
    psi[m + 1:] = pR[m + 1:] * (pL[m] / pR[m])           # glue at node m
    return D, psi, m


def count_nodes(psi):
    p = psi[np.abs(psi) > 1e-6 * np.max(np.abs(psi))]   # ignore numerical tails
    return int(np.sum(p[:-1] * p[1:] < 0))


def eigenstates(nmax=5, L=7.0, N=1400, Emax=None, dE=0.05, tol=1e-14):
    x, h = make_grid(L, N)
    Emax = Emax if Emax else nmax + 1.0
    Egrid = np.arange(0.05, Emax, dE)
    Dg = np.array([shoot(E, x, h)[0] for E in Egrid])
    out = []
    for i in np.where(np.sign(Dg[:-1]) != np.sign(Dg[1:]))[0]:
        a, b, fa = Egrid[i], Egrid[i + 1], Dg[i]
        while b - a > tol * max(1.0, a):                # bisection on D(E)
            mid = 0.5 * (a + b)
            fm = shoot(mid, x, h)[0]
            if np.sign(fm) == np.sign(fa):
                a, fa = mid, fm
            else:
                b = mid
        E = 0.5 * (a + b)
        _, psi, m = shoot(E, x, h)
        psi /= np.sqrt(simpson(psi * psi, h))           # NORMALISATION
        n = count_nodes(psi)
        ex = psi_exact(n, x)
        if np.dot(psi, ex) < 0:
            psi = -psi
        out.append((n, E, psi))
    return x, h, out


def slope(h, e):
    return np.polyfit(np.log(h), np.log(e), 1)[0]


def main():
    print("=" * 78)
    print("Harmonic oscillator  V = x^2/2, box [-7, 7], N = 1400 (h = 0.01)")
    print("=" * 78)
    x, h, st = eigenstates(nmax=5)
    print(f"{'n':>2} {'E_numerov':>18} {'E_exact':>6} {'|dE|':>10} {'norm':>14} "
          f"{'max|psi-ex|':>12} {'<V>':>10} {'E/2':>8}")
    for n, E, psi in st:
        ex = psi_exact(n, x)
        Vexp = simpson(psi * psi * V(x), h)
        print(f"{n:2d} {E:18.12f} {n+0.5:6.1f} {abs(E-n-0.5):10.2e} {simpson(psi*psi,h):14.12f} "
              f"{np.max(abs(psi-ex)):12.2e} {Vexp:10.6f} {E/2:8.4f}")
    P = np.array([s[2] for s in st])
    O = np.array([[simpson(P[i] * P[j], h) for j in range(len(P))] for i in range(len(P))])
    print("max |<psi_i|psi_j> - delta_ij| =", f"{np.max(np.abs(O - np.eye(len(P)))):.2e}")

    # -------------------------------------------- convergence in h (fixed L = 7)
    Ns = [70, 140, 280, 560, 1120, 2240, 4480, 8960]
    hs, errs = [], {0: [], 2: [], 4: []}
    for N in Ns:
        xg, hg, s = eigenstates(nmax=5, L=7.0, N=N)
        hs.append(hg)
        d = {n: E for n, E, _ in s}
        for n in errs:
            errs[n].append(abs(d[n] - (n + 0.5)))
    print("\nEigenvalue error vs step (L = 7):")
    print(f"{'N':>6} {'h':>9}" + "".join(f" {'|dE_'+str(n)+'|':>11} {'ratio':>6}" for n in errs))
    for i, N in enumerate(Ns):
        line = f"{N:6d} {hs[i]:9.5f}"
        for n in errs:
            r = errs[n][i - 1] / errs[n][i] if i else np.nan
            line += f" {errs[n][i]:11.2e} {r:6.1f}"
        print(line)
    for n in errs:
        print(f"  n={n}: observed order p = {slope(hs[:5], errs[n][:5]):.2f}")

    # Richardson extrapolation for E_4
    _, _, s1 = eigenstates(L=7.0, N=140); _, _, s2 = eigenstates(L=7.0, N=280)
    E1 = dict((n, E) for n, E, _ in s1)[4]; E2 = dict((n, E) for n, E, _ in s2)[4]
    ER = (16 * E2 - E1) / 15
    print(f"\nRichardson n=4: E(h=0.1)={E1:.10f}  E(h=0.05)={E2:.10f}  E_R={ER:.10f}  err={abs(ER-4.5):.1e}")

    # -------------------------------------------- truncation error vs L (h = 0.005)
    Ls = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
    errL = {0: [], 2: [], 4: []}
    for L in Ls:
        N = int(round(2 * L / 0.005)); N += N % 2
        _, _, s = eigenstates(nmax=5, L=L, N=N, Emax=6.0)
        d = {n: E for n, E, _ in s}
        for n in errL:
            errL[n].append(abs(d.get(n, np.nan) - (n + 0.5)))
    print("\nBox-truncation error (h = 0.005):")
    print(f"{'L':>5}" + "".join(f" {'|dE_'+str(n)+'|':>11}" for n in errL))
    for i, L in enumerate(Ls):
        print(f"{L:5.1f}" + "".join(f" {errL[n][i]:11.2e}" for n in errL))

    # -------------------------------------------- figures
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    xx = np.linspace(-5, 5, 400)
    ax[0].plot(xx, V(xx), color="#888888", lw=1)
    for i, (n, E, psi) in enumerate(st[:5]):
        sel = (x > -5) & (x < 5)
        ax[0].axhline(E, color="#dddddd", lw=.8)
        ax[0].plot(x[sel], E + 0.6 * psi[sel], color=C[i], lw=2, label=f"n={n}, E={E:.6f}")
        ax[0].plot(x[sel][::40], E + 0.6 * psi_exact(n, x[sel][::40]), "o", ms=3, color="#333333")
    ax[0].set(xlim=(-5, 5), ylim=(0, 7.2), xlabel="x", ylabel=r"$E_n + 0.6\,\psi_n(x)$",
              title="(a) Normalised eigenfunctions (dots: exact)")
    ax[0].legend(frameon=False, fontsize=7, loc="upper center", ncol=2)
    for i, (n, E, psi) in enumerate(st[:5]):
        ax[1].semilogy(x, np.abs(psi - psi_exact(n, x)) + 1e-18, color=C[i], lw=1.2, label=f"n={n}")
    ax[1].set(xlim=(-7, 7), ylim=(1e-14, 1e-5), xlabel="x", ylabel=r"$|\psi_{num}-\psi_{exact}|$",
              title="(b) Pointwise error, h = 0.01")
    ax[1].legend(frameon=False, fontsize=7, ncol=5, loc="upper center")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25)
    fig.tight_layout(); fig.savefig("p3_wavefunctions.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    hs_a = np.array(hs)
    for i, n in enumerate(errs):
        ax[0].loglog(hs_a, errs[n], "o-", color=C[i], lw=2, ms=5, label=f"$E_{n}$")
    ax[0].loglog(hs_a, errs[4][0] * (hs_a / hs_a[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax[0].set(xlabel="step h", ylabel=r"$|E_n^{num}-(n+1/2)|$",
              title="(a) Discretisation error (L = 7)")
    ax[0].legend(frameon=False, fontsize=8)
    for i, n in enumerate(errL):
        ax[1].semilogy(Ls, errL[n], "o-", color=C[i], lw=2, ms=5, label=f"$E_{n}$")
    ax[1].set(xlabel="half-width of box L", ylabel=r"$|E_n^{num}-(n+1/2)|$",
              title="(b) Truncation error of the box (h = 0.005)")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25, which="both")
    fig.tight_layout(); fig.savefig("p3_convergence.png", dpi=170)
    print("\nfigures written: p3_wavefunctions.png, p3_convergence.png")


if __name__ == "__main__":
    main()
