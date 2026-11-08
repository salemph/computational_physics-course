#!/usr/bin/env python3
"""
Problem 1 -- FREE PARTICLE solved with the Numerov method
==========================================================
Units: hbar = m = 1  ->  psi''(x) = -k^2(x) psi(x),  k^2(x) = 2 [E - V(x)],  V = 0

Part A  (initial-value problem, unbounded free particle)
        psi(0) = 0, energy E = k^2/2 given  ->  exact solution psi(x) = sin(kx).
        We integrate with Numerov and measure the global error vs. step h.
        We also show the "starting-value fault": using a Taylor start
        psi(h) = k h (instead of the exact sin(kh)) degrades the global order.

Part B  (free particle confined to a box [0, L] -- hard walls)
        A free particle is NOT normalisable on the whole real line (plane waves),
        so to perform a normalisation check we confine it:
        psi(0) = psi(L) = 0,  E_n = n^2 pi^2 / (2 L^2),  psi_n = sqrt(2/L) sin(n pi x / L).
        Eigenvalues are found by shooting + bisection, eigenfunctions are
        normalised with Simpson's rule and compared with the analytic ones.

Outputs: printed tables + figures  p1_wavefunctions.png, p1_convergence.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]   # fixed series colours


# ---------------------------------------------------------------- Numerov core
def numerov(k2, psi0, psi1, h):
    """Integrate psi'' = -k2(x) psi on a uniform grid.
    k2   : array of k^2 at every grid point
    psi0, psi1 : two starting values
    Returns psi on the whole grid.
      psi_{n+1} (1 + h^2 k2_{n+1}/12) = 2 psi_n (1 - 5 h^2 k2_n/12)
                                       - psi_{n-1} (1 + h^2 k2_{n-1}/12)
    """
    N = len(k2)
    c = 1.0 + h * h * k2 / 12.0                 # Numerov weights
    psi = np.zeros(N, dtype=np.result_type(psi0, psi1, float))
    psi[0], psi[1] = psi0, psi1
    for n in range(1, N - 1):
        psi[n + 1] = ((12.0 - 10.0 * c[n]) * psi[n] - c[n - 1] * psi[n - 1]) / c[n + 1]
    return psi


def simpson(y, h):
    """Composite Simpson rule (needs an even number of intervals)."""
    n = len(y) - 1
    if n % 2:
        raise ValueError("Simpson needs an even number of intervals")
    return h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


# ---------------------------------------------------------------- Part A
def part_A(k=2.0, X=10.0, N=200, exact_start=True):
    h = X / N
    x = np.linspace(0.0, X, N + 1)
    k2 = np.full_like(x, k * k)                        # V = 0
    psi1 = np.sin(k * h) if exact_start else k * h     # starting value
    psi = numerov(k2, 0.0, psi1, h)
    err = np.max(np.abs(psi - np.sin(k * x)))
    return x, psi, err


# ---------------------------------------------------------------- Part B
def shoot_box(E, L, N):
    h = L / N
    k2 = np.full(N + 1, 2.0 * E)
    psi = numerov(k2, 0.0, 1e-3 * h, h)                # psi(0)=0, small slope
    return psi


def box_eigen(n, L=1.0, N=200, tol=1e-13):
    """n-th eigenvalue (n = 1,2,..): bracket by scanning psi(L), then bisection."""
    Eexact = (n * np.pi / L) ** 2 / 2.0
    # bracket: scan an interval around the estimate (generic scan works too)
    Ea, Eb = 0.5 * Eexact * (1 - 1 / n) + 1e-3, Eexact * (1 + 1.0 / n)
    grid = np.linspace(Ea, Eb, 200)
    f = np.array([shoot_box(E, L, N)[-1] for E in grid])
    i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    # pick the root closest to the estimate (the n-th one)
    j = i[np.argmin(np.abs(grid[i] - Eexact))]
    a, b, fa = grid[j], grid[j + 1], f[j]
    while b - a > tol * max(1.0, abs(a)):
        m = 0.5 * (a + b)
        fm = shoot_box(m, L, N)[-1]
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    E = 0.5 * (a + b)
    psi = shoot_box(E, L, N)
    h = L / N
    psi /= np.sqrt(simpson(psi ** 2, h))               # NORMALISATION
    x = np.linspace(0, L, N + 1)
    exact = np.sqrt(2 / L) * np.sin(n * np.pi * x / L)
    if np.dot(psi, exact) < 0:
        psi = -psi                                     # fix arbitrary sign
    return E, Eexact, x, psi, exact


def slope(h, e):
    p = np.polyfit(np.log(h), np.log(e), 1)
    return p[0]


def main():
    print("=" * 70)
    print("PART A : free particle, psi(0)=0, k = 2, exact solution sin(kx)")
    print("=" * 70)
    Ns = [25, 50, 100, 200, 400, 800, 1600, 3200]
    hA, eA, eT = [], [], []
    print(f"{'N':>6} {'h':>10} {'err(exact start)':>18} {'ratio':>7} {'err(Taylor start)':>18} {'ratio':>7}")
    for i, N in enumerate(Ns):
        _, _, e1 = part_A(N=N, exact_start=True)
        _, _, e2 = part_A(N=N, exact_start=False)
        r1 = eA[-1] / e1 if eA else np.nan
        r2 = eT[-1] / e2 if eT else np.nan
        hA.append(10.0 / N); eA.append(e1); eT.append(e2)
        print(f"{N:6d} {10.0/N:10.5f} {e1:18.3e} {r1:7.2f} {e2:18.3e} {r2:7.2f}")
    print(f"observed order (exact start)  p = {slope(hA, eA):.2f}   (theory 4)")
    print(f"observed order (Taylor start) p = {slope(hA, eT):.2f}   (theory 2: start error O(h^3) / h)")

    print("\n" + "=" * 70)
    print("PART B : free particle in a box L = 1  (E_n = n^2 pi^2 / 2)")
    print("=" * 70)
    print(f"{'n':>3} {'E_numerov':>18} {'E_exact':>18} {'rel.err':>11} {'norm':>14} {'max|psi-exact|':>15}")
    states = {}
    for n in (1, 2, 3, 4):
        E, Ex, x, psi, ex = box_eigen(n, N=200)
        states[n] = (x, psi, ex)
        print(f"{n:3d} {E:18.12f} {Ex:18.12f} {abs(E-Ex)/Ex:11.3e} "
              f"{simpson(psi**2, x[1]-x[0]):14.12f} {np.max(abs(psi-ex)):15.3e}")
    # orthogonality check
    O = np.array([[simpson(states[i][1] * states[j][1], states[1][0][1]) for j in (1, 2, 3, 4)]
                  for i in (1, 2, 3, 4)])
    print("overlap matrix <psi_i|psi_j>  (should be identity):")
    print(np.array2string(O, precision=2, suppress_small=False, formatter={'float': lambda v: f"{v: .1e}"}))

    print("\nConvergence of the eigenvalues in the box:")
    NsB = [20, 40, 80, 160, 320, 640, 1280]
    hB = [1.0 / N for N in NsB]
    errB = {}
    for n in (1, 2, 3):
        errB[n] = []
        for N in NsB:
            E, Ex, *_ = box_eigen(n, N=N)
            errB[n].append(abs(E - Ex) / Ex)
        print(f"  n={n}: rel. errors", " ".join(f"{e:.2e}" for e in errB[n]),
              f"-> order {slope(hB[:5], errB[n][:5]):.2f}")

    # Richardson extrapolation of E_1 (error ~ h^4)
    E_h = box_eigen(1, N=40)[0]
    E_h2 = box_eigen(1, N=80)[0]
    E_R = (16 * E_h2 - E_h) / 15
    print(f"\nRichardson (n=1): E(h)={E_h:.12f} E(h/2)={E_h2:.12f} "
          f"E_R={E_R:.12f}  exact={np.pi**2/2:.12f}  err_R={abs(E_R-np.pi**2/2):.2e}")

    # ------------------------------------------------------------ figures
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    x, psi, _ = part_A(N=100)
    xe = np.linspace(0, 10, 1000)
    ax[0].plot(xe, np.sin(2 * xe), color="#888888", lw=1, label="exact  sin(kx)")
    ax[0].plot(x[::3], psi[::3], "o", ms=4, color=C[0], label="Numerov (N=100)")
    ax[0].set(xlabel="x", ylabel=r"$\psi(x)$", title="(a) Free particle, k = 2")
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    for i, n in enumerate((1, 2, 3, 4)):
        x, psi, ex = states[n]
        ax[1].plot(x, ex, color="#888888", lw=1)
        ax[1].plot(x[::8], psi[::8], "o", ms=3.5, color=C[i], label=f"n={n}")
    ax[1].set(xlabel="x", ylabel=r"$\psi_n(x)$", title="(b) Free particle in a box (normalised)")
    ax[1].legend(frameon=False, fontsize=8, ncol=4, loc="lower center")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25)
    fig.tight_layout(); fig.savefig("p1_wavefunctions.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].loglog(hA, eA, "o-", color=C[0], lw=2, label="exact start  sin(kh)")
    ax[0].loglog(hA, eT, "s-", color=C[1], lw=2, label="Taylor start  kh")
    h = np.array(hA)
    ax[0].loglog(h, eA[0] * (h / h[0]) ** 4, "--", color="#888888", lw=1, label=r"$\propto h^4$")
    ax[0].loglog(h, eT[0] * (h / h[0]) ** 2, ":", color="#888888", lw=1, label=r"$\propto h^2$")
    ax[0].set(xlabel="step h", ylabel=r"max $|\psi-\psi_{exact}|$", title="(a) Global error, free particle")
    ax[0].legend(frameon=False, fontsize=8)
    for i, n in enumerate((1, 2, 3)):
        ax[1].loglog(hB, errB[n], "o-", color=C[i], lw=2, label=f"E_{n}")
    hb = np.array(hB)
    ax[1].loglog(hb, errB[1][0] * (hb / hb[0]) ** 4, "--", color="#888888", lw=1, label=r"$\propto h^4$")
    ax[1].set(xlabel="step h", ylabel="relative error of $E_n$", title="(b) Eigenvalue convergence, box")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25, which="both")
    fig.tight_layout(); fig.savefig("p1_convergence.png", dpi=170)
    print("\nfigures written: p1_wavefunctions.png, p1_convergence.png")


if __name__ == "__main__":
    main()
