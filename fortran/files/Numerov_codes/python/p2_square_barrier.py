#!/usr/bin/env python3
"""
Problem 2 -- RECTANGULAR (SQUARE) POTENTIAL BARRIER with the Numerov method
===========================================================================
Units: hbar = m = 1.  V(x) = V0 for 0 < x < a, 0 elsewhere.
psi''(x) = -k^2(x) psi(x),  k^2(x) = 2 [E - V(x)]   (k^2 < 0 inside barrier if E < V0)

Scattering boundary conditions (wave incident from the left):
    x < 0 :  psi = A e^{ikx} + B e^{-ikx}      (incident + reflected)
    x > a :  psi = F e^{ikx}                    (transmitted only)
Numerical strategy ("backward integration"):
    choose F = 1, impose psi = e^{ikx} at the two right-most grid points,
    integrate the COMPLEX Numerov recursion from right to left,
    at the two left-most grid points decompose psi into A, B.
    T = |F/A|^2 = 1/|A|^2,   R = |B/A|^2.
Verification ("normalisation" of a scattering state = probability-flux
conservation):  R + T = 1.   Accuracy: compare T with the analytic formula.

Discontinuity issue: V jumps at x = 0 and x = a.  Numerov's O(h^4) error
assumes a smooth psi''.  We study two grids:
    (i)  edges ON grid nodes, with V(edge) = V0/2 (mean value)  -> O(h^4)
    (ii) edges falling INSIDE a cell (misaligned grid)            -> O(h^2) or worse
Outputs: tables + p2_wavefunction.png, p2_transmission.png, p2_convergence.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
V0, A_WIDTH = 5.0, 1.0           # barrier height and width


def T_exact(E, V0=V0, a=A_WIDTH):
    """Analytic transmission coefficient of a rectangular barrier."""
    if abs(E - V0) < 1e-12:
        return 1.0 / (1.0 + V0 * a * a / 2.0)           # limit E -> V0  (m = hbar = 1)
    if E < V0:
        kap = np.sqrt(2.0 * (V0 - E))
        return 1.0 / (1.0 + V0 ** 2 * np.sinh(kap * a) ** 2 / (4.0 * E * (V0 - E)))
    q = np.sqrt(2.0 * (E - V0))
    return 1.0 / (1.0 + V0 ** 2 * np.sin(q * a) ** 2 / (4.0 * E * (E - V0)))


def potential(x, h, aligned):
    V = np.where((x > 0.0) & (x < A_WIDTH), V0, 0.0)
    edges = []
    if aligned:                                          # mean value on the edge nodes
        for xe in (0.0, A_WIDTH):
            i = np.where(np.isclose(x, xe, atol=1e-3 * h))[0]
            if len(i):
                V[i[0]] = V0 / 2
                edges.append(i[0])
    return V, edges


def scatter(E, h, mode="corrected", Lpad=2.0, return_wave=False):
    """mode = 'misaligned' : barrier edges fall inside grid cells
              'aligned'    : edges on nodes, V(edge) = V0/2
              'corrected'  : aligned + jump correction term delta at edge nodes"""
    k = np.sqrt(2.0 * E)
    aligned = mode != "misaligned"
    shift = 0.0 if aligned else 0.37 * h                 # edges inside cells if shifted
    xl, xr = -Lpad + shift, A_WIDTH + Lpad + shift
    N = int(round((xr - xl) / h))
    x = xl + h * np.arange(N + 1)
    V, edges = potential(x, h, aligned)
    k2 = 2.0 * (E - V)
    c = 1.0 + h * h * k2 / 12.0
    is_edge = np.zeros(N + 1, dtype=bool)
    is_edge[edges] = True
    jump = {}
    if mode == "corrected":
        for i in edges:   # jump of k^2 across the edge node: right limit - left limit
            jump[i] = 2.0 * (E - V[i + 1]) - 2.0 * (E - V[i - 1])
    psi = np.zeros(N + 1, dtype=complex)
    psi[N], psi[N - 1] = np.exp(1j * k * x[N]), np.exp(1j * k * x[N - 1])   # F = 1
    for n in range(N - 1, 0, -1):                       # backward Numerov
        cp, cm = c[n + 1], c[n - 1]
        if mode == "corrected" and not is_edge[n]:
            # a stencil centred OFF the edge must see the ONE-SIDED value of V
            # at a neighbouring edge node (the value on its own side), not V0/2
            if is_edge[n + 1]:
                cp = 1.0 + h * h * 2.0 * (E - V[n]) / 12.0
            if is_edge[n - 1]:
                cm = 1.0 + h * h * 2.0 * (E - V[n]) / 12.0
        rhs = (12.0 - 10.0 * c[n]) * psi[n] - cp * psi[n + 1]
        if n in jump:
            # stencil centred ON the edge: V0/2 at the centre + jump term
            # delta = -h^3/12 [k^2] psi'(x_n), psi' from the known right side (O(h^2))
            dpsi = (-3.0 * psi[n] + 4.0 * psi[n + 1] - psi[n + 2]) / (2.0 * h)
            rhs += -h ** 3 / 12.0 * jump[n] * dpsi
        psi[n - 1] = rhs / cm
    # decomposition  psi = A e^{ikx} + B e^{-ikx}  at x0, x1
    M = np.array([[np.exp(1j * k * x[0]), np.exp(-1j * k * x[0])],
                  [np.exp(1j * k * x[1]), np.exp(-1j * k * x[1])]])
    Acoef, Bcoef = np.linalg.solve(M, psi[:2])
    T = 1.0 / abs(Acoef) ** 2
    R = abs(Bcoef) ** 2 / abs(Acoef) ** 2
    if return_wave:
        return T, R, x, psi / Acoef                      # incident amplitude -> 1
    return T, R


def slope(h, e):
    return np.polyfit(np.log(h), np.log(e), 1)[0]


def main():
    print("=" * 74)
    print(f"Square barrier  V0 = {V0},  a = {A_WIDTH}   (hbar = m = 1)")
    print("=" * 74)
    h = 0.01
    print(f"{'E':>6} {'T_numerov':>14} {'T_exact':>14} {'|dT|':>10} {'R':>12} {'R+T-1':>11}")
    for E in (0.5, 1.0, 2.5, 4.0, 4.9, 5.1, 6.0, 8.0, 9.935, 12.0, 15.0):
        T, R = scatter(E, h)
        print(f"{E:6.3f} {T:14.10f} {T_exact(E):14.10f} {abs(T-T_exact(E)):10.2e} "
              f"{R:12.8f} {R+T-1:11.2e}")

    # ------------------------------------------------ convergence study
    hs = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625])
    Es = (2.0, 8.0)
    modes = ("misaligned", "aligned", "corrected")
    err = {(E, m): [] for E in Es for m in modes}
    flux = {(E, m): [] for E in Es for m in modes}
    for E in Es:
        for m in modes:
            for hh in hs:
                T, R = scatter(E, hh, mode=m)
                err[E, m].append(abs(T - T_exact(E)))
                flux[E, m].append(abs(R + T - 1.0))
    for E in Es:
        print(f"\nConvergence of T at E = {E}   (T_exact = {T_exact(E):.12f})")
        print(f"{'h':>10} {'misaligned':>12} {'aligned V0/2':>13} {'corrected':>12} {'|R+T-1| corr.':>14}")
        for i, hh in enumerate(hs):
            print(f"{hh:10.6f} {err[E,'misaligned'][i]:12.2e} {err[E,'aligned'][i]:13.2e} "
                  f"{err[E,'corrected'][i]:12.2e} {flux[E,'corrected'][i]:14.1e}")
        print("observed order: " + ", ".join(f"{m} = {slope(hs[:5], err[E, m][:5]):.2f}" for m in modes))

    # ------------------------------------------------ figures
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for i, E in enumerate((2.0, 8.0)):
        T, R, x, psi = scatter(E, 0.005, return_wave=True)
        ax.plot(x, abs(psi) ** 2, color=C[i], lw=2, label=f"E = {E}  (T = {T:.4f})")
    ax.axvspan(0, A_WIDTH, color="#cccccc", alpha=.5, lw=0)
    ax.text(A_WIDTH / 2, ax.get_ylim()[1] * .92, f"barrier V0={V0}", ha="center", fontsize=8)
    ax.set(xlabel="x", ylabel=r"$|\psi(x)|^2$  (incident amplitude 1)",
           title="Probability density across the barrier")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig("p2_wavefunction.png", dpi=170)

    Eg = np.linspace(0.05, 20, 400)
    Tn = np.array([scatter(E, 0.01)[0] for E in Eg])
    Te = np.array([T_exact(E) for E in Eg])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(Eg, Te, color="#888888", lw=1, label="analytic")
    ax[0].plot(Eg[::8], Tn[::8], "o", ms=4, color=C[0], label="Numerov h=0.01")
    ax[0].axvline(V0, color="#888888", ls=":", lw=1)
    ax[0].set(xlabel="E", ylabel="T(E)", title="(a) Transmission coefficient")
    ax[0].legend(frameon=False, fontsize=8, loc="lower right")
    ax[1].semilogy(Eg, abs(Tn - Te) + 1e-18, color=C[1], lw=1.5)
    ax[1].set(xlabel="E", ylabel=r"$|T_{num}-T_{exact}|$", title="(b) Absolute error, h = 0.01")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25)
    fig.tight_layout(); fig.savefig("p2_transmission.png", dpi=170)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
    styles = {"misaligned": ("s--", C[1]), "aligned": ("^-.", C[3]), "corrected": ("o-", C[0])}
    labels = {"misaligned": "edges inside cells", "aligned": "edges on nodes, V=V0/2",
              "corrected": "edges on nodes + jump correction"}
    for j, E in enumerate(Es):
        for m in modes:
            st, col = styles[m]
            ax[j].loglog(hs, err[E, m], st, color=col, lw=2, ms=5, label=labels[m])
        e0 = err[E, "corrected"][0]
        ax[j].loglog(hs, e0 * (hs / hs[0]) ** 4, ":", color="#888888", lw=1, label=r"$\propto h^4$")
        ax[j].loglog(hs, err[E, "aligned"][0] * (hs / hs[0]) ** 2, "-", color="#bbbbbb", lw=1, label=r"$\propto h^2$")
        ax[j].set(xlabel="step h", title=f"E = {E}" + ("  (tunnelling)" if E < V0 else "  (above barrier)"))
        ax[j].spines[["top", "right"]].set_visible(False); ax[j].grid(alpha=.25, which="both")
    ax[0].set_ylabel(r"$|T_{num}-T_{exact}|$")
    ax[0].legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig("p2_convergence.png", dpi=170)
    print("\nfigures written: p2_wavefunction.png, p2_transmission.png, p2_convergence.png")


if __name__ == "__main__":
    main()
