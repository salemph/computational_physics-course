#!/usr/bin/env python3
"""Illustrations of the Numerov method used in the course notes."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

# ---------------------------------------------------------- 1. stencil
fig, ax = plt.subplots(figsize=(9, 3.4))
x = np.linspace(-1.6, 1.6, 300)
f = lambda t: 0.6 * np.sin(1.7 * t + 0.4) + 0.15 * t
ax.plot(x, f(x) + 1.3, color="#888888", lw=1.5)
for i, (xi, lab, w) in enumerate(zip((-1, 0, 1), (r"$x_{n-1}$", r"$x_n$", r"$x_{n+1}$"), ("1", "10", "1"))):
    col = C[1] if i == 2 else C[0]
    ax.plot([xi, xi], [0, f(xi) + 1.3], ":", color="#aaaaaa")
    ax.plot(xi, f(xi) + 1.3, "o", ms=11, color=col, mec="white", mew=2, zorder=5)
    ax.text(xi, -0.18, lab, ha="center", va="top", fontsize=12)
    ax.text(xi, f(xi) + 1.62, [r"$\psi_{n-1}$ (known)", r"$\psi_n$ (known)", r"$\psi_{n+1}$ (new)"][i],
            ha="center", fontsize=10)
    ax.text(xi, 0.25, f"weight of $k^2\\psi$:  {w}", ha="center", fontsize=9, color="#555555")
ax.annotate("", xy=(0, 0.02), xytext=(-1, 0.02), arrowprops=dict(arrowstyle="<->", color="#555555"))
ax.text(-0.5, 0.06, "h", ha="center", fontsize=11)
ax.annotate("", xy=(1, 0.02), xytext=(0, 0.02), arrowprops=dict(arrowstyle="<->", color="#555555"))
ax.text(0.5, 0.06, "h", ha="center", fontsize=11)
ax.add_patch(FancyArrowPatch((0.15, f(0) + 1.3), (0.88, f(1) + 1.3), connectionstyle="arc3,rad=-.35",
                             arrowstyle="-|>", mutation_scale=15, color=C[1], lw=1.5))
ax.text(1.75, 1.05, r"$\psi_{n+1}=\dfrac{2\left(1-\frac{5h^2}{12}k_n^2\right)\psi_n-\left(1+\frac{h^2}{12}k_{n-1}^2\right)\psi_{n-1}}{1+\frac{h^2}{12}k_{n+1}^2}$",
        fontsize=12, va="center")
ax.set(xlim=(-1.7, 5.6), ylim=(-0.45, 2.15)); ax.axis("off")
ax.set_title("The three-point Numerov stencil: two known values give the next one", fontsize=11, pad=22)
fig.tight_layout(); fig.savefig("fig_stencil.png", dpi=170)

# ---------------------------------------------------------- helpers
def numerov(k2, p0, p1, h):
    c = 1 + h * h * k2 / 12
    p = np.zeros(len(k2)); p[0], p[1] = p0, p1
    for n in range(1, len(k2) - 1):
        p[n + 1] = ((12 - 10 * c[n]) * p[n] - c[n - 1] * p[n - 1]) / c[n + 1]
    return p

# ---------------------------------------------------------- 2. shooting idea (box)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
N = 400; xg = np.linspace(0, 1, N + 1); h = xg[1]
Ex = np.pi ** 2 / 2
for i, (E, lab) in enumerate(((0.8 * Ex, r"$E<E_1$: misses, $\psi(L)>0$"), (Ex, r"$E=E_1$: hits $\psi(L)=0$"),
                              (1.25 * Ex, r"$E>E_1$: misses, $\psi(L)<0$"))):
    p = numerov(np.full(N + 1, 2 * E), 0, 1e-3 * h, h)
    ax[0].plot(xg, p / 1e-3, color=[C[3], C[2], C[1]][i], lw=2, label=lab)
ax[0].axhline(0, color="#888888", lw=.8)
ax[0].plot([1, 1], [-0.1, 0.1], color="k", lw=3)
ax[0].set(xlabel="x", ylabel=r"$\psi(x)$ (arb.)", title="(a) Shooting: integrate from x=0 with trial E")
ax[0].legend(frameon=False, fontsize=8, loc="lower left")
Es = np.linspace(0.5, 90, 600)
fL = [numerov(np.full(N + 1, 2 * E), 0, 1e-3 * h, h)[-1] / 1e-3 for E in Es]
ax[1].plot(Es, fL, color=C[0], lw=2)
ax[1].axhline(0, color="#888888", lw=.8)
for n in (1, 2, 3, 4):
    ax[1].plot(n * n * Ex, 0, "o", color=C[1], ms=7, mec="white")
ax[1].set(xlabel="trial energy E", ylabel=r"$\psi(L;E)$", title="(b) Eigenvalues = zeros of $\\psi(L;E)$ (bisection)")
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25)
fig.tight_layout(); fig.savefig("fig_shooting.png", dpi=170)

# ---------------------------------------------------------- 3. instability fault (HO)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
L = 6.0; N = 1200; xg = np.linspace(-L, L, N + 1); h = xg[1] - xg[0]
for i, E in enumerate((0.5, 0.5 + 1e-6, 0.5 - 1e-6)):
    p = numerov(2 * E - xg ** 2, 0, 1e-10, h)
    p /= np.max(np.abs(p[: N // 2]))
    ax[0].plot(xg, p, color=[C[2], C[1], C[0]][i], lw=2,
               label=["E = 0.5 (exact, up to round-off)", "E = 0.5 + 1e-6", "E = 0.5 - 1e-6"][i])
ax[0].set(ylim=(-3, 3), xlabel="x", ylabel=r"$\psi$ (scaled)",
          title="(a) One-sided shooting into the forbidden region")
ax[0].legend(frameon=False, fontsize=8, loc="upper left")
# two-sided
E = 1.37
c = 1 + h * h * (2 * E - xg ** 2) / 12
m = np.searchsorted(xg, np.sqrt(2 * E))
pL = numerov(2 * E - xg ** 2, 0, 1e-10, h)[: m + 2]
pR = numerov((2 * E - xg ** 2)[::-1], 0, 1e-10, h)[::-1]
pL /= pL[m]; pR /= pR[m]
ax[1].plot(xg[: m + 2], pL, color=C[0], lw=2, label="outward from -L")
ax[1].plot(xg[m - 1:], pR[m - 1:], color=C[1], lw=2, label="inward from +L")
ax[1].axvline(xg[m], color="#888888", ls=":", lw=1)
ax[1].text(xg[m] + 0.1, -1.6, "matching node m\n(turning point)", fontsize=8)
ax[1].set(ylim=(-2, 2), xlabel="x", title="(b) Cure: two-sided shooting, kink at m unless E = E_n\n(here E = 1.37, not an eigenvalue)")
ax[1].legend(frameon=False, fontsize=8, loc="upper left")
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=.25)
fig.tight_layout(); fig.savefig("fig_instability.png", dpi=170)

# ---------------------------------------------------------- 4. potentials of the 3 problems
fig, ax = plt.subplots(1, 3, figsize=(11, 2.8))
xx = np.linspace(0, 1, 10)
ax[0].plot([0, 0, 1, 1], [3, 0, 0, 3], color="k", lw=2)
ax[0].set(title="1. Free particle (V=0)\nin a box of length L", xticks=[0, 1], xticklabels=["0", "L"], yticks=[])
xb = [-2, 0, 0, 1, 1, 3]; vb = [0, 0, 5, 5, 0, 0]
ax[1].plot(xb, vb, color="k", lw=2)
ax[1].annotate("", xy=(-0.2, 2.5), xytext=(-1.8, 2.5), arrowprops=dict(arrowstyle="-|>", color=C[0]))
ax[1].text(-1.8, 2.9, "incident", color=C[0], fontsize=8)
ax[1].annotate("", xy=(2.9, 2.5), xytext=(1.3, 2.5), arrowprops=dict(arrowstyle="-|>", color=C[2]))
ax[1].text(1.4, 2.9, "transmitted", color=C[2], fontsize=8)
ax[1].set(title="2. Square barrier\n$V_0$ for 0<x<a", xticks=[0, 1], xticklabels=["0", "a"], yticks=[5], yticklabels=["$V_0$"])
xh = np.linspace(-3.5, 3.5, 200)
ax[2].plot(xh, xh ** 2 / 2, color="k", lw=2)
for n in range(4):
    ax[2].axhline(n + .5, xmin=.5 - np.sqrt(2 * n + 1) / 7, xmax=.5 + np.sqrt(2 * n + 1) / 7, color=C[n], lw=1.5)
ax[2].set(title="3. Harmonic oscillator\n$V=x^2/2$, $E_n=n+1/2$", ylim=(0, 5), yticks=[])
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig_potentials.png", dpi=170)
print("done")
