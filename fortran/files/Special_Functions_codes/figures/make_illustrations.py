#!/usr/bin/env python3
"""Illustrations for the special-functions course."""
import numpy as np
from scipy.special import jv, yv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


# ---------------------------------------------------------- 1. the start fault
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
r = np.linspace(1e-4, 1, 600)
k = 2.405
ax[0].plot(r, np.sqrt(r) * jv(0, k * r), color=C[0], lw=2, label=r"regular  $\sqrt{r}J_0(kr)\sim r^{1/2}$")
ax[0].plot(r, np.sqrt(r) * yv(0, k * r), color=C[1], lw=2, label=r"irregular $\sqrt{r}Y_0(kr)\sim r^{1/2}\ln r$")
ax[0].set(ylim=(-1.2, 0.8), xlabel="r", title="(a) The two solutions of the m = 0 equation")
ax[0].legend(frameon=False, fontsize=8, loc="lower right")


def numerov_residual(f, k2, h, j):
    """relative local defect of the Numerov stencil applied to the exact solution at node j"""
    x = h * np.array([j - 1, j, j + 1], float)
    c = 1 + h * h * k2(x) / 12
    res = c[2] * f(x[2]) - (12 - 10 * c[1]) * f(x[1]) + c[0] * f(x[0])
    return abs(res) / abs(f(x[1]))


for i, (m, lab) in enumerate(((0, r"$u=\sqrt{r}J_0$ (circular well, m=0)"), (1, r"$u=\sqrt{r}J_1$ (m=1)"))):
    f = lambda x, m=m: np.sqrt(x) * jv(m, k * x)
    k2 = lambda x, m=m: k * k - (m * m - 0.25) / x ** 2
    for h, mk in ((0.02, "o-"), (0.005, "s--")):
        js = np.arange(1, 21)
        ax[1].semilogy(js, [numerov_residual(f, k2, h, j) for j in js], mk, color=C[i], lw=1.5, ms=4,
                       label=f"{lab}, h={h}")
fsph = lambda x: x * (np.sin(k * x) / (k * x) ** 2 - np.cos(k * x) / (k * x))   # r j_1: analytic
k2s = lambda x: k * k - 2 / x ** 2
js = np.arange(1, 21)
ax[1].semilogy(js, [numerov_residual(fsph, k2s, 0.02, j) for j in js], "^:", color=C[2], lw=1.5, ms=4,
               label=r"$u=r\,j_1$ (spherical well), h=0.02")
ax[1].set(xlabel="node index j  (r = j h)", ylabel="relative local defect",
          title="(b) Numerov defect near r = 0 depends on j, not on h")
ax[1].legend(frameon=False, fontsize=6.5)
# panel c: schematic of fixed-r_s start
x = np.linspace(0, 1, 400)
u = np.sqrt(x) * jv(0, k * x)
ax[2].plot(x, u, color="#888888", lw=1.5)
rs = 0.125; h = 0.025
xs = np.arange(0, rs + h / 2, h)
ax[2].plot(xs, np.sqrt(xs) * jv(0, k * xs), "o", color=C[3], ms=6, label="Frobenius series values")
xn = np.arange(rs + 2 * h, 1.0001, h)
ax[2].plot(xn, np.sqrt(xn) * jv(0, k * xn), "o", color=C[0], ms=4, label="Numerov recursion")
ax[2].axvline(rs, color=C[3], ls=":", lw=1)
ax[2].text(rs + 0.01, 0.42, r"$r_s$ (fixed)", fontsize=9, color=C[3])
ax[2].set(xlabel="r", ylabel="u(r)", title="(c) Cure: series on $[0, r_s]$, Numerov beyond")
ax[2].legend(frameon=False, fontsize=8, loc="lower center")
for a in ax: style(a)
fig.tight_layout(); fig.savefig("fig_start_fault.png", dpi=170)

# ---------------------------------------------------------- 2. stability of recurrences
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
n = np.arange(0, 31)
x = 5.0
ax[0].semilogy(n, abs(jv(n, x)), "o-", color=C[0], lw=2, ms=4, label=r"$|J_n(5)|$ (minimal solution)")
ax[0].semilogy(n, abs(yv(n, x)), "s-", color=C[1], lw=2, ms=4, label=r"$|Y_n(5)|$ (dominant solution)")
ax[0].axvline(x, color="#888888", ls=":", lw=1); ax[0].text(x + 0.4, 1e-20, "n = x", fontsize=8)
ax[0].set(xlabel="n", title=r"(a) Two solutions of $f_{n+1}=\frac{2n}{x}f_n-f_{n-1}$")
ax[0].legend(frameon=False, fontsize=8)
ax[1].axis("off")
txt = ("Rule (Gautschi 1967):\n\n"
       "• A three-term recurrence has two independent\n  solutions. If one decays relative to the other\n"
       "  (minimal vs dominant), round-off excites the\n  dominant one.\n\n"
       "• Compute the DOMINANT solution forwards\n  (Y_n, y_l, P_l, L_n are safe upwards).\n\n"
       "• Compute the MINIMAL solution backwards\n  (J_n, j_l for n > x): Miller's algorithm,\n"
       "  start far above n with arbitrary values,\n  normalise with a known sum or value.")
ax[1].text(0.0, 0.98, txt, va="top", fontsize=9.5, family="DejaVu Sans")
style(ax[0])
fig.tight_layout(); fig.savefig("fig_recurrence_stability.png", dpi=170)

# ---------------------------------------------------------- 3. separation of variables map
fig, ax = plt.subplots(figsize=(10, 3.3))
ax.axis("off")
boxes = [(0.02, "Schrödinger equation\n$-\\frac{1}{2}\\nabla^2\\psi+V\\psi=E\\psi$", "#eef4fc"),
         (0.27, "separation of\nvariables", "#fff4dc"),
         (0.52, "Sturm–Liouville\nODE", "#fff4dc"),
         (0.77, "special function", "#dff3ea")]
for xb, t, col in boxes:
    ax.text(xb + 0.1, 0.85, t, ha="center", va="center", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", fc=col, ec="#2a78d6"))
for xa in (0.2, 0.45, 0.7):
    ax.annotate("", xy=(xa + 0.06, 0.85), xytext=(xa, 0.85), arrowprops=dict(arrowstyle="-|>", color="#555555"))
rows = [("spherical, $V(r)=-1/r$", "radial eq. (Coulomb)", r"$L_{n-l-1}^{(2l+1)}(2r/n)$  Laguerre"),
        ("spherical, any $V(r)$", "angular (polar) eq.", r"$P_l^m(\cos\theta)$  Legendre"),
        ("cylindrical, $V=0$ inside", "radial eq. in $\\rho$", r"$J_m(k\rho)$  Bessel"),
        ("spherical, $V=0$ inside", "radial eq. in $r$", r"$j_l(kr)$  spherical Bessel")]
for i, (a, b, c) in enumerate(rows):
    y = 0.6 - 0.17 * i
    ax.text(0.12, y, a, ha="center", fontsize=8.5)
    ax.text(0.62, y, b, ha="center", fontsize=8.5)
    ax.text(0.87, y, c, ha="center", fontsize=8.5, color=C[i])
ax.set_title("Where the four families come from", fontsize=11)
fig.tight_layout(); fig.savefig("fig_overview.png", dpi=170)
print("done")
