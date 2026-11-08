#!/usr/bin/env python3
"""Solutions of exercises 6-9 of the Numerov course (hbar = m = 1)."""
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal, eigh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "../python")
import p3_harmonic_oscillator as P3          # re-use the course code of Problem 3

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


def numerov(k2, p0, p1, h):
    c = 1.0 + h * h * k2 / 12.0
    p = np.zeros(len(k2)); p[0], p[1] = p0, p1
    for n in range(1, len(k2) - 1):
        p[n + 1] = ((12.0 - 10.0 * c[n]) * p[n] - c[n - 1] * p[n - 1]) / c[n + 1]
    return p


# =================================================================== Ex 6
def T_poschl(E, V0, w):
    k = np.sqrt(2 * E)
    s = 8 * V0 * w * w
    den = np.cosh(np.pi / 2 * np.sqrt(s - 1)) ** 2 if s > 1 else np.cos(np.pi / 2 * np.sqrt(1 - s)) ** 2
    return np.sinh(np.pi * k * w) ** 2 / (np.sinh(np.pi * k * w) ** 2 + den)


def scatter_smooth(E, h, Vf, X):
    k = np.sqrt(2 * E)
    N = int(round(2 * X / h)); x = -X + h * np.arange(N + 1)
    c = 1 + h * h * 2 * (E - Vf(x)) / 12
    psi = np.zeros(N + 1, complex)
    psi[N], psi[N - 1] = np.exp(1j * k * x[N]), np.exp(1j * k * x[N - 1])
    for n in range(N - 1, 0, -1):
        psi[n - 1] = ((12 - 10 * c[n]) * psi[n] - c[n + 1] * psi[n + 1]) / c[n - 1]
    M = np.array([[np.exp(1j * k * x[0]), np.exp(-1j * k * x[0])],
                  [np.exp(1j * k * x[1]), np.exp(-1j * k * x[1])]])
    A, B = np.linalg.solve(M, psi[:2])
    return 1 / abs(A) ** 2, abs(B) ** 2 / abs(A) ** 2


def ex6():
    print("\n### Exercise 6: smooth barrier V0/cosh^2(x/w), V0=5, w=0.5")
    V0, w, X = 5.0, 0.5, 12.0
    Vf = lambda x: V0 / np.cosh(x / w) ** 2
    print(f"  8 V0 w^2 = {8*V0*w*w}; V(+-{X}) = {Vf(X):.1e} (negligible)")
    for E in (1.0, 3.0, 5.0, 8.0):
        T, R = scatter_smooth(E, 0.01, Vf, X)
        print(f"  E={E}: T={T:.10f}  exact={T_poschl(E,V0,w):.10f}  |dT|={abs(T-T_poschl(E,V0,w)):.2e}  R+T-1={R+T-1:.1e}")
    hs = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125])
    errs = {}
    for E in (2.0, 8.0):
        errs[E] = [abs(scatter_smooth(E, h, Vf, X)[0] - T_poschl(E, V0, w)) for h in hs]
        print(f"  E={E}: errors", " ".join(f"{e:.2e}" for e in errs[E]),
              f" ratios", " ".join(f"{errs[E][i]/errs[E][i+1]:.1f}" for i in range(len(hs) - 1)),
              f" order={np.polyfit(np.log(hs[:5]), np.log(errs[E][:5]), 1)[0]:.2f}")
    # square barrier of the same area for comparison (naive scheme, edges on nodes, V0/2)
    Eg = np.linspace(0.05, 15, 300)
    Tn = np.array([scatter_smooth(E, 0.01, Vf, X)[0] for E in Eg])
    Te = np.array([T_poschl(E, V0, w) for E in Eg])
    print(f"  max |dT| over 0<E<15 at h=0.01: {np.max(abs(Tn-Te)):.2e}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
    ax[0].plot(Eg, Te, color="#888888", lw=1, label="exact (Pöschl–Teller)")
    ax[0].plot(Eg[::10], Tn[::10], "o", ms=4, color=C[0], label="Numerov h = 0.01")
    ax[0].set(xlabel="E", ylabel="T(E)", title=r"(a) $V=V_0/\cosh^2(x/w)$, $V_0=5$, $w=0.5$")
    ax[0].legend(frameon=False, fontsize=8)
    for i, E in enumerate((2.0, 8.0)):
        ax[1].loglog(hs, errs[E], "o-", color=C[i], lw=2, label=f"E = {E}")
    ax[1].loglog(hs, errs[8.0][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax[1].set(xlabel="h", ylabel=r"$|T-T_{exact}|$", title="(b) plain Numerov, no correction: order 4")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax: style(a)
    fig.tight_layout(); fig.savefig("ex6.png", dpi=170)


# =================================================================== Ex 7
def ex7():
    print("\n### Exercise 7: one-sided shooting for the oscillator ground state, h = 0.01")
    h = 0.01
    rows = []
    for L in (3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 38):
        N = int(round(2 * L / h)); x = -L + h * np.arange(N + 1)
        def F(E):
            with np.errstate(all="ignore"):
                return numerov(2 * E - x * x, 0.0, 1e-10, h)[-1]
        a, b = 0.3, 0.7
        fa, fb = F(a), F(b)
        if not (np.isfinite(fa) and np.isfinite(fb)) or np.sign(fa) == np.sign(fb):
            rows.append((L, np.nan, "no sign change / overflow"))
            print(f"  L={L:3d}: psi_N(0.3)={fa:.2e}, psi_N(0.7)={fb:.2e} -> no usable bracket")
            continue
        for _ in range(100):
            m = 0.5 * (a + b); fm = F(m)
            if np.sign(fm) == np.sign(fa): a, fa = m, fm
            else: b = m
        E = 0.5 * (a + b)
        # sensitivity: how big is psi_N for dE = 1e-12 compared with its round-off noise?
        rows.append((L, E, ""))
        print(f"  L={L:3d}: E0 = {E:.15f}   |E0-0.5| = {abs(E-0.5):.2e}   e^(L^2/2)={np.exp(L*L/2):.1e}")
    # wavefunction obtained by one-sided shooting at the converged energy, L = 8
    L = 8.0; N = int(round(2 * L / h)); x = -L + h * np.arange(N + 1)
    E8 = [r[1] for r in rows if r[0] == 8][0]
    p = numerov(2 * E8 - x * x, 0.0, 1e-10, h)
    p /= np.max(abs(p[x < 0])) / np.pi ** -0.25          # scale so that psi(0) ~ exact
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    Ls = [r[0] for r in rows if np.isfinite(r[1])]
    es = [max(abs(r[1] - 0.5), 1e-17) for r in rows if np.isfinite(r[1])]
    ax[0].semilogy(Ls, es, "o-", color=C[0], lw=2)
    ax[0].axhline(3.9e-11, color="#888888", ls=":")
    ax[0].text(9, 1e-10, "h-discretisation level (h = 0.01)", fontsize=8, color="#555555")
    ax[0].set(xlabel="L (box [-L, L])", ylabel=r"$|E_0-0.5|$", title="(a) One-sided shooting: eigenvalue error")
    ax[1].semilogy(x, np.abs(p) + 1e-300, color=C[1], lw=2, label="one-sided Numerov at converged $E_0$")
    ax[1].semilogy(x, np.pi ** -0.25 * np.exp(-x * x / 2), ":", color="k", lw=1.5, label=r"exact $\pi^{-1/4}e^{-x^2/2}$")
    ax[1].set(ylim=(1e-16, 1e3), xlabel="x", ylabel=r"$|\psi(x)|$", title="(b) ...but the right tail is destroyed (L = 8)")
    ax[1].legend(frameon=False, fontsize=8, loc="lower center")
    i_bad = np.where((x > 0) & (np.abs(p) > 10 * np.pi ** -0.25 * np.exp(-x * x / 2)))[0]
    print(f"  L=8: one-sided psi departs from exact (factor 10) beyond x = {x[i_bad[0]]:.2f}; psi(L) = {p[-1]:.2e}")
    for a in ax: style(a)
    fig.tight_layout(); fig.savefig("ex7.png", dpi=170)


# =================================================================== Ex 8
def ho_basis_energies(lam, M=160):
    """Independent check: diagonalise H in the harmonic-oscillator basis."""
    n = np.arange(M + 4)
    X = np.diag(np.sqrt((n[1:]) / 2.0), 1); X = X + X.T
    X4 = np.linalg.matrix_power(X, 4)[:M, :M]
    H = np.diag(np.arange(M) + 0.5) + lam * X4
    return np.linalg.eigvalsh(H)


def ex8():
    print("\n### Exercise 8: anharmonic oscillator V = x^2/2 + lambda x^4")
    orig = P3.V
    lam_list = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]
    table = {}
    for lam in lam_list:
        P3.V = lambda x, lam=lam: 0.5 * x * x + lam * x ** 4
        x, h, st = P3.eigenstates(nmax=3, L=6.0, N=1200, Emax=8.0)
        table[lam] = {n: E for n, E, _ in st}
    P3.V = orig
    print(f"{'lambda':>8} {'n':>2} {'E Numerov':>16} {'E basis':>16} {'dE/lambda':>11} {'1st-order coef':>14} {'1st-order E':>12}")
    for lam in lam_list:
        Eb = ho_basis_energies(lam)
        for n in (0, 1, 2):
            E = table[lam][n]
            c1 = 0.75 * (2 * n * n + 2 * n + 1)
            print(f"{lam:8.0e} {n:2d} {E:16.12f} {Eb[n]:16.12f} {(E-n-0.5)/lam:11.6f} {c1:14.4f} {n+0.5+c1*lam:12.6f}")
    lam = 0.1
    print(f"  lambda=0.1: E0 = {table[0.1][0]:.10f}; first order 0.575, 2nd order 0.575-21/8*0.01 = {0.575-21/8*0.01:.6f},"
          f" 3rd order + 333/16*0.001 = {0.575-0.02625+333/16*0.001:.6f}")
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    lams = np.array(lam_list)
    for i, n in enumerate((0, 1, 2)):
        y = [(table[l][n] - n - 0.5) / l for l in lam_list]
        ax.semilogx(lams, y, "o-", color=C[i], lw=2, label=f"n={n} Numerov")
        ax.axhline(0.75 * (2 * n * n + 2 * n + 1), color=C[i], ls=":", lw=1)
    ax.semilogx(lams, 0.75 - 21 / 8 * lams, "--", color="#888888", lw=1, label="n=0: 3/4 − (21/8)λ")
    ax.set(xlabel=r"$\lambda$", ylabel=r"$(E_n-E_n^{(0)})/\lambda$",
           title="Ex. 8: shift / λ → first-order coefficient (dotted) as λ→0")
    ax.legend(frameon=False, fontsize=8); style(ax)
    fig.tight_layout(); fig.savefig("ex8.png", dpi=170)


# =================================================================== Ex 9
V0W, AW = 10.0, 2.0              # well depth and width


def well_exact():
    """Exact bound states from q tan(qa/2) = kappa (even), -q cot(qa/2) = kappa (odd)."""
    res = []
    f_even = lambda E: np.sqrt(2 * (V0W + E)) * np.tan(np.sqrt(2 * (V0W + E)) * AW / 2) - np.sqrt(-2 * E)
    f_odd = lambda E: -np.sqrt(2 * (V0W + E)) / np.tan(np.sqrt(2 * (V0W + E)) * AW / 2) - np.sqrt(-2 * E)
    Eg = np.linspace(-V0W + 1e-9, -1e-9, 200001)
    for f, par in ((f_even, "even"), (f_odd, "odd")):
        with np.errstate(all="ignore"):
            v = f(Eg)
        for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
            if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and abs(v[i]) < 50 and abs(v[i + 1]) < 50:
                res.append((brentq(f, Eg[i], Eg[i + 1], xtol=1e-15), par))
    return sorted(res)


def well_solve(h, mode="corrected", L=8.0, Emin=None):
    """Two-sided shooting for the finite well; mode in {'misaligned','aligned','corrected'}."""
    shift = 0.37 * h if mode == "misaligned" else 0.0
    N = int(round(2 * L / h)); x = -L + shift + h * np.arange(N + 1)
    inside = np.abs(x) < AW / 2
    V = np.where(inside, -V0W, 0.0)
    edges = []
    if mode != "misaligned":
        for xe in (-AW / 2, AW / 2):
            i = int(np.argmin(abs(x - xe))); V[i] = -V0W / 2; edges.append(i)
    Vl = {edges[0]: 0.0, edges[1]: -V0W} if edges else {}
    Vr = {edges[0]: -V0W, edges[1]: 0.0} if edges else {}
    isedge = np.zeros(N + 1, bool); isedge[edges] = True
    m = int(np.argmin(abs(x)))                    # match at the centre of the well
    corr = mode == "corrected"

    def shoot(E):
        k2 = 2 * (E - V); c = 1 + h * h * k2 / 12
        def cc(j, centre):                        # weight of node j as seen from stencil 'centre'
            if corr and isedge[j] and not isedge[centre]:
                return 1 + h * h * 2 * (E - V[centre]) / 12
            return c[j]
        pL = np.zeros(N + 1); pL[1] = 1e-10
        for n in range(1, m + 1):
            rhs = (12 - 10 * c[n]) * pL[n] - cc(n - 1, n) * pL[n - 1]
            if corr and isedge[n]:
                dpsi = (3 * pL[n] - 4 * pL[n - 1] + pL[n - 2]) / (2 * h)   # from the left (known) side
                rhs += -h ** 3 / 12 * (2 * (E - Vr[n]) - 2 * (E - Vl[n])) * dpsi
            pL[n + 1] = rhs / cc(n + 1, n)
        pR = np.zeros(N + 1); pR[N - 1] = 1e-10
        for n in range(N - 1, m - 1, -1):
            rhs = (12 - 10 * c[n]) * pR[n] - cc(n + 1, n) * pR[n + 1]
            if corr and isedge[n]:
                dpsi = (-3 * pR[n] + 4 * pR[n + 1] - pR[n + 2]) / (2 * h)  # from the right side
                rhs += -h ** 3 / 12 * (2 * (E - Vr[n]) - 2 * (E - Vl[n])) * dpsi
            pR[n - 1] = rhs / cc(n - 1, n)
        pL /= np.max(abs(pL[:m + 2])); pR /= np.max(abs(pR[m - 1:]))
        D = pL[m + 1] * pR[m] - pL[m] * pR[m + 1]
        psi = np.concatenate([pL[:m + 1], pR[m + 1:] * pL[m + 1] / pR[m + 1]])  # glue (pR[m] may vanish for odd states)
        return D, psi

    Es = np.arange(-V0W + 0.01, -0.005, 0.01)
    Ds = np.array([shoot(E)[0] for E in Es])
    out = []
    for i in np.where(np.sign(Ds[:-1]) != np.sign(Ds[1:]))[0]:
        a, b, fa = Es[i], Es[i + 1], Ds[i]
        for _ in range(70):
            mm = 0.5 * (a + b); fm = shoot(mm)[0]
            if np.sign(fm) == np.sign(fa): a, fa = mm, fm
            else: b = mm
        E = 0.5 * (a + b)
        psi = shoot(E)[1]
        psi /= np.sqrt(np.sum(psi ** 2) * h)
        out.append((E, x, psi))
    return out


def ex9():
    print(f"\n### Exercise 9: finite square well V=-{V0W} for |x|<{AW/2}")
    ex = well_exact()
    for E, par in ex:
        q, kap = np.sqrt(2 * (V0W + E)), np.sqrt(-2 * E)
        lhs = q * np.tan(q * AW / 2) if par == "even" else -q / np.tan(q * AW / 2)
        print(f"  exact {par:4s}: E = {E:.12f}  q={q:.6f} kappa={kap:.6f}  check {lhs:.10f} = {kap:.10f}")
    st = well_solve(0.01)
    print("  Numerov, corrected scheme, h=0.01:")
    for (E, x, psi), (Ee, par) in zip(st, ex):
        q, kap = np.sqrt(2 * (V0W + E)), np.sqrt(-2 * E)
        lhs = q * np.tan(q * AW / 2) if par == "even" else -q / np.tan(q * AW / 2)
        print(f"   E = {E:.12f} ({par})  |dE| = {abs(E-Ee):.2e}   residual of exact condition {lhs-kap:.2e}"
              f"  parity check psi(x)={'+' if par=='even' else '-'}psi(-x): {np.max(abs(psi - (1 if par=='even' else -1)*psi[::-1])):.1e}")
    hs = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125])
    errs = {m: [] for m in ("misaligned", "aligned", "corrected")}
    for h in hs:
        for mode in errs:
            r = well_solve(h, mode)
            errs[mode].append(max(abs(r[i][0] - ex[i][0]) for i in range(len(ex))))
    print(f"{'h':>9} {'misaligned':>12} {'aligned':>12} {'corrected':>12}   (max |dE| over all bound states)")
    for i, h in enumerate(hs):
        print(f"{h:9.5f} {errs['misaligned'][i]:12.2e} {errs['aligned'][i]:12.2e} {errs['corrected'][i]:12.2e}")
    for mm in errs:
        print(f"  order {mm}: {np.polyfit(np.log(hs[:5]), np.log(errs[mm][:5]), 1)[0]:.2f}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    xx = np.linspace(-3, 3, 7)
    ax[0].plot([-3, -1, -1, 1, 1, 3], [0, 0, -V0W, -V0W, 0, 0], color="k", lw=1.5)
    for i, (E, x, psi) in enumerate(st):
        sel = abs(x) < 3
        ax[0].axhline(E, color="#dddddd", lw=.8)
        ax[0].plot(x[sel], E + 1.5 * psi[sel], color=C[i], lw=2, label=f"E = {E:.6f}")
    ax[0].set(xlabel="x", ylabel=r"$E + 1.5\,\psi(x)$", title="(a) Bound states of the finite well")
    ax[0].set(ylim=(-11, 2.5))
    ax[0].legend(frameon=False, fontsize=7, loc="upper center", ncol=3)
    st_ = {"misaligned": ("s--", C[1]), "aligned": ("^-.", C[3]), "corrected": ("o-", C[0])}
    for mm in errs:
        ax[1].loglog(hs, errs[mm], st_[mm][0], color=st_[mm][1], lw=2, label=mm)
    ax[1].loglog(hs, errs["corrected"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax[1].set(xlabel="h", ylabel=r"max $|E_n-E_n^{exact}|$", title="(b) Convergence: the correction restores order 4")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax: style(a)
    fig.tight_layout(); fig.savefig("ex9.png", dpi=170)


if __name__ == "__main__":
    which = sys.argv[1:] or ["6", "7", "8", "9"]
    for w in which:
        {"6": ex6, "7": ex7, "8": ex8, "9": ex9}[w]()
