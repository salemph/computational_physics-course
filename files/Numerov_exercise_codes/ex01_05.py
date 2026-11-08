#!/usr/bin/env python3
"""Solutions of exercises 1-5 of the Numerov course (hbar = m = 1)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def numerov(k2, p0, p1, h):
    c = 1.0 + h * h * k2 / 12.0
    p = np.zeros(len(k2), dtype=np.result_type(p0, p1, float))
    p[0], p[1] = p0, p1
    for n in range(1, len(k2) - 1):
        p[n + 1] = ((12.0 - 10.0 * c[n]) * p[n] - c[n - 1] * p[n - 1]) / c[n + 1]
    return p


def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")


# =================================================================== Ex 1
def ex1():
    print("\n### Exercise 1: order from two grids (Problem 1, Part A)")
    k, X = 2.0, 10.0
    Ns = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800]
    e = []
    for N in Ns:
        h = X / N; x = np.arange(N + 1) * h
        p = numerov(np.full(N + 1, k * k), 0.0, np.sin(k * h), h)
        e.append(np.max(np.abs(p - np.sin(k * x))))
    print(f"{'N -> 2N':>14} {'eps(h)':>10} {'eps(h/2)':>10} {'p':>6}")
    for i in range(len(Ns) - 1):
        print(f"{Ns[i]:6d}->{Ns[i+1]:6d} {e[i]:10.3e} {e[i+1]:10.3e} {np.log2(e[i]/e[i+1]):6.2f}")
    # round-off estimate: eps_mach * N^2 (second-order recursion)
    for N, ee in zip(Ns, e):
        print(f"N={N:6d}  error={ee:.2e}   truncation est. {e[3]*(200/N)**4:.2e}   roundoff est. eps*N^2={2.2e-16*N*N:.2e}")
    return Ns, e


# =================================================================== Ex 2
def ex2():
    print("\n### Exercise 2: third-order Taylor start psi_1 = kh - (kh)^3/6")
    k, X = 2.0, 10.0
    Ns = [25, 50, 100, 200, 400, 800, 1600, 3200]
    res = {"exact": [], "Taylor 1st": [], "Taylor 3rd": []}
    hs = []
    for N in Ns:
        h = X / N; x = np.arange(N + 1) * h; hs.append(h)
        for lab, p1 in (("exact", np.sin(k * h)), ("Taylor 1st", k * h),
                        ("Taylor 3rd", k * h - (k * h) ** 3 / 6)):
            p = numerov(np.full(N + 1, k * k), 0.0, p1, h)
            res[lab].append(np.max(np.abs(p - np.sin(k * x))))
    print(f"{'N':>6} {'h':>8} {'exact':>10} {'Taylor1':>10} {'Taylor3':>10} {'ratio3':>7}")
    for i, N in enumerate(Ns):
        r = res["Taylor 3rd"][i - 1] / res["Taylor 3rd"][i] if i else np.nan
        print(f"{N:6d} {hs[i]:8.5f} {res['exact'][i]:10.3e} {res['Taylor 1st'][i]:10.3e} "
              f"{res['Taylor 3rd'][i]:10.3e} {r:7.2f}")
    hs = np.array(hs)
    for lab in res:
        print(f"  order {lab}: {np.polyfit(np.log(hs[:6]), np.log(res[lab][:6]), 1)[0]:.2f}")
    # decomposition: Taylor-3 error minus exact-start error = start-value contribution
    d = np.abs(np.array(res["Taylor 3rd"]) - np.array(res["exact"]))
    print("  start-value contribution |e3-e_exact|:", " ".join(f"{v:.1e}" for v in d))
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for i, (lab, st) in enumerate(zip(res, ("o-", "s-", "^--"))):
        ax.loglog(hs, res[lab], st, color=C[i], lw=2, ms=7 - 2 * i, label=f"start: {lab}")
    ax.loglog(hs, res["exact"][0] * (hs / hs[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
    ax.loglog(hs, res["Taylor 1st"][0] * (hs / hs[0]) ** 2, "--", color="#888888", label=r"$\propto h^2$")
    ax.set(xlabel="h", ylabel=r"max $|\psi-\sin kx|$", title="Ex. 2: effect of the starting value")
    ax.legend(frameon=False, fontsize=8); style(ax)
    fig.tight_layout(); fig.savefig("ex2.png", dpi=170)


# =================================================================== Ex 3
def box_E(n, N, L=1.0):
    h = L / N
    f = lambda E: numerov(np.full(N + 1, 2 * E), 0.0, 1e-3 * h, h)[-1]
    Ex = (n * np.pi / L) ** 2 / 2
    a, b = Ex * (1 - 0.5 / n), Ex * (1 + 0.5 / n)          # brackets only the n-th zero
    fa = f(a)
    assert np.sign(fa) != np.sign(f(b))
    for _ in range(80):
        m = 0.5 * (a + b); fm = f(m)
        if np.sign(fm) == np.sign(fa): a, fa = m, fm
        else: b = m
    return 0.5 * (a + b), Ex


def ex3():
    print("\n### Exercise 3: smallest N giving E_10 with rel. error < 1e-6 (box L=1)")
    n = 10
    Ns, errs = [], []
    for N in range(20, 401, 2):
        E, Ex = box_E(n, N)
        Ns.append(N); errs.append(abs(E - Ex) / Ex)
    Ns, errs = np.array(Ns), np.array(errs)
    Nmin = Ns[np.argmax(errs < 1e-6)]
    lam = 2.0 / n
    print(f"  N_min = {Nmin}, rel.err = {errs[Ns == Nmin][0]:.3e}, previous N={Nmin-2}: {errs[Ns == Nmin-2][0]:.3e}")
    print(f"  wavelength lambda = 2L/n = {lam}, points per wavelength = {Nmin*lam:.1f}, k h = {n*np.pi/Nmin:.4f}")
    # theory: rel error of E ~ (k h)^4 / 240 for the box (from the dispersion relation)
    for N in (40, 100, Nmin, 200):
        E, Ex = box_E(n, N); kh = n * np.pi / N
        print(f"  N={N:4d} kh={kh:.3f} rel.err={abs(E-Ex)/Ex:.3e}  (kh)^4/240={kh**4/240:.3e}")
    # other n: N needed for 1e-6
    print("  points per wavelength needed for 1e-6, other states:")
    for nn in (1, 2, 5, 10, 20):
        for N in range(max(8, 20 * nn), 4000, 2):
            E, Ex = box_E(nn, N)
            if abs(E - Ex) / Ex < 1e-6:
                print(f"    n={nn:3d}: N_min={N:5d}  points/wavelength={N*2/nn:.1f}")
                break
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ppw = Ns * lam
    ax.loglog(ppw, errs, color=C[0], lw=2, label=r"Numerov, $E_{10}$")
    ax.loglog(ppw, (2 * np.pi / ppw) ** 4 / 240, "--", color=C[1], label=r"$(kh)^4/240$")
    ax.axhline(1e-6, color="#888888", ls=":"); ax.axvline(Nmin * lam, color="#888888", ls=":")
    ax.text(Nmin * lam * 1.04, 3e-3, f"{Nmin*lam:.0f} points\nper wavelength", fontsize=8)
    ax.set(xlabel="grid points per wavelength", ylabel="relative error of $E_{10}$",
           title="Ex. 3: resolution needed for $10^{-6}$")
    ax.legend(frameon=False, fontsize=8); style(ax)
    fig.tight_layout(); fig.savefig("ex3.png", dpi=170)


# =================================================================== general scattering solver
def scatter(E, h, Vfun, xl, xr, use_kh=False, edges=(), V0jump=None):
    """Backward Numerov for a general potential.  edges: x-positions of jumps
    (assumed on nodes); for those the corrected scheme of Problem 2 is used."""
    k = np.sqrt(2 * E)
    N = int(round((xr - xl) / h)); x = xl + h * np.arange(N + 1)
    V = Vfun(x)
    eidx = [int(round((xe - xl) / h)) for xe in edges]
    Vl = {i: Vfun(np.array([x[i] - 1e-9]))[0] for i in eidx}
    Vr = {i: Vfun(np.array([x[i] + 1e-9]))[0] for i in eidx}
    for i in eidx:
        V[i] = 0.5 * (Vl[i] + Vr[i])
    c = 1 + h * h * 2 * (E - V) / 12
    isedge = np.zeros(N + 1, bool); isedge[eidx] = True
    if use_kh:   # discrete wave number of the Numerov plane wave
        a = (k * h) ** 2
        kk = np.arccos((1 - 5 * a / 12) / (1 + a / 12)) / h
    else:
        kk = k
    psi = np.zeros(N + 1, complex)
    psi[N], psi[N - 1] = np.exp(1j * kk * x[N]), np.exp(1j * kk * x[N - 1])
    for n in range(N - 1, 0, -1):
        cp, cm = c[n + 1], c[n - 1]
        if not isedge[n]:
            if n + 1 < N + 1 and isedge[n + 1]: cp = 1 + h * h * 2 * (E - V[n]) / 12
            if isedge[n - 1]: cm = 1 + h * h * 2 * (E - V[n]) / 12
        rhs = (12 - 10 * c[n]) * psi[n] - cp * psi[n + 1]
        if isedge[n]:
            jump = 2 * (E - Vr[n]) - 2 * (E - Vl[n])
            dpsi = (-3 * psi[n] + 4 * psi[n + 1] - psi[n + 2]) / (2 * h)
            rhs += -h ** 3 / 12 * jump * dpsi
        psi[n - 1] = rhs / cm
    M = np.array([[np.exp(1j * kk * x[0]), np.exp(-1j * kk * x[0])],
                  [np.exp(1j * kk * x[1]), np.exp(-1j * kk * x[1])]])
    A, B = np.linalg.solve(M, psi[:2])
    # flux: incident current is proportional to the (discrete) velocity, identical on both sides
    return 1 / abs(A) ** 2, abs(B) ** 2 / abs(A) ** 2


def T_square(E, V0, a):
    kap = np.sqrt(2 * (V0 - E))
    return 1 / (1 + V0 ** 2 * np.sinh(kap * a) ** 2 / (4 * E * (V0 - E)))


# =================================================================== Ex 4
def ex4():
    print("\n### Exercise 4: plane waves e^{ikx} vs discrete e^{i k_h x}")
    V0, a = 5.0, 1.0
    Vf = lambda x: np.where((x > 0) & (x < a), V0, 0.0)
    hs = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125])
    out = {}
    for E in (2.0, 8.0):
        Tex = T_square(E, V0, a) if E < V0 else 1 / (1 + V0 ** 2 * np.sin(np.sqrt(2 * (E - V0)) * a) ** 2 / (4 * E * (E - V0)))
        for pad in (2.0, 20.0):
            e1, e2 = [], []
            for h in hs:
                T1, _ = scatter(E, h, Vf, -pad, a + pad, False, edges=(0.0, a))
                T2, _ = scatter(E, h, Vf, -pad, a + pad, True, edges=(0.0, a))
                e1.append(abs(T1 - Tex)); e2.append(abs(T2 - Tex))
            out[E, pad] = (e1, e2)
            print(f" E={E}, padding={pad}:")
            for i, h in enumerate(hs):
                print(f"   h={h:8.5f}  |dT| e^ikx = {e1[i]:.3e}   e^ik_hx = {e2[i]:.3e}   ratio {e1[i]/e2[i]:.2f}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for j, E in enumerate((2.0, 8.0)):
        for i, pad in enumerate((2.0, 20.0)):
            e1, e2 = out[E, pad]
            ax[j].loglog(hs, e1, "o-", color=C[2 * i], lw=2, label=f"$e^{{ikx}}$, padding {pad:g}")
            ax[j].loglog(hs, e2, "s--", color=C[2 * i + 1], lw=2, label=f"$e^{{ik_hx}}$, padding {pad:g}")
        ax[j].set(xlabel="h", title=f"E = {E}"); style(ax[j])
    ax[0].set_ylabel(r"$|T-T_{exact}|$"); ax[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Ex. 4: continuum vs discrete plane waves at the boundaries", fontsize=10)
    fig.tight_layout(); fig.savefig("ex4.png", dpi=170)


# =================================================================== Ex 5
def ex5():
    print("\n### Exercise 5: thick barrier V0=5, a=6, E=1")
    V0, a, E = 5.0, 6.0, 1.0
    kap = np.sqrt(2 * (V0 - E)); k = np.sqrt(2 * E)
    Tex = T_square(E, V0, a)
    Twkb = np.exp(-2 * kap * a)
    Tasym = 16 * E * (V0 - E) / V0 ** 2 * np.exp(-2 * kap * a)
    print(f"  kappa={kap:.5f}, kappa*a={kap*a:.4f}")
    print(f"  T exact          = {Tex:.6e}")
    print(f"  T WKB e^-2ka     = {Twkb:.6e}   ratio exact/WKB = {Tex/Twkb:.4f}")
    print(f"  T asymptotic 16E(V0-E)/V0^2 e^-2ka = {Tasym:.6e}  (prefactor {16*E*(V0-E)/V0**2:.3f})")
    Vf = lambda x: np.where((x > 0) & (x < a), V0, 0.0)
    rows = []
    for h in (0.05, 0.025, 0.0125, 0.00625, 0.003125):
        T, R = scatter(E, h, Vf, -2.0, a + 2.0, edges=(0.0, a))
        rows.append((h, T, abs(T - Tex) / Tex, R, R + T - 1))
        print(f"  h={h:8.5f}  T={T:.8e}  rel.err={abs(T-Tex)/Tex:.2e}   R={R:.17f}  R+T-1={R+T-1:.2e}")
    # growth factor inside the barrier
    print(f"  |A| ~ 1/sqrt(T) = {1/np.sqrt(Tex):.3e}  (amplitude grows by e^(kappa a)={np.exp(kap*a):.2e} across barrier)")
    # thickness scan
    print("  thickness scan, h=0.00625:")
    for aa in (1, 2, 4, 6, 8, 10, 12):
        Vf2 = lambda x, aa=aa: np.where((x > 0) & (x < aa), V0, 0.0)
        T, R = scatter(E, 0.00625, Vf2, -2.0, aa + 2.0, edges=(0.0, float(aa)))
        Te = T_square(E, V0, aa)
        print(f"    a={aa:3d}  T={T:.4e}  exact={Te:.4e}  rel.err={abs(T-Te)/Te:.2e}  WKB/exact={np.exp(-2*kap*aa)/Te:.3f}  R+T-1={R+T-1:.1e}")


if __name__ == "__main__":
    ex1(); ex2(); ex3(); ex4(); ex5()
