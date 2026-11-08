#!/usr/bin/env python3
"""Plot the convergence data written by the four Fortran programs (run ./ch1 ... ./ch4 first).
Produces fortran_convergence.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a"]
fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))
D = np.loadtxt("ch1_f_conv.dat")
for i, (n, l) in enumerate(((1, 0), (2, 1), (3, 2))):
    d = D[(D[:, 0] == n) & (D[:, 1] == l)]
    ax[0].loglog(d[:, 2], d[:, 3], "s--", color=C[i], ms=4, lw=1.3, label=f"{n}{'spd'[l]} naive")
    ax[0].loglog(d[:, 2], d[:, 4], "o-", color=C[i], ms=4, lw=2, label=f"{n}{'spd'[l]} Frobenius")
ax[0].set(title="Ch. 1 hydrogen (Fortran)", xlabel="h", ylabel="|ΔE|")
for k, (f, lab, var) in enumerate((("ch2_f_conv.dat", "Ch. 2 angular eq. |Δλ|", "m"),
                                   ("ch3_f_conv.dat", "Ch. 3 circular well |ΔE|", "m"),
                                   ("ch4_f_conv.dat", "Ch. 4 spherical well |ΔE|", "l"))):
    D = np.loadtxt(f)
    a = ax[k + 1]
    names = ("first nodes", "fixed start") if k < 2 else ("naive", "Frobenius")
    for i in range(3):
        d = D[D[:, 0] == i]
        a.loglog(d[:, 2], d[:, 3], "s--", color=C[i], ms=4, lw=1.3, label=f"{var}={i} {names[0]}")
        a.loglog(d[:, 2], d[:, 4], "o-", color=C[i], ms=4, lw=2, label=f"{var}={i} {names[1]}")
    a.set(title=lab + " (Fortran)", xlabel="h")
for a in ax:
    hh = np.array([a.get_xlim()[0] * 1.5, a.get_xlim()[1] / 1.5])
    a.legend(frameon=False, fontsize=6.5)
    a.grid(alpha=.25, which="both"); a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fortran_convergence.png", dpi=170)
print("fortran_convergence.png written")
