#!/usr/bin/env python3
"""Plot the convergence data written by the three Fortran programs.
Run after ./p1 ./p2 ./p3 in the same folder.  Produces fortran_convergence.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a"]
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))

h, eA, eT = np.loadtxt("p1_partA.dat", unpack=True)
ax[0].loglog(h, eA, "o-", color=C[0], lw=2, label="free particle, exact start")
ax[0].loglog(h, eT, "s-", color=C[1], lw=2, label="free particle, Taylor start")
B = np.loadtxt("p1_partB.dat")
ax[0].loglog(B[:, 0], B[:, 1], "^-", color=C[2], lw=2, label="box, rel. error $E_1$")
ax[0].loglog(h, eA[0] * (h / h[0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
ax[0].set(title="Problem 1 (Fortran)", xlabel="h", ylabel="error")

D = np.loadtxt("p2_conv.dat")
for E, mk in ((2.0, "-"), (8.0, "--")):
    d = D[D[:, 0] == E]
    for j, (lab, col) in enumerate(zip(("misaligned", "aligned V0/2", "corrected"), C)):
        ax[1].loglog(d[:, 1], d[:, 2 + j], "o" + mk, color=col, ms=4, lw=1.5,
                     label=f"{lab}" if E == 2.0 else None)
ax[1].set(title="Problem 2 (Fortran): |T-T_exact|\n(solid E=2, dashed E=8)", xlabel="h")

H = np.loadtxt("p3_conv_h.dat")
for j, n in enumerate((0, 2, 4)):
    ax[2].loglog(H[:, 0], H[:, 1 + j], "o-", color=C[j], lw=2, label=f"$E_{n}$")
ax[2].loglog(H[:, 0], H[0, 3] * (H[:, 0] / H[0, 0]) ** 4, ":", color="#888888", label=r"$\propto h^4$")
ax[2].set(title="Problem 3 (Fortran): |E_n - (n+1/2)|", xlabel="h")
for a in ax:
    a.legend(frameon=False, fontsize=7)
    a.grid(alpha=.25, which="both"); a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fortran_convergence.png", dpi=170)
print("fortran_convergence.png written")
