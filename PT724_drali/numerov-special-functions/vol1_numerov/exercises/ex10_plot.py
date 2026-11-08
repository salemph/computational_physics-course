#!/usr/bin/env python3
"""Plot ex10_wp4.dat / ex10_wp8.dat written by ex10_sp and ex10_dp."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
C = ["#2a78d6", "#eb6834"]
fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, (f, lab) in enumerate((("ex10_wp4.dat", "single precision real(4)"), ("ex10_wp8.dat", "double precision real(8)"))):
    h, e = np.loadtxt(f, unpack=True)
    ee = np.where((e > 0) & (e < 0.1), e, np.nan)       # break the line at exact zeros / failures
    ax.loglog(h, ee, "o-", color=C[i], lw=2, ms=5, label=lab)
    z = e == 0
    if z.any():
        ax.plot(h[z], np.full(z.sum(), 2e-12), "v", color=C[i], ms=7, label="error exactly 0 (lucky rounding)")
hh = np.logspace(-3.7, 0, 50)
ax.loglog(hh, 9.6e-7 * (hh / 0.125) ** 4, ":", color="#888888", label=r"truncation $\propto h^4$")
ax.loglog(hh, 3.7e-9 * (hh / 2.44e-4) ** -2, "--", color="#888888", lw=1, label=r"round-off $\propto \epsilon/h^2$ (double)")
ax.loglog(hh, 3.05e-5 * (hh / 0.0625) ** -2, "-.", color="#888888", lw=1, label=r"round-off $\propto \epsilon/h^2$ (single)")
ax.set(ylim=(1e-12, 1e-1), xlabel="step h", ylabel=r"$|E_0-0.5|$",
       title="Ex. 10: optimum step in single vs double precision")
ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=.25, which="both")
fig.tight_layout(); fig.savefig("ex10.png", dpi=170)
