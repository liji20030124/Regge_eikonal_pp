"""Figures 4, 6 and 7 of the paper: fig_profile (profiles of sets A and B and the profile of Eq. (22)
reconstructed from the data), fig_kmr (models M1-M4 of Eqs. (27)-(28)) and fig_universal (Eq. (29)).
Run from this directory after universal_curves.py, profile_inversion.py and kmr_fits.py:
    python3 make_new_figures.py [output_dir, default ../output]"""
import json, os, sys, numpy as np, matplotlib
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, ROOT)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "output")
def o(name): return os.path.join(OUT, name)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from regge_eikonal_pp import params_A, calibrate_B
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "cm", "pdf.fonttype": 42, "ps.fonttype": 42})
HBAR_C = 0.1973269804
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(HERE, name + ".pdf")); fig.savefig(os.path.join(HERE, name + ".png"), dpi=200); plt.close(fig)

# Fig. 4: profiles of sets A and B and the profile reconstructed from the data
HA, lA, nA = params_A(13000.0); lB, nB = calibrate_B(110.6, 0.281)
b = np.linspace(0, 4.5, 400) / HBAR_C
fig, ax = plt.subplots(figsize=(5, 3.4))
for l, n, c, lab in [(lA, nA, "C0", "set A"), (lB, nB, "C3", "set B")]:
    Om = n*np.exp(-b**2/(4*l))
    ax.plot(b*HBAR_C, 1-np.exp(-Om/2), c, label=r"$\Gamma(b)$, " + lab)
    ax.plot(b*HBAR_C, 1-np.exp(-Om), c, ls="--", label=r"$G_{\rm inel}(b)$, " + lab)
P = np.loadtxt(o("profile_13TeV.csv"), delimiter=",", skiprows=1)
ax.plot(P[:, 0], P[:, 1], "k", lw=1.4, label=r"Re$\,\Gamma(b)$, TOTEM")
ax.plot(P[:, 0], P[:, 2], "k", lw=1.4, ls="--", label=r"$G_{\rm inel}(b)$, TOTEM")
ax.axhline(1, color="gray", lw=0.6, ls=":")
ax.set_xlabel(r"$b$ [fm]"); ax.set_ylabel("profile"); ax.set_ylim(0, 1.08); ax.set_xlim(0, 4.5)
ax.legend(fontsize=7.5, frameon=False, loc="upper right")
save(fig, "fig_profile")

# Fig. 6: models M1-M4
KF = os.path.join(ROOT, "kmr_fits.py"); src = open(KF).read().split("#        H   b1")[0]; ns = {"__file__": KF}; exec(src, ns)
R = json.load(open(o("kmr_fits.json"))); TR, Y, ERR = ns["tr"], ns["y"], np.hypot(ns["st"], ns["sy"])
tt = np.linspace(0.001, 4.1, 1500)
sty = {"M1": ("C0", "--"), "M2": ("C2", "-."), "M3": ("C1", ":"), "M4": ("C3", "-")}
lab = {"M1": "M1: Pomeron, one channel", "M2": "M2: Pomeron, two channels",
       "M3": "M3: Pomeron + Odderon, one channel", "M4": "M4: Pomeron + Odderon, two channels"}
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
a1.errorbar(TR, Y, yerr=ERR, fmt=".", color="k", ms=2.5, lw=0.6, label="TOTEM 13 TeV, 290 bins")
for k in ["M1", "M2", "M3", "M4"]:
    p = np.load(o(f"kmr_{k}.npy")); two = R[k]["two"]; c, ls = sty[k]
    a1.semilogy(tt, abs(ns["amp_t"](tt, p, two))**2/(16*np.pi)*ns["GEV2MB"], color=c, ls=ls, lw=1.1, label=lab[k])
    a2.semilogy(TR, ns["binned"](p, two)/Y, ".", color=c, ms=2.5, label=k)
a1.set_ylim(1e-8, 1e3); a1.set_ylabel(r"$d\sigma_{\rm el}/dt$ [mb/GeV$^2$]"); a1.legend(fontsize=7.5, frameon=False)
a2.axhline(1, color="k", lw=0.6); a2.set_ylim(0.1, 10); a2.set_xlim(0, 4.1)
a2.set_xlabel(r"$|t|$ [GeV$^2$]"); a2.set_ylabel("bin-averaged model / data")
a2.legend(fontsize=7, frameon=False, ncol=4, loc="lower left")
save(fig, "fig_kmr")

# Fig. 7: one-parameter curve (r(a), W(a)) of the single-channel Gaussian eikonal
U = np.loadtxt(o("universal_curve.csv"), delimiter=",", skiprows=1); J = json.load(open(o("universal_curves.json")))
fig, ax = plt.subplots(figsize=(5, 3.4))
ax.plot(U[:, 1], U[:, 2], "C0", lw=1.3, label=r"Gaussian eikonal, $0.05\leq a\leq 20$")
for a0 in [0.5, 1, 2, 5, 20]:
    i = np.argmin(abs(U[:, 0]-a0)); ax.plot(U[i, 1], U[i, 2], "o", color="C0", ms=3)
    ax.annotate(f"$a={a0:g}$", (U[i, 1], U[i, 2]), textcoords="offset points", xytext=(-4, -12), fontsize=7)
ax.plot(J["r_meas"][0], J["W_at_a_r"], "s", color="C3", ms=4, mfc="none",
        label="set B")
ax.axhline(J["W_black_disk"], color="gray", lw=0.6, ls=":", label=r"black disk, $j_{1,1}^2/4$")
ax.errorbar([J["r_meas"][0]], [J["W_meas"][0]], xerr=[J["r_meas"][1]], yerr=[J["W_meas"][1]], fmt="o", color="k", ms=4,
            capsize=2, label="TOTEM 13 TeV")
ax.set_xlabel(r"$r=\sigma_{\rm el}/\sigma_{\rm tot}$"); ax.set_ylabel(r"$W=\sigma_{\rm tot}|t_0|/(8\pi)$")
ax.set_xlim(0, 0.45); ax.set_ylim(0, 6); ax.legend(fontsize=7.5, frameon=False, loc="lower right")
save(fig, "fig_universal")
