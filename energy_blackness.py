"""Energy dependence of the central opacity Gamma(b=0) = 1 - e^{-a} of Eq. (12) (Sec. VII).
Set A follows Eq. (7) with the parameters of Refs. [8, 11]; set B is evolved from its 13 TeV
calibration with the same Delta and alpha'_P.  Usage: python3 energy_blackness.py [output_dir]"""
import json, os, sys, numpy as np
OUT = sys.argv[1] if len(sys.argv) > 1 else 'output'; os.makedirs(OUT, exist_ok=True)
from scipy.optimize import brentq
from regge_eikonal_pp import Ein, GEV2_TO_MB, calibrate_B, params_A, SET_A
D, ap = SET_A['Delta'], SET_A['alphaP']
lamB, nuB = calibrate_B(110.6, 0.281)
def parsB(rs):
    L = 2*np.log(rs/13000.); lam = lamB+ap*L; H = 4*np.pi*lamB*nuB*np.exp(D*L); return lam, H/(4*np.pi*lam)
def parsA(rs):
    H, lam, nu = params_A(rs); return lam, nu
def stot(p, rs): lam, nu = p(rs); return 8*np.pi*lam*Ein(nu/2)*GEV2_TO_MB
def slope(p, rs, h=1e-3): return (np.log(stot(p, rs*np.exp(h)))-np.log(stot(p, rs*np.exp(-h))))/(4*h)
out = {}
for name, p in (('A', parsA), ('B', parsB)):
    s99 = np.exp(brentq(lambda x: 1-np.exp(-p(np.exp(x))[1]/2)-0.99, np.log(13e3), np.log(1e14)))
    out[name] = dict(Gamma0_13TeV=1-np.exp(-p(13e3)[1]/2), sqrts_Gamma0_0p99=float(s99),
                     dlnsig_dlns_13TeV=slope(p, 13e3), dlnsig_dlns_1e6=slope(p, 1e6), dlnsig_dlns_s99=slope(p, s99))
json.dump(out, open(os.path.join(OUT, 'energy_blackness.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
