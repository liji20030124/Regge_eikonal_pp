"""One-parameter curve of the single-channel Gaussian eikonal, Eqs. (21) and (29) (Sec. VIII).
Usage: python3 universal_curves.py [output_dir]

At fixed energy the model has two parameters, lambda and a = nu/2, and lambda is
only a length scale.  Two dimensionless combinations of the observables of
Eqs. (16)-(18) therefore depend on a alone:
    r(a) = sigma_el/sigma_tot = 1 - Ein(2a)/(2 Ein(a))              [Eq. (21)]
    W(a) = sigma_tot |t_0| /(8 pi) = Ein(a) Q_0(a),   F(Q_0(a), a) = 0          [Eq. (29)]
The script compares both with the TOTEM values at 13 TeV.
"""
import json, os, sys, numpy as np
OUT = sys.argv[1] if len(sys.argv) > 1 else 'output'; os.makedirs(OUT, exist_ok=True)
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from regge_eikonal_pp import Ein, F, GEV2_TO_MB

def r_of(a): return 1 - Ein(2*a)/(2*Ein(a))
def Q0_of(a):
    qs = np.linspace(0.01, 80, 8000); f = np.array([F(q, a, nmax=400) for q in qs])
    i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    return brentq(lambda q: F(q, a, nmax=400), qs[i[0]], qs[i[0]+1], xtol=1e-13)
def W_of(a): return Ein(a)*Q0_of(a)

# TOTEM 13 TeV: sigma_tot [Ref. 13, Eq. (8)], r [Ref. 13, Table 5], t_dip [Ref. 14]
stot, dstot = 110.6, 3.4
r_m, dr = 0.281, 0.009
tdip, dtdip = 0.47, float(np.hypot(0.004, 0.01))
W_m = stot/GEV2_TO_MB*tdip/(8*np.pi); dW = W_m*np.hypot(dstot/stot, dtdip/tdip)

a_r = brentq(lambda a: r_of(a)-r_m, 1e-3, 60, xtol=1e-12)
from scipy.optimize import minimize_scalar
from scipy.special import jn_zeros
mx = minimize_scalar(lambda la: -W_of(np.exp(la)), bounds=(np.log(2), np.log(30)), method='bounded')
a_max, W_max = float(np.exp(mx.x)), float(-mx.fun)
W_disk = float(jn_zeros(1, 1)[0]**2/4)   # black disk: sigma_tot = 2 pi R^2, |t_0| = j_{1,1}^2/R^2
out = dict(r_meas=[r_m, dr], W_meas=[W_m, dW], a_from_r=a_r, W_at_a_r=W_of(a_r),
           pull_W=(W_of(a_r)-W_m)/dW, a_at_W_max=a_max, W_max=W_max, r_at_W_max=r_of(a_max),
           pull_W_max=(W_max-W_m)/dW, W_black_disk=W_disk,
           W_small_a=W_of(0.05), W_large_a=W_of(20.0), r_small_a=r_of(0.05), r_large_a=r_of(20.0),
           a_band_r=[brentq(lambda a: r_of(a)-(r_m-dr), 1e-3, 60), brentq(lambda a: r_of(a)-(r_m+dr), 1e-3, 60)])
out['W_band_r'] = [W_of(x) for x in out['a_band_r']]
A = np.exp(np.linspace(np.log(0.05), np.log(20), 70))
np.savetxt(os.path.join(OUT, 'universal_curve.csv'), np.array([(a, r_of(a), W_of(a)) for a in A]), delimiter=',',
           header='a,r,W', comments='')
json.dump(out, open(os.path.join(OUT, 'universal_curves.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
