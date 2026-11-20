#!/usr/bin/env python3
"""
Two-channel Good-Walker eikonal (PDG review [Khoze-Ryskin-Tasevsky], Eqs. 20.20-20.26):
  |p> = (|phi_1> + |phi_2>)/sqrt(2),  Omega_ik = gamma_i gamma_k Omega(b),  gamma_{1,2} = 1 +- gamma,
  F_ik = 1 - exp(-Omega_ik/2),  <F> = (1/4) sum_ik F_ik,
  dsigma_tot/d2b = 2 Re<F>,  dsigma_el/d2b = |<F>|^2,  dsigma_(el+low-mass diffr.)/d2b = <|F|^2>.
Omega(b) is the complex single-Gaussian opacity of Eq. (24) (set B_c), so that gamma = 0 reproduces set B_c.
(H, lambda) are calibrated to sigma_tot and r of TOTEM; gamma is fitted (equal weights in ln dsigma/dt)
to the bins with |t| < TFIT.   Usage: python3 goodwalker.py TOTEM13_dsdt.csv outdir
"""
import sys, os, json
import numpy as np
from scipy.special import j0
from scipy.optimize import least_squares, minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extensions import bb, wbb, comps_from, omega, bin_avg_matrix, GEV2_TO_MB, TOTEM13, dip_bump as _db

def profiles(H, lam, g):
    Om = omega(bb, comps_from(H, lam))
    gs = [1 + g, 1 - g]
    F = [[1 - np.exp(-gi * gk * Om / 2) for gk in gs] for gi in gs]
    Fm = 0.25 * (F[0][0] + F[0][1] + F[1][0] + F[1][1])
    F2 = 0.25 * sum(abs(F[i][k])**2 for i in range(2) for k in range(2))
    return Fm, F2

def obs(H, lam, g):
    Fm, F2 = profiles(H, lam, g)
    I = lambda f: 2 * np.pi * np.sum(bb * f * wbb)
    stot = 2 * I(Fm.real); sel = I(abs(Fm)**2); seld = I(F2)
    A0 = 4 * np.pi * np.sum(bb * Fm * wbb)
    mm = Fm.real > 1e-12
    pump = float(np.max(F2[mm] / Fm.real[mm]))               # Pumplin bound: <|F|^2> <= Re<F>
    return dict(stot=stot * GEV2_TO_MB, sel=sel * GEV2_TO_MB, r=sel / stot, s_lowmass_diff=(seld - sel) * GEV2_TO_MB,
                rho=float(-A0.imag / A0.real), pumplin_max_ratio=pump)

_K = {}
def dsdt_t(t, H, lam, g, key=None):
    Fm, _ = profiles(H, lam, g)
    if key is None or key not in _K:
        K = 4 * np.pi * j0(np.sqrt(np.ravel(t))[:, None] * bb[None, :]) * (bb * wbb)[None, :]
        if key is not None: _K[key] = K
    else:
        K = _K[key]
    return np.abs(K @ Fm)**2 / (16 * np.pi) * GEV2_TO_MB

def calib(g, x0):
    def r(p):
        o = obs(p[0], p[1], g)
        return [o['stot'] / TOTEM13['stot'] - 1, o['r'] / TOTEM13['r'] - 1]
    return least_squares(r, x0, xtol=1e-13, ftol=1e-13, bounds=([1, 0.5], [20000, 50])).x

def main(csv, outdir, tfit):
    os.makedirs(outdir, exist_ok=True)
    d = np.loadtxt(csv, delimiter=',', skiprows=1); lo, hi, tr, y, st, sy = d.T
    T, W = bin_avg_matrix(lo, hi); msk = tr < tfit
    state = {'x': np.array([463.7, 7.516])}
    def binned(H, lam, g):
        return (dsdt_t(T, H, lam, g, key='T').reshape(T.shape) * W).sum(axis=1)
    def cost(g):
        H, lam = calib(g, state['x']); state['x'] = np.array([H, lam])
        m = binned(H, lam, g)
        return float(np.sum((np.log(m[msk]) - np.log(y[msk]))**2))
    scan = [(g, cost(g)) for g in np.linspace(0.0, 0.8, 17)]
    g0 = min(scan, key=lambda x: x[1])[0]
    res = minimize_scalar(cost, bounds=(max(0, g0 - 0.05), min(0.95, g0 + 0.05)), method='bounded', options={'xatol': 1e-5})
    g = float(res.x); H, lam = calib(g, state['x'])
    o = obs(H, lam, g); m = binned(H, lam, g)
    tt = np.linspace(0.1, 1.5, 7001); v = dsdt_t(tt, H, lam, g)
    i = np.where((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:]))[0]
    if len(i):
        k = i[0] + 1; j = np.where((v[k+1:-1] > v[k:-2]) & (v[k+1:-1] > v[k+2:]))[0]; kb = k + 1 + j[0]
        db = [float(tt[k]), float(tt[kb]), float(v[kb] / v[k])]
    else:
        db = None
    kd = int(np.argmin(abs(tr - 0.46799)))
    sl = (tr >= 0.04) & (tr <= 0.2)
    R = dict(gamma=g, H=H, lam=lam, nu=H / (4 * np.pi * lam), **o, dipbump=db, ratio_dip=float(m[kd] / y[kd]),
             rms_below1=float(np.sqrt(np.mean(np.log(m[tr < 1] / y[tr < 1])**2))), rms_all=float(np.sqrt(np.mean(np.log(m / y)**2))),
             frac10_below1=float(np.mean(abs(m[tr < 1] / y[tr < 1] - 1) < 0.1)), slope=float(-np.polyfit(tr[sl], np.log(m[sl]), 1)[0]),
             ratio_above1=[float((m / y)[tr > 1].min()), float((m / y)[tr > 1].max())], tfit=tfit, scan=scan)
    json.dump(R, open(os.path.join(outdir, 'goodwalker.json'), 'w'), indent=1)
    tg = np.linspace(0.001, 4.1, 3000)
    np.savetxt(os.path.join(outdir, 'gw_curve.csv'), np.column_stack([tg, dsdt_t(tg, H, lam, g)]), delimiter=',', header='t,D', comments='')
    np.savetxt(os.path.join(outdir, 'gw_bins.csv'), np.column_stack([tr, y, m]), delimiter=',', header='t_repr,TOTEM,D_binavg', comments='')
    return R

def gamma_scan(csv, outdir, gammas=(0.0, 0.2, 0.4, 0.6)):
    """Fixed-gamma scan quoted in Sec. VII: for each gamma, (H, lambda) are calibrated to sigma_tot and r;
    output lambda, nu, low-mass diffractive cross section, rho, Pumplin ratio, position of the first dip,
    model/TOTEM ratio at |t| = 0.468 GeV^2 and the rms of ln(model/TOTEM) for |t| < 1 GeV^2."""
    d = np.loadtxt(csv, delimiter=',', skiprows=1); lo, hi, tr, y, st, sy = d.T
    T, W = bin_avg_matrix(lo, hi); kd = int(np.argmin(abs(tr - 0.46799)))
    tt = np.linspace(0.1, 1.5, 7001); x = np.array([463.7, 7.516]); out = []
    for g in gammas:
        H, lam = calib(g, x); x = np.array([H, lam]); o = obs(H, lam, g)
        m = (dsdt_t(T, H, lam, g).reshape(T.shape) * W).sum(axis=1)
        v = dsdt_t(tt, H, lam, g); i = np.where((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:]))[0]
        out.append(dict(g=float(g), H=float(H), lam=float(lam), nu=float(H / (4 * np.pi * lam)), smd=o['s_lowmass_diff'],
                        rho=o['rho'], pump=o['pumplin_max_ratio'], tmin=float(tt[i[0] + 1]), ratio_dip=float(m[kd] / y[kd]),
                        rms1=float(np.sqrt(np.mean(np.log(m[tr < 1] / y[tr < 1])**2)))))
    json.dump(out, open(os.path.join(outdir, 'gw_gamma_scan.json'), 'w'), indent=1)
    return out

if __name__ == '__main__':
    tfit = float(os.environ.get('TFIT', '1.0'))
    R = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'output', tfit)
    print(json.dumps({k: v for k, v in R.items() if k != 'scan'}, indent=1))
    print('scan', [(round(a, 2), round(b, 2)) for a, b in R['scan']])
    if tfit == 1.0:
        for e in gamma_scan(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'output'):
            print('gamma=%.1f lam=%.3f nu=%.3f s_lowmass=%.3f mb t_min=%.4f ratio(0.468)=%.2f Pumplin=%.3f' %
                  (e['g'], e['lam'], e['nu'], e['smd'], e['tmin'], e['ratio_dip'], e['pump']))
