#!/usr/bin/env python3
"""Bin-averaged comparison of the Gaussian Regge-eikonal model with all 290 TOTEM 13 TeV bins.
Usage: python3 full_comparison.py TOTEM13_dsdt.csv outdir"""
import sys, os, json
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regge_eikonal_pp import params_A, calibrate_B, dsdt, TOTEM13

def bin_average(f, lo, hi, n=24):
    x, w = leggauss(n)
    t = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    return 0.5 * np.sum(w * np.array([f(tt) for tt in t]))

def main(csv, outdir):
    os.makedirs(outdir, exist_ok=True)
    d = np.loadtxt(csv, delimiter=',', skiprows=1)
    lo, hi, tr, y, st, sy = d.T
    err = np.hypot(st, sy)
    HA, lA, nA = params_A(13000.0)
    lB, nB = calibrate_B(TOTEM13['stot'], TOTEM13['r'])
    fA = lambda t: dsdt(t, lA, nA / 2)
    fB = lambda t: dsdt(t, lB, nB / 2)
    out = []
    for i in range(len(d)):
        aA, aB = bin_average(fA, lo[i], hi[i]), bin_average(fB, lo[i], hi[i])
        out.append([lo[i], hi[i], tr[i], y[i], st[i], sy[i], aA, fA(tr[i]), aB, fB(tr[i])])
    out = np.array(out)
    np.savetxt(os.path.join(outdir, 'all_290_bins.csv'), out, delimiter=',', fmt='%.10g',
               header='t_low,t_high,t_repr,TOTEM,stat,syst,A_binavg,A_point,B_binavg,B_point', comments='')
    R = {}
    rA, rB = out[:, 6] / y, out[:, 8] / y
    R['n_bins'] = int(len(d))
    R['max_rel_point_vs_binavg_B'] = float(np.max(np.abs(out[:, 9] / out[:, 8] - 1)))
    # six representative bins
    six = [0.04990, 0.10012, 0.19979, 0.29960, 0.46799, 0.63179]
    R['six'] = []
    for t in six:
        k = int(np.argmin(abs(tr - t)))
        R['six'].append(dict(t=float(tr[k]), lo=float(lo[k]), hi=float(hi[k]), totem=float(y[k]), stat=float(st[k]), syst=float(sy[k]),
                             A_bin=float(out[k, 6]), A_point=float(out[k, 7]), B_bin=float(out[k, 8]), B_point=float(out[k, 9]),
                             A_ratio=float(out[k, 6] / y[k]), B_ratio=float(out[k, 8] / y[k])))
    # log-linear diagnostic slope, equal weights, representative |t| in [0.04, 0.2]
    m = (tr >= 0.04) & (tr <= 0.2)
    R['n_slope_points'] = int(m.sum())
    R['slope_data'] = float(-np.polyfit(tr[m], np.log(y[m]), 1)[0])
    R['slope_B_binavg'] = float(-np.polyfit(tr[m], np.log(out[m, 8]), 1)[0])
    R['slope_A_binavg'] = float(-np.polyfit(tr[m], np.log(out[m, 6]), 1)[0])
    # descriptive agreement statistics (no chi2: correlated syst and normalization not available)
    for key, r in (('A', rA), ('B', rB)):
        R[f'frac_within_10pct_{key}'] = float(np.mean(np.abs(r - 1) < 0.10))
        R[f'max_t_within_5pct_contiguous_{key}'] = float(tr[np.argmax(np.abs(r - 1) >= 0.05)]) if np.any(np.abs(r - 1) >= 0.05) else None
        for a, b in [(0.0, 0.2), (0.2, 0.4), (0.4, 1.0), (1.0, 4.1)]:
            mm = (tr >= a) & (tr < b)
            R[f'ratio_range_{key}_{a}_{b}'] = [float(r[mm].min()), float(r[mm].max())]
    # local data minimum and following maximum
    k_dip = int(np.argmin(np.where((tr > 0.3) & (tr < 0.6), y, np.inf)))
    k_bump = k_dip + int(np.argmax(y[k_dip:k_dip + 40]))
    R['data_min_bin'] = [float(tr[k_dip]), float(y[k_dip])]
    R['data_bump_bin'] = [float(tr[k_bump]), float(y[k_bump])]
    R['data_bump_dip_ratio'] = float(y[k_bump] / y[k_dip])
    json.dump(R, open(os.path.join(outdir, 'full_comparison.json'), 'w'), indent=1)
    return R, out

if __name__ == '__main__':
    R, out = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'output')
    print(json.dumps({k: v for k, v in R.items() if k != 'six'}, indent=1))
    for s in R['six']:
        print('%.5f [%.5f,%.5f] data=%g+-%g+-%g  A bin=%.6g pt=%.6g (%.4f)  B bin=%.6g pt=%.6g (%.4f)' %
              (s['t'], s['lo'], s['hi'], s['totem'], s['stat'], s['syst'], s['A_bin'], s['A_point'], s['A_ratio'], s['B_bin'], s['B_point'], s['B_ratio']))
