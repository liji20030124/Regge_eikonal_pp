#!/usr/bin/env python3
"""
Extensions of the Gaussian Regge-eikonal model (paper, Sec. VII):
  (1) complex amplitude: signature phase restored,  H -> H_c = H (1 - i tan(pi Delta/2)),
      lambda -> lambda_c = lambda - i pi alpha'/2                         [Poghosyan p.22; PDG Eqs. 20.5-20.7]
  (2) two-component opacity: A_B = i s sum_j H_cj exp(lambda_cj t),  j = 1, 2 (model assumption)
The eikonal amplitude is evaluated numerically in impact-parameter space,
  Omega(b) = sum_j H_cj/(4 pi lambda_cj) exp(-b^2/(4 lambda_cj)),
  A/(i s) = 2 * 2 pi int b db J0(q b) [1 - exp(-Omega/2)],
and checked against the closed-form series for a single real Gaussian.
Usage: python3 extensions.py TOTEM13_dsdt.csv outdir
"""
import sys, os, json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0
from scipy.optimize import least_squares, brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regge_eikonal_pp import GEV2_TO_MB, params_A, calibrate_B, F, Ein, TOTEM13

DELTA, ALPHAP = 0.0808, 0.25          # PDG review Eq. (20.7)
TAU = np.tan(np.pi * DELTA / 2)       # Born real-to-imaginary ratio

# ---------------------------------------------------------------- b-space machinery
NB = 3000
xb, wb = leggauss(NB)
BMAX = 60.0                            # GeV^-1  (Omega < 1e-30 there for all models used)
bb = 0.5 * BMAX * (xb + 1); wbb = 0.5 * BMAX * wb

def omega(b, comps):
    """comps: list of (H_c, lambda_c) complex Born components."""
    return sum(H / (4 * np.pi * l) * np.exp(-b**2 / (4 * l)) for H, l in comps)

def amp(t, comps):
    """A(s,t)/(i s) [GeV^-2] for array t >= 0."""
    t = np.atleast_1d(t)
    G = 1 - np.exp(-omega(bb, comps) / 2)
    K = j0(np.sqrt(t)[:, None] * bb[None, :])
    return 4 * np.pi * (K * (bb * G * wbb)[None, :]).sum(axis=1)

def observables(comps):
    G = 1 - np.exp(-omega(bb, comps) / 2)
    Om = omega(bb, comps)
    stot = 2 * 2 * np.pi * np.sum(bb * G.real * wbb)
    sel = 2 * np.pi * np.sum(bb * np.abs(G)**2 * wbb)
    sinel = 2 * np.pi * np.sum(bb * (1 - np.exp(-Om.real)) * wbb)
    A0 = amp(np.array([0.0]), comps)[0]       # A/(i s): Re -> Im A/s, Im -> -Re A/s
    rho = -A0.imag / A0.real
    h = 1e-5
    a1 = amp(np.array([h]), comps)[0]
    B0 = (np.log(abs(A0)**2) - np.log(abs(a1)**2)) / h
    return dict(stot=stot * GEV2_TO_MB, sel=sel * GEV2_TO_MB, sinel=sinel * GEV2_TO_MB,
                r=sel / stot, rho=rho, B0=B0, unitarity=float(abs(stot - sel - sinel) / stot))

def dsdt(t, comps):
    return np.abs(amp(t, comps))**2 / (16 * np.pi) * GEV2_TO_MB      # |A|^2/(16 pi s^2)

def comps_from(H, lam, phase=True):
    if phase:
        return [(H * (1 - 1j * TAU), lam - 1j * np.pi * ALPHAP / 2)]
    return [(H + 0j, lam + 0j)]

# ---------------------------------------------------------------- data
def load(csv):
    d = np.loadtxt(csv, delimiter=',', skiprows=1)
    return d

def bin_avg_matrix(lo, hi, n=8):
    x, w = leggauss(n)
    T = 0.5 * (hi - lo)[:, None] * x[None, :] + 0.5 * (hi + lo)[:, None]
    W = 0.5 * w[None, :] * np.ones_like(T)
    return T, W

_KCACHE = {}
def binned(comps, T, W):
    key = id(T)
    if key not in _KCACHE:
        _KCACHE[key] = 4 * np.pi * j0(np.sqrt(T.ravel())[:, None] * bb[None, :]) * (bb * wbb)[None, :]
    G = 1 - np.exp(-omega(bb, comps) / 2)
    a = _KCACHE[key] @ G
    v = (np.abs(a)**2 / (16 * np.pi) * GEV2_TO_MB).reshape(T.shape)
    return (v * W).sum(axis=1)

def dip_bump(comps):
    """first local minimum of dsigma/dt beyond 0.1 GeV^2 and the following local maximum"""
    t = np.linspace(0.1, 1.5, 7001)
    v = dsdt(t, comps)
    i = np.where((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:]))[0]
    if len(i) == 0:
        return None
    k = i[0] + 1
    j = np.where((v[k+1:-1] > v[k:-2]) & (v[k+1:-1] > v[k+2:]))[0]
    kb = k + 1 + j[0]
    return float(t[k]), float(v[k]), float(t[kb]), float(v[kb]), float(v[kb] / v[k])

def describe(comps, d, T, W):
    lo, hi, tr, y, st, sy = d.T
    m = binned(comps, T, W)
    err = np.hypot(st, sy)
    r = m / y
    sl = (tr >= 0.04) & (tr <= 0.2)
    D = float(np.mean(((np.log(m) - np.log(y)) / (err / y))**2))
    D1 = float(np.mean((((np.log(m) - np.log(y)) / (err / y))**2)[tr < 1.0]))
    return dict(frac10=float(np.mean(abs(r - 1) < 0.1)), frac10_below1=float(np.mean(abs(r[tr < 1] - 1) < 0.1)),
                slope=float(-np.polyfit(tr[sl], np.log(m[sl]), 1)[0]),
                D_all=D, D_below1=D1, ratio_min_above1=float(r[tr > 1].min()), ratio_max_above1=float(r[tr > 1].max())), m

# ---------------------------------------------------------------- main
def main(csv, outdir):
    os.makedirs(outdir, exist_ok=True)
    d = load(csv); lo, hi, tr, y, st, sy = d.T
    T, W = bin_avg_matrix(lo, hi)
    R = {}
    # check against closed form (single real Gaussian, set B)
    lB, nB = calibrate_B(TOTEM13['stot'], TOTEM13['r'])
    HB = 4 * np.pi * lB * nB
    cB = comps_from(HB, lB, phase=False)
    obsB = observables(cB)
    tt = np.array([0.0, 0.1, 0.3, 0.6, 1.5])
    series = np.array([8 * np.pi * lB * F(lB * t, nB / 2) for t in tt])
    R['check_numeric_vs_series_maxrel'] = float(np.max(np.abs(amp(tt, cB).real - series) / np.abs(series)))
    R['check_stot_B'] = obsB['stot']
    # ---- (1a) set A with phase (same H, lambda as set A at 13 TeV)
    HA, lA, nA = params_A(13000.0)
    cAc = comps_from(HA, lA)
    R['Ac'] = observables(cAc); R['Ac']['dipbump'] = dip_bump(cAc)
    R['Ac']['desc'], mAc = describe(cAc, d, T, W)
    # ---- (1b) set Bc: phase from Delta, alpha'; (H, lambda) calibrated to sigma_tot and r
    def resid_B(p):
        o = observables(comps_from(p[0], p[1]))
        return [o['stot'] / TOTEM13['stot'] - 1, o['r'] / TOTEM13['r'] - 1]
    sol = least_squares(resid_B, [HB, lB], xtol=1e-14, ftol=1e-14)
    HBc, lBc = sol.x
    cBc = comps_from(HBc, lBc)
    R['Bc'] = observables(cBc); R['Bc'].update(H=HBc, lam=lBc, nu=HBc / (4 * np.pi * lBc)); R['Bc']['dipbump'] = dip_bump(cBc)
    R['Bc']['desc'], mBc = describe(cBc, d, T, W)
    R['B'] = observables(cB); R['B']['dipbump_real_zero'] = True
    R['B']['desc'], mB = describe(cB, d, T, W)
    # ---- (2) two-component opacity with phase, fitted to the 290 bins (log-space, stat+syst in quadrature)
    err = np.hypot(st, sy)
    # set C: two components with the signature phase.  Four parameters: total coupling H = H1 + H2,
    # width lambda1, fraction f = H2/H, width ratio k = lambda2/lambda1.  (H, lambda1) are fixed by
    # sigma_tot and r exactly as for set B; (f, k) by an equal-weight least-squares fit of ln(dsigma/dt)
    # to the bins with |t| < TFIT (descriptive; correlations not used).
    TFIT = float(os.environ.get("TFIT", "1.0"))
    msk = tr < TFIT
    def comps_C(H, l1, f, k):
        return comps_from(H * (1 - f), l1) + comps_from(H * f, l1 * k)
    cache = {'x': np.array([HBc, lBc])}
    def calib(f, k):
        def r2(p):
            o = observables(comps_C(p[0], p[1], f, k))
            return [o['stot'] / TOTEM13['stot'] - 1, o['r'] / TOTEM13['r'] - 1]
        sol = least_squares(r2, cache['x'], xtol=1e-13, ftol=1e-13, bounds=([1, 0.5], [5000, 50]))
        cache['x'] = sol.x
        return sol.x
    def resid_C(q):
        f, k = 1 / (1 + np.exp(-q[0])), np.exp(q[1])
        H, l1 = calib(f, k)
        m = binned(comps_C(H, l1, f, k), T, W)
        return (np.log(np.maximum(m, 1e-300)) - np.log(y))[msk]
    best = None
    for q0 in [(-1.5, -1.0), (-2.5, -1.5), (-1.0, -0.5), (-3.0, -2.0), (0.0, -0.3), (-0.5, 0.5), (-4.0, -2.5)]:
        cache['x'] = np.array([HBc, lBc])
        sres = least_squares(resid_C, np.array(q0, float), diff_step=1e-4, max_nfev=200)
        if best is None or sres.cost < best.cost:
            best = sres
    fC, kC = 1 / (1 + np.exp(-best.x[0])), np.exp(best.x[1])
    cache['x'] = np.array([HBc, lBc]); HC, l1C = calib(fC, kC)
    cC = comps_C(HC, l1C, fC, kC)
    R['C'] = observables(cC)
    R['C'].update(H=HC, lam1=l1C, f=fC, k=kC, lam2=l1C * kC, H1=HC * (1 - fC), H2=HC * fC,
                  nu1=HC * (1 - fC) / (4 * np.pi * l1C), nu2=HC * fC / (4 * np.pi * l1C * kC), n_params=4, t_fit=TFIT,
                  rms_log_dev_fit=float(np.sqrt(2 * best.cost / msk.sum())))
    R['C']['dipbump'] = dip_bump(cC)
    R['C']['desc'], mC = describe(cC, d, T, W)
    # derivative dispersion relation check for set A family: rho ~ (pi/2) dln sigma/dln s (PDG Eq. 20.19)
    s13 = 13000.0
    st_p = observables(comps_from(*params_A(s13 * 1.01)[:2]))['stot']
    st_m = observables(comps_from(*params_A(s13 / 1.01)[:2]))['stot']
    R['Ac']['rho_dispersion'] = float(np.pi / 2 * (np.log(st_p) - np.log(st_m)) / (4 * np.log(1.01)))
    # outputs
    np.savetxt(os.path.join(outdir, 'extensions_bins.csv'),
               np.column_stack([lo, hi, tr, y, st, sy, mAc, mB, mBc, mC]), delimiter=',', fmt='%.10g',
               header='t_low,t_high,t_repr,TOTEM,stat,syst,Ac_binavg,B_binavg,Bc_binavg,C_binavg', comments='')
    tg = np.linspace(0.001, 4.1, 3000)
    np.savetxt(os.path.join(outdir, 'extensions_curves.csv'),
               np.column_stack([tg, dsdt(tg, cB), dsdt(tg, cBc), dsdt(tg, cC), dsdt(tg, cAc)]), delimiter=',', fmt='%.10g',
               header='t,B,Bc,C,Ac', comments='')
    json.dump(R, open(os.path.join(outdir, 'extensions.json'), 'w'), indent=1, default=float)
    return R

if __name__ == '__main__':
    R = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'output')
    print(json.dumps(R, indent=1, default=float))
