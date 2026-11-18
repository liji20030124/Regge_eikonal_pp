#!/usr/bin/env python3
"""
Single-channel Gaussian Regge-eikonal model for high-energy pp elastic scattering.

Conventions and equation numbers of the paper:
  sigma_tot = Im A(s,0)/s,  dsigma/dt = |A|^2/(16 pi s^2)                      Eq. (1)
  A_B = i s H e^{lambda t},  H = H_ref e^{Delta L},  lambda = lambda_ref + alpha' L,
  L = ln(s/s_ref)                                                               Eq. (7)
  chi = i Omega/2,  Omega = nu exp(-b^2/(4 lambda)),  nu = H/(4 pi lambda)      Eq. (11)
  Gamma = 1 - exp(-Omega/2),  G_inel = 1 - exp(-Omega)                          Eqs. (12), (13)
  A = 8 pi i s lambda F(Q,a),  F = sum (-1)^{n+1} a^n e^{-Q/n}/(n! n),
  a = nu/2,  Q = lambda |t|                                                     Eq. (16)
  sigma_tot = 8 pi lambda Ein(a),  sigma_inel = 4 pi lambda Ein(2a)              Eqs. (17), (18)
  dsigma/dt = 4 pi lambda^2 F^2,  B(0) = 2 lambda I1(a)/Ein(a)                  Eq. (19)
  r(a) = 1 - Ein(2a)/(2 Ein(a))                                                 Eq. (21)

Usage:  python3 regge_eikonal_pp.py [output_dir]
Writes results.json, energy_scan.csv, differential_13TeV.csv, dsdt_curves_13TeV.csv.
Only standard numpy/scipy/mpmath are used.  No experimental number is fitted except
in parameter set B, which is a calibration to two TOTEM central values by construction.
"""
import json, os, sys
import numpy as np
import mpmath as mp
from scipy.special import gammaln, j0, exp1
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

GEV2_TO_MB = 0.3893793721      # 1 GeV^-2 in mb
EULER_GAMMA = 0.5772156649015329

# ---------------------------------------------------------------- inputs
SET_A = dict(sqrts_ref=7000.0, lam_ref=10.0, nu_ref=4.0, Delta=0.0808, alphaP=0.25)
TOTEM13 = dict(stot=110.6, stot_err=3.4, sel=31.0, sel_err=1.7, sinel=79.5, sinel_err=1.8,
               r=0.281, tdip=0.47, tdip_stat=0.004, tdip_syst=0.01)
# six representative bins of TOTEM 13 TeV dsigma/dt [mb/GeV^2] (Ref. [14])
TOTEM13_POINTS = [(0.04990, 227.243), (0.10012, 81.31), (0.19979, 10.595),
                  (0.29960, 1.18712), (0.46799, 0.0272682), (0.63179, 0.0483139)]

# ---------------------------------------------------------------- special functions
def Ein(a):
    """Entire exponential integral, series form, Eq. (17)."""
    return float(mp.nsum(lambda n: (-1)**(n + 1) * mp.mpf(a)**n / (mp.factorial(n) * n), [1, mp.inf]))

def Ein_closed(a):
    """Independent form Ein(a) = gamma_E + ln a + E1(a)."""
    return EULER_GAMMA + np.log(a) + exp1(a)

def I1(a):
    return float(mp.nsum(lambda n: (-1)**(n + 1) * mp.mpf(a)**n / (mp.factorial(n) * n**2), [1, mp.inf]))

def F(Q, a, nmax=250):
    n = np.arange(1, nmax + 1)
    terms = (-1.0)**(n + 1) * np.exp(n * np.log(a) - gammaln(n + 1)) / n * np.exp(-Q / n)
    return float(np.sum(terms))

# ---------------------------------------------------------------- model
def params_A(sqrts, p=SET_A):
    L = 2.0 * np.log(sqrts / p['sqrts_ref'])
    H_ref = 4.0 * np.pi * p['lam_ref'] * p['nu_ref']
    H = H_ref * np.exp(p['Delta'] * L)
    lam = p['lam_ref'] + p['alphaP'] * L
    return H, lam, H / (4.0 * np.pi * lam)

def dsdt(t, lam, a):
    return 4.0 * np.pi * lam**2 * F(lam * t, a)**2 * GEV2_TO_MB

def dsdt_born(t, H, lam):
    return H**2 / (16.0 * np.pi) * np.exp(-2.0 * lam * t) * GEV2_TO_MB

def first_zero(lam, a):
    ts = np.linspace(1e-3, 3.0, 30001)
    f = np.array([F(lam * t, a) for t in ts])
    i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    return brentq(lambda t: F(lam * t, a), ts[i], ts[i + 1], xtol=1e-14)

def second_max(lam, a, t0):
    ts = np.linspace(t0 + 1e-4, 3.0, 3001)
    v = np.array([F(lam * t, a)**2 for t in ts])
    k = int(np.argmax(v))
    r = minimize_scalar(lambda t: -F(lam * t, a)**2, bounds=(ts[max(k - 1, 0)], ts[k + 1]),
                        method='bounded', options={'xatol': 1e-12})
    return float(r.x), dsdt(r.x, lam, a)

def observables(lam, nu, H=None):
    a = nu / 2.0
    i0, i02 = Ein(a), Ein(2 * a)
    st = 8 * np.pi * lam * i0
    si = 4 * np.pi * lam * i02
    t0 = first_zero(lam, a)
    tb, vb = second_max(lam, a, t0)
    out = dict(lam=lam, nu=nu, a=a,
               stot=st * GEV2_TO_MB, sinel=si * GEV2_TO_MB, sel=(st - si) * GEV2_TO_MB,
               ratio_el_tot=(st - si) / st, B0=2 * lam * I1(a) / i0,
               t_zero=t0, t_secondmax=tb, dsdt_secondmax=vb,
               Gamma0=1 - np.exp(-a), Ginel0=1 - np.exp(-nu))
    if H is not None:
        out.update(H=H, stot_born=H * GEV2_TO_MB, B0_born=2 * lam,
                   ratio_born=nu / 8.0, sel_born=H**2 / (32 * np.pi * lam) * GEV2_TO_MB)
    return out

def calibrate_B(stot_mb, r):
    ratio = lambda a: 1 - Ein(2 * a) / (2 * Ein(a))
    a = brentq(lambda x: ratio(x) - r, 1e-3, 50.0, xtol=1e-14)
    lam = (stot_mb / GEV2_TO_MB) / (8 * np.pi * Ein(a))
    return lam, 2 * a

# ---------------------------------------------------------------- numerical checks
def checks(lam, nu, H):
    a = nu / 2
    Om = lambda b: nu * np.exp(-b * b / (4 * lam))
    G = lambda b: 1 - np.exp(-Om(b) / 2)
    bmax = 40 * np.sqrt(lam)
    st_b = 2 * quad(lambda b: 2 * np.pi * b * G(b), 0, bmax, epsabs=1e-13, epsrel=1e-13, limit=500)[0]
    se_b = quad(lambda b: 2 * np.pi * b * G(b)**2, 0, bmax, epsabs=1e-13, epsrel=1e-13, limit=500)[0]
    si_b = quad(lambda b: 2 * np.pi * b * (1 - np.exp(-Om(b))), 0, bmax, epsabs=1e-13, epsrel=1e-13, limit=500)[0]
    st, si = 8 * np.pi * lam * Ein(a), 4 * np.pi * lam * Ein(2 * a)
    se_t = quad(lambda t: 4 * np.pi * lam**2 * F(lam * t, a)**2, 0, 60 / lam, limit=500, epsrel=1e-12)[0]
    hank = []
    for t in [0.0, 0.1, 0.3, 0.6]:
        q = np.sqrt(t)
        num = 2 * quad(lambda b: 2 * np.pi * b * j0(q * b) * G(b), 0, bmax, limit=800, epsabs=1e-13)[0]
        ser = 8 * np.pi * lam * F(lam * t, a)
        hank.append(dict(t=t, series=ser, hankel=num, rel_diff=abs(num - ser) / abs(ser)))
    # local unitarity 2 Re Gamma = |Gamma|^2 + G_inel on a grid
    bs = np.linspace(0, bmax, 2001)
    unit = float(np.max(np.abs(2 * G(bs) - G(bs)**2 - (1 - np.exp(-Om(bs))))))
    # Born limit
    nu_s, lam_s = 1e-4, 10.0
    born_ratio = 8 * np.pi * lam_s * Ein(nu_s / 2) / (4 * np.pi * lam_s * nu_s)
    # forward slope by numerical derivative
    h = 1e-6
    B_num = (np.log(F(0, a)**2) - np.log(F(lam * h, a)**2)) / h
    return dict(
        stot_rel=abs(st_b - st) / st, sel_rel=abs(se_b - (st - si)) / (st - si),
        sinel_rel=abs(si_b - si) / si, sel_tint_rel=abs(se_t - (st - si)) / (st - si),
        hankel=hank, local_unitarity_maxdev=unit,
        born_limit_ratio=born_ratio, born_limit_expected=1 - nu_s / 8 + nu_s**2 / 72,
        Ein_series_vs_closed=abs(Ein(a) - Ein_closed(a)) / Ein(a),
        B0_series=2 * lam * I1(a) / Ein(a), B0_numderiv=B_num)

# ---------------------------------------------------------------- main
def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    R = {}
    # set A at three energies
    for rs in (7000.0, 8000.0, 13000.0):
        H, lam, nu = params_A(rs)
        R[f'A_{int(rs)}GeV'] = observables(lam, nu, H)
    HA, lA, nA = params_A(13000.0)
    R['A_13TeV_vs_TOTEM'] = dict(
        dstot_mb=R['A_13000GeV']['stot'] - TOTEM13['stot'],
        dstot_rel=(R['A_13000GeV']['stot'] - TOTEM13['stot']) / TOTEM13['stot'],
        dsinel_mb=R['A_13000GeV']['sinel'] - TOTEM13['sinel'],
        dsel_mb=R['A_13000GeV']['sel'] - TOTEM13['sel'],
        born_to_eik_discrepancy_factor=(R['A_13000GeV']['stot_born'] - TOTEM13['stot']) /
                                       (R['A_13000GeV']['stot'] - TOTEM13['stot']),
        rho_born=np.tan(np.pi * SET_A['Delta'] / 2))
    # set B calibration
    lB, nB = calibrate_B(TOTEM13['stot'], TOTEM13['r'])
    R['B_13TeV'] = observables(lB, nB)
    # differential comparison
    rows = []
    for t, e in TOTEM13_POINTS:
        mA, mB = dsdt(t, lA, nA / 2), dsdt(t, lB, nB / 2)
        rows.append(dict(t=t, totem=e, setA=mA, A_over_TOTEM=mA / e, setB=mB, B_over_TOTEM=mB / e))
    R['differential_13TeV'] = rows
    # log-linear slope of model in 0.04-0.2 GeV^2 (uniform grid, diagnostic only)
    ts = np.linspace(0.04, 0.2, 401)
    for key, (l, n) in {'A': (lA, nA), 'B': (lB, nB)}.items():
        R[f'loglin_slope_0.04_0.2_{key}'] = -np.polyfit(ts, np.log([dsdt(t, l, n / 2) for t in ts]), 1)[0]
    # numerical checks at 13 TeV, set A
    R['checks_A_13TeV'] = checks(lA, nA, HA)
    # energy scan
    scan = []
    for rs in np.geomspace(1000.0, 1.0e5, 41):
        H, l, n = params_A(rs)
        a = n / 2
        st, si = 8 * np.pi * l * Ein(a), 4 * np.pi * l * Ein(2 * a)
        scan.append([rs / 1000, H * GEV2_TO_MB, st * GEV2_TO_MB, si * GEV2_TO_MB, (st - si) * GEV2_TO_MB,
                     (st - si) / st, n / 8, 2 * l, 2 * l * I1(a) / Ein(a)])
    np.savetxt(os.path.join(outdir, 'energy_scan.csv'), scan, delimiter=',', fmt='%.10g',
               header='sqrts_TeV,stot_Born_mb,stot_mb,sinel_mb,sel_mb,ratio_eik,ratio_Born,B0_Born,B0_eik', comments='')
    with open(os.path.join(outdir, 'differential_13TeV.csv'), 'w') as f:
        f.write('t_GeV2,TOTEM,setA,A_over_TOTEM,setB,B_over_TOTEM\n')
        for r in rows:
            f.write(','.join(f'{r[k]:.10g}' for k in ['t', 'totem', 'setA', 'A_over_TOTEM', 'setB', 'B_over_TOTEM']) + '\n')
    tg = np.linspace(0.001, 1.0, 1000)
    np.savetxt(os.path.join(outdir, 'dsdt_curves_13TeV.csv'),
               [[t, dsdt_born(t, HA, lA), dsdt(t, lA, nA / 2), dsdt(t, lB, nB / 2)] for t in tg],
               delimiter=',', fmt='%.10g', header='t_GeV2,Born_A,eik_A,eik_B', comments='')
    with open(os.path.join(outdir, 'results.json'), 'w') as f:
        json.dump(R, f, indent=1, default=float)
    return R

if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'output'
    R = main(out)
    A, B = R['A_13000GeV'], R['B_13TeV']
    print('Set A, 13 TeV: lam=%.6f nu=%.6f a=%.6f' % (A['lam'], A['nu'], A['a']))
    print('  stot=%.4f sel=%.4f sinel=%.4f mb  Born stot=%.4f mb' % (A['stot'], A['sel'], A['sinel'], A['stot_born']))
    print('  B0=%.4f  B0_Born=%.4f  t0=%.6f  2nd max at %.4f  ratio=%.4f  Born ratio=%.4f' %
          (A['B0'], A['B0_born'], A['t_zero'], A['t_secondmax'], A['ratio_el_tot'], A['ratio_born']))
    print('Set B, 13 TeV: lam=%.8f nu=%.8f a=%.8f' % (B['lam'], B['nu'], B['a']))
    print('  stot=%.4f sel=%.4f sinel=%.4f  B0=%.4f  t0=%.6f  2nd max at %.4f' %
          (B['stot'], B['sel'], B['sinel'], B['B0'], B['t_zero'], B['t_secondmax']))
    for r in R['differential_13TeV']:
        print('  t=%.5f TOTEM=%-10g A=%-10.5g (%.4f)  B=%-10.5g (%.4f)' %
              (r['t'], r['totem'], r['setA'], r['A_over_TOTEM'], r['setB'], r['B_over_TOTEM']))
    c = R['checks_A_13TeV']
    print('Checks: stot %.1e sel %.1e sinel %.1e t-int %.1e unit %.1e Ein %.1e Born %.10f (exp %.10f)' %
          (c['stot_rel'], c['sel_rel'], c['sinel_rel'], c['sel_tint_rel'], c['local_unitarity_maxdev'],
           c['Ein_series_vs_closed'], c['born_limit_ratio'], c['born_limit_expected']))
    for h in c['hankel']:
        print('  Hankel t=%.1f rel diff %.1e' % (h['t'], h['rel_diff']))
