#!/usr/bin/env python3
"""
Compare program output (output/results.json) with
  (i)  every number printed in the manuscript (Regge_eikonal_pp.tex), and
  (ii) the numbers in the handoff summary (强子 Regge-eikonal 研究汇总, 2026-11-18).
For (i) the printed string must occur in the .tex file, and the program value rounded
to the printed number of decimals must reproduce it.  For (ii) the relative difference
is reported.  Usage: python3 compare.py output/results.json Regge_eikonal_pp.tex
"""
import json, sys, re
R = json.load(open(sys.argv[1]))
FC = json.load(open(sys.argv[1].replace('results.json', 'full_comparison.json')))
EX = json.load(open(sys.argv[1].replace('results.json', 'extensions.json')))
EX4 = json.load(open(sys.argv[1].replace('output/results.json', 'output_t4/extensions.json')))
tex = open(sys.argv[2]).read()
A, B, A7 = R['A_13000GeV'], R['B_13TeV'], R['A_7000GeV']
V, C = R['A_13TeV_vs_TOTEM'], R['checks_A_13TeV']
# Table 2 of the paper uses bin-averaged model values (full_comparison.py)
D = [dict(t=s['t'], setA=s['A_bin'], A_over_TOTEM=s['A_ratio'], setB=s['B_bin'], B_over_TOTEM=s['B_ratio']) for s in FC['six']]

# (quantity, program value, string printed in paper, location in paper)
paper = [
 ('A: lambda [GeV^-2]',            A['lam'],        '10.3095',  'Table 1'),
 ('A: nu',                          A['nu'],         '4.28812',  'Table 1'),
 ('A: Gamma(0)',                    A['Gamma0'],     '0.8828',   'Table 1'),
 ('A: H = sigma_tot^Born [mb]',     A['stot_born'],  '216.32',   'Table 1'),
 ('A: sigma_tot [mb]',              A['stot'],       '139.23',   'Table 1'),
 ('A: sigma_tot [mb]',              A['stot'],       '139.2',    'Abstract, Sec. VI'),
 ('A: sigma_el [mb]',               A['sel'],        '36.54',    'Table 1'),
 ('A: sigma_inel [mb]',             A['sinel'],      '102.69',   'Table 1'),
 ('A: B(0) [GeV^-2]',               A['B0'],         '25.53',    'Table 1'),
 ('A: |t_0| [GeV^2]',               A['t_zero'],     '0.2442',   'Table 1'),
 ('A: |t_0| [GeV^2]',               A['t_zero'],     '0.244',    'Abstract, Sec. VI'),
 ('A: 2nd max [GeV^2]',             A['t_secondmax'],'0.357',    'Sec. VI'),
 ('A: sigma_tot - TOTEM [mb]',      V['dstot_mb'],   '28.6',     'Sec. VI'),
 ('A: relative excess [%]',         100*V['dstot_rel'],'25.9',   'Sec. VI'),
 ('A: sigma_inel - TOTEM [mb]',     V['dsinel_mb'],  '23.2',     'Sec. VI'),
 ('A: sigma_el - TOTEM [mb]',       V['dsel_mb'],    '5.5',      'Sec. VI'),
 ('A: (H-exp)/(eik-exp)',           V['born_to_eik_discrepancy_factor'], '3.7', 'Sec. VI'),
 ('A: Born ratio nu/8',             A['ratio_born'], '0.536',    'Sec. VI, VIII'),
 ('A: eikonal sigma_el/sigma_tot',  A['ratio_el_tot'],'0.262',   'Sec. VI'),
 ('A: Gamma(0)',                    A['Gamma0'],     '0.883',    'Sec. VI'),
 ('rho_Born = tan(pi Delta/2)',     V['rho_born'],   '0.13',     'Sec. III, VII'),
 ('B: lambda [GeV^-2]',             B['lam'],        '7.38434',  'Table 1'),
 ('B: nu',                          B['nu'],         '5.06608',  'Table 1'),
 ('B: Gamma(0)',                    B['Gamma0'],     '0.9206',   'Table 1'),
 ('B: Gamma(0)',                    B['Gamma0'],     '0.921',    'Sec. VI'),
 ('B: sigma_tot [mb] (input)',      B['stot'],       '110.60',   'Table 1'),
 ('B: sigma_el [mb] (input)',       B['sel'],        '31.08',    'Table 1'),
 ('B: sigma_inel [mb] (input)',     B['sinel'],      '79.52',    'Table 1'),
 ('B: B(0) [GeV^-2]',               B['B0'],         '18.83',    'Table 1, Sec. VI'),
 ('B: |t_0| [GeV^2]',               B['t_zero'],     '0.3155',   'Table 1'),
 ('B: |t_0| [GeV^2]',               B['t_zero'],     '0.315',    'Abstract, Sec. VI-VIII'),
 ('B: 2nd max [GeV^2]',             B['t_secondmax'],'0.470',    'Sec. VI'),
 ('B/TOTEM at 0.468',               D[4]['B_over_TOTEM'], '14.8', 'Sec. VI'),
 ('B/TOTEM at 0.632',               D[5]['B_over_TOTEM'], '4.0',  'Sec. VI'),
 ('B/TOTEM at 0.2: deficit [%]',    100*(1-D[2]['B_over_TOTEM']), '35', 'Sec. VI'),
 ('A: min G_inel(0)',               min(A['Ginel0'], B['Ginel0']), '0.98', 'Sec. VI (lower bound)'),
 ('Born-limit sigma_tot/H (nu=1e-4)', C['born_limit_ratio'], '0.9999875', 'Sec. V'),
]
for i, (t, d) in enumerate([(0.04990, 0), (0.10012, 1), (0.19979, 2), (0.29960, 3), (0.46799, 4), (0.63179, 5)]):
    pass
tab2 = [('263.95', 'setA'), ('60.93', 'setA'), ('1.033', 'setA'), ('0.2862', 'setA'), ('0.2279', 'setA'), ('0.02594', 'setA'),
        ('1.162', 'A_over_TOTEM'), ('0.749', 'A_over_TOTEM'), ('0.098', 'A_over_TOTEM'), ('0.241', 'A_over_TOTEM'), ('8.36', 'A_over_TOTEM'), ('0.537', 'A_over_TOTEM'),
        ('237.35', 'setB'), ('83.39', 'setB'), ('6.931', 'setB'), ('0.03720', 'setB'), ('0.4025', 'setB'), ('0.1945', 'setB'),
        ('1.044', 'B_over_TOTEM'), ('1.026', 'B_over_TOTEM'), ('0.654', 'B_over_TOTEM'), ('0.031', 'B_over_TOTEM'), ('14.76', 'B_over_TOTEM'), ('4.03', 'B_over_TOTEM')]
paper += [
 ('App.B: data slope 0.04-0.2', FC['slope_data'], '20.39', 'App. B'),
 ('App.B: set B slope (bin avg)', FC['slope_B_binavg'], '23.06', 'App. B'),
 ('App.B: set A slope (bin avg)', FC['slope_A_binavg'], '34.98', 'App. B'),
 ('App.B: n bins in slope fit', FC['n_slope_points'], '94', 'App. B'),
 ('App.B: data bump/dip ratio', FC['data_bump_dip_ratio'], '1.772', 'App. B'),
 ('App.B: B within 10% [%]', 100*FC['frac_within_10pct_B'], '24', 'App. B'),
 ('App.B: A within 10% [%]', 100*FC['frac_within_10pct_A'], '7', 'App. B'),
 ('App.B: B within 5% up to |t|', FC['max_t_within_5pct_contiguous_B'], '0.068', 'App. B'),
 ('App.A: max point-vs-bin, all bins [%]', 100*FC['max_rel_point_vs_binavg_B'], '98', 'App. A'),
]
import numpy as _np
def _binratio(o, col):
    b = _np.loadtxt(sys.argv[1].replace('results.json', '').replace('output/', o + '/') + 'extensions_bins.csv', delimiter=',', skiprows=1)
    k = int(_np.argmin(abs(b[:, 2] - 0.46799))); return b[k, col] / b[k, 3]
def _rms(o, col, below1):
    b = _np.loadtxt(sys.argv[1].replace('results.json', '').replace('output/', o + '/') + 'extensions_bins.csv', delimiter=',', skiprows=1)
    m = b[:, 2] < 1 if below1 else _np.ones(len(b), bool); return float(_np.sqrt(_np.mean(_np.log(b[m, col] / b[m, 3])**2)))
paper += [
 ('Sec.VII: Ac rho', EX['Ac']['rho'], '0.0950', 'Sec. VII'),
 ('Sec.VII: Ac rho (dispersion)', EX['Ac']['rho_dispersion'], '0.0950', 'Sec. VII'),
 ('Tab.ext: Bc lambda', EX['Bc']['lam'], '7.516', 'Table 6'), ('Tab.ext: Bc nu', EX['Bc']['nu'], '4.909', 'Table 6'),
 ('Tab.ext: Bc rho', EX['Bc']['rho'], '0.098', 'Table 6'), ('Tab.ext: Bc B0', EX['Bc']['B0'], '19.06', 'Table 6'),
 ('Tab.ext: Bc slope', EX['Bc']['desc']['slope'], '23.21', 'Table 6'),
 ('Tab.ext: Bc min', EX['Bc']['dipbump'][0], '0.317', 'Table 6'), ('Tab.ext: Bc max', EX['Bc']['dipbump'][2], '0.460', 'Table 6'),
 ('Tab.ext: Bc R', EX['Bc']['dipbump'][4], '4.67', 'Table 6'),
 ('Tab.ext: C lam1', EX['C']['lam1'], '8.947', 'Table 6'), ('Tab.ext: C lam2', EX['C']['lam2'], '1.983', 'Table 6'),
 ('Tab.ext: C nu1', EX['C']['nu1'], '3.248', 'Table 6'), ('Tab.ext: C nu2', EX['C']['nu2'], '6.353', 'Table 6'),
 ('Tab.ext: C rho', EX['C']['rho'], '0.103', 'Table 6'), ('Tab.ext: C B0', EX['C']['B0'], '19.88', 'Table 6'),
 ('Tab.ext: C slope', EX['C']['desc']['slope'], '21.14', 'Table 6'),
 ('Tab.ext: C frac10 <1 [%]', 100*EX['C']['desc']['frac10_below1'], '47', 'Table 6'),
 ('Tab.ext: C\' lam1', EX4['C']['lam1'], '8.784', 'Table 6'), ('Tab.ext: C\' lam2', EX4['C']['lam2'], '2.849', 'Table 6'),
 ('Tab.ext: C\' nu1', EX4['C']['nu1'], '3.279', 'Table 6'), ('Tab.ext: C\' nu2', EX4['C']['nu2'], '3.623', 'Table 6'),
 ('Tab.ext: C\' rho', EX4['C']['rho'], '0.100', 'Table 6'), ('Tab.ext: C\' B0', EX4['C']['B0'], '19.67', 'Table 6'),
 ('Tab.ext: C\' slope', EX4['C']['desc']['slope'], '21.61', 'Table 6'),
 ('Tab.ext: C\' min', EX4['C']['dipbump'][0], '0.397', 'Table 6'), ('Tab.ext: C\' max', EX4['C']['dipbump'][2], '0.531', 'Table 6'),
 ('Tab.ext: C\' R', EX4['C']['dipbump'][4], '1.56', 'Table 6'),
 ('Tab.ext: C\' frac10 <1 [%]', 100*EX4['C']['desc']['frac10_below1'], '31', 'Table 6'),
 ('Tab.ext: Bc frac10 <1 [%]', 100*EX['Bc']['desc']['frac10_below1'], '25', 'Table 6'),
 ('Tab.ext: B frac10 <1 [%]', 100*EX['B']['desc']['frac10_below1'], '26', 'Table 6'),
 ('Tab.ext: B ratio dip', _binratio('output', 7), '14.76', 'Table 6'), ('Tab.ext: Bc ratio dip', _binratio('output', 8), '14.97', 'Table 6'),
 ('Tab.ext: C ratio dip', _binratio('output', 9), '2.04', 'Table 6'), ('Tab.ext: C\' ratio dip', _binratio('output_t4', 9), '3.69', 'Table 6'),
 ('Tab.ext: B rms<1', _rms('output', 7, True), '1.80', 'Table 6'), ('Tab.ext: Bc rms<1', _rms('output', 8, True), '1.20', 'Table 6'),
 ('Tab.ext: C rms<1', _rms('output', 9, True), '0.21', 'Table 6'), ('Tab.ext: C\' rms<1', _rms('output_t4', 9, True), '0.54', 'Table 6'),
 ('Tab.ext: B rms all', _rms('output', 7, False), '1.86', 'Table 6'), ('Tab.ext: Bc rms all', _rms('output', 8, False), '1.25', 'Table 6'),
 ('Tab.ext: C rms all', _rms('output', 9, False), '1.00', 'Table 6'), ('Tab.ext: C\' rms all', _rms('output_t4', 9, False), '0.71', 'Table 6'),
 ('Sec.VII: C overshoot min', EX['C']['desc']['ratio_min_above1'], '1.9', 'Sec. VII'),
 ('Sec.VII: C overshoot max', EX['C']['desc']['ratio_max_above1'], '739', 'Sec. VII'),
 ('Sec.VII: check numeric vs series', EX['check_numeric_vs_series_maxrel']*1e12, '5', 'Sec. VII (x1e-12)'),
]
GW = json.load(open(sys.argv[1].replace('results.json', 'gw_gamma_scan.json')))
GWfit = json.load(open(sys.argv[1].replace('results.json', 'goodwalker.json')))
for g, lam, nu, smd, tmin, rd in [(0.2,'7.35','5.25','0.71','0.311','16.9'), (0.4,'6.87','6.58','2.71','0.298','22.6'), (0.6,'6.06','10.46','5.82','0.283','29.7')]:
    e = [x for x in GW if abs(x['g'] - g) < 1e-9][0]
    paper += [(f'GW g={g} lambda', e['lam'], lam, 'Sec. VII'), (f'GW g={g} nu', e['nu'], nu, 'Sec. VII'),
              (f'GW g={g} sigma_lowmass', e['smd'], smd, 'Sec. VII'), (f'GW g={g} t_min', e['tmin'], tmin, 'Sec. VII'),
              (f'GW g={g} ratio dip', e['ratio_dip'], rd, 'Sec. VII')]
paper += [('GW Pumplin max ratio', max(x['pump'] for x in GW), '0.92', 'Sec. VII') ,
          ('GW best-fit gamma', GWfit['gamma'], '0', 'Sec. VII'),
          ('asymptotic coefficient 8 pi alpha Delta [mb]', 8*3.141592653589793*0.25*0.0808*0.3893793721, '0.198', 'Sec. VII')]
# Secs. VI-VIII: profile reconstruction (Eqs. (22)-(23)), Gamma(0) energy dependence, models M1-M4 (Table 4),
# one-parameter curve W(a) (Eq. (29)); produced by profile_inversion.py, energy_blackness.py, kmr_fits.py,
# universal_curves.py
import os as _os
_od = _os.path.dirname(sys.argv[1])
if _os.path.exists(_os.path.join(_od, 'kmr_fits.json')):
    PI = json.load(open(_os.path.join(_od, 'profile_inversion.json')))
    EB = json.load(open(_os.path.join(_od, 'energy_blackness.json')))
    KF = json.load(open(_os.path.join(_od, 'kmr_fits.json')))
    UC = json.load(open(_os.path.join(_od, 'universal_curves.json')))
    pb, pv = PI['base'], PI['variations'].values()
    paper += [('profile: Re Gamma(0)', pb['ReGamma0'], '1.020', 'Eq. (23)'),
              ('profile: G_inel(0)', pb['Ginel0'], '0.989', 'Eq. (23)'),
              ('profile: max G_inel', pb['Ginel_max'], '0.990', 'Eq. (23)'),
              ('profile: b at max [fm]', pb['b_at_max_fm'], '0.33', 'Eq. (23)'),
              ('profile: sigma_tot [mb]', PI['check']['stot_mb'], '110.4', 'Sec. VI'),
              ('profile: sigma_inel [mb]', PI['check']['sinel_mb'], '79.5', 'Sec. VI'),
              ('profile: min Re Gamma(0) over variations', min(v['ReGamma0'] for v in pv), '1.002', 'Sec. VI'),
              ('profile: max Re Gamma(0) over variations', max(v['ReGamma0'] for v in pv), '1.038', 'Sec. VI'),
              ('A: sqrt(s) [1e6 GeV] at Gamma(0)=0.99', EB['A']['sqrts_Gamma0_0p99'] / 1e6, '8', 'Sec. VI'),
              ('B: sqrt(s) [1e6 GeV] at Gamma(0)=0.99', EB['B']['sqrts_Gamma0_0p99'] / 1e6, '4', 'Sec. VI'),
              ('W measured', UC['W_meas'][0], '5.31', 'Sec. VIII'),
              ('W measured, uncertainty', UC['W_meas'][1], '0.20', 'Sec. VIII'),
              ('W maximum of Gaussian eikonal', UC['W_max'], '3.75', 'Sec. VIII'),
              ('a at W maximum', UC['a_at_W_max'], '6.8', 'Sec. VIII'),
              ('r at W maximum', UC['r_at_W_max'], '0.361', 'Sec. VIII'),
              ('(W_meas - W_max)/dW', -UC['pull_W_max'], '7.7', 'Sec. VIII'),
              ('W black disk j11^2/4', UC['W_black_disk'], '3.67', 'Sec. VIII'),
              ('W at a = 0.05', UC['W_small_a'], '0.44', 'Sec. VIII'),
              ('W at r = 0.281 - 0.009', UC['W_band_r'][0], '3.52', 'Sec. VIII'),
              ('W at r = 0.281 + 0.009', UC['W_band_r'][1], '3.60', 'Sec. VIII'),
              ('W of set B', UC['W_at_a_r'], '3.57', 'Sec. VIII')]
    t4 = {'M1': ['0.100', '106.6', '26.9', '0.102', '20.23', '0.476', '0.659', '1.99', '0.951', '0.119', '0.226', '49', '4569', '4572'],
          'M2': ['0.113', '103.7', '27.5', '0.098', '19.77', '0.469', '0.637', '1.89', '0.971', '0.053', '0.081', '92', '733', '738'],
          'M3': ['0.052', '96.6', '24.6', '0.004', '19.20', '0.468', '0.635', '2.00', '0.934', '0.096', '0.217', '75', '2469', '2560'],
          'M4': ['0.161', '106.1', '29.5', '0.188', '20.08', '0.470', '0.634', '1.75', '1.013', '0.016', '0.122', '100', '71', '169']}
    for m, strs in t4.items():
        k = KF[m]
        vals = [k['params']['Delta'], k['stot'], k['sel'], k['rho'], k['B'], k['tdip'], k['tbump'], k['R'], k['ReGamma0'],
                k['rms_log_lt1'], k['rms_log'], 100 * k['frac10_lt1'], k['chi2_dsdt'], k['chi2_total']]
        names = ['Delta', 'sigma_tot', 'sigma_el', 'rho', 'slope 0.04-0.2', 'first minimum', 'following maximum', 'max/min',
                 'Re Gamma(0)', 'rms ln |t|<1', 'rms ln all', 'bins within 10% |t|<1 [%]', 'chi2_dsigma', 'chi2']
        paper += [(f'{m}: {n}', v, s_, 'Table 4') for n, v, s_ in zip(names, vals, strs)]
    paper += [('M2: gamma', KF['M2']['params']['gamma'], '0.403', 'Table 4'),
              ('M4: gamma', KF['M4']['params']['gamma'], '0.103', 'Table 4'),
              ('M3: H_O [GeV^-2]', KF['M3']['params']['Ho'], '-56.7', 'Table 4'),
              ('M4: H_O [GeV^-2]', KF['M4']['params']['Ho'], '55.4', 'Table 4')]

for k, (s, key) in enumerate(tab2):
    row = D[k % 6]
    paper.append((f'Table 2: {key} at |t|={row["t"]}', row[key], s, 'Table 2'))

def check(val, s, loc):
    dec = len(s.split('.')[1]) if '.' in s else 0
    in_tex = s in tex
    if 'lower bound' in loc:
        ok = val > float(s)
    else:
        ok = abs(round(val, dec) - float(s)) < 0.5 * 10**(-dec) + 1e-12
    return dec, in_tex, ok

rows_p = []
for q, v, s, loc in paper:
    dec, in_tex, ok = check(v, s, loc)
    rows_p.append(dict(quantity=q, program=v, paper=s, location=loc, found_in_tex=in_tex, agrees=ok))

# handoff summary values (quantity, program value, handoff value, note)
handoff = [
 ('A: sigma_tot [mb]', A['stot'], 139.2336, ''), ('A: sigma_el [mb]', A['sel'], 36.5402, ''),
 ('A: sigma_inel [mb]', A['sinel'], 102.6934, ''), ('A: |t_0| [GeV^2]', A['t_zero'], 0.244157, ''),
 ('A: B(0) [GeV^-2]', A['B0'], 25.5336, ''), ('A: excess over TOTEM [%]', 100*V['dstot_rel'], 25.89, ''),
 ('B: a', B['a'], 2.53304174, ''), ('B: nu', B['nu'], 5.06608348, ''), ('B: lambda [GeV^-2]', B['lam'], 7.38433659, ''),
 ('B: sigma_tot [mb]', B['stot'], 110.6000, 'input'), ('B: sigma_el [mb]', B['sel'], 31.0786, ''),
 ('B: sigma_inel [mb]', B['sinel'], 79.5214, ''), ('B: B(0) [GeV^-2]', B['B0'], 18.8314, ''),
 ('B: |t_0| [GeV^2]', B['t_zero'], 0.315492, ''),
 ('B: dsdt at 0.04990', D[0]['setB'], 237.346, 'handoff: bin average; program: point value'),
 ('B: dsdt at 0.10012', D[1]['setB'], 83.3933, 'bin avg vs point'),
 ('B: dsdt at 0.19979', D[2]['setB'], 6.9314, 'bin avg vs point'),
 ('B: dsdt at 0.29960', D[3]['setB'], 0.0371967, 'bin avg vs point'),
 ('B: dsdt at 0.46799', D[4]['setB'], 0.402545, 'bin avg vs point'),
 ('B: dsdt at 0.63179', D[5]['setB'], 0.194537, 'bin avg vs point'),
 ('B: log-lin slope 0.04-0.2', R['loglin_slope_0.04_0.2_B'], 23.0607, 'handoff: 94 data |t| points; program: uniform grid'),
]
rows_h = [dict(quantity=q, program=v, handoff=h, rel_diff=(v - h) / h, note=n) for q, v, h, n in handoff]
json.dump(dict(paper=rows_p, handoff=rows_h), open('comparison.json', 'w'), indent=1, default=float)

print(f'{"quantity":38s} {"program":>14s} {"paper":>10s}  tex  ok   location')
for r in rows_p:
    print(f'{r["quantity"]:38s} {r["program"]:14.7g} {r["paper"]:>10s}  {"Y" if r["found_in_tex"] else "N":>3s}  '
          f'{"Y" if r["agrees"] else "N":>2s}   {r["location"]}')
print(f'\n{"quantity":30s} {"program":>14s} {"handoff":>12s} {"rel.diff":>10s}  note')
for r in rows_h:
    print(f'{r["quantity"]:30s} {r["program"]:14.9g} {r["handoff"]:12.9g} {r["rel_diff"]:10.2e}  {r["note"]}')
nbad = sum(not (r['found_in_tex'] and r['agrees']) for r in rows_p)
print(f'\npaper numbers checked: {len(rows_p)}, mismatches: {nbad}')
