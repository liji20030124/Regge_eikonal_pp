r"""Data-driven impact-parameter profile at 13 TeV (Sec. VI): Eqs. (22)-(23), inversion of Eq. (10)
with the TOTEM dsigma/dt.  Usage: python3 profile_inversion.py [output_dir]
Convention (repo): A = 2 i s \int d^2b e^{iqb} Gamma(b), sigma_tot = ImA(0)/s, dsigma/dt = |A|^2/(16 pi s^2).
Inverse: Gamma(b) = (1/(8 pi)) \int dt J0(b sqrt t) A(t)/(i s).
Im A/s = sign(t) sqrt(16 pi dsdt/(1+rho(t)^2)); sign flips at the dip.
"""
import json, os, sys, numpy as np
OUT = sys.argv[1] if len(sys.argv) > 1 else 'output'; os.makedirs(OUT, exist_ok=True)
from scipy.special import j0
GEV2MB=0.3893793721
d=np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','TOTEM13_dsdt.csv'),delimiter=',',skiprows=1)
lo,hi,tr,y,st,sy=d.T
def profile(stot=110.6,B=20.40,rho0=0.10,tdip=0.47,rho_mode='const',bgrid=np.linspace(0,3.0,301)/0.1973,yscale=1.0,tmax=None):
    # low-t extrapolation 0<t<t1 : exponential with slope B normalised to optical point
    st_g=stot/GEV2MB
    t1=tr[0]
    tt=np.linspace(0,t1,200)
    imA_low=st_g*np.exp(-B*tt/2)            # Im A/s at small t (B defined on dsigma/dt)
    # data region
    ys=y*yscale
    if tmax is not None:
        m=tr<=tmax
    else: m=np.ones_like(tr,bool)
    rho_t=np.full(m.sum(),rho0) if rho_mode=='const' else rho0*np.clip(1-tr[m]/tdip,-1,1)
    mod=np.sqrt(16*np.pi*ys[m]/GEV2MB/(1+rho_t**2))
    sgn=np.where(tr[m]<tdip,1.0,-1.0)
    imA_d=sgn*mod
    # match low-t piece to data continuity (rescale exponential) -- report mismatch
    mis=imA_low[-1]/imA_d[0]
    T=np.concatenate([tt,tr[m]]); IA=np.concatenate([imA_low,imA_d])
    o=np.argsort(T); T=T[o]; IA=IA[o]
    # trapezoid in t
    G=np.array([np.trapezoid(j0(b*np.sqrt(T))*IA,T)/(8*np.pi) for b in bgrid])
    Ginel=2*G-G**2*(1+rho0**2)   # |Gamma|^2 approx (1+rho^2) ReGamma^2
    return dict(b_fm=bgrid*0.1973,ReGamma=G,Ginel=Ginel,lowt_mismatch=mis)
res={}
base=profile()
res['base']=dict(ReGamma0=float(base['ReGamma'][0]),Ginel0=float(base['Ginel'][0]),
  Ginel_max=float(base['Ginel'].max()),b_at_max_fm=float(base['b_fm'][np.argmax(base['Ginel'])]),
  lowt_mismatch=float(base['lowt_mismatch']))
var={}
for name,kw in dict(stot_up=dict(stot=114.0),stot_dn=dict(stot=107.2),rho0=dict(rho0=0.0),rho_lin=dict(rho_mode='lin'),
                    B_up=dict(B=20.60),B_dn=dict(B=20.20),norm_up=dict(yscale=1.055),norm_dn=dict(yscale=0.945),tmax2=dict(tmax=2.0)).items():
    p=profile(**kw)
    var[name]=dict(ReGamma0=float(p['ReGamma'][0]),Ginel0=float(p['Ginel'][0]),Ginel_max=float(p['Ginel'].max()),
                   b_at_max_fm=float(p['b_fm'][np.argmax(p['Ginel'])]),hollow=float(p['Ginel'].max()-p['Ginel'][0]))
res['variations']=var
# integrals check: sigma_tot=2\int d^2b ReGamma, sigma_inel=\int d^2b Ginel
bg=base['b_fm']/0.1973
res['check']=dict(stot_mb=float(2*np.trapezoid(2*np.pi*bg*base['ReGamma'],bg)*GEV2MB),
                  sinel_mb=float(np.trapezoid(2*np.pi*bg*base['Ginel'],bg)*GEV2MB))
np.savetxt(os.path.join(OUT,'profile_13TeV.csv'),np.c_[base['b_fm'],base['ReGamma'],base['Ginel']],delimiter=',',header='b_fm,ReGamma,Ginel',comments='')
json.dump(res,open(os.path.join(OUT,'profile_inversion.json'),'w'),indent=1); print(json.dumps(res,indent=1))
