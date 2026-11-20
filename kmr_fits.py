r"""Models M1-M4 of Sec. VII, Eqs. (26)-(28), and the numbers of Table 4.
Usage: [MODELS=M1,M2,M3,M4] [NSEED=4] python3 kmr_fits.py [output_dir]
KMR-type two-channel eikonal with eigenstate-dependent form factors, real part and Odderon.
Born for eigenstate pair (i,k), |t| written t >= 0:
 A_ik/(is) = g_i g_k [ H e^{-i pi (Delta - aP t)/2} e^{-aP L t} F_i(t) F_k(t)  - i Ho e^{-bO t} ],  g_{1,2} = 1 +- gamma
 F_i(t) = exp(-(b_i (c_i + t))^{d_i} + (b_i c_i)^{d_i})   [form of Ref. KMR]
 Omega_ik(b) = (1/4pi) int dt J0(b sqrt t) A_ik/(is),   <Gamma> = (1/4) sum_ik (1 - e^{-Omega_ik/2})
One-channel models: gamma = 0 and F_1 = F_2.  The normalisation N of dsigma/dt is fixed to 1;
errors: stat (+) syst (+) FLOOR*y.  Extra data: sigma_tot = 110.6 +- 3.4 mb, rho = 0.09 +- 0.01.
"""
import json, os, sys, numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0
from scipy.optimize import least_squares
GEV2MB=0.3893793721; L=np.log(13000.**2)
d=np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','TOTEM13_dsdt.csv'),delimiter=',',skiprows=1)
lo,hi,tr,y,st,sy=d.T; FLOOR=float(os.environ.get('FLOOR','0.03')); err=np.sqrt(st**2+sy**2+(FLOOR*y)**2)
x8,w8=leggauss(6); T=0.5*(hi-lo)[:,None]*x8+0.5*(hi+lo)[:,None]; W=0.5*w8[None,:]*np.ones_like(T)
NB=1000; xb,wb=leggauss(NB); BMAX=35.; bb=0.5*BMAX*(xb+1); wbb=0.5*BMAX*wb
NT=900; xt,wt=leggauss(NT); TMAX=40.; tg=0.5*TMAX*(xt+1); wtg=0.5*TMAX*wt
Kt2b=j0(np.sqrt(tg)[None,:]*bb[:,None])*wtg[None,:]/(4*np.pi)
Kb2t=4*np.pi*j0(np.sqrt(T.ravel())[:,None]*bb[None,:])*(bb*wbb)[None,:]
NAMES=['H','b1','c1','d1','b2','c2','d2','aP','Delta','Ho','bO','gamma','N']
def FF(t,b,c,dd): return np.exp(-(b*(c+t))**dd+(b*c)**dd)
def omegas(p,two):
    H,b1,c1,d1,b2,c2,d2,aP,Dl,Ho,bO,g,N=p
    reg=H*np.exp(-1j*np.pi*(Dl-aP*tg)/2)*np.exp(-aP*L*tg); odd=-1j*Ho*np.exp(-bO*tg)
    F1=FF(tg,b1,c1,d1); F2=FF(tg,b2,c2,d2) if two else F1
    g=g if two else 0.0
    gs=[1+g,1-g]; Fs=[F1,F2]; out=[]
    for i in range(2):
        for k in range(2):
            out.append(gs[i]*gs[k]*(Kt2b@(reg*Fs[i]*Fs[k]+odd)))
    return out
def gamma_b(p,two): return 0.25*sum(1-np.exp(-O/2) for O in omegas(p,two))
def ginel_b(p,two): return 0.25*sum(1-np.exp(-O.real) for O in omegas(p,two))
def amp_t(t,p,two):
    K=4*np.pi*j0(np.sqrt(np.atleast_1d(t))[:,None]*bb[None,:])*(bb*wbb)[None,:]; return K@gamma_b(p,two)
def obs(p,two):
    G=gamma_b(p,two); A0=4*np.pi*np.sum(bb*wbb*G)
    return A0.real*GEV2MB, 2*np.pi*np.sum(bb*wbb*abs(G)**2)*GEV2MB, -A0.imag/A0.real, G
def binned(p,two):
    a=(Kb2t@gamma_b(p,two)).reshape(T.shape); return (abs(a)**2/(16*np.pi)*GEV2MB*W).sum(1)
def resid(q,free,p0,two):
    p=p0.copy(); p[free]=q; m=p[12]*binned(p,two); stot,sel,rho,_=obs(p,two)
    return np.concatenate([(m-y)/err,[(stot-110.6)/3.4,(rho-0.09)/0.01]])
def describe(p,two):
    stot,sel,rho,G=obs(p,two); m=p[12]*binned(p,two)
    tt=np.linspace(0.25,1.0,1501); v=abs(amp_t(tt,p,two))**2/(16*np.pi)*GEV2MB
    i=np.where((v[1:-1]<v[:-2])&(v[1:-1]<v[2:]))[0]
    if len(i): k=i[0]+1; kb=k+int(np.argmax(v[k:])); dip=dict(tdip=float(tt[k]),tbump=float(tt[kb]),R=float(v[kb]/v[k]))
    else: dip=dict(tdip=None,tbump=None,R=None)
    sl=(tr>=0.04)&(tr<=0.2); B=float(-np.polyfit(tr[sl],np.log(m[sl]),1)[0])
    Gi=ginel_b(p,two); pull=(m-y)/err
    reg=lambda a,b_: float(np.sum(pull[(tr>=a)&(tr<b_)]**2))
    return dict(stot=float(stot),sel=float(sel),rho=float(rho),r=float(sel/stot),B=B,**dip,
                chi2_dsdt=float(np.sum(pull**2)),chi2_regions=dict(lt03=reg(0,0.3),dip=reg(0.3,0.8),gt08=reg(0.8,5)),
                ReGamma0=float(G[0].real),Ginel0=float(Gi[0]),Ginel_max=float(Gi.max()),
                rms_log=float(np.sqrt(np.mean(np.log(m/y)**2))),rms_log_lt1=float(np.sqrt(np.mean(np.log(m[tr<1]/y[tr<1])**2))),N=float(p[12]),
                frac10_lt1=float(np.mean(abs(m[tr<1]/y[tr<1]-1)<0.1)),frac10_all=float(np.mean(abs(m/y-1)<0.1)))
#        H   b1  c1    d1    b2  c2   d2   aP    Dl   Ho  bO  g    N    (start: KMR values)
p0=np.array([700,10.,0.233,0.462,4.9,0.52,0.47,0.052,0.13,0.0,6.0,0.56,1.0])
lb=np.array([50,0.5,0.01,0.2,0.5,0.01,0.2,0.0,0.0,-200,0.5,0.0,0.85])
ub=np.array([5000,20,3,1.2,20,3,1.2,0.4,0.3,200,30,0.95,1.15])
models={'M1':('Pomeron, one channel',[0,1,2,3,7,8],False),
        'M2':('Pomeron, two channels',[0,1,2,3,4,5,6,7,8,11],True),
        'M3':('Pomeron + Odderon, one channel',[0,1,2,3,7,8,9,10],False),
        'M4':('Pomeron + Odderon, two channels',[0,1,2,3,4,5,6,7,8,9,10,11],True)}
sel_models=os.environ.get('MODELS','M1,M2,M3,M4').split(',')
OUT=sys.argv[1] if len(sys.argv)>1 else 'output'; os.makedirs(OUT,exist_ok=True)
OJ=os.path.join(OUT,'kmr_fits.json'); out=json.load(open(OJ)) if os.path.exists(OJ) else {}
for key in sel_models:
    lab,free,two=models[key]; best=None
    for seed in range(int(os.environ.get('NSEED','4'))):
        rng=np.random.default_rng(200+seed); q0=p0[free]*(1+(0.15 if seed else 0)*rng.standard_normal(len(free)))
        if 9 in free: q0[free.index(9)]=[0,5,-5,15][seed%4]
        q0=np.clip(q0,lb[free]+1e-6,ub[free]-1e-6)
        r=least_squares(resid,q0,args=(free,p0,two),bounds=(lb[free],ub[free]),x_scale='jac',max_nfev=2500)
        if best is None or r.cost<best.cost: best=r
        print(key,'seed',seed,'chi2',round(2*r.cost,1),flush=True)
    p=p0.copy(); p[free]=best.x
    D=describe(p,two); D.update(chi2_total=float(2*best.cost),nfree=len(free),ndf=len(y)+2-len(free),label=lab,params=dict(zip(NAMES,map(float,p))),two=two)
    out[key]=D; np.save(os.path.join(OUT,f'kmr_{key}.npy'),p); json.dump(out,open(OJ,'w'),indent=1)
    print(key,lab,json.dumps({k:(round(v,4) if isinstance(v,float) else v) for k,v in D.items() if k!='params'}),flush=True)
    print('  params',{k:round(v,4) for k,v in D['params'].items()},flush=True)
