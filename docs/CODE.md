# 代码说明

本文件说明仓库中每个程序的用途、调用方式、输入输出和主要函数。式号、表号和图号与论文
*Regge–eikonal description of high-energy proton–proton elastic scattering* 一致。
新增程序（第 VI–VIII 节）的公式推导与数值细节另见 [new_programs.pdf](new_programs.pdf)。

## 1. 约定

- 振幅归一化（式 (1)）：$\sigma_{\rm tot}={\rm Im}\,A(s,0)/s$，$d\sigma_{\rm el}/dt=|A|^2/(16\pi s^2)$。
- 碰撞参数表示（式 (10)）：$A=2is\int d^2b\,e^{i\mathbf q\cdot\mathbf b}\,\Gamma(s,b)$，$t=-\mathbf q^2$。
- 程序内部单位为 GeV；截面输出前乘 `GEV2_TO_MB = 0.3893793721`（$1\,{\rm GeV}^{-2}$ 换成 mb）。
- 程序中的 `t` 一律表示 $|t|\ge 0$。
- 数据：`data/TOTEM13_dsdt.csv`，TOTEM 13 TeV 的 290 个 bin（EPJC 79, 861 (2019) 表 5–10），
  列为 `t_low, t_high, t_repr, dsdt, stat, syst`，单位 GeV² 与 mb/GeV²。
- 积分截面输入：$\sigma_{\rm tot}=110.6\pm3.4$ mb，$r=\sigma_{\rm el}/\sigma_{\rm tot}=0.281\pm0.009$（EPJC 79, 103 (2019)）。

## 2. 运行

```bash
pip install -r requirements.txt
./run_all.sh                                             # 全部计算与作图，约 25 分钟
python3 compare.py output/results.json Regge_eikonal_pp.tex   # 与论文逐数比对
```

`run_all.sh` 的执行顺序与各步耗时：

| 步骤 | 命令 | 输出目录 | 耗时 |
|---|---|---|---|
| 1 | `regge_eikonal_pp.py output` | `output/` | 约 1 分钟 |
| 2 | `full_comparison.py data/TOTEM13_dsdt.csv output` | `output/` | 数秒 |
| 3 | `extensions.py …`（`TFIT=1.0` 与 `TFIT=4.1`） | `output/`, `output_t4/` | 数分钟 |
| 4 | `goodwalker.py …`（`TFIT=1.0` 与 `TFIT=4.1`） | `output/`, `output_t4/` | 数分钟 |
| 5 | `profile_inversion.py output` | `output/` | 数秒 |
| 6 | `energy_blackness.py output` | `output/` | 数秒 |
| 7 | `MODELS=M1,M2,M3,M4 NSEED=4 kmr_fits.py output` | `output/` | 约 20 分钟 |
| 8 | `universal_curves.py output` | `output/` | 约 1 分钟 |
| 9 | `figures/make_figures.py`，`figures/make_new_figures.py ../output` | `figures/` | 约 1 分钟 |

每一步的屏幕输出保存在 `logs/`（`01_…` 至 `14_…`）。

## 3. 程序与论文的对应

| 程序 | 论文 | 内容 |
|---|---|---|
| `regge_eikonal_pp.py` | 第 III–V 节；表 1、5、6；图 1–3 的数值 | 单道高斯 Regge–eikonal 的闭合式与数值检验 |
| `full_comparison.py` | 表 2、7；附录 B；图 8 | 与 290 个 bin 的 bin 平均比较 |
| `extensions.py` | 第 VII 节式 (24)–(25)；表 3；图 5 | 复振幅（信号相位）与双分量 opacity |
| `goodwalker.py` | 第 VII 节式 (26) | 高斯 profile 上的两道 Good–Walker eikonal |
| `profile_inversion.py` | 第 VI 节式 (22)–(23)；图 4 | 由数据反演碰撞参数分布 |
| `energy_blackness.py` | 第 VI 节末段 | $\Gamma(0)$ 达到 0.99 的能量 |
| `kmr_fits.py` | 第 VII 节式 (27)–(28)；表 4；图 6 | KMR 形状因子、Odderon、两道 eikonal 的拟合 M1–M4 |
| `universal_curves.py` | 第 VIII 节式 (29)；图 7 | 单高斯 eikonal 的单参数曲线 $(r(a),W(a))$ |
| `figures/make_figures.py` | 图 1、2、3、5、8 | 自含数据的作图程序 |
| `figures/make_new_figures.py` | 图 4、6、7 | 读取 `output/` 的新增结果作图 |
| `compare.py` | 全文 | 215 个论文数字与程序输出的比对 |

## 4. 各程序

### 4.1 `regge_eikonal_pp.py`

单道高斯模型：Born 振幅 $A_B=isHe^{-\lambda|t|}$（式 (7)），opacity $\Omega=\nu e^{-b^2/4\lambda}$，
$\nu=H/(4\pi\lambda)$（式 (11)），$\Gamma=1-e^{-\Omega/2}$（式 (12)）。

| 函数 | 作用 | 式号 |
|---|---|---|
| `Ein(a)` | $\mathrm{Ein}(a)=\sum_{n\ge1}(-1)^{n+1}a^n/(n!\,n)$，mpmath 任意精度求和 | (17) |
| `Ein_closed(a)` | 独立形式 $\gamma_E+\ln a+E_1(a)$，用于检验 | 第 IV 节 |
| `I1(a)` | $I_1(a)=\sum(-1)^{n+1}a^n/(n!\,n^2)$ | (19) |
| `F(Q, a, nmax=250)` | $F(Q,a)=\sum(-1)^{n+1}a^ne^{-Q/n}/(n!\,n)$，用 `gammaln` 避免溢出 | (16) |
| `params_A(sqrts)` | set A 的 $H,\lambda,\nu$；$L=\ln(s/s_{\rm ref})$，$\sqrt{s_{\rm ref}}=7$ TeV，$\lambda_{\rm ref}=10$，$\nu_{\rm ref}=4$，$\Delta=0.0808$，$\alpha'=0.25$ | (7) |
| `dsdt(t, lam, a)` | $4\pi\lambda^2F^2$，单位 mb/GeV² | (19) |
| `dsdt_born(t, H, lam)` | $H^2e^{-2\lambda|t|}/(16\pi)$ | (8) |
| `first_zero(lam, a)` | 在 $0.001$–$3$ GeV² 的 30001 点网格上找 $F$ 第一次变号，再用 Brent 法求 $|t_0|$ | 第 IV 节 |
| `second_max(lam, a, t0)` | $|t_0|$ 之后 $F^2$ 的极大（有界一维极小化） | 第 VI 节 |
| `observables(lam, nu, H)` | $\sigma_{\rm tot},\sigma_{\rm el},\sigma_{\rm inel},r,B(0),|t_0|$，第二极大，$\Gamma(0)$，以及 Born 量 | (17)–(19) |
| `calibrate_B(stot_mb, r)` | 由 $r(a)=1-\mathrm{Ein}(2a)/[2\mathrm{Ein}(a)]$ 解 $a$，再由 $\sigma_{\rm tot}=8\pi\lambda\mathrm{Ein}(a)$ 得 $\lambda$；返回 $(\lambda,\nu)$ | (21), (18) |
| `checks(lam, nu, H)` | 表 6 的全部数值检验：$b$ 空间积分、Hankel 变换、局域幺正性、Born 极限、$B(0)$ 有限差分 | 表 6 |
| `main(outdir)` | 计算 7、8、13 TeV 的 set A、13 TeV 的 set B、六个代表 bin 与能量扫描并写文件 | — |

输出：`results.json`（键 `A_7000GeV`、`A_8000GeV`、`A_13000GeV`、`B_13TeV`、`A_13TeV_vs_TOTEM`、
`checks_A_13TeV` 等）；`differential_13TeV.csv`、`dsdt_curves_13TeV.csv`、`energy_scan.csv`（1–100 TeV）。

### 4.2 `full_comparison.py`

`bin_average(f, lo, hi, n=24)`：24 点 Gauss–Legendre 求 $\frac{1}{\delta t}\int f\,dt$（附录 A 式 (A1)）。
`main(csv, outdir)`：对 290 个 bin 计算 set A、B 的 bin 平均与点值、比值分布、10% 以内的比例、
$0.04\le|t|\le0.2$ GeV² 的等权对数线性斜率（数据 20.39，set B 23.06，set A 34.98 GeV⁻²）。
输出 `full_comparison.json`、`all_290_bins.csv`。

### 4.3 `extensions.py`

Born 振幅写成复分量列表 `comps = [(H_c, λ_c), …]`。

| 函数 | 作用 | 式号 |
|---|---|---|
| `comps_from(H, lam, phase=True)` | 信号相位：$H_c=H(1-i\tan\frac{\pi\Delta}{2})$，$\lambda_c=\lambda-i\pi\alpha'/2$ | (24) |
| `omega(b, comps)` | $\Omega(b)=\sum_jH_{cj}/(4\pi\lambda_{cj})\,e^{-b^2/4\lambda_{cj}}$ | (11), (25) |
| `amp(t, comps)` | $A/(is)=4\pi\int b\,db\,J_0(qb)(1-e^{-\Omega/2})$，$b\le60$ GeV⁻¹，3000 点 Gauss–Legendre | (10) |
| `observables(comps)` | $\sigma_{\rm tot},\sigma_{\rm el},\sigma_{\rm inel},r,\rho,B(0)$，幺正性检验 | (1)–(2) |
| `bin_avg_matrix(lo, hi, n=8)`, `binned` | 每个 bin 8 点 Gauss–Legendre 平均 | (A1) |
| `dip_bump(comps)` | 第一极小、其后极大及比值 | 表 3 |
| `describe(…)` | 对数残差 rms、10% 以内比例、0.04–0.2 GeV² 斜率 | 表 3 |
| `main(csv, outdir)` | set A$_c$、B（实）、B$_c$（校准到 $\sigma_{\rm tot}$ 与 $r$）、C/C′（双分量，$f,k$ 对 $|t|<$ `TFIT` 的 bin 做 $\ln$ 等权拟合） | 第 VII 节 |

环境变量 `TFIT`：1.0 给出 set C（写入 `output/`），4.1 给出 set C′（写入 `output_t4/`）。
输出 `extensions.json`、`extensions_bins.csv`、`extensions_curves.csv`；`check_numeric_vs_series_maxrel`
为数值 Hankel 变换与式 (16) 级数的最大相对差（$4.7\times10^{-12}$）。

### 4.4 `goodwalker.py`

两道 Good–Walker：$\Omega_{ik}=\gamma_i\gamma_k\Omega$，$\gamma_{1,2}=1\pm\gamma$，
$\langle F\rangle=\frac14\sum_{ik}(1-e^{-\Omega_{ik}/2})$（式 (26)），$\Omega$ 为 set B$_c$ 的复 opacity。

| 函数 | 作用 |
|---|---|
| `profiles(H, lam, g)` | $\langle F\rangle(b)$ 与 $\langle|F|^2\rangle(b)$ |
| `obs(H, lam, g)` | $\sigma_{\rm tot},\sigma_{\rm el},r$，低质量衍射截面 $\int(\langle|F|^2\rangle-|\langle F\rangle|^2)$，$\rho$，Pumplin 比 $\max_b\langle|F|^2\rangle/{\rm Re}\langle F\rangle$ |
| `dsdt_t(t, H, lam, g)` | 两道 $d\sigma/dt$ |
| `calib(g, x0)` | 固定 $\gamma$，把 $(H,\lambda)$ 校准到 $\sigma_{\rm tot}$ 与 $r$ |
| `main(csv, outdir, tfit)` | 对 $\gamma$ 做 $\ln d\sigma/dt$ 等权拟合（结果 $\gamma=0$） |
| `gamma_scan(…)` | $\gamma=0,0.2,0.4,0.6$ 的校准参数、低质量衍射截面、极小位置与 dip 处比值 |

输出 `goodwalker.json`、`gw_gamma_scan.json`、`gw_bins.csv`、`gw_curve.csv`。

### 4.5 `profile_inversion.py`

由式 (10) 逆变换得式 (22)：$\Gamma(b)=\frac{1}{8\pi}\int_0^\infty d|t|\,J_0(b\sqrt{|t|})\,A/(is)$，
其中 ${\rm Re}\,\Gamma$ 由 ${\rm Im}A/s=\pm[16\pi\,(d\sigma/dt)/(1+\rho^2)]^{1/2}$ 得到（dip 处变号）。

函数 `profile(stot, B, rho0, tdip, rho_mode, bgrid, yscale, tmax)`：

- $|t|<0.0376$ GeV²：外推 ${\rm Im}A/s=\sigma_{\rm tot}e^{-B|t|/2}$（200 点）。
- 数据区：290 个代表 $|t|$。
- 梯形公式积分，$b$ 取 0–3 fm 的 301 个点。
- $G_{\rm inel}=2{\rm Re}\,\Gamma-(1+\rho^2)({\rm Re}\,\Gamma)^2$（式 (13)）。

参数：`rho_mode='const'` 或 `'lin'`（$\rho_0(1-|t|/0.47)$，截断于 $\pm\rho_0$）；`yscale` 为数据整体缩放；
`tmax` 为积分上限。

输出：
- `profile_inversion.json`：
  - `base`：${\rm Re}\,\Gamma(0)=1.020$，$G_{\rm inel}(0)=0.989$，$\max G_{\rm inel}=0.990$，$b=0.33$ fm；
  - `variations`：九种变化；
  - `check`：$\sigma_{\rm tot}=110.4$ mb，$\sigma_{\rm inel}=79.5$ mb。
- `profile_13TeV.csv`（列 `b_fm, ReGamma, Ginel`）。

### 4.6 `energy_blackness.py`

$\Gamma(0)=1-e^{-a}$，$a(s)=H(s)/(8\pi\lambda(s))$，$H,\lambda$ 随能量按式 (7) 变化。
`parsA(rs)` 用 set A，`parsB(rs)` 以 13 TeV 校准值为参考点，同一 $\Delta$、$\alpha'$。
Brent 法在 $\ln\sqrt s$ 上解 $\Gamma(0)=0.99$。

输出 `energy_blackness.json`：
- 13 TeV 的 $\Gamma(0)$：0.883 与 0.921；
- $\Gamma(0)=0.99$ 的 $\sqrt s$：set A $7.9\times10^6$ GeV，set B $4.0\times10^6$ GeV；
- $d\ln\sigma_{\rm tot}/d\ln s$。

### 4.7 `kmr_fits.py`

Born 振幅（式 (27)–(28)），$i,k=1,2$：

$$\frac{A_{B,ik}}{is}=\gamma_i\gamma_k\Big[H\,e^{-i\pi(\Delta-\alpha'|t|)/2}e^{-\alpha' L|t|}F_iF_k-iH_Oe^{-b_O|t|}\Big],\qquad
F_i=\exp[-(b_i(c_i+|t|))^{d_i}+(b_ic_i)^{d_i}]$$

| 函数 | 作用 |
|---|---|
| `FF(t, b, c, d)` | KMR 形状因子 $F_i$ |
| `omegas(p, two)` | 四个 $\Omega_{ik}(b)=\frac{1}{4\pi}\int d|t|\,J_0(b\sqrt{|t|})A_{B,ik}/(is)$；900 点 Gauss–Legendre，$|t|\le40$ GeV² |
| `gamma_b`, `ginel_b` | $\langle\Gamma\rangle=\frac14\sum(1-e^{-\Omega_{ik}/2})$，$G_{\rm inel}=\frac14\sum(1-e^{-{\rm Re}\,\Omega_{ik}})$ |
| `amp_t(t, p, two)` | $A/(is)=4\pi\int b\,db\,J_0\langle\Gamma\rangle$；1000 点，$b\le35$ GeV⁻¹ |
| `obs(p, two)` | $\sigma_{\rm tot},\sigma_{\rm el},\rho$ |
| `binned(p, two)` | 每个 bin 6 点 Gauss–Legendre 平均 |
| `resid(q, free, p0, two)` | 残差：290 个 $(\langle d\sigma/dt\rangle-y)/\epsilon$，加 $(\sigma_{\rm tot}-110.6)/3.4$ 与 $(\rho-0.09)/0.01$；$\epsilon^2=\epsilon_{\rm stat}^2+\epsilon_{\rm syst}^2+(0.03y)^2$ |
| `describe(p, two)` | 表 4 的全部诊断量 |

参数向量 `p = [H, b1, c1, d1, b2, c2, d2, aP, Delta, Ho, bO, gamma, N]`：
- 起点为 KMR 值；边界见程序中 `lb`、`ub`；
- `N`（$d\sigma/dt$ 归一化）固定为 1；单道模型取 `gamma=0`、$F_2=F_1$。

| 模型 | 自由参数 |
|---|---|
| M1 Pomeron，单道 | `H, b1, c1, d1, aP, Delta` |
| M2 Pomeron，两道 | M1 + `b2, c2, d2, gamma` |
| M3 Pomeron+Odderon，单道 | M1 + `Ho, bO` |
| M4 Pomeron+Odderon，两道 | 全部 12 个 |

环境变量：
- `MODELS`：拟合哪些模型，默认全部；
- `NSEED`：起点数，默认 4。第 0 个为 KMR 值，其余为 KMR 值乘 $1+0.15\xi$；Odderon 的起点 $H_O$ 依次取 0、5、−5、15；
- `FLOOR`：模型误差比例，默认 0.03。

输出 `kmr_fits.json` 与 `kmr_M1.npy … kmr_M4.npy`（最优参数向量）。`kmr_fits.json` 中每个模型的键：
- `stot`、`sel`、`rho`、`B`（0.04–0.2 GeV² 斜率）；
- `tdip`、`tbump`、`R`；
- `ReGamma0`、`Ginel0`；
- `rms_log`、`rms_log_lt1`、`frac10_lt1`；
- `chi2_dsdt`、`chi2_regions`、`chi2_total`、`ndf`；
- `params`。

### 4.8 `universal_curves.py`

$W(a)\equiv\sigma_{\rm tot}|t_0|/(8\pi)=\mathrm{Ein}(a)\,Q_0(a)$，$F(Q_0,a)=0$（式 (29)）；$r(a)$ 为式 (21)。

| 函数 | 作用 |
|---|---|
| `r_of(a)` | $1-\mathrm{Ein}(2a)/[2\mathrm{Ein}(a)]$ |
| `Q0_of(a)` | $F(Q,a)$ 第一零点（8000 点网格找变号 + Brent，级数取 400 项） |
| `W_of(a)` | $\mathrm{Ein}(a)Q_0(a)$ |

$a$ 限于 0.05–20，$a$ 更大时交错级数在双精度下相消误差过大。

输出 `universal_curves.json`：
- 测量值 $W=5.31\pm0.20$；
- 极大 $W=3.75$（$a=6.8$，$r=0.361$），与测量值相差 7.7 个标准差；
- 黑盘值 $j_{1,1}^2/4=3.67$；
- set B 的 $W=3.57$。

以及 `universal_curve.csv`（列 `a, r, W`）。

### 4.9 `figures/make_figures.py` 与 `figures/make_new_figures.py`

`make_figures.py` 内嵌全部数据与模型函数，可单独运行，画图 1（`fig_dsdt`）、2（`fig_sigma`）、
3（`fig_ratio_slope`）、5（`fig_ext`）、8（`fig_full290`）。
`make_new_figures.py [output_dir]` 读取 `output/` 中的新增结果，画图 4（`fig_profile`，含数据反演 profile）、
6（`fig_kmr`）、7（`fig_universal`）。每个图都存 PDF 与 PNG 两份。

### 4.10 `compare.py`

`python3 compare.py output/results.json <论文 tex>`。列表 `paper` 中每一项为
（量名，程序值，论文中印出的字符串，位置）。判定一致的条件有两个：
- 程序值按印出的小数位舍入后等于印出值；
- 该字符串出现在 tex 中。

结果写入 `comparison.json` 并打印汇总，当前为 215 个数字、0 处不一致（`logs/14_manuscript_check.log`）。

## 5. 输出文件一览

| 文件 | 生成程序 | 内容 |
|---|---|---|
| `output/results.json` | regge_eikonal_pp.py | 表 1、5、6 的全部数值 |
| `output/differential_13TeV.csv` | regge_eikonal_pp.py | 六个代表 bin（表 2） |
| `output/dsdt_curves_13TeV.csv` | regge_eikonal_pp.py | 图 1 曲线 |
| `output/energy_scan.csv` | regge_eikonal_pp.py | 图 2、3 的 1–100 TeV 扫描 |
| `output/full_comparison.json`, `all_290_bins.csv` | full_comparison.py | 附录 B |
| `output[_t4]/extensions*.json/csv` | extensions.py | 表 3、图 5 |
| `output[_t4]/goodwalker.json`, `gw_*` | goodwalker.py | 第 VII 节两道（高斯） |
| `output/profile_inversion.json`, `profile_13TeV.csv` | profile_inversion.py | 式 (23)、图 4 |
| `output/energy_blackness.json` | energy_blackness.py | 第 VI 节末段 |
| `output/kmr_fits.json`, `kmr_M*.npy` | kmr_fits.py | 表 4、图 6 |
| `output/universal_curves.json`, `universal_curve.csv` | universal_curves.py | 第 VIII 节、图 7 |
| `comparison.json` | compare.py | 论文数字比对 |
