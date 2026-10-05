"""Parameterised re-implementation of the TEHV scenario generator (R1 revision).
Reproduces code/tehv_generate_all_results.py exactly at default settings once the
pulmonary-target parsing bug is fixed (G_pul = 23.37 mmHg instead of 1.40)."""
import numpy as np, pandas as pd

GROUPS = ["baseline", "hydrodynamic_only", "biology_only", "full_data_informed"]
CAL0 = dict(G_pul=23.37, G_ao=32.75, R=11.47, E=2.15, H=58.761932, I12=0.75, Q=4.41)
# endpoint-model coefficients (manuscript Eqs. 4.9-4.19)
COEF0 = dict(g_T=5.5, g_P=4.0, g_M=2.2, g_ao=1.8, e_P=0.35, e_T=0.22, m_M=4.5,
             i_0=-0.2, i_P=1.35, i_S=0.85, i_T=0.25, i_M=0.45, i_H=0.18,
             q_S=0.35, q_P=0.25, r_thin=6.0, r_G=3.0, r_E=2.0,
             f_0=-3.0, f_thin=2.8, f_Mw=1.8, f_G=1.4, f_M=0.9, f_R=0.8)
W0 = dict(I=0.42, Mw=0.18, E=0.15, PG=0.12, F=0.16, M=0.08, PT=0.10)
THR0 = dict(T=0.32, dG=12.0, R=25.0, Mw=65.0, I=0.45, F=0.35)
sig = lambda z: 1/(1+np.exp(-z))

def generate(seed=42, n=1200, cal=None, coef=None, w=None, thr=None,
             equal_noise=False, anchor=True, full_T_lo=0.32, groups=GROUPS):
    cal = {**CAL0, **(cal or {})}; c = {**COEF0, **(coef or {})}
    w = {**W0, **(w or {})}; th = {**THR0, **(thr or {})}
    rng = np.random.default_rng(seed); frames = []
    for g in groups:
        if g == "baseline":
            T=rng.uniform(.10,.95,n); P=rng.uniform(.05,.95,n); H=rng.uniform(8,40,n); M=rng.uniform(.05,.95,n); S=rng.uniform(.05,.95,n)
        elif g == "hydrodynamic_only":
            T=np.clip(rng.normal(.55,.16,n),.18,1.05); P=rng.beta(2,2,n); H=rng.uniform(10,48,n); M=rng.uniform(.05,.95,n); S=rng.uniform(.05,.95,n)
        elif g == "biology_only":
            T=np.clip(rng.normal(.62,.18,n),.22,1.12); P=rng.beta(2.6,1.7,n); H=np.clip(rng.normal(cal["H"],11,n),18,96); M=rng.beta(1.5,4.0,n); S=rng.beta(2.8,1.8,n)
        else:
            T=np.clip(rng.normal(.62,.11,n),full_T_lo,.95); P=rng.beta(2.4,1.9,n); H=np.clip(rng.normal(cal["H"],8,n),30,92); M=rng.beta(1.4,4.8,n); S=rng.beta(2.4,2.1,n)
        full = (g == "full_data_informed") and not equal_noise
        ctx = rng.choice(["pulmonary","aortic"], size=n, p=[.6,.4]); ao = ctx=="aortic"
        Gs = np.where(ao, cal["G_ao"], cal["G_pul"]); Z=(T-.62)/.12
        G = np.clip(Gs + c["g_T"]*Z - c["g_P"]*(P-.55) + c["g_M"]*M + c["g_ao"]*ao + rng.normal(0, 2.0 if full else 3.2, n), 5, 70)
        E = np.clip(cal["E"] + c["e_P"]*(P-.55) - c["e_T"]*Z + rng.normal(0,.08,n), .6, 3.5)
        Mw = np.clip(100*np.exp(-np.log(2)*12/H) - c["m_M"]*M + rng.normal(0,2.5,n), 5, 100)
        I = sig(c["i_0"] + c["i_P"]*P + c["i_S"]*S - c["i_T"]*np.abs(Z) - c["i_M"]*M + c["i_H"]*(H>40) + rng.normal(0, .25 if full else .35, n))
        if g == "full_data_informed":
            eA = rng.normal(0,.03,n)
            if anchor: I = np.clip(.65*I + .35*cal["I12"] + eA, 0, 1)
        Q = np.clip(np.exp(np.log(max(cal["Q"],.2)) + c["q_S"]*S - c["q_P"]*P + rng.normal(0,.35,n)), .2, 12)
        R = np.clip(cal["R"] + c["r_thin"]*(T<.30) + c["r_G"]*(G>Gs+12) - c["r_E"]*(E>cal["E"]) + rng.normal(0,2.5,n), 0, 50)
        F = sig(c["f_0"] + c["f_thin"]*(T<.30) + c["f_Mw"]*(Mw<70) + c["f_G"]*(G>Gs+10) + c["f_M"]*M + c["f_R"]*(R>20) + rng.normal(0,.35,n))
        d = pd.DataFrame(dict(group=g, ctx=ctx, T=T, P=P, H=H, M=M, S=S, Gs=Gs, G=G, E=E, Mw=Mw, I=I, Q=Q, R=R, F=F))
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    return score(d, cal, w, th)

def score(d, cal=None, w=None, th=None):
    cal = {**CAL0, **(cal or {})}; w = {**W0, **(w or {})}; th = {**THR0, **(th or {})}
    d = d.copy()
    PG = np.abs(d.G-d.Gs)/np.maximum(d.Gs,1); PE = np.maximum(0,cal["E"]-d.E)/cal["E"]; PT = np.maximum(0,.32-d["T"])/.32
    d["J"] = np.clip(w["I"]*d.I + w["Mw"]*d.Mw/100 + w["E"]*d.E/3 - w["PG"]*PG - w["F"]*d.F - w["M"]*d.M - w["PT"]*PT, -1, 1)
    c = {"c_T": d["T"]>=th["T"], "c_G": d.G<=d.Gs+th["dG"], "c_R": d.R<=th["R"], "c_Mw": d.Mw>=th["Mw"], "c_I": d.I>=th["I"], "c_F": d.F<=th["F"]}
    for k,v in c.items(): d[k]=v
    d["feasible"] = np.logical_and.reduce(list(c.values())).astype(int)
    return d

def pareto_mask(P):
    eff = np.ones(len(P), bool)
    for i in range(len(P)):
        if eff[i] and np.any(np.all(P>=P[i],1) & np.any(P>P[i],1)): eff[i]=False
    return eff

def pareto(d, top=1800):
    s = d.sort_values("J", ascending=False).head(top).copy()
    cg = (1-np.abs(s.G-s.Gs)/s.Gs).clip(0,1)
    s["pareto"] = pareto_mask(np.c_[s.I, s.Mw/100, cg, 1-s.F])
    return s

def cliffs(x, y):
    from scipy.stats import mannwhitneyu
    U = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2*U/(len(x)*len(y)) - 1
