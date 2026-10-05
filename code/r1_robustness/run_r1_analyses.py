"""R1 robustness / sensitivity / held-out analyses. Writes CSVs to ./tables and a JSON summary."""
import numpy as np, pandas as pd, json, sys
from scipy import stats
sys.path.insert(0,'.'); from tehv_core import *
import os; os.makedirs('tables', exist_ok=True)
OUT = {}
def grp(d): return d.groupby("group", sort=False).agg(J=("J","mean"), feas=("feasible","mean"))
def order_ok(s): return bool(np.all(np.diff(s.loc[GROUPS].values) > 0))

# ---------- 0. corrected baseline run -------------------------------------
d0 = generate(); d0.to_csv('tables/R1_candidates_corrected.csv', index=False)
summ = d0.groupby("group", sort=False).agg(**{f"{m}_{s}":(m,s) for m in ["J","G","E","Mw","I","R","F","feasible"] for s in ["mean","std","median"]})
summ.to_csv('tables/R1_T1_scenario_summary.csv'); OUT["summary"] = summ.round(4).to_dict("index")
rows=[]; rng=np.random.default_rng(42)
for m in ["J","G","I","Mw","R","F","feasible"]:
    H,p = stats.kruskal(*[d0[d0.group==g][m] for g in GROUPS])
    eps2 = (H-3)/(len(d0)-4)
    for a in GROUPS[1:]:
        x=d0[d0.group==a][m].values; y=d0[d0.group=="baseline"][m].values
        bs=[rng.choice(x,len(x)).mean()-rng.choice(y,len(y)).mean() for _ in range(2000)]
        rows.append(dict(endpoint=m, comparison=f"{a} vs baseline", mean_diff=x.mean()-y.mean(), ci_lo=np.percentile(bs,2.5), ci_hi=np.percentile(bs,97.5), cliffs_delta=cliffs(x,y), KW_H=H, KW_p=p, epsilon2=eps2))
eff=pd.DataFrame(rows); eff.to_csv('tables/R1_T2_effect_sizes.csv', index=False)
top=d0[d0.feasible==1].sort_values("J",ascending=False).head(50)
OUT["top50"]=dict(ctx=top.ctx.value_counts().to_dict(), group=top.group.value_counts().to_dict(), means=top[["T","P","H","M","G","E","Mw","I","R","F","J"]].agg(["mean","std"]).round(4).to_dict(), best={k:(round(float(v),4) if not isinstance(v,(str,bool,np.bool_)) else str(v)) for k,v in top.iloc[0].items()})
p=pareto(d0); OUT["pareto"]=dict(n=int(p.pareto.sum()), by_group=p[p.pareto].group.value_counts().to_dict(), in_top1800=p.group.value_counts().to_dict())
feats=["T","P","H","M","S","G","E","Mw","R","F"]
OUT["spearman"]={f: round(float(stats.spearmanr(d0[f],d0.J)[0]),3) for f in feats}
# constraint-wise pass rates (feasibility decomposition)
cw=d0.groupby("group",sort=False)[["c_T","c_G","c_R","c_Mw","c_I","c_F","feasible"]].mean(); cw.to_csv('tables/R1_T3_constraint_pass_rates.csv'); OUT["constraint_pass"]=cw.round(4).to_dict("index")

# ---------- 1. seed robustness --------------------------------------------
S=[]
for s in range(200):
    g=grp(generate(seed=s)); S.append(dict(seed=s, **{f"feas_{k}":g.feas[k] for k in GROUPS}, **{f"J_{k}":g.J[k] for k in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
S=pd.DataFrame(S); S.to_csv('tables/R1_T4_seed_robustness.csv', index=False)
OUT["seeds"]={c: [round(S[c].mean(),4), round(S[c].quantile(.025),4), round(S[c].quantile(.975),4)] for c in S.columns if c[:4] in("feas","J_ba","J_hy","J_bi","J_fu")}
OUT["seeds"]["order_feas_frac"]=S.order_feas.mean(); OUT["seeds"]["order_J_frac"]=S.order_J.mean()

# ---------- 2. one-at-a-time coefficient sensitivity ----------------------
R=[]
for k,v in COEF0.items():
    for f in [.5,.8,1.2,1.5]:
        g=grp(generate(coef={k:v*f}))
        R.append(dict(coef=k, base=v, factor=f, **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
R=pd.DataFrame(R); R.to_csv('tables/R1_T5_coef_OAT.csv', index=False)
b=grp(d0)
oat=R.groupby("coef").agg(full_feas_min=("feas_full_data_informed","min"), full_feas_max=("feas_full_data_informed","max"), base_feas_min=("feas_baseline","min"), base_feas_max=("feas_baseline","max"), J_full_min=("J_full_data_informed","min"), J_full_max=("J_full_data_informed","max"), order_feas=("order_feas","all"), order_J=("order_J","all"))
oat["gap_min"]=R.assign(gap=R.feas_full_data_informed-R.feas_baseline).groupby("coef").gap.min(); oat["gap_max"]=R.assign(gap=R.feas_full_data_informed-R.feas_baseline).groupby("coef").gap.max()
oat.to_csv('tables/R1_T5b_coef_OAT_summary.csv'); OUT["oat"]=dict(all_order_feas=bool(R.order_feas.all()), all_order_J=bool(R.order_J.all()), n=len(R), gap_range=[float(oat.gap_min.min()), float(oat.gap_max.max())], most_influential=oat.assign(span=oat.gap_max-oat.gap_min).sort_values("span",ascending=False).head(6).round(4)[["gap_min","gap_max","full_feas_min","full_feas_max","base_feas_min","base_feas_max"]].to_dict("index"))

# ---------- 3. global (joint) coefficient perturbation --------------------
rng=np.random.default_rng(2026); Gl=[]
for lvl in [.2,.5]:
    for i in range(1000):
        cf={k:v*rng.uniform(1-lvl,1+lvl) for k,v in COEF0.items()}
        g=grp(generate(coef=cf)); Gl.append(dict(level=lvl, **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J), full_gt_base=g.feas[GROUPS[3]]>g.feas[GROUPS[0]]))
Gl=pd.DataFrame(Gl); Gl.to_csv('tables/R1_T6_coef_global.csv', index=False)
OUT["global"]={str(l): dict(order_feas=float(x.order_feas.mean()), order_J=float(x.order_J.mean()), full_gt_base=float(x.full_gt_base.mean()), **{c:[round(x[c].median(),4), round(x[c].quantile(.025),4), round(x[c].quantile(.975),4)] for c in x.columns if c.startswith(("feas_","J_"))}) for l,x in Gl.groupby("level")}

# ---------- 4. objective-weight sensitivity -------------------------------
schemes = {"original":W0,
 "equal":dict(I=1/7,Mw=1/7,E=1/7,PG=1/7,F=1/7,M=1/7,PT=1/7),
 "hydrodynamic_emphasis":dict(I=.20,Mw=.10,E=.25,PG=.25,F=.10,M=.05,PT=.05),
 "biology_emphasis":dict(I=.55,Mw=.20,E=.05,PG=.05,F=.08,M=.04,PT=.03),
 "safety_emphasis":dict(I=.20,Mw=.15,E=.10,PG=.10,F=.30,M=.10,PT=.05),
 "degradation_emphasis":dict(I=.20,Mw=.40,E=.10,PG=.10,F=.10,M=.05,PT=.05)}
tot=sum(W0.values()); J0=d0.J.values; feas=d0.feasible.values==1; top0=set(top.index); W=[]
def wrow(name,w):
    w={k:v*tot/sum(w.values()) for k,v in w.items()}; d=score(d0,w=w); g=d.groupby("group",sort=False).J.mean()
    t=d[d.feasible==1].sort_values("J",ascending=False).head(50)
    return dict(scheme=name, **{f"w_{k}":round(v,3) for k,v in w.items()}, **{f"J_{x}":g[x] for x in GROUPS}, order_J=order_ok(g), rho_all=stats.spearmanr(J0,d.J)[0], rho_feasible=stats.spearmanr(J0[feas],d.J.values[feas])[0], tau_feasible=stats.kendalltau(J0[feas],d.J.values[feas])[0], top50_overlap=len(top0&set(t.index))/50, top50_full=(t.group==GROUPS[3]).mean(), top50_bio=(t.group==GROUPS[2]).mean(), top50_aortic=(t.ctx=="aortic").mean(), top1_same=t.index[0]==top.index[0])
for n,w in schemes.items(): W.append(wrow(n,w))
pd.DataFrame(W).to_csv('tables/R1_T7_weight_schemes.csv', index=False); OUT["weights_named"]=pd.DataFrame(W).round(3).set_index("scheme")[["J_baseline","J_hydrodynamic_only","J_biology_only","J_full_data_informed","order_J","rho_all","rho_feasible","tau_feasible","top50_overlap","top50_full","top50_bio","top50_aortic","top1_same"]].to_dict("index")
rng=np.random.default_rng(7); Wr=[]
for lvl in [.2,.5]:
    for i in range(1000): Wr.append(dict(level=lvl, **wrow("rand",{k:v*rng.uniform(1-lvl,1+lvl) for k,v in W0.items()})))
for i in range(1000): Wr.append(dict(level="dirichlet", **wrow("dir",dict(zip(W0, rng.dirichlet(np.ones(7)))))))
Wr=pd.DataFrame(Wr); Wr.to_csv('tables/R1_T8_weight_random.csv', index=False)
OUT["weights_random"]={str(l): dict(order_J=float(x.order_J.mean()), full_best=float((x.J_full_data_informed>=x[[f"J_{g}" for g in GROUPS]].max(axis=1)).mean()), base_worst=float((x.J_baseline<=x[[f"J_{g}" for g in GROUPS]].min(axis=1)).mean()), **{c:[round(x[c].median(),3), round(x[c].quantile(.025),3), round(x[c].quantile(.975),3)] for c in ["rho_all","rho_feasible","tau_feasible","top50_overlap","top50_full","top50_bio"]}, top1_same=float(x.top1_same.mean())) for l,x in Wr.groupby("level", sort=False)}

# ---------- 5. feasibility-threshold sensitivity --------------------------
Th=[]
for k,v in THR0.items():
    for f in [.8,.9,1.1,1.2]:
        g=grp(score(d0,th={k:v*f})); Th.append(dict(threshold=k, base=v, factor=f, **{f"feas_{x}":g.feas[x] for x in GROUPS}, order=order_ok(g.feas)))
Th=pd.DataFrame(Th); Th.to_csv('tables/R1_T9_threshold_sensitivity.csv', index=False); OUT["thresholds"]=dict(all_order=bool(Th.order.all()), full_range=[float(Th.feas_full_data_informed.min()), float(Th.feas_full_data_informed.max())], base_range=[float(Th.feas_baseline.min()), float(Th.feas_baseline.max())])

# ---------- 6. ablation: what drives the full-data gain? ------------------
A=[]
for name,kw in [("full (as reported)",{}),("no integration anchoring",dict(anchor=False)),("equal noise + no anchoring",dict(equal_noise=True,anchor=False)),("equal noise + no anchoring + T truncated at 0.22",dict(equal_noise=True,anchor=False,full_T_lo=.22))]:
    d=generate(**kw); x=d[d.group==GROUPS[3]]; pp=pareto(d)
    A.append(dict(variant=name, feas=x.feasible.mean(), J=x.J.mean(), I=x.I.mean(), I_sd=x.I.std(), I_p95=x.I.quantile(.95), F=x.F.mean(), pareto_full=int((pp[pp.pareto].group==GROUPS[3]).sum()), pareto_bio=int((pp[pp.pareto].group==GROUPS[2]).sum()), top50_full=int((d[d.feasible==1].sort_values("J",ascending=False).head(50).group==GROUPS[3]).sum())))
A=pd.DataFrame(A); A.to_csv('tables/R1_T10_ablation_full_group.csv', index=False); OUT["ablation"]=A.round(4).set_index("variant").to_dict("index")
OUT["I_tail"]={g: [round(float(d0[d0.group==g].I.quantile(q)),3) for q in (.5,.95,.99)] + [round(float(d0[d0.group==g].I.max()),3)] for g in GROUPS}

# ---------- 7. held-out / leave-one-out checks on calibration -------------
t=np.array([1,3,6,9,12.]); mw=np.array([81.317,75.5861,73.5892,66.792,72.9166]); L=[]
for i in range(5):
    k=np.arange(5)!=i; s,c=np.polyfit(t[k],np.log(mw[k]),1); pred=np.exp(c+s*t[i]); L.append(dict(held_out_month=t[i], observed=mw[i], predicted=pred, abs_pct_err=abs(pred-mw[i])/mw[i]*100, half_life=np.log(2)/abs(s)))
L=pd.DataFrame(L); L.to_csv('tables/R1_T11_LOTO_Mw.csv', index=False)
s,c=np.polyfit(t,np.log(mw),1)
gr=np.array([31,31,36,33.]); LG=[dict(held_out=m, observed=gr[i], predicted=np.delete(gr,i).mean(), abs_err=abs(gr[i]-np.delete(gr,i).mean())) for i,m in enumerate([3,6,9,12])]
pd.DataFrame(LG).to_csv('tables/R1_T12_LOTO_gradient.csv', index=False)
yac=np.array([84.25,86.65,78.82]); pred6=100*np.exp(-np.log(2)*6/CAL0["H"]); H_yac=6*np.log(2)/np.log(100/yac.mean())
OUT["heldout"]=dict(LOTO_Mw_MAPE=float(L.abs_pct_err.mean()), LOTO_Mw_max=float(L.abs_pct_err.max()), LOTO_halflife_range=[float(L.half_life.min()), float(L.half_life.max())], full_fit_R2=float(np.corrcoef(t,np.log(mw))[0,1]**2), LOTO_grad_MAE=float(np.mean([x["abs_err"] for x in LG])), LOTO_grad_max=float(max(x["abs_err"] for x in LG)), yacoub_obs_6m=float(yac.mean()), yacoub_range=[float(yac.min()),float(yac.max())], model_pred_6m=float(pred6), cross_study_err_pp=float(pred6-yac.mean()), H_yacoub_implied=float(H_yac))
# leave-one-study-out on the degradation calibration: recalibrate H from each LOTO fit and from Yacoub only
LS=[]
for name,H in [("Vis GPC (as reported)",CAL0["H"]),("LOTO min half-life",L.half_life.min()),("LOTO max half-life",L.half_life.max()),("Yacoub-only (6-month retention)",H_yac),("Vis Mn-based",42.2)]:
    g=grp(generate(cal=dict(H=H))); d=generate(cal=dict(H=H))
    LS.append(dict(calibration=name, H=H, **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, Mw_full=d[d.group==GROUPS[3]].Mw.mean(), order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
# hydrodynamic targets +-1 SD / LOTO extremes
for name,cal in [("G_pul -1SD",dict(G_pul=23.37-1.11)),("G_pul +1SD",dict(G_pul=23.37+1.11)),("G_ao LOTO min",dict(G_ao=31.67)),("G_ao LOTO max",dict(G_ao=33.33)),("E -1SD",dict(E=2.09)),("E +1SD",dict(E=2.21)),("R -1SD",dict(R=10.89)),("R +1SD",dict(R=12.05)),("I12 = 0.44 (MC 2.5%)",dict(I12=.44)),("I12 = 1.00 (MC 97.5%)",dict(I12=1.0))]:
    d=generate(cal=cal); g=grp(d)
    LS.append(dict(calibration=name, H=CAL0["H"], **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, Mw_full=d[d.group==GROUPS[3]].Mw.mean(), order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
LS=pd.DataFrame(LS); LS.to_csv('tables/R1_T13_calibration_target_sensitivity.csv', index=False); OUT["cal_sens"]=LS.round(4).set_index("calibration").to_dict("index")
json.dump(OUT, open('tables/R1_summary.json','w'), indent=1, default=float)
print(json.dumps(OUT, indent=1, default=float))
