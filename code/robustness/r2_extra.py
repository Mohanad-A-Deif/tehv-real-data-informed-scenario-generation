# ---------- R2 additions (executed inside run_r2_analyses.py) ---------------
# 8. coefficient influence classification (from OAT): largest change in any group's feasibility or mean J
b0=grp(d0); FE=[f"feas_{x}" for x in GROUPS]; JJ=[f"J_{x}" for x in GROUPS]
ref_f=np.array([b0.feas[x] for x in GROUPS]); ref_J=np.array([b0.J[x] for x in GROUPS]); ref_gap=ref_f[3]-ref_f[0]
infl=pd.DataFrame([dict(coef=k, max_dJ_any_group=np.abs(x[JJ].values-ref_J).max(), max_dfeas_any_group=np.abs(x[FE].values-ref_f).max(),
        max_dgap=np.abs((x[FE[3]]-x[FE[0]]).values-ref_gap).max()) for k,x in R.groupby("coef", sort=False)])
STRUCT_INERT=["q_S","q_P"]   # hinge-to-belly model: enters neither J nor the feasibility rule
infl["cls"]=np.select([infl.coef.isin(STRUCT_INERT), (infl.max_dfeas_any_group<=0.005)&(infl.max_dJ_any_group<=0.002),
                       (infl.max_dfeas_any_group<=0.01)&(infl.max_dJ_any_group<=0.01)], ["structurally inert","negligible","weak"], "active")
infl=infl.sort_values(["max_dfeas_any_group","max_dJ_any_group"], ascending=False)
infl.to_csv('tables/R2_T14_coefficient_influence.csv', index=False)
OUT["influence"]=dict(counts=infl.cls.value_counts().to_dict(), **{c:list(infl.coef[infl.cls==c]) for c in ["active","weak","negligible","structurally inert"]})
OUT["influence"]["inert"]=OUT["influence"]["structurally inert"]
# 9. sign reversal of each non-inert coefficient
SF=[]
for k,v in ACTIVE.items():
    if k in OUT["influence"]["inert"]: continue
    g=grp(generate(coef={k:-v})); SF.append(dict(coef=k, base=v, **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J), full_gt_base=g.feas[GROUPS[3]]>g.feas[GROUPS[0]]))
SF=pd.DataFrame(SF); SF.to_csv('tables/R2_T15_sign_reversal.csv', index=False)
OUT["sign_flip"]=dict(n=len(SF), order_feas=int(SF.order_feas.sum()), order_J=int(SF.order_J.sum()), full_gt_base=int(SF.full_gt_base.sum()), broken=SF.loc[~(SF.order_feas&SF.order_J),["coef","feas_baseline","feas_hydrodynamic_only","feas_biology_only","feas_full_data_informed","J_biology_only","J_full_data_informed"]].round(4).to_dict("records"))
# 10. internal cut-points (one at a time and jointly)
CP=[]
for k,v in CUT0.items():
    if k=="Tref": continue
    for f in [.8,.9,1.1,1.2]:
        g=grp(generate(cut={k:v*f})); CP.append(dict(cut=k, base=v, factor=f, **{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
CP=pd.DataFrame(CP); CP.to_csv('tables/R2_T16_cutpoint_OAT.csv', index=False)
rng=np.random.default_rng(77); CJ=[]
for i in range(1000):
    ct={k:v*rng.uniform(.8,1.2) for k,v in CUT0.items() if k!="Tref"}
    g=grp(generate(cut=ct)); CJ.append(dict(**{f"feas_{x}":g.feas[x] for x in GROUPS}, **{f"J_{x}":g.J[x] for x in GROUPS}, order_feas=order_ok(g.feas), order_J=order_ok(g.J)))
CJ=pd.DataFrame(CJ); CJ.to_csv('tables/R2_T17_cutpoint_joint.csv', index=False)
OUT["cutpoints"]=dict(oat_n=len(CP), oat_order_feas=int(CP.order_feas.sum()), oat_order_J=int(CP.order_J.sum()), oat_full_range=[float(CP.feas_full_data_informed.min()),float(CP.feas_full_data_informed.max())], oat_base_range=[float(CP.feas_baseline.min()),float(CP.feas_baseline.max())],
  oat_broken=CP.loc[~(CP.order_feas&CP.order_J),["cut","factor","feas_hydrodynamic_only","feas_biology_only","feas_full_data_informed","J_biology_only","J_full_data_informed"]].round(4).to_dict("records"),
  joint=dict(order_feas=float(CJ.order_feas.mean()), order_J=float(CJ.order_J.mean()), **{c:[round(CJ[c].median(),4),round(CJ[c].quantile(.025),4),round(CJ[c].quantile(.975),4)] for c in ["feas_baseline","feas_full_data_informed","J_baseline","J_full_data_informed"]}))
# 11. prior ladder: how much of the full-data-informed gain do data-derived centres explain?
B=PRIORS["baseline"]
ladder=[("L0 baseline group",None,{}),
 ("L1 data-derived centres for T and H, baseline dispersion; baseline M, P, S; shared noise; no anchoring",{"full_data_informed":dict(T=("N",.62,.245,.10,.95),H=("NH",9.24,8,96),M=B["M"],P=B["P"],S=B["S"])},dict(equal_noise=True,anchor=False)),
 ("L2 full-data-informed T and H priors (narrow, clipped); baseline M, P, S; shared noise; no anchoring",{"full_data_informed":dict(M=B["M"],P=B["P"],S=B["S"])},dict(equal_noise=True,anchor=False)),
 ("L3 all full-data-informed priors; shared noise; no anchoring",None,dict(equal_noise=True,anchor=False)),
 ("L4 full-data-informed as reported (reduced noise and anchoring)",None,{}),
 ("Reviewer variant: full-data-informed T and H priors; baseline M, P, S; reduced noise and anchoring",{"full_data_informed":dict(M=B["M"],P=B["P"],S=B["S"])},{})]
LD=[]
for name,ov,kw in ladder:
    d=generate(prior_override=ov,**kw); x=d[d.group==("baseline" if name.startswith("L0") else GROUPS[3])]
    LD.append(dict(step=name, feas=x.feasible.mean(), J=x.J.mean(), c_T=x.c_T.mean(), c_G=x.c_G.mean(), c_Mw=x.c_Mw.mean(), c_I=x.c_I.mean(), c_F=x.c_F.mean(), I=x.I.mean(), Mw=x.Mw.mean(), F=x.F.mean()))
LD=pd.DataFrame(LD); LD.to_csv('tables/R2_T18_prior_ladder.csv', index=False); OUT["ladder"]=LD.round(4).to_dict("records")
# 12. hydrodynamic-only group: what distinguishes it from baseline
OUT["hydro_vs_base"]={g:{v:[round(float(d0[d0.group==g][v].mean()),3),round(float(d0[d0.group==g][v].std()),3)] for v in "TPHMS"} for g in GROUPS}
# 13. aortic context check after removing the offset
x=d0[d0.group==GROUPS[3]]; OUT["context_check"]={c:[round(float(x[x.ctx==c].G.mean()),2),round(float(x[x.ctx==c].G.std()),2),round(float(x[x.ctx==c].Gs.iloc[0]),2)] for c in ["pulmonary","aortic"]}
