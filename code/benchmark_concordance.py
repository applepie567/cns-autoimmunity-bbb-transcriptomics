"""Matched gene benchmark with stratified subject resampling.

All calculations use archived strict orthologues. Cohorts and gene sets are
fixed before the resampling. No gene bootstrap is used for subject intervals.
"""
from pathlib import Path
import json, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output';OUT.mkdir(exist_ok=True)
SEED=20260907;B=2000
COHORTS=['GSE199460','GSE210776','GSE95401','GSE279183']
strict=pd.read_csv(ROOT/'baseline/BBI_v2.0.0/results/strict_orthology_gene_effects.csv').set_index('gene').sort_index()
GENES=strict.index.tolist()

def hedges(x,case,ctrl):
    a=x[:,case];b=x[:,ctrl];n1=a.shape[1];n0=b.shape[1];df=n1+n0-2
    ma=a.mean(1);mb=b.mean(1)
    s2=(((a-ma[:,None])**2).sum(1)+((b-mb[:,None])**2).sum(1))/df
    # Repeated identical values can have nonzero floating-point residuals.
    scale=np.maximum(np.mean(np.column_stack([a,b])**2,axis=1),1.0)
    s2[s2<=np.finfo(float).eps*scale]=np.nan
    with np.errstate(divide='ignore',invalid='ignore'):
        g=(1-3/(4*df-1))*(ma-mb)/np.sqrt(s2)
    g[~np.isfinite(g)]=np.nan
    v=(n1+n0)/(n1*n0)+g*g/(2*df)
    return g,v

def pooled(g,v):
    w=1/v;fixed=(w*g).sum(1)/w.sum(1)
    q=(w*(g-fixed[:,None])**2).sum(1)
    c=w.sum(1)-(w*w).sum(1)/w.sum(1)
    tau=np.maximum(0,np.divide(q-(g.shape[1]-1),c,out=np.zeros_like(c),where=c>0))
    wr=1/(v+tau[:,None]);return (wr*g).sum(1)/wr.sum(1)

def summarize(g,v):
    # Every comparison in a draw uses exactly the same estimable genes.
    mask=np.isfinite(g).all(1)&np.isfinite(v).all(1)
    if mask.sum()<3:
        keys=[f'{COHORTS[i]}__{COHORTS[j]}' for i,j in itertools.combinations(range(4),2)]
        keys+=['pooled_mouse__GSE279183','mean_within_mouse','mean_individual_mouse_human','within_minus_individual_cross','within_minus_pooled_cross']
        return dict(**{k:np.nan for k in keys},n_genes=int(mask.sum()))
    a=g[mask];vv=v[mask];mp=pooled(a[:,:3],vv[:,:3])
    a=np.column_stack([a,mp])
    r=rankdata(a,axis=0);corr=np.corrcoef(r,rowvar=False)
    vals={f'{COHORTS[i]}__{COHORTS[j]}':float(corr[i,j]) for i,j in itertools.combinations(range(4),2)}
    vals['pooled_mouse__GSE279183']=float(corr[4,3])
    vals['mean_within_mouse']=float(np.mean([corr[0,1],corr[0,2],corr[1,2]]))
    vals['mean_individual_mouse_human']=float(np.mean(corr[:3,3]))
    vals['within_minus_individual_cross']=vals['mean_within_mouse']-vals['mean_individual_mouse_human']
    vals['within_minus_pooled_cross']=vals['mean_within_mouse']-vals['pooled_mouse__GSE279183']
    vals['n_genes']=int(mask.sum())
    return vals

def main():
    payload=[];gr=[];vr=[]
    for cohort in COHORTS:
        z=np.load(ROOT/'derived'/f'{cohort}_payload.npz')
        x=pd.DataFrame(z['expression'],index=z['genes']).groupby(level=0).mean().reindex(GENES).to_numpy()
        groups=z['groups'];case=np.flatnonzero(groups==z['case']);ctrl=np.flatnonzero(groups==z['control'])
        g,v=hedges(x,case,ctrl);gr.append(g);vr.append(v)
        payload.append((x,case,ctrl,z['sample_ids']))
        print(cohort,len(case),len(ctrl),'missing',np.isnan(g).sum(),flush=True)
    g=np.column_stack(gr);v=np.column_stack(vr)
    ref=pd.read_csv(ROOT/'baseline/BBI_v2.0.0/results/acute_EAE_gene_effects_by_dataset.csv').pivot(index='gene',columns='dataset',values='hedges_g').reindex(GENES)[COHORTS[:3]].to_numpy()
    errors={'mouse_max_abs_error':float(np.nanmax(np.abs(g[:,:3]-ref))),
            'human_max_abs_error':float(np.nanmax(np.abs(g[:,3]-strict.human_hedges_g.to_numpy()))),
            'pooled_max_abs_error':float(np.nanmax(np.abs(pooled(g[:,:3],v[:,:3])-strict.mouse_pooled_g.to_numpy())))}
    assert max(errors.values())<1e-7,errors
    print('RECONSTRUCTION',errors,flush=True)
    point=summarize(g,v)
    mat=pd.DataFrame(g,index=GENES,columns=COHORTS)
    mat['pooled_mouse']=pooled(g[:,:3],v[:,:3]);mat.to_csv(OUT/'S19_matched_gene_effects.csv',index_label='gene')
    rng=np.random.default_rng(SEED);boot=[]
    for k in range(B):
        gs=[];vs=[]
        for x,case,ctrl,ids in payload:
            cg,cv=hedges(x,rng.choice(case,len(case),replace=True),rng.choice(ctrl,len(ctrl),replace=True))
            gs.append(cg);vs.append(cv)
        boot.append(summarize(np.column_stack(gs),np.column_stack(vs)))
        if (k+1)%500==0:print('bootstrap',k+1,flush=True)
    bt=pd.DataFrame(boot);bt.to_csv(OUT/'subject_bootstrap_draws.csv',index=False)
    rows=[]
    for key,value in point.items():
        if key=='n_genes':continue
        rows.append(dict(comparison=key,rho_or_difference=value,ci_low=float(bt[key].quantile(.025)),ci_high=float(bt[key].quantile(.975)),n_genes=point['n_genes'],interval='stratified animal and donor percentile bootstrap conditional on estimable draws',draws_attempted=B,draws_estimable=int(bt[key].notna().sum())))
    summary=pd.DataFrame(rows);summary.to_csv(OUT/'S19_concordance_benchmark.csv',index=False)
    jack=[]
    for j,(x,case,ctrl,ids) in enumerate(payload):
        for s in np.r_[case,ctrl]:
            gg=g.copy();vv=v.copy()
            gg[:,j],vv[:,j]=hedges(x,case[case!=s],ctrl[ctrl!=s])
            jack.append(dict(omitted_cohort=COHORTS[j],omitted_subject=ids[s],**summarize(gg,vv)))
    pd.DataFrame(jack).to_csv(OUT/'S19_leave_one_subject_out.csv',index=False)
    valid=bt[bt.mean_within_mouse.notna()]
    info=dict(seed=SEED,draws_attempted=B,draws_estimable=len(valid),draws_unestimable=B-len(valid),reconstruction=errors,point=point,
              bootstrap_gene_count_min=int(valid.n_genes.min()),bootstrap_gene_count_median=float(valid.n_genes.median()),bootstrap_gene_count_max=int(valid.n_genes.max()),
              explanation='Undefined standardized effects caused by zero variance after resampling are excluded jointly from every comparison within the same draw. Draws with fewer than 3 estimable genes are omitted from percentile calculations, without replacing or hiding them. The original target set is fixed at 7705 orthologues. These descriptive intervals condition on the selected cohorts and estimable resamples. They do not have guaranteed nominal coverage with two cases in one cohort.')
    (OUT/'benchmark_details.json').write_text(json.dumps(info,indent=2))
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
