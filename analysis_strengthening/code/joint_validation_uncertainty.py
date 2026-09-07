"""Propagate uncertainty from the mouse cohorts and both human datasets."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from benchmark_concordance import hedges,pooled,COHORTS
R=Path(__file__).resolve().parents[1];O=R/'output';B=2000;SEED=20260909
def load_payload(path,genes):
    z=np.load(path);x=pd.DataFrame(z['expression'],index=z['genes']).groupby(level=0).mean().reindex(genes).to_numpy()
    return x,np.flatnonzero(z['groups']==z['case']),np.flatnonzero(z['groups']==z['control'])
def main():
    summary=pd.read_csv(O/'S20_validation_concordance.csv');out=[]
    for label in summary.subset.drop_duplicates():
        dge=pd.read_csv(O/f'S20_{label}_gene_effects.csv');genes=dge.gene.tolist()
        ps=[load_payload(R/'derived'/f'{c}_payload.npz',genes) for c in COHORTS]
        ps.append(load_payload(R/'derived'/f'Macnair_{label}_payload.npz',genes))
        rng=np.random.default_rng(SEED);draws=[]
        for k in range(B):
            gg=[];vv=[]
            for x,case,ctrl in ps:
                g,v=hedges(x,rng.choice(case,len(case),replace=True),rng.choice(ctrl,len(ctrl),replace=True));gg.append(g);vv.append(v)
            g=np.column_stack(gg);v=np.column_stack(vv);ok=np.isfinite(g).all(1)&np.isfinite(v).all(1)
            if ok.sum()<3:draws.append([np.nan,np.nan,int(ok.sum())]);continue
            g=g[ok];v=v[ok];mp=pooled(g[:,:3],v[:,:3])
            draws.append([spearmanr(g[:,4],mp).statistic,spearmanr(g[:,4],g[:,3]).statistic,int(ok.sum())])
        dr=np.array(draws);valid=np.isfinite(dr[:,0])&np.isfinite(dr[:,1])
        for j,target in enumerate(['pooled_mouse','GSE279183_human']):
            row=summary[(summary.subset==label)&(summary.reference==target)].iloc[0].to_dict()
            row={k:v for k,v in row.items() if not k.startswith('conditional_')}
            row.update(subject_bootstrap_low=float(np.quantile(dr[valid,j],.025)),subject_bootstrap_high=float(np.quantile(dr[valid,j],.975)),draws_attempted=B,draws_estimable=int(valid.sum()),bootstrap_genes_min=int(dr[valid,2].min()),interval_scope='joint animal and donor resampling of both study sides, conditional on estimable draws')
            out.append(row)
        pd.DataFrame(dr,columns=['rho_pooled_mouse','rho_original_human','n_genes']).to_csv(O/f'S20_{label}_joint_bootstrap_draws.csv',index=False)
        print(label,'valid',valid.sum(),'ranges',np.quantile(dr[valid,:2],[.025,.975],axis=0),flush=True)
    pd.DataFrame(out).to_csv(O/'S20_validation_concordance.csv',index=False)
    p=O/'validation_details.json';d=json.loads(p.read_text());d.update(joint_bootstrap_seed=SEED,joint_bootstrap_draws=B,caution='Small healthy donor samples and changing donor eligibility limit precision. Final concordance intervals jointly resample both reference and validation subjects.',interval_scope='All reported concordance intervals propagate subject uncertainty from both sides of the comparison. Degenerate draws are counted and excluded from descriptive percentile ranges.');p.write_text(json.dumps(d,indent=2))
if __name__=='__main__':main()
