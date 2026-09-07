"""Exploratory revision analyses, fixed before effect calculation.

Primary comparison uses the intersection of the archived strict orthologues and
the validation cohort's previously defined expression eligibility (6611 genes).
All cohort pairs use identical genes. All 11 existing focused programs are
retained in the coverage table. A composite requires >=3 genes and >=40% of
its original membership, and nonconstant expression in every cohort. This rule
applies to these new full transcriptome comparisons, not the old MERFISH panel.
An explicitly secondary comparison uses all commonly measured uppercase gene
symbols without claiming these are strict one to one orthologues.
Subject resampling holds original gene scaling fixed and jointly resamples
subjects for all programs. Intervals are descriptive conditional percentiles.
"""
from pathlib import Path
import sys,json,itertools,hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr

R=Path(__file__).resolve().parents[2]
O=R/'fbcns_revision/results';O.mkdir(exist_ok=True,parents=True)
sys.path.insert(0,str(R/'code'))
from benchmark_concordance import hedges,pooled
sys.path.insert(0,str(R/'baseline/BBI_v2.0.0/code'))
from run_bbi_extension import LOCKED_MODULES,bh_fdr
COHORTS=['GSE199460','GSE210776','GSE95401','GSE279183','Macnair_All_EC']
B=2000; SEED=20260910

def exact(x,case,ctrl):
    vals=np.r_[x[case],x[ctrl]];n=len(vals);nc=len(case)
    obs=abs(x[case].mean()-x[ctrl].mean()); count=0;total=0
    for ix in itertools.combinations(range(n),nc):
        mask=np.zeros(n,bool);mask[list(ix)]=True
        count+=abs(vals[mask].mean()-vals[~mask].mean())>=obs-1e-12;total+=1
    return count/total

def summarize(g,v):
    ok=np.isfinite(g).all(1)&np.isfinite(v).all(1)
    if ok.sum()<3:return {'n_genes':int(ok.sum())}
    mat=np.column_stack([g[ok],pooled(g[ok,:3],v[ok,:3])])
    c=np.corrcoef(rankdata(mat,axis=0),rowvar=False)
    d={f'{COHORTS[i]}__{COHORTS[j]}':float(c[i,j]) for i,j in itertools.combinations(range(5),2)}
    d['mean_within_mouse']=float(np.mean([c[0,1],c[0,2],c[1,2]]))
    d['mean_mouse_human_first']=float(c[:3,3].mean())
    d['mean_mouse_human_validation']=float(c[:3,4].mean())
    d['mean_mouse_human_all']=float(c[:3,3:5].mean())
    d['within_minus_all_cross']=d['mean_within_mouse']-d['mean_mouse_human_all']
    d['pooled_mouse__GSE279183']=float(c[5,3]);d['pooled_mouse__Macnair_All_EC']=float(c[5,4])
    d['n_genes']=int(ok.sum());return d

def main():
    genes=sorted(pd.read_csv(R/'output/S20_All_EC_gene_effects.csv').gene)
    assert len(genes)==6611
    data=[]
    for name in COHORTS:
        z=np.load(R/'derived'/f'{name}_payload.npz')
        # Retain only the actual case/control comparison for new program scaling.
        keep=np.isin(z['groups'],[str(z['case']),str(z['control'])])
        x=pd.DataFrame(z['expression'][:,keep],index=z['genes']).groupby(level=0).mean()
        gr=z['groups'][keep]; ids=z['sample_ids'][keep]
        ca=np.flatnonzero(gr==z['case']);ct=np.flatnonzero(gr==z['control'])
        data.append(dict(name=name,x=x,case=ca,ctrl=ct,ids=ids,groups=gr))
    gg=[];vv=[]
    for d in data:
        g,v=hedges(d['x'].reindex(genes).to_numpy(),d['case'],d['ctrl']);gg.append(g);vv.append(v)
    g=np.column_stack(gg);v=np.column_stack(vv)
    assert np.isfinite(g).all()
    mat=pd.DataFrame(g,index=genes,columns=COHORTS);mat['pooled_mouse']=pooled(g[:,:3],v[:,:3])
    archived=pd.read_csv(R/'output/S20_All_EC_gene_effects.csv').set_index('gene').reindex(genes)
    errors={'validation':float(np.max(np.abs(mat.Macnair_All_EC-archived.validation_g))),
            'mouse_pooled':float(np.max(np.abs(mat.pooled_mouse-archived.mouse_pooled_g))),
            'first_human':float(np.max(np.abs(mat.GSE279183-archived.human_hedges_g)))}
    assert max(errors.values())<1e-7,errors
    mat.to_csv(O/'S21_five_cohort_gene_effects.csv',index_label='gene')
    point=summarize(g,v);rng=np.random.default_rng(SEED);draws=[]
    for b in range(B):
        bg=[];bv=[]
        for d in data:
            cg,cv=hedges(d['x'].reindex(genes).to_numpy(),rng.choice(d['case'],len(d['case']),replace=True),rng.choice(d['ctrl'],len(d['ctrl']),replace=True))
            bg.append(cg);bv.append(cv)
        draws.append(summarize(np.column_stack(bg),np.column_stack(bv)))
    bt=pd.DataFrame(draws);bt.to_csv(O/'S21_subject_bootstrap_draws.csv',index=False)
    summary=[]
    for k,p in point.items():
        if k=='n_genes':continue
        a=bt[k].dropna();summary.append(dict(comparison=k,estimate=p,range_low=a.quantile(.025),range_high=a.quantile(.975),genes=len(genes),draws_attempted=B,draws_estimable=len(a)))
    pd.DataFrame(summary).to_csv(O/'S21_five_cohort_concordance.csv',index=False)
    loo=[]
    for j,d in enumerate(data):
        for s in np.r_[d['case'],d['ctrl']]:
            aa=g.copy();av=v.copy();aa[:,j],av[:,j]=hedges(d['x'].reindex(genes).to_numpy(),d['case'][d['case']!=s],d['ctrl'][d['ctrl']!=s])
            loo.append(dict(omitted_cohort=d['name'],omitted_subject=d['ids'][s],**summarize(aa,av)))
    pd.DataFrame(loo).to_csv(O/'S21_leave_one_subject_out.csv',index=False)
    # Common measured genes, excluding nonvarying genes in any comparison.
    common=set.intersection(*(set(d['x'].index) for d in data))
    nonconstant=set.intersection(*(set(d['x'].index[d['x'].std(axis=1,ddof=1)>1e-12]) for d in data))
    allcoverage=[];alltests=[];allscores=[];pairrows=[]
    for mode,universe in [('strict_6611',set(genes)),('common_symbol_sensitivity',common)]:
        panels={}
        for mod,members in LOCKED_MODULES.items():
            used=[a for a in members if a in universe and a in nonconstant]
            eligible=len(used)>=3 and len(used)/len(members)>=.4
            allcoverage.append(dict(matching=mode,program=mod,original_gene_count=len(members),common_gene_count=len(used),coverage=len(used)/len(members),eligible_composite=eligible,used_genes='|'.join(used),excluded_genes='|'.join(a for a in members if a not in used),criterion='at least 3 nonconstant genes and 40 percent coverage in every cohort'))
            if eligible:panels[mod]=used
        modnames=list(panels); effects=[]; mode_rows=[]
        for d in data:
            scores=[]
            for mod,gs in panels.items():
                xx=d['x'].loc[gs].to_numpy();zz=(xx-xx.mean(axis=1,keepdims=True))/xx.std(axis=1,ddof=1,keepdims=True)
                sc=zz.mean(axis=0);scores.append(sc)
                for sid,group,s in zip(d['ids'],d['groups'],sc):allscores.append(dict(matching=mode,cohort=d['name'],program=mod,subject=sid,group=group,score=float(s)))
            ss=np.array(scores);eg,ev=hedges(ss,d['case'],d['ctrl']);effects.append(eg)
            brng=np.random.default_rng(SEED+1);bootstrap=[]
            for _ in range(B):
                bca=brng.choice(d['case'],len(d['case']),replace=True);bct=brng.choice(d['ctrl'],len(d['ctrl']),replace=True)
                bootstrap.append(hedges(ss,bca,bct)[0])
            bs=np.array(bootstrap)
            for i,mod in enumerate(modnames):
                bs0=bs[:,i];valid=bs0[np.isfinite(bs0)]
                mode_rows.append(dict(matching=mode,cohort=d['name'],program=mod,genes=len(panels[mod]),n_case=len(d['case']),n_control=len(d['ctrl']),hedges_g=float(eg[i]),variance=float(ev[i]),exact_p=exact(ss[i],d['case'],d['ctrl']),range_low=float(np.quantile(valid,.025)),range_high=float(np.quantile(valid,.975)),draws_attempted=B,draws_estimable=len(valid)))
        adj=bh_fdr([x['exact_p'] for x in mode_rows])
        for a,q in zip(mode_rows,adj):a['fdr_all_cohorts_and_programs']=float(q);a['fdr_family_size']=len(mode_rows)
        alltests+=mode_rows
        em=np.array(effects).T
        for j,k in itertools.combinations(range(5),2):
            pairrows.append(dict(matching=mode,cohort_1=COHORTS[j],cohort_2=COHORTS[k],programs=len(modnames),direction_agreement=float(np.mean(np.sign(em[:,j])==np.sign(em[:,k]))),note='descriptive across overlapping programs, not independent replication'))
    pd.DataFrame(allcoverage).to_csv(O/'S22_function_gene_coverage.csv',index=False)
    pd.DataFrame(alltests).to_csv(O/'S22_function_effects.csv',index=False)
    pd.DataFrame(allscores).to_csv(O/'S22_subject_function_scores.csv',index=False)
    pd.DataFrame(pairrows).to_csv(O/'S22_function_direction_agreement.csv',index=False)
    audit=dict(seed=SEED,program_seed=SEED+1,genes=len(genes),draws_attempted=B,draws_estimable=int(bt.mean_within_mouse.notna().sum()),resample_gene_min=int(bt.loc[bt.mean_within_mouse.notna(),'n_genes'].min()),resample_gene_median=float(bt.loc[bt.mean_within_mouse.notna(),'n_genes'].median()),resample_gene_max=int(bt.loc[bt.mean_within_mouse.notna(),'n_genes'].max()),reconstruction_errors=errors,cohorts=[dict(cohort=d['name'],cases=len(d['case']),controls=len(d['ctrl']),subjects=list(d['ids'])) for d in data],point=point,
       composition_scope='Published endothelial subtype counts and donor regional scores are retained. No region by cell type count matrix was available in the retained inputs for estimating nonendothelial contributions to GSE279183 spatial pathway scores. No adjusted cell type or causal pathway claim is made.',
       interpretation='Exploratory revision. Functions overlap in gene membership. Bootstrap ranges are conditional descriptive percentiles, not calibrated confidence intervals. The human cohort pair is one pair, not a population estimate of human reproducibility. Strict orthologues retain identical uppercase symbols under the archived mapping rule.')
    (O/'ANALYSIS_AUDIT.json').write_text(json.dumps(audit,indent=2))
    print(pd.DataFrame(summary).to_string(index=False))
    print(pd.DataFrame(alltests).query("matching=='strict_6611'").pivot(index='program',columns='cohort',values='hedges_g').round(3).to_string())
    print('AUDIT',json.dumps({k:v for k,v in audit.items() if k not in ['point','cohorts']},indent=2))

if __name__=='__main__':main()
