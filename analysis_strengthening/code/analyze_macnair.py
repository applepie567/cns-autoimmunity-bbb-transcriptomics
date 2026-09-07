"""Lesion-matched and capillary-restricted validation, without cell-level inference."""
from pathlib import Path
import sys,json,itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1];OUT=R/'output'
sys.path.insert(0,str(R/'baseline/BBI_v2.0.0/code'))
import run_bbi_extension as old
from benchmark_concordance import hedges
B=2000;SEED=20260908

def exact(x,case,ctrl):
    obs=float(x[case].mean()-x[ctrl].mean());allix=np.r_[case,ctrl];extreme=0;n=0
    for inds in itertools.combinations(allix,len(case)):
        a=np.array(inds);b=np.setdiff1d(allix,a)
        extreme+=abs(x[a].mean()-x[b].mean())>=abs(obs)-1e-12;n+=1
    return obs,extreme/n

def subset_counts(counts,m,subtype,minimum,bank):
    keep=m.lesion_type.isin(['CAL','WM'])
    if subtype:keep &= m.type_fine.eq(subtype)
    if bank:keep &= m.sample_source.eq(bank)
    records=[];arrays=[]
    for (donor,lesion),block in m[keep].groupby(['individual_id_anon','lesion_type'],sort=True):
        if block.nuclei.sum()<minimum:continue
        arrays.append(counts[:,block.index].sum(1))
        records.append(dict(donor=donor,lesion=lesion,nuclei=int(block.nuclei.sum()),bank=block.sample_source.iloc[0],age=block.age.iloc[0],sex=block.sex.iloc[0],pmi=block.pmi.iloc[0]))
    return np.column_stack(arrays),pd.DataFrame(records)

def main():
    z=np.load(R/'derived/Macnair_WM_endothelial_counts.npz');m=pd.read_csv(R/'derived/Macnair_pseudobulk_groups.csv')
    counts=pd.DataFrame(z['counts'],index=z['genes']).groupby(level=0,sort=True).sum()
    genes=counts.index.tolist();lookup={g:i for i,g in enumerate(genes)}
    ref=pd.read_csv(R/'baseline/BBI_v2.0.0/results/strict_orthology_gene_effects.csv')[['gene','mouse_pooled_g','human_hedges_g']]
    rng=np.random.default_rng(SEED);module_rows=[];cor_rows=[]
    main_counts,_=subset_counts(counts.to_numpy(),m,None,20,None)
    main_log=old.log_cpm(main_counts)
    center=main_log.mean(1);scale=main_log.std(1,ddof=1);scale[scale==0]=1
    _,capmeta=subset_counts(counts.to_numpy(),m,'Endocyte_Capillary',20,None)
    capkeys=set(zip(capmeta.donor,capmeta.lesion))
    specs=[('All_EC',None,20,None),('Capillary','Endocyte_Capillary',20,None),('All_EC_capillary_donors',None,20,None),('All_EC_same_bank',None,20,'Amsterdam BB'),('All_EC_min10',None,10,None),('Capillary_min10','Endocyte_Capillary',10,None)]
    for label,subtype,minimum,bank in specs:
        c,meta=subset_counts(counts.to_numpy(),m,subtype,minimum,bank)
        if label=='All_EC_capillary_donors':
            keep=np.array([(d,l) in capkeys for d,l in zip(meta.donor,meta.lesion)])
            c=c[:,keep];meta=meta.loc[keep].reset_index(drop=True)
        if label=='All_EC':primary=meta.copy()
        case=np.flatnonzero(meta.lesion.eq('CAL'));ctrl=np.flatnonzero(meta.lesion.eq('WM'))
        assert min(len(case),len(ctrl))>=2
        print(label,len(case),len(ctrl),int(meta.nuclei.sum()),flush=True)
        log=old.log_cpm(c);g,v=hedges(log,case,ctrl)
        cpm=c/np.maximum(c.sum(0),1)[None,:]*1e6
        dge=pd.DataFrame(dict(gene=genes,validation_g=g,validation_var=v,expression_eligible=(cpm>=1).sum(1)>=3)).merge(ref,on='gene',how='inner')
        dge=dge[dge.expression_eligible & np.isfinite(dge.validation_g)].copy()
        dge.to_csv(OUT/f'S20_{label}_gene_effects.csv',index=False)
        # Preserve gene weights when restricting subtypes or changing donors.
        zs=(log-center[:,None])/scale[:,None];st=meta.copy();mr=[]
        for mod,gs in old.STATE_MODULES.items():
            ix=[lookup[g] for g in gs if g in lookup];s=zs[ix].mean(0);st[mod]=s
            mg,mv=hedges(s[None,:],case,ctrl);delta,p=exact(s,case,ctrl);draw=[];loo=[]
            for k in range(B):
                gg,vv=hedges(s[None,:],rng.choice(case,len(case),replace=True),rng.choice(ctrl,len(ctrl),replace=True));draw.append(gg[0])
            a=np.array(draw);a=a[np.isfinite(a)]
            for i in np.r_[case,ctrl]:
                gg,vv=hedges(s[None,:],case[case!=i],ctrl[ctrl!=i]);loo.append(gg[0])
            mr.append(dict(subset=label,module=mod,n_case=len(case),n_control=len(ctrl),n_genes=len(ix),mean_difference=delta,hedges_g=float(mg[0]),p_exact=p,bootstrap_low=float(np.quantile(a,.025)),bootstrap_high=float(np.quantile(a,.975)),bootstrap_valid=len(a),loo_g_min=float(np.nanmin(loo)),loo_g_max=float(np.nanmax(loo))))
        q=old.bh_fdr([r['p_exact'] for r in mr])
        for row,fdr in zip(mr,q):row['FDR_BH']=float(fdr)
        module_rows.extend(mr);st.to_csv(OUT/f'S20_{label}_sample_scores.csv',index=False)
        ix=np.array([lookup[g] for g in dge.gene]);xb=log[ix];references=dge[['mouse_pooled_g','human_hedges_g']].to_numpy();br=[]
        for k in range(B):
            gg,vv=hedges(xb,rng.choice(case,len(case),replace=True),rng.choice(ctrl,len(ctrl),replace=True));ok=np.isfinite(gg)
            br.append([spearmanr(gg[ok],references[ok,j]).statistic for j in range(2)])
        br=np.array(br)
        for j,target in enumerate(['pooled_mouse','GSE279183_human']):
            cor_rows.append(dict(subset=label,reference=target,n_genes=len(dge),rho=float(spearmanr(dge.validation_g,references[:,j]).statistic),same_direction=float((np.sign(dge.validation_g)==np.sign(references[:,j])).mean()),conditional_bootstrap_low=float(np.nanquantile(br[:,j],.025)),conditional_bootstrap_high=float(np.nanquantile(br[:,j],.975)),conditional_bootstrap_valid=int(np.isfinite(br[:,j]).sum()),interval_scope='validation donors only, reference effects fixed',n_case=len(case),n_control=len(ctrl)))
        np.savez_compressed(R/'derived'/f'Macnair_{label}_payload.npz',genes=np.array(genes),expression=log,sample_ids=meta.donor.to_numpy(dtype=str),groups=meta.lesion.to_numpy(dtype=str),case='CAL',control='WM')
    pd.DataFrame(module_rows).to_csv(OUT/'S20_validation_program_effects.csv',index=False)
    pd.DataFrame(cor_rows).to_csv(OUT/'S20_validation_concordance.csv',index=False)
    allg=pd.read_csv(OUT/'S20_All_EC_gene_effects.csv').set_index('gene');capg=pd.read_csv(OUT/'S20_Capillary_gene_effects.csv').set_index('gene');matchg=pd.read_csv(OUT/'S20_All_EC_capillary_donors_gene_effects.csv').set_index('gene');common=allg.index.intersection(capg.index).intersection(matchg.index);cr=[]
    for label,t in [('All_EC',allg),('All_EC_capillary_donors',matchg),('Capillary',capg)]:
        for target in ['mouse_pooled_g','human_hedges_g']:
            cr.append(dict(subset=label,reference=target,n_genes=len(common),rho=float(spearmanr(t.loc[common,'validation_g'],t.loc[common,target]).statistic),same_direction=float((np.sign(t.loc[common,'validation_g'])==np.sign(t.loc[common,target])).mean())))
    pd.DataFrame(cr).to_csv(OUT/'S20_common_gene_subtype_sensitivity.csv',index=False)
    comp=m[m.lesion_type.isin(['CAL','WM'])].pivot_table(index=['individual_id_anon','lesion_type'],columns='type_fine',values='nuclei',aggfunc='sum',fill_value=0)
    idx=pd.MultiIndex.from_frame(primary[['donor','lesion']].rename(columns={'donor':'individual_id_anon','lesion':'lesion_type'}));comp=comp.reindex(idx)
    props=comp.div(comp.sum(1),axis=0);props.reset_index().to_csv(OUT/'S20_endothelial_subtype_proportions.csv',index=False)
    case=np.flatnonzero(primary.lesion.eq('CAL'));ctrl=np.flatnonzero(primary.lesion.eq('WM'));rows=[]
    for s in props:
        x=props[s].to_numpy();delta,p=exact(x,case,ctrl);rows.append(dict(subtype=s,mean_CAL=float(x[case].mean()),mean_control=float(x[ctrl].mean()),mean_difference=delta,p_exact=p,n_case=len(case),n_control=len(ctrl)))
    df=pd.DataFrame(rows);df['FDR_BH']=old.bh_fdr(df.p_exact);df.to_csv(OUT/'S20_composition_comparisons.csv',index=False)
    info=dict(seed=SEED,draws=B,minimum_nuclei_primary=20,minimum_nuclei_sensitivity=10,expression_filter='CPM at least 1 in at least 3 eligible donors',primary_tissue_banks=primary.bank.unique().tolist(),original_human_tissue_bank='UK Multiple Sclerosis Tissue Bank',primary='CAL versus healthy white matter, cohort I',caution='Primary capillary eligibility changes donors, so an additional all EC comparison uses exactly those same donors. Healthy endothelial samples are small. Reference effects remain fixed in conditional donor bootstrap.')
    (OUT/'validation_details.json').write_text(json.dumps(info,indent=2))
    print(pd.DataFrame(module_rows).to_string(index=False));print(pd.DataFrame(cor_rows).to_string(index=False));print(df.to_string(index=False))
if __name__=='__main__':main()
