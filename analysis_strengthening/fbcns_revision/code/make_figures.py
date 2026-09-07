from pathlib import Path
import sys,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib.text import Text
R=Path(__file__).resolve().parents[2];D=R/'fbcns_revision';F=D/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
N=['GSE199460','GSE210776','GSE95401','GSE279183','Macnair_All_EC']
L=['Mouse brain','Mouse spinal 1','Mouse spinal 2','Human 1','Human 2']
P={'Tight_junction':'Tight junctions','BBB_transport_identity':'BBB transport','Caveolae_structural':'Caveolar structure','PLVAP_permeability':'PLVAP associated genes','Adhesion_trafficking':'Adhesion and trafficking','IFN_antigen_presentation':'Interferon and antigen presentation','Wnt_BBB':'Wnt and BBB properties','TGF_response':'TGF response','VEGF_response':'VEGF response','ROS_Src':'Oxidative stress and Src','ECM_protease':'Matrix and proteases'}
def save(fig,n):
    for t in fig.findobj(match=Text):t.set_color('black')
    for ext in ['png','pdf']:fig.savefig(F/(n+'.'+ext),dpi=350,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def workflow():
    fig,ax=plt.subplots(figsize=(8,6));ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
    def box(x,y,w,h,title,body):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012',facecolor='#f3f5f7',edgecolor='#61727d'))
        ax.text(x+.014,y+h-.025,title,fontweight='bold',va='top',fontsize=10)
        ax.text(x+.014,y+h-.09,body,va='top',fontsize=9,linespacing=1.5)
    box(.02,.68,.45,.30,'Three acute EAE cohorts','GSE199460   Brain   2 EAE and 3 controls\nGSE210776   Spinal cord   4 and 3\nGSE95401   Spinal cord   3 and 3\nOne expression profile per animal')
    box(.53,.68,.45,.30,'Two chronic active MS cohorts','Human 1   GSE279183\n4 MS and 6 control donors\nHuman 2   Macnair cohort I   7 and 4\nEndothelial nuclei combined by donor\nChronic active lesions and control WM')
    box(.12,.355,.76,.245,'Comparable effects across five cohorts','6,611 archived strict orthologues and the same effect estimator\nAll ten cohort pairs with animal and donor resampling\nEleven focused functions with explicit gene coverage')
    for x in [.24,.76]:ax.add_patch(FancyArrowPatch((x,.66),(.5,.62),arrowstyle='-|>',mutation_scale=13,color='#61727d'))
    box(.02,.035,.45,.245,'Human endothelial checks','Source annotated endothelial subtypes\nCapillary restriction in the same donors\nSensitivity to tissue bank and cell yield')
    box(.53,.035,.45,.245,'Spatial and functional context','Local regions in GSE279183\nand GSE284005\nLesion categories in GSE208747\nSeparate published barrier measurements')
    save(fig,'Figure_1')
def function_heatmap(ax,mode,title,fig):
    d=pd.read_csv(D/'results/S22_function_effects.csv');c=pd.read_csv(D/'results/S22_function_gene_coverage.csv');c=c[c.matching.eq(mode)].set_index('program')
    d=d[d.matching.eq(mode)].pivot(index='program',columns='cohort',values='hedges_g').reindex(index=P,columns=N)
    # Pale diverging colors keep every numeric label black and readable.
    from matplotlib.colors import LinearSegmentedColormap
    cm=LinearSegmentedColormap.from_list('pale',['#9dbbd0','#ffffff','#e1afa0']);cm.set_bad('#eeeeee')
    im=ax.imshow(d.to_numpy(),cmap=cm,vmin=-5,vmax=5,aspect='auto')
    for i,m in enumerate(P):
        for j in range(5):
            val=d.iloc[i,j];ax.text(j,i,f'{val:.2f}' if np.isfinite(val) else 'Not scored',ha='center',va='center',fontsize=8)
    labs=[f'{P[m]}   ({int(c.loc[m,"common_gene_count"])} of {int(c.loc[m,"original_gene_count"])})' for m in P]
    ax.set_yticks(range(len(P)),labs,fontsize=8.2);ax.set_xticks(range(5),L,fontsize=8)
    ax.set_title(title,loc='left',fontweight='bold',pad=14,fontsize=10)
    ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,11,1),minor=True);ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',length=0)
    cb=fig.colorbar(im,ax=ax,fraction=.023,pad=.025,extend='both');cb.set_label("Hedges' g",fontsize=8)
def figure2():
    e=pd.read_csv(R/'baseline/BBI_v2.0.0/results/dataset_module_effects.csv');e=e[e.contrast.str.startswith('Acute')]
    meta=pd.read_csv(R/'baseline/BBI_v2.0.0/results/acute_EAE_random_effects_meta.csv').set_index('module')
    fig=plt.figure(figsize=(9.5,8.7));gs=fig.add_gridspec(2,3,height_ratios=[1,1.95],hspace=.47,wspace=.60)
    for i,(m,label) in enumerate([('Endothelial_immune_activation','Immune activation'),('BBB_specialization','BBB specialization'),('Structural_ECM_remodeling','Structural and ECM')]):
        a=fig.add_subplot(gs[0,i]);p=meta.loc[m]
        for j,ds in enumerate(['GSE199460_all','GSE210776','GSE95401']):
            r=e[(e.dataset.eq(ds))&e.module.eq(m)].iloc[0]
            a.plot([r.ci_low,r.ci_high],[j,j],color='#607e8d',lw=1.2);a.plot(r.hedges_g,j,'o',color='#37647c',ms=4)
        a.plot([p.mKH_ci_low,p.mKH_ci_high],[3,3],ls=':',color='#37647c',lw=1.5)
        a.plot([p.ci_low,p.ci_high],[3,3],color='#37647c',lw=2);a.plot(p.pooled_g,3,'D',color='#37647c',ms=5)
        a.axvline(0,lw=.7,color='.6');a.set_yticks(range(4),['Brain','Spinal 1','Spinal 2','Pooled']);a.set_ylim(3.6,-.5)
        a.set_xlabel("Hedges' g",fontsize=9);a.set_title(chr(65+i)+'  '+label,loc='left',fontweight='bold',fontsize=9.7,pad=12)
    a=fig.add_subplot(gs[1,:]);function_heatmap(a,'strict_6611','D  Same genes for each function across all five cohorts',fig)
    fig.subplots_adjust(left=.32,right=.97,top=.94,bottom=.065)
    save(fig,'Figure_2')
    fig,ax=plt.subplots(figsize=(9.5,5.4));function_heatmap(ax,'common_symbol_sensitivity','Common gene symbols as a sensitivity analysis',fig);fig.subplots_adjust(left=.39,right=.96,bottom=.12,top=.90);save(fig,'Figure_S10')
def figure4():
    d=pd.read_csv(D/'results/S21_five_cohort_gene_effects.csv').set_index('gene');c=d[N].corr(method='spearman')
    s=pd.read_csv(D/'results/S21_five_cohort_concordance.csv').set_index('comparison')
    fig=plt.figure(figsize=(8.5,9.4));gs=fig.add_gridspec(2,2,height_ratios=[1,1.25],wspace=.45,hspace=.55)
    a=fig.add_subplot(gs[0,0]);a.imshow(c,vmin=-.1,vmax=1,cmap='Blues')
    for i in range(5):
        for j in range(5):a.text(j,i,f'{c.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.7))
    a.set_xticks(range(5),['Brain','Spinal 1','Spinal 2','Human 1','Human 2'],rotation=50,ha='right',fontsize=8);a.set_yticks(range(5),L,fontsize=8)
    a.set_title('A  Same 6,611 genes',loc='left',fontweight='bold',fontsize=10,pad=12)
    b=fig.add_subplot(gs[0,1]);keys=['mean_within_mouse','mean_mouse_human_first','mean_mouse_human_validation','mean_mouse_human_all','within_minus_all_cross']
    labs=['Mean within mice','Mean mice with human 1','Mean mice with human 2','Mean all six human and mouse pairs','Within mice minus all six pairs']
    for i,k in enumerate(keys):
        r=s.loc[k];b.plot([r.range_low,r.range_high],[i,i],color='#476b80',lw=1.8);b.plot(r.estimate,i,'o',ms=4,color='#476b80');b.text(-.14,i-.30,labs[i],fontsize=7.7,va='center')
    b.axvline(0,color='.6',lw=.7);b.set(xlim=(-.15,.50),ylim=(4.6,-.6));b.set_yticks([]);b.set_xlabel('Correlation or difference');b.set_title('B  Descriptive subject ranges',loc='left',fontweight='bold',fontsize=10,pad=12)
    a=fig.add_subplot(gs[1,:]);pairs=[(i,j) for i in range(3) for j in range(i+1,3)]+[(i,j) for j in [3,4] for i in range(3)]+[(3,4)]
    for y,(i,j) in enumerate(pairs):
        r=s.loc[N[i]+'__'+N[j]];col='#38647b' if j<3 else '#966549';a.plot([r.range_low,r.range_high],[y,y],color=col,lw=1.5);a.plot(r.estimate,y,'o',ms=4,color=col)
        a.text(.62,y,f'{r.estimate:.2f}',va='center',ha='right',fontsize=8)
    a.set_yticks(range(10),[L[i]+' with '+L[j].lower() for i,j in pairs],fontsize=8);a.set_ylim(9.7,-.7);a.set_xlim(-.2,.65);a.axvline(0,color='.6',lw=.7);a.set_xlabel("Spearman's correlation");a.set_title('C  All ten cohort pairs',loc='left',fontweight='bold',pad=14,fontsize=10)
    fig.subplots_adjust(left=.28,right=.97,top=.95,bottom=.065);save(fig,'Figure_4')
def supplementary_black_text():
    import importlib.util
    from matplotlib.figure import Figure
    from matplotlib.colors import to_rgba
    out=F/'supplementary';out.mkdir(exist_ok=True)
    original=Figure.savefig
    def wrapped(fig,*args,**kwargs):
        for t in fig.findobj(match=Text):
            try:white=to_rgba(t.get_color())[:3]==(1.,1.,1.)
            except Exception:white=False
            if white:t.set_bbox(dict(facecolor='white',alpha=.85,edgecolor='none',pad=.5))
            t.set_color('black')
            if 'Independent groups;' in t.get_text():t.set_text('Independent groups. Three biological samples per group. Bars show means.')
        return original(fig,*args,**kwargs)
    Figure.savefig=wrapped
    try:
        for fn in ['make_bbi_figures','make_orthogonal_supplementary_figure','make_s9_from_frozen_results']:
            sp=importlib.util.spec_from_file_location(fn,R/'baseline/BBI_v2.0.0/code'/f'{fn}.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
            m.OUT=out;m.SUPP_OUT=out;m.FIGURES=out
            if fn=='make_bbi_figures':
                m.set_style()
                def suppsave(fig,stem,out_dir=None):
                    fig.savefig(out/(stem+'.png'),dpi=350,bbox_inches='tight',facecolor='white')
                    fig.savefig(out/(stem+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
                m.save=suppsave;m.supplementary_figures()
            else:m.main()
    finally:Figure.savefig=original
if __name__=='__main__':
    workflow();figure2();figure4();supplementary_black_text()
    for n,old in [(3,'Figure_3_human_MS_microenvironment_stage_revised'),(5,'Figure_5_human_validation')]:
        for ext in ['png','pdf']:shutil.copy2(R/'output/figures'/f'{old}.{ext}',F/f'Figure_{n}.{ext}')
