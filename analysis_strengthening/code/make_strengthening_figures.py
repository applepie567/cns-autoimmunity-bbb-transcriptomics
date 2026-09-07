from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
R=Path(__file__).resolve().parents[1];O=R/'output';F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def save(fig,name):
    fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)

def study_design():
    fig,ax=plt.subplots(figsize=(8,6.3));ax.set(xlim=(0,1),ylim=(.1,1));ax.axis('off')
    def box(x,y,w,h,title,body):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',facecolor='#f3f5f7',edgecolor='#516473',linewidth=1))
        ax.text(x+.018,y+h-.025,title,va='top',fontsize=10,fontweight='bold')
        ax.text(x+.018,y+h-.078,body,va='top',fontsize=9,linespacing=1.4)
    box(.02,.71,.45,.265,'Acute EAE reference cohorts','GSE199460: brain, 2 EAE and 3 controls\nGSE210776: spinal cord, 4 and 3\nGSE95401: spinal cord, 3 and 3\nOne expression profile per animal')
    box(.53,.71,.45,.265,'Human reference','GSE279183: chronic active MS\n4 MS donors and 6 controls\nEndothelial nuclei combined by donor\nUK MS Tissue Bank brain samples')
    box(.13,.425,.74,.22,'Matched gene benchmark','7,705 strict orthologues with the same effect estimator\nMouse cohort agreement compared with human correspondence\nAnimal and donor resampling and subject omission')
    for x in [.245,.755]:ax.add_patch(FancyArrowPatch((x,.695),(.5,.66),arrowstyle='-|>',mutation_scale=14,color='#516473'))
    box(.02,.135,.45,.235,'Additional human validation','Macnair et al. 2025, cohort I\nChronic active lesions and healthy WM\nSame donors and capillary restriction\nEndothelial composition assessed')
    box(.53,.135,.45,.235,'Spatial and functional context','GSE279183: vascular niches\nGSE208747: lesion categories\nGSE284005: paired MERFISH regions\nPublished barrier assays as support')
    save(fig,'Figure_1_revised_workflow')

def benchmark():
    d=pd.read_csv(O/'S19_matched_gene_effects.csv').set_index('gene')
    names=['GSE199460','GSE210776','GSE95401','GSE279183']
    corr=d[names].corr(method='spearman')
    fig=plt.figure(figsize=(8.6,8.2));gs=fig.add_gridspec(2,2,height_ratios=[1.1,1],hspace=.45,wspace=.55)
    a=fig.add_subplot(gs[0,0]);a.imshow(corr,vmin=0,vmax=1,cmap='Blues')
    for i in range(4):
        for j in range(4):a.text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center',color='black',bbox=dict(facecolor='white',alpha=.72,edgecolor='none',pad=1))
    labels=['Mouse brain','Mouse spinal 1','Mouse spinal 2','Human MS']
    a.set_xticks(range(4),labels,rotation=40,ha='right',fontsize=9);a.set_yticks(range(4),labels,fontsize=9)
    a.set_title('A  Same 7,705 genes',loc='left',fontweight='bold',pad=12)
    b=fig.add_subplot(gs[0,1]);s=pd.read_csv(O/'S19_concordance_benchmark.csv').set_index('comparison')
    keys=['mean_within_mouse','mean_individual_mouse_human','pooled_mouse__GSE279183','within_minus_individual_cross']
    labs=['Mean mouse to mouse','Mean individual mouse to human','Pooled mouse to human','Difference between means']
    for i,k in enumerate(keys):
        r=s.loc[k];b.plot([r.ci_low,r.ci_high],[i,i],color='#416b83',lw=2);b.plot(r.rho_or_difference,i,'o',color='#416b83');b.text(.64,i,f'{r.rho_or_difference:.2f}',va='center',ha='right',fontsize=10)
        b.text(.02,i-.27,labs[i],fontsize=8.5,ha='left',va='center')
    b.set_yticks([]);b.set_ylim(3.45,-.6);b.set_xlim(-.03,.66);b.axvline(0,color='.55',ls=':',lw=.8);b.set_xlabel('Correlation or difference')
    b.set_title('B  Subject bootstrap ranges',loc='left',fontweight='bold',pad=12,fontsize=10.5)
    c=fig.add_subplot(gs[1,:]);pairs=[('GSE199460__GSE210776','Brain with spinal 1'),('GSE199460__GSE95401','Brain with spinal 2'),('GSE210776__GSE95401','Spinal 1 with spinal 2'),('GSE199460__GSE279183','Brain with human'),('GSE210776__GSE279183','Spinal 1 with human'),('GSE95401__GSE279183','Spinal 2 with human')]
    for i,(k,l) in enumerate(pairs):
        r=s.loc[k];col='#376680' if i<3 else '#965a3d';c.plot([r.ci_low,r.ci_high],[i,i],color=col,lw=2);c.plot(r.rho_or_difference,i,'o',color=col)
        c.text(.68,i,f'{r.rho_or_difference:.2f}  ({r.ci_low:.2f} to {r.ci_high:.2f})',va='center',fontsize=9)
    c.set_yticks(range(6),[x[1] for x in pairs]);c.invert_yaxis();c.set_xlim(-.15,.98);c.axvline(0,color='.55',ls=':',lw=.8);c.set_xlabel("Spearman's correlation");c.set_title('C  Individual cohort comparisons',loc='left',fontweight='bold',pad=12)
    fig.subplots_adjust(left=.18,right=.97,bottom=.12,top=.94)
    save(fig,'Figure_4_concordance_benchmark')

def validation():
    prop=pd.read_csv(O/'S20_endothelial_subtype_proportions.csv')
    prop=prop.sort_values(['lesion_type','individual_id_anon'],ascending=[False,True]).reset_index(drop=True)
    fig=plt.figure(figsize=(8.2,8.8));gs=fig.add_gridspec(3,1,height_ratios=[1,1.1,.85],hspace=.9)
    a=fig.add_subplot(gs[0]);bottom=np.zeros(len(prop))
    cats=[('Endocyte_Capillary','Capillary','#63859a'),('Endocyte_Venous','Venous','#c8b089'),('Endocyte_Stressed','Stressed','#b36d51'),('Endocyte_Arterial','Arterial','#b9c5c7')]
    for key,label,col in cats:
        a.bar(np.arange(len(prop)),prop[key],bottom=bottom,label=label,color=col,width=.8,edgecolor='white',linewidth=.3);bottom+=prop[key]
    a.set_xticks(range(len(prop)),[x.replace('Ind','') for x in prop.individual_id_anon],rotation=45,ha='right',fontsize=8)
    a.set_ylabel('Fraction of recovered\nendothelial nuclei');a.set_ylim(0,1.06)
    a.text(1.5,1.10,'Healthy WM',ha='center',fontsize=9);a.text(7,1.10,'Chronic active MS',ha='center',fontsize=9)
    a.set_title('A  Endothelial subtype composition',loc='left',fontweight='bold',y=1.25,fontsize=11)
    a.legend(frameon=False,bbox_to_anchor=(1.01,1.02),loc='upper left',fontsize=9)
    b=fig.add_subplot(gs[1]);rng=np.random.default_rng(53)
    for group,(label,title) in enumerate([('All_EC','All endothelial'),('All_EC_capillary_donors','Same six donors'),('Capillary','Capillary')]):
        d=pd.read_csv(O/f'S20_{label}_sample_scores.csv')
        for j,(condition,color) in enumerate([('WM','#50758a'),('CAL','#a06748')]):
            v=d.loc[d.lesion.eq(condition),'BBB_specialization'].to_numpy();x=group*3+j
            b.scatter(x+rng.uniform(-.12,.12,len(v)),v,color=color,s=26,alpha=.95,label='Healthy WM' if group==0 and j==0 else ('Chronic active MS' if group==0 else None))
            b.plot([x-.25,x+.25],[v.mean(),v.mean()],color='black',lw=1.4)
        b.text(group*3+.5,1.05,title,ha='center',transform=b.get_xaxis_transform(),fontsize=9)
    b.set_xticks([0,1,3,4,6,7],['CTRL','MS']*3,fontsize=9);b.set_xlim(-.6,7.6);b.set_ylabel('BBB specialization score');b.axhline(0,ls=':',lw=.8,color='.55')
    b.set_title('B  Directional BBB program change',loc='left',fontweight='bold',y=1.20,fontsize=11)
    b.legend(frameon=False,bbox_to_anchor=(1.01,1.02),loc='upper left',fontsize=8)
    c=fig.add_subplot(gs[2]);cor=pd.read_csv(O/'S20_common_gene_subtype_sensitivity.csv').set_index(['subset','reference'])
    labels=['All_EC','All_EC_capillary_donors','Capillary']
    for i,lab in enumerate(labels):
        for j,(ref,color,mark) in enumerate([('mouse_pooled_g','#416b83','o'),('human_hedges_g','#965a3d','s')]):
            r=cor.loc[(lab,ref),'rho'];y=i+(j-.5)*.22;c.plot(r,y,marker=mark,color=color,ms=6,ls='none',label=('Pooled mouse' if j==0 else 'Original human') if i==0 else None);c.text(r+.009,y,f'{r:.2f}',va='center',fontsize=8)
    c.set_yticks(range(3),['All EC, full cohort','All EC, same six donors','Capillary, same six donors'],fontsize=9);c.set_ylim(2.6,-.6);c.set_xlim(-.11,.22);c.axvline(0,ls=':',lw=.8,color='.55');c.set_xlabel("Spearman's correlation")
    c.set_title('C  Gene correspondence on 5,811 shared genes',loc='left',fontweight='bold',pad=16,fontsize=11)
    c.legend(frameon=False,bbox_to_anchor=(1.01,1.04),loc='upper left',fontsize=9)
    fig.subplots_adjust(left=.21,right=.78,bottom=.07,top=.90)
    save(fig,'Figure_5_human_validation')

if __name__=='__main__':study_design();benchmark();validation()
