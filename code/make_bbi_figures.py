#!/usr/bin/env python3
"""Generate the five BBI main figures and supporting sensitivity figures."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SOURCE = ROOT / "supplementary_tables"
OUT = ROOT / "figures" / "main"
OUT.mkdir(parents=True, exist_ok=True)
SUPP_OUT = ROOT / "figures" / "supplementary"
SUPP_OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "immune": "#C84C35",
    "barrier": "#2F6F9F",
    "ecm": "#76528B",
    "control": "#7D8790",
    "acute": "#C84C35",
    "bevac": "#2A9D8F",
    "positive": "#B4473D",
    "negative": "#3A78A1",
    "light": "#EFF3F5",
    "dark": "#25313B",
    "gold": "#D29B32",
}

MODULE_LABEL = {
    "Endothelial_immune_activation": "Endothelial immune\nactivation",
    "BBB_specialization": "BBB specialization",
    "Structural_ECM_remodeling": "Structural and ECM\nremodeling",
}
MODULE_COLOR = {
    "Endothelial_immune_activation": COLORS["immune"],
    "BBB_specialization": COLORS["barrier"],
    "Structural_ECM_remodeling": COLORS["ecm"],
}


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


MAIN_PANEL_TITLE_X = -0.065


def panel(ax: plt.Axes, letter: str, title: str, title_x: float = 0.0) -> None:
    """Use one title string so the panel letter is exactly two spaces from the title."""
    ax.set_title(f"{letter}  {title}", loc="left", x=title_x, ha="left",
                 fontweight="bold", pad=8)


def header_panel(ax: plt.Axes, letter: str, title: str, title_x: float = 0.0) -> None:
    """Place the panel letter and title on one baseline with a two-space gap."""
    ax.text(title_x, 0.52, f"{letter}  {title}", transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="center", ha="left")


def box(ax: plt.Axes, xy: tuple[float, float], wh: tuple[float, float], text: str, fc: str, ec: str = "#B7C1C8", fs: float = 7.0) -> None:
    p = FancyBboxPatch(xy, wh[0], wh[1], boxstyle="round,pad=0.018,rounding_size=0.018", fc=fc, ec=ec, lw=0.8)
    ax.add_patch(p)
    ax.text(xy[0] + wh[0] / 2, xy[1] + wh[1] / 2, text, ha="center", va="center", fontsize=fs, color="#000000")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, lw=0.9, color="#72808A"))


def save(fig: plt.Figure, stem: str, out_dir: Path = OUT) -> None:
    png = out_dir / f"{stem}.png"
    tif = out_dir / f"{stem}.tiff"
    png_tmp = out_dir / f".{stem}.tmp.png"
    tif_tmp = out_dir / f".{stem}.tmp.tiff"
    fig.savefig(png_tmp, dpi=350, bbox_inches="tight", format="png")
    fig.savefig(tif_tmp, dpi=600, bbox_inches="tight", format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    os.replace(png_tmp, png)
    os.replace(tif_tmp, tif)
    plt.close(fig)


def fig1_overview() -> None:
    fig = plt.figure(figsize=(7.4, 6.8))
    fig.subplots_adjust(left=.055, right=.985, top=.95, bottom=.06, hspace=.38, wspace=.28)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.08], width_ratios=[1.12, 0.88])

    ax = fig.add_subplot(gs[0, :]); ax.axis("off")
    fig.text(.055, .945, "A  Study logic and analytical contributions", fontsize=9.0,
             fontweight="bold", va="center", ha="left", color="#000000")
    box(ax, (0.02, 0.56), (0.19, 0.23), "Acute EAE\n3 endothelial cohorts", "#FBEDEA", ec="#DBA99F", fs=7.5)
    box(ax, (0.28, 0.56), (0.19, 0.23), "Biological sample scores\nrandom effects\nmeta analysis", "#F4F6F7", fs=7.2)
    box(ax, (0.54, 0.56), (0.19, 0.23), "Human MS\nsingle nucleus RNA\nand spatial data", "#EAF2F8", ec="#A8C5D8", fs=7.2)
    box(ax, (0.80, 0.56), (0.18, 0.23), "Stage and local\nmicroenvironment", "#EEF5F3", ec="#A8CEC5", fs=7.5)
    for x in [(0.225, 0.675, 0.265, 0.675), (0.485, 0.675, 0.525, 0.675),
              (0.745, 0.675, 0.785, 0.675)]:
        arrow(ax, (x[0], x[1]), (x[2], x[3]))
    box(ax, (0.20, 0.13), (0.24, 0.21), "Comparison of disease effects\nbetween mouse and human", "#FFF7E5", ec="#DFC27B", fs=7.2)
    box(ax, (0.56, 0.13), (0.24, 0.21), "VEGF-A blockade\ntranscriptional audit", "#E9F5F2", ec="#9FC9BF", fs=7.5)

    ax = fig.add_subplot(gs[1, 0]); ax.axis("off")
    fig.text(.055, .493, "B  Primary datasets and independent units", fontsize=9.0,
             fontweight="bold", va="center", ha="left", color="#000000")
    rows = [
        ("GSE210776", "Mouse scRNA, EAE and BEVAC groups", "mouse", "meta analysis and audit"),
        ("GSE199460", "Mouse CD31 scRNA, control and acute", "mouse", "independent validation"),
        ("GSE95401", "Mouse EC bulk, EAE time course", "sample", "temporal support"),
        ("GSE279183", "Human snRNA and spatial MS data", "donor", "local MS context"),
        ("GSE208747", "Human Visium, NAWM and lesions", "donor", "lesion stage"),
        ("GSE284005", "Human MERFISH, MS regions", "donor", "vascular audit"),
    ]
    col_x = [0.00, 0.19, 0.68, 0.79]
    headers = ["Dataset", "Data and groups", "Unit", "Role"]
    for x, h in zip(col_x, headers): ax.text(x, 0.92, h, fontweight="bold", fontsize=7.2, color="#000000", transform=ax.transAxes)
    for i, row in enumerate(rows):
        y = 0.81 - i * 0.123
        if i % 2 == 0: ax.add_patch(plt.Rectangle((0, y - .047), 1, .10, transform=ax.transAxes, fc="#F6F8F9", ec="none"))
        for x, val in zip(col_x, row): ax.text(x, y, val, fontsize=5.9, color="#000000", va="center", transform=ax.transAxes)

    ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
    fig.text(.570, .493, "C  Prespecified endothelial states", fontsize=9.0,
             fontweight="bold", va="center", ha="left", color="#000000")
    states = [
        ("Immune activation", "STAT1, STAT3, CXCL9, CXCL10, B2M,\nTAP1, ICAM1, VCAM1", "#FBEDEA", COLORS["immune"]),
        ("BBB specialization", "CLDN5, OCLN, MFSD2A,\nSLC2A1, ABCG2, Wnt-BBB genes", "#EAF2F8", COLORS["barrier"]),
        ("Structural and ECM", "CAV1, CAV2, PLVAP, MMP2, MMP9,\nCOL4A1, COL4A2, FN1, SPARC", "#F1ECF4", COLORS["ecm"]),
    ]
    for i, (name, genes, fc, ec) in enumerate(states):
        y = 0.69 - i * 0.29
        box(ax, (0.04, y), (0.92, 0.20), f"{name}\n{genes}", fc, ec=ec, fs=7.3)
    save(fig, "Figure_1_study_design")


def forest_panel(ax: plt.Axes, block: pd.DataFrame, pooled: pd.Series, module: str, letter: str) -> None:
    short = {"Endothelial_immune_activation":"Immune activation","BBB_specialization":"BBB specialization","Structural_ECM_remodeling":"Structural and ECM"}[module]
    panel(ax, letter, short, title_x=MAIN_PANEL_TITLE_X)
    order = ["GSE210776", "GSE199460_all", "GSE95401"]
    labels = ["GSE210776", "GSE199460", "GSE95401", "Pooled"]
    y = np.arange(4)[::-1]
    for i, ds in enumerate(order):
        r = block[block.dataset == ds].iloc[0]
        ax.plot([r.ci_low, r.ci_high], [y[i], y[i]], color="#65737D", lw=1)
        ax.scatter(r.hedges_g, y[i], s=28, color=MODULE_COLOR[module], edgecolor="white", linewidth=.5, zorder=3)
    yp = y[3]
    ax.plot([pooled.mKH_ci_low, pooled.mKH_ci_high], [yp, yp], color=MODULE_COLOR[module], lw=1.1, ls=(0, (2, 2)))
    ax.plot([pooled.ci_low, pooled.ci_high], [yp, yp], color=MODULE_COLOR[module], lw=2.4)
    ax.scatter(pooled.pooled_g, yp, s=48, marker="D", color=MODULE_COLOR[module], edgecolor="white", linewidth=.5, zorder=3)
    ax.axvline(0, color="#AAB2B7", lw=.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Hedges’ g (EAE − control)")
    ax.set_ylim(-.7, 3.7)
    ax.grid(axis="x", color="#E4E8EA", lw=.6)
    ax.text(.02, .02, f"Wald p={pooled.p:.3g}\nmKH p={pooled.mKH_p:.3g}\nI²={pooled.I2:.1f}%", transform=ax.transAxes, fontsize=5.7, va="bottom")


def fig2_meta() -> None:
    effects = pd.read_csv(RESULTS / "dataset_module_effects.csv")
    effects = effects[effects.dataset.isin(["GSE210776", "GSE199460_all", "GSE95401"]) & effects.contrast.str.startswith("Acute")]
    meta = pd.read_csv(RESULTS / "acute_EAE_random_effects_meta.csv")
    loo = pd.read_csv(RESULTS / "acute_EAE_leave_one_dataset_out.csv")
    sample_scores = pd.read_csv(RESULTS / "sample_level_module_scores.csv")
    fig = plt.figure(figsize=(7.7, 6.7))
    fig.subplots_adjust(left=.09,right=.98,top=.94,bottom=.09,wspace=.55,hspace=.48)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.02, .98])
    modules = list(MODULE_LABEL)
    for i, mod in enumerate(modules):
        ax = fig.add_subplot(gs[0, i])
        forest_panel(ax, effects[effects.module == mod], meta[meta.module == mod].iloc[0], mod, chr(65+i))
        if mod == "Endothelial_immune_activation": ax.set_xlim(-2, 11)
        elif mod == "BBB_specialization": ax.set_xlim(-6, 2)
        else: ax.set_xlim(-2, 4)

    bottom = gs[1, :].subgridspec(1, 2, width_ratios=[2.25, 1], wspace=.58)
    ax = fig.add_subplot(bottom[0, 0]); panel(ax, "D", "Estimates after omitting each dataset", title_x=MAIN_PANEL_TITLE_X)
    mods = ["Endothelial_immune_activation", "BBB_specialization", "Structural_ECM_remodeling"]
    omissions = ["GSE210776", "GSE199460_all", "GSE95401"]
    mat = np.array([[loo[(loo.module==m)&(loo.omitted_dataset==o)].pooled_g.iloc[0] for o in omissions] for m in mods])
    im=ax.imshow(mat, cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]): ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=7, color="white" if abs(mat[i,j])>2.8 else COLORS["dark"])
    ax.set_xticks(range(3), ["omit\nGSE210776", "omit\nGSE199460", "omit\nGSE95401"])
    ax.set_yticks(range(3), ["Immune activation", "BBB specialization", "Structural and ECM"])
    cb=fig.colorbar(im, ax=ax, fraction=.028, pad=.018)
    cb.set_label("Pooled Hedges’ g", fontsize=6, labelpad=2)
    cb.ax.tick_params(labelsize=5.5)
    ax = fig.add_subplot(bottom[0, 1]); panel(ax, "E", "GSE95401 sampled stages", title_x=MAIN_PANEL_TITLE_X)
    stage = sample_scores[sample_scores.dataset.eq("GSE95401")]
    stage_order = ["Control", "Acute", "Subacute", "Chronic"]
    final_labels = {
        "Endothelial_immune_activation": ("Immune", .16),
        "BBB_specialization": ("BBB", -.39),
        "Structural_ECM_remodeling": ("Structural and ECM", .03),
    }
    for mod in MODULE_LABEL:
        means = [stage.loc[stage.group.eq(group), mod].mean() for group in stage_order]
        sems = [stage.loc[stage.group.eq(group), mod].sem() for group in stage_order]
        ax.errorbar(range(4), means, yerr=sems, marker="o", ms=4, lw=1.15, capsize=1.7, color=MODULE_COLOR[mod])
        label, label_y = final_labels[mod]
        ax.text(3.12, label_y, label, color=MODULE_COLOR[mod], fontsize=5.4, va="center")
    ax.axhline(0, color="#D6DCDF", lw=.7)
    ax.set_xlim(-.12, 4.05); ax.set_ylim(-1.02, .88)
    ax.set_xticks(range(4), ["Control", "Acute", "Subacute", "Chronic"], rotation=32, ha="right")
    ax.set_ylabel("Mean score ± SEM", labelpad=2)
    save(fig, "Figure_2_acute_EAE_meta")


def fig3_human_ms() -> None:
    scores = pd.read_csv(RESULTS / "sample_level_module_scores.csv")
    h = scores[scores.dataset == "GSE279183"].copy()
    spatial = pd.read_csv(RESULTS / "GSE279183_inflammatory_vascular_microenvironment_source_effects.csv")
    barrier = pd.read_csv(RESULTS / "GSE279183_inflammatory_vascular_microenvironment_barrier_effects.csv")
    stage = pd.read_csv(RESULTS / "GSE208747_stage_locked_modules_q30_ECadj.csv")
    merfish_scores = pd.read_csv(RESULTS / "GSE284005_donor_region_module_scores.csv")
    merfish_cells = pd.read_csv(RESULTS / "GSE284005_endothelial_cell_counts.csv")
    merfish_summary = pd.read_csv(RESULTS / "GSE284005_paired_summary.csv")
    fig = plt.figure(figsize=(7.7, 9.25))
    fig.subplots_adjust(left=.10,right=.98,top=.96,bottom=.065,hspace=.52,wspace=.60)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.08, .82], width_ratios=[1.04, .96])

    sub = gs[0,0].subgridspec(2,3, height_ratios=[.12,.88], hspace=.18, wspace=.54)
    holder = fig.add_subplot(sub[0,:]); holder.axis("off"); panel(holder, "A", "Endothelial profiles summarized for each donor", title_x=MAIN_PANEL_TITLE_X)
    rng=np.random.default_rng(20260901)
    groups=["Control","Chronic_active","Chronic_inactive"]; glabel=["Control","Chronic\nactive","Chronic\ninactive"]
    for j,mod in enumerate(MODULE_LABEL):
        ax=fig.add_subplot(sub[1,j])
        for i,g in enumerate(groups):
            vals=h.loc[h.group==g,mod].dropna().to_numpy()
            x=i+rng.uniform(-.08,.08,len(vals)); ax.scatter(x,vals,s=24,color=[COLORS['control'],COLORS['immune'],COLORS['ecm']][i],edgecolor='white',lw=.5,zorder=3)
            ax.plot([i-.12,i+.12],[np.mean(vals),np.mean(vals)],color=COLORS['dark'],lw=1.2)
        ax.axhline(0,color='#D6DCDF',lw=.7); ax.set_xticks(range(3),glabel,rotation=30,ha='right'); ax.set_title(MODULE_LABEL[mod],fontsize=7.2)
        if j==0: ax.set_ylabel("Module score")
        ax.grid(axis='y',color='#EDF0F2',lw=.5)

    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Pathways in MS inflammatory microenvironments",title_x=MAIN_PANEL_TITLE_X)
    s=spatial.sort_values('delta_VI_vs_PPWM'); colors=[COLORS['gold'] if p=='WNT' else COLORS['immune'] for p in s.pathway]
    ax.barh(s.pathway,s.delta_VI_vs_PPWM,color=colors,height=.68)
    for i,r in enumerate(s.itertuples()): ax.text(r.delta_VI_vs_PPWM+.035,i,f"FDR {r.FDR_VI_vs_PPWM:.3f}",va='center',fontsize=6.2)
    ax.axvline(0,color='#8F999F',lw=.7); ax.set_xlabel("Paired difference in pathway scores (n=7 donors)"); ax.set_xlim(0,max(s.delta_VI_vs_PPWM)*1.45)

    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","BBB modules in MS inflammatory microenvironments",title_x=MAIN_PANEL_TITLE_X)
    b=barrier.copy(); y=np.arange(len(b)); cols=[COLORS['barrier']]*len(b)
    ax.barh(y,b.delta_microenvironment_vs_control_WM,color=cols,height=.62)
    labels=["Tight junction","BBB transport","Wnt-BBB identity"]
    ax.set_yticks(y,labels); ax.axvline(0,color='#8F999F',lw=.7); ax.set_xlabel("Microenvironment − control white matter")
    ax.invert_yaxis()
    ax.set_xlim(-1.55,.34)
    for i,r in enumerate(b.itertuples()):
        ax.text(.03,i,f"FDR {r.FDR_BH:.3f}",ha='left',va='center',fontsize=6.3,color=COLORS['dark'])
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","GSE208747 lesion stage sensitivity",title_x=MAIN_PANEL_TITLE_X)
    chosen=["IFN_antigen_presentation","PLVAP_permeability","Adhesion_trafficking","TGF_response","Caveolae_structural","ECM_protease","Tight_junction","BBB_transport_identity","Wnt_BBB"]
    contrasts=["Active_vs_NAWM","Mixed_vs_NAWM"]
    mat=np.array([[stage[(stage.module==m)&(stage.contrast==c)].delta_median.iloc[0] for c in contrasts] for m in chosen])
    im=ax.imshow(mat,cmap='RdBu_r',vmin=-1.6,vmax=1.6,aspect='auto')
    for i in range(mat.shape[0]):
        for j in range(2): ax.text(j,i,f"{mat[i,j]:.2f}",ha='center',va='center',fontsize=6,color='white' if abs(mat[i,j])>.85 else COLORS['dark'])
    ax.set_xticks([0,1],["Active − NAWM\n3 vs 3 donors","Mixed − NAWM\n6 vs 3 donors"])
    ax.set_yticks(range(len(chosen)),[m.replace('_',' ') for m in chosen])
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.03); cb.set_label("EC-adjusted median difference")
    ax=fig.add_subplot(gs[2,0]); panel(ax,"E","GSE284005 paired endothelial programs",title_x=MAIN_PANEL_TITLE_X)
    merfish_modules = ["Endothelial immune activation", "Adhesion and trafficking", "IFN and antigen presentation"]
    module_labels = ["Immune", "Adhesion", "Interferon and\nantigen presentation"]
    for j,module in enumerate(merfish_modules):
        wide=merfish_scores[merfish_scores.module.eq(module)].pivot(index='donor',columns='region',values='score')
        base=j*3
        for _,row in wide.iterrows():
            ax.plot([base,base+1],[row['DMWM'],row['Vas_Imm']],color='#BDBDBD',lw=.8,zorder=1)
            ax.scatter(base,row['DMWM'],color=COLORS['control'],s=20,zorder=2)
            ax.scatter(base+1,row['Vas_Imm'],color=COLORS['positive'],s=20,zorder=2)
    ax.axhline(0,color='#D6DCDF',lw=.7)
    ax.set_xticks([j*3+.5 for j in range(3)],module_labels,rotation=20,ha='right')
    ax.set_ylabel('Module score')
    ax.scatter([],[],color=COLORS['control'],label='DMWM'); ax.scatter([],[],color=COLORS['positive'],label='Vascular inflammatory')
    ax.legend(frameon=False,fontsize=5.8,ncol=2,loc='lower left')

    ax=fig.add_subplot(gs[2,1]); panel(ax,"F","GSE284005 endothelial composition",title_x=MAIN_PANEL_TITLE_X)
    frac=merfish_cells.pivot(index='donor',columns='region',values='stress_or_inflammatory_fraction')
    for _,row in frac.iterrows():
        ax.plot([0,1],[row['DMWM'],row['Vas_Imm']],color='#BDBDBD',lw=.9)
        ax.scatter(0,row['DMWM'],color=COLORS['control'],s=24)
        ax.scatter(1,row['Vas_Imm'],color=COLORS['positive'],s=24)
    ax.set_xticks([0,1],['DMWM','Vascular\ninflammatory'])
    ax.set_ylabel('Fraction of stress responsive\nand inflammatory ECs'); ax.set_ylim(0,1)
    save(fig,"Figure_3_human_MS_microenvironment_stage")


def fig4_cross_species(
    stem: str = 'Figure_4_cross_species_concordance',
    background_color: str = '#BFC4C8',
    background_alpha: float = .48,
    background_size: float = 2.1,
) -> None:
    genes=pd.read_csv(RESULTS / 'cross_species_gene_effects.csv')
    summary=pd.read_csv(RESULTS / 'cross_species_concordance_summary.csv')
    orthology=pd.read_csv(RESULTS / 'strict_orthology_sensitivity_summary.csv')
    fig=plt.figure(figsize=(7.55,6.7))
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.09,hspace=.45,wspace=.42)
    gs=fig.add_gridspec(2,2,height_ratios=[1.08,.92],width_ratios=[1.08,.92])

    ax=fig.add_subplot(gs[0,0]); panel(ax,'A','Gene effects in animals and donors across species',title_x=MAIN_PANEL_TITLE_X)
    ax.scatter(genes.mouse_pooled_g,genes.human_hedges_g,s=background_size,color=background_color,alpha=background_alpha,
               linewidth=0,rasterized=True)
    state_genes=genes[genes.prespecified_state_gene].copy()
    for state,color in MODULE_COLOR.items():
        b=state_genes[state_genes.state.str.contains(state,regex=False)]
        ax.scatter(b.mouse_pooled_g,b.human_hedges_g,s=22,color=color,alpha=.82,edgecolor='white',lw=.35,label=MODULE_LABEL[state].replace('\n',' '))
    for gene in ['CXCL9','CXCL10','CCL2','TAP1','MFSD2A','OCLN','TIMP1','VIM','VCAM1','LEF1']:
        b=state_genes[state_genes.gene==gene]
        if b.empty: continue
        r=b.iloc[0]; ax.annotate(gene,(r.mouse_pooled_g,r.human_hedges_g),xytext=(3,3),textcoords='offset points',fontsize=5.4)
    ax.axhline(0,color='#9DA7AD',lw=.7); ax.axvline(0,color='#9DA7AD',lw=.7)
    ax.set_xlabel('Mouse acute EAE pooled Hedges’ g'); ax.set_ylabel('Human chronic active MS Hedges’ g')
    ax.legend(frameon=False,fontsize=5.4,loc='lower right')

    order=['All_common_genes_sensitivity','All_prespecified_state_genes','Endothelial_immune_activation','BBB_specialization','Structural_ECM_remodeling']
    labels=['All common genes','All state genes','Immune activation','BBB specialization','Structural and ECM']
    s=summary.set_index('gene_set').loc[order].reset_index()
    ax=fig.add_subplot(gs[0,1]); panel(ax,'B','Rank correlation of gene effects',title_x=MAIN_PANEL_TITLE_X)
    y=np.arange(len(s))[::-1]
    colors=['#7D8790','#25313B',COLORS['immune'],COLORS['barrier'],COLORS['ecm']]
    for i,r in enumerate(s.itertuples()):
        ax.plot([r.rho_ci_low,r.rho_ci_high],[y[i],y[i]],color=colors[i],lw=1.2)
        ax.scatter(r.spearman_rho,y[i],s=38,color=colors[i],edgecolor='white',lw=.5,zorder=3)
        ax.text(.88,y[i],f"n={int(r.n_genes)}",va='center',fontsize=5.8)
    ax.axvline(0,color='#9DA7AD',lw=.7); ax.set_xlim(-.72,1.02); ax.set_yticks(y,labels); ax.set_xlabel('Spearman ρ (95% CI)'); ax.grid(axis='x',color='#E6EAEC',lw=.5)

    mapping_order = [
        ('Uppercase-symbol match','All shared genes'),
        ('HCOP reciprocal 1:1; Ensembl+NCBI support','All shared genes'),
        ('Uppercase-symbol match','Prespecified focused genes'),
        ('HCOP reciprocal 1:1; Ensembl+NCBI support','Prespecified focused genes'),
    ]
    mapping_labels=['Symbol\nall','Strict 1:1\nall','Symbol\nfocused','Strict 1:1\nfocused']
    o=pd.concat([orthology[(orthology.mapping.eq(m))&(orthology.scope.eq(scope))].iloc[[0]] for m,scope in mapping_order],ignore_index=True)
    mapping_colors=['#688492','#2F769B','#B7833C','#94507F']

    ax=fig.add_subplot(gs[1,0]); panel(ax,'C','Sensitivity to orthology mapping',title_x=MAIN_PANEL_TITLE_X)
    x=np.arange(len(o)); yerr=np.vstack([o.spearman_rho-o.ci_low,o.ci_high-o.spearman_rho])
    ax.bar(x,o.spearman_rho,color=mapping_colors,width=.66)
    ax.errorbar(x,o.spearman_rho,yerr=yerr,fmt='none',ecolor='#25313B',elinewidth=1,capsize=3)
    ax.axhline(0,color='#8F999F',lw=.7); ax.set_ylim(-.15,.76); ax.set_ylabel('Spearman ρ')
    ax.set_xticks(x,mapping_labels)
    for i,r in enumerate(o.itertuples()): ax.text(i,r.ci_high+.035,f"n={int(r.n_genes):,}",ha='center',fontsize=5.7)

    ax=fig.add_subplot(gs[1,1]); panel(ax,'D','Directional agreement after mapping',title_x=MAIN_PANEL_TITLE_X)
    pct=o.same_direction_fraction*100
    bars=ax.bar(x,pct,color=mapping_colors,width=.66)
    ax.axhline(50,color='#4F5559',lw=.8,ls='--'); ax.set_ylim(0,80); ax.set_ylabel('Genes with the same direction (%)')
    ax.set_xticks(x,mapping_labels)
    for bar,r in zip(bars,o.itertuples()): ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+2,f"n={int(r.n_genes):,}",ha='center',fontsize=5.7)
    save(fig,stem)


def fig5_vegf() -> None:
    scores=pd.read_csv(RESULTS/'sample_level_module_scores.csv'); scores=scores[(scores.dataset=='GSE210776')&scores.group.isin(['CFA','Acute','Bevacizumab'])]
    vegf=pd.read_csv(RESULTS/'VEGFA_three_arm_transcriptional_perturbation.csv')
    ven=pd.read_csv(RESULTS/'VEGFA_venous_source_effects.csv'); ven=ven[ven.subtype=='Venous']
    functional=pd.read_csv(RESULTS/'published_source_barrier_and_vascular_summary.csv')
    ihc=pd.read_csv(RESULTS/'published_human_MS_IHC_descriptive_summary.csv')
    fig=plt.figure(figsize=(7.7,7.0)); fig.subplots_adjust(left=.09,right=.985,top=.95,bottom=.085,hspace=.47,wspace=.70); gs=fig.add_gridspec(2,2,height_ratios=[1,.98],width_ratios=[1.07,.93])

    sub=gs[0,0].subgridspec(3,3,height_ratios=[.13,.10,.77],hspace=0,wspace=.58); holder=fig.add_subplot(sub[0,:]); holder.axis('off'); header_panel(holder,'A','All endothelial cells across three groups',title_x=MAIN_PANEL_TITLE_X)
    rng=np.random.default_rng(20260901); groups=['CFA','Acute','Bevacizumab']; gl=['CFA','Acute\nEAE','Bev.']
    for j,mod in enumerate(MODULE_LABEL):
        ax=fig.add_subplot(sub[2,j])
        for i,g in enumerate(groups):
            vals=scores.loc[scores.group==g,mod].to_numpy(); x=i+rng.uniform(-.07,.07,len(vals)); col=[COLORS['control'],COLORS['acute'],COLORS['bevac']][i]
            ax.scatter(x,vals,s=28,color=col,edgecolor='white',lw=.5); ax.plot([i-.12,i+.12],[vals.mean(),vals.mean()],color=COLORS['dark'],lw=1.2)
        ax.axhline(0,color='#D6DCDF',lw=.7); ax.set_xticks(range(3),gl,rotation=25,ha='right'); ax.set_title(MODULE_LABEL[mod],fontsize=7.1); ax.grid(axis='y',color='#EDF0F2',lw=.5)
        if j==0: ax.set_ylabel('Module score')

    sub_b=gs[0,1].subgridspec(2,1,height_ratios=[.13,.87],hspace=.12)
    holder=fig.add_subplot(sub_b[0,0]); holder.axis('off'); header_panel(holder,'B','Disease, treatment and residual vectors',title_x=MAIN_PANEL_TITLE_X)
    ax=fig.add_subplot(sub_b[1,0])
    v=vegf[~vegf.module.isin(MODULE_LABEL)].copy(); mat=v[['acute_minus_CFA','bevacizumab_minus_acute','bevacizumab_minus_CFA']].to_numpy()
    im=ax.imshow(mat,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_xticks(range(3),['Acute − CFA','Bevacizumab\n− acute','Bevacizumab\n− CFA'],rotation=12,ha='right'); ax.set_yticks(range(len(v)),[m.replace('_',' ') for m in v.module])
    for i in range(mat.shape[0]):
        for j in range(3): ax.text(j,i,f"{mat[i,j]:.2f}",ha='center',va='center',fontsize=5.3,color='white' if abs(mat[i,j])>1.15 else COLORS['dark'])
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.03); cb.set_label('Difference in module score')

    sub_c=gs[1,0].subgridspec(2,1,height_ratios=[.13,.87],hspace=.12)
    holder=fig.add_subplot(sub_c[0,0]); holder.axis('off'); header_panel(holder,'C','Venous endothelial source analysis',title_x=MAIN_PANEL_TITLE_X)
    ax=fig.add_subplot(sub_c[1,0])
    features=['Tight_junction','Wnt_BBB','Cldn5','Abcg2','Mfsd2a','IFN_antigen_presentation','Vcam1','Cxcl9']
    s=ven[ven.feature.isin(features)].copy(); s['ord']=s.feature.map({f:i for i,f in enumerate(features)}); s=s.sort_values('ord',ascending=False)
    cols=[COLORS['bevac'] if x>0 else COLORS['immune'] for x in s.delta_bevac_minus_acute]
    ax.barh(s.feature.str.replace('_',' '),s.delta_bevac_minus_acute,color=cols,height=.65); ax.axvline(0,color='#909AA0',lw=.7)
    ax.set_xlabel('Bevacizumab − acute EAE expression score'); ax.set_xlim(-.29,.48)
    for i,r in enumerate(s.itertuples()):
        if r.delta_bevac_minus_acute >= 0:
            xpos, ha = r.delta_bevac_minus_acute + .015, 'left'
        else:
            xpos, ha = -.275, 'left'
        ax.text(xpos,i,f"FDR {r.fdr_within_subtype:.3f}",ha=ha,va='center',fontsize=5.9)
    sub=gs[1,1].subgridspec(2,2,height_ratios=[.13,.87],width_ratios=[1.06,.94],hspace=.12,wspace=.55)
    holder=fig.add_subplot(sub[0,:]); holder.axis('off'); header_panel(holder,'D','Published functional and protein evidence',title_x=MAIN_PANEL_TITLE_X)
    ax=fig.add_subplot(sub[1,0])
    endpoints=['Venous endothelial proliferation','IgG leakage','Fibrinogen leakage']
    f=functional.set_index('endpoint').loc[endpoints].reset_index(); y=np.arange(len(f))[::-1]
    for i,r in enumerate(f.itertuples()):
        col=COLORS['bevac'] if r.endpoint.startswith('Venous') else COLORS['barrier']
        ax.plot([r.g_ci_low,r.g_ci_high],[y[i],y[i]],color=col,lw=1.1)
        ax.scatter(r.hedges_g,y[i],s=30,color=col,edgecolor='white',lw=.5,zorder=3)
        ax.text(-2.95,y[i]-.32,f"exact p={r.p_exact_two_sided:.3f}",fontsize=5.1,color='#66737C')
    ax.axvline(0,color='#9DA7AD',lw=.7); ax.set_xlim(-3.05,1.65); ax.set_ylim(-.75,2.55)
    ax.set_yticks(y,['Venous\nproliferation','IgG\nleakage','Fibrinogen\nleakage']); ax.set_xlabel('Hedges’ g (r84 − IgG)'); ax.grid(axis='x',color='#E6EAEC',lw=.5)

    ax=fig.add_subplot(sub[1,1])
    regions=['Lesion MS','NAWM MS','NAWM HC']; markers=['CD31','EGFL7','MCAM']
    mat=np.array([[ihc[(ihc.marker==m)&(ihc.region==r)]['mean'].iloc[0] for r in regions] for m in markers])
    im=ax.imshow(mat,cmap='YlGnBu',vmin=0,vmax=12,aspect='auto')
    for i in range(3):
        for j in range(3): ax.text(j,i,f"{mat[i,j]:.1f}",ha='center',va='center',fontsize=5.8,color='white' if mat[i,j]>7 else COLORS['dark'])
    ax.set_xticks(range(3),['MS\nlesion','MS\nNAWM','Control\nNAWM'],rotation=38,ha='right',fontsize=5.5); ax.set_yticks(range(3),markers)
    ax.set_title('Vessel-positive fragments',fontsize=6.6,pad=4); cb=fig.colorbar(im,ax=ax,fraction=.055,pad=.04); cb.ax.tick_params(labelsize=5.2)
    save(fig,'Figure_5_VEGFA_perturbation')


def supplementary_figures() -> None:
    qc=pd.read_csv(RESULTS/'reconstruction_qc.csv')
    scores=pd.read_csv(RESULTS/'sample_level_module_scores.csv')
    effects=pd.read_csv(RESULTS/'dataset_module_effects.csv')
    stage_all=pd.read_csv(SOURCE/'S10_GSE208747.csv')
    vegf=pd.read_csv(RESULTS/'VEGFA_three_arm_transcriptional_perturbation.csv')

    # S1: raw-droplet filtering and endothelial-cell yield.
    q=qc[qc.dataset=='GSE210776'].dropna(subset=['endothelial_cells']).copy(); q['x']=np.arange(len(q))
    fig,ax=plt.subplots(figsize=(7.4,4.0),layout='constrained'); panel(ax,'A','GSE210776 raw-droplet filtering and endothelial-cell yield')
    ax.bar(q.x,q.barcodes_after_qc,color='#CCD3D7',label='Barcodes after QC'); ax.bar(q.x,q.endothelial_cells,color='#4C93A8',label='Selected endothelial cells'); ax.set_xticks(q.x,q.sample_id,rotation=55,ha='right'); ax.set_ylabel('Number per mouse')
    ax.legend(frameon=False,ncol=2,loc='upper right',borderaxespad=.2)
    save(fig,'Figure_S1_GSE210776_QC',SUPP_OUT)

    # S2: sample-level acute EAE scores across cohorts.
    d=scores[((scores.dataset=='GSE210776')&scores.group.isin(['CFA','Acute']))|((scores.dataset=='GSE199460')&(scores.subset=='All_endothelial'))|((scores.dataset=='GSE95401')&scores.group.isin(['Control','Acute']))].copy()
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.8),layout='constrained')
    for j,mod in enumerate(MODULE_LABEL):
        ax=axs[j]; panel(ax,chr(65+j),MODULE_LABEL[mod].replace('\n',' '));
        for i,ds in enumerate(['GSE210776','GSE199460','GSE95401']):
            b=d[d.dataset==ds]; case='Acute'; ctrl='CFA' if ds=='GSE210776' else 'Control'
            for off,g,col in [(-.13,ctrl,COLORS['control']),(.13,case,COLORS['acute'])]:
                vals=b[b.group==g][mod].dropna().to_numpy(); ax.scatter(np.full(len(vals),i+off),vals,s=25,color=col,edgecolor='white',lw=.4)
                if len(vals): ax.plot([i+off-.08,i+off+.08],[vals.mean(),vals.mean()],color=COLORS['dark'],lw=1)
        ax.axhline(0,color='#D6DCDF',lw=.7); ax.set_xticks(range(3),['210776','199460','95401']); ax.set_xlabel('GSE accession'); ax.grid(axis='y',color='#EDF0F2',lw=.5)
    fig.legend(handles=[Line2D([0],[0],marker='o',ls='',markerfacecolor=COLORS['control'],markeredgecolor='white',label='Control / CFA'),Line2D([0],[0],marker='o',ls='',markerfacecolor=COLORS['acute'],markeredgecolor='white',label='Acute EAE')],loc='lower center',bbox_to_anchor=(.5,-.08),ncol=2,frameon=False)
    save(fig,'Figure_S2_acute_EAE_sample_scores',SUPP_OUT)

    # S3: all-endothelial vs venous-like GSE199460 sensitivity.
    e=effects[effects.dataset.isin(['GSE199460_all','GSE199460_venous'])].copy(); pivot=e.pivot(index='module',columns='dataset',values='hedges_g').loc[list(MODULE_LABEL)]
    fig,ax=plt.subplots(figsize=(5.4,3.3),layout='constrained'); panel(ax,'A','GSE199460 all-endothelial and venous-like sensitivity')
    im=ax.imshow(pivot,cmap='RdBu_r',vmin=-5,vmax=5,aspect='auto');
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]): ax.text(j,i,f"{pivot.iloc[i,j]:.2f}",ha='center',va='center',color='white' if abs(pivot.iloc[i,j])>2.8 else COLORS['dark'])
    ax.set_xticks([0,1],['All endothelial','Venous-like']); ax.set_yticks(range(3),[MODULE_LABEL[m].replace('\n',' ') for m in pivot.index]); fig.colorbar(im,ax=ax,label='Hedges g')
    save(fig,'Figure_S3_GSE199460_subtype_sensitivity',SUPP_OUT)

    # S4: biological-sample distributions underlying the condensed stage panel
    # in main Figure 2E. Groups are independent and are not joined as a trajectory.
    d=scores[scores.dataset=='GSE95401']; order=['Control','Acute','Subacute','Chronic']
    group_colors=[COLORS['control'],COLORS['acute'],COLORS['gold'],COLORS['ecm']]
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.75),layout='constrained')
    rng=np.random.default_rng(20260901)
    for j,mod in enumerate(MODULE_LABEL):
        ax=axs[j]; panel(ax,chr(65+j),MODULE_LABEL[mod].replace('\n',' '))
        for i,(group,color) in enumerate(zip(order,group_colors)):
            vals=d.loc[d.group.eq(group),mod].dropna().to_numpy()
            x=i+rng.uniform(-.075,.075,len(vals))
            ax.scatter(x,vals,s=31,color=color,edgecolor='white',lw=.5,zorder=3)
            ax.plot([i-.14,i+.14],[vals.mean(),vals.mean()],color=COLORS['dark'],lw=1.2,zorder=4)
        ax.axhline(0,color='#D6DCDF',lw=.7)
        ax.set_xticks(range(4),order,rotation=32,ha='right')
        ax.grid(axis='y',color='#EDF0F2',lw=.5)
        if j==0: ax.set_ylabel('Module score')
    fig.text(.5,-.02,'Independent groups; n=3 biological samples per group. Horizontal bars show group means.',ha='center',fontsize=6.4,color='#66737C')
    save(fig,'Figure_S4_GSE95401_time_course',SUPP_OUT)

    # S5: GSE279183 EC selection QC.
    q=qc[qc.dataset=='GSE279183'].dropna(subset=['endothelial_nuclei']).copy(); q=q.sort_values(['group','donor','tissue']); x=np.arange(len(q))
    fig,ax=plt.subplots(figsize=(7.4,4.0),layout='constrained'); panel(ax,'A','GSE279183 endothelial-nucleus selection by tissue')
    cmap={'Control':COLORS['control'],'Chronic_active':COLORS['immune'],'Chronic_inactive':COLORS['ecm']}; ax.bar(x,q.endothelial_nuclei,color=[cmap[g] for g in q.group]); ax.set_xticks(x,q.tissue,rotation=55,ha='right'); ax.set_ylabel('Selected endothelial nuclei')
    ax.legend(handles=[Patch(facecolor=COLORS['control'],label='Control'),Patch(facecolor=COLORS['immune'],label='Chronic active'),Patch(facecolor=COLORS['ecm'],label='Chronic inactive')],frameon=False,ncol=3,loc='upper right',borderaxespad=.2)
    save(fig,'Figure_S5_GSE279183_EC_QC',SUPP_OUT)

    # S6: threshold sensitivity for selected locked modules.
    chosen=['IFN_antigen_presentation','PLVAP_permeability','Caveolae_structural','ECM_protease','Tight_junction','BBB_transport_identity','Wnt_BBB']; contrasts=['Active_vs_NAWM','Mixed_vs_NAWM']; qs=['q20','q30','q40']
    rows=[]; labels=[]
    for m in chosen:
        for c in contrasts:
            rows.append([stage_all[(stage_all.threshold==q)&(stage_all.value=='module_z_ECadj')&(stage_all.module==m)&(stage_all.contrast==c)].delta_median.iloc[0] for q in qs]); labels.append(f"{m.replace('_',' ')} | {'Active' if c.startswith('Active') else 'Mixed'}")
    mat=np.array(rows); fig,ax=plt.subplots(figsize=(6.3,5.0),layout='constrained'); panel(ax,'A','GSE208747 endothelial-rich spot threshold sensitivity')
    im=ax.imshow(mat,cmap='RdBu_r',vmin=-1.6,vmax=1.6,aspect='auto'); ax.set_xticks(range(3),['top 20%','top 30%','top 40%']); ax.set_yticks(range(len(labels)),labels)
    for i in range(mat.shape[0]):
        for j in range(3): ax.text(j,i,f"{mat[i,j]:.2f}",ha='center',va='center',fontsize=5.7,color='white' if abs(mat[i,j])>.9 else COLORS['dark'])
    fig.colorbar(im,ax=ax,label='EC-adjusted median difference')
    save(fig,'Figure_S6_GSE208747_threshold_sensitivity',SUPP_OUT)

    # S7: three-arm VEGF module vectors.
    v=vegf.set_index('module').loc[list(MODULE_LABEL)+[m for m in vegf.module if m not in MODULE_LABEL]]; mat=v[['acute_minus_CFA','bevacizumab_minus_acute','bevacizumab_minus_CFA']].to_numpy()
    fig,ax=plt.subplots(figsize=(6.2,5.3),layout='constrained'); panel(ax,'A','GSE210776 three-arm transcriptional perturbation across all modules')
    im=ax.imshow(mat,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto'); ax.set_xticks(range(3),['Acute − CFA','BEVAC − acute','BEVAC − CFA']); ax.set_yticks(range(len(v)),[x.replace('_',' ') for x in v.index])
    for i in range(mat.shape[0]):
        for j in range(3): ax.text(j,i,f"{mat[i,j]:.2f}",ha='center',va='center',fontsize=5.5,color='white' if abs(mat[i,j])>1.15 else COLORS['dark'])
    fig.colorbar(im,ax=ax,label='Module-score difference')
    save(fig,'Figure_S7_VEGFA_all_modules',SUPP_OUT)

def main() -> None:
    set_style()
    fig1_overview()
    fig2_meta()
    fig3_human_ms()
    fig4_cross_species()
    fig5_vegf()
    supplementary_figures()
    print(f"Wrote figures to {OUT}")


if __name__ == '__main__':
    main()
