#!/usr/bin/env python3
"""Create a consistent BMC Genomics-ready figure set from frozen and raw-audit tables.

All stochastic operations are disabled. The fixed project seed is recorded for
reproducibility and future extensions that may introduce stochastic procedures.
"""

from __future__ import annotations

from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "analysis_results" / "source_tables"
RAW = ROOT / "analysis_results"
OUT = ROOT / "outputs" / "bmc_submission" / "figures"
MAIN = OUT / "main"
SUPP = OUT / "supplementary"
PREV = OUT / "previews"
for d in [MAIN, SUPP, PREV]:
    d.mkdir(parents=True, exist_ok=True)

SEED = 20260820
NEG, POS, MID = "#2166AC", "#B2182B", "#F2F2F2"
NAVY, TEAL, GOLD, PURPLE, GRAY = "#2C3E50", "#1B9E77", "#D99A00", "#7570B3", "#6B7280"
GROUP_COLORS = {
    "Control white matter": "#7F8C8D", "NAWM": "#4C78A8",
    "Active lesion": "#B2182B", "Mixed active/inactive lesion": "#E08214",
}

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.titlesize": 9.0, "axes.titleweight": "semibold",
    "axes.labelsize": 7.5, "xtick.labelsize": 6.6, "ytick.labelsize": 6.6,
    "legend.fontsize": 6.4, "axes.linewidth": 0.7,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(TAB / f"{name}.csv")


def clean(name: str) -> str:
    short = {
        "Adhesion_trafficking": "Adhesion/trafficking",
        "BBB_transport_identity": "BBB transport",
        "Caveolae_structural": "Caveolae",
        "ECM_protease": "ECM/protease",
        "IFN_antigen_presentation": "IFN/antigen",
        "PLVAP_permeability": "PLVAP/permeability",
        "Tight_junction": "Tight junction",
        "VEGF_response": "VEGF response",
        "Wnt_BBB": "Wnt–BBB",
        "Myeloid_proteolytic": "Myeloid/proteolytic",
    }
    if name in short:
        return short[name]
    return (name.replace("_", " ").replace("BBB", "BBB").replace("Wnt", "Wnt")
            .replace("IFN antigen", "IFN/antigen").replace("PLVAP permeability", "PLVAP/permeability")
            .replace("TNF NFkB", "TNF/NF-κB").replace("IL1 IL6", "IL-1/IL-6")
            .replace("TGF remodeling", "TGF/remodeling").replace("VEGF angiogenic", "VEGF/angiogenic")
            .replace("Myeloid proteolytic", "Myeloid/proteolytic"))


def panel(ax, letter: str, title: str | None = None, letter_x=-0.11, letter_y=1.07):
    ax.text(letter_x, letter_y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", clip_on=False)
    if title:
        ax.set_title(textwrap.fill(title, width=34), loc="left", pad=6)


def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def effect_bar(ax, labels, effects, pvals=None, fdr=None, title=None, xlabel="Effect",
               sig_offset=None, sig_columns=None):
    labels = list(labels); effects = np.asarray(effects, float)
    y = np.arange(len(labels))
    ax.barh(y, effects, color=np.where(effects >= 0, POS, NEG), edgecolor="none", alpha=.92)
    ax.axvline(0, color="#333333", lw=.7)
    ax.set_yticks(y, [clean(x) for x in labels]); ax.invert_yaxis(); ax.set_xlabel(xlabel)
    if title: ax.set_title(title, loc="left", pad=6)
    despine(ax)
    span = max(np.max(np.abs(effects)), .2)
    offset = span*.035 if sig_offset is None else float(sig_offset)
    if pvals is not None:
        for i, (v, p) in enumerate(zip(effects, pvals)):
            q = None if fdr is None else float(np.asarray(fdr)[i])
            mark = "†" if q is not None and q < .1 else ("*" if float(p) < .05 else "")
            if q is not None and q < .05: mark = "**"
            if mark:
                if sig_columns is not None:
                    x_mark = sig_columns[1] if v >= 0 else sig_columns[0]
                    ha = "center"
                else:
                    x_mark = v + np.sign(v if v else 1)*offset
                    ha = "left" if v >= 0 else "right"
                ax.text(x_mark, i, mark, va="center",ha=ha,fontsize=7,
                        fontweight="bold",clip_on=False)


def heatmap(ax, frame: pd.DataFrame, title=None, cbar_label="Standardized effect", annotate=True,
            vlim=None, xrot=25):
    vals = frame.to_numpy(float)
    lim = vlim or max(np.nanmax(np.abs(vals)), .5)
    im = ax.imshow(vals, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), aspect="auto")
    ax.set_xticks(range(frame.shape[1]), frame.columns, rotation=xrot, ha="right")
    ax.set_yticks(range(frame.shape[0]), [clean(x) for x in frame.index])
    if title: ax.set_title(title, loc="left", pad=6)
    if annotate:
        for i in range(frame.shape[0]):
            for j in range(frame.shape[1]):
                v = vals[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.8,
                        color="white" if abs(v) > lim*.62 else "#222222")
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.03)
    cb.set_label(cbar_label, fontsize=6.5); cb.ax.tick_params(labelsize=5.8)
    return im


def save(fig, stem: str, dest: Path):
    fig.subplots_adjust(left=.14, right=.965, top=.90, bottom=.12)
    png = PREV / f"{stem}.png"
    tmp = dest / f".{stem}_tmp.png"
    fig.savefig(tmp, dpi=600, facecolor="white")
    plt.close(fig)
    with Image.open(tmp) as im:
        rgb = im.convert("RGB")
        rgb.save(dest / f"{stem}.tiff", format="TIFF", compression="tiff_lzw", dpi=(600, 600))
        rgb.resize((max(1, rgb.width//4), max(1, rgb.height//4)), Image.Resampling.LANCZOS).save(png, dpi=(150,150))
    tmp.unlink()


def box(ax, xy, wh, text, color="#FFFFFF", edge=NAVY, fontsize=5.8):
    x,y=xy; w,h=wh
    p=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012,rounding_size=0.018",
                     transform=ax.transAxes, facecolor=color, edgecolor=edge, linewidth=1)
    ax.add_patch(p); ax.text(x+w/2,y+h/2,text,ha="center",va="center",transform=ax.transAxes,fontsize=fontsize)


def square_box(ax, xy, wh, text, color="#FFFFFF", edge=NAVY, fontsize=5.0):
    """Draw a closed rectangular workflow box with an inset safety margin."""
    x,y=xy; w,h=wh
    p=Rectangle((x,y),w,h,transform=ax.transAxes,facecolor=color,
                edgecolor=edge,linewidth=1.2,clip_on=False)
    ax.add_patch(p)
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",transform=ax.transAxes,
            fontsize=fontsize,linespacing=1.15)


def arrow(ax, a, b):
    ax.add_patch(FancyArrowPatch(a,b,transform=ax.transAxes,arrowstyle="-|>",mutation_scale=10,
                                 lw=1,color=NAVY))


def fig1():
    s3=load("S3_EAE_modules"); s6=load("S6_GSE95401")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.58,wspace=.62)
    ax=fig.add_subplot(gs[0,0]); ax.axis("off"); panel(ax,"A","Cross-platform study design")
    box(ax,(.10,.56),(.22,.24),"Mouse discovery\nGSE210776",color="#EAF2F8",fontsize=4.5)
    box(ax,(.39,.56),(.22,.24),"Mouse replication\nGSE199460\nGSE95401",color="#E8F6F3",fontsize=4.5)
    box(ax,(.68,.56),(.22,.24),"Human validation\nGSE279183\nGSE208747",color="#FDEDEC",fontsize=4.5)
    arrow(ax,(.32,.68),(.39,.68)); arrow(ax,(.61,.68),(.68,.68))
    box(ax,(.22,.06),(.56,.30),"Biological replicate\nis the inferential unit\nCells/spots are aggregated\nwithin each sample",color="#FBFCFC",fontsize=4.4)
    arrow(ax,(.50,.56),(.50,.36))
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Acute EAE endothelial programs (n=4 vs 3 mice)")
    effect_bar(ax,s3.Feature,s3.Delta,s3.p,s3.FDR,
               xlabel="Median score difference (EAE − CFA)",
               sig_offset=.045)
    ax.set_xlim(-.55,.50)
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Temporal endothelial support (acute stage)")
    genes=["Cxcl9","Esm1","Mfsd2a","Cav1","Cav2","Axin2"]
    d=s6[(s6.Stage=="Acute") & (s6.Gene.isin(genes))].set_index("Gene").reindex(genes).dropna()
    effect_bar(ax,d.index,d.log2FC,d.p,d.FDR,
               xlabel="log2 fold change vs control",
               sig_columns=(-1.10,7.80))
    ax.set_xlim(-2.00,8.40)
    ax=fig.add_subplot(gs[1,1]); ax.axis("off"); panel(ax,"D","Working model")
    box(ax,(.03,.63),(.28,.22),"Inflammatory\ncues",color="#FDEDEC",fontsize=5.0)
    box(ax,(.38,.63),(.24,.22),"Venous EC\nstate shift",color="#F5EEF8")
    box(ax,(.69,.63),(.28,.22),"Barrier\nremodeling",color="#EAF2F8")
    arrow(ax,(.31,.74),(.38,.74)); arrow(ax,(.62,.74),(.69,.74))
    box(ax,(.17,.16),(.66,.26),"Conserved direction is emphasized;\ncross-platform raw values\nare not pooled",color="#FBFCFC",fontsize=5.0)
    arrow(ax,(.50,.63),(.50,.41))
    save(fig,"Figure1_Study_design_and_acute_EAE",MAIN)


def fig2():
    s4=load("S4_EAE_genes"); s5=load("S5_Mouse_sender"); s7=load("S7_BEVAC3")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.58,wspace=.62)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Selected acute-EAE gene effects")
    genes=["H2-Aa","H2-Ab1","Cxcl9","Cxcl10","Cav1","Cav2","Mfsd2a","Cldn5","Plvap","Vegfa"]
    d=s4[s4.Feature.isin(genes)].set_index("Feature").reindex(genes).dropna()
    effect_bar(ax,d.index,d.Delta,d.p,d.FDR,xlabel="Expression-score difference")
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Candidate sender-cell signals")
    hm=s5.set_index("CellType").rename_axis(None)
    heatmap(ax,hm,"",cbar_label="Mean expression shift",annotate=False,xrot=35)
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","VEGF-A perturbation in venous ECs (n=4 vs 3 mice)")
    feats=["Tight_junction","Wnt_BBB","BBB_transport_identity","Caveolae_structural","VEGF_response","IFN_antigen_presentation"]
    d=s7[(s7.subtype=="Venous") & (s7.feature.isin(feats))].set_index("feature").reindex(feats).dropna()
    effect_bar(ax,d.index,d.delta_bevac_minus_acute,d.p_exact,d.fdr_within_subtype,
               xlabel="Score difference (bevacizumab − acute EAE)")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Selected gene shifts after VEGF-A blockade")
    feats=["Cldn5","Abcg2","Kdr","Cav1","Mfsd2a","Plvap","Angpt2","Cxcl9","H2-Aa","H2-Ab1"]
    d=s7[(s7.subtype=="Venous") & (s7.feature.isin(feats))].set_index("feature").reindex(feats).dropna()
    effect_bar(ax,d.index,d.delta_bevac_minus_acute,d.p_exact,d.fdr_within_subtype,
               xlabel="Expression-score difference")
    save(fig,"Figure2_Immune_reprogramming_and_VEGFA",MAIN)


def fig3():
    s8=load("S8_GSE279183_snRNA"); s9=load("S9_GSE279183_spatial")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.60,wspace=.62)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Human MS endothelial nuclei (6 MS vs 6 controls)")
    d=s8[s8.contrast=="MS_vs_CTRL"]
    effect_bar(ax,d.metric_label,d.delta,d.p,d.FDR,xlabel="Mean score difference (MS − control)")
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Vascular-inflammation niches vs paired PPWM (n=7 donors)")
    effect_bar(ax,s9.pathway,s9.delta_VI_vs_PPWM,s9.p_VI_vs_PPWM,s9.FDR_VI_vs_PPWM,
               xlabel="Paired score difference")
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Vascular-inflammation niches vs control WM (7 vs 5 donors)")
    effect_bar(ax,s9.pathway,s9.delta_VI_vs_CTRL,s9.p_VI_vs_CTRL,s9.FDR_VI_vs_CTRL,
               xlabel="Donor-level score difference")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Evidence synthesis")
    frame=pd.DataFrame({"snRNA MS−control":d.set_index("metric_label").delta,
                        "Spatial VI−PPWM":s9.set_index("pathway").delta_VI_vs_PPWM})
    common={"IFN/antigen presentation":"JAK-STAT","TGF response":"TGFb","VEGF response":"VEGF"}
    rows=[]
    for a,b in common.items(): rows.append([d.set_index("metric_label").delta.get(a,np.nan),s9.set_index("pathway").delta_VI_vs_PPWM.get(b,np.nan)])
    h=pd.DataFrame(rows,index=list(common.keys()),columns=["snRNA MS−control","Spatial VI−PPWM"])
    heatmap(ax,h,cbar_label="Within-dataset effect",annotate=True,xrot=20)
    save(fig,"Figure3_Human_MS_vascular_niches",MAIN)


def fig4():
    q=RAW/"gse208747"/"uploaded_archive_qc.csv"; qc=pd.read_csv(q); s10=load("S10_GSE208747")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.62,wspace=.68)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","GSE208747 spatial dataset: 15 sections / 14 donors")
    order=["Control white matter","NAWM","Active lesion","Mixed active/inactive lesion"]
    counts=qc.groupby("group").size().reindex(order)
    ax.bar(range(4),counts,color=[GROUP_COLORS[x] for x in order]); ax.set_xticks(range(4),["Control WM","NAWM","Active","Mixed"],rotation=25,ha="right"); ax.set_ylabel("Sections"); despine(ax)
    for i,v in enumerate(counts): ax.text(i,v+.08,str(int(v)),ha="center",fontsize=7)
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Per-section sequencing QC")
    for g,dd in qc.groupby("group"):
        ax.scatter(dd.median_genes_per_spot,dd.median_UMI_per_spot,s=32,label=g,color=GROUP_COLORS[g],edgecolor="white",linewidth=.4)
        for _,r in dd.iterrows(): ax.text(r.median_genes_per_spot,r.median_UMI_per_spot,r["sample"],fontsize=5.5,ha="left",va="bottom")
    ax.set_xlabel("Median genes per tissue spot"); ax.set_ylabel("Median UMI per tissue spot"); despine(ax); ax.legend(frameon=False,loc="upper left")
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","EC-adjusted lesion vs non-lesional effects (q30)")
    d=s10[(s10.threshold=="q30")&(s10.value=="module_z_ECadj")&(s10.contrast=="Lesion_vs_nonlesional")]
    effect_bar(ax,d.module,d.delta_median,d.p_exact,d.FDR_BH,xlabel="Median donor/section-level difference")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Threshold sensitivity of key pathways")
    d=s10[(s10.threshold.isin(["q20","q30","q40"]))&(s10.value=="module_z_ECadj")&(s10.contrast=="Lesion_vs_nonlesional")]
    keys=["IFN_antigen_presentation","TGF_response","ECM_protease","Tight_junction","BBB_transport_identity","VEGF_response"]
    h=d[d.module.isin(keys)].pivot(index="module",columns="threshold",values="delta_median").reindex(keys)[["q20","q30","q40"]]
    h.columns=["top 20%","top 30%","top 40%"]
    heatmap(ax,h,cbar_label="Median difference",annotate=True,xrot=20)
    save(fig,"Figure4_Independent_MS_spatial_validation",MAIN)


def fig5():
    s11=load("S11_GSE163005"); s12=load("S12_GSE163194"); s13=load("S13_Cross_disease")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.62,wspace=.68)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Viral encephalitis specificity (patient-level)")
    d=s11[s11.contrast.isin(["VE vs IIH","VE vs MS"])]
    h=d.pivot(index="module",columns="contrast",values="delta_median")
    heatmap(ax,h,cbar_label="Median module difference",annotate=True,xrot=20)
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Bacterial meningitis: acute vs recovery (4 vs 9 patients)")
    d=s12[(s12.celltype=="WholeCSF")&(s12.contrast=="Acute_vs_All_recovery")]
    effect_bar(ax,d.module,d.delta_median,d.p_MWU,d.FDR_BH,xlabel="Median module difference")
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Cross-disease effect map")
    h=s13.pivot(index="module",columns="comparison",values="standardized_effect")
    heatmap(ax,h,cbar_label="Standardized effect",annotate=True,xrot=20)
    ax=fig.add_subplot(gs[1,1]); ax.axis("off"); panel(ax,"D","Specificity model")
    box(ax,(.02,.66),(.30,.22),"MS\nlocalized vascular\nremodeling",color="#EAF2F8",fontsize=4.6)
    box(ax,(.35,.66),(.30,.22),"Viral encephalitis\nantiviral\nresponse",color="#FDEDEC",fontsize=4.6)
    box(ax,(.68,.66),(.30,.22),"Bacterial meningitis\nIL-1/chemokine\nresponse",color="#FFF4E6",fontsize=4.6)
    box(ax,(.12,.14),(.76,.28),"Shared genes do not imply\nidentical disease programs;\ncompare sample-level effects\nacross cohorts",color="#FBFCFC",fontsize=4.6)
    arrow(ax,(.17,.66),(.37,.42)); arrow(ax,(.50,.66),(.50,.42)); arrow(ax,(.83,.66),(.63,.42))
    save(fig,"Figure5_Cross_disease_specificity_model",MAIN)


def sfig1():
    sig=load("S2_Locked_signatures"); ds=load("S1_Datasets")
    fig=plt.figure(figsize=(7.5,5.8)); gs=fig.add_gridspec(2,2,hspace=.48,wspace=.35)
    ax=fig.add_subplot(gs[:,0]); ax.axis("off"); panel(ax,"A","Predefined BBB / vascular programs")
    for i,(_,r) in enumerate(sig.iterrows()):
        col=i%2; row=i//2; x=.02+col*.49; y=.84-row*.16
        genes = [g.strip() for g in str(r.Genes).split(";") if g.strip()]
        shown = "; ".join(genes[:5]) + ("; …" if len(genes) > 5 else "")
        box(ax,(x,y),(.45,.12),f"{clean(r.Module)}\n{textwrap.fill(shown, 30)}",color="#FBFCFC",fontsize=4.8)
    ax=fig.add_subplot(gs[0,1]); ax.axis("off"); panel(ax,"B","Inference principles")
    box(ax,(.02,.60),(.28,.25),"Biological\nreplicate",color="#EAF2F8")
    box(ax,(.36,.60),(.28,.25),"Within-sample\naggregation",color="#E8F6F3")
    box(ax,(.70,.60),(.28,.25),"Sample-level\ntests",color="#FDEDEC")
    arrow(ax,(.30,.72),(.36,.72)); arrow(ax,(.64,.72),(.70,.72))
    box(ax,(.18,.15),(.64,.22),"Cross-dataset synthesis uses\ndirection + standardized effect",color="#FBFCFC")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"C","Dataset roles and inferential units")
    y=np.arange(len(ds)); ax.barh(y,np.ones(len(ds)),color=[TEAL if "Mouse" in x else PURPLE for x in ds["Species/tissue"]]); ax.set_xlim(0,1); ax.set_yticks(y,ds.Dataset); ax.invert_yaxis(); ax.set_xticks([])
    for i,r in ds.iterrows(): ax.text(.03,i,str(r["Primary statistical unit"]),va="center",color="white",fontweight="bold",fontsize=6)
    despine(ax)
    save(fig,"Figure_S1_framework_and_inference",SUPP)


def sfig2():
    s3=load("S3_EAE_modules"); s6=load("S6_GSE95401")
    fig,axs=plt.subplots(2,2,figsize=(7.5,5.6));
    panel(axs[0,0],"A","Discovery module effects (GSE210776)"); effect_bar(axs[0,0],s3.Feature,s3.Delta,s3.p,s3.FDR,xlabel="EAE − CFA",sig_offset=.045); axs[0,0].set_xlim(-.55,.50)
    panel(axs[0,1],"B","GSE95401 selected genes across stages")
    genes=["Cxcl9","Esm1","Cav1","Cav2","Mfsd2a","Axin2"]; h=s6[s6.Gene.isin(genes)].pivot(index="Gene",columns="Stage",values="log2FC").reindex(genes)
    heatmap(axs[0,1],h,cbar_label="log2 fold change",annotate=True,xrot=20)
    panel(axs[1,0],"C","Acute-stage directional support")
    d=s6[(s6.Stage=="Acute")&(s6.Gene.isin(genes))].set_index("Gene").reindex(genes).dropna(); effect_bar(axs[1,0],d.index,d.log2FC,d.p,d.FDR,xlabel="log2 fold change",sig_columns=(-1.18,7.72)); axs[1,0].set_xlim(-1.60,8.25)
    panel(axs[1,1],"D","Interpretation boundary"); axs[1,1].axis("off")
    box(axs[1,1],(.08,.62),(.84,.22),"Independent datasets support\nendothelial state direction",color="#E8F6F3")
    box(axs[1,1],(.08,.22),(.84,.22),"No direct pooling of raw expression\nacross platforms or studies",color="#FBFCFC")
    fig.subplots_adjust(hspace=.62, wspace=.62)
    save(fig,"Figure_S2_mouse_replication_and_temporal_support",SUPP)


def sfig3():
    s7=load("S7_BEVAC3")
    fig,axs=plt.subplots(2,2,figsize=(7.5,5.6))
    axs[0,0].axis("off"); panel(axs[0,0],"A","VEGF-A perturbation design")
    box(axs[0,0],(.05,.55),(.32,.25),"Acute EAE\nn=4 mice",color="#FDEDEC"); box(axs[0,0],(.63,.55),(.32,.25),"Bevacizumab EAE\nn=3 mice",color="#E8F6F3"); arrow(axs[0,0],(.37,.67),(.63,.67))
    box(axs[0,0],(.12,.10),(.76,.28),"Primary question:\ndoes VEGF-A blockade shift\nvenous endothelial transcription\ntoward homeostasis?",color="#FBFCFC",fontsize=4.9)
    d=s7[s7.subtype=="Venous"]
    feats=["Tight_junction","BBB_transport_identity","Caveolae_structural","PLVAP_permeability","IFN_antigen_presentation","Wnt_BBB","VEGF_response","ECM_protease"]
    x=d[d.feature.isin(feats)].set_index("feature").reindex(feats).dropna(); panel(axs[0,1],"B","Venous program effects"); effect_bar(axs[0,1],x.index,x.delta_bevac_minus_acute,x.p_exact,x.fdr_within_subtype,xlabel="Bevacizumab − acute EAE")
    feats=["Cldn5","Abcg2","Kdr","Cav1","Cav2","Mfsd2a","Plvap","Angpt2","Cxcl9","Cxcl10","H2-Aa","H2-Ab1"]
    x=d[d.feature.isin(feats)].set_index("feature").reindex(feats).dropna(); panel(axs[1,0],"C","Selected gene effects"); effect_bar(axs[1,0],x.index,x.delta_bevac_minus_acute,x.p_exact,x.fdr_within_subtype,xlabel="Bevacizumab − acute EAE")
    axs[1,1].axis("off"); panel(axs[1,1],"D","Interpretation boundary")
    box(axs[1,1],(.07,.62),(.86,.22),"Supported: partial/directional transcriptomic\nnormalization of selected BBB programs",color="#E8F6F3")
    box(axs[1,1],(.07,.20),(.86,.22),"Not demonstrated: functional BBB restoration,\ntracer rescue, or transport normalization",color="#FBFCFC")
    fig.subplots_adjust(hspace=.62, wspace=.62)
    save(fig,"Figure_S3_VEGFA_perturbation_details",SUPP)


def sfig4():
    qc=pd.read_csv(RAW/"gse279183"/"sample_qc_summary.csv"); scores=pd.read_csv(RAW/"gse279183"/"spot_qc_and_module_scores.csv.gz")
    fig,axs=plt.subplots(2,2,figsize=(7.5,5.6))
    panel(axs[0,0],"A","Uploaded GSE279183 spatial matrices")
    axs[0,0].bar(qc["sample"],qc["spots_under_tissue"],color=[TEAL,PURPLE,GOLD]); axs[0,0].set_ylabel("Tissue spots"); axs[0,0].tick_params(axis='x',rotation=20); despine(axs[0,0])
    panel(axs[0,1],"B","Per-sample spot QC",letter_y=1.14)
    axs[0,1].scatter(qc.median_genes,qc.median_UMI,s=45,color=[TEAL,PURPLE,GOLD])
    for _,r in qc.iterrows(): axs[0,1].text(r.median_genes,r.median_UMI,r["sample"],fontsize=6,va="bottom")
    axs[0,1].set_xlabel("Median genes"); axs[0,1].set_ylabel("Median UMI"); despine(axs[0,1])
    panel(axs[1,0],"C","Spot-level BBB identity vs inflammation")
    for sample,dd in scores.groupby("sample"):
        axs[1,0].scatter(dd.BBB_transport_identity,dd.IFN_antigen_presentation,s=2,alpha=.18,label=sample)
    axs[1,0].set_xlabel("BBB transport identity score"); axs[1,0].set_ylabel("IFN/antigen score"); despine(axs[1,0]); axs[1,0].legend(markerscale=3,frameon=False)
    panel(axs[1,1],"D","Per-sample module medians")
    cols=["IFN_antigen_presentation","ROS_Src","ECM_protease"]
    h=scores.groupby("sample")[cols].median(); h=(h-h.mean())/h.std(ddof=0)
    h.columns=["IFN/antigen","ROS/Src","ECM/protease"]
    heatmap(axs[1,1],h,cbar_label="Across-sample z score",annotate=True,xrot=18)
    fig.subplots_adjust(hspace=.62, wspace=.62)
    save(fig,"Figure_S4_GSE279183_QC_and_robustness",SUPP)


def sfig5():
    qc=pd.read_csv(RAW/"gse208747"/"uploaded_archive_qc.csv"); s10=load("S10_GSE208747")
    fig=plt.figure(figsize=(7.5,8.2)); gs=fig.add_gridspec(3,2,hspace=.85,wspace=.62,height_ratios=[1,1.05,1.05])
    axa=fig.add_subplot(gs[0,0]); panel(axa,"A","Tissue spots per section",letter_y=1.20)
    c=[GROUP_COLORS[g] for g in qc.group]; axa.bar(qc["sample"],qc.spots_under_tissue,color=c); axa.tick_params(axis='x',rotation=45); axa.set_ylabel("Spots under tissue"); despine(axa)
    axb=fig.add_subplot(gs[0,1]); panel(axb,"B","Median molecular complexity",letter_y=1.20)
    for g,d in qc.groupby("group"): axb.scatter(d.median_genes_per_spot,d.median_UMI_per_spot,s=30,label=g,color=GROUP_COLORS[g])
    axb.set_xlabel("Median genes per spot"); axb.set_ylabel("Median UMI per spot"); despine(axb); axb.legend(frameon=False)
    axc=fig.add_subplot(gs[1,:]); panel(axc,"C","Threshold sensitivity: EC-adjusted lesion vs non-lesional",letter_y=1.20)
    d=s10[(s10.value=="module_z_ECadj")&(s10.contrast=="Lesion_vs_nonlesional")&(s10.threshold.isin(["q20","q30","q40"]))]
    h=d.pivot(index="module",columns="threshold",values="delta_median")[["q20","q30","q40"]]; h.columns=["top 20%","top 30%","top 40%"]
    heatmap(axc,h,cbar_label="Median difference",annotate=True,xrot=0)
    axd=fig.add_subplot(gs[2,:]); panel(axd,"D","Exact P and BH-FDR at the primary q30 threshold",letter_y=1.20)
    d=s10[(s10.value=="module_z_ECadj")&(s10.contrast=="Lesion_vs_nonlesional")&(s10.threshold=="q30")].sort_values("FDR_BH")
    y=np.arange(len(d)); axd.scatter(d.p_exact,y,label="Exact P",color=NAVY,s=20); axd.scatter(d.FDR_BH,y,label="BH-FDR",color=POS,s=20,marker="s"); axd.axvline(.05,color=GRAY,lw=.7,ls="--"); axd.set_xscale("log"); axd.set_yticks(y,[clean(x) for x in d.module]); axd.invert_yaxis(); axd.set_xlabel("Probability (log scale)"); despine(axd); axd.legend(frameon=False,loc="upper right",bbox_to_anchor=(1.0,1.02))
    save(fig,"Figure_S5_GSE208747_QC_and_sensitivity",SUPP)


def sfig6():
    qc=pd.read_csv(RAW/"gse163005"/"qc_summary.csv"); recon=pd.read_csv(RAW/"gse163005"/"raw_reconstructed_module_contrasts.csv"); audit=pd.read_csv(RAW/"raw_vs_frozen_direction_audit.csv"); comp=pd.read_csv(RAW/"gse163005"/"patient_cell_composition.csv")
    fig,axs=plt.subplots(2,2,figsize=(7.5,6.2))
    panel(axs[0,0],"A","Cells and patients by diagnosis"); x=np.arange(len(qc)); axs[0,0].bar(x-.18,qc.cells/1000,.36,label="Cells (×10³)",color=TEAL); axs[0,0].bar(x+.18,qc.patients,.36,label="Patients",color=NAVY); axs[0,0].set_xticks(x,qc.diagnosis); axs[0,0].legend(frameon=False); despine(axs[0,0])
    panel(axs[0,1],"B","Patient-level cell composition")
    h=comp.groupby(["diagnosis","celltype"]).fraction.median().unstack(fill_value=0); bottom=np.zeros(len(h))
    for ct in h.columns: axs[0,1].bar(h.index,h[ct],bottom=bottom,label=ct); bottom+=h[ct].values
    axs[0,1].set_ylabel("Median fraction"); despine(axs[0,1]); axs[0,1].legend(frameon=False,ncol=4,fontsize=4.6,loc="upper center",bbox_to_anchor=(.5,-.15),columnspacing=.8,handlelength=1.5)
    panel(axs[1,0],"C","Raw reconstruction: VE vs MS")
    d=recon[recon.contrast=="VE vs MS"]; effect_bar(axs[1,0],d.module,d.delta_median,d.p_exact,d.FDR_BH,xlabel="Median patient-level difference")
    panel(axs[1,1],"D","Frozen vs reconstructed effects")
    d=audit[audit.dataset=="GSE163005"]; colors=np.where(d.direction_concordant,TEAL,POS); axs[1,1].scatter(d.effect_frozen,d.effect_reconstructed,c=colors,s=22,alpha=.8); lo=min(d.effect_frozen.min(),d.effect_reconstructed.min()); hi=max(d.effect_frozen.max(),d.effect_reconstructed.max()); axs[1,1].plot([lo,hi],[lo,hi],color=GRAY,lw=.7,ls="--"); axs[1,1].axhline(0,color=GRAY,lw=.5); axs[1,1].axvline(0,color=GRAY,lw=.5); axs[1,1].set_xlabel("Frozen source effect"); axs[1,1].set_ylabel("Raw reconstructed effect"); despine(axs[1,1])
    legend_handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,markeredgecolor=TEAL,label="Direction concordant",markersize=4.5),Line2D([0],[0],marker="o",color="none",markerfacecolor=POS,markeredgecolor=POS,label="Direction discordant",markersize=4.5)]
    axs[1,1].legend(handles=legend_handles,frameon=False,loc="lower right",fontsize=5.0)
    fig.subplots_adjust(hspace=.92, wspace=.62)
    save(fig,"Figure_S6_GSE163005_specificity_details",SUPP)


def sfig7():
    lib=pd.read_csv(RAW/"gse163194"/"library_qc.csv"); recon=pd.read_csv(RAW/"gse163194"/"raw_reconstructed_module_contrasts.csv"); audit=pd.read_csv(RAW/"raw_vs_frozen_direction_audit.csv"); comp=pd.read_csv(RAW/"gse163194"/"marker_reconstructed_composition.csv")
    fig,axs=plt.subplots(2,2,figsize=(7.5,6.2))
    panel(axs[0,0],"A","Libraries collapsed to biological samples"); g=lib.groupby("patient").agg(libraries=("library","count"),cells=("cells","sum")); axs[0,0].scatter(g.libraries,g.cells,s=22,color=TEAL); axs[0,0].set_xlabel("Libraries per patient"); axs[0,0].set_ylabel("Cells per patient"); despine(axs[0,0])
    panel(axs[0,1],"B","Marker-reconstructed composition")
    h=comp[comp.group.isin(["Acute","Recovery"])].groupby(["group","celltype"]).fraction.median().unstack(fill_value=0); bottom=np.zeros(len(h))
    for ct in h.columns: axs[0,1].bar(h.index,h[ct],bottom=bottom,label=ct); bottom+=h[ct].values
    axs[0,1].set_ylabel("Median fraction"); despine(axs[0,1]); axs[0,1].legend(frameon=False,ncol=5,fontsize=4.3,loc="upper center",bbox_to_anchor=(.5,-.15),columnspacing=.7,handlelength=1.4)
    panel(axs[1,0],"C","Raw reconstruction: acute vs recovery")
    effect_bar(axs[1,0],recon.module,recon.delta_median,recon.p_MWU,recon.FDR_BH,xlabel="Median patient-level difference")
    panel(axs[1,1],"D","Frozen vs reconstructed effects")
    d=audit[audit.dataset=="GSE163194"]; colors=np.where(d.direction_concordant,TEAL,POS); axs[1,1].scatter(d.effect_frozen,d.effect_reconstructed,c=colors,s=30); lo=min(d.effect_frozen.min(),d.effect_reconstructed.min()); hi=max(d.effect_frozen.max(),d.effect_reconstructed.max()); axs[1,1].plot([lo,hi],[lo,hi],color=GRAY,lw=.7,ls="--"); axs[1,1].axhline(0,color=GRAY,lw=.5); axs[1,1].axvline(0,color=GRAY,lw=.5); axs[1,1].set_xlabel("Frozen source effect"); axs[1,1].set_ylabel("Raw reconstructed effect"); despine(axs[1,1])
    legend_handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,markeredgecolor=TEAL,label="Direction concordant",markersize=4.5),Line2D([0],[0],marker="o",color="none",markerfacecolor=POS,markeredgecolor=POS,label="Direction discordant",markersize=4.5)]
    axs[1,1].legend(handles=legend_handles,frameon=False,loc="lower right",fontsize=5.0)
    fig.subplots_adjust(hspace=.92, wspace=.62)
    save(fig,"Figure_S7_GSE163194_bacterial_details",SUPP)


def sfig8():
    s13=load("S13_Cross_disease"); s3=load("S3_EAE_modules"); s8=load("S8_GSE279183_snRNA")
    fig,axs=plt.subplots(2,1,figsize=(7.5,7.2))
    panel(axs[0],"A","Cross-disease evidence map")
    h=s13.pivot(index="module",columns="comparison",values="standardized_effect"); heatmap(axs[0],h,cbar_label="Standardized effect",annotate=True,xrot=25)
    panel(axs[1],"B","BBB-focused cross-cohort map")
    map_eae=s3.set_index("Feature").Delta; ms=s8[s8.contrast=="MS_vs_CTRL"].set_index("metric").delta
    feats=["Tight_junction","BBB_transport_identity","Caveolae_structural","PLVAP_permeability","IFN_antigen_presentation","Wnt_BBB","TGF_response","VEGF_response","ECM_protease"]
    hh=pd.DataFrame({"Acute EAE":map_eae.reindex(feats),"Human MS EC":ms.reindex(feats)},index=feats)
    heatmap(axs[1],hh,cbar_label="Within-cohort effect",annotate=True,xrot=25)
    fig.subplots_adjust(hspace=.62)
    save(fig,"Figure_S8_cross_cohort_evidence_maps",SUPP)


def main():
    np.random.seed(SEED)
    for fn in [fig1,fig2,fig3,fig4,fig5,sfig1,sfig2,sfig3,sfig4,sfig5,sfig6,sfig7,sfig8]:
        fn()
    print(f"Created 13 TIFF figures at 600 dpi in {OUT}")


if __name__ == "__main__":
    main()
