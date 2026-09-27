"""Regenerate selected PAPER003 V3.1 publication figures from frozen public outputs.

Run from repository root:
    python scripts/make_publication_figures_v31.py

The script reads only machine-readable files committed under results/frozen and
results/v31. Molecular depictions use the representative collision-pair CSV.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
RF = ROOT / "results" / "frozen"
RV = ROOT / "results" / "v31"


def save(fig, stem):
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def compactness_generalization():
    d = pd.read_csv(RF / "representation_benchmark_summary.csv")
    c = pd.read_csv(RF / "representation_complexity_performance.csv")
    order = ["Graph19", "RDKit2D", "ECFP4_2048"]
    labels = ["Graph19", "RDKit2D", "ECFP4"]
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 6.2), gridspec_kw={"width_ratios": [.85, 1.25]})
    fig.suptitle("Extreme representation compression trades predictive accuracy for compactness", fontsize=18, fontweight="bold")
    dims = [int(c.loc[c.representation == r, "n_features"].iloc[0]) for r in order]
    y = np.arange(3)
    axs[0].barh(y, dims); axs[0].set_xscale("log"); axs[0].set_yticks(y, labels); axs[0].invert_yaxis()
    axs[0].set_xlabel("Representation dimensionality (log scale)"); axs[0].set_title("A  Compactness", loc="left", fontweight="bold")
    for yi, v in zip(y, dims): axs[0].text(v * 1.08, yi, f"{v:,}", va="center")
    x = np.arange(3); width = .32
    for k, sp in enumerate(["random", "scaffold"]):
        sub = d[d.split_type == sp].set_index("representation").loc[order]
        pos = x + (k - .5) * width
        axs[1].bar(pos, sub["r2__mean"], width, yerr=sub["r2__std"], capsize=4, label=sp.capitalize())
    axs[1].set_xticks(x, labels); axs[1].set_ylabel("Mean R²"); axs[1].set_ylim(0, .84); axs[1].legend(frameon=False)
    axs[1].set_title("B  Random versus scaffold generalization", loc="left", fontweight="bold")
    save(fig, "Fig_compactness_generalization")


def graph19_interpretability():
    a = pd.read_csv(RF / "descriptor_family_ablation_summary.csv")
    imp = pd.read_csv(RF / "interpretability_descriptor_importance.csv").sort_values("permutation_importance_mean_r2_drop", ascending=False).head(10)
    fam = [("degree_only", "Degree"), ("distance_only", "Distance"), ("spectral_only", "Spectral"), ("degree_plus_distance", "Degree + distance"), ("full19", "Full Graph19")]
    fig, axs = plt.subplots(1, 2, figsize=(13.0, 7.2))
    y = np.arange(len(fam)); h = .34
    for k, sp in enumerate(["random", "scaffold"]):
        means=[]; std=[]
        for key, _ in fam:
            row=a[(a.family==key)&(a.split_type==sp)].iloc[0]; means.append(row["r2__mean"]); std.append(row["r2__std"])
        axs[0].barh(y+(k-.5)*h, means, h, xerr=std, capsize=3, label=sp.capitalize())
    axs[0].set_yticks(y, [x[1] for x in fam]); axs[0].invert_yaxis(); axs[0].set_xlabel("Mean R²"); axs[0].legend(frameon=False)
    axs[0].set_title("A  Descriptor-family ablation", loc="left", fontweight="bold")
    ii=np.arange(len(imp)); axs[1].barh(ii, imp["permutation_importance_mean_r2_drop"], xerr=imp["permutation_importance_std"], capsize=3)
    axs[1].set_yticks(ii, imp.descriptor); axs[1].invert_yaxis(); axs[1].set_xlabel("Held-out permutation ΔR²")
    axs[1].set_title("B  Individual descriptor importance", loc="left", fontweight="bold")
    save(fig, "Fig_graph19_interpretability")


def degeneracy_examples():
    s = json.loads((RV / "graph19_ambiguity_summary.json").read_text(encoding="utf-8"))
    pairs = pd.read_csv(RV / "representative_collision_pairs.csv").head(4)
    fig = plt.figure(figsize=(12.5, 8.8)); gs = fig.add_gridspec(3, 4, height_ratios=[.7, 1, 1])
    ax = fig.add_subplot(gs[0, :]); ax.axis("off")
    ax.text(.5,.72,f"{s['n_molecules']:,} molecules → {s['n_unique_graph19_vectors']:,} unique Graph19 vectors",ha="center",fontsize=18,fontweight="bold")
    ax.text(.5,.36,f"{s['n_molecules_in_collision_groups']:,} molecules ({100*s['fraction_molecules_in_collision_groups']:.1f}%) occur in {s['n_collision_groups']:,} exact collision groups; {s['pairs_tanimoto_lt_0_5_and_delta_pic50_gt_2']} discordant pairs satisfy Tanimoto < 0.5 and |ΔpIC50| > 2.",ha="center",fontsize=11)
    for r,(_,row) in enumerate(pairs.iterrows()):
        for c, suf in enumerate(["A","B"]):
            j=r*2+c; a=fig.add_subplot(gs[1+j//4, j%4]); a.axis("off")
            mol=Chem.MolFromSmiles(row[f"SMILES_{suf}"]); a.imshow(Draw.MolToImage(mol,size=(420,260)))
            a.set_title(f"Pair {r+1}{suf} | pIC50 {row[f'pIC50_{suf}']:.2f}\n{row[f'formula_{suf}']}",fontsize=8)
    save(fig, "Fig_graph19_representation_degeneracy")


def chemical_space_uncertainty():
    b=pd.read_csv(RF/"chemical_space_similarity_bins.csv")
    b=b[(b.split_type=="scaffold")&(b.n>=20)&(b.representation.isin(["Graph19","ECFP4_2048"]))]
    c=pd.read_csv(RF/"conformal_metrics.csv")
    fig,axs=plt.subplots(1,2,figsize=(12.4,6.6))
    for rep,label in [("Graph19","Graph19"),("ECFP4_2048","ECFP4")]:
        s=b[b.representation==rep].sort_values("mean_similarity"); axs[0].plot(s.mean_similarity,s.mean_abs_error,marker="o",label=label)
    axs[0].set_xlabel("Nearest-training ECFP4 Tanimoto similarity"); axs[0].set_ylabel("Mean absolute pIC50 error"); axs[0].legend(frameon=False)
    axs[0].set_title("A  Chemical-space distance",loc="left",fontweight="bold")
    reps=["Graph19","ECFP4_2048"]; labels=["Graph19","ECFP4"]
    random=[c[(c.representation==r)&(c.split_type=="random")].coverage.mean() for r in reps]
    scaffold=[c[(c.representation==r)&(c.split_type=="scaffold")].coverage.mean() for r in reps]
    x=np.arange(2); w=.34; axs[1].bar(x-w/2,random,w,label="Random"); axs[1].bar(x+w/2,scaffold,w,label="Scaffold"); axs[1].axhline(.9,ls="--",label="Nominal 90%")
    axs[1].set_xticks(x,labels); axs[1].set_ylim(.78,.94); axs[1].set_ylabel("Empirical coverage"); axs[1].legend(frameon=False)
    axs[1].set_title("B  Split-conformal coverage",loc="left",fontweight="bold")
    save(fig, "Fig_chemical_space_uncertainty")


if __name__ == "__main__":
    compactness_generalization(); graph19_interpretability(); degeneracy_examples(); chemical_space_uncertainty()
    print("Figures written to", OUT)
