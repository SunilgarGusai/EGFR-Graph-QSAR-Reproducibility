from __future__ import annotations
from pathlib import Path
import json, shutil, zipfile, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import ROOT, RESULTS, LOGS, STATE, CFG, sha256_file, json_dump, environment_report, pip_freeze, input_manifest_check


def _base(smoke): return RESULTS/("smoke" if smoke else "full")
def _out(smoke):
    p=_base(smoke)/"07_final"; p.mkdir(parents=True,exist_ok=True); return p

def _savefig(fig,path):
    fig.tight_layout(); fig.savefig(path.with_suffix(".png"),dpi=300,bbox_inches="tight"); fig.savefig(path.with_suffix(".pdf"),bbox_inches="tight"); plt.close(fig)

def run(logger,smoke=False):
    base=_base(smoke); outdir=_out(smoke); figs=outdir/"figures"; tables=outdir/"tables"; figs.mkdir(exist_ok=True); tables.mkdir(exist_ok=True)
    outputs=[]
    bench=pd.read_csv(base/"04_models"/"representation_benchmark_metrics.csv")
    bsum=bench.groupby(["representation","split_type"])["r2"].agg(["mean","std"]).reset_index()
    pivot=bsum.pivot(index="representation",columns="split_type",values="mean")
    ax=pivot.plot(kind="bar",yerr=bsum.pivot(index="representation",columns="split_type",values="std"),capsize=3,figsize=(8,5))
    ax.set_ylabel("R²"); ax.set_xlabel("Molecular representation"); ax.set_title("Representation performance under random and scaffold splits"); ax.axhline(0,linewidth=0.8)
    _savefig(ax.figure,figs/"fig_representation_random_vs_scaffold"); outputs += [figs/"fig_representation_random_vs_scaffold.png",figs/"fig_representation_random_vs_scaffold.pdf"]

    bins=pd.read_csv(base/"05_generalization_ad"/"chemical_space_similarity_bins.csv")
    fig,ax=plt.subplots(figsize=(8,5))
    for (rep,stype),d in bins.groupby(["representation","split_type"]): ax.plot(d["mean_similarity"],d["mean_abs_error"],marker="o",label=f"{rep} / {stype}")
    ax.set_xlabel("Mean nearest-training ECFP4 Tanimoto similarity"); ax.set_ylabel("Mean absolute prediction error"); ax.set_title("Prediction error across chemical-space distance"); ax.legend(fontsize=8)
    _savefig(fig,figs/"fig_error_vs_chemical_similarity"); outputs += [figs/"fig_error_vs_chemical_similarity.png",figs/"fig_error_vs_chemical_similarity.pdf"]

    abl=pd.read_csv(base/"06_strengthening"/"descriptor_family_ablation.csv")
    asu=abl.groupby(["family","split_type"])["r2"].agg(["mean","std"]).reset_index(); piv=asu.pivot(index="family",columns="split_type",values="mean")
    ax=piv.plot(kind="bar",yerr=asu.pivot(index="family",columns="split_type",values="std"),capsize=3,figsize=(9,5)); ax.set_ylabel("R²"); ax.set_xlabel("Descriptor family"); ax.set_title("Ablation of classical graph-descriptor families"); ax.axhline(0,linewidth=0.8)
    _savefig(ax.figure,figs/"fig_descriptor_family_ablation"); outputs += [figs/"fig_descriptor_family_ablation.png",figs/"fig_descriptor_family_ablation.pdf"]

    imp=pd.read_csv(base/"06_strengthening"/"interpretability_descriptor_importance.csv").head(12).sort_values("permutation_importance_mean_r2_drop")
    fig,ax=plt.subplots(figsize=(8,6)); ax.barh(imp["descriptor"],imp["permutation_importance_mean_r2_drop"],xerr=imp["permutation_importance_std"]); ax.set_xlabel("Decrease in held-out R² after permutation"); ax.set_title("Held-out permutation importance of graph descriptors")
    _savefig(fig,figs/"fig_permutation_importance"); outputs += [figs/"fig_permutation_importance.png",figs/"fig_permutation_importance.pdf"]

    conf=pd.read_csv(base/"06_strengthening"/"conformal_metrics.csv"); cs=conf.groupby(["representation","split_type"])["coverage"].agg(["mean","std"]).reset_index(); cp=cs.pivot(index="representation",columns="split_type",values="mean")
    ax=cp.plot(kind="bar",yerr=cs.pivot(index="representation",columns="split_type",values="std"),capsize=3,figsize=(7,5)); ax.axhline(1-float(CFG["conformal_alpha"]),linestyle="--",label="Nominal coverage"); ax.set_ylim(0,1.05); ax.set_ylabel("Empirical coverage"); ax.set_xlabel("Representation"); ax.set_title("Split-conformal coverage under random and scaffold evaluation"); ax.legend()
    _savefig(ax.figure,figs/"fig_conformal_coverage"); outputs += [figs/"fig_conformal_coverage.png",figs/"fig_conformal_coverage.pdf"]

    w=pd.read_csv(base/"05_generalization_ad"/"williams_test_data.csv"); hstar=float(w["h_star"].iloc[0]); fig,ax=plt.subplots(figsize=(8,5)); ax.scatter(w["leverage"],w["standardized_residual"],s=10,alpha=0.6); ax.axvline(hstar,linestyle="--",label=f"h*={hstar:.6f}"); ax.axhline(3,linestyle="--"); ax.axhline(-3,linestyle="--"); ax.set_xlabel("Leverage"); ax.set_ylabel("Standardized residual"); ax.set_title("Williams applicability-domain diagnostic"); ax.legend()
    _savefig(fig,figs/"fig_williams_corrected"); outputs += [figs/"fig_williams_corrected.png",figs/"fig_williams_corrected.pdf"]

    key_tables=[base/"01_audit"/"provenance_summary.json",base/"02_descriptors"/"descriptor_comparison_summary.json",base/"03_representations_splits"/"representation_split_summary.json",base/"04_models"/"assurance_holdout_metrics.csv",base/"04_models"/"historical_10fold_training_cv_summary.json",base/"04_models"/"y_randomization_summary.json",base/"04_models"/"representation_benchmark_summary.csv",base/"05_generalization_ad"/"historical_scaffold_rf600_summary.json",base/"05_generalization_ad"/"williams_summary.json",base/"05_generalization_ad"/"chemical_space_similarity_bins.csv",base/"05_generalization_ad"/"tanimoto_ad_metrics.csv",base/"06_strengthening"/"descriptor_family_ablation_summary.csv",base/"06_strengthening"/"interpretability_descriptor_importance.csv",base/"06_strengthening"/"interpretability_grouped_family_importance.csv",base/"06_strengthening"/"conformal_summary.csv",base/"06_strengthening"/"representation_complexity_performance.csv"]
    for p in key_tables: dest=tables/p.name; shutil.copy2(p,dest); outputs.append(dest)

    ass=pd.read_csv(base/"04_models"/"assurance_holdout_metrics.csv"); rf=ass[ass["model"]=="RF600_postV2_canonical"].iloc[0]; ss=json.loads((base/"05_generalization_ad"/"historical_scaffold_rf600_summary.json").read_text()); ws=json.loads((base/"05_generalization_ad"/"williams_summary.json").read_text()); prov=json.loads((base/"01_audit"/"provenance_summary.json").read_text())
    comp=[
      {"item":"Initial ChEMBL activity records","V2":"17705","V3":prov["raw_rows"],"difference":"none expected","reason":"archived-source assurance"},
      {"item":"Stage-1 records","V2":"17412","V3":prov["stage1_rows"],"difference":"none expected","reason":"historical curation reproduced"},
      {"item":"Unique compounds","V2":"10056","V3":prov["final_unique_compounds"],"difference":"none expected","reason":"InChIKey median aggregation reproduced"},
      {"item":"Classical descriptor count","V2":"19","V3":"19","difference":"none","reason":"scientific identity preserved"},
      {"item":"RF holdout R2","V2":"0.397 (original manuscript)","V3":f"{rf['r2']:.6f}","difference":"post-V2 600-tree canonical assurance configuration","reason":"original V2 RF hyperparameters were not recoverable; CHHELI_GANATRI 600-tree config is provenance-backed"},
      {"item":"Scaffold R2 mean","V2":"post-V2 exploratory 0.210998","V3":f"{ss['r2']['mean']:.6f}","difference":"definitive reproducible rerun","reason":"same Murcko GroupKFold concept promoted to core result"},
      {"item":"Williams h*","V2":"old figure likely ~0.005967 (full cohort denominator)","V3":f"{ws['h_star']:.9f}","difference":"corrected","reason":"threshold uses training n and p=19: 3(p+1)/n_train"},
      {"item":"Modern representation comparison","V2":"absent","V3":"Graph19 vs RDKit2D vs ECFP4","difference":"added","reason":"contemporary benchmark"},
      {"item":"Chemical-space AD","V2":"Williams leverage only","V3":"Williams + nearest-training ECFP4/Tanimoto","difference":"added","reason":"descriptor-space leverage and chemical-space distance"},
      {"item":"Descriptor-family ablation","V2":"absent","V3":"degree, distance, spectral, degree+distance, full19","difference":"added","reason":"tests mathematical representation families"},
      {"item":"Uncertainty","V2":"absent","V3":"90% split conformal + scaffold-shift empirical coverage","difference":"added","reason":"reliability layer"},
      {"item":"Interpretability","V2":"RF impurity importance","V3":"impurity + held-out permutation + grouped-family permutation","difference":"strengthened","reason":"correlation-aware robustness"}
    ]
    cpath=tables/"V2_V3_COMPARISON.csv"; pd.DataFrame(comp).to_csv(cpath,index=False); outputs.append(cpath)

    env=environment_report(); envpath=outdir/"environment.json"; json_dump(env,envpath); outputs.append(envpath); freeze=outdir/"environment_freeze.txt"; pip_freeze(freeze); outputs.append(freeze)
    manifest=[]
    for p in sorted(base.rglob("*")):
        if p.is_file() and p.name!="EGFR_QSAR_REPRODUCIBILITY_RETURN.zip": manifest.append({"path":str(p.relative_to(ROOT)),"sha256":sha256_file(p),"bytes":p.stat().st_size})
    manpath=outdir/"results_manifest.json"; json_dump({"generated_epoch":time.time(),"environment":env,"files":manifest},manpath); outputs.append(manpath)
    qa={"input_manifest_all_ok":all(v["ok"] for v in input_manifest_check().values()),"descriptor_assurance":json.loads((base/"02_descriptors"/"descriptor_comparison_summary.json").read_text()),"scaffold_assurance":ss,"williams":ws,"manuscript_rewrite_status":"NOT STARTED by design; freeze results first, then rewrite V3.","github_status":"Repository-ready outputs produced locally; public repository creation/release is a later manuscript-finalization step."}
    qapath=outdir/"QA_SUMMARY.json"; json_dump(qa,qapath); outputs.append(qapath)

    ret=outdir/"EGFR_QSAR_REPRODUCIBILITY_RETURN.zip"
    with zipfile.ZipFile(ret,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for folder in [base,LOGS,STATE]:
            if folder.exists():
                for p in folder.rglob("*"):
                    if p.is_file() and p!=ret: z.write(p,p.relative_to(ROOT))
        for p in [ROOT/"config.json",ROOT/"INPUT_MANIFEST.json",ROOT/"docs"/"V3_SCIENTIFIC_DESIGN_FREEZE.md",ROOT/"docs"/"NOVELTY_AND_JOURNAL_AUDIT.md",ROOT/"docs"/"V3_GAP_MATRIX.csv"]:
            if p.exists(): z.write(p,p.relative_to(ROOT))
    outputs.append(ret)
    done=outdir/"RUN_COMPLETE.txt"; done.write_text("EGFR Graph QSAR V3 computational assurance/strengthening run completed. Upload EGFR_QSAR_REPRODUCIBILITY_RETURN.zip back to ChatGPT for result audit and manuscript V3 writing.\n",encoding="utf-8"); outputs.append(done)
    logger.info("Final QA complete. Return package: %s",ret)
    return outputs,{"return_package":str(ret),"qa":qa}
