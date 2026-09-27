from __future__ import annotations
from pathlib import Path
import hashlib, time
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from common import RESULTS, CACHE, CFG, jobs_from_cfg, metric_dict, json_dump, experiment_cache_path, cache_load, cache_save, sha256_file
from phase04_models import _rf, _run_cached_fit, _sig


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"05_generalization_ad"; p.mkdir(parents=True,exist_ok=True); return p

def _max_tanimoto_for_split(smiles, train_idx, test_idx, radius, nbits, cache_path):
    if cache_path.exists():
        d=pd.read_csv(cache_path)
        if len(d)==len(test_idx): return d
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=radius,fpSize=nbits,includeChirality=False)
    fps=[gen.GetFingerprint(Chem.MolFromSmiles(str(s))) for s in smiles]
    trainfps=[fps[int(i)] for i in train_idx]
    rows=[]
    for i in test_idx:
        sims=DataStructs.BulkTanimotoSimilarity(fps[int(i)],trainfps)
        rows.append({"row_index":int(i),"max_train_tanimoto":float(max(sims)) if sims else np.nan})
    d=pd.DataFrame(rows); cache_path.parent.mkdir(parents=True,exist_ok=True); d.to_csv(cache_path,index=False)
    return d

def _safe_metrics(df):
    if len(df)<2: return {"r2":np.nan,"rmse":np.nan,"mae":np.nan,"spearman":np.nan,"pearson":np.nan,"n":len(df)}
    return metric_dict(df["y_true"],df["y_pred"])


def run(logger, smoke=False):
    outdir=_out(smoke); base=RESULTS/("smoke" if smoke else "full")
    g=pd.read_csv(base/"02_descriptors"/"graph19_regenerated.csv")
    y=g["pIC50"].to_numpy(float); X=g.iloc[:,3:].to_numpy(float); idx=np.arange(len(g)); primary=int(CFG["primary_seed"])
    scaffold=pd.read_csv(base/"03_representations_splits"/"scaffold_assignments.csv")
    hist_params=CFG["historical_rf"]
    data_sig=_sig(sha256_file(base/"02_descriptors"/"graph19_regenerated.csv"),hist_params,"scaffold-assurance")

    folds=sorted(scaffold["fold"].unique());
    if smoke: folds=folds[:2]
    srows=[]
    for fold in folds:
        te=scaffold.index[scaffold["fold"]==fold].to_numpy(); tr=scaffold.index[scaffold["fold"]!=fold].to_numpy()
        model=Pipeline([("scaler",StandardScaler()),("model",_rf(hist_params,primary,smoke))])
        res,_=_run_cached_fit("historical_scaffold",f"rf600_scaffold_fold{int(fold)}",_sig(data_sig,fold),model,X[tr],y[tr],X[te],y[te],te,smoke)
        srows.append({"fold":int(fold),"train_n":len(tr),"test_n":len(te),**res})
    sdf=pd.DataFrame(srows); scaffold_metrics=outdir/"historical_scaffold_rf600_metrics.csv"; sdf.to_csv(scaffold_metrics,index=False)
    ssum={c:{"mean":float(sdf[c].mean()),"std":float(sdf[c].std(ddof=0))} for c in ["r2","rmse","mae"]}
    if not smoke:
        expected={"r2":0.210998,"rmse":1.196759,"mae":0.928434}
        ssum["historical_expected"]=expected
        ssum["within_0.01_of_handoff"]={k:abs(ssum[k]["mean"]-v)<0.01 for k,v in expected.items()}
    scaffold_summary=outdir/"historical_scaffold_rf600_summary.json"; json_dump(ssum,scaffold_summary)

    tr,te=train_test_split(idx,test_size=float(CFG["test_size"]),random_state=primary,shuffle=True)
    scaler=StandardScaler(); Xtr_s=scaler.fit_transform(X[tr]); Xte_s=scaler.transform(X[te])
    A=np.column_stack([np.ones(len(tr)),Xtr_s]); B=np.column_stack([np.ones(len(te)),Xte_s])
    inv=np.linalg.pinv(A.T@A)
    htrain=np.einsum("ij,jk,ik->i",A,inv,A); htest=np.einsum("ij,jk,ik->i",B,inv,B)
    hstar=3*(X.shape[1]+1)/len(tr)
    model=Pipeline([("scaler",StandardScaler()),("model",_rf(hist_params,primary,smoke))]); model.fit(X[tr],y[tr]); pred=model.predict(X[te])
    resid=y[te]-pred; rsd=float(np.std(resid,ddof=1)); std_resid=resid/rsd if rsd>0 else resid*np.nan
    w=pd.DataFrame({"row_index":te,"InChIKey":g.iloc[te]["InChIKey"].to_numpy(),"leverage":htest,"h_star":hstar,
                    "y_true":y[te],"y_pred":pred,"residual":resid,"standardized_residual":std_resid,
                    "high_leverage":htest>hstar,"residual_outlier":np.abs(std_resid)>3})
    wpath=outdir/"williams_test_data.csv"; w.to_csv(wpath,index=False)
    wsum={"p":int(X.shape[1]),"n_train":int(len(tr)),"n_test":int(len(te)),"h_star":float(hstar),
          "expected_h_star_for_full_run":float(3*20/8044),"high_leverage_test_n":int(w["high_leverage"].sum()),
          "residual_outlier_test_n":int(w["residual_outlier"].sum()),
          "residual_standardization":"test residual divided by sample standard deviation of test residuals; shown for diagnostic Williams plotting",
          "leverage_definition":"x^T (X_train^T X_train)^+ x with intercept; threshold 3(p+1)/n_train"}
    wsum_path=outdir/"williams_summary.json"; json_dump(wsum,wsum_path)

    preds=pd.read_csv(base/"04_models"/"representation_benchmark_predictions.csv")
    radius=int(CFG["morgan"]["radius"]); nbits=int(CFG["morgan"]["n_bits"])
    sim_cache=CACHE/("smoke" if smoke else "full")/"similarity"; sim_cache.mkdir(parents=True,exist_ok=True)
    sim_parts=[]
    trr,ter=train_test_split(idx,test_size=float(CFG["test_size"]),random_state=primary,shuffle=True)
    sim=_max_tanimoto_for_split(g["SMILES"].tolist(),trr,ter,radius,nbits,sim_cache/f"random_seed{primary}.csv")
    sim["split_type"]="random"; sim["split_id"]=primary; sim_parts.append(sim)
    use_folds=folds
    for fold in use_folds:
        tef=scaffold.index[scaffold["fold"]==fold].to_numpy(); trf=scaffold.index[scaffold["fold"]!=fold].to_numpy()
        sim=_max_tanimoto_for_split(g["SMILES"].tolist(),trf,tef,radius,nbits,sim_cache/f"scaffold_fold{int(fold)}.csv")
        sim["split_type"]="scaffold"; sim["split_id"]=int(fold); sim_parts.append(sim)
    sims=pd.concat(sim_parts,ignore_index=True)
    keep_preds=preds[preds["representation"].isin(["Graph19","ECFP4_2048"])].copy()
    chem=keep_preds.merge(sims,on=["row_index","split_type","split_id"],how="inner")
    chem["abs_error"]=np.abs(chem["y_true"]-chem["y_pred"])
    chem_path=outdir/"chemical_space_error_vs_similarity.csv"; chem.to_csv(chem_path,index=False)

    bins=CFG["similarity_bins"]
    chem["similarity_bin"]=pd.cut(chem["max_train_tanimoto"],bins=bins,right=False,include_lowest=True)
    bin_rows=[]
    for keys,d in chem.groupby(["representation","split_type","similarity_bin"],observed=True):
        m=_safe_metrics(d); bin_rows.append({"representation":keys[0],"split_type":keys[1],"similarity_bin":str(keys[2]),
                                            "mean_similarity":float(d["max_train_tanimoto"].mean()),"mean_abs_error":float(d["abs_error"].mean()),**m})
    bins_df=pd.DataFrame(bin_rows); bins_path=outdir/"chemical_space_similarity_bins.csv"; bins_df.to_csv(bins_path,index=False)
    ad_rows=[]
    for (rep,stype),d0 in chem.groupby(["representation","split_type"]):
        for th in CFG["tanimoto_thresholds"]:
            for region,mask in [("inside",d0["max_train_tanimoto"]>=th),("outside",d0["max_train_tanimoto"]<th)]:
                d=d0[mask]; m=_safe_metrics(d)
                ad_rows.append({"representation":rep,"split_type":stype,"threshold":th,"region":region,
                                "coverage_fraction":float(len(d)/len(d0)) if len(d0) else np.nan,**m})
    ad=pd.DataFrame(ad_rows); ad_path=outdir/"tanimoto_ad_metrics.csv"; ad.to_csv(ad_path,index=False)
    chem_summary={"fingerprint":"Morgan/ECFP4 binary radius 2, 2048 bits","purpose":"chemical-space distance only; complements rather than replaces Williams leverage AD",
                  "random_primary_seed":primary,"scaffold_folds_evaluated":[int(x) for x in use_folds]}
    chem_summary_path=outdir/"chemical_space_ad_summary.json"; json_dump(chem_summary,chem_summary_path)
    logger.info("Generalization/AD complete: h*=%.8f; scaffold R2 mean=%.4f",hstar,ssum["r2"]["mean"])
    outputs=[scaffold_metrics,scaffold_summary,wpath,wsum_path,chem_path,bins_path,ad_path,chem_summary_path]
    return outputs,{"scaffold":ssum,"williams":wsum,"chemical_space":chem_summary}
