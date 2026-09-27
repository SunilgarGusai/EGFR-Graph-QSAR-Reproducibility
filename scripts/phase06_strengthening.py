from __future__ import annotations
from pathlib import Path
import hashlib, time, math
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from common import RESULTS, CACHE, CFG, jobs_from_cfg, metric_dict, json_dump, experiment_cache_path, cache_load, cache_save, sha256_file
from phase04_models import _rf, _run_cached_fit, _sig


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"06_strengthening"; p.mkdir(parents=True,exist_ok=True); return p

def _quantile_conformal(residuals,alpha):
    r=np.sort(np.asarray(residuals,float)); n=len(r)
    k=min(n,max(1,int(math.ceil((n+1)*(1-alpha)))))
    return float(r[k-1]),k,n

def _conformal_cached(key,signature,X,y,proper,cal,test,seed,smoke=False):
    cp=experiment_cache_path("conformal",key,smoke); pp=cp.with_suffix(".pred.csv")
    cached=cache_load(cp,signature)
    if cached is not None and pp.exists(): return cached["result"],pd.read_csv(pp)
    model=_rf(CFG["benchmark_rf"],seed,smoke)
    t0=time.perf_counter(); model.fit(X[proper],y[proper]); train_s=time.perf_counter()-t0
    cal_pred=model.predict(X[cal]); absr=np.abs(y[cal]-cal_pred)
    q,k,ncal=_quantile_conformal(absr,float(CFG["conformal_alpha"]))
    t0=time.perf_counter(); pred=model.predict(X[test]); pred_s=time.perf_counter()-t0
    lo=pred-q; hi=pred+q; coverage=float(np.mean((y[test]>=lo)&(y[test]<=hi)))
    result={**metric_dict(y[test],pred),"coverage":coverage,"mean_interval_width":float(2*q),"q":q,
            "alpha":float(CFG["conformal_alpha"]),"calibration_n":int(ncal),"finite_sample_rank":int(k),
            "train_seconds":train_s,"predict_seconds":pred_s}
    pdf=pd.DataFrame({"row_index":test,"y_true":y[test],"y_pred":pred,"lower":lo,"upper":hi,"covered":(y[test]>=lo)&(y[test]<=hi)})
    pdf.to_csv(pp,index=False); cache_save(cp,signature,{"result":result,"prediction_file":str(pp)})
    return result,pdf


def run(logger,smoke=False):
    outdir=_out(smoke); base=RESULTS/("smoke" if smoke else "full")
    g=pd.read_csv(base/"02_descriptors"/"graph19_regenerated.csv"); y=g["pIC50"].to_numpy(float); idx=np.arange(len(g))
    feature_names=list(g.columns[3:]); Xg=g[feature_names].to_numpy(float)
    Xm=np.load(base/"03_representations_splits"/"morgan_ecfp4_2048.npz")["X"].astype(np.float32,copy=False)
    scaffold=pd.read_csv(base/"03_representations_splits"/"scaffold_assignments.csv")
    primary=int(CFG["primary_seed"]); test_size=float(CFG["test_size"]); data_sig=sha256_file(base/"02_descriptors"/"graph19_regenerated.csv")

    groups={
      "degree_only":["M1","M2","Randic","Harmonic","ABC","GA","AZI","F","HM","SCI"],
      "distance_only":["Wiener","Harary","Szeged","Schultz","Gutman","EdgeWiener","VertexEdgeWiener"],
      "spectral_only":["Energy","Estrada"],
      "degree_plus_distance":["M1","M2","Randic","Harmonic","ABC","GA","AZI","F","HM","SCI","Wiener","Harary","Szeged","Schultz","Gutman","EdgeWiener","VertexEdgeWiener"],
      "full19":feature_names
    }
    seeds=CFG["random_seeds"] if not smoke else CFG["random_seeds"][:2]
    folds=sorted(scaffold["fold"].unique())
    if smoke: folds=folds[:2]
    arows=[]
    for gname,cols in groups.items():
        X=g[cols].to_numpy(float)
        for seed in seeds:
            tr,te=train_test_split(idx,test_size=test_size,random_state=int(seed),shuffle=True)
            model=_rf(CFG["benchmark_rf"],int(seed),smoke)
            res,_=_run_cached_fit("ablation",f"{gname}_random_{seed}",_sig(data_sig,"ablation",gname,"random",seed),model,X[tr],y[tr],X[te],y[te],te,smoke)
            arows.append({"family":gname,"n_features":len(cols),"split_type":"random","split_id":int(seed),**res})
        for fold in folds:
            te=scaffold.index[scaffold["fold"]==fold].to_numpy(); tr=scaffold.index[scaffold["fold"]!=fold].to_numpy()
            model=_rf(CFG["benchmark_rf"],primary+int(fold),smoke)
            res,_=_run_cached_fit("ablation",f"{gname}_scaffold_{int(fold)}",_sig(data_sig,"ablation",gname,"scaffold",fold),model,X[tr],y[tr],X[te],y[te],te,smoke)
            arows.append({"family":gname,"n_features":len(cols),"split_type":"scaffold","split_id":int(fold),**res})
    adf=pd.DataFrame(arows); ablation_path=outdir/"descriptor_family_ablation.csv"; adf.to_csv(ablation_path,index=False)
    asum=(adf.groupby(["family","n_features","split_type"])[["r2","rmse","mae","spearman"]].agg(["mean","std"]).reset_index())
    asum.columns=["__".join([str(x) for x in c if x]) if isinstance(c,tuple) else c for c in asum.columns]
    ablation_summary=outdir/"descriptor_family_ablation_summary.csv"; asum.to_csv(ablation_summary,index=False)

    tr,te=train_test_split(idx,test_size=test_size,random_state=primary,shuffle=True)
    hist=CFG["historical_rf"]
    pipe=Pipeline([("scaler",StandardScaler()),("model",_rf(hist,primary,smoke))]); pipe.fit(Xg[tr],y[tr]); pred=pipe.predict(Xg[te])
    impurity=pipe.named_steps["model"].feature_importances_
    nrep=min(3,int(CFG["permutation_repeats"])) if smoke else int(CFG["permutation_repeats"])
    perm=permutation_importance(pipe,Xg[te],y[te],n_repeats=nrep,random_state=primary,scoring="r2",n_jobs=jobs_from_cfg())
    idf=pd.DataFrame({"descriptor":feature_names,"impurity_importance":impurity,
                      "permutation_importance_mean_r2_drop":perm.importances_mean,"permutation_importance_std":perm.importances_std})
    idf=idf.sort_values("permutation_importance_mean_r2_drop",ascending=False)
    importance_path=outdir/"interpretability_descriptor_importance.csv"; idf.to_csv(importance_path,index=False)
    base_r2=float(r2_score(y[te],pred)); rng=np.random.default_rng(primary); grows=[]
    group_cols={"degree":groups["degree_only"],"distance":groups["distance_only"],"spectral":groups["spectral_only"]}
    col_index={c:i for i,c in enumerate(feature_names)}
    for fam,cols in group_cols.items():
        drops=[]; jj=[col_index[c] for c in cols]
        for rep in range(nrep):
            permidx=rng.permutation(len(te)); Xp=Xg[te].copy(); Xp[:,jj]=Xp[permidx][:,jj]
            rp=float(r2_score(y[te],pipe.predict(Xp))); drops.append(base_r2-rp)
        grows.append({"family":fam,"n_features":len(cols),"baseline_r2":base_r2,"mean_group_permutation_r2_drop":float(np.mean(drops)),"std":float(np.std(drops,ddof=0))})
    gimp=pd.DataFrame(grows); group_importance_path=outdir/"interpretability_grouped_family_importance.csv"; gimp.to_csv(group_importance_path,index=False)

    reps={"Graph19":Xg,"ECFP4_2048":Xm}; crows=[]; cpreds=[]
    outer_tr,outer_te=train_test_split(idx,test_size=test_size,random_state=primary,shuffle=True)
    proper,cal=train_test_split(outer_tr,test_size=float(CFG["conformal_calibration_fraction_within_training"]),random_state=primary+1,shuffle=True)
    for repname,X in reps.items():
        sig=_sig(data_sig,"conformal",repname,"random",primary,X.shape,CFG["conformal_alpha"])
        res,pdf=_conformal_cached(f"{repname}_random",sig,X,y,proper,cal,outer_te,primary,smoke)
        crows.append({"representation":repname,"split_type":"random","split_id":primary,"coverage_guarantee_context":"exchangeability approximately assumed under random split",**res})
        pdf["representation"]=repname; pdf["split_type"]="random"; pdf["split_id"]=primary; cpreds.append(pdf)
        for fold in folds:
            test=scaffold.index[scaffold["fold"]==fold].to_numpy(); train=scaffold.index[scaffold["fold"]!=fold].to_numpy()
            proper2,cal2=train_test_split(train,test_size=float(CFG["conformal_calibration_fraction_within_training"]),random_state=primary+100+int(fold),shuffle=True)
            sig=_sig(data_sig,"conformal",repname,"scaffold",fold,X.shape,CFG["conformal_alpha"])
            res,pdf=_conformal_cached(f"{repname}_scaffold_{int(fold)}",sig,X,y,proper2,cal2,test,primary+int(fold),smoke)
            crows.append({"representation":repname,"split_type":"scaffold","split_id":int(fold),
                          "coverage_guarantee_context":"diagnostic only under scaffold shift; ordinary split-conformal exchangeability guarantee does not hold",**res})
            pdf["representation"]=repname; pdf["split_type"]="scaffold"; pdf["split_id"]=int(fold); cpreds.append(pdf)
    cdf=pd.DataFrame(crows); conformal_path=outdir/"conformal_metrics.csv"; cdf.to_csv(conformal_path,index=False)
    cpd=pd.concat(cpreds,ignore_index=True); conformal_pred_path=outdir/"conformal_predictions.csv"; cpd.to_csv(conformal_pred_path,index=False)
    csum=(cdf.groupby(["representation","split_type"])[["coverage","mean_interval_width","r2","rmse","mae"]].agg(["mean","std"]).reset_index())
    csum.columns=["__".join([str(x) for x in c if x]) if isinstance(c,tuple) else c for c in csum.columns]
    conformal_summary=outdir/"conformal_summary.csv"; csum.to_csv(conformal_summary,index=False)

    bench=pd.read_csv(base/"04_models"/"representation_benchmark_metrics.csv")
    dims={"Graph19":Xg.shape[1],"RDKit2D":pd.read_csv(base/"03_representations_splits"/"rdkit2d.csv",nrows=1).shape[1]-3,"ECFP4_2048":Xm.shape[1]}
    comp=[]
    for rep,d in dims.items():
        x=bench[bench["representation"]==rep]
        comp.append({"representation":rep,"n_features":int(d),"mean_train_seconds":float(x["train_seconds"].mean()),
                     "mean_predict_seconds":float(x["predict_seconds"].mean()),"mean_random_r2":float(x[x["split_type"]=="random"]["r2"].mean()),
                     "mean_scaffold_r2":float(x[x["split_type"]=="scaffold"]["r2"].mean())})
    complexity=pd.DataFrame(comp); complexity_path=outdir/"representation_complexity_performance.csv"; complexity.to_csv(complexity_path,index=False)
    meta={"interpretability_note":"Impurity importance retained only as historical/contextual; held-out permutation and grouped-family permutation are the primary complementary analyses because graph descriptors are correlated.",
          "conformal_note":"90% split-conformal intervals are calibrated under random exchangeability; scaffold-shift coverage is reported as empirical stress-test coverage, not as a formal guarantee.",
          "shap_decision":"Not included: no added scientific question beyond permutation/grouped importance, and it would increase decorative complexity.",
          "gnn_decision":"Not included in core V3: recent 2026 EGFR benchmarking already compared a basic GCN with RF+ECFP4 at essentially the same dataset scale; reproducing a GPU/deep-learning track would dilute the classical-graph question."}
    meta_path=outdir/"strengthening_metadata.json"; json_dump(meta,meta_path)
    logger.info("Strengthening phase complete: ablation=%d fits; conformal=%d evaluations",len(adf),len(cdf))
    outputs=[ablation_path,ablation_summary,importance_path,group_importance_path,conformal_path,conformal_pred_path,conformal_summary,complexity_path,meta_path]
    return outputs,{"ablation_fits":len(adf),"conformal_evaluations":len(cdf),"metadata":meta}
