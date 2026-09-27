from __future__ import annotations
from pathlib import Path
import json, time, hashlib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from common import RESULTS, CACHE, CFG, jobs_from_cfg, metric_dict, sha256_file, experiment_cache_path, cache_load, cache_save, json_dump, package_signature


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"04_models"; p.mkdir(parents=True,exist_ok=True); return p

def _sig(*parts):
    h=hashlib.sha256()
    h.update(package_signature().encode())
    for p in parts: h.update(str(p).encode())
    return h.hexdigest()

def _rf(params, seed, smoke=False):
    trees=int(CFG["smoke_test_trees"]) if smoke else int(params["n_estimators"])
    return RandomForestRegressor(n_estimators=trees,max_features=params.get("max_features",1.0),
        min_samples_leaf=int(params.get("min_samples_leaf",1)),random_state=int(seed),n_jobs=jobs_from_cfg())

def _run_cached_fit(category,key,signature,model,Xtr,ytr,Xte,yte,test_idx,smoke=False):
    cp=experiment_cache_path(category,key,smoke); pp=cp.with_suffix(".pred.csv")
    cached=cache_load(cp,signature)
    if cached is not None and pp.exists(): return cached["result"],pd.read_csv(pp)
    t0=time.perf_counter(); model.fit(Xtr,ytr); train_s=time.perf_counter()-t0
    t0=time.perf_counter(); pred=model.predict(Xte); pred_s=time.perf_counter()-t0
    result={**metric_dict(yte,pred),"train_seconds":train_s,"predict_seconds":pred_s}
    pd.DataFrame({"row_index":np.asarray(test_idx,int),"y_true":np.asarray(yte,float),"y_pred":np.asarray(pred,float)}).to_csv(pp,index=False)
    cache_save(cp,signature,{"result":result,"prediction_file":str(pp)})
    return result,pd.read_csv(pp)

def _load_reps(smoke=False):
    base=RESULTS/("smoke" if smoke else "full")
    g=pd.read_csv(base/"02_descriptors"/"graph19_regenerated.csv")
    rd=pd.read_csv(base/"03_representations_splits"/"rdkit2d.csv")
    mo=np.load(base/"03_representations_splits"/"morgan_ecfp4_2048.npz")["X"]
    return g,rd,mo


def run(logger, smoke=False):
    outdir=_out(smoke); g,rd,mo=_load_reps(smoke)
    y=g["pIC50"].to_numpy(float); n=len(y); idx=np.arange(n)
    Xg=g.iloc[:,3:].to_numpy(float); Xrd=rd.iloc[:,3:].to_numpy(float); Xm=mo.astype(np.float32,copy=False)
    primary=int(CFG["primary_seed"]); test_size=float(CFG["test_size"])
    tr,te=train_test_split(idx,test_size=test_size,random_state=primary,shuffle=True)
    hist_params=CFG["historical_rf"]
    base_sig=_sig(sha256_file(RESULTS/("smoke" if smoke else "full")/"02_descriptors"/"graph19_regenerated.csv"),hist_params,primary,smoke)

    assurance=[]
    models={
      "MLR":Pipeline([("scaler",StandardScaler()),("model",LinearRegression())]),
      "SVR_RBF_default":Pipeline([("scaler",StandardScaler()),("model",SVR())]),
      "RF600_postV2_canonical":Pipeline([("scaler",StandardScaler()),("model",_rf(hist_params,primary,smoke))])
    }
    for name,model in models.items():
        key=f"assurance_{name}_seed{primary}".replace("/","_")
        res,pred=_run_cached_fit("assurance",key,_sig(base_sig,name),model,Xg[tr],y[tr],Xg[te],y[te],te,smoke)
        assurance.append({"model":name,"seed":primary,"train_n":len(tr),"test_n":len(te),**res})
    assurance_path=outdir/"assurance_holdout_metrics.csv"; pd.DataFrame(assurance).to_csv(assurance_path,index=False)

    cv_splits=int(CFG["historical_cv_folds"]); cv=KFold(n_splits=cv_splits,shuffle=True,random_state=primary)
    cv_rows=[]
    for fold,(a,b) in enumerate(cv.split(tr)):
        tri=tr[a]; vai=tr[b]; key=f"rf600_trainCV_fold{fold}"
        model=Pipeline([("scaler",StandardScaler()),("model",_rf(hist_params,primary,smoke))])
        res,_=_run_cached_fit("historical_cv",key,_sig(base_sig,"cv",fold),model,Xg[tri],y[tri],Xg[vai],y[vai],vai,smoke)
        cv_rows.append({"fold":fold,**res})
        if smoke and fold>=2: break
    cv_df=pd.DataFrame(cv_rows); cv_path=outdir/"historical_10fold_training_cv.csv"; cv_df.to_csv(cv_path,index=False)
    cv_summary={c:{"mean":float(cv_df[c].mean()),"std":float(cv_df[c].std(ddof=0))} for c in ["r2","rmse","mae"]}
    json_dump(cv_summary,outdir/"historical_10fold_training_cv_summary.json")

    rng=np.random.default_rng(primary); yruns=int(CFG["y_randomization_runs"]); yrows=[]
    if smoke: yruns=min(3,yruns)
    for run_i in range(yruns):
        perm_seed=int(rng.integers(1,2_147_483_647)); yp=np.array(y[tr],copy=True); np.random.default_rng(perm_seed).shuffle(yp)
        model=Pipeline([("scaler",StandardScaler()),("model",_rf(hist_params,primary,smoke))])
        key=f"yrand_{run_i:02d}_{perm_seed}"
        res,_=_run_cached_fit("y_randomization",key,_sig(base_sig,"yrand",run_i,perm_seed),model,Xg[tr],yp,Xg[te],y[te],te,smoke)
        yrows.append({"run":run_i,"shuffle_seed":perm_seed,**res})
    yr=pd.DataFrame(yrows); yr_path=outdir/"y_randomization_metrics.csv"; yr.to_csv(yr_path,index=False)
    yr_summary={c:{"mean":float(yr[c].mean()),"std":float(yr[c].std(ddof=0))} for c in ["r2","rmse","mae"]}
    json_dump(yr_summary,outdir/"y_randomization_summary.json")

    reps={"Graph19":Xg,"RDKit2D":Xrd,"ECFP4_2048":Xm}
    bench_params=CFG["benchmark_rf"]; seeds=CFG["random_seeds"] if not smoke else CFG["random_seeds"][:2]
    scaffold=pd.read_csv(RESULTS/("smoke" if smoke else "full")/"03_representations_splits"/"scaffold_assignments.csv")
    bench_rows=[]; pred_rows=[]
    for repname,X in reps.items():
        for seed in seeds:
            tri,tei=train_test_split(idx,test_size=test_size,random_state=int(seed),shuffle=True)
            model=_rf(bench_params,int(seed),smoke)
            key=f"{repname}_random_seed{seed}"
            sig=_sig(base_sig,"benchmark",bench_params,repname,"random",seed,X.shape)
            res,pred=_run_cached_fit("representation_benchmark",key,sig,model,X[tri],y[tri],X[tei],y[tei],tei,smoke)
            bench_rows.append({"representation":repname,"split_type":"random","split_id":int(seed),"train_n":len(tri),"test_n":len(tei),**res})
            pred["representation"]=repname; pred["split_type"]="random"; pred["split_id"]=int(seed); pred_rows.append(pred)
        folds=sorted(scaffold["fold"].unique())
        if smoke: folds=folds[:2]
        for fold in folds:
            tei=scaffold.index[scaffold["fold"]==fold].to_numpy(); tri=scaffold.index[scaffold["fold"]!=fold].to_numpy()
            model=_rf(bench_params,primary+int(fold),smoke)
            key=f"{repname}_scaffold_fold{int(fold)}"
            sig=_sig(base_sig,"benchmark",bench_params,repname,"scaffold",fold,X.shape)
            res,pred=_run_cached_fit("representation_benchmark",key,sig,model,X[tri],y[tri],X[tei],y[tei],tei,smoke)
            bench_rows.append({"representation":repname,"split_type":"scaffold","split_id":int(fold),"train_n":len(tri),"test_n":len(tei),**res})
            pred["representation"]=repname; pred["split_type"]="scaffold"; pred["split_id"]=int(fold); pred_rows.append(pred)
    bench=pd.DataFrame(bench_rows); bench_path=outdir/"representation_benchmark_metrics.csv"; bench.to_csv(bench_path,index=False)
    preds=pd.concat(pred_rows,ignore_index=True); pred_path=outdir/"representation_benchmark_predictions.csv"; preds.to_csv(pred_path,index=False)
    summary=(bench.groupby(["representation","split_type"])[["r2","rmse","mae","spearman","pearson","train_seconds","predict_seconds"]]
             .agg(["mean","std"]).reset_index())
    summary.columns=["__".join([str(x) for x in c if x]) if isinstance(c,tuple) else c for c in summary.columns]
    summary_path=outdir/"representation_benchmark_summary.csv"; summary.to_csv(summary_path,index=False)
    logger.info("Assurance and representation benchmark phase complete (%d benchmark fits)",len(bench))
    outputs=[assurance_path,cv_path,outdir/"historical_10fold_training_cv_summary.json",yr_path,outdir/"y_randomization_summary.json",bench_path,pred_path,summary_path]
    return outputs,{"assurance":assurance,"cv_summary":cv_summary,"y_randomization":yr_summary,"benchmark_fits":len(bench)}
