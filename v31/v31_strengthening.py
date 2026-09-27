from __future__ import annotations
import os, sys, json, time, zipfile, hashlib, shutil, logging, platform
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import t
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'; LOGS=ROOT/'logs'; STATE=ROOT/'state'; CACHE=ROOT/'cache'
for p in [RESULTS,LOGS,STATE,CACHE]: p.mkdir(parents=True,exist_ok=True)
CFG=json.loads((ROOT/'config.json').read_text())
logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s',handlers=[logging.FileHandler(LOGS/'v31.log',mode='a',encoding='utf-8'),logging.StreamHandler(sys.stdout)])
log=logging.getLogger('v31')

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def jobs():
    n=os.cpu_count() or 4
    return max(1,min(int(CFG['n_jobs_cap']),n-int(CFG['reserve_cpu_cores'])))

def metric(y,p):
    return dict(r2=float(r2_score(y,p)),rmse=float(np.sqrt(mean_squared_error(y,p))),mae=float(mean_absolute_error(y,p)))

def prepare_inputs():
    src=ROOT/'inputs'/'PAPER003_V3_RETURN_PACKAGE.zip'
    if not src.exists(): raise FileNotFoundError(src)
    dest=CACHE/'v3_return'; marker=dest/'.extracted_sha256'; sig=sha256(src)
    if not marker.exists() or marker.read_text().strip()!=sig:
        if dest.exists(): shutil.rmtree(dest)
        dest.mkdir(parents=True)
        with zipfile.ZipFile(src) as z: z.extractall(dest)
        marker.write_text(sig)
    return dest/'results'/'full'

def load(base):
    g=pd.read_csv(base/'02_descriptors'/'graph19_regenerated.csv')
    rd=pd.read_csv(base/'03_representations_splits'/'rdkit2d.csv')
    ec=np.load(base/'03_representations_splits'/'morgan_ecfp4_2048.npz')['X'].astype(np.float32,copy=False)
    sc=pd.read_csv(base/'03_representations_splits'/'scaffold_assignments.csv')
    bm=pd.read_csv(base/'04_models'/'representation_benchmark_metrics.csv')
    return g,rd,ec,sc,bm

def collision_analysis(g,ec):
    out=RESULTS/'01_graph19_ambiguity'; out.mkdir(exist_ok=True)
    desc=list(g.columns[3:]); dec=int(CFG['round_decimals_collision'])
    arr=np.round(g[desc].to_numpy(float),dec)
    key=pd.Series([hashlib.sha256(row.tobytes()).hexdigest() for row in arr],index=g.index,name='graph19_key')
    sizes=key.map(key.value_counts()).to_numpy()
    meta=g[['InChIKey','SMILES','pIC50']].copy(); meta['graph19_key']=key; meta['collision_group_size']=sizes
    formulas=[]
    for sm in meta.SMILES:
        m=Chem.MolFromSmiles(sm); formulas.append(rdMolDescriptors.CalcMolFormula(m) if m else '')
    meta['molecular_formula']=formulas
    agg=meta.groupby('graph19_key').agg(n=('pIC50','size'),n_smiles=('SMILES','nunique'),n_formula=('molecular_formula','nunique'),pIC50_min=('pIC50','min'),pIC50_max=('pIC50','max'),pIC50_mean=('pIC50','mean'),pIC50_median=('pIC50','median'),pIC50_std=('pIC50','std')).reset_index()
    agg['pIC50_range']=agg.pIC50_max-agg.pIC50_min
    dup=agg[agg.n>1].copy().sort_values(['pIC50_range','n'],ascending=False)
    dup.to_csv(out/'graph19_collision_groups.csv',index=False); meta.to_csv(out/'graph19_collision_membership.csv',index=False)
    groups=meta.groupby('graph19_key').indices; rows=[]
    for k,idxs in groups.items():
        if len(idxs)<2: continue
        for i,j in combinations(list(idxs),2):
            a=ec[i].astype(bool,copy=False); b=ec[j].astype(bool,copy=False)
            inter=np.logical_and(a,b).sum(); union=np.logical_or(a,b).sum(); sim=float(inter/union) if union else 1.0
            rows.append(dict(graph19_key=k,row_i=int(i),row_j=int(j),tanimoto_ecfp4=sim,delta_pIC50=float(abs(meta.loc[i,'pIC50']-meta.loc[j,'pIC50'])),InChIKey_A=meta.loc[i,'InChIKey'],formula_A=meta.loc[i,'molecular_formula'],SMILES_A=meta.loc[i,'SMILES'],pIC50_A=float(meta.loc[i,'pIC50']),InChIKey_B=meta.loc[j,'InChIKey'],formula_B=meta.loc[j,'molecular_formula'],SMILES_B=meta.loc[j,'SMILES'],pIC50_B=float(meta.loc[j,'pIC50'])))
    pairs=pd.DataFrame(rows); pairs.to_csv(out/'graph19_collision_pairs.csv',index=False)
    strict=pairs[(pairs.tanimoto_ecfp4<float(CFG['ecfp_tanimoto_threshold']))&(pairs.delta_pIC50>float(CFG['activity_delta_threshold']))].copy().sort_values(['delta_pIC50','tanimoto_ecfp4'],ascending=[False,True])
    strict.to_csv(out/'chemically_discordant_collision_pairs.csv',index=False)
    means=meta.groupby('graph19_key').pIC50.transform('mean'); med=meta.groupby('graph19_key').pIC50.transform('median'); y=meta.pIC50.to_numpy(float)
    sse=float(np.sum((y-means.to_numpy())**2)); sst=float(np.sum((y-y.mean())**2)); collided=sizes>1
    summary={'n_molecules':len(meta),'n_unique_graph19_vectors':int(meta.graph19_key.nunique()),'n_collision_groups':int(len(dup)),'n_molecules_in_collision_groups':int(dup.n.sum()),'fraction_molecules_in_collision_groups':float(dup.n.sum()/len(meta)),'collision_groups_with_distinct_molecular_formula':int((dup.n_formula>1).sum()),'fraction_collision_groups_with_distinct_molecular_formula':float((dup.n_formula>1).mean()),'groups_pic50_range_gt_1':int((dup.pIC50_range>1).sum()),'groups_pic50_range_gt_2':int((dup.pIC50_range>2).sum()),'groups_pic50_range_gt_3':int((dup.pIC50_range>3).sum()),'molecules_in_groups_pic50_range_gt_2':int(dup.loc[dup.pIC50_range>2,'n'].sum()),'pairs_tanimoto_lt_0_5_and_delta_pic50_gt_2':int(((pairs.tanimoto_ecfp4<0.5)&(pairs.delta_pIC50>2)).sum()),'exact_collision_dataset_rmse_floor':float(np.sqrt(np.mean((y-means.to_numpy())**2))),'exact_collision_dataset_mae_floor':float(np.mean(np.abs(y-med.to_numpy()))),'exact_collision_in_sample_r2_ceiling':float(1-sse/sst),'collided_subset_rmse_floor':float(np.sqrt(np.mean((y[collided]-means.to_numpy()[collided])**2))),'collided_subset_mae_floor':float(np.mean(np.abs(y[collided]-med.to_numpy()[collided])))}
    (out/'graph19_ambiguity_summary.json').write_text(json.dumps(summary,indent=2))
    n=min(int(CFG['representative_pairs_to_draw']),len(strict)); strict.head(n).to_csv(out/'representative_collision_pairs.csv',index=False)
    log.info('Graph19 ambiguity: %s',summary); return summary

def paired_effects(bm):
    out=RESULTS/'02_paired_effects'; out.mkdir(exist_ok=True); rows=[]
    comps=[('RDKit2D','Graph19'),('ECFP4_2048','Graph19'),('ECFP4_2048','RDKit2D')]
    for st in ['random','scaffold']:
      for a,b in comps:
        aa=bm[(bm.representation==a)&(bm.split_type==st)].sort_values('split_id'); bb=bm[(bm.representation==b)&(bm.split_type==st)].sort_values('split_id'); assert aa.split_id.tolist()==bb.split_id.tolist()
        for met in ['r2','rmse','mae']:
            d=aa[met].to_numpy()-bb[met].to_numpy(); n=len(d); m=float(d.mean()); sd=float(d.std(ddof=1)); se=sd/np.sqrt(n); q=float(t.ppf(.975,n-1))
            rows.append(dict(split_type=st,representation_A=a,representation_B=b,metric=met,n_pairs=n,mean_delta_A_minus_B=m,sd_delta=sd,ci95_low=m-q*se,ci95_high=m+q*se,n_positive=int((d>0).sum()),n_negative=int((d<0).sum()),deltas=';'.join(f'{x:.8f}' for x in d)))
    pd.DataFrame(rows).to_csv(out/'paired_representation_effect_sizes.csv',index=False); return rows

def fit_one(X,y,tr,te,rep,st,sid,mode,max_features,trees,random_state):
    out=RESULTS/'03_rf_sensitivity'; out.mkdir(exist_ok=True); cache=out/f'{mode}__{rep}__{st}__{sid}.json'
    if cache.exists(): return json.loads(cache.read_text())
    model=RandomForestRegressor(n_estimators=int(trees),max_features=max_features,min_samples_leaf=1,random_state=int(random_state),n_jobs=jobs())
    t0=time.perf_counter(); model.fit(X[tr],y[tr]); train_s=time.perf_counter()-t0; t0=time.perf_counter(); p=model.predict(X[te]); pred_s=time.perf_counter()-t0
    r={'mode':mode,'representation':rep,'split_type':st,'split_id':int(sid),'max_features':max_features,'n_estimators':int(trees),'train_n':int(len(tr)),'test_n':int(len(te)),**metric(y[te],p),'train_seconds':train_s,'predict_seconds':pred_s}
    cache.write_text(json.dumps(r,indent=2)); log.info('FIT %s %s %s %s R2=%.4f train=%.1fs',mode,rep,st,sid,r['r2'],train_s); return r

def rf_sensitivity(g,rd,ec,sc):
    y=g.pIC50.to_numpy(float); idx=np.arange(len(y)); reps={'Graph19':g.iloc[:,3:].to_numpy(float),'RDKit2D':rd.iloc[:,3:].to_numpy(float),'ECFP4_2048':ec}; rows=[]
    frac=float(CFG['rf_fraction_max_features']); trees=int(CFG['rf_trees_fraction_sensitivity'])
    for rep,X in reps.items():
      for seed in CFG['random_seeds']:
        tr,te=train_test_split(idx,test_size=.2,random_state=int(seed),shuffle=True); rows.append(fit_one(X,y,tr,te,rep,'random',seed,'fraction_0.5',frac,trees,int(seed)))
      for fold in sorted(sc['fold'].unique()):
        te=sc.index[sc['fold']==fold].to_numpy(); tr=sc.index[sc['fold']!=fold].to_numpy(); rows.append(fit_one(X,y,tr,te,rep,'scaffold',int(fold),'fraction_0.5',frac,trees,42+int(fold)))
    trees2=int(CFG['rf_trees_full_feature_spotcheck'])
    for rep,X in reps.items():
      for seed in CFG['rf_full_feature_spotcheck_random_seeds']:
        tr,te=train_test_split(idx,test_size=.2,random_state=int(seed),shuffle=True); rows.append(fit_one(X,y,tr,te,rep,'random',seed,'full_feature_spotcheck',1.0,trees2,int(seed)))
      for fold in CFG['rf_full_feature_spotcheck_scaffold_folds']:
        te=sc.index[sc['fold']==fold].to_numpy(); tr=sc.index[sc['fold']!=fold].to_numpy(); rows.append(fit_one(X,y,tr,te,rep,'scaffold',int(fold),'full_feature_spotcheck',1.0,trees2,42+int(fold)))
    df=pd.DataFrame(rows); out=RESULTS/'03_rf_sensitivity'; df.to_csv(out/'rf_sensitivity_metrics.csv',index=False)
    sm=df.groupby(['mode','representation','split_type'])[['r2','rmse','mae','train_seconds']].agg(['mean','std','count']).reset_index(); sm.columns=['__'.join([str(v) for v in c if v]) if isinstance(c,tuple) else c for c in sm.columns]; sm.to_csv(out/'rf_sensitivity_summary.csv',index=False); return df

def package(base,ambiguity):
    env={'python':sys.version,'platform':platform.platform(),'cpu_count':os.cpu_count(),'n_jobs_used':jobs(),'v3_return_sha256':sha256(ROOT/'inputs'/'PAPER003_V3_RETURN_PACKAGE.zip')}; (RESULTS/'environment.json').write_text(json.dumps(env,indent=2)); (RESULTS/'RUN_COMPLETE.txt').write_text('PAPER003 V3.1 strengthening completed successfully.\n')
    manifest=[]
    for p in sorted(RESULTS.rglob('*')):
        if p.is_file() and p.name!='PAPER003_V31_RETURN_PACKAGE.zip': manifest.append({'path':str(p.relative_to(RESULTS)).replace('\\','/'),'sha256':sha256(p),'bytes':p.stat().st_size})
    (RESULTS/'results_manifest.json').write_text(json.dumps(manifest,indent=2))
    zpath=RESULTS/'PAPER003_V31_RETURN_PACKAGE.zip'
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(RESULTS.rglob('*')):
            if p.is_file() and p!=zpath: z.write(p,p.relative_to(RESULTS))
        z.write(LOGS/'v31.log',Path('logs')/'v31.log'); z.write(ROOT/'config.json','config.json'); z.write(ROOT/'README_RUN_FIRST.md','README_RUN_FIRST.md')
    log.info('Return package: %s SHA256=%s',zpath,sha256(zpath))

def main():
    log.info('PAPER003 V3.1 started; CPU=%s jobs=%s',os.cpu_count(),jobs()); base=prepare_inputs(); g,rd,ec,sc,bm=load(base); assert len(g)==10056 and len(rd)==10056 and ec.shape==(10056,2048); amb=collision_analysis(g,ec); paired_effects(bm); rf_sensitivity(g,rd,ec,sc); package(base,amb); log.info('PAPER003 V3.1 COMPLETE')
if __name__=='__main__': main()
