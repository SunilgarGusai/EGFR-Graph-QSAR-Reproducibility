from __future__ import annotations
from pathlib import Path
import math, os
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from scipy.sparse.csgraph import shortest_path
from scipy.linalg import eigvalsh
from common import RESULTS, CACHE, INPUT, CFG, jobs_from_cfg, sha256_file, json_dump

DESC_COLS=["M1","M2","Randic","Harmonic","ABC","GA","AZI","F","HM","SCI",
           "Wiener","Harary","Szeged","Schultz","Gutman","EdgeWiener","VertexEdgeWiener","Energy","Estrada"]


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"02_descriptors"; p.mkdir(parents=True,exist_ok=True); return p


def calculate_one(item):
    idx, inchikey, smiles, y = item
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError(f"Invalid SMILES at row {idx}: {smiles}")
    A=Chem.GetAdjacencyMatrix(mol).astype(float)
    n=A.shape[0]
    if n < 1: raise ValueError(f"Empty graph at row {idx}")
    deg=A.sum(axis=1)
    edges=np.argwhere(np.triu(A,1)>0)
    if n>1 and len(edges)==0: raise ValueError(f"Disconnected/edgeless molecular graph at row {idx}")
    M1=float(np.sum(deg**2)); M2=float(sum(deg[u]*deg[v] for u,v in edges))
    Randic=float(sum(1/np.sqrt(deg[u]*deg[v]) for u,v in edges))
    Harmonic=float(sum(2/(deg[u]+deg[v]) for u,v in edges))
    ABC=float(sum(np.sqrt((deg[u]+deg[v]-2)/(deg[u]*deg[v])) for u,v in edges))
    GA=float(sum(2*np.sqrt(deg[u]*deg[v])/(deg[u]+deg[v]) for u,v in edges))
    AZI=0.0
    for u,v in edges:
        den=deg[u]+deg[v]-2
        if den!=0: AZI += float((deg[u]*deg[v]/den)**3)
    F=float(np.sum(deg**3)); HM=float(sum((deg[u]+deg[v])**2 for u,v in edges))
    SCI=float(sum(1/np.sqrt(deg[u]+deg[v]) for u,v in edges))
    if n==1:
        D=np.zeros((1,1),float)
    else:
        D=shortest_path(A,directed=False,unweighted=True)
    if not np.isfinite(D).all(): raise ValueError(f"Non-connected graph at row {idx}")
    iu=np.triu_indices(n,1)
    Wiener=float(D[iu].sum())
    Harary=float((1/D[iu]).sum()) if len(iu[0]) else 0.0
    Szeged=0.0
    for u,v in edges:
        du=D[:,u]; dv=D[:,v]
        Szeged += float(np.sum(du<dv)*np.sum(dv<du))
    Schultz=0.0; Gutman=0.0
    for i in range(n):
        for j in range(i+1,n):
            dij=D[i,j]
            Schultz += float((deg[i]+deg[j])*dij)
            Gutman += float(deg[i]*deg[j]*dij)
    m=len(edges); EdgeWiener=0.0
    for a in range(m):
        u,v=edges[a]
        for b in range(a+1,m):
            x,z=edges[b]
            EdgeWiener += float(min(D[u,x],D[u,z],D[v,x],D[v,z])+1)
    VertexEdgeWiener=0.0
    for i in range(n):
        for u,v in edges:
            VertexEdgeWiener += float(min(D[i,u],D[i,v]))
    lam=eigvalsh(A)
    Energy=float(np.abs(lam).sum()); Estrada=float(np.exp(lam).sum())
    vals=[M1,M2,Randic,Harmonic,ABC,GA,AZI,F,HM,SCI,Wiener,Harary,Szeged,Schultz,Gutman,EdgeWiener,VertexEdgeWiener,Energy,Estrada]
    return [idx,inchikey,smiles,float(y),*vals]


def _chunk_valid(path: Path, expected_n: int):
    if not path.exists(): return False
    try:
        d=pd.read_csv(path)
        return len(d)==expected_n and list(d.columns)==["row_index","InChIKey","SMILES","pIC50",*DESC_COLS]
    except Exception: return False


def run(logger, smoke=False):
    outdir=_out(smoke)
    frozen=RESULTS/("smoke" if smoke else "full")/"01_audit"/"dataset_frozen.csv"
    df=pd.read_csv(frozen)
    if smoke: df=df.head(min(len(df),CFG["smoke_test_rows"])).copy()
    cache_dir=CACHE/("smoke" if smoke else "full")/"descriptors"; cache_dir.mkdir(parents=True,exist_ok=True)
    chunk_size=CFG["descriptor_chunk_size"] if not smoke else min(80,CFG["descriptor_chunk_size"])
    jobs=jobs_from_cfg("descriptor_workers") if not smoke else min(2,jobs_from_cfg("descriptor_workers"))
    logger.info("Regenerating 19 descriptors for %d compounds using %d workers",len(df),jobs)
    cols=["row_index","InChIKey","SMILES","pIC50",*DESC_COLS]
    chunk_paths=[]
    for start in range(0,len(df),chunk_size):
        stop=min(start+chunk_size,len(df)); cp=cache_dir/f"chunk_{start:06d}_{stop:06d}.csv"; chunk_paths.append(cp)
        if _chunk_valid(cp,stop-start):
            logger.info("Descriptor chunk cached: %d:%d",start,stop); continue
        sub=df.iloc[start:stop]
        tasks=[(int(i),r.InChIKey,str(r.canonical_smiles),float(r.pIC50_median)) for i,r in sub.iterrows()]
        if jobs>1:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                rows=list(ex.map(calculate_one,tasks,chunksize=10))
        else:
            rows=[calculate_one(x) for x in tasks]
        pd.DataFrame(rows,columns=cols).to_csv(cp,index=False)
        logger.info("Descriptor chunk complete: %d:%d",start,stop)
    regen=pd.concat([pd.read_csv(p) for p in chunk_paths],ignore_index=True).sort_values("row_index")
    regen=regen.drop(columns=["row_index"])
    regen_path=outdir/"graph19_regenerated.csv"; regen.to_csv(regen_path,index=False)

    hist=pd.read_csv(INPUT/"egfr_qspr_descriptors.csv")
    if smoke: hist=hist[hist["InChIKey"].isin(regen["InChIKey"])].copy()
    merged=regen.merge(hist,on="InChIKey",suffixes=("_v3","_v2"),how="inner")
    rows=[]; tol=float(CFG["descriptor_tolerance"])
    for c in DESC_COLS:
        d=np.abs(merged[f"{c}_v3"].to_numpy(float)-merged[f"{c}_v2"].to_numpy(float))
        rows.append({"descriptor":c,"n":len(d),"max_abs_diff":float(d.max()),"mean_abs_diff":float(d.mean()),"within_tolerance":bool((d<=tol).all())})
    comp=pd.DataFrame(rows); comp_path=outdir/"descriptor_comparison.csv"; comp.to_csv(comp_path,index=False)
    summary={"n_regenerated":int(len(regen)),"n_matched_to_v2":int(len(merged)),"tolerance":tol,
             "all_19_within_tolerance":bool(comp["within_tolerance"].all()),
             "global_max_abs_diff":float(comp["max_abs_diff"].max()),"descriptor_columns":DESC_COLS}
    summary_path=outdir/"descriptor_comparison_summary.json"; json_dump(summary,summary_path)
    if not smoke and (len(regen)!=10056 or not summary["all_19_within_tolerance"]):
        raise RuntimeError(f"Descriptor assurance gate failed: {summary}")
    logger.info("Descriptor assurance: %s",summary)
    return [regen_path,comp_path,summary_path],summary
