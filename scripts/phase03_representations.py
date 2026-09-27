from __future__ import annotations
from pathlib import Path
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import train_test_split, GroupKFold
from common import RESULTS, CACHE, CFG, jobs_from_cfg, json_dump


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"03_representations_splits"; p.mkdir(parents=True,exist_ok=True); return p

RDKIT_NAMES=[name for name,_ in Descriptors._descList]
RDKIT_FUNCS=dict(Descriptors._descList)

def _rdkit2d_one(item):
    idx,key,smi,y=item
    mol=Chem.MolFromSmiles(smi)
    if mol is None: raise ValueError(f"Invalid SMILES row {idx}")
    vals=[]
    for name in RDKIT_NAMES:
        try:
            v=float(RDKIT_FUNCS[name](mol))
            vals.append(v if np.isfinite(v) else np.nan)
        except Exception:
            vals.append(np.nan)
    return [idx,key,smi,float(y),*vals]


def run(logger, smoke=False):
    outdir=_out(smoke)
    graph_path=RESULTS/("smoke" if smoke else "full")/"02_descriptors"/"graph19_regenerated.csv"
    graph=pd.read_csv(graph_path)
    if smoke: graph=graph.head(min(len(graph),CFG["smoke_test_rows"])).copy()
    n=len(graph); logger.info("Generating contemporary representations for %d compounds",n)

    t0=time.perf_counter()
    jobs=jobs_from_cfg("descriptor_workers") if not smoke else min(2,jobs_from_cfg("descriptor_workers"))
    tasks=[(i,r.InChIKey,str(r.SMILES),float(r.pIC50)) for i,r in graph.iterrows()]
    if jobs>1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            rows=list(ex.map(_rdkit2d_one,tasks,chunksize=20))
    else: rows=[_rdkit2d_one(x) for x in tasks]
    all_cols=["row_index","InChIKey","SMILES","pIC50",*RDKIT_NAMES]
    rd=pd.DataFrame(rows,columns=all_cols).sort_values("row_index").drop(columns="row_index")
    numeric=rd[RDKIT_NAMES].replace([np.inf,-np.inf],np.nan)
    keep=[]; dropped={}
    for c in RDKIT_NAMES:
        miss=float(numeric[c].isna().mean())
        nun=int(numeric[c].nunique(dropna=True))
        if miss>0:
            dropped[c]=f"nonfinite_or_missing_fraction={miss:.6g}"
        elif nun<=1:
            dropped[c]="constant"
        else: keep.append(c)
    rd=rd[["InChIKey","SMILES","pIC50",*keep]]
    rd_path=outdir/"rdkit2d.csv"; rd.to_csv(rd_path,index=False)
    rd_meta={"available_descriptor_names":RDKIT_NAMES,"retained_descriptor_names":keep,"dropped":dropped,
             "n_retained":len(keep),"generation_seconds":time.perf_counter()-t0}
    rd_meta_path=outdir/"rdkit2d_metadata.json"; json_dump(rd_meta,rd_meta_path)

    t0=time.perf_counter(); radius=int(CFG["morgan"]["radius"]); nbits=int(CFG["morgan"]["n_bits"])
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=radius,fpSize=nbits,includeChirality=bool(CFG["morgan"].get("use_chirality",False)))
    sparse_gen=rdFingerprintGenerator.GetMorganGenerator(radius=radius,includeChirality=bool(CFG["morgan"].get("use_chirality",False)))
    X=np.zeros((n,nbits),dtype=np.uint8); coll=[]
    for i,smi in enumerate(graph["SMILES"].astype(str)):
        mol=Chem.MolFromSmiles(smi)
        fp=gen.GetFingerprint(mol)
        arr=np.zeros((nbits,),dtype=np.uint8); DataStructs.ConvertToNumpyArray(fp,arr); X[i]=arr
        ids=list(sparse_gen.GetSparseCountFingerprint(mol).GetNonzeroElements().keys())
        folded={int(x)%nbits for x in ids}
        collision_count=max(0,len(ids)-len(folded))
        coll.append({"row_index":i,"InChIKey":graph.iloc[i]["InChIKey"],"unfolded_unique_identifiers":len(ids),
                     "folded_unique_bins":len(folded),"within_molecule_hash_collisions":collision_count,
                     "collision_fraction":float(collision_count/len(ids)) if ids else 0.0})
    morgan_path=outdir/"morgan_ecfp4_2048.npz"; np.savez_compressed(morgan_path,X=X)
    coll_df=pd.DataFrame(coll); coll_path=outdir/"morgan_collision_audit.csv"; coll_df.to_csv(coll_path,index=False)
    collision_meta={"radius":radius,"n_bits":nbits,"binary":True,"use_chirality":bool(CFG["morgan"].get("use_chirality",False)),
                    "n_compounds":n,"molecules_with_any_hash_collision":int((coll_df["within_molecule_hash_collisions"]>0).sum()),
                    "fraction_molecules_with_any_hash_collision":float((coll_df["within_molecule_hash_collisions"]>0).mean()),
                    "median_collision_fraction":float(coll_df["collision_fraction"].median()),
                    "mean_collision_fraction":float(coll_df["collision_fraction"].mean()),
                    "max_collision_fraction":float(coll_df["collision_fraction"].max()),
                    "generation_seconds":time.perf_counter()-t0,
                    "interpretation":"Collision audit counts unique sparse Morgan identifiers mapping to the same hashed bit within each molecule. It is diagnostic; the 2048-bit ECFP4 representation remains the recognizable literature benchmark."}
    collision_meta_path=outdir/"morgan_metadata.json"; json_dump(collision_meta,collision_meta_path)

    split_rows=[]; indices=np.arange(n)
    for seed in CFG["random_seeds"]:
        tr,te=train_test_split(indices,test_size=float(CFG["test_size"]),random_state=int(seed),shuffle=True)
        split_rows += [{"seed":seed,"row_index":int(i),"InChIKey":graph.iloc[i]["InChIKey"],"split":"train"} for i in tr]
        split_rows += [{"seed":seed,"row_index":int(i),"InChIKey":graph.iloc[i]["InChIKey"],"split":"test"} for i in te]
    random_splits=pd.DataFrame(split_rows); random_path=outdir/"random_split_assignments.csv"; random_splits.to_csv(random_path,index=False)

    scaff=[]
    for i,smi in enumerate(graph["SMILES"].astype(str)):
        scaff.append(MurckoScaffold.MurckoScaffoldSmiles(smiles=smi,includeChirality=False))
    scaffold_df=pd.DataFrame({"row_index":np.arange(n),"InChIKey":graph["InChIKey"],"scaffold":scaff})
    scaffold_df["fold"]=-1
    gkf=GroupKFold(n_splits=int(CFG["scaffold_folds"]) if n>=int(CFG["scaffold_folds"]) else 2)
    dummy=np.zeros((n,1)); y=graph["pIC50"].to_numpy(float)
    for fold,(_,te) in enumerate(gkf.split(dummy,y,groups=np.array(scaff))): scaffold_df.loc[te,"fold"]=fold
    scaffold_path=outdir/"scaffold_assignments.csv"; scaffold_df.to_csv(scaffold_path,index=False)
    scaffold_meta={"n_compounds":n,"unique_murcko_scaffolds":int(scaffold_df["scaffold"].nunique(dropna=False)),
                   "empty_scaffold_count":int((scaffold_df["scaffold"]=="").sum()),
                   "fold_sizes":{str(k):int(v) for k,v in scaffold_df["fold"].value_counts().sort_index().to_dict().items()},
                   "method":"RDKit Bemis-Murcko scaffold (includeChirality=False) + sklearn GroupKFold, no shuffle"}
    scaffold_meta_path=outdir/"scaffold_metadata.json"; json_dump(scaffold_meta,scaffold_meta_path)
    if not smoke and scaffold_meta["unique_murcko_scaffolds"]!=3579:
        raise RuntimeError(f"Scaffold assurance gate failed; expected 3579, got {scaffold_meta['unique_murcko_scaffolds']}")

    rep_summary={"graph19_dimensions":[n,19],"rdkit2d_dimensions":[n,len(keep)],"morgan_dimensions":[n,nbits],
                 "scaffolds":scaffold_meta,"morgan_collision_audit":collision_meta}
    summary_path=outdir/"representation_split_summary.json"; json_dump(rep_summary,summary_path)
    logger.info("Representations ready: graph19=19, RDKit2D=%d, Morgan=%d; scaffolds=%d",len(keep),nbits,scaffold_meta["unique_murcko_scaffolds"])
    return [rd_path,rd_meta_path,morgan_path,coll_path,collision_meta_path,random_path,scaffold_path,scaffold_meta_path,summary_path],rep_summary
