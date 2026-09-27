from __future__ import annotations
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from common import ROOT, INPUT, RESULTS, CFG, sha256_file, json_dump, input_manifest_check


def _out(smoke):
    p=RESULTS/("smoke" if smoke else "full")/"01_audit"
    p.mkdir(parents=True, exist_ok=True)
    return p


def reconstruct_stage1(raw: pd.DataFrame):
    df=raw.copy()
    df["pchembl_value"]=pd.to_numeric(df["pchembl_value"], errors="coerce")
    df["pIC50_from_units"]=pd.to_numeric(df["pIC50_from_units"], errors="coerce")
    df["pIC50_final"]=df["pchembl_value"]
    fb=df["pIC50_final"].isna() & df["pIC50_from_units"].notna()
    df.loc[fb,"pIC50_final"]=df.loc[fb,"pIC50_from_units"]
    df["label_source"]=np.where(df["pchembl_value"].notna(),"pChEMBL",
                         np.where(df["pIC50_from_units"].notna(),"unit-derived","missing"))
    df=df.dropna(subset=["pIC50_final","canonical_smiles"]).copy()
    df=df[~df["canonical_smiles"].astype(str).str.contains(r"\.", regex=True)].copy()
    df=df[np.isfinite(df["pIC50_final"].astype(float)) & (df["pIC50_final"].astype(float)>0)].copy()
    df=df[df["assay_type"]=="B"].copy()
    keep=["molecule_chembl_id","canonical_smiles","pIC50_final","IC50_nM","assay_chembl_id","document_chembl_id"]
    return df, df[keep].copy()


def reconstruct_final(stage1: pd.DataFrame):
    df=stage1.copy()
    def to_key(s):
        try:
            m=Chem.MolFromSmiles(str(s))
            return Chem.MolToInchiKey(m) if m is not None else None
        except Exception: return None
    df["InChIKey"]=df["canonical_smiles"].map(to_key)
    df=df.dropna(subset=["InChIKey"])
    g=(df.groupby("InChIKey", sort=True)
       .agg(canonical_smiles=("canonical_smiles","first"),
            pIC50_median=("pIC50_final","median"),
            pIC50_std=("pIC50_final","std"),
            n_records=("pIC50_final","count")).reset_index())
    g["pIC50_std"]=g["pIC50_std"].fillna(0.0)
    return g


def _compare_frame(a,b,key, numeric_tol=1e-7):
    a=a.sort_values(key).reset_index(drop=True); b=b.sort_values(key).reset_index(drop=True)
    rep={"rows_a":len(a),"rows_b":len(b),"same_row_count":len(a)==len(b)}
    common=[c for c in a.columns if c in b.columns]
    rep["common_columns"]=common
    rep["column_checks"]={}
    for c in common:
        if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            aa=pd.to_numeric(a[c],errors="coerce").to_numpy(float); bb=pd.to_numeric(b[c],errors="coerce").to_numpy(float)
            if len(aa)==len(bb):
                d=np.abs(aa-bb)
                d=d[np.isfinite(d)]
                rep["column_checks"][c]={"max_abs_diff":float(d.max()) if len(d) else 0.0,"within_tol":bool((d<=numeric_tol).all()) if len(d) else True}
        else:
            rep["column_checks"][c]={"exact_equal":bool(a[c].astype(str).equals(b[c].astype(str))) if len(a)==len(b) else False}
    return rep


def standardization_sensitivity(df: pd.DataFrame, outdir: Path, smoke=False):
    if smoke: df=df.head(min(len(df), CFG["smoke_test_rows"]))
    te=rdMolStandardize.TautomerEnumerator()
    uncharger=rdMolStandardize.Uncharger()
    rows=[]
    for i,row in df.iterrows():
        smi=str(row["canonical_smiles"]); oldkey=str(row["InChIKey"])
        rec={"InChIKey":oldkey,"canonical_smiles":smi}
        try:
            m=Chem.MolFromSmiles(smi)
            if m is None: raise ValueError("MolFromSmiles failed")
            clean=rdMolStandardize.Cleanup(m)
            parent=rdMolStandardize.FragmentParent(clean)
            uncharged=uncharger.uncharge(parent)
            taut=te.Canonicalize(uncharged)
            stages={"cleanup":clean,"fragment_parent":parent,"uncharged":uncharged,"canonical_tautomer":taut}
            for name,mol in stages.items():
                ns=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                nk=Chem.MolToInchiKey(mol)
                rec[f"{name}_smiles"]=ns; rec[f"{name}_inchikey"]=nk
                rec[f"{name}_identity_changed"]=(nk!=oldkey)
            rec["error"]=""
        except Exception as e:
            rec["error"]=repr(e)
        rows.append(rec)
    det=pd.DataFrame(rows)
    detail=outdir/"standardization_sensitivity.csv"; det.to_csv(detail,index=False)
    summary={"n":len(det),"errors":int((det.get("error",pd.Series(dtype=str)).astype(str)!="").sum())}
    for name in ["cleanup","fragment_parent","uncharged","canonical_tautomer"]:
        col=f"{name}_identity_changed"
        if col in det:
            summary[f"{name}_changed_n"]=int(det[col].fillna(False).sum())
            summary[f"{name}_changed_fraction"]=float(det[col].fillna(False).mean())
    js=outdir/"standardization_sensitivity_summary.json"; json_dump(summary,js)
    return detail,js,summary


def run(logger, smoke=False):
    outdir=_out(smoke)
    checks=input_manifest_check()
    if not all(v["ok"] for v in checks.values()):
        bad=[k for k,v in checks.items() if not v["ok"]]
        raise RuntimeError(f"Archived V2 input checksum failure: {bad}")
    json_dump(checks,outdir/"input_checksum_audit.json")

    raw=pd.read_csv(INPUT/"egfr_chembl_raw_ic50_binding.csv")
    hist_stage=pd.read_csv(INPUT/"egfr_stage1_clean.csv")
    hist_final=pd.read_csv(INPUT/"egfr_qspr_final.csv")
    if smoke:
        pass
    enriched, rebuilt_stage=reconstruct_stage1(raw)
    rebuilt_final=reconstruct_final(rebuilt_stage)

    rebuilt_stage.to_csv(outdir/"stage1_reconstructed.csv",index=False)
    rebuilt_final.to_csv(outdir/"dataset_frozen.csv",index=False)

    both=(pd.to_numeric(raw["pchembl_value"],errors="coerce")-pd.to_numeric(raw["pIC50_from_units"],errors="coerce")).abs().dropna()
    units=raw["standard_units"].fillna("<missing>").value_counts(dropna=False).to_dict()
    provenance={
        "raw_rows":int(len(raw)),"stage1_rows":int(len(rebuilt_stage)),"final_unique_compounds":int(len(rebuilt_final)),
        "target_ids":sorted(map(str,raw["target_chembl_id"].dropna().unique().tolist())),
        "assay_types":sorted(map(str,raw["assay_type"].dropna().unique().tolist())),
        "standard_types":sorted(map(str,raw["standard_type"].dropna().unique().tolist())),
        "standard_relations":sorted(map(str,raw["standard_relation"].dropna().unique().tolist())),
        "standard_units_counts":{str(k):int(v) for k,v in units.items()},
        "stage1_label_source_counts":{str(k):int(v) for k,v in enriched["label_source"].value_counts().to_dict().items()},
        "records_with_both_pchembl_and_unit_pic50":int(len(both)),
        "pchembl_vs_unit_pic50_absdiff_mean":float(both.mean()),
        "pchembl_vs_unit_pic50_absdiff_median":float(both.median()),
        "pchembl_vs_unit_pic50_absdiff_max":float(both.max()),
        "all_both_sources_agree_within_0.01":bool((both<0.01).all()),
        "endpoint_formula": "pIC50 = -log10(IC50 [M]) = 9 - log10(IC50 [nM])",
        "historical_label_rule": "Use ChEMBL pchembl_value when present; otherwise use pIC50_from_units derived from standard_value/standard_units.",
        "primary_cohort_policy": "V3 assurance preserves the V2-compatible cohort; standardization changes are sensitivity analyses only unless a fatal defect is detected."
    }
    json_dump(provenance,outdir/"provenance_summary.json")

    c1=_compare_frame(rebuilt_stage,hist_stage,"molecule_chembl_id",1e-7)
    c2=_compare_frame(rebuilt_final,hist_final,"InChIKey",1e-7)
    json_dump({"stage1_vs_archived":c1,"final_vs_archived":c2},outdir/"reconstruction_comparison.json")

    enriched[["molecule_chembl_id","canonical_smiles","pIC50_final","label_source","assay_chembl_id","document_chembl_id"]].to_csv(outdir/"label_provenance_stage1.csv",index=False)
    sensitivity_input=rebuilt_final if not smoke else rebuilt_final.head(CFG["smoke_test_rows"])
    sens_detail,sens_json,sens_summary=standardization_sensitivity(sensitivity_input,outdir,smoke=False)

    logger.info("Audit counts: raw=%d stage1=%d final=%d",len(raw),len(rebuilt_stage),len(rebuilt_final))
    logger.info("Label sources after Stage-1: %s", provenance["stage1_label_source_counts"])
    outputs=[outdir/"input_checksum_audit.json",outdir/"stage1_reconstructed.csv",outdir/"dataset_frozen.csv",
             outdir/"provenance_summary.json",outdir/"reconstruction_comparison.json",outdir/"label_provenance_stage1.csv",
             sens_detail,sens_json]
    return outputs,{"counts":provenance,"standardization_sensitivity":sens_summary}
