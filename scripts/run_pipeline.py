from __future__ import annotations
import os, sys, argparse, traceback, shutil, multiprocessing
from pathlib import Path

# Constrain hidden BLAS/OpenMP oversubscription; RF/process parallelism is controlled explicitly.
os.environ.setdefault("PYTHONHASHSEED","42")
os.environ.setdefault("OMP_NUM_THREADS","1")
os.environ.setdefault("MKL_NUM_THREADS","1")
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("NUMEXPR_NUM_THREADS","1")

from common import ROOT, RESULTS, STATE, CFG, setup_logging, phase_is_done, mark_phase_done, json_dump, environment_report, input_manifest_check, pip_freeze, package_signature
import phase01_audit, phase02_descriptors, phase03_representations, phase04_models, phase05_generalization_ad, phase06_strengthening, phase07_outputs_qa

PHASES=[
 ("00_preflight",None),
 ("01_audit",phase01_audit.run),
 ("02_descriptors",phase02_descriptors.run),
 ("03_representations_splits",phase03_representations.run),
 ("04_models",phase04_models.run),
 ("05_generalization_ad",phase05_generalization_ad.run),
 ("06_strengthening",phase06_strengthening.run),
 ("07_outputs_qa",phase07_outputs_qa.run),
]

def preflight(logger,smoke=False):
    base=RESULTS/("smoke" if smoke else "full")/"00_preflight"; base.mkdir(parents=True,exist_ok=True)
    checks=input_manifest_check();
    if not all(v["ok"] for v in checks.values()): raise RuntimeError("Input checksum audit failed. V2 archive has changed or is incomplete.")
    env=environment_report(); disk=shutil.disk_usage(ROOT)
    env["disk_free_gb"]=disk.free/(1024**3); env["package_signature"]=package_signature(); env["smoke_test"]=smoke
    if disk.free < 5*(1024**3): raise RuntimeError("Less than 5 GB free disk space. Free space before running EGFR Graph QSAR V3.")
    p=base/"preflight.json"; json_dump({"environment":env,"inputs":checks},p)
    f=base/"pip_freeze.txt"; pip_freeze(f)
    logger.info("Preflight OK: CPUs=%s, free disk=%.1f GB",env.get("cpu_count"),env["disk_free_gb"])
    return [p,f],env

def main():
    ap=argparse.ArgumentParser(description="EGFR Graph QSAR V3 assurance + strengthening pipeline")
    ap.add_argument("--smoke-test",action="store_true",help="Run reduced non-scientific QA mode; outputs are isolated under results/smoke.")
    ap.add_argument("--force-phase",default=None,help="Force one named phase to rerun by deleting its completion marker.")
    args=ap.parse_args(); smoke=args.smoke_test
    logger=setup_logging(smoke)
    logger.info("EGFR Graph QSAR V3 pipeline starting | smoke=%s | root=%s",smoke,ROOT)
    logger.info("Design: preserve V2 provenance; no manuscript rewrite is performed by this execution package.")
    if args.force_phase:
        from common import phase_state_path
        p=phase_state_path(args.force_phase,smoke)
        if p.exists(): p.unlink(); logger.info("Removed completion marker for forced phase %s",args.force_phase)
    try:
        for name,func in PHASES:
            if phase_is_done(name,smoke):
                logger.info("SKIP %s: validated checkpoint exists",name); continue
            logger.info("START %s",name)
            if name=="00_preflight": outputs,meta=preflight(logger,smoke)
            else: outputs,meta=func(logger,smoke)
            mark_phase_done(name,outputs,smoke,meta)
            logger.info("DONE %s",name)
        logger.info("EGFR Graph QSAR V3 PIPELINE COMPLETED SUCCESSFULLY")
        final=RESULTS/("smoke" if smoke else "full")/"07_final"/"EGFR_QSAR_REPRODUCIBILITY_RETURN.zip"
        if final.exists(): logger.info("UPLOAD BACK TO CHATGPT: %s",final)
    except KeyboardInterrupt:
        logger.warning("Interrupted safely. Completed phase/chunk/experiment checkpoints are preserved; rerun the CMD to resume.")
        raise
    except Exception as e:
        logger.error("PIPELINE FAILED: %s",e)
        logger.error(traceback.format_exc())
        fail=STATE/("smoke_LAST_FAILURE.txt" if smoke else "LAST_FAILURE.txt")
        fail.write_text(traceback.format_exc(),encoding="utf-8")
        raise

if __name__=="__main__":
    multiprocessing.freeze_support()
    main()
