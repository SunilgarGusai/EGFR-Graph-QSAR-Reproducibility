from __future__ import annotations
import os, sys, json, hashlib, logging, platform, subprocess, time
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "inputs" / "v2_archive"
STATE = ROOT / "state"
CACHE = ROOT / "cache"
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"
DOCS = ROOT / "docs"
for p in (STATE, CACHE, RESULTS, LOGS): p.mkdir(parents=True, exist_ok=True)


def load_config():
    return json.loads((ROOT / "config.json").read_text(encoding="utf-8"))

CFG = load_config()


def auto_jobs(reserve: int = 2, cap: int = 10) -> int:
    n = os.cpu_count() or 4
    return max(1, min(cap, n - reserve if n > reserve else 1))


def jobs_from_cfg(key="model_jobs") -> int:
    v = CFG.get(key, "auto")
    return auto_jobs() if v == "auto" else int(v)


def sha256_file(path: Path | str) -> str:
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def json_dump(obj, path: Path | str):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str), encoding="utf-8")


def setup_logging(smoke=False):
    log = LOGS / ("pipeline_smoke.log" if smoke else "pipeline.log")
    logger = logging.getLogger("paper003")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(log, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh); logger.addHandler(sh)
    return logger


def package_signature() -> str:
    h = hashlib.sha256()
    for rel in ["config.json", "requirements.txt", "INPUT_MANIFEST.json"]:
        p = ROOT / rel
        h.update(p.read_bytes())
    for p in sorted((ROOT / "scripts").glob("*.py")):
        h.update(p.name.encode()); h.update(p.read_bytes())
    return h.hexdigest()


def phase_state_path(name: str, smoke=False) -> Path:
    prefix = "smoke_" if smoke else ""
    return STATE / f"{prefix}{name}.done.json"


def outputs_valid(state: dict) -> bool:
    for rec in state.get("outputs", []):
        p = ROOT / rec["path"]
        if not p.exists() or sha256_file(p) != rec.get("sha256"):
            return False
    return state.get("package_signature") == package_signature()


def phase_is_done(name: str, smoke=False) -> bool:
    p = phase_state_path(name, smoke)
    if not p.exists(): return False
    try:
        return outputs_valid(json.loads(p.read_text(encoding="utf-8")))
    except Exception:
        return False


def mark_phase_done(name: str, outputs: Iterable[Path], smoke=False, metadata=None):
    recs=[]
    for p in outputs:
        p=Path(p)
        if not p.exists():
            raise FileNotFoundError(f"Expected phase output missing: {p}")
        recs.append({"path": str(p.relative_to(ROOT)), "sha256": sha256_file(p), "bytes": p.stat().st_size})
    payload={
        "phase": name,
        "completed_epoch": time.time(),
        "package_signature": package_signature(),
        "outputs": recs,
        "metadata": metadata or {}
    }
    json_dump(payload, phase_state_path(name, smoke))


def experiment_cache_path(category: str, key: str, smoke=False) -> Path:
    base = CACHE / ("smoke" if smoke else "full") / category
    base.mkdir(parents=True, exist_ok=True)
    return base / f"{key}.json"


def cache_load(path: Path, signature: str):
    if not path.exists(): return None
    try:
        obj=json.loads(path.read_text(encoding="utf-8"))
        if obj.get("signature") == signature:
            return obj
    except Exception:
        pass
    return None


def cache_save(path: Path, signature: str, result: dict):
    payload={"signature": signature, **result}
    json_dump(payload, path)


def environment_report() -> dict:
    rep={
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "python_executable": sys.executable,
    }
    try:
        import numpy, pandas, scipy, sklearn, matplotlib, psutil, joblib
        from rdkit import rdBase
        rep.update({
            "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__, "rdkit": rdBase.rdkitVersion,
            "matplotlib": matplotlib.__version__, "psutil": psutil.__version__, "joblib": joblib.__version__
        })
    except Exception as e:
        rep["dependency_probe_error"] = repr(e)
    return rep


def pip_freeze(path: Path):
    try:
        out=subprocess.check_output([sys.executable,"-m","pip","freeze"], text=True, stderr=subprocess.STDOUT)
        path.write_text(out, encoding="utf-8")
    except Exception as e:
        path.write_text(f"pip freeze failed: {e}\n", encoding="utf-8")


def input_manifest_check() -> dict:
    expected=json.loads((ROOT/"INPUT_MANIFEST.json").read_text(encoding="utf-8"))
    out={}
    for name, want in expected.items():
        p=INPUT/name
        got=sha256_file(p) if p.exists() else None
        out[name]={"expected":want,"actual":got,"ok":got==want}
    return out


def rmse(y_true, y_pred):
    import numpy as np
    return float(np.sqrt(np.mean((np.asarray(y_true)-np.asarray(y_pred))**2)))


def metric_dict(y_true, y_pred):
    import numpy as np
    from sklearn.metrics import r2_score, mean_absolute_error
    from scipy.stats import spearmanr, pearsonr
    y_true=np.asarray(y_true,float); y_pred=np.asarray(y_pred,float)
    try: sp=float(spearmanr(y_true,y_pred).statistic)
    except Exception: sp=float("nan")
    try: pr=float(pearsonr(y_true,y_pred).statistic)
    except Exception: pr=float("nan")
    return {"r2":float(r2_score(y_true,y_pred)),"rmse":rmse(y_true,y_pred),
            "mae":float(mean_absolute_error(y_true,y_pred)),"spearman":sp,"pearson":pr,"n":int(len(y_true))}
