@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo First run MASTER_RUN_PAPER003_V3.cmd once so the environment is created.
  pause
  exit /b 1
)
call ".venv\Scripts\activate.bat"
set PYTHONHASHSEED=42
set OMP_NUM_THREADS=1
set MKL_NUM_THREADS=1
set OPENBLAS_NUM_THREADS=1
set NUMEXPR_NUM_THREADS=1
python scripts\run_pipeline.py --smoke-test
pause
