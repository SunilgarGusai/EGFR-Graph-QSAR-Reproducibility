@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"
title EGFR Graph QSAR V3 - Assurance and Strengthening

set PYTHONHASHSEED=42
set OMP_NUM_THREADS=1
set MKL_NUM_THREADS=1
set OPENBLAS_NUM_THREADS=1
set NUMEXPR_NUM_THREADS=1

if not exist logs mkdir logs
set "BOOTLOG=logs\bootstrap.log"

echo ============================================================
echo EGFR Graph QSAR V3 - ONE-CLICK LOCAL EXECUTION
echo V2 is preserved. This package performs assurance + strengthening.
echo ============================================================
echo.
echo [%date% %time%] Bootstrap starting>>"%BOOTLOG%"

set "PYBOOT="
py -3.11 -c "import sys; assert sys.version_info >= (3,11)" >nul 2>&1 && set "PYBOOT=py -3.11"
if not defined PYBOOT py -3.12 -c "import sys; assert sys.version_info >= (3,11)" >nul 2>&1 && set "PYBOOT=py -3.12"
if not defined PYBOOT py -3.13 -c "import sys; assert sys.version_info >= (3,11)" >nul 2>&1 && set "PYBOOT=py -3.13"
if not defined PYBOOT python -c "import sys; assert (3,11) <= sys.version_info[:2] <= (3,13)" >nul 2>&1 && set "PYBOOT=python"

if not defined PYBOOT (
  echo ERROR: Python 3.11-3.13 x64 is required but was not found.
  echo Install Python from python.org, enable Add Python to PATH, then double-click this file again.
  pause
  exit /b 1
)

echo Using bootstrap Python: %PYBOOT%
if not exist ".venv\Scripts\python.exe" (
  echo Creating isolated virtual environment...
  call %PYBOOT% -m venv .venv >>"%BOOTLOG%" 2>&1
  if errorlevel 1 goto :FAIL
)
call ".venv\Scripts\activate.bat"
python -m pip install --upgrade pip >>"%BOOTLOG%" 2>&1

set "REQHASHFILE=.venv\.egfr_qsar_requirements_hash.txt"
for /f "usebackq delims=" %%H in (`python -c "import hashlib; print(hashlib.sha256(open('requirements.txt','rb').read()).hexdigest())"`) do set "REQHASH=%%H"
set "OLDHASH="
if exist "%REQHASHFILE%" set /p OLDHASH=<"%REQHASHFILE%"
if /I "!OLDHASH!"=="!REQHASH!" goto :DEPS_READY

echo Installing scientific dependencies. Internet is needed only for first setup or dependency changes.
set /a ATTEMPT=0
:DEP_RETRY
set /a ATTEMPT+=1
echo Dependency install attempt !ATTEMPT! of 3...
python -m pip install -r requirements.txt >>"%BOOTLOG%" 2>&1
if not errorlevel 1 goto :DEP_OK
if !ATTEMPT! GEQ 3 goto :FAIL
timeout /t 10 /nobreak >nul
goto :DEP_RETRY
:DEP_OK
>"%REQHASHFILE%" echo !REQHASH!
:DEPS_READY
echo Dependencies validated.
echo.
echo Starting resumable scientific pipeline...
echo You may safely interrupt; rerunning this CMD resumes from validated checkpoints.
python scripts\run_pipeline.py
if errorlevel 1 goto :FAIL
echo.
echo ============================================================
echo EGFR Graph QSAR V3 RUN COMPLETED
echo Return package: results\full\07_final\EGFR_QSAR_REPRODUCIBILITY_RETURN.zip
echo ============================================================
pause
exit /b 0
:FAIL
echo.
echo EGFR Graph QSAR V3 STOPPED WITH AN ERROR. Valid checkpoints are preserved.
echo Check logs\bootstrap.log, logs\pipeline.log and state\LAST_FAILURE.txt.
pause
exit /b 1
